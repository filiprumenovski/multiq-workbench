//! Run records: every run (live segment, sweep, validation) is written as JSON (+ CSV for series) under `mqviz/runs/`
//! with everything needed to reproduce it: parameters, preset, geometry key, seeds, physics settings, git commit
//! (marked dirty if the tree has uncommitted changes), start time, results with errors, and a headless command line.

use crate::analysis::{self, SweepFit};
use crate::engine::{Params, SpStats, Snapshot, SweepCfg, SweepView, ValResult, DEFAULT_SEED};
use crate::json::J;
use crate::presets::{custom_species, PresetId};
use crate::reality::RealRow;
use crate::stats::{Settle, Val};
use crate::validate;
use crate::world::*;
use std::path::PathBuf;

pub const FORMAT: &str = "mqviz-run-1";

#[derive(Clone, Debug, Default)]
pub struct GitInfo { pub hash: String, pub dirty: bool, pub branch: String }

impl GitInfo {
    /// `git rev-parse HEAD` in the crate directory; "dirty" if `git status --porcelain` lists anything (untracked files count; ignored paths such as runs/, target/, cache/ do not).
    pub fn detect() -> GitInfo {
        let run = |args: &[&str]| -> Option<String> {
            let o = std::process::Command::new("git").args(args).current_dir(env!("CARGO_MANIFEST_DIR")).output().ok()?;
            if o.status.success() { Some(String::from_utf8_lossy(&o.stdout).trim().to_string()) } else { None }
        };
        let hash = run(&["rev-parse", "HEAD"]).unwrap_or_else(|| "unknown".into());
        let dirty = run(&["status", "--porcelain"]).map_or(false, |s| !s.trim().is_empty());
        let branch = run(&["rev-parse", "--abbrev-ref", "HEAD"]).unwrap_or_default();
        GitInfo { hash, dirty, branch }
    }
    pub fn short(&self) -> String { format!("{}{}", self.hash.chars().take(9).collect::<String>(), if self.dirty { "-dirty" } else { "" }) }
}

/// Environment of one program run.
#[derive(Clone, Debug)]
pub struct RunEnv {
    pub git: GitInfo,
    pub geom_key: String,
    pub started: String,
    pub threads: usize,
}

impl RunEnv {
    pub fn new(geom_key: String) -> RunEnv { RunEnv { git: GitInfo::detect(), geom_key, started: iso_now(), threads: default_threads() } }
}

// ---------------------------------------------------------------- time and files

fn civil(days: i64) -> (i64, u32, u32) {
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    (if m <= 2 { y + 1 } else { y }, m, d)
}

/// (iso "2026-10-07T10:43:12Z", compact "20261007T104312Z") for a unix time.
pub fn fmt_time(secs: u64) -> (String, String) {
    let (days, rem) = ((secs / 86400) as i64, secs % 86400);
    let (y, m, d) = civil(days);
    let (hh, mm, ss) = (rem / 3600, rem % 3600 / 60, rem % 60);
    (format!("{y:04}-{m:02}-{d:02}T{hh:02}:{mm:02}:{ss:02}Z"), format!("{y:04}{m:02}{d:02}T{hh:02}{mm:02}{ss:02}Z"))
}
pub fn now_secs() -> u64 { std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.as_secs()) }
pub fn iso_now() -> String { fmt_time(now_secs()).0 }

pub fn runs_dir() -> PathBuf {
    let p = std::env::var_os("MQVIZ_RUNS").map(PathBuf::from).unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("runs"));
    let _ = std::fs::create_dir_all(&p);
    p
}

/// `<runs>/<timestamp>_<kind>.<ext>`, with a numeric suffix if the name exists.
pub fn new_path(kind: &str, ext: &str) -> PathBuf {
    let stamp = fmt_time(now_secs()).1;
    let dir = runs_dir();
    let mut p = dir.join(format!("{stamp}_{kind}.{ext}"));
    let mut i = 2;
    while p.exists() { p = dir.join(format!("{stamp}_{kind}_{i}.{ext}")); i += 1; }
    p
}

pub fn write_json(path: &std::path::Path, j: &J) -> std::io::Result<()> { std::fs::write(path, j.pretty()) }

// ---------------------------------------------------------------- parameters <-> JSON

pub fn preset_name(id: PresetId) -> &'static str { match id { PresetId::S14 => "s14", PresetId::Fig5 => "fig5", PresetId::Custom => "custom" } }
pub fn preset_from(name: &str) -> Option<PresetId> { match name.to_lowercase().as_str() { "s14" => Some(PresetId::S14), "fig5" | "fig.5" | "fig-5" => Some(PresetId::Fig5), "custom" => Some(PresetId::Custom), _ => None } }

pub fn seed_str(s: u64) -> String { format!("0x{:016x}", s) }
pub fn parse_seed(s: &str) -> Option<u64> { let t = s.trim(); if let Some(h) = t.strip_prefix("0x") { u64::from_str_radix(h, 16).ok() } else { t.parse().ok() } }

pub fn sp_j(s: &Sp) -> J { J::obj().with("label", s.label.clone()).with("mass_da", s.mass).with("z", s.z).with("mz", s.mz()).with("ccs_n2_A2", s.ccs).with("weight", s.weight) }

pub fn params_j(p: &Params) -> J {
    J::obj().with("preset", preset_name(p.preset)).with("species", J::Arr(p.species.iter().map(sp_j).collect()))
        .with("rf_amplitude_V", p.v_rf).with("ac_amplitude_V", p.v_ac).with("f_ac_kHz", p.f_ac_khz).with("gate_dc_V", p.gate_v)
        .with("pressure_Pa", p.pressure_pa).with("injection_rate_per_s", p.inj_rate).with("trap_residence_wall_s", p.residence_s)
        .with("recirculate", p.recirculate).with("averaging_window_ms", p.window_ms).with("seed", seed_str(p.seed))
}

/// Parameters from a saved state / record (missing keys keep `base`).
pub fn params_from_j(j: &J, base: &Params) -> Result<Params, String> {
    let pj = j.get("params").unwrap_or(j);
    let mut p = base.clone();
    if let Some(n) = pj.s("preset") { p.preset = preset_from(n).ok_or(format!("unknown preset {n}"))?; }
    if let Some(a) = pj.arr("species") {
        let mut sp = vec![];
        for s in a {
            sp.push(Sp { label: s.s("label").unwrap_or("?").to_string(), mass: s.f("mass_da").ok_or("species without mass_da")?, z: s.f("z").ok_or("species without z")?, ccs: s.f("ccs_n2_A2").ok_or("species without ccs")?, weight: s.f("weight").unwrap_or(1.0) });
        }
        if !sp.is_empty() { p.species = sp; }
    }
    let g = |k: &str, d: f64| pj.f(k).unwrap_or(d);
    p.v_rf = g("rf_amplitude_V", p.v_rf); p.v_ac = g("ac_amplitude_V", p.v_ac); p.f_ac_khz = g("f_ac_kHz", p.f_ac_khz); p.gate_v = g("gate_dc_V", p.gate_v);
    p.pressure_pa = g("pressure_Pa", p.pressure_pa); p.inj_rate = g("injection_rate_per_s", p.inj_rate); p.residence_s = g("trap_residence_wall_s", p.residence_s);
    p.window_ms = g("averaging_window_ms", p.window_ms);
    if let Some(b) = pj.b("recirculate") { p.recirculate = b; }
    if let Some(s) = pj.s("seed").and_then(parse_seed) { p.seed = s; }
    Ok(p)
}

pub fn save_state(p: &Params, sweep: &SweepCfg, env: &RunEnv) -> J {
    J::obj().with("format", "mqviz-state-1").with("saved", iso_now()).with("git", git_j(&env.git))
        .with("params", params_j(p)).with("sweep", sweep_cfg_j(sweep))
}

pub fn sweep_cfg_j(c: &SweepCfg) -> J { J::obj().with("f_lo_kHz", c.f_lo_khz).with("f_hi_kHz", c.f_hi_khz).with("step_kHz", c.step_khz).with("ions_per_point", c.ions_per_point).with("sweeps", c.sweeps) }
pub fn sweep_cfg_from_j(j: &J, base: &SweepCfg) -> SweepCfg {
    let s = j.get("sweep").unwrap_or(j);
    SweepCfg { f_lo_khz: s.f("f_lo_kHz").unwrap_or(base.f_lo_khz), f_hi_khz: s.f("f_hi_kHz").unwrap_or(base.f_hi_khz), step_khz: s.f("step_kHz").unwrap_or(base.step_khz),
               ions_per_point: s.f("ions_per_point").map_or(base.ions_per_point, |x| x as usize), sweeps: s.f("sweeps").map_or(base.sweeps, |x| x as usize) }
}

pub fn git_j(g: &GitInfo) -> J { J::obj().with("commit", g.hash.clone()).with("dirty", g.dirty).with("branch", g.branch.clone()) }

/// What the integrator does, as run (all constants of the engine).
pub fn physics_j(p: &Params) -> J {
    J::obj().with("integrator", "f64 kick-drift leapfrog, full force").with("dt_s", DT).with("dt_definition", "T_RF/50").with("f_rf_Hz", F_RF)
        .with("gas", "N2 (28 u), Maxwellian bath, hard-sphere null-collision, centre-of-mass scattering").with("bath_temperature_K", T_BATH).with("pressure_Pa", p.pressure_pa)
        .with("ccs", "real N2 CCS per species").with("fields", "trilinear interpolation, 0.1 mm grid, calibrated code geometry")
        .with("injection", J::obj().with("model", "effusive thermal, xtalk-trace").with("z_mm", Z_INJECT * 1e3).with("r_disc_mm", R_INJECT * 1e3))
        .with("outcomes", J::obj().with("transmitted_z_gt_mm", (Z_EXIT * 1e4).round() / 10.0).with("returned_z_lt_mm", Z_RETURN * 1e3).with("crash", "conductor node or box wall").with("t_max_ms", T_MAX * 1e3))
        .with("dc_calibrated_V", DC_CAL.to_vec()).with("not_modelled", vec!["space charge in the channel", "trap-side delivery (effusive injection at z = 1 mm, recirculation = delay only)"])
}

pub fn meta_j(kind: &str, env: &RunEnv, p: &Params, extra_cmd: &str) -> J {
    J::obj().with("format", FORMAT).with("kind", kind).with("tool", "mqviz").with("version", env!("CARGO_PKG_VERSION"))
        .with("git", git_j(&env.git)).with("program_started", env.started.clone()).with("written", iso_now())
        .with("preset", preset_name(p.preset)).with("geometry_key", env.geom_key.clone()).with("geometry_hash", format!("{:016x}", fnv(&env.geom_key)))
        .with("seeds", J::obj().with("simulation", seed_str(p.seed)).with("note", "per-ion RNG streams derive from the simulation seed and the ion index"))
        .with("physics", physics_j(p)).with("params", params_j(p)).with("threads", env.threads).with("command", extra_cmd)
}

pub fn fnv(s: &str) -> u64 { s.bytes().fold(0xcbf29ce484222325u64, |h, b| (h ^ b as u64).wrapping_mul(0x100000001b3)) }

// ---------------------------------------------------------------- results

pub fn val_j(v: Val) -> J { J::obj().with("value", v.v).with("sigma", v.e) }
pub fn settle_j(s: &Settle) -> J { J::obj().with("state", match s { Settle::Settled { .. } => "settled", Settle::Accumulating { .. } => "accumulating", Settle::Drifting { .. } => "drifting" }).with("ions", s.n()) }
fn obs_j(o: &analysis::Obs) -> J { val_j(o.val).with("convergence", settle_j(&o.settle)) }

pub fn fit_j(f: &SweepFit) -> J {
    let ft = &f.fit;
    J::obj().with("species", f.label.clone()).with("shape", ft.shape.name()).with("converged", ft.converged)
        .with("baseline", val_j(Val::new(ft.p[0], ft.e[0]))).with("amplitude", val_j(Val::new(ft.p[1], ft.e[1])))
        .with("f0_kHz", val_j(Val::new(ft.p[2], ft.e[2]))).with("fwhm_kHz", val_j(f.fwhm)).with("peak_transmission", val_j(f.peak))
        .with("R_f_over_fwhm", val_j(f.r_f)).with("R_bootstrap", ft.r_boot.map_or(J::Null, |(m, s)| J::obj().with("mean", m).with("sigma", s)))
        .with("dfdmz_kHz_per_mz", f.dfdmz).with("mz", f.mz).with("fwhm_mz", val_j(f.fwhm_mz)).with("R_mz", val_j(f.r_mz))
        .with("chi2", ft.chi2).with("chi2_flat", ft.chi2_flat).with("delta_chi2_vs_flat", ft.delta_chi2()).with("dof", ft.dof as f64).with("chi2_red", ft.chi2_red()).with("error_scale", ft.scale).with("points", ft.n_pts).with("min_ions_per_point", f.min_ions)
}

pub fn sweep_record(env: &RunEnv, p: &Params, sv: &SweepView, fits: &[Option<SweepFit>], cmd: &str) -> J {
    let species = &p.species;
    let pts: Vec<J> = sv.f.iter().enumerate().map(|(i, f)| {
        J::obj().with("f_kHz", *f).with("species", J::Arr(species.iter().enumerate().map(|(s, sp)| {
            let c = sv.t[i][s];
            J::obj().with("label", sp.label.clone()).with("transmitted", c.k).with("finished", c.n).with("T", c.p()).with("sigma", c.sigma())
        }).collect()))
    }).collect();
    meta_j("sweep", env, p, cmd).with("sweep", sv.cfg.as_ref().map_or(J::Null, sweep_cfg_j)).with("sweeps_completed", sv.sweeps_done).with("results", J::obj()
        .with("points", J::Arr(pts)).with("fits", J::Arr(fits.iter().flatten().map(fit_j).collect()))
        .with("error_model", "T = k/n per tagged cohort, sigma = Wilson-score 1-sigma; fit: weighted LM, binomial variances from the model, covariance scaled by sqrt(chi2/dof) if > 1"))
}

pub fn sweep_csv(sv: &SweepView, species: &[Sp]) -> String {
    let mut s = String::from("f_kHz,species,transmitted,finished,T,sigma\n");
    for (i, f) in sv.f.iter().enumerate() {
        for (k, sp) in species.iter().enumerate() {
            let c = sv.t[i][k];
            s.push_str(&format!("{},{},{},{},{},{}\n", f, sp.label, c.k, c.n, c.p(), c.sigma()));
        }
    }
    s
}

pub fn species_stats_j(sp: &Sp, st: &SpStats, dur: f64) -> J {
    let rate = |o: usize| analysis::window_rate(st, o, dur);
    let pa = |v: Val| Val::new(v.v * sp.z * QE * 1e12, v.e * sp.z * QE * 1e12);
    J::obj().with("species", sp_j(sp))
        .with("counts_since_epoch", J::obj().with("transmitted", st.ep[0]).with("crashed", st.ep[1]).with("returned", st.ep[2]).with("timeout", st.ep[3]))
        .with("T_per_pass", obs_j(&analysis::t_pass(st))).with("efficiency_T_over_TC", obs_j(&analysis::efficiency_sp(st)))
        .with("window_rates_per_s", J::obj().with("transmitted", obs_j(&rate(0))).with("crashed", obs_j(&rate(1))).with("returned", obs_j(&rate(2))))
        .with("window_currents_pA", J::obj().with("transmitted", val_j(pa(rate(0).val))).with("crashed", val_j(pa(rate(1).val))).with("returned", val_j(pa(rate(2).val))))
}

pub fn rows_j(rows: &[RealRow]) -> J {
    J::Arr(rows.iter().map(|r| J::obj().with("observable", r.name.clone()).with("model", val_j(r.ours)).with("theirs", r.theirs.map_or(J::Null, val_j)).with("theirs_text", r.theirs_txt.clone())
        .with("pull_sigma", r.pull).with("mark", r.mark()).with("convergence", r.settle.as_ref().map_or(J::Null, settle_j)).with("source", r.src)).collect())
}

pub fn live_record(env: &RunEnv, p: &Params, snap: &Snapshot, rtf: f64, target: usize, rows: &[RealRow], cmd: &str) -> J {
    let sp = &snap.species;
    let eff = analysis::efficiency(&snap.stats);
    let sel = analysis::selectivity(&snap.stats, target);
    let delay_ms = p.residence_s * rtf * 1e3;
    meta_j("live", env, p, cmd).with("segment", J::obj().with("sim_time_ms", snap.t_sim * 1e3).with("epoch_start_ms", snap.epoch_start_ms).with("epoch_ms", snap.epoch_ms)
            .with("window_ms", snap.win_dur * 1e3).with("realtime_factor", rtf).with("mean_trap_delay_sim_ms", delay_ms).with("n_in_flight", snap.n_flight)
            .with("note", "fractions accumulate over the epoch (since the last parameter change); rates over the rolling window; the GUI recirculation delay depends on the wall-clock realtime factor, so a GUI segment reproduces statistically, not bit-wise"))
        .with("results", J::obj().with("species", J::Arr(sp.iter().zip(&snap.stats).map(|(s, st)| species_stats_j(s, st, snap.win_dur)).collect()))
            .with("efficiency_T_over_TC_all", obs_j(&eff)).with("selectivity_target", J::obj().with("species", sp.get(target).map_or("", |s| s.label.as_str())).with("share_of_transmitted", obs_j(&sel)))
            .with("reality", rows_j(rows)))
}

pub fn counts_j(c: &validate::Counts) -> J {
    J::obj().with("n", c.n).with("transmitted", c.trans).with("crashed", c.crash).with("returned", c.ret).with("timeout", c.timeout)
        .with("frac_transmitted", val_j(Val::new(c.frac(c.trans), c.se(c.trans)))).with("frac_crashed", val_j(Val::new(c.frac(c.crash), c.se(c.crash))))
        .with("frac_returned", val_j(Val::new(c.frac(c.ret), c.se(c.ret)))).with("T_over_TC", val_j(Val::new(c.t_over_tc(), c.se_eff())))
}

pub fn validation_record(env: &RunEnv, r: &ValResult, cmd: &str) -> J {
    let p = Params::from_preset(PresetId::S14);
    let rows = |a: &validate::Counts, b: &validate::Counts| J::Arr(validate::rows(a, b).iter().map(|x| J::obj().with("outcome", x.name).with("live", x.live).with("live_se", x.live_se).with("ref", x.refv).with("ref_se", x.ref_se).with("pull_sigma", x.pull())).collect());
    let mut res = J::obj().with("n_live", r.n).with("live_random_rf", counts_j(&r.live)).with("live_rf0", counts_j(&r.live0)).with("xtalk_port_n", r.port_n)
        .with("xtalk_port_random_rf", counts_j(&r.port)).with("xtalk_port_rf0", counts_j(&r.port0))
        .with("live_vs_port_random", rows(&r.live, &r.port)).with("live_rf0_vs_port_rf0", rows(&r.live0, &r.port0));
    if let Some(m) = validate::embedded_metal() { res.set("chsim_metal_n100000", counts_j(&m)); res.set("live_rf0_vs_chsim_metal", rows(&r.live0, &m)); }
    res.set("seconds_live_per_run", r.live_s); res.set("M_ion_steps_per_s", r.live_mstep_s);
    meta_j("validation", env, &p, cmd).with("point", J::obj().with("species", "Bak BH3 2+ m/z 863.47 CCS 425").with("rf_V", validate::CAL_VRF).with("ac_V", validate::CAL_VAC).with("f_ac_kHz", validate::CAL_FAC_KHZ).with("pressure_Pa", validate::CAL_P).with("dc_V", DC_CAL.to_vec()))
        .with("results", res)
}

/// `runs/last_validation.json`: the latest validation numbers and date (read by the Fidelity panel).
pub fn last_validation_path() -> PathBuf { runs_dir().join("last_validation.json") }

pub fn last_validation_summary(r: &ValResult, git: &GitInfo) -> J {
    let rows = validate::rows(&r.live, &r.port);
    J::obj().with("date", iso_now()).with("git", git.short()).with("n_live", r.n).with("n_port", r.port_n)
        .with("rows", J::Arr(rows.iter().map(|x| J::obj().with("outcome", x.name).with("live", x.live).with("live_se", x.live_se).with("ref", x.refv).with("ref_se", x.ref_se).with("pull", x.pull())).collect()))
}

// ---------------------------------------------------------------- time series CSV

#[derive(Clone, Debug)]
pub struct HistPt {
    pub t_ms: f64,
    /// per species: [T, C, R] rates (ions/s) and their Poisson errors
    pub r: Vec<[f64; 3]>,
    pub e: Vec<[f64; 3]>,
}

pub fn series_csv(h: &[HistPt], species: &[Sp]) -> String {
    let mut s = String::from("t_sim_ms");
    for sp in species { for o in ["T", "C", "R"] { s.push_str(&format!(",{l}_{o}_per_s,{l}_{o}_sigma", l = sp.label.replace(' ', "_"))); } }
    s.push('\n');
    for p in h {
        s.push_str(&format!("{}", p.t_ms));
        for i in 0..species.len() { for o in 0..3 { s.push_str(&format!(",{},{}", p.r.get(i).map_or(f64::NAN, |r| r[o]), p.e.get(i).map_or(f64::NAN, |r| r[o]))); } }
        s.push('\n');
    }
    s
}

// ---------------------------------------------------------------- command lines

fn fnum(x: f64) -> String { let t = format!("{}", x); t }

/// Species override for custom presets: "mz:z:ccs".
fn species_arg(p: &Params) -> Option<String> {
    if p.preset != PresetId::Custom { return None; }
    p.species.first().map(|s| format!("{}:{}:{}", s.mz(), s.z, s.ccs))
}

fn common_args(p: &Params) -> String {
    let mut a = format!("{} --rf {} --vac {} --gate {} --pressure {}", preset_name(p.preset), fnum(p.v_rf), fnum(p.v_ac), fnum(p.gate_v), fnum(p.pressure_pa));
    if let Some(s) = species_arg(p) { a.push_str(&format!(" --species {}", s)); }
    if p.seed != DEFAULT_SEED { a.push_str(&format!(" --seed {}", seed_str(p.seed))); }
    a
}

pub fn cmd_sweep(p: &Params, c: &SweepCfg, inj_rate: f64) -> String {
    format!("mqviz --sweep {} --n {} --flo {} --fhi {} --step {} --sweeps {} --inj-rate {} --out sweep.json", common_args(p), c.ions_per_point, fnum(c.f_lo_khz), fnum(c.f_hi_khz), fnum(c.step_khz), c.sweeps, fnum(inj_rate))
}
pub fn cmd_live(p: &Params, sim_ms: f64, rtf: f64) -> String {
    format!("mqviz --live {} --fac {} --sim-ms {} --rtf {:.3e} --inj-rate {} --out live.json", common_args(p), fnum(p.f_ac_khz), fnum((sim_ms * 10.0).round() / 10.0), rtf, fnum(p.inj_rate))
}
pub fn cmd_validate(n: usize, ref_n: usize) -> String { format!("mqviz --validate {} --ref-n {} --out validation.json", n, ref_n) }

// ---------------------------------------------------------------- headless argument parsing

pub fn arg_val(args: &[String], key: &str) -> Option<String> { args.iter().position(|a| a == key).and_then(|i| args.get(i + 1).cloned()) }
pub fn arg_f(args: &[String], key: &str) -> Option<f64> { arg_val(args, key).and_then(|s| s.parse().ok()) }

/// `--sweep|--live <preset> [--rf V] [--vac V] [--fac kHz] [--gate V] [--pressure Pa] [--seed S] [--inj-rate r] [--species mz:z:ccs]`
pub fn params_from_args(args: &[String], mode_flag: &str) -> Result<Params, String> {
    let name = args.iter().skip_while(|a| *a != mode_flag).nth(1).filter(|a| !a.starts_with("--")).cloned().unwrap_or_else(|| "s14".into());
    let id = preset_from(&name).ok_or(format!("unknown preset '{name}' (s14 | fig5 | custom)"))?;
    let mut p = Params::from_preset(id);
    if id == PresetId::Custom {
        if let Some(sp) = arg_val(args, "--species") {
            let v: Vec<f64> = sp.split(':').filter_map(|x| x.parse().ok()).collect();
            if v.len() != 3 { return Err("--species mz:z:ccs".into()); }
            p.species = vec![custom_species(v[0], v[1], v[2])];
        }
    }
    if let Some(x) = arg_f(args, "--rf") { p.v_rf = x; }
    if let Some(x) = arg_f(args, "--vac") { p.v_ac = x; }
    if let Some(x) = arg_f(args, "--fac") { p.f_ac_khz = x; }
    if let Some(x) = arg_f(args, "--gate") { p.gate_v = x; }
    if let Some(x) = arg_f(args, "--pressure") { p.pressure_pa = x; }
    if let Some(x) = arg_f(args, "--inj-rate") { p.inj_rate = x; }
    if let Some(s) = arg_val(args, "--seed") { p.seed = parse_seed(&s).ok_or("bad --seed")?; }
    Ok(p)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn time_formatting() {
        assert_eq!(fmt_time(0).0, "1970-01-01T00:00:00Z");
        assert_eq!(fmt_time(1_791_367_392).0, "2026-10-07T10:03:12Z");
        assert_eq!(fmt_time(951_782_400).0, "2000-02-29T00:00:00Z");
    }
    #[test]
    fn state_roundtrip() {
        let mut p = Params::from_preset(PresetId::Fig5);
        p.v_rf = 57.5; p.seed = 0xDEAD_BEEF_0123_4567; p.recirculate = false;
        let env = RunEnv { git: GitInfo { hash: "abc".into(), dirty: true, branch: "x".into() }, geom_key: "g".into(), started: "t".into(), threads: 3 };
        let j = save_state(&p, &SweepCfg { f_lo_khz: 30.0, f_hi_khz: 70.0, step_khz: 0.5, ions_per_point: 123, sweeps: 2 }, &env);
        let t = j.pretty();
        let q = params_from_j(&J::parse(&t).unwrap(), &Params::from_preset(PresetId::S14)).unwrap();
        assert!(q == p, "state round trip");
        let c = sweep_cfg_from_j(&J::parse(&t).unwrap(), &SweepCfg { f_lo_khz: 0.0, f_hi_khz: 0.0, step_khz: 0.0, ions_per_point: 0, sweeps: 0 });
        assert_eq!((c.ions_per_point, c.sweeps), (123, 2));
    }
    #[test]
    fn command_parses_back_to_the_same_params() {
        let mut p = Params::from_preset(PresetId::S14);
        p.v_ac = 0.3; p.seed = 77;
        let cmd = cmd_sweep(&p, &SweepCfg { f_lo_khz: 70.0, f_hi_khz: 100.0, step_khz: 1.0, ions_per_point: 3000, sweeps: 1 }, 8.0e6);
        let args: Vec<String> = cmd.split_whitespace().skip(1).map(String::from).collect();
        let mut q = params_from_args(&args, "--sweep").unwrap();
        q.inj_rate = p.inj_rate;
        assert!(q == p, "cmd '{cmd}' -> {:?}", q.v_ac);
    }
}
