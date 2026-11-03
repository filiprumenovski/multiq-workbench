//! Validation: the live engine (block stepper, global phase tables) against (a) an independent per-ion port of chsim's
//! CPU f64 `xtalk-trace` loop, and (b) chsim's own Metal engine, at the S14 calibrated point.
use crate::presets::{preset, PresetId};
use crate::world::*;
use rayon::prelude::*;
use std::sync::atomic::{AtomicU32, Ordering};

pub const CAL_VRF: f64 = 70.0;
pub const CAL_VAC: f64 = 0.35;
pub const CAL_FAC_KHZ: f64 = 84.25;
pub const CAL_P: f64 = 0.1;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Counts { pub n: usize, pub trans: usize, pub crash: usize, pub ret: usize, pub timeout: usize, pub steps: u64 }
impl Counts {
    pub fn frac(&self, k: usize) -> f64 { k as f64 / self.n.max(1) as f64 }
    pub fn t_over_tc(&self) -> f64 { self.trans as f64 / (self.trans + self.crash).max(1) as f64 }
    /// Binomial standard error of a count out of n.
    pub fn se(&self, k: usize) -> f64 { let p = self.frac(k); (p * (1.0 - p) / self.n.max(1) as f64).sqrt() }
    /// Standard error of T/(T+C) (binomial on the T+C subset).
    pub fn se_eff(&self) -> f64 { let m = (self.trans + self.crash).max(1) as f64; let p = self.t_over_tc(); (p * (1.0 - p) / m).sqrt() }
}

/// Run `n` single-pass ions to completion with the live-engine kernel (all injected at t = 0, one AC tone).
pub fn run_live(w: &World, pool: &rayon::ThreadPool, sp: &Sp, n: usize, vrf: f64, vac: f64, f_ac: f64, p_pa: f64, seed: u64, rf0: bool, progress: Option<&AtomicU32>) -> Counts {
    let sps = [SpPhys::new(sp)];
    let cfg = StepCfg::new(vrf, vac, p_pa);
    let mut ions: Vec<Ion> = (0..n).map(|i| Ion::inject(splitmix(seed.wrapping_add(i as u64 * 0x9E37)), w.a, 0, &sps[0], i as u32, 0, rf0)).collect();
    let mut clock = Clock::new(1, f_ac);
    let mut tabs = Tabs::new();
    let mut c = Counts { n, ..Default::default() };
    let mut alive = n;
    while alive > 0 {
        clock.fill(&mut tabs);
        advance_all(pool, &mut ions, w, &sps, &cfg, &tabs);
        for ion in ions.iter() {
            c.steps += ion.steps_blk as u64;
            match ion.state { ST_TRANS => c.trans += 1, ST_CRASH => c.crash += 1, ST_RETURN => c.ret += 1, ST_TIMEOUT => c.timeout += 1, _ => {} }
        }
        ions.retain(|i| i.state == ST_ALIVE);
        for i in ions.iter_mut() { i.steps_blk = 0; }
        alive = ions.len();
        if let Some(pr) = progress { pr.store(((n - alive) * 1000 / n.max(1)) as u32, Ordering::Relaxed); }
    }
    c
}

/// Independent port of chsim `xtalk-trace` (per-ion loop, own f64 time accumulator, sin(2 pi f t + phase) each step).
pub fn run_reference(w: &World, sp: &Sp, n: usize, vrf: f64, vac: f64, f_ac: f64, p_pa: f64, seed: u64, rf0: bool) -> Counts {
    let ph = SpPhys::new(sp);
    let (qm, s_ion, mr, sig) = (ph.qm, ph.s_ion, ph.mr, ph.sig);
    let ngas = p_pa / (KB * T_BATH);
    let s_gas = (KB * T_BATH / M_GAS).sqrt();
    let a = w.a;
    let dt = DT;
    let res: Vec<(u8, u64)> = (0..n).into_par_iter().map(|id| {
        let mut st: u64 = (id as u64 + 1).wrapping_mul(0x9E3779B97F4A7C15) ^ seed;
        let mut u = move || { st ^= st << 13; st ^= st >> 7; st ^= st << 17; (st >> 11) as f64 / (1u64 << 53) as f64 };
        let gs = |u1: f64, u2: f64| (-2.0 * u1.max(1e-300).ln()).sqrt() * (TWO_PI * u2).cos();
        let r = R_INJECT * u().sqrt();
        let th = TWO_PI * u();
        let mut x = [a + r * th.cos(), a + r * th.sin(), Z_INJECT];
        let mut v = [s_ion * gs(u(), u()), s_ion * gs(u(), u()), s_ion * (-2.0 * u().max(1e-300).ln()).sqrt()];
        let ph_ac = TWO_PI * u();
        let ph_rf = { let r = TWO_PI * u(); if rf0 { 0.0 } else { r } };
        let (mut t, mut steps) = (0.0f64, 0u64);
        while t < T_MAX {
            let (srf, sac) = ((TWO_PI * F_RF * t + ph_rf).sin() * vrf, (TWO_PI * f_ac * t + ph_ac).sin() * vac);
            let (er, ea, ed) = (w.tri(x, 0), w.tri(x, 1), w.tri(x, 2));
            for d in 0..3 { v[d] += qm * dt * (srf * er[d] + sac * ea[d] + ed[d]); x[d] += v[d] * dt; }
            t += dt;
            steps += 1;
            let vm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            let bound = vm + 6.0 * s_gas;
            if u() < ngas * sig * dt * bound {
                let vg = [s_gas * gs(u(), u()), s_gas * gs(u(), u()), s_gas * gs(u(), u())];
                let gv = [v[0] - vg[0], v[1] - vg[1], v[2] - vg[2]];
                let gl = (gv[0] * gv[0] + gv[1] * gv[1] + gv[2] * gv[2]).sqrt();
                if u() * bound < gl {
                    let c = 2.0 * u() - 1.0; let p2 = TWO_PI * u(); let sn = (1.0 - c * c).max(0.0).sqrt();
                    let nn = [sn * p2.cos(), sn * p2.sin(), c];
                    for d in 0..3 { v[d] = v[d] - mr * gv[d] + mr * gl * nn[d]; }
                }
            }
            if x[2] > Z_EXIT { return (ST_TRANS, steps); }
            if x[2] < Z_RETURN { return (ST_RETURN, steps); }
            if w.is_crash(x) { return (ST_CRASH, steps); }
        }
        (ST_TIMEOUT, steps)
    }).collect();
    let mut c = Counts { n, ..Default::default() };
    for (o, s) in res { c.steps += s; match o { ST_TRANS => c.trans += 1, ST_CRASH => c.crash += 1, ST_RETURN => c.ret += 1, _ => c.timeout += 1 } }
    c
}

/// chsim's Metal engine (f32 GPU kernel, `Engine::run_batch`), same physics and drive. Rebuilds Fields from the cache.
pub fn run_chsim_gpu(n: usize, seed: u64) -> Counts {
    use chsim::channel::*;
    let gargs: Vec<String> = vec!["".into(), "".into(), "".into(), "".into(), "".into(), "code".into(), "7.0".into()];
    let eng = Engine::new(Fields::build(&s14_geom(&gargs)));
    let mut ph = Physics::correct(CAL_P);
    ph.t_max = T_MAX;
    let bak = Species { mass: 2.0 * 863.473, z: 2.0, ccs_n2: 425.0 };
    let drive = RfAc { v_rf: CAL_VRF, f_rf: F_RF, v_ac: CAL_VAC, f_ac: CAL_FAC_KHZ * 1e3 };
    let o = eng.run_batch(&[Job { sp: bak, drive, n, seed }], DC_CAL, &ph)[0];
    eprintln!("  chsim Metal raw: transmitted {} crashed {} returned {} left-box {} timeout {}", o.transmitted, o.crashed, o.returned, o.left, o.timeout);
    Counts { n: o.n, trans: o.transmitted, crash: o.crashed + o.left, ret: o.returned, timeout: o.timeout, steps: 0 }
}

/// chsim Metal engine (`Engine::run_batch`, Physics::correct(0.1), t_max 3 ms) at the S14 calibrated point, n = 100000, seed 4242:
/// (n, transmitted, crashed + left-box, returned, timeout). Regenerate with `mqviz --validate 2000 --gpu --ref-n 100000`.
pub const REF_CHSIM_GPU: Option<(usize, usize, usize, usize, usize)> = Some((100_000, 19_669, 3_412, 75_797, 1_122));

pub struct Row { pub name: &'static str, pub live: f64, pub live_se: f64, pub refv: f64, pub ref_se: f64 }
impl Row {
    pub fn pull(&self) -> f64 { (self.live - self.refv) / (self.live_se * self.live_se + self.ref_se * self.ref_se).sqrt().max(1e-12) }
}

pub fn rows(live: &Counts, r: &Counts) -> Vec<Row> {
    vec![
        Row { name: "transmitted", live: live.frac(live.trans), live_se: live.se(live.trans), refv: r.frac(r.trans), ref_se: r.se(r.trans) },
        Row { name: "crashed", live: live.frac(live.crash), live_se: live.se(live.crash), refv: r.frac(r.crash), ref_se: r.se(r.crash) },
        Row { name: "returned", live: live.frac(live.ret), live_se: live.se(live.ret), refv: r.frac(r.ret), ref_se: r.se(r.ret) },
        Row { name: "T/(T+C)", live: live.t_over_tc(), live_se: live.se_eff(), refv: r.t_over_tc(), ref_se: r.se_eff() },
    ]
}

pub fn print_table(title: &str, a_name: &str, b_name: &str, a: &Counts, b: &Counts) {
    println!("  {title}");
    println!("    {:<12} {:>18} {:>18} {:>8}", "outcome", a_name, b_name, "pull");
    for row in rows(a, b) {
        let ok = if row.pull().abs() <= 2.0 { "ok" } else { "DIFF" };
        println!("    {:<12} {:>9.4} ± {:.4} {:>9.4} ± {:.4} {:>7.2}σ {}", row.name, row.live, row.live_se, row.refv, row.ref_se, row.pull(), ok);
    }
}

pub fn embedded_metal() -> Option<Counts> { REF_CHSIM_GPU.map(|(n, t, c, rt, to)| Counts { n, trans: t, crash: c, ret: rt, timeout: to, steps: 0 }) }

/// `mqviz --validate [N] [--gpu] [--ref-n M]`: headless validation at the S14 calibrated point.
pub fn cli(args: &[String]) {
    let n: usize = args.iter().skip_while(|a| *a != "--validate").nth(1).and_then(|s| s.parse().ok()).unwrap_or(2000);
    let ref_n: usize = args.iter().skip_while(|a| *a != "--ref-n").nth(1).and_then(|s| s.parse().ok()).unwrap_or(20000);
    let gpu = args.iter().any(|a| a == "--gpu");
    let t0 = std::time::Instant::now();
    let w = World::build_calibrated();
    println!("fields ready in {:.1} s", t0.elapsed().as_secs_f64());
    let pool = make_pool(default_threads());
    let bak = preset(PresetId::S14).species[0].clone();
    let fac = CAL_FAC_KHZ * 1e3;
    println!("S14 calibrated point: Bak BH3 2+ (m/z {:.1}, CCS {} A²), RF {} V @ 500 kHz, AC ±{} V @ {} kHz, DC {:?}, {} Pa, dt = T_RF/50, crash = conductor node or box wall",
             bak.mz(), bak.ccs, CAL_VRF, CAL_VAC, CAL_FAC_KHZ, DC_CAL, CAL_P);
    let timed = |f: &dyn Fn() -> Counts| { let t = std::time::Instant::now(); let c = f(); (c, t.elapsed().as_secs_f64()) };
    let (live, tl) = timed(&|| run_live(&w, &pool, &bak, n, CAL_VRF, CAL_VAC, fac, CAL_P, 12345, false, None));
    let (live0, _) = timed(&|| run_live(&w, &pool, &bak, n, CAL_VRF, CAL_VAC, fac, CAL_P, 54321, true, None));
    println!("  mqviz live engine (block stepper, global phase tables), n = {}: {:.2} s, {:.1} M ion-steps/s", n, tl, live.steps as f64 / tl / 1e6);
    let (port, tp) = timed(&|| run_reference(&w, &bak, ref_n, CAL_VRF, CAL_VAC, fac, CAL_P, 777, false));
    let (port0, _) = timed(&|| run_reference(&w, &bak, ref_n, CAL_VRF, CAL_VAC, fac, CAL_P, 778, true));
    println!("  xtalk-trace port (independent per-ion loop), n = {}: {:.2} s, {:.1} M ion-steps/s", ref_n, tp, port.steps as f64 / tp / 1e6);
    println!();
    println!("A. random RF phase at injection (physical for a continuous source; what xtalk-trace does)");
    print_table(&format!("live n={} vs xtalk-trace port n={}", n, ref_n), "live", "xtalk port", &live, &port);
    let metal = if gpu { let t = std::time::Instant::now(); let g = run_chsim_gpu(ref_n.max(20000), 4242); println!("  chsim Metal engine (f32 GPU), n = {}: {:.2} s", g.n, t.elapsed().as_secs_f64()); Some(g) } else { embedded_metal() };
    println!();
    println!("B. RF phase 0 for every ion (the convention of chsim's Metal engine: ChIon aux.w starts at 0)");
    print_table(&format!("live n={} vs xtalk-trace port n={}, RF phase 0", n, ref_n), "live", "xtalk port", &live0, &port0);
    if let Some(m) = &metal {
        print_table(&format!("live n={} vs chsim Metal engine n={}", n, m.n), "live (RF0)", "chsim Metal", &live0, m);
        print_table("xtalk-trace port (RF0) vs chsim Metal engine", "port (RF0)", "chsim Metal", &port0, m);
        print_table("xtalk-trace port (random RF) vs chsim Metal engine: effect of the RF-phase convention", "port (random)", "chsim Metal", &port, m);
    }
    println!();
    println!("  documented chsim numbers (docs/blind_prediction_ac_scan.txt, R28/R31): T_pass ≈ 0.17–0.21, T/(T+C) ≈ 0.82–0.87 at AC ±0.35 V");
    println!("  live (random RF phase): T = {:.3} ± {:.3}, T/(T+C) = {:.3} ± {:.3}", live.frac(live.trans), live.se(live.trans), live.t_over_tc(), live.se_eff());
    // run record (JSON) and the latest-validation summary read by the GUI's Fidelity panel
    let vr = crate::engine::ValResult { n, live, live0, port_n: ref_n, port, port0, live_s: tl, port_s: tp, live_mstep_s: live.steps as f64 / tl / 1e6 };
    let env = crate::record::RunEnv::new(format!("{:?}", w.geom));
    let cmd = crate::record::cmd_validate(n, ref_n);
    let out = crate::record::arg_val(args, "--out").map(std::path::PathBuf::from).unwrap_or_else(|| crate::record::new_path("validation", "json"));
    if let Some(d) = out.parent() { let _ = std::fs::create_dir_all(d); }
    match crate::record::write_json(&out, &crate::record::validation_record(&env, &vr, &cmd)).and_then(|_| crate::record::write_json(&crate::record::last_validation_path(), &crate::record::last_validation_summary(&vr, &env.git))) {
        Ok(_) => println!("  record: {}\n  reproduce: {}", out.display(), cmd),
        Err(e) => eprintln!("  record not written: {e}"),
    }
}
