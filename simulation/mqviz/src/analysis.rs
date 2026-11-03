//! Derived measurements shared by the GUI and the headless paths: resonance fits of the sweep cohorts with R in frequency
//! and in m/z terms, and the live observables with their errors and convergence.

use crate::engine::{SpStats, SweepView};
use crate::fit::{self, ResFit};
use crate::stats::{self, Frac, Settle, Val};
use crate::terrain::{self};
use crate::world::*;

/// d f_sec / d(m/z) [kHz per m/z unit] of species `sp` at the body plateau (z = 6 mm, x mode) from the frequency terrain:
/// central difference of the exact Mathieu secular frequency under m/z -> m/z (1 +- 1 %).
pub fn dfdmz(tt: &TerrainTables, sp: &Sp, v_rf: f64, dc: [f64; 4]) -> f64 {
    let eps = 0.01;
    let mk = |k: f64| Sp { mass: sp.mass * k, ..sp.clone() };
    let t = terrain::compute(tt, &[mk(1.0 + eps), mk(1.0 - eps)], v_rf, dc);
    let (fp, fm) = (t.plateau(0), t.plateau(1));
    (fp - fm) / (2.0 * eps * sp.mz())
}

#[derive(Clone, Debug)]
pub struct SweepFit {
    pub label: String,
    pub mz: f64,
    pub fit: ResFit,
    /// local mapping slope [kHz per m/z] and the FWHM in m/z units
    pub dfdmz: f64,
    pub fwhm: Val,
    pub fwhm_mz: Val,
    pub r_f: Val,
    pub r_mz: Val,
    /// peak (baseline + amplitude) as single-pass transmission
    pub peak: Val,
    pub min_ions: u64,
}

/// Fit every species' single-pass transmission curve. Species whose maximum is below 0.5 % are skipped.
pub fn fit_sweep(sv: &SweepView, species: &[Sp], tt: &TerrainTables, v_rf: f64, dc: [f64; 4], boot: usize, seed: u64) -> Vec<Option<SweepFit>> {
    (0..species.len()).map(|s| {
        let ks: Vec<u64> = sv.t.iter().map(|r| r.get(s).map_or(0, |f| f.k)).collect();
        let ns: Vec<u64> = sv.t.iter().map(|r| r.get(s).map_or(0, |f| f.n)).collect();
        if ks.iter().max().copied().unwrap_or(0) == 0 { return None; }
        let ymax = ks.iter().zip(&ns).filter(|(_, &n)| n > 0).map(|(&k, &n)| k as f64 / n as f64).fold(0.0, f64::max);
        if ymax < 0.005 { return None; }
        let f = fit::fit_binomial_curve(&sv.f, &ks, &ns, boot, seed ^ s as u64)?;
        if !fit::is_significant(&f, &sv.f) { return None; }
        let sp = &species[s];
        let slope = dfdmz(tt, sp, v_rf, dc);
        let (_, _, r_mz, r_mz_e) = fit::r_mz(sp.mz(), slope, f.p[3], f.e[3]);
        let fw_mz = f.p[3] / slope.abs().max(1e-12);
        let min_ions = ns.iter().copied().filter(|&n| n > 0).min().unwrap_or(0);
        Some(SweepFit {
            label: sp.label.clone(), mz: sp.mz(), dfdmz: slope, fwhm: Val::new(f.p[3], f.e[3]), fwhm_mz: Val::new(fw_mz, f.e[3] / slope.abs().max(1e-12)),
            r_f: Val::new(f.r, f.r_err), r_mz: Val::new(r_mz, r_mz_e), peak: Val::new(f.peak(), (f.e[0] * f.e[0] + f.e[1] * f.e[1] + 2.0 * f.cov[0][1]).max(0.0).sqrt()), min_ions, fit: f,
        })
    }).collect()
}

// ---------------------------------------------------------------- live observables

/// One observable: value, error, convergence.
#[derive(Clone, Copy, Debug)]
pub struct Obs { pub val: Val, pub settle: Settle }
impl Obs { pub fn none() -> Obs { Obs { val: Val::nan(), settle: Settle::Accumulating { n: 0 } } } }

/// Sum of the counts of all species (charge-flux weights are not applied: counts of ions).
pub fn sum_counts(st: &[SpStats], f: impl Fn(&SpStats) -> [u64; 4]) -> [u64; 4] { st.iter().fold([0; 4], |a, s| { let c = f(s); [a[0] + c[0], a[1] + c[1], a[2] + c[2], a[3] + c[3]] }) }

/// Steady-state efficiency T/(T+C) of the whole mixture, with convergence of its two halves.
pub fn efficiency(st: &[SpStats]) -> Obs {
    let all = sum_counts(st, |s| s.ep);
    let (a, b) = (sum_counts(st, |s| s.half_a), sum_counts(st, |s| s.half_b));
    let f = SpStats::eff(&all);
    if f.n == 0 { return Obs::none(); }
    Obs { val: f.into(), settle: stats::settle_frac(SpStats::eff(&a), SpStats::eff(&b), 10) }
}

/// Efficiency of one species.
pub fn efficiency_sp(s: &SpStats) -> Obs {
    let f = SpStats::eff(&s.ep);
    if f.n == 0 { return Obs::none(); }
    Obs { val: f.into(), settle: stats::settle_frac(SpStats::eff(&s.half_a), SpStats::eff(&s.half_b), 10) }
}

/// Selectivity: the target species' share of the transmitted ions (charge-weighted current share is `current_share`).
pub fn selectivity(st: &[SpStats], target: usize) -> Obs {
    let share = |c: &dyn Fn(&SpStats) -> [u64; 4]| Frac::new(st.get(target).map_or(0, |s| c(s)[0]), st.iter().map(|s| c(s)[0]).sum());
    let f = share(&|s| s.ep);
    if f.n == 0 { return Obs::none(); }
    Obs { val: f.into(), settle: stats::settle_frac(share(&|s| s.half_a), share(&|s| s.half_b), 10) }
}

/// Per-pass transmitted fraction of one species.
pub fn t_pass(s: &SpStats) -> Obs {
    let f = SpStats::t_pass(&s.ep);
    if f.n == 0 { return Obs::none(); }
    Obs { val: f.into(), settle: stats::settle_frac(SpStats::t_pass(&s.half_a), SpStats::t_pass(&s.half_b), 10) }
}

/// Rolling-window rate [ions/s] of outcome `o` (0 T, 1 C, 2 R) with its Poisson error, convergence from the previous window.
pub fn window_rate(s: &SpStats, o: usize, dur: f64) -> Obs {
    Obs { val: stats::rate(s.win[o], dur), settle: stats::settle_rate(s.win_prev[o], s.win[o]) }
}

pub fn window_rate_total(st: &[SpStats], o: usize, dur: f64, z_weight: Option<&[f64]>) -> Obs {
    // sum of Poisson variances (charge-weighted when z_weight is given: current in elementary charges per second)
    let w = |i: usize| z_weight.map_or(1.0, |z| z[i]);
    let (mut v, mut var, mut a, mut b) = (0.0, 0.0, 0u64, 0u64);
    for (i, s) in st.iter().enumerate() {
        v += w(i) * s.win[o] as f64 / dur;
        var += w(i) * w(i) * stats::poisson_sigma(s.win[o]).powi(2) / (dur * dur);
        a += s.win_prev[o]; b += s.win[o];
    }
    Obs { val: Val::new(v, var.sqrt()), settle: stats::settle_rate(a, b) }
}

/// Live per-pass transmission of ubiquitin charge state z for the spectrum model: Jeffreys estimate (k+1/2)/(n+1) and relative error 1/sqrt(k+1/2)
/// (the convention of chsim's `TModel`), or None if z is not simulated.
pub fn live_t_lookup<'a>(sp: &'a [Sp], st: &'a [SpStats]) -> impl Fn(u32) -> Option<(f64, f64)> + 'a {
    move |z| {
        let i = sp.iter().position(|s| s.z as u32 == z)?;
        let f = SpStats::t_pass(&st.get(i)?.ep);
        if f.n == 0 { return None; }
        Some(((f.k as f64 + 0.5) / (f.n as f64 + 1.0), 1.0 / (f.k as f64 + 0.5).sqrt()))
    }
}
