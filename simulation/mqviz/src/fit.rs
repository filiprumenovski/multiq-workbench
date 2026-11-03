//! Resonance-curve fit: Lorentzian or Gaussian on a flat baseline, y = b + A shape((f - f0)/w), with w the FWHM.
//! Weighted Levenberg-Marquardt (bounded, three starting widths); parameter errors from the weighted covariance
//! (scaled by sqrt(chi2/dof) when chi2/dof > 1); sigma_i of binomial data follow the fitted model (iterated), not the noisy
//! observed fraction, which would bias the fit towards downward fluctuations. R = f0 / FWHM carries the full covariance;
//! a parametric bootstrap gives a second, assumption-light error on R.

use crate::stats::Rng;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Shape { Lorentz, Gauss }
impl Shape {
    pub fn name(&self) -> &'static str { match self { Shape::Lorentz => "Lorentzian", Shape::Gauss => "Gaussian" } }
    /// Peak-normalised shape at offset `d` from the centre for FWHM `w`.
    pub fn g(&self, d: f64, w: f64) -> f64 {
        match self {
            Shape::Lorentz => { let u = 2.0 * d / w; 1.0 / (1.0 + u * u) }
            Shape::Gauss => (-4.0 * std::f64::consts::LN_2 * d * d / (w * w)).exp(),
        }
    }
    /// d g / d d  and  d g / d w
    fn dg(&self, d: f64, w: f64) -> (f64, f64) {
        match self {
            Shape::Lorentz => { let u = 2.0 * d / w; let g = 1.0 / (1.0 + u * u); let g2 = g * g; (-g2 * 2.0 * u * 2.0 / w, g2 * 2.0 * u * u / w) }
            Shape::Gauss => { let g = self.g(d, w); (g * (-8.0 * std::f64::consts::LN_2 * d / (w * w)), g * 8.0 * std::f64::consts::LN_2 * d * d / (w * w * w)) }
        }
    }
}

/// Parameters [b, A, f0, w].
pub type P4 = [f64; 4];

pub fn model(sh: Shape, p: &P4, x: f64) -> f64 { p[0] + p[1] * sh.g(x - p[2], p[3]) }

#[derive(Clone, Debug)]
pub struct ResFit {
    pub shape: Shape,
    pub p: P4,
    /// 1-sigma errors of [b, A, f0, w]
    pub e: P4,
    pub cov: [[f64; 4]; 4],
    pub chi2: f64,
    /// chi2 of the flat (no resonance) model on the same points, for the significance of the peak
    pub chi2_flat: f64,
    pub dof: i64,
    /// factor applied to the covariance (max(1, sqrt(chi2/dof)))
    pub scale: f64,
    pub n_pts: usize,
    /// R = f0 / FWHM and its error (delta method, full covariance)
    pub r: f64,
    pub r_err: f64,
    /// parametric bootstrap: (mean, sigma) of R over the replicas that converged
    pub r_boot: Option<(f64, f64)>,
    pub f0_boot: Option<f64>,
    pub converged: bool,
}

impl ResFit {
    pub fn peak(&self) -> f64 { self.p[0] + self.p[1] }
    /// Improvement over a flat curve; a resonance is called significant above `MIN_DELTA_CHI2`.
    pub fn delta_chi2(&self) -> f64 { self.chi2_flat - self.chi2 }
    pub fn chi2_red(&self) -> f64 { if self.dof > 0 { self.chi2 / self.dof as f64 } else { f64::NAN } }
    /// Value of the fit at x.
    pub fn at(&self, x: f64) -> f64 { model(self.shape, &self.p, x) }
}

fn solve4(a: [[f64; 4]; 4], b: [f64; 4]) -> Option<[f64; 4]> {
    let mut m = [[0.0; 5]; 4];
    for i in 0..4 { for j in 0..4 { m[i][j] = a[i][j]; } m[i][4] = b[i]; }
    for c in 0..4 {
        let piv = (c..4).max_by(|&i, &j| m[i][c].abs().partial_cmp(&m[j][c].abs()).unwrap())?;
        if m[piv][c].abs() < 1e-300 { return None; }
        m.swap(c, piv);
        for r in 0..4 { if r != c { let f = m[r][c] / m[c][c]; for k in c..5 { m[r][k] -= f * m[c][k]; } } }
    }
    Some([m[0][4] / m[0][0], m[1][4] / m[1][1], m[2][4] / m[2][2], m[3][4] / m[3][3]])
}

fn invert4(a: [[f64; 4]; 4]) -> Option<[[f64; 4]; 4]> {
    let mut out = [[0.0; 4]; 4];
    for c in 0..4 { let mut e = [0.0; 4]; e[c] = 1.0; let s = solve4(a, e)?; for r in 0..4 { out[r][c] = s[r]; } }
    Some(out)
}

fn chi2_of(sh: Shape, p: &P4, xs: &[f64], ys: &[f64], sg: &[f64]) -> f64 {
    xs.iter().zip(ys).zip(sg).map(|((&x, &y), &s)| ((y - model(sh, p, x)) / s).powi(2)).sum()
}

struct Bounds { x_lo: f64, x_hi: f64, w_lo: f64, w_hi: f64 }

fn clamp(p: &mut P4, b: &Bounds) {
    p[2] = p[2].clamp(b.x_lo, b.x_hi);
    p[3] = p[3].clamp(b.w_lo, b.w_hi);
}

/// Levenberg-Marquardt from a start; returns (params, chi2, converged).
fn lm(sh: Shape, start: P4, xs: &[f64], ys: &[f64], sg: &[f64], b: &Bounds) -> (P4, f64, bool) {
    let mut p = start;
    clamp(&mut p, b);
    let mut chi = chi2_of(sh, &p, xs, ys, sg);
    let mut lam = 1e-3;
    let mut conv = false;
    for _ in 0..200 {
        let mut jtj = [[0.0; 4]; 4];
        let mut jtr = [0.0; 4];
        for i in 0..xs.len() {
            let d = xs[i] - p[2];
            let (gd, gw) = sh.dg(d, p[3]);
            let g = sh.g(d, p[3]);
            let j = [1.0 / sg[i], g / sg[i], -p[1] * gd / sg[i], p[1] * gw / sg[i]];
            let r = (ys[i] - model(sh, &p, xs[i])) / sg[i];
            for a in 0..4 { jtr[a] += j[a] * r; for c in 0..4 { jtj[a][c] += j[a] * j[c]; } }
        }
        let mut improved = false;
        for _ in 0..12 {
            let mut a = jtj;
            for k in 0..4 { a[k][k] += lam * jtj[k][k].max(1e-12); }
            let Some(dp) = solve4(a, jtr) else { lam *= 10.0; continue };
            let mut q = [p[0] + dp[0], p[1] + dp[1], p[2] + dp[2], p[3] + dp[3]];
            clamp(&mut q, b);
            let c2 = chi2_of(sh, &q, xs, ys, sg);
            if c2.is_finite() && c2 < chi {
                let rel = (chi - c2) / chi.max(1e-300);
                p = q; chi = c2; lam = (lam * 0.3).max(1e-9); improved = true;
                if rel < 1e-10 { conv = true; }
                break;
            }
            lam *= 5.0;
        }
        if !improved { conv = true; }
        if conv { break; }
    }
    (p, chi, conv)
}

fn half_max_width(xs: &[f64], ys: &[f64], base: f64, imax: usize) -> f64 {
    let half = base + 0.5 * (ys[imax] - base);
    let mut lo = xs[0]; let mut hi = xs[xs.len() - 1];
    for i in (1..=imax).rev() { if ys[i - 1] < half { lo = xs[i - 1] + (half - ys[i - 1]) / (ys[i] - ys[i - 1]) * (xs[i] - xs[i - 1]); break; } }
    for i in imax..xs.len() - 1 { if ys[i + 1] < half { hi = xs[i] + (ys[i] - half) / (ys[i] - ys[i + 1]) * (xs[i + 1] - xs[i]); break; } }
    (hi - lo).max(1e-6)
}

fn bounds_of(xs: &[f64]) -> Bounds {
    let (lo, hi) = (xs[0], xs[xs.len() - 1]);
    let dx = (hi - lo) / (xs.len() - 1) as f64;
    Bounds { x_lo: lo - 0.2 * (hi - lo), x_hi: hi + 0.2 * (hi - lo), w_lo: 0.3 * dx, w_hi: 3.0 * (hi - lo) }
}

fn fit_shape(sh: Shape, xs: &[f64], ys: &[f64], sg: &[f64], init: Option<P4>) -> Option<(P4, f64, bool)> {
    let b = bounds_of(xs);
    // smoothed argmax (3-point weighted by 1/sigma^2 is overkill; plain 3-point mean)
    let n = xs.len();
    let sm: Vec<f64> = (0..n).map(|i| { let (a, c) = (i.saturating_sub(1), (i + 1).min(n - 1)); (ys[a] + ys[i] + ys[c]) / 3.0 }).collect();
    let imax = (0..n).max_by(|&i, &j| sm[i].partial_cmp(&sm[j]).unwrap())?;
    let mut base = ys.iter().cloned().fold(f64::INFINITY, f64::min).max(0.0);
    if !base.is_finite() { base = 0.0; }
    let amp = (ys[imax] - base).max(1e-9);
    let w0 = half_max_width(xs, ys, base, imax);
    let mut starts: Vec<P4> = vec![];
    if let Some(i) = init { starts.push(i); }
    for f in [1.0, 0.5, 2.0] { starts.push([base, amp, xs[imax], (w0 * f).clamp(b.w_lo, b.w_hi)]); }
    let mut best: Option<(P4, f64, bool)> = None;
    for s in starts { let r = lm(sh, s, xs, ys, sg, &b); if best.as_ref().map_or(true, |bb| r.1 < bb.1) { best = Some(r); } }
    best
}

fn covariance(sh: Shape, p: &P4, xs: &[f64], sg: &[f64]) -> Option<[[f64; 4]; 4]> {
    let mut jtj = [[0.0; 4]; 4];
    for i in 0..xs.len() {
        let d = xs[i] - p[2];
        let (gd, gw) = sh.dg(d, p[3]);
        let j = [1.0 / sg[i], sh.g(d, p[3]) / sg[i], -p[1] * gd / sg[i], p[1] * gw / sg[i]];
        for a in 0..4 { for c in 0..4 { jtj[a][c] += j[a] * j[c]; } }
    }
    invert4(jtj)
}

fn r_and_err(p: &P4, cov: &[[f64; 4]; 4]) -> (f64, f64) {
    let (f0, w) = (p[2], p[3]);
    let r = f0 / w;
    let (gf, gw) = (1.0 / w, -f0 / (w * w));
    let v = gf * gf * cov[2][2] + gw * gw * cov[3][3] + 2.0 * gf * gw * cov[2][3];
    (r, v.max(0.0).sqrt())
}

/// Fit with given sigma_i (shape fixed).
#[allow(dead_code)]
pub fn fit_shape_sigma(sh: Shape, xs: &[f64], ys: &[f64], sg: &[f64]) -> Option<ResFit> {
    if xs.len() < 6 || xs.len() != ys.len() || xs.len() != sg.len() { return None; }
    let (p, chi2, converged) = fit_shape(sh, xs, ys, sg, None)?;
    finish(sh, p, chi2, converged, xs, sg)
}

fn finish(sh: Shape, p: P4, chi2: f64, converged: bool, xs: &[f64], sg: &[f64]) -> Option<ResFit> {
    let mut cov = covariance(sh, &p, xs, sg)?;
    let dof = xs.len() as i64 - 4;
    let scale = if dof > 0 && chi2 / dof as f64 > 1.0 { (chi2 / dof as f64).sqrt() } else { 1.0 };
    for r in 0..4 { for c in 0..4 { cov[r][c] *= scale * scale; } }
    let e = [cov[0][0].max(0.0).sqrt(), cov[1][1].max(0.0).sqrt(), cov[2][2].max(0.0).sqrt(), cov[3][3].max(0.0).sqrt()];
    let (r, r_err) = r_and_err(&p, &cov);
    Some(ResFit { shape: sh, p, e, cov, chi2, chi2_flat: f64::NAN, dof, scale, n_pts: xs.len(), r, r_err, r_boot: None, f0_boot: None, converged })
}

fn bin_var(p: f64, ni: u64) -> f64 {
    // Jeffreys floor keeps sigma finite where the model predicts ~0 or ~1
    let ni = ni as f64;
    let pj = p.clamp(0.0, 1.0).max(0.5 / (ni + 1.0)).min(1.0 - 0.5 / (ni + 1.0));
    pj * (1.0 - pj) / ni
}

/// Fit one shape to a binomial resonance curve (k successes out of n at each frequency), model-based variances iterated
/// three times. `boot` parametric replicas give a bootstrap error on R and f0.
pub fn fit_binomial_shape(sh: Shape, xs: &[f64], ks: &[u64], ns: &[u64], boot: usize, seed: u64) -> Option<ResFit> {
    let keep: Vec<usize> = (0..xs.len()).filter(|&i| ns[i] > 0).collect();
    if keep.len() < 6 { return None; }
    let x: Vec<f64> = keep.iter().map(|&i| xs[i]).collect();
    let k: Vec<u64> = keep.iter().map(|&i| ks[i]).collect();
    let n: Vec<u64> = keep.iter().map(|&i| ns[i]).collect();
    let y: Vec<f64> = (0..x.len()).map(|i| k[i] as f64 / n[i] as f64).collect();
    // initial sigma from the data (Jeffreys-shifted), then model-based
    let mut sg: Vec<f64> = (0..x.len()).map(|i| bin_var((k[i] as f64 + 0.5) / (n[i] as f64 + 1.0), n[i]).sqrt()).collect();
    let (mut p, mut chi2, mut conv) = fit_shape(sh, &x, &y, &sg, None)?;
    for _ in 0..3 {
        sg = (0..x.len()).map(|i| bin_var(model(sh, &p, x[i]), n[i]).sqrt()).collect();
        let Some(r) = fit_shape(sh, &x, &y, &sg, Some(p)) else { break };
        p = r.0; chi2 = r.1; conv = r.2;
    }
    sg = (0..x.len()).map(|i| bin_var(model(sh, &p, x[i]), n[i]).sqrt()).collect();
    let _ = chi2;
    let chi2 = chi2_of(sh, &p, &x, &y, &sg);
    let mut f = finish(sh, p, chi2, conv, &x, &sg)?;
    {
        let (ks, nn): (f64, f64) = (k.iter().sum::<u64>() as f64, n.iter().sum::<u64>() as f64);
        let pbar = ks / nn;
        f.chi2_flat = (0..x.len()).map(|i| (y[i] - pbar).powi(2) / bin_var(pbar, n[i])).sum();
    }
    if boot > 0 {
        let mut rng = Rng::new(seed ^ 0xB007);
        let (mut rs, mut fs) = (vec![], vec![]);
        let bnds = bounds_of(&x);
        for _ in 0..boot {
            let ys: Vec<f64> = (0..x.len()).map(|i| rng.binomial(n[i], model(sh, &f.p, x[i]).clamp(0.0, 1.0)) as f64 / n[i] as f64).collect();
            let sgb: Vec<f64> = (0..x.len()).map(|i| bin_var(model(sh, &f.p, x[i]), n[i]).sqrt()).collect();
            let (pb, _, _) = lm(sh, f.p, &x, &ys, &sgb, &bnds);
            rs.push(pb[2] / pb[3]); fs.push(pb[2]);
        }
        let sd = |v: &[f64]| { let m = v.iter().sum::<f64>() / v.len() as f64; (m, (v.iter().map(|r| (r - m).powi(2)).sum::<f64>() / (v.len() as f64 - 1.0)).sqrt()) };
        f.r_boot = Some(sd(&rs)); f.f0_boot = Some(sd(&fs).1);
    }
    Some(f)
}

/// Fit both shapes; the lower chi2 wins.
pub fn fit_binomial_curve(xs: &[f64], ks: &[u64], ns: &[u64], boot: usize, seed: u64) -> Option<ResFit> {
    let mut best: Option<ResFit> = None;
    for sh in [Shape::Lorentz, Shape::Gauss] {
        if let Some(f) = fit_binomial_shape(sh, xs, ks, ns, boot, seed) { if best.as_ref().map_or(true, |b| f.chi2 < b.chi2) { best = Some(f); } }
    }
    best
}

/// A fitted resonance is reported only if it beats the flat curve by this much in chi2 (about 4 sigma for 3 extra parameters)
/// and is well-posed (centre inside the scanned range, finite errors, width between half a step and the range).
pub const MIN_DELTA_CHI2: f64 = 16.0;

pub fn is_significant(f: &ResFit, xs: &[f64]) -> bool {
    let (lo, hi) = (xs.iter().cloned().fold(f64::INFINITY, f64::min), xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
    let step = (hi - lo) / (xs.len().max(2) - 1) as f64;
    f.e.iter().all(|e| e.is_finite()) && f.delta_chi2() >= MIN_DELTA_CHI2 && f.p[2] >= lo && f.p[2] <= hi && f.e[2] < 0.5 * (hi - lo) && f.p[3] >= 0.5 * step && f.p[3] <= hi - lo && f.r_err.is_finite()
}

/// FWHM -> resolution in m/z terms: Delta(m/z) = FWHM_f / |df/d(m/z)|, R_mz = (m/z) / Delta(m/z); errors from the FWHM error.
pub fn r_mz(mz: f64, dfdmz: f64, fwhm: f64, fwhm_err: f64) -> (f64, f64, f64, f64) {
    let slope = dfdmz.abs().max(1e-12);
    let (d, de) = (fwhm / slope, fwhm_err / slope);
    (d, de, mz / d, mz / (d * d) * de)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Synthetic binomial Lorentzian scan (as the Sweep produces it).
    fn synth(sh: Shape, truth: P4, rng: &mut Rng, nper: u64) -> (Vec<f64>, Vec<u64>, Vec<u64>) {
        let xs: Vec<f64> = (0..=30).map(|i| 70.0 + i as f64).collect();
        let ns = vec![nper; xs.len()];
        let ks = xs.iter().map(|&x| rng.binomial(nper, model(sh, &truth, x).clamp(0.0, 1.0))).collect();
        (xs, ks, ns)
    }

    #[test]
    fn recovers_noiseless_lorentzian_exactly() {
        let truth: P4 = [0.01, 0.2, 83.0, 9.0];
        let xs: Vec<f64> = (0..=30).map(|i| 70.0 + i as f64).collect();
        let ys: Vec<f64> = xs.iter().map(|&x| model(Shape::Lorentz, &truth, x)).collect();
        let sg = vec![0.01; xs.len()];
        let f = fit_shape_sigma(Shape::Lorentz, &xs, &ys, &sg).unwrap();
        for k in 0..4 { assert!((f.p[k] - truth[k]).abs() < 1e-6 * truth[k].abs().max(1.0), "{:?} vs {:?}", f.p, truth); }
        assert!(f.chi2 < 1e-9);
    }

    #[test]
    fn lorentzian_parameters_and_uncertainties_cover() {
        // 400 synthetic scans (31 points x 700 ions, R = 83/9 = 9.2): the +-1 sigma intervals of every parameter and of R
        // must cover the truth ~68 % of the time, and the bias must be small compared to the error.
        let truth: P4 = [0.01, 0.2, 83.0, 9.0];
        let r_true = truth[2] / truth[3];
        let mut rng = Rng::new(2024);
        let reps = 400;
        let mut hit = [0usize; 5];
        let (mut sum_dev, mut sum_err) = ([0.0f64; 4], [0.0f64; 4]);
        let mut used = 0;
        for _ in 0..reps {
            let (xs, ks, ns) = synth(Shape::Lorentz, truth, &mut rng, 700);
            used += 1;
            let Some(fl) = fit_binomial_shape(Shape::Lorentz, &xs, &ks, &ns, 0, 1) else { continue };
            for k in 0..4 { if (fl.p[k] - truth[k]).abs() <= fl.e[k] { hit[k] += 1; } sum_dev[k] += fl.p[k] - truth[k]; sum_err[k] += fl.e[k]; }
            if (fl.r - r_true).abs() <= fl.r_err { hit[4] += 1; }
        }
        let cov: Vec<f64> = hit.iter().map(|&h| h as f64 / used as f64).collect();
        println!("coverage b/A/f0/w/R = {:.3?}; mean bias/err = {:.2?}", cov, (0..4).map(|k| sum_dev[k] / sum_err[k]).collect::<Vec<_>>());
        for (i, c) in cov.iter().enumerate() { assert!((c - 0.6827).abs() < 0.075, "param {i}: coverage {c}"); }
        for k in 0..4 { assert!((sum_dev[k] / sum_err[k]).abs() < 0.25, "bias of param {k}"); }
    }

    #[test]
    fn gaussian_is_chosen_for_gaussian_data_and_recovers_width() {
        let truth: P4 = [0.0, 0.3, 83.0, 8.0];
        let mut rng = Rng::new(77);
        let mut gauss_wins = 0;
        let mut w_ok = 0;
        let reps = 60;
        for _ in 0..reps {
            let (xs, ks, ns) = synth(Shape::Gauss, truth, &mut rng, 1500);
            let f = fit_binomial_curve(&xs, &ks, &ns, 0, 1).unwrap();
            if f.shape == Shape::Gauss { gauss_wins += 1; if (f.p[3] - truth[3]).abs() < 4.0 * f.e[3] { w_ok += 1; } }
        }
        assert!(gauss_wins as f64 / reps as f64 > 0.8, "gauss chosen {gauss_wins}/{reps}");
        assert!(w_ok as f64 / gauss_wins as f64 > 0.95);
    }

    #[test]
    fn bootstrap_error_on_r_agrees_with_covariance() {
        let truth: P4 = [0.01, 0.2, 83.0, 9.0];
        let mut rng = Rng::new(5);
        let (xs, ks, ns) = synth(Shape::Lorentz, truth, &mut rng, 700);
        let f = fit_binomial_curve(&xs, &ks, &ns, 200, 3).unwrap();
        let (_, sd) = f.r_boot.unwrap();
        println!("R = {:.3} ± {:.3} (cov)  boot sd {:.3}", f.r, f.r_err, sd);
        assert!((sd / f.r_err - 1.0).abs() < 0.35, "boot {sd} vs cov {}", f.r_err);
    }

    #[test]
    fn a_flat_scan_is_not_called_a_resonance_and_a_real_one_is() {
        let mut rng = Rng::new(404);
        let xs: Vec<f64> = (0..=30).map(|i| 35.0 + i as f64).collect();
        let mut false_pos = 0;
        let reps = 300;
        for _ in 0..reps {
            let ks: Vec<u64> = xs.iter().map(|_| rng.binomial(1000, 0.004)).collect();
            if let Some(f) = fit_binomial_curve(&xs, &ks, &vec![1000; xs.len()], 0, 1) { if is_significant(&f, &xs) { false_pos += 1; } }
        }
        println!("false resonances on flat scans: {false_pos}/{reps}");
        assert!(false_pos as f64 / (reps as f64) < 0.03, "{false_pos}");
        // a weak real one (peak 4 % on a 0.4 % floor, FWHM 1.5 kHz, 1 kHz steps, 1000 ions per point)
        let mut found = 0;
        for _ in 0..100 {
            let ks: Vec<u64> = xs.iter().map(|&x| rng.binomial(1000, model(Shape::Lorentz, &[0.004, 0.04, 52.7, 1.5], x).clamp(0.0, 1.0))).collect();
            let f = fit_binomial_curve(&xs, &ks, &vec![1000; xs.len()], 0, 1).unwrap();
            if is_significant(&f, &xs) && (f.p[2] - 52.7).abs() < 0.5 { found += 1; }
        }
        assert!(found >= 95, "found {found}/100");
    }

    #[test]
    fn r_in_mz_terms() {
        // f = K/(m/z): df/dmz = -f/mz, so R_mz = R_f
        let (mz, f0, fw) = (863.5, 83.0, 9.0);
        let (_, _, r, _) = r_mz(mz, -f0 / mz, fw, 0.5);
        assert!((r - f0 / fw).abs() < 1e-9);
    }
}
