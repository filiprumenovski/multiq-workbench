//! Counting statistics for every number mqviz shows: binomial fractions, Poisson rates, ratios, one-sided limits and the
//! convergence ("settled" / "accumulating") test that compares two consecutive windows.
//!
//! Conventions (all 1-sigma unless a name says 95):
//! * a fraction k/n is shown as k/n +- sigma_W, the Wilson-score half-width at z = 1,
//!   sigma_W = sqrt(p(1-p)/n + 1/(4 n^2)) / (1 + 1/n), which equals the Wald error at large n and stays non-zero at k = 0 or k = n;
//! * a Poisson count N has sigma = sqrt(max(N, 1)) (a count of zero is quoted as 0 +- 1);
//! * a ratio of two fractions uses the log-delta method with Jeffreys-shifted counts.

/// Binomial fraction: `k` successes out of `n` trials.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Frac { pub k: u64, pub n: u64 }

impl Frac {
    pub fn new(k: u64, n: u64) -> Frac { Frac { k, n } }
    pub fn p(&self) -> f64 { if self.n == 0 { f64::NAN } else { self.k as f64 / self.n as f64 } }
    /// Wilson-score 1-sigma error (see module docs).
    pub fn sigma(&self) -> f64 {
        if self.n == 0 { return f64::NAN; }
        let (n, p) = (self.n as f64, self.p());
        ((p * (1.0 - p) / n + 0.25 / (n * n)).sqrt()) / (1.0 + 1.0 / n)
    }
    pub fn add(&self, o: &Frac) -> Frac { Frac { k: self.k + o.k, n: self.n + o.n } }
}

/// Poisson count error.
pub fn poisson_sigma(n: u64) -> f64 { (n.max(1) as f64).sqrt() }

/// One-sided 95 % upper limit of a Poisson mean given `k` observed counts (Wilson-Hilferty; 2.97 for k = 0, exact 3.00).
pub fn poisson_ul95(k: u64) -> f64 {
    let a = (k + 1) as f64;
    a * (1.0 - 1.0 / (9.0 * a) + 1.644_854 / (3.0 * a.sqrt())).powi(3)
}

/// A value with its 1-sigma error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Val { pub v: f64, pub e: f64 }
impl Val {
    pub fn new(v: f64, e: f64) -> Val { Val { v, e } }
    pub fn nan() -> Val { Val { v: f64::NAN, e: f64::NAN } }
    pub fn ok(&self) -> bool { self.v.is_finite() && self.e.is_finite() }
}

impl From<Frac> for Val { fn from(f: Frac) -> Val { Val { v: f.p(), e: f.sigma() } } }

/// Rate N/dur from a Poisson count.
pub fn rate(n: u64, dur: f64) -> Val { if dur > 0.0 { Val { v: n as f64 / dur, e: poisson_sigma(n) / dur } } else { Val::nan() } }

/// Ratio a/b of two binomial fractions (independent). Log-delta error with Jeffreys-shifted counts; if the numerator has no
/// counts the value is 0 and the error is the one-count scale (a/b at k = 1), so a pull against it is meaningful.
pub fn ratio(a: Frac, b: Frac) -> Val {
    if a.n == 0 || b.n == 0 || b.k == 0 { return Val::nan(); }
    let (pa, pb) = (a.p(), b.p());
    let vb = (1.0 - pb) / (b.k as f64 + 0.5);
    if a.k == 0 {
        let one = (1.0 / a.n as f64) / pb;
        return Val { v: 0.0, e: one * (1.0 + vb).sqrt() };
    }
    let va = (1.0 - pa) / (a.k as f64 + 0.5);
    let r = pa / pb;
    Val { v: r, e: r * (va + vb).sqrt() }
}

/// 95 % upper limit of a/b (Poisson limit on the numerator, denominator error added in quadrature in log space).
pub fn ratio_ul95(a: Frac, b: Frac) -> f64 {
    if a.n == 0 || b.n == 0 || b.k == 0 { return f64::NAN; }
    let pb = b.p();
    (poisson_ul95(a.k) / a.n as f64) / pb * (1.0 + 1.645 * ((1.0 - pb) / (b.k as f64 + 0.5)).sqrt())
}

/// Pull of `ours` against `theirs` in sigma, errors in quadrature.
pub fn pull(ours: Val, theirs: Val) -> f64 { (ours.v - theirs.v) / (ours.e * ours.e + theirs.e * theirs.e).sqrt().max(1e-300) }

/// Mark like chsim reality: check (|pull| <= 1), tilde (<= 2), cross.
pub fn mark(p: f64) -> &'static str { if !p.is_finite() { "…" } else if p.abs() <= 1.0 { "✓" } else if p.abs() <= 2.0 { "~" } else { "✗" } }

// ---------------------------------------------------------------- convergence

/// State of an estimate that is being accumulated.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Settle {
    /// Two consecutive windows agree (|z| <= 2) and the estimate is precise enough; `n` ions/counts rest behind it.
    Settled { n: u64 },
    /// Not enough counts yet for the precision target.
    Accumulating { n: u64 },
    /// Enough counts in both windows but they disagree beyond 2 sigma (the system is still relaxing or the drive changed).
    Drifting { n: u64, z: f64 },
}
impl Settle {
    pub fn n(&self) -> u64 { match *self { Settle::Settled { n } | Settle::Accumulating { n } | Settle::Drifting { n, .. } => n } }
    pub fn is_settled(&self) -> bool { matches!(self, Settle::Settled { .. }) }
    pub fn label(&self) -> String {
        match *self {
            Settle::Settled { n } => format!("settled · {} ions", short_count(n)),
            Settle::Accumulating { n } => format!("accumulating · {} ions", short_count(n)),
            Settle::Drifting { n, z } => format!("drifting {:+.1}σ · {} ions", z, short_count(n)),
        }
    }
}

pub fn short_count(n: u64) -> String { if n >= 10_000_000 { format!("{:.1}M", n as f64 / 1e6) } else if n >= 10_000 { format!("{:.0}k", n as f64 / 1e3) } else if n >= 1000 { format!("{:.1}k", n as f64 / 1e3) } else { format!("{}", n) } }

/// Precision targets: a fraction is precise enough when sigma/p <= `REL_TARGET` and at least `k_min` successes stand behind it;
/// each window needs `N_MIN` trials.
pub const N_MIN: u64 = 100;
pub const REL_TARGET: f64 = 0.25;
/// |z| above which two consecutive windows are called different.
pub const Z_DRIFT: f64 = 2.0;

/// Convergence of a binomial fraction from its two consecutive windows `a` (earlier) and `b` (later).
/// `k_min` is the number of successes below which a fraction is never called settled (rare-event estimates).
pub fn settle_frac(a: Frac, b: Frac, k_min: u64) -> Settle {
    let all = a.add(&b);
    let n = all.n;
    if a.n < N_MIN || b.n < N_MIN { return Settle::Accumulating { n }; }
    let z = (a.p() - b.p()) / (a.sigma().powi(2) + b.sigma().powi(2)).sqrt().max(1e-300);
    if z.abs() > Z_DRIFT && a.k + b.k >= 10 { return Settle::Drifting { n, z }; }
    if all.k >= k_min && all.sigma() <= REL_TARGET * all.p() { Settle::Settled { n } } else { Settle::Accumulating { n } }
}

/// Convergence of a rate from the counts of two consecutive equal-duration windows.
pub fn settle_rate(a: u64, b: u64) -> Settle {
    let n = a + b;
    if a < 30 || b < 30 { return Settle::Accumulating { n }; }
    let z = (a as f64 - b as f64) / ((a + b) as f64).sqrt();
    if z.abs() > Z_DRIFT { Settle::Drifting { n, z } } else { Settle::Settled { n } }
}

// ---------------------------------------------------------------- display

/// "value ± error" with the error to two significant digits.
pub fn pm(v: f64, e: f64) -> String {
    if !v.is_finite() { return "—".into(); }
    if !e.is_finite() || e <= 0.0 { return format!("{:.3}", v); }
    let d = ((1 - e.log10().floor() as i32).max(0) as usize).min(9);
    format!("{:.*} ± {:.*}", d, v, d, e)
}

/// Like `pm`, switching to a shared exponent "(2.15 ± 0.09)e6" for large or tiny magnitudes.
pub fn pm_auto(v: f64, e: f64) -> String {
    if !v.is_finite() { return "—".into(); }
    if v != 0.0 && (v.abs() >= 1e5 || v.abs() < 1e-3) && e.is_finite() && e > 0.0 {
        let x = v.abs().log10().floor() as i32;
        let (m, me) = (v / 10f64.powi(x), e / 10f64.powi(x));
        return format!("({})e{}", pm(m, me), x);
    }
    pm(v, e)
}

// ---------------------------------------------------------------- small deterministic RNG (tests, bootstrap)

pub struct Rng(pub u64);
impl Rng {
    pub fn new(seed: u64) -> Rng { Rng(seed.wrapping_mul(0x9E3779B97F4A7C15) | 1) }
    pub fn u(&mut self) -> f64 { self.0 ^= self.0 << 13; self.0 ^= self.0 >> 7; self.0 ^= self.0 << 17; (self.0 >> 11) as f64 / (1u64 << 53) as f64 }
    pub fn normal(&mut self) -> f64 { let (a, b) = (self.u().max(1e-300), self.u()); (-2.0 * a.ln()).sqrt() * (std::f64::consts::TAU * b).cos() }
    pub fn binomial(&mut self, n: u64, p: f64) -> u64 {
        if n > 400 && p * (1.0 - p) * n as f64 > 50.0 { // normal approximation is exact enough for the coverage checks at large n
            let x = n as f64 * p + (n as f64 * p * (1.0 - p)).sqrt() * self.normal();
            return x.round().clamp(0.0, n as f64) as u64;
        }
        (0..n).filter(|_| self.u() < p).count() as u64
    }
    #[allow(dead_code)]
    pub fn poisson(&mut self, lam: f64) -> u64 {
        if lam > 60.0 { return (lam + lam.sqrt() * self.normal()).round().max(0.0) as u64; }
        let l = (-lam).exp();
        let (mut k, mut p) = (0u64, 1.0);
        loop { p *= self.u(); if p <= l { return k; } k += 1; }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Coverage of value +- m sigma over many synthetic experiments.
    fn coverage_frac(n: u64, p: f64, reps: usize, m: f64, seed: u64) -> f64 {
        let mut rng = Rng::new(seed);
        let hit = (0..reps).filter(|_| { let f = Frac::new(rng.binomial(n, p), n); (f.p() - p).abs() <= m * f.sigma() }).count();
        hit as f64 / reps as f64
    }

    #[test]
    fn binomial_one_sigma_coverage() {
        // Averaged over the discrete sawtooth of binomial coverage (several n at each p), 1-sigma bands cover 68 +- 4 %,
        // from p = 0.5 down to expected counts of ~10 (the regime in which fractions are shown as measurements).
        let mut worst: f64 = 0.0;
        for &p in &[0.5, 0.2, 0.05, 0.01] {
            let ns = [(10.0 / p) as u64, (14.0 / p) as u64, (20.0 / p) as u64, (30.0 / p) as u64, (50.0 / p) as u64, (100.0 / p) as u64];
            let c1: f64 = ns.iter().enumerate().map(|(i, &n)| coverage_frac(n, p, 4000, 1.0, 11 + i as u64)).sum::<f64>() / ns.len() as f64;
            let c2: f64 = ns.iter().enumerate().map(|(i, &n)| coverage_frac(n, p, 4000, 2.0, 77 + i as u64)).sum::<f64>() / ns.len() as f64;
            println!("p {p}: coverage 1σ {c1:.3}  2σ {c2:.3}");
            assert!((c1 - 0.6827).abs() < 0.04, "1σ coverage {c1} at p={p}");
            assert!((c2 - 0.9545).abs() < 0.03, "2σ coverage {c2} at p={p}");
            worst = worst.max((c1 - 0.6827).abs());
        }
        println!("worst 1σ deviation {worst:.3}");
    }

    #[test]
    fn binomial_error_is_nonzero_at_the_boundary() {
        assert!(Frac::new(0, 50).sigma() > 0.0 && Frac::new(50, 50).sigma() > 0.0);
        // large n: equals the Wald error
        let f = Frac::new(215_000, 1_000_000);
        assert!((f.sigma() - (0.215f64 * 0.785 / 1e6).sqrt()).abs() < 1e-6);
    }

    #[test]
    fn poisson_rate_coverage() {
        let mut rng = Rng::new(5);
        for &lam in &[10.0, 30.0, 100.0, 1000.0] {
            let reps = 20000;
            let hit = (0..reps).filter(|_| { let k = rng.poisson(lam); (k as f64 - lam).abs() <= poisson_sigma(k) }).count() as f64 / reps as f64;
            println!("λ {lam}: coverage {hit:.3}");
            assert!((hit - 0.6827).abs() < 0.06, "λ={lam} cov {hit}");
        }
    }

    #[test]
    fn poisson_upper_limit_is_95_percent() {
        for k in [0u64, 1, 2, 5, 10, 30, 100] {
            // P(N <= k | mu = UL) must be ~ 5 %
            let ul = poisson_ul95(k);
            let mut cdf = 0.0; let mut term = (-ul).exp();
            for i in 0..=k { if i > 0 { term *= ul / i as f64; } cdf += term; }
            assert!((cdf - 0.05).abs() < 0.012, "k {k}: UL {ul} has P(N<=k) = {cdf}");
        }
    }

    #[test]
    fn ratio_coverage_binomial() {
        // T_a / T_b with independent binomial counts, moderate counts
        let mut rng = Rng::new(21);
        let (na, nb, pa, pb) = (20_000u64, 8_000u64, 0.02, 0.2);
        let reps = 6000;
        let hit = (0..reps).filter(|_| { let r = ratio(Frac::new(rng.binomial(na, pa), na), Frac::new(rng.binomial(nb, pb), nb)); (r.v - pa / pb).abs() <= r.e }).count() as f64 / reps as f64;
        println!("ratio coverage {hit:.3}");
        assert!((hit - 0.6827).abs() < 0.05, "{hit}");
    }

    #[test]
    fn plus_minus_formatting() {
        assert_eq!(pm(21.4321, 0.62), "21.43 ± 0.62");
        assert_eq!(pm(0.8549, 0.0157), "0.855 ± 0.016");
        assert_eq!(pm(1234.0, 56.0), "1234 ± 56");
        assert_eq!(pm_auto(2.1523e6, 0.09e6), "(2.152 ± 0.090)e6");
    }

    #[test]
    fn settle_detects_stationary_and_drifting() {
        let mut rng = Rng::new(3);
        let n = 4000u64;
        let a = Frac::new(rng.binomial(n, 0.86), n);
        let b = Frac::new(rng.binomial(n, 0.86), n);
        assert!(settle_frac(a, b, 10).is_settled());
        let c = Frac::new(rng.binomial(n, 0.70), n);
        assert!(matches!(settle_frac(a, c, 10), Settle::Drifting { .. }));
        assert!(matches!(settle_frac(Frac::new(3, 50), Frac::new(4, 50), 10), Settle::Accumulating { .. }));
        assert!(settle_rate(1000, 1030).is_settled());
        assert!(matches!(settle_rate(1000, 1400), Settle::Drifting { .. }));
        // a rare fraction with a handful of counts is never "settled"
        assert!(matches!(settle_frac(Frac::new(2, 3000), Frac::new(1, 3000), 10), Settle::Accumulating { .. }));
    }

    #[test]
    fn false_drift_rate_is_about_five_percent() {
        let mut rng = Rng::new(9);
        let reps = 4000;
        let drift = (0..reps).filter(|_| matches!(settle_frac(Frac::new(rng.binomial(2000, 0.3), 2000), Frac::new(rng.binomial(2000, 0.3), 2000), 10), Settle::Drifting { .. })).count() as f64 / reps as f64;
        println!("false drift rate {drift:.3}");
        assert!(drift < 0.08, "{drift}");
    }
}
