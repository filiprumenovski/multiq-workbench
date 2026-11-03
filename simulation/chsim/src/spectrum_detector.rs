//! Detector view of the spectrum forward model (see chsim/SPECTRUM.md).
//!
//! species -> isotope envelope -> Orbitrap line profile -> list of resolved "lines" (local maxima of the profile)
//! -> pen kernel (how the figure's polyline draws a line on the pixel grid) -> apex height per pixel column.
//!
//! Everything here is CPU-only and deterministic.

pub const PROTON: f64 = 1.007276;
/// mean mass increment of one heavy-isotope step (13C-12C = 1.003355, 2H-1H = 1.006277, 15N-14N = 0.997035, 17/18O ...)
pub const NEUTRON_STEP: f64 = 1.002866;
pub const M_UBQ_AVG: f64 = 8564.8;
const KMAX: usize = 48;

// ---------------------------------------------------------------- isotope envelopes

#[derive(Clone, Debug)]
pub struct Envelope { pub masses: Vec<f64>, pub p: Vec<f64> } // neutral masses, probabilities (sum 1)

#[derive(Clone, Copy)]
struct Bin { p: f64, m: f64 }

type Dist = Vec<Bin>;

fn conv(a: &Dist, b: &Dist) -> Dist {
    let n = (a.len() + b.len() - 1).min(KMAX);
    let mut out = vec![Bin { p: 0.0, m: 0.0 }; n];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            if i + j >= n { continue; }
            let p = x.p * y.p;
            out[i + j].p += p;
            out[i + j].m += p * (x.m + y.m); // probability-weighted mass, divided below
        }
    }
    for o in out.iter_mut() { if o.p > 0.0 { o.m /= o.p; } }
    out
}

fn power(base: &Dist, mut n: u32) -> Dist {
    let mut result: Dist = vec![Bin { p: 1.0, m: 0.0 }];
    let mut b = base.clone();
    while n > 0 {
        if n & 1 == 1 { result = conv(&result, &b); }
        n >>= 1;
        if n > 0 { b = conv(&b, &b); }
    }
    result
}

/// (mass, abundance) per isotope, grouped by nominal shift from the lightest isotope.
const C: &[(u32, f64, f64)] = &[(0, 12.0, 0.9893), (1, 13.003355, 0.0107)];
const H: &[(u32, f64, f64)] = &[(0, 1.007825, 0.999885), (1, 2.014102, 0.000115)];
const N: &[(u32, f64, f64)] = &[(0, 14.003074, 0.99636), (1, 15.000109, 0.00364)];
const O: &[(u32, f64, f64)] = &[(0, 15.994915, 0.99757), (1, 16.999132, 0.00038), (2, 17.999160, 0.00205)];
const S: &[(u32, f64, f64)] = &[(0, 31.972071, 0.9499), (1, 32.971459, 0.0075), (2, 33.967867, 0.0425), (4, 35.967081, 0.0001)];

fn atom(el: &[(u32, f64, f64)]) -> Dist {
    let kmax = el.iter().map(|e| e.0).max().unwrap() as usize;
    let mut d = vec![Bin { p: 0.0, m: 0.0 }; kmax + 1];
    for &(k, m, a) in el { d[k as usize] = Bin { p: a, m }; }
    d
}

/// Exact isotope distribution of an integer formula (binomial convolution; nominal-mass bins with mean masses).
pub fn envelope_formula(c: u32, h: u32, n: u32, o: u32, s: u32) -> Envelope {
    let mut d: Dist = vec![Bin { p: 1.0, m: 0.0 }];
    for (el, k) in [(C, c), (H, h), (N, n), (O, o), (S, s)] { d = conv(&d, &power(&atom(el), k)); }
    let tot: f64 = d.iter().map(|b| b.p).sum();
    let keep: Vec<&Bin> = d.iter().filter(|b| b.p / tot > 1e-5).collect();
    let kt: f64 = keep.iter().map(|b| b.p).sum();
    Envelope { masses: keep.iter().map(|b| b.m).collect(), p: keep.iter().map(|b| b.p / kt).collect() }
}

/// Ubiquitin: C378H629N105O118S1 (average mass 8564.7).
pub fn envelope_ubiquitin() -> Envelope { envelope_formula(378, 629, 105, 118, 1) }

/// Averagine envelope (compound-Poisson, real-valued atom counts) for an organic ion of neutral mass `mass`.
pub fn envelope_averagine(mass: f64) -> Envelope {
    let nav = mass / 111.1254;
    let comp: [(&[(u32, f64, f64)], f64); 5] = [(C, 4.9384), (H, 7.7583), (N, 1.3577), (O, 1.4773), (S, 0.0417)];
    let mut lam = [0.0f64; KMAX];
    for (el, per) in comp { for &(k, _, a) in el.iter().filter(|e| e.0 > 0) { lam[k as usize] += nav * per * a; } }
    let mut p = vec![0.0f64; 24];
    p[0] = (-lam.iter().sum::<f64>()).exp();
    for k in 1..p.len() {
        let mut s = 0.0;
        for sft in 1..=k.min(KMAX - 1) { s += sft as f64 * lam[sft] * p[k - sft]; }
        p[k] = s / k as f64;
    }
    let tot: f64 = p.iter().sum();
    let mean_shift: f64 = p.iter().enumerate().map(|(k, x)| k as f64 * x / tot).sum();
    let mono = mass - mean_shift * NEUTRON_STEP;
    let keep: Vec<usize> = (0..p.len()).filter(|&k| p[k] / tot > 1e-5).collect();
    let kt: f64 = keep.iter().map(|&k| p[k]).sum();
    Envelope { masses: keep.iter().map(|&k| mono + k as f64 * NEUTRON_STEP).collect(), p: keep.iter().map(|&k| p[k] / kt).collect() }
}

// ---------------------------------------------------------------- Orbitrap line profile

#[derive(Clone, Copy, Debug)]
pub struct Orbitrap { pub r200: f64 }
impl Orbitrap {
    /// Resolving power at m/z: R = R200 * sqrt(200 / (m/z)) (Orbitrap: frequency ~ 1/sqrt(m/z), transient length fixed).
    pub fn resolving_power(&self, mz: f64) -> f64 { self.r200 * (200.0 / mz).sqrt() }
    pub fn fwhm(&self, mz: f64) -> f64 { mz / self.resolving_power(mz) }
    pub fn sigma(&self, mz: f64) -> f64 { self.fwhm(mz) / 2.354820045 }
}

/// A species as the detector sees it: resolved lines relative to the apex line.
#[derive(Clone, Debug)]
pub struct Lines {
    /// (m/z offset from apex, relative height 0..1) of every local maximum of the profile
    pub lines: Vec<(f64, f64)>,
    /// apex height of the area-normalised profile (1/m/z units): detector response to unit abundance
    pub apex_gain: f64,
}

/// Profile of an envelope at charge z: sum_i p_i N(m; mz_i, sigma(mz_i)); returns the local maxima.
pub fn detector_lines(env: &Envelope, z: f64, det: &Orbitrap) -> Lines {
    let mzs: Vec<f64> = env.masses.iter().map(|&m| (m + z * PROTON) / z).collect();
    let (lo, hi) = (mzs[0], *mzs.last().unwrap());
    let sig = det.sigma(0.5 * (lo + hi));
    let step = (sig / 4.0).min(0.02).max(sig / 12.0);
    let (g0, g1) = (lo - 5.0 * sig, hi + 5.0 * sig);
    let ns = ((g1 - g0) / step).ceil() as usize + 1;
    let y: Vec<f64> = (0..ns).map(|i| {
        let m = g0 + i as f64 * step;
        env.p.iter().zip(&mzs).map(|(&p, &mu)| { let s = det.sigma(mu); p * (-(m - mu).powi(2) / (2.0 * s * s)).exp() / (s * (2.0 * std::f64::consts::PI).sqrt()) }).sum::<f64>()
    }).collect();
    let ymax = y.iter().cloned().fold(0.0, f64::max);
    let imax = y.iter().position(|&v| v == ymax).unwrap();
    let apex = g0 + imax as f64 * step;
    let mut lines = vec![];
    for i in 1..ns - 1 {
        if y[i] >= y[i - 1] && y[i] > y[i + 1] && y[i] / ymax > 0.004 {
            // parabolic refinement of the maximum position
            let (a, b, c) = (y[i - 1], y[i], y[i + 1]);
            let off = 0.5 * (a - c) / (a - 2.0 * b + c);
            lines.push((g0 + (i as f64 + off) * step - apex, b / ymax));
        }
    }
    Lines { lines, apex_gain: ymax }
}

// ---------------------------------------------------------------- pen kernel (how the figure draws a line)

/// How an isolated line is drawn on the pixel grid: `nbox` columns at full height, then tails (fractions of the
/// plateau height) on the low-m/z (`left`) and high-m/z (`right`) side. Index 0 = first column outside the box.
#[derive(Clone, Debug)]
pub struct Pen { pub nbox: usize, pub left: Vec<f64>, pub right: Vec<f64> }
impl Pen {
    pub fn bare(nbox: usize) -> Pen { Pen { nbox, left: vec![], right: vec![] } }
    /// First box column for a line at fractional column position x (columns are integers, box = [x - nbox/2, x + nbox/2)).
    pub fn p0(&self, x: f64) -> i64 { (x - 0.5 * self.nbox as f64).ceil() as i64 }
    /// (column offset from p0, kernel value) for all columns this pen touches.
    pub fn taps(&self) -> Vec<(i64, f64)> {
        let mut t = vec![];
        for (d, &v) in self.left.iter().enumerate() { t.push((-(d as i64) - 1, v)); }
        for k in 0..self.nbox { t.push((k as i64, 1.0)); }
        for (d, &v) in self.right.iter().enumerate() { t.push((self.nbox as i64 + d as i64, v)); }
        t
    }
}

// ---------------------------------------------------------------- panel geometry + rendering

/// Linear column <-> m/z map of one panel.
#[derive(Clone, Debug)]
pub struct Grid { pub col0: i64, pub mz0: f64, pub dmz: f64, pub n: usize }
impl Grid {
    pub fn x_of(&self, mz: f64) -> f64 { self.col0 as f64 + (mz - self.mz0) / self.dmz }
    pub fn mz_of(&self, i: usize) -> f64 { self.mz0 + i as f64 * self.dmz }
}

/// A species placed on a panel: `apex_mz` is where its apex line sits.
#[derive(Clone, Debug)]
pub struct Placed { pub apex_mz: f64, pub lines: Vec<(f64, f64)> }
impl Placed {
    /// Apex line only. The 5A/5B/5C pixels are 3 m/z wide and drawn with a ~4 px pen, and the data show the pen plateau
    /// follows the envelope apex (no partial-height shoulder columns from satellite isotope lines): collapse the envelope.
    pub fn collapsed(&self) -> Placed { Placed { apex_mz: self.apex_mz, lines: vec![(0.0, 1.0)] } }
}

/// Per-species unit-apex response over panel columns (sparse), for a given calibration shift (m/z).
#[derive(Clone, Debug)]
pub struct Response { pub first: usize, pub v: Vec<f64> }

pub fn response(sp: &Placed, grid: &Grid, pen: &Pen, shift_mz: f64) -> Response {
    let taps = pen.taps();
    // best plateau per box start column, then convolve with the taps (max-plus)
    let mut acc: std::collections::BTreeMap<i64, f64> = std::collections::BTreeMap::new();
    for &(dm, y) in &sp.lines {
        let x = grid.x_of(sp.apex_mz + dm + shift_mz);
        let p0 = pen.p0(x);
        let e = acc.entry(p0).or_insert(0.0);
        if y > *e { *e = y; }
    }
    let mut cols: std::collections::BTreeMap<i64, f64> = std::collections::BTreeMap::new();
    for (&p0, &y) in &acc {
        for &(off, k) in &taps {
            let c = p0 + off;
            let e = cols.entry(c).or_insert(0.0);
            if y * k > *e { *e = y * k; }
        }
    }
    let (c_lo, c_hi) = match (cols.keys().next(), cols.keys().next_back()) { (Some(&a), Some(&b)) => (a, b), _ => return Response { first: 0, v: vec![] } };
    let (c_lo, c_hi) = ((c_lo - grid.col0).max(0), (c_hi - grid.col0).min(grid.n as i64 - 1));
    if c_hi < c_lo { return Response { first: 0, v: vec![] }; }
    Response { first: c_lo as usize, v: (c_lo..=c_hi).map(|c| *cols.get(&(c + grid.col0)).unwrap_or(&0.0)).collect() }
}

/// Max-plus combination of species responses: returns (height, kernel factor of the dominant term, relative error of the
/// dominant species, index of dominant species or usize::MAX).
pub fn combine(resp: &[Response], amp: &[f64], rel: &[f64], n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<usize>) {
    let mut h = vec![0.0; n];
    let mut kf = vec![0.0; n];
    let mut re = vec![0.0; n];
    let mut who = vec![usize::MAX; n];
    for (j, r) in resp.iter().enumerate() {
        let a = amp[j];
        if a <= 0.0 { continue; }
        for (i, &v) in r.v.iter().enumerate() {
            let c = r.first + i;
            let val = a * v;
            if val > h[c] { h[c] = val; kf[c] = v; re[c] = rel[j]; who[c] = j; }
        }
    }
    (h, kf, re, who)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn ubiquitin_envelope_mass_and_normalisation() {
        let e = envelope_ubiquitin();
        let tot: f64 = e.p.iter().sum();
        assert!((tot - 1.0).abs() < 1e-9);
        let avg: f64 = e.masses.iter().zip(&e.p).map(|(m, p)| m * p).sum();
        assert!((avg - 8564.8).abs() < 0.3, "average mass {}", avg);
        assert!((e.masses[0] - 8559.617).abs() < 0.02, "mono {}", e.masses[0]);
        // most-abundant isotope is M+4..M+5 for ubiquitin
        let imax = e.p.iter().cloned().enumerate().fold((0, 0.0), |a, (i, p)| if p > a.1 { (i, p) } else { a }).0;
        assert!(imax == 4 || imax == 5, "imax {}", imax);
    }
    #[test]
    fn averagine_mean_mass_is_preserved() {
        let e = envelope_averagine(402.0);
        let avg: f64 = e.masses.iter().zip(&e.p).map(|(m, p)| m * p).sum::<f64>() / e.p.iter().sum::<f64>();
        assert!((avg - 402.0).abs() < 0.05, "avg {}", avg);
        assert!(e.p[1] > 0.1 && e.p[1] < 0.35); // ~ 17 C atoms -> M+1 about 20-25 % of M
    }
    #[test]
    fn resolved_lines_scale_with_resolution() {
        let env = envelope_ubiquitin();
        let hi = detector_lines(&env, 8.0, &Orbitrap { r200: 140_000.0 });
        let lo = detector_lines(&env, 8.0, &Orbitrap { r200: 3_000.0 });
        assert!(hi.lines.len() > 8, "high R should resolve the isotopes: {}", hi.lines.len());
        assert!(lo.lines.len() <= 2, "low R merges them: {}", lo.lines.len());
        assert!(hi.apex_gain > lo.apex_gain);
    }
    #[test]
    fn pen_box_covers_nbox_columns() {
        let pen = Pen { nbox: 4, left: vec![0.1, 0.05], right: vec![0.3, 0.1, 0.05] };
        let g = Grid { col0: 0, mz0: 0.0, dmz: 3.0, n: 100 };
        for k in 0..20 {
            let mz = 150.0 + 0.7 * k as f64;
            let r = response(&Placed { apex_mz: mz, lines: vec![(0.0, 1.0)] }, &g, &pen, 0.0);
            assert_eq!(r.v.iter().filter(|&&v| v == 1.0).count(), 4);
        }
    }
}
