//! Likelihood side of the spectrum forward model: panels (pixel data + display floor), noise model (Gaussian with
//! Tobit censoring below the display floor), calibration-shift nuisance (marginalised on a grid), kernel fit and
//! source (abundance) fit. CPU only.

use crate::spectrum_detector::*;

// ---------------------------------------------------------------- numerics

/// Standard normal CDF (Hart 1968 / West 2005, double precision).
pub fn ncdf(x: f64) -> f64 {
    let xa = x.abs();
    let c = if xa > 37.0 { 0.0 } else {
        let e = (-0.5 * xa * xa).exp();
        if xa < 7.07106781186547 {
            let mut b = 3.52624965998911e-02 * xa + 0.700383064443688;
            b = b * xa + 6.37396220353165; b = b * xa + 33.912866078383; b = b * xa + 112.079291497871;
            b = b * xa + 221.213596169931; b = b * xa + 220.206867912376;
            let mut c = 8.83883476483184e-02 * xa + 1.75566716318264;
            c = c * xa + 16.064177579207; c = c * xa + 86.7807322029461; c = c * xa + 296.564248779674;
            c = c * xa + 637.333633378831; c = c * xa + 793.826512519948; c = c * xa + 440.413735824752;
            e * b / c
        } else {
            let mut b = xa + 0.65; b = xa + 4.0 / b; b = xa + 3.0 / b; b = xa + 2.0 / b; b = xa + 1.0 / b;
            e / b / 2.506628274631
        }
    };
    if x > 0.0 { 1.0 - c } else { c }
}
pub fn log_ncdf(x: f64) -> f64 {
    if x < -30.0 { -0.5 * x * x - (-x * 2.506628274631).ln() } else { ncdf(x).max(1e-300).ln() }
}
pub fn logsumexp(v: &[f64]) -> f64 {
    let m = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if !m.is_finite() { return m; }
    m + v.iter().map(|x| (x - m).exp()).sum::<f64>().ln()
}

// ---------------------------------------------------------------- panels

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PenKind { Fig5, S14 }

#[derive(Clone, Debug)]
pub struct Panel {
    pub name: String,
    pub grid: Grid,
    pub h: Vec<f64>,
    pub floor: f64,    // display floor, % of base peak: data below it are censored
    pub pxh: f64,      // height of one pixel, % of base peak
    pub sig_cal: f64,  // 1-sigma m/z calibration uncertainty
    pub base_obs: f64, // observed plateau of the base peak
    pub pen: PenKind,
}

impl Panel {
    pub fn mz(&self, i: usize) -> f64 { self.grid.mz_of(i) }
    pub fn from_rows(name: &str, cols: Vec<i64>, mz: Vec<f64>, h: Vec<f64>, floor: f64, pxh: f64, sig_cal: f64, pen: PenKind) -> Panel {
        let n = h.len();
        let dmz = (mz[n - 1] - mz[0]) / (cols[n - 1] - cols[0]) as f64;
        let nb = if pen == PenKind::Fig5 { 4 } else { 1 };
        let mut top = h.clone();
        top.sort_by(|a, b| b.partial_cmp(a).unwrap());
        let base_obs = top.iter().take(nb).sum::<f64>() / nb as f64;
        Panel { name: name.into(), grid: Grid { col0: cols[0], mz0: mz[0], dmz, n }, h, floor, pxh, sig_cal, base_obs, pen }
    }
    /// Same grid / floors, new heights (synthetic panels).
    pub fn with_heights(&self, h: Vec<f64>) -> Panel {
        let mut p = self.clone();
        let nb = if p.pen == PenKind::Fig5 { 4 } else { 1 };
        let mut top = h.clone(); top.sort_by(|a, b| b.partial_cmp(a).unwrap());
        p.base_obs = top.iter().take(nb).sum::<f64>() / nb as f64;
        p.h = h; p
    }
}

/// (name, floor %, pixel %, calibration sigma m/z) from SPECTRA_EXTRACTION.md sections 1-2 and 7.
pub fn panel_meta(name: &str) -> (f64, f64, f64) {
    match name {
        "A1" => (1.27, 0.741, 1.0), "A2" => (1.21, 0.735, 1.0), "A3" => (1.48, 0.733, 2.0), "A4" => (1.18, 0.735, 1.0),
        "B" => (0.52, 0.255, 0.3), "C" => (0.62, 0.265, 0.3),
        "S14I" => (0.50, 0.267, 0.35), "S14II" => (0.80, 0.278, 0.38),
        _ => (1.0, 0.5, 1.0),
    }
}

pub fn load_fig5(dir: &str) -> Result<Vec<Panel>, String> {
    let mut out = vec![];
    for (file, names) in [("fig5A_profiles_native.csv", vec!["A1", "A2", "A3", "A4"]), ("fig5BC_profiles_native.csv", vec!["B", "C"])] {
        let text = std::fs::read_to_string(format!("{}/{}", dir, file)).map_err(|e| format!("{}/{}: {}", dir, file, e))?;
        for nm in names {
            let (mut cols, mut mz, mut h) = (vec![], vec![], vec![]);
            for l in text.lines().filter(|l| !l.starts_with('#') && !l.starts_with("panel")) {
                let p: Vec<&str> = l.split(',').collect();
                if p.len() < 4 || p[0] != nm { continue; }
                cols.push(p[1].parse::<i64>().map_err(|e| e.to_string())?); mz.push(p[2].parse::<f64>().map_err(|e| e.to_string())?); h.push(p[3].parse::<f64>().map_err(|e| e.to_string())?);
            }
            if h.is_empty() { return Err(format!("panel {} not found in {}", nm, file)); }
            let (fl, px, sc) = panel_meta(nm);
            out.push(Panel::from_rows(nm, cols, mz, h, fl, px, sc, PenKind::Fig5));
        }
    }
    Ok(out)
}

pub fn load_s14(dir: &str) -> Result<Vec<Panel>, String> {
    let mut out = vec![];
    for (file, nm) in [("figS14_I_spectrum_profile.csv", "S14I"), ("figS14_II_spectrum_profile.csv", "S14II")] {
        let text = std::fs::read_to_string(format!("{}/{}", dir, file)).map_err(|e| format!("{}/{}: {}", dir, file, e))?;
        let (mut mz, mut h) = (vec![], vec![]);
        for l in text.lines().filter(|l| !l.starts_with('#') && !l.starts_with("mz")) {
            let p: Vec<&str> = l.split(',').collect();
            if p.len() < 2 { continue; }
            mz.push(p[0].parse::<f64>().map_err(|e| e.to_string())?); h.push(p[1].parse::<f64>().map_err(|e| e.to_string())?);
        }
        let cols: Vec<i64> = (0..h.len() as i64).collect();
        let (fl, px, sc) = panel_meta(nm);
        out.push(Panel::from_rows(nm, cols, mz, h, fl, px, sc, PenKind::S14));
    }
    Ok(out)
}

// ---------------------------------------------------------------- noise model and per-panel evaluation

#[derive(Clone, Copy, Debug)]
pub struct Noise {
    pub sigq_px: f64, // quantisation/extraction error in pixel heights (sigma = sigq_px * pixel height)
    pub rho_p: f64,   // relative error of plateau columns (extraction +-1.5 % of height)
    pub rho_t: f64,   // relative error of tail columns (pen-tail scatter between peaks is ~20 %)
}

#[derive(Clone, Debug, Default)]
pub struct ColDetail { pub mz: Vec<f64>, pub data: Vec<f64>, pub model: Vec<f64>, pub nll: Vec<f64>, pub kind: Vec<u8>, pub who: Vec<usize> } // kind 0 censored, 1 gaussian, 2 unmodelled (data above floor, model 0)

#[derive(Clone, Copy, Debug, Default)]
pub struct Split { pub total: f64, pub unmodelled: f64, pub n_unmod: usize, pub n_gauss: usize, pub n_cens: usize }

pub fn panel_nll(p: &Panel, resp: &[Response], amp: &[f64], rel: &[f64], normalise: bool, nz: &Noise, mut det: Option<&mut ColDetail>) -> Split {
    let (h, kf, re, who) = combine(resp, amp, rel, p.grid.n);
    let hmax = h.iter().cloned().fold(0.0, f64::max);
    let scale = if normalise && hmax > 0.0 { p.base_obs / hmax } else { 1.0 };
    let sq = nz.sigq_px * p.pxh;
    let mut s = Split::default();
    let l0 = -log_ncdf(p.floor / sq); // censored column with no model
    if let Some(d) = det.as_deref_mut() { *d = ColDetail::default(); }
    for c in 0..p.grid.n {
        let (dv, m) = (p.h[c], h[c] * scale);
        let l;
        let mut kind = 0u8;
        if m == 0.0 && dv < p.floor { l = l0; s.n_cens += 1; }
        else {
            let rho = if kf[c] >= 0.999 { nz.rho_p } else { nz.rho_t };
            let sg = (sq * sq + (rho * m).powi(2) + (re[c] * m).powi(2)).sqrt();
            if dv < p.floor { l = -log_ncdf((p.floor - m) / sg); s.n_cens += 1; }
            else {
                l = 0.5 * ((dv - m) / sg).powi(2) + sg.ln() + 0.918938533204673;
                if m == 0.0 { kind = 2; s.unmodelled += l; s.n_unmod += 1; } else { kind = 1; s.n_gauss += 1; }
            }
        }
        s.total += l;
        if let Some(d) = det.as_deref_mut() { d.mz.push(p.mz(c)); d.data.push(dv); d.model.push(m); d.nll.push(l); d.kind.push(kind); d.who.push(who[c]); }
    }
    s
}

/// Calibration-shift grid (m/z) and normalised log prior weights.
pub fn shift_grid(sig: f64) -> Vec<(f64, f64)> {
    let ds: Vec<f64> = (-6..=6).map(|i| 0.5 * sig * i as f64).collect();
    let lw: Vec<f64> = ds.iter().map(|d| -0.5 * (d / sig).powi(2)).collect();
    let z = logsumexp(&lw);
    ds.iter().zip(&lw).map(|(&d, &w)| (d, w - z)).collect()
}

#[derive(Clone, Debug)]
pub struct PanelResult { pub nll_marg: f64, pub nll_best: f64, pub delta_best: f64, pub split: Split, pub detail: ColDetail }

pub fn pen_of<'a>(p: &Panel, fig5: &'a Pen, s14: &'a Pen) -> &'a Pen { if p.pen == PenKind::Fig5 { fig5 } else { s14 } }

/// Likelihood of one panel for given species amplitudes (marginalised over the calibration shift).
pub fn eval_panel(p: &Panel, plac: &[Placed], pen: &Pen, amp: &[f64], rel: &[f64], normalise: bool, nz: &Noise) -> PanelResult {
    let grid = shift_grid(p.sig_cal);
    let mut terms = vec![];
    let (mut best, mut db, mut bs) = (f64::INFINITY, 0.0, Split::default());
    for &(d, lw) in &grid {
        let resp: Vec<Response> = plac.iter().map(|s| response(s, &p.grid, pen, d)).collect();
        let s = panel_nll(p, &resp, amp, rel, normalise, nz, None);
        terms.push(lw - s.total);
        if s.total < best { best = s.total; db = d; bs = s; }
    }
    let resp: Vec<Response> = plac.iter().map(|s| response(s, &p.grid, pen, db)).collect();
    let mut det = ColDetail::default();
    panel_nll(p, &resp, amp, rel, normalise, nz, Some(&mut det));
    PanelResult { nll_marg: -logsumexp(&terms), nll_best: best, delta_best: db, split: bs, detail: det }
}

// ---------------------------------------------------------------- pen kernel from the data

/// Measure the pen kernel from every isolated bright peak (>= 12 % of base) of the given panels: box width and tail
/// fractions (median over peaks, only where the tail is at least 2.5x the display floor).
pub fn fit_pen(panels: &[&Panel]) -> (Pen, usize) {
    let mut runs: Vec<(usize, usize, f64, f64, usize)> = vec![]; // (panel idx, p0 col index, plateau mean, floor, nbox)
    for (pi, p) in panels.iter().enumerate() {
        let n = p.grid.n;
        let mut c = 0;
        while c < n {
            let lo = c.saturating_sub(4); let hi = (c + 4).min(n - 1);
            let lm = p.h[lo..=hi].iter().cloned().fold(0.0, f64::max);
            if p.h[c] >= 12.0 && p.h[c] >= 0.93 * lm {
                let mut e = c; while e + 1 < n && p.h[e + 1] >= 0.93 * lm && (p.h[e + 1] - p.h[c]).abs() < 0.07 * p.h[c] { e += 1; }
                let len = e - c + 1;
                if (3..=5).contains(&len) { let mean = p.h[c..=e].iter().sum::<f64>() / len as f64; runs.push((pi, c, mean, p.floor, len)); }
                c = e + 1;
            } else { c += 1; }
        }
    }
    let mut lens: Vec<usize> = runs.iter().map(|r| r.4).collect(); lens.sort();
    let nbox = if lens.is_empty() { 4 } else { lens[lens.len() / 2] };
    let med = |mut v: Vec<f64>| -> Option<f64> { if v.len() < 5 { return None; } v.sort_by(|a, b| a.partial_cmp(b).unwrap()); Some(v[v.len() / 2]) };
    let side = |left: bool| -> Vec<f64> {
        let mut out = vec![];
        for d in 1..=12 {
            let mut v = vec![];
            for &(pi, p0, hm, fl, nb) in runs.iter().filter(|r| r.4 == nbox) {
                let _ = nb;
                let p = panels[pi];
                let idx = if left { p0 as i64 - d } else { (p0 + nbox - 1) as i64 + d };
                if idx < 0 || idx as usize >= p.grid.n { continue; }
                let val = p.h[idx as usize];
                if val >= 2.5 * fl && val < hm { v.push(val / hm); }
            }
            match med(v) { Some(x) => out.push(x), None => break }
        }
        out
    };
    let pen = Pen { nbox, left: side(true), right: side(false) };
    (pen, runs.iter().filter(|r| r.4 == nbox).count())
}

// ---------------------------------------------------------------- source (abundance) fit

#[derive(Clone, Debug)]
pub struct SourceFit { pub h: Vec<f64>, pub delta: f64, pub nll_best: f64, pub nll_marg: f64, pub split: Split }

fn golden(f: &mut dyn FnMut(f64) -> f64, mut a: f64, mut b: f64, it: usize) -> (f64, f64) {
    let g = 0.5 * (5f64.sqrt() - 1.0);
    let (mut x1, mut x2) = (b - g * (b - a), a + g * (b - a));
    let (mut f1, mut f2) = (f(x1), f(x2));
    for _ in 0..it {
        if f1 < f2 { b = x2; x2 = x1; f2 = f1; x1 = b - g * (b - a); f1 = f(x1); }
        else { a = x1; x1 = x2; f1 = f2; x2 = a + g * (b - a); f2 = f(x2); }
    }
    if f1 < f2 { (x1, f1) } else { (x2, f2) }
}

/// Fit plateau heights H_j (% units of the panel) of every species to one panel (no renormalisation).
/// `fixed[j] = Some(h)` pins a species. Coordinate descent, 1-D bracketing + golden section per species.
pub fn fit_source(p: &Panel, plac: &[Placed], pen: &Pen, nz: &Noise, init: &[f64], fixed: &[Option<f64>], sweeps: usize) -> SourceFit {
    let nsp = plac.len();
    let zeros = vec![0.0; nsp];
    let grid = shift_grid(p.sig_cal);
    let mut terms = vec![];
    let mut best: Option<SourceFit> = None;
    let mut h0 = init.to_vec();
    for &(d, lw) in &grid {
        let resp: Vec<Response> = plac.iter().map(|s| response(s, &p.grid, pen, d)).collect();
        let mut h = h0.clone();
        for j in 0..nsp { if let Some(v) = fixed[j] { h[j] = v; } }
        // neighbourhood maximum of the data around each species (upper bound of the search)
        let hmax: Vec<f64> = resp.iter().map(|r| {
            if r.v.is_empty() { return p.floor * 4.0; }
            let m = p.h[r.first..r.first + r.v.len()].iter().cloned().fold(0.0, f64::max);
            (m * 1.5).max(p.floor * 4.0)
        }).collect();
        for _ in 0..sweeps {
            for j in 0..nsp {
                if fixed[j].is_some() || resp[j].v.is_empty() { continue; }
                let mut cands = vec![0.0];
                let mut v = p.floor * 0.25; while v < hmax[j] { cands.push(v); v *= 1.8; }
                cands.push(hmax[j]);
                let mut hh = h.clone();
                let mut ev = |x: f64| -> f64 { hh[j] = x; panel_nll(p, &resp, &hh, &zeros, false, nz, None).total };
                let vals: Vec<f64> = cands.iter().map(|&x| ev(x)).collect();
                let ib = vals.iter().enumerate().fold(0, |a, (i, &x)| if x < vals[a] { i } else { a });
                let (lo, hi) = (cands[ib.saturating_sub(1)], cands[(ib + 1).min(cands.len() - 1)]);
                let (xb, fb) = golden(&mut ev, lo, hi, 22);
                h[j] = if fb < vals[ib] { xb } else { cands[ib] };
            }
        }
        let s = panel_nll(p, &resp, &h, &zeros, false, nz, None);
        terms.push(lw - s.total);
        if best.as_ref().map_or(true, |b| s.total < b.nll_best) {
            best = Some(SourceFit { h: h.clone(), delta: d, nll_best: s.total, nll_marg: 0.0, split: s });
        }
        h0 = h;
    }
    let mut b = best.unwrap();
    b.nll_marg = -logsumexp(&terms);
    b
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn normal_cdf_reference_values() {
        for (x, v) in [(-1.0, 0.15865525393145707), (-2.0, 0.022750131948179195), (-3.0, 0.0013498980316300946), (-5.0, 2.866515718791939e-7), (0.0, 0.5), (1.5, 0.9331927987311419), (-10.0, 7.619853024160527e-24)] {
            let r = ncdf(x);
            assert!(((r - v) / v).abs() < 1e-6, "Phi({}) = {} vs {}", x, r, v);
        }
        assert!((log_ncdf(-40.0) - (-804.608442013754)).abs() < 0.05, "{}", log_ncdf(-40.0));
    }
    #[test]
    fn shift_prior_normalised() {
        let g = shift_grid(1.0);
        assert!((logsumexp(&g.iter().map(|x| x.1).collect::<Vec<_>>())).abs() < 1e-12);
    }
}
