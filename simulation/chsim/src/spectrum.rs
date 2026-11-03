//! Spectrum forward model: `chsim spectrum <config>` (see chsim/SPECTRUM.md).
//!
//! source spectrum (Fig. 5C / S14 Part I, deconvolved to species)  x  channel transmission T_k(species) (cached simulations)
//! x downstream factor D(m/z)  ->  detector view (isotopes, Orbitrap line, pen kernel)  ->  render on the panel's pixel grid
//! (base-peak normalised, display floor censoring)  ->  Gaussian/Tobit likelihood against the published pixels.
//!
//! Subcommands: spectrum <cfg> | spectrum-selftest [cfg] | spectrum-grid <cfg> | spectrum-run <jobs.tsv>   (the last one
//! runs simulations on the GPU; everything else is CPU-only and never writes the result cache).

use crate::spectrum_detector::*;
use crate::spectrum_like::*;
use crate::spectrum_trans::*;
use std::collections::HashMap;
use std::fmt::Write as _;

// ---------------------------------------------------------------- config

pub struct Cfg(pub HashMap<String, String>);
impl Cfg {
    pub fn load(path: &str) -> Result<Cfg, String> {
        let text = std::fs::read_to_string(path).map_err(|e| format!("config {}: {}", path, e))?;
        let mut m = HashMap::new();
        for l in text.lines() {
            let l = l.split('#').next().unwrap().trim();
            if l.is_empty() { continue; }
            let (k, v) = l.split_once('=').ok_or_else(|| format!("bad config line: {}", l))?;
            m.insert(k.trim().to_string(), v.trim().to_string());
        }
        Ok(Cfg(m))
    }
    pub fn empty() -> Cfg { Cfg(HashMap::new()) }
    pub fn s(&self, k: &str, d: &str) -> String { self.0.get(k).cloned().unwrap_or(d.to_string()) }
    pub fn f(&self, k: &str, d: f64) -> f64 { self.0.get(k).and_then(|v| v.parse().ok()).unwrap_or(d) }
    pub fn opt(&self, k: &str) -> Option<f64> { self.0.get(k).and_then(|v| v.parse().ok()) }
}

/// The geometry and physics the whole R28 baseline cache was produced with (code geometry, gate 7 mm, 0.1 Pa, t_max 3 ms).
pub fn cache_keys() -> (String, String) {
    use crate::channel::{Geom, Physics};
    // same geometry as the R28 cache (reality/predict): s14_geom(code, gate 7 mm)
    let gargs: Vec<String> = vec!["".into(), "".into(), "".into(), "".into(), "".into(), "code".into(), "7.0".into()];
    let geom = crate::channel::s14_geom(&gargs); let _ = Geom::default;
    let mut ph = Physics::correct(0.1); ph.t_max = 3e-3;
    (format!("{:?}", geom), format!("{:?}", ph))
}

// ---------------------------------------------------------------- species definitions

pub struct SpDef { pub sim: SimSp, pub parent: Option<usize>, pub plac: Placed, pub gain: f64 }

pub fn mk(det: &Orbitrap, ubq: &Envelope, sim: SimSp, parent: Option<usize>) -> SpDef {
    let env = match sim.kind { Kind::Ubq(_) => ubq.clone(), _ => envelope_averagine((sim.mass - sim.z * PROTON).max(30.0)) };
    let l = detector_lines(&env, sim.z, det);
    SpDef { plac: Placed { apex_mz: sim.mz, lines: l.lines.clone() }, gain: l.apex_gain, sim, parent }
}

/// Fig. 5 species: ubiquitin 5+..12+, the non-ubiquitin peaks seen in 5C or in any 5A panel, and the same-charge adduct series.
pub fn fig5_species(det: &Orbitrap) -> Vec<SpDef> {
    let ubq = envelope_ubiquitin();
    let mut v: Vec<SpDef> = vec![];
    for z in 5..=12u32 { v.push(mk(det, &ubq, SimSp::new(&format!("ubq{}+", z), Kind::Ubq(z), 0.0), None)); }
    let ui = |z: u32| (z - 5) as usize;
    // chemical background (assumed singly charged); m/z = 5C/best-estimate values of SPECTRA_EXTRACTION.md section 6
    for mz in [170.1, 328.8, 402.6, 500.8, 574.6] { v.push(mk(det, &ubq, SimSp::new(&format!("m/z {:.0}", mz), Kind::Chem(1), mz), None)); }
    // same-charge adducts / satellites of ubiquitin (5C positions); parent = the ubiquitin charge state
    for (mz, z) in [(981.4, 9u32), (1003.5, 9), (1103.9, 8), (1125.25, 8), (1160.96, 8), (1261.85, 7), (1270.46, 7), (1276.08, 7), (1286.65, 7), (1416.6, 6)] {
        v.push(mk(det, &ubq, SimSp::new(&format!("ubq{}+ adduct {:.1}", z, mz), Kind::Adduct(z), mz), Some(ui(z))));
    }
    v
}

/// S14 species (peptide mixture). Charges from the printed labels; unlabelled peaks: 2+ if > 5 %, else 1+. Bak BH3 2+ and its
/// Na adduct use the exact simulated species of the scoreboard (mass 2 m/z, CCS 425 / 430).
pub fn s14_species(det: &Orbitrap, dir: &str) -> Result<Vec<SpDef>, String> {
    let ubq = envelope_ubiquitin();
    let mut v: Vec<SpDef> = vec![];
    let bak = SimSp { name: "Bak BH3 2+ 863.47".into(), kind: Kind::Chem(2), mz: 863.473, z: 2.0, mass: 2.0 * 863.473, ccs: 425.0 };
    let na = SimSp { name: "[M+H+Na]2+ 874.47".into(), kind: Kind::Adduct(2), mz: 874.466, z: 2.0, mass: 2.0 * 874.466, ccs: 430.0 };
    let kk = SimSp { name: "[M+H+K]2+ 882.42".into(), kind: Kind::Adduct(2), mz: 882.42, z: 2.0, mass: 2.0 * 882.42, ccs: 435.0 };
    v.push(mk(det, &ubq, bak, None));
    v.push(mk(det, &ubq, na, Some(0)));
    v.push(mk(det, &ubq, kk, Some(0)));
    let text = std::fs::read_to_string(format!("{}/figS14_peak_catalog.csv", dir)).map_err(|e| format!("S14 catalog: {}", e))?;
    for l in text.lines().filter(|l| !l.starts_with('#') && !l.starts_with("part")) {
        let p: Vec<&str> = l.split(',').collect();
        if p.len() < 6 || p[0] != "I" { continue; }
        let (mz_cat, h): (f64, f64) = (p[1].parse().unwrap_or(0.0), p[2].parse().unwrap_or(0.0));
        if (mz_cat - 863.14).abs() < 1.0 || (mz_cat - 874.31).abs() < 1.0 { continue; } // = the Bak species above
        let mz = p[4].parse::<f64>().unwrap_or(mz_cat);
        let z = p[5].parse::<u32>().ok().unwrap_or(if h > 5.0 { 2 } else { 1 });
        let kind = Kind::Chem(z);
        v.push(mk(det, &ubq, SimSp::new(&format!("m/z {:.1} z={}", mz, z), kind, mz), None));
    }
    // species that appear only in Part II (AC on) at < 2.5 %: 854.74 (2+), 530.65 / 728.46 are the catalogued 531 / 728 peaks
    v.push(mk(det, &ubq, SimSp::new("m/z 854.7 z=2", Kind::Chem(2), 854.74), None));
    Ok(v)
}

// ---------------------------------------------------------------- context

pub struct Ctx {
    pub cfg: Cfg,
    pub det: Orbitrap,
    pub nz: Noise,
    pub fig5: Vec<Panel>,   // A1 A2 A3 A4 B C
    pub s14: Vec<Panel>,    // S14I S14II
    pub defs5: Vec<SpDef>,
    pub defs14: Vec<SpDef>,
    pub pen5: Pen,
    pub pen14: Pen,
    pub pen_anchors: usize,
    pub src5: SourceFit,
    pub src14: SourceFit,
    pub tm: TModel,
    pub out: String,
    pub mask_chem: bool, // treat every chemical-background species as MISSING (identical at all grid points)
}

impl Ctx {
    pub fn panel(&self, name: &str) -> &Panel { self.fig5.iter().chain(self.s14.iter()).find(|p| p.name == name).unwrap() }
    pub fn plac5(&self) -> Vec<Placed> { self.defs5.iter().map(|d| d.plac.collapsed()).collect() }
    pub fn plac14(&self) -> Vec<Placed> { self.defs14.iter().map(|d| d.plac.clone()).collect() }
}

pub fn init_heights(p: &Panel, defs: &[SpDef]) -> Vec<f64> {
    defs.iter().map(|d| {
        let x = p.grid.x_of(d.sim.mz).round() as i64 - p.grid.col0;
        let (a, b) = ((x - 2).max(0) as usize, ((x + 2).max(0) as usize).min(p.grid.n - 1));
        p.h[a..=b].iter().cloned().fold(0.0, f64::max).max(0.3 * p.floor)
    }).collect()
}

/// 5C source fit; peaks of 5C that the species list does not explain (data above 1.5x the display floor, model zero) are
/// added as extra species (adduct of the nearest ubiquitin charge state when within 8 % in m/z, else singly charged chemical) and refitted.
pub fn fit_5c(defs5: &mut Vec<SpDef>, pc: &Panel, pen5: &Pen, nz: &Noise, det: &Orbitrap, frac: f64, rounds: usize) -> SourceFit {
    let ubq = envelope_ubiquitin();
    let pinned = |d: &SpDef| matches!(d.sim.kind, Kind::Chem(_)) && [328.8, 500.8, 574.6].iter().any(|m| (d.sim.mz - m).abs() < 0.01);
    let mut fit;
    let mut round = 0;
    loop {
        let plac: Vec<Placed> = defs5.iter().map(|d| d.plac.collapsed()).collect();
        let fixed: Vec<Option<f64>> = defs5.iter().map(|d| if pinned(d) { Some(frac * pc.floor) } else { None }).collect();
        fit = fit_source(pc, &plac, pen5, nz, &init_heights(pc, defs5), &fixed, 3);
        if round >= rounds { break; }
        round += 1;
        let resp: Vec<Response> = plac.iter().map(|s| response(s, &pc.grid, pen5, fit.delta)).collect();
        let zeros = vec![0.0; defs5.len()];
        let mut det_c = ColDetail::default();
        panel_nll(pc, &resp, &fit.h, &zeros, false, nz, Some(&mut det_c));
        // unexplained columns -> peaks: data well above the model (residual >= 1.5x the display floor and >= 40 % of the data), grouped within 3 px
        let resid: Vec<f64> = (0..det_c.mz.len()).map(|i| if det_c.data[i] >= pc.floor { det_c.data[i] - det_c.model[i] } else { 0.0 }).collect();
        let cols: Vec<usize> = (0..det_c.mz.len()).filter(|&i| resid[i] >= 1.5 * pc.floor && resid[i] >= 0.4 * det_c.data[i]).collect();
        let mut groups: Vec<Vec<usize>> = vec![];
        for &c in &cols { match groups.last_mut() { Some(g) if c - g[g.len() - 1] <= 3 => g.push(c), _ => groups.push(vec![c]) } }
        let mut added = 0;
        for g in groups {
            if defs5.len() >= 60 { break; }
            let w: f64 = g.iter().map(|&i| resid[i]).sum();
            let mz = g.iter().map(|&i| det_c.mz[i] * resid[i]).sum::<f64>() / w - fit.delta;
            if defs5.iter().any(|d| (d.sim.mz - mz).abs() < 1.5 * pc.grid.dmz) { continue; }
            let zn = (5..=12u32).min_by(|&a, &b| ((mz / (M_UBQ_AVG / a as f64 + PROTON) - 1.0).abs()).partial_cmp(&((mz / (M_UBQ_AVG / b as f64 + PROTON) - 1.0).abs())).unwrap()).unwrap();
            let mzz = (M_UBQ_AVG + zn as f64 * 1.00728) / zn as f64;
            let near = (mz / mzz - 1.0).abs() < 0.08;
            let sim = if near { SimSp::new(&format!("5C peak {:.1} (ubq{}+ adduct)", mz, zn), Kind::Adduct(zn), mz) } else { SimSp::new(&format!("5C peak {:.1}", mz), Kind::Chem(1), mz) };
            let parent = if near { Some((zn - 5) as usize) } else { None };
            defs5.push(mk(det, &ubq, sim, parent));
            added += 1;
        }
        if added == 0 { break; }
    }
    fit
}

pub fn build_ctx(cfg: Cfg) -> Result<Ctx, String> {
    let mut data = cfg.s("data", "data/their_side/figure_data");
    if !std::path::Path::new(&data).exists() { for alt in ["../data/their_side/figure_data", "../../data/their_side/figure_data"] { if std::path::Path::new(alt).exists() { data = alt.to_string(); break; } } }
    let cache = std::env::var("CHSIM_CACHE").unwrap_or(cfg.s("cache", "cache/results.tsv"));
    let det = Orbitrap { r200: cfg.f("r200", 70000.0) };
    let nz = Noise { sigq_px: cfg.f("sigq", 0.5), rho_p: cfg.f("rho_p", 0.015), rho_t: cfg.f("rho_t", 0.25) };
    let fig5 = load_fig5(&data)?;
    let s14 = load_s14(&data)?;
    let (gk, pk) = cache_keys();
    let tcache = TCache::load(&cache, &gk, &pk)?;
    let mode = TMode::parse(&cfg.s("tmode", "on"))?;
    let down = Downstream::parse(&cfg.s("downstream", "unity"))?;
    let tm = TModel { cache: tcache, mode, down, n_new: cfg.f("n_new", 3000.0) as usize, interp_max_dmz: cfg.f("interp_dmz", 15.0) };
    let mut defs5 = fig5_species(&det);
    let defs14 = s14_species(&det, &data)?;

    let refs: Vec<&Panel> = fig5.iter().collect();
    let (pen5, pen_anchors) = fit_pen(&refs);
    let pen14 = Pen::bare(cfg.f("s14_nbox", 1.0) as usize);

    // ---- source fits: 5C (Fig. 5 control) and S14 Part I (gate open)
    let pc = fig5.iter().find(|p| p.name == "C").unwrap();
    let frac = cfg.f("src_floor_frac", 0.5);
    let src5 = fit_5c(&mut defs5, pc, &pen5, &nz, &det, frac, cfg.f("augment_rounds", 4.0) as usize);

    // S14 Part I: only Bak 2+ and its adducts are fitted (they decide Part II); the other peaks of the peptide mixture are pinned at
    // the apex column of the data (the S14 pen is not shift-invariant: neighbouring columns carry 30-50 % of a peak, see SPECTRUM.md)
    let p1 = s14.iter().find(|p| p.name == "S14I").unwrap();
    let plac14: Vec<Placed> = defs14.iter().map(|d| d.plac.clone()).collect();
    let init14 = init_heights(p1, &defs14);
    let fixed14: Vec<Option<f64>> = defs14.iter().enumerate().map(|(j, d)| {
        if (d.sim.mz - 882.42).abs() < 0.01 || (d.sim.mz - 854.74).abs() < 0.01 { Some(frac * p1.floor) } else if j <= 1 { None } else { Some(init14[j]) }
    }).collect();
    let src14 = fit_source(p1, &plac14, &pen14, &nz, &init14, &fixed14, 3);

    let out = cfg.s("out", "spectrum_out");
    let mask_chem = cfg.f("mask_chem", 0.0) > 0.5;
    Ok(Ctx { cfg, det, nz, fig5, s14, defs5, defs14, pen5, pen14, pen_anchors, src5, src14, tm, out, mask_chem })
}

// ---------------------------------------------------------------- predictions

pub struct Pred { pub amp: Vec<f64>, pub rel: Vec<f64>, pub tv: Vec<Vec<TVal>> /* [channel][species] */ }

fn pen_for<'a>(ctx: &'a Ctx, p: &Panel) -> &'a Pen { pen_of(p, &ctx.pen5, &ctx.pen14) }

/// amp_j = H_j D(mz_j) sum_k T_k(species j); rel = MC relative error of that sum. `freqs` = channel AC frequencies [kHz].
pub fn predict(ctx: &Ctx, defs: &[SpDef], h: &[f64], pt: &Point, freqs: &[f64], mask_chem: bool) -> Pred {
    let mut amp = vec![0.0; defs.len()];
    let mut var = vec![0.0; defs.len()];
    let mut tv: Vec<Vec<TVal>> = vec![];
    for &f in freqs {
        let mut row = vec![];
        for (j, d) in defs.iter().enumerate() {
            let parent = d.parent.map(|i| &defs[i].sim);
            let t = if mask_chem && matches!(d.sim.kind, Kind::Chem(_)) { TVal::missing() } else { ctx.tm.t(&d.sim, pt, f, parent) };
            let a = h[j] * ctx.tm.down.d(d.sim.mz) * t.t;
            amp[j] += a; var[j] += (a * t.rel).powi(2);
            row.push(t);
        }
        tv.push(row);
    }
    let rel: Vec<f64> = amp.iter().zip(&var).map(|(&a, &v)| if a > 0.0 { v.sqrt() / a } else { 0.0 }).collect();
    Pred { amp, rel, tv }
}

/// Channel frequencies: config `f9 f8 f7 f6` (kHz) if given, else the argmax of the cached resonance scan of the channel's charge state.
pub fn channel_freqs(ctx: &Ctx, pt: &Point) -> Vec<f64> {
    let mut fs = vec![];
    for z in [9u32, 8, 7, 6] {
        if let Some(f) = ctx.cfg.opt(&format!("f{}", z)) { fs.push(f); continue; }
        let sp = SimSp::new("", Kind::Ubq(z), 0.0);
        let mut scan = ctx.tm.cache.adaptive_scan_pts(sp.mass, sp.z, sp.ccs, pt.vrf, pt.vac, &pt.dc, ctx.cfg.f("scan_seed_mult", 2000.0) as u64 * z as u64);
        if scan.is_empty() { scan = ctx.tm.cache.freq_scan_best_n(sp.mass, sp.z, sp.ccs, pt.vrf, pt.vac, &pt.dc); }
        let best = scan.iter().fold((f64::NAN, -1.0), |a, (f, c)| { let t = c.transmitted as f64 / c.n as f64; if t > a.1 { (*f, t) } else { a } });
        fs.push(best.0);
    }
    fs
}

pub fn fig5_point(ctx: &Ctx) -> Point {
    Point { vrf: ctx.cfg.f("vrf5", 56.0), frf: 500e3, vac: ctx.cfg.f("vac5", 0.2), dc: [18.0, ctx.cfg.f("gate5", 24.0), 15.0, 20.0] }
}
pub fn s14_point(ctx: &Ctx) -> Point {
    Point { vrf: ctx.cfg.f("vrf14", 70.0), frf: 500e3, vac: ctx.cfg.f("vac14", 0.35), dc: [18.0, ctx.cfg.f("gate14", 24.0), 15.0, 20.0] }
}
pub fn s14_freq(ctx: &Ctx, pt: &Point) -> f64 {
    if let Some(f) = ctx.cfg.opt("f14") { return f; }
    let b = &ctx.defs14[0].sim;
    let scan = ctx.tm.cache.freq_scan_best_n(b.mass, b.z, b.ccs, pt.vrf, pt.vac, &pt.dc);
    scan.iter().fold((83.0, -1.0), |a, (f, c)| { let t = c.transmitted as f64 / c.n as f64; if t > a.1 { (*f, t) } else { a } }).0
}

pub struct Eval { pub name: String, pub res: PanelResult, pub pred: Option<Pred> }

/// A1..A4 (one channel active), B (all four), against the 5C-derived source and the cached transmissions.
pub fn eval_fig5(ctx: &Ctx, pt: &Point, freqs: &[f64]) -> Vec<Eval> {
    let plac = ctx.plac5();
    let mut out = vec![];
    for (i, nm) in ["A1", "A2", "A3", "A4"].iter().enumerate() {
        let p = ctx.panel(nm);
        let pr = predict(ctx, &ctx.defs5, &ctx.src5.h, pt, &freqs[i..=i], ctx.mask_chem);
        let res = eval_panel(p, &plac, pen_for(ctx, p), &pr.amp, &pr.rel, true, &ctx.nz);
        out.push(Eval { name: nm.to_string(), res, pred: Some(pr) });
    }
    let p = ctx.panel("B");
    let pr = predict(ctx, &ctx.defs5, &ctx.src5.h, pt, freqs, ctx.mask_chem);
    let res = eval_panel(p, &plac, pen_for(ctx, p), &pr.amp, &pr.rel, true, &ctx.nz);
    out.push(Eval { name: "B".into(), res, pred: Some(pr) });
    out
}

pub fn eval_s14(ctx: &Ctx, pt: &Point, f: f64) -> Eval {
    let p = ctx.panel("S14II");
    let pr = predict(ctx, &ctx.defs14, &ctx.src14.h, pt, &[f], false);
    let res = eval_panel(p, &ctx.plac14(), pen_for(ctx, p), &pr.amp, &pr.rel, true, &ctx.nz);
    Eval { name: "S14II".into(), res, pred: Some(pr) }
}

// ---------------------------------------------------------------- residual features

pub struct Feature { pub mz_lo: f64, pub mz_hi: f64, pub nll: f64, pub data_max: f64, pub model_max: f64, pub label: String }

pub fn features(det: &ColDetail, defs: &[SpDef], thr: f64) -> Vec<Feature> {
    let n = det.mz.len();
    let mut fs: Vec<Feature> = vec![];
    let mut c = 0;
    while c < n {
        if det.nll[c] > thr {
            let mut e = c; let mut last = c;
            while e + 1 < n && e + 1 <= last + 2 { e += 1; if det.nll[e] > thr { last = e; } }
            let e = last;
            let nll: f64 = det.nll[c..=e].iter().sum();
            let (dmax, mmax) = (det.data[c..=e].iter().cloned().fold(0.0, f64::max), det.model[c..=e].iter().cloned().fold(0.0, f64::max));
            let ic = (c..=e).fold(c, |a, i| if det.data[i].max(det.model[i]) > det.data[a].max(det.model[a]) { i } else { a });
            let mzc = det.mz[ic];
            let near = defs.iter().filter(|d| (d.sim.mz - mzc).abs() < 4.0 * (det.mz.get(1).unwrap_or(&0.0) - det.mz[0]).abs().max(1.0)).map(|d| d.sim.name.clone()).collect::<Vec<_>>();
            let unmod = (c..=e).any(|i| det.kind[i] == 2);
            let what = if mmax <= 0.0 { "UNMODELLED in model" } else if dmax < mmax { "model OVER-predicts" } else { "model UNDER-predicts" };
            fs.push(Feature { mz_lo: det.mz[c], mz_hi: det.mz[e], nll, data_max: dmax, model_max: mmax, label: format!("{}{}{}", what, if near.is_empty() { String::new() } else { format!(" [{}]", near.join(", ")) }, if unmod { "" } else { "" }) });
            c = e + 1;
        } else { c += 1; }
    }
    fs.sort_by(|a, b| b.nll.partial_cmp(&a.nll).unwrap());
    fs
}

// ---------------------------------------------------------------- output helpers

fn write_detail(path: &str, name: &str, det: &ColDetail, defs: &[SpDef]) -> std::io::Result<()> {
    let mut s = format!("# {}: model vs data per pixel column. mz, data (% of base), model (same normalisation), resid = data - model, nll = -log density of the column, kind: censored (data below display floor) | gauss | UNMODELLED (data above floor, model zero); species = dominant model term\nmz,data,model,resid,nll,kind,species\n", name);
    for i in 0..det.mz.len() {
        let k = ["censored", "gauss", "UNMODELLED"][det.kind[i] as usize];
        let sp = if det.who[i] == usize::MAX { String::new() } else { defs[det.who[i]].sim.name.clone() };
        let _ = writeln!(s, "{:.3},{:.4},{:.4},{:.4},{:.4},{},{}", det.mz[i], det.data[i], det.model[i], det.data[i] - det.model[i], det.nll[i], k, sp);
    }
    std::fs::write(path, s)
}

/// Cheap SVG: top = full scale, bottom = zoom to 0-30 % with the display floor drawn; data black, model red.
fn write_svg(path: &str, title: &str, det: &ColDetail, floor: f64) -> std::io::Result<()> {
    let (w, h1, h2, ml) = (1100.0, 220.0, 220.0, 50.0);
    let (x0, x1) = (det.mz[0], *det.mz.last().unwrap());
    let sx = |m: f64| ml + (m - x0) / (x1 - x0) * (w - ml - 10.0);
    let mut s = format!("<svg xmlns='http://www.w3.org/2000/svg' width='{}' height='{}' font-family='monospace' font-size='11'>\n<rect width='100%' height='100%' fill='white'/>\n<text x='{}' y='14'>{} - data (black) vs forward model (red); lower plot zoomed to 0-30 %, dashed = display floor</text>\n", w, h1 + h2 + 70.0, ml, title);
    for (k, (top, hh, ymax)) in [(25.0, h1, 105.0), (25.0 + h1 + 30.0, h2, 30.0)].iter().enumerate() {
        let sy = |v: f64| top + hh - (v.min(*ymax) / ymax) * hh;
        let _ = writeln!(s, "<rect x='{}' y='{}' width='{}' height='{}' fill='none' stroke='#888'/>", ml, top, w - ml - 10.0, hh);
        for t in (0..=((*ymax / 10.0) as i32)).map(|i| i as f64 * if k == 0 { 20.0 } else { 5.0 }).filter(|t| *t <= *ymax) { let _ = writeln!(s, "<text x='4' y='{:.1}'>{}</text><line x1='{}' x2='{}' y1='{:.1}' y2='{:.1}' stroke='#ddd'/>", sy(t) + 4.0, t, ml, w - 10.0, sy(t), sy(t)); }
        for m in (x0 / 200.0).ceil() as i64..=(x1 / 200.0) as i64 { let x = sx(m as f64 * 200.0); let _ = writeln!(s, "<line x1='{:.1}' x2='{:.1}' y1='{}' y2='{}' stroke='#eee'/><text x='{:.1}' y='{}'>{}</text>", x, x, top, top + hh, x - 12.0, top + hh + 12.0, m * 200); }
        if k == 1 { let _ = writeln!(s, "<line x1='{}' x2='{}' y1='{:.1}' y2='{:.1}' stroke='#c60' stroke-dasharray='4 3'/>", ml, w - 10.0, sy(floor), sy(floor)); }
        for (vals, col) in [(&det.data, "#000"), (&det.model, "#d22")] {
            let pts: Vec<String> = (0..det.mz.len()).map(|i| format!("{:.1},{:.1}", sx(det.mz[i]), sy(vals[i]))).collect();
            let _ = writeln!(s, "<polyline fill='none' stroke='{}' stroke-width='{}' points='{}'/>", col, if col == "#000" { 1.2 } else { 1.0 }, pts.join(" "));
        }
    }
    s.push_str("</svg>\n");
    std::fs::write(path, s)
}

pub fn dedupe_needs(v: Vec<Need>) -> Vec<Need> {
    let mut seen: HashMap<String, usize> = HashMap::new();
    let mut out: Vec<Need> = vec![];
    for n in v {
        let k = format!("{:.6}|{}|{:.4}|{}|{}|{}|{:?}", n.mass, n.z, n.ccs, n.vrf, n.vac, n.fac, n.dc);
        if seen.contains_key(&k) { continue; }
        seen.insert(k, out.len()); out.push(n);
    }
    out
}

pub fn write_needs(path: &str, needs: &[Need]) -> std::io::Result<()> {
    let mut s = String::from("# exact Jobs (channel::Job) still needed so every species x channel is an Exact cache hit. Run on the GPU with:\n#   cd chsim && taskpolicy -b ./target/release/chsim spectrum-run <this file> [--dry parses only]      (writes cache/results.tsv through scan::cached_batch)\n# columns: label mass[Da incl. charge protons] z ccs_n2[A^2] v_rf f_rf v_ac f_ac[Hz] n seed dc0 dc1 dc2 dc3 ; geometry = code, gate hole 7 mm; Physics::correct(0.1 Pa), t_max 3 ms (same as the R28 cache)\n");
    for (i, n) in needs.iter().enumerate() {
        let _ = writeln!(s, "{}\t{:.6}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", n.label.replace('\t', " "), n.mass, n.z, n.ccs, n.vrf, n.frf, n.vac, n.fac, n.n, JOB_SEED0 + i as u64, n.dc[0], n.dc[1], n.dc[2], n.dc[3]);
    }
    std::fs::write(path, s)
}

// ---------------------------------------------------------------- the command

pub fn cmd(args: &[String]) {
    let r = match args.get(0).map(|s| s.as_str()) {
        Some("spectrum") => run(args),
        Some("spectrum-selftest") => crate::spectrum_val::selftest(args),
        Some("spectrum-grid") => crate::spectrum_val::grid(args),
        Some("spectrum-res") => crate::spectrum_val::res_scan(args),
        Some("spectrum-run") => run_jobs(args),
        _ => Err("usage: chsim spectrum <config> | spectrum-selftest [config] | spectrum-grid <config> | spectrum-run <jobs.tsv>".into()),
    };
    if let Err(e) = r { eprintln!("error: {}", e); std::process::exit(1); }
}

fn run(args: &[String]) -> Result<(), String> {
    let cfgp = args.get(1).ok_or("usage: chsim spectrum <config>")?;
    let cfg = Cfg::load(cfgp)?;
    let ctx = build_ctx(cfg)?;
    std::fs::create_dir_all(&ctx.out).map_err(|e| e.to_string())?;
    let (pt5, pt14) = (fig5_point(&ctx), s14_point(&ctx));
    let freqs = channel_freqs(&ctx, &pt5);
    let f14 = s14_freq(&ctx, &pt14);

    println!("SPECTRUM FORWARD MODEL  config {}", cfgp);
    println!("  cache: {} usable entries ({} lines in file; {} with a different geometry/physics ignored)", ctx.tm.cache.len(), ctx.tm.cache.n_lines, ctx.tm.cache.n_other_geom);
    println!("  Fig.5 point: RF {} V, AC +-{} V, DC {:?}, channels 9+/8+/7+/6+ at {:?} kHz   S14 point: RF {} V, AC +-{} V, {:.2} kHz", pt5.vrf, pt5.vac, pt5.dc, freqs, pt14.vrf, pt14.vac, f14);
    println!("  T mode {:?}, downstream {:?}, Orbitrap R(200) = {}, noise: sigma_q {} px, rho_plateau {}, rho_tail {}", ctx.tm.mode, ctx.tm.down, ctx.det.r200, ctx.nz.sigq_px, ctx.nz.rho_p, ctx.nz.rho_t);
    println!("\n[1] pen kernel (how the figure draws one line; measured on {} isolated bright peaks of 5A/5B/5C): plateau {} columns", ctx.pen_anchors, ctx.pen5.nbox);
    println!("    low-m/z side (cols outside the box, nearest first): {:?}", ctx.pen5.left.iter().map(|x| (x * 1000.0).round() / 1000.0).collect::<Vec<_>>());
    println!("    high-m/z side:                                      {:?}", ctx.pen5.right.iter().map(|x| (x * 1000.0).round() / 1000.0).collect::<Vec<_>>());
    println!("    S14 (Xcalibur screenshot): box {} column, no tails (isotope structure is modelled explicitly)", ctx.pen14.nbox);

    // ---- source tables
    println!("\n[2] source model: Fig. 5C deconvolved to species (plateau height H in % of the 5C base peak; abundance = H / apex gain, 8+ = 100)");
    println!("    5C fit: calibration shift {:+.2} m/z, -log L {:.1} (marginal {:.1}); unmodelled-data columns: {} (-log L {:.1})", ctx.src5.delta, ctx.src5.nll_best, ctx.src5.nll_marg, ctx.src5.split.n_unmod, ctx.src5.split.unmodelled);
    let ab8 = ctx.src5.h[3] / ctx.defs5[3].gain;
    println!("    {:<26} {:>9} {:>8} {:>10} {:>9}", "species", "m/z", "H [%]", "abund.", "FWHM[m/z]");
    let mut sp_csv = String::from("species,mz,z,H_pct,abundance_rel_8plus,apex_gain,fwhm_mz,kind\n");
    for (j, d) in ctx.defs5.iter().enumerate() {
        let a = 100.0 * (ctx.src5.h[j] / d.gain) / ab8;
        println!("    {:<26} {:>9.2} {:>8.3} {:>10.3} {:>9.4}{}", d.sim.name, d.sim.mz, ctx.src5.h[j], a, ctx.det.fwhm(d.sim.mz), if ctx.src5.h[j] < 1.5 * ctx.panel("C").floor { "   (near/below 5C floor)" } else { "" });
        let _ = writeln!(sp_csv, "{},{:.3},{},{:.4},{:.4},{:.5},{:.5},{:?}", d.sim.name, d.sim.mz, d.sim.z, ctx.src5.h[j], a, d.gain, ctx.det.fwhm(d.sim.mz), d.sim.kind);
    }
    std::fs::write(format!("{}/species_source_5C.csv", ctx.out), sp_csv).map_err(|e| e.to_string())?;

    // ---- predictions
    let plac5 = ctx.plac5();
    let ev5 = eval_fig5(&ctx, &pt5, &freqs);
    let ev14 = eval_s14(&ctx, &pt14, f14);
    let pc = ctx.panel("C");
    let mut c_det = ColDetail::default();
    {
        let resp: Vec<Response> = plac5.iter().map(|s| response(s, &pc.grid, &ctx.pen5, ctx.src5.delta)).collect();
        let zeros = vec![0.0; plac5.len()];
        panel_nll(pc, &resp, &ctx.src5.h, &zeros, false, &ctx.nz, Some(&mut c_det));
    }
    let p1 = ctx.panel("S14I");
    let mut s14i_det = ColDetail::default();
    {
        let pl = ctx.plac14();
        let resp: Vec<Response> = pl.iter().map(|s| response(s, &p1.grid, &ctx.pen14, ctx.src14.delta)).collect();
        let zeros = vec![0.0; pl.len()];
        panel_nll(p1, &resp, &ctx.src14.h, &zeros, false, &ctx.nz, Some(&mut s14i_det));
    }

    println!("\n[3] channel transmission T_k(species) used (per pass; provenance: exact cache pool | interp | proxy | MISSING -> contributes 0)");
    println!("    {:<26} {:>8} {}", "species", "H*D [%]", ["ch9+", "ch8+", "ch7+", "ch6+"].iter().map(|s| format!("{:>22}", s)).collect::<Vec<_>>().join(""));
    let mut tcsv = String::from("species,mz,channel_z,f_kHz,T,rel_err,provenance,k_transmitted,n\n");
    for (j, d) in ctx.defs5.iter().enumerate() {
        let cells: Vec<String> = (0..4).map(|k| { let t = &ev5[k].pred.as_ref().unwrap().tv[0][j]; format!("{:>10.5} {:<11}", t.t, t.prov.tag()) }).collect();
        println!("    {:<26} {:>8.2} {}", d.sim.name, ctx.src5.h[j], cells.join(" "));
        for k in 0..4 { let t = &ev5[k].pred.as_ref().unwrap().tv[0][j]; let _ = writeln!(tcsv, "{},{:.2},{},{:.2},{:.6},{:.4},{},{},{}", d.sim.name, d.sim.mz, [9, 8, 7, 6][k], freqs[k], t.t, t.rel, t.prov.tag(), t.k_on, t.n_on); }
    }
    std::fs::write(format!("{}/transmission_fig5.csv", ctx.out), tcsv).map_err(|e| e.to_string())?;

    println!("\n[4] likelihood (−log L; lower is better; calibration shift marginalised with the stated prior)");
    println!("    {:<8} {:>6} {:>12} {:>12} {:>8} {:>9} {:>9} {:>14}", "panel", "role", "-logL marg", "-logL best", "shift", "n_gauss", "n_unmod", "unmodelled part");
    let row = |nm: &str, role: &str, marg: f64, best: f64, d: f64, s: &Split| println!("    {:<8} {:>6} {:>12.1} {:>12.1} {:>+8.2} {:>9} {:>9} {:>14.1}", nm, role, marg, best, d, s.n_gauss, s.n_unmod, s.unmodelled);
    row("C", "source", ctx.src5.nll_marg, ctx.src5.nll_best, ctx.src5.delta, &ctx.src5.split);
    row("S14I", "source", ctx.src14.nll_marg, ctx.src14.nll_best, ctx.src14.delta, &ctx.src14.split);
    let mut tot = 0.0; let mut tot_um = 0.0;
    for e in ev5.iter().chain(std::iter::once(&ev14)) { row(&e.name, "output", e.res.nll_marg, e.res.nll_best, e.res.delta_best, &e.res.split); tot += e.res.nll_marg; tot_um += e.res.split.unmodelled; }
    println!("    (S14I: only Bak 2+ and its adducts are fitted, the other peptide peaks are pinned at their apex column, so its -logL is not a goodness of fit; the S14 pen is not shift-invariant)");
    println!("    TOTAL (A1-A4, B, S14II) -logL = {:.1}   of which unmodelled-data columns {:.1}   modelled/censored part {:.1}", tot, tot_um, tot - tot_um);

    // ---- residual features
    println!("\n[5] dominant residual features per panel (columns with -log density > 4, merged within 2 px)");
    let mut summary = String::new();
    for e in ev5.iter().chain(std::iter::once(&ev14)) {
        let defs: &[SpDef] = if e.name.starts_with("S14") { &ctx.defs14 } else { &ctx.defs5 };
        let fs = features(&e.res.detail, defs, 4.0);
        let tot_f: f64 = fs.iter().map(|f| f.nll).sum();
        let _ = writeln!(summary, "  {} (-logL {:.1}; features {:.1}):", e.name, e.res.nll_marg, tot_f);
        for f in fs.iter().take(6) { let _ = writeln!(summary, "      m/z {:7.1}-{:7.1}  -logL {:9.1}  data {:6.2} %  model {:6.2} %  {}", f.mz_lo, f.mz_hi, f.nll, f.data_max, f.model_max, f.label); }
    }
    for (nm, det, defs) in [("C (source fit)", &c_det, &ctx.defs5), ("S14I (source fit)", &s14i_det, &ctx.defs14)] {
        let fs = features(det, defs, 4.0);
        let _ = writeln!(summary, "  {} (features {:.1}):", nm, fs.iter().map(|f| f.nll).sum::<f64>());
        for f in fs.iter().take(6) { let _ = writeln!(summary, "      m/z {:7.1}-{:7.1}  -logL {:9.1}  data {:6.2} %  model {:6.2} %  {}", f.mz_lo, f.mz_hi, f.nll, f.data_max, f.model_max, f.label); }
    }
    print!("{}", summary);

    // ---- files
    for e in &ev5 { let p = ctx.panel(&e.name); write_detail(&format!("{}/{}_model_vs_data.csv", ctx.out, e.name), &e.name, &e.res.detail, &ctx.defs5).map_err(|x| x.to_string())?; write_svg(&format!("{}/{}.svg", ctx.out, e.name), &e.name, &e.res.detail, p.floor).map_err(|x| x.to_string())?; }
    write_detail(&format!("{}/S14II_model_vs_data.csv", ctx.out), "S14II", &ev14.res.detail, &ctx.defs14).map_err(|x| x.to_string())?;
    write_svg(&format!("{}/S14II.svg", ctx.out), "S14 Part II", &ev14.res.detail, ctx.panel("S14II").floor).map_err(|x| x.to_string())?;
    write_detail(&format!("{}/C_source_fit.csv", ctx.out), "C (source fit)", &c_det, &ctx.defs5).map_err(|x| x.to_string())?;
    write_svg(&format!("{}/C.svg", ctx.out), "5C source fit", &c_det, pc.floor).map_err(|x| x.to_string())?;
    write_detail(&format!("{}/S14I_source_fit.csv", ctx.out), "S14 Part I (source fit)", &s14i_det, &ctx.defs14).map_err(|x| x.to_string())?;
    std::fs::write(format!("{}/residual_features.txt", ctx.out), &summary).map_err(|e| e.to_string())?;

    // ---- missing jobs: (1) cells with no exact entry, (2) statistics upgrades for exact cells whose leak would be visible but whose MC error is large
    let n_leak = ctx.cfg.f("n_leak", 30000.0) as usize;
    let chan_mz = [952.65, 1071.61, 1224.55, 1428.47];
    // priority buckets: 0 ubiquitin charge states (the cross-talk matrix), 1 chemical background in Fig. 5, 2 S14 Bak 2+ and its adducts,
    // 3 ubiquitin adducts / 5C satellites (proxy = parent until simulated), 4 the other S14 peptides
    let mut bucket: Vec<Vec<Need>> = vec![vec![]; 5];
    for (j, d) in ctx.defs5.iter().enumerate() {
        let b = match d.sim.kind { Kind::Ubq(_) => 0, Kind::Chem(_) => 1, Kind::Adduct(_) => 3 };
        for (k, &f) in freqs.iter().enumerate() {
            let near = (d.sim.mz / chan_mz[k] - 1.0).abs() < 0.08;
            let nn = if near { ctx.tm.n_new } else { n_leak };
            bucket[b].extend(ctx.tm.needs(&d.sim, &pt5, f, &format!("P{} fig5 ch{}+", b, [9, 8, 7, 6][k]), nn));
            // statistics upgrade of an exact cell whose leak would be visible but whose MC error is large
            let pr = ev5[k].pred.as_ref().unwrap();
            let t = &pr.tv[0][j];
            let amax = (0..ctx.defs5.len()).map(|i| pr.amp[i]).fold(0.0, f64::max);
            if t.prov == Prov::Exact && t.rel > 0.25 && amax > 0.0 && pr.amp[j] / amax > 0.002 && t.n_on < n_leak {
                bucket[b].push(Need { label: format!("P{} stat+ fig5 ch{}+ {} (k={} of {})", b, [9, 8, 7, 6][k], d.sim.name, t.k_on, t.n_on), mass: d.sim.mass, z: d.sim.z, ccs: d.sim.ccs, vrf: pt5.vrf, frf: pt5.frf, vac: pt5.vac, fac: f * 1e3, dc: pt5.dc, n: n_leak });
            }
        }
    }
    // S14: every catalogued Part I peak at the S14 channel frequency (Part II shows nothing but 863 / 874 / 854 / 882 / 530 / 728)
    for (j, d) in ctx.defs14.iter().enumerate() { if ctx.src14.h[j] >= 1.5 || j < 3 { let b = if j < 3 { 2 } else { 4 }; bucket[b].extend(ctx.tm.needs(&d.sim, &pt14, f14, &format!("P{} S14", b), if j < 3 { ctx.tm.n_new } else { n_leak })); } }
    let needs: Vec<Need> = bucket.into_iter().flatten().collect();
    let needs = dedupe_needs(needs);
    let np = needs.iter().filter(|n| n.vac > 0.0).count();
    write_needs(&format!("{}/needed_jobs.tsv", ctx.out), &needs).map_err(|e| e.to_string())?;
    println!("\n[6] exact simulations still missing for this operating point: {} jobs ({} AC-on, {} AC-off baselines) x {} ions -> {}/needed_jobs.tsv", needs.len(), np, needs.len() - np, ctx.tm.n_new, ctx.out);
    let by_prov = |pv: Prov| ev5.iter().take(4).map(|e| e.pred.as_ref().unwrap().tv[0].iter().filter(|t| t.prov == pv).count()).sum::<usize>();
    println!("    species x channel cells: exact {}, interp {}, proxy {}, MISSING {} (of {})", by_prov(Prov::Exact), by_prov(Prov::Interp), by_prov(Prov::Proxy), by_prov(Prov::Missing), 4 * ctx.defs5.len());
    println!("    outputs in {}/: *_model_vs_data.csv, *.svg, species_source_5C.csv, transmission_fig5.csv, residual_features.txt, needed_jobs.tsv", ctx.out);
    Ok(())
}

// ---------------------------------------------------------------- GPU runner for the missing-job list (NOT run by the agent that wrote it)

fn run_jobs(args: &[String]) -> Result<(), String> {
    use crate::channel::{Engine, Fields, Geom, Job, Physics, RfAc, Species};
    let path = args.get(1).ok_or("usage: chsim spectrum-run <jobs.tsv>")?;
    let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut groups: Vec<(([f64; 4], f64, f64), Vec<Job>)> = vec![];
    for l in text.lines().filter(|l| !l.starts_with('#') && !l.trim().is_empty()) {
        let p: Vec<&str> = l.split('\t').collect();
        if p.len() < 14 { return Err(format!("bad job line: {}", l)); }
        let f = |i: usize| -> Result<f64, String> { p[i].parse::<f64>().map_err(|e| format!("{}: {}", p[i], e)) };
        let dc = [f(10)?, f(11)?, f(12)?, f(13)?];
        let job = Job { sp: Species { mass: f(1)?, z: f(2)?, ccs_n2: f(3)? }, drive: RfAc { v_rf: f(4)?, f_rf: f(5)?, v_ac: f(6)?, f_ac: f(7)? }, n: f(8)? as usize, seed: f(9)? as u64 };
        let key = (dc, job.drive.v_rf, job.drive.f_rf);
        match groups.iter_mut().find(|g| g.0 == key) { Some(g) => g.1.push(job), None => groups.push((key, vec![job])) }
    }
    let total: usize = groups.iter().map(|g| g.1.len()).sum();
    println!("{} jobs in {} (DC, RF) groups, {} ions", total, groups.len(), groups.iter().flat_map(|g| g.1.iter()).map(|j| j.n).sum::<usize>());
    if args.iter().any(|a| a == "--dry") { return Ok(()); }
    // same geometry as the R28 cache (reality/predict): s14_geom(code, gate 7 mm)
    let gargs: Vec<String> = vec!["".into(), "".into(), "".into(), "".into(), "".into(), "code".into(), "7.0".into()];
    let geom = crate::channel::s14_geom(&gargs); let _ = Geom::default;
    let eng = Engine::new(Fields::build(&geom));
    let mut ph = Physics::correct(0.1); ph.t_max = 3e-3;
    for ((dc, vrf, frf), jobs) in &groups { println!("running {} jobs at DC {:?}, RF {} V @ {} Hz", jobs.len(), dc, vrf, frf); let _ = crate::scan::cached_batch(&eng, jobs, *dc, &ph); }
    Ok(())
}
