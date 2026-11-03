//! Validation of the spectrum forward model: (1) synthetic round trip, (2) ranking of the R28 Fig. 5 grid by the forward-model
//! likelihood vs the old 12-scalar scoreboard (replicated here from the result cache, read-only, no GPU).

use crate::spectrum::*;
use crate::spectrum_detector::*;
use crate::spectrum_like::*;
use crate::spectrum_trans::*;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 { self.0 ^= self.0 << 13; self.0 ^= self.0 >> 7; self.0 ^= self.0 << 17; (self.0 >> 11) as f64 / (1u64 << 53) as f64 }
    fn normal(&mut self) -> f64 { let (a, b) = (self.next().max(1e-12), self.next()); (-2.0 * a.ln()).sqrt() * (2.0 * std::f64::consts::PI * b).cos() }
}

/// Render panel heights from amplitudes (no noise), exactly as the likelihood does.
fn render(p: &Panel, plac: &[Placed], pen: &Pen, amp: &[f64], delta: f64, normalise: bool) -> Vec<f64> {
    let resp: Vec<Response> = plac.iter().map(|s| response(s, &p.grid, pen, delta)).collect();
    let zeros = vec![0.0; plac.len()];
    let (h, _, _, _) = combine(&resp, amp, &zeros, p.grid.n);
    let m = h.iter().cloned().fold(0.0, f64::max);
    let s = if normalise && m > 0.0 { p.base_obs / m } else { 1.0 };
    h.iter().map(|v| v * s).collect()
}

/// Add measurement noise as the data have it: Gaussian (quantisation + relative) and censoring below the display floor.
fn observe(p: &Panel, h: &[f64], nz: &Noise, rng: &mut Rng) -> Panel {
    let sq = nz.sigq_px * p.pxh;
    let d: Vec<f64> = h.iter().map(|&m| { let v = m + rng.normal() * (sq * sq + (nz.rho_p * m).powi(2)).sqrt(); if v < p.floor { 0.0 } else { v } }).collect();
    p.with_heights(d)
}

pub fn selftest(args: &[String]) -> Result<(), String> {
    let cfg = match args.get(1) { Some(p) => Cfg::load(p)?, None => Cfg::empty() };
    let ctx = build_ctx(cfg)?;
    let mut rng = Rng(0x9E3779B97F4A7C15);
    let mut fails = 0;
    println!("SPECTRUM SELF-TEST (synthetic round trip; truth = the fitted 5C source, a made-up R=10 transmission, delta_true per panel)");

    // ---- truth
    let defs = &ctx.defs5;
    let plac = ctx.plac5();
    let h_true: Vec<f64> = ctx.src5.h.clone();
    let chan_mz = [952.64, 1071.6, 1224.55, 1428.47]; // channels 9+ 8+ 7+ 6+
    let eta = [0.10, 0.07, 0.05, 0.03];
    let t_true = |k: usize, j: usize| -> f64 {
        let mz = defs[j].sim.mz;
        let l = 1.0 / (1.0 + (2.0 * 10.0 * (mz / chan_mz[k] - 1.0)).powi(2));
        let mut t = eta[k] * l + 0.0005;
        if k == 3 && (mz - 402.6).abs() < 0.1 { t = 0.09; } // the 403 enrichment seen in 6+
        t
    };
    let amps: Vec<Vec<f64>> = (0..4).map(|k| (0..defs.len()).map(|j| h_true[j] * t_true(k, j)).collect()).collect();
    let delta_true = [("C", 0.15), ("A1", 0.45), ("A2", -0.30), ("A3", 0.60), ("A4", 0.10), ("B", -0.09)];

    // ---- 1. source (abundance) recovery from synthetic 5C
    let pc = ctx.panel("C");
    let dtc = delta_true[0].1;
    let hc = render(pc, &plac, &ctx.pen5, &h_true, dtc, false);
    let obs_c = observe(pc, &hc, &ctx.nz, &mut rng);
    let init: Vec<f64> = h_true.iter().map(|_| 1.0).collect();
    let fixed: Vec<Option<f64>> = defs.iter().enumerate().map(|(j, d)| if matches!(d.sim.kind, Kind::Chem(_)) && [328.8, 500.8, 574.6].iter().any(|m| (d.sim.mz - m).abs() < 0.01) { Some(h_true[j]) } else { None }).collect();
    let fit = fit_source(&obs_c, &plac, &ctx.pen5, &ctx.nz, &init, &fixed, 3);
    println!("\n[1] abundance recovery from a synthetic 5C (noise: {} px quantisation + {} % relative, floor {} %)", ctx.nz.sigq_px, 100.0 * ctx.nz.rho_p, pc.floor);
    println!("    fitted shift {:+.2} m/z (true {:+.2}; grid step {:.2})", fit.delta, dtc, 0.5 * pc.sig_cal);
    println!("    {:<26} {:>9} {:>9} {:>8}", "species", "H true", "H fit", "err");
    let mut worst: f64 = 0.0;
    for (j, d) in defs.iter().enumerate() {
        if h_true[j] < 5.0 * pc.floor { continue; }
        let e = fit.h[j] / h_true[j] - 1.0;
        worst = worst.max(e.abs());
        println!("    {:<26} {:>9.3} {:>9.3} {:>+7.1}%", d.sim.name, h_true[j], fit.h[j], 100.0 * e);
    }
    let ok1 = worst < 0.04 && (fit.delta - dtc).abs() <= 0.5 * pc.sig_cal + 1e-9;
    println!("    worst relative error over species >= 5x floor: {:.1} %   -> {}", 100.0 * worst, if ok1 { "PASS" } else { "FAIL" });
    if !ok1 { fails += 1; }

    // ---- 2. transmission recovery from synthetic channel panels (fit panel heights, compare ratios)
    println!("\n[2] transmission recovery: synthetic A1..A4 (base-peak normalised, own calibration shift, censored below floor) -> fit species heights -> T_k(j)/T_k(main)");
    let mains = [4usize, 3, 2, 1]; // indices of ubq 9+, 8+, 7+, 6+ in defs (index = z - 5)
    let mut worst_t: f64 = 0.0;
    for (k, nm) in ["A1", "A2", "A3", "A4"].iter().enumerate() {
        let p = ctx.panel(nm);
        let dt = delta_true[k + 1].1;
        let h = render(p, &plac, &ctx.pen5, &amps[k], dt, true);
        let obs = observe(p, &h, &ctx.nz, &mut rng);
        let init: Vec<f64> = vec![1.0; defs.len()];
        let fixed = vec![None; defs.len()];
        let f = fit_source(&obs, &plac, &ctx.pen5, &ctx.nz, &init, &fixed, 3);
        let main = mains[k];
        print!("    {} (shift fit {:+.2}, true {:+.2}):", nm, f.delta, dt);
        for (j, d) in defs.iter().enumerate() {
            let truth = amps[k][j] / amps[k][main];
            let rec = f.h[j] / f.h[main];
            // only species clearly above the display floor in this panel are testable
            if amps[k][j] / amps[k][main] * p.base_obs < 4.0 * p.floor { continue; }
            let e = rec / truth - 1.0;
            worst_t = worst_t.max(e.abs());
            print!("  {} {:.3}/{:.3}", d.sim.name.replace("ubq", "u").replace(" adduct", "+"), rec, truth);
        }
        println!();
    }
    let ok2 = worst_t < 0.10;
    println!("    worst relative error of T_k(j)/T_k(main) for species >= 4x floor: {:.1} %   -> {}", 100.0 * worst_t, if ok2 { "PASS" } else { "FAIL" });
    if !ok2 { fails += 1; }

    // ---- 3. the likelihood prefers the truth
    println!("\n[3] likelihood discrimination on synthetic panels at the true shift: truth vs perturbed transmission");
    let p = ctx.panel("A1");
    let h = render(p, &plac, &ctx.pen5, &amps[0], delta_true[1].1, true);
    let obs = observe(p, &h, &ctx.nz, &mut rng);
    let zeros = vec![0.0; defs.len()];
    let truth = eval_panel(&obs, &plac, &ctx.pen5, &amps[0], &zeros, true, &ctx.nz);
    let mut ok3 = true;
    for (label, j, f) in [("10+ leak x2", 5usize, 2.0), ("10+ leak x0.5", 5, 0.5), ("8+ leak x3", 3, 3.0)] {
        let mut a = amps[0].clone(); a[j] *= f;
        let r = eval_panel(&obs, &plac, &ctx.pen5, &a, &zeros, true, &ctx.nz);
        let ok = r.nll_marg > truth.nll_marg;
        ok3 &= ok;
        println!("    {:<16} -logL {:9.1}  vs truth {:9.1}  (+{:.1})  {}", label, r.nll_marg, truth.nll_marg, r.nll_marg - truth.nll_marg, if ok { "ok" } else { "WRONG" });
    }
    let per_col = truth.nll_marg / p.grid.n as f64;
    println!("    truth: -logL/column = {:.3} (expect ~ -ln(censored prob) ~ 0 for censored columns, ~ +ln(sigma)+1.42 for Gaussian columns)", per_col);
    if !ok3 { fails += 1; }
    println!("    -> {}", if ok3 { "PASS" } else { "FAIL" });

    // ---- 4. detector-resolution sensitivity: heights (data units) must not move; abundances (area units) do
    println!("\n[4] Orbitrap resolution sensitivity: apex gain of 8+ and 12+ vs R(200); the fitted H (height units) is independent of R by construction");
    let ubq = envelope_ubiquitin();
    for r200 in [17_500.0, 35_000.0, 70_000.0, 140_000.0] {
        let d = Orbitrap { r200 };
        let (l8, l12, l5) = (detector_lines(&ubq, 8.0, &d), detector_lines(&ubq, 12.0, &d), detector_lines(&ubq, 5.0, &d));
        println!("    R200 {:>7.0}: FWHM(1071.6) {:.3} m/z, lines(8+) {:>2}, gain 8+ {:.2}, 12+/8+ {:.2}, 5+/8+ {:.2}", r200, d.fwhm(1071.6), l8.lines.len(), l8.apex_gain, l12.apex_gain / l8.apex_gain, l5.apex_gain / l8.apex_gain);
    }
    println!("\nSELF-TEST {}", if fails == 0 { "PASSED" } else { "FAILED" });
    if fails > 0 { Err(format!("{} self-test sections failed", fails)) } else { Ok(()) }
}

// ---------------------------------------------------------------- old scoreboard replicated from the cache

pub struct SbRow { pub name: &'static str, pub model: f64, pub meas: f64, pub sigma: f64, pub msig: f64, pub kind: u8, pub pull: f64 }
pub struct Sb { pub loss: f64, pub rows: Vec<SbRow>, pub dead: usize }

fn ratio(a: f64, b: f64, n: usize) -> (f64, f64) {
    let nf = n as f64;
    if b * nf < 10.0 { return (if b > 0.0 { a / b } else { f64::NAN }, f64::NAN); }
    let sa = (a.max(1.0 / nf) * (1.0 - a).max(0.0) / nf).sqrt(); let sb = (b * (1.0 - b) / nf).sqrt();
    let r = a / b; (r, r * ((sa / a.max(1.0 / nf)).powi(2) + (sb / b).powi(2)).sqrt())
}

/// `chsim reality` (REALITY_ONLY=fig5) re-evaluated from the cache: same observables, targets and loss. None if the cache is incomplete.
pub fn scoreboard_fig5(tc: &TCache, vrf: f64, vac: f64, gate: f64, n: usize) -> Option<Sb> {
    let dc = [18.0, gate, 15.0, 20.0];
    let sp = |z: u32| SimSp::new("", Kind::Ubq(z), 0.0);
    let frac = |c: Counts| c.transmitted as f64 / c.n as f64;
    let mut base = vec![];
    for z in 6..=10u32 { let s = sp(z); base.push(frac(tc.exact_seed(s.mass, s.z, s.ccs, vrf, 0.0, Some(50e3), &dc, 50 + z as u64)?)); }
    let mut fch = vec![];
    for z in [6u32, 7, 8, 9] {
        let s = sp(z);
        let scan = tc.adaptive_scan_pts(s.mass, s.z, s.ccs, vrf, vac, &dc, 2000 * z as u64);
        if scan.is_empty() { return None; }
        fch.push(scan.iter().fold((0.0, -1.0), |a, (f, c)| if frac(*c) > a.1 { (*f, frac(*c)) } else { a }).0);
    }
    let mut mtx = vec![];
    for ci in 0..4 {
        let mut row = vec![];
        let zc = [6u32, 7, 8, 9][ci];
        for (i, z) in [6u32, 7, 8, 9, 10].iter().enumerate() { let s = sp(*z); row.push((frac(tc.exact_seed(s.mass, s.z, s.ccs, vrf, vac, Some(fch[ci] * 1e3), &dc, 7000 + 10 * zc as u64 + *z as u64)?) - base[i]).max(0.0)); }
        mtx.push(row);
    }
    let nn = n as f64;
    let mut rows = vec![];
    rows.push(SbRow { name: "Fig5 channel 6+ frequency [kHz]", model: fch[0], meas: 40.0, sigma: 3.0, kind: 0, msig: if mtx[0][0] * nn < 10.0 { f64::NAN } else { 0.5 }, pull: 0.0 });
    rows.push(SbRow { name: "Fig5 channel 9+ frequency [kHz]", model: fch[3], meas: 60.0, sigma: 3.0, kind: 0, msig: if mtx[3][3] * nn < 10.0 { f64::NAN } else { 0.5 }, pull: 0.0 });
    for (name, c, meas, sig) in [("Fig5 efficiency 6+ / 8+", 0usize, 0.21, 0.07), ("Fig5 efficiency 7+ / 8+", 1, 1.15, 0.35), ("Fig5 efficiency 9+ / 8+", 3, 2.5, 0.75)] {
        let (v, e) = ratio(mtx[c][c], mtx[2][2], n); rows.push(SbRow { name, model: v, meas, sigma: sig, kind: 0, msig: e, pull: 0.0 });
    }
    let (v, e) = ratio(mtx[3][4], mtx[3][3], n); rows.push(SbRow { name: "Fig5 ch9+: T(10+)/T(9+)", model: v, meas: 0.272, sigma: 0.046, kind: 0, msig: e, pull: 0.0 });
    let (v, e) = ratio(mtx[3][2], mtx[3][3], n); rows.push(SbRow { name: "Fig5 ch9+: T(8+)/T(9+)", model: v, meas: 0.0030, sigma: 0.0012, kind: 0, msig: e, pull: 0.0 });
    let (v, e) = ratio(mtx[0][1], mtx[0][0], n); rows.push(SbRow { name: "Fig5 ch6+: T(7+)/T(6+)", model: v, meas: 0.0031, sigma: 0.0031, kind: 1, msig: e, pull: 0.0 });
    let (mut loss, mut dead) = (0.0, 0);
    for r in rows.iter_mut() {
        if !r.msig.is_finite() { dead += 1; continue; }
        let s = (r.sigma * r.sigma + r.msig * r.msig).sqrt();
        r.pull = if r.kind == 1 { ((r.model - r.meas) / s).max(0.0) } else { (r.model - r.meas) / s };
        loss += r.pull * r.pull;
    }
    loss += 9.0 * dead as f64;
    Some(Sb { loss, rows, dead })
}

fn ranks(v: &[f64]) -> Vec<f64> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap());
    let mut r = vec![0.0; v.len()];
    let mut i = 0;
    while i < idx.len() {
        let mut j = i; while j + 1 < idx.len() && v[idx[j + 1]] == v[idx[i]] { j += 1; }
        for k in i..=j { r[idx[k]] = 0.5 * (i + j) as f64 + 1.0; }
        i = j + 1;
    }
    r
}
fn spearman(a: &[f64], b: &[f64]) -> f64 {
    let (ra, rb) = (ranks(a), ranks(b));
    let n = a.len() as f64;
    let (ma, mb) = (ra.iter().sum::<f64>() / n, rb.iter().sum::<f64>() / n);
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for i in 0..a.len() { sab += (ra[i] - ma) * (rb[i] - mb); saa += (ra[i] - ma).powi(2); sbb += (rb[i] - mb).powi(2); }
    sab / (saa * sbb).sqrt()
}

pub fn grid(args: &[String]) -> Result<(), String> {
    let cfgp = args.get(1).ok_or("usage: chsim spectrum-grid <config>")?;
    let mut ctx = build_ctx(Cfg::load(cfgp)?)?;
    ctx.mask_chem = ctx.cfg.f("grid_mask_chem", 1.0) > 0.5;
    let vrf = ctx.cfg.f("vrf5", 56.0);
    let n_score = ctx.cfg.f("n_score", 3000.0) as usize;
    let vacs = [0.2, 0.27, 0.35];
    let gates = [24.0, 24.5, 25.0];
    println!("R28 Fig. 5 GRID: RF {} V, AC +-{{0.2,0.27,0.35}} V x gate {{24,24.5,25}} V, code geometry. forward-model -logL (A1-A4 + B, calibration shift marginalised) vs the 12-scalar scoreboard (Fig. 5 rows only, post-R30 targets)", vrf);
    println!("(the scoreboard is replicated from the cache; at +-0.2 V / 24 V R30 reports loss 41.2)\n");
    println!("{:>6} {:>6} | {:>8} {:>5} | {:>12} {:>12} | {:>9} {:>9} {:>9}  channels [kHz]", "AC[V]", "gate", "scoreb", "dead", "-logL (on)", "-logL (sub)", "ch9 T10/9", "ch9 T8/9", "ch6 T7/6");
    let mut pts: Vec<(f64, f64, f64, f64, f64)> = vec![];
    let mut incomplete = vec![];
    let mut per_panel: Vec<(f64, f64, Vec<f64>)> = vec![];
    for &vac in &vacs {
        for &gate in &gates {
            let sb = match scoreboard_fig5(&ctx.tm.cache, vrf, vac, gate, n_score) { Some(s) => s, None => { incomplete.push((vac, gate)); continue; } };
            let pt = Point { vrf, frf: 500e3, vac, dc: [18.0, gate, 15.0, 20.0] };
            let freqs = channel_freqs(&ctx, &pt);
            let mut nll = [0.0f64; 2];
            for (m, mode) in [TMode::On, TMode::Sub].iter().enumerate() {
                ctx.tm.mode = *mode;
                let ev = eval_fig5(&ctx, &pt, &freqs);
                nll[m] = ev.iter().map(|e| e.res.nll_marg).sum();
                if m == 0 { per_panel.push((vac, gate, ev.iter().map(|e| e.res.nll_marg).collect::<Vec<f64>>())); }
            }
            ctx.tm.mode = TMode::parse(&ctx.cfg.s("tmode", "on"))?;
            let g = |name: &str| sb.rows.iter().find(|r| r.name.contains(name)).map(|r| r.model).unwrap_or(f64::NAN);
            println!("{:>6.2} {:>6.1} | {:>8.1} {:>5} | {:>12.1} {:>12.1} | {:>9.4} {:>9.4} {:>9.4}  {:?}", vac, gate, sb.loss, sb.dead, nll[0], nll[1], g("T(10+)"), g("T(8+)"), g("T(7+)"), freqs);
            pts.push((vac, gate, sb.loss, nll[0], nll[1]));
        }
    }
    println!("\nper-panel -logL (on mode) at each grid point:   A1 (9+)  A2 (8+)  A3 (7+)  A4 (6+)  B (all)");
    for (v, g, p) in &per_panel { println!("  AC {:.2} V gate {:4.1} V:  {}", v, g, p.iter().map(|x| format!("{:8.0}", x)).collect::<Vec<_>>().join(" ")); }
    for (v, g) in &incomplete { println!("  (AC {} V, gate {} V: cache incomplete for the scoreboard observables, skipped)", v, g); }
    if pts.len() >= 3 {
        let sbv: Vec<f64> = pts.iter().map(|p| p.2).collect();
        let on: Vec<f64> = pts.iter().map(|p| p.3).collect();
        let sub: Vec<f64> = pts.iter().map(|p| p.4).collect();
        println!("\nSpearman rank correlation over {} grid points: scoreboard vs forward(on) {:+.2}, scoreboard vs forward(sub) {:+.2}, forward(on) vs forward(sub) {:+.2}", pts.len(), spearman(&sbv, &on), spearman(&sbv, &sub), spearman(&on, &sub));
        let best = |v: &[f64]| { let i = v.iter().enumerate().fold(0, |a, (i, &x)| if x < v[a] { i } else { a }); (pts[i].0, pts[i].1) };
        println!("best point: scoreboard AC +-{} V gate {} V | forward(on) AC +-{} V gate {} V | forward(sub) AC +-{} V gate {} V", best(&sbv).0, best(&sbv).1, best(&on).0, best(&on).1, best(&sub).0, best(&sub).1);
        let (rs, ro): (Vec<f64>, Vec<f64>) = (ranks(&sbv), ranks(&on));
        println!("ranks (1 = best):   {}", pts.iter().enumerate().map(|(i, p)| format!("({:.2},{:.1}) sb {} fwd {}", p.0, p.1, rs[i], ro[i])).collect::<Vec<_>>().join(" | "));
    }
    Ok(())
}

/// Is the Orbitrap resolution identifiable from the published pixels? Refit sources and re-evaluate S14 / Fig. 5 for several R(200).
pub fn res_scan(args: &[String]) -> Result<(), String> {
    let cfgp = args.get(1).ok_or("usage: chsim spectrum-res <config>")?;
    println!("Orbitrap resolution scan: R(200) -> detector lines of the S14 main ion (Bak BH3 2+, 863.47), 5C / S14 source fits and the S14 Part II likelihood");
    println!("{:>9} {:>11} {:>8} {:>10} {:>12} {:>12} {:>12} {:>14}", "R(200)", "FWHM(863)", "lines", "gain Bak", "-logL 5C", "-logL S14I", "-logL S14II", "S14: Na/main");
    for r in [300.0, 1000.0, 3000.0, 10000.0, 17500.0, 35000.0, 70000.0, 140000.0] {
        let mut cfg = Cfg::load(cfgp)?;
        cfg.0.insert("r200".into(), format!("{}", r));
        let ctx = build_ctx(cfg)?;
        let (pt14, f14) = { let p = s14_point(&ctx); let f = s14_freq(&ctx, &p); (p, f) };
        let e = eval_s14(&ctx, &pt14, f14);
        let bak = &ctx.defs14[0];
        let pr = e.pred.as_ref().unwrap();
        let na = pr.amp[1] / pr.amp[0];
        println!("{:>9.0} {:>11.4} {:>8} {:>10.2} {:>12.1} {:>12.1} {:>12.1} {:>14.3}", r, ctx.det.fwhm(863.47), bak.plac.lines.len(), bak.gain, ctx.src5.nll_marg, ctx.src14.nll_marg, e.res.nll_marg, na);
    }
    Ok(())
}
