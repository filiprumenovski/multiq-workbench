//! Reality: the instrument's published numbers (digitised from their figures, with the uncertainties stated in
//! data/their_side/figure_data) against the live model, with pulls in sigma. Mirrors `chsim reality` (marks, pull sign).
//!
//! Everything labelled DIGITISED comes from a figure raster (pixel quantisation errors, not statistical CIs); everything
//! labelled MODEL is computed live by this program.

use crate::engine::SpStats;
use crate::fit::ResFit;
use crate::stats::{self, Frac, Settle, Val};
use crate::world::Sp;
use std::path::{Path, PathBuf};

/// Their S14 resonance (fig. S14 caption: AC ON 83 kHz) and the reading uncertainty (task spec and `chsim reality`).
pub const S14_F_RES: Val = Val { v: 83.0, e: 2.0 };
/// Caption: "90-100 %" efficiency in fig. S14.
pub const S14_CAPTION_BAND: (f64, f64) = (0.90, 1.00);

#[derive(Clone, Debug)]
pub struct S14Data {
    /// Part I (gate open) Na-adduct / main apex-height ratio.
    pub adduct_over_main_open: Val,
    /// Part II (AC on, 83 kHz) ratio.
    pub adduct_over_main_on: Val,
    /// Ratio of ratios II/I = T(Na)/T(main) through the AC-on channel.
    pub ratio_of_ratios: Val,
    /// Main-ion counts AC on / gate open (apex height x NL).
    pub main_counts_on_over_open: Val,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum XtKind { Measured, Marginal, UpperLimit, Unconstrained, Other }

#[derive(Clone, Debug)]
pub struct XtRow {
    pub channel: u32,
    pub nb: u32,
    pub kind: XtKind,
    pub value: Option<f64>,
    pub lo: Option<f64>,
    pub hi: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct RealityData {
    pub s14: S14Data,
    pub xt: Vec<XtRow>,
}

fn csv_split(line: &str) -> Vec<String> {
    let (mut out, mut cur, mut q) = (vec![], String::new(), false);
    for c in line.chars() {
        match c { '"' => q = !q, ',' if !q => { out.push(std::mem::take(&mut cur)); } c => cur.push(c) }
    }
    out.push(cur);
    out
}

fn read_rows(path: &Path) -> Result<Vec<Vec<String>>, String> {
    let t = std::fs::read_to_string(path).map_err(|e| format!("{}: {}", path.display(), e))?;
    Ok(t.lines().filter(|l| !l.starts_with('#') && !l.trim().is_empty()).map(csv_split).collect())
}

/// data/their_side/figure_data next to the crate (or $MQVIZ_DATA).
pub fn find_data_dir() -> Option<PathBuf> {
    if let Some(d) = std::env::var_os("MQVIZ_DATA") { let p = PathBuf::from(d); if p.is_dir() { return Some(p); } }
    let man = Path::new(env!("CARGO_MANIFEST_DIR"));
    [man.join("../data/their_side/figure_data"), PathBuf::from("data/their_side/figure_data"), PathBuf::from("../data/their_side/figure_data")]
        .into_iter().find(|p| p.join("fig5A_crosstalk_matrix.csv").exists())
}

impl RealityData {
    pub fn load(dir: &Path) -> Result<RealityData, String> {
        // ---- S14
        let rows = read_rows(&dir.join("figS14_adduct_ratio.csv"))?;
        let get = |name: &str| -> Result<Val, String> {
            let r = rows.iter().find(|r| r.first().map_or(false, |c| c.starts_with(name))).ok_or(format!("figS14_adduct_ratio.csv: no row {name}"))?;
            Ok(Val::new(r[1].parse().map_err(|_| format!("bad value in {name}"))?, r[2].parse().map_err(|_| format!("bad uncertainty in {name}"))?))
        };
        let s14 = S14Data {
            adduct_over_main_open: get("Part_I_adduct_over_main")?,
            adduct_over_main_on: get("Part_II_adduct_over_main")?,
            ratio_of_ratios: get("ratio_of_ratios_II_over_I")?,
            main_counts_on_over_open: get("main_ion_counts_AC_on_over_open")?,
        };
        // ---- Fig. 5A cross-talk matrix
        let rows = read_rows(&dir.join("fig5A_crosstalk_matrix.csv"))?;
        let mut xt = vec![];
        for r in rows.iter().skip(1) {
            if r.len() < 13 { continue; }
            let num = |s: &str| s.trim().parse::<f64>().ok();
            let kind = match r[12].trim() { "measured" => XtKind::Measured, "measured_marginal" => XtKind::Marginal, "upper-limit" => XtKind::UpperLimit, "unconstrained" => XtKind::Unconstrained, _ => XtKind::Other };
            let (Ok(channel), Ok(nb)) = (r[0].trim().parse::<u32>(), r[1].trim().parse::<u32>()) else { continue };
            xt.push(XtRow { channel, nb, kind, value: num(&r[9]), lo: num(&r[10]), hi: num(&r[11]) });
        }
        if xt.is_empty() { return Err("fig5A_crosstalk_matrix.csv: no rows".into()); }
        Ok(RealityData { s14, xt })
    }
    pub fn xt(&self, channel: u32, nb: u32) -> Option<&XtRow> { self.xt.iter().find(|r| r.channel == channel && r.nb == nb) }
}

/// Calibrated model channel frequencies at RF 56 V (R28 resonance maxima, SPECTRUM.md): 6+ 39.0, 7+ 47.0, 8+ 53.0, 9+ 59.25 kHz.
pub const CHANNELS_56V: [(u32, f64); 4] = [(6, 39.0), (7, 47.0), (8, 53.0), (9, 59.25)];

/// Channel (charge state, nominal frequency, distance) nearest to `f_ac`; frequencies scale with the RF amplitude (low-q, f ~ V_rf).
pub fn nearest_channel(f_ac: f64, v_rf: f64) -> (u32, f64, f64) {
    let k = v_rf / 56.0;
    CHANNELS_56V.iter().map(|&(z, f)| (z, f * k, (f * k - f_ac).abs())).min_by(|a, b| a.2.partial_cmp(&b.2).unwrap()).unwrap()
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RowKind { Measured, Limit, Info }

#[derive(Clone, Debug)]
pub struct RealRow {
    pub name: String,
    pub ours: Val,
    pub ours_txt: String,
    pub theirs: Option<Val>,
    pub theirs_txt: String,
    pub kind: RowKind,
    pub pull: f64,
    pub settle: Option<Settle>,
    /// Where their number comes from.
    pub src: &'static str,
}

impl RealRow {
    pub fn mark(&self) -> &'static str { if self.kind == RowKind::Info { "" } else { stats::mark(self.pull) } }
}

fn fmt_val(v: Val, dec: usize) -> String { if v.ok() { format!("{:.*} ± {:.*}", dec, v.v, dec, v.e) } else if v.v.is_finite() { format!("{:.*}", dec, v.v) } else { "—".into() } }
fn fmt_sig(x: f64) -> String { if !x.is_finite() { "—".into() } else if x == 0.0 { "0".into() } else if x.abs() >= 0.1 { format!("{:.3}", x) } else if x.abs() >= 0.01 { format!("{:.4}", x) } else { format!("{:.5}", x) } }
/// two significant digits
fn sig2(x: f64) -> String { if !x.is_finite() { "—".into() } else if x == 0.0 { "0".into() } else { format!("{:.*}", ((1 - x.abs().log10().floor() as i32).max(0)) as usize, x) } }
fn fmt_sigv(v: Val) -> String { if v.ok() { format!("{} ± {}", fmt_sig(v.v), fmt_sig(v.e)) } else { "—".into() } }

fn worst(a: Settle, b: Settle) -> Settle {
    let rank = |s: &Settle| match s { Settle::Settled { .. } => 0, Settle::Accumulating { .. } => 1, Settle::Drifting { .. } => 2 };
    if rank(&a) >= rank(&b) { a } else { b }
}

/// Per-pass transmission of species `s`: T out of every finished ion since the epoch start (and its two halves).
fn tfrac(st: &SpStats) -> (Frac, Settle) {
    let f = SpStats::t_pass(&st.ep);
    (f, stats::settle_frac(SpStats::t_pass(&st.half_a), SpStats::t_pass(&st.half_b), 10))
}

/// Efficiency T/(T+C).
fn eff(st: &SpStats) -> (Frac, Settle) {
    let f = SpStats::eff(&st.ep);
    (f, stats::settle_frac(SpStats::eff(&st.half_a), SpStats::eff(&st.half_b), 10))
}

/// T_a / T_b with its convergence (worst of the two).
pub fn t_ratio(a: &SpStats, b: &SpStats) -> (Val, Frac, Frac, Settle) {
    let (fa, sa) = tfrac(a);
    let (fb, sb) = tfrac(b);
    (stats::ratio(fa, fb), fa, fb, worst(sa, sb))
}

/// Live S14 rows. Species 0 = main Bak BH3 2+, species 1 = Na adduct. `fit` is the latest sweep fit of the main ion.
pub fn s14_rows(d: &RealityData, sp: &[Sp], st: &[SpStats], fit: Option<&ResFit>) -> Vec<RealRow> {
    let mut rows = vec![];
    // resonance
    let (ours, txt) = match fit { Some(f) => (Val::new(f.p[2], f.e[2]), format!("{:.2} ± {:.2} kHz ({} fit, χ²/dof {:.1})", f.p[2], f.e[2], f.shape.name(), f.chi2_red())), None => (Val::nan(), "run a sweep".into()) };
    rows.push(RealRow { name: "resonance f₀ [kHz]".into(), ours, ours_txt: txt, theirs: Some(S14_F_RES), theirs_txt: "83 ± 2".into(), kind: RowKind::Measured,
                        pull: if ours.ok() { stats::pull(ours, S14_F_RES) } else { f64::NAN }, settle: None, src: "fig. S14 caption 83 kHz; ±2 kHz reading uncertainty" });
    if sp.len() < 2 || st.len() < 2 { return rows; }
    // composition of the source in the preset (an input, shown for transparency)
    let (w0, w1) = (sp[0].weight, sp[1].weight);
    let src_in = Val::new(w1 / w0, f64::NAN);
    rows.push(RealRow { name: "Na/main in the source (preset input)".into(), ours: src_in, ours_txt: format!("{:.4} (set by the preset)", src_in.v), theirs: Some(d.s14.adduct_over_main_open), theirs_txt: fmt_val(d.s14.adduct_over_main_open, 4), kind: RowKind::Info,
                        pull: f64::NAN, settle: None, src: "DIGITISED fig. S14 Part I (gate open); the preset uses 5 % Na" });
    // Na adduct transmission relative to main ion
    let (r, fa, fb, settle) = t_ratio(&st[1], &st[0]);
    let th = d.s14.ratio_of_ratios;
    rows.push(RealRow { name: "Na adduct / main  T ratio".into(), ours: r, ours_txt: fmt_val(r, 3), theirs: Some(th), theirs_txt: fmt_val(th, 2), kind: RowKind::Measured,
                        pull: if r.ok() && fb.k >= 10 { stats::pull(r, th) } else { f64::NAN }, settle: Some(settle), src: "DIGITISED fig. S14 II/I ratio of ratios" });
    // output ratio in the AC-on spectrum = digitised source ratio x model T ratio
    let src = d.s14.adduct_over_main_open;
    let out = if r.ok() { Val::new(src.v * r.v, src.v * r.v * ((src.e / src.v).powi(2) + (r.e / r.v.max(1e-300)).powi(2)).sqrt()) } else { Val::nan() };
    let thr = d.s14.adduct_over_main_on;
    rows.push(RealRow { name: "Na/main in AC-on spectrum".into(), ours: out, ours_txt: fmt_val(out, 4), theirs: Some(thr), theirs_txt: fmt_val(thr, 4), kind: RowKind::Measured,
                        pull: if out.ok() && fb.k >= 10 { stats::pull(out, thr) } else { f64::NAN }, settle: Some(settle), src: "DIGITISED fig. S14 II; model = 0.0893 ± 0.0223 (Part I) × model T ratio" });
    let _ = fa;
    // steady-state efficiency of the main ion
    let (e, es) = eff(&st[0]);
    let ev = Val::from(e);
    let th = d.s14.main_counts_on_over_open;
    rows.push(RealRow { name: "efficiency T/(T+C), main ion".into(), ours: ev, ours_txt: fmt_val(ev, 3), theirs: Some(th), theirs_txt: fmt_val(th, 3), kind: RowKind::Measured,
                        pull: if e.k + 0 >= 1 && e.n >= 100 { stats::pull(ev, th) } else { f64::NAN }, settle: Some(es), src: "DIGITISED fig. S14: main-ion counts AC-on/open; model error is statistical only" });
    // caption band
    let below = if e.n >= 100 && ev.e > 0.0 { ((S14_CAPTION_BAND.0 - ev.v) / ev.e).max(0.0) } else { f64::NAN };
    rows.push(RealRow { name: "efficiency within 90–100 %".into(), ours: ev, ours_txt: fmt_val(ev, 3), theirs: Some(Val::new(0.95, 0.05)), theirs_txt: "0.90 – 1.00".into(), kind: RowKind::Limit,
                        pull: below, settle: Some(es), src: "fig. S14 caption" });
    rows
}

/// One cell of the live cross-talk row.
#[derive(Clone, Debug)]
pub struct XtCell {
    pub nb: u32,
    pub ours: Val,
    pub k_num: u64,
    pub ul95: f64,
    pub settle: Settle,
    pub theirs: Option<XtRow>,
    pub pull: f64,
    pub verdict: &'static str,
}

/// Live row of the cross-talk matrix for `channel_z`: T(neighbour)/T(channel ion) per pass, with pulls against their matrix.
pub fn xtalk_row(d: &RealityData, sp: &[Sp], st: &[SpStats], channel_z: u32) -> Vec<XtCell> {
    let Some(ci) = sp.iter().position(|s| s.z as u32 == channel_z) else { return vec![] };
    let mut out = vec![];
    for (j, s) in sp.iter().enumerate() {
        let z = s.z as u32;
        if j == ci || st.len() <= j { continue; }
        let (r, fa, fb, settle) = t_ratio(&st[j], &st[ci]);
        let ul = stats::ratio_ul95(fa, fb);
        let theirs = d.xt(channel_z, z).cloned();
        let mut pull = f64::NAN;
        let mut verdict = "";
        if fb.k < 10 { verdict = "no signal in channel"; }
        else if let Some(t) = &theirs {
            match t.kind {
                XtKind::Measured | XtKind::Marginal => {
                    if let (Some(v), Some(lo), Some(hi)) = (t.value, t.lo, t.hi) { pull = stats::pull(r, Val::new(v, 0.5 * (hi - lo))); }
                }
                XtKind::UpperLimit => {
                    let u = t.hi.unwrap_or(f64::NAN);
                    if fa.k == 0 { if ul <= u { pull = 0.0; } else { verdict = "need more ions"; } } else { pull = ((r.v - u) / r.e).max(0.0); }
                }
                _ => {}
            }
        }
        out.push(XtCell { nb: z, ours: r, k_num: fa.k, ul95: ul, settle, theirs, pull, verdict });
    }
    out.sort_by_key(|c| c.nb);
    out
}

pub fn theirs_text(t: &XtRow) -> String {
    match t.kind {
        XtKind::Measured | XtKind::Marginal => match (t.value, t.lo, t.hi) {
            (Some(v), Some(lo), Some(hi)) => format!("{} [{}–{}]{}", fmt_sig(v), sig2(lo), sig2(hi), if t.kind == XtKind::Marginal { " marginal" } else { "" }),
            _ => "—".into(),
        },
        XtKind::UpperLimit => format!("< {}", fmt_sig(t.hi.unwrap_or(f64::NAN))),
        XtKind::Unconstrained => "unconstrained".into(),
        XtKind::Other => "—".into(),
    }
}

pub fn ours_text(c: &XtCell) -> String {
    let t = if c.k_num == 0 { format!("0 / <{} (95%)", fmt_sig(c.ul95)) } else { fmt_sigv(c.ours) };
    if c.verdict.is_empty() { t } else { format!("{} · {}", t, c.verdict) }
}

pub fn fig5_rows(d: &RealityData, sp: &[Sp], st: &[SpStats], channel_z: u32) -> Vec<RealRow> {
    xtalk_row(d, sp, st, channel_z).into_iter().map(|c| {
        let (tt, kind, th) = match &c.theirs {
            Some(t) => (theirs_text(t), if t.kind == XtKind::UpperLimit { RowKind::Limit } else if matches!(t.kind, XtKind::Measured | XtKind::Marginal) { RowKind::Measured } else { RowKind::Info },
                        match (t.kind, t.value, t.lo, t.hi) { (XtKind::Measured | XtKind::Marginal, Some(v), Some(lo), Some(hi)) => Some(Val::new(v, 0.5 * (hi - lo))), _ => None }),
            None => ("—".into(), RowKind::Info, None),
        };
        RealRow { name: format!("ch {}+ : T({}+)/T({}+)", channel_z, c.nb, channel_z), ours: c.ours, ours_txt: ours_text(&c), theirs: th, theirs_txt: tt, kind, pull: c.pull, settle: Some(c.settle), src: "DIGITISED fig. 5A/5C, ±1 px bounds" }
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn loads_published_numbers_when_the_data_is_there() {
        let Some(dir) = find_data_dir() else { return };
        let d = RealityData::load(&dir).unwrap();
        assert!((d.s14.ratio_of_ratios.v - 0.64).abs() < 1e-9 && (d.s14.ratio_of_ratios.e - 0.186).abs() < 1e-9);
        assert!((d.s14.adduct_over_main_on.v - 0.0571).abs() < 1e-9 && (d.s14.adduct_over_main_on.e - 0.0085).abs() < 1e-9);
        assert!((d.s14.adduct_over_main_open.v - 0.0893).abs() < 1e-9);
        assert!((d.s14.main_counts_on_over_open.v - 1.016).abs() < 1e-9 && (d.s14.main_counts_on_over_open.e - 0.012).abs() < 1e-9);
        let t = d.xt(9, 10).unwrap();
        assert_eq!(t.kind, XtKind::Measured);
        assert!((t.value.unwrap() - 0.2719).abs() < 1e-9 && (t.lo.unwrap() - 0.2294).abs() < 1e-9 && (t.hi.unwrap() - 0.3209).abs() < 1e-9);
        assert_eq!(d.xt(6, 7).unwrap().kind, XtKind::UpperLimit);
        assert!((d.xt(6, 7).unwrap().hi.unwrap() - 0.003139).abs() < 1e-9);
        assert_eq!(d.xt(9, 5).unwrap().kind, XtKind::Unconstrained);
        assert_eq!(d.xt.len(), 24);
    }
    fn sp_of(z: f64, label: &str) -> Sp { Sp { label: label.into(), mass: 8564.8 + z * 1.00728, z, ccs: 2000.0, weight: 1.0 } }
    fn st_of(t: u64, c: u64, rr: u64, to: u64) -> SpStats { SpStats { ep: [t, c, rr, to], half_a: [t / 2, c / 2, rr / 2, to / 2], half_b: [t - t / 2, c - c / 2, rr - rr / 2, to - to / 2], ..Default::default() } }

    #[test]
    fn xtalk_pulls_follow_the_conventions() {
        let Some(dir) = find_data_dir() else { return };
        let d = RealityData::load(&dir).unwrap();
        let sp = vec![sp_of(6.0, "6+"), sp_of(7.0, "7+"), sp_of(8.0, "8+"), sp_of(9.0, "9+"), sp_of(10.0, "10+")];
        // channel 9+: T = 6 % per pass; 6+ and 7+ have no transmitted ions; 8+ at 0.3 %; 10+ at 27 % of the channel ion
        let st = vec![st_of(0, 500, 500, 99_000), st_of(0, 500, 500, 9_000), st_of(60, 500, 500, 19_000), st_of(1200, 5000, 6000, 7800), st_of(325, 2000, 2000, 15_675)];
        let cells = xtalk_row(&d, &sp, &st, 9);
        let by = |z: u32| cells.iter().find(|c| c.nb == z).unwrap();
        // 6+: no counts, n = 100k: 95 % limit 3/1e5/0.06 = 5e-4 < their 0.0124 -> consistent, pull 0
        assert!(by(6).k_num == 0 && by(6).ul95 < 0.0124 && by(6).pull == 0.0, "{:?}", by(6));
        // 7+: n = 10k: limit 3/1e4/0.06 = 0.005 > their 0.0032 limit -> cannot say yet
        assert!(by(7).k_num == 0 && by(7).ul95 > 0.0032 && by(7).pull.is_nan(), "{:?}", by(7));
        // 8+: 60 / 20k = 0.003, ratio 0.05 vs their 0.00299 -> large positive pull
        assert!(by(8).pull > 5.0, "{:?}", by(8));
        // 10+: 325/16.3k = 2 % -> ratio 0.33 vs 0.272 +- 0.046 -> small pull
        assert!(by(10).pull.abs() < 2.0, "{:?}", by(10));
    }

    #[test]
    fn channels_follow_rf() {
        assert_eq!(nearest_channel(59.0, 56.0).0, 9);
        assert_eq!(nearest_channel(40.0, 56.0).0, 6);
        assert_eq!(nearest_channel(49.0, 56.0).0, 7);
        assert_eq!(nearest_channel(73.0, 70.0).0, 9); // 59.25 * 1.25 = 74.06
    }
    #[test]
    fn csv_quotes() {
        assert_eq!(csv_split("a,\"b,c\",d"), vec!["a", "b,c", "d"]);
    }
}
