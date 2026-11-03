//! Predicted single-channel spectrum (Fig. 5A) from the live per-species transmission, using chsim's spectrum forward model
//! (chsim/SPECTRUM.md): source abundances from the Fig. 5C fit x live per-pass transmission -> isotope/Orbitrap line ->
//! pen kernel (how the figure draws a line) -> base-peak normalisation on the panel's pixel grid.
//!
//! What is MODEL: the source heights H_j (fitted to the digitised 5C control), the live T_j (this program), the pen kernel
//! and the calibration shift (fitted to the digitised panels). What is DIGITISED: the panel pixels (data).
//! Not simulated here: ubiquitin 5+, 11+, 12+ and the chemical background ions (T = 0, listed in `missing`); same-charge
//! adducts take the parent charge state's T (`proxied`, SPECTRUM.md calls this an approximation).

use chsim::spectrum::{fig5_species, fit_5c, SpDef};
use chsim::spectrum_detector::{combine, response, Orbitrap, Pen, Response};
use chsim::spectrum_like::{eval_panel, fit_pen, load_fig5, Noise, Panel};
use chsim::spectrum_trans::Kind;
use crate::stats::Rng;
use std::path::Path;

pub struct SpecModel {
    nz: Noise,
    pen: Pen,
    panels: Vec<Panel>,
    defs: Vec<SpDef>,
    h: Vec<f64>,
    pub build_secs: f64,
}

#[derive(Clone, Debug, Default)]
pub struct PanelPred {
    pub panel: String,
    pub channel_z: u32,
    pub mz: Vec<f64>,
    /// DIGITISED panel heights [% of base peak]
    pub data: Vec<f64>,
    /// MODEL, best calibration shift
    pub model: Vec<f64>,
    /// 16th and 84th percentile of the model under the transmission errors
    pub lo: Vec<f64>,
    pub hi: Vec<f64>,
    pub floor: f64,
    pub shift: f64,
    pub nll: f64,
    pub n_gauss: usize,
    /// species that contribute nothing because the live engine does not simulate them
    pub missing: Vec<String>,
    pub proxied: Vec<String>,
}

/// Model columns: max-plus combination of species responses, tallest column scaled to the panel's observed base plateau
/// (exactly the normalisation of chsim's `panel_nll`).
fn render(resp: &[Response], amp: &[f64], rel: &[f64], n: usize, base_obs: f64) -> Vec<f64> {
    let (h, _, _, _) = combine(resp, amp, rel, n);
    let hmax = h.iter().cloned().fold(0.0, f64::max);
    let scale = if hmax > 0.0 { base_obs / hmax } else { 1.0 };
    h.iter().map(|x| x * scale).collect()
}

impl SpecModel {
    pub fn build(dir: &Path) -> Result<SpecModel, String> {
        let t0 = std::time::Instant::now();
        let d = dir.to_string_lossy().to_string();
        let det = Orbitrap { r200: 70000.0 };
        let nz = Noise { sigq_px: 0.5, rho_p: 0.015, rho_t: 0.25 };
        let panels = load_fig5(&d)?;
        let refs: Vec<&Panel> = panels.iter().collect();
        let (pen, _) = fit_pen(&refs);
        let mut defs = fig5_species(&det);
        let pc = panels.iter().find(|p| p.name == "C").ok_or("panel C missing")?;
        let src = fit_5c(&mut defs, pc, &pen, &nz, &det, 0.5, 4);
        Ok(SpecModel { nz, pen, panels, defs, h: src.h.clone(), build_secs: t0.elapsed().as_secs_f64() })
    }

    pub fn panel_name(channel_z: u32) -> Option<&'static str> { match channel_z { 9 => Some("A1"), 8 => Some("A2"), 7 => Some("A3"), 6 => Some("A4"), _ => None } }

    /// `t(z)` = (per-pass T, relative error) of ubiquitin z+ from the live engine, or None if z is not simulated.
    pub fn predict(&self, channel_z: u32, t: &dyn Fn(u32) -> Option<(f64, f64)>) -> Option<PanelPred> {
        let pname = Self::panel_name(channel_z)?;
        let panel = self.panels.iter().find(|p| p.name == pname)?;
        let n = self.defs.len();
        let (mut amp, mut rel) = (vec![0.0; n], vec![0.0; n]);
        let (mut missing, mut proxied) = (vec![], vec![]);
        for (j, d) in self.defs.iter().enumerate() {
            let tv = match d.sim.kind {
                Kind::Ubq(z) => t(z),
                Kind::Adduct(z) => { let v = t(z); if v.is_some() { proxied.push(d.sim.name.clone()); } v }
                Kind::Chem(_) => None,
            };
            match tv { Some((tt, r)) => { amp[j] = self.h[j] * tt; rel[j] = r; } None => if self.h[j] > 0.0 { missing.push(d.sim.name.clone()); } }
        }
        if amp.iter().all(|&a| a <= 0.0) { return None; }
        let plac: Vec<_> = self.defs.iter().map(|d| d.plac.collapsed()).collect();
        let res = eval_panel(panel, &plac, &self.pen, &amp, &rel, true, &self.nz);
        // Monte-Carlo band: lognormal transmission errors, responses at the best shift
        let resp: Vec<Response> = plac.iter().map(|s| response(s, &panel.grid, &self.pen, res.delta_best)).collect();
        let nc = panel.grid.n;
        let mut draws: Vec<Vec<f64>> = vec![vec![]; nc];
        let mut rng = Rng::new(0x5EC7 ^ channel_z as u64);
        for _ in 0..48 {
            let a: Vec<f64> = (0..n).map(|j| if amp[j] > 0.0 { let r = rel[j].min(1.5); amp[j] * (r * rng.normal() - 0.5 * r * r).exp() } else { 0.0 }).collect();
            let m = render(&resp, &a, &rel, nc, panel.base_obs);
            for c in 0..nc { draws[c].push(m[c]); }
        }
        let (mut lo, mut hi) = (vec![0.0; nc], vec![0.0; nc]);
        for c in 0..nc {
            let v = &mut draws[c];
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            lo[c] = v[(v.len() as f64 * 0.16) as usize];
            hi[c] = v[((v.len() as f64 * 0.84) as usize).min(v.len() - 1)];
        }
        let det = &res.detail;
        Some(PanelPred {
            panel: pname.into(), channel_z, mz: det.mz.clone(), data: det.data.clone(), model: det.model.clone(), lo, hi, floor: panel.floor, shift: res.delta_best,
            nll: res.nll_best, n_gauss: res.split.n_gauss, missing, proxied,
        })
    }

}
