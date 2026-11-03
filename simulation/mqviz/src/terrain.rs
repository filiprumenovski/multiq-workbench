//! "Frequency terrain": local secular frequency f(z) of each species along the axis.
//! Method of chsim `xtalk-freq`: principal-axis RF and DC gradients on the axis -> Mathieu (a, q) -> exact beta
//! (continued fraction) -> f = beta * f_RF / 2. Gradients are precomputed per z in `World::terr`; this only recombines them.
use crate::world::*;

/// Mathieu characteristic exponent beta(a, q) (chsim xtalk.rs), valid inside the first stability region.
pub fn beta(a: f64, q: f64) -> f64 {
    let mut b = (a + q * q / 2.0).max(1e-12).sqrt();
    for _ in 0..60 {
        let cf = |b: f64, sgn: f64| {
            let mut s = 0.0;
            for k in (1..12).rev() {
                let kk = b + sgn * 2.0 * k as f64;
                s = q * q / (kk * kk - a - s);
            }
            s
        };
        let nb = (a + cf(b, 1.0) + cf(b, -1.0)).max(1e-12).sqrt();
        if (nb - b).abs() < 1e-12 { b = nb; break; }
        b = 0.5 * (b + nb);
    }
    b
}

/// Secular frequency [kHz] for principal RF gradient `gr`, DC gradient `gd` [V/m^2 per volt], charge-to-mass `qm`.
pub fn sec_freq_khz(gr: f64, gd: f64, qm: f64, vrf: f64) -> f64 {
    let q = 2.0 * qm * gr.abs() * vrf / (OMEGA * OMEGA);
    let a = -4.0 * qm * gd / (OMEGA * OMEGA);
    if a + q * q / 2.0 <= 0.0 || q > 0.95 { f64::NAN } else { beta(a, q) * F_RF / 1e3 / 2.0 }
}

#[derive(Clone)]
pub struct Terrain {
    pub z_mm: Vec<f64>,
    /// [species][z] f_x, f_y in kHz (NaN where the motion is unstable / untrapped).
    pub fx: Vec<Vec<f64>>,
    pub fy: Vec<Vec<f64>>,
}

pub fn compute(t: &TerrainTables, species: &[Sp], vrf: f64, dc: [f64; 4]) -> Terrain {
    let n = t.z_mm.len();
    let mut fx = vec![vec![f64::NAN; n]; species.len()];
    let mut fy = fx.clone();
    for (s, sp) in species.iter().enumerate() {
        let qm = sp.z * QE / (sp.mass * AMU);
        for k in 0..n {
            let (mut gdx, mut gdy) = (0.0, 0.0);
            for b in 0..4 { gdx += dc[b] * t.dc[b][k][0]; gdy += dc[b] * t.dc[b][k][1]; }
            fx[s][k] = sec_freq_khz(t.rf[k][0], gdx, qm, vrf);
            fy[s][k] = sec_freq_khz(t.rf[k][1], gdy, qm, vrf);
        }
    }
    Terrain { z_mm: t.z_mm.clone(), fx, fy }
}

impl Terrain {
    /// z positions [mm] where a curve crosses `f_khz` (linear interpolation between finite neighbours).
    pub fn crossings(z: &[f64], f: &[f64], f_khz: f64) -> Vec<f64> {
        let mut out = vec![];
        for k in 1..f.len() {
            let (a, b) = (f[k - 1] - f_khz, f[k] - f_khz);
            if a.is_finite() && b.is_finite() && (a <= 0.0) != (b <= 0.0) {
                let t = a / (a - b);
                out.push(z[k - 1] + t * (z[k] - z[k - 1]));
            }
        }
        out
    }
    /// Body-plateau frequency at z = 6 mm [kHz] (x mode).
    pub fn plateau(&self, s: usize) -> f64 {
        let k = self.z_mm.iter().position(|&z| z >= 6.0).unwrap_or(0);
        self.fx[s][k]
    }
    /// Index of the species whose plateau frequency is closest to `f_khz`.
    pub fn nearest(&self, f_khz: f64) -> usize {
        (0..self.fx.len()).filter(|&s| self.plateau(s).is_finite())
            .min_by(|&a, &b| (self.plateau(a) - f_khz).abs().partial_cmp(&(self.plateau(b) - f_khz).abs()).unwrap()).unwrap_or(0)
    }
}
