//! Species mixtures and instrument settings for the three presets.
use crate::world::Sp;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PresetId { S14, Fig5, Custom }

#[derive(Clone, Debug)]
pub struct Preset {
    pub id: PresetId,
    pub name: &'static str,
    pub species: Vec<Sp>,
    pub v_rf: f64,
    pub v_ac: f64,
    pub f_ac_khz: f64,
    pub gate_v: f64,
    pub pressure_pa: f64,
    pub sweep: (f64, f64, f64), // lo, hi, step [kHz]
}

/// N2 CCS of ubiquitin charge states [A^2] (chsim/src/reality.rs `ubq()`).
pub fn ubq_ccs(z: u32) -> f64 {
    match z { 6 => 1250.0, 7 => 1262.0, 8 => 1800.0, 9 => 2090.0, _ => 2200.0 }
}

/// Fig. 5C control (channel input) peak heights, % of base peak (data/their_side/figure_data/fig5_labelled_peaks.csv, panel C).
/// 6+ 15.9, 7+ 60.3, 8+ 100.0, 9+ 15.2, 10+ 4.8.
pub const FIG5C_HEIGHT: [(u32, f64); 5] = [(6, 15.9), (7, 60.3), (8, 100.0), (9, 15.2), (10, 4.8)];

pub fn ubq(z: u32, weight: f64) -> Sp {
    Sp { label: format!("{}+", z), mass: 8564.8 + z as f64 * 1.00728, z: z as f64, ccs: ubq_ccs(z), weight }
}

/// Empirical N2 CCS of a compact ion of mass `m` Da [A^2] (2.95 M^(2/3): 425 at 1727 Da, 1235 at 8565 Da).
pub fn ccs_guess(m: f64) -> f64 { 2.95 * m.powf(2.0 / 3.0) }

pub fn custom_species(mz: f64, z: f64, ccs: f64) -> Sp {
    let mass = mz * z;
    Sp { label: format!("{}+ m/z {:.1}", z as i32, mz), mass, z, ccs: if ccs > 0.0 { ccs } else { ccs_guess(mass) }, weight: 1.0 }
}

pub fn preset(id: PresetId) -> Preset {
    match id {
        PresetId::S14 => Preset {
            id, name: "S14 · Bak BH3 2+",
            species: vec![
                Sp { label: "2+ 863.5".into(), mass: 2.0 * 863.473, z: 2.0, ccs: 425.0, weight: 0.95 },
                Sp { label: "2+ Na 874.5".into(), mass: 2.0 * 874.466, z: 2.0, ccs: 430.0, weight: 0.05 },
            ],
            v_rf: 70.0, v_ac: 0.35, f_ac_khz: 84.25, gate_v: 24.0, pressure_pa: 0.1, sweep: (70.0, 100.0, 1.0),
        },
        PresetId::Fig5 => Preset {
            id, name: "Fig. 5 · ubiquitin 6+…10+",
            species: FIG5C_HEIGHT.iter().map(|&(z, h)| ubq(z, h)).collect(),
            v_rf: 56.0, v_ac: 0.20, f_ac_khz: 59.25, gate_v: 24.0, pressure_pa: 0.1, sweep: (35.0, 65.0, 1.0),
        },
        PresetId::Custom => Preset {
            id, name: "Custom m/z",
            species: vec![custom_species(863.5, 2.0, 0.0)],
            v_rf: 70.0, v_ac: 0.35, f_ac_khz: 84.0, gate_v: 24.0, pressure_pa: 0.1, sweep: (70.0, 100.0, 1.0),
        },
    }
}
