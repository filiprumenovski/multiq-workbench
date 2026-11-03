//! m/z-selective ejection channel ("tilted quadrupole") from Krutchinsky's Ion-Ejection-from-MultiQ code.
//! Geometry follows docs/channel_geometry_spec.md exactly (parametric, with their quirks behind switches).

use fieldsolve::{Cell, Grid, Problem};
use metal::*;
use std::cell::RefCell;
use std::collections::HashMap;
use std::ffi::c_void;
use std::time::Instant;

const AMU: f64 = 1.660_539_066_6e-27;
const QE: f64 = 1.602_176_634e-19;
const KB: f64 = 1.380_649e-23;
const POS_SCALE: f64 = 1.6e10;

// ---------------------------------------------------------------- geometry

#[derive(Clone, Debug)]
pub struct Geom {
    pub de: f64,          // rod diameter [mm]
    pub coef: f64,        // gap/De (0.323 => r_rod = 1.148 r0)
    pub q1_start: f64,    // [mm]
    pub q1_len: f64,
    pub alpha: f64,       // squeeze per mm of z (their 0.02)
    pub gap: f64,
    pub plate: f64,
    pub hole_d: f64,
    pub q2_len: f64,
    pub stem_d: f64,
    pub q2_flare_bug: bool, // their code reuses the Q1 squeeze formula in Q2 (negative sq)
    pub grow: f64,        // staircase correction: convex electrodes grown by grow*h (0 = their rasterization; 0.4 recommended)
    pub tilt_mode: u8,    // 0 = their code: all four Q1 rods shifted in y by ±sq; 1 = paper/inset: only the upper pair (LU,RU) converges in x; 2 = symmetric radial squeeze (x/y degenerate)
    pub ac_mode: u8,      // 0 = their code: left rods -1, right rods +1; 1 = paper: AC only across the non-inclined (lower) pair LL(-1), RL(+1)
    pub mis: [[f64; 4]; 4], // per Q1 rod [LL,RL,LU,RU]: [dx mm, dy mm, dx/dz rad, dy/dz rad] about the rod midpoint (tolerance studies)
    pub n_xy: usize,      // 160 nodes at 0.1 mm
    pub h_mm: f64,
}

impl Default for Geom {
    fn default() -> Self {
        Geom { de: 6.35, coef: 0.323, q1_start: 3.0, q1_len: 12.0, alpha: 0.02, gap: 1.0, plate: 1.6, hole_d: 7.0,
               q2_len: 3.0, stem_d: 4.6, q2_flare_bug: true, grow: 0.0, tilt_mode: 0, ac_mode: 0, mis: [[0.0; 4]; 4], n_xy: 160, h_mm: 0.1 }
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum Elec { Q1(u8), Stem3(u8), Plate, Stem5(u8), Q2(u8), Stem7(u8), Slab }

impl Geom {
    fn layers(&self) -> [f64; 8] {
        let l1 = self.q1_start; let l2 = l1 + self.q1_len; let l3 = l2 + self.gap; let l4 = l3 + self.plate;
        let l5 = l4 + self.gap; let l6 = l5 + self.q2_len; let l7 = l6 + self.gap; let l8 = l7 + self.plate;
        [l1, l2, l3, l4, l5, l6, l7, l8]
    }
    pub fn z_len(&self) -> f64 { self.layers()[7] }
    pub fn nz(&self) -> usize { (self.z_len() / self.h_mm + 1e-6) as usize }
    pub fn center(&self) -> f64 { self.coef * self.de / 2.0 + self.de / 2.0 }
    pub fn r0(&self) -> f64 { (self.de + self.coef * self.de) / 2f64.sqrt() - self.de / 2.0 }
    pub fn axis(&self) -> f64 { (self.n_xy / 2) as f64 * self.h_mm }

    /// Classify a node by its indices (exact decimal coordinates like their float grid).
    fn classify(&self, i: usize, j: usize, k: usize) -> Option<Elec> {
        let (x, y, z) = (i as f64 * self.h_mm, j as f64 * self.h_mm, k as f64 * self.h_mm);
        let a = self.axis();
        let quad = |x: f64, y: f64| -> Option<u8> {
            if x < a && y < a { Some(0) } else if x > a && y < a { Some(1) } else if x < a && y > a { Some(2) } else if x > a && y > a { Some(3) } else { None }
        };
        let c = self.center();
        let zmid = self.q1_start + self.q1_len / 2.0;
        let rod_q = |q: u8, sq: f64, d: f64, q1: bool| -> bool {
            let (mut cx, mut cy) = match (q, self.tilt_mode) {
                (0, 0) => (a - c, a - c + sq), (1, 0) => (a + c, a - c + sq), (2, 0) => (a - c, a + c - sq), (_, 0) => (a + c, a + c - sq),
                (0, 2) => (a - c + sq / 2f64.sqrt(), a - c + sq / 2f64.sqrt()), (1, 2) => (a + c - sq / 2f64.sqrt(), a - c + sq / 2f64.sqrt()),
                (2, 2) => (a - c + sq / 2f64.sqrt(), a + c - sq / 2f64.sqrt()), (_, 2) => (a + c - sq / 2f64.sqrt(), a + c - sq / 2f64.sqrt()),
                (0, _) => (a - c, a - c), (1, _) => (a + c, a - c), (2, _) => (a - c + sq, a + c), (_, _) => (a + c - sq, a + c),
            };
            if q1 { let m = self.mis[q as usize]; cx += m[0] + m[2] * (z - zmid); cy += m[1] + m[3] * (z - zmid); }
            ((x - cx).powi(2) + (y - cy).powi(2)).sqrt() < d / 2.0 + self.grow * self.h_mm
        };
        let rod = |q: u8, sq: f64, d: f64| rod_q(q, sq, d, false);
        let [l1, l2, l3, l4, l5, l6, l7, l8] = self.layers();
        let eps = 1e-9;
        if z <= l1 + eps { return None; }
        if z <= l2 + eps {
            let sq = self.alpha * (l2 - z);
            return quad(x, y).filter(|&q| rod_q(q, sq, self.de, true)).map(Elec::Q1);
        }
        if z <= l3 + eps { return quad(x, y).filter(|&q| rod(q, 0.0, self.stem_d)).map(Elec::Stem3); }
        if z <= l4 + eps {
            let r2 = (x - a).powi(2) + (y - a).powi(2);
            return if r2 < (self.hole_d / 2.0 - self.grow * self.h_mm).powi(2) { None } else { Some(Elec::Plate) };
        }
        if z <= l5 + eps { return quad(x, y).filter(|&q| rod(q, 0.0, self.stem_d)).map(Elec::Stem5); }
        if z <= l6 + eps {
            let sq = if self.q2_flare_bug { self.alpha * (l2 - z) } else { 0.0 };
            return quad(x, y).filter(|&q| rod(q, sq, self.de)).map(Elec::Q2);
        }
        if z <= l7 + eps { return quad(x, y).filter(|&q| rod(q, 0.0, self.stem_d)).map(Elec::Stem7); }
        if z <= l8 + eps { return Some(Elec::Slab); }
        None
    }
}

fn rf_sign(q: u8) -> f64 { match q { 0 | 3 => -1.0, _ => 1.0 } }   // LL,RU = -1; RL,LU = +1
fn ac_sign(q: u8) -> f64 { match q { 0 | 2 => -1.0, _ => 1.0 } }   // left -1, right +1

/// Basis identifiers.
#[derive(Clone, Copy, Debug)]
pub enum Basis { Rf, Ac, DcQ1, DcPlate, DcDown, DcBox }

fn basis_value(b: Basis, e: Elec, ac_mode: u8) -> f64 {
    match b {
        Basis::Rf => match e { Elec::Q1(q) | Elec::Stem3(q) | Elec::Stem5(q) | Elec::Q2(q) | Elec::Stem7(q) => rf_sign(q), _ => 0.0 },
        Basis::Ac => match e { Elec::Q1(q) | Elec::Stem3(q) => if ac_mode == 1 && q >= 2 { 0.0 } else { ac_sign(q) }, _ => 0.0 },
        Basis::DcQ1 => match e { Elec::Q1(_) | Elec::Stem3(_) => 1.0, _ => 0.0 },
        Basis::DcPlate => match e { Elec::Plate => 1.0, _ => 0.0 },
        Basis::DcDown => match e { Elec::Stem5(_) | Elec::Q2(_) | Elec::Stem7(_) | Elec::Slab => 1.0, _ => 0.0 },
        Basis::DcBox => 0.0,
    }
}

pub struct Fields {
    pub geom: Geom,
    pub nx: usize, pub ny: usize, pub nz: usize, pub h: f64,
    pub rf: Vec<[f32; 4]>,          // E per volt; w = 1 for conductor nodes (mask)
    pub ac: Vec<[f32; 4]>,
    pub dc: [Vec<[f32; 4]>; 4],     // DcQ1, DcPlate, DcDown, DcBox
    pub phi_rf: Vec<f64>, pub phi_ac: Vec<f64>, pub phi_dc: [Vec<f64>; 4],
}

fn cache_path(g: &Geom, b: Basis) -> std::path::PathBuf {
    let key = format!("{:?}{:?}", g, b);
    let mut h: u64 = 1469598103934665603;
    for c in key.bytes() { h ^= c as u64; h = h.wrapping_mul(1099511628211); }
    cache_dir().join(format!("ch_{:016x}.bin", h))
}

/// Field-cache directory: `$CHSIM_CACHE_DIR` if set, else `./cache` (original behaviour).
pub fn cache_dir() -> std::path::PathBuf {
    std::env::var_os("CHSIM_CACHE_DIR").map(std::path::PathBuf::from).unwrap_or_else(|| std::path::PathBuf::from("cache"))
}

fn solve_basis(g: &Geom, b: Basis, mask_out: &mut Vec<u8>) -> Vec<f64> {
    let (nx, ny, nz) = (g.n_xy, g.n_xy, g.nz());
    let path = cache_path(g, b);
    let grid = Grid { nx, ny, nz, h: g.h_mm * 1e-3, origin: [0.0; 3] };
    let p = Problem::from_fn(grid, |x, y, z| {
        let (i, j, k) = ((x / (g.h_mm * 1e-3)).round() as usize, (y / (g.h_mm * 1e-3)).round() as usize, (z / (g.h_mm * 1e-3)).round() as usize);
        match g.classify(i, j, k) { Some(e) => Cell::Fixed(basis_value(b, e, g.ac_mode)), None => Cell::Free }
    });
    if mask_out.is_empty() { *mask_out = p.mask.clone(); }
    if let Ok(bytes) = std::fs::read(&path) {
        if bytes.len() == nx * ny * nz * 8 {
            return bytes.chunks_exact(8).map(|c| f64::from_le_bytes(c.try_into().unwrap())).collect();
        }
    }
    let outer = if let Basis::DcBox = b { 1.0 } else { 0.0 };
    let t = Instant::now();
    let (sol, st) = fieldsolve::solve(&p, outer, 1e-7);
    eprintln!("  solved {:?}: {} iters, residual {:.1e}, {:.1} s", b, st.iterations, st.final_residual, t.elapsed().as_secs_f64());
    let _ = std::fs::create_dir_all(cache_dir());
    let bytes: Vec<u8> = sol.phi.iter().flat_map(|v| v.to_le_bytes()).collect();
    let _ = std::fs::write(&path, bytes);
    sol.phi
}

/// E = -grad(phi) by central differences (one-sided at the box border), per volt.
fn efield(phi: &[f64], nx: usize, ny: usize, nz: usize, h: f64) -> Vec<[f32; 4]> {
    let idx = |i: usize, j: usize, k: usize| i + nx * (j + ny * k);
    let mut out = vec![[0f32; 4]; phi.len()];
    for k in 0..nz { for j in 0..ny { for i in 0..nx {
        let d = |a: usize, b: usize, da: f64| -(phi[b] - phi[a]) / da;
        let (im, ip) = (i.saturating_sub(1), (i + 1).min(nx - 1));
        let (jm, jp) = (j.saturating_sub(1), (j + 1).min(ny - 1));
        let (km, kp) = (k.saturating_sub(1), (k + 1).min(nz - 1));
        out[idx(i, j, k)] = [
            d(idx(im, j, k), idx(ip, j, k), (ip - im) as f64 * h) as f32,
            d(idx(i, jm, k), idx(i, jp, k), (jp - jm) as f64 * h) as f32,
            d(idx(i, j, km), idx(i, j, kp), (kp - km) as f64 * h) as f32,
            0.0,
        ];
    }}}
    out
}

impl Fields {
    /// Fields from THEIR exported potential arrays (float32 RF/DC/AC + masks; see legacy_java/results/potentials).
    /// Their DC basis (plate +1, downstream -1) goes into the plate slot: drive dc = [0, V_DC, 0, 0].
    pub fn from_theirs(dir: &str) -> Fields {
        let g = Geom::default();
        let (nx, ny, nz, h) = (160usize, 160usize, 242usize, 1e-4);
        let rd = |f: &str| -> Vec<f64> { std::fs::read(format!("{}/{}", dir, f)).unwrap().chunks_exact(4).map(|c| f32::from_le_bytes(c.try_into().unwrap()) as f64).collect() };
        let rm = |f: &str| -> Vec<u8> { std::fs::read(format!("{}/{}", dir, f)).unwrap() };
        let (phi_rf, phi_dc_t, phi_ac) = (rd("eject_RF.f32"), rd("eject_DC.f32"), rd("eject_AC.f32"));
        let (m_rf, m_ac) = (rm("eject_RF_electrode.u8"), rm("eject_AC_electrode.u8"));
        let mut rf = efield(&phi_rf, nx, ny, nz, h);
        for n in 0..rf.len() { rf[n][3] = (m_rf[n] | m_ac[n]) as f32; }
        let zero = vec![0.0f64; nx * ny * nz];
        let zf = vec![[0f32; 4]; nx * ny * nz];
        Fields { geom: g, nx, ny, nz, h, ac: efield(&phi_ac, nx, ny, nz, h), dc: [zf.clone(), efield(&phi_dc_t, nx, ny, nz, h), zf.clone(), zf],
                 rf, phi_rf, phi_ac, phi_dc: [zero.clone(), phi_dc_t, zero.clone(), zero] }
    }

    pub fn build(g: &Geom) -> Fields {
        let (nx, ny, nz, h) = (g.n_xy, g.n_xy, g.nz(), g.h_mm * 1e-3);
        eprintln!("channel grid {}x{}x{} @ {} mm", nx, ny, nz, g.h_mm);
        let mut mask = Vec::new();
        let phi_rf = solve_basis(g, Basis::Rf, &mut mask);
        let phi_ac = solve_basis(g, Basis::Ac, &mut mask);
        let phi_dc = [solve_basis(g, Basis::DcQ1, &mut mask), solve_basis(g, Basis::DcPlate, &mut mask),
                      solve_basis(g, Basis::DcDown, &mut mask), solve_basis(g, Basis::DcBox, &mut mask)];
        let mut rf = efield(&phi_rf, nx, ny, nz, h);
        for (e, m) in rf.iter_mut().zip(mask.iter()) { e[3] = *m as f32; }
        let ac = efield(&phi_ac, nx, ny, nz, h);
        let dc = [efield(&phi_dc[0], nx, ny, nz, h), efield(&phi_dc[1], nx, ny, nz, h), efield(&phi_dc[2], nx, ny, nz, h), efield(&phi_dc[3], nx, ny, nz, h)];
        Fields { geom: g.clone(), nx, ny, nz, h, rf, ac, dc, phi_rf, phi_ac, phi_dc }
    }
    fn idx(&self, i: usize, j: usize, k: usize) -> usize { i + self.nx * (j + self.ny * k) }
    /// Combined static DC field for electrode voltages [q1, plate, down, box].
    pub fn dc_field(&self, v: [f64; 4]) -> Vec<[f32; 4]> {
        (0..self.rf.len()).map(|n| {
            let mut e = [0f32; 4];
            for b in 0..4 { for d in 0..3 { e[d] += (v[b] as f32) * self.dc[b][n][d]; } }
            e
        }).collect()
    }
    pub fn dc_phi_axis(&self, v: [f64; 4]) -> Vec<(f64, f64)> {
        let a = self.nx / 2;
        (0..self.nz).map(|k| { let n = self.idx(a, a, k); (k as f64 * self.geom.h_mm, (0..4).map(|b| v[b] * self.phi_dc[b][n]).sum()) }).collect()
    }
}

// ---------------------------------------------------------------- species, drive, physics

#[derive(Clone, Copy, Debug)]
pub struct Species { pub mass: f64, pub z: f64, pub ccs_n2: f64 }
impl Species {
    pub fn mz(&self) -> f64 { self.mass / self.z }
    fn ccs_he_poly(&self) -> f64 { let m = self.mass; 2.81e-9 * m * m * m - 3.55e-5 * m * m + 0.232 * m + 41.91 }
}

#[derive(Clone, Copy, Debug)]
pub struct Drive { pub v_rf: f64, pub f_rf: f64, pub v_ac: f64, pub f_ac: f64, pub dc: [f64; 4] }

/// RF + AC drive of one job (DC is shared by a whole batch).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RfAc { pub v_rf: f64, pub f_rf: f64, pub v_ac: f64, pub f_ac: f64 }
impl Drive { pub fn rfac(&self) -> RfAc { RfAc { v_rf: self.v_rf, f_rf: self.f_rf, v_ac: self.v_ac, f_ac: self.f_ac } } }

/// One scan point of a batch: species, RF/AC drive, ion count, injection seed.
#[derive(Clone, Copy, Debug)]
pub struct Job { pub sp: Species, pub drive: RfAc, pub n: usize, pub seed: u64 }

#[derive(Clone, Copy, Debug)]
pub struct Physics {
    pub legacy_force: bool, pub legacy_collisions: bool, pub legacy_fields: bool, pub legacy_injection: bool,
    pub legacy_criteria: bool, pub pressure_pa: f64, pub temp_k: f64, pub dt: f64, pub t_max: f64,
    pub r_inject_mm: f64,
}
impl Physics {
    pub fn legacy(pressure_pa: f64) -> Physics {
        Physics { legacy_force: true, legacy_collisions: true, legacy_fields: true, legacy_injection: true, legacy_criteria: true,
                  pressure_pa, temp_k: 300.0, dt: 1e-7, t_max: 1.0, r_inject_mm: 0.0 }
    }
    /// Correct-physics defaults. dt = T_RF/50 (4e-8 s at 500 kHz): converged vs T/100 within statistical error on the R28
    /// workloads (chsim/JUMP.md, `chsim dtconv`); `CHSIM_DT_DIV=100` selects the original reference step T/100 (2e-8 s).
    /// Validated for q_Mathieu <~ 0.45; use T/100 for larger q.
    pub fn correct(pressure_pa: f64) -> Physics {
        let div: f64 = std::env::var("CHSIM_DT_DIV").ok().and_then(|v| v.parse().ok()).unwrap_or(50.0);  // T/50 default: peak bias +0.028 kHz (+0.04%) = leapfrog dt^2 law (chsim peakbias, JUMP.md); CHSIM_DT_DIV=100 = reference
        Physics { legacy_force: false, legacy_collisions: false, legacy_fields: false, legacy_injection: false, legacy_criteria: false,
                  pressure_pa, temp_k: 300.0, dt: 2e-6 / div, t_max: 0.02, r_inject_mm: 2.0 }
    }
}

/// Batch-wide kernel parameters (layout mirrors `ChParams` in channel.metal).
#[repr(C)]
#[derive(Clone, Copy, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub(crate) struct ChParams {
    pub(crate) dt: f32, pub(crate) v_rf: f32, pub(crate) w_rf: f32, inv_h: f32, h: f32, pos_scale: f32,
    z_return: f32, z_exit: f32, t_max: f32, crash_z0: f32, crash_z1: f32,
    pub(crate) nx: i32, pub(crate) ny: i32, pub(crate) nz: i32, pub(crate) nsteps: i32, pub(crate) nearest: i32, pub(crate) legacy: i32,
}

/// Per-job kernel parameters (layout mirrors `ChJob`).
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, bytemuck::Pod, bytemuck::Zeroable)]
pub(crate) struct ChJob {
    qm_kick: f32, v_ac: f32, w_ac: f32, nsig_dt: f32, s_gas: f32, v_bound: f32, mr: f32, legacy_floor: f32,
    legacy_mfp_mean: f32, pad: [f32; 3],
}

#[repr(C)]
#[derive(Clone, Copy, Debug, Default, bytemuck::Pod, bytemuck::Zeroable)]
pub struct ChIon { pub(crate) p: [i32; 4], pub(crate) v: [f32; 4], pub(crate) aux: [f32; 4], pub(crate) rng: [u32; 4] }

#[derive(Clone, Copy, Debug, Default)]
pub struct Outcome { pub n: usize, pub crashed: usize, pub returned: usize, pub transmitted: usize, pub left: usize, pub timeout: usize, pub mean_coll: f64, pub mean_time: f64 }
impl Outcome {
    pub fn frac(&self, k: usize) -> f64 { k as f64 / self.n as f64 }
}

/// Execution path of the channel kernel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Path { Tiled, Direct }

/// Process-wide default path (`--untiled` / CHSIM_UNTILED=1 selects the direct per-ion kernel).
pub static DEFAULT_DIRECT: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

pub struct Engine {
    pub(crate) device: Device, pub(crate) queue: CommandQueue, pub(crate) pso: ComputePipelineState,
    pub(crate) rf: Buffer, pub(crate) ac: Buffer, pub fields: Fields,
    pub path: Path,
    dc_cache: RefCell<HashMap<[u64; 4], Buffer>>,
    pub(crate) bricks: RefCell<Option<([u64; 5], Buffer)>>,
    pub tile: RefCell<crate::tiled::TileCfg>,
    pub(crate) tpipes: RefCell<[Option<std::rc::Rc<crate::tiled::TPipes>>; 2]>,
    pub(crate) verbose: bool,
}

pub(crate) fn shared_buf(dev: &Device, bytes: usize) -> Buffer { dev.new_buffer(bytes.max(16) as u64, MTLResourceOptions::StorageModeShared) }
pub(crate) fn buf_from<T: bytemuck::Pod>(dev: &Device, d: &[T]) -> Buffer {
    if d.is_empty() { return shared_buf(dev, 16); }
    dev.new_buffer_with_data(d.as_ptr() as *const c_void, std::mem::size_of_val(d) as u64, MTLResourceOptions::StorageModeShared)
}

impl Engine {
    pub fn new(fields: Fields) -> Engine {
        let device = Device::system_default().expect("no Metal device");
        let opts = CompileOptions::new();
        opts.set_fast_math_enabled(false); // accuracy-first for validation
        let lib = device.new_library_with_source(include_str!("channel.metal"), &opts).expect("MSL compile failed");
        let pso = device.new_compute_pipeline_state_with_function(&lib.get_function("channel", None).unwrap()).unwrap();
        let queue = device.new_command_queue();
        if std::env::var("CHSIM_VERBOSE").is_ok() { eprintln!("  direct kernel max threads/TG: {}", pso.max_total_threads_per_threadgroup()); }
        let mk = |d: &[[f32; 4]]| device.new_buffer_with_data(d.as_ptr() as *const c_void, (d.len() * 16) as u64, MTLResourceOptions::StorageModeShared);
        let rf = mk(&fields.rf); let ac = mk(&fields.ac);
        let direct = DEFAULT_DIRECT.load(std::sync::atomic::Ordering::Relaxed) || std::env::var("CHSIM_UNTILED").map_or(false, |v| v == "1");
        Engine { device, queue, pso, rf, ac, fields, path: if direct { Path::Direct } else { Path::Tiled },
                 dc_cache: RefCell::new(HashMap::new()), bricks: RefCell::new(None), tile: RefCell::new(crate::tiled::TileCfg::from_env()), tpipes: RefCell::new([None, None]), verbose: std::env::var("CHSIM_VERBOSE").is_ok() }
    }

    fn inject(&self, sp: &Species, ph: &Physics, n: usize, seed: u64) -> Vec<ChIon> {
        let mut st = seed.wrapping_mul(0x9E3779B97F4A7C15) ^ 0xD1B54A32D192ED03;
        let mut next = move || { st ^= st << 13; st ^= st >> 7; st ^= st << 17; st };
        let m = sp.mass * AMU; let g = &self.fields.geom; let a = g.axis() * 1e-3;
        let legacy_mean_mfp = self.legacy_mfp_mean(sp, ph);
        (0..n).map(|_| {
            let mut u = || (next() >> 11) as f64 / (1u64 << 53) as f64;
            let (pos, vel);
            if ph.legacy_injection {
                let kin = 0.1 * (0.75 + 0.75 * u()) * sp.z;                 // eV (their VIN*random(0.75,1.5), q = z e)
                let vin = (2.0 * kin * QE / m).sqrt();
                let ang = (-10.0 + 20.0 * u()).to_radians();
                let vx = 0.707 * vin * ang.sin();
                pos = [a, a, 1.0e-3]; vel = [vx, vx, vin * ang.cos()];
            } else {
                let s = (KB * ph.temp_k / m).sqrt();
                let r = ph.r_inject_mm * 1e-3 * u().sqrt(); let th = 2.0 * std::f64::consts::PI * u();
                let gs = |u1: f64, u2: f64| (-2.0 * u1.max(1e-300).ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
                let vz = s * (-2.0 * u().max(1e-300).ln()).sqrt();          // flux-weighted (effusive) normal component
                pos = [a + r * th.cos(), a + r * th.sin(), 1.0e-3]; vel = [s * gs(u(), u()), s * gs(u(), u()), vz];
            }
            let phase = 2.0 * std::f64::consts::PI * u();
            let vmag = (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]).sqrt();
            let first_mfp = -u().max(1e-7).ln() * legacy_mean_mfp / 10.0;  // their first free path has an extra /10
            let steps0 = (first_mfp / vmag / ph.dt).floor() + 1.0;
            ChIon {
                p: [(pos[0] * POS_SCALE) as i32, (pos[1] * POS_SCALE) as i32, (pos[2] * POS_SCALE) as i32, 0],
                v: [vel[0] as f32, vel[1] as f32, vel[2] as f32, 0.0],
                aux: [(phase / (2.0 * std::f64::consts::PI)) as f32, 0.0, steps0 as f32, 0.0],
                rng: [next() as u32 | 1, next() as u32, next() as u32, next() as u32],
            }
        }).collect()
    }

    fn legacy_mfp_mean(&self, sp: &Species, ph: &Physics) -> f64 {
        let p_torr = ph.pressure_pa / 133.322 + 1e-10;
        let np = 2.68e5 * (p_torr / 760.0);
        let ccs = sp.ccs_he_poly() * 10.0 * sp.z.sqrt();
        1.0 / (np * ccs)
    }

    /// Combined static DC field for the electrode voltages, uploaded once per DC configuration.
    pub(crate) fn dc_buffer(&self, v: [f64; 4]) -> Buffer {
        let key = [v[0].to_bits(), v[1].to_bits(), v[2].to_bits(), v[3].to_bits()];
        if let Some(b) = self.dc_cache.borrow().get(&key) { return b.clone(); }
        let t0 = Instant::now();
        let dc = self.fields.dc_field(v);
        let b = buf_from(&self.device, &dc);
        if self.verbose { eprintln!("  combined DC field built+uploaded in {:.3} s", t0.elapsed().as_secs_f64()); }
        let mut cache = self.dc_cache.borrow_mut();
        if cache.len() >= 2 { cache.clear(); }       // keep at most two DC configurations resident (99 MB each)
        cache.insert(key, b.clone());
        b
    }

    pub(crate) fn batch_params(&self, ph: &Physics, rf: &RfAc) -> ChParams {
        let g = &self.fields.geom;
        ChParams {
            dt: ph.dt as f32, v_rf: rf.v_rf as f32, w_rf: (rf.f_rf * ph.dt) as f32,
            inv_h: (1.0 / self.fields.h) as f32, h: self.fields.h as f32, pos_scale: POS_SCALE as f32,
            z_return: if ph.legacy_criteria { 1.0e-3 } else { 0.5e-3 },
            z_exit: if ph.legacy_criteria { 23.7e-3 } else { ((g.layers()[4]) * 1e-3) as f32 },
            t_max: ph.t_max as f32,
            crash_z0: if ph.legacy_criteria { (g.layers()[0] * 1e-3) as f32 } else { -1.0 },
            crash_z1: if ph.legacy_criteria { (g.layers()[1] * 1e-3) as f32 } else { 1.0 },
            nx: self.fields.nx as i32, ny: self.fields.ny as i32, nz: self.fields.nz as i32,
            nsteps: ((2e-4 / ph.dt) as i32).max(1000), nearest: ph.legacy_fields as i32, legacy: ph.legacy_collisions as i32,
        }
    }

    fn job_params(&self, j: &Job, ph: &Physics) -> ChJob {
        let sp = &j.sp;
        let m = sp.mass * AMU; let mg = 28.0 * AMU;
        let qm = sp.z * QE / m;
        let n_gas = ph.pressure_pa / (KB * ph.temp_k);
        let s_gas = (KB * ph.temp_k / mg).sqrt();
        ChJob {
            qm_kick: (if ph.legacy_force { 0.5 } else { 1.0 } * qm) as f32,
            v_ac: j.drive.v_ac as f32, w_ac: (j.drive.f_ac * ph.dt) as f32,
            nsig_dt: (n_gas * sp.ccs_n2 * 1e-20 * ph.dt) as f32, s_gas: s_gas as f32, v_bound: (6.0 * s_gas) as f32,
            mr: (mg / (m + mg)) as f32,
            legacy_floor: if ph.legacy_collisions { (2.0 * 0.035 * 1.6e-19 / (sp.mass * 1.66e-27)).sqrt() as f32 } else { -1.0 },
            legacy_mfp_mean: self.legacy_mfp_mean(sp, ph) as f32, pad: [0.0; 3],
        }
    }

    /// Run one scan point (a batch of one job).
    #[allow(dead_code)]
    pub fn run(&self, sp: &Species, dr: &Drive, ph: &Physics, n: usize, seed: u64) -> Outcome {
        self.run_batch(&[Job { sp: *sp, drive: dr.rfac(), n, seed }], dr.dc, ph)[0]
    }

    /// Run many scan points together in one launch (all ions in one buffer, per-ion parameters from a job table).
    /// Every job gets exactly the ions/RNG streams it would get from `run` alone, so results are bitwise identical
    /// per ion. RF (v_rf, f_rf), DC and `ph` are shared by the batch. Returns one `Outcome` per job.
    pub fn run_batch(&self, jobs: &[Job], dc: [f64; 4], ph: &Physics) -> Vec<Outcome> {
        let (ions, jid) = self.run_batch_ions(jobs, dc, ph, self.path);
        summarize(&ions, &jid, jobs.len())
    }

    /// As `run_batch` but returns final per-ion states (concatenated in job order) and the job id of each ion.
    pub(crate) fn prepare(&self, jobs: &[Job], ph: &Physics) -> (Vec<ChIon>, Vec<u32>, Vec<ChJob>, ChParams) {
        assert!(!jobs.is_empty());
        let rf0 = jobs[0].drive;
        assert!(jobs.iter().all(|j| j.drive.v_rf == rf0.v_rf && j.drive.f_rf == rf0.f_rf), "RF drive must be shared by a batch");
        let pr = self.batch_params(ph, &rf0);
        let jp: Vec<ChJob> = jobs.iter().map(|j| self.job_params(j, ph)).collect();
        let mut ions: Vec<ChIon> = Vec::new(); let mut jid: Vec<u32> = Vec::new();
        for (k, j) in jobs.iter().enumerate() {
            let v = self.inject(&j.sp, ph, j.n, j.seed);
            jid.extend(std::iter::repeat(k as u32).take(v.len()));
            ions.extend(v);
        }
        (ions, jid, jp, pr)
    }

    pub fn run_batch_ions(&self, jobs: &[Job], dc: [f64; 4], ph: &Physics, path: Path) -> (Vec<ChIon>, Vec<u32>) {
        let (ions, jid, jp, pr) = self.prepare(jobs, ph);
        let dcb = self.dc_buffer(dc);
        let t_run = Instant::now();
        let ib = buf_from(&self.device, &ions); let jb = buf_from(&self.device, &jid); let jpb = buf_from(&self.device, &jp);
        match path {
            Path::Direct => self.run_direct(&ib, &jb, &jpb, &dcb, &pr, ions.len()),
            Path::Tiled => { self.run_tiled(&ib, &jb, &jpb, &dcb, dc, &pr, ions.len()); }
        }
        let out = unsafe { std::slice::from_raw_parts(ib.contents() as *const ChIon, ions.len()) };
        if self.verbose {
            let steps: f64 = out.iter().map(|i| i.v[3] as f64 / ph.dt).sum();
            let tmax = out.iter().map(|i| i.v[3]).fold(0f32, f32::max);
            eprintln!("  batch: {} ions, {:.3e} ion-steps, longest ion {:.1} ms ({:.2e} steps), wall {:.3} s => {:.3} G ion-steps/s",
                out.len(), steps, tmax * 1e3, tmax as f64 / ph.dt, t_run.elapsed().as_secs_f64(), steps / t_run.elapsed().as_secs_f64() / 1e9);
        }
        (out.to_vec(), jid)
    }

    /// The original per-ion kernel: one thread per ion, fields read from global memory.
    pub(crate) fn run_direct(&self, ib: &Buffer, jb: &Buffer, jpb: &Buffer, dcb: &Buffer, pr0: &ChParams, n: usize) {
        let pr = *pr0;
        loop {
            let cmd = self.queue.new_command_buffer();
            let enc = cmd.new_compute_command_encoder();
            enc.set_compute_pipeline_state(&self.pso);
            enc.set_buffer(0, Some(ib), 0); enc.set_buffer(1, Some(&self.rf), 0);
            enc.set_buffer(2, Some(dcb), 0); enc.set_buffer(3, Some(&self.ac), 0);
            enc.set_bytes(4, std::mem::size_of::<ChParams>() as u64, &pr as *const ChParams as *const c_void);
            enc.set_buffer(5, Some(jpb), 0); enc.set_buffer(6, Some(jb), 0);
            enc.dispatch_threads(MTLSize::new(n as u64, 1, 1), MTLSize::new(128, 1, 1));
            enc.end_encoding(); cmd.commit(); cmd.wait_until_completed();
            let out = unsafe { std::slice::from_raw_parts(ib.contents() as *const ChIon, n) };
            if out.iter().all(|i| i.p[3] != 0) { break; }
        }
    }
}

pub fn summarize(out: &[ChIon], jid: &[u32], njobs: usize) -> Vec<Outcome> {
    let mut os = vec![Outcome::default(); njobs];
    for (i, &j) in out.iter().zip(jid) {
        let o = &mut os[j as usize];
        o.n += 1;
        match i.p[3] { 1 => o.crashed += 1, 2 => o.returned += 1, 3 => o.transmitted += 1, 4 => o.left += 1, _ => o.timeout += 1 }
        o.mean_coll += i.aux[1] as f64; o.mean_time += i.v[3] as f64;
    }
    for o in os.iter_mut() { if o.n > 0 { o.mean_coll /= o.n as f64; o.mean_time /= o.n as f64; } }
    os
}

// ---------------------------------------------------------------- validation commands

/// Fit the near-axis RF potential phi = -2 X Y / r0^2 (+ higher order) per z-slice.
/// FWHM of a peaked curve (linear interpolation), returns (peak x, peak y, fwhm).
pub fn fwhm(xs: &[f64], ys: &[f64]) -> (f64, f64, f64) {
    let (mut im, mut ym) = (0, f64::MIN);
    for (i, &y) in ys.iter().enumerate() { if y > ym { ym = y; im = i; } }
    let half = ym / 2.0;
    let mut lo = xs[0]; let mut hi = xs[xs.len() - 1];
    for i in (1..=im).rev() { if ys[i - 1] < half { lo = xs[i - 1] + (half - ys[i - 1]) / (ys[i] - ys[i - 1]) * (xs[i] - xs[i - 1]); break; } }
    for i in im..xs.len() - 1 { if ys[i + 1] < half { hi = xs[i] + (ys[i] - half) / (ys[i] - ys[i + 1]) * (xs[i + 1] - xs[i]); break; } }
    (xs[im], ym, hi - lo)
}

fn r0_profile(f: &Fields) {
    let g = &f.geom; let a = f.nx / 2;
    println!("RF field check (ideal r0 = {:.4} mm):", g.r0());
    println!("   z[mm]   r0_eff[mm]   E_ac_x(axis)[V/m per V]");
    for &zmm in &[3.5, 5.0, 7.0, 9.0, 11.0, 13.0, 14.5] {
        let k = (zmm / g.h_mm).round() as usize;
        // least squares on points |X|,|Y| <= 0.5 mm: phi ≈ c XY
        let (mut sxy, mut sxx) = (0.0, 0.0);
        for dj in -5i64..=5 { for di in -5i64..=5 {
            if di == 0 || dj == 0 { continue; }
            let (x, y) = (di as f64 * g.h_mm, dj as f64 * g.h_mm);
            let phi = f.phi_rf[f.idx((a as i64 + di) as usize, (a as i64 + dj) as usize, k)];
            sxy += phi * x * y; sxx += (x * y) * (x * y);
        }}
        let c = sxy / sxx; // = -2/r0^2
        let r0 = (-2.0 / c).sqrt();
        let eac = f.ac[f.idx(a, a, k)][0];
        println!("   {:5.1}   {:8.4}     {:8.1}", zmm, r0, eac);
    }
}

/// f64 CPU trajectory in RF-only field (no gas, no DC, no AC); returns secular frequency from zero crossings.
fn secular_freq(f: &Fields, sp: &Species, v_rf: f64, f_rf: f64, z_mm: f64, half_force: bool) -> f64 {
    let g = &f.geom; let a = g.axis() * 1e-3; let h = f.h;
    let qm = sp.z * QE / (sp.mass * AMU) * if half_force { 0.5 } else { 1.0 };
    let w = 2.0 * std::f64::consts::PI * f_rf;
    let tri = |x: [f64; 3]| -> [f64; 3] {
        let gx = [x[0] / h, x[1] / h, x[2] / h];
        let (i, j, k) = (gx[0] as usize, gx[1] as usize, gx[2] as usize);
        let (tx, ty, tz) = (gx[0] - i as f64, gx[1] - j as f64, gx[2] - k as f64);
        let mut e = [0.0; 3];
        for (di, wx) in [(0, 1.0 - tx), (1, tx)] { for (dj, wy) in [(0, 1.0 - ty), (1, ty)] { for (dk, wz) in [(0, 1.0 - tz), (1, tz)] {
            let v = f.rf[f.idx(i + di, j + dj, k + dk)];
            for d in 0..3 { e[d] += wx * wy * wz * v[d] as f64; }
        }}}
        e
    };
    // small displacement along the u = (1,-1)/sqrt2 principal direction, no initial velocity
    let amp = 0.2e-3;
    let mut x = [a + amp / 2f64.sqrt(), a - amp / 2f64.sqrt(), z_mm * 1e-3];
    let mut v = [0.0; 3];
    let dt = 1.0 / f_rf / 200.0; let n = (4e-3 / dt) as usize;
    let mut crossings = Vec::new(); let mut prev = (x[0] - a) - (x[1] - a);
    let mut t = 0.0;
    for _ in 0..n {
        let e = tri(x); let s = v_rf * (w * t).sin();
        for d in 0..2 { v[d] += qm * s * e[d] * dt; }      // transverse only (z frozen at slice)
        for d in 0..2 { x[d] += v[d] * dt; }
        t += dt;
        let u = (x[0] - a) - (x[1] - a);
        if prev < 0.0 && u >= 0.0 { crossings.push(t); }
        prev = u;
        if (x[0] - a).abs() > 2.5e-3 || (x[1] - a).abs() > 2.5e-3 { return f64::NAN; }
    }
    if crossings.len() < 3 { return f64::NAN; }
    (crossings.len() - 1) as f64 / (crossings[crossings.len() - 1] - crossings[0])
}

/// RF-only solve at arbitrary resolution/tilt; prints r0_eff(z) from a near-axis fit. Light on memory.
fn field_check(h_mm: f64, alpha: f64) {
    let g = Geom { h_mm, n_xy: (16.0 / h_mm).round() as usize, alpha, ..Geom::default() };
    let (nx, nz) = (g.n_xy, g.nz());
    let mut mask = Vec::new();
    let t = Instant::now();
    let phi = solve_basis(&g, Basis::Rf, &mut mask);
    let a = nx / 2;
    print!("h={:.3} mm alpha={:.3}: r0_eff(z) =", h_mm, alpha);
    for &zmm in &[3.5, 5.0, 7.0, 9.0, 11.0, 13.0, 14.5] {
        let k = (zmm / h_mm).round() as usize;
        let w = (0.5 / h_mm).round() as i64;
        let (mut sxy, mut sxx) = (0.0, 0.0);
        for dj in -w..=w { for di in -w..=w {
            if di == 0 || dj == 0 { continue; }
            let (x, y) = (di as f64 * h_mm, dj as f64 * h_mm);
            let p = phi[(a as i64 + di) as usize + nx * ((a as i64 + dj) as usize + nx * k)];
            sxy += p * x * y; sxx += (x * y) * (x * y);
        }}
        print!("  z{:.1}:{:.4}", zmm, (-2.0 * sxx / sxy).sqrt());
    }
    println!("   ({:.1} s, {}x{}x{})", t.elapsed().as_secs_f64(), nx, nx, nz);
}

/// S14 geometry variants (see docs/01_results_log.md R16, R18). args[5] = variant, args[6] = gate hole override [mm].
pub fn s14_geom(args: &[String]) -> Geom {
    let variant = args.get(5).map(|s| s.as_str()).unwrap_or("code");
    let mut geom = match variant {
        "inset" => Geom { grow: 0.4, q2_flare_bug: false, tilt_mode: 1, ac_mode: 1, alpha: 0.525 / 12.0, ..Geom::default() },
        "paper" => Geom { grow: 0.4, q2_flare_bug: false, tilt_mode: 1, ac_mode: 1, alpha: 0.01, ..Geom::default() },
        // as-built from photo forensics (data/their_side/deep/figures_forensics.md): D 6.34 on 8.30 pitch, 0.021 rad/rod (4 rods)
        // or 0.042 rad on one pair, gate aperture ~5.75 mm, AC on two rods
        "asbuilt4" => Geom { grow: 0.4, q2_flare_bug: false, de: 6.34, coef: (8.30 - 6.34) / 6.34, alpha: 0.021, hole_d: 5.75, ac_mode: 1, ..Geom::default() },
        "asbuilt2" => Geom { grow: 0.4, q2_flare_bug: false, de: 6.34, coef: (8.30 - 6.34) / 6.34, alpha: 0.042, tilt_mode: 1, hole_d: 5.75, ac_mode: 1, ..Geom::default() },
        _ => Geom { grow: 0.4, q2_flare_bug: false, ..Geom::default() },
    };
    if let Some(hd) = args.get(6).and_then(|s| s.parse::<f64>().ok()) { geom.hole_d = hd; }
    println!("geometry variant: {} (tilt_mode {}, ac_mode {}, alpha {:.4} per rod, gate hole {:.2} mm)", variant, geom.tilt_mode, geom.ac_mode, geom.alpha, geom.hole_d);
    geom
}

/// S14 corrected geometry: staircase-compensated rods, Q2 parallel. Tilt as in their code (as-built geometry unknown).
pub(crate) fn s14_engine() -> Engine {
    let geom = Geom { grow: 0.4, q2_flare_bug: false, ..Geom::default() };
    Engine::new(Fields::build(&geom))
}

pub(crate) fn bak2() -> Species { Species { mass: 2.0 * 862.974, z: 2.0, ccs_n2: 425.0 } }

/// DC configs of fig. S14: (closed, open). With a barrier the gate voltage is chosen so the on-axis saddle sits `barrier` V above the trap.
pub(crate) fn s14_dc(eng: &Engine, barrier: Option<f64>, verbose: bool) -> ([f64; 4], [f64; 4]) {
    let mut closed = [18.0, 24.0, 15.0, 20.0]; // [channel quads, gate, downstream, trap/box]  (fig. S14)
    if let Some(b) = barrier {
        let saddle = |vg: f64| eng.fields.dc_phi_axis([18.0, vg, 15.0, 20.0]).iter().filter(|(z, _)| *z > 14.0 && *z < 19.0).map(|&(_, p)| p).fold(f64::MIN, f64::max);
        let (mut lo, mut hi) = (20.0, 2000.0);
        for _ in 0..60 { let mid = 0.5 * (lo + hi); if saddle(mid) < 20.0 + b { lo = mid } else { hi = mid } }
        closed[1] = hi;
        if verbose { println!("barrier-matched gate: {:.1} V gives on-axis saddle {:.2} V ({} V above trap)", hi, saddle(hi), b); }
    }
    (closed, [18.0, 19.0, 15.0, 20.0])
}

fn s14(args: &[String]) {
    let eng = Engine::new(Fields::build(&s14_geom(args)));
    let pa: f64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(0.1);
    let n: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(2000);
    let mut ph = Physics::correct(pa); ph.t_max = 5e-3;
    let bak2 = bak2();
    let barrier: Option<f64> = args.get(3).and_then(|s| s.parse().ok()).filter(|b: &f64| *b > -50.0);
    let vac: f64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(3.0);
    let (closed, open) = s14_dc(&eng, barrier, true);
    let ax = eng.fields.dc_phi_axis(closed);
    let (zmax, vmax) = ax.iter().filter(|(z, _)| *z > 12.0 && *z < 19.0).fold((0.0, f64::MIN), |a, &(z, p)| if p > a.1 { (z, p) } else { a });
    let vmin = ax.iter().filter(|(z, _)| *z > 4.0 && *z < 15.0).map(|&(_, p)| p).fold(f64::MAX, f64::min);
    println!("S14 config on axis: trap side 20 V, channel well min {:.2} V, gate max {:.2} V at z={:.1} mm (=> {:.2} V above trap)", vmin, vmax, zmax, vmax - 20.0);
    let rf = |v_ac: f64| RfAc { v_rf: 75.0, f_rf: 500e3, v_ac, f_ac: 83e3 };
    let o_open = eng.run_batch(&[Job { sp: bak2, drive: rf(0.0), n, seed: 7 }], open, &ph)[0];
    let (mut xs, mut ys) = (vec![], vec![]);
    let fs: Vec<f64> = (0..).map(|i| 70.0 + i as f64).take_while(|&f| f <= 110.0).collect();
    // the gate-closed, no-AC reference rides in the same batch (same DC): job 0
    let mut jobs: Vec<Job> = vec![Job { sp: bak2, drive: rf(0.0), n: if std::env::var("CHSIM_S14_NONOAC").is_ok() { 1 } else { n }, seed: 8 }];
    jobs.extend(fs.iter().map(|&f| Job { sp: bak2, drive: RfAc { v_rf: 75.0, f_rf: 500e3, v_ac: vac, f_ac: f * 1e3 }, n, seed: 100 + f as u64 }));
    let t0 = Instant::now();
    if std::env::var("CHSIM_WARM").is_ok() { let t = Instant::now(); let _ = eng.run_batch(&jobs, closed, &ph); eprintln!("warm-up batch: {:.3} s", t.elapsed().as_secs_f64()); }
    let t0 = Instant::now();
    let outs_all = eng.run_batch(&jobs, closed, &ph);
    let dt_scan = t0.elapsed().as_secs_f64();
    let (o_noac, outs) = (outs_all[0], &outs_all[1..]);
    println!("Bak BH3 2+ (m/z 863.0), {} Pa, {} ions/pt, per pass from trap side:", pa, n);
    println!("  gate OPEN, no AC:    transmitted {:.3}  returned {:.3}  crashed {:.3}  in-channel@{}ms {:.3}", o_open.frac(o_open.transmitted), o_open.frac(o_open.returned), o_open.frac(o_open.crashed), ph.t_max * 1e3, o_open.frac(o_open.timeout));
    println!("  gate CLOSED, no AC:  transmitted {:.3}  returned {:.3}  crashed {:.3}  in-channel {:.3}", o_noac.frac(o_noac.transmitted), o_noac.frac(o_noac.returned), o_noac.frac(o_noac.crashed), o_noac.frac(o_noac.timeout));
    println!("  gate CLOSED, AC {} V (per side) scan:", vac);
    println!("   f[kHz]  transm  return  crash  in-ch  <coll>  <t>[ms]   eff_ss=T/(T+C)");
    for (&f, o) in fs.iter().zip(outs) {
        let (tr, cr) = (o.frac(o.transmitted), o.frac(o.crashed));
        println!("   {:5.1}   {:.3}   {:.3}   {:.3}  {:.3}  {:6.1}  {:6.3}    {:.3}", f, tr, o.frac(o.returned), cr, o.frac(o.timeout), o.mean_coll, o.mean_time * 1e3, if tr + cr > 0.0 { tr / (tr + cr) } else { 0.0 });
        xs.push(f); ys.push(tr);
    }
    let (pk, py, w) = fwhm(&xs, &ys);
    println!("resonance: peak {:.1}% per pass at {:.1} kHz (their S14 setting: 83 kHz); FWHM {:.1} kHz => f/df = {:.1}", py * 100.0, pk, w, pk / w);
    println!("scan time {:.3} s ({} jobs x {} ions = {} ions in one batch, {} path)", dt_scan, jobs.len(), n, jobs.len() * n, if eng.path == Path::Tiled { "tiled" } else { "direct" });
}

fn compare_theirs(dir: &str, nions: usize) {
    let th = Fields::from_theirs(dir);
    let ours = Fields::build(&Geom::default()); // same rasterization as theirs (grow 0, Q2 flare quirk), converged
    let ours_dc: Vec<f64> = (0..ours.phi_rf.len()).map(|n| ours.phi_dc[1][n] - ours.phi_dc[2][n]).collect();
    let stats = |a: &[f64], b: &[f64]| { let mut mx = 0f64; let mut s = 0f64; for (x, y) in a.iter().zip(b) { let d = (x - y).abs(); mx = mx.max(d); s += d * d; } (mx, (s / a.len() as f64).sqrt()) };
    for (name, a, b) in [("RF", &th.phi_rf, &ours.phi_rf), ("AC", &th.phi_ac, &ours.phi_ac), ("DC", &th.phi_dc[1], &ours_dc)] {
        let (mx, rms) = stats(a, b);
        println!("{} basis: max |their - converged| = {:.4} V per V, rms {:.4}", name, mx, rms);
    }
    let a = 80usize;
    println!("on-axis DC (x5 V): z[mm]  theirs   converged");
    for k in (120..200).step_by(8) { let n = a + 160 * (a + 160 * k); println!("               {:5.1}  {:7.3}  {:7.3}", k as f64 * 0.1, 5.0 * th.phi_dc[1][n], 5.0 * ours_dc[n]); }
    let eng = Engine::new(th);
    let ph = Physics::legacy(1e-3 * 133.322);
    let dr = Drive { v_rf: 100.0, f_rf: 500e3, v_ac: 0.5, f_ac: 100e3, dc: [0.0, 5.0, 0.0, 0.0] };
    let mzs: Vec<f64> = (0..).map(|i| 485.0 + 2.0 * i as f64).take_while(|&m| m <= 515.0).collect();
    let jobs: Vec<Job> = mzs.iter().map(|&mz| Job { sp: Species { mass: mz * 5.0, z: 5.0, ccs_n2: 0.0 }, drive: dr.rfac(), n: nions, seed: 1000 + mz as u64 }).collect();
    let outs = eng.run_batch(&jobs, dr.dc, &ph);
    let (mut xs, mut ys) = (vec![], vec![]);
    for (&mz, o) in mzs.iter().zip(&outs) {
        xs.push(mz); ys.push(o.frac(o.transmitted));
        if (mz as i64) % 4 == 1 || mz == 499.0 { println!("  m/z {:5.0}: transmitted {:.3} returned {:.3} crashed {:.3} left {:.3}", mz, o.frac(o.transmitted), o.frac(o.returned), o.frac(o.crashed), o.frac(o.left)); }
    }
    let (pk, py, w) = fwhm(&xs, &ys);
    println!("OUR legacy engine on THEIR potentials: peak {:.1}% at m/z {:.0}, FWHM {:.2}, R {:.1}   (their code: 74.2%, 12.09, 41.3; ours on converged fields: 55.2%, 12.3, 41)", py * 100.0, pk, w, pk / w);
}

pub fn cmd(args: &[String]) {
    let mut args: Vec<String> = args.to_vec();
    if let Some(i) = args.iter().position(|a| a == "--untiled") { args.remove(i); DEFAULT_DIRECT.store(true, std::sync::atomic::Ordering::Relaxed); }
    let args = &args[..];
    let sub = args.get(0).map(|s| s.as_str()).unwrap_or("validate");
    if sub == "compare-theirs" { compare_theirs(args.get(1).map(|s| s.as_str()).unwrap_or("."), args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1000)); return; }
    if sub == "s14" { s14(args); return; }
    if sub == "hash" { crate::check::hash_jobs(args); return; }
    if sub == "check-tiled" { crate::check::check_tiled(args); return; }
    if sub == "field-check" {
        for &(h, al) in &[(0.1, 0.02), (0.05, 0.02), (0.1, 0.0), (0.05, 0.0)] { field_check(h, al); }
        return;
    }
    let geom = Geom::default();
    let fields = Fields::build(&geom);
    match sub {
        "field" => {
            r0_profile(&fields);
            let ax = fields.dc_phi_axis([0.0, 5.0, -5.0, 0.0]);
            println!("Fig.4 DC config on axis (Q1 0, plate +5, downstream -5, box 0):");
            for (z, p) in ax.iter().step_by(10) { println!("   z {:5.1} mm  phi {:7.3} V", z, p); }
        }
        "secular" => {
            let bak2 = Species { mass: 1725.95, z: 2.0, ccs_n2: 425.0 };
            for &(zmm, label) in &[(5.0, "near entrance"), (9.0, "mid"), (13.0, "near exit")] {
                let fc = secular_freq(&fields, &bak2, 75.0, 500e3, zmm, false);
                let fl = secular_freq(&fields, &bak2, 75.0, 500e3, zmm, true);
                println!("m/z 863 @ 75 V, z = {:4.1} mm ({}): f_sec correct {:6.2} kHz | legacy half-force {:6.2} kHz | S14 used 83 kHz", zmm, label, fc / 1e3, fl / 1e3);
            }
            let m500 = Species { mass: 500.0, z: 1.0, ccs_n2: 230.0 };
            let fc = secular_freq(&fields, &m500, 100.0, 500e3, 9.0, false);
            let fl = secular_freq(&fields, &m500, 100.0, 500e3, 9.0, true);
            println!("m/z 500 @ 100 V (Fig. 4 point), mid: correct {:.2} kHz (expected unstable) | legacy {:.2} kHz (their 100 kHz)", fc / 1e3, fl / 1e3);
        }
        "legacy-fig4" => {
            let eng = Engine::new(fields);
            let ph = Physics::legacy(1e-3 * 133.322);
            let dr = Drive { v_rf: 100.0, f_rf: 500e3, v_ac: 0.5, f_ac: 100e3, dc: [0.0, 5.0, -5.0, 0.0] };
            let nions: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(1000);
            println!("Legacy replication of Fig. 4 (m/z scan, {} ions/point, 5+ ions of m/z = mass/5):", nions);
            println!("  m/z    returned  transmitted  collided  left  timeout   <coll>   <t>[ms]");
            let mzs: Vec<f64> = (0..).map(|i| 475.0 + 2.0 * i as f64).take_while(|&m| m <= 525.0).collect();
            let jobs: Vec<Job> = mzs.iter().map(|&mz| Job { sp: Species { mass: mz * 5.0, z: 5.0, ccs_n2: 0.0 }, drive: dr.rfac(), n: nions, seed: 1000 + mz as u64 }).collect();
            let t0 = Instant::now();
            let outs = eng.run_batch(&jobs, dr.dc, &ph);
            let dt_scan = t0.elapsed().as_secs_f64();
            let (mut xs, mut ys) = (vec![], vec![]);
            for (&mz, o) in mzs.iter().zip(&outs) {
                xs.push(mz); ys.push(o.frac(o.transmitted));
                println!("  {:5.0}  {:7.3}   {:7.3}     {:7.3}  {:5.3}  {:5.3}   {:7.1}  {:7.3}", mz, o.frac(o.returned), o.frac(o.transmitted), o.frac(o.crashed), o.frac(o.left), o.frac(o.timeout), o.mean_coll, o.mean_time * 1e3);
            }
            let (pk, py, w) = fwhm(&xs, &ys);
            println!("peak {:.1}% at m/z {:.0}; FWHM {:.1} m/z => R = {:.0}   (their Fig. 4: ~80%, R ~ 50)", py * 100.0, pk, w, pk / w);
            println!("scan time {:.3} s ({} path)", dt_scan, if eng.path == Path::Tiled { "tiled" } else { "direct" });
        }
        _ => eprintln!("channel subcommands: field | secular | legacy-fig4 [nions] | s14 [Pa n barrier vac] | check-tiled [n] | compare-theirs   (global flag: --untiled)"),
    }
}
