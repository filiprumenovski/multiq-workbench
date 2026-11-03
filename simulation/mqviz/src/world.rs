//! Packed field world and the f64 ion integrator (a faithful port of chsim `xtalk.rs` trace physics).
//!
//! Integrator (identical to xtalk-trace): f64 kick-drift leapfrog, dt = T_RF/50, trilinear fields (RF, AC, DC),
//! hard-sphere collisions against a Maxwellian N2 bath with real N2 CCS (null-collision with bound v + 6 s_gas),
//! crash = nearest node is conductor (or box exit), transmitted at z > 18.6 mm, returned at z < 0.5 mm.
//! Differences from xtalk.rs are bookkeeping only: the RF/AC phase of every ion is a per-ion offset on global
//! phase accumulators (never accumulated in f32), so a live frequency change is phase-continuous and one
//! sin/cos per block replaces one per ion-step.

use chsim::channel::{s14_geom, Fields, Geom};
use rayon::prelude::*;

pub const AMU: f64 = 1.660_539_066_6e-27;
pub const QE: f64 = 1.602_176_634e-19;
pub const KB: f64 = 1.380_649e-23;
pub const TWO_PI: f64 = std::f64::consts::TAU;

pub const F_RF: f64 = 500e3;
pub const OMEGA: f64 = TWO_PI * F_RF;
pub const DT: f64 = 1.0 / F_RF / 50.0; // T_RF/50 = 4e-8 s
pub const Z_EXIT: f64 = 18.6e-3;
pub const Z_RETURN: f64 = 0.5e-3;
pub const Z_INJECT: f64 = 1.0e-3;
pub const R_INJECT: f64 = 2.0e-3; // chsim Physics::correct(..).r_inject_mm
pub const T_MAX: f64 = 3.0e-3;
pub const T_BATH: f64 = 300.0;
pub const M_GAS: f64 = 28.0 * AMU;
pub const R_EXCITED: f64 = 1.2e-3; // |r| beyond which an ion counts as resonantly excited
/// Steps per parallel block (5.12 us of simulated time).
pub const BLOCK: usize = 128;
pub const NTRAIL: usize = 4;
/// Calibrated DC set [channel quads, gate plate, downstream, trap/box] (fig. S14).
pub const DC_CAL: [f64; 4] = [18.0, 24.0, 15.0, 20.0];

pub const ST_ALIVE: u8 = 0;
pub const ST_TRANS: u8 = 1;
pub const ST_CRASH: u8 = 2;
pub const ST_RETURN: u8 = 3;
pub const ST_TIMEOUT: u8 = 4;

// ------------------------------------------------------------------ packed field node

/// One grid node, one cache line: RF E per volt (+ conductor mask in [3]), AC E per volt, combined DC E.
#[repr(C, align(64))]
#[derive(Clone, Copy, Default)]
pub struct Node {
    pub rf: [f32; 4],
    pub ac: [f32; 4],
    pub dc: [f32; 4],
    pub pad: [f32; 4],
}

/// Tables for the topographic view (x-z slice through the axis) and the end view.
pub struct VizTables {
    pub nx: usize,
    pub nz: usize,
    pub ny: usize,
    pub h_mm: f64,
    pub axis_mm: f64,
    /// Conductor mask on the y = axis plane, index i + nx*k.
    pub cut: Vec<bool>,
    /// Projection over y of conductor nodes (rods only meaningful; plate holes are hidden by it).
    pub shadow: Vec<bool>,
    /// |E_rf|^2 per V^2 on the y = axis plane.
    pub e2: Vec<f32>,
    /// DC basis potentials per volt on the plane.
    pub phi_dc: [Vec<f32>; 4],
    /// Conductor mask in the x-y plane at the end-view z, index i + nx*j.
    pub end_mask: Vec<bool>,
    pub end_z_mm: f64,
    pub layers: [f64; 8],
}

/// Per-z principal-axis gradient tables on the axis (Jacobian method of xtalk-freq with the principal-axis fix).
pub struct TerrainTables {
    pub z_mm: Vec<f64>,
    /// RF principal gradients [gx, gy] per z (V/m^2 per volt).
    pub rf: Vec<[f64; 2]>,
    /// DC principal gradients per basis, per z: [basis][z] = [gx, gy] (V/m^2 per volt of that basis).
    pub dc: [Vec<[f64; 2]>; 4],
}

pub struct World {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub h: f64,
    pub nodes: Vec<Node>,
    pub dc_bases: [Vec<[f32; 4]>; 4],
    pub a: f64,
    pub lim: [f64; 3],
    pub dc_volts: [f64; 4],
    pub viz: std::sync::Arc<VizTables>,
    pub terr: std::sync::Arc<TerrainTables>,
    pub geom: Geom,
}

/// Combined static DC field node value (same f32 formula and order as `Fields::dc_field`).
#[inline]
fn dc_combine(bases: &[Vec<[f32; 4]>; 4], v: [f64; 4], n: usize) -> [f32; 4] {
    let mut e = [0f32; 4];
    for b in 0..4 {
        for d in 0..3 {
            e[d] += (v[b] as f32) * bases[b][n][d];
        }
    }
    e
}

/// Trilinear interpolation exactly as xtalk.rs `tri` (clamped cell index, unclamped weights).
#[inline(always)]
fn tri_one(nodes: &[Node], nx: usize, ny: usize, nz: usize, h: f64, x: [f64; 3], sel: u8) -> [f64; 3] {
    let g = [x[0] / h, x[1] / h, x[2] / h];
    let i = (g[0].max(0.0) as usize).min(nx - 2);
    let j = (g[1].max(0.0) as usize).min(ny - 2);
    let k = (g[2].max(0.0) as usize).min(nz - 2);
    let (tx, ty, tz) = (g[0] - i as f64, g[1] - j as f64, g[2] - k as f64);
    let mut e = [0.0; 3];
    for (di, wx) in [(0usize, 1.0 - tx), (1, tx)] {
        for (dj, wy) in [(0usize, 1.0 - ty), (1, ty)] {
            for (dk, wz) in [(0usize, 1.0 - tz), (1, tz)] {
                let c = &nodes[(i + di) + nx * ((j + dj) + ny * (k + dk))];
                let src = match sel { 0 => &c.rf, 1 => &c.ac, _ => &c.dc };
                for d in 0..3 {
                    e[d] += wx * wy * wz * src[d] as f64;
                }
            }
        }
    }
    e
}

impl World {
    /// Solve (or load from the field cache) the calibrated S14 geometry and pack it.
    pub fn build_calibrated() -> World {
        let gargs: Vec<String> = vec!["".into(), "".into(), "".into(), "".into(), "".into(), "code".into(), "7.0".into()];
        let g = s14_geom(&gargs);
        let f = Fields::build(&g);
        World::from_fields(f, DC_CAL)
    }

    pub fn from_fields(f: Fields, dc: [f64; 4]) -> World {
        let (nx, ny, nz, h) = (f.nx, f.ny, f.nz, f.h);
        let a = f.geom.axis() * 1e-3;
        let lim = [(nx - 1) as f64 * h, (ny - 1) as f64 * h, (nz - 1) as f64 * h];
        let n = nx * ny * nz;
        // ---- viz tables (before the heavy arrays are dropped)
        let j0 = ny / 2;
        let idx = |i: usize, j: usize, k: usize| i + nx * (j + ny * k);
        let mut cut = vec![false; nx * nz];
        let mut shadow = vec![false; nx * nz];
        let mut e2 = vec![0f32; nx * nz];
        let mut phi_dc: [Vec<f32>; 4] = [vec![0f32; nx * nz], vec![0f32; nx * nz], vec![0f32; nx * nz], vec![0f32; nx * nz]];
        for k in 0..nz {
            for i in 0..nx {
                let c = f.rf[idx(i, j0, k)];
                cut[i + nx * k] = c[3] > 0.5;
                e2[i + nx * k] = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
                for b in 0..4 {
                    phi_dc[b][i + nx * k] = f.phi_dc[b][idx(i, j0, k)] as f32;
                }
                shadow[i + nx * k] = (0..ny).any(|j| f.rf[idx(i, j, k)][3] > 0.5);
            }
        }
        let layers = {
            let g = &f.geom;
            let l1 = g.q1_start;
            let l2 = l1 + g.q1_len;
            let l3 = l2 + g.gap;
            let l4 = l3 + g.plate;
            let l5 = l4 + g.gap;
            let l6 = l5 + g.q2_len;
            let l7 = l6 + g.gap;
            let l8 = l7 + g.plate;
            [l1, l2, l3, l4, l5, l6, l7, l8]
        };
        let end_z_mm = 9.0;
        let ke = (end_z_mm / (h * 1e3)).round() as usize;
        let mut end_mask = vec![false; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                end_mask[i + nx * j] = f.rf[idx(i, j, ke)][3] > 0.5;
            }
        }
        let viz = VizTables { nx, nz, ny, h_mm: h * 1e3, axis_mm: f.geom.axis(), cut, shadow, e2, phi_dc, end_mask, end_z_mm, layers };

        // ---- terrain tables (principal-axis Jacobians on the axis)
        let terr = {
            let hh = 0.2e-3;
            let dcs: Vec<Vec<[f32; 4]>> = (0..4).map(|b| { let mut w = [0.0; 4]; w[b] = 1.0; f.dc_field(w) }).collect();
            let trif = |v: &[[f32; 4]], x: [f64; 3]| -> [f64; 3] {
                let gg = [x[0] / h, x[1] / h, x[2] / h];
                let i = (gg[0].max(0.0) as usize).min(nx - 2);
                let j = (gg[1].max(0.0) as usize).min(ny - 2);
                let k = (gg[2].max(0.0) as usize).min(nz - 2);
                let (tx, ty, tz) = (gg[0] - i as f64, gg[1] - j as f64, gg[2] - k as f64);
                let mut e = [0.0; 3];
                for (di, wx) in [(0, 1.0 - tx), (1, tx)] { for (dj, wy) in [(0, 1.0 - ty), (1, ty)] { for (dk, wz) in [(0, 1.0 - tz), (1, tz)] {
                    let c = v[(i + di) + nx * ((j + dj) + ny * (k + dk))];
                    for d in 0..3 { e[d] += wx * wy * wz * c[d] as f64; } } } }
                e
            };
            let mut z_mm = vec![];
            let mut rf = vec![];
            let mut dcv: [Vec<[f64; 2]>; 4] = [vec![], vec![], vec![], vec![]];
            let mut z = 0.3e-3;
            while z < (nz - 2) as f64 * h && z < 20.0e-3 {
                let jac = |v: &[[f32; 4]]| {
                    let mut jm = [[0.0f64; 2]; 2];
                    for dj in 0..2 {
                        let mut p = [a, a, z];
                        let mut m = p;
                        p[dj] += hh;
                        m[dj] -= hh;
                        let (ep, em) = (trif(v, p), trif(v, m));
                        for di in 0..2 { jm[di][dj] = (ep[di] - em[di]) / (2.0 * hh); }
                    }
                    jm
                };
                let jr = jac(&f.rf);
                let (s11, s12, s22) = (jr[0][0], 0.5 * (jr[0][1] + jr[1][0]), jr[1][1]);
                let th = 0.5 * (2.0 * s12).atan2(s11 - s22);
                let (c, sn) = (th.cos(), th.sin());
                let proj = |j: [[f64; 2]; 2], ux: f64, uy: f64| ux * (j[0][0] * ux + j[0][1] * uy) + uy * (j[1][0] * ux + j[1][1] * uy);
                z_mm.push(z * 1e3);
                rf.push([proj(jr, c, sn), proj(jr, -sn, c)]);
                for b in 0..4 {
                    let jd = jac(&dcs[b]);
                    dcv[b].push([proj(jd, c, sn), proj(jd, -sn, c)]);
                }
                z += 0.1e-3;
            }
            TerrainTables { z_mm, rf, dc: dcv }
        };

        // ---- pack the nodes (parallel), keep only the DC bases for later gate changes
        let mut nodes = vec![Node::default(); n];
        let dc_bases = f.dc;
        {
            let (rf, ac) = (&f.rf, &f.ac);
            nodes.par_iter_mut().enumerate().for_each(|(i, nd)| {
                nd.rf = rf[i];
                nd.ac = ac[i];
                nd.dc = dc_combine(&dc_bases, dc, i);
            });
        }
        let geom = f.geom.clone();
        drop(f.rf);
        drop(f.ac);
        drop(f.phi_rf);
        drop(f.phi_ac);
        drop(f.phi_dc);
        World { nx, ny, nz, h, nodes, dc_bases, a, lim, dc_volts: dc, viz: std::sync::Arc::new(viz), terr: std::sync::Arc::new(terr), geom }
    }

    /// Re-combine the static DC field for new electrode voltages (parallel; ~0.1 s).
    pub fn set_dc(&mut self, v: [f64; 4]) {
        if v == self.dc_volts { return; }
        let bases = &self.dc_bases;
        self.nodes.par_iter_mut().enumerate().for_each(|(i, nd)| nd.dc = dc_combine(bases, v, i));
        self.dc_volts = v;
    }

    #[inline(always)]
    pub fn fields_at(&self, x: [f64; 3]) -> ([f64; 3], [f64; 3], [f64; 3]) {
        // One pass over the 8 corner nodes for all three fields (same arithmetic order as three `tri` calls).
        let h = self.h;
        let g = [x[0] / h, x[1] / h, x[2] / h];
        let i = (g[0].max(0.0) as usize).min(self.nx - 2);
        let j = (g[1].max(0.0) as usize).min(self.ny - 2);
        let k = (g[2].max(0.0) as usize).min(self.nz - 2);
        let (tx, ty, tz) = (g[0] - i as f64, g[1] - j as f64, g[2] - k as f64);
        let (mut er, mut ea, mut ed) = ([0.0; 3], [0.0; 3], [0.0; 3]);
        for (di, wx) in [(0usize, 1.0 - tx), (1, tx)] {
            for (dj, wy) in [(0usize, 1.0 - ty), (1, ty)] {
                for (dk, wz) in [(0usize, 1.0 - tz), (1, tz)] {
                    let c = &self.nodes[(i + di) + self.nx * ((j + dj) + self.ny * (k + dk))];
                    let w = wx * wy * wz;
                    for d in 0..3 {
                        er[d] += w * c.rf[d] as f64;
                        ea[d] += w * c.ac[d] as f64;
                        ed[d] += w * c.dc[d] as f64;
                    }
                }
            }
        }
        (er, ea, ed)
    }

    /// Gate-voltage adjustment for a blocked port without duplicating the channel grid.
    pub fn gate_delta(&self,x:[f64;3],volts:f64)->[f64;3] {
        let g=x.map(|v|v/self.h); let dims=[self.nx,self.ny,self.nz];
        let c:[usize;3]=std::array::from_fn(|d|(g[d].max(0.) as usize).min(dims[d]-2));
        let t:[f64;3]=std::array::from_fn(|d|g[d]-c[d] as f64);
        let mut out=[0.;3];
        for i in 0..2 {for j in 0..2 {for k in 0..2 {
          let q=c[0]+i+self.nx*(c[1]+j+self.ny*(c[2]+k));
          let w=volts * (if i==1{t[0]}else{1.-t[0]}) * (if j==1{t[1]}else{1.-t[1]}) * (if k==1{t[2]}else{1.-t[2]});
          for d in 0..3 {out[d]+=w*self.dc_bases[1][q][d] as f64;}
        }}}
        out
    }
    /// Single-field trilinear (0 = RF, 1 = AC, 2 = DC), for diagnostics.
    pub fn tri(&self, x: [f64; 3], sel: u8) -> [f64; 3] {
        tri_one(&self.nodes, self.nx, self.ny, self.nz, self.h, x, sel)
    }

    #[inline(always)]
    pub fn is_crash(&self, x: [f64; 3]) -> bool {
        if (0..3).any(|d| x[d] < 0.0 || x[d] > self.lim[d]) { return true; }
        let k = ((x[0] / self.h).round() as usize).min(self.nx - 1)
            + self.nx * (((x[1] / self.h).round() as usize).min(self.ny - 1) + self.ny * ((x[2] / self.h).round() as usize).min(self.nz - 1));
        self.nodes[k].rf[3] > 0.5
    }
}

// ------------------------------------------------------------------ species and per-step constants

/// A simulated species: neutral mass [Da] (includes charging protons), charge number, N2 CCS [A^2].
#[derive(Clone, Debug, PartialEq)]
pub struct Sp {
    pub label: String,
    pub mass: f64,
    pub z: f64,
    pub ccs: f64,
    /// Relative source abundance (normalised by the sampler).
    pub weight: f64,
}
impl Sp {
    pub fn mz(&self) -> f64 { self.mass / self.z }
}

#[derive(Clone, Copy, Debug)]
pub struct SpPhys {
    pub qm: f64,
    pub s_ion: f64,
    pub mr: f64,
    pub sig: f64,
    pub mass_kg: f64,
    pub z: f64,
}
impl SpPhys {
    pub fn new(sp: &Sp) -> SpPhys {
        let m = sp.mass * AMU;
        SpPhys { qm: sp.z * QE / m, s_ion: (KB * T_BATH / m).sqrt(), mr: M_GAS / (m + M_GAS), sig: sp.ccs * 1e-20, mass_kg: m, z: sp.z }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct StepCfg {
    pub dt: f64,
    pub vrf: f64,
    pub vac: f64,
    pub ngas: f64,
    pub s_gas: f64,
    pub tmax_steps: u32,
}
impl StepCfg {
    pub fn new(vrf: f64, vac: f64, pressure_pa: f64) -> StepCfg {
        StepCfg { dt: DT, vrf, vac, ngas: pressure_pa / (KB * T_BATH), s_gas: (KB * T_BATH / M_GAS).sqrt(), tmax_steps: (T_MAX / DT).round() as u32 }
    }
}

// ------------------------------------------------------------------ ion

#[derive(Clone)]
pub struct Ion {
    pub x: [f64; 3],
    pub v: [f64; 3],
    pub rng: u64,
    pub rs: f64, // sin / cos of the RF phase offset
    pub rc: f64,
    pub as_: f64, // sin / cos of the AC phase offset
    pub ac: f64,
    pub age: u32,
    pub id: u32,
    pub sp: u8,
    pub tab: u8,
    pub state: u8,
    pub fin_step: u32, // step index within the last block at which the ion finished
    pub steps_blk: u32,
    pub trail: [[f32; 3]; NTRAIL],
}

#[inline(always)]
fn xs_next(st: &mut u64) -> f64 {
    *st ^= *st << 13;
    *st ^= *st >> 7;
    *st ^= *st << 17;
    (*st >> 11) as f64 / (1u64 << 53) as f64
}
#[inline(always)]
fn gs(u1: f64, u2: f64) -> f64 {
    (-2.0 * u1.max(1e-300).ln()).sqrt() * (TWO_PI * u2).cos()
}

pub fn splitmix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E3779B97F4A7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

impl Ion {
    /// Effusive thermal injection as xtalk-trace: uniform disc r <= 2 mm at z = 1 mm, Maxwellian transverse,
    /// flux-weighted (Rayleigh) axial velocity, random RF and AC phase (`rf_phase0`: RF phase 0 for every ion, the convention of chsim's Metal engine).
    pub fn inject(seed: u64, a: f64, sp_idx: u8, sp: &SpPhys, id: u32, tab: u8, rf_phase0: bool) -> Ion {
        let mut st = seed | 1;
        let mut u = || xs_next(&mut st);
        let r = R_INJECT * u().sqrt();
        let th = TWO_PI * u();
        let x = [a + r * th.cos(), a + r * th.sin(), Z_INJECT];
        let v = [sp.s_ion * gs(u(), u()), sp.s_ion * gs(u(), u()), sp.s_ion * (-2.0 * u().max(1e-300).ln()).sqrt()];
        let ph = TWO_PI * u();
        let ph_rf = TWO_PI * u();
        let ph_rf = if rf_phase0 { 0.0 } else { ph_rf };
        Ion {
            x, v, rng: st | 1, rs: ph_rf.sin(), rc: ph_rf.cos(), as_: ph.sin(), ac: ph.cos(),
            age: 0, id, sp: sp_idx, tab, state: ST_ALIVE, fin_step: 0, steps_blk: 0,
            trail: [[x[0] as f32 * 1e3, x[1] as f32 * 1e3, x[2] as f32 * 1e3]; NTRAIL],
        }
    }

    #[inline]
    pub fn r_mm(&self, a: f64) -> f64 { (((self.x[0] - a).powi(2) + (self.x[1] - a).powi(2)).sqrt()) * 1e3 }
}

/// Phase tables for one block: sin/cos of the global RF phase and of every AC tone for steps 0..BLOCK.
pub struct Tabs {
    pub rf_s: [f64; BLOCK],
    pub rf_c: [f64; BLOCK],
    pub ac_s: Vec<f64>, // [tab * BLOCK + s]
    pub ac_c: Vec<f64>,
}

/// Global phase accumulators (f64, wrapped each block).
pub struct Clock {
    pub rf: f64,
    pub ac: Vec<f64>,
    pub f_ac: Vec<f64>, // Hz per tone (tab 0 = live)
}
impl Clock {
    pub fn new(ntab: usize, f0: f64) -> Clock { Clock { rf: 0.0, ac: vec![0.0; ntab], f_ac: vec![f0; ntab] } }
    /// Fill `tabs` for the next block starting at the current phases, then advance the phases by one block.
    pub fn fill(&mut self, tabs: &mut Tabs) {
        let w = OMEGA * DT;
        for s in 0..BLOCK {
            let p = self.rf + w * s as f64;
            tabs.rf_s[s] = p.sin();
            tabs.rf_c[s] = p.cos();
        }
        self.rf = (self.rf + w * BLOCK as f64) % TWO_PI;
        let nt = self.ac.len();
        tabs.ac_s.resize(nt * BLOCK, 0.0);
        tabs.ac_c.resize(nt * BLOCK, 0.0);
        for t in 0..nt {
            let wa = TWO_PI * self.f_ac[t] * DT;
            for s in 0..BLOCK {
                let p = self.ac[t] + wa * s as f64;
                tabs.ac_s[t * BLOCK + s] = p.sin();
                tabs.ac_c[t * BLOCK + s] = p.cos();
            }
            self.ac[t] = (self.ac[t] + wa * BLOCK as f64) % TWO_PI;
        }
    }
}
impl Tabs {
    pub fn new() -> Tabs { Tabs { rf_s: [0.0; BLOCK], rf_c: [0.0; BLOCK], ac_s: vec![], ac_c: vec![] } }
}

/// Shared nitrogen bath; nsd = number density × cross section × time step.
#[inline]
pub fn collide(v: &mut [f64;3], rng: &mut u64, sp: &SpPhys, cfg: &StepCfg, nsd: f64) {
        let vm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        let bound = vm + 6.0 * cfg.s_gas;
        if xs_next(rng) < nsd * bound {
            let vg = [cfg.s_gas * gs(xs_next(rng), xs_next(rng)), cfg.s_gas * gs(xs_next(rng), xs_next(rng)), cfg.s_gas * gs(xs_next(rng), xs_next(rng))];
            let gv = [v[0] - vg[0], v[1] - vg[1], v[2] - vg[2]];
            let gl = (gv[0] * gv[0] + gv[1] * gv[1] + gv[2] * gv[2]).sqrt();
            if xs_next(rng) * bound < gl {
                let c = 2.0 * xs_next(rng) - 1.0;
                let p2 = TWO_PI * xs_next(rng);
                let sn = (1.0 - c * c).max(0.0).sqrt();
                let nn = [sn * p2.cos(), sn * p2.sin(), c];
                for d in 0..3 { v[d] = v[d] - sp.mr * gv[d] + sp.mr * gl * nn[d]; }
            }
        }
}

/// Advance one ion by up to `nsteps` steps (<= BLOCK) using the phase tables of the current block.
/// Sets `state` and `fin_step` when the ion finishes. This is the hot loop.
#[inline(never)]
pub fn advance(ion: &mut Ion, w: &World, sp: &SpPhys, cfg: &StepCfg, tabs: &Tabs) {
    advance_with_field(ion, w, sp, cfg, tabs, |x| w.fields_at(x))
}

/// Replay with an independently refined force grid, retaining the original collision mask and box.
pub fn advance_with_field(ion: &mut Ion, w: &World, sp: &SpPhys, cfg: &StepCfg, tabs: &Tabs,
    fields: impl Fn([f64; 3]) -> ([f64; 3], [f64; 3], [f64; 3])) {
    let (dt, qm) = (cfg.dt, sp.qm);
    let nsd = cfg.ngas * sp.sig * dt;
    let a = w.a;
    let off = ion.tab as usize * BLOCK;
    let use_ac = cfg.vac != 0.0;
    let mut x = ion.x;
    let mut v = ion.v;
    let mut rng = ion.rng;
    let mut age = ion.age;
    let mut done = 0u8;
    let mut fin = BLOCK as u32;
    let mut steps = 0u32;
    for s in 0..BLOCK {
        if age >= cfg.tmax_steps { done = ST_TIMEOUT; fin = s as u32; break; }
        let srf = cfg.vrf * (tabs.rf_s[s] * ion.rc + tabs.rf_c[s] * ion.rs);
        let sac = if use_ac { cfg.vac * (tabs.ac_s[off + s] * ion.ac + tabs.ac_c[off + s] * ion.as_) } else { 0.0 };
        let (er, ea, ed) = fields(x);
        for d in 0..3 {
            v[d] += qm * dt * (srf * er[d] + sac * ea[d] + ed[d]);
            x[d] += v[d] * dt;
        }
        age += 1;
        steps += 1;
        collide(&mut v, &mut rng, sp, cfg, nsd);
        if x[2] > Z_EXIT { done = ST_TRANS; fin = s as u32 + 1; break; }
        if x[2] < Z_RETURN { done = ST_RETURN; fin = s as u32 + 1; break; }
        if w.is_crash(x) { done = ST_CRASH; fin = s as u32 + 1; break; }
    }
    ion.x = x;
    ion.v = v;
    ion.rng = rng;
    ion.age = age;
    ion.state = done;
    ion.fin_step = fin;
    ion.steps_blk = steps;
    let _ = a;
}

/// Run a block over all ions (rayon, small chunks for work stealing across P and E cores).
pub fn advance_all(pool: &rayon::ThreadPool, ions: &mut [Ion], w: &World, sps: &[SpPhys], cfg: &StepCfg, tabs: &Tabs) {
    pool.install(|| {
        ions.par_chunks_mut(16).for_each(|c| {
            for ion in c.iter_mut() {
                if ion.state == ST_ALIVE {
                    advance(ion, w, &sps[ion.sp as usize], cfg, tabs);
                }
            }
        });
    });
}

pub fn make_pool(n: usize) -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new().num_threads(n.max(1)).thread_name(|i| format!("mqsim-{i}")).build().unwrap()
}

pub fn default_threads() -> usize {
    let n = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(8);
    if n > 2 { n - 1 } else { n }
}
