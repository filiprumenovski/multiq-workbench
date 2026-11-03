//! Field solves for the cube: one RF basis (+-1 V checkerboard) and three DC bases (wall plane, entrance rods, exit rods).
//! The post basis is implicit: the four DC bases sum to 1 everywhere (outer box value belongs to the wall basis), so
//! Phi_dc = U0 + dU_wall*phi_wall + (q_in-U0)*phi_ent + (q_out-U0)*phi_ext  (no large-offset cancellation in float32).

use crate::geometry::*;
use fieldsolve::{solve, Grid, Problem};
use rayon::prelude::*;
use std::time::Instant;

pub struct Fields {
    pub geom: Geometry,
    pub h_mm: f64,
    pub grow: f64,            // electrode growth in units of h
    pub lo: [f64; 3],         // grid origin in cube coordinates [mm]
    pub n: [usize; 3],
    pub rf: Vec<[f32; 4]>,    // (Ex,Ey,Ez [V/m per V], id as f32)
    pub phi_rf: Vec<f32>,
    pub dc: [Vec<[f32; 4]>; 3], // wall, ent, ext bases: (E, phi)
    pub id: Vec<u8>,
    pub sign: Vec<i8>,
    pub solve_s: f64,
}

impl Fields {
    pub fn len(&self) -> usize { self.n[0] * self.n[1] * self.n[2] }
    pub fn idx(&self, i: usize, j: usize, k: usize) -> usize { i + self.n[0] * (j + self.n[1] * k) }

    pub fn build(geom: &Geometry, h_mm: f64, grow: f64, tol: f64) -> Fields {
        let (lo, hi) = geom.domain(geom.cap_gap + 1.0);
        let n = [0, 1, 2].map(|a| ((hi[a] - lo[a]) / h_mm).round() as usize + 1);
        let grid = Grid { nx: n[0], ny: n[1], nz: n[2], h: h_mm * 1e-3, origin: [lo[0] * 1e-3, lo[1] * 1e-3, lo[2] * 1e-3] };
        let t0 = Instant::now();
        // classification once
        let len = grid.len();
        let mut id = vec![0u8; len];
        let mut sign = vec![0i8; len];
        let plane = n[0] * n[1];
        id.par_chunks_mut(plane).zip(sign.par_chunks_mut(plane)).enumerate().for_each(|(k, (idc, sgc))| {
            let z = lo[2] + h_mm * k as f64;
            for j in 0..n[1] { for i in 0..n[0] {
                let (a, s) = geom.classify([lo[0] + h_mm * i as f64, lo[1] + h_mm * j as f64, z], grow * h_mm);
                idc[i + n[0] * j] = a; sgc[i + n[0] * j] = s;
            }}
        });
        eprintln!("  classified {} nodes ({}x{}x{} @ {} mm) in {:.2}s", len, n[0], n[1], n[2], h_mm, t0.elapsed().as_secs_f64());
        let mk = |val: &(dyn Fn(u8, i8) -> Option<f64> + Sync)| -> Problem {
            let mut mask = vec![0u8; len];
            let mut fixed = vec![0.0f64; len];
            for q in 0..len {
                if let Some(v) = val(id[q], sign[q]) { mask[q] = 1; fixed[q] = v; }
            }
            Problem { grid, mask, fixed }
        };
        let rf_p = mk(&|i, s| match i { ID_POST | ID_ENT | ID_EXT => Some(s as f64), ID_PLANE => Some(0.0), _ => None });
        let wall_p = mk(&|i, _| if i == ID_FREE { None } else { Some(if i == ID_PLANE { 1.0 } else { 0.0 }) });
        let ent_p = mk(&|i, _| if i == ID_FREE { None } else { Some(if i == ID_ENT { 1.0 } else { 0.0 }) });
        let ext_p = mk(&|i, _| if i == ID_FREE { None } else { Some(if i == ID_EXT { 1.0 } else { 0.0 }) });
        let solve_one = |p: &Problem, outer: f64, name: &str| {
            let t = Instant::now();
            let (sol, st) = solve(p, outer, tol);
            eprintln!("  solved {:<5} {} cycles, residual {:.1e}, {:.2}s", name, st.iterations, st.final_residual, t.elapsed().as_secs_f64());
            assert!(st.converged, "{name} field did not converge");
            sol
        };
        let sol_rf = solve_one(&rf_p, 0.0, "rf");
        let mut rf = sol_rf.efield_f32();
        let phi_rf: Vec<f32> = rf.iter().map(|e| e[3]).collect();
        rf.par_iter_mut().zip(id.par_iter()).for_each(|(e, &i)| e[3] = i as f32);
        drop(sol_rf);
        let wall = solve_one(&wall_p, 1.0, "wall").efield_f32();
        let ent = solve_one(&ent_p, 0.0, "ent").efield_f32();
        let ext = solve_one(&ext_p, 0.0, "ext").efield_f32();
        Fields { geom: geom.clone(), h_mm, grow, lo, n, rf, phi_rf, dc: [wall, ent, ext], id, sign, solve_s: t0.elapsed().as_secs_f64() }
    }

    /// Combined DC field for post potential u0 (irrelevant for the force), wall potential, entrance and exit rod potentials.
    pub fn dc_combined(&self, u0: f64, uwall: f64, q_in: f64, q_out: f64) -> Vec<[f32; 4]> {
        let w = [uwall - u0, q_in - u0, q_out - u0];
        let mut out = vec![[0f32; 4]; self.len()];
        out.par_iter_mut().enumerate().for_each(|(q, o)| {
            let mut e = [0f64; 4];
            for b in 0..3 { for c in 0..4 { e[c] += w[b] * self.dc[b][q][c] as f64; } }
            *o = [e[0] as f32, e[1] as f32, e[2] as f32, (e[3] + u0) as f32];
        });
        out
    }

    /// Trilinear sample of a node-array component at cube coordinates (mm).
    pub fn grid_coord(&self, p: [f64; 3]) -> ([usize; 3], [f64; 3]) {
        let mut c = [0usize; 3]; let mut t = [0.0; 3];
        for a in 0..3 {
            let f = ((p[a] - self.lo[a]) / self.h_mm).clamp(0.0, (self.n[a] - 1) as f64 - 1e-9);
            c[a] = f.floor() as usize; t[a] = f - c[a] as f64;
        }
        (c, t)
    }
    pub fn trilinear(&self, arr: &dyn Fn(usize) -> f64, p: [f64; 3]) -> f64 {
        let (c, t) = self.grid_coord(p);
        let mut s = 0.0;
        for dk in 0..2 { for dj in 0..2 { for di in 0..2 {
            let w = (if di == 1 { t[0] } else { 1.0 - t[0] }) * (if dj == 1 { t[1] } else { 1.0 - t[1] }) * (if dk == 1 { t[2] } else { 1.0 - t[2] });
            s += w * arr(self.idx(c[0] + di, c[1] + dj, c[2] + dk));
        }}}
        s
    }
    /// RF |E| per volt [V/m] at point (trilinear nodal field).
    pub fn rf_e(&self, p: [f64; 3]) -> [f64; 3] {
        [0, 1, 2].map(|c| self.trilinear(&|q| self.rf[q][c] as f64, p))
    }
    pub fn rf_phi(&self, p: [f64; 3]) -> f64 { self.trilinear(&|q| self.phi_rf[q] as f64, p) }
    /// DC potential relative to the posts (V per volt of each basis)
    pub fn dc_phi(&self, b: usize, p: [f64; 3]) -> f64 { self.trilinear(&|q| self.dc[b][q][3] as f64, p) }
    pub fn dc_e(&self, b: usize, p: [f64; 3]) -> [f64; 3] { [0, 1, 2].map(|c| self.trilinear(&|q| self.dc[b][q][c] as f64, p)) }
}
