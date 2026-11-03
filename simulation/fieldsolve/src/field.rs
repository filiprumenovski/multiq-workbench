//! Electric field evaluation from a solved potential.
//!
//! `E = -grad(phi)` is computed per node with second-order differences:
//!
//! * **central** differences when both axis neighbours are usable,
//! * **one-sided** differences (second order when two usable nodes exist on
//!   that side, otherwise first order) at box borders and next to electrode
//!   nodes, so the stencil never reaches into an electrode,
//! * electrode nodes themselves use only their non-electrode neighbours (this
//!   gives the surface field estimate) and return 0 along an axis with none.
//!
//! "Usable" means in range and not an electrode node (box-face nodes that are
//! not masked are ordinary Dirichlet values and are usable).

use crate::Solution;
use rayon::prelude::*;

impl Solution {
    #[inline]
    fn is_electrode(&self, idx: usize) -> bool {
        !self.mask.is_empty() && self.mask[idx] != 0
    }

    /// `d phi / d(axis)` at node `idx`, whose index along the axis is `i` of `n` and
    /// whose linear-index stride along the axis is `s`.
    #[inline]
    fn axis_deriv(&self, idx: usize, s: usize, i: usize, n: usize) -> f64 {
        let phi = &self.phi;
        let h = self.grid.h;
        let has_m = i > 0 && !self.is_electrode(idx - s);
        let has_p = i + 1 < n && !self.is_electrode(idx + s);
        if has_m && has_p {
            (phi[idx + s] - phi[idx - s]) / (2.0 * h)
        } else if has_p {
            if i + 2 < n && !self.is_electrode(idx + 2 * s) {
                (-3.0 * phi[idx] + 4.0 * phi[idx + s] - phi[idx + 2 * s]) / (2.0 * h)
            } else {
                (phi[idx + s] - phi[idx]) / h
            }
        } else if has_m {
            if i >= 2 && !self.is_electrode(idx - 2 * s) {
                (3.0 * phi[idx] - 4.0 * phi[idx - s] + phi[idx - 2 * s]) / (2.0 * h)
            } else {
                (phi[idx] - phi[idx - s]) / h
            }
        } else if self.is_electrode(idx) {
            0.0
        } else {
            // Free node squeezed between electrodes (or a border node whose only
            // neighbour is an electrode): fall back to whatever neighbours exist.
            match (i > 0, i + 1 < n) {
                (true, true) => (phi[idx + s] - phi[idx - s]) / (2.0 * h),
                (false, true) => (phi[idx + s] - phi[idx]) / h,
                (true, false) => (phi[idx] - phi[idx - s]) / h,
                _ => 0.0,
            }
        }
    }

    /// `E = -grad(phi)` (V/m) at node `(i,j,k)`.
    pub fn node_field(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        let g = &self.grid;
        let idx = g.idx(i, j, k);
        [
            -self.axis_deriv(idx, 1, i, g.nx),
            -self.axis_deriv(idx, g.nx, j, g.ny),
            -self.axis_deriv(idx, g.nx * g.ny, k, g.nz),
        ]
    }

    /// `E = -grad(phi)` at every node as `[Ex, Ey, Ez, phi]` (`f32`), index
    /// `i + nx*(j + ny*k)`.  Units: V/m for `phi` in V (i.e. V/m per volt of basis).
    pub fn efield_f32(&self) -> Vec<[f32; 4]> {
        let g = self.grid;
        let mut out = vec![[0.0f32; 4]; g.len()];
        let plane = g.nx * g.ny;
        out.par_chunks_mut(plane).enumerate().for_each(|(k, chunk)| {
            for j in 0..g.ny {
                for i in 0..g.nx {
                    let e = self.node_field(i, j, k);
                    let idx = g.idx(i, j, k);
                    chunk[i + g.nx * j] = [e[0] as f32, e[1] as f32, e[2] as f32, self.phi[idx] as f32];
                }
            }
        });
        out
    }

    /// Tri-linear potential and tri-linearly interpolated nodal field at an
    /// arbitrary point.  Points outside the grid are clamped onto it.
    pub fn sample(&self, x: f64, y: f64, z: f64) -> (f64, [f64; 3]) {
        let g = &self.grid;
        let locate = |p: f64, o: f64, n: usize| {
            let f = ((p - o) / g.h).clamp(0.0, (n - 1) as f64);
            let i = (f.floor() as usize).min(n - 2);
            (i, f - i as f64)
        };
        let (i, tx) = locate(x, g.origin[0], g.nx);
        let (j, ty) = locate(y, g.origin[1], g.ny);
        let (k, tz) = locate(z, g.origin[2], g.nz);
        let mut phi = 0.0;
        let mut e = [0.0; 3];
        for dk in 0..2 {
            for dj in 0..2 {
                for di in 0..2 {
                    let w = (if di == 1 { tx } else { 1.0 - tx })
                        * (if dj == 1 { ty } else { 1.0 - ty })
                        * (if dk == 1 { tz } else { 1.0 - tz });
                    if w == 0.0 {
                        continue;
                    }
                    phi += w * self.phi[g.idx(i + di, j + dj, k + dk)];
                    let en = self.node_field(i + di, j + dj, k + dk);
                    for c in 0..3 {
                        e[c] += w * en[c];
                    }
                }
            }
        }
        (phi, e)
    }
}
