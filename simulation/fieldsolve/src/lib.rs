//! # fieldsolve
//!
//! Fast 3D electrostatic (Laplace) field solver for ion-optics simulation.
//!
//! The potential `phi` obeys `laplace(phi) = 0` on a regular Cartesian grid.
//! Electrodes are Dirichlet nodes (a `mask` plus a potential per node); the
//! outer faces of the box are Dirichlet nodes at a constant potential unless
//! they are masked themselves.  The 7-point finite-difference system is solved
//! by geometric multigrid ([`mg`] module) to a *true-residual* tolerance.
//!
//! ```
//! use fieldsolve::*;
//! // a 1 V cube electrode in the middle of a grounded box
//! let grid = Grid { nx: 33, ny: 33, nz: 33, h: 1e-3, origin: [0.0; 3] };
//! let p = Problem::from_fn(grid, |x, y, z| {
//!     let c = 16e-3;
//!     if (x - c).abs() < 2.1e-3 && (y - c).abs() < 2.1e-3 && (z - c).abs() < 2.1e-3 {
//!         Cell::Fixed(1.0)
//!     } else {
//!         Cell::Free
//!     }
//! });
//! let (sol, stats) = solve(&p, 0.0, 1e-9);
//! assert!(stats.converged);
//! let (phi, e) = sol.sample(16e-3, 16e-3, 20e-3);
//! assert!(phi > 0.0 && phi < 1.0);
//! assert!(e[2] > 0.0); // field points away from the electrode
//! ```
//!
//! # Accuracy notes
//!
//! * The scheme is second order for the interior operator, but electrodes are node masks, so a
//!   curved electrode is a *staircase*: the effective Dirichlet surface sits on the outermost fixed
//!   nodes, on average about `h/2` inside the electrode.  This gives an `O(h)` shift of curved
//!   electrode surfaces (see `VERIFICATION.md`).  Marking nodes within `R + 0.4 h` of a convex
//!   electrode of radius `R` (and `R - 0.4 h` for a concave one) as fixed removes most of it.
//! * `Cell::Fixed` may vary from node to node, so Dirichlet data can follow a known analytic
//!   solution exactly (used by the verification tests to show second-order convergence).
//!
//! # Memory and speed
//!
//! At 400^3 the whole pipeline (problem, solver, solution) peaks at about 2.2 GB; a four-rod
//! quadrupole solves to 1e-8 in about 0.65 s at 256^3 and 2.1 s at 384^3 on an Apple M4.
//!
//! See `tests/` and `VERIFICATION.md` for the accuracy and performance checks.

pub mod field;
pub mod mg;
mod par;

use rayon::prelude::*;

pub use mg::{Options, Poisson, SolveStats, solve_with};

/// Regular grid: node `(i,j,k)` sits at `origin + h*(i,j,k)` (metres).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Grid {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub h: f64,
    pub origin: [f64; 3],
}

impl Grid {
    /// Total number of nodes.
    #[inline]
    pub fn len(&self) -> usize {
        self.nx * self.ny * self.nz
    }
    /// True if the grid has no nodes.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Linear index `i + nx*(j + ny*k)`.
    #[inline]
    pub fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.nx * (j + self.ny * k)
    }
    /// Position of node `(i,j,k)` in metres.
    #[inline]
    pub fn pos(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        [
            self.origin[0] + self.h * i as f64,
            self.origin[1] + self.h * j as f64,
            self.origin[2] + self.h * k as f64,
        ]
    }
}

/// A node is either free (solved for) or a Dirichlet electrode node.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Cell {
    Free,
    Fixed(f64),
}

/// Boundary-value problem: grid plus electrode mask and potentials.
///
/// `mask[idx] == 1` marks electrode nodes with potential `fixed[idx]`;
/// `idx = i + nx*(j + ny*k)`.  Entries of `fixed` at unmasked nodes are ignored.
#[derive(Clone, Debug)]
pub struct Problem {
    pub grid: Grid,
    pub mask: Vec<u8>,
    pub fixed: Vec<f64>,
}

impl Problem {
    /// Build a problem by evaluating `f(x, y, z)` (metres) at every node, in parallel.
    pub fn from_fn(grid: Grid, f: impl Fn(f64, f64, f64) -> Cell + Sync) -> Problem {
        let n = grid.len();
        let mut mask = vec![0u8; n];
        let mut fixed = vec![0.0f64; n];
        let plane = grid.nx * grid.ny;
        if plane > 0 {
            mask.par_chunks_mut(plane)
                .zip(fixed.par_chunks_mut(plane))
                .enumerate()
                .for_each(|(k, (m, fx))| {
                    for j in 0..grid.ny {
                        for i in 0..grid.nx {
                            let [x, y, z] = grid.pos(i, j, k);
                            if let Cell::Fixed(v) = f(x, y, z) {
                                m[i + grid.nx * j] = 1;
                                fx[i + grid.nx * j] = v;
                            }
                        }
                    }
                });
        }
        Problem { grid, mask, fixed }
    }
}

/// Potential on the grid together with the electrode mask used to produce it
/// (needed by the field evaluation to avoid differencing across electrodes).
#[derive(Clone, Debug)]
pub struct Solution {
    pub grid: Grid,
    pub phi: Vec<f64>,
    /// Electrode mask copied from the [`Problem`] (`mask[idx] == 1` for electrode nodes).
    /// May be empty, in which case no node is treated as an electrode.
    pub mask: Vec<u8>,
}

impl Solution {
    /// Wrap an existing potential array (no electrode information).
    pub fn new(grid: Grid, phi: Vec<f64>) -> Solution {
        assert_eq!(phi.len(), grid.len());
        Solution { grid, phi, mask: Vec::new() }
    }
}

/// Solve the Laplace problem.
///
/// * Outer box faces are Dirichlet nodes at `outer_bc` volts unless masked.
/// * `tol` is the stopping tolerance on the **true** max-norm residual of the
///   finest-level iterate, defined as `max_free |h^2 * laplace_h(phi)| / max(|fixed|, |outer_bc|)`
///   (dimensionless; `h^2 * laplace_h(phi) = sum of 6 neighbours - 6 phi`, in volts).
///
/// Uses [`Options::default`]; see [`solve_with`] to tune the cycle.
pub fn solve(p: &Problem, outer_bc: f64, tol: f64) -> (Solution, SolveStats) {
    solve_with(p, outer_bc, tol, &Options::default())
}
