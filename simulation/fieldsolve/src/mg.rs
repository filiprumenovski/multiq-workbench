//! Geometric-algebraic multigrid for the 7-point Laplacian with masked Dirichlet nodes.
//!
//! # Why Galerkin coarse operators
//!
//! Electrodes are arbitrary node masks (curved staircases, one-node-thick
//! plates, ...).  The textbook choice "coarse node fixed iff the coincident
//! fine node is fixed, coarse operator = re-discretised Laplacian" is *not*
//! robust for such geometry: a plate that sits on odd fine planes disappears on
//! the coarse grid and the coarse correction then diverges, and curved
//! staircase surfaces give convergence factors around 0.7.  Instead we use
//!
//! * **masked tri-linear interpolation** `P'`: `P` is the tri-linear
//!   interpolation (weights 1, 1/2, 1/4, 1/8) with the rows of fixed fine nodes
//!   removed, so every coarse-grid function vanishes on all electrodes;
//! * restriction `R = P'^T`;
//! * coarse operator `A_c = P'^T A P'` (Galerkin), a 27-point stencil.
//!
//! The coarse space is a subspace of the fine functions that satisfy the
//! homogeneous Dirichlet condition, and the coarse operator is the exact
//! projection of the fine one, so the two-grid correction is energy
//! non-increasing for any mask.  A coarse node is "fixed" iff the coincident
//! fine node is fixed (its unknown is eliminated; couplings to it vanish).
//!
//! Away from electrodes and box faces all coarse stencils on a level are
//! identical, so only the *deviating* ("dirty") nodes store their own 27
//! coefficients (`sid` index into `dst`); the rest use the shared stencil
//! `cst`.  The finest level is never stored: it is the constant 7-point stencil
//! `(6, -1, ..., -1)` (the Laplacian scaled by `h^2`), solved with fast
//! red-black Gauss-Seidel.  Coarse levels use an 8-colour Gauss-Seidel (a
//! 27-point stencil couples nodes of equal red-black colour).
//!
//! Grid sizes need not be `2^k+1`; coarse node `I` coincides with fine node
//! `2I`, plus (for an even node count) the last coarse node coincides with the
//! last fine node.  An axis stops coarsening once it has fewer than 5 nodes.
//!
//! # Cycle
//!
//! Optional full-multigrid start (nested iteration on the residual equation of
//! the initial guess), then V-cycles with red-black/8-colour GS smoothing.
//! Convergence is judged on the *true* residual of the finest-level iterate.

#![allow(clippy::needless_range_loop)]

use crate::par::{SyncPtr, par_max, par_ranges};
use crate::{Problem, Solution};
use std::time::Instant;

/// Tunable parameters of the multigrid solver.
#[derive(Clone, Debug)]
pub struct Options {
    /// Pre-smoothing sweeps on every level.  With `0` the finest-level
    /// residual pass that checks convergence is reused as the first step of the
    /// cycle (post-smoothing of the previous cycle acts as pre-smoothing).
    pub nu_pre: usize,
    /// Post-smoothing sweeps on every level.
    pub nu_post: usize,
    /// Start from a nested-iteration (full multigrid) initial guess.  Off by default: for
    /// Laplace problems with electrode-dominated errors the first V-cycles are just as effective
    /// and the nested start costs about as much as it saves.
    pub fmg: bool,
    /// V-cycles used on each coarser level during the FMG start.
    pub fmg_cycles: usize,
    /// Maximum number of fine-level V-cycles before giving up.
    pub max_cycles: usize,
    /// Allow the red-black SOR fallback if multigrid stagnates or diverges.
    pub allow_sor_fallback: bool,
    /// Maximum number of SOR sweeps in the fallback.
    pub max_sor_sweeps: usize,
    /// Print a line per cycle to stderr.
    pub verbose: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            nu_pre: 0,
            nu_post: 2,
            fmg: false,
            fmg_cycles: 1,
            max_cycles: 60,
            allow_sor_fallback: true,
            max_sor_sweeps: 20000,
            verbose: false,
        }
    }
}

/// Diagnostics returned by [`solve`](crate::solve).
#[derive(Clone, Debug)]
pub struct SolveStats {
    /// Number of finest-level V-cycles (plus SOR sweeps if the fallback was used).
    pub iterations: usize,
    /// True max-norm residual of the returned solution: `max_free |sum of 6 neighbours - 6 phi|`
    /// divided by `max(|fixed|, |outer_bc|)`.
    pub final_residual: f64,
    /// Wall-clock seconds for the whole `solve` call (setup + iterations).
    pub seconds: f64,
    /// `final_residual <= tol` was reached.
    pub converged: bool,
    /// The red-black SOR fallback had to be used (multigrid misbehaved).
    pub used_sor_fallback: bool,
    /// Number of multigrid levels.
    pub levels: usize,
    /// Geometric-mean residual reduction per V-cycle.
    pub convergence_factor: f64,
    /// Seconds spent building the multigrid hierarchy (Galerkin operators).
    pub setup_seconds: f64,
}

// ---------------------------------------------------------------------------
// Transfer operators (1D, tensor-product)
// ---------------------------------------------------------------------------

/// Fine <-> coarse transfer along one axis (tri-linear interpolation).
struct Transfer {
    /// Whether this axis is coarsened.
    coarsened: bool,
    /// Coarse index -> fine index.
    cmap: Vec<usize>,
    /// Fine node -> lower/upper coarse node of the interpolation pair.
    pc0: Vec<u32>,
    pc1: Vec<u32>,
    /// Weight of `pc1` (0 for coincident nodes, 1/2 otherwise).
    pw: Vec<f64>,
    /// Restriction (`P^T`) weights of the fine nodes `cmap-1`, `cmap`, `cmap+1`.
    rwm: Vec<f64>,
    rwp: Vec<f64>,
}

impl Transfer {
    fn build(n: usize, coarsen: bool) -> (Transfer, usize) {
        if !coarsen {
            let t = Transfer {
                coarsened: false,
                cmap: (0..n).collect(),
                pc0: (0..n as u32).collect(),
                pc1: (0..n as u32).collect(),
                pw: vec![0.0; n],
                rwm: vec![0.0; n],
                rwp: vec![0.0; n],
            };
            return (t, n);
        }
        let nc = if n % 2 == 1 { n.div_ceil(2) } else { n / 2 + 1 };
        let mut cmap: Vec<usize> = (0..nc).map(|i| 2 * i).collect();
        cmap[nc - 1] = n - 1;
        let mut is_coarse = vec![false; n];
        for &f in &cmap {
            is_coarse[f] = true;
        }
        let mut pc0 = vec![0u32; n];
        let mut pc1 = vec![0u32; n];
        let mut pw = vec![0.0; n];
        for f in 0..n {
            if is_coarse[f] {
                let c = if f == n - 1 { nc - 1 } else { f / 2 } as u32;
                pc0[f] = c;
                pc1[f] = c;
            } else {
                let c0 = ((f - 1) / 2) as u32;
                pc0[f] = c0;
                pc1[f] = c0 + 1;
                pw[f] = 0.5;
            }
        }
        let mut rwm = vec![0.0; nc];
        let mut rwp = vec![0.0; nc];
        for c in 0..nc {
            let f = cmap[c];
            if f >= 1 && !is_coarse[f - 1] {
                rwm[c] = 0.5;
            }
            if f + 1 < n && !is_coarse[f + 1] {
                rwp[c] = 0.5;
            }
        }
        (Transfer { coarsened: true, cmap, pc0, pc1, pw, rwm, rwp }, nc)
    }
}

// ---------------------------------------------------------------------------
// Level
// ---------------------------------------------------------------------------

const S27: usize = 27;
const DIAG: usize = 13;

struct Level {
    /// Finest level: constant 7-point stencil, solved for the potential itself.
    is_fine: bool,
    nx: usize,
    ny: usize,
    nz: usize,
    /// 1 = Dirichlet node (electrode, outer face; on coarse levels: coincident fine node fixed).
    fix: Vec<u8>,
    /// Iterate (finest level: potential incl. Dirichlet values; coarse levels: correction, zero at fixed nodes).
    u: Vec<f64>,
    /// Right-hand side (coarse levels only; the finest level solves `A u = 0`).
    b: Vec<f64>,
    /// Residual `b - A u` (zero at fixed nodes).
    r: Vec<f64>,
    /// Shared stencil of "clean" nodes (27 coefficients, index `(dx+1) + 3 (dy+1) + 9 (dz+1)`).
    cst: [f64; S27],
    /// Per-node stencil id: 0 = shared stencil, else 1 + index into `dst`.  Empty on the finest level.
    sid: Vec<u32>,
    /// Stencils of the deviating nodes.
    dst: Vec<[f64; S27]>,
    /// Transfer operators to the next coarser level (x, y, z).
    tr: Option<[Transfer; 3]>,
    /// Finest level only, Poisson mode: right-hand side `h^2 rho / eps0` (empty = Laplace mode).
    fb: Vec<f64>,
    /// Finest level only, Poisson mode: reaction term `h^2 kappa^2 >= 0` added to the diagonal (empty = none).
    fc: Vec<f64>,
}

impl Level {
    fn nlines(&self) -> usize {
        (self.ny - 2) * (self.nz - 2)
    }

    #[inline]
    fn stencil(&self, idx: usize) -> &[f64; S27] {
        if self.sid.is_empty() {
            return &self.cst;
        }
        let id = self.sid[idx] as usize;
        if id == 0 { &self.cst } else { &self.dst[id - 1] }
    }

    fn offsets(&self) -> [isize; S27] {
        let (sy, sz) = (self.nx as isize, (self.nx * self.ny) as isize);
        let mut off = [0isize; S27];
        for t in 0..S27 {
            let (dx, dy, dz) = ((t % 3) as isize - 1, ((t / 3) % 3) as isize - 1, (t / 9) as isize - 1);
            off[t] = dx + dy * sy + dz * sz;
        }
        off
    }
}

/// Box dilation (radius 2 in each axis) of a byte mask: 1 where any node within +-2 is non-zero.
fn dilate2(src: &[u8], nx: usize, ny: usize, nz: usize) -> Vec<u8> {
    let n = nx * ny * nz;
    let mut a = vec![0u8; n];
    let mut b = vec![0u8; n];
    // x pass: src -> a
    {
        let ap = SyncPtr::new(&mut a);
        par_ranges(ny * nz, nx, |s, e| {
            for line in s..e {
                let base = line * nx;
                let row = &src[base..base + nx];
                // SAFETY: disjoint lines.
                let out = unsafe { std::slice::from_raw_parts_mut(ap.get().add(base), nx) };
                for i in 0..nx {
                    let (lo, hi) = (i.saturating_sub(2), (i + 2).min(nx - 1));
                    out[i] = row[lo..=hi].iter().fold(0u8, |m, &v| m | v);
                }
            }
        });
    }
    // y pass: a -> b
    {
        let bp = SyncPtr::new(&mut b);
        let a = &a;
        par_ranges(nz, nx * ny, |s, e| {
            for k in s..e {
                for j in 0..ny {
                    let base = (k * ny + j) * nx;
                    // SAFETY: disjoint planes.
                    let out = unsafe { std::slice::from_raw_parts_mut(bp.get().add(base), nx) };
                    let (lo, hi) = (j.saturating_sub(2), (j + 2).min(ny - 1));
                    for jj in lo..=hi {
                        let row = &a[(k * ny + jj) * nx..][..nx];
                        for (o, r) in out.iter_mut().zip(row) {
                            *o |= *r;
                        }
                    }
                }
            }
        });
    }
    // z pass: b -> a
    {
        let ap = SyncPtr::new(&mut a);
        let b = &b;
        let plane = nx * ny;
        par_ranges(nz, plane, |s, e| {
            for k in s..e {
                // SAFETY: disjoint planes.
                let out = unsafe { std::slice::from_raw_parts_mut(ap.get().add(k * plane), plane) };
                out.fill(0);
                let (lo, hi) = (k.saturating_sub(2), (k + 2).min(nz - 1));
                for kk in lo..=hi {
                    let row = &b[kk * plane..][..plane];
                    for (o, r) in out.iter_mut().zip(row) {
                        *o |= *r;
                    }
                }
            }
        });
    }
    a
}

/// Per-axis support weight of a fine node (offset `a` from the coincident node) in a coarse node.
#[inline]
fn supp_weight(coarsened: bool, a: i32) -> f64 {
    match (coarsened, a) {
        (_, 0) => 1.0,
        (true, -1) | (true, 1) => 0.5,
        _ => 0.0,
    }
}

/// Coarse nodes (offset from the coarse node `I`, weight) that a fine node at offset `d`
/// from the coincident node `2I` contributes to under interpolation.
#[inline]
fn interp_targets(coarsened: bool, d: i32) -> ([(i32, f64); 2], usize) {
    if !coarsened {
        ([(d, 1.0), (0, 0.0)], 1)
    } else if d % 2 == 0 {
        ([(d / 2, 1.0), (0, 0.0)], 1)
    } else {
        ([((d - 1).div_euclid(2), 0.5), ((d + 1).div_euclid(2), 0.5)], 2)
    }
}

/// Galerkin stencil `(P'^T A P')` row of the coarse node whose coincident fine node is `f0 = [i,j,k]`.
fn galerkin_row(lv: &Level, f0: [usize; 3], coarsened: [bool; 3]) -> [f64; S27] {
    let mut out = [0.0f64; S27];
    let (nx, ny) = (lv.nx, lv.ny);
    for fz in -1i32..=1 {
        let wz = supp_weight(coarsened[2], fz);
        if wz == 0.0 {
            continue;
        }
        for fy in -1i32..=1 {
            let wy = supp_weight(coarsened[1], fy);
            if wy == 0.0 {
                continue;
            }
            for fx in -1i32..=1 {
                let wx = supp_weight(coarsened[0], fx);
                if wx == 0.0 {
                    continue;
                }
                let wf = wx * wy * wz;
                let f = [(f0[0] as i32 + fx) as usize, (f0[1] as i32 + fy) as usize, (f0[2] as i32 + fz) as usize];
                let fidx = (f[2] * ny + f[1]) * nx + f[0];
                if lv.fix[fidx] != 0 {
                    continue; // row of a fixed node is removed by the masked prolongation
                }
                let mut rowv = *lv.stencil(fidx);
                if !lv.fc.is_empty() {
                    rowv[DIAG] += lv.fc[fidx];
                }
                let row = &rowv;
                for t in 0..S27 {
                    let a = row[t];
                    if a == 0.0 {
                        continue;
                    }
                    let (dx, dy, dz) = ((t % 3) as i32 - 1, ((t / 3) % 3) as i32 - 1, (t / 9) as i32 - 1);
                    let g = [
                        (f[0] as i32 + dx) as usize,
                        (f[1] as i32 + dy) as usize,
                        (f[2] as i32 + dz) as usize,
                    ];
                    if lv.fix[(g[2] * ny + g[1]) * nx + g[0]] != 0 {
                        continue; // fixed fine node: zero correction
                    }
                    let (tx, ntx) = interp_targets(coarsened[0], fx + dx);
                    let (ty, nty) = interp_targets(coarsened[1], fy + dy);
                    let (tz, ntz) = interp_targets(coarsened[2], fz + dz);
                    let wa = wf * a;
                    for &(oz, w3) in &tz[..ntz] {
                        for &(oy, w2) in &ty[..nty] {
                            for &(ox, w1) in &tx[..ntx] {
                                out[(ox + 1) as usize + 3 * (oy + 1) as usize + 9 * (oz + 1) as usize] +=
                                    wa * w1 * w2 * w3;
                            }
                        }
                    }
                }
            }
        }
    }
    out
}

/// Build the next coarser level (and store the transfer operators in `fine`).
fn coarsen(fine: &mut Level) -> Level {
    let (tx, ncx) = Transfer::build(fine.nx, fine.nx >= 5);
    let (ty, ncy) = Transfer::build(fine.ny, fine.ny >= 5);
    let (tz, ncz) = Transfer::build(fine.nz, fine.nz >= 5);
    let coarsened = [tx.coarsened, ty.coarsened, tz.coarsened];
    let n = ncx * ncy * ncz;
    let (nxf, nyf, nzf) = (fine.nx, fine.ny, fine.nz);

    // coarse fixed mask by injection
    let mut fix = vec![0u8; n];
    {
        let fixp = SyncPtr::new(&mut fix);
        let ffix = &fine.fix;
        par_ranges(ncy * ncz, ncx, |a, b| {
            for line in a..b {
                let (j, k) = (line % ncy, line / ncy);
                let fbase = (tz.cmap[k] * nyf + ty.cmap[j]) * nxf;
                let cbase = (k * ncy + j) * ncx;
                for i in 0..ncx {
                    // SAFETY: disjoint line.
                    unsafe { *fixp.get().add(cbase + i) = ffix[fbase + tx.cmap[i]] };
                }
            }
        });
    }

    // nodes whose Galerkin row can differ from the shared one: anything within +-2 of a
    // fixed node or of a node with a non-shared stencil
    let dirty_src: Vec<u8> = if fine.sid.is_empty() {
        if fine.fc.is_empty() {
            fine.fix.clone()
        } else {
            fine.fix.iter().zip(&fine.fc).map(|(&f, &c)| (f != 0 || c != 0.0) as u8).collect()
        }
    } else {
        fine.fix.iter().zip(&fine.sid).map(|(&f, &s)| (f != 0 || s != 0) as u8).collect()
    };
    let dil = dilate2(&dirty_src, nxf, nyf, nzf);
    drop(dirty_src);

    // count dirty free coarse nodes per line, assign ids by prefix sum
    let nl = ncy * ncz;
    let mut counts = vec![0u32; nl + 1];
    {
        let cp = SyncPtr::new(&mut counts);
        let (fix, dil) = (&fix, &dil);
        par_ranges(nl, ncx, |a, b| {
            for line in a..b {
                let (j, k) = (line % ncy, line / ncy);
                let fbase = (tz.cmap[k] * nyf + ty.cmap[j]) * nxf;
                let cbase = line * ncx;
                let mut c = 0u32;
                for i in 1..ncx - 1 {
                    if fix[cbase + i] == 0 && dil[fbase + tx.cmap[i]] != 0 {
                        c += 1;
                    }
                }
                // SAFETY: one slot per line.
                unsafe { *cp.get().add(line + 1) = c };
            }
        });
    }
    for l in 0..nl {
        counts[l + 1] += counts[l];
    }
    let ndirty = counts[nl] as usize;
    let mut sid = vec![0u32; n];
    let mut dst = vec![[0.0f64; S27]; ndirty];
    {
        let sp = SyncPtr::new(&mut sid);
        let dp = SyncPtr::new(&mut dst);
        let (fix, dil, counts) = (&fix, &dil, &counts);
        let finer: &Level = fine;
        par_ranges(nl, ncx * 8, |a, b| {
            for line in a..b {
                let (j, k) = (line % ncy, line / ncy);
                let fbase = (tz.cmap[k] * nyf + ty.cmap[j]) * nxf;
                let cbase = line * ncx;
                let mut next = counts[line] as usize;
                for i in 1..ncx - 1 {
                    if fix[cbase + i] == 0 && dil[fbase + tx.cmap[i]] != 0 {
                        let mut row = galerkin_row(finer, [tx.cmap[i], ty.cmap[j], tz.cmap[k]], coarsened);
                        // drop couplings to fixed coarse nodes (their correction is zero)
                        for t in 0..S27 {
                            let (dx, dy, dz) = ((t % 3) as isize - 1, ((t / 3) % 3) as isize - 1, (t / 9) as isize - 1);
                            let (ii, jj, kk) = (i as isize + dx, j as isize + dy, k as isize + dz);
                            if fix[(kk as usize * ncy + jj as usize) * ncx + ii as usize] != 0 {
                                row[t] = 0.0;
                            }
                        }
                        // SAFETY: each line owns the id range [counts[line], counts[line+1]).
                        unsafe {
                            *dp.get().add(next) = row;
                            *sp.get().add(cbase + i) = next as u32 + 1;
                        }
                        next += 1;
                    }
                }
            }
        });
    }

    // shared stencil from any clean free node (single-threaded search; cheap)
    let mut cst = [0.0f64; S27];
    'search: for k in 1..ncz - 1 {
        for j in 1..ncy - 1 {
            for i in 1..ncx - 1 {
                let c = (k * ncy + j) * ncx + i;
                if fix[c] == 0 && sid[c] == 0 {
                    cst = galerkin_row(fine, [tx.cmap[i], ty.cmap[j], tz.cmap[k]], coarsened);
                    break 'search;
                }
            }
        }
    }

    fine.tr = Some([tx, ty, tz]);
    Level {
        is_fine: false,
        nx: ncx,
        ny: ncy,
        nz: ncz,
        fix,
        u: vec![0.0; n],
        b: vec![0.0; n],
        r: vec![0.0; n],
        cst,
        sid,
        dst,
        tr: None,
        fb: Vec::new(),
        fc: Vec::new(),
    }
}

// ---------------------------------------------------------------------------
// Kernels: finest level (7-point, scaled Laplacian: 6 u - sum of neighbours)
// ---------------------------------------------------------------------------

/// One red-black half sweep of (over-relaxed) Gauss-Seidel on the finest level.
fn smooth7_half(lv: &mut Level, color: usize, omega: f64) {
    if !lv.fc.is_empty() {
        return smooth7c_half(lv, color, omega);
    }
    let (nx, ny) = (lv.nx, lv.ny);
    let (sy, sz) = (nx, nx * ny);
    let up = SyncPtr::new(&mut lv.u);
    let fix = &lv.fix[..];
    let sixth = 1.0 / 6.0;
    par_ranges(lv.nlines(), nx, |a, b| {
        let u = up.get();
        for line in a..b {
            let j = 1 + line % (ny - 2);
            let k = 1 + line / (ny - 2);
            let base = k * sz + j * sy;
            let mut i = if ((1 + j + k) & 1) == color { 1 } else { 2 };
            while i < nx - 1 {
                let idx = base + i;
                // SAFETY: interior node; the six neighbours have the other colour, so
                // no concurrent task writes them.
                unsafe {
                    if *fix.get_unchecked(idx) == 0 {
                        let s = *u.add(idx - 1)
                            + *u.add(idx + 1)
                            + *u.add(idx - sy)
                            + *u.add(idx + sy)
                            + *u.add(idx - sz)
                            + *u.add(idx + sz);
                        let old = *u.add(idx);
                        *u.add(idx) = old + omega * (s * sixth - old);
                    }
                }
                i += 2;
            }
        }
    });
}

/// Residual `r = -(6u - sum of neighbours)` on free nodes of the finest level; returns max |r|.
fn residual7(lv: &mut Level, store: bool) -> f64 {
    if !lv.fc.is_empty() {
        return residual7c(lv, store);
    }
    let (nx, ny) = (lv.nx, lv.ny);
    let (sy, sz) = (nx, nx * ny);
    let rp = SyncPtr::new(&mut lv.r);
    let fix = &lv.fix[..];
    let u = &lv.u[..];
    par_max(lv.nlines(), nx, |a, b| {
        let mut local = 0.0f64;
        for line in a..b {
            let j = 1 + line % (ny - 2);
            let k = 1 + line / (ny - 2);
            let base = k * sz + j * sy;
            let c = &u[base..base + nx];
            let ym = &u[base - sy..base - sy + nx];
            let yp = &u[base + sy..base + sy + nx];
            let zm = &u[base - sz..base - sz + nx];
            let zp = &u[base + sz..base + sz + nx];
            let fx = &fix[base..base + nx];
            if store {
                // SAFETY: lines are disjoint across tasks.
                let rl = unsafe { std::slice::from_raw_parts_mut(rp.get().add(base), nx) };
                for i in 1..nx - 1 {
                    let s = c[i - 1] + c[i + 1] + ym[i] + yp[i] + zm[i] + zp[i];
                    let r = if fx[i] != 0 { 0.0 } else { s - 6.0 * c[i] };
                    local = local.max(r.abs());
                    rl[i] = r;
                }
            } else {
                for i in 1..nx - 1 {
                    let s = c[i - 1] + c[i + 1] + ym[i] + yp[i] + zm[i] + zp[i];
                    let r = if fx[i] != 0 { 0.0 } else { s - 6.0 * c[i] };
                    local = local.max(r.abs());
                }
            }
        }
        local
    })
}


/// Poisson-mode red-black half sweep on the finest level: `(6 + c_i) u_i - sum nb = b_i`.
fn smooth7c_half(lv: &mut Level, color: usize, omega: f64) {
    let (nx, ny) = (lv.nx, lv.ny);
    let (sy, sz) = (nx, nx * ny);
    let up = SyncPtr::new(&mut lv.u);
    let fix = &lv.fix[..];
    let (fb, fc) = (&lv.fb[..], &lv.fc[..]);
    par_ranges(lv.nlines(), nx, |a, b| {
        let u = up.get();
        for line in a..b {
            let j = 1 + line % (ny - 2);
            let k = 1 + line / (ny - 2);
            let base = k * sz + j * sy;
            let mut i = if ((1 + j + k) & 1) == color { 1 } else { 2 };
            while i < nx - 1 {
                let idx = base + i;
                // SAFETY: as in `smooth7_half`.
                unsafe {
                    if *fix.get_unchecked(idx) == 0 {
                        let s = *u.add(idx - 1)
                            + *u.add(idx + 1)
                            + *u.add(idx - sy)
                            + *u.add(idx + sy)
                            + *u.add(idx - sz)
                            + *u.add(idx + sz);
                        let old = *u.add(idx);
                        let new = (s + *fb.get_unchecked(idx)) / (6.0 + *fc.get_unchecked(idx));
                        *u.add(idx) = old + omega * (new - old);
                    }
                }
                i += 2;
            }
        }
    });
}

/// Poisson-mode residual `r = b + sum nb - (6 + c) u` on free nodes of the finest level; returns max |r|.
fn residual7c(lv: &mut Level, store: bool) -> f64 {
    let (nx, ny) = (lv.nx, lv.ny);
    let (sy, sz) = (nx, nx * ny);
    let rp = SyncPtr::new(&mut lv.r);
    let fix = &lv.fix[..];
    let u = &lv.u[..];
    let (fb, fc) = (&lv.fb[..], &lv.fc[..]);
    par_max(lv.nlines(), nx, |a, b| {
        let mut local = 0.0f64;
        for line in a..b {
            let j = 1 + line % (ny - 2);
            let k = 1 + line / (ny - 2);
            let base = k * sz + j * sy;
            for i in 1..nx - 1 {
                let idx = base + i;
                let s = u[idx - 1] + u[idx + 1] + u[idx - sy] + u[idx + sy] + u[idx - sz] + u[idx + sz];
                let r = if fix[idx] != 0 { 0.0 } else { fb[idx] + s - (6.0 + fc[idx]) * u[idx] };
                local = local.max(r.abs());
                if store {
                    // SAFETY: lines are disjoint across tasks.
                    unsafe { *rp.get().add(idx) = r };
                }
            }
        }
        local
    })
}

// ---------------------------------------------------------------------------
// Kernels: coarse levels (27-point stencils)
// ---------------------------------------------------------------------------

/// One colour class (parity of i, j, k) of an 8-colour Gauss-Seidel sweep.
fn smooth27_color(lv: &mut Level, cx: usize, cy: usize, cz: usize) {
    let (nx, ny, nz) = (lv.nx, lv.ny, lv.nz);
    let (i0, j0, k0) = (if cx == 1 { 1 } else { 2 }, if cy == 1 { 1 } else { 2 }, if cz == 1 { 1 } else { 2 });
    if i0 > nx - 2 || j0 > ny - 2 || k0 > nz - 2 {
        return;
    }
    let nj = (ny - 2 - j0) / 2 + 1;
    let nk = (nz - 2 - k0) / 2 + 1;
    let off = lv.offsets();
    let up = SyncPtr::new(&mut lv.u);
    let (fix, b, sid, dst, cst) = (&lv.fix[..], &lv.b[..], &lv.sid[..], &lv.dst[..], &lv.cst);
    let cinv = 1.0 / cst[DIAG];
    par_ranges(nj * nk, nx * 4, |s, e| {
        let u = up.get();
        for line in s..e {
            let j = j0 + 2 * (line % nj);
            let k = k0 + 2 * (line / nj);
            let base = (k * ny + j) * nx;
            let mut i = i0;
            while i <= nx - 2 {
                let idx = base + i;
                // SAFETY: interior node; all 26 neighbours belong to other colour classes, which
                // no concurrent task writes during this half sweep.
                unsafe {
                    if *fix.get_unchecked(idx) == 0 {
                        let id = *sid.get_unchecked(idx) as usize;
                        let (coef, dinv) = if id == 0 {
                            (cst, cinv)
                        } else {
                            let c = dst.get_unchecked(id - 1);
                            (c, 1.0 / c[DIAG])
                        };
                        let mut acc = *b.get_unchecked(idx);
                        let p = u.add(idx);
                        for t in 0..S27 {
                            if t != DIAG {
                                acc -= coef[t] * *p.offset(off[t]);
                            }
                        }
                        *u.add(idx) = acc * dinv;
                    }
                }
                i += 2;
            }
        }
    });
}

/// Residual `r = b - A u` (zero at fixed nodes) on a coarse level.
fn residual27(lv: &mut Level) {
    let (nx, ny) = (lv.nx, lv.ny);
    let off = lv.offsets();
    let rp = SyncPtr::new(&mut lv.r);
    let (fix, u, b, sid, dst, cst) = (&lv.fix[..], &lv.u[..], &lv.b[..], &lv.sid[..], &lv.dst[..], &lv.cst);
    par_ranges(lv.nlines(), nx * 4, |s, e| {
        for line in s..e {
            let j = 1 + line % (ny - 2);
            let k = 1 + line / (ny - 2);
            let base = (k * ny + j) * nx;
            for i in 1..nx - 1 {
                let idx = base + i;
                let r = if fix[idx] != 0 {
                    0.0
                } else {
                    let coef = if sid[idx] == 0 { cst } else { &dst[sid[idx] as usize - 1] };
                    let mut acc = b[idx];
                    for t in 0..S27 {
                        // SAFETY: interior node, offsets stay inside the array.
                        acc -= coef[t] * unsafe { *u.get_unchecked((idx as isize + off[t]) as usize) };
                    }
                    acc
                };
                // SAFETY: disjoint line per task.
                unsafe { *rp.get().add(idx) = r };
            }
        }
    });
}

/// `n` Gauss-Seidel sweeps on any level (red-black on the finest, 8-colour on coarse levels).
fn smooth(lv: &mut Level, sweeps: usize, omega: f64) {
    for _ in 0..sweeps {
        if lv.is_fine {
            smooth7_half(lv, 0, omega);
            smooth7_half(lv, 1, omega);
        } else {
            for c in 0..8 {
                smooth27_color(lv, c & 1, (c >> 1) & 1, (c >> 2) & 1);
            }
        }
    }
}

/// Residual of a level (finest: `-A u`; coarse: `b - A u`), stored in `lv.r`.
fn residual_store(lv: &mut Level) {
    if lv.is_fine {
        residual7(lv, true);
    } else {
        residual27(lv);
    }
}

// ---------------------------------------------------------------------------
// Transfers
// ---------------------------------------------------------------------------

/// `coarse_dst = P'^T src` (zero on coarse fixed nodes).
fn restrict(src: &[f64], fine: &Level, coarse_dst: &mut [f64], coarse: (usize, usize, usize), cfix: &[u8]) {
    let [tx, ty, tz] = fine.tr.as_ref().expect("transfer operators");
    let (nxf, nyf) = (fine.nx, fine.ny);
    let (ncx, ncy, ncz) = coarse;
    let dp = SyncPtr::new(coarse_dst);
    par_ranges((ncy - 2) * (ncz - 2), nxf * 4, |a, b| {
        let mut t = vec![0.0f64; nxf];
        for line in a..b {
            let jj = 1 + line % (ncy - 2);
            let kk = 1 + line / (ncy - 2);
            let (fj, fk) = (ty.cmap[jj], tz.cmap[kk]);
            t.fill(0.0);
            let wzs = [(fk - 1, tz.rwm[kk]), (fk, 1.0), (fk + 1, tz.rwp[kk])];
            let wys = [(fj - 1, ty.rwm[jj]), (fj, 1.0), (fj + 1, ty.rwp[jj])];
            for &(kf, wz) in &wzs {
                if wz == 0.0 {
                    continue;
                }
                for &(jf, wy) in &wys {
                    if wy == 0.0 {
                        continue;
                    }
                    let w = wz * wy;
                    let row = &src[(kf * nyf + jf) * nxf..][..nxf];
                    for (ti, ri) in t.iter_mut().zip(row) {
                        *ti += w * ri;
                    }
                }
            }
            let cbase = (kk * ncy + jj) * ncx;
            for ii in 1..ncx - 1 {
                let f = tx.cmap[ii];
                let v = t[f] + tx.rwm[ii] * t[f - 1] + tx.rwp[ii] * t[f + 1];
                // SAFETY: disjoint line per (jj,kk).
                unsafe {
                    *dp.get().add(cbase + ii) = if cfix[cbase + ii] != 0 { 0.0 } else { v };
                }
            }
        }
    });
}

/// `fine.u[free] += P coarse.u` (fixed fine nodes are left alone).
fn prolong_add(coarse: &Level, fine: &mut Level) {
    let [tx, ty, tz] = fine.tr.as_ref().expect("transfer operators");
    let (nx, ny, nz) = (fine.nx, fine.ny, fine.nz);
    let (ncx, ncy) = (coarse.nx, coarse.ny);
    let up = SyncPtr::new(&mut fine.u);
    let fix = &fine.fix[..];
    let cu = &coarse.u[..];
    par_ranges((ny - 2) * (nz - 2), nx * 2, |a, b| {
        let mut t = vec![0.0f64; ncx];
        for line in a..b {
            let j = 1 + line % (ny - 2);
            let k = 1 + line / (ny - 2);
            let (j0, j1, wy) = (ty.pc0[j] as usize, ty.pc1[j] as usize, ty.pw[j]);
            let (k0, k1, wz) = (tz.pc0[k] as usize, tz.pc1[k] as usize, tz.pw[k]);
            let row = |kk: usize, jj: usize| &cu[(kk * ncy + jj) * ncx..][..ncx];
            let (r00, r10, r01, r11) = (row(k0, j0), row(k0, j1), row(k1, j0), row(k1, j1));
            let (w00, w10, w01, w11) = ((1.0 - wy) * (1.0 - wz), wy * (1.0 - wz), (1.0 - wy) * wz, wy * wz);
            for i in 0..ncx {
                t[i] = w00 * r00[i] + w10 * r10[i] + w01 * r01[i] + w11 * r11[i];
            }
            let base = (k * ny + j) * nx;
            for i in 1..nx - 1 {
                let idx = base + i;
                if fix[idx] == 0 {
                    let w = tx.pw[i];
                    let v = (1.0 - w) * t[tx.pc0[i] as usize] + w * t[tx.pc1[i] as usize];
                    // SAFETY: disjoint line per task.
                    unsafe {
                        *up.get().add(idx) += v;
                    }
                }
            }
        }
    });
}

fn fill_zero(v: &mut [f64]) {
    const BLK: usize = 4096;
    let n = v.len();
    let p = SyncPtr::new(v);
    par_ranges(n.div_ceil(BLK), BLK, |a, b| {
        let lo = a * BLK;
        let hi = (b * BLK).min(n);
        // SAFETY: disjoint ranges.
        unsafe { std::slice::from_raw_parts_mut(p.get().add(lo), hi - lo) }.fill(0.0);
    });
}

// ---------------------------------------------------------------------------
// Cycles
// ---------------------------------------------------------------------------

/// Coarsest level (every axis <= 4 nodes, <= 8 unknowns): converge with plain GS.
fn coarsest_solve(lv: &mut Level) {
    smooth(lv, 80, 1.0);
}

/// One V-cycle on `levels[0]`.
///
/// * `zero_start`: `levels[0].u == 0` on entry (so the residual is just `b`);
/// * `residual_ready`: `levels[0].r` already holds the current residual (finest level, `nu_pre == 0`).
fn vcycle(levels: &mut [Level], zero_start: bool, residual_ready: bool, o: &Options) {
    let (cur, rest) = levels.split_first_mut().expect("levels");
    if rest.is_empty() {
        coarsest_solve(cur);
        return;
    }
    if o.nu_pre > 0 {
        smooth(cur, o.nu_pre, 1.0);
    }
    let coarse_dims = (rest[0].nx, rest[0].ny, rest[0].nz);
    let c0 = &mut rest[0];
    if o.nu_pre == 0 && zero_start && !cur.is_fine {
        restrict(&cur.b, cur, &mut c0.b, coarse_dims, &c0.fix); // u = 0  =>  r = b
    } else {
        if !(o.nu_pre == 0 && residual_ready) {
            residual_store(cur);
        }
        restrict(&cur.r, cur, &mut c0.b, coarse_dims, &c0.fix);
    }
    fill_zero(&mut c0.u);
    vcycle(rest, true, false, o);
    prolong_add(&rest[0], cur);
    smooth(cur, o.nu_post, 1.0);
}

fn build_level0(p: &Problem, outer_bc: f64) -> (Level, f64) {
    let g = &p.grid;
    let (nx, ny, nz) = (g.nx, g.ny, g.nz);
    let n = nx * ny * nz;
    let mut fix = vec![0u8; n];
    let mut u = vec![0.0f64; n];
    let maxfix;
    {
        let fixp = SyncPtr::new(&mut fix);
        let up = SyncPtr::new(&mut u);
        let mask = &p.mask;
        let fixed = &p.fixed;
        let plane = nx * ny;
        maxfix = par_max(nz, plane, |a, b| {
            let mut m = 0.0f64;
            for k in a..b {
                for j in 0..ny {
                    let face_jk = k == 0 || k == nz - 1 || j == 0 || j == ny - 1;
                    let base = (k * ny + j) * nx;
                    for i in 0..nx {
                        let idx = base + i;
                        let face = face_jk || i == 0 || i == nx - 1;
                        let (fx, val) = if mask[idx] != 0 {
                            m = m.max(fixed[idx].abs());
                            (1u8, fixed[idx])
                        } else if face {
                            (1u8, outer_bc)
                        } else {
                            (0u8, outer_bc)
                        };
                        // SAFETY: disjoint planes per task.
                        unsafe {
                            *fixp.get().add(idx) = fx;
                            *up.get().add(idx) = val;
                        }
                    }
                }
            }
            m
        });
    }
    let scale = {
        let s = maxfix.max(outer_bc.abs());
        if s > 0.0 { s } else { 1.0 }
    };
    let mut cst = [0.0f64; S27];
    cst[DIAG] = 6.0;
    for t in [4usize, 10, 12, 14, 16, 22] {
        cst[t] = -1.0;
    }
    let lv = Level {
        is_fine: true,
        nx,
        ny,
        nz,
        fix,
        u,
        b: Vec::new(),
        r: vec![0.0; n],
        cst,
        sid: Vec::new(),
        dst: Vec::new(),
        tr: None,
        fb: Vec::new(),
        fc: Vec::new(),
    };
    (lv, scale)
}

/// Solve with explicit [`Options`].  See [`crate::solve`].
pub fn solve_with(p: &Problem, outer_bc: f64, tol: f64, o: &Options) -> (Solution, SolveStats) {
    let t0 = Instant::now();
    let g = p.grid;
    assert!(g.nx >= 3 && g.ny >= 3 && g.nz >= 3, "grid must have at least 3 nodes per axis");
    let n = g.nx * g.ny * g.nz;
    assert_eq!(p.mask.len(), n, "mask length");
    assert_eq!(p.fixed.len(), n, "fixed length");

    let (lv0, scale) = build_level0(p, outer_bc);
    let mut levels = vec![lv0];
    while {
        let l = levels.last().unwrap();
        l.nx >= 5 || l.ny >= 5 || l.nz >= 5
    } {
        let c = coarsen(levels.last_mut().unwrap());
        levels.push(c);
    }
    let nlev = levels.len();
    let setup_seconds = t0.elapsed().as_secs_f64();
    if o.verbose {
        for (l, lv) in levels.iter().enumerate() {
            eprintln!("  level {l}: {}x{}x{}  dirty stencils {}", lv.nx, lv.ny, lv.nz, lv.dst.len());
        }
        eprintln!("  setup {setup_seconds:.3}s");
    }
    let to_scaled = |r: f64| r / scale;

    // --- full-multigrid start (nested iteration on the residual equation) ---------------
    if o.fmg && nlev > 1 {
        residual7(&mut levels[0], true);
        {
            let (a, b) = levels.split_at_mut(1);
            let c = &mut b[0];
            let dims = (c.nx, c.ny, c.nz);
            restrict(&a[0].r, &a[0], &mut c.b, dims, &c.fix);
        }
        for l in 1..nlev - 1 {
            let (a, b) = levels.split_at_mut(l + 1);
            let c = &mut b[0];
            let dims = (c.nx, c.ny, c.nz);
            restrict(&a[l].b, &a[l], &mut c.b, dims, &c.fix);
        }
        fill_zero(&mut levels[nlev - 1].u);
        coarsest_solve(&mut levels[nlev - 1]);
        for l in (0..nlev - 1).rev() {
            {
                let (a, b) = levels.split_at_mut(l + 1);
                prolong_add(&b[0], &mut a[l]);
            }
            if l > 0 {
                for _ in 0..o.fmg_cycles {
                    vcycle(&mut levels[l..], false, false, o);
                }
            }
        }
    }

    // --- fine-level V-cycles ---------------------------------------------------
    let mut cycles = 0usize;
    let mut hist: Vec<f64> = Vec::new();
    let mut converged = false;
    let mut need_fallback = false;
    let mut final_res;
    let reuse = o.nu_pre == 0;
    loop {
        final_res = to_scaled(residual7(&mut levels[0], reuse));
        hist.push(final_res);
        if o.verbose {
            eprintln!("cycle {cycles:3}  residual {final_res:.3e}  t={:.3}s", t0.elapsed().as_secs_f64());
        }
        if final_res <= tol {
            converged = true;
            break;
        }
        if !final_res.is_finite() || cycles >= o.max_cycles {
            need_fallback = true;
            break;
        }
        // stagnation: less than a factor 2 reduction over 5 cycles (and not at the round-off floor)
        if cycles >= 5 && final_res > 1e-12 && final_res > 0.5 * hist[cycles - 5] {
            need_fallback = true;
            break;
        }
        vcycle(&mut levels, false, reuse, o);
        cycles += 1;
    }
    let mg_factor = if cycles > 0 && hist[0] > 0.0 && final_res > 0.0 {
        (final_res / hist[0]).powf(1.0 / cycles as f64)
    } else {
        0.0
    };

    // --- SOR fallback -----------------------------------------------------------
    let mut used_sor = false;
    let mut iterations = cycles;
    // (not when the residual sits at the round-off floor: the requested tolerance is simply too tight)
    if need_fallback && o.allow_sor_fallback && final_res > 1e-11 {
        used_sor = true;
        let lv = &mut levels[0];
        let pi = std::f64::consts::PI;
        let rho = ((pi / (lv.nx - 1) as f64).cos() + (pi / (lv.ny - 1) as f64).cos() + (pi / (lv.nz - 1) as f64).cos()) / 3.0;
        let omega = 2.0 / (1.0 + (1.0 - rho * rho).max(0.0).sqrt());
        if o.verbose {
            eprintln!("multigrid failed to converge: falling back to red-black SOR, omega = {omega:.5}");
        }
        let mut sweeps = 0usize;
        while sweeps < o.max_sor_sweeps {
            smooth(lv, 10, omega);
            sweeps += 10;
            final_res = to_scaled(residual7(lv, false));
            if final_res <= tol {
                converged = true;
                break;
            }
        }
        iterations += sweeps;
    }

    let lv0 = levels.swap_remove(0);
    let seconds = t0.elapsed().as_secs_f64();
    let sol = Solution { grid: g, phi: lv0.u, mask: p.mask.clone() };
    let stats = SolveStats {
        iterations,
        final_residual: final_res,
        seconds,
        converged,
        used_sor_fallback: used_sor,
        levels: nlev,
        convergence_factor: mg_factor,
        setup_seconds,
    };
    (sol, stats)
}


// ---------------------------------------------------------------------------
// Poisson / screened-Poisson solver (reusable hierarchy)
// ---------------------------------------------------------------------------

/// Solver for `-laplace(phi) + kappa2 phi = rho/eps0` on a masked Dirichlet grid, with the same Galerkin
/// multigrid hierarchy as [`solve_with`].  The hierarchy is built once per `(mask, kappa2)` and can be reused
/// for several right-hand sides (Newton steps, Woodbury solves).
///
/// Units: callers pass `react[i] = h^2 kappa2_i` (dimensionless, >= 0) and `rhs[i] = h^2 rho_i / eps0` (volts).
/// Dirichlet values (electrodes and box faces, `outer_bc`) are those of the [`Problem`].
pub struct Poisson {
    levels: Vec<Level>,
    pub grid: crate::Grid,
    mask: Vec<u8>,
    pub setup_seconds: f64,
}

impl Poisson {
    pub fn new(p: &Problem, outer_bc: f64, react: Option<&[f64]>) -> Poisson {
        let t0 = Instant::now();
        let n = p.grid.len();
        let (mut lv0, _) = build_level0(p, outer_bc);
        lv0.fb = vec![0.0; n];
        lv0.fc = match react {
            Some(c) => {
                assert_eq!(c.len(), n);
                c.to_vec()
            }
            None => vec![0.0; n],
        };
        let mut levels = vec![lv0];
        while {
            let l = levels.last().unwrap();
            l.nx >= 5 || l.ny >= 5 || l.nz >= 5
        } {
            let c = coarsen(levels.last_mut().unwrap());
            levels.push(c);
        }
        Poisson { levels, grid: p.grid, mask: p.mask.clone(), setup_seconds: t0.elapsed().as_secs_f64() }
    }

    /// Number of multigrid levels.
    pub fn n_levels(&self) -> usize {
        self.levels.len()
    }

    /// Solve `A phi = rhs` to a max-norm residual `tol_abs` (volts, i.e. `h^2 |laplace phi + rho/eps0| <= tol_abs`).
    /// `x0` (potential at all nodes; Dirichlet entries are ignored) is the initial guess; `None` = the previous
    /// solution of this object (or `outer_bc` on first use).
    pub fn solve(&mut self, rhs: &[f64], x0: Option<&[f64]>, tol_abs: f64, max_cycles: usize) -> (&[f64], SolveStats) {
        let t0 = Instant::now();
        let o = Options::default();
        {
            let lv = &mut self.levels[0];
            assert_eq!(rhs.len(), lv.fb.len());
            lv.fb.copy_from_slice(rhs);
            if let Some(x) = x0 {
                for i in 0..lv.u.len() {
                    if lv.fix[i] == 0 {
                        lv.u[i] = x[i];
                    }
                }
            }
        }
        let mut cycles = 0usize;
        let mut hist: Vec<f64> = Vec::new();
        let mut converged = false;
        let mut final_res;
        loop {
            final_res = residual7(&mut self.levels[0], true);
            hist.push(final_res);
            if final_res <= tol_abs {
                converged = true;
                break;
            }
            if !final_res.is_finite() || cycles >= max_cycles {
                break;
            }
            vcycle(&mut self.levels, false, true, &o);
            cycles += 1;
        }
        let mg_factor = if cycles > 0 && hist[0] > 0.0 && final_res > 0.0 {
            (final_res / hist[0]).powf(1.0 / cycles as f64)
        } else {
            0.0
        };
        let stats = SolveStats {
            iterations: cycles,
            final_residual: final_res,
            seconds: t0.elapsed().as_secs_f64(),
            converged,
            used_sor_fallback: false,
            levels: self.levels.len(),
            convergence_factor: mg_factor,
            setup_seconds: self.setup_seconds,
        };
        (&self.levels[0].u, stats)
    }

    /// Electrode/box-face mask in the sense of the Problem (1 = electrode node, from the Problem's mask).
    pub fn mask(&self) -> &[u8] {
        &self.mask
    }
}
