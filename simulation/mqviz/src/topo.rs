//! Topographic helpers: effective potential on the x-z slice, marching-squares contours (chained polylines),
//! index-contour label placement, mask outlines and merged row-run rectangles for solid electrodes.
use crate::world::*;
use std::collections::HashMap;

/// Effective potential [V per unit charge] on the y = axis plane: DC potential + RF pseudopotential
/// Psi = (z e / m) V_rf^2 |E_rf|^2 / (4 Omega^2).
pub fn eff_potential(viz: &VizTables, sp: &Sp, v_rf: f64, dc: [f64; 4]) -> Vec<f32> {
    let m = sp.mass * AMU;
    let k = sp.z * QE / (4.0 * m * OMEGA * OMEGA) * v_rf * v_rf;
    let n = viz.nx * viz.nz;
    (0..n).map(|i| {
        let mut u = k * viz.e2[i] as f64;
        for b in 0..4 { u += dc[b] * viz.phi_dc[b][i] as f64; }
        u as f32
    }).collect()
}

pub struct Level { pub n: i32, pub lines: Vec<Vec<[f32; 2]>> }

#[inline]
fn seg_table(case: u8, center_hi: bool) -> &'static [(u8, u8)] {
    match case {
        1 | 14 => &[(3, 0)],
        2 | 13 => &[(0, 1)],
        3 | 12 => &[(3, 1)],
        4 | 11 => &[(1, 2)],
        6 | 9 => &[(0, 2)],
        7 | 8 => &[(3, 2)],
        5 => if center_hi { &[(0, 1), (2, 3)] } else { &[(3, 0), (1, 2)] },
        10 => if center_hi { &[(3, 0), (1, 2)] } else { &[(0, 1), (2, 3)] },
        _ => &[],
    }
}

/// Marching squares for all levels n*step with n in [n_lo, n_hi] over a nx x nz grid (index i + nx*k).
/// Cells touching an invalid node (conductor) are skipped. Coordinates in grid units (i, k).
pub fn isolines(v: &[f32], valid: &[bool], nx: usize, nz: usize, step: f32, n_lo: i32, n_hi: i32) -> Vec<Level> {
    if n_hi < n_lo { return vec![]; }
    let nl = (n_hi - n_lo + 1) as usize;
    // per level segments: (edge a, edge b, p a, p b)
    let mut segs: Vec<Vec<(u32, u32, [f32; 2], [f32; 2])>> = vec![vec![]; nl];
    for k in 0..nz - 1 {
        for i in 0..nx - 1 {
            let (i00, i10, i01, i11) = (i + nx * k, i + 1 + nx * k, i + nx * (k + 1), i + 1 + nx * (k + 1));
            if !(valid[i00] && valid[i10] && valid[i01] && valid[i11]) { continue; }
            let c = [v[i00], v[i10], v[i11], v[i01]];
            let lo = c.iter().cloned().fold(f32::INFINITY, f32::min);
            let hi = c.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
            let a = ((lo / step).ceil() as i32).max(n_lo);
            let b = ((hi / step).floor() as i32).min(n_hi);
            for n in a..=b {
                let l = n as f32 * step;
                let case = (c[0] >= l) as u8 | ((c[1] >= l) as u8) << 1 | ((c[2] >= l) as u8) << 2 | ((c[3] >= l) as u8) << 3;
                if case == 0 || case == 15 { continue; }
                let center_hi = (c[0] + c[1] + c[2] + c[3]) * 0.25 >= l;
                let (fi, fk) = (i as f32, k as f32);
                let pt = |e: u8| -> (u32, [f32; 2]) {
                    match e {
                        0 => { let t = (l - c[0]) / (c[1] - c[0]); ((2 * (k * nx + i)) as u32, [fi + t, fk]) }
                        1 => { let t = (l - c[1]) / (c[2] - c[1]); ((2 * (k * nx + i + 1) + 1) as u32, [fi + 1.0, fk + t]) }
                        2 => { let t = (l - c[3]) / (c[2] - c[3]); ((2 * ((k + 1) * nx + i)) as u32, [fi + t, fk + 1.0]) }
                        _ => { let t = (l - c[0]) / (c[3] - c[0]); ((2 * (k * nx + i) + 1) as u32, [fi, fk + t]) }
                    }
                };
                for &(e0, e1) in seg_table(case, center_hi) {
                    let (ea, pa) = pt(e0);
                    let (eb, pb) = pt(e1);
                    segs[(n - n_lo) as usize].push((ea, eb, pa, pb));
                }
            }
        }
    }
    segs.into_iter().enumerate().filter(|(_, s)| !s.is_empty()).map(|(li, s)| Level { n: n_lo + li as i32, lines: chain(&s) }).collect()
}

fn chain(segs: &[(u32, u32, [f32; 2], [f32; 2])]) -> Vec<Vec<[f32; 2]>> {
    let mut adj: HashMap<u32, [u32; 2]> = HashMap::with_capacity(segs.len() * 2);
    const NONE: u32 = u32::MAX;
    for (si, s) in segs.iter().enumerate() {
        for e in [s.0, s.1] {
            let ent = adj.entry(e).or_insert([NONE, NONE]);
            if ent[0] == NONE { ent[0] = si as u32 } else { ent[1] = si as u32 }
        }
    }
    let mut used = vec![false; segs.len()];
    let mut out = vec![];
    let next = |adj: &HashMap<u32, [u32; 2]>, edge: u32, cur: usize, used: &Vec<bool>| -> Option<usize> {
        adj.get(&edge).and_then(|a| a.iter().find(|&&s| s != NONE && s as usize != cur && !used[s as usize]).map(|&s| s as usize))
    };
    for s0 in 0..segs.len() {
        if used[s0] { continue; }
        used[s0] = true;
        let mut pts: std::collections::VecDeque<[f32; 2]> = std::collections::VecDeque::new();
        pts.push_back(segs[s0].2);
        pts.push_back(segs[s0].3);
        // forward from edge b
        let (mut cur, mut edge) = (s0, segs[s0].1);
        while let Some(n) = next(&adj, edge, cur, &used) {
            used[n] = true;
            let s = &segs[n];
            if s.0 == edge { pts.push_back(s.3); edge = s.1; } else { pts.push_back(s.2); edge = s.0; }
            cur = n;
        }
        // backward from edge a
        let (mut cur, mut edge) = (s0, segs[s0].0);
        while let Some(n) = next(&adj, edge, cur, &used) {
            used[n] = true;
            let s = &segs[n];
            if s.0 == edge { pts.push_front(s.3); edge = s.1; } else { pts.push_front(s.2); edge = s.0; }
            cur = n;
        }
        out.push(pts.into_iter().collect());
    }
    out
}

/// 3x3 box blur of a bool mask to f32 (smooth electrode outlines).
pub fn blur_mask(m: &[bool], nx: usize, ny: usize) -> Vec<f32> {
    let mut out = vec![0f32; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let mut s = 0.0;
            let mut c = 0.0;
            for dj in -1i32..=1 { for di in -1i32..=1 {
                let (ii, jj) = (i as i32 + di, j as i32 + dj);
                if ii < 0 || jj < 0 || ii >= nx as i32 || jj >= ny as i32 { s += m[(i.min(nx - 1)) + nx * j.min(ny - 1)] as i32 as f32; } else { s += m[ii as usize + nx * jj as usize] as i32 as f32; }
                c += 1.0;
            } }
            out[i + nx * j] = s / c;
        }
    }
    out
}

/// Merge a bool grid into rectangles: runs along i per row k, merged across identical consecutive rows.
/// Returns (i0, i1_exclusive, k0, k1_exclusive).
pub fn rects(m: &[bool], nx: usize, nz: usize) -> Vec<(usize, usize, usize, usize)> {
    let mut open: Vec<(usize, usize, usize)> = vec![]; // (i0, i1, k0)
    let mut out = vec![];
    for k in 0..=nz {
        let mut runs: Vec<(usize, usize)> = vec![];
        if k < nz {
            let mut i = 0;
            while i < nx {
                if m[i + nx * k] { let s = i; while i < nx && m[i + nx * k] { i += 1; } runs.push((s, i)); } else { i += 1; }
            }
        }
        let mut still = vec![];
        for o in open.drain(..) {
            if runs.contains(&(o.0, o.1)) { still.push(o); } else { out.push((o.0, o.1, o.2, k)); }
        }
        for r in runs { if !still.iter().any(|o| (o.0, o.1) == r) { still.push((r.0, r.1, k)); } }
        open = still;
    }
    out
}

/// "Nice" contour spacing [V] giving at most `max_levels` levels over a span.
pub fn nice_step(span: f64, max_levels: f64) -> f64 {
    let raw = span / max_levels;
    let p = 10f64.powf(raw.log10().floor());
    for m in [1.0, 2.0, 5.0, 10.0] { if m * p >= raw { return m * p; } }
    10.0 * p
}
