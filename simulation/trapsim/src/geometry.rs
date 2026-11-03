//! Parametric geometry of the 486-quadrupole MultiQ-IT cube (all lengths in mm).
//!
//! Lattice: N x N x N lattice points at c0 + pitch*i; a post sits on every *surface* lattice point (488 for N=10),
//! its axis normal to the face it belongs to (z faces first, then y, then x: 10x10 + 10x10 + 10x8 + 10x8 + 8x8 + 8x8).
//! Each post is a cylinder (diameter d, height `height`) standing `gap` above a conducting wall plane of thickness `t_pcb`;
//! a stem (radius stem_r) carries it through a clearance hole (radius clear_r) in the plane (their code: 2 mm stems,
//! clearance = d/2).  RF sign = (-1)^(i+j+k).  The six wall planes form a closed shell, so the cube interior is only
//! open through the ports.  A port is one quadrupole cell of a face: its four posts are replaced by rods of the same
//! diameter that continue `l_ext` mm outside the wall through enlarged clearance holes, and the plane is cut away in a
//! disk of radius `ap_r` about the cell axis (so ions can pass).  Ports carry their own DC (q_in, q_out).

pub const ID_FREE: u8 = 0;
pub const ID_POST: u8 = 1; // ordinary posts and stems (basis: U0)
pub const ID_PLANE: u8 = 2; // wall planes (basis: Uwall)
pub const ID_ENT: u8 = 3; // entrance-port rods (q_in)
pub const ID_EXT: u8 = 4; // exit-port rods (q_out)

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PortSpec {
    pub axis: usize,   // face normal axis 0..3
    pub side: usize,   // 0: face at coordinate 0, 1: face at coordinate L
    pub cell: [usize; 2], // quadrupole cell indices (0..n-2) along in-plane axes (b, c) = ((a+1)%3, (a+2)%3)
}

#[derive(Clone, Debug)]
pub struct Geometry {
    pub n: usize,
    pub d: f64,
    pub pitch: f64,
    pub height: f64,
    pub t_pcb: f64,
    pub gap: f64,
    pub stem_r: f64,
    pub clear_r: f64,
    pub l_ext: f64,    // port extension length outside the wall
    pub ap_r: f64,     // aperture radius cut in the wall plane at a port
    pub rod_clear_extra: f64, // port rod clearance hole radius = d/2 + this
    pub cap_gap: f64,  // distance from the rod ends to the end cap (potential of the downstream/upstream guide)
    pub cap_r: f64,    // end-cap radius
    pub entrance: PortSpec,
    pub exit: PortSpec,
    pub array: bool,
}

impl Geometry {
    /// As-built (photo forensics, SM Fig. S1): D 6.34 mm on 8.30 mm pitch.
    pub fn as_built() -> Geometry {
        let d = 6.34;
        let pitch = 8.30;
        Geometry {
            array: false, n: 10, d, pitch, height: d, t_pcb: 1.5, gap: (pitch - d) / 2.0, stem_r: 2.0, clear_r: d / 2.0,
            l_ext: 6.0, ap_r: 4.0, rod_clear_extra: 0.5, cap_gap: 2.0, cap_r: 12.0,
            entrance: PortSpec { axis: 1, side: 1, cell: [4, 4] },
            exit: PortSpec { axis: 2, side: 1, cell: [4, 4] },
        }
    }
    /// Their code's dimensions: De = 6.04686, pitch 8.0, 83 mm box.
    pub fn code_dims() -> Geometry {
        let de = (83.0 - 3.0) / 10.0 / 1.323;
        let pitch = de * 1.323;
        Geometry { d: de, pitch, height: de, gap: (pitch - de) / 2.0, clear_r: de / 2.0, ..Geometry::as_built() }
    }
    /// Snap the geometry to a grid of target spacing `h_target` [mm]: the spacing becomes pitch/(2*round(pitch/(2 h))) so that
    /// the quadrupole axes and post axes fall exactly on grid nodes (symmetric staircases), and the plane thickness and
    /// port extension become whole numbers of cells.  Returns (geometry, h).
    pub fn snap(&self, h_target: f64) -> (Geometry, f64) {
        let nh = (self.pitch / (2.0 * h_target)).round().max(1.0);
        let h = self.pitch / (2.0 * nh);
        let mut g = self.clone();
        g.t_pcb = (self.t_pcb / h).round() * h;
        g.l_ext = (self.l_ext / h).round() * h;
        (g, h)
    }
    pub fn c0(&self) -> f64 { self.t_pcb + self.gap + self.d / 2.0 }
    pub fn size(&self) -> f64 { 2.0 * self.c0() + (self.n - 1) as f64 * self.pitch }
    pub fn lat(&self, i: usize) -> f64 { self.c0() + self.pitch * i as f64 }
    /// ideal quadrupole r0 of four posts at the corners of a square of side `pitch`
    pub fn r0(&self) -> f64 { self.pitch / 2f64.sqrt() - self.d / 2.0 }
    pub fn post_top_depth(&self) -> f64 { self.t_pcb + self.gap + self.height }

    /// Does lattice point `idx` carry a post belonging to face axis `a` (priority z > y > x)?
    pub fn post_valid(&self, a: usize, idx: [usize; 3]) -> bool {
        let n = self.n;
        if idx[a] != 0 && idx[a] != n - 1 { return false; }
        for a2 in (a + 1)..3 { if idx[a2] == 0 || idx[a2] == n - 1 { return false; } }
        true
    }
    pub fn n_posts(&self) -> usize {
        let mut c = 0;
        for i in 0..self.n { for j in 0..self.n { for k in 0..self.n {
            let idx = [i, j, k];
            if (0..3).any(|a| self.post_valid(a, idx)) { c += 1; }
        }}}
        c
    }
    /// Port rod lattice index sets.
    fn port_posts(&self, p: &PortSpec) -> [[usize; 3]; 4] {
        let a = p.axis; let b = (a + 1) % 3; let c = (a + 2) % 3;
        let mut out = [[0usize; 3]; 4];
        let mut m = 0;
        for db in 0..2 { for dc in 0..2 {
            let mut idx = [0usize; 3];
            idx[a] = if p.side == 0 { 0 } else { self.n - 1 };
            idx[b] = p.cell[0] + db; idx[c] = p.cell[1] + dc;
            out[m] = idx; m += 1;
        }}
        out
    }
    /// Centre of the port quadrupole axis in cube coordinates (the axial coordinate is the outer face).
    pub fn port_axis_xy(&self, p: &PortSpec) -> (usize, usize, f64, f64) {
        let a = p.axis; let b = (a + 1) % 3; let c = (a + 2) % 3;
        (b, c, 0.5 * (self.lat(p.cell[0]) + self.lat(p.cell[0] + 1)), 0.5 * (self.lat(p.cell[1]) + self.lat(p.cell[1] + 1)))
    }
    pub fn port(&self, index: usize) -> PortSpec {
        if !self.array { return self.exit; }
        PortSpec { cell: [2 + 2 * (index % 3), 2 + 2 * (index / 3)], ..self.exit }
    }
    // The nine mouths are separated by two pitches. Choose the nearest mouth
    // only for classifying the local metal; this never chooses an ion's destination.
    pub fn nearest_port(&self, x: f64, y: f64) -> usize {
        if !self.array { return 0; }
        let first = (self.lat(2) + self.lat(3)) / 2.;
        let ix = |x: f64| ((x-first)/(2.*self.pitch)).round().clamp(0.,2.) as usize;
        ix(x) + 3*ix(y)
    }
    /// Grid-domain bounds (cube coordinates) [lo, hi] per axis: the cube plus the port extensions (+ margin).
    pub fn domain(&self, margin: f64) -> ([f64; 3], [f64; 3]) {
        let l = self.size();
        let mut lo = [0.0; 3]; let mut hi = [l; 3];
        for p in [&self.entrance, &self.exit] {
            if p.side == 0 { lo[p.axis] = -(self.l_ext + margin); } else { hi[p.axis] = l + self.l_ext + margin; }
        }
        (lo, hi)
    }

    /// Classify a point (cube coordinates, mm). `g` = electrode growth: convex electrodes grow by g, holes shrink by g.
    /// Returns (id, rf_sign).
    pub fn classify(&self, p: [f64; 3], g: f64) -> (u8, i8) {
        let l = self.size();
        let n = self.n;
        let half_d = self.d / 2.0;
        let mut found: (u8, i8) = (ID_FREE, 0);
        let mut plane = false;
        let inside = (0..3).all(|a| p[a] >= -1e-12 && p[a] <= l + 1e-12);
        for a in 0..3 {
            let b = (a + 1) % 3; let c = (a + 2) % 3;
            for side in 0..2usize {
                let u = if side == 0 { p[a] } else { l - p[a] };
                let selected;
                let port_id = if self.entrance.axis == a && self.entrance.side == side { Some((ID_ENT, &self.entrance)) }
                    else if self.exit.axis == a && self.exit.side == side { selected = self.port(self.nearest_port(p[b],p[c])); Some((ID_EXT, &selected)) } else { None };
                // outside the cube: only port rods
                if u < 0.0 || (port_id.is_some() && u < 0.0 + 1e-12) {
                    if let Some((id, ps)) = port_id {
                        // end cap: continues the port rods' potential beyond the rods (downstream/upstream guide)
                        if u <= -(self.l_ext + self.cap_gap) {
                            let (_, _, cb, cc) = self.port_axis_xy(ps);
                            if ((p[b] - cb).powi(2) + (p[c] - cc).powi(2)).sqrt() < self.cap_r { return (id, 0); }
                        }
                        if u >= -(self.l_ext + g) {
                            for idx in self.port_posts(ps) {
                                let dd = ((p[b] - self.lat(idx[b])).powi(2) + (p[c] - self.lat(idx[c])).powi(2)).sqrt();
                                if dd < half_d + g { return (id, if (idx[0] + idx[1] + idx[2]) % 2 == 0 { 1 } else { -1 }); }
                            }
                        }
                    }
                    continue;
                }
                if !inside || u > self.post_top_depth() + g { continue; }
                // nearest lattice posts in-plane
                let ib = (((p[b] - self.c0()) / self.pitch).round().max(0.0) as usize).min(n - 1);
                let ic = (((p[c] - self.c0()) / self.pitch).round().max(0.0) as usize).min(n - 1);
                let mut idx = [0usize; 3];
                idx[a] = if side == 0 { 0 } else { n - 1 }; idx[b] = ib; idx[c] = ic;
                let valid = self.post_valid(a, idx);
                let dd = ((p[b] - self.lat(ib)).powi(2) + (p[c] - self.lat(ic)).powi(2)).sqrt();
                let sign: i8 = if (idx[0] + idx[1] + idx[2]) % 2 == 0 { 1 } else { -1 };
                let mut port_post: Option<u8> = None;
                let mut aperture = false;
                if let Some((id, ps)) = port_id {
                    if valid && self.port_posts(ps).contains(&idx) { port_post = Some(id); }
                    let (_, _, cb, cc) = self.port_axis_xy(ps);
                    if ((p[b] - cb).powi(2) + (p[c] - cc).powi(2)).sqrt() < self.ap_r - g { aperture = true; }
                }
                if let Some(id) = port_post {
                    // rod: full radius from the wall (depth 0) to the post top; replaces the stem
                    if u <= self.post_top_depth() + g && dd < half_d + g { return (id, sign); }
                }
                // wall plane
                if u <= self.t_pcb + g {
                    let hole_r = if port_post.is_some() { half_d + self.rod_clear_extra } else { self.clear_r };
                    let in_hole = aperture || (valid && dd < hole_r - g);
                    if !in_hole { plane = true; }
                }
                // stem + post
                if valid && port_post.is_none() {
                    if u <= self.t_pcb + self.gap + g && dd < self.stem_r + g { found = (ID_POST, sign); }
                    if u >= self.t_pcb + self.gap - g && u <= self.post_top_depth() + g && dd < half_d + g { found = (ID_POST, sign); }
                }
            }
        }
        if found.0 != ID_FREE { return found; }
        if plane { return (ID_PLANE, 0); }
        (ID_FREE, 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn post_counts() {
        let g = Geometry::as_built();
        assert_eq!(g.n_posts(), 488);
        assert!((g.size() - 86.0).abs() < 1e-9);
        let c = Geometry::code_dims();
        assert!((c.size() - 83.0).abs() < 1e-9, "{}", c.size());
        assert!((c.r0() - 2.6334).abs() < 1e-3);
        assert!((g.r0() - 2.699).abs() < 1e-3);
    }
    #[test]
    fn quad_axis_is_open_until_plane() {
        let g = Geometry::code_dims();
        // a wall quadrupole cell not used by a port: axis between lattice (1,1) and (2,2) on the -z face
        let (x, y) = (0.5 * (g.lat(1) + g.lat(2)), 0.5 * (g.lat(1) + g.lat(2)));
        assert_eq!(g.classify([x, y, 6.0], 0.0).0, ID_FREE);
        assert_eq!(g.classify([x, y, 1.0], 0.0).0, ID_PLANE);
        assert_eq!(g.classify([g.lat(1), g.lat(1), 5.0], 0.0).0, ID_POST);
        assert_eq!(g.classify([g.lat(1), g.lat(1), 1.0], 0.0).0, ID_POST); // stem through the plane
        assert_eq!(g.classify([g.lat(1) + 2.5, g.lat(1), 1.0], 0.0).0, ID_FREE); // clearance annulus
    }
}
