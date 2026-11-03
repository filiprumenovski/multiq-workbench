//! Painters for the three spatial views: main x-z topographic side view, end view, frequency terrain strip.
use crate::engine::Snapshot;
use crate::terrain::Terrain;
use crate::theme::{mono, sans, species_color, with_alpha, Palette};
use crate::topo::{self, Level};
use crate::world::*;
use egui::{self, Align2, Color32, Painter, Pos2, Rect, Shape, Stroke, Vec2};

pub const Z0: f32 = -0.5;
pub const Z1: f32 = 20.7;

/// Cached solid-electrode geometry (grid coordinates).
pub struct Statics {
    pub cut_rects: Vec<(usize, usize, usize, usize)>,
    pub cut_lines: Vec<Vec<[f32; 2]>>,
    pub ghost_rects: Vec<(usize, usize, usize, usize)>,
    pub end_rects: Vec<(usize, usize, usize, usize)>,
    pub end_lines: Vec<Vec<[f32; 2]>>,
}

pub fn build_statics(v: &VizTables) -> Statics {
    let (nx, nz, ny) = (v.nx, v.nz, v.ny);
    let all = vec![true; nx * nz.max(ny)];
    let lines = |m: &[bool], n1: usize, n2: usize| -> Vec<Vec<[f32; 2]>> {
        let b = topo::blur_mask(m, n1, n2);
        topo::isolines(&b, &all, n1, n2, 0.5, 1, 1).into_iter().flat_map(|l| l.lines).collect()
    };
    let rods = |k: usize| { let z = k as f64 * v.h_mm; (z >= v.layers[0] && z <= v.layers[1]) || (z >= v.layers[4] && z <= v.layers[5]) };
    let ghost: Vec<bool> = (0..nx * nz).map(|n| v.shadow[n] && !v.cut[n] && rods(n / nx)).collect();
    Statics {
        cut_rects: topo::rects(&v.cut, nx, nz),
        cut_lines: lines(&v.cut, nx, nz),
        ghost_rects: topo::rects(&ghost, nx, nz),
        end_rects: topo::rects(&v.end_mask, nx, ny),
        end_lines: lines(&v.end_mask, nx, ny),
    }
}

/// Contour cache for the selected species' effective potential.
pub struct Topo {
    pub key: (u64, u64, usize, u64),
    pub eff: Vec<f32>,
    pub valid: Vec<bool>,
    pub step: f32,
    pub levels: Vec<Level>,
}

pub fn build_topo(v: &VizTables, sp: &Sp, v_rf: f64, dc: [f64; 4]) -> Topo {
    let eff = topo::eff_potential(v, sp, v_rf, dc);
    let valid: Vec<bool> = v.cut.iter().map(|c| !c).collect();
    // contour spacing from the ion-accessible valley (|x| < 1.6 mm) so the gate barrier and axial steps resolve;
    // levels run on up the flanks to the summits (|x| < 4.6 mm), capped at 100 levels
    let (mut lo, mut hi, mut hi_all) = (f32::INFINITY, f32::NEG_INFINITY, f32::NEG_INFINITY);
    let a = (v.nx / 2) as i64;
    let (w1, w2) = ((1.6 / v.h_mm) as i64, (4.6 / v.h_mm) as i64);
    for k in 0..v.nz {
        for di in -w2..=w2 {
            let i = (a + di) as usize;
            if !valid[i + v.nx * k] { continue; }
            let u = eff[i + v.nx * k];
            hi_all = hi_all.max(u);
            if di.abs() <= w1 { lo = lo.min(u); hi = hi.max(u); }
        }
    }
    if !lo.is_finite() { lo = 0.0; hi = 1.0; hi_all = 2.0; }
    let step = topo::nice_step((hi - lo) as f64, 14.0).max(0.1) as f32;
    let n_lo = (lo / step).floor() as i32 - 1;
    let n_hi = ((hi_all / step).ceil() as i32).min(n_lo + 100);
    let levels = topo::isolines(&eff, &valid, v.nx, v.nz, step, n_lo, n_hi);
    Topo { key: (0, 0, 0, 0), eff, valid, step, levels }
}

/// Screen mapping of the z axis shared by the main view and the terrain strip.
#[derive(Clone, Copy)]
pub struct Frame { pub rect: Rect, pub s: f32 }
impl Frame {
    pub fn new(rect: Rect, min_xspan: f32) -> Frame {
        let s = (rect.width() / (Z1 - Z0)).min(rect.height() / min_xspan);
        Frame { rect, s }
    }
    fn zspan_px(&self) -> f32 { (Z1 - Z0) * self.s }
    pub fn left(&self) -> f32 { self.rect.center().x - self.zspan_px() * 0.5 }
    pub fn zx(&self, z_mm: f32) -> f32 { self.left() + (z_mm - Z0) * self.s }
    pub fn xy(&self, x_mm: f32) -> f32 { self.rect.center().y - x_mm * self.s }
    pub fn z_of(&self, px: f32) -> f32 { (px - self.left()) / self.s + Z0 }
    pub fn x_of(&self, py: f32) -> f32 { (self.rect.center().y - py) / self.s }
    pub fn grid(&self, v: &VizTables, i: f32, k: f32) -> Pos2 {
        Pos2::new(self.zx(k * v.h_mm as f32), self.xy(i * v.h_mm as f32 - v.axis_mm as f32))
    }
}


/// Cheap batched ion drawing: one mesh for all dots and trails (no per-shape tessellation).
pub struct IonMesh { pub m: egui::Mesh }
impl IonMesh {
    pub fn new() -> IonMesh { IonMesh { m: egui::Mesh::default() } }
    fn vtx(&mut self, p: Pos2, c: Color32) -> u32 {
        self.m.vertices.push(egui::epaint::Vertex { pos: p, uv: egui::epaint::WHITE_UV, color: c });
        (self.m.vertices.len() - 1) as u32
    }
    pub fn dot(&mut self, c: Pos2, r: f32, col: Color32) {
        let i0 = self.vtx(c, col);
        let mut prev = 0;
        for k in 0..=6 {
            let a = k as f32 * std::f32::consts::TAU / 6.0;
            let i = if k == 6 { i0 + 1 } else { self.vtx(c + Vec2::new(r * a.cos(), r * a.sin()), col) };
            if k > 0 { self.m.indices.extend_from_slice(&[i0, prev, i]); }
            prev = i;
        }
    }
    pub fn seg(&mut self, a: Pos2, b: Pos2, w: f32, ca: Color32, cb: Color32) {
        let d = b - a;
        let l = d.length();
        if l < 1e-3 { return; }
        let n = Vec2::new(-d.y, d.x) / l * (w * 0.5);
        let (i0, i1, i2, i3) = (self.vtx(a + n, ca), self.vtx(a - n, ca), self.vtx(b - n, cb), self.vtx(b + n, cb));
        self.m.indices.extend_from_slice(&[i0, i1, i2, i0, i2, i3]);
    }
    pub fn finish(self, p: &Painter) { if !self.m.is_empty() { p.add(Shape::mesh(self.m)); } }
}

pub struct MainOut { pub clicked_species: Option<usize>, pub hover: Option<(f32, f32, f32)> }

pub struct Scene<'a> {
    pub pal: &'a Palette,
    pub viz: &'a VizTables,
    pub st: &'a Statics,
    pub topo: &'a Topo,
    pub snap: &'a Snapshot,
    pub species: &'a [Sp],
    pub target: usize,
    pub eff_sp: usize,
    pub f_ac: f64,
    pub terr: &'a Terrain,
    pub show_trails: bool,
    pub show_ghost: bool,
    pub show_ions: bool,
}

fn rect_px(f: &Frame, v: &VizTables, r: (usize, usize, usize, usize)) -> Rect {
    // r = (i0, i1, k0, k1) in nodes; cells are centred on nodes
    let (a, b) = (f.grid(v, r.0 as f32 - 0.5, r.2 as f32 - 0.5), f.grid(v, r.1 as f32 - 0.5, r.3 as f32 - 0.5));
    Rect::from_two_pos(a, b)
}

fn zticks(p: &Painter, f: &Frame, pal: &Palette, y: f32, labels: bool) {
    let mut z = 0.0;
    while z <= 20.01 {
        let x = f.zx(z);
        p.line_segment([Pos2::new(x, y - 3.0), Pos2::new(x, y)], Stroke::new(1.0, pal.faint));
        if labels { p.text(Pos2::new(x, y - 5.0), Align2::CENTER_BOTTOM, format!("{}", z as i32), mono(9.5), pal.dim); }
        z += 2.0;
    }
}

pub fn draw_main(ui: &mut egui::Ui, rect: Rect, sc: &Scene) -> MainOut {
    let pal = sc.pal;
    let f = Frame::new(rect, 13.0);
    let resp = ui.interact(rect, ui.id().with("main"), egui::Sense::click());
    let p = ui.painter_at(rect);
    p.rect_filled(rect, 0.0, pal.bg);
    let v = sc.viz;
    // faint z grid
    let mut z = 0.0;
    while z <= 20.01 { p.line_segment([Pos2::new(f.zx(z), rect.top()), Pos2::new(f.zx(z), rect.bottom())], Stroke::new(1.0, with_alpha(pal.line, 90))); z += 2.0; }
    p.line_segment([Pos2::new(rect.left(), f.xy(0.0)), Pos2::new(rect.right(), f.xy(0.0))], Stroke::new(1.0, with_alpha(pal.line, 140)));

    // projected rods
    if sc.show_ghost {
        for r in &sc.st.ghost_rects { p.rect_filled(rect_px(&f, v, *r), 0.0, pal.ghost); }
    }
    // topographic contours (thin) with every 5th an index contour
    for lv in &sc.topo.levels {
        let idx = lv.n.rem_euclid(5) == 0;
        let stroke = if idx { Stroke::new(1.3, pal.index) } else { Stroke::new(0.8, pal.contour) };
        for l in &lv.lines {
            if l.len() < 2 { continue; }
            let pts: Vec<Pos2> = l.iter().map(|q| f.grid(v, q[0], q[1])).collect();
            p.add(Shape::line(pts, stroke));
        }
    }
    // solid electrodes from the conductor mask
    for r in &sc.st.cut_rects { p.rect_filled(rect_px(&f, v, *r), 0.0, pal.elec_fill); }
    for l in &sc.st.cut_lines {
        let pts: Vec<Pos2> = l.iter().map(|q| f.grid(v, q[0], q[1])).collect();
        p.add(Shape::line(pts, Stroke::new(1.2, pal.elec_line)));
    }
    // index labels, knocked out of the line work
    let mut placed: Vec<Rect> = vec![];
    for lv in sc.topo.levels.iter().filter(|l| l.n.rem_euclid(5) == 0) {
        let val = lv.n as f32 * sc.topo.step;
        let txt = if sc.topo.step < 0.1 { format!("{:.2}", val) } else if sc.topo.step < 1.0 { format!("{:.1}", val) } else { format!("{:.0}", val) };
        let mut cands: Vec<(f32, Pos2)> = vec![];
        for l in &lv.lines {
            let mut acc = 0.0;
            let mut prev: Option<Pos2> = None;
            for q in l {
                let pt = f.grid(v, q[0], q[1]);
                if let Some(pp) = prev { acc += pt.distance(pp); }
                prev = Some(pt);
                if acc > 40.0 { acc = 0.0; cands.push(((pt.y - (rect.bottom() - 26.0)).abs(), pt)); }
            }
        }
        cands.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        for (_, pt) in cands {
            let g = p.layout_no_wrap(txt.clone(), mono(9.5), pal.index);
            let r = Rect::from_center_size(pt, g.size() + Vec2::new(8.0, 2.0));
            if !rect.shrink(6.0).contains_rect(r) || placed.iter().any(|o| o.expand(8.0).intersects(r)) { continue; }
            p.rect_filled(r, 2.0, pal.bg);
            p.galley(r.center() - g.size() * 0.5, g, pal.index);
            placed.push(r);
            break;
        }
    }
    // reference planes
    let dash = |x: f32, col: Color32| { for s in Shape::dashed_line(&[Pos2::new(x, rect.top() + 22.0), Pos2::new(x, rect.bottom() - 4.0)], Stroke::new(1.0, col), 4.0, 4.0) { p.add(s); } };
    dash(f.zx(Z_EXIT as f32 * 1e3), with_alpha(pal.dim, 160));
    dash(f.zx(Z_RETURN as f32 * 1e3), with_alpha(pal.faint, 200));
    let lab = |z: f32, t: &str| { p.text(Pos2::new(f.zx(z), rect.top() + 6.0), Align2::CENTER_TOP, t, mono(9.0), pal.faint); };
    lab(1.2, "TRAP");
    if sc.show_ghost { p.text(Pos2::new(f.zx(3.4), f.xy(-6.6)), Align2::LEFT_BOTTOM, "rods, projected", mono(8.5), pal.faint); } lab(9.0, "Q1 CHANNEL"); lab(16.8, "GATE"); lab(19.7, "TO 2ND TRAP");
    zticks(&p, &f, pal, rect.bottom(), true);
    for xm in [-6.0f32, -4.0, -2.0, 2.0, 4.0, 6.0] {
        let y = f.xy(xm);
        if y > rect.top() + 40.0 && y < rect.bottom() - 14.0 {
            p.line_segment([Pos2::new(rect.left(), y), Pos2::new(rect.left() + 4.0, y)], Stroke::new(1.0, pal.faint));
            p.text(Pos2::new(rect.left() + 7.0, y), Align2::LEFT_CENTER, format!("{:+}", xm as i32), mono(9.0), pal.faint);
        }
    }
    p.text(Pos2::new(rect.left() + 7.0, f.xy(0.0) - 2.0), Align2::LEFT_BOTTOM, "x [mm]", mono(9.0), pal.dim);
    p.text(Pos2::new(rect.right() - 6.0, rect.bottom() - 16.0), Align2::RIGHT_BOTTOM, "z [mm]", mono(9.0), pal.dim);

    // excitation planes: where each species' secular frequency crosses the waterline
    for (s, _) in sc.species.iter().enumerate() {
        let col = species_color(pal.dark, s);
        for zc in Terrain::crossings(&sc.terr.z_mm, &sc.terr.fx[s], sc.f_ac) {
            let x = f.zx(zc as f32);
            p.line_segment([Pos2::new(x, rect.top() + 20.0), Pos2::new(x, rect.top() + 28.0)], Stroke::new(1.5, col));
        }
    }

    // ions: trails then dots, batched into one mesh
    let snap = sc.snap;
    let ions: &[crate::engine::IonVis] = if sc.show_ions { &snap.ions } else { &[] };
    let pos = |q: [f32; 3]| Pos2::new(f.zx(q[2]), f.xy(q[0]));
    let mut mesh = IonMesh::new();
    if sc.show_trails {
        for ion in ions {
            let col = if ion.hot { pal.accent } else { species_color(pal.dark, ion.sp as usize) };
            let mut prev = pos(ion.p);
            for (k, t) in ion.trail.iter().enumerate() {
                let pt = pos(*t);
                if pt.distance(prev) > 0.3 && pt.distance(prev) < 120.0 {
                    let (a0, a1) = ([130u8, 84, 50, 26][k.min(3)], [84u8, 50, 26, 0][k.min(3)]);
                    mesh.seg(prev, pt, 1.0, with_alpha(col, a0), with_alpha(col, a1));
                }
                prev = pt;
            }
        }
    }
    for ion in ions {
        let c = pos(ion.p);
        if !rect.contains(c) { continue; }
        if ion.hot {
            mesh.dot(c, 5.2, with_alpha(pal.accent, 30));
            mesh.dot(c, 2.2, pal.accent);
        } else {
            mesh.dot(c, 1.7, with_alpha(species_color(pal.dark, ion.sp as usize), 235));
        }
    }
    mesh.finish(&p);

    // legend (click = show this species' potential and make it the selectivity target)
    let mut clicked = None;
    let mut x = rect.left() + 12.0;
    let y = rect.top() + 26.0 + 14.0;
    for (s, sp) in sc.species.iter().enumerate() {
        let col = species_color(pal.dark, s);
        let txt = format!("{}", sp.label);
        let g = p.layout_no_wrap(txt, mono(10.5), pal.text);
        let r = Rect::from_min_size(Pos2::new(x, y - 8.0), Vec2::new(g.size().x + 22.0, 16.0));
        let sel = s == sc.eff_sp;
        p.rect_filled(r, 3.0, if sel { with_alpha(col, 40) } else { with_alpha(pal.card, 200) });
        if sel { p.rect_stroke(r, 3.0, Stroke::new(1.0, with_alpha(col, 160)), egui::StrokeKind::Inside); }
        p.circle_filled(Pos2::new(r.left() + 9.0, y), 3.0, col);
        p.galley(Pos2::new(r.left() + 17.0, y - g.size().y * 0.5), g, pal.text);
        if s == sc.target { p.text(Pos2::new(r.right() - 3.0, y), Align2::RIGHT_CENTER, "◆", mono(7.0), pal.accent); }
        if resp.clicked() { if let Some(pp) = resp.interact_pointer_pos() { if r.contains(pp) { clicked = Some(s); } } }
        x = r.right() + 6.0;
    }
    p.text(Pos2::new(rect.right() - 10.0, rect.top() + 40.0), Align2::RIGHT_CENTER, "x–z slice · contours: effective potential [V/charge] of the ◆ species", mono(9.0), pal.faint);
    let _ = x;

    // hover readout: potential and position
    let mut hover = None;
    if let Some(hp) = resp.hover_pos().or_else(|| ui.input(|i| i.pointer.hover_pos())) {
        if rect.contains(hp) {
            let (zmm, xmm) = (f.z_of(hp.x), f.x_of(hp.y));
            let (fi, fk) = ((xmm as f64 + v.axis_mm) / v.h_mm, zmm as f64 / v.h_mm);
            if fi >= 0.0 && fk >= 0.0 && (fi as usize) < v.nx - 1 && (fk as usize) < v.nz - 1 {
                let (i, k) = (fi as usize, fk as usize);
                let (tx, tk) = ((fi - i as f64) as f32, (fk - k as f64) as f32);
                let g = |ii: usize, kk: usize| sc.topo.eff[ii + v.nx * kk];
                let u = (1.0 - tx) * (1.0 - tk) * g(i, k) + tx * (1.0 - tk) * g(i + 1, k) + (1.0 - tx) * tk * g(i, k + 1) + tx * tk * g(i + 1, k + 1);
                let inside = v.cut[i + v.nx * k];
                hover = Some((zmm, xmm, u));
                p.line_segment([Pos2::new(hp.x, rect.top()), Pos2::new(hp.x, rect.bottom())], Stroke::new(1.0, with_alpha(pal.dim, 60)));
                p.line_segment([Pos2::new(rect.left(), hp.y), Pos2::new(rect.right(), hp.y)], Stroke::new(1.0, with_alpha(pal.dim, 60)));
                let txt = if inside { format!("z {:5.2} mm   x {:+5.2} mm   electrode", zmm, xmm) } else { format!("z {:5.2} mm   x {:+5.2} mm   U {:6.2} V", zmm, xmm, u) };
                let g = p.layout_no_wrap(txt, mono(10.5), pal.text);
                let mut r = Rect::from_min_size(hp + Vec2::new(14.0, 14.0), g.size() + Vec2::new(12.0, 8.0));
                if r.right() > rect.right() - 4.0 { r = r.translate(Vec2::new(-(r.width() + 28.0), 0.0)); }
                if r.bottom() > rect.bottom() - 4.0 { r = r.translate(Vec2::new(0.0, -(r.height() + 28.0))); }
                p.rect_filled(r, 3.0, with_alpha(pal.panel, 235));
                p.rect_stroke(r, 3.0, Stroke::new(1.0, pal.line), egui::StrokeKind::Inside);
                p.galley(r.min + Vec2::new(6.0, 4.0), g, pal.text);
            }
        }
    }
    MainOut { clicked_species: clicked, hover }
}

/// End view (x-y): electrodes at z = 9 mm, ions projected, ring at the excitation threshold.
pub fn draw_end(ui: &mut egui::Ui, rect: Rect, sc: &Scene) {
    let pal = sc.pal;
    let v = sc.viz;
    let p = ui.ctx().layer_painter(ui.layer_id()).with_clip_rect(rect);
    p.rect_filled(rect, 4.0, with_alpha(pal.panel, 245));
    p.rect_stroke(rect, 4.0, Stroke::new(1.0, pal.line), egui::StrokeKind::Inside);
    let half = 7.0f32; // mm shown each side
    let s = rect.width().min(rect.height()) / (2.0 * half);
    let c = rect.center();
    let map = |x: f32, y: f32| Pos2::new(c.x + x * s, c.y - y * s);
    let node = |i: f32, j: f32| map(i * v.h_mm as f32 - v.axis_mm as f32, j * v.h_mm as f32 - v.axis_mm as f32);
    for r in &sc.st.end_rects {
        let rr = Rect::from_two_pos(node(r.0 as f32 - 0.5, r.2 as f32 - 0.5), node(r.1 as f32 - 0.5, r.3 as f32 - 0.5));
        p.rect_filled(rr, 0.0, pal.elec_fill);
    }
    for l in &sc.st.end_lines {
        let pts: Vec<Pos2> = l.iter().map(|q| node(q[0], q[1])).collect();
        p.add(Shape::line(pts, Stroke::new(1.0, pal.elec_line)));
    }
    p.line_segment([Pos2::new(rect.left() + 4.0, c.y), Pos2::new(rect.right() - 4.0, c.y)], Stroke::new(1.0, with_alpha(pal.line, 160)));
    p.line_segment([Pos2::new(c.x, rect.top() + 4.0), Pos2::new(c.x, rect.bottom() - 4.0)], Stroke::new(1.0, with_alpha(pal.line, 160)));
    // threshold ring r = 1.2 mm
    let n = 48;
    let ring: Vec<Pos2> = (0..=n).map(|k| { let a = k as f32 / n as f32 * std::f32::consts::TAU; map(R_EXCITED as f32 * 1e3 * a.cos(), R_EXCITED as f32 * 1e3 * a.sin()) }).collect();
    for sh in Shape::dashed_line(&ring, Stroke::new(1.0, with_alpha(pal.accent, 200)), 3.0, 3.0) { p.add(sh); }
    let mut mesh = IonMesh::new();
    if sc.show_ions {
        for ion in &sc.snap.ions {
            let pt = map(ion.p[0], ion.p[1]);
            if !rect.contains(pt) { continue; }
            mesh.dot(pt, 1.4, if ion.hot { pal.accent } else { with_alpha(species_color(pal.dark, ion.sp as usize), 200) });
        }
    }
    mesh.finish(&p);
    p.text(rect.left_top() + Vec2::new(6.0, 4.0), Align2::LEFT_TOP, "END VIEW x–y · z 9 mm", mono(8.5), pal.faint);
}

/// Frequency terrain strip. Returns Some(new f_AC kHz) while dragging.
pub fn draw_terrain(ui: &mut egui::Ui, rect: Rect, main_rect: Rect, sc: &Scene) -> Option<f64> {
    let pal = sc.pal;
    let f = Frame::new(main_rect, 13.0);
    let id = ui.id().with("terrain");
    let resp = ui.interact(rect, id, egui::Sense::click_and_drag());
    let p = ui.painter_at(rect);
    p.rect_filled(rect, 0.0, pal.bg);
    p.line_segment([rect.left_top(), rect.right_top()], Stroke::new(1.0, pal.line));
    // y range from the data (+ waterline)
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for s in 0..sc.species.len() {
        for k in 0..sc.terr.z_mm.len() {
            let (z, fx) = (sc.terr.z_mm[k], sc.terr.fx[s][k]);
            if z > 2.0 && z < 17.0 && fx.is_finite() { lo = lo.min(fx); hi = hi.max(fx); }
        }
    }
    lo = lo.min(sc.f_ac);
    hi = hi.max(sc.f_ac);
    if !lo.is_finite() || hi <= lo { lo = 20.0; hi = 120.0; }
    let pad = (hi - lo) * 0.18 + 2.0;
    let (flo, fhi) = (((lo - pad) / 5.0).floor() * 5.0, ((hi + pad) / 5.0).ceil() * 5.0);
    let (top, bot) = (rect.top() + 20.0, rect.bottom() - 20.0);
    let fy = |fk: f64| bot - ((fk - flo) / (fhi - flo)) as f32 * (bot - top);
    let f_of_y = |y: f32| flo + ((bot - y) / (bot - top)) as f64 * (fhi - flo);
    // grid
    let span = fhi - flo;
    let gstep = if span > 80.0 { 20.0 } else if span > 40.0 { 10.0 } else { 5.0 };
    let mut g = (flo / gstep).ceil() * gstep;
    while g <= fhi {
        let y = fy(g);
        p.line_segment([Pos2::new(rect.left(), y), Pos2::new(rect.right(), y)], Stroke::new(1.0, with_alpha(pal.line, 110)));
        p.text(Pos2::new(rect.left() + 8.0, y - 2.0), Align2::LEFT_BOTTOM, format!("{:.0}", g), mono(9.5), pal.dim);
        g += gstep;
    }
    let mut z = 0.0;
    while z <= 20.01 { p.line_segment([Pos2::new(f.zx(z), rect.top()), Pos2::new(f.zx(z), rect.bottom())], Stroke::new(1.0, with_alpha(pal.line, 70))); z += 2.0; }
    p.text(Pos2::new(rect.left() + 40.0, rect.top() + 5.0), Align2::LEFT_TOP, "FREQUENCY TERRAIN · secular f(z) [kHz] · solid f_x, dotted f_y · drag the waterline to set f_AC", mono(9.0), pal.faint);
    // curves
    for (s, _) in sc.species.iter().enumerate() {
        let col = species_color(pal.dark, s);
        for (ser, w, a) in [(&sc.terr.fy[s], 1.0f32, 90u8), (&sc.terr.fx[s], 1.6f32, 235u8)] {
            let mut run: Vec<Pos2> = vec![];
            let flush = |run: &mut Vec<Pos2>| {
                if run.len() > 1 {
                    if w < 1.2 { for sh in Shape::dashed_line(run, Stroke::new(w, with_alpha(col, a)), 2.0, 3.0) { p.add(sh); } }
                    else { p.add(Shape::line(std::mem::take(run), Stroke::new(w, with_alpha(col, a)))); }
                }
                run.clear();
            };
            for k in 0..sc.terr.z_mm.len() {
                let fv = ser[k];
                if fv.is_finite() { run.push(Pos2::new(f.zx(sc.terr.z_mm[k] as f32), fy(fv.clamp(flo - 5.0, fhi + 5.0)))); } else { flush(&mut run); }
            }
            flush(&mut run);
        }
    }
    // waterline
    let hovered = resp.hovered() || resp.dragged();
    let yw = fy(sc.f_ac);
    p.line_segment([Pos2::new(rect.left(), yw), Pos2::new(rect.right(), yw)], Stroke::new(if hovered { 1.8 } else { 1.2 }, sc.pal.accent));
    p.text(Pos2::new(rect.right() - 8.0, yw - 3.0), Align2::RIGHT_BOTTOM, format!("f_AC {:.2} kHz", sc.f_ac), mono(10.5), pal.accent);
    // crossings
    for (s, _) in sc.species.iter().enumerate() {
        let col = species_color(pal.dark, s);
        for zc in Terrain::crossings(&sc.terr.z_mm, &sc.terr.fx[s], sc.f_ac) {
            let c = Pos2::new(f.zx(zc as f32), yw);
            p.circle_filled(c, 4.5, pal.bg);
            p.circle_stroke(c, 4.5, Stroke::new(1.6, col));
        }
    }
    // gate / rod end markers
    for (z, t) in [(sc.viz.layers[0], "rods in"), (sc.viz.layers[1], "rods out")] {
        let x = f.zx(z as f32);
        p.line_segment([Pos2::new(x, bot), Pos2::new(x, bot + 4.0)], Stroke::new(1.0, pal.faint));
        p.text(Pos2::new(x, rect.bottom() - 3.0), Align2::CENTER_BOTTOM, t, mono(8.5), pal.faint);
    }
    let mut out = None;
    if resp.dragged() || resp.clicked() {
        if let Some(pp) = resp.interact_pointer_pos() { out = Some(f_of_y(pp.y).clamp(flo.max(5.0), fhi)); }
    }
    if hovered { ui.ctx().set_cursor_icon(egui::CursorIcon::ResizeVertical); }
    out
}

pub fn spark(p: &Painter, rect: Rect, data: &[f32], max: f32, col: Color32) {
    if data.len() < 2 || max <= 0.0 { return; }
    let n = data.len();
    let pts: Vec<Pos2> = data.iter().enumerate().map(|(k, &y)| Pos2::new(rect.left() + rect.width() * k as f32 / (n - 1) as f32, rect.bottom() - rect.height() * (y / max).clamp(0.0, 1.0))).collect();
    p.add(Shape::line(pts, Stroke::new(1.1, col)));
}

/// Sparkline with 1-sigma whiskers (every `every`-th sample) on the same scale.
pub fn spark_err(p: &Painter, rect: Rect, data: &[f32], err: &[f32], max: f32, col: Color32, every: usize) {
    if data.len() < 2 || max <= 0.0 { return; }
    let n = data.len();
    let x = |k: usize| rect.left() + rect.width() * k as f32 / (n - 1) as f32;
    let y = |v: f32| rect.bottom() - rect.height() * (v / max).clamp(0.0, 1.0);
    let wc = with_alpha(col, 90);
    let mut k = every / 2;
    while k < n {
        let e = err.get(k).copied().unwrap_or(0.0);
        let (y0, y1) = (y(data[k] - e), y(data[k] + e));
        p.line_segment([Pos2::new(x(k), y0), Pos2::new(x(k), y1)], Stroke::new(1.0, wc));
        k += every.max(1);
    }
    spark(p, rect, data, max, col);
}

pub fn _unused(_: &Painter) { let _ = sans(10.0); }
