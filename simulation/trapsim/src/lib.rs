//! Dilute, domain-coupled tandem transport. SI units except the inherited cube geometry (mm).
pub mod fields;
pub mod geometry;
use fields::Fields;
use geometry::Geometry;
use mqviz::{
    presets::{preset, PresetId},
    world::*,
};
use serde::Serialize;
use std::f64::consts::TAU;

pub const DEFAULT_TONES: [f64; 9] = [39., 0., 47., 0., 0., 0., 53., 0., 59.25]; // reconstructed kHz
#[derive(Clone, Copy, Debug, PartialEq)]
enum Region {
    Split,
    Port(usize),
    Combine,
    Lost,
    Back,
    Detected,
}
struct Track {
    ion: Ion,
    region: Region,
    last: Option<usize>,
}
#[derive(Clone, Default, Serialize)]
pub struct Counts {
    pub entered: u32,
    pub returned: u32,
    pub forwarded: u32,
    pub revisits: u32,
    pub from_other: u32,
    pub retry_from: [u32; 9],
    pub reverse_entries: u32,
    pub lost: u32,
    pub by_charge: [u32; 5],
}
#[derive(Serialize)]
pub struct Frame {
    pub ms: f64,
    pub ports: [Counts; 9],
    pub split: usize,
    pub channel: usize,
    pub combine: usize,
    pub lost: usize,
    pub back: usize,
    pub detected: usize,
    pub points: Vec<[f64; 5]>,
    pub size_mm: f64,
    pub seam_mm: f64,
}
pub struct Model {
    cube: Fields,
    dc: [Vec<[f32; 4]>; 2],
    channel: World,
    mouth: [[f64; 2]; 9],
    size: f64,
    seam: f64,
    cfg: StepCfg,
    species: Vec<SpPhys>,
    tracks: Vec<Track>,
    counts: [Counts; 9],
    step: u64,
    tones: [f64; 9],
}
/// Shared by the live cube and the reference-field replay (SI coordinates).
#[inline]
pub fn cube_field(cube: &Fields, dc: &[[f32; 4]], x: [f64; 3], rf: f64) -> [f64; 3] {
    let p = x.map(|v| v * 1000.);
    let (c, w) = cube.grid_coord(p);
    let mut e = [0.; 3];
    for k in 0..2 {
        for j in 0..2 {
            for i in 0..2 {
                let q = cube.idx(c[0] + i, c[1] + j, c[2] + k);
                let wt = if i == 1 { w[0] } else { 1. - w[0] }
                    * if j == 1 { w[1] } else { 1. - w[1] }
                    * if k == 1 { w[2] } else { 1. - w[2] };
                for d in 0..3 {
                    e[d] += wt * (rf * cube.rf[q][d] as f64 + dc[q][d] as f64);
                }
            }
        }
    }
    e
}
#[inline]
pub fn transport_step(
    x: &mut [f64; 3],
    v: &mut [f64; 3],
    rng: &mut u64,
    sp: &SpPhys,
    cfg: &StepCfg,
    e: [f64; 3],
) {
    for d in 0..3 {
        v[d] += sp.qm * cfg.dt * e[d];
        x[d] += v[d] * cfg.dt;
    }
    collide(v, rng, sp, cfg, cfg.ngas * sp.sig * cfg.dt);
}

impl Model {
    pub fn new(tones: [f64; 9], n: usize, seed: u64, fine: bool, mouth_start: bool) -> Self {
        let mut g = Geometry::as_built();
        g.array = true;
        let (g, h) = g.snap(if fine { 0.52 } else { 1.04 });
        let mouth = std::array::from_fn(|i| {
            let (_, _, x, y) = g.port_axis_xy(&g.port(i));
            [x * 0.001, y * 0.001]
        });
        let size = g.size() * 0.001;
        let mut cube = Fields::build(&g, h, 0.4, 1e-7);
        let dc = [
            cube.dc_combined(20., 20.1, 20.5, 18.),
            cube.dc_combined(15., 15.1, 14.8, 15.),
        ];
        // Sampling only needs the RF grid and the two combined DC grids.
        cube.dc = Default::default();
        cube.id = vec![];
        cube.sign = vec![];
        cube.phi_rf = vec![];
        let channel = World::build_calibrated();
        let seam = size - (g.post_top_depth() + channel.geom.q1_start - Z_RETURN * 1000.) * 0.001;
        let preset = preset(PresetId::Fig5);
        let cfg = StepCfg::new(preset.v_rf, preset.v_ac, preset.pressure_pa);
        let species: Vec<_> = preset.species.iter().map(SpPhys::new).collect();
        let total: f64 = preset.species.iter().map(|s| s.weight).sum();
        let tracks = (0..n)
            .map(|i| {
                let key = splitmix(seed.wrapping_add(i as u64));
                let mut r = (key >> 11) as f64 / (1u64 << 53) as f64 * total;
                let sp = preset
                    .species
                    .iter()
                    .position(|s| {
                        r -= s.weight;
                        r <= 0.
                    })
                    .unwrap_or(4);
                let mut ion =
                    Ion::inject(splitmix(key), 0., sp as u8, &species[sp], i as u32, 0, true);
                let [x, y, _] = ion.x;
                if mouth_start {
                    let c = mouth[i % 9];
                    ion.x = [c[0] + x, c[1] + y, seam - 0.0001];
                } else {
                    let (_, _, b, c) = g.port_axis_xy(&g.entrance);
                    ion.x = [c * 0.001 + x * 0.1, size - 0.0001, b * 0.001 + y * 0.1];
                    ion.v = [ion.v[0], -ion.v[2], ion.v[1]];
                }
                Track {
                    ion,
                    region: Region::Split,
                    last: None,
                }
            })
            .collect();
        Self {
            cube,
            dc,
            channel,
            mouth,
            size,
            seam,
            cfg,
            species,
            tracks,
            counts: Default::default(),
            step: 0,
            tones,
        }
    }
    pub fn advance(&mut self, steps: usize) {
        for _ in 0..steps {
            let t = self.step as f64 * DT;
            let rf = self.cfg.vrf * (TAU * F_RF * t).sin();
            let ac = self
                .tones
                .map(|f| self.cfg.vac * (TAU * f * 1000. * t).sin());
            for tr in &mut self.tracks {
                let reg = tr.region;
                let ion = &mut tr.ion;
                let sp = &self.species[ion.sp as usize];
                let mut x = ion.x;
                let mut v = ion.v;
                let e = match reg {
                    Region::Split | Region::Combine => cube_field(
                        &self.cube,
                        &self.dc[(reg == Region::Combine) as usize],
                        x,
                        rf,
                    ),
                    Region::Port(i) => {
                        let (r, a, d) = self.channel.fields_at(x);
                        // Inactive ports have no excitation and a +4 V gate increment (assumed blocking bias).
                        let bias = if self.tones[i] == 0. {
                            self.channel.gate_delta(x, 4.)
                        } else {
                            [0.; 3]
                        };
                        std::array::from_fn(|k| rf * r[k] + ac[i] * a[k] + d[k] + bias[k])
                    }
                    _ => continue,
                };
                transport_step(&mut x, &mut v, &mut ion.rng, sp, &self.cfg, e);
                ion.x = x;
                ion.v = v;
                match reg {
                    Region::Split | Region::Combine => {
                        let mm = x.map(|v| v * 1000.);
                        if self.cube.geom.classify(mm, 0.).0 != 0 {
                            tr.region = Region::Lost;
                            continue;
                        }
                        if x[2] > self.seam {
                            let i = self.cube.geom.nearest_port(mm[0], mm[1]);
                            let c = self.mouth[i];
                            if (x[0] - c[0]).hypot(x[1] - c[1]) < self.cube.geom.ap_r * 0.001 {
                                enter(
                                    tr,
                                    i,
                                    reg,
                                    self.seam,
                                    c,
                                    self.channel.a,
                                    &mut self.counts[i],
                                );
                                continue;
                            }
                            if x[2] > self.size {
                                tr.region = Region::Lost;
                                continue;
                            }
                        }
                        if x[1] > self.size {
                            let (_, _, b, c) =
                                self.cube.geom.port_axis_xy(&self.cube.geom.entrance);
                            if (mm[0] - c).hypot(mm[2] - b) >= self.cube.geom.ap_r {
                                tr.region = Region::Lost;
                                continue;
                            }
                        }
                        if x[1] > self.size + self.cube.geom.l_ext * 0.001 {
                            tr.region = if reg == Region::Split {
                                Region::Back
                            } else {
                                Region::Detected
                            };
                        } else if x.iter().any(|p| *p < 0.) || x[0] > self.size {
                            tr.region = Region::Lost;
                        }
                    }
                    Region::Port(i) => {
                        if self.channel.is_crash(x) {
                            tr.region = Region::Lost;
                            self.counts[i].lost += 1;
                        } else if x[2] < Z_RETURN || x[2] > Z_EXIT {
                            leave(
                                tr,
                                i,
                                self.seam,
                                self.mouth[i],
                                self.channel.a,
                                &mut self.counts[i],
                            );
                        }
                    }
                    _ => unreachable!(),
                }
            }
            self.step += 1;
        }
    }
    pub fn frame(&self) -> Frame {
        let mut f = Frame {
            ms: self.step as f64 * DT * 1000.,
            ports: self.counts.clone(),
            split: 0,
            channel: 0,
            combine: 0,
            lost: 0,
            back: 0,
            detected: 0,
            points: vec![],
            size_mm: self.size * 1000.,
            seam_mm: self.seam * 1000.,
        };
        for tr in &self.tracks {
            let [x, y, z] = tr.ion.x;
            let (x, z) = match tr.region {
                Region::Lost => {
                    f.lost += 1;
                    continue;
                }
                Region::Back => {
                    f.back += 1;
                    continue;
                }
                Region::Detected => {
                    f.detected += 1;
                    continue;
                }
                Region::Split => {
                    f.split += 1;
                    (x, z)
                }
                Region::Combine => {
                    f.combine += 1;
                    (x, 2. * self.seam + Z_EXIT - Z_RETURN - z)
                }
                Region::Port(i) => {
                    f.channel += 1;
                    (
                        self.mouth[i][0] + x - self.channel.a,
                        self.seam + z - Z_RETURN,
                    )
                }
            };
            let y = if let Region::Port(i) = tr.region {
                self.mouth[i][1] + y - self.channel.a
            } else {
                y
            };
            f.points.push([
                x * 1000.,
                z * 1000.,
                tr.ion.sp as f64,
                tr.ion.id as f64,
                y * 1000.,
            ]);
        }
        f
    }
}
fn enter(
    tr: &mut Track,
    i: usize,
    from: Region,
    size: f64,
    c: [f64; 2],
    axis: f64,
    count: &mut Counts,
) {
    let x = &mut tr.ion.x;
    x[0] += axis - c[0];
    x[1] += axis - c[1];
    let depth = x[2] - size;
    if from == Region::Split {
        x[2] = Z_RETURN + depth;
        count.entered += 1;
        if let Some(last) = tr.last {
            count.revisits += 1;
            count.retry_from[last] += 1;
            if last != i {
                count.from_other += 1;
            }
        }
    } else {
        x[2] = Z_EXIT - depth;
        tr.ion.v[2] *= -1.;
        count.reverse_entries += 1;
    }
    tr.region = Region::Port(i);
}
fn leave(tr: &mut Track, i: usize, size: f64, c: [f64; 2], axis: f64, count: &mut Counts) {
    let x = &mut tr.ion.x;
    x[0] += c[0] - axis;
    x[1] += c[1] - axis;
    if x[2] < Z_RETURN {
        x[2] = size + x[2] - Z_RETURN;
        tr.region = Region::Split;
        tr.last = Some(i);
        count.returned += 1;
    } else {
        x[2] = size - (x[2] - Z_EXIT);
        tr.ion.v[2] *= -1.;
        tr.region = Region::Combine;
        tr.last = None;
        count.forwarded += 1;
        count.by_charge[tr.ion.sp as usize] += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn nine_open_mouths_with_solid_board_between() {
        let mut g = Geometry::as_built();
        g.array = true;
        for i in 0..9 {
            let (_, _, x, y) = g.port_axis_xy(&g.port(i));
            assert_eq!(g.nearest_port(x, y), i);
            assert_eq!(g.classify([x, y, g.size() - 0.2], 0.).0, 0);
        }
        let (_, _, x, y) = g.port_axis_xy(&g.port(0));
        assert_eq!(
            g.classify([x + g.pitch, y, g.size() - 0.2], 0.).0,
            geometry::ID_PLANE
        );
    }
    #[test]
    fn rejection_preserves_phase_space_and_allows_a_different_port() {
        let phys = SpPhys::new(&mqviz::presets::ubq(7, 1.));
        let mut tr = Track {
            ion: Ion::inject(111, 0., 0, &phys, 11, 0, true),
            region: Region::Split,
            last: None,
        };
        let mut counts: [Counts; 9] = Default::default();
        let size = 0.086;
        let c = [0.024, 0.042];
        let axis = 0.008;
        let pos = [0.0242, 0.0419, size + 1e-6];
        tr.ion.x = pos;
        let rng = tr.ion.rng;
        enter(&mut tr, 0, Region::Split, size, c, axis, &mut counts[0]);
        tr.ion.v = [3., 7., -20.];
        tr.ion.x[2] = Z_RETURN - 1e-6;
        leave(&mut tr, 0, size, c, axis, &mut counts[0]);
        assert_eq!(tr.region, Region::Split);
        assert!((tr.ion.x[0] - pos[0]).abs() < 1e-14);
        assert!((tr.ion.x[2] - (size - 1e-6)).abs() < 1e-14);
        assert_eq!(tr.ion.v, [3., 7., -20.]);
        assert_eq!(tr.ion.rng, rng);
        tr.ion.x = [0.04, 0.042, size + 1e-6];
        enter(
            &mut tr,
            1,
            Region::Split,
            size,
            [0.04, 0.042],
            axis,
            &mut counts[1],
        );
        assert_eq!(counts[1].from_other, 1);
        tr.ion.x[2] = Z_EXIT + 1e-6;
        tr.ion.v[2] = 30.;
        leave(&mut tr, 1, size, [0.04, 0.042], axis, &mut counts[1]);
        assert_eq!(tr.region, Region::Combine);
        assert_eq!(tr.ion.v[2], -30.);
    }
}

#[cfg(test)]
mod transport_regression {
    use super::*;
    // Explicit opt-in: uses the real channel cache and solves both trap field settings.
    #[test]
    #[ignore]
    fn physical_routing() {
        let fine = std::env::var_os("MULTIQ_FINE").is_some();
        let inlet = std::env::var_os("MULTIQ_INLET").is_some();
        let mut m = Model::new(DEFAULT_TONES, 180, 486, fine, !inlet);
        // A track inside an active port must advance like Explore, including
        // collisions and RNG consumption. Isolate a block before testing routing.
        let ion = Ion::inject(42, m.channel.a, 0, &m.species[0], 0, 0, true);
        let mut old = ion.clone();
        old.as_ = 0.;
        old.ac = 1.;
        let mut tabs = Tabs::new();
        Clock::new(1, DEFAULT_TONES[0] * 1000.).fill(&mut tabs);
        let tracks = std::mem::replace(
            &mut m.tracks,
            vec![Track {
                ion: old.clone(),
                region: Region::Port(0),
                last: None,
            }],
        );
        advance(&mut old, &m.channel, &m.species[0], &m.cfg, &tabs);
        m.advance(BLOCK);
        assert_eq!(old.state, ST_ALIVE);
        assert_eq!(m.tracks[0].region, Region::Port(0));
        for j in 0..3 {
            assert!((old.x[j] - m.tracks[0].ion.x[j]).abs() < 1e-12);
            assert!((old.v[j] - m.tracks[0].ion.v[j]).abs() < 1e-7);
        }
        assert_eq!(old.rng, m.tracks[0].ion.rng);
        m.tracks = tracks;
        m.step = 0;
        m.advance(250_000);
        let f = m.frame();
        let sum = |fun: fn(&Counts) -> u32| f.ports.iter().map(fun).sum::<u32>();
        println!("inlet={inlet} fine={fine} entries={} returned={} other={} forwarded={} reverse={} collected={} retained={} lost={} source_back={}",sum(|p|p.entered),sum(|p|p.returned),sum(|p|p.from_other),sum(|p|p.forwarded),sum(|p|p.reverse_entries),f.detected,f.split+f.combine+f.channel,f.lost,f.back);
        assert_eq!(
            f.split + f.combine + f.channel + f.detected + f.lost + f.back,
            180
        );
        assert_eq!(
            sum(|p| p.entered + p.reverse_entries),
            sum(|p| p.returned + p.forwarded + p.lost) + f.channel as u32
        );
        assert!(sum(|p| p.from_other) > 0);
        assert!(sum(|p| p.forwarded) > 0);
        if !inlet {
            assert!(f.ports.iter().all(|p| p.entered > 0));
            assert!(sum(|p| p.returned) > 0);
        }
    }
}
