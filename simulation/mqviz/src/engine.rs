//! The simulation thread: continuous injection, recirculation, rolling currents, resonance sweeps, validation.
//! It owns the `World`, steps blocks of BLOCK steps with rayon over ions, and publishes `Snapshot`s (Arc swap) to the UI.

use crate::json::J;
use crate::presets::{preset, PresetId};
use crate::stats::Frac;
use crate::validate::{self, Counts};
use crate::world::*;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, VecDeque};
use std::sync::atomic::{AtomicBool, AtomicU32, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

pub const MAXSP: usize = 8;
pub const MAXION: usize = 12_000;
const BIN_STEPS: u64 = 512; // rolling-rate bin = 20.48 us
const TRAIL_EVERY: u64 = 8; // blocks between trail samples (41 us)
const WARM_S: f64 = 3.5e-3; // simulated seconds of pre-roll so the beam is already flowing
const SNAP_EVERY: Duration = Duration::from_millis(16);

#[derive(Clone, PartialEq)]
pub struct Params {
    pub preset: PresetId,
    pub species: Vec<Sp>,
    pub v_rf: f64,
    pub v_ac: f64,
    pub f_ac_khz: f64,
    pub gate_v: f64,
    pub pressure_pa: f64,
    /// Fresh ions injected per simulated second.
    pub inj_rate: f64,
    /// Mean trap residence in wall-clock seconds as seen on screen; converted to simulated time with the realtime factor.
    pub residence_s: f64,
    pub recirculate: bool,
    /// Maximum simulated seconds per wall second (INFINITY = unthrottled).
    pub speed_cap: f64,
    pub paused: bool,
    /// Averaging window of the currents [simulated ms].
    pub window_ms: f64,
    /// Seed of the simulation RNG (injection times, species draw, recirculation delays; per-ion streams derive from it).
    pub seed: u64,
}

/// Default simulation seed (the constant the live engine has always used).
pub const DEFAULT_SEED: u64 = 0x1234_5678_9ABC_DEF1;

impl Params {
    pub fn from_preset(id: PresetId) -> Params {
        let p = preset(id);
        Params {
            preset: id, species: p.species, v_rf: p.v_rf, v_ac: p.v_ac, f_ac_khz: p.f_ac_khz, gate_v: p.gate_v, pressure_pa: p.pressure_pa,
            inj_rate: 3.0e6, residence_s: 0.2, recirculate: true, speed_cap: 6.0e-4, paused: false, window_ms: 1.0, seed: DEFAULT_SEED,
        }
    }
}

#[derive(Clone, Debug)]
pub struct SweepCfg { pub f_lo_khz: f64, pub f_hi_khz: f64, pub step_khz: f64, pub ions_per_point: usize, pub sweeps: usize }

pub enum Cmd { Reset, Step, SweepStart(SweepCfg), SweepStop, Validate { n: usize, n_port: usize } }

#[derive(Clone)]
pub struct IonVis { pub p: [f32; 3], pub trail: [[f32; 3]; NTRAIL], pub sp: u8, pub hot: bool }

#[derive(Clone, Default)]
pub struct SweepView {
    pub f: Vec<f64>,
    /// [point][species] = transmitted ions out of finished ions (single pass); value and error from `stats::Frac`.
    pub t: Vec<Vec<Frac>>,
    /// The configuration that produced it (record / command line).
    pub cfg: Option<SweepCfg>,
    pub tag: u64,
    pub cur: usize,
    pub sweeps_done: usize,
    pub sweeps_target: usize,
    pub running: bool,
    pub frac: f64,
}

#[derive(Clone)]
pub struct ValResult { pub n: usize, pub live: Counts, pub live0: Counts, pub port_n: usize, pub port: Counts, pub port0: Counts, pub live_s: f64, pub port_s: f64, pub live_mstep_s: f64 }

#[derive(Clone, Default)]
pub struct Snapshot {
    pub seq: u64,
    pub t_sim: f64,
    pub rtf: f64,
    pub steps_per_s: f64,
    pub n_flight: usize,
    pub warm: f64,
    pub ions: Vec<IonVis>,
    pub species: Vec<Sp>,
    /// Per species: transmitted, crashed, returned [ions/s] over the averaging window.
    pub rates: Vec<[f64; 3]>,
    pub cum: Vec<[u64; 3]>,
    pub fresh_rate: f64,
    pub reentry_rate: f64,
    pub f_ac_eff_khz: f64,
    pub sweep: Option<SweepView>,
    pub validating: bool,
    pub rebuilding_dc: bool,
    pub dropped: u64,
    pub timeouts: u64,
    pub gate_applied: f64,
    /// Per species: raw counts behind every displayed number (epoch totals, two consecutive halves, rolling windows).
    pub stats: Vec<SpStats>,
    /// Duration of one averaging window [simulated s] (the rates and `SpStats::win*` refer to it).
    pub win_dur: f64,
    /// Simulated time the fractions have been counted over [ms]: counting starts SETTLE_S after the last parameter change / reset.
    pub epoch_ms: f64,
    pub epoch_start_ms: f64,
}

/// Outcome counts of one species: [transmitted, crashed, returned, timeout].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SpStats {
    /// Ions born after the epoch start, all finished outcomes.
    pub ep: [u64; 4],
    /// First and second half of the epoch (equal duration), for the convergence test.
    pub half_a: [u64; 4],
    pub half_b: [u64; 4],
    /// Rolling window [T, C, R] counts, the window before it, and the one before that.
    pub win: [u64; 3],
    pub win_prev: [u64; 3],
}

impl SpStats {
    /// Per-pass transmitted fraction: T out of every finished ion.
    pub fn t_pass(c: &[u64; 4]) -> Frac { Frac::new(c[0], c.iter().sum()) }
    /// Steady-state efficiency T / (T + C).
    pub fn eff(c: &[u64; 4]) -> Frac { Frac::new(c[0], c[0] + c[1]) }
}

pub struct Shared {
    pub params: Mutex<Params>,
    pub cmds: Mutex<Vec<Cmd>>,
    pub snap: Mutex<Arc<Snapshot>>,
    pub val_progress: AtomicU32,
    pub val_result: Mutex<Option<ValResult>>,
    pub quit: AtomicBool,
}

impl Shared {
    pub fn new(p: Params) -> Arc<Shared> {
        Arc::new(Shared {
            params: Mutex::new(p), cmds: Mutex::new(vec![]), snap: Mutex::new(Arc::new(Snapshot::default())),
            val_progress: AtomicU32::new(0), val_result: Mutex::new(None), quit: AtomicBool::new(false),
        })
    }
    pub fn send(&self, c: Cmd) { self.cmds.lock().unwrap().push(c); }
    pub fn snapshot(&self) -> Arc<Snapshot> { self.snap.lock().unwrap().clone() }
}

// ------------------------------------------------------------------ rolling tallies

#[derive(Clone, Copy, Default)]
struct Bin { idx: u64, c: [[u32; 3]; MAXSP], fresh: u32, reentry: u32 }

struct Tally {
    bins: VecDeque<Bin>,
    cum: [[u64; 3]; MAXSP],
    timeouts: u64,
    dropped: u64,
}
impl Tally {
    fn new() -> Tally { Tally { bins: VecDeque::new(), cum: [[0; 3]; MAXSP], timeouts: 0, dropped: 0 } }
    fn bin(&mut self, step: u64) -> &mut Bin {
        let idx = step / BIN_STEPS;
        if self.bins.back().map_or(true, |b| b.idx != idx) {
            self.bins.push_back(Bin { idx, ..Default::default() });
            while self.bins.len() > 800 { self.bins.pop_front(); }
        }
        self.bins.back_mut().unwrap()
    }
    fn record(&mut self, step: u64, sp: usize, o: usize) { self.bin(step).c[sp][o] += 1; self.cum[sp][o] += 1; }
    fn window(&self, now_step: u64, w_steps: u64, nsp: usize) -> (Vec<[f64; 3]>, f64, f64) {
        let now_idx = now_step / BIN_STEPS;
        let nb = (w_steps / BIN_STEPS).max(1);
        let elapsed_bins = (now_idx + 1).min(nb) as f64;
        let dur = elapsed_bins * BIN_STEPS as f64 * DT;
        let mut r = vec![[0.0; 3]; nsp];
        let (mut fr, mut re) = (0.0, 0.0);
        for b in self.bins.iter().rev() {
            if b.idx + nb <= now_idx { break; }
            for s in 0..nsp { for o in 0..3 { r[s][o] += b.c[s][o] as f64; } }
            fr += b.fresh as f64;
            re += b.reentry as f64;
        }
        for s in 0..nsp { for o in 0..3 { r[s][o] /= dur; } }
        (r, fr / dur, re / dur)
    }
    /// Raw [T, C, R] counts of the window `back` windows before the current one (0 = current), same length as `window`.
    fn window_counts(&self, now_step: u64, w_steps: u64, nsp: usize, back: u64) -> Vec<[u64; 3]> {
        let now_idx = now_step / BIN_STEPS;
        let nb = (w_steps / BIN_STEPS).max(1);
        let hi = now_idx.saturating_sub(back * nb); // newest bin of this window
        let mut r = vec![[0u64; 3]; nsp];
        for b in self.bins.iter().rev() {
            if b.idx > hi { continue; }
            if b.idx + nb <= hi { break; }
            for s in 0..nsp { for o in 0..3 { r[s][o] += b.c[s][o] as u64; } }
        }
        r
    }
}

/// Outcome counts of ions born after the epoch start, kept in equal-width time batches so that the first and the second
/// half of the epoch can be compared (convergence test). Batches merge pairwise when there are more than `MAXB`.
struct Epoch {
    /// Step of the last parameter change: ions born before it are never counted.
    start_step: u64,
    /// Counting begins `SETTLE_STEPS` later: every ion that finishes after that was born after `start_step`
    /// (lifetimes are capped at T_MAX), so the fractions carry no right-censoring bias from the changeover.
    count0: u64,
    bw: u64,
    batches: Vec<[[u32; 4]; MAXSP]>,
    total: [[u64; 4]; MAXSP],
}
const MAXB: usize = 64;
/// Settling time after a parameter change before the fractions start counting [simulated s] = the longest ion lifetime.
pub const SETTLE_S: f64 = T_MAX;
const SETTLE_STEPS: u64 = (SETTLE_S / DT) as u64;
impl Epoch {
    fn new(start_step: u64) -> Epoch { Epoch { start_step, count0: start_step + SETTLE_STEPS, bw: 4096, batches: vec![], total: [[0; 4]; MAXSP] } }
    fn record(&mut self, event_step: u64, birth_step: u64, sp: usize, o: usize) {
        if birth_step < self.start_step || event_step < self.count0 { return; }
        let mut i = ((event_step - self.count0) / self.bw) as usize;
        while i >= MAXB {
            // merge pairs, double the width
            let mut m: Vec<[[u32; 4]; MAXSP]> = Vec::with_capacity(self.batches.len() / 2 + 1);
            for ch in self.batches.chunks(2) {
                let mut b = ch[0];
                if ch.len() > 1 { for s in 0..MAXSP { for k in 0..4 { b[s][k] += ch[1][s][k]; } } }
                m.push(b);
            }
            self.batches = m;
            self.bw *= 2;
            i = ((event_step - self.count0) / self.bw) as usize;
        }
        if self.batches.len() <= i { self.batches.resize(i + 1, [[0; 4]; MAXSP]); }
        self.batches[i][sp][o] += 1;
        self.total[sp][o] += 1;
    }
    /// (first half, second half) counts of species `sp` over equal durations: only complete batches are used (the newest is still
    /// filling), and with an odd number of them the middle one is left out.
    fn halves(&self, sp: usize) -> ([u64; 4], [u64; 4]) {
        let m = self.batches.len().saturating_sub(1);
        let (mut a, mut b) = ([0u64; 4], [0u64; 4]);
        if m < 2 { return (a, b); }
        let h = m / 2;
        for i in 0..h { for k in 0..4 { a[k] += self.batches[i][sp][k] as u64; } }
        for i in (m - h)..m { for k in 0..4 { b[k] += self.batches[i][sp][k] as u64; } }
        (a, b)
    }
}

// ------------------------------------------------------------------ sweep

struct SweepRun {
    cfg: SweepCfg,
    f: Vec<f64>,
    point: usize,
    inj_in_point: usize,
    acc: Vec<Vec<[u32; 4]>>, // [point][sp] trans, crash, ret, timeout
    sweeps_done: usize,
    draining: bool,
    running: bool,
}

// ------------------------------------------------------------------ simulation state

struct Sim {
    w: World,
    pool: rayon::ThreadPool,
    p: Params,
    sps: Vec<SpPhys>,
    cdf: Vec<f64>,
    cfg: StepCfg,
    clock: Clock,
    tabs: Tabs,
    ions: Vec<Ion>,
    step: u64,
    block_n: u64,
    rng: u64,
    next_id: u32,
    returns: BinaryHeap<Reverse<(u64, u8)>>,
    tally: Tally,
    epoch: Epoch,
    phys_key: [u64; 6],
    rtf_fixed: Option<f64>,
    sweep_id: u64,
    sweep: Option<SweepRun>,
    // pacing / metrics
    rtf: f64,
    sps_ema: f64,
    t_snap: Instant,
    step_snap: u64,
    ion_steps_snap: u64,
    ion_steps_acc: u64,
    seq: u64,
    anchor_wall: Instant,
    anchor_step: u64,
    // gate debounce
    gate_seen: f64,
    gate_changed: Instant,
    rebuilding: bool,
    warm_until: u64,
    validating: bool,
}

/// Everything that changes the per-pass physics: the epoch (fraction statistics) restarts when it changes.
fn phys_key_of(p: &Params) -> [u64; 6] {
    let sp = p.species.iter().fold(7u64, |h, s| h.wrapping_mul(1_000_003) ^ s.mass.to_bits() ^ s.z.to_bits().rotate_left(7) ^ s.ccs.to_bits().rotate_left(13) ^ s.weight.to_bits().rotate_left(29));
    [p.v_rf.to_bits(), p.v_ac.to_bits(), p.f_ac_khz.to_bits(), p.gate_v.to_bits(), p.pressure_pa.to_bits(), sp]
}

fn poisson(rng: &mut u64, lam: f64) -> usize {
    let mut u = || { *rng ^= *rng << 13; *rng ^= *rng >> 7; *rng ^= *rng << 17; (*rng >> 11) as f64 / (1u64 << 53) as f64 };
    if lam < 30.0 {
        let l = (-lam).exp();
        let (mut k, mut p) = (0usize, 1.0);
        loop { p *= u(); if p <= l { return k; } k += 1; }
    } else {
        let g = (-2.0 * u().max(1e-300).ln()).sqrt() * (TWO_PI * u()).cos();
        (lam + lam.sqrt() * g).round().max(0.0) as usize
    }
}

impl Sim {
    fn new(w: World, p: Params) -> Sim { Self::with_threads(w, p, default_threads()) }
    fn with_threads(w: World, p: Params, threads: usize) -> Sim {
        let ntab = 1;
        let mut s = Sim {
            pool: make_pool(threads), sps: vec![], cdf: vec![], cfg: StepCfg::new(p.v_rf, p.v_ac, p.pressure_pa),
            clock: Clock::new(ntab, p.f_ac_khz * 1e3), tabs: Tabs::new(), ions: Vec::with_capacity(MAXION), step: 0, block_n: 0, rng: p.seed | 1,
            next_id: 0, returns: BinaryHeap::new(), tally: Tally::new(), epoch: Epoch::new(0), phys_key: phys_key_of(&p), rtf_fixed: None, sweep_id: 0, sweep: None, rtf: 0.002, sps_ema: 0.0,
            t_snap: Instant::now(), step_snap: 0, ion_steps_snap: 0, ion_steps_acc: 0, seq: 0, anchor_wall: Instant::now(), anchor_step: 0,
            gate_seen: p.gate_v, gate_changed: Instant::now(), rebuilding: false, warm_until: (WARM_S / DT) as u64, validating: false, p: p.clone(), w,
        };
        s.set_species(&p.species);
        s
    }

    fn set_species(&mut self, sp: &[Sp]) {
        self.sps = sp.iter().map(SpPhys::new).collect();
        let tot: f64 = sp.iter().map(|s| s.weight).sum::<f64>().max(1e-12);
        let mut acc = 0.0;
        self.cdf = sp.iter().map(|s| { acc += s.weight / tot; acc }).collect();
    }

    fn reset(&mut self) {
        self.ions.clear();
        self.returns.clear();
        self.tally = Tally::new();
        self.step = 0;
        self.block_n = 0;
        self.clock = Clock::new(1 + self.sweep.as_ref().map_or(0, |s| s.f.len()), self.p.f_ac_khz * 1e3);
        if let Some(s) = &mut self.sweep { for (i, f) in s.f.iter().enumerate() { self.clock.f_ac[i + 1] = *f * 1e3; } }
        self.warm_until = (WARM_S / DT) as u64;
        self.anchor_wall = Instant::now();
        self.anchor_step = 0;
        self.epoch = Epoch::new(0);
    }

    fn pick_species(&mut self) -> u8 {
        self.rng ^= self.rng << 13; self.rng ^= self.rng >> 7; self.rng ^= self.rng << 17;
        let u = (self.rng >> 11) as f64 / (1u64 << 53) as f64;
        self.cdf.iter().position(|&c| u <= c).unwrap_or(self.cdf.len() - 1) as u8
    }

    fn spawn_ion(&mut self, sp: u8, reentry: bool) {
        if self.ions.len() >= MAXION { self.tally.dropped += 1; return; }
        let id = self.next_id;
        self.next_id = self.next_id.wrapping_add(1);
        let tab = match &self.sweep { Some(s) if s.running && !s.draining => (s.point + 1) as u8, _ => 0 };
        let seed = splitmix(self.rng ^ (id as u64).wrapping_mul(0x9E3779B97F4A7C15) ^ self.step);
        let ion = Ion::inject(seed, self.w.a, sp, &self.sps[sp as usize], id, tab, false);
        self.ions.push(ion);
        let b = self.tally.bin(self.step);
        if reentry { b.reentry += 1 } else { b.fresh += 1 }
    }

    fn apply_params(&mut self, np: &Params) {
        if np.seed != self.p.seed {
            self.rng = np.seed | 1;
            self.p.seed = np.seed;
            self.sweep = None;
            self.reset();
        }
        if np.species != self.p.species {
            self.set_species(&np.species);
            self.p.species = np.species.clone();
            self.p.preset = np.preset;
            self.sweep = None;
            self.reset();
        }
        self.cfg = StepCfg::new(np.v_rf, np.v_ac, np.pressure_pa);
        if !self.sweep.as_ref().map_or(false, |s| s.running) { self.clock.f_ac[0] = np.f_ac_khz * 1e3; }
        if np.paused != self.p.paused || np.speed_cap != self.p.speed_cap { self.anchor_wall = Instant::now(); self.anchor_step = self.step; }
        // gate: debounce the 0.1 s field recombination until the slider settles
        if np.gate_v != self.gate_seen { self.gate_seen = np.gate_v; self.gate_changed = Instant::now(); }
        if self.w.dc_volts[1] != np.gate_v && self.gate_changed.elapsed() > Duration::from_millis(140) {
            let mut v = DC_CAL;
            v[1] = np.gate_v;
            self.rebuilding = true;
            self.w.set_dc(v);
            self.rebuilding = false;
            self.anchor_wall = Instant::now();
            self.anchor_step = self.step;
        }
        self.p = np.clone();
        // the fraction statistics restart whenever the per-pass physics changes (gate: the field actually applied)
        let mut kp = np.clone();
        kp.gate_v = self.w.dc_volts[1];
        let key = phys_key_of(&kp);
        if key != self.phys_key { self.phys_key = key; self.epoch = Epoch::new(self.step); }
    }

    fn start_sweep(&mut self, c: SweepCfg) {
        self.sweep_id += 1;
        let mut f = vec![];
        let n = ((c.f_hi_khz - c.f_lo_khz) / c.step_khz).round() as usize;
        for i in 0..=n.min(100) { f.push(c.f_lo_khz + i as f64 * c.step_khz); }
        let nsp = self.sps.len();
        self.clock = Clock { rf: self.clock.rf, ac: { let mut a = self.clock.ac.clone(); a.resize(1 + f.len(), 0.0); a }, f_ac: { let mut a = vec![self.p.f_ac_khz * 1e3]; a.extend(f.iter().map(|x| x * 1e3)); a } };
        self.clock.f_ac[0] = f[0] * 1e3;
        self.sweep = Some(SweepRun { cfg: c, acc: vec![vec![[0; 4]; nsp]; f.len()], f, point: 0, inj_in_point: 0, sweeps_done: 0, draining: false, running: true });
    }

    fn stop_sweep(&mut self) {
        if let Some(s) = &mut self.sweep { s.draining = true; s.running = true; }
    }

    fn sweep_view(&self) -> Option<SweepView> {
        let s = self.sweep.as_ref()?;
        let nsp = self.sps.len();
        let t: Vec<Vec<Frac>> = s.acc.iter().map(|row| (0..nsp).map(|k| {
            let c = row[k];
            Frac::new(c[0] as u64, (c[0] + c[1] + c[2] + c[3]) as u64)
        }).collect()).collect();
        let tot = (s.cfg.sweeps * s.f.len()).max(1) as f64;
        Some(SweepView { f: s.f.clone(), t, cfg: Some(s.cfg.clone()), tag: self.sweep_id, cur: s.point, sweeps_done: s.sweeps_done, sweeps_target: s.cfg.sweeps, running: s.running,
                         frac: ((s.sweeps_done * s.f.len() + s.point) as f64 / tot).min(1.0) })
    }

    /// One block: advance all ions, account outcomes, inject.
    fn block(&mut self) {
        self.clock.fill(&mut self.tabs);
        advance_all(&self.pool, &mut self.ions, &self.w, &self.sps, &self.cfg, &self.tabs);
        let step0 = self.step;
        let step_end = step0 + BLOCK as u64;
        let mut finished: Vec<(u64, u8, u8, u8, u64)> = vec![]; // (event step, sp, outcome state, tab, birth step)
        let mut steps_done = 0u64;
        for ion in self.ions.iter_mut() {
            steps_done += ion.steps_blk as u64;
            ion.steps_blk = 0;
            if ion.state != ST_ALIVE { let ev = step0 + ion.fin_step as u64; finished.push((ev, ion.sp, ion.state, ion.tab, ev.saturating_sub(ion.age as u64))); }
        }
        self.ion_steps_acc += steps_done;
        if !finished.is_empty() { self.ions.retain(|i| i.state == ST_ALIVE); }
        for (st, sp, state, tab, birth) in finished {
            let spi = sp as usize;
            if tab == 0 { self.epoch.record(st, birth, spi, match state { ST_TRANS => 0, ST_CRASH => 1, ST_RETURN => 2, _ => 3 }); }
            match state {
                ST_TRANS => self.tally.record(st, spi, 0),
                ST_CRASH => self.tally.record(st, spi, 1),
                ST_RETURN => {
                    self.tally.record(st, spi, 2);
                    if self.p.recirculate {
                        self.rng ^= self.rng << 13; self.rng ^= self.rng >> 7; self.rng ^= self.rng << 17;
                        let u = ((self.rng >> 11) as f64 / (1u64 << 53) as f64).max(1e-12);
                        let mean_sim = self.p.residence_s * self.rtf_fixed.unwrap_or(self.rtf).max(1e-5);
                        let delay = (-u.ln() * mean_sim / DT).ceil() as u64;
                        self.returns.push(Reverse((st + delay, sp)));
                    }
                }
                _ => self.tally.timeouts += 1,
            }
            if tab > 0 {
                if let Some(s) = &mut self.sweep {
                    let o = match state { ST_TRANS => 0, ST_CRASH => 1, ST_RETURN => 2, _ => 3 };
                    if let Some(row) = s.acc.get_mut(tab as usize - 1) { row[spi][o] += 1; }
                }
            }
        }
        // trails
        self.block_n += 1;
        if self.block_n % TRAIL_EVERY == 0 {
            for ion in self.ions.iter_mut() {
                for k in (1..NTRAIL).rev() { ion.trail[k] = ion.trail[k - 1]; }
                ion.trail[0] = [ion.x[0] as f32 * 1e3, ion.x[1] as f32 * 1e3, ion.x[2] as f32 * 1e3];
            }
        }
        self.step = step_end;
        // re-entries
        while let Some(&Reverse((due, sp))) = self.returns.peek() {
            if due > step_end { break; }
            self.returns.pop();
            self.spawn_ion(sp, true);
        }
        // fresh injection (stops while a sweep drains)
        let draining = self.sweep.as_ref().map_or(false, |s| s.draining);
        if !draining {
            let lam = self.p.inj_rate * BLOCK as f64 * DT;
            let n = poisson(&mut self.rng, lam);
            for _ in 0..n {
                let sp = self.pick_species();
                self.spawn_ion(sp, false);
                if let Some(s) = &mut self.sweep {
                    if s.running {
                        s.inj_in_point += 1;
                        if s.inj_in_point >= s.cfg.ions_per_point {
                            s.inj_in_point = 0;
                            s.point += 1;
                            if s.point >= s.f.len() {
                                s.sweeps_done += 1;
                                s.point = 0;
                                if s.sweeps_done >= s.cfg.sweeps { s.draining = true; }
                            }
                            let f = s.f[s.point] * 1e3;
                            self.clock.f_ac[0] = f;
                        }
                    }
                }
            }
        }
        // sweep finished when the last tagged ion is gone
        if let Some(s) = &mut self.sweep {
            if s.draining && s.running && !self.ions.iter().any(|i| i.tab > 0) {
                s.running = false;
                s.draining = false; // the live beam resumes (previously fresh injection stayed off after a sweep until reset)
                s.point = s.f.len().saturating_sub(1);
                self.clock.f_ac[0] = self.p.f_ac_khz * 1e3;
                self.epoch = Epoch::new(step_end);
            }
        }
    }

    fn snapshot(&mut self) -> Snapshot {
        let now = Instant::now();
        let dtw = now.duration_since(self.t_snap).as_secs_f64().max(1e-6);
        let rtf_i = (self.step - self.step_snap) as f64 * DT / dtw;
        let sps_i = self.ion_steps_acc as f64 / dtw;
        if !self.p.paused || self.step != self.step_snap {
            self.rtf = if self.rtf_valid() { 0.7 * self.rtf + 0.3 * rtf_i } else { rtf_i };
            self.sps_ema = 0.7 * self.sps_ema + 0.3 * sps_i;
        }
        self.t_snap = now;
        self.step_snap = self.step;
        self.ion_steps_snap += self.ion_steps_acc;
        self.ion_steps_acc = 0;
        let a = self.w.a;
        let ions: Vec<IonVis> = self.ions.iter().map(|i| {
            let tr = |p: [f32; 3]| [p[0] - (a * 1e3) as f32, p[1] - (a * 1e3) as f32, p[2]];
            IonVis { p: [(i.x[0] - a) as f32 * 1e3, (i.x[1] - a) as f32 * 1e3, i.x[2] as f32 * 1e3],
                     trail: [tr(i.trail[0]), tr(i.trail[1]), tr(i.trail[2]), tr(i.trail[3])], sp: i.sp, hot: i.r_mm(a) > R_EXCITED * 1e3 && i.x[2] > 3.5e-3 }
        }).collect();
        let nsp = self.sps.len();
        let w_steps = (self.p.window_ms * 1e-3 / DT) as u64;
        let (rates, fr, re) = self.tally.window(self.step, w_steps, nsp);
        let (w0, w1) = (self.tally.window_counts(self.step, w_steps, nsp, 0), self.tally.window_counts(self.step, w_steps, nsp, 1));
        let stats: Vec<SpStats> = (0..nsp).map(|s| {
            let (a, b) = self.epoch.halves(s);
            SpStats { ep: self.epoch.total[s], half_a: a, half_b: b, win: w0[s], win_prev: w1[s] }
        }).collect();
        self.seq += 1;
        let warm = (self.step as f64 / self.warm_until.max(1) as f64).min(1.0);
        Snapshot {
            seq: self.seq, t_sim: self.step as f64 * DT, rtf: if self.p.paused { 0.0 } else { self.rtf }, steps_per_s: if self.p.paused { 0.0 } else { self.sps_ema },
            n_flight: self.ions.len(), warm, ions, species: self.p.species.clone(), rates,
            cum: (0..nsp).map(|s| self.tally.cum[s]).collect(), fresh_rate: fr, reentry_rate: re,
            f_ac_eff_khz: self.clock.f_ac[0] / 1e3, sweep: self.sweep_view(), validating: self.validating, rebuilding_dc: self.rebuilding,
            dropped: self.tally.dropped, timeouts: self.tally.timeouts, gate_applied: self.w.dc_volts[1],
            stats, win_dur: (((w_steps / BIN_STEPS).max(1)) * BIN_STEPS) as f64 * DT, epoch_ms: (self.step.saturating_sub(self.epoch.count0)) as f64 * DT * 1e3,
            epoch_start_ms: self.epoch.start_step as f64 * DT * 1e3,
        }
    }
    fn rtf_valid(&self) -> bool { self.sps_ema > 0.0 }
}

// Scheduling only: the desktop frontend bounds work without changing integration or counts.
fn desktop_idle(active: bool, sweeping: bool, step: bool) -> bool { !active && !sweeping && !step }
fn desktop_rest(work: Duration, sweeping: bool) -> Duration { work.mul_f64(if sweeping { 1.0 } else { 3.0 }).min(Duration::from_millis(100)) }
fn run(w: World, sh: Arc<Shared>, ctx: Option<egui::Context>, desktop: Option<(Arc<AtomicBool>, Arc<AtomicBool>)>) {
    let p0 = sh.params.lock().unwrap().clone();
    let mut sim = Sim::with_threads(w, p0, if desktop.is_some() { 2 } else { default_threads() });
    let mut step_req = false;
    let mut last_pub = Instant::now();
    while !sh.quit.load(Ordering::Relaxed) {
        let cycle = Instant::now();
        let np = sh.params.lock().unwrap().clone();
        sim.apply_params(&np);
        // commands (after params, so a preset change and a sweep start sent together stay in order)
        let cmds: Vec<Cmd> = std::mem::take(&mut *sh.cmds.lock().unwrap());
        for c in cmds {
            match c {
                Cmd::Reset => { sim.sweep = None; sim.reset(); }
                Cmd::Step => step_req = true,
                Cmd::SweepStart(cfg) => { sim.start_sweep(cfg); }
                Cmd::SweepStop => sim.stop_sweep(),
                Cmd::Validate { n, n_port } => {
                    sim.validating = true;
                    *sh.snap.lock().unwrap() = Arc::new(sim.snapshot());
                    let was = sim.w.dc_volts;
                    if was != DC_CAL { sim.w.set_dc(DC_CAL); }
                    let bak = preset(PresetId::S14).species[0].clone();
                    let t = Instant::now();
                    let fac = validate::CAL_FAC_KHZ * 1e3;
                    let live = validate::run_live(&sim.w, &sim.pool, &bak, n, validate::CAL_VRF, validate::CAL_VAC, fac, validate::CAL_P, 0xC0FFEE, false, Some(&sh.val_progress));
                    let live0 = validate::run_live(&sim.w, &sim.pool, &bak, n, validate::CAL_VRF, validate::CAL_VAC, fac, validate::CAL_P, 0xBEEF, true, None);
                    let live_s = t.elapsed().as_secs_f64() / 2.0;
                    let t = Instant::now();
                    let port = validate::run_reference(&sim.w, &bak, n_port, validate::CAL_VRF, validate::CAL_VAC, fac, validate::CAL_P, 777, false);
                    let port0 = validate::run_reference(&sim.w, &bak, n_port, validate::CAL_VRF, validate::CAL_VAC, fac, validate::CAL_P, 778, true);
                    let port_s = t.elapsed().as_secs_f64() / 2.0;
                    if was != DC_CAL { sim.w.set_dc(was); }
                    *sh.val_result.lock().unwrap() = Some(ValResult { n, live, live0, port_n: n_port, port, port0, live_s, port_s, live_mstep_s: live.steps as f64 / live_s / 1e6 });
                    sim.validating = false;
                    sim.anchor_wall = Instant::now();
                    sim.anchor_step = sim.step;
                }
            }
        }
        let sweeping = sim.sweep.as_ref().map_or(false, |s| s.running);
        let studying = sweeping || desktop.as_ref().is_some_and(|(_,study)| study.load(Ordering::Relaxed));
        if desktop.as_ref().is_some_and(|(active,_)| desktop_idle(active.load(Ordering::Relaxed), studying, step_req)) {
            // Exclude background time from pacing/RTF; do not advance or discard ions.
            sim.anchor_wall = Instant::now(); sim.anchor_step = sim.step;
            sim.t_snap = Instant::now(); sim.step_snap = sim.step; sim.ion_steps_acc = 0;
            std::thread::sleep(Duration::from_millis(100));
            continue;
        }
        let warm = sim.step < sim.warm_until;
        if sim.p.paused && !step_req && !warm {
            std::thread::sleep(Duration::from_millis(if desktop.is_some() { 50 } else { 8 }));
            if last_pub.elapsed() > Duration::from_millis(50) {
                *sh.snap.lock().unwrap() = Arc::new(sim.snapshot());
                last_pub = Instant::now();
                if let Some(c) = &ctx { c.request_repaint(); }
            }
            continue;
        }
        let nblocks = if step_req { 8 } else { 1 };
        for _ in 0..nblocks { sim.block(); }
        step_req = false;
        if desktop.is_some() { std::thread::sleep(desktop_rest(cycle.elapsed(), studying)); }
        // pacing
        let sweeping = sim.sweep.as_ref().map_or(false, |s| s.running);
        if !warm && !sweeping && sim.p.speed_cap.is_finite() {
            let sim_elapsed = (sim.step - sim.anchor_step) as f64 * DT;
            let want = sim_elapsed / sim.p.speed_cap;
            let have = sim.anchor_wall.elapsed().as_secs_f64();
            if want > have { std::thread::sleep(Duration::from_secs_f64((want - have).min(0.02))); }
            // never accumulate a deficit: if the machine cannot keep up, run at its own pace
            if have > want + 0.1 { sim.anchor_wall = Instant::now(); sim.anchor_step = sim.step; }
        } else {
            sim.anchor_wall = Instant::now();
            sim.anchor_step = sim.step;
        }
        if last_pub.elapsed() >= if desktop.is_some() { Duration::from_millis(100) } else { SNAP_EVERY } {
            *sh.snap.lock().unwrap() = Arc::new(sim.snapshot());
            last_pub = Instant::now();
            if let Some(c) = &ctx { c.request_repaint(); }
        }
    }
}

pub fn spawn(w: World, sh: Arc<Shared>, ctx: Option<egui::Context>) -> std::thread::JoinHandle<()> {
    std::thread::Builder::new().name("mqsim-main".into()).spawn(move || run(w, sh, ctx, None)).unwrap()
}

/// Bounded desktop work: two workers, 25% duty live / 50% during explicit sweeps.
pub fn spawn_desktop(w: World, sh: Arc<Shared>, active: Arc<AtomicBool>, study: Arc<AtomicBool>) -> std::thread::JoinHandle<()> {
    std::thread::Builder::new().name("mqsim-main".into()).spawn(move || run(w, sh, None, Some((active, study)))).unwrap()
}

/// Headless throughput benchmark of the live engine: `mqviz --bench [seconds]`.
pub fn bench(secs: f64, inj_rate: f64) {
    let w = World::build_calibrated();
    let mut p = Params::from_preset(PresetId::S14);
    p.inj_rate = inj_rate;
    let mut sim = Sim::new(w, p);
    let t0 = Instant::now();
    let mut samples: Vec<(f64, usize, f64, f64)> = vec![];
    let mut last = Instant::now();
    let (mut s_last, mut steps_last) = (0u64, 0u64);
    while t0.elapsed().as_secs_f64() < secs {
        sim.block();
        if last.elapsed() >= Duration::from_millis(500) {
            let dtw = last.elapsed().as_secs_f64();
            let rtf = (sim.step - s_last) as f64 * DT / dtw;
            let sps = (sim.ion_steps_acc - steps_last.min(sim.ion_steps_acc)) as f64 / dtw;
            samples.push((t0.elapsed().as_secs_f64(), sim.ions.len(), rtf, sps));
            println!("  t {:5.1} s  sim {:7.3} ms  in flight {:5}  RTF {:.5}  {:6.1} M ion-steps/s", t0.elapsed().as_secs_f64(), sim.step as f64 * DT * 1e3, sim.ions.len(), rtf, sps / 1e6);
            last = Instant::now();
            s_last = sim.step;
            steps_last = 0;
            sim.ion_steps_acc = 0;
        }
    }
    let (tr, cr, rt): (u64, u64, u64) = (0..sim.sps.len()).fold((0, 0, 0), |a, s| (a.0 + sim.tally.cum[s][0], a.1 + sim.tally.cum[s][1], a.2 + sim.tally.cum[s][2]));
    println!("cumulative: transmitted {tr}, crashed {cr}, returned {rt}, timeouts {}, dropped {}", sim.tally.timeouts, sim.tally.dropped);
}

/// Headless resonance scan using the same multiplexed-tag machinery as the GUI sweep: `mqviz --scan [s14|fig5] [ions_per_point]`.
/// (Kept as it was; `--sweep` is the full-featured path with fits and a run record.)
pub fn scan_cli(args: &[String]) {
    let which = args.iter().skip_while(|a| *a != "--scan").nth(1).cloned().unwrap_or_else(|| "s14".into());
    let n: usize = args.iter().skip_while(|a| *a != "--scan").nth(2).and_then(|s| s.parse().ok()).unwrap_or(2000);
    let id = if which == "fig5" { PresetId::Fig5 } else { PresetId::S14 };
    let pr = preset(id);
    let w = World::build_calibrated();
    let mut p = Params::from_preset(id);
    p.inj_rate = 8.0e6;
    let mut sim = Sim::new(w, p);
    sim.warm_until = 0;
    sim.start_sweep(SweepCfg { f_lo_khz: pr.sweep.0, f_hi_khz: pr.sweep.1, step_khz: pr.sweep.2, ions_per_point: n, sweeps: 1 });
    let t = Instant::now();
    while sim.sweep.as_ref().map_or(false, |s| s.running) { sim.block(); }
    println!("scan {} · RF {} V · AC ±{} V · gate {} V · {} Pa · {} ions/point (tagged cohorts) · {:.1} s wall", pr.name, pr.v_rf, pr.v_ac, pr.gate_v, pr.pressure_pa, n, t.elapsed().as_secs_f64());
    let sv = sim.sweep_view().unwrap();
    print!("  f[kHz] ");
    for s in &sim.p.species { print!(" {:>9}", s.label); }
    println!();
    for (f, row) in sv.f.iter().zip(&sv.t) {
        print!("  {:6.2} ", f);
        for c in row { print!(" {:>8.2}%", 100.0 * c.p()); }
        println!();
    }
    for (si, s) in sim.p.species.iter().enumerate() {
        let (xs, ys): (Vec<f64>, Vec<f64>) = sv.f.iter().zip(&sv.t).filter_map(|(f, r)| if r[si].n > 0 { Some((*f, r[si].p())) } else { None }).unzip();
        if xs.len() > 3 { let (fp, yp, w) = chsim::channel::fwhm(&xs, &ys); println!("  {}: peak {:.1} % at {:.2} kHz, FWHM {:.2} kHz, f/Δf {:.1}", s.label, 100.0 * yp, fp, w, fp / w.max(1e-9)); }
    }
}

fn die(msg: String) -> ! { eprintln!("error: {msg}"); std::process::exit(2) }

/// World with the DC set the parameters ask for (the GUI debounces this; headless applies it directly).
fn headless_world(p: &Params) -> World {
    let mut w = World::build_calibrated();
    let mut dc = DC_CAL;
    dc[1] = p.gate_v;
    w.set_dc(dc);
    w
}

fn fmt_pct(v: crate::stats::Val) -> String { format!("{:6.2} ± {:4.2} %", 100.0 * v.v, 100.0 * v.e) }

/// `mqviz --sweep [s14|fig5|custom] [--rf V] [--vac V] [--gate V] [--pressure Pa] [--n ions/pt] [--flo kHz] [--fhi kHz] [--step kHz]
///  [--sweeps K] [--seed S] [--inj-rate r] [--boot B] [--out file.json]`: resonance sweep with per-point binomial errors, Lorentzian/Gaussian
/// fits (R in f/Δf and m/z), and a run record (JSON + CSV). Deterministic given the seed.
pub fn sweep_cli(args: &[String]) {
    use crate::record::{self, arg_f, arg_val};
    let mut p = record::params_from_args(args, "--sweep").unwrap_or_else(|e| die(e));
    if arg_f(args, "--inj-rate").is_none() { p.inj_rate = 8.0e6; }
    let pr = preset(p.preset);
    let cfg = SweepCfg {
        f_lo_khz: arg_f(args, "--flo").unwrap_or(pr.sweep.0), f_hi_khz: arg_f(args, "--fhi").unwrap_or(pr.sweep.1), step_khz: arg_f(args, "--step").unwrap_or(pr.sweep.2),
        ions_per_point: arg_f(args, "--n").map_or(2000, |x| x as usize), sweeps: arg_f(args, "--sweeps").map_or(1, |x| x as usize),
    };
    let boot = arg_f(args, "--boot").map_or(200, |x| x as usize);
    let w = headless_world(&p);
    let env = record::RunEnv::new(format!("{:?}", w.geom));
    let mut sim = Sim::new(w, p.clone());
    sim.warm_until = 0;
    sim.start_sweep(cfg.clone());
    let t = Instant::now();
    while sim.sweep.as_ref().map_or(false, |s| s.running) { sim.block(); }
    let wall = t.elapsed().as_secs_f64();
    let sv = sim.sweep_view().unwrap();
    let dc = { let mut d = DC_CAL; d[1] = p.gate_v; d };
    let fits = crate::analysis::fit_sweep(&sv, &p.species, &sim.w.terr, p.v_rf, dc, boot, p.seed);
    println!("sweep {} · RF {} V · AC ±{} V · gate {} V · {} Pa · {} ions/point × {} sweep(s) (tagged cohorts, seed {}) · {:.1} s wall · git {}",
             pr.name, p.v_rf, p.v_ac, p.gate_v, p.pressure_pa, cfg.ions_per_point, cfg.sweeps, record::seed_str(p.seed), wall, env.git.short());
    print!("  f[kHz] ");
    for s in &p.species { print!(" {:>18}", s.label); }
    println!();
    for (f, row) in sv.f.iter().zip(&sv.t) {
        print!("  {:6.2} ", f);
        for c in row { print!(" {:>18}", if c.n > 0 { fmt_pct((*c).into()) } else { "—".into() }); }
        println!();
    }
    for f in fits.iter().flatten() {
        println!("  {}: {} fit  f0 = {:.2} ± {:.2} kHz  FWHM = {:.2} ± {:.2} kHz  peak {:.1} ± {:.1} %  χ²/dof {:.2} ({} pts)", f.label, f.fit.shape.name(), f.fit.p[2], f.fit.e[2], f.fit.p[3], f.fit.e[3], 100.0 * f.peak.v, 100.0 * f.peak.e, f.fit.chi2_red(), f.fit.n_pts);
        let boot_txt = f.fit.r_boot.map_or(String::new(), |(m, s)| format!(" (bootstrap {:.2} ± {:.2})", m, s));
        println!("      R = f/Δf = {:.2} ± {:.2}{}   |   df/d(m/z) = {:.4} kHz per m/z → FWHM = {:.2} ± {:.2} m/z, R = (m/z)/Δ(m/z) = {:.2} ± {:.2}",
                 f.r_f.v, f.r_f.e, boot_txt, f.dfdmz, f.fwhm_mz.v, f.fwhm_mz.e, f.r_mz.v, f.r_mz.e);
    }
    if p.preset == PresetId::S14 {
        if let (Some(f), Some(dir)) = (fits.get(0).and_then(|f| f.as_ref()), crate::reality::find_data_dir()) {
            if let Ok(d) = crate::reality::RealityData::load(&dir) {
                let r = &crate::reality::s14_rows(&d, &p.species, &[], Some(&f.fit))[0];
                println!("  reality: {}  model {}  vs theirs {}  pull {:+.2}σ {}", r.name, r.ours_txt, r.theirs_txt, r.pull, r.mark());
            }
        }
    }
    let cmd = record::cmd_sweep(&p, &cfg, p.inj_rate);
    let mut j = record::sweep_record(&env, &p, &sv, &fits, &cmd);
    if let J::Obj(o) = &mut j { o.push(("wall_seconds".into(), wall.into())); }
    let out = arg_val(args, "--out").map(std::path::PathBuf::from).unwrap_or_else(|| record::new_path("sweep", "json"));
    if let Some(d) = out.parent() { let _ = std::fs::create_dir_all(d); }
    record::write_json(&out, &j).unwrap_or_else(|e| die(format!("{}: {e}", out.display())));
    let csv = out.with_extension("csv");
    std::fs::write(&csv, record::sweep_csv(&sv, &p.species)).unwrap_or_else(|e| die(format!("{}: {e}", csv.display())));
    println!("  record: {}\n          {}\n  reproduce: {}", out.display(), csv.display(), cmd);
}

/// `mqviz --live [s14|fig5|custom] [--fac kHz] [--rf V] [--vac V] [--gate V] [--pressure Pa] [--sim-ms T] [--warm-ms W] [--rtf x] [--seed S]
///  [--inj-rate r] [--out file.json]`: run the live engine headless (fixed realtime factor for the recirculation delay) and report every observable
/// with its error and convergence, plus the reality rows against their published numbers.
pub fn live_cli(args: &[String]) {
    use crate::record::{self, arg_f, arg_val};
    use crate::reality;
    let p = record::params_from_args(args, "--live").unwrap_or_else(|e| die(e));
    let sim_ms = arg_f(args, "--sim-ms").unwrap_or(40.0);
    let warm_ms = arg_f(args, "--warm-ms").unwrap_or(WARM_S * 1e3);
    let rtf = arg_f(args, "--rtf").unwrap_or(1.5e-3);
    let w = headless_world(&p);
    let env = record::RunEnv::new(format!("{:?}", w.geom));
    let mut sim = Sim::new(w, p.clone());
    sim.rtf_fixed = Some(rtf);
    sim.rtf = rtf;
    sim.warm_until = (warm_ms * 1e-3 / DT) as u64;
    let t = Instant::now();
    while sim.step < sim.warm_until { sim.block(); }
    sim.tally = Tally::new();
    let end = sim.warm_until + (sim_ms * 1e-3 / DT) as u64;
    while sim.step < end { sim.block(); }
    let wall = t.elapsed().as_secs_f64();
    let snap = sim.snapshot();
    let pr = preset(p.preset);
    println!("live {} · RF {} V · AC ±{} V @ {} kHz · gate {} V · {} Pa · {} ms after {} ms warm-up · rtf {:.2e} · seed {} · {:.1} s wall · git {}",
             pr.name, p.v_rf, p.v_ac, p.f_ac_khz, p.gate_v, p.pressure_pa, sim_ms, warm_ms, rtf, record::seed_str(p.seed), wall, env.git.short());
    println!("  {:<14} {:>9} {:>9} {:>9} {:>8}   {:>20}  {:>20}  convergence", "species", "T", "C", "R", "timeout", "T per pass", "T/(T+C)");
    for (s, st) in snap.species.iter().zip(&snap.stats) {
        let (tp, ef) = (crate::analysis::t_pass(st), crate::analysis::efficiency_sp(st));
        println!("  {:<14} {:>9} {:>9} {:>9} {:>8}   {:>20}  {:>20}  {} / {}", s.label, st.ep[0], st.ep[1], st.ep[2], st.ep[3],
                 if tp.val.ok() { fmt_pct(tp.val) } else { "—".into() }, if ef.val.ok() { fmt_pct(ef.val) } else { "—".into() }, tp.settle.label(), ef.settle.label());
    }
    let eff = crate::analysis::efficiency(&snap.stats);
    println!("  mixture T/(T+C) = {}  ({})", fmt_pct(eff.val), eff.settle.label());
    let target = 0usize;
    // reality rows
    let mut rows = vec![];
    let mut spec_j = J::Null;
    let data = reality::find_data_dir().and_then(|d| reality::RealityData::load(&d).ok());
    if let Some(d) = &data {
        match p.preset {
            PresetId::S14 => rows = reality::s14_rows(d, &snap.species, &snap.stats, None),
            PresetId::Fig5 => {
                let (cz, cf, dist) = reality::nearest_channel(p.f_ac_khz, p.v_rf);
                println!("  channel: {}+ (nominal {:.2} kHz at RF {} V; f_AC is {:.2} kHz away)", cz, cf, p.v_rf, dist);
                rows = reality::fig5_rows(d, &snap.species, &snap.stats, cz);
                if let Some(dir) = reality::find_data_dir() {
                    match crate::spectrum_model::SpecModel::build(&dir) {
                        Ok(m) => {
                            let tl = crate::analysis::live_t_lookup(&snap.species, &snap.stats);
                            if let Some(pp) = m.predict(cz, &tl) {
                                let resid: Vec<f64> = pp.data.iter().zip(&pp.model).map(|(d, m)| d - m).collect();
                                let rms = (resid.iter().map(|r| r * r).sum::<f64>() / resid.len() as f64).sqrt();
                                println!("  predicted {} spectrum: -logL {:.0} over {} columns, calibration shift {:+.2} m/z, rms residual {:.2} % (model build {:.1} s); not simulated: {} species",
                                         pp.panel, pp.nll, pp.mz.len(), pp.shift, rms, m.build_secs, pp.missing.len());
                                spec_j = J::obj().with("panel", pp.panel.clone()).with("channel_z", pp.channel_z).with("shift_mz", pp.shift).with("neg_log_likelihood", pp.nll).with("mz", pp.mz.clone()).with("data_pct", pp.data.clone())
                                    .with("model_pct", pp.model.clone()).with("model_lo_pct", pp.lo.clone()).with("model_hi_pct", pp.hi.clone()).with("display_floor_pct", pp.floor)
                                    .with("not_simulated", pp.missing.clone()).with("adducts_with_parent_T", pp.proxied.clone());
                            }
                        }
                        Err(e) => eprintln!("  spectrum model: {e}"),
                    }
                }
            }
            PresetId::Custom => {}
        }
        for r in &rows {
            println!("  reality  {:<34} model {:>24}   theirs {:>20}   pull {:>6}  {}   {}", r.name, r.ours_txt, r.theirs_txt, if r.pull.is_finite() { format!("{:+.2}σ", r.pull) } else { "—".into() }, r.mark(), r.settle.map_or(String::new(), |s| s.label()));
        }
    } else { println!("  (their published data not found: set MQVIZ_DATA to data/their_side/figure_data)"); }
    let cmd = record::cmd_live(&p, sim_ms, rtf);
    let mut j = record::live_record(&env, &p, &snap, rtf, target, &rows, &cmd);
    if let J::Obj(o) = &mut j { o.push(("predicted_spectrum".into(), spec_j)); o.push(("wall_seconds".into(), wall.into())); }
    let out = arg_val(args, "--out").map(std::path::PathBuf::from).unwrap_or_else(|| record::new_path("live", "json"));
    if let Some(d) = out.parent() { let _ = std::fs::create_dir_all(d); }
    record::write_json(&out, &j).unwrap_or_else(|e| die(format!("{}: {e}", out.display())));
    println!("  record: {}\n  reproduce: {}", out.display(), cmd);
}

/// Headless benchmark through the same thread + snapshot path as the GUI: `mqviz --bench-thread [secs]`.
pub fn bench_thread(secs: f64) {
    let w = World::build_calibrated();
    let sh = Shared::new(Params::from_preset(PresetId::S14));
    let h = spawn(w, sh.clone(), None);
    let t0 = Instant::now();
    while t0.elapsed().as_secs_f64() < secs {
        std::thread::sleep(Duration::from_millis(500));
        let s = sh.snapshot();
        println!("  t {:5.1} s  warm {:.2}  sim {:7.3} ms  in flight {:5}  RTF {:.5}  {:6.1} M ion-steps/s", t0.elapsed().as_secs_f64(), s.warm, s.t_sim * 1e3, s.n_flight, s.rtf, s.steps_per_s / 1e6);
    }
    sh.quit.store(true, Ordering::Relaxed);
    let _ = h.join();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn epoch_counts_only_ions_born_after_the_change_and_after_the_settling_time() {
        let mut e = Epoch::new(1000);
        e.record(1000 + SETTLE_STEPS - 1, 1500, 0, 0); // finishes before counting starts
        e.record(1000 + SETTLE_STEPS + 10, 900, 0, 0); // born before the change
        assert_eq!(e.total[0][0], 0);
        e.record(1000 + SETTLE_STEPS + 10, 1200, 0, 0);
        e.record(1000 + SETTLE_STEPS + 5000, 1300, 1, 2);
        assert_eq!((e.total[0][0], e.total[1][2]), (1, 1));
    }

    #[test]
    fn epoch_batches_merge_without_losing_counts_and_halves_are_ordered() {
        let mut e = Epoch::new(0);
        let n = 200_000u64;
        for i in 0..n { e.record(SETTLE_STEPS + i * 97, 0, 0, (i % 3) as usize); }
        assert!(e.batches.len() <= MAXB);
        let tot: u64 = e.total[0].iter().sum();
        assert_eq!(tot, n);
        let in_batches: u64 = e.batches.iter().map(|b| b[0].iter().map(|&x| x as u64).sum::<u64>()).sum();
        assert_eq!(in_batches, n);
        let (a, b) = e.halves(0);
        let (sa, sb): (u64, u64) = (a.iter().sum(), b.iter().sum());
        assert!(sa + sb <= n && sa + sb > n * 9 / 10, "halves {sa} {sb}");
        // uniform arrival: the halves hold the same number of events within a few %
        assert!((sa as f64 / sb as f64 - 1.0).abs() < 0.05, "{sa} vs {sb}");
    }

    #[test]
    fn window_counts_match_the_window_rates() {
        let mut t = Tally::new();
        for i in 0..5000u64 { t.record(i * 40, 0, (i % 3) as usize); }
        let (step, w) = (5000 * 40, 20 * BIN_STEPS);
        let (rates, _, _) = t.window(step, w, 1);
        let c0 = t.window_counts(step, w, 1, 0);
        let dur = 20.0 * BIN_STEPS as f64 * DT;
        for o in 0..3 { assert!((rates[0][o] * dur - c0[0][o] as f64).abs() < 1e-6); }
        // consecutive windows tile the history: the sum of the first k windows back equals the cumulative count
        let all: u64 = (0..((step / BIN_STEPS) / 20 + 1)).map(|k| t.window_counts(step, w, 1, k)[0].iter().sum::<u64>()).sum();
        assert_eq!(all, 5000);
    }
}

#[cfg(test)]
mod desktop_scheduling_tests {
 use super::*;
 #[test] fn background_stops_live_work_but_not_requested_studies_or_steps() {
  assert!(desktop_idle(false,false,false));assert!(!desktop_idle(false,true,false));assert!(!desktop_idle(false,false,true));assert!(!desktop_idle(true,false,false));
 }
 #[test] fn work_budget_leaves_cpu_headroom_and_keeps_commands_responsive() {
  assert_eq!(desktop_rest(Duration::from_millis(4),false),Duration::from_millis(12));
  assert_eq!(desktop_rest(Duration::from_millis(4),true),Duration::from_millis(4));
  assert_eq!(desktop_rest(Duration::from_secs(1),false),Duration::from_millis(100));
 }
 #[test]
 #[ignore = "builds the calibrated world; run explicitly with the field cache"]
 fn desktop_background_lifecycle() {
  let shared=Shared::new(Params::from_preset(PresetId::S14));
  let active=Arc::new(AtomicBool::new(false));
  let study=Arc::new(AtomicBool::new(false));
  let thread=spawn_desktop(World::build_calibrated(),shared.clone(),active.clone(),study.clone());
  let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
   let wait_for=|predicate:&dyn Fn()->bool| {
    let deadline=Instant::now()+Duration::from_secs(30);
    while !predicate() { assert!(Instant::now()<deadline,"desktop engine did not progress");std::thread::sleep(Duration::from_millis(50)); }
   };
   std::thread::sleep(Duration::from_millis(250));
   assert_eq!(shared.snapshot().t_sim,0.);
   // Study preparation must advance while unfocused, before SweepStart exists.
   study.store(true,Ordering::Relaxed);
   wait_for(&||shared.snapshot().t_sim>0.);
   shared.send(Cmd::SweepStart(SweepCfg {f_lo_khz:83.,f_hi_khz:83.,step_khz:0.5,ions_per_point:500,sweeps:1}));
   wait_for(&||shared.snapshot().sweep.is_some());
   let before=shared.snapshot().t_sim;
   wait_for(&||shared.snapshot().t_sim>before);
   shared.send(Cmd::Reset);study.store(false,Ordering::Relaxed);
   std::thread::sleep(Duration::from_millis(300));
   let held=shared.snapshot().t_sim;
   std::thread::sleep(Duration::from_millis(300));
   assert_eq!(shared.snapshot().t_sim,held,"cancelled background work must stop");
   active.store(true,Ordering::Relaxed);
   wait_for(&||shared.snapshot().t_sim>held);
  }));
  shared.quit.store(true,Ordering::Relaxed);thread.join().unwrap();
  result.unwrap();
 }

}
