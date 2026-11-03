//! Tiled / batched execution path of the channel kernel (`Path::Tiled`, the default). See channel_tiled.metal and PERF.md.
//! Per-ion arithmetic is identical to the direct reference kernel (`Engine::run_direct`); `chsim check-tiled` verifies bitwise equality.

use crate::channel::*;
use metal::*;
use std::ffi::c_void;
use std::rc::Rc;
use std::time::Instant;

fn tail_tg() -> usize { std::env::var("CHSIM_TAILTG").ok().and_then(|v| v.parse().ok()).unwrap_or(128) }

#[derive(Clone, Copy, Debug)]
pub struct TileCfg {
    pub tg: usize,             // threads (ions) per threadgroup of the tiled rounds (fixed at pipeline build)
    pub chunk: i32,            // max steps per ion per tiled round (CHSIM_CHUNK)
    pub rounds_per_cb: usize,  // tiled rounds encoded per command buffer (CHSIM_RPC)
    /// Run threadgroup-memory tiled rounds while at least this many ions are alive, then hand over to the worker kernel.
    /// 0 = never (default: measured no faster than the register-cached brick workers, see PERF.md). CHSIM_ROUNDS.
    pub rounds_min_live: usize,
    pub tail_chunk: i32,       // steps per ion in the first worker launch (doubles each launch up to tail_cap) CHSIM_TCHUNK
    pub tail_cap: i32,         // CHSIM_TCAP
}
impl TileCfg {
    pub fn from_env() -> TileCfg {
        let g = |k: &str, d: usize| std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(d);
        TileCfg { tg: g("CHSIM_TG", 256), chunk: g("CHSIM_CHUNK", 64) as i32, rounds_per_cb: g("CHSIM_RPC", 8),
                  rounds_min_live: g("CHSIM_ROUNDS", 0), tail_chunk: g("CHSIM_TCHUNK", 20000) as i32, tail_cap: g("CHSIM_TCAP", 200000) as i32 }
    }
}

pub struct TPipes {
    brick: ComputePipelineState, clear: ComputePipelineState, count: ComputePipelineState, scan: ComputePipelineState,
    scatter: ComputePipelineState, push: ComputePipelineState, tail: ComputePipelineState, filter: ComputePipelineState,
    zero_args: ComputePipelineState, size_args: ComputePipelineState, tail_leg: ComputePipelineState,
}

#[repr(C)]
#[derive(Clone, Copy, Debug, bytemuck::Pod, bytemuck::Zeroable)]
struct TParams { nbx: i32, nby: i32, nbz: i32, nblocks: i32, chunk: i32, n_ions: i32, pad: [i32; 2] }

impl Engine {
    fn tpipes(&self, nearest: bool) -> Rc<TPipes> {
        let k = nearest as usize;
        if let Some(p) = &self.tpipes.borrow()[k] { return p.clone(); }
        let t_compile = Instant::now();
        let tg = self.tile.borrow().tg;
        // developer hook: CHSIM_DEFS="NAME=VALUE,NAME2" injects #defines (e.g. TAIL_LSPAN=32, TAIL_WMAX=8192)
        let extra: String = std::env::var("CHSIM_DEFS").unwrap_or_default().split(',').filter(|d| !d.is_empty())
            .map(|d| format!("#define {}\n", d.replace('=', " "))).collect();
        let src = format!("{extra}#define TILE_TG {tg}\n#define TAIL_TG {}\n#define CH_NEAREST {}\n{}{}\n{}", tail_tg(), nearest as i32,
            if self.verbose { "#define CH_STATS 1\n" } else { "" }, include_str!("channel.metal"), include_str!("channel_tiled.metal"));
        let opts = CompileOptions::new();
        opts.set_fast_math_enabled(false); // same arithmetic as the direct kernel (bitwise parity)
        let lib = self.device.new_library_with_source(&src, &opts).expect("MSL compile failed (tiled)");
        let mk = |name: &str| self.device.new_compute_pipeline_state_with_function(&lib.get_function(name, None).unwrap()).unwrap();
        let p = Rc::new(TPipes { brick: mk("tc_brick"), clear: mk("tc_clear"), count: mk("tc_count"), scan: mk("tc_scan"), scatter: mk("tc_scatter"),
                                 push: mk("tc_push"), tail: mk("tc_tail"), filter: mk("tc_filter"), zero_args: mk("tc_zero_args"), size_args: mk("tc_size_args"), tail_leg: mk("tc_tail_leg") });
        if self.verbose { eprintln!("  tiled pipelines (nearest={}) compiled in {:.3} s; max threads/TG: tail {} tail_leg {} push {} (exec width {})", nearest, t_compile.elapsed().as_secs_f64(),
            p.tail.max_total_threads_per_threadgroup(), p.tail_leg.max_total_threads_per_threadgroup(), p.push.max_total_threads_per_threadgroup(), p.tail.thread_execution_width()); }
        self.tpipes.borrow_mut()[k] = Some(p.clone());
        p
    }

    /// Brick-layout field buffer (RF | combined DC | AC per node, 7-cell bricks) for this DC configuration; built on the GPU, last config cached.
    fn brick_buffer(&self, tp: &TPipes, dc: [f64; 4], dcb: &Buffer, pr: &ChParams) -> Buffer {
        let key = [dc[0].to_bits(), dc[1].to_bits(), dc[2].to_bits(), dc[3].to_bits(), pr.nearest as u64];
        if let Some((k, b)) = &*self.bricks.borrow() { if *k == key { return b.clone(); } }
        *self.bricks.borrow_mut() = None; // free the old layout before allocating the new one
        let (nx, ny, nz) = (pr.nx as usize, pr.ny as usize, pr.nz as usize);
        let nbricks = ((nx - 1 + 6) / 7) * ((ny - 1 + 6) / 7) * ((nz - 1 + 6) / 7);
        let t0 = Instant::now();
        let out = shared_buf(&self.device, nbricks * 512 * 48);
        let cmd = self.queue.new_command_buffer();
        let enc = cmd.new_compute_command_encoder();
        enc.set_compute_pipeline_state(&tp.brick);
        enc.set_buffer(1, Some(&self.rf), 0); enc.set_buffer(2, Some(dcb), 0); enc.set_buffer(3, Some(&self.ac), 0);
        enc.set_bytes(4, std::mem::size_of::<ChParams>() as u64, pr as *const ChParams as *const c_void);
        enc.set_buffer(22, Some(&out), 0);
        enc.dispatch_threads(MTLSize::new((nbricks * 512) as u64, 1, 1), MTLSize::new(256, 1, 1));
        enc.end_encoding(); cmd.commit(); cmd.wait_until_completed();
        if self.verbose { eprintln!("  brick layout built in {:.3} s", t0.elapsed().as_secs_f64()); }
        *self.bricks.borrow_mut() = Some((key, out.clone()));
        out
    }

    /// Advance every ion in `ib` until it has finished (status != 0). `dcb` is the combined DC field of configuration `dc`.
    pub(crate) fn run_tiled(&self, ib: &Buffer, jb: &Buffer, jpb: &Buffer, dcb: &Buffer, dc: [f64; 4], pr: &ChParams, n: usize) {
        let nearest = pr.nearest != 0;
        let cfg = *self.tile.borrow();
        let tp = self.tpipes(nearest);
        let fb = self.brick_buffer(&tp, dc, dcb, pr);
        let dev = &self.device;
        // dense live-ion lists (ping-pong) with their counters, and indirect-dispatch arguments
        let act = [shared_buf(dev, n * 4), shared_buf(dev, n * 4)]; let ctr = shared_buf(dev, 32); let args = shared_buf(dev, 64);
        unsafe {
            let a = act[0].contents() as *mut u32; for i in 0..n { *a.add(i) = i as u32; }
            let c = ctr.contents() as *mut u32; *c = n as u32; *c.add(2) = 0;
        }
        let mut round = 0usize;
        if cfg.rounds_min_live > 0 && n >= cfg.rounds_min_live {
            round = self.tiled_rounds(&tp, &cfg, ib, &fb, jb, jpb, pr, n, &act, &ctr, &args);
        }
        self.run_tail(&tp, &cfg, ib, &fb, jb, jpb, pr, &act, &ctr, &args, round, n);
    }

    /// Threadgroup-memory tiled rounds (bin ions by block, stage the block, advance, park). Returns the round counter.
    fn tiled_rounds(&self, tp: &TPipes, cfg: &TileCfg, ib: &Buffer, fb: &Buffer, jb: &Buffer, jpb: &Buffer, pr: &ChParams, n: usize,
                    act: &[Buffer; 2], ctr: &Buffer, args: &Buffer) -> usize {
        let (nx, ny, nz) = (pr.nx as usize, pr.ny as usize, pr.nz as usize);
        let (nbx, nby, nbz) = ((nx - 1 + 6) / 7, (ny - 1 + 6) / 7, (nz - 1 + 6) / 7);
        let nblocks = nbx * nby * nbz;
        let dev = &self.device;
        let tparams = TParams { nbx: nbx as i32, nby: nby as i32, nbz: nbz as i32, nblocks: nblocks as i32, chunk: cfg.chunk, n_ions: n as i32, pad: [0; 2] };
        let blk_of = shared_buf(dev, n * 4); let rank_of = shared_buf(dev, n * 4);
        let cnt = shared_buf(dev, nblocks * 4); let offs = shared_buf(dev, nblocks * 4);
        unsafe { std::ptr::write_bytes(cnt.contents() as *mut u8, 0, nblocks * 4); }
        let max_tg = n / cfg.tg + nblocks.min(n) + 2;
        let tgb = shared_buf(dev, max_tg * 4); let sorted = shared_buf(dev, max_tg * cfg.tg * 4);
        let stats_b = shared_buf(dev, 64);
        unsafe { std::ptr::write_bytes(stats_b.contents() as *mut u8, 0, 64); }
        let (mut visits, mut tsteps) = (0u64, 0u64);
        let psz = std::mem::size_of::<ChParams>() as u64; let tsz = std::mem::size_of::<TParams>() as u64;
        let pp = pr as *const ChParams as *const c_void; let tpp = &tparams as *const TParams as *const c_void;
        let t_start = Instant::now();
        let mut round = 0usize;
        loop {
            let cmd = self.queue.new_command_buffer();
            let enc = cmd.new_compute_command_encoder();
            enc.set_buffer(0, Some(ib), 0); enc.set_buffer(1, Some(fb), 0);
            enc.set_bytes(4, psz, pp); enc.set_buffer(5, Some(jpb), 0); enc.set_buffer(6, Some(jb), 0); enc.set_bytes(7, tsz, tpp);
            enc.set_buffer(8, Some(&cnt), 0); enc.set_buffer(9, Some(&blk_of), 0); enc.set_buffer(10, Some(&rank_of), 0);
            enc.set_buffer(11, Some(&offs), 0); enc.set_buffer(12, Some(&sorted), 0); enc.set_buffer(13, Some(&tgb), 0);
            enc.set_buffer(18, Some(args), 0); enc.set_buffer(19, Some(&stats_b), 0);
            let one = MTLSize::new(1, 1, 1); let w256 = MTLSize::new(256, 1, 1);
            for _ in 0..cfg.rounds_per_cb {
                let (ain, aout) = (&act[round % 2], &act[(round + 1) % 2]);
                let (c_in, c_out) = ((round % 2) as u64 * 8, ((round + 1) % 2) as u64 * 8);
                enc.set_buffer(14, Some(ain), 0); enc.set_buffer(15, Some(aout), 0);
                enc.set_buffer(16, Some(ctr), c_in); enc.set_buffer(17, Some(ctr), c_out);
                enc.set_compute_pipeline_state(&tp.clear);
                enc.dispatch_thread_groups(MTLSize::new(((nblocks + 255) / 256) as u64, 1, 1), w256);
                enc.set_compute_pipeline_state(&tp.count); enc.dispatch_thread_groups_indirect(args, 0, w256);
                enc.set_compute_pipeline_state(&tp.scan); enc.dispatch_thread_groups(one, MTLSize::new(1024, 1, 1));
                enc.set_compute_pipeline_state(&tp.scatter); enc.dispatch_thread_groups_indirect(args, 16, w256);
                enc.set_compute_pipeline_state(&tp.push); enc.dispatch_thread_groups_indirect(args, 32, MTLSize::new(cfg.tg as u64, 1, 1));
                round += 1;
            }
            enc.end_encoding(); cmd.commit(); cmd.wait_until_completed();
            let alive = unsafe { *(ctr.contents() as *const u32).add((round % 2) * 2) } as usize;
            unsafe { let sp = stats_b.contents() as *mut u32; visits += *sp.add(2) as u64; *sp.add(2) = 0; *sp.add(3) = 0; tsteps += *sp.add(4) as u64; *sp.add(4) = 0; }
            if alive == 0 || alive < cfg.rounds_min_live { break; }
        }
        if self.verbose { eprintln!("  tiled rounds: {} rounds, {:.3} s, {} ion-steps ({:.1}/visit)", round, t_start.elapsed().as_secs_f64(), tsteps, tsteps as f64 / visits.max(1) as f64); }
        round
    }

    /// Persistent-worker phase: dense list of live ions, advanced by the brick-layout direct kernel until all finish.
    /// `act[round % 2][0..ctr[round % 2])` is a superset of the live ions (same ping-pong as the tiled rounds).
    fn run_tail(&self, tp: &TPipes, cfg: &TileCfg, ib: &Buffer, fb: &Buffer, jb: &Buffer, jpb: &Buffer, pr0: &ChParams,
                act: &[Buffer; 2], ctr: &Buffer, args: &Buffer, mut round: usize, n_all: usize) {
        let last_tot = std::cell::Cell::new(0f64);
        let head = shared_buf(&self.device, 16); let cfgb = shared_buf(&self.device, 16);
        let mut pr = *pr0; let mut chunk = cfg.tail_chunk;
        let mut iter = 0; let mut prev_alive = n_all;
        let envn = |k: &str| std::env::var(k).ok().and_then(|v| v.parse::<usize>().ok());
        let bench_at = envn("CHSIM_BENCH_AT");
        // one compaction + worker launch of `chunk` steps per ion; returns (live ions in list, seconds)
        let batched = pr0.legacy != 0 && pr0.nearest != 0 && envn("CHSIM_NOBATCHCOLL").is_none();
        let batch_min = envn("CHSIM_BCMIN").unwrap_or(3000);
        let launch = |round: usize, chunk: i32, pr: &mut ChParams, use_leg: bool| -> (usize, f64) {
            let (ain, aout) = (&act[round % 2], &act[(round + 1) % 2]);
            let (c_in, c_out) = ((round % 2) as u64 * 8, ((round + 1) % 2) as u64 * 8);
            pr.nsteps = chunk;
            let cmd = self.queue.new_command_buffer();
            let enc = cmd.new_compute_command_encoder();
            enc.set_buffer(0, Some(ib), 0); enc.set_buffer(1, Some(fb), 0);
            enc.set_bytes(4, std::mem::size_of::<ChParams>() as u64, pr as *const ChParams as *const c_void);
            enc.set_buffer(5, Some(jpb), 0); enc.set_buffer(6, Some(jb), 0);
            enc.set_buffer(20, Some(&head), 0); enc.set_buffer(21, Some(&cfgb), 0);
            enc.set_buffer(14, Some(ain), 0); enc.set_buffer(15, Some(aout), 0);
            enc.set_buffer(16, Some(ctr), c_in); enc.set_buffer(17, Some(ctr), c_out); enc.set_buffer(18, Some(args), 0);
            enc.set_compute_pipeline_state(&tp.zero_args); enc.dispatch_thread_groups(MTLSize::new(1, 1, 1), MTLSize::new(32, 1, 1));
            enc.set_compute_pipeline_state(&tp.filter); enc.dispatch_thread_groups_indirect(args, 0, MTLSize::new(256, 1, 1));
            enc.set_compute_pipeline_state(&tp.size_args); enc.dispatch_thread_groups(MTLSize::new(1, 1, 1), MTLSize::new(32, 1, 1));
            enc.set_buffer(14, Some(aout), 0);      // the worker kernel reads the compacted list
            enc.set_compute_pipeline_state(if use_leg { &tp.tail_leg } else { &tp.tail }); enc.dispatch_thread_groups_indirect(args, 16, MTLSize::new(tail_tg() as u64, 1, 1));
            enc.end_encoding();
            let tt = Instant::now();
            cmd.commit(); cmd.wait_until_completed();
            let alive = unsafe { *((ctr.contents() as *const u8).add(c_out as usize) as *const u32) } as usize;
            (alive, tt.elapsed().as_secs_f64())
        };
        loop {
            if bench_at == Some(iter) && n_all >= envn("CHSIM_BENCH_MINN").unwrap_or(0) {
                // developer hook: replay this launch from a snapshot to compare kernel variants on realistic ion states
                let snap: Vec<ChIon> = unsafe { std::slice::from_raw_parts(ib.contents() as *const ChIon, n_all) }.to_vec();
                let bc = envn("CHSIM_BENCH_CHUNK").unwrap_or(20000) as i32;
                for rep in 0..5 {
                    unsafe { std::ptr::copy_nonoverlapping(snap.as_ptr(), ib.contents() as *mut ChIon, n_all); }
                    let (alive, dt) = launch(round, bc, &mut pr, batched && n_all >= batch_min);
                    println!("  bench rep {rep}: {} live ions, {} steps/ion in {:.4} s => {:.3} us/step/ion", alive, bc, dt, dt / bc as f64 * 1e6);
                }
                std::process::exit(0);
            }
            let (alive, dt_it) = launch(round, chunk, &mut pr, batched && prev_alive >= batch_min);
            prev_alive = alive;
            if self.verbose {
                let all = unsafe { std::slice::from_raw_parts(ib.contents() as *const ChIon, n_all) };
                let tot: f64 = all.iter().map(|i| i.v[3] as f64).sum::<f64>() / pr.dt as f64;
                eprintln!("    worker launch {iter}: {} steps/ion in {:.4} s, {} ions live; cumulative steps {:.3e} (+{:.3e} => {:.2} G ion-steps/s)",
                          chunk, dt_it, alive, tot, tot - last_tot.get(), (tot - last_tot.get()) / dt_it / 1e9);
                last_tot.set(tot);
            }
            chunk = (chunk.saturating_mul(2)).min(cfg.tail_cap);
            iter += 1; round += 1;
            if alive == 0 { break; }
        }
    }
}
