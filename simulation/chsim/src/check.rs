//! `chsim check-tiled [n]`: bitwise comparison of the tiled/batched path against the direct per-ion kernel run job by job.
use crate::channel::*;
use std::time::Instant;

#[derive(Default)]
struct Diff { n: usize, differ: usize, p: usize, v: usize, aux: usize, rng: usize, status: usize, zero_sign_only: usize }

fn compare(a: &[ChIon], b: &[ChIon]) -> Diff {
    let mut d = Diff { n: a.len(), ..Default::default() };
    for (x, y) in a.iter().zip(b) {
        let pd = x.p[..3] != y.p[..3]; let sd = x.p[3] != y.p[3];
        let vd = (0..4).any(|k| x.v[k].to_bits() != y.v[k].to_bits());
        let ad = (0..4).any(|k| x.aux[k].to_bits() != y.aux[k].to_bits());
        let rd = x.rng != y.rng;
        d.p += pd as usize; d.status += sd as usize; d.v += vd as usize; d.aux += ad as usize; d.rng += rd as usize;
        d.differ += (pd || sd || vd || ad || rd) as usize;
        if !(pd || sd || rd) && (vd || ad) && (0..4).all(|k| x.v[k] == y.v[k]) && (0..4).all(|k| x.aux[k] == y.aux[k]) { d.zero_sign_only += 1; }
    }
    d
}

/// Returns number of differing ions. Reference = the direct per-ion kernel run job by job (each job alone). Compared against:
/// the direct kernel run as one batch, and the tiled engine in three configurations (workers only [default]; threadgroup-memory
/// tiled rounds until every ion is done; tiled rounds while >= 1000 ions live, then workers).
fn check_case(label: &str, eng: &mut Engine, jobs: &[Job], dc: [f64; 4], ph: &Physics, rounds_stride: usize) -> usize {
    let t = Instant::now();
    let mut alone: Vec<ChIon> = Vec::new();
    for j in jobs { alone.extend(eng.run_batch_ions(std::slice::from_ref(j), dc, ph, Path::Direct).0); }
    let t_alone = t.elapsed().as_secs_f64();
    let t = Instant::now();
    let (id, jt) = eng.run_batch_ions(jobs, dc, ph, Path::Direct);
    let t_dbatch = t.elapsed().as_secs_f64();
    let o_t = summarize(&alone, &jt, jobs.len());
    let tot = |f: fn(&Outcome) -> usize| o_t.iter().map(f).sum::<usize>();
    println!("[{label}] {} jobs, {} ions; outcomes: returned {} transmitted {} crashed {} left {} timeout {}", jobs.len(), alone.len(),
        tot(|o| o.returned), tot(|o| o.transmitted), tot(|o| o.crashed), tot(|o| o.left), tot(|o| o.timeout));
    let d2 = compare(&id, &alone);
    println!("   direct batch vs direct per-job              : {} / {} ions differ   [direct per-job {:.2} s, direct batch {:.2} s]", d2.differ, d2.n, t_alone, t_dbatch);
    let mut bad = d2.differ;
    // job subsets for the (slow) threadgroup-tiled-round variants: every `stride`-th job, compared against the matching slices
    let offs: Vec<usize> = jobs.iter().scan(0usize, |a, j| { let o = *a; *a += j.n; Some(o) }).collect();
    let pick = |stride: usize| -> (Vec<Job>, Vec<ChIon>) {
        let mut js = Vec::new(); let mut al = Vec::new();
        for k in (0..jobs.len()).step_by(stride) { js.push(jobs[k]); al.extend_from_slice(&alone[offs[k]..offs[k] + jobs[k].n]); }
        (js, al)
    };
    for (name, rounds, stride) in [("tiled engine (default: workers)    ", 0usize, 1usize), ("tiled rounds only (TG blocks)        ", 1, rounds_stride), ("tiled rounds >= 1000 live, then workers", 1000, 3)] {
        eng.tile.borrow_mut().rounds_min_live = rounds;
        let (js, al) = pick(stride);
        let t = Instant::now();
        let (it, _) = eng.run_batch_ions(&js, dc, ph, Path::Tiled);
        let dt = t.elapsed().as_secs_f64();
        let d = compare(&it, &al);
        println!("   {name} vs direct per-job: {} / {} ions differ (pos {} vel/t {} aux {} rng {} status {}; zero-sign-only {})   [{} jobs, {:.2} s]",
            d.differ, d.n, d.p, d.v, d.aux, d.rng, d.status, d.zero_sign_only, js.len(), dt);
        bad += d.differ;
    }
    eng.tile.borrow_mut().rounds_min_live = 0;
    bad
}

pub fn check_tiled(args: &[String]) {
    let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(300);
    let mut bad = 0;
    // (a) legacy Fig. 4 conditions: legacy physics, nearest-node fields, dt 1e-7
    {
        let eng = Engine::new(Fields::build(&Geom::default()));
        let mut eng = eng;
        let ph = Physics::legacy(1e-3 * 133.322);
        let dr = Drive { v_rf: 100.0, f_rf: 500e3, v_ac: 0.5, f_ac: 100e3, dc: [0.0, 5.0, -5.0, 0.0] };
        let jobs: Vec<Job> = (0..26).map(|i| 475.0 + 2.0 * i as f64).map(|mz| Job { sp: Species { mass: mz * 5.0, z: 5.0, ccs_n2: 0.0 }, drive: dr.rfac(), n, seed: 1000 + mz as u64 }).collect();
        bad += check_case("legacy Fig.4", &mut eng, &jobs, dr.dc, &ph, 9);
        // mixed switch combinations (legacy collisions with interpolated fields; thermal gas with nearest-node fields; thermal gas, AC on)
        let mut ph2 = Physics::legacy(1e-3 * 133.322); ph2.legacy_fields = false;
        bad += check_case("legacy collisions + trilinear fields", &mut eng, &jobs[8..14], dr.dc, &ph2, 6);
        let mut ph3 = Physics::correct(0.1); ph3.legacy_fields = true; ph3.t_max = 2e-3;
        let jobs3: Vec<Job> = (0..6).map(|i| Job { sp: Species { mass: 2.0 * 862.974, z: 2.0, ccs_n2: 425.0 }, drive: RfAc { v_rf: 75.0, f_rf: 500e3, v_ac: 1.0, f_ac: (70.0 + 4.0 * i as f64) * 1e3 }, n, seed: 50 + i }).collect();
        bad += check_case("thermal gas + nearest-node fields", &mut eng, &jobs3, [18.0, 24.0, 15.0, 20.0], &ph3, 2);
    }
    // (b) correct-physics S14: trilinear fields, thermal null-collision gas (0.1 Pa, 5 ms, AC 1.5 V)
    {
        let mut eng = s14_engine();
        let mut ph = Physics::correct(0.1); ph.t_max = 5e-3;
        let (closed, open) = s14_dc(&eng, None, false);
        let jobs: Vec<Job> = (70..=96).map(|f| Job { sp: bak2(), drive: RfAc { v_rf: 75.0, f_rf: 500e3, v_ac: 1.5, f_ac: f as f64 * 1e3 }, n, seed: 100 + f as u64 }).collect();
        bad += check_case("S14 closed gate, AC 1.5 V scan", &mut eng, &jobs, closed, &ph, 4);
        let jobs: Vec<Job> = vec![Job { sp: bak2(), drive: RfAc { v_rf: 75.0, f_rf: 500e3, v_ac: 0.0, f_ac: 83e3 }, n: n * 4, seed: 7 }];
        bad += check_case("S14 gate open, no AC", &mut eng, &jobs, open, &ph, 1);
    }
    println!("{}", if bad == 0 { "check-tiled: PASS (bitwise identical: position, velocity, time, aux, rng, status)" } else { "check-tiled: FAIL" });
    if bad != 0 { std::process::exit(1); }
}

fn fnv(ions: &[ChIon]) -> u64 {
    let mut h: u64 = 1469598103934665603;
    for b in unsafe { std::slice::from_raw_parts(ions.as_ptr() as *const u8, ions.len() * 64) } { h ^= *b as u64; h = h.wrapping_mul(1099511628211); }
    h
}

/// `chsim hash [n] [direct|tiled]`: FNV hash of the full final state of every ion for a fixed set of single jobs. The same jobs
/// run through the ORIGINAL (pre-batching) engine print identical hashes -- see PERF.md for how that comparison was made.
pub fn hash_jobs(args: &[String]) {
    let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(300);
    let path = if args.get(2).map(|s| s.as_str()) == Some("direct") { Path::Direct } else { Path::Tiled };
    let eng = Engine::new(Fields::build(&Geom::default()));
    let ph = Physics::legacy(1e-3 * 133.322);
    let dr = Drive { v_rf: 100.0, f_rf: 500e3, v_ac: 0.5, f_ac: 100e3, dc: [0.0, 5.0, -5.0, 0.0] };
    for mz in [475.0, 491.0, 499.0, 507.0, 523.0] {
        println!("legacy mz {}", mz);
        let (ions, _) = eng.run_batch_ions(&[Job { sp: Species { mass: mz * 5.0, z: 5.0, ccs_n2: 0.0 }, drive: dr.rfac(), n, seed: 1000 + mz as u64 }], dr.dc, &ph, path);
        println!("HASH {:016x} n={}", fnv(&ions), n);
    }
    let eng = s14_engine();
    let mut ph = Physics::correct(0.1); ph.t_max = 5e-3;
    let (closed, _) = s14_dc(&eng, None, false);
    for f in [75.0, 81.0, 83.0, 90.0] {
        println!("s14 f {}", f);
        let (ions, _) = eng.run_batch_ions(&[Job { sp: bak2(), drive: RfAc { v_rf: 75.0, f_rf: 500e3, v_ac: 1.5, f_ac: f * 1e3 }, n, seed: 100 + f as u64 }], closed, &ph, path);
        println!("HASH {:016x} n={}", fnv(&ions), n);
    }
}
