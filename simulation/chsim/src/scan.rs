//! Compute-frugal helpers: a persistent result cache (never rerun an identical job) and adaptive two-stage
//! frequency scans (coarse low-statistics pass, then high statistics only near resonances).
use crate::channel::{Engine, Job, Outcome, Physics, RfAc, Species};
use std::collections::HashMap;
use std::io::Write;

const CACHE: &str = "cache/results.tsv";

fn key(eng: &Engine, j: &Job, dc: [f64; 4], ph: &Physics) -> String {
    format!("{:?}|{:?}|{:?}|{}|{}|{:?}|{:?}", eng.fields.geom, j.sp, j.drive, j.n, j.seed, dc, ph)
}

fn load() -> HashMap<String, Outcome> {
    let mut m = HashMap::new();
    if let Ok(s) = std::fs::read_to_string(CACHE) {
        for line in s.lines() {
            let mut it = line.split('\t');
            if let (Some(k), Some(v)) = (it.next(), it.next()) {
                let f: Vec<f64> = v.split(',').filter_map(|x| x.parse().ok()).collect();
                if f.len() == 8 {
                    m.insert(k.to_string(), Outcome { n: f[0] as usize, crashed: f[1] as usize, returned: f[2] as usize, transmitted: f[3] as usize,
                                                      left: f[4] as usize, timeout: f[5] as usize, mean_coll: f[6], mean_time: f[7] });
                }
            }
        }
    }
    m
}

/// run_batch with a persistent cache keyed on (geometry, species, drive, n, seed, dc, physics).
pub fn cached_batch(eng: &Engine, jobs: &[Job], dc: [f64; 4], ph: &Physics) -> Vec<Outcome> {
    let mut cache = load();
    let keys: Vec<String> = jobs.iter().map(|j| key(eng, j, dc, ph)).collect();
    let miss: Vec<usize> = (0..jobs.len()).filter(|&i| !cache.contains_key(&keys[i])).collect();
    if !miss.is_empty() {
        let mj: Vec<Job> = miss.iter().map(|&i| jobs[i]).collect();
        let out = eng.run_batch(&mj, dc, ph);
        let _ = std::fs::create_dir_all("cache");
        let mut fh = std::fs::OpenOptions::new().create(true).append(true).open(CACHE).unwrap();
        for (&i, o) in miss.iter().zip(&out) {
            let _ = writeln!(fh, "{}\t{},{},{},{},{},{},{},{}", keys[i], o.n, o.crashed, o.returned, o.transmitted, o.left, o.timeout, o.mean_coll, o.mean_time);
            cache.insert(keys[i].clone(), *o);
        }
    }
    eprintln!("  [cache: {} hits, {} computed]", jobs.len() - miss.len(), miss.len());
    keys.iter().map(|k| cache[k]).collect()
}

/// Adaptive frequency scan: stage 1 at `f_step` with n_lo ions; stage 2 at `f_step/4` with n_hi ions only where
/// stage-1 transmission exceeds baseline by > 3 sigma (or within one step of the maximum).
/// Returns sorted (f_kHz, transmitted fraction, crashed fraction, n).
pub fn adaptive_scan(eng: &Engine, sp: Species, rf: RfAc, dc: [f64; 4], ph: &Physics, f_lo: f64, f_hi: f64, f_step: f64,
                     n_lo: usize, n_hi: usize, base: f64, seed0: u64) -> Vec<(f64, f64, f64, usize)> {
    let fs: Vec<f64> = (0..).map(|i| f_lo + f_step * i as f64).take_while(|&f| f <= f_hi + 1e-9).collect();
    let j1: Vec<Job> = fs.iter().enumerate().map(|(i, &f)| Job { sp, drive: RfAc { f_ac: f * 1e3, ..rf }, n: n_lo, seed: seed0 + i as u64 }).collect();
    let o1 = cached_batch(eng, &j1, dc, ph);
    let t1: Vec<f64> = o1.iter().map(|o| o.frac(o.transmitted)).collect();
    let tmax = t1.iter().cloned().fold(0.0, f64::max);
    let mut fine: Vec<f64> = vec![];
    for (i, &f) in fs.iter().enumerate() {
        // refine only the peak region(s): coarse points within 35% of the maximum, and only if clearly above baseline
        let sig = (t1[i].max(1.0 / n_lo as f64) / n_lo as f64).sqrt();
        if t1[i] >= 0.35 * tmax && t1[i] - base > 4.0 * sig {
            let mut g = f - f_step; while g <= f + f_step + 1e-9 { fine.push((g * 1e4).round() / 1e4); g += f_step / 4.0; }
        }
    }
    fine.sort_by(|a, b| a.partial_cmp(b).unwrap()); fine.dedup();
    let j2: Vec<Job> = fine.iter().enumerate().map(|(i, &f)| Job { sp, drive: RfAc { f_ac: f * 1e3, ..rf }, n: n_hi, seed: seed0 + 100_000 + i as u64 }).collect();
    let o2 = cached_batch(eng, &j2, dc, ph);
    let mut pts: Vec<(f64, f64, f64, usize)> = fs.iter().zip(&o1).filter(|(f, _)| !fine.iter().any(|g| (*g - **f).abs() < 1e-6))
        .map(|(&f, o)| (f, o.frac(o.transmitted), o.frac(o.crashed), n_lo)).collect();
    pts.extend(fine.iter().zip(&o2).map(|(&f, o)| (f, o.frac(o.transmitted), o.frac(o.crashed), n_hi)));
    pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let cost_full = (fs.len() as f64 * 4.0) * n_hi as f64;
    let cost = fs.len() as f64 * n_lo as f64 + fine.len() as f64 * n_hi as f64;
    eprintln!("  [adaptive scan: {} coarse x {} + {} fine x {} ions = {:.0}% of a uniform fine scan]", fs.len(), n_lo, fine.len(), n_hi, 100.0 * cost / cost_full);
    pts
}
