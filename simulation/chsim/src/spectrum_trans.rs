//! Channel transmission interface for the spectrum forward model: T_k(species) from the cached simulation results
//! (`cache/results.tsv`, read-only: this module never launches a simulation), with an explicit provenance for every
//! value and an exact list of the jobs that would turn missing entries into exact ones.
//!
//! Lookup order: Exact (pooled over every cached seed with identical species/drive/DC) > Interp (neighbouring m/z of the
//! same charge in a cached m/z scan) > Proxy (adducts: parent charge state) > Missing (T = 0, flagged in the report).

use std::collections::HashMap;

#[derive(Clone, Copy, Debug, Default)]
pub struct Counts { pub n: usize, pub crashed: usize, pub transmitted: usize }
impl Counts { fn add(&mut self, o: &Counts) { self.n += o.n; self.crashed += o.crashed; self.transmitted += o.transmitted; } }

#[derive(Clone, Debug)]
#[allow(dead_code)]
struct Entry { mass: f64, z: f64, ccs: f64, vrf: f64, frf: f64, vac: f64, fac: f64, dc: [f64; 4], n: usize, crashed: usize, transmitted: usize, seed: u64 }

fn num(s: &str, key: &str) -> Option<f64> {
    let i = s.find(key)? + key.len();
    let t: String = s[i..].chars().take_while(|c| c.is_ascii_digit() || matches!(c, '.' | '-' | 'e' | 'E' | '+')).collect();
    t.parse().ok()
}

pub struct TCache { e: Vec<Entry>, by_z: HashMap<i32, Vec<usize>>, pub n_lines: usize, pub n_other_geom: usize }

impl TCache {
    /// Load only the entries whose geometry / physics debug strings equal the given ones.
    pub fn load(path: &str, geom: &str, phys: &str) -> Result<TCache, String> {
        let text = std::fs::read_to_string(path).map_err(|e| format!("cannot read cache {}: {}", path, e))?;
        let (mut e, mut n_lines, mut other) = (vec![], 0, 0);
        for line in text.lines() {
            n_lines += 1;
            let mut it = line.split('\t');
            let (Some(k), Some(v)) = (it.next(), it.next()) else { continue };
            let p: Vec<&str> = k.split('|').collect();
            if p.len() != 7 { continue; }
            if p[0] != geom || p[6] != phys { other += 1; continue; }
            let f: Vec<f64> = v.split(',').filter_map(|x| x.parse().ok()).collect();
            if f.len() != 8 { continue; }
            let dc: Vec<f64> = p[5].trim_matches(|c| c == '[' || c == ']').split(',').filter_map(|x| x.trim().parse().ok()).collect();
            if dc.len() != 4 { continue; }
            let (Some(mass), Some(z), Some(ccs)) = (num(p[1], "mass: "), num(p[1], "z: "), num(p[1], "ccs_n2: ")) else { continue };
            let (Some(vrf), Some(frf), Some(vac), Some(fac)) = (num(p[2], "v_rf: "), num(p[2], "f_rf: "), num(p[2], "v_ac: "), num(p[2], "f_ac: ")) else { continue };
            e.push(Entry { mass, z, ccs, vrf, frf, vac, fac, dc: [dc[0], dc[1], dc[2], dc[3]], n: f[0] as usize, crashed: f[1] as usize, transmitted: f[3] as usize, seed: p[4].parse().unwrap_or(0) });
        }
        let mut by_z: HashMap<i32, Vec<usize>> = HashMap::new();
        for (i, x) in e.iter().enumerate() { by_z.entry(x.z.round() as i32).or_default().push(i); }
        Ok(TCache { e, by_z, n_lines, n_other_geom: other })
    }
    pub fn len(&self) -> usize { self.e.len() }

    fn matches_sp(x: &Entry, mass: f64, z: f64, ccs: f64) -> bool { (x.mass - mass).abs() <= 1e-6 * mass && (x.z - z).abs() < 1e-9 && (x.ccs - ccs).abs() <= 1e-6 * ccs.max(1.0) }
    fn same_dc(a: &[f64; 4], b: &[f64; 4]) -> bool { a.iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-9) }

    /// Pooled counts of identical (species, RF, AC amplitude, AC frequency, DC) jobs over every cached seed.
    /// `fac = None` pools every AC frequency (used for the AC-off baseline, where the frequency is irrelevant).
    pub fn pooled(&self, mass: f64, z: f64, ccs: f64, vrf: f64, vac: f64, fac: Option<f64>, dc: &[f64; 4]) -> Option<Counts> {
        let mut c = Counts::default();
        let mut hit = false;
        for &i in self.by_z.get(&(z.round() as i32))? {
            let x = &self.e[i];
            if !Self::matches_sp(x, mass, z, ccs) || (x.vrf - vrf).abs() > 1e-9 || (x.vac - vac).abs() > 1e-9 || !Self::same_dc(&x.dc, dc) { continue; }
            if let Some(f) = fac { if (x.fac - f).abs() > 1e-6 { continue; } }
            c.add(&Counts { n: x.n, crashed: x.crashed, transmitted: x.transmitted }); hit = true;
        }
        if hit { Some(c) } else { None }
    }

    /// Entries of one (z, ccs-matching-rule) m/z scan at fixed drive: (mass/z, counts), for interpolation.
    fn scan_points(&self, z: f64, vrf: f64, vac: f64, fac: Option<f64>, dc: &[f64; 4], ccs_rule: &dyn Fn(f64, f64) -> f64) -> Vec<(f64, Counts)> {
        let mut m: std::collections::BTreeMap<i64, (f64, Counts)> = Default::default();
        if let Some(ix) = self.by_z.get(&(z.round() as i32)) {
            for &i in ix {
                let x = &self.e[i];
                let mz = x.mass / x.z;
                if (x.vrf - vrf).abs() > 1e-9 || (x.vac - vac).abs() > 1e-9 || !Self::same_dc(&x.dc, dc) { continue; }
                if (x.ccs - ccs_rule(mz, z)).abs() > 1e-6 * x.ccs { continue; }
                if let Some(f) = fac { if (x.fac - f).abs() > 1e-6 { continue; } }
                let e = m.entry((mz * 1e4).round() as i64).or_insert((mz, Counts::default()));
                e.1.add(&Counts { n: x.n, crashed: x.crashed, transmitted: x.transmitted });
            }
        }
        m.into_values().collect()
    }

    /// Per-(species, f) entry with the largest n (what the adaptive scans use for their argmax): (f_kHz, counts).
    pub fn freq_scan_best_n(&self, mass: f64, z: f64, ccs: f64, vrf: f64, vac: f64, dc: &[f64; 4]) -> Vec<(f64, Counts)> {
        let mut m: std::collections::BTreeMap<i64, (f64, Counts, u64)> = Default::default();
        if let Some(ix) = self.by_z.get(&(z.round() as i32)) {
            for &i in ix {
                let x = &self.e[i];
                if !Self::matches_sp(x, mass, z, ccs) || (x.vrf - vrf).abs() > 1e-9 || (x.vac - vac).abs() > 1e-9 || !Self::same_dc(&x.dc, dc) { continue; }
                let k = (x.fac * 1e-2).round() as i64;
                let c = Counts { n: x.n, crashed: x.crashed, transmitted: x.transmitted };
                match m.get(&k) { Some(old) if old.1.n >= c.n => {}, _ => { m.insert(k, (x.fac * 1e-3, c, x.seed)); } }
            }
        }
        m.into_values().map(|(f, c, _)| (f, c)).collect()
    }
    /// The exact job (one seed) `chsim reality` used: pooled over nothing, the entry with this seed.
    pub fn exact_seed(&self, mass: f64, z: f64, ccs: f64, vrf: f64, vac: f64, fac: Option<f64>, dc: &[f64; 4], seed: u64) -> Option<Counts> {
        for &i in self.by_z.get(&(z.round() as i32))? {
            let x = &self.e[i];
            if x.seed != seed || !Self::matches_sp(x, mass, z, ccs) || (x.vrf - vrf).abs() > 1e-9 || (x.vac - vac).abs() > 1e-9 || !Self::same_dc(&x.dc, dc) { continue; }
            if let Some(f) = fac { if (x.fac - f).abs() > 1e-6 { continue; } }
            return Some(Counts { n: x.n, crashed: x.crashed, transmitted: x.transmitted });
        }
        None
    }
    /// The points of one `scan::adaptive_scan` (coarse seeds seed0.., fine seeds seed0+100000..): fine where it exists, else coarse.
    pub fn adaptive_scan_pts(&self, mass: f64, z: f64, ccs: f64, vrf: f64, vac: f64, dc: &[f64; 4], seed0: u64) -> Vec<(f64, Counts)> {
        let mut m: std::collections::BTreeMap<i64, (f64, Counts, bool)> = Default::default();
        if let Some(ix) = self.by_z.get(&(z.round() as i32)) {
            for &i in ix {
                let x = &self.e[i];
                if !Self::matches_sp(x, mass, z, ccs) || (x.vrf - vrf).abs() > 1e-9 || (x.vac - vac).abs() > 1e-9 || !Self::same_dc(&x.dc, dc) { continue; }
                let (coarse, fine) = (x.seed >= seed0 && x.seed < seed0 + 100_000, x.seed >= seed0 + 100_000 && x.seed < seed0 + 200_000);
                if !coarse && !fine { continue; }
                let k = (x.fac * 1e-2).round() as i64;
                let c = Counts { n: x.n, crashed: x.crashed, transmitted: x.transmitted };
                match m.get(&k) { Some(old) if old.2 || !fine => {}, _ => { m.insert(k, (x.fac * 1e-3, c, fine)); } }
            }
        }
        m.into_values().map(|(f, c, _)| (f, c)).collect()
    }
}

// ---------------------------------------------------------------- species -> simulated Species

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Kind { Ubq(u32), Chem(u32), Adduct(u32) }

/// N2 CCS (A^2) used when a species has to be simulated. Ubiquitin 6..10 are the values in the cache (docs/ion_parameters.md);
/// 5, 11, 12 are NOT in the cache and are estimates (flagged in the job list).
pub fn ubq_ccs(z: u32) -> f64 { match z { 0..=6 => 1250.0, 7 => 1262.0, 8 => 1800.0, 9 => 2090.0, 10 => 2200.0, 11 => 2300.0, _ => 2400.0 } }
/// Small-molecule rule used by `chsim predict` (so its cached m/z scans are found): 1+: 9.41 m^0.52; z>1: 425 (m/1726.9)^(2/3).
pub fn chem_ccs(mz: f64, z: f64) -> f64 { let m = mz * z; if z < 1.5 { 9.41 * m.powf(0.52) } else { 425.0 * (m / 1726.9).powf(2.0 / 3.0) } }

#[derive(Clone, Debug)]
pub struct SimSp { pub name: String, pub kind: Kind, pub mz: f64, pub z: f64, pub mass: f64, pub ccs: f64 }

impl SimSp {
    pub fn new(name: &str, kind: Kind, mz: f64) -> SimSp {
        match kind {
            Kind::Ubq(z) => { let zf = z as f64; let mass = crate::spectrum_detector::M_UBQ_AVG + zf * 1.00728; SimSp { name: name.into(), kind, mz: mass / zf, z: zf, mass, ccs: ubq_ccs(z) } }
            Kind::Chem(z) => { let zf = z as f64; SimSp { name: name.into(), kind, mz, z: zf, mass: mz * zf, ccs: chem_ccs(mz, zf) } }
            Kind::Adduct(z) => {
                let zf = z as f64; let mass = mz * zf; let pm = crate::spectrum_detector::M_UBQ_AVG + zf * 1.00728;
                SimSp { name: name.into(), kind, mz, z: zf, mass, ccs: ubq_ccs(z) * (mass / pm).powf(2.0 / 3.0) }
            }
        }
    }
}

// ---------------------------------------------------------------- downstream factor D(m/z)

#[derive(Clone, Debug)]
pub enum Downstream { Unity, /// logistic low-mass cutoff: D = 1/(1+exp(-(mz-mz0)/w))
    Cutoff { mz0: f64, w: f64 }, Table(Vec<(f64, f64)>) }
impl Downstream {
    pub fn parse(s: &str) -> Result<Downstream, String> {
        let p: Vec<&str> = s.trim().split(':').collect();
        match p[0] {
            "" | "unity" | "1" => Ok(Downstream::Unity),
            "cutoff" if p.len() == 3 => Ok(Downstream::Cutoff { mz0: p[1].parse().map_err(|e| format!("{}", e))?, w: p[2].parse().map_err(|e| format!("{}", e))? }),
            "file" if p.len() == 2 => {
                let t = std::fs::read_to_string(p[1]).map_err(|e| format!("downstream file: {}", e))?;
                let mut v: Vec<(f64, f64)> = t.lines().filter(|l| !l.starts_with('#')).filter_map(|l| { let q: Vec<f64> = l.split(|c| c == ',' || c == ' ' || c == '\t').filter_map(|x| x.parse().ok()).collect(); if q.len() >= 2 { Some((q[0], q[1])) } else { None } }).collect();
                v.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap()); Ok(Downstream::Table(v))
            }
            _ => Err(format!("downstream = unity | cutoff:<mz0>:<width> | file:<path>  (got {})", s)),
        }
    }
    pub fn d(&self, mz: f64) -> f64 {
        match self {
            Downstream::Unity => 1.0,
            Downstream::Cutoff { mz0, w } => 1.0 / (1.0 + (-(mz - mz0) / w).exp()),
            Downstream::Table(t) => {
                if t.is_empty() { return 1.0; }
                if mz <= t[0].0 { return t[0].1; }
                for w in t.windows(2) { if mz <= w[1].0 { let f = (mz - w[0].0) / (w[1].0 - w[0].0); return w[0].1 + f * (w[1].1 - w[0].1); } }
                t.last().unwrap().1
            }
        }
    }
}

// ---------------------------------------------------------------- transmission values

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Prov { Exact, Interp, Proxy, Missing }
impl Prov { pub fn tag(&self) -> &'static str { match self { Prov::Exact => "exact", Prov::Interp => "interp", Prov::Proxy => "proxy", Prov::Missing => "MISSING" } } }

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TMode { On, Sub, Steady }
impl TMode { pub fn parse(s: &str) -> Result<TMode, String> { match s { "on" => Ok(TMode::On), "sub" => Ok(TMode::Sub), "steady" => Ok(TMode::Steady), _ => Err(format!("tmode = on | sub | steady (got {})", s)) } } }

#[derive(Clone, Copy, Debug)]
#[allow(dead_code)]
pub struct TVal { pub t: f64, pub rel: f64, pub prov: Prov, pub t_on: f64, pub t_base: f64, pub k_on: usize, pub n_on: usize }
impl TVal { pub fn missing() -> TVal { TVal { t: 0.0, rel: 1.0, prov: Prov::Missing, t_on: 0.0, t_base: 0.0, k_on: 0, n_on: 0 } } }

/// Operating point of the channels (everything but the AC frequency).
#[derive(Clone, Copy, Debug)]
pub struct Point { pub vrf: f64, pub frf: f64, pub vac: f64, pub dc: [f64; 4] }

/// Seeds for newly requested jobs start here (cache keys include the seed; unique and deterministic).
pub const JOB_SEED0: u64 = 8_000_000;

#[derive(Clone, Debug)]
pub struct Need { pub label: String, pub mass: f64, pub z: f64, pub ccs: f64, pub vrf: f64, pub frf: f64, pub vac: f64, pub fac: f64, pub dc: [f64; 4], pub n: usize }

pub struct TModel { pub cache: TCache, pub mode: TMode, pub down: Downstream, pub n_new: usize, pub interp_max_dmz: f64 }

/// Jeffreys estimate of a binomial fraction: (k + 1/2) / (n + 1) with relative error 1/sqrt(k + 1/2) (never a claim of exactly zero).
fn jeffreys(k: usize, n: usize) -> (f64, f64) { ((k as f64 + 0.5) / (n as f64 + 1.0), 1.0 / (k as f64 + 0.5).sqrt()) }

impl TModel {
    fn counts_to(&self, on: Counts, base: Option<Counts>, prov: Prov) -> TVal {
        let n = on.n.max(1);
        let (t_on, r_on) = jeffreys(on.transmitted, n);
        let (t_base, r_base) = base.map_or((0.0, 0.0), |b| jeffreys(b.transmitted, b.n.max(1)));
        let (t, rel) = match self.mode {
            TMode::On => (t_on, r_on),
            TMode::Sub => {
                let t = (t_on - t_base).max(0.5 / (n as f64 + 1.0));
                let var = (t_on * r_on).powi(2) + (t_base * r_base).powi(2);
                (t, (var.sqrt() / t).min(1.5))
            }
            TMode::Steady => { let d = on.transmitted + on.crashed; let (f, r) = jeffreys(on.transmitted, d.max(1)); (f, r) }
        };
        TVal { t, rel: rel.min(1.0), prov, t_on, t_base, k_on: on.transmitted, n_on: on.n }
    }

    fn exact(&self, sp: &SimSp, pt: &Point, fac_khz: f64) -> Option<TVal> {
        let on = self.cache.pooled(sp.mass, sp.z, sp.ccs, pt.vrf, pt.vac, Some(fac_khz * 1e3), &pt.dc)?;
        let base = self.cache.pooled(sp.mass, sp.z, sp.ccs, pt.vrf, 0.0, None, &pt.dc);
        Some(self.counts_to(on, base, Prov::Exact))
    }

    fn interp(&self, sp: &SimSp, pt: &Point, fac_khz: f64) -> Option<TVal> {
        let rule: &dyn Fn(f64, f64) -> f64 = &chem_ccs;
        if !matches!(sp.kind, Kind::Chem(_)) { return None; }
        let on = self.cache.scan_points(sp.z, pt.vrf, pt.vac, Some(fac_khz * 1e3), &pt.dc, rule);
        let off = self.cache.scan_points(sp.z, pt.vrf, 0.0, None, &pt.dc, rule);
        let bracket = |pts: &[(f64, Counts)]| -> Option<(f64, Counts, f64, Counts)> {
            let lo = pts.iter().filter(|p| p.0 <= sp.mz + 1e-9).last()?; let hi = pts.iter().find(|p| p.0 >= sp.mz - 1e-9)?;
            if sp.mz - lo.0 > self.interp_max_dmz || hi.0 - sp.mz > self.interp_max_dmz { return None; }
            Some((lo.0, lo.1, hi.0, hi.1))
        };
        let (l0, c0, l1, c1) = bracket(&on)?;
        let w = if l1 > l0 { (sp.mz - l0) / (l1 - l0) } else { 0.0 };
        let mix = |a: Counts, b: Counts| -> Counts {
            // interpolate the fractions, keep the smaller sample size for the error estimate
            let (fa, fb) = (a.transmitted as f64 / a.n as f64, b.transmitted as f64 / b.n as f64);
            let (ca, cb) = (a.crashed as f64 / a.n as f64, b.crashed as f64 / b.n as f64);
            let n = a.n.min(b.n);
            Counts { n, transmitted: (((1.0 - w) * fa + w * fb) * n as f64).round() as usize, crashed: (((1.0 - w) * ca + w * cb) * n as f64).round() as usize }
        };
        let onc = mix(c0, c1);
        let basec = bracket(&off).map(|(_, a, _, b)| mix(a, b));
        Some(self.counts_to(onc, basec, Prov::Interp))
    }

    /// T of species `sp` through a channel tuned to `fac_khz` at operating point `pt`; `parent` = the parent species of an
    /// adduct (its T is used as a proxy when the adduct itself is not cached).
    pub fn t(&self, sp: &SimSp, pt: &Point, fac_khz: f64, parent: Option<&SimSp>) -> TVal {
        if let Some(v) = self.exact(sp, pt, fac_khz) { return v; }
        if let Some(v) = self.interp(sp, pt, fac_khz) { return v; }
        if let (Kind::Adduct(_), Some(p)) = (sp.kind, parent) {
            if let Some(mut v) = self.exact(p, pt, fac_khz) { v.prov = Prov::Proxy; return v; }
        }
        TVal::missing()
    }

    /// The exact jobs that would make `t(...)` Exact (AC-on at each channel frequency, plus the AC-off baseline when the mode needs it).
    pub fn needs(&self, sp: &SimSp, pt: &Point, fac_khz: f64, label: &str, n: usize) -> Vec<Need> {
        let mut v = vec![];
        if self.cache.pooled(sp.mass, sp.z, sp.ccs, pt.vrf, pt.vac, Some(fac_khz * 1e3), &pt.dc).is_none() {
            v.push(Need { label: format!("{} {} ch@{:.2}kHz", label, sp.name, fac_khz), mass: sp.mass, z: sp.z, ccs: sp.ccs, vrf: pt.vrf, frf: pt.frf, vac: pt.vac, fac: fac_khz * 1e3, dc: pt.dc, n });
        }
        if self.cache.pooled(sp.mass, sp.z, sp.ccs, pt.vrf, 0.0, None, &pt.dc).is_none() {
            v.push(Need { label: format!("{} {} AC-off baseline", label, sp.name), mass: sp.mass, z: sp.z, ccs: sp.ccs, vrf: pt.vrf, frf: pt.frf, vac: 0.0, fac: 50e3, dc: pt.dc, n });
        }
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn num_parser() {
        assert_eq!(num("Species { mass: 8570.84368, z: 6.0, ccs_n2: 1250.0 }", "z: "), Some(6.0));
        assert_eq!(num("RfAc { v_rf: 60.0, f_rf: 500000.0, v_ac: 0.27, f_ac: 1e5 }", "f_ac: "), Some(1e5));
    }
    #[test]
    fn downstream_parse_and_eval() {
        let d = Downstream::parse("cutoff:444:30").unwrap();
        assert!(d.d(300.0) < 0.01 && d.d(700.0) > 0.99 && (d.d(444.0) - 0.5).abs() < 1e-12);
        assert_eq!(Downstream::parse("unity").unwrap().d(5.0), 1.0);
    }
}
