use crate::simulation::Simulation;
use mqviz::{analysis, engine::{Cmd, Params, SweepCfg, SweepView}, presets::PresetId, reality::{self, RealityData}, record, stats, world::TerrainTables};
use serde_json::{json, Value};
use std::{path::PathBuf, sync::{Mutex, OnceLock, atomic::{AtomicBool, Ordering}}, time::{Duration, Instant}};
use tauri::{Manager, State};

pub struct Research {
    pub data: OnceLock<RealityData>,
    pub state: Mutex<Value>,
    cancel: AtomicBool,
    directory: OnceLock<PathBuf>,
}
impl Research {
    pub fn new() -> Self { Self { data:OnceLock::new(), state:Mutex::new(json!({"running":false,"message":"Ready to compare with the paper.","rows":[]})), cancel:AtomicBool::new(false), directory:OnceLock::new() } }
    pub fn running(&self) -> bool {self.state.lock().unwrap()["running"]==true}
    pub fn init(&self, dir:PathBuf) -> Result<(),String> {
        std::fs::create_dir_all(&dir).map_err(|e|e.to_string())?;
        let data=dataset();
        for (name,file) in data["files"].as_object().ok_or("Missing embedded measurements")? {
            std::fs::write(dir.join(name),file["text"].as_str().unwrap()).map_err(|e|e.to_string())?;
        }
        self.data.set(RealityData::load(&dir)?).map_err(|_|"Already loaded")?;
        if let Ok(bytes)=std::fs::read(dir.join("last-experiment.json")) {
            if let Ok(mut v)=serde_json::from_slice::<Value>(&bytes) {
                if v["running"]==true {v["message"]=json!("Previous run interrupted; completed curves retained.");v["running"]=json!(false);}
                if let Some(rows)=v["rows"].as_array_mut(){for row in rows {annotate_leakage(row);}}
                if v["confirmation"].is_object(){annotate_leakage(&mut v["confirmation"]);}
                *self.state.lock().unwrap()=v;
            }
        }
        self.directory.set(dir).map_err(|_|"Already initialized")?;
        Ok(())
    }
    fn save(&self) -> Result<(),String> {
        let dir=self.directory.get().ok_or("Research storage unavailable")?;
        let bytes=serde_json::to_vec_pretty(&*self.state.lock().unwrap()).map_err(|e|e.to_string())?;
        let tmp=dir.join("last-experiment.tmp");std::fs::write(&tmp,bytes).map_err(|e|e.to_string())?;
        std::fs::rename(tmp,dir.join("last-experiment.json")).map_err(|e|e.to_string())
    }
}
fn dataset() -> Value {serde_json::from_str(include_str!("../../data/research.json")).expect("Bundled published data")}

pub fn reference_params(id:PresetId, d:&RealityData) -> Params {
    let mut p=Params::from_preset(id);
    if id==PresetId::S14 {p.species[0].weight=1.;p.species[1].weight=d.s14.adduct_over_main_open.v;p.f_ac_khz=83.;}
    p
}

#[tauri::command]
pub fn research_data() -> Value {dataset()}

#[tauri::command]
pub fn research_status(sim:State<'_,Simulation>, research:State<'_,Research>) -> Value {
    let s=sim.shared.snapshot();let p=sim.shared.params.lock().unwrap().clone();
    let d=research.data.get().unwrap();
    let (rows,context)=if p.preset==PresetId::S14 {
        (reality::s14_rows(d,&s.species,&s.stats,None),"S14 · reference drive 83 kHz; calibrated model RF ±70 V / AC ±0.35 V")
    } else {
        let (ch,_,_)=reality::nearest_channel(p.f_ac_khz,p.v_rf);
        (reality::fig5_rows(d,&s.species,&s.stats,ch),"Fig. 5 · nearest calibrated channel; source abundances from Fig. 5C")
    };
    let comparisons:Vec<_>=rows.iter().map(|r|json!({"name":if r.name=="efficiency T/(T+C), main ion" {"Ideal repeat-pass proxy vs measured signal ratio"}else{&r.name},"model":r.ours_txt,"published":r.theirs_txt,"source":if r.name.contains("source (preset input)") {"DIGITISED S14 Part I; source mixture set to measured Na/main ratio"}else{r.src},"settle":r.settle.map(|s|s.label()),"value":r.ours.v,"error":r.ours.e})).collect();
    let baseline=if p.preset==PresetId::S14 {(p.v_rf-70.).abs()<1e-6&&(p.v_ac-0.35).abs()<1e-6&&(p.f_ac_khz-83.).abs()<0.01&&(p.gate_v-24.).abs()<1e-6&&(p.pressure_pa-0.1).abs()<1e-6} else {let (_,f,delta)=reality::nearest_channel(p.f_ac_khz,p.v_rf);let _=f;(p.v_rf-56.).abs()<1e-6&&(p.v_ac-0.2).abs()<1e-6&&delta<0.05&&(p.gate_v-24.).abs()<1e-6&&(p.pressure_pa-0.1).abs()<1e-6};
    json!({"experiment":research.state.lock().unwrap().clone(),"comparisons":comparisons,"reference_settings":baseline,"context":context})
}

fn fit_quality(f:&analysis::SweepFit, sv:&SweepView) -> Option<String> {
    let step=sv.f[1]-sv.f[0];let lo=sv.f[0];let hi=*sv.f.last().unwrap();
    if !f.fit.converged || !f.r_mz.v.is_finite() || !f.r_mz.e.is_finite() {return Some("Fit did not converge".into());}
    if f.fwhm.v < 4.*step {return Some("Peak too narrow for this frequency grid".into());}
    if f.fit.p[2]-f.fwhm.v<lo || f.fit.p[2]+f.fwhm.v>hi {return Some("Peak wings extend beyond scan".into());}
    if f.r_mz.e > 0.3*f.r_mz.v {return Some("Resolution uncertainty is too large".into());}
    if f.fit.chi2_red()>10. {return Some("Single-peak fit misses the curve shape".into());}
    None
}
fn summarize(p:&Params, sv:&SweepView, tt:&TerrainTables, stage:&str) -> Value {
    let fits=analysis::fit_sweep(sv,&p.species,tt,p.v_rf,[18.,p.gate_v,15.,20.],100,p.seed);
    let points:Vec<_>=sv.f.iter().zip(&sv.t).map(|(f,r)|json!({"f":f,"species":r.iter().map(|s|json!({"k":s.k,"n":s.n,"t":s.p(),"sigma":s.sigma()})).collect::<Vec<_>>()})).collect();
    let mut row=json!({"ac":p.v_ac,"stage":stage,"seed":record::seed_str(p.seed),"points":points,"issue":"No significant peak","parameters":serde_json::from_str::<Value>(&record::params_j(p).pretty()).unwrap(),"scan":{"lo":sv.f[0],"hi":sv.f.last(),"step":sv.f[1]-sv.f[0],"ions_per_point":sv.cfg.as_ref().unwrap().ions_per_point},"sampling":"Fresh injections per frequency, shared across species; no recirculation. Raw finished counts retained."});
    if let Some(f)=&fits[0] {
        let i=sv.f.iter().enumerate().min_by(|a,b|(a.1-f.fit.p[2]).abs().total_cmp(&(b.1-f.fit.p[2]).abs())).unwrap().0;
        let t=sv.t[i][0];let leak=stats::ratio(sv.t[i][1],t);
        row["fit"]=json!({"f":f.fit.p[2],"f_error":f.fit.e[2],"width":f.fwhm.v,"width_error":f.fwhm.e,"r_mass":f.r_mz.v,"r_mass_error":f.r_mz.e,"r_frequency":f.r_f.v,"chi2":f.fit.chi2_red(),"shape":f.fit.shape.name(),"slope":f.dfdmz});
        row["frequency"]=json!(sv.f[i]);row["transmission"]=json!(t.p());row["transmission_error"]=json!(t.sigma());
        row["leakage"]=json!(leak.v);row["leakage_error"]=json!(leak.e);
        row["issue"]=json!(fit_quality(f,sv));
    }
    annotate_leakage(&mut row);
    row
}
// A zero transmitted count is an upper limit, never evidence of perfect rejection.
fn annotate_leakage(row:&mut Value) {
    let Some(f)=row["frequency"].as_f64() else{return};
    let Some(point)=row["points"].as_array().and_then(|p|p.iter().find(|p|p["f"].as_f64().is_some_and(|v|(v-f).abs()<1e-8))) else{return};
    let frac=|i:usize|stats::Frac::new(point["species"][i]["k"].as_u64().unwrap_or(0),point["species"][i]["n"].as_u64().unwrap_or(0));
    let (main,na)=(frac(0),frac(1));
    row["leakage_count"]=json!(na.k);row["leakage_upper95"]=json!(stats::ratio_ul95(na,main));
}
fn eligible(r:&Value,min_r:f64)->bool {
    r["issue"].is_null() && r["fit"]["r_mass"].as_f64().zip(r["fit"]["r_mass_error"].as_f64()).is_some_and(|(v,e)|v-2.*e>=min_r)
}
fn best(rows:&[Value],min_r:f64)->Option<usize> {
    rows.iter().enumerate().filter(|(_,r)|eligible(r,min_r)).max_by(|(_,a),(_,b)|{
        let score=|r:&Value|r["transmission"].as_f64().unwrap_or(0.)-2.*r["transmission_error"].as_f64().unwrap_or(1.);
        score(a).total_cmp(&score(b))
    }).map(|(i,_)|i)
}

fn sweep(sim:&Simulation,r:&Research,p:&Params,cfg:SweepCfg,label:&str)->Result<SweepView,String> {
    let old=sim.shared.snapshot().seq;
    *sim.shared.params.lock().unwrap()=p.clone();
    sim.shared.send(Cmd::Reset);
    let start=Instant::now();
    loop {
        if r.cancel.load(Ordering::Relaxed)||sim.shared.quit.load(Ordering::Relaxed) {return Err("Cancelled; live settings restored.".into());}
        let s=sim.shared.snapshot();
        // Gate rebuilding is debounced in mqviz. Wait for the fresh beam and the applied field.
        if s.seq>old+2 && s.warm>=1. && (s.gate_applied-p.gate_v).abs()<1e-6 && !s.rebuilding_dc && s.sweep.is_none() {break;}
        if start.elapsed()>Duration::from_secs(180) {return Err("Timed out preparing the beam".into());}
        std::thread::sleep(Duration::from_millis(100));
    }
    sim.shared.send(Cmd::SweepStart(cfg));
    loop {
        if r.cancel.load(Ordering::Relaxed)||sim.shared.quit.load(Ordering::Relaxed) {return Err("Cancelled; live settings restored.".into());}
        if start.elapsed()>Duration::from_secs(600) {return Err("Scan timed out; partial run retained".into());}
        if let Some(s)=sim.shared.snapshot().sweep.as_ref() {
            let mut state=r.state.lock().unwrap();state["message"]=json!(format!("{label} · {:.0}% · draining finished ions before fitting",s.frac*100.));state["progress"]=json!(s.frac);
            if !s.running {return Ok(s.clone());}
        }
        std::thread::sleep(Duration::from_millis(150));
    }
}
fn experiment(sim:&Simulation,r:&Research,min_r:f64,n:usize)->Result<(),String> {
    let tt=sim.tables.lock().unwrap().as_ref().unwrap().1.clone();
    let mut rows=vec![];
    for (i,ac) in [0.35,0.2625,0.175,0.0875].iter().enumerate() {
        let mut p=reference_params(PresetId::S14,r.data.get().unwrap());p.v_ac=*ac;p.recirculate=false;p.speed_cap=f64::INFINITY;p.seed=0xA110_0101+i as u64*2;p.inj_rate=3e6;
        let label=format!("Drive {} / 4 · AC ±{ac:.4} V",i+1);
        r.state.lock().unwrap()["message"]=json!(format!("{label} · preparing beam"));
        let cfg=SweepCfg{f_lo_khz:70.,f_hi_khz:100.,step_khz:0.5,ions_per_point:n,sweeps:1};
        let sv=sweep(sim,r,&p,cfg,&label)?;let coarse=summarize(&p,&sv,&tt,"exploratory");
        let mut row=coarse.clone();
        if let (Some(center),Some(width))=(row["fit"]["f"].as_f64(),row["fit"]["width"].as_f64()) {
            if width<2. {
                let half=(width*3.).clamp(2.,6.);p.seed+=100;
                let lo=(center-half).max(70.);let hi=(center+half).min(100.);
                let cfg=SweepCfg{f_lo_khz:lo,f_hi_khz:hi,step_khz:(hi-lo)/60.,ions_per_point:n*2,sweeps:1};
                let sv=sweep(sim,r,&p,cfg,&format!("{label} · resolving narrow peak"))?;
                row=summarize(&p,&sv,&tt,"refined");row["coarse"]=coarse;
            }
        }
        rows.push(row);r.state.lock().unwrap()["rows"]=json!(rows);r.save()?;
    }
    if let Some(i)=best(&rows,min_r) {
        let row=&rows[i];let mut p=reference_params(PresetId::S14,r.data.get().unwrap());p.v_ac=row["ac"].as_f64().unwrap();p.recirculate=false;p.speed_cap=f64::INFINITY;p.seed=0xC0FF_EE13;
        let cfg=SweepCfg{f_lo_khz:row["scan"]["lo"].as_f64().unwrap(),f_hi_khz:row["scan"]["hi"].as_f64().unwrap(),step_khz:row["scan"]["step"].as_f64().unwrap(),ions_per_point:n*3,sweeps:1};
        let sv=sweep(sim,r,&p,cfg,"Independent confirmation of leading candidate")?;
        let confirmed=summarize(&p,&sv,&tt,"confirmation");let passed=eligible(&confirmed,min_r);
        let mut s=r.state.lock().unwrap();s["confirmation"]=confirmed;s["selected"]=if passed {json!(i)} else {Value::Null};
        s["message"]=json!(if passed {"Candidate confirmed within this model. Compare the baseline mismatch before planning a hardware test."} else {"Leading candidate did not meet the resolution threshold on confirmation. No recommendation."});
    } else {r.state.lock().unwrap()["message"]=json!("No resolved candidate meets the conservative resolution threshold. Curves retained; lower the threshold or increase sampling.");}
    Ok(())
}

#[tauri::command]
pub fn research_start(app:tauri::AppHandle,sim:State<'_,Simulation>,r:State<'_,Research>,min_r:f64,ions:usize)->Result<(),String> {
    if !min_r.is_finite()||!(5.0..=150.0).contains(&min_r)||!(500..=5000).contains(&ions) {return Err("Choose resolving power 5–150 and 500–5000 ions per coarse point".into());}
    if sim.tables.lock().unwrap().is_none() {return Err("Wait for the fields to load".into());}
    let mut state=r.state.lock().unwrap();if state["running"]==true {return Err("An experiment is already running".into());}
    let saved=sim.shared.params.lock().unwrap().clone();
    *state=json!({"format":"multiq-experiment-1","running":true,"started":record::iso_now(),"message":"Preparing S14 reference scan…","rows":[],"min_r":min_r,"ions":ions,"objective":"Maximize target per-pass transmission minus 2 Monte Carlo standard errors, subject to fitted mass resolving power minus 2 standard errors >= requested minimum. Four AC drives; independent confirmation. Statistical screening, not a hardware optimum.","limits":"Calibrated code geometry, inferred 7 mm aperture, RF ±70 V / AC reference ±0.35 V. Paper reports AC 3 V; drive convention/transfer unverified. No space charge or trap-side delivery. R(m/z) uses local Mathieu slope; model uncertainty excluded. Seeds recorded; mqviz live ion IDs are not rewound, so reruns are statistical replicates, not bitwise replays."});
    r.cancel.store(false,Ordering::Relaxed);sim.study.store(true,Ordering::Relaxed);drop(state);
    std::thread::spawn(move||{
        let sim=app.state::<Simulation>();let r=app.state::<Research>();
        let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||experiment(&sim,&r,min_r,ions)));
        *sim.shared.params.lock().unwrap()=saved;sim.shared.send(Cmd::Reset);
        let mut state=r.state.lock().unwrap();state["running"]=json!(false);state["finished"]=json!(record::iso_now());
        match result {Ok(Ok(()))=>{},Ok(Err(e))=>{state["message"]=json!(e);state["selected"]=Value::Null;},Err(_)=>{state["message"]=json!("Experiment failed; live settings restored");state["selected"]=Value::Null;}}
        drop(state);sim.study.store(false,Ordering::Relaxed);if let Err(e)=r.save(){r.state.lock().unwrap()["storage_error"]=json!(e);}
    });
    Ok(())
}
#[tauri::command]
pub fn research_cancel(r:State<'_,Research>) {r.cancel.store(true,Ordering::Relaxed);}
#[tauri::command]
pub fn research_apply(sim:State<'_,Simulation>,r:State<'_,Research>,index:usize)->Result<(),String> {
    let s=r.state.lock().unwrap();if s["running"]==true {return Err("Wait for the scan or cancel it first".into());}
    let original=s["rows"].get(index).ok_or("Unknown candidate")?;
    let row=if s["selected"].as_u64()==Some(index as u64)&&s["confirmation"].is_object(){&s["confirmation"]}else{original};
    let mut p=reference_params(PresetId::S14,r.data.get().unwrap());p.v_ac=row["ac"].as_f64().ok_or("No drive")?;p.f_ac_khz=row["frequency"].as_f64().ok_or("No resolved peak")?;
    *sim.shared.params.lock().unwrap()=p;sim.shared.send(Cmd::Reset);Ok(())
}
#[tauri::command]
pub fn research_export(app:tauri::AppHandle,r:State<'_,Research>)->Result<String,String> {
    let mut value=r.state.lock().unwrap().clone();value["published_data"]=dataset();
    let dir=app.path().download_dir().map_err(|e|e.to_string())?.join("MultiQ Workbench");std::fs::create_dir_all(&dir).map_err(|e|e.to_string())?;
    for suffix in 0..1000 {
        let path=dir.join(format!("multiq-{}-{suffix}.json",record::fmt_time(record::now_secs()).1));
        match std::fs::OpenOptions::new().write(true).create_new(true).open(&path) {
            Ok(mut f)=>{use std::io::Write;f.write_all(&serde_json::to_vec_pretty(&value).unwrap()).map_err(|e|e.to_string())?;return Ok(path.to_string_lossy().into());},
            Err(e) if e.kind()==std::io::ErrorKind::AlreadyExists=>continue,Err(e)=>return Err(e.to_string())
        }
    }
    Err("Could not choose a unique export name".into())
}
#[cfg(test)] mod tests {
 use super::*;
 #[test] fn ranking_rejects_uncertain_and_invalid_peaks() {
  let rows=vec![json!({"issue":null,"fit":{"r_mass":55.,"r_mass_error":4.},"transmission":0.9,"transmission_error":0.01}),json!({"issue":null,"fit":{"r_mass":65.,"r_mass_error":4.},"transmission":0.2,"transmission_error":0.02}),json!({"issue":"under sampled","fit":{"r_mass":100.,"r_mass_error":1.},"transmission":0.99,"transmission_error":0.01})];
  assert_eq!(best(&rows,50.),Some(1));assert_eq!(best(&rows,80.),None);
 }
 #[test] fn embedded_measurements_have_provenance() {let d=dataset();assert!(d["files"].as_object().unwrap().len()>40);assert!(d["files"]["figS14_adduct_ratio.csv"]["text"].as_str().unwrap().contains("0.64,0.186"));}
 #[test] fn zero_adduct_counts_are_reported_as_a_limit() {let mut row=json!({"frequency":83.,"points":[{"f":83.,"species":[{"k":200,"n":1000},{"k":0,"n":100}]}]});annotate_leakage(&mut row);assert_eq!(row["leakage_count"],0);assert!(row["leakage_upper95"].as_f64().unwrap()>0.15);}
}
