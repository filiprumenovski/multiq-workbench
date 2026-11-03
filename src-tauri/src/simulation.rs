use mqviz::{analysis, engine::{self, Cmd, Params, Shared}, presets::PresetId, terrain, view, world::{World, VizTables, TerrainTables}};
use serde::Deserialize;
use serde_json::{json, Value};
use std::sync::{Arc, Mutex, OnceLock, atomic::{AtomicBool, Ordering}};
use tauri::{Manager, State};

pub struct Simulation {
    pub shared: Arc<Shared>,
    pub active: Arc<AtomicBool>,
    pub study: Arc<AtomicBool>,
    pub tables: Mutex<Option<(Arc<VizTables>, Arc<TerrainTables>)>>,
    geometry: OnceLock<Value>,
    error: Mutex<Option<String>>,
}
impl Simulation {
    pub fn new() -> Self {
        Self { study:Arc::new(AtomicBool::new(false)), active:Arc::new(AtomicBool::new(true)), shared: Shared::new(Params::from_preset(PresetId::S14)), tables: Mutex::new(None), geometry:OnceLock::new(), error: Mutex::new(None) }
    }
    pub fn stop(&self) { self.shared.quit.store(true, Ordering::Relaxed); }
}

pub fn start(app: tauri::AppHandle) {
    std::thread::spawn(move || {
        let sim = app.state::<Simulation>();
        let built = std::panic::catch_unwind(World::build_calibrated);
        match built {
            Ok(world) => {
                let g=&world.geom;
                let _=sim.geometry.set(json!({"diameter":g.de,"center":g.center(),"start":g.q1_start,"length":g.q1_len,"tilt":g.alpha,"gap":g.gap,"plate":g.plate,"hole":g.hole_d,"q2_length":g.q2_len,"stem":g.stem_d,"axis":g.axis(),"grid":g.h_mm,"grow":g.grow,"layers":world.viz.layers}));
                *sim.tables.lock().unwrap() = Some((world.viz.clone(), world.terr.clone()));
                engine::spawn_desktop(world, sim.shared.clone(), sim.active.clone(), sim.study.clone());
            }
            Err(_) => *sim.error.lock().unwrap() = Some("The field solver failed to initialize. Restart the app to retry.".into()),
        }
    });
}

fn observation(o: analysis::Obs) -> Value {
    json!({"value":o.val.v,"sigma":o.val.e,"status":o.settle.label()})
}
fn params_json(p: &Params) -> Value {
    json!({"preset":match p.preset {PresetId::S14=>"s14",PresetId::Fig5=>"fig5",_=>"custom"},"rf":p.v_rf,"ac":p.v_ac,"frequency":p.f_ac_khz,"gate":p.gate_v,"pressure":p.pressure_pa,"paused":p.paused,"recirculate":p.recirculate})
}

#[tauri::command]
pub fn simulation_snapshot(sim: State<'_, Simulation>, r: State<'_, crate::research::Research>, positions: Option<bool>) -> Value {
    if let Some(error) = sim.error.lock().unwrap().as_ref() { return json!({"error":error}); }
    let ready = sim.tables.lock().unwrap().is_some();
    let snap = sim.shared.snapshot();
    let p = sim.shared.params.lock().unwrap().clone();
    let species: Vec<_> = snap.species.iter().enumerate().map(|(i,sp)| {
        let st = snap.stats.get(i).copied().unwrap_or_default();
        let rate = analysis::window_rate(&st, 0, snap.win_dur);
        json!({"label":sp.label,"mz":sp.mass/sp.z,"charge":sp.z,"transmission":observation(analysis::t_pass(&st)),"efficiency":observation(analysis::efficiency_sp(&st)),"current":{"value":rate.val.v*sp.z*1.602176634e-7,"sigma":rate.val.e*sp.z*1.602176634e-7},"counts":st.ep})
    }).collect();
    // Transfer visible positions only; the mqviz engine owns all integration and counting.
    let ions: Vec<_> = snap.ions.iter().take(if positions.unwrap_or(true) { usize::MAX } else { 0 }).map(|i| json!({"p":i.p,"trail":i.trail,"sp":i.sp,"hot":i.hot})).collect();
    json!({"background_suspended":!sim.active.load(Ordering::Relaxed)&&!r.running(),"ready":ready,"seq":snap.seq,"time_ms":snap.t_sim*1e3,"warm":snap.warm,"n_flight":snap.n_flight,"rtf":snap.rtf,"ions":ions,"species":species,"params":params_json(&p),"effective_frequency":snap.f_ac_eff_khz,"experiment_running":r.running(),"gate_applied":snap.gate_applied,"rebuilding":snap.rebuilding_dc,"epoch_ms":snap.epoch_ms,"dropped":snap.dropped,"timeouts":snap.timeouts})
}

#[derive(Deserialize, Default)]
#[serde(deny_unknown_fields)]
pub struct Controls {
    preset: Option<String>, rf: Option<f64>, ac: Option<f64>, frequency: Option<f64>,
    gate: Option<f64>, pressure: Option<f64>, paused: Option<bool>, recirculate: Option<bool>,
    action: Option<String>,
}
fn apply(p: &mut Params, c: &Controls) -> Result<Option<Cmd>, String> {
    let mut next = p.clone();
    if let Some(id) = &c.preset {
        next = Params::from_preset(match id.as_str() {"s14"=>PresetId::S14,"fig5"=>PresetId::Fig5,_=>return Err("Unknown preset".into())});
    }
    for (value,low,high,name,dest) in [
        (c.rf,40.,100.,"RF",&mut next.v_rf),(c.ac,0.,1.,"AC",&mut next.v_ac),
        (c.frequency,20.,120.,"frequency",&mut next.f_ac_khz),(c.gate,18.,30.,"gate",&mut next.gate_v),
        (c.pressure,0.02,0.3,"pressure",&mut next.pressure_pa),
    ] {
        if let Some(v)=value { if !v.is_finite() || v<low || v>high {return Err(format!("{name} must be between {low} and {high}"));} *dest=v; }
    }
    if let Some(v)=c.paused {next.paused=v;}
    if let Some(v)=c.recirculate {next.recirculate=v;}
    let cmd = match c.action.as_deref() {None=>None,Some("reset")=>Some(Cmd::Reset),Some("step")=>{next.paused=true;Some(Cmd::Step)},_=>return Err("Unknown action".into())};
    *p=next;
    Ok(cmd)
}

#[tauri::command]
pub fn simulation_control(sim: State<'_, Simulation>, r:State<'_,crate::research::Research>, controls: Controls) -> Result<Value, String> {
    let guard=r.state.lock().unwrap();
    if guard["running"]==true {return Err("The experiment owns the controls. Cancel it to tune manually.".into());}
    let mut p=sim.shared.params.lock().map_err(|_|"Simulation settings unavailable")?;
    let cmd=apply(&mut p,&controls)?;
    if controls.preset.is_some() {let reference=crate::research::reference_params(p.preset,r.data.get().unwrap());p.species=reference.species;if controls.frequency.is_none(){p.f_ac_khz=reference.f_ac_khz;}}
    let result=params_json(&p);
    if let Some(cmd)=cmd {sim.shared.send(cmd);}
    Ok(result)
}

#[tauri::command]
pub async fn simulation_scene(sim: State<'_, Simulation>, focus: usize) -> Result<Value, String> {
    let (v,t)=sim.tables.lock().map_err(|_|"Field data unavailable")?.as_ref().cloned().ok_or("Fields are loading")?;
    let p=sim.shared.params.lock().map_err(|_|"Settings unavailable")?.clone();
    let geometry=sim.geometry.get().cloned();
    if focus>=p.species.len() {return Err("Unknown species".into());}
    tauri::async_runtime::spawn_blocking(move || {
        let dc=[18.,p.gate_v,15.,20.];
        let st=view::build_statics(&v);
        let topo=view::build_topo(&v,&p.species[focus],p.v_rf,dc);
        let terrain=terrain::compute(&t,&p.species,p.v_rf,dc);
        let contours:Vec<_>=topo.levels.iter().map(|l|json!({"value":l.n as f32*topo.step,"index":l.n.rem_euclid(5)==0,"lines":l.lines})).collect();
        json!({"h":v.h_mm,"axis":v.axis_mm,"geometry":geometry,"cut":st.cut_rects,"ghost":st.ghost_rects,"end":st.end_rects,"contours":contours,"terrain":{"z":terrain.z_mm,"fx":terrain.fx,"fy":terrain.fy},"focus":focus,"params":params_json(&p)})
    }).await.map_err(|e|e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn control_validation_is_atomic_and_preserves_physics_defaults() {
        let mut p=Params::from_preset(PresetId::S14);
        let original=p.clone();
        let bad=Controls {rf:Some(80.),pressure:Some(-1.),..Default::default()};
        assert!(apply(&mut p,&bad).is_err());assert!(p==original);
        assert!(apply(&mut p,&Controls {frequency:Some(f64::NAN),..Default::default()}).is_err());
        assert!(apply(&mut p,&Controls {frequency:Some(82.5),..Default::default()}).is_ok());
        assert_eq!(p.f_ac_khz,82.5);assert_eq!(p.seed,original.seed);
        assert!(apply(&mut p,&Controls {action:Some("step".into()),..Default::default()}).unwrap().is_some());
        assert!(p.paused);
    }
}
