use serde_json::{json, Value};
use std::{
    sync::{
        atomic::{AtomicBool, Ordering},
        Arc, Mutex,
    },
    time::Instant,
};
use tauri::{Manager, State};
use trapsim::Model;
#[derive(Clone, Default)]
pub struct Tandem {
    busy: Arc<AtomicBool>,
    stop: Arc<AtomicBool>,
    state: Arc<Mutex<Value>>,
}
#[tauri::command]
pub fn tandem_status(t: State<Tandem>) -> Value {
    t.state.lock().unwrap().clone()
}
#[tauri::command]
pub fn tandem_stop(t: State<Tandem>) {
    t.stop.store(true, Ordering::Release);
}
#[tauri::command]
pub fn tandem_start(
    t: State<Tandem>,
    tones: [f64; 9],
    ions: usize,
    ms: u32,
    seed: u32,
    fine: bool,
    mouth: bool,
) -> Result<(), String> {
    if tones
        .iter()
        .any(|v| !v.is_finite() || (*v != 0. && !(20.0..=120.).contains(v)))
        || !(1..=512).contains(&ions)
        || !(1..=500).contains(&ms)
    {
        return Err("Use 1–512 ions, 1–500 ms, and 0 (off) or 20–120 kHz for each port.".into());
    }
    let t = t.inner().clone();
    if t.busy.swap(true, Ordering::AcqRel) {
        return Err("A tandem run is already active.".into());
    }
    t.stop.store(false, Ordering::Release);
    let settings = json!({"tones_khz":tones,"ions":ions,"duration_ms":ms,"seed":seed,"fine_grid":fine,"mouth_start":mouth,"dt_ns":40,"pressure_pa":0.1,"rf_v":56,"active_ac_v":0.2,"inactive_gate_v":28});
    *t.state.lock().unwrap() =
        json!({"running":true,"status":"Solving trap and channel fields…","settings":settings});
    std::thread::spawn(move || {
        let result = std::panic::catch_unwind(|| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(2)
                .build()
                .unwrap();
            let mut m = pool.install(|| Model::new(tones, ions, seed as u64, fine, mouth));
            let limit = ms as usize * 25000;
            let mut steps = 0;
            let mut last = Instant::now();
            while steps < limit && !t.stop.load(Ordering::Acquire) {
                let k = (limit - steps).min(128);
                m.advance(k);
                steps += k;
                if last.elapsed().as_millis() > 250 {
                    let mut s = t.state.lock().unwrap();
                    s["frame"] = json!(m.frame());
                    s["status"] = json!("Tracking ions");
                    last = Instant::now();
                }
            }
            let mut s = t.state.lock().unwrap();
            s["frame"] = json!(m.frame());
            s["status"] = json!(if steps == limit {
                "Complete"
            } else {
                "Stopped"
            });
        });
        let mut s = t.state.lock().unwrap();
        if result.is_err() {
            s["status"] = json!("Field solve or transport failed; result incomplete.");
        }
        s["running"] = json!(false);
        t.busy.store(false, Ordering::Release);
    });
    Ok(())
}
#[tauri::command]
pub fn original_sketch(app: tauri::AppHandle, channel: bool) -> Result<(), String> {
    let name = if channel {
        "MQ_MZ_ejection_singleQ"
    } else {
        "Demiurge_MultiQ_486"
    };
    let path = app
        .path()
        .resource_dir()
        .map_err(|e| e.to_string())?
        .join("reference")
        .join(name);
    let dest = app
        .path()
        .app_data_dir()
        .map_err(|e| e.to_string())?
        .join("original-sketches")
        .join(name);
    std::fs::create_dir_all(&dest).map_err(|e| e.to_string())?;
    for entry in std::fs::read_dir(path).map_err(|e| e.to_string())? {
        let f = entry.map_err(|e| e.to_string())?;
        let local = dest.join(f.file_name());
        if !local.exists() {
            std::fs::copy(f.path(), local).map_err(|e| e.to_string())?;
        }
    }
    let path = dest.join(format!("{name}.pde"));
    let status = std::process::Command::new("open")
        .arg("-a")
        .arg("Processing")
        .arg(path)
        .status()
        .map_err(|e| e.to_string())?;
    if status.success() {
        Ok(())
    } else {
        Err("Install Processing 4 to run the original sketch.".into())
    }
}
#[tauri::command]
pub fn tandem_export(t: State<Tandem>, app: tauri::AppHandle) -> Result<String, String> {
    let s = t.state.lock().unwrap().clone();
    if s["running"] == true || s["frame"].is_null() {
        return Err("Finish or stop a run before exporting.".into());
    }
    let dir = app
        .path()
        .download_dir()
        .map_err(|e| e.to_string())?
        .join("MultiQ Workbench");
    std::fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
    let path = dir.join(format!(
        "tandem-{}.json",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::write(
        &path,
        serde_json::to_vec_pretty(&s).map_err(|e| e.to_string())?,
    )
    .map_err(|e| e.to_string())?;
    Ok(path.display().to_string())
}
