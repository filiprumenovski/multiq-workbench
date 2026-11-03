#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]
mod tandem;
mod simulation;
mod research;
use tauri::Manager;

#[tauri::command]
fn smoke_report(app: tauri::AppHandle, success: bool, message: String, fixture: Option<serde_json::Value>) {
    if !std::env::args().any(|a|a=="--smoke-test"||a=="--research-test"||a=="--tandem-test") {return;}
    if let Some(data)=fixture {std::fs::write("/tmp/multiq-native-fixture.json",serde_json::to_vec(&data).unwrap()).unwrap();}
    println!("NATIVE SMOKE {}: {}",if success {"PASS"}else{"FAIL"},message);
    app.exit(if success {0}else{1});
}

fn main() {
    tauri::Builder::default()
        .manage(tandem::Tandem::default())
        .manage(simulation::Simulation::new())
        .manage(research::Research::new())
        .invoke_handler(tauri::generate_handler![tandem::tandem_export,tandem::original_sketch,tandem::tandem_start,tandem::tandem_stop,tandem::tandem_status,simulation::simulation_snapshot, simulation::simulation_control, simulation::simulation_scene, research::research_data, research::research_status, research::research_start, research::research_cancel, research::research_apply, research::research_export, smoke_report])
        .on_window_event(|window,event| {
            if let tauri::WindowEvent::Focused(active)=event {
                if !std::env::args().any(|a|a=="--smoke-test"||a=="--research-test"||a=="--tandem-test") {
                    window.state::<simulation::Simulation>().active.store(*active,std::sync::atomic::Ordering::Relaxed);
                }
            }
        })
        .on_page_load(|webview,payload| {
            if matches!(payload.event(),tauri::webview::PageLoadEvent::Finished) {
                if std::env::args().any(|a|a=="--tandem-test") {webview.eval(include_str!("../tests/tandem-smoke.js")).expect("Tandem smoke script");}
                else if std::env::args().any(|a|a=="--research-test") {webview.eval(include_str!("../tests/research-smoke.js")).expect("Research smoke script");}
                else if std::env::args().any(|a|a=="--smoke-test") {webview.eval(include_str!("../tests/native-smoke.js")).expect("Native smoke script");}
            }
        })
        .setup(|app| {
            let r=app.state::<research::Research>();
            r.init(app.path().app_data_dir()?.join("research")).map_err(std::io::Error::other)?;
            *app.state::<simulation::Simulation>().shared.params.lock().unwrap()=research::reference_params(mqviz::presets::PresetId::S14,r.data.get().unwrap());
            let cache=app.path().app_cache_dir()?.join("fields");
            std::fs::create_dir_all(&cache)?;
            if std::env::var_os("CHSIM_CACHE_DIR").is_none() {std::env::set_var("CHSIM_CACHE_DIR",cache);}
            simulation::start(app.handle().clone());
            Ok(())
        })
        .build(tauri::generate_context!())
        .expect("Could not start MultiQ Workbench")
        .run(|app,event| {if matches!(event,tauri::RunEvent::Exit) {app.state::<simulation::Simulation>().stop();}});
}
