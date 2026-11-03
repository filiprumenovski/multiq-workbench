(async()=>{
 const invoke=window.__TAURI__.core.invoke,wait=ms=>new Promise(r=>setTimeout(r,ms));
 const assert=(ok,msg)=>{if(!ok)throw new Error(msg);};
 async function until(fn,ms=180000){const end=Date.now()+ms;while(Date.now()<end){if(await fn())return;await wait(600);}throw new Error('Timed out: '+document.getElementById('research-status')?.textContent);}
 try{
  await until(()=>document.getElementById('trap-live')?.dataset.geometry==='native'&&Number(document.getElementById('trap-live').dataset.ions)>0);
  const data=await invoke('research_data');assert(Object.keys(data.files).length>40,'Published dataset must be bundled');
  await invoke('simulation_control',{controls:{preset:'s14',paused:true}});
  await until(async()=>(await invoke('simulation_snapshot')).warm===1);
  const saved=(await invoke('simulation_snapshot')).params;
  let rejected=false;try{await invoke('research_start',{minR:50,ions:1});}catch{rejected=true;}assert(rejected,'Reject unbounded/invalid requests');
  await invoke('research_start',{minR:50,ions:1000});
  rejected=false;try{await invoke('simulation_control',{controls:{ac:.1}});}catch{rejected=true;}assert(rejected,'Scan must own the controls');
  rejected=false;try{await invoke('research_start',{minR:50,ions:1000});}catch{rejected=true;}assert(rejected,'No overlapping scans');
  await invoke('research_cancel');await until(async()=>!(await invoke('research_status')).experiment.running);
  const restored=(await invoke('simulation_snapshot')).params;assert(JSON.stringify(saved)===JSON.stringify(restored),'Cancel must restore live settings');
  await invoke('simulation_control',{controls:{paused:false}});
  await invoke('research_start',{minR:50,ions:1000});
  await until(async()=>!(await invoke('research_status')).experiment.running,1500000);
  const result=await invoke('research_status'),ex=result.experiment;
  assert(ex.rows.length===4,'All four AC settings must produce a scan: '+ex.message);
  assert(ex.rows.every(r=>r.points.length>=50&&r.points.every(p=>p.species.length===2&&p.species.every(s=>s.n>0))),'Tagged species counts retained');
  if(ex.selected!=null){assert(ex.confirmation.stage==='confirmation'&&ex.confirmation.issue==null,'Recommendation requires independent confirmation');assert(ex.confirmation.fit.r_mass-2*ex.confirmation.fit.r_mass_error>=50,'Confirmed conservative R threshold');}
  const snapshot=await invoke('simulation_snapshot');await wait(2000);
  assert(!document.getElementById('research-error').textContent,'No research UI errors');
  assert(Number(document.getElementById('trap-live').dataset.entranceIons)>0,'Real ions must appear at the trap exit');
  const path=await invoke('research_export');assert(path.endsWith('.json'),'Export native data file');
  const fixture={snapshot,scene:await invoke('simulation_scene',{focus:0}),research:result};
  await invoke('smoke_report',{success:true,message:`Published data, 3D ions, exclusive scans, cancellation, 4 AC scans, confirmation screening and export verified. ${ex.message} Export: ${path}`,fixture});
 }catch(e){await invoke('smoke_report',{success:false,message:String(e),fixture:null});}
})();
