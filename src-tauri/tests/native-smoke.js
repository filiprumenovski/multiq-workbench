(async()=>{
 const invoke=window.__TAURI__.core.invoke,wait=ms=>new Promise(r=>setTimeout(r,ms));
 const assert=(ok,message)=>{if(!ok)throw new Error(message);};
 async function until(fn){const deadline=Date.now()+120000;while(Date.now()<deadline){if(await fn())return;await wait(200);}throw new Error('Timed out: '+document.getElementById('live-error')?.textContent+' '+document.getElementById('live-status')?.textContent);}
 try{
  await until(()=>document.getElementById('chait-multiq')?.dataset.liveScene==='ready'&&!document.getElementById('live-controls').disabled);
  const a=await invoke('simulation_snapshot');await wait(400);const b=await invoke('simulation_snapshot');
  assert(b.time_ms>a.time_ms&&b.ions.length>0,'Live ions must advance');
  assert(b.ions.every(i=>i.p.length===3&&i.trail.length===4),'Real mqviz ion positions and trails');
  const compact=await invoke('simulation_snapshot',{positions:false});assert(compact.ions.length===0&&compact.species.length===b.species.length&&compact.n_flight>0,'Hidden views retain statistics without transferring positions');
  await until(()=>Number(document.getElementById('chait-multiq').dataset.frames)>0);
  assert(Number(document.getElementById('chait-multiq').dataset.frames)>0,'3D renderer must draw');
  const cloudRoot=document.getElementById('chait-multiq');const cloudTime=Number(cloudRoot.dataset.cloudTime);await wait(250);
  if(Number(cloudRoot.dataset.clouds)>0&&!matchMedia('(prefers-reduced-motion: reduce)').matches)assert(Number(cloudRoot.dataset.cloudTime)>cloudTime,'Illustrative cloud animates in the native webview');
  const l=document.getElementById('mq-render').getBoundingClientRect(),r=document.getElementById('mq-results').getBoundingClientRect();
  assert(r.left>=l.right,'Both panels must be side by side');
  document.getElementById('tab-compare').click();await until(()=>document.getElementById('research-spectrum').width>300);
  assert(document.getElementById('panel-explore').hidden,'Explore hides when comparing');
  document.getElementById('tab-optimize').click();assert(!document.getElementById('panel-optimize').hidden,'Optimize navigation');
  document.getElementById('tab-explore').click();await until(()=>document.getElementById('live-side').width>200);
  const fixture={snapshot:b,scene:await invoke('simulation_scene',{focus:0}),research:await invoke('research_status')};
  await until(()=>Number(document.getElementById('trap-live').dataset.ions)>0);
  document.getElementById('live-pause').click();await until(async()=>(await invoke('simulation_snapshot')).params.paused);await wait(150);
  const p=await invoke('simulation_snapshot');await wait(200);const q=await invoke('simulation_snapshot');assert(p.time_ms===q.time_ms,'Pause must stop ion integration');
  document.getElementById('live-step').click();await until(async()=>(await invoke('simulation_snapshot')).time_ms>q.time_ms);
  const f=document.getElementById('live-frequency-slider');f.value='82.5';f.dispatchEvent(new Event('input',{bubbles:true}));assert(document.getElementById('live-frequency').value==='82.5','Slider updates its precise readout');f.dispatchEvent(new Event('change',{bubbles:true}));
  await until(async()=>(await invoke('simulation_snapshot')).params.frequency===82.5);
  let rejected=false;try{await invoke('simulation_control',{controls:{rf:500}});}catch{rejected=true;}assert(rejected,'Invalid parameters must be rejected');
  await invoke('simulation_control',{controls:{paused:false}});await wait(300);
  const preset=document.getElementById('live-preset');preset.value='fig5';preset.dispatchEvent(new Event('change',{bubbles:true}));
  await until(async()=>{const s=await invoke('simulation_snapshot');return s.params.preset==='fig5'&&s.species.length===5&&s.warm===1;});
  assert((await invoke('simulation_snapshot')).ions.length>0,'Fig. 5 must produce live ions');
  await invoke('simulation_control',{controls:{preset:'s14'}});
  await until(async()=>{const s=await invoke('simulation_snapshot');return s.species.length===2&&s.warm===1;});
  await invoke('simulation_control',{controls:{preset:'s14',paused:true,frequency:82.5}});
  const combined=(await invoke('simulation_snapshot')).params;assert(combined.paused&&combined.frequency===82.5,'Preset plus additional controls preserve all supplied values');
  const study=(await invoke('research_status')).experiment;
  if(study.selected!=null){await invoke('research_apply',{index:study.selected});const p=(await invoke('simulation_snapshot')).params;assert(p.ac===study.confirmation.ac&&p.frequency===study.confirmation.frequency,'Preview must use independently confirmed settings');}
  await invoke('simulation_control',{controls:{preset:'s14'}});
  const exportPath=await invoke('research_export');
  assert(!document.getElementById('live-error').textContent,'No UI errors');
  await invoke('smoke_report',{success:true,message:`${b.ions.length} computed ions; S14 and Fig. 5, advancing time, pause, step, UI tuning, validation, contours, both panels and confirmed-candidate preview verified. Export: ${exportPath}`,fixture});
 }catch(e){await invoke('smoke_report',{success:false,message:String(e),fixture:null});}
})();
