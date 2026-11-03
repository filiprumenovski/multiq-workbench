(() => {
'use strict';
const $=id=>document.getElementById(id), invoke=window.__TAURI__?.core?.invoke;
if(!invoke){$('live-status').textContent='Open MultiQ Workbench to run the live Rust simulator.';return;}
let timer=null,polling=false,renderKey='',terrainKey='',snapshot=null, scene=null, focus=0, sceneKey='', sceneBusy=false, controlsBusy=false, queuedControls=null, disposed=false;
let background=document.createElement('canvas'), backgroundDirty=true, palette;
const side=$('live-side'), end=$('live-end'), terrain=$('live-terrain');
function colors(){
 const s=getComputedStyle(document.documentElement);
 // Resolve light-dark() as an actual color for Canvas, just as the instrument renderer does.
 const probe=document.createElement('span');document.body.append(probe);
 const read=name=>{probe.style.color=`var(${name})`;return getComputedStyle(probe).color;};
 palette={bg:read('--background'),fg:read('--foreground'),muted:read('--muted-foreground'),border:read('--border'),solid:read('--card'),hot:read('--orange'),species:['#58a5c9','#c89b62','#a98bd0','#70af92','#d27b88']};probe.remove();backgroundDirty=true;
}
function fit(canvas){const r=canvas.getBoundingClientRect(),d=Math.min(devicePixelRatio,1.25);const w=Math.round(r.width*d),h=Math.round(r.height*d);if(canvas.width!==w||canvas.height!==h){canvas.width=w;canvas.height=h;backgroundDirty=true;}const c=canvas.getContext('2d');c.setTransform(d,0,0,d,0,0);return {c,w:r.width,h:r.height,d};}
function line(c,points,map,color,width=1,dash=[]){c.beginPath();let started=false;for(const p of points){if(p==null||p.some(v=>v==null||!Number.isFinite(v))){started=false;continue;}const [x,y]=map(p);if(!started){c.moveTo(x,y);started=true;}else c.lineTo(x,y);}c.strokeStyle=color;c.lineWidth=width;c.setLineDash(dash);c.stroke();c.setLineDash([]);}
function label(c,text,x,y,align='left'){c.font='11px -apple-system, sans-serif';c.fillStyle=palette.muted;c.textAlign=align;c.fillText(text,x,y);}
function mapping(w,h){const x0=40,y0=24,pw=w-54,ph=h-50;return {x0,y0,pw,ph,map:p=>[x0+(p[0]+.5)/21.2*pw,y0+(6.8-p[1])/13.6*ph]};}
function rects(c,items,map,kind){
 c.fillStyle=kind==='ghost'?palette.border:palette.solid;c.globalAlpha=kind==='ghost'?.25:1;
 for(const [i0,i1,k0,k1] of items){const a=map([i0-.5,k0-.5]),b=map([i1-.5,k1-.5]);c.fillRect(Math.min(a[0],b[0]),Math.min(a[1],b[1]),Math.abs(b[0]-a[0]),Math.abs(b[1]-a[1]));}c.globalAlpha=1;
}
function sideBackground(w,h,d,m){
 background.width=Math.round(w*d);background.height=Math.round(h*d);
 const c=background.getContext('2d');c.scale(d,d);c.fillStyle=palette.bg;c.fillRect(0,0,w,h);
 const grid=p=>m.map([p[1]*scene.h,p[0]*scene.h-scene.axis]);
 c.save();c.beginPath();c.rect(m.x0,m.y0,m.pw,m.ph);c.clip();
 for(let z=0;z<=20;z+=2)line(c,[[z,-7],[z,7]],m.map,palette.border,.5);
 line(c,[[-.5,0],[20.7,0]],m.map,palette.border,.7);
 rects(c,scene.ghost,grid,'ghost');
 for(const level of scene.contours){c.globalAlpha=level.index?.6:.25;for(const points of level.lines)line(c,points,grid,palette.muted,level.index?1:.55);}
 c.globalAlpha=1;rects(c,scene.cut,grid,'solid');
 for(const z of [.5,18.6])line(c,[[z,-7],[z,7]],m.map,palette.muted,1,[3,4]);
 c.restore();
 for(let z=0;z<=20;z+=5)label(c,String(z),m.map([z,0])[0],h-10,'center');
 for(const x of [-6,-3,0,3,6])label(c,String(x),30,m.map([0,x])[1]+4,'right');
 label(c,'x (mm)',4,12);label(c,'z (mm)',w-3,14,'right');
 for(const [z,name] of [[1,'Trap'],[9,'Q1 channel'],[16.8,'Gate']])label(c,name,m.map([z,0])[0],14,'center');
 backgroundDirty=false;
}
function drawSide(){
 const {c,w,h,d}=fit(side),m=mapping(w,h);if(!scene||!snapshot)return;
 if(backgroundDirty)sideBackground(w,h,d,m);
 c.clearRect(0,0,w,h);c.drawImage(background,0,0,w,h);
 c.save();c.beginPath();c.rect(m.x0,m.y0,m.pw,m.ph);c.clip();
 const map=p=>m.map([p[2],p[0]]),trails=$('live-trails').checked;
 for(let s=0;s<snapshot.species.length;s++){
  const color=palette.species[s%palette.species.length];
  if(trails){c.globalAlpha=.25;c.strokeStyle=color;c.lineWidth=.75;c.beginPath();for(const ion of snapshot.ions){if(ion.sp!==s)continue;const pts=[...ion.trail].reverse().concat([ion.p]);let first=true;for(const p of pts){const [x,y]=map(p);if(first){c.moveTo(x,y);first=false;}else c.lineTo(x,y);}}c.stroke();}
  c.globalAlpha=.9;c.fillStyle=color;c.beginPath();for(const ion of snapshot.ions){if(ion.sp!==s)continue;const [x,y]=map(ion.p);c.moveTo(x+1.2,y);c.arc(x,y,1.2,0,2*Math.PI);}c.fill();
 }
 c.fillStyle=palette.hot;c.beginPath();for(const ion of snapshot.ions){if(!ion.hot)continue;const [x,y]=map(ion.p);c.moveTo(x+1.6,y);c.arc(x,y,1.6,0,2*Math.PI);}c.fill();c.restore();
}
function drawEnd(){
 const {c,w,h}=fit(end);c.fillStyle=palette.bg;c.fillRect(0,0,w,h);if(!scene||!snapshot)return;
 const scale=(Math.min(w,h)-25)/14,map=p=>[w/2+p[0]*scale,h/2-p[1]*scale];
 rects(c,scene.end,p=>map([p[0]*scene.h-scene.axis,p[1]*scene.h-scene.axis]),'solid');
 line(c,[[-7,0],[7,0]],map,palette.border);line(c,[[0,-7],[0,7]],map,palette.border);
 c.strokeStyle=palette.hot;c.setLineDash([3,3]);c.beginPath();c.arc(w/2,h/2,1.2*scale,0,2*Math.PI);c.stroke();c.setLineDash([]);
 for(let s=0;s<snapshot.species.length;s++){c.fillStyle=palette.species[s%palette.species.length];c.globalAlpha=.65;c.beginPath();for(const ion of snapshot.ions){if(ion.sp!==s)continue;const [x,y]=map(ion.p);c.moveTo(x+1,y);c.arc(x,y,1,0,2*Math.PI);}c.fill();}c.globalAlpha=1;
 label(c,'±7 mm',w-4,h-4,'right');
}
let terrainRange=[0,120];
function drawTerrain(){
 const {c,w,h}=fit(terrain);c.fillStyle=palette.bg;c.fillRect(0,0,w,h);if(!scene||!snapshot)return;
 const t=scene.terrain,values=[snapshot.params.frequency,...t.fx.flat(),...t.fy.flat()].filter(v=>v!=null&&Number.isFinite(v));
 const hi=Math.max(120,Math.ceil(Math.max(...values)/20)*20);terrainRange=[0,hi];
 const map=p=>[40+(p[0]+.5)/21.2*(w-50),h-25-p[1]/hi*(h-45)];
 for(let f=0;f<=hi;f+=hi/3){line(c,[[-.5,f],[20.7,f]],map,palette.border,.5);label(c,f.toFixed(0),34,map([0,f])[1]+3,'right');}
 for(let z=0;z<=20;z+=5)label(c,String(z),map([z,0])[0],h-8,'center');
 for(let s=0;s<t.fx.length;s++)for(const [m,dash]of [[t.fx,[]],[t.fy,[3,3]]])line(c,t.z.map((z,i)=>[z,m[s][i]]),map,palette.species[s%palette.species.length],s===focus?1.7:1,dash);
 line(c,[[-.5,snapshot.params.frequency],[20.7,snapshot.params.frequency]],map,palette.hot,1.5,[7,4]);
 label(c,'kHz',5,12);label(c,'z (mm)',w-2,12,'right');
}
function redraw(){if(!side.clientWidth||document.hidden||!document.hasFocus())return;if(!palette)colors();drawSide();drawEnd();const key=sceneKey+':'+snapshot?.params.frequency;if(key!==terrainKey){drawTerrain();terrainKey=key;}}
function invalidatePlots(){renderKey='';terrainKey='';backgroundDirty=true;redraw();}
const observer=new ResizeObserver(invalidatePlots);[side,end,terrain].forEach(c=>observer.observe(c));
matchMedia('(prefers-color-scheme: dark)').addEventListener('change',()=>{colors();invalidatePlots();});
function showError(error){$('live-error').hidden=false;$('live-error').textContent=String(error);}
const parameterKeys=['frequency','rf','ac','gate','pressure'];
function sliderValue(key,value){const slider=$('live-'+key+'-slider');slider.value=value;slider.style.setProperty('--fill',100*(slider.value-slider.min)/(slider.max-slider.min)+'%');slider.setAttribute('aria-valuetext',value+' '+slider.dataset.unit);}
function inputs(p){
 if(!controlsBusy){
  if(document.activeElement!==$('live-preset'))$('live-preset').value=p.preset;
  for(const key of parameterKeys){const number=$('live-'+key),slider=$('live-'+key+'-slider');if(document.activeElement!==number&&document.activeElement!==slider){number.value=p[key];sliderValue(key,p[key]);}}
 }
 $('live-pause').textContent=p.paused?'Resume':'Pause';$('live-recirculate').checked=p.recirculate;
}
async function control(changes){
 if(snapshot?.experiment_running)return;
 if(controlsBusy){if(Object.keys(changes).every(k=>parameterKeys.includes(k)))queuedControls={...queuedControls,...changes};return;}
 controlsBusy=true;
 try{const p=await invoke('simulation_control',{controls:changes});if(snapshot)snapshot.params=p;$('live-error').hidden=true;sceneKey='';}
 catch(e){showError(e);}
 finally{controlsBusy=false;if(queuedControls){const next=queuedControls;queuedControls=null;control(next);}else{if(snapshot)inputs(snapshot.params);refresh();}}
}
for(const key of parameterKeys){
 const number=$('live-'+key),slider=$('live-'+key+'-slider');sliderValue(key,number.value);
 number.addEventListener('change',()=>{if(!number.checkValidity()||number.value===''){number.reportValidity();return;}sliderValue(key,number.value);control({[key]:Number(number.value)});});
 slider.addEventListener('input',()=>{number.value=slider.value;sliderValue(key,slider.value);});
 // Commit on release (or a keyboard step), avoiding expensive field rebuilds during a drag.
 slider.addEventListener('change',()=>control({[key]:Number(slider.value)}));
}
$('live-preset').addEventListener('change',e=>{focus=0;control({preset:e.target.value});});
$('live-pause').addEventListener('click',()=>control({paused:!snapshot.params.paused}));
$('live-step').addEventListener('click',()=>control({action:'step'}));
$('live-reset').addEventListener('click',()=>control({action:'reset'}));
$('live-recirculate').addEventListener('change',e=>control({recirculate:e.target.checked}));
$('live-trails').addEventListener('change',invalidatePlots);
terrain.addEventListener('click',e=>{if(!snapshot)return;const r=terrain.getBoundingClientRect();const f=(r.height-25-(e.clientY-r.top))/(r.height-45)*terrainRange[1];control({frequency:Math.round(Math.max(20,Math.min(120,f))*10)/10});});
let speciesKey='',lastStats=0;
function speciesButtons(){
 const key=snapshot.species.map(s=>s.label).join('|')+focus;if(key===speciesKey)return;speciesKey=key;
 $('live-species').replaceChildren(...snapshot.species.map((s,i)=>{const b=document.createElement('button');b.className='btn';b.type='button';b.textContent=s.label+' · '+s.mz.toFixed(1);b.setAttribute('aria-pressed',String(i===focus));b.style.borderBottomColor=palette.species[i%palette.species.length];b.addEventListener('click',()=>{focus=i;sceneKey='';speciesKey='';speciesButtons();});return b;}));
}
function fmt(value,sigma,mult=1){return value==null||sigma==null?'accumulating':`${(value*mult).toFixed(2)} ± ${(sigma*mult).toFixed(2)}`;}
function stats(){
 if(snapshot.experiment_running){$('live-currents').textContent='Scan mode: ions are independent frequency cohorts, not one physical beam. Live current and paper comparisons resume after the scan.';return;}
 const table=document.createElement('table'),head=document.createElement('tr');
 for(const t of ['Species','Transmitted (pA)','Per pass (%)','T/(T+C) (%)']){const th=document.createElement('th');th.textContent=t;head.append(th);}table.append(head);
 for(const s of snapshot.species){const tr=document.createElement('tr');for(const [i,t]of [s.label,fmt(s.current.value,s.current.sigma),fmt(s.transmission.value,s.transmission.sigma,100),fmt(s.efficiency.value,s.efficiency.sigma,100)].entries()){const td=document.createElement('td');td.textContent=t;if(i===2||i===3){const note=document.createElement('div');note.className='text-muted';note.textContent=(i===2?s.transmission:s.efficiency).status;td.append(note);}tr.append(td);}table.append(tr);}
 $('live-currents').replaceChildren(table);
}
async function updateScene(){
 if(sceneBusy||!snapshot.ready||!snapshot.species.length||snapshot.rebuilding)return;
 const p=snapshot.params,key=JSON.stringify([p.preset,p.rf,p.gate,focus]);if(key===sceneKey)return;sceneBusy=true;
 try{const next=await invoke('simulation_scene',{focus});scene=next;sceneKey=key;backgroundDirty=true;renderKey='';terrainKey='';document.getElementById('chait-multiq').dataset.liveScene='ready';}
 catch(e){showError(e);}finally{sceneBusy=false;}
}
function refresh(){clearTimeout(timer);if(!polling)poll();}
window.addEventListener('focus',()=>{invalidatePlots();refresh();});
window.addEventListener('multiq-viewchange',()=>{invalidatePlots();refresh();});
async function poll(){
 if(disposed||polling)return;polling=true;
 const active=!document.hidden&&document.hasFocus(),positions=active&&document.getElementById('chait-multiq').dataset.workspace==='explore';
 try{
  const next=await invoke('simulation_snapshot',{positions});if(next.error)throw new Error(next.error);snapshot=next;
  $('live-controls').disabled=!snapshot.ready||snapshot.warm<1||snapshot.experiment_running;
  const status=!snapshot.ready?'Solving channel fields…':snapshot.warm<1?`Preparing ion beam · ${Math.round(snapshot.warm*100)}%`:snapshot.rebuilding?'Rebuilding gate field…':snapshot.experiment_running?`Scanning independent frequency cohorts · current point ${snapshot.effective_frequency.toFixed(2)} kHz`:`${snapshot.params.paused?'Paused':'Running'} · ${snapshot.n_flight.toLocaleString()} ions · ${snapshot.time_ms.toFixed(2)} ms simulated`;
  $('live-status').textContent=status;document.getElementById('chait-multiq').dataset.liveSeq=String(snapshot.seq);
  inputs(snapshot.params);if(!palette)colors();speciesButtons();if(positions)await updateScene();const key=snapshot.time_ms+':'+snapshot.params.frequency+':'+positions;if(positions&&key!==renderKey){redraw();renderKey=key;}
  window.dispatchEvent(new CustomEvent('multiq-frame',{detail:{snapshot,scene,positions}}));
  if(performance.now()-lastStats>500){stats();lastStats=performance.now();}
 }catch(e){showError(e);$('live-status').textContent='Simulation connection stopped';$('live-controls').disabled=true;}
 finally{polling=false;}
 timer=setTimeout(poll,!active||!positions?1000:snapshot?.params.paused?500:100);
}
window.addEventListener('beforeunload',()=>{disposed=true;clearTimeout(timer);observer.disconnect();});poll();
})();
