(() => {
'use strict';
const root=document.getElementById('chait-multiq'),tabs=[...document.querySelectorAll('[role=tab][data-view]')];
function show(view,focus=false){
 const selected=tabs.find(t=>t.dataset.view===view);if(!selected)return;
 for(const tab of tabs){const active=tab===selected;tab.setAttribute('aria-selected',String(active));tab.tabIndex=active?0:-1;document.getElementById(tab.getAttribute('aria-controls')).hidden=!active;}
 root.dataset.workspace=view;if(focus)selected.focus();
 window.dispatchEvent(new CustomEvent('multiq-viewchange',{detail:view}));
}
for(const tab of tabs){tab.addEventListener('click',()=>show(tab.dataset.view));tab.addEventListener('keydown',e=>{const i=tabs.indexOf(tab);let next;if(e.key==='ArrowRight')next=(i+1)%tabs.length;if(e.key==='ArrowLeft')next=(i+tabs.length-1)%tabs.length;if(e.key==='Home')next=0;if(e.key==='End')next=tabs.length-1;if(next!=null){e.preventDefault();show(tabs[next].dataset.view,true);}});}
for(const button of document.querySelectorAll('[data-open-view]'))button.addEventListener('click',()=>show(button.dataset.openView,true));
window.addEventListener('multiq-navigate',e=>{show(e.detail,true);window.scrollTo({top:0,behavior:'instant'});});
let last=0;
window.addEventListener('multiq-frame',e=>{if(performance.now()-last<300)return;last=performance.now();const s=e.detail.snapshot,el=document.getElementById('app-run-state');const state=s.experiment_running?'scan':!s.ready||s.warm<1?'loading':s.params.paused?'paused':'running';const label={scan:'Study running',loading:'Preparing',paused:'Paused',running:'Live simulation'}[state];if(el.dataset.state!==state){el.dataset.state=state;el.replaceChildren(Object.assign(document.createElement('i'),{ariaHidden:'true'}),document.createTextNode(label));}document.getElementById('settings-lock').hidden=!s.experiment_running;});
window.addEventListener('multiq-research',e=>{const ex=e.detail,badge=document.getElementById('study-count');badge.hidden=!ex.rows?.length;badge.textContent=String(ex.rows?.length||0);});
show('explore');
})();

(()=>{
 const el=id=>document.getElementById(id), ipc=(name,args)=>window.__TAURI__.core.invoke(name,args);
 const tones=[39,0,47,0,0,0,53,0,59.25],inputs=tones.map((hz,i)=>{
  const l=document.createElement('label'),n=document.createElement('input');l.textContent='P'+(i+1);n.type='number';n.min=0;n.max=120;n.step='any';n.required=true;n.value=hz;l.append(n);el('tandem-tones').append(l);return n;
 });
 let timer,visible=false;
 function render(s){
  if(!s)return;const f=s.frame;
  el('tandem-controls').disabled=el('tandem-run').disabled=!!s.running;el('tandem-stop').disabled=!s.running;el('tandem-save').disabled=!f||s.running;
  el('tandem-status').textContent=s.status+(f?' · '+f.ms.toFixed(3)+' ms':'');
  if(!f){el('tandem-counts').textContent='';el('tandem-table').tBodies[0].replaceChildren();el('tandem-map').getContext('2d').clearRect(0,0,880,260);return}
  el('tandem-counts').textContent='Splitter '+f.split+' · In channels '+f.channel+' · Combiner '+f.combine+' · Collected '+f.detected+' · Source returns '+f.back+' · Losses '+f.lost;
  const body=el('tandem-table').tBodies[0];body.replaceChildren();
  f.ports.forEach((p,i)=>{const r=body.insertRow();for(const v of [i+1,s.settings.tones_khz[i]||'Blocked',p.entered,p.returned,p.from_other,p.forwarded,p.lost,p.by_charge.join(" / ")])r.insertCell().textContent=v;});
  const ctx=el('tandem-map').getContext('2d'),size=f.size_mm,seam=f.seam_mm,gap=18.1,scale=210/size,offset=35;
  ctx.clearRect(0,0,880,260);ctx.strokeStyle='#686877';ctx.fillStyle='#d6d6e3';ctx.font='13px system-ui';
  for(const [z,text] of [[0,'Splitter'],[seam+gap,'Combiner']]){ctx.strokeRect(offset+z*scale,24,seam*scale,size*scale);ctx.fillText(text,offset+z*scale,17);}
  const pitch=8.3;for(const k of [-1,0,1]){const y=24+scale*(size/2+k*2*pitch);ctx.strokeRect(offset+seam*scale,y-5,gap*scale,10);}
  ctx.fillText('x–z projection',offset,252);const faceX=555;ctx.strokeRect(faceX,24,size*scale,size*scale);ctx.fillText('End view · all 9 ports',faceX,17);for(let i=0;i<9;i++){const px=faceX+(size/2+(i%3-1)*2*pitch)*scale,py=24+(size/2+(Math.floor(i/3)-1)*2*pitch)*scale;ctx.strokeStyle=s.settings.tones_khz[i]?'#a391ff':'#535360';ctx.strokeRect(px-9,py-9,18,18);ctx.fillText(String(i+1),px+10,py-8);}for(const [x,z,sp,,y] of f.points){ctx.fillStyle=['#aa99ff','#54c7ff','#56dba5','#f2ca70','#ef80a5'][sp];ctx.fillRect(offset+z*scale,24+x*scale,3,3);ctx.fillRect(faceX+x*scale,24+y*scale,2,2);}
 }
 async function poll(){clearTimeout(timer);if(!visible)return;try{render(await ipc('tandem_status'));}catch(e){el('tandem-status').textContent=String(e)} if(visible)timer=setTimeout(poll,500);}
 window.addEventListener('multiq-viewchange',e=>{visible=e.detail==='tandem';if(visible)poll();else clearTimeout(timer);});
 el('tandem-form').onsubmit=async e=>{e.preventDefault();try{await ipc('tandem_start',{tones:inputs.map(n=>n.valueAsNumber),ions:el('tandem-n').valueAsNumber,ms:el('tandem-ms').valueAsNumber,seed:el('tandem-seed').valueAsNumber,fine:el('tandem-grid').value==='fine',mouth:el('tandem-source').value==='mouth'});await poll();}catch(err){el('tandem-status').textContent=String(err);}};
 el('tandem-stop').onclick=()=>ipc('tandem_stop');
 el('tandem-save').onclick=async()=>{try{const path=await ipc('tandem_export');el('research-export-note').textContent='Saved '+path;}catch(e){el('tandem-status').textContent=String(e)}};
 for(const [id,channel] of [['original-trap',false],['original-channel',true]])el(id).onclick=async()=>{try{await ipc('original_sketch',{channel})}catch(e){el('tandem-status').textContent=String(e)}};
})();
