(() => {
'use strict';
const $=id=>document.getElementById(id),invoke=window.__TAURI__?.core?.invoke;if(!invoke)return;
let data=null,status=null,snapshot=null,rowsKey='',comparisonKey='',disposed=false;
const series=['#c89b62','#70af92','#a98bd0','#58a5c9'];
function error(e){$('research-error').hidden=false;$('research-error').textContent=String(e);}
function num(v,d=2){return Number.isFinite(v)?v.toFixed(d):'—';}
function csv(name){const lines=data.files[name].text.split(/\r?\n/).filter(l=>l&&!l.startsWith('#'));const keys=lines.shift().split(',');return lines.map(l=>Object.fromEntries(l.split(',').map((v,i)=>[keys[i],v])));}
function table(headers,rows){const t=document.createElement('table'),head=document.createElement('tr');for(const h of headers){const th=document.createElement('th');th.textContent=h;head.append(th);}t.append(head);for(const row of rows){const tr=document.createElement('tr');for(const v of row){const td=document.createElement('td');if(v instanceof Node)td.append(v);else td.textContent=v;tr.append(td);}t.append(tr);}return t;}
function canvas(id,xlo,xhi,yhi,xname){const el=$(id),r=el.getBoundingClientRect(),d=Math.min(devicePixelRatio,2);el.width=r.width*d;el.height=r.height*d;const c=el.getContext('2d');c.scale(d,d);const w=r.width,h=r.height;const fg=getComputedStyle(document.body).color;c.font='11px -apple-system,sans-serif';c.fillStyle=fg;c.strokeStyle=fg;c.globalAlpha=.4;c.beginPath();c.moveTo(40,15);c.lineTo(40,h-28);c.lineTo(w-10,h-28);c.stroke();c.globalAlpha=1;const map=(x,y)=>[40+(x-xlo)/(xhi-xlo)*(w-55),h-28-y/yhi*(h-48)];for(let i=0;i<=4;i++){const x=xlo+(xhi-xlo)*i/4;const [px]=map(x,0);c.textAlign='center';c.fillText(num(x,xhi-xlo>100?0:1),px,h-12);}c.textAlign='left';c.fillText(xname,w-60,12);c.fillText(num(yhi,0)+'%',2,15);c.fillText('0',20,h-27);return {c,map,w,h};}
function line(c,map,points,color){c.strokeStyle=color;c.lineWidth=1.4;c.beginPath();for(let i=0;i<points.length;i++){const p=map(...points[i]);if(i)c.lineTo(...p);else c.moveTo(...p);}c.stroke();}
function spectrum(){
 if(!data||!snapshot||!$('research-spectrum').clientWidth)return;
 const s14=snapshot.params.preset==='s14';let measured,xlo,xhi;
 if(s14){measured=csv('figS14_II_spectrum_profile.csv').map(r=>[+r.mz,+r.rel_abundance_pct]);xlo=840;xhi=900;}
 else {const f=snapshot.params.frequency/(snapshot.params.rf/56);const channels=[[9,59.25,'A1'],[8,53,'A2'],[7,47,'A3'],[6,39,'A4']].sort((a,b)=>Math.abs(a[1]-f)-Math.abs(b[1]-f));measured=csv('fig5A_profiles_native.csv').filter(r=>r.panel===channels[0][2]).map(r=>[+r.mz,+r.height_pct_of_axis100]);xlo=750;xhi=1500;}
 const {c,map,w,h}=canvas('research-spectrum',xlo,xhi,105,'m/z');c.save();c.beginPath();c.rect(40,16,w-50,h-44);c.clip();line(c,map,measured.filter(p=>p[0]>=xlo&&p[0]<=xhi),'#58a5c9');
 const weights=s14?[1,.0893]:[15.9,60.3,100,15.2,4.8];const amp=snapshot.species.map((s,i)=>(s.transmission.value||0)*weights[i]);const scale=s14?amp[0]:Math.max(...amp);
 if(scale>0&&!snapshot.experiment_running){c.strokeStyle='#c89b62';c.setLineDash([4,3]);snapshot.species.forEach((s,i)=>{const p=map(s.mz,100*amp[i]/scale),base=map(s.mz,0);c.beginPath();c.moveTo(...base);c.lineTo(...p);c.stroke();c.beginPath();c.arc(...p,3,0,Math.PI*2);c.stroke();});c.setLineDash([]);}c.restore();
 $('research-spectrum').dataset.ready='true';
 $('research-spectrum-note').textContent=s14?'S14 Part II · blue: digitized spectrum; amber: two computed envelope heights, normalized to the main ion. Isotope structure and detector response are not simulated here.':'Fig. 5A · blue: digitized selected-channel spectrum; amber: computed 6+–10+ peak heights, normalized to their own maximum. Chemical background and adduct envelopes are omitted. Use the ratio table for quantitative comparison.';
}
function curves(){
 const rows=status?.experiment.rows||[],el=$('research-curves');el.hidden=!rows.length;$('research-empty').hidden=!!rows.length;
 const legend=$('research-curve-legend');legend.replaceChildren();legend.hidden=!rows.length;
 rows.forEach((r,i)=>{const item=document.createElement('span'),swatch=document.createElement('i');swatch.style.background=series[i];item.append(swatch,document.createTextNode('±'+num(r.ac,4)+' V'));legend.append(item);});
 if(!rows.length||!el.clientWidth)return;
 const ymax=Math.max(10,...rows.flatMap(r=>[...(r.coarse?.points||[]),...r.points].map(p=>p.species[0].t*100)));
 const {c,map}=canvas('research-curves',70,100,Math.ceil(ymax*1.15/5)*5,'kHz');
 rows.forEach((r,i)=>{const p=r.coarse?.points||r.points;line(c,map,p.map(p=>[p.f,100*p.species[0].t]),series[i]);if(r.coarse)line(c,map,r.points.map(p=>[p.f,100*p.species[0].t]),series[i]);});
}
function outcome(ex){
 const el=$('research-outcome'),r=ex.confirmation;el.hidden=ex.running||ex.selected==null||!r?.fit||!!r.issue;el.replaceChildren();if(el.hidden)return;
 const intro=document.createElement('div'),title=document.createElement('strong'),caption=document.createElement('p');title.textContent='Confirmed in model';caption.textContent='Independent scan · sampling errors shown';intro.append(title,caption);el.append(intro);
 for(const [label,value,unit] of [['Mass resolving power',num(r.fit.r_mass,1),'± '+num(r.fit.r_mass_error,1)],['Target transmission',num(r.transmission*100,1)+'%','± '+num(r.transmission_error*100,1)+'%'],['Candidate drive','±'+num(r.ac,4)+' V',num(r.frequency,3)+' kHz']]){const cell=document.createElement('div'),heading=document.createElement('div'),metric=document.createElement('div'),detail=document.createElement('small');heading.className='outcome-label';heading.textContent=label;metric.className='outcome-value';metric.append(document.createTextNode(value+' '));detail.textContent=unit;metric.append(detail);cell.append(heading,metric);el.append(cell);}
}
function comparisons(){const key=JSON.stringify([status.comparisons,status.experiment.running,status.reference_settings]);if(key===comparisonKey)return;comparisonKey=key;
 $('research-comparison').replaceChildren(table(['Quantity','Their measurement','Live model'],status.comparisons.filter(r=>!r.name.startsWith('resonance')).map(r=>{const name=document.createElement('span');name.textContent=r.name;name.title=r.source;return[name,r.published,status.experiment.running?'Comparison paused during scan':r.model+(r.settle?' · '+r.settle:'')];})));
 $('research-context').textContent=(status.reference_settings?'Reference model settings · ':'Exploratory settings · ')+status.context;
}
function results(){const ex=status.experiment,key=JSON.stringify([ex.rows,ex.selected,ex.confirmation,ex.running]);if(key===rowsKey)return;rowsKey=key;
 $('research-results').replaceChildren(table(['AC ±V','R(m/z) ±σ','Per pass','Na/main T','Fit / action'],(ex.rows||[]).map((r,i)=>{const actions=document.createElement('div');const note=document.createElement('div');note.textContent=r.issue||`χ²/dof ${num(r.fit?.chi2,1)} · ${r.stage}`;actions.append(note);if(r.frequency!=null){const b=document.createElement('button');b.className='btn';b.textContent=ex.selected===i?'Preview confirmed':'Preview';b.disabled=ex.running;b.onclick=()=>act('research_apply',{index:i});actions.append(b);}return[num(r.ac,4),r.fit?`${num(r.fit.r_mass,1)} ± ${num(r.fit.r_mass_error,1)}`:'unresolved',num(100*r.transmission,1)+'%',r.leakage_count===0?`< ${num(r.leakage_upper95,3)} (95%)`:`${num(r.leakage)} ± ${num(r.leakage_error)}`,actions];})));
 if(ex.selected!=null)$('research-results').querySelectorAll('tr')[ex.selected+1]?.classList.add('confirmed-row');outcome(ex);
 const baseline=ex.rows?.[0],p=baseline?.points.find(p=>Math.abs(p.f-83)<.01);let note='';
 if(p){const [main,na]=p.species;note=`Baseline at the paper’s 83 kHz: model Na/main transmission ${num(na.t/main.t)}; published ${num(.64)} ± ${num(.186,3)} (digitization error). `;}
 if(baseline?.fit)note+=`Baseline resonance ${num(baseline.fit.f)} ± ${num(baseline.fit.f_error)} kHz; paper 83 ± 2 kHz reading range. `;
 if(ex.confirmation){const r=ex.confirmation;note+=`Confirmation: R(m/z) ${num(r.fit?.r_mass,1)} ± ${num(r.fit?.r_mass_error,1)}, transmitted ${num(r.transmission*100,1)} ± ${num(r.transmission_error*100,1)}%. `;}
 $('research-result-note').textContent=note+'Uncertainties here describe ion sampling and curve fitting, not geometry, calibration or whole-trap physics. All completed curves and counts are retained in the export.';curves();
}
async function act(name,args={}){try{$('research-error').hidden=true;const result=await invoke(name,args);if(name==='research_export')$('research-export-note').textContent='Saved '+result;else {await pollOnce();if(name==='research_apply')window.dispatchEvent(new CustomEvent('multiq-navigate',{detail:'explore'}));}return result;}catch(e){error(e);}}
$('research-start').onclick=()=>{if(!$('research-resolution').reportValidity())return;act('research_start',{minR:Number($('research-resolution').value),ions:Number($('research-ions').value)});};
$('research-cancel').onclick=()=>act('research_cancel');$('research-export').onclick=()=>act('research_export');$('research-reference').onclick=()=>act('simulation_control',{controls:{preset:'s14'}});
async function pollOnce(){status=await invoke('research_status');const ex=status.experiment;$('research-start').disabled=!!ex.running||!snapshot?.ready;$('research-cancel').disabled=!ex.running;$('research-reference').disabled=!!ex.running;$('research-resolution').disabled=!!ex.running;$('research-ions').disabled=!!ex.running;$('research-status').textContent=ex.message;$('research-progress').hidden=!ex.running;$('research-progress').value=ex.progress||0;if(ex.storage_error)error(ex.storage_error);comparisons();results();spectrum();window.dispatchEvent(new CustomEvent('multiq-research',{detail:ex}));}
window.addEventListener('multiq-frame',e=>{const first=!snapshot;snapshot=e.detail.snapshot;if(first&&data)spectrum();});
const resize=new ResizeObserver(()=>{spectrum();curves();});resize.observe($('research-spectrum'));resize.observe($('research-curves'));window.addEventListener('multiq-viewchange',()=>{spectrum();curves();});matchMedia('(prefers-color-scheme: dark)').addEventListener('change',()=>{spectrum();curves();});
async function poll(){if(disposed)return;try{await pollOnce();}catch(e){error(e);}setTimeout(poll,document.hidden||!document.hasFocus()?5000:1500);}
(async()=>{try{data=await invoke('research_data');$('research-sources').textContent=`${Object.keys(data.files).length} source files bundled with SHA-256 checksums. Sources: ${data.source}; SPECTRA_EXTRACTION.md. Export includes the original CSV text, extraction notes, parameters, seeds and raw counts. Simulator ${data.simulator_commit.slice(0,9)}${data.simulator_dirty?' (research checkout has local changes)':''}.`;poll();}catch(e){error(e);}})();
window.addEventListener('beforeunload',()=>{disposed=true;resize.disconnect();});
})();
