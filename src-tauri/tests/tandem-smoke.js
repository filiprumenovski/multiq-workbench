(async()=>{
 const wait=ms=>new Promise(r=>setTimeout(r,ms)),invoke=window.__TAURI__.core.invoke;
 try{
 document.getElementById('tab-tandem').click();
 document.getElementById('tandem-source').value='inlet';
 document.getElementById('tandem-n').value=180;
 document.getElementById('tandem-ms').value=10;
 document.getElementById('tandem-run').click();
 const until=Date.now()+240000;let s;
 do{await wait(250);s=await invoke('tandem_status');}while((!s||s.running)&&Date.now()<until);
 const f=s?.frame;if(s?.status!=='Complete'||!f||f.ports.length!==9||f.ports.some(p=>!p.entered)||f.split+f.combine+f.channel+f.detected+f.lost+f.back!==180||!f.ports.some(p=>p.from_other)||!f.ports.some(p=>p.forwarded))throw Error(JSON.stringify(s).slice(0,500));
 await invoke('smoke_report',{success:true,message:'Tandem: '+JSON.stringify({ms:f.ms,split:f.split,channel:f.channel,combine:f.combine,lost:f.lost,detected:f.detected,ports:f.ports}),fixture:null});
 }catch(e){await invoke('smoke_report',{success:false,message:String(e),fixture:null});}
})();
