(() => {
'use strict';
const host=document.getElementById('trap-live'),note=document.getElementById('trap-live-note');
if(!host||!window.THREE)return;
const scene=new THREE.Scene(),camera=new THREE.PerspectiveCamera(38,1,.03,200),renderer=new THREE.WebGLRenderer({antialias:true,alpha:true,powerPreference:'low-power'});
renderer.setPixelRatio(Math.min(devicePixelRatio,1.25));host.append(renderer.domElement);
const controls=new THREE.OrbitControls(camera,renderer.domElement);controls.enableDamping=false;controls.minDistance=3;controls.maxDistance=65;
scene.add(new THREE.HemisphereLight(0xffffff,0x8394a3,.65));const light=new THREE.DirectionalLight(0xffffff,.7);light.position.set(-8,14,20);scene.add(light);
const hardware=new THREE.Group();scene.add(hardware);
const positions=new Float32Array(12000*3),colors=new Float32Array(12000*3),geo=new THREE.BufferGeometry();geo.setAttribute('position',new THREE.BufferAttribute(positions,3).setUsage(THREE.DynamicDrawUsage));geo.setAttribute('color',new THREE.BufferAttribute(colors,3).setUsage(THREE.DynamicDrawUsage));geo.setDrawRange(0,0);
const points=new THREE.Points(geo,new THREE.PointsMaterial({size:.09,vertexColors:true,transparent:true,opacity:.95,depthWrite:false}));points.frustumCulled=false;scene.add(points);
const palette=['#58a5c9','#c89b62','#a98bd0','#70af92','#d27b88'].map(c=>new THREE.Color(c));
let frameKey='',snapshot=null,built=false,view='channel';
function draw(){if(!host.clientWidth||document.hidden||!document.hasFocus())return;renderer.render(scene,camera);}
function positionView(){if(view==='entrance'){controls.target.set(-7.7,0,0);camera.position.set(-12,4,7);}else{controls.target.set(0,0,0);camera.position.set(-19,17,25);}controls.update();draw();}
function rod(a,b,r,mat){const start=new THREE.Vector3(...a),end=new THREE.Vector3(...b),delta=end.clone().sub(start),mesh=new THREE.Mesh(new THREE.CylinderGeometry(r,r,delta.length(),24),mat);mesh.position.copy(start.add(end).multiplyScalar(.5));mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0,1,0),delta.normalize());hardware.add(mesh);}
function build(g){
 const mat=new THREE.MeshPhongMaterial({color:0x9ea9b1,transparent:true,opacity:.2,depthWrite:false,side:THREE.DoubleSide});
 const [l1,l2,l3,l4,l5,l6,l7,l8]=g.layers,c=g.center,r=g.diameter/2+g.grow*g.grid;
 // Coordinates map exactly: native (x,y,z) mm -> display (z-9,x,y). No invented full-cube paths.
 for(const sx of [-1,1])for(const sy of [-1,1]){
  rod([l1-9,sx*c,sy*(c-g.tilt*(l2-l1))],[l2-9,sx*c,sy*c],r,mat);
  for(const [a,b,rad] of [[l2,l3,g.stem/2],[l4,l5,g.stem/2],[l5,l6,r],[l6,l7,g.stem/2]])rod([a-9,sx*c,sy*c],[b-9,sx*c,sy*c],rad,mat);
 }
 const shape=new THREE.Shape();const a=g.axis;shape.moveTo(-a,-a);shape.lineTo(a,-a);shape.lineTo(a,a);shape.lineTo(-a,a);shape.closePath();const hole=new THREE.Path();hole.absarc(0,0,g.hole/2-g.grow*g.grid,0,Math.PI*2,true);shape.holes.push(hole);
 const gate=new THREE.Mesh(new THREE.ExtrudeGeometry(shape,{depth:l4-l3,bevelEnabled:false}),new THREE.MeshPhongMaterial({color:0xc89b62,transparent:true,opacity:.2,depthWrite:false,side:THREE.DoubleSide}));gate.rotation.y=Math.PI/2;gate.position.x=l3-9;hardware.add(gate);
 const boundary=new THREE.Mesh(new THREE.CircleGeometry(2,48),new THREE.MeshBasicMaterial({color:0x58a5c9,transparent:true,opacity:.1,side:THREE.DoubleSide,depthWrite:false}));boundary.rotation.y=Math.PI/2;boundary.position.x=-8;hardware.add(boundary);
 const box=new THREE.Box3(new THREE.Vector3(-8.5,-2.5,-2.5),new THREE.Vector3(-6,2.5,2.5));hardware.add(new THREE.Box3Helper(box,0x58a5c9));
 built=true;positionView();host.dataset.geometry='native';
}
function update(){if(!snapshot||!host.clientWidth||document.hidden||!document.hasFocus())return;const key=[snapshot.time_ms,snapshot.ions.length,snapshot.params.paused,snapshot.experiment_running].join(':');if(key===frameKey)return;frameKey=key;let n=0,entrance=0;for(const ion of snapshot.ions){if(n>=12000)break;const p=ion.p;positions.set([p[2]-9,p[0],p[1]],n*3);const c=palette[ion.sp%palette.length];colors.set([c.r,c.g,c.b],n*3);if(p[2]<3)entrance++;n++;}geo.setDrawRange(0,n);geo.attributes.position.needsUpdate=true;geo.attributes.color.needsUpdate=true;host.dataset.ions=String(n);host.dataset.entranceIons=String(entrance);host.dataset.seq=String(snapshot.seq);note.textContent=`${entrance.toLocaleString()} visible ions at the trap exit · ${n.toLocaleString()} in the computed volume · ${snapshot.experiment_running?'independent scan cohorts':snapshot.params.paused?'paused':snapshot.time_ms.toFixed(2)+' ms simulated'}`;draw();}
window.addEventListener('multiq-frame',e=>{if(e.detail.positions===false)return;snapshot=e.detail.snapshot;if(!built&&e.detail.scene?.geometry)build(e.detail.scene.geometry);update();});
document.getElementById('trap-view').addEventListener('change',e=>{view=e.target.value;positionView();});
controls.addEventListener('change',draw);window.addEventListener('focus',draw);
const resize=new ResizeObserver(()=>{const r=host.getBoundingClientRect();if(!r.width||!r.height)return;renderer.setSize(r.width,r.height,false);camera.aspect=r.width/r.height;camera.updateProjectionMatrix();draw();});resize.observe(host);positionView();
window.addEventListener('beforeunload',()=>{resize.disconnect();controls.dispose();scene.traverse(o=>{o.geometry?.dispose();o.material?.dispose();});renderer.dispose();});
})();
