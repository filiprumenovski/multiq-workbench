import { readFile, writeFile, mkdir, copyFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import { resolve, dirname } from 'node:path';

const root = fileURLToPath(new URL('../', import.meta.url));
const dist = resolve(root, 'dist');
await mkdir(dist, { recursive: true });
let fragment = await readFile(resolve(root, 'src/instrument.html'), 'utf8');
fragment = fragment.replace('<!-- LIVE SIMULATOR -->', await readFile(resolve(root, 'src/live-panel.html'), 'utf8'));
fragment = fragment.replace('<!-- RESEARCH -->', await readFile(resolve(root, 'src/research-panel.html'), 'utf8'));
const styles = fragment.match(/<style>([\s\S]*?)<\/style>/);
if (!styles) throw new Error('Missing renderer styles');
await writeFile(resolve(dist, 'viewer.css'), styles[1]);
fragment = fragment.replace(styles[0], '');
const prefix = 'https://cdn.jsdelivr.net/npm/three@0.128.0/';
const files = [...fragment.matchAll(/<script src="([^"]+)"><\/script>/g)];
for (const [, url] of files) {
  if (!url.startsWith(prefix)) throw new Error(`Unexpected script source: ${url}`);
  const relative = url.slice(prefix.length);
  const destination = resolve(dist, 'vendor/three', relative);
  await mkdir(dirname(destination), { recursive: true });
  await copyFile(resolve(root, 'node_modules/three', relative), destination);
  fragment = fragment.replace(url, `vendor/three/${relative}`);
}
await copyFile(resolve(root, 'node_modules/three/LICENSE'), resolve(dist, 'vendor/three/LICENSE'));
const template = await readFile(resolve(root, 'src/index.html'), 'utf8');
if (!template.includes('<!-- INSTRUMENT -->')) throw new Error('Missing instrument mount');
await writeFile(resolve(dist, 'index.html'), template.replace('<!-- INSTRUMENT -->', fragment));
await copyFile(resolve(root, 'src/app.css'), resolve(dist, 'app.css'));
await copyFile(resolve(root, 'src/live.js'), resolve(dist, 'live.js'));
for (const name of ['research.js','trap.js','layout.js']) await copyFile(resolve(root,'src',name),resolve(dist,name));
console.log(`Built offline frontend with ${files.length} local renderer scripts.`);
