import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile, access } from 'node:fs/promises';
import { resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = fileURLToPath(new URL('../', import.meta.url));
test('desktop bundle contains both panels and all renderer scripts offline', async () => {
  const html = await readFile(resolve(root, 'dist/index.html'), 'utf8');
  assert.match(html, /id="mq-render"/);
  assert.match(html, /id="mq-results"/);
  assert.doesNotMatch(html, /window\.openai|https:\/\/cdn|<!-- INSTRUMENT -->/);
  const scripts = [...html.matchAll(/<script src="([^"]+)"[^>]*>/g)];
  assert.equal(scripts.length, 14);
  assert.match(html,/id="trap-live"/);
  assert.match(html,/id="research-start"/);
  assert.match(html, /id="live-side"/);
  assert.match(html, /id="live-terrain"/);
  for (const [, src] of scripts) await access(resolve(root, 'dist', src));
  const data = JSON.parse(html.match(/id="mq-recorded-data">([\s\S]*?)<\/script>/)[1]);
  assert.equal(data.rows.length, 6);
  assert.equal(data.rows.find(row => row.rel_AC === 1).R_f, 9.4);
});
