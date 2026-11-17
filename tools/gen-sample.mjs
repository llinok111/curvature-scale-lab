// tools/gen-sample.mjs — deterministic sample corpus into data/meshes.json
import { writeFileSync, mkdirSync } from 'node:fs';
import { sampleMeshes } from '../src/shapes.mjs';
import { analyze, hashFloats } from '../src/core.mjs';

const SCALES = [4, 8, 16];
const median = (a) => { const s = [...a].sort((x, y) => x - y); const m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; };

export function buildCorpus() {
  return sampleMeshes().map((m) => ({
    name: m.name,
    closed: m.closed,
    positions: m.positions,
    indices: m.indices,
    analytic: m.analytic,
    interior: m.interior || null,
    note: m.note || null,
  }));
}

export function corpusText() {
  return JSON.stringify({ generator: 'tools/gen-sample.mjs', seed: 20261117, meshes: buildCorpus() }, null, 0);
}

if (import.meta.url === `file://${process.argv[1]}`) {
  mkdirSync(new URL('../data/', import.meta.url), { recursive: true });
  const text = corpusText();
  writeFileSync(new URL('../data/meshes.json', import.meta.url), text);
  const meshes = buildCorpus();
  let verts = 0, tris = 0;
  for (const m of meshes) { verts += m.positions.length / 3; tris += m.indices.length / 3; }
  console.log(`wrote data/meshes.json: ${meshes.length} meshes, ${verts} vertices, ${tris} triangles, ${text.length} bytes`);

  // node-side reference hashes; the browser recomputes with the same core and compares.
  const reference = { scales: SCALES, hash: 'FNV-1a-64 over Float64 K*diag^2', meshes: {} };
  for (const m of meshes) {
    const a = analyze(m, { scales: SCALES });
    reference.meshes[m.name] = {
      hash: hashFloats(a.Kn),
      certified: Array.from(a.stable).reduce((s, v) => s + v, 0),
      vertices: a.topo.n,
      medianKn: median(Array.from(a.Kn)),
    };
  }
  writeFileSync(new URL('../data/reference.json', import.meta.url), JSON.stringify(reference, null, 1));
  console.log(`wrote data/reference.json: ${Object.keys(reference.meshes).length} node-side hashes`);
}
