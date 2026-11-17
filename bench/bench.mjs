// bench/bench.mjs — cold start + >=5 samples, median/p95/stddev -> bench/result.json
import { analyze, buildTopology, gaussianCurvature, hashFloats } from '../src/core.mjs';
import { icosphere } from '../src/shapes.mjs';
import { singleScaleSigns, uniformAreaKn } from '../baseline/single-scale.mjs';
import { writeFileSync } from 'node:fs';

const SCALES = [4, 8, 16];

const mesh = icosphere(5, 1); // 10,242 vertices / 20,480 triangles
const S = 7;
function stats(a) {
  const v = [...a].sort((x, y) => x - y);
  const med = (w) => (w.length % 2 ? w[w.length >> 1] : (w[(w.length >> 1) - 1] + w[w.length >> 1]) / 2);
  const mean = v.reduce((x, y) => x + y, 0) / v.length;
  const sd = Math.sqrt(v.reduce((x, y) => x + (y - mean) ** 2, 0) / (v.length - 1));
  return { min: v[0], median: med(v), p95: v[Math.min(v.length - 1, Math.floor(0.95 * v.length))], max: v[v.length - 1], stddev: sd, relStddev: sd / mean };
}

const jobs = {
  gaussian_curvature_only: () => { const t = buildTopology(mesh); return hashFloats(gaussianCurvature(t, 'mixed')); },
  single_scale_analysis: () => hashFloats(singleScaleSigns(mesh).Kn),
  scale_space_certificate: () => hashFloats(analyze(mesh, { scales: SCALES }).Kn),
  baseline_uniform_area: () => hashFloats(uniformAreaKn(mesh)),
};
const out = { mesh: { name: mesh.name, vertices: mesh.positions.length / 3, triangles: mesh.indices.length / 3 }, scales: SCALES, samples: S, node: process.version, ops: {} };
for (const [name, fn] of Object.entries(jobs)) {
  const coldStart = (() => { const t0 = performance.now(); const h = fn(); return { ms: performance.now() - t0, hash: h }; })();
  fn(); // warmup
  const times = [];
  let hash = null;
  for (let i = 0; i < S; i++) { const t0 = performance.now(); hash = fn(); times.push(performance.now() - t0); }
  out.ops[name] = { coldStartMs: coldStart.ms, ...stats(times), samples: times.map((x) => +x.toFixed(4)), hash };
}
out.ratio_certificate_over_single = out.ops.scale_space_certificate.median / out.ops.single_scale_analysis.median;
const serialized = JSON.stringify(out, null, 2);
writeFileSync(new URL('./result.json', import.meta.url), serialized + '\n');
console.log(serialized);
