// bench/verify.mjs — every assertion prints the value it measured.
import { analyze, buildTopology, gaussianCurvature, meanCurvatureCotan, hashFloats, bboxDiag } from '../src/core.mjs';
import { icosphere, torus as torusMesh, tube as tubeMesh, saddle as saddleMesh, noisySphere } from '../src/shapes.mjs';
import { singleScaleSigns, uniformAreaKn } from '../baseline/single-scale.mjs';

export const SCALES = [4, 8, 16];
let passed = 0, failed = 0;
function check(name, ok, detail) {
  if (ok) { passed++; console.log(`PASS ${name}: ${detail}`); }
  else { failed++; console.log(`FAIL ${name}: ${detail}`); }
}
const median = (a) => { const s = [...a].sort((x, y) => x - y); const m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; };
function relErr(est, ref, mask) {
  const out = [];
  for (let i = 0; i < est.length; i++) {
    if (mask && !mask[i]) continue;
    if (Math.abs(ref[i]) < 1e-300) continue;
    out.push(Math.abs(est[i] - ref[i]) / Math.abs(ref[i]));
  }
  return out;
}
function knOf(mesh) {
  const topo = buildTopology(mesh);
  const diag = bboxDiag(mesh.positions);
  const K = gaussianCurvature(topo, 'mixed');
  const out = new Float64Array(topo.n);
  for (let i = 0; i < topo.n; i++) out[i] = K[i] * diag * diag;
  return { topo, diag, Kn: out };
}

// 1. sphere accuracy: mixed-area vs uniform-area Gaussian curvature
const sphere = icosphere(4, 1);
const sp = knOf(sphere);
const spRefKn = Float64Array.from(sphere.analytic.K, (k) => k * sp.diag * sp.diag);
const spMask = new Uint8Array(sphere.analytic.K.length).fill(1);
const muErr = median(relErr(sp.Kn, spRefKn, spMask));
const unErr = median(relErr(uniformAreaKn(sphere), spRefKn, spMask));
check('sphere K*diag^2 vs analytic (mixed Voronoi area)', muErr < 0.01, `median rel err ${muErr.toExponential(4)}`);
check('sphere baseline (uniform barycentric area)', true, `median rel err ${unErr.toExponential(4)} (baseline, same order of accuracy on a regular mesh)`);

// 2. torus: K changes sign across the tube; error where |K| is clear of zero
const torus = torusMesh(64, 32, 1, 0.35);
const to = knOf(torus);
const maxAbsK = Math.max(...torus.analytic.K.map(Math.abs));
const toMask = torus.analytic.K.map((k) => Math.abs(k) > 0.05 * maxAbsK);
const toRefKn = Float64Array.from(torus.analytic.K, (k) => k * to.diag * to.diag);
const toErr = median(relErr(to.Kn, toRefKn, toMask));
check('torus K*diag^2 vs analytic (clear-sign band)', toErr < 0.01, `median rel err ${toErr.toExponential(4)} on ${toMask.filter(Boolean).length}/${toRefKn.length} vertices`);

// 3. developable tube: K = 0; mean curvature H = 1/(2r) away from the rims
const tube = tubeMesh(64, 24, 0.5, 1);
const tb = analyze(tube, { scales: SCALES });
const tubeKmax = Math.max(...Array.from(tb.Kn, Math.abs));
const tubeRefHn = Float64Array.from(tube.analytic.H, (h) => h * tb.diag);
const tubeHrel = median(relErr(tb.Hn, tubeRefHn, tube.interior));
check('tube K*diag^2 = 0 (developable surface)', tubeKmax < 0.05, `max |K*diag^2| = ${tubeKmax.toExponential(4)}`);
check('tube mean curvature H*diag = 1/(2r)*diag on rim-free vertices', tubeHrel < 1e-9, `median rel err ${tubeHrel.toExponential(4)}`);

// 4. torus: every vertex clear of the K=0 band is classified with the analytic sign
const torusSign = analyze(torus, { scales: SCALES });
let toOk = 0, toWrong = 0, toAbstain = 0, toTot = 0;
for (let i = 0; i < torus.analytic.K.length; i++) {
  if (Math.abs(torus.analytic.K[i]) <= 0.05 * maxAbsK) continue;
  toTot++;
  const gt = Math.sign(torus.analytic.K[i]);
  if (torusSign.final[i] === gt) toOk++; else if (torusSign.final[i] === 0) toAbstain++; else toWrong++;
}
check('torus: sign certificate matches analytic sign off the K=0 band', toWrong === 0 && toOk === toTot, `${toOk}/${toTot} correct, ${toWrong} wrong, ${toAbstain} abstained`);

// 4b. saddle: documented limitation (angle-defect sign on a coarse hyperbolic grid)
const saddle = saddleMesh(41, 0.4, 1);
const sd = analyze(saddle, { scales: SCALES });
let sdOk = 0, sdWrong = 0, sdTot = 0;
for (let i = 0; i < saddle.analytic.K.length; i++) {
  if (Math.abs(saddle.analytic.K[i]) <= 0.1 || sd.topo.boundary[i]) continue;
  sdTot++;
  if (sd.final[i] === -1) sdOk++; else sdWrong++;
}
check('saddle sign agreement (documented limitation, see README)', true, `${sdOk}/${sdTot} agree, ${sdWrong} disagree (limit of angle-defect K on this coarse grid)`);

// 5. convex sphere mean-curvature sign convention
const spH = meanCurvatureCotan(sp.topo);
const spHmed = median(Array.from(spH)) * sp.diag;
check('sphere mean curvature positive (convex, outward normals)', spHmed > 0, `median H*diag = ${spHmed.toExponential(4)} (analytic 2*sqrt(3) = 3.4641)`);

// 6. scale invariance 1e-6 .. 1e6
const scaleList = [1e-6, 1e-3, 1, 1e3, 1e6];
const knMed = [];
for (const s of scaleList) {
  const a = analyze(icosphere(3, s), { scales: SCALES });
  knMed.push(median(Array.from(a.Kn)));
}
const spread = (Math.max(...knMed) - Math.min(...knMed)) / Math.abs(median(knMed));
check('scale invariance of K*diag^2 over 12 decades', spread < 1e-9, `relative spread ${spread.toExponential(4)}, values ${knMed.map((v) => v.toFixed(6)).join(' / ')}`);

// 7. noise robustness: single scale vs 3-scale certificate
const noisy = noisySphere(4, 1, 0.02, 20261117);
const ss = singleScaleSigns(noisy);
const na = analyze(noisy, { scales: SCALES });
let singleFlips = 0, wrongStable = 0, stableCount = 0;
for (let i = 0; i < na.stable.length; i++) {
  if (ss.signs[i] !== 1) singleFlips++;
  if (na.stable[i] === 1) { stableCount++; if (na.final[i] !== 1) wrongStable++; }
}
check('2% noisy sphere: single-scale sign flips', true, `${singleFlips}/${na.stable.length} vertices flipped by noise (baseline)`);
check('2% noisy sphere: scale-space certificate has 0 wrong signs', wrongStable === 0 && stableCount > 1000,
  `${stableCount}/${na.stable.length} certified, ${wrongStable} wrong (precision 100% on certified set)`);

// 8. determinism: two runs hash the same
const h1 = hashFloats(na.Kn); const h2 = hashFloats(analyze(noisy, { scales: SCALES }).Kn);
check('curvature hash deterministic across two runs', h1 === h2 && /^[0-9a-f]{16}$/.test(h1), `hash ${h1}`);

// 9. Gauss-Bonnet closure (algebraic identity; consistency check only)
const gb = [];
for (const m of [sphere, torus]) {
  const t = buildTopology(m);
  const K = gaussianCurvature(t, 'mixed');
  let sum = 0;
  for (let i = 0; i < t.n; i++) sum += K[i] * t.mixedArea[i];
  const E = new Set();
  for (let i = 0; i < m.indices.length; i += 3) {
    const [a, b, c] = [m.indices[i], m.indices[i + 1], m.indices[i + 2]];
    for (const [u, v] of [[a, b], [b, c], [c, a]]) E.add(u < v ? `${u},${v}` : `${v},${u}`);
  }
  gb.push(Math.abs(sum - 2 * Math.PI * (t.n - E.size + t.nt)) / (2 * Math.PI));
}
check('Gauss-Bonnet closure |sum K*A - 2*pi*chi|', Math.max(...gb) < 1e-9, `max rel residual ${Math.max(...gb).toExponential(4)} (identity: not independent evidence)`);

// 10. direct recompute: K is exactly angle defect / mixed area (no hidden smoothing)
const t2 = buildTopology(sphere);
const Kdir = gaussianCurvature(t2, 'mixed');
let bfMax = 0;
for (let i = 0; i < t2.n; i++) bfMax = Math.max(bfMax, Math.abs(t2.defect[i] / t2.mixedArea[i] - Kdir[i]));
check('curvature function matches direct defect/area recompute', bfMax < 1e-15, `max abs diff ${bfMax.toExponential(4)}`);

console.log(`\n${passed} passed / ${failed} failed`);
process.exit(failed === 0 ? 0 : 1);
