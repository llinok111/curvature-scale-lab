// baseline/single-scale.mjs — the baseline we beat: classify curvature sign at
// the original resolution only (no scale-space). Also exposes uniform-area
// Gaussian curvature, the classic naive area weighting.
import { buildTopology, gaussianCurvature, classifySign, bboxDiag } from '../src/core.mjs';

export function singleScaleSigns(mesh, relThreshold = 1e-3) {
  const topo = buildTopology(mesh);
  const diag = bboxDiag(mesh.positions);
  const K = gaussianCurvature(topo, 'mixed');
  const Kn = new Float64Array(topo.n);
  for (let i = 0; i < topo.n; i++) Kn[i] = K[i] * diag * diag;
  return { topo, signs: classifySign(Kn, relThreshold), Kn };
}

export function uniformAreaKn(mesh) {
  const topo = buildTopology(mesh);
  const diag = bboxDiag(mesh.positions);
  const K = gaussianCurvature(topo, 'uniform');
  const out = new Float64Array(topo.n);
  for (let i = 0; i < topo.n; i++) out[i] = K[i] * diag * diag;
  return out;
}
