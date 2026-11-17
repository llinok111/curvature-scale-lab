// src/core.mjs — pure curvature analysis. No three.js, no DOM, no I/O.
// All functions operate on plain {positions:number[], indices:number[]} meshes.
// Curvature is reported in scale-normalised form: K * diag^2 (dimensionless).

const TWO_PI = Math.PI * 2;

export function bboxDiag(positions) {
  let mn = [Infinity, Infinity, Infinity];
  let mx = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < positions.length; i += 3) {
    for (let c = 0; c < 3; c++) {
      const v = positions[i + c];
      if (v < mn[c]) mn[c] = v;
      if (v > mx[c]) mx[c] = v;
    }
  }
  let s = 0;
  for (let c = 0; c < 3; c++) { const d = mx[c] - mn[c]; s += d * d; }
  return Math.sqrt(s);
}

// Vertex/edge/triangle topology with per-vertex areas and normals.
export function buildTopology(mesh) {
  const pos = mesh.positions;
  const idx = mesh.indices;
  const n = pos.length / 3;
  const nt = idx.length / 3;
  const nbSets = new Array(n);
  const vTris = new Array(n);
  for (let i = 0; i < n; i++) { nbSets[i] = new Set(); vTris[i] = []; }
  const triArea = new Float64Array(nt);
  const faceN = new Float64Array(nt * 3);
  const edgeKey = new Map();
  const cotW = new Map(); // edge -> sum of cot of opposite angles

  const p = (i) => [pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]];
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const len = (a) => Math.sqrt(dot(a, a));

  for (let t = 0; t < nt; t++) {
    const [ia, ib, ic] = [idx[3 * t], idx[3 * t + 1], idx[3 * t + 2]];
    const a = p(ia), b = p(ib), c = p(ic);
    const ab = sub(b, a), ac = sub(c, a);
    const cr = cross(ab, ac);
    const area = 0.5 * len(cr);
    triArea[t] = area;
    for (let k = 0; k < 3; k++) faceN[3 * t + k] = cr[k];
    vTris[ia].push(t); vTris[ib].push(t); vTris[ic].push(t);
    nbSets[ia].add(ib); nbSets[ia].add(ic);
    nbSets[ib].add(ia); nbSets[ib].add(ic);
    nbSets[ic].add(ia); nbSets[ic].add(ib);
    // edge -> boundary count
    for (const [u, v] of [[ia, ib], [ib, ic], [ic, ia]]) {
      const key = u < v ? `${u},${v}` : `${v},${u}`;
      edgeKey.set(key, (edgeKey.get(key) || 0) + 1);
    }
    // cotangent weights: for edge (i,j) the opposite angle is at k
    const pairs = [[ia, ib, ic], [ib, ic, ia], [ic, ia, ib]];
    for (const [u, v, w] of pairs) {
      const pu = p(u), pv = p(v), pw = p(w);
      const e1 = sub(pu, pw), e2 = sub(pv, pw);
      const cot = dot(e1, e2) / Math.max(len(cross(e1, e2)), 1e-300);
      const key = u < v ? `${u},${v}` : `${v},${u}`;
      cotW.set(key, (cotW.get(key) || 0) + cot);
    }
  }

  const neighbors = nbSets.map((s) => Int32Array.from(s));
  const boundary = new Uint8Array(n);
  for (const [key, count] of edgeKey) {
    if (count === 1) {
      const [u, v] = key.split(',').map(Number);
      boundary[u] = 1; boundary[v] = 1;
    }
  }

  // angle defect and mixed/uniform areas
  const defect = new Float64Array(n);
  defect.fill(0);
  const mixed = new Float64Array(n);
  const uniform = new Float64Array(n);
  for (let t = 0; t < nt; t++) {
    const ia = idx[3 * t], ib = idx[3 * t + 1], ic = idx[3 * t + 2];
    const a = p(ia), b = p(ib), c = p(ic);
    const area = triArea[t];
    const ang = (u, v, w) => { // angle at u between uv, uw
      const e1 = sub(v, u), e2 = sub(w, u);
      return Math.acos(Math.max(-1, Math.min(1, dot(e1, e2) / Math.max(len(e1) * len(e2), 1e-300))));
    };
    const A = ang(a, b, c), B = ang(b, c, a), C = ang(c, a, b);
    defect[ia] += A; defect[ib] += B; defect[ic] += C;
    uniform[ia] += area / 3; uniform[ib] += area / 3; uniform[ic] += area / 3;
    const la = len(sub(b, c)), lb = len(sub(c, a)), lc = len(sub(a, b));
    const obtuse = [A, B, C].findIndex((x) => x > Math.PI / 2);
    if (obtuse >= 0) {
      const vs = [ia, ib, ic];
      for (let k = 0; k < 3; k++) mixed[vs[k]] += (k === obtuse ? area / 2 : area / 4);
    } else {
      const cotA = 1 / Math.tan(A), cotB = 1 / Math.tan(B), cotC = 1 / Math.tan(C);
      mixed[ia] += (lb * lb * cotB + lc * lc * cotC) / 8;
      mixed[ib] += (lc * lc * cotC + la * la * cotA) / 8;
      mixed[ic] += (la * la * cotA + lb * lb * cotB) / 8;
    }
  }
  for (let i = 0; i < n; i++) {
    defect[i] = (boundary[i] ? Math.PI : TWO_PI) - defect[i];
    if (mixed[i] <= 0) mixed[i] = Math.max(uniform[i], 1e-300);
    if (uniform[i] <= 0) uniform[i] = 1e-300;
  }

  const normals = new Float64Array(n * 3);
  for (let t = 0; t < nt; t++) {
    for (let k = 0; k < 3; k++) {
      const v = idx[3 * t + k];
      for (let c = 0; c < 3; c++) normals[3 * v + c] += faceN[3 * t + c];
    }
  }
  for (let i = 0; i < n; i++) {
    const x = normals[3 * i], y = normals[3 * i + 1], z = normals[3 * i + 2];
    const l = Math.sqrt(x * x + y * y + z * z) || 1;
    normals[3 * i] = x / l; normals[3 * i + 1] = y / l; normals[3 * i + 2] = z / l;
  }

  return { n, nt, neighbors, vTris, triArea, mixedArea: mixed, uniformArea: uniform, defect, normals, cotW, boundary, positions: pos, indices: idx };
}

// Gaussian curvature (angle defect / area). weight: 'mixed' | 'uniform'.
export function gaussianCurvature(topo, weight = 'mixed') {
  const area = weight === 'uniform' ? topo.uniformArea : topo.mixedArea;
  const out = new Float64Array(topo.n);
  for (let i = 0; i < topo.n; i++) out[i] = topo.defect[i] / area[i];
  return out;
}

// Cotangent mean curvature (signed w.r.t. vertex normal). Positive = convex.
export function meanCurvatureCotan(topo) {
  const pos = topo.positions, n = topo.n;
  const out = new Float64Array(n);
  const key = (u, v) => (u < v ? `${u},${v}` : `${v},${u}`);
  // Accumulate sum w_ij (p_j - p_i) per vertex, reuse cotW.
  const acc = new Float64Array(n * 3);
  for (const [k, w] of topo.cotW) {
    const [i, j] = k.split(',').map(Number);
    for (let c = 0; c < 3; c++) {
      const d = pos[3 * j + c] - pos[3 * i + c];
      acc[3 * i + c] += w * d;
      acc[3 * j + c] -= w * d;
    }
  }
  for (let i = 0; i < n; i++) {
    const A = topo.mixedArea[i];
    const dx = acc[3 * i] / (2 * A), dy = acc[3 * i + 1] / (2 * A), dz = acc[3 * i + 2] / (2 * A);
    const nd = dx * topo.normals[3 * i] + dy * topo.normals[3 * i + 1] + dz * topo.normals[3 * i + 2];
    out[i] = -0.5 * nd;
  }
  return out;
}

// Uniform (umbrella) mean curvature, calibrated by mean incident edge length^2.
export function meanCurvatureUniform(topo) {
  const pos = topo.positions, n = topo.n;
  const out = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const nb = topo.neighbors[i];
    if (nb.length === 0) continue;
    let sx = 0, sy = 0, sz = 0, ls = 0;
    for (const j of nb) {
      sx += pos[3 * j] - pos[3 * i];
      sy += pos[3 * j + 1] - pos[3 * i + 1];
      sz += pos[3 * j + 2] - pos[3 * i + 2];
      const dx = pos[3 * j] - pos[3 * i], dy = pos[3 * j + 1] - pos[3 * i + 1], dz = pos[3 * j + 2] - pos[3 * i + 2];
      ls += Math.sqrt(dx * dx + dy * dy + dz * dz);
    }
    const k = nb.length;
    const h2 = (ls / k) * (ls / k);
    const mx = sx / k, my = sy / k, mz = sz / k;
    const nd = mx * topo.normals[3 * i] + my * topo.normals[3 * i + 1] + mz * topo.normals[3 * i + 2];
    out[i] = -0.5 * nd / h2;
  }
  return out;
}

// Scalar diffusion on the mesh graph: f <- (1-lambda) f + lambda * mean(neighbours).
export function diffuse(topo, field, steps, lambda = 0.5) {
  let cur = Float64Array.from(field);
  const n = topo.n, pos = topo.positions, nbAll = topo.neighbors;
  for (let s = 0; s < steps; s++) {
    const next = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const nb = nbAll[i];
      if (nb.length === 0) { next[i] = cur[i]; continue; }
      let sum = 0;
      for (let q = 0; q < nb.length; q++) sum += cur[nb[q]];
      next[i] = (1 - lambda) * cur[i] + lambda * (sum / nb.length);
    }
    cur = next;
  }
  return cur;
}

// sign with a dead-band: -1 / 0 / +1
export function classifySign(values, threshold) {
  const out = new Int8Array(values.length);
  for (let i = 0; i < values.length; i++) {
    const v = values[i];
    out[i] = v > threshold ? 1 : v < -threshold ? -1 : 0;
  }
  return out;
}

// Multi-scale stability: a vertex is stable if all scales agree on a non-zero sign.
export function stability(topo, field, scales, threshold) {
  const signs = scales.map((s) => classifySign(s === 0 ? field : diffuse(topo, field, s), threshold));
  const stable = new Int8Array(topo.n);
  const final = new Int8Array(topo.n);
  for (let i = 0; i < topo.n; i++) {
    let ok = 1;
    const s0 = signs[0][i];
    for (let k = 1; k < signs.length; k++) if (signs[k][i] !== s0 || s0 === 0) ok = 0;
    stable[i] = ok;
    final[i] = signs[signs.length - 1][i];
  }
  return { signs, stable, final };
}

export function signFlips(classification, analytic) {
  let flips = 0;
  for (let i = 0; i < classification.length; i++) {
    const a = analytic[i];
    if (a === 0) continue;
    if (classification[i] !== 0 && classification[i] !== Math.sign(a)) flips++;
    if (classification[i] === 0) flips++;
  }
  return flips;
}

// Deterministic 64-bit FNV-1a over the raw bytes of a Float64Array.
export function hashFloats(arr) {
  const bytes = new Uint8Array(arr.buffer, arr.byteOffset, arr.byteLength);
  let h = 0xcbf29ce484222325n;
  const prime = 0x100000001b3n;
  const mask = 0xffffffffffffffffn;
  for (let i = 0; i < bytes.length; i++) {
    h ^= BigInt(bytes[i]);
    h = (h * prime) & mask;
  }
  return h.toString(16).padStart(16, '0');
}

// Full analysis: dimensionless curvature K*diag^2, H*diag, scale-space stability.
export function analyze(mesh, opts = {}) {
  const scales = opts.scales || [0, 3, 6];
  const relThreshold = opts.relThreshold ?? 1e-3;
  const topo = buildTopology(mesh);
  const diag = bboxDiag(mesh.positions);
  const K = gaussianCurvature(topo, 'mixed');
  const H = meanCurvatureCotan(topo);
  const Kn = new Float64Array(topo.n), Hn = new Float64Array(topo.n);
  for (let i = 0; i < topo.n; i++) { Kn[i] = K[i] * diag * diag; Hn[i] = H[i] * diag; }
  const thr = relThreshold;
  const st = stability(topo, Kn, scales, thr);
  return { topo, diag, K, H, Kn, Hn, stable: st.stable, final: st.final, signs: st.signs, threshold: thr, scales };
}
