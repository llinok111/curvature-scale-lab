// src/shapes.mjs — deterministic analytic meshes with per-vertex ground truth.
// Pure: no three.js, no DOM, no I/O.

export function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const ICO_T = (1 + Math.sqrt(5)) / 2;
function icosahedron() {
  const v = [[-1, ICO_T, 0], [1, ICO_T, 0], [-1, -ICO_T, 0], [1, -ICO_T, 0],
    [0, -1, ICO_T], [0, 1, ICO_T], [0, -1, -ICO_T], [0, 1, -ICO_T],
    [ICO_T, 0, -1], [ICO_T, 0, 1], [-ICO_T, 0, -1], [-ICO_T, 0, 1]];
  const f = [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
    [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
    [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
    [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]];
  return { v, f };
}

// Icosphere: subdiv midpoints projected to radius r. Outward CCW winding.
export function icosphere(subdiv = 3, radius = 1) {
  const { v, f } = icosahedron();
  let verts = v.map((p) => { const l = Math.hypot(...p); return p.map((x) => x / l); });
  let faces = f.map((t) => t.slice());
  const cache = new Map();
  const mid = (a, b) => {
    const key = a < b ? `${a},${b}` : `${b},${a}`;
    if (cache.has(key)) return cache.get(key);
    const m = [(verts[a][0] + verts[b][0]) / 2, (verts[a][1] + verts[b][1]) / 2, (verts[a][2] + verts[b][2]) / 2];
    const l = Math.hypot(...m);
    verts.push([m[0] / l, m[1] / l, m[2] / l]);
    cache.set(key, verts.length - 1);
    return verts.length - 1;
  };
  for (let s = 0; s < subdiv; s++) {
    const nf = [];
    for (const [a, b, c] of faces) {
      const ab = mid(a, b), bc = mid(b, c), ca = mid(c, a);
      nf.push([a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]);
    }
    faces = nf;
    cache.clear();
  }
  const positions = [], indices = [], K = [], H = [];
  for (const p of verts) {
    positions.push(p[0] * radius, p[1] * radius, p[2] * radius);
    K.push(1 / (radius * radius));
    H.push(1 / radius);
  }
  for (const [a, b, c] of faces) indices.push(a, b, c);
  return { name: `icosphere-s${subdiv}`, positions, indices, analytic: { K, H }, closed: true };
}

export function noisySphere(subdiv = 4, radius = 1, amplitude = 0.02, seed = 20261117) {
  const base = icosphere(subdiv, radius);
  const rnd = mulberry32(seed);
  const positions = base.positions.slice();
  for (let i = 0; i < positions.length; i += 3) {
    const s = 1 + amplitude * (rnd() * 2 - 1);
    positions[i] *= s; positions[i + 1] *= s; positions[i + 2] *= s;
  }
  return { name: `noisy-sphere-${(amplitude * 100).toFixed(0)}pct`, positions, indices: base.indices.slice(), analytic: null, closed: true, note: 'ground truth = smooth sphere: K>0' };
}

// Open tube (no caps). K = 0, H = 1/(2r) on vertices away from the two rims.
export function tube(nu = 64, nz = 24, radius = 0.5, height = 1) {
  const positions = [], indices = [], K = [], H = [], interior = [];
  for (let j = 0; j <= nz; j++) {
    const z = -height / 2 + (height * j) / nz;
    for (let i = 0; i < nu; i++) {
      const u = (2 * Math.PI * i) / nu;
      positions.push(radius * Math.cos(u), z, radius * Math.sin(u));
      K.push(0); H.push(1 / (2 * radius));
      interior.push(j > 0 && j < nz);
    }
  }
  for (let j = 0; j < nz; j++) {
    for (let i = 0; i < nu; i++) {
      const i2 = (i + 1) % nu;
      const a = j * nu + i, b = j * nu + i2, c = (j + 1) * nu + i, d = (j + 1) * nu + i2;
      // du x dz points outward
      indices.push(a, c, b, b, c, d);
    }
  }
  return { name: `tube-nu${nu}`, positions, indices, analytic: { K, H }, interior, closed: false };
}

export function torus(nu = 64, nv = 32, R = 1, r = 0.35) {
  const positions = [], indices = [], K = [], H = [];
  for (let i = 0; i < nu; i++) {
    const u = (2 * Math.PI * i) / nu;
    for (let j = 0; j < nv; j++) {
      const v = (2 * Math.PI * j) / nv;
      const rr = R + r * Math.cos(v);
      positions.push(rr * Math.cos(u), r * Math.sin(v), rr * Math.sin(u));
      K.push(Math.cos(v) / (r * rr));
      H.push((R + 2 * r * Math.cos(v)) / (2 * r * rr));
    }
  }
  for (let i = 0; i < nu; i++) {
    for (let j = 0; j < nv; j++) {
      const i2 = (i + 1) % nu, j2 = (j + 1) % nv;
      const a = i * nv + j, b = i2 * nv + j, c = i * nv + j2, d = i2 * nv + j2;
      indices.push(a, b, c, b, d, c);
    }
  }
  return { name: `torus-${nu}x${nv}`, positions, indices, analytic: { K, H }, closed: true };
}

// Hyperbolic paraboloid z = a(x^2 - y^2); K < 0 everywhere, H = 0 at the origin.
export function saddle(n = 41, a = 0.4, half = 1) {
  const positions = [], indices = [], K = [], H = [];
  for (let j = 0; j < n; j++) {
    const y = -half + (2 * half * j) / (n - 1);
    for (let i = 0; i < n; i++) {
      const x = -half + (2 * half * i) / (n - 1);
      positions.push(x, y, a * (x * x - y * y));
      const s = 1 + 4 * a * a * (x * x + y * y);
      K.push((-4 * a * a) / (s * s));
      H.push((4 * a * a * a * (y * y - x * x)) / Math.pow(s, 1.5));
    }
  }
  for (let j = 0; j < n - 1; j++) {
    for (let i = 0; i < n - 1; i++) {
      const a0 = j * n + i, b0 = j * n + i + 1, c0 = (j + 1) * n + i, d0 = (j + 1) * n + i + 1;
      indices.push(a0, c0, b0, b0, c0, d0);
    }
  }
  return { name: `saddle-${n}x${n}`, positions, indices, analytic: { K, H }, closed: false };
}

export function sampleMeshes() {
  return [icosphere(4, 1), torus(64, 32, 1, 0.35), tube(64, 24, 0.5, 1), saddle(41, 0.4, 1), noisySphere(4, 1, 0.02, 20261117)];
}
