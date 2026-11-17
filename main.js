// main.js — rendering + interaction only. All curvature maths comes from src/core.mjs.
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { analyze, hashFloats } from './src/core.mjs';

const SCALES = [4, 8, 16];
const view = document.getElementById('view');
const panel = {
  mesh: document.getElementById('mesh'), toggle: document.getElementById('toggle'),
  verts: document.getElementById('verts'), certified: document.getElementById('certified'),
  mediankn: document.getElementById('mediankn'), bhash: document.getElementById('bhash'),
  nhash: document.getElementById('nhash'), match: document.getElementById('match'), note: document.getElementById('note'),
};

const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
renderer.setSize(view.clientWidth, view.clientHeight);
view.appendChild(renderer.domElement);
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x0e1116);
const camera = new THREE.PerspectiveCamera(45, view.clientWidth / view.clientHeight, 0.01, 5000);
camera.position.set(2.4, 1.8, 3.2);
const controls = new OrbitControls(camera, renderer.domElement);
scene.add(new THREE.AmbientLight(0xffffff, 0.65));
const dir = new THREE.DirectionalLight(0xffffff, 1.1); dir.position.set(3, 5, 4); scene.add(dir);
addEventListener('resize', () => {
  camera.aspect = view.clientWidth / view.clientHeight; camera.updateProjectionMatrix();
  renderer.setSize(view.clientWidth, view.clientHeight);
});

const [corpus, reference] = await Promise.all([
  fetch('./data/meshes.json').then((r) => r.json()),
  fetch('./data/reference.json').then((r) => r.json()),
]);
let current = null, showUnstable = false;

function colorFor(v) { // diverging: red K>0, blue K<0, dark near 0
  const k = new THREE.Color();
  if (v >= 0) k.setRGB(0.15 + 0.85 * Math.min(1, v), 0.12 + 0.35 * Math.min(1, v), 0.16);
  else k.setRGB(0.16, 0.2 + 0.25 * Math.min(1, -v), 0.18 + 0.82 * Math.min(1, -v));
  return k;
}

function load(name) {
  if (current) { scene.remove(current.group); current.group.geometry.dispose(); current.group.material.dispose(); }
  const mesh = corpus.meshes.find((m) => m.name === name);
  const a = analyze(mesh, { scales: SCALES });
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(mesh.positions, 3));
  g.setIndex(mesh.indices);
  g.computeVertexNormals();
  const scale = Math.max(...Array.from(a.Kn, Math.abs)) || 1;
  const colors = new Float32Array(a.topo.n * 3);
  for (let i = 0; i < a.topo.n; i++) {
    const c = showUnstable && a.stable[i] === 0 ? new THREE.Color(0.45, 0.45, 0.45) : colorFor(a.Kn[i] / scale);
    colors[3 * i] = c.r; colors[3 * i + 1] = c.g; colors[3 * i + 2] = c.b;
  }
  g.setAttribute('color', new THREE.BufferAttribute(colors, 3));
  const mat = new THREE.MeshStandardMaterial({ vertexColors: true, flatShading: false, side: THREE.DoubleSide, roughness: 0.75 });
  const group = new THREE.Mesh(g, mat);
  scene.add(group);
  current = { group, a, name };
  const ref = reference.meshes[name];
  panel.verts.textContent = String(a.topo.n);
  panel.certified.textContent = `${Array.from(a.stable).reduce((s, v) => s + v, 0)} / ${a.topo.n}`;
  panel.mediankn.textContent = a.Kn.length ? median(Array.from(a.Kn)).toFixed(4) : '-';
  panel.bhash.textContent = hashFloats(a.Kn);
  panel.nhash.textContent = ref.hash;
  const ok = hashFloats(a.Kn) === ref.hash;
  panel.match.textContent = ok ? 'node == browser' : 'MISMATCH';
  panel.match.className = ok ? 'v ok' : 'v warn';
  panel.note.textContent = `threshold ${a.threshold}; mesh=${mesh.name}; geometry from tools/gen-sample.mjs (no mesh is bundled by hand)`;
}
function median(a) { const s = [...a].sort((x, y) => x - y); const m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; }

for (const m of corpus.meshes) {
  const o = document.createElement('option'); o.value = m.name; o.textContent = `${m.name} (${m.positions.length / 3} v)`;
  panel.mesh.appendChild(o);
}
panel.mesh.addEventListener('change', () => load(panel.mesh.value));
panel.toggle.addEventListener('click', () => { showUnstable = !showUnstable; panel.toggle.textContent = showUnstable ? 'hide unstable vertices' : 'show unstable vertices'; load(panel.mesh.value); });
load(corpus.meshes[0].name);

renderer.setAnimationLoop(() => { controls.update(); renderer.render(scene, camera); });
