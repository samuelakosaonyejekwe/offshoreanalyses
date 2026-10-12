import test from 'node:test';
import assert from 'node:assert/strict';
import {
  mat4Identity, mat4Multiply, mat4Perspective, mat4Ortho, mat4FromQuat, projectPoint,
  quatMultiply, quatFromAxisAngle, quatRotate, quatFromBasis, quatNormalize,
  arcballVector, arcballRotation, viewOrientation, defaultCamera, cameraAxes, cameraMatrices, fitDistance, fitCamera,
  orbitArcball, orbitTurntable, panTarget, zoomCamera, screenRay,
  colourMap, colourMapTexture, colourMapStops, luminance, valueRange, normaliseValues, parseColour, formatNumber,
  vertexNormals, triEdgePairs, uniqueEdges, extractBoundary, decimateIndices, autoExaggeration, gridMesh,
  buildItem, applyDeformed, itemEdges, itemBounds, unionBounds, thickLineBuffers, expandIndexed, highlightGeometry,
  pickTriangles, pickScreen, dimensionText, describeScene, shaderSource,
} from '../../src/ui/viewer3d.js';

const close = (a, b, tol = 1e-9, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} expected ${b}, got ${a}`);
const closeV = (a, b, tol = 1e-9) => { assert.equal(a.length, b.length); a.forEach((v, i) => close(v, b[i], tol, `component ${i}`)); };

const cube = () => ({
  kind: 'surface', name: 'cube',
  positions: Float32Array.from([0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1]),
  indices: Uint32Array.from([0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6, 3, 0, 4, 3, 4, 7]),
});
const hexCell = (ox = 0) => [0, 1, 2, 3, 4, 5, 6, 7].map((i) => i + ox);
const twoHexes = () => {
  // two unit cubes sharing the face x = 1
  const P = [];
  for (const z of [0, 1]) for (const [x, y] of [[0, 0], [1, 0], [1, 1], [0, 1]]) P.push(x, y, z);
  for (const z of [0, 1]) for (const [x, y] of [[2, 0], [2, 1]]) P.push(x, y, z);
  // nodes 8: (2,0,0) 9: (2,1,0) 10: (2,0,1) 11: (2,1,1)
  return { kind: 'volume', name: 'two', positions: Float32Array.from(P), cells: [{ type: 'hex8', conn: Uint32Array.from([...hexCell(), 1, 8, 9, 2, 5, 10, 11, 6]) }] };
};

test('matrix product: identity is neutral and the product is associative', () => {
  const a = mat4Perspective(0.7, 1.5, 0.1, 100), b = mat4FromQuat(quatFromAxisAngle([1, 2, 3], 0.8)), c = mat4Ortho(3, 2, -1, 9);
  closeV(Array.from(mat4Multiply(a, mat4Identity())), Array.from(a));
  closeV(Array.from(mat4Multiply(mat4Multiply(a, b), c)), Array.from(mat4Multiply(a, mat4Multiply(b, c))), 1e-12);
});

test('perspective projection maps the near and far planes to z = −1 and +1 and the frustum edge to x, y = ±1', () => {
  const fov = Math.PI / 3, aspect = 2, n = 0.5, f = 80;
  const m = mat4Perspective(fov, aspect, n, f);
  close(projectPoint(m, 0, 0, -n)[2], -1, 1e-12);
  close(projectPoint(m, 0, 0, -f)[2], 1, 1e-12);
  const t = Math.tan(fov / 2) * 10;
  const p = projectPoint(m, t * aspect, t, -10);
  close(p[0], 1, 1e-12); close(p[1], 1, 1e-12); close(p[3], 10, 1e-12);
  assert.throws(() => mat4Perspective(1, 1, 0, 10), RangeError);
});

test('orthographic projection maps the view box to the unit cube', () => {
  const m = mat4Ortho(4, 1.5, 2, 12);
  closeV(projectPoint(m, 6, 4, -2), [1, 1, -1, 1], 1e-12);
  closeV(projectPoint(m, -6, -4, -12), [-1, -1, 1, 1], 1e-12);
});

test('quaternions: axis-angle rotation, composition order and matrix agree', () => {
  const qz = quatFromAxisAngle([0, 0, 1], Math.PI / 2);
  closeV(quatRotate(qz, [1, 0, 0]), [0, 1, 0], 1e-12);
  const qx = quatFromAxisAngle([1, 0, 0], Math.PI / 2);
  // a ⊗ b applies b first: x → (z-rotation) y → (x-rotation) z
  closeV(quatRotate(quatMultiply(qx, qz), [1, 0, 0]), [0, 0, 1], 1e-12);
  const q = quatNormalize([0.3, -0.5, 0.2, 0.7]);
  const m = mat4FromQuat(q);
  const v = [0.4, -1.2, 2.5];
  const r = quatRotate(q, v);
  closeV(projectPoint(m, v[0], v[1], v[2]).slice(0, 3), r, 1e-12);
  close(Math.hypot(...r), Math.hypot(...v), 1e-12);
});

test('quatFromBasis sends the given world axes to the screen axes', () => {
  for (const name of ['iso', 'top', 'front', 'side']) {
    const q = viewOrientation(name);
    const ax = cameraAxes({ orientation: q });
    closeV(quatRotate(q, ax.right), [1, 0, 0], 1e-12);
    closeV(quatRotate(q, ax.up), [0, 1, 0], 1e-12);
    closeV(quatRotate(q, ax.back), [0, 0, 1], 1e-12);
  }
  // z is up: the world z axis points up on screen in the front, side and isometric views
  for (const name of ['front', 'side', 'iso']) assert.ok(quatRotate(viewOrientation(name), [0, 0, 1])[1] > 0.5, name);
  closeV(quatRotate(viewOrientation('top'), [0, 1, 0]), [0, 1, 0], 1e-12);
  closeV(quatRotate(viewOrientation('front'), [1, 0, 0]), [1, 0, 0], 1e-12);
  closeV(quatRotate(viewOrientation('side'), [0, 1, 0]), [1, 0, 0], 1e-12);
  assert.throws(() => viewOrientation('upside-down'), RangeError);
});

test('arcball: points map to the unit sphere, continuously across the rim', () => {
  closeV(arcballVector(0, 0), [0, 0, 1]);
  for (const [x, y] of [[0.2, 0.1], [0.7, 0.0], [0.5, 0.5], [1.5, -0.4], [3, 3]]) close(Math.hypot(...arcballVector(x, y)), 1, 1e-12);
  const r = Math.SQRT1_2;
  closeV(arcballVector(r - 1e-9, 0), arcballVector(r + 1e-9, 0), 1e-8);
});

test('arcball rotation takes the first vector to the second (Shoemake 1992) and the drag composes on the left', () => {
  const a = arcballVector(0.1, 0.2), b = arcballVector(-0.4, 0.3);
  const q = arcballRotation(a, b);
  closeV(quatRotate(q, a), b, 1e-12);
  closeV(arcballRotation(a, a), [0, 0, 0, 1], 1e-12);
  // dragging the front point of the ball to the right turns the scene about the screen y axis
  const q0 = viewOrientation('front');
  const q1 = orbitArcball(q0, arcballVector(0, 0), arcballVector(0.5, 0));
  const world = cameraAxes({ orientation: q0 }).back; // the world point facing the camera
  const after = quatRotate(q1, world);
  assert.ok(after[0] > 0.4 && Math.abs(after[1]) < 1e-12);
  close(Math.hypot(...q1), 1, 1e-12);
});

test('turntable orbit keeps the world z axis upright on screen', () => {
  let q = viewOrientation('iso');
  for (let i = 0; i < 12; i++) {
    q = orbitTurntable(q, 0.3, 0);
    const z = quatRotate(q, [0, 0, 1]);
    close(z[0], 0, 1e-12, 'z axis has no sideways component');
  }
});

test('camera matrices: target projects to the centre, eye sits at the stated distance', () => {
  const cam = { ...defaultCamera(), target: [10, -5, 3], distance: 40 };
  const m = cameraMatrices(cam, 1.6, 12);
  const p = projectPoint(m.viewProj, 10, -5, 3);
  close(p[0], 0, 1e-12); close(p[1], 0, 1e-12); close(p[3], 40, 1e-9);
  close(Math.hypot(m.eye[0] - 10, m.eye[1] + 5, m.eye[2] - 3), 40, 1e-9);
  assert.ok(m.near > 0 && m.near <= 28 && m.far >= 52);
  // the same point expressed relative to a model origin gives the same clip coordinates
  const o = [9, -4, 1];
  const m2 = cameraMatrices(cam, 1.6, 12, o);
  closeV(projectPoint(m2.viewProj, 12 - o[0], -7 - o[1], 2 - o[2]), projectPoint(m.viewProj, 12, -7, 2), 1e-9);
});

test('camera precision: a model 600 km from the origin projects exactly through its local origin', () => {
  const far = [612345.25, 6543210.5, -80];
  const cam = { ...defaultCamera(), target: far, distance: 5 };
  const m = cameraMatrices(cam, 1, 2, far);
  const ax = cameraAxes(cam);
  const p = projectPoint(m.viewProj, ax.right[0] * 0.001, ax.right[1] * 0.001, ax.right[2] * 0.001);
  // 1 mm to the right at 5 m distance
  close(p[0], 0.001 / (5 * Math.tan(cam.fovY / 2)), 1e-12);
});

test('fit: the bounding sphere touches the frustum in both wide and tall viewports', () => {
  const cam = defaultCamera();
  for (const aspect of [2, 1, 0.5]) {
    const d = fitDistance(3, cam.fovY, aspect);
    const half = Math.min(cam.fovY / 2, Math.atan(Math.tan(cam.fovY / 2) * aspect));
    close(d * Math.sin(half), 3, 1e-12);
  }
  const box = { min: [0, 0, 0], max: [4, 2, 1] };
  const c = fitCamera(cam, box, 1.5, 1);
  closeV(c.target, [2, 1, 0.5]);
  const m = cameraMatrices(c, 1.5, Math.hypot(4, 2, 1) / 2);
  for (const x of [0, 4]) for (const y of [0, 2]) for (const z of [0, 1]) {
    const p = projectPoint(m.viewProj, x, y, z);
    assert.ok(Math.abs(p[0]) <= 1 + 1e-9 && Math.abs(p[1]) <= 1 + 1e-9 && p[2] > -1 && p[2] < 1, `corner ${x},${y},${z} inside the view`);
  }
});

test('pan moves the target in the screen plane by the pixel size at the target', () => {
  const cam = { ...defaultCamera(), orientation: viewOrientation('top'), distance: 10 };
  const h = 500;
  const t = panTarget(cam, 50, 0, h);
  close(t[0], -50 * 2 * 10 * Math.tan(cam.fovY / 2) / h, 1e-12);
  close(t[1], 0, 1e-12); close(t[2], 0, 1e-12);
  const t2 = panTarget(cam, 0, 50, h);
  assert.ok(t2[1] > 0 && Math.abs(t2[0]) < 1e-12);
});

test('zoom keeps the world point under the cursor fixed', () => {
  const cam = { ...defaultCamera(), orientation: viewOrientation('top'), distance: 20, target: [1, 2, 0] };
  const aspect = 1.5, nx = 0.4, ny = -0.6;
  const h = 20 * Math.tan(cam.fovY / 2);
  const world = [1 + nx * aspect * h, 2 + ny * h, 0];
  const z = zoomCamera(cam, 0.5, nx, ny, aspect);
  close(z.distance, 10);
  const p = projectPoint(cameraMatrices(z, aspect, 5).viewProj, ...world);
  close(p[0], nx, 1e-12); close(p[1], ny, 1e-12);
});

test('screen ray passes through the eye and the unprojected point (perspective and orthographic)', () => {
  for (const ortho of [false, true]) {
    const cam = { ...defaultCamera(), ortho, target: [3, 1, 2], distance: 15 };
    const aspect = 1.3;
    const world = [4, 0.5, 2.5];
    const p = projectPoint(cameraMatrices(cam, aspect, 5).viewProj, ...world);
    const ray = screenRay(cam, aspect, p[0], p[1]);
    const v = [world[0] - ray.origin[0], world[1] - ray.origin[1], world[2] - ray.origin[2]];
    const t = v[0] * ray.dir[0] + v[1] * ray.dir[1] + v[2] * ray.dir[2];
    close(Math.hypot(v[0] - t * ray.dir[0], v[1] - t * ray.dir[1], v[2] - t * ray.dir[2]), 0, 1e-9, ortho ? 'ortho' : 'perspective');
    assert.ok(t > 0);
  }
});

test('colour map: end colours and monotonically rising luminance', () => {
  closeV(colourMap(0).map((v) => Math.round(v * 255)), [68, 1, 84]);
  closeV(colourMap(1).map((v) => Math.round(v * 255)), [253, 231, 37]);
  closeV(colourMap(-3), colourMap(0)); closeV(colourMap(7), colourMap(1)); closeV(colourMap(NaN), colourMap(0));
  let prev = -1;
  for (let i = 0; i <= 100; i++) {
    const c = colourMap(i / 100);
    const L = luminance(c);
    assert.ok(L > prev, `luminance rises at t = ${i / 100}`);
    prev = L;
    assert.ok(c.every((v) => v >= 0 && v <= 1));
  }
  // lightness alone separates the ends by a wide margin, so the ordering survives any colour-vision deficiency
  assert.ok(luminance(colourMap(1)) / luminance(colourMap(0)) > 20);
  const tex = colourMapTexture(256);
  assert.equal(tex.length, 1024);
  closeV(Array.from(tex.slice(0, 4)), [68, 1, 84, 255]);
  closeV(Array.from(tex.slice(1020)), [253, 231, 37, 255]);
  assert.equal(colourMapStops(9).length, 9);
  assert.equal(colourMapStops(9)[4], 'rgb(33,144,141)');
});

test('value range and normalisation handle NaN, clamping and a flat field', () => {
  const v = Float32Array.from([2, NaN, 6, 4, Infinity]);
  assert.deepEqual(valueRange(v), { min: 2, max: 6 });
  closeV(Array.from(normaliseValues(v, { min: 2, max: 6 })), [0, -1, 1, 0.5, -1]);
  closeV(Array.from(normaliseValues([0, 10], { min: 2, max: 6 })), [0, 1]);
  closeV(Array.from(normaliseValues([5, 5], { min: 5, max: 5 })), [0.5, 0.5]);
  assert.ok(Number.isNaN(valueRange([NaN]).min));
});

test('CSS colours and number labels', () => {
  closeV(parseColour('#0b6e8f'), [11 / 255, 110 / 255, 143 / 255]);
  closeV(parseColour(' #fff '), [1, 1, 1]);
  closeV(parseColour('rgb(255, 0, 51)'), [1, 0, 0.2]);
  closeV(parseColour('rgba(0 128 255 / 0.5)'), [0, 128 / 255, 1]);
  closeV(parseColour('', [0.1, 0.2, 0.3]), [0.1, 0.2, 0.3]);
  assert.equal(formatNumber(1234.5678), '1235');
  assert.equal(formatNumber(0.012345), '0.01235');
  assert.equal(formatNumber(2.5e7), '2.50e7');
  assert.equal(formatNumber(0), '0');
  assert.equal(formatNumber(NaN), '–');
  assert.equal(dimensionText({ min: [0, 0, -10], max: [120, 40, 15] }), '120 × 40 × 25 m');
});

test('unique edges of a closed cube: 18 for the triangulation (12 cube edges + 6 diagonals)', () => {
  const c = cube();
  const e = uniqueEdges(triEdgePairs(c.indices), 8);
  assert.equal(e.length / 2, 18);
  for (let i = 0; i < e.length; i += 2) assert.ok(e[i] < e[i + 1]);
  assert.equal(uniqueEdges(Uint32Array.from([3, 3, 1, 2, 2, 1]), 4).length, 2);
});

test('vertex normals of a cube point away from its centre and are unit length', () => {
  const c = cube();
  const n = vertexNormals(c.positions, c.indices);
  for (let v = 0; v < 8; v++) {
    close(Math.hypot(n[3 * v], n[3 * v + 1], n[3 * v + 2]), 1, 1e-6);
    const d = [c.positions[3 * v] - 0.5, c.positions[3 * v + 1] - 0.5, c.positions[3 * v + 2] - 0.5];
    assert.ok(n[3 * v] * d[0] + n[3 * v + 1] * d[1] + n[3 * v + 2] * d[2] > 0);
  }
});

test('boundary extraction: two hexahedra sharing a face give 10 quadrilaterals = 20 triangles and 20 outline edges', () => {
  const b = extractBoundary(twoHexes());
  assert.equal(b.tris.length / 3, 20);
  assert.equal(b.edges.length / 2, 20); // 12 + 12 − 4 shared, diagonals excluded
  assert.deepEqual([...new Set(b.triCell)].sort(), [0, 1]);
  // the interior face (nodes 1, 2, 5, 6) is not on the boundary
  const inner = new Set([1, 2, 5, 6]);
  for (let t = 0; t < b.tris.length; t += 3) assert.ok(!(inner.has(b.tris[t]) && inner.has(b.tris[t + 1]) && inner.has(b.tris[t + 2])));
  // closed boundary: every triangle edge is used by exactly two triangles
  const count = new Map();
  for (let t = 0; t < b.tris.length; t += 3) for (let k = 0; k < 3; k++) {
    const a = b.tris[t + k], c = b.tris[t + (k + 1) % 3];
    const key = a < c ? `${a},${c}` : `${c},${a}`;
    count.set(key, (count.get(key) || 0) + 1);
  }
  assert.ok([...count.values()].every((n) => n === 2));
});

test('boundary extraction: a single tetrahedron keeps its four faces', () => {
  const m = { kind: 'volume', positions: Float32Array.from([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1]), cells: [{ type: 'tet4', conn: Uint32Array.from([0, 1, 2, 3]) }] };
  const b = extractBoundary(m);
  assert.equal(b.tris.length / 3, 4);
  assert.equal(b.edges.length / 2, 6);
});

test('decimation keeps the requested number of evenly spread items', () => {
  assert.deepEqual(Array.from(decimateIndices(5, 10)), [0, 1, 2, 3, 4]);
  const d = decimateIndices(1000, 100);
  assert.equal(d.length, 100);
  assert.equal(d[0], 0); assert.equal(d[1], 10); assert.equal(d[99], 990);
});

test('grid mesh: node positions, vertical exaggeration, no-data cells skipped, decimation stride', () => {
  const g = { kind: 'grid', name: 'bathy', nx: 3, ny: 3, x0: 100, y0: 200, dx: 10, dy: 5, z: Float32Array.from([-10, -11, -12, -10, NaN, -12, -9, -10, -11]) };
  const m = gridMesh(g, { zScale: 2 });
  assert.equal(m.tris.length, 0); // every cell touches the missing centre node
  closeV(Array.from(m.coords.slice(0, 3)), [100, 200, -20]);
  closeV(Array.from(m.coords.slice(24, 27)), [120, 210, -22]);
  g.z[4] = -11;
  const m2 = gridMesh(g);
  assert.equal(m2.tris.length / 3, 8);
  // upward-facing triangles for positive spacings
  const [a, b, c] = m2.tris;
  const ux = m2.coords[3 * b] - m2.coords[3 * a], uy = m2.coords[3 * b + 1] - m2.coords[3 * a + 1], vx = m2.coords[3 * c] - m2.coords[3 * a], vy = m2.coords[3 * c + 1] - m2.coords[3 * a + 1];
  assert.ok(ux * vy - uy * vx > 0);
  const big = { nx: 1001, ny: 501, x0: 0, y0: 0, dx: 1, dy: 1, z: new Float32Array(1001 * 501) };
  const m3 = gridMesh(big, { maxNodes: 101 });
  assert.equal(m3.stride, 10); assert.equal(m3.nx, 101); assert.equal(m3.ny, 51);
  assert.equal(m3.source[m3.nx + 1], 10 * 1001 + 10);
  assert.throws(() => gridMesh({ nx: 2, ny: 2, z: new Float32Array(3) }), RangeError);
});

test('automatic exaggeration: flat seabed is stretched, steep terrain is not', () => {
  const flat = { nx: 101, ny: 101, dx: 100, dy: 100, x0: 0, y0: 0, z: Float32Array.from({ length: 10201 }, (_, i) => -50 - (i % 101) * 0.1) };
  assert.equal(autoExaggeration(flat), 50); // 10 km wide, 10 m relief
  const steep = { ...flat, z: Float32Array.from({ length: 10201 }, (_, i) => -(i % 101) * 30) };
  assert.equal(autoExaggeration(steep), 1);
  assert.equal(autoExaggeration({ ...flat, z: new Float32Array(10201) }), 1);
});

test('buildItem renders every model kind of the contract', () => {
  const s = buildItem(cube());
  assert.equal(s.tris.length, 36); assert.equal(s.nVerts, 8);
  closeV(s.origin, [0.5, 0.5, 0.5]);
  closeV(Array.from(s.positions.slice(0, 3)), [-0.5, -0.5, -0.5]);
  assert.equal(itemEdges(s).length / 2, 18);
  assert.deepEqual(itemBounds(s), { min: [0, 0, 0], max: [1, 1, 1] });

  const v = buildItem(twoHexes());
  assert.equal(v.tris.length / 3, 20); assert.equal(itemEdges(v).length / 2, 20); assert.equal(v.triElement.length, 20);
  assert.equal(v.counts.cells, 2);

  const f = buildItem({ kind: 'frame', name: 'jacket', nodes: Float64Array.from([0, 0, 0, 0, 0, 10, 5, 0, 10]), members: [{ a: 0, b: 1 }, { a: 1, b: 2 }] });
  assert.deepEqual(Array.from(f.wires), [0, 1, 1, 2]); assert.equal(f.tris, null);

  const p = buildItem({ kind: 'polyline', name: 'route', points: Float64Array.from([0, 0, 0, 10, 0, -1, 20, 5, -2]) });
  assert.deepEqual(Array.from(p.wires), [0, 1, 1, 2]);

  const n = 5000;
  const cloud = new Float32Array(3 * n);
  for (let i = 0; i < n; i++) { cloud[3 * i] = i; cloud[3 * i + 2] = -i / 100; }
  const pts = buildItem({ kind: 'points', name: 'soundings', positions: cloud }, { maxPoints: 500 });
  assert.equal(pts.nVerts, 500); assert.equal(pts.sourceCount, n); assert.ok(pts.isPoints);
  assert.equal(pts.source[1], 10);
  assert.equal(pts.valueLabel, 'Elevation');
  close(pts.trueBounds.max[0], n - 1); // bounds are those of the full cloud, not of the drawn subset
  assert.equal(buildItem({ kind: 'points', positions: cloud }, { colourBy: 'none' }).values, null);

  const g = buildItem({ kind: 'grid', name: 'bathy', nx: 4, ny: 3, x0: 0, y0: 0, dx: 1000, dy: 1000, z: Float32Array.from([-20, -21, -22, -23, -20, -21, -22, -23, -20, -21, -22, -23]) });
  assert.equal(g.tris.length / 3, 12); assert.equal(g.zScale, 50); assert.equal(g.valueUnit, 'm');
  assert.deepEqual(g.range, { min: -23, max: -20 });
  assert.deepEqual(g.trueBounds, { min: [0, 0, -23], max: [3000, 2000, -20] });
  close(itemBounds(g).min[2], -23 * 50, 1e-3);
  assert.equal(buildItem({ kind: 'grid', nx: 4, ny: 3, x0: 0, y0: 0, dx: 1000, dy: 1000, z: g.values }, { zScale: 3 }).zScale, 3);

  const sec = buildItem({ kind: 'section', name: 'tube', loops: [Float64Array.from([0, 0, 2, 0, 2, 1, 0, 1]), Float64Array.from([0.5, 0.25, 1.5, 0.25, 1, 0.75])] });
  assert.equal(sec.wires.length / 2, 7);
  assert.deepEqual(Array.from(sec.wires.slice(6, 8)), [3, 0]); // closing segment of the first loop
  assert.equal(sec.trueBounds.max[2], 0);

  // a route laid over an exaggerated grid is stretched by the same factor but reports its true size
  const over = buildItem({ kind: 'polyline', points: Float64Array.from([0, 0, -20, 10, 0, -22]) }, { zScale: 50 });
  close(itemBounds(over).min[2], -1100, 1e-3); assert.equal(over.trueBounds.min[2], -22);

  assert.throws(() => buildItem({ kind: 'hologram' }), RangeError);
  assert.throws(() => buildItem(null), RangeError);
});

test('buildItem: model3d plot with nested arrays, values and deformed shape', () => {
  const plot = {
    kind: 'model3d', title: 'Mode 1', nodes: [[0, 0, 0], [0, 0, 10], [0, 0, 20]], lines: [[0, 1], [1, 2]],
    nodeValues: [0, 0.4, 1], valueLabel: 'Utilisation', valueUnit: '-', deformed: [[0, 0, 0], [0.5, 0, 10], [2, 0, 20]],
  };
  const it = buildItem(plot);
  assert.equal(it.name, 'Mode 1');
  assert.deepEqual(Array.from(it.wires), [0, 1, 1, 2]);
  assert.deepEqual(it.range, { min: 0, max: 1 });
  assert.ok(it.deformed);
  close(it.positions[6] - it.base[6], 2, 1e-6);
  assert.deepEqual(it.trueBounds, { min: [0, 0, 0], max: [0, 0, 20] });
  // only nodes → drawn as points; indices out of range are dropped rather than drawn
  assert.ok(buildItem({ kind: 'model3d', nodes: [[0, 0, 0], [1, 1, 1]] }).isPoints);
  const bad = buildItem({ kind: 'model3d', nodes: [[0, 0, 0], [1, 0, 0], [0, 1, 0]], tris: [[0, 1, 2], [0, 1, 9]] });
  assert.equal(bad.tris.length, 3);
  // value array of the wrong length is ignored
  assert.equal(buildItem({ kind: 'model3d', nodes: [[0, 0, 0], [1, 1, 1]], lines: [[0, 1]], nodeValues: [1] }).values, null);
});

test('deformed shape: scale factor, automatic scale, restore and size check', () => {
  const it = buildItem({ kind: 'frame', nodes: Float64Array.from([0, 0, 0, 0, 0, 10]), members: [{ a: 0, b: 1 }] });
  const k = applyDeformed(it, [[0, 0, 0], [0.01, 0, 10]], 100);
  assert.equal(k, 100);
  close(it.positions[3], 1, 1e-5);
  const ka = applyDeformed(it, Float64Array.from([0, 0, 0, 0.01, 0, 10]), 'auto');
  close(ka, 0.05 * 10 / 0.01, 1e-3 * ka); // 5 % of the 10 m diagonal (single-precision buffers)
  applyDeformed(it, null);
  assert.equal(it.positions, it.base); assert.equal(it.deformed, false);
  assert.throws(() => applyDeformed(it, [[0, 0, 0]]), RangeError);
});

test('thick-line streams: six corners per segment, sides alternate across the line', () => {
  const P = Float32Array.from([0, 0, 0, 1, 0, 0, 1, 2, 0]);
  const w = thickLineBuffers(P, Uint32Array.from([0, 1, 1, 2]), Float32Array.from([0, 0.5, 1]));
  assert.equal(w.count, 12);
  closeV(Array.from(w.pos.slice(0, 3)), [0, 0, 0]); closeV(Array.from(w.other.slice(0, 3)), [1, 0, 0]);
  closeV(Array.from(w.pos.slice(6, 9)), [1, 0, 0]); closeV(Array.from(w.other.slice(6, 9)), [0, 0, 0]);
  // reconstruct the offset corners as the shader does (2-D, direction to the other end rotated by 90°)
  const corner = (v) => {
    const dx = w.other[3 * v] - w.pos[3 * v], dy = w.other[3 * v + 1] - w.pos[3 * v + 1], l = Math.hypot(dx, dy);
    return [w.pos[3 * v] - dy / l * w.side[v], w.pos[3 * v + 1] + dx / l * w.side[v]];
  };
  const quad = [0, 1, 2, 3, 4, 5].map(corner);
  const area = (a, b, c) => Math.abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2;
  // the two triangles tile a 1 × 2 rectangle around the first segment
  close(area(quad[0], quad[1], quad[2]) + area(quad[3], quad[4], quad[5]), 2, 1e-6);
  closeV(Array.from(w.t.slice(0, 6)), [0, 0, 0.5, 0, 0.5, 0.5]);
  closeV(Array.from(expandIndexed(P, 3, Uint32Array.from([2, 0]))), [1, 2, 0, 0, 0, 0]);
});

test('highlight geometry: edges, faces, volume cells, frame members and decimated vertices', () => {
  const s = buildItem(cube());
  const h = highlightGeometry(s, { edges: [[0, 1], [0, 99]], elements: [0, 5], vertices: [7] });
  assert.equal(h.segments.length, 6); assert.equal(h.triangles.length, 18); assert.equal(h.points.length, 3);
  const v = buildItem(twoHexes());
  assert.equal(highlightGeometry(v, { elements: [1] }).triangles.length / 9, 10); // five boundary quads of cell 1
  const f = buildItem({ kind: 'frame', nodes: Float64Array.from([0, 0, 0, 0, 0, 10, 5, 0, 10]), members: [{ a: 0, b: 1 }, { a: 1, b: 2 }] });
  const hf = highlightGeometry(f, { elements: [1] });
  closeV(Array.from(hf.segments), [-2.5, 0, 5, 2.5, 0, 5]);
  const cloud = new Float32Array(300);
  const pts = buildItem({ kind: 'points', positions: cloud }, { maxPoints: 10 });
  assert.equal(highlightGeometry(pts, { vertices: [10, 11] }).points.length, 3); // node 11 is not drawn
  assert.equal(highlightGeometry(s).segments.length, 0);
});

test('ray picking finds the nearest face of a cube and the hit point', () => {
  const c = cube();
  const hit = pickTriangles(c.positions, c.indices, [0.25, 0.3, 5], [0, 0, -1]);
  close(hit.t, 4, 1e-9); closeV(hit.point, [0.25, 0.3, 1], 1e-9);
  assert.ok(hit.tri === 2 || hit.tri === 3);
  assert.equal(pickTriangles(c.positions, c.indices, [3, 3, 5], [0, 0, -1]), null);
  // from inside the cube the nearest wall is hit, not the one behind the origin
  close(pickTriangles(c.positions, c.indices, [0.5, 0.5, 0.5], [1, 0, 0]).t, 0.5, 1e-9);
});

test('screen picking of segments and points honours the tolerance and returns the 3-D position', () => {
  const cam = { ...defaultCamera(), orientation: viewOrientation('top'), distance: 10 };
  const mvp = cameraMatrices(cam, 1, 5).viewProj;
  const P = Float32Array.from([-1, 0, 0, 1, 0, 0, 0, 2, 0]);
  const pairs = Uint32Array.from([0, 1]);
  const at = projectPoint(mvp, 0.5, 0, 0);
  const hit = pickScreen(P, pairs, mvp, at[0], at[1] + 0.01, 0.02, 0.02);
  assert.equal(hit.index, 0); assert.equal(hit.vertex, 1);
  closeV(hit.point, [0.5, 0, 0], 1e-6);
  assert.equal(pickScreen(P, pairs, mvp, at[0], at[1] + 0.05, 0.02, 0.02), null);
  const top = projectPoint(mvp, 0, 2, 0);
  assert.equal(pickScreen(P, null, mvp, top[0] + 0.005, top[1], 0.02, 0.02).vertex, 2);
  assert.equal(pickScreen(P, null, mvp, 0.9, 0.9, 0.02, 0.02), null);
});

test('scene bounds and description', () => {
  assert.equal(unionBounds([]), null);
  assert.deepEqual(unionBounds([{ min: [0, 0, 0], max: [1, 1, 1] }, null, { min: [-2, 0.5, 0], max: [0, 3, 0.5] }]), { min: [-2, 0, 0], max: [1, 3, 1] });
  assert.equal(describeScene([]), 'Empty 3-D view.');
  const text = describeScene([buildItem(cube()), buildItem({ kind: 'model3d', title: 'Tower', nodes: [[0, 0, 0], [0, 0, 80]], lines: [[0, 1]], nodeValues: [0.2, 0.9], valueLabel: 'Utilisation', valueUnit: '-' })]);
  assert.match(text, /^3-D view of 2 models\./);
  assert.match(text, /cube: surface with 8 vertices, 12 triangles/);
  assert.match(text, /Tower: model with 2 nodes, 1 lines; coloured by Utilisation \[-\] from 0\.2 to 0\.9/);
  assert.match(text, /Overall size 1 × 1 × 80 m/);
});

test('shader source is adapted to the context version', () => {
  const fs2 = shaderSource('void main() { gl_FragColor = vec4(1.0); }', 'fragment', true, true);
  assert.ok(fs2.startsWith('#version 300 es\n'));
  assert.match(fs2, /#define DERIV/); assert.match(fs2, /out vec4 fragColour/);
  const fs1 = shaderSource('void main() {}', 'fragment', false, true);
  assert.ok(fs1.startsWith('#extension GL_OES_standard_derivatives : enable'));
  assert.ok(!/DERIV/.test(shaderSource('void main() {}', 'fragment', false, false)));
  assert.ok(shaderSource('attribute vec3 aPos;', 'vertex', true).includes('#define attribute in'));
  assert.equal(shaderSource('attribute vec3 aPos;', 'vertex', false), 'attribute vec3 aPos;');
});
