/**
 * Dependency-free 3-D viewer (WebGL2, WebGL1 fallback, 2-D canvas wireframe when WebGL is unavailable).
 *
 * World axes: z is up. Renders every geometry model kind of CONTRACT §4 and the `model3d` plot type of §2.2.
 * The camera, arcball, colour-map, picking and buffer-building math are exported as pure functions that run
 * without a DOM (they are tested under Node); only `createViewer` touches the document.
 *
 * Conventions
 *  - Matrices are column-major Float64Array(16), as WebGL expects after conversion to Float32.
 *  - Quaternions are [x, y, z, w]; the camera orientation q rotates world vectors into camera space
 *    (camera looks down −z, +y is up on screen). View matrix V = T(0, 0, −distance) · R(q) · T(−target).
 *  - Vertex buffers are stored relative to a per-model origin (its bounding-box centre) so that models far
 *    from the global origin (UTM coordinates) keep their precision in single-precision buffers.
 *
 * References: Shoemake, "ARCBALL: a user interface for specifying three-dimensional orientation using a mouse",
 * Graphics Interface '92; Möller & Trumbore, "Fast, minimum storage ray/triangle intersection" (1997);
 * colour map after van der Walt & Smith, "A better default colormap for Matplotlib" (viridis, 2015).
 */
import { boundaryFaces } from '../mesh/topology.js';

/* ------------------------------------------------------------------------------------------------ matrices */

/** 4×4 identity (column-major). */
export function mat4Identity() {
  const m = new Float64Array(16);
  m[0] = m[5] = m[10] = m[15] = 1;
  return m;
}

/** Matrix product a · b (column-major). */
export function mat4Multiply(a, b) {
  const o = new Float64Array(16);
  for (let c = 0; c < 4; c++) {
    for (let r = 0; r < 4; r++) {
      o[4 * c + r] = a[r] * b[4 * c] + a[4 + r] * b[4 * c + 1] + a[8 + r] * b[4 * c + 2] + a[12 + r] * b[4 * c + 3];
    }
  }
  return o;
}

/**
 * Perspective projection (OpenGL clip space, z in [−1, 1]).
 * @param {number} fovY vertical field of view [rad]
 * @param {number} aspect width / height
 * @param {number} near near plane distance [m], > 0
 * @param {number} far far plane distance [m], > near
 */
export function mat4Perspective(fovY, aspect, near, far) {
  if (!(near > 0) || !(far > near)) throw new RangeError('near/far: need 0 < near < far');
  const f = 1 / Math.tan(fovY / 2);
  const m = new Float64Array(16);
  m[0] = f / aspect; m[5] = f; m[10] = (far + near) / (near - far); m[11] = -1; m[14] = 2 * far * near / (near - far);
  return m;
}

/** Symmetric orthographic projection with half-height `halfH` [m]. */
export function mat4Ortho(halfH, aspect, near, far) {
  if (!(halfH > 0) || !(far > near)) throw new RangeError('halfH/near/far: need halfH > 0 and near < far');
  const m = new Float64Array(16);
  m[0] = 1 / (halfH * aspect); m[5] = 1 / halfH; m[10] = -2 / (far - near); m[14] = -(far + near) / (far - near); m[15] = 1;
  return m;
}

/** Rotation matrix of a unit quaternion. */
export function mat4FromQuat(q) {
  const [x, y, z, w] = q;
  const m = mat4Identity();
  m[0] = 1 - 2 * (y * y + z * z); m[4] = 2 * (x * y - w * z); m[8] = 2 * (x * z + w * y);
  m[1] = 2 * (x * y + w * z); m[5] = 1 - 2 * (x * x + z * z); m[9] = 2 * (y * z - w * x);
  m[2] = 2 * (x * z - w * y); m[6] = 2 * (y * z + w * x); m[10] = 1 - 2 * (x * x + y * y);
  return m;
}

/** Transform a point by a 4×4 matrix; returns [x, y, z, w] after the perspective divide (w is the clip w). */
export function projectPoint(m, x, y, z) {
  const w = m[3] * x + m[7] * y + m[11] * z + m[15];
  const iw = w !== 0 ? 1 / w : 0;
  return [(m[0] * x + m[4] * y + m[8] * z + m[12]) * iw, (m[1] * x + m[5] * y + m[9] * z + m[13]) * iw, (m[2] * x + m[6] * y + m[10] * z + m[14]) * iw, w];
}

/* ---------------------------------------------------------------------------------------------- quaternions */

export function quatNormalize(q) {
  const l = Math.hypot(q[0], q[1], q[2], q[3]);
  return l > 0 ? [q[0] / l, q[1] / l, q[2] / l, q[3] / l] : [0, 0, 0, 1];
}

/** Hamilton product a ⊗ b: the rotation b followed by the rotation a. */
export function quatMultiply(a, b) {
  return [
    a[3] * b[0] + a[0] * b[3] + a[1] * b[2] - a[2] * b[1],
    a[3] * b[1] - a[0] * b[2] + a[1] * b[3] + a[2] * b[0],
    a[3] * b[2] + a[0] * b[1] - a[1] * b[0] + a[2] * b[3],
    a[3] * b[3] - a[0] * b[0] - a[1] * b[1] - a[2] * b[2],
  ];
}

export function quatFromAxisAngle(axis, angle) {
  const l = Math.hypot(axis[0], axis[1], axis[2]) || 1;
  const s = Math.sin(angle / 2) / l;
  return [axis[0] * s, axis[1] * s, axis[2] * s, Math.cos(angle / 2)];
}

/** Rotate vector v by unit quaternion q. */
export function quatRotate(q, v) {
  const [x, y, z, w] = q;
  const tx = 2 * (y * v[2] - z * v[1]), ty = 2 * (z * v[0] - x * v[2]), tz = 2 * (x * v[1] - y * v[0]);
  return [v[0] + w * tx + y * tz - z * ty, v[1] + w * ty + z * tx - x * tz, v[2] + w * tz + x * ty - y * tx];
}

/**
 * Camera orientation from the world-space directions of the screen axes: `right` (+x on screen), `up`
 * (+y on screen) and `back` (from the target towards the eye). The three must be orthonormal and right-handed.
 */
export function quatFromBasis(right, up, back) {
  const m00 = right[0], m01 = right[1], m02 = right[2];
  const m10 = up[0], m11 = up[1], m12 = up[2];
  const m20 = back[0], m21 = back[1], m22 = back[2];
  const tr = m00 + m11 + m22;
  let q;
  if (tr > 0) { const s = Math.sqrt(tr + 1) * 2; q = [(m21 - m12) / s, (m02 - m20) / s, (m10 - m01) / s, s / 4]; }
  else if (m00 > m11 && m00 > m22) { const s = Math.sqrt(1 + m00 - m11 - m22) * 2; q = [s / 4, (m01 + m10) / s, (m02 + m20) / s, (m21 - m12) / s]; }
  else if (m11 > m22) { const s = Math.sqrt(1 + m11 - m00 - m22) * 2; q = [(m01 + m10) / s, s / 4, (m12 + m21) / s, (m02 - m20) / s]; }
  else { const s = Math.sqrt(1 + m22 - m00 - m11) * 2; q = [(m02 + m20) / s, (m12 + m21) / s, s / 4, (m10 - m01) / s]; }
  return quatNormalize(q);
}

/* -------------------------------------------------------------------------------------------------- arcball */

/**
 * Map a pointer position to the arcball. Coordinates are normalised so that the shorter canvas side spans
 * [−1, 1], y up. Inside the ball the point lies on the unit sphere; outside it lies on the hyperbolic sheet
 * z = 1 / (2r) (Bell's blend), so dragging past the rim stays smooth.
 * @returns {number[]} unit vector in camera space
 */
export function arcballVector(nx, ny) {
  const r2 = nx * nx + ny * ny;
  const z = r2 <= 0.5 ? Math.sqrt(1 - r2) : 0.5 / Math.sqrt(r2);
  const l = Math.sqrt(r2 + z * z);
  return [nx / l, ny / l, z / l];
}

/** Shortest rotation taking unit vector `from` to unit vector `to`. */
export function arcballRotation(from, to) {
  const d = from[0] * to[0] + from[1] * to[1] + from[2] * to[2];
  if (d < -0.999999) return [0, 1, 0, 0];
  return quatNormalize([from[1] * to[2] - from[2] * to[1], from[2] * to[0] - from[0] * to[2], from[0] * to[1] - from[1] * to[0], 1 + d]);
}

/* --------------------------------------------------------------------------------------------------- camera */

const S2 = Math.SQRT1_2, S3 = 1 / Math.sqrt(3), S6 = 1 / Math.sqrt(6);

/**
 * Orientation of a named view (z up).
 *  iso   — from (+x, −y, +z);  top — looking down, x right and y up on screen;
 *  front — looking along +y, x right;  side — looking along −x, y right.
 */
export function viewOrientation(name) {
  switch (name) {
    case 'top': return quatFromBasis([1, 0, 0], [0, 1, 0], [0, 0, 1]);
    case 'front': return quatFromBasis([1, 0, 0], [0, 0, 1], [0, -1, 0]);
    case 'side': return quatFromBasis([0, 1, 0], [0, 0, 1], [1, 0, 0]);
    case 'iso': return quatFromBasis([S2, S2, 0], [-S6, S6, 2 * S6], [S3, -S3, S3]);
    default: throw new RangeError(`view: unknown view '${name}'`);
  }
}

/** Default camera state. */
export function defaultCamera() {
  return { target: [0, 0, 0], distance: 10, orientation: viewOrientation('iso'), fovY: 35 * Math.PI / 180, ortho: false };
}

/** World-space screen axes of a camera: { right, up, back } (unit vectors). */
export function cameraAxes(cam) {
  const c = [-cam.orientation[0], -cam.orientation[1], -cam.orientation[2], cam.orientation[3]];
  return { right: quatRotate(c, [1, 0, 0]), up: quatRotate(c, [0, 1, 0]), back: quatRotate(c, [0, 0, 1]) };
}

/**
 * View and projection matrices.
 * @param {object} cam { target [m], distance [m], orientation, fovY [rad], ortho }
 * @param {number} aspect viewport width / height
 * @param {number} radius radius of the sphere around `sceneCentre` that must stay inside the depth range [m]
 * @param {number[]} [origin] model origin subtracted from the buffers [m] (default [0, 0, 0])
 * @param {number[]} [sceneCentre] centre of the scene bounds (default: the target)
 * @returns {{ view: Float64Array, proj: Float64Array, viewProj: Float64Array, eye: number[], near: number, far: number }}
 */
export function cameraMatrices(cam, aspect, radius, origin = [0, 0, 0], sceneCentre = cam.target) {
  const rot = mat4FromQuat(cam.orientation);
  const o = [origin[0] - cam.target[0], origin[1] - cam.target[1], origin[2] - cam.target[2]];
  const view = rot;
  view[12] = rot[0] * o[0] + rot[4] * o[1] + rot[8] * o[2];
  view[13] = rot[1] * o[0] + rot[5] * o[1] + rot[9] * o[2];
  view[14] = rot[2] * o[0] + rot[6] * o[1] + rot[10] * o[2] - cam.distance;
  const back = [rot[2], rot[6], rot[10]];
  // depth of the scene centre in front of the eye
  const dc = cam.distance - ((sceneCentre[0] - cam.target[0]) * back[0] + (sceneCentre[1] - cam.target[1]) * back[1] + (sceneCentre[2] - cam.target[2]) * back[2]);
  const r = Math.max(radius, 1e-9);
  const far = Math.max(dc + r * 1.05, cam.distance * 1e-3, 1e-9);
  let proj, near;
  if (cam.ortho) {
    near = Math.min(dc - r * 1.05, far - 2.1 * r);
    proj = mat4Ortho(cam.distance * Math.tan(cam.fovY / 2), aspect, near, far);
  } else {
    near = Math.max(dc - r * 1.05, far * 2e-4);
    proj = mat4Perspective(cam.fovY, aspect, near, far);
  }
  const eye = [cam.target[0] + back[0] * cam.distance, cam.target[1] + back[1] * cam.distance, cam.target[2] + back[2] * cam.distance];
  return { view, proj, viewProj: mat4Multiply(proj, view), eye, near, far };
}

/** Camera distance at which a sphere of the given radius fills the view (both projections). */
export function fitDistance(radius, fovY, aspect) {
  const half = fovY / 2;
  const halfX = Math.atan(Math.tan(half) * aspect);
  return Math.max(radius, 1e-9) / Math.sin(Math.min(half, halfX));
}

/** Camera framing an axis-aligned box { min, max }, keeping the current orientation. */
export function fitCamera(cam, box, aspect, margin = 1.12) {
  const c = [(box.min[0] + box.max[0]) / 2, (box.min[1] + box.max[1]) / 2, (box.min[2] + box.max[2]) / 2];
  const r = Math.hypot(box.max[0] - box.min[0], box.max[1] - box.min[1], box.max[2] - box.min[2]) / 2 || 1;
  return { ...cam, target: c, distance: fitDistance(r, cam.fovY, aspect) * margin };
}

/** Orbit by arcball drag between two arcball vectors, starting from orientation q0. */
export function orbitArcball(q0, from, to) {
  return quatNormalize(quatMultiply(arcballRotation(from, to), q0));
}

/** Turntable orbit: `dAz` about the world z axis and `dEl` about the screen x axis [rad]. */
export function orbitTurntable(q0, dAz, dEl) {
  return quatNormalize(quatMultiply(quatFromAxisAngle([1, 0, 0], dEl), quatMultiply(q0, quatFromAxisAngle([0, 0, 1], dAz))));
}

/** Pan by a screen displacement [px] (x right, y down) in a viewport `height` px tall. Returns the new target. */
export function panTarget(cam, dxPx, dyPx, height) {
  const { right, up } = cameraAxes(cam);
  const s = 2 * cam.distance * Math.tan(cam.fovY / 2) / Math.max(1, height);
  return [cam.target[0] + (-dxPx * right[0] + dyPx * up[0]) * s, cam.target[1] + (-dxPx * right[1] + dyPx * up[1]) * s, cam.target[2] + (-dxPx * right[2] + dyPx * up[2]) * s];
}

/**
 * Zoom by `factor` (< 1 moves in) keeping the point under the cursor fixed.
 * @param {number} nx cursor x in normalised device coordinates [−1, 1]
 * @param {number} ny cursor y in normalised device coordinates [−1, 1], up positive
 */
export function zoomCamera(cam, factor, nx = 0, ny = 0, aspect = 1) {
  const d = Math.min(Math.max(cam.distance * factor, 1e-9), 1e12);
  const { right, up } = cameraAxes(cam);
  const h = Math.tan(cam.fovY / 2) * (cam.distance - d);
  const ox = nx * aspect * h, oy = ny * h;
  return { ...cam, distance: d, target: [cam.target[0] + right[0] * ox + up[0] * oy, cam.target[1] + right[1] * ox + up[1] * oy, cam.target[2] + right[2] * ox + up[2] * oy] };
}

/** Pick ray through a point of the viewport given in normalised device coordinates. */
export function screenRay(cam, aspect, nx, ny) {
  const { right, up, back } = cameraAxes(cam);
  const t = Math.tan(cam.fovY / 2);
  const eye = [cam.target[0] + back[0] * cam.distance, cam.target[1] + back[1] * cam.distance, cam.target[2] + back[2] * cam.distance];
  if (cam.ortho) {
    const h = cam.distance * t;
    return { origin: [eye[0] + right[0] * nx * aspect * h + up[0] * ny * h, eye[1] + right[1] * nx * aspect * h + up[1] * ny * h, eye[2] + right[2] * nx * aspect * h + up[2] * ny * h], dir: [-back[0], -back[1], -back[2]] };
  }
  const d = [right[0] * nx * aspect * t + up[0] * ny * t - back[0], right[1] * nx * aspect * t + up[1] * ny * t - back[1], right[2] * nx * aspect * t + up[2] * ny * t - back[2]];
  const l = Math.hypot(d[0], d[1], d[2]);
  return { origin: eye, dir: [d[0] / l, d[1] / l, d[2] / l] };
}

/* --------------------------------------------------------------------------------------------------- colour */

// Anchor colours of the viridis map at equal steps of 1/8 (sRGB, 0–255). Lightness rises monotonically and the
// hue path avoids the red–green axis, so the map reads correctly with any colour-vision deficiency and in print.
const MAP_ANCHORS = [[68, 1, 84], [71, 44, 122], [59, 81, 139], [44, 113, 142], [33, 144, 141], [39, 173, 129], [92, 200, 99], [170, 220, 50], [253, 231, 37]];

/**
 * Sequential colour map, perceptually uniform and colour-blind safe.
 * @param {number} t position in [0, 1] (clamped)
 * @returns {number[]} [r, g, b] in 0–1
 */
export function colourMap(t) {
  const u = Math.min(1, Math.max(0, Number.isFinite(t) ? t : 0)) * (MAP_ANCHORS.length - 1);
  const i = Math.min(MAP_ANCHORS.length - 2, Math.floor(u));
  const f = u - i;
  const a = MAP_ANCHORS[i], b = MAP_ANCHORS[i + 1];
  return [(a[0] + f * (b[0] - a[0])) / 255, (a[1] + f * (b[1] - a[1])) / 255, (a[2] + f * (b[2] - a[2])) / 255];
}

/** RGBA lookup table of the colour map (n texels, 4 bytes each) for use as a texture. */
export function colourMapTexture(n = 256) {
  const out = new Uint8Array(4 * n);
  for (let i = 0; i < n; i++) {
    const c = colourMap(i / (n - 1));
    out[4 * i] = Math.round(c[0] * 255); out[4 * i + 1] = Math.round(c[1] * 255); out[4 * i + 2] = Math.round(c[2] * 255); out[4 * i + 3] = 255;
  }
  return out;
}

/** CSS colour stops of the map, low to high, for a legend gradient. */
export function colourMapStops(n = 9) {
  const out = [];
  for (let i = 0; i < n; i++) { const c = colourMap(i / (n - 1)); out.push(`rgb(${Math.round(c[0] * 255)},${Math.round(c[1] * 255)},${Math.round(c[2] * 255)})`); }
  return out;
}

/** Relative luminance (WCAG) of an sRGB colour given as [r, g, b] in 0–1. */
export function luminance(c) {
  const f = (v) => (v <= 0.03928 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4);
  return 0.2126 * f(c[0]) + 0.7152 * f(c[1]) + 0.0722 * f(c[2]);
}

/** Finite minimum and maximum of a value array: { min, max } (both NaN when nothing is finite). */
export function valueRange(values) {
  let lo = Infinity, hi = -Infinity;
  for (let i = 0; i < values.length; i++) { const v = values[i]; if (v < lo) lo = v; if (v > hi) hi = v; }
  if (!(lo <= hi) || !Number.isFinite(lo) || !Number.isFinite(hi)) {
    lo = Infinity; hi = -Infinity;
    for (let i = 0; i < values.length; i++) { const v = values[i]; if (Number.isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; } }
  }
  return lo <= hi ? { min: lo, max: hi } : { min: NaN, max: NaN };
}

/**
 * Map values to colour-map coordinates: (v − min) / (max − min) clamped to [0, 1]; non-finite values give −1
 * (drawn in the no-data colour). A zero-width range maps everything to 0.5.
 */
export function normaliseValues(values, range, out = new Float32Array(values.length)) {
  const lo = range.min, span = range.max - range.min;
  for (let i = 0; i < values.length; i++) {
    const v = values[i];
    out[i] = !Number.isFinite(v) ? -1 : span > 0 ? Math.min(1, Math.max(0, (v - lo) / span)) : 0.5;
  }
  return out;
}

/** Parse a CSS colour (#rgb, #rrggbb, rgb(), rgba()) to [r, g, b] in 0–1; returns `fallback` if not understood. */
export function parseColour(str, fallback = [0.5, 0.5, 0.5]) {
  if (Array.isArray(str)) return str;
  const s = String(str || '').trim();
  let m = /^#([0-9a-f]{3})$/i.exec(s);
  if (m) return [0, 1, 2].map((i) => parseInt(m[1][i] + m[1][i], 16) / 255);
  m = /^#([0-9a-f]{6})([0-9a-f]{2})?$/i.exec(s);
  if (m) return [0, 2, 4].map((i) => parseInt(m[1].slice(i, i + 2), 16) / 255);
  m = /^rgba?\(\s*([\d.]+)[\s,]+([\d.]+)[\s,]+([\d.]+)/i.exec(s);
  if (m) return [Number(m[1]) / 255, Number(m[2]) / 255, Number(m[3]) / 255];
  return fallback;
}

const mix = (a, b, t) => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
const css = (c, a = 1) => `rgba(${Math.round(c[0] * 255)},${Math.round(c[1] * 255)},${Math.round(c[2] * 255)},${a})`;

/** Compact number for labels: four significant figures, exponent form outside 1e-3 … 1e6. */
export function formatNumber(v, sig = 4) {
  if (!Number.isFinite(v)) return '–';
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1e6 || a < 1e-3) return v.toExponential(Math.max(0, sig - 2)).replace('e+', 'e');
  return String(Number(v.toPrecision(sig)));
}

/* ------------------------------------------------------------------------------------------ buffer building */

/** Area-weighted vertex normals of an indexed triangle mesh (Float32Array, 3 per vertex). */
export function vertexNormals(positions, tris) {
  const n = new Float32Array(positions.length);
  for (let f = 0; f < tris.length; f += 3) {
    const a = 3 * tris[f], b = 3 * tris[f + 1], c = 3 * tris[f + 2];
    const ux = positions[b] - positions[a], uy = positions[b + 1] - positions[a + 1], uz = positions[b + 2] - positions[a + 2];
    const vx = positions[c] - positions[a], vy = positions[c + 1] - positions[a + 1], vz = positions[c + 2] - positions[a + 2];
    const x = uy * vz - uz * vy, y = uz * vx - ux * vz, z = ux * vy - uy * vx;
    n[a] += x; n[a + 1] += y; n[a + 2] += z; n[b] += x; n[b + 1] += y; n[b + 2] += z; n[c] += x; n[c + 1] += y; n[c + 2] += z;
  }
  for (let i = 0; i < n.length; i += 3) {
    const l = Math.hypot(n[i], n[i + 1], n[i + 2]);
    if (l > 0) { n[i] /= l; n[i + 1] /= l; n[i + 2] /= l; } else n[i + 2] = 1;
  }
  return n;
}

/** Vertex pairs of the three edges of every triangle (with repeats). */
export function triEdgePairs(tris) {
  const out = new Uint32Array(tris.length * 2);
  for (let f = 0, o = 0; f < tris.length; f += 3) {
    out[o++] = tris[f]; out[o++] = tris[f + 1]; out[o++] = tris[f + 1]; out[o++] = tris[f + 2]; out[o++] = tris[f + 2]; out[o++] = tris[f];
  }
  return out;
}

/** Unique undirected edges of a list of vertex pairs (collapsed pairs dropped), sorted by lower vertex id. */
export function uniqueEdges(pairs, nVerts) {
  const m = pairs.length >> 1;
  const keys = new Float64Array(m);
  let k = 0;
  for (let i = 0; i < m; i++) {
    const a = pairs[2 * i], b = pairs[2 * i + 1];
    if (a !== b) keys[k++] = a < b ? a * nVerts + b : b * nVerts + a;
  }
  const s = keys.subarray(0, k).sort();
  const out = new Uint32Array(2 * k);
  let c = 0, prev = -1;
  for (let i = 0; i < k; i++) {
    if (s[i] === prev) continue;
    prev = s[i];
    const lo = Math.floor(prev / nVerts);
    out[2 * c] = lo; out[2 * c + 1] = prev - lo * nVerts; c++;
  }
  return out.slice(0, 2 * c);
}

/**
 * Boundary of a volume mesh ready for drawing: cell faces owned by a single cell, quadrilaterals split into two
 * triangles, plus the face outline edges (the quadrilateral diagonals are not edges).
 * @param {object} model CONTRACT §4 volume model
 * @returns {{ tris: Uint32Array, triCell: Uint32Array, edges: Uint32Array }} `triCell[t]` is the running cell
 *   index owning triangle t.
 */
export function extractBoundary(model) {
  const bf = boundaryFaces(model);
  const nt = bf.tris.length / 3, nq = bf.quads.length / 4;
  const tris = new Uint32Array(3 * (nt + 2 * nq));
  const triCell = new Uint32Array(nt + 2 * nq);
  const pairs = new Uint32Array(6 * nt + 8 * nq);
  tris.set(bf.tris);
  triCell.set(bf.triCell);
  pairs.set(triEdgePairs(bf.tris));
  for (let q = 0; q < nq; q++) {
    const a = bf.quads[4 * q], b = bf.quads[4 * q + 1], c = bf.quads[4 * q + 2], d = bf.quads[4 * q + 3];
    tris.set([a, b, c, a, c, d], 3 * (nt + 2 * q));
    triCell[nt + 2 * q] = triCell[nt + 2 * q + 1] = bf.quadCell[q];
    pairs.set([a, b, b, c, c, d, d, a], 6 * nt + 8 * q);
  }
  return { tris, triCell, edges: uniqueEdges(pairs, (model.positions.length / 3) | 0) };
}

/** Evenly spaced subset of `count` items, at most `max` of them: returns the kept indices. */
export function decimateIndices(count, max) {
  if (!(max > 0) || count <= max) { const all = new Uint32Array(count); for (let i = 0; i < count; i++) all[i] = i; return all; }
  const out = new Uint32Array(max);
  const step = count / max;
  for (let i = 0; i < max; i++) out[i] = Math.min(count - 1, Math.floor(i * step));
  return out;
}

/**
 * Vertical exaggeration that makes the relief of a grid readable: 1 when the relief already exceeds about a
 * seventh of the plan extent, otherwise the largest of 2, 5, 10, 20, 50 that keeps it below a fifth.
 */
export function autoExaggeration(grid) {
  const r = valueRange(grid.z);
  const relief = r.max - r.min;
  const ext = Math.max(Math.abs((grid.nx - 1) * grid.dx), Math.abs((grid.ny - 1) * grid.dy));
  if (!(relief > 0) || !(ext > 0)) return 1;
  const e = 0.2 * ext / relief;
  let best = 1;
  for (const s of [2, 5, 10, 20, 50]) if (s <= e) best = s;
  return best;
}

/**
 * Height-field mesh of a grid, decimated to at most `maxNodes` nodes along each side. Cells touching a no-data
 * node are left out.
 * @returns {{ coords: Float64Array, tris: Uint32Array, values: Float32Array, source: Uint32Array, nx: number, ny: number, stride: number }}
 *   absolute coordinates with z multiplied by `zScale`, true elevations in `values`, and for every mesh vertex
 *   the index of its grid node in `source`.
 */
export function gridMesh(grid, { maxNodes = 400, zScale = 1 } = {}) {
  const { nx, ny } = grid;
  if (!(nx > 0) || !(ny > 0) || !grid.z || grid.z.length < nx * ny) throw new RangeError('grid: nx·ny must match z.length');
  const s = Math.max(1, Math.ceil((Math.max(nx, ny) - 1) / Math.max(1, maxNodes - 1)));
  const mx = Math.floor((nx - 1) / s) + 1, my = Math.floor((ny - 1) / s) + 1;
  const coords = new Float64Array(3 * mx * my), values = new Float32Array(mx * my), source = new Uint32Array(mx * my);
  for (let j = 0; j < my; j++) for (let i = 0; i < mx; i++) {
    const g = j * s * nx + i * s, v = j * mx + i;
    const z = grid.z[g];
    coords[3 * v] = grid.x0 + i * s * grid.dx; coords[3 * v + 1] = grid.y0 + j * s * grid.dy; coords[3 * v + 2] = Number.isFinite(z) ? z * zScale : 0;
    values[v] = z; source[v] = g;
  }
  const tris = new Uint32Array(6 * Math.max(0, mx - 1) * Math.max(0, my - 1));
  const ccw = grid.dx * grid.dy > 0;
  let o = 0;
  for (let j = 0; j + 1 < my; j++) for (let i = 0; i + 1 < mx; i++) {
    const a = j * mx + i, b = a + 1, c = a + mx + 1, d = a + mx;
    if (!Number.isFinite(values[a] + values[b] + values[c] + values[d])) continue;
    if (ccw) { tris[o++] = a; tris[o++] = b; tris[o++] = c; tris[o++] = a; tris[o++] = c; tris[o++] = d; } else { tris[o++] = a; tris[o++] = c; tris[o++] = b; tris[o++] = a; tris[o++] = d; tris[o++] = c; }
  }
  return { coords, tris: tris.slice(0, o), values, source, nx: mx, ny: my, stride: s };
}

function flat3(a, name) {
  if (!a) return new Float64Array(0);
  if (ArrayBuffer.isView(a)) return a;
  if (a.length && typeof a[0] === 'number') return Float64Array.from(a);
  const out = new Float64Array(3 * a.length);
  for (let i = 0; i < a.length; i++) {
    const p = a[i];
    if (!p || p.length < 2) throw new RangeError(`${name}[${i}]: expected [x, y, z]`);
    out[3 * i] = p[0]; out[3 * i + 1] = p[1]; out[3 * i + 2] = p[2] || 0;
  }
  return out;
}

function flatIdx(a, per) {
  if (!a) return null;
  if (ArrayBuffer.isView(a)) return a instanceof Uint32Array ? a : Uint32Array.from(a);
  if (a.length && typeof a[0] === 'number') return Uint32Array.from(a);
  const out = new Uint32Array(per * a.length);
  for (let i = 0; i < a.length; i++) for (let k = 0; k < per; k++) out[per * i + k] = a[i][k];
  return out;
}

function boundsOf(coords) {
  const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i + 2 < coords.length; i += 3) {
    const x = coords[i], y = coords[i + 1], z = coords[i + 2];
    if (!Number.isFinite(x + y + z)) continue;
    if (x < min[0]) min[0] = x; if (x > max[0]) max[0] = x;
    if (y < min[1]) min[1] = y; if (y > max[1]) max[1] = y;
    if (z < min[2]) min[2] = z; if (z > max[2]) max[2] = z;
  }
  return min[0] <= max[0] ? { min, max } : null;
}

function dropBadIndices(idx, per, nVerts) {
  if (!idx) return null;
  let ok = true;
  for (let i = 0; i < idx.length; i++) if (!(idx[i] < nVerts)) { ok = false; break; }
  if (ok && idx.length % per === 0) return idx;
  const out = [];
  for (let i = 0; i + per <= idx.length; i += per) {
    let good = true;
    for (let k = 0; k < per; k++) if (!(idx[i + k] < nVerts)) good = false;
    if (good) for (let k = 0; k < per; k++) out.push(idx[i + k]);
  }
  return Uint32Array.from(out);
}

/**
 * Normalise any model into the arrays the renderers draw.
 * @param {object} model CONTRACT §4 model or a `model3d` plot object
 * @param {object} [style] { zScale (vertical exaggeration: automatic for grids, 1 for other kinds), maxPoints = 300000,
 *   maxGridNodes = 400, colourBy: 'z' | 'none' }
 * @returns {object} render item:
 *   kind, name, origin [m], nVerts, positions (Float32Array, relative to origin, current shape), base (same,
 *   undeformed), tris / triElement, wires / wireElement (line segments as vertex pairs: frame members, polylines,
 *   outlines), isPoints, presetEdges (volume face outlines), source (render vertex → model node, when decimated),
 *   sourceCount, values (Float32Array per render vertex or null), valueLabel, valueUnit, range, zScale,
 *   trueBounds (bounds of the model as given, without exaggeration or deformation), counts.
 */
export function buildItem(model, style = {}) {
  if (!model || typeof model !== 'object') throw new RangeError('model: a geometry model is required');
  const it = {
    kind: model.kind, name: model.name || model.title || model.kind, tris: null, triElement: null, wires: null, wireElement: null,
    isPoints: false, presetEdges: null, source: null, values: null, valueLabel: '', valueUnit: '', range: null, zScale: 1, counts: {},
  };
  let coords, trueBounds = null, deformed = null;
  const colourZ = style.colourBy !== 'none';
  switch (model.kind) {
    case 'surface':
      coords = model.positions; it.tris = flatIdx(model.indices, 3);
      it.counts = { vertices: coords.length / 3, triangles: it.tris.length / 3 };
      break;
    case 'volume': {
      coords = model.positions;
      const b = extractBoundary(model);
      it.tris = b.tris; it.triElement = b.triCell; it.presetEdges = b.edges;
      it.counts = { nodes: coords.length / 3, cells: (model.cells || []).reduce((s, c) => s + c.conn.length / ({ tet4: 4, tet10: 10, hex8: 8, hex20: 20, wedge6: 6, pyr5: 5 }[c.type] || 1), 0), 'boundary triangles': b.tris.length / 3 };
      break;
    }
    case 'frame': {
      coords = model.nodes;
      const m = model.members || [];
      it.wires = new Uint32Array(2 * m.length);
      for (let i = 0; i < m.length; i++) { it.wires[2 * i] = m[i].a; it.wires[2 * i + 1] = m[i].b; }
      it.counts = { nodes: coords.length / 3, members: m.length };
      break;
    }
    case 'polyline': {
      coords = model.points;
      const n = (coords.length / 3) | 0;
      it.wires = new Uint32Array(2 * Math.max(0, n - 1));
      for (let i = 0; i + 1 < n; i++) { it.wires[2 * i] = i; it.wires[2 * i + 1] = i + 1; }
      it.counts = { points: n };
      break;
    }
    case 'points': {
      const n = (model.positions.length / 3) | 0;
      const keep = decimateIndices(n, style.maxPoints ?? 300000);
      if (keep.length < n) {
        coords = new Float64Array(3 * keep.length);
        for (let i = 0; i < keep.length; i++) { coords[3 * i] = model.positions[3 * keep[i]]; coords[3 * i + 1] = model.positions[3 * keep[i] + 1]; coords[3 * i + 2] = model.positions[3 * keep[i] + 2]; }
        it.source = keep; trueBounds = boundsOf(model.positions);
      } else coords = model.positions;
      it.isPoints = true; it.sourceCount = n;
      it.counts = { points: n, ...(keep.length < n ? { 'points drawn': keep.length } : {}) };
      if (colourZ) {
        it.values = new Float32Array(keep.length);
        for (let i = 0; i < keep.length; i++) it.values[i] = coords[3 * i + 2];
        it.valueLabel = 'Elevation'; it.valueUnit = 'm';
      }
      break;
    }
    case 'grid': {
      it.zScale = style.zScale > 0 ? style.zScale : autoExaggeration(model);
      const g = gridMesh(model, { maxNodes: style.maxGridNodes ?? 400, zScale: it.zScale });
      coords = g.coords; it.tris = g.tris; it.source = g.source; it.sourceCount = model.nx * model.ny;
      if (colourZ) { it.values = g.values; it.valueLabel = 'Elevation'; it.valueUnit = 'm'; }
      const r = valueRange(model.z);
      const x1 = model.x0 + (model.nx - 1) * model.dx, y1 = model.y0 + (model.ny - 1) * model.dy;
      trueBounds = Number.isFinite(r.min) ? { min: [Math.min(model.x0, x1), Math.min(model.y0, y1), r.min], max: [Math.max(model.x0, x1), Math.max(model.y0, y1), r.max] } : null;
      it.counts = { nodes: model.nx * model.ny, columns: model.nx, rows: model.ny, ...(g.stride > 1 ? { 'drawn every': g.stride } : {}) };
      break;
    }
    case 'section': {
      const loops = model.loops || [];
      let n = 0;
      for (const L of loops) n += (L.length / 2) | 0;
      coords = new Float64Array(3 * n);
      const w = [];
      let o = 0;
      for (const L of loops) {
        const m = (L.length / 2) | 0;
        for (let i = 0; i < m; i++) { coords[3 * (o + i)] = L[2 * i]; coords[3 * (o + i) + 1] = L[2 * i + 1]; if (m > 1) w.push(o + i, o + (i + 1) % m); }
        o += m;
      }
      it.wires = Uint32Array.from(w);
      it.counts = { loops: loops.length, points: n };
      break;
    }
    case 'model3d': {
      coords = flat3(model.nodes, 'nodes');
      it.tris = flatIdx(model.tris, 3);
      it.wires = flatIdx(model.lines, 2);
      if (it.tris && !it.tris.length) it.tris = null;
      if (it.wires && !it.wires.length) it.wires = null;
      if (!it.tris && !it.wires) it.isPoints = true;
      if (model.nodeValues && model.nodeValues.length === coords.length / 3) {
        it.values = Float32Array.from(model.nodeValues);
        it.valueLabel = model.valueLabel || 'Value'; it.valueUnit = model.valueUnit || '';
      }
      if (model.deformed && model.deformed.length) deformed = flat3(model.deformed, 'deformed');
      it.counts = { nodes: coords.length / 3, ...(it.wires ? { lines: it.wires.length / 2 } : {}), ...(it.tris ? { triangles: it.tris.length / 3 } : {}) };
      break;
    }
    default: throw new RangeError(`model.kind: '${model.kind}' cannot be displayed`);
  }
  if (model.kind !== 'grid' && style.zScale > 0 && style.zScale !== 1) {
    // overlay on an exaggerated grid: stretch this model by the same factor, report its true size
    trueBounds = trueBounds || boundsOf(coords);
    const scaled = new Float64Array(coords.length);
    for (let i = 0; i < coords.length; i += 3) { scaled[i] = coords[i]; scaled[i + 1] = coords[i + 1]; scaled[i + 2] = coords[i + 2] * style.zScale; }
    coords = scaled; it.zScale = style.zScale;
  }
  const nVerts = (coords.length / 3) | 0;
  it.nVerts = nVerts;
  if (it.sourceCount === undefined) it.sourceCount = nVerts;
  it.tris = dropBadIndices(it.tris, 3, nVerts);
  it.wires = dropBadIndices(it.wires, 2, nVerts);
  const bb = boundsOf(coords) || { min: [0, 0, 0], max: [0, 0, 0] };
  it.origin = [(bb.min[0] + bb.max[0]) / 2, (bb.min[1] + bb.max[1]) / 2, (bb.min[2] + bb.max[2]) / 2];
  it.trueBounds = trueBounds || bb;
  const base = new Float32Array(3 * nVerts);
  for (let i = 0; i < nVerts; i++) {
    base[3 * i] = coords[3 * i] - it.origin[0]; base[3 * i + 1] = coords[3 * i + 1] - it.origin[1]; base[3 * i + 2] = coords[3 * i + 2] - it.origin[2];
  }
  it.base = base;
  it.positions = base;
  it.deformed = false;
  if (it.values) it.range = valueRange(it.values);
  if (deformed) applyDeformed(it, deformed, 1);
  return it;
}

/**
 * Set the displayed shape of an item to base + scale · (deformed − base).
 * @param {object} item render item from `buildItem`
 * @param {ArrayLike<number>|number[][]|null} deformed deformed node coordinates [m] (absolute, one per model
 *   node), or null to restore the undeformed shape
 * @param {number|'auto'} [scale] displacement magnification; 'auto' makes the largest displacement 5 % of the
 *   bounding-box diagonal
 * @returns {number} the scale applied
 */
export function applyDeformed(item, deformed, scale = 1) {
  if (!deformed) { item.positions = item.base; item.deformed = false; return 1; }
  const d = flat3(deformed, 'positions');
  if (d.length !== 3 * item.sourceCount) throw new RangeError(`positions: expected ${item.sourceCount} nodes, got ${d.length / 3}`);
  const n = item.nVerts, o = item.origin, src = item.source, zs = item.zScale;
  const disp = new Float32Array(3 * n);
  let dmax = 0;
  for (let i = 0; i < n; i++) {
    const s = src ? src[i] : i;
    const dx = d[3 * s] - o[0] - item.base[3 * i], dy = d[3 * s + 1] - o[1] - item.base[3 * i + 1], dz = d[3 * s + 2] * zs - o[2] - item.base[3 * i + 2];
    disp[3 * i] = dx; disp[3 * i + 1] = dy; disp[3 * i + 2] = dz;
    const m = dx * dx + dy * dy + dz * dz;
    if (m > dmax) dmax = m;
  }
  let k = scale;
  if (scale === 'auto') {
    const b = item.trueBounds;
    const diag = Math.hypot(b.max[0] - b.min[0], b.max[1] - b.min[1], b.max[2] - b.min[2]);
    k = dmax > 0 ? 0.05 * diag / Math.sqrt(dmax) : 1;
  }
  if (!Number.isFinite(k)) throw new RangeError('scale: must be a finite number or "auto"');
  const pos = new Float32Array(3 * n);
  for (let i = 0; i < 3 * n; i++) pos[i] = item.base[i] + k * disp[i];
  item.positions = pos;
  item.deformed = true;
  return k;
}

/** Mesh edges of an item, computed once on demand (volume face outlines, or the unique triangle edges). */
export function itemEdges(item) {
  if (item._edges === undefined) item._edges = item.presetEdges || (item.tris && item.tris.length ? uniqueEdges(triEdgePairs(item.tris), item.nVerts) : null);
  return item._edges;
}

/** Bounds of an item as displayed (with exaggeration and deformation), in world coordinates. */
export function itemBounds(item) {
  const b = boundsOf(item.positions);
  if (!b) return null;
  return { min: [b.min[0] + item.origin[0], b.min[1] + item.origin[1], b.min[2] + item.origin[2]], max: [b.max[0] + item.origin[0], b.max[1] + item.origin[1], b.max[2] + item.origin[2]] };
}

/** Union of boxes { min, max } (nulls ignored); null when there is nothing. */
export function unionBounds(list) {
  let out = null;
  for (const b of list) {
    if (!b) continue;
    if (!out) { out = { min: b.min.slice(), max: b.max.slice() }; continue; }
    for (let k = 0; k < 3; k++) { if (b.min[k] < out.min[k]) out.min[k] = b.min[k]; if (b.max[k] > out.max[k]) out.max[k] = b.max[k]; }
  }
  return out;
}

/**
 * Vertex streams for lines of constant screen width: every segment becomes two triangles whose corners are
 * pushed sideways in the vertex shader.
 * @param {Float32Array} positions packed vertex coordinates
 * @param {Uint32Array} pairs segment end vertex ids
 * @param {Float32Array} [t] colour-map coordinate per vertex
 * @returns {{ pos: Float32Array, other: Float32Array, side: Float32Array, t: Float32Array|null, count: number }}
 *   six vertices per segment: own position, the other end, side (±1) and colour coordinate.
 */
export function thickLineBuffers(positions, pairs, t = null) {
  const n = pairs.length >> 1;
  const pos = new Float32Array(18 * n), other = new Float32Array(18 * n), side = new Float32Array(6 * n);
  const tt = t ? new Float32Array(6 * n) : null;
  // corner order A+, A−, B+ / A−, B−, B+ ; at the B end the direction is reversed, so the side sign flips
  const END = [0, 0, 1, 0, 1, 1], SIDE = [1, -1, -1, -1, 1, -1];
  for (let s = 0; s < n; s++) {
    const a = pairs[2 * s], b = pairs[2 * s + 1];
    for (let k = 0; k < 6; k++) {
      const p = END[k] ? b : a, q = END[k] ? a : b, v = 6 * s + k;
      pos[3 * v] = positions[3 * p]; pos[3 * v + 1] = positions[3 * p + 1]; pos[3 * v + 2] = positions[3 * p + 2];
      other[3 * v] = positions[3 * q]; other[3 * v + 1] = positions[3 * q + 1]; other[3 * v + 2] = positions[3 * q + 2];
      side[v] = SIDE[k];
      if (tt) tt[v] = t[p];
    }
  }
  return { pos, other, side, t: tt, count: 6 * n };
}

/** Copy the indexed entries of a packed array into a flat, non-indexed one (`comps` values per vertex). */
export function expandIndexed(arr, comps, indices) {
  const out = new Float32Array(indices.length * comps);
  for (let i = 0; i < indices.length; i++) for (let k = 0; k < comps; k++) out[comps * i + k] = arr[comps * indices[i] + k];
  return out;
}

/**
 * Geometry of a highlight request, in item (origin-relative) coordinates.
 * @param {object} item render item
 * @param {{ edges?: number[][], elements?: number[], vertices?: number[] }} what vertex-id pairs, element ids
 *   (faces of a surface, cells of a volume, members of a frame, segments of a polyline) and vertex ids, all
 *   numbered as in the model given to the viewer
 * @returns {{ segments: Float32Array, triangles: Float32Array, points: Float32Array }} packed coordinates
 */
export function highlightGeometry(item, what = {}) {
  const P = item.positions;
  let lookup = null;
  if (item.source) { lookup = new Map(); for (let i = 0; i < item.source.length; i++) lookup.set(item.source[i], i); }
  const vid = (v) => (lookup ? lookup.get(v) : v < item.nVerts ? v : undefined);
  const seg = [], tri = [], pts = [];
  const pushV = (arr, v) => arr.push(P[3 * v], P[3 * v + 1], P[3 * v + 2]);
  for (const e of what.edges || []) {
    const a = vid(e[0]), b = vid(e[1]);
    if (a !== undefined && b !== undefined) { pushV(seg, a); pushV(seg, b); }
  }
  const els = what.elements ? new Set(what.elements) : null;
  if (els && els.size) {
    if (item.tris) {
      for (let t = 0; t < item.tris.length / 3; t++) {
        if (!els.has(item.triElement ? item.triElement[t] : t)) continue;
        pushV(tri, item.tris[3 * t]); pushV(tri, item.tris[3 * t + 1]); pushV(tri, item.tris[3 * t + 2]);
      }
    }
    if (item.wires && !item.tris) {
      for (let s = 0; s < item.wires.length / 2; s++) {
        if (!els.has(item.wireElement ? item.wireElement[s] : s)) continue;
        pushV(seg, item.wires[2 * s]); pushV(seg, item.wires[2 * s + 1]);
      }
    }
  }
  for (const v of what.vertices || []) { const i = vid(v); if (i !== undefined) pushV(pts, i); }
  return { segments: Float32Array.from(seg), triangles: Float32Array.from(tri), points: Float32Array.from(pts) };
}

/* -------------------------------------------------------------------------------------------------- picking */

/**
 * Nearest intersection of a ray with an indexed triangle mesh (both faces hit).
 * @returns {{ tri: number, t: number, point: number[] }|null} triangle index, ray parameter and hit point
 */
export function pickTriangles(positions, tris, origin, dir) {
  let best = Infinity, hit = -1;
  const ox = origin[0], oy = origin[1], oz = origin[2], dx = dir[0], dy = dir[1], dz = dir[2];
  for (let f = 0; f < tris.length; f += 3) {
    const a = 3 * tris[f], b = 3 * tris[f + 1], c = 3 * tris[f + 2];
    const e1x = positions[b] - positions[a], e1y = positions[b + 1] - positions[a + 1], e1z = positions[b + 2] - positions[a + 2];
    const e2x = positions[c] - positions[a], e2y = positions[c + 1] - positions[a + 1], e2z = positions[c + 2] - positions[a + 2];
    const px = dy * e2z - dz * e2y, py = dz * e2x - dx * e2z, pz = dx * e2y - dy * e2x;
    const det = e1x * px + e1y * py + e1z * pz;
    if (det === 0) continue;
    const inv = 1 / det;
    const tx = ox - positions[a], ty = oy - positions[a + 1], tz = oz - positions[a + 2];
    const u = (tx * px + ty * py + tz * pz) * inv;
    if (!(u >= 0 && u <= 1)) continue;
    const qx = ty * e1z - tz * e1y, qy = tz * e1x - tx * e1z, qz = tx * e1y - ty * e1x;
    const v = (dx * qx + dy * qy + dz * qz) * inv;
    if (!(v >= 0 && u + v <= 1)) continue;
    const t = (e2x * qx + e2y * qy + e2z * qz) * inv;
    if (t > 0 && t < best) { best = t; hit = f / 3; }
  }
  return hit < 0 ? null : { tri: hit, t: best, point: [ox + best * dx, oy + best * dy, oz + best * dz] };
}

/**
 * Nearest line segment or point to a screen position.
 * @param {Float32Array} positions packed coordinates
 * @param {Uint32Array|null} pairs segment vertex pairs, or null to test every vertex as a point
 * @param {Float64Array} mvp model-view-projection matrix
 * @param {number} nx cursor x in normalised device coordinates
 * @param {number} ny cursor y in normalised device coordinates
 * @param {number} tolX pick tolerance in NDC units along x
 * @param {number} tolY pick tolerance in NDC units along y
 * @returns {{ index: number, vertex: number, depth: number, point: number[] }|null} segment (or point) index,
 *   nearest end vertex, clip-space w of the hit (distance in front of the eye) and the hit position
 */
export function pickScreen(positions, pairs, mvp, nx, ny, tolX, tolY) {
  let best = 1, hit = null;
  const proj = (v) => projectPoint(mvp, positions[3 * v], positions[3 * v + 1], positions[3 * v + 2]);
  if (!pairs) {
    for (let v = 0; v < positions.length / 3; v++) {
      const p = proj(v);
      if (!(p[3] > 0)) continue;
      const d = ((p[0] - nx) / tolX) ** 2 + ((p[1] - ny) / tolY) ** 2;
      if (d < best) { best = d; hit = { index: v, vertex: v, depth: p[3], point: [positions[3 * v], positions[3 * v + 1], positions[3 * v + 2]] }; }
    }
    return hit;
  }
  for (let s = 0; s < pairs.length / 2; s++) {
    const a = pairs[2 * s], b = pairs[2 * s + 1];
    const p = proj(a), q = proj(b);
    if (!(p[3] > 0) || !(q[3] > 0)) continue;
    const ax = (p[0] - nx) / tolX, ay = (p[1] - ny) / tolY, bx = (q[0] - nx) / tolX, by = (q[1] - ny) / tolY;
    const ex = bx - ax, ey = by - ay;
    const l2 = ex * ex + ey * ey;
    const u = l2 > 0 ? Math.min(1, Math.max(0, -(ax * ex + ay * ey) / l2)) : 0;
    const cx = ax + u * ex, cy = ay + u * ey;
    const d = cx * cx + cy * cy;
    if (d >= best) continue;
    best = d;
    // perspective-correct parameter along the 3-D segment
    const w = u * p[3] / ((1 - u) * q[3] + u * p[3]) || 0;
    hit = {
      index: s, vertex: u < 0.5 ? a : b, depth: 1 / ((1 - u) / p[3] + u / q[3]),
      point: [positions[3 * a] + w * (positions[3 * b] - positions[3 * a]), positions[3 * a + 1] + w * (positions[3 * b + 1] - positions[3 * a + 1]), positions[3 * a + 2] + w * (positions[3 * b + 2] - positions[3 * a + 2])],
    };
  }
  return hit;
}

/* ------------------------------------------------------------------------------------------------- summary */

/** Bounding-box dimension text, e.g. "120 × 40 × 25 m". */
export function dimensionText(box) {
  if (!box) return '';
  return `${formatNumber(box.max[0] - box.min[0])} × ${formatNumber(box.max[1] - box.min[1])} × ${formatNumber(box.max[2] - box.min[2])} m`;
}

/**
 * Plain-language description of what is on screen, for assistive technology.
 * @param {object[]} items render items
 * @returns {string}
 */
export function describeScene(items) {
  if (!items.length) return 'Empty 3-D view.';
  const parts = items.map((it) => {
    const counts = Object.entries(it.counts || {}).map(([k, v]) => `${formatNumber(v, 6)} ${k}`).join(', ');
    let s = `${it.name}: ${it.kind === 'model3d' ? 'model' : it.kind}${counts ? ` with ${counts}` : ''}`;
    if (it.values && it.range && Number.isFinite(it.range.min)) s += `; coloured by ${it.valueLabel || 'value'}${it.valueUnit ? ` [${it.valueUnit}]` : ''} from ${formatNumber(it.range.min)} to ${formatNumber(it.range.max)}`;
    if (it.zScale !== 1) s += `; vertical scale exaggerated ${formatNumber(it.zScale)} times`;
    if (it.deformed) s += '; deformed shape shown';
    return s;
  });
  const box = unionBounds(items.map((it) => it.trueBounds));
  return `3-D view of ${items.length} model${items.length === 1 ? '' : 's'}. ${parts.join('. ')}. Overall size ${dimensionText(box)} (x × y × z, z up).`;
}

/* ------------------------------------------------------------------------------------------------- shaders */

const VS_MESH = `
attribute vec3 aPos; attribute vec3 aNormal; attribute float aT;
uniform mat4 uMVP; uniform mat4 uMV;
varying vec3 vView; varying vec3 vN; varying float vT;
void main() {
  vec4 mv = uMV * vec4(aPos, 1.0);
  vView = mv.xyz; vN = (uMV * vec4(aNormal, 0.0)).xyz; vT = aT;
  gl_Position = uMVP * vec4(aPos, 1.0);
}`;
const FS_MESH = `
uniform vec3 uColor; uniform vec3 uNoData; uniform float uUseMap; uniform float uFlat; uniform float uOpacity; uniform sampler2D uMap;
varying vec3 vView; varying vec3 vN; varying float vT;
void main() {
  vec3 n = vN;
#ifdef DERIV
  if (uFlat > 0.5) n = cross(dFdx(vView), dFdy(vView));
#endif
  float l = length(n);
  float d = l > 0.0 ? abs(dot(n / l, vec3(0.3345, 0.4301, 0.8385))) : 1.0;
  vec3 base = uColor;
  if (uUseMap > 0.5) base = vT < -0.5 ? uNoData : texture2D(uMap, vec2(clamp(vT, 0.0, 1.0), 0.5)).rgb;
  gl_FragColor = vec4(base * (0.42 + 0.58 * d), uOpacity);
}`;
const VS_LINE = `
attribute vec3 aPos;
uniform mat4 uMVP; uniform float uPointSize;
void main() { gl_Position = uMVP * vec4(aPos, 1.0); gl_PointSize = uPointSize; }`;
const FS_LINE = `
uniform vec4 uColor;
void main() { gl_FragColor = uColor; }`;
const VS_WIRE = `
attribute vec3 aPos; attribute vec3 aOther; attribute float aSide; attribute float aT;
uniform mat4 uMVP; uniform vec2 uViewport; uniform float uWidth;
varying float vT;
void main() {
  vec4 p = uMVP * vec4(aPos, 1.0);
  vec4 q = uMVP * vec4(aOther, 1.0);
  // keep the far end in front of the eye so the screen direction stays defined
  float eps = 1e-4 * max(abs(p.w), 1e-6);
  if (q.w < eps && p.w > eps) q = mix(p, q, (p.w - eps) / (p.w - q.w));
  vec2 d = (q.xy / max(q.w, eps) - p.xy / max(p.w, eps)) * uViewport;
  float l = length(d);
  d = l > 1e-6 ? d / l : vec2(1.0, 0.0);
  p.xy += vec2(-d.y, d.x) * aSide * 2.0 * uWidth / uViewport * p.w;
  vT = aT;
  gl_Position = p;
}`;
const FS_WIRE = `
uniform vec4 uColor; uniform vec3 uNoData; uniform float uUseMap; uniform sampler2D uMap;
varying float vT;
void main() {
  vec4 c = uColor;
  if (uUseMap > 0.5) c.rgb = vT < -0.5 ? uNoData : texture2D(uMap, vec2(clamp(vT, 0.0, 1.0), 0.5)).rgb;
  gl_FragColor = c;
}`;
const VS_POINT = `
attribute vec3 aPos; attribute float aT;
uniform mat4 uMVP; uniform float uPointSize;
varying float vT;
void main() { gl_Position = uMVP * vec4(aPos, 1.0); gl_PointSize = uPointSize; vT = aT; }`;
const FS_POINT = `
uniform vec4 uColor; uniform vec3 uNoData; uniform float uUseMap; uniform sampler2D uMap;
varying float vT;
void main() {
  vec2 c = gl_PointCoord - 0.5;
  if (dot(c, c) > 0.25) discard;
  vec4 col = uColor;
  if (uUseMap > 0.5) col.rgb = vT < -0.5 ? uNoData : texture2D(uMap, vec2(clamp(vT, 0.0, 1.0), 0.5)).rgb;
  gl_FragColor = col;
}`;

/** Adapt GLSL ES 1.00 source to the context version; `deriv` enables screen-space derivatives. */
export function shaderSource(src, type, webgl2, deriv) {
  const prec = '#ifdef GL_FRAGMENT_PRECISION_HIGH\nprecision highp float;\n#else\nprecision mediump float;\n#endif\n';
  if (webgl2) {
    if (type === 'vertex') return `#version 300 es\n#define attribute in\n#define varying out\n${src}`;
    return `#version 300 es\n${prec}${deriv ? '#define DERIV\n' : ''}#define varying in\n#define texture2D texture\nout vec4 fragColour;\n#define gl_FragColor fragColour\n${src}`;
  }
  if (type === 'vertex') return src;
  return `${deriv ? '#extension GL_OES_standard_derivatives : enable\n#define DERIV\n' : ''}${prec}${src}`;
}

/* -------------------------------------------------------------------------------------------------- viewer */

const STYLE_ID = 'v3d-style';
const STYLE = `
.v3d{position:absolute;inset:0;overflow:hidden;font:11px/1.3 system-ui,-apple-system,"Segoe UI",sans-serif;color:var(--ink,#111)}
.v3d canvas.v3d-c{position:absolute;inset:0;width:100%;height:100%;display:block;touch-action:none;cursor:grab}
.v3d canvas.v3d-c:active{cursor:grabbing}
.v3d canvas.v3d-c:focus-visible{outline:2px solid var(--focus,var(--accent,#06c));outline-offset:-2px}
.v3d-bar{position:absolute;top:6px;left:6px;display:flex;gap:4px;flex-wrap:wrap;max-width:calc(100% - 90px)}
.v3d-bar button{font:600 11px/1 inherit;font-family:inherit;padding:0 8px;min-width:34px;height:28px;border:1px solid var(--line,rgba(128,128,128,.5));border-radius:6px;background:var(--surface,#fff);color:var(--ink,#111);cursor:pointer}
.v3d-bar button[aria-pressed="true"]{background:var(--accent,#06c);color:var(--accent-ink,#fff);border-color:var(--accent,#06c)}
.v3d-dims{position:absolute;left:8px;bottom:6px;right:84px;pointer-events:none;color:var(--ink-2,var(--ink,#111));text-shadow:0 0 3px var(--surface-2,#fff)}
.v3d-legend{position:absolute;right:8px;top:8px;display:flex;flex-direction:column;align-items:flex-end;gap:2px;pointer-events:none;text-align:right;text-shadow:0 0 3px var(--surface-2,#fff)}
.v3d-legend i{display:block;width:12px;height:110px;border:1px solid var(--line,rgba(128,128,128,.5));border-radius:2px}
.v3d-legend b{font-weight:600;max-width:140px}
.v3d-triad{position:absolute;right:6px;bottom:6px;width:64px;height:64px;pointer-events:none}
.v3d-sr{position:absolute;width:1px;height:1px;margin:-1px;overflow:hidden;clip-path:inset(50%);white-space:nowrap}
.v3d-note{position:absolute;left:8px;bottom:22px;color:var(--ink-3,var(--ink,#111));pointer-events:none}
`;

let uid = 0;

function readTheme(el, win) {
  let cs = null;
  try { cs = win.getComputedStyle(el); } catch { cs = null; }
  const get = (name, fb) => parseColour(cs ? cs.getPropertyValue(name) : '', fb);
  const bg = get('--surface-2', [0.93, 0.95, 0.96]);
  const ink = get('--ink', [0.04, 0.09, 0.13]);
  const accent = get('--accent', [0.04, 0.43, 0.56]);
  const grid = get('--grid', [0.85, 0.88, 0.91]);
  const dark = luminance(bg) < 0.3;
  return {
    bg, ink, accent, grid, dark,
    // neutral body colour that keeps shading contrast on both themes
    mesh: mix(mix(bg, ink, dark ? 0.55 : 0.28), accent, 0.22),
    edge: mix(bg, ink, dark ? 0.9 : 0.78),
    wire: mix(ink, accent, 0.35),
    noData: mix(bg, ink, 0.35),
    highlight: dark ? [1, 0.62, 0.2] : [0.85, 0.2, 0.1],
  };
}

/**
 * Create a viewer inside `container` (the container is given `position: relative` when it is static and a
 * minimum height when it has none).
 *
 * @param {HTMLElement} container host element
 * @param {object} [opts] { edges?: boolean (default: on up to 150 000 triangles), view?: 'iso'|'top'|'front'|'side',
 *   ortho?: boolean, label?: string (accessible name), controls?: boolean (toolbar, default true),
 *   renderer?: 'auto'|'2d' (force the canvas fallback), maxPixelRatio?: number (default 2, never more than 2) }
 * @returns {object} viewer API:
 *   setScene({ models, styles? }) -> ids[]     replace everything (model: §4 model or model3d plot object)
 *   addModel(model, style?) -> id              style: { color, opacity, edges, lineWidth, pointSize, zScale,
 *                                              maxPoints, maxGridNodes, colourBy, visible, range, smooth }
 *   removeModel(id)
 *   setNodeValues(id, values, { label, unit, range: [min, max] }?)   values per model node; null clears
 *   setDeformed(id, positions, scale = 1) -> applied scale           positions per model node [m]; null restores
 *   highlight(id, { edges, elements, vertices }?)                    omit the second argument to clear
 *   clear(), fit(), setView(name), setEdges(on), setProjection('perspective'|'orthographic'),
 *   setVerticalScale(k) (same exaggeration for every model, e.g. a route over bathymetry), setStyle(id, style)
 *   onPick(cb) -> unsubscribe                  cb({ id, kind, name, point, element, node, value } | null)
 *   snapshot() -> PNG data URL, resize(), dispose(), getCamera(), setCamera(cam), summary() -> string
 *   canvas, mode ('webgl2' | 'webgl' | '2d')
 */
export function createViewer(container, opts = {}) {
  if (!container || !container.ownerDocument) throw new RangeError('container: an element is required');
  const doc = container.ownerDocument;
  const win = doc.defaultView;
  if (!doc.getElementById(STYLE_ID)) {
    const st = doc.createElement('style');
    st.id = STYLE_ID; st.textContent = STYLE;
    doc.head.append(st);
  }
  const el = (tag, cls, props) => { const e = doc.createElement(tag); if (cls) e.className = cls; if (props) Object.assign(e, props); return e; };
  if (win.getComputedStyle(container).position === 'static') container.style.position = 'relative';
  if (!container.clientHeight) container.style.minHeight = '280px';

  const myId = ++uid;
  const root = el('div', 'v3d');
  const canvas = el('canvas', 'v3d-c');
  const sr = el('p', 'v3d-sr', { id: `v3d-sum-${myId}` });
  canvas.tabIndex = 0;
  canvas.setAttribute('role', 'application');
  canvas.setAttribute('aria-label', opts.label || container.getAttribute('aria-label') || '3-D model view');
  canvas.setAttribute('aria-describedby', sr.id);
  canvas.setAttribute('aria-keyshortcuts', 'ArrowLeft ArrowRight ArrowUp ArrowDown + - F E 1 2 3 4');
  const bar = el('div', 'v3d-bar');
  const dims = el('div', 'v3d-dims');
  const note = el('div', 'v3d-note');
  const legend = el('div', 'v3d-legend');
  const SVG = 'http://www.w3.org/2000/svg';
  const triad = doc.createElementNS(SVG, 'svg');
  triad.setAttribute('class', 'v3d-triad'); triad.setAttribute('viewBox', '-32 -32 64 64'); triad.setAttribute('aria-hidden', 'true');
  root.append(canvas, bar, legend, dims, note, triad, sr);
  container.append(root);

  let theme = readTheme(container, win);
  let cam = defaultCamera();
  cam.ortho = !!opts.ortho;
  if (opts.view) cam.orientation = viewOrientation(opts.view);
  const entries = new Map();
  let nextId = 1;
  let edgesOn = opts.edges;
  let raf = 0, disposed = false, lost = false;
  let pxW = 1, pxH = 1, dpr = 1;
  const pickers = new Set();
  const cancelFrame = win.cancelAnimationFrame ? win.cancelAnimationFrame.bind(win) : clearTimeout;
  const requestFrame = win.requestAnimationFrame ? win.requestAnimationFrame.bind(win) : (f) => setTimeout(f, 16);

  /* ---- toolbar ---- */
  const buttons = {};
  const addButton = (key, text, title, fn) => {
    const b = el('button', '', { type: 'button', textContent: text, title });
    b.setAttribute('aria-label', title);
    b.addEventListener('click', () => { fn(); canvas.focus({ preventScroll: true }); });
    bar.append(b); buttons[key] = b;
  };
  if (opts.controls !== false) {
    addButton('fit', 'Fit', 'Fit the model to the view (F)', () => api.fit());
    addButton('iso', 'Iso', 'Isometric view (1)', () => api.setView('iso'));
    addButton('top', 'Top', 'Top view, looking down (2)', () => api.setView('top'));
    addButton('front', 'Front', 'Front view, looking along +y (3)', () => api.setView('front'));
    addButton('side', 'Side', 'Side view, looking along −x (4)', () => api.setView('side'));
    addButton('edges', 'Edges', 'Show or hide mesh edges (E)', () => api.setEdges(!edgesShown()));
  }

  /* ---- renderer ---- */
  let gl = null, mode = '2d', ctx2d = null, G = null;
  const glAttrs = { antialias: true, alpha: false, depth: true, stencil: false, premultipliedAlpha: false, preserveDrawingBuffer: false, powerPreference: 'default' };
  if (opts.renderer !== '2d') {
    try { gl = canvas.getContext('webgl2', glAttrs); if (gl) mode = 'webgl2'; } catch { gl = null; }
    if (!gl) {
      try { gl = canvas.getContext('webgl', glAttrs) || canvas.getContext('experimental-webgl', glAttrs); if (gl) mode = 'webgl'; } catch { gl = null; }
    }
  }
  if (gl) {
    try { initGl(); } catch {
      // the context exists but cannot be used (driver blacklist, shader failure): start again on a fresh canvas
      root.remove();
      return createViewer(container, { ...opts, renderer: '2d' });
    }
  } else {
    ctx2d = canvas.getContext('2d');
    note.textContent = 'Wireframe view (3-D graphics acceleration is not available).';
  }

  function compile(type, src) {
    const s = gl.createShader(type === 'vertex' ? gl.VERTEX_SHADER : gl.FRAGMENT_SHADER);
    gl.shaderSource(s, shaderSource(src, type, mode === 'webgl2', G.deriv));
    gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS) && !gl.isContextLost()) throw new Error(`shader: ${gl.getShaderInfoLog(s)}`);
    return s;
  }
  function program(vs, fs, attrs, uniforms) {
    const p = gl.createProgram();
    gl.attachShader(p, compile('vertex', vs)); gl.attachShader(p, compile('fragment', fs));
    gl.bindAttribLocation(p, 0, 'aPos');
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS) && !gl.isContextLost()) throw new Error(`program: ${gl.getProgramInfoLog(p)}`);
    const a = {}, u = {};
    for (const n of attrs) a[n] = gl.getAttribLocation(p, n);
    for (const n of uniforms) u[n] = gl.getUniformLocation(p, n);
    return { p, a, u };
  }
  function initGl() {
    G = { deriv: mode === 'webgl2' || !!gl.getExtension('OES_standard_derivatives'), uint: mode === 'webgl2' || !!gl.getExtension('OES_element_index_uint') };
    G.mesh = program(VS_MESH, FS_MESH, ['aPos', 'aNormal', 'aT'], ['uMVP', 'uMV', 'uColor', 'uNoData', 'uUseMap', 'uFlat', 'uOpacity', 'uMap']);
    G.line = program(VS_LINE, FS_LINE, ['aPos'], ['uMVP', 'uColor', 'uPointSize']);
    G.wire = program(VS_WIRE, FS_WIRE, ['aPos', 'aOther', 'aSide', 'aT'], ['uMVP', 'uViewport', 'uWidth', 'uColor', 'uNoData', 'uUseMap', 'uMap']);
    G.point = program(VS_POINT, FS_POINT, ['aPos', 'aT'], ['uMVP', 'uPointSize', 'uColor', 'uNoData', 'uUseMap', 'uMap']);
    G.map = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, G.map);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 256, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, colourMapTexture(256));
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    G.scratch = gl.createBuffer();
    G.box = gl.createBuffer();
  }
  const abuf = (data) => { const b = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, b); gl.bufferData(gl.ARRAY_BUFFER, data, gl.STATIC_DRAW); return b; };
  const ibuf = (data) => { const b = gl.createBuffer(); gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, b); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, data, gl.STATIC_DRAW); return b; };
  const smallIndex = (idx) => (G.uint ? { data: idx, type: gl.UNSIGNED_INT } : { data: Uint16Array.from(idx), type: gl.UNSIGNED_SHORT });

  function tNorm(e) {
    const it = e.item;
    if (!it.values) return null;
    const r = e.style.range ? { min: e.style.range[0], max: e.style.range[1] } : it.range;
    return normaliseValues(it.values, r);
  }

  function upload(e) {
    const it = e.item;
    const g = { bufs: [] };
    const keep = (b) => { g.bufs.push(b); return b; };
    const t = tNorm(e);
    // without 32-bit indices a large mesh is drawn from expanded, non-indexed streams
    const expand = !G.uint && it.nVerts > 65535;
    const needN = !G.deriv || e.style.smooth;
    if (it.tris && it.tris.length) {
      const nrm = needN ? vertexNormals(it.positions, it.tris) : null;
      if (expand) {
        g.triPos = keep(abuf(expandIndexed(it.positions, 3, it.tris)));
        if (nrm) g.triN = keep(abuf(expandIndexed(nrm, 3, it.tris)));
        if (t) g.triT = keep(abuf(expandIndexed(t, 1, it.tris)));
      } else {
        g.pos = keep(abuf(it.positions));
        g.triPos = g.pos;
        if (nrm) g.triN = keep(abuf(nrm));
        if (t) g.triT = keep(abuf(t));
        const ix = smallIndex(it.tris);
        g.triIdx = keep(ibuf(ix.data)); g.idxType = ix.type;
      }
      g.triCount = it.tris.length;
    }
    if (it.wires && it.wires.length) {
      const w = thickLineBuffers(it.positions, it.wires, t);
      g.wire = { pos: keep(abuf(w.pos)), other: keep(abuf(w.other)), side: keep(abuf(w.side)), t: w.t ? keep(abuf(w.t)) : null, count: w.count };
    }
    if (it.isPoints) {
      g.ptPos = keep(abuf(it.positions));
      if (t) g.ptT = keep(abuf(t));
      g.ptCount = it.nVerts;
    }
    if (it.deformed) {
      // undeformed outline drawn faintly behind the deformed shape
      const pairs = it.wires && it.wires.length ? it.wires : itemEdges(it);
      if (pairs && pairs.length) { g.ghost = keep(abuf(expandIndexed(it.base, 3, pairs))); g.ghostCount = pairs.length; }
    }
    g.expand = expand;
    e.gpu = g;
  }
  function uploadEdges(e) {
    const g = e.gpu, it = e.item;
    if (g.edgeReady) return;
    g.edgeReady = true;
    const ed = itemEdges(it);
    if (!ed || !ed.length) return;
    if (g.expand || !g.pos) { g.edgePos = abuf(expandIndexed(it.positions, 3, ed)); g.bufs.push(g.edgePos); g.edgeCount = ed.length; g.edgeIndexed = false; }
    else { const ix = smallIndex(ed); g.edgeIdx = ibuf(ix.data); g.bufs.push(g.edgeIdx); g.edgeCount = ed.length; g.edgeIndexed = true; }
  }
  function uploadHighlight(e) {
    const h = e.hl, g = { bufs: [] };
    if (h.segments.length) {
      const n = h.segments.length / 6;
      const pairs = new Uint32Array(2 * n);
      for (let i = 0; i < 2 * n; i++) pairs[i] = i;
      const w = thickLineBuffers(h.segments, pairs, null);
      g.wire = { pos: abuf(w.pos), other: abuf(w.other), side: abuf(w.side), t: null, count: w.count };
      g.bufs.push(g.wire.pos, g.wire.other, g.wire.side);
    }
    if (h.triangles.length) { g.tri = abuf(h.triangles); g.triCount = h.triangles.length / 3; g.bufs.push(g.tri); }
    if (h.points.length) { g.pts = abuf(h.points); g.ptCount = h.points.length / 3; g.bufs.push(g.pts); }
    e.hlGpu = g;
  }
  function freeGpu(e) {
    for (const key of ['gpu', 'hlGpu']) {
      if (e[key] && gl && !lost) for (const b of e[key].bufs) gl.deleteBuffer(b);
      e[key] = null;
    }
  }

  let enabled = [];
  const attr = (loc, buf, size) => {
    if (loc < 0) return;
    gl.bindBuffer(gl.ARRAY_BUFFER, buf);
    gl.enableVertexAttribArray(loc);
    gl.vertexAttribPointer(loc, size, gl.FLOAT, false, 0, 0);
    enabled.push(loc);
  };
  const doneAttrs = () => { for (const l of enabled) gl.disableVertexAttribArray(l); enabled = []; };
  const f32 = (m) => new Float32Array(m);

  function sceneBox() { return unionBounds([...entries.values()].filter((e) => e.style.visible !== false).map((e) => itemBounds(e.item))); }
  function trueBox() { return unionBounds([...entries.values()].map((e) => e.item.trueBounds)); }
  function triTotal() { let n = 0; for (const e of entries.values()) if (e.item.tris) n += e.item.tris.length / 3; return n; }
  function edgesShown() { return edgesOn ?? triTotal() <= 150000; }
  function edgesFor(e) { return e.style.edges ?? (edgesShown() && e.item.kind !== 'grid'); }

  function frameState() {
    const box = sceneBox();
    const aspect = pxW / pxH;
    const centre = box ? [(box.min[0] + box.max[0]) / 2, (box.min[1] + box.max[1]) / 2, (box.min[2] + box.max[2]) / 2] : cam.target;
    const radius = box ? Math.hypot(box.max[0] - box.min[0], box.max[1] - box.min[1], box.max[2] - box.min[2]) / 2 || 1 : 1;
    return { box, aspect, centre, radius };
  }

  function drawGl() {
    const fs = frameState();
    gl.viewport(0, 0, pxW, pxH);
    gl.clearColor(theme.bg[0], theme.bg[1], theme.bg[2], 1);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.enable(gl.DEPTH_TEST);
    gl.depthFunc(gl.LEQUAL);
    gl.disable(gl.CULL_FACE);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, G.map);
    const list = [...entries.values()].filter((e) => e.style.visible !== false);
    for (const e of list) {
      if (!e.gpu) upload(e);
      const m = cameraMatrices(cam, fs.aspect, fs.radius, e.item.origin, fs.centre);
      e.mvp = m.viewProj; e.mv32 = f32(m.view); e.mvp32 = f32(m.viewProj);
    }
    const meshPass = (transparent) => {
      gl.useProgram(G.mesh.p);
      gl.uniform1i(G.mesh.u.uMap, 0);
      gl.uniform3fv(G.mesh.u.uNoData, theme.noData);
      for (const e of list) {
        const g = e.gpu, op = e.style.opacity ?? 1;
        if (!g.triCount || (op < 1) !== transparent) continue;
        gl.uniformMatrix4fv(G.mesh.u.uMVP, false, e.mvp32);
        gl.uniformMatrix4fv(G.mesh.u.uMV, false, e.mv32);
        gl.uniform3fv(G.mesh.u.uColor, e.style.color ? parseColour(e.style.color, theme.mesh) : theme.mesh);
        gl.uniform1f(G.mesh.u.uUseMap, g.triT ? 1 : 0);
        gl.uniform1f(G.mesh.u.uFlat, g.triN && e.style.smooth ? 0 : 1);
        gl.uniform1f(G.mesh.u.uOpacity, op);
        attr(G.mesh.a.aPos, g.triPos, 3);
        if (g.triN) attr(G.mesh.a.aNormal, g.triN, 3); else if (G.mesh.a.aNormal >= 0) gl.vertexAttrib3f(G.mesh.a.aNormal, 0, 0, 1);
        if (g.triT) attr(G.mesh.a.aT, g.triT, 1); else if (G.mesh.a.aT >= 0) gl.vertexAttrib1f(G.mesh.a.aT, 0);
        if (g.triIdx) { gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, g.triIdx); gl.drawElements(gl.TRIANGLES, g.triCount, g.idxType, 0); }
        else gl.drawArrays(gl.TRIANGLES, 0, g.triCount);
        doneAttrs();
      }
    };
    // faces are pushed back slightly so that edges and highlights drawn on them are not z-fighting
    gl.enable(gl.POLYGON_OFFSET_FILL);
    gl.polygonOffset(1, 1);
    gl.disable(gl.BLEND);
    meshPass(false);
    gl.disable(gl.POLYGON_OFFSET_FILL);

    gl.useProgram(G.line.p);
    gl.uniform1f(G.line.u.uPointSize, 1);
    const lines = (mvp32, buf, count, colour, idx, idxType) => {
      gl.uniformMatrix4fv(G.line.u.uMVP, false, mvp32);
      gl.uniform4f(G.line.u.uColor, colour[0], colour[1], colour[2], colour[3] ?? 1);
      attr(G.line.a.aPos, buf, 3);
      if (idx) { gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, idx); gl.drawElements(gl.LINES, count, idxType, 0); } else gl.drawArrays(gl.LINES, 0, count);
      doneAttrs();
    };
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
    for (const e of list) {
      const g = e.gpu;
      if (g.triCount && edgesFor(e)) {
        uploadEdges(e);
        const c = e.style.edgeColor ? parseColour(e.style.edgeColor, theme.edge) : theme.edge;
        if (g.edgeCount) lines(e.mvp32, g.edgeIndexed ? g.pos : g.edgePos, g.edgeCount, [c[0], c[1], c[2], g.triT ? 0.35 : 0.55], g.edgeIndexed ? g.edgeIdx : null, g.idxType);
      }
      if (g.ghost) lines(e.mvp32, g.ghost, g.ghostCount, [theme.ink[0], theme.ink[1], theme.ink[2], 0.28]);
    }
    if (fs.box && entries.size) {
      const b = fs.box, o = fs.centre;
      const x0 = b.min[0] - o[0], y0 = b.min[1] - o[1], z0 = b.min[2] - o[2], x1 = b.max[0] - o[0], y1 = b.max[1] - o[1], z1 = b.max[2] - o[2];
      const c = [[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0], [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]];
      const seg = [0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6, 6, 7, 7, 4, 0, 4, 1, 5, 2, 6, 3, 7];
      const data = new Float32Array(72);
      seg.forEach((v, i) => data.set(c[v], 3 * i));
      gl.bindBuffer(gl.ARRAY_BUFFER, G.box);
      gl.bufferData(gl.ARRAY_BUFFER, data, gl.DYNAMIC_DRAW);
      const m = cameraMatrices(cam, fs.aspect, fs.radius, o, fs.centre);
      const gc = mix(theme.grid, theme.ink, 0.25);
      lines(f32(m.viewProj), G.box, 24, [gc[0], gc[1], gc[2], 1]);
    }

    const wirePass = (e, w, colour, width, useMap) => {
      gl.uniformMatrix4fv(G.wire.u.uMVP, false, e.mvp32);
      gl.uniform2f(G.wire.u.uViewport, pxW, pxH);
      gl.uniform1f(G.wire.u.uWidth, width * dpr);
      gl.uniform4f(G.wire.u.uColor, colour[0], colour[1], colour[2], colour[3] ?? 1);
      gl.uniform1f(G.wire.u.uUseMap, useMap && w.t ? 1 : 0);
      attr(G.wire.a.aPos, w.pos, 3); attr(G.wire.a.aOther, w.other, 3); attr(G.wire.a.aSide, w.side, 1);
      if (w.t && useMap) attr(G.wire.a.aT, w.t, 1); else if (G.wire.a.aT >= 0) gl.vertexAttrib1f(G.wire.a.aT, 0);
      gl.drawArrays(gl.TRIANGLES, 0, w.count);
      doneAttrs();
    };
    gl.useProgram(G.wire.p);
    gl.uniform1i(G.wire.u.uMap, 0);
    gl.uniform3fv(G.wire.u.uNoData, theme.noData);
    for (const e of list) {
      if (!e.gpu.wire) continue;
      const c = e.style.color ? parseColour(e.style.color, theme.wire) : theme.wire;
      wirePass(e, e.gpu.wire, [c[0], c[1], c[2], e.style.opacity ?? 1], (e.style.lineWidth ?? (e.item.tris ? 1.25 : 1.75)) / 2 + 0.25, true);
    }

    gl.useProgram(G.point.p);
    gl.uniform1i(G.point.u.uMap, 0);
    gl.uniform3fv(G.point.u.uNoData, theme.noData);
    for (const e of list) {
      const g = e.gpu;
      if (!g.ptCount) continue;
      const c = e.style.color ? parseColour(e.style.color, theme.wire) : theme.wire;
      gl.uniformMatrix4fv(G.point.u.uMVP, false, e.mvp32);
      gl.uniform1f(G.point.u.uPointSize, (e.style.pointSize ?? (g.ptCount > 20000 ? 2 : g.ptCount > 2000 ? 3 : 5)) * dpr);
      gl.uniform4f(G.point.u.uColor, c[0], c[1], c[2], 1);
      gl.uniform1f(G.point.u.uUseMap, g.ptT ? 1 : 0);
      attr(G.point.a.aPos, g.ptPos, 3);
      if (g.ptT) attr(G.point.a.aT, g.ptT, 1); else if (G.point.a.aT >= 0) gl.vertexAttrib1f(G.point.a.aT, 0);
      gl.drawArrays(gl.POINTS, 0, g.ptCount);
      doneAttrs();
    }

    gl.depthMask(false);
    meshPass(true);
    gl.depthMask(true);

    // highlights are drawn over everything so that problems hidden inside the model are still found
    gl.disable(gl.DEPTH_TEST);
    const hc = theme.highlight;
    for (const e of list) {
      if (!e.hl) continue;
      if (!e.hlGpu) uploadHighlight(e);
      const h = e.hlGpu;
      if (h.triCount) {
        gl.useProgram(G.line.p);
        gl.uniformMatrix4fv(G.line.u.uMVP, false, e.mvp32);
        gl.uniform4f(G.line.u.uColor, hc[0], hc[1], hc[2], 0.55);
        attr(G.line.a.aPos, h.tri, 3);
        gl.drawArrays(gl.TRIANGLES, 0, h.triCount);
        doneAttrs();
      }
      if (h.wire) { gl.useProgram(G.wire.p); wirePass(e, h.wire, [hc[0], hc[1], hc[2], 1], 1.75, false); }
      if (h.ptCount) {
        gl.useProgram(G.point.p);
        gl.uniformMatrix4fv(G.point.u.uMVP, false, e.mvp32);
        gl.uniform1f(G.point.u.uPointSize, 9 * dpr);
        gl.uniform4f(G.point.u.uColor, hc[0], hc[1], hc[2], 1);
        gl.uniform1f(G.point.u.uUseMap, 0);
        attr(G.point.a.aPos, h.pts, 3);
        if (G.point.a.aT >= 0) gl.vertexAttrib1f(G.point.a.aT, 0);
        gl.drawArrays(gl.POINTS, 0, h.ptCount);
        doneAttrs();
      }
    }
    gl.enable(gl.DEPTH_TEST);
    gl.disable(gl.BLEND);
  }

  function draw2d() {
    const fs = frameState();
    const c = ctx2d;
    c.setTransform(1, 0, 0, 1, 0, 0);
    c.fillStyle = css(theme.bg);
    c.fillRect(0, 0, pxW, pxH);
    c.lineCap = 'round';
    const MAX_SEG = 40000;
    const path = (mvp, P, pairs, stride) => {
      c.beginPath();
      for (let s = 0; s < pairs.length / 2; s += stride) {
        const a = pairs[2 * s], b = pairs[2 * s + 1];
        const p = projectPoint(mvp, P[3 * a], P[3 * a + 1], P[3 * a + 2]), q = projectPoint(mvp, P[3 * b], P[3 * b + 1], P[3 * b + 2]);
        if (!(p[3] > 0) || !(q[3] > 0)) continue;
        c.moveTo((p[0] + 1) / 2 * pxW, (1 - p[1]) / 2 * pxH);
        c.lineTo((q[0] + 1) / 2 * pxW, (1 - q[1]) / 2 * pxH);
      }
      c.stroke();
    };
    for (const e of entries.values()) {
      if (e.style.visible === false) continue;
      const it = e.item;
      const m = cameraMatrices(cam, fs.aspect, fs.radius, it.origin, fs.centre);
      e.mvp = m.viewProj;
      const ed = itemEdges(it);
      if (ed && ed.length) {
        c.strokeStyle = css(theme.edge, 0.6); c.lineWidth = dpr;
        path(m.viewProj, it.positions, ed, Math.max(1, Math.ceil(ed.length / 2 / MAX_SEG)));
      }
      if (it.wires && it.wires.length) {
        c.strokeStyle = css(e.style.color ? parseColour(e.style.color, theme.wire) : theme.wire); c.lineWidth = 1.75 * dpr;
        path(m.viewProj, it.positions, it.wires, Math.max(1, Math.ceil(it.wires.length / 2 / MAX_SEG)));
      }
      if (it.isPoints) {
        const st = Math.max(1, Math.ceil(it.nVerts / MAX_SEG));
        const sz = Math.max(2, 2 * dpr);
        const t = tNorm(e);
        c.fillStyle = css(theme.wire);
        for (let v = 0; v < it.nVerts; v += st) {
          const p = projectPoint(m.viewProj, it.positions[3 * v], it.positions[3 * v + 1], it.positions[3 * v + 2]);
          if (!(p[3] > 0)) continue;
          if (t) c.fillStyle = t[v] < 0 ? css(theme.noData) : css(colourMap(t[v]));
          c.fillRect((p[0] + 1) / 2 * pxW - sz / 2, (1 - p[1]) / 2 * pxH - sz / 2, sz, sz);
        }
      }
      if (e.hl) {
        c.strokeStyle = css(theme.highlight); c.fillStyle = css(theme.highlight, 0.5); c.lineWidth = 3 * dpr;
        const seq = (n) => { const a = new Uint32Array(n); for (let i = 0; i < n; i++) a[i] = i; return a; };
        if (e.hl.segments.length) path(m.viewProj, e.hl.segments, seq(e.hl.segments.length / 3), 1);
        if (e.hl.triangles.length) path(m.viewProj, e.hl.triangles, triEdgePairs(seq(e.hl.triangles.length / 3)), 1);
        for (let v = 0; v < e.hl.points.length / 3; v++) {
          const p = projectPoint(m.viewProj, e.hl.points[3 * v], e.hl.points[3 * v + 1], e.hl.points[3 * v + 2]);
          if (p[3] > 0) { c.beginPath(); c.arc((p[0] + 1) / 2 * pxW, (1 - p[1]) / 2 * pxH, 4 * dpr, 0, 2 * Math.PI); c.fill(); }
        }
      }
    }
  }

  function drawTriad() {
    const ax = cameraAxes(cam);
    const items = [['X', [1, 0, 0], '#d55e00'], ['Y', [0, 1, 0], '#009e73'], ['Z', [0, 0, 1], '#0072b2']].map(([n, v, col]) => ({
      n, col, x: v[0] * ax.right[0] + v[1] * ax.right[1] + v[2] * ax.right[2], y: v[0] * ax.up[0] + v[1] * ax.up[1] + v[2] * ax.up[2], z: v[0] * ax.back[0] + v[1] * ax.back[1] + v[2] * ax.back[2],
    })).sort((a, b) => a.z - b.z);
    while (triad.firstChild) triad.firstChild.remove();
    for (const a of items) {
      const ln = doc.createElementNS(SVG, 'line');
      ln.setAttribute('x1', 0); ln.setAttribute('y1', 0); ln.setAttribute('x2', (a.x * 20).toFixed(1)); ln.setAttribute('y2', (-a.y * 20).toFixed(1));
      ln.setAttribute('stroke', a.col); ln.setAttribute('stroke-width', 2.5); ln.setAttribute('stroke-linecap', 'round');
      const tx = doc.createElementNS(SVG, 'text');
      tx.setAttribute('x', (a.x * 27).toFixed(1)); tx.setAttribute('y', (-a.y * 27 + 3.5).toFixed(1));
      tx.setAttribute('text-anchor', 'middle'); tx.setAttribute('font-size', 10); tx.setAttribute('font-weight', 700);
      tx.setAttribute('fill', css(theme.ink));
      tx.textContent = a.n;
      triad.append(ln, tx);
    }
  }

  function render() {
    raf = 0;
    if (disposed) return;
    syncSize();
    if (gl) { if (lost || gl.isContextLost()) return; drawGl(); } else draw2d();
    drawTriad();
  }
  function invalidate() { if (!raf && !disposed) raf = requestFrame(render); }

  function syncSize() {
    const c = canvas;
    dpr = Math.min(win.devicePixelRatio || 1, Math.min(2, opts.maxPixelRatio || 2));
    const w = Math.max(1, Math.round(c.clientWidth * dpr)), h = Math.max(1, Math.round(c.clientHeight * dpr));
    if (c.width !== w || c.height !== h) { c.width = w; c.height = h; }
    pxW = w; pxH = h;
  }

  /* ---- overlay text ---- */
  function legendEntry() {
    for (const e of entries.values()) if (e.item.values && e.style.visible !== false) return e;
    return null;
  }
  function legendInfo() {
    const e = legendEntry();
    if (!e) return null;
    const r = e.style.range ? { min: e.style.range[0], max: e.style.range[1] } : e.item.range;
    if (!r || !Number.isFinite(r.min)) return null;
    return { label: e.item.valueLabel || 'Value', unit: e.item.valueUnit || '', min: r.min, max: r.max };
  }
  function updateOverlay() {
    const tb = trueBox();
    const items = [...entries.values()].map((e) => e.item);
    const ex = items.find((it) => it.zScale !== 1);
    dims.textContent = tb && entries.size ? `${dimensionText(tb)}${ex ? ` · vertical ×${formatNumber(ex.zScale)}` : ''}` : '';
    dims.title = 'Bounding box: x × y × z';
    while (legend.firstChild) legend.firstChild.remove();
    const li = legendInfo();
    if (li) {
      const title = el('b', '', { textContent: li.unit ? `${li.label} [${li.unit}]` : li.label });
      const hi = el('span', '', { textContent: formatNumber(li.max) });
      const grad = el('i');
      grad.style.background = `linear-gradient(to top, ${colourMapStops(9).join(',')})`;
      const lo = el('span', '', { textContent: formatNumber(li.min) });
      legend.append(title, hi, grad, lo);
    }
    sr.textContent = describeScene(items);
    if (buttons.edges) buttons.edges.setAttribute('aria-pressed', String(edgesShown()));
  }

  /* ---- interaction ---- */
  const live = canvas;
  const pointers = new Map();
  let drag = null;
  const rel = (ev) => { const r = live.getBoundingClientRect(); return { x: ev.clientX - r.left, y: ev.clientY - r.top, w: r.width || 1, h: r.height || 1 }; };
  const ball = (p) => { const s = Math.min(p.w, p.h); return arcballVector((2 * p.x - p.w) / s, (p.h - 2 * p.y) / s); };
  const two = () => { const [a, b] = [...pointers.values()]; return { cx: (a.x + b.x) / 2, cy: (a.y + b.y) / 2, d: Math.hypot(a.x - b.x, a.y - b.y) || 1 }; };

  function onDown(ev) {
    const p = rel(ev);
    pointers.set(ev.pointerId, p);
    try { live.setPointerCapture(ev.pointerId); } catch { /* synthetic events have no capture */ }
    if (pointers.size === 1) {
      const pan = ev.button === 1 || ev.button === 2 || ev.shiftKey || ev.ctrlKey;
      drag = { mode: pan ? 'pan' : 'orbit', q0: cam.orientation, v0: ball(p), x: p.x, y: p.y, x0: p.x, y0: p.y, moved: 0, t: Date.now(), single: true };
    } else if (pointers.size === 2) {
      drag = { mode: 'pinch', ...two(), moved: 99, single: false };
    }
    live.focus({ preventScroll: true });
    ev.preventDefault();
  }
  function onMove(ev) {
    if (!pointers.has(ev.pointerId) || !drag) return;
    const p = rel(ev);
    pointers.set(ev.pointerId, p);
    if (drag.mode === 'pinch' && pointers.size >= 2) {
      const t = two();
      cam = zoomCamera(cam, drag.d / t.d, (2 * t.cx - p.w) / p.w, (p.h - 2 * t.cy) / p.h, p.w / p.h);
      cam.target = panTarget(cam, t.cx - drag.cx, t.cy - drag.cy, p.h);
      drag.cx = t.cx; drag.cy = t.cy; drag.d = t.d;
    } else if (drag.mode === 'orbit') {
      drag.moved = Math.max(drag.moved, Math.hypot(p.x - drag.x0, p.y - drag.y0));
      if (drag.moved > 3) cam.orientation = orbitArcball(drag.q0, drag.v0, ball(p));
    } else if (drag.mode === 'pan') {
      drag.moved = Math.max(drag.moved, Math.hypot(p.x - drag.x0, p.y - drag.y0));
      cam.target = panTarget(cam, p.x - drag.x, p.y - drag.y, p.h);
      drag.x = p.x; drag.y = p.y;
    }
    invalidate();
  }
  function onUp(ev) {
    if (!pointers.has(ev.pointerId)) return;
    const p = pointers.get(ev.pointerId);
    pointers.delete(ev.pointerId);
    if (drag && drag.single && drag.moved <= 3 && Date.now() - drag.t < 600 && ev.type === 'pointerup' && pickers.size) {
      const hit = pickAt(p.x, p.y, p.w, p.h);
      for (const cb of pickers) cb(hit);
    }
    if (pointers.size === 1) {
      const q = [...pointers.values()][0];
      drag = { mode: 'orbit', q0: cam.orientation, v0: ball(q), x: q.x, y: q.y, x0: q.x, y0: q.y, moved: 99, t: 0, single: false };
    } else if (pointers.size === 0) drag = null;
  }
  function onWheel(ev) {
    ev.preventDefault();
    const p = rel(ev);
    const dy = ev.deltaMode === 1 ? ev.deltaY * 16 : ev.deltaMode === 2 ? ev.deltaY * p.h : ev.deltaY;
    cam = zoomCamera(cam, Math.exp(Math.max(-1, Math.min(1, dy * 0.0015))), (2 * p.x - p.w) / p.w, (p.h - 2 * p.y) / p.h, p.w / p.h);
    invalidate();
  }
  function onKey(ev) {
    if (ev.altKey || ev.metaKey || ev.ctrlKey) return;
    const step = 10 * Math.PI / 180;
    const h = live.clientHeight || 1;
    let used = true;
    const arrow = (dx, dy) => {
      if (ev.shiftKey) cam.target = panTarget(cam, dx * h * 0.1, dy * h * 0.1, h);
      else cam.orientation = orbitTurntable(cam.orientation, -dx * step, -dy * step);
    };
    switch (ev.key) {
      case 'ArrowLeft': arrow(-1, 0); break;
      case 'ArrowRight': arrow(1, 0); break;
      case 'ArrowUp': arrow(0, -1); break;
      case 'ArrowDown': arrow(0, 1); break;
      case '+': case '=': cam = zoomCamera(cam, 1 / 1.2); break;
      case '-': case '_': cam = zoomCamera(cam, 1.2); break;
      case 'f': case 'F': api.fit(); break;
      case 'e': case 'E': api.setEdges(!edgesShown()); break;
      case '1': case 'Home': api.setView('iso'); break;
      case '2': api.setView('top'); break;
      case '3': api.setView('front'); break;
      case '4': api.setView('side'); break;
      default: used = false;
    }
    if (used) { ev.preventDefault(); invalidate(); }
  }
  const noMenu = (ev) => ev.preventDefault();
  function onLost(ev) { ev.preventDefault(); lost = true; if (raf) { cancelFrame(raf); raf = 0; } note.textContent = '3-D graphics were interrupted; restoring…'; }
  function onRestored() {
    lost = false;
    for (const e of entries.values()) { e.gpu = null; e.hlGpu = null; }
    try { initGl(); note.textContent = ''; } catch { note.textContent = '3-D graphics could not be restored. Reload the page.'; return; }
    invalidate();
  }
  function onTheme() { theme = readTheme(container, win); invalidate(); }

  live.addEventListener('pointerdown', onDown);
  live.addEventListener('pointermove', onMove);
  live.addEventListener('pointerup', onUp);
  live.addEventListener('pointercancel', onUp);
  live.addEventListener('wheel', onWheel, { passive: false });
  live.addEventListener('keydown', onKey);
  live.addEventListener('contextmenu', noMenu);
  live.addEventListener('webglcontextlost', onLost);
  live.addEventListener('webglcontextrestored', onRestored);
  win.addEventListener('themechange', onTheme);
  let mq = null;
  try { mq = win.matchMedia('(prefers-color-scheme: dark)'); mq.addEventListener('change', onTheme); } catch { mq = null; }
  let ro = null;
  if (typeof win.ResizeObserver === 'function') { ro = new win.ResizeObserver(() => invalidate()); ro.observe(root); }
  else win.addEventListener('resize', invalidate);

  function pickAt(x, y, w, h) {
    const fs = frameState();
    const aspect = w / h;
    const nx = (2 * x - w) / w, ny = (h - 2 * y) / h;
    const ray = screenRay(cam, aspect, nx, ny);
    const { back } = cameraAxes(cam);
    let best = null;
    for (const [id, e] of entries) {
      if (e.style.visible === false) continue;
      const it = e.item;
      const o = [ray.origin[0] - it.origin[0], ray.origin[1] - it.origin[1], ray.origin[2] - it.origin[2]];
      const mvp = cameraMatrices(cam, aspect, fs.radius, it.origin, fs.centre).viewProj;
      // distance in front of the eye plane, comparable between faces, lines and points in both projections
      const depthOf = (pt) => -(pt[0] - o[0]) * back[0] - (pt[1] - o[1]) * back[1] - (pt[2] - o[2]) * back[2];
      const consider = (depth, res) => { if (!best || depth < best.depth) best = { depth, id, e, ...res }; };
      if (it.tris && it.tris.length) {
        const hit = pickTriangles(it.positions, it.tris, o, ray.dir);
        if (hit) {
          let v = it.tris[3 * hit.tri], dmin = Infinity;
          for (let k = 0; k < 3; k++) {
            const q = it.tris[3 * hit.tri + k];
            const d = Math.hypot(it.positions[3 * q] - hit.point[0], it.positions[3 * q + 1] - hit.point[1], it.positions[3 * q + 2] - hit.point[2]);
            if (d < dmin) { dmin = d; v = q; }
          }
          consider(depthOf(hit.point), { point: hit.point, element: it.triElement ? it.triElement[hit.tri] : hit.tri, vertex: v });
        }
      }
      const tolX = 16 / w, tolY = 16 / h;
      if (it.wires && it.wires.length) {
        const hit = pickScreen(it.positions, it.wires, mvp, nx, ny, tolX, tolY);
        if (hit) consider(depthOf(hit.point) * 0.999, { point: hit.point, element: it.wireElement ? it.wireElement[hit.index] : hit.index, vertex: hit.vertex });
      }
      if (it.isPoints) {
        const hit = pickScreen(it.positions, null, mvp, nx, ny, tolX, tolY);
        if (hit) consider(depthOf(hit.point) * 0.999, { point: hit.point, element: hit.index, vertex: hit.vertex });
      }
    }
    if (!best) return null;
    const it = best.e.item;
    const node = it.source ? it.source[best.vertex] : best.vertex;
    return {
      id: best.id, kind: it.kind, name: it.name,
      point: [best.point[0] + it.origin[0], best.point[1] + it.origin[1], (best.point[2] + it.origin[2]) / it.zScale],
      element: best.element, node, value: it.values ? it.values[best.vertex] : undefined,
    };
  }

  const need = (id) => { const e = entries.get(id); if (!e) throw new RangeError(`id: no model with id ${id}`); return e; };

  const api = {
    get canvas() { return live; },
    get mode() { return mode; },
    setScene(scene) {
      api.clear();
      const models = (scene && scene.models) || [];
      const ids = models.map((m, i) => api.addModel(m, scene.styles ? scene.styles[i] : undefined));
      const box = trueBox();
      // a scene with no height (cross-sections, plans) is best seen from above
      if (!opts.view && box && box.max[2] - box.min[2] <= 1e-9 * Math.max(box.max[0] - box.min[0], box.max[1] - box.min[1])) cam.orientation = viewOrientation('top');
      else if (!opts.view) cam.orientation = viewOrientation('iso');
      api.fit();
      return ids;
    },
    addModel(model, style = {}) {
      const id = nextId++;
      entries.set(id, { model, style: { ...style }, item: buildItem(model, style), gpu: null, hl: null, hlGpu: null });
      updateOverlay(); invalidate();
      return id;
    },
    removeModel(id) { const e = entries.get(id); if (!e) return; freeGpu(e); entries.delete(id); updateOverlay(); invalidate(); },
    setStyle(id, style) {
      const e = need(id);
      const rebuild = ['zScale', 'maxPoints', 'maxGridNodes', 'colourBy'].some((k) => k in style && style[k] !== e.style[k]);
      Object.assign(e.style, style);
      if (rebuild) e.item = buildItem(e.model, e.style);
      freeGpu(e); updateOverlay(); invalidate();
    },
    setNodeValues(id, values, info = {}) {
      const e = need(id), it = e.item;
      if (values == null) { it.values = null; it.range = null; e.style.range = undefined; }
      else {
        if (values.length !== it.sourceCount) throw new RangeError(`values: expected ${it.sourceCount} values (one per node), got ${values.length}`);
        const v = new Float32Array(it.nVerts);
        for (let i = 0; i < it.nVerts; i++) v[i] = values[it.source ? it.source[i] : i];
        it.values = v; it.range = valueRange(v);
        it.valueLabel = info.label ?? 'Value'; it.valueUnit = info.unit ?? '';
        const r = info.range;
        e.style.range = r ? (Array.isArray(r) ? [r[0], r[1]] : [r.min, r.max]) : undefined;
      }
      freeGpu(e); updateOverlay(); invalidate();
    },
    setDeformed(id, positions, scale = 1) {
      const e = need(id);
      const k = applyDeformed(e.item, positions, scale);
      if (e.hlSpec) e.hl = highlightGeometry(e.item, e.hlSpec);
      freeGpu(e); updateOverlay(); invalidate();
      return k;
    },
    highlight(id, what) {
      const e = need(id);
      if (e.hlGpu && gl && !lost) for (const b of e.hlGpu.bufs) gl.deleteBuffer(b);
      e.hlGpu = null;
      e.hlSpec = what && ((what.edges && what.edges.length) || (what.elements && what.elements.length) || (what.vertices && what.vertices.length)) ? what : null;
      e.hl = e.hlSpec ? highlightGeometry(e.item, e.hlSpec) : null;
      invalidate();
    },
    clear() { for (const e of entries.values()) freeGpu(e); entries.clear(); updateOverlay(); invalidate(); },
    fit() {
      syncSize();
      const box = sceneBox();
      if (box) cam = fitCamera(cam, box, pxW / pxH);
      invalidate();
    },
    setView(name) { cam.orientation = viewOrientation(name); api.fit(); },
    setEdges(on) { edgesOn = !!on; if (buttons.edges) buttons.edges.setAttribute('aria-pressed', String(edgesOn)); invalidate(); },
    setVerticalScale(k) {
      if (!(k > 0)) throw new RangeError('k: vertical scale must be > 0');
      for (const e of entries.values()) {
        e.style.zScale = k;
        const keep = e.item;
        e.item = buildItem(e.model, e.style);
        if (keep.values && !e.item.values) { e.item.values = keep.values; e.item.range = keep.range; }
        if (keep.values) { e.item.valueLabel = keep.valueLabel; e.item.valueUnit = keep.valueUnit; }
        if (e.hlSpec) e.hl = highlightGeometry(e.item, e.hlSpec);
        freeGpu(e);
      }
      updateOverlay(); api.fit();
    },
    setProjection(p) { cam.ortho = p === 'orthographic'; invalidate(); },
    getCamera() { return { ...cam, target: cam.target.slice(), orientation: cam.orientation.slice() }; },
    setCamera(c) { cam = { ...cam, ...c }; cam.orientation = quatNormalize(cam.orientation); invalidate(); },
    onPick(cb) { pickers.add(cb); return () => pickers.delete(cb); },
    pick(x, y) { const r = live.getBoundingClientRect(); return pickAt(x, y, r.width || 1, r.height || 1); },
    summary() { return sr.textContent; },
    snapshot() {
      if (raf) { cancelFrame(raf); raf = 0; }
      render();
      const out = doc.createElement('canvas');
      out.width = pxW; out.height = pxH;
      const c = out.getContext('2d');
      c.drawImage(live, 0, 0);
      c.font = `${11 * dpr}px system-ui, sans-serif`;
      c.fillStyle = css(theme.ink);
      c.textBaseline = 'alphabetic';
      if (dims.textContent) c.fillText(dims.textContent, 8 * dpr, pxH - 8 * dpr);
      const li = legendInfo();
      if (li) {
        const bw = 12 * dpr, bh = 110 * dpr, x = pxW - 8 * dpr - bw, y = 36 * dpr;
        const grd = c.createLinearGradient(0, y + bh, 0, y);
        colourMapStops(9).forEach((s, i) => grd.addColorStop(i / 8, s));
        c.fillStyle = grd; c.fillRect(x, y, bw, bh);
        c.fillStyle = css(theme.ink); c.textAlign = 'right';
        c.fillText(li.unit ? `${li.label} [${li.unit}]` : li.label, pxW - 8 * dpr, 16 * dpr);
        c.fillText(formatNumber(li.max), pxW - 8 * dpr, 30 * dpr);
        c.fillText(formatNumber(li.min), pxW - 8 * dpr, y + bh + 14 * dpr);
      }
      return out.toDataURL('image/png');
    },
    resize() { syncSize(); invalidate(); },
    invalidate,
    dispose() {
      if (disposed) return;
      disposed = true;
      if (raf) cancelFrame(raf);
      live.removeEventListener('pointerdown', onDown); live.removeEventListener('pointermove', onMove);
      live.removeEventListener('pointerup', onUp); live.removeEventListener('pointercancel', onUp);
      live.removeEventListener('wheel', onWheel); live.removeEventListener('keydown', onKey);
      live.removeEventListener('contextmenu', noMenu);
      live.removeEventListener('webglcontextlost', onLost); live.removeEventListener('webglcontextrestored', onRestored);
      win.removeEventListener('themechange', onTheme);
      if (mq) mq.removeEventListener('change', onTheme);
      if (ro) ro.disconnect(); else win.removeEventListener('resize', invalidate);
      for (const e of entries.values()) freeGpu(e);
      entries.clear(); pickers.clear();
      if (gl) { const x = gl.getExtension('WEBGL_lose_context'); if (x) x.loseContext(); }
      root.remove();
    },
  };
  updateOverlay();
  invalidate();
  return api;
}
