// Hydrostatics of closed triangulated surfaces clipped at an arbitrary flotation plane, closed-form primitives,
// equilibrium solver and free-surface corrections.
//
// Body frame: right-handed, x forward, y to port, z up; z = 0 is normally the keel/baseline.
// Floating position (draft, heel, trim):  p_earth = Ry(trim) · Rx(heel) · p_body - (0, 0, draft)
//   draft: depth of the body origin below the waterline, m
//   heel:  rotation about x, rad, positive = port side up (starboard down)
//   trim:  rotation about y, rad, positive = bow down
// The earth frame therefore has its origin on the waterline vertically above the (rotated) body origin.

import { solveDense } from '../math/linalg.js';
import { G, needPos } from './_wavehelp.js';

let _buf = new Float64Array(0);

function rotation(heel, trim) {
  const cf = Math.cos(heel); const sf = Math.sin(heel);
  const ct = Math.cos(trim); const st = Math.sin(trim);
  return [ct, st * sf, st * cf, 0, cf, -sf, -st, ct * sf, ct * cf];
}

/**
 * Body → earth transformation of a point for a floating position.
 * @param {number[]} p body coordinates, m
 * @param {number} draft m
 * @param {number} heel rad
 * @param {number} trim rad
 * @returns {number[]} earth coordinates, m (z relative to the waterline)
 */
export function bodyToEarth(p, draft, heel, trim) {
  const R = rotation(heel, trim);
  return [
    R[0] * p[0] + R[1] * p[1] + R[2] * p[2],
    R[3] * p[0] + R[4] * p[1] + R[5] * p[2],
    R[6] * p[0] + R[7] * p[1] + R[8] * p[2] - draft,
  ];
}

/**
 * Earth → body transformation (inverse of bodyToEarth).
 * @returns {number[]} body coordinates, m
 */
export function earthToBody(p, draft, heel, trim) {
  const R = rotation(heel, trim);
  const z = p[2] + draft;
  return [
    R[0] * p[0] + R[3] * p[1] + R[6] * z,
    R[1] * p[0] + R[4] * p[1] + R[7] * z,
    R[2] * p[0] + R[5] * p[1] + R[8] * z,
  ];
}

const KEYS = ['V', 'Mx', 'My', 'Mz', 'Aw', 'Sx', 'Sy', 'Jxx', 'Jyy', 'Jxy', 'wet'];

function emptyRaw() {
  return { V: 0, Mx: 0, My: 0, Mz: 0, Aw: 0, Sx: 0, Sy: 0, Jxx: 0, Jyy: 0, Jxy: 0, wet: 0 };
}

/**
 * Raw immersed-volume integrals of a closed triangulated surface below the flotation plane, in the earth
 * frame. Volume integrals use the divergence theorem with the field (0, 0, z·f), which vanishes on the
 * waterplane cap, so no explicit waterline polygon is needed; waterplane integrals follow from
 * ∮ f(x, y) n_z dS = 0 over the closed immersed surface. Exact for the faceted geometry.
 * The mesh must be closed; triangle orientation may be consistently outward or inward.
 * @param {{positions:ArrayLike<number>, indices:ArrayLike<number>}} mesh contract §4 surface, m
 * @param {number} draft m
 * @param {number} heel rad
 * @param {number} trim rad
 * @returns {{V:number, Mx:number, My:number, Mz:number, Aw:number, Sx:number, Sy:number, Jxx:number, Jyy:number,
 *   Jxy:number, wet:number}} volume m³; first moments of volume ∫x dV, ∫y dV, ∫z dV (m⁴); waterplane area m²;
 *   waterplane first moments ∫x dA, ∫y dA (m³); second moments Jxx = ∫y² dA, Jyy = ∫x² dA, Jxy = ∫xy dA about
 *   the earth origin (m⁴); wetted surface area m²
 */
export function clipIntegrals(mesh, draft, heel, trim) {
  const { positions, indices } = mesh;
  if (!positions || !indices || indices.length % 3 !== 0) throw new RangeError('mesh must have positions and triangle indices');
  const nv = positions.length / 3;
  if (_buf.length < positions.length) _buf = new Float64Array(positions.length);
  const P = _buf;
  const R = rotation(heel, trim);
  for (let i = 0; i < nv; i++) {
    const x = positions[3 * i]; const y = positions[3 * i + 1]; const z = positions[3 * i + 2];
    P[3 * i] = R[0] * x + R[1] * y + R[2] * z;
    P[3 * i + 1] = R[3] * x + R[4] * y + R[5] * z;
    P[3 * i + 2] = R[6] * x + R[7] * y + R[8] * z - draft;
  }
  let V = 0; let Mx = 0; let My = 0; let Mz = 0; let Aw = 0; let Sx = 0; let Sy = 0;
  let Jxx = 0; let Jyy = 0; let Jxy = 0; let wet = 0;
  const acc = (ax, ay, az, bx, by, bz, cx, cy, cz) => {
    const e1x = bx - ax; const e1y = by - ay; const e1z = bz - az;
    const e2x = cx - ax; const e2y = cy - ay; const e2z = cz - az;
    const nx = e1y * e2z - e1z * e2y; const ny = e1z * e2x - e1x * e2z; const nz = e1x * e2y - e1y * e2x;
    const dA = 0.5 * nz;
    const sx = ax + bx + cx; const sy = ay + by + cy; const sz = az + bz + cz;
    wet += 0.5 * Math.sqrt(nx * nx + ny * ny + nz * nz);
    V += (dA * sz) / 3;
    Mx += (dA * (ax * az + bx * bz + cx * cz + sx * sz)) / 12;
    My += (dA * (ay * az + by * bz + cy * cz + sy * sz)) / 12;
    Mz += (dA * (az * az + bz * bz + cz * cz + sz * sz)) / 24;
    Aw -= dA;
    Sx -= (dA * sx) / 3;
    Sy -= (dA * sy) / 3;
    Jxx -= (dA * (ay * ay + by * by + cy * cy + sy * sy)) / 12;
    Jyy -= (dA * (ax * ax + bx * bx + cx * cx + sx * sx)) / 12;
    Jxy -= (dA * (ax * ay + bx * by + cx * cy + sx * sy)) / 12;
  };
  const q = new Float64Array(12);
  for (let f = 0; f < indices.length; f += 3) {
    const ia = 3 * indices[f]; const ib = 3 * indices[f + 1]; const ic = 3 * indices[f + 2];
    const az = P[ia + 2]; const bz = P[ib + 2]; const cz = P[ic + 2];
    if (az > 0 && bz > 0 && cz > 0) continue;
    if (az <= 0 && bz <= 0 && cz <= 0) {
      acc(P[ia], P[ia + 1], az, P[ib], P[ib + 1], bz, P[ic], P[ic + 1], cz);
      continue;
    }
    // Sutherland–Hodgman clip of the triangle against z <= 0, keeping the winding.
    let n = 0;
    let pi = ic;
    for (let e = 0; e < 3; e++) {
      const ci = e === 0 ? ia : e === 1 ? ib : ic;
      const zp = P[pi + 2]; const zc = P[ci + 2];
      const inP = zp <= 0; const inC = zc <= 0;
      if (inP !== inC) {
        const s = zp / (zp - zc);
        q[3 * n] = P[pi] + s * (P[ci] - P[pi]);
        q[3 * n + 1] = P[pi + 1] + s * (P[ci + 1] - P[pi + 1]);
        q[3 * n + 2] = 0;
        n++;
      }
      if (inC) { q[3 * n] = P[ci]; q[3 * n + 1] = P[ci + 1]; q[3 * n + 2] = zc; n++; }
      pi = ci;
    }
    for (let k = 1; k + 1 < n; k++) {
      acc(q[0], q[1], q[2], q[3 * k], q[3 * k + 1], q[3 * k + 2], q[3 * k + 3], q[3 * k + 4], q[3 * k + 5]);
    }
  }
  if (V < 0) {
    V = -V; Mx = -Mx; My = -My; Mz = -Mz; Aw = -Aw; Sx = -Sx; Sy = -Sy; Jxx = -Jxx; Jyy = -Jyy; Jxy = -Jxy;
  }
  return { V, Mx, My, Mz, Aw, Sx, Sy, Jxx, Jyy, Jxy, wet };
}

/**
 * Immersed integrals of a body with lost buoyancy: intact hull minus permeability-weighted damaged
 * compartments, all clipped at the same flotation plane (lost-buoyancy method).
 * @param {object} mesh hull surface
 * @param {{mesh:object, permeability?:number}[]} [lost] damaged compartments (closed meshes in body coordinates)
 * @returns {object} raw integrals as clipIntegrals
 */
export function bodyIntegrals(mesh, lost, draft, heel, trim) {
  const r = clipIntegrals(mesh, draft, heel, trim);
  if (!lost || lost.length === 0) return r;
  const out = { ...r };
  for (const c of lost) {
    const mu = c.permeability === undefined ? 1 : c.permeability;
    if (!(mu >= 0 && mu <= 1)) throw new RangeError('permeability must be within 0..1');
    const rc = clipIntegrals(c.mesh, draft, heel, trim);
    for (const key of KEYS) if (key !== 'wet') out[key] -= mu * rc[key];
  }
  return out;
}

/**
 * Hydrostatic particulars of a closed triangulated surface floating at (draft, heel, trim).
 * KB, BM, KM, GM are the classical upright definitions generalised to the instantaneous waterplane:
 * BM = I/V with I about the waterplane centroidal axes parallel to earth x (transverse stability) and y
 * (longitudinal); GM = BM + z_B - z_G in the earth frame. Stiffness is the 6×6 linear hydrostatic restoring
 * matrix about the earth origin (surge, sway, heave, roll, pitch, yaw):
 *   C33 = rho g Aw, C34 = rho g ∫y dA, C35 = -rho g ∫x dA, C44 = rho g ∫y² dA + rho g V z_B - m g z_G,
 *   C45 = -rho g ∫xy dA, C55 = rho g ∫x² dA + rho g V z_B - m g z_G, C46 = -rho g V x_B + m g x_G,
 *   C56 = -rho g V y_B + m g y_G.
 * Reference: Newman (1977) Marine Hydrodynamics ch. 6 (hydrostatic restoring coefficients); Biran & López-Pulido (2014) Ship Hydrostatics and
 * Stability ch. 2. Exact for the faceted mesh; wall-sided linearisation for the stiffness.
 * @param {object} o
 * @param {{positions:ArrayLike<number>, indices:ArrayLike<number>}} o.mesh closed surface, m (body frame)
 * @param {number} o.draft m
 * @param {number} [o.heel=0] rad
 * @param {number} [o.trim=0] rad
 * @param {number} [o.rho=1025] kg/m³
 * @param {number} [o.g=9.80665] m/s²
 * @param {number[]} [o.cog] centre of gravity in body coordinates, m (enables GM and gravity stiffness terms)
 * @param {number} [o.mass] kg (default rho·V)
 * @param {number} [o.length] reference length for MCT, m (default: mesh extent in x)
 * @param {{mesh:object, permeability?:number}[]} [o.lost] damaged compartments (lost buoyancy)
 * @returns {object} volume m³; displacement kg; buoyancy N; cob (earth) and cobBody [x,y,z] m; wettedArea m²;
 *   waterplaneArea m²; cof (centre of flotation, earth) and cofBody m; Ixx, Iyy, Ixy centroidal waterplane
 *   second moments m⁴ (Ixx = ∫(y-yc)² dA governs heel, Iyy = ∫(x-xc)² dA governs trim); I1 >= I2 principal
 *   moments m⁴ and principalAngle rad (axis of I2 from earth x); KB m (height of B above the body origin level);
 *   BMt, BMl, KMt, KMl m; GMt, GMl m (null without cog); TPC t/cm; massPerDraft kg/m; MCT N·m per metre of
 *   trim over `length` and MCTcm t·m/cm (null without cog); stiffness number[6][6] in N/m, N, N·m
 */
export function meshHydrostatics({ mesh, draft, heel = 0, trim = 0, rho = 1025, g = G, cog, mass, length, lost }) {
  needPos(rho, 'rho');
  if (!Number.isFinite(draft)) throw new RangeError('draft must be finite');
  const r = bodyIntegrals(mesh, lost, draft, heel, trim);
  const V = r.V;
  if (!(V > 0)) throw new RangeError('draft gives no immersed volume');
  const cob = [r.Mx / V, r.My / V, r.Mz / V];
  const Aw = r.Aw;
  const hasWp = Aw > 1e-12 * Math.cbrt(V) ** 2;
  const xc = hasWp ? r.Sx / Aw : 0; const yc = hasWp ? r.Sy / Aw : 0;
  const Ixx = hasWp ? r.Jxx - Aw * yc * yc : 0;
  const Iyy = hasWp ? r.Jyy - Aw * xc * xc : 0;
  const Ixy = hasWp ? r.Jxy - Aw * xc * yc : 0;
  const avg = 0.5 * (Ixx + Iyy); const rad = Math.hypot(0.5 * (Ixx - Iyy), Ixy);
  const m = mass === undefined ? rho * V : mass;
  const cg = cog ? bodyToEarth(cog, draft, heel, trim) : null;
  const BMt = Ixx / V; const BMl = Iyy / V;
  const KB = cob[2] + draft;
  let L = length;
  if (L === undefined) {
    let lo = Infinity; let hi = -Infinity;
    for (let i = 0; i < mesh.positions.length; i += 3) { lo = Math.min(lo, mesh.positions[i]); hi = Math.max(hi, mesh.positions[i]); }
    L = hi - lo;
  }
  const GMt = cg ? BMt + cob[2] - cg[2] : null;
  const GMl = cg ? BMl + cob[2] - cg[2] : null;
  const C = [];
  for (let i = 0; i < 6; i++) C.push([0, 0, 0, 0, 0, 0]);
  const rg = rho * g;
  const zg = cg ? cg[2] : cob[2]; const xg = cg ? cg[0] : cob[0]; const yg = cg ? cg[1] : cob[1];
  C[2][2] = rg * Aw;
  C[2][3] = C[3][2] = rg * r.Sy;
  C[2][4] = C[4][2] = -rg * r.Sx;
  C[3][3] = rg * r.Jxx + rg * r.Mz - m * g * zg;
  C[4][4] = rg * r.Jyy + rg * r.Mz - m * g * zg;
  C[3][4] = C[4][3] = -rg * r.Jxy;
  C[3][5] = -rg * r.Mx + m * g * xg;
  C[4][5] = -rg * r.My + m * g * yg;
  return {
    volume: V, displacement: rho * V, buoyancy: rg * V, cob, cobBody: earthToBody(cob, draft, heel, trim),
    wettedArea: r.wet, waterplaneArea: Aw, cof: [xc, yc, 0], cofBody: earthToBody([xc, yc, 0], draft, heel, trim),
    Ixx, Iyy, Ixy, I1: avg + rad, I2: avg - rad, principalAngle: 0.5 * Math.atan2(-2 * Ixy, Ixx - Iyy),
    KB, BMt, BMl, KMt: KB + BMt, KMl: KB + BMl, GMt, GMl,
    TPC: (rho * Aw) / 1e5, massPerDraft: rho * Aw,
    MCT: GMl === null ? null : (rg * V * GMl) / L, MCTcm: GMl === null ? null : (rho * V * GMl) / L / 1e5,
    stiffness: C,
  };
}

function finishMesh(pos, idx, name) {
  return { kind: 'surface', name, positions: Float32Array.from(pos), indices: Uint32Array.from(idx) };
}

/**
 * Closed, outward-oriented triangulated box.
 * @param {{xMin:number,xMax:number,yMin:number,yMax:number,zMin:number,zMax:number}} o extents, m
 * @returns {{kind:'surface', name:string, positions:Float32Array, indices:Uint32Array}}
 */
export function boxMesh({ xMin, xMax, yMin, yMax, zMin, zMax }) {
  if (!(xMax > xMin && yMax > yMin && zMax > zMin)) throw new RangeError('box extents must satisfy min < max');
  const pos = [];
  for (const z of [zMin, zMax]) for (const y of [yMin, yMax]) for (const x of [xMin, xMax]) pos.push(x, y, z);
  // vertex index = ix + 2 iy + 4 iz
  const quads = [[0, 2, 3, 1], [4, 5, 7, 6], [0, 1, 5, 4], [2, 6, 7, 3], [0, 4, 6, 2], [1, 3, 7, 5]];
  const idx = [];
  for (const [a, b, c, d] of quads) idx.push(a, b, c, a, c, d);
  return finishMesh(pos, idx, 'box');
}

/**
 * Closed, outward-oriented triangulated circular cylinder with axis along z, centred at (x0, y0).
 * The polygonal section is scaled so that its area equals the exact circle area.
 * @param {{D:number, zMin:number, zMax:number, x0?:number, y0?:number, nTheta?:number}} o m; nTheta segments (>= 8)
 * @returns {{kind:'surface', name:string, positions:Float32Array, indices:Uint32Array}}
 */
export function cylinderMesh({ D, zMin, zMax, x0 = 0, y0 = 0, nTheta = 48 }) {
  needPos(D, 'D');
  if (!(zMax > zMin)) throw new RangeError('zMax must exceed zMin');
  if (!(nTheta >= 8)) throw new RangeError('nTheta must be >= 8');
  const n = Math.round(nTheta);
  const rp = (D / 2) * Math.sqrt((2 * Math.PI) / (n * Math.sin((2 * Math.PI) / n)));
  const pos = [];
  for (const z of [zMin, zMax]) for (let i = 0; i < n; i++) pos.push(x0 + rp * Math.cos((2 * Math.PI * i) / n), y0 + rp * Math.sin((2 * Math.PI * i) / n), z);
  pos.push(x0, y0, zMin, x0, y0, zMax);
  const idx = [];
  for (let i = 0; i < n; i++) {
    const j = (i + 1) % n;
    idx.push(i, j, n + j, i, n + j, n + i);
    idx.push(2 * n, j, i);
    idx.push(2 * n + 1, n + i, n + j);
  }
  return finishMesh(pos, idx, 'cylinder');
}

/**
 * Concatenates closed surfaces into one mesh (the union must not overlap for hydrostatics to be additive).
 * @param {object[]} meshes
 * @returns {{kind:'surface', name:string, positions:Float32Array, indices:Uint32Array}}
 */
export function mergeMeshes(meshes) {
  let np = 0; let ni = 0;
  for (const m of meshes) { np += m.positions.length; ni += m.indices.length; }
  const positions = new Float32Array(np);
  const indices = new Uint32Array(ni);
  let op = 0; let oi = 0;
  for (const m of meshes) {
    positions.set(m.positions, op);
    for (let i = 0; i < m.indices.length; i++) indices[oi + i] = m.indices[i] + op / 3;
    op += m.positions.length; oi += m.indices.length;
  }
  return { kind: 'surface', name: 'merged', positions, indices };
}

/**
 * Closed-form upright hydrostatics of a rectangular box barge L × B at draft T:
 * V = L B T, KB = T/2, BMt = B²/(12 T), BMl = L²/(12 T), Aw = L B.
 * Reference: Biran & López-Pulido (2014) ch. 2 (wall-sided box).
 * @param {{L:number,B:number,T:number,KG?:number,rho?:number,g?:number}} o m, m, m, m, kg/m³, m/s²
 * @returns {{volume:number, displacement:number, waterplaneArea:number, KB:number, BMt:number, BMl:number,
 *   KMt:number, KMl:number, GMt:number|null, GMl:number|null, TPC:number, MCT:number|null, C33:number,
 *   C44:number|null, C55:number|null}} SI; TPC t/cm; MCT N·m per metre of trim; C33 N/m; C44, C55 N·m/rad
 */
export function boxHydrostatics({ L, B, T, KG, rho = 1025, g = G }) {
  needPos(L, 'L'); needPos(B, 'B'); needPos(T, 'T');
  const V = L * B * T;
  const KB = T / 2; const BMt = (B * B) / (12 * T); const BMl = (L * L) / (12 * T);
  const has = KG !== undefined;
  const GMt = has ? KB + BMt - KG : null; const GMl = has ? KB + BMl - KG : null;
  return {
    volume: V, displacement: rho * V, waterplaneArea: L * B, KB, BMt, BMl, KMt: KB + BMt, KMl: KB + BMl, GMt, GMl,
    TPC: (rho * L * B) / 1e5, MCT: has ? (rho * g * V * GMl) / L : null, C33: rho * g * L * B,
    C44: has ? rho * g * V * GMt : null, C55: has ? rho * g * V * GMl : null,
  };
}

/**
 * Closed-form upright hydrostatics of a vertical circular cylinder (spar, buoy) of diameter D at draft T:
 * V = pi D² T/4, KB = T/2, BM = D²/(16 T).
 * Reference: first principles (I = pi D⁴/64).
 * @param {{D:number,T:number,KG?:number,rho?:number,g?:number}} o m, m, m, kg/m³, m/s²
 * @returns {{volume:number, displacement:number, waterplaneArea:number, KB:number, BM:number, KM:number,
 *   GM:number|null, C33:number, C44:number|null}} SI
 */
export function verticalCylinderHydrostatics({ D, T, KG, rho = 1025, g = G }) {
  needPos(D, 'D'); needPos(T, 'T');
  const Aw = (Math.PI * D * D) / 4;
  const V = Aw * T;
  const KB = T / 2; const BM = (D * D) / (16 * T);
  const GM = KG === undefined ? null : KB + BM - KG;
  return {
    volume: V, displacement: rho * V, waterplaneArea: Aw, KB, BM, KM: KB + BM, GM,
    C33: rho * g * Aw, C44: GM === null ? null : rho * g * V * GM,
  };
}

/**
 * Closed-form upright hydrostatics of a multi-column body (semi-submersible, TLP hull): vertical
 * surface-piercing columns (circular D, or rectangular a × b with a along x) plus fully submerged volumes
 * (pontoons, bracings). Waterplane inertia by the parallel-axis theorem about the centre of flotation.
 * Reference: first principles; e.g. Biran & López-Pulido (2014) ch. 2.
 * @param {object} o
 * @param {{x:number, y:number, D?:number, a?:number, b?:number, zBottom:number}[]} o.columns plan position m,
 *   section m, elevation of the column base relative to the keel m (the column is wet from zBottom to the draft)
 * @param {{volume:number, x?:number, y?:number, z:number}[]} [o.submerged] volumes m³ with centroids (z above keel) m
 * @param {number} o.T draft above keel, m
 * @param {number} [o.KG] m
 * @param {number} [o.rho=1025] kg/m³
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{volume:number, displacement:number, waterplaneArea:number, cof:number[], cob:number[], KB:number,
 *   Ixx:number, Iyy:number, BMt:number, BMl:number, KMt:number, KMl:number, GMt:number|null, GMl:number|null,
 *   C33:number, C44:number|null, C55:number|null}} SI (Ixx governs heel, Iyy trim)
 */
export function multiColumnHydrostatics({ columns, submerged = [], T, KG, rho = 1025, g = G }) {
  needPos(T, 'T');
  if (!Array.isArray(columns) || columns.length === 0) throw new RangeError('columns must be a non-empty array');
  let V = 0; let Mx = 0; let My = 0; let Mz = 0; let Aw = 0; let Sx = 0; let Sy = 0; let Jxx = 0; let Jyy = 0;
  for (const c of columns) {
    let A; let ix; let iy;
    if (c.D) { A = (Math.PI * c.D * c.D) / 4; ix = iy = (Math.PI * c.D ** 4) / 64; } else {
      needPos(c.a, 'column a'); needPos(c.b, 'column b');
      A = c.a * c.b; ix = (c.a * c.b ** 3) / 12; iy = (c.b * c.a ** 3) / 12;
    }
    const h = T - c.zBottom;
    if (!(h > 0)) throw new RangeError('draft T must exceed every column zBottom');
    V += A * h; Mx += A * h * c.x; My += A * h * c.y; Mz += A * h * (c.zBottom + h / 2);
    Aw += A; Sx += A * c.x; Sy += A * c.y; Jxx += ix + A * c.y * c.y; Jyy += iy + A * c.x * c.x;
  }
  for (const s of submerged) {
    V += s.volume; Mx += s.volume * (s.x || 0); My += s.volume * (s.y || 0); Mz += s.volume * s.z;
  }
  const xc = Sx / Aw; const yc = Sy / Aw;
  const Ixx = Jxx - Aw * yc * yc; const Iyy = Jyy - Aw * xc * xc;
  const KB = Mz / V; const BMt = Ixx / V; const BMl = Iyy / V;
  const has = KG !== undefined;
  const GMt = has ? KB + BMt - KG : null; const GMl = has ? KB + BMl - KG : null;
  return {
    volume: V, displacement: rho * V, waterplaneArea: Aw, cof: [xc, yc], cob: [Mx / V, My / V, KB], KB, Ixx, Iyy,
    BMt, BMl, KMt: KB + BMt, KMl: KB + BMl, GMt, GMl, C33: rho * g * Aw,
    C44: has ? rho * g * V * GMt : null, C55: has ? rho * g * V * GMl : null,
  };
}

function zExtent(mesh) {
  let lo = Infinity; let hi = -Infinity; let r = 0;
  const p = mesh.positions;
  for (let i = 0; i < p.length; i += 3) {
    lo = Math.min(lo, p[i + 2]); hi = Math.max(hi, p[i + 2]);
    r = Math.max(r, Math.abs(p[i]), Math.abs(p[i + 1]), Math.abs(p[i + 2]));
  }
  return { lo, hi, r };
}

/**
 * Draft at which the immersed volume equals a target, for fixed heel and trim (bisection on the monotonic
 * volume–draft curve).
 * @param {object} o
 * @param {object} o.mesh closed surface
 * @param {number} o.volume target immersed volume, m³
 * @param {number} [o.heel=0] rad
 * @param {number} [o.trim=0] rad
 * @param {object[]} [o.lost] damaged compartments
 * @returns {number} draft, m
 */
export function draftForVolume({ mesh, volume, heel = 0, trim = 0, lost }) {
  needPos(volume, 'volume');
  const e = zExtent(mesh);
  let lo = -2 * e.r - 1; let hi = 2 * e.r + 1;
  const vol = (T) => bodyIntegrals(mesh, lost, T, heel, trim).V;
  if (vol(hi) < volume) throw new RangeError('volume exceeds the total (intact) volume of the body: it sinks');
  for (let i = 0; i < 80; i++) {
    const mid = 0.5 * (lo + hi);
    if (vol(mid) < volume) lo = mid; else hi = mid;
    if (hi - lo < 1e-12 * (1 + e.r)) break;
  }
  return 0.5 * (lo + hi);
}

/**
 * Free-floating equilibrium (draft, heel, trim) of a body of given mass and centre of gravity: Newton
 * iteration on rho·V = m and on the horizontal coincidence of B and G in the earth frame, with a
 * finite-difference Jacobian and step limiting. Heel and/or trim can be held fixed.
 * Reference: Archimedes' principle and moment equilibrium; Biran & López-Pulido (2014) ch. 2.
 * Finds the equilibrium nearest the starting guess (upright by default); it does not search for capsized ones.
 * @param {object} o
 * @param {object} o.mesh closed surface, m
 * @param {number} o.mass kg
 * @param {number[]} o.cog centre of gravity, body coordinates, m
 * @param {number} [o.rho=1025] kg/m³
 * @param {{heel?:number, trim?:number, draft?:number}} [o.guess] starting angles rad and draft m (the draft
 *   defaults to the one matching the volume at the starting angles)
 * @param {boolean} [o.fixHeel=false]
 * @param {boolean} [o.fixTrim=false]
 * @param {object[]} [o.lost] damaged compartments (lost buoyancy)
 * @param {number} [o.tol=1e-9] convergence tolerance on the scaled residual
 * @param {number} [o.maxIter=60]
 * @returns {{draft:number, heel:number, trim:number, iterations:number, converged:boolean, residual:number}}
 *   m, rad, rad
 */
export function solveEquilibrium({ mesh, mass, cog, rho = 1025, guess = {}, fixHeel = false, fixTrim = false, lost, tol = 1e-9, maxIter = 60 }) {
  needPos(mass, 'mass'); needPos(rho, 'rho');
  if (!cog || cog.length !== 3) throw new RangeError('cog must be [x, y, z]');
  const V0 = mass / rho;
  const ext = zExtent(mesh);
  const Lc = ext.r || 1;
  let heel = guess.heel || 0; let trim = guess.trim || 0;
  let draft = Number.isFinite(guess.draft) ? guess.draft : draftForVolume({ mesh, volume: V0, heel, trim, lost });
  const free = [true, !fixHeel, !fixTrim];
  const res = (T, h, t) => {
    const r = bodyIntegrals(mesh, lost, T, h, t);
    const cg = bodyToEarth(cog, T, h, t);
    const out = [(r.V - V0) / V0];
    if (free[1]) out.push(r.V > 0 ? (r.My / r.V - cg[1]) / Lc : 1);
    if (free[2]) out.push(r.V > 0 ? (r.Mx / r.V - cg[0]) / Lc : 1);
    return out;
  };
  const norm = (v) => Math.sqrt(v.reduce((s, x) => s + x * x, 0));
  const pack = () => { const x = [draft]; if (free[1]) x.push(heel); if (free[2]) x.push(trim); return x; };
  const unpack = (x) => { let i = 1; return [x[0], free[1] ? x[i++] : heel, free[2] ? x[i] : trim]; };
  let it = 0; let rn = Infinity; let converged = false;
  for (; it < maxIter; it++) {
    const x = pack();
    const f = res(draft, heel, trim);
    rn = norm(f);
    if (rn < tol) { converged = true; break; }
    const n = x.length;
    const J = [];
    for (let i = 0; i < n; i++) J.push(new Array(n).fill(0));
    for (let j = 0; j < n; j++) {
      const h = j === 0 ? 1e-5 * Lc : 1e-6;
      const xp = x.slice(); xp[j] += h;
      const xm = x.slice(); xm[j] -= h;
      const fp = res(...unpack(xp)); const fm = res(...unpack(xm));
      for (let i = 0; i < n; i++) J[i][j] = (fp[i] - fm[i]) / (2 * h);
    }
    let dx;
    try { dx = solveDense(J, f.map((v) => -v)); } catch (err) { break; }
    if (!dx.every(Number.isFinite)) break;
    let scale = 1;
    scale = Math.min(scale, (0.25 * (ext.hi - ext.lo)) / Math.max(Math.abs(dx[0]), 1e-300));
    for (let j = 1; j < n; j++) scale = Math.min(scale, 0.2 / Math.max(Math.abs(dx[j]), 1e-300));
    let ok = false;
    for (let ls = 0; ls < 12; ls++) {
      const xt = x.map((v, j) => v + scale * dx[j]);
      const ft = res(...unpack(xt));
      if (norm(ft) < rn || scale * norm(dx) < 1e-14) { [draft, heel, trim] = unpack(xt); ok = true; break; }
      scale *= 0.5;
    }
    if (!ok) break;
  }
  if (!converged) { rn = norm(res(draft, heel, trim)); converged = rn < Math.max(tol, 1e-7); }
  return { draft, heel, trim, iterations: it, converged, residual: rn };
}

/**
 * Free-surface correction of slack tanks: virtual rise of G, GG' = sum(rho_t i_t) / (rho V), where i_t is the
 * second moment of the liquid surface about its own centroidal axis parallel to the heel axis.
 * Rectangular tanks: i = l b³/12 (b athwartships). Mesh tanks: the tank surface (body coordinates, upright)
 * is clipped at the liquid level `fillLevel` (z, m) and i taken from the cut plane.
 * Reference: Biran & López-Pulido (2014) (free-surface effect); IMO 2008 IS Code Part B ch. 3 (free-surface moment).
 * Small-angle (wall-sided tank) correction; tanks that are full or empty carry no correction.
 * @param {object} o
 * @param {({l:number, b:number, rho?:number} | {mesh:object, fillLevel:number, rho?:number})[]} o.tanks
 *   rectangular: free-surface length (x) and breadth (y), m; liquid density kg/m³ (default 1025)
 * @param {number} o.volume displaced volume of the vessel, m³
 * @param {number} [o.rho=1025] sea-water density, kg/m³
 * @returns {{GGt:number, GGl:number, fsmT:number, fsmL:number, perTank:{iT:number, iL:number, fsmT:number,
 *   fsmL:number}[]}} transverse and longitudinal rise of G m; free-surface moments rho_t·i in kg·m (t·m × 1000)
 */
export function freeSurfaceCorrection({ tanks, volume, rho = 1025 }) {
  needPos(volume, 'volume');
  if (!Array.isArray(tanks)) throw new RangeError('tanks must be an array');
  let fsmT = 0; let fsmL = 0;
  const perTank = tanks.map((tk) => {
    const rt = tk.rho === undefined ? 1025 : tk.rho;
    let iT; let iL;
    if (tk.mesh) {
      const r = clipIntegrals(tk.mesh, tk.fillLevel, 0, 0);
      if (r.Aw > 0) {
        iT = r.Jxx - (r.Sy * r.Sy) / r.Aw;
        iL = r.Jyy - (r.Sx * r.Sx) / r.Aw;
      } else { iT = 0; iL = 0; }
    } else {
      needPos(tk.l, 'tank l'); needPos(tk.b, 'tank b');
      iT = (tk.l * tk.b ** 3) / 12; iL = (tk.b * tk.l ** 3) / 12;
    }
    fsmT += rt * iT; fsmL += rt * iL;
    return { iT, iL, fsmT: rt * iT, fsmL: rt * iL };
  });
  return { GGt: fsmT / (rho * volume), GGl: fsmL / (rho * volume), fsmT, fsmL, perTank };
}
