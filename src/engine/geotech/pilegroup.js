// Pile groups under a rigid cap: three-dimensional stiffness method and group-effect helpers.
//
// Global axes: x, y horizontal, z UP; loads and displacements are given at the reference point (0, 0, 0) of the
// cap. A downward load is therefore a negative Fz. Moments and rotations follow the right-hand rule.
// Pile local axes: e3 along the pile pointing down, e1 horizontal-most direction in the plane of the batter,
// e2 = e3 × e1. A vertical pile has e1 = +x, e2 = −y, e3 = −z.
// Pile-head stiffness terms (KHH, KHM, KMM) use the convention of pile.js: KHM is negative.

import { solveDense } from '../math/linalg.js';

const RAD = Math.PI / 180;

function pos(name, v, allowZero = false) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < 0 || (!allowZero && v === 0)) {
    throw new RangeError(`${name}: must be a finite number ${allowZero ? '>= 0' : '> 0'} (got ${v})`);
  }
}

function zeros6() {
  return Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]);
}

/**
 * 6×6 pile-head stiffness matrix in pile local axes, degrees of freedom [u1, u2, u3, r1, r2, r3]
 * (translations along and rotations about e1, e2, e3). Axisymmetric pile: the same KHH, KHM, KMM in both
 * lateral directions. With theta = −dy/dz of pile.js: theta(direction 1) = −r2 and theta(direction 2) = +r1,
 * hence K[u1][r2] = −KHM and K[u2][r1] = +KHM.
 * A pinned head carries no moment: the lateral stiffness is condensed to KHH − KHM²/KMM.
 * Reference: stiffness method for pile groups, e.g. Poulos & Davis (1980) ch. 9; first-principles assembly.
 * Validity: linear (use secant stiffnesses at the working load for nonlinear soil and iterate if necessary).
 * @param {{ KHH: number, KHM?: number, KMM?: number, Kaxial: number, Ktorsion?: number, head?: string }} o
 *   KHH [N/m]; KHM [N] (<= 0); KMM [N·m/rad]; Kaxial [N/m]; Ktorsion [N·m/rad] (default 0); head 'fixed'|'pinned'
 * @returns {number[][]} local stiffness matrix (SI units by block: N/m, N, N·m/rad)
 */
export function pileLocalStiffness({ KHH, KHM = 0, KMM = 0, Kaxial, Ktorsion = 0, head = 'fixed' }) {
  pos('KHH', KHH, true); pos('KMM', KMM, true); pos('Kaxial', Kaxial); pos('Ktorsion', Ktorsion, true);
  if (head !== 'fixed' && head !== 'pinned') throw new RangeError("head: must be 'fixed' or 'pinned'");
  if (KHM > 0) throw new RangeError('KHM: must be <= 0 in the pile.js sign convention');
  const K = zeros6();
  K[2][2] = Kaxial; K[5][5] = Ktorsion;
  if (head === 'pinned') {
    const kl = KMM > 0 ? KHH - (KHM * KHM) / KMM : KHH;
    K[0][0] = kl; K[1][1] = kl;
    return K;
  }
  K[0][0] = KHH; K[1][1] = KHH; K[3][3] = KMM; K[4][4] = KMM;
  K[0][4] = -KHM; K[4][0] = -KHM; K[1][3] = KHM; K[3][1] = KHM;
  return K;
}

function pileAxes(batterDeg, azimuthDeg) {
  const b = batterDeg * RAD, a = azimuthDeg * RAD;
  const e3 = [Math.sin(b) * Math.cos(a), Math.sin(b) * Math.sin(a), -Math.cos(b)];
  const e1 = [Math.cos(b) * Math.cos(a), Math.cos(b) * Math.sin(a), Math.sin(b)];
  const e2 = [e3[1] * e1[2] - e3[2] * e1[1], e3[2] * e1[0] - e3[0] * e1[2], e3[0] * e1[1] - e3[1] * e1[0]];
  return [e1, e2, e3];
}

// B = Q·T: maps cap displacement (6) to pile-head displacement in local axes (6).
function kinematics(p) {
  const ax = pileAxes(p.batterDeg ?? 0, p.azimuthDeg ?? 0);
  const r = [p.x, p.y, p.z ?? 0];
  const B = zeros6();
  for (let i = 0; i < 3; i++) {
    const e = ax[i];
    // translation: e·(u + θ × r) = e·u + θ·(r × e)
    const rxe = [r[1] * e[2] - r[2] * e[1], r[2] * e[0] - r[0] * e[2], r[0] * e[1] - r[1] * e[0]];
    for (let j = 0; j < 3; j++) { B[i][j] = e[j]; B[i][3 + j] = rxe[j]; B[3 + i][3 + j] = e[j]; }
  }
  return B;
}

/**
 * Rigid-cap pile group in three dimensions by the stiffness method: K = Σ Bᵀ K_local B, cap displacement
 * d = K⁻¹ F, pile-head forces f = K_local B d. Raked piles are defined by batter and azimuth.
 * Reference: stiffness (displacement) method for pile groups — Poulos & Davis (1980) ch. 9; Reese & Van Impe
 * (2011) ch. 5. Validity: rigid cap, linear pile-head stiffnesses, pile–soil–pile interaction only through the
 * stiffnesses supplied (apply group reduction factors beforehand).
 * @param {{ piles: object[], load: object }} o
 *   piles: [{ x, y, z? [m] head position; batterDeg? [deg] rake from vertical (default 0); azimuthDeg? [deg]
 *   direction of the rake measured from +x towards +y; stiffness: options of pileLocalStiffness, or
 *   K: number[][] a ready 6×6 local matrix; id? }];
 *   load: { Fx, Fy, Fz [N], Mx, My, Mz [N·m] } at the cap reference point (z up; missing terms are 0)
 * @returns {{ displacement: object, K: number[][], piles: object[], equilibriumError: number }}
 *   displacement { ux, uy, uz [m], rx, ry, rz [rad] }; K 6×6 group stiffness at the reference point;
 *   piles: [{ id, axial [N] (compression +), shear1, shear2, shear [N], moment1, moment2, moment [N·m] (about
 *   e1, e2 and resultant), torsion [N·m], local [6] head displacement in local axes }];
 *   equilibriumError [-] = |Σ pile forces − load| / |load|
 */
export function analysePileGroup({ piles, load }) {
  if (!Array.isArray(piles) || piles.length === 0) throw new RangeError('piles: must be a non-empty array');
  const F = ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'].map((k) => {
    const v = load?.[k] ?? 0;
    if (typeof v !== 'number' || !Number.isFinite(v)) throw new RangeError(`load.${k}: must be a finite number`);
    return v;
  });
  const K = zeros6();
  const data = piles.map((p, idx) => {
    if (![p.x, p.y].every(Number.isFinite)) throw new RangeError(`piles[${idx}].x: pile head coordinates must be finite`);
    const Kl = p.K ?? pileLocalStiffness(p.stiffness ?? {});
    const B = kinematics(p);
    for (let a = 0; a < 6; a++) for (let b = 0; b < 6; b++) {
      let s = 0;
      for (let i = 0; i < 6; i++) { const Bia = B[i][a]; if (Bia === 0) continue; for (let j = 0; j < 6; j++) s += Bia * Kl[i][j] * B[j][b]; }
      K[a][b] += s;
    }
    return { Kl, B };
  });
  let d;
  try { d = Array.from(solveDense(K, F)); } catch (err) { throw new Error('Pile group is a mechanism for the applied loading (singular group stiffness): add lateral, torsional or raked-pile stiffness.'); }
  if (!d.every(Number.isFinite)) throw new Error('Pile group is a mechanism for the applied loading (singular group stiffness).');
  const sum = [0, 0, 0, 0, 0, 0];
  const out = data.map(({ Kl, B }, idx) => {
    const ul = B.map((row) => row.reduce((s, v, j) => s + v * d[j], 0));
    const f = Kl.map((row) => row.reduce((s, v, j) => s + v * ul[j], 0));
    for (let a = 0; a < 6; a++) for (let i = 0; i < 6; i++) sum[a] += B[i][a] * f[i];
    return {
      id: piles[idx].id ?? idx, axial: f[2], shear1: f[0], shear2: f[1], shear: Math.hypot(f[0], f[1]),
      moment1: f[3], moment2: f[4], moment: Math.hypot(f[3], f[4]), torsion: f[5], local: ul,
    };
  });
  const nF = Math.hypot(...F) || 1;
  return {
    displacement: { ux: d[0], uy: d[1], uz: d[2], rx: d[3], ry: d[4], rz: d[5] }, K, piles: out,
    equilibriumError: Math.hypot(...sum.map((v, i) => v - F[i])) / nF,
  };
}

/**
 * Example p-multiplier table by row position and centre-to-centre spacing in the direction of loading:
 * spacing 3 D → 0.80 / 0.40 / 0.30 and spacing 5 D → 1.00 / 0.85 / 0.70 for row 1 (leading) / row 2 / row 3 and
 * higher. Supply a project table to pMultiplier where another source governs.
 * Reference: AASHTO LRFD Bridge Design Specifications (2012 and later), Table 10.7.2.4-1.
 */
export const P_MULTIPLIER_TABLE_AASHTO = Object.freeze({ spacingOverD: [3, 5], rows: [[0.8, 0.4, 0.3], [1.0, 0.85, 0.7]] });

/**
 * p-multiplier for a pile row from a spacing/row table, linearly interpolated in spacing and held constant
 * outside the tabulated spacings (above the largest tabulated spacing the last row of the table applies; set
 * the multiplier to 1 yourself where the source allows it).
 * Reference: the table supplied (default P_MULTIPLIER_TABLE_AASHTO); concept of Brown, Morrison & Reese (1988).
 * Validity: as stated by the source of the table.
 * @param {{ spacingOverD: number, row: number, table?: { spacingOverD: number[], rows: number[][] } }} o
 *   spacingOverD [-] centre-to-centre spacing / diameter in the load direction; row [-] 1 = leading row;
 *   table: spacings ascending and one array of multipliers (by row, last value used for further rows) per spacing
 * @returns {number} p-multiplier [-]
 */
export function pMultiplier({ spacingOverD, row, table = P_MULTIPLIER_TABLE_AASHTO }) {
  pos('spacingOverD', spacingOverD);
  if (!(Number.isInteger(row) && row >= 1)) throw new RangeError('row: must be an integer >= 1');
  const s = table.spacingOverD;
  const at = (i) => { const r = table.rows[i]; return r[Math.min(row, r.length) - 1]; };
  if (spacingOverD <= s[0]) return at(0);
  for (let i = 1; i < s.length; i++) {
    if (spacingOverD <= s[i]) return at(i - 1) + ((at(i) - at(i - 1)) * (spacingOverD - s[i - 1])) / (s[i] - s[i - 1]);
  }
  return at(s.length - 1);
}

/**
 * Axial pile-to-pile interaction factor from the logarithmic shaft displacement field of Randolph & Wroth:
 * α = ln(rm/s) / ln(rm/r0) for s < rm, else 0, with rm = 2.5 ρ (1 − ν) L.
 * Base interaction and pile compressibility are neglected (rigid, shaft-dominated piles).
 * Reference: Randolph & Wroth (1979) Géotechnique 29(4) 423–439. Validity: floating piles, linear elastic soil.
 * @param {{ s: number, D: number, L: number, nu: number, rho?: number }} o s [m] centre-to-centre spacing;
 *   D [m]; L [m]; nu [-]; rho [-] G_avg/G_L (default 1)
 * @returns {number} interaction factor α [-]
 */
export function axialInteractionFactor({ s, D, L, nu, rho = 1 }) {
  pos('s', s); pos('D', D); pos('L', L); pos('rho', rho);
  const r0 = D / 2, rm = 2.5 * rho * (1 - nu) * L;
  if (!(rm > r0)) throw new RangeError('L: too short for the Randolph–Wroth displacement field (rm <= r0)');
  if (s < D) throw new RangeError('s: spacing must be at least one diameter');
  return s >= rm ? 0 : Math.log(rm / s) / Math.log(rm / r0);
}

/**
 * Lateral interaction factor for deflection between two flexible piles, Randolph (1981):
 * α = c ρc (Ep/Gc)^(1/7) (r0/s) (1 + cos²ψ); where this exceeds 1/3 it is replaced by 1 − 2/√(27 α).
 * c = 0.6 for fixed-head piles. For free-head piles under force loading Randolph indicates a lower coefficient
 * (about 0.4); pass it explicitly through coefficient.
 * Reference: Randolph (1981) Géotechnique 31(2) 247–259. Validity: flexible piles (L >= lc), elastic soil.
 * @param {{ s: number, D: number, EI: number, Gc: number, rhoc?: number, psiDeg?: number, coefficient?: number }} o
 *   s [m]; D [m]; EI [N·m²]; Gc [Pa] characteristic modulus (see randolphLateral in pile.js); rhoc [-]
 *   (default 1); psiDeg [deg] angle between the line joining the piles and the load direction (default 0);
 *   coefficient [-] (default 0.6)
 * @returns {number} interaction factor α [-]
 */
export function lateralInteractionFactor({ s, D, EI, Gc, rhoc = 1, psiDeg = 0, coefficient = 0.6 }) {
  pos('s', s); pos('D', D); pos('EI', EI); pos('Gc', Gc); pos('rhoc', rhoc); pos('coefficient', coefficient);
  const r0 = D / 2, Ep = EI / ((Math.PI * r0 ** 4) / 4);
  const c = Math.cos(psiDeg * RAD);
  const a = coefficient * rhoc * Math.pow(Ep / Gc, 1 / 7) * (r0 / s) * (1 + c * c);
  return a > 1 / 3 ? 1 - 2 / Math.sqrt(27 * a) : a;
}

/**
 * Group stiffness under a rigid cap from pairwise interaction factors: displacement of pile i is
 * w_i = (1/K1) Σ_j α_ij P_j with α_ii = 1; for a common displacement w the pile loads follow from the linear
 * system and K_group = Σ P / w. Works for axial (alphaFn from axialInteractionFactor) or lateral translation.
 * Reference: interaction-factor superposition, Poulos (1968) Géotechnique 18(4); Poulos & Davis (1980) ch. 6.
 * Validity: identical piles, linear response.
 * @param {{ piles: { x: number, y: number }[], K1: number, alphaFn: Function }} o piles [m] plan positions;
 *   K1 [N/m] single-pile stiffness; alphaFn(s, psiDeg) → α for spacing s [m] and bearing psiDeg [deg] of the
 *   neighbour measured from +x
 * @returns {{ Kgroup: number, efficiency: number, loadShare: number[], settlementRatio: number }} Kgroup [N/m];
 *   efficiency = Kgroup/(n K1) [-]; loadShare = P_i / mean P [-]; settlementRatio = 1/efficiency [-]
 */
export function groupStiffnessFromInteraction({ piles, K1, alphaFn }) {
  pos('K1', K1);
  const n = piles.length;
  if (n < 1) throw new RangeError('piles: must not be empty');
  const A = Array.from({ length: n }, (_, i) => piles.map((pj, j) => {
    if (i === j) return 1;
    const dx = pj.x - piles[i].x, dy = pj.y - piles[i].y;
    return alphaFn(Math.hypot(dx, dy), Math.atan2(dy, dx) / RAD);
  }));
  const P = Array.from(solveDense(A, new Array(n).fill(K1))); // loads for unit cap displacement
  const total = P.reduce((s, v) => s + v, 0);
  return { Kgroup: total, efficiency: total / (n * K1), loadShare: P.map((v) => (v * n) / total), settlementRatio: (n * K1) / total };
}

/**
 * Converse–Labarre axial group efficiency: Eg = 1 − θ [(n − 1) m + (m − 1) n] / (90 m n), θ = atan(D/s) in
 * degrees, for m rows of n piles.
 * Reference: Converse–Labarre formula as given in Bowles, Foundation Analysis and Design (5th ed.), pile-group efficiency.
 * Validity: friction piles in clay in a rectangular pattern; empirical.
 * @param {{ m: number, n: number, D: number, s: number }} o m [-] rows; n [-] piles per row; D [m]; s [m] spacing
 * @returns {number} efficiency Eg [-]
 */
export function converseLabarre({ m, n, D, s }) {
  if (!(Number.isInteger(m) && m >= 1)) throw new RangeError('m: must be an integer >= 1');
  if (!(Number.isInteger(n) && n >= 1)) throw new RangeError('n: must be an integer >= 1');
  pos('D', D); pos('s', s);
  const theta = Math.atan(D / s) / RAD;
  return 1 - (theta * ((n - 1) * m + (m - 1) * n)) / (90 * m * n);
}

/**
 * Ultimate torsional resistance of a pile group about a vertical axis: each pile contributes its ultimate
 * lateral resistance times its lever arm plus its own ultimate torsional resistance, T = Σ (Hu_i r_i + Tu_i).
 * (The elastic torsional stiffness of the group is K[5][5] of analysePileGroup.)
 * Reference: statics (upper estimate: all piles reach Hu simultaneously, perpendicular to their radius).
 * Validity: rigid cap; lateral–torsional interaction in each pile neglected.
 * @param {{ piles: { x: number, y: number, Hu: number, Tu?: number }[], centre?: { x: number, y: number } }} o
 *   piles: plan position [m], Hu [N] ultimate lateral resistance, Tu [N·m] ultimate torsion of the pile itself;
 *   centre [m] axis of rotation (default: centroid of the pile heads)
 * @returns {{ Tult: number, centre: { x: number, y: number } }} Tult [N·m]
 */
export function groupTorsionalCapacity({ piles, centre }) {
  if (!Array.isArray(piles) || piles.length === 0) throw new RangeError('piles: must be a non-empty array');
  const c = centre ?? { x: piles.reduce((s, p) => s + p.x, 0) / piles.length, y: piles.reduce((s, p) => s + p.y, 0) / piles.length };
  let T = 0;
  for (const p of piles) { pos('Hu', p.Hu, true); T += p.Hu * Math.hypot(p.x - c.x, p.y - c.y) + (p.Tu ?? 0); }
  return { Tult: T, centre: c };
}

/**
 * Utilisation table for the piles of a group: axial (compression or tension), lateral and bending demand
 * against the capacities supplied.
 * Reference: demand / capacity ratios; the capacities and their safety format are the caller's.
 * @param {{ piles: object[], capacity: { compression: number, tension: number, lateral?: number,
 *   moment?: number } }} o piles: the piles array returned by analysePileGroup; capacity: design resistances
 *   [N], [N], [N], [N·m]
 * @returns {{ rows: object[], maxUtilisation: number, governing: string|number }} rows: [{ id, axial [N], shear
 *   [N], moment [N·m], ucAxial, ucLateral, ucMoment, uc [-] }]
 */
export function pileUtilisationTable({ piles, capacity }) {
  pos('capacity.compression', capacity?.compression); pos('capacity.tension', capacity?.tension);
  let max = 0, gov = null;
  const rows = piles.map((p) => {
    const ucAxial = p.axial >= 0 ? p.axial / capacity.compression : -p.axial / capacity.tension;
    const ucLateral = capacity.lateral ? p.shear / capacity.lateral : 0;
    const ucMoment = capacity.moment ? p.moment / capacity.moment : 0;
    const uc = Math.max(ucAxial, ucLateral, ucMoment);
    if (uc > max) { max = uc; gov = p.id; }
    return { id: p.id, axial: p.axial, shear: p.shear, moment: p.moment, ucAxial, ucLateral, ucMoment, uc };
  });
  return { rows, maxUtilisation: max, governing: gov };
}
