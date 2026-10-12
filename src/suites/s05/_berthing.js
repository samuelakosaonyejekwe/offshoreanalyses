// Berthing energy by the kinetic-energy method and fender performance curves.
// Practice followed: PIANC WG33 "Guidelines for the Design of Fender Systems" (2002) and PIANC WG211 (2024),
// BS 6349-4:2014, ROM 2.0-11. Every formula is derived or stated with its origin below. SI units throughout.
//
// Coefficients and tables in this file that originate from published guidance are INDICATIVE DEFAULTS exposed to
// the user through module inputs; they must be confirmed against the PIANC edition in force on the project.
import { erfc } from '../../engine/math/special.js';

/* ------------------------------------------------------------------------------------------ berthing energy */

/**
 * Block coefficient Cb = M / (ρ L B D): displaced volume over the circumscribing box of the underwater hull.
 * M displacement [kg], L length between perpendiculars [m], B beam [m], D draft [m], rho [kg/m³].
 */
export function blockCoefficient({ M, L, B, D, rho = 1025 }) {
  return M / (rho * L * B * D);
}

/**
 * Radius of gyration of the vessel about the vertical axis through its centre of mass, K = (0.19 Cb + 0.11) L.
 * Empirical fit quoted in PIANC WG33 (2002) §4.2.4 and BS 6349-4; about 0.25 L for a full-form tanker (Cb ≈ 0.8)
 * and 0.2 L for a fine-form vessel. Returned in metres.
 */
export function radiusOfGyration({ Cb, L }) {
  return (0.19 * Cb + 0.11) * L;
}

/**
 * Eccentricity factor from rigid-body mechanics.
 *
 * A vessel of mass M and yaw radius of gyration K moves with velocity v. It touches the fender at a point P at
 * distance R from its centre of mass; γ is the angle between the velocity vector and the line from P to the
 * centre of mass. Assume the contact point is brought to rest without sliding (the fender arrests the hull at P)
 * and the vessel is otherwise free. The fender force passes through P, so the angular momentum about P is
 * conserved during the impact:
 *     before:  M v · (R sin γ)                 (moment arm of the linear momentum about P)
 *     after:   M (K² + R²) ω                   (rotation about P, parallel-axis theorem)
 *     ⇒  ω = v R sin γ / (K² + R²).
 * The kinetic energy left in the vessel is ½ M (K² + R²) ω² = ½ M v² · R² sin²γ / (K² + R²), so the fraction that
 * must be absorbed by fender and structure is
 *     Ce = 1 − R² sin²γ / (K² + R²) = (K² + R² cos²γ) / (K² + R²).
 * Limits: R = 0 (midships contact) → Ce = 1; γ = 90° → Ce = K² / (K² + R²).
 * K, R in metres, gamma in radians.
 */
export function eccentricityFactor({ K, R, gamma }) {
  const c = Math.cos(gamma);
  return (K * K + R * R * c * c) / (K * K + R * R);
}

/**
 * Geometry of the contact point for a side berthing: the contact is at distance x from the bow along the hull
 * (x = L/4 quarter-point, L/3 third-point, L/2 midships), on the hull side at half the beam from the centreline:
 *     R = √[(L/2 − x)² + (B/2)²]
 *     γ = 90° − α − asin(B / (2R))
 * where α is the berthing angle between the hull and the berthing line and the velocity is normal to the
 * berthing line (PIANC WG33 §4.2.4). Returns R [m] and γ [rad].
 */
export function contactGeometry({ L, B, x, alpha }) {
  const R = Math.hypot(L / 2 - x, B / 2);
  const gamma = Math.PI / 2 - alpha - Math.asin(Math.min(B / (2 * R), 1));
  return { R, gamma };
}

/**
 * Added (hydrodynamic) mass factor Cm = (M + M_added)/M for sway.
 *  'vasco-costa': Vasco Costa (1964), Cm = 1 + 2D/B — an empirical expression in which the entrained water grows
 *     with the draft-to-beam ratio. Stated valid for under-keel clearance ≥ 0.1 D and velocity ≥ 0.08 m/s.
 *  'ueda': Ueda (1981) from model tests and field observation, Cm = 1 + π D / (2 Cb B): the added mass is that of a
 *     half cylinder of water of radius D per unit length, ρ π D² L / 2, divided by the vessel mass ρ Cb L B D.
 *  'pianc': PIANC WG33 (2002) §4.2.3 values by under-keel clearance Kc: Cm = 1.8 for Kc/D ≤ 0.1,
 *     1.875 − 0.75 Kc/D for 0.1 < Kc/D < 0.5, 1.5 for Kc/D ≥ 0.5.
 *  'user': value supplied.
 */
export function addedMassFactor({ method, D, B, Cb, Kc, user }) {
  if (method === 'vasco-costa') return 1 + 2 * D / B;
  if (method === 'ueda') return 1 + Math.PI * D / (2 * Cb * B);
  if (method === 'pianc') { const r = Kc / D; return r <= 0.1 ? 1.8 : r >= 0.5 ? 1.5 : 1.875 - 0.75 * r; }
  return user;
}

/**
 * Berthing energy by the kinetic-energy method:
 *     E_N = ½ M v² · Ce · Cm · Cs · Cc           (normal berthing energy, J)
 *     E_A = F_ab · E_N                           (abnormal berthing energy)
 * ½ M v² is the translational kinetic energy of the vessel normal to the berth; Cm adds the entrained water;
 * Ce removes the energy left in yaw rotation (derived above); Cs (softness) removes the part absorbed by elastic
 * deformation of the hull (1.0 for soft fenders deflecting more than about 150 mm, 0.9 for hard fenders); Cc
 * (berth configuration) removes the cushioning of water trapped against a closed quay (1.0 for open piled
 * structures such as dolphins). F_ab is the abnormal-berthing factor of PIANC WG33 Table 4.2.5.
 */
export function berthingEnergy({ M, v, Ce, Cm, Cs = 1, Cc = 1, Fab = 1 }) {
  const EN = 0.5 * M * v * v * Ce * Cm * Cs * Cc;
  return { EN, EA: Fab * EN };
}

/**
 * Indicative berthing-velocity guidance (Brolsma-type curves as tabulated in fender design manuals after PIANC
 * WG33 Table 4.2.1): velocity normal to the berth [m/s] against deadweight [t] for five navigation conditions
 *   a easy berthing, sheltered;  b difficult berthing, sheltered;  c easy berthing, exposed;
 *   d good berthing, exposed;    e difficult berthing, exposed.
 * Indicative only: WG211 (2024) revises the velocity basis using measured berthings — confirm with the edition in
 * force and with pilot / port data.
 */
export const VELOCITY_TABLE_DEFAULT = [
  { dwt: 1000, a: 0.179, b: 0.343, c: 0.517, d: 0.669, e: 0.865 },
  { dwt: 2000, a: 0.151, b: 0.296, c: 0.445, d: 0.577, e: 0.726 },
  { dwt: 5000, a: 0.117, b: 0.236, c: 0.352, d: 0.459, e: 0.558 },
  { dwt: 10000, a: 0.094, b: 0.192, c: 0.287, d: 0.377, e: 0.448 },
  { dwt: 20000, a: 0.074, b: 0.153, c: 0.228, d: 0.303, e: 0.355 },
  { dwt: 50000, a: 0.052, b: 0.110, c: 0.164, d: 0.221, e: 0.258 },
  { dwt: 100000, a: 0.039, b: 0.083, c: 0.126, d: 0.171, e: 0.201 },
  { dwt: 200000, a: 0.028, b: 0.062, c: 0.095, d: 0.131, e: 0.158 },
  { dwt: 500000, a: 0.017, b: 0.041, c: 0.064, d: 0.090, e: 0.115 },
];

/** Velocity from a guidance table by interpolation in log(deadweight); clamped at the ends. Returns m/s. */
export function velocityFromTable(table, dwt, cond) {
  const rows = (table || []).filter((r) => r.dwt > 0 && r[cond] > 0).sort((p, q) => p.dwt - q.dwt);
  if (rows.length < 1) throw new Error(`The berthing-velocity table has no usable row for condition "${cond}".`);
  if (dwt <= rows[0].dwt) return rows[0][cond];
  const last = rows[rows.length - 1];
  if (dwt >= last.dwt) return last[cond];
  for (let i = 1; i < rows.length; i++) {
    if (dwt <= rows[i].dwt) {
      const f = Math.log(dwt / rows[i - 1].dwt) / Math.log(rows[i].dwt / rows[i - 1].dwt);
      return rows[i - 1][cond] + f * (rows[i][cond] - rows[i - 1][cond]);
    }
  }
  return last[cond];
}

/** Standard normal distribution function. */
export const normCdf = (x) => 0.5 * erfc(-x / Math.SQRT2);

/**
 * Probability that a lognormally distributed berthing velocity exceeds vCap, given the characteristic velocity vk
 * taken as the fractile pk of the distribution with coefficient of variation cov:
 *   σ_ln = √ln(1 + cov²),  μ_ln = ln vk − z_k σ_ln,  P(v > vCap) = 1 − Φ[(ln vCap − μ_ln)/σ_ln].
 * Energy scales with v², so P(E > E_cap) = P(v > v_cap) with v_cap = vk √(E_cap/E_k).
 */
export function exceedanceLognormal({ vk, zk, cov, vCap }) {
  const s = Math.sqrt(Math.log(1 + cov * cov));
  const mu = Math.log(vk) - zk * s;
  return 1 - normCdf((Math.log(vCap) - mu) / s);
}

/* ----------------------------------------------------------------------------------------- fender curves */

/**
 * GENERIC normalised performance curves: deflection as % of fender height (or diameter) against reaction as a
 * fraction of the rated reaction. The shapes reproduce the character of each family — buckling fenders (cone,
 * cell, arch) reach the rated reaction early, hold a plateau with a shallow dip and stiffen beyond the rated
 * deflection; cylindrical and pneumatic fenders stiffen progressively. They are not any manufacturer's data: a
 * certified curve entered in the module's table replaces them.
 */
export const GENERIC_FENDERS = Object.freeze({
  cone: { label: 'Cone (buckling)', rated: 72, pts: [[0, 0], [5, 0.19], [10, 0.39], [15, 0.59], [20, 0.75], [25, 0.89], [30, 0.97], [35, 1.0], [40, 0.98], [45, 0.92], [50, 0.87], [55, 0.84], [60, 0.84], [65, 0.89], [70, 0.97], [72, 1.0], [75, 1.1]] },
  cell: { label: 'Cell (buckling)', rated: 52.5, pts: [[0, 0], [5, 0.25], [10, 0.5], [15, 0.72], [20, 0.89], [25, 0.98], [27.5, 1.0], [30, 0.99], [35, 0.95], [40, 0.92], [45, 0.92], [50, 0.97], [52.5, 1.0], [55, 1.08]] },
  arch: { label: 'Arch / V (buckling)', rated: 52.5, pts: [[0, 0], [5, 0.22], [10, 0.45], [15, 0.66], [20, 0.83], [25, 0.94], [30, 1.0], [35, 0.98], [40, 0.94], [45, 0.93], [50, 0.97], [52.5, 1.0], [55, 1.07]] },
  cylindrical: { label: 'Cylindrical (side-loaded)', rated: 50, pts: [[0, 0], [10, 0.08], [20, 0.2], [30, 0.38], [40, 0.63], [50, 1.0], [55, 1.3]] },
  pneumatic: { label: 'Pneumatic (floating)', rated: 60, pts: [[0, 0], [10, 0.05], [20, 0.14], [30, 0.28], [40, 0.46], [50, 0.7], [60, 1.0], [65, 1.2]] },
});

/**
 * Dimensional fender curve for one unit.
 *   type: key of GENERIC_FENDERS; H fender height / diameter [m]; RR rated reaction [N];
 *   table: optional certified curve rows { defl [% of H], reaction [kN] } replacing the generic shape; the rated
 *     point is at ratedDefl [%], or at the last row when ratedDefl is not given.
 * Returns { d [m], R [N], E [J] (E(δ) = ∫₀^δ R dδ by the trapezoidal rule on a refined grid — the energy is the area
 * under the force–deflection curve), dRated, RRated, ERated, generic }.
 */
export function fenderCurve({ type, H, RR, table, ratedDefl }) {
  let pts, rated, generic = true;
  const rows = (table || []).filter((r) => Number.isFinite(r.defl) && Number.isFinite(r.reaction));
  if (rows.length >= 2) {
    const s = rows.slice().sort((p, q) => p.defl - q.defl);
    if (s[0].defl !== 0 || s[0].reaction !== 0) s.unshift({ defl: 0, reaction: 0 });
    for (let i = 1; i < s.length; i++) if (!(s[i].defl > s[i - 1].defl) || s[i].reaction < 0) throw new Error('The manufacturer fender curve must have strictly increasing deflection and non-negative reaction.');
    pts = s.map((r) => [r.defl, r.reaction * 1e3]);
    rated = ratedDefl > 0 ? Math.min(ratedDefl, s[s.length - 1].defl) : s[s.length - 1].defl;
    generic = false;
  } else {
    const g = GENERIC_FENDERS[type];
    if (!g) throw new Error(`Unknown fender type "${type}".`);
    pts = g.pts.map(([x, r]) => [x, r * RR]); rated = g.rated;
  }
  // refine to 0.5 % steps by linear interpolation so that the energy integral and inverse look-ups are smooth
  const dMaxPct = pts[pts.length - 1][0];
  const n = Math.max(2, Math.round(dMaxPct / 0.5)) + 1;
  const d = new Array(n), R = new Array(n), E = new Array(n);
  let j = 0;
  for (let i = 0; i < n; i++) {
    const x = dMaxPct * i / (n - 1);
    while (j < pts.length - 2 && x > pts[j + 1][0]) j++;
    const f = (x - pts[j][0]) / (pts[j + 1][0] - pts[j][0]);
    d[i] = x / 100 * H; R[i] = pts[j][1] + f * (pts[j + 1][1] - pts[j][1]);
    E[i] = i === 0 ? 0 : E[i - 1] + 0.5 * (R[i] + R[i - 1]) * (d[i] - d[i - 1]);
  }
  const at = (x) => interp(d, R, x), eAt = (x) => interp(d, E, x);
  const dRated = rated / 100 * H;
  return { d, R, E, dRated, RRated: at(dRated), ERated: eAt(dRated), dMax: d[n - 1], generic, H };
}

export function interp(xs, ys, x) {
  const n = xs.length;
  if (x <= xs[0]) return ys[0];
  if (x >= xs[n - 1]) return ys[n - 1];
  let lo = 0, hi = n - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
  return ys[lo] + (ys[hi] - ys[lo]) * (x - xs[lo]) / (xs[hi] - xs[lo]);
}

/**
 * State of a fender (or of nUnits identical units sharing the impact equally) that has absorbed energy E [J]:
 * deflection where E(δ)·fE = E, the reaction there and the largest reaction passed on the way (a buckling fender
 * peaks before the rated deflection). fE and fR are the combined correction factors on energy and reaction
 * (velocity, temperature, angle, tolerance) applied to the catalogue curve.
 */
export function fenderAtEnergy(curve, E, { fE = 1, fR = 1 } = {}) {
  const Ecap = curve.E[curve.E.length - 1] * fE;
  const exceeded = E > Ecap;
  let dfl = curve.dMax;
  if (!exceeded) {
    const Es = curve.E.map((e) => e * fE);
    // E(δ) is monotone: invert by interpolation
    dfl = interp(Es, curve.d, E);
  }
  let Rpeak = 0;
  for (let i = 0; i < curve.d.length && curve.d[i] <= dfl + 1e-12; i++) Rpeak = Math.max(Rpeak, curve.R[i]);
  const Rat = interp(curve.d, curve.R, dfl);
  return { deflection: dfl, reaction: Rat * fR, peakReaction: Math.max(Rpeak, Rat) * fR, exceeded, energyCapacity: Ecap };
}
