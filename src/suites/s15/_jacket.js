// Planar rigid-body strip model of a jacket for the launch and upending simulations.
//
// Frame: x horizontal (positive astern / in the launch direction), z up, still water at z = 0. The jacket is a
// rigid line through its centre of gravity c = (x, z); θ is the angle of its axis below the horizontal, so the
// axis unit vector is u = (cos θ, −sin θ) and the normal is n = (sin θ, cos θ). A point at axial station ξ
// (positive towards the leading / lower end) is p = c + ξ·u, moves with v − ξ·ω·n, and a force F applied there
// does virtual work F·δc − ξ·(F·n)·δθ. Buoyancy, drag and added mass are integrated over strips with a linear
// buoyant-volume distribution fixed by the total volume and its centroid. Added mass enters as inertia only
// (no slam or Munk terms). This is a reduced-order model: see the module limits.
import { G, need } from './_common.js';

export const jacketInputs = ({ M = 12000, L = 140, xiG = 0.55 * L, xiB = 0.52 * L, reserve = 13, group = 'Jacket' } = {}) => [
  { key: 'M', label: 'Jacket mass (incl. launch/upending appurtenances)', unit: 't', type: 'number', default: M, min: 50, max: 80000, group, role: 'load', uncertain: { dist: 'normal', cov: 0.03 } },
  { key: 'L', label: 'Jacket length along its axis', unit: 'm', type: 'number', default: L, min: 10, max: 500, group, role: 'geometry' },
  { key: 'xiG', label: 'Centre of gravity from the trailing (upper) end', unit: 'm', type: 'number', default: xiG, min: 1, max: 499, group, role: 'geometry', uncertain: { dist: 'normal', cov: 0.02 }, help: 'Measured along the axis from the end that leaves the barge last (launch) or ends up on top (upending).' },
  { key: 'kGyr', label: 'Pitch radius of gyration about the CoG', unit: 'm', type: 'number', default: 0.3 * L, min: 1, max: 250, group, role: 'geometry', help: 'About 0.29·L for a uniform bar; jackets are typically 0.27–0.33·L.' },
  { key: 'reserve', label: 'Reserve buoyancy (buoyant volume × ρ / mass − 1)', unit: '%', type: 'number', default: reserve, min: -30, max: 100, group, role: 'load', design: true, uncertain: { dist: 'normal', cov: 0.15 }, help: 'Total sealed buoyancy of legs, braces and auxiliary tanks over the mass. Launch practice asks for roughly 10–15 % after allowing for a damaged compartment.' },
  { key: 'xiB', label: 'Centre of buoyancy from the trailing (upper) end', unit: 'm', type: 'number', default: xiB, min: 1, max: 499, group, role: 'geometry', help: 'Centroid of the total buoyant volume along the axis. Must lie within the middle third for the linear volume distribution used here.' },
  { key: 'Dn', label: 'Broadside projected width (drag area per unit length)', unit: 'm', type: 'number', default: 18, min: 0.2, max: 200, group: 'Hydrodynamics', role: 'geometry', help: 'Sum of the projected member widths seen normal to the axis, averaged along the length.' },
  { key: 'Cdn', label: 'Normal drag coefficient', unit: '-', type: 'number', default: 1.0, min: 0.2, max: 3, group: 'Hydrodynamics', role: 'load', uncertain: { dist: 'normal', cov: 0.2 } },
  { key: 'Aax', label: 'End-on projected area', unit: 'm²', type: 'number', default: 350, min: 0, max: 20000, group: 'Hydrodynamics', role: 'geometry', help: 'Projected area of braces and legs seen along the axis; gives the axial drag.' },
  { key: 'Cdt', label: 'Axial drag coefficient', unit: '-', type: 'number', default: 1.0, min: 0, max: 3, group: 'Hydrodynamics', role: 'load' },
  { key: 'Ca', label: 'Added-mass coefficient (normal to the axis)', unit: '-', type: 'number', default: 1.0, min: 0, max: 2.5, group: 'Hydrodynamics', role: 'load', uncertain: { dist: 'normal', cov: 0.2 }, help: 'Applied to the displaced volume of each submerged strip; 1.0 for circular members.' },
];

/** Build the strip model from the jacket inputs (t → kg, % → fraction). */
export function buildJacket(inp, rho, nStrip = 40) {
  const M = inp.M * 1e3, L = inp.L, La = inp.xiG, Lf = L - inp.xiG;
  need(La > 0 && Lf > 0, 'The centre of gravity must lie within the jacket length.');
  const V = M * (1 + inp.reserve / 100) / rho, vBar = V / L;
  const off = inp.xiB - L / 2;
  need(Math.abs(off) <= L / 6 + 1e-9, `The centre of buoyancy must lie within the middle third of the length (${(L / 3).toFixed(1)}–${(2 * L / 3).toFixed(1)} m from the trailing end) for the linear buoyancy distribution.`);
  const slope = 12 * vBar * off / (L * L);
  const dxi = L / nStrip, xi = new Float64Array(nStrip), vol = new Float64Array(nStrip), rad = new Float64Array(nStrip);
  for (let i = 0; i < nStrip; i++) {
    const s = (i + 0.5) * dxi; // from the trailing end
    xi[i] = s - La;
    const a = Math.max(vBar + slope * (s - L / 2), 0);
    vol[i] = a * dxi; rad[i] = Math.sqrt(Math.max(a, 1e-6) / Math.PI);
  }
  return { M, I: M * inp.kGyr ** 2, L, La, Lf, V, n: nStrip, dxi, xi, vol, rad, rho, Dn: inp.Dn, Cdn: inp.Cdn, Aax: inp.Aax, Cdt: inp.Cdt, Ca: inp.Ca };
}

/**
 * Hydrostatic and hydrodynamic generalised forces on the jacket at state (x, z, θ, vx, vz, ω).
 * Adds to Q (length 3) and to the symmetric mass matrix Mm (3×3 as 9 numbers); returns the submerged volume.
 */
export function hydro(j, z, th, vx, vz, om, Q, Mm) {
  const s = Math.sin(th), c = Math.cos(th), ux = c, uz = -s, nx = s, nz = c;
  const half = 0.5 * j.dxi * Math.abs(s);
  let Vsub = 0;
  for (let i = 0; i < j.n; i++) {
    const xi = j.xi[i], zi = z + xi * uz;
    const h = j.rad[i] * Math.abs(c) + half; // half-extent of the strip in z
    const f = zi <= -h ? 1 : zi >= h ? 0 : (h - zi) / (2 * h);
    if (f <= 0) continue;
    const vi = j.vol[i] * f; Vsub += vi;
    const vn = vx * nx + vz * nz - xi * om, vt = vx * ux + vz * uz;
    const Fb = j.rho * G * vi;
    const Fn = -0.5 * j.rho * j.Cdn * j.Dn * j.dxi * f * Math.abs(vn) * vn;
    const Ft = -0.5 * j.rho * j.Cdt * (j.Aax / j.n) * f * Math.abs(vt) * vt;
    const Fx = Fn * nx + Ft * ux, Fz = Fb + Fn * nz + Ft * uz;
    Q[0] += Fx; Q[1] += Fz; Q[2] += -xi * (Fx * nx + Fz * nz);
    const ma = j.rho * j.Ca * vi, J0 = nx, J1 = nz, J2 = -xi;
    Mm[0] += ma * J0 * J0; Mm[1] += ma * J0 * J1; Mm[2] += ma * J0 * J2;
    Mm[4] += ma * J1 * J1; Mm[5] += ma * J1 * J2; Mm[8] += ma * J2 * J2;
  }
  return Vsub;
}

/** Add a force (Fx, Fz) acting at axial station xi. */
export function addForce(Q, th, xi, Fx, Fz) {
  Q[0] += Fx; Q[1] += Fz; Q[2] += -xi * (Fx * Math.sin(th) + Fz * Math.cos(th));
}
/** Add a point mass m at axial station xi to the mass matrix, with its weight and centripetal term in Q. */
export function addPointMass(Q, Mm, th, om, xi, m) {
  const nx = Math.sin(th), nz = Math.cos(th), ux = nz, uz = -nx;
  Mm[0] += m; Mm[4] += m; Mm[2] += -m * xi * nx; Mm[5] += -m * xi * nz; Mm[8] += m * xi * xi;
  Q[0] += m * xi * om * om * ux; Q[1] += m * xi * om * om * uz - m * G;
  Q[2] += xi * m * G * nz;
}
/** Solve the symmetric 3×3 system Mm·a = Q (upper triangle filled). */
export function solve3(Mm, Q) {
  const a = Mm[0], b = Mm[1], c = Mm[2], d = Mm[4], e = Mm[5], f = Mm[8];
  const det = a * (d * f - e * e) - b * (b * f - e * c) + c * (b * e - d * c);
  return [
    (Q[0] * (d * f - e * e) - b * (Q[1] * f - e * Q[2]) + c * (Q[1] * e - d * Q[2])) / det,
    (a * (Q[1] * f - e * Q[2]) - Q[0] * (b * f - e * c) + c * (b * Q[2] - Q[1] * c)) / det,
    (a * (d * Q[2] - Q[1] * e) - b * (b * Q[2] - Q[1] * c) + Q[0] * (b * e - d * c)) / det,
  ];
}
