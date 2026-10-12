// Intact and damage stability: GZ curves by re-clipping the hull mesh, IMO / MODU criteria, lost-buoyancy
// damage cases, progressive flooding, ballast planning and inclining-test reduction.
// Frames and sign conventions: see hydrostatics.js (heel positive = port side up, trim positive = bow down).

import { G, needPos, boundedLeastSquares } from './_wavehelp.js';
import { bodyIntegrals, clipIntegrals, bodyToEarth, earthToBody, draftForVolume, solveEquilibrium, boxMesh } from './hydrostatics.js';

const DEG = Math.PI / 180;

function compartmentMeshes(compartments) {
  return (compartments || []).map((c) => {
    const mesh = c.mesh || (c.box ? boxMesh(c.box) : null);
    if (!mesh) throw new RangeError('each compartment needs a mesh or a box');
    return { mesh, permeability: c.permeability === undefined ? 1 : c.permeability };
  });
}

/**
 * Righting-lever (GZ) curve at constant displacement, by re-clipping the hull at each heel angle.
 * With `freeTrim` the trim is re-balanced at each heel (longitudinal B–G coincidence); otherwise the trim is
 * held. GZ = side·(y_G - y_B) in the earth frame, positive when righting.
 * Reference: Biran & López-Pulido (2014) (statical stability at large angles of heel); IMO 2008 IS Code Part A §2.2 (curve definition).
 * Hydrostatic (calm-water) levers; no wave, wind or dynamic effects.
 * @param {object} o
 * @param {object} o.mesh closed hull surface, body coordinates, m
 * @param {number} o.mass kg
 * @param {number[]} o.cog centre of gravity, body coordinates, m
 * @param {number} [o.rho=1025] kg/m³
 * @param {number[]} [o.heelsDeg] heel magnitudes, deg, ascending (default 0…70 step 2.5)
 * @param {number} [o.side=1] +1 heels port side up, -1 starboard side up
 * @param {boolean} [o.freeTrim=false]
 * @param {number} [o.trim=0] fixed trim, rad (start value when freeTrim)
 * @param {{mesh?:object, box?:object, permeability?:number}[]} [o.lost] damaged compartments (lost buoyancy)
 * @param {{name?:string, point:number[]}[]} [o.openings] downflooding points, body coordinates, m
 * @returns {{heelDeg:number[], gz:number[], draft:number[], trimDeg:number[], openingHeight:number[],
 *   downfloodDeg:number|null}} deg, m, m, deg; minimum height of the openings above the waterline m; first heel
 *   at which an opening immerses (null if none within range)
 */
export function gzCurve({ mesh, mass, cog, rho = 1025, heelsDeg, side = 1, freeTrim = false, trim = 0, lost, openings }) {
  needPos(mass, 'mass');
  if (side !== 1 && side !== -1) throw new RangeError('side must be 1 or -1');
  const heels = heelsDeg || Array.from({ length: 29 }, (_, i) => 2.5 * i);
  const lostM = compartmentMeshes(lost);
  const out = { heelDeg: [], gz: [], draft: [], trimDeg: [], openingHeight: [], downfloodDeg: null };
  let tr = trim; let prevDraft;
  for (const hd of heels) {
    if (!(hd >= 0 && hd < 180)) throw new RangeError('heelsDeg must be within 0..180');
    const heel = side * hd * DEG;
    let draft;
    if (freeTrim) {
      let eq = solveEquilibrium({ mesh, mass, cog, rho, guess: { heel, trim: tr, draft: prevDraft }, fixHeel: true, lost: lostM });
      if (!eq.converged) eq = solveEquilibrium({ mesh, mass, cog, rho, guess: { heel, trim: tr }, fixHeel: true, lost: lostM });
      if (!eq.converged) throw new Error(`Free-trim equilibrium did not converge at heel ${hd} deg`);
      draft = eq.draft; tr = eq.trim; prevDraft = draft;
    } else {
      draft = draftForVolume({ mesh, volume: mass / rho, heel, trim: tr, lost: lostM });
    }
    const r = bodyIntegrals(mesh, lostM, draft, heel, tr);
    const cg = bodyToEarth(cog, draft, heel, tr);
    out.heelDeg.push(hd);
    out.gz.push(side * (cg[1] - r.My / r.V));
    out.draft.push(draft);
    out.trimDeg.push(tr / DEG);
    let hmin = Infinity;
    for (const op of openings || []) hmin = Math.min(hmin, bodyToEarth(op.point, draft, heel, tr)[2]);
    out.openingHeight.push(hmin);
  }
  if (openings && openings.length) {
    const h = out.openingHeight;
    if (h[0] <= 0) out.downfloodDeg = out.heelDeg[0];
    else {
      for (let i = 1; i < h.length; i++) {
        if (h[i] <= 0) {
          out.downfloodDeg = out.heelDeg[i - 1] + ((out.heelDeg[i] - out.heelDeg[i - 1]) * h[i - 1]) / (h[i - 1] - h[i]);
          break;
        }
      }
    }
  }
  return out;
}

function lin(xs, ys, x) {
  if (x <= xs[0]) return ys[0];
  const n = xs.length;
  if (x >= xs[n - 1]) return ys[n - 1];
  let i = 1;
  while (xs[i] < x) i++;
  return ys[i - 1] + ((ys[i] - ys[i - 1]) * (x - xs[i - 1])) / (xs[i] - xs[i - 1]);
}

/**
 * Area under a piecewise-linear lever curve between two heel angles.
 * @param {number[]} heelDeg deg, ascending
 * @param {number[]} gz m
 * @param {number} aDeg lower limit, deg
 * @param {number} bDeg upper limit, deg (clamped to the data range)
 * @returns {number} area, m·rad
 */
export function areaUnderGz(heelDeg, gz, aDeg, bDeg) {
  const a = Math.max(aDeg, heelDeg[0]);
  const b = Math.min(bDeg, heelDeg[heelDeg.length - 1]);
  if (!(b > a)) return 0;
  const xs = [a];
  for (const h of heelDeg) if (h > a && h < b) xs.push(h);
  xs.push(b);
  let s = 0;
  for (let i = 1; i < xs.length; i++) s += 0.5 * (lin(heelDeg, gz, xs[i - 1]) + lin(heelDeg, gz, xs[i])) * (xs[i] - xs[i - 1]);
  return s * DEG;
}

/**
 * Characteristic values of a GZ curve: maximum lever and its angle, range of positive stability, areas,
 * initial GM (slope at the first interval) and static equilibrium angle when GZ(0) < 0.
 * Reference: IMO 2008 IS Code Part A §2.2 (quantities used by the criteria).
 * Piecewise-linear interpolation of the supplied points; use a heel step of 2.5° or finer.
 * @param {{heelDeg:number[], gz:number[], downfloodDeg?:number|null}} o deg, m, deg
 * @returns {{gzMax:number, heelAtGzMaxDeg:number, vanishingDeg:number|null, equilibriumDeg:number, rangeDeg:number|null,
 *   area0to30:number, area0to40:number, area30to40:number, gzAt30:number, gm0:number, limitDeg:number}}
 *   m, deg, deg (null when GZ stays positive to the end of the data), deg, deg, m·rad ×3 (the 40° limit is
 *   replaced by the downflooding angle when smaller: limitDeg), m, m, deg
 */
export function gzProperties({ heelDeg, gz, downfloodDeg = null }) {
  if (!heelDeg || heelDeg.length < 3 || heelDeg.length !== gz.length) throw new RangeError('heelDeg and gz must be equal-length arrays (>= 3 points)');
  let gzMax = -Infinity; let iMax = 0;
  for (let i = 0; i < gz.length; i++) if (gz[i] > gzMax) { gzMax = gz[i]; iMax = i; }
  let heelAtMax = heelDeg[iMax];
  if (iMax > 0 && iMax < gz.length - 1) {
    // parabolic refinement through the three points around the maximum
    const x0 = heelDeg[iMax - 1]; const x1 = heelDeg[iMax]; const x2 = heelDeg[iMax + 1];
    const y0 = gz[iMax - 1]; const y1 = gz[iMax]; const y2 = gz[iMax + 1];
    const d = (x0 - x1) * (x0 - x2) * (x1 - x2);
    const A = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / d;
    const B = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / d;
    if (A < 0) {
      const xv = -B / (2 * A);
      if (xv > x0 && xv < x2) { heelAtMax = xv; gzMax = Math.max(gzMax, y1 - A * (x1 - xv) ** 2); }
    }
  }
  let eq = heelDeg[0];
  let i0 = 0;
  if (gz[0] < 0) {
    for (let i = 1; i < gz.length; i++) {
      if (gz[i] >= 0) { eq = heelDeg[i - 1] + ((heelDeg[i] - heelDeg[i - 1]) * -gz[i - 1]) / (gz[i] - gz[i - 1]); i0 = i; break; }
    }
  }
  let vanishing = null;
  for (let i = Math.max(i0, 1); i < gz.length; i++) {
    if (gz[i] < 0 && gz[i - 1] >= 0 && heelDeg[i - 1] >= eq) {
      vanishing = heelDeg[i - 1] + ((heelDeg[i] - heelDeg[i - 1]) * gz[i - 1]) / (gz[i - 1] - gz[i]);
      break;
    }
  }
  const limit = downfloodDeg !== null && downfloodDeg !== undefined && downfloodDeg < 40 ? downfloodDeg : 40;
  return {
    gzMax, heelAtGzMaxDeg: heelAtMax, vanishingDeg: vanishing, equilibriumDeg: eq,
    rangeDeg: vanishing === null ? null : vanishing - eq,
    area0to30: areaUnderGz(heelDeg, gz, 0, 30), area0to40: areaUnderGz(heelDeg, gz, 0, limit),
    area30to40: areaUnderGz(heelDeg, gz, 30, limit), gzAt30: lin(heelDeg, gz, 30),
    gm0: (gz[1] - gz[0]) / ((heelDeg[1] - heelDeg[0]) * DEG), limitDeg: limit,
  };
}

/**
 * General intact stability criteria for ships: area 0–30° >= 0.055 m·rad; area 0–40° (or downflooding angle if
 * less) >= 0.090 m·rad; area 30–40° (or downflooding) >= 0.030 m·rad; GZ >= 0.20 m at an angle >= 30°; maximum
 * GZ at an angle >= 25°; initial GM >= 0.15 m.
 * Reference: IMO 2008 IS Code (MSC.267(85)) Part A §2.2.1–2.2.4. The weather criterion (§2.3) is not included.
 * @param {{heelDeg:number[], gz:number[], gm0?:number, downfloodDeg?:number|null}} o deg, m, m (initial GM corrected
 *   for free surfaces; default: slope of the curve), deg
 * @returns {{pass:boolean, checks:{name:string, value:number, limit:number, unit:string, ok:boolean, ref:string}[]}}
 */
export function imoIntactCriteria({ heelDeg, gz, gm0, downfloodDeg = null }) {
  const p = gzProperties({ heelDeg, gz, downfloodDeg });
  let gz30plus = -Infinity;
  for (let i = 0; i < gz.length; i++) if (heelDeg[i] >= 30) gz30plus = Math.max(gz30plus, gz[i]);
  gz30plus = Math.max(gz30plus, p.gzAt30);
  const gm = gm0 === undefined ? p.gm0 : gm0;
  const mk = (name, value, limit, unit, ref) => ({ name, value, limit, unit, ok: value >= limit, ref });
  const checks = [
    mk('Area under GZ, 0–30°', p.area0to30, 0.055, 'm·rad', 'IS Code 2008 Part A 2.2.1'),
    mk(`Area under GZ, 0–${p.limitDeg.toFixed(1)}°`, p.area0to40, 0.09, 'm·rad', 'IS Code 2008 Part A 2.2.1'),
    mk(`Area under GZ, 30–${p.limitDeg.toFixed(1)}°`, p.area30to40, 0.03, 'm·rad', 'IS Code 2008 Part A 2.2.1'),
    mk('GZ at an angle of heel >= 30°', gz30plus, 0.2, 'm', 'IS Code 2008 Part A 2.2.2'),
    mk('Angle of maximum GZ', p.heelAtGzMaxDeg, 25, 'deg', 'IS Code 2008 Part A 2.2.3'),
    mk('Initial metacentric height GM0', gm, 0.15, 'm', 'IS Code 2008 Part A 2.2.4'),
  ];
  return { pass: checks.every((c) => c.ok), checks };
}

/**
 * Wind heeling force, moment and upright lever from projected areas:
 * F = ½ Cs Ch rho_air V² A summed over the exposed blocks; lever = sum(F_i (z_i - z_lat)) / (displacement·g).
 * Reference: IMO MODU Code 2009 §3.2.3 (rho_air = 1.222 kg/m³; 36 m/s normal operating, 51.5 m/s severe storm).
 * Shape coefficient Cs (e.g. 0.4 spherical, 0.5 cylindrical, 1.0 large flat surface) and height coefficient Ch
 * (1.0 up to 15.3 m, increasing with height) come from the tables of that clause and are inputs here.
 * @param {object} o
 * @param {number} o.V wind speed, m/s
 * @param {{A:number, z:number, Cs?:number, Ch?:number}[]} o.areas projected area m², centroid height above the
 *   waterline m, coefficients (default 1)
 * @param {number} o.zLateralResistance elevation of the centre of lateral resistance of the underwater body
 *   relative to the waterline, m (negative below)
 * @param {number} o.displacement kg
 * @param {number} [o.rhoAir=1.222] kg/m³
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{force:number, moment:number, lever:number}} N, N·m, m
 */
export function windHeelingLever({ V, areas, zLateralResistance, displacement, rhoAir = 1.222, g = G }) {
  if (!(V >= 0)) throw new RangeError('V must be >= 0');
  needPos(displacement, 'displacement');
  if (!Array.isArray(areas)) throw new RangeError('areas must be an array');
  if (!Number.isFinite(zLateralResistance)) throw new RangeError('zLateralResistance must be finite');
  let force = 0; let moment = 0;
  for (const a of areas) {
    const f = 0.5 * (a.Cs === undefined ? 1 : a.Cs) * (a.Ch === undefined ? 1 : a.Ch) * rhoAir * V * V * a.A;
    force += f; moment += f * (a.z - zLateralResistance);
  }
  return { force, moment, lever: moment / (displacement * g) };
}

/**
 * Intact wind-heeling stability criterion for mobile offshore units: the area under the righting-lever curve
 * up to the second intercept with the wind heeling lever (or the downflooding angle, whichever is less) must be
 * at least 1.4 times (surface and self-elevating units) or 1.3 times (column-stabilised units) the area under
 * the wind heeling lever curve to the same angle, and the righting lever must be positive from upright to the
 * second intercept.
 * Reference: IMO MODU Code 2009 §3.3.1–3.3.2. The variation of the heeling lever with heel is an input:
 * 'constant' (conservative default), 'cos' or 'cos2' (lever0·cosⁿ heel).
 * @param {object} o
 * @param {number[]} o.heelDeg deg, ascending from 0
 * @param {number[]} o.gz m
 * @param {number} o.lever0 upright wind heeling lever, m
 * @param {string} [o.variation='constant']
 * @param {string} [o.unitType='surface'] 'surface' | 'self-elevating' | 'column-stabilised'
 * @param {number|null} [o.downfloodDeg]
 * @returns {{firstInterceptDeg:number|null, secondInterceptDeg:number|null, limitDeg:number, areaRighting:number,
 *   areaHeeling:number, ratio:number, required:number, positiveRange:boolean, pass:boolean}} deg, m·rad, -
 */
export function moduWindCriterion({ heelDeg, gz, lever0, variation = 'constant', unitType = 'surface', downfloodDeg = null }) {
  if (!(lever0 >= 0)) throw new RangeError('lever0 must be >= 0');
  const pow = { constant: 0, cos: 1, cos2: 2 }[variation];
  if (pow === undefined) throw new RangeError("variation must be 'constant', 'cos' or 'cos2'");
  const required = { surface: 1.4, 'self-elevating': 1.4, 'column-stabilised': 1.3 }[unitType];
  if (required === undefined) throw new RangeError("unitType must be 'surface', 'self-elevating' or 'column-stabilised'");
  const hl = (d) => lever0 * Math.cos(d * DEG) ** pow;
  const end = heelDeg[heelDeg.length - 1];
  const n = 20 * (heelDeg.length - 1);
  const xs = []; const dr = []; const hs = []; const gs = [];
  for (let i = 0; i <= n; i++) {
    const x = heelDeg[0] + ((end - heelDeg[0]) * i) / n;
    xs.push(x); gs.push(lin(heelDeg, gz, x)); hs.push(hl(x)); dr.push(gs[i] - hs[i]);
  }
  let first = dr[0] >= 0 ? xs[0] : null; let second = null;
  for (let i = 1; i <= n; i++) {
    if (first === null && dr[i - 1] < 0 && dr[i] >= 0) first = xs[i - 1] + ((xs[i] - xs[i - 1]) * -dr[i - 1]) / (dr[i] - dr[i - 1]);
    else if (first !== null && second === null && dr[i - 1] >= 0 && dr[i] < 0) second = xs[i - 1] + ((xs[i] - xs[i - 1]) * dr[i - 1]) / (dr[i - 1] - dr[i]);
  }
  let limit = second === null ? end : second;
  if (downfloodDeg !== null && downfloodDeg !== undefined) limit = Math.min(limit, downfloodDeg);
  const areaRighting = areaUnderGz(xs, gs, heelDeg[0], limit);
  const areaHeeling = areaUnderGz(xs, hs, heelDeg[0], limit);
  let positive = true;
  for (let i = 1; i <= n; i++) if (xs[i] <= limit && gs[i] < 0) positive = false;
  const ratio = areaHeeling > 0 ? areaRighting / areaHeeling : Infinity;
  return {
    firstInterceptDeg: first, secondInterceptDeg: second, limitDeg: limit, areaRighting, areaHeeling, ratio, required,
    positiveRange: positive, pass: first !== null && positive && ratio >= required,
  };
}

/**
 * Damage stability by the lost-buoyancy (constant-displacement) method: the buoyancy of the flooded
 * compartments (times permeability) is removed from the hull, mass and centre of gravity are unchanged.
 * Returns the damaged equilibrium, the residual GZ curve towards the damaged side and the margin-line check.
 * Reference: Biran & López-Pulido (2014) ch. 11; SOLAS II-1 Part B (lost-buoyancy convention).
 * Final-stage flooding only (compartments open to the sea and fully equalised); no intermediate stages,
 * cross-flooding or cargo shift.
 * @param {object} o
 * @param {object} o.mesh intact hull, closed surface, body coordinates, m
 * @param {number} o.mass kg
 * @param {number[]} o.cog body coordinates, m
 * @param {number} [o.rho=1025] kg/m³
 * @param {{mesh?:object, box?:{xMin:number,xMax:number,yMin:number,yMax:number,zMin:number,zMax:number},
 *   permeability?:number}[]} o.compartments flooded spaces (inside the hull), permeability 0..1 (default 1)
 * @param {number[]} [o.heelsDeg] heel magnitudes for the residual GZ curve, deg
 * @param {number[][]} [o.marginLine] points of the margin line, body coordinates, m
 * @param {{name?:string, point:number[]}[]} [o.openings] downflooding points
 * @returns {{equilibrium:{draft:number, heel:number, trim:number, heelDeg:number, trimDeg:number, converged:boolean},
 *   side:number, curve:object, properties:object, marginLine:{minHeight:number, submerged:boolean}|null,
 *   lostVolume:number}} m, rad, deg; curve as gzCurve; properties as gzProperties; margin-line height above the
 *   damaged waterline m; buoyant volume lost at the damaged waterline m³
 */
export function damageStability({ mesh, mass, cog, rho = 1025, compartments, heelsDeg, marginLine, openings }) {
  const lost = compartmentMeshes(compartments);
  if (lost.length === 0) throw new RangeError('compartments must list at least one flooded space');
  const eq = solveEquilibrium({ mesh, mass, cog, rho, lost });
  if (!eq.converged) throw new Error('Damaged equilibrium did not converge (vessel may capsize or sink in this damage case)');
  const side = eq.heel < -1e-9 ? -1 : 1;
  const curve = gzCurve({ mesh, mass, cog, rho, heelsDeg, side, freeTrim: true, trim: eq.trim, lost, openings });
  const properties = gzProperties({ heelDeg: curve.heelDeg, gz: curve.gz, downfloodDeg: curve.downfloodDeg });
  let ml = null;
  if (marginLine && marginLine.length) {
    let h = Infinity;
    for (const p of marginLine) h = Math.min(h, bodyToEarth(p, eq.draft, eq.heel, eq.trim)[2]);
    ml = { minHeight: h, submerged: h < 0 };
  }
  let lostVolume = 0;
  for (const c of lost) lostVolume += c.permeability * clipIntegrals(c.mesh, eq.draft, eq.heel, eq.trim).V;
  return {
    equilibrium: { draft: eq.draft, heel: eq.heel, trim: eq.trim, heelDeg: eq.heel / DEG, trimDeg: eq.trim / DEG, converged: eq.converged },
    side, curve, properties, marginLine: ml, lostVolume,
  };
}

// Internal water level (earth z) and centroid of floodwater of given volume inside a compartment.
function floodLevel(cmesh, volume, draft, heel, trim) {
  let lo = Infinity; let hi = -Infinity;
  const p = cmesh.positions;
  for (let i = 0; i < p.length; i += 3) {
    const z = bodyToEarth([p[i], p[i + 1], p[i + 2]], draft, heel, trim)[2];
    lo = Math.min(lo, z); hi = Math.max(hi, z);
  }
  const top = hi;
  for (let i = 0; i < 60; i++) {
    const mid = 0.5 * (lo + hi);
    if (clipIntegrals(cmesh, draft + mid, heel, trim).V < volume) lo = mid; else hi = mid;
  }
  const level = 0.5 * (lo + hi);
  const r = clipIntegrals(cmesh, draft + level, heel, trim);
  const c = r.V > 0 ? [r.Mx / r.V, r.My / r.V, r.Mz / r.V + level] : [0, 0, level];
  return { level: Math.min(level, top), centroidBody: earthToBody(c, draft, heel, trim) };
}

/**
 * Progressive flooding of one compartment through an opening to the sea, with quasi-static re-equilibration
 * (added-weight method): Q = Cd A sqrt(2 g |Δh|) (Torricelli), Δh = external head minus internal head at the
 * opening; at each step the floodwater mass and its centroid are added to the ship and draft, heel and trim are
 * re-solved. Explicit Euler in time.
 * Reference: Bernoulli's equation for an orifice with a discharge coefficient (Torricelli's law) and
 * hydrostatic equilibrium. Quasi-static: no sloshing, air compression, inertia of the flow or transient roll;
 * a single compartment and opening.
 * End states. The march stops when the levels equalise (`equalised`) or when no equilibrium is found for the
 * flooded ship. `sunk` is set in every such case, as it always was; the cause is reported separately:
 *  - `lostBuoyancy`: ship plus floodwater weigh more than the whole closed hull can displace — foundering;
 *  - `capsized`: buoyancy is sufficient but the only equilibrium found (searched from heel guesses of ±45°, ±90°,
 *    ±135° and 180° about the last state) lies beyond 40° of heel or trim from upright;
 *  - `solverFailed`: buoyancy is sufficient and no equilibrium was found at all, or one exists near the last
 *    attitude that the step-to-step solve missed — a numerical failure, not a prediction: shorten dt.
 * Exactly one of the three is true when `sunk` is true.
 * @param {object} o
 * @param {object} o.mesh intact hull, closed surface, m
 * @param {number} o.mass kg (intact)
 * @param {number[]} o.cog body coordinates, m
 * @param {{mesh?:object, box?:object, permeability?:number}} o.compartment flooded space
 * @param {{point:number[], area:number, Cd?:number}} o.opening centre (body coordinates, m), area m², discharge
 *   coefficient (typical 0.6 for a sharp-edged orifice; input)
 * @param {number} o.dt time step, s
 * @param {number} o.tMax simulated duration, s
 * @param {number} [o.rho=1025] kg/m³
 * @param {number} [o.g=9.80665] m/s²
 * @returns {{t:number[], floodMass:number[], draft:number[], heelDeg:number[], trimDeg:number[], level:number[],
 *   head:number[], flow:number[], equalised:boolean, sunk:boolean, lostBuoyancy:boolean, capsized:boolean,
 *   solverFailed:boolean, capsizeHeelDeg:number|null, timeToEqualise:number|null}} s, kg, m, deg, deg,
 *   internal water level relative to the sea surface m, driving head m, inflow m³/s; sunk = no equilibrium found
 *   for the flooded ship (any cause), split into lostBuoyancy / capsized / solverFailed (see above);
 *   capsizeHeelDeg = heel of the capsized equilibrium, deg
 */
export function progressiveFlooding({ mesh, mass, cog, compartment, opening, dt, tMax, rho = 1025, g = G }) {
  needPos(mass, 'mass'); needPos(dt, 'dt'); needPos(tMax, 'tMax');
  if (!opening || !(opening.area > 0)) throw new RangeError('opening.area must be > 0');
  const Cd = opening.Cd === undefined ? 0.6 : opening.Cd;
  if (!(Cd > 0 && Cd <= 1)) throw new RangeError('opening.Cd must be within (0, 1]');
  const [c] = compartmentMeshes([compartment]);
  const vMax = c.permeability * clipIntegrals(c.mesh, 1e5, 0, 0).V;
  const out = { t: [], floodMass: [], draft: [], heelDeg: [], trimDeg: [], level: [], head: [], flow: [], equalised: false, sunk: false,
    lostBuoyancy: false, capsized: false, solverFailed: false, capsizeHeelDeg: null, timeToEqualise: null };
  const hullVolume = clipIntegrals(mesh, 1e5, 0, 0).V;
  // why did the equilibrium solve fail for total mass mt with centre of gravity cg, coming from state s0?
  const classify = (mt, cg, s0) => {
    out.sunk = true;
    if (mt >= rho * hullVolume * (1 - 1e-6)) { out.lostBuoyancy = true; return; }
    const wrapPi = (a) => Math.atan2(Math.sin(a), Math.cos(a));
    let near = false;
    // towards the side the ship is already heeling to first
    const sg = s0.heel < 0 ? -1 : 1;
    for (const dh of [sg * Math.PI / 4, sg * Math.PI / 2, sg * 0.75 * Math.PI, Math.PI, -sg * Math.PI / 4, -sg * Math.PI / 2, -sg * 0.75 * Math.PI]) {
      let eq;
      try { eq = solveEquilibrium({ mesh, mass: mt, cog: cg, rho, guess: { heel: s0.heel + dh, trim: s0.trim } }); } catch (err) { continue; }
      if (!eq || !eq.converged) continue;
      const heel = wrapPi(eq.heel); const trim = wrapPi(eq.trim);
      if (Math.abs(heel) > 40 * DEG || Math.abs(trim) > 40 * DEG) { out.capsized = true; out.capsizeHeelDeg = heel / DEG; return; }
      near = true;
    }
    // an upright-range equilibrium exists (near) or nothing was found: either way the march failed numerically
    out.solverFailed = true;
    return near;
  };
  let mw = 0;
  let st = solveEquilibrium({ mesh, mass, cog, rho });
  if (!st.converged) throw new Error('Intact equilibrium did not converge');
  let cw = [opening.point[0], opening.point[1], opening.point[2]];
  let prevSign = 0;
  const nSteps = Math.ceil(tMax / dt);
  for (let n = 0; n <= nSteps; n++) {
    if (mw > 0) {
      for (let pass = 0; pass < 2; pass++) {
        const mt = mass + mw;
        const cg = [0, 1, 2].map((i) => (mass * cog[i] + mw * cw[i]) / mt);
        let eq;
        try {
          eq = solveEquilibrium({ mesh, mass: mt, cog: cg, rho, guess: { heel: st.heel, trim: st.trim, draft: st.draft } });
          if (!eq.converged) eq = solveEquilibrium({ mesh, mass: mt, cog: cg, rho, guess: { heel: st.heel, trim: st.trim } });
        } catch (err) { classify(mt, cg, st); break; }
        if (!eq.converged) { classify(mt, cg, st); break; }
        st = eq;
        cw = floodLevel(c.mesh, mw / (rho * c.permeability), st.draft, st.heel, st.trim).centroidBody;
      }
      if (out.sunk) break;
    }
    const zo = bodyToEarth(opening.point, st.draft, st.heel, st.trim)[2];
    const level = mw > 0 ? floodLevel(c.mesh, mw / (rho * c.permeability), st.draft, st.heel, st.trim).level : -Infinity;
    const hOut = Math.max(-zo, 0);
    const hIn = Math.max(level - zo, 0);
    const dh = hOut - hIn;
    const full = mw >= rho * vMax * (1 - 1e-9);
    const sign = Math.sign(dh);
    const Q = full && dh > 0 ? 0 : sign * Cd * opening.area * Math.sqrt(2 * g * Math.abs(dh));
    out.t.push(n * dt); out.floodMass.push(mw); out.draft.push(st.draft); out.heelDeg.push(st.heel / DEG);
    out.trimDeg.push(st.trim / DEG); out.level.push(level); out.head.push(dh); out.flow.push(Q);
    if (mw > 0 && (Q === 0 || (prevSign !== 0 && sign !== prevSign))) {
      out.equalised = true; out.timeToEqualise = n * dt; break;
    }
    prevSign = sign;
    mw = Math.min(Math.max(mw + rho * Q * dt, 0), rho * vMax);
  }
  return out;
}

/**
 * Ballast plan: tank volume changes that bring the vessel from its current draft/heel/trim to target values,
 * by bounded linear least squares on the small-angle hydrostatic model
 *   d(draft at LCF) = sum(dm)/(rho Aw),  d(heel) = -sum(dm·y)/(Δ GMt),  d(trim) = sum(dm·(x - x_F))/(Δ GMl),
 * with each tank limited to 0 … capacity, and the pumping time from the volume moved and the pump capacity.
 * 'transfer' keeps the total ballast constant (internal transfer); 'fill-discharge' exchanges with the sea.
 * Reference: small-angle hydrostatics, Biran & López-Pulido (2014) ch. 7 (weight shifts, trim). The vertical
 * shift of G, free-surface changes and the change of hydrostatics with draft are neglected: re-check the final
 * condition with meshHydrostatics / solveEquilibrium.
 * @param {object} o
 * @param {{name:string, x:number, y:number, capacity:number, volume:number, rho?:number}[]} o.tanks centroid
 *   plan position (body frame) m, capacity m³, current content m³, liquid density kg/m³ (default 1025)
 * @param {{displacement:number, waterplaneArea:number, lcf:number, GMt:number, GMl:number, rho?:number}} o.hydro
 *   kg, m², x of the centre of flotation m, m, m, sea-water density kg/m³
 * @param {{draft:number, heelDeg:number, trimDeg:number}} o.current m (at LCF), deg, deg
 * @param {{draft?:number, heelDeg?:number, trimDeg?:number}} o.target omitted quantities are left free
 * @param {string} [o.mode='transfer']
 * @param {number} [o.pumpRate] capacity of one pump, m³/s
 * @param {number} [o.nPumps=1]
 * @returns {{tanks:{name:string, deltaVolume:number, newVolume:number, deltaMass:number}[], achieved:{draft:number,
 *   heelDeg:number, trimDeg:number}, error:{draft:number|null, heelDeg:number|null, trimDeg:number|null},
 *   volumeMoved:number, pumpTime:number|null, feasible:boolean}} m³, kg, m, deg, m³, s; feasible when every target
 *   is met within 1 mm / 0.01°
 */
export function ballastPlan({ tanks, hydro, current, target, mode = 'transfer', pumpRate, nPumps = 1 }) {
  if (!Array.isArray(tanks) || tanks.length === 0) throw new RangeError('tanks must be a non-empty array');
  if (mode !== 'transfer' && mode !== 'fill-discharge') throw new RangeError("mode must be 'transfer' or 'fill-discharge'");
  for (const key of ['displacement', 'waterplaneArea', 'GMt', 'GMl']) needPos(hydro[key], `hydro.${key}`);
  const rho = hydro.rho === undefined ? 1025 : hydro.rho;
  const n = tanks.length;
  const rt = tanks.map((t) => (t.rho === undefined ? 1025 : t.rho));
  // rows scaled so that 1 mm of draft and 0.01° of heel/trim have equal weight
  const sD = 1 / 1e-3; const sA = 1 / (0.01 * DEG);
  const rowD = tanks.map((t, i) => rt[i] / (rho * hydro.waterplaneArea));
  const rowH = tanks.map((t, i) => (-rt[i] * t.y) / (hydro.displacement * hydro.GMt));
  const rowT = tanks.map((t, i) => (rt[i] * (t.x - hydro.lcf)) / (hydro.displacement * hydro.GMl));
  const A = []; const b = [];
  if (target.draft !== undefined) { A.push(rowD.map((v) => v * sD)); b.push((target.draft - current.draft) * sD); }
  if (target.heelDeg !== undefined) { A.push(rowH.map((v) => v * sA)); b.push((target.heelDeg - current.heelDeg) * DEG * sA); }
  if (target.trimDeg !== undefined) { A.push(rowT.map((v) => v * sA)); b.push((target.trimDeg - current.trimDeg) * DEG * sA); }
  if (A.length === 0) throw new RangeError('target must specify at least one of draft, heelDeg, trimDeg');
  if (mode === 'transfer') { A.push(rowD.map((v) => v * sD * 30)); b.push(0); }
  const lo = []; const hi = [];
  for (const t of tanks) {
    if (!(t.capacity > 0) || !(t.volume >= 0 && t.volume <= t.capacity * (1 + 1e-12))) throw new RangeError(`tank ${t.name}: need capacity > 0 and 0 <= volume <= capacity`);
    lo.push(-t.volume); hi.push(t.capacity - t.volume);
  }
  let tr = 0;
  for (const row of A) for (const v of row) tr += v * v;
  const sol = boundedLeastSquares(A, b, lo, hi, { lambda: (1e-10 * tr) / n, maxIter: 200000, tol: 1e-13 });
  const dv = sol.x;
  const sum = (row) => row.reduce((s, v, i) => s + v * dv[i], 0);
  const achieved = {
    draft: current.draft + sum(rowD), heelDeg: current.heelDeg + sum(rowH) / DEG, trimDeg: current.trimDeg + sum(rowT) / DEG,
  };
  const error = {
    draft: target.draft === undefined ? null : achieved.draft - target.draft,
    heelDeg: target.heelDeg === undefined ? null : achieved.heelDeg - target.heelDeg,
    trimDeg: target.trimDeg === undefined ? null : achieved.trimDeg - target.trimDeg,
  };
  let pos = 0; let neg = 0;
  for (const v of dv) { if (v > 0) pos += v; else neg -= v; }
  const volumeMoved = mode === 'transfer' ? Math.max(pos, neg) : pos + neg;
  const feasible = (error.draft === null || Math.abs(error.draft) < 1e-3)
    && (error.heelDeg === null || Math.abs(error.heelDeg) < 0.01) && (error.trimDeg === null || Math.abs(error.trimDeg) < 0.01);
  return {
    tanks: tanks.map((t, i) => ({ name: t.name, deltaVolume: dv[i], newVolume: t.volume + dv[i], deltaMass: rt[i] * dv[i] })),
    achieved, error, volumeMoved, pumpTime: pumpRate > 0 ? volumeMoved / (pumpRate * nPumps) : null, feasible,
  };
}

/**
 * Inclining-test reduction: GM from transverse weight shifts, w d = Δ GM tan(phi), with tan(phi) from pendulum
 * deflection / pendulum length; several shifts are combined by least squares through the origin.
 * Reference: IMO 2008 IS Code Part B Annex 1 (inclining test); Biran & López-Pulido (2014) ch. 7.
 * Valid for small heel (typically < 4°) with wall-sided waterlines; correct the result for free surfaces.
 * @param {object} o
 * @param {number} o.displacement as-inclined displacement, kg
 * @param {{mass:number, distance:number, deflection:number, pendulumLength:number}[]} o.shifts shifted mass kg,
 *   transverse distance m, pendulum deflection m, pendulum length m
 * @param {number} [o.KM] transverse KM at the test draft, m (gives KG)
 * @returns {{GM:number, KG:number|null, perShift:{tanPhi:number, heelDeg:number, GM:number}[], spread:number}}
 *   m; spread = (max - min)/mean of the single-shift GM values, -
 */
export function incliningTest({ displacement, shifts, KM }) {
  needPos(displacement, 'displacement');
  if (!Array.isArray(shifts) || shifts.length === 0) throw new RangeError('shifts must be a non-empty array');
  let num = 0; let den = 0;
  const perShift = shifts.map((s) => {
    needPos(s.pendulumLength, 'pendulumLength');
    const tanPhi = s.deflection / s.pendulumLength;
    if (tanPhi === 0) throw new RangeError('deflection must be non-zero');
    const moment = s.mass * s.distance;
    num += moment * tanPhi; den += tanPhi * tanPhi;
    return { tanPhi, heelDeg: Math.atan(tanPhi) / DEG, GM: moment / (displacement * tanPhi) };
  });
  const GM = num / (displacement * den);
  const gms = perShift.map((p) => p.GM);
  return { GM, KG: KM === undefined ? null : KM - GM, perShift, spread: (Math.max(...gms) - Math.min(...gms)) / GM };
}
