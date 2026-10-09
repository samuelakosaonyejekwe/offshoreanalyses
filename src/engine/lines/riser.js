/**
 * Riser configuration builders and global analyses on the planar rod solver (rod.js): steel catenary risers,
 * lazy/steep/pliant-wave risers and dynamic cables, top-tensioned and drilling risers, free-standing hybrid
 * risers, vessel-offset and quasi-static heave cases, analytic boundary-layer cross-checks, effective-weight
 * bookkeeping and wall stresses.
 *
 * Conventions: x horizontal, z vertical upward; risers are described from the top (hang-off, s = 0) to the
 * seabed end. Angles named `...Deg` are in degrees from the downward vertical, positive towards +x.
 * All analyses are in-plane (2-D); "cross" offsets are not modelled.
 */

import { solveRodStatic } from './rod.js';
import { brent } from '../math/roots.js';

const G0 = 9.80665;
const DEG = Math.PI / 180;

function num(v, name, { min = -Infinity, max = Infinity, strictMin = false } = {}) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < min || v > max || (strictMin && v <= min)) {
    throw new RangeError(`riser: argument ${name} = ${v} is outside its valid range`);
  }
  return v;
}

/* ------------------------------------------------------------------------------------------------------- */
/* Section properties, effective tension, wall stress                                                       */

/**
 * Cross-section and weight bookkeeping of a pipe with contents, coatings and optional buoyancy modules.
 *
 * Apparent (submerged) weight per unit length = (steel + coating + module + contents mass)·g − ρext·g·(displaced
 * area) (Archimedes on the closed displaced volume; Sparks, "Fundamentals of Marine Riser Mechanics", Ch. 2).
 * Validity: concentric circular layers; bending and axial stiffness from the steel wall only.
 *
 * @param {object} o
 * @param {number} o.D steel outer diameter [m]
 * @param {number} o.t steel wall thickness [m]
 * @param {number} [o.E=207e9] Young's modulus [Pa]
 * @param {number} [o.rhoSteel=7850] steel density [kg/m³]
 * @param {number} [o.rhoInt=0] contents density [kg/m³]
 * @param {number} [o.rhoExt=1025] sea-water density [kg/m³]
 * @param {Array<{thickness: number, rho: number}>} [o.coatings] layers from the steel outwards [m, kg/m³]
 * @param {{ Dout: number, rho: number, fill?: number }} [o.buoyancy] module outer diameter [m], material density
 *   [kg/m³] and fraction of the length covered (default 1) [-]
 * @param {number} [o.Ca=1] normal added-mass coefficient on the outer displaced area [-]
 * @param {number} [o.g=9.80665] [m/s²]
 * @returns {{ As: number, Ai: number, Ae: number, I: number, EA: number, EI: number, Dhydro: number,
 *   massSteel: number, massCoating: number, massBuoyancy: number, massContents: number, m: number, ma: number,
 *   wDry: number, w: number, displacedArea: number }} areas [m²] (Ai, Ae are the steel-bore and steel-outer areas
 *   used for pressure end-cap terms), I [m⁴], EA [N], EI [N·m²], Dhydro [m], masses [kg/m], added mass ma [kg/m],
 *   dry weight and apparent weight w [N/m, + down], displaced area [m²]
 */
export function pipeProperties(o) {
  const { D, t, E = 207e9, rhoSteel = 7850, rhoInt = 0, rhoExt = 1025, coatings = [], buoyancy, Ca = 1, g = G0 } = o ?? {};
  num(D, 'D', { min: 0, strictMin: true }); num(t, 't', { min: 0, strictMin: true, max: D / 2 });
  num(E, 'E', { min: 0, strictMin: true }); num(rhoInt, 'rhoInt', { min: 0 }); num(rhoExt, 'rhoExt', { min: 0 });
  const Di = D - 2 * t;
  const Ae = Math.PI * D * D / 4, Ai = Math.PI * Di * Di / 4, As = Ae - Ai;
  const I = Math.PI * (D ** 4 - Di ** 4) / 64;
  let Do = D, massCoating = 0;
  coatings.forEach((c, k) => {
    num(c.thickness, `coatings[${k}].thickness`, { min: 0 }); num(c.rho, `coatings[${k}].rho`, { min: 0 });
    const Dn = Do + 2 * c.thickness;
    massCoating += c.rho * Math.PI * (Dn * Dn - Do * Do) / 4;
    Do = Dn;
  });
  let massBuoyancy = 0, dispArea = Math.PI * Do * Do / 4, Dhydro = Do;
  if (buoyancy) {
    num(buoyancy.Dout, 'buoyancy.Dout', { min: Do }); num(buoyancy.rho, 'buoyancy.rho', { min: 0 });
    const fill = num(buoyancy.fill ?? 1, 'buoyancy.fill', { min: 0, max: 1 });
    const Am = fill * Math.PI * (buoyancy.Dout ** 2 - Do * Do) / 4;
    massBuoyancy = buoyancy.rho * Am;
    dispArea += Am;
    // equal-volume hydrodynamic diameter for partially covered sections
    Dhydro = Math.sqrt(4 * dispArea / Math.PI);
  }
  const massSteel = rhoSteel * As, massContents = rhoInt * Ai;
  const m = massSteel + massCoating + massBuoyancy + massContents;
  return {
    As, Ai, Ae, I, EA: E * As, EI: E * I, Dhydro, massSteel, massCoating, massBuoyancy, massContents, m,
    ma: Ca * rhoExt * dispArea, wDry: (m - massContents) * g, w: (m - rhoExt * dispArea) * g, displacedArea: dispArea,
  };
}

/**
 * Effective tension from true wall tension: Te = Tw − pi·Ai + pe·Ae (Sparks 2007, Eq. 2.5; DNV-ST-F201 §4).
 * @param {object} o `{ Tw [N], pi [Pa], pe [Pa], Ai [m²], Ae [m²] }`
 * @returns {number} effective tension [N]
 */
export function effectiveTension({ Tw, pi = 0, pe = 0, Ai, Ae }) {
  num(Tw, 'Tw'); num(Ai, 'Ai', { min: 0 }); num(Ae, 'Ae', { min: 0 });
  return Tw - pi * Ai + pe * Ae;
}

/**
 * True wall tension from effective tension: Tw = Te + pi·Ai − pe·Ae (Sparks 2007, Eq. 2.5).
 * @param {object} o `{ Te [N], pi [Pa], pe [Pa], Ai [m²], Ae [m²] }`
 * @returns {number} true wall tension [N]
 */
export function trueWallTension({ Te, pi = 0, pe = 0, Ai, Ae }) {
  num(Te, 'Te'); num(Ai, 'Ai', { min: 0 }); num(Ae, 'Ae', { min: 0 });
  return Te + pi * Ai - pe * Ae;
}

/**
 * Wall stresses of a pressurised pipe under effective tension and bending: thick-wall Lamé hoop/radial stress,
 * axial stress from true wall tension ± bending, von Mises equivalent at the inner and outer fibres.
 * Reference: Lamé (Timoshenko & Goodier, "Theory of Elasticity", §28); Sparks (2007) Ch. 4.
 * Validity: linear elastic, plane sections, no shear or torsion.
 *
 * @param {object} o `{ Te [N], M [N·m], pi [Pa], pe [Pa], D [m], t [m] }`
 * @returns {{ Tw: number, axial: number, bendingOuter: number, hoopInner: number, hoopOuter: number,
 *   radialInner: number, radialOuter: number, vonMisesInner: number, vonMisesOuter: number, vonMises: number }}
 *   Tw [N]; stresses [Pa] (axial = Tw/As; vonMises = governing of the four fibre points)
 */
export function wallStresses({ Te, M = 0, pi = 0, pe = 0, D, t }) {
  num(Te, 'Te'); num(M, 'M'); num(D, 'D', { min: 0, strictMin: true }); num(t, 't', { min: 0, strictMin: true, max: D / 2 });
  const ro = D / 2, ri = ro - t;
  const Ae = Math.PI * ro * ro, Ai = Math.PI * ri * ri, As = Ae - Ai, I = Math.PI * (ro ** 4 - ri ** 4) / 4;
  const Tw = Te + pi * Ai - pe * Ae;
  const axial = Tw / As;
  const A = (pi * ri * ri - pe * ro * ro) / (ro * ro - ri * ri);
  const B = (pi - pe) * ri * ri * ro * ro / (ro * ro - ri * ri);
  const hoopInner = A + B / (ri * ri), hoopOuter = A + B / (ro * ro);
  const radialInner = -pi, radialOuter = -pe;
  const vm = (sa, sh, sr) => Math.sqrt(0.5 * ((sa - sh) ** 2 + (sh - sr) ** 2 + (sr - sa) ** 2));
  const bo = Math.abs(M) * ro / I, bi = Math.abs(M) * ri / I;
  const vonMisesInner = Math.max(vm(axial + bi, hoopInner, radialInner), vm(axial - bi, hoopInner, radialInner));
  const vonMisesOuter = Math.max(vm(axial + bo, hoopOuter, radialOuter), vm(axial - bo, hoopOuter, radialOuter));
  return { Tw, axial, bendingOuter: bo, hoopInner, hoopOuter, radialInner, radialOuter, vonMisesInner, vonMisesOuter, vonMises: Math.max(vonMisesInner, vonMisesOuter) };
}

/* ------------------------------------------------------------------------------------------------------- */
/* Analytic catenary and boundary layers                                                                     */

/**
 * Inextensible free-hanging catenary riser on a flat seabed (no bending stiffness): geometry and tensions from
 * the hang-off height and one of hang-off angle, top tension or horizontal tension.
 * Reference: Irvine, "Cable Structures", §1.2; Sparks (2007) App. catenary equations.
 * Validity: uniform apparent weight, no current, horizontal tangent at touchdown.
 *
 * @param {object} o
 * @param {number} o.h hang-off height above the seabed [m]
 * @param {number} o.w apparent weight per unit length [N/m]
 * @param {number} [o.hangoffAngleDeg] top angle from the vertical [deg]
 * @param {number} [o.topTension] top effective tension [N] (must exceed w·h)
 * @param {number} [o.H] horizontal tension [N]
 * @returns {{ H: number, a: number, topTension: number, V: number, suspendedLength: number,
 *   horizontalSpan: number, hangoffAngleDeg: number, tdpCurvature: number, topCurvature: number }}
 *   H, topTension, V (top vertical force) [N]; catenary parameter a, lengths [m]; curvature [1/m]
 */
export function catenaryRiser({ h, w, hangoffAngleDeg, topTension, H }) {
  num(h, 'h', { min: 0, strictMin: true }); num(w, 'w', { min: 0, strictMin: true });
  let a;
  if (H !== undefined) a = num(H, 'H', { min: 0, strictMin: true }) / w;
  else if (topTension !== undefined) a = (num(topTension, 'topTension', { min: w * h, strictMin: true }) - w * h) / w;
  else {
    const sa = Math.sin(num(hangoffAngleDeg, 'hangoffAngleDeg', { min: 0, strictMin: true, max: 89.999 }) * DEG);
    a = h * sa / (1 - sa);
  }
  const S = Math.sqrt(h * (h + 2 * a));
  const Hh = w * a, T = Hh + w * h;
  return {
    H: Hh, a, topTension: T, V: w * S, suspendedLength: S, horizontalSpan: a * Math.asinh(S / a),
    hangoffAngleDeg: Math.atan2(a, S) / DEG, tdpCurvature: 1 / a, topCurvature: a / (a + h) ** 2,
  };
}

/**
 * Boundary-layer (stiffened-catenary) bending solution at the touchdown point of a rod laid on a rigid flat
 * seabed: curvature rises from zero at touchdown to the catenary value over the decay length λ = sqrt(EI/H),
 * κ(ξ) = (w/H)(1 − e^(−ξ/λ)), with a concentrated seabed reaction w·λ at the touchdown point.
 * Reference: Plunkett (1967), "Static bending stresses in catenaries and drill strings", J. Eng. Ind. 89;
 * Dixon & Rutledge (1968); Sparks (2007) Ch. 13. Validity: λ much shorter than the suspended length, rigid seabed.
 *
 * @param {object} o `{ EI [N·m²], H [N horizontal (TDP) tension], w [N/m], xi? [m distances above touchdown] }`
 * @returns {{ lambda: number, catenaryCurvature: number, pointReaction: number, tdpShear: number,
 *   curvature: number[] }} lambda [m]; curvature [1/m] (at each xi); pointReaction, tdpShear [N]
 */
export function boundaryLayerTdp({ EI, H, w, xi = [] }) {
  num(EI, 'EI', { min: 0 }); num(H, 'H', { min: 0, strictMin: true }); num(w, 'w', { min: 0, strictMin: true });
  const lambda = Math.sqrt(EI / H), k0 = w / H;
  return {
    lambda, catenaryCurvature: k0, pointReaction: w * lambda, tdpShear: w * lambda,
    curvature: Array.from(xi, (x) => (lambda > 0 ? k0 * (1 - Math.exp(-x / lambda)) : k0)),
  };
}

/**
 * Boundary-layer bending at a hang-off (or any support) where the imposed end angle differs from the cable
 * (catenary) angle by Δθ: θ(s) − θcat = Δθ0·e^(−s/λ), λ = sqrt(EI/T); with a flex joint of rotational stiffness
 * kRot in series, end moment M0 = Δθ / (1/kRot + 1/sqrt(EI·T)).
 * Reference: Sparks (2007) Ch. 13 (tensioned-beam end rotation stiffness sqrt(EI·T)).
 * Validity: T > 0 and λ much shorter than the riser, small Δθ.
 *
 * @param {object} o `{ EI [N·m²], T [N effective tension at the end], deltaAngle [rad], kRot? [N·m/rad, omit for a
 *   rigid (clamped) end] }`
 * @returns {{ lambda: number, endMoment: number, endCurvature: number, jointRotation: number,
 *   rotationalStiffness: number }} lambda [m]; endMoment [N·m]; endCurvature [1/m]; jointRotation [rad] taken by the
 *   flex joint; rotationalStiffness sqrt(EI·T) [N·m/rad]
 */
export function boundaryLayerHangoff({ EI, T, deltaAngle, kRot }) {
  num(EI, 'EI', { min: 0, strictMin: true }); num(T, 'T', { min: 0, strictMin: true }); num(deltaAngle, 'deltaAngle');
  const kb = Math.sqrt(EI * T);
  const M0 = kRot === undefined ? deltaAngle * kb : deltaAngle / (1 / num(kRot, 'kRot', { min: 0, strictMin: true }) + 1 / kb);
  return { lambda: Math.sqrt(EI / T), endMoment: M0, endCurvature: M0 / EI, jointRotation: kRot === undefined ? 0 : M0 / kRot, rotationalStiffness: kb };
}

/**
 * Idealised seabed trench profile for touchdown-zone studies: a raised-cosine depression of given length and
 * maximum depth (shape is an assumption — measured trench profiles should be supplied as a custom function).
 * @param {object} o `{ zSeabed [m], xStart [m], length [m], depth [m] }`
 * @returns {(x: number) => number} seabed elevation [m] at x [m]
 */
export function trenchProfile({ zSeabed = 0, xStart, length, depth }) {
  num(xStart, 'xStart'); num(length, 'length', { min: 0, strictMin: true }); num(depth, 'depth', { min: 0 });
  return (x) => {
    const f = (x - xStart) / length;
    return f <= 0 || f >= 1 ? zSeabed : zSeabed - 0.5 * depth * (1 - Math.cos(2 * Math.PI * f));
  };
}

/* ------------------------------------------------------------------------------------------------------- */
/* Generic riser solve                                                                                       */

function endFromCondition(pos, cond) {
  const e = { x: pos.x, z: pos.z };
  if (cond && cond.type && cond.type !== 'pinned') {
    e.angle = num(cond.angleDeg, 'condition.angleDeg') * DEG - Math.PI / 2;
    if (cond.type === 'flex') e.kRot = num(cond.kRot, 'condition.kRot', { min: 0, strictMin: true });
    else if (cond.type !== 'clamped') throw new RangeError("riser: argument condition.type must be 'pinned', 'clamped' or 'flex'");
  }
  return e;
}

function summarise(r, cfg) {
  const n = r.s.length - 1;
  let maxK = 0, sMaxK = 0, maxT = -Infinity, minT = Infinity, maxM = 0;
  for (let k = 0; k <= n; k++) {
    const c = Math.abs(r.curvature[k]);
    if (c > maxK) { maxK = c; sMaxK = r.s[k]; }
    maxT = Math.max(maxT, r.Te[k]); minT = Math.min(minT, r.Te[k]);
    maxM = Math.max(maxM, Math.abs(r.moment[k]));
  }
  let tdp = null;
  if (cfg.seabed) {
    for (let k = 1; k <= n; k++) {
      if (r.reaction[k] > 0) {
        // touchdown located between the last free node and the first reacting node
        tdp = { s: 0.5 * (r.s[k - 1] + r.s[k]), x: 0.5 * (r.x[k - 1] + r.x[k]), node: k };
        break;
      }
    }
  }
  // sag (local minimum) and hog (local maximum) bends of wave configurations, from the top downwards
  let sag = null, hog = null;
  const kEnd = tdp ? tdp.node - 1 : n - 1;
  for (let k = 1; k < kEnd; k++) {
    if (!sag && r.z[k] < r.z[k - 1] && r.z[k] <= r.z[k + 1]) sag = { s: r.s[k], x: r.x[k], z: r.z[k] };
    else if (sag && !hog && r.z[k] > r.z[k - 1] && r.z[k] >= r.z[k + 1]) hog = { s: r.s[k], x: r.x[k], z: r.z[k] };
  }
  if (!hog) sag = null;
  return {
    ...r,
    topTension: r.Te[0], bottomTension: r.Te[n], maxTension: maxT, minTension: minT,
    hangoffAngleDeg: (r.angle[0] + Math.PI / 2) / DEG, bottomAngleDeg: (r.angle[n] + Math.PI / 2) / DEG,
    maxCurvature: maxK, sMaxCurvature: sMaxK, minBendRadius: maxK > 0 ? 1 / maxK : Infinity, maxMoment: maxM,
    tdp, sag, hog, config: cfg,
  };
}

/**
 * Static configuration of a riser, dynamic cable or umbilical between a hang-off point and a seabed end.
 *
 * Method: rod.js corotational beam model started from a lumped catenary; sections are listed from the top.
 * Reference: see solveRodStatic. Validity: in-plane, linear elastic sections, linear seabed springs.
 *
 * @param {object} cfg
 * @param {Array<object>} cfg.sections from top to seabed end, fields as in solveRodStatic (length, EA, EI, w, …)
 * @param {{x: number, z: number}} cfg.top hang-off position [m]
 * @param {{x: number, z: number}} cfg.bottom seabed-end position [m]
 * @param {object} [cfg.topCondition] `{ type: 'pinned' | 'clamped' | 'flex', angleDeg [deg from the downward
 *   vertical, + towards +x], kRot [N·m/rad] }` (default pinned)
 * @param {object} [cfg.bottomCondition] same form (e.g. steep-wave riser base clamped at angleDeg = 0)
 * @param {object} [cfg.seabed] `{ z, kn, mu?, kt?, profile? }` as in solveRodStatic
 * @param {object} [cfg.current] @param {object} [cfg.pressure] @param {Array<object>} [cfg.pointLoads]
 *   as in solveRodStatic (pointLoads[].s measured from the top; tethers and clumps of pliant-wave risers)
 * @param {number} [cfg.nElem=200] @param {object} [cfg.initial] previous result on the same mesh (warm start)
 * @returns {object} the solveRodStatic result plus `topTension`, `bottomTension`, `maxTension`, `minTension` [N],
 *   `hangoffAngleDeg`, `bottomAngleDeg` [deg], `maxCurvature` [1/m], `sMaxCurvature` [m], `minBendRadius` [m],
 *   `maxMoment` [N·m], `tdp` `{ s, x [m], node }` or null (resolution: one element), `sag`, `hog` `{ s, x, z [m] }`
 *   or null, and `config` (the input, for offset studies)
 */
export function solveRiser(cfg) {
  if (!cfg || !cfg.top || !cfg.bottom) throw new RangeError('riser: arguments top and bottom are required');
  num(cfg.top.x, 'top.x'); num(cfg.top.z, 'top.z'); num(cfg.bottom.x, 'bottom.x'); num(cfg.bottom.z, 'bottom.z');
  const model = {
    sections: cfg.sections, nElem: cfg.nElem ?? 200,
    endA: endFromCondition(cfg.top, cfg.topCondition), endB: endFromCondition(cfg.bottom, cfg.bottomCondition),
    seabed: cfg.seabed, current: cfg.current, pressure: cfg.pressure, pointLoads: cfg.pointLoads, supports: cfg.supports,
    initial: cfg.initial ?? 'cable', nSteps: cfg.nSteps, tol: cfg.tol,
  };
  return summarise(solveRodStatic(model), cfg);
}

/**
 * Steel catenary riser (or any free-hanging catenary line) with bending stiffness, flex/stress joint at the top
 * and a touchdown zone on linear seabed springs (optional trench profile).
 *
 * The configuration is sized from the inextensible catenary (catenaryRiser) for the requested hang-off angle or
 * top tension, then solved with the rod model; the mesh is refined at the top and through the touchdown zone.
 * For a tapered stress joint pass `topSections` (stiffer sections listed from the top) with a clamped top.
 * Reference: API RP 2RD §5 (global analysis); Sparks (2007) Ch. 13. Validity: as solveRiser.
 *
 * @param {object} o
 * @param {object} o.pipe section fields `{ EA, EI, w, D?, Cdn?, Cdt?, Ai?, Ae?, m?, ma?, contactOffset? }`
 * @param {number} o.h hang-off height above the seabed [m]
 * @param {number} [o.hangoffAngleDeg] nominal catenary top angle from vertical [deg] (or give topTension / H)
 * @param {number} [o.topTension] [N] @param {number} [o.H] horizontal tension [N]
 * @param {number} [o.groundedLength] pipe length kept on the seabed beyond the nominal touchdown [m]
 *   (default 30 % of the suspended length)
 * @param {object} [o.seabed] `{ kn [Pa] = 1e6, mu?, kt?, profile? }` (seabed at z = 0)
 * @param {object} [o.topCondition] as solveRiser (angleDeg defaults to the catenary angle)
 * @param {Array<object>} [o.topSections] extra sections at the top (stress joint), each with its own length
 * @param {number} [o.nElem=240] @param {object} [o.current] @param {object} [o.pressure]
 * @returns {object} solveRiser result plus `catenary` (the catenaryRiser estimate)
 */
export function solveScr(o) {
  const { pipe, h, groundedLength, nElem = 240 } = o ?? {};
  if (!pipe) throw new RangeError('riser: argument pipe is required');
  const cat = catenaryRiser({ h, w: pipe.w, hangoffAngleDeg: o.hangoffAngleDeg, topTension: o.topTension, H: o.H });
  const S = cat.suspendedLength;
  const Lg = num(groundedLength ?? 0.3 * S, 'groundedLength', { min: 0, strictMin: true });
  const lamT = Math.sqrt(pipe.EI / cat.H), lamTop = Math.sqrt(pipe.EI / cat.topTension);
  const tdz = Math.min(0.25 * S, Lg * 0.8, Math.max(15 * lamT, 0.08 * S));
  const topLen = Math.min(0.15 * S, Math.max(12 * lamTop, 0.03 * S));
  const extra = (o.topSections ?? []).reduce((a, s) => a + s.length, 0);
  const upLen = S - tdz - topLen - extra;
  if (!(upLen > 0)) throw new RangeError('riser: argument topSections is longer than the suspended length allows');
  const nTop = Math.max(8, Math.round(0.15 * nElem)), nTdz = Math.max(16, Math.round(0.45 * nElem));
  const nUp = Math.max(8, Math.round(0.28 * nElem)), nGr = Math.max(4, nElem - nTop - nTdz - nUp);
  const sections = [
    ...(o.topSections ?? []).map((s) => ({ ...pipe, ...s })),
    { ...pipe, length: topLen, nElem: nTop },
    { ...pipe, length: upLen, nElem: nUp },
    { ...pipe, length: 2 * tdz, nElem: nTdz },
    { ...pipe, length: Lg - tdz, nElem: nGr },
  ];
  const topCondition = o.topCondition ? { angleDeg: cat.hangoffAngleDeg, ...o.topCondition } : undefined;
  const res = solveRiser({
    sections, top: { x: 0, z: h }, bottom: { x: cat.horizontalSpan + Lg, z: 0 },
    seabed: { z: 0, kn: 1e6, ...(o.seabed ?? {}) }, topCondition, current: o.current, pressure: o.pressure,
  });
  res.catenary = cat;
  return res;
}

/**
 * Lazy-wave riser / dynamic cable / umbilical: upper catenary, buoyancy section and lower catenary to a
 * touchdown on the seabed. Steep-wave and pliant-wave variants use `bottomCondition` (riser base clamped
 * vertically, no grounded length) and `pointLoads` (tether or clump at the touchdown end of the buoyant
 * section) respectively.
 * Reference: API RP 17B §7 (configurations); Bai & Bai, "Subsea Pipelines and Risers", Ch. 23.
 * Validity: as solveRiser.
 *
 * @param {object} o
 * @param {object} o.pipe bare-section fields `{ EA, EI, w, … }`
 * @param {object} o.buoyant buoyancy-section fields; `w` is the net apparent weight [N/m], negative when buoyant
 * @param {number} o.upperLength @param {number} o.buoyancyLength @param {number} o.lowerLength [m] (lowerLength
 *   includes the grounded part)
 * @param {{x: number, z: number}} o.top @param {{x: number, z: number}} o.bottom [m]
 * @param {object} [o.seabed] default `{ z: bottom.z, kn: 1e6 }`; pass `null` for none (steep wave)
 * @param {number} [o.nElem=240] and the optional solveRiser fields
 * @returns {object} solveRiser result plus `sagClearance` (sag-bend height above the seabed, 0 when no wave has
 *   formed) and `hogElevation` (hog-bend height above the seabed) [m], `netBuoyancy` [N] of the buoyant section
 */
export function solveWaveRiser(o) {
  const { pipe, buoyant, upperLength, buoyancyLength, lowerLength, nElem = 240 } = o ?? {};
  if (!pipe || !buoyant) throw new RangeError('riser: arguments pipe and buoyant are required');
  num(upperLength, 'upperLength', { min: 0, strictMin: true }); num(buoyancyLength, 'buoyancyLength', { min: 0, strictMin: true });
  num(lowerLength, 'lowerLength', { min: 0, strictMin: true });
  const L = upperLength + buoyancyLength + lowerLength;
  const ne = (len) => Math.max(8, Math.round(nElem * len / L));
  const sections = [
    { ...pipe, length: upperLength, nElem: ne(upperLength) },
    { ...buoyant, length: buoyancyLength, nElem: Math.max(16, ne(buoyancyLength)) },
    { ...pipe, length: lowerLength, nElem: ne(lowerLength) },
  ];
  const seabed = o.seabed === null ? undefined : { z: o.bottom.z, kn: 1e6, ...(o.seabed ?? {}) };
  const res = solveRiser({ ...o, sections, seabed });
  const zs = o.bottom.z;
  res.sagClearance = res.sag ? res.sag.z - zs : 0;
  res.hogElevation = res.hog ? res.hog.z - zs : 0;
  res.netBuoyancy = -buoyant.w * buoyancyLength;
  return res;
}

/**
 * Buoyancy-section length of a lazy-wave configuration that gives a target sag-bend clearance above the seabed,
 * keeping the upper length and the total length fixed (the lower length absorbs the change).
 * Method: Brent root search on solveWaveRiser. Validity: the bracket must contain a sign change of
 * (clearance − target).
 *
 * @param {object} o solveWaveRiser arguments (buoyancyLength is ignored) plus
 * @param {number} o.targetSagClearance [m] @param {number} o.minLength @param {number} o.maxLength bracket [m]
 * @param {number} [o.tol=0.01] tolerance on the length [m]
 * @returns {{ buoyancyLength: number, lowerLength: number, result: object }} lengths [m] and the converged
 *   solveWaveRiser result
 */
export function solveBuoyancyLength(o) {
  const { targetSagClearance, minLength, maxLength, tol = 0.01 } = o ?? {};
  num(targetSagClearance, 'targetSagClearance', { min: 0 });
  num(minLength, 'minLength', { min: 0, strictMin: true }); num(maxLength, 'maxLength', { min: minLength, strictMin: true });
  const tail = o.buoyancyLength + o.lowerLength;
  if (!(maxLength < tail)) throw new RangeError('riser: argument maxLength must be shorter than buoyancyLength + lowerLength');
  const run = (Lb) => solveWaveRiser({ ...o, buoyancyLength: Lb, lowerLength: tail - Lb });
  const f = (Lb) => run(Lb).sagClearance - targetSagClearance;
  const Lb = brent(f, minLength, maxLength, { tol: tol / maxLength, maxIter: 60 });
  return { buoyancyLength: Lb, lowerLength: tail - Lb, result: run(Lb) };
}

/**
 * Vessel-offset cases of a solved riser: the hang-off is moved in-plane and the configuration re-solved from
 * the previous shape. Positive offset = "far" (away from the seabed end), negative = "near". Cross offsets are
 * out of plane and not modelled.
 * Reference: API RP 2RD §5.4 (near/far/cross load cases). Validity: as solveRiser.
 *
 * @param {object} base result of solveRiser / solveScr / solveWaveRiser
 * @param {number[]} offsets horizontal hang-off offsets [m]
 * @param {object} [opts] `{ heave?: number [m, vertical hang-off shift applied to every case] }`
 * @returns {Array<{ offset: number, topTension: number, bottomTension: number, hangoffAngleDeg: number,
 *   maxCurvature: number, tdpX: number, tdpS: number, tdpExcursion: number, minTension: number, result: object }>}
 *   offset, tdp positions and excursion relative to the base case [m]; tensions [N]; angle [deg]; curvature [1/m]
 */
export function riserOffsetCases(base, offsets, { heave = 0 } = {}) {
  if (!base || !base.config) throw new RangeError('riser: argument base must be a solveRiser result');
  const cfg = base.config;
  const away = cfg.bottom.x >= cfg.top.x ? -1 : 1;
  const sorted = offsets.map((v, i) => ({ v: num(v, `offsets[${i}]`), i }));
  const out = new Array(offsets.length);
  // march outwards from zero so that each solve starts from a nearby shape
  for (const sign of [1, -1]) {
    let prev = base;
    for (const o of sorted.filter((q) => (sign > 0 ? q.v >= 0 : q.v < 0)).sort((p, q) => sign * (p.v - q.v))) {
      const c = { ...cfg, top: { x: cfg.top.x + away * o.v, z: cfg.top.z + heave }, initial: prev };
      const r = solveRiser(c);
      r.config = { ...c, initial: undefined };
      out[o.i] = {
        offset: o.v, topTension: r.topTension, bottomTension: r.bottomTension, hangoffAngleDeg: r.hangoffAngleDeg,
        maxCurvature: r.maxCurvature, minTension: r.minTension, tdpX: r.tdp ? r.tdp.x : NaN, tdpS: r.tdp ? r.tdp.s : NaN,
        tdpExcursion: r.tdp && base.tdp ? r.tdp.x - base.tdp.x : NaN, result: r,
      };
      prev = r;
    }
  }
  return out;
}

/**
 * Quasi-static estimate of the dynamic response to hang-off heave: the riser is re-solved statically with the
 * hang-off raised and lowered by the heave amplitude (= RAO × wave amplitude). Inertia and drag are ignored, so
 * the result is a lower bound on dynamic tension for fast motions (use rodHarmonicResponse for the
 * frequency-domain estimate).
 * Reference: API RP 2RD §5.5 (quasi-static methods). Validity: motion periods well above the riser natural periods.
 *
 * @param {object} base solved riser @param {object} o `{ heaveAmp? [m] | raoHeave [m/m] and waveAmp [m] }`
 * @returns {{ heaveAmp: number, topTensionMax: number, topTensionMin: number, tensionAmplitude: number,
 *   tdpRange: number, maxCurvature: number, minTension: number }} [m], [N], [N], [N], [m], [1/m], [N]
 */
export function quasiStaticHeave(base, o) {
  const amp = o.heaveAmp ?? num(o.raoHeave, 'raoHeave', { min: 0 }) * num(o.waveAmp, 'waveAmp', { min: 0 });
  num(amp, 'heaveAmp', { min: 0 });
  const up = riserOffsetCases(base, [0], { heave: amp })[0], dn = riserOffsetCases(base, [0], { heave: -amp })[0];
  return {
    heaveAmp: amp, topTensionMax: Math.max(up.topTension, dn.topTension), topTensionMin: Math.min(up.topTension, dn.topTension),
    tensionAmplitude: 0.5 * Math.abs(up.topTension - dn.topTension), tdpRange: Math.abs(up.tdpX - dn.tdpX),
    maxCurvature: Math.max(up.maxCurvature, dn.maxCurvature), minTension: Math.min(up.minTension, dn.minTension),
  };
}

/* ------------------------------------------------------------------------------------------------------- */
/* Top-tensioned, drilling and free-standing risers                                                          */

/**
 * Top-tensioned (production or drilling) riser as a vertical tensioned beam: lower end at the wellhead with a
 * flex/ball joint, upper end offset with the vessel and held by a tensioner; current drag along the length.
 *
 * Method: rod.js beam with apparent weight; tensioner modelled as constant tension plus optional linear
 * stiffness on stroke (fixed-point iteration on the top set-down). Reference: API RP 16Q §3 (riser analysis),
 * API RP 2RD §5; Sparks (2007) Ch. 7–8. Validity: in-plane, quasi-static, top tension high enough to keep the
 * whole string in positive effective tension.
 *
 * @param {object} o
 * @param {number} o.length riser length wellhead to tensioner ring [m]
 * @param {object} o.pipe section fields `{ EA, EI, w, D?, Cdn?, … }`
 * @param {number} [o.topTension] nominal top tension [N] or @param {number} [o.tensionFactor] top tension /
 *   total apparent weight [-]
 * @param {number} [o.offset=0] vessel offset [m] @param {object} [o.current] as solveRodStatic (z from the wellhead)
 * @param {number} [o.kFlexLower] lower flex-joint stiffness [N·m/rad]; omit for a free ball joint
 * @param {number} [o.kFlexUpper] upper flex-joint stiffness [N·m/rad]; omit for a free joint
 * @param {number} [o.tensionerStiffness=0] [N/m] tension increase per metre of stroke-out
 * @param {number} [o.nElem=100] @param {object} [o.pressure]
 * @returns {{ topTension: number, bottomTension: number, tensionFactor: number, lowerAngleDeg: number,
 *   upperAngleDeg: number, stroke: number, maxMoment: number, maxCurvature: number, maxDeflection: number,
 *   lowerJointMoment: number, upperJointMoment: number, rod: object }} tensions [N]; angles from vertical [deg,
 *   + leaning towards +x]; stroke [m, top set-down relative to the zero-offset, zero-current position];
 *   moments [N·m]; curvature [1/m]; maxDeflection [m] from the straight chord
 */
export function solveTtr(o) {
  const { length, pipe, offset = 0, tensionerStiffness = 0, nElem = 100 } = o ?? {};
  num(length, 'length', { min: 0, strictMin: true });
  if (!pipe) throw new RangeError('riser: argument pipe is required');
  const W = pipe.w * length;
  const T0 = o.topTension ?? num(o.tensionFactor, 'tensionFactor', { min: 0, strictMin: true }) * W;
  num(T0, 'topTension', { min: 0, strictMin: true });
  num(tensionerStiffness, 'tensionerStiffness', { min: 0 });
  const up = Math.PI / 2;
  const run = (T, off, current) => solveRodStatic({
    sections: [{ ...pipe, length }], nElem,
    endA: { x: 0, z: 0, ...(o.kFlexLower !== undefined ? { angle: up, kRot: o.kFlexLower } : {}) },
    endB: { x: off, Fz: T, ...(o.kFlexUpper !== undefined ? { angle: up, kRot: o.kFlexUpper } : {}) },
    current, pressure: o.pressure, initialAngle: up, nSteps: 6,
  });
  let T = T0, r, stroke = 0;
  if (tensionerStiffness > 0) {
    const zRef = run(T0, 0, undefined).z[nElem];
    for (let it = 0; it < 60; it++) {
      r = run(T, offset, o.current);
      const sNew = zRef - r.z[nElem];
      const Tn = T0 + tensionerStiffness * sNew;
      const done = Math.abs(Tn - T) <= 1e-8 * T0;
      stroke = sNew; T = 0.5 * (T + Tn);
      if (done) break;
    }
  } else {
    r = run(T0, offset, o.current);
    stroke = run(T0, 0, undefined).z[nElem] - r.z[nElem];
  }
  const n = nElem;
  let maxM = 0, maxK = 0, maxD = 0;
  for (let k = 0; k <= n; k++) {
    maxM = Math.max(maxM, Math.abs(r.moment[k])); maxK = Math.max(maxK, Math.abs(r.curvature[k]));
    maxD = Math.max(maxD, Math.abs(r.x[k] - offset * r.z[k] / r.z[n]));
  }
  return {
    topTension: r.Te[n], bottomTension: r.Te[0], tensionFactor: r.Te[n] / W,
    lowerAngleDeg: (up - r.angle[0]) / DEG, upperAngleDeg: (up - r.angle[n]) / DEG, stroke,
    maxMoment: maxM, maxCurvature: maxK, maxDeflection: maxD,
    lowerJointMoment: Math.abs(r.moment[0]), upperJointMoment: Math.abs(r.moment[n]), rod: r,
  };
}

/**
 * Drilling-riser quasi-static offset sweep: lower/upper flex-joint angles against operating limits and the
 * offsets at which the limits are reached (watch-circle radii).
 * Reference: API RP 16Q (marine drilling riser operating guidelines: lower flex-joint mean and maximum angle
 * limits — the commonly applied 2° mean / 4° maximum while drilling are the argument defaults and must be
 * confirmed against the governing edition and the equipment ratings). Validity: as solveTtr.
 *
 * @param {object} o solveTtr arguments plus
 * @param {number[]} o.offsets ascending vessel offsets [m]
 * @param {number} [o.meanLimitDeg=2] @param {number} [o.maxLimitDeg=4] lower flex-joint angle limits [deg]
 * @returns {{ rows: Array<{ offset: number, lowerAngleDeg: number, upperAngleDeg: number, topTension: number,
 *   stroke: number, ucMean: number, ucMax: number }>, offsetAtMeanLimit: number, offsetAtMaxLimit: number }}
 *   offsets [m] (NaN when the limit is not reached inside the sweep), angles [deg], tension [N], stroke [m]
 */
export function drillingRiserOffsets(o) {
  const { offsets, meanLimitDeg = 2, maxLimitDeg = 4 } = o ?? {};
  if (!Array.isArray(offsets) || offsets.length < 2) throw new RangeError('riser: argument offsets needs at least two values');
  num(meanLimitDeg, 'meanLimitDeg', { min: 0, strictMin: true }); num(maxLimitDeg, 'maxLimitDeg', { min: 0, strictMin: true });
  const rows = offsets.map((off) => {
    const r = solveTtr({ ...o, offset: off });
    return { offset: off, lowerAngleDeg: r.lowerAngleDeg, upperAngleDeg: r.upperAngleDeg, topTension: r.topTension, stroke: r.stroke,
      ucMean: Math.abs(r.lowerAngleDeg) / meanLimitDeg, ucMax: Math.abs(r.lowerAngleDeg) / maxLimitDeg };
  });
  const cross = (lim) => {
    for (let k = 1; k < rows.length; k++) {
      const a0 = Math.abs(rows[k - 1].lowerAngleDeg), a1 = Math.abs(rows[k].lowerAngleDeg);
      if ((a0 - lim) * (a1 - lim) <= 0 && a1 !== a0) return rows[k - 1].offset + (lim - a0) * (rows[k].offset - rows[k - 1].offset) / (a1 - a0);
    }
    return NaN;
  };
  return { rows, offsetAtMeanLimit: cross(meanLimitDeg), offsetAtMaxLimit: cross(maxLimitDeg) };
}

/**
 * Drift-off of a vessel after loss of position under a constant net environmental force with quadratic
 * hydrodynamic resistance, and the watch-circle (alarm) offsets that leave time for the emergency disconnect
 * sequence before the physical limit is reached.
 *
 * Closed form of M·ẍ = F − c·ẋ|ẋ| from rest: v(t) = v∞·tanh(t/τ), x(t) = x0 + v∞·τ·ln cosh(t/τ) with v∞ = sqrt(F/c),
 * τ = M/sqrt(F·c) (first principles). Red alarm offset = position reached EDS-time before the limit; yellow =
 * a further reaction time earlier. Validity: constant force and resistance coefficient, single degree of freedom.
 *
 * @param {object} o
 * @param {number} o.mass vessel mass including surge/sway added mass [kg]
 * @param {number} o.force net mean environmental force [N]
 * @param {number} o.dragCoefficient resistance c in F = c·v² [N·s²/m²]
 * @param {number} o.limitOffset offset at which the riser/joint limit is reached [m]
 * @param {number} [o.x0=0] initial offset [m] @param {number} [o.edsTime=60] disconnect sequence duration [s]
 * @param {number} [o.reactionTime=60] time between yellow and red alarms [s]
 * @param {number} [o.nPoints=50]
 * @returns {{ terminalVelocity: number, timeConstant: number, timeToLimit: number, redOffset: number,
 *   yellowOffset: number, t: number[], x: number[], v: number[] }} [m/s], [s], [s], [m], [m], curve [s, m, m/s]
 */
export function driftOff(o) {
  const { mass, force, dragCoefficient, limitOffset, x0 = 0, edsTime = 60, reactionTime = 60, nPoints = 50 } = o ?? {};
  num(mass, 'mass', { min: 0, strictMin: true }); num(force, 'force', { min: 0, strictMin: true });
  num(dragCoefficient, 'dragCoefficient', { min: 0, strictMin: true }); num(limitOffset, 'limitOffset', { min: x0, strictMin: true });
  num(edsTime, 'edsTime', { min: 0 }); num(reactionTime, 'reactionTime', { min: 0 });
  const vInf = Math.sqrt(force / dragCoefficient), tau = mass / Math.sqrt(force * dragCoefficient);
  const lncosh = (y) => y + Math.log1p(Math.exp(-2 * y)) - Math.LN2;
  const xOf = (t) => (t <= 0 ? x0 : x0 + vInf * tau * lncosh(t / tau));
  const tLim = tau * Math.acosh(Math.exp((limitOffset - x0) / (vInf * tau)));
  const t = [], x = [], v = [];
  for (let k = 0; k < nPoints; k++) { const tk = tLim * k / (nPoints - 1); t.push(tk); x.push(xOf(tk)); v.push(vInf * Math.tanh(tk / tau)); }
  return { terminalVelocity: vInf, timeConstant: tau, timeToLimit: tLim, redOffset: xOf(tLim - edsTime), yellowOffset: xOf(tLim - edsTime - reactionTime), t, x, v };
}

function catenaryTwoPoint(L, X, h, w) {
  // inextensible catenary through two points: horizontal tension and end vertical forces (Irvine §1.2)
  if (L <= Math.hypot(X, h)) throw new RangeError('riser: jumper length must exceed the distance between its ends');
  const rhs = Math.sqrt(L * L - h * h) / X;
  const u = brent((q) => Math.sinh(q) / q - rhs, 1e-6, 60, { tol: 1e-13 });
  const a = X / (2 * u), H = w * a;
  return { H, V1: 0.5 * w * (L - h / Math.tanh(u)), V2: 0.5 * w * (L + h / Math.tanh(u)) };
}

/**
 * Free-standing hybrid riser (riser tower with buoyancy can) and flexible jumper: documented simple model.
 *
 * Model: the tower is a vertical tensioned beam pinned (or on a flex joint) at the seabed, tensioned at the top
 * by the net uplift of the buoyancy can; the jumper is an inextensible catenary between the tower top and the
 * vessel, applying its horizontal tension and its tower-end weight share to the tower top. Tower-top offset and
 * jumper span are iterated to consistency. Current acts on the tower (and optionally as a point drag on the can).
 * Reference: first-principles statics (tensioned beam + catenary; Irvine §1.2); configuration per API RP 2RD
 * Annex (hybrid risers). Validity: in-plane, jumper bending stiffness and current load on the jumper neglected,
 * can treated as a point load.
 *
 * @param {object} o
 * @param {number} o.length tower length seabed to can connection [m] @param {object} o.pipe section fields
 * @param {number} o.canNetUplift net buoyancy of the can (buoyancy − weight) [N]
 * @param {object} o.jumper `{ length [m], w [N/m apparent], vesselX [m], vesselZ [m] }` vessel end position relative
 *   to the tower base
 * @param {number} [o.kFlexBase] base flex-joint stiffness [N·m/rad] @param {object} [o.current]
 * @param {number} [o.canDragForce=0] horizontal current force on the can [N] @param {number} [o.nElem=80]
 * @returns {{ topOffset: number, topSetDown: number, baseAngleDeg: number, topTension: number, baseTension: number,
 *   jumperH: number, jumperTowerV: number, jumperVesselV: number, jumperVesselTension: number, maxMoment: number,
 *   iterations: number, rod: object }} offsets [m]; angle from vertical [deg]; forces [N]; moment [N·m]
 */
export function solveFreeStandingRiser(o) {
  const { length, pipe, canNetUplift, jumper, canDragForce = 0, nElem = 80 } = o ?? {};
  num(length, 'length', { min: 0, strictMin: true }); num(canNetUplift, 'canNetUplift', { min: 0, strictMin: true });
  if (!pipe || !jumper) throw new RangeError('riser: arguments pipe and jumper are required');
  num(jumper.length, 'jumper.length', { min: 0, strictMin: true }); num(jumper.w, 'jumper.w', { min: 0, strictMin: true });
  const up = Math.PI / 2;
  let xt = 0, zt = length, r, jc, it = 0;
  for (; it < 80; it++) {
    const X = jumper.vesselX - xt;
    jc = catenaryTwoPoint(jumper.length, Math.abs(X), jumper.vesselZ - zt, jumper.w);
    r = solveRodStatic({
      sections: [{ ...pipe, length }], nElem, initialAngle: up, nSteps: 6,
      endA: { x: 0, z: 0, ...(o.kFlexBase !== undefined ? { angle: up, kRot: o.kFlexBase } : {}) },
      endB: { Fx: Math.sign(X) * jc.H + canDragForce, Fz: canNetUplift - jc.V1 },
      current: o.current,
    });
    const dx = r.x[nElem] - xt, dz = r.z[nElem] - zt;
    xt += 0.7 * dx; zt += 0.7 * dz;
    if (Math.hypot(dx, dz) < 1e-7 * length) break;
  }
  let maxM = 0;
  for (const m of r.moment) maxM = Math.max(maxM, Math.abs(m));
  return {
    topOffset: r.x[nElem], topSetDown: length - r.z[nElem], baseAngleDeg: (up - r.angle[0]) / DEG,
    topTension: r.Te[nElem], baseTension: r.Te[0], jumperH: jc.H, jumperTowerV: jc.V1, jumperVesselV: jc.V2,
    jumperVesselTension: Math.hypot(jc.H, jc.V2), maxMoment: maxM, iterations: it + 1, rod: r,
  };
}

/**
 * Wall stress along a solved rod/riser: von Mises stress at each node from effective tension, bending moment and
 * hydrostatic internal/external pressure (wallStresses at every node).
 * @param {object} result rod or riser result (s, z, Te, moment)
 * @param {object} o `{ D [m], t [m], pIntRef? [Pa at zIntRef], zIntRef? [m], rhoInt? [kg/m³], rhoExt? [kg/m³],
 *   zSurface? [m], g? }`
 * @returns {{ s: Float64Array, vonMises: Float64Array, axial: Float64Array, hoop: Float64Array, max: number,
 *   sMax: number }} s [m]; stresses [Pa] (hoop at the inner fibre); max von Mises [Pa] at sMax [m]
 */
export function stressAlong(result, { D, t, pIntRef = 0, zIntRef = 0, rhoInt = 0, rhoExt = 1025, zSurface = 0, g = G0 }) {
  const n = result.s.length;
  const vonMises = new Float64Array(n), axial = new Float64Array(n), hoop = new Float64Array(n);
  let max = 0, sMax = 0;
  for (let k = 0; k < n; k++) {
    const pi = Math.max(0, pIntRef + rhoInt * g * (zIntRef - result.z[k]));
    const pe = Math.max(0, rhoExt * g * (zSurface - result.z[k]));
    const ws = wallStresses({ Te: result.Te[k], M: result.moment[k], pi, pe, D, t });
    vonMises[k] = ws.vonMises; axial[k] = ws.axial; hoop[k] = ws.hoopInner;
    if (ws.vonMises > max) { max = ws.vonMises; sMax = result.s[k]; }
  }
  return { s: Float64Array.from(result.s), vonMises, axial, hoop, max, sMax };
}
