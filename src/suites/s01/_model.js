// Shared jacket model of the fixed-platform suite.
//  * geometry: parametric jacket (builders.jacketModel, X braces noded at their crossing) or an imported frame;
//  * foundation: fixed, pinned, linear 6×6 pile-head springs (given, or computed with the p–y / t–z pile engine);
//  * environmental loading: Morison force on every wetted member from a regular wave (Airy, Wheeler-stretched
//    Airy, Stokes 5th or stream function) plus a power-law current, stepped through the wave phase, and wind and
//    interface loads on the topsides;
//  * post-processing: member forces, tubular code checks, reactions and model3d views.
// Model axes: z up, z = 0 at the still-water level (before surge), mudline at z = −depth.
import { jacketModel, fromImportedFrame } from '../../engine/structural/builders.js';
import { memberAxes, analyzeFrame, memberDiagram, sectionStress, refineModel, buildFrameSystem, factorFrame, solveFrame, resolveLoads, loadParts, elementMatVec } from '../../engine/structural/frame3d.js';
import { frameGeometricStiffness } from '../../engine/structural/framebuckling.js';
import { eigGenSym } from '../../engine/math/eigen.js';
import { resolveSection } from '../../engine/structural/sections.js';
import { tubularMemberCheck } from '../../engine/structural/tubular.js';
import { airyKinematics, stokes5Wave, stokes5Kinematics, selectWaveTheory, wavenumber } from '../../engine/metocean/waves.js';
import { streamFunctionWave, streamFunctionKinematics, streamFunctionEta } from '../../engine/metocean/streamfunction.js';
import { currentPowerLaw } from '../../engine/metocean/current.js';
import { windForce, heightCoefficient } from '../../engine/metocean/wind.js';
import { pileHeadStiffness, axialLoadSettlement } from '../../engine/geotech/pile.js';
import { createProfile } from '../../engine/geotech/soilprofile.js';
import { G, DEG, need, run, round, fmt, cap, CONFIRM } from './_common.js';

/* ------------------------------------------------------------------------------------------------ inputs */

export function geomInputs(o = {}) {
  const g = 'Jacket geometry', m = 'Member sizes';
  return [
    { key: 'geom', label: 'Imported frame model (optional)', type: 'geometry', kinds: ['frame'], default: null, group: g, role: 'geometry',
      help: 'SDNF, structured JSON / CSV or FE beam model. Its lowest nodes are taken as the pile heads at the mudline and its highest nodes as the topsides interface. When empty, the parametric jacket below is analysed.' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: o.depth ?? 80, min: 5, max: 400, group: g, role: 'environment', bind: 'site.depth' },
    { key: 'topEl', label: 'Jacket top elevation above still water', unit: 'm', type: 'number', default: 15, min: 2, max: 40, group: g, role: 'geometry', help: 'Work point where the topsides legs stab in.' },
    { key: 'topW', label: 'Leg spacing at the top', unit: 'm', type: 'number', default: 18, min: 4, max: 80, group: g, role: 'geometry', design: true },
    { key: 'batter', label: 'Leg batter (1 horizontal : n vertical)', unit: '-', type: 'number', default: 10, min: 3, max: 1000, group: g, role: 'geometry', help: 'Both faces battered; 1000 gives practically vertical legs.' },
    { key: 'bays', label: 'Bracing bays', type: 'int', default: 4, min: 1, max: 8, group: g, role: 'geometry' },
    { key: 'bracing', label: 'Face bracing pattern', type: 'select', default: 'X', group: g, role: 'geometry',
      options: [{ value: 'X', label: 'X bracing (noded at the crossing)' }, { value: 'K', label: 'K bracing' }, { value: 'diagonal', label: 'Single diagonal' }] },
    { key: 'legD', label: 'Leg diameter', unit: 'mm', type: 'number', default: 1500, min: 200, max: 6000, group: m, role: 'geometry', design: true },
    { key: 'legT', label: 'Leg wall thickness', unit: 'mm', type: 'number', default: 40, min: 6, max: 150, group: m, role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'brD', label: 'Diagonal brace diameter', unit: 'mm', type: 'number', default: 700, min: 100, max: 4000, group: m, role: 'geometry', design: true },
    { key: 'brT', label: 'Diagonal brace wall thickness', unit: 'mm', type: 'number', default: 20, min: 5, max: 120, group: m, role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'hzD', label: 'Horizontal brace diameter', unit: 'mm', type: 'number', default: 800, min: 100, max: 4000, group: m, role: 'geometry', design: true },
    { key: 'hzT', label: 'Horizontal brace wall thickness', unit: 'mm', type: 'number', default: 20, min: 5, max: 120, group: m, role: 'geometry', design: true },
    { key: 'mat', label: 'Steel', type: 'select', source: 'materials:steel', default: 's355', group: m, role: 'material' },
  ];
}

export function envInputs(o = {}) {
  const g = o.group ?? 'Environment', h = 'Hydrodynamic coefficients';
  const out = [
    { key: 'H', label: o.Hlabel ?? 'Design wave height', unit: 'm', type: 'number', default: o.H ?? 18, min: 0, max: 40, group: g, role: 'bc', bind: o.Hbind === null ? undefined : (o.Hbind ?? 'site.Hmax100'), uncertain: { dist: 'lognormal', cov: 0.1 },
      help: 'Individual (maximum) wave height of the design sea state, crest to trough.' },
    { key: 'T', label: 'Wave period', unit: 's', type: 'number', default: o.T ?? 13.5, min: 3, max: 25, group: g, role: 'bc', bind: o.Tbind ?? 'site.Tp100', uncertain: { dist: 'normal', cov: 0.07 } },
    { key: 'Uc', label: 'Surface current speed', unit: 'm/s', type: 'number', default: o.Uc ?? 1.0, min: 0, max: 4, group: g, role: 'bc', bind: o.Ucbind ?? 'site.Uc100', uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'In line with the wave; 1/7 power-law profile to the seabed.' },
    { key: 'Uw', label: '1-hour mean wind speed at 10 m', unit: 'm/s', type: 'number', default: o.Uw ?? 35, min: 0, max: 80, group: g, role: 'bc', bind: o.Uwbind ?? 'site.Uw100', uncertain: { dist: 'lognormal', cov: 0.1 } },
    { key: 'surge', label: 'Storm surge and tide above still water', unit: 'm', type: 'number', default: o.surge ?? 1.0, min: 0, max: 8, group: g, role: 'bc', bind: 'site.surge100' },
  ];
  if (o.heading !== false) out.push({ key: 'heading', label: 'Wave, current and wind heading', unit: 'deg', type: 'number', default: o.heading ?? 0, min: 0, max: 360, group: g, role: 'bc', help: 'Direction of travel measured from +x towards +y: 0 is broadside on, 45 is the diagonal of a square jacket.' });
  out.push(
    { key: 'theory', label: 'Wave theory', type: 'select', default: o.theory ?? 'auto', group: g, role: 'control',
      options: [{ value: 'auto', label: 'Automatic (applicability chart)' }, { value: 'airy', label: 'Airy, loaded to still water' }, { value: 'wheeler', label: 'Airy with Wheeler stretching' }, { value: 'stokes5', label: 'Stokes 5th order' }, { value: 'stream', label: 'Stream function' }],
      help: 'Automatic follows the approximated applicability chart of API RP 2A / ISO 19901-1 and falls back on the stream function near breaking.' },
    { key: 'Cd', label: 'Drag coefficient Cd', unit: '-', type: 'number', default: 1.05, min: 0, max: 2.5, group: h, role: 'criteria', uncertain: { dist: 'normal', cov: 0.1 }, help: `0.65 smooth, 1.05 rough (marine growth) in API RP 2A. ${CONFIRM}` },
    { key: 'Cm', label: 'Inertia coefficient Cm', unit: '-', type: 'number', default: 1.2, min: 0, max: 3, group: h, role: 'criteria', help: '1.6 smooth, 1.2 rough in API RP 2A. Morison loading is valid for slender members only (D/λ < 0.2).' },
    { key: 'mg', label: 'Marine growth thickness', unit: 'mm', type: 'number', default: 50, min: 0, max: 300, group: h, role: 'environment', help: 'Added to the radius between the water surface and the depth below.' },
    { key: 'mgDepth', label: 'Marine growth down to depth', unit: 'm', type: 'number', default: 50, min: 0, max: 400, group: h, role: 'environment' },
    { key: 'kinFac', label: 'Wave kinematics factor', unit: '-', type: 'number', default: 0.9, min: 0.5, max: 1, group: h, role: 'criteria', help: `Reduction of horizontal particle velocity and acceleration for directional spreading (0.85–0.95 for storms). ${CONFIRM}` },
    { key: 'block', label: 'Current blockage factor', unit: '-', type: 'number', default: 0.8, min: 0.5, max: 1, group: h, role: 'criteria', help: `0.80 end-on / 0.85 diagonal / 0.80 broadside for four legs in API RP 2A. ${CONFIRM}` },
    { key: 'appF', label: 'Wave-load allowance for conductors and appurtenances', unit: '-', type: 'number', default: 1.5, min: 1, max: 4, group: h, role: 'bc', uncertain: { dist: 'normal', cov: 0.1 },
      help: 'Factor on the Morison load of the modelled members for conductors, risers, caissons, boat landings and anodes that are not in the frame model (typically 1.3–2.0). Use 1.0 when they are modelled as members.' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: h, role: 'environment', bind: 'site.rhoW' },
  );
  return out;
}

export function topsidesInputs(o = {}) {
  const g = 'Topsides and gravity loads';
  return [
    { key: 'Wtop', label: 'Topsides weight', unit: 'MN', type: 'number', default: o.Wtop ?? 40, min: 0, max: 2000, group: g, role: 'bc', bind: 'topsides.weight', uncertain: { dist: 'normal', cov: 0.05 }, help: 'Dead, equipment and variable deck load at the jacket interface, shared equally by the leg tops.' },
    { key: 'FxTop', label: 'Topsides horizontal interface load', unit: 'kN', type: 'number', default: 0, min: 0, max: 1e5, group: g, role: 'bc', bind: 'topsides.Fx', help: 'Published by the topsides suite (wind, crane, drilling). Applied along the wave heading. Set the wind area to zero when this already contains wind.' },
    { key: 'MyTop', label: 'Topsides overturning moment at the interface', unit: 'MN·m', type: 'number', default: 0, min: 0, max: 1e4, group: g, role: 'bc', bind: 'topsides.My', help: 'Applied about the horizontal axis normal to the wave heading as a push–pull on the leg tops.' },
    { key: 'windA', label: 'Topsides wind area', unit: 'm²', type: 'number', default: o.windA ?? 600, min: 0, max: 2e4, group: g, role: 'bc', help: 'Projected area normal to the wind.' },
    { key: 'windCs', label: 'Wind shape coefficient', unit: '-', type: 'number', default: 1.0, min: 0.3, max: 2, group: g, role: 'criteria', help: '1.0 for the overall projected area of a platform, 1.5 for beams and sides of buildings (API RP 2A).' },
    { key: 'windZ', label: 'Wind area centroid above still water', unit: 'm', type: 'number', default: 30, min: 5, max: 150, group: g, role: 'bc' },
    { key: 'swCont', label: 'Self-weight contingency', unit: '-', type: 'number', default: 1.1, min: 0, max: 1.5, group: g, role: 'bc', help: 'Factor on modelled steel weight for joint cans, anodes, appurtenances and weld metal; 0 switches self-weight off.' },
  ];
}

export function foundationInputs(o = {}) {
  const g = 'Foundation (boundary condition)', sp = { key: 'found', in: ['springs'] }, pl = { key: 'found', in: ['pile'] };
  return [
    { key: 'found', label: 'Pile-head condition', type: 'select', default: o.found ?? 'springs', group: g, role: 'bc',
      options: [{ value: 'fixed', label: 'Fixed at the mudline' }, { value: 'pinned', label: 'Pinned at the mudline' }, { value: 'springs', label: 'Linear pile-head springs (given)' }, { value: 'pile', label: 'Pile-head springs computed (p–y / t–z)' }],
      help: 'Springs are a linearisation of the pile and soil at the working load: 6×6 per leg with lateral–rocking coupling.' },
    { key: 'Kh', label: 'Lateral pile-head stiffness', unit: 'MN/m', type: 'number', default: 100, min: 1, max: 1e5, group: g, role: 'bc', bind: 'foundation.Kh', showIf: sp, uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'Kr', label: 'Rotational pile-head stiffness', unit: 'GN·m/rad', type: 'number', default: 6.8, min: 0.01, max: 1e4, group: g, role: 'bc', bind: 'foundation.Kr', showIf: sp },
    { key: 'Khr', label: 'Lateral–rotational coupling (magnitude)', unit: 'MN', type: 'number', default: 690, min: 0, max: 1e6, group: g, role: 'bc', bind: 'foundation.Khr', showIf: sp, help: 'Entered positive; the coupling term is negative in the pile sign convention. Kh·Kr must exceed Khr².' },
    { key: 'Kv', label: 'Axial pile-head stiffness', unit: 'MN/m', type: 'number', default: 2000, min: 1, max: 1e6, group: g, role: 'bc', bind: 'foundation.Kv', showIf: sp, uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'pileD', label: 'Pile diameter', unit: 'mm', type: 'number', default: 1829, min: 300, max: 5000, group: g, role: 'geometry', showIf: pl },
    { key: 'pileT', label: 'Pile wall thickness', unit: 'mm', type: 'number', default: 60, min: 10, max: 150, group: g, role: 'geometry', showIf: pl },
    { key: 'pileL', label: 'Pile penetration', unit: 'm', type: 'number', default: 70, min: 10, max: 200, group: g, role: 'geometry', showIf: pl },
    { key: 'soilType', label: 'Soil', type: 'select', default: 'clay', group: g, role: 'soil', showIf: pl, options: [{ value: 'clay', label: 'Clay (su profile)' }, { value: 'sand', label: 'Sand (φ′)' }] },
    { key: 'su0', label: 'Undrained shear strength at mudline', unit: 'kPa', type: 'number', default: 50, min: 1, max: 500, group: g, role: 'soil', bind: 'site.su0', showIf: pl },
    { key: 'suGrad', label: 'Strength gradient', unit: 'kPa/m', type: 'number', default: 4, min: 0, max: 50, group: g, role: 'soil', bind: 'site.suGrad', showIf: pl },
    { key: 'phiDeg', label: 'Sand friction angle φ′', unit: 'deg', type: 'number', default: 35, min: 28, max: 42, group: g, role: 'soil', bind: 'site.phiDeg', showIf: pl },
    { key: 'gammaSub', label: 'Submerged unit weight', unit: 'kN/m³', type: 'number', default: 8, min: 3, max: 13, group: g, role: 'soil', bind: 'site.gammaSub', showIf: pl },
  ];
}

export function criteriaInputs(o = {}) {
  const g = 'Design criteria';
  return [
    { key: 'code', label: 'Member check', type: 'select', default: 'ISO', group: g, role: 'criteria', options: [{ value: 'ISO', label: 'ISO 19902 (partial factors)' }, { value: 'API-LRFD', label: 'API RP 2A-LRFD' }, { value: 'API-WSD', label: 'API RP 2A-WSD (use action factors of 1.0)' }],
      help: 'Equations as recalled from the standards; confirm against the governing edition.' },
    { key: 'gG', label: 'Action factor — permanent and variable', unit: '-', type: 'number', default: o.gG ?? 1.1, min: 0.8, max: 2, group: g, role: 'criteria', help: `ISO 19902 extreme condition: 1.1 G + 1.35 E; operating: 1.3 G + 1.5 Q + 0.9·1.35 E. ${CONFIRM}` },
    { key: 'gE', label: 'Action factor — environmental', unit: '-', type: 'number', default: o.gE ?? 1.35, min: 0.5, max: 2.5, group: g, role: 'criteria', help: CONFIRM },
    { key: 'gRt', label: 'Resistance factor — tension γR,t', unit: '-', type: 'number', default: 1.05, min: 1, max: 1.5, group: g, role: 'criteria', help: `ISO 19902 only. ${CONFIRM}` },
    { key: 'gRc', label: 'Resistance factor — compression γR,c', unit: '-', type: 'number', default: 1.18, min: 1, max: 1.6, group: g, role: 'criteria', help: `ISO 19902 only: 1.18 in the 2007 edition, 1.10 in the 2020 edition — the default is the more conservative of the two. ${CONFIRM}` },
    { key: 'gRb', label: 'Resistance factor — bending γR,b', unit: '-', type: 'number', default: 1.05, min: 1, max: 1.5, group: g, role: 'criteria', help: `ISO 19902 only. ${CONFIRM}` },
    { key: 'Kleg', label: 'Effective length factor — legs', unit: '-', type: 'number', default: 1.0, min: 0.5, max: 2.5, group: g, role: 'criteria' },
    { key: 'Kbr', label: 'Effective length factor — braces', unit: '-', type: 'number', default: 0.8, min: 0.5, max: 2.5, group: g, role: 'criteria', help: '0.7–0.8 for braces in ISO 19902 Table 13.5-1 (node to node). X braces are noded at their crossing.' },
  ];
}

export function meshInputs(o = {}) {
  const g = 'Solver';
  const out = [{ key: 'nSub', label: 'Elements per member', type: 'int', default: o.nSub ?? 2, min: 1, max: 8, group: g, role: 'control', help: 'Sub-division of every member: refines the distributed wave load, member bending and mode shapes.' }];
  if (o.phase !== false) out.push({ key: 'nPhase', label: 'Wave phase steps per period', type: 'int', default: o.nPhase ?? 16, min: 6, max: 144, group: g, role: 'control', help: 'The wave is stepped past the structure in this many positions to find the largest base shear and overturning moment.' });
  return out;
}

/* ---------------------------------------------------------------------------------------------- foundation */

const kCache = new Map();
/** Pile-head stiffness (SI) for the chosen foundation option, or null for rigid supports. */
export function foundationK(inp, found = inp.found) {
  if (found === 'springs') {
    const K = { Kh: inp.Kh * 1e6, Kr: inp.Kr * 1e9, Khr: inp.Khr * 1e6, Kv: inp.Kv * 1e6, source: 'given' };
    need(K.Kh * K.Kr > K.Khr * K.Khr * 1.0001, 'The pile-head stiffness matrix is not positive definite: Kh·Kr must exceed Khr². Reduce the coupling term or use consistent values from the pile module.');
    return K;
  }
  if (found !== 'pile') return null;
  const key = [inp.pileD, inp.pileT, inp.pileL, inp.soilType, inp.su0, inp.suGrad, inp.phiDeg, inp.gammaSub].join('|');
  if (kCache.has(key)) return kCache.get(key);
  const D = inp.pileD / 1e3, t = inp.pileT / 1e3, L = inp.pileL;
  need(t < D / 2, 'Pile wall thickness must be smaller than the pile radius.');
  const thick = L + 20;
  const layer = inp.soilType === 'clay'
    ? { type: 'clay', top: 0, thickness: thick, gammaSub: inp.gammaSub * 1e3, suTop: inp.su0 * 1e3, suGrad: inp.suGrad * 1e3, eps50: inp.su0 + inp.suGrad * 10 > 100 ? 0.005 : 0.01 }
    : { type: 'sand', top: 0, thickness: thick, gammaSub: inp.gammaSub * 1e3, phiDeg: inp.phiDeg, Dr: Math.min(0.95, Math.max(0.36, 0.04 * (inp.phiDeg - 20))) };
  const profile = run(() => createProfile([layer]), 'Soil profile');
  // Secant stiffness at a representative working load (fixed head, as for a pile grouted in the leg).
  const hs = run(() => pileHeadStiffness({ L, D, t, profile, loading: 'cyclic', nElem: 40, head: { H: 1.0e6, fixedRotation: true } }), 'Pile-head stiffness (p–y)');
  const ax = run(() => axialLoadSettlement({ profile, L, D, t, nElem: 30 }, { steps: 30 }), 'Pile axial stiffness (t–z)');
  const K = { Kh: hs.secant.KHH, Kr: hs.secant.KMM, Khr: Math.abs(hs.secant.KHM), Kv: ax.headStiffness, source: 'computed', Qc: ax.capacity };
  need([K.Kh, K.Kr, K.Khr, K.Kv].every((v) => Number.isFinite(v) && v >= 0) && K.Kh > 0 && K.Kv > 0, 'The pile-head stiffness could not be computed for this pile and soil: check penetration and soil strength.');
  kCache.set(key, K);
  return K;
}

/** 6×6 support spring in global axes for a vertical pile. Torsional head stiffness is taken as 25 % of the rocking stiffness. */
export function springMatrix(K) {
  const k = Array.from({ length: 6 }, () => new Array(6).fill(0));
  k[0][0] = k[1][1] = K.Kh; k[2][2] = K.Kv; k[3][3] = k[4][4] = K.Kr; k[5][5] = 0.25 * K.Kr;
  // A head pushed towards +x leans towards +x (rotation +ry) when free: the ux–ry term is negative, uy–rx positive.
  k[0][4] = k[4][0] = -K.Khr; k[1][3] = k[3][1] = K.Khr;
  return k;
}

/* ------------------------------------------------------------------------------------------------ geometry */

// X braces generated by the builder cross without a node; real X joints are welded, so the crossing is noded.
function splitX(model) {
  const out = [], P = model.nodes;
  for (let i = 0; i < model.members.length; i++) {
    const m1 = model.members[i], m2 = model.members[i + 1];
    if (m1.group === 'brace' && m2 && m2.group === 'brace') {
      const a = P[m1.a], d1 = P[m1.b].map((v, k) => v - a[k]), c = P[m2.a], d2 = P[m2.b].map((v, k) => v - c[k]);
      const r = c.map((v, k) => v - a[k]);
      const a11 = d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2], a22 = d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2], a12 = -(d1[0] * d2[0] + d1[1] * d2[1] + d1[2] * d2[2]);
      const b1 = d1[0] * r[0] + d1[1] * r[1] + d1[2] * r[2], b2 = -(d2[0] * r[0] + d2[1] * r[1] + d2[2] * r[2]);
      const det = a11 * a22 - a12 * a12, s = (b1 * a22 - a12 * b2) / det, u = (a11 * b2 - a12 * b1) / det;
      if (det > 1e-9 && s > 0.05 && s < 0.95 && u > 0.05 && u < 0.95) {
        P.push([a[0] + s * d1[0], a[1] + s * d1[1], a[2] + s * d1[2]]);
        const x = P.length - 1;
        out.push({ ...m1, b: x }, { ...m1, a: x }, { ...m2, b: x }, { ...m2, a: x });
        i++; continue;
      }
    }
    out.push(m1);
  }
  model.members = out;
}

const GROUP_LABEL = { leg: 'Legs', brace: 'Diagonal braces', horizontal: 'Horizontals', plan: 'Plan braces', pile: 'Piles' };

/**
 * Build the analysis model.
 * @param {object} inp module inputs (engineering units)
 * @param {object} ctx module context
 * @param {{nSub?:number, found?:string, wallLoss?:Object<string,number>, remove?:number[], fy?:number, layDown?:boolean}} [o] wallLoss in m per group key;
 *   layDown turns the parametric jacket onto one face (transport position)
 * @returns {object} J: { model (refined frame3d model), base (one element per member), memberMap, orig[], elems[], topNodes, mudNodes,
 *   steel, steelT, zMud, zTop, centre, imported, K }
 */
export function buildJacket(inp, ctx, o = {}) {
  const mat = ctx.material(inp.mat);
  need(mat, `Steel '${inp.mat}' is not in the materials library.`);
  const steel = { E: mat.E, nu: mat.nu ?? 0.3, rho: mat.rho ?? 7850, alpha: mat.alpha ?? 1.2e-5, fy: o.fy ?? mat.fy, fu: mat.fu ?? 1.2 * mat.fy };
  const found = o.found ?? inp.found ?? 'fixed';
  let base, topNodes, mudNodes, imported = false;
  if (inp.geom) {
    need(inp.geom.kind === 'frame', 'The imported model must be a frame (nodes and members).');
    base = run(() => fromImportedFrame(inp.geom, { defaultSection: { type: 'tube', D: inp.brD / 1e3, t: inp.brT / 1e3 } }), 'Imported frame');
    const zs = base.nodes.map((n) => n[2]), z0 = Math.min(...zs), z1 = Math.max(...zs), tol = Math.max(1e-3, 1e-3 * (z1 - z0));
    // Datum shift: the lowest nodes are the mudline.
    base.nodes = base.nodes.map((n) => [n[0], n[1], n[2] - z0 - inp.depth]);
    mudNodes = base.supports.map((s) => s.node);
    const used = new Set(base.members.flatMap((m) => [m.a, m.b]));
    topNodes = base.nodes.map((n, i) => i).filter((i) => used.has(i) && zs[i] >= z1 - tol && !mudNodes.includes(i));
    if (!topNodes.length) topNodes = base.nodes.map((n, i) => i).filter((i) => used.has(i) && !mudNodes.includes(i));
    need(topNodes.length, 'The imported frame has no free node above its supports: nothing can be loaded.');
    base.members.forEach((m) => { m.group = m.group || 'member'; });
    imported = true;
  } else {
    for (const [d, t, n] of [[inp.legD, inp.legT, 'leg'], [inp.brD, inp.brT, 'diagonal brace'], [inp.hzD, inp.hzT, 'horizontal brace']]) need(t < d / 2, `Wall thickness of the ${n} must be smaller than its radius.`);
    need(inp.brD <= inp.legD && inp.hzD <= inp.legD, 'Brace diameters must not exceed the leg diameter (β ≤ 1 at the joints).');
    base = run(() => jacketModel({ waterDepth: inp.depth, topElevation: inp.topEl, topWidthX: inp.topW, batter: inp.batter, bays: inp.bays, bracing: inp.bracing,
      legSection: { D: inp.legD / 1e3, t: inp.legT / 1e3 }, braceSection: { D: inp.brD / 1e3, t: inp.brT / 1e3 }, horizontalSection: { D: inp.hzD / 1e3, t: inp.hzT / 1e3 }, piles: { mode: 'fixed' }, floodedLegs: true }), 'Jacket geometry');
    if (inp.bracing === 'X') splitX(base);
    topNodes = base.meta.topNodes.slice(); mudNodes = base.meta.mudlineNodes.slice();
    if (o.layDown) {
      // Transport position: the jacket lies on the face with the lowest x, rotated about y until that face is level.
      // New axes: x along the face from base to top (barge longitudinal), z up from the face; supports at every
      // framing level of the two legs of that face.
      const N = base.nodes, cx0 = mudNodes.reduce((q, n) => q + N[n][0], 0) / mudNodes.length, legs = mudNodes.map((n, k) => k).filter((k) => N[mudNodes[k]][0] < cx0);
      const p0 = N[mudNodes[legs[0]]], p1 = N[topNodes[legs[0]]], l1 = Math.hypot(p1[0] - p0[0], p1[2] - p0[2]), e1 = [(p1[0] - p0[0]) / l1, (p1[2] - p0[2]) / l1];
      base.nodes = N.map((q) => { const dx = q[0] - p0[0], dz = q[2] - p0[2]; return [dx * e1[0] + dz * e1[1], q[1], dx * e1[1] - dz * e1[0]]; });
      mudNodes = base.meta.legNodes.flatMap((lv) => legs.map((k) => lv[k]));
      topNodes = base.meta.legNodes[base.meta.legNodes.length - 1].slice();
    }
  }
  if (o.remove && o.remove.length) {
    const rm = new Set(o.remove);
    for (const i of rm) need(Number.isInteger(i) && i >= 0 && i < base.members.length, `Member index ${i} does not exist (the model has members 0 to ${base.members.length - 1}).`);
    base.members = base.members.filter((m, i) => !rm.has(i));
  }
  const K = foundationK(inp, found);
  base.supports = mudNodes.map((n) => (K ? { node: n, spring: springMatrix(K) } : { node: n, fix: found === 'pinned' ? 'pinned' : 'fixed' }));
  base.materials = { steel }; base.loadCases = []; delete base.nodeMasses;
  const P = base.nodes;
  const orig = base.members.map((m, i) => {
    let sdef = m.section;
    const loss = o.wallLoss && (o.wallLoss[m.group] || 0);
    if (loss && sdef.type === 'tube') { need(sdef.t - loss > 0.002, `The wall loss leaves less than 2 mm of wall in the ${m.group} members.`); sdef = { type: 'tube', D: sdef.D, t: sdef.t - loss }; m.section = sdef; }
    const sec = resolveSection(sdef), a = P[m.a], b = P[m.b], L = Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
    m.material = 'steel';
    const Kf = m.group === 'leg' ? (inp.Kleg ?? 1) : (inp.Kbr ?? 0.8);
    m.K = Kf; m.bucklingLength = L;   // read by the pushover engine; survives refinement
    return { i, a: m.a, b: m.b, sec, L, K: Kf, group: m.group, label: GROUP_LABEL[m.group] || m.group, flooded: !!m.flooded, tube: sec.type === 'tube',
      D: sec.D || sec.h || 2 * Math.sqrt(sec.A / Math.PI), zmid: 0.5 * (a[2] + b[2]) };
  });
  const nSub = o.nSub ?? inp.nSub ?? 1;
  const ref = nSub > 1 ? run(() => refineModel(base, nSub), 'Mesh refinement') : { model: base, memberMap: base.members.map((m, i) => [i]) };
  const model = ref.model, Q = model.nodes;
  const elems = [];
  ref.memberMap.forEach((list, i) => list.forEach((e) => {
    const m = model.members[e], a = Q[m.a], b = Q[m.b], L = Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
    elems[e] = { e, o: i, a, b, L, t: [(b[0] - a[0]) / L, (b[1] - a[1]) / L, (b[2] - a[2]) / L], D: orig[i].D };
  }));
  const steelT = orig.reduce((s, m) => s + m.sec.A * m.L * steel.rho, 0) / 1e3;
  const zMud = Math.min(...mudNodes.map((n) => P[n][2])), zTop = Math.max(...topNodes.map((n) => P[n][2]));
  const centre = [mudNodes.reduce((s, n) => s + P[n][0], 0) / mudNodes.length, mudNodes.reduce((s, n) => s + P[n][1], 0) / mudNodes.length];
  const size = Math.max(zTop - zMud, ...Q.map((n) => Math.max(Math.abs(n[0] - centre[0]), Math.abs(n[1] - centre[1]))));
  return { model, base, memberMap: ref.memberMap, orig, elems, topNodes, mudNodes, steel, steelT, zMud, zTop, centre, size, imported, K, found, nSub };
}

/* ---------------------------------------------------------------------------------------------- wave field */

const ZERO = { u: 0, w: 0, ax: 0, az: 0 };
/**
 * Regular-wave kinematics provider. kin(x, z, t) returns { u, w, ax, az } along the wave direction (z from the
 * storm water level), surf(x, t) the elevation used to truncate the loaded length.
 */
export function waveField(inp, o = {}) {
  const H = o.H ?? inp.H, T = o.T ?? inp.T, d = inp.depth + (inp.surge || 0);
  if (!(H > 1e-4)) return { H: 0, T, d, theory: 'none', label: 'no wave', kin: () => ZERO, surf: () => 0, crest: 0, L: Infinity, k: 0 };
  const k0 = run(() => wavenumber({ T, d }), 'Wave dispersion'), L0 = 2 * Math.PI / k0;
  let theory = o.theory ?? inp.theory ?? 'auto', sel = null, order = 10;
  if (theory === 'auto') {
    sel = selectWaveTheory({ H, T, d });
    need(sel.theory !== 'breaking', `A ${fmt(H)} m wave of ${fmt(T)} s period breaks in ${fmt(d)} m of water (H/Hb = ${fmt(sel.HoverHb)}): reduce the height to the breaking limit or use a breaking-wave load model.`);
    theory = sel.theory === 'airy' ? 'wheeler' : sel.theory === 'stokes5' ? 'stokes5' : 'stream';
    order = Math.max(sel.streamFunctionOrder, 5);
  }
  if (theory === 'stokes5') {
    let sol = null;
    try { sol = stokes5Wave({ H, d, T }); } catch { sol = null; }
    if (sol && sol.valid) {
      const eta = (x, t) => { const th = sol.k * (x - sol.c * t); let e = 0; for (let j = 1; j <= 5; j++) e += sol.etaAmp[j] * Math.cos(j * th); return e; };
      return { H, T, d, theory, label: 'Stokes 5th order', L: sol.L, k: sol.k, crest: sol.crest, surf: eta,
        kin: (x, z, t) => { const q = stokes5Kinematics(sol, { x, z, t }); return { u: q.u, w: q.w, ax: q.axTotal, az: q.azTotal }; } };
    }
    need(o.theory !== 'stokes5' && inp.theory !== 'stokes5' || o.fallback !== false, 'Stokes 5th-order theory is outside its range for this wave.');
    theory = 'stream';
  }
  if (theory === 'stream') {
    const sol = run(() => streamFunctionWave({ H, d, T, order }), 'Stream-function wave (the wave may exceed the breaking limit)');
    return { H, T, d, theory, label: `Stream function (order ${order})`, L: sol.L, k: sol.k, crest: sol.crest, surf: (x, t) => streamFunctionEta(sol, x, t),
      kin: (x, z, t) => { const q = streamFunctionKinematics(sol, { x, z, t }); return { u: q.u, w: q.w, ax: q.axTotal, az: q.azTotal }; } };
  }
  const a = H / 2, om = 2 * Math.PI / T;
  if (theory === 'wheeler') {
    return { H, T, d, theory, label: 'Airy with Wheeler stretching', L: L0, k: k0, crest: a, surf: (x, t) => a * Math.cos(k0 * x - om * t),
      kin: (x, z, t) => airyKinematics({ H, T, d, x, z, t, k: k0, stretching: 'wheeler', clipAboveSurface: false }) };
  }
  return { H, T, d, theory: 'airy', label: 'Airy, loaded to still water', L: L0, k: k0, crest: a, surf: () => 0,
    kin: (x, z, t) => airyKinematics({ H, T, d, x, z, t, k: k0, stretching: 'none', clipAboveSurface: false }) };
}

/**
 * Morison loading of the whole model for one heading. at(t) integrates drag + inertia on the wetted part of every
 * element as piecewise-linear distributed loads (returned as frame3d member loads on request) and their
 * resultants about the mudline centre; scan() steps the wave through one period.
 * Fh is the force along the heading, Mh the overturning moment about the horizontal axis normal to it.
 */
export function envLoader(J, inp, wf, betaDeg, o = {}) {
  const cb = Math.cos(betaDeg * DEG), sb = Math.sin(betaDeg * DEG), rho = inp.rhoW ?? 1025, sur = inp.surge || 0, d = wf.d;
  const Cd = o.Cd ?? inp.Cd, Cm = o.Cm ?? inp.Cm, kf = inp.kinFac ?? 1, tmg = (inp.mg || 0) / 1e3, zmg = -(inp.mgDepth ?? 0);
  const Uc0 = (o.Uc ?? inp.Uc ?? 0) * (inp.block ?? 1), nInt = Math.max(1, o.nInt ?? inp.nInt ?? (J.nSub >= 2 ? 2 : 3)), fac = o.factor ?? inp.appF ?? 1;
  const cur = (zw) => (Uc0 > 0 ? currentPowerLaw(Math.min(0, Math.max(zw, -d)), { v0: Uc0, d }) : 0);
  const cx = J.centre[0], cy = J.centre[1], zM = J.zMud;
  const f = [new Float64Array(3), new Float64Array(3)];
  const pt = (el, s, t, out) => {
    const x = el.a[0] + s * (el.b[0] - el.a[0]), y = el.a[1] + s * (el.b[1] - el.a[1]), z = el.a[2] + s * (el.b[2] - el.a[2]);
    const xp = (x - cx) * cb + (y - cy) * sb;
    let zw = z - sur;
    const top = wf.surf(xp, t);
    if (zw > top) zw = top;
    if (zw < -d) zw = -d;
    const q = wf.kin(xp, zw, t), uh = kf * q.u + cur(zw), ah = kf * q.ax;
    const vx = uh * cb, vy = uh * sb, vz = q.w, ax = ah * cb, ay = ah * sb, az = q.az;
    const vt = vx * el.t[0] + vy * el.t[1] + vz * el.t[2], at = ax * el.t[0] + ay * el.t[1] + az * el.t[2];
    const nx = vx - vt * el.t[0], ny = vy - vt * el.t[1], nz = vz - vt * el.t[2], vm = Math.sqrt(nx * nx + ny * ny + nz * nz);
    const Dm = el.D + (zw >= zmg ? 2 * tmg : 0), kD = 0.5 * rho * Cd * Dm * fac, kI = rho * Cm * Math.PI * Dm * Dm / 4 * fac;
    out[0] = kD * vm * nx + kI * (ax - at * el.t[0]); out[1] = kD * vm * ny + kI * (ay - at * el.t[1]); out[2] = kD * vm * nz + kI * (az - at * el.t[2]);
  };
  const at = (t, want = false) => {
    let Fx = 0, Fy = 0, Fz = 0, Mx = 0, My = 0;
    const items = want ? [] : null, elemF = want ? new Float64Array(J.elems.length) : null;
    if (wf.H > 0 || Uc0 > 0) for (const el of J.elems) {
      const xa = (el.a[0] - cx) * cb + (el.a[1] - cy) * sb, xb = (el.b[0] - cx) * cb + (el.b[1] - cy) * sb;
      const g0 = el.a[2] - sur - wf.surf(xa, t), g1 = el.b[2] - sur - wf.surf(xb, t);
      if (g0 > 0 && g1 > 0) continue;
      let s0 = 0, s1 = 1;
      if (g0 > 0) s0 = g0 / (g0 - g1); else if (g1 > 0) s1 = g0 / (g0 - g1);
      if (s1 - s0 < 1e-6) continue;
      const ds = (s1 - s0) / nInt;
      pt(el, s0, t, f[0]);
      for (let j = 0; j < nInt; j++) {
        const sa = s0 + j * ds, sb2 = sa + ds, fa = f[j % 2], fb = f[(j + 1) % 2];
        pt(el, sb2, t, fb);
        const Lp = ds * el.L, hx = 0.5 * Lp * (fa[0] + fb[0]), hy = 0.5 * Lp * (fa[1] + fb[1]), hz = 0.5 * Lp * (fa[2] + fb[2]);
        Fx += hx; Fy += hy; Fz += hz;
        // ∫ r × f over a piece with linear r and linear f: L [ra×fa + (ra×Δf + Δr×fa)/2 + Δr×Δf/3]
        const rax = el.a[0] + sa * (el.b[0] - el.a[0]) - cx, ray = el.a[1] + sa * (el.b[1] - el.a[1]) - cy, raz = el.a[2] + sa * (el.b[2] - el.a[2]) - zM;
        const drx = ds * (el.b[0] - el.a[0]), dry = ds * (el.b[1] - el.a[1]), drz = ds * (el.b[2] - el.a[2]);
        const dfx = fb[0] - fa[0], dfy = fb[1] - fa[1], dfz = fb[2] - fa[2];
        Mx += Lp * ((ray * fa[2] - raz * fa[1]) + 0.5 * ((ray * dfz - raz * dfy) + (dry * fa[2] - drz * fa[1])) + (dry * dfz - drz * dfy) / 3);
        My += Lp * ((raz * fa[0] - rax * fa[2]) + 0.5 * ((raz * dfx - rax * dfz) + (drz * fa[0] - drx * fa[2])) + (drz * dfx - drx * dfz) / 3);
        if (want) {
          items.push({ member: el.e, type: 'linear', q1: [fa[0], fa[1], fa[2]], q2: [fb[0], fb[1], fb[2]], axes: 'global', s1: sa * el.L, s2: Math.min(sb2 * el.L, el.L) });
          elemF[el.e] += hx * cb + hy * sb;
        }
      }
    }
    return { Fx, Fy, Fz, Fh: Fx * cb + Fy * sb, Mh: -Mx * sb + My * cb, items, elemF };
  };
  const scan = (nPhase = inp.nPhase ?? 24) => {
    const T = wf.T, t = [], Fh = [], Mh = [], Fz = [];
    let iF = 0, iM = 0;
    const n = wf.H > 0 ? nPhase : 1;
    for (let i = 0; i < n; i++) {
      const r = at(i * T / n);
      t.push(i * T / n); Fh.push(r.Fh); Mh.push(r.Mh); Fz.push(r.Fz);
      if (Math.abs(r.Fh) > Math.abs(Fh[iF])) iF = i;
      if (Math.abs(r.Mh) > Math.abs(Mh[iM])) iM = i;
    }
    return { t, Fh, Mh, Fz, iF, iM, FhMax: Fh[iF], MhMax: Mh[iM], range: Math.max(...Fh) - Math.min(...Fh) };
  };
  return { at, scan, cb, sb, betaDeg };
}

/** Wind force on the topsides area (API RP 2A Eq. 5.6 with a power-law height coefficient), N. */
export function windLoad(inp, Uw = inp.Uw) {
  if (!(Uw > 0) || !(inp.windA > 0)) return 0;
  return windForce({ U: Uw, A: inp.windA, Cs: inp.windCs ?? 1, Ch: heightCoefficient(Math.max(inp.windZ ?? 10, 1), 10, 0.125) }).F;
}

/**
 * Nodal loads on the leg tops: topsides weight, horizontal interface + wind force along the heading and the
 * overturning moment as a push–pull. Returns the loads and their resultants about the mudline centre.
 */
export function topLoads(J, inp, betaDeg, o = {}) {
  const fG = o.fG ?? 1, fE = o.fE ?? 1, cb = Math.cos(betaDeg * DEG), sb = Math.sin(betaDeg * DEG), P = J.model.nodes, n = J.topNodes.length;
  const W = (o.W ?? inp.Wtop * 1e6) * fG, Fw = windLoad(inp, o.Uw ?? inp.Uw) * fE, Fi = (inp.FxTop || 0) * 1e3 * fE, Hx = Fw + Fi;
  const M = ((inp.MyTop || 0) * 1e6 + Fw * Math.max((inp.windZ ?? 0) - J.zTop, 0)) * fE;
  const xi = J.topNodes.map((k) => (P[k][0] - J.centre[0]) * cb + (P[k][1] - J.centre[1]) * sb), xm = xi.reduce((s, v) => s + v, 0) / n;
  const s2 = xi.reduce((s, v) => s + (v - xm) ** 2, 0);
  const nodal = J.topNodes.map((k, j) => {
    const F = [Hx * cb / n, Hx * sb / n, -W / n, 0, 0, 0];
    // push–pull couple (down-wind legs pushed down); a single top node takes the moment directly
    if (s2 > 1e-9) F[2] -= M * (xi[j] - xm) / s2; else { F[3] = -M * sb / n; F[4] = M * cb / n; }
    return { node: k, F };
  });
  return { nodal, W, Fwind: Fw, Hx, M, Mh: M + Hx * (J.zTop - J.zMud) };
}

/** Self-weight with buoyancy (marine method) and marine growth, as a frame3d load-case fragment. */
export function selfWeight(J, inp, factor = 1) {
  const f = factor * (inp.swCont ?? 1);
  if (!(f > 0)) return null;
  const sw = { factor: f, buoyancy: { waterline: inp.surge || 0, rhoW: inp.rhoW ?? 1025 } };
  if (inp.mg > 0) sw.marineGrowth = { thickness: inp.mg / 1e3, rho: 1325, zTop: inp.surge || 0, zBottom: -(inp.mgDepth ?? 0) };
  return sw;
}

/** One complete load case: factored gravity + topsides + environmental member loads. */
export function loadCase(name, J, inp, { betaDeg = 0, fG = 1, fE = 1, env = null, top = true, sw = true, W, Uw } = {}) {
  const lc = { name, nodal: [], member: [] };
  if (top) lc.nodal = topLoads(J, inp, betaDeg, { fG, fE, W, Uw }).nodal;
  if (sw) { const s = selfWeight(J, inp, fG); if (s) lc.selfWeight = s; }
  if (env && env.items) lc.member = fE === 1 ? env.items : env.items.map((it) => ({ ...it, q1: it.q1.map((v) => v * fE), q2: it.q2.map((v) => v * fE) }));
  return lc;
}

export function solve(J, cases, hint = 'Frame analysis') {
  J.model.loadCases = cases;
  const sol = run(() => analyzeFrame(J.model), `${hint} (check that the structure is connected and supported)`);
  for (const c of cases) {
    const r = sol.results[c.name];
    need(r && r.displacements && Number.isFinite(r.displacements[0]), 'The frame analysis did not return a finite solution: the structure is a mechanism or unsupported.');
  }
  return sol.results;
}

/* ----------------------------------------------------------------------------------------- post-processing */

/** Extreme section forces of every original member (signed axial of largest magnitude; other terms absolute). */
export function memberForces(J, res) {
  const nSt = J.nSub >= 2 ? 3 : 5;   // stations per element: five or more along every member
  return J.orig.map((m, i) => {
    let N = 0, My = 0, Mz = 0, V = 0, T = 0;
    for (const k of J.memberMap[i]) {
      const d = memberDiagram(res, k, nSt), n1 = d.N[0], n2 = d.N[d.N.length - 1];
      if (Math.abs(n1) > Math.abs(N)) N = n1;
      if (Math.abs(n2) > Math.abs(N)) N = n2;
      My = Math.max(My, d.maxAbs.My); Mz = Math.max(Mz, d.maxAbs.Mz); V = Math.max(V, Math.hypot(d.maxAbs.Vy, d.maxAbs.Vz)); T = Math.max(T, d.maxAbs.T);
    }
    return { N, My, Mz, V, T };
  });
}

/** Code check of one member for given forces. Non-tubular imported sections fall back on von Mises against fy/γ. */
export function checkMember(J, m, f, inp, o = {}) {
  const fy = J.steel.fy;
  const ss = sectionStress({ section: m.sec, N: f.N, Vy: f.V, T: f.T, My: f.My, Mz: f.Mz });
  let uc = ss.vonMises / (fy / 1.15), ref = 'von Mises, γM = 1.15', strengths = null;
  if (m.tube) {
    const code = inp.code ?? 'ISO', opt = { Ky: m.K, Kz: m.K, E: J.steel.E };
    if (code === 'ISO') opt.gamma = { t: inp.gRt ?? 1.05, c: inp.gRc ?? 1.18, b: inp.gRb ?? 1.05 };
    if (o.p > 0 && code === 'ISO') opt.p = o.p;
    const c = tubularMemberCheck({ code, section: m.sec, fy, L: m.L, forces: { N: f.N, My: f.My, Mz: f.Mz, Vy: f.V, Vz: 0, T: f.T }, options: opt });
    uc = c.governing.uc; ref = c.governing.name; strengths = c.strengths;
  }
  return { uc: cap(uc), ref, vm: ss.vonMises, sigA: ss.sigmaAxial, sigB: ss.sigmaBending, tau: ss.tau, strengths };
}

/** Member rows (forces + check) for one solved case. */
export function memberRows(J, res, inp, o = {}) {
  const F = memberForces(J, res);
  return J.orig.map((m, i) => {
    const f = F[i], c = checkMember(J, m, f, inp, o.hydro && !m.flooded && m.zmid < 0 ? { p: (inp.rhoW ?? 1025) * G * (-(m.zmid) + (inp.surge || 0)) * (inp.gG ?? 1) } : {});
    return { i, group: m.group, label: m.label, L: m.L, D: m.D, N: f.N, M: Math.hypot(f.My, f.Mz), My: f.My, Mz: f.Mz, V: f.V, T: f.T, vm: c.vm, sigA: c.sigA, sigB: c.sigB, tau: c.tau, uc: c.uc, ref: c.ref, case: o.name ?? '' };
  });
}

/** Envelope of member rows over several cases: the row with the largest utilisation per member. */
export function envelope(list) {
  return list[0].map((r, i) => list.reduce((a, rows) => (rows[i].uc > a.uc ? rows[i] : a), r));
}

export function byGroup(rows) {
  const g = {};
  for (const r of rows) { const q = g[r.group] || (g[r.group] = { group: r.group, label: r.label, n: 0, uc: -1 }); q.n++; if (r.uc > q.uc) Object.assign(q, { uc: r.uc, N: r.N, M: r.M, vm: r.vm, ref: r.ref, i: r.i, case: r.case }); }
  return Object.values(g);
}

/** Support reactions as pile-head forces: axial (compression positive), shear and moment per pile. */
export function pileForces(J, res) {
  const piles = res.reactions.map((x) => ({ node: x.node, axial: x.R[2], shear: Math.hypot(x.R[0], x.R[1]), moment: Math.hypot(x.R[3], x.R[4]), Rx: x.R[0], Ry: x.R[1] }));
  return { piles, axialMax: Math.max(...piles.map((p) => p.axial)), tensionMax: Math.max(0, ...piles.map((p) => -p.axial)), shearMax: Math.max(...piles.map((p) => p.shear)),
    momentMax: Math.max(...piles.map((p) => p.moment)), sumV: piles.reduce((s, p) => s + p.axial, 0), sumH: Math.hypot(piles.reduce((s, p) => s + p.Rx, 0), piles.reduce((s, p) => s + p.Ry, 0)) };
}

export function displacements(J, res) {
  const u = res.displacements, n = J.model.nodes.length, mag = new Array(n);
  let max = 0, rot = 0;
  for (let i = 0; i < n; i++) { mag[i] = Math.hypot(u[6 * i], u[6 * i + 1], u[6 * i + 2]); if (mag[i] > max) max = mag[i]; rot = Math.max(rot, Math.hypot(u[6 * i + 3], u[6 * i + 4], u[6 * i + 5])); }
  const top = J.topNodes.reduce((s, k) => Math.max(s, Math.hypot(u[6 * k], u[6 * k + 1])), 0);
  const mud = J.mudNodes.reduce((s, k) => Math.max(s, Math.hypot(u[6 * k], u[6 * k + 1])), 0);
  need(Number.isFinite(max), 'The frame is unstable (mechanism): check supports and connectivity.');
  return { u, mag, max, rotMax: rot, top, mud };
}

/** Utilisation of the members meeting at each node of the refined model (for model3d colouring). */
export function nodeValuesFromMembers(J, values) {
  const out = new Array(J.model.nodes.length).fill(0);
  J.memberMap.forEach((list, i) => list.forEach((e) => { const m = J.model.members[e]; out[m.a] = Math.max(out[m.a], values[i]); out[m.b] = Math.max(out[m.b], values[i]); }));
  return out;
}

/** model3d plot of the refined model with nodal values and an exaggerated deformed shape (u: 6 dof per node). */
export function plot3d(J, { title, values, label, unit, u, explain, shape }) {
  const P = J.model.nodes, p = { type: 'model3d', title, nodes: P.map((n) => [round(n[0]), round(n[1]), round(n[2])]), lines: J.model.members.map((m) => [m.a, m.b]), nodeValues: values.map((v) => round(v)), valueLabel: label, valueUnit: unit };
  const d = u || shape;
  let sc = 0;
  if (d) {
    let mx = 0;
    for (let i = 0; i < P.length; i++) mx = Math.max(mx, Math.hypot(d[6 * i], d[6 * i + 1], d[6 * i + 2]));
    sc = mx > 0 ? 0.06 * J.size / mx : 0;
    p.deformed = P.map((n, i) => [round(n[0] + sc * d[6 * i]), round(n[1] + sc * d[6 * i + 1]), round(n[2] + sc * d[6 * i + 2])]);
  }
  p.explain = `${explain}${u ? ` The deformed shape is exaggerated ${fmt(sc)} times.` : ''}`;
  return p;
}

export const memberCols = [{ key: 'i', label: 'Member' }, { key: 'group', label: 'Group' }, { key: 'L', label: 'Length', unit: 'm' }, { key: 'N', label: 'Axial (tension +)', unit: 'kN' }, { key: 'M', label: 'Moment', unit: 'kN·m' }, { key: 'V', label: 'Shear', unit: 'kN' }, { key: 'T', label: 'Torsion', unit: 'kN·m' }, { key: 'vm', label: 'von Mises', unit: 'MPa' }, { key: 'uc', label: 'Utilisation' }, { key: 'ref', label: 'Governing check' }, { key: 'case', label: 'Load case' }];
export const memberTableRows = (rows, n = 15) => rows.slice().sort((a, b) => b.uc - a.uc).slice(0, n).map((r) => ({ i: r.i, group: r.label, L: round(r.L), N: round(r.N / 1e3), M: round(r.M / 1e3), V: round(r.V / 1e3), T: round(r.T / 1e3), vm: round(r.vm / 1e6), uc: round(r.uc), ref: r.ref, case: r.case }));

/**
 * Brace-to-chord connections of the model: at every node the member of largest diameter is the chord and each
 * smaller tubular member ending there is a brace. Equal-diameter crossings (X-brace nodes) are not listed.
 */
export function jointList(J) {
  const at = new Map();
  J.orig.forEach((m) => { for (const [n, end] of [[m.a, 'a'], [m.b, 'b']]) { if (!at.has(n)) at.set(n, []); at.get(n).push({ m, end }); } });
  const P = J.base.nodes, out = [];
  const dir = (m) => { const a = P[m.a], b = P[m.b]; return [(b[0] - a[0]) / m.L, (b[1] - a[1]) / m.L, (b[2] - a[2]) / m.L]; };
  for (const [node, list] of at) {
    const tubes = list.filter((x) => x.m.tube);
    if (tubes.length < 2) continue;
    const chord = tubes.reduce((a, x) => (x.m.sec.D > a.m.sec.D || (x.m.sec.D === a.m.sec.D && x.m.sec.t > a.m.sec.t) ? x : a));
    const tc = dir(chord.m);
    for (const x of tubes) {
      if (x === chord || x.m.sec.D >= chord.m.sec.D * 0.999) continue;
      const tb = dir(x.m), c = Math.abs(tb[0] * tc[0] + tb[1] * tc[1] + tb[2] * tc[2]);
      out.push({ node, brace: x.m.i, end: x.end, chord: chord.m.i, chordEnd: chord.end, thetaDeg: Math.acos(Math.min(1, c)) / DEG, z: P[node][2] });
    }
  }
  return out;
}

/** End forces {N, Vy, Vz, T, My, Mz} of an original member at end 'a' or 'b' for a solved case. */
export function endForces(J, res, i, end) {
  const list = J.memberMap[i], mf = res.memberForces[end === 'a' ? list[0] : list[list.length - 1]];
  return end === 'a' ? mf.endA : mf.endB;
}

/** Governing regular-wave environment for one heading: scan, then the load state at the phase of maximum base shear. */
export function designWave(J, inp, betaDeg, o = {}) {
  const wf = o.wf || waveField(inp, o), ld = envLoader(J, inp, wf, betaDeg, o), sc = ld.scan(o.nPhase ?? inp.nPhase);
  const top = topLoads(J, inp, betaDeg, {});
  return { wf, ld, sc, top, baseShear: Math.abs(sc.FhMax) + top.Hx, otm: Math.abs(sc.MhMax) + top.Mh };
}

/**
 * Load cases of one design condition (wave + current + wind + gravity) for a list of headings: the wave is scanned
 * through its phase per heading and a case is generated at the phase of maximum base shear and, when it differs by
 * more than one step, at the phase of maximum overturning moment.
 * @returns {{wf:object, dirs:Array, cases:Array}} dirs[i] = { beta, sc, top, baseShear, otm (unfactored, N, N·m), names }
 */
export function conditionCases(J, inp, c) {
  const wf = waveField(inp, { H: c.H, T: c.T, theory: c.theory }), cases = [], dirs = [];
  for (const beta of c.headings) {
    const ld = envLoader(J, inp, wf, beta, { Uc: c.Uc }), sc = ld.scan(c.nPhase ?? inp.nPhase);
    const top = topLoads(J, inp, beta, { Uw: c.Uw, W: c.W });
    const sgn = Math.sign(sc.FhMax) || 1;
    const d = { beta, sc, top, baseShear: Math.abs(sc.FhMax + sgn * top.Hx), otm: Math.abs(sc.MhMax + (Math.sign(sc.MhMax) || 1) * top.Mh), names: [], ld };
    const n = sc.t.length, gap = Math.min(Math.abs(sc.iF - sc.iM), n - Math.abs(sc.iF - sc.iM));
    for (const [k, i] of [['F', sc.iF], ['M', sc.iM]]) {
      if (k === 'M' && gap <= 1) continue;
      const name = `${c.tag} ${fmt(beta)}° ${k === 'F' ? 'max shear' : 'max OTM'}`;
      cases.push(loadCase(name, J, inp, { betaDeg: beta, fG: c.fG, fE: c.fE, env: ld.at(sc.t[i], true), W: c.W, Uw: c.Uw }));
      d.names.push(name);
    }
    dirs.push(d);
  }
  return { wf, dirs, cases };
}

export const headingList = (v) => (String(v) === '8' ? [0, 45, 90, 135, 180, 225, 270, 315] : String(v) === '1' ? [0] : [0, 45, 90]);
export const headingInput = (def = '3') => ({ key: 'nDir', label: 'Wave headings analysed', type: 'select', default: def, group: 'Environment', role: 'bc',
  options: [{ value: '1', label: 'Broadside only (0°)' }, { value: '3', label: '0°, 45°, 90° (symmetric jacket)' }, { value: '8', label: 'Every 45° (eight headings)' }],
  help: 'A doubly symmetric jacket needs only 0°, 45° and 90°; use eight headings for an imported or asymmetric model.' });

/** Factored storm load case of one heading at the wave phase of maximum base shear ('F') or overturning moment ('M'). */
export function stormCase(J, inp, { name = 'Storm', beta = inp.heading ?? 0, fG = inp.gG ?? 1, fE = inp.gE ?? 1, at = 'F', top = true, sw = true, W } = {}) {
  const dw = designWave(J, inp, beta), env = dw.ld.at(dw.sc.t[at === 'M' ? dw.sc.iM : dw.sc.iF], true);
  return { ...dw, env, name, beta, lc: loadCase(name, J, inp, { betaDeg: beta, fG, fE, env, top, sw, W }) };
}

/**
 * Linear (eigenvalue) buckling of the model for a reference load case: (K + λ Kg) φ = 0. Subspace iteration on
 * (−Kg) φ = θ (K + σ Kg) φ with a shift σ placed just below the first factor after a few unshifted sweeps, so the
 * closely spaced brace modes of a jacket separate quickly and the reversed-load (negative) factors drop out.
 * @returns {{factors:number[], modes:Float64Array[], axial:number[], iterations:number}} factors ascending; modes over all dofs, largest component 1
 */
export function linearBuckling(J, caseName, nModes = 3) {
  const model = J.model, sys = run(() => buildFrameSystem(model, {}), 'Buckling model');
  let fac = run(() => factorFrame(sys), 'Buckling analysis (check supports and connectivity)');
  const ref = solveFrame(sys, fac, resolveLoads(sys, loadParts(model, caseName)));
  const axial = ref.memberForces.map((m) => 0.5 * (m.endA.N + m.endB.N));
  need(Math.min(...axial) < 0, 'No member is in compression in the reference load case: nothing can buckle. Add gravity or environmental load.');
  const kg = frameGeometricStiffness(sys, axial).global, n = sys.nFree, p = Math.min(nModes, n), q = Math.min(n, Math.max(2 * p, p + 6));
  const neg = (x) => { const y = elementMatVec(sys, kg, x); for (let i = 0; i < n; i++) y[i] = -y[i]; return y; };
  // deterministic start vectors (linear congruential sequence)
  let seed = 12345; const rnd = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648 - 0.5; };
  let X = Array.from({ length: q }, () => Float64Array.from({ length: n }, rnd));
  let sigma = 0, lam = [], prev = [], it = 0, shifted = false;
  const dot = (a, b) => { let s = 0; for (let k = 0; k < n; k++) s += a[k] * b[k]; return s; };
  while (it < 60) {
    it++;
    // Z = (K + σKg)⁻¹ (−Kg) X, so (K + σKg) Z = Y and the reduced matrices need no product with K
    const Y = X.map(neg), Z = Y.map((y) => fac.solve(y));
    Z.forEach((z, j) => { const sc = 1 / Math.sqrt(dot(z, z) || 1), y = Y[j]; for (let k = 0; k < n; k++) { z[k] *= sc; y[k] *= sc; } });
    const Wz = Z.map(neg), A = Array.from({ length: q }, () => new Array(q)), B = Array.from({ length: q }, () => new Array(q));
    for (let i = 0; i < q; i++) for (let j = i; j < q; j++) { A[i][j] = A[j][i] = 0.5 * (dot(Z[i], Wz[j]) + dot(Z[j], Wz[i])); B[i][j] = B[j][i] = 0.5 * (dot(Z[i], Y[j]) + dot(Z[j], Y[i])); }
    let eg;
    try { eg = eigGenSym(A, B); } catch { const tr = B.reduce((s2, r, i) => s2 + r[i], 0) / q; for (let i = 0; i < q; i++) B[i][i] += 1e-10 * tr; eg = eigGenSym(A, B); }
    const order = eg.values.map((v, i) => i).sort((i, j) => eg.values[j] - eg.values[i]);
    X = order.map((c) => { const x = new Float64Array(n), v = eg.vectors[c]; for (let j = 0; j < q; j++) { const zj = v[j], z = Z[j]; for (let k = 0; k < n; k++) x[k] += zj * z[k]; } return x; });
    lam = order.map((c) => eg.values[c]).filter((v) => v > 1e-12).map((th) => sigma + 1 / th);
    need(lam.length, 'No positive buckling factor was found: the reference load case does not compress any member critically.');
    if ((shifted || it > 25) && lam.length >= Math.min(p, prev.length || p) && lam.slice(0, p).every((v, c) => Math.abs(v - prev[c]) < 1e-6 * v)) break;
    prev = lam.slice(0, p);
    if (!shifted && it >= 4) {
      // place the shift below the current (upper-bound) estimate; back off while K + σKg is not positive definite
      let s = 0.9 * lam[0];
      for (let k = 0; k < 8; k++) {
        try { const f2 = factorFrame(sys, { add: kg, addFactor: s }); if (f2.negativePivots === 0) { fac = f2; sigma = s; shifted = true; prev = []; break; } } catch { /* singular at the shift: move it */ }
        s *= 0.8;
      }
      if (!shifted) shifted = true;
    }
  }
  const factors = lam.slice(0, p), modes = factors.map((f, c) => {
    const full = new Float64Array(sys.nDof); let big = 0;
    for (let k = 0; k < n; k++) if (Math.abs(X[c][k]) > Math.abs(big)) big = X[c][k];
    for (let k = 0; k < n; k++) full[sys.eqDof[k]] = X[c][k] / (big || 1);
    return full;
  });
  return { factors, modes, axial, iterations: it, reference: ref };
}

/** Small imported-frame models used by the verification cases (closed-form benchmarks): kind 'frame', metres. */
export function testFrame(kind, o = {}) {
  const sec = { type: 'tube', D: o.D ?? 1.5, t: o.t ?? 0.045 }, H = o.H ?? 95, s = o.span ?? 20;
  if (kind === 'portal') return { kind: 'frame', name: 'portal', nodes: new Float64Array([0, 0, 0, s, 0, 0, 0, 0, H, s, 0, H]), members: [{ a: 0, b: 2, section: sec, group: 'leg' }, { a: 1, b: 3, section: sec, group: 'leg' }, { a: 2, b: 3, section: sec, group: 'horizontal' }] };
  if (kind === 'tjoint') {
    // 20 m chord 1000 × 25 with a 500 × 12.5 brace at mid-height: horizontal stub (θ = 90°) or rising at 45° to the loaded top node
    const ch = { type: 'tube', D: 1.0, t: 0.025 }, br = { type: 'tube', D: 0.5, t: 0.0125 }, tip = o.thetaDeg === 45 ? [20, 0, 30] : [5, 0, 10];
    return { kind: 'frame', name: 'tjoint', nodes: new Float64Array([0, 0, 0, 0, 0, 10, 0, 0, 20, ...tip]), members: [{ a: 0, b: 1, section: ch, group: 'leg' }, { a: 1, b: 2, section: ch, group: 'leg' }, { a: 1, b: 3, section: br, group: 'brace' }] };
  }
  return { kind: 'frame', name: 'column', nodes: new Float64Array([0, 0, 0, 0, 0, H]), members: [{ a: 0, b: 1, section: sec, group: 'leg' }] };
}

/** Full input set of a storm analysis on the shared model (geometry, gravity, environment, foundation, criteria, solver). */
export const stdInputs = (o = {}) => [...geomInputs(o.geom), ...topsidesInputs(o.top), ...envInputs(o.env), ...foundationInputs(o.found), ...criteriaInputs(o.crit), ...meshInputs(o.mesh)];
export const steelMass = (J, inp) => J.steelT * (inp.swCont > 0 ? inp.swCont : 1);
export const groupSteel = (J, group) => J.orig.filter((m) => m.group === group).reduce((s, m) => s + m.sec.A * m.L * J.steel.rho, 0) / 1e3;
export const foundNote = (J) => `Foundation: ${J.K ? `pile-head springs (${J.K.source}) Kh ${fmt(J.K.Kh / 1e6)} MN/m, Kr ${fmt(J.K.Kr / 1e9)} GN·m/rad, Kv ${fmt(J.K.Kv / 1e6)} MN/m` : `${J.found} supports at the mudline`}.`;
export const maxLegD = (J) => Math.max(...J.orig.map((m) => m.D));

/** One member-strength check per member group from byGroup() output. */
export const groupChecks = (groups, inp, tag, explain) => groups.map((g) => ({ name: `Member strength${tag ? ` (${tag})` : ''} — ${g.label}`, demand: round(g.uc), capacity: 1, unit: 'interaction', uc: g.uc,
  ref: `${inp.code ?? 'ISO'} tubular member check: ${g.ref}`, explain: explain ?? 'Governing tension, compression, bending, shear and beam-column interaction of the most loaded member of this group. Above 1.0 the member yields or buckles before the factored load is reached.' }));

/** Loads of a table [{ z, Fx, Fy, Fz }] (kN) placed on the leg nodes nearest to elevation z and shared equally between them. */
export function tableLoads(J, rows, factor = 1) {
  const P = J.base.nodes, cand = [...new Set(J.orig.filter((m) => m.group === 'leg' || J.imported).flatMap((m) => [m.a, m.b]))].filter((n) => !J.mudNodes.includes(n));
  const out = [];
  for (const r of rows || []) {
    const F = [Number(r.Fx) || 0, Number(r.Fy) || 0, Number(r.Fz) || 0].map((v) => v * 1e3 * factor);
    if (!F.some((v) => v !== 0)) continue;
    need(Number.isFinite(Number(r.z)), 'Every load row needs an elevation.');
    need(cand.length, 'The model has no free node to carry the tabulated loads.');
    const dz = Math.min(...cand.map((n) => Math.abs(P[n][2] - r.z))), at = cand.filter((n) => Math.abs(P[n][2] - r.z) <= dz + 0.25);
    for (const n of at) out.push({ node: n, F: [F[0] / at.length, F[1] / at.length, F[2] / at.length, 0, 0, 0] });
  }
  return out;
}

/** Gravity (held) and environmental (scaled) load cases for an event-to-event pushover of one heading. */
export function pushCases(J, inp, st) {
  const cases = [loadCase('Gravity', J, inp, { betaDeg: st.beta, fG: inp.gG ?? 1, fE: 0 }), loadCase('Env', J, inp, { betaDeg: st.beta, fG: 0, fE: inp.gE ?? 1, env: st.env, sw: false, W: 0 })];
  J.model.loadCases = cases;
  return cases;
}
/** Decimated pushover curve as plot series values (displacement mm, base shear MN, load factor). */
export function pushCurve(po) {
  const n = po.curve.loadFactor.length, step = Math.max(1, Math.ceil(n / 400)), d = [], v = [], lf = [];
  for (let i = 0; i < n; i += step) { d.push(round(Math.abs(po.curve.displacement[i]) * 1e3)); v.push(round(Math.abs(po.curve.baseShear[i]) / 1e6)); lf.push(round(po.curve.loadFactor[i])); }
  return { d, v, lf };
}

const dot3 = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross3 = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const unit3 = (a) => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };

/**
 * Geometry of every brace-to-chord connection for joint, hot-spot and fatigue checks: brace and chord dimensions
 * (the chord wall is the joint-can thickness on legs), brace angle, the joint plane and the group of coplanar
 * braces at the node used to classify K / X / Y action.
 */
export function jointData(J, inp) {
  const P = J.base.nodes, canT = (inp.canT || 0) / 1e3;
  const list = jointList(J).map((q, k) => {
    const b = J.orig[q.brace], c = J.orig[q.chord], ax = memberAxes(P[b.a], P[b.b]).R, tc = memberAxes(P[c.a], P[c.b]).R[0];
    const away = q.end === 'a' ? ax[0] : ax[0].map((v) => -v);   // brace axis pointing away from the joint
    let n = cross3(away, tc); const ln = Math.hypot(n[0], n[1], n[2]);
    n = ln > 1e-6 ? n.map((v) => v / ln) : ax[1];
    const e = unit3(away.map((v, i) => v - dot3(away, tc) * tc[i]));   // in-plane normal to the chord, towards the brace
    return { k, ...q, b, c, d: b.sec.D, t: b.sec.t, D: c.sec.D, T: c.group === 'leg' && canT > c.sec.t ? canT : c.sec.t, Lc: c.L, thetaDeg: Math.max(q.thetaDeg, 1), n, e, ey: ax[1], ez: ax[2], eo: cross3(n, away), tc, plane: -1, side: 1 };
  });
  // coplanar braces at the same node form one planar joint
  let np = 0;
  for (const j of list) {
    if (j.plane >= 0) continue;
    j.plane = np++; j.side = 1;
    for (const q of list) if (q.plane < 0 && q.node === j.node && Math.abs(dot3(q.n, j.n)) > 0.95) { q.plane = j.plane; q.side = dot3(q.e, j.e) >= 0 ? 1 : -1; }
  }
  return list;
}

/** Brace and chord loads of every joint for a solved case: P (tension +), in-plane and out-of-plane moments, nominal stresses. */
export function jointLoads(J, res, joints) {
  return joints.map((j) => {
    const f = endForces(J, res, j.brace, j.end), fc = endForces(J, res, j.chord, j.chordEnd);
    const Mipb = f.My * dot3(j.ey, j.n) + f.Mz * dot3(j.ez, j.n), Mopb = f.My * dot3(j.ey, j.eo) + f.Mz * dot3(j.ez, j.eo), W = j.b.sec.Wely ?? j.b.sec.Iy / (j.d / 2);
    return { P: f.N, Mipb, Mopb, sx: f.N / j.b.sec.A, sip: Mipb / W, sop: Mopb / W, chord: { P: fc.N, Mipb: Math.abs(fc.My), Mopb: Math.abs(fc.Mz) } };
  });
}
