/**
 * Pipeline crossings and supports: a beam on discrete supports and a unilateral seabed (reactions, bending
 * stress, spans created, separation), support quantities, impact energies and protection capacity, trawl
 * pull-over loading and rock armour sizing.
 *
 * Design notes for crossings: the crossing and crossed lines must stay electrically separated (cathodic
 * protection systems isolated by the physical gap — mattress, sleeper or rock — typically not less than
 * 0.3 m), and support settlement should be included as a reduced support height.
 * All quantities SI. In the beam model loads act downward-positive and deflections are upward-positive.
 */

import { solveBandedSym } from '../math/linalg.js';
import { coverQuantities } from './pipesoil.js';

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function nonNeg(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}
function fin(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v))) throw new RangeError(`${name} must be a finite number`);
  return v;
}

/**
 * Euler–Bernoulli beam on discrete supports and, optionally, a unilateral elastic seabed, under uniform
 * weight and point loads. Hermite cubic finite elements on a uniform mesh; supports and point loads act at
 * the nearest node. Supports and the seabed can be unilateral (the pipe may lift off; found by an
 * active-set iteration). Rigid supports are represented by a spring 10⁷ times the element stiffness.
 *
 * Typical use: a pipe laid over a crossing — give the sleeper / mattress / berm crest heights as supports
 * above a seabed at elevation zero; the solution returns the touchdown points, the spans created, the
 * support reactions and the bending moments. Settlement of a support is modelled by lowering its height.
 *
 * Reference: beam finite-element method (exact at the nodes for uniform load); contact by active set.
 * Validity: small deflections, no axial force (tension stiffening and effective compression are not
 * included), linear elastic pipe.
 *
 * @param {object} o
 * @param {number} o.L modelled length [m]
 * @param {number} o.EI bending stiffness [N·m²]
 * @param {number} [o.w=0] uniform downward load (submerged weight) [N/m]
 * @param {Array<{x:number, height?:number, stiffness?:number, unilateral?:boolean}>} [o.supports=[]]
 *   position [m], top elevation [m] (default 0), vertical stiffness [N/m] (default rigid), lift-off allowed (default true)
 * @param {{stiffness:number, elevation?:number}} [o.seabed] Winkler seabed: stiffness [N/m²] per unit length,
 *   elevation [m] (default 0); contact only in compression
 * @param {Array<{x:number, F:number}>} [o.pointLoads=[]] downward point loads [N]
 * @param {'free'|'pinned'|'fixed'} [o.ends='free'] condition at both ends (pinned / fixed at elevation 0)
 * @param {number} [o.nElem=200] number of elements (4–5000) [-]
 * @param {number} [o.D] steel outer diameter [m] (with E, gives the bending stress)
 * @param {number} [o.E=207e9] [Pa]
 * @returns {{x:Float64Array, deflection:Float64Array, rotation:Float64Array, moment:Float64Array, shear:Float64Array,
 *   reactions:Array<{x:number, R:number, active:boolean}>, seabedReaction:Float64Array, spans:Array<{from:number, to:number, length:number}>,
 *   maxMoment:number, xMaxMoment:number, maxStress:(number|null), totalReaction:number, totalLoad:number, iterations:number}}
 *   deflection [m] (elevation of the pipe axis relative to its seabed-resting position); moment [N·m]
 *   (positive sagging); shear = dM/dx [N] (averaged across a node where a reaction acts); R upward support reactions [N]; seabedReaction [N/m]; spans between
 *   contact points [m]; maxStress [Pa]
 */
export function beamOnSupports({ L, EI, w = 0, supports = [], seabed, pointLoads = [], ends = 'free', nElem = 200, D, E = 207e9 }) {
  pos('L', L); pos('EI', EI); fin('w', w);
  const n = Math.floor(nElem);
  if (!(n >= 4 && n <= 5000)) throw new RangeError('nElem must be within 4–5000');
  const h = L / n;
  const nn = n + 1;
  const nd = 2 * nn;
  const kRef = (12 * EI) / h ** 3;
  const kRigid = 1e7 * kRef;
  const nodeOf = (x, name) => {
    if (!(x >= 0 && x <= L)) throw new RangeError(`${name} must lie within [0, L]`);
    return Math.round(x / h);
  };
  const sup = supports.map((s, i) => ({
    node: nodeOf(s.x, `supports[${i}].x`), height: s.height === undefined ? 0 : fin(`supports[${i}].height`, s.height),
    k: s.stiffness === undefined || s.stiffness === Infinity ? kRigid : pos(`supports[${i}].stiffness`, s.stiffness),
    unilateral: s.unilateral !== false, active: true, x: s.x,
  }));
  const bedK = seabed ? pos('seabed.stiffness', seabed.stiffness) : 0;
  const bedZ = seabed && seabed.elevation !== undefined ? seabed.elevation : 0;
  const bedActive = new Uint8Array(nn).fill(seabed ? 1 : 0);
  const trib = (i) => (i === 0 || i === n ? h / 2 : h);

  // beam stiffness in symmetric band storage (band[i][k] = K[i][i+k], k = 0..3) and load vector
  const band = Array.from({ length: nd }, () => new Float64Array(4));
  const f0 = new Float64Array(nd);
  const c = EI / h ** 3;
  const ke = [[12, 6 * h, -12, 6 * h], [6 * h, 4 * h * h, -6 * h, 2 * h * h], [-12, -6 * h, 12, -6 * h], [6 * h, 2 * h * h, -6 * h, 4 * h * h]];
  const q = -w;                                   // upward load intensity
  const fe = [(q * h) / 2, (q * h * h) / 12, (q * h) / 2, (-q * h * h) / 12];
  for (let e = 0; e < n; e++) {
    for (let a = 0; a < 4; a++) {
      f0[2 * e + a] += fe[a];
      for (let b = a; b < 4; b++) band[2 * e + a][b - a] += c * ke[a][b];
    }
  }
  pointLoads.forEach((p, i) => { f0[2 * nodeOf(p.x, `pointLoads[${i}].x`)] -= fin(`pointLoads[${i}].F`, p.F); });
  if (ends !== 'free' && ends !== 'pinned' && ends !== 'fixed') throw new RangeError('ends must be "free", "pinned" or "fixed"');
  // nodal force the beam exerts on its restraints at a dof: (K_beam u − f)_dof, i.e. the reaction, evaluated
  // from equilibrium rather than from the tiny penalty displacement
  const restraintForce = (u, dof) => {
    let s = -f0[dof];
    for (let j = Math.max(0, dof - 3); j <= Math.min(nd - 1, dof + 3); j++) s += (j >= dof ? band[dof][j - dof] : band[j][dof - j]) * u[j];
    return s;
  };

  let u = null;
  let iterations = 0;
  for (; iterations < 200; iterations++) {
    const K = band.map((row) => Float64Array.from(row));
    const f = Float64Array.from(f0);
    if (ends !== 'free') {
      for (const node of [0, n]) {
        K[2 * node][0] += kRigid;
        if (ends === 'fixed') K[2 * node + 1][0] += kRigid * h * h;
      }
    }
    for (const s of sup) if (s.active) { K[2 * s.node][0] += s.k; f[2 * s.node] += s.k * s.height; }
    if (seabed) for (let i = 0; i < nn; i++) if (bedActive[i]) { const k = bedK * trib(i); K[2 * i][0] += k; f[2 * i] += k * bedZ; }
    try {
      u = solveBandedSym(K, f);
    } catch {
      throw new Error('beam is not sufficiently supported (mechanism): add supports, a seabed or end restraints');
    }
    for (let i = 0; i < nd; i++) if (!Number.isFinite(u[i]) || Math.abs(u[i]) > 1e6 * L) throw new Error('beam is not sufficiently supported (mechanism): add supports, a seabed or end restraints');
    let changed = false;
    for (const s of sup) {
      if (!s.unilateral) continue;
      const v = u[2 * s.node];
      if (s.active && s.k * (s.height - v) < 0) { s.active = false; changed = true; }
      else if (!s.active && v < s.height) { s.active = true; changed = true; }
    }
    if (seabed) {
      for (let i = 0; i < nn; i++) {
        const below = u[2 * i] < bedZ;
        if (bedActive[i] && !below) { bedActive[i] = 0; changed = true; } else if (!bedActive[i] && below) { bedActive[i] = 1; changed = true; }
      }
    }
    if (!changed) break;
  }
  if (iterations >= 200) throw new Error('contact iteration did not converge');

  const x = new Float64Array(nn);
  const deflection = new Float64Array(nn);
  const rotation = new Float64Array(nn);
  const moment = new Float64Array(nn);
  const shear = new Float64Array(nn);
  const seabedReaction = new Float64Array(nn);
  const cnt = new Uint8Array(nn);
  for (let i = 0; i < nn; i++) { x[i] = i * h; deflection[i] = u[2 * i]; rotation[i] = u[2 * i + 1]; }
  for (let e = 0; e < n; e++) {
    const v1 = u[2 * e];
    const t1 = u[2 * e + 1];
    const v2 = u[2 * e + 2];
    const t2 = u[2 * e + 3];
    // curvature from the shape functions plus the fixed-end particular solution of the uniform load
    const m1 = EI * ((-6 * v1) / (h * h) - (4 * t1) / h + (6 * v2) / (h * h) - (2 * t2) / h) + (q * h * h) / 12;
    const m2 = EI * ((6 * v1) / (h * h) + (2 * t1) / h - (6 * v2) / (h * h) + (4 * t2) / h) + (q * h * h) / 12;
    moment[e] += m1; cnt[e]++;
    moment[e + 1] += m2; cnt[e + 1]++;
    const vMid = (m2 - m1) / h;                    // dM/dx at mid-element
    shear[e] += vMid; shear[e + 1] += vMid;
  }
  let maxMoment = 0;
  let xMax = 0;
  // dV/dx = q, so the end-node shear is the mid-element value shifted by half an element
  shear[0] -= (q * h) / 2;
  shear[n] += (q * h) / 2;
  for (let i = 0; i < nn; i++) {
    moment[i] /= cnt[i]; shear[i] /= cnt[i];
    if (Math.abs(moment[i]) > maxMoment) { maxMoment = Math.abs(moment[i]); xMax = x[i]; }
  }
  const contact = new Uint8Array(nn);
  let total = 0;
  if (seabed) {
    for (let i = 0; i < nn; i++) if (bedActive[i]) {
      seabedReaction[i] = bedK * (bedZ - u[2 * i]);
      total += seabedReaction[i] * trib(i);
      contact[i] = 1;
    }
  }
  const reactions = sup.map((s) => {
    // a rigid support takes whatever the node needs beyond the seabed spring at the same node
    const R = !s.active ? 0 : s.k === kRigid ? restraintForce(u, 2 * s.node) - seabedReaction[s.node] * trib(s.node) : s.k * (s.height - u[2 * s.node]);
    total += R;
    if (s.active) contact[s.node] = 1;
    return { x: s.x, R, active: s.active };
  });
  if (ends !== 'free') {
    for (const node of [0, n]) {
      let R = restraintForce(u, 2 * node) - seabedReaction[node] * trib(node);
      for (const s of sup) if (s.active && s.node === node) R -= s.k * (s.height - u[2 * node]);
      total += R;
      contact[node] = 1;
    }
  }
  const spans = [];
  let last = -1;
  for (let i = 0; i < nn; i++) {
    if (!contact[i]) continue;
    if (last >= 0 && i - last > 1) spans.push({ from: x[last], to: x[i], length: x[i] - x[last] });
    last = i;
  }
  let totalLoad = w * L;
  for (const p of pointLoads) totalLoad += p.F;
  return { x, deflection, rotation, moment, shear, reactions, seabedReaction, spans, maxMoment, xMaxMoment: xMax,
    maxStress: D === undefined ? null : (maxMoment * pos('D', D) * E) / (2 * EI), totalReaction: total, totalLoad, iterations };
}

/**
 * Vertical separation between a crossing pipe and the crossed line, against the required minimum.
 * separation = (elevation of the crossing pipe axis − D_cross/2) − (top of the existing line incl. coating).
 *
 * Reference: DNV-ST-F101 (2021) §5 (crossings: minimum separation 0.3 m unless otherwise agreed);
 * the limit is an argument.
 *
 * @param {object} o
 * @param {number} o.axisElevation elevation of the crossing pipe axis at the crossing point [m]
 * @param {number} o.Dcross overall outer diameter of the crossing pipe [m]
 * @param {number} o.topOfExisting elevation of the top of the existing line [m] (same datum)
 * @param {number} [o.minSeparation=0.3] [m]
 * @param {number} [o.settlement=0] expected settlement of the crossing supports [m]
 * @returns {{separation:number, minSeparation:number, ok:boolean, margin:number}} [m]
 */
export function crossingSeparation({ axisElevation, Dcross, topOfExisting, minSeparation = 0.3, settlement = 0 }) {
  fin('axisElevation', axisElevation); pos('Dcross', Dcross); fin('topOfExisting', topOfExisting); nonNeg('minSeparation', minSeparation); nonNeg('settlement', settlement);
  const separation = axisElevation - Dcross / 2 - settlement - topOfExisting;
  return { separation, minSeparation, ok: separation >= minSeparation, margin: separation - minSeparation };
}

/**
 * Quantity take-off for a crossing: concrete mattresses and rock berms.
 *
 * Reference: geometry (rock via pipesoil.coverQuantities). Typical articulated concrete mattress:
 * 6 m × 3 m × 0.15–0.5 m, concrete density about 2400 kg/m³ — dimensions are arguments.
 *
 * @param {object} o
 * @param {Array<{count:number, length:number, width:number, thickness:number, density?:number, solidFraction?:number}>} [o.mattresses=[]]
 *   number, plan dimensions [m], thickness [m], concrete density [kg/m³] (default 2400) and solid fraction of
 *   the envelope volume (default 1)
 * @param {Array<object>} [o.rockBerms=[]] arguments of pipesoil.coverQuantities for each berm
 * @param {number} [o.rhoWater=1025] [kg/m³]
 * @returns {{mattressCount:number, mattressVolume:number, mattressMass:number, mattressSubmergedWeight:number,
 *   mattressArea:number, rockVolume:number, rockMass:number}} volumes [m³]; masses [kg]; weight [N]; area [m²]
 */
export function crossingQuantities({ mattresses = [], rockBerms = [], rhoWater = 1025 }) {
  let count = 0;
  let vol = 0;
  let mass = 0;
  let area = 0;
  let wsub = 0;
  mattresses.forEach((m, i) => {
    nonNeg(`mattresses[${i}].count`, m.count); pos(`mattresses[${i}].length`, m.length); pos(`mattresses[${i}].width`, m.width); pos(`mattresses[${i}].thickness`, m.thickness);
    const rho = m.density === undefined ? 2400 : pos(`mattresses[${i}].density`, m.density);
    const sf = m.solidFraction === undefined ? 1 : pos(`mattresses[${i}].solidFraction`, m.solidFraction);
    const v = m.count * m.length * m.width * m.thickness * sf;
    count += m.count; vol += v; mass += v * rho; area += m.count * m.length * m.width; wsub += v * (rho - rhoWater) * 9.80665;
  });
  let rockVolume = 0;
  let rockMass = 0;
  for (const b of rockBerms) {
    const r = coverQuantities(b);
    rockVolume += r.volume; rockMass += r.mass;
  }
  return { mattressCount: count, mattressVolume: vol, mattressMass: mass, mattressSubmergedWeight: wsub, mattressArea: area, rockVolume, rockMass };
}

/**
 * Dropped object falling through water: terminal velocity and impact energy.
 * v_T = √( 2 g (m − ρ_w V)/(ρ_w C_D A) ), E_T = ½ m v_T², effective impact energy E_E = ½ (m + m_a) v_T²,
 * m_a = ρ_w C_a V.
 *
 * Reference: DNV-RP-F107 (2019) §5 (object excursion and hit energy): force balance of weight, buoyancy and drag. Typical
 * coefficients there: slender shapes C_D 0.7–1.5, C_a 0.1–1.0; box shapes C_D 1.2–1.3, C_a 0.6–1.5;
 * complex (spherical to complex) shapes C_D 0.6–2.0, C_a 1.0–2.0 — arguments.
 *
 * @param {object} o
 * @param {number} o.mass [kg]
 * @param {number} o.volume displaced volume [m³]
 * @param {number} o.CD drag coefficient [-]
 * @param {number} o.area projected area in the direction of fall [m²]
 * @param {number} [o.Ca=0] added-mass coefficient [-]
 * @param {number} [o.rhoW=1025] [kg/m³]
 * @param {number} [o.g=9.80665]
 * @returns {{terminalVelocity:number, kineticEnergy:number, impactEnergy:number, addedMass:number}} [m/s], [J], [J], [kg]
 */
export function droppedObjectImpact({ mass, volume, CD, area, Ca = 0, rhoW = 1025, g = 9.80665 }) {
  pos('mass', mass); nonNeg('volume', volume); pos('CD', CD); pos('area', area); nonNeg('Ca', Ca); pos('rhoW', rhoW);
  const net = mass - rhoW * volume;
  if (!(net > 0)) throw new RangeError('the object is buoyant (mass must exceed rhoW·volume)');
  const vT = Math.sqrt((2 * g * net) / (rhoW * CD * area));
  const ma = rhoW * Ca * volume;
  return { terminalVelocity: vT, kineticEnergy: 0.5 * mass * vT * vT, impactEnergy: 0.5 * (mass + ma) * vT * vT, addedMass: ma };
}

/**
 * Kinetic energy of an impacting body with hydrodynamic added mass: E = ½ (m + m_a) v² (anchors, trawl
 * boards, clump weights).
 *
 * Reference: classical mechanics; format of DNV-RP-F111 (2014) §3 (the code's empirical shape and velocity
 * coefficients are applied by the caller through m_a and v).
 *
 * @param {object} o
 * @param {number} o.mass [kg]
 * @param {number} [o.addedMass=0] [kg]
 * @param {number} o.velocity [m/s]
 * @returns {number} energy [J]
 */
export function impactEnergy({ mass, addedMass = 0, velocity }) {
  pos('mass', mass); nonNeg('addedMass', addedMass); nonNeg('velocity', velocity);
  return 0.5 * (mass + addedMass) * velocity * velocity;
}

/**
 * Energy absorbed by a bare steel pipe in forming a dent of depth δ under a knife-edge load:
 * E = 16 (2π/9)^½ m_p (D/t)^½ D (δ/D)^(3/2), m_p = ¼ σ_y t².
 *
 * Reference: DNV-RP-F107 (2019) §4.3 (steel pipe capacity; after Wierzbicki & Suh 1988). Damage classes in the RP are stated in
 * terms of δ/D (e.g. < 5 % minor, > 20 % rupture likely).
 * Validity: empty pipe, knife-edge indenter transverse to the axis; conservative for larger contact areas.
 *
 * @param {object} o
 * @param {number} o.D steel outer diameter [m]
 * @param {number} o.t wall thickness [m]
 * @param {number} o.fy yield stress [Pa]
 * @param {number} [o.dentRatio] δ/D [-] (gives the energy)
 * @param {number} [o.energy] impact energy [J] (gives the dent)
 * @returns {{energy:(number|null), dentRatio:(number|null), mp:number}} mp plastic moment capacity of the wall [N]
 */
export function steelPipeDentEnergy({ D, t, fy, dentRatio, energy }) {
  pos('D', D); pos('t', t); pos('fy', fy);
  const mp = 0.25 * fy * t * t;
  const c = 16 * Math.sqrt((2 * Math.PI) / 9) * mp * Math.sqrt(D / t) * D;
  return { energy: dentRatio === undefined ? null : c * nonNeg('dentRatio', dentRatio) ** 1.5,
    dentRatio: energy === undefined ? null : (nonNeg('energy', energy) / c) ** (2 / 3), mp };
}

/**
 * Energy absorbed by crushing of a concrete weight coating.
 * 'edge' (indenter edge across the pipe): E_k = Y · b · h · x_0
 * 'flat' (indenter face of breadth b along the pipe): E_k = Y · b · (4/3) √(D x_0³)
 * Y crushing strength, b breadth of the impacting object, h depth of the indenting edge, x_0 penetration
 * (at most the coating thickness).
 *
 * Reference: DNV-RP-F107 (2019) §4.4 (concrete coating; as recalled from the RP: confirm the two expressions before use).
 * Y is given there as 3 to 5 times the cube strength for normal-density concrete; argument.
 *
 * @param {object} o
 * @param {'edge'|'flat'} [o.mode='edge']
 * @param {number} o.Y crushing strength [Pa]
 * @param {number} o.b breadth of the impacting object [m]
 * @param {number} [o.h] depth of the indenting edge [m] ('edge')
 * @param {number} [o.D] outer diameter over the coating [m] ('flat')
 * @param {number} o.x0 penetration [m]
 * @returns {number} absorbed energy [J]
 */
export function concreteCoatingEnergy({ mode = 'edge', Y, b, h, D, x0 }) {
  pos('Y', Y); pos('b', b); nonNeg('x0', x0);
  if (mode === 'edge') return Y * b * pos('h', h) * x0;
  if (mode === 'flat') return Y * b * (4 / 3) * Math.sqrt(pos('D', D) * x0 ** 3);
  throw new RangeError('mode must be "edge" or "flat"');
}

/**
 * Required protection against an impact: energy left for the steel after the protection layers have
 * absorbed their share, the resulting dent, and the acceptance against a dent limit.
 *
 * Reference: energy balance per DNV-RP-F107 (2019) §4 (capacities of concrete coating, polymer coating,
 * gravel cover etc. are supplied by the caller as energies).
 *
 * @param {object} o
 * @param {number} o.impactEnergy [J]
 * @param {number[]} [o.absorbed=[]] energy absorbed by each protection layer [J]
 * @param {number} o.D [m]
 * @param {number} o.t [m]
 * @param {number} o.fy [Pa]
 * @param {number} [o.dentLimit=0.05] acceptable δ/D [-]
 * @returns {{energyToSteel:number, dentRatio:number, ok:boolean, additionalProtection:number}}
 *   additionalProtection: energy still to be absorbed by extra protection to meet the dent limit [J]
 */
export function impactProtection({ impactEnergy: Ei, absorbed = [], D, t, fy, dentLimit = 0.05 }) {
  nonNeg('impactEnergy', Ei); pos('dentLimit', dentLimit);
  let a = 0;
  for (const e of absorbed) a += nonNeg('absorbed', e);
  const energyToSteel = Math.max(0, Ei - a);
  const dent = steelPipeDentEnergy({ D, t, fy, energy: energyToSteel }).dentRatio;
  const cap = steelPipeDentEnergy({ D, t, fy, dentRatio: dentLimit }).energy;
  return { energyToSteel, dentRatio: dent, ok: dent <= dentLimit, additionalProtection: Math.max(0, energyToSteel - cap) };
}

/**
 * Trawl pull-over: maximum horizontal force F_p = C_F · V · √(m_t · k_w) and the quasi-static response of a
 * pipe section treated as a beam of length L between points of fixity, loaded at mid-length, with the
 * lateral soil resistance q acting against the movement.
 *   fixed ends:  M_end = F L/8 − q L²/12, M_mid = F L/8 − q L²/24, δ = F L³/(192 EI) − q L⁴/(384 EI)
 *   pinned ends: M_mid = F L/4 − q L²/8, δ = F L³/(48 EI) − 5 q L⁴/(384 EI)
 *
 * Reference: force format of DNV-RP-F111 (2014) §4.3 (trawl boards). The empirical coefficient C_F depends
 * on board type and on the span height through the RP's expressions, and the warp stiffness k_w on wire
 * size and length; neither is built in — supply C_F and k_w (or the pull-over force directly).
 * Beam formulas: elementary statics. A dynamic or non-linear pipe–soil analysis is needed when the
 * quasi-static response is close to the limit.
 *
 * @param {object} o
 * @param {number} [o.Fp] pull-over force [N] (given directly), or computed from the next four
 * @param {number} [o.CF] empirical pull-over coefficient [-]
 * @param {number} [o.V] trawling velocity [m/s]
 * @param {number} [o.mass] steel mass of the board / beam / clump [kg]
 * @param {number} [o.kWarp] warp line stiffness [N/m]
 * @param {number} o.L length between points of fixity [m]
 * @param {number} o.EI [N·m²]
 * @param {number} [o.q=0] lateral soil resistance [N/m]
 * @param {'fixed'|'pinned'} [o.ends='fixed']
 * @param {number} [o.D] steel outer diameter [m] (for the stress)
 * @param {number} [o.E=207e9] [Pa]
 * @returns {{Fp:number, momentMax:number, deflection:number, stress:(number|null), soilHoldsLoad:boolean}}
 */
export function trawlPullover({ Fp, CF, V, mass, kWarp, L, EI, q = 0, ends = 'fixed', D, E = 207e9 }) {
  pos('L', L); pos('EI', EI); nonNeg('q', q);
  const F = Fp === undefined ? pos('CF', CF) * pos('V', V) * Math.sqrt(pos('mass', mass) * pos('kWarp', kWarp)) : pos('Fp', Fp);
  if (q * L >= F) return { Fp: F, momentMax: 0, deflection: 0, stress: D === undefined ? null : 0, soilHoldsLoad: true };
  let momentMax;
  let deflection;
  if (ends === 'fixed') {
    momentMax = Math.max(Math.abs((F * L) / 8 - (q * L * L) / 12), Math.abs((F * L) / 8 - (q * L * L) / 24));
    deflection = (F * L ** 3) / (192 * EI) - (q * L ** 4) / (384 * EI);
  } else if (ends === 'pinned') {
    momentMax = (F * L) / 4 - (q * L * L) / 8;
    deflection = (F * L ** 3) / (48 * EI) - (5 * q * L ** 4) / (384 * EI);
  } else throw new RangeError('ends must be "fixed" or "pinned"');
  return { Fp: F, momentMax, deflection, stress: D === undefined ? null : (momentMax * pos('D', D) * E) / (2 * EI), soilHoldsLoad: false };
}

/**
 * Rock armour size for hydraulic stability from the Shields criterion:
 * D_50 = τ_b/(θ_cr (ρ_s − ρ_w) g), with the bed shear stress τ_b = ½ ρ_w f U² when a velocity and friction
 * factor are given instead of τ_b.
 *
 * Reference: Shields (1936); CIRIA/CUR/CETMEF Rock Manual (2007) §5.2.1. The critical Shields parameter is
 * an argument: about 0.03–0.035 for no movement of rock armour, 0.05–0.055 for the onset of general
 * transport. Slope and wave–current interaction corrections are not applied (hook for a project method).
 *
 * @param {object} o
 * @param {number} o.thetaCr critical Shields parameter [-]
 * @param {number} [o.tau] bed shear stress [Pa]
 * @param {number} [o.U] near-bed velocity [m/s]
 * @param {number} [o.f] friction factor in τ_b = ½ ρ_w f U² [-]
 * @param {number} [o.rhoRock=2650] [kg/m³]
 * @param {number} [o.rhoW=1025] [kg/m³]
 * @param {number} [o.g=9.80665]
 * @returns {{D50:number, tau:number}} D50 [m]; tau [Pa]
 */
export function rockSizeShields({ thetaCr, tau, U, f, rhoRock = 2650, rhoW = 1025, g = 9.80665 }) {
  pos('thetaCr', thetaCr); pos('rhoRock', rhoRock); pos('rhoW', rhoW);
  if (rhoRock <= rhoW) throw new RangeError('rhoRock must exceed rhoW');
  const tb = tau === undefined ? 0.5 * rhoW * pos('f', f) * pos('U', U) ** 2 : nonNeg('tau', tau);
  return { D50: tb / (thetaCr * (rhoRock - rhoW) * g), tau: tb };
}
