/**
 * Accidental impact: ship collision and dropped objects.
 *
 *  - Collision kinetic energy with hydrodynamic added mass and the strain energy to be dissipated for fixed and
 *    compliant installations (conservation of momentum; DNV-RP-C204 ch. 3, NORSOK N-004 Annex A.3).
 *  - Local denting of tubular members (Furnes & Amdahl 1980 form used in API RP 2A and DNV-RP-C204) and its energy.
 *  - Rigid-plastic three-hinge beam mechanism with membrane action for axially restrained members, derived from the
 *    normality rule for a thin-walled tube and a rectangular section (Søreide 1981; Haythornthwaite 1957; Jones,
 *    Structural Impact, 1989, ch. 7).
 *  - Energy sharing between striking and struck body from their force–deformation curves, design-category label.
 *  - Dropped objects: fall in air, fall through water with drag and added mass (terminal velocity, ODE), plastic
 *    energy absorption of a plate (yield-line mechanism + membrane), dent energy of a steel pipe (DNV-RP-F107).
 *  - Vessel impact on a single pile: elastic energy plus plastic hinge rotation.
 *
 * SI units (m, kg, s, N, J, Pa, rad). Plastic moment per unit width of a wall or plate: mp = fy t²/4.
 */
import { integrate } from '../math/quad.js';
import { brent } from '../math/roots.js';

const G0 = 9.80665;

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function nonNeg(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}

/* ------------------------------------------------------------------ collision energy */

/**
 * Kinetic energy of a vessel including hydrodynamic added mass: E = ½ (m + a) v² = ½ m Cm v².
 * The added-mass factor Cm = (m + a)/m is an input; DNV-RP-C204 ch. 3 and NORSOK N-003 quote a = 0.1 m for bow or stern
 * impact (Cm = 1.1, surge) and a = 0.4 m for broadside impact (Cm = 1.4, sway) of supply-type vessels.
 * @param {object} o
 * @param {number} o.mass vessel displacement mass m [kg]
 * @param {number} o.velocity impact speed v [m/s]
 * @param {number} o.addedMassFactor Cm = 1 + a/m [-] (>= 1; 1.1 surge, 1.4 sway as above)
 * @returns {{energy:number, addedMass:number, effectiveMass:number, momentum:number}} E [J], a [kg], m + a [kg], (m + a) v [N·s]
 */
export function collisionEnergy({ mass, velocity, addedMassFactor }) {
  pos('mass', mass); nonNeg('velocity', velocity);
  if (!(addedMassFactor >= 1)) throw new RangeError('addedMassFactor must be >= 1');
  const me = mass * addedMassFactor;
  return { energy: 0.5 * me * velocity * velocity, addedMass: me - mass, effectiveMass: me, momentum: me * velocity };
}

/**
 * Strain energy to be dissipated in a central collision. Both bodies move with a common velocity at the end of
 * the impact (perfectly plastic collision); conservation of momentum gives
 *   Es = ½ (ms + as) vs² (1 − vi/vs)² / (1 + (ms + as)/(mi + ai))   (compliant installation, as in DNV-RP-C204 ch. 3)
 *   Es = ½ (ms + as) vs²                                            (fixed installation, mi → ∞).
 * @param {object} o
 * @param {number} o.shipMass ms + as, ship mass including added mass [kg]
 * @param {number} o.shipVelocity vs [m/s]
 * @param {number} [o.installationMass=Infinity] mi + ai, installation mass including added mass [kg] (Infinity = fixed)
 * @param {number} [o.installationVelocity=0] vi, velocity of the installation in the direction of vs [m/s]
 * @returns {{strainEnergy:number, commonVelocity:number, kineticEnergyBefore:number}} Es [J], velocity after impact [m/s], total kinetic
 *   energy before impact [J]
 */
export function collisionStrainEnergy({ shipMass, shipVelocity, installationMass = Infinity, installationVelocity = 0 }) {
  pos('shipMass', shipMass);
  if (!Number.isFinite(shipVelocity) || !Number.isFinite(installationVelocity)) throw new RangeError('velocities must be finite');
  if (!(installationMass > 0)) throw new RangeError('installationMass must be positive');
  const dv = shipVelocity - installationVelocity;
  if (installationMass === Infinity) {
    return { strainEnergy: 0.5 * shipMass * dv * dv, commonVelocity: installationVelocity, kineticEnergyBefore: 0.5 * shipMass * shipVelocity ** 2 };
  }
  const vc = (shipMass * shipVelocity + installationMass * installationVelocity) / (shipMass + installationMass);
  return {
    strainEnergy: 0.5 * shipMass * dv * dv / (1 + shipMass / installationMass), commonVelocity: vc,
    kineticEnergyBefore: 0.5 * shipMass * shipVelocity ** 2 + 0.5 * installationMass * installationVelocity ** 2,
  };
}

/* ------------------------------------------------------------------ local denting of tubes */

/**
 * Resistance of a tubular member to local denting by a lateral concentrated load, and the energy absorbed:
 *   R/Rc = c (wd/D)^n,  Rc = fy t²/4 · √(D/t),  E = ∫R dwd = Rc D c/(n + 1) · (wd/D)^(n+1).
 * Defaults c = 22, n = 0.5. The same relation appears with other constants: API RP 2A-WSD commentary (Furnes & Amdahl)
 * Pd = 15 Mp (D/t)^½ (X/R)^½ with Mp = fy t²/4 and R = D/2 corresponds to c = 15√2 = 21.2, n = 0.5 (energy
 * 14.14 Mp X^1.5/t^0.5); DNV-RP-C204 (denting of tubular members) gives c = 22 + 1.2 B/D and n = 1.925/(3.5 + B/D) for a contact width B, with a
 * further reduction for axial load. Pass the constants of the code being applied.
 * Validity: the published relations are calibrated for dents up to roughly wd/D ≈ 0.5 and D/t of offshore tubulars.
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.t wall thickness [m]
 * @param {number} o.fy yield strength [Pa]
 * @param {number} o.dent dent depth wd [m] (0 ≤ wd < D)
 * @param {number} [o.coefficient=22] c [-]
 * @param {number} [o.exponent=0.5] n [-]
 * @returns {{Rc:number, R:number, energy:number, dentRatio:number}} characteristic resistance [N], denting force [N], absorbed energy [J], wd/D [-]
 */
export function tubeDentResistance({ D, t, fy, dent, coefficient = 22, exponent = 0.5 }) {
  pos('D', D); pos('t', t); pos('fy', fy); nonNeg('dent', dent); pos('coefficient', coefficient); pos('exponent', exponent);
  if (dent >= D) throw new RangeError('dent must be smaller than D');
  const Rc = fy * t * t / 4 * Math.sqrt(D / t);
  const x = dent / D;
  return { Rc, R: Rc * coefficient * x ** exponent, energy: Rc * D * coefficient / (exponent + 1) * x ** (exponent + 1), dentRatio: x };
}

/**
 * Dent depth of a tubular member that absorbs a given energy or resists a given force (inverse of tubeDentResistance).
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.t wall thickness [m]
 * @param {number} o.fy yield strength [Pa]
 * @param {number} [o.energy] energy absorbed by denting [J], or
 * @param {number} [o.force] lateral force [N]
 * @param {number} [o.coefficient=22] c [-]
 * @param {number} [o.exponent=0.5] n [-]
 * @returns {{dent:number, dentRatio:number, R:number, energy:number}} wd [m], wd/D [-], force [N] and energy [J] at that dent
 */
export function tubeDentDepth({ D, t, fy, energy, force, coefficient = 22, exponent = 0.5 }) {
  pos('D', D); pos('t', t); pos('fy', fy); pos('coefficient', coefficient); pos('exponent', exponent);
  const Rc = fy * t * t / 4 * Math.sqrt(D / t);
  let x;
  if (energy !== undefined) x = (nonNeg('energy', energy) * (exponent + 1) / (Rc * D * coefficient)) ** (1 / (exponent + 1));
  else if (force !== undefined) x = (nonNeg('force', force) / (Rc * coefficient)) ** (1 / exponent);
  else throw new RangeError('energy or force is required');
  return { dent: x * D, dentRatio: x, R: Rc * coefficient * x ** exponent, energy: Rc * D * coefficient / (exponent + 1) * x ** (exponent + 1) };
}

/* ------------------------------------------------------------------ three-hinge mechanism with membrane action */

/*
 * Rigid-plastic beam of span L loaded at midspan, ends fully restrained against axial displacement.
 * Kinematics for central deflection w (half span l = L/2): rotation rate of each half ẇ/l, elongation rate of the beam
 * 2wẇ/l. Clamped ends: hinges at both supports and at midspan — total rotation rate 4ẇ/l, so every hinge has
 * ė/θ̇ = w/2. Pinned ends: one hinge at midspan with rotation rate 2ẇ/l, so ė/θ̇ = w. Write ė/θ̇ = w/c1 with c1 = 2
 * (clamped) or 1 (pinned). Normality to the yield surface M = M(N) requires ė/θ̇ = −dM/dN, and the work equation gives
 *   R = (2 c1/l) (M + N w/c1) = R0 [M/Mp + (N/Mp) w/c1],   R0 = 4 c1 Mp/L.
 * Thin-walled tube, M/Mp = cos(πn/2), n = N/Np, Mp/Np = Dm/π:  −dM/dN = (Dm/2) sin(πn/2) = w/c1 → sin(πn/2) = w̄,
 *   w̄ = 2w/(c1 Dm):   R/R0 = √(1 − w̄²) + w̄ asin(w̄)  (w̄ ≤ 1),   R/R0 = (π/2) w̄  (w̄ > 1, pure membrane N = Np).
 * Rectangular section (depth H), M/Mp = 1 − n², Mp/Np = H/4:  −dM/dN = n H/2 = w/c1 → n = w̄, w̄ = 2w/(c1 H):
 *   R/R0 = 1 + w̄²  (w̄ ≤ 1),   R/R0 = 2 w̄  (w̄ > 1).
 * These are the classical results (clamped rectangular beam: P/Pc = 1 + (W/H)², 2W/H; the tube result is the curve given in DNV-RP-C204 for tubular beams with axial restraint).
 */
function mechanismRatio(type, wb) {
  if (type === 'tube') return wb <= 1 ? { r: Math.sqrt(1 - wb * wb) + wb * Math.asin(wb), n: 2 / Math.PI * Math.asin(wb) } : { r: Math.PI / 2 * wb, n: 1 };
  return wb <= 1 ? { r: 1 + wb * wb, n: wb } : { r: 2 * wb, n: 1 };
}

function mechanismSetup(o, fn) {
  const { section, fy, L, ends = 'clamped' } = o;
  pos('fy', fy); pos('L', L);
  if (ends !== 'clamped' && ends !== 'pinned') throw new RangeError(`${fn}: ends must be "clamped" or "pinned"`);
  const c1 = ends === 'clamped' ? 2 : 1;
  let Mp; let Np; let depth;
  if (section === 'tube') {
    pos('D', o.D); pos('t', o.t);
    if (o.t > o.D / 2) throw new RangeError('t must not exceed D/2');
    const Dm = o.D - o.t;
    Mp = fy * Dm * Dm * o.t; Np = fy * Math.PI * Dm * o.t; depth = Dm;
  } else if (section === 'rect') {
    pos('H', o.H); pos('B', o.B);
    Mp = fy * o.B * o.H * o.H / 4; Np = fy * o.B * o.H; depth = o.H;
  } else throw new RangeError(`${fn}: section must be "tube" or "rect"`);
  return { c1, Mp, Np, depth, R0: 4 * c1 * Mp / L };
}

/**
 * Load–deflection relation of an axially restrained beam deforming in a three-hinge plastic mechanism under a
 * midspan point load, including the membrane force that develops with deflection (rigid–perfectly-plastic material,
 * moderately large deflections, no local denting, no rupture limit). Derivation in the comment above:
 * R0 = 4 c1 Mp/L (c1 = 2 clamped, 1 pinned); thin-walled tube R/R0 = √(1 − w̄²) + w̄ asin w̄ for w̄ ≤ 1 and (π/2) w̄ beyond,
 * w̄ = 2w/(c1 (D − t)); rectangular section R/R0 = 1 + w̄² and 2w̄, w̄ = 2w/(c1 H).
 * Elastic axial flexibility of the supports delays the membrane action and is not included; use the bending value R0
 * when the ends can pull in.
 * @param {object} o
 * @param {'tube'|'rect'} o.section section type
 * @param {number} [o.D] tube outer diameter [m]
 * @param {number} [o.t] tube wall thickness [m]
 * @param {number} [o.H] depth of the rectangular section in the load direction [m]
 * @param {number} [o.B] width of the rectangular section [m]
 * @param {number} o.fy yield strength [Pa]
 * @param {number} o.L span [m]
 * @param {'clamped'|'pinned'} [o.ends='clamped'] rotational fixity of the (axially fixed) ends
 * @param {number} o.w midspan deflection [m]
 * @returns {{R0:number, R:number, ratio:number, wBar:number, axialRatio:number, Mp:number, Np:number, energy:number}} plastic collapse load in
 *   bending [N], resistance at w [N], R/R0 [-], w̄ [-], membrane force N/Np [-], plastic moment [N·m], plastic axial force [N], energy
 *   absorbed up to w [J]
 */
export function beamMechanismResistance(o) {
  const s = mechanismSetup(o, 'beamMechanismResistance');
  nonNeg('w', o.w);
  const scale = s.c1 * s.depth / 2;
  const wb = o.w / scale;
  const m = mechanismRatio(o.section, wb);
  // ∫ R dw = R0·scale·∫ ratio dw̄, split at the kink w̄ = 1
  const f = (x) => mechanismRatio(o.section, x).r;
  const e = wb <= 1 ? integrate(f, 0, wb, { tol: 1e-12 }) : integrate(f, 0, 1, { tol: 1e-12 }) + integrate(f, 1, wb, { tol: 1e-12 });
  return { R0: s.R0, R: s.R0 * m.r, ratio: m.r, wBar: wb, axialRatio: m.n, Mp: s.Mp, Np: s.Np, energy: s.R0 * scale * e };
}

/**
 * Midspan deflection at which the three-hinge mechanism of beamMechanismResistance has absorbed a given energy.
 * @param {object} o same section, fy, L, ends arguments as beamMechanismResistance, plus
 * @param {number} o.energy energy to absorb [J]
 * @returns {{w:number, R:number, ratio:number, wBar:number, axialRatio:number}} deflection [m], resistance [N], R/R0, w̄, N/Np
 */
export function beamMechanismDeflection(o) {
  const s = mechanismSetup(o, 'beamMechanismDeflection');
  nonNeg('energy', o.energy);
  if (o.energy === 0) return { w: 0, R: s.R0, ratio: 1, wBar: 0, axialRatio: 0 };
  const g = (w) => beamMechanismResistance({ ...o, w }).energy - o.energy;
  // the bending-only deflection energy/R0 is an upper bound because the resistance never falls below R0
  const w = brent(g, 0, o.energy / s.R0 * (1 + 1e-9), { tol: 1e-12 });
  const r = beamMechanismResistance({ ...o, w });
  return { w, R: r.R, ratio: r.ratio, wBar: r.wBar, axialRatio: r.axialRatio };
}

/* ------------------------------------------------------------------ energy sharing */

function curveCheck(name, c) {
  if (!c || !c.d || !c.F || c.d.length !== c.F.length || c.d.length < 2) throw new RangeError(`${name} must be { d: [...], F: [...] } with at least two points`);
  if (c.d[0] !== 0) throw new RangeError(`${name}.d must start at 0`);
  for (let i = 1; i < c.d.length; i++) {
    if (!(c.d[i] > c.d[i - 1])) throw new RangeError(`${name}.d must be strictly ascending`);
    if (c.F[i] < c.F[i - 1]) throw new RangeError(`${name}.F must be non-decreasing (use the envelope of a softening curve)`);
  }
}

/* Deformation and absorbed energy when a non-decreasing piecewise-linear curve first reaches force F (whole curve if F > Fmax). */
function upToForce(c, F) {
  let e = 0;
  const n = c.d.length;
  if (F <= c.F[0]) return { d: 0, e: 0, end: false };
  for (let i = 1; i < n; i++) {
    if (c.F[i] >= F) {
      const f = (F - c.F[i - 1]) / (c.F[i] - c.F[i - 1]);
      const d = c.d[i - 1] + f * (c.d[i] - c.d[i - 1]);
      return { d, e: e + 0.5 * (c.F[i - 1] + F) * (d - c.d[i - 1]), end: false };
    }
    e += 0.5 * (c.F[i] + c.F[i - 1]) * (c.d[i] - c.d[i - 1]);
  }
  return { d: c.d[n - 1], e, end: true };
}

/* Energy of a flat segment (constant force F, within tol) that starts at deformation d on curve c; 0 if there is none. */
function flatAt(c, d, F, tol) {
  for (let i = 1; i < c.d.length; i++) {
    if (Math.abs(c.d[i - 1] - d) <= 1e-9 * Math.max(1, c.d[c.d.length - 1]) && Math.abs(c.F[i] - F) <= tol && Math.abs(c.F[i - 1] - F) <= tol) {
      // consecutive flat segments at the same force are merged
      let j = i;
      while (j + 1 < c.d.length && Math.abs(c.F[j + 1] - F) <= tol) j++;
      return { length: c.d[j] - c.d[i - 1], force: c.F[i] };
    }
  }
  return null;
}

/**
 * Distribution of the collision strain energy between the striking ship and the struck structure from their
 * force–deformation curves: the bodies act in series, so both carry the same contact force and the energy is the sum
 * of the areas under the two curves up to that force (DNV-RP-C204 ch. 3; NORSOK N-004 Annex A.3). A body is rigid below
 * the first ordinate of its curve; when the force reaches a plateau of one curve, that body takes the further
 * deformation at constant force until the plateau ends (if both curves have a plateau at the same force the split is
 * indeterminate and the ship is deformed first).
 * Curves are piecewise linear with non-decreasing force (use the envelope of a softening curve); the last point of a
 * curve is its deformation limit.
 * @param {object} o
 * @param {{d:number[], F:number[]}} o.shipCurve ship force–indentation curve [m], [N], starting at d = 0
 * @param {{d:number[], F:number[]}} o.structureCurve structure force–deformation curve [m], [N]
 * @param {number} o.energy strain energy to dissipate [J]
 * @returns {{force:number, shipEnergy:number, structureEnergy:number, shipDeformation:number, structureDeformation:number,
 *   structureFraction:number, exceeded:boolean, capacity:number}} contact force at the end of impact [N], energies [J], deformations [m],
 *   structure share of the dissipated energy [-]; `exceeded` is true when the energy cannot be absorbed before the weaker curve ends
 *   (capacity [J] is the energy that can be absorbed up to that point)
 */
export function sharedEnergy({ shipCurve, structureCurve, energy }) {
  curveCheck('shipCurve', shipCurve); curveCheck('structureCurve', structureCurve); nonNeg('energy', energy);
  const fTop = Math.min(shipCurve.F[shipCurve.F.length - 1], structureCurve.F[structureCurve.F.length - 1]);
  const tol = 1e-7 * fTop;
  const total = (F) => upToForce(shipCurve, F).e + upToForce(structureCurve, F).e;
  // largest force whose "first reach" energy does not exceed the demand (total is non-decreasing, with jumps at plateaus)
  let lo = 0; let hi = fTop;
  if (total(fTop) <= energy) lo = fTop;
  else for (let it = 0; it < 200 && hi - lo > 1e-13 * fTop; it++) { const mid = 0.5 * (lo + hi); if (total(mid) <= energy) lo = mid; else hi = mid; }
  // state when the force first reaches F, then plateaus at F taken up to the remaining energy
  const walk = (F, E) => {
    const a = upToForce(shipCurve, F); const b = upToForce(structureCurve, F);
    const out = { force: F, shipEnergy: a.e, structureEnergy: b.e, shipDeformation: a.d, structureDeformation: b.d };
    let rest = E - a.e - b.e;
    for (const [c, eKey, dKey] of [[shipCurve, 'shipEnergy', 'shipDeformation'], [structureCurve, 'structureEnergy', 'structureDeformation']]) {
      if (!(rest > 0)) break;
      const fl = flatAt(c, out[dKey], F, tol);
      if (!fl) continue;
      const take = Math.min(rest, fl.force * fl.length);
      out[eKey] += take; out[dKey] += take / fl.force; rest -= take;
    }
    out.rest = rest;
    return out;
  };
  const out = walk(lo, energy);
  const absorbed = out.shipEnergy + out.structureEnergy;
  out.exceeded = lo >= fTop * (1 - 1e-12) && out.rest > 1e-9 * Math.max(1, energy);
  delete out.rest;
  out.structureFraction = absorbed > 0 ? out.structureEnergy / absorbed : 0;
  // capacity: everything that can be absorbed once the force has reached the lower of the two curve maxima
  const top = walk(fTop, Infinity);
  out.capacity = top.shipEnergy + top.structureEnergy;
  return out;
}

/**
 * Design-category label for a collision from the share of the strain energy dissipated by the installation
 * (terminology of NORSOK N-004 Annex A.3 and DNV-RP-C204): 'strength' — the ship dissipates most of the energy,
 * 'ductile' — the installation dissipates most of it, 'shared' otherwise. The standards define the categories
 * qualitatively; the numerical limits are project choices supplied as arguments.
 * @param {object} o
 * @param {number} o.structureEnergy energy dissipated by the installation [J]
 * @param {number} o.shipEnergy energy dissipated by the ship [J]
 * @param {number} [o.strengthLimit=0.1] structure share at or below which the design is labelled 'strength' [-]
 * @param {number} [o.ductileLimit=0.9] structure share at or above which the design is labelled 'ductile' [-]
 * @returns {{category:'strength'|'shared'|'ductile', structureFraction:number}}
 */
export function classifyCollisionDesign({ structureEnergy, shipEnergy, strengthLimit = 0.1, ductileLimit = 0.9 }) {
  nonNeg('structureEnergy', structureEnergy); nonNeg('shipEnergy', shipEnergy);
  if (!(strengthLimit >= 0 && ductileLimit <= 1 && strengthLimit < ductileLimit)) throw new RangeError('limits must satisfy 0 <= strengthLimit < ductileLimit <= 1');
  const sum = structureEnergy + shipEnergy;
  if (!(sum > 0)) throw new RangeError('structureEnergy + shipEnergy must be positive');
  const f = structureEnergy / sum;
  return { category: f <= strengthLimit ? 'strength' : (f >= ductileLimit ? 'ductile' : 'shared'), structureFraction: f };
}

/* ------------------------------------------------------------------ dropped objects */

/**
 * Dropped object falling in air (drag neglected): v = √(2 g h), E = m g h.
 * @param {object} o
 * @param {number} o.mass [kg]
 * @param {number} o.height drop height [m]
 * @param {number} [o.g=9.80665] [m/s²]
 * @returns {{velocity:number, energy:number, time:number}} impact velocity [m/s], kinetic energy [J], fall time [s]
 */
export function dropInAir({ mass, height, g = G0 }) {
  pos('mass', mass); nonNeg('height', height); pos('g', g);
  return { velocity: Math.sqrt(2 * g * height), energy: mass * g * height, time: Math.sqrt(2 * height / g) };
}

/**
 * Terminal velocity of an object sinking through water, from the balance of submerged weight and drag
 * (m − ρw V) g = ½ ρw Cd A vT², and the kinetic energy including added mass E = ½ (m + ma) vT², ma = ρw Ca V
 * (DNV-RP-F107 §5, hit energy of dropped objects). Cd and Ca are shape-dependent inputs.
 * @param {object} o
 * @param {number} o.mass [kg]
 * @param {number} o.volume displaced volume V [m³]
 * @param {number} o.Cd drag coefficient [-]
 * @param {number} o.area projected area in the direction of fall [m²]
 * @param {number} [o.Ca=0] added mass coefficient [-]
 * @param {number} [o.rhoW=1025] water density [kg/m³]
 * @param {number} [o.g=9.80665] [m/s²]
 * @returns {{velocity:number, energy:number, addedMass:number, submergedWeight:number}} vT [m/s] (0 for a buoyant object), E [J], ma [kg], W' [N]
 */
export function terminalVelocityWater({ mass, volume, Cd, area, Ca = 0, rhoW = 1025, g = G0 }) {
  pos('mass', mass); nonNeg('volume', volume); pos('Cd', Cd); pos('area', area); nonNeg('Ca', Ca); pos('rhoW', rhoW); pos('g', g);
  const W = (mass - rhoW * volume) * g;
  const ma = rhoW * Ca * volume;
  const v = W > 0 ? Math.sqrt(2 * W / (rhoW * Cd * area)) : 0;
  return { velocity: v, energy: 0.5 * (mass + ma) * v * v, addedMass: ma, submergedWeight: W };
}

/**
 * Fall of a dropped object through water from the surface to a given depth: (m + ma) dv/dt = (m − ρw V) g − ½ ρw Cd A v|v|,
 * dz/dt = v, integrated by the classical fourth-order Runge–Kutta method from the entry velocity (splash losses at the
 * surface are not modelled; pass a reduced `vEntry` to account for them). The velocity tends to the terminal velocity;
 * for constant coefficients v² = vT² + (v0² − vT²) exp(−ρw Cd A z/(m + ma)).
 * @param {object} o
 * @param {number} o.mass [kg]
 * @param {number} o.volume displaced volume [m³]
 * @param {number} o.Cd drag coefficient [-]
 * @param {number} o.area projected area [m²]
 * @param {number} [o.Ca=0] added mass coefficient [-]
 * @param {number} o.depth water depth to the target [m]
 * @param {number} [o.vEntry=0] velocity at the water surface [m/s] (e.g. dropInAir(...).velocity)
 * @param {number} [o.rhoW=1025] [kg/m³]
 * @param {number} [o.g=9.80665] [m/s²]
 * @param {number} [o.nOut=50] number of stored profile points
 * @returns {{reached:boolean, velocity:number, energy:number, time:number, terminalVelocity:number, maxDepth:number,
 *   profile:{t:number[], z:number[], v:number[]}}} whether the object reaches the depth (a buoyant object may stop), velocity [m/s] and
 *   kinetic energy ½(m + ma)v² [J] at the depth, travel time [s], terminal velocity [m/s], deepest point reached [m], sampled history
 */
export function dropThroughWater({ mass, volume, Cd, area, Ca = 0, depth, vEntry = 0, rhoW = 1025, g = G0, nOut = 50 }) {
  const tv = terminalVelocityWater({ mass, volume, Cd, area, Ca, rhoW, g });
  pos('depth', depth); nonNeg('vEntry', vEntry);
  const me = mass + tv.addedMass; const c = 0.5 * rhoW * Cd * area; const W = tv.submergedWeight;
  if (W <= 0 && vEntry === 0) return { reached: false, velocity: 0, energy: 0, time: 0, terminalVelocity: 0, maxDepth: 0, profile: { t: [0], z: [0], v: [0] } };
  const acc = (v) => (W - c * v * Math.abs(v)) / me;
  // time scales: velocity relaxation me/(c v) and transit depth/v
  const vRef = Math.max(vEntry, tv.velocity, Math.sqrt(2 * Math.abs(W) / me * depth) * 0.1, 1e-6);
  const dt = Math.min(me / (c * vRef), depth / vRef) / 400;
  let t = 0; let z = 0; let v = vEntry;
  const ts = [0]; const zs = [0]; const vs = [v];
  let reached = false; let steps = 0;
  const maxSteps = 5e6;
  while (steps++ < maxSteps) {
    const k1 = acc(v); const k2 = acc(v + 0.5 * dt * k1); const k3 = acc(v + 0.5 * dt * k2); const k4 = acc(v + dt * k3);
    const vn = v + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    const zn = z + dt * (v + dt / 6 * (k1 + k2 + k3));
    if (zn >= depth) {
      const f = (depth - z) / (zn - z);
      t += f * dt; v += f * (vn - v); z = depth; reached = true;
      break;
    }
    if (vn <= 0) { t += dt; z = zn; v = 0; break; }
    t += dt; z = zn; v = vn;
    ts.push(t); zs.push(z); vs.push(v);
  }
  ts.push(t); zs.push(z); vs.push(v);
  const stride = Math.max(1, Math.floor(ts.length / nOut));
  const pick = (a) => a.filter((_, i) => i % stride === 0 || i === a.length - 1);
  return {
    reached, velocity: reached ? v : 0, energy: reached ? 0.5 * me * v * v : 0, time: t, terminalVelocity: tv.velocity, maxDepth: z,
    profile: { t: pick(ts), z: pick(zs), v: pick(vs) },
  };
}

/**
 * Plastic resistance and energy absorption of a plate panel struck at its centre by a blunt object, idealised as a
 * circular rigid–perfectly-plastic plate of radius Rp loaded through a rigid punch of radius a.
 * Bending (yield-line / upper-bound mechanism): a conical field between the punch and the edge with a hinge circle at the
 * punch, and at the edge when clamped. Internal work for a punch displacement δ: hoop curvature 2π mp δ, hinge circle at the punch
 * 2π a mp δ/(Rp − a), edge hinge 2π Rp mp δ/(Rp − a); external work P δ, hence
 *   Pc = 2π mp Rp/(Rp − a)  (simply supported),   Pc = 4π mp Rp/(Rp − a)  (clamped),   mp = fy t²/4
 * (a → 0 recovers the classical 2π mp and 4π mp).
 * Membrane: with the radial membrane force at its plastic value N0 = fy t and small slopes, vertical equilibrium of the plate
 * outside the punch gives dw/dr = P/(2π N0 r), i.e. Pm(δ) = 2π N0 δ / ln(Rp/a) for edges restrained in the plane.
 * The resistance is taken as the envelope max(Pc, Pm(δ)) — bending first, membrane once it exceeds the bending load — and
 * the energy as its integral. No rupture criterion is included: `deflection` must be limited by the caller (strain or
 * punching limits), and a rectangular panel is represented by its inscribed circle.
 * @param {object} o
 * @param {number} o.t plate thickness [m]
 * @param {number} o.fy yield strength [Pa]
 * @param {number} o.radius plate radius Rp [m] (half the shorter panel side)
 * @param {number} o.punchRadius contact radius a [m] (0 < a < Rp)
 * @param {'clamped'|'simple'} [o.edge='clamped'] rotational edge condition
 * @param {boolean} [o.membrane=true] include membrane action (edges restrained against pull-in)
 * @param {number} o.deflection punch displacement δ [m]
 * @returns {{Pc:number, membraneLoad:number, resistance:number, energy:number, transitionDeflection:number, mp:number}} bending collapse load [N],
 *   membrane load at δ [N], resistance at δ [N], absorbed energy [J], deflection at which membrane action takes over [m], mp [N·m/m]
 */
export function plateImpactCapacity({ t, fy, radius, punchRadius, edge = 'clamped', membrane = true, deflection }) {
  pos('t', t); pos('fy', fy); pos('radius', radius); pos('punchRadius', punchRadius); nonNeg('deflection', deflection);
  if (punchRadius >= radius) throw new RangeError('punchRadius must be smaller than radius');
  if (edge !== 'clamped' && edge !== 'simple') throw new RangeError('edge must be "clamped" or "simple"');
  const mp = fy * t * t / 4;
  const Pc = (edge === 'clamped' ? 4 : 2) * Math.PI * mp * radius / (radius - punchRadius);
  const km = membrane ? 2 * Math.PI * fy * t / Math.log(radius / punchRadius) : 0;
  const wt = membrane ? Pc / km : Infinity;
  const Pm = km * deflection;
  const energy = deflection <= wt ? Pc * deflection : Pc * wt + 0.5 * km * (deflection * deflection - wt * wt);
  return { Pc, membraneLoad: Pm, resistance: Math.max(Pc, Pm), energy, transitionDeflection: wt, mp };
}

/**
 * Energy absorbed by a dent in a steel pipe wall struck by a knife-edge load perpendicular to the pipe
 * (DNV-RP-F107 §4, steel pipe dent–energy relation): E = 16 (2π/9)^½ mp (D/t)^½ D (δ/D)^(3/2), mp = fy t²/4. Coating and concrete
 * contributions are not included.
 * @param {object} o
 * @param {number} o.D steel outer diameter [m]
 * @param {number} o.t wall thickness [m]
 * @param {number} o.fy yield stress [Pa]
 * @param {number} o.dent dent depth δ [m]
 * @returns {{energy:number, dentRatio:number, mp:number}} E [J], δ/D [-], mp [N·m/m]
 */
export function pipeDentEnergy({ D, t, fy, dent }) {
  pos('D', D); pos('t', t); pos('fy', fy); nonNeg('dent', dent);
  if (dent >= D) throw new RangeError('dent must be smaller than D');
  const mp = fy * t * t / 4;
  return { energy: 16 * Math.sqrt(2 * Math.PI / 9) * mp * Math.sqrt(D / t) * D * (dent / D) ** 1.5, dentRatio: dent / D, mp };
}

/**
 * Dent depth of a steel pipe for a given absorbed impact energy (inverse of pipeDentEnergy).
 * @param {object} o
 * @param {number} o.D steel outer diameter [m]
 * @param {number} o.t wall thickness [m]
 * @param {number} o.fy yield stress [Pa]
 * @param {number} o.energy absorbed energy [J]
 * @returns {{dent:number, dentRatio:number}} δ [m], δ/D [-] (a ratio above 1 means the relation is outside its range)
 */
export function pipeDentFromEnergy({ D, t, fy, energy }) {
  pos('D', D); pos('t', t); pos('fy', fy); nonNeg('energy', energy);
  const mp = fy * t * t / 4;
  const x = (energy / (16 * Math.sqrt(2 * Math.PI / 9) * mp * Math.sqrt(D / t) * D)) ** (2 / 3);
  return { dent: x * D, dentRatio: x };
}

/* ------------------------------------------------------------------ vessel impact on a pile */

/**
 * Energy absorption of a single free-head pile (or dolphin pile) struck at a height h above its point of fixity.
 * Elastic stage (cantilever): stiffness 3EI/h³, hinge forms at the fixity point when F h = Mp, elastic energy
 * Ee = ½ Fc δe = Mp² h/(6 EI). Plastic stage: the hinge rotates at constant moment, Epl = Mp θp. The demand is compared with
 * the capacity Ee + Mp θcap, where the plastic rotation capacity θcap depends on section class, D/t and soil support and
 * is supplied by the caller. Shear deformation, soil flexibility above the fixity point and P–Δ effects are neglected.
 * @param {object} o
 * @param {number} o.energy impact energy to be absorbed by the pile [J]
 * @param {number} o.Mp plastic moment of the pile section [N·m]
 * @param {number} o.EI bending stiffness [N·m²]
 * @param {number} o.height height of the impact point above the point of fixity [m]
 * @param {number} o.rotationCapacity plastic hinge rotation capacity θcap [rad]
 * @returns {{collapseLoad:number, elasticDeflection:number, elasticEnergy:number, plasticRotation:number, deflection:number, force:number,
 *   capacity:number, uc:number, elastic:boolean}} hinge-forming load Mp/h [N], deflection at first hinge [m], elastic energy capacity [J],
 *   required plastic rotation [rad], deflection of the impact point at the end of impact [m], peak contact force [N], energy capacity
 *   Ee + Mp θcap [J], utilisation energy/capacity [-], and whether the pile remains elastic
 */
export function pileImpact({ energy, Mp, EI, height, rotationCapacity }) {
  nonNeg('energy', energy); pos('Mp', Mp); pos('EI', EI); pos('height', height); nonNeg('rotationCapacity', rotationCapacity);
  const Fc = Mp / height;
  const k = 3 * EI / height ** 3;
  const de = Fc / k;
  const Ee = 0.5 * Fc * de;
  const capacity = Ee + Mp * rotationCapacity;
  if (energy <= Ee) {
    const d = Math.sqrt(2 * energy / k);
    return { collapseLoad: Fc, elasticDeflection: de, elasticEnergy: Ee, plasticRotation: 0, deflection: d, force: k * d, capacity, uc: energy / capacity, elastic: true };
  }
  const th = (energy - Ee) / Mp;
  return {
    collapseLoad: Fc, elasticDeflection: de, elasticEnergy: Ee, plasticRotation: th, deflection: de + th * height, force: Fc, capacity,
    uc: energy / capacity, elastic: false,
  };
}
