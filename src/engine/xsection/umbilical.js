/**
 * Umbilical / dynamic-cable cross-section assembler and capacity checks (ISO 13628-5 / API Spec 17E scope).
 *
 * Model. Components (steel tubes, hoses, electrical cables, fibre units, fillers, armour wires) are placed as
 * helical elements on pitch-circle layers around an optional straight centre element, inside concentric sheaths.
 * Each element is a helical tendon with its own axial stiffness EA; the bundle is taken as radially rigid
 * (no diameter contraction), so with the kinematics of helix.js
 *     F_el = EA (cos^2(a) eps + lay R sin(a) cos(a) tau)
 *     K11 = sum EA cos^3 a,  K12 = sum lay EA R sin a cos^2 a,  K22 = sum EA R^2 sin^2 a cos a
 * plus the straight centre element and the sheaths. The elements' own torsion/bending stiffness is neglected in
 * K22 and included in the bending stiffness (slip state). Radial contraction, which softens real bundles with
 * polymer fillers, is not modelled: measured axial stiffness should be used to calibrate where available.
 */

import { helixBendingStiffness } from './helix.js';
import { lameStresses, vonMises } from './tubes.js';
import { minerDamage } from './flexfatigue.js';

const DEG = Math.PI / 180;

function pos(v, name) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function fin(v, name) {
  if (!(typeof v === 'number' && Number.isFinite(v))) throw new RangeError(`${name} must be a finite number`);
  return v;
}
function nonneg(v, name) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}

/**
 * Utilisation factors (allowable von Mises stress / yield) by load condition used as defaults for the capacity
 * curve. They follow the pattern of ISO 13628-5 / API 17E load-condition dependent utilisation (lower for normal
 * operation, higher for abnormal and accidental conditions); the numerical values are placeholders to be replaced
 * by those of the governing specification edition.
 */
export const UMBILICAL_UTILISATION = Object.freeze({ normal: 0.67, abnormal: 0.8, accidental: 1.0 });

/**
 * Normalise one component (tube, hose, cable, fibre unit, filler, wire).
 * Derived where not supplied: metal/solid area A = pi (od^2 - id^2)/4, EA = E A, EI = E pi (od^4 - id^4)/64,
 * GJ = G pi (od^4 - id^4)/32 with G = E/2.6, mass = rho A, bore area pi id^2/4. Composite items (cables, hoses)
 * should be given by their measured EA, EI, GJ and mass.
 * @param {object} c { kind, name, od [m], id [m], E [Pa], G [Pa], rho [kg/m^3], EA [N], EI [N m^2], GJ [N m^2],
 *   mass [kg/m], fluidDensity [kg/m^3], pressure [Pa], yieldStress [Pa] }
 * @returns {object} component with A [m^2], EA [N], EI, GJ [N m^2], mass [kg/m], boreArea [m^2], I [m^4]
 */
export function umbilicalComponent(c) {
  if (!c || typeof c !== 'object') throw new RangeError('component must be an object');
  const od = pos(c.od, 'component.od'), id = nonneg(c.id ?? 0, 'component.id');
  if (!(id < od)) throw new RangeError('component.id must be smaller than component.od');
  const A = Math.PI * (od * od - id * id) / 4, I = Math.PI * (od ** 4 - id ** 4) / 64;
  const E = c.E ?? 0;
  const G = c.G ?? E / 2.6;
  const out = {
    ...c, kind: c.kind ?? 'filler', name: c.name ?? c.kind ?? 'component', od, id, A, I, E,
    EA: c.EA ?? E * A, EI: c.EI ?? E * I, GJ: c.GJ ?? G * 2 * I,
    mass: c.mass ?? (c.rho ?? 0) * A, boreArea: Math.PI * id * id / 4, fluidDensity: c.fluidDensity ?? 0, pressure: c.pressure ?? 0,
  };
  for (const k of ['EA', 'EI', 'GJ', 'mass', 'fluidDensity', 'pressure']) nonneg(out[k], `component.${k}`);
  return out;
}

/**
 * Assemble an umbilical cross-section.
 * @param {object} def
 * @param {object} [def.centre] centre component (straight)
 * @param {{R:number, alphaDeg:number, lay?:number, elements:{component:object, count:number}[]}[]} def.layers
 *   helical layers inner to outer: pitch radius m, lay angle deg (0 < alpha < 90), hand +1/-1
 * @param {{ri:number, ro:number, E:number, nu?:number, rho?:number}[]} [def.sheaths] concentric polymer sheaths
 * @param {number} [def.outerDiameter] overall diameter, m (default: largest of sheath and layer envelopes)
 * @param {object} [o]
 * @param {number} [o.rhoSeawater=1025] kg/m^3
 * @param {number} [o.g=9.80665] m/s^2
 * @returns {object} { outerDiameter [m], massEmpty, massFilled [kg/m], weightAir, submergedEmpty,
 *   submergedFilled [N/m] (displacement based on the overall diameter: interstices taken as filled/sealed),
 *   K11 [N], K12 [N m], K22 [N m^2], EA (torque-free) [N], EAfixed [N], GJ [N m^2] (tension-free), GJfixed,
 *   EIstick, EIslip [N m^2], torquePerTension [m] (reaction torque per unit tension with twist restrained,
 *   K12/K11), twistPerTension [rad/m/N] (free end), torqueBalanced (|K12| < 1e-3 sqrt(K11 K22)),
 *   elements: flattened list [{ name, kind, layer, count, R, alpha, lay, component }] }
 */
export function assembleUmbilical(def, { rhoSeawater = 1025, g = 9.80665 } = {}) {
  if (!def || !Array.isArray(def.layers)) throw new RangeError('def.layers must be an array');
  pos(rhoSeawater, 'rhoSeawater'); pos(g, 'g');
  let K11 = 0, K12 = 0, K22 = 0, EIstick = 0, EIslip = 0, mE = 0, mF = 0, env = 0;
  const elements = [];
  if (def.centre) {
    const c = umbilicalComponent(def.centre);
    K11 += c.EA; K22 += c.GJ; EIstick += c.EI; EIslip += c.EI;
    mE += c.mass; mF += c.mass + c.fluidDensity * c.boreArea; env = Math.max(env, c.od);
    elements.push({ name: c.name, kind: c.kind, layer: -1, count: 1, R: 0, alpha: 0, lay: 1, component: c });
  }
  def.layers.forEach((L, li) => {
    const R = pos(L.R, `layers[${li}].R`);
    if (!(L.alphaDeg > 0 && L.alphaDeg < 90)) throw new RangeError(`layers[${li}].alphaDeg must be in (0, 90)`);
    const lay = L.lay ?? 1;
    if (lay !== 1 && lay !== -1) throw new RangeError(`layers[${li}].lay must be +1 or -1`);
    if (!Array.isArray(L.elements) || L.elements.length === 0) throw new RangeError(`layers[${li}].elements must be a non-empty array`);
    const a = L.alphaDeg * DEG, c = Math.cos(a), s = Math.sin(a);
    for (const el of L.elements) {
      const comp = umbilicalComponent(el.component);
      const n = pos(el.count, `layers[${li}].elements.count`);
      K11 += n * comp.EA * c ** 3;
      K12 += lay * n * comp.EA * R * s * c * c;
      K22 += n * comp.EA * R * R * s * s * c;
      if (comp.EA > 0) {
        const b = helixBendingStiffness({ n, E: 1, A: comp.EA, R, alpha: a });
        EIstick += b.stick;
      }
      // own bending/torsion of the element following the loxodromic path (helix.js energy expression, EI isotropic)
      const own = n / (2 * c) * (comp.EI * c ** 8 + comp.EI * (1 + s * s) ** 2 * c * c + comp.GJ * s * s * c ** 6);
      EIstick += own; EIslip += own;
      mE += n * comp.mass / c;
      mF += n * (comp.mass + comp.fluidDensity * comp.boreArea) / c;
      env = Math.max(env, 2 * R + comp.od);
      elements.push({ name: comp.name, kind: comp.kind, layer: li, count: n, R, alpha: a, lay, component: comp });
    }
  });
  for (const [i, sh] of (def.sheaths ?? []).entries()) {
    pos(sh.ri, `sheaths[${i}].ri`); pos(sh.ro, `sheaths[${i}].ro`); pos(sh.E, `sheaths[${i}].E`);
    if (!(sh.ro > sh.ri)) throw new RangeError(`sheaths[${i}].ro must exceed ri`);
    const A = Math.PI * (sh.ro ** 2 - sh.ri ** 2), I = Math.PI * (sh.ro ** 4 - sh.ri ** 4) / 4;
    K11 += sh.E * A; K22 += sh.E / (2 * (1 + (sh.nu ?? 0.4))) * 2 * I; EIstick += sh.E * I; EIslip += sh.E * I;
    mE += (sh.rho ?? 0) * A; mF += (sh.rho ?? 0) * A; env = Math.max(env, 2 * sh.ro);
  }
  const OD = def.outerDiameter ?? env;
  pos(OD, 'outerDiameter');
  const disp = rhoSeawater * Math.PI * OD * OD / 4;
  return {
    outerDiameter: OD, massEmpty: mE, massFilled: mF, weightAir: mE * g, submergedEmpty: (mE - disp) * g, submergedFilled: (mF - disp) * g,
    K11, K12, K22,
    EA: K22 > 0 ? K11 - K12 * K12 / K22 : K11, EAfixed: K11, GJ: K22 - K12 * K12 / K11, GJfixed: K22,
    EIstick, EIslip,
    torquePerTension: K12 / K11, twistPerTension: K22 > 0 ? -K12 / (K11 * K22 - K12 * K12) : 0,
    torqueBalanced: Math.abs(K12) < 1e-3 * Math.sqrt(K11 * K22),
    elements,
  };
}

/**
 * Sharing of tension (and torque) among the elements, and the resulting contact line loads.
 * (eps, tau) from the 2x2 stiffness: twist restrained (tau = 0) or torque-free. Element force
 * F = EA (cos^2 a eps + lay R sin a cos a tau). A tensioned helical element presses inward with the line load
 * w = F sin^2(a)/R per unit element length (helix.js). The contact line load on the inner side of an element is
 * its own w plus the radial load of all outer layers shared equally among the elements of its layer; on its outer
 * side only the latter:
 *     q_out,i = (sum_{j>i} n_j w_j / cos a_j) cos(a_i) / n_i,   q_in,i = w_i + q_out,i.
 * Reference: helix statics as derived in helix.js; ISO 13628-5 Annex (load sharing analysis) in intent.
 * Validity: radially rigid bundle, equal sharing within a layer, compressive element forces give no contact.
 * @param {object} section result of {@link assembleUmbilical}
 * @param {object} o
 * @param {number} o.tension N
 * @param {number} [o.torque=0] N m (used when twist is free)
 * @param {boolean} [o.twistRestrained=true]
 * @returns {{axialStrain:number, twist:number, torque:number, elements:{name:string, kind:string, layer:number,
 *   force:number, share:number, stress:number, lineLoad:number, contactInner:number, contactOuter:number}[]}}
 *   force N per element; share = fraction of the tension carried by all elements of that entry; stress = force
 *   over the component metal area, Pa; line loads N/m per unit element length
 */
export function umbilicalLoadSharing(section, { tension, torque = 0, twistRestrained = true }) {
  fin(tension, 'tension'); fin(torque, 'torque');
  const { K11, K12, K22 } = section;
  let eps, tau;
  if (twistRestrained || K22 === 0) { eps = tension / K11; tau = 0; } else {
    const det = K11 * K22 - K12 * K12;
    eps = (K22 * tension - K12 * torque) / det; tau = (K11 * torque - K12 * tension) / det;
  }
  const rows = section.elements.map((e) => {
    const c = Math.cos(e.alpha), s = Math.sin(e.alpha);
    const force = e.component.EA * (c * c * eps + e.lay * e.R * s * c * tau);
    const lineLoad = e.R > 0 ? Math.max(force, 0) * s * s / e.R : 0;
    return { name: e.name, kind: e.kind, layer: e.layer, force, share: tension !== 0 ? e.count * force * c / tension : 0, stress: e.component.A > 0 ? force / e.component.A : 0, lineLoad, contactInner: 0, contactOuter: 0 };
  });
  // radial load per metre of umbilical produced by each layer, then accumulate from the outside inwards
  const nLayers = section.elements.reduce((m, e) => Math.max(m, e.layer + 1), 0);
  const radial = new Array(nLayers).fill(0), wireLen = new Array(nLayers).fill(0);
  section.elements.forEach((e, i) => {
    if (e.layer < 0) return;
    radial[e.layer] += e.count * rows[i].lineLoad / Math.cos(e.alpha);
    wireLen[e.layer] += e.count / Math.cos(e.alpha);
  });
  let outer = 0;
  const fromOutside = new Array(nLayers).fill(0);
  for (let l = nLayers - 1; l >= 0; l--) { fromOutside[l] = outer; outer += radial[l]; }
  section.elements.forEach((e, i) => {
    if (e.layer < 0) return;
    rows[i].contactOuter = fromOutside[e.layer] / wireLen[e.layer];
    rows[i].contactInner = rows[i].lineLoad + rows[i].contactOuter;
  });
  return { axialStrain: eps, twist: tau, torque: twistRestrained ? K12 * eps : torque, elements: rows };
}

function elementUtilisation(e, row1, tension, kappa, mu, eta, pExt) {
  // row1: load-sharing row for unit tension (linear in tension)
  const c = e.component, cs = Math.cos(e.alpha), sn = Math.sin(e.alpha);
  const F = row1.force * tension;
  const f = mu * (row1.contactInner + row1.contactOuter) * Math.max(tension, 0);
  const stick = e.R > 0 ? (c.EA / c.A) * cs * cs * e.R * kappa : 0;
  const sf = e.R > 0 ? Math.min(stick, Math.PI * e.R * f / (2 * sn * c.A)) : 0;
  const ro = c.od / 2, ri = c.id / 2;
  const Emat = c.E || c.EA / c.A;
  const sb = Emat * ro * kappa;
  let vm;
  if (ri > 0) {
    vm = 0;
    for (const r of [ri, ro]) {
      const s = lameStresses({ pInt: c.pressure, pExt, ri, ro, r, axialForce: F });
      for (const sign of [1, -1]) vm = Math.max(vm, vonMises({ s1: s.radial, s2: s.hoop, s3: s.axial + sign * (sb * r / ro + sf) }).vonMises);
    }
  } else vm = Math.abs(F / c.A) + sb + sf;
  return { vm, utilisation: vm / (eta * c.yieldStress), frictionStress: sf, bendingStress: sb, axialStress: F / c.A };
}

/**
 * Utilisation of every stress-limited element (steel tubes: von Mises including internal pressure; armour wires
 * and strength members: axial stress) under combined tension and curvature.
 * Stress components: tension share F/A (plus closed-end pressure stresses by Lame, tubes.js); bending of the
 * element about its own axis E (d/2) kappa (the D/2R rule); friction stress from the stick-slip model of helix.js,
 * min(E cos^2 a R kappa, pi R f/(2 sin a A)) with f = mu (q_in + q_out). Only elements with `yieldStress` are
 * checked.
 * Reference: ISO 13628-5 §7 / API 17E (capacity based on von Mises stress of the most loaded tube); Saevik
 * (2011) for the friction stress.
 * @param {object} section assembled section
 * @param {object} o
 * @param {number} o.tension N
 * @param {number} o.curvature 1/m
 * @param {number} [o.mu=0.2] friction coefficient between elements (argument; no normative value)
 * @param {number} [o.eta=0.67] allowable utilisation factor
 * @param {number} [o.pExt=0] external pressure on the tubes, Pa
 * @param {boolean} [o.twistRestrained=true]
 * @returns {{utilisation:number, governing:string|null, elements:{name:string, vonMises:number,
 *   utilisation:number, axialStress:number, bendingStress:number, frictionStress:number}[]}} Pa, -
 */
export function umbilicalUtilisation(section, { tension, curvature, mu = 0.2, eta = 0.67, pExt = 0, twistRestrained = true }) {
  fin(tension, 'tension'); nonneg(curvature, 'curvature'); nonneg(mu, 'mu'); pos(eta, 'eta'); nonneg(pExt, 'pExt');
  const unit = umbilicalLoadSharing(section, { tension: 1, twistRestrained });
  let worst = 0, governing = null;
  const rows = [];
  section.elements.forEach((e, i) => {
    if (!e.component.yieldStress || !(e.component.A > 0)) return;
    const u = elementUtilisation(e, unit.elements[i], tension, curvature, mu, eta, pExt);
    rows.push({ name: e.name, vonMises: u.vm, utilisation: u.utilisation, axialStress: u.axialStress, bendingStress: u.bendingStress, frictionStress: u.frictionStress });
    if (u.utilisation > worst) { worst = u.utilisation; governing = e.name; }
  });
  if (rows.length === 0) throw new RangeError('no element carries a yieldStress: nothing to check');
  return { utilisation: worst, governing, elements: rows };
}

/**
 * Capacity curve: allowable tension versus curvature at a given utilisation factor (bisection on tension for
 * utilisation = 1 at each curvature), and the limiting curvature at zero tension.
 * Reference: ISO 13628-5 / API 17E capacity curves (normal, abnormal, accidental utilisation, see
 * {@link UMBILICAL_UTILISATION}).
 * @param {object} section assembled section
 * @param {object} [o]
 * @param {number} [o.eta=0.67]
 * @param {number} [o.mu=0.2]
 * @param {number} [o.pExt=0] Pa
 * @param {number} [o.nPoints=21] number of curvature points (>= 2)
 * @param {number} [o.kappaMax] largest curvature of the curve, 1/m (default: the zero-tension limit)
 * @param {boolean} [o.twistRestrained=true]
 * @returns {{curvature:number[], allowableTension:number[], governing:(string|null)[], kappaLimit:number,
 *   tensionLimit:number}} 1/m, N; kappaLimit = curvature at which the allowable tension vanishes; tensionLimit =
 *   allowable straight tension
 */
export function umbilicalCapacityCurve(section, { eta = 0.67, mu = 0.2, pExt = 0, nPoints = 21, kappaMax, twistRestrained = true } = {}) {
  if (!(Number.isInteger(nPoints) && nPoints >= 2)) throw new RangeError('nPoints must be an integer >= 2');
  const ut = (tension, curvature) => umbilicalUtilisation(section, { tension, curvature, mu, eta, pExt, twistRestrained });
  if (ut(0, 0).utilisation >= 1) throw new RangeError('section is over-utilised by pressure alone at this eta');
  const solveT = (kappa) => {
    if (ut(0, kappa).utilisation >= 1) return 0;
    let hi = 1e3;
    for (let i = 0; i < 60 && ut(hi, kappa).utilisation < 1; i++) hi *= 2;
    let lo = 0;
    for (let i = 0; i < 70; i++) { const mid = 0.5 * (lo + hi); if (ut(mid, kappa).utilisation < 1) lo = mid; else hi = mid; }
    return 0.5 * (lo + hi);
  };
  let kHi = 1e-3;
  for (let i = 0; i < 60 && ut(0, kHi).utilisation < 1; i++) kHi *= 2;
  let kLo = 0;
  for (let i = 0; i < 70; i++) { const mid = 0.5 * (kLo + kHi); if (ut(0, mid).utilisation < 1) kLo = mid; else kHi = mid; }
  const kappaLimit = 0.5 * (kLo + kHi);
  const kEnd = kappaMax === undefined ? kappaLimit : pos(kappaMax, 'kappaMax');
  const curvature = [], allowableTension = [], governing = [];
  for (let i = 0; i < nPoints; i++) {
    const k = kEnd * i / (nPoints - 1);
    const T = solveT(k);
    curvature.push(k); allowableTension.push(T); governing.push(ut(T, k).governing);
  }
  return { curvature, allowableTension, governing, kappaLimit, tensionLimit: allowableTension[0] };
}

/**
 * Minimum bend radius at a given tension: the bend radius at which the utilisation reaches 1 (bisection on
 * curvature). Storage MBR is evaluated at zero tension (and zero pressure if the components carry none) with the
 * storage utilisation; dynamic MBR at the operating tension with the normal-operation utilisation.
 * Reference: ISO 13628-5 / API 17E (MBR for storage, installation and operation derived from the capacity curve).
 * @param {object} section assembled section
 * @param {object} [o]
 * @param {number} [o.tension=0] N
 * @param {number} [o.eta=0.67]
 * @param {number} [o.mu=0.2]
 * @param {number} [o.pExt=0] Pa
 * @param {boolean} [o.twistRestrained=true]
 * @returns {{mbr:number, curvature:number}} m, 1/m (mbr = Infinity if no curvature is admissible)
 */
export function umbilicalMbr(section, { tension = 0, eta = 0.67, mu = 0.2, pExt = 0, twistRestrained = true } = {}) {
  const ut = (k) => umbilicalUtilisation(section, { tension, curvature: k, mu, eta, pExt, twistRestrained }).utilisation;
  if (ut(0) >= 1) return { mbr: Infinity, curvature: 0 };
  let hi = 1e-3;
  for (let i = 0; i < 60 && ut(hi) < 1; i++) hi *= 2;
  let lo = 0;
  for (let i = 0; i < 70; i++) { const mid = 0.5 * (lo + hi); if (ut(mid) < 1) lo = mid; else hi = mid; }
  const k = 0.5 * (lo + hi);
  return { mbr: 1 / k, curvature: k };
}

/**
 * Maximum bending moment in a thin circular ring loaded by n equal radial line loads P equally spaced round the
 * circumference (per unit length of tube). At the load points
 *     M = (P R / 2)(1/theta - cot(theta)),  theta = pi/n,
 * which for two opposing loads is M = P R/pi = 0.3183 P R exactly; between the loads
 * M = -(P R/2)(1/sin(theta) - 1/theta).
 * Reference: Roark's Formulas for Stress and Strain, Table 9.2 (ring under equal radial loads equally spaced);
 * Timoshenko, Strength of Materials Part I (ring under two opposing forces).
 * Validity: thin elastic ring (R/t > ~5), small deflection.
 * @param {object} o
 * @param {number} o.P line load at each point, N/m
 * @param {number} o.R mean ring radius, m
 * @param {number} [o.nLoads=2] number of equally spaced loads (>= 2)
 * @returns {{momentAtLoad:number, momentBetween:number, coefficient:number}} N m/m, N m/m, M/(P R)
 */
export function ringRadialLoadMoment({ P, R, nLoads = 2 }) {
  fin(P, 'P'); pos(R, 'R');
  if (!(Number.isInteger(nLoads) && nLoads >= 2)) throw new RangeError('nLoads must be an integer >= 2');
  const th = Math.PI / nLoads;
  const coefficient = 0.5 * (1 / th - 1 / Math.tan(th));
  return { momentAtLoad: coefficient * P * R, momentBetween: -0.5 * (1 / Math.sin(th) - 1 / th) * P * R, coefficient };
}

/**
 * Crush check of a steel tube inside the umbilical under caterpillar-tensioner squeeze or chute/sheave line load.
 *  - Tensioner with nTracks pads, each squeezing with q (N per metre of umbilical): by symmetry every track
 *    load passes radially through the bundle; the tube under a pad is squeezed between the pad side and the core
 *    by a line load P = q x (1/elementsUnderPad) x cos(alpha) per unit tube length (two opposing loads).
 *    `elementsUnderPad` is the number of outer-layer elements sharing one pad (1 = conservative: no spreading by
 *    sheath and fillers); the load distribution through the layers is this single documented factor.
 *  - Chute / sheave: line load T / Rchute on the umbilical, treated as a single pad load.
 *  - Tube capacity: ring bending stress 6 M/t^2 with M = 0.3183 P Rm ({@link ringRadialLoadMoment}), limited to
 *    allowableStress: P_allow = allowableStress t^2 / (6 x 0.3183 Rm).
 *  - Grip: friction capacity mu x nTracks x q x L must hold the tension: L_req = SF T/(mu nTracks q).
 * Reference: ring formula as above; ISO 13628-5 crush (lay-up and installation) checks in intent.
 * Validity: elastic first-yield ring criterion (conservative: a ductile tube tolerates local yielding), no
 * interaction with internal pressure or tension.
 * @param {object} o
 * @param {number} o.tubeOd m
 * @param {number} o.tubeWall m
 * @param {number} o.allowableStress Pa
 * @param {number} [o.padLoad] q, squeeze load per track per metre of umbilical, N/m
 * @param {number} [o.nTracks=2] 2, 3 or 4
 * @param {number} [o.elementsUnderPad=1]
 * @param {number} [o.alpha=0] lay angle of the tube, rad
 * @param {number} [o.tension] line tension, N
 * @param {number} [o.chuteRadius] m
 * @param {number} [o.mu] friction coefficient pad / sheath
 * @param {number} [o.safetyFactor=1] on the grip length
 * @returns {{tubeLineLoadCapacity:number, allowablePadLoad:number, tubeLineLoadTensioner:number|null,
 *   utilisationTensioner:number|null, chuteLineLoad:number|null, utilisationChute:number|null,
 *   requiredTensionerLength:number|null, minimumChuteRadius:number|null}} N/m, N/m, N/m, -, N/m, -, m, m
 */
export function umbilicalCrushCheck({ tubeOd, tubeWall, allowableStress, padLoad, nTracks = 2, elementsUnderPad = 1, alpha = 0, tension, chuteRadius, mu, safetyFactor = 1 }) {
  pos(tubeOd, 'tubeOd'); pos(tubeWall, 'tubeWall'); pos(allowableStress, 'allowableStress'); pos(elementsUnderPad, 'elementsUnderPad');
  if (!(tubeWall < tubeOd / 2)) throw new RangeError('tubeWall must be smaller than tubeOd/2');
  if (![2, 3, 4].includes(nTracks)) throw new RangeError('nTracks must be 2, 3 or 4');
  if (!(alpha >= 0 && alpha < Math.PI / 2)) throw new RangeError('alpha must be in [0, pi/2)');
  const Rm = (tubeOd - tubeWall) / 2;
  const cap = allowableStress * tubeWall ** 2 / (6 * ringRadialLoadMoment({ P: 1, R: Rm }).coefficient * Rm);
  const spread = Math.cos(alpha) / elementsUnderPad;
  const out = {
    tubeLineLoadCapacity: cap, allowablePadLoad: cap / spread, tubeLineLoadTensioner: null, utilisationTensioner: null,
    chuteLineLoad: null, utilisationChute: null, requiredTensionerLength: null, minimumChuteRadius: null,
  };
  if (padLoad !== undefined) {
    out.tubeLineLoadTensioner = pos(padLoad, 'padLoad') * spread;
    out.utilisationTensioner = out.tubeLineLoadTensioner / cap;
    if (tension !== undefined && mu !== undefined) out.requiredTensionerLength = pos(safetyFactor, 'safetyFactor') * nonneg(tension, 'tension') / (pos(mu, 'mu') * nTracks * padLoad);
  }
  if (tension !== undefined) {
    out.minimumChuteRadius = nonneg(tension, 'tension') * spread / cap;
    if (chuteRadius !== undefined) {
      out.chuteLineLoad = tension / pos(chuteRadius, 'chuteRadius');
      out.utilisationChute = out.chuteLineLoad * spread / cap;
    }
  }
  return out;
}

/**
 * Steady temperature of electrical cores inside the bundle by a series thermal-resistance network (per unit
 * length): each core has losses W_i and its own insulation resistance ln(d_ins/d_cond)/(2 pi k_ins); the total
 * heat then crosses the bundle (filler) resistance supplied by the caller, the sheath ln(ro/ri)/(2 pi k) and the
 * external film 1/(pi D h):
 *     T_core,i = T_amb + W_i R_ins,i + (sum W)(R_bundle + R_sheath + R_ext).
 * Reference: thermal-resistance method of IEC 60287 in its simplest form (series resistances); the
 * bundle resistance depends on lay-up and fill and must be estimated or measured (argument).
 * Validity: steady state, cores close together so that all losses cross the same outer resistances, no
 * solar/soil effects.
 * @param {object} o
 * @param {{loss:number, conductorDiameter:number, insulationDiameter:number, kInsulation:number}[]} o.cores W/m,
 *   m, m, W/m/K
 * @param {number} o.ambientC ambient temperature, degC
 * @param {number} [o.bundleResistance=0] K m/W
 * @param {{ri:number, ro:number, k:number}} [o.sheath] m, m, W/m/K
 * @param {number} [o.hExternal=Infinity] external film coefficient, W/m^2/K
 * @param {number} [o.outerDiameter] m (needed with hExternal)
 * @returns {{coreTempsC:number[], sheathInnerTempC:number, surfaceTempC:number, totalLoss:number,
 *   outerResistance:number}} degC, degC, degC, W/m, K m/W
 */
export function umbilicalCoreTemperature({ cores, ambientC, bundleResistance = 0, sheath, hExternal = Infinity, outerDiameter }) {
  if (!Array.isArray(cores) || cores.length === 0) throw new RangeError('cores must be a non-empty array');
  fin(ambientC, 'ambientC'); nonneg(bundleResistance, 'bundleResistance');
  const W = cores.reduce((a, c, i) => a + nonneg(c.loss, `cores[${i}].loss`), 0);
  const rExt = Number.isFinite(hExternal) ? 1 / (Math.PI * pos(outerDiameter, 'outerDiameter') * pos(hExternal, 'hExternal')) : 0;
  const rSh = sheath ? Math.log(pos(sheath.ro, 'sheath.ro') / pos(sheath.ri, 'sheath.ri')) / (2 * Math.PI * pos(sheath.k, 'sheath.k')) : 0;
  const surface = ambientC + W * rExt, shIn = surface + W * rSh, bundle = shIn + W * bundleResistance;
  return {
    coreTempsC: cores.map((c, i) => {
      if (!(c.insulationDiameter > c.conductorDiameter && c.conductorDiameter > 0)) throw new RangeError(`cores[${i}] diameters are inconsistent`);
      return bundle + c.loss * Math.log(c.insulationDiameter / c.conductorDiameter) / (2 * Math.PI * pos(c.kInsulation, `cores[${i}].kInsulation`));
    }),
    sheathInnerTempC: shIn, surfaceTempC: surface, totalLoss: W, outerResistance: rExt + rSh + bundleResistance,
  };
}

/**
 * Fatigue damage of one tube (or wire) element from tension and curvature ranges.
 * Stress range per bin = (element force per unit tension / A) dT + E (d/2) dKappa + friction range
 * min(E cos^2 a R dKappa, 2 pi R f/(2 sin a A)) with f = mu (q_in + q_out) at the bin mean tension; damage by
 * Miner's rule with the caller's S-N curve (e.g. a seam-welded super-duplex tube curve from qualification tests).
 * Reference: ISO 13628-5 fatigue analysis of steel tubes (S-N + Miner); flexfatigue.js `minerDamage`.
 * @param {object} section assembled section
 * @param {object} o
 * @param {string} o.element name of the element to assess
 * @param {{tensionRange:number, curvatureRange:number, cycles:number, meanTension?:number}[]} o.bins N, 1/m, -
 * @param {{a:number, m:number}} o.curve S-N curve N = a S^-m, S in Pa (bilinear/limit fields as in snCycles)
 * @param {number} [o.mu=0.2]
 * @param {number} [o.scf=1] stress concentration factor (welds)
 * @param {number} [o.dff=10] design fatigue factor
 * @param {number} [o.exposureYears=1]
 * @param {boolean} [o.twistRestrained=true]
 * @returns {object} minerDamage result plus ranges: [{ tensionPart, bendingPart, frictionPart, total }] in Pa
 */
export function umbilicalTubeFatigue(section, { element, bins, curve, mu = 0.2, scf = 1, dff = 10, exposureYears = 1, twistRestrained = true }) {
  const idx = section.elements.findIndex((e) => e.name === element);
  if (idx < 0) throw new RangeError(`element '${element}' not found in the section`);
  if (!Array.isArray(bins)) throw new RangeError('bins must be an array');
  const e = section.elements[idx], c = e.component;
  if (!(c.A > 0)) throw new RangeError('element has no metal area');
  const unit = umbilicalLoadSharing(section, { tension: 1, twistRestrained }).elements[idx];
  const cs = Math.cos(e.alpha), sn = Math.sin(e.alpha), Emat = c.E || c.EA / c.A;
  const ranges = [];
  const sBins = bins.map((b, i) => {
    nonneg(b.tensionRange, `bins[${i}].tensionRange`); nonneg(b.curvatureRange, `bins[${i}].curvatureRange`);
    const tensionPart = Math.abs(unit.force) * b.tensionRange / c.A;
    const bendingPart = Emat * c.od / 2 * b.curvatureRange;
    const f = mu * (unit.contactInner + unit.contactOuter) * Math.max(b.meanTension ?? 0, 0);
    const frictionPart = e.R > 0 ? Math.min((c.EA / c.A) * cs * cs * e.R * b.curvatureRange, 2 * Math.PI * e.R * f / (2 * sn * c.A)) : 0;
    const total = tensionPart + bendingPart + frictionPart;
    ranges.push({ tensionPart, bendingPart, frictionPart, total });
    return { stressRange: total, cycles: b.cycles };
  });
  return { ...minerDamage({ bins: sBins, curve, scf, dff, exposureYears }), ranges };
}

/**
 * Darcy-Weisbach pressure drop in a tube: dp = f (L/d) rho v^2/2, with f = 64/Re for Re < 2300 and the Haaland
 * explicit formula 1/sqrt(f) = -1.8 log10[(eps/(3.7 d))^1.11 + 6.9/Re] otherwise.
 * Reference: Haaland (1983) J. Fluids Eng. 105, 89-90; Darcy-Weisbach.
 * Validity: fully developed single-phase Newtonian flow; transition region 2300 < Re < 4000 treated as turbulent.
 * @param {object} o
 * @param {number} o.flow volumetric flow rate, m^3/s
 * @param {number} o.d bore diameter, m
 * @param {number} o.length m
 * @param {number} o.rho fluid density, kg/m^3
 * @param {number} o.viscosity dynamic viscosity, Pa s
 * @param {number} [o.roughness=0] absolute roughness, m
 * @returns {{dp:number, velocity:number, reynolds:number, frictionFactor:number, regime:string}} Pa, m/s, -, -
 */
export function tubePressureDrop({ flow, d, length, rho, viscosity, roughness = 0 }) {
  nonneg(flow, 'flow'); pos(d, 'd'); nonneg(length, 'length'); pos(rho, 'rho'); pos(viscosity, 'viscosity'); nonneg(roughness, 'roughness');
  const v = flow / (Math.PI * d * d / 4), Re = rho * v * d / viscosity;
  if (Re === 0) return { dp: 0, velocity: 0, reynolds: 0, frictionFactor: 0, regime: 'static' };
  const laminar = Re < 2300;
  const f = laminar ? 64 / Re : 1 / (1.8 * Math.log10((roughness / (3.7 * d)) ** 1.11 + 6.9 / Re)) ** 2;
  return { dp: f * length / d * rho * v * v / 2, velocity: v, reynolds: Re, frictionFactor: f, regime: laminar ? 'laminar' : 'turbulent' };
}

/**
 * Hydraulic response time of a long control line closed at the far end, linear (laminar) distributed model.
 * Per unit length: resistance r = 128 mu/(pi d^4) (Hagen-Poiseuille), capacitance c = A/K + ve (fluid
 * compressibility plus hose volumetric expansion ve in m^3/(m Pa)). The pressure obeys the diffusion equation
 * p_t = p_xx/(r c); for a step at the inlet and a closed end the far-end response is dominated by the first mode
 *     p_end(t)/p_step ~ 1 - (4/pi) exp(-t/tau1),  tau1 = (4/pi^2) r c L^2,
 * so the time to reach a fraction x is tau1 ln(4/(pi (1 - x))) (valid for x > ~0.5).
 * Reference: diffusion solution (Carslaw & Jaeger, Conduction of Heat in Solids, §3.3-3.4), applied to hydraulic
 * lines as in umbilical response-time analysis.
 * Validity: laminar flow, small pressure steps (constant properties), actuator volume not included.
 * @param {object} o
 * @param {number} o.length m
 * @param {number} o.d bore diameter, m
 * @param {number} o.viscosity Pa s
 * @param {number} o.bulkModulus fluid bulk modulus K, Pa
 * @param {number} [o.volumetricExpansion=0] ve, m^3/(m Pa)
 * @param {number} [o.fraction=0.9] target fraction of the step at the far end
 * @returns {{tau1:number, time:number, resistance:number, capacitance:number, storedVolumePerPa:number}} s, s,
 *   Pa s/m^3 (total), m^3/Pa (total), m^3/Pa
 */
export function hydraulicResponseTime({ length, d, viscosity, bulkModulus, volumetricExpansion = 0, fraction = 0.9 }) {
  pos(length, 'length'); pos(d, 'd'); pos(viscosity, 'viscosity'); pos(bulkModulus, 'bulkModulus'); nonneg(volumetricExpansion, 'volumetricExpansion');
  if (!(fraction >= 0.5 && fraction < 1)) throw new RangeError('fraction must be in [0.5, 1)');
  const r = 128 * viscosity / (Math.PI * d ** 4), c = Math.PI * d * d / 4 / bulkModulus + volumetricExpansion;
  const tau1 = 4 / Math.PI ** 2 * r * c * length * length;
  return { tau1, time: tau1 * Math.log(4 / (Math.PI * (1 - fraction))), resistance: r * length, capacitance: c * length, storedVolumePerPa: c * length };
}
