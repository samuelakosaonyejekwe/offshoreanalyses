// Pile drivability: static resistance to driving (SRD), Smith one-dimensional wave equation, blow counts.
// Residual stresses between blows are not modelled (each blow starts from an unstressed pile).

import { unitResistanceAt } from './axialcapacity.js';
import { propertiesAt } from './soilprofile.js';

const G = 9.80665;

function pos(name, v, allowZero = false) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < 0 || (!allowZero && v === 0)) {
    throw new RangeError(`${name}: must be a finite number ${allowZero ? '>= 0' : '> 0'} (got ${v})`);
  }
}

/**
 * Generic hydraulic impact hammers for demonstration. These are NOT manufacturer data: rated energy and ram
 * mass are round numbers chosen so that rated energy = ram weight × ~2.04 m equivalent stroke, and the
 * efficiency is an editable default. Replace with the data sheet of the actual hammer.
 * Fields: ratedEnergy [J], ramMass [kg], efficiency [-] (energy at impact / rated energy).
 */
export const GENERIC_HAMMERS = Object.freeze({
  'generic-500': Object.freeze({ label: 'Generic hydraulic hammer 500 kJ', ratedEnergy: 500e3, ramMass: 25e3, efficiency: 0.9 }),
  'generic-800': Object.freeze({ label: 'Generic hydraulic hammer 800 kJ', ratedEnergy: 800e3, ramMass: 40e3, efficiency: 0.9 }),
  'generic-1200': Object.freeze({ label: 'Generic hydraulic hammer 1200 kJ', ratedEnergy: 1200e3, ramMass: 60e3, efficiency: 0.9 }),
  'generic-2000': Object.freeze({ label: 'Generic hydraulic hammer 2000 kJ', ratedEnergy: 2000e3, ramMass: 100e3, efficiency: 0.9 }),
  'generic-3000': Object.freeze({ label: 'Generic hydraulic hammer 3000 kJ', ratedEnergy: 3000e3, ramMass: 150e3, efficiency: 0.9 }),
  'generic-4000': Object.freeze({ label: 'Generic hydraulic hammer 4000 kJ', ratedEnergy: 4000e3, ramMass: 200e3, efficiency: 0.9 }),
});

/**
 * Ram impact velocity from the delivered kinetic energy: v0 = sqrt(2 η E / m).
 * Reference: energy balance. Validity: η lumps all losses before impact.
 * @param {{ ratedEnergy: number, ramMass: number, efficiency?: number }} hammer ratedEnergy [J]; ramMass [kg];
 *   efficiency [-] in (0, 1] (default 1)
 * @returns {number} impact velocity [m/s]
 */
export function impactVelocity({ ratedEnergy, ramMass, efficiency = 1 }) {
  pos('ratedEnergy', ratedEnergy); pos('ramMass', ramMass);
  if (!(efficiency > 0 && efficiency <= 1)) throw new RangeError('efficiency: must be in (0, 1]');
  return Math.sqrt((2 * efficiency * ratedEnergy) / ramMass);
}

/**
 * One blow by the Smith wave-equation model: ram (rigid mass) – cushion spring – helmet mass – pile discretised
 * into lumped masses and springs, with elastic–perfectly-plastic soil springs (quake q, ultimate Ru) and Smith
 * damping along the shaft and at the toe (toe in compression only). The helmet mass is
 * lumped on the pile head. With a cushion the ram acts through a compression-only spring; without one the ram
 * is in rigid unilateral contact with the pile head (common velocity from momentum conservation at impact, so a
 * fraction m_head/(M + m_head) of the impact energy is lost, vanishing with mesh refinement).
 * Smith's expression is applied as R = Rs + J |Rs| v, identical for Rs >= 0 and dissipative during rebound.
 * Explicit (leapfrog) time integration of the pile chain with the step dtFactor times the stability limit ΔL/c
 * (reduced where a cushion spring is stiffer). Soil springs and dampers are integrated semi-implicitly and do
 * not restrict the step. At dtFactor = 1 a uniform chain carries a wave front exactly, but the scheme is only
 * marginally stable and can diverge with stiff toe resistance; the default 0.9 is robust. Below the limit the
 * lumped model disperses a sharp front: for a rigid ram striking without a cushion the peak segment stress is
 * overestimated by up to about 20 % (the analytic peak is ρ c v0), while the head force, the transferred energy
 * and the blow count are not affected appreciably. A cushion removes the sharp front and the overshoot.
 * pileDamping adds stiffness-proportional internal damping c = η k ΔL/c (default 0), which suppresses that
 * overshoot at the price of mesh-dependent energy loss.
 * Set = maximum toe displacement − toe quake (Smith's definition); blow count = 1/set.
 * Reference: Smith (1960) "Bearing capacity of piles by the wave equation", J. Soil Mech. Found. Div. ASCE
 * 86(SM4); Smith's recommended parameters q = 2.54 mm, J(shaft) = 0.16 s/m, J(toe) = 0.49 s/m (0.1 in,
 * 0.05 and 0.15 s/ft) are the defaults. Validity: one-dimensional stress waves in a prismatic or stepped pile;
 * no residual stress analysis; no hammer assembly dynamics beyond a rigid ram.
 * @param {object} o
 *   hammer: { ramMass [kg], ratedEnergy [J], efficiency [-] } or { ramMass, impactVelocity [m/s] };
 *   cushionStiffness [N/m] (optional); helmetMass [kg] (optional);
 *   pile: { length [m], area [m²] (or D, t [m]), E [Pa] (210e9), rho [kg/m³] (7850), nSegments [-] (100) };
 *   shaftRu: number[] ultimate static shaft resistance per pile segment from the head down [N] (default zeros);
 *   toeRu [N] (default 0); shaftQuake, toeQuake [m] (default 0.00254), number or (shaft) array per segment;
 *   shaftDamping, toeDamping [s/m] (defaults 0.164, 0.492), number or (shaft) array per segment;
 *   dtFactor [-] in (0, 1] (default 0.9); pileDamping η [-] (default 0); duration [s] (default 6 L/c + 0.02); history [bool] (default false)
 * @returns {object} { set [m] permanent set per blow; blowsPerMetre [1/m] (Infinity at refusal);
 *   maxCompressiveStress, maxTensileStress [Pa] in the pile (both >= 0); maxToeDisplacement [m];
 *   enthru [J] maximum energy transferred to the pile head; impactEnergy [J] ram kinetic energy at impact;
 *   energyDrift [-] max |E(t) − E0| / E0 of kinetic + strain energy (meaningful without soil and while
 *   the ram stays in contact or on its cushion); waveSpeed [m/s];
 *   dt [s]; steps [-]; impedance [N·s/m] = ρ c A; history { t [s], headForce [N], headVelocity [m/s],
 *   toeDisplacement [m] } when requested }
 */
export function waveEquation(o) {
  const { hammer, pile } = o;
  if (!hammer || !pile) throw new RangeError('hammer: hammer and pile are required');
  pos('hammer.ramMass', hammer.ramMass);
  const v0 = hammer.impactVelocity ?? impactVelocity(hammer);
  pos('hammer.impactVelocity', v0);
  const E = pile.E ?? 210e9, rho = pile.rho ?? 7850, nSeg = pile.nSegments ?? 100;
  pos('pile.length', pile.length); pos('pile.E', E); pos('pile.rho', rho);
  if (!(Number.isInteger(nSeg) && nSeg >= 2)) throw new RangeError('pile.nSegments: must be an integer >= 2');
  let area = pile.area;
  if (area === undefined) { pos('pile.D', pile.D); pos('pile.t', pile.t); area = (Math.PI * (pile.D ** 2 - (pile.D - 2 * pile.t) ** 2)) / 4; }
  pos('pile.area', area);
  const dL = pile.length / nSeg, c = Math.sqrt(E / rho);
  const mSeg = rho * area * dL, kSeg = (E * area) / dL;
  const dtFactor = o.dtFactor ?? 0.9;
  if (!(dtFactor > 0 && dtFactor <= 1)) throw new RangeError('dtFactor: must be in (0, 1]');
  const arr = (v, def, name) => {
    const a = new Float64Array(nSeg);
    for (let i = 0; i < nSeg; i++) { a[i] = Array.isArray(v) || ArrayBuffer.isView(v) ? v[i] : (v ?? def); pos(name, a[i], true); }
    return a;
  };
  const Ru = arr(o.shaftRu ?? 0, 0, 'shaftRu');
  if ((Array.isArray(o.shaftRu) || ArrayBuffer.isView(o.shaftRu)) && o.shaftRu.length !== nSeg) throw new RangeError('shaftRu: needs one value per pile segment');
  const qS = arr(o.shaftQuake, 0.00254, 'shaftQuake'), jS = arr(o.shaftDamping, 0.164, 'shaftDamping');
  const toeRu = o.toeRu ?? 0, qT = o.toeQuake ?? 0.00254, jT = o.toeDamping ?? 0.492;
  pos('toeRu', toeRu, true); pos('toeQuake', qT); pos('toeDamping', jT, true);
  for (let i = 0; i < nSeg; i++) if (Ru[i] > 0 && !(qS[i] > 0)) throw new RangeError('shaftQuake: must be > 0');

  // Pile nodes at the segment ends (head = node 0, toe = node nSeg) with half masses at the two ends, so that
  // the chain has the uniform stability limit ΔL/c and propagates a wave front without numerical dispersion
  // when the time step is close to that limit.
  const helmet = o.helmetMass ?? 0;
  pos('helmetMass', helmet, true);
  const kc = o.cushionStiffness ?? null;
  if (kc !== null) pos('cushionStiffness', kc);
  const n = nSeg + 1, M = hammer.ramMass;
  const mass = new Float64Array(n).fill(mSeg);
  mass[0] = 0.5 * mSeg + helmet; mass[n - 1] = 0.5 * mSeg;
  const RuN = new Float64Array(n), qN = new Float64Array(n), jN = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    RuN[i] = 0.5 * ((i > 0 ? Ru[i - 1] : 0) + (i < nSeg ? Ru[i] : 0));
    const j = Math.min(i, nSeg - 1);
    qN[i] = qS[j]; jN[i] = jS[j];
  }
  // Soil springs and Smith damping are integrated semi-implicitly (they only touch the diagonal), so the step
  // is governed by the pile chain (ΔL/c) and, where present, by the cushion spring.
  let dtCrit = dL / c;
  if (kc !== null) dtCrit = Math.min(dtCrit, Math.sqrt((2 * mass[0]) / (kSeg + kc)), Math.sqrt((2 * M) / kc));
  const eta = o.pileDamping ?? 0;
  pos('pileDamping', eta, true);
  // stiffness-proportional damping c = η k ΔL/c shortens the stable step by the factor sqrt(1 + η²) − η
  dtCrit *= Math.sqrt(1 + eta * eta) - eta;
  const dt = dtFactor * dtCrit;
  const cPile = (eta * kSeg * dL) / c;
  const duration = o.duration ?? (6 * pile.length) / c + 0.02;
  pos('duration', duration);
  const steps = Math.ceil(duration / dt);

  const u = new Float64Array(n), v = new Float64Array(n), vOld = new Float64Array(n), f = new Float64Array(n);
  const up = new Float64Array(n); // plastic slip of the shaft springs
  let upToe = 0, uRam = 0, vRam = v0, contact = kc === null;
  const E0 = 0.5 * M * v0 * v0;
  // rigid contact (no cushion): ram and pile head move together from the instant of impact (momentum conserved)
  if (contact) { v[0] = (M * v0) / (M + mass[0]); vRam = v[0]; }
  let maxC = 0, maxT = 0, maxToe = 0, enthru = 0, work = 0, drift = 0;
  const hist = o.history ? { t: [], headForce: [], headVelocity: [], toeDisplacement: [] } : null;
  // v⁺ = (v + dt (F − Rs)/m) / (1 + dt c/m + dt² k/m): static soil force Rs now, tangent k and damper c implicit
  const advance = (vi, F, m, Rs, k, cd) => (vi + (dt * (F - Rs)) / m) / (1 + (dt * cd) / m + (dt * dt * k) / m);
  for (let step = 0; step < steps; step++) {
    f.fill(0);
    let strain = 0;
    for (let sgm = 0; sgm < nSeg; sgm++) {
      const el = kSeg * (u[sgm] - u[sgm + 1]); // compression positive
      const force = el + cPile * (v[sgm] - v[sgm + 1]);
      f[sgm] -= force; f[sgm + 1] += force;
      strain += (0.5 * el * el) / kSeg;
      if (force > maxC) maxC = force; else if (-force > maxT) maxT = -force;
    }
    let Fc = 0;
    if (kc !== null) {
      Fc = Math.max(kc * (uRam - u[0]), 0);
      strain += (0.5 * Fc * Fc) / kc;
      f[0] += Fc;
    } else if (!contact && uRam >= u[0] && vRam > v[0]) {
      contact = true; v[0] = (M * vRam + mass[0] * v[0]) / (M + mass[0]); vRam = v[0]; uRam = u[0];
    }
    vOld.set(v);
    const vRamOld = vRam;
    for (let i = 0; i < n; i++) {
      let Rs = 0, k = 0, cd = 0;
      if (RuN[i] > 0) {
        const kq = RuN[i] / qN[i];
        let rs = kq * (u[i] - up[i]);
        if (rs > RuN[i]) { rs = RuN[i]; up[i] = u[i] - qN[i]; } else if (rs < -RuN[i]) { rs = -RuN[i]; up[i] = u[i] + qN[i]; }
        Rs = rs; cd = jN[i] * Math.abs(rs);
        // no further stiffness while the spring keeps yielding in the direction of motion
        if (!(Math.abs(rs) >= RuN[i] && rs * vOld[i] > 0)) k = kq;
      }
      let RsT = 0, kT = 0, cT = 0;
      if (i === n - 1 && toeRu > 0) {
        const kq = toeRu / qT;
        let rs = kq * (u[i] - upToe);
        if (rs > toeRu) { rs = toeRu; upToe = u[i] - qT; } else if (rs < 0) { rs = 0; upToe = u[i]; }
        RsT = rs; cT = jT * rs;
        if (!(rs >= toeRu && vOld[i] > 0) && !(rs === 0 && vOld[i] <= 0)) kT = kq;
      }
      const m = i === 0 && contact ? mass[0] + M : mass[i];
      let vn = advance(vOld[i], f[i], m, Rs + RsT, k + kT, cd + cT);
      // Smith: the toe cannot pull on the pile, so the toe term is dropped when 1 + J v < 0
      if (RsT > 0 && 1 + jT * vn < 0) vn = advance(vOld[i], f[i], m, Rs, k, cd);
      if (i === 0 && contact) {
        Fc = (-M * (vn - vOld[0])) / dt; // force decelerating the ram
        if (Fc < 0) { Fc = 0; contact = false; vn = advance(vOld[0], f[0], mass[0], Rs, k, cd); } else vRam = vn;
      }
      v[i] = vn;
    }
    if (kc !== null) vRam -= (Fc / M) * dt;
    if (Fc > maxC) maxC = Fc;
    // kinetic energy at the time of the displacements: the product of the two staggered velocities is the
    // quantity conserved by the leapfrog scheme
    let kin = 0.5 * M * vRamOld * vRam;
    for (let i = 0; i < n; i++) kin += 0.5 * mass[i] * vOld[i] * v[i];
    const d = Math.abs(kin + strain - E0) / E0;
    if (d > drift) drift = d;
    const vHead = 0.5 * (vOld[0] + v[0]);
    work += Fc * vHead * dt;
    if (work > enthru) enthru = work;
    uRam += vRam * dt;
    for (let i = 0; i < n; i++) u[i] += v[i] * dt;
    if (u[n - 1] > maxToe) maxToe = u[n - 1];
    if (hist) { hist.t.push(step * dt); hist.headForce.push(Fc); hist.headVelocity.push(vHead); hist.toeDisplacement.push(u[n - 1]); }
  }
  const set = Math.max(maxToe - qT, 0);
  const out = {
    set, blowsPerMetre: set > 0 ? 1 / set : Infinity, maxCompressiveStress: maxC / area, maxTensileStress: maxT / area,
    maxToeDisplacement: maxToe, enthru, impactEnergy: E0, energyDrift: drift, waveSpeed: c, dt, steps, impedance: rho * c * area,
  };
  if (hist) out.history = hist;
  return out;
}

/**
 * Analytic force at the head of a semi-infinite elastic rod struck by a rigid mass: F(t) = Z v0 exp(−Z t / M),
 * Z = ρ c A, c = sqrt(E/ρ); peak stress ρ c v0. Valid until the reflection from the far end returns (t < 2L/c).
 * Reference: Timoshenko & Goodier, Theory of Elasticity (longitudinal impact of bars); Smith (1960) benchmark.
 * @param {{ ramMass: number, impactVelocity: number, area: number, E?: number, rho?: number, t: number }} o
 *   ramMass [kg]; impactVelocity [m/s]; area [m²]; E [Pa] (210e9); rho [kg/m³] (7850); t [s] time after impact
 * @returns {{ force: number, stress: number, peakStress: number, impedance: number }} [N], [Pa], [Pa], [N·s/m]
 */
export function rodImpactForce({ ramMass, impactVelocity: v0, area, E = 210e9, rho = 7850, t }) {
  pos('ramMass', ramMass); pos('impactVelocity', v0); pos('area', area); pos('E', E); pos('rho', rho); pos('t', t, true);
  const Z = rho * Math.sqrt(E / rho) * area;
  const force = Z * v0 * Math.exp((-Z * t) / ramMass);
  return { force, stress: force / area, peakStress: rho * Math.sqrt(E / rho) * v0, impedance: Z };
}

/**
 * Pile capacity factor for unit skin friction in clay during continuous driving: Fp = 0.5 (OCR)^0.3.
 * Reference: Semple & Gemeinhardt (1981) OTC 4148, as used by Stevens, Wiltsie & Turton (1982) OTC 4205.
 * Validity: clays; OCR >= 1.
 * @param {number} OCR overconsolidation ratio [-]
 * @returns {number} Fp [-]
 */
export function clayPileCapacityFactor(OCR) {
  if (!(OCR >= 1)) throw new RangeError('OCR: must be >= 1');
  return 0.5 * Math.pow(OCR, 0.3);
}

const SRD_CASES = ['coringLower', 'coringUpper', 'pluggedLower', 'pluggedUpper'];

/**
 * Static resistance to driving (SRD) versus depth, after Stevens, Wiltsie & Turton (1982), built on the API-type
 * unit resistances of axialcapacity.js:
 *   clay: unit friction = Fp · f_static with Fp = 0.5 OCR^0.3 (layer.OCR, default 1), unit end bearing 9 su;
 *   sand (and other layers): static unit friction and end bearing with their limits (Fp = 1).
 *   coringLower: outer friction + 50 % of that unit friction on the inner wall + end bearing on the steel annulus;
 *   coringUpper: outer + full inner friction + end bearing on the annulus;
 *   pluggedLower: outer friction + end bearing on the gross area;
 *   pluggedUpper: as pluggedLower with unit friction and end bearing multiplied by the upperPlugged factors.
 * The upper-bound plugged factors are arguments. Defaults: sand friction × 1.3 and sand end bearing × 1.5 as
 * commonly quoted from Stevens et al.; clay factors default to 1.0 and must be set by the user (the original
 * clay upper bound is not reproduced here).
 * Reference: Stevens, Wiltsie & Turton (1982) "Evaluating pile drivability for hard clay, very dense sand and
 * rock", OTC 4205. Validity: continuous driving of open-ended pipe piles; no set-up after delays.
 * @param {{ profile: object, D: number, t: number, Lmax?: number, dz?: number, upperPlugged?: object,
 *   frictionFactorAt?: Function }} o profile; D [m]; t [m]; Lmax [m] (default profile depth); dz [m] depth step
 *   (default 0.5); upperPlugged { sandFriction (1.3), sandBearing (1.5), clayFriction (1), clayBearing (1) } [-];
 *   frictionFactorAt(z, props) → factor [-] overrides the pile capacity factor on unit friction
 * @returns {object} { z [m] depths (0 … Lmax); fUnit [Pa] driving unit friction at each depth; qUnit [Pa] unit
 *   end bearing; isClay [bool[]]; cases: { coringLower | coringUpper | pluggedLower | pluggedUpper:
 *   { total [N][], shaft [N][], toe [N][], shaftPerLength [N/m][] (resistance per metre of pile at each depth
 *   for that case) } } }
 */
export function srdProfile(o) {
  const { profile, D, t } = o;
  pos('D', D); pos('t', t);
  if (t >= D / 2) throw new RangeError('t: must be less than D/2');
  const Lmax = o.Lmax ?? profile.depth, dz = o.dz ?? 0.5;
  pos('Lmax', Lmax); pos('dz', dz);
  if (Lmax > profile.depth + 1e-9) throw new RangeError('Lmax: exceeds the profile depth');
  const up = { sandFriction: 1.3, sandBearing: 1.5, clayFriction: 1, clayBearing: 1, ...(o.upperPlugged ?? {}) };
  for (const k of Object.keys(up)) pos(`upperPlugged.${k}`, up[k]);
  const n = Math.ceil(Lmax / dz - 1e-9) + 1;
  const Di = D - 2 * t, Aann = (Math.PI * (D * D - Di * Di)) / 4, Agross = (Math.PI * D * D) / 4;
  const z = [], fUnit = [], qUnit = [], isClay = [];
  for (let i = 0; i < n; i++) {
    // sample just inside the layer below an interface so that the step in resistance is captured
    const zz = Math.min((i * Lmax) / (n - 1), Lmax);
    const r = unitResistanceAt(profile, zz);
    const pr = propertiesAt(profile, zz);
    const clay = r.type === 'clay';
    const Fp = o.frictionFactorAt ? o.frictionFactorAt(zz, pr) : clay ? clayPileCapacityFactor(pr.OCR ?? 1) : 1;
    z.push(zz); fUnit.push(r.f * Fp); qUnit.push(r.q); isClay.push(clay);
  }
  const perim = {
    coringLower: () => Math.PI * D + 0.5 * Math.PI * Di,
    coringUpper: () => Math.PI * (D + Di),
    pluggedLower: () => Math.PI * D,
    pluggedUpper: (i) => Math.PI * D * (isClay[i] ? up.clayFriction : up.sandFriction),
  };
  const toeOf = {
    coringLower: (i) => qUnit[i] * Aann,
    coringUpper: (i) => qUnit[i] * Aann,
    pluggedLower: (i) => qUnit[i] * Agross,
    pluggedUpper: (i) => qUnit[i] * Agross * (isClay[i] ? up.clayBearing : up.sandBearing),
  };
  const cases = {};
  for (const name of SRD_CASES) {
    const shaftPerLength = z.map((_, i) => fUnit[i] * perim[name](i));
    const shaft = [0];
    for (let i = 1; i < n; i++) shaft.push(shaft[i - 1] + 0.5 * (shaftPerLength[i - 1] + shaftPerLength[i]) * (z[i] - z[i - 1]));
    const toe = z.map((_, i) => toeOf[name](i));
    cases[name] = { total: shaft.map((s, i) => s + toe[i]), shaft, toe, shaftPerLength };
  }
  return { z, fUnit, qUnit, isClay, cases };
}

function interpAt(xs, ys, x) {
  if (x <= xs[0]) return ys[0];
  const n = xs.length;
  if (x >= xs[n - 1]) return ys[n - 1];
  let i = 1;
  while (xs[i] < x) i++;
  return ys[i - 1] + ((ys[i] - ys[i - 1]) * (x - xs[i - 1])) / (xs[i] - xs[i - 1]);
}

/**
 * Bearing graph: blow count, driving stresses and transferred energy for a range of total SRD values, the SRD
 * being split between shaft (distributed uniformly or triangularly over the embedded length) and toe.
 * Reference: Smith (1960); see waveEquation. Validity: as waveEquation.
 * @param {object} o waveEquation options (hammer, pile, cushionStiffness, helmetMass, quakes, damping …) plus
 *   srd: number[] total static resistances [N]; shaftFraction [-] in [0, 1] (default 0.7); penetration [m]
 *   embedded length (default pile length); distribution 'uniform'|'triangular' (default 'triangular',
 *   increasing with depth)
 * @returns {{ srd: number[], blowsPerMetre: number[], blowsPer025m: number[], set: number[],
 *   maxCompressiveStress: number[], maxTensileStress: number[], enthru: number[] }} SRD [N]; blow counts
 *   [1/m], [1/0.25 m]; set [m]; stresses [Pa]; ENTHRU [J]
 */
export function bearingGraph(o) {
  const { srd, shaftFraction = 0.7, distribution = 'triangular' } = o;
  if (!Array.isArray(srd) || srd.length === 0) throw new RangeError('srd: must be a non-empty array');
  if (!(shaftFraction >= 0 && shaftFraction <= 1)) throw new RangeError('shaftFraction: must be within 0..1');
  if (distribution !== 'uniform' && distribution !== 'triangular') throw new RangeError("distribution: must be 'uniform' or 'triangular'");
  const nSeg = o.pile.nSegments ?? 100, Lp = o.pile.length, pen = o.penetration ?? Lp;
  pos('penetration', pen);
  if (pen > Lp) throw new RangeError('penetration: cannot exceed the pile length');
  const dL = Lp / nSeg, w = new Float64Array(nSeg);
  let sum = 0;
  for (let j = 0; j < nSeg; j++) {
    const a = Math.max(j * dL - (Lp - pen), 0), b = Math.max((j + 1) * dL - (Lp - pen), 0); // embedded depth range of the segment
    w[j] = distribution === 'uniform' ? b - a : 0.5 * (b * b - a * a);
    sum += w[j];
  }
  const out = { srd: [], blowsPerMetre: [], blowsPer025m: [], set: [], maxCompressiveStress: [], maxTensileStress: [], enthru: [] };
  for (const R of srd) {
    pos('srd', R, true);
    const r = waveEquation({ ...o, shaftRu: Array.from(w, (x) => (R * shaftFraction * x) / sum), toeRu: R * (1 - shaftFraction), history: false });
    out.srd.push(R); out.blowsPerMetre.push(r.blowsPerMetre); out.blowsPer025m.push(r.blowsPerMetre / 4); out.set.push(r.set);
    out.maxCompressiveStress.push(r.maxCompressiveStress); out.maxTensileStress.push(r.maxTensileStress); out.enthru.push(r.enthru);
  }
  return out;
}

/**
 * Drivability study: for each penetration depth the SRD of the chosen case is distributed along the embedded
 * pile segments, one blow is analysed with the wave equation, and the blow count, driving stresses and ENTHRU
 * are reported; blows are integrated over depth, refusal is flagged against a configurable limit, and a list of
 * (stress range, number of blows) blocks is produced for a driving-fatigue assessment (stress range taken as
 * the maximum compressive plus maximum tensile driving stress anywhere in the pile — an upper estimate).
 * Soil damping may differ between clay and sand segments; the Smith values are the defaults for both.
 * Reference: Smith (1960); Stevens et al. (1982) for the SRD; see waveEquation and srdProfile.
 * Validity: continuous driving, no set-up, no residual stresses; refusal criteria are contractual — the default
 * 250 blows per 0.25 m is an example and must be confirmed for the project.
 * @param {object} o profile; hammer; pile { length [m], D [m], t [m], E, rho, nSegments }; cushionStiffness;
 *   helmetMass; srdCase 'coringLower'|'coringUpper'|'pluggedLower'|'pluggedUpper' (default 'coringUpper');
 *   depths: number[] penetrations [m] (default every 1 m from 1 m to min(pile length, profile depth));
 *   upperPlugged, frictionFactorAt as srdProfile; quake { shaft, toe } [m]; damping { shaftClay, shaftSand,
 *   toeClay, toeSand } [s/m]; refusalBlowsPer025m [-] (default 250); dtFactor
 * @returns {object} { depth [m][], srd [N][], blowsPerMetre [1/m][], blowsPer025m [][], maxCompressiveStress
 *   [Pa][], maxTensileStress [Pa][], enthru [J][], refusal [bool][], refusalDepth [m]|null, totalBlows [-]
 *   (to the last depth or to refusal), fatigueBlocks [{ depthFrom [m], depthTo [m], stressRange [Pa],
 *   blows [-] }] }
 */
export function drivabilityAnalysis(o) {
  const { profile, hammer, pile } = o;
  const srdCase = o.srdCase ?? 'coringUpper';
  if (!SRD_CASES.includes(srdCase)) throw new RangeError(`srdCase: must be one of ${SRD_CASES.join(', ')}`);
  pos('pile.length', pile?.length); pos('pile.D', pile?.D); pos('pile.t', pile?.t);
  const limit = o.refusalBlowsPer025m ?? 250;
  pos('refusalBlowsPer025m', limit);
  const Lmax = Math.min(pile.length, profile.depth);
  const depths = o.depths ?? Array.from({ length: Math.floor(Lmax) }, (_, i) => i + 1);
  const srd = srdProfile({ profile, D: pile.D, t: pile.t, Lmax, dz: o.dz ?? 0.25, upperPlugged: o.upperPlugged, frictionFactorAt: o.frictionFactorAt });
  const cs = srd.cases[srdCase];
  const nSeg = pile.nSegments ?? 100, dL = pile.length / nSeg;
  const q = { shaft: 0.00254, toe: 0.00254, ...(o.quake ?? {}) };
  const dmp = { shaftClay: 0.164, shaftSand: 0.164, toeClay: 0.492, toeSand: 0.492, ...(o.damping ?? {}) };
  const out = { depth: [], srd: [], blowsPerMetre: [], blowsPer025m: [], maxCompressiveStress: [], maxTensileStress: [], enthru: [], refusal: [], refusalDepth: null, totalBlows: 0, fatigueBlocks: [] };
  let prevDepth = 0, prevBpm = null;
  for (const pen of depths) {
    pos('depths', pen);
    if (pen > Lmax + 1e-9) throw new RangeError(`depths: ${pen} m exceeds the pile length or profile depth`);
    const shaftRu = new Array(nSeg).fill(0), jS = new Array(nSeg).fill(dmp.shaftSand);
    for (let j = 0; j < nSeg; j++) {
      const a = Math.max(j * dL - (pile.length - pen), 0), b = Math.max((j + 1) * dL - (pile.length - pen), 0);
      if (b <= a) continue;
      shaftRu[j] = interpAt(srd.z, cs.shaft, b) - interpAt(srd.z, cs.shaft, a);
      jS[j] = propertiesAt(profile, Math.min(0.5 * (a + b), profile.depth)).type === 'clay' ? dmp.shaftClay : dmp.shaftSand;
    }
    const toeClay = propertiesAt(profile, pen).type === 'clay';
    const r = waveEquation({
      hammer, pile: { ...pile, nSegments: nSeg }, cushionStiffness: o.cushionStiffness, helmetMass: o.helmetMass, dtFactor: o.dtFactor,
      shaftRu, toeRu: interpAt(srd.z, cs.toe, pen), shaftQuake: q.shaft, toeQuake: q.toe, shaftDamping: jS, toeDamping: toeClay ? dmp.toeClay : dmp.toeSand,
    });
    const b025 = r.blowsPerMetre / 4, refused = b025 > limit;
    out.depth.push(pen); out.srd.push(interpAt(srd.z, cs.total, pen)); out.blowsPerMetre.push(r.blowsPerMetre); out.blowsPer025m.push(b025);
    out.maxCompressiveStress.push(r.maxCompressiveStress); out.maxTensileStress.push(r.maxTensileStress); out.enthru.push(r.enthru); out.refusal.push(refused);
    if (refused) { if (out.refusalDepth === null) out.refusalDepth = pen; } else if (out.refusalDepth === null) {
      // blows over the interval: trapezoid on blows per metre (first interval uses the value at its end)
      const blows = (prevBpm === null ? r.blowsPerMetre : 0.5 * (prevBpm + r.blowsPerMetre)) * (pen - prevDepth);
      out.totalBlows += blows;
      out.fatigueBlocks.push({ depthFrom: prevDepth, depthTo: pen, stressRange: r.maxCompressiveStress + r.maxTensileStress, blows });
    }
    prevDepth = pen; prevBpm = r.blowsPerMetre;
  }
  return out;
}

/**
 * Self-weight penetration and pile-run screening from an SRD profile.
 * Self-weight penetration: the shallowest depth at which the SRD equals the submerged weight of the pile alone,
 * and of the pile plus hammer. Pile-run risk: any depth interval below the pile-plus-hammer self-weight
 * penetration where the SRD falls below that weight again (a strong layer over a weaker one), so that the pile
 * can drop without hammer blows.
 * Reference: vertical equilibrium, SRD from srdProfile. Validity: quasi-static screening; SRD during a run is
 * rate dependent and usually lower than the continuous-driving SRD — treat results as indicative.
 * @param {{ z: number[], srd: number[], pileWeight: number, hammerWeight?: number }} o z [m] depths (ascending);
 *   srd [N] total static resistance at each depth (e.g. srdProfile(...).cases.coringLower.total); pileWeight [N]
 *   submerged weight of the pile; hammerWeight [N] weight of hammer and sleeve resting on the pile (default 0)
 * @returns {{ swpPile: number|null, swpPileHammer: number|null, runIntervals: { from: number, to: number,
 *   minSrd: number }[], runRisk: boolean }} depths [m] (null if the weight is never resisted); intervals [m];
 *   minSrd [N]
 */
export function selfWeightPenetration({ z, srd, pileWeight, hammerWeight = 0 }) {
  if (!Array.isArray(z) || !Array.isArray(srd) || z.length !== srd.length || z.length < 2) throw new RangeError('z: z and srd must be arrays of equal length >= 2');
  pos('pileWeight', pileWeight); pos('hammerWeight', hammerWeight, true);
  const first = (W) => {
    if (srd[0] >= W) return z[0];
    for (let i = 1; i < z.length; i++) if (srd[i] >= W) return z[i - 1] + ((z[i] - z[i - 1]) * (W - srd[i - 1])) / (srd[i] - srd[i - 1]);
    return null;
  };
  const W = pileWeight + hammerWeight;
  const swpPile = first(pileWeight), swpPileHammer = first(W);
  const runIntervals = [];
  if (swpPileHammer !== null) {
    let cur = null;
    for (let i = 0; i < z.length; i++) {
      if (z[i] <= swpPileHammer) continue;
      if (srd[i] < W) { if (!cur) cur = { from: z[i], to: z[i], minSrd: srd[i] }; cur.to = z[i]; cur.minSrd = Math.min(cur.minSrd, srd[i]); } else if (cur) { runIntervals.push(cur); cur = null; }
    }
    if (cur) runIntervals.push(cur);
  }
  return { swpPile, swpPileHammer, runIntervals, runRisk: runIntervals.length > 0 };
}

/**
 * Submerged weight of a steel pipe pile: W = (ρs − ρw) g A L.
 * Reference: Archimedes. Validity: fully submerged, free-flooding pile.
 * @param {{ D: number, t: number, length: number, rho?: number, rhoWater?: number }} o D [m]; t [m]; length [m];
 *   rho [kg/m³] (7850); rhoWater [kg/m³] (1025)
 * @returns {number} submerged weight [N]
 */
export function pileSubmergedWeight({ D, t, length, rho = 7850, rhoWater = 1025 }) {
  pos('D', D); pos('t', t); pos('length', length); pos('rho', rho); pos('rhoWater', rhoWater, true);
  return ((rho - rhoWater) * G * Math.PI * (D * D - (D - 2 * t) ** 2) * length) / 4;
}
