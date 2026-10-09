/**
 * Reduced-order aero-hydro-servo-elastic time-domain simulation of bottom-fixed and floating offshore wind turbines.
 *
 * Model
 *  - rotor/drivetrain: rigid rotor speed DOF (optional torsional shaft DOF), generator torque and collective pitch
 *    from control.js, second-order pitch actuator;
 *  - aerodynamics: thrust and torque from pre-computed Ct/Cq(tip-speed ratio, pitch) surfaces evaluated with the
 *    rotor-effective wind speed RELATIVE to the moving hub, so aerodynamic damping emerges from the model; optional
 *    reduced-order dynamic inflow;
 *  - support structure: N fore-aft modal coordinates (supplied, or from the beam model of buildSupportStructure);
 *  - floater: rigid-body surge, heave and pitch with constant mass + added mass, stiffness, linear and quadratic
 *    damping, and first-order wave force RAOs;
 *  - hydrodynamics: linear irregular waves (JONSWAP, seeded), Morison loads on the submerged stations with relative
 *    velocity in the drag term, optional MacCamy-Fuchs diffraction correction of the inertia term;
 *  - integration: classical RK4 with fixed step, controller and actuator sampled once per step.
 *
 * Modelling limits (state them with any result): rigid blades; modal (linear) support structure in one plane;
 * rotor-effective wind (no blade-resolved turbulence, no 1P/3P excitation); quasi-steady rotor aerodynamics from
 * steady BEM surfaces; linear (Airy) waves without stretching, no slamming or breaking-wave loads; frequency-
 * independent added mass and damping for floaters, no second-order wave forces; side-side response not computed.
 */

import { eigGenSym } from '../math/eigen.js';
import { ifft, psdWelch } from '../math/fft.js';
import { makeRng } from '../math/rng.js';
import { besselJ0, besselJ1, besselY0, besselY1 } from '../math/special.js';
import { surfaceLookup, inductionFromCt } from './bem.js';
import { createController, pitchActuatorStep } from './control.js';

const G = 9.80665;
const D2R = Math.PI / 180;
const TWO_PI = 2 * Math.PI;

function need(cond, name, msg) {
  if (!cond) throw new RangeError(`${name}: ${msg}`);
}

/**
 * Wave number from the linear dispersion relation omega^2 = g k tanh(k d) (Newton iteration).
 * @param {number} omega angular frequency [rad/s] (> 0) @param {number} depth water depth [m] (> 0)
 * @returns {number} k [rad/m]
 */
export function waveNumber(omega, depth) {
  need(omega > 0 && depth > 0, 'omega/depth', 'must be > 0');
  const k0 = omega * omega / G;
  let k = k0 / Math.sqrt(Math.tanh(k0 * depth));
  for (let i = 0; i < 50; i++) {
    const th = Math.tanh(k * depth);
    const f = k * th - k0;
    const dk = f / (th + k * depth * (1 - th * th));
    k -= dk;
    if (Math.abs(dk) < 1e-13 * k) break;
  }
  return k;
}

/**
 * JONSWAP spectral shape S(f) ~ f^-5 exp(-1.25 (f/fp)^-4) gamma^exp(-(f - fp)^2 / (2 s^2 fp^2)), s = 0.07 (f <= fp)
 * or 0.09 (f > fp) (Hasselmann et al. 1973; DNV-RP-C205 3.5.5), returned UN-normalised: scale to m0 = Hs^2/16.
 * @param {number} f frequency [Hz] (> 0) @param {number} Tp peak period [s] @param {number} [gamma=3.3] peak enhancement [-]
 * @returns {number} relative spectral density [-]
 */
export function jonswapShape(f, Tp, gamma = 3.3) {
  const fp = 1 / Tp;
  const s = f <= fp ? 0.07 : 0.09;
  const r = Math.exp(-((f - fp) ** 2) / (2 * s * s * fp * fp));
  return Math.pow(f, -5) * Math.exp(-1.25 * Math.pow(f / fp, -4)) * Math.pow(gamma, r);
}

/**
 * Complex effective inertia coefficient of a surface-piercing vertical cylinder from MacCamy-Fuchs linear
 * diffraction theory, for the exp(+i omega t) time convention:
 *   CM_eff(ka) = 4 (Y1'(ka) - i J1'(ka)) / (pi (ka)^2 (J1'(ka)^2 + Y1'(ka)^2)),
 * i.e. magnitude 4 / (pi (ka)^2 sqrt(J1'^2 + Y1'^2)) with the force lagging the undisturbed acceleration at the
 * cylinder axis; CM_eff -> 2 as ka -> 0. Reference: MacCamy, R.C. and Fuchs, R.A. (1954), "Wave forces on piles: a
 * diffraction theory", Beach Erosion Board TM-69; derivation via the Wronskian of J1 and H1.
 * @param {number} ka wave number times cylinder radius [-] (> 0)
 * @returns {{re:number, im:number, magnitude:number, phaseLag:number}} complex coefficient [-], |CM| [-], lag [rad]
 */
export function macCamyFuchsCm(ka) {
  need(ka > 0, 'ka', 'must be > 0');
  if (ka < 1e-3) return { re: 2, im: 0, magnitude: 2, phaseLag: 0 };
  const j1p = besselJ0(ka) - besselJ1(ka) / ka;
  const y1p = besselY0(ka) - besselY1(ka) / ka;
  const den = Math.PI * ka * ka * (j1p * j1p + y1p * y1p);
  const re = 4 * y1p / den; const im = -4 * j1p / den;
  return { re, im, magnitude: Math.hypot(re, im), phaseLag: Math.atan2(-im, re) };
}

/**
 * Linear irregular long-crested wave: surface elevation and horizontal particle velocity/acceleration at given
 * elevations on the vertical axis x = 0, synthesised by inverse FFT from a JONSWAP spectrum normalised so that the
 * represented band carries m0 = Hs^2/16. Kinematics from Airy theory, evaluated at min(z, 0) (no stretching).
 * @param {object} o
 * @param {number} o.Hs significant wave height [m] @param {number} o.Tp peak period [s] @param {number} [o.gamma=3.3]
 * @param {number} o.depth water depth [m]
 * @param {number} o.dt sample interval [s] @param {number} o.n number of samples (made even)
 * @param {number} o.seed integer seed
 * @param {number[]} [o.z=[]] elevations, 0 at still water level, negative below [m]
 * @param {number[]} [o.diameters] cylinder diameter at each elevation [m]; when given, `accMF` is returned
 * @param {number} [o.fMax] cut-off frequency [Hz]; default max(0.5, 4/Tp), capped at Nyquist
 * @returns {{t:Float64Array, eta:Float64Array, u:Float64Array[], acc:Float64Array[], accMF:(Float64Array[]|null), m0:number, dt:number}}
 *   time [s], elevation [m], horizontal velocity [m/s] and acceleration [m/s2] per elevation; `accMF` is the
 *   acceleration filtered by the MacCamy-Fuchs complex coefficient (so that inertia force per length =
 *   rho pi D^2/4 * accMF, replacing CM * acc); zeroth spectral moment [m2]
 */
export function linearIrregularWave(o) {
  need(o && o.Hs >= 0 && o.Tp > 0, 'Hs/Tp', 'Hs >= 0 and Tp > 0 required');
  need(o.depth > 0, 'depth', 'must be > 0');
  need(o.dt > 0 && o.n >= 4, 'dt/n', 'dt > 0 and n >= 4 required');
  const n = o.n % 2 ? o.n + 1 : o.n;
  const T = n * o.dt; const df = 1 / T;
  const fMax = Math.min(0.5 / o.dt, o.fMax ?? Math.max(0.5, 4 / o.Tp));
  const mMax = Math.min(n / 2 - 1, Math.floor(fMax / df));
  const zs = o.z || [];
  const rng = makeRng(o.seed ?? 1);
  const S = new Float64Array(mMax + 1);
  let sum = 0;
  for (let m = 1; m <= mMax; m++) { S[m] = jonswapShape(m * df, o.Tp, o.gamma ?? 3.3); sum += S[m] * df; }
  const scale = sum > 0 ? o.Hs * o.Hs / 16 / sum : 0;
  const eRe = new Float64Array(n); const eIm = new Float64Array(n);
  const uRe = zs.map(() => new Float64Array(n)); const uIm = zs.map(() => new Float64Array(n));
  const aRe = zs.map(() => new Float64Array(n)); const aIm = zs.map(() => new Float64Array(n));
  const mf = !!o.diameters;
  const mRe = mf ? zs.map(() => new Float64Array(n)) : null; const mIm = mf ? zs.map(() => new Float64Array(n)) : null;
  for (let m = 1; m <= mMax; m++) {
    const amp = Math.sqrt(2 * S[m] * scale * df);
    const ph = TWO_PI * rng.next();
    const cr = 0.5 * amp * n * Math.cos(ph); const ci = 0.5 * amp * n * Math.sin(ph);
    eRe[m] = cr; eIm[m] = ci; eRe[n - m] = cr; eIm[n - m] = -ci;
    const w = TWO_PI * m * df;
    const k = waveNumber(w, o.depth);
    const kd = k * o.depth;
    for (let s = 0; s < zs.length; s++) {
      const z = Math.min(0, zs[s]);
      if (z < -o.depth - 1e-9) continue;
      // cosh(k (z + d)) / sinh(k d), overflow-safe
      const Gz = kd > 30 ? Math.exp(k * z) : Math.cosh(k * (z + o.depth)) / Math.sinh(kd);
      const ur = w * Gz * cr; const ui = w * Gz * ci;
      uRe[s][m] = ur; uIm[s][m] = ui; uRe[s][n - m] = ur; uIm[s][n - m] = -ui;
      const ar = -w * ui; const ai = w * ur;            // i omega * u
      aRe[s][m] = ar; aIm[s][m] = ai; aRe[s][n - m] = ar; aIm[s][n - m] = -ai;
      if (mf) {
        const c = macCamyFuchsCm(0.5 * k * o.diameters[s]);
        const xr = c.re * ar - c.im * ai; const xi = c.re * ai + c.im * ar;
        mRe[s][m] = xr; mIm[s][m] = xi; mRe[s][n - m] = xr; mIm[s][n - m] = -xi;
      }
    }
  }
  const real = (re, im) => Float64Array.from(ifft(re, im).re);
  const t = new Float64Array(n);
  for (let i = 0; i < n; i++) t[i] = i * o.dt;
  return {
    t, eta: real(eRe, eIm), u: zs.map((_, s) => real(uRe[s], uIm[s])), acc: zs.map((_, s) => real(aRe[s], aIm[s])),
    accMF: mf ? zs.map((_, s) => real(mRe[s], mIm[s])) : null, m0: sum * scale, dt: o.dt,
  };
}

/**
 * Beam (stick) model of a monopile + tower support structure and its fore-aft modes.
 * Two-node Euler-Bernoulli elements (cubic Hermite, consistent mass) with element-mean section properties from
 * linearly tapered tubular segments; lumped rotor-nacelle mass (and optional rotary inertia) at the tower top;
 * foundation represented at the lowest node (mudline) by a coupled lateral/rocking spring matrix [[Kh, Khr],
 * [Khr, Kr]] or a clamp; submerged elements carry hydrodynamic added mass rho_w Ca pi D^2/4 and, optionally,
 * the enclosed water. Shear deformation, axial load (P-delta) and soil damping are not modelled.
 * @param {object} o
 * @param {Array<{zBase:number,zTop:number,dBase:number,dTop:number,tBase:number,tTop:number,density?:number,E?:number}>} o.segments
 *   contiguous tubular segments ordered upward, z = 0 at still water level [m]; the lowest zBase is the mudline
 * @param {number} o.rnaMass rotor-nacelle-assembly mass at the tower top [kg]
 * @param {number} [o.rnaInertia=0] RNA rotary inertia about the lateral axis [kg m2]
 * @param {number} o.zHub hub elevation [m] (hub displacement = top displacement + top rotation x (zHub - zTop))
 * @param {number} [o.zTowerBase] elevation of the tower-base output section [m]; default the first segment joint
 * @param {object|'fixed'} [o.foundation='fixed'] { Kh [N/m], Kr [N m/rad], Khr [N] (default 0) }
 * @param {number} [o.nElem=40] number of beam elements (approximately; segment joints are kept as nodes)
 * @param {number} [o.nModes=2] modes returned
 * @param {number} [o.zeta=0.01] structural damping ratio assigned to every mode [-]
 * @param {number} [o.rhoWater=1025] [kg/m3] @param {number} [o.addedMassCa=1] [-] @param {boolean} [o.internalWater=true]
 * @param {number} [o.cd=1.0] Morison drag coefficient of the wet stations [-] (0.65 smooth - 1.05 rough, DNV-RP-C205)
 * @param {number} [o.cm=2.0] Morison inertia coefficient of the wet stations [-]
 * @param {number} [o.cdAir=0.7] wind drag coefficient of the dry stations [-] (0.5 - 1.2 depending on Reynolds number)
 * @returns {object} structure for simulateCoupled: { zHub, zMudline, zTowerBase, zTop, modes: [{ frequency [Hz],
 *   mass [kg], stiffness [N/m], zeta, phiHub, slopeHub [1/m], phiStations[], momentBase, shearBase, momentMudline,
 *   shearMudline (section forces per unit modal coordinate, N m and N), shape: Float64Array (lateral displacement at
 *   nodes) }], stations: [{ z, D, dz, medium, cd, cm }], nodes: Float64Array, mass (structure without RNA) [kg] };
 *   modes are scaled to unit hub displacement
 */
export function buildSupportStructure(o) {
  need(o && Array.isArray(o.segments) && o.segments.length >= 1, 'segments', 'at least one segment required');
  need(o.rnaMass >= 0, 'rnaMass', 'must be >= 0');
  const segs = o.segments;
  for (let i = 1; i < segs.length; i++) need(Math.abs(segs[i].zBase - segs[i - 1].zTop) < 1e-6, `segments[${i}]`, 'must start where the previous segment ends');
  const zMud = segs[0].zBase;
  const zTop = segs[segs.length - 1].zTop;
  need(o.zHub >= zTop - 1e-9, 'zHub', 'must not be below the tower top');
  const total = zTop - zMud;
  const nTarget = o.nElem ?? 40;
  const rhoW = o.rhoWater ?? 1025;
  const Ca = o.addedMassCa ?? 1;
  const zNodes = [zMud];
  const elems = [];
  for (const sg of segs) {
    const L = sg.zTop - sg.zBase;
    need(L > 0, 'segments', 'zTop must exceed zBase');
    // keep the still water level as a node so that wet/dry elements are distinct
    const cuts = [sg.zBase, sg.zTop];
    if (sg.zBase < -1e-9 && sg.zTop > 1e-9) cuts.splice(1, 0, 0);
    for (let c = 0; c < cuts.length - 1; c++) {
      const Lc = cuts[c + 1] - cuts[c];
      const ne = Math.max(1, Math.round(nTarget * Lc / total));
      for (let e = 0; e < ne; e++) {
        const za = cuts[c] + Lc * e / ne; const zb = cuts[c] + Lc * (e + 1) / ne;
        const s = (0.5 * (za + zb) - sg.zBase) / L;
        const D = sg.dBase + (sg.dTop - sg.dBase) * s;
        const t = sg.tBase + (sg.tTop - sg.tBase) * s;
        const A = Math.PI * (D - t) * t;
        const I = Math.PI / 64 * (D ** 4 - (D - 2 * t) ** 4);
        let m = (sg.density ?? 7850) * A;
        const wet = 0.5 * (za + zb) < 0;
        if (wet) m += rhoW * Ca * Math.PI * D * D / 4 + ((o.internalWater ?? true) ? rhoW * Math.PI * (D - 2 * t) ** 2 / 4 : 0);
        elems.push({ za, zb, D, EI: (sg.E ?? 210e9) * I, m, mSteel: (sg.density ?? 7850) * A, wet });
        zNodes.push(zb);
      }
    }
  }
  const nn = zNodes.length;
  const nd = 2 * nn;
  const K = []; const M = [];
  for (let i = 0; i < nd; i++) { K.push(new Array(nd).fill(0)); M.push(new Array(nd).fill(0)); }
  let massStruct = 0;
  elems.forEach((el, e) => {
    const L = el.zb - el.za;
    const k = el.EI / L ** 3;
    const ke = [[12 * k, 6 * L * k, -12 * k, 6 * L * k], [6 * L * k, 4 * L * L * k, -6 * L * k, 2 * L * L * k],
      [-12 * k, -6 * L * k, 12 * k, -6 * L * k], [6 * L * k, 2 * L * L * k, -6 * L * k, 4 * L * L * k]];
    const c = el.m * L / 420;
    const me = [[156 * c, 22 * L * c, 54 * c, -13 * L * c], [22 * L * c, 4 * L * L * c, 13 * L * c, -3 * L * L * c],
      [54 * c, 13 * L * c, 156 * c, -22 * L * c], [-13 * L * c, -3 * L * L * c, -22 * L * c, 4 * L * L * c]];
    for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) { K[2 * e + a][2 * e + b] += ke[a][b]; M[2 * e + a][2 * e + b] += me[a][b]; }
    massStruct += el.mSteel * L;
  });
  M[nd - 2][nd - 2] += o.rnaMass;
  M[nd - 1][nd - 1] += o.rnaInertia ?? 0;
  const fixed = o.foundation === undefined || o.foundation === 'fixed';
  // K0 keeps the unrestrained stiffness for the free-body section-force recovery
  const K0 = K.map((r) => r.slice());
  let off = 0;
  if (!fixed) {
    need(o.foundation.Kh > 0 && o.foundation.Kr > 0, 'foundation', 'Kh and Kr must be > 0');
    K[0][0] += o.foundation.Kh; K[1][1] += o.foundation.Kr;
    K[0][1] += o.foundation.Khr ?? 0; K[1][0] += o.foundation.Khr ?? 0;
  } else off = 2;
  const Kr = K.slice(off).map((r) => r.slice(off));
  const Mr = M.slice(off).map((r) => r.slice(off));
  const nModes = Math.min(o.nModes ?? 2, nd - off);
  const eig = eigGenSym(Kr, Mr, { nModes });
  const zTB = o.zTowerBase ?? (segs.length > 1 ? segs[0].zTop : zMud);
  const iTB = zNodes.findIndex((z) => Math.abs(z - zTB) < 1e-6);
  need(iTB >= 0, 'zTowerBase', 'must coincide with a segment joint');
  const stations = elems.map((el) => ({ z: 0.5 * (el.za + el.zb), D: el.D, dz: el.zb - el.za, medium: el.wet ? 'water' : 'air',
    cd: el.wet ? (o.cd ?? 1.0) : (o.cdAir ?? 0.7), cm: o.cm ?? 2.0 }));
  const modes = [];
  for (let q = 0; q < nModes; q++) {
    const lam = eig.values[q];
    const v = new Float64Array(nd);
    for (let i = off; i < nd; i++) v[i] = eig.vectors[q][i - off];
    let hub = v[nd - 2] + v[nd - 1] * (o.zHub - zTop);
    if (Math.abs(hub) < 1e-9) { hub = 0; for (let i = 0; i < nd; i += 2) if (Math.abs(v[i]) > Math.abs(hub)) hub = v[i]; }
    for (let i = 0; i < nd; i++) v[i] /= hub;
    let mg = 0;
    for (let i = 0; i < nd; i++) { let s = 0; for (let j = 0; j < nd; j++) s += M[i][j] * v[j]; mg += v[i] * s; }
    const f = new Float64Array(nd);
    for (let i = 0; i < nd; i++) { let s = 0; for (let j = 0; j < nd; j++) s += K0[i][j] * v[j]; f[i] = s; }
    const section = (iNode) => {
      let V = 0; let Mo = 0;
      for (let j = iNode + 1; j < nn; j++) { V += f[2 * j]; Mo += f[2 * j] * (zNodes[j] - zNodes[iNode]) + f[2 * j + 1]; }
      return [Mo, V];
    };
    const sm = section(0); const sb = section(iTB);
    const shape = new Float64Array(nn);
    for (let i = 0; i < nn; i++) shape[i] = v[2 * i];
    modes.push({
      frequency: Math.sqrt(Math.max(0, lam)) / TWO_PI, mass: mg, stiffness: lam * mg, zeta: o.zeta ?? 0.01,
      phiHub: v[nd - 2] + v[nd - 1] * (o.zHub - zTop), slopeHub: v[nd - 1],
      phiStations: elems.map((el, e) => 0.5 * (v[2 * e] + v[2 * e + 2]) + (el.zb - el.za) / 8 * (v[2 * e + 1] - v[2 * e + 3])),
      momentMudline: sm[0], shearMudline: sm[1], momentBase: sb[0], shearBase: sb[1], shape,
    });
  }
  return { zHub: o.zHub, zMudline: zMud, zTowerBase: zTB, zTop, modes, stations, nodes: Float64Array.from(zNodes), mass: massStruct };
}

/**
 * Support-structure segments for a rotor's tower on a monopile (convenience for buildSupportStructure): monopile
 * from the mudline to the transition elevation, tower from there to the tower top.
 * @param {object} o
 * @param {object} o.tower rotor.tower record (height, baseDiameter, topDiameter, baseThickness, topThickness, density, E)
 * @param {number} o.waterDepth [m] @param {number} o.transitionZ elevation of the tower base above still water [m]
 * @param {number} o.monopileD [m] @param {number} o.monopileT wall thickness [m]
 * @param {number} [o.monopileDensity=7850] [kg/m3]
 * @returns {Array<object>} segments
 */
export function monopileTowerSegments(o) {
  need(o && o.tower && o.waterDepth > 0, 'tower/waterDepth', 'required');
  const t = o.tower;
  return [
    { zBase: -o.waterDepth, zTop: o.transitionZ, dBase: o.monopileD, dTop: o.monopileD, tBase: o.monopileT, tTop: o.monopileT, density: o.monopileDensity ?? 7850, E: t.E ?? 210e9 },
    { zBase: o.transitionZ, zTop: o.transitionZ + t.height, dBase: t.baseDiameter, dTop: t.topDiameter, tBase: t.baseThickness, tTop: t.topThickness, density: t.density ?? 7850, E: t.E ?? 210e9 },
  ];
}

/**
 * Mean, standard deviation, minimum and maximum of a series (optionally from index i0).
 * @param {Float64Array|number[]} x series @param {number} [i0=0] first index used
 * @returns {{mean:number, std:number, min:number, max:number, absMax:number}} in the units of x
 */
export function seriesStats(x, i0 = 0) {
  let s = 0; let mn = Infinity; let mx = -Infinity;
  const n = x.length - i0;
  need(n >= 1, 'x', 'empty series');
  for (let i = i0; i < x.length; i++) { s += x[i]; if (x[i] < mn) mn = x[i]; if (x[i] > mx) mx = x[i]; }
  const mean = s / n;
  let v = 0;
  for (let i = i0; i < x.length; i++) v += (x[i] - mean) ** 2;
  return { mean, std: Math.sqrt(v / n), min: mn, max: mx, absMax: Math.max(Math.abs(mn), Math.abs(mx)) };
}

/**
 * Rainflow cycle counting (four-point algorithm, ASTM E1049-85 equivalent); the residue is counted as half cycles.
 * @param {Float64Array|number[]} x load series
 * @returns {{ranges:number[], means:number[], counts:number[]}} cycle ranges and means (units of x), counts (1 or 0.5)
 */
export function rainflowCount(x) {
  const tp = [];
  for (let i = 0; i < x.length; i++) {
    const v = x[i];
    const m = tp.length;
    if (m >= 1 && v === tp[m - 1]) continue;
    if (m >= 2 && (tp[m - 1] - tp[m - 2]) * (v - tp[m - 1]) > 0) tp[m - 1] = v; else tp.push(v);
  }
  const ranges = []; const means = []; const counts = [];
  const st = [];
  for (const p of tp) {
    st.push(p);
    while (st.length >= 4) {
      const m = st.length;
      const a = st[m - 4]; const b = st[m - 3]; const c = st[m - 2]; const d = st[m - 1];
      const inner = Math.abs(c - b);
      if (inner <= Math.abs(b - a) && inner <= Math.abs(d - c)) {
        ranges.push(inner); means.push(0.5 * (b + c)); counts.push(1);
        st.splice(m - 3, 2);
      } else break;
    }
  }
  for (let i = 0; i < st.length - 1; i++) { ranges.push(Math.abs(st[i + 1] - st[i])); means.push(0.5 * (st[i] + st[i + 1])); counts.push(0.5); }
  return { ranges, means, counts };
}

/**
 * Damage-equivalent load: the constant range that, applied nEq times, gives the same Miner damage on an S-N curve
 * of slope m:  DEL = ( sum n_i R_i^m / nEq )^(1/m).
 * @param {Float64Array|number[]} x load series
 * @param {object} o
 * @param {number} o.m Woehler exponent [-] (steel welds 3 - 5, composites 9 - 12)
 * @param {number} o.nEq number of equivalent cycles [-] (e.g. duration x 1 Hz)
 * @param {function(Float64Array|number[]):{ranges:number[],counts:number[]}} [o.rainflow] cycle counter; default rainflowCount
 * @returns {{del:number, sumNRm:number, nCycles:number}} DEL (units of x), sum n R^m, number of counted cycles
 */
export function damageEquivalentLoad(x, o) {
  need(o && o.m > 0 && o.nEq > 0, 'm/nEq', 'must be > 0');
  const rf = (o.rainflow || rainflowCount)(x);
  let s = 0; let nc = 0;
  for (let i = 0; i < rf.ranges.length; i++) { s += rf.counts[i] * Math.pow(rf.ranges[i], o.m); nc += rf.counts[i]; }
  return { del: Math.pow(s / o.nEq, 1 / o.m), sumNRm: s, nCycles: nc };
}

/**
 * Damping ratio and damped frequency of a decaying oscillation from the logarithmic decrement of successive
 * positive peaks about the mean of the last quarter of the record: delta = ln(x_k / x_{k+n}) / n,
 * zeta = delta / sqrt(4 pi^2 + delta^2).
 * @param {Float64Array|number[]} x response series @param {number} dt [s]
 * @param {object} [o] @param {number} [o.maxPeaks=8] peaks used
 * @returns {{zeta:number, frequency:number, nPeaks:number}} damping ratio [-], damped frequency [Hz]
 */
export function logDecrement(x, dt, o = {}) {
  const n = x.length;
  let mean = 0; let cnt = 0;
  for (let i = Math.floor(0.75 * n); i < n; i++) { mean += x[i]; cnt++; }
  mean /= cnt;
  const pk = []; const tk = [];
  for (let i = 1; i < n - 1 && pk.length < (o.maxPeaks ?? 8); i++) {
    const a = x[i - 1] - mean; const b = x[i] - mean; const c = x[i + 1] - mean;
    if (b > a && b >= c && b > 0) {
      // parabolic refinement of the peak
      const den = a - 2 * b + c;
      const sft = den !== 0 ? 0.5 * (a - c) / den : 0;
      pk.push(b - 0.25 * (a - c) * sft); tk.push((i + sft) * dt);
    }
  }
  need(pk.length >= 2, 'x', 'fewer than two peaks found');
  const m = pk.length - 1;
  const delta = Math.log(pk[0] / pk[m]) / m;
  return { zeta: delta / Math.sqrt(4 * Math.PI * Math.PI + delta * delta), frequency: m / (tk[m] - tk[0]), nPeaks: pk.length };
}

/**
 * Identify spectral peaks near named frequencies: for each candidate the largest PSD value within +-`band` of the
 * frequency and its ratio to the median PSD level.
 * @param {Float64Array|number[]} f frequencies [Hz] @param {Float64Array|number[]} psd spectrum
 * @param {Array<{name:string, f:number}>} candidates e.g. 1P, 3P, natural frequencies [Hz]
 * @param {number} [band=0.1] relative half-width of the search band [-]
 * @returns {Array<{name:string, f:number, fPeak:number, psd:number, prominence:number}>} peak frequency [Hz], PSD value,
 *   ratio to the median level [-]
 */
export function identifyPeaks(f, psd, candidates, band = 0.1) {
  const sorted = Float64Array.from(psd).sort();
  const med = sorted[Math.floor(sorted.length / 2)] || 1e-300;
  return candidates.map((c) => {
    let best = -Infinity; let fb = NaN;
    for (let i = 0; i < f.length; i++) if (f[i] >= c.f * (1 - band) && f[i] <= c.f * (1 + band) && psd[i] > best) { best = psd[i]; fb = f[i]; }
    return { name: c.name, f: c.f, fPeak: fb, psd: best, prominence: best / med };
  });
}

// Allocation-free bilinear Ct/Cq lookup (same interpolation as bem.js surfaceLookup); result in LK.
const LK = { ct: 0, cq: 0 };
function cellOf(xs, x) {
  const n = xs.length;
  if (x <= xs[0]) return 0;
  if (x >= xs[n - 1]) return n - 2;
  let lo = 0; let hi = n - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
  return lo;
}
function lookupCtCq(sf, tsr, pitchDeg) {
  const T = sf.tsr; const P = sf.pitchDeg;
  const i = cellOf(T, tsr); const j = cellOf(P, pitchDeg);
  let s = (tsr - T[i]) / (T[i + 1] - T[i]); let t = (pitchDeg - P[j]) / (P[j + 1] - P[j]);
  s = s < 0 ? 0 : (s > 1 ? 1 : s); t = t < 0 ? 0 : (t > 1 ? 1 : t);
  const w00 = (1 - s) * (1 - t); const w10 = s * (1 - t); const w01 = (1 - s) * t; const w11 = s * t;
  const c0 = sf.ct[j]; const c1 = sf.ct[j + 1]; const q0 = sf.cq[j]; const q1 = sf.cq[j + 1];
  LK.ct = w00 * c0[i] + w10 * c0[i + 1] + w01 * c1[i] + w11 * c1[i + 1];
  LK.cq = w00 * q0[i] + w10 * q0[i + 1] + w01 * q1[i] + w11 * q1[i + 1];
}

function sampleInput(v, k, t, dt) {
  if (typeof v === 'number') return v;
  if (typeof v === 'function') return v(t);
  const x = t / dt;
  const i = Math.min(v.length - 1, Math.max(0, Math.floor(x)));
  const j = Math.min(v.length - 1, i + 1);
  return v[i] + (x - i) * (v[j] - v[i]);
}

function inv3(m) {
  const [a, b, c] = m[0]; const [d, e, f] = m[1]; const [g, h, i] = m[2];
  const det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
  need(Math.abs(det) > 0, 'floater.mass', 'singular mass matrix');
  return [[(e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det],
    [(f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det],
    [(d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det]];
}

function as3(m, name) {
  if (m === undefined) return [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  if (typeof m[0] === 'number') return [[m[0], 0, 0], [0, m[1], 0], [0, 0, m[2]]];
  need(m.length === 3 && m[0].length === 3, name, '3 diagonal terms or a 3x3 matrix required');
  return m;
}

/**
 * Coupled time-domain simulation (see the module header for the model and its limits).
 * @param {object} o
 * @param {object} o.surface cpSurface result (bem.js) for the rotor
 * @param {object} o.turbine { inertia [kg m2] (drivetrain on the low-speed shaft), hubHeight [m] (above still water
 *   level), shaft?: { stiffness [N m/rad], damping [N m s/rad], generatorInertiaLss [kg m2] } (torsional DOF; the
 *   rotor inertia is then inertia - generatorInertiaLss) }
 * @param {object|null} o.controller controller configuration (control.js); null holds rotor speed and pitch at their
 *   initial values (prescribed operation, parked or idling)
 * @param {object} [o.structure] { zHub, zMudline, zTowerBase, modes: [{ mass, stiffness, zeta, phiHub, slopeHub,
 *   phiStations[], momentBase, shearBase, momentMudline, shearMudline }], stations: [{ z, D, dz, medium, cd, cm }] }
 *   (buildSupportStructure result or caller-supplied). Modal mass must include hydrodynamic added mass. Section
 *   forces per unit modal coordinate may be omitted (the equivalent-tip-force estimate stiffness/phiHub x lever is used).
 * @param {object} [o.support] arguments for buildSupportStructure, used when `structure` is not supplied
 * @param {object} [o.floater] { mass (3 diagonal terms or 3x3 for surge, heave, pitch; kg, kg, kg m2),
 *   addedMass, stiffness (hydrostatic + mooring; N/m, N/m, N m/rad), damping (linear), quadDamping (N s2/m2, ...,
 *   N m s2/rad2; diagonal), refZ (elevation of the motion reference point) [m], towerBaseZ [m], rnaMass [kg],
 *   forceRao?: { f: [Hz], surge: { amp: [N/m], phaseDeg: [] }, heave: {..}, pitch: { amp: [N m/m], .. } } }
 * @param {number|number[]|Float64Array|function(number):number} o.wind rotor-effective wind speed at hub [m/s]
 *   (series sampled at o.dt)
 * @param {object} [o.waves] { Hs [m], Tp [s], gamma, seed, depth [m], macCamyFuchs (bool), current [m/s] (uniform) }
 * @param {number} o.dt time step [s] @param {number} o.duration [s]
 * @param {number} [o.tTransient=0] start of the statistics window [s]
 * @param {number} [o.rho=1.225] air density @param {number} [o.rhoWater=1025]
 * @param {number} [o.yawErrorDeg=0] steady yaw misalignment (wind component cos(yaw) used) [deg]
 * @param {number} [o.shearExp=0] wind shear exponent for wind drag on dry stations [-]
 * @param {boolean} [o.towerWindDrag=false] apply wind drag to the dry stations
 * @param {{ct:number, cq?:number}} [o.aeroOverride] constant thrust/torque coefficients replacing the surface (parked rotor)
 * @param {boolean} [o.dynamicInflow=false] reduced-order dynamic inflow: the mean induced velocity v = a V (a from
 *   Ct by momentum theory with the Buhl closure) is lagged with the Oye time constant tau1 and the surfaces are
 *   evaluated at the equivalent wind speed V + (v_qs - v_lag)/(1 - a)
 * @param {object} [o.initial] { omega [rad/s], pitchDeg, q[], qd[], platform[3], platformVel[3], mode }
 * @param {function(number):object} [o.events] t [s] -> { command?: 'start'|'stop'|'estop', gridLoss?: boolean }
 * @param {object} [o.actuator] pitch actuator constants (omegaAct [rad/s], zetaAct)
 * @param {object} [o.del] { m: Woehler exponent (default 4), nEq: equivalent cycles (default 1 Hz x window), rainflow: hook }
 * @param {number} [o.psdNfft=1024] Welch segment length for the load spectrum
 * @returns {object} { t, wind, waveElevation, rotorSpeed [rad/s], genSpeed [rad/s], pitchDeg, genTorque [N m],
 *   power [W] (electrical), thrust [N], aeroTorque [N m], towerTopDisp [m], towerTopVel [m/s], towerTopAcc [m/s2],
 *   towerBaseMoment [N m], towerBaseShear [N], mudlineMoment [N m], mudlineShear [N], modal: Float64Array[]
 *   (coordinates), platform: { surge [m], heave [m], pitch [rad] } | null, mode: string[], stats: { channel:
 *   { mean, std, min, max, absMax } }, del: { towerBaseMoment, mudlineMoment } (N m), psd: { f [Hz],
 *   towerBaseMoment [(N m)2/Hz], peaks: [...] }, final: state for restart (pass as `initial`), structure }
 */
export function simulateCoupled(o) {
  need(o && o.dt > 0 && o.duration > 0, 'dt/duration', 'must be > 0');
  const sf = o.surface;
  need(sf && sf.R > 0, 'surface', 'cpSurface result required');
  const tb = o.turbine;
  need(tb && tb.inertia > 0, 'turbine.inertia', 'must be > 0');
  need(tb.hubHeight > 0, 'turbine.hubHeight', 'must be > 0');
  const rho = o.rho ?? 1.225;
  const rhoW = o.rhoWater ?? 1025;
  const dt = o.dt;
  const n = Math.round(o.duration / dt);
  const cfg = o.controller || null;
  const N = cfg ? (cfg.gearboxRatio ?? 1) : 1;
  const eta = cfg ? (cfg.efficiency ?? 1) : 1;
  const R = sf.R; const area = Math.PI * R * R;
  const cosYaw = Math.cos((o.yawErrorDeg ?? 0) * D2R);
  const fl = o.floater || null;
  let str = o.structure || null;
  if (!str && o.support) str = buildSupportStructure(o.support);
  const modes = str ? str.modes : [];
  const nm = modes.length;
  const stations = str ? (str.stations || []) : [];
  const ns = stations.length;
  const zHub = str ? (str.zHub ?? tb.hubHeight) : tb.hubHeight;
  const zMud = str ? (str.zMudline ?? 0) : 0;
  const zTB = str ? (str.zTowerBase ?? zMud) : (fl ? (fl.towerBaseZ ?? 0) : 0);
  const mK = new Float64Array(nm); const mM = new Float64Array(nm); const mC = new Float64Array(nm); const mPhi = new Float64Array(nm);
  let wMax = 0;
  for (let i = 0; i < nm; i++) {
    const md = modes[i];
    need(md.mass > 0 && md.stiffness > 0, `modes[${i}]`, 'mass and stiffness must be > 0');
    mK[i] = md.stiffness; mM[i] = md.mass; mPhi[i] = md.phiHub ?? 1;
    mC[i] = 2 * (md.zeta ?? 0) * Math.sqrt(md.stiffness * md.mass);
    if (ns) need(md.phiStations && md.phiStations.length === ns, `modes[${i}].phiStations`, 'must match stations');
    wMax = Math.max(wMax, Math.sqrt(md.stiffness / md.mass));
  }
  const phiS = new Float64Array(nm * ns);
  for (let i = 0; i < nm; i++) for (let s = 0; s < ns; s++) phiS[i * ns + s] = modes[i].phiStations[s];
  const stZ = Float64Array.from(stations, (s) => s.z);
  need(wMax * dt < 2.0, 'dt', `too large for the highest structural mode (omega dt = ${(wMax * dt).toFixed(2)} >= 2); reduce dt or the number of modes`);
  // section-force participation per unit modal coordinate
  const sMB = new Float64Array(nm); const sVB = new Float64Array(nm); const sMM = new Float64Array(nm); const sVM = new Float64Array(nm);
  for (let i = 0; i < nm; i++) {
    const md = modes[i];
    const tip = md.stiffness / (mPhi[i] || 1);
    sMM[i] = md.momentMudline ?? tip * (zHub - zMud); sVM[i] = md.shearMudline ?? tip;
    sMB[i] = md.momentBase ?? tip * (zHub - zTB); sVB[i] = md.shearBase ?? tip;
  }
  // ---- waves ----
  const wv = o.waves || null;
  const wetIdx = []; const dryIdx = [];
  stations.forEach((s, i) => { if ((s.medium ?? (s.z < 0 ? 'water' : 'air')) === 'water') wetIdx.push(i); else dryIdx.push(i); });
  let n2 = 2;                                  // power of two >= 2n + 2 samples at dt/2 (fast FFT path; the wave record is not periodic within the run)
  while (n2 < 2 * n + 2) n2 *= 2;
  let wave = null; let raoF = null;
  if (wv) {
    need(wv.depth > 0, 'waves.depth', 'must be > 0');
    wave = linearIrregularWave({ Hs: wv.Hs, Tp: wv.Tp, gamma: wv.gamma, depth: wv.depth, dt: dt / 2, n: n2, seed: wv.seed ?? 1,
      z: wetIdx.map((i) => stations[i].z), diameters: wv.macCamyFuchs ? wetIdx.map((i) => stations[i].D) : undefined, fMax: wv.fMax });
    if (fl && fl.forceRao) raoF = raoForces(fl.forceRao, wv, dt / 2, n2);
  }
  const current = wv ? (wv.current ?? 0) : 0;
  const waveAcc = wave ? (wave.accMF || wave.acc) : null;
  const shearFac = Float64Array.from(stations, (s) => Math.pow(Math.max(0.01, s.z / tb.hubHeight), o.shearExp ?? 0));
  const cIn = new Float64Array(ns); const cDr = new Float64Array(ns);
  for (let s = 0; s < ns; s++) {
    const stn = stations[s];
    const wet = wetIdx.includes(s);
    cIn[s] = wet ? rhoW * Math.PI * stn.D * stn.D / 4 * stn.dz * (wv && wv.macCamyFuchs ? 1 : (stn.cm ?? 2)) : 0;
    cDr[s] = 0.5 * (wet ? rhoW * (stn.cd ?? 1) : rho * (stn.cd ?? 0.7)) * stn.D * stn.dz;
  }
  // ---- floater ----
  let Minv = null; let Kp = null; let Bp = null; let Bq = null;
  if (fl) {
    const Mt = as3(fl.mass, 'floater.mass'); const Am = as3(fl.addedMass, 'floater.addedMass');
    Minv = inv3(Mt.map((r, i) => r.map((v, j) => v + Am[i][j])));
    Kp = as3(fl.stiffness, 'floater.stiffness'); Bp = as3(fl.damping, 'floater.damping');
    Bq = fl.quadDamping ? (typeof fl.quadDamping[0] === 'number' ? fl.quadDamping : [fl.quadDamping[0][0], fl.quadDamping[1][1], fl.quadDamping[2][2]]) : [0, 0, 0];
  }
  const hHubP = fl ? tb.hubHeight - (fl.refZ ?? 0) : 0;
  // ---- state ----
  const shaft = tb.shaft || null;
  const init = o.initial || {};
  const omega0 = init.omega ?? (cfg ? cfg.ratedGenSpeed / N : 0);
  const act = { theta: init.pitchDeg !== undefined ? init.pitchDeg * D2R : (cfg ? cfg.pitchMin : 0), rate: 0 };
  const ctrl = cfg ? createController(cfg, { mode: init.mode ?? 'production', pitch: act.theta, genSpeed: omega0 * N }) : null;
  // layout: [omega, omegaGenLss, twist, q(nm), qd(nm), xp(3), xpd(3), vInd]
  const iQ = 3; const iQd = 3 + nm; const iP = 3 + 2 * nm; const iPd = iP + 3; const iV = iPd + 3;
  const nState = iV + 1;
  const y = new Float64Array(nState);
  y[0] = omega0; y[1] = omega0; y[2] = 0;
  for (let i = 0; i < nm; i++) { y[iQ + i] = init.q ? init.q[i] : 0; y[iQd + i] = init.qd ? init.qd[i] : 0; }
  for (let i = 0; i < 3; i++) { y[iP + i] = init.platform ? init.platform[i] : 0; y[iPd + i] = init.platformVel ? init.platformVel[i] : 0; }
  const dynIn = !!o.dynamicInflow;
  y[iV] = init.vInd ?? -1;
  const k1 = new Float64Array(nState); const k2 = new Float64Array(nState); const k3 = new Float64Array(nState); const k4 = new Float64Array(nState);
  const yt = new Float64Array(nState);
  const Fst = new Float64Array(ns);
  const fGen = new Float64Array(nm);
  const ax = { thrust: 0, torque: 0, hubAcc: 0, hubVel: 0, waveShear: 0, waveMomMud: 0, waveMomTB: 0, waveShearTB: 0 };
  let genTq = 0; let brakeTq = 0; let pitch = act.theta;
  const Ir = shaft ? tb.inertia - shaft.generatorInertiaLss : tb.inertia;
  if (shaft) need(Ir > 0 && shaft.generatorInertiaLss > 0 && shaft.stiffness > 0, 'turbine.shaft', 'stiffness and inertias must be > 0');
  const over = o.aeroOverride || null;

  function deriv(yy, dy, j2, V) {
    const om = yy[0];
    let hubVel = 0;
    for (let i = 0; i < nm; i++) hubVel += mPhi[i] * yy[iQd + i];
    if (fl) hubVel += yy[iPd] + hHubP * yy[iPd + 2];
    const Vrel = V * cosYaw - hubVel;
    const aV = Math.abs(Vrel);
    let T = 0; let Q = 0; let vqs = 0; let aInd = 0;
    if (over) {
      const qd = 0.5 * rho * area * Vrel * aV;
      T = over.ct * qd; Q = (over.cq ?? 0) * qd * R;
    } else if (aV > 1e-3) {
      let Ve = Vrel;
      if (dynIn) {
        lookupCtCq(sf, om * R / aV, pitch / D2R);
        aInd = inductionFromCt(Math.min(2, Math.max(0, LK.ct)));
        vqs = aInd * aV;
        if (yy[iV] < 0) yy[iV] = vqs;
        Ve = Vrel + Math.sign(Vrel) * (vqs - yy[iV]) / Math.max(0.3, 1 - aInd);
      }
      const ae = Math.abs(Ve);
      lookupCtCq(sf, om * R / Math.max(ae, 1e-3), pitch / D2R);
      const qd = 0.5 * rho * area * Ve * ae;
      T = LK.ct * qd; Q = LK.cq * qd * R * Math.sign(Ve);
    }
    // drivetrain
    if (!ctrl) { dy[0] = 0; dy[1] = 0; dy[2] = 0; } else if (shaft) {
      const tw = shaft.stiffness * yy[2] + shaft.damping * (yy[0] - yy[1]);
      dy[0] = (Q - tw) / Ir;
      dy[1] = (tw - N * (genTq + (yy[1] > 0 ? brakeTq : 0))) / shaft.generatorInertiaLss;
      dy[2] = yy[0] - yy[1];
      if (yy[1] <= 0 && dy[1] < 0) dy[1] = 0;
    } else {
      dy[0] = (Q - N * (genTq + (om > 0 ? brakeTq : 0))) / tb.inertia;
      if (om <= 0 && dy[0] < 0) dy[0] = 0;
      dy[1] = dy[0]; dy[2] = 0;
    }
    // station loads
    let wSh = 0; let wMm = 0; let wMb = 0; let wSb = 0;
    for (let i = 0; i < nm; i++) fGen[i] = mPhi[i] * T;
    if (ns) {
      for (let w = 0; w < wetIdx.length; w++) {
        const s = wetIdx[w];
        let xs = 0;
        for (let i = 0; i < nm; i++) xs += phiS[i * ns + s] * yy[iQd + i];
        const u = (wave ? wave.u[w][j2] : 0) + current - xs;
        Fst[s] = (wave ? cIn[s] * waveAcc[w][j2] : 0) + cDr[s] * Math.abs(u) * u;
      }
      for (let d = 0; d < dryIdx.length; d++) {
        const s = dryIdx[d];
        let F = 0;
        if (o.towerWindDrag) {
          let xs = 0;
          for (let i = 0; i < nm; i++) xs += phiS[i * ns + s] * yy[iQd + i];
          const u = V * cosYaw * shearFac[s] - xs;
          F = cDr[s] * Math.abs(u) * u;
        }
        Fst[s] = F;
      }
      for (let s = 0; s < ns; s++) {
        const F = Fst[s]; const z = stZ[s];
        if (F === 0) continue;
        for (let i = 0; i < nm; i++) fGen[i] += phiS[i * ns + s] * F;
        wSh += F; wMm += F * (z - zMud);
        if (z > zTB) { wSb += F; wMb += F * (z - zTB); }
      }
    }
    let hubAcc = 0;
    for (let i = 0; i < nm; i++) {
      dy[iQ + i] = yy[iQd + i];
      dy[iQd + i] = (fGen[i] - mC[i] * yy[iQd + i] - mK[i] * yy[iQ + i]) / mM[i];
      hubAcc += mPhi[i] * dy[iQd + i];
    }
    if (fl) {
      const f0 = T + (raoF ? raoF[0][j2] : 0); const f1 = raoF ? raoF[1][j2] : 0; const f2 = T * hHubP + (raoF ? raoF[2][j2] : 0);
      const r = [f0, f1, f2];
      for (let a = 0; a < 3; a++) {
        const v = yy[iPd + a];
        r[a] -= Bq[a] * Math.abs(v) * v;
        for (let b = 0; b < 3; b++) r[a] -= Kp[a][b] * yy[iP + b] + Bp[a][b] * yy[iPd + b];
      }
      for (let a = 0; a < 3; a++) {
        dy[iP + a] = yy[iPd + a];
        dy[iPd + a] = Minv[a][0] * r[0] + Minv[a][1] * r[1] + Minv[a][2] * r[2];
      }
      hubAcc += dy[iPd] + hHubP * dy[iPd + 2];
    } else { for (let a = 0; a < 6; a++) dy[iP + a] = 0; }
    if (dynIn && !over && aV > 1e-3) {
      const tau = 1.1 / (1 - 1.3 * Math.min(0.5, aInd)) * R / Math.max(aV, 1);
      dy[iV] = (vqs - yy[iV]) / tau;
    } else dy[iV] = 0;
    ax.thrust = T; ax.torque = Q; ax.hubAcc = hubAcc; ax.hubVel = hubVel;
    ax.waveShear = wSh; ax.waveMomMud = wMm; ax.waveMomTB = wMb; ax.waveShearTB = wSb;
  }

  const N1 = n + 1;
  const F64 = () => new Float64Array(N1);
  const out = {
    t: F64(), wind: F64(), waveElevation: F64(), rotorSpeed: F64(), genSpeed: F64(), pitchDeg: F64(), genTorque: F64(), power: F64(),
    thrust: F64(), aeroTorque: F64(), towerTopDisp: F64(), towerTopVel: F64(), towerTopAcc: F64(),
    towerBaseMoment: F64(), towerBaseShear: F64(), mudlineMoment: F64(), mudlineShear: F64(),
    modal: modes.map(() => F64()), platform: fl ? { surge: F64(), heave: F64(), pitch: F64() } : null, mode: [],
  };
  const actArgs = { omegaAct: TWO_PI, zetaAct: 0.7, ...(o.actuator || {}), cmd: 0, dt, rateLimit: cfg ? Math.max(cfg.pitchRate, cfg.estopPitchRate ?? 0) : 1, min: cfg ? Math.min(0, cfg.pitchMin) : -Math.PI, max: cfg ? cfg.pitchMax : Math.PI };
  const cin = { dt, genSpeed: 0, pitch: 0, towerTopVel: 0, windSpeed: undefined, command: undefined, gridLoss: false };
  const rnaM = fl ? (fl.rnaMass ?? 0) : 0;
  for (let k = 0; k <= n; k++) {
    const t = k * dt;
    const V = sampleInput(o.wind, k, t, dt);
    need(Number.isFinite(V), 'wind', `non-finite value at t = ${t}`);
    pitch = act.theta;
    if (ctrl) {
      const ev = o.events ? (o.events(t) || {}) : {};
      let hv = 0;
      for (let i = 0; i < nm; i++) hv += mPhi[i] * y[iQd + i];
      cin.genSpeed = (shaft ? y[1] : y[0]) * N; cin.pitch = act.theta; cin.towerTopVel = hv; cin.command = ev.command; cin.gridLoss = !!ev.gridLoss;
      const c = ctrl.step(cin);
      genTq = c.genTorque; brakeTq = c.brakeTorque;
      actArgs.cmd = c.pitchCmd;
      out.mode.push(c.mode);
    } else out.mode.push('prescribed');
    deriv(y, k1, 2 * k, V);
    // ---- outputs at the start of the step ----
    out.t[k] = t; out.wind[k] = V; out.waveElevation[k] = wave ? wave.eta[2 * k] : 0;
    out.rotorSpeed[k] = y[0]; out.genSpeed[k] = (shaft ? y[1] : y[0]) * N; out.pitchDeg[k] = act.theta / D2R; out.genTorque[k] = genTq;
    out.power[k] = eta * genTq * out.genSpeed[k]; out.thrust[k] = ax.thrust; out.aeroTorque[k] = ax.torque;
    let disp = 0; let dMB = 0; let dVB = 0; let dMM = 0; let dVM = 0;
    for (let i = 0; i < nm; i++) {
      const q = y[iQ + i];
      out.modal[i][k] = q;
      disp += mPhi[i] * q;
      const dq = q - fGen[i] / mK[i];           // mode-acceleration: dynamic part beyond the static response
      dMB += sMB[i] * dq; dVB += sVB[i] * dq; dMM += sMM[i] * dq; dVM += sVM[i] * dq;
    }
    if (fl) {
      out.platform.surge[k] = y[iP]; out.platform.heave[k] = y[iP + 1]; out.platform.pitch[k] = y[iP + 2];
      disp += y[iP] + hHubP * y[iP + 2];
      // rigid tower on the platform: RNA inertia and weight offset about the tower base
      const hT = tb.hubHeight - zTB;
      out.towerBaseShear[k] = ax.thrust - rnaM * ax.hubAcc;
      out.towerBaseMoment[k] = (ax.thrust - rnaM * ax.hubAcc) * hT + rnaM * G * hT * Math.sin(y[iP + 2]);
      out.mudlineMoment[k] = out.towerBaseMoment[k]; out.mudlineShear[k] = out.towerBaseShear[k];
    } else {
      out.towerBaseMoment[k] = ax.thrust * (zHub - zTB) + ax.waveMomTB + dMB;
      out.towerBaseShear[k] = ax.thrust + ax.waveShearTB + dVB;
      out.mudlineMoment[k] = ax.thrust * (zHub - zMud) + ax.waveMomMud + dMM;
      out.mudlineShear[k] = ax.thrust + ax.waveShear + dVM;
    }
    out.towerTopDisp[k] = disp; out.towerTopVel[k] = ax.hubVel; out.towerTopAcc[k] = ax.hubAcc;
    if (k === n) break;
    // ---- RK4 ----
    const Vh = sampleInput(o.wind, k, t + 0.5 * dt, dt); const Vn = sampleInput(o.wind, k + 1, t + dt, dt);
    for (let i = 0; i < nState; i++) yt[i] = y[i] + 0.5 * dt * k1[i];
    deriv(yt, k2, 2 * k + 1, Vh);
    for (let i = 0; i < nState; i++) yt[i] = y[i] + 0.5 * dt * k2[i];
    deriv(yt, k3, 2 * k + 1, Vh);
    for (let i = 0; i < nState; i++) yt[i] = y[i] + dt * k3[i];
    deriv(yt, k4, 2 * k + 2, Vn);
    for (let i = 0; i < nState; i++) y[i] += dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    if (y[0] < 0) y[0] = 0;
    if (y[1] < 0) y[1] = 0;
    if (ctrl) pitchActuatorStep(act, actArgs);
  }
  // ---- post-processing ----
  const i0 = Math.min(n, Math.round((o.tTransient ?? 0) / dt));
  const stats = {};
  for (const key of ['wind', 'waveElevation', 'rotorSpeed', 'pitchDeg', 'power', 'thrust', 'towerTopDisp', 'towerTopAcc', 'towerBaseMoment', 'towerBaseShear', 'mudlineMoment', 'mudlineShear']) {
    stats[key] = seriesStats(out[key], i0);
  }
  if (fl) for (const key of ['surge', 'heave', 'pitch']) stats[`platform.${key}`] = seriesStats(out.platform[key], i0);
  const win = (n - i0) * dt;
  const dOpt = o.del || {};
  const dArgs = { m: dOpt.m ?? 4, nEq: dOpt.nEq ?? Math.max(1, win), rainflow: dOpt.rainflow };
  out.del = {
    towerBaseMoment: damageEquivalentLoad(out.towerBaseMoment.subarray(i0), dArgs).del,
    mudlineMoment: damageEquivalentLoad(out.mudlineMoment.subarray(i0), dArgs).del,
    m: dArgs.m, nEq: dArgs.nEq,
  };
  out.stats = stats;
  const nfft = Math.min(o.psdNfft ?? 1024, n - i0);
  if (nfft >= 64) {
    const sp = psdWelch(out.towerBaseMoment.subarray(i0), 1 / dt, { nfft });
    const cands = [];
    const om = stats.rotorSpeed.mean;
    if (om > 0.01) { cands.push({ name: '1P', f: om / TWO_PI }); cands.push({ name: '3P', f: 3 * om / TWO_PI }); }
    modes.forEach((md, i) => cands.push({ name: `mode ${i + 1}`, f: Math.sqrt(md.stiffness / md.mass) / TWO_PI }));
    if (wv) cands.push({ name: 'wave peak', f: 1 / wv.Tp });
    out.psd = { f: sp.f, towerBaseMoment: sp.psd, peaks: identifyPeaks(sp.f, sp.psd, cands) };
  } else out.psd = null;
  out.final = {
    omega: y[0], pitchDeg: act.theta / D2R, q: Array.from(y.subarray(iQ, iQ + nm)), qd: Array.from(y.subarray(iQd, iQd + nm)),
    platform: Array.from(y.subarray(iP, iP + 3)), platformVel: Array.from(y.subarray(iPd, iPd + 3)), vInd: y[iV],
    mode: ctrl ? ctrl.state.mode : 'prescribed',
  };
  out.structure = str;
  return out;
}

function raoForces(rao, wv, dt, n) {
  need(Array.isArray(rao.f) && rao.f.length >= 2, 'floater.forceRao.f', 'at least two frequencies required');
  const T = n * dt; const df = 1 / T;
  const fMax = Math.min(0.5 / dt, wv.fMax ?? Math.max(0.5, 4 / wv.Tp));
  const mMax = Math.min(n / 2 - 1, Math.floor(fMax / df));
  // identical component amplitudes and phases as linearIrregularWave for the same seed
  const rng = makeRng(wv.seed ?? 1);
  let sum = 0;
  const S = new Float64Array(mMax + 1);
  for (let m = 1; m <= mMax; m++) { S[m] = jonswapShape(m * df, wv.Tp, wv.gamma ?? 3.3); sum += S[m] * df; }
  const scale = sum > 0 ? wv.Hs * wv.Hs / 16 / sum : 0;
  const lin = (xs, ys, x) => {
    if (x <= xs[0]) return ys[0];
    if (x >= xs[xs.length - 1]) return ys[ys.length - 1];
    let i = 0;
    while (xs[i + 1] < x) i++;
    return ys[i] + (ys[i + 1] - ys[i]) * (x - xs[i]) / (xs[i + 1] - xs[i]);
  };
  const keys = ['surge', 'heave', 'pitch'];
  const re = keys.map(() => new Float64Array(n)); const im = keys.map(() => new Float64Array(n));
  for (let m = 1; m <= mMax; m++) {
    const amp = Math.sqrt(2 * S[m] * scale * df);
    const ph = TWO_PI * rng.next();
    keys.forEach((key, a) => {
      const r = rao[key];
      if (!r) return;
      const A = lin(rao.f, r.amp, m * df);
      const p = ph + (r.phaseDeg ? lin(rao.f, r.phaseDeg, m * df) * D2R : 0);
      const cr = 0.5 * amp * A * n * Math.cos(p); const ci = 0.5 * amp * A * n * Math.sin(p);
      re[a][m] = cr; im[a][m] = ci; re[a][n - m] = cr; im[a][n - m] = -ci;
    });
  }
  return keys.map((_, a) => Float64Array.from(ifft(re[a], im[a]).re));
}

/**
 * Aerodynamic damping of a support-structure mode with the rotor operating.
 * Two estimates:
 *  - `zetaAnalytic`: from the thrust slope at frozen rotor speed and pitch, c_aero = phiHub^2 dT/dV, zeta =
 *    c_aero / (2 sqrt(k m)) (Salzmann and van der Tempel 2005 closed-loop effects excluded);
 *  - `zetaDecay`: from a free-decay simulation: the model is run to steady state in constant wind, the mode is
 *    displaced by `perturbation` and the logarithmic decrement of the decay gives the total damping; the structural
 *    damping ratio of the mode is subtracted.
 * @param {object} o simulateCoupled arguments with constant `wind` and no waves
 * @param {object} [p]
 * @param {number} [p.modeIndex=0] @param {number} [p.perturbation=0.1] hub displacement step [m]
 * @param {number} [p.settle=120] settling time before the perturbation [s] @param {number} [p.decay] decay record length [s]
 * @returns {{zetaDecay:number, zetaTotal:number, zetaStructural:number, zetaAnalytic:number, frequency:number, dTdV:number}}
 *   damping ratios [-], damped frequency [Hz], thrust slope [N s/m]
 */
export function aeroDampingRatio(o, p = {}) {
  need(typeof o.wind === 'number', 'wind', 'constant wind speed required');
  const mi = p.modeIndex ?? 0;
  const base = { ...o, waves: null };
  const s1 = simulateCoupled({ ...base, duration: p.settle ?? 120 });
  const md = s1.structure.modes[mi];
  need(md, 'modeIndex', 'mode not found');
  const fn = Math.sqrt(md.stiffness / md.mass) / TWO_PI;
  const init = { ...s1.final, q: s1.final.q.slice() };
  init.q[mi] += (p.perturbation ?? 0.1) / (md.phiHub || 1);
  const s2 = simulateCoupled({ ...base, structure: s1.structure, support: undefined, initial: init, duration: p.decay ?? 12 / fn });
  const ld = logDecrement(s2.modal[mi], o.dt);
  // thrust slope at the settled operating point
  const V = o.wind; const om = s1.final.omega; const th = s1.final.pitchDeg;
  const A = Math.PI * o.surface.R ** 2; const rho = o.rho ?? 1.225;
  const Tof = (v) => (o.aeroOverride ? o.aeroOverride.ct : surfaceLookup(o.surface, om * o.surface.R / v, th).ct) * 0.5 * rho * A * v * v;
  const h = 0.05 * V;
  const dTdV = (Tof(V + h) - Tof(V - h)) / (2 * h);
  return {
    zetaTotal: ld.zeta, zetaStructural: md.zeta ?? 0, zetaDecay: ld.zeta - (md.zeta ?? 0),
    zetaAnalytic: (md.phiHub ?? 1) ** 2 * dTdV / (2 * Math.sqrt(md.stiffness * md.mass)), frequency: ld.frequency, dTdV,
  };
}
