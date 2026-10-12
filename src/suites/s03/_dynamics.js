// Shared structural-dynamics helpers of the suite: the stick model of monopile + tower with its fore-aft modes,
// thrust slopes of the rotor, and the frequency-domain response of the first mode to turbulent thrust and to
// linear wave inertia loading (used by the dynamic-response and fatigue modules).
import { macCamyFuchsCm, waveNumber, jonswapShape } from '../../engine/aero/coupled.js';
import { kaimalParams, kaimalSpectrum, rotorAveragingFilter } from '../../engine/aero/turbulence.js';
import { surfaceLookup } from '../../engine/aero/bem.js';
import { run, need, segments, foundationOf, operatingPoint } from './_common.js';

const TWO_PI = 2 * Math.PI;

/** Hydrodynamic and damping inputs shared by the dynamic modules. */
export function dynamicsInputs({ group = 'Dynamics' } = {}) {
  return [
    { key: 'zetaS', label: 'Structural + soil + hydrodynamic damping', unit: '%', type: 'number', default: 1.2, min: 0.1, max: 10, group, role: 'ic', uncertain: { dist: 'lognormal', cov: 0.3 },
      help: 'Damping ratio of the first fore-aft mode without the rotor: steel about 0.2–0.3 %, soil 0.4–1.0 %, hydrodynamic and sloshing 0.1–0.3 % (measured on monopiles in rotor-stop tests). Aerodynamic damping is added by the model.' },
    { key: 'Ca', label: 'Added-mass coefficient of the wet pile', unit: '-', type: 'number', default: 1, min: 0, max: 2, group, role: 'environment' },
    { key: 'flooded', label: 'Monopile flooded inside', type: 'bool', default: true, group, role: 'ic', help: 'Enclosed water moves with the pile and adds mass below the still-water level.' },
    { key: 'Cm', label: 'Morison inertia coefficient', unit: '-', type: 'number', default: 2, min: 1, max: 2.5, group, role: 'environment', help: 'Used when the MacCamy–Fuchs diffraction correction is switched off.' },
    { key: 'Cd', label: 'Morison drag coefficient', unit: '-', type: 'number', default: 0.9, min: 0.4, max: 1.6, group, role: 'environment', help: '0.65 smooth to 1.05 rough cylinder (DNV-RP-C205); marine growth raises it.' },
    { key: 'mcf', label: 'MacCamy–Fuchs diffraction correction', type: 'bool', default: true, group, role: 'environment', help: 'Reduces the inertia load of short waves on a large-diameter pile (D / wavelength > 0.2).' },
  ];
}

/**
 * Stick model of monopile + tower and its lowest fore-aft modes, in the structure format of the coupled simulator
 * (same idealisation as the engine's buildSupportStructure: two-node Euler–Bernoulli elements with consistent mass,
 * lumped rotor–nacelle mass, coupled lateral/rocking springs or a clamp at the mudline, added mass and enclosed
 * water on wet elements). The eigenproblem is solved by inverse iteration on the Cholesky-factored banded stiffness,
 * which costs a few milliseconds; the dense QL path of the engine builder takes between 0.02 s and 2 s on these
 * badly scaled beam matrices and is avoided here.
 * @param {object} o segments (SI, z up from still water), rnaMass kg, zHub m, foundation ('fixed' | { Kh, Kr, Khr }),
 *   nElem, nModes, zeta, rhoWater, addedMassCa, internalWater, cd, cm, cdAir
 * @returns {object} { zHub, zMudline, zTowerBase, zTop, modes[], stations[], nodes, mass } as buildSupportStructure
 */
export function stickModel(o) {
  const segs = o.segments;
  const zMud = segs[0].zBase, zTop = segs[segs.length - 1].zTop, total = zTop - zMud;
  const rhoW = o.rhoWater ?? 1025, Ca = o.addedMassCa ?? 1;
  const zN = [zMud], els = [];
  for (const sg of segs) {
    const L = sg.zTop - sg.zBase;
    const cuts = [sg.zBase, sg.zTop];
    if (sg.zBase < -1e-9 && sg.zTop > 1e-9) cuts.splice(1, 0, 0); // still-water level kept as a node
    for (let c = 0; c < cuts.length - 1; c++) {
      const Lc = cuts[c + 1] - cuts[c];
      const ne = Math.max(1, Math.round((o.nElem * Lc) / total));
      for (let e = 0; e < ne; e++) {
        const za = cuts[c] + (Lc * e) / ne, zb = cuts[c] + (Lc * (e + 1)) / ne;
        const q = (0.5 * (za + zb) - sg.zBase) / L;
        const D = sg.dBase + (sg.dTop - sg.dBase) * q, t = sg.tBase + (sg.tTop - sg.tBase) * q;
        const A = Math.PI * (D - t) * t, I = (Math.PI / 64) * (D ** 4 - (D - 2 * t) ** 4);
        const wet = 0.5 * (za + zb) < 0;
        let m = sg.density * A;
        if (wet) m += rhoW * Ca * Math.PI * D * D / 4 + ((o.internalWater ?? true) ? rhoW * Math.PI * (D - 2 * t) ** 2 / 4 : 0);
        els.push({ za, zb, D, EI: sg.E * I, m, mSteel: sg.density * A, wet });
        zN.push(zb);
      }
    }
  }
  const nn = zN.length, nd = 2 * nn;
  const K0 = new Float64Array(nd * nd), M = new Float64Array(nd * nd);
  let massStruct = 0;
  els.forEach((el, e) => {
    const L = el.zb - el.za, k = el.EI / L ** 3, c = (el.m * L) / 420;
    const ke = [12, 6 * L, -12, 6 * L, 6 * L, 4 * L * L, -6 * L, 2 * L * L, -12, -6 * L, 12, -6 * L, 6 * L, 2 * L * L, -6 * L, 4 * L * L];
    const me = [156, 22 * L, 54, -13 * L, 22 * L, 4 * L * L, 13 * L, -3 * L * L, 54, 13 * L, 156, -22 * L, -13 * L, -3 * L * L, -22 * L, 4 * L * L];
    for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) { K0[(2 * e + a) * nd + 2 * e + b] += k * ke[4 * a + b]; M[(2 * e + a) * nd + 2 * e + b] += c * me[4 * a + b]; }
    massStruct += el.mSteel * L;
  });
  M[(nd - 2) * nd + nd - 2] += o.rnaMass;
  M[(nd - 1) * nd + nd - 1] += o.rnaInertia ?? 0;
  const fixed = o.foundation === undefined || o.foundation === 'fixed';
  const K = Float64Array.from(K0);
  if (!fixed) { const f = o.foundation; K[0] += f.Kh; K[nd + 1] += f.Kr; K[1] += f.Khr ?? 0; K[nd] += f.Khr ?? 0; }
  const off = fixed ? 2 : 0, n = nd - off;
  // Cholesky factor of the restrained stiffness (dense lower triangle; n is a few hundred at most)
  const Lc = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      let sum = K[(i + off) * nd + j + off];
      for (let q = Math.max(0, i - 4, j - 4); q < j; q++) sum -= Lc[i * n + q] * Lc[j * n + q]; // half-bandwidth 3
      if (i === j) { need(sum > 0, 'The stick model is not restrained (stiffness matrix not positive definite): check the foundation stiffness and section sizes.'); Lc[i * n + i] = Math.sqrt(sum); } else Lc[i * n + j] = sum / Lc[j * n + j];
    }
  }
  const buf = new Float64Array(n), mb = new Float64Array(n);
  const solveK = (b) => {
    for (let i = 0; i < n; i++) { let sum = b[i]; for (let q = Math.max(0, i - 3); q < i; q++) sum -= Lc[i * n + q] * buf[q]; buf[i] = sum / Lc[i * n + i]; }
    for (let i = n - 1; i >= 0; i--) { let sum = buf[i]; for (let q = i + 1; q < Math.min(n, i + 4); q++) sum -= Lc[q * n + i] * buf[q]; buf[i] = sum / Lc[i * n + i]; }
    return buf;
  };
  const mulM = (x) => {
    for (let i = 0; i < n; i++) { let sum = 0; const r = (i + off) * nd + off; for (let j = Math.max(0, i - 3); j < Math.min(n, i + 4); j++) sum += M[r + j] * x[j]; mb[i] = sum; }
    return mb;
  };
  // Lowest modes by inverse iteration with M-orthogonal deflation: the beam spectrum is well separated, so each
  // mode converges in a few tens of sweeps of two banded triangular solves.
  const nModes = Math.min(o.nModes ?? 2, n);
  const eig = { values: [], vectors: [], converged: true };
  const mv = [];
  for (let q = 0; q < nModes; q++) {
    let x = Float64Array.from({ length: n }, (_, i) => 1 + ((i * 7919 + q * 104729) % 97) / 97);
    let lam = 0, done = false;
    for (let it = 0; it < 400 && !done; it++) {
      for (let p = 0; p < q; p++) { let d = 0; for (let i = 0; i < n; i++) d += mv[p][i] * x[i]; for (let i = 0; i < n; i++) x[i] -= d * eig.vectors[p][i]; }
      const y = Float64Array.from(solveK(mulM(x)));
      for (let p = 0; p < q; p++) { let d = 0; for (let i = 0; i < n; i++) d += mv[p][i] * y[i]; for (let i = 0; i < n; i++) y[i] -= d * eig.vectors[p][i]; }
      const My = mulM(y);
      let yMy = 0, yKy = 0;
      for (let i = 0; i < n; i++) yMy += y[i] * My[i];
      const Mx = mulM(x);
      for (let i = 0; i < n; i++) yKy += y[i] * Mx[i]; // yᵀK y = yᵀM x because K y = M x
      const lamNew = yKy / yMy, sc = 1 / Math.sqrt(yMy);
      for (let i = 0; i < n; i++) x[i] = y[i] * sc;
      done = it > 2 && Math.abs(lamNew - lam) <= 1e-10 * Math.abs(lamNew);
      lam = lamNew;
    }
    if (!done) eig.converged = false;
    eig.values.push(lam); eig.vectors.push(x); mv.push(Float64Array.from(mulM(x)));
  }
  const zTB = o.zTowerBase ?? (segs.length > 1 ? segs[0].zTop : zMud);
  const iTB = zN.findIndex((z) => Math.abs(z - zTB) < 1e-6);
  const stations = els.map((el) => ({ z: 0.5 * (el.za + el.zb), D: el.D, dz: el.zb - el.za, medium: el.wet ? 'water' : 'air', cd: el.wet ? (o.cd ?? 1.0) : (o.cdAir ?? 0.7), cm: o.cm ?? 2.0 }));
  const modes = [];
  for (let q = 0; q < nModes; q++) {
    const v = new Float64Array(nd);
    for (let i = off; i < nd; i++) v[i] = eig.vectors[q][i - off];
    const hub = v[nd - 2] + v[nd - 1] * (o.zHub - zTop);
    need(Math.abs(hub) > 1e-12, 'A mode with no hub displacement was found: reduce the number of modes.');
    for (let i = 0; i < nd; i++) v[i] /= hub;
    let mg = 0, lr = 0;
    const f = new Float64Array(nd);
    for (let i = 0; i < nd; i++) {
      let sm = 0, sk = 0;
      for (let j = Math.max(0, i - 3); j < Math.min(nd, i + 4); j++) { sm += M[i * nd + j] * v[j]; sk += K0[i * nd + j] * v[j]; }
      mg += v[i] * sm; f[i] = sk;
      if (i % 2 === 0) lr += sm; // φᵀ M r with r = unit lateral translation of every node
    }
    // section forces per unit modal coordinate from the free-body above the section (elastic nodal forces)
    const section = (iNode) => { let V = 0, Mo = 0; for (let j = iNode + 1; j < nn; j++) { V += f[2 * j]; Mo += f[2 * j] * (zN[j] - zN[iNode]) + f[2 * j + 1]; } return [Mo, V]; };
    const sm0 = section(0), sb = section(Math.max(iTB, 0));
    const lam = eig.values[q];
    modes.push({ frequency: Math.sqrt(Math.max(0, lam)) / TWO_PI, mass: mg, stiffness: lam * mg, zeta: o.zeta ?? 0.01, phiHub: 1, slopeHub: v[nd - 1],
      phiStations: els.map((el, e) => 0.5 * (v[2 * e] + v[2 * e + 2]) + ((el.zb - el.za) / 8) * (v[2 * e + 1] - v[2 * e + 3])),
      momentMudline: sm0[0], shearMudline: sm0[1], momentBase: sb[0], shearBase: sb[1], participation: lr / mg, effectiveMass: (lr * lr) / mg, shape: Float64Array.from({ length: nn }, (_, i) => v[2 * i]) });
  }
  let mTot = o.rnaMass;
  els.forEach((el) => { mTot += el.m * (el.zb - el.za); });
  return { zHub: o.zHub, zMudline: zMud, zTowerBase: zTB, zTop, modes, stations, nodes: Float64Array.from(zN), mass: massStruct, totalMass: mTot, converged: eig.converged };
}

/** Stick model from the shared inputs. */
export function supportStructure(inp, tb, steel, { nElem = 30, nModes = 2 } = {}) {
  return run(() => stickModel({
    segments: segments(inp, tb, steel), rnaMass: tb.rnaMass, zHub: inp.hubHeight, foundation: foundationOf(inp), nElem, nModes,
    zeta: (inp.zetaS ?? 1.2) / 100, rhoWater: inp.rhoW ?? 1025, addedMassCa: inp.Ca ?? 1, internalWater: inp.flooded ?? true, cd: inp.Cd ?? 0.9, cm: inp.Cm ?? 2,
  }), 'Support-structure model');
}

/** Thrust slopes dT/dV [N·s/m] at mean wind V: frozen rotor speed and pitch (fast changes) and along the schedule (slow changes). */
export function thrustSlopes(tb, V, parked = false) {
  const op = parked ? { thrust: 0, omega: 0, pitchDeg: 90, ct: 0, power: 0, parked: true } : operatingPoint(tb, V);
  if (op.parked) {
    // parked rotor: T = Ct ½ρAV², slope 2T/V in both limits
    const s = tb.parkedCt * tb.rhoAir * tb.area * V;
    return { frozen: s, schedule: s, op: { ...op, thrust: 0.5 * s * V } };
  }
  const T = (v) => surfaceLookup(tb.surface, (op.omega * tb.R) / v, op.pitchDeg).ct * 0.5 * tb.rhoAir * tb.area * v * v;
  const h = 0.03 * V;
  const frozen = (T(V + h) - T(V - h)) / (2 * h);
  const a = operatingPoint(tb, Math.min(V + 0.5, tb.schedule.U[tb.schedule.U.length - 1])), b = operatingPoint(tb, Math.max(V - 0.5, tb.schedule.U[0]));
  const schedule = (a.thrust - b.thrust) / Math.max(1e-9, Math.min(V + 0.5, tb.schedule.U[tb.schedule.U.length - 1]) - Math.max(V - 0.5, tb.schedule.U[0]));
  return { frozen: Math.max(frozen, 0), schedule, op };
}

// The rotor-averaging filter depends on f/V only (for a given rotor radius and coherence scale): it is tabulated
// once per turbine on a logarithmic reduced-frequency grid and interpolated, which removes its cost from sweeps.
const H2CACHE = new Map();
function rotorFilter(tb, Lc, f, V) {
  const key = `${tb.R}|${Lc}`;
  if (!H2CACHE.has(key)) {
    const xs = Float64Array.from({ length: 241 }, (_, i) => 10 ** (-5 + (5.3 * i) / 240));
    H2CACHE.set(key, rotorAveragingFilter(Array.from(xs, (x) => x * 10), { R: tb.R, Vhub: 10, Lc }));
  }
  const h = H2CACHE.get(key);
  const out = new Float64Array(f.length);
  for (let i = 0; i < f.length; i++) {
    const q = (Math.log10(Math.max(f[i] / V, 1e-5)) + 5) / (5.3 / 240);
    const j = Math.min(Math.max(Math.floor(q), 0), 239), w = Math.min(Math.max(q - j, 0), 1);
    out[i] = h[j] + w * (h[j + 1] - h[j]);
  }
  return out;
}

/**
 * Frequency grid [Hz] from 1/600 Hz (10-minute record) to fMax, graded quadratically towards low frequency where
 * the turbulence spectrum is steep, plus a cluster across the resonance peak so that light damping is resolved.
 */
export function frequencyGrid(n, fMax, f1, zeta) {
  const set = [];
  const f0 = 1 / 600;
  const nu = Math.max(8, Math.round(0.75 * n)), nc = Math.max(8, n - nu);
  for (let i = 0; i <= nu; i++) set.push(f0 + (fMax - f0) * (i / nu) ** 2);
  const half = Math.max(8 * zeta, 0.05) * f1;
  for (let i = 0; i < nc; i++) set.push(f1 - half + (2 * half * (i + 0.5)) / nc);
  return Float64Array.from(set.filter((f) => f >= f0).sort((p, q) => p - q));
}
const trapz = (f, S) => { let s = 0; for (let i = 1; i < f.length; i++) s += 0.5 * (S[i] + S[i - 1]) * (f[i] - f[i - 1]); return s; };
export { trapz };

/**
 * First-mode response spectra under turbulent thrust and linear wave inertia loading.
 * Section forces use the mode-acceleration form: quasi-static part from the external load directly plus the
 * dynamic part of the modal coordinate. Wave drag is not linearised (inertia-dominated large pile).
 * @returns {object} f [Hz]; SM, SMwind, SMwave mudline-moment PSD [(N·m)²/Hz]; Sx hub displacement PSD [m²/Hz];
 *   Sa hub acceleration PSD [(m/s²)²/Hz]; sigma values; zetaAero, zetaTotal [-]; f1 [Hz]; mean thrust [N] and mean
 *   mudline moment [N·m]
 */
export function modalSpectra({ str, tb, inp, V, sigma1, Hs, Tp, gamma = 3.3, nFreq = 400, parked = false }) {
  const md = str.modes[0];
  const f1 = md.frequency, k = md.stiffness, m = md.mass;
  const sl = thrustSlopes(tb, V, parked);
  const zetaAero = sl.op.parked ? 0 : (md.phiHub ** 2 * sl.frozen) / (2 * Math.sqrt(k * m));
  const zeta = md.zeta + zetaAero;
  const depth = inp.depth, rhoW = inp.rhoW ?? 1025;
  const lever = str.zHub - str.zMudline;
  const f = frequencyGrid(nFreq, Math.max(1.0, 2.5 * f1), f1, zeta);
  const n = f.length;
  // wind
  const kp = kaimalParams({ hubHeight: inp.hubHeight, sigma1 });
  const H2 = sigma1 > 0 ? rotorFilter(tb, kp.Lc, f, V) : new Float64Array(n);
  const fc = 0.6 / TWO_PI; // pitch-loop design frequency of the baseline controller, 0.6 rad/s
  // waves
  const wet = [];
  str.stations.forEach((s, i) => { if (s.medium === 'water') wet.push(i); });
  let shape = 0;
  const Sj = new Float64Array(n);
  if (Hs > 0) { for (let i = 0; i < n; i++) Sj[i] = jonswapShape(f[i], Tp, gamma); shape = trapz(f, Sj); }
  const scale = shape > 0 ? (Hs * Hs) / 16 / shape : 0;
  const out = { f, SM: new Float64Array(n), SMwind: new Float64Array(n), SMwave: new Float64Array(n), Sx: new Float64Array(n), Sa: new Float64Array(n), SV: new Float64Array(n) };
  for (let i = 0; i < n; i++) {
    const w = TWO_PI * f[i], r = f[i] / f1;
    const den = (1 - r * r) ** 2 + (2 * zeta * r) ** 2;
    const hr = (1 - r * r) / den, hi = (-2 * zeta * r) / den;
    const resp = (Fq, Mext) => {
      const c = (md.momentMudline * Fq) / k;
      return { M2: (Mext + c * (hr - 1)) ** 2 + (c * hi) ** 2, x2: ((md.phiHub * Fq) / k) ** 2 * (hr * hr + hi * hi) };
    };
    if (sigma1 > 0) {
      const x = (f[i] / fc) ** 2;
      const slope = sl.schedule + (sl.frozen - sl.schedule) * (x / (1 + x));
      const Su = kaimalSpectrum(f[i], { Vhub: V, sigma: sigma1, L: kp.u.L }) * H2[i];
      const q = resp(slope * md.phiHub, slope * lever);
      out.SMwind[i] = q.M2 * Su; out.Sx[i] += q.x2 * Su;
      out.SV[i] += slope * slope * Su;
    }
    if (scale > 0) {
      const kw = waveNumber(w, depth), kd = kw * depth;
      let Fq = 0, Mext = 0, Vext = 0, dLast = -1, cmLast = 2;
      for (const s of wet) {
        const st = str.stations[s];
        const z = Math.min(0, st.z);
        const Gz = kd > 30 ? Math.exp(kw * z) : Math.cosh(kw * (z + depth)) / Math.sinh(kd);
        if (st.D !== dLast) { dLast = st.D; cmLast = inp.mcf ? macCamyFuchsCm(0.5 * kw * st.D).magnitude : (inp.Cm ?? 2); }
        const cm = cmLast;
        const fz = rhoW * cm * (Math.PI * st.D * st.D) / 4 * w * w * Gz * st.dz; // force amplitude per unit wave amplitude
        Fq += md.phiStations[s] * fz; Mext += fz * (st.z - str.zMudline); Vext += fz;
      }
      const Se = Sj[i] * scale;
      const q = resp(Fq, Mext);
      out.SMwave[i] = q.M2 * Se; out.Sx[i] += q.x2 * Se;
      out.SV[i] += Vext * Vext * Se;
    }
    out.SM[i] = out.SMwind[i] + out.SMwave[i];
    out.Sa[i] = out.Sx[i] * w ** 4;
  }
  need(Number.isFinite(out.SM[n - 1]), 'The response spectrum could not be evaluated: check the wind and wave inputs.');
  const meanThrust = sl.op.thrust;
  return {
    ...out, f1, zetaAero, zetaTotal: zeta, slopes: sl, meanThrust, meanMoment: meanThrust * lever,
    sigmaM: Math.sqrt(trapz(f, out.SM)), sigmaMwind: Math.sqrt(trapz(f, out.SMwind)), sigmaMwave: Math.sqrt(trapz(f, out.SMwave)),
    sigmaX: Math.sqrt(trapz(f, out.Sx)), sigmaA: Math.sqrt(trapz(f, out.Sa)), sigmaV: Math.sqrt(trapz(f, out.SV)),
    nu0: Math.sqrt(trapz(f, Float64Array.from(out.SM, (s, i) => s * f[i] * f[i])) / Math.max(trapz(f, out.SM), 1e-300)),
  };
}
