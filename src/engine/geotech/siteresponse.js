/**
 * One-dimensional site response: vertically propagating SH waves in a layered visco-elastic column,
 * equivalent-linear (SHAKE-type) strain-compatible iteration, response spectra and intensity measures.
 *
 * SI units: m, s, kg/m^3, Pa, m/s, m/s^2, Hz; damping ratios and strains as fractions.
 *
 * References
 *   Schnabel P.B., Lysmer J. & Seed H.B. (1972) SHAKE, Report EERC 72-12, UC Berkeley.
 *   Idriss I.M. & Sun J.I. (1992) SHAKE91 user's manual.
 *   Kramer S.L. (1996) Geotechnical Earthquake Engineering, §7.2 (transfer functions, eqs 7.8–7.34).
 *   Nigam N.C. & Jennings P.C. (1969) Calculation of response spectra from strong-motion earthquake
 *     records, BSSA 59(2); Chopra A.K., Dynamics of Structures, Table 5.2.1.
 *   Arias A. (1970) A measure of earthquake intensity; Trifunac M.D. & Brady A.G. (1975) BSSA 65(3).
 */
import { fft, ifft } from '../math/fft.js';
import { modulusReduction } from './soilmodels.js';

const G0 = 9.80665;

function need(name, v, min = -Infinity, max = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < min || v > max) {
    throw new RangeError(`${name} must be a number in [${min}, ${max}] (got ${v})`);
  }
  return v;
}

function prepLayers(layers) {
  if (!Array.isArray(layers) || layers.length === 0) throw new RangeError('layers must be a non-empty array');
  return layers.map((l, i) => {
    need(`layers[${i}].h`, l.h, 1e-9);
    need(`layers[${i}].rho`, l.rho, 1e-9);
    const Gmax = l.Gmax !== undefined ? need(`layers[${i}].Gmax`, l.Gmax, 1e-9) : l.rho * need(`layers[${i}].Vs`, l.Vs, 1e-9) ** 2;
    return { h: l.h, rho: l.rho, Gmax, G: Gmax, damping: l.damping === undefined ? 0 : need(`layers[${i}].damping`, l.damping, 0, 0.99), curve: l.curve };
  });
}

/*
 * Wave amplitudes for unit A1 = B1 = 1 (Kramer eqs 7.31–7.34) at circular frequency w.
 * Complex shear-wave velocity Vs* = Vs √(1 + 2iξ); k* = w / Vs*; impedance ratio α = ρm Vs*m / (ρm+1 Vs*m+1).
 * Fills A, B (re/im pairs) for layers 0..n (n = half-space when present) and k* per layer.
 */
function amplitudes(w, L, rock, ws) {
  const n = L.length;
  ws.Ar[0] = 1; ws.Ai[0] = 0; ws.Br[0] = 1; ws.Bi[0] = 0;
  for (let m = 0; m < n; m++) {
    const lay = L[m];
    const vs = Math.sqrt(lay.G / lay.rho);
    // s = √(1 + 2iξ)
    const mod = Math.sqrt(Math.sqrt(1 + 4 * lay.damping * lay.damping));
    const arg = 0.5 * Math.atan2(2 * lay.damping, 1);
    const sr = mod * Math.cos(arg);
    const si = mod * Math.sin(arg);
    const s2 = sr * sr + si * si;
    const kr = (w * sr) / (vs * s2);
    const ki = (-w * si) / (vs * s2);
    ws.kr[m] = kr; ws.ki[m] = ki;
    // impedance Z = ρ Vs s
    const zr = lay.rho * vs * sr;
    const zi = lay.rho * vs * si;
    let ar = 0;
    let ai = 0;
    const next = m + 1 < n ? L[m + 1] : rock;
    if (next) {
      const vn = Math.sqrt(next.G / next.rho);
      const modn = Math.sqrt(Math.sqrt(1 + 4 * next.damping * next.damping));
      const argn = 0.5 * Math.atan2(2 * next.damping, 1);
      const nr = next.rho * vn * modn * Math.cos(argn);
      const ni = next.rho * vn * modn * Math.sin(argn);
      const d = nr * nr + ni * ni;
      ar = (zr * nr + zi * ni) / d;
      ai = (zi * nr - zr * ni) / d;
    }
    // e^{+ik h} and e^{−ik h}
    const ep = Math.exp(-ki * lay.h);
    const c = Math.cos(kr * lay.h);
    const sn = Math.sin(kr * lay.h);
    const pr = ep * c, pi = ep * sn;
    const mr = c / ep, mi = -sn / ep;
    // A e^{ikh}, B e^{−ikh}
    const aer = ws.Ar[m] * pr - ws.Ai[m] * pi;
    const aei = ws.Ar[m] * pi + ws.Ai[m] * pr;
    const ber = ws.Br[m] * mr - ws.Bi[m] * mi;
    const bei = ws.Br[m] * mi + ws.Bi[m] * mr;
    ws.Er[m] = aer + ber; ws.Ei[m] = aei + bei; // displacement at the layer base
    // A_{m+1} = ½[(1+α) Ae + (1−α) Be],  B_{m+1} = ½[(1−α) Ae + (1+α) Be]
    const sumr = aer + ber, sumi = aei + bei;
    const difr = aer - ber, difi = aei - bei;
    const adr = ar * difr - ai * difi;
    const adi = ar * difi + ai * difr;
    ws.Ar[m + 1] = 0.5 * (sumr + adr); ws.Ai[m + 1] = 0.5 * (sumi + adi);
    ws.Br[m + 1] = 0.5 * (sumr - adr); ws.Bi[m + 1] = 0.5 * (sumi - adi);
  }
}

function workspace(n) {
  const f = () => new Float64Array(n + 2);
  return { Ar: f(), Ai: f(), Br: f(), Bi: f(), kr: f(), ki: f(), Er: f(), Ei: f() };
}

function prepRock(bedrock) {
  if (!bedrock || bedrock.rigid) return null;
  need('bedrock.rho', bedrock.rho, 1e-9);
  need('bedrock.Vs', bedrock.Vs, 1e-9);
  return { rho: bedrock.rho, G: bedrock.rho * bedrock.Vs ** 2, damping: bedrock.damping || 0 };
}

/**
 * Frequency-domain transfer function of a layered visco-elastic column (surface motion / input motion).
 *   'within'  : input is the motion at the top of the bedrock (base of the column). For a rigid base this
 *               is the only meaningful choice; a uniform undamped layer gives |F| = 1/|cos(ωH/Vs)| with
 *               resonances at f_n = (2n − 1) Vs / 4H.
 *   'outcrop' : input is the motion the bedrock would have as a free outcrop (2 × incident wave).
 * @param {object} o
 * @param {{h:number,rho:number,Vs?:number,Gmax?:number,damping?:number}[]} o.layers top to bottom
 * @param {{rho:number,Vs:number,damping?:number}|{rigid:true}} [o.bedrock] elastic half-space or rigid base
 * @param {ArrayLike<number>} o.freqs frequencies [Hz]
 * @param {'outcrop'|'within'} [o.input] default 'outcrop' for elastic bedrock, 'within' for rigid
 * @returns {{f:number[],amp:number[],re:number[],im:number[]}} complex transfer function and modulus [-]
 */
export function transferFunction({ layers, bedrock, freqs, input }) {
  const L = prepLayers(layers);
  const rock = prepRock(bedrock);
  const mode = input || (rock ? 'outcrop' : 'within');
  if (!rock && mode === 'outcrop') throw new RangeError("input 'outcrop' needs an elastic bedrock");
  const ws = workspace(L.length);
  const n = L.length;
  const out = { f: [], amp: [], re: [], im: [] };
  for (const f of freqs) {
    let re = 1;
    let im = 0;
    if (f > 0) {
      amplitudes(2 * Math.PI * f, L, rock, ws);
      // surface displacement = A1 + B1 = 2
      const dr = mode === 'outcrop' ? 2 * ws.Ar[n] : ws.Er[n - 1];
      const di = mode === 'outcrop' ? 2 * ws.Ai[n] : ws.Ei[n - 1];
      const d = dr * dr + di * di;
      re = (2 * dr) / d;
      im = (-2 * di) / d;
    }
    out.f.push(f); out.re.push(re); out.im.push(im); out.amp.push(Math.hypot(re, im));
  }
  return out;
}

/**
 * Equivalent-linear one-dimensional site response (SHAKE procedure). The accelerogram is transformed to
 * the frequency domain, propagated through the column with complex moduli G* = G (1 + 2iξ), and the
 * strain-compatible G and ξ of each layer are iterated with effective strain γeff = strainRatio · γmax at
 * mid-layer until the properties change by less than `tol`.
 * Assumptions: vertically propagating SH waves, horizontal layers, frequency-independent hysteretic
 * damping, total-stress analysis (no pore-pressure generation). Equivalent-linear analysis becomes
 * unreliable for effective strains above roughly 0.5–1 %.
 * @param {object} o
 * @param {{h:number,rho:number,Vs?:number,Gmax?:number,damping?:number,curve?:object}[]} o.layers top to
 *   bottom; `curve` is a modulus-reduction model for soilmodels.modulusReduction (omit for a linear layer
 *   with the given damping)
 * @param {{rho:number,Vs:number,damping?:number}|{rigid:true}} [o.bedrock]
 * @param {ArrayLike<number>} o.acc input accelerogram [m/s^2]
 * @param {number} o.dt time step [s]
 * @param {'outcrop'|'within'} [o.input] how the input motion is specified (see transferFunction)
 * @param {number} [o.strainRatio=0.65] effective / maximum shear strain [-] ((M − 1)/10 is also used)
 * @param {number} [o.tol=0.01] relative convergence tolerance on G and ξ [-]
 * @param {number} [o.maxIter=15]
 * @param {boolean} [o.layerHistories=false] also return the PGA at every layer top (one FFT per layer)
 * @returns {{surfaceAcc:Float64Array,t:Float64Array,pgaInput:number,pgaSurface:number,amplification:number,
 *   layers:{G:number,damping:number,ratio:number,gammaMax:number,gammaEff:number,zMid:number,pgaTop?:number}[],
 *   iterations:number,converged:boolean,tf:{f:Float64Array,amp:Float64Array}}}
 *   surfaceAcc [m/s^2]; PGA values [m/s^2]; per-layer strain-compatible properties; transfer function
 *   modulus (surface / input) up to the Nyquist frequency
 */
export function equivalentLinearResponse({ layers, bedrock, acc, dt, input, strainRatio = 0.65, tol = 0.01, maxIter = 15, layerHistories = false }) {
  const L = prepLayers(layers);
  const rock = prepRock(bedrock);
  need('dt', dt, 1e-9);
  need('strainRatio', strainRatio, 0.1, 1);
  if (!acc || acc.length < 4) throw new RangeError('acc must contain at least four samples');
  const mode = input || (rock ? 'outcrop' : 'within');
  if (!rock && mode === 'outcrop') throw new RangeError("input 'outcrop' needs an elastic bedrock");
  const nt = acc.length;
  let N = 1;
  while (N < nt) N *= 2;
  const x = new Float64Array(N);
  for (let i = 0; i < nt; i++) x[i] = acc[i];
  const X = fft(x);
  const nL = L.length;
  const ws = workspace(nL);
  const half = N / 2;
  const dw = (2 * Math.PI) / (N * dt);
  const sr = new Float64Array(N);
  const si = new Float64Array(N);
  // initial damping for nonlinear layers: small-strain value of the curve
  for (const l of L) if (l.curve) { const c = modulusReduction(1e-6, l.curve); l.G = l.Gmax * c.ratio; l.damping = c.damping; }
  const gammaMax = new Float64Array(nL);
  let iterations = 0;
  let converged = false;
  const tfAmp = new Float64Array(half + 1);
  const Yr = new Float64Array(N);
  const Yi = new Float64Array(N);
  const anyNonlinear = L.some((l) => l.curve);

  // surface acceleration spectrum for the current properties
  const solveSurface = () => {
    Yr[0] = X.re[0]; Yi[0] = X.im[0]; tfAmp[0] = 1;
    for (let j = 1; j <= half; j++) {
      amplitudes(j * dw, L, rock, ws);
      const dr = mode === 'outcrop' ? 2 * ws.Ar[nL] : ws.Er[nL - 1];
      const di = mode === 'outcrop' ? 2 * ws.Ai[nL] : ws.Ei[nL - 1];
      const d = dr * dr + di * di;
      const tr = (2 * dr) / d;
      const ti = (-2 * di) / d;
      tfAmp[j] = Math.hypot(tr, ti);
      Yr[j] = X.re[j] * tr - X.im[j] * ti;
      Yi[j] = X.re[j] * ti + X.im[j] * tr;
      if (j < half) { Yr[N - j] = Yr[j]; Yi[N - j] = -Yi[j]; } else Yi[j] = 0;
    }
  };
  // time history of (cA·A_m e^{ikz} + cB·B_m e^{−ikz}) × factor(ω) × surface spectrum / 2
  const layerQuantity = (m, z, strain) => {
    sr[0] = strain ? 0 : Yr[0]; si[0] = 0;
    for (let j = 1; j <= half; j++) {
      const w = j * dw;
      amplitudes(w, L, rock, ws);
      const kr = ws.kr[m], ki = ws.ki[m];
      const ep = Math.exp(-ki * z);
      const c = Math.cos(kr * z);
      const sn = Math.sin(kr * z);
      const aer = (ws.Ar[m] * c - ws.Ai[m] * sn) * ep;
      const aei = (ws.Ar[m] * sn + ws.Ai[m] * c) * ep;
      const ber = (ws.Br[m] * c + ws.Bi[m] * sn) / ep;
      const bei = (-ws.Br[m] * sn + ws.Bi[m] * c) / ep;
      let qr, qi;
      if (strain) {
        // γ = ik (A e^{ikz} − B e^{−ikz}) · U_surface/2,  U_surface = −Acc_surface/ω²
        const dr = aer - ber, di = aei - bei;
        const tr = -(kr * di + ki * dr);
        const ti = kr * dr - ki * di;
        const sc = -0.5 / (w * w);
        qr = tr * sc; qi = ti * sc;
      } else {
        qr = 0.5 * (aer + ber); qi = 0.5 * (aei + bei);
      }
      sr[j] = Yr[j] * qr - Yi[j] * qi;
      si[j] = Yr[j] * qi + Yi[j] * qr;
      if (j < half) { sr[N - j] = sr[j]; si[N - j] = -si[j]; } else si[j] = 0;
    }
    const y = ifft(sr, si).re;
    let mx = 0;
    for (let i = 0; i < N; i++) { const a = Math.abs(y[i]); if (a > mx) mx = a; }
    return mx;
  };

  for (let it = 0; it < Math.max(maxIter, 1); it++) {
    iterations = it + 1;
    solveSurface();
    let change = 0;
    for (let m = 0; m < nL; m++) {
      gammaMax[m] = layerQuantity(m, L[m].h / 2, true);
      if (L[m].curve) {
        const c = modulusReduction(Math.max(strainRatio * gammaMax[m], 1e-9), L[m].curve);
        const Gn = L[m].Gmax * c.ratio;
        change = Math.max(change, Math.abs(Gn - L[m].G) / L[m].G, Math.abs(c.damping - L[m].damping) / Math.max(L[m].damping, 0.01));
        L[m].G = Gn;
        L[m].damping = c.damping;
      }
    }
    if (!anyNonlinear || change < tol) { converged = true; if (anyNonlinear) { solveSurface(); for (let m = 0; m < nL; m++) gammaMax[m] = layerQuantity(m, L[m].h / 2, true); } break; }
  }
  const ys = ifft(Yr, Yi).re;
  const surfaceAcc = new Float64Array(nt);
  const t = new Float64Array(nt);
  let pgaS = 0;
  let pgaI = 0;
  for (let i = 0; i < nt; i++) { surfaceAcc[i] = ys[i]; t[i] = i * dt; pgaS = Math.max(pgaS, Math.abs(ys[i])); pgaI = Math.max(pgaI, Math.abs(acc[i])); }
  // the zero-padded tail can carry the response peak of lightly damped columns
  for (let i = nt; i < N; i++) pgaS = Math.max(pgaS, Math.abs(ys[i]));
  const outLayers = [];
  let zTop = 0;
  for (let m = 0; m < nL; m++) {
    const o = { G: L[m].G, damping: L[m].damping, ratio: L[m].G / L[m].Gmax, gammaMax: gammaMax[m], gammaEff: strainRatio * gammaMax[m], zMid: zTop + L[m].h / 2 };
    if (layerHistories) o.pgaTop = m === 0 ? pgaS : layerQuantity(m, 0, false);
    outLayers.push(o);
    zTop += L[m].h;
  }
  const f = new Float64Array(half + 1);
  for (let j = 0; j <= half; j++) f[j] = (j * dw) / (2 * Math.PI);
  return { surfaceAcc, t, pgaInput: pgaI, pgaSurface: pgaS, amplification: pgaI > 0 ? pgaS / pgaI : NaN, layers: outLayers, iterations, converged, tf: { f, amp: tfAmp } };
}

/**
 * Elastic response spectrum of an accelerogram by the exact piecewise-linear recurrence of Nigam &
 * Jennings (1969) (coefficients as in Chopra, Table 5.2.1, for unit mass and p = −a_g):
 *   u_{i+1} = A u_i + B v_i + C p_i + D p_{i+1},   v_{i+1} = A′ u_i + B′ v_i + C′ p_i + D′ p_{i+1}.
 * The recurrence is exact for ground acceleration varying linearly between samples, for any Δt/T.
 * @param {object} o
 * @param {ArrayLike<number>} o.acc ground acceleration [m/s^2]
 * @param {number} o.dt [s]
 * @param {ArrayLike<number>} o.periods natural periods [s], > 0
 * @param {number} [o.damping=0.05] damping ratio [-], 0 ≤ ξ < 1
 * @param {number} [o.freeVibrationCycles=0] cycles of free vibration appended after the record [-]
 * @returns {{T:number[],Sd:number[],Sv:number[],Sa:number[],SaAbs:number[],PSV:number[],PSA:number[]}}
 *   Sd relative displacement [m]; Sv relative velocity [m/s]; SaAbs absolute acceleration [m/s^2];
 *   PSV = ω Sd [m/s]; PSA = ω² Sd [m/s^2]; Sa = PSA (alias)
 */
export function responseSpectrum({ acc, dt, periods, damping = 0.05, freeVibrationCycles = 0 }) {
  need('dt', dt, 1e-12);
  need('damping', damping, 0, 0.999);
  if (!acc || acc.length < 2) throw new RangeError('acc must contain at least two samples');
  const out = { T: [], Sd: [], Sv: [], Sa: [], SaAbs: [], PSV: [], PSA: [] };
  const xi = damping;
  const n = acc.length;
  for (const T of periods) {
    need('period', T, 1e-9);
    const w = (2 * Math.PI) / T;
    const k = w * w;
    const r = Math.sqrt(1 - xi * xi);
    const wd = w * r;
    const e = Math.exp(-xi * w * dt);
    const s = Math.sin(wd * dt);
    const c = Math.cos(wd * dt);
    const A = e * ((xi / r) * s + c);
    const B = (e * s) / wd;
    const C = (1 / k) * ((2 * xi) / (w * dt) + e * (((1 - 2 * xi * xi) / (wd * dt) - xi / r) * s - (1 + (2 * xi) / (w * dt)) * c));
    const D = (1 / k) * (1 - (2 * xi) / (w * dt) + e * (((2 * xi * xi - 1) / (wd * dt)) * s + ((2 * xi) / (w * dt)) * c));
    const A1 = -e * ((w / r) * s);
    const B1 = e * (c - (xi / r) * s);
    const C1 = (1 / k) * (-1 / dt + e * ((w / r + xi / (dt * r)) * s + c / dt));
    const D1 = (1 / (k * dt)) * (1 - e * ((xi / r) * s + c));
    let u = 0;
    let v = 0;
    let sd = 0;
    let sv = 0;
    let sa = 0;
    const extra = Math.ceil((freeVibrationCycles * T) / dt);
    for (let i = 0; i < n - 1 + extra; i++) {
      const p0 = i < n ? -acc[i] : 0;
      const p1 = i + 1 < n ? -acc[i + 1] : 0;
      const un = A * u + B * v + C * p0 + D * p1;
      const vn = A1 * u + B1 * v + C1 * p0 + D1 * p1;
      u = un; v = vn;
      const au = Math.abs(u);
      if (au > sd) sd = au;
      const av = Math.abs(v);
      if (av > sv) sv = av;
      const aa = Math.abs(2 * xi * w * v + k * u);
      if (aa > sa) sa = aa;
    }
    out.T.push(T); out.Sd.push(sd); out.Sv.push(sv); out.SaAbs.push(sa); out.PSV.push(w * sd); out.PSA.push(k * sd); out.Sa.push(k * sd);
  }
  return out;
}

/**
 * Arias intensity Ia = (π / 2g) ∫ a² dt and significant duration (time between 5 % and 95 % of Ia by
 * default; Trifunac & Brady 1975), trapezoidal integration.
 * @param {object} o
 * @param {ArrayLike<number>} o.acc [m/s^2]  @param {number} o.dt [s]
 * @param {number} [o.lower=0.05]  @param {number} [o.upper=0.95] Husid fractions [-]
 * @returns {{arias:number,significantDuration:number,t5:number,t95:number,pga:number,husid:Float64Array}}
 *   arias [m/s]; durations [s]; husid = normalised cumulative intensity [-]
 */
export function ariasIntensity({ acc, dt, lower = 0.05, upper = 0.95 }) {
  need('dt', dt, 1e-12);
  need('lower', lower, 0, 1);
  need('upper', upper, lower, 1);
  if (!acc || acc.length < 2) throw new RangeError('acc must contain at least two samples');
  const n = acc.length;
  const cum = new Float64Array(n);
  let pga = Math.abs(acc[0]);
  for (let i = 1; i < n; i++) {
    cum[i] = cum[i - 1] + 0.5 * (acc[i - 1] * acc[i - 1] + acc[i] * acc[i]) * dt;
    pga = Math.max(pga, Math.abs(acc[i]));
  }
  const total = cum[n - 1];
  const husid = new Float64Array(n);
  if (total > 0) for (let i = 0; i < n; i++) husid[i] = cum[i] / total;
  const cross = (frac) => {
    if (!(total > 0)) return 0;
    for (let i = 1; i < n; i++) {
      if (husid[i] >= frac) {
        const d = husid[i] - husid[i - 1];
        return (i - 1 + (d > 0 ? (frac - husid[i - 1]) / d : 0)) * dt;
      }
    }
    return (n - 1) * dt;
  };
  const t5 = cross(lower);
  const t95 = cross(upper);
  return { arias: (Math.PI / (2 * G0)) * total, significantDuration: t95 - t5, t5, t95, pga, husid };
}
