/**
 * Slug-flow characteristics: translational velocity, slug-body holdup, frequency, length statistics, unit-cell
 * model, synthetic slug trains at the outlet, hydrodynamic slug initiation (Kelvin–Helmholtz), riser-induced
 * severe slugging criteria and cycle estimate, surge-volume and pigging-slug sizing.
 *
 * theta is the pipe inclination from horizontal [rad], positive for upward flow.
 */
import { brent } from '../math/roots.js';
import { makeRng } from '../math/rng.js';
import { dist } from '../math/stats.js';
import { G, checkArg, bendiksen, darcyFrictionFactor } from './driftflux.js';
import { stratifiedGeometry, stratifiedForces, levelFromHoldup } from './flowregime.js';

/** Universal gas constant [J/(mol·K)]. */
export const R_GAS = 8.314462618;

const area = (D) => Math.PI * D * D / 4;

/**
 * Slug (Taylor-bubble) translational velocity vt = C0 vm + vd.
 * Reference: Nicklin, Wilkes & Davidson (1962), Trans. Inst. Chem. Eng. 40, 61–68 (vertical: C0 = 1.2,
 * vd = 0.35 √(gD)); Bendiksen (1984), Int. J. Multiphase Flow 10(4), 467–483 (all inclinations, default).
 * @param {object} o
 * @param {number} o.vm mixture velocity [m/s]
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.theta=0] inclination [rad]
 * @param {number} [o.C0] override of the distribution coefficient [-]
 * @param {number} [o.vd] override of the drift velocity [m/s]
 * @returns {{vt: number, C0: number, vd: number}} vt, vd [m/s]; C0 [-]
 */
export function slugTranslationalVelocity({ vm, D, theta = 0, C0, vd }) {
  const b = bendiksen({ vm, D, theta });
  const c0 = C0 !== undefined ? checkArg('C0', C0, { min: 0.5, max: 2 }) : b.C0;
  const v0 = vd !== undefined ? checkArg('vd', vd) : b.vd;
  return { vt: c0 * vm + v0, C0: c0, vd: v0 };
}

/**
 * Liquid holdup in the slug body, HLS = 1 / (1 + (vm/8.66)^1.39), vm in m/s.
 * Reference: Gregory, Nicholson & Aziz (1978), Int. J. Multiphase Flow 4(1), 33–39. Validity: horizontal
 * air–light-oil data, D = 25.8 and 51.2 mm, vm < ~10 m/s.
 * @param {number} vm mixture velocity [m/s], >= 0
 * @returns {number} slug-body liquid holdup [-]
 */
export function slugBodyHoldupGregory(vm) {
  checkArg('vm', vm, { min: 0 });
  return 1 / (1 + (vm / 8.66) ** 1.39);
}

/**
 * Slug-body holdup from the Barnea & Brauner balance between turbulent break-up and coalescence:
 *   αs = 0.058 { dc [ (2 fs/D) vm³ ]^(2/5) (ρL/σ)^(3/5) − 0.725 }²,  dc = 2 [0.4 σ / ((ρL − ρG) g)]^(1/2),
 * with fs the Fanning friction factor of the liquid at the mixture velocity. HLS = 1 − αs, limited to
 * 0 ≤ αs ≤ 0.52 (maximum packing of dispersed bubbles).
 * Reference: Barnea & Brauner (1985), Int. J. Multiphase Flow 11(1), 43–49.
 * @param {object} o
 * @param {number} o.vm mixture velocity [m/s]
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³]
 * @param {number} o.muL liquid viscosity [Pa·s]
 * @param {number} o.sigma surface tension [N/m]
 * @param {number} [o.roughness=0] wall roughness [m]
 * @returns {number} slug-body liquid holdup [-]
 */
export function slugBodyHoldupBarneaBrauner({ vm, D, rhoL, rhoG, muL, sigma, roughness = 0 }) {
  checkArg('vm', vm, { min: 0 });
  checkArg('D', D, { min: 0, exclusiveMin: true });
  checkArg('rhoL', rhoL, { min: 0, exclusiveMin: true });
  checkArg('rhoG', rhoG, { min: 0, max: rhoL });
  checkArg('sigma', sigma, { min: 0, exclusiveMin: true });
  const fs = darcyFrictionFactor(rhoL * vm * D / muL, roughness / D) / 4;
  const dc = 2 * Math.sqrt(0.4 * sigma / ((rhoL - rhoG) * G));
  const term = dc * (2 * fs / D * vm ** 3) ** 0.4 * (rhoL / sigma) ** 0.6 - 0.725;
  const alpha = term > 0 ? Math.min(0.52, 0.058 * term * term) : 0;
  return 1 - alpha;
}

/**
 * Slug frequency, Gregory & Scott: fs = 0.0226 [ vsl/(gD) · (19.75/vm + vm) ]^1.2 (SI; 19.75 in m²/s²).
 * Reference: Gregory & Scott (1969), AIChE J. 15(6), 933–935. Validity: horizontal CO2–water, D = 19 mm.
 * @param {object} o
 * @param {number} o.vsl superficial liquid velocity [m/s], >= 0
 * @param {number} o.vm mixture velocity [m/s], > 0
 * @param {number} o.D pipe inner diameter [m]
 * @returns {number} slug frequency [1/s]
 */
export function slugFrequencyGregoryScott({ vsl, vm, D }) {
  checkArg('vsl', vsl, { min: 0 });
  checkArg('vm', vm, { min: 0, exclusiveMin: true });
  checkArg('D', D, { min: 0, exclusiveMin: true });
  return 0.0226 * (vsl / (G * D) * (19.75 / vm + vm)) ** 1.2;
}

/**
 * Slug frequency with Zabaras' inclination correction:
 *   fs = 0.0226 [vsl/(gD)]^1.2 [212.6/vm + vm]^1.2 [0.836 + 2.75 sin^0.25 θ]    (ft, ft/s, g = 32.2 ft/s²).
 * The bracketed group vsl (c/vm + vm)/(gD) is dimensionless, so the SI form only changes the constant
 * 212.6 ft²/s² = 19.75 m²/s² (the Gregory–Scott constant); that is the form coded here.
 * Reference: Zabaras (2000), SPE Journal 5(3), 252–258. Validity: 0 ≤ θ ≤ 11°, D = 25–200 mm in the data bank.
 * @param {object} o
 * @param {number} o.vsl superficial liquid velocity [m/s], >= 0
 * @param {number} o.vm mixture velocity [m/s], > 0
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.theta=0] upward inclination [rad], 0 … π/2
 * @returns {number} slug frequency [1/s]
 */
export function slugFrequencyZabaras({ vsl, vm, D, theta = 0 }) {
  checkArg('vsl', vsl, { min: 0 });
  checkArg('vm', vm, { min: 0, exclusiveMin: true });
  checkArg('D', D, { min: 0, exclusiveMin: true });
  checkArg('theta', theta, { min: 0, max: Math.PI / 2 });
  return 0.0226 * (vsl / (G * D)) ** 1.2 * (19.75 / vm + vm) ** 1.2 * (0.836 + 2.75 * Math.sin(theta) ** 0.25);
}

/**
 * Mean slug length. 'diameters': Ls = factor · D (≈ 32 D horizontal after Dukler & Hubbard 1975 / Barnea &
 * Brauner 1985; 16–20 D vertical after Taitel, Barnea & Dukler 1980). 'scott': large-diameter Prudhoe Bay
 * correlation ln(Ls) = −25.4 + 28.5 [ln(D)]^0.1 with D in inches and Ls in feet.
 * Reference: Scott, Shoham & Brill (1989), SPE Prod. Eng. 4(3), 335–340 (valid for D ≳ 0.1 m; field data
 * 0.3–0.6 m).
 * @param {object} o
 * @param {number} o.D pipe inner diameter [m]
 * @param {'diameters'|'scott'} [o.method='diameters']
 * @param {number} [o.factor=32] slug length in diameters [-] (16–20 vertical, ~32 horizontal)
 * @returns {number} mean slug length [m]
 */
export function meanSlugLength({ D, method = 'diameters', factor = 32 }) {
  checkArg('D', D, { min: 0, exclusiveMin: true });
  if (method === 'diameters') return checkArg('factor', factor, { min: 0, exclusiveMin: true }) * D;
  if (method === 'scott') {
    const dIn = D / 0.0254;
    checkArg('D (inches, Scott et al. correlation)', dIn, { min: 1, exclusiveMin: true });
    return 0.3048 * Math.exp(-25.4 + 28.5 * Math.log(dIn) ** 0.1);
  }
  throw new RangeError(`method = ${method} is not a known slug-length method`);
}

/**
 * Lognormal slug-length distribution with a stated mean and log-standard deviation.
 * Brill et al. (1981, SPE J. 21(3), 363–378) found slug lengths at Prudhoe Bay lognormal with σ(ln Ls) ≈ 0.5;
 * the 0.1 % exceedance length is the usual design maximum.
 * @param {object} o
 * @param {number} o.meanLength mean slug length [m]
 * @param {number} [o.sigmaLn=0.5] standard deviation of ln(Ls) [-] (typ. 0.3–0.6)
 * @param {number} [o.exceedance=0.001] exceedance probability of the design maximum [-]
 * @returns {{mu: number, sigmaLn: number, median: number, maxLength: number, quantile: (p: number) => number,
 *   sample: (rng: {normal: () => number}) => number}} mu = mean of ln(Ls [m]) [-]; median, maxLength [m];
 *   quantile(p) [m] for cumulative probability p; sample(rng) draws a length [m]
 */
export function slugLengthDistribution({ meanLength, sigmaLn = 0.5, exceedance = 0.001 }) {
  checkArg('meanLength', meanLength, { min: 0, exclusiveMin: true });
  checkArg('sigmaLn', sigmaLn, { min: 0, max: 3 });
  checkArg('exceedance', exceedance, { min: 0, max: 1, exclusiveMin: true });
  const mu = Math.log(meanLength) - sigmaLn * sigmaLn / 2;
  const std = dist({ dist: 'normal', mean: 0, std: 1 });
  const quantile = (p) => Math.exp(mu + sigmaLn * std.icdf(p));
  return {
    mu, sigmaLn, median: Math.exp(mu), maxLength: quantile(1 - exceedance), quantile,
    sample: (rng) => Math.exp(mu + sigmaLn * rng.normal()),
  };
}

/**
 * Slug unit-cell model for horizontal and near-horizontal pipes (equilibrium-film version).
 *   translational velocity   vt = C0 vm + vd (Bendiksen 1984)
 *   slug body                holdup HLS (Gregory et al. 1978 unless given), no slip: vLS = vGS = vm
 *   continuity across cell   (vt − vm) HLS = (vt − vLF) HLF,  (vt − vm)(1 − HLS) = (vt − vGF)(1 − HLF)
 *   film zone                uniform stratified film satisfying the combined momentum balance
 *   liquid balance           vsl = vm HLS β + vLF HLF (1 − β),  β = Ls/Lu
 *   unit holdup              HLU = β HLS + (1 − β) HLF  [= (vt HLS + vm (1 − HLS) − vsg)/vt]
 *   pressure gradient        global momentum balance over the unit: slug-body wall friction over β, film-zone
 *                            wall friction over (1 − β), gravity on the unit-average density
 * Reference: Dukler & Hubbard (1975), Ind. Eng. Chem. Fundam. 14(4), 337–347; Taitel & Barnea (1990), Adv.
 * Heat Transfer 20, 83–132 (equilibrium-film simplification, global momentum balance).
 * Validity: |θ| ≲ 15°, fully developed slug flow; `feasible` is false when the balances admit no film
 * thinner than the slug body (the operating point is then not slug flow).
 * @param {object} o
 * @param {number} o.vsl superficial liquid velocity [m/s], > 0
 * @param {number} o.vsg superficial gas velocity [m/s], > 0
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.theta=0] inclination [rad]
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³]
 * @param {number} o.muL liquid viscosity [Pa·s]
 * @param {number} o.muG gas viscosity [Pa·s]
 * @param {number} [o.roughness=0] wall roughness [m]
 * @param {number} [o.slugLength=32 D] slug-body length Ls [m]
 * @param {number} [o.slugHoldup] slug-body holdup HLS [-] (default Gregory et al.)
 * @param {number} [o.C0] distribution coefficient override [-]
 * @param {number} [o.vd] drift velocity override [m/s]
 * @param {'smooth'|'constant'|'andritsos'} [o.interfacial='smooth'] film interfacial friction (see stratifiedForces)
 * @param {number} [o.fi] constant interfacial Fanning factor [-]
 * @returns {{feasible: boolean, vm: number, vt: number, C0: number, vd: number, slugHoldup: number,
 *   filmHoldup: number, filmLevel: number, vLF: number, vGF: number, slugFraction: number, slugLength: number,
 *   filmLength: number, unitLength: number, frequency: number, holdup: number, gradFrictionSlug: number,
 *   gradFrictionFilm: number, gradFriction: number, gradGravity: number, gradTotal: number,
 *   unitPressureDrop: number, residual: number}} velocities [m/s]; holdups, filmLevel = hF/D and
 *   slugFraction = Ls/Lu [-]; lengths [m]; frequency = vt/Lu [1/s]; holdup = unit-average liquid holdup [-];
 *   gradients = pressure loss per metre [Pa/m] (unit averages except the two zone values);
 *   unitPressureDrop over one slug unit [Pa]; residual of the film momentum balance [Pa/m]
 */
export function slugUnitCell({ vsl, vsg, D, theta = 0, rhoL, rhoG, muL, muG, roughness = 0, slugLength, slugHoldup,
  C0, vd, interfacial = 'smooth', fi }) {
  checkArg('vsl', vsl, { min: 0, exclusiveMin: true });
  checkArg('vsg', vsg, { min: 0, exclusiveMin: true });
  checkArg('D', D, { min: 0, exclusiveMin: true });
  const vm = vsl + vsg;
  const tv = slugTranslationalVelocity({ vm, D, theta, C0, vd });
  const vt = tv.vt;
  const HLS = slugHoldup !== undefined ? checkArg('slugHoldup', slugHoldup, { min: 0, max: 1, exclusiveMin: true })
    : slugBodyHoldupGregory(vm);
  const Ls = slugLength !== undefined ? checkArg('slugLength', slugLength, { min: 0, exclusiveMin: true }) : 32 * D;
  const fluid = { D, theta, rhoL, rhoG, muL, muG, roughness, interfacial, fi };
  const film = (h) => {
    const HLF = stratifiedGeometry(h).holdup;
    const vLF = vt - (vt - vm) * HLS / HLF;
    const vGF = vt - (vt - vm) * (1 - HLS) / (1 - HLF);
    return { HLF, vLF, vGF, forces: stratifiedForces({ ...fluid, h, vL: vLF, vG: vGF }) };
  };
  const rhoS = HLS * rhoL + (1 - HLS) * rhoG;
  const muS = HLS * muL + (1 - HLS) * muG;
  const gradFrictionSlug = darcyFrictionFactor(rhoS * vm * D / muS, roughness / D) * rhoS * vm * vm / (2 * D);
  const sinT = Math.sin(theta);

  // Search downward from the slug-body level for the first film level that balances shear and gravity.
  const hTop = levelFromHoldup(Math.min(HLS, 1 - 1e-9)) * (1 - 1e-9);
  let hRoot = null;
  if (vt > vm) {
    let hPrev = hTop;
    let rPrev = film(hPrev).forces.residual;
    if (rPrev < 0) {
      for (let k = 1; k <= 240; k++) {
        const h = hTop * (1e-5) ** (k / 240);
        const r = film(h).forces.residual;
        if (r >= 0) { hRoot = brent((x) => film(x).forces.residual, h, hPrev, { tol: 1e-13 }); break; }
        hPrev = h; rPrev = r;
      }
    }
  }
  if (hRoot === null) {
    // No bubble zone: the cell degenerates to a continuous aerated slug.
    const gradGravity = rhoS * G * sinT;
    return {
      feasible: false, vm, vt, C0: tv.C0, vd: tv.vd, slugHoldup: HLS, filmHoldup: HLS, filmLevel: hTop, vLF: vm, vGF: vm,
      slugFraction: 1, slugLength: Ls, filmLength: 0, unitLength: Ls, frequency: vt / Ls, holdup: HLS,
      gradFrictionSlug, gradFrictionFilm: gradFrictionSlug, gradFriction: gradFrictionSlug, gradGravity,
      gradTotal: gradFrictionSlug + gradGravity, unitPressureDrop: (gradFrictionSlug + gradGravity) * Ls, residual: NaN,
    };
  }
  const f = film(hRoot);
  let beta = (vsl - f.vLF * f.HLF) / (vm * HLS - f.vLF * f.HLF);
  const feasible = beta > 0 && beta <= 1;
  beta = Math.min(1, Math.max(1e-9, beta));
  const Lu = Ls / beta;
  const holdup = beta * HLS + (1 - beta) * f.HLF;
  const gradFrictionFilm = f.forces.gradFriction;
  const gradFriction = beta * gradFrictionSlug + (1 - beta) * gradFrictionFilm;
  const gradGravity = (holdup * rhoL + (1 - holdup) * rhoG) * G * sinT;
  return {
    feasible, vm, vt, C0: tv.C0, vd: tv.vd, slugHoldup: HLS, filmHoldup: f.HLF, filmLevel: hRoot, vLF: f.vLF, vGF: f.vGF,
    slugFraction: beta, slugLength: Ls, filmLength: Lu - Ls, unitLength: Lu, frequency: vt / Lu, holdup,
    gradFrictionSlug, gradFrictionFilm, gradFriction, gradGravity, gradTotal: gradFriction + gradGravity,
    unitPressureDrop: (gradFriction + gradGravity) * Lu, residual: f.forces.residual,
  };
}

/**
 * Synthetic slug train at the pipe outlet: lognormally distributed slug lengths (seeded), each followed by a
 * film/bubble zone whose length keeps the unit-cell liquid balance, convected past the outlet at the
 * translational velocity. While a slug body passes, qL = A vm HLS and qG = A vm (1 − HLS); while a film zone
 * passes, qL = A vLF HLF and qG = A vGF (1 − HLF).
 * Reference: unit-cell kinematics of Taitel & Barnea (1990); lognormal lengths after Brill et al. (1981).
 * @param {object} o
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} o.vm mixture velocity [m/s]
 * @param {number} o.vt translational velocity [m/s], > 0
 * @param {number} o.slugHoldup HLS [-]
 * @param {number} o.filmHoldup HLF [-]
 * @param {number} o.vLF film liquid velocity [m/s]
 * @param {number} o.vGF gas velocity above the film [m/s]
 * @param {number} o.slugFraction Ls/Lu [-], 0 < β ≤ 1
 * @param {number} o.meanSlugLength mean slug length [m]
 * @param {number} [o.sigmaLn=0.5] standard deviation of ln(Ls) [-]
 * @param {number} o.duration length of the history [s]
 * @param {number} [o.nSamples=1000] samples in the returned rate history [-]
 * @param {number} [o.seed=1] generator seed [-]
 * @returns {{slugs: {tFront: number, tTail: number, length: number, volume: number}[], t: Float64Array,
 *   qL: Float64Array, qG: Float64Array, meanQL: number, meanQG: number, peakQL: number, maxSlugVolume: number,
 *   meanFrequency: number}} slug arrival/departure times [s], lengths [m] and liquid volumes [m³]; sampled
 *   history t [s], qL, qG [m³/s] (interval averages); time-mean and peak rates [m³/s]; largest slug liquid
 *   volume [m³]; observed frequency [1/s]
 */
export function slugTrain({ D, vm, vt, slugHoldup, filmHoldup, vLF, vGF, slugFraction, meanSlugLength: Lmean, sigmaLn = 0.5,
  duration, nSamples = 1000, seed = 1 }) {
  checkArg('D', D, { min: 0, exclusiveMin: true });
  checkArg('vt', vt, { min: 0, exclusiveMin: true });
  checkArg('slugFraction', slugFraction, { min: 0, max: 1, exclusiveMin: true });
  checkArg('duration', duration, { min: 0, exclusiveMin: true });
  checkArg('nSamples', nSamples, { min: 2, max: 200000 });
  const A = area(D);
  const d = slugLengthDistribution({ meanLength: Lmean, sigmaLn });
  const rng = makeRng(seed);
  const qLs = A * vm * slugHoldup;
  const qGs = A * vm * (1 - slugHoldup);
  const qLf = A * vLF * filmHoldup;
  const qGf = A * vGF * (1 - filmHoldup);
  const slugs = [];
  let t = 0;
  let slugTime = 0;
  while (t < duration) {
    const L = d.sample(rng);
    t += L * (1 - slugFraction) / slugFraction / vt;
    if (t >= duration) break;
    const tTail = Math.min(duration, t + L / vt);
    slugs.push({ tFront: t, tTail, length: L, volume: A * L * slugHoldup });
    slugTime += tTail - t;
    t = tTail;
  }
  const dt = duration / nSamples;
  const ts = new Float64Array(nSamples);
  const qL = new Float64Array(nSamples);
  const qG = new Float64Array(nSamples);
  let k = 0;
  for (let i = 0; i < nSamples; i++) {
    const t0 = i * dt;
    const t1 = t0 + dt;
    ts[i] = t0 + dt / 2;
    while (k < slugs.length && slugs[k].tTail <= t0) k++;
    let inSlug = 0;
    for (let j = k; j < slugs.length && slugs[j].tFront < t1; j++) {
      inSlug += Math.max(0, Math.min(t1, slugs[j].tTail) - Math.max(t0, slugs[j].tFront));
    }
    const w = inSlug / dt;
    qL[i] = w * qLs + (1 - w) * qLf;
    qG[i] = w * qGs + (1 - w) * qGf;
  }
  const ws = slugTime / duration;
  let maxVol = 0;
  for (const s of slugs) if (s.volume > maxVol) maxVol = s.volume;
  return {
    slugs, t: ts, qL, qG, meanQL: ws * qLs + (1 - ws) * qLf, meanQG: ws * qGs + (1 - ws) * qGf, peakQL: Math.max(qLs, qLf),
    maxSlugVolume: maxVol, meanFrequency: slugs.length / duration,
  };
}

/**
 * Kelvin–Helmholtz stability of stratified flow from the one-dimensional two-fluid model.
 * With rk = ρk/αk and Gt = (ρL − ρG) g cosθ · A/(dAL/dhL), the dynamic (characteristic) wave speeds are the
 * roots of rL (c − vL)² + rG (c − vG)² = Gt.
 *   inviscid (IKH): the roots are real, i.e. (vG − vL)² ≤ (αL/ρL + αG/ρG) Gt;
 *   viscous (VKH):  the kinematic wave speed CV of the steady momentum balance lies between them, i.e.
 *                   rL (CV − vL)² + rG (CV − vG)² ≤ Gt, with CV = −(∂F/∂HL)/(∂F/∂vsl − ∂F/∂vsg) evaluated by
 *                   central differences on the shear balance F (stratifiedForces residual) at constant vm.
 * VKH unstable marks the onset of slug/roll-wave growth; IKH unstable marks ill-posedness of the two-fluid
 * model. Reference: Barnea & Taitel (1993), Int. J. Multiphase Flow 19(4), 639–649; Lin & Hanratty (1986),
 * Int. J. Multiphase Flow 12(1), 79–98; Whitham (1974), Linear and Nonlinear Waves, §10 (wave hierarchy).
 * @param {object} o
 * @param {number} o.h liquid level hL/D [-]
 * @param {number} o.vsl superficial liquid velocity [m/s]
 * @param {number} o.vsg superficial gas velocity [m/s]
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.theta=0] inclination [rad], |θ| < π/2
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.rhoG gas density [kg/m³]
 * @param {number} o.muL liquid viscosity [Pa·s]
 * @param {number} o.muG gas viscosity [Pa·s]
 * @param {number} [o.roughness=0] wall roughness [m]
 * @param {string} [o.interfacial='smooth'] interfacial friction option (see stratifiedForces)
 * @param {number} [o.fi] constant interfacial Fanning factor [-]
 * @returns {{ikhStable: boolean, vkhStable: boolean, relativeVelocity: number, ikhLimit: number,
 *   kinematicVelocity: number, cMinus: number, cPlus: number, ikhRatio: number, vkhRatio: number}}
 *   relativeVelocity = vG − vL and ikhLimit (critical |vG − vL|) [m/s]; kinematic and dynamic wave speeds
 *   [m/s] (NaN when IKH unstable); ratios = (destabilising term)/(Gt) [-], unstable above 1
 */
export function kelvinHelmholtzStability({ h, vsl, vsg, D, theta = 0, rhoL, rhoG, muL, muG, roughness = 0,
  interfacial = 'smooth', fi }) {
  checkArg('h', h, { min: 0, max: 1, exclusiveMin: true });
  checkArg('theta', theta, { min: -Math.PI / 2, max: Math.PI / 2 });
  const fluid = { D, theta, rhoL, rhoG, muL, muG, roughness, interfacial, fi };
  const F = (hh, sl, sg) => {
    const H = stratifiedGeometry(hh).holdup;
    return { H, r: stratifiedForces({ ...fluid, h: hh, vL: sl / H, vG: sg / (1 - H) }).residual };
  };
  const g = stratifiedGeometry(h);
  const aL = g.holdup;
  const aG = 1 - aL;
  const vL = vsl / aL;
  const vG = vsg / aG;
  const Gt = (rhoL - rhoG) * G * Math.cos(theta) * (Math.PI * D / 4) / g.dAL;
  const rL = rhoL / aL;
  const rG = rhoG / aG;
  const dh = Math.min(1e-5, 0.1 * h, 0.1 * (1 - h));
  const fp = F(h + dh, vsl, vsg);
  const fm = F(h - dh, vsl, vsg);
  const dFdH = (fp.r - fm.r) / (fp.H - fm.H);
  const dv = 1e-6 * Math.max(Math.abs(vsl), Math.abs(vsg), 1e-3);
  const dFdSl = (F(h, vsl + dv, vsg).r - F(h, vsl - dv, vsg).r) / (2 * dv);
  const dFdSg = (F(h, vsl, vsg + dv).r - F(h, vsl, vsg - dv).r) / (2 * dv);
  const CV = -dFdH / (dFdSl - dFdSg);
  const ikhLimit = Math.sqrt((aL / rhoL + aG / rhoG) * Gt);
  const rel = vG - vL;
  const ikhRatio = rel * rel / (ikhLimit * ikhLimit);
  const vkhRatio = (rL * (CV - vL) ** 2 + rG * (CV - vG) ** 2) / Gt;
  const disc = (rL + rG) * Gt - rL * rG * rel * rel;
  const cMid = (rL * vL + rG * vG) / (rL + rG);
  const half = disc >= 0 ? Math.sqrt(disc) / (rL + rG) : NaN;
  return {
    ikhStable: ikhRatio <= 1, vkhStable: vkhRatio <= 1, relativeVelocity: rel, ikhLimit, kinematicVelocity: CV,
    cMinus: cMid - half, cPlus: cMid + half, ikhRatio, vkhRatio,
  };
}

/**
 * Minimum-slip (slug-growth) criterion: a slug of holdup HLS travelling at vm over a stratified layer
 * (holdup HLF, liquid velocity vLF) has a front velocity vF = (HLS vm − HLF vLF)/(HLS − HLF) from the liquid
 * balance across the front; it survives when the front outruns its tail, vF ≥ vt = C0 vm + vd.
 * Reference: Bendiksen & Espedal (1992), Int. J. Multiphase Flow 18(2), 237–247.
 * @param {object} o
 * @param {number} o.vm mixture velocity [m/s]
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.theta=0] inclination [rad]
 * @param {number} o.filmHoldup stratified-layer holdup ahead of the slug HLF [-]
 * @param {number} o.vLF stratified-layer liquid velocity [m/s]
 * @param {number} [o.slugHoldup] HLS [-] (default Gregory et al.)
 * @param {number} [o.C0] override [-]
 * @param {number} [o.vd] override [m/s]
 * @returns {{frontVelocity: number, tailVelocity: number, growthRate: number, slugSurvives: boolean}}
 *   velocities and growthRate = vF − vt (rate of change of slug length) [m/s]
 */
export function slugGrowthCriterion({ vm, D, theta = 0, filmHoldup, vLF, slugHoldup, C0, vd }) {
  checkArg('filmHoldup', filmHoldup, { min: 0, max: 1 });
  const HLS = slugHoldup !== undefined ? checkArg('slugHoldup', slugHoldup, { min: 0, max: 1, exclusiveMin: true })
    : slugBodyHoldupGregory(Math.abs(vm));
  if (!(HLS > filmHoldup)) throw new RangeError(`filmHoldup = ${filmHoldup} must be below the slug holdup ${HLS}`);
  const vt = slugTranslationalVelocity({ vm, D, theta, C0, vd }).vt;
  const vF = (HLS * vm - filmHoldup * vLF) / (HLS - filmHoldup);
  return { frontVelocity: vF, tailVelocity: vt, growthRate: vF - vt, slugSurvives: vF >= vt };
}

/**
 * Bøe criterion for riser-induced (severe) slugging. With the riser base blocked, the hydrostatic head in the
 * riser grows at ρL g vsl sinφ while the gas trapped in the flowline (volume α L A, isothermal) is compressed
 * at P vsg/(α L). Severe slugging is possible when the head grows faster:
 *   vsl ≥ P / (ρL g α L sinφ) · vsg.
 * Reference: Bøe (1981), Severe slugging characteristics, Selected Topics in Two-Phase Flow, NTH Trondheim.
 * Necessary, not sufficient: also requires stratified flow in the flowline; see taitelStability.
 * @param {object} o
 * @param {number} o.vsl superficial liquid velocity in the flowline [m/s]
 * @param {number} o.vsg superficial gas velocity at flowline pressure [m/s]
 * @param {number} o.P flowline (riser-base) absolute pressure [Pa]
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.alphaPipe gas area fraction in the flowline [-], 0 < α ≤ 1
 * @param {number} o.L flowline length upstream of the riser base [m]
 * @param {number} [o.riserAngle=π/2] riser inclination from horizontal φ [rad]
 * @returns {{vslCritical: number, piSS: number, severeSlugging: boolean}} vslCritical [m/s];
 *   piSS = vslCritical/vsl [-] (severe slugging when < 1)
 */
export function boeCriterion({ vsl, vsg, P, rhoL, alphaPipe, L, riserAngle = Math.PI / 2 }) {
  checkArg('vsl', vsl, { min: 0, exclusiveMin: true });
  checkArg('vsg', vsg, { min: 0 });
  checkArg('P', P, { min: 0, exclusiveMin: true });
  checkArg('rhoL', rhoL, { min: 0, exclusiveMin: true });
  checkArg('alphaPipe', alphaPipe, { min: 0, max: 1, exclusiveMin: true });
  checkArg('L', L, { min: 0, exclusiveMin: true });
  checkArg('riserAngle', riserAngle, { min: 0, max: Math.PI / 2, exclusiveMin: true });
  const vslCritical = P * vsg / (rhoL * G * alphaPipe * L * Math.sin(riserAngle));
  return { vslCritical, piSS: vslCritical / vsl, severeSlugging: vsl >= vslCritical };
}

/**
 * Pots severe-slugging number ΠSS = (Z R T / M) WG / (g αF LF WL): ratio of the rate of gas-pressure build-up
 * in the flowline to the rate of hydrostatic-head build-up in the riser (identical to the Bøe ratio for an
 * ideal-gas law with compressibility Z). Severe slugging for ΠSS < 1.
 * Reference: Pots, Bromilow & Konijn (1987), SPE Prod. Eng. 2(4), 319–324.
 * @param {object} o
 * @param {number} o.massRateG gas mass rate WG [kg/s]
 * @param {number} o.massRateL liquid mass rate WL [kg/s], > 0
 * @param {number} o.T gas temperature [K]
 * @param {number} o.M gas molar mass [kg/mol]
 * @param {number} [o.Z=1] gas compressibility factor [-]
 * @param {number} o.alphaPipe gas area fraction in the flowline αF [-]
 * @param {number} o.L flowline length LF [m]
 * @returns {{piSS: number, severeSlugging: boolean}} piSS [-]
 */
export function potsNumber({ massRateG, massRateL, T, M, Z = 1, alphaPipe, L }) {
  checkArg('massRateG', massRateG, { min: 0 });
  checkArg('massRateL', massRateL, { min: 0, exclusiveMin: true });
  checkArg('T', T, { min: 0, exclusiveMin: true });
  checkArg('M', M, { min: 0, exclusiveMin: true });
  checkArg('Z', Z, { min: 0, exclusiveMin: true });
  checkArg('alphaPipe', alphaPipe, { min: 0, max: 1, exclusiveMin: true });
  checkArg('L', L, { min: 0, exclusiveMin: true });
  const piSS = Z * R_GAS * T / M * massRateG / (G * alphaPipe * L * massRateL);
  return { piSS, severeSlugging: piSS < 1 };
}

/**
 * Taitel stability criterion for the flowline–riser system. A gas bubble penetrating a liquid-filled riser
 * (height h) by y lightens the column by ρL g α' y while the flowline gas (length-equivalent α L) expands and
 * loses pressure (Ps + ρL g h) α' y/(α L). The steady state is stable when the net force decreases with y:
 *   Ps/P0 > [ (α/α') L − h ] / (P0/(ρL g)).
 * Reference: Taitel (1986), Int. J. Multiphase Flow 12(2), 203–217, Eq. (6).
 * @param {object} o
 * @param {number} o.Psep separator (riser-top) absolute pressure Ps [Pa]
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.alphaPipe gas area fraction in the flowline α [-]
 * @param {number} o.L flowline length [m]
 * @param {number} o.riserHeight riser height h [m]
 * @param {number} [o.alphaRiser=0.89] gas fraction of the bubble penetrating the riser α' [-] (Taitel used
 *   0.89 for a Taylor bubble in a vertical riser; typ. 0.8–0.9)
 * @param {number} [o.P0=101325] reference (atmospheric) pressure [Pa]
 * @returns {{stable: boolean, lhs: number, rhs: number, PsepMin: number}} lhs = Ps/P0, rhs [-]; PsepMin =
 *   lowest stable separator pressure [Pa] (may be negative: always stable)
 */
export function taitelStability({ Psep, rhoL, alphaPipe, L, riserHeight, alphaRiser = 0.89, P0 = 101325 }) {
  checkArg('Psep', Psep, { min: 0, exclusiveMin: true });
  checkArg('rhoL', rhoL, { min: 0, exclusiveMin: true });
  checkArg('alphaPipe', alphaPipe, { min: 0, max: 1 });
  checkArg('L', L, { min: 0, exclusiveMin: true });
  checkArg('riserHeight', riserHeight, { min: 0, exclusiveMin: true });
  checkArg('alphaRiser', alphaRiser, { min: 0, max: 1, exclusiveMin: true });
  checkArg('P0', P0, { min: 0, exclusiveMin: true });
  const PsepMin = rhoL * G * (alphaPipe / alphaRiser * L - riserHeight);
  return { stable: Psep > PsepMin, lhs: Psep / P0, rhs: PsepMin / P0, PsepMin };
}

/**
 * Severe-slugging cycle estimate from quasi-static volume and gas-mass balances (isothermal gas, liquid of
 * constant density, flowline gas fraction α constant, riser-base blockage for the whole build-up).
 *   build-up     the riser (height h, residual liquid fraction r after blow-down) fills and liquid backs x into
 *                the flowline while the gas is compressed from P0 = Ps + ρL g h r to P1 = Ps + ρL g h:
 *                  liquid  vsl t = h (1 − r)/sinφ + α x,     gas  P1 α (L − x) = P0 α L + (WG RT/A) t
 *                → t_build = [ (P1 − P0) α L + P1 h (1 − r)/sinφ ] / (WG RT/A + P1 vsl)
 *   production   the slug tail travels x at vsg1/α (vsg1 = WG RT/(P1 A)) while liquid leaves at vsl + vsg1
 *   blow-out     quasi-steady balance of the unbalanced head ρL g (h − y) against wall friction on the remaining
 *                liquid length y, evaluated at y = blowoutFraction·h:  v = √(2 g D (h − y) sinφ/(f y)) — an
 *                order-of-magnitude peak velocity; refine with the transient solver.
 * A negative x means the gas pressure reaches P1 before the riser is full (no full blockage cycle, cf. Bøe);
 * the result is then flagged `severeSlugging: false` and x is reported as 0.
 * Reference: derived from the cycle description of Schmidt, Brill & Beggs (1980), SPE J. 20(5), 407–414 and
 * Taitel (1986); the balances are first-principles, the blow-out velocity is an estimate.
 * @param {object} o
 * @param {number} o.vsl superficial liquid velocity [m/s]
 * @param {number} o.massRateG gas mass rate [kg/s]
 * @param {number} o.gasRT specific gas constant term Z R T/M [J/kg] (= P/ρG)
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} o.L flowline length [m]
 * @param {number} o.riserHeight riser height h [m]
 * @param {number} [o.riserAngle=π/2] riser inclination φ [rad]
 * @param {number} o.Psep separator pressure [Pa]
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.alphaPipe flowline gas fraction α [-]
 * @param {number} [o.residualHoldup=0] liquid fraction left in the riser after blow-down r [-]
 * @param {number} [o.frictionFactor=0.02] Darcy friction factor for the blow-out estimate [-] (0.012–0.03)
 * @param {number} [o.blowoutFraction=0.5] remaining liquid length y/h at which the blow-out velocity is
 *   evaluated [-]
 * @returns {{severeSlugging: boolean, buildUpTime: number, productionTime: number, blowoutTime: number,
 *   period: number, penetrationLength: number, slugLength: number, slugVolume: number, peakPressure: number,
 *   productionLiquidRate: number, blowoutVelocity: number, peakLiquidRate: number, meanLiquidRate: number}}
 *   times [s]; lengths [m]; slugVolume = liquid produced per cycle [m³]; peakPressure = riser-base pressure
 *   P1 [Pa]; rates [m³/s]; blowoutVelocity [m/s]
 */
export function severeSluggingCycle({ vsl, massRateG, gasRT, D, L, riserHeight, riserAngle = Math.PI / 2, Psep, rhoL,
  alphaPipe, residualHoldup = 0, frictionFactor = 0.02, blowoutFraction = 0.5 }) {
  checkArg('vsl', vsl, { min: 0, exclusiveMin: true });
  checkArg('massRateG', massRateG, { min: 0 });
  checkArg('gasRT', gasRT, { min: 0, exclusiveMin: true });
  checkArg('D', D, { min: 0, exclusiveMin: true });
  checkArg('L', L, { min: 0, exclusiveMin: true });
  checkArg('riserHeight', riserHeight, { min: 0, exclusiveMin: true });
  checkArg('riserAngle', riserAngle, { min: 0, max: Math.PI / 2, exclusiveMin: true });
  checkArg('Psep', Psep, { min: 0, exclusiveMin: true });
  checkArg('rhoL', rhoL, { min: 0, exclusiveMin: true });
  checkArg('alphaPipe', alphaPipe, { min: 0, max: 1, exclusiveMin: true });
  checkArg('residualHoldup', residualHoldup, { min: 0, max: 1 });
  checkArg('frictionFactor', frictionFactor, { min: 0, exclusiveMin: true });
  checkArg('blowoutFraction', blowoutFraction, { min: 0, max: 1, exclusiveMin: true });
  const A = area(D);
  const sinP = Math.sin(riserAngle);
  const riserLength = riserHeight / sinP;
  const P0 = Psep + rhoL * G * riserHeight * residualHoldup;
  const P1 = Psep + rhoL * G * riserHeight;
  const gasTerm = massRateG * gasRT / A;
  const fill = riserLength * (1 - residualHoldup);
  const buildUpTime = ((P1 - P0) * alphaPipe * L + P1 * fill) / (gasTerm + P1 * vsl);
  let x = (vsl * buildUpTime - fill) / alphaPipe;
  const severeSlugging = x >= 0;
  if (x < 0) x = 0;
  const vsg1 = gasTerm / P1;
  const productionTime = vsg1 > 0 ? alphaPipe * x / vsg1 : Infinity;
  const y = blowoutFraction * riserLength;
  const blowoutVelocity = Math.sqrt(2 * G * D * (riserLength - y) * sinP / (frictionFactor * y));
  const blowoutTime = riserLength / blowoutVelocity;
  const slugVolume = A * (fill + alphaPipe * x);
  const period = buildUpTime + productionTime + blowoutTime;
  const productionLiquidRate = A * (vsl + vsg1);
  return {
    severeSlugging, buildUpTime, productionTime, blowoutTime, period, penetrationLength: x,
    slugLength: riserLength + x, slugVolume, peakPressure: P1, productionLiquidRate, blowoutVelocity,
    peakLiquidRate: Math.max(productionLiquidRate, A * blowoutVelocity), meanLiquidRate: A * vsl,
  };
}

/**
 * Effect of topside choking and riser-base gas lift on the severe-slugging criteria — a simple, documented
 * stability shift, not a dynamic model:
 *   choking   raises the riser-top pressure seen by the system by the choke pressure drop: Ps → Ps + ΔPchoke in
 *             the Taitel criterion and in the flowline pressure of the Bøe criterion (Schmidt et al. 1980;
 *             Taitel 1986 discuss choking as back-pressure);
 *   gas lift  aerates the riser: with the lift gas alone the riser void is αr = vsgLift/(C0 (vsl + vsgLift) + vd)
 *             (drift flux, vd = 0.35 √(gD) vertical), so the head that can build up, and its growth rate, are
 *             multiplied by (1 − αr): ρL → ρL (1 − αr) in both criteria.
 * @param {object} o
 * @param {number} o.vsl superficial liquid velocity [m/s]
 * @param {number} o.vsg superficial gas velocity in the flowline at flowline pressure without mitigation [m/s]
 * @param {number} o.Psep separator pressure [Pa]
 * @param {number} o.rhoL liquid density [kg/m³]
 * @param {number} o.alphaPipe flowline gas fraction [-]
 * @param {number} o.L flowline length [m]
 * @param {number} o.riserHeight riser height [m]
 * @param {number} o.D riser inner diameter [m]
 * @param {number} [o.chokeDp=0] choke pressure drop [Pa]
 * @param {number} [o.vsgLift=0] lift-gas superficial velocity in the riser [m/s]
 * @param {number} [o.C0=1.2] riser distribution coefficient [-]
 * @param {number} [o.alphaRiser=0.89] Taitel's α' [-]
 * @returns {{base: {boe: object, taitel: object}, mitigated: {boe: object, taitel: object},
 *   riserVoidLift: number, chokeDpForStability: number}} criteria before/after (see boeCriterion,
 *   taitelStability); riserVoidLift αr [-]; chokeDpForStability = choke ΔP that just satisfies Taitel's
 *   criterion with the stated gas lift [Pa] (0 when already stable)
 */
export function severeSluggingMitigation({ vsl, vsg, Psep, rhoL, alphaPipe, L, riserHeight, D, chokeDp = 0, vsgLift = 0,
  C0 = 1.2, alphaRiser = 0.89 }) {
  checkArg('chokeDp', chokeDp, { min: 0 });
  checkArg('vsgLift', vsgLift, { min: 0 });
  checkArg('D', D, { min: 0, exclusiveMin: true });
  const evalCase = (ps, rho) => {
    const P = ps + rho * G * riserHeight;
    const Pbase = Psep + rhoL * G * riserHeight;
    return {
      // vsg is rescaled isothermally from the unmitigated flowline pressure to the mitigated one.
      boe: boeCriterion({ vsl, vsg: vsg * Pbase / P, P, rhoL: rho, alphaPipe, L }),
      taitel: taitelStability({ Psep: ps, rhoL: rho, alphaPipe, L, riserHeight, alphaRiser }),
    };
  };
  const riserVoidLift = vsgLift > 0 ? Math.min(1, vsgLift / (C0 * (vsl + vsgLift) + 0.35 * Math.sqrt(G * D))) : 0;
  const rhoEff = Math.max(rhoL * (1 - riserVoidLift), 1e-6);
  const base = evalCase(Psep, rhoL);
  const mitigated = evalCase(Psep + chokeDp, rhoEff);
  const need = taitelStability({ Psep, rhoL: rhoEff, alphaPipe, L, riserHeight, alphaRiser }).PsepMin - Psep;
  return { base, mitigated, riserVoidLift, chokeDpForStability: Math.max(0, need) };
}

/**
 * Surge volume a slug catcher or separator must absorb: liquid arrives at qIn for the slug duration while the
 * vessel is drained at qDrain, V = max(0, Vslug − qDrain · tSlug) = Vslug (1 − qDrain/qIn).
 * Reference: liquid volume balance (e.g. GPSA Engineering Data Book, §17 slug-catcher sizing practice).
 * @param {object} o
 * @param {number} o.slugVolume liquid volume of the design slug [m³]
 * @param {number} o.inflowRate liquid arrival rate during the slug qIn [m³/s] (A vm HLS for a hydrodynamic slug)
 * @param {number} o.drainRate liquid drain (processing) rate [m³/s]
 * @param {number} [o.margin=1] design margin applied to the surge volume [-]
 * @returns {{surgeVolume: number, slugDuration: number, drainedDuringSlug: number}} volumes [m³]; duration [s]
 */
export function surgeVolume({ slugVolume, inflowRate, drainRate, margin = 1 }) {
  checkArg('slugVolume', slugVolume, { min: 0 });
  checkArg('inflowRate', inflowRate, { min: 0, exclusiveMin: true });
  checkArg('drainRate', drainRate, { min: 0 });
  checkArg('margin', margin, { min: 1 });
  const slugDuration = slugVolume / inflowRate;
  const drained = Math.min(slugVolume, drainRate * slugDuration);
  return { surgeVolume: margin * (slugVolume - drained), slugDuration, drainedDuringSlug: drained };
}

/**
 * Pigging slug volume from the steady holdup: a pig travelling at the mixture velocity sweeps the liquid
 * inventory A ∫ HL ds ahead of it, while liquid keeps leaving the outlet at vsl during the transit; the liquid
 * delivered as a slug is the difference, V = A ∫ (HL − λL) ds with λL = vsl/vm the no-slip holdup.
 * Reference: volume balance; cf. Barua (1982) / McDonald & Baker (1964) sphering model.
 * @param {object} o
 * @param {number} o.D pipe inner diameter [m]
 * @param {number|number[]} o.holdup steady liquid holdup, uniform or per segment [-]
 * @param {number|number[]} o.lambdaL no-slip holdup, uniform or per segment [-]
 * @param {number|number[]} o.length line length, or the length of each segment [m]
 * @param {number} [o.vm] mixture velocity for the arrival duration [m/s]
 * @param {number} [o.slugHoldup=1] holdup of the pig-generated slug [-]
 * @returns {{slugVolume: number, inventory: number, slugLength: number, slugDuration: number}} volumes [m³];
 *   slugLength = V/(A·slugHoldup) [m]; slugDuration = slugLength/vm [s] (NaN without vm)
 */
export function pigSlugVolume({ D, holdup, lambdaL, length, vm, slugHoldup = 1 }) {
  checkArg('D', D, { min: 0, exclusiveMin: true });
  checkArg('slugHoldup', slugHoldup, { min: 0, max: 1, exclusiveMin: true });
  const A = area(D);
  const Ls = Array.isArray(length) ? length : [length];
  const at = (v, i) => (Array.isArray(v) ? v[i] : v);
  let inv = 0;
  let vol = 0;
  for (let i = 0; i < Ls.length; i++) {
    const H = checkArg('holdup', at(holdup, i), { min: 0, max: 1 });
    const lam = checkArg('lambdaL', at(lambdaL, i), { min: 0, max: 1 });
    checkArg('length', Ls[i], { min: 0 });
    inv += A * H * Ls[i];
    vol += A * Math.max(0, H - lam) * Ls[i];
  }
  const slugLength = vol / (A * slugHoldup);
  return { slugVolume: vol, inventory: inv, slugLength, slugDuration: vm ? slugLength / vm : NaN };
}
