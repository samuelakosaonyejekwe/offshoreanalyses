/**
 * Constitutive soil models at a stress point: laboratory-test simulation and calibration.
 *
 * SI units: stresses and moduli Pa, strains as fractions (not %), angles rad. Compression positive.
 * Triaxial invariants: p′ = (σ′a + 2σ′r)/3, q = σa − σr, εv = εa + 2εr, εq = 2(εa − εr)/3.
 *
 * References
 *   Roscoe K.H. & Burland J.B. (1968) On the generalised stress–strain behaviour of 'wet' clay;
 *   Wood D.M. (1990) Soil Behaviour and Critical State Soil Mechanics, CUP (Modified Cam-Clay).
 *   Drucker D.C. & Prager W. (1952) Soil mechanics and plastic analysis or limit design, Q. Appl. Math. 10.
 *   Simo J.C. & Ortiz M. (1985) cutting-plane return algorithm; Sloan S.W. et al. (2001) Eng. Comput. 18.
 *   Duncan J.M. & Chang C.Y. (1970) Nonlinear analysis of stress and strain in soils, JSMFD 96(SM5).
 *   Hardin B.O. & Drnevich V.P. (1972) Shear modulus and damping in soils, JSMFD 98(SM7).
 *   Darendeli M.B. (2001) Development of a new family of normalized modulus reduction and material damping
 *     curves, PhD thesis, Univ. Texas at Austin.
 *   Ishibashi I. & Zhang X. (1993) Unified dynamic shear moduli and damping ratios of sand and clay,
 *     Soils Found. 33(1).
 *   Idriss I.M., Dobry R. & Singh R.D. (1978) Nonlinear behavior of soft clays during cyclic loading, JGED 104.
 *   Janbu N. (1963); Hardin B.O. (1978) The nature of stress–strain behavior for soils.
 */

const PI = Math.PI;
const PATM = 101325;

function need(name, v, min = -Infinity, max = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < min || v > max) {
    throw new RangeError(`${name} must be a number in [${min}, ${max}] (got ${v})`);
  }
  return v;
}

/* ------------------------------------------------------------------ generic p–q integrator */

// One strain sub-increment: elastic predictor followed by a cutting-plane return to the yield surface.
function subStep(model, s, dev, deq) {
  const { K, G } = model.elastic(s);
  s.p += K * dev;
  s.q += 3 * G * deq;
  if (model.afterElastic) model.afterElastic(s, dev);
  let f = model.f(s);
  if (f <= model.tol(s)) { if (model.cutoff) model.cutoff(s); return; }
  for (let it = 0; it < 50; it++) {
    const d = model.derivs(s);
    const den = d.H + d.np * K * d.mp + d.nq * 3 * G * d.mq;
    if (!(den > 0)) throw new Error('constitutive integration failed: non-positive plastic modulus (snap-back)');
    const dl = f / den;
    s.p -= dl * K * d.mp;
    s.q -= dl * 3 * G * d.mq;
    s.epv += dl * d.mp;
    s.epq += dl * d.mq;
    if (model.harden) model.harden(s, dl * d.mp, dl * d.mq);
    f = model.f(s);
    if (Math.abs(f) <= model.tol(s)) break;
  }
  if (model.cutoff) model.cutoff(s);
}

function clone(s) { return { ...s }; }

function runTriaxial(model, s0, { drainage, epsAMax, nSteps, subSteps }, record) {
  need('epsAMax', epsAMax, 1e-9, 1);
  if (drainage !== 'drained' && drainage !== 'undrained') throw new RangeError(`drainage must be 'drained' or 'undrained' (got ${drainage})`);
  const s = clone(s0);
  const out = { epsA: [0], epsV: [0], epsQ: [0], p: [s.p], q: [s.q], u: [0] };
  const p0 = s.p;
  const q0 = s.q;
  const dea = epsAMax / nSteps;
  let epsV = 0;
  let epsQ = 0;
  let devGuess = 0;
  const advance = (state, dev, deq) => {
    for (let k = 0; k < subSteps; k++) subStep(model, state, dev / subSteps, deq / subSteps);
  };
  for (let n = 1; n <= nSteps; n++) {
    if (drainage === 'undrained') {
      // εr = −εa/2  →  dεq = dεa, dεv = 0
      advance(s, 0, dea);
      epsQ += dea;
    } else {
      // constant radial stress: find dεv such that q − q0 = 3 (p − p0); dεq = dεa − dεv/3
      const resid = (dev) => {
        const t = clone(s);
        advance(t, dev, dea - dev / 3);
        return { r: t.q - q0 - 3 * (t.p - p0), t };
      };
      let a = devGuess;
      let ra = resid(a);
      let b = a + (Math.abs(a) > 0 ? 0.1 * Math.abs(a) : 1e-3 * dea);
      let rb = resid(b);
      const scale = Math.max(Math.abs(p0), 1);
      for (let it = 0; it < 60 && Math.abs(rb.r) > 1e-10 * scale; it++) {
        const denom = rb.r - ra.r;
        if (denom === 0) break;
        const c = b - (rb.r * (b - a)) / denom;
        a = b; ra = rb;
        b = c; rb = resid(b);
      }
      if (!(Math.abs(rb.r) <= 1e-6 * scale)) throw new Error('drained triaxial step did not converge');
      Object.assign(s, rb.t);
      devGuess = b;
      epsV += b;
      epsQ += dea - b / 3;
    }
    out.epsA.push(n * dea);
    out.epsV.push(epsV);
    out.epsQ.push(epsQ);
    out.p.push(s.p);
    out.q.push(s.q);
    // excess pore pressure for a conventional compression test (Δσr = 0): Δp_total = Δq/3
    out.u.push(drainage === 'undrained' ? p0 + (s.q - q0) / 3 - s.p : 0);
    if (record) record(s);
  }
  out.final = s;
  return out;
}

/* ------------------------------------------------------------------ Modified Cam-Clay */

/**
 * Modified Cam-Clay triaxial compression test from an isotropic state (single element, strain controlled).
 * Model: yield f = q² + M² p′(p′ − p′c), associated flow, hardening dp′c = p′c v/(λ − κ) dεv^p, elastic
 * K = v p′/κ and G from Poisson's ratio (or constant G). Integration: explicit sub-stepping with an elastic
 * predictor and cutting-plane correction back to the yield surface.
 * Analytical checks (Wood 1990): undrained strength from an isotropic state with over-consolidation ratio
 * R = p′c0/p′0:  p′f = p′0 (R/2)^Λ, qf = M p′f, su = qf/2, Λ = (λ − κ)/λ;
 * drained compression ends on the critical-state line at p′f = 3 p′0/(3 − M), v_f = Γ − λ ln p′f,
 * Γ = N − (λ − κ) ln 2.
 * Heavily over-consolidated drained states soften after peak; the strain-controlled integration follows
 * the softening branch as long as the plastic modulus stays positive.
 * @param {object} o
 * @param {number} o.M critical-state stress ratio [-]
 * @param {number} o.lambda slope of the normal compression line in v–ln p′ [-]
 * @param {number} o.kappa slope of the swelling line [-]
 * @param {number} o.N specific volume on the isotropic NCL at p′ = pRef [-]
 * @param {number} [o.pRef=1000] reference pressure for N [Pa] (1 kPa)
 * @param {number} o.p0 initial mean effective stress [Pa]
 * @param {number} [o.ocr=1] isotropic over-consolidation ratio p′c0/p′0 [-]
 * @param {number} [o.nu=0.3] Poisson's ratio [-] (ignored when G is given)
 * @param {number} [o.G] constant shear modulus [Pa]
 * @param {'drained'|'undrained'} o.drainage
 * @param {number} [o.epsAMax=0.2] final axial strain [-]
 * @param {number} [o.nSteps=400] output steps
 * @param {number} [o.subSteps=10] strain sub-increments per step
 * @returns {{epsA:number[],epsV:number[],epsQ:number[],p:number[],q:number[],u:number[],v:number[],
 *   pc:number[],su:number,pf:number,qf:number,vf:number,eta:number}}
 *   stress paths [Pa], excess pore pressure u [Pa], specific volume v [-], preconsolidation pc [Pa];
 *   su = q_final/2 [Pa]; final p′, q, v; eta = qf/pf [-]
 */
export function mccTriaxial({ M, lambda, kappa, N, pRef = 1000, p0, ocr = 1, nu = 0.3, G, drainage, epsAMax = 0.2, nSteps = 400, subSteps = 10 }) {
  need('M', M, 1e-6, 3);
  need('lambda', lambda, 1e-9);
  need('kappa', kappa, 1e-9, lambda * (1 - 1e-9));
  need('p0', p0, 1e-9);
  need('ocr', ocr, 1);
  need('nu', nu, 0, 0.499);
  const pc0 = ocr * p0;
  const v0 = N - lambda * Math.log(pc0 / pRef) + kappa * Math.log(pc0 / p0);
  if (!(v0 > 1)) throw new RangeError('N, lambda and p0 give a specific volume ≤ 1');
  const model = {
    elastic: (s) => {
      const K = (s.v * s.p) / kappa;
      return { K, G: G !== undefined ? G : (3 * K * (1 - 2 * nu)) / (2 * (1 + nu)) };
    },
    afterElastic: (s, dev) => { s.v -= s.v * dev; },
    f: (s) => s.q * s.q + M * M * s.p * (s.p - s.pc),
    tol: (s) => 1e-10 * M * M * s.pc * s.pc,
    derivs: (s) => {
      const np = M * M * (2 * s.p - s.pc);
      const nq = 2 * s.q;
      return { np, nq, mp: np, mq: nq, H: (M * M * s.p * s.pc * s.v * np) / (lambda - kappa) };
    },
    harden: (s, depv) => { s.pc *= Math.exp((s.v * depv) / (lambda - kappa)); },
  };
  const vHist = [v0];
  const pcHist = [pc0];
  const s0 = { p: p0, q: 0, pc: pc0, v: v0, epv: 0, epq: 0 };
  const res = runTriaxial(model, s0, { drainage, epsAMax, nSteps, subSteps }, (s) => { vHist.push(s.v); pcHist.push(s.pc); });
  const pf = res.p[res.p.length - 1];
  const qf = res.q[res.q.length - 1];
  return { epsA: res.epsA, epsV: res.epsV, epsQ: res.epsQ, p: res.p, q: res.q, u: res.u, v: vHist, pc: pcHist, su: qf / 2, pf, qf, vf: vHist[vHist.length - 1], eta: qf / pf };
}

/**
 * Analytical Modified Cam-Clay undrained strength ratio in triaxial compression from an isotropic state:
 *   su / p′0 = (M/2) (R/2)^Λ,   Λ = (λ − κ)/λ,   R = isotropic over-consolidation ratio.
 * @param {object} o M, lambda, kappa, ocr (default 1)
 * @returns {number} su/p′0 [-]
 */
export function mccUndrainedStrengthRatio({ M, lambda, kappa, ocr = 1 }) {
  need('M', M, 1e-6, 3);
  need('lambda', lambda, 1e-9);
  need('kappa', kappa, 0, lambda);
  need('ocr', ocr, 1);
  return (M / 2) * (ocr / 2) ** ((lambda - kappa) / lambda);
}

/* ------------------------------------------------------------------ Mohr–Coulomb / Drucker–Prager */

/**
 * Triaxial-plane parameters equivalent to Mohr–Coulomb strength: q = M p′ + kq with
 *   compression: M = 6 sin φ/(3 − sin φ), kq = 6 c cos φ/(3 − sin φ);
 *   extension:   M = 6 sin φ/(3 + sin φ), kq = 6 c cos φ/(3 + sin φ).
 * @param {object} o phi [rad], c [Pa], mode 'compression' | 'extension'
 * @returns {{M:number,kq:number}}
 */
export function mohrCoulombTriaxialParams({ phi, c = 0, mode = 'compression' }) {
  need('phi', phi, 0, PI / 2 - 1e-9);
  need('c', c, 0);
  const s = Math.sin(phi);
  const den = mode === 'extension' ? 3 + s : 3 - s;
  return { M: (6 * s) / den, kq: (6 * c * Math.cos(phi)) / den };
}

/**
 * Drucker–Prager parameters √J2 = α I1 + k (compression positive) matched to Mohr–Coulomb:
 *   'compression' cone:  α = 2 sin φ / (√3 (3 − sin φ)),  k = 6 c cos φ / (√3 (3 − sin φ));
 *   'extension' cone:    α = 2 sin φ / (√3 (3 + sin φ)),  k = 6 c cos φ / (√3 (3 + sin φ));
 *   'plane-strain' (associated flow):  α = tan φ / √(9 + 12 tan²φ),  k = 3 c / √(9 + 12 tan²φ).
 * In the triaxial plane √J2 = q/√3 and I1 = 3p′, so q = 3√3 α p′ + √3 k.
 * @param {object} o phi [rad], c [Pa], match
 * @returns {{alpha:number,k:number,M:number,kq:number}} α [-], k [Pa]; triaxial slope M [-] and intercept kq [Pa]
 */
export function druckerPragerParams({ phi, c = 0, match = 'compression' }) {
  need('phi', phi, 0, PI / 2 - 1e-9);
  need('c', c, 0);
  const s = Math.sin(phi);
  const r3 = Math.sqrt(3);
  let alpha;
  let k;
  if (match === 'compression') { alpha = (2 * s) / (r3 * (3 - s)); k = (6 * c * Math.cos(phi)) / (r3 * (3 - s)); }
  else if (match === 'extension') { alpha = (2 * s) / (r3 * (3 + s)); k = (6 * c * Math.cos(phi)) / (r3 * (3 + s)); }
  else if (match === 'plane-strain') { const t = Math.tan(phi); const d = Math.sqrt(9 + 12 * t * t); alpha = t / d; k = (3 * c) / d; }
  else throw new RangeError(`match must be 'compression', 'extension' or 'plane-strain' (got ${match})`);
  return { alpha, k, M: 3 * r3 * alpha, kq: r3 * k };
}

/**
 * Linear-elastic perfectly-plastic frictional model in triaxial compression: Mohr–Coulomb (default) or
 * Drucker–Prager, with non-associated flow through the dilation angle and an optional tension cut-off on
 * p′. Yield q = M p′ + kq; plastic potential q = Md p′ (Md from ψ in the same way as M from φ).
 * Checks: drained failure q_f = 3 (M p′0 + kq)/(3 − M), i.e. σ′1 = σ′3 tan²(45° + φ/2) + 2c tan(45° + φ/2);
 * undrained with ψ = 0: p′ stays at p′0 after yield and q_f = M p′0 + kq.
 * @param {object} o
 * @param {'mohr-coulomb'|'drucker-prager'} [o.model='mohr-coulomb']
 * @param {number} o.E Young's modulus [Pa]  @param {number} o.nu Poisson's ratio [-]
 * @param {number} o.phi friction angle [rad]  @param {number} [o.c=0] cohesion [Pa]
 * @param {number} [o.psi=0] dilation angle [rad]
 * @param {'compression'|'extension'|'plane-strain'} [o.match='compression'] Drucker–Prager matching
 * @param {number} [o.tensionCutoff] limiting tensile mean stress (positive number) [Pa]
 * @param {number} o.p0 initial isotropic effective stress [Pa]
 * @param {'drained'|'undrained'} o.drainage
 * @param {number} [o.epsAMax=0.05]  @param {number} [o.nSteps=200]  @param {number} [o.subSteps=5]
 * @returns {{epsA:number[],epsV:number[],epsQ:number[],p:number[],q:number[],u:number[],qf:number,pf:number,
 *   M:number,kq:number}}
 */
export function frictionalTriaxial({ model = 'mohr-coulomb', E, nu, phi, c = 0, psi = 0, match = 'compression', tensionCutoff, p0, drainage, epsAMax = 0.05, nSteps = 200, subSteps = 5 }) {
  need('E', E, 1e-9);
  need('nu', nu, 0, 0.499);
  need('psi', psi, 0, phi);
  need('p0', p0, 0);
  let M;
  let kq;
  let Md;
  if (model === 'mohr-coulomb') {
    ({ M, kq } = mohrCoulombTriaxialParams({ phi, c }));
    Md = mohrCoulombTriaxialParams({ phi: psi }).M;
  } else if (model === 'drucker-prager') {
    ({ M, kq } = druckerPragerParams({ phi, c, match }));
    Md = druckerPragerParams({ phi: psi, match }).M;
  } else throw new RangeError(`model must be 'mohr-coulomb' or 'drucker-prager' (got ${model})`);
  const K = E / (3 * (1 - 2 * nu));
  const G = E / (2 * (1 + nu));
  const scale = Math.max(M * p0 + kq, 1);
  const mdl = {
    elastic: () => ({ K, G }),
    f: (s) => s.q - M * s.p - kq,
    tol: () => 1e-10 * scale,
    derivs: () => ({ np: -M, nq: 1, mp: -Md, mq: 1, H: 0 }),
    cutoff: tensionCutoff === undefined ? undefined : (s) => {
      if (s.p < -tensionCutoff) { s.p = -tensionCutoff; s.q = Math.max(Math.min(s.q, M * s.p + kq), 0); }
    },
  };
  const res = runTriaxial(mdl, { p: p0, q: 0, epv: 0, epq: 0 }, { drainage, epsAMax, nSteps, subSteps });
  return { epsA: res.epsA, epsV: res.epsV, epsQ: res.epsQ, p: res.p, q: res.q, u: res.u, qf: res.q[res.q.length - 1], pf: res.p[res.p.length - 1], M, kq };
}

/**
 * Direct simple shear of an elastic perfectly-plastic Mohr–Coulomb element, assuming the horizontal plane
 * is the plane of maximum stress obliquity (τ = c + σ′n tan φ on it):
 *   drained (σ′n constant): τ = G γ up to τf; plastic dilation dεn = −tan ψ dγp (expansion negative);
 *   undrained (constant volume): dilation is suppressed elastically, dσ′n = Moed tan ψ dγp, so the
 *   strength rises along the failure line: dγp = G dγ / (G + Moed tan φ tan ψ).
 * @param {object} o
 * @param {number} o.G shear modulus [Pa]  @param {number} o.Moed constrained modulus [Pa]
 * @param {number} o.phi [rad]  @param {number} [o.c=0] [Pa]  @param {number} [o.psi=0] [rad]
 * @param {number} o.sigmaN0 initial normal effective stress [Pa]
 * @param {'drained'|'undrained'} o.drainage
 * @param {number} [o.gammaMax=0.1] final shear strain [-]  @param {number} [o.nSteps=200]
 * @returns {{gamma:number[],tau:number[],sigmaN:number[],epsN:number[],u:number[]}} u = σn0 − σ′n [Pa]
 */
export function simpleShearMohrCoulomb({ G, Moed, phi, c = 0, psi = 0, sigmaN0, drainage, gammaMax = 0.1, nSteps = 200 }) {
  need('G', G, 1e-9);
  need('Moed', Moed, 1e-9);
  need('phi', phi, 0, PI / 2 - 1e-9);
  need('sigmaN0', sigmaN0, 0);
  if (drainage !== 'drained' && drainage !== 'undrained') throw new RangeError(`drainage must be 'drained' or 'undrained' (got ${drainage})`);
  const tf = Math.tan(phi);
  const td = Math.tan(psi);
  const out = { gamma: [0], tau: [0], sigmaN: [sigmaN0], epsN: [0], u: [0] };
  let tau = 0;
  let sn = sigmaN0;
  let en = 0;
  const dg = gammaMax / nSteps;
  for (let i = 1; i <= nSteps; i++) {
    const tauLim = c + sn * tf;
    const trial = tau + G * dg;
    if (trial <= tauLim) tau = trial;
    else {
      const dgEl = Math.max((tauLim - tau) / G, 0);
      const dgRem = dg - dgEl;
      tau = tauLim;
      if (drainage === 'drained') en -= td * dgRem;
      else {
        const dgp = (G * dgRem) / (G + Moed * tf * td);
        sn += Moed * td * dgp;
        tau = c + sn * tf;
      }
    }
    out.gamma.push(i * dg);
    out.tau.push(tau);
    out.sigmaN.push(sn);
    out.epsN.push(en);
    out.u.push(sigmaN0 - sn);
  }
  return out;
}

/* ------------------------------------------------------------------ nonlinear elastic stiffness */

/**
 * Duncan & Chang (1970) hyperbolic model:
 *   Ei = K pa (σ3/pa)^n,   (σ1 − σ3)f = (2 c cos φ + 2 σ3 sin φ)/(1 − sin φ),
 *   (σ1 − σ3) = ε / [1/Ei + Rf ε/(σ1 − σ3)f],   Et = Ei [1 − Rf (σ1 − σ3)/(σ1 − σ3)f]²,
 *   unload–reload Eur = Kur pa (σ3/pa)^n.
 * @param {object} o
 * @param {number} o.K modulus number [-]  @param {number} o.n modulus exponent [-]
 * @param {number} o.Rf failure ratio [-] (0.6–0.95)
 * @param {number} o.phi [rad]  @param {number} [o.c=0] [Pa]
 * @param {number} o.sigma3 confining effective stress [Pa]
 * @param {number} [o.epsA] axial strain at which to evaluate [-]
 * @param {number} [o.q] deviator stress at which to evaluate Et [Pa] (default: from epsA)
 * @param {number} [o.Kur] unload–reload modulus number [-]
 * @param {number} [o.pa=101325] [Pa]
 * @returns {{Ei:number,qf:number,qUlt:number,q:number|undefined,Et:number|undefined,Eur:number|undefined}} [Pa]
 */
export function duncanChang({ K, n, Rf, phi, c = 0, sigma3, epsA, q, Kur, pa = PATM }) {
  need('K', K, 1e-9);
  need('n', n, 0, 2);
  need('Rf', Rf, 0.1, 1);
  need('phi', phi, 0, PI / 2 - 1e-9);
  need('sigma3', sigma3, 1e-9);
  const Ei = K * pa * (sigma3 / pa) ** n;
  const s = Math.sin(phi);
  const qf = (2 * c * Math.cos(phi) + 2 * sigma3 * s) / (1 - s);
  let qq = q;
  if (qq === undefined && epsA !== undefined) qq = need('epsA', epsA, 0) / (1 / Ei + (Rf * epsA) / qf);
  const Et = qq !== undefined ? Ei * Math.max(1 - (Rf * qq) / qf, 0) ** 2 : undefined;
  return { Ei, qf, qUlt: qf / Rf, q: qq, Et, Eur: Kur !== undefined ? Kur * pa * (sigma3 / pa) ** n : undefined };
}

/**
 * Janbu (1963) stress-dependent constrained modulus: M = m σr (σ′/σr)^(1 − j).
 * @param {object} o m modulus number [-], j exponent [-], sigma effective stress [Pa], sigmaRef [Pa] (100 kPa)
 * @returns {number} M [Pa]
 */
export function janbuModulus({ m, j, sigma, sigmaRef = 100e3 }) {
  need('m', m, 1e-9);
  need('j', j, 0, 1);
  need('sigma', sigma, 1e-9);
  return m * sigmaRef * (sigma / sigmaRef) ** (1 - j);
}

/**
 * Small-strain shear modulus (Hardin 1978): Gmax = 625 F(e) OCR^k √(pa σ′m),  F(e) = 1/(0.3 + 0.7 e²);
 * k from 0 (PI = 0) to 0.5 (PI ≥ 100) – argument.
 * @param {object} o e void ratio [-], sigmaM mean effective stress [Pa], ocr [-], k [-], pa [Pa]
 * @returns {number} Gmax [Pa]
 */
export function gmaxHardin({ e, sigmaM, ocr = 1, k = 0, pa = PATM }) {
  need('e', e, 0.05, 5);
  need('sigmaM', sigmaM, 1e-9);
  need('ocr', ocr, 1);
  return (625 / (0.3 + 0.7 * e * e)) * ocr ** k * Math.sqrt(pa * sigmaM);
}

/* ------------------------------------------------------------------ modulus reduction and damping */

/**
 * Hardin & Drnevich (1972) hyperbolic model: G/Gmax = 1/(1 + γ/γr),  D = Dmax (γ/γr)/(1 + γ/γr),
 * γr = τmax/Gmax the reference strain.
 * @param {number} gamma shear strain [-]
 * @param {{gammaRef:number,Dmax?:number,Dmin?:number}} o reference strain [-], limiting damping ratios [-]
 * @returns {{ratio:number,damping:number}} G/Gmax [-], damping ratio [-]
 */
export function hardinDrnevich(gamma, { gammaRef, Dmax = 0.25, Dmin = 0 }) {
  need('gamma', gamma, 0);
  need('gammaRef', gammaRef, 1e-12);
  const x = gamma / gammaRef;
  return { ratio: 1 / (1 + x), damping: Dmin + (Dmax * x) / (1 + x) };
}

/**
 * Darendeli (2001) modulus-reduction and damping curves (mean values):
 *   G/Gmax = 1/[1 + (γ/γr)^a],  a = 0.9190,
 *   γr [%] = (0.0352 + 0.0010 PI OCR^0.3246) σ′0^0.3483   (σ′0 in atm),
 *   Dmin [%] = (0.8005 + 0.0129 PI OCR^−0.1069) σ′0^−0.2889 [1 + 0.2919 ln f],
 *   D = b (G/Gmax)^0.1 D_Masing + Dmin,  b = 0.6329 − 0.0057 ln N,
 *   D_Masing = c1 D1 + c2 D1² + c3 D1³,  D1 [%] = (100/π) {4 [γ − γr ln((γ + γr)/γr)] / [γ²/(γ + γr)] − 2},
 *   c1 = −1.1143a² + 1.8618a + 0.2523,  c2 = 0.0805a² − 0.0710a − 0.0095,  c3 = −0.0005a² + 0.0002a + 0.0003.
 * Calibrated for strains 1e-6 to about 1e-2, σ′0 0.25–16 atm, PI 0–100 %.
 * @param {number} gamma shear strain [-]
 * @param {object} o
 * @param {number} o.sigmaM mean effective confining stress [Pa]
 * @param {number} [o.PI=0] plasticity index [%]  @param {number} [o.ocr=1] [-]
 * @param {number} [o.freq=1] loading frequency [Hz]  @param {number} [o.N=10] number of cycles [-]
 * @returns {{ratio:number,damping:number,gammaRef:number}} G/Gmax [-], damping ratio [-], γr [-]
 */
export function darendeli(gamma, { sigmaM, PI: pIdx = 0, ocr = 1, freq = 1, N = 10 }) {
  need('gamma', gamma, 0);
  need('sigmaM', sigmaM, 1e-9);
  need('PI', pIdx, 0);
  need('ocr', ocr, 1);
  const a = 0.919;
  const s0 = sigmaM / PATM;
  const gr = (0.0352 + 0.001 * pIdx * ocr ** 0.3246) * s0 ** 0.3483; // %
  const g = gamma * 100; // %
  const ratio = 1 / (1 + (g / gr) ** a);
  const dmin = (0.8005 + 0.0129 * pIdx * ocr ** -0.1069) * s0 ** -0.2889 * (1 + 0.2919 * Math.log(freq));
  let dm = 0;
  if (g > 1e-9 * gr) {
    // series form of the bracket for small γ/γr avoids cancellation
    const x = g / gr;
    const d1 = x < 1e-3 ? (100 / PI) * ((2 * x) / 3 - x * x / 3) : (100 / PI) * ((4 * (g - gr * Math.log((g + gr) / gr))) / ((g * g) / (g + gr)) - 2);
    const c1 = -1.1143 * a * a + 1.8618 * a + 0.2523;
    const c2 = 0.0805 * a * a - 0.071 * a - 0.0095;
    const c3 = -0.0005 * a * a + 0.0002 * a + 0.0003;
    dm = c1 * d1 + c2 * d1 * d1 + c3 * d1 ** 3;
  }
  const b = 0.6329 - 0.0057 * Math.log(N);
  return { ratio, damping: (b * ratio ** 0.1 * dm + dmin) / 100, gammaRef: gr / 100 };
}

/**
 * Ishibashi & Zhang (1993) unified curves:
 *   G/Gmax = K(γ, PI) σ′m^(m(γ, PI) − m0)   (σ′m in kPa),
 *   K = 0.5 {1 + tanh[ln((0.000102 + n(PI))/γ)^0.492]},
 *   m − m0 = 0.272 {1 − tanh[ln((0.000556/γ)^0.4)]} exp(−0.0145 PI^1.3),
 *   n(PI) = 0 (PI = 0); 3.37e-6 PI^1.404 (0 < PI ≤ 15); 7.0e-7 PI^1.976 (15 < PI ≤ 70); 2.7e-5 PI^1.115 (PI > 70),
 *   D = 0.333 [1 + exp(−0.0145 PI^1.3)]/2 · [0.586 (G/Gmax)² − 1.547 (G/Gmax) + 1].
 * G/Gmax is capped at 1.
 * @param {number} gamma shear strain [-], > 0
 * @param {{sigmaM:number,PI?:number}} o mean effective stress [Pa], plasticity index [%]
 * @returns {{ratio:number,damping:number}}
 */
export function ishibashiZhang(gamma, { sigmaM, PI: pIdx = 0 }) {
  need('gamma', gamma, 0);
  need('sigmaM', sigmaM, 1e-9);
  need('PI', pIdx, 0);
  const g = Math.max(gamma, 1e-9);
  let n = 0;
  if (pIdx > 70) n = 2.7e-5 * pIdx ** 1.115;
  else if (pIdx > 15) n = 7.0e-7 * pIdx ** 1.976;
  else if (pIdx > 0) n = 3.37e-6 * pIdx ** 1.404;
  const K = 0.5 * (1 + Math.tanh(Math.log(((0.000102 + n) / g) ** 0.492)));
  const mm = 0.272 * (1 - Math.tanh(Math.log((0.000556 / g) ** 0.4))) * Math.exp(-0.0145 * pIdx ** 1.3);
  const ratio = Math.min(K * (sigmaM / 1000) ** mm, 1);
  const damping = ((0.333 * (1 + Math.exp(-0.0145 * pIdx ** 1.3))) / 2) * (0.586 * ratio * ratio - 1.547 * ratio + 1);
  return { ratio, damping };
}

/**
 * Uniform access to the modulus-reduction and damping models, used by the site-response engine.
 * @param {number} gamma shear strain [-]
 * @param {{model:'linear'|'hardin-drnevich'|'darendeli'|'ishibashi-zhang'|'table', [k:string]:any}} curve
 *   model parameters as for the individual functions; 'linear' takes {damping}; 'table' takes
 *   {strain:number[], ratio:number[], damping:number[]} interpolated linearly in log10(strain)
 * @returns {{ratio:number,damping:number}}
 */
export function modulusReduction(gamma, curve) {
  if (!curve || typeof curve.model !== 'string') throw new RangeError('curve.model is required');
  switch (curve.model) {
    case 'linear': return { ratio: 1, damping: curve.damping || 0 };
    case 'hardin-drnevich': return hardinDrnevich(gamma, curve);
    case 'darendeli': return darendeli(gamma, curve);
    case 'ishibashi-zhang': return ishibashiZhang(gamma, curve);
    case 'table': {
      const { strain, ratio, damping } = curve;
      const n = strain.length;
      if (gamma <= strain[0]) return { ratio: ratio[0], damping: damping[0] };
      if (gamma >= strain[n - 1]) return { ratio: ratio[n - 1], damping: damping[n - 1] };
      let i = 1;
      while (strain[i] < gamma) i++;
      const w = (Math.log10(gamma) - Math.log10(strain[i - 1])) / (Math.log10(strain[i]) - Math.log10(strain[i - 1]));
      return { ratio: ratio[i - 1] + w * (ratio[i] - ratio[i - 1]), damping: damping[i - 1] + w * (damping[i] - damping[i - 1]) };
    }
    default: throw new RangeError(`unknown modulus-reduction model '${curve.model}'`);
  }
}

/**
 * Cyclic degradation index of clay under strain-controlled cycling (Idriss, Dobry & Singh 1978):
 *   δ = G_N / G_1 = N^(−t),
 * t the degradation parameter (depends on cyclic strain amplitude, PI and OCR; e.g. Vucetic & Dobry 1988).
 * @param {number} N number of cycles [-], ≥ 1
 * @param {number} t degradation parameter [-], ≥ 0
 * @returns {number} δ [-]
 */
export function cyclicDegradationIndex(N, t) {
  need('N', N, 1);
  need('t', t, 0);
  return N ** -t;
}
