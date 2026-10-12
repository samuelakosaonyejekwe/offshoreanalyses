/**
 * Shells of revolution: pressure stresses in cylinders and spheres, elastic buckling of cylindrical shells,
 * collapse of pipes under external pressure, ovalisation in bending and equivalent cylinders for cones.
 *
 * Sign convention: stresses are positive in tension. Internal pressure is positive `pi`, external `po`.
 *
 * References: Timoshenko & Goodier, "Theory of Elasticity" §28, §136 (Lamé); Timoshenko & Gere, "Theory of
 * Elastic Stability" ch. 7, 11; DNV-RP-C202 (Buckling Strength of Shells) §3.2–3.4, §6; DNV-ST-F101 §5.4.4
 * (system collapse); Brazier (1927) Proc. R. Soc. A 116, 104–114.
 */

function pos(v, name, fn) {
  if (!(typeof v === 'number' && v > 0 && Number.isFinite(v))) throw new RangeError(`${fn}: ${name} must be a positive number`);
}
function num(v, name, fn) {
  if (!(typeof v === 'number' && Number.isFinite(v))) throw new RangeError(`${fn}: ${name} must be a finite number`);
}
function poisson(nu, fn) {
  if (!(typeof nu === 'number' && nu >= 0 && nu < 0.5)) throw new RangeError(`${fn}: nu must lie in [0, 0.5)`);
}

/**
 * Lamé solution for a thick-walled cylinder under internal and external pressure.
 * σr = A − B/r², σθ = A + B/r² with A = (pi ri² − po ro²)/(ro² − ri²), B = (pi − po) ri² ro²/(ro² − ri²).
 * Axial stress: 'closed' ends σz = A (end-cap force carried by the wall), 'open' σz = 0, 'planeStrain'
 * σz = ν (σr + σθ).
 * @param {object} o
 * @param {number} o.ri inner radius [m]
 * @param {number} o.ro outer radius [m]
 * @param {number} [o.pi] internal pressure [Pa] (default 0)
 * @param {number} [o.po] external pressure [Pa] (default 0)
 * @param {number} [o.r] radius of evaluation [m] (default ri)
 * @param {'closed'|'open'|'planeStrain'} [o.ends] axial condition (default 'closed')
 * @param {number} [o.E] Young's modulus [Pa] (needed for the displacement)
 * @param {number} [o.nu] Poisson's ratio (default 0.3)
 * @returns {{sigmaR: number, sigmaTheta: number, sigmaZ: number, vonMises: number, tresca: number,
 *   u: number|null}} stresses at r [Pa] and radial displacement [m] (null without E)
 */
export function lameCylinder({ ri, ro, pi = 0, po = 0, r, ends = 'closed', E, nu = 0.3 }) {
  const fn = 'lameCylinder';
  pos(ri, 'ri', fn);
  if (!(ro > ri)) throw new RangeError(`${fn}: ro must exceed ri`);
  num(pi, 'pi', fn); num(po, 'po', fn); poisson(nu, fn);
  const rr = r === undefined ? ri : r;
  if (!(rr >= ri * (1 - 1e-12) && rr <= ro * (1 + 1e-12))) throw new RangeError(`${fn}: r must lie within the wall`);
  const d = ro * ro - ri * ri;
  const A = (pi * ri * ri - po * ro * ro) / d;
  const B = ((pi - po) * ri * ri * ro * ro) / d;
  const sr = A - B / (rr * rr);
  const st = A + B / (rr * rr);
  let sz;
  if (ends === 'closed') sz = A;
  else if (ends === 'open') sz = 0;
  else if (ends === 'planeStrain') sz = nu * (sr + st);
  else throw new RangeError(`${fn}: ends must be 'closed', 'open' or 'planeStrain'`);
  const vm = Math.sqrt(0.5 * ((sr - st) ** 2 + (st - sz) ** 2 + (sz - sr) ** 2));
  const tr = Math.max(sr, st, sz) - Math.min(sr, st, sz);
  let u = null;
  if (E !== undefined) { pos(E, 'E', fn); u = (rr / E) * (st - nu * (sr + sz)); }
  return { sigmaR: sr, sigmaTheta: st, sigmaZ: sz, vonMises: vm, tresca: tr, u };
}

/**
 * Thin-walled cylinder membrane stresses under net internal pressure: σθ = p r / t, σz = p r / (2 t) (closed
 * ends), with r the mean radius. Adequate for r/t ≳ 10 (error of the hoop stress at the bore ≈ t/(2r)).
 * @param {object} o
 * @param {number} o.r mean radius [m]
 * @param {number} o.t wall thickness [m]
 * @param {number} o.p net internal pressure pi − po [Pa]
 * @param {boolean} [o.closedEnds] include the end-cap axial stress (default true)
 * @returns {{sigmaTheta: number, sigmaZ: number, vonMises: number}} [Pa]
 */
export function thinCylinder({ r, t, p, closedEnds = true }) {
  const fn = 'thinCylinder';
  pos(r, 'r', fn); pos(t, 't', fn); num(p, 'p', fn);
  const st = (p * r) / t;
  const sz = closedEnds ? st / 2 : 0;
  return { sigmaTheta: st, sigmaZ: sz, vonMises: Math.sqrt(st * st - st * sz + sz * sz) };
}

/**
 * Lamé solution for a thick-walled sphere under internal and external pressure:
 * σr = (pi ri³ − po ro³)/(ro³ − ri³) − (pi − po) ri³ ro³ / ((ro³ − ri³) r³),
 * σθ = (pi ri³ − po ro³)/(ro³ − ri³) + (pi − po) ri³ ro³ / (2 (ro³ − ri³) r³).
 * The thin-wall limit is σθ = p r / (2 t).
 * @param {object} o
 * @param {number} o.ri inner radius [m]
 * @param {number} o.ro outer radius [m]
 * @param {number} [o.pi] internal pressure [Pa]
 * @param {number} [o.po] external pressure [Pa]
 * @param {number} [o.r] radius of evaluation [m] (default ri)
 * @returns {{sigmaR: number, sigmaTheta: number, vonMises: number}} [Pa]
 */
export function lameSphere({ ri, ro, pi = 0, po = 0, r }) {
  const fn = 'lameSphere';
  pos(ri, 'ri', fn);
  if (!(ro > ri)) throw new RangeError(`${fn}: ro must exceed ri`);
  const rr = r === undefined ? ri : r;
  if (!(rr >= ri * (1 - 1e-12) && rr <= ro * (1 + 1e-12))) throw new RangeError(`${fn}: r must lie within the wall`);
  const d = ro ** 3 - ri ** 3;
  const A = (pi * ri ** 3 - po * ro ** 3) / d;
  const B = ((pi - po) * ri ** 3 * ro ** 3) / d;
  const sr = A - B / rr ** 3;
  const st = A + B / (2 * rr ** 3);
  return { sigmaR: sr, sigmaTheta: st, vonMises: Math.abs(st - sr) };
}

/**
 * Classical (perfect-shell) elastic buckling loads, to be combined with a knock-down factor chosen by the user:
 *  - cylinder, axial compression: σ_cl = E t / (r √(3 (1 − ν²))) (= 0.605 E t/r for ν = 0.3)
 *  - long cylinder, external lateral pressure: p_cl = E (t/r)³ / (4 (1 − ν²))
 *  - complete sphere, external pressure: p_cl = 2 E (t/r)² / √(3 (1 − ν²))
 * (Timoshenko & Gere §11.1, §7.5, §11.6.) Typical knock-down factors from tests: 0.2–0.5 for axially
 * compressed cylinders (decreasing with r/t), 0.7–0.9 for external pressure on cylinders, 0.15–0.3 for spheres.
 * @param {object} o
 * @param {number} o.r mid-surface radius [m]
 * @param {number} o.t thickness [m]
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} [o.nu] Poisson's ratio (default 0.3)
 * @param {number} [o.knockDownAxial] factor applied to the axial result (default 1)
 * @param {number} [o.knockDownPressure] factor applied to the cylinder pressure result (default 1)
 * @param {number} [o.knockDownSphere] factor applied to the sphere result (default 1)
 * @returns {{sigmaAxial: number, pCylinder: number, pSphere: number}} [Pa]
 */
export function classicalShellBuckling({ r, t, E, nu = 0.3, knockDownAxial = 1, knockDownPressure = 1, knockDownSphere = 1 }) {
  const fn = 'classicalShellBuckling';
  pos(r, 'r', fn); pos(t, 't', fn); pos(E, 'E', fn); poisson(nu, fn);
  pos(knockDownAxial, 'knockDownAxial', fn); pos(knockDownPressure, 'knockDownPressure', fn); pos(knockDownSphere, 'knockDownSphere', fn);
  const s3 = Math.sqrt(3 * (1 - nu * nu));
  return {
    sigmaAxial: (knockDownAxial * E * t) / (r * s3),
    pCylinder: (knockDownPressure * E * (t / r) ** 3) / (4 * (1 - nu * nu)),
    pSphere: (knockDownSphere * 2 * E * (t / r) ** 2) / s3,
  };
}

/**
 * Elastic buckling strengths of an unstiffened circular cylindrical shell (or the shell between ring frames)
 * per DNV-RP-C202 §3.4.2: f_E = C π² E / (12 (1 − ν²)) (t/l)², C = ψ √(1 + (ρ ξ / ψ)²), with the curvature
 * parameter Z_l = (l² / (r t)) √(1 − ν²) and Table 3-2:
 *   axial stress       ψ = 1     ξ = 0.702 Z_l       ρ = 0.5 (1 + r/(150 t))^-0.5
 *   bending            ψ = 1     ξ = 0.702 Z_l       ρ = 0.5 (1 + r/(300 t))^-0.5
 *   torsion and shear  ψ = 5.34  ξ = 0.856 Z_l^(3/4) ρ = 0.6
 *   lateral pressure   ψ = 4     ξ = 1.04 √Z_l       ρ = 0.6
 *   hydrostatic press. ψ = 2     ξ = 1.04 √Z_l       ρ = 0.6
 * For long shells, l/r > 3.85 √(r/t), the code replaces the pressure values by f_Eh = 0.25 E (t/r)² and the
 * torsion value by f_Eτ = 0.25 E (t/r)^(3/2).
 * @param {object} o
 * @param {number} o.r mid-surface radius [m]
 * @param {number} o.t thickness [m]
 * @param {number} o.l length between effective supports (ring frames, bulkheads) [m]
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} [o.nu] Poisson's ratio (default 0.3)
 * @returns {{Zl: number, long: boolean, fEa: number, fEm: number, fEtau: number, fEhLateral: number,
 *   fEhHydrostatic: number, C: {axial: number, bending: number, torsion: number, lateral: number,
 *   hydrostatic: number}}} curvature parameter [-], long-shell flag, elastic buckling strengths [Pa] for axial
 *   compression, bending, torsion/shear, lateral pressure and hydrostatic pressure (hoop stress at buckling),
 *   and the reduced buckling coefficients [-]
 */
export function c202ElasticBuckling({ r, t, l, E, nu = 0.3 }) {
  const fn = 'c202ElasticBuckling';
  pos(r, 'r', fn); pos(t, 't', fn); pos(l, 'l', fn); pos(E, 'E', fn); poisson(nu, fn);
  const Zl = ((l * l) / (r * t)) * Math.sqrt(1 - nu * nu);
  const base = ((Math.PI * Math.PI * E) / (12 * (1 - nu * nu))) * (t / l) ** 2;
  const coef = (psi, xi, rho) => psi * Math.sqrt(1 + ((rho * xi) / psi) ** 2);
  const C = {
    axial: coef(1, 0.702 * Zl, 0.5 / Math.sqrt(1 + r / (150 * t))),
    bending: coef(1, 0.702 * Zl, 0.5 / Math.sqrt(1 + r / (300 * t))),
    torsion: coef(5.34, 0.856 * Zl ** 0.75, 0.6),
    lateral: coef(4, 1.04 * Math.sqrt(Zl), 0.6),
    hydrostatic: coef(2, 1.04 * Math.sqrt(Zl), 0.6),
  };
  const long = l / r > 3.85 * Math.sqrt(r / t);
  return {
    Zl, long,
    fEa: C.axial * base,
    fEm: C.bending * base,
    fEtau: long ? 0.25 * E * (t / r) ** 1.5 : C.torsion * base,
    fEhLateral: long ? 0.25 * E * (t / r) ** 2 : C.lateral * base,
    fEhHydrostatic: long ? 0.25 * E * (t / r) ** 2 : C.hydrostatic * base,
    C,
  };
}

/**
 * Shell buckling check of an unstiffened circular cylinder under combined loading per DNV-RP-C202 §3.2:
 *   λs² = (fy / σj) [σa0/f_Ea + σm0/f_Em + σh0/f_Eh + τ/f_Eτ],  f_ks = fy / √(1 + λs⁴),
 *   σj = √((σa + σm)² − (σa + σm) σh + σh² + 3 τ²),
 * where σa0, σm0, σh0 are the compressive parts (zero when the component is tensile), and the design
 * strength f_ksd = f_ks / γM with γM = 1.15 (λs < 0.5), 0.85 + 0.60 λs (0.5 ≤ λs ≤ 1.0), 1.45 (λs > 1.0).
 * Requirement: σj ≤ f_ksd. Column buckling of the whole cylinder (§3.8) is not covered here.
 * @param {object} o
 * @param {number} o.r mid-surface radius [m]
 * @param {number} o.t thickness [m]
 * @param {number} o.l length between supports [m]
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} [o.nu] Poisson's ratio (default 0.3)
 * @param {number} o.fy yield strength [Pa]
 * @param {number} [o.sigmaA] design axial membrane stress [Pa], tension positive
 * @param {number} [o.sigmaM] design bending stress [Pa] (the compressive extreme fibre is checked: enter the
 *   magnitude as a negative number)
 * @param {number} [o.sigmaH] design hoop stress [Pa], tension positive (external pressure: −p r / t)
 * @param {number} [o.tau] design shear stress from torsion and/or shear force [Pa] (magnitude)
 * @param {boolean} [o.hydrostatic] hoop stress comes from hydrostatic (end-capped) pressure (default false:
 *   lateral pressure)
 * @returns {{sigmaJ: number, lambdaS: number, fks: number, gammaM: number, fksd: number, utilisation: number,
 *   elastic: object}} equivalent stress [Pa]; reduced slenderness [-]; characteristic and design buckling
 *   strengths [Pa]; material factor [-]; σj / f_ksd [-]; the elastic strengths from `c202ElasticBuckling`
 */
export function c202ShellCheck(o) {
  const fn = 'c202ShellCheck';
  const { fy, sigmaA = 0, sigmaM = 0, sigmaH = 0, tau = 0, hydrostatic = false } = o;
  pos(fy, 'fy', fn);
  for (const [v, n] of [[sigmaA, 'sigmaA'], [sigmaM, 'sigmaM'], [sigmaH, 'sigmaH'], [tau, 'tau']]) num(v, n, fn);
  const el = c202ElasticBuckling(o);
  const sa0 = sigmaA < 0 ? -sigmaA : 0;
  const sm0 = sigmaM < 0 ? -sigmaM : 0;
  const sh0 = sigmaH < 0 ? -sigmaH : 0;
  const s = sigmaA + sigmaM;
  const sj = Math.sqrt(s * s - s * sigmaH + sigmaH * sigmaH + 3 * tau * tau);
  if (!(sj > 0)) throw new RangeError(`${fn}: at least one stress component must be non-zero`);
  const fEh = hydrostatic ? el.fEhHydrostatic : el.fEhLateral;
  const lam2 = (fy / sj) * (sa0 / el.fEa + sm0 / el.fEm + sh0 / fEh + Math.abs(tau) / el.fEtau);
  const lambdaS = Math.sqrt(lam2);
  const fks = fy / Math.sqrt(1 + lam2 * lam2);
  let gammaM;
  if (lambdaS < 0.5) gammaM = 1.15; else if (lambdaS <= 1) gammaM = 0.85 + 0.6 * lambdaS; else gammaM = 1.45;
  const fksd = fks / gammaM;
  return { sigmaJ: sj, lambdaS, fks, gammaM, fksd, utilisation: sj / fksd, elastic: el };
}

/**
 * Equivalent cylinder for the buckling check of an unstiffened conical shell (DNV-RP-C202 §6.2): the cone
 * between two supports is treated as a cylinder with r_e = (r1 + r2) / (2 cos α) and l_e = l / cos α, where l
 * is the axial length and α the half apex angle. Membrane stresses in the cone wall follow from equilibrium:
 * meridional σ = N / (2π r t cos α) for an axial force N, hoop σ = p r / (t cos α) for lateral pressure.
 * Stated validity of the code approach: α ≤ 30°.
 * @param {object} o
 * @param {number} o.r1 mid-surface radius at one end [m]
 * @param {number} o.r2 mid-surface radius at the other end [m]
 * @param {number} o.l axial length between the supports [m]
 * @returns {{alpha: number, alphaDeg: number, rEquivalent: number, lEquivalent: number, withinCode: boolean}}
 *   half apex angle [rad, deg]; equivalent radius and length [m]; whether α ≤ 30°
 */
export function coneEquivalentCylinder({ r1, r2, l }) {
  const fn = 'coneEquivalentCylinder';
  pos(r1, 'r1', fn); pos(r2, 'r2', fn); pos(l, 'l', fn);
  const alpha = Math.atan2(Math.abs(r2 - r1), l);
  const c = Math.cos(alpha);
  return { alpha, alphaDeg: (alpha * 180) / Math.PI, rEquivalent: (r1 + r2) / (2 * c), lEquivalent: l / c, withinCode: alpha <= Math.PI / 6 + 1e-12 };
}

/**
 * Collapse pressure of a pipe under external pressure including initial ovality — the cubic of
 * DNV-ST-F101 Eq. (5.11) (Haagsma & Schaap 1981):
 *   (pc − pel) (pc² − pp²) = pc pel pp f0 D/t,
 *   pel = 2 E (t/D)³ / (1 − ν²),  pp = fy αfab 2 t / D,  f0 = (Dmax − Dmin)/D.
 * The physically relevant root (0 < pc ≤ min(pel, pp)) is returned. The code requires f0 ≥ 0.005.
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.t wall thickness [m]
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} [o.nu] Poisson's ratio (default 0.3)
 * @param {number} o.fy yield strength [Pa]
 * @param {number} [o.f0] ovality (Dmax − Dmin)/D [-] (default 0.005)
 * @param {number} [o.alphaFab] fabrication factor (default 1.0; DNV-ST-F101 Table 5-5: seamless 1.00, UO/TRB/ERW
 *   0.93, UOE 0.85)
 * @returns {{pc: number, pel: number, pp: number}} collapse, elastic collapse and plastic collapse pressures [Pa]
 */
export function pipeCollapsePressure({ D, t, E, nu = 0.3, fy, f0 = 0.005, alphaFab = 1 }) {
  const fn = 'pipeCollapsePressure';
  pos(D, 'D', fn); pos(t, 't', fn); pos(E, 'E', fn); poisson(nu, fn); pos(fy, 'fy', fn); pos(alphaFab, 'alphaFab', fn);
  if (!(t < D / 2)) throw new RangeError(`${fn}: t must be smaller than D/2`);
  if (!(f0 >= 0 && f0 < 0.2)) throw new RangeError(`${fn}: f0 must lie in [0, 0.2)`);
  const pel = (2 * E * (t / D) ** 3) / (1 - nu * nu);
  const pp = (fy * alphaFab * 2 * t) / D;
  const g = (p) => (p - pel) * (p * p - pp * pp) - (p * pel * pp * f0 * D) / t;
  let lo = 0;
  let hi = Math.min(pel, pp);
  if (f0 === 0) return { pc: hi, pel, pp };
  // g(0) = pel pp² > 0 and g(min(pel, pp)) < 0: a single sign change in between.
  for (let i = 0; i < 200; i++) {
    const mid = 0.5 * (lo + hi);
    if (g(mid) > 0) lo = mid; else hi = mid;
    if (hi - lo <= 1e-14 * hi) break;
  }
  return { pc: 0.5 * (lo + hi), pel, pp };
}

/**
 * External pressure at first yield of an initially oval ring or long tube (Timoshenko & Gere §7.4): with
 * m = r/t, n = w0/r (initial radial out-of-roundness amplitude over radius) and p_cr = E (t/r)³ / (4 (1 − ν²)),
 * the extreme-fibre stress p m + 6 p m² n / (1 − p/p_cr) reaches fy when
 *   p² − [fy/m + (1 + 6 m n) p_cr] p + fy p_cr / m = 0   (smaller root).
 * Use E/(1 − ν²) → E (set nu = 0) for a ring rather than a long tube.
 * @param {object} o
 * @param {number} o.r mean radius [m]
 * @param {number} o.t wall thickness [m]
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} [o.nu] Poisson's ratio (default 0.3)
 * @param {number} o.fy yield strength [Pa]
 * @param {number} o.w0 initial radial deviation amplitude from the mean circle [m] (ovality f0 = 4 w0 / D)
 * @returns {{pYield: number, pCr: number, pSquash: number}} first-yield pressure, elastic buckling pressure
 *   and membrane yield pressure fy t / r [Pa]
 */
export function ovalTubeYieldPressure({ r, t, E, nu = 0.3, fy, w0 }) {
  const fn = 'ovalTubeYieldPressure';
  pos(r, 'r', fn); pos(t, 't', fn); pos(E, 'E', fn); poisson(nu, fn); pos(fy, 'fy', fn);
  if (!(w0 >= 0)) throw new RangeError(`${fn}: w0 must be >= 0`);
  const m = r / t;
  const n = w0 / r;
  const pCr = (E * (t / r) ** 3) / (4 * (1 - nu * nu));
  const b = fy / m + (1 + 6 * m * n) * pCr;
  const c = (fy * pCr) / m;
  // Smaller root in the numerically stable form.
  const disc = Math.sqrt(Math.max(0, b * b - 4 * c));
  return { pYield: (2 * c) / (b + disc), pCr, pSquash: fy / m };
}

/**
 * Brazier ovalisation instability of a long thin tube in pure bending. With the flattening parameter
 * ζ = κ² r⁴ (1 − ν²) / t², the moment is M = π E r³ t κ (1 − 3ζ/2); it peaks at ζ = 2/9:
 *   M_B = (2 √2 / 9) π E r t² / √(1 − ν²),  κ_B = (√2 / 3) t / (r² √(1 − ν²)).
 * Local (bifurcation) buckling on the compression side normally occurs somewhat earlier; M_B is an upper
 * bound for elastic tubes and irrelevant for thick tubes that yield first.
 * @param {object} o
 * @param {number} o.r mean radius [m]
 * @param {number} o.t wall thickness [m]
 * @param {number} o.E Young's modulus [Pa]
 * @param {number} [o.nu] Poisson's ratio (default 0.3)
 * @param {number} [o.kappa] curvature at which to evaluate the moment [1/m] (optional)
 * @returns {{momentLimit: number, curvatureLimit: number, flatteningLimit: number, moment: number|null,
 *   flattening: number|null}} limit moment [N·m], curvature at the limit [1/m], ζ at the limit (2/9) [-];
 *   moment [N·m] and ζ [-] at the given curvature (valid up to the limit)
 */
export function brazierMoment({ r, t, E, nu = 0.3, kappa }) {
  const fn = 'brazierMoment';
  pos(r, 'r', fn); pos(t, 't', fn); pos(E, 'E', fn); poisson(nu, fn);
  const s = Math.sqrt(1 - nu * nu);
  const out = {
    momentLimit: ((2 * Math.SQRT2) / 9) * Math.PI * E * r * t * t / s,
    curvatureLimit: ((Math.SQRT2 / 3) * t) / (r * r * s),
    flatteningLimit: 2 / 9, moment: null, flattening: null,
  };
  if (kappa !== undefined) {
    if (!(kappa >= 0)) throw new RangeError(`${fn}: kappa must be >= 0`);
    const z = (kappa * kappa * r ** 4 * (1 - nu * nu)) / (t * t);
    out.flattening = z;
    out.moment = Math.PI * E * r ** 3 * t * kappa * (1 - 1.5 * z);
  }
  return out;
}
