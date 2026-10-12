/**
 * Scour at offshore and coastal structures, and scour-protection sizing.
 *
 * SI units: m, s, m/s, Pa, kg/m^3. KC = Um T / D (Um orbital velocity amplitude at the bed).
 *
 * References
 *   Sumer B.M., Fredsøe J. & Christiansen N. (1992) Scour around vertical pile in waves, JWPCOE 118(1).
 *   Sumer B.M., Christiansen N. & Fredsøe J. (1992) Time scale of scour around a vertical pile, ISOPE.
 *   Sumer B.M. & Fredsøe J. (2001) Scour around pile in combined waves and current, J. Hydraul. Eng. 127(5).
 *   Sumer B.M. & Fredsøe J. (2002) The Mechanics of Scour in the Marine Environment, World Scientific.
 *   Breusers H.N.C., Nicollet G. & Shen H.W. (1977) Local scour around cylindrical piers, J. Hydraul. Res. 15.
 *   FHWA HEC-18 (5th ed., 2012) Evaluating Scour at Bridges, ch. 7.
 *   Sumer B.M. & Fredsøe J. (1990) Scour below pipelines in waves, JWPCOE 116(3).
 *   Sumer B.M. et al. (2001) Onset of scour below pipelines and self-burial, Coastal Eng. 42.
 *   Rudolph D. & Bos K.J. (2006) Scour around a monopile under combined wave-current conditions, ICSE-3.
 *   Soulsby R.L. (1997) Dynamics of Marine Sands, Thomas Telford.
 *   DNV-ST-0126 (2021) 7.7 / DNV-RP-0618; API RP 2A-WSD 6.3 commentary.
 *   Terzaghi K. & Peck R.B. (1948) filter rules; CIRIA C683 The Rock Manual 5.4.
 */

const G = 9.80665;

function need(name, v, min = -Infinity, max = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < min || v > max) {
    throw new RangeError(`${name} must be a number in [${min}, ${max}] (got ${v})`);
  }
  return v;
}

/**
 * Keulegan–Carpenter number KC = Um T / D.
 * @param {number} Um near-bed orbital velocity amplitude [m/s]
 * @param {number} T wave period [s]
 * @param {number} D diameter [m]
 * @returns {number} KC [-]
 */
export function keuleganCarpenter(Um, T, D) {
  need('Um', Um, 0);
  need('T', T, 1e-9);
  need('D', D, 1e-9);
  return (Um * T) / D;
}

/**
 * Equilibrium scour depth at a slender vertical pile.
 *   Waves only (Sumer et al. 1992):          S/D = 1.3 {1 − exp[−0.03 (KC − 6)]},  KC ≥ 6, else 0.
 *   Waves + current (Sumer & Fredsøe 2001):  S/D = 1.3 {1 − exp[−A (KC − B)]},  KC ≥ B,
 *       A = 0.03 + 0.75 Ucw^2.6,  B = 6 exp(−4.7 Ucw),  Ucw = Uc / (Uc + Um).
 *   Steady current (Ucw = 1 or KC → ∞):      S/D = 1.3 with standard deviation σ(S/D) = 0.7; the design
 *       value is 1.3 + nSigma·0.7 (DNV-ST-0126 quotes 1.3 D as the default when no other data exist).
 * Range of the combined formula. It was fitted to tests at KC of about 4 to 26 and expresses how a current helps
 * the wave-driven horseshoe vortex; it does not contain the scour the current produces on its own. Outside
 * that range it under-predicts badly: a large monopile in a current with small waves has KC ≪ 1, and the
 * expression returns S of a few per cent of D where the current alone gives 1.3 D. The function therefore
 * selects a regime and reports it:
 *   'steady-current'     KC = ∞ or Ucw = 1: S/D = 1.3 (+ nSigma·0.7);
 *   'current-dominated'  Ucw ≥ currentDominatedUcw (default 0.7): the steady-current value (Sumer & Fredsøe
 *                        2001, J. Hydraul. Eng. 127(5): for Ucw above about 0.7 the scour depth approaches its
 *                        steady-current value);
 *   'current-below-wave-threshold'  KC < B with a current present (0 < Ucw < currentDominatedUcw): the waves are
 *                        below their scour threshold, the steady-current value is returned;
 *   'no-wave-scour'      KC < 6 without current (Ucw = 0): S = 0;
 *   'wave'               Ucw = 0, KC ≥ 6: Sumer et al. (1992);
 *   'wave-current'       otherwise: Sumer & Fredsøe (2001), limited from above by the steady-current value.
 * The two current regimes assume the current alone keeps the bed live (θ > θcr): check the mobility, and use
 * `fallback: false` for the bare formula when the current is too weak to move the bed.
 * Live-bed conditions (θ > θcr), D/L < ~0.1 (no diffraction).
 * @param {object} o
 * @param {number} [o.KC=Infinity] Keulegan–Carpenter number [-]
 * @param {number} [o.Ucw=0] current ratio Uc/(Uc + Um) [-], 0 (waves) … 1 (current)
 * @param {number} [o.nSigma=0] number of standard deviations added to the mean equilibrium value [-]
 * @param {boolean} [o.fallback=true] apply the steady-current value in the two current regimes above; false
 *   evaluates the exponential formula for every KC and Ucw (regime 'formula')
 * @param {number} [o.currentDominatedUcw=0.7] current ratio from which the flow counts as current-dominated [-]
 * @returns {{SoverD:number,A:number,B:number,regime:string,formulaSoverD:number}} S/D [-], coefficients [-], the
 *   regime used, and the value of the exponential formula itself for comparison [-]
 */
export function pileScour({ KC = Infinity, Ucw = 0, nSigma = 0, fallback = true, currentDominatedUcw = 0.7 } = {}) {
  need('KC', KC, 0, Infinity);
  need('Ucw', Ucw, 0, 1);
  need('nSigma', nSigma, 0, 5);
  need('currentDominatedUcw', currentDominatedUcw, 0, 1);
  const A = 0.03 + 0.75 * Ucw ** 2.6;
  const B = 6 * Math.exp(-4.7 * Ucw);
  const eq = 1.3 + 0.7 * nSigma;
  if (KC === Infinity || Ucw === 1) return { SoverD: eq, A, B, regime: 'steady-current', formulaSoverD: eq };
  const formula = KC > B ? eq * (1 - Math.exp(-A * (KC - B))) : 0;
  let SoverD = formula;
  let regime;
  if (!fallback) regime = 'formula';
  else if (Ucw >= currentDominatedUcw) { regime = 'current-dominated'; SoverD = eq; }
  else if (!(KC > B)) {
    if (Ucw > 0) { regime = 'current-below-wave-threshold'; SoverD = eq; } else regime = 'no-wave-scour';
  } else regime = Ucw > 0 ? 'wave-current' : 'wave';
  return { SoverD, A, B, regime, formulaSoverD: formula };
}

/**
 * Coefficients proposed by Rudolph & Bos (2006) for the same exponential form at low KC (1 < KC < 10),
 * typical of large monopiles:  A = 0.03 + 1.5 Ucw⁴,  B = 6 exp(−5 Ucw).
 * Quoted from the conference paper; confirm before use in final design.
 * @param {object} o KC [-], Ucw [-]
 * @returns {{SoverD:number,A:number,B:number}} [-]
 */
export function pileScourRudolphBos({ KC, Ucw }) {
  need('KC', KC, 0);
  need('Ucw', Ucw, 0, 1);
  const A = 0.03 + 1.5 * Ucw ** 4;
  const B = 6 * Math.exp(-5 * Ucw);
  return { SoverD: KC > B ? 1.3 * (1 - Math.exp(-A * (KC - B))) : 0, A, B };
}

/**
 * Breusers, Nicollet & Shen (1977) pier scour in a steady current: S = 1.5 K D tanh(h/D), K the product of
 * the shape, alignment, velocity-intensity and grouping factors (1 for an isolated circular pile in
 * live-bed flow). A design coefficient of 2.0 instead of 1.5 is often adopted (argument).
 * @param {object} o
 * @param {number} o.D pile diameter [m]
 * @param {number} o.h water depth [m]
 * @param {number} [o.K=1] product of correction factors [-]
 * @param {number} [o.coefficient=1.5] [-]
 * @returns {number} scour depth S [m]
 */
export function scourBreusers({ D, h, K = 1, coefficient = 1.5 }) {
  need('D', D, 1e-9);
  need('h', h, 0);
  need('K', K, 0);
  return coefficient * K * D * Math.tanh(h / D);
}

/**
 * HEC-18 (CSU) pier-scour equation (5th ed. eq. 7.1/7.3):
 *   ys / a = 2.0 K1 K2 K3 (y1/a)^0.35 Fr1^0.43,   Fr1 = V1 / √(g y1),
 *   K1 nose shape: square 1.1, round 1.0, circular cylinder 1.0, group of cylinders 1.0, sharp 0.9;
 *   K2 = (cos θ + (L/a) sin θ)^0.65 for flow skew θ (L/a ≤ 12);
 *   K3 bed condition: clear-water 1.1, plane bed/antidunes 1.1, small dunes 1.1, medium dunes 1.1–1.2,
 *      large dunes 1.3;  K4 (4th-ed. armouring factor, ≥ 0.4) optional.
 * Round-nosed piers aligned with the flow: ys ≤ 2.4 a (Fr ≤ 0.8) or 3.0 a (Fr > 0.8).
 * @param {object} o
 * @param {number} o.a pier width [m]
 * @param {number} o.y1 approach flow depth [m]
 * @param {number} o.V1 approach velocity [m/s]
 * @param {'square'|'round'|'circular'|'group'|'sharp'} [o.nose='circular']
 * @param {number} [o.skew=0] angle of attack θ [rad]
 * @param {number} [o.L=o.a] pier length [m]
 * @param {number} [o.K3=1.1] bed-condition factor [-]
 * @param {number} [o.K4=1] armouring factor [-]
 * @param {boolean} [o.applyLimit=true] apply the 2.4a / 3.0a cap (round/circular, zero skew)
 * @returns {{ys:number,Fr:number,K1:number,K2:number}} ys scour depth [m]
 */
export function scourHec18({ a, y1, V1, nose = 'circular', skew = 0, L, K3 = 1.1, K4 = 1, applyLimit = true }) {
  need('a', a, 1e-9);
  need('y1', y1, 1e-9);
  need('V1', V1, 0);
  need('skew', skew, 0, Math.PI / 2);
  const K1 = { square: 1.1, round: 1.0, circular: 1.0, group: 1.0, sharp: 0.9 }[nose];
  if (K1 === undefined) throw new RangeError(`nose must be square, round, circular, group or sharp (got ${nose})`);
  const La = Math.min((L === undefined ? a : L) / a, 12);
  const K2 = (Math.cos(skew) + La * Math.sin(skew)) ** 0.65;
  const Fr = V1 / Math.sqrt(G * y1);
  let ys = 2.0 * (skew > (5 * Math.PI) / 180 ? 1 : K1) * K2 * K3 * K4 * a * (y1 / a) ** 0.35 * Fr ** 0.43;
  if (applyLimit && skew === 0 && (nose === 'round' || nose === 'circular')) ys = Math.min(ys, (Fr <= 0.8 ? 2.4 : 3.0) * a);
  return { ys, Fr, K1, K2 };
}

/**
 * Sheppard–Melville pier-scour method as given in HEC-18 5th ed. §7.3 (eqs 7.5–7.9), ae = effective pier
 * width, Vc = critical velocity for the bed sediment (argument):
 *   f1 = tanh[(y1/ae)^0.4],  f2 = 1 − 1.2 [ln(V1/Vc)]²,
 *   f3 = (ae/D50) / [0.4 (ae/D50)^1.2 + 10.6 (ae/D50)^−0.13];
 *   0.4 ≤ V1/Vc < 1:        ys/ae = 2.5 f1 f2 f3;
 *   1 ≤ V1/Vc ≤ Vlp/Vc:     ys/ae = f1 [2.2 (V1/Vc − 1)/(Vlp/Vc − 1) + 2.5 f3 (Vlp/Vc − V1/Vc)/(Vlp/Vc − 1)];
 *   V1/Vc > Vlp/Vc:         ys/ae = 2.2 f1;
 *   Vlp = max(5 Vc, 0.6 √(g y1)).
 * @param {object} o a [m] effective pier width, y1 [m], V1 [m/s], Vc [m/s], D50 [m]
 * @returns {{ys:number,f1:number,f2:number,f3:number,Vlp:number}} ys [m]
 */
export function scourSheppardMelville({ a, y1, V1, Vc, D50 }) {
  need('a', a, 1e-9);
  need('y1', y1, 1e-9);
  need('V1', V1, 0);
  need('Vc', Vc, 1e-9);
  need('D50', D50, 1e-9);
  const r = V1 / Vc;
  const f1 = Math.tanh((y1 / a) ** 0.4);
  const f2 = r > 0 ? 1 - 1.2 * Math.log(r) ** 2 : 0;
  const ad = a / D50;
  const f3 = ad / (0.4 * ad ** 1.2 + 10.6 * ad ** -0.13);
  const Vlp = Math.max(5 * Vc, 0.6 * Math.sqrt(G * y1));
  const rl = Vlp / Vc;
  let ys;
  if (r < 0.4) ys = 0;
  else if (r < 1) ys = 2.5 * f1 * f2 * f3 * a;
  else if (r <= rl) ys = f1 * ((2.2 * (r - 1)) / (rl - 1) + (2.5 * f3 * (rl - r)) / (rl - 1)) * a;
  else ys = 2.2 * f1 * a;
  return { ys, f1, f2, f3, Vlp };
}

/**
 * Time scale of scour at a pile (Sumer, Christiansen & Fredsøe 1992; Sumer & Fredsøe 2002 §3.3) and the
 * exponential development S(t) = S∞ [1 − exp(−t/T)]:
 *   current:  T* = (δ/D) θ^−2.2 / 2000   (δ = boundary-layer thickness, usually the water depth);
 *   waves:    T* = 10^−6 (KC/θ)³;
 *   T = T* D² / √(g (s − 1) d³).
 * θ is the undisturbed Shields parameter (live bed).
 * @param {object} o
 * @param {'current'|'waves'} o.mode
 * @param {number} o.D pile diameter [m]
 * @param {number} o.theta Shields parameter [-]
 * @param {number} o.d grain size d50 [m]
 * @param {number} [o.s=2.65] relative density of the sediment [-]
 * @param {number} [o.delta] boundary-layer thickness [m] (current)
 * @param {number} [o.KC] [-] (waves)
 * @returns {{Tstar:number,T:number}} T [s]
 */
export function scourTimeScale({ mode, D, theta, d, s = 2.65, delta, KC }) {
  need('D', D, 1e-9);
  need('theta', theta, 1e-9);
  need('d', d, 1e-9);
  need('s', s, 1 + 1e-9);
  let Tstar;
  if (mode === 'current') Tstar = (need('delta', delta, 1e-9) / D / 2000) * theta ** -2.2;
  else if (mode === 'waves') Tstar = 1e-6 * (need('KC', KC, 0) / theta) ** 3;
  else throw new RangeError(`mode must be 'current' or 'waves' (got ${mode})`);
  return { Tstar, T: (Tstar * D * D) / Math.sqrt(G * (s - 1) * d ** 3) };
}

/**
 * Scour depth at time t: S(t) = S∞ [1 − exp(−t/T)].
 * @param {number} Seq equilibrium depth [m]
 * @param {number} t time [s]
 * @param {number} T time scale [s]
 * @returns {number} S(t) [m]
 */
export function scourDevelopment(Seq, t, T) {
  need('t', t, 0);
  need('T', T, 1e-12);
  return Seq * (1 - Math.exp(-t / T));
}

/**
 * Scour below a pipeline.
 *   Waves (Sumer & Fredsøe 1990):  S/D = 0.1 √KC  (2 < KC < 300, live bed);
 *   Current (live bed; Sumer & Fredsøe 2002 §2.3):  S/D = 0.6 ± 0.1, and with an initial gap e between the
 *   pipe and the bed  S/D = 0.625 exp(−0.6 e/D)  (Hansen et al. 1986).
 * @param {object} o
 * @param {'waves'|'current'} o.mode
 * @param {number} [o.KC] [-] (waves)
 * @param {number} [o.gapRatio=0] e/D [-] (current)
 * @returns {number} S/D [-]
 */
export function pipelineScour({ mode, KC, gapRatio = 0 }) {
  if (mode === 'waves') return 0.1 * Math.sqrt(need('KC', KC, 0));
  if (mode === 'current') return 0.625 * Math.exp(-0.6 * need('gapRatio', gapRatio, 0));
  throw new RangeError(`mode must be 'waves' or 'current' (got ${mode})`);
}

/**
 * Onset of scour (piping) below a partly embedded pipeline in a current (Sumer et al. 2001):
 *   Ucr² / [g D (1 − n)(s − 1)] = 0.025 exp[9 (e/D)^0.5],   e = embedment depth.
 * @param {object} o
 * @param {number} o.D pipe diameter [m]
 * @param {number} o.e embedment [m], 0 ≤ e ≤ D/2
 * @param {number} [o.n=0.4] bed porosity [-]
 * @param {number} [o.s=2.65] relative grain density [-]
 * @param {number} [o.U] current velocity at the pipe [m/s] for the verdict
 * @returns {{Ucr:number,onset:boolean|undefined}} critical velocity [m/s]
 */
export function pipelineScourOnset({ D, e, n = 0.4, s = 2.65, U }) {
  need('D', D, 1e-9);
  need('e', e, 0, D);
  need('n', n, 0, 0.9);
  const Ucr = Math.sqrt(G * D * (1 - n) * (s - 1) * 0.025 * Math.exp(9 * Math.sqrt(e / D)));
  return { Ucr, onset: U === undefined ? undefined : U > Ucr };
}

/**
 * Scour at the edge of large structures (gravity bases, mudmat edges, large caissons). No closed formula of
 * general validity exists: in the diffraction regime (D/L > ~0.1, KC < ~1) scour is governed by steady
 * streaming and by the edge geometry. This helper applies a user-supplied non-dimensional depth,
 *   S = kS · Lref,
 * with kS taken from model tests or published guidance for the configuration (e.g. Rudolph & Bos 2006 or
 * Khalfin 1983 for cylinders; Whitehouse 1998 for gravity bases), and reports the regime indicators.
 * @param {object} o
 * @param {number} o.Lref reference dimension (diameter, skirt depth or edge height) [m]
 * @param {number} o.kS non-dimensional scour depth S/Lref [-]
 * @param {number} [o.KC] [-]
 * @param {number} [o.DoverL] diameter over wavelength [-]
 * @returns {{S:number,regime:string}} S [m]
 */
export function largeStructureScour({ Lref, kS, KC, DoverL }) {
  need('Lref', Lref, 0);
  need('kS', kS, 0);
  let regime = 'unspecified';
  if (DoverL !== undefined || KC !== undefined) {
    regime = (DoverL !== undefined && DoverL > 0.1) || (KC !== undefined && KC < 1) ? 'diffraction / streaming' : KC !== undefined && KC < 6 ? 'low-KC slender' : 'vortex shedding (slender)';
  }
  return { S: kS * Lref, regime };
}

/**
 * Global (general) scour around a multi-leg structure – documented rule: the global depression is
 * expressed as a fraction of a reference dimension supplied by the user (site experience, model tests or
 * the operator's standard), S_global = kG · Dref, and added to the local scour at each leg.
 * @param {object} o
 * @param {number} o.Dref reference dimension (leg or pile diameter) [m]
 * @param {number} o.kG global scour ratio [-]
 * @param {number} [o.localScour=0] local scour depth at the leg [m]
 * @returns {{global:number,total:number}} [m]
 */
export function globalScour({ Dref, kG, localScour = 0 }) {
  need('Dref', Dref, 0);
  need('kG', kG, 0);
  need('localScour', localScour, 0);
  return { global: kG * Dref, total: kG * Dref + localScour };
}

/**
 * Dimensionless grain size D* = d [g (s − 1)/ν²]^(1/3) and critical Shields parameter (Soulsby 1997 eq. 77):
 *   θcr = 0.30/(1 + 1.2 D*) + 0.055 [1 − exp(−0.020 D*)].
 * @param {object} o
 * @param {number} o.d grain or stone size [m]
 * @param {number} [o.s=2.65] relative density ρs/ρ [-]
 * @param {number} [o.nu=1.36e-6] kinematic viscosity [m^2/s] (sea water, 10 °C)
 * @returns {{Dstar:number,thetaCr:number,tauCr:number}} D* [-]; θcr [-]; tauCr = θcr (s − 1) g d is the
 *   kinematic critical shear stress τcr/ρ [m^2/s^2] (multiply by the water density for Pa)
 */
export function criticalShields({ d, s = 2.65, nu = 1.36e-6 }) {
  need('d', d, 1e-12);
  need('s', s, 1 + 1e-9);
  need('nu', nu, 1e-12);
  const Dstar = d * ((G * (s - 1)) / (nu * nu)) ** (1 / 3);
  const thetaCr = 0.3 / (1 + 1.2 * Dstar) + 0.055 * (1 - Math.exp(-0.02 * Dstar));
  return { Dstar, thetaCr, tauCr: thetaCr * (s - 1) * G * d };
}

/**
 * Bed shear stress under current, waves and their combination on a flat rough bed (Soulsby 1997):
 *   z0 = ks/30, ks = 2.5 d;  current  τc = ρ Cd U²,  Cd = [0.40 / (ln(h/z0) − 1)]²  (eq. 37);
 *   waves  τw = ½ ρ fw Uw²,  fw = 1.39 (A/z0)^−0.52,  A = Uw T / 2π  (eq. 62a);
 *   combined  τm = τc [1 + 1.2 (τw/(τc + τw))^3.2],  τmax = √[(τm + τw cos φ)² + (τw sin φ)²]  (eqs 69–70).
 * @param {object} o
 * @param {number} o.d bed grain / stone size [m]
 * @param {number} o.h water depth [m]
 * @param {number} [o.U=0] depth-averaged current [m/s]
 * @param {number} [o.Uw=0] orbital velocity amplitude [m/s]
 * @param {number} [o.T] wave period [s] (needed when Uw > 0)
 * @param {number} [o.phi=0] angle between current and waves [rad]
 * @param {number} [o.rho=1025] water density [kg/m^3]
 * @returns {{tauC:number,tauW:number,tauMean:number,tauMax:number,Cd:number,fw:number}} stresses [Pa]
 */
export function bedShearStress({ d, h, U = 0, Uw = 0, T, phi = 0, rho = 1025 }) {
  need('d', d, 1e-12);
  need('h', h, 1e-9);
  const z0 = (2.5 * d) / 30;
  if (h <= Math.E * z0) throw new RangeError('h must be much larger than the bed roughness');
  const Cd = (0.4 / (Math.log(h / z0) - 1)) ** 2;
  const tauC = rho * Cd * U * U;
  let tauW = 0;
  let fw = 0;
  if (Uw > 0) {
    const A = (Uw * need('T', T, 1e-9)) / (2 * Math.PI);
    fw = 1.39 * (A / z0) ** -0.52;
    tauW = 0.5 * rho * fw * Uw * Uw;
  }
  const tauMean = tauC > 0 ? tauC * (1 + 1.2 * (tauW / (tauC + tauW)) ** 3.2) : 0;
  const tauMax = Math.hypot(tauMean + tauW * Math.cos(phi), tauW * Math.sin(phi));
  return { tauC, tauW, tauMean, tauMax, Cd, fw };
}

/**
 * Stone size for statically stable scour protection by the critical-shear (Shields) approach: smallest d50
 * for which  α · τmax(d) ≤ θcr(d) (ρs − ρ) g d,  with α the amplification of bed shear stress next to the
 * structure (argument; about 4 for a circular pile in current, lower in waves at small KC – Sumer &
 * Fredsøe 2002; Whitehouse 1998 quotes design values in the range 2–4).
 * @param {object} o as bedShearStress (h, U, Uw, T, phi, rho) plus
 * @param {number} o.amplification α [-]
 * @param {number} [o.rhoS=2650] stone density [kg/m^3]
 * @param {number} [o.nu=1.36e-6] [m^2/s]
 * @returns {{d50:number,thetaCr:number,tauMax:number,tauCr:number}} d50 [m]; stresses [Pa]
 */
export function rockSizeShields({ h, U = 0, Uw = 0, T, phi = 0, rho = 1025, amplification, rhoS = 2650, nu = 1.36e-6 }) {
  need('amplification', amplification, 1e-9);
  need('rhoS', rhoS, rho * (1 + 1e-9));
  const s = rhoS / rho;
  const margin = (d) => {
    const tb = bedShearStress({ d, h, U, Uw, T, phi, rho }).tauMax * amplification;
    return criticalShields({ d, s, nu }).tauCr * rho - tb;
  };
  let lo = 1e-5;
  let hi = Math.min(h / 10, 5);
  if (margin(hi) < 0) throw new RangeError('no stable stone size below h/10: loading too severe for the Shields approach');
  if (margin(lo) >= 0) hi = lo;
  else for (let i = 0; i < 80; i++) { const mid = Math.sqrt(lo * hi); if (margin(mid) < 0) lo = mid; else hi = mid; }
  const cs = criticalShields({ d: hi, s, nu });
  return { d50: hi, thetaCr: cs.thetaCr, tauMax: bedShearStress({ d: hi, h, U, Uw, T, phi, rho }).tauMax * amplification, tauCr: cs.tauCr * rho };
}

/**
 * Geometrically closed granular filter rules (Terzaghi; CIRIA C683 §5.4.3.6):
 *   retention      D15f / D85b < 5;
 *   permeability   D15f / D15b > 4 (limit argument, 4–5 quoted);
 *   uniformity     D50f / D50b < 25.
 * @param {object} o
 * @param {number} o.D15f @param {number} o.D50f filter sizes [m]
 * @param {number} o.D15b @param {number} o.D50b @param {number} o.D85b base sizes [m]
 * @param {number} [o.permeabilityLimit=4] [-]
 * @returns {{retention:{ratio:number,ok:boolean},permeability:{ratio:number,ok:boolean},
 *   uniformity:{ratio:number,ok:boolean},ok:boolean}}
 */
export function filterCriteria({ D15f, D50f, D15b, D50b, D85b, permeabilityLimit = 4 }) {
  for (const [k, v] of Object.entries({ D15f, D50f, D15b, D50b, D85b })) need(k, v, 1e-12);
  const r1 = D15f / D85b;
  const r2 = D15f / D15b;
  const r3 = D50f / D50b;
  const res = {
    retention: { ratio: r1, ok: r1 < 5 },
    permeability: { ratio: r2, ok: r2 > permeabilityLimit },
    uniformity: { ratio: r3, ok: r3 < 25 },
  };
  res.ok = res.retention.ok && res.permeability.ok && res.uniformity.ok;
  return res;
}

/**
 * Quantity take-off for an annular scour-protection layer around a pile: outer radius = extentFactor · D
 * from the pile axis, thickness t (commonly ≥ 2–3 Dn50), bulk porosity n.
 * @param {object} o
 * @param {number} o.D pile diameter [m]
 * @param {number} o.extentFactor outer radius / D [-] (> 0.5)
 * @param {number} o.thickness [m]
 * @param {number} [o.porosity=0.4] [-]
 * @param {number} [o.rhoS=2650] stone density [kg/m^3]
 * @returns {{area:number,volume:number,mass:number,outerRadius:number}} [m^2, m^3 bulk, kg, m]
 */
export function armourTakeoff({ D, extentFactor, thickness, porosity = 0.4, rhoS = 2650 }) {
  need('D', D, 1e-9);
  need('extentFactor', extentFactor, 0.5);
  need('thickness', thickness, 0);
  need('porosity', porosity, 0, 0.9);
  const ro = extentFactor * D;
  const area = Math.PI * (ro * ro - (D * D) / 4);
  return { area, volume: area * thickness, mass: area * thickness * (1 - porosity) * rhoS, outerRadius: ro };
}

/**
 * Effect of scour on a pile or skirted foundation: reduced embedment and effective overburden profile.
 * General scour Sg lowers the mudline uniformly. Local scour Sl removes soil at the foundation only; its
 * influence on vertical effective stress is taken to fade linearly to zero over a depth `reductionDepth`
 * below the base of the local hole (API RP 2A-WSD C6.3 illustrates 6 D; argument), below which the
 * overburden is that of the general-scour mudline.
 * @param {object} o
 * @param {number} o.embedment original embedded length [m]
 * @param {number} [o.generalScour=0] Sg [m]
 * @param {number} [o.localScour=0] Sl [m]
 * @param {number} o.gammaSub submerged unit weight [N/m^3]
 * @param {number} o.reductionDepth depth below the local scour base over which the stress recovers [m]
 * @returns {{embedment:number,scourBase:number,sigmaV:(z:number)=>number}} reduced embedment [m]; depth of
 *   the scour base below the original mudline [m]; σ′v(z) [Pa] with z below the original mudline
 */
export function scourEffect({ embedment, generalScour = 0, localScour = 0, gammaSub, reductionDepth }) {
  need('embedment', embedment, 0);
  need('generalScour', generalScour, 0);
  need('localScour', localScour, 0);
  need('gammaSub', gammaSub, 0);
  need('reductionDepth', reductionDepth, 0);
  const base = generalScour + localScour;
  const zFull = base + reductionDepth;
  const sigmaV = (z) => {
    if (z <= base) return 0;
    const full = gammaSub * (z - generalScour);
    if (z >= zFull || reductionDepth === 0) return full;
    return (gammaSub * (zFull - generalScour) * (z - base)) / reductionDepth;
  };
  return { embedment: Math.max(embedment - base, 0), scourBase: base, sigmaV };
}
