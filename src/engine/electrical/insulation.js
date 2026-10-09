/**
 * Cable insulation: electric stress, DC field distribution, dielectric quantities, ageing and breakdown
 * statistics, and insulation co-ordination levels.
 *
 * References
 *  - Coaxial field E(r) = U/(r·ln(R/r0)) (Gauss's law).
 *  - Eoll (1975), IEEE Trans. Electr. Insul. EI-10; Jeroense & Morshuis (1998), IEEE Trans. Diel. Electr. Insul. 5:
 *    DC field with σ = σ0·exp(α·T + β·E).
 *  - IEC 60216-1 (thermal endurance, Arrhenius); Montsinger (1930) halving-interval rule; Dakin (1948).
 *  - IEC 60505 / inverse power law for electrical ageing; Simoni (1981) multi-stress models (the simple
 *    multiplicative form is used here).
 *  - IEC 62539 (Weibull statistics of breakdown data); volume/length enlargement law.
 *  - IEC 60071-1:2019 Table 2 and 3, IEC 60502-2, IEC 60840, IEC 62067 (rated impulse withstand voltages).
 *
 * Notes not covered by formulas: (1) water treeing — wet-design XLPE insulation without a radial water barrier
 * ages by water trees; the design stress is then kept low and the design is qualified by long-term wet ageing
 * tests (CIGRE TB 722, IEC 63026); no life model is given here, set the stress limits of stressCheck accordingly.
 * (2) polarity reversal of DC cables superposes the capacitive (Laplace) field of twice the voltage on the
 * resistive field carrying space charge; dcPolarityReversalStress gives that first-instant estimate.
 *
 * SI units: V, m, V/m, S/m, F/m, s; temperatures in °C (names ending in C) or K where stated.
 */

const EPS0 = 8.8541878128e-12;
const K_BOLTZMANN_EV = 8.617333262e-5; // eV/K

function pos(v, name) {
  if (!(v > 0) || !Number.isFinite(v)) throw new RangeError(`${name}: must be a positive finite number`);
}
function nonNeg(v, name) {
  if (!(v >= 0) || !Number.isFinite(v)) throw new RangeError(`${name}: must be a non-negative finite number`);
}

/**
 * AC (capacitive) electric stress in coaxial insulation: E(r) = U0/(r·ln(ro/ri)); maximum at the conductor
 * screen (r = ri), minimum at the insulation screen (r = ro).
 * @param {object} o
 * @param {number} o.U0 voltage across the insulation (phase to earth, rms or peak consistently), V
 * @param {number} o.ri radius over the conductor screen, m
 * @param {number} o.ro radius over the insulation, m
 * @param {number} [o.r] radius at which the stress is wanted, m (default: none)
 * @returns {{Emax:number, Emin:number, Emean:number, E:?number}} V/m (E at r when given)
 */
export function acStress({ U0, ri, ro, r }) {
  nonNeg(U0, 'U0'); pos(ri, 'ri'); pos(ro, 'ro');
  if (ro <= ri) throw new RangeError('ro: must exceed ri');
  const ln = Math.log(ro / ri);
  if (r !== undefined && !(r >= ri && r <= ro)) throw new RangeError('r: must lie within the insulation');
  return { Emax: U0 / (ri * ln), Emin: U0 / (ro * ln), Emean: U0 / (ro - ri), E: r !== undefined ? U0 / (r * ln) : null };
}

/**
 * Design-stress check: maximum and minimum AC stress against limits supplied by the designer
 * (indicatively, extruded HV/EHV cables are designed for several kV/mm up to roughly 15 kV/mm at the conductor
 * screen at the highest voltages, and wet-design array cables for lower values; use the manufacturer's limits).
 * @param {object} o arguments of acStress plus
 * @param {number} o.limitInner permissible stress at the conductor screen, V/m
 * @param {number} [o.limitOuter] permissible stress at the insulation screen, V/m
 * @returns {{Emax:number, Emin:number, utilisationInner:number, utilisationOuter:number, ok:boolean,
 *            thicknessRequired:number}} V/m; thicknessRequired is the insulation thickness (m) that would meet
 *          limitInner on the same conductor screen radius, t = ri·(exp(U0/(ri·limitInner)) − 1)
 */
export function stressCheck({ U0, ri, ro, limitInner, limitOuter }) {
  pos(limitInner, 'limitInner');
  const s = acStress({ U0, ri, ro });
  const utilisationInner = s.Emax / limitInner;
  const utilisationOuter = limitOuter ? s.Emin / limitOuter : 0;
  return { Emax: s.Emax, Emin: s.Emin, utilisationInner, utilisationOuter, ok: utilisationInner <= 1 && utilisationOuter <= 1,
    thicknessRequired: ri * (Math.exp(U0 / (ri * limitInner)) - 1) };
}

/**
 * Capacitance, insulation resistance and dielectric loss of a coaxial core:
 * C = 2π·ε0·εr/ln(ro/ri) per metre; Rins = ρ/(2π)·ln(ro/ri) for one metre (Ω·m; divide by length);
 * loss angle δ = atan(tanδ); loss Wd = ω·C·U0²·tanδ; the self-discharge time constant is ρ·ε0·εr.
 * @param {object} o
 * @param {number} o.ri m @param {number} o.ro m @param {number} o.epsR relative permittivity
 * @param {number} [o.resistivity] volume resistivity of the insulation, Ω·m
 * @param {number} [o.tanDelta=0] @param {number} [o.U0] V @param {number} [o.f] Hz
 * @returns {{C:number, Rins:number, lossAngle:number, Wd:number, timeConstant:number}} F/m, Ω·m (resistance of a
 *          1 m length), rad, W/m, s (NaN where inputs are missing)
 */
export function dielectricProperties({ ri, ro, epsR, resistivity, tanDelta = 0, U0, f }) {
  pos(ri, 'ri'); pos(ro, 'ro'); pos(epsR, 'epsR'); nonNeg(tanDelta, 'tanDelta');
  if (ro <= ri) throw new RangeError('ro: must exceed ri');
  const ln = Math.log(ro / ri);
  const C = (2 * Math.PI * EPS0 * epsR) / ln;
  return { C, Rins: resistivity ? (resistivity * ln) / (2 * Math.PI) : NaN, lossAngle: Math.atan(tanDelta),
    Wd: U0 !== undefined && f !== undefined ? 2 * Math.PI * f * C * U0 * U0 * tanDelta : NaN,
    timeConstant: resistivity ? resistivity * EPS0 * epsR : NaN };
}

/**
 * Steady-state DC field in coaxial insulation with temperature- and field-dependent conductivity
 * σ(T, E) = σ0·exp(α·T + β·E). Current continuity gives r·σ·E = constant; the temperature falls logarithmically
 * from the conductor, T(r) = Tc − ΔT·ln(r/ri)/ln(ro/ri). The field is found by Newton iteration at each radius
 * for a trial leakage constant and bisection on that constant so that ∫E·dr = U.
 * For β = 0 the closed form E(r) = δ·U·r^(δ−1)/(ro^δ − ri^δ), δ = α·ΔT/ln(ro/ri), is recovered; stress inversion
 * (maximum at the insulation screen) occurs for δ > 1 (δ including the field term, Eoll's approximation:
 * δ = (α·ΔT/ln(ro/ri) + β·U/(ro − ri))/(1 + β·U/(ro − ri)), returned for reference).
 * Typical coefficients: mass-impregnated paper α ≈ 0.1 1/K, β ≈ 0.03 mm/kV; DC-XLPE α ≈ 0.05–0.1 1/K,
 * β ≈ 0.03–0.13 mm/kV (material-specific, from the manufacturer).
 * Reference: Eoll (1975); Jeroense & Morshuis (1998).
 * @param {object} o
 * @param {number} o.U DC voltage across the insulation, V
 * @param {number} o.ri m @param {number} o.ro m
 * @param {number} o.alpha temperature coefficient of conductivity, 1/K
 * @param {number} o.beta field coefficient of conductivity, m/V
 * @param {number} [o.deltaT=0] temperature drop across the insulation (conductor − screen), K
 * @param {number} [o.n=201] radial points
 * @returns {{r:Float64Array, E:Float64Array, Einner:number, Eouter:number, Emax:number, rMax:number,
 *            inverted:boolean, deltaEoll:number, Elaplace:number}} radii (m), field (V/m), field at the two
 *          screens, maximum and its radius, inversion flag (Eouter > Einner), Eoll's exponent and the
 *          capacitive (no-load, Laplace) field at the conductor screen for comparison
 */
export function dcField({ U, ri, ro, alpha, beta, deltaT = 0, n = 201 }) {
  pos(U, 'U'); pos(ri, 'ri'); pos(ro, 'ro'); nonNeg(alpha, 'alpha'); nonNeg(beta, 'beta'); nonNeg(deltaT, 'deltaT');
  if (ro <= ri) throw new RangeError('ro: must exceed ri');
  if (!(n >= 11)) throw new RangeError('n: at least 11 points');
  const ln = Math.log(ro / ri);
  const r = new Float64Array(n); const E = new Float64Array(n); const g = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    // Logarithmic spacing resolves the steep region near the conductor.
    r[i] = ri * Math.exp((ln * i) / (n - 1));
    // g = ln of the relative resistivity factor from temperature: σ ∝ exp(−α·ΔT·ln(r/ri)/ln(ro/ri)).
    g[i] = (alpha * deltaT * Math.log(r[i] / ri)) / ln;
  }
  const solveE = (rhs) => {
    // E·exp(β·E) = rhs, Newton on f(E) = ln E + β·E − ln rhs (monotonic, concave).
    if (beta === 0) return rhs;
    const lr = Math.log(rhs);
    let e = Math.min(rhs, lr > 0 ? lr / beta + 1 / beta : rhs);
    if (!(e > 0)) e = rhs;
    for (let it = 0; it < 60; it++) {
      const fv = Math.log(e) + beta * e - lr;
      const d = 1 / e + beta;
      let en = e - fv / d;
      if (!(en > 0)) en = e / 10;
      if (Math.abs(en - e) < 1e-13 * e) { e = en; break; }
      e = en;
    }
    return e;
  };
  const integral = (lnK) => {
    for (let i = 0; i < n; i++) E[i] = solveE(Math.exp(lnK + g[i]) / r[i]);
    let s = 0;
    for (let i = 1; i < n; i++) s += 0.5 * (E[i] + E[i - 1]) * (r[i] - r[i - 1]);
    return s;
  };
  let lo = Math.log(U / ln) - 60; let hi = Math.log(U / ln) + 60;
  while (integral(hi) < U) hi += 60;
  while (integral(lo) > U) lo -= 60;
  for (let it = 0; it < 200 && hi - lo > 1e-14; it++) {
    const mid = 0.5 * (lo + hi);
    if (integral(mid) < U) lo = mid; else hi = mid;
  }
  integral(0.5 * (lo + hi));
  let Emax = 0; let rMax = ri;
  for (let i = 0; i < n; i++) if (E[i] > Emax) { Emax = E[i]; rMax = r[i]; }
  const bu = (beta * U) / (ro - ri);
  return { r, E, Einner: E[0], Eouter: E[n - 1], Emax, rMax, inverted: E[n - 1] > E[0],
    deltaEoll: ((alpha * deltaT) / ln + bu) / (1 + bu), Elaplace: U / (ri * ln) };
}

/**
 * First-instant stress at the conductor screen after a fast polarity reversal of a loaded DC cable: the
 * resistive field (with its space charge) of the old polarity plus the capacitive field of the voltage step 2U,
 * E = |Edc(ri) − 2·U/(ri·ln(ro/ri))| at the inner screen and |Edc(ro) − 2·U/(ro·ln(ro/ri))| at the outer.
 * Reference: Jeroense & Morshuis (1998); CIGRE TB 496 commentary.
 * @param {object} o
 * @param {number} o.U voltage before the reversal, V @param {number} o.ri m @param {number} o.ro m
 * @param {number} o.Einner resistive field at the conductor screen before the reversal, V/m
 * @param {number} o.Eouter resistive field at the insulation screen before the reversal, V/m
 * @returns {{Einner:number, Eouter:number}} V/m (magnitudes just after the reversal)
 */
export function dcPolarityReversalStress({ U, ri, ro, Einner, Eouter }) {
  pos(U, 'U'); pos(ri, 'ri'); pos(ro, 'ro');
  if (ro <= ri) throw new RangeError('ro: must exceed ri');
  const ln = Math.log(ro / ri);
  return { Einner: Math.abs(Einner - (2 * U) / (ri * ln)), Eouter: Math.abs(Eouter - (2 * U) / (ro * ln)) };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Ageing                                                                                                       */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Thermal life by the Arrhenius relation, L(T) = L0·exp[(Ea/k)·(1/T − 1/T0)] (absolute temperatures), or, when
 * `halvingIntervalK` is given instead of the activation energy, by Montsinger's rule L = L0·2^(−(θ − θ0)/Δ).
 * Reference: IEC 60216-1; Dakin (1948); Montsinger (1930). Typical values: Ea ≈ 1.0–1.3 eV for XLPE
 * (material-specific); halving interval 8–10 K.
 * @param {object} o
 * @param {number} o.thetaC operating temperature, °C
 * @param {number} o.thetaRefC reference temperature at which the life is lifeRef, °C
 * @param {number} o.lifeRef life at the reference temperature (any time unit)
 * @param {number} [o.activationEnergyEv] activation energy, eV
 * @param {number} [o.halvingIntervalK] temperature rise that halves the life, K
 * @returns {number} life at thetaC, same unit as lifeRef
 */
export function thermalLife({ thetaC, thetaRefC, lifeRef, activationEnergyEv, halvingIntervalK }) {
  pos(lifeRef, 'lifeRef');
  if (!Number.isFinite(thetaC) || !Number.isFinite(thetaRefC)) throw new RangeError('thetaC, thetaRefC: must be finite');
  if (activationEnergyEv !== undefined) {
    pos(activationEnergyEv, 'activationEnergyEv');
    return lifeRef * Math.exp((activationEnergyEv / K_BOLTZMANN_EV) * (1 / (thetaC + 273.15) - 1 / (thetaRefC + 273.15)));
  }
  if (halvingIntervalK === undefined) throw new RangeError('activationEnergyEv or halvingIntervalK: one is required');
  pos(halvingIntervalK, 'halvingIntervalK');
  return lifeRef * Math.pow(2, -(thetaC - thetaRefC) / halvingIntervalK);
}

/**
 * Electrical life by the inverse power law, L(E) = L0·(E/E0)^(−n).
 * Reference: IEC 60505. The exponent is material- and design-specific (values of roughly 9–15 are commonly
 * assumed for extruded AC insulation and about 10 for extruded DC insulation, CIGRE TB 496); it is an argument.
 * @param {object} o
 * @param {number} o.E operating stress (or voltage), any unit @param {number} o.Eref reference stress, same unit
 * @param {number} o.lifeRef life at the reference stress @param {number} o.n life exponent
 * @returns {number} life, same unit as lifeRef
 */
export function electricalLife({ E, Eref, lifeRef, n }) {
  pos(E, 'E'); pos(Eref, 'Eref'); pos(lifeRef, 'lifeRef'); pos(n, 'n');
  return lifeRef * Math.pow(E / Eref, -n);
}

/**
 * Combined electro-thermal life as the product of the two acceleration factors (no synergy term):
 * L = L0·(E/E0)^(−n)·AF_T⁻¹, with the thermal part from thermalLife. A screening model; multi-stress models with
 * interaction (Simoni, Crine) need material test data.
 * @param {object} o arguments of thermalLife and electricalLife (lifeRef is the life at Eref and thetaRefC)
 * @returns {number} life, same unit as lifeRef
 */
export function combinedLife({ E, Eref, n, thetaC, thetaRefC, lifeRef, activationEnergyEv, halvingIntervalK }) {
  const thermal = thermalLife({ thetaC, thetaRefC, lifeRef, activationEnergyEv, halvingIntervalK });
  return electricalLife({ E, Eref, lifeRef: thermal, n });
}

/**
 * Loss of life accumulated over a temperature history (Miner-type summation of life fractions):
 * LoL = Σ Δt_i/L(θ_i), optionally with a constant electrical stress factor. Takes the conductor temperature
 * history from thermaltransient.temperatureHistory directly.
 * @param {object} o
 * @param {ArrayLike<number>} o.thetaC temperature of each step, °C
 * @param {number} o.dt step length, in the time unit of lifeRef (e.g. hours if lifeRef is in hours)
 * @param {number} o.thetaRefC °C @param {number} o.lifeRef life at the reference conditions
 * @param {number} [o.activationEnergyEv] @param {number} [o.halvingIntervalK]
 * @param {number} [o.E] @param {number} [o.Eref] @param {number} [o.n] optional electrical stress terms
 * @returns {{lossOfLife:number, equivalentAgeingFactor:number, lifeAtThisDuty:number}} fraction of life consumed,
 *          mean ageing rate relative to the reference condition, and the life if this duty repeated indefinitely
 */
export function lossOfLife({ thetaC, dt, thetaRefC, lifeRef, activationEnergyEv, halvingIntervalK, E, Eref, n }) {
  const N = thetaC.length;
  if (!(N >= 1)) throw new RangeError('thetaC: must be a non-empty array');
  pos(dt, 'dt');
  const elec = E !== undefined ? Math.pow(E / Eref, n) : 1;
  let sum = 0;
  for (let i = 0; i < N; i++) {
    sum += (dt * elec) / thermalLife({ thetaC: thetaC[i], thetaRefC, lifeRef, activationEnergyEv, halvingIntervalK });
  }
  const duration = N * dt;
  return { lossOfLife: sum, equivalentAgeingFactor: (sum * lifeRef) / duration, lifeAtThisDuty: duration / sum };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Breakdown statistics                                                                                         */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Two-parameter Weibull breakdown probability with time and size terms:
 * P = 1 − exp[−(E/E0)^b·(t/t0)^a·(L/L0)], where E0 is the 63.2 % breakdown stress of a sample of length (or
 * volume) L0 tested for time t0.
 * Reference: IEC 62539; Occhini (1971) for the enlargement law applied to cables.
 * @param {object} o
 * @param {number} o.E stress, V/m @param {number} o.E0 scale stress of the reference sample, V/m
 * @param {number} o.b shape parameter for stress (extruded cables typically 8–20)
 * @param {number} [o.sizeRatio=1] L/L0 (length or stressed volume ratio)
 * @param {number} [o.timeRatio=1] t/t0 @param {number} [o.a=1] shape parameter for time
 * @returns {number} probability of breakdown (0–1)
 */
export function weibullBreakdownProbability({ E, E0, b, sizeRatio = 1, timeRatio = 1, a = 1 }) {
  nonNeg(E, 'E'); pos(E0, 'E0'); pos(b, 'b'); pos(sizeRatio, 'sizeRatio'); nonNeg(timeRatio, 'timeRatio'); pos(a, 'a');
  return -Math.expm1(-Math.pow(E / E0, b) * Math.pow(timeRatio, a) * sizeRatio);
}

/**
 * Enlargement law: stress giving the same breakdown probability for a larger insulation (route length against
 * test length) and longer time, E_L = E0·(L0/L)^(1/b)·(t0/t)^(a/b); also the stress for a target probability,
 * E_p = E0·[−ln(1 − p)/((L/L0)·(t/t0)^a)]^(1/b).
 * Reference: IEC 62539; Weibull weakest-link scaling.
 * @param {object} o
 * @param {number} o.E0 63.2 % stress of the reference sample, V/m @param {number} o.b stress shape parameter
 * @param {number} o.sizeRatio L/L0 @param {number} [o.timeRatio=1] @param {number} [o.a=1]
 * @param {number} [o.probability] target breakdown probability (0–1)
 * @returns {{E632:number, Eprobability:?number}} 63.2 % stress of the enlarged insulation and stress at the
 *          target probability, V/m
 */
export function weibullEnlargement({ E0, b, sizeRatio, timeRatio = 1, a = 1, probability }) {
  pos(E0, 'E0'); pos(b, 'b'); pos(sizeRatio, 'sizeRatio'); pos(timeRatio, 'timeRatio'); pos(a, 'a');
  const k = sizeRatio * Math.pow(timeRatio, a);
  let Ep = null;
  if (probability !== undefined) {
    if (!(probability > 0 && probability < 1)) throw new RangeError('probability: must be in (0, 1)');
    Ep = E0 * Math.pow(-Math.log1p(-probability) / k, 1 / b);
  }
  return { E632: E0 * Math.pow(1 / k, 1 / b), Eprobability: Ep };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Insulation co-ordination                                                                                     */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Rated lightning-impulse (BIL) and switching-impulse withstand voltages by highest voltage for equipment Um,
 * as used for cable systems: IEC 60502-2 (Um ≤ 36 kV), IEC 60840 (36 < Um ≤ 170 kV), IEC 62067 (Um > 170 kV);
 * values are members of the IEC 60071-1 standard series. Where IEC 60071-1 offers several levels for one Um the
 * value normally specified for cables is listed. Editable: pass your own table to impulseLevels.
 * Um and withstand voltages in volts (peak for impulses).
 */
export const IMPULSE_LEVELS = Object.freeze([
  Object.freeze({ Um: 7.2e3, lightning: 60e3, switching: null }),
  Object.freeze({ Um: 12e3, lightning: 75e3, switching: null }),
  Object.freeze({ Um: 17.5e3, lightning: 95e3, switching: null }),
  Object.freeze({ Um: 24e3, lightning: 125e3, switching: null }),
  Object.freeze({ Um: 36e3, lightning: 170e3, switching: null }),
  Object.freeze({ Um: 52e3, lightning: 250e3, switching: null }),
  Object.freeze({ Um: 72.5e3, lightning: 325e3, switching: null }),
  Object.freeze({ Um: 123e3, lightning: 550e3, switching: null }),
  Object.freeze({ Um: 145e3, lightning: 650e3, switching: null }),
  Object.freeze({ Um: 170e3, lightning: 750e3, switching: null }),
  Object.freeze({ Um: 245e3, lightning: 1050e3, switching: null }),
  Object.freeze({ Um: 300e3, lightning: 1050e3, switching: 850e3 }),
  Object.freeze({ Um: 362e3, lightning: 1175e3, switching: 950e3 }),
  Object.freeze({ Um: 420e3, lightning: 1425e3, switching: 1050e3 }),
  Object.freeze({ Um: 550e3, lightning: 1550e3, switching: 1175e3 }),
]);

/**
 * Impulse withstand levels for a system voltage and a simple co-ordination check against the protective level
 * of the surge arrester: margin = withstand/(protective level) must exceed the co-ordination/safety factor
 * (IEC 60071-2 uses a safety factor of 1.15 for internal insulation on top of the co-ordination withstand
 * voltage; the required ratio is an argument).
 * @param {object} o
 * @param {number} o.Um highest voltage for equipment (line-to-line, rms), V
 * @param {{Um:number, lightning:number, switching:?number}[]} [o.table=IMPULSE_LEVELS] editable level table
 * @param {number} [o.arresterLightningProtectiveLevel] V (peak)
 * @param {number} [o.arresterSwitchingProtectiveLevel] V (peak)
 * @param {number} [o.requiredRatio=1.15]
 * @returns {{Um:number, lightning:number, switching:?number, lightningMargin:?number, switchingMargin:?number,
 *            ok:?boolean}} the row for the smallest tabulated Um ≥ the requested one, with margins when arrester
 *          levels are given
 */
export function impulseLevels({ Um, table = IMPULSE_LEVELS, arresterLightningProtectiveLevel, arresterSwitchingProtectiveLevel, requiredRatio = 1.15 }) {
  pos(Um, 'Um');
  const sorted = [...table].sort((a, b) => a.Um - b.Um);
  const row = sorted.find((x) => x.Um >= Um * (1 - 1e-9));
  if (!row) throw new RangeError('Um: above the highest voltage in the table');
  const lightningMargin = arresterLightningProtectiveLevel ? row.lightning / arresterLightningProtectiveLevel : null;
  const switchingMargin = arresterSwitchingProtectiveLevel && row.switching ? row.switching / arresterSwitchingProtectiveLevel : null;
  let ok = null;
  if (lightningMargin !== null || switchingMargin !== null) {
    ok = (lightningMargin === null || lightningMargin >= requiredRatio) && (switchingMargin === null || switchingMargin >= requiredRatio);
  }
  return { Um: row.Um, lightning: row.lightning, switching: row.switching, lightningMargin, switchingMargin, ok };
}
