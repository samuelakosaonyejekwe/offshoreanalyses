// Oscillating-water-column wave-energy converter: rigid-piston water column coupled to a linearised isentropic
// air chamber and an air turbine (Wells: linear; impulse: quadratic), regular and irregular waves, optimum
// turbine damping, power matrix and annual energy.
//
// Units: SI. Chamber pressure p is gauge (Pa), volume flow q in m³/s, turbine conductance G_t = q/p in m³/(s·Pa).

import { cylinderHeave, waveSpectrum, matrixEnergy, regularWavePower } from './wec.js';

const G = 9.80665;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}
function finitePos(v, name) {
  need(Number.isFinite(v) && v > 0, `${name} must be a positive finite number`);
}

/**
 * Natural (piston-mode) period of the water column, T = 2π √((d + l_a)/g), with d the column draft and l_a an
 * added length representing the entrained water below the mouth. For a circular mouth of radius a the added
 * length lies between 4a/(3π) (half the added mass of an unbaffled disc) and 8a/(3π) (piston in an infinite
 * baffle); it is an argument. Reference: Evans (1978) J. Inst. Maths Applics 22:423–433; Falcão & Henriques
 * (2016) Renewable Energy 85:1391–1424.
 * @param {object} o
 * @param {number} o.draft column draft d (m)
 * @param {number} [o.addedLength=0] added length l_a (m)
 * @returns {{period:number, omega:number}} s, rad/s
 */
export function columnResonance({ draft, addedLength = 0 }) {
  finitePos(draft, 'draft');
  need(addedLength >= 0, 'addedLength must be >= 0');
  const omega = Math.sqrt(G / (draft + addedLength));
  return { period: (2 * Math.PI) / omega, omega };
}

/**
 * Linearised air-chamber compliance for isentropic compression, C = V0 / (γ p_atm) (m³/Pa): the chamber mass
 * balance S ẋ = q_turbine + C ṗ. Reference: Sarmento & Falcão (1985) J. Fluid Mech. 150:467–485.
 * @param {object} o
 * @param {number} o.volume still-water chamber air volume V0 (m³)
 * @param {number} [o.patm=101325] atmospheric pressure (Pa)
 * @param {number} [o.gamma=1.4] ratio of specific heats (-)
 * @returns {number} C (m³/Pa)
 */
export function chamberCompliance({ volume, patm = 101325, gamma = 1.4 }) {
  finitePos(volume, 'volume');
  finitePos(patm, 'patm');
  finitePos(gamma, 'gamma');
  return volume / (gamma * patm);
}

/**
 * Coefficients of a circular water column treated as a rigid piston: mass ρ S (d + l_a), hydrostatic stiffness
 * ρ g S, Froude–Krylov excitation on the column mouth and Haskind radiation damping (see `cylinderHeave`).
 * @param {object} o
 * @param {number} o.radius column radius (m)
 * @param {number} o.draft column draft (m)
 * @param {number} [o.addedLength=0] added length (m)
 * @param {number} [o.depth=Infinity] water depth (m)
 * @param {number} [o.rho=1025]
 * @param {number} [o.excitationFactor=1] multiplier on the Froude–Krylov excitation (-)
 * @param {number} [o.Bvisc=0] linear loss damping, e.g. entry losses (N·s/m)
 * @returns {{area:number, column:(omega:number)=>{M:number, B:number, K:number, Fe:number}}}
 */
export function columnBody({ radius, draft, addedLength = 0, depth = Infinity, rho = 1025, excitationFactor = 1, Bvisc = 0 }) {
  finitePos(radius, 'radius');
  finitePos(draft, 'draft');
  need(addedLength >= 0, 'addedLength must be >= 0');
  const S = Math.PI * radius * radius;
  const M = rho * S * (draft + addedLength);
  return {
    area: S,
    column: (omega) => {
      const c = cylinderHeave({ radius, draft, omega, depth, rho, diffractionFactor: excitationFactor });
      return { M, B: c.B + Bvisc, K: c.K, Fe: c.Fe };
    },
  };
}

// complex transfer functions per unit wave amplitude at one frequency
function transfer(col, S, omega, Gt, C) {
  // air impedance seen by the piston: i ω S² / (G_t + i ω C)
  const d = Gt * Gt + omega * omega * C * C;
  const airRe = (omega * omega * C * S * S) / d; // stiffness-like part
  const airIm = (omega * S * S * Gt) / d; // damping-like part (× 1, already includes ω)
  const re = col.K - omega * omega * col.M + airRe;
  const im = omega * col.B + airIm;
  const x = col.Fe / Math.hypot(re, im);
  const p = (omega * S * x) / Math.sqrt(d);
  return { x, p };
}

function meanEfficiency(eta, sigmaOrAmp, harmonic) {
  if (typeof eta !== 'function') return eta;
  // power-weighted mean of η(|q|): weight q² over a sinusoid (regular) or a Gaussian (irregular)
  let num = 0;
  let den = 0;
  const n = 200;
  for (let i = 0; i < n; i++) {
    let q;
    let w;
    if (harmonic) {
      q = sigmaOrAmp * Math.sin((Math.PI * (i + 0.5)) / n);
      w = q * q;
    } else {
      q = (5 * sigmaOrAmp * (i + 0.5)) / n;
      w = q * q * Math.exp((-0.5 * q * q) / (sigmaOrAmp * sigmaOrAmp));
    }
    num += eta(Math.abs(q)) * w;
    den += w;
  }
  return den > 0 ? num / den : 0;
}

/**
 * OWC in a regular wave.
 * Water column (rigid piston, displacement x positive upward): M ẍ + B ẋ + ρ g S x = F_e ζ − p S.
 * Chamber: S ẋ = q + C ṗ. Turbine: Wells — linear, q = G_t p (G_t = 1/damping); impulse — quadratic,
 * p = k_i q |q|, replaced by the energy-equivalent linear turbine G_t = 1 / ((8/(3π)) k_i |q̂|) (iterated).
 * Pneumatic power available to the turbine P = ½ G_t |p̂|²; electrical power = η P with the turbine efficiency
 * either a number or a function of the instantaneous flow rate |q| (power-weighted over the cycle).
 * @param {object} o
 * @param {(omega:number)=>{M:number, B:number, K:number, Fe:number}} o.column from `columnBody`
 * @param {number} o.area column cross-section S (m²)
 * @param {number} o.omega angular frequency (rad/s)
 * @param {number} o.amplitude wave amplitude (m)
 * @param {number} o.compliance chamber compliance C (m³/Pa), see `chamberCompliance`; 0 = incompressible
 * @param {'wells'|'impulse'} [o.turbine='wells']
 * @param {number} [o.conductance] Wells turbine G_t = q/p (m³/(s·Pa))
 * @param {number} [o.ki] impulse-turbine coefficient k_i in p = k_i q|q| (Pa·s²/m⁶)
 * @param {number|((q:number)=>number)} [o.efficiency=1] turbine–generator efficiency (-) or η(|q|) with q in m³/s
 * @returns {{pneumaticPower:number, electricalPower:number, pressure:number, flow:number, stroke:number,
 *   conductance:number, efficiency:number}} W, W, Pa amplitude, m³/s amplitude, m amplitude
 */
export function owcRegular({ column, area, omega, amplitude, compliance, turbine = 'wells', conductance, ki, efficiency = 1 }) {
  finitePos(area, 'area');
  finitePos(omega, 'omega');
  need(amplitude >= 0, 'amplitude must be >= 0');
  need(compliance >= 0, 'compliance must be >= 0');
  const col = column(omega);
  let Gt = conductance;
  let r;
  if (turbine === 'wells') {
    finitePos(Gt, 'conductance');
    r = transfer(col, area, omega, Gt, compliance);
  } else {
    need(turbine === 'impulse', "turbine must be 'wells' or 'impulse'");
    finitePos(ki, 'ki');
    if (amplitude === 0) return { pneumaticPower: 0, electricalPower: 0, pressure: 0, flow: 0, stroke: 0, conductance: Infinity, efficiency: 0 };
    // solve q̂ = G_t(q̂) p̂(G_t) by bisection on ln q̂ (residual is monotonic)
    const f = (lq) => {
      const qh = Math.exp(lq);
      const g = 1 / ((8 / (3 * Math.PI)) * ki * qh);
      return g * transfer(col, area, omega, g, compliance).p * amplitude - qh;
    };
    let lo = -30;
    let hi = 30;
    for (let i = 0; i < 200; i++) {
      const mid = 0.5 * (lo + hi);
      if (f(mid) > 0) lo = mid;
      else hi = mid;
    }
    Gt = 1 / ((8 / (3 * Math.PI)) * ki * Math.exp(0.5 * (lo + hi)));
    r = transfer(col, area, omega, Gt, compliance);
  }
  const p = r.p * amplitude;
  const q = Gt * p;
  const P = 0.5 * Gt * p * p;
  const eta = meanEfficiency(efficiency, q, true);
  return { pneumaticPower: P, electricalPower: eta * P, pressure: p, flow: q, stroke: r.x * amplitude, conductance: Gt, efficiency: eta };
}

/**
 * Wells-turbine conductance maximising the pneumatic power in a regular wave (golden-section search on ln G_t),
 * with the corresponding capture width.
 * @param {object} o as for `owcRegular` without `conductance`, plus
 * @param {number} [o.depth=Infinity] water depth (m) for the incident power
 * @param {number} [o.rho=1025]
 * @returns {{conductance:number, damping:number, pneumaticPower:number, captureWidth:number}} G_t (m³/(s·Pa)),
 *   1/G_t (Pa·s/m³), W, m
 */
export function optimalTurbineDamping({ depth = Infinity, rho = 1025, ...o }) {
  const f = (lg) => owcRegular({ ...o, turbine: 'wells', conductance: Math.exp(lg), efficiency: 1 }).pneumaticPower;
  let a = Math.log(1e-8);
  let b = Math.log(1e4);
  // coarse scan to bracket the maximum, then golden section
  let bi = 0;
  let best = -1;
  const n = 120;
  for (let i = 0; i <= n; i++) {
    const v = f(a + ((b - a) * i) / n);
    if (v > best) {
      best = v;
      bi = i;
    }
  }
  const step = (b - a) / n;
  const c0 = a + step * bi;
  a = c0 - step;
  b = c0 + step;
  const gr = (Math.sqrt(5) - 1) / 2;
  for (let i = 0; i < 80; i++) {
    const c = b - gr * (b - a);
    const d = a + gr * (b - a);
    if (f(c) > f(d)) b = d;
    else a = c;
  }
  const Gt = Math.exp(0.5 * (a + b));
  const P = owcRegular({ ...o, turbine: 'wells', conductance: Gt, efficiency: 1 }).pneumaticPower;
  const J = regularWavePower({ amplitude: o.amplitude, omega: o.omega, depth, rho });
  return { conductance: Gt, damping: 1 / Gt, pneumaticPower: P, captureWidth: J > 0 ? P / J : NaN };
}

/**
 * OWC in an irregular sea state by spectral superposition. Wells turbine: linear. Impulse turbine: statistical
 * linearisation of p = k_i q|q| for Gaussian flow, G_t = 1 / (√(8/π) k_i σ_q) (iterated). Mean pneumatic power
 * P̄ = G_t σ_p².
 * @param {object} o
 * @param {(omega:number)=>{M:number, B:number, K:number, Fe:number}} o.column
 * @param {number} o.area column cross-section (m²)
 * @param {number[]|Float64Array} o.omega frequency grid (rad/s), ascending
 * @param {number[]|Float64Array} o.S wave spectrum (m²·s/rad)
 * @param {number} o.compliance chamber compliance (m³/Pa)
 * @param {'wells'|'impulse'} [o.turbine='wells']
 * @param {number} [o.conductance] Wells G_t (m³/(s·Pa))
 * @param {number} [o.ki] impulse coefficient (Pa·s²/m⁶)
 * @param {number|((q:number)=>number)} [o.efficiency=1]
 * @returns {{pneumaticPower:number, electricalPower:number, sigmaP:number, sigmaQ:number, sigmaX:number,
 *   conductance:number, efficiency:number}}
 */
export function owcIrregular({ column, area, omega, S, compliance, turbine = 'wells', conductance, ki, efficiency = 1 }) {
  need(omega && S && omega.length === S.length && omega.length > 2, 'omega and S must be arrays of equal length > 2');
  finitePos(area, 'area');
  need(compliance >= 0, 'compliance must be >= 0');
  const cols = Array.from(omega, (w) => (w > 0 ? column(w) : null));
  const stats = (Gt) => {
    let vp = 0;
    let vx = 0;
    let pPrev = 0;
    let xPrev = 0;
    for (let i = 0; i < omega.length; i++) {
      let pp = 0;
      let xx = 0;
      if (cols[i] && S[i] > 0) {
        const r = transfer(cols[i], area, omega[i], Gt, compliance);
        pp = r.p * r.p * S[i];
        xx = r.x * r.x * S[i];
      }
      if (i > 0) {
        vp += 0.5 * (pp + pPrev) * (omega[i] - omega[i - 1]);
        vx += 0.5 * (xx + xPrev) * (omega[i] - omega[i - 1]);
      }
      pPrev = pp;
      xPrev = xx;
    }
    return { sp: Math.sqrt(vp), sx: Math.sqrt(vx) };
  };
  let Gt = conductance;
  if (turbine === 'wells') finitePos(Gt, 'conductance');
  else {
    need(turbine === 'impulse', "turbine must be 'wells' or 'impulse'");
    finitePos(ki, 'ki');
    const c = Math.sqrt(8 / Math.PI) * ki;
    const f = (lq) => {
      const sq = Math.exp(lq);
      const g = 1 / (c * sq);
      return g * stats(g).sp - sq;
    };
    if (stats(1).sp === 0) return { pneumaticPower: 0, electricalPower: 0, sigmaP: 0, sigmaQ: 0, sigmaX: 0, conductance: Infinity, efficiency: 0 };
    let lo = -30;
    let hi = 30;
    for (let i = 0; i < 100; i++) {
      const mid = 0.5 * (lo + hi);
      if (f(mid) > 0) lo = mid;
      else hi = mid;
    }
    Gt = 1 / (c * Math.exp(0.5 * (lo + hi)));
  }
  const s = stats(Gt);
  const P = Gt * s.sp * s.sp;
  const eta = meanEfficiency(efficiency, Gt * s.sp, false);
  return { pneumaticPower: P, electricalPower: eta * P, sigmaP: s.sp, sigmaQ: Gt * s.sp, sigmaX: s.sx, conductance: Gt, efficiency: eta };
}

/**
 * OWC power matrix over (Hs, Tp) and annual energy with a scatter diagram.
 * @param {object} o
 * @param {number[]} o.hs significant wave heights (m)
 * @param {number[]} o.tp peak periods (s)
 * @param {(omega:number)=>object} o.column
 * @param {number} o.area (m²)
 * @param {number} o.compliance (m³/Pa)
 * @param {'wells'|'impulse'} [o.turbine='wells']
 * @param {number} [o.conductance] Wells G_t (m³/(s·Pa))
 * @param {number} [o.ki] impulse coefficient
 * @param {number|Function} [o.efficiency=1]
 * @param {number} [o.gamma=1] spectral peak enhancement (-)
 * @param {number} [o.ratedPower=Infinity] electrical rating (W)
 * @param {number[][]} [o.prob] occurrence prob[iHs][iTp]; when given the annual energy is returned
 * @param {number} [o.availability=1]
 * @param {number[]|Float64Array} [o.omega] frequency grid (rad/s); default 0.15–3.0 rad/s, 200 points
 * @returns {{power:number[][], aepWh:number, meanPower:number, capacityFactor:number}} electrical power
 *   power[iHs][iTp] (W); energy fields are NaN without `prob`
 */
export function owcPowerMatrix({ hs, tp, column, area, compliance, turbine = 'wells', conductance, ki, efficiency = 1, gamma = 1, ratedPower = Infinity, prob, availability = 1, omega }) {
  need(Array.isArray(hs) && Array.isArray(tp) && hs.length > 0 && tp.length > 0, 'hs and tp must be non-empty arrays');
  let w = omega;
  if (!w) {
    w = new Float64Array(200);
    for (let i = 0; i < w.length; i++) w[i] = 0.15 + (i * 2.85) / (w.length - 1);
  }
  const power = hs.map((H) =>
    tp.map((T) => {
      if (H === 0) return 0;
      const r = owcIrregular({ column, area, omega: w, S: waveSpectrum({ Hs: H, Tp: T, gamma, omega: w }), compliance, turbine, conductance, ki, efficiency });
      return Math.min(ratedPower, r.electricalPower);
    }),
  );
  if (!prob) return { power, aepWh: NaN, meanPower: NaN, capacityFactor: NaN };
  const e = matrixEnergy({ power, prob, ratedPower: Number.isFinite(ratedPower) ? ratedPower : undefined, availability });
  return { power, aepWh: e.aepWh, meanPower: e.meanPower, capacityFactor: e.capacityFactor };
}

