/**
 * Thermal / pressure expansion of pipelines: restrained axial force, friction build-up from free ends,
 * virtual anchors, end expansion with spool / PLET flexibility, spool stresses, and axial ratcheting
 * (pipeline walking) both by closed-form expressions and by a 1-D elastic bar with Coulomb friction that is
 * cycled through heat-up and cool-down.
 *
 * Sign convention: effective axial force S is tension-positive; x runs from the hot (inlet) end.
 * All quantities SI.
 */

import { solveTridiagonal } from '../math/linalg.js';

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function nonNeg(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}
function fin(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v))) throw new RangeError(`${name} must be a finite number`);
  return v;
}

/**
 * Effective axial force of a fully restrained pipe: S_0 = H − Δp_i A_i (1 − 2ν) − A_s E α ΔT.
 *
 * Reference: DNV-RP-F110 (2018) Eq. 6 / DNV-ST-F101 (2021) Eq. 4.11; thin-wall, linear elastic.
 * H is the residual lay tension, Δp_i the internal pressure difference relative to laying, ΔT the
 * temperature difference relative to laying.
 *
 * @param {object} o
 * @param {number} [o.H=0] residual lay tension [N]
 * @param {number} o.deltaPi internal pressure increase since installation [Pa]
 * @param {number} o.Ai bore area [m²]
 * @param {number} o.As steel area [m²]
 * @param {number} o.deltaT temperature increase since installation [K]
 * @param {number} [o.E=207e9] [Pa]
 * @param {number} [o.alpha=11.7e-6] [1/K]
 * @param {number} [o.nu=0.3] [-]
 * @returns {{S0:number, pressureTerm:number, thermalTerm:number, freeStrain:number}}
 *   S0 [N] (negative = compression); freeStrain = unrestrained axial strain of a capped pipe [-]
 */
export function restrainedAxialForce({ H = 0, deltaPi, Ai, As, deltaT, E = 207e9, alpha = 11.7e-6, nu = 0.3 }) {
  fin('H', H); fin('deltaPi', deltaPi); pos('Ai', Ai); pos('As', As); fin('deltaT', deltaT); pos('E', E); pos('alpha', alpha);
  const pressureTerm = -deltaPi * Ai * (1 - 2 * nu);
  const thermalTerm = -As * E * alpha * deltaT;
  return { S0: H + pressureTerm + thermalTerm, pressureTerm, thermalTerm, freeStrain: -(pressureTerm + thermalTerm) / (E * As) };
}

/**
 * Steady temperature profile along a line with heat loss to the surroundings:
 * T(x) = T_amb + (T_in − T_amb) e^(−x/λ), decay length λ = ṁ c_p/(U π D) when not given directly.
 *
 * Reference: steady energy balance with a constant overall heat-transfer coefficient (no Joule–Thomson term).
 *
 * @param {object} o
 * @param {number|ArrayLike<number>} o.x position(s) from the inlet [m]
 * @param {number} o.TinC inlet temperature [°C]
 * @param {number} o.TambC ambient temperature [°C]
 * @param {number} [o.lambda] decay length [m]
 * @param {number} [o.massFlow] [kg/s]
 * @param {number} [o.cp] specific heat of the fluid [J/(kg·K)]
 * @param {number} [o.U] overall heat-transfer coefficient referred to diameter Dref [W/(m²·K)]
 * @param {number} [o.Dref] reference diameter for U [m]
 * @returns {number|Float64Array} temperature(s) [°C]
 */
export function temperatureProfile({ x, TinC, TambC, lambda, massFlow, cp, U, Dref }) {
  fin('TinC', TinC); fin('TambC', TambC);
  const lam = lambda === undefined ? (pos('massFlow', massFlow) * pos('cp', cp)) / (pos('U', U) * Math.PI * pos('Dref', Dref)) : pos('lambda', lambda);
  const f = (xx) => TambC + (TinC - TambC) * Math.exp(-xx / lam);
  if (typeof x === 'number') return f(x);
  const out = new Float64Array(x.length);
  for (let i = 0; i < x.length; i++) out[i] = f(x[i]);
  return out;
}

/**
 * Internal pressure along the line: p(x) = p_in − (dp/dx)_friction · x + ρ g (z_in − z(x)).
 *
 * Reference: steady single-phase momentum balance with a constant frictional gradient.
 *
 * @param {object} o
 * @param {ArrayLike<number>} o.x positions [m]
 * @param {number} o.pIn inlet pressure [Pa]
 * @param {number} [o.frictionGradient=0] frictional pressure drop per metre [Pa/m]
 * @param {number} [o.rho=0] fluid density [kg/m³]
 * @param {ArrayLike<number>} [o.z] elevation of the pipe at x, positive upwards [m]
 * @param {number} [o.g=9.80665]
 * @returns {Float64Array} pressure [Pa]
 */
export function pressureProfile({ x, pIn, frictionGradient = 0, rho = 0, z, g = 9.80665 }) {
  fin('pIn', pIn);
  const out = new Float64Array(x.length);
  for (let i = 0; i < x.length; i++) out[i] = pIn - frictionGradient * x[i] + (z ? rho * g * (z[0] - z[i]) : 0);
  return out;
}

function sampler(v, name) {
  if (typeof v === 'number') return () => v;
  if (typeof v === 'function') return v;
  throw new RangeError(`${name} must be a number or a function of x`);
}

/**
 * Effective axial force and end expansion of a pipeline that is free to slide at both ends against axial
 * soil friction.
 *
 * The fully restrained compressive force is P(x) = −S_0(x). Friction builds up from each end,
 * F_hot(x) = R_hot + ∫_0^x f dx and F_cold(x) = R_cold + ∫_x^L f dx. If the two slip zones do not meet (long
 * pipeline) there are two virtual anchors with a fully restrained length between; otherwise (short pipeline)
 * the whole line slides about a single virtual anchor where F_hot = F_cold. End expansion is the integral of
 * the released strain (P − C)/(E A_s) from the anchor to the end. End restraints may be forces (R) or linear
 * springs (spool / PLET flexibility), in which case R = k δ is solved by bisection.
 *
 * Reference: DNV-RP-F110 (2018) §6 (effective force build-up); Palmer & King, Subsea Pipeline Engineering,
 * ch. 12. Closed-form checks: anchor length P/f and end expansion P²/(2 f E A_s) for uniform P, long line.
 * Validity: first load (monotonic heat-up from the as-laid state), P(x) ≥ 0; in the restrained length the
 * gradient of P must not exceed f (flagged by `gradientExceedsFriction`).
 *
 * @param {object} o
 * @param {number} o.L pipeline length [m]
 * @param {number} o.As steel area [m²]
 * @param {number} o.Ai bore area [m²]
 * @param {number|function(number):number} o.friction axial friction per unit length f [N/m] (μ·w_s), may vary with x
 * @param {number|function(number):number} o.deltaT temperature rise above installation [K]
 * @param {number|function(number):number} [o.deltaPi=0] internal pressure rise above installation [Pa]
 * @param {number} [o.H=0] residual lay tension [N]
 * @param {number} [o.E=207e9] [Pa]
 * @param {number} [o.alpha=11.7e-6] [1/K]
 * @param {number} [o.nu=0.3] [-]
 * @param {number} [o.Rhot=0] resisting force at the hot end [N]
 * @param {number} [o.Rcold=0] resisting force at the cold end [N]
 * @param {number} [o.kHot] spring stiffness of the hot-end tie-in [N/m] (overrides Rhot)
 * @param {number} [o.kCold] spring stiffness of the cold-end tie-in [N/m] (overrides Rcold)
 * @param {number} [o.n=2000] number of intervals [-]
 * @returns {{x:Float64Array, S0:Float64Array, S:Float64Array, regime:'long'|'short', anchorHot:number,
 *   anchorCold:number, expansionHot:number, expansionCold:number, Rhot:number, Rcold:number,
 *   maxCompression:number, gradientExceedsFriction:boolean}}
 *   x [m]; S0, S [N] tension-positive; anchor positions [m] (equal for a short line); expansions [m] (outward
 *   positive); R* end reactions [N]; maxCompression [N]
 */
export function endExpansion(o) {
  const { L, As, Ai, H = 0, E = 207e9, alpha = 11.7e-6, nu = 0.3, kHot, kCold } = o;
  pos('L', L); pos('As', As); pos('Ai', Ai); pos('E', E);
  const n = Math.max(20, Math.floor(o.n === undefined ? 2000 : o.n));
  const fF = sampler(o.friction, 'friction');
  const fT = sampler(o.deltaT, 'deltaT');
  const fP = sampler(o.deltaPi === undefined ? 0 : o.deltaPi, 'deltaPi');
  const h = L / n;
  const x = new Float64Array(n + 1);
  const P = new Float64Array(n + 1);
  const cumF = new Float64Array(n + 1);
  let fPrev = 0;
  for (let i = 0; i <= n; i++) {
    x[i] = i * h;
    P[i] = -(H - fP(x[i]) * Ai * (1 - 2 * nu) - As * E * alpha * fT(x[i]));
    if (P[i] < -1e-9 * E * As) throw new RangeError('the restrained force must be compressive along the line (deltaT, deltaPi, H)');
    const fi = fF(x[i]);
    if (!(fi > 0)) throw new RangeError('friction must be positive along the line');
    if (i > 0) cumF[i] = cumF[i - 1] + 0.5 * (fi + fPrev) * h;
    fPrev = fi;
  }
  const EA = E * As;
  const C = new Float64Array(n + 1);

  const solve = (R0, RL) => {
    // first crossings of each friction line with the restrained force
    let ia = -1;
    for (let i = 0; i <= n; i++) if (R0 + cumF[i] >= P[i]) { ia = i; break; }
    let ib = -1;
    for (let i = n; i >= 0; i--) if (RL + (cumF[n] - cumF[i]) >= P[i]) { ib = i; break; }
    const long = ia >= 0 && ib >= 0 && ia <= ib;
    let xa;
    let xb;
    if (long) {
      for (let i = 0; i <= n; i++) C[i] = i < ia ? R0 + cumF[i] : i > ib ? RL + cumF[n] - cumF[i] : P[i];
      const cross = (i0, g0, g1) => x[i0] + (h * g0) / (g0 - g1);   // linear root of g between nodes i0 and i0+1
      xa = ia === 0 ? 0 : cross(ia - 1, R0 + cumF[ia - 1] - P[ia - 1], R0 + cumF[ia] - P[ia]);
      xb = ib === n ? L : cross(ib, RL + cumF[n] - cumF[ib] - P[ib], RL + cumF[n] - cumF[ib + 1] - P[ib + 1]);
    } else {
      // single virtual anchor where the two friction lines meet
      let k = 0;
      const g = (i) => R0 + cumF[i] - (RL + cumF[n] - cumF[i]);
      if (g(0) >= 0) xa = 0;
      else if (g(n) <= 0) xa = L;
      else {
        while (g(k + 1) < 0) k++;
        xa = x[k] + (h * g(k)) / (g(k) - g(k + 1));
      }
      xb = xa;
      for (let i = 0; i <= n; i++) C[i] = Math.min(R0 + cumF[i], RL + cumF[n] - cumF[i], P[i]);
    }
    // released strain integrated from each end to its anchor (trapezoid with the partial end cell)
    const integ = (from, to) => {
      let s = 0;
      for (let i = 0; i < n; i++) {
        const a = Math.max(x[i], from);
        const b = Math.min(x[i + 1], to);
        if (b <= a) continue;
        const e0 = (P[i] - C[i]) / EA;
        const e1 = (P[i + 1] - C[i + 1]) / EA;
        const ea = e0 + ((e1 - e0) * (a - x[i])) / h;
        const eb = e0 + ((e1 - e0) * (b - x[i])) / h;
        s += 0.5 * (ea + eb) * (b - a);
      }
      return s;
    };
    return { long, xa, xb, dHot: integ(0, xa), dCold: integ(xb, L) };
  };

  let R0 = kHot === undefined ? nonNeg('Rhot', o.Rhot === undefined ? 0 : o.Rhot) : 0;
  let RL = kCold === undefined ? nonNeg('Rcold', o.Rcold === undefined ? 0 : o.Rcold) : 0;
  let Pmax = 0;
  for (let i = 0; i <= n; i++) if (P[i] > Pmax) Pmax = P[i];
  let r = solve(R0, RL);
  if (kHot !== undefined || kCold !== undefined) {
    if (kHot !== undefined) pos('kHot', kHot);
    if (kCold !== undefined) pos('kCold', kCold);
    for (let outer = 0; outer < 60; outer++) {
      const prev0 = R0;
      const prevL = RL;
      if (kHot !== undefined) {
        let lo = 0;
        let hi = Pmax;
        for (let it = 0; it < 60; it++) {
          const mid = 0.5 * (lo + hi);
          if (solve(mid, RL).dHot * kHot > mid) lo = mid; else hi = mid;
        }
        R0 = 0.5 * (lo + hi);
      }
      if (kCold !== undefined) {
        let lo = 0;
        let hi = Pmax;
        for (let it = 0; it < 60; it++) {
          const mid = 0.5 * (lo + hi);
          if (solve(R0, mid).dCold * kCold > mid) lo = mid; else hi = mid;
        }
        RL = 0.5 * (lo + hi);
      }
      if (Math.abs(R0 - prev0) + Math.abs(RL - prevL) <= 1e-10 * (Pmax + 1)) break;
    }
    r = solve(R0, RL);
  }
  const S = new Float64Array(n + 1);
  const S0 = new Float64Array(n + 1);
  let maxC = 0;
  let steep = false;
  for (let i = 0; i <= n; i++) {
    S[i] = -C[i]; S0[i] = -P[i];
    if (C[i] > maxC) maxC = C[i];
    if (r.long && i > 0 && x[i] > r.xa && x[i] < r.xb && Math.abs(P[i] - P[i - 1]) > (cumF[i] - cumF[i - 1]) * (1 + 1e-9)) steep = true;
  }
  return { x, S0, S, regime: r.long ? 'long' : 'short', anchorHot: r.xa, anchorCold: r.xb, expansionHot: r.dHot, expansionCold: r.dCold,
    Rhot: R0, Rcold: RL, maxCompression: maxC, gradientExceedsFriction: steep };
}

/**
 * Closed-form anchor length and end expansion for a uniform restrained force P and uniform friction f.
 * Long line (L ≥ 2 L_a): L_a = (P − R)/f, δ = (P − R)²/(2 f E A_s).
 * Short line (L < 2 P/f, R = 0): anchor at mid-length, δ = (P L/2 − f L²/8)/(E A_s) per end.
 *
 * Reference: integration of the released strain (P − R − f x)/(E A_s); Palmer & King ch. 12.
 *
 * @param {object} o
 * @param {number} o.P fully restrained compressive force [N]
 * @param {number} o.f friction [N/m]
 * @param {number} o.EA axial stiffness E·A_s [N]
 * @param {number} [o.L=Infinity] pipeline length [m]
 * @param {number} [o.R=0] end resisting force [N] (long-line case only)
 * @returns {{anchorLength:number, expansion:number, regime:'long'|'short'}} [m]
 */
export function endExpansionUniform({ P, f, EA, L = Infinity, R = 0 }) {
  pos('P', P); pos('f', f); pos('EA', EA); nonNeg('R', R);
  const La = Math.max(0, (P - R) / f);
  if (L >= 2 * La) return { anchorLength: La, expansion: (P - R) ** 2 / (2 * f * EA) * (P > R ? 1 : 0), regime: 'long' };
  if (R !== 0) throw new RangeError('R is only supported for the long-line case');
  return { anchorLength: L / 2, expansion: ((P * L) / 2 - (f * L * L) / 8) / EA, regime: 'short' };
}

/**
 * Tie-in spool / expansion-loop leg by the guided-cantilever method: a leg of length L_leg perpendicular to
 * the pipeline absorbs the end displacement δ with both ends held against rotation.
 * F = 12 E I δ/L³, M = 6 E I δ/L², σ_b = M D/(2I) = 3 E D δ/L²; flexibility 1/k = L³/(12 E I);
 * leg length for an allowable stress: L = √(3 E D δ/σ_allow).
 *
 * Reference: guided cantilever method (M.W. Kellogg, Design of Piping Systems); elastic, no elbow
 * flexibility, no soil resistance on the leg (conservative for stress, unconservative for the reaction if the
 * leg is buried).
 *
 * @param {object} o
 * @param {number} o.delta imposed end displacement [m]
 * @param {number} o.D outer diameter [m]
 * @param {number} o.I second moment of area [m⁴]
 * @param {number} [o.Lleg] leg length [m]
 * @param {number} [o.sigmaAllow] allowable bending stress range [Pa] (gives the required leg)
 * @param {number} [o.E=207e9] [Pa]
 * @param {number} [o.sif=1] stress intensification factor [-]
 * @returns {{force:(number|null), moment:(number|null), stress:(number|null), stiffness:(number|null), requiredLeg:(number|null)}}
 *   force [N]; moment [N·m]; stress [Pa]; stiffness [N/m]; requiredLeg [m]
 */
export function spoolGuidedCantilever({ delta, D, I, Lleg, sigmaAllow, E = 207e9, sif = 1 }) {
  nonNeg('delta', delta); pos('D', D); pos('I', I); pos('sif', sif);
  const out = { force: null, moment: null, stress: null, stiffness: null, requiredLeg: null };
  if (Lleg !== undefined) {
    pos('Lleg', Lleg);
    out.stiffness = (12 * E * I) / Lleg ** 3;
    out.force = out.stiffness * delta;
    out.moment = (6 * E * I * delta) / (Lleg * Lleg);
    out.stress = (sif * out.moment * D) / (2 * I);
  }
  if (sigmaAllow !== undefined) out.requiredLeg = Math.sqrt((3 * E * D * delta * sif) / pos('sigmaAllow', sigmaAllow));
  return out;
}

/**
 * Walking per thermal cycle of a short (fully mobilised) pipeline on a seabed slope:
 * Δ = [ |ΔP| + W L |sin φ| − μ W L cos φ ] · L tan φ/(E A μ), down the slope.
 *
 * Derivation (reproduced in the tests against the numerical model): the heat-up and cool-down virtual
 * anchors lie L tan φ/μ apart and the effective force change between them is uniform.
 * Reference: Carr, Sinclair & Bruton (2006), "Pipeline Walking — Understanding the Field Layout
 * Challenges and Analytical Solutions Developed for the SAFEBUCK JIP", OTC 17945.
 * Validity: fully mobilised over the whole length in both heat-up and cool-down (the bracket is positive,
 * i.e. |ΔP| > μ W L cos φ − W L |sin φ|) and tan φ < μ; returns 0 when the bracket is not positive.
 *
 * @param {object} o
 * @param {number} o.deltaP change in fully constrained force between cold and hot conditions [N]
 * @param {number} o.W submerged weight [N/m]
 * @param {number} o.L length [m]
 * @param {number} o.mu axial friction coefficient [-]
 * @param {number} o.slopeRad seabed slope φ [rad]
 * @param {number} o.EA axial stiffness [N]
 * @returns {{walkPerCycle:number, anchorSeparation:number, fullyMobilised:boolean}} [m]
 */
export function walkingSlope({ deltaP, W, L, mu, slopeRad, EA }) {
  fin('deltaP', deltaP); pos('W', W); pos('L', L); pos('mu', mu); fin('slopeRad', slopeRad); pos('EA', EA);
  const phi = Math.abs(slopeRad);
  if (Math.tan(phi) >= mu) throw new RangeError('slope exceeds the friction angle: the pipe slides without cycling');
  const bracket = Math.abs(deltaP) + W * L * Math.sin(phi) - mu * W * L * Math.cos(phi);
  const sep = (L * Math.tan(phi)) / mu;
  return { walkPerCycle: bracket > 0 ? (bracket * sep) / EA : 0, anchorSeparation: sep, fullyMobilised: bracket > 0 };
}

/**
 * Walking per thermal cycle of a short pipeline with a constant tension S_R applied at one end (e.g. a
 * steel catenary riser): Δ = ( |ΔP| + S_R − μ W L ) · S_R/(E A μ W), towards the tensioned end.
 *
 * Reference: Carr, Sinclair & Bruton (2006), OTC 17945 (the anchors lie S_R/(μ W) apart).
 * Validity: fully mobilised (positive bracket) and S_R < μ W L; returns 0 when the bracket is not positive.
 *
 * @param {object} o
 * @param {number} o.deltaP change in fully constrained force between cold and hot [N]
 * @param {number} o.SR end tension [N]
 * @param {number} o.W submerged weight [N/m]
 * @param {number} o.L [m]
 * @param {number} o.mu [-]
 * @param {number} o.EA [N]
 * @returns {{walkPerCycle:number, anchorSeparation:number, fullyMobilised:boolean}} [m]
 */
export function walkingEndTension({ deltaP, SR, W, L, mu, EA }) {
  fin('deltaP', deltaP); nonNeg('SR', SR); pos('W', W); pos('L', L); pos('mu', mu); pos('EA', EA);
  const f = mu * W;
  if (SR >= f * L) throw new RangeError('SR exceeds the total friction capacity: the pipe slides without cycling');
  const bracket = Math.abs(deltaP) + SR - f * L;
  return { walkPerCycle: bracket > 0 ? (bracket * SR) / (EA * f) : 0, anchorSeparation: SR / f, fullyMobilised: bracket > 0 };
}

/**
 * Constant hold-back tension that arrests slope- and tension-driven walking of a fully mobilised pipeline,
 * and the load range seen by a rigid anchor at that end.
 * A constant tension T applied at the upper end moves the heat-up and cool-down anchors to coincide when
 * T = W L sin φ + S_R (derived from the same force profiles as walkingSlope / walkingEndTension).
 * A rigid end anchor instead sees, each cycle, compression up to (μ cos φ − sin φ) W L − S_R on heat-up and
 * tension up to (μ cos φ + sin φ) W L + S_R on cool-down (limited by |ΔP| for partly mobilised lines).
 *
 * Reference: first-principles equilibrium of the fully mobilised friction profiles (cf. Carr et al. 2006).
 *
 * @param {object} o
 * @param {number} o.W [N/m]
 * @param {number} o.L [m]
 * @param {number} o.mu [-]
 * @param {number} [o.slopeRad=0] [rad]
 * @param {number} [o.SR=0] tension at the lower / far end [N]
 * @returns {{holdBackTension:number, rigidAnchorTension:number, rigidAnchorCompression:number}} [N]
 */
export function walkingArrestForce({ W, L, mu, slopeRad = 0, SR = 0 }) {
  pos('W', W); pos('L', L); pos('mu', mu); nonNeg('SR', SR);
  const phi = Math.abs(slopeRad);
  return {
    holdBackTension: W * L * Math.sin(phi) + SR,
    rigidAnchorTension: (mu * Math.cos(phi) + Math.sin(phi)) * W * L + SR,
    rigidAnchorCompression: Math.max(0, (mu * Math.cos(phi) - Math.sin(phi)) * W * L - SR),
  };
}

/**
 * Load steps for one heat-up / cool-down cycle with a thermal front travelling from the hot end.
 * Heat-up step j (1..nHeat): T = T_amb + (T_in − T_amb) e^(−x/λ) behind the front at x_f = L j/nHeat and
 * ambient ahead of it; pressure strain is applied in full with the first step. Cool-down: the whole line
 * returns to ambient uniformly in nCool equal decrements (as assumed by Carr et al. 2006).
 *
 * Reference: idealised start-up transient of Carr, Sinclair & Bruton (2006), OTC 17945.
 *
 * @param {object} o
 * @param {number} o.L [m]
 * @param {number} o.nElem number of elements [-]
 * @param {number} o.alphaDeltaT thermal strain at the inlet α (T_in − T_amb) [-]
 * @param {number} [o.lambda=Infinity] decay length of the steady profile [m]
 * @param {number} [o.pressureStrain=0] free axial strain from pressure Δp A_i (1 − 2ν)/(E A_s) [-]
 * @param {number} [o.nHeat=20] number of front positions [-]
 * @param {number} [o.nCool=4] [-]
 * @returns {Float64Array[]} free (unrestrained) axial strain per element for each load step
 */
export function thermalTransientSteps({ L, nElem, alphaDeltaT, lambda = Infinity, pressureStrain = 0, nHeat = 20, nCool = 4 }) {
  pos('L', L); pos('alphaDeltaT', alphaDeltaT);
  const n = Math.floor(pos('nElem', nElem));
  const h = L / n;
  const steady = new Float64Array(n);
  for (let e = 0; e < n; e++) steady[e] = alphaDeltaT * Math.exp(-((e + 0.5) * h) / lambda);
  const steps = [];
  for (let j = 1; j <= nHeat; j++) {
    const xf = (L * j) / nHeat;
    const s = new Float64Array(n);
    for (let e = 0; e < n; e++) {
      const a = e * h;
      const frac = Math.max(0, Math.min(1, (xf - a) / h));   // part of the element behind the front
      s[e] = pressureStrain + steady[e] * frac;
    }
    steps.push(s);
  }
  for (let j = nCool - 1; j >= 0; j--) {
    const s = new Float64Array(n);
    for (let e = 0; e < n; e++) s[e] = ((pressureStrain + steady[e]) * j) / nCool;
    steps.push(s);
  }
  return steps;
}

/**
 * Load steps for a cycle of uniform heating and cooling (no thermal gradient): free strain ramps to
 * `freeStrain` in nUp steps and back to zero in nDown steps.
 *
 * @param {object} o
 * @param {number} o.nElem [-]
 * @param {number} o.freeStrain unrestrained axial strain when hot [-]
 * @param {number} [o.nUp=4] [-]
 * @param {number} [o.nDown=4] [-]
 * @returns {Float64Array[]} free strain per element for each step
 */
export function uniformCycleSteps({ nElem, freeStrain, nUp = 4, nDown = 4 }) {
  const n = Math.floor(pos('nElem', nElem));
  const steps = [];
  for (let j = 1; j <= nUp; j++) steps.push(new Float64Array(n).fill((freeStrain * j) / nUp));
  for (let j = nDown - 1; j >= 0; j--) steps.push(new Float64Array(n).fill((freeStrain * j) / nDown));
  return steps;
}

/**
 * Numerical axial walking model: an elastic bar (axial stiffness EA) on a Coulomb-friction seabed, loaded
 * by a sequence of free-strain fields (thermal + pressure) and cycled. Friction is elastic–perfectly
 * plastic with mobilisation displacement uMob; the seabed slope adds a distributed force W sin φ along +x
 * (downhill) and scales the friction limit to μ W cos φ; a constant end tension and a linear anchor spring
 * may be applied. Each load step is a convex minimisation, solved on the nodal equilibrium equations
 * (tridiagonal system, return mapping of the friction sliders) by Newton steps that are accepted only when
 * they lower the potential energy, with a secant-stiffness step as the always-descending fallback; a step
 * that still fails to converge is halved.
 *
 * Walking emerges from first principles: with symmetric loading, no slope, no end tension and no thermal
 * gradient the net displacement per cycle is zero; a slope, end tension or travelling thermal front gives a
 * steady axial shift per cycle. The walk per cycle is the change of the mean pipe displacement between
 * successive cold states.
 *
 * Reference: 1-D axial pipe–soil model as in Carr, Sinclair & Bruton (2006), OTC 17945 (their analytical
 * results are limiting cases); Tørnes et al. (2000), OMAE PIPE-5038.
 * Validity: straight pipe, no lateral buckling, linear elastic steel, rate-independent friction.
 *
 * @param {object} o
 * @param {number} o.L length [m]
 * @param {number} o.EA axial stiffness [N]
 * @param {number} o.W submerged weight [N/m]
 * @param {number} o.mu axial friction coefficient [-]
 * @param {Float64Array[]} o.steps free-strain field per element for each load step of one cycle (the last
 *   step should return to the cold state); see uniformCycleSteps / thermalTransientSteps
 * @param {number} [o.nElem=steps[0].length] number of elements [-]
 * @param {number} [o.cycles=5] number of cycles [-]
 * @param {number} [o.slopeRad=0] seabed slope, positive when the seabed falls towards +x [rad]
 * @param {number} [o.tensionAtEnd=0] constant tension applied at x = L (pulling towards +x) [N]
 * @param {number} [o.tensionAtStart=0] constant tension applied at x = 0 (pulling towards −x) [N]
 * @param {number} [o.uMob=0.002] friction mobilisation displacement [m]
 * @param {{x:number, stiffness:number}} [o.anchor] linear anchor spring to ground at position x [m], [N/m]
 * @param {number} [o.tol=1e-9] equilibrium tolerance relative to μ W L [-]
 * @returns {{x:Float64Array, u:Float64Array, force:Float64Array, meanDisplacement:number[], walkPerCycle:number[],
 *   walkRate:number, endDisplacementStart:number[], endDisplacementEnd:number[], hotForce:Float64Array,
 *   coldForce:Float64Array, anchorForceMin:number, anchorForceMax:number, iterations:number}}
 *   x node positions [m]; u final (cold) displacements [m]; force element axial force at the end, tension
 *   positive [N]; meanDisplacement after each cycle [m]; walkPerCycle [m] (positive towards +x); walkRate
 *   last-cycle value [m/cycle]; hotForce / coldForce element forces at the most compressive step and at the
 *   end of the last cycle [N]; anchorForce* extreme forces in the anchor spring over the last cycle [N]
 *   (positive when the anchor pulls the pipe towards −x)
 */
export function axialWalkingModel(o) {
  const { L, EA, W, mu, steps, cycles = 5, slopeRad = 0, tensionAtEnd = 0, tensionAtStart = 0, uMob = 0.002, anchor, tol = 1e-9 } = o;
  pos('L', L); pos('EA', EA); pos('W', W); pos('mu', mu); pos('uMob', uMob); fin('slopeRad', slopeRad);
  if (!Array.isArray(steps) || steps.length < 2) throw new RangeError('steps must hold at least two load steps');
  const n = Math.floor(o.nElem === undefined ? steps[0].length : o.nElem);
  if (!(n >= 4)) throw new RangeError('nElem must be at least 4');
  for (const s of steps) if (s.length !== n) throw new RangeError('every load step must have nElem entries');
  if (!(cycles >= 1)) throw new RangeError('cycles must be at least 1');
  const h = L / n;
  const kb = EA / h;
  const fLim = mu * W * Math.cos(slopeRad);
  const grav = W * Math.sin(slopeRad);
  const nn = n + 1;
  const trib = new Float64Array(nn);
  for (let i = 0; i < nn; i++) trib[i] = i === 0 || i === n ? h / 2 : h;
  const u = new Float64Array(nn);
  const up = new Float64Array(nn);          // plastic slip of each friction slider
  const q = new Float64Array(nn);           // friction force on each node (resisting, acts −sign(u − up))
  const slid = new Uint8Array(nn);
  const res = new Float64Array(nn);
  const lower = new Float64Array(nn);
  const diag = new Float64Array(nn);
  const diagSec = new Float64Array(nn);
  const upper = new Float64Array(nn);
  const eps = new Float64Array(n);
  const epsFrom = new Float64Array(n);
  const epsTo = new Float64Array(n);
  let iAnchor = -1;
  let kAnchor = 0;
  if (anchor) {
    kAnchor = pos('anchor.stiffness', anchor.stiffness);
    iAnchor = Math.max(0, Math.min(n, Math.round(anchor.x / h)));
  }
  const scale = fLim * L;
  let iterations = 0;

  const assemble = () => {
    let norm = 0;
    for (let i = 0; i < nn; i++) {
      const ks = (fLim * trib[i]) / uMob;
      let qi = ks * (u[i] - up[i]);
      const lim = fLim * trib[i];
      if (qi > lim) { qi = lim; slid[i] = 1; } else if (qi < -lim) { qi = -lim; slid[i] = 1; } else slid[i] = 0;
      q[i] = qi;
      let r = grav * trib[i] - qi;
      if (i < n) r += kb * (u[i + 1] - u[i]) - EA * eps[i];
      if (i > 0) r -= kb * (u[i] - u[i - 1]) - EA * eps[i - 1];
      if (i === n) r += tensionAtEnd;
      if (i === 0) r -= tensionAtStart;
      if (i === iAnchor) r -= kAnchor * u[i];
      res[i] = r;
      norm = Math.max(norm, Math.abs(r));
      // two tangents: the consistent one (sliding nodes carry no stiffness) for fast Newton steps, and the
      // secant one, whose quadratic model majorises the convex energy so that its step always lowers it
      const base = (i < n ? kb : 0) + (i > 0 ? kb : 0) + (i === iAnchor ? kAnchor : 0);
      diag[i] = base + (slid[i] ? ks * 1e-9 : ks);
      diagSec[i] = base + (slid[i] ? Math.abs(qi) / Math.abs(u[i] - up[i]) : ks);
      lower[i] = -kb; upper[i] = -kb;
    }
    return norm;
  };

  // total potential energy (convex): bar strain energy + friction potential − work of the external loads
  const energy = () => {
    let e = 0;
    for (let i = 0; i < nn; i++) {
      const lim = fLim * trib[i];
      const ks = lim / uMob;
      const d = u[i] - up[i];
      e += Math.abs(ks * d) <= lim ? 0.5 * ks * d * d : lim * Math.abs(d) - (lim * lim) / (2 * ks);
      e -= grav * trib[i] * u[i];
      if (i < n) { const s = u[i + 1] - u[i] - eps[i] * h; e += 0.5 * kb * s * s; }
    }
    e -= tensionAtEnd * u[n] - tensionAtStart * u[0];
    if (iAnchor >= 0) e += 0.5 * kAnchor * u[iAnchor] * u[iAnchor];
    return e;
  };

  const uSave = new Float64Array(nn);
  const u0 = new Float64Array(nn);
  const maxIter = 400 + 20 * nn;
  const equilibrate = () => {
    uSave.set(u);
    for (let it = 0; it < maxIter; it++) {
      const norm = assemble();
      if (norm <= tol * scale) {
        for (let i = 0; i < nn; i++) up[i] = u[i] - (q[i] * uMob) / (fLim * trib[i]);
        return true;
      }
      iterations++;
      // Newton step, kept only if it lowers the energy; otherwise the (always descending) secant step
      const e0 = energy();
      u0.set(u);
      const duN = solveTridiagonal(lower, diag, upper, res);
      for (let i = 0; i < nn; i++) u[i] = u0[i] + duN[i];
      if (!(energy() < e0)) {
        const duS = solveTridiagonal(lower, diagSec, upper, res);
        for (let i = 0; i < nn; i++) u[i] = u0[i] + duS[i];
      }
    }
    u.set(uSave);
    return false;
  };

  const advance = (from, to, depth) => {
    eps.set(to);
    if (equilibrate()) return;
    if (depth > 12) throw new Error('axial walking model did not converge; reduce the load-step size or check that friction can hold the pipe');
    const mid = new Float64Array(n);
    for (let e = 0; e < n; e++) mid[e] = 0.5 * (from[e] + to[e]);
    advance(from, mid, depth + 1);
    advance(mid, to, depth + 1);
  };

  const mean = () => {
    let s = 0;
    for (let i = 0; i < nn; i++) s += u[i] * trib[i];
    return s / L;
  };
  const forces = () => {
    const f = new Float64Array(n);
    for (let e = 0; e < n; e++) f[e] = kb * (u[e + 1] - u[e]) - EA * eps[e];
    return f;
  };

  // as-laid state under the permanent loads (slope, end tensions)
  advance(epsFrom, epsFrom, 0);
  const meanDisplacement = [];
  const walkPerCycle = [];
  const endStart = [];
  const endEnd = [];
  let prevMean = mean();
  let hotForce = forces();
  let anchorMin = 0;
  let anchorMax = 0;
  for (let c = 0; c < cycles; c++) {
    let minSum = Infinity;
    if (c === cycles - 1 && iAnchor >= 0) { anchorMin = Infinity; anchorMax = -Infinity; }
    for (let s = 0; s < steps.length; s++) {
      epsTo.set(steps[s]);
      advance(Float64Array.from(epsFrom), Float64Array.from(epsTo), 0);
      epsFrom.set(epsTo);
      if (c === cycles - 1) {
        const f = forces();
        let sum = 0;
        for (let e = 0; e < n; e++) sum += f[e];
        if (sum < minSum) { minSum = sum; hotForce = f; }
        if (iAnchor >= 0) {
          const fa = kAnchor * u[iAnchor];
          anchorMin = Math.min(anchorMin, fa); anchorMax = Math.max(anchorMax, fa);
        }
      }
    }
    const m = mean();
    meanDisplacement.push(m);
    walkPerCycle.push(m - prevMean);
    prevMean = m;
    endStart.push(u[0]);
    endEnd.push(u[n]);
  }
  const x = new Float64Array(nn);
  for (let i = 0; i < nn; i++) x[i] = i * h;
  const coldForce = forces();
  return { x, u: Float64Array.from(u), force: coldForce, meanDisplacement, walkPerCycle, walkRate: walkPerCycle[walkPerCycle.length - 1],
    endDisplacementStart: endStart, endDisplacementEnd: endEnd, hotForce, coldForce, anchorForceMin: anchorMin, anchorForceMax: anchorMax, iterations };
}
