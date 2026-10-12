/**
 * One-dimensional transient isothermal gas–liquid drift-flux solver along an elevation profile (finite
 * volume, explicit), scenario helpers (ramp-up, shutdown/restart, riser-induced severe slugging and its
 * mitigation), surge detection, and a simple Lagrangian slug tracker.
 *
 * Model (Evje & Fjelde 2002, J. Comput. Phys. 175, 674–701, Eqs. 1–3):
 *   ∂t(αL ρL) + ∂s(αL ρL vL) = 0
 *   ∂t(αG ρG) + ∂s(αG ρG vG) = qG
 *   ∂t(αL ρL vL + αG ρG vG) + ∂s(αL ρL vL² + αG ρG vG² + p) = −f ρm vm|vm|/(2D) − ρm g sinθ
 * with ρL = ρL0 + (p − p0)/aL² (weakly compressible liquid), ρG = p/(Z R T/M) (isothermal gas) and the slip
 * law vG = C0 vm + vd, vm = αL vL + αG vG. s is the chainage along the pipe, θ the local inclination
 * (positive upward in the +s direction).
 *
 * Numerical flux: hybrid flux-vector splitting of the AUSM family introduced for this model by Evje & Fjelde
 * (2003), Computers & Fluids 32, 1497–1530 ("AUSMV"/"AUSMD"). Phase mass and convective momentum fluxes are
 * split with the velocity polynomials V± (Liou & Steffen 1993; Wada & Liou 1997) using the contact weights
 * χL = 2 mR/(mL + mR), χR = 2 mL/(mL + mR) per phase (m = α ρ), which make the mass flux pure upwinding
 * whenever the two sides share a velocity, so volume-fraction fronts are not smeared at acoustic speed. These
 * weights follow the Wada–Liou construction; they are stated here rather than transcribed from the 2003 paper.
 * The interface pressure is the linearised (acoustic) Riemann value with impedances z = ρ c,
 *   p* = [zR pL + zL pR + zL zR (wL − wR)]/(zL + zR)      (Toro 2009, Riemann Solvers, §9.3),
 * blended into the upwinding P± pressure splitting for mixture Mach numbers between 0.5 and 1. The plain
 * average of the P± splitting is unstable at a liquid/gas interface at rest under gravity (impedance ratio
 * ~10³); the impedance weighting removes that. The matching pressure-diffusion velocity 2 Kp (pL − pR)/(zL + zR)
 * is added to the mass flux (Kp = ½ is the acoustic Riemann value; cf. Liou 2006, J. Comput. Phys. 214,
 * 137–170, AUSM⁺-up) and damps odd–even pressure modes at low Mach number.
 * Space: MUSCL reconstruction of (αG, p, vL) with a slope limiter; each cell extrapolates its pressure to the
 * faces along its own hydrostatic profile (hydrostatic reconstruction, cf. Audusse et al. 2004, SIAM J. Sci.
 * Comput. 25, 2050–2065), so a fluid at rest in discrete hydrostatic balance is preserved. Time: two-stage SSP Runge–Kutta (Shu & Osher 1988) under a CFL limit on |v| + c,
 * where c is the mixture sound speed — the time step is acoustic (liquid-filled cells: Δt ≈ CFL Δs/aL).
 */
import { G, checkArg, darcyFrictionFactorFast } from './driftflux.js';
import { boeCriterion, taitelStability } from './slug.js';

const P_FLOOR = 1;

function asFunction(v, name) {
  if (typeof v === 'function') return v;
  if (typeof v === 'number') { checkArg(name, v); return () => v; }
  if (v && (Array.isArray(v.t) || ArrayBuffer.isView(v.t)) && v.v && v.v.length === v.t.length && v.t.length > 0) {
    const ts = v.t;
    const vs = v.v;
    return (t) => {
      if (t <= ts[0]) return vs[0];
      const n = ts.length;
      if (t >= ts[n - 1]) return vs[n - 1];
      let i = 1;
      while (ts[i] < t) i++;
      const w = (t - ts[i - 1]) / (ts[i] - ts[i - 1]);
      return vs[i - 1] * (1 - w) + vs[i] * w;
    };
  }
  throw new RangeError(`${name} must be a number, a function of time or a series { t: [], v: [] }`);
}

function limiterFunction(name) {
  switch (name) {
    case 'first': return () => 0;
    case 'none': return (a, b) => 0.5 * (a + b);
    case 'minmod': return (a, b) => (a * b <= 0 ? 0 : (Math.abs(a) < Math.abs(b) ? a : b));
    case 'vanleer': return (a, b) => (a * b <= 0 ? 0 : 2 * a * b / (a + b));
    case 'mc': return (a, b) => {
      if (a * b <= 0) return 0;
      const s = a > 0 ? 1 : -1;
      return s * Math.min(2 * Math.abs(a), 2 * Math.abs(b), 0.5 * Math.abs(a + b));
    };
    default: throw new RangeError(`limiter = ${name} is not one of first, minmod, vanleer, mc, none`);
  }
}

/**
 * Uniform finite-volume grid along a route.
 * @param {{x: number, z: number}[]} profile route points: horizontal distance x [m], elevation z [m, up]
 * @param {number} nCells number of cells [-]
 * @returns {{n: number, ds: number, length: number, s: Float64Array, x: Float64Array, z: Float64Array,
 *   sinTheta: Float64Array, theta: Float64Array}} cell size ds and total length [m]; cell-centre chainage s,
 *   position x, z [m]; inclination of the route segment holding the cell centre [rad] and its sine [-]
 */
export function makeGrid(profile, nCells) {
  if (!Array.isArray(profile) || profile.length < 2) throw new RangeError('profile must list at least two {x, z} points');
  checkArg('nCells', nCells, { min: 2, max: 100000 });
  const n = Math.round(nCells);
  const segs = [];
  let length = 0;
  for (let i = 1; i < profile.length; i++) {
    const dx = profile[i].x - profile[i - 1].x;
    const dz = profile[i].z - profile[i - 1].z;
    const len = Math.hypot(dx, dz);
    if (!(len > 0)) throw new RangeError(`profile[${i}] coincides with the previous point`);
    segs.push({ s0: length, len, dx, dz, x0: profile[i - 1].x, z0: profile[i - 1].z });
    length += len;
  }
  const ds = length / n;
  const s = new Float64Array(n);
  const x = new Float64Array(n);
  const z = new Float64Array(n);
  const sinTheta = new Float64Array(n);
  const theta = new Float64Array(n);
  let k = 0;
  for (let j = 0; j < n; j++) {
    s[j] = (j + 0.5) * ds;
    while (k < segs.length - 1 && s[j] > segs[k].s0 + segs[k].len) k++;
    const g = segs[k];
    const w = (s[j] - g.s0) / g.len;
    x[j] = g.x0 + w * g.dx;
    z[j] = g.z0 + w * g.dz;
    sinTheta[j] = g.dz / g.len;
    theta[j] = Math.asin(Math.max(-1, Math.min(1, sinTheta[j])));
  }
  return { n, ds, length, s, x, z, sinTheta, theta };
}

/**
 * Route of a straight flowline followed by a straight riser.
 * @param {object} o
 * @param {number} o.flowlineLength flowline length along the pipe [m]
 * @param {number} [o.flowlineAngleDeg=-2] flowline inclination [deg] (negative = downward towards the riser)
 * @param {number} o.riserHeight riser height [m]
 * @param {number} [o.riserAngleDeg=90] riser inclination [deg]
 * @returns {{x: number, z: number}[]} three route points [m]
 */
export function flowlineRiserProfile({ flowlineLength, flowlineAngleDeg = -2, riserHeight, riserAngleDeg = 90 }) {
  checkArg('flowlineLength', flowlineLength, { min: 0, exclusiveMin: true });
  checkArg('riserHeight', riserHeight, { min: 0, exclusiveMin: true });
  checkArg('flowlineAngleDeg', flowlineAngleDeg, { min: -89, max: 89 });
  checkArg('riserAngleDeg', riserAngleDeg, { min: 1, max: 90 });
  const a = flowlineAngleDeg * Math.PI / 180;
  const r = riserAngleDeg * Math.PI / 180;
  const x1 = flowlineLength * Math.cos(a);
  const z1 = flowlineLength * Math.sin(a);
  return [{ x: 0, z: 0 }, { x: x1, z: z1 }, { x: x1 + (riserAngleDeg === 90 ? 0 : riserHeight / Math.tan(r)), z: z1 + riserHeight }];
}

/**
 * Transient drift-flux simulation.
 *
 * Boundary conditions (the model has two inflow and one outflow characteristic at a subsonic inlet, and the
 * reverse at the outlet, so the inlet takes two conditions and the outlet one):
 *   inlet   { type: 'massflow', massRateL, massRateG }  phase mass rates [kg/s]: number, function of t, or series
 *           { type: 'closed' }                           shut-in valve
 *   outlet  { type: 'pressure', pressure, backflowGasFraction = 1 }   pressure [Pa]: number, function or series
 *           { type: 'closed' }
 *           { type: 'choke', pressure, opening }  both phases leave at v = opening √(2 (p − pdown)/ρm) referred
 *             to the pipe area; opening = Cd·(choke area)/(pipe area) (number, function or series, 0 … 1) and
 *             `pressure` is the pressure downstream of the choke. Sub-critical, incompressible orifice equation
 *             on the mixture density — a simple documented choke, not a critical-flow model.
 *   both    { type: 'periodic' } (verification only)
 * Optional gas-lift sources: sources = [{ s, massRateG }] inject gas mass [kg/s] in the cell containing s.
 *
 * Initial conditions:
 *   { type: 'uniform', P, alphaG, vL = 0 }
 *   { type: 'profile', P, alphaG, vL }   arrays per cell or functions of chainage s
 *   { type: 'steady', massRateL, massRateG, pOut }   steady drift-flux profile marched upstream from the outlet
 *   { type: 'settled', liquidFraction, pRef, refIndex = last }   shut-in state: liquid at rest in the lowest part
 *             of the route (common level), gas above, in discrete hydrostatic balance
 *
 * @param {object} o
 * @param {{x: number, z: number}[]} o.profile route points [m]
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.roughness=0] wall roughness [m]
 * @param {number} o.nCells number of cells [-]
 * @param {number} o.tEnd simulated time [s]
 * @param {object} o.fluid
 * @param {number} o.fluid.rhoL liquid density at p0 [kg/m³]
 * @param {number} [o.fluid.p0=101325] reference pressure of rhoL [Pa]
 * @param {number} [o.fluid.aL=1000] liquid sound speed [m/s]; lower values (≥ ~10× the flow velocity) relax the
 *   acoustic time-step limit at the cost of a more compressible liquid
 * @param {number} o.fluid.gasRT Z R T/M of the gas [J/kg] (ρG = p/gasRT), e.g. 287·293 for air at 20 °C
 * @param {number} [o.fluid.muL=1e-3] liquid viscosity [Pa·s]
 * @param {number} [o.fluid.muG=1.8e-5] gas viscosity [Pa·s]
 * @param {object} [o.slip]
 * @param {number} [o.slip.C0=1.2] distribution coefficient [-] (1 ≤ C0 < 1/alphaBlend)
 * @param {number|'vertical'|'bendiksen'|((theta: number) => number)} [o.slip.vd='vertical'] drift velocity
 *   [m/s] in the +s direction: number; 'vertical' = 0.35 √(gD) sinθ; 'bendiksen' = (0.35 sinθ + 0.54 cosθ) √(gD)
 *   (low-Froude branch of Bendiksen 1984); or a function of the cell inclination [rad]
 * @param {number} [o.slip.alphaBlend=0.75] gas fraction above which C0 − 1 and vd are scaled by
 *   [(1 − αG)/(1 − alphaBlend)]², so the slip velocity vanishes and the slip law stays non-singular as αG → 1
 * @param {boolean} [o.friction=true] include wall friction
 * @param {object} o.inlet inlet boundary (above)
 * @param {object} o.outlet outlet boundary (above)
 * @param {{s: number, massRateG: number|Function|object}[]} [o.sources] gas-lift sources
 * @param {object} o.initial initial condition (above)
 * @param {'first'|'minmod'|'vanleer'|'mc'|'none'} [o.limiter='vanleer'] MUSCL slope limiter ('first' = first order)
 * @param {number} [o.cfl=0.45] Courant number on |v| + c [-], 0 < cfl ≤ 0.9
 * @param {number} [o.Kp=0.25] pressure-diffusion coefficient of the mass flux [-] (0 switches it off and
 *   leaves odd–even pressure modes undamped; 0.5 is the acoustic Riemann value)
 * @param {number} [o.nProfiles=50] number of stored profiles (plus the initial one) [-]
 * @param {number} [o.nHistory=1000] number of history intervals [-]
 * @param {number} [o.maxSteps=5e6] abort threshold on the number of time steps [-]
 * @returns {{grid: object, profiles: {t: number[], P: Float64Array[], holdup: Float64Array[], vL: Float64Array[],
 *   vG: Float64Array[]}, history: {t: Float64Array, outletLiquidRate: Float64Array, outletGasRate: Float64Array,
 *   inletPressure: Float64Array, outletPressure: Float64Array, liquidInventory: Float64Array,
 *   maxPressure: Float64Array}, final: {P: Float64Array, alphaG: Float64Array, holdup: Float64Array,
 *   vL: Float64Array, vG: Float64Array}, conservation: {massL0: number, massG0: number, massL: number,
 *   massG: number, netInflowL: number, netInflowG: number, residualL: number, residualG: number,
 *   clippedL: number, clippedG: number}, steps: number, time: number, dtMin: number, dtMax: number}}
 *   profiles at the stored times t [s]: P [Pa], liquid holdup [-], phase velocities [m/s] per cell; history on
 *   nHistory intervals: mid-interval time [s], interval-mean outlet mass rates [kg/s] (positive out of the
 *   pipe), end-of-interval first/last-cell pressures and maximum pressure [Pa] and liquid inventory [m³];
 *   conservation: phase masses in the pipe [kg] at start and end, net boundary + source inflow [kg], relative
 *   residual (mass − mass0 − netInflow)/max(mass0, |netInflow|) [-] (round-off for a correct run), and any
 *   mass removed by the positivity clip [kg] (0 in a healthy run)
 */
export function simulateTransient(o) {
  const { profile, D, roughness = 0, nCells, tEnd, fluid, slip = {}, friction = true, inlet, outlet, sources = [], initial,
    limiter = 'vanleer', cfl = 0.45, Kp = 0.25, nProfiles = 50, nHistory = 1000, maxSteps = 5e6 } = o;
  checkArg('D', D, { min: 0, exclusiveMin: true });
  checkArg('roughness', roughness, { min: 0 });
  checkArg('tEnd', tEnd, { min: 0, exclusiveMin: true });
  checkArg('cfl', cfl, { min: 0, max: 0.9, exclusiveMin: true });
  checkArg('Kp', Kp, { min: 0, max: 1 });
  checkArg('nProfiles', nProfiles, { min: 1, max: 2000 });
  checkArg('nHistory', nHistory, { min: 1, max: 200000 });
  if (!fluid) throw new RangeError('fluid is required');
  const rhoL0 = checkArg('fluid.rhoL', fluid.rhoL, { min: 0, exclusiveMin: true });
  const p0 = checkArg('fluid.p0', fluid.p0 ?? 101325, { min: 0 });
  const aL = checkArg('fluid.aL', fluid.aL ?? 1000, { min: 0, exclusiveMin: true });
  const aG2 = checkArg('fluid.gasRT', fluid.gasRT, { min: 0, exclusiveMin: true });
  const muL = checkArg('fluid.muL', fluid.muL ?? 1e-3, { min: 0, exclusiveMin: true });
  const muG = checkArg('fluid.muG', fluid.muG ?? 1.8e-5, { min: 0, exclusiveMin: true });
  const aL2 = aL * aL;
  const B = rhoL0 * aL2 - p0;
  if (!(B > 0)) throw new RangeError('fluid.aL is too low: rhoL·aL² must exceed p0');
  const C00 = checkArg('slip.C0', slip.C0 ?? 1.2, { min: 1, max: 2 });
  const ab = checkArg('slip.alphaBlend', slip.alphaBlend ?? 0.75, { min: 0, max: 0.999 });
  if (C00 * ab >= 1) throw new RangeError('slip.C0 · slip.alphaBlend must be below 1');
  const lim = limiterFunction(limiter);
  if (!inlet || !outlet || !initial) throw new RangeError('inlet, outlet and initial are required');

  const grid = makeGrid(profile, nCells);
  const N = grid.n;
  const ds = grid.ds;
  const A = Math.PI * D * D / 4;
  const sinT = grid.sinTheta;
  const sgd = Math.sqrt(G * D);
  const vdSpec = slip.vd ?? 'vertical';
  const vd0 = new Float64Array(N);
  for (let j = 0; j < N; j++) {
    const th = grid.theta[j];
    if (typeof vdSpec === 'number') vd0[j] = vdSpec;
    else if (typeof vdSpec === 'function') vd0[j] = vdSpec(th);
    else if (vdSpec === 'vertical') vd0[j] = 0.35 * sgd * Math.sin(th);
    else if (vdSpec === 'bendiksen') vd0[j] = (0.35 * Math.sin(th) + 0.54 * Math.cos(th)) * sgd;
    else if (vdSpec === 'buoyant') vd0[j] = th === 0 ? 0 : Math.sign(th) * (0.35 * Math.abs(Math.sin(th)) + 0.54 * Math.cos(th)) * sgd;
    else throw new RangeError(`slip.vd = ${vdSpec} is not a number, a function, 'vertical', 'bendiksen' or 'buoyant'`);
    checkArg('slip.vd', vd0[j]);
  }

  // --- boundary specification -------------------------------------------------------------------------
  const periodic = inlet.type === 'periodic' || outlet.type === 'periodic';
  if (periodic && inlet.type !== outlet.type) throw new Error('Inconsistent boundary conditions: a periodic boundary must be set at both ends.');
  let inL = null;
  let inG = null;
  if (inlet.type === 'massflow') {
    if (inlet.pressure !== undefined) throw new Error('Over-constrained inlet: phase mass rates and a pressure cannot both be imposed at a subsonic inlet.');
    inL = asFunction(inlet.massRateL, 'inlet.massRateL');
    inG = asFunction(inlet.massRateG, 'inlet.massRateG');
  } else if (inlet.type !== 'closed' && !periodic) throw new RangeError(`inlet.type = ${inlet.type} is not 'massflow', 'closed' or 'periodic'`);
  let outP = null;
  let outOpening = null;
  const backGas = checkArg('outlet.backflowGasFraction', outlet.backflowGasFraction ?? 1, { min: 0, max: 1 });
  if (outlet.type === 'pressure' || outlet.type === 'choke') {
    if (outlet.massRateL !== undefined || outlet.massRateG !== undefined) throw new Error('Over-constrained outlet: a subsonic outlet takes one condition (pressure or a valve law), not flow rates as well.');
    outP = asFunction(outlet.pressure, 'outlet.pressure');
    if (outlet.type === 'choke') outOpening = asFunction(outlet.opening, 'outlet.opening');
  } else if (outlet.type !== 'closed' && !periodic) throw new RangeError(`outlet.type = ${outlet.type} is not 'pressure', 'choke', 'closed' or 'periodic'`);
  const src = sources.map((q, i) => {
    checkArg(`sources[${i}].s`, q.s, { min: 0, max: grid.length });
    return { cell: Math.min(N - 1, Math.floor(q.s / ds)), rate: asFunction(q.massRateG, `sources[${i}].massRateG`) };
  });

  // --- state ------------------------------------------------------------------------------------------
  const m1 = new Float64Array(N);
  const m2 = new Float64Array(N);
  const mom = new Float64Array(N);
  const a1 = new Float64Array(N);
  const a2 = new Float64Array(N);
  const a3 = new Float64Array(N);
  const d1 = new Float64Array(N);
  const d2 = new Float64Array(N);
  const d3 = new Float64Array(N);
  const pr = new Float64Array(N);
  const ag = new Float64Array(N);
  const vl = new Float64Array(N);
  const vg = new Float64Array(N);
  const rhoGrav = new Float64Array(N);
  const gf = new Float64Array(N + 1);
  const f1 = new Float64Array(N + 1);
  const f2 = new Float64Array(N + 1);
  const f3 = new Float64Array(N + 1);
  const sAg = new Float64Array(N);
  const sP = new Float64Array(N);
  const sV = new Float64Array(N);
  const SL = { ag: 0, p: 0, vl: 0, vg: 0, m1: 0, m2: 0, c: 0 };
  const SR = { ag: 0, p: 0, vl: 0, vg: 0, m1: 0, m2: 0, c: 0 };
  const flux = { f1: 0, f2: 0, f3: 0 };
  let maxSpeed = 0;
  let clippedL = 0;
  let clippedG = 0;

  const fill = (st, agv, pv, vlv, vd) => {
    const a = agv < 0 ? 0 : (agv > 1 ? 1 : agv);
    const p = pv > P_FLOOR ? pv : P_FLOOR;
    const al = 1 - a;
    const rl = (p + B) / aL2;
    const rg = p / aG2;
    let C0 = C00;
    let v0 = vd;
    if (a > ab) { const w = ((1 - a) / (1 - ab)) ** 2; C0 = 1 + (C00 - 1) * w; v0 = vd * w; }
    const den = 1 - C0 * a;
    st.ag = a; st.p = p; st.vl = vlv;
    st.vg = den > 1e-12 ? (C0 * al * vlv + v0) / den : vlv;
    st.m1 = al * rl; st.m2 = a * rg;
    st.c = 1 / Math.sqrt((rl * den + rg * a) * (a / p + al / (rl * aL2)));
  };

  const primitives = (u1, u2, u3) => {
    maxSpeed = 0;
    for (let j = 0; j < N; j++) {
      const x1 = u1[j];
      const x2 = u2[j];
      const b = B - x1 * aL2 - x2 * aG2;
      const cc = x2 * aG2 * B;
      const root = Math.sqrt(b * b + 4 * cc);
      let p = b <= 0 ? 0.5 * (root - b) : 2 * cc / (b + root);
      if (!(p > P_FLOOR)) p = P_FLOOR;
      const rl = (p + B) / aL2;
      let a = 1 - x1 / rl;
      if (a < 0) a = 0; else if (a > 1) a = 1;
      const al = 1 - a;
      let C0 = C00;
      let v0 = vd0[j];
      if (a > ab) { const w = ((1 - a) / (1 - ab)) ** 2; C0 = 1 + (C00 - 1) * w; v0 *= w; }
      const den = 1 - C0 * a;
      let vL;
      let vG;
      if (al > 1e-9 && den > 1e-12) {
        const k = x1 + x2 * C0 * al / den;
        vL = k > 1e-300 ? (u3[j] - x2 * v0 / den) / k : 0;
        vG = (C0 * al * vL + v0) / den;
      } else { vL = x1 + x2 > 1e-300 ? u3[j] / (x1 + x2) : 0; vG = vL; }
      pr[j] = p; ag[j] = a; vl[j] = vL; vg[j] = vG;
      rhoGrav[j] = (x1 + x2) * sinT[j];
      const c = 1 / Math.sqrt((rl * den + x2) * (a / p + al / (rl * aL2)));
      const sp = Math.max(Math.abs(vL), Math.abs(vG)) + c;
      if (sp > maxSpeed) maxSpeed = sp;
    }
  };

  const ausmv = (L, R) => {
    const c = L.c > R.c ? L.c : R.c;
    let F1 = 0;
    let F2 = 0;
    let F3 = 0;
    // liquid
    let sum = L.m1 + R.m1;
    if (sum > 1e-300) {
      const chiL = 2 * R.m1 / sum;
      const chiR = 2 * L.m1 / sum;
      const vp = L.vl;
      const vm = R.vl;
      const up = 0.5 * (vp + Math.abs(vp));
      const um = 0.5 * (vm - Math.abs(vm));
      const Vp = Math.abs(vp) <= c ? chiL * (vp + c) * (vp + c) / (4 * c) + (1 - chiL) * up : up;
      const Vm = Math.abs(vm) <= c ? -chiR * (vm - c) * (vm - c) / (4 * c) + (1 - chiR) * um : um;
      F1 = L.m1 * Vp + R.m1 * Vm;
      F3 += L.m1 * vp * Vp + R.m1 * vm * Vm;
    }
    // gas
    sum = L.m2 + R.m2;
    if (sum > 1e-300) {
      const chiL = 2 * R.m2 / sum;
      const chiR = 2 * L.m2 / sum;
      const vp = L.vg;
      const vm = R.vg;
      const up = 0.5 * (vp + Math.abs(vp));
      const um = 0.5 * (vm - Math.abs(vm));
      const Vp = Math.abs(vp) <= c ? chiL * (vp + c) * (vp + c) / (4 * c) + (1 - chiL) * up : up;
      const Vm = Math.abs(vm) <= c ? -chiR * (vm - c) * (vm - c) / (4 * c) + (1 - chiR) * um : um;
      F2 = L.m2 * Vp + R.m2 * Vm;
      F3 += L.m2 * vp * Vp + R.m2 * vm * Vm;
    }
    const rL = L.m1 + L.m2;
    const rR = R.m1 + R.m2;
    const zL = rL * L.c;
    const zR = rR * R.c;
    if (Kp > 0) {
      // Pressure-diffusion velocity (pL − pR)/(zL + zR) of the acoustic Riemann problem, scaled by 2 Kp; it
      // vanishes for hydrostatically balanced face pressures.
      const u = 2 * Kp * (L.p - R.p) / (zL + zR);
      if (u > 0) { F1 += L.m1 * u; F2 += L.m2 * u; F3 += (L.m1 * L.vl + L.m2 * L.vg) * u; }
      else { F1 += R.m1 * u; F2 += R.m2 * u; F3 += (R.m1 * R.vl + R.m2 * R.vg) * u; }
    }
    const wL = (L.m1 * L.vl + L.m2 * L.vg) / rL;
    const wR = (R.m1 * R.vl + R.m2 * R.vg) / rR;
    // Interface pressure: impedance-weighted acoustic value at low Mach number, blended into the upwinding
    // P± splitting as the mixture Mach number approaches one.
    const pAc = (zR * L.p + zL * R.p + zL * zR * (wL - wR)) / (zL + zR);
    const M2 = Math.max(wL * wL, wR * wR) / (c * c);
    let pFace = pAc;
    if (M2 > 0.25) {
      const Pp = Math.abs(wL) <= c ? (wL + c) * (wL + c) * (2 - wL / c) / (4 * c * c) : (wL > 0 ? 1 : 0);
      const Pm = Math.abs(wR) <= c ? (wR - c) * (wR - c) * (2 + wR / c) / (4 * c * c) : (wR < 0 ? 1 : 0);
      const f = M2 >= 1 ? 0 : (1 - M2) / 0.75;
      pFace = f * pAc + (1 - f) * (Pp * L.p + Pm * R.p);
    }
    flux.f1 = F1; flux.f2 = F2; flux.f3 = F3 + pFace;
  };

  // Phase state carrying stated mass fluxes at a given pressure (inlet boundary and steady initial profile).
  const fromFluxes = (Gl, Gg, p, vd, out) => {
    const rl = (p + B) / aL2;
    const rg = p / aG2;
    let a = 0;
    let vG = 0;
    if (Gg !== 0) {
      let C0 = C00;
      let v0 = vd;
      for (let it = 0; it < 30; it++) {
        vG = C0 * (Gl / rl + Gg / rg) + v0;
        let an = vG > 1e-9 ? Gg / (rg * vG) : 1;
        if (an > 1) an = 1; else if (an < 0) an = 0;
        const rel = it < 4 ? 1 : 0.5;
        const next = a + rel * (an - a);
        if (Math.abs(next - a) < 1e-14) { a = next; break; }
        a = next;
        if (a > ab) { const w = ((1 - a) / (1 - ab)) ** 2; C0 = 1 + (C00 - 1) * w; v0 = vd * w; } else { C0 = C00; v0 = vd; }
      }
      a = Math.min(a, 1 - 1e-9);
      vG = a > 0 ? Gg / (rg * a) : 0;
    }
    out.ag = a; out.p = p; out.vl = Gl / (rl * (1 - a)); out.vg = vG; out.m1 = (1 - a) * rl; out.m2 = a * rg;
  };

  const relRough = roughness / D;
  const frictionGrad = (rho, vm, a) => {
    if (!friction || vm === 0) return 0;
    const mu = (1 - a) * muL + a * muG;
    return darcyFrictionFactorFast(rho * Math.abs(vm) * D / mu, relRough) * rho * vm * Math.abs(vm) / (2 * D);
  };

  // Boundary and source mass rates per unit area of the last residual evaluation (for conservation bookkeeping).
  const bnd = { inL: 0, inG: 0, outL: 0, outG: 0, srcG: 0 };

  const residual = (u1, u2, u3, t, r1, r2, r3) => {
    primitives(u1, u2, u3);
    // Face-averaged ρ sinθ: the hydrostatic pressure jump between neighbouring cell centres is g Δs gf.
    gf[0] = periodic ? 0.5 * (rhoGrav[N - 1] + rhoGrav[0]) : rhoGrav[0];
    gf[N] = periodic ? gf[0] : rhoGrav[N - 1];
    for (let j = 1; j < N; j++) gf[j] = 0.5 * (rhoGrav[j - 1] + rhoGrav[j]);
    const gds = G * ds;
    for (let j = 0; j < N; j++) {
      const jm = j > 0 ? j - 1 : (periodic ? N - 1 : -1);
      const jp = j < N - 1 ? j + 1 : (periodic ? 0 : -1);
      if (jm < 0 || jp < 0) { sAg[j] = 0; sP[j] = 0; sV[j] = 0; continue; }
      sAg[j] = lim(ag[j] - ag[jm], ag[jp] - ag[j]);
      sV[j] = lim(vl[j] - vl[jm], vl[jp] - vl[j]);
      sP[j] = lim(pr[j] - pr[jm] + gds * gf[j], pr[jp] - pr[j] + gds * gf[j + 1]);
    }
    const first = periodic ? 0 : 1;
    const last = periodic ? N : N - 1;
    for (let f = first; f <= last; f++) {
      const jl = f === 0 ? N - 1 : f - 1;
      const jr = f === N ? 0 : f;
      // Each cell extrapolates its own hydrostatic profile to the face (hydrostatic reconstruction).
      fill(SL, ag[jl] + 0.5 * sAg[jl], pr[jl] + 0.5 * sP[jl] - 0.5 * gds * rhoGrav[jl], vl[jl] + 0.5 * sV[jl], vd0[jl]);
      fill(SR, ag[jr] - 0.5 * sAg[jr], pr[jr] - 0.5 * sP[jr] + 0.5 * gds * rhoGrav[jr], vl[jr] - 0.5 * sV[jr], vd0[jr]);
      ausmv(SL, SR);
      f1[f] = flux.f1; f2[f] = flux.f2; f3[f] = flux.f3;
    }
    if (periodic) { f1[N] = f1[0]; f2[N] = f2[0]; f3[N] = f3[0]; bnd.inL = bnd.inG = bnd.outL = bnd.outG = 0; } else {
      // inlet face
      const pIn = pr[0] + 0.5 * gds * gf[0];
      if (inL) {
        const Gl = inL(t) / A;
        const Gg = inG(t) / A;
        if (Gl < 0 || Gg < 0) throw new RangeError('inlet mass rates must be non-negative');
        if (Gl === 0 && Gg === 0) { f1[0] = 0; f2[0] = 0; f3[0] = pIn; } else {
          fromFluxes(Gl, Gg, pIn > P_FLOOR ? pIn : P_FLOOR, vd0[0], SL);
          f1[0] = Gl; f2[0] = Gg; f3[0] = Gl * SL.vl + Gg * SL.vg + pIn;
        }
      } else { f1[0] = 0; f2[0] = 0; f3[0] = pIn; }
      bnd.inL = f1[0]; bnd.inG = f2[0];
      // outlet face
      const e = N - 1;
      const pFace = pr[e] - 0.5 * gds * gf[N];
      if (!outP) { f1[N] = 0; f2[N] = 0; f3[N] = pFace; } else if (outOpening) {
        const pd = outP(t);
        const z = Math.min(1, Math.max(0, outOpening(t)));
        const rho = u1[e] + u2[e];
        if (pFace >= pd) {
          const v = z * Math.sqrt(2 * (pFace - pd) / rho);
          f1[N] = u1[e] * v; f2[N] = u2[e] * v; f3[N] = u1[e] * v * vl[e] + u2[e] * v * vg[e] + pFace;
        } else {
          const rg = pd / aG2;
          const v = -z * Math.sqrt(2 * (pd - pFace) / rg);
          f1[N] = 0; f2[N] = rg * v; f3[N] = rg * v * v + pFace;
        }
      } else {
        fill(SL, ag[e], pFace, vl[e], vd0[e]);
        // Liquid cannot return from the receiving vessel: when the liquid at the outlet falls back, the ghost
        // state is the vessel atmosphere (gas fraction backflowGasFraction).
        fill(SR, vl[e] >= 0 ? ag[e] : Math.max(ag[e], backGas), outP(t), vl[e], vd0[e]);
        ausmv(SL, SR);
        f1[N] = flux.f1; f2[N] = flux.f2; f3[N] = flux.f3;
      }
      bnd.outL = f1[N]; bnd.outG = f2[N];
    }
    for (let j = 0; j < N; j++) {
      const rho = u1[j] + u2[j];
      const vm = (1 - ag[j]) * vl[j] + ag[j] * vg[j];
      r1[j] = -(f1[j + 1] - f1[j]) / ds;
      r2[j] = -(f2[j + 1] - f2[j]) / ds;
      r3[j] = -(f3[j + 1] - f3[j]) / ds - G * rhoGrav[j] - frictionGrad(rho, vm, ag[j]);
    }
    bnd.srcG = 0;
    for (const q of src) { const rate = q.rate(t); r2[q.cell] += rate / (A * ds); bnd.srcG += rate / A; }
  };

  // --- initial condition ------------------------------------------------------------------------------
  const setCell = (j, p, a, vL) => {
    fill(SL, a, p, vL, vd0[j]);
    m1[j] = SL.m1; m2[j] = SL.m2; mom[j] = SL.m1 * SL.vl + SL.m2 * SL.vg;
  };
  const cellValue = (v, j, name) => {
    if (typeof v === 'function') return v(grid.s[j]);
    if (typeof v === 'number') return v;
    if (v && v.length === N) return v[j];
    throw new RangeError(`initial.${name} must be a number, a function of chainage or an array with one value per cell`);
  };
  if (initial.type === 'uniform' || initial.type === 'profile') {
    for (let j = 0; j < N; j++) {
      const p = checkArg('initial.P', cellValue(initial.P, j, 'P'), { min: 0, exclusiveMin: true });
      const a = checkArg('initial.alphaG', cellValue(initial.alphaG, j, 'alphaG'), { min: 0, max: 1 });
      setCell(j, p, a, initial.vL === undefined ? 0 : cellValue(initial.vL, j, 'vL'));
    }
  } else if (initial.type === 'steady') {
    const Gl = checkArg('initial.massRateL', initial.massRateL, { min: 0 }) / A;
    const Gg = checkArg('initial.massRateG', initial.massRateG, { min: 0 }) / A;
    let p = checkArg('initial.pOut', initial.pOut, { min: 0, exclusiveMin: true });
    const grad = (pp, j) => {
      fromFluxes(Gl, Gg, pp, vd0[j], SL);
      const rho = SL.m1 + SL.m2;
      return frictionGrad(rho, (1 - SL.ag) * SL.vl + SL.ag * SL.vg, SL.ag) + rho * G * sinT[j];
    };
    for (let j = N - 1; j >= 0; j--) {
      const h = j === N - 1 ? 0.5 * ds : ds;
      const jj = j === N - 1 ? j : j + 1;
      const k1 = grad(p, jj);
      const k2 = grad(p + h * k1, j);
      p += 0.5 * h * (k1 + k2);
      if (!(p > 0) || !Number.isFinite(p)) throw new Error('Steady initial profile failed: pressure not finite.');
      fromFluxes(Gl, Gg, p, vd0[j], SL);
      m1[j] = SL.m1; m2[j] = SL.m2; mom[j] = SL.m1 * SL.vl + SL.m2 * SL.vg;
    }
  } else if (initial.type === 'settled') {
    const frac = checkArg('initial.liquidFraction', initial.liquidFraction, { min: 0, max: 1 });
    const pRef = checkArg('initial.pRef', initial.pRef, { min: 0, exclusiveMin: true });
    const ref = initial.refIndex === undefined ? N - 1 : checkArg('initial.refIndex', initial.refIndex, { min: 0, max: N - 1 });
    const order = Array.from({ length: N }, (_, j) => j).sort((i, k) => grid.z[i] - grid.z[k] || i - k);
    const nLiq = Math.round(frac * N);
    const a = new Float64Array(N).fill(1);
    for (let k = 0; k < nLiq; k++) a[order[k]] = 0;
    const p = new Float64Array(N);
    const rhoS = (j, pj) => ((1 - a[j]) * (pj + B) / aL2 + a[j] * pj / aG2) * sinT[j];
    p[ref] = pRef;
    // Discrete hydrostatic balance p[j+1] − p[j] = −g Δs ½(ρ sinθ|j + ρ sinθ|j+1), solved by fixed point per face.
    for (let j = ref; j < N - 1; j++) {
      let pn = p[j];
      for (let it = 0; it < 60; it++) pn = p[j] - 0.5 * G * ds * (rhoS(j, p[j]) + rhoS(j + 1, pn));
      p[j + 1] = pn;
    }
    for (let j = ref; j > 0; j--) {
      let pn = p[j];
      for (let it = 0; it < 60; it++) pn = p[j] + 0.5 * G * ds * (rhoS(j, p[j]) + rhoS(j - 1, pn));
      p[j - 1] = pn;
    }
    for (let j = 0; j < N; j++) {
      if (!(p[j] > 0)) throw new Error('Settled initial state failed: hydrostatic pressure not positive.');
      setCell(j, p[j], a[j], 0);
    }
  } else throw new RangeError(`initial.type = ${initial.type} is not 'uniform', 'profile', 'steady' or 'settled'`);

  // --- time integration -------------------------------------------------------------------------------
  const total = (u) => { let s = 0; for (let j = 0; j < N; j++) s += u[j]; return s * A * ds; };
  const massL0 = total(m1);
  const massG0 = total(m2);
  let netL = 0;
  let netG = 0;
  const hT = new Float64Array(nHistory);
  const hL = new Float64Array(nHistory);
  const hG = new Float64Array(nHistory);
  const hPin = new Float64Array(nHistory);
  const hPout = new Float64Array(nHistory);
  const hPmax = new Float64Array(nHistory);
  const hInv = new Float64Array(nHistory);
  const hDt = tEnd / nHistory;
  let hIdx = 0;
  let accL = 0;
  let accG = 0;
  let accT = 0;
  const prof = { t: [], P: [], holdup: [], vL: [], vG: [] };
  const store = (t) => {
    primitives(m1, m2, mom);
    prof.t.push(t); prof.P.push(Float64Array.from(pr));
    prof.holdup.push(Float64Array.from(ag, (a) => 1 - a)); prof.vL.push(Float64Array.from(vl)); prof.vG.push(Float64Array.from(vg));
  };
  store(0);
  const pDt = tEnd / nProfiles;
  let pIdx = 1;
  let t = 0;
  let steps = 0;
  let dtMin = Infinity;
  let dtMax = 0;
  const clip = (u, which) => {
    for (let j = 0; j < N; j++) {
      if (u[j] < 0) { if (which === 1) clippedL -= u[j] * A * ds; else clippedG -= u[j] * A * ds; u[j] = 0; }
    }
  };
  while (t < tEnd * (1 - 1e-13)) {
    if (steps >= maxSteps) throw new Error(`Transient run stopped after ${maxSteps} steps at t = ${t.toFixed(3)} s; raise maxSteps or lower fluid.aL.`);
    residual(m1, m2, mom, t, d1, d2, d3);
    if (!Number.isFinite(maxSpeed) || !(maxSpeed > 0)) throw new Error(`Transient run failed at t = ${t.toFixed(4)} s (non-finite state).`);
    let dt = cfl * ds / maxSpeed;
    const tNext = Math.min(tEnd, hIdx < nHistory ? (hIdx + 1) * hDt : tEnd);
    if (t + dt > tNext) dt = tNext - t;
    let bL = bnd.inL - bnd.outL;
    let bG = bnd.inG - bnd.outG + bnd.srcG;
    let oL = bnd.outL;
    let oG = bnd.outG;
    for (let j = 0; j < N; j++) { a1[j] = m1[j] + dt * d1[j]; a2[j] = m2[j] + dt * d2[j]; a3[j] = mom[j] + dt * d3[j]; }
    clip(a1, 1); clip(a2, 2);
    residual(a1, a2, a3, t + dt, d1, d2, d3);
    bL = 0.5 * (bL + bnd.inL - bnd.outL);
    bG = 0.5 * (bG + bnd.inG - bnd.outG + bnd.srcG);
    oL = 0.5 * (oL + bnd.outL);
    oG = 0.5 * (oG + bnd.outG);
    for (let j = 0; j < N; j++) {
      m1[j] = 0.5 * m1[j] + 0.5 * (a1[j] + dt * d1[j]);
      m2[j] = 0.5 * m2[j] + 0.5 * (a2[j] + dt * d2[j]);
      mom[j] = 0.5 * mom[j] + 0.5 * (a3[j] + dt * d3[j]);
    }
    clip(m1, 1); clip(m2, 2);
    netL += bL * A * dt;
    netG += bG * A * dt;
    accL += oL * A * dt;
    accG += oG * A * dt;
    accT += dt;
    t += dt;
    steps++;
    if (dt < dtMin) dtMin = dt;
    if (dt > dtMax) dtMax = dt;
    if (hIdx < nHistory && t >= (hIdx + 1) * hDt * (1 - 1e-12)) {
      primitives(m1, m2, mom);
      let inv = 0;
      let pm = 0;
      for (let j = 0; j < N; j++) { inv += 1 - ag[j]; if (pr[j] > pm) pm = pr[j]; }
      hT[hIdx] = t - 0.5 * accT; hL[hIdx] = accL / accT; hG[hIdx] = accG / accT;
      hPin[hIdx] = pr[0]; hPout[hIdx] = pr[N - 1]; hPmax[hIdx] = pm; hInv[hIdx] = inv * A * ds;
      hIdx++; accL = 0; accG = 0; accT = 0;
    }
    if (pIdx <= nProfiles && t >= pIdx * pDt * (1 - 1e-12)) { store(t); pIdx++; }
  }
  primitives(m1, m2, mom);
  const massL = total(m1);
  const massG = total(m2);
  return {
    grid, profiles: prof,
    history: {
      t: hT.subarray(0, hIdx), outletLiquidRate: hL.subarray(0, hIdx), outletGasRate: hG.subarray(0, hIdx),
      inletPressure: hPin.subarray(0, hIdx), outletPressure: hPout.subarray(0, hIdx), maxPressure: hPmax.subarray(0, hIdx),
      liquidInventory: hInv.subarray(0, hIdx),
    },
    final: { P: Float64Array.from(pr), alphaG: Float64Array.from(ag), holdup: Float64Array.from(ag, (a) => 1 - a), vL: Float64Array.from(vl), vG: Float64Array.from(vg) },
    conservation: {
      massL0, massG0, massL, massG, netInflowL: netL, netInflowG: netG,
      residualL: (massL - massL0 - netL) / Math.max(massL0, Math.abs(netL), 1e-300),
      residualG: (massG - massG0 - netG) / Math.max(massG0, Math.abs(netG), 1e-300), clippedL, clippedG,
    },
    steps, time: t, dtMin, dtMax,
  };
}

/**
 * Oscillation statistics of a signal after an initial transient: mean, range and period from upward
 * crossings of the mean.
 * @param {object} o
 * @param {ArrayLike<number>} o.t times [s], ascending
 * @param {ArrayLike<number>} o.y signal [any unit]
 * @param {number} [o.discardFraction=0.5] leading fraction of the record ignored [-]
 * @returns {{mean: number, min: number, max: number, relativeAmplitude: number, cycles: number, period: number}}
 *   mean, min, max [unit of y]; relativeAmplitude = (max − min)/|mean| [-]; cycles = number of complete
 *   cycles [-]; period = mean spacing of upward mean-crossings [s] (NaN with fewer than two crossings)
 */
export function analyseOscillation({ t, y, discardFraction = 0.5 }) {
  checkArg('discardFraction', discardFraction, { min: 0, max: 0.99 });
  const n = y.length;
  const i0 = Math.floor(discardFraction * n);
  if (n - i0 < 3) throw new RangeError('y is too short for the stated discardFraction');
  let mean = 0;
  let min = Infinity;
  let max = -Infinity;
  for (let i = i0; i < n; i++) { mean += y[i]; if (y[i] < min) min = y[i]; if (y[i] > max) max = y[i]; }
  mean /= n - i0;
  const cross = [];
  // Hysteresis of 10 % of the half-range keeps noise near the mean from counting as cycles.
  const band = 0.1 * 0.5 * (max - min);
  let armed = false;
  for (let i = i0; i < n; i++) {
    if (y[i] < mean - band) armed = true;
    else if (armed && y[i] > mean + band) { cross.push(t[i]); armed = false; }
  }
  const cycles = Math.max(0, cross.length - 1);
  return {
    mean, min, max, relativeAmplitude: mean !== 0 ? (max - min) / Math.abs(mean) : Infinity, cycles,
    period: cycles > 0 ? (cross[cross.length - 1] - cross[0]) / cycles : NaN,
  };
}

/**
 * Liquid-surge (slug-arrival) statistics of an outlet liquid-rate history: a surge is a run of samples above
 * threshold = thresholdFactor × time-mean rate.
 * @param {object} o
 * @param {ArrayLike<number>} o.t sample times [s] (uniform spacing assumed for volumes)
 * @param {ArrayLike<number>} o.rate liquid rate [kg/s or m³/s]
 * @param {number} [o.thresholdFactor=1.5] surge threshold relative to the mean rate [-]
 * @param {number} [o.discardFraction=0] leading fraction of the record ignored [-]
 * @returns {{surges: {tStart: number, tEnd: number, peak: number, amount: number}[], count: number,
 *   meanRate: number, peakRate: number, meanPeriod: number, periodCov: number, maxAmount: number,
 *   surgeFraction: number}} surges with start/end [s], peak rate and amount = ∫ rate dt [kg or m³]; meanPeriod
 *   between surge starts [s] and its coefficient of variation [-] (NaN with fewer than two / three surges);
 *   surgeFraction = share of the liquid delivered in surges [-]
 */
export function detectSlugArrivals({ t, rate, thresholdFactor = 1.5, discardFraction = 0 }) {
  checkArg('thresholdFactor', thresholdFactor, { min: 1, exclusiveMin: true });
  checkArg('discardFraction', discardFraction, { min: 0, max: 0.99 });
  const n = rate.length;
  const i0 = Math.floor(discardFraction * n);
  if (n - i0 < 3) throw new RangeError('rate is too short');
  const dt = (t[n - 1] - t[i0]) / (n - 1 - i0);
  let mean = 0;
  let peakRate = -Infinity;
  for (let i = i0; i < n; i++) { mean += rate[i]; if (rate[i] > peakRate) peakRate = rate[i]; }
  mean /= n - i0;
  const thr = thresholdFactor * mean;
  const surges = [];
  let cur = null;
  let inSurges = 0;
  for (let i = i0; i < n; i++) {
    if (rate[i] > thr && mean > 0) {
      if (!cur) cur = { tStart: t[i] - dt / 2, tEnd: 0, peak: 0, amount: 0 };
      cur.amount += rate[i] * dt;
      if (rate[i] > cur.peak) cur.peak = rate[i];
      cur.tEnd = t[i] + dt / 2;
    } else if (cur) { surges.push(cur); inSurges += cur.amount; cur = null; }
  }
  if (cur) { surges.push(cur); inSurges += cur.amount; }
  const gaps = [];
  for (let i = 1; i < surges.length; i++) gaps.push(surges[i].tStart - surges[i - 1].tStart);
  const mp = gaps.length ? gaps.reduce((a, b) => a + b, 0) / gaps.length : NaN;
  const cov = gaps.length > 1 ? Math.sqrt(gaps.reduce((a, b) => a + (b - mp) ** 2, 0) / (gaps.length - 1)) / mp : NaN;
  let maxAmount = 0;
  for (const s of surges) if (s.amount > maxAmount) maxAmount = s.amount;
  return {
    surges, count: surges.length, meanRate: mean, peakRate, meanPeriod: mp, periodCov: cov, maxAmount,
    surgeFraction: mean > 0 ? inSurges / (mean * dt * (n - i0)) : 0,
  };
}

function ramp(t0, dur, a, b) {
  return (t) => (t <= t0 ? a : (t >= t0 + dur ? b : a + (b - a) * (t - t0) / dur));
}

/**
 * Production ramp-up: the line starts from the steady state at the initial rates and the inlet rates are ramped
 * linearly to the final rates; returns the transient result with the surge statistics of the outlet liquid.
 * @param {object} o simulateTransient arguments without inlet/initial, plus
 * @param {number} o.pOut outlet pressure [Pa] (used when no outlet is given)
 * @param {{massRateL: number, massRateG: number}} o.from initial rates [kg/s]
 * @param {{massRateL: number, massRateG: number}} o.to final rates [kg/s]
 * @param {number} [o.rampStart=0] start of the ramp [s]
 * @param {number} o.rampDuration ramp duration [s]
 * @returns {object} simulateTransient result plus `surge` (detectSlugArrivals of the outlet liquid rate) and
 *   liquidSwept = initial − final liquid inventory [m³]
 */
export function rampUp({ pOut, from, to, rampStart = 0, rampDuration, outlet, ...rest }) {
  checkArg('rampDuration', rampDuration, { min: 0, exclusiveMin: true });
  const r = simulateTransient({
    ...rest,
    inlet: { type: 'massflow', massRateL: ramp(rampStart, rampDuration, from.massRateL, to.massRateL), massRateG: ramp(rampStart, rampDuration, from.massRateG, to.massRateG) },
    outlet: outlet ?? { type: 'pressure', pressure: pOut },
    initial: { type: 'steady', massRateL: from.massRateL, massRateG: from.massRateG, pOut },
  });
  const inv = r.history.liquidInventory;
  return {
    ...r, surge: detectSlugArrivals({ t: r.history.t, rate: r.history.outletLiquidRate, thresholdFactor: 1.2 }),
    liquidSwept: (r.profiles.holdup[0].reduce((a, b) => a + b, 0) * r.grid.ds * Math.PI * rest.D * rest.D / 4) - inv[inv.length - 1],
  };
}

/**
 * Shutdown and restart: steady production, inlet rates ramped to zero at tShutdown (the outlet stays on
 * pressure control so the line depressurises to the outlet pressure and the liquid settles), then ramped
 * back at tRestart.
 * @param {object} o simulateTransient arguments without inlet/initial, plus
 * @param {number} o.pOut outlet pressure [Pa]
 * @param {number} o.massRateL liquid rate [kg/s]
 * @param {number} o.massRateG gas rate [kg/s]
 * @param {number} o.tShutdown time the inlet starts closing [s]
 * @param {number} o.tRestart time the inlet starts re-opening [s], > tShutdown + valveTime
 * @param {number} [o.valveTime=10] closing / opening ramp time [s]
 * @returns {object} simulateTransient result plus `settledInventory` (liquid in the line just before restart
 *   [m³]) and `restartSurge` (detectSlugArrivals over the part of the record after tRestart)
 */
export function shutdownRestart({ pOut, massRateL, massRateG, tShutdown, tRestart, valveTime = 10, outlet, ...rest }) {
  checkArg('tShutdown', tShutdown, { min: 0 });
  checkArg('valveTime', valveTime, { min: 0, exclusiveMin: true });
  checkArg('tRestart', tRestart, { min: tShutdown + valveTime });
  const f = (q) => (t) => (t < tRestart ? ramp(tShutdown, valveTime, q, 0)(t) : ramp(tRestart, valveTime, 0, q)(t));
  const r = simulateTransient({
    ...rest, inlet: { type: 'massflow', massRateL: f(massRateL), massRateG: f(massRateG) },
    outlet: outlet ?? { type: 'pressure', pressure: pOut }, initial: { type: 'steady', massRateL, massRateG, pOut },
  });
  const h = r.history;
  let k = 0;
  while (k < h.t.length - 1 && h.t[k] < tRestart) k++;
  const frac = Math.min(0.98, k / h.t.length);
  let restartSurge = null;
  if (h.t.length - k >= 3) restartSurge = detectSlugArrivals({ t: h.t, rate: h.outletLiquidRate, thresholdFactor: 1.2, discardFraction: frac });
  return { ...r, settledInventory: h.liquidInventory[Math.max(0, k - 1)], restartSurge };
}

/**
 * Riser-induced (terrain) slugging run in a flowline + riser at constant inlet mass rates, with optional
 * mitigation by topside choking (outlet choke) or riser-base gas lift. Returns the transient result, the
 * oscillation statistics of the riser-base pressure and outlet liquid rate over the last part of the run, and
 * the Bøe and Taitel criteria evaluated with the time-mean flowline state of the run.
 * @param {object} o
 * @param {number} o.flowlineLength flowline length [m]
 * @param {number} [o.flowlineAngleDeg=-2] flowline inclination [deg]
 * @param {number} o.riserHeight riser height [m]
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} [o.roughness=0] wall roughness [m]
 * @param {number} o.massRateL liquid mass rate [kg/s]
 * @param {number} o.massRateG gas mass rate [kg/s]
 * @param {number} o.pSep separator pressure [Pa]
 * @param {object} o.fluid fluid description (see simulateTransient)
 * @param {object} [o.slip] slip closure (see simulateTransient)
 * @param {number} [o.chokeOpening] outlet choke opening Cd·Ac/A [-]; omitted = no choke (pressure outlet)
 * @param {number} [o.gasLiftRate=0] gas-lift mass rate injected at the riser base [kg/s]
 * @param {number} [o.upstreamGasVolume=0] extra gas buffer volume upstream of the flowline [m³], modelled as an
 *   added horizontal pipe length of the same diameter (laboratory rigs use such a vessel to represent a long line)
 * @param {number} o.nCells number of cells [-]
 * @param {number} o.tEnd simulated time [s]
 * @param {number} [o.discardFraction=0.5] part of the record treated as start-up [-]
 * @param {number} [o.unstableAmplitude=0.2] relative peak-to-peak amplitude of the riser-base pressure head
 *   above which the run is classed as slugging [-] (relative to the full-riser liquid head)
 * @returns {object} simulateTransient result plus {riserBaseIndex, riserBasePressure: Float64Array [Pa] (stored
 *   profiles), pressure: analyseOscillation of the inlet pressure, liquid: analyseOscillation of the outlet
 *   liquid rate, headAmplitude [-] = (pmax − pmin)/(ρL g h), slugging: boolean, period [s] (NaN when not slugging), boe, taitel,
 *   flowlineGasFraction [-], flowlinePressure [Pa]}
 */
export function severeSluggingCase({ flowlineLength, flowlineAngleDeg = -2, riserHeight, D, roughness = 0, massRateL, massRateG,
  pSep, fluid, slip, chokeOpening, gasLiftRate = 0, upstreamGasVolume = 0, nCells, tEnd, discardFraction = 0.5,
  unstableAmplitude = 0.2, ...numerics }) {
  checkArg('massRateL', massRateL, { min: 0, exclusiveMin: true });
  checkArg('massRateG', massRateG, { min: 0, exclusiveMin: true });
  checkArg('pSep', pSep, { min: 0, exclusiveMin: true });
  checkArg('gasLiftRate', gasLiftRate, { min: 0 });
  checkArg('upstreamGasVolume', upstreamGasVolume, { min: 0 });
  const A = Math.PI * D * D / 4;
  const buffer = upstreamGasVolume / A;
  const base = flowlineRiserProfile({ flowlineLength, flowlineAngleDeg, riserHeight });
  const profile = buffer > 0 ? [{ x: -buffer, z: 0 }, ...base] : base;
  const sBase = buffer + flowlineLength;
  const r = simulateTransient({
    ...numerics, profile, D, roughness, nCells, tEnd, fluid, slip,
    inlet: { type: 'massflow', massRateL, massRateG },
    outlet: chokeOpening === undefined ? { type: 'pressure', pressure: pSep } : { type: 'choke', pressure: pSep, opening: chokeOpening },
    sources: gasLiftRate > 0 ? [{ s: sBase + 0.02 * riserHeight, massRateG: gasLiftRate }] : [],
    initial: { type: 'steady', massRateL, massRateG, pOut: pSep },
  });
  const g = r.grid;
  const kBase = Math.min(g.n - 1, Math.max(0, Math.floor(sBase / g.ds) - 1));
  const pressure = analyseOscillation({ t: r.history.t, y: r.history.maxPressure, discardFraction });
  const liquid = analyseOscillation({ t: r.history.t, y: r.history.outletLiquidRate, discardFraction });
  const headAmplitude = (pressure.max - pressure.min) / (fluid.rhoL * G * riserHeight);
  // Time-mean flowline state over the analysed part of the run (stored profiles).
  const k0 = Math.floor(discardFraction * r.profiles.t.length);
  let aSum = 0;
  let pSum = 0;
  let cnt = 0;
  const j0 = Math.min(kBase, Math.floor(buffer / g.ds));
  for (let k = k0; k < r.profiles.t.length; k++) {
    for (let j = j0; j <= kBase; j++) { aSum += 1 - r.profiles.holdup[k][j]; pSum += r.profiles.P[k][j]; cnt++; }
  }
  const alphaPipe = Math.min(1, Math.max(1e-6, aSum / cnt));
  const pPipe = pSum / cnt;
  const vsl = massRateL / (fluid.rhoL * A);
  const vsg = massRateG * fluid.gasRT / (pPipe * A);
  const Leff = flowlineLength + buffer / alphaPipe;
  return {
    ...r, riserBaseIndex: kBase, riserBasePressure: Float64Array.from(r.profiles.P, (p) => p[kBase]), pressure, liquid,
    headAmplitude, slugging: headAmplitude > unstableAmplitude, period: headAmplitude > unstableAmplitude ? pressure.period : NaN,
    boe: boeCriterion({ vsl, vsg, P: pPipe, rhoL: fluid.rhoL, alphaPipe, L: Leff }),
    taitel: taitelStability({ Psep: pSep, rhoL: fluid.rhoL, alphaPipe, L: Leff, riserHeight }),
    flowlineGasFraction: alphaPipe, flowlinePressure: pPipe,
  };
}

/**
 * Lagrangian slug tracking over a uniform stratified film: every slug is an object with a tail and a front
 * position. Tail (bubble-nose) velocity vt = C0 vm + vd, optionally increased for short slugs by the wake
 * factor 1 + a·exp(−b Ls/D); front velocity from the liquid balance across the front (pickup of the film
 * ahead), vF = (HLS vm − HLF vLF)/(HLS − HLF). A slug grows when vF > vt (pickup exceeds shedding), decays
 * otherwise and is removed at zero length; a slug whose front reaches the tail of the slug ahead merges with
 * it. Positions advance explicitly (exact for constant velocities).
 * Reference: Taitel & Barnea (1998), Int. J. Multiphase Flow 24(8), 1235–1257 (tracking concept); Nydal &
 * Banerjee (1996), Chem. Eng. Commun. 141–142, 13–39; wake law of the form proposed by Moissis & Griffith
 * (1962), J. Heat Transfer 84, 29–39 (their vertical-flow coefficients are a = 8, b = 1.06).
 * Simplifications: film holdup and velocity are prescribed and uniform (no film dynamics, no gas
 * compressibility), so liquid is not conserved between slugs.
 * @param {object} o
 * @param {number} o.length pipe length [m]
 * @param {number} o.D pipe inner diameter [m]
 * @param {number} o.vm mixture velocity [m/s]
 * @param {number} o.vt tail velocity of a long slug C0 vm + vd [m/s]
 * @param {number} o.slugHoldup HLS [-]
 * @param {number} o.filmHoldup HLF [-], < HLS
 * @param {number} o.vLF film liquid velocity [m/s]
 * @param {{tail: number, front: number}[]} [o.slugs=[]] initial slugs, positions [m]
 * @param {number} [o.inletPeriod] period at which new slugs of length inletSlugLength are created at s = 0 [s]
 * @param {number} [o.inletSlugLength] initial length of inlet slugs [m]
 * @param {number} [o.wakeA=0] wake coefficient a [-] (0 = no wake effect)
 * @param {number} [o.wakeB=1.06] wake decay coefficient b [-]
 * @param {number} o.tEnd simulated time [s]
 * @param {number} o.dt time step [s]
 * @returns {{t: Float64Array, count: Float64Array, slugs: {tail: number, front: number, length: number}[],
 *   arrivals: {t: number, length: number}[], merges: {t: number, s: number}[], dissipated: number,
 *   frontVelocity: number}} t [s]; number of slugs in the pipe at each step [-]; slugs in the pipe at tEnd [m];
 *   arrivals at the outlet (front arrival time [s], slug length then [m]); merge events; number of slugs that
 *   decayed [-]; front velocity [m/s]
 */
export function trackSlugs({ length, D, vm, vt, slugHoldup, filmHoldup, vLF, slugs = [], inletPeriod, inletSlugLength, wakeA = 0,
  wakeB = 1.06, tEnd, dt }) {
  checkArg('length', length, { min: 0, exclusiveMin: true });
  checkArg('D', D, { min: 0, exclusiveMin: true });
  checkArg('slugHoldup', slugHoldup, { min: 0, max: 1, exclusiveMin: true });
  checkArg('filmHoldup', filmHoldup, { min: 0, max: slugHoldup });
  if (filmHoldup >= slugHoldup) throw new RangeError('filmHoldup must be below slugHoldup');
  checkArg('tEnd', tEnd, { min: 0, exclusiveMin: true });
  checkArg('dt', dt, { min: 0, exclusiveMin: true });
  checkArg('wakeA', wakeA, { min: 0 });
  const vF = (slugHoldup * vm - filmHoldup * vLF) / (slugHoldup - filmHoldup);
  // Ordered from the outlet backwards: list[0] is the slug furthest downstream.
  let list = slugs.map((s) => {
    if (!(s.front > s.tail)) throw new RangeError('each initial slug needs front > tail');
    return { tail: s.tail, front: s.front, arrived: false };
  }).sort((a, b) => b.front - a.front);
  const n = Math.ceil(tEnd / dt - 1e-9);
  const ts = new Float64Array(n + 1);
  const count = new Float64Array(n + 1);
  const arrivals = [];
  const merges = [];
  let dissipated = 0;
  let nextInlet = inletPeriod ? inletPeriod : Infinity;
  if (inletPeriod !== undefined) { checkArg('inletPeriod', inletPeriod, { min: 0, exclusiveMin: true }); checkArg('inletSlugLength', inletSlugLength, { min: 0, exclusiveMin: true }); }
  count[0] = list.length;
  for (let k = 1; k <= n; k++) {
    const t = Math.min(tEnd, k * dt);
    const h = t - ts[k - 1];
    ts[k] = t;
    for (const s of list) {
      const Ls = s.front - s.tail;
      const tailV = vt * (1 + (wakeA > 0 ? wakeA * Math.exp(-wakeB * Ls / D) : 0));
      const f0 = s.front;
      s.front += vF * h;
      s.tail += tailV * h;
      if (!s.arrived && s.front >= length) {
        s.arrived = true;
        const tau = (length - f0) / (vF * h);
        arrivals.push({ t: ts[k - 1] + tau * h, length: Ls + (vF - tailV) * tau * h });
      }
    }
    // merge: front of an upstream slug reaches the tail of the one ahead
    for (let i = list.length - 1; i > 0; i--) {
      if (list[i].front >= list[i - 1].tail) {
        merges.push({ t, s: list[i - 1].tail });
        list[i - 1].tail = list[i].tail;
        list.splice(i, 1);
      }
    }
    const before = list.length;
    list = list.filter((s) => s.front > s.tail || s.tail >= length);
    dissipated += before - list.length;
    list = list.filter((s) => s.tail < length);
    while (t >= nextInlet) { list.push({ tail: 0, front: inletSlugLength, arrived: false }); nextInlet += inletPeriod; }
    count[k] = list.length;
  }
  return {
    t: ts, count, slugs: list.map((s) => ({ tail: s.tail, front: Math.min(s.front, length), length: Math.min(s.front, length) - s.tail })),
    arrivals, merges, dissipated, frontVelocity: vF,
  };
}
