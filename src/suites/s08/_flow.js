// Shared flow-assurance helpers: the produced-fluid description (black-oil stream or constant properties), the
// property callback handed to the hydraulic engines, the thermal build-up of the pipe wall and coatings, and the
// coupled steady pressure–temperature march used by the hydraulics, slugging, hydrate, wax and corrosion modules.
import { fluidAt } from '../../engine/flow/blackoil.js';
import { marchPipeline } from '../../engine/flow/multiphase.js';
import { overallU, externalFilmCoefficient, burialFilmCoefficient, temperatureProfile } from '../../engine/flow/thermalpipe.js';
import { num, sel, need, run, clamp, fmt, routeLength, KELVIN, BAR } from './_common.js';

const G_F = 'Produced fluid';
const G_T = 'Thermal design';
const DAY = 86400;

/** Fluid description inputs: black-oil stream by default, or constant phase properties. */
export function fluidInputs({ qo = 3000, GOR = 120, WC = 20 } = {}) {
  const bo = { showIf: { key: 'fluidModel', in: ['blackoil'] } }, cst = { showIf: { key: 'fluidModel', in: ['constant'] } };
  return [
    sel('fluidModel', 'Fluid model', 'blackoil', [['blackoil', 'Black oil (oil, gas and water from correlations)'], ['constant', 'Constant phase properties']], G_F, 'control',
      'Black oil derives phase rates and properties from stock-tank rates and gravities at every pressure and temperature. Constant properties is for single-phase or laboratory fluids.'),
    num('qo', 'Stock-tank oil rate', 'Sm³/d', qo, 0, 200000, G_F, 'bc', 'Oil rate at standard conditions. With the gas–oil ratio and water cut it defines all three phase rates.', { ...bo, design: true, uncertain: { dist: 'normal', cov: 0.1 } }),
    num('GOR', 'Produced gas–oil ratio', 'Sm³/Sm³', GOR, 0, 20000, G_F, 'bc', 'Total produced gas per stock-tank oil volume.', { ...bo, uncertain: { dist: 'lognormal', cov: 0.15 } }),
    num('WC', 'Water cut', '%', WC, 0, 98, G_F, 'bc', 'Water fraction of the produced liquid at stock-tank conditions.', { ...bo, uncertain: { dist: 'normal', cov: 0.2 } }),
    num('api', 'Oil gravity', '°API', 35, 10, 60, G_F, 'material', 'Stock-tank oil gravity.', bo),
    num('gammaG', 'Gas gravity (air = 1)', '-', 0.75, 0.55, 1.5, G_F, 'material', 'Specific gravity of the produced gas.', bo),
    num('salinity', 'Water salinity', 'g/kg', 35, 0, 250, G_F, 'material', 'Dissolved solids of the produced water.', bo),
    num('yCO2', 'CO₂ in the gas', 'mol %', 2, 0, 50, G_F, 'material', 'Mole fraction of carbon dioxide in the gas phase.', bo),
    num('yH2S', 'H₂S in the gas', 'mol %', 0, 0, 30, G_F, 'material', 'Mole fraction of hydrogen sulphide in the gas phase.', bo),
    num('mRate', 'Total mass rate', 'kg/s', 60, 0.01, 5000, G_F, 'bc', 'Total mass flow rate of the stream.', cst),
    num('xGas', 'Gas mass fraction', '-', 0, 0, 1, G_F, 'bc', 'Mass fraction of gas; zero gives single-phase liquid.', cst),
    num('rhoLc', 'Liquid density', 'kg/m³', 850, 300, 1300, G_F, 'material', 'Constant liquid density.', cst),
    num('muLc', 'Liquid viscosity', 'cP', 5, 0.05, 5000, G_F, 'material', 'Constant liquid dynamic viscosity.', cst),
    num('rhoGc', 'Gas density', 'kg/m³', 60, 0.5, 400, G_F, 'material', 'Constant gas density at line conditions.', cst),
    num('muGc', 'Gas viscosity', 'cP', 0.014, 0.005, 0.1, G_F, 'material', 'Constant gas viscosity.', cst),
    num('sigmac', 'Gas–liquid interfacial tension', 'mN/m', 20, 1, 80, G_F, 'material', 'Constant interfacial tension.', cst),
    num('cpMix', 'Stream specific heat', 'J/kgK', 2300, 800, 4500, G_F, 'material', 'Mass-averaged heat capacity for the temperature profile.'),
  ];
}

/**
 * Fluid record: property callback for the hydraulic engines and helpers.
 * @returns {{ at(P,T): object, state(P,T,D): object, massRate: number, black: boolean, scaled(f): object }}
 *   at → { rhoL, rhoG, muL, muG, sigma, quality }; state → full black-oil record (or the constant equivalent) with
 *   superficial velocities for inner diameter D; massRate in kg/s.
 */
export function fluidFrom(inp, rateFactor = 1) {
  const black = inp.fluidModel !== 'constant';
  if (black) {
    const qo = inp.qo / DAY * rateFactor;
    need(qo > 0, 'The oil rate must be positive for the black-oil stream. For a gas–water system use the constant-property model.');
    const base = { qo, GOR: inp.GOR, WC: Math.min(inp.WC / 100, 0.98), gammaG: inp.gammaG, api: inp.api, salinity: inp.salinity, yCO2: (inp.yCO2 ?? 0) / 100, yH2S: (inp.yH2S ?? 0) / 100 };
    const raw = (P, T, D) => fluidAt({ ...base, P: clamp(P, 1.2e5, 1.4e8), T: clamp(T, 255, 470), ...(D ? { D } : {}) });
    const ref = run(() => raw(50e5, 320), 'Fluid properties');
    const at = (P, T) => { const f = raw(P, T); return { rhoL: f.rhoL, rhoG: Math.max(f.rhoG, 0.5), muL: f.muL, muG: f.muG, sigma: Math.max(f.sigmaL, 1e-3), quality: clamp(f.mGas / f.mTotal, 0, 1) }; };
    return { black, at, state: raw, massRate: ref.mTotal, base, scaled: (f) => fluidFrom(inp, rateFactor * f) };
  }
  const m = inp.mRate * rateFactor, x = inp.xGas;
  const props = { rhoL: inp.rhoLc, rhoG: inp.rhoGc, muL: inp.muLc / 1e3, muG: inp.muGc / 1e3, sigma: inp.sigmac / 1e3, quality: x };
  const state = (P, T, D) => {
    const qL = m * (1 - x) / props.rhoL, qG = m * x / props.rhoG, A = D ? Math.PI * D * D / 4 : NaN;
    return { P, T, rhoL: props.rhoL, rhoG: props.rhoG, rhoO: props.rhoL, rhoW: props.rhoL, muL: props.muL, muG: props.muG, muO: props.muL, sigmaL: props.sigma, mTotal: m, mGas: m * x, mOil: m * (1 - x), mWater: 0,
      qOil: qL, qWater: 0, qGas: qG, qLiquid: qL, qTotal: qL + qG, lambdaL: qL / Math.max(qL + qG, 1e-30), waterFraction: 0, Z: 1, Pb: 0, Rs: 0, qgFreeSc: 0, vsl: qL / A, vsg: qG / A, vm: (qL + qG) / A, warnings: [] };
  };
  return { black, at: () => props, state, massRate: m, base: null, scaled: (f) => fluidFrom(inp, rateFactor * f) };
}

/** Thermal design inputs: coating conductivities, burial and the external flow, or a typed U-value. */
export function thermalInputs({ Tin = 70 } = {}) {
  return [
    num('Tin', 'Inlet temperature', '°C', Tin, -20, 200, 'Boundary conditions', 'bc', 'Fluid temperature at the pipeline inlet.', { uncertain: { dist: 'normal', cov: 0.05 } }),
    num('Tamb', 'Ambient sea-water temperature', '°C', 4, -2, 35, 'Boundary conditions', 'bc', 'Sea-water (or seabed) temperature along the route.', { bind: 'site.Tbed' }),
    num('Ugiven', 'Overall heat-transfer coefficient (0 = calculate)', 'W/m²K', 0, 0, 500, G_T, 'operation', 'U-value referred to the steel outer diameter. Zero builds it from the wall, coatings, burial and external flow below.'),
    num('kCoat', 'Coating / insulation conductivity', 'W/mK', 0.22, 0.02, 5, G_T, 'material', 'Solid polypropylene 0.22; PP foam / syntactic 0.15–0.17; FBE 0.3.'),
    num('kConc', 'Concrete conductivity', 'W/mK', 2.0, 0.5, 4, G_T, 'material', 'Thermal conductivity of the weight coating.'),
    num('burial', 'Burial depth to pipe axis (0 = exposed)', 'm', 0, 0, 10, G_T, 'geometry', 'Depth of the pipe centreline below the seabed; zero for a pipe lying on the seabed.'),
    num('kSoil', 'Soil conductivity', 'W/mK', 1.5, 0.5, 4, G_T, 'soil', 'Thermal conductivity of the backfill for a buried pipe.'),
    num('Usea', 'Sea-water velocity across the pipe', 'm/s', 0.2, 0.01, 3, G_T, 'environment', 'Sets the external film coefficient of an exposed pipe.'),
  ];
}

/**
 * Overall U on the steel outer diameter.
 * @param {object} inp thermal inputs @param {object} p pipe record from pipeFrom
 * @returns {{ U:number, hExt:number, Do:number, layers:object[] }} U [W/m²K on π·D_steel]
 */
export function thermalFrom(inp, p) {
  const layers = [{ t: p.t, k: 45, rho: p.rhoS, cp: 500 }];
  if ((inp.tCoat ?? 0) > 0) layers.push({ t: inp.tCoat / 1e3, k: inp.kCoat, rho: inp.rhoCoat ?? 900, cp: 2000 });
  if ((inp.tConc ?? 0) > 0) layers.push({ t: inp.tConc / 1e3, k: inp.kConc, rho: p.concRho, cp: 880 });
  const Do = p.Do;
  const hSea = externalFilmCoefficient({ v: Math.max(inp.Usea, 0.01), D: Do, rho: p.rhoW, mu: 1.5e-3, k: 0.58, cp: 3990 }).h;
  const hExt = inp.burial > Do / 2 ? burialFilmCoefficient({ D: Do, H: inp.burial, kSoil: inp.kSoil, hSea: 200 }).h : hSea;
  if (inp.Ugiven > 0) return { U: inp.Ugiven, hExt, Do, layers, given: true };
  const u = run(() => overallU({ Di: p.Di, layers, hi: 1500, ho: hExt, Dref: p.D }), 'Overall heat-transfer coefficient');
  return { U: u.U, hExt, Do, layers, given: false };
}

/** Resample a route to segments no longer than `ds` [m]. */
export function refineRoute(pts, ds) {
  const out = [pts[0]];
  for (let i = 1; i < pts.length; i++) {
    const a = pts[i - 1], b = pts[i], n = Math.max(1, Math.ceil(Math.hypot(b.x - a.x, b.z - a.z) / ds));
    for (let k = 1; k <= n; k++) out.push({ x: a.x + (b.x - a.x) * k / n, z: a.z + (b.z - a.z) * k / n });
  }
  return out;
}

/**
 * Steady temperature along the route (exponential decay per segment) as a function of chainage along the pipe.
 * @returns {{ s:number[], T:number[], at(s):number, Tarrival:number, heatLoss:number }} s [m], T [K]
 */
export function temperatureAlong({ route, Tin, Tamb, mDot, cp, U, Dref }) {
  const fine = refineRoute(route, Math.max(routeLength(route) / 200, 5));
  const r = run(() => temperatureProfile({ Tin: Tin + KELVIN, mDot, cp, U: Math.max(U, 1e-6), D: Dref, profile: fine, Tamb: Tamb + KELVIN }), 'Temperature profile');
  const s = [0];
  for (let i = 1; i < fine.length; i++) s.push(s[i - 1] + Math.hypot(fine[i].x - fine[i - 1].x, fine[i].z - fine[i - 1].z));
  const at = (q) => {
    if (q <= 0) return r.T[0];
    if (q >= s[s.length - 1]) return r.T[r.T.length - 1];
    let lo = 0, hi = s.length - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (s[m] <= q) lo = m; else hi = m; }
    const w = (q - s[lo]) / (s[hi] - s[lo]);
    return r.T[lo] * (1 - w) + r.T[hi] * w;
  };
  return { s, T: r.T, at, Tarrival: r.Tarrival, heatLoss: r.heatLoss };
}

/** Hydraulic boundary-condition inputs: any two of inlet pressure, outlet pressure and flow rate. */
export function boundaryInputs({ pIn = 120, pOut = 30 } = {}) {
  const G_B = 'Boundary conditions';
  return [
    { key: 'usePin', label: 'Inlet pressure is specified', type: 'bool', default: false, group: G_B, role: 'bc', help: 'Specify exactly two of inlet pressure, outlet pressure and flow rate; the third is solved.' },
    num('pIn', 'Inlet pressure', 'bar', pIn, 1, 1400, G_B, 'bc', 'Absolute pressure at the pipeline inlet (wellhead or manifold).', { showIf: { key: 'usePin', in: [true] } }),
    { key: 'usePout', label: 'Outlet pressure is specified', type: 'bool', default: true, group: G_B, role: 'bc', help: 'Arrival (separator or riser-top) pressure.' },
    num('pOut', 'Outlet pressure', 'bar', pOut, 1, 1400, G_B, 'bc', 'Absolute pressure at the pipeline outlet.', { showIf: { key: 'usePout', in: [true] } }),
    { key: 'useRate', label: 'Flow rate is specified', type: 'bool', default: true, group: G_B, role: 'bc', help: 'When off, the rate in the fluid section is only the starting guess and the rate is solved from the two pressures.' },
  ];
}

/**
 * Check the hydraulic boundary set and run the steady march.
 * @returns {{ res: object, fluid: object, solvedFor: string, rateFactor: number }}
 */
export function solveHydraulics({ inp, route, D, roughness, fluid, temperature, method, maxStep, extra = {}, bracket = [0.03, 8] }) {
  const n = [inp.usePin, inp.usePout, inp.useRate].filter(Boolean).length;
  if (n === 3) throw new Error('Over-constrained hydraulic boundary conditions: inlet pressure, outlet pressure and flow rate are all specified. Specify exactly two; the third follows from the pressure-drop solution.');
  if (n < 2) throw new Error(`Under-constrained hydraulic boundary conditions: ${n === 0 ? 'none' : 'only one'} of inlet pressure, outlet pressure and flow rate is specified. Specify exactly two.`);
  const common = { profile: route, D, roughness, fluidAt: fluid.at, temperature, method, maxStep, ...extra };
  if (inp.useRate) {
    const args = { ...common, massRate: fluid.massRate };
    if (inp.usePin) args.pIn = inp.pIn * BAR; else args.pOut = inp.pOut * BAR;
    const res = run(() => marchPipeline(args), 'Steady pressure march');
    need(res.pOut > 1.05e5 && res.pIn > 1.05e5, `The pipeline cannot deliver this rate: the pressure falls to atmospheric before the outlet (inlet ${fmt(res.pIn / BAR)} bar, outlet ${fmt(res.pOut / BAR)} bar). Raise the inlet pressure, lower the rate or increase the diameter.`);
    return { res, fluid, solvedFor: res.solvedFor, rateFactor: 1 };
  }
  need(inp.pIn > inp.pOut - 400, 'The inlet pressure is far below the outlet pressure: no forward flow is possible.');
  const m0 = fluid.massRate;
  // the root search runs on a coarse step (the fourth-order march is already converged there); the profile is then re-marched
  const coarse = run(() => marchPipeline({ ...common, maxStep: Math.max(maxStep, routeLength(route) / 10), pIn: inp.pIn * BAR, pOut: inp.pOut * BAR, rateBracket: [bracket[0] * m0, bracket[1] * m0] }), 'Flow-rate solution for the stated pressures (no rate between 3 % and 800 % of the starting guess satisfies them)');
  const f = coarse.massRate / m0;
  const res = run(() => marchPipeline({ ...common, massRate: coarse.massRate, pIn: inp.pIn * BAR }), 'Steady pressure march');
  res.multipleSolutions = coarse.multipleSolutions;
  return { res, fluid: fluid.scaled(f), solvedFor: 'massRate', rateFactor: f };
}
