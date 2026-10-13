// Shared helpers of the offshore renewable-energy suite: unit handling, advice/check records, screening
// sustainability factors, the turbine model (NREL 5-MW reference or a geometrically scaled rotor with its Cp/Ct
// surfaces, steady operating schedule and baseline controller), the monopile + tower support-structure description
// and the caveat texts that every module carries into its limits and advice.
// Engineering units are converted to SI here so that every module calls the engines in SI (docs/CONTRACT.md §2.1).
import { nrel5mw, scaleRotor } from '../../engine/aero/rotors.js';
import { cpSurface, maxCp, surfaceLookup } from '../../engine/aero/bem.js';
import { operatingSchedule, controllerFromRotor, piPitchGains } from '../../engine/aero/control.js';

export const G = 9.80665;
export const DEG = Math.PI / 180;
export const RPM = Math.PI / 30;
export const HOURS_YEAR = 8760;
// Screening emission factors (order-of-magnitude, cradle-to-gate): not project-specific EPD values.
export const CO2_STEEL = 1.9;      // t CO2e per t steel
export const CO2_CONCRETE = 0.3;   // t CO2e per m3 concrete
export const CO2_MGO = 3.206;      // t CO2 per t marine gas oil burnt
export const MGO_T_PER_DAY = 30;   // typical jack-up / heavy-lift installation vessel fuel burn, t/day (screening)
export const CO2_GRID = 0.4;       // t CO2 per MWh displaced from a gas-dominated grid (screening)
export const CO2_COPPER = 4.0;     // t CO2e per t refined copper (screening)
export const CO2_ALUMINIUM = 12;   // t CO2e per t primary aluminium (screening)
export const CO2_POLYMER = 2.5;    // t CO2e per t polyethylene (screening)
export const CO2_ROCK = 0.01;      // t CO2e per t quarried rock at the quay (screening)
export const CO2_GFRP = 3.0;       // t CO2e per t glass-fibre composite (screening)

export const need = (cond, msg) => { if (!cond) throw new Error(msg); };

/** Run an engine call and turn range/convergence errors into a message a user can act on. */
export function run(fn, hint) {
  try { return fn(); } catch (e) { throw new Error(`${hint}: ${e.message}`); }
}
export function finite(v, what) {
  if (!Number.isFinite(v)) throw new Error(`${what} could not be evaluated (non-finite result). Check that the inputs describe a physically possible system.`);
  return v;
}
export const round = (v, n = 4) => (Number.isFinite(v) ? Number(v.toPrecision(n)) : v);
export const fmt = (v, n = 3) => (Number.isFinite(v) ? String(Number(v.toPrecision(n))) : '–');
export const clamp = (v, a, b) => Math.min(b, Math.max(a, v));

/** Decimate parallel arrays to at most n points (keeps first and last). */
export function thin(arrs, n = 400) {
  const len = arrs[0].length;
  if (len <= n) return arrs.map((a) => Array.from(a));
  const idx = []; for (let i = 0; i < n; i++) idx.push(Math.round(i * (len - 1) / (n - 1)));
  return arrs.map((a) => idx.map((i) => a[i]));
}
export const scaled = (a, f) => Array.from(a, (v) => round(v * f, 5));
export const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));
export const interp = (xs, ys, x) => {
  if (x <= xs[0]) return ys[0];
  const n = xs.length;
  if (x >= xs[n - 1]) return ys[n - 1];
  let i = 1; while (xs[i] < x) i++;
  return ys[i - 1] + (ys[i] - ys[i - 1]) * (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
};

/* ------------------------------------------------------------------------------------------------ advice */

export const adv = (kind, severity, text, action, basis, sustainability) => (
  sustainability ? { kind, severity, text, action, basis, sustainability } : { kind, severity, text, action, basis });

export const check = (name, demand, capacity, unit, ref, explain) => {
  const uc = capacity > 0 ? Math.abs(demand) / capacity : Infinity;
  return { name, demand: round(demand, 5), capacity: round(capacity, 5), unit, uc: Number.isFinite(uc) ? round(uc, 5) : 99, ref, explain };
};
export const worst = (checks) => checks.reduce((a, c) => (c.uc > a.uc ? c : a), checks[0]);
export const sev = (uc, warnAt = 0.9) => (uc > 1 ? 'fail' : uc > warnAt ? 'warn' : 'ok');

export const vesselCo2 = (days, tPerDay = MGO_T_PER_DAY) => days * tPerDay * CO2_MGO;
export const vesselCo2Text = (days, tPerDay = MGO_T_PER_DAY) => `${fmt(days)} vessel day(s) ≈ ${fmt(days * tPerDay)} t marine gas oil ≈ ${fmt(vesselCo2(days, tPerDay))} t CO₂ (screening factors: ${tPerDay} t/day fuel, ${CO2_MGO} t CO₂ per t fuel).`;
export const steelCo2Text = (t) => `${fmt(t)} t of steel ≈ ${fmt(t * CO2_STEEL)} t CO₂e (screening factor ${CO2_STEEL} t CO₂e/t; use the mill EPD for the project).`;
export const energyCo2Text = (mwh) => `${fmt(mwh)} MWh/yr ≈ ${fmt(mwh * CO2_GRID)} t CO₂ a year displaced from (or made up by) gas-fired generation (screening factor ${CO2_GRID} t CO₂/MWh — use the grid factor of the project).`;

/* --------------------------------------------------------------------------------------------- caveat texts */

export const LIMIT_POLAR = 'Airfoil polars are generic parametric polars by relative thickness, and rotors other than the NREL 5-MW reference are geometric scalings of it (not published reference machines): performance and loads are preliminary and must be replaced by supplier data.';
export const LIMIT_COUPLED = 'The coupled simulator is reduced-order: rigid blades, a modal (linear, fore-aft) support structure, rotor-effective wind without blade-resolved turbulence or 1P/3P excitation, quasi-steady rotor aerodynamics from BEM surfaces and linear (Airy) waves without slamming or second-order forces.';
export const LIMIT_WAKE = 'Wake deficits come from engineering models (Jensen/Park top-hat, Bastankhah–Porté-Agel Gaussian, Frandsen/Crespo added turbulence): no CFD or large-eddy simulation, no dynamic wake meandering, flat uniform inflow per direction sector.';
export const LIMIT_WEC = 'Hydrodynamic coefficients are closed-form approximations for a truncated vertical cylinder in heave (or imported coefficient tables): there is no boundary-element (diffraction–radiation) solver in the module.';
export const LIMIT_TIDAL = 'Channel blockage uses one-dimensional actuator-disc (linear momentum) theory; the tide is a harmonic synthesis without nodal corrections and turbulence is not resolved.';
export const LIMIT_PV = 'Photovoltaic yield uses published transposition (Hay–Davies / isotropic), Erbs decomposition and single-diode models with the floating-specific thermal coefficients as inputs; no ray-traced shading or measured module-temperature data.';
/** Indicative unit-rate input (role 'cost'): never a quoted price, only a way to express a quantity as money in the advice. */
export const COST_HELP = 'Indicative placeholder — replace with current quotations.';
export const costInput = (key, label, unit, def, max, note = '') => ({ key, label, unit, type: 'number', default: def, min: 0, max, group: 'Indicative unit rates', role: 'cost', help: `${note ? `${note} ` : ''}${COST_HELP}` });
export const money = (v) => (Math.abs(v) >= 1e6 ? `${fmt(v / 1e6)} M USD` : `${fmt(v / 1e3)} k USD`);
export const LIMIT_CONTROL = 'Control is limited to the classical baseline loops named in the models; model-predictive control, Kalman filtering (state estimation) and state-space (LQR/LQG or other full-state feedback) control are not implemented.';
export const LIMIT_CODE = 'Code factors and loss assumptions are overridable inputs recalled from the documents and must be confirmed against the governing editions (IEC 61400-1, IEC 61400-3-1 / -3-2, IEC TS 62600, IEC 61853).';
export const CONFIRM = 'Overridable default recalled from the document — confirm against the governing edition.';
export const PXX = 'Pxx is an exceedance value: the annual energy that is exceeded with xx % probability (P90 < P50; P50 is the median estimate).';
export const SUPPLIER = 'Replace the generic rotor and polar data with the turbine supplier\'s certified power, thrust and load documents before design use.';
export const CFD_OUT = 'Three-dimensional CFD (RANS / LES / actuator-line), vortex-lattice aerodynamics and nonlinear finite elements are outside the module.';

/* -------------------------------------------------------------------------------------------- turbine model */

// Tip-speed-ratio and pitch grid of the Cp/Ct/Cq surfaces: dense around the optimum, sparse towards feather.
const TSR_GRID = [1, 2, 3, 4, 5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 10, 11, 12, 14];
const PITCH_GRID = [-2, -1, 0, 1, 2, 4, 6, 8, 10, 12, 15, 18, 21, 24, 27, 30, 40, 60, 90];
let baseCache = null;
const turbineCache = new Map();

/** Turbine selection inputs: the NREL 5-MW reference or a geometric scaling of it to another rating. */
export function turbineInputs({ group = 'Turbine', hub = true } = {}) {
  return [
    { key: 'turbine', label: 'Rotor model', type: 'select', default: 'nrel5', group, role: 'geometry',
      options: [{ value: 'nrel5', label: 'NREL 5-MW reference (126 m rotor)' }, { value: 'scaled', label: 'Generic rotor scaled from the NREL 5-MW reference' }],
      help: `Blade planform (span, chord, twist, airfoil thickness) of the NREL 5-MW reference (Jonkman et al. 2009) with generic parametric polars. ${SUPPLIER}` },
    { key: 'ratedMW', label: 'Rated power of the scaled rotor', unit: 'MW', type: 'number', default: 15, min: 1, max: 25, group, role: 'geometry', showIf: { key: 'turbine', in: ['scaled'] }, design: true,
      help: 'Lengths scale with √(P/5 MW) at constant tip speed and specific power; masses with the cube of the length scale (upper bound). Not the IEA 15-MW or DTU 10-MW reference machines.' },
    ...(hub ? [{ key: 'hubHeight', label: 'Hub height above still water', unit: 'm', type: 'number', default: 90, min: 30, max: 220, group, role: 'geometry',
      help: 'NREL 5-MW reference: 90 m. Must clear the blade tip above the highest wave crest (air gap). Scale with the rotor.' }] : []),
  ];
}

/**
 * Turbine record for the selected rotor: { rotor, surface, schedule, controller, scale }. The Cp/Ct/Cq surfaces
 * depend on tip-speed ratio and pitch only (Reynolds-independent polars), so the scaled rotor reuses the reference
 * surfaces with its own radius and swept area. Results are cached per rating (pure function of the inputs).
 */
export function getTurbine(inp, { rho = 1.225 } = {}) {
  const scaledSel = inp.turbine === 'scaled';
  const P = scaledSel ? inp.ratedMW * 1e6 : 5e6;
  need(P > 0, 'Rated power must be positive.');
  const key = `${scaledSel ? 's' : 'r'}:${P}:${rho}`;
  if (turbineCache.has(key)) return turbineCache.get(key);
  if (!baseCache) {
    const rotor = nrel5mw();
    baseCache = { rotor, surface: run(() => cpSurface({ rotor, tsr: TSR_GRID, pitchDeg: PITCH_GRID }), 'Rotor Cp/Ct surface') };
  }
  let rotor = baseCache.rotor, surface = baseCache.surface, scale = 1;
  if (scaledSel) {
    rotor = run(() => scaleRotor(baseCache.rotor, { ratedPower: P }), 'Rotor scaling');
    scale = rotor.scaleFactor;
    surface = { ...baseCache.surface, R: baseCache.surface.R * scale, area: baseCache.surface.area * scale * scale };
  }
  const schedule = run(() => operatingSchedule({ surface, ratedPower: rotor.ratedPower, ratedRotorSpeed: rotor.ratedRotorSpeed, minRotorSpeed: rotor.minRotorSpeed,
    efficiency: rotor.drivetrainEfficiency, gearboxRatio: rotor.gearboxRatio, rho, cutIn: rotor.cutIn, cutOut: rotor.cutOut }), 'Operating schedule');
  const controller = tunedController({ rotor, surface, schedule, rho });
  const iRated = schedule.thrust.indexOf(Math.max(...schedule.thrust));
  const rec = { rotor, surface, schedule, controller, scale, ratedPower: P, thrustRated: schedule.thrust[iRated], ratedWind: schedule.ratedWind, D: 2 * rotor.tipRadius };
  if (turbineCache.size > 12) turbineCache.clear();
  turbineCache.set(key, rec);
  return rec;
}

/** Turbine record for an arbitrary rotor definition (own Cp/Ct surfaces, not cached): used for user-entered blade geometry. */
export function buildTurbine(rotor, { rho = 1.225 } = {}) {
  const surface = run(() => cpSurface({ rotor, tsr: TSR_GRID, pitchDeg: PITCH_GRID }), 'Rotor Cp/Ct surface');
  const schedule = run(() => operatingSchedule({ surface, ratedPower: rotor.ratedPower, ratedRotorSpeed: rotor.ratedRotorSpeed, minRotorSpeed: rotor.minRotorSpeed,
    efficiency: rotor.drivetrainEfficiency, gearboxRatio: rotor.gearboxRatio, rho, cutIn: rotor.cutIn, cutOut: rotor.cutOut }), 'Operating schedule');
  need(Number.isFinite(schedule.ratedWind), `The rotor does not reach its rated power of ${fmt(rotor.ratedPower / 1e6)} MW below cut-out: reduce the rating or enlarge the rotor.`);
  const controller = tunedController({ rotor, surface, schedule, rho });
  const thrustRated = Math.max(...schedule.thrust);
  return { rotor, surface, schedule, controller, scale: rotor.tipRadius / 63, ratedPower: rotor.ratedPower, thrustRated, ratedWind: schedule.ratedWind, D: 2 * rotor.tipRadius };
}

/**
 * Baseline controller with pitch-loop gains designed on the region-3 pitch sensitivity of the rotor.
 * The engine's controllerFromRotor linearises dP/dθ at fine pitch, i.e. at the Cp optimum where the sensitivity
 * vanishes; with these polars that gives PI gains about eight times the NREL baseline and a sustained pitch limit
 * cycle just above rated wind speed. Here dP/dθ is evaluated along the above-rated schedule and represented, as in
 * Jonkman et al. (2009) §7.3, by dP/dθ(θ) = s0·(1 + θ/θK): s0 is floored at 20 % of the cut-out sensitivity and θK
 * matches the mid-schedule point. For the reference rotor this gives Kp ≈ 0.026 s, Ki ≈ 0.011, θK ≈ 7° (NREL
 * baseline with the measured polars: 0.0188 s, 0.0081, 6.3°).
 */
export function tunedController({ rotor, surface, schedule, rho = 1.225, omegaN = 0.6, zeta = 0.7 }) {
  const cfg = run(() => controllerFromRotor({ rotor, surface, rho, omegaN, zeta }), 'Baseline controller (the blade geometry must give power that falls as the blade pitches to feather, and the rating must be reachable)');
  const pts = [];
  for (let i = 0; i < schedule.U.length; i++) {
    if (schedule.region[i] !== '3') continue;
    const lk = surfaceLookup(surface, schedule.omega[i] * surface.R / schedule.U[i], schedule.pitchDeg[i]);
    pts.push([schedule.pitchDeg[i] * DEG, 0.5 * rho * surface.area * schedule.U[i] ** 3 * lk.dCpDPitchDeg / DEG]);
  }
  if (pts.length < 3) return cfg;
  const last = pts[pts.length - 1], mid = pts[Math.floor(pts.length / 2)];
  const s0 = 0.2 * last[1];
  if (!(s0 < 0) || !(mid[1] / s0 > 1.05)) return cfg;
  const g = piPitchGains({ inertia: rotor.drivetrainInertia, omegaRated: rotor.ratedRotorSpeed, dPdTheta: s0, gearboxRatio: rotor.gearboxRatio, zeta, omegaN });
  return { ...cfg, kp: g.kp, ki: g.ki, thetaK: mid[0] / (mid[1] / s0 - 1), design: { ...cfg.design, dPdTheta0: s0, zeta, omegaN } };
}

/** Maximum power coefficient and optimum tip-speed ratio / pitch of the reference planform (independent of scale and density). */
export function optimumCp() {
  getTurbine({ turbine: 'nrel5' });
  if (!baseCache.opt) baseCache.opt = run(() => maxCp({ rotor: baseCache.rotor }), 'Optimum power coefficient');
  return baseCache.opt;
}

/** Steady thrust [N], electrical power [W], rotor speed [rad/s] and pitch [deg] at hub wind speed U from the schedule (zero outside cut-in..cut-out). */
export function steadyAt(tb, U) {
  const s = tb.schedule;
  if (U < s.U[0] || U > s.U[s.U.length - 1]) return { thrust: 0, power: 0, omega: 0, pitchDeg: 90, ct: 0 };
  return { thrust: interp(s.U, s.thrust, U), power: interp(s.U, s.power, U), omega: interp(s.U, s.omega, U), pitchDeg: interp(s.U, s.pitchDeg, U), ct: interp(s.U, s.ct, U) };
}

/* ------------------------------------------------------------------------- monopile + tower support structure */

/** Tubular support-structure inputs (monopile to the transition, tapered tower above). Defaults: NREL 5-MW tower on a 6.5 m monopile in 30 m of water. */
export function supportInputs({ group = 'Support structure', depthBind = true, steel = true } = {}) {
  return [
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 30, min: 5, max: 70, group: 'Site', role: 'environment', ...(depthBind ? { bind: 'site.depth' } : {}),
      help: 'Still-water depth at the turbine. Monopiles are usual to about 50–60 m; beyond that the range of the form is exceeded and a jacket or a floater is the practical concept.' },
    { key: 'monoD', label: 'Monopile outer diameter', unit: 'm', type: 'number', default: 6.5, min: 2, max: 14, group, role: 'geometry', design: true, help: 'Constant from the mudline to the transition. 5-MW class: 5–7 m; 15-MW class: 9–11 m.' },
    { key: 'monoT', label: 'Monopile wall thickness', unit: 'mm', type: 'number', default: 80, min: 20, max: 150, group, role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.02 }, help: 'API RP 2A minimum for driving is 6.35 mm + D/100.' },
    { key: 'transZ', label: 'Tower-base (transition) elevation above still water', unit: 'm', type: 'number', default: 17, min: 3, max: 30, group, role: 'geometry', help: 'Interface level of tower and substructure: above the extreme wave crest plus surge, tide and air gap (15–22 m in the North Sea).' },
    { key: 'towerBaseD', label: 'Tower base diameter', unit: 'm', type: 'number', default: 6, min: 2, max: 12, group, role: 'geometry', design: true },
    { key: 'towerTopD', label: 'Tower top diameter', unit: 'm', type: 'number', default: 3.87, min: 1.5, max: 9, group, role: 'geometry' },
    { key: 'towerBaseT', label: 'Tower base wall thickness', unit: 'mm', type: 'number', default: 35, min: 10, max: 120, group, role: 'geometry', design: true, help: 'NREL 5-MW reference: 27 mm nominal, raised by 30 % here for the offshore monopile variant.' },
    { key: 'towerTopT', label: 'Tower top wall thickness', unit: 'mm', type: 'number', default: 25, min: 8, max: 80, group, role: 'geometry' },
    { key: 'towerMassFactor', label: 'Tower mass factor for flanges, internals and paint', unit: '-', type: 'number', default: 1.08, min: 1, max: 1.5, group, role: 'material', help: 'Multiplier on the shell mass of the tower (the NREL 5-MW definition uses an effective density of 8500 kg/m³ = 1.083 × steel).' },
    ...(steel ? [{ key: 'steel', label: 'Structural steel', type: 'select', source: 'materials:steel', default: 's355', group, role: 'material', help: 'Yield strength and modulus come from the materials library (nominal, t ≤ 16 mm: reduce for thick plate per EN 10025).' }] : []),
  ];
}

/** Segments for the beam model of coupled.js (z = 0 at still water, mudline at −depth). Tower density carries 8 % for flanges, paint and internals. */
export function supportSegments(inp, { zHub, E = 210e9, rhoSteel = 7850, hubOffset = 2.4 } = {}) {
  const zTop = zHub - hubOffset;   // tower top below the hub centre (NREL 5-MW: 2.4 m)
  need(inp.transZ < zTop - 10, `The transition elevation (${inp.transZ} m) must be at least 10 m below the tower top (${fmt(zTop)} m): raise the hub height or lower the transition.`);
  need(inp.monoT / 1e3 < inp.monoD / 4 && inp.towerBaseT / 1e3 < inp.towerBaseD / 4 && inp.towerTopT / 1e3 < inp.towerTopD / 4, 'Wall thickness must be below a quarter of the diameter (thin-walled tubes).');
  need(inp.towerTopD <= inp.towerBaseD, 'The tower top diameter cannot exceed the base diameter.');
  return [
    { zBase: -inp.depth, zTop: inp.transZ, dBase: inp.monoD, dTop: inp.monoD, tBase: inp.monoT / 1e3, tTop: inp.monoT / 1e3, density: rhoSteel, E },
    { zBase: inp.transZ, zTop, dBase: inp.towerBaseD, dTop: inp.towerTopD, tBase: inp.towerBaseT / 1e3, tTop: inp.towerTopT / 1e3, density: rhoSteel * (inp.towerMassFactor ?? 1.08), E },
  ];
}

/** Foundation boundary condition at the mudline: clamped or a coupled lateral / rocking spring matrix (bound to the pile-head stiffness of the foundation module). */
export function foundationInputs({ group = 'Foundation (mudline boundary)' } = {}) {
  const sp = { key: 'foundation', in: ['springs'] };
  return [
    { key: 'foundation', label: 'Mudline boundary condition', type: 'select', default: 'springs', group, role: 'bc', options: [{ value: 'springs', label: 'Coupled lateral / rocking springs (pile–soil)' }, { value: 'fixed', label: 'Clamped at the mudline' }],
      help: 'Pile–soil interaction represented at the mudline. A clamp over-predicts the natural frequency by 5–15 % for monopiles.' },
    { key: 'Kh', label: 'Lateral stiffness K_LL', unit: 'GN/m', type: 'number', default: 1.5, min: 0.01, max: 50, group, role: 'bc', bind: 'foundation.Kh', showIf: sp, uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'Kr', label: 'Rocking stiffness K_RR', unit: 'GN·m/rad', type: 'number', default: 120, min: 0.5, max: 5000, group, role: 'bc', bind: 'foundation.Kr', showIf: sp, uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'Khr', label: 'Coupling stiffness K_LR', unit: 'MN', type: 'number', default: -10000, min: -500000, max: 0, group, role: 'bc', bind: 'foundation.Khr', showIf: sp, help: 'Negative for a pile (a lateral head force rotates the head with the load). K_LL·K_RR must exceed K_LR².' },
  ];
}
export function foundationOf(inp) {
  if (inp.foundation !== 'springs') return 'fixed';
  const Kh = inp.Kh * 1e9, Kr = inp.Kr * 1e9, Khr = inp.Khr * 1e6;
  need(Kh * Kr > Khr * Khr * 1.02, `The foundation stiffness matrix is not positive definite (K_LL·K_RR = ${fmt(Kh * Kr)} must exceed K_LR² = ${fmt(Khr * Khr)}): reduce the coupling term or use the values of the foundation module.`);
  return { Kh, Kr, Khr };
}

/** Steel mass [kg] of the segments (frustum shell volume × density). */
export function segmentsMass(segs) {
  let m = 0;
  for (const s of segs) {
    const L = s.zTop - s.zBase, n = 20;
    for (let i = 0; i < n; i++) {
      const x = (i + 0.5) / n, D = s.dBase + (s.dTop - s.dBase) * x, t = s.tBase + (s.tTop - s.tBase) * x;
      m += Math.PI * (D - t) * t * (L / n) * s.density;
    }
  }
  return m;
}

/** Diameter and wall thickness [m] of the support structure at elevation z (still-water datum). */
export function sectionAt(segs, z) {
  const s = segs.find((q) => z >= q.zBase - 1e-9 && z <= q.zTop + 1e-9) || (z < segs[0].zBase ? segs[0] : segs[segs.length - 1]);
  const x = clamp((z - s.zBase) / (s.zTop - s.zBase), 0, 1);
  return { D: s.dBase + (s.dTop - s.dBase) * x, t: s.tBase + (s.tTop - s.tBase) * x };
}

/* ------------------------------------------------------------------------------------------- site inputs */

export const SITE = {
  A: { key: 'weibA', label: 'Weibull scale A at hub height', unit: 'm/s', type: 'number', default: 10.6, min: 4, max: 16, group: 'Wind resource', role: 'environment', bind: 'site.windWeibullA', uncertain: { dist: 'normal', cov: 0.05 }, help: 'Long-term 10-minute mean wind-speed distribution at hub height; mean = A·Γ(1 + 1/k).' },
  k: { key: 'weibK', label: 'Weibull shape k', unit: '-', type: 'number', default: 2.2, min: 1.2, max: 3.5, group: 'Wind resource', role: 'environment', bind: 'site.windWeibullK', help: 'Offshore North Sea sites: about 2.0–2.4.' },
  Hs1: { key: 'Hs1', label: '1-year significant wave height', unit: 'm', type: 'number', default: 5.5, min: 0.5, max: 14, group: 'Sea state', role: 'bc', bind: 'site.Hs1' },
  Tp1: { key: 'Tp1', label: '1-year peak period', unit: 's', type: 'number', default: 10, min: 3, max: 20, group: 'Sea state', role: 'bc', bind: 'site.Tp1' },
  Hs100: { key: 'Hs100', label: 'Extreme significant wave height (50/100-year)', unit: 'm', type: 'number', default: 8.5, min: 1, max: 18, group: 'Sea state', role: 'bc', bind: 'site.Hs100', uncertain: { dist: 'normal', cov: 0.08 }, help: 'IEC 61400-3-1 asks for the 50-year sea state; the 100-year value pulled from the site record is slightly conservative.' },
  Tp100: { key: 'Tp100', label: 'Extreme peak period', unit: 's', type: 'number', default: 12.5, min: 4, max: 22, group: 'Sea state', role: 'bc', bind: 'site.Tp100' },
  Uc1: { key: 'Uc1', label: '1-year surface current', unit: 'm/s', type: 'number', default: 0.8, min: 0, max: 5, group: 'Current', role: 'bc', bind: 'site.Uc1' },
  Uc100: { key: 'Uc100', label: 'Extreme surface current', unit: 'm/s', type: 'number', default: 1.2, min: 0, max: 6, group: 'Current', role: 'bc', bind: 'site.Uc100' },
};
export const site = (name, over = {}) => ({ ...SITE[name], ...over });
