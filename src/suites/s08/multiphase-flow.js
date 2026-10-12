// Steady multiphase flow and pressure drop along the route: coupled temperature and pressure profiles, phase
// fractions, velocities, liquid holdup and flow regime at every station, with the turndown curve of the line.
import { turndownCurve } from '../../engine/flow/multiphase.js';
import { regimeMap } from '../../engine/flow/flowregime.js';
import { pipeInputs, pipeFrom, routeInputs, routeFrom, routeLength, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, sev, worst, linspace,
  steelCo2, LIMIT_NO_CFD, ECON_BASIS, BAR, KELVIN } from './_common.js';
import { fluidInputs, fluidFrom, thermalInputs, thermalFrom, temperatureAlong, boundaryInputs, solveHydraulics } from './_flow.js';

const G_N = 'Hydraulic model';
const G_C = 'Acceptance criteria';
const METHODS = [['beggs-brill', 'Beggs & Brill with Payne corrections'], ['mechanistic', 'Mechanistic (Taitel–Dukler map, two-fluid stratified, slug unit cell)'], ['drift-flux', 'Drift flux (Bendiksen closure)'], ['homogeneous', 'Homogeneous (no slip)']];

export default {
  id: 's08-multiphase-flow',
  suite: 8,
  title: 'Multiphase flow and pressure drop',
  summary: 'Computes pressure, temperature, holdup, phase velocities and flow regime along the pipeline, and the inlet pressure or rate the line needs.',
  covers: [14],
  fidelity: 'numerical',
  standards: ['Beggs & Brill (1973) with Payne et al. (1979)', 'Taitel & Dukler (1976) flow-regime transitions', 'API RP 14E (1991, R2019) §2.5 erosional velocity', 'Brill & Mukherjee, Multiphase Flow in Wells (SPE Monograph 17)'],
  models: ['Conservation of mass, momentum and energy (steady, one-dimensional)', 'Beggs–Brill pressure-drop correlation', 'Drift-flux / Zuber–Findlay holdup', 'Homogeneous equilibrium model', 'Taitel–Dukler transition model', 'Darcy–Weisbach friction with Colebrook–White', 'Black-oil PVT correlations', 'Newton cooling with an overall heat-transfer coefficient'],
  limits: `Steady one-dimensional flow marched with a point pressure-gradient model and black-oil (or constant) properties; temperature is solved first from an exponential heat-loss balance without Joule–Thomson cooling, then the pressure. There is no compositional tracking, no mechanistic annular-flow model and no emulsion viscosity model. Transient behaviour is in the transient module. ${LIMIT_NO_CFD}`,
  publishes: ['pipeline.dP', 'pipeline.Tout', 'pipeline.pIn', 'pipeline.holdup', 'pipeline.vMix'],
  searchTerms: 'multiphase flow pressure drop Beggs Brill flowline holdup flow regime',
  inputs: [
    ...routeInputs(),
    ...pipeInputs({ omit: ['rhoCont', 'wsGiven', 'rhoW'] }),
    num('rough', 'Wall roughness', 'mm', 0.045, 0.001, 3, 'Pipe', 'geometry', 'Absolute roughness of the bore: 0.045 mm new carbon steel, 0.005 mm flow-coated, up to 0.5 mm corroded.', { uncertain: { dist: 'lognormal', cov: 0.4 } }),
    ...fluidInputs(),
    ...boundaryInputs(),
    ...thermalInputs(),
    sel('method', 'Pressure-gradient method', 'beggs-brill', METHODS, G_N, 'control', 'Beggs & Brill is the conventional choice; compare with the mechanistic and drift-flux methods to see the closure sensitivity.'),
    { key: 'nSeg', label: 'Integration segments along the route', type: 'int', default: 30, min: 8, max: 1000, group: G_N, role: 'control', help: 'Largest integration step = route length / segments (fourth-order Runge–Kutta, so few are needed).' },
    num('pInMax', 'Maximum allowable inlet pressure', 'bar', 200, 5, 1400, G_C, 'criteria', 'Design pressure or wellhead flowing pressure available at the inlet.'),
    num('pOutMin', 'Minimum arrival pressure', 'bar', 20, 1, 1000, G_C, 'criteria', 'Lowest pressure acceptable at the outlet (separator inlet).'),
    num('erosC', 'API RP 14E erosional constant C', '-', 100, 50, 300, G_C, 'criteria', '100 continuous / 125 intermittent solids-free service; up to 150–200 with corrosion-resistant alloy.'),
  ],
  mesh: { key: 'nSeg', label: 'Integration segments', levels: [10, 20, 40, 80], qoi: 'dP' },
  compute(inp, ctx) {
    const p = pipeFrom({ ...inp, rhoCont: 0 }, ctx);
    const route = routeFrom(inp), Ltot = routeLength(route);
    let fluid = fluidFrom(inp);
    const th = thermalFrom(inp, p);
    const solve = (fl, bracket) => {
      const T = temperatureAlong({ route, Tin: inp.Tin, Tamb: inp.Tamb, mDot: fl.massRate, cp: inp.cpMix, U: th.U, Dref: p.D });
      const h = solveHydraulics({ inp, route, D: p.Di, roughness: inp.rough / 1e3, fluid: fl, temperature: T.at, method: inp.method, maxStep: Ltot / inp.nSeg, extra: { erosionalC: inp.erosC }, bracket });
      return { T, ...h };
    };
    let sol = solve(fluid);
    // rate solved from the two pressures: the temperature profile depends on the rate, so repeat once with the solved rate
    if (!inp.useRate) { sol = solve(sol.fluid, [0.6, 1.6]); fluid = sol.fluid; }
    const r = sol.res, st = r.stations;
    const sIn = fluid.state(r.pIn, st[0].T, p.Di), sOut = fluid.state(r.pOut, st[st.length - 1].T, p.Di);
    const vmMax = Math.max(...st.map((q) => q.vm));
    const TarrC = st[st.length - 1].T - KELVIN;
    // regime transitions along the route
    const trans = []; let prev = null;
    st.forEach((q) => { if (q.regime !== prev) { trans.push({ kp: round(q.s / 1e3, 4), regime: q.regime, holdup: round(q.holdup, 3), vsl: round(q.vsl, 3), vsg: round(q.vsg, 3), angle: round(q.theta * 180 / Math.PI, 3) }); prev = q.regime; } });
    // turndown curve at the solved outlet pressure
    const m0 = r.massRate, fr = [0.2, 0.45, 0.7, 1, 1.5];
    let td = null;
    if (inp.method !== 'mechanistic') try {
      td = turndownCurve({ profile: route, D: p.Di, roughness: inp.rough / 1e3, fluidAt: fluid.at, temperature: sol.T.at, method: inp.method, maxStep: Ltot / 12, pOut: r.pOut, massRates: fr.map((f) => f * m0) });
    } catch { td = null; }
    const mid = st[Math.floor(st.length / 2)];
    let map = null;
    try {
      const fm = fluid.at(mid.P, mid.T);
      map = regimeMap({ model: 'taitel-dukler', D: p.Di, theta: 0, rhoL: fm.rhoL, rhoG: fm.rhoG, muL: fm.muL, muG: fm.muG, sigma: fm.sigma, nx: 20, ny: 16, vsgRange: [0.02, 60], vslRange: [0.005, 8], point: { vsg: Math.max(mid.vsg, 0.021), vsl: Math.max(mid.vsl, 0.0051) } });
    } catch { map = null; }
    const checks = [
      chk('Inlet pressure against the available pressure', r.pIn, inp.pInMax * BAR, 'Pa', 'Hydraulic design basis (design pressure / wellhead flowing pressure)',
        'Guards against the line needing more inlet pressure than the source or the design pressure can supply; if exceeded the rate cannot be delivered.'),
      chk('Minimum arrival pressure', inp.pOutMin * BAR, r.pOut, 'Pa', 'Hydraulic design basis (separator operating pressure)',
        'Guards against the arrival pressure falling below what the receiving facility needs.'),
      chk('Erosional velocity', r.maxErosionalRatio, 1, '-', 'API RP 14E (1991, R2019) §2.5.a',
        'Guards against the mixture velocity exceeding the empirical erosional limit C/√ρ, above which erosion-corrosion and noise become likely.'),
    ];
    const gov = worst(checks);
    const km = st.map((q) => round(q.s / 1e3, 4));
    const steelT = p.steelKgM * Ltot / 1e3;
    const fricShare = r.dP !== 0 ? r.dPFriction / Math.abs(r.dP) : 0;
    const unstable = td && td.hasInteriorMinimum && td.minimum.massRate > m0;
    return {
      kpis: [
        kpi('dP', 'Total pressure drop', r.dP / BAR, 'bar', 'Inlet minus outlet pressure: friction plus elevation plus acceleration. It sets the wellhead back-pressure or the boosting requirement.'),
        kpi('pIn', 'Inlet pressure', r.pIn / BAR, 'bar', inp.usePin ? 'Specified inlet pressure.' : 'Required inlet pressure to deliver the rate against the outlet pressure.', inp.pInMax, 'max'),
        kpi('pOut', 'Outlet pressure', r.pOut / BAR, 'bar', inp.usePout ? 'Specified outlet pressure.' : 'Arrival pressure for the stated inlet pressure and rate.', inp.pOutMin, 'min'),
        kpi('massRate', 'Total mass rate', r.massRate, 'kg/s', inp.useRate ? 'Mass flow rate of the stream (oil, gas and water).' : 'Mass rate the line delivers between the two stated pressures.'),
        kpi('dPfric', 'Frictional pressure drop', r.dPFriction / BAR, 'bar', 'Part of the pressure drop from wall friction; grows with the square of the rate.'),
        kpi('dPgrav', 'Elevation pressure drop', r.dPGravity / BAR, 'bar', 'Part from lifting the mixture up the route profile (negative for a net downhill line); dominates at low rate where liquid accumulates.'),
        kpi('gradient', 'Mean hydraulic gradient', r.dP / BAR / (Ltot / 1e3), 'bar/km', 'Average pressure loss per kilometre of pipe.'),
        kpi('Tarrival', 'Arrival temperature', TarrC, '°C', 'Fluid temperature at the outlet; compare with the hydrate and wax appearance temperatures.'),
        kpi('Uvalue', 'Overall heat-transfer coefficient', th.U, 'W/m²K', 'U-value on the steel outer diameter used for the temperature profile.'),
        kpi('holdup', 'Mean liquid holdup', r.meanHoldup, '-', 'Length-averaged fraction of the pipe cross-section occupied by liquid, including slip.'),
        kpi('inventory', 'Liquid inventory', r.liquidInventory, 'm³', 'Liquid volume standing in the line at steady state; the volume swept out by a pig or a rate increase.'),
        kpi('vmMax', 'Maximum mixture velocity', vmMax, 'm/s', 'Highest sum of the gas and liquid superficial velocities along the route (at the low-pressure end).'),
        kpi('erosRatio', 'Maximum velocity / erosional velocity', r.maxErosionalRatio, '-', 'Mixture velocity over the API RP 14E limit. Keep below 1.0 for continuous service.', 1, 'max'),
        kpi('gvfIn', 'Gas volume fraction at the inlet', 1 - sIn.lambdaL, '-', 'No-slip gas volume fraction at inlet conditions.'),
        kpi('gvfOut', 'Gas volume fraction at the outlet', 1 - sOut.lambdaL, '-', 'No-slip gas volume fraction at outlet conditions; it rises along the line as pressure falls and gas evolves.'),
        kpi('waterFrac', 'Water fraction of the liquid at the outlet', sOut.waterFraction, '-', 'In-situ water volume fraction of the liquid phase.'),
      ],
      checks,
      plots: [
        plot('line', 'Pressure and temperature along the route', 'Chainage', 'km', 'Pressure', 'bar', [
          { name: 'Pressure', x: km, y: st.map((q) => round(q.P / BAR)) }, { name: 'Temperature (°C)', x: km, y: st.map((q) => round(q.T - KELVIN)) },
        ], 'Pressure falls along the line by friction and elevation; temperature decays towards ambient. The arrival values at the right-hand end are the ones the facility sees.'),
        plot('line', 'Liquid holdup and elevation', 'Chainage', 'km', 'Holdup', '-', [
          { name: 'Liquid holdup', x: km, y: st.map((q) => round(q.holdup)) }, { name: 'No-slip liquid fraction', x: km, y: st.map((q) => round(q.lambdaL)) },
          { name: 'Elevation (÷ 1000 m)', x: km, y: st.map((q) => round(q.z / 1000, 4)), dash: true },
        ], 'Holdup rises in uphill sections where liquid slips back and falls downhill. The gap to the no-slip line is the liquid accumulated by slip.'),
        plot('line', 'Phase velocities', 'Chainage', 'km', 'Velocity', 'm/s', [
          { name: 'Superficial liquid', x: km, y: st.map((q) => round(q.vsl)) }, { name: 'Superficial gas', x: km, y: st.map((q) => round(q.vsg)) },
          { name: 'Actual liquid', x: km, y: st.map((q) => round(q.vL)) }, { name: 'Actual gas', x: km, y: st.map((q) => round(q.vG)) },
        ], 'Gas accelerates along the line as it expands. The actual phase velocities set residence time, heat loss and erosion.'),
        plot('line', 'Hydraulic gradient components', 'Chainage', 'km', 'Pressure loss per metre', 'Pa/m', [
          { name: 'Friction', x: km, y: st.map((q) => round(q.gradFriction)) }, { name: 'Gravity', x: km, y: st.map((q) => round(q.gradGravity)) }, { name: 'Total', x: km, y: st.map((q) => round(q.gradTotal)) },
        ], 'Local pressure gradient split into friction and elevation. Negative gravity terms are downhill sections recovering pressure.'),
        ...(td ? [plot('line', 'Turndown curve', 'Total mass rate', 'kg/s', 'Pressure drop', 'bar', [
          { name: 'Pressure drop', x: td.massRate.map((v) => round(v)), y: td.dP.map((v) => round(v / BAR)) }, { name: 'Friction part', x: td.massRate.map((v) => round(v)), y: td.dPFriction.map((v) => round(v / BAR)) },
          { name: 'Operating point', x: [round(m0)], y: [round(r.dP / BAR)] },
        ], 'Pressure drop against rate at the outlet pressure. To the left of the minimum the line is gravity dominated: liquid accumulates and operation is unstable and prone to terrain slugging.')] : []),
        ...(map ? [{ type: 'heatmap', title: 'Flow-regime map at mid-route conditions (horizontal)', x: { label: 'Superficial gas velocity', unit: 'm/s' }, y: { label: 'Superficial liquid velocity', unit: 'm/s' },
          xs: map.xs.map((v) => round(v)), ys: map.ys.map((v) => round(v)), z: map.z, zLabel: `Regime index: ${map.legend.map((l, i) => `${i} = ${l}`).join(', ')}`, zUnit: '-',
          explain: `Taitel–Dukler flow-regime map for the fluid at mid-route. The operating point there (v_sg ${fmt(mid.vsg)} m/s, v_sl ${fmt(mid.vsl)} m/s) is ${map.point ? map.point.regime : mid.regime}; its distance from a boundary shows how robust the regime is to rate changes.` }] : []),
      ],
      tables: [
        tbl('Flow regime along the route', [col('kp', 'From chainage', 'km'), col('regime', 'Flow regime'), col('holdup', 'Liquid holdup', '-'), col('vsl', 'v_sl', 'm/s'), col('vsg', 'v_sg', 'm/s'), col('angle', 'Inclination', 'deg')], trans.slice(0, 80),
          'Each row marks where the predicted flow regime changes along the route, with the local holdup, superficial velocities and pipe angle. Intermittent rows are where slugging loads should be checked.'),
        tbl('Phase split at inlet and outlet', [col('q', 'Quantity'), col('inlet', 'Inlet'), col('outlet', 'Outlet'), col('u', 'Unit')], [
          { q: 'Oil volumetric rate', inlet: round(sIn.qOil * 3600), outlet: round(sOut.qOil * 3600), u: 'm³/h' }, { q: 'Free gas volumetric rate', inlet: round(sIn.qGas * 3600), outlet: round(sOut.qGas * 3600), u: 'm³/h' },
          { q: 'Water volumetric rate', inlet: round(sIn.qWater * 3600), outlet: round(sOut.qWater * 3600), u: 'm³/h' }, { q: 'Liquid density', inlet: round(sIn.rhoL), outlet: round(sOut.rhoL), u: 'kg/m³' },
          { q: 'Gas density', inlet: round(sIn.rhoG), outlet: round(sOut.rhoG), u: 'kg/m³' }, { q: 'Liquid viscosity', inlet: round(sIn.muL * 1e3), outlet: round(sOut.muL * 1e3), u: 'cP' },
          { q: 'Bubble-point pressure', inlet: round((sIn.Pb || 0) / BAR), outlet: round((sOut.Pb || 0) / BAR), u: 'bar' },
        ], 'In-situ phase rates and properties at the two ends of the line: gas comes out of solution and expands towards the outlet.'),
      ],
      notes: [`Boundary set: ${[inp.usePin && 'inlet pressure', inp.usePout && 'outlet pressure', inp.useRate && 'flow rate'].filter(Boolean).join(' + ')}; solved for ${sol.solvedFor === 'massRate' ? 'the flow rate' : sol.solvedFor === 'pIn' ? 'the inlet pressure' : 'the outlet pressure'}.`,
        `Friction is ${fmt(fricShare * 100)} % of the pressure drop; U-value ${fmt(th.U)} W/m²K (${th.given ? 'typed' : 'built from the wall and coatings'}).`,
        r.multipleSolutions ? 'Two rates satisfy the stated pressures; the friction-dominated (higher) one is reported.' : (unstable ? 'The operating rate lies on the gravity-dominated branch of the turndown curve: expect liquid accumulation and terrain slugging.' : 'The operating rate lies on the friction-dominated (stable) branch of the turndown curve.'),
        ...((sIn.warnings || []).slice(0, 2).map((w) => `Fluid correlation range: ${w}`))],
      advice: [
        adv('engineering', gov.uc > 1 ? 'fail' : unstable ? 'warn' : sev(gov.uc), `Pressure drop ${fmt(r.dP / BAR)} bar over ${fmt(Ltot / 1e3)} km (inlet ${fmt(r.pIn / BAR)} bar, arrival ${fmt(r.pOut / BAR)} bar, ${fmt(TarrC)} °C); velocity is ${fmt(r.maxErosionalRatio * 100)} % of the erosional limit.`,
          gov.uc > 1 ? (gov.name.startsWith('Erosional') ? `Reduce the rate or go up one pipe size; the bore must grow by about ${fmt((Math.sqrt(r.maxErosionalRatio) - 1) * 100)} % to meet the velocity limit.` : 'Increase the diameter, lower the arrival pressure or add boosting; friction falls with roughly the fifth power of the bore.')
            : unstable ? `Operate above about ${fmt(td.minimum.massRate)} kg/s or plan for slug handling; below it the line is gravity dominated.` : 'The line delivers the rate within the pressure limits. Check the arrival temperature against hydrate and wax limits and repeat with a second pressure-drop method for closure sensitivity.',
          'Beggs & Brill (1973) / Payne et al. (1979); API RP 14E §2.5; flow-assurance practice per ISO 13628-1 / API RP 17A'),
        adv('economic', 'info', `Friction takes ${fmt(r.dPFriction / BAR)} bar of the ${fmt(r.dP / BAR)} bar pressure drop; the ${fmt(inp.D)} mm line contains ${fmt(steelT)} t of steel.`,
          fricShare > 0.7 && r.maxErosionalRatio > 0.6 ? 'The line is friction limited: the next pipe size up adds roughly 20 % steel but cuts friction by more than half, which may defer boosting.' : 'The line is not friction limited at this rate; a smaller bore would save steel but raise velocity and liquid surges — check the late-life rates before downsizing.',
          ECON_BASIS, `${steelCo2(steelT)} Back-pressure on the wells reduces recovery, so every bar of friction avoided is production gained without added energy.`),
      ],
      publish: { 'pipeline.dP': r.dP, 'pipeline.Tout': TarrC, 'pipeline.pIn': r.pIn, 'pipeline.holdup': r.meanHoldup, 'pipeline.vMix': vmMax },
      quantities: { pipe_km: round(Ltot / 1e3), steel_t: round(steelT) },
    };
  },
  verify: [{
    // Single-phase liquid limit on a horizontal line: Darcy–Weisbach with the Colebrook–White friction factor.
    // D_i = 0.3 m, ρ = 850 kg/m³, μ = 5 cP, 60 kg/s → v = 60/(850·π·0.15²) = 0.9986 m/s, Re = 850·0.9986·0.3/0.005 = 50 930.
    // Smooth pipe: 1/√f = −2 log10(2.51/(Re√f)) → f = 0.02081. ΔP = f (L/D) ρ v²/2 = 0.02081·(10 000/0.3)·850·0.9986²/2 = 2.940 bar.
    name: 'Single-phase liquid limit: Darcy–Weisbach with Colebrook–White',
    source: 'Colebrook (1939); Moody (1944); hand calculation',
    inputs: { D: 320, t: 10, rough: 0.001, fluidModel: 'constant', mRate: 60, xGas: 0, rhoLc: 850, muLc: 5, usePin: false, usePout: true, useRate: true, pOut: 30,
      profile: [{ kp: 0, z: -100 }, { kp: 10, z: -100 }], method: 'beggs-brill', nSeg: 40 },
    expect: [{ kpi: 'dP', value: 2.940, tol: 0.02 }, { kpi: 'holdup', value: 1, tol: 0.001 }],
  }, {
    // Static liquid column: negligible flow in a vertical rise of 100 m gives ΔP = ρ g h = 850·9.80665·100 = 8.336 bar.
    name: 'Hydrostatic limit: near-zero flow up a 100 m rise',
    source: 'Hydrostatics, ΔP = ρ g h',
    inputs: { D: 320, t: 10, fluidModel: 'constant', mRate: 0.05, xGas: 0, rhoLc: 850, muLc: 5, usePin: false, usePout: true, useRate: true, pOut: 30,
      profile: [{ kp: 0, z: -100 }, { kp: 1, z: 0 }], method: 'homogeneous', nSeg: 20 },
    expect: [{ kpi: 'dP', value: 8.336, tol: 0.003 }],
  }],
};
