// Hydrate assessment: equilibrium curve against the flowing pressure–temperature profile, subcooling and stability
// margin, thermodynamic inhibitor requirement, screening of formation, slurry loading and plugging risk, shutdown
// cooldown (no-touch) time and depressurisation.
import { hydrateTemperature, hydratePressure, INHIBITORS, hammerschmidtDepression, nielsenBucklinDepression, massToMoleFraction, saltDepression, inhibitorInjection,
  subcoolingProfile, inductionTime, hydrateFractionAlongPipe, slurryViscosityKriegerDougherty, pluggingRiskIndex, depressurisation } from '../../engine/flow/hydrate.js';
import { radialTransient, timeToTemperature, cooldownLumped } from '../../engine/flow/thermalpipe.js';
import { pipeInputs, pipeFrom, routeInputs, routeFrom, routeLength, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, sev, worst, linspace, interp, clamp,
  LIMIT_NO_CFD, ECON_BASIS, CO2_MGO, BAR, KELVIN } from './_common.js';
import { fluidInputs, fluidFrom, thermalInputs, thermalFrom, temperatureAlong, boundaryInputs, solveHydraulics } from './_flow.js';

const G_H = 'Hydrate equilibrium';
const G_I = 'Inhibition';
const G_K = 'Formation kinetics (user constants)';
const G_S = 'Shutdown';
const G_C = 'Acceptance criteria';

export default {
  id: 's08-hydrate',
  suite: 8,
  title: 'Hydrate formation, inhibition and cooldown',
  summary: 'Shows where the line enters the hydrate region, how much inhibitor keeps it out, how long it can stay shut in before hydrates form, and a screening view of plugging risk.',
  covers: [15],
  fidelity: 'screening',
  standards: ['Katz (1945) gas-gravity hydrate chart as fitted by Motiee (1991) and Towler & Mokhatab (2005)', 'Hammerschmidt (1939); Nielsen & Bucklin (1983) inhibitor depression', 'GPSA Engineering Data Book §20', 'Sloan & Koh, Clathrate Hydrates of Natural Gases (2008)', 'ISO 13628-1 / API RP 17A (flow assurance)'],
  models: ['Gas-gravity hydrate equilibrium correlations (or a user equilibrium curve)', 'Inhibitor thermodynamics (Hammerschmidt, Nielsen–Bucklin, colligative salt effect)', 'Hydrate nucleation induction time (classical nucleation theory form)', 'Subcooling-driven hydrate growth bookkeeping and slurry viscosity (Krieger–Dougherty)', 'Hydrate plugging-risk index (screening heuristic)', 'Transient radial heat conduction (shutdown cooldown)'],
  limits: `Hydrate curves are gas-gravity correlations for sweet gas (±1–3 K, worse with CO₂, H₂S or N₂); the engine holds a van der Waals–Platteeuw framework but no parameter set, so no statistical-thermodynamic equilibrium is computed — supply a laboratory or simulator curve to override. Induction-time and growth constants are user arguments with illustrative defaults, and the plugging-risk index is a screening heuristic, not a validated criterion. Inhibitor losses to the vapour and condensate phases are not included. Hydrate deposition, adhesion and agglomeration models and CFD–DEM are outside this module. ${LIMIT_NO_CFD}`,
  publishes: ['pipeline.hydrateMargin', 'pipeline.cooldownTime'],
  searchTerms: 'gas hydrate formation flowline subcooling MEG methanol injection cooldown time',
  inputs: [
    ...routeInputs(),
    ...pipeInputs({ omit: ['rhoCont', 'wsGiven', 'rhoW'], defaults: { tCoat: 50, rhoCoat: 750, tConc: 0 } }),
    num('rough', 'Wall roughness', 'mm', 0.045, 0.001, 3, 'Pipe', 'geometry', 'Absolute roughness of the bore.'),
    ...fluidInputs(),
    ...boundaryInputs({ pOut: 40 }),
    ...thermalInputs({ Tin: 65 }).map((i) => (i.key === 'kCoat' ? { ...i, default: 0.17 } : i)),
    sel('hydMethod', 'Hydrate curve', 'motiee', [['motiee', 'Motiee (1991) gas-gravity fit'], ['towler', 'Towler & Mokhatab (2005)'], ['hammerschmidt', 'Hammerschmidt (1934), no gravity dependence']], G_H, 'control', 'All are fits of the Katz gas-gravity chart for sweet gas with free water.'),
    { key: 'hydCurve', label: 'Measured or simulated hydrate curve (optional)', type: 'series', default: { x: [], y: [] }, xLabel: 'Pressure', xUnit: 'bar', yLabel: 'Hydrate temperature', yUnit: '°C', group: G_H, role: 'material',
      help: 'At least three points with increasing pressure replace the correlation (use a PVT laboratory or compositional simulator curve for design).' },
    sel('inhibitor', 'Thermodynamic inhibitor', 'MEG', Object.keys(INHIBITORS).map((k) => [k, k]), G_I, 'operation', 'Inhibitor injected at the inlet.'),
    num('wInh', 'Inhibitor in the water phase (present dosing)', 'wt %', 0, 0, 80, G_I, 'ic', 'Concentration currently achieved in the aqueous phase; zero for an uninhibited line.'),
    sel('inhMethod', 'Depression model', 'hammerschmidt', [['hammerschmidt', 'Hammerschmidt'], ['nielsenBucklin', 'Nielsen–Bucklin']], G_I, 'control', 'Hammerschmidt up to about 25 wt % methanol / 60 wt % glycol; Nielsen–Bucklin for concentrated solutions.'),
    num('purity', 'Inhibitor purity of the injected solution', '-', 0.9, 0.5, 1, G_I, 'operation', 'Mass fraction of inhibitor in the injected stream (lean MEG ≈ 0.9).'),
    num('kForm', 'Formation rate constant', '1/(s·K)', 2e-5, 0, 0.01, G_K, 'material', 'Lumped first-order conversion constant dX/dt = k·ΔT_sub·(1 − X). Illustrative default — fit to flow-loop data.'),
    num('dTonset', 'Subcooling needed for nucleation', 'K', 3, 0, 15, G_K, 'material', 'Formation is assumed to start only beyond this subcooling. Illustrative default.'),
    num('Kind', 'Induction-time pre-exponential K', 's', 60, 0.01, 1e7, G_K, 'material', 't_ind = K·exp(B/(T·ΔT²)). System specific; illustrative default.'),
    num('Bind', 'Induction-time constant B', 'K³', 5e4, 0, 1e8, G_K, 'material', 'System specific; illustrative default.'),
    num('pShut', 'Settle-out pressure after shutdown', 'bar', 80, 2, 1400, G_S, 'ic', 'Pressure in the line once shut in; it sets the hydrate temperature to which the line may cool.'),
    num('tSim', 'Cooldown time simulated', 'h', 36, 1, 240, G_S, 'control', 'Length of the shutdown cooldown history.'),
    num('pBlow', 'Blow-down pressure', 'bar', 8, 1.1, 300, G_S, 'operation', 'Pressure reached by depressurising the line to melt or avoid hydrates.'),
    num('margin', 'Required temperature margin to the hydrate curve', 'K', 3, 0, 15, G_C, 'criteria', 'Design margin covering the uncertainty of the curve and of the temperature prediction.'),
    num('tNoTouch', 'Required cooldown (no-touch) time', 'h', 8, 0, 72, G_C, 'criteria', 'Time the operators need before inhibitor injection or depressurisation must start.'),
    { key: 'nSeg', label: 'Integration segments along the route', type: 'int', default: 24, min: 8, max: 400, group: 'Numerics', role: 'control', help: 'Steps of the pressure march and resolution of the subcooling profile.' },
  ],
  mesh: { key: 'nSeg', label: 'Integration segments', levels: [12, 24, 48], qoi: 'subMax' },
  compute(inp, ctx) {
    const p = pipeFrom({ ...inp, rhoCont: 0 }, ctx);
    const route = routeFrom(inp), Ltot = routeLength(route);
    const fluid0 = fluidFrom(inp), th = thermalFrom(inp, p);
    const Tprof0 = temperatureAlong({ route, Tin: inp.Tin, Tamb: inp.Tamb, mDot: fluid0.massRate, cp: inp.cpMix, U: th.U, Dref: p.D });
    const sol = solveHydraulics({ inp, route, D: p.Di, roughness: inp.rough / 1e3, fluid: fluid0, temperature: Tprof0.at, method: 'beggs-brill', maxStep: Ltot / inp.nSeg });
    const fluid = sol.fluid, st = sol.res.stations;
    // equilibrium curve
    const hc = inp.hydCurve && inp.hydCurve.x && inp.hydCurve.x.length >= 3 ? inp.hydCurve : null;
    if (hc) for (let i = 1; i < hc.x.length; i++) need(hc.x[i] > hc.x[i - 1], 'The hydrate curve pressures must increase.');
    const gG = inp.gammaG ?? 0.7;
    const Teq = (P) => (hc ? interp(hc.x.map((v) => v * BAR), hc.y.map((v) => v + KELVIN), P) : hydrateTemperature({ P: clamp(P, 2e5, 4e8), gammaG: clamp(gG, 0.55, 1.0), method: inp.hydMethod }));
    const inh = INHIBITORS[inp.inhibitor];
    const dTinh = inp.wInh > 0 ? (inp.inhMethod === 'hammerschmidt' ? hammerschmidtDepression({ w: inp.wInh, inhibitor: inp.inhibitor }) : nielsenBucklinDepression(massToMoleFraction(inp.wInh / 100, inh.M))) : 0;
    const dTsalt = inp.fluidModel !== 'constant' && inp.salinity > 0 ? saltDepression({ w: inp.salinity / 10 }) : 0;
    const prof = run(() => subcoolingProfile({ s: st.map((q) => q.s), P: st.map((q) => q.P), T: st.map((q) => q.T), inhibitorDT: dTinh + dTsalt, curve: Teq }), 'Subcooling profile');
    const subMax = prof.maxSubcooling;
    // inhibitor requirement to hold the margin with no inhibitor credit
    const noInh = subcoolingProfile({ s: st.map((q) => q.s), P: st.map((q) => q.P), T: st.map((q) => q.T), inhibitorDT: dTsalt, curve: Teq });
    const dTreq = Math.max(noInh.maxSubcooling + inp.margin, 0);
    const sOut = fluid.state(sol.res.pOut, st[st.length - 1].T, p.Di);
    const waterRate = Math.max(sOut.mWater, 0);
    let inj = null;
    if (dTreq > 0 && waterRate > 0) { try { inj = inhibitorInjection({ dT: dTreq, inhibitor: inp.inhibitor, method: inp.inhMethod, waterRate, purity: inp.purity }); } catch { inj = null; } }
    // formation bookkeeping along the line
    const vL = st.reduce((a, q) => a + q.vL, 0) / st.length;
    const wf = clamp(sOut.waterFraction, 0, 1);
    const sub = prof.subcooling;
    const form = wf > 0 && inp.kForm > 0 && vL > 0 ? hydrateFractionAlongPipe({ x: st.map((q) => q.s), subcooling: sub, v: Math.max(vL, 0.01), k: inp.kForm, waterFraction: wf, dTonset: inp.dTonset }) : { phi: st.map(() => 0), conversion: st.map(() => 0), onset: null };
    const phiEnd = form.phi[form.phi.length - 1];
    const muRel = slurryViscosityKriegerDougherty({ phi: Math.min(phiEnd, 0.6) });
    const risk = pluggingRiskIndex({ subcooling: Math.max(subMax, 0), phi: phiEnd, velocity: Math.max(vL, 0.01) });
    const iMax = sub.indexOf(Math.max(...sub));
    const tInd = subMax > 0 ? inductionTime({ K: inp.Kind, B: inp.Bind, T: st[iMax].T, dT: subMax }) : Infinity;
    const tRes = prof.onset != null ? (Ltot - prof.onset) / Math.max(vL, 0.01) : 0;
    // shutdown cooldown at the coldest flowing point (the outlet end)
    const iCold = st.reduce((a, q, i) => (q.T < st[a].T ? i : a), 0);
    const Thyd = Teq(inp.pShut * BAR) - dTinh - dTsalt;
    const fS = fluid.at(inp.pShut * BAR, st[iCold].T), hold = clamp(sol.res.meanHoldup, 0.02, 1);
    const rhoC = hold * fS.rhoL + (1 - hold) * fS.rhoG, cpC = inp.cpMix;
    const cd = run(() => radialTransient({ Di: p.Di, contents: { rho: rhoC, cp: cpC }, layers: th.layers, hi: 200, ho: th.hExt, Tamb: inp.Tamb + KELVIN, T0: st[iCold].T, mode: 'cooldown', tEnd: inp.tSim * 3600, dt: Math.max(inp.tSim * 3600 / 400, 30), cellsPerLayer: 3, theta: 1 }), 'Cooldown');
    const tCool = timeToTemperature({ t: cd.t, T: cd.Tcontents, threshold: Thyd });
    const cdIn = radialTransient({ Di: p.Di, contents: { rho: rhoC, cp: cpC }, layers: th.layers, hi: 200, ho: th.hExt, Tamb: inp.Tamb + KELVIN, T0: st[0].T, mode: 'cooldown', tEnd: inp.tSim * 3600, dt: Math.max(inp.tSim * 3600 / 400, 30), cellsPerLayer: 3, theta: 1 });
    const tCoolH = Number.isFinite(tCool) ? tCool / 3600 : inp.tSim;
    const never = Thyd <= inp.Tamb + KELVIN;
    const dep = depressurisation({ P: inp.pBlow * BAR, Tamb: inp.Tamb + KELVIN, gammaG: clamp(gG, 0.55, 1), method: inp.hydMethod, inhibitorDT: dTinh + dTsalt });
    const marginK = prof.margin;
    const checks = [
      chk('Flowing temperature margin to the hydrate curve', inp.margin, Math.max(marginK, 1e-6), 'K', 'Flow-assurance design basis (ISO 13628-1 / API RP 17A)',
        'Guards against the flowing fluid entering the hydrate stability region anywhere along the route, where free water and gas can form a hydrate plug. A ratio above 1.0 means the required margin is not held.'),
      chk('Shutdown cooldown (no-touch) time', inp.tNoTouch, never ? 1e6 : Math.max(tCoolH, 1e-6), 'h', 'Flow-assurance design basis; transient conduction',
        'Guards against the shut-in line cooling into the hydrate region before the operators can inject inhibitor or depressurise.'),
    ];
    const gov = worst(checks);
    const km = st.map((q) => round(q.s / 1e3, 4));
    const ps = linspace(5, Math.max(1.3 * sol.res.pIn / BAR, inp.pShut * 1.2, 60), 30);
    const hT = cd.t.map((t) => round(t / 3600, 3));
    const dayMass = inj ? inj.injectionRate * 86400 / 1e3 : 0;
    const inside = prof.inHydrateRegion;
    return {
      kpis: [
        kpi('margin', 'Hydrate stability margin while flowing', marginK, 'K', 'Smallest distance between the fluid temperature and the (inhibited) hydrate temperature along the route. Negative means the line runs inside the hydrate region.', inp.margin, 'min'),
        kpi('subMax', 'Maximum subcooling', subMax, 'K', 'Largest amount by which the fluid is colder than the hydrate temperature; positive values drive formation.'),
        kpi('onset', 'Chainage where the hydrate region is entered', prof.onset != null ? prof.onset / 1e3 : Ltot / 1e3, 'km', 'First point along the route at which the fluid is colder than the hydrate temperature. The route length means it is never entered.'),
        kpi('TeqShut', 'Hydrate temperature at settle-out pressure', Teq(inp.pShut * BAR) - KELVIN, '°C', 'Equilibrium temperature of the uninhibited fluid at the shut-in pressure.'),
        kpi('dTinhib', 'Depression from the present inhibitor dose', dTinh, 'K', 'Lowering of the hydrate temperature by the inhibitor concentration now in the water phase.'),
        kpi('dTsalt', 'Depression from produced-water salinity', dTsalt, 'K', 'Lowering of the hydrate temperature by dissolved salt (ideal-solution estimate).'),
        kpi('wReq', 'Inhibitor concentration required in the water phase', inj ? inj.w * 100 : 0, 'wt %', 'Concentration that keeps the whole route outside the hydrate region with the design margin. Zero when no inhibitor is needed.'),
        kpi('injRate', 'Inhibitor solution injection rate', inj ? inj.injectionRate * 3600 : 0, 'kg/h', 'Mass rate of the injected solution at the stated purity; losses to vapour and condensate are not included.'),
        kpi('cooldown', 'Cooldown time to the hydrate temperature', never ? inp.tSim : tCoolH, 'h', 'Time for the shut-in contents at the coldest point to cool to the hydrate temperature at settle-out pressure. Equal to the simulated time when it is not reached.', inp.tNoTouch, 'min'),
        kpi('phiHyd', 'Hydrate volume fraction at the outlet (screening)', phiEnd, '-', 'Hydrate fraction of the liquid slurry from the lumped conversion model with the user rate constant.'),
        kpi('muRel', 'Slurry relative viscosity', Number.isFinite(muRel) ? muRel : 999, '-', 'Viscosity of the hydrate slurry over the carrier liquid (hard spheres); 999 denotes a non-flowing slurry.'),
        kpi('risk', 'Plugging-risk index (screening heuristic)', risk.index, '-', 'Heuristic 0–1 index combining subcooling, hydrate loading and velocity. Below 0.2 low, 0.2–0.5 moderate, above 0.5 high; not a validated criterion.'),
        kpi('tInd', 'Induction time at the coldest point', Number.isFinite(tInd) ? Math.min(tInd / 3600, 1e6) : 1e6, 'h', 'Characteristic nucleation delay from the user constants; 1e6 means no driving force. Compare with the residence time inside the hydrate region.'),
        kpi('pDiss', 'Dissociation pressure at ambient temperature', dep.Pdissociation / BAR, 'bar', 'Pressure below which hydrates are unstable at sea-bed temperature; depressurise below it to melt a plug.'),
      ],
      checks,
      plots: [
        plot('line', 'Operating path against the hydrate curve', 'Temperature', '°C', 'Pressure', 'bar', [
          { name: 'Hydrate curve (no inhibitor)', x: ps.map((q) => round(Teq(q * BAR) - KELVIN)), y: ps.map((q) => round(q)) },
          { name: 'Hydrate curve with inhibitor and salt', x: ps.map((q) => round(Teq(q * BAR) - KELVIN - dTinh - dTsalt)), y: ps.map((q) => round(q)), dash: true },
          { name: 'Flowing P–T path (inlet to outlet)', x: st.map((q) => round(q.T - KELVIN)), y: st.map((q) => round(q.P / BAR)) },
          { name: 'Shut-in at ambient', x: [inp.Tamb], y: [inp.pShut] },
        ], 'Hydrates are stable to the left of the curve. The flowing path should stay to the right with the design margin; the shut-in point shows where the line ends up if left to cool.'),
        plot('line', 'Subcooling along the route', 'Chainage', 'km', 'Subcooling', 'K', [
          { name: 'Subcooling (positive = inside the hydrate region)', x: km, y: sub.map((v) => round(v)) }, { name: 'Required margin', x: [0, round(Ltot / 1e3)], y: [-inp.margin, -inp.margin], dash: true },
          { name: 'Hydrate fraction × 100', x: km, y: form.phi.map((v) => round(v * 100)) },
        ], 'Subcooling grows as the fluid cools along the line. Where it crosses zero the hydrate region is entered; it must stay below the dashed line to hold the margin.'),
        plot('line', 'Shutdown cooldown', 'Time since shutdown', 'h', 'Contents temperature', '°C', [
          { name: 'Coldest point of the line', x: hT, y: cd.Tcontents.map((v) => round(v - KELVIN)) }, { name: 'Inlet end', x: cdIn.t.map((t) => round(t / 3600, 3)), y: cdIn.Tcontents.map((v) => round(v - KELVIN)) },
          { name: 'Hydrate temperature at settle-out pressure', x: [0, inp.tSim], y: [round(Thyd - KELVIN), round(Thyd - KELVIN)], dash: true },
        ], 'Cooling of the shut-in contents. The time at which the curve crosses the dashed hydrate temperature is the no-touch time available.'),
      ],
      tables: [tbl('Inhibition and remediation', [col('q', 'Quantity'), col('v', 'Value'), col('u', 'Unit')], [
        { q: 'Depression required (maximum subcooling + margin)', v: round(dTreq), u: 'K' }, { q: `${inp.inhibitor} required in the water phase`, v: inj ? round(inj.w * 100) : 0, u: 'wt %' },
        { q: 'Free water rate at the outlet', v: round(waterRate * 3600), u: 'kg/h' }, { q: 'Injected solution', v: round(dayMass), u: 't/day' },
        { q: 'Residence time inside the hydrate region', v: round(tRes / 3600, 3), u: 'h' }, { q: 'Plugging-risk class (heuristic)', v: risk.class, u: '' },
        { q: 'Hydrate temperature at the blow-down pressure', v: round(dep.Teq - KELVIN), u: '°C' }, { q: 'Ice risk during blow-down', v: dep.iceRisk ? 'yes — stay above ' + fmt(dep.PiceLimit / BAR) + ' bar' : 'no', u: '' },
        { q: 'Lumped cooldown time constant', v: round(cooldownLumped({ T0: st[iCold].T, Tamb: inp.Tamb + KELVIN, U: th.U, D: p.D, heatCapacity: rhoC * cpC * p.Ai + p.rhoS * 500 * p.As }).tau / 3600), u: 'h' },
      ], 'Inhibitor mass balance and remediation data: the dose that protects the line, the time the fluid spends in the hydrate region, and what depressurisation achieves.')],
      notes: [`Hydrate curve: ${hc ? 'user-supplied' : `${inp.hydMethod} correlation, gas gravity ${fmt(gG)}`}; U-value ${fmt(th.U)} W/m²K; arrival ${fmt(st[st.length - 1].T - KELVIN)} °C at ${fmt(sol.res.pOut / BAR)} bar.`,
        inside ? `The line enters the hydrate region at ${fmt(prof.onset / 1e3)} km with up to ${fmt(subMax)} K of subcooling.` : `The line stays ${fmt(marginK)} K outside the hydrate region while flowing.`,
        never ? 'At settle-out pressure the hydrate temperature is below ambient: the shut-in line never enters the hydrate region.' : (Number.isFinite(tCool) ? `After shutdown the coldest point reaches the hydrate temperature in ${fmt(tCoolH)} h.` : 'The hydrate temperature is not reached within the simulated cooldown time.'),
        'Kinetic constants are illustrative user inputs and the plugging-risk index is a screening heuristic.'],
      advice: [
        adv('engineering', gov.uc > 1 ? 'fail' : risk.index >= 0.5 ? 'warn' : sev(gov.uc), `Flowing margin to the hydrate curve ${fmt(marginK)} K (required ${fmt(inp.margin)} K); cooldown time ${fmt(never ? inp.tSim : tCoolH)} h (required ${fmt(inp.tNoTouch)} h).`,
          gov.uc > 1 ? (gov.name.startsWith('Flowing') ? (inj ? `Inject ${fmt(inj.injectionRate * 3600)} kg/h of ${fmt(inp.purity * 100)} % ${inp.inhibitor} to reach ${fmt(inj.w * 100)} wt % in the water phase, or add insulation to raise the arrival temperature by ${fmt(inp.margin - marginK)} K.` : `Raise the arrival temperature by ${fmt(inp.margin - marginK)} K with insulation or a higher rate.`)
            : `Increase the insulation or bury the line to extend the cooldown from ${fmt(tCoolH)} to ${fmt(inp.tNoTouch)} h, or plan automatic inhibitor injection or blow-down to below ${fmt(dep.Pdissociation / BAR)} bar on shutdown.`)
            : 'The line is protected while flowing and for the required no-touch time. Confirm the hydrate curve with a compositional simulator or laboratory data, and re-check at turndown rates where arrival temperature falls.',
          'ISO 13628-1 / API RP 17A flow-assurance practice; GPSA Engineering Data Book §20; Sloan & Koh (2008)'),
        adv('economic', 'info', inj ? `Continuous inhibition needs ${fmt(dayMass)} t/day of ${inp.inhibitor} solution (${fmt(dayMass * 365)} t/year circulating or consumed).` : 'No continuous inhibitor injection is needed at this operating point.',
          inj ? 'Compare the inhibitor supply, regeneration and umbilical capacity with the one-off cost of more insulation; insulation that lifts the arrival temperature above the curve removes the operating cost for the life of the field.' : 'Operating cost is limited to shutdown and restart inhibitor; size storage for the settle-out volume of the line rather than for continuous dosing.',
          ECON_BASIS, `MEG regeneration and methanol make-up carry an energy and emissions burden for every tonne injected (${fmt(dayMass * 365)} t/year here); supply by vessel is about ${CO2_MGO} t CO₂ per tonne of fuel burned. Passive insulation avoids both and reduces the chemical discharge risk.`),
      ],
      publish: { 'pipeline.hydrateMargin': marginK, 'pipeline.cooldownTime': (never ? inp.tSim : tCoolH) * 3600 },
    };
  },
  verify: [{
    // Hammerschmidt (1939): ΔT(°F) = K_H w/(M (100 − w)) with K_H = 2335 for methanol, M = 32.04 g/mol, w = 20 wt %:
    // ΔT = 2335·20/(32.04·80) = 18.22 °F = 10.12 K.
    // Hammerschmidt (1934) hydrate curve: T(°F) = 8.9 P^0.285 with P = 1000 psia (68.948 bar): 8.9·7.161 = 63.74 °F = 17.63 °C.
    name: 'Hammerschmidt depression for 20 wt % methanol and hydrate temperature at 1000 psia',
    source: 'Hammerschmidt (1934, 1939) as given in the GPSA Engineering Data Book §20; hand calculation',
    inputs: { inhibitor: 'methanol', wInh: 20, inhMethod: 'hammerschmidt', hydMethod: 'hammerschmidt', pShut: 68.948 },
    expect: [{ kpi: 'dTinhib', value: 10.12, tol: 0.003 }, { kpi: 'TeqShut', value: 17.63, abs: 0.1 }],
  }],
};
