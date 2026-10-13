// Transient gas–liquid flow in a flowline and riser: production ramp-up, shutdown and restart, and riser-induced
// slugging with choking or gas lift, by a one-dimensional isothermal drift-flux finite-volume solver.
import { flowlineRiserProfile, rampUp, shutdownRestart, severeSluggingCase, analyseOscillation, detectSlugArrivals } from '../../engine/flow/transient.js';
import { pipeInputs, pipeFrom, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, sev, worst, thin, arr,
  LIMIT_NO_CFD, ECON_BASIS, CO2_STEEL, BAR, KELVIN, G } from './_common.js';

const G_R = 'Flowline and riser';
const G_F = 'Fluids';
const G_S = 'Scenario';
const G_N = 'Numerics';
const R_GAS = 8.314462618;

export default {
  id: 's08-transient-flow',
  suite: 8,
  title: 'Transient multiphase flow: ramp-up, shutdown, restart and riser slugging',
  summary: 'Simulates how pressure, holdup and outlet rates respond in time to a rate change, a shutdown and restart, or riser-base liquid blockage.',
  covers: [14, 15],
  fidelity: 'numerical',
  standards: ['Evje & Fjelde (2002, 2003) drift-flux model and AUSMV-type flux splitting', 'Zuber & Findlay (1965) drift-flux slip relation', 'Bøe (1981); Taitel (1986) severe-slugging criteria', 'ISO 13628-1 / API RP 17A (flow assurance, transient operations)'],
  models: ['One-dimensional drift-flux model (two mass equations and a mixture momentum equation)', 'Zuber–Findlay slip law', 'Finite-volume solver with flux-vector splitting and MUSCL reconstruction', 'Slug front propagation as holdup waves; slug arrival detection', 'Terrain-induced severe slugging with choking and gas-lift mitigation'],
  limits: `Transient multiphase flow is a one-dimensional isothermal drift-flux model with constant liquid properties, an ideal-gas law with a fixed compressibility factor and no interphase mass transfer; the liquid sound speed is a numerical parameter (lower values run faster but make the liquid more compressible). Hydrodynamic slugs shorter than a few cells are not resolved, and there is no thermal transient, hydrate or wax coupling. Full two-fluid six/seven-equation models, slug tracking with film dynamics and 3-D CFD are outside this module. ${LIMIT_NO_CFD}`,
  publishes: ['pipeline.surgeVolume', 'pipeline.transientPmax'],
  searchTerms: 'transient multiphase flow drift flux ramp-up shutdown restart severe slugging simulation',
  inputs: [
    ...pipeInputs({ omit: ['steel', 'tCoat', 'rhoCoat', 'tConc', 'concrete', 'rhoCont', 'rhoW', 'wsGiven'] }),
    num('rough', 'Wall roughness', 'mm', 0.045, 0.001, 3, 'Pipe', 'geometry', 'Absolute roughness of the bore.'),
    num('Lflow', 'Flowline length', 'm', 900, 50, 50000, G_R, 'geometry', 'Length of the flowline up to the riser base. Long lines need many cells and long run times; choose the cell count so that a cell boundary falls at the riser base.'),
    num('angle', 'Flowline inclination (negative = downhill to the riser)', 'deg', -1, -10, 10, G_R, 'geometry', 'Constant slope of the flowline.'),
    num('hRiser', 'Riser height', 'm', 100, 1, 2000, G_R, 'geometry', 'Vertical riser height to the outlet.', { bind: 'site.depth' }),
    num('rhoL', 'Liquid density', 'kg/m³', 820, 400, 1300, G_F, 'material', 'Density of the liquid phase (oil–water mixture).'),
    num('muL', 'Liquid viscosity', 'cP', 2, 0.1, 2000, G_F, 'material', 'Viscosity of the liquid phase.'),
    num('Mgas', 'Gas molar mass', 'kg/kmol', 20, 16, 60, G_F, 'material', 'Molar mass of the gas (0.69 gravity ≈ 20 kg/kmol).'),
    num('Zgas', 'Gas compressibility factor', '-', 0.92, 0.5, 1.2, G_F, 'material', 'Constant Z at the mean line pressure and temperature.'),
    num('Tline', 'Line temperature', '°C', 40, -10, 150, G_F, 'bc', 'Uniform fluid temperature (isothermal model).'),
    sel('scenario', 'Scenario', 'rampup', [['rampup', 'Production ramp-up'], ['shutdown', 'Shutdown and restart'], ['severe', 'Riser-induced slugging at constant rate']], G_S, 'operation', 'Each scenario sets the inlet boundary history and the initial state.'),
    num('mL', 'Liquid mass rate', 'kg/s', 20, 0.05, 2000, G_S, 'bc', 'Inlet liquid rate at full production (inlet boundary: both phase rates are prescribed).', { design: true }),
    num('mG', 'Gas mass rate', 'kg/s', 1.0, 0.001, 500, G_S, 'bc', 'Inlet gas rate at full production.'),
    num('pSep', 'Outlet (separator) pressure', 'bar', 20, 1.5, 300, G_S, 'bc', 'Pressure boundary at the riser top (outlet boundary: one condition).'),
    num('frac0', 'Initial rate as a fraction of full rate', '-', 0.4, 0.05, 1, G_S, 'ic', 'Ramp-up starts from the steady state at this fraction of the full rates.', { showIf: { key: 'scenario', in: ['rampup'] } }),
    num('tRamp', 'Ramp duration', 's', 60, 1, 7200, G_S, 'operation', 'Time over which the inlet rates are raised linearly to full rate.', { showIf: { key: 'scenario', in: ['rampup'] } }),
    num('tShut', 'Time of shutdown', 's', 40, 1, 86400, G_S, 'operation', 'When the inlet starts closing; the line starts from steady production.', { showIf: { key: 'scenario', in: ['shutdown'] } }),
    num('tRestart', 'Time of restart', 's', 150, 20, 172800, G_S, 'operation', 'When the inlet starts re-opening; must be after the shutdown plus the valve time.', { showIf: { key: 'scenario', in: ['shutdown'] } }),
    num('choke', 'Outlet choke opening (0 = no choke)', '-', 0, 0, 1, G_S, 'bc', 'Cd × choke area / pipe area for the riser-slugging scenario; small openings stabilise the flow at the cost of back-pressure.', { showIf: { key: 'scenario', in: ['severe'] } }),
    num('gasLift', 'Riser-base gas-lift rate', 'kg/s', 0, 0, 100, G_S, 'bc', 'Gas injected at the riser base to prevent liquid blockage.', { showIf: { key: 'scenario', in: ['severe'] } }),
    num('tEnd', 'Simulated time', 's', 200, 20, 86400, G_N, 'control', 'Length of the simulation. Riser-slugging cycles need several cycle periods.'),
    { key: 'nCells', label: 'Finite-volume cells', type: 'int', default: 40, min: 12, max: 800, group: G_N, role: 'control', help: 'Cells along the flowline and riser; cost grows with the square of the count.' },
    num('aL', 'Liquid sound speed (numerical)', 'm/s', 300, 50, 1500, G_N, 'control', 'Real liquids are near 1000–1300 m/s; a lower value relaxes the acoustic time-step limit. Keep it at least ten times the flow velocity.'),
    num('C0', 'Drift-flux distribution coefficient', '-', 1.2, 1, 1.3, G_N, 'operation', 'Slip closure v_G = C0·v_m + v_d; 1.2 for slug and bubbly flow, 1.0 for no slip.'),
    num('pDesign', 'Design pressure', 'bar', 100, 5, 1400, 'Acceptance criteria', 'criteria', 'Limit for the highest pressure reached anywhere in the line during the transient.'),
    num('vSurge', 'Available surge volume', 'm³', 10, 0.1, 5000, 'Acceptance criteria', 'criteria', 'Liquid surge the receiving separator can absorb.'),
  ],
  mesh: { key: 'nCells', label: 'Finite-volume cells', levels: [20, 40, 80], qoi: 'pInMean' },
  compute(inp, ctx) {
    const p = pipeFrom({ ...inp, steel: 'api5l-x65' }, ctx);
    const D = p.Di, A = Math.PI * D * D / 4;
    const T = inp.Tline + KELVIN;
    const gasRT = inp.Zgas * R_GAS * T / (inp.Mgas / 1e3);
    const fluid = { rhoL: inp.rhoL, aL: inp.aL, gasRT, muL: inp.muL / 1e3, muG: 1.4e-5 };
    const profile = flowlineRiserProfile({ flowlineLength: inp.Lflow, flowlineAngleDeg: inp.angle, riserHeight: inp.hRiser });
    const base = { profile, D, roughness: inp.rough / 1e3, nCells: inp.nCells, tEnd: inp.tEnd, fluid, slip: { C0: inp.C0 }, nProfiles: 20, nHistory: 400, maxSteps: 4e5 };
    const pOut = inp.pSep * BAR;
    need(inp.pSep * BAR > 1.2e5, 'The outlet pressure must be above atmospheric.');
    let r, extra = {};
    if (inp.scenario === 'rampup') {
      need(inp.tRamp < inp.tEnd, 'The ramp must finish within the simulated time.');
      r = run(() => rampUp({ ...base, pOut, from: { massRateL: inp.frac0 * inp.mL, massRateG: inp.frac0 * inp.mG }, to: { massRateL: inp.mL, massRateG: inp.mG }, rampStart: 0.05 * inp.tEnd, rampDuration: inp.tRamp }), 'Ramp-up simulation');
      extra = { surge: r.surge, swept: r.liquidSwept };
    } else if (inp.scenario === 'shutdown') {
      need(inp.tRestart > inp.tShut + 12, 'Restart must come at least 12 s after the shutdown starts (valve closing time 10 s).');
      need(inp.tRestart < inp.tEnd, 'The restart must happen within the simulated time.');
      r = run(() => shutdownRestart({ ...base, pOut, massRateL: inp.mL, massRateG: inp.mG, tShutdown: inp.tShut, tRestart: inp.tRestart, valveTime: 10 }), 'Shutdown / restart simulation');
      extra = { surge: r.restartSurge, settled: r.settledInventory };
    } else {
      const { profile: _p, ...b2 } = base;
      r = run(() => severeSluggingCase({ ...b2, flowlineLength: inp.Lflow, flowlineAngleDeg: inp.angle, riserHeight: inp.hRiser, massRateL: inp.mL, massRateG: inp.mG, pSep: pOut,
        ...(inp.choke > 0 ? { chokeOpening: inp.choke } : {}), gasLiftRate: inp.gasLift }), 'Riser slugging simulation');
      extra = { surge: detectSlugArrivals({ t: r.history.t, rate: r.history.outletLiquidRate, discardFraction: 0.3 }), slugging: r.slugging, period: r.period, head: r.headAmplitude, boe: r.boe, taitel: r.taitel };
    }
    const h = r.history, c = r.conservation;
    const oscP = analyseOscillation({ t: h.t, y: h.inletPressure, discardFraction: 0.5 });
    let pMax = 0, qLmax = 0; for (let i = 0; i < h.t.length; i++) { pMax = Math.max(pMax, h.maxPressure[i]); qLmax = Math.max(qLmax, h.outletLiquidRate[i]); }
    const nH = h.t.length, tail = Math.max(1, Math.floor(nH / 10));
    let qLend = 0; for (let i = nH - tail; i < nH; i++) qLend += h.outletLiquidRate[i] / tail;
    const surge = extra.surge || { count: 0, maxAmount: 0, meanPeriod: NaN, peakRate: qLmax, meanRate: qLend };
    // surge volume: liquid delivered above the full-rate drain capacity during the largest surge
    let surgeVol = 0, run0 = 0; const dtH = h.t[1] - h.t[0];
    for (let i = 0; i < nH; i++) { const ex = (h.outletLiquidRate[i] - inp.mL) / inp.rhoL * dtH; run0 = ex > 0 ? run0 + ex : 0; surgeVol = Math.max(surgeVol, run0); }
    const inv0 = h.liquidInventory[0], invE = h.liquidInventory[nH - 1];
    // below 1e-9 of the mass in the pipe the imbalance is arithmetic round-off: reported as zero, otherwise to three figures
    const resRaw = Math.max(Math.abs(c.residualL), Math.abs(c.residualG)), residual = resRaw < 1e-9 ? 0 : Number(resRaw.toPrecision(3));
    const checks = [
      chk('Peak pressure during the transient', pMax, inp.pDesign * BAR, 'Pa', 'Design pressure (DNV-ST-F101 §3.4 pressure control)',
        'Guards against the transient pushing the pressure anywhere in the line above the design pressure, for example when gas compresses behind a liquid blockage.'),
      chk('Liquid surge against the separator capacity', surgeVol, inp.vSurge, 'm³', 'Liquid volume balance at the receiving separator',
        'Guards against the liquid arriving faster than the full-rate drain capacity for long enough to overfill the separator.'),
      chk('Mass conservation of the solution', residual, 1e-6, '-', 'Finite-volume conservation (numerical verification)',
        'Confirms that the liquid and gas masses in the pipe change only by what crosses the boundaries; a large residual would mean the numerical solution cannot be trusted.'),
    ];
    const gov = worst(checks);
    const [tt, qL, qG, pInH, invH] = thin([h.t, h.outletLiquidRate, h.outletGasRate, h.inletPressure, h.liquidInventory], 400);
    const s = r.grid.s, pr = r.profiles, nP = pr.t.length;
    const snaps = [0, Math.floor(nP / 3), Math.floor(2 * nP / 3), nP - 1].filter((v, i, a) => a.indexOf(v) === i);
    const slugging = inp.scenario === 'severe' ? !!extra.slugging : (oscP.relativeAmplitude > 0.1 && oscP.cycles >= 2);
    return {
      kpis: [
        kpi('pInMean', 'Mean inlet pressure (second half of the run)', oscP.mean / BAR, 'bar', 'Time-averaged pressure at the flowline inlet after the initial transient; the back-pressure the wells see.'),
        kpi('pInSwing', 'Inlet pressure swing', (oscP.max - oscP.min) / BAR, 'bar', 'Peak-to-peak variation of the inlet pressure over the second half of the run. A large persistent swing indicates slugging.'),
        kpi('pMax', 'Highest pressure in the line', pMax / BAR, 'bar', 'Maximum pressure reached in any cell at any time.', inp.pDesign, 'max'),
        kpi('qLpeak', 'Peak outlet liquid rate', qLmax, 'kg/s', 'Largest interval-mean liquid mass rate leaving the riser top.'),
        kpi('qLend', 'Outlet liquid rate at the end', qLend, 'kg/s', 'Mean liquid rate over the last tenth of the run; it equals the inlet rate once a steady state is reached.'),
        kpi('surgeVol', 'Liquid surge volume above full-rate drain', surgeVol, 'm³', 'Largest volume of liquid delivered in excess of the full-rate liquid capacity in one continuous surge.', inp.vSurge, 'max'),
        kpi('surges', 'Liquid surges detected', surge.count ?? 0, '-', 'Number of periods in which the outlet liquid rate exceeds 1.5 times its mean.'),
        kpi('period', 'Oscillation period', Number.isFinite(oscP.period) ? oscP.period : 0, 's', 'Mean period of the inlet pressure oscillation; zero when the flow is steady.'),
        kpi('slugging', 'Slugging detected (1 = yes)', slugging ? 1 : 0, '-', 'One when the run shows sustained riser-base pressure cycling; zero for stable flow.'),
        kpi('inv0', 'Liquid inventory at the start', inv0, 'm³', 'Liquid volume in the line in the initial state.'),
        kpi('invEnd', 'Liquid inventory at the end', invE, 'm³', 'Liquid volume in the line at the end of the run; a ramp-up sweeps the difference out as a surge.'),
        kpi('massResidual', 'Mass conservation residual', residual, '-', 'Relative imbalance between the change of mass in the pipe and the net boundary flow. It should vanish for a correct run; an imbalance below 1e-9 of the mass is arithmetic round-off and is reported as 0.', 1e-6, 'max'),
        kpi('steps', 'Time steps taken', r.steps, '-', 'Number of explicit time steps; the step is limited by the acoustic Courant condition.'),
        kpi('dtMin', 'Smallest time step', r.dtMin * 1e3, 'ms', 'Smallest explicit time step used in the run; it shrinks where gas-free liquid makes the mixture sound speed high.'),
      ],
      checks,
      plots: [
        plot('line', 'Outlet mass rates', 'Time', 's', 'Mass rate', 'kg/s', [
          { name: 'Liquid', x: tt.map((v) => round(v)), y: qL.map((v) => round(v)) }, { name: 'Gas × 10', x: tt.map((v) => round(v)), y: qG.map((v) => round(10 * v)) },
          { name: 'Full-rate liquid', x: [0, round(inp.tEnd)], y: [inp.mL, inp.mL], dash: true },
        ], 'Liquid and gas leaving the riser top. Spikes above the dashed line are surges the separator must absorb; gaps are periods of no liquid production.'),
        plot('line', 'Inlet pressure', 'Time', 's', 'Pressure', 'bar', [{ name: 'Inlet pressure', x: tt.map((v) => round(v)), y: pInH.map((v) => round(v / BAR)) }],
          'Pressure at the flowline inlet. A slow sawtooth is the signature of riser-base blockage; a smooth approach to a constant value is a stable transient.'),
        plot('line', 'Liquid holdup along the line', 'Chainage', 'm', 'Liquid holdup', '-',
          snaps.map((i) => ({ name: `t = ${fmt(pr.t[i])} s`, x: arr(s), y: arr(pr.holdup[i]) })),
          'Holdup profiles at four times. Liquid collecting at the riser base (right of the flowline) and travelling fronts show how the inventory redistributes.'),
        plot('line', 'Liquid inventory', 'Time', 's', 'Liquid in the line', 'm³', [{ name: 'Inventory', x: tt.map((v) => round(v)), y: invH.map((v) => round(v)) }],
          'Total liquid volume in the flowline and riser. A fall is liquid swept out; a rise is accumulation that will come out later as a surge.'),
      ],
      tables: [tbl('Run summary', [col('q', 'Quantity'), col('v', 'Value'), col('u', 'Unit')], [
        { q: 'Scenario', v: inp.scenario, u: '' }, { q: 'Inlet boundary', v: 'liquid and gas mass rates prescribed', u: '' }, { q: 'Outlet boundary', v: inp.scenario === 'severe' && inp.choke > 0 ? 'choke to separator pressure' : 'pressure prescribed', u: '' },
        { q: 'Initial condition', v: inp.scenario === 'rampup' ? `steady state at ${fmt(inp.frac0 * 100)} % rate` : 'steady state at full rate', u: '' },
        { q: 'Cell size', v: round(r.grid.ds), u: 'm' }, { q: 'Liquid mass in / out balance residual', v: Math.abs(c.residualL) < 1e-9 ? 0 : round(c.residualL, 3), u: '-' }, { q: 'Gas mass balance residual', v: Math.abs(c.residualG) < 1e-9 ? 0 : round(c.residualG, 3), u: '-' },
        { q: 'Mass removed by the positivity limiter (liquid / gas)', v: `${fmt(c.clippedL)} / ${fmt(c.clippedG)}`, u: 'kg' },
        ...(extra.settled != null ? [{ q: 'Liquid inventory just before restart', v: round(extra.settled), u: 'm³' }] : []),
        ...(extra.swept != null ? [{ q: 'Liquid swept out by the ramp-up', v: round(extra.swept), u: 'm³' }] : []),
        ...(extra.boe ? [{ q: 'Bøe number at the mean flowline state', v: round(extra.boe.piSS), u: '-' }, { q: 'Taitel criterion', v: extra.taitel.stable ? 'stable' : 'unstable', u: '' }, { q: 'Pressure amplitude / riser liquid head', v: round(extra.head), u: '-' }] : []),
      ], 'The exact boundary and initial conditions of the run with the numerical health indicators: cell size, mass-balance residuals and any mass removed by the positivity limiter.')],
      notes: [`${r.steps} time steps of ${fmt(r.dtMin * 1e3)}–${fmt(r.dtMax * 1e3)} ms on ${inp.nCells} cells of ${fmt(r.grid.ds)} m; gas RT/M = ${fmt(gasRT)} J/kg.`,
        slugging ? `Sustained oscillation: inlet pressure swings ${fmt((oscP.max - oscP.min) / BAR)} bar with a period of about ${fmt(oscP.period)} s.` : 'No sustained slugging in the second half of the run.',
        'Isothermal drift-flux model with a reduced liquid sound speed: use it for trends and surge sizes and confirm critical cases with a full transient flow-assurance simulator.'],
      advice: [
        adv('engineering', sev(gov.uc), `${inp.scenario === 'rampup' ? 'Ramp-up' : inp.scenario === 'shutdown' ? 'Shutdown and restart' : 'Constant-rate riser flow'}: peak pressure ${fmt(pMax / BAR)} bar, peak liquid rate ${fmt(qLmax)} kg/s (${fmt(qLmax / inp.mL)} × full rate), surge volume ${fmt(surgeVol)} m³.`,
          gov.uc > 1 ? (gov.name.startsWith('Liquid surge') ? `Lengthen the ramp (now ${fmt(inp.tRamp)} s), restart at reduced rate or provide ${fmt(surgeVol)} m³ of surge volume.` : gov.name.startsWith('Peak') ? 'Limit the restart rate or open the outlet before the inlet so gas is not compressed behind the liquid.' : 'Refine the mesh or raise the liquid sound speed: the solution is not conserving mass.')
            : slugging ? 'The flow cycles: test a smaller choke opening or riser-base gas lift in the riser-slugging scenario until the oscillation dies out.' : 'The transient stays within the pressure and surge limits. Check the mesh-converged value and the slowest credible ramp before fixing the operating procedure.',
          'ISO 13628-1 / API RP 17A flow-assurance design; drift-flux model after Evje & Fjelde (2002)'),
        adv('economic', 'info', `The transient delivers a ${fmt(surgeVol)} m³ liquid surge against ${fmt(inp.vSurge)} m³ of separator surge volume.`,
          surgeVol > inp.vSurge ? 'A slower ramp costs minutes of deferred production per restart; extra separator volume costs topsides weight for the whole field life — the procedure change is nearly always cheaper.' : 'The existing surge volume is enough: restart procedures can be set by production considerations rather than by slug handling.',
          ECON_BASIS, `Each avoided separator trip saves a restart and its flaring. Extra surge volume of ${fmt(Math.max(surgeVol - inp.vSurge, 0))} m³ would be about ${fmt(Math.max(surgeVol - inp.vSurge, 0) * 0.35)} t of vessel steel (≈ ${fmt(Math.max(surgeVol - inp.vSurge, 0) * 0.35 * CO2_STEEL)} t CO₂e, screening factor).`),
      ],
      publish: { 'pipeline.surgeVolume': surgeVol, 'pipeline.transientPmax': pMax },
    };
  },
  verify: [{
    // Constant inlet rates: after the start-up adjustment of the discrete initial state has died out (about 250 s here) the outlet liquid rate
    // must equal the inlet rate (20 kg/s), and the mass in the pipe may change only by what crossed the boundaries (round-off residual).
    name: 'Constant boundary conditions: outlet liquid rate tends to the inlet rate and mass is conserved',
    source: 'Conservation of mass (limiting case of constant boundary conditions)',
    inputs: { scenario: 'rampup', frac0: 1, tRamp: 10, tEnd: 300, nCells: 30 },
    expect: [{ kpi: 'qLend', value: 20, tol: 0.02 }, { kpi: 'massResidual', value: 0, abs: 1e-9 }],
  }],
};
