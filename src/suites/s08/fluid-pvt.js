// Compositional fluid behaviour for flow assurance: two-phase flash with a cubic equation of state, phase densities
// and viscosities, the phase envelope with the operating path of the line, and the liquid drop-out curve.
import { makeMixture, flashPT, properties, phaseEnvelope, upperSaturationPressure, vapourPressure, gasViscosityLGE, viscosityLBC, surfaceTension, waterContentBukacek } from '../../engine/flow/eos.js';
import { COMPONENT_IDS } from '../../engine/flow/components.js';
import { num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, linspace, ECON_BASIS, CO2_STEEL, BAR, KELVIN } from './_common.js';

const G_F = 'Fluid composition';
const G_O = 'Operating conditions';
const COMP = [
  { id: 'N2', mol: 0.6 }, { id: 'CO2', mol: 2.0 }, { id: 'C1', mol: 78.0 }, { id: 'C2', mol: 7.5 }, { id: 'C3', mol: 4.2 }, { id: 'nC4', mol: 2.4 },
  { id: 'nC5', mol: 1.5 }, { id: 'nC6', mol: 1.3 }, { id: 'nC7', mol: 1.2 }, { id: 'nC10', mol: 1.3 },
];
const R = 8.314462618;

export default {
  id: 's08-fluid-pvt',
  suite: 8,
  title: 'Fluid phase behaviour (equation of state)',
  summary: 'Flashes the produced fluid at line conditions, draws its phase envelope with the operating path, and shows how much liquid drops out as pressure falls.',
  covers: [14],
  fidelity: 'analytical',
  standards: ['Peng & Robinson (1976); Soave (1972) cubic equations of state', 'Michelsen (1980, 1982) stability analysis, flash and phase-envelope algorithms', 'Lohrenz, Bray & Clark (1964) viscosity; Lee, Gonzalez & Eakin (1966) gas viscosity', 'GPSA Engineering Data Book §20 (water content of gas)'],
  models: ['Peng–Robinson equation of state', 'Soave–Redlich–Kwong equation of state', 'Phase-equilibrium PT flash with stability test (Rachford–Rice)', 'Phase envelope by Newton continuation', 'Residual-property phase densities with volume translation', 'Parachor interfacial tension'],
  limits: 'Two fluid phases only (vapour and hydrocarbon liquid): free water, hydrates, wax and asphaltene solids are not part of the flash. Database components up to n-decane with default binary interaction coefficients; heavier fractions need characterisation and tuning to laboratory PVT data before the envelope is used for design. The critical point of the envelope is interpolated (about ±1 K, ±1 bar).',
  publishes: ['pipeline.vapourFraction', 'pipeline.cricondenbar'],
  searchTerms: 'phase envelope Peng-Robinson flash gas condensate liquid dropout flow assurance',
  inputs: [
    { key: 'comp', label: 'Composition', type: 'table', group: G_F, role: 'material', columns: [{ key: 'id', label: `Component (${COMPONENT_IDS.slice(0, 16).join(', ')}, …)`, type: 'text' }, { key: 'mol', label: 'Amount', unit: 'mol %', type: 'number' }],
      default: COMP, help: 'Hydrocarbon stream composition; amounts are normalised. Water is handled separately as the saturation water content of the gas.' },
    sel('eos', 'Equation of state', 'PR', [['PR', 'Peng–Robinson'], ['SRK', 'Soave–Redlich–Kwong']], G_F, 'control', 'Both use van der Waals mixing rules with default interaction coefficients and volume translation for density.'),
    num('pIn', 'Inlet pressure', 'bar', 150, 0.1, 1400, G_O, 'bc', 'Pressure at the pipeline inlet; the fluid is flashed here.', { bind: 'pipeline.pIn' }),
    num('Tin', 'Inlet temperature', '°C', 70, -100, 250, G_O, 'bc', 'Temperature at the inlet.'),
    num('pOut', 'Outlet pressure', 'bar', 60, 0.1, 1400, G_O, 'bc', 'Pressure at the outlet, the second end of the operating path.'),
    num('Tout', 'Outlet temperature', '°C', 10, -100, 250, G_O, 'bc', 'Temperature at the outlet.', { bind: 'pipeline.Tout' }),
    num('nPath', 'Points along the operating path', '-', 12, 2, 60, G_O, 'control', 'Flash points between inlet and outlet (linear in pressure and temperature).'),
  ],
  compute(inp) {
    const rows = (inp.comp || []).map((r) => ({ id: String(r.id || '').trim(), mol: Number(r.mol) })).filter((r) => r.id && r.mol > 0);
    need(rows.length > 0, 'The composition is empty: enter at least one component with a positive amount.');
    for (const r of rows) need(COMPONENT_IDS.includes(r.id), `Component '${r.id}' is not in the database. Available: ${COMPONENT_IDS.join(', ')}.`);
    need(new Set(rows.map((r) => r.id)).size === rows.length, 'A component is listed twice in the composition.');
    const tot = rows.reduce((s, r) => s + r.mol, 0), z = rows.map((r) => r.mol / tot), ids = rows.map((r) => r.id);
    const mix = run(() => makeMixture(ids, { eos: inp.eos, volumeShift: true }), 'Mixture model');
    const single = ids.length === 1;
    const state = (Pbar, TC) => {
      const P = Pbar * BAR, T = TC + KELVIN;
      let beta, x = z, y = z;
      if (single) {
        let ps = Infinity; try { ps = vapourPressure({ component: ids[0], T, eos: inp.eos }); } catch { ps = Infinity; }
        beta = P < ps || !Number.isFinite(ps) ? 1 : 0;
        if (!Number.isFinite(ps)) beta = properties({ mix, T, P, z }).phase === 'liquid' ? 0 : 1;
      } else {
        const f = run(() => flashPT({ mix, T, P, z }), `Flash at ${fmt(Pbar)} bar, ${fmt(TC)} °C`);
        beta = f.phases === 2 ? f.beta : (f.phase === 'vapour' ? 1 : 0);
        if (f.phases === 2) { x = f.x; y = f.y; }
      }
      const V = beta > 0 ? properties({ mix, T, P, z: y, phase: beta < 1 ? 'vapour' : 'stable' }) : null;
      const Lq = beta < 1 ? properties({ mix, T, P, z: x, phase: beta > 0 ? 'liquid' : 'stable' }) : null;
      const vV = V ? beta * V.v : 0, vL = Lq ? (1 - beta) * Lq.v : 0;
      const mV = V ? beta * V.M : 0, mL = Lq ? (1 - beta) * Lq.M : 0;
      return { P, T, beta, x, y, V, L: Lq, liqVol: vL / (vV + vL), gasMass: mV / (mV + mL), rhoMix: (mV + mL) / (vV + vL) };
    };
    const sIn = state(inp.pIn, inp.Tin), sOut = state(inp.pOut, inp.Tout);
    const nP = Math.round(inp.nPath);
    const path = linspace(0, 1, nP).map((f) => state(inp.pIn + (inp.pOut - inp.pIn) * f, inp.Tin + (inp.Tout - inp.Tin) * f));
    let env = null;
    if (!single) { try { env = phaseEnvelope({ mix, z, maxPoints: 250 }); } catch { env = null; } }
    let pSat = 0, satType = 'none';
    if (single) { try { pSat = vapourPressure({ component: ids[0], T: inp.Tin + KELVIN, eos: inp.eos }); satType = 'vapour pressure'; } catch { pSat = 0; } }
    else { try { const u = upperSaturationPressure({ mix, T: inp.Tin + KELVIN, z }); pSat = u.P; satType = u.type === 'bubble' ? 'bubble point' : 'dew point'; } catch { pSat = 0; } }
    // liquid drop-out with falling pressure at the outlet temperature
    const pTop = Math.max(env ? env.cricondenbar.P / BAR * 1.05 : 0, inp.pIn, inp.pOut) || inp.pIn;
    const pd = linspace(Math.max(1, 0.02 * pTop), pTop, 24);
    // single-phase states (above the envelope) carry no dropped-out liquid whatever the phase label of the dense fluid
    const drop = pd.map((q) => { try { const s = state(q, inp.Tout); return s.beta > 0 && s.beta < 1 ? s.liqVol * 100 : 0; } catch { return null; } });
    const o = sOut;
    const muG = o.V ? gasViscosityLGE({ T: o.T, rho: o.V.rho, M: o.V.M }) : 0;
    let muL = 0; if (o.L) { try { muL = viscosityLBC({ mix, T: o.T, z: o.x, rhoMolar: o.L.rhoMolar }); } catch { muL = 0; } }
    let ift = 0; if (o.L && o.V) { try { ift = surfaceTension({ mix, x: o.x, y: o.y, rhoMolarL: o.L.rhoMolar, rhoMolarV: o.V.rhoMolar }); } catch { ift = 0; } }
    const wc = waterContentBukacek({ T: Math.min(Math.max(o.T, 233), 500), P: Math.max(o.P, 1e5) });
    const twoPhaseOut = o.beta > 0 && o.beta < 1, twoPhaseIn = sIn.beta > 0 && sIn.beta < 1;
    const ccb = env ? env.cricondenbar.P / BAR : 0, cct = env ? env.cricondentherm.T - KELVIN : 0;
    const M = z.reduce((s, zi, i) => s + zi * mix.M[i], 0);
    const maxDrop = Math.max(0, ...drop.filter((v) => v !== null));
    return {
      kpis: [
        kpi('betaOut', 'Vapour mole fraction at the outlet', o.beta, '-', 'Fraction of the stream in the vapour phase at outlet conditions: 1 is all gas, 0 all liquid.'),
        kpi('liqVolOut', 'Liquid volume fraction at the outlet', o.liqVol * 100, '%', 'No-slip volume fraction of hydrocarbon liquid at outlet conditions; the liquid loading the line and the slug catcher must handle.'),
        kpi('betaIn', 'Vapour mole fraction at the inlet', sIn.beta, '-', 'Phase state at inlet conditions.'),
        kpi('rhoGas', 'Gas density at the outlet', o.V ? o.V.rho : 0, 'kg/m³', 'Density of the vapour phase at outlet conditions (zero when no vapour is present).'),
        kpi('rhoLiq', 'Liquid density at the outlet', o.L ? o.L.rho : 0, 'kg/m³', 'Density of the liquid phase with volume translation (zero when no liquid is present).'),
        kpi('rhoIn', 'Mixture density at the inlet', sIn.rhoMix, 'kg/m³', 'No-slip density of the whole stream at inlet conditions.'),
        kpi('Zgas', 'Gas compressibility factor at the outlet', o.V ? o.V.Z : 0, '-', 'Z of the vapour phase from the equation of state.'),
        kpi('muGas', 'Gas viscosity at the outlet', muG * 1e3, 'cP', 'Lee–Gonzalez–Eakin viscosity of the vapour phase.'),
        kpi('muLiq', 'Liquid viscosity at the outlet', muL * 1e3, 'cP', 'Lohrenz–Bray–Clark viscosity of the liquid phase (zero when no liquid).'),
        kpi('ift', 'Interfacial tension at the outlet', ift * 1e3, 'mN/m', 'Gas–liquid interfacial tension from parachors; it falls towards zero near the critical point.'),
        kpi('pSat', `Saturation pressure at the inlet temperature`, pSat / BAR, 'bar', `Upper phase-boundary pressure at the inlet temperature (${satType}); above it the fluid is single phase at that temperature.`),
        kpi('cricondenbar', 'Cricondenbar', ccb, 'bar', 'Highest pressure on the phase envelope: above it the fluid is single phase at every temperature (dense-phase transport).'),
        kpi('cricondentherm', 'Cricondentherm', cct, '°C', 'Highest temperature on the phase envelope: above it no liquid can form at any pressure.'),
        kpi('dropMax', 'Maximum liquid drop-out at the outlet temperature', maxDrop, 'vol %', 'Largest liquid volume fraction reached as pressure falls at the outlet temperature (retrograde condensation peak for a gas condensate).'),
        kpi('water', 'Saturation water content of the gas at the outlet', wc.kgPerSm3 * 1e6, 'mg/Sm³', 'Water the gas can carry at outlet conditions; water above this condenses and is available for hydrates and corrosion.'),
        kpi('M', 'Molar mass of the stream', M * 1e3, 'kg/kmol', 'Mole-weighted molar mass of the feed.'),
      ],
      checks: [],
      plots: [
        ...(env ? [plot('line', 'Phase envelope and operating path', 'Temperature', '°C', 'Pressure', 'bar', [
          { name: 'Bubble line', x: env.bubble.map((q) => round(q.T - KELVIN)), y: env.bubble.map((q) => round(q.P / BAR)) }, { name: 'Dew line', x: env.dew.map((q) => round(q.T - KELVIN)), y: env.dew.map((q) => round(q.P / BAR)) },
          ...(env.critical ? [{ name: 'Critical point', x: [round(env.critical.T - KELVIN)], y: [round(env.critical.P / BAR)] }] : []),
          { name: 'Operating path (inlet to outlet)', x: path.map((q) => round(q.T - KELVIN)), y: path.map((q) => round(q.P / BAR)) },
        ], 'Inside the envelope gas and liquid coexist. Where the operating path enters it, liquid starts to form in the line; a path that stays above the cricondenbar is dense-phase throughout.')] : []),
        plot('line', 'Phase split along the operating path', 'Pressure', 'bar', 'Fraction', '-', [
          { name: 'Vapour mole fraction', x: path.map((q) => round(q.P / BAR)), y: path.map((q) => round(q.beta)) }, { name: 'Liquid volume fraction', x: path.map((q) => round(q.P / BAR)), y: path.map((q) => round(q.liqVol)) },
        ], 'Vapour fraction and no-slip liquid volume fraction from inlet to outlet. A liquid fraction that rises as pressure falls is retrograde condensation.'),
        plot('line', 'Liquid drop-out at the outlet temperature', 'Pressure', 'bar', 'Liquid volume', '%', [
          { name: `At ${fmt(inp.Tout)} °C`, x: pd.filter((_, i) => drop[i] !== null).map((q) => round(q)), y: drop.filter((v) => v !== null).map((v) => round(v)) },
        ], 'Constant-composition liquid drop-out as pressure is reduced. The peak is the worst liquid loading; operating to either side of it reduces liquid handling.'),
      ],
      tables: [tbl('Phase compositions at the outlet', [col('id', 'Component'), col('z', 'Feed', 'mol %'), col('x', 'Liquid', 'mol %'), col('y', 'Vapour', 'mol %'), col('K', 'K-value', '-')],
        ids.map((id, i) => ({ id, z: round(z[i] * 100), x: round(o.x[i] * 100), y: round(o.y[i] * 100), K: round(o.x[i] > 0 ? o.y[i] / o.x[i] : 0) })),
        'Feed, liquid and vapour compositions at outlet conditions with the equilibrium ratio of each component: light components concentrate in the vapour (K above 1), heavy ones in the liquid.')],
      notes: [`${inp.eos === 'PR' ? 'Peng–Robinson' : 'Soave–Redlich–Kwong'} with default interaction coefficients and volume translation; ${ids.length} components.`,
        `Inlet: ${twoPhaseIn ? 'two-phase' : sIn.beta >= 1 ? 'single-phase vapour or dense phase' : 'single-phase liquid'}; outlet: ${twoPhaseOut ? `two-phase with ${fmt(o.liqVol * 100)} % liquid by volume` : o.beta >= 1 ? 'single-phase vapour or dense phase' : 'single-phase liquid'}.`,
        env ? (env.complete ? `Critical point about ${fmt(env.critical ? env.critical.T - KELVIN : 0)} °C, ${fmt(env.critical ? env.critical.P / BAR : 0)} bar (interpolated).` : 'The phase envelope could not be closed at the critical point; the traced branches are shown.') : 'Single component or envelope not available: no phase envelope is drawn.'],
      advice: [
        adv('engineering', twoPhaseOut ? 'warn' : 'ok', twoPhaseOut ? `The stream is two-phase at the outlet: ${fmt(o.liqVol * 100)} % liquid by volume at ${fmt(inp.pOut)} bar and ${fmt(inp.Tout)} °C; cricondenbar ${fmt(ccb)} bar.` : `The stream is single phase at the outlet; the upper saturation pressure at the inlet temperature is ${fmt(pSat / BAR)} bar and the cricondenbar ${fmt(ccb)} bar.`,
          twoPhaseOut ? 'Run the multiphase and slugging modules with this liquid loading, and size liquid handling for the drop-out peak; tune the equation of state to laboratory PVT data before design.' : `Keep the arrival pressure above ${fmt(ccb)} bar if dense-phase transport is intended; otherwise check turndown, when lower pressure and temperature may bring the path into the envelope.`,
          'Michelsen & Mollerup, Thermodynamic Models (2007); flow-assurance practice per ISO 13628-1 / API RP 17A'),
        adv('economic', 'info', `Maximum liquid drop-out at the outlet temperature is ${fmt(maxDrop)} % by volume; the operating pressure sets how much of it the line carries.`,
          twoPhaseOut ? `Operating above the cricondenbar (${fmt(ccb)} bar) removes liquid handling and slug-catcher volume but needs a thicker wall and more compression; compare the two concepts on steel tonnes and compressor energy.` : 'Single-phase transport avoids slug catchers and pigging for liquid management; the cost is carried in wall thickness and compression energy, which should be checked against a two-phase concept at lower pressure.',
          ECON_BASIS, `Dense-phase operation trades steel (about ${CO2_STEEL} t CO₂e per tonne, screening factor) and compressor fuel against avoided liquid handling and flaring during slug-induced trips; the phase envelope is the basis for that trade.`),
      ],
      publish: { 'pipeline.vapourFraction': o.beta, 'pipeline.cricondenbar': ccb * BAR },
    };
  },
  verify: [{
    // Ideal-gas limit: pure methane at 0.1 bar and 300 K has Z = 1 to four figures, so ρ = P M/(R T) = 1e4·0.016043/(8.31446·300) = 0.06432 kg/m³.
    name: 'Ideal-gas limit of the equation of state (methane at 0.1 bar)',
    source: 'Ideal-gas law (limiting case of any cubic equation of state)',
    inputs: { comp: [{ id: 'C1', mol: 100 }], pIn: 0.1, Tin: 26.85, pOut: 0.1, Tout: 26.85 },
    expect: [{ kpi: 'rhoGas', value: 0.06432, tol: 0.002 }, { kpi: 'Zgas', value: 1, tol: 0.001 }],
  }, {
    // Vapour pressure of propane at 300 K: 9.98 bar (NIST Chemistry WebBook, saturation properties of propane).
    name: 'Vapour pressure of propane at 300 K',
    source: 'NIST Chemistry WebBook, propane saturation pressure at 300 K = 0.998 MPa',
    inputs: { comp: [{ id: 'C3', mol: 100 }], pIn: 5, Tin: 26.85, pOut: 5, Tout: 26.85, eos: 'PR' },
    expect: [{ kpi: 'pSat', value: 9.98, tol: 0.02 }],
  }],
};
