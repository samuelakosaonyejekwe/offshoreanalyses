// Thermal stress and end expansion: temperature and pressure profiles, fully restrained and friction-released
// effective force along the line, virtual anchors, end expansion against tie-in spool flexibility, and spool stress.
import { restrainedAxialForce, temperatureProfile, pressureProfile, endExpansion, spoolGuidedCantilever } from '../../engine/pipeline/expansion.js';
import { equivalentStressCheck } from '../../engine/pipeline/strength.js';
import { pipeInputs, pipeFrom, lengthInput, depthInput, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, sev, thin, worst, interp,
  steelCo2, CONFIRM, LIMIT_NO_FE, ECON_BASIS, G, BAR, MPA } from './_common.js';

const G_O = 'Operating conditions';
const G_E = 'End conditions';

export default {
  id: 's08-thermal-expansion',
  suite: 8,
  title: 'Thermal stress and end expansion',
  summary: 'Computes the temperature-driven axial force along the line, where the pipe is anchored by friction, and how far each end moves into its tie-in spool.',
  covers: [13],
  fidelity: 'analytical',
  standards: ['DNV-RP-F110 (2018) §6 (effective axial force and friction build-up)', 'DNV-ST-F101 (2021) §4.7.3 (effective axial force)', 'Palmer & King, Subsea Pipeline Engineering (2008) ch. 12', 'ASME B31.8 (2022) §833 (expansion stress)'],
  models: ['Thermal expansion and Poisson / end-cap pressure strain', 'Steady exponential temperature decay (Newton cooling)', 'Axial force equilibrium with Coulomb friction (virtual anchors)', 'Linear spool spring end restraint', 'Guided-cantilever spool stress', 'von Mises equivalent stress'],
  limits: `First heat-up from the as-laid state on a straight line with uniform or prescribed friction; no lateral buckling (which releases force — see the global buckling module), no route curvature and no cyclic effects (see the walking module). The spool is an elastic guided cantilever without elbow flexibility or soil resistance. ${LIMIT_NO_FE}`,
  publishes: ['pipeline.endExpansion', 'pipeline.endForce', 'pipeline.effectiveForce'],
  searchTerms: 'pipeline end expansion virtual anchor thermal stress tie-in spool',
  inputs: [
    ...pipeInputs(),
    lengthInput(10, 'Pipe'),
    depthInput(120, G_O),
    num('Tin', 'Inlet temperature', '°C', 90, -20, 200, G_O, 'bc', 'Fluid (wall) temperature at the hot end.', { design: true, uncertain: { dist: 'normal', cov: 0.05 } }),
    num('Tamb', 'Ambient sea-water temperature', '°C', 5, -2, 35, G_O, 'bc', 'Temperature the pipe decays towards along the route.', { bind: 'site.Tbed' }),
    num('Tinst', 'Installation temperature', '°C', 5, -10, 40, G_O, 'ic', 'Pipe temperature when laid; reference for thermal strain.'),
    num('U', 'Overall heat-transfer coefficient', 'W/m²K', 4, 0, 200, G_O, 'operation', 'Referred to the steel outer diameter. Zero gives a uniform temperature along the line.'),
    num('mdot', 'Mass flow rate', 'kg/s', 60, 0.1, 5000, G_O, 'bc', 'Total mass rate, with the heat capacity it sets the temperature decay length.'),
    num('cp', 'Fluid specific heat', 'J/kgK', 2300, 500, 5000, G_O, 'operation', 'Mass-averaged heat capacity of the stream.'),
    { key: 'Tprofile', label: 'Measured or simulated temperature profile (optional)', type: 'series', default: { x: [], y: [] }, xLabel: 'Chainage', xUnit: 'km', yLabel: 'Temperature', yUnit: '°C', group: G_O, role: 'bc',
      help: 'When at least two points are given they replace the exponential profile (e.g. the result of the multiphase flow module).' },
    num('pIn', 'Inlet pressure', 'bar', 200, 0, 1500, G_O, 'bc', 'Internal pressure at the hot end; pressure at installation is taken as zero.', { uncertain: { dist: 'normal', cov: 0.03 } }),
    num('dpdx', 'Frictional pressure gradient', 'bar/km', 1.5, 0, 50, G_O, 'operation', 'Pressure loss along the line.'),
    num('Hlay', 'Residual lay tension', 'kN', 100, 0, 5000, G_O, 'ic', 'Effective tension left after laying.'),
    num('muA', 'Axial friction coefficient', '-', 0.5, 0.02, 1.5, 'Pipe–soil friction', 'soil', 'Equivalent axial friction (low estimate for end expansion, high for anchor force).', { bind: 'pipeline.axialFriction', uncertain: { dist: 'lognormal', cov: 0.25 } }),
    sel('hotEnd', 'Hot-end restraint', 'spool', [['free', 'Free'], ['spool', 'Tie-in spool (guided cantilever)'], ['spring', 'Spring of stated stiffness']], G_E, 'bc', 'Restraint offered by the tie-in at the inlet end.'),
    sel('coldEnd', 'Cold-end restraint', 'free', [['free', 'Free'], ['spool', 'Tie-in spool (guided cantilever)'], ['spring', 'Spring of stated stiffness']], G_E, 'bc', 'Restraint offered by the tie-in at the outlet end.'),
    num('Lleg', 'Spool leg length', 'm', 40, 3, 200, G_E, 'geometry', 'Length of the spool leg perpendicular to the pipeline that absorbs the expansion.', { design: true }),
    num('kEnd', 'End spring stiffness', 'kN/m', 500, 0, 1e6, G_E, 'bc', 'Used for ends set to "spring" (PLET sliding foundation, riser base).'),
    num('sigSpool', 'Allowable spool bending stress range', 'MPa', 300, 50, 600, 'Acceptance criteria', 'criteria', `Allowable for the displacement-induced bending stress in the spool (e.g. 0.72–0.9 SMYS less pressure stresses). ${CONFIRM}`),
    num('eta', 'Usage factor η on yield (pipeline)', '-', 0.9, 0.5, 1, 'Acceptance criteria', 'criteria', `For the equivalent stress along the restrained pipeline. ${CONFIRM}`),
    { key: 'nInt', label: 'Integration intervals along the line', type: 'int', default: 1000, min: 100, max: 8000, group: 'Numerics', role: 'control', help: 'Number of intervals for the friction build-up and strain integration.' },
  ],
  mesh: { key: 'nInt', label: 'Integration intervals', levels: [250, 500, 1000, 2000], qoi: 'expHot' },
  compute(inp, ctx) {
    const p = pipeFrom(inp, ctx);
    need(p.ws > 0, 'The pipe is buoyant in operation: there is no seabed friction to anchor it.');
    const L = inp.Lkm * 1e3, pe = p.rhoW * G * inp.depth;
    const series = inp.Tprofile && inp.Tprofile.x && inp.Tprofile.x.length >= 2 ? inp.Tprofile : null;
    if (series) for (let i = 1; i < series.x.length; i++) need(series.x[i] > series.x[i - 1], 'The temperature profile chainage must increase.');
    const lambda = inp.U > 0 ? inp.mdot * inp.cp / (inp.U * Math.PI * p.D) : Infinity;
    const Tat = (x) => (series ? interp(series.x.map((v) => v * 1e3), series.y, x) : (Number.isFinite(lambda) ? temperatureProfile({ x, TinC: inp.Tin, TambC: inp.Tamb, lambda }) : inp.Tin));
    const pAt = (x) => Math.max(inp.pIn * BAR - inp.dpdx * BAR / 1e3 * x, 0);
    const f = inp.muA * p.ws;
    const kSpool = spoolGuidedCantilever({ delta: 1, D: p.D, I: p.I, Lleg: inp.Lleg, E: p.E }).stiffness;
    const kOf = (end) => (end === 'spool' ? kSpool : end === 'spring' ? inp.kEnd * 1e3 : undefined);
    const args = { L, As: p.As, Ai: p.Ai, friction: f, deltaT: (x) => Tat(x) - inp.Tinst, deltaPi: pAt, H: inp.Hlay * 1e3, E: p.E, alpha: p.alpha, nu: p.nu, n: inp.nInt };
    const kH = kOf(inp.hotEnd), kC = kOf(inp.coldEnd);
    if (kH > 0) args.kHot = kH;
    if (kC > 0) args.kCold = kC;
    const ex = run(() => endExpansion(args), 'End expansion');
    const free = run(() => endExpansion({ ...args, kHot: undefined, kCold: undefined }), 'Free-end expansion');
    const r0 = restrainedAxialForce({ H: inp.Hlay * 1e3, deltaPi: pAt(0), Ai: p.Ai, As: p.As, deltaT: Tat(0) - inp.Tinst, E: p.E, alpha: p.alpha, nu: p.nu });
    const sigTh = p.E * p.alpha * (Tat(0) - inp.Tinst);
    // equivalent stress along the line with the released effective force
    const [xs, S0s, Ss] = thin([ex.x, ex.S0, ex.S], 200);
    const vm = xs.map((x, i) => { const pi = pAt(x); return equivalentStressCheck({ D: p.D, t: p.t, fy: p.smys, eta: inp.eta, pi, pe, N: Ss[i] + pi * p.Ai - pe * p.Ae }).sigmaE; });
    const iMax = vm.indexOf(Math.max(...vm));
    const spools = [];
    for (const [end, d] of [['hot', ex.expansionHot], ['cold', ex.expansionCold]]) {
      if (inp[`${end}End`] !== 'spool') continue;
      const s = spoolGuidedCantilever({ delta: Math.abs(d), D: p.D, I: p.I, Lleg: inp.Lleg, sigmaAllow: inp.sigSpool * MPA, E: p.E });
      spools.push({ end, d, ...s });
    }
    const checks = [
      chk('Equivalent stress in the pipeline', vm[iMax], inp.eta * p.smys, 'Pa', 'DNV-ST-F101 (2021) Sec. 5 ASD format; ASME B31.8 §833.4',
        'Guards against yielding of the restrained hot pipe under hoop tension combined with axial compression. Exceeding it means the wall yields where the pipe is anchored by friction.'),
      ...spools.map((s) => chk(`Spool bending stress, ${s.end} end`, s.stress, inp.sigSpool * MPA, 'Pa', 'Guided cantilever method (Kellogg); ASME B31.8 §833.8',
        'Guards against overstressing the tie-in spool as it bends to absorb the pipeline end movement; exceeding it risks yielding or fatigue of the spool and overloading the connector.')),
    ];
    const gov = worst(checks);
    const legReq = Math.max(0, ...spools.map((s) => s.requiredLeg));
    const expMax = Math.max(Math.abs(ex.expansionHot), Math.abs(ex.expansionCold));
    const dLeg = legReq - inp.Lleg;
    const spoolSteel = (len) => 2.2 * len * p.steelKgM / 1e3; // Z- or L-spool: about 2.2 × the leg length of pipe
    return {
      kpis: [
        kpi('expHot', 'Hot-end expansion', ex.expansionHot, 'm', 'Outward movement of the inlet end against its tie-in restraint. The spool, jumper or sliding PLET must accommodate it.'),
        kpi('expCold', 'Cold-end expansion', ex.expansionCold, 'm', 'Outward movement of the outlet end.'),
        kpi('expFree', 'Hot-end expansion with free ends', free.expansionHot, 'm', 'Upper bound of the hot-end movement if the tie-in offered no resistance.'),
        kpi('anchorHot', 'Virtual anchor from the hot end', ex.anchorHot, 'm', 'Distance from the inlet at which friction has built up to the fully restrained force; beyond it the pipe does not move.'),
        kpi('anchorCold', 'Virtual anchor from the hot end (cold side)', ex.anchorCold, 'm', 'Position of the cold-side anchor. Equal to the hot one for a short pipeline that never reaches full restraint.'),
        kpi('Smax', 'Maximum effective compression', ex.maxCompression / 1e3, 'kN', 'Largest compressive effective force in the line after friction release; the driver for global buckling.'),
        kpi('S0in', 'Fully restrained force at the inlet', r0.S0 / 1e3, 'kN', 'Effective force the inlet section would carry if it could not move (negative = compression).'),
        kpi('Rhot', 'Hot-end reaction on the tie-in', ex.Rhot / 1e3, 'kN', 'Axial force the pipeline applies to the spool or structure at the inlet end.'),
        kpi('sigTh', 'Restrained thermal stress at the inlet', sigTh / MPA, 'MPa', 'E·α·ΔT: the axial compressive stress produced by fully restraining the inlet temperature rise.'),
        kpi('vmMax', 'Maximum equivalent stress', vm[iMax] / MPA, 'MPa', 'Largest von Mises stress along the pipeline with the friction-released axial force.', inp.eta * p.smys / MPA, 'max'),
        kpi('spoolStress', 'Spool bending stress', Math.max(0, ...spools.map((s) => s.stress)) / MPA, 'MPa', 'Bending stress range in the spool leg from the end movement (zero when no end has a spool).', inp.sigSpool, 'max'),
        kpi('legReq', 'Spool leg length required', legReq, 'm', 'Leg length at which the spool bending stress equals the allowable for the computed expansion.'),
        kpi('lambda', 'Temperature decay length', Number.isFinite(lambda) ? lambda / 1e3 : 9999, 'km', 'Distance over which the temperature excess falls to 37 %; 9999 denotes a uniform temperature.'),
        kpi('Tout', 'Temperature at the outlet', Tat(L), '°C', 'Wall temperature at the cold end from the profile used.'),
      ],
      checks,
      plots: [
        plot('line', 'Effective axial force along the pipeline', 'Chainage', 'km', 'Effective force (tension +)', 'kN', [
          { name: 'Fully restrained', x: xs.map((x) => round(x / 1e3)), y: S0s.map((v) => round(v / 1e3)) }, { name: 'With friction release', x: xs.map((x) => round(x / 1e3)), y: Ss.map((v) => round(v / 1e3)) },
        ], 'The restrained force follows the temperature and pressure profile; near each end the actual force ramps up at the friction gradient. Where the two curves meet are the virtual anchors.'),
        plot('line', 'Temperature and pressure profile', 'Chainage', 'km', 'Temperature', '°C', [
          { name: 'Temperature', x: xs.map((x) => round(x / 1e3)), y: xs.map((x) => round(Tat(x))) }, { name: 'Pressure (bar ÷ 10)', x: xs.map((x) => round(x / 1e3)), y: xs.map((x) => round(pAt(x) / BAR / 10)) },
        ], 'The loading that drives expansion: temperature decays towards ambient and pressure falls with friction along the line.'),
        plot('line', 'Equivalent stress along the pipeline', 'Chainage', 'km', 'von Mises stress', 'MPa', [
          { name: 'von Mises', x: xs.map((x) => round(x / 1e3)), y: vm.map((v) => round(v / MPA)) }, { name: 'Allowable', x: [0, round(L / 1e3)], y: [round(inp.eta * p.smys / MPA), round(inp.eta * p.smys / MPA)], dash: true },
        ], 'Combined hoop and axial stress. It peaks at the hot-side virtual anchor where compression is greatest.'),
      ],
      tables: [tbl('End conditions', [col('end', 'End'), col('bc', 'Restraint'), col('k', 'Stiffness', 'kN/m'), col('d', 'Expansion', 'm'), col('R', 'Reaction', 'kN')], [
        { end: 'Hot (inlet)', bc: inp.hotEnd, k: round((kH || 0) / 1e3), d: round(ex.expansionHot), R: round(ex.Rhot / 1e3) },
        { end: 'Cold (outlet)', bc: inp.coldEnd, k: round((kC || 0) / 1e3), d: round(ex.expansionCold), R: round(ex.Rcold / 1e3) },
      ], 'Boundary condition used at each end with the resulting movement and the axial load passed to the tie-in.')],
      notes: [`${ex.regime === 'long' ? 'Long' : 'Short'} pipeline: ${ex.regime === 'long' ? `fully restrained between ${fmt(ex.anchorHot / 1e3)} and ${fmt(ex.anchorCold / 1e3)} km` : `it slides over its whole length about a single anchor at ${fmt(ex.anchorHot / 1e3)} km`}.`,
        `Axial friction ${fmt(f)} N/m; spool stiffness ${fmt(kSpool / 1e3)} kN/m for a ${fmt(inp.Lleg)} m leg.`,
        ex.gradientExceedsFriction ? 'The restrained force changes faster than friction can follow in part of the line: the anchored length is partly slipping.' : 'End expansion is governed by the low-estimate friction; anchor force and buckling by the high estimate.'],
      advice: [
        adv('engineering', sev(gov.uc), `The hot end moves ${fmt(ex.expansionHot)} m and the cold end ${fmt(ex.expansionCold)} m; governing utilisation ${fmt(gov.uc)} (${gov.name}).`,
          spools.some((s) => s.stress > inp.sigSpool * MPA) ? `Lengthen the spool leg to at least ${fmt(legReq)} m, add a second leg (Z-spool), or reduce the expansion with a sliding PLET or mid-line anchor.` : `Design the tie-ins for ${fmt(expMax)} m of movement with the low-estimate friction; check global buckling for the ${fmt(ex.maxCompression / 1e3)} kN compression.`,
          'DNV-RP-F110 (2018) §6; DNV-ST-F101 (2021) §4.7.3; guided cantilever method for the spool'),
        adv('economic', 'info', `A spool leg of ${fmt(legReq)} m is needed against ${fmt(inp.Lleg)} m specified; each metre of leg is about ${fmt(2.2 * p.steelKgM)} kg of spool pipe.`,
          dLeg > 0 ? `The longer spool adds about ${fmt(spoolSteel(dLeg))} t of fabricated pipe and may exceed the lift and deck limits of the installation vessel; a sliding PLET is the usual alternative.` : `The spool could be shortened by about ${fmt(-dLeg)} m (${fmt(spoolSteel(-dLeg))} t) if metrology and installation tolerances allow.`,
          ECON_BASIS, `${steelCo2(spoolSteel(Math.abs(dLeg)))} Insulation that lowers the U-value raises the far-end temperature and the expansion; flow assurance and expansion design need to be traded together.`),
      ],
      publish: { 'pipeline.endExpansion': ex.expansionHot, 'pipeline.endForce': ex.Rhot, 'pipeline.effectiveForce': -ex.maxCompression },
    };
  },
  verify: [{
    // Uniform temperature rise ΔT = 50 K, no pressure, free ends, long line. P = A_s E α ΔT = 0.013909·207e9·1.17e-5·50 = 1.6843 MN;
    // friction f = μ w = 0.5·1000 = 500 N/m. Anchor length L_a = P/f = 3368.7 m; expansion δ = P²/(2 f E A_s)
    // = (1.6843e6)²/(2·500·2.8792e9) = 0.9853 m.
    name: 'Long pipeline with uniform temperature rise: anchor length P/f and end expansion P²/(2fEA)',
    source: 'Palmer & King, Subsea Pipeline Engineering (2008) ch. 12; hand calculation',
    inputs: { D: 323.9, t: 14.3, wsGiven: 1000, Lkm: 20, Tin: 55, Tamb: 5, Tinst: 5, U: 0, pIn: 0, dpdx: 0, Hlay: 0, muA: 0.5, hotEnd: 'free', coldEnd: 'free', depth: 0 },
    expect: [{ kpi: 'expHot', value: 0.9853, tol: 0.005 }, { kpi: 'anchorHot', value: 3368.7, tol: 0.005 }],
  }],
};
