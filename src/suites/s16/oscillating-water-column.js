// Oscillating-water-column wave-energy converter: the water column as a rigid piston coupled to a compressible air
// chamber and an air turbine (Wells: linear; impulse: quadratic), in regular and irregular waves — natural period,
// optimum turbine damping, capture width against the axisymmetric limit, chamber pressure and column stroke, the
// power matrix with the site scatter diagram and the response in the extreme sea state.
// Engineering units (kW, kPa, MWh) are converted to SI at the top of compute.
import { columnResonance, chamberCompliance, columnBody, owcRegular, optimalTurbineDamping, owcIrregular, owcPowerMatrix } from '../../engine/aero/owc.js';
import { waveSpectrum, regularWavePower, wavePower } from '../../engine/aero/wec.js';
import { site, adv, check, worst, sev, fmt, round, run, need, linspace, G, LIMIT_CODE, LIMIT_CONTROL, CONFIRM, CO2_CONCRETE, energyCo2Text, costInput, money } from './_common.js';

const TPS = [6, 8, 10, 12, 14];
// Indicative exposed North-Atlantic scatter diagram (occurrence %, rows Hs, columns Tp): replace with the site hindcast.
const SCATTER = [[0.5, 4, 6, 3, 1, 0], [1.5, 6, 14, 10, 4, 1], [2.5, 1, 9, 11, 5, 1], [3.5, 0, 3, 7, 4, 1], [4.5, 0, 1, 3, 2, 1], [5.5, 0, 0, 1, 1, 0]];
// Mean zero-crossing period of the response taken as 0.75 Tp (between 0.71 for a Bretschneider and 0.78 for a JONSWAP γ = 3.3 sea).
const TZ_OVER_TP = 0.75;

export default {
  id: 's16-oscillating-water-column', suite: 16, title: 'Oscillating water column (air-turbine power take-off)',
  summary: 'Water-column resonance, optimum air-turbine damping, capture width, chamber pressure and stroke of an oscillating-water-column converter, with its power matrix, annual energy and extreme-sea response.',
  covers: [2, 12], fidelity: 'analytical',
  standards: ['IEC TS 62600-100:2012 (power performance of wave energy converters)', 'IEC TS 62600-2:2019 (design requirements)', 'DNV-RP-C205 (2021) §3.5.5 (wave spectra)', 'Evans (1978) J. Inst. Maths Applics 22:423–433', 'Sarmento & Falcão (1985) J. Fluid Mech. 150:467–485', 'Falcão & Henriques (2016) Renewable Energy 85:1391–1424'],
  models: ['Rigid-piston water column (single degree of freedom, frequency domain)', 'Froude–Krylov excitation and Haskind radiation damping on the column mouth', 'Linearised isentropic air-chamber compressibility', 'Wells turbine (linear pressure–flow)', 'Impulse turbine (quadratic pressure–flow, harmonic and statistical linearisation)', 'Optimum turbine damping (golden-section search)', 'JONSWAP / Bretschneider spectra', 'Power matrix × scatter diagram', 'Rayleigh most-probable-maximum response'],
  limits: `The water column is a rigid piston with closed-form Froude–Krylov excitation and Haskind radiation damping of a circular mouth: there is no boundary-element (diffraction–radiation) solution, no sloshing mode and no structure motion (a fixed or very large host structure is assumed). The air chamber is linearised about atmospheric pressure (valid while the pressure amplitude stays below about 20 % of atmospheric) and the turbine has a constant mean efficiency: stall of a Wells turbine at high flow, relief-valve dynamics and rotational-speed control are not modelled. Linear theory overstates the response in steep and extreme seas. ${LIMIT_CONTROL} ${LIMIT_CODE}`,
  publishes: ['owc.ratedPower', 'owc.aep', 'owc.capacityFactor', 'owc.turbineDamping'],
  searchTerms: 'oscillating water column OWC Wells turbine impulse turbine air chamber capture width',
  inputs: [
    { key: 'radius', label: 'Water-column radius', unit: 'm', type: 'number', default: 5, min: 0.5, max: 15, group: 'Chamber geometry', role: 'geometry', design: true, help: 'Circular column of equal plan area for other shapes.' },
    { key: 'draft', label: 'Column draft (lip submergence)', unit: 'm', type: 'number', default: 10, min: 1, max: 60, group: 'Chamber geometry', role: 'geometry', design: true, help: 'Depth of the column mouth below still water. With the added length it sets the natural period: T = 2π√((d + lₐ)/g).' },
    { key: 'laFactor', label: 'Added length / radius', unit: '-', type: 'number', default: 0.6, min: 0, max: 1.2, group: 'Chamber geometry', role: 'bc', uncertain: { dist: 'normal', cov: 0.15 }, help: 'Entrained water below the mouth as a length: between 4/(3π) = 0.42 (unbaffled disc, half added mass) and 8/(3π) = 0.85 (piston in an infinite baffle).' },
    { key: 'airHeight', label: 'Still-water air height in the chamber', unit: 'm', type: 'number', default: 12, min: 0.5, max: 40, group: 'Chamber geometry', role: 'ic', design: true, help: 'Initial condition: chamber air volume = column area × this height. It is also the freeboard between the water surface and the turbine duct.' },
    { key: 'compress', label: 'Include air compressibility', type: 'bool', default: true, group: 'Chamber geometry', role: 'bc', help: 'Off = incompressible air (model scale, or the limiting case of a very small chamber).' },
    { key: 'fExc', label: 'Factor on the Froude–Krylov excitation', unit: '-', type: 'number', default: 1, min: 0.3, max: 2, group: 'Hydrodynamics', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Diffraction by the host structure (a breakwater front wall roughly doubles the incident amplitude; an open spar leaves it near 1).' },
    { key: 'Bvisc', label: 'Linearised entry-loss damping', unit: 'kN·s/m', type: 'number', default: 20, min: 0, max: 5000, group: 'Hydrodynamics', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.4 }, help: 'Vortex shedding at the lip, as an equivalent linear damper on the column; calibrate against model tests.' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 50, min: 5, max: 500, group: 'Site', role: 'environment', bind: 'site.depth' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 995, max: 1035, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'patm', label: 'Atmospheric pressure', unit: 'kPa', type: 'number', default: 101.325, min: 85, max: 108, group: 'Site', role: 'environment' },
    { key: 'Hs', label: 'Design operating significant wave height', unit: 'm', type: 'number', default: 2.5, min: 0.2, max: 8, group: 'Operating sea state', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Incident-wave boundary condition for the headline response.' },
    { key: 'Tp', label: 'Design operating peak period', unit: 's', type: 'number', default: 9, min: 3, max: 18, group: 'Operating sea state', role: 'bc' },
    { key: 'gamma', label: 'JONSWAP peak-enhancement factor', unit: '-', type: 'number', default: 2, min: 1, max: 7, group: 'Operating sea state', role: 'bc', help: '1 = Bretschneider, 3.3 = fetch-limited wind sea.' },
    { key: 'scatter', label: 'Scatter diagram (occurrence, any scale)', type: 'table', group: 'Wave resource', role: 'bc',
      columns: [{ key: 'hs', label: 'Hs', unit: 'm', type: 'number' }, ...TPS.map((t) => ({ key: `t${t}`, label: `Tp ${t} s`, unit: '%', type: 'number' }))],
      default: SCATTER.map(([hs, ...p]) => ({ hs, ...Object.fromEntries(TPS.map((t, i) => [`t${t}`, p[i]])) })), help: 'Joint occurrence of Hs and Tp (rows: bin-centre Hs). Normalised internally.' },
    site('Hs100', { group: 'Survival sea state' }), site('Tp100', { group: 'Survival sea state' }),
    { key: 'turbine', label: 'Air turbine', type: 'select', default: 'wells', group: 'Air turbine', role: 'geometry', options: [{ value: 'wells', label: 'Wells (flow proportional to pressure)' }, { value: 'impulse', label: 'Impulse (pressure proportional to flow squared)' }] },
    { key: 'dampMode', label: 'Wells turbine damping', type: 'select', default: 'optimal', group: 'Air turbine', role: 'control', showIf: { key: 'turbine', in: ['wells'] }, options: [{ value: 'optimal', label: 'Optimum for the design peak period' }, { value: 'fixed', label: 'Entered value' }],
      help: 'The damping p/q of a Wells turbine follows from its diameter and rotational speed; the optimum is found for a regular wave at the design peak period and then held for every sea state.' },
    { key: 'damping', label: 'Turbine damping p/q', unit: 'Pa·s/m³', type: 'number', default: 60, min: 0.5, max: 1e5, group: 'Air turbine', role: 'bc', showIf: { key: 'dampMode', in: ['fixed'] }, design: true },
    { key: 'ki', label: 'Impulse-turbine coefficient k in p = k·q·|q|', unit: 'Pa·s²/m⁶', type: 'number', default: 0.5, min: 1e-4, max: 1e4, group: 'Air turbine', role: 'bc', showIf: { key: 'turbine', in: ['impulse'] }, design: true },
    { key: 'etaTurb', label: 'Mean turbine–generator efficiency', unit: '-', type: 'number', default: 0.5, min: 0.1, max: 0.9, group: 'Air turbine', role: 'operation', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Pneumatic to electrical, averaged over the cycle: about 0.4–0.6 for Wells turbines, 0.35–0.5 for impulse turbines.' },
    { key: 'ratedKW', label: 'Rated electrical power per chamber', unit: 'kW', type: 'number', default: 100, min: 5, max: 5000, group: 'Air turbine', role: 'geometry', design: true },
    { key: 'pRelief', label: 'Relief-valve (design) chamber pressure', unit: 'kPa', type: 'number', default: 35, min: 1, max: 60, group: 'Acceptance', role: 'criteria', design: true, help: 'Gauge pressure amplitude the chamber and turbine are designed for; also keeps the linearised air model within its range.' },
    { key: 'survValve', label: 'Survival-mode conductance / operating conductance', unit: '-', type: 'number', default: 0.8, min: 0.01, max: 100, group: 'Air turbine', role: 'control', help: 'Valve setting in the extreme sea state: above 1 vents the chamber through a by-pass (lower pressure, larger stroke); below 1 throttles it (smaller stroke, higher pressure).' },
    { key: 'survFactor', label: 'Load factor on the survival response', unit: '-', type: 'number', default: 1.35, min: 1, max: 2, group: 'Acceptance', role: 'criteria', help: `IEC TS 62600-2 normal design situation. ${CONFIRM}` },
    { key: 'nUnits', label: 'Chambers in the plant', type: 'int', default: 8, min: 1, max: 200, group: 'Plant', role: 'geometry' },
    { key: 'avail', label: 'Availability', unit: '%', type: 'number', default: 92, min: 30, max: 100, group: 'Plant', role: 'ic' },
    { key: 'wallConcrete', label: 'Structural concrete per chamber', unit: 'm³', type: 'number', default: 900, min: 0, max: 50000, group: 'Plant', role: 'geometry', help: 'For the embodied-carbon statement only (caisson or breakwater chamber walls).' },
    costInput('energyPrice', 'Value of delivered energy', 'USD/MWh', 150, 2000, 'Used only to express energy as revenue in the economic advice.'),
    { key: 'nOmega', label: 'Frequency bins', type: 'int', default: 200, min: 25, max: 1600, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nOmega', label: 'Frequency bins', levels: [50, 100, 200, 400], qoi: 'powerSea' },
  compute(inp) {
    const rho = inp.rhoW, a = inp.radius, d = inp.draft, la = inp.laFactor * a, patm = inp.patm * 1e3;   // kPa → Pa
    need(inp.depth > d + 1, `The column draft (${d} m) leaves no water under the mouth in ${inp.depth} m depth.`);
    const body = run(() => columnBody({ radius: a, draft: d, addedLength: la, depth: inp.depth, rho, excitationFactor: inp.fExc, Bvisc: inp.Bvisc * 1e3 }), 'Water column');
    const S = body.area, V0 = S * inp.airHeight;
    const C = inp.compress ? chamberCompliance({ volume: V0, patm }) : 0;
    const res = columnResonance({ draft: d, addedLength: la });
    const wp = 2 * Math.PI / inp.Tp, amp = inp.Hs / 2, wells = inp.turbine === 'wells';
    const reg0 = { column: body.column, area: S, omega: wp, amplitude: amp, compliance: C };
    // Wells damping: optimum for a regular wave at the design peak period, or the entered value
    const opt = run(() => optimalTurbineDamping({ ...reg0, depth: inp.depth, rho }), 'Optimum turbine damping');
    const Gt = wells ? (inp.dampMode === 'optimal' ? opt.conductance : 1 / inp.damping) : undefined;
    const turb = wells ? { turbine: 'wells', conductance: Gt } : { turbine: 'impulse', ki: inp.ki };
    const reg = run(() => owcRegular({ ...reg0, ...turb, efficiency: inp.etaTurb }), 'Regular-wave response');
    const J = regularWavePower({ amplitude: amp, omega: wp, depth: inp.depth, rho });
    const cw = J > 0 ? reg.pneumaticPower / J : 0;
    const colP = body.column(wp), cwLimit = G / (wp * wp);   // axisymmetric heave limit λ/2π in deep water

    // irregular sea
    const n = Math.round(inp.nOmega), om = Float64Array.from(linspace(0.15, 3.0, n));
    const spec = (Hs, Tp) => waveSpectrum({ Hs, Tp, gamma: inp.gamma, omega: om });
    const sea = run(() => owcIrregular({ column: body.column, area: S, omega: om, S: spec(inp.Hs, inp.Tp), compliance: C, ...turb, efficiency: inp.etaTurb }), 'Irregular-sea response');
    const Jsea = wavePower({ Hm0: inp.Hs, Te: 0.9 * inp.Tp, depth: inp.depth, rho });
    const mpm = (sigma, Tp, hours = 3) => sigma * Math.sqrt(2 * Math.log(Math.max(hours * 3600 / (TZ_OVER_TP * Tp), 2)));
    const pMax = mpm(sea.sigmaP, inp.Tp), xMax = mpm(sea.sigmaX, inp.Tp);
    const rated = inp.ratedKW * 1e3;

    // power matrix and annual energy
    const rows = (inp.scatter || []).filter((r) => r.hs > 0);
    need(rows.length >= 1, 'The scatter diagram needs at least one row with a positive Hs.');
    const hs = rows.map((r) => Number(r.hs)), prob = rows.map((r) => TPS.map((t) => Math.max(0, Number(r[`t${t}`]) || 0)));
    need(prob.flat().some((p) => p > 0), 'The scatter diagram has no occurrence above zero.');
    const pm = run(() => owcPowerMatrix({ hs, tp: TPS, column: body.column, area: S, compliance: C, ...turb, efficiency: inp.etaTurb, gamma: inp.gamma, ratedPower: rated, prob, availability: inp.avail / 100, omega: om }), 'Power matrix');
    const aepMWh = pm.aepWh * inp.nUnits / 1e6, cf = pm.capacityFactor;

    // extreme sea state with the survival valve setting (linear model: an upper bound on the motion)
    const surv = run(() => owcIrregular({ column: body.column, area: S, omega: om, S: spec(inp.Hs100, inp.Tp100), compliance: C, ...(wells ? { turbine: 'wells', conductance: Gt * inp.survValve } : { turbine: 'impulse', ki: inp.ki / (inp.survValve * inp.survValve) }), efficiency: inp.etaTurb }), 'Survival response');
    const xSurv = mpm(surv.sigmaX, inp.Tp100) * inp.survFactor, pSurv = mpm(surv.sigmaP, inp.Tp100) * inp.survFactor;

    const checks = [
      check('Chamber pressure in the operating sea (3-hour most probable maximum)', pMax / 1e3, inp.pRelief, 'kPa', 'IEC TS 62600-2:2019 §7 (design loads); linearised chamber model valid to about 0.2 p_atm',
        'Largest pressure amplitude expected in three hours of the design operating sea against the relief setting. Above it the relief valve opens (energy is vented) and the chamber shell and turbine see more than their design pressure.'),
      check('Column stroke in the operating sea against the lip submergence', xMax, d, 'm', 'Falcão & Henriques (2016) §3 (air venting below the lip)',
        'Largest fall of the internal water surface against the depth of the column mouth. If exceeded the mouth uncovers, air escapes under the lip and the pneumatic power collapses for that wave.'),
      check('Survival water rise against the chamber air height (factored)', xSurv, inp.airHeight, 'm', 'IEC TS 62600-2:2019 §6 (survival design situation); Rayleigh most probable maximum',
        'Factored largest rise of the internal water surface in the extreme sea against the height to the turbine duct. If exceeded, green water reaches the turbine and the chamber roof takes a slam load.'),
      check('Survival chamber pressure against the relief setting (factored)', pSurv / 1e3, inp.pRelief, 'kPa', 'IEC TS 62600-2:2019 §6; relief-valve sizing',
        'Factored largest chamber pressure in the extreme sea with the survival valve setting. Above the relief pressure the by-pass area is too small to protect the chamber and turbine.'),
    ];
    const gov = worst(checks);

    // curves
    const Ts = linspace(4, 16, 49);
    const curve = Ts.map((T) => { const w = 2 * Math.PI / T; const r = owcRegular({ ...reg0, omega: w, ...turb, efficiency: inp.etaTurb }); const j = regularWavePower({ amplitude: amp, omega: w, depth: inp.depth, rho }); return { T, cw: r.pneumaticPower / j, lim: G / (w * w), x: r.stroke / amp, p: r.pressure / amp / 1e3 }; });
    const gs = linspace(Math.log(opt.conductance / 20), Math.log(opt.conductance * 20), 41).map(Math.exp);
    const sweep = gs.map((g) => owcIrregular({ column: body.column, area: S, omega: om, S: spec(inp.Hs, inp.Tp), compliance: C, conductance: g }).pneumaticPower / 1e3);
    const revenue = aepMWh * inp.energyPrice, concrete = inp.wallConcrete * inp.nUnits;
    const tuned = Math.abs(res.period - inp.Tp) / inp.Tp;

    return {
      kpis: [
        { key: 'Tn', label: 'Natural period of the water column', value: round(res.period, 5), unit: 's', explain: 'Piston-mode resonance, 2π√((d + lₐ)/g). Capture is highest where this lies inside the band of energetic wave periods; a deeper lip lengthens it.' },
        { key: 'compliance', label: 'Air-chamber compliance', value: round(C * 1e3, 5), unit: 'm³/kPa', explain: 'Air volume change per unit pressure, V₀/(γ·p_atm). It acts as a spring in series with the turbine: a large chamber softens the air and lowers the power at short periods.' },
        { key: 'dampOpt', label: 'Optimum Wells damping p/q at the design period', value: round(opt.damping, 4), unit: 'Pa·s/m³', explain: 'Turbine pressure-to-flow ratio that maximises pneumatic power in a regular wave at the peak period. Select the turbine diameter and speed to give about this value.' },
        { key: 'dampUsed', label: 'Turbine damping used (equivalent linear)', value: round(1 / sea.conductance, 4), unit: 'Pa·s/m³', explain: 'Pressure over flow of the turbine in the design sea state: the entered or optimum value for a Wells turbine, the statistically linearised value for an impulse turbine.' },
        { key: 'powerReg', label: 'Pneumatic power, regular wave (H = Hs, T = Tp)', value: round(reg.pneumaticPower / 1e3, 5), unit: 'kW', explain: 'Mean power delivered by the air to the turbine in a regular wave of the design height and period.' },
        { key: 'cwRegular', label: 'Capture width, regular wave', value: round(cw, 5), unit: 'm', explain: `Pneumatic power over the incident power per metre of crest. The limit for an axisymmetric heaving absorber in deep water is λ/2π = ${fmt(cwLimit)} m at this period; the column is ${fmt(2 * a)} m wide.` },
        { key: 'cwRatio', label: 'Capture-width ratio (regular wave)', value: round(cw / (2 * a), 4), unit: '-', explain: 'Capture width over the column diameter. Above 1 the chamber absorbs more than the energy passing its own width (possible near resonance for a small device).' },
        { key: 'powerSea', label: 'Mean pneumatic power, design sea state', value: round(sea.pneumaticPower / 1e3, 6), unit: 'kW', explain: `Mean air power into the turbine in the irregular design sea, P = G·σp². Incident power is ${fmt(Jsea / 1e3)} kW per metre of crest, so the sea-state capture width is ${fmt(sea.pneumaticPower / Jsea)} m.` },
        { key: 'powerElec', label: 'Mean electrical power, design sea state', value: round(Math.min(sea.electricalPower, rated) / 1e3, 5), unit: 'kW', explain: 'Pneumatic power times the mean turbine–generator efficiency, capped at the rating.' },
        { key: 'pSig', label: 'Chamber pressure, 3-hour most probable maximum', value: round(pMax / 1e3, 4), unit: 'kPa', limit: inp.pRelief, sense: 'max', explain: `Largest gauge-pressure amplitude expected in three hours of the design sea (standard deviation ${fmt(sea.sigmaP / 1e3)} kPa), against the relief setting.` },
        { key: 'strokeMax', label: 'Column stroke, 3-hour most probable maximum', value: round(xMax, 4), unit: 'm', limit: d, sense: 'max', explain: 'Largest excursion of the internal water surface in the design sea, against the lip submergence: beyond it the mouth uncovers and the chamber vents.' },
        { key: 'aep', label: 'Net annual energy of the plant', value: round(aepMWh, 5), unit: 'MWh/yr', explain: `Power matrix times the scatter diagram for ${inp.nUnits} chamber(s) at ${inp.avail} % availability, each capped at ${inp.ratedKW} kW.` },
        { key: 'cf', label: 'Capacity factor', value: round(cf * 100, 4), unit: '%', explain: 'Mean electrical power over the rating. Wave plants usually land between 15 and 35 %; a low value means the rating is oversized for the resource, a high one that energy is being clipped.' },
        { key: 'strokeSurv', label: 'Survival water rise (factored)', value: round(xSurv, 4), unit: 'm', limit: inp.airHeight, sense: 'max', explain: 'Factored most probable maximum rise of the water column in the extreme sea with the survival valve setting, against the height to the turbine duct. Linear theory: an upper-bound screening value.' },
        { key: 'pSurv', label: 'Survival chamber pressure (factored)', value: round(pSurv / 1e3, 4), unit: 'kPa', limit: inp.pRelief, sense: 'max', explain: 'Factored largest chamber pressure in the extreme sea with the by-pass open to the stated extent.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Capture width against wave period', x: { label: 'Wave period', unit: 's' }, y: { label: 'Capture width', unit: 'm' },
          series: [{ name: 'This chamber and turbine', x: Ts.map((t) => round(t)), y: curve.map((c) => round(c.cw)) }, { name: 'Axisymmetric limit λ/2π', x: Ts.map((t) => round(t)), y: curve.map((c) => round(c.lim)), dash: true }],
          explain: 'Pneumatic capture width in regular waves. The peak sits at the natural period of the column; compare its position with the periods that carry the site energy and its height with the dashed theoretical limit.' },
        { type: 'line', title: 'Response per metre of wave amplitude', x: { label: 'Wave period', unit: 's' }, y: { label: 'Stroke (m/m) and pressure (kPa/m)', unit: '' },
          series: [{ name: 'Column stroke per unit amplitude', x: Ts.map((t) => round(t)), y: curve.map((c) => round(c.x)) }, { name: 'Chamber pressure per unit amplitude (kPa/m)', x: Ts.map((t) => round(t)), y: curve.map((c) => round(c.p)) }],
          explain: 'Transfer functions of the water-column motion and the chamber pressure. A stroke ratio above 1 marks resonance; the pressure curve shows which periods load the chamber and turbine.' },
        { type: 'line', title: 'Mean pneumatic power against Wells turbine damping (design sea state)', x: { label: 'Turbine damping p/q', unit: 'Pa·s/m³' }, y: { label: 'Pneumatic power', unit: 'kW' },
          series: [{ name: 'Wells turbine', x: gs.map((g) => round(1 / g)).reverse(), y: sweep.map((v) => round(v)).reverse() }],
          explain: 'Power in the irregular design sea as the turbine damping is varied. The flat top shows how tolerant the plant is to off-optimum damping; too little damping gives large motion with little pressure, too much stalls the column.' },
        { type: 'heatmap', title: 'Electrical power matrix per chamber', xs: TPS, ys: hs.map((h) => round(h)), z: pm.power.map((r) => r.map((p) => round(p / 1e3, 4))), zLabel: 'Power', zUnit: 'kW', x: { label: 'Peak period Tp', unit: 's' }, y: { label: 'Significant wave height Hs', unit: 'm' },
          explain: 'Mean electrical power in each sea state, capped at the rating. Overlay it mentally on the scatter diagram: the plant should produce most where the sea spends most time.' },
      ],
      tables: [
        { title: 'Power matrix and occurrence', columns: [{ key: 'hs', label: 'Hs', unit: 'm' }, ...TPS.flatMap((t) => [{ key: `p${t}`, label: `P at Tp ${t} s`, unit: 'kW' }, { key: `o${t}`, label: `Occurrence at ${t} s`, unit: '%' }])],
          rows: hs.map((h, i) => { const tot = prob.flat().reduce((s, v) => s + v, 0); return { hs: round(h), ...Object.fromEntries(TPS.flatMap((t, j) => [[`p${t}`, round(pm.power[i][j] / 1e3, 4)], [`o${t}`, round(100 * prob[i][j] / tot, 3)]])) }; }),
          explain: 'One row per wave-height bin: electrical power per chamber and the share of the year in each sea state. Their product summed over the table is the annual energy.' },
        { title: 'Chamber and turbine summary', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['Column area', S, 'm²'], ['Added length', la, 'm'], ['Chamber air volume', V0, 'm³'], ['Column mass incl. added length', colP.M / 1e3, 't'], ['Radiation damping at Tp', colP.B / 1e3, 'kN·s/m'], ['Excitation per unit amplitude at Tp', colP.Fe / 1e3, 'kN/m'], ['Regular-wave stroke amplitude', reg.stroke, 'm'], ['Regular-wave pressure amplitude', reg.pressure / 1e3, 'kPa'], ['Regular-wave flow amplitude', reg.flow, 'm³/s'], ['Design-sea flow, standard deviation', sea.sigmaQ, 'm³/s'], ['Survival stroke, standard deviation', surv.sigmaX, 'm']].map(([q, v, u]) => ({ q, v: round(v, 5), u })),
          explain: 'Intermediate quantities for checking by hand and for sizing the turbine (flow and pressure amplitudes give its diameter and speed).' },
      ],
      notes: [
        `Natural period ${fmt(res.period)} s against a design peak period of ${fmt(inp.Tp)} s (${fmt(tuned * 100)} % apart); ${gov.name.toLowerCase()} governs at ${fmt(gov.uc)}.`,
        wells ? `Wells turbine damping ${fmt(1 / Gt)} Pa·s/m³ (${inp.dampMode === 'optimal' ? 'optimum for the design period' : `entered; the optimum is ${fmt(opt.damping)}`}).` : `Impulse turbine linearised to ${fmt(1 / sea.conductance)} Pa·s/m³ in the design sea (optimum linear damping would be ${fmt(opt.damping)}).`,
        'Extreme-sea values come from the linear model with the survival valve setting and overstate resonant motion.',
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name}: ${fmt(gov.demand)} against ${fmt(gov.capacity)} ${gov.unit} (utilisation ${fmt(gov.uc)}). Capture width ${fmt(cw)} m in the design regular wave, ${fmt(cw / cwLimit * 100)} % of the axisymmetric limit.`,
          gov.uc > 1 ? (gov === checks[0] || gov === checks[3] ? `Lower the turbine damping or enlarge the by-pass (survival conductance ratio above ${fmt(inp.survValve)}), or raise the design pressure to about ${fmt(gov.demand * 1.1)} kPa.` : gov === checks[1] ? `Deepen the lip to about ${fmt(xMax * 1.15)} m or raise the turbine damping to limit the stroke.` : `Raise the chamber air height to about ${fmt(xSurv * 1.1)} m, or throttle the chamber in survival mode.`)
            : tuned > 0.25 ? `The column is tuned ${fmt(tuned * 100)} % away from the design period: a lip submergence of about ${fmt(Math.max(G * (inp.Tp / (2 * Math.PI)) ** 2 - la, 1))} m would bring the resonance onto it (check the stroke and the structure cost first).` : 'The chamber is tuned to the design sea and within its limits: confirm the entry losses and the excitation factor with model tests, and the turbine efficiency curve (stall) with the supplier.',
          'IEC TS 62600-2:2019 (design situations) and IEC TS 62600-100:2012 (power performance); Falcão & Henriques (2016) for chamber and turbine matching.',
          `${concrete > 0 ? `${fmt(concrete)} m³ of chamber concrete ≈ ${fmt(concrete * CO2_CONCRETE)} t CO₂e (screening factor ${CO2_CONCRETE} t/m³), repaid in about ${fmt(concrete * CO2_CONCRETE / Math.max(aepMWh * 0.4, 1e-9))} years of production at 0.4 t CO₂/MWh displaced. ` : ''}Building the chambers into a breakwater that is needed anyway shares the concrete and the seabed footprint between coastal protection and generation.`),
        adv('economic', cf < 0.15 ? 'warn' : 'info', `The plant yields ${fmt(aepMWh)} MWh a year at a capacity factor of ${fmt(cf * 100)} %; at the indicative ${fmt(inp.energyPrice)} USD/MWh that is ${money(revenue)} a year, and one point of turbine efficiency is worth ${money(revenue * 0.01 / inp.etaTurb)} a year.`,
          cf < 0.15 ? `The rating of ${inp.ratedKW} kW per chamber is high for this resource: a turbine–generator of about ${fmt(pm.meanPower / (inp.avail / 100) / 0.25 / 1e3)} kW would run at 25 % capacity factor with little lost energy and lower electrical cost.` : 'Size the turbine and generator on energy captured per unit of rating, not on peak power: clipping the rare large seas costs little energy and saves generator, converter and cable capacity.',
          'LCOE convention (IEA-OES / IRENA): annual net energy as denominator; IEC TS 62600-100 power matrix × IEC TS 62600-101 resource for the energy estimate.',
          energyCo2Text(aepMWh)),
      ],
      publish: { 'owc.ratedPower': rated * inp.nUnits, 'owc.aep': pm.aepWh * inp.nUnits, 'owc.capacityFactor': cf, 'owc.turbineDamping': 1 / sea.conductance },
      quantities: { concrete_m3: round(concrete) },
    };
  },
  verify: [
    { // T = 2π√((d + lₐ)/g) with d = 4 m, lₐ = 1 m (radius 3 m × 1/3): 2π√(5/9.80665) = 4.4865 s.
      // C = V₀/(γ p_atm): V₀ = π·3²·17.684 = 500.0 m³ → 500/(1.4·101325) = 3.5247e-3 m³/Pa = 3.5247 m³/kPa.
      name: 'Water-column natural period T = 2π√((d + lₐ)/g) and isentropic chamber compliance V₀/(γ·p_atm)', source: 'Evans (1978); Sarmento & Falcão (1985); hand calculation',
      inputs: { radius: 3, draft: 4, laFactor: 1 / 3, airHeight: 17.684, Tp: 6, Hs: 1 },
      expect: [{ kpi: 'Tn', value: 4.4865, tol: 2e-4 }, { kpi: 'compliance', value: 3.5247, tol: 5e-4 }] },
    { // At resonance with optimum damping, no losses and incompressible air, a heaving axisymmetric absorber captures
      // the energy of a crest length λ/2π = g/ω₀² (Budal & Falnes; Evans 1976). With ω₀² = g/(d + lₐ) this is d + lₐ = 5 m.
      name: 'Optimum damping at resonance reaches the axisymmetric limit: capture width = λ/2π = d + lₐ', source: 'Evans (1976) J. Fluid Mech. 77:1–25; Falnes (2002) Ocean Waves and Oscillating Systems §6.1',
      inputs: { radius: 3, draft: 4, laFactor: 1 / 3, compress: false, Bvisc: 0, depth: 400, Tp: 2 * Math.PI * Math.sqrt(5 / 9.80665), Hs: 1, turbine: 'wells', dampMode: 'optimal' },
      expect: [{ kpi: 'cwRegular', value: 5.0, tol: 2e-3 }] },
  ],
};
