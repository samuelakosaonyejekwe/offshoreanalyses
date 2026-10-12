// Wave-energy converter hydrodynamics: heaving point absorber in the frequency domain — hydrodynamic coefficients
// (truncated-cylinder approximations or an imported table), response and absorbed power in regular and irregular
// waves, capture width against the theoretical limits, the power matrix with the site scatter diagram and the
// survival response in the extreme sea state.
import { dispersion, wavePower, waveSpectrum, spectralPeriods, scatterResource, cylinderHeave, pointAbsorber, captureWidth, heaveLimits, cylinderBody, spectralResponse, powerMatrix, matrixEnergy, survivalResponse } from '../../engine/aero/wec.js';
import { site, adv, check, worst, sev, fmt, round, run, need, linspace, interp, G, LIMIT_WEC, LIMIT_CODE, CONFIRM, steelCo2Text, energyCo2Text } from './_common.js';

const TPS = [6, 8, 10, 12, 14];
// Indicative exposed North-Atlantic scatter diagram (occurrence %, rows Hs, columns Tp): replace with the site hindcast.
const SCATTER = [[0.5, 4, 6, 3, 1, 0], [1.5, 6, 14, 10, 4, 1], [2.5, 1, 9, 11, 5, 1], [3.5, 0, 3, 7, 4, 1], [4.5, 0, 1, 3, 2, 1], [5.5, 0, 0, 1, 1, 0]];

export default {
  id: 's16-wec-hydrodynamics', suite: 16, title: 'Wave-energy converter hydrodynamics (point absorber)',
  summary: 'Heave response, absorbed power and capture width of a point absorber in regular and irregular waves, its power matrix and annual energy, and the survival stroke in the extreme sea state.',
  covers: [2], fidelity: 'analytical',
  standards: ['IEC TS 62600-100:2012 (power performance of wave energy converters)', 'IEC TS 62600-101:2015 (wave energy resource assessment)', 'IEC TS 62600-2:2019 (design requirements)', 'DNV-RP-C205 (2021) §3.5.5 (wave spectra)', 'Falnes (2002) Ocean Waves and Oscillating Systems'],
  models: ['Linear wave radiation and diffraction (single-body, frequency domain)', 'Hydrodynamic added mass and radiation damping (Haskind relation)', 'Froude–Krylov excitation on a truncated cylinder', 'Airy wave theory and linear dispersion', 'JONSWAP / Bretschneider spectra', 'Optimal passive damping and reactive (complex-conjugate) control', 'Optimal energy absorption limits (Budal & Falnes)', 'Power matrix × scatter diagram', 'Rayleigh most-probable-maximum survival response'],
  limits: `${LIMIT_WEC} One rigid body in heave with a linear power take-off; viscous loss is a linearised damping; array interaction is a user factor. Linear theory overstates resonant motion in steep and extreme seas (no nonlinear Froude–Krylov force, slamming or overtopping here — the time-domain PTO module adds end stops and nonlinear PTO). ${LIMIT_CODE}`,
  publishes: ['wec.ratedPower', 'wec.aep', 'wec.BptoOpt', 'wec.capacityFactor', 'wec.maxStroke'],
  searchTerms: 'point absorber wave energy converter capture width power matrix Budal limit',
  inputs: [
    { key: 'radius', label: 'Float radius', unit: 'm', type: 'number', default: 5, min: 0.5, max: 20, group: 'Device geometry', role: 'geometry', design: true },
    { key: 'draft', label: 'Float draft', unit: 'm', type: 'number', default: 6, min: 0.3, max: 40, group: 'Device geometry', role: 'geometry', design: true, help: 'Initial hydrostatic equilibrium: the float sits at this draft in still water.' },
    { key: 'massRatio', label: 'Body mass / displaced mass', unit: '-', type: 'number', default: 1.5, min: 0.5, max: 3, group: 'Device geometry', role: 'ic', help: '1 = freely floating at the draft. Above 1 represents supplementary inertia (a submerged mass or an inerter) used to tune the natural period to the waves.' },
    { key: 'coefSource', label: 'Hydrodynamic coefficients', type: 'select', default: 'cylinder', group: 'Hydrodynamics', role: 'bc', options: [{ value: 'cylinder', label: 'Truncated-cylinder approximations' }, { value: 'table', label: 'Imported coefficient table' }],
      help: 'Radiation and excitation boundary condition. Import the heave added mass, radiation damping and excitation force from a diffraction–radiation analysis or basin test when available.' },
    { key: 'coefTable', label: 'Heave coefficients', type: 'table', group: 'Hydrodynamics', role: 'bc', showIf: { key: 'coefSource', in: ['table'] },
      columns: [{ key: 'omega', label: 'Frequency', unit: 'rad/s', type: 'number' }, { key: 'A', label: 'Added mass', unit: 't', type: 'number' }, { key: 'B', label: 'Radiation damping', unit: 'kN·s/m', type: 'number' }, { key: 'Fe', label: 'Excitation per unit amplitude', unit: 'kN/m', type: 'number' }],
      default: [], help: 'At least four ascending frequencies covering 0.2–2.5 rad/s. Values outside the table are held at the end rows.' },
    { key: 'caAdd', label: 'Added-mass coefficient c in A = c·ρ·D³', unit: '-', type: 'number', default: 0.1667, min: 0.05, max: 0.4, group: 'Hydrodynamics', role: 'bc', showIf: { key: 'coefSource', in: ['cylinder'] }, uncertain: { dist: 'normal', cov: 0.15 }, help: '1/6 = half the added mass of a disc (McCormick 1981); π/24 for a hemisphere.' },
    { key: 'fDiff', label: 'Diffraction factor on the Froude–Krylov force', unit: '-', type: 'number', default: 1, min: 0.5, max: 1.5, group: 'Hydrodynamics', role: 'bc', showIf: { key: 'coefSource', in: ['cylinder'] }, uncertain: { dist: 'normal', cov: 0.1 } },
    { key: 'Bvisc', label: 'Linearised viscous damping', unit: 'kN·s/m', type: 'number', default: 20, min: 0, max: 2000, group: 'Hydrodynamics', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.4 }, help: 'Equivalent linear damping of vortex shedding at the float edges, (8/3π)·½ρCdA|v|; calibrate against decay tests.' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 50, min: 10, max: 500, group: 'Site', role: 'environment', bind: 'site.depth' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 995, max: 1035, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'Hs', label: 'Design operating significant wave height', unit: 'm', type: 'number', default: 2.5, min: 0.2, max: 8, group: 'Operating sea state', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Incident-wave boundary condition for the headline response.' },
    { key: 'Tp', label: 'Design operating peak period', unit: 's', type: 'number', default: 9, min: 3, max: 18, group: 'Operating sea state', role: 'bc' },
    { key: 'gamma', label: 'JONSWAP peak-enhancement factor', unit: '-', type: 'number', default: 2, min: 1, max: 7, group: 'Operating sea state', role: 'bc', help: '1 = Bretschneider (open-ocean swell-dominated), 3.3 = fetch-limited wind sea.' },
    { key: 'scatter', label: 'Scatter diagram (occurrence, any scale)', type: 'table', group: 'Wave resource', role: 'bc',
      columns: [{ key: 'hs', label: 'Hs', unit: 'm', type: 'number' }, ...TPS.map((t) => ({ key: `t${t}`, label: `Tp ${t} s`, unit: '%', type: 'number' }))],
      default: SCATTER.map(([hs, ...p]) => ({ hs, ...Object.fromEntries(TPS.map((t, i) => [`t${t}`, p[i]])) })), help: 'Joint occurrence of Hs and Tp (rows: bin-centre Hs). Normalised internally.' },
    site('Hs100', { group: 'Survival sea state' }), site('Tp100', { group: 'Survival sea state' }),
    { key: 'control', label: 'Power take-off control', type: 'select', default: 'passive', group: 'Power take-off', role: 'control', options: [{ value: 'passive', label: 'Optimal passive damping per sea state' }, { value: 'reactive', label: 'Reactive (complex-conjugate) control in regular waves' }, { value: 'fixed', label: 'Fixed PTO damping' }],
      help: 'Reactive control is evaluated for the regular-wave capture width only; irregular-sea results use passive damping tuned to the peak period (or the fixed value).' },
    { key: 'Bpto', label: 'Fixed PTO damping', unit: 'kN·s/m', type: 'number', default: 600, min: 1, max: 20000, group: 'Power take-off', role: 'bc', showIf: { key: 'control', in: ['fixed'] }, design: true },
    { key: 'strokeLim', label: 'Available stroke (amplitude)', unit: 'm', type: 'number', default: 3.5, min: 0.2, max: 1000, group: 'Power take-off', role: 'bc', design: true, help: 'End-stop limit of the PTO, each way from mid-stroke.' },
    { key: 'forceLim', label: 'PTO force limit (amplitude)', unit: 'kN', type: 'number', default: 2500, min: 10, max: 1e6, group: 'Power take-off', role: 'bc', design: true },
    { key: 'survDamp', label: 'Survival-mode PTO damping / operating damping', unit: '-', type: 'number', default: 10, min: 0, max: 200, group: 'Power take-off', role: 'control', help: 'Braking applied in the extreme sea state: 1 = operating setting, large values approach a locked float, 0 = PTO disconnected.' },
    { key: 'ratedKW', label: 'Rated power per device', unit: 'kW', type: 'number', default: 150, min: 5, max: 5000, group: 'Power take-off', role: 'geometry', design: true },
    { key: 'effPto', label: 'PTO conversion efficiency (absorbed to electrical)', unit: '-', type: 'number', default: 0.8, min: 0.3, max: 1, group: 'Power take-off', role: 'operation', uncertain: { dist: 'normal', cov: 0.05 } },
    { key: 'nUnits', label: 'Devices in the array', type: 'int', default: 10, min: 1, max: 500, group: 'Array', role: 'geometry' },
    { key: 'qArray', label: 'Array interaction factor q', unit: '-', type: 'number', default: 0.95, min: 0.5, max: 1.3, group: 'Array', role: 'operation', help: 'Array power over N isolated devices; 0.9–1.0 for typical spacings.' },
    { key: 'avail', label: 'Availability', unit: '%', type: 'number', default: 90, min: 30, max: 100, group: 'Array', role: 'ic' },
    { key: 'survFactor', label: 'Load factor on the survival response', unit: '-', type: 'number', default: 1.35, min: 1, max: 2, group: 'Acceptance', role: 'criteria', help: `IEC TS 62600-2 normal design situation. ${CONFIRM}` },
    { key: 'nOmega', label: 'Frequency bins', type: 'int', default: 240, min: 30, max: 1920, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nOmega', label: 'Frequency bins', levels: [60, 120, 240, 480], qoi: 'powerSea' },
  compute(inp) {
    const rho = inp.rhoW, a = inp.radius, d = inp.draft;
    need(inp.depth > d + 2, `The draft (${d} m) leaves no water under the float in ${inp.depth} m depth.`);
    const disp = rho * Math.PI * a * a * d, mass = disp * inp.massRatio, Kh = rho * G * Math.PI * a * a;
    const Bv = inp.Bvisc * 1e3;
    let body, coefAt;
    if (inp.coefSource === 'table') {
      const rows = (inp.coefTable || []).filter((r) => r.omega > 0);
      need(rows.length >= 4, 'The coefficient table needs at least four rows (frequency, added mass, radiation damping, excitation).');
      for (let i = 1; i < rows.length; i++) need(rows[i].omega > rows[i - 1].omega, `Coefficient table: frequencies must ascend (row ${i + 1}).`);
      need(rows.every((r) => r.B >= 0 && r.Fe >= 0 && r.A > -mass / 1e3), 'Coefficient table: damping and excitation must be non-negative and the added mass physical.');
      const w = rows.map((r) => r.omega), A = rows.map((r) => r.A * 1e3), Bt = rows.map((r) => r.B * 1e3), F = rows.map((r) => r.Fe * 1e3);
      coefAt = (om) => ({ A: interp(w, A, om), B: interp(w, Bt, om), Fe: interp(w, F, om), K: Kh });
      body = (om) => { const c = coefAt(om); return { M: mass + c.A, B: c.B + Bv, K: Kh, Fe: c.Fe }; };
    } else {
      const geom = { radius: a, draft: d, depth: inp.depth, rho, addedMassCoeff: inp.caAdd, diffractionFactor: inp.fDiff };
      coefAt = (om) => cylinderHeave({ ...geom, omega: om });
      body = cylinderBody({ ...geom, mass, Bvisc: Bv });
    }
    const wp = 2 * Math.PI / inp.Tp, cp = coefAt(wp);
    const Tn = 2 * Math.PI * Math.sqrt((mass + cp.A) / Kh);
    const omega = linspace(0.15, 3.0, inp.nOmega);
    const S = waveSpectrum({ Hs: inp.Hs, Tp: inp.Tp, gamma: inp.gamma, omega });
    const per = spectralPeriods({ omega, S });
    const J = wavePower({ Hm0: inp.Hs, Te: per.Te, depth: inp.depth, rho });
    // PTO damping for the sea state: passive optimum at the peak frequency, or the fixed value
    const bp = body(wp);
    const BptoOpt = Math.hypot(bp.B, wp * bp.M - bp.K / wp);
    const Bsea = inp.control === 'fixed' ? inp.Bpto * 1e3 : BptoOpt;
    const sea = run(() => spectralResponse({ omega, S, body, Bpto: Bsea }), 'Irregular-sea response');
    const powerSea = sea.power;
    const cw = captureWidth({ power: powerSea, J, width: 2 * a });
    // regular wave of equal energy at the peak period (H = Hs/√2) with the selected control and the PTO limits
    const amp = inp.Hs / (2 * Math.SQRT2);
    const lim = { stroke: inp.strokeLim, force: inp.forceLim * 1e3 };
    const reg = run(() => pointAbsorber({ mass, A: cp.A, B: cp.B, K: Kh, Fe: cp.Fe, omega: wp, amplitude: amp, control: inp.control, Bpto: inp.Bpto * 1e3, Bvisc: Bv, limits: lim }), 'Regular-wave response');
    const k = dispersion(wp, inp.depth);
    const Jreg = 0.5 * rho * G * amp * amp * k.cg;
    const cwReg = reg.power / Jreg;
    const hl = heaveLimits({ omega: wp, H: 2 * amp, volume: Math.PI * a * a * d, depth: inp.depth, rho });
    // power matrix and annual energy
    const hsBins = (inp.scatter || []).map((r) => r.hs);
    need(hsBins.length >= 1 && hsBins.every((h) => h > 0), 'The scatter diagram needs at least one row with a positive Hs.');
    const prob = inp.scatter.map((r) => TPS.map((t) => Math.max(0, Number(r[`t${t}`]) || 0)));
    need(prob.flat().some((p) => p > 0), 'The scatter diagram has no occurrence.');
    const omM = linspace(0.15, 3.0, Math.min(inp.nOmega, 240));
    const pm = run(() => powerMatrix({ hs: hsBins, tp: TPS, body, Bpto: inp.control === 'fixed' ? inp.Bpto * 1e3 : 'optimal', gamma: inp.gamma, ratedPower: inp.ratedKW * 1e3, efficiency: inp.effPto, omega: omM }), 'Power matrix');
    const en = matrixEnergy({ power: pm.power, prob, ratedPower: inp.ratedKW * 1e3, availability: inp.avail / 100, arrayFactor: inp.qArray, devices: inp.nUnits });
    // resource: Te ≈ 0.86–0.90 Tp depending on γ; use the ratio of the design spectrum
    const teRatio = per.Te / inp.Tp;
    const res = scatterResource({ hs: hsBins, te: TPS.map((t) => t * teRatio), prob, depth: inp.depth, rho });
    // survival in the extreme sea with the survival-mode damping (braked or locked PTO)
    const Sx = waveSpectrum({ Hs: inp.Hs100, Tp: inp.Tp100, gamma: 3.3, omega });
    const surv = run(() => survivalResponse({ omega, S: Sx, body, Bpto: Bsea * inp.survDamp, duration: 10800, strokeLimit: inp.strokeLim }), 'Survival response');
    const sigStroke = 2 * sea.sigmaX, sigForce = 2 * sea.sigmaF;
    const checks = [
      check('Significant stroke in the operating sea state', sigStroke, inp.strokeLim, 'm', 'IEC TS 62600-2:2019 §7 (serviceability of the PTO); Rayleigh statistics', 'Significant stroke amplitude (2σ) against the available stroke. If it is exceeded the float reaches the end stops in normal operation: lost energy and impact loads on the PTO.'),
      check('Significant PTO force in the operating sea state', sigForce / 1e3, inp.forceLim, 'kN', 'IEC TS 62600-2:2019; PTO rating', 'Significant PTO force amplitude (2σ) against the force rating. Above it the PTO saturates and the linear power estimate is optimistic.'),
      check('Survival stroke (3-h most probable maximum, factored)', surv.maxStroke * inp.survFactor, inp.strokeLim, 'm', 'IEC TS 62600-2:2019 (ULS, extreme sea state); DNV-RP-C205 §3.5', 'Factored most probable maximum stroke in the extreme sea against the available stroke. Above 1 the device needs a survival mode (locking, submergence or detuning) or end stops designed for the impact.'),
    ];
    const gov = worst(checks);
    // curves versus period
    const Ts = linspace(4, 16, 49);
    const rao = Ts.map((T) => { const om = 2 * Math.PI / T, b = body(om); return b.Fe / Math.hypot(b.K - om * om * b.M, om * (b.B + Bsea)); });
    const raoFree = Ts.map((T) => { const om = 2 * Math.PI / T, b = body(om); return Math.min(b.Fe / Math.hypot(b.K - om * om * b.M, om * b.B), 20); });
    const cwPass = [], cwReact = [], cwLim = [];
    for (const T of Ts) {
      const om = 2 * Math.PI / T, c = coefAt(om), kk = dispersion(om, inp.depth), Jr = 0.5 * rho * G * amp * amp * kk.cg;
      const base = { mass, A: c.A, B: c.B, K: Kh, Fe: c.Fe, omega: om, amplitude: amp, Bvisc: Bv, limits: lim };
      cwPass.push(pointAbsorber({ ...base, control: 'passive' }).power / Jr);
      cwReact.push(pointAbsorber({ ...base, control: 'reactive' }).power / Jr);
      cwLim.push(1 / kk.k);
    }
    const wc = linspace(0.3, 2.5, 45).map((om) => ({ om, ...coefAt(om) }));
    const respS = Array.from(omega, (om, i) => { const b = body(om); const H = b.Fe / Math.hypot(b.K - om * om * b.M, om * (b.B + Bsea)); return H * H * S[i]; });
    const idx = linspace(0, omega.length - 1, Math.min(omega.length, 160)).map(Math.round);
    const steelT = (2 * Math.PI * a * a + 2 * Math.PI * a * (d + 2)) * 0.015 * 7850 * 1.5 / 1e3;   // 15 mm shell with 50 % stiffening, per float
    const aepMwh = en.aepWh / 1e6;
    return {
      kpis: [
        { key: 'Tn', label: 'Heave natural period', value: round(Tn, 5), unit: 's', explain: 'Undamped natural period 2π√((m + A)/ρgA_w). Absorption peaks when it is close to the wave period; small floats are usually stiffer (shorter period) than the waves.' },
        { key: 'powerSea', label: 'Mean absorbed power in the operating sea state', value: round(powerSea / 1e3, 5), unit: 'kW', explain: 'Mean hydrodynamic power absorbed by the PTO damping in the design sea state (before conversion losses and the rating cap).' },
        { key: 'cwr', label: 'Capture-width ratio in the operating sea state', value: round(cw.ratio, 5), unit: '-', explain: 'Absorbed power over the wave power passing a width equal to the float diameter. Point absorbers reach 0.1–0.4 in irregular seas.' },
        { key: 'J', label: 'Incident wave power in the operating sea state', value: round(J / 1e3, 5), unit: 'kW/m', explain: 'Energy flux per metre of wave crest, ρg·Hs²·c_g/16 at the energy period.' },
        { key: 'sigStroke', label: 'Significant stroke amplitude', value: round(sigStroke, 4), unit: 'm', limit: inp.strokeLim, sense: 'max', explain: 'Twice the standard deviation of heave relative to the reference in the operating sea state.' },
        { key: 'sigForce', label: 'Significant PTO force amplitude', value: round(sigForce / 1e3, 5), unit: 'kN', limit: inp.forceLim, sense: 'max', explain: 'Twice the standard deviation of the PTO force: the load the PTO sees routinely.' },
        { key: 'BptoOpt', label: 'Optimal passive PTO damping at the peak period', value: round(BptoOpt / 1e3, 5), unit: 'kN·s/m', explain: 'Damping that maximises absorbed power for a sinusoidal wave at the peak period, √(B² + (ωM − K/ω)²).' },
        { key: 'cwRegular', label: 'Capture width in the equal-energy regular wave', value: round(cwReg, 5), unit: 'm', explain: `Absorbed power over the flux of a regular wave of period Tp and height Hs/√2 with the selected control and PTO limits${reg.constrained ? ' (limits active)' : ''}. The theoretical maximum for heave is λ/2π.` },
        { key: 'cwLimit', label: 'Heave capture-width limit λ/2π', value: round(hl.maxCaptureWidth, 5), unit: 'm', explain: 'Largest capture width an axisymmetric heaving body can reach at this period (Budal & Falnes; Evans; Newman).' },
        { key: 'budal', label: 'Budal upper bound at the peak period', value: round(hl.budalLimit / 1e3, 5), unit: 'kW', explain: 'Volume-limited bound (π/4)ρgVH/T: a float of this size cannot absorb more whatever the control.' },
        { key: 'aep', label: 'Annual electrical energy of the array', value: round(en.aepWh / 1e9, 5), unit: 'GWh/yr', explain: 'Power matrix times scatter diagram, with PTO efficiency, rating cap, array interaction and availability.' },
        { key: 'cf', label: 'Capacity factor', value: round(en.capacityFactor * 100, 4), unit: '%', explain: 'Mean electrical power over rated power. Wave devices typically achieve 20–35 %; a high value suggests the rating is too small for the site.' },
        { key: 'resource', label: 'Mean annual wave power of the site', value: round(res.meanJ / 1e3, 5), unit: 'kW/m', explain: 'Occurrence-weighted energy flux from the scatter diagram. Sites above 20–25 kW/m are considered commercially attractive.' },
        { key: 'survStroke', label: 'Survival stroke, 3-h most probable maximum', value: round(surv.maxStroke, 4), unit: 'm', explain: 'Linear estimate of the largest heave excursion in the extreme sea state, σ√(2 ln N). Linear theory is conservative here but shows whether a survival mode is needed.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Heave response amplitude operator', x: { label: 'Wave period', unit: 's' }, y: { label: 'Heave / wave amplitude', unit: 'm/m' },
          series: [{ name: 'With PTO damping', x: Ts.map((v) => round(v, 3)), y: rao.map((v) => round(v, 4)) }, { name: 'PTO disconnected (capped at 20)', x: Ts.map((v) => round(v, 3)), y: raoFree.map((v) => round(v, 4)) }, { name: 'Natural period', x: [Tn, Tn], y: [0, Math.max(...rao)], dash: true }],
          explain: 'Heave motion per unit wave amplitude against period. The peak sits at the natural period; the PTO damping flattens it. Long waves give a ratio of 1 (the float follows the surface).' },
        { type: 'line', title: 'Capture width versus wave period (regular waves)', x: { label: 'Wave period', unit: 's' }, y: { label: 'Capture width', unit: 'm' },
          series: [{ name: 'Passive damping', x: Ts.map((v) => round(v, 3)), y: cwPass.map((v) => round(v, 4)) }, { name: 'Reactive control (with limits)', x: Ts.map((v) => round(v, 3)), y: cwReact.map((v) => round(v, 4)) }, { name: 'Limit λ/2π', x: Ts.map((v) => round(v, 3)), y: cwLim.map((v) => round(v, 4)), dash: true }, { name: 'Float diameter', x: [Ts[0], Ts[Ts.length - 1]], y: [2 * a, 2 * a], dash: true }],
          explain: 'Capture width with optimal passive damping and with reactive control under the stroke and force limits, against the theoretical limit. The gap between passive and reactive is what advanced control could recover.' },
        { type: 'heatmap', title: 'Electrical power matrix per device', xs: TPS, ys: hsBins, z: pm.power.map((row) => row.map((v) => round(v / 1e3, 4))), zLabel: 'Power', zUnit: 'kW', x: { label: 'Peak period', unit: 's' }, y: { label: 'Hs', unit: 'm' },
          explain: 'Mean electrical power in each sea state (IEC TS 62600-100 format), capped at the rating. Multiply by the scatter diagram to see which sea states deliver the annual energy.' },
        { type: 'line', title: 'Hydrodynamic coefficients', x: { label: 'Frequency', unit: 'rad/s' }, y: { label: 'Normalised', unit: '-' },
          series: [{ name: 'Added mass / displaced mass', x: wc.map((q) => round(q.om, 4)), y: wc.map((q) => round(q.A / disp, 5)) }, { name: 'Radiation damping / (ρgA_w/ω)', x: wc.map((q) => round(q.om, 4)), y: wc.map((q) => round(q.B * q.om / Kh, 5)) }, { name: 'Excitation / ρgA_w', x: wc.map((q) => round(q.om, 4)), y: wc.map((q) => round(q.Fe / Kh, 5)) }],
          explain: 'Heave added mass, radiation damping and excitation force against frequency, normalised. Excitation tends to the hydrostatic value at low frequency and decays as the float draft shields the bottom.' },
        { type: 'line', title: 'Wave and heave response spectra', x: { label: 'Frequency', unit: 'rad/s' }, y: { label: 'Spectral density', unit: 'm²·s/rad' },
          series: [{ name: 'Wave', x: idx.map((i) => round(omega[i], 4)), y: idx.map((i) => round(S[i], 5)) }, { name: 'Heave response', x: idx.map((i) => round(omega[i], 4)), y: idx.map((i) => round(respS[i], 5)) }],
          explain: 'Design sea-state spectrum and the heave response spectrum. Power is absorbed where the two overlap; a response peak far from the wave peak means the device is detuned.' },
      ],
      tables: [
        { title: 'Power matrix and occurrence', columns: [{ key: 'hs', label: 'Hs', unit: 'm' }, ...TPS.flatMap((t) => [{ key: `p${t}`, label: `P at Tp ${t} s`, unit: 'kW' }, { key: `o${t}`, label: 'occurrence', unit: '%' }])],
          rows: hsBins.map((h, i) => ({ hs: h, ...Object.fromEntries(TPS.flatMap((t, j) => [[`p${t}`, round(pm.power[i][j] / 1e3, 4)], [`o${t}`, round(prob[i][j] / prob.flat().reduce((x, y) => x + y, 0) * 100, 3)]])) })),
          explain: 'Electrical power per device and occurrence for each Hs–Tp cell. Cells with high occurrence and mid-range power carry the yield; the highest cells are capped by the rating.' },
        { title: 'Response summary', columns: [{ key: 'q', label: 'Quantity' }, { key: 'op', label: 'Operating sea' }, { key: 'sv', label: 'Extreme sea' }, { key: 'u', label: 'Unit' }],
          rows: [{ q: 'Hs / Tp', op: `${inp.Hs} / ${inp.Tp}`, sv: `${inp.Hs100} / ${inp.Tp100}`, u: 'm / s' }, { q: 'Stroke standard deviation', op: round(sea.sigmaX, 4), sv: round(surv.sigmaX, 4), u: 'm' }, { q: 'Stroke maximum (2σ operating, 3-h MPM extreme)', op: round(sigStroke, 4), sv: round(surv.maxStroke, 4), u: 'm' },
            { q: 'Velocity (2σ operating, MPM extreme)', op: round(2 * sea.sigmaV, 4), sv: round(surv.maxVelocity, 4), u: 'm/s' }, { q: 'PTO force (2σ operating, MPM extreme)', op: round(sigForce / 1e3, 4), sv: round(surv.maxPtoForce / 1e3, 4), u: 'kN' }],
          explain: 'Displacement, velocity and PTO force of the float in the operating and extreme sea states. The extreme column is the linear upper screening value used to decide on a survival mode.' },
      ],
      notes: [
        `Float: radius ${fmt(a)} m, draft ${fmt(d)} m, mass ${fmt(mass / 1e3)} t (displaced ${fmt(disp / 1e3)} t); natural period ${fmt(Tn)} s against a design peak period of ${inp.Tp} s.`,
        `Operating sea: ${fmt(powerSea / 1e3)} kW absorbed from ${fmt(J / 1e3)} kW/m (capture-width ratio ${fmt(cw.ratio)}), PTO damping ${fmt(Bsea / 1e3)} kN·s/m.`,
        `Regular wave at Tp: ${fmt(reg.power / 1e3)} kW with ${inp.control} control${reg.constrained ? ', PTO limits active' : ''}${reg.feasible ? '' : ' — no damping satisfies the limits'}; unconstrained optimum ${fmt(reg.powerMax / 1e3)} kW.`,
        `Governing check: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `The float absorbs ${fmt(powerSea / 1e3)} kW in the operating sea (capture-width ratio ${fmt(cw.ratio)}); the factored survival stroke is ${fmt(surv.maxStroke * inp.survFactor)} m against ${inp.strokeLim} m available.`,
          checks[2].uc > 1 ? `Provide a survival mode: lock or submerge the float, or raise the survival damping — the linear estimate needs ${fmt(surv.maxStroke * inp.survFactor)} m of stroke. Then verify end-stop loads in the time-domain PTO module.` : (Tn < 0.8 * inp.Tp ? `The device is stiff for these waves (natural period ${fmt(Tn)} s, waves ${inp.Tp} s): added inertia (mass ratio about ${fmt(Math.min(3, inp.massRatio * (inp.Tp / Tn) ** 2 * 0.8))}) or reactive control would raise capture; confirm coefficients by a diffraction–radiation analysis or basin test.` : 'The response is within the PTO limits: confirm the hydrodynamic coefficients by a diffraction–radiation analysis or basin test and the viscous damping by decay tests.'),
          'IEC TS 62600-100:2012 (power matrix reporting) and IEC TS 62600-2:2019 (design); linear frequency-domain practice after Falnes (2002).',
          `Float shells: roughly ${fmt(steelT * inp.nUnits)} t of steel for ${inp.nUnits} devices. ${steelCo2Text(steelT * inp.nUnits)} Wave energy has no fuel and little visual impact; moorings and cables determine the seabed footprint.`),
        adv('economic', 'info', `The array delivers ${fmt(aepMwh, 4)} MWh/yr at a capacity factor of ${fmt(en.capacityFactor * 100)} % from a site resource of ${fmt(res.meanJ / 1e3)} kW/m.`,
          `Size the rating to the matrix: at ${inp.ratedKW} kW the capacity factor is ${fmt(en.capacityFactor * 100)} % — below about 20 % the PTO is oversized for the float, above 40 % energy is being clipped. Each point of availability is worth ${fmt(aepMwh * 0.01 / (inp.avail / 100))} MWh/yr.`,
          'LCOE conventions for ocean energy (IEA-OES 2015 International LCOE study): energy per tonne of structure and PTO rating are the cost drivers; whole-life costing per ISO 15686-5.',
          energyCo2Text(aepMwh)),
      ],
      publish: { 'wec.ratedPower': inp.ratedKW * 1e3 * inp.nUnits, 'wec.aep': en.aepWh, 'wec.BptoOpt': BptoOpt, 'wec.capacityFactor': en.capacityFactor, 'wec.maxStroke': surv.maxStroke },
      quantities: { steel_t: round(steelT * inp.nUnits, 4) },
    };
  },
  verify: [
    { name: 'Heave natural period of a truncated cylinder: T = 2π√((m + ρD³/6)/(ρgπa²))', source: 'Hand calculation with the McCormick (1981) added mass A = 0.167ρD³',
      // m = 1025·π·25·6 = 483 020 kg; A = 1025·10³/6 = 170 833 kg; K = 1025·9.80665·78.5398 = 789 461 N/m → T = 2π·√(653 853/789 461) = 5.7181 s
      inputs: { radius: 5, draft: 6, massRatio: 1, caAdd: 1 / 6 }, expect: [{ kpi: 'Tn', value: 5.7181, tol: 0.001 }] },
    { name: 'Reactive control without losses or limits: capture width = λ/2π = g/ω² in deep water', source: 'Budal & Falnes (1975); Evans (1976); Newman (1976): maximum capture width of a heaving axisymmetric body',
      // T = 8 s: ω = 0.785398 rad/s, k = ω²/g = 0.062902 rad/m → 1/k = 15.898 m (500 m depth is deep water at this period)
      inputs: { control: 'reactive', Tp: 8, depth: 500, Bvisc: 0, strokeLim: 1000, forceLim: 1e6, massRatio: 1 }, expect: [{ kpi: 'cwRegular', value: 15.898, tol: 0.002 }, { kpi: 'cwLimit', value: 15.898, tol: 0.002 }] },
  ],
};
