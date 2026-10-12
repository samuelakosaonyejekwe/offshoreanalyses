// Natural-frequency and resonance analysis: fore-aft modes of the monopile–tower–RNA system on foundation springs,
// the 1P / 3P excitation bands of the variable-speed rotor, blade flapwise and edgewise frequencies with centrifugal
// stiffening (Campbell diagram), the drivetrain torsional mode and the sensitivity to scour and soil stiffness.
import { buildSupportStructure } from '../../engine/aero/coupled.js';
import { genericBladeDistributions, campbellData } from '../../engine/aero/blade.js';
import { turbineInputs, getTurbine, supportInputs, supportSegments, foundationInputs, foundationOf, segmentsMass, site, adv, check, worst, sev, fmt, round, run, need, linspace, RPM,
  LIMIT_POLAR, LIMIT_CODE, CONFIRM, steelCo2Text } from './_common.js';

export default {
  id: 's16-natural-frequency', suite: 16, title: 'Natural frequencies and resonance margins',
  summary: 'First natural frequencies of the support structure, blades and drivetrain against the 1P and 3P rotor excitation bands and the wave frequencies, with a Campbell diagram and sensitivity to scour and foundation stiffness.',
  covers: [11], fidelity: 'numerical',
  standards: ['DNV-ST-0126 (2021) §4 and §7 (support structures: natural frequency and resonance)', 'IEC 61400-1:2019 §7.5 (resonance and stability)', 'IEC 61400-3-1:2019', 'DNV-ST-0376 (2015) Rotor blades', 'Arany, Bhattacharya, Macdonald & Hogan (2017) Soil Dyn. Earthq. Eng. 92 (monopile design steps)'],
  models: ['Euler–Bernoulli beam finite elements with consistent mass', 'Modal analysis (generalised symmetric eigenproblem)', 'Coupled rotor–tower–foundation dynamics: lumped RNA and three-spring foundation', 'Hydrodynamic added mass and enclosed water', 'Rotating-beam centrifugal stiffening (Southwell coefficients)', 'Campbell diagram', 'Two-inertia torsional drivetrain model'],
  limits: `Fore-aft bending in one plane: side-side modes (within a few percent for an axisymmetric structure), torsion, P-delta softening from the RNA weight and soil damping are not modelled; blade distributions are generic shapes scaled to the stated mass and standstill frequencies, without flap–edge–torsion coupling. ${LIMIT_POLAR} ${LIMIT_CODE}`,
  publishes: ['owt.f1', 'owt.f2'],
  searchTerms: 'offshore wind turbine monopile natural frequency 1P 3P soft-stiff Campbell diagram',
  inputs: [
    ...turbineInputs(),
    ...supportInputs({ steel: false }),
    ...foundationInputs(),
    { key: 'scour', label: 'Scour depth at the pile', unit: 'm', type: 'number', default: 0, min: 0, max: 20, group: 'Foundation (mudline boundary)', role: 'ic', bind: 'scour.depth', help: 'Initial seabed state: local scour lengthens the free pile below the original mudline (the springs act at the scoured level).' },
    { key: 'growth', label: 'Marine-growth thickness on the submerged pile', unit: 'mm', type: 'number', default: 50, min: 0, max: 300, group: 'Support structure', role: 'ic', help: 'Adds mass (1325 kg/m³) over the wet length; 50–100 mm is usual in the North Sea.' },
    { key: 'flapHz', label: 'Blade first flapwise frequency at standstill', unit: 'Hz', type: 'number', default: 0.70, min: 0.2, max: 3, group: 'Blade and drivetrain', role: 'material', uncertain: { dist: 'normal', cov: 0.05 }, help: 'NREL 5-MW reference blade: about 0.70 Hz. Scale with 1/length for a scaled rotor (about 0.40 Hz at 15 MW).' },
    { key: 'edgeHz', label: 'Blade first edgewise frequency at standstill', unit: 'Hz', type: 'number', default: 1.08, min: 0.3, max: 5, group: 'Blade and drivetrain', role: 'material', uncertain: { dist: 'normal', cov: 0.05 }, help: 'NREL 5-MW reference blade: about 1.08 Hz.' },
    { key: 'kShaft', label: 'Drivetrain torsional stiffness (low-speed shaft)', unit: 'MN·m/rad', type: 'number', default: 867.6, min: 10, max: 1e5, group: 'Blade and drivetrain', role: 'material', help: 'NREL 5-MW reference: 867.6 MN·m/rad equivalent on the low-speed shaft.' },
    site('Tp1', { group: 'Excitation', label: 'Peak period of the frequent sea states', default: 7, help: 'Wave energy near the first natural frequency drives fatigue; typical operational seas have Tp of 5–9 s.', bind: undefined }),
    { key: 'margin', label: 'Required separation from 1P and 3P', unit: '%', type: 'number', default: 10, min: 0, max: 25, group: 'Acceptance', role: 'criteria', help: `DNV-ST-0126: natural frequency outside ±10 % of the rotor and blade-passing frequency ranges (±5 % with justification). ${CONFIRM}` },
    { key: 'bladeMargin', label: 'Required separation of blade modes from 3P and 6P at rated speed', unit: '%', type: 'number', default: 5, min: 0, max: 25, group: 'Acceptance', role: 'criteria' },
    { key: 'nElem', label: 'Beam elements', type: 'int', default: 40, min: 6, max: 160, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Beam elements', levels: [10, 20, 40, 80], qoi: 'f1' },
  compute(inp) {
    const tb = getTurbine(inp), rotor = tb.rotor, R = rotor.tipRadius;
    need(inp.hubHeight > R + 5, `The hub height (${inp.hubHeight} m) leaves no clearance under the blade tip (rotor radius ${fmt(R)} m).`);
    const fnd = foundationOf(inp);
    const build = (over = {}, f = fnd, nElem = inp.nElem) => {
      const q = { ...inp, ...over };
      const segs = supportSegments({ ...q, depth: q.depth + q.scour }, { zHub: q.hubHeight, hubOffset: 2.4 * tb.scale });
      // marine growth as extra mass on the wet monopile: equivalent density increase of the shell
      if (q.growth > 0) { const tg = q.growth / 1e3, s0 = segs[0]; s0.density += 1325 * Math.PI * (s0.dBase + tg) * tg / (Math.PI * (s0.dBase - s0.tBase) * s0.tBase) * Math.min(1, q.depth / (q.depth + q.scour + q.transZ)); }
      return { segs, str: run(() => buildSupportStructure({ segments: segs, rnaMass: rotor.rnaMass, zHub: q.hubHeight, foundation: f, nElem, nModes: 3, zeta: 0.01 }), 'Support-structure model') };
    };
    const { segs, str } = build();
    const [f1, f2, f3] = str.modes.map((m) => m.frequency);
    const p1 = [rotor.minRotorSpeed / (2 * Math.PI), rotor.ratedRotorSpeed / (2 * Math.PI)], p3 = p1.map((v) => v * rotor.blades);
    const mg = inp.margin / 100;
    // position of f1: soft-soft (< 1P), soft-stiff (between), stiff-stiff (> 3P)
    const region = f1 < p1[0] ? 'soft-soft (below 1P)' : f1 < p1[1] ? 'inside the 1P band' : f1 < p3[0] ? 'soft-stiff (between 1P and 3P)' : f1 < p3[1] ? 'inside the 3P band' : 'stiff-stiff (above 3P)';
    const sep1 = f1 > p1[1] ? f1 / p1[1] - 1 : f1 < p1[0] ? 1 - f1 / p1[0] : 0;
    const sep3 = f1 < p3[0] ? 1 - f1 / p3[0] : f1 > p3[1] ? f1 / p3[1] - 1 : 0;
    // blade Campbell data
    const gb = run(() => genericBladeDistributions({ length: rotor.bladeLength, mass: rotor.bladeMass, flapFrequency: inp.flapHz, edgeFrequency: inp.edgeHz }), 'Blade distributions');
    const oms = linspace(0, 1.25 * rotor.ratedRotorSpeed, 11);
    const cb = run(() => campbellData({ x: gb.x, EIflap: gb.EIflap, EIedge: gb.EIedge, massPerLength: gb.massPerLength, omegas: oms, hubRadius: rotor.hubRadius, nModes: 2, harmonics: [1, 3, 6, 9] }), 'Campbell diagram');
    const at = (arr) => { const x = rotor.ratedRotorSpeed; let i = 1; while (oms[i] < x) i++; return arr[i - 1] + (arr[i] - arr[i - 1]) * (x - oms[i - 1]) / (oms[i] - oms[i - 1]); };
    const flapR = at(cb.flap[0]), edgeR = at(cb.edge[0]);
    const bladeSep = Math.min(...[flapR, edgeR].flatMap((f) => [3, 6].map((n) => Math.abs(f - n * p1[1]) / (n * p1[1]))));
    const inRange = cb.crossings.filter((c) => c.omega >= rotor.minRotorSpeed && c.omega <= rotor.ratedRotorSpeed);
    // drivetrain: free-free torsional frequency of rotor and generator inertia on the low-speed shaft
    const Jg = rotor.generatorInertia * rotor.gearboxRatio ** 2, Jr = rotor.rotorInertia;
    const fDrive = Math.sqrt(inp.kShaft * 1e6 * (1 / Jr + 1 / Jg)) / (2 * Math.PI);
    const checks = [
      check('Separation of the first mode from the 1P band', mg, Math.max(sep1, 1e-6), '-', 'DNV-ST-0126 (2021) §4.5 (resonance, ±10 % of 1P)', 'Demand is the required relative separation, capacity the actual distance of the first natural frequency from the rotor-speed range. Inside the margin, rotor imbalance excites the tower at resonance for long periods: fatigue and nacelle vibration.'),
      check('Separation of the first mode from the 3P band', mg, Math.max(sep3, 1e-6), '-', 'DNV-ST-0126 (2021) §4.5 (resonance, ±10 % of 3P)', 'Required against actual separation from the blade-passing frequency range. Blade passage (tower shadow, shear) loads the tower three times per revolution; resonance with it must be passed through quickly or avoided by a speed exclusion zone.'),
      check('Blade modes clear of 3P and 6P at rated speed', inp.bladeMargin / 100, Math.max(bladeSep, 1e-6), '-', 'DNV-ST-0376 (2015) §3 (blade resonance); IEC 61400-1:2019 §7.5', 'Required against actual separation of the first flapwise and edgewise blade frequencies from the 3P and 6P harmonics at rated speed. Coincidence amplifies blade loads and can couple with tower modes.'),
    ];
    const gov = worst(checks);
    // sensitivities
    const scs = linspace(0, Math.max(1.5 * inp.monoD, inp.scour), 7), fScour = scs.map((s) => build({ scour: s }, fnd, 20).str.modes[0].frequency);
    const mult = [0.25, 0.5, 1, 2, 4], fStiff = fnd === 'fixed' ? mult.map(() => f1) : mult.map((m) => build({}, { Kh: fnd.Kh * m, Kr: fnd.Kr * m, Khr: fnd.Khr * m }, 20).str.modes[0].frequency);
    const fFixed = build({}, 'fixed', 20).str.modes[0].frequency;
    const steel = segmentsMass(supportSegments(inp, { zHub: inp.hubHeight, hubOffset: 2.4 * tb.scale })) / 1e3;
    // stiffness needed to restore the margin: f ∝ √EI ∝ D^1.5 at constant D/t
    const fTarget = p1[1] * (1 + mg), dScale = f1 < fTarget ? (fTarget / f1) ** (2 / 3) : 1;
    const rpm = oms.map((o) => round(o / RPM, 4));
    const z = Array.from(str.nodes, (v) => round(v, 3));
    return {
      kpis: [
        { key: 'f1', label: 'First fore-aft natural frequency', value: round(f1, 5), unit: 'Hz', explain: `First bending mode of the whole system. It is ${region}; most monopile turbines are designed soft-stiff, between the 1P and 3P bands.` },
        { key: 'f2', label: 'Second fore-aft natural frequency', value: round(f2, 5), unit: 'Hz', explain: 'Second bending mode (tower bending with the RNA nearly stationary). It can be excited by 6P and by breaking-wave impact.' },
        { key: 'p1max', label: '1P at rated rotor speed', value: round(p1[1], 5), unit: 'Hz', explain: 'Rotor rotational frequency at rated speed: upper end of the 1P band.' },
        { key: 'p3min', label: '3P at minimum rotor speed', value: round(p3[0], 5), unit: 'Hz', explain: 'Blade-passing frequency at cut-in speed: lower end of the 3P band.' },
        { key: 'sep1', label: 'Separation from the 1P band', value: round(sep1 * 100, 4), unit: '%', limit: inp.margin, sense: 'min', explain: 'Distance of the first frequency from the nearest edge of the rotor-speed range. 0 means it lies inside the band.' },
        { key: 'sep3', label: 'Separation from the 3P band', value: round(sep3 * 100, 4), unit: '%', limit: inp.margin, sense: 'min', explain: 'Distance of the first frequency from the nearest edge of the blade-passing range.' },
        { key: 'fFixed', label: 'First frequency with a clamped mudline', value: round(fFixed, 5), unit: 'Hz', explain: 'Upper bound with an infinitely stiff foundation. The drop to the spring-supported value is the effect of soil flexibility (typically 5–15 %).' },
        { key: 'waveRatio', label: 'First frequency × wave peak period', value: round(f1 * inp.Tp1, 4), unit: '-', explain: 'Ratio of the natural frequency to the wave peak frequency. Close to 1 means direct wave resonance; values of 1.5–2.5 are typical and still give noticeable dynamic amplification.' },
        { key: 'flapRated', label: 'Blade first flapwise frequency at rated speed', value: round(flapR, 5), unit: 'Hz', explain: 'Standstill flap frequency raised by centrifugal stiffening. Compare with 3P and 6P on the Campbell diagram.' },
        { key: 'edgeRated', label: 'Blade first edgewise frequency at rated speed', value: round(edgeR, 5), unit: 'Hz', explain: 'Edgewise frequency at rated speed; lightly damped, so resonance with a harmonic is more damaging than for flap.' },
        { key: 'fDrive', label: 'Drivetrain torsional frequency (free–free)', value: round(fDrive, 5), unit: 'Hz', explain: 'Rotor inertia against generator inertia on the shaft stiffness. It must stay clear of 3P and 6P and is damped actively through generator torque.' },
        { key: 'modalMass', label: 'Modal mass of the first mode (unit hub displacement)', value: round(str.modes[0].mass / 1e3, 5), unit: 't', explain: 'Effective mass of the first mode referred to the hub: RNA plus the participating tower, pile, added water and marine growth.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Frequency diagram: tower modes against 1P and 3P', x: { label: 'Rotor speed', unit: 'rpm' }, y: { label: 'Frequency', unit: 'Hz' },
          series: [{ name: '1P', x: rpm, y: oms.map((o) => round(o / (2 * Math.PI), 5)) }, { name: '3P', x: rpm, y: oms.map((o) => round(rotor.blades * o / (2 * Math.PI), 5)) }, { name: 'Tower mode 1', x: [rpm[0], rpm[rpm.length - 1]], y: [f1, f1] }, { name: 'Tower mode 2', x: [rpm[0], rpm[rpm.length - 1]], y: [f2, f2] },
            { name: 'Minimum rotor speed', x: [rotor.minRotorSpeed / RPM, rotor.minRotorSpeed / RPM], y: [0, f2], dash: true }, { name: 'Rated rotor speed', x: [rotor.ratedRotorSpeed / RPM, rotor.ratedRotorSpeed / RPM], y: [0, f2], dash: true }],
          explain: 'Excitation lines 1P and 3P against rotor speed with the tower frequencies as horizontal lines. A crossing between the two dashed speed limits is a resonance in the operating range.' },
        { type: 'line', title: 'Blade Campbell diagram', x: { label: 'Rotor speed', unit: 'rpm' }, y: { label: 'Frequency', unit: 'Hz' },
          series: [{ name: 'Flap 1', x: rpm, y: cb.flap[0].map((v) => round(v, 5)) }, { name: 'Edge 1', x: rpm, y: cb.edge[0].map((v) => round(v, 5)) }, { name: 'Flap 2', x: rpm, y: cb.flap[1].map((v) => round(v, 5)) }, ...cb.excitation.map((e) => ({ name: `${e.n}P`, x: rpm, y: e.f.map((v) => round(v, 5)), dash: true }))],
          explain: 'Blade frequencies rising with rotor speed (centrifugal stiffening) and the nP excitation rays. Intersections inside the operating speed range are potential blade resonances.' },
        { type: 'line', title: 'Fore-aft mode shapes', x: { label: 'Lateral displacement (unit at hub)', unit: '-' }, y: { label: 'Elevation above still water', unit: 'm' }, series: str.modes.map((m, i) => ({ name: `Mode ${i + 1}: ${fmt(m.frequency)} Hz`, x: Array.from(m.shape, (v) => round(v, 4)), y: z })),
          explain: 'Bending mode shapes from the mudline to the tower top. Displacement and slope at the mudline show how much the foundation participates.' },
        { type: 'line', title: 'First frequency versus scour depth', x: { label: 'Scour depth', unit: 'm' }, y: { label: 'Frequency', unit: 'Hz' }, series: [{ name: 'f1', x: scs.map((v) => round(v, 3)), y: fScour.map((v) => round(v, 5)) }, { name: '1P at rated + margin', x: [scs[0], scs[scs.length - 1]], y: [fTarget, fTarget], dash: true }],
          explain: 'Loss of natural frequency as scour deepens the free length of the pile. Where the curve meets the dashed line, scour protection or monitoring becomes a requirement.' },
        { type: 'line', title: 'First frequency versus foundation stiffness', x: { label: 'Multiplier on K_LL, K_LR and K_RR', unit: '-' }, y: { label: 'Frequency', unit: 'Hz' }, series: [{ name: 'f1', x: mult, y: fStiff.map((v) => round(v, 5)) }, { name: 'Clamped', x: [mult[0], mult[mult.length - 1]], y: [fFixed, fFixed], dash: true }],
          explain: 'Sensitivity of the first frequency to the soil stiffness, which is the most uncertain input. A flat curve means the design is robust to soil variability; a steep one calls for site-specific cyclic stiffness.' },
      ],
      tables: [
        { title: 'Natural frequencies and excitation bands', columns: [{ key: 'item', label: 'Item' }, { key: 'lo', label: 'From', unit: 'Hz' }, { key: 'hi', label: 'To', unit: 'Hz' }, { key: 'note', label: 'Remark' }],
          rows: [{ item: '1P (rotor speed)', lo: round(p1[0], 4), hi: round(p1[1], 4), note: `${fmt(rotor.minRotorSpeed / RPM)}–${fmt(rotor.ratedRotorSpeed / RPM)} rpm` }, { item: `${rotor.blades}P (blade passing)`, lo: round(p3[0], 4), hi: round(p3[1], 4), note: '' },
            { item: 'Tower mode 1', lo: round(f1, 4), hi: round(f1, 4), note: region }, { item: 'Tower mode 2', lo: round(f2, 4), hi: round(f2, 4), note: '' }, { item: 'Tower mode 3', lo: round(f3, 4), hi: round(f3, 4), note: '' },
            { item: 'Blade flap 1', lo: round(cb.flap[0][0], 4), hi: round(flapR, 4), note: 'standstill → rated' }, { item: 'Blade edge 1', lo: round(cb.edge[0][0], 4), hi: round(edgeR, 4), note: 'standstill → rated' }, { item: 'Drivetrain torsion', lo: round(fDrive, 4), hi: round(fDrive, 4), note: 'free–free' }, { item: 'Wave peak', lo: round(1 / inp.Tp1, 4), hi: round(1 / inp.Tp1, 4), note: `Tp = ${inp.Tp1} s` }],
          explain: 'All natural frequencies and excitation ranges side by side. Overlaps between a structural frequency and an excitation band are the resonances to design out or to manage by control.' },
        { title: 'Blade resonance crossings', columns: [{ key: 'mode', label: 'Blade mode' }, { key: 'n', label: 'Harmonic' }, { key: 'rpm', label: 'Rotor speed', unit: 'rpm' }, { key: 'f', label: 'Frequency', unit: 'Hz' }, { key: 'op', label: 'In operating range' }],
          rows: cb.crossings.map((c) => ({ mode: c.mode, n: `${c.harmonic}P`, rpm: round(c.omega / RPM, 4), f: round(c.f, 4), op: c.omega >= rotor.minRotorSpeed && c.omega <= rotor.ratedRotorSpeed ? 'yes' : 'no' })),
          explain: 'Rotor speeds at which a blade frequency equals a harmonic of the rotor speed. Crossings marked "yes" lie between minimum and rated speed and need damping or a speed exclusion window.' },
      ],
      notes: [
        `${rotor.name}: f1 = ${fmt(f1, 4)} Hz is ${region}; 1P band ${fmt(p1[0])}–${fmt(p1[1])} Hz, 3P band ${fmt(p3[0])}–${fmt(p3[1])} Hz.`,
        `Foundation flexibility lowers f1 by ${fmt((1 - f1 / fFixed) * 100)} % from the clamped value of ${fmt(fFixed, 4)} Hz; scour of ${fmt(scs[scs.length - 1])} m would bring it to ${fmt(fScour[fScour.length - 1], 4)} Hz.`,
        `Governing check: ${gov.name} (utilisation ${fmt(gov.uc)}). ${inRange.length ? `${inRange.length} blade crossing(s) inside the operating speed range.` : 'No blade crossing with 1P–9P inside the operating speed range.'}`,
        'Blade stiffness and mass distributions are generic shapes scaled to the stated mass and standstill frequencies.',
      ],
      advice: [
        adv('engineering', sev(gov.uc), `The first natural frequency is ${fmt(f1, 4)} Hz (${region}), ${fmt(sep1 * 100)} % from the 1P band and ${fmt(sep3 * 100)} % from the 3P band; required separation ${inp.margin} %.`,
          sep1 < mg && f1 >= p1[0] ? `Raise the first frequency to at least ${fmt(fTarget, 4)} Hz: scale the monopile and tower diameters by about ${fmt((dScale - 1) * 100)} % at constant D/t, or stiffen the foundation (longer or larger pile); otherwise define a rotor-speed exclusion zone in the controller.` : sep3 < mg ? `The frequency is too close to the 3P band: reduce stiffness or mass placement so that f1 ≤ ${fmt(p3[0] * (1 - mg), 4)} Hz, or apply a speed exclusion window around ${fmt(f1 / rotor.blades / RPM * 2 * Math.PI, 3)} rpm.` : `Margins are met. Confirm the cyclic foundation stiffness (the frequency moves from ${fmt(fStiff[0], 4)} to ${fmt(fStiff[fStiff.length - 1], 4)} Hz over a factor 16 in stiffness) and plan frequency monitoring against scour.`,
          'DNV-ST-0126 (2021) §4.5 and §7 (frequency window, ±10 % from 1P and 3P); design steps after Arany et al. (2017).',
          `The frequency window, not strength, often sets monopile diameter: the support structure here is ${fmt(steel)} t of steel. ${steelCo2Text(steel)}`),
        adv('economic', 'info', `Stiffness is bought with steel: at constant D/t the first frequency scales with about D^1.5, so +10 % in frequency costs roughly +${fmt((1.1 ** (2 / 3) * 1.1 ** (2 / 3) - 1) * 100)} % in steel mass (${fmt(steel * (1.1 ** (4 / 3) - 1))} t here).`,
          'Before adding steel, check cheaper measures: a speed exclusion zone (a small energy loss), scour protection to keep the frequency, or measured soil stiffness to remove conservatism. Carry frequency monitoring into the O&M plan.',
          'Whole-life costing per ISO 15686-5: structural mass, energy lost in exclusion zones and monitoring cost are compared over the design life; steel tonnage is the main CAPEX quantity of the foundation.'),
      ],
      publish: { 'owt.f1': f1, 'owt.f2': f2 },
      quantities: { steel_t: round(steel, 4) },
    };
  },
  verify: [
    { name: 'Uniform cantilever with tip mass: f₁ = (1/2π)·√(3EI / ((M + 33/140·mL)·L³))', source: 'Rayleigh quotient with the static deflection shape; Blevins, Formulas for Natural Frequency and Mode Shape, Table 8-8',
      // D = 6 m, t = 50 mm: I = 4.1363 m⁴, m = 7336.8 kg/m, L = 5 + 87.6 = 92.6 m; k = 3EI/L³ = 3.2819 MN/m; M_eff = 350 000 + 0.2357·679 388 = 510 130 kg → 0.4037 Hz
      inputs: { foundation: 'fixed', depth: 5, monoD: 6, monoT: 50, towerBaseD: 6, towerTopD: 6, towerBaseT: 50, towerTopT: 50, towerMassFactor: 1, hubHeight: 90, transZ: 10, scour: 0, growth: 0 },
      expect: [{ kpi: 'f1', value: 0.4037, tol: 0.02 }] },
    { name: 'Drivetrain free–free torsion: f = (1/2π)·√(k·(1/J_r + 1/(N²·J_g)))', source: 'Two-inertia torsional oscillator with NREL 5-MW data (Jonkman et al. 2009, Table 5-1): hand calculation',
      // J_r = 38 759 228 kg·m² (hub + 3 blades), J_g·N² = 534.116·97² = 5 025 497 kg·m², k = 867.637 MN·m/rad
      // ω² = 867.637e6·(2.58003e-8 + 1.98985e-7) = 195.03 → ω = 13.965 rad/s → 2.2226 Hz; 1P at 12.1 rpm = 0.20167 Hz
      inputs: { kShaft: 867.637 }, expect: [{ kpi: 'fDrive', value: 2.2226, tol: 0.003 }, { kpi: 'p1max', value: 0.201667, tol: 1e-4 }] },
  ],
};
