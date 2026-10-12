// Free vibration and dynamic amplification of the jacket: natural frequencies, mode shapes, modal mass participation
// and the dynamic amplification of the wave loading, with added mass, entrained water and marine growth.
import { G, need, run, fmt, round, adv, check, linspace, steelCo2, CONFIRM, BASIS_COST } from './_common.js';
import { geomInputs, foundationInputs, meshInputs, buildJacket, plot3d, steelMass, foundNote, testFrame } from './_model.js';
import { modalAnalysis, rayleighDamping } from '../../engine/structural/framemodal.js';
import { dynamicAmplification, harmonicResponse } from '../../engine/structural/framedynamics.js';

export default {
  id: 's01-dynamics', suite: 1, title: 'Natural frequencies and dynamic amplification',
  summary: 'Natural periods, mode shapes and effective modal masses of the platform, the dynamic amplification factor at the wave period and the frequency response of the deck.',
  covers: [9], fidelity: 'numerical',
  standards: ['ISO 19902:2020 §12.3.6 (dynamic analysis, natural period above 2.5–3 s)', 'API RP 2A-WSD 22nd ed. §5.3.3 (dynamic wave analysis)', 'DNV-RP-C205 (2019) §6.9 (added mass)', 'ISO 19901-3:2014 (machinery vibration interface)'],
  models: ['Eigenvalue analysis by subspace iteration', 'Consistent mass matrix with hydrodynamic added mass', 'Modal participation and effective modal mass', 'Single-degree-of-freedom dynamic amplification', 'Rayleigh damping', 'Harmonic response by modal superposition'],
  limits: 'Linear, undamped modes of the frame on linearised pile-head springs; the topsides is a lumped mass at the leg tops (no deck flexibility or rotary inertia). The amplification factor is the single-mode harmonic value at the wave period — irregular-sea dynamic analysis, nonlinear soil damping, diffraction and CFD are outside the module.',
  publishes: ['jacket.f1'],
  searchTerms: 'jacket natural period dynamic amplification factor modal analysis',
  inputs: [
    ...geomInputs(),
    { key: 'Wtop', label: 'Topsides weight', unit: 'MN', type: 'number', default: 40, min: 0, max: 2000, group: 'Mass', role: 'ic', bind: 'topsides.weight', uncertain: { dist: 'normal', cov: 0.05 }, help: 'Converted to mass and lumped equally on the leg tops.' },
    { key: 'Ca', label: 'Added mass coefficient of submerged members', unit: '-', type: 'number', default: 1.0, min: 0, max: 2, group: 'Mass', role: 'criteria', help: `1.0 for circular cylinders (DNV-RP-C205). ${CONFIRM}` },
    { key: 'mg', label: 'Marine growth thickness', unit: 'mm', type: 'number', default: 50, min: 0, max: 300, group: 'Mass', role: 'environment' },
    { key: 'mgDepth', label: 'Marine growth down to depth', unit: 'm', type: 'number', default: 50, min: 0, max: 400, group: 'Mass', role: 'environment' },
    { key: 'surge', label: 'Water level above still water', unit: 'm', type: 'number', default: 0, min: -5, max: 8, group: 'Mass', role: 'ic' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Mass', role: 'environment', bind: 'site.rhoW' },
    { key: 'T', label: 'Wave period of the governing sea state', unit: 's', type: 'number', default: 13.5, min: 2, max: 25, group: 'Excitation', role: 'bc', bind: 'site.Tp100' },
    { key: 'Tfat', label: 'Shortest wave period with fatigue energy', unit: 's', type: 'number', default: 4, min: 1.5, max: 12, group: 'Excitation', role: 'bc', help: 'Short operating sea states are closest to resonance and drive fatigue amplification.' },
    { key: 'zeta', label: 'Damping ratio', unit: '-', type: 'number', default: 0.02, min: 0.002, max: 0.2, group: 'Excitation', role: 'criteria', uncertain: { dist: 'lognormal', cov: 0.3 }, help: `2 % of critical for extreme waves and 1–2 % for fatigue in API RP 2A / ISO 19902 (structural + hydrodynamic + soil). ${CONFIRM}` },
    { key: 'Tlimit', label: 'Period above which dynamic analysis is required', unit: 's', type: 'number', default: 3, min: 1, max: 6, group: 'Excitation', role: 'criteria', help: `3 s in API RP 2A, 2.5 s in ISO 19902. ${CONFIRM}` },
    ...foundationInputs(),
    { key: 'nModes', label: 'Modes', type: 'int', default: 6, min: 1, max: 20, group: 'Solver', role: 'control' },
    ...meshInputs({ nSub: 1, phase: false }),
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'f1' },
  compute(inp, ctx) {
    // SI conversion: MN → kg through g; mm → m.
    const J = buildJacket(inp, ctx), mTop = inp.Wtop * 1e6 / G, n = J.topNodes.length;
    if (mTop > 0) J.model.nodeMasses = J.topNodes.map((k) => ({ node: k, m: mTop / n }));
    const opt = { nModes: inp.nModes, waterline: inp.surge, rhoW: inp.rhoW, Ca: inp.Ca };
    if (inp.mg > 0) opt.marineGrowth = { thickness: inp.mg / 1e3, rho: 1325, zTop: inp.surge, zBottom: -inp.mgDepth };
    const mo = run(() => modalAnalysis(J.model, opt), 'Modal analysis (check supports and that the model has mass)');
    need(mo.frequencies.length && mo.frequencies.every((f) => Number.isFinite(f) && f > 0), 'The modal analysis returned no positive frequency: the structure is a mechanism or has no mass.');
    const f = mo.frequencies, f1 = f[0], T1 = 1 / f1, daf = (T) => dynamicAmplification({ ratio: T1 / T, zeta: inp.zeta }).daf, dafX = daf(inp.T), dafF = daf(inp.Tfat);
    const cum = { x: 0, y: 0, z: 0 }, rowsM = f.map((fr, i) => { for (const k of ['x', 'y', 'z']) cum[k] += mo.massFraction[i][k]; return { k: i + 1, f: round(fr), T: round(1 / fr), mx: round(mo.massFraction[i].x * 100, 3), my: round(mo.massFraction[i].y * 100, 3), mz: round(mo.massFraction[i].z * 100, 3), cx: round(cum.x * 100, 3), cy: round(cum.y * 100, 3) }; });
    const ray = f.length > 1 ? rayleighDamping({ f1, f2: f[Math.min(2, f.length - 1)], zeta1: inp.zeta, zeta2: inp.zeta }) : null;
    // deck receptance: response to a 1 MN harmonic force shared by the leg tops
    const fr = linspace(0.2 * f1, Math.min(3 * f1, 1.5 * f[f.length - 1]), 90);
    const hr = run(() => harmonicResponse(mo, { frequencies: fr, forces: J.topNodes.map((k) => ({ node: k, dof: 'ux', F: 1e6 / n, amplitude: 1e6 / n })), outputs: [{ node: J.topNodes[0], dof: 'ux' }], zeta: inp.zeta }), 'Harmonic response');
    const phi = mo.modes[0], amp = Array.from({ length: J.model.nodes.length }, (_, i) => Math.hypot(phi[6 * i], phi[6 * i + 1], phi[6 * i + 2])), aMax = Math.max(...amp);
    const Ts = linspace(Math.max(1.5, 0.4 * T1), 20, 80), steel = steelMass(J, inp);
    const checks = [
      check('Natural period against the quasi-static limit', T1, inp.Tlimit, 's', 'API RP 2A-WSD 22nd ed. §5.3.3 (3 s); ISO 19902:2020 §12.3.6', 'Above this period wave energy excites the first mode noticeably and a quasi-static in-place analysis underestimates member forces; a dynamic analysis or an explicit amplification factor is then required.'),
      check('Dynamic amplification in the extreme wave', dafX, 1.1, '-', 'Single-mode DAF = 1/√((1 − r²)² + (2ζr)²); 10 % taken as the threshold for neglecting dynamics', 'Ratio of dynamic to static response at the extreme wave period. Above about 1.1 inertia loads must be added to the in-place analysis.'),
    ];
    const ok = checks.every((c) => c.uc <= 1);
    return {
      kpis: [
        { key: 'f1', label: 'First natural frequency', unit: 'Hz', value: f1, explain: 'Lowest natural frequency of the platform, normally the first sway mode. It should stay well above the frequencies where the waves carry energy.' },
        { key: 'T1', label: 'First natural period', unit: 's', value: T1, limit: inp.Tlimit, sense: 'max', explain: 'Fundamental period. Below about 3 s a jacket responds quasi-statically to storm waves; longer periods need dynamic analysis.' },
        { key: 'f2', label: 'Second natural frequency', unit: 'Hz', value: f[1] ?? f1, explain: 'Usually the sway mode in the other direction; a square jacket has two nearly equal sway frequencies.' },
        { key: 'f3', label: 'Third natural frequency', unit: 'Hz', value: f[2] ?? f1, explain: 'Typically the first torsional mode of the platform or a local brace mode.' },
        { key: 'dafExtreme', label: 'Dynamic amplification at the extreme wave period', unit: '-', value: dafX, limit: 1.1, sense: 'max', explain: 'Single-mode amplification of the response to a harmonic load at the wave period. Multiply the quasi-static base shear by this factor as a first estimate of the dynamic base shear.' },
        { key: 'dafFatigue', label: 'Dynamic amplification at the shortest fatigue wave period', unit: '-', value: dafF, explain: 'Amplification for the short operating waves that occur most often. It scales fatigue stress ranges and can dominate fatigue life when the period ratio approaches 1.' },
        { key: 'massTotal', label: 'Total vibrating mass', unit: 't', value: mo.totalMass.x / 1e3, explain: 'Structural steel, topsides, contained water, marine growth and hydrodynamic added mass moving with the structure in the x direction.' },
        { key: 'addedMass', label: 'Hydrodynamic added mass', unit: 't', value: mo.addedMass / 1e3, explain: 'Water accelerated with the submerged members. It lowers the natural frequencies and is absent in air (transport, lift).' },
        { key: 'massPart', label: 'Mass participation of the modes used (x)', unit: '%', value: cum.x * 100, explain: 'Share of the total mass captured by the computed modes in the x direction. Response-spectrum analysis needs about 90 %.' },
      ],
      checks,
      plots: [
        plot3d(J, { title: `Mode 1: ${fmt(f1)} Hz (${fmt(T1)} s)`, values: amp.map((v) => v / (aMax || 1)), label: 'Modal displacement (normalised)', unit: '-', shape: phi, explain: 'First mode shape, normalised to 1. A straight-line sway with movement at the mudline shows foundation flexibility; curvature in the legs shows frame flexibility.' }),
        { type: 'line', title: 'Dynamic amplification against wave period', x: { label: 'Wave period', unit: 's' }, y: { label: 'DAF', unit: '-' }, series: [{ name: `ζ = ${inp.zeta}`, x: Ts.map((v) => round(v)), y: Ts.map((t) => round(Math.min(daf(t), 60))) }, { name: 'Extreme wave', x: [inp.T], y: [round(dafX)] }, { name: 'Shortest fatigue wave', x: [inp.Tfat], y: [round(Math.min(dafF, 60))] }],
          explain: 'Single-mode amplification over wave period with the two operating points. The resonance peak sits at the natural period; the further to the right the design waves lie, the more quasi-static the response.' },
        { type: 'line', title: 'Deck frequency response to a 1 MN harmonic force', x: { label: 'Frequency', unit: 'Hz' }, y: { label: 'Deck displacement amplitude', unit: 'mm' }, series: [{ name: 'ux at leg top', x: fr.map((v) => round(v)), y: Array.from(hr.amplitude[0], (v) => round(v * 1e3)) }],
          explain: 'Amplitude of the deck displacement for a unit harmonic force at deck level. Peaks are the modes that the force can excite; the low-frequency plateau is the static flexibility.' },
        { type: 'bar', title: 'Effective modal mass by mode', categories: f.map((v, i) => `${i + 1}: ${fmt(v)} Hz`), series: [{ name: 'x (%)', y: rowsM.map((r) => r.mx) }, { name: 'y (%)', y: rowsM.map((r) => r.my) }, { name: 'z (%)', y: rowsM.map((r) => r.mz) }], y: { label: 'Effective mass', unit: '%' },
          explain: 'Share of the total mass that each mode mobilises in each direction. Modes with large effective mass control base shear; modes with none are local or torsional.' },
      ],
      tables: [
        { title: 'Natural frequencies and modal mass participation', columns: [{ key: 'k', label: 'Mode' }, { key: 'f', label: 'Frequency', unit: 'Hz' }, { key: 'T', label: 'Period', unit: 's' }, { key: 'mx', label: 'Mass x', unit: '%' }, { key: 'my', label: 'Mass y', unit: '%' }, { key: 'mz', label: 'Mass z', unit: '%' }, { key: 'cx', label: 'Cumulative x', unit: '%' }, { key: 'cy', label: 'Cumulative y', unit: '%' }], rows: rowsM,
          explain: 'One row per mode with its effective mass fraction per direction and the running total. Use it to decide how many modes a response-spectrum analysis needs.' },
        { title: 'Mass budget and damping', columns: [{ key: 'c', label: 'Item' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [{ c: 'Jacket steel (modelled)', v: round(J.steelT), u: 't' }, { c: 'Topsides', v: round(mTop / 1e3), u: 't' }, { c: 'Structural mass incl. contained water and growth', v: round(mo.structuralMass / 1e3), u: 't' }, { c: 'Hydrodynamic added mass', v: round(mo.addedMass / 1e3), u: 't' }, ...(ray ? [{ c: 'Rayleigh α (mass proportional)', v: round(ray.alpha), u: '1/s' }, { c: 'Rayleigh β (stiffness proportional)', v: round(ray.beta), u: 's' }] : [])],
          explain: 'Masses included in the eigenvalue problem and the Rayleigh coefficients that give the chosen damping ratio at modes 1 and 3, for use in a time-domain analysis.' },
      ],
      notes: [`First mode ${fmt(f1)} Hz (${fmt(T1)} s); period ratio to the extreme wave ${fmt(T1 / inp.T)}, DAF ${fmt(dafX)}; to the ${fmt(inp.Tfat)} s fatigue wave ${fmt(T1 / inp.Tfat)}, DAF ${fmt(dafF)}.`,
        `Eigen-solution ${mo.converged ? 'converged' : 'did not fully converge'} in ${mo.iterations} iterations; ${fmt(cum.x * 100)} % of the mass participates in x with ${f.length} modes.`, foundNote(J)],
      advice: [
        ok ? adv('engineering', 'ok', `The platform is quasi-static in the design storm: T1 = ${fmt(T1)} s (limit ${fmt(inp.Tlimit)} s), DAF ${fmt(dafX)} at ${fmt(inp.T)} s.`, dafF > 1.15 ? `Fatigue waves of ${fmt(inp.Tfat)} s are amplified by ${fmt(dafF)}: include this factor on fatigue stress ranges, and confirm the pile-head stiffness, which moves T1 most.` : 'No dynamic correction is needed for strength. Re-check the period whenever topsides weight, water depth assumptions or pile-head stiffness change.', 'API RP 2A-WSD 22nd ed. §5.3.3; ISO 19902:2020 §12.3.6', `Stiffness comes mostly from leg spacing and foundation, not from steel mass: widening the base lowers the period without adding much steel. The jacket holds ${steelCo2(steel)}`)
          : adv('engineering', 'warn', `Dynamic effects are significant: T1 = ${fmt(T1)} s against ${fmt(inp.Tlimit)} s, DAF ${fmt(dafX)} at the extreme wave period.`, `Apply the DAF to the wave load in the in-place analysis or run a time-domain dynamic analysis; to shorten the period, stiffen the foundation, widen the base or reduce topsides mass (T ∝ √(mass/stiffness): 20 % less mass gives ${fmt((1 - Math.sqrt(0.8)) * 100)} % shorter period).`, 'API RP 2A-WSD 22nd ed. §5.3.3; ISO 19902:2020 §12.3.6', `Inertia loads grow steel demand roughly in proportion to the DAF. The jacket holds ${steelCo2(steel)}`),
        adv('economic', 'info', `Topsides mass is ${fmt(mTop / (mo.totalMass.x || 1) * 100)} % of the vibrating mass; each 10 % of topsides weight moves the period by about ${fmt(0.5 * mTop / (mo.totalMass.x || 1) * 10)} %.`, 'Weight growth during topsides design is the usual reason a quasi-static jacket becomes dynamically sensitive: hold a weight budget with contingency, because a late dynamic re-analysis and stiffening costs more than the allowance.', `${BASIS_COST}; ISO 19901-5:2021 (weight control)`),
      ],
      publish: { 'jacket.f1': f1 },
      quantities: { steel_t: round(steel) },
    };
  },
  verify: [
    // Uniform cantilever in air: f1 = (1.8751²/2π)·√(EI/(m L⁴)); EI = 1.14418e10 N·m², m = 7850·0.20570 = 1614.7 kg/m, L = 95 m:
    // f1 = 0.55959·√(1.14418e10/(1614.7·95⁴)) = 0.55959·0.29496 = 0.16506 Hz.
    { name: 'Uniform cantilever: f1 = (1.8751²/2π)√(EI/mL⁴)', source: 'Blevins, Formulas for Natural Frequency and Mode Shape (1979) Table 8-1',
      inputs: { geom: testFrame('column'), Wtop: 0, Ca: 0, mg: 0, found: 'fixed', nSub: 8, nModes: 2 },
      expect: [{ kpi: 'f1', value: 0.16506, tol: 0.005 }] },
  ],
};
