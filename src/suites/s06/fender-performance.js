// Fender system performance: reaction, deflection, absorbed energy, correction factors, hull pressure and the
// force and moment handed to the dolphin (PIANC WG33 2002 §4–6, BS 6349-4:2014).
import { GENERIC_FENDERS, fenderCurve, fenderAtEnergy } from './_berthing.js';
import { need, fmt, round, adv, check, CONFIRM, NOT_HIFI, ECON_BASIS } from './_common.js';

export default {
  id: 's06-fender-performance', suite: 6, title: 'Fender system performance',
  summary: 'Deflection, reaction and absorbed energy of the fender for normal and abnormal berthing, with velocity, temperature, angle and tolerance corrections, hull pressure and the loads passed to the dolphin.',
  covers: [10], fidelity: 'analytical',
  standards: ['PIANC WG33 Guidelines for the Design of Fender Systems (2002) §4.3, §6 and Appendix A', 'PIANC WG211 (2024)', 'BS 6349-4:2014 §6', 'ASTM F2192 (fender performance testing)'],
  models: ['Nonlinear fender force–deflection curve', 'Energy as the area under the force–deflection curve', 'Rate (velocity) and temperature correction factors', 'Angular compression factor', 'Hysteresis energy dissipation', 'Hull contact pressure', 'Coulomb friction on the fender face'],
  limits: 'The built-in curves are generic normalised shapes for each fender family, not manufacturer data: enter the certified curve for design. Rubber hyperelasticity, viscoelastic rate dependence and hysteresis are represented through the measured curve, the correction factors and a dissipation fraction, not through a constitutive model. Units in contact share the impact equally. ' + NOT_HIFI,
  publishes: ['dolphin.energyCapacity', 'dolphin.fenderReaction', 'dolphin.fenderMoment', 'dolphin.fenderDeflection'],
  searchTerms: 'rubber fender performance curve reaction energy correction factor hull pressure',
  inputs: [
    { key: 'EN', label: 'Normal berthing energy', unit: 'kJ', type: 'number', default: 410, min: 0.1, max: 1e5, group: 'Demand', role: 'bc', bind: 'dolphin.berthEnergyNormal', uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'EA', label: 'Abnormal berthing energy', unit: 'kJ', type: 'number', default: 515, min: 0.1, max: 1e5, group: 'Demand', role: 'bc', bind: 'dolphin.berthEnergyAbnormal' },
    { key: 'delta0', label: 'Fender compression before impact', unit: 'mm', type: 'number', default: 0, min: 0, max: 500, group: 'Demand', role: 'ic', help: 'Zero for a vessel approaching a free fender; pre-compression applies to a vessel already leaning on the fender.' },
    { key: 'type', label: 'Fender type', type: 'select', default: 'cone', group: 'Fender', role: 'geometry', options: Object.entries(GENERIC_FENDERS).map(([value, g]) => ({ value, label: `${g.label} — generic curve` })) },
    { key: 'H', label: 'Fender height (or diameter)', unit: 'm', type: 'number', default: 1.3, min: 0.1, max: 4.5, group: 'Fender', role: 'geometry', design: true },
    { key: 'RR', label: 'Rated reaction per unit', unit: 'kN', type: 'number', default: 1300, min: 5, max: 20000, group: 'Fender', role: 'geometry', design: true, help: 'Catalogue reaction at rated deflection, reference conditions (constant velocity, 23 °C, 0°). Used with the generic curves only.' },
    { key: 'ERcat', label: 'Rated energy per unit from the catalogue (0 = area under the curve)', unit: 'kJ', type: 'number', default: 0, min: 0, max: 1e5, group: 'Fender', role: 'geometry', help: 'When given, the energy curve is scaled to this value at rated deflection; the reaction curve is unchanged.' },
    { key: 'nUnits', label: 'Units sharing the impact', unit: '-', type: 'int', default: 1, min: 1, max: 12, group: 'Fender', role: 'geometry', design: true },
    { key: 'curve', label: 'Certified manufacturer curve (replaces the generic shape)', type: 'table', group: 'Fender', role: 'geometry',
      columns: [{ key: 'defl', label: 'Deflection', unit: '% of height', type: 'number' }, { key: 'reaction', label: 'Reaction', unit: 'kN', type: 'number' }], default: [],
      help: 'Reaction against deflection of one unit from the type-approval test report. Leave empty to use the generic curve.' },
    { key: 'ratedDefl', label: 'Rated deflection of the certified curve (0 = last row)', unit: '%', type: 'number', default: 0, min: 0, max: 90, group: 'Fender', role: 'geometry' },
    { key: 'VF', label: 'Velocity factor', unit: '-', type: 'number', default: 1.0, min: 0.8, max: 1.5, group: 'Correction factors', role: 'criteria', help: `Ratio of performance at the real compression rate to the test rate; depends on compound and strain rate (manufacturer data). ${CONFIRM}` },
    { key: 'TFlow', label: 'Temperature factor at the lowest temperature (reaction)', unit: '-', type: 'number', default: 1.1, min: 0.9, max: 2.5, group: 'Correction factors', role: 'criteria', help: 'Rubber stiffens when cold: governs the largest reaction.' },
    { key: 'TFhigh', label: 'Temperature factor at the highest temperature (energy)', unit: '-', type: 'number', default: 0.95, min: 0.5, max: 1.1, group: 'Correction factors', role: 'criteria', help: 'Rubber softens when hot: governs the smallest energy.' },
    { key: 'AF', label: 'Angle factor on energy', unit: '-', type: 'number', default: 0.95, min: 0.4, max: 1.05, group: 'Correction factors', role: 'criteria', help: 'Loss of energy capacity under angular compression (berthing angle, flare, bow radius); from the manufacturer table for the compression angle.' },
    { key: 'tol', label: 'Manufacturing tolerance', unit: '%', type: 'number', default: 10, min: 0, max: 20, group: 'Correction factors', role: 'criteria', help: '±10 % on reaction and energy is the usual performance tolerance.' },
    { key: 'hyst', label: 'Hysteresis: fraction of absorbed energy dissipated', unit: '-', type: 'number', default: 0.5, min: 0, max: 0.95, group: 'Correction factors', role: 'criteria', help: 'About 0.5–0.6 for buckling rubber fenders, 0.1–0.3 for pneumatic and foam units; from the unloading branch of the test curve.' },
    { key: 'mu', label: 'Friction coefficient of the fender face', unit: '-', type: 'number', default: 0.2, min: 0, max: 1, group: 'Contact', role: 'bc', help: '0.2 for UHMW-PE facing on steel, 0.5 or more for bare rubber.' },
    { key: 'panelW', label: 'Contact width on the hull', unit: 'm', type: 'number', default: 3.0, min: 0.2, max: 15, group: 'Contact', role: 'geometry', design: true },
    { key: 'panelH', label: 'Contact height on the hull', unit: 'm', type: 'number', default: 4.0, min: 0.2, max: 15, group: 'Contact', role: 'geometry', design: true },
    { key: 'pHull', label: 'Permissible hull pressure', unit: 'kPa', type: 'number', default: 200, min: 20, max: 1500, group: 'Contact', role: 'criteria', help: `Indicative: 150–200 kPa for large tankers and bulk carriers, 200–400 kPa for container ships, below 200 kPa for gas carriers. ${CONFIRM}` },
    { key: 'eArm', label: 'Fender centreline above the pile-head level', unit: 'm', type: 'number', default: 1.5, min: -10, max: 15, group: 'Contact', role: 'bc', help: 'Lever arm that turns the reaction into a moment on the pile cap.' },
    { key: 'Rallow', label: 'Design lateral load of the dolphin', unit: 'kN', type: 'number', default: 2500, min: 1, max: 1e5, group: 'Acceptance', role: 'criteria', bind: 'dolphin.lateralDesignLoad', help: 'Reaction the structure behind the fender has been designed for.' },
  ],
  compute(inp) {
    const n = inp.nUnits, tol = inp.tol / 100;
    // SI: kN → N, kJ → J, mm → m
    const c = fenderCurve({ type: inp.type, H: inp.H, RR: inp.RR * 1e3, table: inp.curve, ratedDefl: inp.ratedDefl });
    const eScale = inp.ERcat > 0 ? inp.ERcat * 1e3 / c.ERated : 1;
    if (eScale !== 1) { for (let i = 0; i < c.E.length; i++) c.E[i] *= eScale; c.ERated *= eScale; }
    need(c.ERated > 0, 'The fender curve encloses no energy: check the reaction values.');
    const fMin = (1 - tol) * inp.TFhigh * inp.AF * inp.VF, fMax = (1 + tol) * inp.TFlow * inp.VF;
    // Energy already stored by the pre-compression is not available for the impact.
    const d0 = inp.delta0 / 1e3;
    need(d0 < c.dRated, `The pre-compression (${fmt(inp.delta0)} mm) reaches the rated deflection (${fmt(c.dRated * 1e3)} mm): no energy capacity is left.`);
    let i0 = 0; while (i0 < c.d.length - 1 && c.d[i0 + 1] <= d0) i0++;
    const E0 = c.E[i0] + (c.E[Math.min(i0 + 1, c.E.length - 1)] - c.E[i0]) * (d0 - c.d[i0]) / Math.max(c.d[Math.min(i0 + 1, c.d.length - 1)] - c.d[i0], 1e-12);
    const state = (EkJ, f) => fenderAtEnergy(c, EkJ * 1e3 / n + E0 * f, { fE: f, fR: f });
    const nSoft = state(inp.EN, fMin), nHard = state(inp.EN, fMax), aSoft = state(inp.EA, fMin), aHard = state(inp.EA, fMax);
    const Ecap = n * (c.ERated - E0) * fMin / 1e3;            // kJ, at rated deflection, minimum performance
    const Rn = n * Math.max(nSoft.peakReaction, nHard.peakReaction) / 1e3, Ra = n * Math.max(aSoft.peakReaction, aHard.peakReaction) / 1e3;
    const area = n * inp.panelW * inp.panelH;
    const pN = Rn / area, pA = Ra / area;
    const eff = c.ERated / (c.RRated * c.H);
    const Erebound = inp.EN * (1 - inp.hyst);
    const bottomed = aSoft.exceeded;
    const checks = [
      check('Abnormal energy against capacity at rated deflection (minimum performance)', inp.EA, Ecap, 'kJ', 'PIANC WG33 (2002) §4.3 and App. A; BS 6349-4:2014 §6.5', 'Guards against the fender being compressed beyond its rated deflection. Beyond it the reaction climbs steeply, the rubber is overstrained and the surplus energy goes into hull and piles.'),
      check('Normal-berthing deflection against rated deflection', nSoft.deflection * 1e3, c.dRated * 1e3, 'mm', 'PIANC WG33 (2002) §4.3', 'Keeps routine berthings within the working range of the fender, with the softest credible unit (hot, angled, minimum tolerance). Repeated over-compression shortens the life of the rubber.'),
      check('Hull pressure in abnormal berthing', pA, inp.pHull, 'kPa', 'PIANC WG33 (2002) Table 4.4.1', 'Reaction over the contact area against the pressure the hull plating between frames can take. If exceeded the shell plating may be set in.'),
      check('Reaction on the dolphin (maximum performance)', Ra, inp.Rallow, 'kN', 'BS 6349-4:2014 §6.6', 'The stiffest credible fender (cold, maximum tolerance) delivers the largest reaction for the same energy; it must not exceed the lateral load the dolphin is designed for.'),
    ];
    const ok = checks.every((q) => q.uc <= 1);
    const w = checks.reduce((a, q) => (q.uc > a.uc ? q : a));
    // size hint: scale rated reaction and height together (energy ∝ RR·H)
    const need1 = inp.EA / Ecap;
    const rubberT = (h, units) => units * 0.9 * h ** 3 * 1.2; // rough rubber volume ≈ 0.9·H³ per buckling unit at 1.2 t/m³ (screening)
    const dPct = c.d.map((d) => round(d / c.H * 100, 4));
    const advice = [
      adv('engineering', ok ? 'ok' : 'fail', ok ? `The fender absorbs the abnormal energy ${fmt(inp.EA)} kJ within its capacity ${fmt(Ecap)} kJ; design reaction ${fmt(Ra)} kN, hull pressure ${fmt(pA)} kPa.` : `${w.name}: ${fmt(w.demand)} against ${fmt(w.capacity)} ${w.unit}.`,
        ok ? `Energy reserve ${fmt((1 - need1) * 100)} %. Confirm the correction factors with the manufacturer's type-approval data for the selected compound and repeat for the lowest and highest water levels.`
          : (inp.EA > Ecap ? `Increase the rated energy by a factor of ${fmt(need1)} (larger unit, harder compound or ${Math.ceil(n * need1)} units), keeping the reaction below ${fmt(inp.Rallow)} kN.` : pA > inp.pHull ? `Enlarge the panel to at least ${fmt(Ra / inp.pHull)} m² of contact area.` : `Select a softer or taller fender: the same energy at a lower reaction, or raise the design lateral load of the dolphin to ${fmt(Ra)} kN.`),
        'PIANC WG33 (2002) §4.3–4.4 and Appendix A; PIANC WG211 (2024); BS 6349-4:2014',
        `Rubber mass grows roughly with the cube of fender height: one ${fmt(inp.H)} m unit is of the order of ${fmt(rubberT(inp.H, 1))} t of polymer; over-sizing by 10 % in height adds about 33 % of that.`),
      adv('economic', ok ? 'info' : 'warn', `The dolphin sees ${fmt(Ra)} kN and ${fmt(Ra * inp.eArm)} kN·m per abnormal berthing; the energy-to-reaction efficiency of this fender is ${fmt(eff)} (E/(R·H)).`,
        'Compare fender options on reaction per unit energy, not on unit price: a lower reaction for the same energy reduces pile diameter and wall, which usually outweighs the fender cost. Allow for replacement of rubber units once in a 25–30 year berth life.',
        `${ECON_BASIS}; PIANC WG33 §7 whole-life considerations`),
    ];
    const Rcat = c.R.map((r) => round(r * n / 1e3));
    return {
      kpis: [
        { key: 'Ecap', label: 'Energy capacity (minimum performance)', unit: 'kJ', value: Ecap, limit: inp.EA, sense: 'min', explain: 'Energy absorbed at rated deflection by all units with every adverse factor applied (hot, angled, minimum tolerance). It must cover the abnormal berthing energy.' },
        { key: 'ERated', label: 'Rated energy per unit', unit: 'kJ', value: c.ERated / 1e3, explain: 'Area under the catalogue force–deflection curve up to rated deflection, at reference conditions. The basis to which the correction factors are applied.' },
        { key: 'deflN', label: 'Deflection in normal berthing', unit: 'mm', value: nSoft.deflection * 1e3, limit: c.dRated * 1e3, sense: 'max', explain: 'Compression of the softest credible fender under the normal berthing energy. It should stay below the rated deflection.' },
        { key: 'deflA', label: 'Deflection in abnormal berthing', unit: 'mm', value: aSoft.deflection * 1e3, explain: 'Compression under the abnormal energy. At or beyond rated deflection the fender has no further useful capacity.' },
        { key: 'deflNpct', label: 'Normal deflection / height', unit: '%', value: nSoft.deflection / c.H * 100, explain: 'Normal-berthing compression as a percentage of the fender height, for comparison with the catalogue rated deflection.' },
        { key: 'RN', label: 'Peak reaction in normal berthing', unit: 'kN', value: Rn, explain: 'Largest reaction reached on the way to the normal-berthing deflection, for the stiffest credible fender. A buckling fender reaches its plateau early, so this is close to the rated reaction even in light berthings.' },
        { key: 'RA', label: 'Design reaction (abnormal berthing)', unit: 'kN', value: Ra, limit: inp.Rallow, sense: 'max', explain: 'Largest reaction handed to the dolphin: abnormal energy on the stiffest credible fender. This is the lateral design load of the structure.' },
        { key: 'MA', label: 'Moment at pile-head level', unit: 'kN·m', value: Math.abs(Ra * inp.eArm), explain: 'Design reaction times the height of the fender centreline above the pile heads: the moment transmitted to the pile cap.' },
        { key: 'shear', label: 'Friction force along the berth', unit: 'kN', value: inp.mu * Ra, explain: 'Shear dragged along the fender face by vessel movement, μ·R. It loads the dolphin sideways and in torsion.' },
        { key: 'pHull', label: 'Hull pressure', unit: 'kPa', value: pA, limit: inp.pHull, sense: 'max', explain: 'Design reaction over the panel contact area. Above the permissible pressure the hull plating between frames may yield.' },
        { key: 'eff', label: 'Energy efficiency E/(R·H)', unit: '-', value: eff, explain: 'Rated energy over rated reaction times height. 0.5 or more marks an efficient buckling fender; a linear spring is 0.5 of its own deflection, about 0.25 on height.' },
        { key: 'Ediss', label: 'Energy dissipated per normal berthing', unit: 'kJ', value: inp.EN * inp.hyst, explain: 'Part of the absorbed energy turned to heat in the rubber (hysteresis). The remainder is returned to the vessel as rebound.' },
        { key: 'Ereb', label: 'Rebound energy', unit: 'kJ', value: Erebound, explain: 'Energy given back to the vessel on unloading. High rebound pushes the vessel off the berth and loads the mooring lines.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Fender reaction against deflection', x: { label: 'Deflection', unit: '% of height' }, y: { label: 'Reaction (all units)', unit: 'kN' },
          series: [{ name: c.generic ? 'Catalogue (generic shape)' : 'Catalogue (certified curve)', x: dPct, y: Rcat }, { name: 'Maximum performance', x: dPct, y: Rcat.map((r) => round(r * fMax)), dash: true }, { name: 'Minimum performance', x: dPct, y: Rcat.map((r) => round(r * fMin)), dash: true },
            { name: 'Normal berthing', x: [round(nSoft.deflection / c.H * 100)], y: [round(n * nSoft.reaction / 1e3)] }, { name: 'Abnormal berthing', x: [round(aSoft.deflection / c.H * 100)], y: [round(n * aSoft.reaction / 1e3)] }],
          explain: 'Force–deflection curve with the band between softest and stiffest credible performance. The operating points should sit before the final steep rise; the plateau height is the reaction the dolphin must carry.' },
        { type: 'line', title: 'Absorbed energy against deflection', x: { label: 'Deflection', unit: '% of height' }, y: { label: 'Energy (all units)', unit: 'kJ' },
          series: [{ name: 'Catalogue', x: dPct, y: c.E.map((e) => round(e * n / 1e3)) }, { name: 'Minimum performance', x: dPct, y: c.E.map((e) => round(e * n * fMin / 1e3)), dash: true }, { name: 'Normal energy', x: [0, dPct[dPct.length - 1]], y: [inp.EN, inp.EN], dash: true }, { name: 'Abnormal energy', x: [0, dPct[dPct.length - 1]], y: [inp.EA, inp.EA], dash: true }],
          explain: 'Area under the reaction curve. Where the demand lines cross the minimum-performance curve is the working deflection; a crossing beyond the rated deflection means the fender is too small.' },
        { type: 'line', title: 'Loading and unloading in a normal berthing (hysteresis)', x: { label: 'Deflection', unit: '% of height' }, y: { label: 'Reaction (all units)', unit: 'kN' },
          series: (() => { const k = dPct.findIndex((d) => d > nSoft.deflection / c.H * 100); const m = k < 0 ? dPct.length : Math.max(k, 2); return [{ name: 'Loading', x: dPct.slice(0, m), y: Rcat.slice(0, m).map((r) => round(r * fMin)) }, { name: 'Unloading', x: dPct.slice(0, m), y: Rcat.slice(0, m).map((r) => round(r * fMin * (1 - inp.hyst))) }]; })(),
          explain: 'The enclosed area between loading and unloading is the energy dissipated as heat in one berthing; the area under the unloading branch returns to the vessel.' },
      ],
      tables: [{ title: 'Fender states', columns: [{ key: 'case', label: 'Case' }, { key: 'd', label: 'Deflection', unit: 'mm' }, { key: 'pct', label: 'Deflection', unit: '%' }, { key: 'R', label: 'Peak reaction', unit: 'kN' }, { key: 'over', label: 'Beyond curve' }],
        rows: [['Normal, minimum performance', nSoft], ['Normal, maximum performance', nHard], ['Abnormal, minimum performance', aSoft], ['Abnormal, maximum performance', aHard]].map(([name, s]) => ({ case: name, d: round(s.deflection * 1e3), pct: round(s.deflection / c.H * 100), R: round(n * s.peakReaction / 1e3), over: s.exceeded ? 'yes' : 'no' })),
        explain: 'Each row is one combination of energy and fender performance. Minimum performance governs deflection and energy; maximum performance governs reaction.' }],
      notes: [c.generic ? `Generic ${GENERIC_FENDERS[inp.type].label} curve used — replace with the certified manufacturer curve for design.` : 'Certified manufacturer curve used.',
        `Combined factors: energy × ${fmt(fMin)} (minimum), reaction × ${fmt(fMax)} (maximum).`, ...(bottomed ? ['The abnormal energy exceeds the end of the curve: the deflection is reported at the last tabulated point and the true reaction is higher.'] : []),
        ...(eScale !== 1 ? [`Energy curve scaled by ${fmt(eScale)} to the catalogue rated energy.`] : [])],
      advice,
      publish: { 'dolphin.energyCapacity': Ecap * 1e3, 'dolphin.fenderReaction': Ra * 1e3, 'dolphin.fenderMoment': Math.abs(Ra * inp.eArm) * 1e3, 'dolphin.fenderDeflection': aSoft.deflection },
      quantities: { polymer_t: round(rubberT(inp.H, n), 3) },
    };
  },
  verify: [
    // k = 1000 kN / 0.5 m = 2000 kN/m; E = ½ k δ² → δ = √(2·100/2000) = 0.31623 m; R = k δ = 632.46 kN; rated E = ½·1000·0.5 = 250 kJ.
    { name: 'Linear fender: δ = √(2E/k), R = kδ', source: 'Closed form for a linear spring, prescribed force–deflection curve 0–1000 kN over 0–50 % of 1 m',
      inputs: { curve: [{ defl: 0, reaction: 0 }, { defl: 50, reaction: 1000 }], H: 1, nUnits: 1, EN: 100, EA: 100, VF: 1, TFlow: 1, TFhigh: 1, AF: 1, tol: 0, ERcat: 0, delta0: 0 },
      expect: [{ kpi: 'deflN', value: 316.23, tol: 0.002 }, { kpi: 'RN', value: 632.46, tol: 0.002 }, { kpi: 'ERated', value: 250, tol: 1e-6 }, { kpi: 'eff', value: 0.25, tol: 1e-6 }] },
    // Constant reaction 500 kN reached at 0.5 % and held to 60 % of 2 m: E = 500·(1.2 − 0.01) + ½·500·0.01 = 597.5 kJ.
    { name: 'Constant-reaction fender: E = R·δ', source: 'Area under a prescribed rectangular force–deflection test curve',
      inputs: { curve: [{ defl: 0, reaction: 0 }, { defl: 0.5, reaction: 500 }, { defl: 60, reaction: 500 }], H: 2, nUnits: 1, EN: 100, EA: 100, VF: 1, TFlow: 1, TFhigh: 1, AF: 1, tol: 0 },
      expect: [{ kpi: 'ERated', value: 597.5, tol: 1e-4 }, { kpi: 'RN', value: 500, tol: 1e-6 }] },
  ],
};
