// Axial tension and bending stress analysis: stresses in tubes, armour wires and strength members under combined
// tension, curvature and pressure, and the tension–curvature capacity curves of the cross-section.
import { umbilicalUtilisation, umbilicalCapacityCurve, umbilicalLoadSharing } from '../../engine/xsection/umbilical.js';
import { sectionInputs, criteriaInputs, buildSection, sectionMasses, etaOf, adv, check, worst, sev, fmt, round, run, need, G, CO2_STEEL, LIMIT_RIGID, STD_UMB, CONFIRM } from './_common.js';

export default {
  id: 's12-tension-bending', suite: 12, title: 'Axial tension and bending stress',
  summary: 'von Mises utilisation of every tube and wire under combined tension, bending and pressure, with the allowable tension–curvature envelope.',
  covers: [3], fidelity: 'analytical',
  standards: [`${STD_UMB} §7.2–7.3 (stress utilisation by load condition, capacity curves)`, 'Timoshenko & Goodier, Theory of Elasticity §28 (Lamé)'],
  models: ['Axial elasticity and generalised Hooke’s law', 'Lamé thick-wall stresses with von Mises criterion', 'Helical armour stress formulations (tension share, local bending, friction stress)', 'Euler–Bernoulli bending of each element about its own axis', 'Tension–torsion coupling'],
  limits: `${LIMIT_RIGID} Elastic stresses only: plastic bending during reeling, residual stresses and weld details are not modelled. Utilisation factors are overridable defaults, not values taken from a specific edition.`,
  publishes: ['umbilical.allowTension', 'umbilical.tensionCapacity', 'umbilical.stressUtilisation'],
  searchTerms: 'umbilical capacity curve tension curvature von Mises steel tube',
  inputs: [
    ...sectionInputs(),
    { key: 'T', label: 'Effective tension', unit: 'kN', type: 'number', default: 120, min: 0, max: 5000, group: 'Loading', role: 'bc', bind: 'umbilical.topTension', uncertain: { dist: 'lognormal', cov: 0.1 },
      help: 'Maximum tension at the section (hang-off value including dynamics).' },
    { key: 'radius', label: 'Bend radius at the section', unit: 'm', type: 'number', default: 15, min: 0.5, max: 1e5, group: 'Loading', role: 'bc', bind: 'umbilical.minBendRadius', design: true, uncertain: { dist: 'lognormal', cov: 0.15 },
      help: 'Tightest radius occurring together with the tension (bend stiffener, sag bend, chute). Enter a very large value for a straight section.' },
    { key: 'zSection', label: 'Water depth at the section', unit: 'm', type: 'number', default: 10, min: 0, max: 3500, group: 'Loading', role: 'bc', help: 'External hydrostatic pressure on the tubes (bundle interstices flooded).' },
    { key: 'twist', label: 'End rotation', type: 'select', default: 'restrained', group: 'Loading', role: 'bc', options: [{ value: 'restrained', label: 'Twist restrained' }, { value: 'free', label: 'Free to twist' }] },
    ...criteriaInputs(),
    { key: 'sheathStrain', label: 'Allowable sheath bending strain', unit: '%', type: 'number', default: 7.7, min: 0.5, max: 15, group: 'Acceptance criteria', role: 'criteria', help: 'Polyethylene / polyamide 7.7 % as in API Spec 17J §5.3.1 (recalled) — confirm for the sheath compound.' },
    { key: 'lengthKm', label: 'Manufactured length', unit: 'km', type: 'number', default: 5, min: 0.05, max: 200, group: 'Quantities', role: 'operation', bind: 'cable.length' },
  ],
  compute(inp, ctx) {
    const b = buildSection(inp, ctx), s = b.sec;
    need(s.elements.some((e) => e.component.yieldStress), 'No component has a yield stress: nothing can be stress-checked. Enter the yield of the tubes or armour wires.');
    const T = inp.T * 1e3, kappa = 1 / inp.radius, eta = etaOf(inp), pExt = inp.rhoW * G * inp.zSection, restrained = inp.twist === 'restrained';
    const args = { mu: inp.muInt, pExt, twistRestrained: restrained };
    const u = run(() => umbilicalUtilisation(s, { tension: T, curvature: kappa, eta, ...args }), 'Stress utilisation');
    const ls = umbilicalLoadSharing(s, { tension: T, twistRestrained: restrained });
    const curves = {};
    for (const [name, e] of [['normal', inp.etaNormal], ['abnormal', inp.etaAbnormal], ['accidental', inp.etaAccidental]]) {
      curves[name] = run(() => umbilicalCapacityCurve(s, { eta: e, nPoints: 13, ...args }), `Capacity curve (${name}): the tubes are over-utilised by pressure alone — reduce pressure or increase wall`);
    }
    const cur = curves[inp.condition];
    // allowable tension at the operating curvature by interpolation on the capacity curve
    const interp = (c, k) => { if (k >= c.kappaLimit) return 0; const n = c.curvature.length; for (let i = 1; i < n; i++) if (k <= c.curvature[i]) { const f = (k - c.curvature[i - 1]) / (c.curvature[i] - c.curvature[i - 1]); return c.allowableTension[i - 1] + f * (c.allowableTension[i] - c.allowableTension[i - 1]); } return 0; };
    const Tallow = interp(cur, kappa);
    const els = u.elements.map((e) => {
      const el = s.elements.find((q) => q.name === e.name), c = el.component;
      return { name: e.name, kind: el.kind, vm: e.vonMises, fy: c.yieldStress, uc: e.utilisation, axial: e.axialStress, bend: e.bendingStress, fric: e.frictionStress, strainB: c.od / 2 * kappa };
    });
    const checks = els.map((e) => check(`${e.name}: von Mises`, e.vm / 1e6, eta * e.fy / 1e6, 'MPa', `${STD_UMB} §7.2 (utilisation ${eta} — ${inp.condition})`,
      e.kind === 'tube' ? 'Combined pressure, tension-share, own-axis bending and friction stress in the tube wall against the allowable fraction of yield. Exceeding it yields the tube; repeated yielding leads to low-cycle fatigue and leakage.'
        : 'Axial stress from the tension share plus local bending and friction stress in the wire against the allowable fraction of yield. Exceeding it permanently stretches the wire and unbalances the armour.'));
    const sheathEps = b.OD / 2 * kappa;
    checks.push(check('Outer sheath bending strain', sheathEps * 100, inp.sheathStrain, '%', 'API Spec 17J (2014) §5.3.1 strain limit applied to the sheath (recalled value, overridable)', 'Outer-fibre strain of the polymer sheath at the bend radius. Beyond the limit the sheath crazes or cracks and lets sea water reach the armour.'));
    const gov = worst(checks), govEl = els.reduce((a, e) => (e.uc > a.uc ? e : a), els[0]);
    const masses = sectionMasses(b), steelT = masses.steel * inp.lengthKm;
    const armourRows = els.filter((e) => e.kind === 'armour' || e.kind === 'strength');
    const armUc = armourRows.length ? Math.max(...armourRows.map((e) => e.uc)) : 0;
    const kpis = [
      { key: 'utilisation', label: 'Governing stress utilisation', value: round(u.utilisation), unit: '-', limit: 1, sense: 'max', explain: `von Mises stress of the most loaded element divided by ${eta} × yield (${inp.condition} condition). Above 1.0 the element is overstressed.` },
      { key: 'vmGov', label: 'von Mises stress, governing element', value: round(govEl.vm / 1e6), unit: 'MPa', explain: `Peak combined stress, here in "${govEl.name}". It includes pressure, tension share, bending about the element axis and friction stress.` },
      { key: 'Tallow', label: 'Allowable tension at this bend radius', value: round(Tallow / 1e3), unit: 'kN', explain: 'Tension at which the governing element reaches the allowable stress at the specified curvature and pressure — the point on the capacity curve.' },
      { key: 'Tstraight', label: 'Allowable tension, straight', value: round(cur.tensionLimit / 1e3), unit: 'kN', explain: 'Capacity-curve intercept at zero curvature for the selected load condition: the quoted maximum working tension.' },
      { key: 'mbrZeroT', label: 'Bend radius limit at zero tension', value: round(1 / cur.kappaLimit), unit: 'm', explain: 'Capacity-curve intercept on the curvature axis: the tightest radius the pressurised section accepts with no tension.' },
      { key: 'axialStrain', label: 'Axial strain of the bundle', value: round(ls.axialStrain * 100), unit: '%', explain: 'Elongation under the applied tension. Compare with the strain-free window of the fibre-optic units and the allowable conductor strain.' },
      { key: 'bendStrain', label: 'Bending strain, governing element', value: round(govEl.strainB * 100), unit: '%', explain: 'Outer-fibre strain of the governing element bent about its own axis (d/2R). Above about 0.2–0.3 % a steel tube yields in bending.' },
      { key: 'sheathStrain', label: 'Outer sheath bending strain', value: round(sheathEps * 100), unit: '%', limit: inp.sheathStrain, sense: 'max', explain: 'Outer-fibre strain of the sheath (D/2R) against the polymer strain limit.' },
      { key: 'torsion', label: restrained ? 'Reaction torque' : 'Twist', value: restrained ? round(ls.torque) : round(ls.twist * 180 / Math.PI, 3), unit: restrained ? 'N·m' : 'deg/m', explain: restrained ? 'Torque induced by the tension in a twist-restrained umbilical; it loads the termination and its anti-rotation features.' : 'Rotation per metre induced by the tension when the end is free to rotate.' },
      { key: 'armourUc', label: 'Armour-wire utilisation', value: round(armUc), unit: '-', explain: 'Highest utilisation among armour wires and strength members; a low value means the armour is there for weight and protection rather than strength.' },
    ];
    const cats = els.map((e) => e.name);
    return {
      kpis, checks,
      plots: [
        { type: 'line', title: 'Capacity curves with the operating point', x: { label: 'Curvature', unit: '1/m' }, y: { label: 'Allowable effective tension', unit: 'kN' },
          series: [...Object.entries(curves).map(([n, c]) => ({ name: `${n} (η = ${({ normal: inp.etaNormal, abnormal: inp.etaAbnormal, accidental: inp.etaAccidental })[n]})`, x: c.curvature, y: c.allowableTension.map((v) => v / 1e3), dash: n !== inp.condition })),
            { name: 'Operating point', x: [kappa], y: [inp.T] }],
          explain: 'Interaction diagram: any tension–curvature combination below the curve of the applicable load condition is acceptable. The operating point must lie inside; its distance to the curve is the reserve.' },
        { type: 'bar', title: 'Stress components by element', categories: cats, series: [{ name: 'Axial (tension share) [MPa]', y: els.map((e) => round(e.axial / 1e6)) }, { name: 'Own-axis bending [MPa]', y: els.map((e) => round(e.bend / 1e6)) }, { name: 'Friction [MPa]', y: els.map((e) => round(e.fric / 1e6)) }, { name: 'von Mises total [MPa]', y: els.map((e) => round(e.vm / 1e6)) }],
          explain: 'What each stress-checked element is made to carry. Large bending bars point to the bend radius as the lever; large axial bars to tension; the von Mises total of tubes also contains the pressure stresses.' },
        { type: 'bar', title: 'Utilisation by element', categories: cats, series: [{ name: 'Utilisation [-]', y: els.map((e) => round(e.uc)) }],
          explain: 'von Mises stress over allowable for each element. The tallest bar governs the capacity curve at this operating point; bars above 1.0 fail.' },
      ],
      tables: [{ title: 'Element stresses', explain: 'One row per stress-checked component: tension-share stress, own-axis bending stress, friction stress, the combined von Mises stress (tubes include pressure) and utilisation against the allowable fraction of yield.',
        columns: [{ key: 'name', label: 'Component' }, { key: 'axial', label: 'Axial', unit: 'MPa' }, { key: 'bend', label: 'Bending', unit: 'MPa' }, { key: 'fric', label: 'Friction', unit: 'MPa' }, { key: 'vm', label: 'von Mises', unit: 'MPa' }, { key: 'allow', label: 'Allowable', unit: 'MPa' }, { key: 'uc', label: 'Utilisation', unit: '-' }],
        rows: els.map((e) => ({ name: e.name, axial: round(e.axial / 1e6), bend: round(e.bend / 1e6), fric: round(e.fric / 1e6), vm: round(e.vm / 1e6), allow: round(eta * e.fy / 1e6), uc: round(e.uc) })) }],
      notes: [
        `Governing: ${gov.name} at utilisation ${fmt(gov.uc)} (${inp.condition} condition, η = ${eta}).`,
        `At a bend radius of ${fmt(inp.radius)} m the section accepts ${fmt(Tallow / 1e3)} kN against ${fmt(inp.T)} kN applied; straight capacity is ${fmt(cur.tensionLimit / 1e3)} kN.`,
        `Utilisation factors ${inp.etaNormal} / ${inp.etaAbnormal} / ${inp.etaAccidental} are overridable defaults. ${CONFIRM}`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}; bending contributes ${fmt(govEl.bend / 1e6)} MPa and the tension share ${fmt(govEl.axial / 1e6)} MPa of ${fmt(govEl.vm / 1e6)} MPa von Mises.`,
          gov.uc > 1 ? `Increase the bend radius to at least ${fmt(inp.radius * Math.max(gov.uc, 1.05))} m (first estimate, bending-governed) with a longer bend stiffener, or lower the tension to ${fmt(Tallow / 1e3)} kN; if pressure dominates, increase the tube wall.`
            : `Reserve at this point is ${fmt((1 - gov.uc) * 100)} %. Check every tension–curvature pair of the global analysis against the same curve, and confirm the utilisation factors against the governing edition.`,
          `${STD_UMB} §7.2–7.3 (edition to be confirmed by the purchaser)`),
        adv('economic', 'info', `Armour wires are utilised to ${fmt(armUc * 100)} % and tubes to ${fmt(Math.max(...els.filter((e) => e.kind === 'tube').map((e) => e.uc), 0) * 100)} % of their allowable stress.`,
          armUc < 0.3 && armourRows.length ? 'The armour is lightly stressed: if weight-to-diameter ratio and impact protection allow, a smaller wire or a single armour layer cuts steel tonnage and reel payload; otherwise keep it and widen the operating envelope instead.' : 'Material is well used at this operating point; savings would come from relaxing the bend radius demand rather than from the cross-section.',
          'Value engineering against quantity (AACE RP 18R-97 Class 3); whole-life costing ISO 15686-5',
          `Steel in the lay-up: ${fmt(steelT)} t over ${inp.lengthKm} km ≈ ${fmt(steelT * CO2_STEEL)} t CO₂e at the screening factor of ${CO2_STEEL} t CO₂e/t (super-duplex tube is higher). Every 10 % of armour removed saves about ${fmt(0.1 * steelT * 0.7 * CO2_STEEL)} t CO₂e.`),
      ],
      publish: { 'umbilical.allowTension': Tallow, 'umbilical.tensionCapacity': cur.tensionLimit, 'umbilical.stressUtilisation': u.utilisation },
      quantities: { steel_t: steelT, cable_km: inp.lengthKm },
    };
  },
  verify: [
    { name: 'Single unpressurised tube in tension: σ = T/A, utilisation = σ/(η·fy)',
      source: 'Axial elasticity. A = π(31.8² − 26.6²)/4 = 238.51 mm²; σ = 50 kN / 238.51 mm² = 209.63 MPa; UC = 209.63/(0.67 × 550) = 0.5689',
      inputs: { sheathT: 0, layers: [], T: 50, radius: 1e5, zSection: 0, components: [{ name: 'Tube', kind: 'tube', layer: 0, count: 1, od: 31.8, wall: 2.6, E: 200, rho: 7800, fy: 550, p: 0, rhoF: 0 }] },
      expect: [{ kpi: 'vmGov', value: 209.63, tol: 0.002 }, { kpi: 'utilisation', value: 0.5689, tol: 0.002 }] },
    { name: 'Single tube in pure bending: σ = E·(D/2)/R',
      source: 'Euler–Bernoulli beam theory. σ = 200 GPa × 15.9 mm / 20 m = 159.0 MPa',
      inputs: { sheathT: 0, layers: [], T: 0, radius: 20, zSection: 0, components: [{ name: 'Tube', kind: 'tube', layer: 0, count: 1, od: 31.8, wall: 2.6, E: 200, rho: 7800, fy: 550, p: 0, rhoF: 0 }] },
      expect: [{ kpi: 'vmGov', value: 159.0, tol: 0.002 }] },
    { name: 'Closed-end tube under internal pressure: bore von Mises = √3·p·ro²/(ro² − ri²)',
      source: 'Lamé (Timoshenko & Goodier §28). ro = 8.25, ri = 6.35 mm, p = 69 MPa: √3 × 69 × 68.0625/27.74 = 293.2 MPa',
      inputs: { sheathT: 0, layers: [], T: 0, radius: 1e5, zSection: 0, components: [{ name: 'Tube', kind: 'tube', layer: 0, count: 1, od: 16.5, wall: 1.9, E: 200, rho: 7800, fy: 550, p: 690, rhoF: 0 }] },
      expect: [{ kpi: 'vmGov', value: 293.2, tol: 0.003 }] },
  ],
};
