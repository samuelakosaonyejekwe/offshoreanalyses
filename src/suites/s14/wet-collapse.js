// Carcass collapse and wet collapse: with a breached outer sheath the sea pressure reaches the pressure sheath
// and the carcass alone resists it. Compares dry and wet capacity, the flooding scenario and its side effects.
import { collapseAssessment, ovalisedRingCollapse, ringElasticCollapse } from '../../engine/xsection/flexcollapse.js';
import { annulusFreeVolume, annulusFlooding } from '../../engine/xsection/annulus.js';
import { stackInputs, stackFrom, kpi, chk, adv, plot, bar, tbl, col, run, fmt, round, linspace, sev, CONFIRM, UNBONDED, GENERIC, STD_17J, STD_17B, CO2_MGO, MGO_T_PER_DAY, BAR, G } from './_common.js';
import { ringInputs, ringsFrom } from './_collapse.js';

export default {
  id: 's14-wet-collapse',
  suite: 14,
  title: 'Carcass collapse and wet collapse',
  summary: 'Checks the carcass alone against sea pressure acting on the pressure sheath after the outer sheath is breached, and quantifies the flooded-annulus scenario.',
  covers: [7],
  fidelity: 'analytical',
  standards: ['API RP 17B (5th ed.) §8.4 (collapse with flooded annulus)', 'API Spec 17J (4th ed.) Table 6'],
  models: ['Carcass collapse as an oval ring at first yield (Timoshenko)', 'Wet-collapse model: carcass alone, pressure on the pressure sheath', 'Interlocked-profile equivalent ring stiffness (API RP 17B)', 'Annulus flooding volume balance'],
  limits: 'Ring model at first yield with an equivalent ring for the interlocked carcass: no profile contact mechanics, no sheath–carcass gap or sheath-collapse mode, no plasticity or post-buckling. Wet-collapse capacity must be confirmed by flooded-annulus collapse tests on the actual carcass. ' + UNBONDED + ' ' + GENERIC,
  publishes: ['flexible.wetCollapseDepth'],
  searchTerms: 'flexible riser wet collapse carcass flooded annulus',
  inputs: [
    ...stackInputs(),
    { key: 'depth', label: 'Water depth at the damage', unit: 'm', type: 'number', default: 1000, min: 10, max: 3500, group: 'External pressure', role: 'bc', bind: 'site.depth', design: true, help: 'Depth at which sea water enters the annulus and acts on the pressure sheath.' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'External pressure', role: 'environment', bind: 'site.rhoW' },
    { key: 'pMin', label: 'Minimum bore pressure', unit: 'bar', type: 'number', default: 0, min: 0, max: 500, group: 'External pressure', role: 'bc', help: 'Bore pressure that can be relied on. 0 = depressurised or empty.' },
    { key: 'annulusState', label: 'Annulus at the start of the event', type: 'select', default: 'dry', options: [{ value: 'dry', label: 'Initially dry, floods through a breach' }, { value: 'flooded', label: 'Already flooded' }], group: 'External pressure', role: 'ic',
      help: 'Both cases end with the carcass carrying the full pressure; an initially dry annulus additionally takes on water mass.' },
    { key: 'floodLength', label: 'Flooded length', unit: 'm', type: 'number', default: 500, min: 1, max: 20000, group: 'External pressure', role: 'ic', help: 'Length of annulus that fills with sea water.' },
    ...ringInputs(),
  ],
  compute(inp, ctx) {
    const S = stackFrom(inp, ctx);
    const R = ringsFrom(inp, S);
    const pDesign = Math.max(inp.rhoW * G * inp.depth - inp.pMin * BAR, 0);
    const base = { carcass: R.carcass, pressureArmour: R.pa ?? undefined, curvatureFactor: inp.curvFactor, factor: R.factor, rho: inp.rhoW, ovality: inp.ovality / 100 };
    const wet = run(() => collapseAssessment({ ...base, annulus: 'flooded', pDesign }), 'Wet-collapse assessment');
    const dry = collapseAssessment({ ...base, annulus: 'dry', pDesign });
    const pCr = ringElasticCollapse({ EI: R.carcass.EI, R: R.carcass.R }).pCr;
    const ring = ovalisedRingCollapse({ pCr, yieldStress: R.carcass.yieldStress, R: R.carcass.R, thickness: R.carcass.thickness, sectionModulus: R.carcass.sectionModulus, ovality: inp.ovality / 100 });
    // annulus free volume between the two seals
    const between = S.stack.layers.filter((L) => L.index > S.stack.sealInner && L.index < S.stack.sealOuter);
    const vol = between.length ? annulusFreeVolume({ layers: between.map((L) => ({ rIn: L.rIn, rOut: L.rOut, fill: L.kind === 'helical' ? Math.min(1, L.nA / Math.cos(L.alpha) / (Math.PI * (L.rOut ** 2 - L.rIn ** 2))) : 1 })), length: inp.floodLength }) : { volume: 0, voidArea: 0 };
    const flood = vol.volume > 0 ? annulusFlooding({ freeVolume: vol.volume, floodedLength: inp.floodLength, breachDepth: inp.depth, rho: inp.rhoW }) : { waterVolume: 0, addedMassPerLength: 0, annulusPressure: 101325, consequences: [] };
    const ovs = linspace(0, Math.max(2, 2 * inp.ovality), 21);
    const pOv = ovs.map((o) => collapseAssessment({ ...base, ovality: o / 100, annulus: 'flooded' }).pCollapse);
    const ds = linspace(0, Math.max(1.3 * dry.collapseDepth, 1.2 * inp.depth), 25);
    const pOf = (d) => Math.max(inp.rhoW * G * d - inp.pMin * BAR, 0);
    const uc = wet.utilisation;
    const tReq = R.tEq * Math.cbrt(Math.max(uc, 1e-6));
    const pBoreReq = Math.max(0, pDesign - R.factor * wet.pCollapse);
    return {
      kpis: [
        kpi('pWet', 'Wet-collapse pressure', wet.pCollapse / BAR, 'bar', 'Pressure on the pressure sheath at which the carcass alone collapses (first yield of the oval ring, with the curvature knock-down).'),
        kpi('pDry', 'Dry-collapse pressure', dry.pCollapse / BAR, 'bar', 'Collapse pressure of the intact pipe, with the pressure armour helping. The difference to the wet value is what a breach costs.'),
        kpi('pElasticCarcass', 'Carcass elastic critical pressure', pCr / BAR, 'bar', 'Ring buckling pressure 3EI/R³ of the carcass alone for a perfectly round ring.'),
        kpi('pYieldCarcass', 'Carcass membrane yield pressure', ring.pYield / BAR, 'bar', 'Pressure at which the carcass hoop compression alone reaches yield; an upper bound for thick, stocky carcasses.'),
        kpi('wetDepth', 'Wet-collapse water depth', wet.collapseDepth, 'm', 'Depth at which a flooded annulus would collapse the empty pipe.'),
        kpi('wetAllowDepth', 'Allowable depth, flooded annulus', wet.allowableDepth, 'm', 'Wet-collapse depth multiplied by the permissible utilisation.'),
        kpi('ucWet', 'Wet-collapse utilisation', uc, '-', 'Net external pressure over the permissible fraction of the wet-collapse pressure. Above 1.0 the pipe cannot tolerate a breached outer sheath at this depth when depressurised.', 1, 'max'),
        kpi('ucDry', 'Dry-collapse utilisation', dry.utilisation, '-', 'Same ratio for the intact pipe, for comparison.'),
        kpi('wetLoss', 'Capacity lost on flooding', (1 - wet.pCollapse / dry.pCollapse) * 100, '%', 'Reduction of collapse pressure when the annulus floods and the pressure armour no longer supports the carcass.'),
        kpi('floodWater', 'Sea water entering the annulus', flood.waterVolume * 1e3, 'L', 'Volume of water that fills the annulus voids over the flooded length.'),
        kpi('floodMass', 'Added mass from flooding', flood.addedMassPerLength, 'kg/m', 'Extra mass in air per metre; it matters for recovery loads and for the dynamic response of the riser.'),
      ],
      checks: [
        chk('Wet collapse: carcass alone', pDesign / BAR, R.factor * wet.pCollapse / BAR, 'bar', `${STD_17B} §8.4; ${STD_17J} Table 6 (factor ${fmt(R.factor)})`,
          'Guards against collapse of the carcass when sea water behind a damaged outer sheath presses directly on the pressure sheath. Collapse blocks the bore and usually means replacing the pipe.'),
        chk('Dry collapse: carcass + pressure armour', pDesign / BAR, R.factor * dry.pCollapse / BAR, 'bar', `${STD_17B} §8.4`,
          'The same failure mode for the intact pipe, where the pressure armour shares the ovalisation load with the carcass.'),
      ],
      plots: [
        plot('line', 'Wet-collapse pressure vs ovality', 'Ovality', '%', 'Collapse pressure', 'bar', [
          { name: 'Carcass alone', x: ovs.map((o) => round(o)), y: pOv.map((p) => round(p / BAR)) },
          { name: 'Design external pressure', x: [0, round(ovs[ovs.length - 1])], y: [round(pDesign / BAR), round(pDesign / BAR)], dash: true },
        ], 'Imperfection sensitivity of the carcass. A carcass dented or ovalised during installation moves to the right on this curve; where it meets the dashed line the margin is gone.'),
        plot('line', 'Utilisation vs depth, dry and flooded annulus', 'Water depth', 'm', 'Utilisation', '-', [
          { name: 'Flooded annulus', x: ds.map((d) => round(d)), y: ds.map((d) => round(pOf(d) / (R.factor * wet.pCollapse))) },
          { name: 'Dry annulus', x: ds.map((d) => round(d)), y: ds.map((d) => round(pOf(d) / (R.factor * dry.pCollapse))) },
          { name: 'Limit', x: [0, round(ds[ds.length - 1])], y: [1, 1], dash: true },
        ], 'The two lines cross the limit at the allowable depths with and without a flooded annulus. The gap between them is the depth range in which outer-sheath integrity is safety critical.'),
        bar('Collapse pressure by condition', ['Dry annulus', 'Flooded annulus', 'Design pressure'], 'Pressure', 'bar', [{ name: 'Pressure', y: [round(dry.pCollapse / BAR), round(wet.pCollapse / BAR), round(pDesign / BAR)] }],
          'Collapse resistance intact and flooded against the pressure to be resisted. Both bars, multiplied by the permissible factor, must exceed the third.'),
      ],
      tables: [tbl('Flooding scenario', [col('item', 'Item'), col('value', 'Value'), col('unit', 'Unit')], [
        { item: 'Annulus free volume over the flooded length', value: round(vol.volume * 1e3), unit: 'L' },
        { item: 'Annulus pressure at the breach (absolute)', value: round(flood.annulusPressure / BAR), unit: 'bar' },
        { item: 'Added mass in air', value: round(flood.addedMassPerLength), unit: 'kg/m' },
        { item: 'Bore pressure needed to pass the wet-collapse check', value: round(pBoreReq / BAR), unit: 'bar' },
        ...flood.consequences.map((t, i) => ({ item: `Consequence ${i + 1}`, value: t, unit: '' })),
      ], 'What a breach of the outer sheath changes: water volume, annulus pressure, weight, and the design checks that must be re-run for a flooded annulus.')],
      notes: [`Wet collapse at ${fmt(wet.pCollapse / BAR)} bar against ${fmt(dry.pCollapse / BAR)} bar dry: flooding removes ${fmt((1 - wet.pCollapse / dry.pCollapse) * 100)} % of the capacity.`,
        `Annulus ${inp.annulusState === 'dry' ? 'initially dry: flooding adds ' + fmt(flood.addedMassPerLength) + ' kg/m' : 'already flooded'}; utilisation ${fmt(uc)} with factor ${fmt(R.factor)}.`, GENERIC],
      advice: [
        adv('engineering', sev(uc), `With a flooded annulus the carcass alone is ${fmt(uc * 100)} % utilised at ${fmt(inp.depth)} m; allowable depth ${fmt(wet.allowableDepth)} m.`,
          uc > 1 ? `Either hold a minimum bore pressure of ${fmt(pBoreReq / BAR)} bar whenever the annulus may be flooded, or stiffen the carcass to an equivalent thickness of about ${fmt(tReq * 1e3)} mm.` : `The pipe tolerates a breached outer sheath when depressurised. Still treat annulus flooding as an integrity event: it changes the fatigue environment of the armour.`,
          `${STD_17B} §8.4; ${STD_17J} Table 6. ${CONFIRM}`),
        adv('economic', uc > 1 ? 'warn' : 'info', `Flooding costs ${fmt((1 - wet.pCollapse / dry.pCollapse) * 100)} % of the collapse capacity, so outer-sheath integrity carries real value here.`,
          `Budget periodic annulus vacuum or volume testing and a sheath-repair clamp: a test takes a few hours from the platform, while a collapsed riser means a replacement campaign of the order of 10–20 vessel days plus deferred production.`,
          'Risk-based inspection per DNV-RP-G101 / API RP 580; consequence in vessel days and downtime',
          `A replacement campaign of 15 vessel days burns about ${fmt(15 * MGO_T_PER_DAY)} t of marine gas oil ≈ ${fmt(15 * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening factors: ${MGO_T_PER_DAY} t/day, ${CO2_MGO} t CO₂ per t fuel), besides scrapping the pipe — monitoring the annulus avoids both.`),
      ],
      publish: { 'flexible.wetCollapseDepth': wet.collapseDepth },
    };
  },
  verify: [{
    // Round ring (ovality → 0): collapse = min(elastic buckling, membrane yield). Carcass E = 193 GPa, t_eq = 6 mm,
    // R = 79.7 mm: p_cr = 3·(193e9·0.006³/12)/0.0797³ = 20.586 MPa. Membrane yield p_y = σ_y·a/R with a = fill × t =
    // 0.55·7 mm = 3.85 mm and σ_y = 500 MPa: 500e6·0.00385/0.0797 = 24.153 MPa. The smaller root is p_cr = 205.86 bar.
    name: 'Perfectly round carcass: collapse equals the Levy ring-buckling pressure',
    source: 'Levy (1884); Timoshenko & Gere, Theory of Elastic Stability, Ch. 7; hand calculation',
    inputs: { ringMethod: 'equivalent', tEq: 6, ovality: 0, curvFactor: 1 },
    expect: [{ kpi: 'pWet', value: 205.86, tol: 0.002 }, { kpi: 'pYieldCarcass', value: 241.53, tol: 0.002 }],
  }],
};
