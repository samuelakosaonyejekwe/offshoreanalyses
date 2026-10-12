// External-pressure collapse of the intact pipe (dry annulus): the carcass and the pressure armour ovalise
// together and their ring stiffnesses add; first-yield collapse with initial ovality.
import { collapseAssessment, combinedElasticCollapse } from '../../engine/xsection/flexcollapse.js';
import { stackInputs, stackFrom, kpi, chk, adv, plot, tbl, col, run, fmt, round, linspace, sev, CONFIRM, UNBONDED, GENERIC, STD_17J, STD_17B, CO2_STEEL, BAR, G } from './_common.js';
import { ringInputs, ringsFrom } from './_collapse.js';

export default {
  id: 's14-collapse',
  suite: 14,
  title: 'External pressure collapse (intact pipe)',
  summary: 'Calculates the hydrostatic collapse pressure and depth of the pipe with a dry annulus, where carcass and pressure armour resist together.',
  covers: [6],
  fidelity: 'analytical',
  standards: ['API RP 17B (5th ed.) §8.4 (collapse)', 'API Spec 17J (4th ed.) Table 6 (carcass buckling factor)'],
  models: ['Elastic ring buckling p = 3EI/R³ (Levy; Timoshenko & Gere)', 'Summed ring stiffness of carcass and pressure armour', 'Imperfection-sensitive first-yield collapse of an oval ring (Timoshenko)', 'Equivalent ring stiffness of an interlocked profile'],
  limits: 'Ring model at first yield: the interlocked carcass and pressure armour are equivalent rings, with no profile contact mechanics, plasticity or post-buckling. It does not replace the manufacturer\'s validated collapse methodology and hydrostatic collapse tests. The bending knock-down is a user factor. ' + UNBONDED + ' ' + GENERIC,
  publishes: ['flexible.collapseDepth'],
  searchTerms: 'flexible pipe carcass collapse external pressure ovality',
  inputs: [
    ...stackInputs(),
    { key: 'depth', label: 'Maximum water depth', unit: 'm', type: 'number', default: 1000, min: 10, max: 3500, group: 'External pressure', role: 'bc', bind: 'site.depth', design: true, help: 'Sets the external hydrostatic pressure on the outer sheath.' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'External pressure', role: 'environment', bind: 'site.rhoW' },
    { key: 'pMin', label: 'Minimum bore pressure', unit: 'bar', type: 'number', default: 0, min: 0, max: 500, group: 'External pressure', role: 'bc', help: 'Bore pressure that can be relied on to oppose the external pressure. 0 = empty pipe (installation or depressurised), the conservative case.' },
    ...ringInputs(),
  ],
  compute(inp, ctx) {
    const S = stackFrom(inp, ctx);
    const R = ringsFrom(inp, S);
    const pDesign = Math.max(inp.rhoW * G * inp.depth - inp.pMin * BAR, 0);
    const base = { carcass: R.carcass, pressureArmour: R.pa ?? undefined, annulus: 'dry', curvatureFactor: inp.curvFactor, factor: R.factor, rho: inp.rhoW };
    const c = run(() => collapseAssessment({ ...base, ovality: inp.ovality / 100, pDesign }), 'Collapse assessment');
    const rings = [{ EI: R.carcass.EI, R: R.carcass.R }, ...(R.pa ? [R.pa] : [])];
    const el = combinedElasticCollapse({ rings });
    const ovs = linspace(0, Math.max(2, 2 * inp.ovality), 21);
    const pOv = ovs.map((o) => collapseAssessment({ ...base, ovality: o / 100 }).pCollapse);
    const ds = linspace(0, Math.max(1.3 * c.collapseDepth, inp.depth * 1.2), 25);
    const ucD = ds.map((d) => Math.max(inp.rhoW * G * d - inp.pMin * BAR, 0) / (R.factor * c.pCollapse));
    const uc = c.utilisation;
    // ring stiffness scales with t_eq³: stiffness needed to bring the utilisation to 1 (elastic-governed estimate)
    const tReq = R.tEq * Math.cbrt(Math.max(uc, 1e-6));
    const carSteel = R.C.rho * R.C.nA / Math.cos(R.C.alpha); // kg/m = t/km
    return {
      kpis: [
        kpi('pCollapse', 'Collapse pressure (dry annulus)', c.pCollapse / BAR, 'bar', 'External pressure at which the oval carcass ring first yields with the pressure armour sharing the load. Collapse of the bore follows.'),
        kpi('pElastic', 'Elastic critical pressure', c.pElastic / BAR, 'bar', 'Sum of the ring buckling pressures 3EI/R³ of carcass and pressure armour — the upper bound for a perfectly round pipe.'),
        kpi('collapseDepth', 'Collapse water depth', c.collapseDepth, 'm', 'Depth of sea water giving the collapse pressure on an empty pipe.'),
        kpi('allowableDepth', 'Allowable water depth', c.allowableDepth, 'm', 'Collapse depth multiplied by the permissible utilisation; the deepest water the empty pipe may see.'),
        kpi('uc', 'Collapse utilisation', uc, '-', 'Net external pressure over the permissible fraction of the collapse pressure. Above 1.0 the pipe is not acceptable at this depth.', 1, 'max'),
        kpi('factor', 'Permissible utilisation', R.factor, '-', 'Fraction of the collapse pressure that may be used; API 17J makes it depth dependent (recalled default — confirm).'),
        kpi('paShare', 'Pressure-armour share of elastic resistance', R.pa ? el.contributions[1] / el.pCr * 100 : 0, '%', 'Part of the elastic collapse resistance provided by the pressure armour. It is lost entirely if the annulus floods.'),
        kpi('tEq', 'Carcass equivalent ring thickness', R.tEq * 1e3, 'mm', 'Thickness of the solid ring with the same bending stiffness as the interlocked carcass.'),
      ],
      checks: [chk('Hydrostatic collapse, dry annulus', pDesign / BAR, R.factor * c.pCollapse / BAR, 'bar', `${STD_17B} §8.4; ${STD_17J} Table 6 (factor ${fmt(R.factor)})`,
        'Guards against the bore collapsing under external pressure when the pipe is empty or depressurised. A collapsed carcass blocks the bore and cannot be recovered.')],
      plots: [
        plot('line', 'Collapse pressure vs initial ovality', 'Ovality', '%', 'Collapse pressure', 'bar', [
          { name: 'Collapse pressure', x: ovs.map((o) => round(o)), y: pOv.map((p) => round(p / BAR)) },
          { name: 'Design external pressure', x: [0, round(ovs[ovs.length - 1])], y: [round(pDesign / BAR), round(pDesign / BAR)], dash: true },
          { name: 'This pipe', x: [inp.ovality], y: [round(c.pCollapse / BAR)] },
        ], 'Imperfection sensitivity: collapse pressure drops quickly with the first fractions of a per cent of ovality. The curve must stay above the dashed design pressure with the required margin.'),
        plot('line', 'Collapse utilisation vs water depth', 'Water depth', 'm', 'Utilisation', '-', [
          { name: 'Utilisation', x: ds.map((d) => round(d)), y: ucD.map((u) => round(u)) },
          { name: 'Limit', x: [0, round(ds[ds.length - 1])], y: [1, 1], dash: true },
        ], 'Utilisation against depth with the present permissible factor. The crossing of the limit line is the allowable depth.'),
      ],
      tables: [tbl('Ring contributions', [col('ring', 'Layer'), col('EI', 'Ring stiffness', 'N·m'), col('R', 'Mean radius', 'mm'), col('pcr', 'Elastic critical pressure', 'bar')],
        rings.map((r, i) => ({ ring: i === 0 ? 'Carcass' : 'Pressure armour', EI: round(r.EI), R: round(r.R * 1e3), pcr: round(el.contributions[i] / BAR) })),
        'Elastic buckling resistance of each supporting ring; with a dry annulus they add because the layers ovalise together.')],
      notes: [`Dry-annulus collapse at ${fmt(c.pCollapse / BAR)} bar (${fmt(c.collapseDepth)} m); design external pressure ${fmt(pDesign / BAR)} bar; utilisation ${fmt(uc)} with factor ${fmt(R.factor)}.`,
        'If the outer sheath is breached the pressure armour no longer helps — run the wet-collapse module.', GENERIC],
      advice: [
        adv('engineering', sev(uc), `Collapse utilisation is ${fmt(uc)} at ${fmt(inp.depth)} m; allowable depth ${fmt(c.allowableDepth)} m.`,
          uc > 1 ? `Increase the carcass ring stiffness to an equivalent thickness of about ${fmt(tReq * 1e3)} mm (now ${fmt(R.tEq * 1e3)} mm), or keep a minimum bore pressure of ${fmt((pDesign - R.factor * c.pCollapse) / BAR)} bar at all times.` : `Reserve of ${fmt((1 - uc) * 100)} %. Confirm the profile factor K and ovality with the manufacturer's collapse test at the minimum bend radius.`,
          `${STD_17B} §8.4; ${STD_17J} Table 6. ${CONFIRM}`),
        adv('economic', 'info', `The carcass is ${fmt(carSteel)} t of stainless strip per kilometre and is ${fmt(uc * 100)} % utilised against collapse.`,
          uc < 0.6 ? `Collapse capacity scales with t_eq³: an equivalent thickness of ${fmt(tReq * 1e3)} mm would suffice, a strip-mass saving of roughly ${fmt(carSteel * (1 - tReq / R.tEq))} t/km — worth a manufacturer optimisation if depth is fixed.` : 'The carcass is matched to the water depth; the remaining lever is operational — guaranteeing a minimum bore pressure avoids a heavier carcass.',
          'Whole-life costing per ISO 15686-5; AACE Class 4 quantity-based estimate',
          `Stainless carcass strip is the most carbon-intensive steel in the pipe: ${fmt(carSteel)} t/km ≈ ${fmt(carSteel * CO2_STEEL)} t CO₂e/km at the carbon-steel screening factor of ${CO2_STEEL} t CO₂e/t (stainless grades are typically two to three times higher).`),
      ],
      publish: { 'flexible.collapseDepth': c.collapseDepth },
    };
  },
  verify: [{
    // Elastic critical pressure = Σ 3EI/R³. Carcass: E = 193 GPa, t_eq = 6 mm → EI = 193e9·0.006³/12 = 3 474 N·m,
    // R = 76.2 + 3.5 = 79.7 mm → 3·3474/0.0797³ = 20.586 MPa. Pressure armour: E = 207 GPa, t_eq = 6.4 mm →
    // EI = 207e9·0.0064³/12 = 4 521.98 N·m, R = 95.2 mm → 3·4521.98/0.0952³ = 15.723 MPa. Sum = 36.309 MPa = 363.1 bar.
    name: 'Elastic collapse of two concentric rings, p = Σ 3EI/R³',
    source: 'Timoshenko & Gere, Theory of Elastic Stability, Ch. 7 (Levy ring); hand calculation',
    inputs: { ringMethod: 'equivalent', tEq: 6, paTeq: 6.4, ovality: 0 },
    expect: [{ kpi: 'pElastic', value: 363.1, tol: 0.002 }],
  }],
};
