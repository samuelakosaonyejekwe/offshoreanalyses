// Annulus environment and pressure: gas permeation through the pressure sheath, annulus free volume, pressure
// build-up and venting, sweet/sour classification, condensed-water pH, outer-sheath loading by annulus pressure
// and rapid-gas-decompression screening.
import { permeabilityFromCm3, annulusPermeation, annulusFreeVolume, annulusPressureBuildUp, annulusClassification, annulusPH, rgdScreening, H2S_SOUR_LIMIT } from '../../engine/xsection/annulus.js';
import { stackInputs, stackFrom, kpi, chk, adv, plot, bar, tbl, col, run, need, fmt, round, linspace, sev, worst, UNBONDED, GENERIC, STD_17J, STD_17B, CONFIRM, BAR } from './_common.js';

export default {
  id: 's14-annulus',
  suite: 14,
  title: 'Annulus pressure and environment',
  summary: 'Estimates gas permeation into the annulus, the vent rate and interval, whether the annulus is sweet or sour, and the loading the annulus pressure puts on the outer sheath.',
  covers: [12],
  fidelity: 'analytical',
  standards: ['API RP 17B (5th ed.) §6 (annulus permeation), §8', 'API Spec 17J (4th ed.) §5.3 (venting, blistering qualification)', 'ISO 15156-2 / NACE MR0175 (sour-service threshold)'],
  models: ['Steady solution–diffusion permeation through a cylindrical sheath', 'Ideal-gas annulus pressure build-up and venting', 'Annulus free-volume balance', 'Carbonic-acid pH of condensed water', 'Rapid gas decompression screening against the qualified envelope'],
  limits: 'Steady-state permeation with constant coefficients and ideal-gas annulus; shielding by the pressure armour is a user factor. Permeability coefficients, the pH confinement shift and the blistering envelope are test data — defaults are indicative. Corrosion-fatigue of the armour in the resulting environment needs S–N data from qualified tests. ' + UNBONDED + ' ' + GENERIC,
  publishes: ['flexible.annulusH2S', 'flexible.ventInterval'],
  searchTerms: 'flexible pipe annulus permeation vent sour service armour wire',
  inputs: [
    ...stackInputs(),
    { key: 'length', label: 'Pipe length draining to one vent', unit: 'm', type: 'number', default: 1500, min: 1, max: 50000, group: 'Pipe construction', role: 'geometry' },
    { key: 'pInt', label: 'Bore pressure', unit: 'bar', type: 'number', default: 200, min: 1, max: 1500, group: 'Internal fluid', role: 'bc', uncertain: { dist: 'normal', cov: 0.05 } },
    { key: 'Tsheath', label: 'Pressure sheath temperature', unit: '°C', type: 'number', default: 60, min: 0, max: 150, group: 'Internal fluid', role: 'bc', help: 'Permeability rises steeply with temperature; enter coefficients valid at this temperature.' },
    { key: 'gas', label: 'Bore gas composition and sheath permeability', type: 'table', group: 'Internal fluid', role: 'load',
      columns: [{ key: 'name', label: 'Component (CH4, CO2, H2S, H2O …)', type: 'text' }, { key: 'y', label: 'Mole fraction', unit: '%', type: 'number' }, { key: 'Pe', label: 'Permeability', unit: 'cm³(STP)·cm/(cm²·s·bar)', type: 'number' }],
      default: [{ name: 'CH4', y: 88, Pe: 2e-9 }, { name: 'CO2', y: 10, Pe: 1.5e-8 }, { name: 'H2S', y: 0.01, Pe: 4e-8 }, { name: 'H2O', y: 1.99, Pe: 3e-7 }],
      help: 'Composition of the gas phase in contact with the sheath and the permeability coefficient of the sheath polymer for each component at the sheath temperature (indicative PA11 values — replace with supplier data).' },
    { key: 'areaFactor', label: 'Effective-area factor (armour shielding)', unit: '-', type: 'number', default: 1, min: 0.05, max: 1, group: 'Internal fluid', role: 'material', help: '1 ignores shielding by the pressure armour (conservative for vent rate).' },
    { key: 'Tann', label: 'Annulus gas temperature', unit: '°C', type: 'number', default: 30, min: 1, max: 100, group: 'Annulus boundaries', role: 'ic' },
    { key: 'pVent', label: 'Vent-valve opening pressure (absolute)', unit: 'bar', type: 'number', default: 3, min: 1.05, max: 30, group: 'Annulus boundaries', role: 'bc', help: 'Annulus pressure at which the vent valve opens.' },
    { key: 'pReseat', label: 'Annulus pressure after venting (absolute)', unit: 'bar', type: 'number', default: 1.013, min: 0.5, max: 29, group: 'Annulus boundaries', role: 'ic' },
    { key: 'pAmb', label: 'External pressure at the outer sheath (absolute)', unit: 'bar', type: 'number', default: 1.013, min: 1, max: 400, group: 'Annulus boundaries', role: 'bc', help: 'Topside section above water = atmospheric (worst case for the outer sheath).' },
    { key: 'sourWire', label: 'Armour wires qualified for sour service', type: 'bool', default: false, group: 'Environmental degradation', role: 'material', help: 'Sour-qualified wires are of lower strength; high-strength sweet wires must not see H₂S above the threshold.' },
    { key: 'sourLimit', label: 'Sour threshold, H₂S partial pressure', unit: 'Pa', type: 'number', default: H2S_SOUR_LIMIT, min: 10, max: 5000, group: 'Acceptance criteria', role: 'criteria', help: 'ISO 15156-2: 0.3 kPa (0.05 psi).' },
    { key: 'phShift', label: 'pH rise from iron saturation in confined water', unit: '-', type: 'number', default: 0, min: 0, max: 3, group: 'Environmental degradation', role: 'material', help: 'Confined annulus water saturates with iron and its pH rises; take the shift from annulus-environment tests (0 = pure condensed water).' },
    { key: 'dpBlow', label: 'Bore depressurisation step', unit: 'bar', type: 'number', default: 150, min: 0, max: 1500, group: 'Pressure history', role: 'operation', help: 'Pressure drop during an emergency blow-down.' },
    { key: 'tBlow', label: 'Depressurisation time', unit: 'min', type: 'number', default: 30, min: 0.1, max: 10000, group: 'Pressure history', role: 'operation' },
    { key: 'rgdRate', label: 'Qualified decompression rate', unit: 'bar/min', type: 'number', default: 70, min: 0.1, max: 1000, group: 'Acceptance criteria', role: 'criteria', help: 'Rate used in the sheath blistering qualification test (API 17J). ' + CONFIRM },
  ],
  compute(inp, ctx) {
    const S = stackFrom(inp, ctx);
    need(S.idx.sheath >= 0, 'The stack has no pressure sheath: permeation cannot be assessed.');
    need(S.stack.sealOuter > S.stack.sealInner, 'The stack has no outer sheath, so there is no closed annulus.');
    need(inp.pReseat < inp.pVent, 'The annulus pressure after venting must be below the vent opening pressure.');
    const rowsIn = (inp.gas || []).filter((g) => Number(g.y) > 0 && Number(g.Pe) > 0);
    need(rowsIn.length > 0, 'Enter at least one gas component with a positive mole fraction and permeability.');
    const ySum = rowsIn.reduce((a, g) => a + Number(g.y), 0);
    need(ySum <= 100.5, `Gas mole fractions add up to ${fmt(ySum)} % — they must not exceed 100 %.`);
    const Sh = S.stack.layers[S.idx.sheath];
    const comps = rowsIn.map((g) => ({ name: String(g.name).trim().toUpperCase(), Pe: permeabilityFromCm3(Number(g.Pe)).Pe, fugacityBore: Number(g.y) / 100 * inp.pInt * BAR }));
    const perm = run(() => annulusPermeation({ components: comps, ri: Sh.rIn, ro: Sh.rOut, length: inp.length, areaFactor: inp.areaFactor }), 'Permeation');
    const between = S.stack.layers.filter((L) => L.index > S.stack.sealInner && L.index < S.stack.sealOuter);
    const vol = annulusFreeVolume({ layers: between.map((L) => ({ rIn: L.rIn, rOut: L.rOut, fill: L.kind === 'helical' ? Math.min(1, L.nA / Math.cos(L.alpha) / (Math.PI * (L.rOut ** 2 - L.rIn ** 2))) : 1 })), length: inp.length });
    need(vol.volume > 0, 'The annulus has no free volume: add the armour layers between the sheaths.');
    const T = inp.Tann + 273.15;
    const bu = run(() => annulusPressureBuildUp({ molarRate: perm.totalMolarRate, volume: vol.volume, T, p0: inp.pReseat * BAR, pVent: inp.pVent * BAR }), 'Annulus pressure build-up');
    const cls = annulusClassification({ moleFractions: perm.moleFractions, pAnnulus: inp.pVent * BAR, TC: inp.Tann, sourLimit: inp.sourLimit });
    const pH2S = cls.partialPressures.H2S ?? 0, pCO2 = cls.partialPressures.CO2 ?? 0;
    const ph = pCO2 > 0 ? annulusPH({ pCO2, confinementShift: inp.phShift }).pH : 7;
    const O = S.stack.layers[S.stack.sealOuter];
    const hoop = Math.max(0, (inp.pVent - inp.pAmb) * BAR) * O.R / O.t;
    const rgd = inp.dpBlow > 0 ? rgdScreening({ pressureDrop: inp.dpBlow * BAR, duration: inp.tBlow * 60, qualifiedRate: inp.rgdRate * BAR / 60 }) : { rateUtilisation: 0, rate: 0 };
    const checks = [
      chk('Outer sheath hoop stress at vent pressure', hoop / 1e6, (O.fy ?? 20e6) / 1e6, 'MPa', `${STD_17J} §5.3 (vent system must open before the outer sheath is overloaded)`,
        'Guards against the annulus gas pressure bursting the outer sheath from inside when the vent valve is slow or blocked. A split outer sheath floods the annulus.'),
      chk('Depressurisation rate vs blistering qualification', rgd.rate * 60 / BAR, inp.rgdRate, 'bar/min', `${STD_17J} (blistering / rapid gas decompression qualification)`,
        'Guards against blistering of the pressure sheath when dissolved gas cannot escape during a fast blow-down.'),
    ];
    if (!inp.sourWire) checks.push(chk('Annulus H₂S partial pressure vs sour threshold (sweet-service wires)', pH2S, inp.sourLimit, 'Pa', 'ISO 15156-2 §7 / NACE MR0175 (0.3 kPa threshold)',
      'Guards against sulphide stress cracking and hydrogen-induced cracking of high-strength armour wires. Above the threshold the wires must be sour-qualified.'));
    const gov = worst(checks);
    const ts = linspace(0, Math.min(bu.ventInterval, 3.2e9), 2);
    const dayMol = 86400;
    const cat = perm.rates.map((q) => q.name);
    return {
      kpis: [
        kpi('permRate', 'Gas permeation into the annulus', perm.totalMolarRate * dayMol, 'mol/day', 'Total gas passing through the pressure sheath into the annulus over the drained length.'),
        kpi('ventRate', 'Average vent flow', bu.ventRateStd * dayMol * 1e3, 'L/day (std)', 'Time-averaged gas volume released at the vent at standard conditions; a measured vent rate well above this points to a leaking sheath.'),
        kpi('ventInterval', 'Time between vent openings', bu.ventInterval / 3600, 'h', 'Time for permeated gas to raise the annulus from the re-seat pressure to the vent opening pressure.'),
        kpi('freeVolume', 'Annulus free volume', vol.volume * 1e3, 'L', 'Gas volume in the gaps of the armour layers over the drained length — the reference for annulus volume tests.'),
        kpi('pH2S', 'H₂S partial pressure in the annulus', pH2S, 'Pa', `Partial pressure of H₂S at vent pressure. At or above ${inp.sourLimit} Pa the annulus is sour and the armour wires must be sour-qualified.`),
        kpi('sourRatio', 'H₂S partial pressure / sour threshold', pH2S / inp.sourLimit, '-', 'Above 1.0 the annulus is classed as sour service under ISO 15156-2.'),
        kpi('pCO2', 'CO₂ partial pressure in the annulus', pCO2 / 1e3, 'kPa', 'Partial pressure of CO₂ at vent pressure; with condensed water it makes the annulus corrosive (sweet).'),
        kpi('pH', 'pH of condensed annulus water', ph, '-', 'pH from carbonic-acid equilibrium at the CO₂ partial pressure, plus the entered confinement shift. Lower pH means more severe corrosion fatigue.'),
        kpi('sheathUc', 'Outer sheath utilisation at vent pressure', checks[0].uc, '-', 'Hoop stress in the outer sheath from the annulus over-pressure over the polymer yield strength.', 1, 'max'),
        kpi('rgdUc', 'Decompression-rate utilisation', rgd.rateUtilisation, '-', 'Blow-down rate over the rate covered by the blistering qualification. Above 1.0 the event is outside the qualified envelope.', 1, 'max'),
      ],
      checks,
      plots: [
        bar('Permeation by component', cat, 'Molar rate', 'mol/day', [{ name: 'Permeation', y: perm.rates.map((q) => round(q.molarRate * dayMol)) }],
          'Which gases reach the annulus. Water vapour and CO₂ permeate far faster than methane, so the annulus gas differs from the bore gas.'),
        bar('Annulus gas composition', cat, 'Mole fraction', '%', [{ name: 'Annulus', y: cat.map((n) => round((perm.moleFractions[n] ?? 0) * 100)) }, { name: 'Bore', y: comps.map((c) => round(c.fugacityBore / (inp.pInt * BAR) * 100)) }],
          'Steady annulus composition (permeate stream) against the bore gas. H₂S and CO₂ are enriched in the annulus relative to the bore.'),
        plot('line', 'Annulus pressure between vents', 'Time', 'h', 'Annulus pressure', 'bar', [{ name: 'Pressure', x: ts.map((t) => round(t / 3600)), y: [inp.pReseat, round(inp.pReseat + bu.pressureRate * ts[1] / BAR)] }, { name: 'Vent opens', x: [0, round(ts[1] / 3600)], y: [inp.pVent, inp.pVent], dash: true }],
          'Linear pressure build-up from the re-seat pressure to the vent opening pressure. A faster rise in service indicates a leak or a flooded (smaller) free volume.'),
      ],
      tables: [tbl('Annulus gas', [col('gas', 'Component'), col('pBore', 'Bore partial pressure', 'bar'), col('rate', 'Permeation', 'mol/day'), col('y', 'Annulus mole fraction', '%'), col('pp', 'Annulus partial pressure', 'Pa')],
        perm.rates.map((q, i) => ({ gas: q.name, pBore: round(comps[i].fugacityBore / BAR), rate: round(q.molarRate * dayMol), y: round((perm.moleFractions[q.name] ?? 0) * 100), pp: round(cls.partialPressures[q.name] ?? 0) })),
        'One row per gas component: driving partial pressure in the bore, permeation rate and the resulting annulus composition at vent pressure.')],
      notes: [`Annulus classified as ${cls.classification}${cls.waterCondenses ? ', with water condensing' : ''}; H₂S ${fmt(pH2S)} Pa against ${inp.sourLimit} Pa.`, `Vent opens about every ${fmt(bu.ventInterval / 3600)} h; average vent flow ${fmt(bu.ventRateStd * dayMol * 1e3)} L/day (std).`,
        'There are no built-in sweet or sour armour-wire S–N curves: fatigue in this environment needs qualified test data.', GENERIC],
      advice: [
        adv('engineering', cls.classification === 'sour' && !inp.sourWire ? 'fail' : sev(gov.uc), `Annulus is ${cls.classification}: H₂S ${fmt(pH2S)} Pa (${fmt(pH2S / inp.sourLimit * 100)} % of the sour threshold), CO₂ ${fmt(pCO2 / 1e3)} kPa, condensed-water pH ${fmt(ph)}.`,
          cls.classification === 'sour' && !inp.sourWire ? 'Specify sour-qualified armour wire, or a sheath material / H₂S scavenging layer that keeps the annulus below the threshold, and use sour S–N data in the fatigue assessment.' : 'Use armour-wire S–N data qualified for this annulus environment in the fatigue modules, and record the vent rate at commissioning as the integrity baseline.',
          `ISO 15156-2 (2020) §7; ${STD_17B} §6; ${STD_17J} §5.3`),
        adv('economic', 'info', `Expected vent flow is ${fmt(bu.ventRateStd * dayMol * 1e3)} L/day (std) with an opening every ${fmt(bu.ventInterval / 3600)} h.`,
          'Fit vent-gas flow monitoring and trend it: a step change detects a sheath leak or annulus flooding months before a failure, at the cost of a flow meter rather than a riser replacement and its downtime.',
          'Risk-based inspection per DNV-RP-G101 / API RP 580; value-of-information of condition monitoring',
          `Vented permeate is mostly methane and CO₂: about ${fmt((perm.rates.find((q) => q.name === 'CH4')?.molarRate ?? 0) * 16.04e-3 * 86400 * 365)} kg of methane per year from this length — small, but routing vents to a closed drain or flare avoids the release.`),
      ],
      publish: { 'flexible.annulusH2S': pH2S, 'flexible.ventInterval': bu.ventInterval },
    };
  },
  verify: [{
    // One component: n = Pe·Δp·2πL/ln(ro/ri). Pe = 1e-8 cm³(STP)·cm/(cm²·s·bar) = 1e-8·(1/22414 mol)·(1e-2 m)/(1e-4 m²·1e5 Pa)
    // = 4.4615e-16 mol·m/(m²·s·Pa). Pure CO₂ at 200 bar (Δp = 2e7 Pa), L = 1000 m, sheath radii 83.2 → 91.2 mm:
    // ln(91.2/83.2) = 0.091808. n = 4.4615e-16·2e7·2π·1000/0.091808 = 6.1067e-4 mol/s = 52.76 mol/day.
    name: 'Steady permeation of one gas through a cylindrical sheath',
    source: 'Crank, The Mathematics of Diffusion, §5.2 (steady radial diffusion); hand calculation',
    inputs: { length: 1000, pInt: 200, areaFactor: 1, gas: [{ name: 'CO2', y: 100, Pe: 1e-8 }] },
    expect: [{ kpi: 'permRate', value: 52.76, tol: 0.003 }],
  }],
};
