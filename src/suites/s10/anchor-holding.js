// Anchor holding-capacity check for the common anchor families: drag-embedment (design-chart power law),
// vertically loaded plate anchors and gravity (deadweight) anchors, with set-up, cyclic and proof-load factors.
import { DRAG_ANCHOR_INDICATIVE, dragAnchorUhc, plateAnchorCapacity, gravityAnchorCapacity } from '../../engine/geotech/anchors.js';
import { DEG, G, need, run, round, fmt, adv, check, worst, sev, linspace, steelCo2Text, CO2_STEEL } from './_common.js';
import { anchorLoadInputs, clayInputs, capacityFactorInputs, requiredSf, anchorLoad, designCheck, LIMIT_ANCHOR } from './_anchor.js';

const SOILS = Object.keys(DRAG_ANCHOR_INDICATIVE);
const on = (...v) => ({ key: 'anchorType', in: v });

export default {
  id: 's10-anchor-holding',
  suite: 10,
  title: 'Anchor holding capacity',
  summary: 'Compares the design line load at the anchor with the holding capacity of a drag, plate or gravity anchor and returns the required anchor size or embedment.',
  covers: [9],
  fidelity: 'screening',
  standards: ['API RP 2SK (3rd ed., 2005) §6.5, Annex on anchor capacity', 'DNV-RP-E301 (2017) fluke anchors', 'DNV-RP-E302 (2017) plate anchors', 'ISO 19901-7 (2013) §10.4'],
  models: ['Drag-anchor design-chart power law UHC = a·W^b (NCEL 1987)', 'Plate anchor bearing capacity Nc·su·A (Rowe & Davis 1982; Martin & Randolph 2001)', 'Gravity anchor sliding and uplift statics (Coulomb friction / undrained base shear)', 'Consolidation set-up and cyclic degradation factors', 'Installation proof-load criterion'],
  limits: `${LIMIT_ANCHOR} The drag-anchor power-law coefficients are indicative values of the kind printed on manufacturers' charts and must be replaced by the chart of the actual anchor; plate capacity assumes deep embedment in clay.`,
  publishes: ['mooring.anchorCapacity'],
  searchTerms: 'anchor holding capacity drag embedment plate anchor gravity anchor',
  inputs: [
    ...anchorLoadInputs(),
    { key: 'anchorType', label: 'Anchor type', type: 'select', default: 'drag', group: 'Anchor', role: 'geometry',
      options: [{ value: 'drag', label: 'Drag-embedment (fluke) anchor' }, { value: 'plate', label: 'Plate anchor (VLA / SEPLA)' }, { value: 'gravity', label: 'Gravity (deadweight) anchor' }],
      help: 'Drag anchors take horizontal load only; plate anchors take load normal to the fluke at any angle; gravity anchors rely on weight and base friction.' },
    { key: 'mass', label: 'Anchor mass', unit: 't', type: 'number', default: 28, min: 0.5, max: 2000, group: 'Anchor', role: 'geometry', design: true, showIf: on('drag', 'gravity'),
      help: 'Mass in air. For a gravity anchor this is the total ballasted mass.' },
    { key: 'dragSoil', label: 'Drag-anchor soil class', type: 'select', default: SOILS[1], options: SOILS.map((v) => ({ value: v, label: v })), group: 'Anchor', role: 'soil', showIf: on('drag') },
    { key: 'aCoef', label: 'Power-law coefficient a (0 = soil class value)', unit: 'kN/t^b', type: 'number', default: 0, min: 0, max: 2000, group: 'Anchor', role: 'soil', showIf: on('drag'),
      help: 'UHC[kN] = a·W[t]^b. Indicative defaults: 161 very soft clay, 229 medium clay, 324 sand/hard clay with b = 0.92. Replace with the manufacturer chart.' },
    { key: 'bExp', label: 'Power-law exponent b', unit: '-', type: 'number', default: 0.92, min: 0.5, max: 1.2, group: 'Anchor', role: 'soil', showIf: on('drag') },
    { key: 'plateA', label: 'Plate (fluke) area', unit: 'm²', type: 'number', default: 14, min: 0.5, max: 200, group: 'Anchor', role: 'geometry', design: true, showIf: on('plate') },
    { key: 'plateB', label: 'Plate width', unit: 'm', type: 'number', default: 3, min: 0.3, max: 15, group: 'Anchor', role: 'geometry', showIf: on('plate') },
    { key: 'zInstall', label: 'Installed depth of plate centre', unit: 'm', type: 'number', default: 22, min: 0, max: 100, group: 'Anchor', role: 'ic', design: true, showIf: on('plate'),
      help: 'Initial condition: embedment reached at the end of installation, before keying.' },
    { key: 'keying', label: 'Embedment lost during keying', unit: 'm', type: 'number', default: 1.5, min: 0, max: 20, group: 'Anchor', role: 'ic', showIf: on('plate'), help: 'Typically 0.25–1.0 plate heights for SEPLAs.' },
    { key: 'NcPlate', label: 'Plate bearing factor Nc', unit: '-', type: 'number', default: 11.42, min: 5, max: 15, group: 'Anchor', role: 'soil', showIf: on('plate'), help: '11.42 strip (exact), 12.42–13.11 circular plate (smooth–rough), deep embedment.' },
    { key: 'gravMode', label: 'Gravity-anchor base response', type: 'select', default: 'drained', options: [{ value: 'drained', label: 'Drained (sand): friction' }, { value: 'undrained', label: 'Undrained (clay): base shear' }], group: 'Anchor', role: 'soil', showIf: on('gravity') },
    { key: 'rhoAnchor', label: 'Gravity-anchor bulk density', unit: 'kg/m³', type: 'number', default: 4000, min: 1500, max: 7850, group: 'Anchor', role: 'material', showIf: on('gravity'), help: 'Concrete 2400, iron-ore ballasted box ≈ 3500–4500, solid steel 7850. Sets the submerged weight.' },
    { key: 'deltaDeg', label: 'Base interface friction angle', unit: 'deg', type: 'number', default: 28, min: 5, max: 45, group: 'Anchor', role: 'soil', bind: 'site.phiDeg', showIf: on('gravity') },
    { key: 'baseA', label: 'Base area', unit: 'm²', type: 'number', default: 80, min: 1, max: 2000, group: 'Anchor', role: 'geometry', showIf: on('gravity') },
    { key: 'Hkey', label: 'Shear-key / skirt passive resistance', unit: 'kN', type: 'number', default: 0, min: 0, max: 50000, group: 'Anchor', role: 'soil', showIf: on('gravity') },
    ...clayInputs(),
    ...capacityFactorInputs({ sf: 1.5 }),
    { key: 'proofFactor', label: 'Installation proof load / design load', unit: '-', type: 'number', default: 1.0, min: 0.5, max: 2, group: 'Capacity factors', role: 'operation',
      help: 'Drag anchors are proof-tensioned at installation; API RP 2SK asks for the maximum intact design load for permanent moorings in soft clay unless set-up is demonstrated.' },
    { key: 'proofAvail', label: 'Available installation tension', unit: 'kN', type: 'number', default: 3000, min: 0, max: 60000, group: 'Capacity factors', role: 'operation', help: 'Bollard pull or winch tension that the installation spread can deliver at the anchor.' },
  ],
  compute(inp) {
    const ld = anchorLoad(inp);
    const type = inp.anchorType;
    const sfKey = type === 'plate' ? 'plateAnchor' : type === 'drag' ? 'dragAnchor' : 'pileOrSuctionLateral';
    const sf = requiredSf(inp, sfKey);
    const su0 = inp.su0 * 1e3, k = inp.suGrad * 1e3;
    let cap, capLabel, sizeKpi, plot, steel = 0, concrete = 0, sizeAdvice;
    const need0 = ld.T * sf / (inp.setUp * inp.cyclic);            // static capacity needed
    if (type === 'drag') {
      const c = DRAG_ANCHOR_INDICATIVE[inp.dragSoil];
      const a = inp.aCoef > 0 ? inp.aCoef : c.a, b = inp.bExp;
      const r = run(() => dragAnchorUhc({ mass: inp.mass * 1e3, a, b }), 'Drag anchor');
      cap = r.UHC; capLabel = `power law a = ${a}, b = ${b}`;
      const req = (need0 / 1e3 / a) ** (1 / b);                    // tonnes
      sizeKpi = { key: 'sizeReq', label: 'Required anchor mass', value: round(req, 4), unit: 't', limit: inp.mass, sense: 'max', explain: 'Anchor mass at which the chart capacity, reduced by set-up and cyclic factors, equals the design load times the safety factor.' };
      const ms = linspace(Math.max(0.5, 0.3 * Math.min(inp.mass, req)), 1.8 * Math.max(inp.mass, req), 30);
      plot = { type: 'line', title: 'Holding capacity vs anchor mass', x: { label: 'Anchor mass', unit: 't' }, y: { label: 'Capacity / demand', unit: 'kN' },
        series: [{ name: 'Design capacity', x: ms.map((v) => round(v)), y: ms.map((m) => round(a * m ** b * inp.setUp * inp.cyclic)) }, { name: 'Demand × safety factor', x: [ms[0], ms[ms.length - 1]], y: [ld.T * sf / 1e3, ld.T * sf / 1e3], dash: true }, { name: 'Selected anchor', x: [inp.mass], y: [round(cap * inp.setUp * inp.cyclic / 1e3)] }],
        explain: 'Chart capacity of the anchor family against mass. The crossing with the dashed demand line is the smallest anchor that satisfies the safety factor.' };
      steel = inp.mass;
      sizeAdvice = { req, unit: 't', now: inp.mass, word: 'anchor mass' };
    } else if (type === 'plate') {
      need(inp.keying < inp.zInstall, 'The keying loss exceeds the installed depth: the plate would surface. Increase the installed depth.');
      const r = run(() => plateAnchorCapacity({ A: inp.plateA, B: inp.plateB, zInstall: inp.zInstall, keyingLoss: inp.keying, su0, k, Nc: inp.NcPlate }), 'Plate anchor');
      cap = r.T; capLabel = `Nc = ${inp.NcPlate}, su = ${fmt(r.su / 1e3)} kPa at ${fmt(r.zFinal)} m`;
      need(k > 0 || su0 > 0, 'The soil has no strength: enter a mudline strength or a gradient.');
      const suReq = need0 / (inp.NcPlate * inp.plateA);
      const zReq = k > 0 ? Math.max(0, (suReq - su0) / k) + inp.keying : (suReq <= su0 ? inp.keying : Infinity);
      sizeKpi = { key: 'sizeReq', label: 'Required installed embedment', value: round(Number.isFinite(zReq) ? zReq : 999, 4), unit: 'm', limit: inp.zInstall, sense: 'max', explain: 'Installed depth of the plate centre at which the bearing capacity after keying, set-up and cyclic factors equals the design load times the safety factor.' };
      const zs = linspace(inp.keying + 0.5, 1.8 * Math.max(inp.zInstall, Number.isFinite(zReq) ? zReq : inp.zInstall), 30);
      plot = { type: 'line', title: 'Holding capacity vs installed depth', x: { label: 'Installed depth of plate centre', unit: 'm' }, y: { label: 'Capacity / demand', unit: 'kN' },
        series: [{ name: 'Design capacity', x: zs.map((v) => round(v)), y: zs.map((z) => round(inp.NcPlate * (su0 + k * (z - inp.keying)) * inp.plateA * inp.setUp * inp.cyclic / 1e3)) }, { name: 'Demand × safety factor', x: [zs[0], zs[zs.length - 1]], y: [ld.T * sf / 1e3, ld.T * sf / 1e3], dash: true }, { name: 'Selected depth', x: [inp.zInstall], y: [round(cap * inp.setUp * inp.cyclic / 1e3)] }],
        explain: 'Plate capacity grows with depth because clay strength does. The crossing with the dashed line is the embedment the installation must achieve.' };
      steel = inp.plateA * 0.35;                                    // screening: ≈ 0.35 t of fabricated steel per m² of fluke
      sizeAdvice = { req: zReq, unit: 'm', now: inp.zInstall, word: 'installed depth' };
    } else {
      const vol = inp.mass * 1e3 / inp.rhoAnchor, Wsub = (inp.mass * 1e3 - 1025 * vol) * G;
      const r = run(() => gravityAnchorCapacity({ W: Wsub, angle: ld.angle, mode: inp.gravMode, delta: inp.deltaDeg * DEG, su: su0, A: inp.baseA, Hkey: inp.Hkey * 1e3, T: ld.T }), 'Gravity anchor');
      cap = r.Tult; capLabel = `${inp.gravMode}, submerged weight ${fmt(Wsub / 1e3)} kN`;
      const capOf = (m) => { const W = (m * 1e3 - 1025 * m * 1e3 / inp.rhoAnchor) * G; return gravityAnchorCapacity({ W, angle: ld.angle, mode: inp.gravMode, delta: inp.deltaDeg * DEG, su: su0, A: inp.baseA, Hkey: inp.Hkey * 1e3 }).Tult; };
      let lo = 0.1, hi = Math.max(inp.mass, 1);
      for (let i = 0; i < 40 && capOf(hi) < need0; i++) hi *= 1.6;
      const reachable = capOf(hi) >= need0;
      if (reachable) for (let i = 0; i < 60; i++) { const mid = 0.5 * (lo + hi); if (capOf(mid) < need0) lo = mid; else hi = mid; }
      const req = reachable ? hi : 99999;
      sizeKpi = { key: 'sizeReq', label: 'Required anchor mass', value: round(req, 4), unit: 't', limit: inp.mass, sense: 'max', explain: 'Total anchor mass at which sliding and uplift resistance, with set-up and cyclic factors, equal the design load times the safety factor. 99999 means weight alone cannot do it (undrained base shear governs: enlarge the base).' };
      const ms = linspace(0.3 * inp.mass, 1.8 * Math.max(inp.mass, Math.min(req, 10 * inp.mass)), 30);
      plot = { type: 'line', title: 'Holding capacity vs anchor mass', x: { label: 'Anchor mass', unit: 't' }, y: { label: 'Capacity / demand', unit: 'kN' },
        series: [{ name: 'Design capacity', x: ms.map((v) => round(v)), y: ms.map((m) => round(capOf(m) * inp.setUp * inp.cyclic / 1e3)) }, { name: 'Demand × safety factor', x: [ms[0], ms[ms.length - 1]], y: [ld.T * sf / 1e3, ld.T * sf / 1e3], dash: true }, { name: 'Selected anchor', x: [inp.mass], y: [round(cap * inp.setUp * inp.cyclic / 1e3)] }],
        explain: 'Gravity-anchor capacity against mass for the line inclination entered. A curve that flattens means base shear on clay, not weight, limits the capacity.' };
      if (inp.rhoAnchor > 5000) steel = inp.mass; else concrete = vol;
      sizeAdvice = { req, unit: 't', now: inp.mass, word: 'anchor mass' };
    }
    const dc = designCheck({ staticCapacity: cap, setUpFactor: inp.setUp, cyclicFactor: inp.cyclic, tension: ld.T, safetyFactor: sf, proofLoad: inp.proofAvail * 1e3, proofFactor: inp.proofFactor });
    const checks = [
      check('Holding capacity', ld.T * sf / 1e3, dc.resistance / 1e3, 'kN', `API RP 2SK (3rd ed.) anchor safety factors (${type}, ${inp.service}, ${inp.condition})`,
        'Guards against the anchor dragging or pulling out in the design storm: the factored line load must not exceed the capacity after set-up and cyclic effects.'),
      check('Installation proof tension available', inp.proofFactor * ld.T / 1e3, inp.proofAvail * 1e3 / 1e3, 'kN', 'API RP 2SK (3rd ed.) §6.5 (anchor test load); DNV-RP-E301 (2017) §5',
        'Guards against an anchor that was never pulled to its service load: without the proof tension the embedment and capacity are unverified and the anchor may drag on first storm loading.'),
    ];
    if (type === 'drag') checks.push({ name: 'Uplift on a drag anchor', demand: round(ld.angleDeg, 4), capacity: 0, unit: 'deg', uc: ld.V > 1e-6 * ld.T ? 99 : 0, ref: 'API RP 2SK (3rd ed.) §6.5',
      explain: 'A drag anchor must be pulled horizontally; any uplift load reduces embedment and capacity, so the line must keep grounded length in all design conditions.' });
    const g = worst(checks);
    const quantities = {};
    if (steel > 0) quantities.steel_t = round(steel);
    if (concrete > 0) quantities.concrete_m3 = round(concrete);
    const over = sizeAdvice.req < sizeAdvice.now;
    return {
      kpis: [
        { key: 'demand', label: 'Design line load at anchor', value: round(ld.T / 1e3, 5), unit: 'kN', explain: 'Resultant of the horizontal and vertical line load at the anchor, before the safety factor.' },
        { key: 'angle', label: 'Load inclination above horizontal', value: round(ld.angleDeg, 4), unit: 'deg', explain: 'Uplift angle of the line at the anchor; decides which anchor types are admissible.' },
        { key: 'capacity', label: 'Static holding capacity', value: round(cap / 1e3, 5), unit: 'kN', explain: 'Reference ultimate capacity of the selected anchor in the given soil, before set-up and cyclic factors.' },
        { key: 'resistance', label: 'Design resistance', value: round(dc.resistance / 1e3, 5), unit: 'kN', explain: 'Static capacity × set-up × cyclic factor: what the anchor can be counted on to hold in the design storm.' },
        { key: 'fs', label: 'Factor of safety', value: round(dc.fs, 4), unit: '-', limit: sf, sense: 'min', explain: 'Design resistance over design load. Compare with the required factor for the anchor type and condition.' },
        sizeKpi,
        { key: 'efficiency', label: 'Capacity / anchor weight', value: round(type === 'plate' ? cap / (steel * 1e3 * G) : cap / (inp.mass * 1e3 * G), 4), unit: '-', explain: 'Holding efficiency: ultimate capacity divided by anchor weight in air. Modern fluke anchors reach 20–50, gravity anchors stay below 1.' },
        { key: 'proof', label: 'Required proof tension', value: round(inp.proofFactor * ld.T / 1e3, 5), unit: 'kN', limit: inp.proofAvail, sense: 'max', explain: 'Tension to apply and hold at installation to demonstrate the anchor; must be within the installation spread\'s capability.' },
      ],
      checks,
      plots: [plot,
        { type: 'bar', title: 'Demand and capacity', categories: ['Design load', 'Load × safety factor', 'Design resistance', 'Static capacity'], series: [{ name: 'kN', y: [ld.T, ld.T * sf, dc.resistance, cap].map((v) => round(v / 1e3)) }],
          explain: 'The four numbers of the check side by side. The design resistance bar must be at least as tall as the factored load bar.' }],
      tables: [
        { title: 'Indicative drag-anchor capacity by soil class', columns: [{ key: 'soil', label: 'Soil class' }, { key: 'a', label: 'a', unit: 'kN/t^b' }, { key: 'uhc', label: `UHC of a ${inp.mass} t anchor`, unit: 'kN' }, { key: 'req', label: 'Mass required for this load', unit: 't' }],
          rows: SOILS.map((s) => { const c = DRAG_ANCHOR_INDICATIVE[s]; return { soil: s, a: c.a, uhc: round(c.a * inp.mass ** c.b), req: round((need0 / 1e3 / c.a) ** (1 / c.b)) }; }),
          explain: 'What a drag anchor of the entered mass would hold in each indicative soil class, and the mass each class would need for this load. Use it to see how sensitive the anchor size is to soil classification.' },
      ],
      notes: [`${type} anchor (${capLabel}); required safety factor ${sf} (${inp.sfRequired > 0 ? 'user value' : 'API RP 2SK table'}).`,
        `Set-up ${inp.setUp} and cyclic ${inp.cyclic} factors scale the static capacity to ${fmt(dc.resistance / 1e3)} kN.`,
        'Drag-anchor coefficients are indicative; use the manufacturer\'s chart or a site-specific trajectory analysis for design.'],
      advice: [
        adv('engineering', sev(g.uc), `Factor of safety ${fmt(dc.fs)} against ${sf} required on a ${fmt(ld.T / 1e3)} kN load at ${fmt(ld.angleDeg)}°; governing: ${g.name} (${fmt(g.uc)}).`,
          g.uc > 1 ? (g.name.startsWith('Holding') ? `Increase the ${sizeAdvice.word} from ${fmt(sizeAdvice.now)} to at least ${fmt(sizeAdvice.req)} ${sizeAdvice.unit}, or improve the soil factors by site testing.` : g.name.startsWith('Uplift') ? 'Remove the uplift (more ground chain) or change to a plate, suction or pile anchor.' : `Provide ${fmt(inp.proofFactor * ld.T / 1e3)} kN at the anchor during installation (larger vessel, or cross-tension opposite lines with the platform winches).`)
            : `The anchor is adequate; the smallest admissible ${sizeAdvice.word} is ${fmt(sizeAdvice.req)} ${sizeAdvice.unit}. Confirm the soil class or strength profile at each anchor location before fixing the size.`,
          'API RP 2SK (3rd ed., 2005) §6.5 and anchor safety-factor tables; DNV-RP-E301 / E302 (2017).'),
        adv('economic', 'info', over ? `The selected ${sizeAdvice.word} exceeds the requirement by ${fmt(100 * (sizeAdvice.now / sizeAdvice.req - 1))} %.` : `The ${sizeAdvice.word} must rise by ${fmt(100 * (sizeAdvice.req / sizeAdvice.now - 1))} % to pass.`,
          type === 'plate' ? 'Each extra metre of embedment is installation time, not steel: weigh deeper installation against a larger plate.' : `Per anchor the difference is ${fmt(Math.abs(sizeAdvice.now - sizeAdvice.req))} t; multiply by the number of anchors, and note that anchor-handling vessel size (deck load, bollard pull ${fmt(inp.proofFactor * ld.T / 1e3 / 9.81)} t for proof loading) steps up with anchor mass.`,
          'Quantity-based Class 4 estimate per AACE RP 18R-97; marine spread selection per DNV-ST-N001 (2021).',
          steel > 0 ? `${steelCo2Text(type === 'plate' ? steel : Math.abs(sizeAdvice.now - sizeAdvice.req) || steel)} High-efficiency anchors hold the same load with a fraction of the steel of a gravity anchor.` : `A ${fmt(concrete)} m³ ballasted anchor ≈ ${fmt(concrete * 0.3)} t CO₂e at a screening factor of 0.3 t CO₂e per m³ of concrete; an embedded anchor of a few tens of tonnes of steel (${CO2_STEEL} t CO₂e/t) holds the same load with a smaller seabed footprint.`),
      ],
      publish: { 'mooring.anchorCapacity': dc.resistance },
      quantities,
    };
  },
  verify: [
    { // UHC = a·W^b = 229·15^0.92; 15^0.92 = exp(0.92·ln 15) = exp(2.4914) = 12.078 → 2765.9 kN.
      name: 'Drag-anchor power law, 15 t anchor in medium clay',
      source: 'NCEL Techdata Sheet 83-08R (1987) power-law form; coefficients as printed in the module (a = 229, b = 0.92)',
      inputs: { anchorType: 'drag', mass: 15, dragSoil: 'high-efficiency / medium clay', aCoef: 0, bExp: 0.92 },
      expect: [{ kpi: 'capacity', value: 2765.9, tol: 0.002 }] },
    { // Deep strip plate: T = Nc·su·A with su = 5 + 1.5·20 = 35 kPa at the final depth 20 m, A = 10 m², Nc = 11.42 → 3997 kN.
      name: 'Plate anchor in clay, T = Nc·su·A',
      source: 'Rowe & Davis (1982) Géotechnique 32(1): Nc = 3π + 2 for a deep strip anchor',
      inputs: { anchorType: 'plate', plateA: 10, plateB: 2, zInstall: 21, keying: 1, NcPlate: 11.42, su0: 5, suGrad: 1.5 },
      expect: [{ kpi: 'capacity', value: 3997, tol: 0.002 }] },
    { // Gravity anchor on sand, horizontal pull: Tult = W'·tan δ. 200 t at 4000 kg/m³ → W' = 200e3·9.80665·(1 − 1025/4000) = 1458.7 kN;
      // tan 30° = 0.57735 → 842.2 kN.
      name: 'Gravity anchor sliding, T = W′·tan δ',
      source: 'Coulomb friction statics (hand calculation)',
      inputs: { anchorType: 'gravity', mass: 200, rhoAnchor: 4000, gravMode: 'drained', deltaDeg: 30, anchorH: 500, anchorV: 0, Hkey: 0 },
      expect: [{ kpi: 'capacity', value: 842.2, tol: 0.002 }] },
  ],
};
