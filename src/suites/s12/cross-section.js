// Cross-sectional structural response: helical-tendon assembly of the lay-up into equivalent axial, torsional and
// bending stiffness, load sharing between elements, contact line loads and the stick–slip moment–curvature loop.
import { umbilicalLoadSharing } from '../../engine/xsection/umbilical.js';
import { criticalCurvature, momentCurvatureLoop } from '../../engine/xsection/helix.js';
import { sectionInputs, buildSection, sectionMasses, adv, check, fmt, round, run, need, G, CO2_STEEL, CO2_POLYMER, LIMIT_RIGID, STD_UMB } from './_common.js';

export default {
  id: 's12-cross-section', suite: 12, title: 'Cross-sectional structural response',
  summary: 'Equivalent axial, bending and torsional stiffness, mass, tension sharing, contact loads and stick–slip bending of the umbilical lay-up.',
  covers: [5], fidelity: 'reduced-order',
  standards: [`${STD_UMB} §7 (design analysis: load-sharing and stiffness)`, 'API RP 17B (2014) §8 (helical layer mechanics)'],
  models: ['Helical wire mechanics (Féret & Bournazel 1987; Witz & Tan 1992)', 'Equivalent axial, bending and torsional stiffness', 'Layer interaction: contact line loads and Coulomb stick–slip', 'Bilinear moment–curvature hysteresis'],
  limits: `${LIMIT_RIGID} Elements in a layer share radial load equally; local contact deformation, fillers' creep and a multilayer finite-element cross-section model are outside this module.`,
  publishes: ['umbilical.OD', 'umbilical.massPerM', 'umbilical.EA', 'umbilical.EI', 'umbilical.EIstick', 'umbilical.GJ', 'umbilical.wSub', 'umbilical.torquePerTension'],
  searchTerms: 'umbilical cross-section helical stiffness load sharing stick slip bending',
  inputs: [
    ...sectionInputs(),
    { key: 'T', label: 'Applied tension', unit: 'kN', type: 'number', default: 120, min: 0, max: 5000, group: 'Loading', role: 'bc', bind: 'umbilical.topTension',
      uncertain: { dist: 'lognormal', cov: 0.1 }, help: 'Effective tension at the section considered (hang-off value from the global analysis).' },
    { key: 'kappaAmp', label: 'Curvature amplitude of the bending cycle', unit: '1/m', type: 'number', default: 0.02, min: 1e-4, max: 0.5, group: 'Loading', role: 'bc',
      help: 'Amplitude used for the hysteresis loop (1/radius).' },
    { key: 'twist', label: 'End rotation', type: 'select', default: 'restrained', group: 'Loading', role: 'bc',
      options: [{ value: 'restrained', label: 'Twist restrained (terminations fixed against rotation)' }, { value: 'free', label: 'Free to twist (torque-free, e.g. swivel)' }] },
    { key: 'mu', label: 'Friction coefficient between elements', unit: '-', type: 'number', default: 0.2, min: 0, max: 0.6, group: 'Loading', role: 'material', uncertain: { dist: 'lognormal', cov: 0.3 },
      help: 'No normative value; 0.1–0.3 for steel on polymer. Calibrate against cyclic bending tests.' },
    { key: 'balanceLimit', label: 'Torque-balance limit |K12|/√(K11·K22)', unit: '-', type: 'number', default: 0.05, min: 0.001, max: 0.5, group: 'Acceptance criteria', role: 'criteria',
      help: 'Project criterion for a torque-balanced design (no normative value).' },
    { key: 'lengthKm', label: 'Manufactured length', unit: 'km', type: 'number', default: 5, min: 0.05, max: 200, group: 'Quantities', role: 'operation', bind: 'cable.length' },
  ],
  compute(inp, ctx) {
    const b = buildSection(inp, ctx), s = b.sec;
    const T = inp.T * 1e3, restrained = inp.twist === 'restrained';
    const ls = run(() => umbilicalLoadSharing(s, { tension: T, twistRestrained: restrained }), 'Load sharing');
    const coupling = Math.abs(s.K12) / Math.sqrt(s.K11 * (s.K22 || 1));
    // stick–slip per helical entry: available friction f = μ(q_in + q_out) per unit element length (helix.js shear-lag model)
    let Mf = 0, kCr = Infinity;
    const rows = s.elements.map((e, i) => {
      const r = ls.elements[i], c = e.component;
      let kc = null;
      if (e.R > 0 && c.EA > 0 && T > 0 && inp.mu > 0) {
        const f = inp.mu * (r.contactInner + r.contactOuter);
        if (f > 0) { const cc = criticalCurvature({ E: c.EA / c.A, A: c.A, R: e.R, alpha: e.alpha, n: e.count, frictionPerLength: f }); Mf += cc.frictionMoment; kc = cc.kappaCr; kCr = Math.min(kCr, kc); }
      }
      return { name: e.name, kind: e.kind, layer: e.layer + 1, count: e.count, R: round(e.R * 1e3), alpha: round(e.alpha * 180 / Math.PI), force: round(r.force / 1e3), share: round(r.share * 100), stress: round(r.stress / 1e6),
        qIn: round(r.contactInner), qOut: round(r.contactOuter), kCr: kc == null ? null : round(kc, 3) };
    });
    const loop = run(() => momentCurvatureLoop({ EIstick: s.EIstick, EIslip: Math.min(s.EIslip, 0.999 * s.EIstick), frictionMoment: Mf, kappaAmp: inp.kappaAmp }), 'Moment–curvature loop');
    const masses = sectionMasses(b), L = inp.lengthKm * 1e3;
    const steelT = masses.steel * L / 1e3, polyT = masses.polymer * L / 1e3;
    const shareSteelTubes = ls.elements.reduce((a, r, i) => a + (s.elements[i].kind === 'tube' ? r.share : 0), 0);
    const circle = (x0, y0, r, n = 24) => ({ x: Array.from({ length: n + 1 }, (_, k) => (x0 + r * Math.cos(2 * Math.PI * k / n)) * 1e3), y: Array.from({ length: n + 1 }, (_, k) => (y0 + r * Math.sin(2 * Math.PI * k / n)) * 1e3) });
    // layout drawing: elements spread evenly round each pitch circle in table order
    const byKind = {};
    const perLayer = {};
    for (const e of s.elements) (perLayer[e.layer] ||= []).push(e);
    for (const [, els] of Object.entries(perLayer)) {
      const n = els.reduce((a, e) => a + e.count, 0); let k = 0;
      for (const e of els) for (let j = 0; j < e.count; j++, k++) {
        const th = 2 * Math.PI * k / n;
        const sr = (byKind[e.kind] ||= { name: e.kind, x: [], y: [] });
        sr.x.push(round(e.R * Math.cos(th) * 1e3)); sr.y.push(round(e.R * Math.sin(th) * 1e3));
      }
    }
    const outline = circle(0, 0, b.OD / 2, 72);
    const kpis = [
      { key: 'OD', label: 'Overall diameter', value: round(b.OD * 1e3), unit: 'mm', explain: 'Outside diameter over the sheath: sets buoyancy, drag and the handling radii quoted as multiples of diameter.' },
      { key: 'massFilled', label: 'Mass in air, tubes filled', value: round(s.massFilled), unit: 'kg/m', explain: 'Mass per metre including the helix over-length of each element and the fluid in the tubes; drives top tension and reel payload.' },
      { key: 'wSub', label: 'Submerged weight, tubes filled', value: round(s.submergedFilled), unit: 'N/m', explain: 'Weight in sea water with interstices sealed. A weight-to-diameter ratio that is too low makes a dynamic umbilical wander in current; too high raises top tension.' },
      { key: 'wdRatio', label: 'Submerged weight / diameter', value: round(s.submergedFilled / b.OD / 1e3), unit: 'kN/m²', explain: 'Stability indicator of a dynamic or seabed umbilical: values above roughly 1–2 kN/m² resist current-induced excursion; lighter sections need ballast or armour.' },
      { key: 'EA', label: 'Axial stiffness, torque-free', value: round(s.EA / 1e6), unit: 'MN', explain: 'Tension per unit strain with the ends free to rotate. Upper estimate: radial contraction of fillers, not modelled, lowers the measured value.' },
      { key: 'EAfixed', label: 'Axial stiffness, twist restrained', value: round(s.EAfixed / 1e6), unit: 'MN', explain: 'Tension per unit strain with rotation prevented (K11 = Σ n·EA·cos³α); the value that applies between two fixed terminations.' },
      { key: 'GJ', label: 'Torsional stiffness (tension-free)', value: round(s.GJ / 1e3), unit: 'kN·m²', explain: 'Torque per unit twist. Low values mean installation twist is easily induced and must be controlled by the lay spread.' },
      { key: 'EIstick', label: 'Bending stiffness, stick', value: round(s.EIstick / 1e3), unit: 'kN·m²', explain: 'Stiffness while friction locks the helical elements together (very small curvature only); it governs small-amplitude wave-frequency bending.' },
      { key: 'EIslip', label: 'Bending stiffness, slip', value: round(s.EIslip / 1e3), unit: 'kN·m²', explain: 'Stiffness after the elements have slipped: sum of each element bending about its own axis. Use it for global shape and large-curvature bending.' },
      { key: 'coupling', label: 'Torque balance |K12|/√(K11·K22)', value: round(coupling), unit: '-', limit: inp.balanceLimit, sense: 'max', explain: 'Normalised tension–torsion coupling. Zero is a perfectly torque-balanced lay-up; large values rotate the umbilical or load the terminations in torsion when tension changes.' },
      { key: 'torque', label: restrained ? 'Reaction torque at the terminations' : 'Twist under tension', value: restrained ? round(ls.torque / 1e3) : round(ls.twist * 180 / Math.PI), unit: restrained ? 'kN·m' : 'deg/m', explain: restrained ? 'Torque the terminations must react at the applied tension because the lay-up is not perfectly balanced.' : 'Rotation per metre of a torque-free umbilical at the applied tension.' },
      { key: 'strain', label: 'Axial strain at the applied tension', value: round(ls.axialStrain * 100), unit: '%', explain: 'Elongation of the bundle; compare with the allowable strain of fibre-optic units (strain-free window) and copper conductors.' },
      { key: 'steelShare', label: 'Tension carried by steel tubes', value: round(shareSteelTubes * 100), unit: '%', explain: 'Fraction of the tension taken by the hydraulic tubes; the rest goes to armour, cables and sheath. Tubes carrying most of the tension leave less stress available for pressure.' },
      { key: 'kappaCr', label: 'Slip-onset curvature', value: Number.isFinite(kCr) ? round(kCr, 3) : 0, unit: '1/m', explain: 'Curvature at which the first helical element overcomes friction and starts to slide; below it the section bends with the stick stiffness.' },
      { key: 'damping', label: 'Hysteretic damping ratio of the cycle', value: round(loop.dampingRatio, 3), unit: '-', explain: 'Equivalent viscous damping from the friction loop at the chosen curvature amplitude; it damps wave-frequency bending response near the hang-off.' },
    ];
    const checks = [
      check('Torque balance', coupling, inp.balanceLimit, '-', 'Project criterion; ISO 13628-5 §7.2 (torque-balanced design in intent)', 'Guards against rotation of the umbilical and torsional loading of terminations when tension varies; an unbalanced lay-up can unlay, loop or bird-cage during installation.'),
    ];
    const bal = checks[0].uc <= 1;
    const tubeShareTxt = `${fmt(shareSteelTubes * 100)} % of the tension is carried by the steel tubes and ${fmt((1 - shareSteelTubes) * 100)} % by armour, cables and sheath`;
    return {
      kpis, checks,
      plots: [
        { type: 'scatter', title: 'Cross-section lay-out (element centres)', x: { label: 'y', unit: 'mm' }, y: { label: 'z', unit: 'mm' },
          series: [{ name: 'Outer sheath', ...outline }, ...Object.values(byKind)],
          explain: 'Centres of every element, spread evenly round each pitch circle in table order, with the outer sheath boundary. Check that the layers sit where the manufacturer drawing shows them and that the sheath encloses the outer layer.' },
        { type: 'bar', title: 'Tension share by component', categories: rows.map((r) => r.name), series: [{ name: 'Share of tension [%]', y: rows.map((r) => r.share) }],
          explain: 'How the applied tension divides among the component groups in proportion to n·EA·cos³α. The dominant bar is the load path whose stiffness uncertainty matters most.' },
        { type: 'bar', title: 'Contact line loads per element', categories: rows.filter((r) => r.layer > 0).map((r) => r.name), series: [{ name: 'Inner side [N/m]', y: rows.filter((r) => r.layer > 0).map((r) => r.qIn) }, { name: 'Outer side [N/m]', y: rows.filter((r) => r.layer > 0).map((r) => r.qOut) }],
          explain: 'Radial squeeze on each element from its own helix tension and from the layers outside it. Inner layers accumulate the load of everything wound over them: look for the largest bar when judging tube crushing and friction.' },
        { type: 'line', title: 'Moment–curvature hysteresis', x: { label: 'Curvature', unit: '1/m' }, y: { label: 'Bending moment', unit: 'kN·m' },
          series: [{ name: 'Stick–slip loop', x: loop.kappa, y: loop.moment.map((m) => m / 1e3) }],
          explain: 'Bilinear loop: steep stick branch up to the friction moment, then the slip stiffness. The enclosed area is the energy dissipated per cycle; a fat loop means strong hysteretic damping but also friction stress cycling in the elements.' },
      ],
      tables: [{ title: 'Element forces and contact loads', explain: 'One row per component type: helix geometry, force per element, share of the total tension, axial stress, contact line loads on its inner and outer sides, and the curvature at which it starts to slip. Use it to identify the most loaded element.',
        columns: [{ key: 'name', label: 'Component' }, { key: 'layer', label: 'Layer' }, { key: 'count', label: 'No.' }, { key: 'R', label: 'Pitch radius', unit: 'mm' }, { key: 'alpha', label: 'Lay angle', unit: 'deg' }, { key: 'force', label: 'Force per element', unit: 'kN' }, { key: 'share', label: 'Share', unit: '%' }, { key: 'stress', label: 'Axial stress', unit: 'MPa' }, { key: 'qIn', label: 'Contact inner', unit: 'N/m' }, { key: 'qOut', label: 'Contact outer', unit: 'N/m' }, { key: 'kCr', label: 'Slip curvature', unit: '1/m' }],
        rows: rows.map((r) => ({ ...r, layer: r.layer === 0 ? 'centre' : r.layer })) }],
      notes: [
        `${tubeShareTxt}.`,
        `Stick bending stiffness is ${fmt(s.EIstick / s.EIslip)} × the slip stiffness; friction moment at ${fmt(inp.T)} kN tension is ${fmt(Mf)} N·m and slip starts at a bend radius of about ${Number.isFinite(kCr) ? fmt(1 / kCr) : '∞'} m.`,
        bal ? 'The lay-up is torque balanced within the project limit.' : `The lay-up is not torque balanced: reaction torque per unit tension is ${fmt(s.torquePerTension * 1e3)} mm.`,
      ],
      advice: [
        adv('engineering', bal ? 'ok' : 'warn', bal ? `Tension–torsion coupling ${fmt(coupling)} is inside the limit of ${inp.balanceLimit}.` : `Tension–torsion coupling ${fmt(coupling)} exceeds the limit of ${inp.balanceLimit} (${fmt(Math.abs(ls.torque) / 1e3)} kN·m at ${fmt(inp.T)} kN).`,
          bal ? 'Keep the contra-wound lay; confirm EA, EI and GJ with a full-section tensile, bending and torsion test and replace the assembler values with the measured ones in the global analyses.' : 'Rebalance by adjusting the wire count or lay angle of the outer armour layer (Σ lay·n·EA·R·sinα·cos²α → 0), then repeat this run.',
          `${STD_UMB} §7 design analysis and verification tests (tensile, bend stiffness, torque balance)`),
        adv('engineering', 'info', `Axial stiffness ${fmt(s.EA / 1e6)} MN is an upper estimate (radially rigid bundle); ${tubeShareTxt}.`,
          'Use a reduced EA (typically the measured value) when the global analysis is sensitive to axial stiffness, e.g. snatch loads during lay.', 'Witz & Tan (1992); calibration to full-section tests per ISO 13628-5 §7 / §11'),
        adv('economic', 'info', `The lay-up holds ${fmt(masses.steel)} kg/m of steel and ${fmt(masses.polymer)} kg/m of polymer: ${fmt(steelT)} t steel and ${fmt(polyT)} t polymer over ${inp.lengthKm} km.`,
          `Material quantity scales directly with armour wire count: each armour wire removed saves ${fmt(L * 7850 * Math.PI * 0.0025 ** 2 / 1e3, 2)} t over the length if the tension capacity and weight-to-diameter ratio still pass. Carry the tonnages into the cost estimate as quantities, priced at the manufacturer's rates.`,
          'AACE RP 18R-97 (quantity-based Class 3–4 estimate); whole-life costing ISO 15686-5',
          `Embodied carbon ≈ ${fmt(steelT * CO2_STEEL)} t CO₂e for steel (screening ${CO2_STEEL} t/t; super-duplex tube is higher) plus ≈ ${fmt(polyT * CO2_POLYMER)} t CO₂e for polymers (screening ${CO2_POLYMER} t/t). Thermoplastic sheath and steel are separable and recyclable at decommissioning if the umbilical is recovered.`),
      ],
      publish: { 'umbilical.OD': b.OD, 'umbilical.massPerM': s.massFilled, 'umbilical.EA': s.EA, 'umbilical.EI': s.EIslip, 'umbilical.EIstick': s.EIstick, 'umbilical.GJ': s.GJ, 'umbilical.wSub': s.submergedFilled, 'umbilical.torquePerTension': s.torquePerTension },
      quantities: { steel_t: steelT, polymer_t: polyT, cable_km: inp.lengthKm },
    };
  },
  verify: [
    { name: 'Single straight tube: EA = E·π(D² − d²)/4',
      source: 'Axial elasticity (Hooke). A = π(31.8² − 26.6²)/4 = 238.51 mm²; EA = 200 GPa × 238.51 mm² = 47.70 MN',
      inputs: { sheathT: 0, layers: [], components: [{ name: 'Tube', kind: 'tube', layer: 0, count: 1, od: 31.8, wall: 2.6, E: 200, rho: 7800, fy: 550, p: 0, rhoF: 0 }] },
      expect: [{ kpi: 'EAfixed', value: 47.70, tol: 0.002 }, { kpi: 'EIslip', value: 5.1228, tol: 0.003 }] },
    // I = π(31.8⁴ − 26.6⁴)/64 = π(1 022 600 − 500 635)/64 = 25 622 mm⁴ → EI = 200e9 × 2.5622e-8 = 5.124 kN·m²
    { name: 'Single helical wire layer: K11 = n·EA·cos³α, stick EI = n·EA·R²·cos³α/2',
      source: 'Witz & Tan (1992) Marine Structures 5. 68 wires Ø5 mm (A = 19.635 mm²) at 15°, R = 59 mm: K11 = 68 × 3.927 MN × 0.90122 = 240.66 MN',
      inputs: { sheathT: 0, layers: [{ layer: 1, R: 59, alpha: 15, lay: 1 }], components: [{ name: 'Wire', kind: 'armour', layer: 1, count: 68, od: 5, wall: 0, E: 200, rho: 7850, fy: 1100, p: 0, rhoF: 0 }] },
      expect: [{ kpi: 'EAfixed', value: 240.66, tol: 0.002 }] },
  ],
};
