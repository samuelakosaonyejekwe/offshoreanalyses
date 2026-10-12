// Global static configuration: free-hanging catenary or lazy wave between a floater hang-off and a seabed
// termination, with in-plane current and near / far vessel offsets, on the 2-D rod solver.
import { riserOffsetCases } from '../../engine/lines/riser.js';
import { hangOffClamp } from '../../engine/xsection/ancillaries.js';
import { lineInputs, configInputs, currentInputs, solveConfig, applyRoute, rodCurrent, shape3d, adv, check, worst, sev, fmt, round, run, thin, scaled, G, CO2_POLYMER, LIMIT_2D, STD_UMB, vesselCo2Text } from './_common.js';

export default {
  id: 's12-static-configuration', suite: 12, title: 'Global static configuration',
  summary: 'Shape, effective tension, curvature, clearances and end loads of a dynamic umbilical in its mean, near and far positions.',
  covers: [1], fidelity: 'numerical',
  standards: [`${STD_UMB} §7.3 (global analysis)`, 'API RP 17B (2014) §7 (configurations)', 'DNV-RP-C205 (2021) §6.2 (Morison drag)'],
  models: ['Elastic rod theory with geometric nonlinearity (corotational Euler–Bernoulli beam)', 'Cable catenary equations', 'Effective tension', 'Morison drag in steady current', 'Linear seabed contact springs with axial friction'],
  limits: `${LIMIT_2D} A 3-D Cosserat rod and coupled vessel analysis are outside this module; buoyancy modules are smeared along the section.`,
  publishes: ['umbilical.topTension', 'umbilical.hangoffAngleDeg', 'umbilical.tdpTension', 'umbilical.terminationH', 'umbilical.minBendRadius', 'cable.length'],
  searchTerms: 'dynamic umbilical lazy wave static configuration',
  inputs: [
    ...lineInputs(), ...configInputs(),
    { key: 'route', label: 'Imported centreline / route (optional)', type: 'geometry', kinds: ['polyline'], default: null, group: 'Configuration', role: 'geometry',
      help: 'CSV / XYZ / DXF / GeoJSON polyline, z up, first point at the hang-off and last at the termination. Only the end positions are used (they replace the horizontal distance and the hang-off height); the shape between them is solved. Leave empty to use the form values.' },
    ...currentInputs(),
    { key: 'offsetPct', label: 'Vessel offset (near and far)', unit: '% depth', type: 'number', default: 8, min: 0, max: 20, group: 'Vessel offset', role: 'bc', uncertain: { dist: 'normal', cov: 0.15 },
      help: 'Static offset of the hang-off towards (near) and away from (far) the termination, from the mooring analysis. Cross offsets are out of plane and not modelled.' },
    { key: 'Tallow', label: 'Allowable tension', unit: 'kN', type: 'number', default: 350, min: 5, max: 5000, group: 'Acceptance criteria', role: 'criteria', bind: 'umbilical.allowTension',
      help: 'From the capacity curve at the governing curvature (tension-and-bending module) or the manufacturer data sheet.' },
    { key: 'mbrAllow', label: 'Allowable minimum bend radius (dynamic)', unit: 'm', type: 'number', default: 11, min: 0.5, max: 60, group: 'Acceptance criteria', role: 'criteria', bind: 'umbilical.mbrDynamic' },
    { key: 'sagMin', label: 'Minimum sag-bend clearance above the seabed', unit: 'm', type: 'number', default: 10, min: 0, max: 200, group: 'Acceptance criteria', role: 'criteria', showIf: { key: 'config', in: ['lazy-wave'] } },
    { key: 'hogMin', label: 'Minimum depth of the hog bend below still water', unit: 'm', type: 'number', default: 30, min: 0, max: 500, group: 'Acceptance criteria', role: 'criteria', showIf: { key: 'config', in: ['lazy-wave'] },
      help: 'Keeps the buoyancy section out of wave action and clear of vessel traffic.' },
    { key: 'termRating', label: 'Termination pull rating (horizontal)', unit: 'kN', type: 'number', default: 60, min: 1, max: 5000, group: 'Acceptance criteria', role: 'criteria',
      help: 'Allowable horizontal load on the subsea termination / its foundation.' },
    { key: 'nBolts', label: 'Hang-off clamp bolts', type: 'int', default: 8, min: 2, max: 64, group: 'Hang-off clamp', role: 'geometry' },
    { key: 'boltPreload', label: 'Bolt preload', unit: 'kN', type: 'number', default: 80, min: 1, max: 2000, group: 'Hang-off clamp', role: 'bc', design: true },
    { key: 'clampMu', label: 'Clamp / sheath friction coefficient', unit: '-', type: 'number', default: 0.2, min: 0.03, max: 0.8, group: 'Hang-off clamp', role: 'material', uncertain: { dist: 'lognormal', cov: 0.25 }, help: 'Lined clamp on a polyethylene sheath: 0.15–0.3; armour-pot terminations are rated by test instead.' },
    { key: 'clampL', label: 'Clamp length', unit: 'm', type: 'number', default: 0.6, min: 0.05, max: 5, group: 'Hang-off clamp', role: 'geometry' },
    { key: 'clampSF', label: 'Slip safety factor', unit: '-', type: 'number', default: 2, min: 1, max: 5, group: 'Hang-off clamp', role: 'criteria' },
    { key: 'clampP', label: 'Allowable clamp contact pressure', unit: 'MPa', type: 'number', default: 12, min: 0.5, max: 200, group: 'Hang-off clamp', role: 'criteria', help: 'Radial pressure the lay-up tolerates under the clamp (from the crush analysis or crush test).' },
    { key: 'nElem', label: 'Rod elements', type: 'int', default: 120, min: 30, max: 1200, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Rod elements', levels: [60, 120, 240, 480], qoi: 'topTension' },
  compute(raw, ctx) {
    const inp = applyRoute(raw);
    const current = rodCurrent(inp, inp.rhoW);
    const { res, line, h, buoy, buoySec } = solveConfig(inp, { nElem: inp.nElem, current });
    const off = inp.offsetPct / 100 * inp.depth;
    const cases = off > 0 ? run(() => riserOffsetCases(res, [-off, off]), 'Offset cases did not solve (the line may be too short for the far offset or pile up in the near offset)') : [];
    const all = [{ label: 'Near', r: cases[0]?.result }, { label: 'Mean', r: res }, { label: 'Far', r: cases[1]?.result }].filter((c) => c.r);
    const wave = inp.config === 'lazy-wave';
    const row = (c) => {
      const r = c.r, n = r.s.length - 1;
      const sag = wave && r.sag && r.hog ? r.sag.z : null, hog = wave && r.hog ? r.hog.z : null;
      return { case: c.label, topTension: r.topTension, angle: r.hangoffAngleDeg, mbr: r.minBendRadius, sMbr: r.sMaxCurvature, minT: r.minTension, sag, hog, tdpX: r.tdp ? r.tdp.x : null, tdpS: r.tdp ? r.tdp.s : null,
        termH: Math.abs(r.endB.Fx), termV: r.endB.Fz, botT: r.Te[n] };
    };
    const rows = all.map(row), mean = rows.find((r) => r.case === 'Mean');
    const maxT = Math.max(...rows.map((r) => r.topTension)), minMbr = Math.min(...rows.map((r) => r.mbr)), minT = Math.min(...rows.map((r) => r.minT));
    const maxTermH = Math.max(...rows.map((r) => r.termH));
    const tdpXs = rows.map((r) => r.tdpX).filter((v) => v != null);
    const tdpTravel = tdpXs.length ? Math.max(...tdpXs) - Math.min(...tdpXs) : 0;
    const sagMinV = wave ? Math.min(...rows.map((r) => (r.sag == null ? 0 : r.sag))) : null;
    const hogMaxV = wave ? Math.max(...rows.map((r) => (r.hog == null ? 0 : r.hog))) : null;
    if (wave && !(sagMinV > 0)) throw new Error('No lazy wave forms in at least one position: the buoyancy section rests on the seabed or the wave is pulled flat. Increase buoyancy, shorten the upper catenary or reduce the offset.');
    const checks = [
      check('Maximum effective tension', maxT / 1e3, inp.Tallow, 'kN', `${STD_UMB} §7.3 (tension within the capacity curve)`, 'Guards against overstressing tubes and armour in tension at the hang-off; exceeding it yields the most loaded tube or wire.'),
      check('Minimum bend radius', inp.mbrAllow, minMbr, 'm', `${STD_UMB} §7.3 (dynamic MBR)`, 'Demand is the allowable radius, capacity the tightest radius found along the line in any position; tighter bending than allowed overstresses tubes and can kink the bundle.'),
      check('Termination horizontal load', maxTermH / 1e3, inp.termRating, 'kN', 'Project termination / foundation rating', 'Horizontal pull the subsea termination must react. Exceeding the rating drags the termination unit or overloads its connector.'),
    ];
    if (wave) {
      checks.push(check('Sag-bend seabed clearance', inp.sagMin, Math.max(sagMinV, 1e-9), 'm', 'API RP 17B §7 (configuration clearances); project criterion', 'Demand is the required clearance, capacity the lowest sag-bend elevation. A sag bend touching the seabed abrades the sheath and changes the touchdown behaviour.'));
      checks.push(check('Hog-bend submergence', inp.hogMin, Math.max(inp.depth - hogMaxV, 1e-9), 'm', 'API RP 17B §7; project criterion', 'Demand is the required depth of the buoyancy section below the surface, capacity the actual depth of the hog bend in its highest position. Too shallow exposes the modules to wave loading and vessels.'));
    }
    // hang-off clamp: friction grip and contact pressure at the maximum top tension
    const clamp = run(() => hangOffClamp({ tension: maxT, nBolts: inp.nBolts, boltPreload: inp.boltPreload * 1e3, mu: inp.clampMu, diameter: line.D, clampLength: inp.clampL, safetyFactor: inp.clampSF, allowablePressure: inp.clampP * 1e6 }), 'Hang-off clamp');
    checks.push(check('Hang-off clamp slip', maxT / 1e3, clamp.slipCapacity / 1e3, 'kN', `${STD_UMB} §9 (terminations and hang-off); friction grip μ·N·faces/SF`, 'Maximum top tension against the factored friction capacity of the hang-off clamp. If it slips the umbilical drops through the hang-off and the tubes take the load at the topside termination.'));
    checks.push(check('Hang-off clamp contact pressure', clamp.contactPressure / 1e6, inp.clampP, 'MPa', `${STD_UMB} §7 (crush load at clamps)`, 'Radial pressure of the bolted clamp on the umbilical against the allowable crush pressure of the lay-up. Too much preload ovalises the tubes under the clamp.'));
    const comp = minT < 0;
    const gov = worst(checks);
    const [s, x, z, te, kap] = thin([res.s, res.x, res.z, scaled(res.Te, 1e-3), Array.from(res.curvature, Math.abs)], 300);
    const foamT = wave ? inp.foamRho * buoy.foamVolume * inp.buoyLength / inp.modPitch / 1e3 : 0;
    const nMod = wave ? Math.ceil(inp.buoyLength / inp.modPitch) : 0;
    const kpis = [
      { key: 'topTension', label: 'Top effective tension (mean position)', value: round(res.topTension / 1e3, 7), unit: 'kN', explain: 'Effective tension at the hang-off in the mean position with current: the load on the hang-off clamp and the start of the tension profile.' },
      { key: 'maxTension', label: 'Maximum top tension over offsets', value: round(maxT / 1e3), unit: 'kN', limit: inp.Tallow, sense: 'max', explain: 'Largest hang-off tension of the near, mean and far positions; it must stay inside the capacity curve of the cross-section.' },
      { key: 'hangoffAngle', label: 'Hang-off angle from vertical', value: round(mean.angle), unit: 'deg', explain: 'Static departure angle at the hang-off; with the dynamic angle range it sizes the bend stiffener or bellmouth.' },
      { key: 'minBendRadius', label: 'Minimum bend radius along the line', value: round(minMbr), unit: 'm', limit: inp.mbrAllow, sense: 'min', explain: 'Tightest static radius in any position (usually the sag, hog or touchdown bend). It must stay above the allowable dynamic MBR.' },
      { key: 'minTension', label: 'Minimum effective tension', value: round(minT / 1e3), unit: 'kN', limit: 0, sense: 'min', explain: 'Lowest effective tension anywhere; negative values mean effective compression, which a flexible umbilical tolerates only briefly before it buckles into loops.' },
      { key: 'tdpTension', label: 'Tension at touchdown / termination', value: round(mean.botT / 1e3), unit: 'kN', explain: 'Effective tension at the seabed end in the mean position; low values mean the touchdown zone is governed by bending rather than tension.' },
      { key: 'termH', label: 'Maximum horizontal load on the termination', value: round(maxTermH / 1e3), unit: 'kN', limit: inp.termRating, sense: 'max', explain: 'Support reaction the subsea termination must resist horizontally (far position usually governs).' },
      { key: 'hangoffV', label: 'Hang-off vertical reaction', value: round(Math.abs(res.endA.Fz) / 1e3), unit: 'kN', explain: 'Vertical support reaction at the hang-off in the mean position: the weight of the suspended line less the buoyancy-section uplift and the seabed share.' },
      { key: 'hangoffH', label: 'Hang-off horizontal reaction', value: round(Math.abs(res.endA.Fx) / 1e3), unit: 'kN', explain: 'Horizontal support reaction at the hang-off in the mean position; it loads the I-tube or porch laterally.' },
      { key: 'clampUc', label: 'Hang-off clamp utilisation', value: round(Math.max(clamp.slipUtilisation, clamp.pressureUtilisation)), unit: '-', limit: 1, sense: 'max', explain: 'Larger of the slip and contact-pressure utilisations of the hang-off clamp at the maximum top tension.' },
      { key: 'tdpTravel', label: 'Touchdown travel near → far', value: round(tdpTravel), unit: 'm', explain: 'Distance the touchdown point moves along the seabed between the near and far positions: the length of sheath exposed to abrasion and trenching.' },
      ...(wave ? [
        { key: 'sagClearance', label: 'Lowest sag-bend clearance', value: round(sagMinV), unit: 'm', limit: inp.sagMin, sense: 'min', explain: 'Height of the sag bend above the seabed in its lowest position. It must stay clear of the seabed with margin for buoyancy loss and marine growth.' },
        { key: 'hogDepth', label: 'Shallowest hog-bend depth', value: round(inp.depth - hogMaxV), unit: 'm', limit: inp.hogMin, sense: 'min', explain: 'Water depth above the crest of the wave in its highest position.' },
        { key: 'netBuoyancy', label: 'Net buoyancy of the buoyant section', value: round(res.netBuoyancy / 1e3), unit: 'kN', explain: 'Uplift of the modules minus the weight of the umbilical they carry; it sets the size of the wave and decouples the touchdown from vessel motion.' },
      ] : []),
    ];
    const plots = [
      shape3d('Configuration (mean position)', [{ x: res.x, z: res.z, v: scaled(res.Te, 1e-3) }], { valueLabel: 'Effective tension', valueUnit: 'kN', depth: inp.depth,
        explain: 'Static shape in the vertical plane through hang-off and termination (z = 0 at still water), coloured by effective tension. The rod model is planar, so the third coordinate is zero.' }),
      { type: 'line', title: 'Near, mean and far shapes', x: { label: 'Horizontal distance from hang-off (mean)', unit: 'm' }, y: { label: 'Elevation above seabed', unit: 'm' },
        series: all.map((c, i) => { const shift = i === 0 && all.length === 3 ? off : i === 2 ? -off : 0; const [xx, zz] = thin([c.r.x, c.r.z], 200); return { name: c.label, x: xx.map((v) => v + (c.label === 'Mean' ? 0 : shift)), y: zz, dash: c.label !== 'Mean' }; }),
        explain: 'The three static shapes in a common frame. Look at how far the touchdown point and the sag bend move and whether the far position pulls the wave flat.' },
      { type: 'line', title: 'Effective tension along the line', x: { label: 'Arc length from hang-off', unit: 'm' }, y: { label: 'Tension', unit: 'kN' }, series: [{ name: 'Effective tension', x: s, y: te }, { name: 'True wall tension (Te − pe·Ae)', x: s, y: te.map((v, i) => v - line.rhoW * G * Math.max(0, inp.depth - z[i]) * line.Ae / 1e3), dash: true }],
        explain: 'Effective tension governs shape and buckling: it falls from the hang-off, is lowest at the sag or touchdown and rises over the buoyant section; any part below zero is in effective compression. The dashed true wall tension (external pressure end-cap removed, tube pressures excluded) is what strains the section and can be negative at depth without any buckling tendency.' },
      { type: 'line', title: 'Curvature along the line', x: { label: 'Arc length from hang-off', unit: 'm' }, y: { label: 'Curvature', unit: '1/m' },
        series: [{ name: 'Mean position', x: s, y: kap }, { name: 'Allowable (1/MBR)', x: [s[0], s[s.length - 1]], y: [1 / inp.mbrAllow, 1 / inp.mbrAllow], dash: true }],
        explain: 'Peaks mark the sag, hog and touchdown bends. Every peak must stay below the dashed allowable curvature; the hang-off itself needs the bend-stiffener analysis.' },
    ];
    return {
      kpis, checks, plots,
      tables: [{ title: 'Offset cases', explain: 'One row per vessel position: top tension, hang-off angle, tightest bend radius and where it occurs, lowest tension, wave elevations, touchdown position and termination loads. The governing value of each column feeds the checks.',
        columns: [{ key: 'case', label: 'Position' }, { key: 'T', label: 'Top tension', unit: 'kN' }, { key: 'angle', label: 'Hang-off angle', unit: 'deg' }, { key: 'mbr', label: 'Min bend radius', unit: 'm' }, { key: 'sMbr', label: 'at arc length', unit: 'm' }, { key: 'minT', label: 'Min tension', unit: 'kN' }, { key: 'sag', label: 'Sag elevation', unit: 'm' }, { key: 'hog', label: 'Hog elevation', unit: 'm' }, { key: 'tdpX', label: 'Touchdown x', unit: 'm' }, { key: 'termH', label: 'Termination H', unit: 'kN' }],
        rows: rows.map((r) => ({ case: r.case, T: round(r.topTension / 1e3), angle: round(r.angle), mbr: round(r.mbr), sMbr: round(r.sMbr), minT: round(r.minT / 1e3), sag: r.sag == null ? null : round(r.sag), hog: r.hog == null ? null : round(r.hog), tdpX: r.tdpX == null ? null : round(r.tdpX), termH: round(r.termH / 1e3) })) }],
      notes: [
        `Governing check: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        ...(raw.route?.points ? [`End positions taken from the imported route: horizontal distance ${fmt(inp.span, 5)} m, hang-off ${fmt(h, 5)} m above the termination.`] : []),
        `Apparent weight ${fmt(line.w)} N/m${wave ? `; buoyant section ${fmt(buoySec.w)} N/m (net uplift) over ${inp.buoyLength} m with ${nMod} modules` : ''}. Hang-off is ${fmt(h, 4)} m above the seabed.`,
        comp ? 'Effective compression is predicted: the configuration is not acceptable as it stands.' : 'The whole line stays in positive effective tension in all three positions.',
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at utilisation ${fmt(gov.uc)}; top tension ${fmt(maxT / 1e3)} kN, tightest radius ${fmt(minMbr)} m, touchdown travel ${fmt(tdpTravel)} m.`,
          gov.uc > 1 ? (wave ? 'Re-tune the wave: lengthen the buoyancy section or reduce the module pitch to raise the sag bend, or lengthen the line to relieve the far position; then re-run near and far.' : 'Lengthen the line or move the termination closer to raise the bend radius at touchdown, or adopt a lazy-wave configuration to decouple the touchdown from offsets.')
            : 'Carry the top tension and hang-off angle into the bend-stiffener and vessel-motion modules, and confirm the result with cross-current and 3-D offsets in a dedicated global analysis before freezing the configuration.',
          `${STD_UMB} §7.3; API RP 17B (2014) §7`),
        adv('economic', 'info', wave ? `The lazy wave needs ${nMod} buoyancy modules (${fmt(foamT)} t of syntactic foam) on a ${fmt(inp.length)} m dynamic section.` : `The free-hanging catenary uses ${fmt(inp.length)} m of dynamic umbilical and no buoyancy hardware.`,
          wave ? `Each 10 m removed from the buoyancy section saves ${Math.round(10 / inp.modPitch)} modules and their offshore fitting time, provided the sag clearance (${fmt(sagMinV)} m against ${inp.sagMin} m) is kept; module fitting typically adds about a vessel day per 100–150 modules.` : 'A catenary is the lowest-cost configuration but transfers vessel motion to the touchdown; compare the fatigue life at touchdown before choosing it over a wave.',
          'AACE RP 18R-97 (quantity-driven estimate); whole-life costing ISO 15686-5',
          wave ? `Foam embodied carbon ≈ ${fmt(foamT * CO2_POLYMER)} t CO₂e (screening ${CO2_POLYMER} t CO₂e/t polymer). ${vesselCo2Text(Math.max(0.5, nMod / 120))}` : 'No buoyancy hardware: lowest material footprint, but a shorter fatigue life at touchdown can mean earlier replacement — compare on a whole-life basis.'),
      ],
      publish: { 'umbilical.topTension': maxT, 'umbilical.hangoffAngleDeg': mean.angle, 'umbilical.tdpTension': mean.botT, 'umbilical.terminationH': maxTermH, 'umbilical.minBendRadius': minMbr, 'cable.length': inp.length },
      quantities: { cable_km: inp.length / 1e3, polymer_t: foamT },
    };
  },
  verify: [
    // Inextensible catenary with seabed contact (Irvine 1981 §1.2): a = H/w, suspended length S = sqrt(h(h + 2a)),
    // horizontal projection a·asinh(S/a); closure X = a·asinh(S/a) + (L − S) with h = 340 m, L = 600 m, X = 330 m
    // gives a = 32.00 m (bisection). w = (40.9 − 1025·π·0.145²/4)·9.80665 = 235.11 N/m → Ttop = w(a + h) = 235.11 × 372.00 = 87.46 kN.
    { name: 'Free-hanging catenary: Ttop = w·(a + h), a from X = a·asinh(S/a) + L − S',
      source: 'Irvine, Cable Structures (1981) §1.2; hand solution a = 32.00 m for h = 340 m, L = 600 m, X = 330 m',
      inputs: { config: 'catenary', length: 600, span: 330, Us: 0, Ub: 0, offsetPct: 0, EI: 200, muAx: 0, nElem: 240, mbrAllow: 1, Tallow: 5000 },
      expect: [{ kpi: 'topTension', value: 87.46, tol: 0.01 }] },
  ],
};
