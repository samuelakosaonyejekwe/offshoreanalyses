// Static configuration of a dynamic power cable between a floating platform and the seabed: free-hanging catenary
// or lazy wave, with in-plane current and near / far platform offsets, on the 2-D geometrically nonlinear rod solver.
import { riserOffsetCases } from '../../engine/lines/riser.js';
import { lineInputs, configInputs, currentInputs, routeInput, solveConfig, applyRoute, rodCurrent, shape3d, adv, check, worst, sev, fmt, round, run, thin, scaled, G,
  CO2_POLYMER, LIMIT_2D, STD_CABLE, vesselCo2Text } from './_common.js';

export default {
  id: 's13-static-configuration', suite: 13, title: 'Static cable configuration',
  summary: 'Installed shape, effective tension, curvature, clearances and end reactions of a dynamic cable in its mean, near and far positions.',
  covers: [1], fidelity: 'numerical',
  standards: [`${STD_CABLE} (dynamic cable global analysis)`, 'CIGRE TB 862 (2022) Recommendations for mechanical testing of submarine cables for dynamic applications', 'DNV-RP-C205 (2021) §6.2 (Morison drag)', 'API RP 17B (2014) §7 (compliant configurations)'],
  models: ['Elastic rod theory with geometric nonlinearity (corotational Euler–Bernoulli beam)', 'Nonlinear catenary equations', 'Effective tension and apparent weight', 'Morison drag in steady current', 'Seabed contact springs with axial Coulomb friction'],
  limits: `${LIMIT_2D} A 3-D Cosserat rod, torsion and the coupled platform analysis are outside this module; buoyancy modules are smeared along their section and marine growth is entered through mass, diameter and drag.`,
  publishes: ['cable.topTension', 'cable.hangoffAngleDeg', 'cable.tdpTension', 'cable.minBendRadius', 'cable.dynamicLength'],
  searchTerms: 'dynamic power cable lazy wave floating wind static configuration',
  inputs: [
    ...lineInputs(), ...configInputs(), routeInput(), ...currentInputs(),
    { key: 'offsetPct', label: 'Platform offset (near and far)', unit: '% depth', type: 'number', default: 8, min: 0, max: 30, group: 'Platform offset', role: 'bc', uncertain: { dist: 'normal', cov: 0.15 },
      help: 'Static excursion of the hang-off towards (near) and away from (far) the anchor point, from the mooring analysis. Cross offsets are out of plane and not modelled.' },
    { key: 'Tallow', label: 'Allowable tension', unit: 'kN', type: 'number', default: 250, min: 2, max: 5000, group: 'Acceptance criteria', role: 'criteria', bind: 'cable.allowTension',
      help: 'From the tension–bending capacity of the cross-section or the manufacturer data sheet.' },
    { key: 'mbrAllow', label: 'Allowable minimum bend radius (dynamic)', unit: 'm', type: 'number', default: 3.2, min: 0.3, max: 60, group: 'Acceptance criteria', role: 'criteria', bind: 'cable.mbrDynamic',
      help: 'Typically about 20 × overall diameter under tension — manufacturer specific; confirm against the type test.' },
    { key: 'sagMin', label: 'Minimum sag-bend clearance above the seabed', unit: 'm', type: 'number', default: 5, min: 0, max: 200, group: 'Acceptance criteria', role: 'criteria', showIf: { key: 'config', in: ['lazy-wave'] } },
    { key: 'hogMin', label: 'Minimum depth of the hog bend below still water', unit: 'm', type: 'number', default: 30, min: 0, max: 500, group: 'Acceptance criteria', role: 'criteria', showIf: { key: 'config', in: ['lazy-wave'] },
      help: 'Keeps the buoyancy section below wave action and vessel draught.' },
    { key: 'nElem', label: 'Rod elements', type: 'int', default: 80, min: 30, max: 800, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Rod elements', levels: [40, 80, 160, 320], qoi: 'topTension' },
  compute(raw) {
    const inp = applyRoute(raw);
    const current = rodCurrent(inp, inp.rhoW);
    const { res, line, h, buoy, buoySec } = solveConfig(inp, { nElem: inp.nElem, current });
    const off = inp.offsetPct / 100 * inp.depth;
    const cases = off > 0 ? run(() => riserOffsetCases(res, [-off, off]), 'Offset cases did not solve (the line may be too short for the far offset or pile up in the near offset)') : [];
    const all = [{ label: 'Near', r: cases[0]?.result }, { label: 'Mean', r: res }, { label: 'Far', r: cases[1]?.result }].filter((c) => c.r);
    const wave = inp.config === 'lazy-wave';
    const rows = all.map((c) => {
      const r = c.r, n = r.s.length - 1;
      return { case: c.label, T: r.topTension, angle: r.hangoffAngleDeg, mbr: r.minBendRadius, sMbr: r.sMaxCurvature, minT: r.minTension,
        sag: wave && r.sag ? r.sag.z : null, hog: wave && r.hog ? r.hog.z : null, tdpX: r.tdp ? r.tdp.x : null, botT: r.Te[n], endH: Math.abs(r.endB.Fx) };
    });
    const mean = rows.find((r) => r.case === 'Mean');
    const maxT = Math.max(...rows.map((r) => r.T)), minMbr = Math.min(...rows.map((r) => r.mbr)), minT = Math.min(...rows.map((r) => r.minT));
    const tdpXs = rows.map((r) => r.tdpX).filter((v) => v != null);
    const tdpTravel = tdpXs.length > 1 ? Math.max(...tdpXs) - Math.min(...tdpXs) : 0;
    const sagMinV = wave ? Math.min(...rows.map((r) => r.sag ?? 0)) : null, hogMaxV = wave ? Math.max(...rows.map((r) => r.hog ?? 0)) : null;
    if (wave && !(sagMinV > 0)) throw new Error('No lazy wave forms in at least one position: the buoyancy section rests on the seabed or the wave is pulled flat. Increase buoyancy, shorten the upper catenary or reduce the offset.');
    const checks = [
      check('Maximum effective tension', maxT / 1e3, inp.Tallow, 'kN', `${STD_CABLE}; CIGRE TB 862 (tensile capacity from type test)`, 'Largest hang-off tension of the three positions against the allowable tension of the cable. Above it the armour wires and conductors are overstrained.'),
      check('Minimum bend radius', inp.mbrAllow, minMbr, 'm', `${STD_CABLE} (MBR under tension); CIGRE TB 862`, 'Demand is the allowable radius, capacity the tightest radius found along the cable in any position. Tighter bending damages the lead sheath and insulation screens and can birdcage the armour.'),
      check('No effective compression', Math.max(0, -minT) / 1e3 + 1e-9, Math.max(line.w * line.D * 50, 1) / 1e3, 'kN', 'DNV-RP-0360 (compression in dynamic cables to be avoided); Euler screening', 'Compressive effective tension against a small tolerance (apparent weight of 50 diameters). A cable carries almost no compression: it loops or kinks at touchdown or in the sag bend.'),
    ];
    if (wave) {
      checks.push(check('Sag-bend seabed clearance', inp.sagMin, Math.max(sagMinV, 1e-9), 'm', 'API RP 17B §7; project criterion', 'Demand is the required clearance, capacity the lowest sag-bend elevation. A sag bend touching the seabed abrades the serving and shifts the touchdown fatigue hot spot.'));
      checks.push(check('Hog-bend submergence', inp.hogMin, Math.max(inp.depth - hogMaxV, 1e-9), 'm', 'API RP 17B §7; project criterion', 'Demand is the required depth of the buoyancy section below the surface, capacity the actual depth of the hog bend in its highest position. Too shallow exposes the modules to wave loading and vessels.'));
    }
    const gov = worst(checks);
    const [s, , z, te, kap] = thin([res.s, res.x, res.z, scaled(res.Te, 1e-3), Array.from(res.curvature, Math.abs)], 300);
    const foamT = wave ? inp.foamRho * buoy.foamVolume * inp.buoyLength / inp.modPitch / 1e3 : 0;
    const nMod = wave ? Math.ceil(inp.buoyLength / inp.modPitch) : 0;
    const sh = (c) => { const [xx, zz] = thin([c.r.x, c.r.z], 200); return { name: c.label, x: xx, y: zz, dash: c.label !== 'Mean' }; };
    return {
      kpis: [
        { key: 'topTension', label: 'Top effective tension (mean position)', value: round(res.topTension / 1e3, 7), unit: 'kN', explain: 'Effective tension at the hang-off in the mean position with current: the load on the hang-off and the start of the tension profile.' },
        { key: 'maxTension', label: 'Maximum top tension over offsets', value: round(maxT / 1e3), unit: 'kN', limit: inp.Tallow, sense: 'max', explain: 'Largest hang-off tension of the near, mean and far positions; it must stay below the allowable tension of the cable.' },
        { key: 'hangoffAngle', label: 'Hang-off angle from vertical', value: round(mean.angle), unit: 'deg', explain: 'Static departure angle at the hang-off; with the dynamic angle range it sizes the bend stiffener.' },
        { key: 'minBendRadius', label: 'Minimum bend radius along the cable', value: round(minMbr), unit: 'm', limit: inp.mbrAllow, sense: 'min', explain: 'Tightest static radius in any position (sag, hog or touchdown bend). It must stay above the allowable dynamic bend radius.' },
        { key: 'minTension', label: 'Minimum effective tension', value: round(minT / 1e3), unit: 'kN', explain: 'Lowest effective tension anywhere along the cable; a negative value is effective compression, which a cable tolerates only marginally before it loops.' },
        { key: 'tdpTension', label: 'Tension at the seabed end', value: round(mean.botT / 1e3), unit: 'kN', explain: 'Effective tension at the anchor point in the mean position; it is the residual tension the seabed tether or the buried section must hold.' },
        { key: 'tdpX', label: 'Touchdown distance from the hang-off', value: round(mean.tdpX ?? inp.span), unit: 'm', explain: 'Horizontal position where the cable first meets the seabed in the mean position (the anchor point when the whole lower catenary is suspended).' },
        { key: 'tdpTravel', label: 'Touchdown travel near → far', value: round(tdpTravel), unit: 'm', explain: 'Distance the touchdown point moves between the near and far positions: the length of serving exposed to abrasion and the zone to keep free of rock and burial.' },
        { key: 'hangoffV', label: 'Hang-off vertical reaction', value: round(Math.abs(res.endA.Fz) / 1e3), unit: 'kN', explain: 'Vertical support reaction at the hang-off in the mean position: suspended cable weight less the buoyancy uplift and the seabed share.' },
        { key: 'anchorH', label: 'Maximum horizontal load at the seabed anchor point', value: round(Math.max(...rows.map((r) => r.endH)) / 1e3), unit: 'kN', explain: 'Horizontal support reaction the seabed tether, clump weight or burial must resist (far position usually governs).' },
        ...(wave ? [
          { key: 'sagClearance', label: 'Lowest sag-bend clearance', value: round(sagMinV), unit: 'm', limit: inp.sagMin, sense: 'min', explain: 'Height of the sag bend above the seabed in its lowest position, to be kept with margin for buoyancy loss and marine growth.' },
          { key: 'hogDepth', label: 'Shallowest hog-bend depth', value: round(inp.depth - hogMaxV), unit: 'm', limit: inp.hogMin, sense: 'min', explain: 'Water depth above the crest of the wave in its highest position.' },
          { key: 'netBuoyancy', label: 'Net buoyancy of the buoyant section', value: round(res.netBuoyancy / 1e3), unit: 'kN', explain: 'Uplift of the modules minus the weight of the cable they carry; it sets the size of the wave that decouples the touchdown from platform motion.' },
        ] : []),
      ],
      checks,
      plots: [
        shape3d('Installed configuration (mean position)', [{ x: res.x, z: res.z, v: scaled(res.Te, 1e-3) }], { valueLabel: 'Effective tension', valueUnit: 'kN', depth: inp.depth,
          explain: 'Static shape in the vertical plane through hang-off and anchor point (z = 0 at still water), coloured by effective tension. The rod model is planar, so the third coordinate is zero.' }),
        { type: 'line', title: 'Near, mean and far shapes', x: { label: 'Horizontal distance from the mean hang-off', unit: 'm' }, y: { label: 'Elevation above seabed', unit: 'm' }, series: all.map(sh),
          explain: 'The three static shapes in a common frame. Look at how far the touchdown point and the sag bend move and whether the far position pulls the wave flat.' },
        { type: 'line', title: 'Effective tension along the cable', x: { label: 'Arc length from hang-off', unit: 'm' }, y: { label: 'Tension', unit: 'kN' },
          series: [{ name: 'Effective tension', x: s, y: te }, { name: 'True tension (Te − pe·Ae)', x: s, y: te.map((v, i) => v - line.rhoW * G * Math.max(0, inp.depth - z[i]) * line.Ae / 1e3), dash: true }],
          explain: 'Effective tension governs shape and buckling: it falls from the hang-off and is lowest at the sag bend or touchdown; any part below zero is in effective compression. The dashed true tension is what strains the armour.' },
        { type: 'line', title: 'Curvature along the cable', x: { label: 'Arc length from hang-off', unit: 'm' }, y: { label: 'Curvature', unit: '1/m' },
          series: [{ name: 'Mean position', x: s, y: kap }, { name: 'Allowable (1/MBR)', x: [s[0], s[s.length - 1]], y: [1 / inp.mbrAllow, 1 / inp.mbrAllow], dash: true }],
          explain: 'Peaks mark the sag, hog and touchdown bends. Every peak must stay below the dashed allowable curvature; the hang-off itself needs a bend-stiffener analysis.' },
      ],
      tables: [{ title: 'Offset cases', explain: 'One row per platform position: top tension, hang-off angle, tightest bend radius and where it occurs, lowest tension, wave elevations and touchdown position. The governing value of each column feeds the checks.',
        columns: [{ key: 'case', label: 'Position' }, { key: 'T', label: 'Top tension', unit: 'kN' }, { key: 'angle', label: 'Hang-off angle', unit: 'deg' }, { key: 'mbr', label: 'Min bend radius', unit: 'm' }, { key: 'sMbr', label: 'at arc length', unit: 'm' }, { key: 'minT', label: 'Min tension', unit: 'kN' }, { key: 'sag', label: 'Sag elevation', unit: 'm' }, { key: 'hog', label: 'Hog elevation', unit: 'm' }, { key: 'tdpX', label: 'Touchdown x', unit: 'm' }],
        rows: rows.map((r) => ({ case: r.case, T: round(r.T / 1e3), angle: round(r.angle), mbr: round(r.mbr), sMbr: round(r.sMbr), minT: round(r.minT / 1e3), sag: r.sag == null ? null : round(r.sag), hog: r.hog == null ? null : round(r.hog), tdpX: r.tdpX == null ? null : round(r.tdpX) })) }],
      notes: [
        `Governing check: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        ...(raw.route?.points ? [`End positions taken from the imported route: horizontal distance ${fmt(inp.span, 5)} m, hang-off ${fmt(h, 5)} m above the anchor point.`] : []),
        `Apparent weight ${fmt(line.w)} N/m${wave ? `; buoyant section ${fmt(buoySec.w)} N/m over ${inp.buoyLength} m with ${nMod} modules` : ''}. Hang-off is ${fmt(h, 4)} m above the seabed.`,
        minT < 0 ? 'Effective compression is predicted in at least one position: the configuration is not acceptable as it stands.' : 'The whole cable stays in positive effective tension in all positions.',
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at utilisation ${fmt(gov.uc)}; top tension ${fmt(maxT / 1e3)} kN, tightest radius ${fmt(minMbr)} m, touchdown travel ${fmt(tdpTravel)} m.`,
          gov.uc > 1 ? (wave ? 'Re-tune the wave: lengthen the buoyancy section or reduce the module pitch to raise the sag bend, or lengthen the cable to relieve the far position; then re-run near and far.' : 'Lengthen the cable or move the anchor point closer to open the touchdown radius, or adopt a lazy wave to decouple the touchdown from platform offsets.')
            : 'Carry the top tension and hang-off angle into the dynamic-response and storm modules, and confirm cross-current and 3-D offsets in a dedicated global analysis before freezing the configuration.',
          `${STD_CABLE}; CIGRE TB 862 (2022); API RP 17B (2014) §7`),
        adv('economic', 'info', wave ? `The lazy wave needs ${nMod} buoyancy modules (${fmt(foamT)} t of syntactic foam) on ${fmt(inp.length)} m of dynamic cable.` : `The free-hanging catenary uses ${fmt(inp.length)} m of dynamic cable and no buoyancy hardware.`,
          wave ? `Each 10 m removed from the buoyancy section saves ${Math.round(10 / inp.modPitch)} modules and their offshore fitting time, provided the sag clearance (${fmt(sagMinV)} m against ${inp.sagMin} m) is kept; fitting typically adds about one vessel day per 100–150 modules.` : 'A catenary is the cheapest configuration but transfers platform motion to the touchdown; compare the fatigue life at touchdown before choosing it over a wave.',
          'AACE RP 18R-97 (quantity-driven estimate); whole-life costing ISO 15686-5',
          wave ? `Foam embodied carbon ≈ ${fmt(foamT * CO2_POLYMER)} t CO₂e (screening ${CO2_POLYMER} t CO₂e/t polymer). ${vesselCo2Text(Math.max(0.5, nMod / 120))}` : 'No buoyancy hardware: lowest material footprint, but a shorter fatigue life at touchdown can mean an earlier cable replacement — compare on a whole-life basis.'),
      ],
      publish: { 'cable.topTension': maxT, 'cable.hangoffAngleDeg': mean.angle, 'cable.tdpTension': mean.botT, 'cable.minBendRadius': minMbr, 'cable.dynamicLength': inp.length },
      quantities: { cable_km: inp.length / 1e3, polymer_t: round(foamT) },
    };
  },
  verify: [
    // Inextensible catenary with seabed contact (Irvine 1981 §1.2): a = H/w, suspended length S = sqrt(h(h + 2a)),
    // closure X = a·asinh(S/a) + (L − S). h = 138 m, L = 300 m, X = 200 m: a = 21.6 m gives S = 157.9 m,
    // a·asinh(S/a) = 21.6 × 2.688 = 58.1 m and 58.1 + 142.1 = 200.2 m, so a ≈ 21.5 m.
    // w = (63 − 1025·π·0.16²/4)·9.80665 = 415.7 N/m → Ttop = w·(a + h) = 415.7 × 159.5 = 66.3 kN.
    { name: 'Free-hanging catenary: Ttop = w·(a + h), a from X = a·asinh(S/a) + L − S',
      source: 'Irvine, Cable Structures (1981) §1.2; hand solution a ≈ 21.5 m for h = 138 m, L = 300 m, X = 200 m',
      inputs: { config: 'catenary', length: 300, span: 200, Us: 0, Ub: 0, offsetPct: 0, EI: 200, muAx: 0, nElem: 160, mbrAllow: 0.5, Tallow: 5000, OD: 160, mass: 63 },
      expect: [{ kpi: 'topTension', value: 66.3, tol: 0.015 }] },
  ],
};
