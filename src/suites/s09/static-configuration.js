// Global static configuration of a steel catenary or lazy-wave riser: converged equilibrium shape under apparent
// weight, buoyancy, internal fluid, current and vessel offset, with the force and stress distributions along it.
import { solveWaveRiser, riserOffsetCases, catenaryRiser } from '../../engine/lines/riser.js';
import {
  G, DEG, need, run, fmt, round, thin, scaled, pipeFrom, currentInputs, steelTonnes, adv, check, kpi, plot, meshInput, finish,
  handWeight, stations, steelCo2Text, CO2_STEEL, LIMIT_2D, CONFIRM,
} from './_common.js';
import { scrInputs, criteriaInputs, buildScr, stressOf, currentOf, tdpNode } from './_scr.js';

const lazy = { key: 'config', in: ['lazy'] };

function solveLazy(inp, ctx, nElem, current) {
  const pp = pipeFrom(inp, ctx);
  const outer = pp.D + 2 * inp.coatT / 1e3;
  need(inp.buoyD / 1e3 > outer, `Buoyancy module diameter ${inp.buoyD} mm must exceed the coated pipe diameter ${fmt(outer * 1e3)} mm.`);
  const pb = pipeFrom(inp, ctx, { buoyancy: { Dout: inp.buoyD / 1e3, rho: inp.buoyRho, fill: inp.buoyFill } });
  need(pb.section.w < 0, `The buoyancy section is not buoyant (net weight ${Math.round(pb.section.w)} N/m). Increase module diameter or coverage, or use lighter foam.`);
  const h = inp.depth - inp.hangoffDepth, L = inp.upperLength + inp.buoyLength + inp.lowerLength;
  need(L > Math.hypot(inp.lazySpan, h) * 1.02, `Total riser length ${fmt(L)} m is too short to reach from the hang-off to the seabed end (${fmt(Math.hypot(inp.lazySpan, h))} m apart). Lengthen a section or reduce the span.`);
  const pressure = { rhoInt: inp.rhoInt, pIntRef: inp.pInt * 1e5, zIntRef: h, rhoExt: pp.rhoW, zSurface: inp.depth, nu: pp.mat.nu };
  const base = run(() => solveWaveRiser({
    pipe: pp.section, buoyant: pb.section, upperLength: inp.upperLength, buoyancyLength: inp.buoyLength, lowerLength: inp.lowerLength,
    top: { x: 0, z: h }, bottom: { x: inp.lazySpan, z: 0 }, seabed: { kn: inp.kn * 1e3, mu: inp.mu }, current, pressure, nElem,
  }), 'The lazy-wave riser did not solve');
  need(base.converged, 'The lazy-wave solution did not converge. Adjust the section lengths or the horizontal span.');
  return { pp, pb, h, base, pressure, length: L };
}

export default {
  id: 's09-static-configuration', suite: 9, title: 'Global static configuration',
  summary: 'Converged static shape, effective tension, moment, shear and curvature of a steel catenary or lazy-wave riser for mean, near and far vessel positions.',
  covers: [1], fidelity: 'numerical',
  standards: ['API RP 2RD (1998, R2013) §5 global analysis', 'DNV-ST-F201 (2018) Sec. 4', 'API RP 17B (2014) §7 riser configurations'],
  models: ['Corotational Euler–Bernoulli beam (geometrically nonlinear finite elements)', 'Catenary equations', 'Effective tension formulation', 'Winkler seabed springs with Coulomb axial friction', 'Morison current drag', 'Buoyancy-controlled lazy-wave configuration'],
  limits: `${LIMIT_2D} Torsional response is therefore not produced. Flexible-pipe layer behaviour, bend stiffeners and 3-D riser arrays need a dedicated 3-D nonlinear FE model.`,
  publishes: ['riser.topTension', 'riser.hangoffAngleDeg', 'riser.tdpTension', 'riser.reactionH', 'riser.reactionV', 'riser.suspendedLength', 'riser.length', 'riser.tdpCurvature', 'riser.maxMoment'],
  searchTerms: 'steel catenary riser static configuration lazy wave',
  inputs: [
    { key: 'config', label: 'Riser configuration', type: 'select', default: 'scr', group: 'Configuration', role: 'geometry',
      options: [{ value: 'scr', label: 'Steel catenary riser (free-hanging)' }, { value: 'lazy', label: 'Lazy-wave riser (distributed buoyancy)' }],
      help: 'Top-tensioned and free-standing hybrid risers are in "Effective tension and bending moment"; drilling risers in "Extreme storm and disconnect".' },
    ...scrInputs(),
    { key: 'upperLength', label: 'Upper catenary length', unit: 'm', type: 'number', default: 1150, min: 50, max: 5000, group: 'Lazy-wave sections', role: 'geometry', showIf: lazy, design: true },
    { key: 'buoyLength', label: 'Buoyancy section length', unit: 'm', type: 'number', default: 600, min: 20, max: 2000, group: 'Lazy-wave sections', role: 'geometry', showIf: lazy, design: true,
      help: 'Longer buoyancy raises the hog bend and the sag-bend clearance.' },
    { key: 'lowerLength', label: 'Lower length (incl. pipe on the seabed)', unit: 'm', type: 'number', default: 900, min: 50, max: 5000, group: 'Lazy-wave sections', role: 'geometry', showIf: lazy },
    { key: 'lazySpan', label: 'Horizontal distance hang-off to seabed end', unit: 'm', type: 'number', default: 1900, min: 50, max: 8000, group: 'Lazy-wave sections', role: 'bc', showIf: lazy },
    { key: 'buoyD', label: 'Buoyancy module outer diameter', unit: 'mm', type: 'number', default: 700, min: 150, max: 2500, group: 'Buoyancy system', role: 'geometry', showIf: lazy },
    { key: 'buoyRho', label: 'Module (syntactic foam) density', unit: 'kg/m³', type: 'number', default: 450, min: 250, max: 800, group: 'Buoyancy system', role: 'material', showIf: lazy,
      help: 'Depth-rated syntactic foam: about 400–550 kg/m³ for 1000–2000 m.' },
    { key: 'buoyFill', label: 'Fraction of the section length covered', unit: '-', type: 'number', default: 0.8, min: 0.1, max: 1, group: 'Buoyancy system', role: 'geometry', showIf: lazy },
    ...currentInputs(),
    { key: 'curDir', label: 'Current direction (in plane)', type: 'select', default: 'far', group: 'Current', role: 'bc',
      options: [{ value: 'far', label: 'Towards the touchdown side' }, { value: 'near', label: 'Towards the vessel side' }] },
    { key: 'offset', label: 'Vessel offset (near / far)', unit: 'm', type: 'number', default: 40, min: 0, max: 300, group: 'Top end (hang-off)', role: 'bc', bind: 'vessel.offsetMax',
      uncertain: { dist: 'normal', cov: 0.2 }, help: 'In-plane excursion of the hang-off from the mean position; both the near and the far case are solved.' },
    ...criteriaInputs(),
    { key: 'centreline', label: 'Surveyed / imported riser centreline (optional)', type: 'geometry', kinds: ['polyline'], default: null, group: 'Configuration', role: 'geometry',
      help: 'A 3-D polyline (CSV, XYZ, DXF, GeoJSON) is projected on its own vertical plane, overlaid on the computed shape and the largest deviation reported. Leave empty to use the parametric configuration.' },
    meshInput(160),
  ],
  mesh: { key: 'nElem', label: 'Beam elements', levels: [80, 160, 320, 640], qoi: 'topTension' },

  compute(inp, ctx) {
    // Form values are in engineering units (mm, bar, kPa, kN·m/deg); the helpers convert to SI.
    const sign = inp.curDir === 'far' ? 1 : -1;
    const current = currentOf(inp, sign);
    const isLazy = inp.config === 'lazy';
    const m = isLazy ? solveLazy(inp, ctx, inp.nElem, current) : buildScr(inp, ctx, { nElem: inp.nElem, current });
    const { pp, base, h, pressure } = m;
    need(inp.offset <= 0.25 * h, `A vessel offset of ${inp.offset} m is more than 25 % of the hang-off height (${fmt(h)} m): outside any practical riser design envelope. Reduce the offset.`);
    const cases = [{ name: 'Mean', res: base, offset: 0 }];
    if (inp.offset > 0) {
      const oc = run(() => riserOffsetCases(base, [-inp.offset, inp.offset]), 'Vessel-offset cases did not solve');
      oc.forEach((o) => { need(o.result.converged, `The ${o.offset < 0 ? 'near' : 'far'} offset case did not converge; reduce the offset.`); });
      cases.push({ name: `Near (−${inp.offset} m)`, res: oc[0].result, offset: -inp.offset }, { name: `Far (+${inp.offset} m)`, res: oc[1].result, offset: inp.offset });
    }
    const fy = pp.mat.fy, allow = inp.Cf * inp.CaAllow * fy;
    for (const c of cases) { c.stress = stressOf(c.res, pp, pressure); c.k = tdpNode(c.res); }
    const gov = cases.reduce((a, c) => (c.stress.max > a.stress.max ? c : a), cases[0]);
    const kT = tdpNode(base), n = base.s.length - 1;
    const length = base.s[n];
    const minTe = Math.min(...cases.map((c) => c.res.minTension));
    const maxCurv = Math.max(...cases.map((c) => c.res.maxCurvature));
    const tdpX = cases.map((c) => (c.res.tdp ? c.res.tdp.x - (c.res.config.top.x) : NaN));
    const tdpRange = cases.length > 1 && tdpX.every(Number.isFinite) ? Math.max(...tdpX) - Math.min(...tdpX) : 0;
    const strain = maxCurv * pp.D / 2;
    need(strain < 0.03, `The solution bends the pipe to ${fmt(strain * 100)} % strain, far beyond the elastic range this model assumes. The configuration is not viable: increase the hang-off angle or reduce the vessel offset.`);

    const kpis = [
      kpi('topTension', 'Top effective tension (mean position)', base.topTension / 1e3, 'kN', 'Effective tension the hang-off must carry in the mean position. It sizes the porch, the flex joint and the hull payload; it grows with water depth, contents density and hang-off angle.'),
      kpi('topTensionMax', 'Largest top tension over the offset cases', Math.max(...cases.map((c) => c.res.topTension)) / 1e3, 'kN', 'Highest hang-off tension among mean, near and far positions; the far case normally governs because the catenary is pulled flatter.'),
      kpi('hangoffAngle', 'Hang-off angle from vertical', base.hangoffAngleDeg, 'deg', 'Angle of the riser axis at the porch in the mean position. The difference from the porch (receptacle) angle is what the flex joint or stress joint has to accommodate.'),
      kpi('tdpTension', 'Effective tension at touchdown', base.Te[kT] / 1e3, 'kN', 'Horizontal tension left in the pipe where it meets the seabed. Low values mean a tight sag bend; it must stay positive in every case to avoid compression.'),
      kpi('minTe', 'Minimum effective tension (all cases)', minTe / 1e3, 'kN', 'Lowest effective tension anywhere along the riser. Negative effective tension means compression and a risk of global buckling at the touchdown or sag bend.', 0, 'min'),
      kpi('reactionH', 'Horizontal hang-off reaction', Math.abs(base.endA.Fx) / 1e3, 'kN', 'Horizontal pull of the riser on the vessel in the mean position; it enters the mooring or station-keeping load balance.'),
      kpi('reactionV', 'Vertical hang-off reaction', Math.abs(base.endA.Fz) / 1e3, 'kN', 'Vertical load of the riser on the porch: the suspended apparent weight. It is hull payload that the host must carry for the field life.'),
      kpi('suspendedLength', 'Suspended length', base.tdp ? base.tdp.s : length, 'm', 'Arc length from the hang-off to the touchdown point in the mean position.'),
      kpi('tdpX', 'Touchdown distance from hang-off', base.tdp ? base.tdp.x : inp.lazySpan ?? 0, 'm', 'Horizontal distance from the hang-off to the touchdown point; it fixes the seabed layout and the riser corridor.'),
      kpi('tdpRange', 'Touchdown movement near → far', tdpRange, 'm', 'Length of seabed swept by the touchdown point between the near and far offsets: the zone where trenching and fatigue concentrate.'),
      kpi('minBendRadius', 'Minimum bend radius', 1 / maxCurv, 'm', 'Tightest radius anywhere in the riser over all cases, normally in the sag bend near touchdown (near case) or in the hog bend of a lazy wave.'),
      kpi('maxStrain', 'Maximum bending strain', strain * 100, '%', 'Outer-fibre bending strain κ·D/2 at the tightest bend. Steel risers are normally kept elastic, well below 0.2 %.'),
      kpi('maxMoment', 'Maximum bending moment', Math.max(...cases.map((c) => c.res.maxMoment)) / 1e3, 'kN·m', 'Largest bending moment along the riser over all cases; it combines with tension and pressure in the wall-stress check.'),
      kpi('vmMax', 'Maximum von Mises stress', gov.stress.max / 1e6, 'MPa', `Largest combined stress from tension, bending and pressure along the riser (${gov.name} case), compared with the API RP 2RD allowable Cf·Ca·σy.`, allow / 1e6, 'max'),
    ];
    if (isLazy) {
      kpis.push(kpi('sagClearance', 'Sag-bend clearance above the seabed', base.sagClearance, 'm', 'Height of the lowest point of the sag bend above the seabed. Zero means no wave has formed; too little clearance lets the sag bend touch down in the near or heavy-contents case.'),
        kpi('hogElevation', 'Hog-bend elevation above the seabed', base.hogElevation, 'm', 'Height of the top of the buoyant arch. It must stay clear of the surface zone, neighbouring lines and the vessel keel.'),
        kpi('netBuoyancy', 'Net uplift of the buoyancy section', base.netBuoyancy / 1e3, 'kN', 'Net upward force provided by the modules; foam water absorption over life reduces it and lowers the arch.'));
    }
    let overlay = null;
    if (inp.centreline?.points?.length >= 6) {
      const p = inp.centreline.points, np = p.length / 3;
      let zMin = Infinity; for (let i = 0; i < np; i++) zMin = Math.min(zMin, p[3 * i + 2]);
      // project on the vertical plane through the first point, highest point first (hang-off)
      const first = p[2] >= p[3 * (np - 1) + 2] ? 0 : np - 1;
      const xs = [], zs = []; let dev = 0;
      for (let i = 0; i < np; i++) { xs.push(Math.hypot(p[3 * i] - p[3 * first], p[3 * i + 1] - p[3 * first + 1])); zs.push(p[3 * i + 2] - zMin); }
      for (let i = 0; i < np; i++) { let d = Infinity; for (let k = 0; k <= n; k++) d = Math.min(d, Math.hypot(xs[i] - base.x[k], zs[i] - base.z[k])); dev = Math.max(dev, d); }
      overlay = { xs, zs };
      kpis.push(kpi('centrelineDev', 'Largest deviation from the imported centreline', dev, 'm', 'Greatest distance between the imported centreline points and the computed mean shape. A large value points to different contents, offset, length or seabed level than assumed.'));
    }

    const checks = [
      check('Combined stress (von Mises)', gov.stress.max / 1e6, allow / 1e6, 'MPa', 'API RP 2RD (1998) §5.2.3, Table 2', 'Guards against yielding of the pipe wall under tension, bending and pressure together. Above the allowable the wall has less than the code margin against gross plastic deformation.'),
    ];
    if (inp.config === 'scr' && inp.topType === 'flex') {
      const rot = Math.max(...cases.map((c) => Math.abs(c.res.hangoffAngleDeg - base.catenary.hangoffAngleDeg)));
      checks.push(check('Flex-joint rotation range', rot, inp.flexAllow, 'deg', 'Flex-joint supplier rating (project data)', 'Rotation of the riser relative to the porch between the offset cases. Beyond the rated cone angle the elastomer is over-strained and the joint body can contact.'));
    }

    const xz = cases.map((c) => { const [x, z] = thin([c.res.x, c.res.z], 250); return { name: c.name, x, y: z, ...(c.offset ? { dash: true } : {}) }; });
    const xmax = Math.max(...cases.map((c) => c.res.x[n])), xmin = Math.min(0, ...cases.map((c) => c.res.x[0]));
    xz.push({ name: 'Seabed', x: [xmin, xmax], y: [0, 0], dash: true });
    if (overlay) xz.push({ name: 'Imported centreline', x: overlay.xs, y: overlay.zs });
    const along = (f, scale) => cases.map((c) => { const [s, y] = thin([c.res.s, f(c)], 300); return { name: c.name, x: s, y: scaled(y, scale) }; });
    const [x3, z3, te3] = thin([base.x, base.z, base.Te], 200);
    const plots = [
      plot('Static configuration', 'Horizontal distance from mean hang-off', 'm', 'Height above seabed', 'm', xz, 'Riser shape for the mean, near and far vessel positions. Look at how far the touchdown point moves and how the sag bend tightens in the near case.'),
      plot('Effective tension along the riser', 'Arc length from hang-off', 'm', 'Effective tension', 'kN', along((c) => c.res.Te, 1e-3), 'Effective tension from the hang-off (left) to the seabed end. It falls with the apparent weight carried below each point; watch for values approaching zero at the touchdown.'),
      plot('Bending moment along the riser', 'Arc length from hang-off', 'm', 'Bending moment', 'kN·m', along((c) => c.res.moment, 1e-3), 'Bending moment along the riser. Peaks appear at the hang-off (if rotationally restrained) and in the touchdown zone, where the curvature changes from catenary to straight.'),
      plot('Shear force along the riser', 'Arc length from hang-off', 'm', 'Shear force', 'kN', along((c) => c.res.shear, 1e-3), 'Transverse shear in the pipe. It is small except in the boundary layers at the hang-off and at touchdown, where the seabed reaction is picked up.'),
      plot('Curvature along the riser', 'Arc length from hang-off', 'm', 'Curvature', '1/km', along((c) => c.res.curvature, 1e3), 'Curvature (1/bend radius). The maximum sets the minimum bend radius and the bending strain; compare the near case with the mean.'),
      plot('von Mises stress along the riser', 'Arc length from hang-off', 'm', 'von Mises stress', 'MPa', [...cases.map((c) => { const [s, v] = thin([c.stress.s, c.stress.vonMises], 300); return { name: c.name, x: s, y: scaled(v, 1e-6) }; }), { name: 'Allowable', x: [0, length], y: [allow / 1e6, allow / 1e6], dash: true }], 'Combined wall stress against the allowable line. The governing location is where the curve comes closest to the dashed limit.'),
      { type: 'model3d', title: 'Riser model (coloured by effective tension)', nodes: x3.map((x, i) => [x, 0, z3[i]]), lines: x3.slice(1).map((_, i) => [i, i + 1]), nodeValues: scaled(te3, 1e-3), valueLabel: 'Effective tension', valueUnit: 'kN',
        explain: 'Three-dimensional view of the mean configuration in its vertical plane, coloured by effective tension from the hang-off to the seabed end.' },
    ];
    const tables = [
      { title: 'Stations along the riser (mean position)', columns: [{ key: 's', label: 'Arc length', unit: 'm' }, { key: 'x', label: 'x', unit: 'm' }, { key: 'z', label: 'Height', unit: 'm' }, { key: 'Te', label: 'Effective tension', unit: 'kN' }, { key: 'Tw', label: 'True wall tension', unit: 'kN' }, { key: 'M', label: 'Moment', unit: 'kN·m' }, { key: 'V', label: 'Shear', unit: 'kN' }, { key: 'kappa', label: 'Curvature', unit: '1/km' }, { key: 'vm', label: 'von Mises', unit: 'MPa' }],
        rows: stations({ s: base.s, x: base.x, z: base.z, Te: scaled(base.Te, 1e-3), Tw: scaled(base.Tw, 1e-3), M: scaled(base.moment, 1e-3), V: scaled(base.shear, 1e-3), kappa: scaled(base.curvature, 1e3), vm: scaled(cases[0].stress.vonMises, 1e-6) }),
        explain: 'Spatial coordinates and section forces at stations from the hang-off to the seabed end, for the report and for seeding local analyses.' },
      { title: 'Load cases', columns: [{ key: 'case', label: 'Case' }, { key: 'T', label: 'Top tension', unit: 'kN' }, { key: 'ang', label: 'Hang-off angle', unit: 'deg' }, { key: 'tdp', label: 'Touchdown x', unit: 'm' }, { key: 'mbr', label: 'Min bend radius', unit: 'm' }, { key: 'vm', label: 'Max von Mises', unit: 'MPa' }, { key: 'uc', label: 'Utilisation', unit: '-' }],
        rows: cases.map((c, i) => ({ case: c.name, T: round(c.res.topTension / 1e3, 5), ang: round(c.res.hangoffAngleDeg, 4), tdp: round(tdpX[i], 5), mbr: round(c.res.minBendRadius, 4), vm: round(c.stress.max / 1e6, 4), uc: round(c.stress.max / allow, 3) })),
        explain: 'One row per vessel position. Use it to see which offset governs tension (far) and which governs bending and touchdown stress (near).' },
    ];

    const steel = steelTonnes(pp.D, pp.t, length, pp.mat.rho);
    const uc = gov.stress.max / allow;
    const advice = [];
    advice.push(adv('engineering', uc > 1 ? 'fail' : uc > 0.9 ? 'warn' : 'ok',
      `Combined stress utilisation is ${fmt(uc)} in the ${gov.name} case at s = ${fmt(gov.stress.sMax)} m; minimum effective tension is ${fmt(minTe / 1e3)} kN.`,
      uc > 1 ? `Reduce bending at the governing location: ${isLazy ? 'lengthen the buoyancy section or reduce module lift per metre' : `raise the hang-off angle above ${inp.hangoffAngle}° or soften the top restraint`}; a wall of about ${fmt(pp.t * 1e3 * uc / 0.9)} mm would bring the stress to 0.9 if pressure and tension govern.`
        : 'Confirm the result with the dynamic and fatigue analyses before fixing the configuration: static stress rarely governs a deep-water riser, fatigue at the touchdown and hang-off usually does.',
      `API RP 2RD (1998, R2013) §5.2.3 working-stress format with Cf = ${inp.Cf}, Ca = ${inp.CaAllow}. ${CONFIRM}`));
    if (minTe <= 0) advice.push(adv('engineering', 'fail', `Effective compression of ${fmt(-minTe / 1e3)} kN occurs along the riser.`, 'Increase the hang-off angle or reduce the near offset so that the touchdown zone stays in tension; compression leads to global buckling of the sag bend.', 'DNV-ST-F201 (2018) Sec. 5 — effective tension to remain positive unless buckling is explicitly assessed.'));
    if (!isLazy) {
      // sensitivity of payload and length to the hang-off angle, from the closed-form catenary
      const w = pp.section.w, c0 = catenaryRiser({ h, w, hangoffAngleDeg: inp.hangoffAngle }), c1 = catenaryRiser({ h, w, hangoffAngleDeg: Math.min(89, inp.hangoffAngle + 1) });
      const dS = c1.suspendedLength - c0.suspendedLength, dT = c1.topTension - c0.topTension;
      advice.push(adv('economic', 'info',
        `The riser needs ${fmt(length)} m of pipe, ${fmt(steel)} t of steel, and loads the host with ${fmt(Math.abs(base.endA.Fz) / 1e3 / G)} t of vertical payload. One more degree of hang-off angle adds ${fmt(dS)} m of suspended pipe (${fmt(steelTonnes(pp.D, pp.t, dS, pp.mat.rho))} t) and ${fmt(dT / 1e3)} kN of top tension.`,
        'Fix the hang-off angle at the smallest value that passes touchdown fatigue and near-offset compression: every degree saved reduces pipe, welds, installation time and hull payload.',
        'Whole-life costing (ISO 15686-5) with quantities from the configuration; estimate class 4 per AACE RP 18R-97 until the dynamic design is closed.',
        `Steel in the riser: ${steelCo2Text(steel)}`));
    } else {
      const foam = Math.PI / 4 * ((inp.buoyD / 1e3) ** 2 - (pp.D + 2 * inp.coatT / 1e3) ** 2) * inp.buoyFill * inp.buoyLength * inp.buoyRho / 1e3;
      advice.push(adv('economic', 'info',
        `The lazy wave uses ${fmt(length)} m of pipe (${fmt(steel)} t of steel) and ${fmt(foam)} t of buoyancy foam giving ${fmt(m.base.netBuoyancy / 1e3)} kN of net uplift; top tension is ${fmt(base.topTension / 1e3)} kN.`,
        'Trade buoyancy length against sag-bend clearance: remove modules until the clearance in the heaviest-contents, near-offset case is just acceptable. Modules and their clamps are bought, installed and inspected one by one.',
        'Whole-life costing (ISO 15686-5); quantities are the basis for an AACE class 4 estimate.',
        `Steel: ${steelCo2Text(steel)} Syntactic foam is a thermoset composite with no established recycling route, so each tonne avoided also reduces end-of-life waste.`));
    }
    advice.push(adv('engineering', 'info', 'The model is two-dimensional: only in-plane current and offsets are represented.', 'Run the cross-current and cross-offset cases, and any riser-array interference, in a 3-D model before detailed design.', 'API RP 2RD §5.4 (near, far and cross load cases).'));

    const quantities = { steel_t: steel, pipe_km: length / 1e3 };
    if (isLazy) quantities.polymer_t = Math.PI / 4 * ((inp.buoyD / 1e3) ** 2 - (pp.D + 2 * inp.coatT / 1e3) ** 2) * inp.buoyFill * inp.buoyLength * inp.buoyRho / 1e3;
    return finish({
      kpis, checks, plots, tables,
      notes: [
        `${isLazy ? 'Lazy-wave' : 'Steel catenary'} riser solved with ${n} beam elements; residual force imbalance ${fmt(base.residual)} N.`,
        `Apparent weight ${fmt(pp.section.w)} N/m; governing stress in the ${gov.name} case. ${current ? `In-plane current ${inp.Us} m/s at the surface, ${inp.curDir === 'far' ? 'towards the touchdown side' : 'towards the vessel side'}.` : 'No current applied.'}`,
        `Embodied carbon factor used in advice: ${CO2_STEEL} t CO₂e per tonne of steel (screening).`,
      ],
      advice,
      publish: {
        'riser.topTension': base.topTension, 'riser.hangoffAngleDeg': base.hangoffAngleDeg, 'riser.tdpTension': base.Te[kT], 'riser.reactionH': Math.abs(base.endA.Fx), 'riser.reactionV': Math.abs(base.endA.Fz),
        'riser.suspendedLength': base.tdp ? base.tdp.s : length, 'riser.length': length, 'riser.tdpCurvature': base.maxCurvature, 'riser.maxMoment': base.maxMoment,
      },
      quantities,
    });
  },

  verify: [
    {
      // Inextensible catenary: a = h·sinθ/(1 − sinθ), top tension T = w·(h + a) = w·h/(1 − sinθ).
      // Defaults: h = 1000 − 15 = 985 m, θ = 12°, w from the hand weight formula ≈ 955.9 N/m → T ≈ 1188.7 kN.
      name: 'Free-hanging catenary top tension T = w·h/(1 − sin θ)', source: 'Irvine, Cable Structures (1981) §1.2; Sparks, Fundamentals of Marine Riser Mechanics (2007) catenary equations',
      inputs: { config: 'scr', topType: 'pinned', Us: 0, Ub: 0, offset: 0, pInt: 0, nElem: 320 },
      expect: [{ kpi: 'topTension', value: handWeight() * 985 / (1 - Math.sin(12 * DEG)) / 1e3, tol: 0.005 },
        // vertical equilibrium: porch carries the suspended apparent weight w·S, S = sqrt(h(h + 2a))
        { kpi: 'reactionV', value: handWeight() * Math.sqrt(985 * (985 + 2 * 985 * Math.sin(12 * DEG) / (1 - Math.sin(12 * DEG)))) / 1e3, tol: 0.01 }],
    },
  ],
};
