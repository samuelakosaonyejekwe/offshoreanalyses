// End-fitting and termination loads: reactions transferred through the end fitting (axial, shear, moment, torque,
// pressure end cap), bend-stiffener response at the termination, armour-wire anchorage capacity and the connector
// rating check.
import { axisymmetricResponse, bendingProperties, minimumBendRadius } from '../../engine/xsection/flexible.js';
import { bendStiffener } from '../../engine/xsection/ancillaries.js';
import { endFittingAnchorage } from '../../engine/xsection/flexfatigue.js';
import { stackInputs, stackFrom, loadInputs, loadsFrom, muInput, kpi, chk, adv, plot, bar, tbl, col, run, need, fmt, round, thin, sev, worst, LIMIT_AXI, GENERIC, STD_17J, STD_17B, CONFIRM, DEG, BAR, CO2_MGO, MGO_T_PER_DAY } from './_common.js';

export default {
  id: 's14-end-fitting',
  suite: 14,
  title: 'End-fitting and termination loads',
  summary: 'Resolves the loads an end fitting must transfer, the curvature inside the bend stiffener, the armour anchorage capacity and the connector utilisation.',
  covers: [13],
  fidelity: 'numerical',
  standards: ['API Spec 17J (4th ed.) §5.4 (end fittings)', 'API RP 17B (5th ed.) §8', 'API Spec 17L1 / RP 17L2 (ancillary equipment: bend stiffeners)'],
  models: ['Elastica of a tapered polyurethane bend stiffener (Boef & Out 1990)', 'Armour-wire anchorage: bond and cone-wedge statics', 'Pressure end-cap reaction and true wall force', 'Linear tension–moment connector interaction'],
  limits: LIMIT_AXI + ' The bend stiffener is a planar elastic beam with a linear taper and a curvature-independent pipe stiffness; polyurethane non-linearity and gaps are not modelled. Anchorage is uniform bond plus wedge statics with test-derived coefficients; the epoxy stress field, crimped seals and the end-fitting body are not analysed (no local finite-element model). ' + GENERIC,
  publishes: ['flexible.endFittingMoment', 'flexible.endFittingShear'],
  searchTerms: 'flexible riser end fitting bend stiffener armour anchorage epoxy',
  inputs: [
    ...stackInputs(),
    ...loadInputs({ Te: 800, depth: 0, torque: true }).map((i) => (i.key === 'Te' ? { ...i, label: 'Effective tension at the termination', min: 1 } : i)),
    { key: 'angle', label: 'Angle between tension and end-fitting axis', unit: 'deg', type: 'number', default: 12, min: 0.05, max: 60, group: 'Loads and pressure boundaries', role: 'bc', design: true,
      uncertain: { dist: 'normal', cov: 0.15 }, help: 'Largest deviation of the riser tangent from the end-fitting axis (vessel motion plus static angle).' },
    muInput(),
    { key: 'bsLength', label: 'Bend stiffener length', unit: 'm', type: 'number', default: 4, min: 0.3, max: 15, group: 'Bend stiffener', role: 'geometry', design: true },
    { key: 'bsRootD', label: 'Stiffener root outer diameter', unit: 'mm', type: 'number', default: 700, min: 30, max: 3000, group: 'Bend stiffener', role: 'geometry', design: true, help: 'Set root and tip diameter equal to the pipe diameter to analyse a termination without a stiffener.' },
    { key: 'bsTipD', label: 'Stiffener tip outer diameter', unit: 'mm', type: 'number', default: 280, min: 30, max: 3000, group: 'Bend stiffener', role: 'geometry' },
    { key: 'matPu', label: 'Stiffener material', type: 'select', source: 'materials:polymer', default: 'pu-stiffener', group: 'Bend stiffener', role: 'material', help: 'Supplies the secant modulus of the polyurethane; it is strain-, rate- and temperature-dependent.' },
    { key: 'puStrainLimit', label: 'Allowable polyurethane strain', unit: '%', type: 'number', default: 10, min: 1, max: 30, group: 'Acceptance criteria', role: 'criteria', help: 'Supplier design value for the load case.' },
    { key: 'EIoverride', label: 'Pipe bending stiffness (0 = slip stiffness from the layer stack)', unit: 'kN·m²', type: 'number', default: 0, min: 0, max: 1e5, group: 'Bend stiffener', role: 'geometry', help: 'Data-sheet bending stiffness of the pipe at the operating state.' },
    { key: 'bondStress', label: 'Epoxy bond stress on the wire', unit: 'MPa', type: 'number', default: 8, min: 0, max: 60, group: 'Armour anchorage', role: 'material', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Average pull-out bond stress from anchorage tests on the actual wire and resin.' },
    { key: 'embedLength', label: 'Embedded wire length', unit: 'mm', type: 'number', default: 350, min: 20, max: 2000, group: 'Armour anchorage', role: 'geometry', design: true },
    { key: 'coneAngle', label: 'Cone half-angle of the epoxy plug (0 = no wedge action)', unit: 'deg', type: 'number', default: 8, min: 0, max: 45, group: 'Armour anchorage', role: 'geometry' },
    { key: 'muCone', label: 'Friction epoxy / end-fitting body', unit: '-', type: 'number', default: 0.1, min: 0, max: 0.6, group: 'Armour anchorage', role: 'material' },
    { key: 'muWire', label: 'Friction wire / epoxy', unit: '-', type: 'number', default: 0.1, min: 0, max: 0.6, group: 'Armour anchorage', role: 'material' },
    { key: 'sfAnchor', label: 'Safety factor on anchorage capacity', unit: '-', type: 'number', default: 1.5, min: 1, max: 4, group: 'Acceptance criteria', role: 'criteria', help: 'Project value; the anchorage must in any case be stronger than the wire (efficiency ≥ 1).' },
    { key: 'ratedTension', label: 'Connector rated tension', unit: 'kN', type: 'number', default: 3000, min: 10, max: 50000, group: 'Connector rating', role: 'criteria', help: 'Mechanical load rating of the flange or hub with zero moment.' },
    { key: 'ratedMoment', label: 'Connector rated bending moment', unit: 'kN·m', type: 'number', default: 600, min: 1, max: 50000, group: 'Connector rating', role: 'criteria' },
    { key: 'ratedPressure', label: 'End-fitting pressure rating', unit: 'bar', type: 'number', default: 345, min: 1, max: 2000, group: 'Connector rating', role: 'criteria' },
    { key: 'nSteps', label: 'Integration steps along the stiffener', type: 'int', default: 200, min: 20, max: 3200, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nSteps', label: 'Integration steps', levels: [50, 100, 200, 400], qoi: 'maxCurvature' },
  compute(inp, ctx) {
    const S = stackFrom(inp, ctx);
    need(S.idx.ta.length > 0, 'The stack has no tensile armour layer to anchor.');
    const loads = loadsFrom(inp);
    const r = run(() => axisymmetricResponse(S.stack, loads), 'Axisymmetric solution');
    const bend = bendingProperties(S.stack, { response: r, mu: inp.mu });
    const EIpipe = inp.EIoverride > 0 ? inp.EIoverride * 1e3 : bend.EIslip;
    const od = S.od, root = inp.bsRootD / 1e3, tip = inp.bsTipD / 1e3;
    need(tip >= od - 1e-9, `Stiffener tip diameter ${fmt(inp.bsTipD)} mm is smaller than the pipe outer diameter ${fmt(od * 1e3)} mm.`);
    need(root >= tip, 'Stiffener root diameter must not be smaller than the tip diameter.');
    const pu = ctx.material(inp.matPu);
    need(pu && pu.E > 0, `Material '${inp.matPu}' has no modulus: pick a stiffener polyurethane.`);
    let mbrDyn = null;
    try { mbrDyn = minimumBendRadius(S.stack).dynamic; } catch { mbrDyn = null; }
    const bs = run(() => bendStiffener({ tension: loads.tension, angle: inp.angle * DEG, length: inp.bsLength, rootDiameter: Math.max(root, tip, od), tipDiameter: Math.max(tip, od), innerDiameter: od, Epu: pu.E, EIpipe, nSteps: inp.nSteps,
      mbr: mbrDyn ?? undefined, strainLimit: inp.puStrainLimit / 100 }), 'Bend stiffener');
    const anch = S.idx.ta.map((i) => {
      const L = S.stack.layers[i];
      const t = L.wireThickness ?? L.t, b = L.wireWidth ?? 3 * t;
      const n = L.n ?? Math.max(1, Math.round(L.fill * 2 * Math.PI * L.R * Math.cos(L.alpha) / b));
      const A = L.nA / n;
      need(L.uts > 0, `${L.name}: enter the wire tensile strength in the layer table.`);
      const a = run(() => endFittingAnchorage({ bondStress: inp.bondStress * 1e6, perimeter: 2 * (b + t), embeddedLength: inp.embedLength / 1e3, wireArea: A, uts: L.uts, nWires: n,
        coneHalfAngle: inp.coneAngle > 0 ? inp.coneAngle * DEG : undefined, muCone: inp.muCone, muWire: inp.muWire }), `${L.name} anchorage`);
      // friction from wedging can exceed the pull (self-locking plug): the wire then governs, cap the anchorage at a large multiple
      if (!Number.isFinite(a.capacity)) { a.capacity = 99 * a.wireStrength; a.efficiency = 99; }
      return { L, n, A, a, force: Math.abs(r.layers[i].wireStress) * A };
    });
    const T = loads.tension, M = bs.rootMoment;
    const ucConn = T / (inp.ratedTension * 1e3) + M / (inp.ratedMoment * 1e3);
    const checks = [
      ...anch.map((y) => chk(`${y.L.name}: wire force vs anchorage capacity`, y.force / 1e3, y.a.capacity / inp.sfAnchor / 1e3, 'kN', `${STD_17J} §5.4 (end-fitting anchorage); bond and wedge statics, safety factor ${inp.sfAnchor}`,
        'Guards against the armour wires pulling out of the epoxy. Pull-out releases the tensile armour and the pipe parts at the end fitting.')),
      ...anch.map((y) => chk(`${y.L.name}: anchorage stronger than the wire`, y.a.wireStrength / 1e3, y.a.capacity / 1e3, 'kN', `${STD_17J} §5.4 (termination shall not be the weak link)`,
        'Checks that the wire, not its anchorage, would fail first: the anchorage capacity must exceed the wire breaking load.')),
      chk('Connector: tension + moment interaction', ucConn, 1, '-', 'Manufacturer rating envelope (linear T/T_r + M/M_r); API Spec 17J §5.4',
        'Guards against overloading the flange or hub connection by combined tension and bending; leakage or bolt overstress follows.'),
      chk('End-fitting pressure rating', inp.pInt, inp.ratedPressure, 'bar', `${STD_17J} §5.4`, 'The bore pressure must not exceed the pressure rating of the end fitting and its seals.'),
      chk('Polyurethane strain in the bend stiffener', bs.maxPuStrain * 100, inp.puStrainLimit, '%', 'API Spec 17L1 (bend stiffener design); supplier allowable strain',
        'Guards against tearing or fatigue of the stiffener cone at its most strained fibre.'),
    ];
    if (mbrDyn) checks.push(chk('Pipe curvature at the termination vs dynamic MBR', bs.maxCurvature, 1 / mbrDyn, '1/m', `${STD_17J} §5.3.1 (dynamic MBR) — ${CONFIRM}`,
      'Guards against over-bending the pipe where it leaves the end fitting, the classic location of armour fatigue and sheath failure.'));
    const gov = worst(checks);
    const [ss, kk, mm] = thin([bs.s, bs.curvature, bs.moment], 200);
    const minEff = Math.min(...anch.map((y) => y.a.efficiency));
    const needLen = Math.max(...anch.map((y) => inp.embedLength * Math.max(y.force * inp.sfAnchor, y.a.wireStrength) / y.a.capacity));
    const puVol = Math.PI / 12 * inp.bsLength * (root * root + root * tip + tip * tip) - Math.PI / 4 * od * od * inp.bsLength; // frustum minus bore, m³
    return {
      kpis: [
        kpi('endCapForce', 'Pressure end-cap force', r.endCapForce / 1e3, 'kN', 'Net pressure force on the sealing layers that the end fitting must react in addition to the effective tension.'),
        kpi('wallForce', 'Axial force through the armour anchorage', r.wallForce / 1e3, 'kN', 'True wall force transferred from the layers into the end-fitting body (effective tension plus end cap).'),
        kpi('rootAxial', 'Axial reaction at the end fitting', bs.rootAxial / 1e3, 'kN', 'Component of the effective tension along the end-fitting axis.'),
        kpi('rootShear', 'Shear reaction at the end fitting', bs.rootShear / 1e3, 'kN', 'Component of the tension across the end-fitting axis; it is carried through the stiffener interface structure.'),
        kpi('rootMoment', 'Bending moment at the end fitting', M / 1e3, 'kN·m', 'Total moment reacted at the stiffener root; almost all of it goes through the stiffener cone into the support, not through the pipe.'),
        kpi('torque', 'Torque at the end fitting', r.torque / 1e3, 'kN·m', 'Torque transferred through the armour anchorage.'),
        kpi('maxCurvature', 'Largest pipe curvature at the termination', bs.maxCurvature, '1/m', `Peak curvature of the pipe along the stiffener, at ${fmt(bs.maxCurvatureAt)} m from the root. It drives armour fatigue at the hang-off.`),
        kpi('curvPerAngle', 'Curvature per degree of angle', bs.curvaturePerAngle * DEG, '1/m per deg', 'Transfer function from end angle to peak curvature, for converting angle ranges from global analysis into curvature ranges for fatigue.'),
        kpi('puStrain', 'Peak polyurethane strain', bs.maxPuStrain * 100, '%', 'Outer-fibre strain of the stiffener cone.', inp.puStrainLimit, 'max'),
        kpi('wireForce', 'Largest armour wire force', Math.max(...anch.map((y) => y.force)) / 1e3, 'kN', 'Force in the most loaded tensile armour wire, to be anchored in the epoxy.'),
        kpi('anchorEfficiency', 'Anchorage efficiency', minEff, '-', 'Anchorage capacity over wire breaking load. At or above 1.0 the wire breaks before it pulls out, as required.', 1, 'min'),
        kpi('connectorUc', 'Connector utilisation', ucConn, '-', 'Tension over rated tension plus moment over rated moment. Above 1.0 the connection is outside its rating envelope.', 1, 'max'),
        kpi('maxUc', 'Governing utilisation', gov.uc, '-', `Highest utilisation of all termination checks (${gov.name}).`, 1, 'max'),
      ],
      checks,
      plots: [
        plot('line', 'Pipe curvature along the bend stiffener', 'Distance from root', 'm', 'Curvature', '1/m', [{ name: 'Curvature', x: ss.map((v) => round(v)), y: kk.map((v) => round(v)) },
          ...(mbrDyn ? [{ name: '1 / dynamic MBR', x: [0, round(inp.bsLength)], y: [round(1 / mbrDyn), round(1 / mbrDyn)], dash: true }] : [])],
        'Curvature of the pipe inside the stiffener. A well-proportioned stiffener gives a flat distribution; a peak at the tip means it is too short or too stiff, a peak at the root that it is too soft.'),
        plot('line', 'Bending moment along the bend stiffener', 'Distance from root', 'm', 'Moment', 'kN·m', [{ name: 'Moment', x: ss.map((v) => round(v)), y: mm.map((v) => round(v / 1e3)) }],
          'Total bending moment in stiffener plus pipe. The root value is the moment the support structure must be designed for.'),
        bar('Armour wire force and anchorage capacity', anch.map((y) => y.L.name), 'Force per wire', 'kN', [
          { name: 'Wire force', y: anch.map((y) => round(y.force / 1e3)) }, { name: 'Factored anchorage capacity', y: anch.map((y) => round(y.a.capacity / inp.sfAnchor / 1e3)) }, { name: 'Wire breaking load', y: anch.map((y) => round(y.a.wireStrength / 1e3)) }],
        'Per wire: the force to be anchored, the factored pull-out capacity and the wire breaking load. Capacity must exceed both the demand and the breaking load.'),
      ],
      tables: [tbl('Armour anchorage', [col('layer', 'Layer'), col('n', 'Wires'), col('F', 'Wire force', 'kN'), col('bond', 'Bond capacity', 'kN'), col('cap', 'Capacity with wedge', 'kN'), col('str', 'Wire breaking load', 'kN'), col('eff', 'Efficiency', '-'), col('rel', 'Wedge self-releasing')],
        anch.map((y) => ({ layer: y.L.name, n: y.n, F: round(y.force / 1e3), bond: round(y.a.bondCapacity / 1e3), cap: round(y.a.capacity / 1e3), str: round(y.a.wireStrength / 1e3), eff: round(y.a.efficiency, 3), rel: y.a.selfReleasing === null ? 'no wedge' : (y.a.selfReleasing ? 'yes' : 'no') })),
        'One row per tensile armour layer: wire count, force per wire and the pull-out capacity from bond and cone wedging.')],
      notes: [`${gov.name} governs (${fmt(gov.uc)}).`, `Reactions: axial ${fmt(bs.rootAxial / 1e3)} kN, shear ${fmt(bs.rootShear / 1e3)} kN, moment ${fmt(M / 1e3)} kN·m, end cap ${fmt(r.endCapForce / 1e3)} kN.`,
        `Pipe bending stiffness used: ${fmt(EIpipe / 1e3)} kN·m² (${inp.EIoverride > 0 ? 'entered' : 'slip stiffness of the stack'}).`, GENERIC],
      advice: [
        adv('engineering', sev(gov.uc), `Governing termination utilisation ${fmt(gov.uc)} (${gov.name}); peak pipe curvature ${fmt(bs.maxCurvature)} 1/m (radius ${fmt(1 / Math.max(bs.maxCurvature, 1e-9))} m) at ${fmt(inp.angle)}°.`,
          gov.uc > 1 ? (minEff < 1 || anch.some((y) => y.force * inp.sfAnchor > y.a.capacity) ? `Increase the embedded wire length to about ${fmt(needLen)} mm, and/or re-proportion the stiffener.` : 'Lengthen the stiffener or enlarge its root diameter to spread the curvature, or raise the connector rating.') : 'Use the curvature-per-angle transfer function for the armour fatigue load bins, and confirm bond stress and wedge friction from anchorage pull-out tests.',
          `${STD_17J} §5.4; API Spec 17L1 (2021); Boef & Out (1990) OTC 6436`),
        adv('economic', 'info', `The stiffener cone is about ${fmt(puVol)} m³ of polyurethane (${fmt(puVol * (pu.rho ?? 1150) / 1e3)} t).`,
          bs.strainUtilisation !== null && bs.strainUtilisation < 0.5 && gov.uc < 0.7 ? `The cone is only ${fmt(bs.strainUtilisation * 100)} % strained: a smaller root diameter would cut polyurethane volume and casting cost — run the design sweep on root diameter and length.` : 'The stiffener is working: a re-termination offshore to cure a fatigue problem here costs several vessel days, so the casting is cheap insurance.',
          'Whole-life costing per ISO 15686-5; quantity-based estimate at AACE Class 4',
          `Polyurethane is not readily recycled: ${fmt(puVol * (pu.rho ?? 1150) / 1e3)} t per stiffener ends as waste at decommissioning. A re-termination campaign of about 5 vessel days emits roughly ${fmt(5 * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening: ${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂ per t).`),
      ],
      publish: { 'flexible.endFittingMoment': M, 'flexible.endFittingShear': bs.rootShear },
      quantities: { polymer_t: puVol * (pu.rho ?? 1150) / 1e3 },
    };
  },
  verify: [
    {
      // No stiffener (root = tip = pipe diameter): a uniform beam pulled at angle θ. First integral of the elastica,
      // EI·κ²/2 = T(1 − cos θ) → root curvature κ = 2·√(T/EI)·sin(θ/2). T = 100 kN, EI = 50 kN·m², θ = 2°:
      // 2·√2·sin(1°) = 2.82843·0.0174524 = 0.04936 1/m. Shear = T·sin θ = 3.490 kN.
      name: 'Uniform pipe pulled at an angle: root curvature from the elastica first integral',
      source: 'Love, Mathematical Theory of Elasticity, Ch. XIX (elastica); Boef & Out (1990)',
      inputs: { Te: 100, angle: 2, EIoverride: 50, bsRootD: 236.4, bsTipD: 236.4, bsLength: 6, pInt: 0, pAnn: 0, nSteps: 400 },
      expect: [{ kpi: 'maxCurvature', value: 0.04936, tol: 0.005 }, { kpi: 'rootShear', value: 3.490, tol: 0.002 }],
    },
    {
      // End-cap force with atmospheric outside and an empty annulus: p·π·r² on the inner face of the pressure sheath,
      // r = 83.2 mm: 20e6·π·0.0832² = 434.94 kN.
      name: 'Pressure end-cap force on the pressure sheath',
      source: 'Statics, F = p·A; Sparks (2007) Ch. 2',
      inputs: { pInt: 200, depth: 0, pAnn: 0 },
      expect: [{ kpi: 'endCapForce', value: 434.94, tol: 0.001 }],
    },
  ],
};
