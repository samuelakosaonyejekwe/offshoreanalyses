// Mechanical crushing and external pressure: hydrostatic pressure on the cable and its rigid accessories, stability of
// the metallic sheath, radial squeeze of the core by tensioned armour, and the crushing line loads of handling and
// protection — tensioner pads, chute sidewall pressure, storage stacking, rock overburden and falling-rock impact.
import { cableHandlingLimits, cableImpactClass } from '../../engine/xsection/cablemech.js';
import { helixContactPressure, wireStrain } from '../../engine/xsection/helix.js';
import { solidRingStiffness, ringElasticCollapse, ovalisedRingCollapse } from '../../engine/xsection/flexcollapse.js';
import { collapsePressure } from '../../engine/pipeline/strength.js';
import { droppedObjectImpact } from '../../engine/pipeline/crossing.js';
import { cableInputs, buildMech, adv, check, worst, sev, fmt, round, run, need, linspace, G, STD_CABLE, LIMIT_GENERIC, CONFIRM, CO2_STEEL, CO2_ROCK } from './_common.js';

const NU_SHEATH = 0.42;       // Poisson's ratio of lead alloy (0.40–0.44)
const DEFAULT_BANDS = [
  { label: 'No damage expected', maxEnergy: 500 },
  { label: 'Serving / armour damage, core intact', maxEnergy: 2000 },
  { label: 'Core damage possible — inspect and test', maxEnergy: 5000 },
];

export default {
  id: 's13-crushing-external-pressure', suite: 13, title: 'Mechanical crushing and external pressure',
  summary: 'Hydrostatic pressure integrity of the cable and its joint housing, and the crushing line loads from tensioner, chute, storage stack, rock cover and falling rock against the tested crush resistance.',
  covers: [11], fidelity: 'analytical',
  standards: ['CIGRE TB 623 (2015) §5 (crush, sidewall-pressure and impact tests)', 'CIGRE TB 490 (2012) (water-pressure withstand of submarine cables and joints)', `${STD_CABLE}`, 'DNV-ST-F101 (2021) §5.4.4 (external-pressure collapse, used for the rigid joint housing)', 'DNV-RP-F107 (2019) §5 (falling-object impact energy)', 'Timoshenko & Gere, Theory of Elastic Stability ch. 7 (ring under external pressure)'],
  models: ['Hydrostatic pressure p = ρ·g·d', 'Elastic ring buckling 3EI/R³ and first-yield collapse of an oval ring (metallic sheath)', 'Bulk compression of the extruded core', 'Armour contact pressure from the helical-tendon model', 'Characteristic collapse pressure of a steel cylinder (joint housing)', 'Tensioner grip and squeeze statics, sidewall pressure T/R', 'Stacking line load in a hexagonally packed coil', 'Overburden line load under a rock berm', 'Terminal velocity and impact energy of a falling rock', 'Chord-width indentation estimate of the flat-plate crush load'],
  limits: `${LIMIT_GENERIC} Crush, sidewall-pressure and impact capacities are type-test values: the defaults here are indicative and overridable, and the indentation estimate is a screening figure, not a substitute for the test. Local 3-D contact between cores, fillers and armour wires, creep of the insulation under sustained squeeze and water ingress along the conductor after damage need finite-element analysis or tests and are outside the module.`,
  publishes: ['cable.crushLimit', 'cable.waterPressure'],
  searchTerms: 'submarine cable crush test sidewall pressure tensioner squeeze hydrostatic pressure joint housing',
  inputs: [
    ...cableInputs({ mech: true }),
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'External pressure', role: 'environment', bind: 'site.rhoW' },
    { key: 'depth', label: 'Maximum water depth on the route', unit: 'm', type: 'number', default: 150, min: 1, max: 3500, group: 'External pressure', role: 'bc', bind: 'site.depth', design: true, help: 'Pressure boundary condition: hydrostatic head at the deepest point, including tide and surge.' },
    { key: 'support', label: 'Support of the metallic sheath', type: 'select', default: 'solid', group: 'External pressure', role: 'ic', options: [{ value: 'solid', label: 'Extruded core fills the sheath (normal construction)' }, { value: 'void', label: 'Sheath unsupported over a gap (joint sleeve, corrugated or loose sheath)' }],
      help: 'Initial condition of the construction. A sheath extruded tight on the core cannot buckle inwards; an unsupported sheath is checked as a free ring.' },
    { key: 'sheathYield', label: 'Yield stress of the sheath metal', unit: 'MPa', type: 'number', default: 12, min: 3, max: 400, group: 'External pressure', role: 'material', uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Lead alloy E about 10–15 MPa (creeping material); aluminium 60–100; copper 70–200.' },
    { key: 'ovality', label: 'Sheath ovality (Dmax − Dmin)/(Dmax + Dmin)', unit: '%', type: 'number', default: 0.5, min: 0.05, max: 5, group: 'External pressure', role: 'ic', showIf: { key: 'support', in: ['void'] } },
    { key: 'bulkK', label: 'Bulk modulus of the insulation', unit: 'GPa', type: 'number', default: 2.5, min: 0.5, max: 8, group: 'External pressure', role: 'material', help: 'Polyethylene 2–3 GPa; EPR about 2 GPa.' },
    { key: 'jointD', label: 'Rigid joint housing outside diameter', unit: 'mm', type: 'number', default: 400, min: 0, max: 1500, group: 'Joint housing', role: 'geometry', help: '0 = no rigid housing on the route (flexible factory joints only).' },
    { key: 'jointT', label: 'Housing wall thickness', unit: 'mm', type: 'number', default: 12, min: 2, max: 80, group: 'Joint housing', role: 'geometry', design: true },
    { key: 'jointFy', label: 'Housing yield strength', unit: 'MPa', type: 'number', default: 355, min: 150, max: 900, group: 'Joint housing', role: 'material' },
    { key: 'jointOval', label: 'Housing ovality (Dmax − Dmin)/D', unit: '%', type: 'number', default: 1.0, min: 0.5, max: 5, group: 'Joint housing', role: 'ic' },
    { key: 'jointL', label: 'Housing length', unit: 'm', type: 'number', default: 3, min: 0.3, max: 15, group: 'Joint housing', role: 'geometry' },
    { key: 'nJoints', label: 'Number of rigid joints', type: 'int', default: 2, min: 0, max: 200, group: 'Joint housing', role: 'geometry' },
    { key: 'tension', label: 'Cable tension in the tensioner and over the chute', unit: 'kN', type: 'number', default: 90, min: 0, max: 5000, group: 'Handling loads', role: 'bc', bind: 'cable.layTension', design: true, uncertain: { dist: 'lognormal', cov: 0.1 } },
    { key: 'squeeze', label: 'Tensioner squeeze per track', unit: 'kN/m', type: 'number', default: 40, min: 1, max: 500, group: 'Handling loads', role: 'bc', design: true, help: 'Radial pad force per metre of each track.' },
    { key: 'nTracks', label: 'Tensioner tracks', type: 'select', default: '2', group: 'Handling loads', role: 'operation', options: [{ value: '2', label: 'Two tracks' }, { value: '3', label: 'Three tracks' }, { value: '4', label: 'Four tracks' }] },
    { key: 'padLength', label: 'Tensioner contact length', unit: 'm', type: 'number', default: 8, min: 0.5, max: 30, group: 'Handling loads', role: 'operation' },
    { key: 'muPad', label: 'Pad / serving friction coefficient', unit: '-', type: 'number', default: 0.3, min: 0.05, max: 0.8, group: 'Handling loads', role: 'material', uncertain: { dist: 'lognormal', cov: 0.25 }, help: 'Wet polypropylene-yarn serving on polyurethane pads: 0.2–0.4. Bitumen-coated serving can be lower.' },
    { key: 'gripSf', label: 'Safety factor on grip', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Handling loads', role: 'criteria' },
    { key: 'chuteR', label: 'Chute / sheave radius', unit: 'm', type: 'number', default: 4, min: 0.5, max: 30, group: 'Handling loads', role: 'bc', design: true },
    { key: 'stackH', label: 'Storage stack height (carousel or tank)', unit: 'm', type: 'number', default: 4, min: 0, max: 15, group: 'Handling loads', role: 'operation', help: 'Height of cable coiled above the bottom layer.' },
    { key: 'rockCover', label: 'Rock cover above the cable', unit: 'm', type: 'number', default: 1.0, min: 0, max: 5, group: 'Rock protection', role: 'bc', help: '0 = no rock placement.' },
    { key: 'rockBulk', label: 'Rock bulk density', unit: 'kg/m³', type: 'number', default: 1600, min: 1200, max: 2200, group: 'Rock protection', role: 'material' },
    { key: 'rockD', label: 'Largest rock size falling on the cable', unit: 'm', type: 'number', default: 0.25, min: 0.02, max: 1.2, group: 'Rock protection', role: 'load', help: 'Upper grading size (about D90–D100) of the armour or filter layer placed directly on the cable.' },
    { key: 'rockRho', label: 'Rock solid density', unit: 'kg/m³', type: 'number', default: 2650, min: 2200, max: 3200, group: 'Rock protection', role: 'material' },
    { key: 'rockCd', label: 'Rock drag coefficient', unit: '-', type: 'number', default: 1.0, min: 0.4, max: 2, group: 'Rock protection', role: 'load', help: 'DNV-RP-F107: 0.6–2.0 for spherical to complex shapes.' },
    { key: 'crushLimit', label: 'Allowable crush line load (flat plates)', unit: 'kN/m', type: 'number', default: 60, min: 2, max: 1000, group: 'Acceptance criteria', role: 'criteria', help: `From the crush test of the cable (load on a plate of stated length without damage). ${CONFIRM}` },
    { key: 'swpLimit', label: 'Allowable sidewall pressure', unit: 'kN/m', type: 'number', default: 50, min: 2, max: 500, group: 'Acceptance criteria', role: 'criteria', help: `From the sidewall-pressure / tensile-bending test. ${CONFIRM}` },
    { key: 'squeezeLimit', label: 'Allowable tensioner squeeze per track', unit: 'kN/m', type: 'number', default: 60, min: 2, max: 1000, group: 'Acceptance criteria', role: 'criteria', help: `From the manufacturer's caterpillar squeeze test. ${CONFIRM}` },
    { key: 'impactRated', label: 'Rated impact energy', unit: 'J', type: 'number', default: 500, min: 10, max: 50000, group: 'Acceptance criteria', role: 'criteria', help: `Impact-test energy without loss of function. ${CONFIRM}` },
    { key: 'bands', label: 'Impact classification bands', type: 'table', group: 'Acceptance criteria', role: 'criteria', columns: [{ key: 'label', label: 'Consequence', type: 'text' }, { key: 'maxEnergy', label: 'Up to', unit: 'J', type: 'number' }], default: DEFAULT_BANDS.map((r) => ({ ...r })), help: 'Ascending energy bands from the impact type test or a dropped-object study (project data).' },
    { key: 'sigmaIns', label: 'Allowable compressive stress on the insulation', unit: 'MPa', type: 'number', default: 8, min: 1, max: 40, group: 'Acceptance criteria', role: 'criteria', help: `Sustained radial stress the extruded insulation takes without measurable indentation; of the order of half the compressive yield. ${CONFIRM}` },
    { key: 'indent', label: 'Permissible flattening in the crush estimate', unit: '% of diameter', type: 'number', default: 0.5, min: 0.1, max: 10, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'gammaP', label: 'Safety factor on external-pressure collapse', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Acceptance criteria', role: 'criteria', help: `DNV-ST-F101 uses γm·γSC = 1.15 × 1.04–1.26 on collapse; 1.5 is a common value for qualification by calculation. ${CONFIRM}` },
  ],
  compute(inp, ctx) {
    const b = buildMech(inp, ctx, { rhoW: inp.rhoW }), cab = b.cab, OD = cab.outerDiameter;
    const p = inp.rhoW * G * inp.depth;                                   // Pa
    const T = inp.tension * 1e3;
    // metallic sheath as a free ring (plane strain): lower bound that applies only where the core does not support it
    const sh = cab.layers.find((L) => L.type === 'sheath');
    let ring = null;
    if (sh) {
      const t = sh.ro - sh.ri, R = 0.5 * (sh.ro + sh.ri), E = b.m.sh.E ?? 16e9;
      const pCr = ringElasticCollapse({ EI: solidRingStiffness({ E, thickness: t, nu: NU_SHEATH }).EI, R }).pCr;
      const ov = ovalisedRingCollapse({ pCr, yieldStress: inp.sheathYield * 1e6, R, thickness: t, ovality: inp.ovality / 100 });
      ring = { t, R, pCr, pCollapse: ov.pCollapse, pYield: ov.pYield, hoop: p * R / t };
    }
    need(!(inp.support === 'void' && !ring), 'Incompatible combination: an unsupported sheath is selected but the construction has no metallic sheath (thickness 0).');
    // armour squeeze on the core under tension (helical-tendon model, ends restrained against rotation)
    const eps = T / cab.EAfixed;
    const qArm = b.armourLayers.reduce((a, L) => a + helixContactPressure({ n: L.nWires, wireForce: Math.max(L.E * wireStrain({ alpha: L.alpha, R: L.R, lay: L.lay, epsAxial: eps, twist: 0 }).strain, 0) * L.wireArea, R: L.R, alpha: L.alpha }).pressure, 0);
    // the tensioner pad load is carried across the cable on a chord: mean bearing stress = q / chord width
    const delta = inp.indent / 100 * OD, chord = 2 * Math.sqrt(delta * (OD - delta));
    const crushEst = inp.sigmaIns * 1e6 * chord;                         // N/m, screening
    const volStrain = p / (inp.bulkK * 1e9);
    // rigid joint housing (steel cylinder)
    const joint = inp.jointD > 0 && inp.nJoints > 0 ? (() => {
      const D = inp.jointD / 1e3, t = inp.jointT / 1e3;
      need(t < D / 4, 'The joint-housing wall is thicker than a quarter of its diameter: check the housing dimensions.');
      const c = run(() => collapsePressure({ D, t, fy: inp.jointFy * 1e6, f0: inp.jointOval / 100 }), 'Joint-housing collapse pressure');
      // wall for which pc/γ = p, by bisection on the monotonic collapse pressure
      let lo = 1e-4, hi = D / 4.01;
      for (let i = 0; i < 60; i++) { const mid = 0.5 * (lo + hi); if (collapsePressure({ D, t: mid, fy: inp.jointFy * 1e6, f0: inp.jointOval / 100 }).pc / inp.gammaP >= p) hi = mid; else lo = mid; }
      const mass = (tt) => Math.PI * (D - tt) * tt * 7850 * inp.jointL * inp.nJoints / 1e3;   // t of steel, shells only
      return { D, t, ...c, tReq: hi, steel: mass(t), steelReq: mass(hi) };
    })() : null;
    // handling: tensioner, chute, stack
    const nTracks = Number(inp.nTracks);
    const hand = cableHandlingLimits({ outerDiameter: OD, tension: T, bendRadius: inp.chuteR, sidewallLimit: inp.swpLimit * 1e3, squeeze: inp.squeeze * 1e3, squeezeLimit: inp.squeezeLimit * 1e3, nTracks, mu: inp.muPad, safetyFactor: inp.gripSf });
    const gripCap = inp.muPad * nTracks * inp.squeeze * 1e3 * inp.padLength / inp.gripSf;   // N the tensioner can hold
    // hexagonal close packing: each cable rests on two below at 30° from the vertical → N = W/(2·cos30°) per contact line
    const layerPitch = OD * Math.sqrt(3) / 2, nAbove = inp.stackH > 0 ? Math.floor(inp.stackH / layerPitch) : 0;
    const stackLine = nAbove * cab.mass * G / (2 * Math.cos(Math.PI / 6));
    const stackMax = inp.crushLimit * 1e3 * 2 * Math.cos(Math.PI / 6) / (cab.mass * G) * layerPitch;   // m of stack at the crush limit
    // rock: static overburden on the projected width and falling-stone impact at terminal velocity
    const rockLine = inp.rockCover > 0 ? (inp.rockBulk - inp.rhoW * (1 - 0.4)) * G * inp.rockCover * OD : 0;   // submerged bulk weight, porosity 0.4 water-filled
    const vol = Math.PI * inp.rockD ** 3 / 6, rockMass = inp.rockRho * vol;
    const imp = inp.rockCover > 0 ? run(() => droppedObjectImpact({ mass: rockMass, volume: vol, CD: inp.rockCd, area: Math.PI * inp.rockD ** 2 / 4, Ca: 0.5, rhoW: inp.rhoW }), 'Falling-rock impact') : { terminalVelocity: 0, impactEnergy: 0 };
    const bandRows = (inp.bands || []).filter((r) => r.maxEnergy > 0).map((r) => ({ label: String(r.label || ''), maxEnergy: Number(r.maxEnergy) })).sort((x, y) => x.maxEnergy - y.maxEnergy);
    need(bandRows.length > 0, 'The impact classification table needs at least one band with a positive energy.');
    const cls = run(() => cableImpactClass({ energy: imp.impactEnergy, bands: bandRows }), 'Impact classification (energies must be distinct)');
    const dRockMax = (() => {            // largest stone whose terminal-velocity impact stays within the rated energy (E ∝ d⁴)
      const e1 = droppedObjectImpact({ mass: inp.rockRho * Math.PI / 6, volume: Math.PI / 6, CD: inp.rockCd, area: Math.PI / 4, Ca: 0.5, rhoW: inp.rhoW }).impactEnergy;
      return (inp.impactRated / e1) ** 0.25;
    })();
    const checks = [
      ...(ring && inp.support === 'void' ? [check('Unsupported metallic sheath under water pressure', p * inp.gammaP / 1e6, ring.pCollapse / 1e6, 'MPa', 'Timoshenko & Gere ch. 7 (oval ring, first yield); CIGRE TB 490', 'Factored hydrostatic pressure against the collapse pressure of the sheath acting as a free ring. Beyond it the sheath buckles inwards into the gap, creases and loses its water-tightness.')] : []),
      check('Radial squeeze of the core by the armour', qArm / 1e6, inp.sigmaIns, 'MPa', 'Helical-tendon contact pressure (Féret & Bournazel 1987); CIGRE TB 623 §5', 'Inward pressure of the tensioned armour wires on the core against the sustained compressive stress the insulation tolerates. Too much squeeze indents the insulation screen under the wires.'),
      ...(joint ? [check('Joint housing collapse under water pressure', p * inp.gammaP / 1e6, joint.pc / 1e6, 'MPa', 'DNV-ST-F101 (2021) Eq. 5.10 (characteristic collapse pressure); CIGRE TB 490', 'Factored hydrostatic pressure against the collapse pressure of the steel joint casing with its ovality. A collapsed housing crushes the joint insulation and floods the joint.')] : []),
      check('Tensioner squeeze per track', inp.squeeze, inp.squeezeLimit, 'kN/m', 'CIGRE TB 623 §5 (squeeze / crush test); manufacturer handling limits', 'Pad force per metre of track against the tested squeeze. Excess squeeze ovalises the cable and flattens the cores against each other.'),
      check('Tensioner grip', T * 1e-3, Math.max(gripCap * 1e-3, 1e-9), 'kN', 'Friction statics T ≤ μ·n·q·L / SF', 'Cable tension against the factored friction capacity of the tensioner. Without enough grip the cable slips through the pads and runs away over the chute.'),
      check('Sidewall pressure on the chute', hand.sidewallPressure / 1e3, inp.swpLimit, 'kN/m', 'CIGRE TB 623 §5.4; SWP = T/R', 'Radial line load where the tensioned cable bears on the chute. Above the tested value the insulation is indented on the inside of the bend.'),
      ...(nAbove > 0 ? [check('Storage stack on the bottom layer', stackLine / 1e3, inp.crushLimit, 'kN/m', 'Hexagonal packing statics; CIGRE TB 623 crush test', 'Contact line load on a cable of the bottom layer from the stack above against the crush resistance. Long storage under excess load leaves permanent ovality.')] : []),
      ...(inp.rockCover > 0 ? [
        check('Rock overburden', rockLine / 1e3, inp.crushLimit, 'kN/m', 'Overburden statics (submerged bulk weight × cover × diameter)', 'Static weight of the rock column on the cable against the crush resistance; normally small, it matters for deep berms and soft-bedded cables.'),
        check('Falling-rock impact', imp.impactEnergy, inp.impactRated, 'J', 'DNV-RP-F107 (2019) §5 (terminal velocity and effective impact energy); CIGRE TB 623 impact test', 'Energy of the largest stone arriving at terminal velocity against the rated impact energy. Larger stones dent the armour and bruise the core: use a finer filter layer first.'),
      ] : []),
    ];
    const gov = worst(checks);
    const depths = linspace(0, Math.max(2 * inp.depth, 300), 30), pD = depths.map((d) => inp.rhoW * G * d / 1e6);
    const pLim = joint ? joint.pc / inp.gammaP / 1e6 : null, pRing = ring ? ring.pCollapse / inp.gammaP / 1e6 : null;
    const Ts = linspace(0, Math.max(2 * inp.tension, 50), 25);
    const dSave = joint ? joint.steel - joint.steelReq : 0;
    return {
      kpis: [
        { key: 'pExt', label: 'Hydrostatic pressure at maximum depth', value: round(p / 1e6, 5), unit: 'MPa', explain: 'External water pressure ρ·g·d acting on the serving, the sheath and every rigid accessory at the deepest point of the route.' },
        { key: 'volStrain', label: 'Volumetric compression of the insulation', value: round(volStrain * 100, 4), unit: '%', explain: 'Bulk strain p/K of the extruded core under the water pressure. A supported solid core only shrinks slightly; it is the voids and interfaces (joints, strand interstices) that pressure tests probe.' },
        { key: 'pSheathEl', label: 'Elastic buckling pressure of the sheath as a free ring', value: round(ring ? ring.pCr / 1e6 : 0, 5), unit: 'MPa', explain: 'E·t³/(4(1 − ν²)R³) for the metallic sheath with no support from the core: the pressure it could carry on its own over a gap.' },
        { key: 'pSheath', label: 'Collapse pressure of the unsupported sheath', value: round(ring ? ring.pCollapse / 1e6 : 0, 5), unit: 'MPa', explain: 'First-yield collapse of the oval free ring. It is a lower bound that governs only where the sheath is not supported by the core.' },
        { key: 'depthSheath', label: 'Water depth the unsupported sheath can take', value: round(ring ? ring.pCollapse / inp.gammaP / (inp.rhoW * G) : 0), unit: 'm', explain: 'Depth at which the factored pressure equals the free-ring collapse pressure. Deeper than this the sheath relies entirely on the core behind it.' },
        { key: 'qArmour', label: 'Armour contact pressure on the core', value: round(qArm / 1e6, 4), unit: 'MPa', limit: inp.sigmaIns, sense: 'max', explain: 'Radial pressure from the helical armour wires at the applied tension, smeared over the circumference.' },
        ...(joint ? [
          { key: 'jointPel', label: 'Joint housing elastic collapse pressure', value: round(joint.pel / 1e6, 5), unit: 'MPa', explain: '2E(t/D)³/(1 − ν²): the buckling pressure of the perfect casing.' },
          { key: 'jointPc', label: 'Joint housing characteristic collapse pressure', value: round(joint.pc / 1e6, 5), unit: 'MPa', explain: 'Collapse pressure of the casing with yielding and ovality combined; divide by the safety factor for the permissible pressure.' },
          { key: 'jointDepth', label: 'Rated water depth of the joint housing', value: round(joint.pc / inp.gammaP / (inp.rhoW * G)), unit: 'm', limit: inp.depth, sense: 'min', explain: 'Depth at which the factored water pressure reaches the collapse pressure of the casing; it must exceed the route depth.' },
          { key: 'jointTreq', label: 'Required housing wall thickness', value: round(joint.tReq * 1e3), unit: 'mm', explain: 'Wall that just satisfies the collapse check at the route depth for this diameter, grade and ovality.' },
        ] : []),
        { key: 'gripLength', label: 'Required tensioner contact length', value: round(hand.gripLength), unit: 'm', limit: inp.padLength, sense: 'max', explain: 'SF·T/(μ·n·q): track length needed to hold the tension at the chosen squeeze. It must fit within the tensioner.' },
        { key: 'swp', label: 'Sidewall pressure on the chute', value: round(hand.sidewallPressure / 1e3), unit: 'kN/m', limit: inp.swpLimit, sense: 'max', explain: 'T/R on the chute or sheave at the applied tension.' },
        { key: 'chuteRmin', label: 'Minimum chute radius for the sidewall limit', value: round(hand.minRadiusForSidewall), unit: 'm', explain: 'Radius at which the sidewall pressure equals its allowable at this tension.' },
        { key: 'stackLoad', label: 'Line load on the bottom layer of the stack', value: round(stackLine / 1e3, 4), unit: 'kN/m', explain: 'Contact force per metre on one cable at the bottom of the storage stack (hexagonal packing).' },
        { key: 'stackMax', label: 'Stack height at the crush limit', value: round(stackMax), unit: 'm', explain: 'Stack height at which the bottom-layer contact load equals the crush resistance; compare with the carousel or tank depth.' },
        { key: 'rockLoad', label: 'Rock overburden line load', value: round(rockLine / 1e3, 4), unit: 'kN/m', explain: 'Submerged weight of the rock column above the cable per metre.' },
        { key: 'vRock', label: 'Terminal velocity of the largest stone', value: round(imp.terminalVelocity), unit: 'm/s', explain: 'Fall velocity in water at which weight, buoyancy and drag balance.' },
        { key: 'impact', label: 'Impact energy of the largest stone', value: round(imp.impactEnergy, 4), unit: 'J', limit: inp.impactRated, sense: 'max', explain: 'Effective energy including added mass at terminal velocity; compare with the impact-test energy of the cable.' },
        { key: 'crushEst', label: 'Indentation estimate of the crush load', value: round(crushEst / 1e3), unit: 'kN/m', explain: 'Allowable bearing stress × chord width at the permissible flattening. A screening cross-check of the test value entered as the crush limit, not a replacement for it.' },
        { key: 'uc', label: 'Governing utilisation', value: round(gov.uc), unit: '-', limit: 1, sense: 'max', explain: 'Largest demand / capacity ratio of the crushing and pressure checks.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Water pressure against depth and pressure limits', x: { label: 'Water depth', unit: 'm' }, y: { label: 'Pressure', unit: 'MPa' },
          series: [{ name: 'Hydrostatic pressure', x: depths, y: pD }, ...(pLim != null ? [{ name: 'Joint housing permissible', x: [depths[0], depths[29]], y: [pLim, pLim], dash: true }] : []), ...(pRing != null ? [{ name: 'Unsupported sheath permissible', x: [depths[0], depths[29]], y: [pRing, pRing], dash: true }] : []), { name: 'Route depth', x: [inp.depth, inp.depth], y: [0, pD[29]], dash: true }],
          explain: 'Where the pressure line crosses a dashed limit is the rated depth of that component. The route-depth marker must lie to the left of the joint-housing crossing.' },
        { type: 'line', title: 'Handling loads against tension', x: { label: 'Tension', unit: 'kN' }, y: { label: 'Utilisation', unit: '-' },
          series: [{ name: 'Sidewall pressure on the chute', x: Ts, y: Ts.map((t) => t / inp.chuteR / inp.swpLimit) }, { name: 'Tensioner grip', x: Ts, y: Ts.map((t) => t * 1e3 / Math.max(gripCap, 1e-9)) }, { name: 'Limit', x: [0, Ts[24]], y: [1, 1], dash: true }, { name: 'Applied tension', x: [inp.tension, inp.tension], y: [0, 1], dash: true }],
          explain: 'Both handling utilisations grow linearly with tension. The first line to reach 1.0 gives the tension at which the chute radius or the tensioner length runs out.' },
        { type: 'bar', title: 'Crushing line loads against the crush limit', categories: ['Tensioner squeeze', 'Chute sidewall', 'Storage stack', 'Rock overburden'], series: [{ name: 'Demand (kN/m)', y: [inp.squeeze, hand.sidewallPressure / 1e3, stackLine / 1e3, rockLine / 1e3].map((v) => round(v, 4)) }, { name: 'Allowable (kN/m)', y: [inp.squeezeLimit, inp.swpLimit, inp.crushLimit, inp.crushLimit] }],
          explain: 'Each radial line load beside its own allowable. The pair with the smallest gap is the handling step that needs the tightest control offshore.' },
      ],
      tables: [
        { title: 'Crushing and pressure summary', explain: 'One row per load: demand, capacity and where the capacity comes from. Rows marked "test" need the manufacturer value before design use.',
          columns: [{ key: 'load', label: 'Load' }, { key: 'd', label: 'Demand' }, { key: 'c', label: 'Capacity' }, { key: 'u', label: 'Unit' }, { key: 'src', label: 'Capacity basis' }],
          rows: [
            { load: 'Water pressure on the cable', d: round(p / 1e6), c: '—', u: 'MPa', src: 'water-pressure withstand test (CIGRE TB 490)' },
            ...(ring ? [{ load: 'Sheath as a free ring', d: round(p * inp.gammaP / 1e6), c: round(ring.pCollapse / 1e6), u: 'MPa', src: inp.support === 'void' ? 'calculation (checked)' : 'calculation (information: core supports the sheath)' }] : []),
            ...(joint ? [{ load: 'Joint housing', d: round(p * inp.gammaP / 1e6), c: round(joint.pc / 1e6), u: 'MPa', src: 'calculation' }] : []),
            { load: 'Armour squeeze on the core', d: round(qArm / 1e6), c: inp.sigmaIns, u: 'MPa', src: 'material data' },
            { load: 'Tensioner squeeze', d: inp.squeeze, c: inp.squeezeLimit, u: 'kN/m', src: 'test' },
            { load: 'Chute sidewall pressure', d: round(hand.sidewallPressure / 1e3), c: inp.swpLimit, u: 'kN/m', src: 'test' },
            { load: 'Storage stack', d: round(stackLine / 1e3), c: inp.crushLimit, u: 'kN/m', src: 'test' },
            { load: 'Rock overburden', d: round(rockLine / 1e3), c: inp.crushLimit, u: 'kN/m', src: 'test' },
            { load: 'Falling rock', d: round(imp.impactEnergy), c: inp.impactRated, u: 'J', src: 'test' },
          ] },
      ],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        `Water pressure ${fmt(p / 1e6)} MPa at ${inp.depth} m; the core compresses by ${fmt(volStrain * 100)} % in volume.${ring ? ` The sheath alone would collapse at ${fmt(ring.pCollapse / 1e6)} MPa (${fmt(ring.pCollapse / (inp.rhoW * G))} m of water before the safety factor): ${inp.support === 'void' ? 'it is treated as unsupported.' : 'it relies on the extruded core, which is the normal construction.'}` : ''}`,
        `Storage: ${nAbove} layer(s) above the bottom one at ${fmt(layerPitch * 1e3)} mm layer pitch.`,
        inp.rockCover > 0 ? `Falling rock of ${inp.rockD} m (${fmt(rockMass)} kg) arrives at ${fmt(imp.terminalVelocity)} m/s with ${fmt(imp.impactEnergy)} J: class '${cls.label}'. The rated energy allows stones up to ${fmt(dRockMax)} m.` : 'No rock placement on this cable.',
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}. Water pressure ${fmt(p / 1e6)} MPa, sidewall pressure ${fmt(hand.sidewallPressure / 1e3)} kN/m, tensioner needs ${fmt(hand.gripLength)} m of ${inp.padLength} m, stack load ${fmt(stackLine / 1e3)} kN/m, rock impact ${fmt(imp.impactEnergy)} J.`,
          gov.uc > 1 ? (gov.name.startsWith('Joint') ? `Increase the housing wall to at least ${fmt(joint.tReq * 1e3)} mm, use a higher grade or pressure-compensate the joint.` : gov.name.startsWith('Tensioner grip') ? `Use ${fmt(hand.gripLength)} m of track, add a track or raise the squeeze towards its limit of ${inp.squeezeLimit} kN/m; do not rely on friction above ${inp.muPad} without a pull test.` : gov.name.startsWith('Sidewall') ? `Increase the chute radius to at least ${fmt(hand.minRadiusForSidewall)} m or lower the lay tension.` : gov.name.startsWith('Storage') ? `Limit the stack to ${fmt(stackMax)} m or interleave load-spreading layers.` : gov.name.startsWith('Falling') ? `Place a filter layer with stones no larger than ${fmt(dRockMax)} m directly on the cable before the armour rock.` : gov.name.startsWith('Unsupported') ? 'Support the sheath (fill the gap, shorter sleeve) or use a thicker or stronger sheath; a free lead sheath carries very little external pressure.' : 'Reduce the governing load or obtain a higher tested capacity for it.')
            : 'Replace the indicative crush, squeeze, sidewall and impact limits with the type-test values of the actual cable and repeat; record the tensioner settings and stack height as installation hold points.',
          `CIGRE TB 623 (2015) §5; CIGRE TB 490 (2012); ${STD_CABLE}; DNV-RP-F107 (2019)`),
        adv('economic', 'info', joint ? `The ${inp.nJoints} joint housing(s) use ${fmt(joint.steel)} t of steel shell at ${inp.jointT} mm wall; ${fmt(joint.tReq * 1e3)} mm is the wall the pressure needs.` : `Handling governs: the tensioner needs ${fmt(hand.gripLength)} m of contact for ${inp.tension} kN at ${inp.squeeze} kN/m squeeze.`,
          joint && dSave > 0 ? `Up to ${fmt(dSave)} t of housing steel could be removed on pressure grounds — small beside the cost of a failed joint, so keep the margin unless weight on the lay spread is critical. A crushed or slipped cable offshore costs a repair joint and of the order of 10–20 vessel days, far more than a longer tensioner or a larger chute.`
            : 'A crushed or slipped cable offshore costs a repair joint and of the order of 10–20 vessel days; a longer tensioner, a larger chute or a finer filter layer are cheap by comparison. Fix these in the installation engineering rather than in the cable design.',
          'AACE RP 18R-97 (quantity-based estimating); risk-based contingency per AACE RP 44R-08',
          `${joint ? `Housing steel ≈ ${fmt(joint.steel * CO2_STEEL)} t CO₂e (screening ${CO2_STEEL} t CO₂e/t steel). ` : ''}${inp.rockCover > 0 ? `A 0.3 m filter layer adds roughly ${fmt(inp.rockBulk * 0.3 * (2 * 3 * 0.3 + OD) / 1e3)} t of rock per metre of berm (≈ ${fmt(inp.rockBulk * 0.3 * (2 * 3 * 0.3 + OD) / 1e3 * CO2_ROCK * 1e3)} kg CO₂e, screening ${CO2_ROCK} t CO₂e/t) but avoids impact damage and an early repair. ` : ''}Avoiding one offshore repair saves the manufacture of a spare length and a repair-vessel campaign.`),
      ],
      publish: { 'cable.crushLimit': inp.crushLimit * 1e3, 'cable.waterPressure': p },
      quantities: joint ? { steel_t: round(joint.steel) } : {},
    };
  },
  verify: [
    // Hydrostatics: p = ρ·g·d = 1025 × 9.80665 × 1000 = 10.052 MPa.
    // Single-core default build-up: conductor radius sqrt(630e-6/(0.9π)) = 14.927 mm, + 1 + 9 + 1 = 25.927 mm under the lead sheath,
    // t = 1.8 mm → R = 26.827 mm. Free ring in plane strain (Levy / Timoshenko): p_cr = 3·[E·t³/(12(1 − ν²))]/R³
    // = 16e9 × (1.8e-3)³ / (4 × (1 − 0.42²) × (26.827e-3)³) = 93.312 / (4 × 0.8236 × 1.93071e-5) = 1.467 MPa.
    // Joint housing D = 400 mm, t = 12 mm: p_el = 2E(t/D)³/(1 − ν²) = 2 × 207e9 × 0.03³/0.91 = 12.284 MPa.
    { name: 'Hydrostatic pressure, free-ring buckling of the sheath and elastic collapse of the joint housing',
      source: 'Hydrostatics; Levy (1884) / Timoshenko & Gere ch. 7 (p_cr = 3EI/R³); DNV-ST-F101 Eq. 5.11 (hand calculation)',
      inputs: { cores: '1', depth: 1000 },
      expect: [{ kpi: 'pExt', value: 10.052, tol: 0.001 }, { kpi: 'pSheathEl', value: 1.467, tol: 0.003 }, { kpi: 'jointPel', value: 12.284, tol: 0.002 }] },
    // Tensioner: L = SF·T/(μ·n·q) = 1.5 × 100e3/(0.3 × 2 × 40e3) = 6.25 m. Chute: SWP = T/R = 100/4 = 25 kN/m.
    // Falling stone d = 0.25 m: V = π·d³/6 = 8.1812e-3 m³, m = 2650·V = 21.680 kg, A = π·d²/4 = 0.049087 m²;
    // v_T = sqrt(2g(m − ρV)/(ρ·Cd·A)) = sqrt(2 × 9.80665 × (21.680 − 8.3857)/(1025 × 0.049087)) = sqrt(5.1824) = 2.2765 m/s;
    // E = ½(m + 0.5·ρ·V)·v² = 0.5 × (21.680 + 4.1929) × 5.1824 = 67.04 J.
    { name: 'Tensioner grip length, sidewall pressure and terminal-velocity impact energy',
      source: 'Friction statics; SWP = T/R (CIGRE TB 623); DNV-RP-F107 §5 (hand calculation)',
      inputs: { tension: 100 },
      expect: [{ kpi: 'gripLength', value: 6.25, tol: 0.001 }, { kpi: 'swp', value: 25, tol: 0.001 }, { kpi: 'vRock', value: 2.2765, tol: 0.002 }, { kpi: 'impact', value: 67.04, tol: 0.003 }] },
  ],
};
