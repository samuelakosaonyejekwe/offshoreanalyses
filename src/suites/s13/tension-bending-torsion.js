// Cross-section mechanics: axial, bending and torsional stiffness of the armoured cable, load sharing under
// combined tension, curvature and torque, armour-wire and conductor stresses, allowable tension, bend radius
// limits and coiling ability.
import { cablePullingTension, cableHandlingLimits, cableCoiling } from '../../engine/xsection/cablemech.js';
import { criticalCurvature, wireBendingStress, helixContactPressure, momentCurvatureLoop, wireStrain } from '../../engine/xsection/helix.js';
import { cableInputs, buildMech, adv, check, worst, sev, fmt, round, run, need, linspace, DEG, G, STD_CABLE, LIMIT_GENERIC, CONFIRM, CO2_STEEL, CO2_COPPER, CO2_ALUMINIUM, CO2_LEAD, CO2_POLYMER } from './_common.js';

export default {
  id: 's13-tension-bending-torsion', suite: 13, title: 'Axial tension, bending and torsion',
  summary: 'Stiffnesses of the cable cross-section and its utilisation under combined tension, curvature and torque, with allowable tension and bend radius.',
  covers: [3], fidelity: 'analytical',
  standards: ['CIGRE TB 623 (2015) Recommendations for mechanical testing of submarine cables', 'CIGRE TB 862 (2022) (dynamic cables)', `${STD_CABLE}`, 'IEC 63026 (2019) (submarine cables 6–60 kV: mechanical tests)'],
  models: ['Axial tension–bending–torsion relationships of helically armoured cables (helical-tendon model)', 'Euler–Bernoulli bending with stick–slip armour (bilinear moment–curvature)', 'Coulomb friction between armour layers', 'Torque balance of contra-wound layers', 'Capstan / sidewall bearing pressure'],
  limits: `${LIMIT_GENERIC} The section is radially rigid (no radial contraction of the core), so the axial stiffness is an upper estimate; polymer moduli are room-temperature values, and bitumen-bound armour, temperature dependence and local 3-D contact (finite-element cross-section models) are not represented. Timoshenko shear deformation is negligible for these slender sections and is not included.`,
  publishes: ['cable.OD', 'cable.massPerM', 'cable.wSub', 'cable.EA', 'cable.EI', 'cable.GJ', 'cable.allowTension', 'cable.mbrDynamic', 'cable.mbrStatic'],
  searchTerms: 'submarine cable armour tension bending torsion stiffness CIGRE TB 623',
  inputs: [
    ...cableInputs({ mech: true }),
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'tension', label: 'Applied tension', unit: 'kN', type: 'number', default: 90, min: 0, max: 5000, group: 'Applied loads', role: 'bc', bind: 'cable.topTension', design: true, uncertain: { dist: 'lognormal', cov: 0.1 }, help: 'Governing tension from the configuration, laying or pull-in analysis.' },
    { key: 'bendRadius', label: 'Bend radius at that tension', unit: 'm', type: 'number', default: 4, min: 0.3, max: 1e4, group: 'Applied loads', role: 'bc', design: true, help: 'Chute, sheave, bend-stiffener or sag-bend radius acting together with the tension.' },
    { key: 'torque', label: 'Applied torque', unit: 'kN·m', type: 'number', default: 0, min: -200, max: 200, group: 'Applied loads', role: 'bc' },
    { key: 'ends', label: 'Torsional end condition', type: 'select', default: 'fixed', group: 'Applied loads', role: 'bc', options: [{ value: 'fixed', label: 'Both ends restrained against rotation (installed cable)' }, { value: 'free', label: 'One end free to rotate (swivel / free end)' }],
      help: 'A free end with an applied torque of zero is the torque-free case; fixed ends with a non-zero applied torque is over-constrained and rejected.' },
    { key: 'length', label: 'Cable length for quantities and twist', unit: 'm', type: 'number', default: 1000, min: 1, max: 5e5, group: 'Applied loads', role: 'geometry', bind: 'cable.dynamicLength' },
    { key: 'mu', label: 'Friction coefficient armour / bedding', unit: '-', type: 'number', default: 0.2, min: 0.02, max: 0.6, group: 'Applied loads', role: 'material', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'No normative value: 0.1–0.3 for steel on polymer yarn. Calibrate against the bending-stiffness test.' },
    { key: 'coilRadius', label: 'Coiling / turntable radius', unit: 'm', type: 'number', default: 4, min: 0.5, max: 30, group: 'Applied loads', role: 'operation' },
    { key: 'etaArmour', label: 'Allowable armour stress / tensile strength', unit: '-', type: 'number', default: 0.4, min: 0.1, max: 1, group: 'Acceptance criteria', role: 'criteria', help: `CIGRE TB 623 leaves the permissible armour stress to the tensile-bending type test. ${CONFIRM}` },
    { key: 'condStress', label: 'Allowable conductor stress', unit: 'MPa', type: 'number', default: 50, min: 5, max: 150, group: 'Acceptance criteria', role: 'criteria', help: `Pulling practice: 50 MPa copper, 30 MPa aluminium. ${CONFIRM}` },
    { key: 'mbrDynFactor', label: 'Minimum bend radius under tension / diameter', unit: '-', type: 'number', default: 20, min: 5, max: 40, group: 'Acceptance criteria', role: 'criteria', help: `Typically 15–20 × OD during installation and for dynamic service. ${CONFIRM}` },
    { key: 'mbrStatFactor', label: 'Minimum bend radius installed / diameter', unit: '-', type: 'number', default: 12, min: 4, max: 30, group: 'Acceptance criteria', role: 'criteria', help: `Typically 10–12 × OD without tension. ${CONFIRM}` },
    { key: 'swpLimit', label: 'Allowable sidewall pressure', unit: 'kN/m', type: 'number', default: 50, min: 2, max: 500, group: 'Acceptance criteria', role: 'criteria', help: `From the crush / sidewall-pressure type test of the cable. ${CONFIRM}` },
    { key: 'sheathStrain', label: 'Allowable bending strain in the metallic sheath', unit: '%', type: 'number', default: 1.5, min: 0.1, max: 5, group: 'Acceptance criteria', role: 'criteria', help: 'Single-event limit for lead alloy; cyclic strain is checked in the fatigue module.' },
    { key: 'twistLimit', label: 'Allowable twist', unit: 'deg/m', type: 'number', default: 2, min: 0.05, max: 30, group: 'Acceptance criteria', role: 'criteria', help: 'Manufacturer limit; unbalanced single-armoured cables tolerate more in the lay-opening direction.' },
  ],
  compute(inp, ctx) {
    const b = buildMech(inp, ctx, { rhoW: inp.rhoW }), cab = b.cab;
    const T = inp.tension * 1e3, kappa = 1 / inp.bendRadius, Mt = inp.torque * 1e3;
    need(!(inp.ends === 'fixed' && Mt !== 0), 'Over-constrained torsion: with both ends restrained against rotation the torque is a reaction, not an input. Set the applied torque to 0 or free one end.');
    const K11 = cab.EAfixed, K12 = cab.coupling, K22 = cab.GJfixed;
    // [T; Mt] = [K11 K12; K12 K22]·[ε; τ]
    let eps, twist, torqueOut;
    if (inp.ends === 'fixed') { eps = T / K11; twist = 0; torqueOut = K12 * eps; }
    else { const det = K11 * K22 - K12 * K12; need(det > 0, 'The section has no torsional stiffness: add armour or a sheath.'); eps = (K22 * T - K12 * Mt) / det; twist = (K11 * Mt - K12 * T) / det; torqueOut = Mt; }
    const arm = b.armourLayers;
    // armour wires: axial strain from tension and twist, local bending and friction stress from curvature
    const wires = arm.map((L, i) => {
      const st = wireStrain({ alpha: L.alpha, R: L.R, lay: L.lay, epsAxial: eps, twist }).strain;
      const sigAx = L.E * st, d = Math.sqrt(4 * L.wireArea / Math.PI);
      const q = helixContactPressure({ n: L.nWires, wireForce: Math.max(sigAx, 0) * L.wireArea, R: L.R, alpha: L.alpha });
      const qOut = arm.slice(i + 1).reduce((a, Lo) => a + helixContactPressure({ n: Lo.nWires, wireForce: Math.max(Lo.E * wireStrain({ alpha: Lo.alpha, R: Lo.R, lay: Lo.lay, epsAxial: eps, twist }).strain, 0) * Lo.wireArea, R: Lo.R, alpha: Lo.alpha }).pressure, 0);
      const cc = criticalCurvature({ E: L.E, A: L.wireArea, R: L.R, alpha: L.alpha, n: L.nWires, mu: inp.mu, qOuter: qOut, qInner: qOut + q.pressure });
      const wb = wireBendingStress({ E: L.E, thickness: d, width: d, alpha: L.alpha, kappa });
      const sigFr = kappa >= cc.kappaFullSlip ? cc.sigmaFrictionMax : cc.sigmaFrictionMax * Math.min(1, kappa / Math.max(cc.kappaFullSlip, 1e-12));
      const total = Math.abs(sigAx) + Math.abs(wb.sigmaNormal) + sigFr;
      return { L, d, strain: st, sigAx, sigBend: Math.abs(wb.sigmaNormal), sigFr, total, cc, q: q.pressure, share: L.nWires * sigAx * L.wireArea * Math.cos(L.alpha) / Math.max(T, 1) };
    });
    const cond = cab.layers.find((L) => L.type === 'conductor');
    const cosC = cond.layAngle ? Math.cos(cond.layAngle) : 1;
    const sigCond = cond.E * eps * cosC * cosC;
    const pull = run(() => cablePullingTension(cab, { conductor: b.m.cond.id === 'al' ? 'aluminium' : 'copper', conductorStress: inp.condStress * 1e6, ...(arm.length ? { armourUtilisation: inp.etaArmour } : {}) }), 'Pulling-tension limits');
    // allowable tension: armour stress limit less the bending and friction stress at the applied curvature; conductor stress limit
    const wOut = wires[wires.length - 1];
    const perT = (w) => (T > 0 ? Math.abs(w.sigAx) / T : 0);
    const TallowArm = wires.length ? Math.min(...wires.map((w) => (perT(w) > 0 ? Math.max(0, inp.etaArmour * w.L.uts - w.sigBend - w.sigFr) / perT(w) : pull.armourLimit))) : Infinity;
    const TallowCond = sigCond > 0 ? inp.condStress * 1e6 / (sigCond / T) : pull.conductorLimit;
    const Tallow = Math.min(TallowArm, TallowCond, 1e12);
    need(Tallow > 0, 'The bend radius alone exhausts the allowable armour stress: increase the bend radius.');
    const hand = cableHandlingLimits({ outerDiameter: cab.outerDiameter, installFactor: inp.mbrDynFactor, staticFactor: inp.mbrStatFactor, tension: T, bendRadius: inp.bendRadius, sidewallLimit: inp.swpLimit * 1e3 });
    const coil = cableCoiling(cab, { coilRadius: inp.coilRadius });
    // bending strain of a core's metallic sheath: each helically laid core bends about its own axis
    const sheath = cab.layers.find((L) => L.type === 'sheath'), rSheath = sheath ? sheath.ro : cab.coreDiameter / 2;
    const epsSheath = kappa * rSheath;
    const Mf = wires.reduce((a, w) => a + w.cc.frictionMoment, 0);
    const EIstick = cab.EImax, EIslip = cab.EImin;
    const loop = run(() => momentCurvatureLoop({ EIstick, EIslip, frictionMoment: Math.max(Mf, 1e-6), kappaAmp: Math.max(kappa, 1e-6), nPoints: 21 }), 'Moment–curvature loop');
    const checks = [
      ...(wires.length ? [check('Armour wire stress (tension + bending + friction)', Math.max(...wires.map((w) => w.total / w.L.uts)) * 1, inp.etaArmour, '× UTS', 'CIGRE TB 623 §5 (tensile-bending test rationale); helical-tendon model', 'Largest wire stress as a fraction of the wire tensile strength against the allowable fraction. Above it wires yield at the bend and the armour loosens or birdcages on unloading.')] : []),
      check('Conductor stress', sigCond / 1e6, inp.condStress, 'MPa', 'CIGRE TB 623 (conductor elongation); cable-pulling practice', 'Axial stress carried by the conductors at the shared strain. Excess stress stretches the strands permanently and disturbs the conductor screen.'),
      check('Bend radius under tension', hand.mbrInstall, inp.bendRadius, 'm', `${STD_CABLE} (MBR × OD); CIGRE TB 623`, 'Demand is the minimum bend radius under tension, capacity the actual radius. A tighter bend overstrains the sheath and screens.'),
      check('Sidewall pressure', hand.sidewallPressure / 1e3, inp.swpLimit, 'kN/m', 'CIGRE TB 623 §5.4 (sidewall pressure / crush test); SWP = T/R', 'Radial line load on the cable where it bears on a chute or sheave under tension. Above the tested value the insulation is indented and the cores are flattened.'),
      check('Metallic-sheath bending strain', epsSheath * 100, inp.sheathStrain, '%', 'CIGRE TB 623; lead-sheath strain practice', 'Outer-fibre strain of the lead sheath of one core bending about its own axis at the applied radius. Excess strain wrinkles or cracks the radial water barrier.'),
      check('Twist', Math.abs(twist) / DEG, inp.twistLimit, 'deg/m', 'Manufacturer twist limit; CIGRE TB 623 (coiling)', 'Rotation per metre caused by tension–torsion coupling or applied torque. Too much twist opens one armour layer and tightens the other, risking birdcaging or core damage.'),
    ];
    const gov = worst(checks);
    // capacity curve: allowable tension against curvature (armour stress limit)
    const kap = linspace(0, Math.min(2 / (inp.mbrStatFactor * cab.outerDiameter), 1.5), 30);
    const capT = kap.map((k) => {
      if (!wires.length) return TallowCond / 1e3;
      const lim = Math.min(...wires.map((w) => {
        const wb = wireBendingStress({ E: w.L.E, thickness: w.d, width: w.d, alpha: w.L.alpha, kappa: k });
        return perT(w) > 0 ? Math.max(0, inp.etaArmour * w.L.uts - Math.abs(wb.sigmaNormal)) / perT(w) : 0;
      }));
      return Math.min(lim, TallowCond) / 1e3;
    });
    const L = inp.length, mk = (pred) => cab.layers.filter(pred).reduce((a, q) => a + q.mass * (q.count ?? 1) / (q.layAngle ? Math.cos(q.layAngle) : 1), 0);
    const mCond = mk((q) => q.type === 'conductor'), mSheath = mk((q) => q.type === 'sheath'), mArm = mk((q) => q.type === 'armour');
    const mPoly = cab.mass - mCond - mSheath - mArm, al = b.m.cond.id === 'al';
    const co2 = (mCond * (al ? CO2_ALUMINIUM : CO2_COPPER) + mSheath * CO2_LEAD + mArm * CO2_STEEL + mPoly * CO2_POLYMER) / 1e3;   // t CO2e per metre × 1e0
    const armEA = arm.reduce((a, q) => a + q.EA, 0);
    return {
      kpis: [
        { key: 'OD', label: 'Overall diameter', value: round(cab.outerDiameter * 1e3), unit: 'mm', explain: 'Diameter over the serving from the layer build-up; it sets buoyancy, drag and the bend-radius limits.' },
        { key: 'mass', label: 'Mass in air', value: round(cab.mass), unit: 'kg/m', explain: 'Sum of all layers including fillers; compare with the data sheet to validate the build-up.' },
        { key: 'wSub', label: 'Submerged weight', value: round(cab.submergedWeight), unit: 'N/m', explain: 'Weight in sea water (flooded interstices not credited): the load that drives catenary tension and on-bottom stability.' },
        { key: 'EA', label: 'Axial stiffness (torque-free)', value: round(cab.EA / 1e6), unit: 'MN', explain: 'Tension per unit strain when the cable is free to rotate; the restrained value is slightly higher. Calculated for a radially rigid core, so an upper estimate.' },
        { key: 'armourEA', label: 'Armour contribution to axial stiffness', value: round(armEA / 1e6, 5), unit: 'MN', explain: 'Σ n·E·A·cos³α of the armour layers: the part of the axial stiffness (and of the tension) carried by the armour.' },
        { key: 'EIslip', label: 'Bending stiffness, full slip', value: round(cab.EImin / 1e3), unit: 'kN·m²', explain: 'Lower bound: every strand and wire bends about its own axis. It governs the global shape at large curvature.' },
        { key: 'EIstick', label: 'Bending stiffness, full stick', value: round(cab.EImax / 1e3), unit: 'kN·m²', explain: 'Upper bound: solid conductors and armour locked by friction. It applies only at very small curvature (and to cold, bitumen-bound cables).' },
        { key: 'GJ', label: 'Torsional stiffness (tension-free)', value: round(cab.GJ / 1e3), unit: 'kN·m²', explain: 'Torque per unit twist with no axial restraint.' },
        { key: 'torquePerTension', label: 'Torque per unit tension (ends fixed)', value: round(cab.torquePerTension * 1e3, 3), unit: 'N·m/kN', explain: 'Torque generated by tension when rotation is prevented. Near zero for contra-wound armour (torque-balanced); a single armour layer unwinds under tension.' },
        { key: 'strain', label: 'Axial strain at the applied tension', value: round(eps * 100, 4), unit: '%', explain: 'Common strain of conductors and armour under the applied tension and end condition.' },
        { key: 'twist', label: 'Twist', value: round(twist / DEG, 4), unit: 'deg/m', limit: inp.twistLimit, sense: 'max', explain: 'Rotation per metre under the applied loads (zero when both ends are fixed). Multiply by the free length for the end rotation.' },
        { key: 'armourStress', label: 'Largest armour wire stress', value: round(wires.length ? Math.max(...wires.map((w) => w.total)) / 1e6 : 0), unit: 'MPa', explain: 'Axial stress from tension plus local wire bending and friction stress at the applied curvature, in the most loaded layer.' },
        { key: 'condStressVal', label: 'Conductor stress', value: round(sigCond / 1e6), unit: 'MPa', limit: inp.condStress, sense: 'max', explain: 'Stress in the conductor metal at the shared axial strain.' },
        { key: 'Tallow', label: 'Allowable tension at this bend radius', value: round(Tallow / 1e3), unit: 'kN', explain: 'Smaller of the armour-stress and conductor-stress limits with the bending and friction stresses at the applied radius deducted. Published as the allowable tension for the global analyses.' },
        { key: 'mbrDyn', label: 'Minimum bend radius under tension', value: round(hand.mbrInstall), unit: 'm', explain: 'Bend-radius multiple × overall diameter for installation and dynamic service.' },
        { key: 'mbrStat', label: 'Minimum bend radius installed', value: round(hand.mbrStatic), unit: 'm', explain: 'Bend-radius multiple × overall diameter for the installed cable without tension.' },
        { key: 'uc', label: 'Governing utilisation', value: round(gov.uc), unit: '-', limit: 1, sense: 'max', explain: 'Largest demand / capacity ratio of the checks below.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Allowable tension against curvature', x: { label: 'Curvature', unit: '1/m' }, y: { label: 'Tension', unit: 'kN' },
          series: [{ name: 'Capacity (armour and conductor limits)', x: kap, y: capT }, { name: 'Applied load', x: [kappa], y: [inp.tension] }, { name: 'MBR under tension', x: [1 / hand.mbrInstall, 1 / hand.mbrInstall], y: [0, Math.max(...capT)], dash: true }],
          explain: 'Interaction diagram of tension and curvature: the operating point must lie below the curve and left of the dashed bend-radius limit. The curve falls as wire bending stress uses up the allowable armour stress.' },
        { type: 'line', title: 'Moment–curvature hysteresis at the applied curvature', x: { label: 'Curvature', unit: '1/m' }, y: { label: 'Bending moment', unit: 'kN·m' }, series: [{ name: 'Stick–slip loop', x: loop.kappa, y: loop.moment.map((v) => v / 1e3) }],
          explain: 'Bilinear bending response: stiff until the armour wires slip, then the slip stiffness. The enclosed area is the energy dissipated per cycle (structural damping); the knee is the slip curvature.' },
        { type: 'bar', title: 'Share of tension by component', categories: ['Conductors', ...wires.map((w) => w.L.name), 'Sheaths and polymers'],
          series: [{ name: 'Share of tension', y: (() => { const c = T > 0 ? cab.cores * cond.EA * cosC ** 3 * eps / T : 0, a = wires.map((w) => w.share); return [c, ...a, Math.max(0, 1 - c - a.reduce((p, q) => p + q, 0))].map((v) => round(v * 100)); })() }],
          explain: 'Percentage of the applied tension carried by each component at the shared strain. The armour should carry most of it; a large conductor share signals an under-armoured design.' },
      ],
      tables: [
        { title: 'Layer build-up and properties', explain: 'One row per layer from the conductor outwards: radii, mass and axial stiffness contribution. Check the build-up against the manufacturer drawing.',
          columns: [{ key: 'name', label: 'Layer' }, { key: 'group', label: 'In' }, { key: 'ri', label: 'Inner radius', unit: 'mm' }, { key: 'ro', label: 'Outer radius', unit: 'mm' }, { key: 'mass', label: 'Mass', unit: 'kg/m' }, { key: 'EA', label: 'EA', unit: 'MN' }],
          rows: cab.layers.map((q) => ({ name: q.name, group: q.group === 'core' ? `core × ${q.count ?? 1}` : 'outer', ri: round(q.ri * 1e3), ro: round(q.ro * 1e3), mass: round(q.mass * (q.count ?? 1)), EA: round(q.EA * (q.count ?? 1) / 1e6) })) },
        ...(wires.length ? [{ title: 'Armour wire stresses', explain: 'One row per armour layer: wire count, axial stress from tension and twist, local bending stress, friction stress, slip curvature and contact pressure on the layer below.',
          columns: [{ key: 'name', label: 'Layer' }, { key: 'n', label: 'Wires' }, { key: 'ax', label: 'Axial', unit: 'MPa' }, { key: 'bend', label: 'Bending', unit: 'MPa' }, { key: 'fr', label: 'Friction', unit: 'MPa' }, { key: 'tot', label: 'Total', unit: 'MPa' }, { key: 'kcr', label: 'Slip curvature', unit: '1/m' }, { key: 'q', label: 'Contact pressure', unit: 'MPa' }],
          rows: wires.map((w) => ({ name: w.L.name, n: w.L.nWires, ax: round(w.sigAx / 1e6), bend: round(w.sigBend / 1e6), fr: round(w.sigFr / 1e6), tot: round(w.total / 1e6), kcr: round(w.cc.kappaCr), q: round(w.q / 1e6) })) }] : []),
      ],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        cab.torqueBalanced ? 'The armour is torque-balanced: tension produces negligible torque or twist.' : `The section is not torque-balanced (${fmt(cab.torquePerTension * 1e3)} N·m per kN): ${inp.ends === 'fixed' ? `the fixed ends react ${fmt(torqueOut / 1e3)} kN·m` : `a free end rotates ${fmt(Math.abs(twist) * L / (2 * Math.PI))} turns over ${fmt(L)} m`}.`,
        coil.note,
        `Specific gravity ${fmt(cab.specificGravity)}; submerged weight ${fmt(cab.submergedWeight / G)} kg/m.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}. Allowable tension at ${inp.bendRadius} m radius is ${fmt(Tallow / 1e3)} kN against ${inp.tension} kN applied; minimum bend radius under tension ${fmt(hand.mbrInstall)} m.`,
          gov.uc > 1 ? (gov.name.startsWith('Armour') ? `Increase the bend radius, raise the armour wire diameter or strength, or add a second armour layer; the tension must drop below ${fmt(Tallow / 1e3)} kN at this radius.` : gov.name.startsWith('Sidewall') ? `Increase the chute or sheave radius to at least ${fmt(hand.minRadiusForSidewall)} m or reduce the tension.` : 'Increase the bend radius or relieve the governing load, then re-check.')
            : 'Confirm EA, EI and the allowable tension against the manufacturer tensile-bending and bending-stiffness tests (CIGRE TB 623 / TB 862) before using the published values in the global analyses.',
          `CIGRE TB 623 (2015); CIGRE TB 862 (2022); ${STD_CABLE}`),
        adv('economic', 'info', `The section weighs ${fmt(cab.mass)} kg/m: ${fmt(mCond)} kg/m conductor, ${fmt(mArm)} kg/m armour steel, ${fmt(mSheath)} kg/m sheath metal and ${fmt(mPoly)} kg/m polymers — ${fmt(cab.mass * L / 1e3)} t over ${fmt(L)} m.`,
          wires.length && gov.uc < 0.7 ? `The armour is utilised to ${fmt(Math.max(...wires.map((w) => w.total / w.L.uts)) / inp.etaArmour * 100)} % of its allowable: a smaller wire or a single layer (where torque balance is not needed) would save up to ${fmt(mArm / Math.max(b.nArm, 1))} kg/m of steel and reduce turntable loading.` : 'Armour and conductor metal dominate cable cost and carousel capacity; optimise wire size against the governing installation tension rather than a generic rule.',
          'AACE RP 18R-97 (quantity-based estimating); whole-life costing ISO 15686-5',
          `Embodied carbon ≈ ${fmt(co2 * L)} t CO₂e over ${fmt(L)} m (screening factors: copper ${CO2_COPPER}, aluminium ${CO2_ALUMINIUM}, lead ${CO2_LEAD}, steel ${CO2_STEEL}, polymer ${CO2_POLYMER} t CO₂e/t). Lead sheaths are under regulatory pressure (REACH): a lead-free barrier removes ${fmt(mSheath * L / 1e3)} t of lead.`),
      ],
      publish: { 'cable.OD': cab.outerDiameter, 'cable.massPerM': cab.mass, 'cable.wSub': cab.submergedWeight, 'cable.EA': cab.EA, 'cable.EI': cab.EImin, 'cable.GJ': cab.GJ, 'cable.allowTension': Tallow, 'cable.mbrDynamic': hand.mbrInstall, 'cable.mbrStatic': hand.mbrStatic },
      quantities: { cable_km: L / 1e3, copper_t: round(al ? 0 : mCond * L / 1e3), aluminium_t: round(al ? mCond * L / 1e3 : 0), lead_t: round(b.m.sh.id === 'lead-sheath' ? mSheath * L / 1e3 : 0), steel_t: round(mArm * L / 1e3), polymer_t: round(mPoly * L / 1e3) },
    };
  },
  verify: [
    // Single-core cable, one armour layer of 5 mm wires at 15°: conductor radius sqrt(630e-6/(0.9π)) = 14.927 mm;
    // + 1 + 9 + 1 + 1.8 + 2 + 2 (bedding) = 31.727 mm under the armour → pitch radius R = 34.227 mm.
    // Wires: floor(0.95 × 2π × 34.227 × cos15° / 5) = floor(39.47) = 39. Axial stiffness of the layer
    // n·E·A·cos³α = 39 × 200e9 × 19.635e-6 × 0.9012 = 138.0 MN.
    { name: 'Helical armour layer: axial stiffness n·E·A·cos³α',
      source: 'Witz & Tan (1992) / Féret & Bournazel (1987) helical-tendon stiffness; hand calculation',
      inputs: { cores: '1', nArmour: '1', armAngle: 15, tension: 50 },
      expect: [{ kpi: 'armourEA', value: 138.0, tol: 0.005 }] },
    // Torque-free and torsionally fixed limits bracket the strain: fixed ends ε = T/K11. With cores = 1 and no armour
    // (dWire = 0) the section is concentric: EA = Σ E·A. Conductor 110e9 × 630e-6 = 69.3 MN dominates; XLPE
    // (0.2 GPa) rings 14.927→25.927 mm: π(25.927² − 14.927²)·1e-6 × 0.2e9 = 0.2824 MN; lead 16 GPa × π(27.727² − 25.927²)e-6
    // = 4.855 MN; jacket 0.8e9 × π(29.727² − 27.727²)e-6 = 0.2888 MN; bedding 0.05e9 × π(31.727² − 29.727²)e-6 = 0.0193 MN;
    // serving 0.8e9 × π(35.727² − 31.727²)e-6 = 0.678 MN → EA = 75.42 MN; strain at 50 kN = 0.0663 %.
    { name: 'Concentric unarmoured core: EA = Σ E·A and ε = T/EA',
      source: 'Composite bar in tension (strength of materials); hand calculation',
      inputs: { cores: '1', dWire: 0, tension: 50, condStress: 150, bendRadius: 50 },
      expect: [{ kpi: 'EA', value: 75.42, tol: 0.005 }, { kpi: 'strain', value: 0.0663, tol: 0.005 }] },
  ],
};
