// Cable crossing and interference: the cable laid over an existing pipeline or cable on a separation support —
// profile, lifted length, curvature and contact loads from a beam on unilateral supports — with the separation
// clearance, the protection quantities, the thermal interference between the two assets and the spacing to
// parallel circuits.
import { beamOnSupports, crossingSeparation, crossingQuantities } from '../../engine/pipeline/crossing.js';
import { rockBerm, rockSizeShields, mattressCount } from '../../engine/xsection/burial.js';
import { mutualHeatingRise } from '../../engine/electrical/ampacity.js';
import { lineInputs, lineFrom, adv, check, worst, sev, fmt, round, run, need, linspace, thin, G, STD_CABLE, CONFIRM, CO2_ROCK, vesselCo2Text } from './_common.js';

const CO2_CONCRETE = 0.3;     // t CO2e per m³ concrete (screening)

export default {
  id: 's13-crossing-interference', suite: 13, title: 'Cable crossing and interference',
  summary: 'Shape, bending, contact loads and separation of a cable crossing over an existing pipeline or cable, the rock and mattress quantities, the mutual heating at the crossing and the spacing to neighbouring circuits.',
  covers: [13], fidelity: 'numerical',
  standards: [`${STD_CABLE} (crossings and proximity)`, 'ICPC Recommendation No. 2 (cable routing and crossing criteria) and No. 3 (crossings of telecommunication cables and pipelines / power cables)', 'DNV-ST-F101 (2021) §5 (crossing separation, used as the default clearance)', 'IEC 60287-2-1 §4.2.3 and IEC 60287-3-3 (cables crossing external heat sources)', 'CIRIA C683 The Rock Manual (2007) (berm geometry and rock stability)', 'The Crown Estate / ESCA Guideline No. 6 (proximity of offshore renewable installations and submarine cables)'],
  models: ['Euler–Bernoulli beam on unilateral supports and a Winkler seabed (lift-off by active set)', 'Seabed contact mechanics with touchdown reactions', 'Fourier heat conduction: Kennelly image method for mutual heating', 'Rating reduction from an external temperature rise (IEC 60287-3-3 form)', 'Trapezoidal rock-berm geometry', 'Shields criterion for rock stability', 'Repair-bight spacing rule for parallel circuits'],
  limits: 'The crossing beam has no axial force: residual lay tension stiffens the real cable and shortens the lifted length, so curvature here is on the high side for a tensioned lay and on the low side where slack is laid in. The mutual heating assumes parallel line sources — an upper bound for a crossing, whose heated length is finite (IEC 60287-3-3 gives the exact crossing solution). Electromagnetic interference on the crossed asset is treated in the electromagnetic module; trawl-gear interaction with the berm, scour around it and 3-D berm stability are outside the module. Clearances and spacing factors are overridable defaults to confirm against the crossing agreement.',
  publishes: ['cable.crossingSeparation', 'cable.crossingDerating', 'cable.crossingRock'],
  searchTerms: 'submarine cable pipeline crossing design mattress rock berm separation ICPC',
  inputs: [
    ...lineInputs({ hydro: false }),
    { key: 'asset', label: 'Crossed asset', type: 'select', default: 'pipeline', group: 'Crossed asset', role: 'geometry', options: [{ value: 'pipeline', label: 'Pipeline' }, { value: 'cable', label: 'Power or telecom cable' }] },
    { key: 'crossedD', label: 'Overall diameter of the crossed asset', unit: 'm', type: 'number', default: 0.5, min: 0.02, max: 2, group: 'Crossed asset', role: 'geometry' },
    { key: 'topExisting', label: 'Top of the crossed asset relative to the seabed', unit: 'm', type: 'number', default: 0.5, min: -5, max: 2.5, group: 'Crossed asset', role: 'ic', help: 'As-found position from the pre-lay survey: positive = proud of the seabed, negative = buried.' },
    { key: 'angle', label: 'Crossing angle', unit: 'deg', type: 'number', default: 80, min: 10, max: 90, group: 'Crossed asset', role: 'geometry', design: true, help: 'Angle between the two routes; as near to 90° as the route allows.' },
    { key: 'nCross', label: 'Number of similar crossings', type: 'int', default: 1, min: 1, max: 100, group: 'Crossed asset', role: 'geometry', help: 'Multiplies the quantities only.' },
    { key: 'sepT', label: 'Separation layer thickness above the crossed asset', unit: 'm', type: 'number', default: 0.3, min: 0, max: 2, group: 'Separation support', role: 'bc', design: true, help: 'Concrete mattress or rock pre-lay berm between the two assets: it is the support the cable rests on.' },
    { key: 'supportW', label: 'Width of the support along the cable', unit: 'm', type: 'number', default: 3, min: 0, max: 30, group: 'Separation support', role: 'bc', help: 'Crest length on which the cable bears (mattress width, or pre-lay berm crest). 0 = knife-edge support.' },
    { key: 'settlement', label: 'Expected settlement of the support', unit: 'm', type: 'number', default: 0.05, min: 0, max: 1, group: 'Separation support', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.5 } },
    { key: 'kn', label: 'Seabed vertical stiffness', unit: 'kPa', type: 'number', default: 500, min: 5, max: 1e6, group: 'Separation support', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.5 }, help: 'Force per metre of cable per metre of penetration. Soft clay 20–200, sand 1000 and above.' },
    { key: 'modelL', label: 'Modelled cable length', unit: 'm', type: 'number', default: 80, min: 10, max: 600, group: 'Separation support', role: 'geometry', help: 'Centred on the crossing; it must reach beyond both touchdown points (free ends on the seabed).' },
    { key: 'matL', label: 'Mattress length', unit: 'm', type: 'number', default: 6, min: 1, max: 12, group: 'Protection', role: 'geometry' },
    { key: 'matW', label: 'Mattress width', unit: 'm', type: 'number', default: 3, min: 1, max: 6, group: 'Protection', role: 'geometry' },
    { key: 'matT', label: 'Mattress thickness', unit: 'm', type: 'number', default: 0.3, min: 0.1, max: 0.6, group: 'Protection', role: 'geometry' },
    { key: 'nMatUnder', label: 'Separation mattresses under the cable', type: 'int', default: 1, min: 0, max: 20, group: 'Protection', role: 'geometry', help: '0 when the separation is a rock pre-lay berm.' },
    { key: 'protect', label: 'Post-lay protection', type: 'select', default: 'rock', group: 'Protection', role: 'bc', options: [{ value: 'rock', label: 'Rock berm over the crossing' }, { value: 'mattress', label: 'Concrete mattresses over the cable' }, { value: 'none', label: 'None (cable left exposed)' }] },
    { key: 'bermL', label: 'Protected length along the cable', unit: 'm', type: 'number', default: 40, min: 2, max: 500, group: 'Protection', role: 'geometry', design: true },
    { key: 'cover', label: 'Rock cover above the cable', unit: 'm', type: 'number', default: 0.6, min: 0.2, max: 3, group: 'Protection', role: 'geometry', design: true, showIf: { key: 'protect', in: ['rock'] } },
    { key: 'crestW', label: 'Berm crest width', unit: 'm', type: 'number', default: 1.0, min: 0.3, max: 10, group: 'Protection', role: 'geometry', showIf: { key: 'protect', in: ['rock'] } },
    { key: 'slope', label: 'Berm side slope (horizontal per vertical)', unit: '-', type: 'number', default: 3, min: 1.5, max: 6, group: 'Protection', role: 'geometry', showIf: { key: 'protect', in: ['rock'] } },
    { key: 'rockBulk', label: 'Rock bulk density', unit: 'kg/m³', type: 'number', default: 1600, min: 1200, max: 2200, group: 'Protection', role: 'material', showIf: { key: 'protect', in: ['rock'] } },
    { key: 'd50', label: 'Rock median size D50', unit: 'm', type: 'number', default: 0.08, min: 0.01, max: 0.6, group: 'Protection', role: 'material', design: true, showIf: { key: 'protect', in: ['rock'] } },
    { key: 'Ubed', label: 'Near-bed design velocity (current + wave)', unit: 'm/s', type: 'number', default: 1.2, min: 0, max: 5, group: 'Protection', role: 'environment', bind: 'site.UcBed100', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'fw', label: 'Bed friction factor on the berm', unit: '-', type: 'number', default: 0.03, min: 0.003, max: 0.3, group: 'Protection', role: 'environment', help: 'τ = ½·ρ·f·U². Rough rock in combined waves and current: 0.02–0.05 (Soulsby 1997); raise it on the crest where flow accelerates.' },
    { key: 'shields', label: 'Critical Shields parameter', unit: '-', type: 'number', default: 0.035, min: 0.02, max: 0.06, group: 'Protection', role: 'criteria', help: `0.03–0.035 for no movement, 0.05–0.055 for limited movement (CIRIA C683). ${CONFIRM}` },
    { key: 'daysPerCrossing', label: 'Vessel time per crossing', unit: 'days', type: 'number', default: 1.5, min: 0.1, max: 20, group: 'Protection', role: 'operation', help: 'Survey, mattress placement and rock placement, excluding transit and weather.' },
    { key: 'Wext', label: 'Heat released by the crossed asset', unit: 'W/m', type: 'number', default: 30, min: 0, max: 600, group: 'Thermal interference', role: 'bc', help: 'Losses of a crossed power cable (20–80 W/m) or heat loss of a warm pipeline; 0 for a telecom cable or a cold line.' },
    { key: 'Wown', label: 'Losses of this cable', unit: 'W/m', type: 'number', default: 75, min: 0, max: 600, group: 'Thermal interference', role: 'bc', bind: 'cable.lossPerM' },
    { key: 'rhoCover', label: 'Thermal resistivity of the cover and soil at the crossing', unit: 'K·m/W', type: 'number', default: 1.0, min: 0.2, max: 4, group: 'Thermal interference', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.25 }, help: 'Water-filled rock berm 0.6–1.2; convection in coarse rock lowers the effective value.' },
    { key: 'ambient', label: 'Ambient temperature', unit: '°C', type: 'number', default: 12, min: -2, max: 35, group: 'Thermal interference', role: 'bc', bind: 'site.Tbed' },
    { key: 'thetaMax', label: 'Maximum conductor temperature', unit: '°C', type: 'number', default: 90, min: 40, max: 110, group: 'Thermal interference', role: 'criteria' },
    { key: 'ampacity', label: 'Continuous rating away from the crossing', unit: 'A', type: 'number', default: 630, min: 10, max: 5000, group: 'Thermal interference', role: 'criteria', bind: 'cable.ampacity' },
    { key: 'I', label: 'Load current', unit: 'A', type: 'number', default: 580, min: 0, max: 5000, group: 'Thermal interference', role: 'bc', bind: 'cable.loadCurrent' },
    { key: 'spacing', label: 'Spacing to the nearest parallel cable', unit: 'm', type: 'number', default: 100, min: 1, max: 5000, group: 'Parallel circuits', role: 'geometry', design: true },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 40, min: 2, max: 3000, group: 'Parallel circuits', role: 'environment', bind: 'site.depth' },
    { key: 'spacingFactor', label: 'Required spacing / water depth', unit: '-', type: 'number', default: 2, min: 0.5, max: 6, group: 'Acceptance criteria', role: 'criteria', help: `Room to lay a repair bight beside the cable without crossing its neighbour: 2–3 × depth is common practice. ${CONFIRM}` },
    { key: 'spacingMin', label: 'Minimum spacing regardless of depth', unit: 'm', type: 'number', default: 50, min: 1, max: 1000, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'sepMin', label: 'Minimum vertical separation', unit: 'm', type: 'number', default: 0.3, min: 0.05, max: 2, group: 'Acceptance criteria', role: 'criteria', help: `0.3 m is the usual figure in crossing agreements and DNV-ST-F101. ${CONFIRM}` },
    { key: 'angleMin', label: 'Minimum crossing angle', unit: 'deg', type: 'number', default: 60, min: 10, max: 90, group: 'Acceptance criteria', role: 'criteria', help: `ICPC recommends crossings as near to a right angle as practicable. ${CONFIRM}` },
    { key: 'mbrAllow', label: 'Allowable bend radius (installed)', unit: 'm', type: 'number', default: 1.9, min: 0.2, max: 60, group: 'Acceptance criteria', role: 'criteria', bind: 'cable.mbrStatic' },
    { key: 'spanAllow', label: 'Allowable free-span length', unit: 'm', type: 'number', default: 12, min: 0.5, max: 200, group: 'Acceptance criteria', role: 'criteria', bind: 'cable.spanAllowable', help: 'From the free-span and vibration module; applies while the crossing is not yet covered.' },
    { key: 'loadAllow', label: 'Allowable load on the crossed asset', unit: 'kN', type: 'number', default: 30, min: 0.5, max: 2000, group: 'Acceptance criteria', role: 'criteria', help: `From the owner of the crossed asset (crossing agreement). ${CONFIRM}` },
    { key: 'dThetaAllow', label: 'Allowable temperature rise on the crossed asset', unit: 'K', type: 'number', default: 5, min: 0.5, max: 40, group: 'Acceptance criteria', role: 'criteria', help: `Set in the crossing agreement where the crossed asset is a rated power cable. ${CONFIRM}` },
    { key: 'nElem', label: 'Beam elements', type: 'int', default: 320, min: 40, max: 4000, group: 'Solver', role: 'control', help: 'Supports act at the nearest node: choose a count for which the support edges fall on nodes (element length dividing half the support width).' },
  ],
  mesh: { key: 'nElem', label: 'Beam elements', levels: [160, 320, 640, 1280], qoi: 'reaction' },
  compute(inp) {
    const line = lineFrom(inp), D = line.D, L = inp.modelL;
    need(line.w > 0, `The cable is buoyant in sea water (apparent weight ${fmt(line.w)} N/m): it cannot be laid over a crossing without ballast.`);
    const hSup = inp.topExisting + inp.sepT;            // support crest above the seabed
    // a buried asset whose separation layer stays below the seabed needs no support: the cable lies flat and the soil cover is the separation
    const flat = hSup <= 0.01;
    need(inp.supportW < L / 3, 'The support is wider than a third of the modelled length: increase the modelled cable length.');
    const xs = inp.supportW > 0 ? [L / 2 - inp.supportW / 2, L / 2, L / 2 + inp.supportW / 2] : [L / 2];
    const beam = flat ? (() => { const nn = inp.nElem + 1; return { x: Float64Array.from({ length: nn }, (_, i) => i * L / inp.nElem), deflection: new Float64Array(nn), moment: new Float64Array(nn), seabedReaction: new Float64Array(nn).fill(line.w), reactions: [], maxMoment: 0 }; })()
      : run(() => beamOnSupports({ L, EI: line.EI, w: line.w, supports: xs.map((x) => ({ x, height: hSup })), seabed: { stiffness: inp.kn * 1e3 }, nElem: inp.nElem }), 'The crossing profile did not solve');
    const n = beam.x.length - 1;
    // lifted length = stretch around the crossing where the cable is clear of the seabed
    let i0 = Math.round(n / 2), i1 = i0;
    while (i0 > 0 && beam.seabedReaction[i0] <= 0) i0--;
    while (i1 < n && beam.seabedReaction[i1] <= 0) i1++;
    need(i0 > 1 && i1 < n - 1, `The cable is lifted over the whole modelled length of ${L} m: increase the modelled cable length so that both touchdown points fall inside it.`);
    // touchdown points interpolated between the last node in contact and the first node clear of the seabed
    const cross = (a, c) => { const ua = beam.deflection[a], uc = beam.deflection[c]; return ua < 0 && uc > 0 ? beam.x[a] + (beam.x[c] - beam.x[a]) * (-ua) / (uc - ua) : beam.x[a]; };
    const xA = cross(i0, i0 + 1), xB = cross(i1, i1 - 1);
    const lifted = xB - xA, spanSide = Math.max(0, L / 2 - inp.supportW / 2 - xA, xB - L / 2 - inp.supportW / 2);
    const kMax = beam.maxMoment / line.EI, R = beam.reactions.filter((r) => r.active).reduce((a, r) => a + r.R, 0);
    const tdR = Math.max(...Array.from(beam.seabedReaction));                   // N/m peak seabed reaction at touchdown
    const iMid = Math.round(n / 2), axis = beam.deflection[iMid] + D / 2;
    const sep = crossingSeparation({ axisElevation: axis, Dcross: D, topOfExisting: inp.topExisting, minSeparation: inp.sepMin, settlement: inp.settlement });
    // protection quantities
    const rock = inp.protect === 'rock' ? run(() => rockBerm({ productDiameter: D, cover: inp.cover + Math.max(hSup, 0), crestWidth: inp.crestW, slope: inp.slope, length: inp.bermL, bulkDensity: inp.rockBulk }), 'Rock berm') : null;
    const matOver = inp.protect === 'mattress' ? mattressCount({ length: inp.bermL, mattressLength: inp.matL, coverWidth: inp.matW, mattressWidth: inp.matW }) : null;
    const nMat = inp.nMatUnder + (matOver ? matOver.count : 0);
    const q = crossingQuantities({ mattresses: nMat > 0 ? [{ count: nMat, length: inp.matL, width: inp.matW, thickness: inp.matT, solidFraction: 0.6 }] : [], rhoWater: line.rhoW });
    const tau = 0.5 * line.rhoW * inp.fw * inp.Ubed ** 2;
    const stab = rock ? (tau > 0 ? rockSizeShields({ bedShear: tau, shieldsCritical: inp.shields, rockDensity: 2650, rhoW: line.rhoW, d50: inp.d50 }) : { requiredD50: 0, utilisation: 0 }) : null;
    // load on the crossed asset: cable reaction through the support + the rock column over its width along the crossing length
    const rockOnAsset = rock ? (inp.rockBulk - line.rhoW * 0.6) * G * (inp.cover + D) * inp.crossedD * Math.max(inp.supportW, D) : 0;
    const matOnAsset = inp.nMatUnder > 0 ? q.mattressSubmergedWeight / nMat * Math.min(1, inp.crossedD / inp.matL * 2) : 0;   // share of one mattress bearing on the asset (it bridges to the seabed either side)
    const loadAsset = inp.topExisting > 0 ? R + rockOnAsset + matOnAsset : (R + rockOnAsset) * Math.min(1, inp.crossedD / (inp.crossedD + 2 * Math.abs(inp.topExisting)));   // buried asset: 1:1 load spread through the soil cover
    // thermal interference (parallel line sources, conservative): depths below the top of the cover
    const hasCover = inp.protect !== 'none';
    const topCover = hasCover ? (rock ? inp.cover : inp.matT) : 0.05;            // exposed cable: thin boundary layer, heat goes to the water
    const zOwn = topCover + D / 2, zExt = topCover + D + sep.separation + inp.settlement + inp.crossedD / 2;
    const dThOwn = inp.Wext > 0 ? mutualHeatingRise({ rhoT: inp.rhoCover, target: { x: 0, depth: zOwn }, sources: [{ x: 0, depth: zExt, W: inp.Wext }] }) : 0;
    const dThExt = inp.Wown > 0 ? mutualHeatingRise({ rhoT: inp.rhoCover, target: { x: 0, depth: zExt }, sources: [{ x: 0, depth: zOwn, W: inp.Wown }] }) : 0;
    const dThMax = inp.thetaMax - inp.ambient;
    need(dThMax > 5, 'The maximum conductor temperature must exceed the ambient temperature by at least 5 K.');
    // IEC 60287-3-3 form with dielectric loss neglected: I/I0 = sqrt(1 − Δθ_int/Δθ_max)
    const DF = dThOwn < dThMax ? Math.sqrt(1 - dThOwn / dThMax) : 0;
    const spacingReq = Math.max(inp.spacingMin, inp.spacingFactor * inp.depth);
    const checks = [
      check('Vertical separation after settlement', inp.sepMin, Math.max(sep.separation, 1e-9), 'm', 'DNV-ST-F101 (2021) §5 (crossings); ICPC Recommendation No. 3; crossing agreement', 'Demand is the required clearance, capacity the gap between the underside of the cable and the top of the crossed asset after the support has settled. Contact abrades both assets and couples them electrically and thermally.'),
      check('Bending at the crossing', kMax, 1 / inp.mbrAllow, '1/m', `${STD_CABLE} (installed MBR); CIGRE TB 623`, 'Largest curvature of the cable over the support and at the touchdown points against the allowable installed curvature. A support that is too high or too narrow kinks the cable over its edge.'),
      check('Free span each side of the support (before cover)', spanSide, inp.spanAllow, 'm', 'DNV-RP-F105 principles applied to cables; free-span module of this suite', 'Unsupported length between the support and the touchdown point against the allowable span. Until the rock is placed the span can vibrate in current and is exposed to fishing gear.'),
      check('Load on the crossed asset', loadAsset / 1e3, inp.loadAllow, 'kN', 'Crossing agreement (owner of the crossed asset); statics', 'Cable reaction, mattress share and rock column bearing on the crossed asset against the load its owner accepts. Excess load dents a pipeline coating or crushes a crossed cable.'),
      check('Crossing angle', inp.angleMin, inp.angle, 'deg', 'ICPC Recommendation No. 2 / No. 3 (as near to 90° as practicable)', 'Demand is the minimum acceptable angle, capacity the actual one. Shallow crossings lengthen the zone where neither asset can be repaired or trenched and increase induction and heating.'),
      ...(rock ? [check('Rock stability on the berm', stab.requiredD50, inp.d50, 'm', 'Shields (1936); CIRIA C683 §5.2', 'Rock size needed for the bed shear stress against the size provided. Undersized rock is winnowed away in storms and the cable is re-exposed.'),
        check('Lifted length covered by the berm', lifted, inp.bermL, 'm', 'Project criterion: no residual free span at a crossing', 'Demand is the length over which the cable is clear of the seabed, capacity the protected length. Any lifted length left outside the berm remains a free span.')] : []),
      check('Load current against the rating at the crossing', inp.I, Math.max(inp.ampacity * DF, 1e-9), 'A', 'IEC 60287-3-3 (crossing heat sources); IEC 60287-2-1 §4.2.3.1 (mutual heating)', 'Load current against the continuous rating reduced by the temperature rise the crossed asset imposes on this cable. A hot crossing is a thermal hotspot that limits the whole circuit.'),
      ...(inp.asset === 'cable' && inp.Wown > 0 ? [check('Temperature rise imposed on the crossed cable', dThExt, inp.dThetaAllow, 'K', 'IEC 60287-2-1 §4.2.3.1; crossing agreement', 'Warming of the crossed cable by the losses of this one. It takes rating away from the third-party cable, which is why crossing agreements limit it.')] : []),
      check('Spacing to the parallel cable', spacingReq, inp.spacing, 'm', 'ESCA / The Crown Estate proximity guidance; ICPC Recommendation No. 2', 'Demand is the spacing needed to lay a repair bight without crossing the neighbour, capacity the actual spacing. Too little spacing means a fault on one cable cannot be repaired without disturbing the other.'),
    ];
    const gov = worst(checks);
    const [xp, yp, mp] = thin([Array.from(beam.x, (x) => x - L / 2), Array.from(beam.deflection), Array.from(beam.moment, (m) => m / line.EI)], 300);
    const seps = linspace(0.1, 2, 20);
    const dThS = seps.map((s) => (inp.Wext > 0 ? mutualHeatingRise({ rhoT: inp.rhoCover, target: { x: 0, depth: zOwn }, sources: [{ x: 0, depth: topCover + D + s + inp.crossedD / 2, W: inp.Wext }] }) : 0));
    const rockT = (rock ? rock.tonnes : 0) * inp.nCross, concrete = q.mattressVolume * inp.nCross, days = inp.daysPerCrossing * inp.nCross;
    const co2 = rockT * CO2_ROCK + concrete * CO2_CONCRETE;
    return {
      kpis: [
        { key: 'crestHeight', label: 'Cable elevation over the crossing', value: round(beam.deflection[iMid]), unit: 'm', explain: 'Height of the underside of the cable above the seabed at the crossing point: the support crest less any sag between support edges.' },
        { key: 'separation', label: 'Vertical separation after settlement', value: round(sep.separation), unit: 'm', limit: inp.sepMin, sense: 'min', explain: 'Gap between the underside of the cable and the top of the crossed asset once the support has settled.' },
        { key: 'liftedLength', label: 'Lifted length of cable', value: round(lifted, 5), unit: 'm', explain: 'Length between the two touchdown points over which the cable is clear of the seabed. It is the minimum length the protection must cover.' },
        { key: 'spanSide', label: 'Free span each side of the support', value: round(spanSide), unit: 'm', limit: inp.spanAllow, sense: 'max', explain: 'Unsupported length from the edge of the support to the touchdown point, before any rock is placed.' },
        { key: 'curvMax', label: 'Largest curvature', value: round(kMax, 5), unit: '1/m', limit: round(1 / inp.mbrAllow, 5), sense: 'max', explain: 'Peak curvature M/EI over the support (hogging) or at touchdown (sagging); its reciprocal is the tightest bend radius at the crossing.' },
        { key: 'mbr', label: 'Tightest bend radius', value: round(Math.min(1 / Math.max(kMax, 1e-9), 1e4)), unit: 'm', limit: inp.mbrAllow, sense: 'min', explain: 'Reciprocal of the largest curvature; it must stay above the installed minimum bend radius. Shown as 10 000 m for a cable lying flat.' },
        { key: 'reaction', label: 'Cable reaction on the support', value: round(R / 1e3, 5), unit: 'kN', explain: 'Vertical force the cable puts on the separation support: the weight of the lifted length less what the touchdown points carry.' },
        { key: 'tdReaction', label: 'Peak seabed reaction at touchdown', value: round(tdR / 1e3, 4), unit: 'kN/m', explain: 'Concentration of seabed contact force where the cable lands either side of the crossing; compare with the bearing and crush resistance.' },
        { key: 'loadAsset', label: 'Load on the crossed asset', value: round(loadAsset / 1e3), unit: 'kN', limit: inp.loadAllow, sense: 'max', explain: 'Cable reaction plus the share of mattress and rock bearing on the crossed asset.' },
        { key: 'rock', label: 'Rock for the crossing berm(s)', value: round(rockT), unit: 't', explain: 'Trapezoidal berm over the protected length, all crossings; add the placement tolerance of the fall-pipe vessel.' },
        { key: 'mattresses', label: 'Concrete mattresses', value: nMat * inp.nCross, unit: '-', explain: 'Separation mattresses under the cable plus any laid over it, all crossings.' },
        { key: 'd50Req', label: 'Rock size required for stability', value: round(stab ? stab.requiredD50 : 0, 4), unit: 'm', explain: 'D50 at which the design bed shear stress equals the critical Shields stress of the rock.' },
        { key: 'dThetaOwn', label: 'Temperature rise imposed on this cable', value: round(dThOwn, 4), unit: 'K', explain: 'Warming of this cable at the crossing by the heat of the crossed asset (parallel-source upper bound).' },
        { key: 'derating', label: 'Rating factor at the crossing', value: round(DF, 5), unit: '-', explain: 'Fraction of the continuous rating left at the crossing: √(1 − Δθ/Δθmax). Below 1 the crossing is the thermal hotspot of the route.' },
        { key: 'dThetaExt', label: 'Temperature rise imposed on the crossed asset', value: round(dThExt, 4), unit: 'K', explain: 'Warming of the crossed asset by the losses of this cable.' },
        { key: 'spacingReq', label: 'Required spacing to a parallel cable', value: round(spacingReq), unit: 'm', explain: 'Larger of the minimum spacing and the depth multiple needed for a repair bight.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Cable profile over the crossing', x: { label: 'Distance along the cable from the crossing', unit: 'm' }, y: { label: 'Elevation of the cable underside above the seabed', unit: 'm' },
          series: [{ name: 'Cable', x: xp, y: yp }, { name: 'Support crest', x: [-inp.supportW / 2 - 0.01, inp.supportW / 2 + 0.01], y: [hSup, hSup], dash: true }, { name: 'Top of crossed asset', x: [-inp.crossedD / 2, inp.crossedD / 2], y: [inp.topExisting, inp.topExisting], dash: true }],
          explain: 'Installed shape: the cable rides over the support and lands on the seabed either side. The distance between the landing points is the lifted length; the gap to the dashed asset line is the separation.' },
        { type: 'line', title: 'Curvature along the cable', x: { label: 'Distance along the cable from the crossing', unit: 'm' }, y: { label: 'Curvature (sagging positive)', unit: '1/m' },
          series: [{ name: 'Curvature', x: xp, y: mp }, { name: 'Allowable', x: [xp[0], xp[xp.length - 1]], y: [1 / inp.mbrAllow, 1 / inp.mbrAllow], dash: true }, { name: 'Allowable (hogging)', x: [xp[0], xp[xp.length - 1]], y: [-1 / inp.mbrAllow, -1 / inp.mbrAllow], dash: true }],
          explain: 'Hogging peaks at the support edges and sagging peaks near the touchdown points. Both must stay inside the dashed band; a wider or lower support flattens the hogging peak.' },
        { type: 'line', title: 'Mutual heating against separation', x: { label: 'Vertical separation', unit: 'm' }, y: { label: 'Temperature rise on this cable', unit: 'K' },
          series: [{ name: 'Rise from the crossed asset', x: seps, y: dThS }, { name: 'As designed', x: [sep.separation, sep.separation], y: [0, Math.max(...dThS, dThOwn, 0.1)], dash: true }],
          explain: 'Temperature rise imposed on this cable as the separation is increased. The curve is flat: separation helps little, so a hot crossing is mitigated by lower thermal resistivity cover or a larger conductor on the crossing length.' },
      ],
      tables: [
        { title: 'Support reactions', explain: 'One row per support point of the separation layer: position along the cable and the upward reaction. Inactive points have lifted off.',
          columns: [{ key: 'x', label: 'Position', unit: 'm' }, { key: 'R', label: 'Reaction', unit: 'kN' }, { key: 'active', label: 'In contact' }],
          rows: beam.reactions.map((r) => ({ x: round(r.x - L / 2), R: round(r.R / 1e3), active: r.active ? 'yes' : 'no' })) },
        { title: 'Crossing quantities', explain: 'Material take-off for all crossings of this type, for the estimate and the installation plan.',
          columns: [{ key: 'item', label: 'Item' }, { key: 'v', label: 'Quantity' }, { key: 'u', label: 'Unit' }],
          rows: [{ item: 'Rock', v: round(rockT), u: 't' }, { item: 'Rock volume', v: round(rock ? rock.volume * inp.nCross : 0), u: 'm³' }, { item: 'Berm height / base width', v: rock ? `${fmt(rock.height)} / ${fmt(rock.baseWidth)}` : '—', u: 'm' }, { item: 'Mattresses', v: nMat * inp.nCross, u: '-' }, { item: 'Concrete', v: round(concrete), u: 'm³' }, { item: 'Vessel time', v: round(days), u: 'days' }] },
      ],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        flat ? `The crossed asset is buried and the separation layer does not rise above the seabed: the cable lies flat and the soil cover gives a separation of ${fmt(sep.separation)} m after ${inp.settlement} m of settlement.` : `The cable is lifted over ${fmt(lifted)} m and rides ${fmt(beam.deflection[iMid])} m above the seabed; separation ${fmt(sep.separation)} m after ${inp.settlement} m of settlement.`,
        `Thermal: +${fmt(dThOwn)} K on this cable (rating factor ${fmt(DF)}), +${fmt(dThExt)} K on the crossed asset; parallel-source upper bound.`,
        `Parallel circuits: ${inp.spacing} m provided against ${fmt(spacingReq)} m required in ${inp.depth} m of water.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}. ${flat ? 'The cable lies flat over the buried asset' : `Lifted length ${fmt(lifted)} m, tightest radius ${fmt(1 / Math.max(kMax, 1e-9))} m against ${inp.mbrAllow} m`}, separation ${fmt(sep.separation)} m, load on the crossed asset ${fmt(loadAsset / 1e3)} kN, rating factor ${fmt(DF)}.`,
          gov.uc > 1 ? (gov.name.startsWith('Bending') ? 'Widen the support crest or grade its edges with a pre-lay rock ramp so that the cable leaves the support tangentially; a lower separation layer also reduces curvature.' : gov.name.startsWith('Vertical') ? `Thicken the separation layer to at least ${fmt(inp.sepMin + inp.settlement)} m above the crossed asset or found it on a stiffer base to cut settlement.` : gov.name.startsWith('Free span') || gov.name.startsWith('Lifted') ? `Extend the berm to cover the full lifted length of ${fmt(lifted)} m with overrun, or place grout bags under the spans before the rock.` : gov.name.startsWith('Load on') ? 'Bridge the crossed asset: support the cable on mattresses or sleepers founded either side of it so that the reaction bypasses the asset.' : gov.name.startsWith('Rock') ? `Use rock with D50 of at least ${fmt(stab.requiredD50)} m as armour over a finer filter layer.` : gov.name.startsWith('Load current') || gov.name.startsWith('Temperature') ? 'Use a larger conductor or a low-resistivity cover on the crossing length, or increase the separation with a thicker pre-lay berm; confirm with the exact crossing solution of IEC 60287-3-3.' : gov.name.startsWith('Spacing') ? `Open the spacing to ${fmt(spacingReq)} m or agree a reduced repair method (in-line joint, shallow-water lay-down) with both owners.` : 'Re-route to cross nearer to a right angle.')
            : 'Fix the support dimensions, berm length and tolerances in the crossing agreement, survey the crossed asset immediately before the lay, and verify the as-built separation and cover by post-lay survey.',
          `${STD_CABLE}; ICPC Recommendations No. 2 and No. 3; IEC 60287-3-3; CIRIA C683 (2007)`),
        adv('economic', 'info', `${inp.nCross} crossing(s) need ${fmt(rockT)} t of rock, ${nMat * inp.nCross} mattress(es) (${fmt(concrete)} m³ of concrete) and about ${fmt(days)} vessel day(s).`,
          `Each metre taken off the protected length saves ${fmt(rock ? rock.tonnes / inp.bermL : 0)} t of rock, but the length cannot go below the lifted length of ${fmt(lifted)} m plus overrun. Crossings are priced per event: route selection that removes a crossing saves the full ${fmt(days / inp.nCross)} vessel day(s), the crossing agreement and a permanent derating or hotspot.`,
          'AACE RP 18R-97 (quantity-based estimating); whole-life costing ISO 15686-5',
          `Rock and concrete ≈ ${fmt(co2)} t CO₂e (screening ${CO2_ROCK} t CO₂e/t rock, ${CO2_CONCRETE} t CO₂e/m³ concrete). ${vesselCo2Text(days)} The berm adds ${fmt(rock ? rock.baseWidth * inp.bermL * inp.nCross : nMat * inp.nCross * inp.matL * inp.matW)} m² of hard substrate to the seabed — a permanent footprint to declare in the environmental statement.`),
      ],
      publish: { 'cable.crossingSeparation': sep.separation, 'cable.crossingDerating': DF, 'cable.crossingRock': rockT * 1e3 },
      quantities: { rock_t: round(rockT), concrete_m3: round(concrete), vessel_days: { crossings: round(days) } },
    };
  },
  verify: [
    // Heavy beam on a rigid floor lifted by a knife edge to height h (classical lift-off problem): with x from the touchdown point,
    // y = (w/EI)(L·x³/18 − x⁴/24), zero slope at the support by symmetry → touchdown shear V0 = wL/3 and
    // h = w·L⁴/(72·EI) → L = (72·EI·h/w)^¼; support reaction R = 2wL − 2V0 = 4wL/3; moment at the support M = wL²/6.
    // w = (63 − 1025·π·0.16²/4)·9.80665 = 415.7 N/m, EI = 12 000 N·m², h = 1.0 m:
    // L = (72 × 12000/415.7)^¼ = 2078.4^¼ = 6.752 m → lifted length 2L = 13.50 m; R = 4 × 415.7 × 6.752/3 = 3.742 kN;
    // curvature = M/EI = 415.7 × 6.752²/(6 × 12000) = 0.2632 1/m.
    { name: 'Cable lifted off a stiff seabed by a knife-edge support: L = (72·EI·h/w)^¼, R = 4wL/3, κ = wL²/(6·EI)',
      source: 'Classical beam lift-off solution (e.g. Timoshenko, Strength of Materials I, beam resting on a rigid plane); hand calculation',
      inputs: { topExisting: 0.7, sepT: 0.3, supportW: 0, kn: 1e6, nElem: 800, modelL: 40, settlement: 0 },
      expect: [{ kpi: 'liftedLength', value: 13.50, tol: 0.02 }, { kpi: 'reaction', value: 3.742, tol: 0.02 }, { kpi: 'curvMax', value: 0.2632, tol: 0.02 }, { kpi: 'separation', value: 0.3, tol: 0.01 }] },
    // Mutual heating by the image method: Δθ = ρ·W/(2π)·ln(d′/d). Knife-edge support (separation = 0.3 m exactly), rock cover 0.6 m, cable OD 0.16 m → own axis at 0.68 m;
    // crossed asset axis at 0.6 + 0.16 + 0.3 + 0.25 = 1.31 m (settlement 0, separation 0.3 m, asset diameter 0.5 m).
    // d = 0.63 m, d′ = 1.99 m; ρ = 1.0, W = 30: Δθ = 30/(2π) × ln(3.1587) = 4.7746 × 1.15022 = 5.492 K;
    // rating factor = sqrt(1 − 5.492/78) = 0.9642.
    { name: 'Mutual heating Δθ = ρ·W/(2π)·ln(d′/d) and the rating factor √(1 − Δθ/Δθmax)',
      source: 'IEC 60287-2-1 §4.2.3.1 (Kennelly image method); IEC 60287-3-3 (rating reduction form); hand calculation',
      inputs: { settlement: 0, supportW: 0 },
      expect: [{ kpi: 'dThetaOwn', value: 5.492, tol: 0.01 }, { kpi: 'derating', value: 0.9642, tol: 0.002 }] },
  ],
};
