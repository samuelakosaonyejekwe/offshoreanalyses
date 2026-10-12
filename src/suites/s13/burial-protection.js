// Cable burial and protection: risk-based depth of lowering against anchors (CBRA structure) and fishing gear,
// burial protection index, trenching-tool progress and feasibility, seabed mobility (sandwave) exposure, rock-berm
// and mattress quantities for unburied lengths, and the thermal consequence of the burial depth.
import { burialProtectionIndex, anchorPenetration, cbraStrikeProbability, depthOfLowering, fishingGearPenetration, GENERIC_FISHING_PENETRATION, trenchingRate, ploughTowForce,
  naturalInfillTime, rockBerm, mattressCount, burialThermalRecord, sandwaveExposure } from '../../engine/xsection/burial.js';
import { polylineLength, adv, check, worst, sev, fmt, round, run, need, linspace, CONFIRM, CO2_ROCK, vesselCo2Text, MGO_T_PER_DAY } from './_common.js';

const KN = 0.514444;
const DEFAULT_VESSELS = [
  { name: 'Fishing < 500 t', crossings: 2500, speed: 6, hazard: 600, pInc: 0.05, anchor: 400 },
  { name: 'Coaster 500–5 000 t', crossings: 900, speed: 10, hazard: 900, pInc: 0.04, anchor: 1500 },
  { name: 'Cargo 5 000–30 000 t', crossings: 500, speed: 12, hazard: 1200, pInc: 0.03, anchor: 5000 },
  { name: 'Tanker / bulk > 30 000 t', crossings: 120, speed: 13, hazard: 1500, pInc: 0.03, anchor: 12000 },
];
const SOIL_KEY = { sand: 'sand', mud: 'mud', stiffClay: 'stiffClay', gravel: 'gravel' };

export default {
  id: 's13-burial-protection', suite: 13, title: 'Cable burial and protection',
  summary: 'Required depth of lowering from anchor-strike and fishing risk, burial feasibility and progress of the trenching tool, exposure by seabed mobility, and rock or mattress protection quantities.',
  covers: [7], fidelity: 'screening',
  standards: ['Carbon Trust (2015) Cable Burial Risk Assessment Methodology (CBRA)', 'DNV-RP-0360 (2016) §5 (cable protection and burial)', 'Mole, Featherstone & Winter (1997) / Allan (1998) Burial Protection Index', 'CIRIA C683 The Rock Manual (2007)', 'IEC 60287-2-1 (external thermal resistance of a buried cable)'],
  models: ['Risk-based cable protection: probabilistic anchor-strike rate (CBRA)', 'Anchor fluke penetration geometry', 'Fishing-gear penetration table', 'Burial protection index', 'Jetting / cutting specific-energy progress model', 'Plough tow force (Cathie & Wintgens form; Rankine passive wedge)', 'Sandwave migration exposure (kinematic bedform)', 'Sediment-continuity natural backfill', 'Rock-berm and mattress take-off', 'Kennelly external thermal resistance'],
  limits: `Screening level: CBRA modifiers, anchor fluke factors, fishing-gear penetration depths and tool energy factors are overridable defaults recalled from the documents — confirm against the governing edition, AIS traffic data and tool trials. The trenching model is power-limited steady progress, not a jet-fluidisation or plough-dynamics simulation; scour around structures and liquefaction are not modelled.`,
  publishes: ['cable.burialDepth', 'cable.strikeRate', 'quantities.rock_t'],
  searchTerms: 'cable burial risk assessment depth of lowering anchor fishing trenching',
  inputs: [
    { key: 'route', label: 'Imported route (optional)', type: 'geometry', kinds: ['polyline'], default: null, group: 'Route section', role: 'geometry', help: 'Route polyline (CSV / GeoJSON / DXF); its length replaces the section length. Bathymetry grids and sediment surfaces are evaluated per section — enter the governing section here.' },
    { key: 'L', label: 'Section length', unit: 'km', type: 'number', default: 25, min: 0.05, max: 800, group: 'Route section', role: 'geometry', bind: 'cable.length' },
    { key: 'OD', label: 'Cable overall diameter', unit: 'mm', type: 'number', default: 150, min: 15, max: 400, group: 'Route section', role: 'geometry', bind: 'cable.OD' },
    { key: 'depthLow', label: 'Planned depth of lowering (to top of cable)', unit: 'm', type: 'number', default: 1.5, min: 0, max: 6, group: 'Route section', role: 'ic', design: true, help: 'Initial trench condition: target burial depth below the mean seabed.' },
    { key: 'soil', label: 'Seabed soil', type: 'select', default: 'sand', group: 'Soil', role: 'soil', options: [{ value: 'sand', label: 'Sand' }, { value: 'mud', label: 'Soft mud / clay' }, { value: 'stiffClay', label: 'Stiff clay' }, { value: 'gravel', label: 'Gravel' }] },
    { key: 'su', label: 'Undrained shear strength (clays)', unit: 'kPa', type: 'number', default: 20, min: 0.5, max: 400, group: 'Soil', role: 'soil', bind: 'site.su0', showIf: { key: 'soil', in: ['mud', 'stiffClay'] }, uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'phi', label: 'Friction angle (granular)', unit: 'deg', type: 'number', default: 33, min: 20, max: 45, group: 'Soil', role: 'soil', bind: 'site.phiDeg', showIf: { key: 'soil', in: ['sand', 'gravel'] } },
    { key: 'gammaSub', label: 'Submerged unit weight', unit: 'kN/m³', type: 'number', default: 9.5, min: 2, max: 13, group: 'Soil', role: 'soil', bind: 'site.gammaSub' },
    { key: 'rhoSoil', label: 'Native soil thermal resistivity', unit: 'K·m/W', type: 'number', default: 0.7, min: 0.2, max: 3, group: 'Soil', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'rhoBackfill', label: 'Backfill thermal resistivity', unit: 'K·m/W', type: 'number', default: 0.9, min: 0.2, max: 4, group: 'Soil', role: 'soil', help: 'Jetted or ploughed backfill is looser than native soil and usually has a higher resistivity.' },
    { key: 'vessels', label: 'Vessel traffic crossing the section', type: 'table', group: 'Anchor threat (CBRA)', role: 'load',
      columns: [{ key: 'name', label: 'Vessel class', type: 'text' }, { key: 'crossings', label: 'Crossings per year', type: 'number' }, { key: 'speed', label: 'Speed', unit: 'kn', type: 'number' }, { key: 'hazard', label: 'Hazard distance', unit: 'm', type: 'number' }, { key: 'pInc', label: 'Incident probability per vessel-year', type: 'number' }, { key: 'anchor', label: 'Anchor mass', unit: 'kg', type: 'number' }],
      default: DEFAULT_VESSELS.map((r) => ({ ...r })), help: 'From AIS data. Hazard distance = track length within which a dragged anchor can reach the cable. Indicative demonstration values.' },
    { key: 'flukeCoef', label: 'Fluke length coefficient', unit: 'm/kg^⅓', type: 'number', default: 0.1, min: 0.03, max: 0.3, group: 'Anchor threat (CBRA)', role: 'criteria', help: `Fluke length = coefficient × anchor mass^⅓ (geometric similarity, fitted to stockless-anchor catalogues). ${CONFIRM}` },
    { key: 'flukeAngle', label: 'Fluke opening angle', unit: 'deg', type: 'number', default: 45, min: 25, max: 60, group: 'Anchor threat (CBRA)', role: 'criteria', help: CONFIRM },
    { key: 'soilFactor', label: 'Soil penetration factor', unit: '-', type: 'number', default: 1.0, min: 0.3, max: 5, group: 'Anchor threat (CBRA)', role: 'criteria', help: `About 1 in sand and stiff clay, larger in soft clay where the anchor sinks. ${CONFIRM}` },
    { key: 'pTraffic', label: 'Traffic modifier P_traffic', unit: '-', type: 'number', default: 1, min: 0, max: 1, group: 'Anchor threat (CBRA)', role: 'criteria', help: CONFIRM },
    { key: 'pWd', label: 'Water-depth modifier P_wd', unit: '-', type: 'number', default: 0.5, min: 0, max: 1, group: 'Anchor threat (CBRA)', role: 'criteria', help: `Feasibility of anchoring at the water depth of the section. ${CONFIRM}` },
    { key: 'pTarget', label: 'Acceptable strike probability', unit: '1/yr', type: 'number', default: 1e-4, min: 1e-7, max: 1, group: 'Anchor threat (CBRA)', role: 'criteria', help: 'Project risk criterion per section (commonly a 1-in-10 000-year to 1-in-100 000-year event).' },
    { key: 'gear', label: 'Fishing gear', type: 'select', default: 'otterTrawl', group: 'Fishing threat', role: 'load', options: [{ value: 'otterTrawl', label: 'Otter trawl' }, { value: 'beamTrawl', label: 'Beam trawl' }, { value: 'scallopDredge', label: 'Scallop dredge' }, { value: 'hydraulicDredge', label: 'Hydraulic dredge' }, { value: 'custom', label: 'Custom penetration' }] },
    { key: 'gearPen', label: 'Custom gear penetration', unit: 'm', type: 'number', default: 0.3, min: 0, max: 2, group: 'Fishing threat', role: 'criteria', showIf: { key: 'gear', in: ['custom'] } },
    { key: 'gearSF', label: 'Safety factor on gear penetration', unit: '-', type: 'number', default: 2, min: 1, max: 5, group: 'Fishing threat', role: 'criteria', help: `Generic penetration table (indicative, below 0.3 m). ${CONFIRM}` },
    { key: 'margin', label: 'Margin on the governing penetration', unit: 'm', type: 'number', default: 0.3, min: 0, max: 2, group: 'Fishing threat', role: 'criteria' },
    { key: 'bpiRef', label: 'Burial depth giving BPI = 1 in this soil', unit: 'm', type: 'number', default: 0.6, min: 0.1, max: 5, group: 'Fishing threat', role: 'criteria', help: `Read from the published BPI chart for the soil strength (larger in soft clay). ${CONFIRM}` },
    { key: 'bpiTarget', label: 'Target burial protection index', unit: '-', type: 'number', default: 2, min: 0.5, max: 3, group: 'Fishing threat', role: 'criteria' },
    { key: 'tool', label: 'Burial tool', type: 'select', default: 'jet', group: 'Burial tool', role: 'operation', options: [{ value: 'jet', label: 'Jetting ROV / sled' }, { value: 'cutter', label: 'Mechanical chain cutter' }, { value: 'plough', label: 'Cable plough' }] },
    { key: 'power', label: 'Installed jetting / cutting power', unit: 'kW', type: 'number', default: 600, min: 20, max: 5000, group: 'Burial tool', role: 'operation', design: true, showIf: { key: 'tool', in: ['jet', 'cutter'] } },
    { key: 'eff', label: 'Fraction of power reaching the soil', unit: '-', type: 'number', default: 0.3, min: 0.02, max: 1, group: 'Burial tool', role: 'operation', showIf: { key: 'tool', in: ['jet', 'cutter'] }, uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'eSpec', label: 'Specific energy of the soil', unit: 'kJ/m³', type: 'number', default: 2500, min: 5, max: 50000, group: 'Burial tool', role: 'soil', showIf: { key: 'tool', in: ['jet', 'cutter'] }, uncertain: { dist: 'lognormal', cov: 0.4 }, help: 'Tool energy delivered to the soil per cubic metre of trench (fluidising, cutting and transporting the spoil): of the order of 1 000–5 000 for jetting in sand and higher in stiff clay or weak rock. No normative value — calibrate against trenching trials and installation records.' },
    { key: 'trenchW', label: 'Trench width', unit: 'm', type: 'number', default: 0.5, min: 0.1, max: 5, group: 'Burial tool', role: 'geometry' },
    { key: 'toolMax', label: 'Maximum burial depth of the tool', unit: 'm', type: 'number', default: 3, min: 0.3, max: 10, group: 'Burial tool', role: 'operation' },
    { key: 'ploughW', label: 'Plough submerged weight', unit: 'kN', type: 'number', default: 250, min: 10, max: 2000, group: 'Burial tool', role: 'operation', showIf: { key: 'tool', in: ['plough'] } },
    { key: 'ploughV', label: 'Plough speed', unit: 'm/h', type: 'number', default: 300, min: 10, max: 2000, group: 'Burial tool', role: 'operation', showIf: { key: 'tool', in: ['plough'] } },
    { key: 'Cw', label: 'Plough friction coefficient Cw', unit: '-', type: 'number', default: 0.4, min: 0.1, max: 1, group: 'Burial tool', role: 'soil', showIf: { key: 'tool', in: ['plough'] } },
    { key: 'Cs', label: 'Plough static coefficient Cs (sand) / Cc (clay)', unit: '-', type: 'number', default: 10, min: 0.5, max: 60, group: 'Burial tool', role: 'soil', showIf: { key: 'tool', in: ['plough'] }, help: 'Cathie & Wintgens (2001) form; plough- and soil-specific.' },
    { key: 'CdP', label: 'Plough rate coefficient Cd (sand)', unit: 'kN·s/m³', type: 'number', default: 100, min: 0, max: 2000, group: 'Burial tool', role: 'soil', showIf: { key: 'tool', in: ['plough'] } },
    { key: 'bollard', label: 'Available tow force', unit: 'kN', type: 'number', default: 1200, min: 50, max: 5000, group: 'Burial tool', role: 'operation', showIf: { key: 'tool', in: ['plough'] } },
    { key: 'operability', label: 'Weather / tool operability', unit: '-', type: 'number', default: 0.7, min: 0.1, max: 1, group: 'Burial tool', role: 'operation' },
    { key: 'qSed', label: 'Sediment transport rate across the trench', unit: 'm²/day', type: 'number', default: 0.05, min: 0, max: 20, group: 'Seabed mobility', role: 'environment', help: 'Volume per metre width per day from a sediment-transport study; 0 = no natural backfill.' },
    { key: 'swH', label: 'Sandwave height (trough to crest)', unit: 'm', type: 'number', default: 2, min: 0, max: 20, group: 'Seabed mobility', role: 'environment' },
    { key: 'swL', label: 'Sandwave wavelength', unit: 'm', type: 'number', default: 200, min: 5, max: 2000, group: 'Seabed mobility', role: 'environment' },
    { key: 'swC', label: 'Sandwave migration rate', unit: 'm/yr', type: 'number', default: 2, min: 0, max: 100, group: 'Seabed mobility', role: 'environment', uncertain: { dist: 'lognormal', cov: 0.4 } },
    { key: 'swPhase', label: 'Position on the sandwave at installation', type: 'select', default: 'flank', group: 'Seabed mobility', role: 'ic', options: [{ value: 'crest', label: 'Under a crest (worst case)' }, { value: 'flank', label: 'Mid flank' }, { value: 'trough', label: 'In a trough (pre-swept / best case)' }] },
    { key: 'life', label: 'Design life', unit: 'years', type: 'number', default: 30, min: 1, max: 60, group: 'Seabed mobility', role: 'criteria', bind: 'life.designYears' },
    { key: 'rockLen', label: 'Length protected by rock (crossings, unburied)', unit: 'm', type: 'number', default: 300, min: 0, max: 1e5, group: 'Rock and mattresses', role: 'geometry' },
    { key: 'rockCover', label: 'Rock cover above the cable', unit: 'm', type: 'number', default: 0.6, min: 0.1, max: 3, group: 'Rock and mattresses', role: 'geometry', design: true },
    { key: 'crest', label: 'Berm crest width', unit: 'm', type: 'number', default: 1, min: 0.2, max: 10, group: 'Rock and mattresses', role: 'geometry' },
    { key: 'bermSlope', label: 'Berm side slope (horizontal : 1 vertical)', unit: '-', type: 'number', default: 3, min: 1, max: 6, group: 'Rock and mattresses', role: 'geometry' },
    { key: 'matLen', label: 'Length protected by mattresses', unit: 'm', type: 'number', default: 60, min: 0, max: 5000, group: 'Rock and mattresses', role: 'geometry' },
  ],
  compute(inp) {
    const D = inp.OD / 1e3, L = polylineLength(inp.route) ?? inp.L * 1e3;
    const vs = (inp.vessels || []).filter((v) => v.crossings > 0);
    need(vs.length > 0, 'The traffic table is empty: add at least one vessel class with crossings per year (or a class with a negligible count if there is no traffic).');
    for (const v of vs) need(v.speed > 0 && v.hazard > 0 && v.pInc >= 0 && v.anchor > 0, `Vessel class '${v.name}': speed, hazard distance and anchor mass must be positive.`);
    const classes = vs.map((v) => ({ name: v.name, crossings: v.crossings, speed: v.speed * KN, hazardDistance: v.hazard, pIncident: v.pInc,
      penetration: anchorPenetration({ anchorMass: v.anchor, flukeCoefficient: inp.flukeCoef, flukeAngleDeg: inp.flukeAngle, soilFactor: inp.soilFactor }).penetration }));
    const fish = inp.gear === 'custom' ? { penetration: inp.gearPen, design: inp.gearPen * inp.gearSF } : fishingGearPenetration({ gear: inp.gear, soil: SOIL_KEY[inp.soil], safetyFactor: inp.gearSF, table: GENERIC_FISHING_PENETRATION });
    const mods = { pTraffic: inp.pTraffic, pWd: inp.pWd };
    const at = run(() => cbraStrikeProbability({ classes, burialDepth: inp.depthLow, ...mods }), 'Anchor-strike probability');
    const dol = run(() => depthOfLowering({ classes, targetProbability: inp.pTarget, fishingPenetration: fish.design, margin: inp.margin, ...mods }), 'Depth of lowering');
    const bpi = burialProtectionIndex({ referenceDepth: inp.bpiRef, burialDepth: inp.depthLow, targetBpi: inp.bpiTarget });
    // trenching progress
    let rate, towUc = 0, tow = null;
    const trenchDepth = inp.depthLow + D;
    if (inp.tool === 'plough') {
      const clay = inp.soil === 'mud' || inp.soil === 'stiffClay';
      tow = ploughTowForce({ soil: clay ? 'clay' : 'sand', depth: trenchDepth, ploughWeight: inp.ploughW * 1e3, gammaSub: inp.gammaSub * 1e3, speed: inp.ploughV / 3600, Cw: inp.Cw, ...(clay ? { Cc: inp.Cs, su: inp.su * 1e3 } : { Cs: inp.Cs, Cd: inp.CdP * 1e3, phiDeg: inp.phi }), shareWidth: inp.trenchW, mu: inp.Cw });
      rate = inp.ploughV; towUc = tow.empirical / (inp.bollard * 1e3);
    } else {
      rate = run(() => trenchingRate({ power: inp.power * 1e3, efficiency: inp.eff, specificEnergy: inp.eSpec * 1e3, trenchWidth: inp.trenchW, trenchDepth }), 'Trenching rate').ratePerHour;
    }
    const burialDays = L / (rate * 24) / inp.operability;
    const trenchArea = inp.trenchW * trenchDepth;
    const infill = inp.qSed > 0 ? naturalInfillTime({ trenchArea, transportRate: inp.qSed / 86400 }).timeDays : Infinity;
    // sandwave exposure: cover must not fall below the fishing design penetration within the design life
    const phase = { crest: 0, flank: Math.PI / 2, trough: Math.PI }[inp.swPhase];
    const sw = inp.swH > 0 && inp.swC > 0 ? sandwaveExposure({ height: inp.swH, wavelength: inp.swL, migrationRate: inp.swC, initialCover: inp.depthLow, targetCover: fish.design, phase }) : { time: Infinity, minimumCover: inp.depthLow, exposed: false, period: Infinity };
    const coverAt = (t) => (inp.swH > 0 ? inp.depthLow + inp.swH / 2 * (Math.cos(phase - 2 * Math.PI * inp.swC * t / inp.swL) - Math.cos(phase)) : inp.depthLow);
    const ts = linspace(0, inp.life, 80), covers = ts.map(coverAt), minCoverLife = Math.min(...covers);
    const th = burialThermalRecord({ depthToTop: Math.max(inp.depthLow, 0.01), outerDiameter: D, soilResistivity: inp.rhoSoil, backfillResistivity: inp.rhoBackfill });
    const berm = inp.rockLen > 0 ? rockBerm({ productDiameter: D, cover: inp.rockCover, crestWidth: inp.crest, slope: inp.bermSlope, length: inp.rockLen, overdump: 0.15, requiredAbsorptionCover: fish.design }) : { tonnes: 0, height: 0, baseWidth: 0, coverUtilisation: 0, volume: 0 };
    const mats = inp.matLen > 0 ? mattressCount({ length: inp.matLen, mattressLength: 6, overlap: 0, endOverrun: 3, coverWidth: 3, mattressWidth: 3, mattressMass: 9000 }) : { count: 0, totalMass: 0 };
    const rockDays = berm.tonnes / 8000 + (berm.tonnes > 0 ? 2 : 0);    // fall-pipe vessel: ≈ 8 000 t placed per day + transit/loading (screening)
    const req = Math.max(dol.recommendedDepth, bpi.requiredDepth ?? 0);
    const checks = [
      check('Depth of lowering against the risk-based requirement', dol.recommendedDepth, Math.max(inp.depthLow, 1e-6), 'm', 'Carbon Trust CBRA (2015) §5 (depth of lowering specification)', 'Demand is the smallest depth meeting the acceptable anchor-strike probability and the factored fishing penetration, plus margin; capacity is the planned depth. A shallower cable is hooked or struck more often than the project accepts.'),
      check('Anchor-strike probability at the planned depth', at.pStrike, inp.pTarget, '1/yr', 'Carbon Trust CBRA (2015) probabilistic assessment', 'Annual probability that a vessel in distress drops an anchor that penetrates below the planned burial depth on this section, against the acceptable value.'),
      check('Fishing-gear penetration', fish.design, Math.max(inp.depthLow, 1e-6), 'm', 'DNV-RP-0360 §5; fishing-gear interaction studies (generic table)', 'Factored penetration of the gear in this soil against the planned cover. Gear reaching the cable snags it or strips the serving.'),
      check('Burial protection index', inp.bpiTarget, Math.max(bpi.bpi, 1e-6), '-', 'Mole et al. (1997); Allan (1998)', 'Target index against the index achieved by the planned depth in this soil (BPI 1 fishing gear, 2 anchors to about 2 t, 3 larger anchors).'),
      check('Tool burial capability', inp.depthLow, inp.toolMax, 'm', 'DNV-RP-0360 §5 (burial assessment study); tool specification', 'Planned depth against the maximum depth the tool can reach. Beyond it the section needs a second pass, a different tool or rock.'),
      ...(inp.tool === 'plough' ? [check('Plough tow force', tow.empirical / 1e3, inp.bollard, 'kN', 'Cathie & Wintgens (2001) OTC 13145', 'Predicted steady tow force at the planned depth and speed against the available tow force of the vessel or winch.')] : []),
      check('Cover retained under seabed mobility', fish.design, Math.max(minCoverLife, 1e-6), 'm', 'DNV-RP-0360 §5 (mobile seabed, reference seabed level)', 'Demand is the cover needed against fishing gear, capacity the smallest cover reached within the design life as sandwaves migrate past. Burial should refer to the non-mobile reference level.'),
      ...(inp.rockLen > 0 ? [check('Rock cover for gear absorption', berm.coverUtilisation, 1, '-', 'CIRIA C683 (2007); DNV-RP-0360 §5', 'Required rock thickness over the cable against the cover provided on the unburied length.')] : []),
    ];
    const gov = worst(checks);
    const ds = linspace(0, Math.max(...classes.map((c) => c.penetration), inp.depthLow) * 1.2 + 0.2, 60);
    const pOfD = ds.map((d) => cbraStrikeProbability({ classes, burialDepth: d, ...mods }).pStrike);
    const t4OfD = ds.map((d) => burialThermalRecord({ depthToTop: Math.max(d, 0.02), outerDiameter: D, soilResistivity: inp.rhoSoil, backfillResistivity: inp.rhoBackfill }).T4);
    const rockT = berm.tonnes, extra = Math.max(0, req - inp.depthLow);
    return {
      kpis: [
        { key: 'recommendedDepth', label: 'Required depth of lowering', value: round(dol.recommendedDepth), unit: 'm', explain: `Smallest burial depth meeting the acceptable strike probability and the fishing penetration, with margin. Governed by ${dol.governing}.` },
        { key: 'pStrike', label: 'Anchor-strike probability at the planned depth', value: round(at.pStrike, 4), unit: '1/yr', limit: inp.pTarget, sense: 'max', explain: 'Annual frequency of an anchor reaching the cable on this section at the planned depth of lowering.' },
        { key: 'returnPeriod', label: 'Strike return period', value: round(Math.min(at.returnPeriod, 1e9)), unit: 'years', explain: 'Reciprocal of the strike probability (10⁹ = no vessel class can reach the cable).' },
        { key: 'pUnburied', label: 'Strike probability if left unburied', value: round(dol.pStrikeUnburied, 4), unit: '1/yr', explain: 'Frequency with the cable on the seabed: the risk that burial removes.' },
        { key: 'fishDepth', label: 'Design fishing-gear penetration', value: round(fish.design), unit: 'm', explain: 'Tabulated penetration of the selected gear in this soil times the safety factor.' },
        { key: 'maxAnchor', label: 'Deepest anchor penetration in the traffic', value: round(Math.max(...classes.map((c) => c.penetration))), unit: 'm', explain: 'Fluke-tip depth of the largest anchor class; burying below it removes the anchor risk entirely.' },
        { key: 'bpi', label: 'Burial protection index achieved', value: round(bpi.bpi), unit: '-', limit: inp.bpiTarget, sense: 'min', explain: 'Planned depth divided by the depth giving BPI = 1 in this soil.' },
        { key: 'rate', label: 'Burial progress rate', value: round(rate), unit: 'm/h', explain: 'Steady progress of the tool at the planned depth: power-limited for jetting and cutting, the entered speed for ploughing.' },
        { key: 'burialDays', label: 'Burial duration incl. operability', value: round(burialDays), unit: 'days', explain: 'Section length over progress rate, divided by the operability.' },
        { key: 'infillDays', label: 'Natural backfill time', value: round(Math.min(infill, 1e6)), unit: 'days', explain: 'Time for sediment transport to fill the open trench (10⁶ = none). Until then the cable relies on the trench walls alone.' },
        { key: 'minCover', label: 'Minimum cover within the design life', value: round(minCoverLife), unit: 'm', explain: 'Smallest cover as sandwaves migrate past the cable. Negative means the cable is exposed and may span.' },
        { key: 'exposureTime', label: 'Time until cover falls to the fishing depth', value: round(Math.min(sw.time, 1e6)), unit: 'years', explain: 'First time the cover drops to the fishing-gear design penetration (10⁶ = never).' },
        { key: 'T4', label: 'External thermal resistance at the planned depth', value: round(th.T4), unit: 'K·m/W', explain: 'Soil thermal resistance seen by a single cable at this depth with the backfill resistivity: deeper burial protects better but derates the cable.' },
        { key: 'rockT', label: 'Rock for unburied lengths', value: round(rockT), unit: 't', explain: 'Berm tonnage including 15 % placement allowance.' },
        { key: 'mattresses', label: 'Concrete mattresses', value: mats.count, unit: '-', explain: 'Number of 6 m × 3 m mattresses for the mattress-protected length.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Anchor-strike probability against depth of lowering', x: { label: 'Depth of lowering', unit: 'm' }, y: { label: 'log₁₀ strike probability', unit: '1/yr' },
          series: [{ name: 'CBRA strike rate', x: ds, y: pOfD.map((p) => Math.log10(Math.max(p, 1e-12))) }, { name: 'Acceptable', x: [0, ds[59]], y: [Math.log10(inp.pTarget), Math.log10(inp.pTarget)], dash: true }, { name: 'Planned depth', x: [inp.depthLow, inp.depthLow], y: [-12, Math.log10(Math.max(pOfD[0], 1e-12))], dash: true }],
          explain: 'The rate falls in steps as the burial depth passes the fluke penetration of each vessel class. The required depth is where the curve first drops below the acceptable line.' },
        { type: 'bar', title: 'Strike contribution by vessel class (planned depth)', categories: at.classes.map((c) => c.name), series: [{ name: 'Contribution', y: at.classes.map((c) => c.contribution) }],
          explain: 'Which traffic drives the residual risk at the planned depth. Classes with zero contribution cannot reach the cable.' },
        { type: 'line', title: 'Cover over the design life (seabed mobility)', x: { label: 'Time', unit: 'years' }, y: { label: 'Cover over the cable', unit: 'm' },
          series: [{ name: 'Cover', x: ts, y: covers }, { name: 'Fishing design penetration', x: [0, inp.life], y: [fish.design, fish.design], dash: true }],
          explain: 'Cover as the sandwave migrates past. The crossing of the dashed line is when re-burial or rock becomes necessary; cover below zero means exposure.' },
        { type: 'line', title: 'Thermal resistance against depth of lowering', x: { label: 'Depth of lowering', unit: 'm' }, y: { label: 'T4', unit: 'K·m/W' }, series: [{ name: 'Single buried cable', x: ds, y: t4OfD }],
          explain: 'The external thermal resistance grows with the logarithm of depth: every extra metre of protection costs ampacity. Use it with the thermal-rating module.' },
      ],
      tables: [{ title: 'Vessel classes', explain: 'One row per vessel class: anchor fluke penetration, exposure hours per year over the cable, whether the anchor can reach the planned depth and its share of the strike rate.',
        columns: [{ key: 'name', label: 'Class' }, { key: 'pen', label: 'Anchor penetration', unit: 'm' }, { key: 'exp', label: 'Exposure', unit: 'h/yr' }, { key: 'thr', label: 'Reaches the cable' }, { key: 'p', label: 'Contribution', unit: '1/yr' }],
        rows: at.classes.map((c, i) => ({ name: c.name, pen: round(classes[i].penetration), exp: round(c.exposureHours), thr: c.threatens ? 'yes' : 'no', p: round(c.contribution, 3) })) }],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}). Required depth ${fmt(req)} m against ${inp.depthLow} m planned.`,
        `Trench ${inp.trenchW} m × ${fmt(trenchDepth)} m in ${inp.soil}: ${fmt(rate)} m/h${inp.tool === 'plough' ? `, tow force ${fmt(tow.empirical / 1e3)} kN (Rankine static estimate ${fmt((tow.rankine ?? 0) / 1e3)} kN)` : ''}.`,
        sw.exposed ? `Sandwave migration brings the cover down to the fishing depth after ${fmt(sw.time)} years (period ${fmt(sw.period)} years).` : 'Seabed mobility does not reduce the cover to the fishing depth.',
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}. Strike rate ${fmt(at.pStrike)} per year (1 in ${fmt(Math.min(at.returnPeriod, 1e9))} years) at ${inp.depthLow} m; required depth ${fmt(req)} m.`,
          gov.uc > 1 ? (extra > 0 ? `Lower the cable by a further ${fmt(extra)} m (to ${fmt(req)} m), or protect the high-traffic stretch with rock; where sandwaves govern, pre-sweep to the trough level before burial.` : gov.name.includes('Tool') ? 'Select a tool rated for the planned depth or split the section between burial and rock placement.' : 'Increase the rock cover or the tow / jetting capacity as indicated by the governing check.')
            : 'Confirm traffic with 12 months of AIS data, the soil units with the burial assessment study, and specify the depth of lowering relative to the non-mobile reference seabed level; verify as-buried depth by survey.',
          'Carbon Trust CBRA (2015); DNV-RP-0360 (2016) §5'),
        adv('economic', 'info', `Burial of ${fmt(L / 1e3)} km takes about ${fmt(burialDays)} days at ${fmt(rate)} m/h; rock for ${inp.rockLen} m is ${fmt(rockT)} t (≈ ${fmt(rockDays)} fall-pipe vessel days) and ${mats.count} mattresses.`,
          `Depth is paid twice: each extra 0.5 m slows the tool by ${fmt((1 - (inp.depthLow + D) / (inp.depthLow + D + 0.5)) * 100)} % and raises T4 from ${fmt(th.T4)} to ${fmt(burialThermalRecord({ depthToTop: inp.depthLow + 0.5, outerDiameter: D, soilResistivity: inp.rhoSoil, backfillResistivity: inp.rhoBackfill }).T4)} K·m/W (lower rating). Bury to the risk-based depth section by section rather than to one blanket depth.`,
          'Risk-based cable protection (Carbon Trust CBRA); AACE RP 18R-97 quantity-based estimate; whole-life costing ISO 15686-5',
          `Rock ${fmt(rockT)} t ≈ ${fmt(rockT * CO2_ROCK)} t CO₂e at the quarry (screening ${CO2_ROCK} t CO₂e/t) and a permanent hard-substrate footprint of ${fmt(berm.baseWidth * inp.rockLen)} m²; trenching disturbs about ${fmt(inp.trenchW * L)} m² of seabed that recovers. ${vesselCo2Text(burialDays + rockDays)}`),
      ],
      publish: { 'cable.burialDepth': inp.depthLow, 'cable.strikeRate': at.pStrike, 'quantities.rock_t': rockT },
      quantities: { cable_km: round(L / 1e3, 6), rock_t: round(rockT), concrete_m3: round(mats.count * 9000 / 2400), vessel_days: { trenching: round(burialDays), rockPlacement: round(rockDays) } },
    };
  },
  verify: [
    // One vessel class: 1000 crossings/yr at 5 m/s over a 1000 m hazard distance, P_incident = 0.1 per vessel-year,
    // modifiers 1: exposure per crossing 200 s = 0.05556 h → P = 1000 × 0.05556/8760 × 0.1 = 6.342e-4 per year (unburied).
    // Anchor 5000 kg, coefficient 0.1 → fluke 1.710 m × sin45° = 1.209 m. BPI = 0 / 0.6 = 0 (unburied) — not checked.
    { name: 'CBRA strike rate P = Σ N·D/(V·8760 h)·P_incident and fluke penetration L·sin(angle)',
      source: 'Carbon Trust (2015) Cable Burial Risk Assessment Methodology, probabilistic method; hand calculation',
      inputs: { depthLow: 0, pWd: 1, pTraffic: 1, swH: 0, vessels: [{ name: 'A', crossings: 1000, speed: 5 / 0.514444, hazard: 1000, pInc: 0.1, anchor: 5000 }] },
      expect: [{ kpi: 'pStrike', value: 6.342e-4, tol: 0.003 }, { kpi: 'maxAnchor', value: 1.209, tol: 0.003 }] },
    // Jetting: 600 kW × 0.3 = 180 kW into a 0.5 m × (1.5 + 0.15) m trench at 2500 kJ/m³ → v = 180e3/(2500e3 × 0.825) = 0.08727 m/s = 314.2 m/h.
    // Sandwave under a crest: H = 2 m, λ = 200 m, c = 2 m/yr, cover 1.5 → 0.4 m: t = λ/(2πc)·acos(1 − 2·1.1/2) = 15.915 × 1.6710 = 26.59 yr.
    { name: 'Specific-energy trenching rate and sandwave exposure time',
      source: 'Energy balance v = ηP/(e·A) (Teale 1965 concept); kinematic bedform migration (closed form)',
      inputs: { swPhase: 'crest' },
      expect: [{ kpi: 'rate', value: 314.2, tol: 0.003 }, { kpi: 'exposureTime', value: 26.59, tol: 0.005 }] },
  ],
};
