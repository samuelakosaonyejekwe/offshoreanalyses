// What can honestly be said about soil strength before a site investigation exists: published typical ranges for
// the kind of seabed identified, an opt-in screening profile taken from the cautious end of those ranges, and
// what investigation is normally expected for the structure in hand. Nothing here is site data. Pure functions.

const REFS = {
  ncClay: 'Normally consolidated clays have an undrained strength of about 0.20 to 0.30 times the vertical effective stress (Skempton 1957; Ladd & Foott 1974; Mesri 1975: about 0.22 as mobilised in the field).',
  softClay: 'Soft marine clays weigh about 4 to 8 kN/m³ submerged near the mudline, so the strength rises by roughly 1 to 2 kPa per metre; deep-water clays commonly show a few kilopascals at the mudline (Randolph & Gourvenec 2011, Offshore Geotechnical Engineering; Lunne, Robertson & Powell 1997).',
  sandPhi: 'Friction angle of silica sand by density: very loose below about 29°, loose 29–30°, medium dense 30–36°, dense 36–41°, very dense above 41° (Peck, Hanson & Thornburn 1974; Meyerhof 1956).',
  sandGamma: 'Saturated silica sand weighs about 8 to 11 kN/m³ submerged.',
  codes: 'ISO 19901-4 and API RP 2GEO require characteristic soil parameters from a site-specific investigation; they give no default strengths.',
  carbonate: 'ISO 19901-4 and API RP 2GEO single out carbonate soils: the grains crush and cemented layers break up during pile driving, so shaft friction can be a small fraction of that in silica sand of the same friction angle, and laboratory friction angles (often 35–45°) say little about foundation capacity.',
  opal: 'Diatom and radiolarian sediments have very high water content and low unit weight, and the correlations established for ordinary clays do not hold for them.',
};

// families of seabed classes with the same screening treatment
const FAMILY = {
  mud: 'clay', sandyMud: 'clay', clay: 'clay', sand: 'sand', gravel: 'sand', muddySand: 'siltySand',
  silt: 'silt', lithogenous: 'silt', mixed: 'mixed', volcanic: 'volcanic', biogenic: 'carbonate', calcOoze: 'carbonate', calcMud: 'carbonate', mixedOoze: 'carbonate',
  silOoze: 'opal', silMud: 'opal', rock: 'rock',
};

const PROFILES = {
  clay: {
    title: 'Soft, normally consolidated clay',
    ranges: [['Undrained shear strength at the mudline', '1 to 5 kPa'], ['Strength gain with depth', '1.2 to 2 kPa per metre (0.20 to 0.30 × the vertical effective stress)'], ['Submerged unit weight', '4 to 8 kN/m³']],
    refs: [REFS.ncClay, REFS.softClay, REFS.codes],
    values: { 'site.su0': 2000, 'site.suGrad': 1200, 'site.gammaSub': 5000 },
    chosen: 'Mudline strength 2 kPa, strength gain 1.2 kPa/m, submerged unit weight 5 kN/m³: the weak end of the ranges.',
    caution: 'The weak end is the cautious side for bearing, holding and pile capacity. It is the wrong side for anything that must be pushed or driven into the ground (skirt and anchor penetration, pile driving, trenching) and for pipeline embedment checks that rely on soil resistance: run those with the strong end as well. An over-consolidated or glacial clay is many times stronger than this profile.',
  },
  sand: {
    title: 'Silica sand of unknown density',
    ranges: [['Friction angle', '29 to 41° from loose to dense'], ['Submerged unit weight', '8 to 11 kN/m³']],
    refs: [REFS.sandPhi, REFS.sandGamma, REFS.codes],
    values: { 'site.phiDeg': 30, 'site.gammaSub': 9000 },
    chosen: 'Friction angle 30° (loose to medium dense) and submerged unit weight 9 kN/m³: the weak end of the ranges.',
    caution: 'The weak end is the cautious side for bearing and sliding capacity. Dense sand is the governing case for pile driving, skirt penetration and trenching: run those with 40° as well. If the sand is rich in shell or coral, do not use this profile — carbonate sand needs its own tests.',
  },
  siltySand: {
    title: 'Muddy (silty or clayey) sand',
    ranges: [['Friction angle', '27 to 36°: fines lower the friction angle of a sand and slow its drainage'], ['Submerged unit weight', '8 to 11 kN/m³']],
    refs: [REFS.sandPhi, REFS.sandGamma, REFS.codes],
    values: { 'site.phiDeg': 29, 'site.gammaSub': 9000 },
    chosen: 'Friction angle 29° and submerged unit weight 9 kN/m³: the weak end of the ranges.',
    caution: 'A muddy sand may respond partly undrained under wave and earthquake loading; a friction angle alone does not describe that. Dense material governs driving and penetration: check those with 36° as well.',
  },
};
const WITHHELD = {
  silt: ['Silt lies between drained and undrained behaviour: neither a clay strength line nor a sand friction angle describes it, and its response to cyclic loading has to be tested.', [REFS.codes]],
  mixed: ['A mixed or volcanic bed changes over short distances, so no single profile can stand for it.', [REFS.codes]],
  volcanic: ['A mixed or volcanic bed changes over short distances, so no single profile can stand for it.', [REFS.codes]],
  carbonate: ['Carbonate sediment: generic strength values are known to mislead here.', [REFS.carbonate, REFS.codes]],
  opal: ['Siliceous (opal-rich) sediment: generic strength values are known to mislead here.', [REFS.opal, REFS.codes]],
  rock: ['Rock: strength depends on the rock type, weathering and jointing, which only coring can show.', [REFS.codes]],
};

/** Bus keys a screening profile may set; all are cleared before a new profile is applied. */
export const SCREENING_KEYS = ['site.su0', 'site.suGrad', 'site.phiDeg', 'site.gammaSub'];
/** The wording shown wherever a screening value is used. */
export const SCREENING_FLAG = 'screening value — no site investigation yet';

/**
 * Screening soil profile for a seabed class.
 * @param {string|null} cls key of SEABED_CLASSES (null when the seabed type is unknown)
 * @param {{ carbonate?: number|null }} [o] measured carbonate content of nearby samples [%]
 * @returns {{ offered: boolean, title?: string, ranges?: string[][], refs: string[], values?: object, chosen?: string,
 *   caution?: string, why?: string }} `offered` false with the reason in `why` for rock, unknown ground and the
 *   classes whose strength generic values cannot represent
 */
export function screeningProfile(cls, { carbonate = null } = {}) {
  if (!cls) return { offered: false, why: 'The kind of seabed is not known, so there is nothing to base a profile on.', refs: [REFS.codes] };
  const fam = carbonate >= 50 ? 'carbonate' : FAMILY[cls];
  if (PROFILES[fam]) return { offered: true, ...PROFILES[fam] };
  const [why, refs] = WITHHELD[fam] || WITHHELD.mixed;
  return { offered: false, why: `${why} No screening profile is offered.`, refs };
}

/** Structure types the investigation planner knows. */
export const STRUCTURES = {
  jacket: { label: 'Piled jacket', depth: 'to below the planned pile tip by several pile diameters — the ground under the tip carries the end bearing', where: 'one location per leg group for a single platform; at every foundation for a wind farm' },
  monopile: { label: 'Monopile', depth: 'to below the planned pile tip with a margin of the order of the pile diameter: the whole embedded length resists the overturning moment', where: 'at every foundation location' },
  gravity: { label: 'Gravity base or mudmat', depth: 'to a depth of about the width of the base or more: the ground is loaded to that depth, and soft layers within it govern settlement and sliding', where: 'several points across the footprint, to show how the ground varies under the base' },
  suction: { label: 'Suction anchor or suction bucket', depth: 'to below the planned skirt tip by about one diameter', where: 'at every anchor or bucket: one hard layer or boulder bed stops installation' },
  drag: { label: 'Drag or plate anchor', depth: 'over the depth the anchor is expected to reach — a few metres in sand and stiff clay, tens of metres in soft clay', where: 'along the anchor drag path, not only at the touchdown point' },
  route: { label: 'Pipeline or cable route', depth: 'the top few metres: the burial depth plus a metre or two', where: 'at regular spacing along the route — of the order of one test per kilometre, closer where the geophysical survey shows a change of ground — with deeper tests at crossings and landfalls' },
};

const TESTING = {
  clay: ['Cone penetration tests with pore-pressure measurement (CPTU); in very soft clay a T-bar or ball penetrometer and in-situ vane measure the strength better than the cone.', 'Undisturbed piston or push samples from boreholes.', 'Laboratory: water content, plasticity, unit weight; triaxial and direct simple shear strength, consolidation (oedometer), and cyclic and sensitivity tests where waves, earthquakes or installation remoulding matter.'],
  sand: ['Cone penetration tests are the main tool: they give density and friction angle where samples cannot be taken undisturbed.', 'Vibrocores or borehole samples for grading.', 'Laboratory: particle-size distribution, minimum and maximum density, triaxial or direct shear, and interface shear against steel for piles and skirts.'],
  mixed: ['Cone penetration tests with pore pressure at close spacing to map the changes, with boreholes placed to sample each kind of ground the cones show.', 'Laboratory: the clay and the sand programmes, each on its own material.'],
  carbonate: ['Boreholes with continuous sampling: cemented layers matter and the cone alone does not show what they are made of.', 'Cone penetration tests, expecting refusal on cemented layers.', 'Laboratory: carbonate content, grain crushing and compressibility, cementation, and constant-normal-stiffness interface shear for pile friction.'],
  rock: ['Rotary cored boreholes: recovery, fracture spacing (RQD) and weathering.', 'Laboratory: unconfined compressive strength and point-load tests; abrasivity where drilling or trenching is planned.', 'A cone test serves only to find the depth to rock.'],
  coarse: ['Cone tests may refuse on gravel and cobbles: plan boreholes, and rely on the geophysical survey for layer thickness.', 'Laboratory: particle-size distribution; large-diameter samples are needed to be representative.'],
};
const TEST_FAMILY = { mud: 'clay', sandyMud: 'clay', clay: 'clay', silt: 'clay', lithogenous: 'mixed', silMud: 'clay', silOoze: 'clay', sand: 'sand', muddySand: 'sand', gravel: 'coarse', mixed: 'mixed', volcanic: 'mixed', biogenic: 'carbonate', calcOoze: 'carbonate', calcMud: 'carbonate', mixedOoze: 'carbonate', rock: 'rock' };

/**
 * What a site investigation for this seabed and structure normally contains. Planning guidance in general terms,
 * to be checked against the documents named in `basis`.
 * @param {string|null} cls key of SEABED_CLASSES
 * @param {string} structure key of STRUCTURES
 * @returns {{ structure: string, depth: string, where: string, testing: string[], first: string, basis: string }}
 */
export function investigationPlan(cls, structure, { carbonate = null } = {}) {
  const s = STRUCTURES[structure] || STRUCTURES.jacket;
  const fam = carbonate >= 50 ? 'carbonate' : TEST_FAMILY[cls] || 'mixed';
  return { structure: s.label, depth: s.depth, where: s.where, testing: TESTING[fam],
    first: 'A geophysical survey comes first (bathymetry, sub-bottom profiler, side-scan sonar): it shows where the ground changes and where to place the tests.',
    basis: 'ISO 19901-8 (marine soil investigations: planning, drilling and sampling, in-situ and laboratory testing, reporting); DNV-RP-C212 (offshore soil mechanics and geotechnical engineering); DNV-ST-0126 for wind turbine structures; the SUT Offshore Site Investigation and Geotechnics Committee guidance notes on planning geophysical and geotechnical investigations. The scope, depths and spacing for a project are set by the governing document and the certifier.' };
}
