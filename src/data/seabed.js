// Seabed-surface sediment: the classes the app works with, their measured grain sizes, and the functions that
// turn a description, a laboratory analysis or a map legend code into one of those classes. Pure functions.

/** Water depth (m) taken as the edge of the continental shelf: shallower, a mapped class needs a second line of evidence before a grain size is passed on. */
export const SHELF_DEPTH = 200;

// Where the grain sizes below come from. Each class stores the typical value (d50), the range between the 10th
// and 90th percentile of the measurements (d50Range) and this note (d50Source); nothing here is a chosen figure.
const LAB = (n) => `Measured mean grain size of ${n.toLocaleString('en')} seabed-surface samples of this class (median, with the 10th to 90th percentile as the range): laboratory analyses of the Geological Survey of Canada (Open Government Licence – Canada) and the Geoscience Australia Marine Sediments database (CC BY 4.0), computed October 2026.`;
const SVALNOV = 'Sval’nov & Alekseeva (2005), Grain size composition of ocean and sea bottom sediments, PANGAEA, doi:10.1594/PANGAEA.746414 (CC BY 3.0; sieve and pipette analyses of deep-sea cores)';
const PIPETTE = 'The size fractions of the same tables place the median one size class finer than the published median column, so the lower end of the range is set at 1 µm.';

/**
 * Seabed-surface classes with a measured median grain size and what each usually means for foundations, scour
 * and burial. Grain size is a property of the class; strength is not — no strength is ever derived from a
 * surface description. `d50` is null where no measurements support a single figure.
 */
export const SEABED_CLASSES = {
  rock: { label: 'Rock or boulders', d50: null, d50Range: null, d50Source: 'Not a granular bed.', means: 'Driven piles may refuse and trenching is slow: expect drilled or rock-socketed foundations and surface-laid, protected cables or pipelines.' },
  gravel: { label: 'Gravel or coarse-grained sediment', d50: 1.8e-3, d50Range: [5.9e-4, 1.2e-2], d50Source: `${LAB(3772)} The class holds every sample with at least 30 % gravel, so its mean lies in very coarse sand.`, means: 'Free-draining and resistant to scour; hard on ploughs and jetting tools, and suction installation may not seal.' },
  mixed: { label: 'Mixed sediment', d50: null, d50Range: null, d50Source: 'No single grain size describes a mixed bed.', means: 'Variable over short distances: expect both granular and cohesive behaviour and plan the site investigation to map the boundaries.' },
  sand: { label: 'Sand', d50: 2.6e-4, d50Range: [1.0e-4, 1.0e-3], d50Source: LAB(7275), means: 'Drained behaviour and good bearing, but mobile under waves and current: scour protection and burial depth need attention.' },
  muddySand: { label: 'Muddy sand', d50: 1.0e-4, d50Range: [4.7e-5, 4.0e-4], d50Source: LAB(2372), means: 'Mostly granular with fines: partly drained response, moderate scour potential.' },
  sandyMud: { label: 'Sandy mud', d50: 2.0e-5, d50Range: [8.0e-6, 5.6e-5], d50Source: LAB(2228), means: 'Mostly cohesive: undrained strength governs and consolidation settlement should be expected.' },
  mud: { label: 'Mud (silt and clay)', d50: 7.8e-6, d50Range: [2.7e-6, 1.5e-5], d50Source: LAB(3074), means: 'Soft cohesive seabed: low bearing capacity near the surface, large embedment of pipelines and mudmats, favourable for suction and drag anchors.' },
  biogenic: { label: 'Shell, coral or other biogenic material', d50: null, d50Range: null, d50Source: 'No grain-size measurements were found for this class: shell and coral debris ranges from sand to cobble size.', means: 'Carbonate material crushes under load: pile shaft friction can be far below that of silica sand, and coral is usually an environmental constraint.' },
  // classes below occur on continental margins and in the deep sea
  silt: { label: 'Silt', d50: null, d50Range: [1e-5, 5e-5], d50Source: 'No measured median was found for this class. The range is the definition of “silt” in the legend of the Shom worldwide seabed sediment map: at least half of the grains between 0.01 and 0.05 mm.', means: 'Fine, low-plasticity sediment between drained and undrained behaviour: easily eroded and jetted, sensitive to cyclic loading, and slow to regain strength after disturbance around anchors and trenches.' },
  clay: { label: 'Clay (terrigenous or deep-sea pelagic clay)', d50: 4.5e-6, d50Range: [1e-6, 6.9e-6], d50Source: `Median grain size of 119 samples from the top 0.3 m of pelagic clay (median 4.5 µm, 10th–90th percentile 3.7–5.1 µm; 104 further samples of eupelagic clay reach 6.9 µm): ${SVALNOV}, tables 4 and 5. ${PIPETTE}`, means: 'Very soft, highly plastic clay at the mudline with strength rising slowly with depth: suited to suction, drag-embedment and plate anchors; expect deep pipeline and mudmat embedment, long consolidation times and low axial pipe-soil resistance.' },
  calcOoze: { label: 'Calcareous ooze', d50: 6.9e-6, d50Range: [1e-6, 3.1e-4], d50Source: `Median grain size of 120 samples from the top 0.3 m of coccolith ooze (median 6.9 µm, 10th–90th percentile 5.8–8.3 µm); foraminiferal ooze is a mud with sand-sized shells and reaches 0.31 mm (90th percentile of 34 samples): ${SVALNOV}, tables 12 and 13. ${PIPETTE}`, means: 'Soft carbonate sediment of hollow microfossil shells: the grains crush under load, so pile shaft friction and anchor capacity can be far below those of a silica soil of the same grading, and it is compressible and often lightly cemented in layers — design with carbonate-specific methods and laboratory tests.' },
  silOoze: { label: 'Siliceous ooze (diatoms, radiolarians or sponge spicules)', d50: 5.1e-6, d50Range: [1e-6, 2.4e-5], d50Source: `Median grain size of 221 samples from the top 0.3 m of radiolarian ooze (median 5.1 µm, 10th–90th percentile 4.1–5.8 µm) and 35 samples of diatom ooze (4.7 µm): ${SVALNOV}, tables 8 and 9. ${PIPETTE} Laser diffraction gives 15–24 µm for the top centimetre of abyssal Pacific sediment of the same province (Macheriotou et al. 2022, PANGAEA, doi:10.1594/PANGAEA.941216 and doi:10.1594/PANGAEA.952416, CC BY 4.0), which sets the upper end: the answer depends on the method.`, means: 'Very porous, light sediment of fragile opal skeletons: unit weight is exceptionally low and water content high, the grains crush, and index-test correlations for ordinary clays do not apply — anchors and mudmats need site-specific strength testing.' },
  mixedOoze: { label: 'Mixed calcareous and siliceous ooze', d50: null, d50Range: null, d50Source: 'No grain-size measurements were found for this class.', means: 'Soft biogenic sediment with crushable carbonate and opal grains: low unit weight, high compressibility and low pile or anchor friction; treat as a carbonate soil until tests show otherwise.' },
  volcanic: { label: 'Volcanic ash, sand or gravel', d50: null, d50Range: [8.8e-6, 7.4e-4], d50Source: `Two populations, so no single figure: 17 samples of acid ash (median 33 µm, 10th–90th percentile 8.8–196 µm) and 6 of basaltic volcanic sand (median 0.30 mm, 0.07–0.74 mm): ${SVALNOV}, tables 18 and 19.`, means: 'Angular, vesicular grains that crush and may be altered to clay: friction angles look high in the laboratory but pile and anchor friction is unreliable, the ground is abrasive to trenching tools, and the area is likely to be seismically active.' },
  silMud: { label: 'Siliceous mud', d50: 5.2e-6, d50Range: [1e-6, 6.0e-6], d50Source: `Median grain size of 28 samples from the top 0.3 m of radiolarian clayey mud (median 5.2 µm, 10th–90th percentile 4.5–5.4 µm; 6.0 µm at depth): ${SVALNOV}, table 10. ${PIPETTE}`, means: 'Soft cohesive mud rich in opal microfossils: lower unit weight and higher water content than ordinary clay of the same strength, large embedment of pipelines and mudmats, favourable for suction and drag anchors once the strength profile is measured.' },
  calcMud: { label: 'Fine-grained calcareous sediment', d50: 5.2e-6, d50Range: [1e-6, 6.9e-6], d50Source: `Median grain size of 13 samples from the top 0.3 m of calcareous clayey mud (median 5.2 µm, 10th–90th percentile 4.3–6.9 µm): ${SVALNOV}, table 15. ${PIPETTE}`, means: 'Carbonate mud or marl: soft and compressible near the surface, with crushable grains and patchy cementation — pile friction is low and variable, and hard cemented layers can stop ploughs and suction-anchor skirts.' },
  lithogenous: { label: 'Land-derived (lithogenous) deep-sea sediment', d50: null, d50Range: null, d50Source: 'No measured median was found for this class: surface samples of continental margins below 1,000 m hold about half silt, with sand between 1 and 75 % (Lamy et al. 1998, doi:10.1594/PANGAEA.58251; Frenz et al. 2003, doi:10.1594/PANGAEA.70693).', means: 'Mud, silt and sand washed off the continents — fans, slopes and turbidite plains: layering of soft clay with sand beds is usual, so anchors and piles meet abrupt changes of strength, and slopes may be unstable.' },
};

const SEABED_WORDS = [
  ['rock', /\b(rock|rk|r|rky|boulders?|blds?|bo|bedrock|ledge|lava|hard)\b/], ['biogenic', /\b(coral|co|crl|shells?|sh|oysters?|oys|maerl|sponge|weed|grass|kelp)\b/],
  ['gravel', /\b(gravel|gravelly|g|grv|pebbles?|p|cobbles?|cb|stones?|st|shingle|coarse-grained)\b/], ['mixed', /\bmixed\b/],
  ['sand', /\b(sand|sandy|sd|s)\b/], ['mud', /\b(mud|muddy|md|m|clay|cy|cl|silt|si|ooze|oz|marl|mrl)\b/],
];
/** Class of one seabed description: a Folk class name ("muddy Sand"), chart abbreviations ("SD,BRN,CRS") or plain words. */
export function classifySeabed(text) {
  // a bracketed adjective means "slightly" on the Folk scheme ("(gravelly) Sand"): it does not change the class
  // what follows "with" is a minor constituent ("sand with shells")
  const t = String(text || '').toLowerCase().replace(/\bwith\b[^,;/]*/g, ' ').replace(/\((slightly\s+)?(gravelly|sandy|muddy)\)/g, ' ').replace(/\bslightly (gravelly|sandy|muddy)\b/g, ' ').replace(/[;/]/g, ',').replace(/\s+/g, ' ').trim();
  if (!t || /no data|unknown/.test(t)) return null;
  // Folk names: the noun (last word) is the main fraction, an adjective before it the secondary one
  if (/\bmud to muddy sand\b/.test(t)) return 'sandyMud';
  if (/\b(sandy|gravelly sandy) mud\b/.test(t)) return 'sandyMud';
  if (/\b(muddy|gravelly muddy) sand\b/.test(t)) return 'muddySand';
  if (/\b(sandy|muddy sandy|muddy) gravel\b/.test(t)) return 'gravel';
  if (/\bgravelly sand\b/.test(t)) return 'sand';
  if (/\bgravelly mud\b/.test(t)) return 'sandyMud';
  // otherwise the first recognised term of the description is the main constituent (chart convention)
  for (const part of t.split(',').map((x) => x.trim()).filter(Boolean)) {
    for (const [k, re] of SEABED_WORDS) if (re.test(part)) return k;
  }
  return null;
}
/** Class of a laboratory grain-size analysis: gravel, sand and mud fractions (any common scale), or failing those a Folk code such as "(g)mS". */
export function classifyGrainSize({ gravel, sand, mud, folk } = {}) {
  const g = Math.max(0, gravel || 0), s = Math.max(0, sand || 0), m = Math.max(0, mud || 0);
  if (g + s + m > 0) {
    if (g / (g + s + m) >= 0.3 || s + m === 0) return 'gravel';
    const r = s / (s + m);
    return r >= 0.9 ? 'sand' : r >= 0.5 ? 'muddySand' : r >= 0.1 ? 'sandyMud' : 'mud';
  }
  // Folk code: the capital is the main fraction, lower-case letters the modifiers, brackets "slightly"
  const code = /^((?:\(?[gms]\)?)*)([GSMZC])$/.exec(String(folk || '').trim());
  if (!code) return null;
  const mod = code[1].replace(/\([gms]\)/g, '');
  return code[2] === 'G' ? 'gravel' : code[2] === 'S' ? (mod.includes('m') ? 'muddySand' : 'sand') : /[sg]/.test(mod) ? 'sandyMud' : 'mud';
}
/** Grain-size analyses in the shape reduceSeabedSamples expects: a class and a one-line description each. */
export function gradeSeabedSamples(samples) {
  const pc = (v) => `${Math.round(v || 0)} %`;
  return samples.map((sm) => ({ ...sm, cls: classifyGrainSize(sm), text: sm.sand != null || sm.mud != null || sm.gravel != null ? `gravel ${pc(sm.gravel)}, sand ${pc(sm.sand)}, mud ${pc(sm.mud)}` : `Folk class ${sm.folk}` }));
}
/**
 * Reduce nearby described samples to a class: the commonest class among the nearest few, with distances. Where
 * the samples of that class carry a measured mean grain size (`meanSize`, m), their median is returned as `d50`
 * with the smallest and largest value as `d50Range`.
 */
export function reduceSeabedSamples(samples, lat, lon, { maxKm = 15, take = 7 } = {}) {
  const R = Math.PI / 180;
  const near = samples.map((sm) => ({ ...sm, cls: sm.cls || classifySeabed(sm.text), km: 6371 * Math.hypot((sm.lon - lon) * R * Math.cos(lat * R), (sm.lat - lat) * R) }))
    .filter((sm) => sm.cls && sm.km <= maxKm).sort((a, b) => a.km - b.km).slice(0, take);
  if (!near.length) return null;
  const votes = {};
  for (const sm of near) votes[sm.cls] = (votes[sm.cls] || 0) + 1 / (1 + sm.km);
  const cls = Object.entries(votes).sort((a, b) => b[1] - a[1])[0][0];
  const sizes = near.filter((sm) => sm.cls === cls && sm.meanSize > 0).map((sm) => sm.meanSize).sort((a, b) => a - b);
  const d50 = sizes.length ? (sizes.length % 2 ? sizes[(sizes.length - 1) / 2] : Math.sqrt(sizes[sizes.length / 2 - 1] * sizes[sizes.length / 2])) : null;
  return { cls, n: near.length, nearestKm: near[0].km, farthestKm: near[near.length - 1].km, agree: near.filter((sm) => sm.cls === cls).length, years: [Math.min(...near.map((sm) => sm.year || 9999)), Math.max(...near.map((sm) => sm.year || 0))],
    d50, d50Range: sizes.length ? [sizes[0], sizes[sizes.length - 1]] : null, nSized: sizes.length, samples: near };
}

// ---- legend of the Shom worldwide seabed sediment map --------------------------------------------------------
// A code is "NF" followed by the main sediment and up to two admixtures, coarsest names first in the legend:
// Roche, C (pebbles, over 20 mm), G (gravel, 2–20 mm), S (sand, 0.5–2 mm), SF (fine sand, 0.05–0.5 mm),
// V (mud: 20 % or more finer than 0.05 mm), Si (silt, 0.01–0.05 mm), A (clay, finer than 0.01 mm).
const SHOM_MAIN = {
  Roche: { cls: 'rock', band: null, def: 'rock outcrop' },
  C: { cls: 'gravel', band: [2e-2, 6.4e-2], def: 'at least half of the grains coarser than 20 mm' },
  G: { cls: 'gravel', band: [2e-3, 2e-2], def: 'at least half of the grains between 2 and 20 mm' },
  S: { cls: 'sand', band: [5e-4, 2e-3], def: 'at least half of the grains between 0.5 and 2 mm' },
  SF: { cls: 'sand', band: [5e-5, 5e-4], def: 'at least half of the grains between 0.05 and 0.5 mm' },
  V: { cls: 'mud', band: null, def: '20 % or more of the grains finer than 0.05 mm' },
  Si: { cls: 'silt', band: [1e-5, 5e-5], def: 'at least half of the grains between 0.01 and 0.05 mm' },
  A: { cls: 'clay', band: [1e-6, 1e-5], def: 'at least half of the grains finer than 0.01 mm' },
};
const SHOM_NAMES = {
  NFRoche: 'Rock', NFC: 'Pebble', NFCG: 'Gravelly pebble', NFCS: 'Sandy pebble', NFCV: 'Muddy pebble', NFGC: 'Pebbly gravel', NFG: 'Gravel', NFGS: 'Sandy gravel', NFGSV: 'Gravel and muddy sand',
  NFGVS: 'Gravel and sandy mud', NFGV: 'Muddy gravel', NFSC: 'Pebbly sand', NFSG: 'Gravelly sand', NFSGV: 'Muddy gravelly sand', NFS: 'Sand', NFSSi: 'Silty sand', NFSVG: 'Sand and gravelly mud',
  NFSV: 'Muddy sand', NFSFC: 'Pebbly fine sand', NFSFG: 'Fine sand and gravel', NFSF: 'Fine sand', NFSFSi: 'Silty fine sand', NFSFV: 'Muddy fine sand', NFVC: 'Pebbly mud', NFVG: 'Gravelly mud',
  NFVSG: 'Mud and sandy gravel', NFVS: 'Sandy mud', NFVSF: 'Mud with fine sand', NFV: 'Mud', NFSi: 'Silt', NFSiA: 'Clayey silt', NFASi: 'Silty clay', NFA: 'Clay', NFNoP: 'Polymetallic nodules', NFCo: 'Coral',
};
/**
 * Decode a legend code of the Shom map ("NFSFV") into a seabed class, the English legend name, the grain-size
 * band the legend defines for the main sediment and that definition in words. Null for an unknown code.
 * @returns {null|{ code: string, cls: string, name: string, band: number[]|null, definition: string }}
 */
export function shomClass(code) {
  const c = String(code || '').trim();
  if (c === 'NFCo') return { code: c, cls: 'biogenic', name: SHOM_NAMES[c], band: null, definition: 'coral' };
  const parts = /^NF(Roche|SF|Si|[CGSVA])((?:SF|Si|[CGSVA])*)$/.exec(c);
  if (!parts || !SHOM_NAMES[c]) return null;
  const main = SHOM_MAIN[parts[1]], rest = parts[2];
  // a sand with mud or silt is a muddy sand; a mud with sand, gravel or pebbles is a sandy mud (Folk groups)
  const cls = main.cls === 'sand' && /V|Si/.test(rest) ? 'muddySand' : main.cls === 'mud' && rest ? 'sandyMud' : main.cls;
  return { code: c, cls, name: SHOM_NAMES[c], band: main.band, definition: main.def };
}

/**
 * A seabed observation entered by the user — a description, measured fractions or a measured median grain size —
 * as a class with the grain size to pass on. A measured d50 always wins; fractions outrank a description.
 * @param {{ text?: string, gravel?: number, sand?: number, mud?: number, d50mm?: number }} o
 * @returns {null|{ cls: string, d50: number|null, measured: boolean, how: string }}
 */
export function classifyObservation({ text, gravel, sand, mud, d50mm } = {}) {
  const d = Number(d50mm) > 0 ? Number(d50mm) / 1000 : null;
  const fr = [gravel, sand, mud].map((v) => (Number(v) > 0 ? Number(v) : 0));
  let cls = fr[0] + fr[1] + fr[2] > 0 ? classifyGrainSize({ gravel: fr[0], sand: fr[1], mud: fr[2] }) : null, how = cls ? 'from the gravel, sand and mud percentages' : '';
  if (!cls && text) { cls = classifySeabed(text); how = cls ? 'from the description' : ''; }
  if (!cls && d) { cls = d >= 2e-3 ? 'gravel' : d >= 6.3e-5 ? 'sand' : 'mud'; how = 'from the measured median grain size (Wentworth scale: gravel from 2 mm, sand from 0.063 mm)'; }
  if (!cls) return null;
  return { cls, d50: d ?? SEABED_CLASSES[cls].d50, measured: d != null, how };
}
