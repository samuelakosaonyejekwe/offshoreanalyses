import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { classifySeabed, classifyGrainSize, gradeSeabedSamples, reduceSeabedSamples, SEABED_CLASSES, SHELF_DEPTH, worldSeabed } from '../src/data/site.js';
import { shomClass, classifyObservation } from '../src/data/seabed.js';
import { SEABED_GRID_CLASSES, DEEPSEA_GRID_CLASSES, decodeSeabedGrid, seabedGridAt, seabedGridCodes, seabedGridSize } from '../src/data/seabedgrid.js';
import { bedStress, largestMobileGrain, hydroSorting, combineSeabedEvidence, seabedConfidence, FINEST_SAND, MUD_DEPOSITION_STRESS } from '../src/data/seabedevidence.js';
import { screeningProfile, investigationPlan, STRUCTURES, SCREENING_KEYS } from '../src/data/soilscreening.js';

test('seabed descriptions: Folk class names', () => {
  for (const [t, c] of [['Sand', 'sand'], ['Mud', 'mud'], ['sandy Mud', 'sandyMud'], ['muddy Sand', 'muddySand'], ['sandy Gravel', 'gravel'], ['Mud to muddy Sand', 'sandyMud'],
    ['Coarse-grained sediment', 'gravel'], ['Mixed sediment', 'mixed'], ['Rock & boulders', 'rock'], ['gravelly muddy Sand', 'muddySand'], ['gravelly Sand', 'sand'],
    // a bracketed or "slightly" adjective does not change the class
    ['(gravelly) Sand', 'sand'], ['fine sand with shells', 'sand'], ['SLIGHTLY GRAVELLY SAND', 'sand'], ['(slightly gravelly) muddy Sand', 'muddySand'], ['(sandy) Mud', 'mud']]) assert.equal(classifySeabed(t), c, t);
});
test('seabed descriptions: chart abbreviations take the first constituent', () => {
  for (const [t, c] of [['SD,BRN,CRS,BRKN SHELLS', 'sand'], ['M', 'mud'], ['Cy,S', 'mud'], ['RK', 'rock'], ['SH,SD', 'biogenic'], ['G,S', 'gravel'], ['BRN,FNE,SD', 'sand'], ['hard', 'rock']]) assert.equal(classifySeabed(t), c, t);
  for (const t of ['', null, 'No data at this level', 'BRN,FNE']) assert.equal(classifySeabed(t), null);
});
test('seabed samples: nearest described samples decide, far ones are ignored', () => {
  const at = (dLat, text, year = 1990) => ({ lat: 30 + dLat, lon: -80, text, year });
  const r = reduceSeabedSamples([at(0.01, 'SD,GY'), at(0.02, 'SD'), at(0.03, 'M'), at(0.5, 'RK'), at(0.04, 'BRN')], 30, -80);
  assert.equal(r.cls, 'sand'); assert.equal(r.n, 3); assert.equal(r.agree, 2);
  assert.ok(Math.abs(r.nearestKm - 1.112) < 0.01);
  assert.equal(reduceSeabedSamples([at(0.5, 'RK')], 30, -80), null);
});
test('every class has a label, a meaning and a measured grain size with its range and source', () => {
  for (const [k, c] of Object.entries(SEABED_CLASSES)) {
    assert.ok(c.label && c.means.length > 40, k); assert.ok(c.d50Source.length > 15, `${k} says where its grain size comes from`);
    if (c.d50 != null) { assert.ok(c.d50 > 1e-6 && c.d50 < 0.1, k); assert.ok(c.d50Range[0] <= c.d50 && c.d50 <= c.d50Range[1], `${k}: typical value inside its range`); }
    if (c.d50Range) assert.ok(c.d50Range[0] > 0 && c.d50Range[0] < c.d50Range[1], k);
  }
  // Wentworth limits: mud below 0.063 mm, sand up to 2 mm; the gravel class holds every sample with 30 % gravel or more
  assert.ok(SEABED_CLASSES.mud.d50 < 6.3e-5 && SEABED_CLASSES.sand.d50 > 6.3e-5 && SEABED_CLASSES.sand.d50 < 2e-3 && SEABED_CLASSES.gravel.d50 > SEABED_CLASSES.sand.d50);
  assert.ok(SEABED_CLASSES.mud.d50 < SEABED_CLASSES.sandyMud.d50 && SEABED_CLASSES.sandyMud.d50 < SEABED_CLASSES.muddySand.d50 && SEABED_CLASSES.muddySand.d50 < SEABED_CLASSES.sand.d50);
  // no measurements, no figure
  for (const k of ['rock', 'mixed', 'biogenic', 'silt', 'mixedOoze', 'volcanic', 'lithogenous']) assert.equal(SEABED_CLASSES[k].d50, null, k);
});
test('seabed samples: a measured mean grain size of the winning class becomes the grain size', () => {
  const at = (dLat, o) => ({ lat: -20 + dLat, lon: 116, year: 2017, ...o });
  // three sand samples 0.2, 0.4 and 0.3 mm (median 0.3 mm); the muddy sample's 0.05 mm does not count
  const r = reduceSeabedSamples(gradeSeabedSamples([at(0.01, { sand: 95, mud: 5, meanSize: 2e-4 }), at(0.02, { sand: 96, mud: 4, meanSize: 4e-4 }), at(0.03, { sand: 97, mud: 3, meanSize: 3e-4 }), at(0.04, { sand: 20, mud: 80, meanSize: 5e-5 }), at(0.05, { sand: 99, mud: 1 })]), -20, 116);
  assert.equal(r.cls, 'sand'); assert.equal(r.nSized, 3); assert.ok(Math.abs(r.d50 - 3e-4) < 1e-12); assert.deepEqual(r.d50Range, [2e-4, 4e-4]);
  // an even number: the geometric mean of the middle two, sqrt(0.2 × 0.4) = 0.2828 mm
  const e = reduceSeabedSamples(gradeSeabedSamples([at(0.01, { sand: 95, mud: 5, meanSize: 2e-4 }), at(0.02, { sand: 96, mud: 4, meanSize: 4e-4 })]), -20, 116);
  assert.ok(Math.abs(e.d50 - Math.sqrt(8e-8)) < 1e-12);
  assert.equal(reduceSeabedSamples(gradeSeabedSamples([at(0.01, { sand: 95, mud: 5 })]), -20, 116).d50, null);
});
test('own observation: a measured grain size outranks percentages, which outrank the description', () => {
  assert.deepEqual(classifyObservation({ text: 'grey fine sand with shells' }), { cls: 'sand', d50: SEABED_CLASSES.sand.d50, measured: false, how: 'from the description' });
  const r = classifyObservation({ text: 'sand', gravel: 0, sand: 30, mud: 70, d50mm: 0.04 });
  assert.equal(r.cls, 'sandyMud'); assert.ok(Math.abs(r.d50 - 4e-5) < 1e-15); assert.equal(r.measured, true);
  assert.equal(classifyObservation({ d50mm: 0.3 }).cls, 'sand'); assert.equal(classifyObservation({ d50mm: 5 }).cls, 'gravel'); assert.equal(classifyObservation({ d50mm: 0.01 }).cls, 'mud');
  assert.equal(classifyObservation({ text: 'rock' }).d50, null);
  for (const o of [{}, { text: 'brown' }, { d50mm: -1 }, { gravel: 0, sand: 0, mud: 0 }]) assert.equal(classifyObservation(o), null);
});

// ---- laboratory grain-size samples (Australia, Canada) ----
test('grain-size analyses: Folk class from the measured fractions, or from a Folk code', () => {
  for (const [g, s, m, c] of [[0, 95, 5, 'sand'], [4, 78, 18, 'muddySand'], [0, 30, 70, 'sandyMud'], [0, 3, 97, 'mud'], [45, 50, 5, 'gravel'], [10, 85, 5, 'sand'], [100, 0, 0, 'gravel'], [0, 0.65, 0.35, 'muddySand']])
    assert.equal(classifyGrainSize({ gravel: g, sand: s, mud: m }), c, `${g}/${s}/${m}`);
  for (const [f, c] of [['S', 'sand'], ['gS', 'sand'], ['(g)S', 'sand'], ['mS', 'muddySand'], ['(g)mS', 'muddySand'], ['gmS', 'muddySand'], ['sM', 'sandyMud'], ['(g)sM', 'sandyMud'], ['gM', 'sandyMud'], ['M', 'mud'], ['(g)M', 'mud'], ['sG', 'gravel'], ['msG', 'gravel'], ['G', 'gravel'], ['Z', 'mud']])
    assert.equal(classifyGrainSize({ folk: f }), c, f);
  for (const x of [{}, { folk: 'unknown' }, { gravel: null, sand: null, mud: null }, { gravel: 0, sand: 0, mud: 0 }]) assert.equal(classifyGrainSize(x), null);
  assert.equal(classifyGrainSize({ gravel: 0, sand: 90, mud: 10, folk: 'M' }), 'sand', 'measured fractions outrank the code');
});
test('grain-size analyses reduce like described samples and keep a readable description', () => {
  const at = (dLat, o) => ({ lat: -20 + dLat, lon: 116, year: 2017, ...o });
  const graded = gradeSeabedSamples([at(0.01, { gravel: 4, sand: 78, mud: 18, carbonate: 88 }), at(0.02, { folk: '(g)mS' }), at(0.03, { gravel: 0, sand: 97, mud: 3 }), at(0.6, { folk: 'G' })]);
  assert.deepEqual(graded.map((x) => x.cls), ['muddySand', 'muddySand', 'sand', 'gravel']);
  assert.equal(graded[0].text, 'gravel 4 %, sand 78 %, mud 18 %'); assert.equal(graded[1].text, 'Folk class (g)mS');
  const r = reduceSeabedSamples(graded, -20, 116);
  assert.equal(r.cls, 'muddySand'); assert.equal(r.n, 3); assert.equal(r.agree, 2); assert.equal(r.samples[0].carbonate, 88);
});

// ---- legend of the Shom worldwide seabed sediment map ----
test('Shom legend codes decode to a class, the legend name and the size band of the main sediment', () => {
  for (const [code, cls, name] of [['NFRoche', 'rock', 'Rock'], ['NFC', 'gravel', 'Pebble'], ['NFG', 'gravel', 'Gravel'], ['NFGS', 'gravel', 'Sandy gravel'], ['NFS', 'sand', 'Sand'], ['NFSG', 'sand', 'Gravelly sand'],
    ['NFSV', 'muddySand', 'Muddy sand'], ['NFSSi', 'muddySand', 'Silty sand'], ['NFSF', 'sand', 'Fine sand'], ['NFSFV', 'muddySand', 'Muddy fine sand'], ['NFSFSi', 'muddySand', 'Silty fine sand'], ['NFV', 'mud', 'Mud'],
    ['NFVS', 'sandyMud', 'Sandy mud'], ['NFVSF', 'sandyMud', 'Mud with fine sand'], ['NFVG', 'sandyMud', 'Gravelly mud'], ['NFSi', 'silt', 'Silt'], ['NFSiA', 'silt', 'Clayey silt'], ['NFA', 'clay', 'Clay'], ['NFASi', 'clay', 'Silty clay'], ['NFCo', 'biogenic', 'Coral']]) {
    const r = shomClass(code); assert.equal(r?.cls, cls, code); assert.equal(r.name, name); assert.ok(SEABED_CLASSES[r.cls]);
  }
  assert.deepEqual(shomClass('NFSF').band, [5e-5, 5e-4]); assert.deepEqual(shomClass('NFS').band, [5e-4, 2e-3]); assert.deepEqual(shomClass('NFA').band, [1e-6, 1e-5]); assert.equal(shomClass('NFV').band, null);
  for (const c of ['', null, 'NFX', 'SF', 'NFSFX', 'NFNoP']) assert.equal(shomClass(c), null, String(c));
});

// ---- bundled 0.1° grids ----
const SED = new Uint8Array(await readFile(new URL('../assets/seabed-sediment.bin', import.meta.url)));
const DEEP = new Uint8Array(await readFile(new URL('../assets/seabed-deepsea.bin', import.meta.url)));
const deepAt = (lat, lon) => seabedGridAt(DEEP, lat, lon, { reach: 1, classes: DEEPSEA_GRID_CLASSES });
// run-length writer mirroring tools/seabed-grid.py, for the round trip
function encodeGrid(nx, ny, cells) {
  const out = [0x53, 0x42, 0x4c, 0x31, nx & 255, nx >> 8, ny & 255, ny >> 8];
  for (let i = 0; i < cells.length;) {
    let n = 1; while (i + n < cells.length && cells[i + n] === cells[i]) n++;
    if (n <= 15) out.push(cells[i] | (n - 1) << 4);
    else { out.push(cells[i] | 0xf0); let v = n - 16; do { out.push((v & 127) | (v > 127 ? 128 : 0)); v = Math.floor(v / 128); } while (v); }
    i += n;
  }
  return new Uint8Array(out);
}
test('seabed grids: encode → decode round trip, including runs of 15, 16 and many thousand nodes', () => {
  const nx = 400, ny = 201, cells = new Uint8Array(nx * ny);
  let seed = 7; const rnd = () => (seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648;
  for (let i = 0, c = 0; i < cells.length;) { const n = [1, 2, 15, 16, 17, 143, 144, 2063, 20000][Math.floor(rnd() * 9)]; c = (c + 1 + Math.floor(rnd() * 12)) % 14; cells.fill(c, i, i + n); i += n; }
  const buf = encodeGrid(nx, ny, cells), back = decodeSeabedGrid(buf);
  assert.equal(back.nx, nx); assert.equal(back.ny, ny); assert.deepEqual(back.cells, cells);
  const idx = [0, 1, 14, 15, 16, 399, 400, 40199, 80399]; assert.deepEqual(seabedGridCodes(buf, idx), idx.map((i) => cells[i]));
  assert.deepEqual(seabedGridCodes(buf.buffer, [nx * ny - 1, nx * ny]), [cells[nx * ny - 1], -1], 'an ArrayBuffer works too; past the end is -1');
  assert.throws(() => decodeSeabedGrid(buf.subarray(0, buf.length - 3)), /truncated/); assert.throws(() => seabedGridSize(new Uint8Array(16)), /Not a seabed/);
});
test('seabed grids: the bundled files are complete 0.1° grids with the node counts of the packing run', () => {
  const count = (buf) => { const g = decodeSeabedGrid(buf); assert.equal(g.nx, 3600); assert.equal(g.ny, 1801); const c = new Array(14).fill(0); for (const x of g.cells) c[x]++; return c; };
  // nodes per class printed by tools/seabed-grid.py for the two sources
  assert.deepEqual(count(SED), [2858372, 30874, 88845, 39259, 1239, 7074, 166569, 78998, 114655, 9820, 17622, 697961, 926535, 1445777]);
  assert.deepEqual(count(DEEP), [2742573, 1115953, 1030978, 543515, 499567, 551014, 0, 0, 0, 0, 0, 0, 0, 0]);
});
test('sediment grid: nodes equal the answer of the Shom map service at the same points', () => {
  // [lat, lon, legend code returned by the map service for the polygon at the point, group it belongs to]
  const known = [[4.2, 6.5, 'NFS', 'NFS'] /* Gulf of Guinea */, [54, 3, 'NFS', 'NFS'] /* southern North Sea */, [58.5, 4.5, 'NFV', 'NFV'] /* Norwegian Trench */, [30, -160, 'NFA', 'NFA'] /* North Pacific clay */,
    [15, 88, 'NFSiA', 'NFSi'] /* Bay of Bengal */, [-23, -41, 'NFASi', 'NFA'] /* off Brazil */, [26.5, 52, 'NFSF', 'NFSF'] /* Persian Gulf */, [1.2, 103.9, 'NFS', 'NFS'] /* Singapore Strait */, [-20, 116, 'NFSF', 'NFSF'], [44, -60, 'NFS', 'NFS'],
    [54.05, 7.95, 'NFSFV', 'NFSFV'], [43.3, 4.7, 'NFV', 'NFV'], [51.4, -6, 'NFSV', 'NFSV'], [46.2, -124.2, 'NFS', 'NFS'], [53.9, -5.6, 'NFVS', 'NFVS']];
  for (const [lat, lon, code, group] of known) {
    const r = seabedGridAt(SED, lat, lon);
    assert.equal(r?.code, group, `${lat}, ${lon}`); assert.equal(r.shifted, false); assert.ok(r.offsetKm < 8); assert.equal(r.stepDeg, 0.1);
    assert.equal(shomClass(r.code).cls, shomClass(code).cls, 'the group keeps the class of the legend code');
  }
  assert.equal(seabedGridAt(SED, 54.04, 3.04).code, 'NFS'); assert.equal(seabedGridAt(SED, 30, 200).code, 'NFA', 'longitudes wrap');
});
test('sediment grid: land gives null, a coastal gap takes the closest node with data', () => {
  for (const [lat, lon] of [[23, 10] /* Sahara */, [-80, -60] /* beyond the map */, [85, 0], [90, 0], [NaN, 0], [0, Infinity], [95, 0]]) assert.equal(seabedGridAt(SED, lat, lon), null, `${lat}, ${lon}`);
  // Lagos shoreline: the node itself is on land, its neighbour 0.1° away is sand
  const r = seabedGridAt(SED, 6.5, 3.4);
  assert.equal(r.code, 'NFS'); assert.equal(r.shifted, true); assert.ok(r.offsetKm > 10.5 && r.offsetKm < 11.5);
  assert.equal(seabedGridAt(SED, 6.5, 3.4, { reach: 0 }), null);
});
test('deep-sea grid: sediment types of the source, nothing on shelves', () => {
  assert.equal(deepAt(30, -160).name, 'Clay'); assert.equal(deepAt(30, -160).cls, 'clay'); assert.equal(deepAt(-20, -12).cls, 'calcOoze'); assert.equal(deepAt(50, -30).name, 'Calcareous sediment');
  assert.equal(deepAt(5, -140).name, 'Radiolarian ooze'); assert.equal(deepAt(-58, 0).cls, 'lithogenous');
  for (const [lat, lon] of [[4.2, 6.5], [54, 3], [26.5, 52], [23, 10]]) assert.equal(deepAt(lat, lon), null, `${lat}, ${lon}`);
});
test('every grid class maps to a seabed class with an engineering meaning', () => {
  assert.deepEqual(Object.keys(SEABED_GRID_CLASSES).map(Number), [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]);
  for (const [code, c] of Object.entries(SEABED_GRID_CLASSES)) assert.ok(SEABED_CLASSES[shomClass(c.code)?.cls], `${code} → ${c.code}`);
  for (const c of Object.values(DEEPSEA_GRID_CLASSES)) assert.ok(c.name && SEABED_CLASSES[c.cls], c.name);
  for (const k of ['calcOoze', 'mixedOoze', 'calcMud', 'biogenic']) assert.match(SEABED_CLASSES[k].means, /crush/, `${k} warns about crushable carbonate grains`);
  assert.ok(SHELF_DEPTH >= 100 && SHELF_DEPTH <= 300);
});

// ---- hydrodynamic sorting and the combination of evidence ----
test('bed shear stress under waves and current: hand calculation', () => {
  // 30 m of water, Hs 2 m, Tp 7 s, depth-mean current 0.3 m/s, 0.2 mm sand: z0 = d/12 = 1.6667e-5 m.
  // Linear theory: kh = 2.4983, Uw = π·2/(7·sinh kh) = 0.14862 m/s, A = Uw·T/2π = 0.16558 m.
  // fw = 1.39 (A/z0)^-0.52 = 1.39 × 9934.5^-0.52 = 0.011601; τw = ½·1025·fw·Uw² = 0.13133 Pa.
  // Cd = [0.4/(ln(30/z0) − 1)]² = (0.4/13.403)² = 8.9066e-4; τc = 1025·Cd·0.3² = 0.08216 Pa.
  // τm = τc [1 + 1.2 (τw/(τc+τw))^3.2] = 0.08216 × 1.2535 = 0.10299 Pa; in line: τmax = τm + τw = 0.23432 Pa.
  const s = bedStress({ depth: 30, Hs: 2, Tp: 7, U: 0.3, d: 2e-4 });
  assert.ok(Math.abs(s.Uw - 0.14862) < 1e-4); assert.ok(Math.abs(s.tauW - 0.13133) < 2e-4); assert.ok(Math.abs(s.tauC - 0.08216) < 1e-4); assert.ok(Math.abs(s.tauMax - 0.23432) < 3e-4);
  assert.equal(bedStress({ depth: 30, d: 2e-4 }).tauMax, 0);
  // a very rough bed under a small orbit: the wave friction factor is capped at 0.3, τw = ½·1025·0.3·Uw²
  const r = bedStress({ depth: 10, Hs: 0.3, Tp: 3, d: 0.05 }); assert.ok(r.Uw > 0 && r.Uw < 0.02); assert.ok(Math.abs(r.tauW - 0.5 * 1025 * 0.3 * r.Uw ** 2) < 1e-12);
});
test('largest mobile grain: moves at the size returned, not just above it', () => {
  const c = { depth: 30, Hs: 2, Tp: 7, U: 0.3 }, d = largestMobileGrain(c);
  assert.ok(d > 5e-4 && d < 2e-3, `about 0.9 mm, got ${d}`);
  // threshold of 0.2 mm sand in sea water is 0.176 Pa (Soulsby & Whitehouse: D* = 4.07, θcr = 0.0553), below the 0.234 Pa above
  assert.equal(largestMobileGrain({ depth: 150, Hs: 1, Tp: 5, U: 0.02 }), null);
  assert.equal(largestMobileGrain({ depth: 3, Hs: 2.2, Tp: 9, U: 1.5 }), 0.064, 'capped at the coarsest size examined');
  assert.ok(FINEST_SAND === 6.3e-5 && MUD_DEPOSITION_STRESS === 0.05);
});
test('hydrodynamic setting: worked, calm or between', () => {
  const worked = hydroSorting({ depth: 30, waves: { p90: 2, hs1: 5, tp1: 10 }, current: { p90: 0.4, max: 0.9 } });
  assert.equal(worked.setting, 'worked'); assert.ok(worked.frequent.tau > 0.1763); assert.ok(worked.dStorm > worked.dFrequent); assert.match(worked.text, /worked all the time/);
  assert.ok(Math.abs(worked.frequent.Tp - 4.4 * Math.sqrt(2)) < 1e-12, 'Tp = 4.4 √Hs when no factor is given'); assert.equal(worked.storm.Tp, 10); assert.ok(Math.abs(worked.frequent.U - 0.35) < 1e-12, 'depth mean = 7/8 of the surface speed');
  const calm = hydroSorting({ depth: 150, waves: { p90: 1.5, hs1: 4 }, current: { p90: 0.05, max: 0.2 } });
  assert.equal(calm.setting, 'calm'); assert.ok(calm.frequent.tau < 0.05); assert.match(calm.text, /mud settles and stays/);
  const mid = hydroSorting({ depth: 60, waves: { p90: 2, hs1: 6 }, current: { p90: 0.35, max: 0.6 } });
  assert.equal(mid.setting, 'between'); assert.ok(mid.frequent.tau >= 0.05 && mid.frequent.tau < 0.1763);
  assert.equal(hydroSorting({ depth: 30 }), null); assert.equal(hydroSorting({ depth: -5, waves: { p90: 2, hs1: 5 } }), null);
});
test('evidence: agreement passes a grain size on, a contradiction or an undecided check passes none', () => {
  const worked = { setting: 'worked' }, calm = { setting: 'calm' }, between = { setting: 'between' };
  const sand = { cls: 'sand', d50: SEABED_CLASSES.sand.d50, d50Range: SEABED_CLASSES.sand.d50Range }, mud = { cls: 'mud', d50: SEABED_CLASSES.mud.d50, d50Range: SEABED_CLASSES.mud.d50Range };
  const a = combineSeabedEvidence({ ...sand, hydro: worked }); assert.equal(a.verdict, 'agree'); assert.equal(a.d50, 2.6e-4); assert.deepEqual(a.range, [1e-4, 1e-3]);
  const b = combineSeabedEvidence({ ...mud, hydro: calm }); assert.equal(b.verdict, 'agree'); assert.equal(b.d50, 7.8e-6);
  const c = combineSeabedEvidence({ ...mud, hydro: worked }); assert.equal(c.verdict, 'conflict'); assert.equal(c.d50, null); assert.ok(c.check.some((x) => /river/.test(x)));
  const d = combineSeabedEvidence({ ...sand, hydro: calm }); assert.equal(d.verdict, 'conflict'); assert.equal(d.d50, null); assert.ok(d.check.some((x) => /[Rr]elict/.test(x)));
  for (const h of [between, null]) { const o = combineSeabedEvidence({ ...sand, hydro: h }); assert.equal(o.verdict, 'open'); assert.equal(o.d50, null); }
  for (const cls of ['rock', 'mixed', 'volcanic']) assert.equal(combineSeabedEvidence({ cls, hydro: worked }).verdict, 'open', cls);
  // a muddy sand in a worked setting: never finer than very fine sand
  assert.equal(combineSeabedEvidence({ cls: 'muddySand', d50: 5e-5, d50Range: [4e-5, 4e-4], hydro: worked }).d50, 6.3e-5);
});
test('world map tier: on a shelf a grain size needs the hydrodynamics to agree; in deep water the sediment type names the class', () => {
  const worked = hydroSorting({ depth: 30, waves: { p90: 2, hs1: 5, tp1: 10 }, current: { p90: 0.4, max: 0.9 } });
  const a = worldSeabed({ depth: 30, live: { code: 'NFSF', sheet: 8 }, grid: null, hydro: worked });
  assert.equal(a.cls, 'sand'); assert.equal(a.basis, 'worldmap'); assert.equal(a.described, 'Fine sand'); assert.equal(a.sheet, 8); assert.equal(a.shelf, true); assert.equal(a.evidence.verdict, 'agree'); assert.equal(a.d50, 2.6e-4); assert.equal(seabedConfidence(a).key, 'consistent');
  const b = worldSeabed({ depth: 30, live: { code: 'NFV', sheet: 1 }, hydro: worked });
  assert.equal(b.cls, 'mud'); assert.equal(b.evidence.verdict, 'conflict'); assert.equal(b.d50, null); assert.equal(b.d50Class, 7.8e-6); assert.equal(seabedConfidence(b).key, 'conflict');
  const c = worldSeabed({ depth: 30, live: null, grid: { sediment: seabedGridAt(SED, 4.2, 6.5), deepsea: null }, hydro: null });
  assert.equal(c.cls, 'sand'); assert.equal(c.basis, 'grid'); assert.equal(c.d50, null); assert.ok(Math.abs(c.resolutionKm - 11.12) < 0.01); assert.equal(seabedConfidence(c).key, 'single');
  // 5,000 m of water in the North Pacific: clay on both; off Brazil the map says silty clay, the deep-sea grid calcareous sediment
  const d = worldSeabed({ depth: 5000, live: { code: 'NFA', sheet: 1 }, grid: { sediment: null, deepsea: deepAt(30, -160) } });
  assert.equal(d.cls, 'clay'); assert.equal(d.shelf, false); assert.equal(d.evidence, undefined); assert.equal(d.d50, 4.5e-6); assert.equal(d.composition, 'Clay');
  const e = worldSeabed({ depth: 3000, live: { code: 'NFASi', sheet: 1 }, grid: { sediment: null, deepsea: deepAt(-20, -12) } });
  assert.equal(e.cls, 'calcOoze'); assert.equal(e.described, 'Silty clay'); assert.equal(e.composition, 'Calcareous sediment'); assert.equal(e.d50, 6.9e-6);
  // land-derived sediment is not a class of its own where the map gives a grain-size class; a deep sand stays sand
  assert.equal(worldSeabed({ depth: 3000, live: { code: 'NFSi' }, grid: { sediment: null, deepsea: deepAt(-58, 0) } }).cls, 'silt');
  assert.equal(worldSeabed({ depth: 3000, live: { code: 'NFS' }, grid: { sediment: null, deepsea: deepAt(-20, -12) } }).cls, 'sand');
  // the deep-sea grid is not used above 500 m; nothing anywhere gives null
  assert.equal(worldSeabed({ depth: 300, live: null, grid: { sediment: null, deepsea: deepAt(30, -160) } }), null);
  assert.equal(worldSeabed({ depth: 30, live: { code: 'NFX' }, grid: { sediment: null, deepsea: null } }), null);
});
test('confidence grades', () => {
  for (const [sb, key] of [[null, 'none'], [{ cls: null, basis: 'none' }, 'none'], [{ cls: 'sand', basis: 'map' }, 'mapped'], [{ cls: 'sand', basis: 'samples' }, 'measured'], [{ cls: 'sand', basis: 'grainsize' }, 'measured'],
    [{ cls: 'sand', basis: 'worldmap', evidence: { verdict: 'open' } }, 'single'], [{ cls: 'clay', basis: 'worldmap' }, 'single'], [{ cls: 'sand', basis: 'observation' }, 'observed']]) {
    const c = seabedConfidence(sb); assert.equal(c.key, key); assert.ok(c.label && c.raise.length > 30);
  }
});

// ---- screening soil profile and investigation planner ----
test('screening profile: offered for ordinary clay and sand only, with values from the weak end of the published ranges', () => {
  for (const cls of ['mud', 'sandyMud', 'clay']) {
    const p = screeningProfile(cls); assert.equal(p.offered, true, cls);
    // 0.20 × σ′v with γ′ = 6 kN/m³ gives 1.2 kPa/m, the low end of the 1.2–2 kPa/m range; mudline 2 kPa inside 1–5 kPa
    assert.deepEqual(p.values, { 'site.su0': 2000, 'site.suGrad': 1200, 'site.gammaSub': 5000 }); assert.equal(p.values['site.phiDeg'], undefined);
  }
  assert.deepEqual(screeningProfile('sand').values, { 'site.phiDeg': 30, 'site.gammaSub': 9000 }); assert.deepEqual(screeningProfile('gravel').values, screeningProfile('sand').values);
  assert.equal(screeningProfile('muddySand').values['site.phiDeg'], 29);
  for (const cls of [null, 'rock', 'silt', 'mixed', 'volcanic', 'biogenic', 'calcOoze', 'calcMud', 'mixedOoze', 'silOoze', 'silMud', 'lithogenous']) {
    const p = screeningProfile(cls); assert.equal(p.offered, false, String(cls)); assert.equal(p.values, undefined); assert.ok(p.why.length > 30 && p.refs.length);
  }
  // a sand that the samples show to be carbonate gets no profile
  assert.equal(screeningProfile('sand', { carbonate: 85 }).offered, false); assert.equal(screeningProfile('sand', { carbonate: 20 }).offered, true);
  for (const k of Object.keys({ ...screeningProfile('mud').values, ...screeningProfile('sand').values })) assert.ok(SCREENING_KEYS.includes(k), k);
  assert.ok(!SCREENING_KEYS.includes('site.d50'));
});
test('investigation planner: every structure and every class gives a programme', () => {
  for (const s of Object.keys(STRUCTURES)) for (const cls of [null, ...Object.keys(SEABED_CLASSES)]) {
    const p = investigationPlan(cls, s); assert.ok(p.structure && p.depth && p.where && p.testing.length >= 2 && /ISO 19901-8/.test(p.basis), `${s} / ${cls}`);
  }
  assert.match(investigationPlan('rock', 'jacket').testing[0], /cored/); assert.match(investigationPlan('sand', 'route', { carbonate: 90 }).testing.join(' '), /carbonate content/);
});
