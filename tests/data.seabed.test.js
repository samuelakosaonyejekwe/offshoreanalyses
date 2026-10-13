import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { classifySeabed, classifyGrainSize, gradeSeabedSamples, reduceSeabedSamples, SEABED_CLASSES, SHELF_DEPTH } from '../src/data/site.js';
import { SEABED_GRID_CLASSES, decodeSeabedGrid, seabedGridAt, seabedGridCodes, seabedGridSize } from '../src/data/seabedgrid.js';

test('seabed descriptions: Folk class names', () => {
  for (const [t, c] of [['Sand', 'sand'], ['Mud', 'mud'], ['sandy Mud', 'sandyMud'], ['muddy Sand', 'muddySand'], ['sandy Gravel', 'gravel'], ['Mud to muddy Sand', 'sandyMud'],
    ['Coarse-grained sediment', 'gravel'], ['Mixed sediment', 'mixed'], ['Rock & boulders', 'rock'], ['gravelly muddy Sand', 'muddySand'], ['gravelly Sand', 'sand']]) assert.equal(classifySeabed(t), c, t);
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
test('every class has a label, a meaning and a grain size in the range of its name', () => {
  for (const [k, c] of Object.entries(SEABED_CLASSES)) { assert.ok(c.label && c.means.length > 40, k); if (c.d50 != null) assert.ok(c.d50 > 1e-6 && c.d50 < 0.1, k); }
  assert.ok(SEABED_CLASSES.mud.d50 < 6.3e-5 && SEABED_CLASSES.sand.d50 > 6.3e-5 && SEABED_CLASSES.sand.d50 < 2e-3 && SEABED_CLASSES.gravel.d50 > 2e-3);
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

// ---- bundled global lithology grid ----
const GRID = new Uint8Array(await readFile(new URL('../assets/seabed-lithology.bin', import.meta.url)));
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
test('lithology grid: encode → decode round trip, including runs of 15, 16 and many thousand nodes', () => {
  const nx = 400, ny = 201, cells = new Uint8Array(nx * ny);
  let seed = 7; const rnd = () => (seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648;
  for (let i = 0, c = 0; i < cells.length;) { const n = [1, 2, 15, 16, 17, 143, 144, 2063, 20000][Math.floor(rnd() * 9)]; c = (c + 1 + Math.floor(rnd() * 12)) % 14; cells.fill(c, i, i + n); i += n; }
  const buf = encodeGrid(nx, ny, cells), back = decodeSeabedGrid(buf);
  assert.equal(back.nx, nx); assert.equal(back.ny, ny); assert.deepEqual(back.cells, cells);
  const idx = [0, 1, 14, 15, 16, 399, 400, 40199, 80399]; assert.deepEqual(seabedGridCodes(buf, idx), idx.map((i) => cells[i]));
  assert.deepEqual(seabedGridCodes(buf.buffer, [nx * ny - 1, nx * ny]), [cells[nx * ny - 1], -1], 'an ArrayBuffer works too; past the end is -1');
  assert.throws(() => decodeSeabedGrid(buf.subarray(0, buf.length - 3)), /truncated/); assert.throws(() => seabedGridSize(new Uint8Array(16)), /Not a seabed/);
});
test('lithology grid: the bundled file is the complete 0.1° source grid', () => {
  const g = decodeSeabedGrid(GRID);
  assert.equal(g.nx, 3600); assert.equal(g.ny, 1801);
  const count = new Array(14).fill(0); for (const c of g.cells) count[c]++;
  // node counts per class in seabed_lithology_v1.nc without its repeated +180° column
  assert.deepEqual(count, [2706651, 35544, 134146, 65012, 1521435, 899980, 59384, 260535, 6287, 61090, 9164, 16361, 282407, 425604]);
});
test('lithology grid: known nodes read from the source file', () => {
  // [lat, lon, code in seabed_lithology_v1.nc]
  const known = [[4.2, 6.5, 4] /* Gulf of Guinea */, [15, 88, 4] /* Bay of Bengal */, [-23, -41, 5] /* offshore Brazil */, [-58, 0, 7], [-60, -120, 7] /* Southern Ocean diatom belt */,
    [30, -160, 4] /* North Pacific red clay */, [-20, -12, 5], [30, -42, 5] /* Mid-Atlantic Ridge */, [0, -110, 5], [-10, 60, 5], [55, 170, 12] /* Bering Sea */, [52, -175, 2] /* Aleutian arc */,
    [20, -90.5, 13], [13, 145, 13], [47.3, -22.7, 13], [0, 180, 4], [0, -180, 4], [-52.8, -144.3, 1], [-44.7, 160.2, 3], [-50.4, 123.3, 6], [-56.5, 163.7, 8], [-10.1, -96.5, 9], [-21.5, 153, 10], [-16.2, 168.8, 11]];
  for (const [lat, lon, code] of known) {
    const r = seabedGridAt(GRID, lat, lon);
    assert.equal(r?.code, code, `${lat}, ${lon}`); assert.equal(r.shifted, false); assert.ok(r.offsetKm < 1e-6); assert.equal(r.stepDeg, 0.1);
    assert.equal(r.name, SEABED_GRID_CLASSES[code].name); assert.equal(r.cls, SEABED_GRID_CLASSES[code].cls);
  }
  // off the nodes the nearest one answers, and longitudes wrap
  assert.equal(seabedGridAt(GRID, -23.04, -40.96).code, 5); assert.equal(seabedGridAt(GRID, -58, 360).code, 7); assert.equal(seabedGridAt(GRID, 0, 179.97).code, 4);
});
test('lithology grid: land and unmapped sea give null, a coastal gap takes the closest node with data', () => {
  for (const [lat, lon] of [[23, 10] /* Sahara */, [-75, 0], [-80, -60] /* Antarctica */, [85, 0] /* unmapped Arctic */, [90, 0], [NaN, 0], [0, Infinity], [95, 0]]) assert.equal(seabedGridAt(GRID, lat, lon), null, `${lat}, ${lon}`);
  // Lagos shoreline: the node itself is empty in the source, its neighbour 0.1° away is clay
  const r = seabedGridAt(GRID, 6.5, 3.4);
  assert.equal(r.code, 4); assert.equal(r.shifted, true); assert.ok(r.offsetKm > 10.5 && r.offsetKm < 11.5);
  assert.equal(seabedGridAt(GRID, 6.5, 3.4, { reach: 0 }), null);
});
test('lithology grid: every class maps to a seabed class with an engineering meaning', () => {
  assert.deepEqual(Object.keys(SEABED_GRID_CLASSES).map(Number), [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]);
  for (const [code, c] of Object.entries(SEABED_GRID_CLASSES)) { assert.ok(c.name, code); assert.ok(SEABED_CLASSES[c.cls], `${code} → ${c.cls}`); }
  for (const k of ['calcOoze', 'mixedOoze', 'calcMud', 'biogenic']) assert.match(SEABED_CLASSES[k].means, /crush/, `${k} warns about crushable carbonate grains`);
  assert.ok(SEABED_CLASSES.clay.d50 < SEABED_CLASSES.silt.d50 && SEABED_CLASSES.silt.d50 < 6.3e-5 && SEABED_CLASSES.volcanic.d50 > 6.3e-5);
  assert.ok(SHELF_DEPTH >= 100 && SHELF_DEPTH <= 300);
});
