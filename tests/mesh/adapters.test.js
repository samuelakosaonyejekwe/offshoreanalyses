import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile, readdir } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import * as occt from '../../src/geometry/adapters/occt.js';
import * as hdf5 from '../../src/geometry/adapters/hdf5.js';
import * as ifc from '../../src/geometry/adapters/web-ifc.js';
import * as laz from '../../src/geometry/adapters/laz.js';
import * as sqlite from '../../src/geometry/adapters/sqlite.js';
import { adapterStatus } from '../../src/geometry/registry.js';
import { inspect } from '../../src/mesh/topology.js';
import { inertiaSolid } from '../../src/mesh/measure.js';
import { repair } from '../../src/mesh/repair.js';
import { buildItem } from '../../src/ui/viewer3d.js';
import { close, rel, closeV } from './helpers.js';

const fixture = (name) => readFile(new URL(`./fixtures/${name}`, import.meta.url));
const vendor = (name) => readFile(new URL(`../../vendor/occt/${name}`, import.meta.url));

test('vendored CAD kernel matches the record in docs/vendor.md (size and SHA-256)', async () => {
  const doc = await readFile(new URL('../../docs/vendor.md', import.meta.url), 'utf8');
  let total = 0;
  for (const f of ['occt-import-js.js', 'occt-import-js.wasm']) {
    const buf = await vendor(f);
    total += buf.length;
    const sha = createHash('sha256').update(buf).digest('hex');
    assert.ok(doc.includes(sha), `${f}: SHA-256 ${sha} is recorded`);
    assert.ok(doc.includes(String(buf.length)), `${f}: size ${buf.length} is recorded`);
  }
  assert.equal(occt.sizeBytes, total);
  assert.ok(doc.includes(occt.library.version) && doc.includes(occt.library.licence));
  assert.ok((await vendor('license.occt.txt')).length > 1000 && (await vendor('license.occt-import-js.txt')).length > 1000);
  assert.equal(occt.available(), true);
});

test('STEP block 2000 × 3000 × 4000 mm → one watertight solid of 2 × 3 × 4 m with six face groups', async () => {
  const bytes = new Uint8Array(await fixture('block.step'));
  assert.equal(occt.sniff({ name: 'unknown.dat', bytes }), 'step');
  const r = await occt.parse({ name: 'block.step', bytes });
  assert.equal(r.models.length, 1);
  assert.deepEqual(r.warnings, []);
  assert.equal(r.meta.format, 'step');
  assert.deepEqual(r.meta.units, { unit: 'm', scale: 1, basis: 'file' });
  assert.deepEqual(r.meta.bbox, { min: [0, 0, 0], max: [2, 3, 4] });
  assert.equal(r.meta.assembly.children[0].name, 'block');
  const m = r.models[0];
  assert.equal(m.kind, 'surface'); assert.equal(m.name, 'block');
  assert.ok(m.positions instanceof Float32Array && m.indices instanceof Uint32Array);
  assert.equal(m.indices.length / 3, 12);
  assert.equal(m.groups.length, 6);
  assert.equal(m.groups.reduce((s, g) => s + g.count, 0), m.indices.length);
  const i = inspect(m);
  assert.ok(i.ok && i.watertight && i.stats.outward);
  assert.equal(i.stats.eulerCharacteristic, 2);
  rel(i.stats.volume, 24, 1e-6); rel(i.stats.area, 52, 1e-6);
  closeV(inertiaSolid(m.positions, m.indices).centroid, [1, 1.5, 2], 1e-5);
  // faces are tessellated separately; welding gives the eight corners back
  assert.equal(m.positions.length / 3, 24);
  assert.equal(repair(m).model.positions.length / 3, 8);
  assert.equal(buildItem(m).tris.length, 36); // displays in the viewer as it is
});

test('import options: text input, y-up conversion, explicit format; unreadable content fails clearly', async () => {
  const text = (await fixture('block.step')).toString('utf8');
  const y = await occt.parse({ name: 'block.stp', text }, { upAxis: 'y', quality: 'coarse', format: 'step' });
  const i = inspect(y.models[0]);
  rel(i.stats.volume, 24, 1e-6);
  // file z (4 m) is not the vertical any more
  close(i.stats.bbox.max[2] - i.stats.bbox.min[2], 3, 1e-5);
  assert.equal(y.meta.tessellation.linearDeflection, 0.005);
  await assert.rejects(() => occt.parse({ name: 'notes.txt', text: 'hello' }), /not recognised/);
  await assert.rejects(() => occt.parse({ name: 'broken.step', text: 'ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\nENDSEC;\nEND-ISO-10303-21;\n' }), /could not be read|no surfaces/);
  assert.equal(occt.sniff({ name: 'a.igs', text: 'x' }), 'iges');
  assert.equal(occt.sniff({ name: 'a.brep', text: 'DBRep_DrawableShape\n' }), 'brep');
  assert.equal(occt.sniff({ name: 'a.xyz', text: 'x' }), null);
});

test('every adapter of this build is real: common interface, runtime on disk and recorded in docs/vendor.md', async () => {
  const doc = await readFile(new URL('../../docs/vendor.md', import.meta.url), 'utf8');
  const all = { occt, hdf5, 'web-ifc': ifc, laz, sqlite };
  for (const [name, a] of Object.entries(all)) {
    assert.equal(a.id, name);
    assert.equal(typeof a.title, 'string');
    assert.ok(Array.isArray(a.formats) && a.formats.length && Array.isArray(a.extensions) && a.extensions.length, name);
    assert.equal(typeof a.library.name, 'string');
    assert.equal(a.available(), true, `${name} is available`);
    assert.equal(typeof a.parse, 'function');
    // input without content is refused with the file name in front, not with a crash inside a runtime
    await assert.rejects(() => a.parse({ name: 'survey.bin' }), /^Error: survey\.bin: /, name);
    await assert.rejects(() => a.parse({ name: 'survey.bin', bytes: new Uint8Array(4) }), /^Error: survey\.bin: .*not (an? |recognised)/, name);
  }
  // runtimes that are downloaded on first use state their size, version and licence, as recorded
  for (const a of [occt, hdf5, ifc, laz]) {
    assert.ok(a.sizeBytes > 100000, a.id);
    assert.match(a.library.version, /^\d+\.\d+\.\d+$/);
    assert.ok(doc.includes(`| ${a.library.version} |`) && doc.includes(a.library.licence), `${a.id}: version and licence are in docs/vendor.md`);
    assert.ok(doc.includes(a.sizeBytes.toLocaleString('en-US').replace(/,/g, ' ')), `${a.id}: download size ${a.sizeBytes} is in docs/vendor.md`);
    assert.equal(typeof a.load, 'function');
  }
  // the GeoPackage reader is part of the application: nothing to download
  assert.equal(sqlite.sizeBytes, 0);
  // what to convert a file to when a runtime cannot start (no WebAssembly) stays available as advice
  for (const a of [hdf5, ifc, laz, sqlite]) assert.ok(typeof a.convertTo === 'string' && a.convertTo.length > 10, a.id);
  const st = await adapterStatus();
  assert.deepEqual(Object.fromEntries(Object.entries(st).map(([k, v]) => [k, v.available])), { occt: true, 'web-ifc': true, hdf5: true, laz: true, sqlite: true });
});

test('no vendored file is left unrecorded: every file under vendor/ has its SHA-256 in docs/vendor.md', async () => {
  const doc = await readFile(new URL('../../docs/vendor.md', import.meta.url), 'utf8');
  const root = new URL('../../vendor/', import.meta.url);
  const walk = async (dir) => (await Promise.all((await readdir(dir, { withFileTypes: true })).map((e) => (e.isDirectory() ? walk(new URL(`${e.name}/`, dir)) : [new URL(e.name, dir)])))).flat();
  const files = await walk(root);
  assert.deepEqual(files.map((u) => u.href.slice(root.href.length)).sort(), [
    'h5wasm/LICENSE.txt', 'h5wasm/hdf5_hl.js', 'h5wasm/hdf5_util.js',
    'laz-perf/LICENSE.txt', 'laz-perf/laz-perf.js', 'laz-perf/laz-perf.wasm', 'laz-perf/package.json',
    'occt/license.occt-import-js.txt', 'occt/license.occt.txt', 'occt/occt-import-js.js', 'occt/occt-import-js.wasm', 'occt/package.json',
    'web-ifc/LICENSE.md', 'web-ifc/node/package.json', 'web-ifc/node/web-ifc-api-node.js', 'web-ifc/node/web-ifc-node.wasm', 'web-ifc/web-ifc-api.js', 'web-ifc/web-ifc.wasm',
  ]);
  for (const u of files) {
    const buf = await readFile(u);
    assert.ok(doc.includes(createHash('sha256').update(buf).digest('hex')), `${u.href.slice(root.href.length)}: SHA-256 is recorded`);
    assert.ok(doc.includes(`| ${buf.length} |`), `${u.href.slice(root.href.length)}: size ${buf.length} is recorded`);
  }
});
