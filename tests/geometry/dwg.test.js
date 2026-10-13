import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, readdirSync, statSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import * as dwg from '../../src/geometry/adapters/dwg.js';
import { importFile, formatFor, adapterStatus, GeometryFormatError, FORMATS } from '../../src/geometry/registry.js';
import { enc, close, assertContract } from './fixtures/helpers.js';

// The DWG library itself is fetched from its CDN by a browser worker and cannot run here, so these tests cover
// what belongs to the application: the consent gate, the message protocol of the worker, the conversion of the
// drawing description into the DXF reader's group codes, and the separation itself (no library code in the
// repository, no import of it from application modules).

const root = fileURLToPath(new URL('../../', import.meta.url));
const DWG = (v = 'AC1032') => enc(`${v}\0\0\0\0\0\0binary`);

/** A stand-in for the worker: answers the "read" request with a prepared drawing description. */
function fakeWorker(drawing, log = []) {
  return () => {
    const w = {
      terminated: false,
      postMessage(m, transfer) {
        log.push({ cmd: m.cmd, consent: m.consent, url: m.url, bytes: m.bytes.byteLength, transferred: Array.isArray(transfer) && transfer[0] === m.bytes });
        queueMicrotask(() => w.onmessage({ data: typeof drawing === 'string' ? { id: m.id, ok: false, error: drawing } : { id: m.id, ok: true, drawing } }));
      },
      terminate() { w.terminated = true; log.push('terminated'); },
    };
    return w;
  };
}

const P3 = (x, y, z = 0) => ({ x, y, z });
const DRAWING = {
  header: { INSUNITS: 4, ACADVER: 'AC1032' },
  layers: ['0', 'LEGS', 'PLATING'],
  blocks: [
    { name: '*Model_Space', basePoint: P3(0, 0), entities: [] },
    { name: 'BRACE', basePoint: P3(0, 0), entities: [{ type: 'LINE', layer: '0', startPoint: P3(0, 0, 0), endPoint: P3(1000, 0, 1000) }] },
  ],
  entities: [
    { type: 'LINE', layer: 'LEGS', startPoint: P3(0, 0, 0), endPoint: P3(0, 0, 10000) },
    { type: 'LINE', layer: 'LEGS', startPoint: P3(8000, 0, 0), endPoint: P3(8000, 0, 10000) },
    { type: 'INSERT', layer: 'BRACES', name: 'BRACE', insertionPoint: P3(0, 0, 2000), xScale: 8, yScale: 1, zScale: 1, rotation: 0, columnCount: 1, rowCount: 1, columnSpacing: 0, rowSpacing: 0, extrusionDirection: P3(0, 0, 1) },
    { type: '3DFACE', layer: 'PLATING', corner1: P3(0, 0, 10000), corner2: P3(8000, 0, 10000), corner3: P3(8000, 4000, 10000), corner4: P3(0, 4000, 10000), flag: 0 },
    // quarter arc of radius 2000 about (0, 0, 12000): angles in radians, as the library gives them
    { type: 'ARC', layer: 'DECK', center: P3(0, 0, 12000), radius: 2000, startAngle: 0, endAngle: Math.PI / 2, extrusionDirection: P3(0, 0, 1) },
    { type: 'LWPOLYLINE', layer: 'DECK', flag: 1, elevation: 12000, vertices: [{ x: 0, y: 0, bulge: 0 }, { x: 4000, y: 0, bulge: 0 }, { x: 4000, y: 3000, bulge: 0 }, { x: 0, y: 3000, bulge: 0 }], extrusionDirection: P3(0, 0, 1) },
    { type: 'POLYLINE3D', layer: 'ROUTE', flag: 0, vertices: [{ x: 0, y: 0, z: 0, flag: 32 }, { x: 0, y: 5000, z: -1000, flag: 32 }, { x: 0, y: 9000, z: -3000, flag: 32 }] },
    { type: 'LINE', layer: 'SHEET', isInPaperSpace: true, startPoint: P3(0, 0, 0), endPoint: P3(1e6, 0, 0) },
    { type: 'MTEXT', layer: '0', text: 'note' },
  ],
};

test('DWG adapter: nothing is read without the consent of the user, and the refusal says what would happen', async () => {
  assert.equal(dwg.id, 'dwg');
  assert.deepEqual(dwg.formats, ['dwg']);
  assert.equal(dwg.library.licence, 'GPL-3.0');
  assert.equal(dwg.library.bundled, false);
  assert.deepEqual(dwg.requirements, { consent: true, network: 'first use', environment: 'browser (module worker)' });
  // the library address is one constant: its own distribution on a public CDN, at an exact version
  assert.match(dwg.LIBRARY_URL, /^https:\/\/cdn\.jsdelivr\.net\/npm\/@mlightcad\/libredwg-web@\d+\.\d+\.\d+\/dist\/libredwg-web\.js$/);
  assert.ok(dwg.LIBRARY_URL.includes(`@${dwg.library.version}/`));
  const log = [];
  const input = { name: 'jacket.dwg', bytes: DWG() };
  for (const opts of [{}, { consent: false }, { consent: 'yes' }, { consent: 1 }]) {
    await assert.rejects(dwg.parse(input, { ...opts, worker: fakeWorker(DRAWING, log) }), /jacket\.dwg: Reading a DWG drawing downloads the LibreDWG library \(GNU GPL version 3, about 11 MB\) from the public CDN cdn\.jsdelivr\.net into your browser and runs it there in a separate worker\. The library is not part of this application\. Confirm to continue \(pass consent: true\), or save the drawing as DXF instead\./);
  }
  assert.deepEqual(log, [], 'no worker was started and nothing was sent');
  // outside a browser there is no worker to load the library into
  assert.equal(dwg.available(), false);
  assert.match(dwg.unavailableReason(), /runs in a web browser only/);
  await assert.rejects(dwg.parse(input, { consent: true }), /jacket\.dwg: DWG reading runs in a web browser only .* Save the drawing as DXF/);
});

test('DWG adapter: with consent the bytes go to the worker and the drawing that comes back is read like a DXF drawing', async () => {
  const log = [];
  const bytes = DWG('AC1027');
  const r = await dwg.parse({ name: 'jacket.dwg', bytes }, { consent: true, worker: fakeWorker(DRAWING, log) });
  assert.deepEqual(log, [{ cmd: 'read', consent: true, url: dwg.LIBRARY_URL, bytes: bytes.length, transferred: true }, 'terminated']);
  assertContract(r, ['surface', 'frame', 'polyline', 'section']);
  assert.equal(r.meta.format, 'dwg');
  assert.equal(r.meta.dwgVersion, 'AC1027');
  assert.deepEqual(r.meta.units, { unit: 'mm', scale: 0.001, basis: 'file' });
  assert.deepEqual(r.meta.library, { name: dwg.library.name, version: dwg.library.version, licence: 'GPL-3.0', loadedFrom: dwg.LIBRARY_URL });
  assert.deepEqual(r.meta.elementTypes, { LINE: 3, INSERT: 1, '3DFACE': 1, ARC: 1, LWPOLYLINE: 1, POLYLINE: 1 });
  assert.deepEqual(r.meta.skipped, { MTEXT: 1 });
  assert.deepEqual(r.meta.layers, ['0', 'LEGS', 'PLATING']);
  // two legs and the brace of the block, stretched to the leg spacing by the insert; paper space is left out
  const frame = r.models.find((m) => m.kind === 'frame');
  assert.equal(frame.members.length, 3);
  assert.deepEqual(frame.members.map((m) => m.group).sort(), ['BRACES', 'LEGS', 'LEGS']);
  close(r.meta.bbox.max[0], 8, 1e-9); close(r.meta.bbox.max[2], 12, 1e-9); close(r.meta.bbox.min[2], -3, 1e-9);
  const deck = r.models.find((m) => m.kind === 'surface');
  assert.equal(deck.indices.length / 3, 2);
  // the arc sweeps a quarter turn: it ends on the y axis, 2 m from the centre
  const arc = r.models.find((m) => m.kind === 'polyline' && m.entity === 'ARC');
  const n = arc.points.length;
  close(arc.points[0], 2, 1e-9); close(arc.points[1], 0, 1e-9);
  close(arc.points[n - 3], 0, 1e-9); close(arc.points[n - 2], 2, 1e-9); close(arc.points[n - 1], 12, 1e-9);
  assert.ok(r.models.some((m) => m.kind === 'section' && m.loops.length === 1));
  assert.ok(r.warnings.some((w) => /read with @mlightcad\/libredwg-web .* \(GPL-3\.0\), downloaded into the browser at your request; it is not part of this application/.test(w)));
  // the same through the registry: the format is an adapter format and the consent option is passed on
  const viaRegistry = await importFile({ name: 'jacket.dwg', bytes: DWG() }, { consent: true, adapters: { dwg: { parse: (i, o) => dwg.parse(i, { ...o, worker: fakeWorker(DRAWING) }) } } });
  assert.equal(viaRegistry.format, 'dwg');
  assert.equal(viaRegistry.meta.via, 'adapter:dwg');
  assert.equal(viaRegistry.meta.detected.reason, 'DWG version string');
});

test('DWG adapter: failures of the library or of the network are reported with what to do', async () => {
  const input = { name: 'plan.dwg', bytes: DWG() };
  await assert.rejects(dwg.parse(input, { consent: true, worker: fakeWorker('Failed to fetch dynamically imported module: https://cdn.jsdelivr.net/…') }), /The DWG library could not be loaded from https:\/\/cdn\.jsdelivr\.net\/npm\/@mlightcad\/libredwg-web@.* A network connection is needed the first time; alternatively save the drawing as DXF\./);
  await assert.rejects(dwg.parse(input, { consent: true, worker: fakeWorker('the library could not open the data as a DWG drawing') }), /The DWG library could not read the drawing \(the library could not open the data as a DWG drawing\)\. Save it as DXF/);
  await assert.rejects(dwg.parse(input, { consent: true, worker: () => { throw new Error('blocked'); } }), /The DWG worker could not be started \(blocked\)/);
  // a worker that never answers is given up on and terminated
  const log = [];
  const silent = () => ({ postMessage() { log.push('sent'); }, terminate() { log.push('terminated'); } });
  await assert.rejects(dwg.parse(input, { consent: true, worker: silent, timeoutMs: 20 }), /did not answer in time/);
  assert.deepEqual(log, ['sent', 'terminated']);
  await assert.rejects(dwg.parse({ name: 'x.dwg', bytes: enc('not a drawing') }, { consent: true, worker: fakeWorker(DRAWING) }), /x\.dwg: not a DWG drawing \(the "AC10xx" version string at the start is missing\)/);
  await assert.rejects(dwg.parse(input, { consent: true, worker: fakeWorker({ header: {}, layers: [], blocks: [], entities: [] }) }), /plan\.dwg: the drawing holds no entities/);
});

test('DWG in the registry: an adapter format whose library is not in the build, and says so', async () => {
  const f = formatFor('plan.dwg', DWG());
  assert.equal(f.id, 'dwg');
  assert.equal(f.status, 'adapter');
  assert.equal(f.via, 'adapter:dwg');
  assert.match(f.notes, /LibreDWG library \(GNU GPL 3\), which is NOT included in this application/);
  assert.match(f.notes, /Needs a network connection the first time/);
  assert.match(f.notes, /opts\.consent = true/);
  assert.match(f.missing, /save the drawing as ASCII or binary DXF/);
  // not one of the readers the build contains
  assert.ok(!('dwg' in (await adapterStatus())));
  // in Node the adapter is unavailable: the registry says why and what to do instead
  await assert.rejects(importFile({ name: 'plan.dwg', bytes: DWG() }, { consent: true }), (e) => e instanceof GeometryFormatError && e.reason === 'adapter-missing' && /AutoCAD DWG needs the "dwg" reader .* which runs in a web browser only .* save the drawing as ASCII or binary DXF/.test(e.message));
  assert.ok(FORMATS.find((x) => x.id === 'dwg').kinds.includes('frame'));
});

test('DWG separation: the repository holds no copy of the library and no application module imports it', () => {
  // the worker is a message loop: one dynamic import of the address it is handed, nothing else from outside
  const worker = readFileSync(`${root}src/geometry/adapters/dwg.worker.js`, 'utf8');
  assert.equal((worker.match(/\bimport\s*\(/g) || []).length, 1);
  assert.match(worker, /await import\(m\.url\)/);
  assert.ok(!/^\s*import\s/m.test(worker), 'no static import');
  assert.ok(!/libredwg|jsdelivr/i.test(worker.replace(/^\/\/.*$/gm, '').replace(/lib\.LibreDwg|lib\.Dwg_File_Type/g, '')), 'the worker does not name the library address itself');
  assert.match(worker, /if \(m\.consent !== true\) throw/);
  // the address appears in exactly one place of the source tree
  const walk = (dir) => readdirSync(dir).flatMap((n) => { const p = `${dir}/${n}`; return statSync(p).isDirectory() ? walk(p) : [p]; });
  const src = walk(`${root}src`).filter((p) => /\.(js|mjs|html|json)$/.test(p));
  const naming = src.filter((p) => /cdn\.jsdelivr\.net\/npm\/@mlightcad/.test(readFileSync(p, 'utf8')));
  assert.deepEqual(naming.map((p) => p.slice(root.length)), ['src/geometry/adapters/dwg.js']);
  // nothing of the library is stored: no vendored directory, no file carrying its name
  const everything = walk(`${root}src`).concat(walk(`${root}vendor`));
  assert.deepEqual(everything.filter((p) => /libredwg/i.test(p)), []);
  for (const p of src) {
    const text = readFileSync(p, 'utf8');
    assert.ok(!/from\s+['"][^'"]*libredwg[^'"]*['"]/i.test(text), `${p} imports the library`);
  }
});
