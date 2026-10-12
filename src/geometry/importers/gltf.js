// glTF 2.0 (.gltf JSON with embedded or side-car buffers, and binary .glb) and 3MF packages.
//
// glTF: the default scene is flattened — node transforms are applied, every mesh primitive becomes a group
// of one surface model. glTF is metres, y up, by specification; the result is rotated to z up.
// Not read: Draco / meshopt compressed primitives, morph targets, skins, animations, materials, textures.
//
// 3MF: the core specification — objects with meshes, components, build items with transforms, the `unit`
// attribute (millimetre by default) — plus production-extension `path` references to other model parts.
// Not read: beam lattices, slices, materials, volumetric data.

import {
  inputText, baseName, Grow, finish, surfaceModel, pointsModel, M4, findFile, base64Decode,
  decodeText, cleanName, groupsFromBuckets, tally,
} from './common.js';
import { resolveUnits, scaleArray, yUpToZUp } from '../units.js';
import { isZip, readZip } from '../zip.js';
import { parseXml, kids, kid, attr, numberList } from '../xml.js';

const COMP = { 5120: [Int8Array, 1], 5121: [Uint8Array, 1], 5122: [Int16Array, 2], 5123: [Uint16Array, 2], 5125: [Uint32Array, 4], 5126: [Float32Array, 4] };
const NCOMP = { SCALAR: 1, VEC2: 2, VEC3: 3, VEC4: 4, MAT2: 4, MAT3: 9, MAT4: 16 };
const NORM = { 5120: (v) => Math.max(v / 127, -1), 5121: (v) => v / 255, 5122: (v) => Math.max(v / 32767, -1), 5123: (v) => v / 65535 };
const GETTER = { 5120: 'getInt8', 5121: 'getUint8', 5122: 'getInt16', 5123: 'getUint16', 5125: 'getUint32', 5126: 'getFloat32' };

/** glTF 2.0 JSON or GLB. */
export function parseGltf(input, opts = {}) {
  const name = baseName(input.name);
  const warnings = [];
  let json, binChunk = null, isGlb = false;
  const bytes = input.bytes || null;
  if (bytes && bytes.length >= 12 && bytes[0] === 0x67 && bytes[1] === 0x6c && bytes[2] === 0x54 && bytes[3] === 0x46) {
    const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
    const version = dv.getUint32(4, true);
    isGlb = true;
    if (version !== 2) throw new Error(`${input.name || 'GLB'}: binary glTF version ${version} is not supported (expected version 2; re-export as glTF 2.0).`);
    let p = 12;
    while (p + 8 <= bytes.length) {
      const len = dv.getUint32(p, true), type = dv.getUint32(p + 4, true);
      const data = bytes.subarray(p + 8, p + 8 + len);
      if (p + 8 + len > bytes.length) throw new Error(`${input.name || 'GLB'}: a chunk of ${len} bytes runs past the end of the file (${bytes.length} bytes) — the GLB is truncated.`);
      if (type === 0x4e4f534a) { try { json = JSON.parse(decodeText(data)); } catch (e) { throw new Error(`${input.name || 'GLB'}: the JSON chunk of the binary glTF is not valid JSON (${e.message}).`); } } else if (type === 0x004e4942 && !binChunk) binChunk = data;
      p += 8 + len;
    }
    if (!json) throw new Error(`${input.name || 'GLB'}: no JSON chunk found in the binary glTF container.`);
  } else {
    try { json = JSON.parse(inputText(input)); } catch (e) { throw new Error(`${input.name || 'glTF'}: not valid JSON (${e.message}); expected a glTF 2.0 document with "asset", "meshes" and "accessors".`); }
  }
  if (!json || !json.asset) throw new Error(`${input.name || 'glTF'}: missing the "asset" member (expected a glTF 2.0 document).`);
  if (json.asset.version && !String(json.asset.version).startsWith('2')) throw new Error(`${input.name || 'glTF'}: glTF version ${json.asset.version} is not supported (expected 2.x).`);
  const required = json.extensionsRequired || [];
  const blocked = required.filter((e) => /draco|meshopt/i.test(e));

  const buffers = (json.buffers || []).map((b, i) => {
    if (b.uri === undefined) {
      if (i === 0 && binChunk) return binChunk;
      throw new Error(`${input.name || 'glTF'}: buffer ${i} has no uri and the file has no binary chunk.`);
    }
    if (b.uri.startsWith('data:')) { const c = b.uri.indexOf(','); return base64Decode(b.uri, c + 1); }
    const f = findFile(input.files, decodeURIComponent(b.uri));
    if (!f) throw new Error(`${input.name || 'glTF'}: buffer file "${b.uri}" was not supplied. Upload the .gltf together with its .bin files (or zipped), or export a single .glb.`);
    return f;
  });

  const accessorCache = new Map();
  const readAccessor = (ai) => {
    if (accessorCache.has(ai)) return accessorCache.get(ai);
    const a = json.accessors && json.accessors[ai];
    if (!a) throw new Error(`${input.name || 'glTF'}: accessor ${ai} does not exist.`);
    const ct = COMP[a.componentType], nc = NCOMP[a.type];
    if (!ct || !nc) throw new Error(`${input.name || 'glTF'}: accessor ${ai} has an unknown componentType/type (${a.componentType}/${a.type}).`);
    const out = new Float64Array(a.count * nc);
    const fill = (view, byteOffset, count, dst, dstIndex, ctype, ncomp) => {
      const bv = json.bufferViews[view];
      const buf = buffers[bv.buffer];
      const size = COMP[ctype][1];
      const stride = bv.byteStride || size * ncomp;
      const base = (bv.byteOffset || 0) + (byteOffset || 0);
      if (base + (count - 1) * stride + size * ncomp > buf.length) throw new Error(`${input.name || 'glTF'}: accessor data runs past the end of buffer ${bv.buffer} (the file is truncated or a side-car .bin is the wrong one).`);
      const dv = new DataView(buf.buffer, buf.byteOffset, buf.byteLength);
      const get = GETTER[ctype];
      for (let i = 0; i < count; i++) for (let c = 0; c < ncomp; c++) dst[dstIndex(i) * ncomp + c] = dv[get](base + i * stride + c * size, true);
    };
    if (a.bufferView !== undefined) fill(a.bufferView, a.byteOffset, a.count, out, (i) => i, a.componentType, nc);
    if (a.sparse) {
      const s = a.sparse;
      const ids = new Float64Array(s.count);
      fill(s.indices.bufferView, s.indices.byteOffset, s.count, ids, (i) => i, s.indices.componentType, 1);
      fill(s.values.bufferView, s.values.byteOffset, s.count, out, (i) => ids[i], a.componentType, nc);
    }
    if (a.normalized && NORM[a.componentType]) { const f = NORM[a.componentType]; for (let i = 0; i < out.length; i++) out[i] = f(out[i]); }
    accessorCache.set(ai, out);
    return out;
  };

  const coords = new Grow(Float64Array, 3 * 4096);
  const buckets = new Map();
  const polylines = [], pointSets = [];
  const hist = {}, skipped = {};
  const uniq = (base, map) => { let n = base, k = 2; while (map.has(n)) n = `${base} (${k++})`; return n; };

  const emitMesh = (mi, M, nodeName) => {
    const mesh = json.meshes && json.meshes[mi];
    if (!mesh) return;
    (mesh.primitives || []).forEach((prim, pi) => {
      if (prim.extensions && (prim.extensions.KHR_draco_mesh_compression || prim.extensions.EXT_meshopt_compression)) { tally(skipped, 'compressedPrimitive'); return; }
      if (!prim.attributes || prim.attributes.POSITION === undefined) { tally(skipped, 'primitiveWithoutPositions'); return; }
      let pos;
      try { pos = readAccessor(prim.attributes.POSITION); } catch (e) { if (blocked.length) { tally(skipped, 'compressedPrimitive'); return; } throw e; }
      const n = pos.length / 3;
      const idx = prim.indices !== undefined ? readAccessor(prim.indices) : null;
      const count = idx ? idx.length : n;
      const at = (i) => (idx ? idx[i] : i);
      const mode = prim.mode === undefined ? 4 : prim.mode;
      const world = new Float64Array(pos.length);
      for (let i = 0; i < n; i++) M4.apply(M, pos[3 * i], pos[3 * i + 1], pos[3 * i + 2], world, 3 * i);
      const label = cleanName(nodeName || mesh.name, `mesh ${mi}`) + (mesh.primitives.length > 1 ? ` [${pi}]` : '');
      if (mode >= 4 && mode <= 6) {
        const base = coords.n / 3;
        coords.ensure(world.length);
        coords.a.set(world, coords.n); coords.n += world.length;
        const g = new Grow(Uint32Array, Math.max(3, count));
        // a negative determinant mirrors the part: flip the winding so that outward normals stay outward
        const det = M[0] * (M[5] * M[10] - M[6] * M[9]) - M[4] * (M[1] * M[10] - M[2] * M[9]) + M[8] * (M[1] * M[6] - M[2] * M[5]);
        const tri = (a, b, c) => { if (a === b || b === c || a === c) return; g.push(base + a); if (det < 0) { g.push(base + c); g.push(base + b); } else { g.push(base + b); g.push(base + c); } };
        if (mode === 4) for (let i = 0; i + 2 < count; i += 3) tri(at(i), at(i + 1), at(i + 2));
        else if (mode === 5) for (let i = 0; i + 2 < count; i++) { if (i & 1) tri(at(i + 1), at(i), at(i + 2)); else tri(at(i), at(i + 1), at(i + 2)); } else for (let i = 1; i + 1 < count; i++) tri(at(0), at(i), at(i + 1));
        tally(hist, 'tri3', g.n / 3);
        buckets.set(uniq(label, buckets), g);
      } else if (mode === 0) {
        const p = new Float64Array(count * 3);
        for (let i = 0; i < count; i++) { const k = at(i); p[3 * i] = world[3 * k]; p[3 * i + 1] = world[3 * k + 1]; p[3 * i + 2] = world[3 * k + 2]; }
        pointSets.push({ name: label, p });
        tally(hist, 'point1', count);
      } else if (mode >= 1 && mode <= 3) {
        const seq = (list) => { const p = new Float64Array(list.length * 3); list.forEach((k, i) => { p[3 * i] = world[3 * k]; p[3 * i + 1] = world[3 * k + 1]; p[3 * i + 2] = world[3 * k + 2]; }); polylines.push({ name: label, p }); };
        if (mode === 1) for (let i = 0; i + 1 < count; i += 2) seq([at(i), at(i + 1)]);
        else { const l = []; for (let i = 0; i < count; i++) l.push(at(i)); if (mode === 2 && l.length) l.push(l[0]); seq(l); }
        tally(hist, 'line2', mode === 1 ? count >> 1 : count - (mode === 3 ? 1 : 0));
      } else tally(skipped, `primitiveMode${mode}`);
    });
  };

  const nodes = json.nodes || [];
  const seen = new Set();
  const walk = (ni, parent, depth) => {
    const nd = nodes[ni];
    if (!nd || depth > 64) return;
    seen.add(ni);
    const local = nd.matrix ? Float64Array.from(nd.matrix) : M4.fromTRS(nd.translation, nd.rotation, nd.scale);
    const M = M4.mul(parent, local);
    if (nd.mesh !== undefined) emitMesh(nd.mesh, M, nd.name);
    for (const c of nd.children || []) walk(c, M, depth + 1);
  };
  const scene = json.scenes && json.scenes[json.scene ?? 0];
  const I = M4.identity();
  if (scene && scene.nodes) for (const r of scene.nodes) walk(r, I, 0);
  else if (nodes.length) { const child = new Set(); for (const nd of nodes) for (const c of nd.children || []) child.add(c); nodes.forEach((_, i) => { if (!child.has(i)) walk(i, I, 0); }); } else (json.meshes || []).forEach((_, i) => emitMesh(i, I, null));

  if (skipped.compressedPrimitive) warnings.push(`${skipped.compressedPrimitive} primitive(s) use Draco or meshopt compression, which is not decoded here; re-export the model without mesh compression.`);
  for (const k in skipped) if (k !== 'compressedPrimitive') warnings.push(`Skipped ${skipped[k]} ${k}.`);
  if (!coords.n && !polylines.length && !pointSets.length) {
    if (skipped.compressedPrimitive) throw new Error(`${input.name || 'glTF'}: every mesh primitive is Draco/meshopt compressed, which is not decoded. Re-export without mesh compression (for example "gltf-transform copy" without --draco, or Blender with Compression off).`);
    throw new Error(`${input.name || 'glTF'}: the scene contains no mesh geometry (expected at least one mesh primitive with a POSITION accessor).`);
  }
  const units = resolveUnits(opts, 'm', 'm');
  const fix = (c) => { scaleArray(c, units.scale); return opts.upAxis === 'z' ? c : yUpToZUp(c); };
  const models = [];
  if (coords.n) {
    const { indices, groups } = groupsFromBuckets(buckets);
    models.push(surfaceModel(name, fix(coords.finish()), indices, groups.length > 1 ? groups : null));
  }
  for (const l of polylines) models.push({ kind: 'polyline', name: l.name, points: fix(l.p) });
  for (const s of pointSets) models.push(pointsModel(s.name, fix(s.p)));
  return finish(models, warnings, {
    format: isGlb ? 'glb' : 'gltf', units, upAxis: { source: opts.upAxis === 'z' ? 'z' : 'y', basis: opts.upAxis ? 'option' : 'file' }, crs: null,
    elementTypes: hist, skipped, generator: json.asset.generator || null, nodes: seen.size, extensionsUsed: json.extensionsUsed || [],
  });
}

// ---------------------------------------------------------------------------------------------------------
// 3MF
// ---------------------------------------------------------------------------------------------------------

function tmf(s) {
  // "m00 m01 m02 m10 m11 m12 m20 m21 m22 m30 m31 m32": row vectors, so x' = x·m00 + y·m10 + z·m20 + m30
  const v = s ? numberList(s) : null;
  if (!v || v.length < 12) return M4.identity();
  const m = M4.identity();
  m[0] = v[0]; m[1] = v[1]; m[2] = v[2];
  m[4] = v[3]; m[5] = v[4]; m[6] = v[5];
  m[8] = v[6]; m[9] = v[7]; m[10] = v[8];
  m[12] = v[9]; m[13] = v[10]; m[14] = v[11];
  return m;
}

/** 3MF package (ZIP) or a bare 3D model part (.model XML). */
export async function parse3mf(input, opts = {}) {
  const name = baseName(input.name);
  const warnings = [];
  let files = input.files || null;
  let rootText = null;
  if (input.bytes && isZip(input.bytes)) files = await readZip(input.bytes, { filter: (n) => /\.(model|rels)$/i.test(n) });
  else if (input.bytes || typeof input.text === 'string') rootText = inputText(input);
  const parts = new Map(); // path → { objects: Map(id → def), build: [] , unit }
  const loadPart = (path, text) => {
    const root = parseXml(text);
    if (root.name !== 'model') throw new Error(`${input.name || '3MF'}: part "${path}" is not a 3MF model (expected root element <model>, found <${root.name}>).`);
    const part = { unit: attr(root, 'unit', 'millimeter'), objects: new Map(), build: [] };
    const res = kid(root, 'resources');
    for (const o of res ? kids(res, 'object') : []) {
      const def = { id: attr(o, 'id'), name: attr(o, 'name', ''), type: attr(o, 'type', 'model'), mesh: null, components: [] };
      const mesh = kid(o, 'mesh');
      if (mesh) {
        const vs = kid(mesh, 'vertices'), ts = kid(mesh, 'triangles');
        const vl = vs ? kids(vs, 'vertex') : [], tl = ts ? kids(ts, 'triangle') : [];
        const pos = new Float64Array(vl.length * 3), tri = new Uint32Array(tl.length * 3);
        vl.forEach((v, i) => { pos[3 * i] = +attr(v, 'x'); pos[3 * i + 1] = +attr(v, 'y'); pos[3 * i + 2] = +attr(v, 'z'); });
        tl.forEach((t, i) => { tri[3 * i] = +attr(t, 'v1'); tri[3 * i + 1] = +attr(t, 'v2'); tri[3 * i + 2] = +attr(t, 'v3'); });
        def.mesh = { pos, tri };
        for (const c of mesh.children) if (!['vertices', 'triangles'].includes(c.name)) warnings.push(`3MF mesh extension <${c.tag}> is not read.`);
      }
      const comps = kid(o, 'components');
      for (const c of comps ? kids(comps, 'component') : []) def.components.push({ objectid: attr(c, 'objectid'), path: attr(c, 'path', null), M: tmf(attr(c, 'transform')) });
      part.objects.set(def.id, def);
    }
    const build = kid(root, 'build');
    for (const it of build ? kids(build, 'item') : []) part.build.push({ objectid: attr(it, 'objectid'), path: attr(it, 'path', null), M: tmf(attr(it, 'transform')) });
    parts.set(path, part);
    return part;
  };
  const norm = (p) => '/' + String(p).replace(/\\/g, '/').replace(/^\/+/, '');
  let rootPath = '/3D/3dmodel.model';
  if (rootText !== null) loadPart(rootPath, rootText);
  else {
    if (!files || !files.size) throw new Error(`${input.name || '3MF'}: expected a 3MF package (ZIP archive holding 3D/3dmodel.model) or a 3MF model XML file.`);
    const rels = findFile(files, '_rels/.rels');
    if (rels) {
      const m = /Target\s*=\s*"([^"]+\.model)"/i.exec(decodeText(rels));
      if (m) rootPath = norm(m[1]);
    }
    let found = false;
    for (const [k, v] of files) if (/\.model$/i.test(k)) { loadPart(norm(k), decodeText(v)); if (norm(k).toLowerCase() === rootPath.toLowerCase()) { rootPath = norm(k); found = true; } }
    if (!parts.size) throw new Error(`${input.name || '3MF'}: the package holds no .model part (expected 3D/3dmodel.model).`);
    if (!found) rootPath = [...parts.keys()][0];
  }
  const root = parts.get(rootPath);
  const coords = new Grow(Float64Array, 3 * 4096);
  const buckets = new Map();
  let missing = 0, badTris = 0, nTri = 0;
  const emit = (path, id, M, label, depth) => {
    const part = parts.get(path) || [...parts.entries()].find(([k]) => k.toLowerCase() === String(path).toLowerCase())?.[1];
    const def = part && part.objects.get(id);
    if (!def || depth > 32) { missing++; return; }
    const nm = cleanName(label || def.name, `object ${id}`);
    if (def.mesh) {
      const base = coords.n / 3, { pos, tri } = def.mesh, nv = pos.length / 3;
      coords.ensure(pos.length);
      for (let i = 0; i < nv; i++) { M4.apply(M, pos[3 * i], pos[3 * i + 1], pos[3 * i + 2], coords.a, coords.n); coords.n += 3; }
      let g = buckets.get(nm);
      if (!g) { g = new Grow(Uint32Array, tri.length || 3); buckets.set(nm, g); }
      for (let i = 0; i + 2 < tri.length; i += 3) {
        if (tri[i] >= nv || tri[i + 1] >= nv || tri[i + 2] >= nv) { badTris++; continue; }
        g.push(base + tri[i]); g.push(base + tri[i + 1]); g.push(base + tri[i + 2]); nTri++;
      }
    }
    for (const c of def.components) emit(c.path ? norm(c.path) : path, c.objectid, M4.mul(M, c.M), nm, depth + 1);
  };
  const items = root.build.length ? root.build : [...root.objects.values()].filter((o) => o.type === 'model').map((o) => ({ objectid: o.id, path: null, M: M4.identity() }));
  if (!root.build.length && items.length) warnings.push('The model has no <build> section; every object is placed at its own coordinates.');
  items.forEach((it, i) => {
    const p = it.path ? norm(it.path) : rootPath;
    const def = (parts.get(p) || root).objects.get(it.objectid);
    emit(p, it.objectid, it.M, items.length > 1 ? `${cleanName(def && def.name, `object ${it.objectid}`)} #${i + 1}` : null, 0);
  });
  if (!coords.n) throw new Error(`${input.name || '3MF'}: no mesh objects were found (expected <object><mesh><vertices>/<triangles> in the model part).`);
  if (missing) warnings.push(`${missing} component/build reference(s) point to objects that are not in the package.`);
  if (badTris) warnings.push(`${badTris} triangle(s) reference vertices that do not exist and were skipped.`);
  const units = resolveUnits(opts, root.unit === 'micron' ? 'um' : root.unit, 'mm');
  const c = scaleArray(coords.finish(), units.scale);
  const { indices, groups } = groupsFromBuckets(buckets);
  return finish([surfaceModel(name, c, indices, groups.length > 1 ? groups : null)], warnings, {
    format: '3mf', units, upAxis: { source: 'z', basis: 'file' }, crs: null, elementTypes: { tri3: nTri },
    skipped: { missingReferences: missing, invalidTriangles: badTris }, objects: root.objects.size, buildItems: root.build.length,
  });
}
