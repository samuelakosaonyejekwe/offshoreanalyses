/**
 * Mesh repair and frame healing. Every function returns a new model and leaves its input untouched;
 * `repair(model, opts)` runs the appropriate sequence for the model kind and lists what it did.
 *
 * Surface groups ({ name, start, count } in index units) and per-face grouping survive every operation.
 * What is not attempted: resolving self-intersections, closing large or non-planar-friendly gaps by surface
 * reconstruction, and untangling edges shared by more than two genuinely distinct sheets — these are reported
 * by `topology.inspect` and left to the user.
 *
 * References: Botsch et al., "Polygon Mesh Processing" (2010) ch. 8 (model repair); Attene, Campen & Kobbelt,
 * "Polygon mesh repairing: an application perspective", ACM Comput. Surv. 45(2) 2013.
 */
import {
  weldVertices, buildEdges, boundaryLoops, cornerFans, connectedComponents, duplicateFaces, degenerateFaces,
  orientConsistently, faceGroupIds, groupsFromIds, bounds, CELL_INFO, cellVolume,
} from './topology.js';

function withFaces(model, keep) {
  const F = (model.indices.length / 3) | 0;
  const g = faceGroupIds(model);
  const idx = [], ids = [];
  for (let f = 0; f < F; f++) {
    if (!keep[f]) continue;
    idx.push(model.indices[3 * f], model.indices[3 * f + 1], model.indices[3 * f + 2]);
    ids.push(g.ids[f]);
  }
  const out = { ...model, indices: Uint32Array.from(idx) };
  delete out.faceCell;
  if (model.groups) { const r = groupsFromIds(out.indices, Int32Array.from(ids), g.names); out.indices = r.indices; out.groups = mergeGroupMeta(model.groups, r.groups); }
  return out;
}

function mergeGroupMeta(oldGroups, newGroups) {
  const byName = new Map(oldGroups.map((g) => [g.name, g]));
  return newGroups.map((g) => ({ ...(byName.get(g.name) || {}), ...g }));
}

/** Remove zero-area faces (repeated or collinear corners). Returns { model, removed }. */
export function removeDegenerateFaces(model, relTol = 1e-7) {
  const bad = degenerateFaces(model.positions, model.indices, relTol);
  if (!bad.length) return { model, removed: 0 };
  const keep = new Uint8Array(model.indices.length / 3).fill(1);
  for (const f of bad) keep[f] = 0;
  return { model: withFaces(model, keep), removed: bad.length };
}

/** Remove faces that repeat an earlier face (same three vertices). Returns { model, removed }. */
export function removeDuplicateFaces(model) {
  const bad = duplicateFaces(model.indices, model.positions.length / 3);
  if (!bad.length) return { model, removed: 0 };
  const keep = new Uint8Array(model.indices.length / 3).fill(1);
  for (const f of bad) keep[f] = 0;
  return { model: withFaces(model, keep), removed: bad.length };
}

/** Drop vertices no face uses and renumber. Returns { model, removed, map } (map old→new, 0xFFFFFFFF if dropped). */
export function removeUnreferencedVertices(model) {
  const nV = (model.positions.length / 3) | 0;
  const map = new Uint32Array(nV).fill(0xffffffff);
  let n = 0;
  for (let i = 0; i < model.indices.length; i++) { const v = model.indices[i]; if (map[v] === 0xffffffff) map[v] = n++; }
  if (n === nV) return { model, removed: 0, map };
  const pos = new model.positions.constructor(3 * n);
  for (let v = 0; v < nV; v++) { const m = map[v]; if (m !== 0xffffffff) { pos[3 * m] = model.positions[3 * v]; pos[3 * m + 1] = model.positions[3 * v + 1]; pos[3 * m + 2] = model.positions[3 * v + 2]; } }
  const idx = new Uint32Array(model.indices.length);
  for (let i = 0; i < idx.length; i++) idx[i] = map[model.indices[i]];
  return { model: { ...model, positions: pos, indices: idx }, removed: nV - n, map };
}

/**
 * Close boundary loops of up to `maxEdges` edges: three- and four-edge holes directly, larger ones with a fan
 * from the loop centroid (adequate for small, nearly planar holes; a strongly non-planar hole gets a cone).
 * New faces are wound to match their neighbours and collected in a group named 'filled-holes' when the model
 * has groups.
 * @param {{ maxEdges?: number }} [opts] largest hole to fill (default 32 edges)
 * @returns {{ model: object, filled: number, skipped: number, addedFaces: number }}
 */
export function fillHoles(model, { maxEdges = 32 } = {}) {
  const nV = (model.positions.length / 3) | 0;
  const bl = boundaryLoops(model.indices, nV);
  const add = [], extra = [];
  let filled = 0, skipped = 0;
  bl.loops.forEach((loop, i) => {
    const n = loop.length;
    if (!bl.closed[i] || n < 3 || n > maxEdges || new Set(loop).size !== n) { skipped++; return; }
    filled++;
    if (n === 3) add.push(loop[2], loop[1], loop[0]);
    else if (n === 4) add.push(loop[1], loop[0], loop[3], loop[3], loop[2], loop[1]);
    else {
      const c = [0, 0, 0];
      for (const v of loop) for (let k = 0; k < 3; k++) c[k] += model.positions[3 * v + k] / n;
      const ci = nV + extra.length / 3;
      extra.push(c[0], c[1], c[2]);
      for (let q = 0; q < n; q++) add.push(loop[(q + 1) % n], loop[q], ci);
    }
  });
  if (!add.length) return { model, filled: 0, skipped, addedFaces: 0 };
  const pos = new model.positions.constructor(model.positions.length + extra.length);
  pos.set(model.positions); pos.set(extra, model.positions.length);
  const idx = new Uint32Array(model.indices.length + add.length);
  idx.set(model.indices); idx.set(add, model.indices.length);
  const out = { ...model, positions: pos, indices: idx };
  delete out.faceCell;
  if (model.groups) {
    const last = model.groups.find((g) => g.name === 'filled-holes');
    if (last && last.start + last.count === model.indices.length) out.groups = model.groups.map((g) => (g === last ? { ...g, count: g.count + add.length } : g));
    else out.groups = [...model.groups, { name: 'filled-holes', start: model.indices.length, count: add.length }];
  }
  return { model: out, filled, skipped, addedFaces: add.length / 3 };
}

/** Consistent, outward winding (see topology.orientConsistently). Returns { model, flipped, orientable }. */
export function fixOrientation(model) {
  const r = orientConsistently(model.positions, model.indices);
  return { model: r.flipped ? { ...model, indices: r.indices } : model, flipped: r.flipped, orientable: r.orientable };
}

/**
 * Give every edge-connected fan of faces around a bow-tie vertex its own copy of the vertex, so that the
 * surface becomes vertex-manifold. Returns { model, split } with the number of vertices added.
 */
export function splitNonManifoldVertices(model) {
  const nV = (model.positions.length / 3) | 0;
  const fans = cornerFans(model.indices, nV);
  const firstFan = new Int32Array(nV).fill(-1);
  const idx = Uint32Array.from(model.indices);
  const extra = [];
  const copies = new Map();
  for (let c = 0; c < idx.length; c++) {
    const v = model.indices[c];
    if (firstFan[v] < 0) { firstFan[v] = fans[c]; continue; }
    if (firstFan[v] === fans[c]) continue;
    let id = copies.get(fans[c]);
    if (id === undefined) {
      id = nV + extra.length / 3;
      extra.push(model.positions[3 * v], model.positions[3 * v + 1], model.positions[3 * v + 2]);
      copies.set(fans[c], id);
    }
    idx[c] = id;
  }
  if (!extra.length) return { model, split: 0 };
  const pos = new model.positions.constructor(model.positions.length + extra.length);
  pos.set(model.positions); pos.set(extra, model.positions.length);
  return { model: { ...model, positions: pos, indices: idx }, split: extra.length / 3 };
}

/**
 * Remove small disconnected parts (debris left by scanning or CAD export).
 * @param {{ minFaces?: number, minAreaFraction?: number }} [opts] a component is dropped when it has fewer than
 *   `minFaces` faces (default 4) or less than `minAreaFraction` of the total area (default 1e-4)
 * @returns {{ model: object, removedComponents: number, removedFaces: number }}
 */
export function dropSmallComponents(model, { minFaces = 4, minAreaFraction = 1e-4 } = {}) {
  const nV = (model.positions.length / 3) | 0;
  const F = (model.indices.length / 3) | 0;
  const cc = connectedComponents(model.indices, nV, 'vertex');
  if (cc.count <= 1) return { model, removedComponents: 0, removedFaces: 0 };
  const area = new Float64Array(cc.count);
  const P = model.positions;
  let total = 0;
  for (let f = 0; f < F; f++) {
    const a = 3 * model.indices[3 * f], b = 3 * model.indices[3 * f + 1], c = 3 * model.indices[3 * f + 2];
    const ux = P[b] - P[a], uy = P[b + 1] - P[a + 1], uz = P[b + 2] - P[a + 2], vx = P[c] - P[a], vy = P[c + 1] - P[a + 1], vz = P[c + 2] - P[a + 2];
    const A = Math.hypot(uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx) / 2;
    area[cc.faceComponent[f]] += A; total += A;
  }
  const drop = new Uint8Array(cc.count);
  let nd = 0;
  for (let c = 0; c < cc.count; c++) if (cc.sizes[c] < minFaces || area[c] < minAreaFraction * total) { drop[c] = 1; nd++; }
  if (!nd || nd === cc.count) return { model, removedComponents: 0, removedFaces: 0 };
  const keep = new Uint8Array(F);
  let rf = 0;
  for (let f = 0; f < F; f++) { keep[f] = drop[cc.faceComponent[f]] ? 0 : 1; if (!keep[f]) rf++; }
  return { model: withFaces(model, keep), removedComponents: nd, removedFaces: rf };
}

/**
 * Frame healing: merge near-coincident nodes, split members at nodes lying on their span, remove zero-length
 * and duplicate members, and flag what an engineer should look at.
 * @param {object} frame frame model
 * @param {{ tolerance?: number, split?: boolean, removeUnusedNodes?: boolean }} [opts] merge / on-member
 *   distance [m] (default 1e-4 × bounding-box diagonal, at least 1 mm); split = false keeps members whole
 * @returns {{ model: object, actions: Array<{ code, count, message }>, flags: { freeEndNodes: number[],
 *   overlappingMembers: number[][], shortMembers: number[], isolatedParts: number } }} `overlappingMembers` are
 *   pairs of input member indices that lie on top of each other over some length (the duplicate part is removed);
 *   `freeEndNodes` and `shortMembers` (shorter than 10 × tolerance) index the healed model.
 */
export function healFrame(frame, { tolerance, split = true, removeUnusedNodes = true } = {}) {
  const bb = bounds(frame.nodes);
  const tol = tolerance ?? Math.max(1e-3, 1e-4 * bb.diagonal);
  if (!(tol >= 0)) throw new RangeError('tolerance: must be >= 0');
  const actions = [];
  const w = weldVertices(frame.nodes, tol);
  if (w.merged) actions.push({ code: 'merged-nodes', count: w.merged, message: `Merged ${w.merged} node(s) lying within ${tol.toPrecision(2)} m of another node.` });
  let N = Array.from(w.positions);
  const nN = w.count;
  let members = frame.members.map((m, i) => ({ ...m, a: w.map[m.a], b: w.map[m.b], _src: i }));
  let nSplit = 0;
  if (split) {
    const out = [];
    for (const m of members) {
      const ax = N[3 * m.a], ay = N[3 * m.a + 1], az = N[3 * m.a + 2];
      const dx = N[3 * m.b] - ax, dy = N[3 * m.b + 1] - ay, dz = N[3 * m.b + 2] - az;
      const L2 = dx * dx + dy * dy + dz * dz;
      if (!(L2 > tol * tol)) { out.push(m); continue; }
      const L = Math.sqrt(L2);
      const on = [];
      for (let i = 0; i < nN; i++) {
        if (i === m.a || i === m.b) continue;
        const px = N[3 * i] - ax, py = N[3 * i + 1] - ay, pz = N[3 * i + 2] - az;
        const t = (px * dx + py * dy + pz * dz) / L2;
        if (t * L <= tol || (1 - t) * L <= tol) continue;
        const qx = px - t * dx, qy = py - t * dy, qz = pz - t * dz;
        if (qx * qx + qy * qy + qz * qz <= tol * tol) on.push([t, i]);
      }
      if (!on.length) { out.push(m); continue; }
      on.sort((p, q) => p[0] - q[0]);
      let prev = m.a;
      for (const [, i] of on) { out.push({ ...m, a: prev, b: i }); prev = i; }
      out.push({ ...m, a: prev, b: m.b });
      nSplit += on.length;
    }
    members = out;
    if (nSplit) actions.push({ code: 'split-members', count: nSplit, message: `Split members at ${nSplit} intermediate node(s) so that crossing members are connected.` });
  }
  const seen = new Map();
  const kept = [];
  const overlaps = [];
  let zero = 0, dup = 0;
  for (const m of members) {
    if (m.a === m.b) { zero++; continue; }
    const key = m.a < m.b ? `${m.a},${m.b}` : `${m.b},${m.a}`;
    const prev = seen.get(key);
    if (prev !== undefined) { dup++; if (prev !== m._src) overlaps.push([prev, m._src]); continue; }
    seen.set(key, m._src);
    kept.push(m);
  }
  if (zero) actions.push({ code: 'removed-zero-length', count: zero, message: `Removed ${zero} zero-length member(s).` });
  if (dup) actions.push({ code: 'removed-duplicates', count: dup, message: `Removed ${dup} duplicate or overlapping member segment(s).` });
  let nodeMap = null;
  let nodes = Float64Array.from(N);
  if (removeUnusedNodes) {
    const used = new Int32Array(nN).fill(-1);
    let n = 0;
    for (const m of kept) { if (used[m.a] < 0) used[m.a] = n++; if (used[m.b] < 0) used[m.b] = n++; }
    if (n < nN) {
      const nn = new Float64Array(3 * n);
      for (let i = 0; i < nN; i++) if (used[i] >= 0) { nn[3 * used[i]] = N[3 * i]; nn[3 * used[i] + 1] = N[3 * i + 1]; nn[3 * used[i] + 2] = N[3 * i + 2]; }
      actions.push({ code: 'removed-unused-nodes', count: nN - n, message: `Removed ${nN - n} node(s) not used by any member.` });
      nodes = nn; nodeMap = used;
      for (const m of kept) { m.a = used[m.a]; m.b = used[m.b]; }
    }
  }
  N = null;
  const origin = kept.map((m) => m._src);
  for (const m of kept) delete m._src;
  const model = { ...frame, nodes, members: kept, memberOrigin: Uint32Array.from(origin) };
  delete model.nodeIds;
  if (frame.sets) {
    const sets = {};
    if (frame.sets.nodes) {
      sets.nodes = {};
      for (const [k, ids] of Object.entries(frame.sets.nodes)) {
        const s = new Set();
        for (const i of ids) { let j = w.map[i]; if (nodeMap) j = nodeMap[j]; if (j >= 0) s.add(j); }
        sets.nodes[k] = Uint32Array.from(s);
      }
    }
    if (frame.sets.members) {
      sets.members = {};
      for (const [k, ids] of Object.entries(frame.sets.members)) {
        const s = new Set(ids);
        sets.members[k] = Uint32Array.from(origin.map((o, i) => (s.has(o) ? i : -1)).filter((i) => i >= 0));
      }
    }
    model.sets = sets;
  }
  const nn = nodes.length / 3;
  const deg = new Uint32Array(nn);
  const p = new Uint32Array(nn);
  for (let i = 0; i < nn; i++) p[i] = i;
  const find = (i) => { while (p[i] !== i) { p[i] = p[p[i]]; i = p[i]; } return i; };
  const short = [];
  kept.forEach((m, i) => {
    deg[m.a]++; deg[m.b]++; p[find(m.a)] = find(m.b);
    if (Math.hypot(nodes[3 * m.a] - nodes[3 * m.b], nodes[3 * m.a + 1] - nodes[3 * m.b + 1], nodes[3 * m.a + 2] - nodes[3 * m.b + 2]) < 10 * tol) short.push(i);
  });
  const free = [];
  const roots = new Set();
  for (let i = 0; i < nn; i++) { if (deg[i] === 1) free.push(i); if (deg[i] > 0) roots.add(find(i)); }
  return { model, actions, flags: { freeEndNodes: free, overlappingMembers: overlaps, shortMembers: short, isolatedParts: roots.size } };
}

function repairSurface(model, opts) {
  const actions = [];
  const act = (code, count, message) => { if (count) actions.push({ code, count, message }); };
  let m = model;
  const bb = bounds(m.positions);
  const tol = opts.tolerance ?? 1e-6 * bb.diagonal;
  if (opts.weld !== false) {
    const w = weldVertices(m.positions, tol);
    if (w.merged) {
      const idx = new Uint32Array(m.indices.length);
      for (let i = 0; i < idx.length; i++) idx[i] = w.map[m.indices[i]];
      m = { ...m, positions: w.positions, indices: idx };
      act('welded-vertices', w.merged, `Welded ${w.merged} coincident vertices (tolerance ${tol.toExponential(1)} m).`);
    }
  }
  let r = removeDegenerateFaces(m); m = r.model; act('removed-degenerate-faces', r.removed, `Removed ${r.removed} zero-area face(s).`);
  r = removeDuplicateFaces(m); m = r.model; act('removed-duplicate-faces', r.removed, `Removed ${r.removed} duplicate face(s).`);
  if (opts.dropSmallComponents !== false) {
    r = dropSmallComponents(m, { minFaces: opts.minComponentFaces, minAreaFraction: opts.minComponentAreaFraction });
    m = r.model; act('dropped-small-components', r.removedComponents, `Dropped ${r.removedComponents} tiny disconnected part(s) (${r.removedFaces} faces).`);
  }
  if (opts.splitNonManifoldVertices !== false) {
    r = splitNonManifoldVertices(m); m = r.model; act('split-non-manifold-vertices', r.split, `Separated ${r.split} bow-tie vertex connection(s).`);
  }
  if (opts.fillHoles !== false) {
    r = fillHoles(m, { maxEdges: opts.maxHoleEdges ?? 32 }); m = r.model;
    act('filled-holes', r.filled, `Filled ${r.filled} hole(s) with ${r.addedFaces} new face(s).`);
    if (r.skipped) actions.push({ code: 'holes-left-open', count: r.skipped, message: `${r.skipped} boundary loop(s) were left open (larger than ${opts.maxHoleEdges ?? 32} edges or not a simple loop).` });
  }
  if (opts.fixOrientation !== false) {
    r = fixOrientation(m); m = r.model; act('flipped-faces', r.flipped, `Reversed ${r.flipped} face(s) to give consistent outward normals.`);
    if (!r.orientable) actions.push({ code: 'not-orientable', count: 1, message: 'A consistent orientation does not exist (non-manifold or one-sided surface); normals were made consistent as far as possible.' });
  }
  r = removeUnreferencedVertices(m); m = r.model; act('removed-unreferenced-vertices', r.removed, `Removed ${r.removed} unused vertices.`);
  const ed = buildEdges(m.indices, m.positions.length / 3);
  let nm = 0;
  for (let e = 0; e < ed.edgeCount; e++) if (ed.edgeFaceStart[e + 1] - ed.edgeFaceStart[e] > 2) nm++;
  if (nm) actions.push({ code: 'non-manifold-edges-remain', count: nm, message: `${nm} edge(s) are still shared by more than two faces; this needs a modelling decision and was not changed.` });
  return { model: m, actions };
}

const TET10_SWAP = [0, 2, 1, 3, 6, 5, 4, 7, 9, 8];

function repairVolume(model) {
  const actions = [];
  const nV = (model.positions.length / 3) | 0;
  let flipped = 0, dropped = 0, unfixable = 0;
  const cells = [];
  const cellMap = [];
  let gi = 0, gn = 0;
  for (const block of model.cells) {
    const nn = CELL_INFO[block.type].nodes;
    const nc = (block.conn.length / nn) | 0;
    const conn = new Uint32Array(block.conn.length);
    const keepIds = [];
    let n = 0;
    for (let c = 0; c < nc; c++, gi++) {
      const o = c * nn;
      const v = cellVolume(model.positions, block.conn, o, block.type);
      const loc = block.conn.subarray ? block.conn.subarray(o, o + nn) : block.conn.slice(o, o + nn);
      if (v === 0 || !Number.isFinite(v)) { dropped++; cellMap.push(-1); continue; }
      const dst = n * nn;
      conn.set(loc, dst);
      if (v < 0) {
        if (block.type === 'tet4') { conn[dst + 1] = loc[2]; conn[dst + 2] = loc[1]; flipped++; }
        else if (block.type === 'tet10') { for (let k = 0; k < 10; k++) conn[dst + k] = loc[TET10_SWAP[k]]; flipped++; }
        else if (block.type === 'hex8') { for (let k = 0; k < 4; k++) { conn[dst + k] = loc[k + 4]; conn[dst + k + 4] = loc[k]; } flipped++; }
        else if (block.type === 'wedge6') { for (let k = 0; k < 3; k++) { conn[dst + k] = loc[k + 3]; conn[dst + k + 3] = loc[k]; } flipped++; }
        else if (block.type === 'pyr5') { conn[dst + 1] = loc[3]; conn[dst + 3] = loc[1]; flipped++; }
        else unfixable++;
      }
      cellMap.push(gn++); keepIds.push(c); n++;
    }
    const nb = { ...block, conn: conn.slice(0, n * nn) };
    for (const k of ['ids', 'tags']) if (block[k]) nb[k] = block[k].constructor.from(keepIds, (c) => block[k][c]);
    if (n) cells.push(nb);
  }
  if (flipped) actions.push({ code: 'reordered-inverted-cells', count: flipped, message: `Reordered the nodes of ${flipped} inverted cell(s) to give positive volume.` });
  if (dropped) actions.push({ code: 'removed-zero-volume-cells', count: dropped, message: `Removed ${dropped} zero-volume cell(s).` });
  if (unfixable) actions.push({ code: 'inverted-cells-remain', count: unfixable, message: `${unfixable} inverted hex20 cell(s) were left as they are (mid-side node ordering is format dependent).` });
  const map = new Int32Array(nV).fill(-1);
  let n = 0;
  for (const b of cells) for (let i = 0; i < b.conn.length; i++) if (map[b.conn[i]] < 0) map[b.conn[i]] = n++;
  let positions = model.positions;
  const out = { ...model, cells };
  if (n < nV) {
    positions = new model.positions.constructor(3 * n);
    for (let v = 0; v < nV; v++) if (map[v] >= 0) { positions[3 * map[v]] = model.positions[3 * v]; positions[3 * map[v] + 1] = model.positions[3 * v + 1]; positions[3 * map[v] + 2] = model.positions[3 * v + 2]; }
    for (const b of cells) for (let i = 0; i < b.conn.length; i++) b.conn[i] = map[b.conn[i]];
    actions.push({ code: 'removed-unreferenced-nodes', count: nV - n, message: `Removed ${nV - n} unused node(s).` });
    out.positions = positions;
    if (model.nodeIds) out.nodeIds = remapPick(model.nodeIds, map, n);
  }
  if (model.sets && (n < nV || dropped)) {
    const remapNodes = (ids) => Uint32Array.from(Array.from(ids, (i) => map[i]).filter((i) => i >= 0));
    const sets = {};
    for (const [k, s] of Object.entries(model.sets)) {
      if (k === 'nodes' && isBag(s)) { sets.nodes = {}; for (const [nm, ids] of Object.entries(s)) sets.nodes[nm] = n < nV ? remapNodes(ids) : ids; }
      else if (k === 'cells' && isBag(s)) { sets.cells = {}; for (const [nm, ids] of Object.entries(s)) sets.cells[nm] = Uint32Array.from(Array.from(ids, (i) => cellMap[i]).filter((i) => i >= 0)); }
      else if (s && s.nodes && n < nV) {
        const fs = { ...s, nodes: remapNodes(s.nodes) };
        for (const key of ['tri3', 'quad4', 'tri6']) if (s[key]) fs[key] = Uint32Array.from(s[key], (i) => map[i]);
        sets[k] = fs;
      } else sets[k] = s;
    }
    out.sets = sets;
  }
  return { model: out, actions };
}

const isBag = (s) => s && typeof s === 'object' && !ArrayBuffer.isView(s) && !Array.isArray(s) && !ArrayBuffer.isView(s.nodes);

function remapPick(arr, map, n) {
  const out = new arr.constructor(n);
  for (let v = 0; v < map.length; v++) if (map[v] >= 0) out[map[v]] = arr[v];
  return out;
}

function repairSimple(model, opts) {
  const actions = [];
  if (model.kind === 'polyline') {
    const P = model.points, n = (P.length / 3) | 0;
    const keep = [];
    for (let i = 0; i < n; i++) {
      const x = P[3 * i], y = P[3 * i + 1], z = P[3 * i + 2];
      if (!Number.isFinite(x + y + z)) continue;
      const l = keep.length ? keep[keep.length - 1] : -1;
      if (l >= 0 && P[3 * l] === x && P[3 * l + 1] === y && P[3 * l + 2] === z) continue;
      keep.push(i);
    }
    if (keep.length === n) return { model, actions };
    const pts = new P.constructor(3 * keep.length);
    keep.forEach((i, j) => { pts[3 * j] = P[3 * i]; pts[3 * j + 1] = P[3 * i + 1]; pts[3 * j + 2] = P[3 * i + 2]; });
    const out = { ...model, points: pts };
    if (model.attrs) { out.attrs = {}; for (const [k, a] of Object.entries(model.attrs)) out.attrs[k] = a && a.length === n ? a.constructor.from(keep, (i) => a[i]) : a; }
    actions.push({ code: 'removed-points', count: n - keep.length, message: `Removed ${n - keep.length} repeated or invalid point(s).` });
    return { model: out, actions };
  }
  if (model.kind === 'points') {
    const P = model.positions, n = (P.length / 3) | 0;
    const keep = [];
    const seen = opts.removeDuplicates === false ? null : new Set();
    for (let i = 0; i < n; i++) {
      const x = P[3 * i], y = P[3 * i + 1], z = P[3 * i + 2];
      if (!Number.isFinite(x + y + z)) continue;
      if (seen) { const k = `${x},${y},${z}`; if (seen.has(k)) continue; seen.add(k); }
      keep.push(i);
    }
    if (keep.length === n) return { model, actions };
    const pts = new P.constructor(3 * keep.length);
    keep.forEach((i, j) => { pts[3 * j] = P[3 * i]; pts[3 * j + 1] = P[3 * i + 1]; pts[3 * j + 2] = P[3 * i + 2]; });
    const out = { ...model, positions: pts };
    if (model.attrs) { out.attrs = {}; for (const [k, a] of Object.entries(model.attrs)) out.attrs[k] = a && a.length === n ? a.constructor.from(keep, (i) => a[i]) : a; }
    actions.push({ code: 'removed-points', count: n - keep.length, message: `Removed ${n - keep.length} duplicate or invalid point(s).` });
    return { model: out, actions };
  }
  if (model.kind === 'grid') {
    // fill isolated no-data nodes from the mean of valid 8-neighbours, a few passes for small gaps
    const { nx, ny } = model;
    let z = Float32Array.from(model.z);
    const passes = opts.maxGap ?? 2;
    let filled = 0;
    for (let p = 0; p < passes; p++) {
      const next = z.slice();
      let did = 0;
      for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
        if (Number.isFinite(z[j * nx + i])) continue;
        let s = 0, c = 0;
        for (let dj = -1; dj <= 1; dj++) for (let di = -1; di <= 1; di++) {
          const a = i + di, b = j + dj;
          if (a < 0 || b < 0 || a >= nx || b >= ny) continue;
          const v = z[b * nx + a];
          if (Number.isFinite(v)) { s += v; c++; }
        }
        if (c >= 3) { next[j * nx + i] = s / c; did++; }
      }
      z = next; filled += did;
      if (!did) break;
    }
    if (!filled) return { model, actions };
    actions.push({ code: 'filled-no-data', count: filled, message: `Filled ${filled} no-data node(s) from neighbouring values (gaps up to ${passes} node(s) wide from each side).` });
    return { model: { ...model, z }, actions };
  }
  if (model.kind === 'section') {
    let fixed = 0, flipped = 0;
    const loops = model.loops.map((L, li) => {
      const pts = [];
      const n = (L.length / 2) | 0;
      for (let i = 0; i < n; i++) {
        const x = L[2 * i], y = L[2 * i + 1];
        if (!Number.isFinite(x + y)) { fixed++; continue; }
        if (pts.length && pts[pts.length - 2] === x && pts[pts.length - 1] === y) { fixed++; continue; }
        pts.push(x, y);
      }
      if (pts.length >= 4 && pts[0] === pts[pts.length - 2] && pts[1] === pts[pts.length - 1]) { pts.length -= 2; fixed++; }
      let a = 0;
      const m = pts.length / 2;
      for (let i = 0; i < m; i++) { const j = (i + 1) % m; a += pts[2 * i] * pts[2 * j + 1] - pts[2 * j] * pts[2 * i + 1]; }
      // outer boundary counter-clockwise, further loops (holes) clockwise
      if ((li === 0 && a < 0) || (li > 0 && a > 0)) {
        const rev = [];
        for (let i = m - 1; i >= 0; i--) rev.push(pts[2 * i], pts[2 * i + 1]);
        flipped++;
        return Float64Array.from(rev);
      }
      return Float64Array.from(pts);
    }).filter((L) => L.length >= 6);
    if (fixed) actions.push({ code: 'cleaned-loop-points', count: fixed, message: `Removed ${fixed} repeated, closing or invalid loop point(s).` });
    if (flipped) actions.push({ code: 'reversed-loops', count: flipped, message: `Reversed ${flipped} loop(s): outer boundary counter-clockwise, holes clockwise.` });
    if (loops.length !== model.loops.length) actions.push({ code: 'removed-degenerate-loops', count: model.loops.length - loops.length, message: `Removed ${model.loops.length - loops.length} loop(s) with fewer than three points.` });
    return { model: actions.length ? { ...model, loops } : model, actions };
  }
  throw new RangeError(`model.kind: '${model.kind}' cannot be repaired`);
}

/**
 * Repair a model of any kind and report each action taken.
 *
 * surface: weld → remove degenerate → remove duplicate → drop tiny components → split bow-tie vertices →
 *          fill small holes → consistent outward orientation → drop unused vertices
 * volume:  drop zero-volume cells, reorder inverted cells, drop unused nodes (sets remapped)
 * frame:   healFrame
 * polyline / points: drop invalid and repeated points; grid: fill small no-data gaps; section: clean and orient loops
 *
 * @param {object} model CONTRACT §4 model
 * @param {object} [opts] surface: { tolerance, weld, fillHoles, maxHoleEdges = 32, fixOrientation,
 *   splitNonManifoldVertices, dropSmallComponents, minComponentFaces, minComponentAreaFraction } (booleans default
 *   true); frame: healFrame options; grid: { maxGap = 2 }; points: { removeDuplicates = true }
 * @returns {{ model: object, actions: Array<{ code: string, count: number, message: string }>, flags?: object }}
 */
export function repair(model, opts = {}) {
  if (!model || typeof model !== 'object') throw new RangeError('model: a geometry model is required');
  if (model.kind === 'surface') return repairSurface(model, opts);
  if (model.kind === 'volume') return repairVolume(model, opts);
  if (model.kind === 'frame') return healFrame(model, opts);
  return repairSimple(model, opts);
}
