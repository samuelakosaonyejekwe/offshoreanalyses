// Hand-off to the finite-element solver CalculiX (ccx) and back.
//
//   writeCalculixCase   frame / surface / volume model → Abaqus-style input deck (.inp) with material, supports,
//                       loads and one analysis step
//   readCalculixFrd     ASCII results file (.frd): nodal DISP, STRESS … blocks → arrays in the app's node order
//   readCalculixDat     printed tables (.dat): displacements, stresses at integration points, reaction forces
//
// Target: CalculiX CrunchiX 2.20 / 2.21 (Dhondt, "CalculiX CrunchiX User's Manual version 2.21", 2023). The deck
// keeps to keywords shared with Abaqus/Standard, so the app's own Abaqus importer reads the mesh back.
//
// Numbering: node id = app node index + 1; element id = member / triangle / cell index + 1 (cells counted
// through the blocks of a volume model in order; a 5-node pyramid, which CalculiX does not have, becomes two
// tetrahedra and shifts the following ids). Units are SI throughout (m, N, Pa, kg/m³).
//
// Element choice: frame → B31 (2-node) or B32 (3-node, mid-side nodes generated) beams; surface → S3 shells;
// volume → C3D4, C3D10, C3D8, C3D20, C3D6. CalculiX expands beams and shells into solid elements internally;
// the output requests use OUTPUT=2D so that results come back on the original nodes.

import { worldCoords } from '../geometry/exporters.js';

const num = (v) => { if (!Number.isFinite(v)) throw new RangeError('calculix: a non-finite number reached the deck'); return String(Number(v.toPrecision(12))); };
const ANALYSES = ['static', 'static-nlgeom', 'buckle', 'frequency', 'dynamic'];
const CELL = { tet4: ['C3D4', 4], tet10: ['C3D10', 10], hex8: ['C3D8', 8], hex20: ['C3D20', 20], wedge6: ['C3D6', 6] };

function idLines(ids, perLine = 12) {
  const out = [];
  for (let i = 0; i < ids.length; i += perLine) out.push(ids.slice(i, i + perLine).join(', '));
  return out;
}

function resolveNodes(spec, model, nNodes, what) {
  let list = spec;
  if (typeof spec === 'string') {
    const sets = (model.sets && model.sets.nodes) || {};
    const key = Object.keys(sets).find((k) => k.toUpperCase() === spec.toUpperCase());
    if (!key) throw new RangeError(`${what}: node set "${spec}" is not in the model`);
    list = Array.from(sets[key]);
  }
  if (!list || !list.length) throw new RangeError(`${what}: nodes must be a non-empty list of node indices or the name of a node set`);
  const out = Array.from(list);
  for (const i of out) if (!Number.isInteger(i) || i < 0 || i >= nNodes) throw new RangeError(`${what}: node index ${i} is outside 0…${nNodes - 1}`);
  return out;
}

/** Rectangle a × b with the same second moments as a generic section (Iy about local 1, Iz about local 2). */
function equivalentRect(sec) {
  const k = Math.sqrt(sec.Iy / sec.Iz);
  const a = (12 * sec.Iy / k ** 3) ** 0.25;
  return { a, b: k * a };
}

function sectionCard(sec, mode, notes, label) {
  if (!sec) throw new RangeError(`calculix: member group "${label}" has no section`);
  if (sec.type === 'tube') {
    if (!(sec.D > 0) || !(sec.t > 0) || !(sec.t < sec.D / 2)) throw new RangeError(`calculix: tube section of "${label}" needs D > 2 t > 0`);
    return { kw: '*BEAM SECTION', section: 'PIPE', data: `${num(sec.D / 2)}, ${num(sec.t)}` };
  }
  if (sec.shape === 'circle' && sec.D > 0) return { kw: '*BEAM SECTION', section: 'CIRC', data: `${num(sec.D / 2)}, ${num(sec.D / 2)}` };
  if (sec.shape === 'rect' && sec.b > 0 && sec.h > 0) return { kw: '*BEAM SECTION', section: 'RECT', data: `${num(sec.b)}, ${num(sec.h)}` };
  if (sec.shape === 'box' && sec.b > 0 && sec.h > 0 && sec.tw > 0 && sec.tf > 0) return { kw: '*BEAM SECTION', section: 'BOX', data: `${num(sec.b)}, ${num(sec.h)}, ${num(sec.tw)}, ${num(sec.tf)}, ${num(sec.tw)}, ${num(sec.tf)}` };
  if (!(sec.A > 0) || !(sec.Iy > 0) || !(sec.Iz > 0)) throw new RangeError(`calculix: section of "${label}" needs { type: 'tube', D, t } or { A, Iy, Iz, J }`);
  if (mode === 'general') {
    notes.push(`Section "${label}" is written as *BEAM GENERAL SECTION (A, I11, I12, I22, J): Abaqus reads it, CalculiX does not — choose genericSection: 'rect-equivalent' to run the deck in ccx.`);
    return { kw: '*BEAM GENERAL SECTION', section: 'GENERAL', data: `${num(sec.A)}, ${num(sec.Iy)}, 0., ${num(sec.Iz)}, ${num(sec.J > 0 ? sec.J : sec.Iy + sec.Iz)}` };
  }
  const r = equivalentRect(sec);
  notes.push(`Section "${label}" (A = ${num(sec.A)} m², Iy = ${num(sec.Iy)}, Iz = ${num(sec.Iz)} m⁴) is written as the rectangle ${num(Number(r.a.toPrecision(6)))} × ${num(Number(r.b.toPrecision(6)))} m with the same Iy and Iz: its area is ${num(Number((r.a * r.b / sec.A).toPrecision(4)))} × A, so axial stiffness, shear stiffness, torsion and mass differ — CalculiX has no general beam section.`);
  return { kw: '*BEAM SECTION', section: 'RECT', data: `${num(r.a)}, ${num(r.b)}` };
}

/** Direction of the local 1-axis of a beam: the member's own `orient`, else (0, 0, -1), else (0, 1, 0) for vertical members. */
function beamOrient(m, nodes) {
  if (m.orient && m.orient.length >= 3) return m.orient.slice(0, 3);
  const dx = nodes[3 * m.b] - nodes[3 * m.a]; const dy = nodes[3 * m.b + 1] - nodes[3 * m.a + 1]; const dz = nodes[3 * m.b + 2] - nodes[3 * m.a + 2];
  const l = Math.hypot(dx, dy, dz);
  if (!(l > 0)) throw new RangeError('calculix: a member has zero length');
  return Math.abs(dz) / l > 0.99 ? [0, 1, 0] : [0, 0, -1];
}

/**
 * Writes a CalculiX input deck for a `frame`, `surface` or `volume` model (CONTRACT §4).
 *
 * Keywords written: *HEADING, *NODE, *ELEMENT, *NSET, *ELSET, *MATERIAL with *ELASTIC, *DENSITY and (for
 * static, static-nlgeom and dynamic analyses with a yield stress) *PLASTIC, *BEAM SECTION / *SHELL SECTION /
 * *SOLID SECTION, *BOUNDARY, one *STEP with *STATIC, *BUCKLE, *FREQUENCY or *DYNAMIC, *CLOAD, *DLOAD (gravity,
 * pressure), *NODE FILE, *EL FILE, *NODE PRINT, *EL PRINT.
 * Reference: CalculiX CrunchiX User's Manual v2.21, ch. 6 (element types) and ch. 7 (input deck format).
 * Validity: small-strain metal plasticity with isotropic hardening; one material for the whole model; loads
 * applied in a single step (ramped over the step for static analyses, applied as a step for *DYNAMIC unless an
 * amplitude is added by hand).
 *
 * @param {object} model frame, surface or volume model, coordinates in m
 * @param {object} o
 * @param {{E:number, nu:number, rho?:number, fy?:number, hardening?:{Et?:number, epsMax?:number,
 *   table?:number[][]}, name?:string}} o.material Young's modulus Pa, Poisson's ratio, density kg/m³, yield
 *   stress Pa; hardening: tangent modulus Et (Pa, bilinear, written up to plastic strain epsMax, default 0.2)
 *   or a table of [true stress Pa, plastic strain] pairs starting at zero plastic strain
 * @param {{nodes:number[]|string, dofs?:'fixed'|'pinned'|number[], value?:number}[]} o.supports node indices
 *   (0-based) or the name of a node set of the model; dofs 'fixed' (all), 'pinned' (translations) or a list of
 *   degrees of freedom 1…6; value: prescribed displacement, m or rad (default 0)
 * @param {object[]} [o.loads] each one of
 *   `{ type: 'force', nodes, F: [Fx, Fy, Fz] (N), M?: [Mx, My, Mz] (N·m), share?: boolean }` — the load is applied
 *   to every listed node, or divided equally among them when `share` is true;
 *   `{ type: 'gravity', g?: 9.80665 (m/s²), direction?: [0, 0, -1] }`;
 *   `{ type: 'pressure', p (Pa), elements?: number[] (0-based; default all), face?: 1…6 (volume models) }` —
 *   on shells a positive pressure acts against the element normal;
 *   `{ type: 'line', members: number[], q: [qx, qy, qz] (N/m) }` — frame only, lumped to the member end nodes
 *   (q L / 2 each; the fixed-end moments of a consistent load vector are not applied)
 * @param {'static'|'static-nlgeom'|'buckle'|'frequency'|'dynamic'} [o.analysis='static']
 * @param {number|{nModes?:number, increments?:number, initialIncrement?:number, minIncrement?:number,
 *   maxIncrement?:number, timeStep?:number, endTime?:number, nlgeom?:boolean}} [o.steps] a number sets the
 *   mode count (buckle, frequency) or the increment limit (other analyses); time values in s (dynamic) or as
 *   fractions of the step (static)
 * @param {1|2} [o.elementOrder=1] frame: 1 → B31, 2 → B32 (mid-side nodes are generated)
 * @param {number} [o.thickness] shell thickness, m (surface models; a group's own `thickness` wins)
 * @param {'rect-equivalent'|'general'} [o.genericSection='rect-equivalent'] how { A, Iy, Iz, J } sections are written
 * @param {string} [o.name] job name (file name without extension)
 * @returns {{files:{path:string, text:string}[], notes:string[], counts:{nodes:number, elements:number,
 *   appNodes:number}}} the deck `<name>.inp`; counts.nodes includes generated mid-side nodes
 */
export function writeCalculixCase(model, { material, supports, loads = [], analysis = 'static', steps = {}, elementOrder = 1, thickness, genericSection = 'rect-equivalent', name } = {}) {
  if (!model || !['frame', 'surface', 'volume'].includes(model.kind)) throw new RangeError("model must be a 'frame', 'surface' or 'volume' model");
  if (!ANALYSES.includes(analysis)) throw new RangeError(`analysis must be one of ${ANALYSES.join(', ')}`);
  if (!material || !(material.E > 0)) throw new RangeError('material.E must be > 0');
  if (!(material.nu > -1 && material.nu < 0.5)) throw new RangeError('material.nu must lie in (-1, 0.5)');
  const needsMass = analysis === 'frequency' || analysis === 'dynamic' || loads.some((l) => l && l.type === 'gravity');
  if (needsMass && !(material.rho > 0)) throw new RangeError('material.rho must be > 0 for frequency, dynamic and gravity analyses');
  if (material.fy !== undefined && !(material.fy > 0)) throw new RangeError('material.fy must be > 0');
  if (![1, 2].includes(elementOrder)) throw new RangeError('elementOrder must be 1 or 2');
  if (!['rect-equivalent', 'general'].includes(genericSection)) throw new RangeError("genericSection must be 'rect-equivalent' or 'general'");
  if (!Array.isArray(supports)) throw new RangeError('supports must be a list');
  if (supports.length === 0 && analysis !== 'frequency') throw new RangeError('supports must not be empty (the model would be free to move as a rigid body)');
  const st = typeof steps === 'number' ? (analysis === 'buckle' || analysis === 'frequency' ? { nModes: steps } : { increments: steps }) : (steps || {});
  const job = String(name || model.name || 'model').replace(/[^A-Za-z0-9_-]+/g, '_').replace(/^_+|_+$/g, '') || 'model';
  const mat = String(material.name || 'MAT').replace(/[^A-Za-z0-9_-]+/g, '_').toUpperCase();
  const notes = [];
  const X = worldCoords(model);
  const nApp = X.length / 3;
  const L = ['*HEADING', `${job} - ${model.kind} model, ${analysis} analysis, SI units (m, N, Pa, kg)`, '*NODE, NSET=NALL'];
  for (let i = 0; i < nApp; i++) L.push(`${i + 1}, ${num(X[3 * i])}, ${num(X[3 * i + 1])}, ${num(X[3 * i + 2])}`);
  let nNodes = nApp; let nElem = 0;
  const sections = [];
  const rotational = model.kind !== 'volume';
  // element id(s) of each app element, for load sets
  const elemIds = [];

  if (model.kind === 'frame') {
    const mem = model.members || [];
    if (!mem.length) throw new RangeError('model.members is empty');
    const groups = new Map();
    mem.forEach((m, i) => {
      if (!Number.isInteger(m.a) || !Number.isInteger(m.b) || m.a < 0 || m.b < 0 || m.a >= nApp || m.b >= nApp) throw new RangeError(`member ${i} refers to a node outside the model`);
      const o = beamOrient(m, X);
      const key = JSON.stringify([m.section || null, o.map((v) => Number(v.toPrecision(9)))]);
      if (!groups.has(key)) groups.set(key, { section: m.section, orient: o, label: m.group || `section ${groups.size + 1}`, ids: [] });
      groups.get(key).ids.push(i + 1);
      elemIds.push([i + 1]);
    });
    if (elementOrder === 2) {
      for (const m of mem) L.push(`${++nNodes}, ${num((X[3 * m.a] + X[3 * m.b]) / 2)}, ${num((X[3 * m.a + 1] + X[3 * m.b + 1]) / 2)}, ${num((X[3 * m.a + 2] + X[3 * m.b + 2]) / 2)}`);
      notes.push(`${mem.length} mid-side nodes (ids ${nApp + 1}…${nNodes}) were generated for the B32 elements; results on them are ignored by the readers.`);
    }
    L.push(`*ELEMENT, TYPE=${elementOrder === 2 ? 'B32' : 'B31'}, ELSET=EALL`);
    mem.forEach((m, i) => L.push(elementOrder === 2 ? `${i + 1}, ${m.a + 1}, ${nApp + 1 + i}, ${m.b + 1}` : `${i + 1}, ${m.a + 1}, ${m.b + 1}`));
    nElem = mem.length;
    let g = 0;
    for (const grp of groups.values()) {
      const set = `BEAM${++g}`;
      L.push(`*ELSET, ELSET=${set}`, ...idLines(grp.ids));
      const card = sectionCard(grp.section, genericSection, notes, grp.label);
      sections.push(`${card.kw}, ELSET=${set}, MATERIAL=${mat}, SECTION=${card.section}`, card.data, grp.orient.map(num).join(', '));
    }
    notes.push('Beam local 1-axis: the member\'s own orientation when the model carries one, otherwise (0, 0, -1), or (0, 1, 0) for members within 8° of vertical.');
  } else if (model.kind === 'surface') {
    const idx = model.indices;
    if (!idx || idx.length < 3) throw new RangeError('model.indices is empty');
    L.push('*ELEMENT, TYPE=S3, ELSET=EALL');
    for (let f = 0; f < idx.length / 3; f++) { L.push(`${f + 1}, ${idx[3 * f] + 1}, ${idx[3 * f + 1] + 1}, ${idx[3 * f + 2] + 1}`); elemIds.push([f + 1]); }
    nElem = idx.length / 3;
    const groups = model.groups && model.groups.length ? model.groups : [{ name: 'all', start: 0, count: idx.length }];
    groups.forEach((grp, k) => {
      const t = grp.thickness > 0 ? grp.thickness : thickness;
      if (!(t > 0)) throw new RangeError(`thickness must be > 0 (group "${grp.name}" carries none)`);
      const set = `SHELL${k + 1}`;
      const ids = [];
      for (let f = grp.start / 3; f < (grp.start + grp.count) / 3; f++) ids.push(f + 1);
      L.push(`*ELSET, ELSET=${set}`, ...idLines(ids));
      sections.push(`*SHELL SECTION, ELSET=${set}, MATERIAL=${mat}`, num(t));
    });
  } else {
    const cells = model.cells || [];
    if (!cells.length) throw new RangeError('model.cells is empty');
    let pyr = 0;
    for (const blk of cells) {
      const def = blk.type === 'pyr5' ? ['C3D4', 5] : CELL[blk.type];
      if (!def) throw new RangeError(`cell type "${blk.type}" cannot be written for CalculiX`);
      const n = def[1]; const count = blk.conn.length / n;
      L.push(`*ELEMENT, TYPE=${def[0]}, ELSET=EALL`);
      for (let e = 0; e < count; e++) {
        const c = Array.from(blk.conn.subarray ? blk.conn.subarray(e * n, (e + 1) * n) : blk.conn.slice(e * n, (e + 1) * n), (v) => v + 1);
        if (blk.type === 'pyr5') {
          pyr++;
          L.push(`${++nElem}, ${c[0]}, ${c[1]}, ${c[2]}, ${c[4]}`, `${++nElem}, ${c[0]}, ${c[2]}, ${c[3]}, ${c[4]}`);
          elemIds.push([nElem - 1, nElem]);
        } else {
          // long connectivities wrap: a trailing comma continues the record
          const first = [++nElem, ...c];
          for (let q = 0; q < first.length; q += 11) L.push(first.slice(q, q + 11).join(', ') + (q + 11 < first.length ? ',' : ''));
          elemIds.push([nElem]);
        }
      }
    }
    if (pyr) notes.push(`${pyr} pyramid cell(s) were split into two C3D4 tetrahedra each (CalculiX has no pyramid element); element ids after them are shifted.`);
    sections.push(`*SOLID SECTION, ELSET=EALL, MATERIAL=${mat}`);
  }

  // supports
  const bc = [];
  supports.forEach((s, k) => {
    const nodes = resolveNodes(s.nodes, model, nApp, `supports[${k}]`);
    const set = `SUP${k + 1}`;
    L.push(`*NSET, NSET=${set}`, ...idLines(nodes.map((i) => i + 1)));
    let dofs = s.dofs === undefined || s.dofs === 'fixed' ? (rotational ? [1, 2, 3, 4, 5, 6] : [1, 2, 3]) : s.dofs === 'pinned' ? [1, 2, 3] : s.dofs;
    if (!Array.isArray(dofs) || !dofs.length || dofs.some((d) => !Number.isInteger(d) || d < 1 || d > 6)) throw new RangeError(`supports[${k}].dofs must be 'fixed', 'pinned' or a list of degrees of freedom 1…6`);
    if (!rotational && dofs.some((d) => d > 3)) { dofs = dofs.filter((d) => d <= 3); notes.push(`supports[${k}]: rotational degrees of freedom do not exist on solid elements and were dropped.`); }
    const val = s.value === undefined ? 0 : s.value;
    if (val !== 0 && (analysis === 'buckle' || analysis === 'frequency')) throw new RangeError(`supports[${k}].value: prescribed displacements are not meaningful in a ${analysis} step`);
    for (const d of dofs) bc.push({ set, dof: d, value: val });
  });
  if (supports.length) L.push('*NSET, NSET=NSUP', supports.map((_, k) => `SUP${k + 1}`).join(', '));

  // loads → sets and step cards
  const cload = []; const dload = [];
  const lumped = new Map();
  const addNodal = (node, dof, v) => { if (v === 0) return; const key = `${node},${dof}`; lumped.set(key, (lumped.get(key) || 0) + v); };
  loads.forEach((ld, k) => {
    if (!ld || !ld.type) throw new RangeError(`loads[${k}].type is required`);
    if (ld.type === 'force') {
      const nodes = resolveNodes(ld.nodes, model, nApp, `loads[${k}]`);
      const F = ld.F || [0, 0, 0]; const M = ld.M || [0, 0, 0];
      if (F.length !== 3 || M.length !== 3 || [...F, ...M].some((v) => !Number.isFinite(v))) throw new RangeError(`loads[${k}]: F and M must be [x, y, z]`);
      if (!rotational && M.some((v) => v !== 0)) throw new RangeError(`loads[${k}].M: moments cannot be applied to the nodes of solid elements`);
      const f = ld.share ? 1 / nodes.length : 1;
      if (nodes.length === 1 || ld.share) { for (const i of nodes) { for (let d = 0; d < 3; d++) { addNodal(i + 1, d + 1, f * F[d]); addNodal(i + 1, d + 4, f * M[d]); } } } else {
        const set = `LOAD${k + 1}`;
        L.push(`*NSET, NSET=${set}`, ...idLines(nodes.map((i) => i + 1)));
        for (let d = 0; d < 3; d++) { if (F[d] !== 0) cload.push(`${set}, ${d + 1}, ${num(F[d])}`); if (M[d] !== 0) cload.push(`${set}, ${d + 4}, ${num(M[d])}`); }
      }
    } else if (ld.type === 'gravity') {
      const g = ld.g === undefined ? 9.80665 : ld.g; const d = ld.direction || [0, 0, -1];
      const l = Math.hypot(d[0], d[1], d[2]);
      if (!(g > 0) || !(l > 0)) throw new RangeError(`loads[${k}]: gravity needs g > 0 and a non-zero direction`);
      dload.push(`EALL, GRAV, ${num(g)}, ${num(d[0] / l)}, ${num(d[1] / l)}, ${num(d[2] / l)}`);
    } else if (ld.type === 'pressure') {
      if (model.kind === 'frame') throw new RangeError(`loads[${k}]: pressure applies to surface and volume models (use a 'line' load on frames)`);
      if (!Number.isFinite(ld.p)) throw new RangeError(`loads[${k}].p must be a pressure in Pa`);
      let set = 'EALL';
      if (ld.elements) {
        set = `PRES${k + 1}`;
        const ids = [];
        for (const e of ld.elements) { if (!elemIds[e]) throw new RangeError(`loads[${k}].elements: index ${e} is outside the model`); ids.push(...elemIds[e]); }
        L.push(`*ELSET, ELSET=${set}`, ...idLines(ids));
      }
      if (model.kind === 'volume') {
        if (!Number.isInteger(ld.face) || ld.face < 1 || ld.face > 6) throw new RangeError(`loads[${k}].face must be the element face number 1…6 for a volume model`);
        dload.push(`${set}, P${ld.face}, ${num(ld.p)}`);
      } else dload.push(`${set}, P, ${num(ld.p)}`);
    } else if (ld.type === 'line') {
      if (model.kind !== 'frame') throw new RangeError(`loads[${k}]: line loads apply to frame models`);
      const q = ld.q || [];
      if (q.length !== 3 || q.some((v) => !Number.isFinite(v))) throw new RangeError(`loads[${k}].q must be [qx, qy, qz] in N/m`);
      const list = ld.members || model.members.map((_, i) => i);
      for (const i of list) {
        const m = model.members[i];
        if (!m) throw new RangeError(`loads[${k}].members: index ${i} is outside the model`);
        const len = Math.hypot(X[3 * m.b] - X[3 * m.a], X[3 * m.b + 1] - X[3 * m.a + 1], X[3 * m.b + 2] - X[3 * m.a + 2]);
        for (let d = 0; d < 3; d++) { addNodal(m.a + 1, d + 1, q[d] * len / 2); addNodal(m.b + 1, d + 1, q[d] * len / 2); }
      }
      notes.push(`loads[${k}]: the line load is lumped to the member end nodes (q L / 2 each); refine the members where local bending between nodes matters.`);
    } else throw new RangeError(`loads[${k}].type "${ld.type}" is not one of force, gravity, pressure, line`);
  });
  const nodal = Array.from(lumped, ([key, v]) => { const [node, dof] = key.split(',').map(Number); return { node, dof, v }; }).sort((p, q) => p.node - q.node || p.dof - q.dof);
  for (const c of nodal) cload.push(`${c.node}, ${c.dof}, ${num(c.v)}`);
  if ((analysis === 'static' || analysis === 'static-nlgeom' || analysis === 'buckle') && !cload.length && !dload.length && !bc.some((b) => b.value !== 0)) throw new RangeError(`loads must not be empty for a ${analysis} analysis`);

  // material
  L.push(`*MATERIAL, NAME=${mat}`, '*ELASTIC', `${num(material.E)}, ${num(material.nu)}`);
  if (material.rho > 0) L.push('*DENSITY', num(material.rho));
  const plastic = material.fy > 0 && (analysis === 'static' || analysis === 'static-nlgeom' || analysis === 'dynamic');
  if (plastic) {
    L.push('*PLASTIC');
    const h = material.hardening || {};
    if (h.table) {
      if (!h.table.length || h.table[0][1] !== 0) throw new RangeError('material.hardening.table must start at zero plastic strain');
      for (let i = 0; i < h.table.length; i++) {
        if (i > 0 && !(h.table[i][1] > h.table[i - 1][1])) throw new RangeError('material.hardening.table: plastic strain must increase');
        L.push(`${num(h.table[i][0])}, ${num(h.table[i][1])}`);
      }
    } else if (h.Et !== undefined) {
      if (!(h.Et >= 0 && h.Et < material.E)) throw new RangeError('material.hardening.Et must lie in [0, E)');
      // bilinear: plastic modulus H = E Et / (E - Et), stress against plastic strain
      const H = material.E * h.Et / (material.E - h.Et); const eps = h.epsMax > 0 ? h.epsMax : 0.2;
      L.push(`${num(material.fy)}, 0.`, `${num(material.fy + H * eps)}, ${num(eps)}`);
    } else L.push(`${num(material.fy)}, 0.`);
  } else if (material.fy > 0) notes.push(`The yield stress is not written: a ${analysis} step is a linear perturbation.`);
  L.push(...sections);

  if (bc.length) { L.push('*BOUNDARY'); for (const b of bc) L.push(b.value === 0 ? `${b.set}, ${b.dof}, ${b.dof}` : `${b.set}, ${b.dof}, ${b.dof}, ${num(b.value)}`); }

  // step
  const inc = Number.isInteger(st.increments) && st.increments > 0 ? st.increments : 200;
  const nl = analysis === 'static-nlgeom' || (analysis === 'dynamic' && st.nlgeom === true);
  if (analysis === 'static') L.push(`*STEP${plastic ? `, INC=${inc}` : ''}`, '*STATIC');
  else if (analysis === 'static-nlgeom') {
    const d0 = st.initialIncrement > 0 ? st.initialIncrement : 0.1; const dmin = st.minIncrement > 0 ? st.minIncrement : 1e-6; const dmax = st.maxIncrement > 0 ? st.maxIncrement : Math.max(d0, 0.25);
    L.push(`*STEP, NLGEOM, INC=${inc}`, '*STATIC', `${num(d0)}, 1., ${num(dmin)}, ${num(dmax)}`);
  } else if (analysis === 'buckle') {
    const n = Number.isInteger(st.nModes) && st.nModes > 0 ? st.nModes : 5;
    L.push('*STEP', '*BUCKLE', `${n}`);
    notes.push('The buckling factors multiply the loads of this step; a factor below 1 means the applied load exceeds the elastic critical load.');
  } else if (analysis === 'frequency') {
    const n = Number.isInteger(st.nModes) && st.nModes > 0 ? st.nModes : 10;
    L.push('*STEP', '*FREQUENCY', `${n}`);
    if (cload.length || dload.length) notes.push('Loads are not written in a frequency step.');
  } else {
    const dt = st.timeStep; const tEnd = st.endTime;
    if (!(dt > 0) || !(tEnd > dt)) throw new RangeError('steps.timeStep and steps.endTime (s, endTime > timeStep > 0) are required for a dynamic analysis');
    const nInc = Number.isInteger(st.increments) && st.increments > 0 ? st.increments : Math.ceil(tEnd / dt) + 10;
    L.push(`*STEP${nl ? ', NLGEOM' : ''}, INC=${nInc}`, '*DYNAMIC, DIRECT', `${num(dt)}, ${num(tEnd)}`);
    notes.push('Dynamic step: implicit direct integration with a fixed time step; the loads are applied at full value from t = 0 (add an *AMPLITUDE by hand for a time history).');
  }
  if (analysis !== 'frequency') {
    if (cload.length) L.push('*CLOAD', ...cload);
    if (dload.length) L.push('*DLOAD', ...dload);
  }
  const twoD = model.kind === 'volume' ? '' : ', OUTPUT=2D';
  L.push(`*NODE FILE${twoD}`, analysis === 'static' || analysis === 'static-nlgeom' || analysis === 'dynamic' ? 'U, RF' : 'U');
  L.push(`*EL FILE${twoD}`, plastic ? 'S, E, PEEQ' : 'S, E');
  if (analysis === 'static' || analysis === 'static-nlgeom' || analysis === 'dynamic') {
    L.push('*NODE PRINT, NSET=NALL', 'U');
    if (supports.length) L.push('*NODE PRINT, NSET=NSUP, TOTALS=YES', 'RF');
    L.push('*EL PRINT, ELSET=EALL', 'S');
  }
  L.push('*END STEP', '');
  if (model.kind !== 'volume') notes.push('Beams and shells are expanded to solid elements inside CalculiX; OUTPUT=2D returns the results on the original nodes. Stresses of expanded beams are section-point values, not member forces.');
  notes.push(`Run:  ccx ${job}   (reads ${job}.inp, writes ${job}.frd and ${job}.dat).`);
  return { files: [{ path: `${job}.inp`, text: L.join('\n') }], notes, counts: { nodes: nNodes, elements: nElem, appNodes: nApp } };
}

function flatPositions(positions) {
  if (!positions) return null;
  if (positions.kind) return worldCoords(positions);
  return positions;
}

function vonMises(s, o) {
  const xx = s[o]; const yy = s[o + 1]; const zz = s[o + 2]; const xy = s[o + 3]; const yz = s[o + 4]; const zx = s[o + 5];
  return Math.sqrt(0.5 * ((xx - yy) ** 2 + (yy - zz) ** 2 + (zz - xx) ** 2) + 3 * (xy * xy + yz * yz + zx * zx));
}

function viewerFields(out, pos, n, scale) {
  if (out.stress) {
    out.vonMises = new Float64Array(n);
    for (let i = 0; i < n; i++) out.vonMises[i] = vonMises(out.stress, 6 * i);
    out.nodeValues = Array.from(out.vonMises); out.valueLabel = 'von Mises stress'; out.valueUnit = 'Pa';
  } else if (out.displacement) {
    out.nodeValues = [];
    for (let i = 0; i < n; i++) out.nodeValues.push(Math.hypot(out.displacement[3 * i], out.displacement[3 * i + 1], out.displacement[3 * i + 2]));
    out.valueLabel = 'Displacement magnitude'; out.valueUnit = 'm';
  } else { out.nodeValues = null; out.valueLabel = null; out.valueUnit = null; }
  out.deformed = null;
  if (pos && out.displacement) {
    if (pos.length < 3 * n) throw new RangeError('positions must hold three coordinates per node');
    out.deformed = [];
    for (let i = 0; i < n; i++) out.deformed.push([pos[3 * i] + scale * out.displacement[3 * i], pos[3 * i + 1] + scale * out.displacement[3 * i + 1], pos[3 * i + 2] + scale * out.displacement[3 * i + 2]]);
  }
}

/**
 * Reads an ASCII CalculiX results file (.frd): every nodal result block (DISP, STRESS, TOSTRAIN, FORC, PE …)
 * of every step, mapped to the app's node numbering (node id − 1). Nodes with ids above `nodeCount`
 * (generated mid-side nodes, nodes of expanded elements) are ignored.
 * Reference: CalculiX cgx manual, "Result file format (frd)": record keys 1C/2C/3C, 100C, -4, -5, -1, -2, -3;
 * fixed-width data records (1X, I2, I10, 6E12.5).
 * Validity: ASCII .frd (binary .frd written with *NODE FILE, OUTPUT=3D … and "frdbin" is not read).
 * @param {string} text file content
 * @param {object} [o]
 * @param {number} [o.nodeCount] number of app nodes (default: largest node id in the file)
 * @param {ArrayLike<number>|object} [o.positions] undeformed coordinates (flat xyz, m) or the model itself,
 *   needed for `deformed`
 * @param {number|'last'} [o.step='last'] index into `steps` used for the summary fields
 * @param {number} [o.scale=1] magnification of the displacements in `deformed` (-)
 * @returns {{steps:{index:number, time:number, fields:Object<string,{components:string[], values:Float64Array}>}[],
 *   nodeCount:number, displacement:Float64Array|null, stress:Float64Array|null, vonMises?:Float64Array,
 *   nodeValues:number[]|null, valueLabel:string|null, valueUnit:string|null, deformed:number[][]|null,
 *   time:number}} steps: every result set in file order (values: components × node, row per node; NaN where
 *   the file has no value); displacement 3 per node (m); stress 6 per node (xx, yy, zz, xy, yz, zx; Pa);
 *   nodeValues: von Mises stress (Pa) when the step has stresses, otherwise displacement magnitude (m);
 *   deformed: [x, y, z] per node (m); time: step time or frequency as written by ccx
 */
export function readCalculixFrd(text, { nodeCount, positions, step = 'last', scale = 1 } = {}) {
  const lines = String(text).split(/\r?\n/);
  let maxId = 0;
  const raw = []; // { time, name, components, ids[], rows[] }
  let i = 0; let time = 0; let stepNo = 0;
  const parseData = (line) => {
    const id = Number(line.slice(3, 13));
    const vals = [];
    for (let c = 13; c + 12 <= line.length || (c < line.length && line.slice(c).trim()); c += 12) vals.push(Number(line.slice(c, c + 12)));
    return { id, vals };
  };
  while (i < lines.length) {
    const line = lines[i];
    if (/^\s{4}2C/.test(line)) {
      for (i++; i < lines.length && !/^\s-3/.test(lines[i]); i++) if (/^\s-1/.test(lines[i])) { const id = Number(lines[i].slice(3, 13)); if (id > maxId) maxId = id; }
    } else if (/^\s{2}100C/.test(line)) {
      time = Number(line.slice(12, 24));
      const sn = Number(line.slice(58, 63));
      stepNo = Number.isFinite(sn) && sn > 0 ? sn : stepNo + 1;
    } else if (/^\s-4/.test(line)) {
      const tk = line.trim().split(/\s+/);
      const blk = { time, stepNo, name: tk[1], components: [], ids: [], rows: [] };
      for (i++; i < lines.length && /^\s-5/.test(lines[i]); i++) blk.components.push(lines[i].trim().split(/\s+/)[1]);
      // the trailing "ALL" pseudo-component of vector blocks carries no data column
      if (blk.components[blk.components.length - 1] === 'ALL') blk.components.pop();
      let cur = null;
      for (; i < lines.length && !/^\s-3/.test(lines[i]); i++) {
        if (/^\s-1/.test(lines[i])) { cur = parseData(lines[i]); blk.ids.push(cur.id); blk.rows.push(cur.vals); if (cur.id > maxId) maxId = cur.id; } else if (/^\s-2/.test(lines[i]) && cur) for (const v of parseData(lines[i]).vals) cur.vals.push(v);
      }
      raw.push(blk);
    }
    i++;
  }
  if (!raw.length) throw new RangeError('frd: no result blocks (-4 records) found; is this an ASCII .frd file?');
  const n = nodeCount === undefined ? maxId : nodeCount;
  if (!Number.isInteger(n) || n < 1) throw new RangeError('nodeCount must be a positive integer');
  const steps = [];
  for (const blk of raw) {
    let s = steps[steps.length - 1];
    if (!s || s.stepNo !== blk.stepNo || s.time !== blk.time || s.fields[blk.name]) { s = { index: steps.length, stepNo: blk.stepNo, time: blk.time, fields: {} }; steps.push(s); }
    const nc = blk.components.length;
    const values = new Float64Array(nc * n).fill(NaN);
    blk.ids.forEach((id, r) => { if (id >= 1 && id <= n) for (let c = 0; c < nc; c++) values[(id - 1) * nc + c] = blk.rows[r][c]; });
    s.fields[blk.name] = { components: blk.components, values };
  }
  for (const s of steps) delete s.stepNo;
  const k = step === 'last' ? steps.length - 1 : step;
  if (!Number.isInteger(k) || k < 0 || k >= steps.length) throw new RangeError(`step must be 'last' or an index 0…${steps.length - 1}`);
  const f = steps[k].fields;
  const out = { steps, nodeCount: n, time: steps[k].time, displacement: null, stress: null };
  if (f.DISP) { out.displacement = new Float64Array(3 * n); for (let q = 0; q < n; q++) for (let c = 0; c < 3; c++) out.displacement[3 * q + c] = f.DISP.values[q * f.DISP.components.length + c]; }
  if (f.STRESS && f.STRESS.components.length >= 6) { out.stress = new Float64Array(6 * n); for (let q = 0; q < n; q++) for (let c = 0; c < 6; c++) out.stress[6 * q + c] = f.STRESS.values[q * f.STRESS.components.length + c]; }
  viewerFields(out, flatPositions(positions), n, scale);
  return out;
}

/**
 * Reads the printed tables of a CalculiX .dat file: nodal displacements (*NODE PRINT U), reaction forces
 * (*NODE PRINT RF, with totals) and stresses at the integration points (*EL PRINT S).
 * Reference: CalculiX CrunchiX User's Manual v2.21, *NODE PRINT and *EL PRINT (table headers
 * "displacements (vx,vy,vz) for set … and time …", "forces (fx,fy,fz) …", "total force (fx,fy,fz) …",
 * "stresses (elem, integ.pnt.,sxx,syy,szz,sxy,sxz,syz) …").
 * Validity: tables in the global coordinate system; one row per node or integration point.
 * @param {string} text file content
 * @param {object} [o]
 * @param {number} [o.nodeCount] number of app nodes (default: largest node id in the displacement tables)
 * @param {ArrayLike<number>|object} [o.positions] undeformed coordinates (flat xyz, m) or the model
 * @param {number} [o.scale=1] magnification of the displacements in `deformed` (-)
 * @returns {{tables:{kind:string, set:string, time:number, rows:number[][]}[], nodeCount:number,
 *   displacement:Float64Array|null, reactions:{node:number, F:number[]}[], totalForce:number[]|null,
 *   elementStress:{element:number[], vonMises:number[]}|null, nodeValues:number[]|null, valueLabel:string|null,
 *   valueUnit:string|null, deformed:number[][]|null, time:number|null}} tables in file order; the summary fields
 *   come from the last table of each kind: displacement 3 per node (m, NaN where not printed), reactions (N, node
 *   as 0-based app index), total reaction (N), mean von Mises stress per element over its integration points
 *   (Pa, element as 0-based index), nodeValues = displacement magnitude (m)
 */
export function readCalculixDat(text, { nodeCount, positions, scale = 1 } = {}) {
  const tables = [];
  let cur = null;
  for (const line of String(text).split(/\r?\n/)) {
    const h = /^\s*([A-Za-z][A-Za-z .]*?)\s*\(([^)]*)\)\s*for set\s+(\S+)\s+and time\s+(\S+)/i.exec(line);
    if (h) { cur = { kind: h[1].trim().toLowerCase(), columns: h[2], set: h[3], time: Number(h[4]), rows: [] }; tables.push(cur); continue; }
    if (!cur) continue;
    const t = line.trim();
    if (!t) continue;
    const v = t.split(/\s+/).map(Number);
    if (v.every((x) => Number.isFinite(x))) cur.rows.push(v); else cur = null;
  }
  if (!tables.length) throw new RangeError('dat: no result tables found (expected lines such as "displacements (vx,vy,vz) for set NALL and time …")');
  const last = (kind) => { for (let k = tables.length - 1; k >= 0; k--) if (tables[k].kind === kind) return tables[k]; return null; };
  const d = last('displacements');
  let n = nodeCount;
  if (n === undefined) { n = 0; if (d) for (const r of d.rows) n = Math.max(n, r[0]); }
  const out = { tables, nodeCount: n, displacement: null, stress: null, reactions: [], totalForce: null, elementStress: null, time: d ? d.time : null };
  if (d && n > 0) {
    out.displacement = new Float64Array(3 * n).fill(NaN);
    for (const r of d.rows) if (r[0] >= 1 && r[0] <= n) for (let c = 0; c < 3; c++) out.displacement[3 * (r[0] - 1) + c] = r[1 + c];
  }
  const rf = last('forces');
  if (rf) out.reactions = rf.rows.map((r) => ({ node: r[0] - 1, F: [r[1], r[2], r[3]] }));
  const tf = last('total force');
  if (tf && tf.rows.length) out.totalForce = tf.rows[0].slice(0, 3);
  const s = last('stresses');
  if (s) {
    const acc = new Map();
    for (const r of s.rows) {
      if (r.length < 8) continue;
      // .dat order is sxx, syy, szz, sxy, sxz, syz; the invariant does not depend on the order of the shears
      const vm = vonMises(r, 2);
      const a = acc.get(r[0]) || [0, 0];
      a[0] += vm; a[1]++; acc.set(r[0], a);
    }
    const ids = Array.from(acc.keys()).sort((a, b) => a - b);
    out.elementStress = { element: ids.map((e) => e - 1), vonMises: ids.map((e) => acc.get(e)[0] / acc.get(e)[1]) };
  }
  viewerFields(out, flatPositions(positions), n, scale);
  delete out.stress;
  return out;
}
