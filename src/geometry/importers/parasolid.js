// Parasolid XT, text transmit files (.x_t, .xmt_txt) — boundary representation.
//
// Layout (Siemens, "Parasolid XT Format Reference"): a header of "**" lines up to "**END_OF_HEADER", then one
// stream in which line breaks mean nothing:
//   "T", the modeller version and the schema name as counted strings, and the user-field size;
//   nodes  "<node type> [<length of the variable field>] <index> <fields …> [<user field>]" in any order;
//   the terminator "1 0".
// Integers and reals are followed by one blank; characters and logicals (T / F) are not; "?" is an unset value
// (one "?" for a whole vector). Pointers are node indices, 0 = null.
//
// The order and types of the fields of every node type are given by the schema. Files of Parasolid 14 and later
// name a base schema (SCH_<version>_<schema>_13006) and carry, in front of the first node of each type, either
// the flag 255 ("as in the base schema") or an edit sequence against it: C copy, D delete, I insert and A append
// a field (each with its name, pointer class, element count and type), Z end; a node type the base schema does
// not have is described in full. The table BASE below is that base schema, entered from the reference; fields
// added by later releases are therefore read (and passed over) without knowing them.
// Files of Parasolid 9 – 13 (and later files saved without the base-schema suffix) carry no such description;
// they are read with the layouts listed under VARIANT, tried in turn, and are refused when the stream does not
// come out consistent under any of them.
//
// Read: bodies (solid, sheet, wire and general) through region → shell → face → loop → fin → edge → vertex → point;
// part files holding several bodies; assemblies, with the transforms of their instances.
//   surfaces  PLANE, CYLINDER, CONE, SPHERE, TORUS (apple included), B_SURFACE (NURBS, trimmed by the face
//             boundary), SWEPT_SURF over a line or B-curve
//   curves    LINE, CIRCLE, ELLIPSE, B_CURVE, TRIMMED_CURVE over one of these, INTERSECTION (the chart of points
//             the file stores for it, as a polyline), SP_CURVE on a surface that is read
// Each face is triangulated in the parameter space of its surface between its loops (brep.js), with the points of
// the shared edges, so a solid whose faces are all read comes out watertight. Coordinates are metres: the
// Parasolid kernel works in metres whatever unit the originating program showed.
//
// Conventions relied on:
//   * the face normal is the natural normal ∂u × ∂v of its surface when the senses of face and surface are equal,
//     and its opposite otherwise; a loop runs with the face on its left seen against the face normal;
//   * an edge runs along its curve in the direction of increasing parameter when the curve sense is "+";
//     edge.fin is a fin along the edge, and fin.vertex is the vertex a fin runs towards;
//   * CONE: R(u, v) = P + v·A + (X cos u + Y sin u)(r + v·tan a) — the radius grows along the axis, and the
//     natural normal points away from it. (The reference prints "P − v·A"; every cone face of the files checked
//     has its boundary on the surface only with "+", and is oriented consistently with its neighbours only so.)
//   * rational control points are stored multiplied by their weight; a B-surface lists its vertices with v
//     running fastest.
//
// Not read: binary transmit files (.x_b); faces on BLENDED_EDGE, OFFSET_SURF, SPUN_SURF, foreign (PE) surfaces
// and SWEPT_SURF over other curves — counted per type and reported; periodic B-splines whose knots are not
// stored in the unwrapped form; attributes (names, colours), groups.
// A node that cannot be interpreted stops the file: nothing is guessed.

import { inputText, baseName, finish, surfaceModel, Welder, V3, tally, cleanName, Warnings } from './common.js';
import { sampleCurve, findSpan, basisFuns } from './nurbs.js';
import { stepper, planeSurface, coneSurface, sphereSurface, torusSurface, splineSurface, tessellateFace, BodyMesh, CROSSING } from './brep.js';

const TAU = 2 * Math.PI, D2R = Math.PI / 180;
const dist3 = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);

// ---------------------------------------------------------------------------------------------------------
// Base schema SCH_13006: node type → [name, fields]. A field is "name:type"; "type*" is the variable-length
// field, "typeN" an array of N. Types: d int, n short, u byte, w unicode, t tag, p pointer, f double, c char,
// l logical, v vector, i interval, b box, h point of an intersection chart.
// ---------------------------------------------------------------------------------------------------------

const CURVE_HEAD = 'node_id:d attributes_groups:p owner:p next:p previous:p geometric_owner:p sense:c';
const BASE = {
  10: ['ASSEMBLY', 'highest_node_id:d attributes_groups:p attribute_chains:p list:p surface:p curve:p point:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p type:u sub_instance:p'],
  11: ['INSTANCE', 'node_id:d attributes_groups:p type:u part:p transform:p assembly:p next_in_part:p prev_in_part:p next_of_part:p prev_of_part:p'],
  12: ['BODY', 'highest_node_id:d attributes_groups:p attribute_chains:p surface:p curve:p point:p key:p res_size:f res_linear:f ref_instance:p next:p previous:p state:u owner:p body_type:u nom_geom_state:u shell:p boundary_surface:p boundary_curve:p boundary_point:p region:p edge:p vertex:p'],
  13: ['SHELL', 'node_id:d attributes_groups:p body:p next:p face:p edge:p vertex:p region:p front_face:p'],
  14: ['FACE', 'node_id:d attributes_groups:p tolerance:f next:p previous:p loop:p shell:p surface:p sense:c next_on_surface:p previous_on_surface:p next_front:p previous_front:p front_shell:p'],
  15: ['LOOP', 'node_id:d attributes_groups:p fin:p face:p next:p'],
  16: ['EDGE', 'node_id:d attributes_groups:p tolerance:f fin:p previous:p next:p curve:p next_on_curve:p previous_on_curve:p owner:p'],
  17: ['FIN', 'attributes_groups:p loop:p forward:p backward:p vertex:p other:p edge:p curve:p next_at_vx:p sense:c'],
  18: ['VERTEX', 'node_id:d attributes_groups:p fin:p previous:p next:p point:p tolerance:f owner:p'],
  19: ['REGION', 'node_id:d attributes_groups:p body:p next:p previous:p shell:p type:c'],
  29: ['POINT', 'node_id:d attributes_groups:p owner:p next:p previous:p pvec:v'],
  30: ['LINE', `${CURVE_HEAD} pvec:v direction:v`],
  31: ['CIRCLE', `${CURVE_HEAD} centre:v normal:v x_axis:v radius:f`],
  32: ['ELLIPSE', `${CURVE_HEAD} centre:v normal:v x_axis:v major_radius:f minor_radius:f`],
  38: ['INTERSECTION', `${CURVE_HEAD} surface:p2 chart:p start:p end:p`],
  40: ['CHART', 'base_parameter:f base_scale:f chart_count:d chordal_error:f angular_error:f parameter_error:f2 hvec:h*'],
  41: ['LIMIT', 'type:c hvec:h*'],
  45: ['BSPLINE_VERTICES', 'vertices:f*'],
  50: ['PLANE', `${CURVE_HEAD} pvec:v normal:v x_axis:v`],
  51: ['CYLINDER', `${CURVE_HEAD} pvec:v axis:v radius:f x_axis:v`],
  52: ['CONE', `${CURVE_HEAD} pvec:v axis:v radius:f sin_half_angle:f cos_half_angle:f x_axis:v`],
  53: ['SPHERE', `${CURVE_HEAD} centre:v radius:f axis:v x_axis:v`],
  54: ['TORUS', `${CURVE_HEAD} centre:v axis:v major_radius:f minor_radius:f x_axis:v`],
  56: ['BLENDED_EDGE', `${CURVE_HEAD} blend_type:c surface:p2 spine:p range:f2 thumb_weight:f2 boundary:p2 start:p end:p`],
  59: ['BLEND_BOUND', `${CURVE_HEAD} boundary:n blend:p`],
  60: ['OFFSET_SURF', `${CURVE_HEAD} check:c true_offset:l surface:p offset:f scale:f`],
  67: ['SWEPT_SURF', `${CURVE_HEAD} section:p sweep:v scale:f`],
  68: ['SPUN_SURF', `${CURVE_HEAD} profile:p base:v axis:v start:v end:v start_param:f end_param:f x_axis:v scale:f`],
  70: ['LIST', 'node_id:d owner:p next:p previous:p list_type:d list_length:d block_length:d size_of_entry:d list_block:p finger_block:p finger_index:d notransmit:l'],
  74: ['POINTER_LIS_BLOCK', 'n_entries:d next_block:p entries:p*'],
  79: ['ATT_DEF_ID', 'string:c*'],
  80: ['ATTRIB_DEF', 'next:p identifier:p type_id:d actions:u8 field_names:p legal_owners:l14 fields:u*'],
  81: ['ATTRIBUTE', 'node_id:d definition:p owner:p next:p previous:p next_of_type:p previous_of_type:p fields:p*'],
  82: ['INT_VALUES', 'values:d*'],
  83: ['REAL_VALUES', 'values:f*'],
  84: ['CHAR_VALUES', 'values:c*'],
  85: ['POINT_VALUES', 'values:v*'],
  86: ['VECTOR_VALUES', 'values:v*'],
  87: ['AXIS_VALUES', 'values:v*'],
  88: ['TAG_VALUES', 'values:t*'],
  89: ['DIRECTION_VALUES', 'values:v*'],
  90: ['GROUP', 'node_id:d attributes_groups:p owner:p next:p previous:p type:u first_member:p'],
  91: ['MEMBER_OF_GROUP', 'dummy_node_id:d owning_group:p owner:p next:p previous:p next_member:p previous_member:p'],
  98: ['UNICODE_VALUES', 'values:w*'],
  99: ['FIELD_NAMES', 'names:p*'],
  100: ['TRANSFORM', 'node_id:d owner:p next:p previous:p rotation_matrix:f9 translation_vector:v scale:f flag:d perspective_vector:v'],
  101: ['WORLD', 'assembly:p attribute:p body:p transform:p surface:p curve:p point:p alive:l attrib_def:p highest_id:d current_id:d'],
  102: ['KEY', 'string:c*'],
  120: ['PE_SURF', `${CURVE_HEAD} type:c data:p tf:p internal_geom:p*`],
  121: ['INT_PE_DATA', 'geom_type:d real_array:p int_array:p'],
  122: ['EXT_PE_DATA', 'key:p real_array:p int_array:p'],
  124: ['B_SURFACE', `${CURVE_HEAD} nurbs:p data:p`],
  125: ['SURFACE_DATA', 'original_uint:i original_vint:i extended_uint:i extended_vint:i self_int:u original_u_start:c original_u_end:c original_v_start:c original_v_end:c extended_u_start:c extended_u_end:c extended_v_start:c extended_v_end:c analytic_form_type:c swept_form_type:c spun_form_type:c blend_form_type:c analytic_form:p swept_form:p spun_form:p blend_form:p'],
  126: ['NURBS_SURF', 'u_periodic:l v_periodic:l u_degree:n v_degree:n n_u_vertices:d n_v_vertices:d u_knot_type:u v_knot_type:u n_u_knots:d n_v_knots:d rational:l u_closed:l v_closed:l surface_form:u vertex_dim:n bspline_vertices:p u_knot_mult:p v_knot_mult:p u_knots:p v_knots:p'],
  127: ['KNOT_MULT', 'mult:n*'],
  128: ['KNOT_SET', 'knots:f*'],
  130: ['PE_CURVE', `${CURVE_HEAD} type:c data:p tf:p internal_geom:p*`],
  133: ['TRIMMED_CURVE', `${CURVE_HEAD} basis_curve:p point_1:v point_2:v parm_1:f parm_2:f`],
  134: ['B_CURVE', `${CURVE_HEAD} nurbs:p data:p`],
  135: ['CURVE_DATA', 'self_int:u analytic_form:p'],
  136: ['NURBS_CURVE', 'degree:n n_vertices:d vertex_dim:n n_knots:d knot_type:u periodic:l closed:l rational:l curve_form:u bspline_vertices:p knot_mult:p knots:p'],
  137: ['SP_CURVE', `${CURVE_HEAD} surface:p b_curve:p original:p tolerance_to_original:f`],
  141: ['GEOMETRIC_OWNER', 'owner:p next:p previous:p shared_geometry:p'],
  163: ['HELIX_CU_FORM', 'axis_pt:v axis_dir:v point:v hand:c turns:i pitch:f tol:f'],
  184: ['HELIX_SU_FORM', 'axis_pt:v axis_dir:v hand:c turns:i pitch:f gap:f tol:f'],
};

const SURFACE_TYPES = new Set([50, 51, 52, 53, 54, 56, 59, 60, 67, 68, 120, 124]);
const CURVE_TYPES = new Set([30, 31, 32, 38, 130, 133, 134, 137]);
const parseFields = (text) => text.split(' ').map((f) => {
  const [name, t] = f.split(':');
  const type = t[0], rest = t.slice(1);
  return { name, type, n: rest === '*' ? 1 : rest ? Number(rest) : 0 };
});
const baseFields = new Map();
const baseOf = (type) => {
  if (!BASE[type]) return null;
  if (!baseFields.has(type)) baseFields.set(type, parseFields(BASE[type][1]));
  return baseFields.get(type);
};

/**
 * Schemas other than the base one, for files that do not describe their schema themselves (everything before
 * Parasolid 14, and later files saved "as version n" without the base-schema suffix). Each entry edits the base
 * field list of one node type. Which of them apply to a schema number is known only in part, so the candidates
 * below are tried in turn and the first under which the whole node stream reads consistently — every node type
 * known, every index unique, every value of its type, the terminator at the end — is taken; if none does, the
 * file is refused.
 *   old-attdef   attribute definitions before 12.1: no field names, 13 owner classes          (reference, 12006)
 *   no-nom-geom  bodies of Parasolid 9: no nominal-geometry state                              (files, 9003 / 9008)
 *   bare-fin     fins of Parasolid 9: no attribute chain                                       (files, 9003 / 9008)
 *   new-list     list header and pointer block as written since 14 (embedded schemas 19008 …)
 *   body-maps    the four index-map fields of a body (embedded schemas 19008, 20000)
 */
const VARIANT = {
  'old-attdef': [80, (f) => f.filter((x) => x.name !== 'field_names').map((x) => (x.name === 'legal_owners' ? { ...x, n: 13 } : x))],
  'no-nom-geom': [12, (f) => f.filter((x) => x.name !== 'nom_geom_state')],
  'bare-fin': [17, (f) => f.filter((x) => x.name !== 'attributes_groups')],
  'new-list': [70, () => parseFields('node_id:d list_type:u notransmit:l owner:p next:p previous:p list_length:d block_length:d finger_index:d finger_block:p list_block:p')],
  'new-block': [74, () => parseFields('n_entries:d index_map_offset:d next_block:p entries:p*')],
  'body-maps': [12, (f) => f.concat(parseFields('index_map_offset:d index_map:p node_id_index_map:p schema_embedding_map:p'))],
  // later releases, as their embedded schemas describe them (number of the first schema seen with the change)
  'body-25001': [12, (f) => f.concat(parseFields('child:p lowest_node_id:d'))],
  'region-25001': [19, (f) => f.concat(parseFields('owner:p'))],
  'body-26105': [12, (f) => { const k = f.findIndex((x) => x.name === 'index_map_offset'); return f.slice(0, k).concat(parseFields('boundary_mesh:p'), f.slice(k)); }],
  'body-28101': [12, (f) => { const a = f.findIndex((x) => x.name === 'key'), b = f.findIndex((x) => x.name === 'region'); return f.slice(0, a).concat(parseFields('mesh:p polyline:p'), f.slice(a, b), parseFields('boundary_polyline:p'), f.slice(b)); }],
  'limit-28101': [41, () => parseFields('type:c term_use:c hvec:h*')],
  'assembly-28101': [10, (f) => { const a = f.findIndex((x) => x.name === 'key'); return f.slice(0, a).concat(parseFields('mesh:p polyline:p'), f.slice(a)); }],
  'intersection-30000': [38, (f) => f.concat(parseFields('intersection_data:p'))],
  'intersection-data': [204, () => parseFields('uv_type:u values:f*')],
  'body-31001': [12, (f) => f.concat(parseFields('mesh_offset_data:p'))],
  'assembly-31001': [10, (f) => f.concat(parseFields('mesh_offset_data:p'))],
};
const LATER = [[25001, ['body-25001', 'region-25001']], [26105, ['body-26105']], [28101, ['body-28101', 'limit-28101', 'assembly-28101']], [30000, ['intersection-30000', 'intersection-data']], [31001, ['body-31001', 'assembly-31001']]];
function candidates(schemaNumber) {
  if (schemaNumber >= 13006) {
    const full = ['new-list', 'new-block', 'body-maps'];
    for (const [from, names] of LATER) if (schemaNumber >= from) full.push(...names);
    return [full, ['new-list', 'new-block', 'body-maps'], ['new-list', 'new-block'], ['new-list'], []];
  }
  if (schemaNumber >= 12100) return [[], ['old-attdef']];
  if (schemaNumber >= 10000) return [['old-attdef'], [], ['old-attdef', 'bare-fin'], ['old-attdef', 'no-nom-geom', 'bare-fin']];
  return [['old-attdef', 'no-nom-geom', 'bare-fin'], ['old-attdef', 'bare-fin'], ['old-attdef', 'no-nom-geom'], ['old-attdef']];
}

class Stop extends Error {}

/**
 * Decode a text transmit file into its nodes.
 * @returns {{ header: object, nodes: Map<number, object>, types: object, userField: number }} every node is
 *   `{ index, type, name, …fields }` with pointers left as indices
 */
export function readXt(text, fname = 'x_t') {
  const fail = (why) => new Error(`${fname}: ${why}`);
  const mark = text.indexOf('**END_OF_HEADER');
  if (!text.startsWith('**') || mark < 0) throw fail('not a Parasolid text transmit file (the "**…" header ending in "**END_OF_HEADER" is missing).');
  const head = text.slice(0, mark).replace(/[\r\n]/g, '');
  const key = (k) => { const m = new RegExp(`[;*]${k}=([^;]*);`).exec(head); return m ? m[1] : null; };
  let nl = text.indexOf('\n', mark);
  if (nl < 0) nl = text.length;
  const s = text.slice(nl + 1).replace(/[\r\n]/g, '');
  let p = 0;
  const n = s.length;
  const where = () => `near character ${p} of the data ("${s.slice(Math.max(0, p - 12), p + 12)}")`;
  const end = (what) => new Stop(`the file ends inside ${what}`);
  const int = (what) => {
    if (p >= n) throw end(what);
    if (s[p] === '?') { p++; return null; }
    const a = p;
    let q = p;
    if (s[q] === '-' || s[q] === '+') q++;
    while (q < n && s.charCodeAt(q) >= 48 && s.charCodeAt(q) <= 57) q++;
    if (q === a || (q < n && s[q] !== ' ') || ((s[a] === '-' || s[a] === '+') && q === a + 1)) throw new Stop(`an integer was expected for ${what} ${where()}`);
    p = q + 1;
    return Number(s.slice(a, q));
  };
  const real = (what) => {
    if (p >= n) throw end(what);
    if (s[p] === '?') { p++; return null; }
    let q = s.indexOf(' ', p);
    if (q < 0) q = n;
    const w = s.slice(p, q);
    if (!/^[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$/.test(w)) throw new Stop(`a number was expected for ${what} ${where()}`);
    p = q + 1;
    return Number(w);
  };
  const vec = (what) => { if (s[p] === '?') { p++; return null; } return [real(what), real(what), real(what)]; };
  const chr = (what) => { if (p >= n) throw end(what); return s[p++]; };
  const str = (what) => { const len = int(what); if (!(len >= 0) || p + len > n) throw new Stop(`a counted string was expected for ${what} ${where()}`); const v = s.slice(p, p + len); p += len; return v; };

  // ---- flag sequence, versions, user-field size ----
  if (s[p] !== 'T') throw fail(s[p] === 'B' || s.startsWith('PS') ? 'this is a binary Parasolid transmit file, which is not read. Save it as text (.x_t) or export STEP.' : 'the data after the header do not start with the text flag "T".');
  p++;
  let modeller, schema, userField, header;
  try {
    modeller = str('the modeller version');
    schema = str('the schema name');
    const m = /^SCH_(\d+)_(\d+)(?:_(\d+))?$/.exec(schema);
    if (!m) throw new Stop(`the schema name "${schema}" is not of the form SCH_<version>_<schema>`);
    if (m[3]) int('the number of node types');
    userField = int('the user-field size');
    header = { modeller: /version (\d+)/.test(modeller) ? Number(/version (\d+)/.exec(modeller)[1]) : null, schema, schemaNumber: Number(m[2]), baseSchema: m[3] ? Number(m[3]) : null, application: key('APPL'), date: key('DATE'), format: key('FORMAT'), partition: /partition/.test(modeller) };
  } catch (e) { if (e instanceof Stop) throw fail(`the start of the data could not be read: ${e.message}.`); throw e; }
  if (header.baseSchema !== null && header.baseSchema !== 13006) throw fail(`the file is written against the base schema ${header.baseSchema}; only 13006 is known to this reader.`);
  if (!(userField >= 0 && userField <= 16)) throw fail(`the user-field size ${userField} is outside 0 … 16.`);
  const embedded = header.baseSchema !== null;

  // ---- schema in force for this file ----
  const fieldsOf = new Map();
  const names = new Map();
  const legacy = (type, variant) => {
    let f = baseOf(type) || (variant.some((v) => VARIANT[v][0] === type) ? [] : null);
    if (!f) return null;
    for (const v of variant) if (VARIANT[v][0] === type) f = VARIANT[v][1](f);
    return f;
  };
  const descriptor = () => {
    const name = str('a field name');
    const cls = int('a pointer class'), cnt = int('an element count');
    const type = cls ? 'p' : str('a field type');
    const xmt = cnt === 1 ? chr('a transmit flag') : 'T';
    if (xmt !== 'T' && xmt !== 'F') throw new Stop(`the transmit flag of field "${name}" is "${xmt}" ${where()}`);
    return { name, type: type[0] || '?', n: cnt, skip: xmt === 'F' };
  };
  const readSchema = (type) => {
    const count = int('the schema of a node type');
    const base = baseOf(type);
    if (count === 255) { if (!base) throw new Stop(`node type ${type} is declared equal to the base schema, which has no such type`); return base; }
    const out = [];
    if (!base) {
      names.set(type, str('a node name'));
      str('a node description');
      for (let k = 0; k < count; k++) out.push(descriptor());
      return out;
    }
    let b = 0;
    for (;;) {
      const op = chr('a schema edit');
      if (op === 'Z') break;
      if (op === 'C') { if (b >= base.length) throw new Stop(`the schema edit of node type ${type} copies more fields than the base schema has`); out.push(base[b++]); } else if (op === 'D') b++;
      else if (op === 'I' || op === 'A') out.push(descriptor());
      else throw new Stop(`unknown schema edit "${op}" for node type ${type} ${where()}`);
    }
    if (out.length !== count) throw new Stop(`the schema edit of node type ${type} gives ${out.length} fields where ${count} are announced`);
    return out;
  };

  const value = (type, what) => {
    switch (type) {
      case 'd': case 'n': case 'u': case 'w': case 't': case 'p': return int(what);
      case 'f': return real(what);
      case 'c': return chr(what);
      case 'l': { const c = chr(what); if (c !== 'T' && c !== 'F') throw new Stop(`a logical (T / F) was expected for ${what} ${where()}`); return c === 'T'; }
      case 'v': case 'h': return vec(what);
      case 'i': return [real(what), real(what)];
      case 'b': return [real(what), real(what), real(what), real(what), real(what), real(what)];
      default: throw new Stop(`field type "${type}" of ${what} is not known`);
    }
  };

  const start = p;
  /** Read the node stream under one schema variant; throws Stop where it does not fit. */
  const run = (variant) => {
    p = start;
    fieldsOf.clear(); names.clear();
    const nodes = new Map(), types = {};
    let done = false;
    try {
        while (p < n) {
          const type = int('a node type');
          if (type === 1) { const z = int('the terminator'); if (z !== 0) throw new Stop(`the terminator is followed by ${z} instead of 0`); done = true; break; }
          if (!(type > 1 && type < 100000)) throw new Stop(`${type} is not a node type ${where()}`);
          let fields = fieldsOf.get(type);
          if (!fields) {
            fields = embedded ? readSchema(type) : legacy(type, variant);
            if (!fields) throw new Stop(`node type ${type} is not in the schema of this reader and the file does not describe it`);
            fieldsOf.set(type, fields);
          }
          const name = BASE[type] ? BASE[type][0] : names.get(type) || (type === 204 ? 'INTERSECTION_DATA' : `type ${type}`);
          const variable = fields.some((f) => f.n === 1);
          const len = variable ? int(`the length of a ${name} node`) : 0;
          if (variable && !(len >= 0 && len <= n)) throw new Stop(`the variable length ${len} of a ${name} node is not plausible ${where()}`);
          const index = int(`the index of a ${name} node`);
          if (!(index > 0)) throw new Stop(`a ${name} node has the index ${index} ${where()}`);
          if (nodes.has(index)) throw new Stop(`the index ${index} is used by two nodes (${nodes.get(index).name} and ${name})`);
          const node = { index, type, name };
          for (const f of fields) {
            if (f.skip) continue;
            const what = `${name}.${f.name}`;
          // fields called like the bookkeeping of this reader are stored under another name
          const key = f.name === 'type' ? 'kind' : f.name === 'index' || f.name === 'name' ? `${f.name}_` : f.name;
            if (f.n === 0) node[key] = value(f.type, what);
            else {
              const cnt = f.n === 1 ? len : f.n;
              if (f.type === 'c') {
                // character arrays are written as text; "\9" stands for nine blanks
                let v = '';
                while (v.length < cnt) { const c = chr(what); if (c === '\\' && s[p] === '9') { p++; v += '         '; } else if (c === '\\' && p < n) { const e = s[p++]; v += e === 'n' ? '\n' : e === 'r' ? '\r' : e === '0' ? '\0' : e; } else v += c; }
                node[key] = v;
              } else {
                const a = new Array(cnt);
                for (let k = 0; k < cnt; k++) a[k] = value(f.type, what);
                node[key] = a;
              }
            }
          }
          // the user field follows every node (the reference speaks of the visible ones only; files written by NX
          // with a user-field size of 1 have it after all of them)
          for (let k = 0; k < userField; k++) int(`the user field of a ${name} node`);
          nodes.set(index, node);
          types[name] = (types[name] || 0) + 1;
        }
    } catch (e) {
      if (e instanceof Stop) e.count = nodes.size;
      throw e;
    }
    return { nodes, types, done };
  };
  let res = null, worst = null;
  const tried = embedded ? [[]] : candidates(header.schemaNumber);
  for (const variant of tried) {
    try { res = run(variant); res.variant = variant; break; } catch (e) {
      if (!(e instanceof Stop)) throw e;
      if (!worst || e.count > worst.count) worst = e;
    }
  }
  if (!res) throw fail(`the node stream could not be interpreted (schema ${schema}${embedded ? '' : `, which the file does not describe; ${tried.length} known layouts were tried`}): after ${worst.count} nodes, ${worst.message}. Nothing was read from the file; export STEP from the originating program.`);
  const { nodes, types, done } = res;
  if (!done) throw fail(`the file is truncated: the node stream ends after ${nodes.size} nodes without its terminator.`);
  if (!nodes.size) throw fail('the file holds a header but no nodes.');
  return { header, nodes, types, userField, schema: fieldsOf, typeNames: names, layout: embedded ? 'described by the file' : res.variant.length ? res.variant.join(', ') : 'base schema' };
}

// ---------------------------------------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------------------------------------


class Unreadable extends Error {}

/** Knot vector written out from distinct values and multiplicities. */
function knotVector(knots, mult) {
  if (!knots || !mult || knots.length !== mult.length) return null;
  const U = [];
  for (let i = 0; i < knots.length; i++) { if (!(mult[i] >= 1) || knots[i] === null || (i && knots[i] < knots[i - 1])) return null; for (let k = 0; k < mult[i]; k++) U.push(knots[i]); }
  return U;
}

/** Control points of a B-spline node: `dim` values per vertex, the last one a weight (points stored multiplied by it) when rational. */
function controlNet(values, count, dim, rational) {
  const geo = rational ? dim - 1 : dim;
  if (!values || values.length < count * dim || geo < 2 || geo > 3) return null;
  const P = new Float64Array(count * 3), W = rational ? new Float64Array(count) : null;
  for (let i = 0; i < count; i++) {
    const w = rational ? values[i * dim + geo] : 1;
    if (!(w > 0)) return null;
    for (let d = 0; d < geo; d++) { const x = values[i * dim + d]; if (x === null) return null; P[3 * i + d] = x / w; }
    if (W) W[i] = w;
  }
  return { P, W, geo };
}

/**
 * Parasolid XT (text).
 * @param {{ name?: string, text?: string, bytes?: Uint8Array }} input
 * @param {{ deflection?: number, angularStep?: number }} [opts] `deflection`: largest distance between a curved
 *   edge or surface and its facets [m] (default 1 % of the bounding-box diagonal of each body); `angularStep`:
 *   largest angle one facet spans [deg], default 15
 */
export function parseParasolid(input, opts = {}) {
  const fname = input.name || 'x_t';
  if (input.bytes && input.bytes.length > 4 && input.bytes[0] === 0x50 && input.bytes[1] === 0x53 && input.bytes[2] === 0 && input.bytes[3] === 0) throw new Error(`${fname}: this is a binary Parasolid transmit file (.x_b), which is not read. Save it as text (.x_t) in the originating program, or export STEP.`);
  if (opts.deflection !== undefined && !(opts.deflection > 0)) throw new RangeError('deflection: must be a positive length in metres.');
  if (opts.angularStep !== undefined && !(opts.angularStep > 0 && opts.angularStep <= 90)) throw new RangeError('angularStep: must be within 0…90 degrees.');
  const xt = readXt(inputText(input), fname);
  const N = xt.nodes;
  const name = baseName(input.name);
  const W = new Warnings();
  const extra = [];
  const hist = {}, skipped = {}, faceSkips = {};
  const angStep = (opts.angularStep ?? 15) * D2R;
  const get = (i) => (i ? N.get(i) || null : null);
  const need = (i, ok, what) => { const q = get(i); if (!q || !ok(q)) throw new Unreadable(`${what} points to ${q ? `a ${q.name} node` : 'no node'}`); return q; };
  const isType = (...t) => (q) => t.includes(q.type);
  const vecOk = (v) => Array.isArray(v) && v.length === 3 && v.every((x) => typeof x === 'number' && Number.isFinite(x));
  const V = (node, field) => { const v = node[field]; if (!vecOk(v)) throw new Unreadable(`${node.name}.${field} is not set`); return v; };
  const F = (node, field) => { const v = node[field]; if (!(typeof v === 'number' && Number.isFinite(v))) throw new Unreadable(`${node.name}.${field} is not set`); return v; };

  // ---- evaluation of surfaces at Parasolid parameters, for SP-curves ----
  const frame = (axis, x) => { const A = V3.norm(axis), X = V3.norm(x); return { A, X, Y: V3.cross(A, X) }; };
  const nurbsSurfOf = (sn) => {
    if (sn.net !== undefined) return sn.net;
    const nb = get(sn.nurbs);
    if (!nb || nb.type !== 126) return (sn.net = null);
    const U = knotVector((get(nb.u_knots) || {}).knots, (get(nb.u_knot_mult) || {}).mult), Vk = knotVector((get(nb.v_knots) || {}).knots, (get(nb.v_knot_mult) || {}).mult);
    const nu = nb.n_u_vertices, nv = nb.n_v_vertices;
    const cn = controlNet((get(nb.bspline_vertices) || {}).vertices, nu * nv, nb.vertex_dim, nb.rational);
    if (!U || !Vk || !cn || cn.geo !== 3 || U.length !== nu + nb.u_degree + 1 || Vk.length !== nv + nb.v_degree + 1) return (sn.net = null);
    // the vertices run with v fastest: vertex (i, j) is entry i·nv + j; brep.js wants u fastest
    const P = new Float64Array(cn.P.length), Wt = cn.W ? new Float64Array(cn.W.length) : null;
    for (let iu = 0; iu < nu; iu++) for (let iv = 0; iv < nv; iv++) {
      const from = iv + nv * iu, to = iu + nu * iv;
      P[3 * to] = cn.P[3 * from]; P[3 * to + 1] = cn.P[3 * from + 1]; P[3 * to + 2] = cn.P[3 * from + 2];
      if (Wt) Wt[to] = cn.W[from];
    }
    return (sn.net = { pu: nb.u_degree, pv: nb.v_degree, U, V: Vk, nu, nv, P, W: Wt });
  };
  const netPoint = (net, u, v) => {
    const { pu, pv, U, V: Vk, nu, nv, P, W: Wt } = net;
    u = Math.max(U[pu], Math.min(U[nu], u)); v = Math.max(Vk[pv], Math.min(Vk[nv], v));
    const Nu = new Float64Array(pu + 1), Nv = new Float64Array(pv + 1), L = new Float64Array(Math.max(pu, pv) + 1), R = new Float64Array(Math.max(pu, pv) + 1);
    const su = findSpan(nu - 1, pu, u, U); basisFuns(su, u, pu, U, Nu, L, R);
    const sv = findSpan(nv - 1, pv, v, Vk); basisFuns(sv, v, pv, Vk, Nv, L, R);
    let x = 0, y = 0, z = 0, w = 0;
    for (let b = 0; b <= pv; b++) for (let a = 0; a <= pu; a++) {
      const i = su - pu + a + nu * (sv - pv + b), c = Nu[a] * Nv[b] * (Wt ? Wt[i] : 1);
      x += c * P[3 * i]; y += c * P[3 * i + 1]; z += c * P[3 * i + 2]; w += c;
    }
    if (!Wt || w === 0) w = 1;
    return [x / w, y / w, z / w];
  };
  /** (u, v) → point on a surface node in the parametrisation of the format, or null when the type is not evaluated. */
  const paramPoint = (sn) => {
    if (sn.at !== undefined) return sn.at;
    let at = null;
    const add = (o, ...terms) => { const p = o.slice(); for (const [k, d] of terms) { p[0] += k * d[0]; p[1] += k * d[1]; p[2] += k * d[2]; } return p; };
    try {
      if (sn.type === 50) { const O = V(sn, 'pvec'), Nn = V3.norm(V(sn, 'normal')), X = V3.norm(V(sn, 'x_axis')), Y = V3.cross(Nn, X); at = (u, v) => add(O, [u, X], [v, Y]); } else if (sn.type === 51) { const O = V(sn, 'pvec'), { A, X, Y } = frame(V(sn, 'axis'), V(sn, 'x_axis')), r = F(sn, 'radius'); at = (u, v) => add(O, [r * Math.cos(u), X], [r * Math.sin(u), Y], [v, A]); } else if (sn.type === 52) { const O = V(sn, 'pvec'), { A, X, Y } = frame(V(sn, 'axis'), V(sn, 'x_axis')), r = F(sn, 'radius'), t = F(sn, 'sin_half_angle') / F(sn, 'cos_half_angle'); at = (u, v) => add(O, [v, A], [(r + v * t) * Math.cos(u), X], [(r + v * t) * Math.sin(u), Y]); } else if (sn.type === 53) { const O = V(sn, 'centre'), { A, X, Y } = frame(V(sn, 'axis'), V(sn, 'x_axis')), r = F(sn, 'radius'); at = (u, v) => add(O, [r * Math.cos(v) * Math.cos(u), X], [r * Math.cos(v) * Math.sin(u), Y], [r * Math.sin(v), A]); } else if (sn.type === 54) { const O = V(sn, 'centre'), { A, X, Y } = frame(V(sn, 'axis'), V(sn, 'x_axis')), a = F(sn, 'major_radius'), b = F(sn, 'minor_radius'); at = (u, v) => add(O, [(a + b * Math.cos(v)) * Math.cos(u), X], [(a + b * Math.cos(v)) * Math.sin(u), Y], [b * Math.sin(v), A]); } else if (sn.type === 124) { const net = nurbsSurfOf(sn); if (net) at = (u, v) => netPoint(net, u, v); }
    } catch (e) { if (!(e instanceof Unreadable)) throw e; }
    return (sn.at = at);
  };

  // ---- curves: a polyline along the curve in the direction of its sense, and whether it closes ----
  const bcurveOf = (cn) => {
    const nb = get(cn.nurbs);
    if (!nb || nb.type !== 136) return null;
    const U = knotVector((get(nb.knots) || {}).knots, (get(nb.knot_mult) || {}).mult);
    const net = controlNet((get(nb.bspline_vertices) || {}).vertices, nb.n_vertices, nb.vertex_dim, nb.rational);
    if (!U || !net || U.length !== nb.n_vertices + nb.degree + 1) return null;
    return { deg: nb.degree, U, P: net.P, W: net.W, geo: net.geo };
  };
  const toPts = (flat) => { const o = []; for (let i = 0; i < flat.length; i += 3) o.push([flat[i], flat[i + 1], flat[i + 2]]); return o; };
  /** { pts, closed, exact } of a curve node, running with its sense; `why` when it cannot be followed. */
  const polylineOf = (cn, stepFor) => {
    if (cn.poly !== undefined) return cn.poly;
    const rev = cn.sense === '-';
    const done = (pts, closed, note) => (cn.poly = { pts: rev ? pts.slice().reverse() : pts, closed, note });
    if (cn.type === 134) {
      const bs = bcurveOf(cn);
      if (!bs || bs.geo !== 3) return (cn.poly = { why: 'a B-curve whose data are inconsistent' });
      const flat = sampleCurve(bs.deg, bs.U, bs.P, bs.W, 8);
      if (!flat) return (cn.poly = { why: 'a B-curve whose data are inconsistent' });
      const pts = toPts(flat);
      return done(pts, dist3(pts[0], pts[pts.length - 1]) <= 1e-9 * (1 + V3.len(pts[0])) + 1e-10);
    }
    if (cn.type === 38) {
      const chart = get(cn.chart);
      const pts = chart && chart.type === 40 ? chart.hvec.filter(vecOk) : [];
      if (pts.length < 2) return (cn.poly = { why: 'an intersection curve without chart points' });
      const a = get(cn.start), b = get(cn.end);
      return done(pts, !!(a && a.type === 41 && a.kind === 'H') || dist3(pts[0], pts[pts.length - 1]) <= 1e-12 * (1 + V3.len(pts[0])), 'chart');
    }
    if (cn.type === 133) {
      const basis = get(cn.basis_curve);
      if (!basis || !CURVE_TYPES.has(basis.type) || basis.type === 133) return (cn.poly = { why: 'a trimmed curve without a basis curve' });
      const inner = basis.type === 30 || analyticOf(basis) ? null : polylineOf(basis, stepFor);
      if (inner && inner.why) return (cn.poly = inner);
      return (cn.poly = { basis, inner });
    }
    if (cn.type === 137) {
      const sn = get(cn.surface), bc = get(cn.b_curve);
      const at = sn && SURFACE_TYPES.has(sn.type) ? paramPoint(sn) : null;
      const bs = bc && bc.type === 134 ? bcurveOf(bc) : null;
      if (!at) return (cn.poly = { why: `a curve in the parameter space of a ${sn ? sn.name : 'missing'} surface, which is not evaluated` });
      if (!bs || bs.geo !== 2) return (cn.poly = { why: 'a parameter-space curve whose 2-D B-curve is inconsistent' });
      const flat = sampleCurve(bs.deg, bs.U, bs.P, bs.W, 8);
      if (!flat) return (cn.poly = { why: 'a parameter-space curve whose 2-D B-curve is inconsistent' });
      const uv = toPts(flat);
      const pts = (bc.sense === '-' ? uv.reverse() : uv).map((q) => at(q[0], q[1]));
      return done(pts, dist3(pts[0], pts[pts.length - 1]) <= 1e-9 * (1 + V3.len(pts[0])) + 1e-10, 'sp');
    }
    return (cn.poly = { why: `a ${cn.name} curve, which is not read` });
  };
  /** Circle / ellipse data of a curve node, or null. */
  function analyticOf(cn) {
    if (cn.type !== 31 && cn.type !== 32) return null;
    if (cn.ana) return cn.ana;
    const C = V(cn, 'centre'), Nn = V3.norm(V(cn, 'normal')), X = V3.norm(V(cn, 'x_axis')), Y = V3.cross(Nn, X);
    const a = cn.type === 31 ? F(cn, 'radius') : F(cn, 'major_radius'), b = cn.type === 31 ? a : F(cn, 'minor_radius');
    if (!(a > 0 && b > 0)) throw new Unreadable(`${cn.name} with a radius that is not positive`);
    return (cn.ana = {
      a, b, dir: cn.sense === '-' ? -1 : 1,
      at: (t) => [C[0] + a * Math.cos(t) * X[0] + b * Math.sin(t) * Y[0], C[1] + a * Math.cos(t) * X[1] + b * Math.sin(t) * Y[1], C[2] + a * Math.cos(t) * X[2] + b * Math.sin(t) * Y[2]],
      ang: (p) => { const d = V3.sub(p, C); return Math.atan2(V3.dot(d, Y) / b, V3.dot(d, X) / a); },
      off: (p) => { const d = V3.sub(p, C), x = V3.dot(d, X) / a, y = V3.dot(d, Y) / b; return Math.hypot((Math.hypot(x, y) - 1) * Math.min(a, b), V3.dot(d, Nn)); },
    });
  }

  /** Position (segment index + fraction) of the point of a polyline nearest to p, and its distance. */
  const locate = (pts, p) => {
    let best = 0, bd = Infinity;
    for (let i = 0; i + 1 < pts.length; i++) {
      const a = pts[i], b = pts[i + 1], e = V3.sub(b, a), l2 = V3.dot(e, e);
      const f = l2 > 0 ? Math.max(0, Math.min(1, V3.dot(V3.sub(p, a), e) / l2)) : 0;
      const d = dist3(p, [a[0] + f * e[0], a[1] + f * e[1], a[2] + f * e[2]]);
      if (d < bd) { bd = d; best = i + f; }
    }
    return { s: best, d: bd };
  };

  // ---- topology ----
  const pointOfVertex = (v) => (v ? V(need(v.point, isType(29), 'VERTEX.point'), 'pvec') : null);
  /** Start and end vertices of an edge: the vertex of a fin against the edge, and of a fin along it. */
  const edgeEnds = (e) => {
    if (e.ends) return e.ends;
    let start = null, end = null;
    const seen = new Set();
    for (let f = get(e.fin); f && !seen.has(f); f = get(f.other)) {
      seen.add(f);
      if (f.type !== 17) throw new Unreadable('EDGE.fin does not lead to fins');
      const v = get(f.vertex);
      if (v && v.type !== 18) throw new Unreadable('FIN.vertex points to another kind of node');
      if (f.sense === '+') end = end || v; else start = start || v;
    }
    return (e.ends = { start, end, fins: [...seen] });
  };

  const bodies = [...N.values()].filter((q) => q.type === 12);
  const facesOf = new Map(bodies.map((b) => [b, []]));
  const edgesOf = new Map(bodies.map((b) => [b, []]));
  for (const q of N.values()) {
    if (q.type === 14) {
      const sh = get(q.shell) || get(q.front_shell), rg = sh && get(sh.region), b = rg && get(rg.body);
      if (b && facesOf.has(b)) facesOf.get(b).push(q); else W.add('A face that belongs to no body was left out');
    } else if (q.type === 16) {
      let o = get(q.owner);
      if (o && o.type === 13) { const rg = get(o.region); o = rg && get(rg.body); }
      if (o && edgesOf.has(o)) edgesOf.get(o).push(q);
    }
  }

  const results = new Map();
  let facesRead = 0, faceTotal = 0, chord = 0, chartEdges = 0, wireEdges = 0;
  bodies.forEach((body, bi) => {
    const bodyName = bodies.length > 1 ? `${name} body ${bi + 1}` : name;
    const faces = facesOf.get(body), edges = edgesOf.get(body);
    faceTotal += faces.length;
    try {
      const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
      const see = (p, r = 0) => { for (let k = 0; k < 3; k++) { if (p[k] - r < lo[k]) lo[k] = p[k] - r; if (p[k] + r > hi[k]) hi[k] = p[k] + r; } };
      for (const e of edges) {
        const { start, end } = edgeEnds(e);
        if (start) see(pointOfVertex(start));
        if (end) see(pointOfVertex(end));
        const cv = get(e.curve);
        if (cv && !start && !end && (cv.type === 31 || cv.type === 32)) see(V(cv, 'centre'), cv.type === 31 ? F(cv, 'radius') : F(cv, 'major_radius'));
      }
      for (const f of faces) { const sn = get(f.surface); if (sn && sn.type === 53) see(V(sn, 'centre'), Math.abs(F(sn, 'radius'))); else if (sn && sn.type === 54) see(V(sn, 'centre'), Math.abs(F(sn, 'major_radius')) + Math.abs(F(sn, 'minor_radius'))); }
      const diag = lo[0] === Infinity ? 0 : Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
      if (!(diag > 0)) { if (faces.length || edges.length) W.add('A body without any extent was left out'); return; }
      const first = opts.deflection !== undefined ? opts.deflection : 0.01 * diag;
      const tol = diag * 1e-9;

      const attempt = (defl) => {
        const loc = { W: new Warnings(), hist: {}, skips: {}, read: 0, charts: 0, crossing: 0, skip(why) { tally(this.skips, why); } };
        for (const e of edges) delete e.pts;
        for (const q of N.values()) if (q.surf !== undefined) delete q.surf;
        const stepFor = stepper(defl, angStep);
        // a vertex must lie this close to the curve of its edge (tolerant vertices: their own tolerance)
        const meet = Math.max(1e-5 * diag, 2 * defl * 1e-2);

        const edgePts = (e) => {
          if (e.pts) return e.pts;
          const { start, end, fins } = edgeEnds(e);
          const Sp = pointOfVertex(start), Ep = pointOfVertex(end);
          let cv = get(e.curve);
          // a tolerant edge has no curve of its own: the curve of its first fin that has one is used for all faces
          if (!cv) for (const f of fins) { const c = get(f.curve); if (c) { cv = c; if (f.sense === '-') cv = { ...c, sense: c.sense === '-' ? '+' : '-', poly: undefined, flipped: true }; break; } }
          const slack = Math.max(meet, (start && start.tolerance) || 0, (end && end.tolerance) || 0, e.tolerance || 0);
          const chordPts = (why) => { if (why) loc.W.add(why); return Sp && Ep && Sp !== Ep ? [Sp, Ep] : Sp ? [Sp] : []; };
          if (!cv) return (e.pts = chordPts(Sp && Ep && dist3(Sp, Ep) > tol ? 'An edge without a curve was replaced by a straight chord' : ''));
          if (!CURVE_TYPES.has(cv.type)) throw new Unreadable(`EDGE.curve points to a ${cv.name} node`);
          if (cv.type === 30) {
            if (!Sp || !Ep) throw new Unreadable('an edge on a straight line has no vertices');
            return (e.pts = [Sp, Ep]);
          }
          let trimmed = null;
          if (cv.type === 133) { const t = polylineOf(cv, stepFor); if (t.why) return (e.pts = chordPts(`An edge on ${t.why} was replaced by a straight chord`)); trimmed = t; }
          const base = trimmed ? trimmed.basis : cv;
          if (base.type === 30) { if (!Sp || !Ep) throw new Unreadable('an edge on a straight line has no vertices'); return (e.pts = [Sp, Ep]); }
          const an = analyticOf(base);
          if (an) {
            const same = !Sp || !Ep || Sp === Ep || dist3(Sp, Ep) <= tol;
            if (Sp && an.off(Sp) > slack + 1e-6 * an.a) return (e.pts = chordPts('An edge whose vertex does not lie on its circle or ellipse was replaced by a straight chord'));
            const a0 = Sp ? an.ang(Sp) : 0;
            let sweep = an.dir * TAU;
            if (!same) { sweep = an.ang(Ep) - a0; if (an.dir > 0) { while (sweep <= 1e-12) sweep += TAU; } else { while (sweep >= -1e-12) sweep -= TAU; } }
            const k = Math.max(same ? 3 : 1, Math.ceil(Math.abs(sweep) / stepFor(Math.max(an.a, an.b)) - 1e-9));
            const pts = [];
            for (let i = 0; i <= k; i++) {
              if (i === 0 && Sp) pts.push(Sp);
              else if (i === k && (same ? Sp : Ep)) pts.push(same ? Sp : Ep);
              else pts.push(an.at(a0 + (sweep * i) / k));
            }
            if (same && !Sp) pts[k] = pts[0];
            return (e.pts = pts);
          }
          const pl = trimmed ? trimmed.inner : polylineOf(cv, stepFor);
          if (pl.why) return (e.pts = chordPts(`An edge on ${pl.why} was replaced by a straight chord`));
          if (pl.note === 'chart') loc.charts++;
          const D = pl.pts;
          if (!Sp || !Ep) { const ring = D.slice(); if (pl.closed) ring[ring.length - 1] = ring[0]; return (e.pts = ring); }
          const a = locate(D, Sp), b = locate(D, Ep);
          // sampled curves are met within the sag of their own facets
          let sag = 0;
          for (let i = 0; i + 1 < D.length; i++) sag = Math.max(sag, dist3(D[i], D[i + 1]));
          if (a.d > slack + 0.05 * sag || b.d > slack + 0.05 * sag) return (e.pts = chordPts('An edge whose vertices do not lie on its curve was replaced by a straight chord'));
          const same = Sp === Ep || dist3(Sp, Ep) <= tol;
          let mid;
          if (!same && b.s > a.s) mid = D.slice(Math.floor(a.s) + 1, Math.ceil(b.s));
          else if (pl.closed) mid = D.slice(Math.floor(a.s) + 1, D.length - 1).concat(D.slice(0, Math.ceil(b.s)));
          else return (e.pts = chordPts('An edge whose vertices lie in the opposite order on its curve was replaced by a straight chord'));
          const keep = mid.filter((q) => dist3(q, Sp) > tol && dist3(q, Ep) > tol);
          return (e.pts = [Sp, ...keep, same ? Sp : Ep]);
        };

        const used = new Map();
        const loopsOf = (face) => {
          const loops = [];
          const seenLoops = new Set();
          for (let lp = get(face.loop); lp && !seenLoops.has(lp); lp = get(lp.next)) {
            if (lp.type !== 15) throw new Unreadable(`FACE.loop leads to a ${lp.name} node`);
            seenLoops.add(lp);
            const pts = [];
            const seen = new Set();
            for (let f = get(lp.fin); f && !seen.has(f); f = get(f.forward)) {
              if (f.type !== 17) throw new Unreadable(`LOOP.fin leads to a ${f.name} node`);
              seen.add(f);
              const e = get(f.edge);
              if (!e) { const p = pointOfVertex(get(f.vertex)); if (p && (!pts.length || pts[pts.length - 1] !== p)) pts.push(p); continue; }
              if (e.type !== 16) throw new Unreadable(`FIN.edge points to a ${e.name} node`);
              used.set(e, (used.get(e) || 0) + 1);
              const ep = edgePts(e);
              const seq = f.sense === '-' ? ep.slice().reverse() : ep;
              for (let i = 0; i + 1 < seq.length; i++) if (!pts.length || pts[pts.length - 1] !== seq[i]) pts.push(seq[i]);
              if (seq.length === 1 && (!pts.length || pts[pts.length - 1] !== seq[0])) pts.push(seq[0]);
            }
            while (pts.length > 1 && pts[0] === pts[pts.length - 1]) pts.pop();
            if (pts.length) loops.push(pts);
          }
          return loops;
        };

        /** brep.js surface of a surface node (or of both readings of a B-surface net), or { why }. */
        const surfaceOf = (sn, loops) => {
          if (sn.surf !== undefined) return sn.surf;
          let s = null, why = `${sn.name} surface, which is not read`;
          switch (sn.type) {
            case 50: s = planeSurface(V(sn, 'pvec'), V3.norm(V(sn, 'normal')), V(sn, 'x_axis')); break;
            case 51: s = coneSurface(V(sn, 'pvec'), V3.norm(V(sn, 'axis')), V(sn, 'x_axis'), F(sn, 'radius'), 1, 0, false, stepFor); break;
            case 52: {
              const sin = F(sn, 'sin_half_angle'), cos = F(sn, 'cos_half_angle');
              if (!(Math.abs(cos) > 1e-12) || !(F(sn, 'radius') >= 0)) { why = 'CONE with unusable values'; break; }
              // radius r + v·tan a at the distance v along the axis from pvec; ∂u × ∂v points away from the axis
              const r0 = F(sn, 'radius'), A = V3.norm(V(sn, 'axis'));
              s = r0 > 0 ? coneSurface(V(sn, 'pvec'), A, V(sn, 'x_axis'), r0, 1, sin / cos, false, stepFor) : sin / cos > 0 ? coneSurface(V3.add(V(sn, 'pvec'), A), A, V(sn, 'x_axis'), sin / cos, 1, sin / cos, false, stepFor) : null;
              if (!s) why = 'CONE with unusable values';
              break;
            }
            case 53: s = sphereSurface(V(sn, 'centre'), F(sn, 'radius'), V3.norm(V(sn, 'axis')), V(sn, 'x_axis'), false, stepFor); break;
            case 54: s = F(sn, 'major_radius') > 0 ? torusSurface(V(sn, 'centre'), V3.norm(V(sn, 'axis')), F(sn, 'major_radius'), F(sn, 'minor_radius'), V(sn, 'x_axis'), false, stepFor) : null; if (!s) why = 'TORUS of the "lemon" kind (negative major radius), which is not read'; break;
            case 124: {
              const net = nurbsSurfOf(sn);
              s = net ? splineSurface(net, defl, angStep) : null;
              if (!s) why = 'B-surface whose data are inconsistent or collapsed';
              break;
            }
            case 67: {
              // S(u, v) = C(u) + v·sweep
              const sec = get(sn.section), sw = V(sn, 'sweep');
              if (sec && sec.type === 30) { const d = V3.norm(V(sec, 'direction')); const nrm = V3.cross(sec.sense === '-' ? V3.scale(d, -1) : d, V3.norm(sw)); if (V3.len(nrm) > 1e-9) s = planeSurface(V(sec, 'pvec'), V3.norm(nrm), d); break; }
              const bs = sec && sec.type === 134 ? bcurveOf(sec) : null;
              if (!bs || bs.geo !== 3) { why = sec && sec.type === 134 ? 'SWEPT_SURF whose section B-curve is inconsistent' : `SWEPT_SURF over a ${sec ? sec.name : 'missing'} curve, which is not read`; break; }
              // the sweep is unbounded: take the stretch the face boundary needs
              const dir = V3.norm(sw);
              let c0 = Infinity, c1 = -Infinity, p0 = Infinity, p1 = -Infinity;
              for (let i = 0; i < bs.P.length; i += 3) { const h = bs.P[i] * dir[0] + bs.P[i + 1] * dir[1] + bs.P[i + 2] * dir[2]; if (h < c0) c0 = h; if (h > c1) c1 = h; }
              for (const L of loops) for (const p of L) { const h = V3.dot(p, dir); if (h < p0) p0 = h; if (h > p1) p1 = h; }
              if (!(p1 >= p0)) { why = 'SWEPT_SURF face without a boundary'; break; }
              const v0 = p0 - c1 - 1e-3 * diag, v1 = p1 - c0 + 1e-3 * diag, m = bs.P.length / 3;
              const P = new Float64Array(2 * bs.P.length), Wt = bs.W ? new Float64Array(2 * m) : null;
              for (let j = 0; j < 2; j++) for (let i = 0; i < m; i++) { const o = i + m * j, h = j ? v1 : v0; P[3 * o] = bs.P[3 * i] + h * dir[0]; P[3 * o + 1] = bs.P[3 * i + 1] + h * dir[1]; P[3 * o + 2] = bs.P[3 * i + 2] + h * dir[2]; if (Wt) Wt[o] = bs.W[i]; }
              s = splineSurface({ pu: bs.deg, pv: 1, U: bs.U, V: [v0, v0, v1, v1], nu: m, nv: 2, P, W: Wt }, defl, angStep);
              if (s && sec.sense === '-') s = { ...s, inward: !s.inward };
              if (!s) why = 'SWEPT_SURF whose section is inconsistent';
              return s ? { s, once: true } : { why };
            }
            default: if (!SURFACE_TYPES.has(sn.type)) throw new Unreadable(`FACE.surface points to a ${sn.name} node`);
          }
          return (sn.surf = s ? { s } : { why });
        };

        const mesh = new BodyMesh(diag);
        let faceNo = 0;
        for (const face of faces) {
          const label = `face ${++faceNo}`;
          const sn = get(face.surface);
          if (!sn) { loc.skip('face without a surface'); continue; }
          const loops = loopsOf(face);
          const so = surfaceOf(sn, loops);
          if (so.why) { loc.skip(so.why); continue; }
          const surf = so.s;
          // aligned when the senses of face and surface agree; `inward`: the natural normal opposes ∂u × ∂v of brep.js
          const flip = surf.inward !== (face.sense !== sn.sense);
          let res;
          try { res = tessellateFace(surf, loops, flip, Math.max(2 * defl, 1e-3 * diag)); } catch (e) { res = `triangulation failed (${e.message})`; }
          if (typeof res === 'string') { if (res === CROSSING) loc.crossing++; loc.skip(`${surf.type === 'spline' ? sn.name : surf.type}: ${res}`); continue; }
          if (!mesh.addFace(res.pts, res.indices, flip, label)) { loc.skip(`${surf.type}: its boundary is degenerate`); continue; }
          loc.read++;
          tally(loc.hist, sn.type === 124 ? 'B-surface' : sn.type === 67 ? 'swept surface' : surf.type);
        }
        return { loc, mesh, used, edgePts };
      };

      let defl = first, pass = attempt(defl);
      for (let k = 0; k < 3 && pass.loc.crossing; k++) { defl /= 4; pass = attempt(defl); }
      if (defl !== first) extra.push(`${bodyName}: tessellated with a chord error of ${defl.toPrecision(3)} m instead of ${first.toPrecision(3)} m, because the coarser facets of neighbouring boundary loops crossed.`);
      chord = Math.max(chord, defl);
      const { loc, mesh, used, edgePts } = pass;
      const out = { name: cleanName(bodyName), mesh: null, wires: [], faces: faces.length, facesRead: loc.read, freeEdges: 0, nonManifoldEdges: 0 };
      if (!mesh.empty) {
        out.mesh = mesh.finish(null);
        for (const [, n] of used) { if (n === 1) out.freeEdges++; else if (n > 2) out.nonManifoldEdges++; }
        if (out.mesh.reoriented) extra.push(`${bodyName}: the faces enclose a negative volume as stored, so their orientation was reversed to point outwards.`);
      }
      // wireframe edges
      for (const e of edges) {
        if (used.has(e)) continue;
        const pts = edgePts(e);
        if (pts.length > 1) { wireEdges++; for (let i = 0; i + 1 < pts.length; i++) out.wires.push([pts[i], pts[i + 1]]); }
      }
      facesRead += loc.read; chartEdges += loc.charts;
      for (const k in loc.hist) tally(hist, k, loc.hist[k]);
      for (const k in loc.skips) tally(faceSkips, k, loc.skips[k]);
      for (const [t, n] of loc.W.m) for (let q = 0; q < n; q++) W.add(t);
      results.set(body, out);
    } catch (e) {
      if (!(e instanceof Unreadable)) throw e;
      tally(skipped, 'body that could not be interpreted');
      extra.push(`${bodyName} was not read: ${e.message}. Nothing of this body is shown, because its data could not be interpreted with certainty.`);
    }
  });

  // ---- placements: bodies as they are, or once per instance of the assemblies that hold them ----
  const ident = { R: [1, 0, 0, 0, 1, 0, 0, 0, 1], t: [0, 0, 0], s: 1 };
  const compose = (a, b) => {
    // apply b first, then a:  x' = (Ra·((Rb·x + tb)·sb) + ta)·sa = (Ra·Rb·x + Ra·tb + ta / sb)·sa·sb
    const R = new Array(9);
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) R[3 * i + j] = a.R[3 * i] * b.R[j] + a.R[3 * i + 1] * b.R[3 + j] + a.R[3 * i + 2] * b.R[6 + j];
    return { R, t: [0, 1, 2].map((i) => a.R[3 * i] * b.t[0] + a.R[3 * i + 1] * b.t[1] + a.R[3 * i + 2] * b.t[2] + a.t[i] / b.s), s: a.s * b.s };
  };
  const placements = [];
  const instanced = new Set();
  const assemblies = [...N.values()].filter((q) => q.type === 10);
  if (assemblies.length) {
    const child = new Set();
    for (const q of N.values()) if (q.type === 11) { const pt = get(q.part); if (pt) child.add(pt); }
    const walk = (asm, M, depth) => {
      if (depth > 32 || placements.length > 20000) return;
      const seen = new Set();
      for (let ins = get(asm.sub_instance); ins && !seen.has(ins); ins = get(ins.next_in_part)) {
        seen.add(ins);
        if (ins.type !== 11) break;
        const part = get(ins.part), tr = get(ins.transform);
        let T = ident;
        if (tr && tr.type === 100 && Array.isArray(tr.rotation_matrix) && tr.rotation_matrix.every((x) => typeof x === 'number')) T = { R: tr.rotation_matrix, t: vecOk(tr.translation_vector) ? tr.translation_vector : [0, 0, 0], s: typeof tr.scale === 'number' && tr.scale > 0 ? tr.scale : 1 };
        const here = compose(M, T);
        if (part && part.type === 12) { placements.push([part, here]); instanced.add(part); } else if (part && part.type === 10) walk(part, here, depth + 1);
      }
    };
    for (const a of assemblies) if (!child.has(a)) walk(a, ident, 0);
  }
  for (const b of bodies) if (!instanced.has(b)) placements.push([b, ident]);

  const models = [];
  const frameSegs = [];
  const perBody = new Map();
  let closedBodies = 0, orientedBodies = 0, solidBodies = 0;
  const volumes = [];
  for (const [body, M] of placements) {
    const r = results.get(body);
    if (!r) continue;
    const k = (perBody.get(body) || 0) + 1;
    perBody.set(body, k);
    const nm = k > 1 ? `${r.name} (${k})` : r.name;
    const put = (p) => { const x = M.R[0] * p[0] + M.R[1] * p[1] + M.R[2] * p[2] + M.t[0], y = M.R[3] * p[0] + M.R[4] * p[1] + M.R[5] * p[2] + M.t[1], z = M.R[6] * p[0] + M.R[7] * p[1] + M.R[8] * p[2] + M.t[2]; return [x * M.s, y * M.s, z * M.s]; };
    if (r.mesh) {
      const c = new Float64Array(r.mesh.coords.length);
      for (let i = 0; i < c.length; i += 3) { const q = put([r.mesh.coords[i], r.mesh.coords[i + 1], r.mesh.coords[i + 2]]); c[i] = q[0]; c[i + 1] = q[1]; c[i + 2] = q[2]; }
      const det = M.R[0] * (M.R[4] * M.R[8] - M.R[5] * M.R[7]) - M.R[1] * (M.R[3] * M.R[8] - M.R[5] * M.R[6]) + M.R[2] * (M.R[3] * M.R[7] - M.R[4] * M.R[6]);
      let idx = r.mesh.indices;
      if (det < 0) { idx = Uint32Array.from(idx); for (let i = 0; i < idx.length; i += 3) { const x = idx[i + 1]; idx[i + 1] = idx[i + 2]; idx[i + 2] = x; } }
      models.push(surfaceModel(nm, c, idx, r.mesh.groups.map((g) => ({ ...g }))));
      if (r.mesh.closed) closedBodies++;
      if (r.mesh.closed && r.mesh.oriented) orientedBodies++;
      if (!r.freeEdges && !r.nonManifoldEdges) solidBodies++;
      volumes.push(r.mesh.volume === null ? null : r.mesh.volume * M.s ** 3);
    }
    for (const [a, b] of r.wires) frameSegs.push([put(a), put(b), nm]);
  }
  if (frameSegs.length) {
    let lo = Infinity, hi = -Infinity;
    for (const sg of frameSegs) for (const p of [sg[0], sg[1]]) for (const v of p) { if (v < lo) lo = v; if (v > hi) hi = v; }
    const weld = new Welder(Math.max(hi - lo, 1e-12) * 1e-9, frameSegs.length * 2);
    const members = [];
    for (const [a, b, g] of frameSegs) { const i = weld.add(a[0], a[1], a[2]), j = weld.add(b[0], b[1], b[2]); if (i !== j) members.push({ a: i, b: j, group: g }); }
    if (members.length) { models.push({ kind: 'frame', name: cleanName(models.length ? `${name} (wires)` : name), nodes: weld.finish(), members }); tally(hist, 'wire edge', wireEdges); }
  }

  const warnings = [];
  const nSkip = Object.values(faceSkips).reduce((a, b) => a + b, 0);
  if (nSkip) warnings.push(`${nSkip} of ${faceTotal} face${faceTotal > 1 ? 's' : ''} could not be read (${Object.entries(faceSkips).map(([k, v]) => `${v} × ${k}`).join('; ')}) — export those bodies as STEP to get them complete.`);
  if (chartEdges) warnings.push(`${chartEdges} edge${chartEdges > 1 ? 's follow' : ' follows'} the chart of points stored for an intersection curve, joined by straight segments.`);
  for (const w of W.list()) warnings.push(`${w}.`);
  warnings.push(...extra);
  if (!models.length) throw new Error(`${fname}: no geometry could be read from the Parasolid file.${warnings.length ? ` ${warnings.join(' ')}` : ''}`);
  for (const k in faceSkips) skipped[`face: ${k}`] = faceSkips[k];
  return finish(models, warnings, {
    format: 'parasolid', units: { unit: 'm', scale: 1, basis: 'file' }, crs: null, elementTypes: hist, skipped,
    schema: xt.header.schema, schemaLayout: xt.layout, modellerVersion: xt.header.modeller, application: xt.header.application, saved: xt.header.date,
    bodies: bodies.length, instances: placements.length, solidBodies, closedBodies, orientedBodies,
    faces: facesRead, faceTotal, facesSkipped: nSkip, volumes, nodeTypes: xt.types,
    tessellation: { chord, angularStepDeg: opts.angularStep ?? 15 },
  });
}
