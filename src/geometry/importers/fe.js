// Shared wrap-up for the finite-element deck readers: unit resolution, model building and the report every
// importer returns (counts, bounding box, element-type histogram, skipped records, set names).

import { finish, unitAdvice } from './common.js';
import { resolveUnits, suggestUnits } from '../units.js';

/**
 * @param {import('./common.js').FeMesh} fe collected mesh
 * @param {string} format registry id
 * @param {object} opts importer options (`units`, `scale`)
 * @param {object} o `{ name, unitHint?, fallback?, warnings?, extra?, expected }` — `expected` is the text of
 *   the error raised when the file held no nodes at all
 */
export function feResult(fe, format, opts, o) {
  if (!fe.nodeCount) throw new Error(`${o.name}: no nodes were found (${o.expected}).`);
  const units = resolveUnits(opts, o.unitHint || null, o.fallback || 'm');
  const b = fe.build({ scale: units.scale, name: fe.name });
  const r = finish(b.models, (o.warnings || []).concat(b.warnings), {
    format, units, crs: null, elementTypes: b.elementTypes, skipped: { ...b.skipped, ...(o.skipped || {}) },
    sourceNodes: fe.nodeCount, nodeSets: [...fe.nodeSets.keys()], elementSets: [...fe.elemSets.keys()], ...(o.extra || {}),
  });
  unitAdvice(units, r.meta.bbox, suggestUnits, r.warnings);
  return r;
}

/** Scale the dimensions of a section record (lengths ×s, areas ×s², second moments ×s⁴). */
export function scaleSection(sec, s) {
  if (!sec || s === 1) return sec;
  const out = { ...sec };
  for (const k of ['D', 't', 'h', 'b', 'tw', 'tf', 'r']) if (typeof out[k] === 'number') out[k] *= s;
  if (typeof out.A === 'number') out.A *= s * s;
  for (const k of ['Iy', 'Iz', 'J']) if (typeof out[k] === 'number') out[k] *= s ** 4;
  return out;
}
