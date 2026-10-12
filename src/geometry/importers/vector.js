// Vector geospatial data: the feature collector shared by every vector reader (GeoJSON, WKT, Shapefile, GML,
// KML, GPX, LandXML) and the GeoJSON and WKT readers themselves.
//
// Mapping to model kinds:
//   points            → one `points` model; feature properties become per-point attributes
//   lines             → one `polyline` per line string, `properties` attached
//   polygons, flat    → one `section` per polygon (outer loop first, then holes), `properties` attached
//   polygons, 3-D     → closed `polyline` per ring (pass opts.polygons = 'polyline' to force this)
//   TIN / polyhedral  → `surface`
//
// Coordinates in degrees are projected to metres with units.geoProjector (local tangent plane about
// `opts.origin` or the data centre by default; `opts.projection = 'utm' | 'none'`).

import { inputText, baseName, Grow, finish, pointsModel, surfaceModel, cleanName, tally, triangulate3D, closeLoop } from './common.js';
import { geoProjector, looksGeographic, epsgInfo, resolveUnits } from '../units.js';

/** Interpret a CRS identifier ("EPSG:32631", "urn:ogc:def:crs:EPSG::4326", "…/EPSG/0/25832", "CRS84"). */
export function crsFromName(s) {
  if (!s) return null;
  const str = String(s);
  if (/CRS84|CRS:84/i.test(str)) return { epsg: 4326, name: 'WGS 84 (CRS84)', kind: 'geographic', unit: 'deg', latFirst: false, source: 'EPSG:4326' };
  const m = /EPSG(?::+|\/\d*\.?\d*\/|:\d*\.?\d*:)(\d{4,6})/i.exec(str) || /^(\d{4,6})$/.exec(str.trim());
  if (!m) return { epsg: null, name: str, kind: 'unknown', unit: null, source: str };
  const code = Number(m[1]);
  const info = epsgInfo(code);
  // the URN / URL spellings carry the authority's axis order (latitude first for geographic systems)
  const urn = /^urn:|^https?:/i.test(str);
  if (info) return { ...info, latFirst: info.kind === 'geographic' && urn, source: `EPSG:${code}` };
  return { epsg: code, name: `EPSG:${code}`, kind: 'unknown', unit: null, source: `EPSG:${code}` };
}

/** Collects point / line / polygon / mesh features and builds the models. */
export class Features {
  constructor(name) {
    this.name = name;
    this.pts = new Grow(Float64Array, 3 * 256);
    this.ptProps = [];
    this.lines = [];
    this.polys = [];
    this.meshes = [];
    this.hist = {};
    this.skipped = {};
    this.b = [Infinity, Infinity, -Infinity, -Infinity];
  }
  see(x, y) { const b = this.b; if (x < b[0]) b[0] = x; if (x > b[2]) b[2] = x; if (y < b[1]) b[1] = y; if (y > b[3]) b[3] = y; }
  addPoint(x, y, z, props, name) {
    if (!(x === x && y === y)) { tally(this.skipped, 'pointWithoutCoordinates'); return; }
    this.pts.push3(x, y, z === z && z !== undefined ? z : 0);
    this.ptProps.push(name !== undefined && name !== null ? { ...(props || {}), name } : props || null);
    this.see(x, y);
    tally(this.hist, 'point');
  }
  /** @param {ArrayLike<number>} c flat x, y, z */
  addLine(c, props, name) {
    if (c.length < 6) { tally(this.skipped, 'lineWithFewerThanTwoPoints'); return; }
    for (let i = 0; i < c.length; i += 3) this.see(c[i], c[i + 1]);
    this.lines.push({ c: c instanceof Float64Array ? c : Float64Array.from(c), props: props || null, name });
    tally(this.hist, 'lineString');
  }
  /** @param {ArrayLike<number>[]} rings flat x, y, z per ring, outer ring first */
  addPolygon(rings, props, name) {
    const ok = rings.filter((r) => r.length >= 9);
    if (!ok.length || rings[0].length < 9) { tally(this.skipped, 'polygonWithFewerThanThreePoints'); return; }
    for (const r of ok) for (let i = 0; i < r.length; i += 3) this.see(r[i], r[i + 1]);
    this.polys.push({ rings: ok.map((r) => (r instanceof Float64Array ? r : Float64Array.from(r))), props: props || null, name });
    tally(this.hist, 'polygon');
  }
  addMesh(coords, indices, name, groups) {
    for (let i = 0; i < coords.length; i += 3) this.see(coords[i], coords[i + 1]);
    this.meshes.push({ coords, indices, name, groups });
    tally(this.hist, 'tri3', indices.length / 3);
  }
  get empty() { return !this.pts.n && !this.lines.length && !this.polys.length && !this.meshes.length; }

  /**
   * @param {object} opts importer options
   * @param {object} o `{ format, crs?: crsFromName()/parsePrj() style info, geographic?: boolean, unitHint?,
   *   warnings?, extra?, expected }`. `geographic` forces the interpretation; otherwise the CRS decides, and
   *   without a CRS the coordinates are taken as projected.
   */
  build(opts, o) {
    const warnings = o.warnings || [];
    if (this.empty) throw new Error(`${o.fileName || this.name}: no geometry was found (${o.expected}).`);
    let crs = o.crs || null;
    if (opts.crs) { const c = /^geographic$/i.test(String(opts.crs)) ? crsFromName('EPSG:4326') : crsFromName(opts.crs); if (c) crs = { ...c, basis: 'option' }; }
    const b = this.b;
    let geographic = o.geographic;
    if (crs && crs.kind === 'geographic') geographic = true; else if (crs && crs.kind === 'projected') geographic = false;
    if (geographic === undefined || geographic === null) geographic = false;
    if (geographic && !looksGeographic(b[0], b[2], b[1], b[3])) {
      if (crs && crs.kind === 'geographic' && crs.basis !== 'default') throw new Error(`${o.fileName || this.name}: the coordinate system is geographic (${crs.source || crs.name}) but the coordinates (${b[0]} … ${b[2]}, ${b[1]} … ${b[3]}) are not longitude/latitude in degrees. If they are projected, re-import with crs set to their EPSG code.`);
      geographic = false;
      warnings.push('The coordinates are outside ±180° / ±90°, so they were read as projected metres rather than longitude/latitude.');
    }
    let units, fwd, crsMeta;
    if (geographic) {
      const P = geoProjector(opts, (b[0] + b[2]) / 2, (b[1] + b[3]) / 2);
      fwd = P.forward;
      units = { unit: 'deg', scale: null, basis: 'file' };
      crsMeta = { source: (crs && crs.source) || 'EPSG:4326', kind: 'geographic', epsg: (crs && crs.epsg) || 4326, name: (crs && crs.name) || 'WGS 84', basis: (crs && crs.basis) || 'file', ...P.crs };
    } else {
      units = crs && crs.unitScale && crs.unitScale !== 1 && !opts.units ? { unit: crs.unit || 'custom', scale: crs.unitScale, basis: 'file' } : resolveUnits(opts, o.unitHint || (crs && crs.kind === 'projected' ? 'm' : null), 'm');
      const s = units.scale;
      fwd = s === 1 ? null : (x, y) => [x * s, y * s];
      crsMeta = crs ? { source: crs.source || (crs.epsg ? `EPSG:${crs.epsg}` : crs.name), kind: crs.kind, epsg: crs.epsg ?? null, name: crs.name || null, basis: crs.basis || 'file', units: 'm' } : null;
      if (!crs && looksGeographic(b[0], b[2], b[1], b[3]) && o.warnIfDegrees) warnings.push('The file names no coordinate system and its coordinates would fit longitude/latitude in degrees; they were read as projected metres. Re-import with crs = "EPSG:4326" if they are geographic.');
    }
    if (opts.vertical && opts.vertical !== 'depth' && opts.vertical !== 'elevation') throw new RangeError(`vertical: opts.vertical "${opts.vertical}" is not one of elevation, depth.`);
    const zs = (opts.vertical === 'depth' ? -1 : 1) * (geographic ? 1 : (o.zScale ?? units.scale ?? 1));
    const conv = (src) => {
      const out = new Float64Array(src.length);
      for (let i = 0; i < src.length; i += 3) {
        if (fwd) { const p = fwd(src[i], src[i + 1]); out[i] = p[0]; out[i + 1] = p[1]; } else { out[i] = src[i]; out[i + 1] = src[i + 1]; }
        out[i + 2] = src[i + 2] * zs;
      }
      return out;
    };
    const models = [];
    for (const m of this.meshes) models.push(surfaceModel(cleanName(m.name, this.name), conv(m.coords), m.indices, m.groups || null));
    if (this.pts.n) {
      const n = this.pts.n / 3;
      const attrs = {};
      const keys = new Set();
      for (const p of this.ptProps) if (p) for (const k in p) keys.add(k);
      for (const k of keys) {
        let numeric = true, any = false;
        for (const p of this.ptProps) { const v = p ? p[k] : undefined; if (v === undefined || v === null || v === '') continue; any = true; if (typeof v !== 'number') { numeric = false; break; } }
        if (!any) continue;
        const key = k === 'name' ? 'names' : k;
        if (numeric) { const a = new Float64Array(n).fill(NaN); this.ptProps.forEach((p, i) => { if (p && typeof p[k] === 'number') a[i] = p[k]; }); attrs[key] = a; } else attrs[key] = this.ptProps.map((p) => (p && p[k] !== undefined && p[k] !== null ? (typeof p[k] === 'object' ? JSON.stringify(p[k]) : String(p[k])) : ''));
      }
      models.push(pointsModel(this.meshes.length || this.lines.length || this.polys.length ? `${this.name} (points)` : this.name, conv(this.pts.view()), attrs));
    }
    const label = (f, i, n, what) => cleanName(f.name ?? (f.props && (f.props.name ?? f.props.Name ?? f.props.NAME ?? f.props.id ?? f.props.ID)), n > 1 ? `${this.name} ${what} ${i + 1}` : this.name);
    this.lines.forEach((l, i) => {
      const m = { kind: 'polyline', name: label(l, i, this.lines.length, 'line'), points: conv(l.c) };
      if (l.props && Object.keys(l.props).length) m.properties = l.props;
      models.push(m);
    });
    let flat = 0, spatial = 0;
    this.polys.forEach((pg, i) => {
      const nm = label(pg, i, this.polys.length, 'polygon');
      const rings = pg.rings.map(conv);
      let zlo = Infinity, zhi = -Infinity;
      for (const r of rings) for (let k = 2; k < r.length; k += 3) { if (r[k] < zlo) zlo = r[k]; if (r[k] > zhi) zhi = r[k]; }
      const isFlat = zhi - zlo <= 1e-9 * Math.max(1, Math.abs(zhi));
      if (isFlat && opts.polygons !== 'polyline') {
        const loops = [];
        for (const r of rings) { const l2 = new Float64Array((r.length / 3) * 2); for (let k = 0; k < r.length / 3; k++) { l2[2 * k] = r[3 * k]; l2[2 * k + 1] = r[3 * k + 1]; } const c = closeLoop(l2); if (c) loops.push(c); }
        if (!loops.length) return;
        const m = { kind: 'section', name: nm, loops };
        if (zlo !== 0 && Number.isFinite(zlo)) m.elevation = zlo;
        if (pg.props && Object.keys(pg.props).length) m.properties = pg.props;
        models.push(m);
        flat++;
      } else {
        rings.forEach((r, k) => {
          const m = { kind: 'polyline', name: k ? `${nm} (hole ${k})` : nm, points: r, closed: true };
          if (pg.props && Object.keys(pg.props).length) m.properties = pg.props;
          models.push(m);
        });
        spatial++;
      }
    });
    if (!models.length) throw new Error(`${o.fileName || this.name}: no usable geometry was found (${o.expected}).`);
    for (const k in this.skipped) warnings.push(`Skipped ${this.skipped[k]} ${k}.`);
    return finish(models, warnings, {
      format: o.format, units, crs: crsMeta, vertical: { source: opts.vertical === 'depth' ? 'depth-positive-down' : 'elevation-up', sign: opts.vertical === 'depth' ? -1 : 1, basis: opts.vertical ? 'option' : 'assumed' },
      elementTypes: this.hist, skipped: this.skipped, features: { points: this.pts.n / 3, lines: this.lines.length, polygons: this.polys.length, polygonsAsSections: flat, polygonsAsPolylines: spatial, meshes: this.meshes.length },
      ...(o.extra || {}),
    });
  }
}

/** Triangulate polygon patches (each `[ring, ring…]` of flat xyz) into one welded-by-patch mesh. */
export function patchesToMesh(patches) {
  const coords = new Grow(Float64Array, 3072), idx = new Grow(Uint32Array, 3072);
  for (const rings of patches) {
    const total = rings.reduce((s, r) => s + r.length, 0);
    const flat = new Float64Array(total), holes = [];
    let o = 0;
    rings.forEach((r, k) => { if (k) holes.push(o / 3); flat.set(r, o); o += r.length; });
    const base = coords.n / 3;
    const t = triangulate3D(flat, holes.length ? holes : undefined);
    if (!t.length) continue;
    coords.ensure(total);
    coords.a.set(flat, coords.n); coords.n += total;
    for (const k of t) idx.push(base + k);
  }
  return { coords: coords.finish(), indices: idx.finish() };
}

const dropClosing = (c) => { const n = c.length; return n >= 6 && c[0] === c[n - 3] && c[1] === c[n - 2] && c[2] === c[n - 1] ? c.subarray(0, n - 3) : c; };

// ---------------------------------------------------------------------------------------------------------
// GeoJSON
// ---------------------------------------------------------------------------------------------------------

function flatCoords(list) {
  const out = new Float64Array(list.length * 3);
  for (let i = 0; i < list.length; i++) { const p = list[i]; out[3 * i] = +p[0]; out[3 * i + 1] = +p[1]; out[3 * i + 2] = p.length > 2 && Number.isFinite(+p[2]) ? +p[2] : 0; }
  return out;
}

function addGeoJsonGeometry(F, g, props, name) {
  if (!g || !g.type) { tally(F.skipped, 'featureWithoutGeometry'); return; }
  const c = g.coordinates;
  switch (g.type) {
    case 'Point': F.addPoint(+c[0], +c[1], +c[2], props, name); break;
    case 'MultiPoint': for (const p of c) F.addPoint(+p[0], +p[1], +p[2], props, name); break;
    case 'LineString': F.addLine(flatCoords(c), props, name); break;
    case 'MultiLineString': for (const l of c) F.addLine(flatCoords(l), props, name); break;
    case 'Polygon': F.addPolygon(c.map((r) => dropClosing(flatCoords(r))), props, name); break;
    case 'MultiPolygon': for (const pg of c) F.addPolygon(pg.map((r) => dropClosing(flatCoords(r))), props, name); break;
    case 'GeometryCollection': for (const s of g.geometries || []) addGeoJsonGeometry(F, s, props, name); break;
    default: tally(F.skipped, `geometryType:${g.type}`);
  }
}

/** GeoJSON (RFC 7946), including newline-delimited feature files and the legacy `crs` member. */
export function parseGeoJson(input, opts = {}) {
  const text = inputText(input);
  const name = baseName(input.name);
  const F = new Features(name);
  let roots;
  try { roots = [JSON.parse(text)]; } catch (e) {
    roots = [];
    for (const l of text.split(/\r?\n/)) { const t = l.trim().replace(/^\u001e/, ''); if (!t) continue; try { roots.push(JSON.parse(t)); } catch { throw new Error(`${input.name || 'GeoJSON'}: not valid JSON (${e.message}).`); } }
  }
  let crs = null;
  const visit = (o) => {
    if (!o || typeof o !== 'object') return;
    if (o.crs && o.crs.properties && o.crs.properties.name && !crs) crs = crsFromName(o.crs.properties.name);
    if (o.type === 'FeatureCollection') { for (const f of o.features || []) visit(f); } else if (o.type === 'Feature') addGeoJsonGeometry(F, o.geometry, o.properties || null, o.properties ? undefined : o.id);
    else if (o.type === 'Topology') throw new Error(`${input.name || 'GeoJSON'}: TopoJSON is not read. Convert it to GeoJSON first (for example with "topo2geo" or mapshaper).`);
    else if (o.type) addGeoJsonGeometry(F, o, null);
    else throw new Error(`${input.name || 'GeoJSON'}: not GeoJSON (expected an object with "type": "FeatureCollection", "Feature" or a geometry type).`);
  };
  for (const r of roots) visit(r);
  const b = F.b;
  // RFC 7946: coordinates are WGS 84 longitude/latitude unless a (legacy) crs member says otherwise
  const looks = !F.empty && looksGeographic(b[0], b[2], b[1], b[3]);
  const warnings = [];
  if (!crs && !looks && !F.empty) warnings.push('The coordinates are not longitude/latitude in degrees (RFC 7946) and the file has no crs member; they were read as projected metres.');
  if (crs && crs.kind === 'unknown') { warnings.push(`Coordinate system ${crs.source} is not one this application can classify; coordinates were read as projected metres.`); crs.kind = 'projected'; }
  return F.build(opts, { format: 'geojson', fileName: input.name, crs: crs || (looks ? { ...crsFromName('EPSG:4326'), basis: 'default' } : null), warnings, expected: 'expected features with Point, LineString or Polygon geometries' });
}

// ---------------------------------------------------------------------------------------------------------
// WKT
// ---------------------------------------------------------------------------------------------------------

/** Parse one WKT geometry starting at `pos`. Returns `{ type, dim, coords, next }` or null. */
export function parseWktGeometry(s, pos = 0) {
  const re = /\s*([A-Za-z]+)\s*(ZM|Z|M)?\s*/y;
  re.lastIndex = pos;
  const m = re.exec(s);
  if (!m) return null;
  const type = m[1].toUpperCase();
  const flag = (m[2] || '').toUpperCase();
  let i = re.lastIndex;
  if (/^EMPTY/i.test(s.slice(i, i + 5))) return { type, flag, coords: [], next: i + 5 };
  if (s[i] !== '(') return null;
  if (type === 'GEOMETRYCOLLECTION') {
    const list = [];
    i++;
    for (;;) {
      const g = parseWktGeometry(s, i);
      if (!g) break;
      list.push(g);
      i = g.next;
      while (i < s.length && /[\s,]/.test(s[i])) i++;
      if (s[i] === ')') { i++; break; }
    }
    return { type, flag, geometries: list, next: i };
  }
  const nested = () => {
    // at '(': either a list of parenthesised items or a flat run of "x y z, x y z"
    i++;
    while (i < s.length && /\s/.test(s[i])) i++;
    if (s[i] === '(') {
      const out = [];
      for (;;) {
        while (i < s.length && /[\s,]/.test(s[i])) i++;
        if (s[i] === '(') out.push(nested()); else { if (s[i] === ')') i++; return out; }
      }
    }
    const e = s.indexOf(')', i);
    const pts = pointsOf(s.slice(i, e < 0 ? s.length : e));
    i = e < 0 ? s.length : e + 1;
    return pts;
  };
  const pointsOf = (t) => t.split(',').map((p) => p.trim().split(/\s+/).map(Number)).filter((p) => p.length >= 2 && p[0] === p[0]);
  const coords = nested();
  return { type, flag, coords, next: i };
}

function wktFlat(pts, flag) {
  // "M" alone means the third ordinate is a measure, not an elevation
  const out = new Float64Array(pts.length * 3);
  pts.forEach((p, k) => { out[3 * k] = p[0]; out[3 * k + 1] = p[1]; out[3 * k + 2] = flag === 'M' ? 0 : (p[2] ?? 0); });
  return out;
}

function addWkt(F, g, name) {
  const isPts = (c) => Array.isArray(c) && (!c.length || typeof c[0][0] === 'number');
  const fl = (p) => wktFlat(p, g.flag);
  switch (g.type) {
    case 'POINT': if (g.coords.length) { const p = g.coords[0]; F.addPoint(p[0], p[1], g.flag === 'M' ? 0 : p[2], null, name); } break;
    case 'MULTIPOINT': for (const c of isPts(g.coords) ? [g.coords] : g.coords) for (const p of c) F.addPoint(p[0], p[1], g.flag === 'M' ? 0 : p[2], null, name); break;
    case 'LINESTRING': case 'LINEARRING': F.addLine(fl(g.coords), null, name); break;
    case 'MULTILINESTRING': for (const l of g.coords) F.addLine(fl(l), null, name); break;
    case 'POLYGON': case 'TRIANGLE': F.addPolygon(g.coords.map((r) => dropClosing(fl(r))), null, name); break;
    case 'MULTIPOLYGON': for (const pg of g.coords) F.addPolygon(pg.map((r) => dropClosing(fl(r))), null, name); break;
    case 'TIN': case 'POLYHEDRALSURFACE': { const m = patchesToMesh(g.coords.map((pg) => pg.map((r) => dropClosing(fl(r))))); if (m.indices.length) F.addMesh(m.coords, m.indices, name); break; }
    case 'GEOMETRYCOLLECTION': for (const s of g.geometries) addWkt(F, s, name); break;
    default: tally(F.skipped, `geometryType:${g.type}`);
  }
}

/** Well-known text (one geometry per line or separated by semicolons; EWKT "SRID=…;" prefixes honoured). */
export function parseWkt(input, opts = {}) {
  const text = inputText(input);
  const name = baseName(input.name);
  const F = new Features(name);
  let crs = null, pos = 0, n = 0;
  const srid = /SRID\s*=\s*(\d+)\s*;/i.exec(text);
  if (srid) crs = crsFromName(`EPSG:${srid[1]}`);
  const body = text.replace(/SRID\s*=\s*\d+\s*;/gi, ' ');
  while (pos < body.length) {
    while (pos < body.length && /[\s;,]/.test(body[pos])) pos++;
    if (pos >= body.length) break;
    const g = parseWktGeometry(body, pos);
    if (!g) { const e = body.indexOf('\n', pos); tally(F.skipped, 'unparsedLine'); pos = e < 0 ? body.length : e + 1; continue; }
    addWkt(F, g, undefined);
    n++;
    pos = g.next;
  }
  if (!n) throw new Error(`${input.name || 'WKT'}: no WKT geometry found (expected text such as "LINESTRING Z (x y z, x y z)" or "POLYGON ((x y, …))").`);
  if (crs && crs.kind === 'unknown') crs.kind = 'projected';
  return F.build(opts, { format: 'wkt', fileName: input.name, crs, warnIfDegrees: true, expected: 'expected POINT, LINESTRING, POLYGON, MULTI…, TIN or POLYHEDRALSURFACE geometries' });
}
