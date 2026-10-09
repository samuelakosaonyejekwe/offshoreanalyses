// Units, axis conventions and coordinate reference handling for imported geometry.
//
// All models leave the importers in metres, z up. This module holds the pieces every importer shares:
// length-unit names and scale factors, y-up → z-up rotation, geographic → metric projection (local tangent
// plane or UTM), EPSG sniffing from ESRI .prj WKT and GeoTIFF GeoKeys, and the vertical sign convention.
//
// What this module does NOT do: datum transformations (no Helmert shifts, no grid shifts, no geoid models).
// Coordinates are projected on the ellipsoid of their own datum; WGS 84 / GRS 80 parameters are used for the
// metric projection unless the WKT supplies a spheroid. Vertical datums are reported, never converted.

/** Metres per unit. US survey foot = 1200/3937 m exactly. */
export const UNIT_SCALE = Object.freeze({
  m: 1, mm: 1e-3, cm: 1e-2, dm: 0.1, km: 1e3, um: 1e-6,
  in: 0.0254, ft: 0.3048, 'us-ft': 1200 / 3937, yd: 0.9144, mi: 1609.344, nmi: 1852, mil: 2.54e-5,
});

const UNIT_ALIASES = {
  m: 'm', meter: 'm', metre: 'm', meters: 'm', metres: 'm', mtr: 'm', mtrs: 'm',
  mm: 'mm', millimeter: 'mm', millimetre: 'mm', millimeters: 'mm', millimetres: 'mm', milli: 'mm',
  cm: 'cm', centimeter: 'cm', centimetre: 'cm', centimeters: 'cm', centimetres: 'cm',
  dm: 'dm', decimeter: 'dm', decimetre: 'dm', decimeters: 'dm', decimetres: 'dm',
  km: 'km', kilometer: 'km', kilometre: 'km', kilometers: 'km', kilometres: 'km',
  um: 'um', micron: 'um', microns: 'um', micrometer: 'um', micrometre: 'um', micrometers: 'um', micrometres: 'um', 'µm': 'um',
  in: 'in', inch: 'in', inches: 'in', '"': 'in',
  ft: 'ft', foot: 'ft', feet: 'ft', "'": 'ft', internationalfoot: 'ft', ifoot: 'ft',
  'us-ft': 'us-ft', usft: 'us-ft', ussurveyfoot: 'us-ft', surveyfoot: 'us-ft', footus: 'us-ft', ftus: 'us-ft', ussurveyfeet: 'us-ft', surveyfeet: 'us-ft',
  yd: 'yd', yard: 'yd', yards: 'yd',
  mi: 'mi', mile: 'mi', miles: 'mi',
  nmi: 'nmi', nauticalmile: 'nmi', nauticalmiles: 'nmi',
  mil: 'mil', mils: 'mil', thou: 'mil',
};

/**
 * Canonical unit key for a free-form unit name ('MILLIMETRE', 'Inches', 'US survey foot', 'ft' …).
 * @param {string} name
 * @returns {string|null} key of UNIT_SCALE, or null when not recognised
 */
export function normaliseUnit(name) {
  if (name === undefined || name === null) return null;
  const raw = String(name).trim();
  if (raw === '"' || raw === "'") return UNIT_ALIASES[raw];
  const k = raw.toLowerCase().replace(/[\s_.]/g, '').replace(/-(?!ft)/g, '');
  if (Object.prototype.hasOwnProperty.call(UNIT_ALIASES, k)) return UNIT_ALIASES[k];
  if (k === 'us-ft') return 'us-ft';
  return null;
}

/**
 * Metres per unit for a unit name.
 * @param {string} name any spelling accepted by normaliseUnit
 * @returns {number} scale to metres
 * @throws {RangeError} for an unknown unit
 */
export function unitScale(name) {
  const k = normaliseUnit(name);
  if (!k) throw new RangeError(`units: unknown length unit "${name}" (use one of ${Object.keys(UNIT_SCALE).join(', ')}).`);
  return UNIT_SCALE[k];
}

/**
 * Decide the length unit of a file. Priority: `opts.units` (user override) → the unit found in the file
 * (`hint`) → `fallback`. The result is what importers report as `meta.units`.
 * @param {object} opts importer options; `opts.units` is a unit name, `opts.scale` an explicit metres-per-unit
 * @param {string|null} hint unit name found in the file, if any
 * @param {string} [fallback='m'] unit assumed when neither is present
 * @returns {{ unit: string, scale: number, basis: 'option'|'file'|'assumed' }}
 */
export function resolveUnits(opts, hint, fallback = 'm') {
  if (opts && Number.isFinite(opts.scale) && opts.scale > 0 && !opts.units) return { unit: 'custom', scale: opts.scale, basis: 'option' };
  if (opts && opts.units) {
    const u = normaliseUnit(opts.units);
    if (!u) throw new RangeError(`units: unknown length unit "${opts.units}" (use one of ${Object.keys(UNIT_SCALE).join(', ')}).`);
    return { unit: u, scale: UNIT_SCALE[u], basis: 'option' };
  }
  const h = normaliseUnit(hint);
  if (h) return { unit: h, scale: UNIT_SCALE[h], basis: 'file' };
  const f = normaliseUnit(fallback) || 'm';
  return { unit: f, scale: UNIT_SCALE[f], basis: 'assumed' };
}

/**
 * Plausibility hint for unit-less files: offshore structures are metres to hundreds of metres across, so an
 * extent in the thousands with no unit declared usually means millimetres. Advisory only — never applied
 * automatically.
 * @param {{min:number[], max:number[]}} bbox in file units
 * @returns {{ unit: string, reason: string }|null}
 */
export function suggestUnits(bbox) {
  const ext = Math.max(bbox.max[0] - bbox.min[0], bbox.max[1] - bbox.min[1], bbox.max[2] - bbox.min[2]);
  if (!(ext > 0)) return null;
  if (ext > 5000) return { unit: 'mm', reason: `largest extent is ${ext.toPrecision(4)} file units; as metres that is implausibly large for a structure, so the file is probably in millimetres` };
  if (ext < 0.05) return { unit: 'km', reason: `largest extent is ${ext.toPrecision(3)} file units; the file may be in kilometres or normalised coordinates` };
  return null;
}

/** Multiply a coordinate array by a scalar in place. */
export function scaleArray(arr, s) {
  if (s !== 1) for (let i = 0; i < arr.length; i++) arr[i] *= s;
  return arr;
}

/** Rotate y-up coordinates to z-up in place: (x, y, z) → (x, −z, y). Right-handedness is preserved. */
export function yUpToZUp(arr) {
  for (let i = 0; i + 2 < arr.length; i += 3) { const y = arr[i + 1]; arr[i + 1] = -arr[i + 2]; arr[i + 2] = y; }
  return arr;
}

/** Inverse of yUpToZUp: (x, y, z) → (x, z, −y). */
export function zUpToYUp(arr) {
  for (let i = 0; i + 2 < arr.length; i += 3) { const y = arr[i + 1]; arr[i + 1] = arr[i + 2]; arr[i + 2] = -y; }
  return arr;
}

// ---------------------------------------------------------------------------------------------------------
// Ellipsoid and projections
// ---------------------------------------------------------------------------------------------------------

/** WGS 84 ellipsoid (a in metres). GRS 80 differs in the 10th significant digit of 1/f and is treated alike. */
export const WGS84 = Object.freeze({ a: 6378137, invF: 298.257223563 });

const D2R = Math.PI / 180;

/** Meridional (M) and prime-vertical (N) radii of curvature at latitude `latDeg`, metres. */
export function radii(latDeg, ell = WGS84) {
  const f = 1 / ell.invF, e2 = f * (2 - f);
  const s = Math.sin(latDeg * D2R);
  const w = Math.sqrt(1 - e2 * s * s);
  return { M: (ell.a * (1 - e2)) / (w * w * w), N: ell.a / w };
}

/**
 * Local tangent-plane metres from longitude/latitude: equirectangular about (lon0, lat0) using the ellipsoid's
 * radii of curvature at lat0:  x = N(φ0)·cos φ0·Δλ (east),  y = M(φ0)·Δφ (north).
 * Distortion grows with distance from the origin: about 0.01 % at 10 km and 1 % at 100 km of north–south
 * extent in mid latitudes — fine for a field layout, not for a basin. Use UTM beyond a few tens of km.
 * @returns {[number, number]} [east, north] in metres
 */
export function lonLatToLocal(lonDeg, latDeg, lon0, lat0, ell = WGS84) {
  const { M, N } = radii(lat0, ell);
  let dl = lonDeg - lon0;
  if (dl > 180) dl -= 360; else if (dl < -180) dl += 360;
  return [N * Math.cos(lat0 * D2R) * dl * D2R, M * (latDeg - lat0) * D2R];
}

/** Inverse of lonLatToLocal. @returns {[number, number]} [lon, lat] in degrees */
export function localToLonLat(x, y, lon0, lat0, ell = WGS84) {
  const { M, N } = radii(lat0, ell);
  return [lon0 + x / (N * Math.cos(lat0 * D2R)) / D2R, lat0 + y / M / D2R];
}

/** UTM zone number (1–60) for a longitude, with the Norway (32V) and Svalbard (31X–37X) exceptions. */
export function utmZone(lonDeg, latDeg = 0) {
  let lon = ((((lonDeg + 180) % 360) + 360) % 360) - 180;
  let zone = Math.floor((lon + 180) / 6) + 1;
  if (zone > 60) zone = 60;
  if (latDeg >= 56 && latDeg < 64 && lon >= 3 && lon < 12) zone = 32;
  if (latDeg >= 72 && latDeg < 84) {
    if (lon >= 0 && lon < 9) zone = 31; else if (lon >= 9 && lon < 21) zone = 33; else if (lon >= 21 && lon < 33) zone = 35; else if (lon >= 33 && lon < 42) zone = 37;
  }
  return zone;
}

function krueger(ell) {
  const f = 1 / ell.invF;
  const n = f / (2 - f);
  const n2 = n * n, n3 = n2 * n, n4 = n2 * n2;
  return {
    e: Math.sqrt(f * (2 - f)),
    A: (ell.a / (1 + n)) * (1 + n2 / 4 + n4 / 64),
    alpha: [n / 2 - (2 * n2) / 3 + (5 * n3) / 16, (13 * n2) / 48 - (3 * n3) / 5, (61 * n3) / 240],
    beta: [n / 2 - (2 * n2) / 3 + (37 * n3) / 96, n2 / 48 + n3 / 15, (17 * n3) / 480],
  };
}

/**
 * Transverse Mercator forward (Krüger n-series to third order, as used for UTM; error below a millimetre
 * within a zone and below a centimetre to ±12° from the central meridian).
 * Reference: L. Krüger (1912); the series coefficients are those tabulated by Karney, J. Geodesy 85 (2011), eq. 35.
 * @param {number} latDeg @param {number} lonDeg
 * @param {object} p `{ lon0 (deg), k0, fe (false easting, m), fn (false northing, m), lat0?: 0 }`
 * @returns {{ easting: number, northing: number }} metres
 */
export function tmForward(latDeg, lonDeg, p, ell = WGS84) {
  const K = krueger(ell);
  const phi = latDeg * D2R;
  let dl = lonDeg - p.lon0;
  if (dl > 180) dl -= 360; else if (dl < -180) dl += 360;
  const lam = dl * D2R;
  const s = Math.sin(phi);
  const t = Math.sinh(Math.atanh(s) - K.e * Math.atanh(K.e * s));
  const xi0 = Math.atan2(t, Math.cos(lam));
  const eta0 = Math.atanh(Math.sin(lam) / Math.sqrt(1 + t * t));
  let xi = xi0, eta = eta0;
  for (let j = 1; j <= 3; j++) {
    xi += K.alpha[j - 1] * Math.sin(2 * j * xi0) * Math.cosh(2 * j * eta0);
    eta += K.alpha[j - 1] * Math.cos(2 * j * xi0) * Math.sinh(2 * j * eta0);
  }
  const k0 = p.k0 ?? 0.9996;
  return { easting: (p.fe ?? 500000) + k0 * K.A * eta, northing: (p.fn ?? 0) + k0 * K.A * xi };
}

/**
 * Transverse Mercator inverse (Krüger series, third order).
 * @returns {{ lat: number, lon: number }} degrees
 */
export function tmInverse(easting, northing, p, ell = WGS84) {
  const K = krueger(ell);
  const k0 = p.k0 ?? 0.9996;
  const xi = (northing - (p.fn ?? 0)) / (k0 * K.A);
  const eta = (easting - (p.fe ?? 500000)) / (k0 * K.A);
  let xi0 = xi, eta0 = eta;
  for (let j = 1; j <= 3; j++) {
    xi0 -= K.beta[j - 1] * Math.sin(2 * j * xi) * Math.cosh(2 * j * eta);
    eta0 -= K.beta[j - 1] * Math.cos(2 * j * xi) * Math.sinh(2 * j * eta);
  }
  // Conformal latitude χ, then geodetic latitude by fixed-point iteration on the isometric latitude.
  const chi = Math.asin(Math.sin(xi0) / Math.cosh(eta0));
  const psi = Math.atanh(Math.sin(chi));
  let s = Math.sin(chi);
  for (let it = 0; it < 12; it++) {
    const sNew = Math.tanh(psi + K.e * Math.atanh(K.e * s));
    if (Math.abs(sNew - s) < 1e-15) { s = sNew; break; }
    s = sNew;
  }
  const lam = Math.atan2(Math.sinh(eta0), Math.cos(xi0));
  return { lat: Math.asin(s) / D2R, lon: p.lon0 + lam / D2R };
}

/**
 * Geographic → UTM.
 * @param {number} latDeg −80…84 @param {number} lonDeg
 * @param {number} [zone] force a zone (1–60); default is the zone containing the point
 * @returns {{ easting: number, northing: number, zone: number, south: boolean }}
 */
export function utmForward(latDeg, lonDeg, zone, ell = WGS84) {
  if (!(latDeg >= -90 && latDeg <= 90)) throw new RangeError(`utmForward: latDeg ${latDeg} is outside −90…90.`);
  const z = zone || utmZone(lonDeg, latDeg);
  const south = latDeg < 0;
  const r = tmForward(latDeg, lonDeg, { lon0: z * 6 - 183, k0: 0.9996, fe: 500000, fn: south ? 1e7 : 0 }, ell);
  return { easting: r.easting, northing: r.northing, zone: z, south };
}

/**
 * UTM → geographic.
 * @param {number} easting m @param {number} northing m @param {number} zone 1–60 @param {boolean} [south=false]
 * @returns {{ lat: number, lon: number }} degrees
 */
export function utmInverse(easting, northing, zone, south = false, ell = WGS84) {
  if (!(zone >= 1 && zone <= 60)) throw new RangeError(`utmInverse: zone ${zone} is outside 1…60.`);
  return tmInverse(easting, northing, { lon0: zone * 6 - 183, k0: 0.9996, fe: 500000, fn: south ? 1e7 : 0 }, ell);
}

/** True when an x/y extent can only be degrees of longitude/latitude (and is not a tiny local model). */
export function looksGeographic(xmin, xmax, ymin, ymax) {
  return xmin >= -180.000001 && xmax <= 360.000001 && ymin >= -90.000001 && ymax <= 90.000001 && (xmax - xmin) <= 360.000001;
}

/**
 * Projector from longitude/latitude (degrees) to metres, shared by every geospatial importer so that
 * datasets imported with the same `origin` line up.
 *
 * `opts.projection`:
 *  - 'local' (default): equirectangular tangent plane about `opts.origin` `{ lon, lat }`; when no origin is
 *    given the supplied data centre is used (two files then get two different origins — pass the project
 *    site as origin to overlay them);
 *  - 'utm': UTM easting/northing in `opts.utmZone` or the zone of the data centre;
 *  - 'none': leave degrees untouched (models are then NOT in metres; flagged in the description).
 * @param {object} opts @param {number} lonC data-centre longitude @param {number} latC data-centre latitude
 * @returns {{ forward: (lon:number, lat:number) => [number, number], crs: object }}
 */
export function geoProjector(opts = {}, lonC = 0, latC = 0) {
  const mode = opts.projection || 'local';
  if (mode === 'none') return { forward: (lon, lat) => [lon, lat], crs: { projection: 'none', units: 'deg', note: 'coordinates left in degrees; lengths are not metres' } };
  if (mode === 'utm') {
    const zone = opts.utmZone || utmZone(lonC, latC);
    const south = opts.utmSouth ?? latC < 0;
    const p = { lon0: zone * 6 - 183, k0: 0.9996, fe: 500000, fn: south ? 1e7 : 0 };
    return {
      forward: (lon, lat) => { const r = tmForward(lat, lon, p); return [r.easting, r.northing]; },
      crs: { projection: 'utm', zone, south, units: 'm', epsg: (south ? 32700 : 32600) + zone },
    };
  }
  if (mode !== 'local') throw new RangeError(`geoProjector: opts.projection "${mode}" is not one of local, utm, none.`);
  const lon0 = opts.origin && Number.isFinite(opts.origin.lon) ? opts.origin.lon : lonC;
  const lat0 = opts.origin && Number.isFinite(opts.origin.lat) ? opts.origin.lat : latC;
  const { M, N } = radii(lat0);
  const kx = N * Math.cos(lat0 * D2R) * D2R, ky = M * D2R;
  return {
    forward: (lon, lat) => { let dl = lon - lon0; if (dl > 180) dl -= 360; else if (dl < -180) dl += 360; return [kx * dl, ky * (lat - lat0)]; },
    crs: { projection: 'local-tangent-plane', method: 'equirectangular', origin: { lon: lon0, lat: lat0 }, originBasis: opts.origin ? 'option' : 'data centre', units: 'm' },
  };
}

// ---------------------------------------------------------------------------------------------------------
// CRS identification
// ---------------------------------------------------------------------------------------------------------

/**
 * What can be said about an EPSG code without a database: the common geographic CRSs and the UTM families.
 * @param {number} code
 * @returns {{ epsg: number, name: string, kind: 'geographic'|'projected', unit: 'deg'|'m', utm?: { zone: number, south: boolean } }|null}
 */
export function epsgInfo(code) {
  const c = Number(code);
  if (!Number.isInteger(c)) return null;
  const geo = { 4326: 'WGS 84', 4979: 'WGS 84 (3D)', 4258: 'ETRS89', 4269: 'NAD83', 4267: 'NAD27', 4230: 'ED50', 4283: 'GDA94', 7844: 'GDA2020', 4674: 'SIRGAS 2000', 4322: 'WGS 72', 4277: 'OSGB36', 4619: 'SWEREF99', 4617: 'NAD83(CSRS)', 6318: 'NAD83(2011)' };
  if (geo[c]) return { epsg: c, name: geo[c], kind: 'geographic', unit: 'deg' };
  const utm = (name, zone, south) => ({ epsg: c, name: `${name} / UTM zone ${zone}${south ? 'S' : 'N'}`, kind: 'projected', unit: 'm', utm: { zone, south } });
  if (c >= 32601 && c <= 32660) return utm('WGS 84', c - 32600, false);
  if (c >= 32701 && c <= 32760) return utm('WGS 84', c - 32700, true);
  if (c >= 25828 && c <= 25838) return utm('ETRS89', c - 25800, false);
  if (c >= 23028 && c <= 23038) return utm('ED50', c - 23000, false);
  if (c >= 26901 && c <= 26923) return utm('NAD83', c - 26900, false);
  if (c >= 26701 && c <= 26722) return utm('NAD27', c - 26700, false);
  if (c >= 32201 && c <= 32260) return utm('WGS 72', c - 32200, false);
  if (c >= 32301 && c <= 32360) return utm('WGS 72', c - 32300, true);
  if (c >= 31965 && c <= 31976) return utm('SIRGAS 2000', c - 31954, false);
  if (c >= 31977 && c <= 31985) return utm('SIRGAS 2000', c - 31960, true);
  if (c >= 28348 && c <= 28358) return utm('GDA94 / MGA', c - 28300, true);
  if (c >= 7846 && c <= 7859) return utm('GDA2020 / MGA', c - 7800, true);
  const proj = { 3857: 'WGS 84 / Pseudo-Mercator', 3395: 'WGS 84 / World Mercator', 27700: 'OSGB36 / British National Grid', 28992: 'Amersfoort / RD New', 2154: 'RGF93 / Lambert-93', 3035: 'ETRS89-extended / LAEA Europe', 3006: 'SWEREF99 TM', 25832: 'ETRS89 / UTM zone 32N' };
  if (proj[c]) return { epsg: c, name: proj[c], kind: 'projected', unit: 'm' };
  return null;
}

function parseWktTree(wkt) {
  let i = 0;
  const n = wkt.length;
  const ws = () => { while (i < n && /\s/.test(wkt[i])) i++; };
  const node = () => {
    ws();
    let s = i;
    while (i < n && /[A-Za-z0-9_]/.test(wkt[i])) i++;
    const kw = wkt.slice(s, i).toUpperCase();
    ws();
    const args = [];
    if (wkt[i] === '[' || wkt[i] === '(') {
      i++;
      for (;;) {
        ws();
        if (i >= n) break;
        const ch = wkt[i];
        if (ch === ']' || ch === ')') { i++; break; }
        if (ch === ',') { i++; continue; }
        if (ch === '"') {
          let e = i + 1, str = '';
          for (; e < n; e++) { if (wkt[e] === '"') { if (wkt[e + 1] === '"') { str += '"'; e++; } else break; } else str += wkt[e]; }
          args.push(str);
          i = e + 1;
        } else if (/[-+.\d]/.test(ch)) {
          s = i;
          while (i < n && /[-+.\deE]/.test(wkt[i])) i++;
          args.push(Number(wkt.slice(s, i)));
        } else if (/[A-Za-z_]/.test(ch)) args.push(node());
        else i++;
      }
    }
    return { kw, args };
  };
  return node();
}

/**
 * Read what matters from a coordinate-system WKT string (ESRI .prj, OGC WKT 1, WKT 2 keywords accepted).
 * @param {string} wkt
 * @returns {{ name: string, kind: 'geographic'|'projected'|'local'|'unknown', epsg: number|null, unit: string|null,
 *   unitScale: number|null, projection: string|null, params: Object<string, number>, datum: string|null,
 *   spheroid: { a: number, invF: number }|null, utm: { zone: number, south: boolean }|null, vertical: string|null }}
 *   `unit` is the canonical linear unit for projected systems ('m', 'ft', 'us-ft') or 'deg' for geographic.
 */
export function parsePrj(wkt) {
  const out = { name: '', kind: 'unknown', epsg: null, unit: null, unitScale: null, projection: null, params: {}, datum: null, spheroid: null, utm: null, vertical: null };
  if (!wkt || !/[\[(]/.test(wkt)) return out;
  let root;
  try { root = parseWktTree(wkt.trim()); } catch { return out; }
  const sub = (nd, ...kws) => nd.args.find((a) => a && typeof a === 'object' && kws.includes(a.kw));
  const subs = (nd, ...kws) => nd.args.filter((a) => a && typeof a === 'object' && kws.includes(a.kw));
  const authority = (nd) => {
    const a = sub(nd, 'AUTHORITY', 'ID');
    if (a && /epsg/i.test(String(a.args[0]))) { const c = Number(a.args[1]); return Number.isFinite(c) ? c : null; }
    return null;
  };
  let horiz = root;
  if (root.kw === 'COMPD_CS' || root.kw === 'COMPOUNDCRS') {
    const v = sub(root, 'VERT_CS', 'VERTCRS', 'VERTICALCRS');
    if (v) out.vertical = String(v.args[0] ?? '');
    horiz = sub(root, 'PROJCS', 'PROJCRS', 'PROJECTEDCRS', 'GEOGCS', 'GEOGCRS', 'GEODCRS', 'GEOGRAPHICCRS') || root;
  }
  out.name = typeof horiz.args[0] === 'string' ? horiz.args[0] : '';
  out.epsg = authority(horiz);
  const isProj = ['PROJCS', 'PROJCRS', 'PROJECTEDCRS'].includes(horiz.kw);
  const isGeog = ['GEOGCS', 'GEOGCRS', 'GEODCRS', 'GEOGRAPHICCRS', 'BASEGEOGCRS', 'BASEGEODCRS'].includes(horiz.kw);
  const geog = isProj ? sub(horiz, 'GEOGCS', 'GEOGCRS', 'BASEGEOGCRS', 'BASEGEODCRS', 'GEODCRS') : horiz;
  if (geog) {
    const d = sub(geog, 'DATUM', 'GEODETICDATUM', 'TRF');
    if (d) {
      out.datum = String(d.args[0] ?? '');
      const sp = sub(d, 'SPHEROID', 'ELLIPSOID');
      if (sp && Number.isFinite(sp.args[1]) && Number.isFinite(sp.args[2])) out.spheroid = { a: sp.args[1], invF: sp.args[2] };
    }
  }
  if (isProj) {
    out.kind = 'projected';
    const conv = sub(horiz, 'CONVERSION');
    const pj = sub(horiz, 'PROJECTION') || (conv && sub(conv, 'METHOD'));
    if (pj) out.projection = String(pj.args[0] ?? '');
    for (const p of [...subs(horiz, 'PARAMETER'), ...(conv ? subs(conv, 'PARAMETER') : [])]) {
      out.params[String(p.args[0]).toLowerCase().replace(/[\s]+/g, '_')] = Number(p.args[1]);
    }
    const units = subs(horiz, 'UNIT', 'LENGTHUNIT');
    const cs = sub(horiz, 'CS');
    const u = units[units.length - 1] || (cs && sub(cs, 'UNIT', 'LENGTHUNIT'));
    if (u) {
      const f = Number(u.args[1]);
      out.unit = normaliseUnit(u.args[0]);
      if (!out.unit && Number.isFinite(f)) out.unit = Math.abs(f - 1) < 1e-9 ? 'm' : Math.abs(f - 0.3048) < 1e-9 ? 'ft' : Math.abs(f - 1200 / 3937) < 1e-9 ? 'us-ft' : null;
      out.unitScale = Number.isFinite(f) ? f : out.unit ? UNIT_SCALE[out.unit] : null;
    } else { out.unit = 'm'; out.unitScale = 1; }
    const zm = /UTM[\s_]*zone[\s_]*(\d{1,2})\s*([NS])?/i.exec(out.name);
    const cm = out.params.central_meridian ?? out.params.longitude_of_natural_origin ?? out.params.longitude_of_origin;
    const k0 = out.params.scale_factor ?? out.params.scale_factor_at_natural_origin;
    const fe = out.params.false_easting, fn = out.params.false_northing ?? 0;
    const tm = out.projection && /transverse[_ ]mercator/i.test(out.projection);
    if (tm && Number.isFinite(cm) && Math.abs((k0 ?? 0) - 0.9996) < 1e-9 && Math.abs((fe ?? 0) * (out.unitScale || 1) - 500000) < 1 && Math.abs((((cm + 183) % 6) + 6) % 6) < 1e-9) {
      out.utm = { zone: Math.round((cm + 183) / 6), south: Math.abs(fn * (out.unitScale || 1) - 1e7) < 1 };
    } else if (zm && !out.projection) out.utm = { zone: Number(zm[1]), south: (zm[2] || 'N').toUpperCase() === 'S' };
    if (!out.epsg && out.utm && /WGS[\s_]*(19)?84/i.test(`${out.name} ${out.datum}`)) out.epsg = (out.utm.south ? 32700 : 32600) + out.utm.zone;
    if (!out.epsg && out.utm && /ETRS[\s_]*(19)?89/i.test(`${out.name} ${out.datum}`) && out.utm.zone >= 28 && out.utm.zone <= 38 && !out.utm.south) out.epsg = 25800 + out.utm.zone;
  } else if (isGeog) {
    out.kind = 'geographic';
    out.unit = 'deg';
    if (!out.epsg) {
      const d = `${out.name} ${out.datum}`;
      if (/WGS[\s_]*(19)?84/i.test(d)) out.epsg = 4326; else if (/ETRS[\s_]*(19)?89/i.test(d)) out.epsg = 4258; else if (/North_American_1983|NAD[\s_]*83/i.test(d)) out.epsg = 4269;
    }
  } else if (['LOCAL_CS', 'ENGCRS', 'ENGINEERINGCRS'].includes(horiz.kw)) {
    out.kind = 'local';
    const u = sub(horiz, 'UNIT', 'LENGTHUNIT');
    if (u) { out.unit = normaliseUnit(u.args[0]); out.unitScale = Number(u.args[1]) || (out.unit ? UNIT_SCALE[out.unit] : null); }
  }
  return out;
}

/**
 * Interpret a GeoTIFF GeoKey directory (key id → value) as far as the importers need.
 * Keys used: 1024 GTModelType (1 projected, 2 geographic), 1025 GTRasterType (1 area, 2 point),
 * 2048 GeographicType, 3072 ProjectedCSType, 3076 ProjLinearUnits (9001 m, 9002 ft, 9003 US survey ft),
 * 4096 VerticalCSType, 4099 VerticalUnits, 1026/2049/3073 citations.
 * @param {Object<number, number|string>} keys
 * @returns {{ kind: 'geographic'|'projected'|'unknown', epsg: number|null, name: string, unit: string, unitScale: number,
 *   rasterType: 'area'|'point', utm: object|null, verticalEpsg: number|null, verticalUnit: string }}
 */
export function geoKeysToCrs(keys) {
  const model = keys[1024];
  const user = 32767;
  const linear = { 9001: 'm', 9002: 'ft', 9003: 'us-ft', 9036: 'km' };
  const out = { kind: 'unknown', epsg: null, name: '', unit: 'm', unitScale: 1, rasterType: keys[1025] === 2 ? 'point' : 'area', utm: null, verticalEpsg: null, verticalUnit: 'm' };
  const pcs = keys[3072], gcs = keys[2048];
  if (model === 1 || (pcs && pcs !== user && model !== 2)) {
    out.kind = 'projected';
    if (pcs && pcs !== user) out.epsg = pcs;
    out.unit = linear[keys[3076]] || 'm';
    out.unitScale = UNIT_SCALE[out.unit];
    out.name = String(keys[3073] || keys[1026] || '');
  } else if (model === 2 || gcs) {
    out.kind = 'geographic';
    if (gcs && gcs !== user) out.epsg = gcs;
    out.unit = 'deg';
    out.unitScale = 1;
    out.name = String(keys[2049] || keys[1026] || '');
  } else out.name = String(keys[1026] || '');
  const info = out.epsg ? epsgInfo(out.epsg) : null;
  if (info) { if (!out.name) out.name = info.name; out.utm = info.utm || null; }
  if (keys[4096] && keys[4096] !== user) out.verticalEpsg = keys[4096];
  if (keys[4099]) out.verticalUnit = linear[keys[4099]] || 'm';
  return out;
}

/** "EPSG:32631"-style label for a parsePrj()/geoKeysToCrs() result, falling back to the CRS name, else null. */
export function crsLabel(info) {
  if (!info) return null;
  if (info.epsg) return `EPSG:${info.epsg}`;
  return info.name ? info.name : null;
}

// ---------------------------------------------------------------------------------------------------------
// Vertical convention
// ---------------------------------------------------------------------------------------------------------

/**
 * Decide whether a vertical column is an elevation (positive up) or a depth (positive down).
 * `opts.vertical` ('elevation' | 'depth') wins; otherwise a column/variable name is used when it is
 * unambiguous ('depth', 'water_depth', 'wd', 'sounding' → depth; 'z', 'elevation', 'height', 'altitude' →
 * elevation); otherwise elevation is assumed and the basis says so.
 * @returns {{ source: 'elevation-up'|'depth-positive-down', sign: 1|-1, basis: 'option'|'name'|'assumed' }}
 */
export function verticalConvention(opts, columnName) {
  const o = opts && opts.vertical;
  if (o === 'depth') return { source: 'depth-positive-down', sign: -1, basis: 'option' };
  if (o === 'elevation') return { source: 'elevation-up', sign: 1, basis: 'option' };
  if (o) throw new RangeError(`vertical: opts.vertical "${o}" is not one of elevation, depth.`);
  const n = String(columnName || '').toLowerCase().replace(/[^a-z]/g, '');
  if (/^(depth|waterdepth|wd|sounding|soundings|dep|depthm|tvd|bsl|depthbelowsealevel)$/.test(n)) return { source: 'depth-positive-down', sign: -1, basis: 'name' };
  if (/^(z|elev|elevation|height|altitude|alt|level|zm|elevationm|rl|ele)$/.test(n)) return { source: 'elevation-up', sign: 1, basis: 'name' };
  return { source: 'elevation-up', sign: 1, basis: 'assumed' };
}
