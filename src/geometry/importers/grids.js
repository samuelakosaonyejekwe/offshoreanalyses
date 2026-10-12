// Gridded elevation / bathymetry: ESRI ASCII grid (.asc) and Golden Software Surfer grids (.grd: DSAA text,
// DSBB Surfer 6 binary, DSRR Surfer 7 binary), plus the raster → `grid` wrap-up shared with GeoTIFF and NetCDF.
//
// Grid convention (CONTRACT §4 and src/mesh): z[iy·nx + ix] is the value at the node x = x0 + ix·dx,
// y = y0 + iy·dy with dx, dy > 0, so row 0 is the southern row. Cell-registered rasters are stored at their
// cell centres. No-data cells are NaN. z is positive up; depth rasters are negated (opts.vertical = 'depth'
// or a variable name that says depth).
//
// Neither format stores a coordinate reference system. A .prj side-car supplied with an .asc is honoured;
// otherwise coordinates are taken as projected metres unless `opts.crs` says 'EPSG:4326' (or 'geographic').

import { inputBytes, lineCursor, baseName, finish, findFile, decodeText, extOf } from './common.js';
import { resolveUnits, verticalConvention, geoProjector, looksGeographic, parsePrj, epsgInfo, crsLabel, UNIT_SCALE, normaliseUnit } from '../units.js';

/** Block-average a row-major raster by an integer factor; NaN cells are ignored within each block. */
export function blockAverage(z, nx, ny, f) {
  if (f <= 1) return { z, nx, ny };
  const mx = Math.ceil(nx / f), my = Math.ceil(ny / f);
  const sum = new Float64Array(mx * my), cnt = new Uint32Array(mx * my);
  for (let iy = 0; iy < ny; iy++) {
    const row = Math.floor(iy / f) * mx;
    for (let ix = 0; ix < nx; ix++) { const v = z[iy * nx + ix]; if (v === v) { const k = row + Math.floor(ix / f); sum[k] += v; cnt[k]++; } }
  }
  const out = new Float32Array(mx * my);
  for (let k = 0; k < out.length; k++) out[k] = cnt[k] ? sum[k] / cnt[k] : NaN;
  return { z: out, nx: mx, ny: my };
}

/** Integer decimation factor that brings nx × ny under `maxCells`. */
export function decimationFactor(nx, ny, maxCells) {
  return maxCells > 0 && nx * ny > maxCells ? Math.ceil(Math.sqrt((nx * ny) / maxCells)) : 1;
}

export const DEFAULT_MAX_CELLS = 4000000;

/**
 * Turn a decoded raster into the importer result.
 * @param {object} g `{ name, nx, ny, x0, y0, dx, dy, z }` — node coordinates in source units, row 0 south
 * @param {object} opts importer options: `maxCells`, `vertical`, `units`, `zUnits`, `crs`, `projection`, `origin`
 * @param {object} info `{ format, crs?: parsePrj()/geoKeysToCrs() result, zName?, zUnit?, warnings?, extra?, decimated?: factor already applied }`
 */
export function finishGrid(g, opts, info) {
  const warnings = info.warnings || [];
  let crs = info.crs && info.crs.kind !== 'unknown' ? info.crs : null;
  if (opts.crs) {
    const e = /(\d{4,5})/.exec(String(opts.crs));
    const known = e ? epsgInfo(Number(e[1])) : null;
    if (/^geographic$/i.test(String(opts.crs))) crs = { kind: 'geographic', epsg: 4326, name: 'WGS 84', unit: 'deg', basis: 'option' };
    else if (known) crs = { ...known, unitScale: 1, basis: 'option' };
    else crs = { kind: 'projected', epsg: e ? Number(e[1]) : null, name: String(opts.crs), unit: 'm', unitScale: 1, basis: 'option' };
  }
  const f = decimationFactor(g.nx, g.ny, opts.maxCells ?? DEFAULT_MAX_CELLS);
  let { z, nx, ny } = blockAverage(g.z, g.nx, g.ny, f);
  let x0 = g.x0 + ((f - 1) / 2) * g.dx, y0 = g.y0 + ((f - 1) / 2) * g.dy, dx = g.dx * f, dy = g.dy * f;
  if (z === g.z && !(z instanceof Float32Array)) z = Float32Array.from(z);
  const x1 = g.x0 + (g.nx - 1) * g.dx, y1 = g.y0 + (g.ny - 1) * g.dy;
  let units, crsMeta;
  if (crs && crs.kind === 'geographic') {
    if (!looksGeographic(g.x0, x1, g.y0, y1)) throw new Error(`${g.name}: the coordinate system is geographic but the grid extent (${g.x0} … ${x1}, ${g.y0} … ${y1}) is not within ±180° / ±90°.`);
    const mode = opts.projection || 'local';
    if (mode === 'none') { units = { unit: 'deg', scale: 1, basis: 'file' }; crsMeta = { projection: 'none', units: 'deg', note: 'grid left in degrees; dx and dy are not metres' }; } else {
      if (mode === 'utm') warnings.push('A regular longitude/latitude grid is not regular in UTM, so the grid is placed on the local tangent plane instead (opts.projection = "utm" is not available for grids).');
      // equirectangular about the grid centre is linear in longitude and latitude, so the lattice stays regular
      const P = geoProjector({ ...opts, projection: 'local' }, (g.x0 + x1) / 2, (g.y0 + y1) / 2);
      const o = P.forward(x0, y0), e = P.forward(x0 + dx, y0 + dy);
      dx = e[0] - o[0]; dy = e[1] - o[1]; x0 = o[0]; y0 = o[1];
      units = { unit: 'deg', scale: null, basis: 'file' };
      crsMeta = P.crs;
    }
  } else {
    units = crs && crs.unitScale && !opts.units ? { unit: crs.unit || 'm', scale: crs.unitScale, basis: 'file' } : resolveUnits(opts, info.unitHint || null, 'm');
    x0 *= units.scale; y0 *= units.scale; dx *= units.scale; dy *= units.scale;
    if (!crs && looksGeographic(g.x0, x1, g.y0, y1) && Math.max(g.dx, g.dy) < 0.5) warnings.push('The file names no coordinate system and its extent would fit longitude/latitude in degrees; it was read as projected metres. Re-import with crs = "EPSG:4326" if the grid is geographic.');
  }
  const vert = verticalConvention(opts, info.zName || '');
  let zs = vert.sign;
  const zu = opts.zUnits ? normaliseUnit(opts.zUnits) : info.zUnit ? normaliseUnit(info.zUnit) : null;
  if (opts.zUnits && !zu) throw new RangeError(`zUnits: unknown length unit "${opts.zUnits}".`);
  if (zu) zs *= UNIT_SCALE[zu]; else if (!(crs && crs.kind === 'geographic') && units.scale) zs *= units.scale;
  if (zs !== 1) for (let i = 0; i < z.length; i++) z[i] *= zs;
  if (!zu && units.scale && units.scale !== 1 && !(crs && crs.kind === 'geographic')) warnings.push(`Elevations are assumed to be in the horizontal unit (${units.unit}); set zUnits if they are not.`);
  let valid = 0;
  for (let i = 0; i < z.length; i++) if (z[i] === z[i]) valid++;
  if (!valid) warnings.push('Every cell of the grid is no-data.');
  const model = { kind: 'grid', name: g.name, nx, ny, x0, y0, dx, dy, z };
  const label = crs ? crsLabel(crs) : null;
  if (label) model.crs = label;
  const factor = f * (info.decimated || 1);
  if (factor > 1) warnings.push(`The raster (${info.sourceSize ? info.sourceSize.join(' × ') : `${g.nx} × ${g.ny}`} cells) was block-averaged by ${factor} in each direction to ${nx} × ${ny} to stay under maxCells.`);
  return finish([model], warnings, {
    format: info.format, units, vertical: vert,
    crs: crs ? { source: label, kind: crs.kind, epsg: crs.epsg ?? null, name: crs.name || null, basis: crs.basis || 'file', ...(crsMeta || {}) } : null,
    elementTypes: { gridCell: nx * ny }, skipped: { noDataCells: z.length - valid },
    decimation: { factor, source: info.sourceSize || [g.nx, g.ny] }, ...(info.extra || {}),
  });
}

/** CRS from a .prj side-car in a multi-file input, or null. */
export function sidecarCrs(input) {
  if (!input.files) return null;
  const base = baseName(input.name);
  let prj = findFile(input.files, `${base}.prj`);
  if (!prj) for (const [k, v] of input.files) if (extOf(k) === 'prj') { prj = v; break; }
  if (!prj) return null;
  const info = parsePrj(decodeText(prj));
  return info.kind === 'unknown' ? null : info;
}

/** ESRI ASCII grid (Arc/Info ASCII raster). */
export function parseEsriAscii(input, opts = {}) {
  const cur = lineCursor(input);
  const name = baseName(input.name);
  const h = {};
  let line = cur.nextNonBlank();
  while (line !== null) {
    const m = /^([A-Za-z_]+)\s+(\S+)\s*$/.exec(line);
    if (!m) break;
    h[m[1].toLowerCase()] = Number(m[2]);
    line = cur.nextNonBlank();
  }
  const nx = h.ncols, ny = h.nrows;
  const dx = h.cellsize ?? h.dx, dy = h.cellsize ?? h.dy;
  if (!(nx > 0 && ny > 0) || !(dx > 0 && dy > 0)) throw new Error(`${input.name || 'ASC'}: not an ESRI ASCII grid (expected header lines ncols, nrows, xllcorner, yllcorner, cellsize before the values).`);
  if (nx * ny > 2 ** 28) throw new Error(`${input.name || 'ASC'}: ${nx} × ${ny} cells is too large to load; resample the raster first.`);
  const centre = h.xllcenter !== undefined;
  const xll = h.xllcenter ?? h.xllcorner ?? 0, yll = h.yllcenter ?? h.yllcorner ?? 0;
  const nodata = h.nodata_value ?? -9999;
  const z = new Float32Array(nx * ny).fill(NaN);
  let k = 0;
  const total = nx * ny;
  while (line !== null && k < total) {
    const p = line.split(/\s+/);
    for (let i = 0; i < p.length && k < total; i++, k++) {
      const v = Number(p[i]);
      // file rows run north → south; the grid stores south → north
      if (v !== nodata && v === v) z[(ny - 1 - Math.floor(k / nx)) * nx + (k % nx)] = v;
    }
    line = cur.nextNonBlank();
  }
  const warnings = [];
  if (k < total) warnings.push(`The header declares ${total} cells but the file holds ${k}; the remainder is no-data.`);
  return finishGrid({ name, nx, ny, x0: centre ? xll : xll + dx / 2, y0: centre ? yll : yll + dy / 2, dx, dy, z }, opts, {
    format: 'esri-ascii', crs: sidecarCrs(input), warnings, extra: { registration: 'cell centres', nodata },
  });
}

const SURFER_BLANK = 1.70141e38;

/** Surfer grid: DSAA (ASCII), DSBB (Surfer 6 binary) or DSRR (Surfer 7 binary). */
export function parseSurferGrd(input, opts = {}) {
  const bytes = inputBytes(input);
  const name = baseName(input.name);
  const magic = decodeText(bytes.subarray(0, 4));
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  let nx, ny, x0, y0, dx, dy, z, variant;
  if (magic === 'DSAA') {
    variant = 'DSAA (ASCII)';
    const cur = lineCursor({ bytes });
    cur.next();
    const need = () => { const l = cur.nextNonBlank(); if (l === null) throw new Error(`${input.name || 'GRD'}: Surfer ASCII grid header is incomplete (expected nx ny, xlo xhi, ylo yhi, zlo zhi).`); return l.split(/\s+/).map(Number); };
    [nx, ny] = need();
    const [xlo, xhi] = need(), [ylo, yhi] = need();
    need();
    if (!(nx > 1 && ny > 1)) throw new Error(`${input.name || 'GRD'}: Surfer grid needs at least 2 × 2 nodes (found ${nx} × ${ny}).`);
    x0 = xlo; y0 = ylo; dx = (xhi - xlo) / (nx - 1); dy = (yhi - ylo) / (ny - 1);
    z = new Float32Array(nx * ny).fill(NaN);
    let k = 0;
    for (let l = cur.nextNonBlank(); l !== null && k < nx * ny; l = cur.nextNonBlank()) for (const s of l.split(/\s+/)) { if (k >= nx * ny) break; const v = Number(s); z[k++] = v >= SURFER_BLANK ? NaN : v; }
  } else if (magic === 'DSBB') {
    variant = 'DSBB (Surfer 6 binary)';
    if (bytes.length < 56) throw new Error(`${input.name || 'GRD'}: Surfer 6 binary grid is truncated (header needs 56 bytes).`);
    nx = dv.getInt16(4, true); ny = dv.getInt16(6, true);
    const xlo = dv.getFloat64(8, true), xhi = dv.getFloat64(16, true), ylo = dv.getFloat64(24, true), yhi = dv.getFloat64(32, true);
    x0 = xlo; y0 = ylo; dx = (xhi - xlo) / (nx - 1); dy = (yhi - ylo) / (ny - 1);
    if (56 + nx * ny * 4 > bytes.length) throw new Error(`${input.name || 'GRD'}: Surfer 6 binary grid declares ${nx} × ${ny} nodes but the file is too short.`);
    z = new Float32Array(nx * ny);
    for (let i = 0; i < nx * ny; i++) { const v = dv.getFloat32(56 + 4 * i, true); z[i] = v >= SURFER_BLANK ? NaN : v; }
  } else if (magic === 'DSRR') {
    variant = 'DSRR (Surfer 7 binary)';
    let p = 0, blank = SURFER_BLANK;
    while (p + 8 <= bytes.length) {
      const tag = decodeText(bytes.subarray(p, p + 4)), size = dv.getInt32(p + 4, true);
      p += 8;
      if (tag === 'GRID') {
        ny = dv.getInt32(p, true); nx = dv.getInt32(p + 4, true);
        x0 = dv.getFloat64(p + 8, true); y0 = dv.getFloat64(p + 16, true); dx = dv.getFloat64(p + 24, true); dy = dv.getFloat64(p + 32, true);
        blank = dv.getFloat64(p + 64, true);
      } else if (tag === 'DATA') {
        if (!nx) throw new Error(`${input.name || 'GRD'}: Surfer 7 grid has a DATA section before its GRID section.`);
        if (p + nx * ny * 8 > bytes.length) throw new Error(`${input.name || 'GRD'}: Surfer 7 grid declares ${nx} × ${ny} nodes but the DATA section is too short.`);
        z = new Float32Array(nx * ny);
        for (let i = 0; i < nx * ny; i++) { const v = dv.getFloat64(p + 8 * i, true); z[i] = v >= blank || v >= SURFER_BLANK ? NaN : v; }
      }
      p += size;
    }
    if (!z) throw new Error(`${input.name || 'GRD'}: Surfer 7 grid has no DATA section.`);
  } else throw new Error(`${input.name || 'GRD'}: not a Surfer grid (expected the file to start with DSAA, DSBB or DSRR). Other programs also use .grd — export the grid as ESRI ASCII (.asc) or GeoTIFF.`);
  return finishGrid({ name, nx, ny, x0, y0, dx, dy, z }, opts, { format: 'surfer-grd', crs: sidecarCrs(input), extra: { variant, registration: 'grid nodes' } });
}
