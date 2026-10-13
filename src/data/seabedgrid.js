// Reader for the two bundled 0.1° seabed grids written by tools/seabed-grid.py: assets/seabed-sediment.bin (the
// Shom worldwide seabed sediment map, Open Licence 2.0) and assets/seabed-deepsea.bin (Diesing 2020 deep-sea
// sediment types, CC BY 4.0). Pure functions on a byte buffer, so the same code runs in the browser and under Node.

/**
 * Classes of the sediment grid: groups of the legend of the Shom map. `code` is the legend code that stands for
 * the group (decoded by shomClass in seabed.js, which also gives the grain-size band of the legend definition).
 */
export const SEABED_GRID_CLASSES = {
  1: { code: 'NFRoche' }, 2: { code: 'NFC' }, 3: { code: 'NFG' }, 4: { code: 'NFGS' }, 5: { code: 'NFSG' }, 6: { code: 'NFS' }, 7: { code: 'NFSV' },
  8: { code: 'NFSF' }, 9: { code: 'NFSFV' }, 10: { code: 'NFVS' }, 11: { code: 'NFV' }, 12: { code: 'NFSi' }, 13: { code: 'NFA' },
};
/** Classes of the deep-sea grid (what the sediment is made of), each mapped to a key of SEABED_CLASSES. */
export const DEEPSEA_GRID_CLASSES = {
  1: { name: 'Calcareous sediment', cls: 'calcOoze' }, 2: { name: 'Clay', cls: 'clay' }, 3: { name: 'Diatom ooze', cls: 'silOoze' },
  4: { name: 'Lithogenous sediment', cls: 'lithogenous' }, 5: { name: 'Radiolarian ooze', cls: 'silOoze' },
};

const bytesOf = (buf) => (buf instanceof Uint8Array ? buf : new Uint8Array(buf));
/** Grid size from the 8-byte header ("SBL1", nx, ny). */
export function seabedGridSize(buf) {
  const b = bytesOf(buf);
  if (b.length < 8 || String.fromCharCode(b[0], b[1], b[2], b[3]) !== 'SBL1') throw new Error('Not a seabed grid file');
  return { nx: b[4] | b[5] << 8, ny: b[6] | b[7] << 8 };
}

/** Codes at the given node indices (row-major, ascending) in one pass over the runs; -1 past the end of the data. */
export function seabedGridCodes(buf, indices) {
  const b = bytesOf(buf), out = new Array(indices.length).fill(-1);
  let p = 8, end = 0, k = 0;
  while (p < b.length && k < indices.length) {
    const code = b[p] & 15;
    let n = (b[p++] >> 4) + 1;
    if (n === 16) { let shift = 0, v; do { v = b[p++]; n += (v & 127) * 2 ** shift; shift += 7; } while (v & 128); }
    end += n;
    while (k < indices.length && indices[k] < end) out[k++] = code;
  }
  return out;
}

/** Every node of the grid as one byte each (for checks and tools; the app itself only reads single points). */
export function decodeSeabedGrid(buf) {
  const b = bytesOf(buf), { nx, ny } = seabedGridSize(b), cells = new Uint8Array(nx * ny);
  let p = 8, at = 0;
  while (p < b.length) {
    const code = b[p] & 15;
    let n = (b[p++] >> 4) + 1;
    if (n === 16) { let shift = 0, v; do { v = b[p++]; n += (v & 127) * 2 ** shift; shift += 7; } while (v & 128); }
    cells.fill(code, at, at + n); at += n;
  }
  if (at !== cells.length) throw new Error('Seabed grid file is truncated');
  return { nx, ny, cells };
}

/**
 * Class at a point: the nearest grid node, or — where that node holds no data, as it does along coasts drawn
 * coarser than the real shoreline — the closest node with data within `reach` nodes. Null on land and wherever
 * the source has no value. `offsetKm` is the distance from the point to the node that answered.
 */
export function seabedGridAt(buf, lat, lon, { reach = 2, classes = SEABED_GRID_CLASSES } = {}) {
  if (!Number.isFinite(lat) || !Number.isFinite(lon) || Math.abs(lat) > 90) return null;
  const { nx, ny } = seabedGridSize(buf), step = 360 / nx, R = Math.PI / 180;
  const x = (((lon + 180) % 360) + 360) % 360 / step, y = (90 - lat) / step;
  const i0 = Math.round(x), j0 = Math.round(y), cand = [];
  for (let dj = -reach; dj <= reach; dj++) for (let di = -reach; di <= reach; di++) {
    const j = j0 + dj, i = ((i0 + di) % nx + nx) % nx;
    if (j < 0 || j >= ny) continue;
    const nLat = 90 - j * step, km = 6371 * R * Math.hypot((i0 + di - x) * step * Math.cos(lat * R), nLat - lat);
    cand.push({ index: j * nx + i, km, direct: !di && !dj, lat: nLat, lon: i * step - 180 });
  }
  cand.sort((a, b) => a.index - b.index);
  const codes = seabedGridCodes(buf, cand.map((c) => c.index));
  const hit = cand.map((c, k) => ({ ...c, code: codes[k] })).filter((c) => c.code > 0).sort((a, b) => b.direct - a.direct || a.km - b.km)[0];
  return hit ? { code: hit.code, ...classes[hit.code], lat: hit.lat, lon: hit.lon, offsetKm: hit.km, shifted: !hit.direct, stepDeg: step } : null;
}
