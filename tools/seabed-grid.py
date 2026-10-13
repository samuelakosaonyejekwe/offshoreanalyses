"""Pack the two seabed grids the app reads offline.

1. Worldwide seabed sediment map (grain-size classes, shelves and deep sea)
   Source: Shom, 2026. Carte sédimentaire mondiale / Worldwide seabed sediment map, 1:500 000,
   https://dx.doi.org/10.17183/HOM_GEOL_SEDIM_MONDIALE — "Licence Ouverte version 2.0 / Open Licence version 2.0".
   Download (shapefile inside a 7z archive, link from the metadata record HOM_GEOL_SEDIM_MONDIALE):
   https://services.data.shom.fr/INSPIRE/telechargement/prepackageGroup/SEDIM_MONDIALE_PACK_DIFF_DL/prepackage/SEDIM_MONDIALE/file/SEDIM_MONDIALE.7z

       python3 -I tools/seabed-grid.py sediment path/to/SHAPE/CarteSedimMonde_2026.shp assets/seabed-sediment.bin

   The polygons are sampled at the nodes of a 0.1 degree grid (the polygon that contains the node decides) and
   the 30 classes of the map legend are merged into the 13 groups of GROUPS below.

2. Deep-sea sediment type (what the sediment is made of, below about 500 m of water)
   Source: Diesing, M. (2020). Deep-sea sediments of the global ocean mapped with Random Forest machine learning
   algorithm. PANGAEA, https://doi.org/10.1594/PANGAEA.911692 — "Creative Commons Attribution 4.0 International".
   File lithology_classes.tif (10 km cells, Wagner IV equal-area projection) from Deep-sea_sediments_5_classes.zip.

       python3 -I tools/seabed-grid.py deepsea path/to/lithology_classes.tif assets/seabed-deepsea.bin

   Each node of a 0.1 degree grid takes the class of the source cell it falls in; nothing is interpolated.

Needs numpy and Pillow, and pyshp for the first grid; not part of the build or the tests.

Format of both outputs, decoded by src/data/seabedgrid.js:
    bytes 0-3  "SBL1"
    bytes 4-5  nx, bytes 6-7 ny (unsigned 16-bit, little-endian)
    then runs over the nodes in row order, north to south, each row west to east:
        one byte  class + 16 * (length - 1)    for runs of 1 to 15 nodes
        one byte  class + 240, then (length - 16) as a little-endian base-128 varint, for longer runs
Node (i, j) lies at longitude -180 + 360 * i / nx and latitude 90 - 180 * j / (ny - 1). Class 0 is "no data".
"""
import math
import sys

import numpy as np

NX, NY, STEP = 3600, 1801, 0.1

# legend code of the Shom map (without the "NF" prefix) -> class of the packed grid
GROUPS = {
    1: ["Roche"], 2: ["C", "CG", "CS", "CV"], 3: ["G", "GC"], 4: ["GS", "GSV", "GVS", "GV"], 5: ["SG", "SGV", "SC"],
    6: ["S"], 7: ["SV", "SVG", "SSi"], 8: ["SF", "SFG", "SFC"], 9: ["SFV", "SFSi"], 10: ["VS", "VSF", "VG", "VC", "VSG"],
    11: ["V"], 12: ["Si", "SiA"], 13: ["A", "ASi"],
}
CODE = {f"NF{c}": k for k, v in GROUPS.items() for c in v}


def pack(grid, out):
    ny, nx = grid.shape
    flat = grid.ravel()
    starts = np.concatenate(([0], np.flatnonzero(np.diff(flat)) + 1))
    lengths = np.diff(np.concatenate((starts, [flat.size])))
    buf = bytearray(b"SBL1") + nx.to_bytes(2, "little") + ny.to_bytes(2, "little")
    for cls, n in zip(flat[starts].tolist(), lengths.tolist()):
        if n <= 15:
            buf.append(cls | (n - 1) << 4)
            continue
        buf.append(cls | 0xF0)
        n -= 16
        while True:
            low, n = n & 0x7F, n >> 7
            buf.append(low | (0x80 if n else 0))
            if not n:
                break
    with open(out, "wb") as f:
        f.write(buf)
    share = np.bincount(flat, minlength=14) / flat.size
    print(f"{out}: {nx} x {ny} nodes, {len(starts)} runs, {len(buf)} bytes")
    print("nodes per class 0..13:", " ".join(str(int(c)) for c in np.bincount(flat, minlength=14)))
    print("share of nodes per class 0..13:", " ".join(f"{s:.3f}" for s in share))


def sediment(src, out):
    import shapefile

    grid = np.zeros((NY, NX + 1), dtype=np.uint8)  # one extra column for +180, folded onto -180 at the end
    skipped = {}
    with shapefile.Reader(src) as r:
        for sr in r.iterShapeRecords():
            code = CODE.get(str(sr.record[0]).strip())
            if not code:
                skipped[sr.record[0]] = skipped.get(sr.record[0], 0) + 1
                continue
            sh = sr.shape
            pts = np.asarray(sh.points, dtype=float)
            xs, ys = (pts[:, 0] + 180) / STEP, (90 - pts[:, 1]) / STEP  # nodes sit at whole numbers
            x0, y0 = max(0, int(math.floor(xs.min()))), max(0, int(math.ceil(ys.min())))
            x1, y1 = min(NX, int(math.floor(xs.max()))) + 1, min(NY - 1, int(math.floor(ys.max()))) + 1
            if x1 <= x0 or y1 <= y0:
                continue
            # exact test at the nodes: every outline edge that crosses a node row flips "inside" for the nodes to
            # its east (even-odd rule), so holes and islands inside holes come out right whatever the winding
            flips = np.zeros((y1 - y0, x1 - x0 + 1), dtype=np.int32)
            for a, b in zip(sh.parts, list(sh.parts[1:]) + [len(pts)]):
                ax, ay, bx, by = xs[a:b - 1], ys[a:b - 1], xs[a + 1:b], ys[a + 1:b]
                lo, hi = np.ceil(np.minimum(ay, by)).astype(int), np.ceil(np.maximum(ay, by)).astype(int)  # rows lo..hi-1
                n = hi - lo
                keep = n > 0
                if not keep.any():
                    continue
                idx = np.repeat(np.flatnonzero(keep), n[keep])
                rows = np.concatenate([np.arange(l, h) for l, h in zip(lo[keep], hi[keep])])
                xc = ax[idx] + (rows - ay[idx]) * (bx[idx] - ax[idx]) / (by[idx] - ay[idx])
                cols = np.clip(np.ceil(xc).astype(int), x0, x1) - x0
                ok = (rows >= y0) & (rows < y1)
                np.add.at(flips, (rows[ok] - y0, cols[ok]), 1)
            m = (np.cumsum(flips, axis=1)[:, :-1] % 2).astype(bool)
            grid[y0:y1, x0:x1][m] = code
    first = grid[:, 0]
    first[first == 0] = grid[:, NX][first == 0]
    if skipped:
        print("legend codes without a class (left as no data):", skipped)
    pack(np.ascontiguousarray(grid[:, :NX]), out)


def deepsea(src, out):
    from PIL import Image

    Image.MAX_IMAGE_PIXELS = None
    im = Image.open(src)
    scale, tie = im.tag_v2[33550], im.tag_v2[33922]  # ModelPixelScale, ModelTiepoint
    z = np.asarray(im)
    print("source", z.shape, z.dtype, "pixel", scale[:2], "origin", tie[3:5], "values", np.unique(z[np.isfinite(z)])[:12])
    lon = np.radians(-180 + STEP * np.arange(NX))
    lat = np.radians(90 - STEP * np.arange(NY))
    # Wagner IV forward projection on the WGS 84 authalic sphere: 2θ + sin 2θ = ((4π + 3√3)/6) sin φ
    R = 6371007.181
    theta = lat.copy()
    k = (4 * math.pi + 3 * math.sqrt(3)) / 6 * np.sin(lat)
    for _ in range(60):
        theta -= (2 * theta + np.sin(2 * theta) - k) / (2 + 2 * np.cos(2 * theta) + 1e-12)
    theta = np.clip(theta, -math.pi / 3, math.pi / 3)
    x = 0.86310 * R * lon[None, :] * np.cos(theta)[:, None]
    y = 1.56548 * R * np.sin(theta)[:, None] * np.ones((1, NX))
    col = np.floor((x - tie[3]) / scale[0]).astype(int)
    row = np.floor((tie[4] - y) / scale[1]).astype(int)
    ok = (col >= 0) & (col < z.shape[1]) & (row >= 0) & (row < z.shape[0])
    v = np.zeros((NY, NX))
    v[ok] = z[row[ok], col[ok]]
    v[~np.isfinite(v) | (v < 1) | (v > 13)] = 0
    pack(v.astype(np.uint8), out)


if __name__ == "__main__":
    {"sediment": sediment, "deepsea": deepsea}[sys.argv[1]](sys.argv[2], sys.argv[3])
