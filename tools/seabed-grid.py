"""Pack the global seafloor-lithology grid into the compact file the app reads offline.

Source: Dutkiewicz, Müller, O'Callaghan & Jónasson (2015), "Census of seafloor sediments in the world's ocean",
Geology 43(9), 795-798, doi:10.1130/G36883.1. Grid `seabed_lithology_v1.nc` (0.1 degree, 13 classes) from
https://www.earthbyte.org/webdav/ftp/papers/Dutkiewicz_etal_seafloor_lithology/ ; the same files are archived at
PANGAEA, doi:10.1594/PANGAEA.932412, under CC BY-SA 4.0. The packed file is an adaptation and carries that licence.

Usage (needs numpy and scipy; not part of the build or the tests):
    python3 -I tools/seabed-grid.py path/to/seabed_lithology_v1.nc assets/seabed-lithology.bin

Format of the output, decoded by src/data/seabedgrid.js:
    bytes 0-3  "SBL1"
    bytes 4-5  nx, bytes 6-7 ny (unsigned 16-bit, little-endian)
    then runs over the nodes in row order, north to south, each row west to east:
        one byte  class + 16 * (length - 1)    for runs of 1 to 15 nodes
        one byte  class + 240, then (length - 16) as a little-endian base-128 varint, for longer runs
Node (i, j) lies at longitude -180 + 360 * i / nx and latitude 90 - 180 * j / (ny - 1). Class 0 is "no data"
(land, ice shelves, unmapped Arctic); 1 to 13 are the lithologies of the source, unchanged. Nothing is resampled:
the duplicate column at +180 degrees is dropped, the rest is stored node for node.
"""
import sys

import numpy as np
from scipy.io import netcdf_file

src, out = sys.argv[1], sys.argv[2]

with netcdf_file(src, "r", mmap=False) as f:
    z = np.array(f.variables["z"].data, dtype=np.int16)
    lon = np.array(f.variables["lon"].data)
    lat = np.array(f.variables["lat"].data)

assert z.shape == (1801, 3601) and lon[0] == -180 and lon[-1] == 180 and lat[0] == -90 and lat[-1] == 90
assert (z[:, 0] == z[:, -1]).all(), "the +180 column is expected to repeat the -180 column"
z[z < 0] = 0  # fill value
assert z.max() <= 13
grid = z[::-1, :-1]  # north first, drop the repeated meridian
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
print("share of nodes per class 0..13:", " ".join(f"{s:.3f}" for s in share))
