# Vendored third-party runtimes

Everything under `vendor/` is third-party code shipped unmodified so that the application works offline. Nothing
is fetched from a CDN at run time. Each runtime is loaded lazily, relative to the module that uses it, the first
time a file of the corresponding format is opened. Only the files that are actually loaded are kept: no source
maps, type definitions, multi-threaded builds or bundles for other loaders.

| Directory | Library | Version | Licence | Used by | Browser download | Directory |
| --- | --- | --- | --- | --- | ---: | ---: |
| `vendor/occt/` | occt-import-js (OpenCASCADE Technology compiled to WebAssembly) | 0.0.23 | LGPL-2.1 (library and OCCT, with the Open CASCADE exception) | `src/geometry/adapters/occt.js` — STEP, IGES, BREP | 7 700 902 bytes | 7 754 891 bytes |
| `vendor/h5wasm/` | h5wasm (HDF5 2.0.0 compiled to WebAssembly) | 0.10.3 | NIST software licence; HDF5 library BSD-3-Clause | `src/geometry/adapters/hdf5.js` — NetCDF-4, MED, CGNS, Exodus II on NetCDF-4, generic HDF5 | 4 198 695 bytes | 4 206 272 bytes |
| `vendor/web-ifc/` | web-ifc (IFC geometry engine compiled to WebAssembly) | 0.0.78 | MPL-2.0 | `src/geometry/adapters/web-ifc.js` — IFC2x3 / IFC4 / IFC4x3 full geometry | 7 504 318 bytes | 14 973 803 bytes |
| `vendor/laz-perf/` | laz-perf (LASzip decoder compiled to WebAssembly) | 0.0.7 | Apache-2.0 | `src/geometry/adapters/laz.js` — LAZ | 303 406 bytes | 314 776 bytes |

"Browser download" is what a browser fetches the first time the format is opened (it equals `sizeBytes` of the
adapter); "Directory" is everything stored in the repository, including licence texts and, for web-ifc, the
separate build that Node needs.

`src/geometry/adapters/sqlite.js` (GeoPackage) has no vendored runtime: it reads the SQLite file format itself.

Each adapter has a test that checks its runtime files against the sizes and digests recorded here:
`tests/mesh/adapters.test.js` (occt), `adapters.hdf5.test.js`, `adapters.ifc.test.js` and `adapters.laz.test.js`.

## vendor/occt

Source: npm package `occt-import-js` 0.0.23, <https://github.com/kovacsv/occt-import-js>, which builds
OpenCASCADE Technology (<https://dev.opencascade.org>). The files are the unmodified contents of the package's
`dist/` directory.

| File | Bytes | SHA-256 |
| --- | ---: | --- |
| `occt-import-js.js` | 96871 | `3fb44ce11d00611f9b3f3c5775d520ebab48930c1f08279b7b1316f05f0d3379` |
| `occt-import-js.wasm` | 7604031 | `33391fc9d94ea5c869a6718488bf0a9a464222bac9bdc764dfe1690cef281952` |
| `license.occt-import-js.txt` | 27030 | `7ffe1954587c77dfba1cf8eb9b2ea743671fa6e63f9e7a2f258119d42e14eefe` |
| `license.occt.txt` | 26936 | `090f8eb6ba63e72887c1f89b009fb1656cfda02d83b6e9805b20468520a5d77f` |
| `package.json` | 23 | `d18f120db5a99372058690f86a1332484253e3cdbdbddbc0f5effa6a392b26e0` (local, one line: marks the directory as CommonJS so that Node can `require` the UMD glue) |

Licence notes. Both the wrapper and OpenCASCADE are under the GNU Lesser General Public License 2.1; the two
licence texts are shipped beside the binaries. The library is used as a separate, dynamically loaded file and is
not modified or linked into application code, so a user can replace it with another build of the same interface
(drop-in files of the same names). The corresponding source is available from the upstream repository at the
tag of the version above.

Loaded files: `occt-import-js.js` and `occt-import-js.wasm` in the browser and in Node (plus `package.json`
under Node).

## vendor/h5wasm

Source: npm package `h5wasm` 0.10.3, <https://github.com/usnistgov/h5wasm>, built on the HDF5 library 2.0.0
(<https://github.com/HDFGroup/hdf5>) through <https://github.com/usnistgov/libhdf5-wasm>. The two modules are the
unmodified `dist/esm/hdf5_hl.js` and `dist/esm/hdf5_util.js` of the package; `LICENSE.txt` is the package's.

| File | Bytes | SHA-256 |
| --- | ---: | --- |
| `hdf5_hl.js` | 48561 | `fe46fdd39690a348a61561d627f7d9feb81b92d805cde992cfa23866f1c1a8e1` |
| `hdf5_util.js` | 4150134 | `c41874d94e9523f4c7d498c951d7707f278f62040155ee8525aafce46ce1e2c8` |
| `LICENSE.txt` | 7577 | `f3ba6b8afe2a0d6f482f29a46672f88668ae02b16dfc4a2e878fd50a4f34fa6a` |

`hdf5_hl.js` is the JavaScript API (File, Group, Dataset, Attribute); it imports `hdf5_util.js`, the compiled
library with the WebAssembly binary embedded in the module, so there is no separate `.wasm` file. The same two
ES modules are loaded in the browser, in a module worker and in Node. Not vendored because nothing loads them:
the package's `dist/node` build, the `dist/iife` bundle, `dist/esm/file_handlers.js` and the TypeScript sources.

Licence notes. h5wasm is NIST-developed software provided as a public service: it may be used, copied and
distributed provided the notice in `LICENSE.txt` is kept intact, and it is not subject to copyright protection
within the United States. The HDF5 library inside it is under the HDF Group's three-clause BSD-style licence,
reproduced in the same file. Binary redistribution requires both notices to accompany the files, which
`LICENSE.txt` does.

Compression: the deflate (gzip) and shuffle filters are part of the build and are exercised by the tests.
Filters that HDF5 loads as plug-ins (Blosc, Zstandard, LZF, bzip2 …) are not shipped, so datasets compressed
with them cannot be read.

## vendor/web-ifc

Source: npm package `web-ifc` 0.0.78, <https://github.com/ThatOpen/engine_web-ifc>. Files are unmodified from
the package root.

| File | Bytes | SHA-256 |
| --- | ---: | --- |
| `web-ifc-api.js` | 5909050 | `1edb1dd8e8dba63757932f62001f1b2ee812d86594b6ffc5e1e8beebe1f16690` |
| `web-ifc.wasm` | 1595268 | `1fbd30bd5515ff6ad15268e87aa26e41d57b29411c8461618b63bee92735d349` |
| `node/web-ifc-api-node.js` | 5929555 | `592491e48f22ba5dd2f1d33a6485469bed63ac1b5f48547a80b146536556ff6f` |
| `node/web-ifc-node.wasm` | 1523182 | `36b10f4f9fada12bb89ed664f9e61eae1c82b2ce5f04ee4175205149516498f5` |
| `node/package.json` | 23 | `d18f120db5a99372058690f86a1332484253e3cdbdbddbc0f5effa6a392b26e0` (local, one line: marks the directory as CommonJS) |
| `LICENSE.md` | 16725 | `1f256ecad192880510e84ad60474eab7589218784b9a50bc7ceee34c2b91f1d5` |

Why there are two builds. The package ships one build per environment and neither runs in the other: the
ES-module build (`web-ifc-api.js` + `web-ifc.wasm`) refuses to start outside a window or worker, and the
CommonJS build (`web-ifc-api-node.js` + `web-ifc-node.wasm`) needs Node's file system. The browser and module
workers load the first pair (7 504 318 bytes); Node — the test suite and any command-line use — loads the pair
in `node/` (7 452 737 bytes) and never the first. The two JavaScript files are large because each carries the
IFC schema tables for IFC2x3, IFC4 and IFC4x3.

The files under `node/` are never requested by a browser. `tools/precache.mjs` lists everything under
`vendor/`, so unless `vendor/web-ifc/node/` is excluded there, the offline cache also stores those 7.5 MB.

Not vendored because nothing loads them: `web-ifc-mt.wasm` (multi-threaded build; it needs cross-origin
isolation, and the adapter always starts the single-threaded engine), `web-ifc-api-iife.js`, the `.d.ts` type
definitions (`ifc-schema.d.ts` alone is 2 MB) and `helpers/`.

Licence notes. Mozilla Public License 2.0, text in `LICENSE.md`. The files are distributed unmodified in
executable form; as section 3.2 of the licence asks, the source is available from the upstream repository at
the tag of the version above.

## vendor/laz-perf

Source: npm package `laz-perf` 0.0.7, <https://github.com/hobuinc/laz-perf>. `laz-perf.js` and `laz-perf.wasm`
are the unmodified `lib/laz-perf.js` and `lib/laz-perf.wasm` of the package (the build that runs in browsers,
workers and Node alike). The npm package contains no licence text; `LICENSE.txt` is the `COPYING` file of the
upstream repository.

| File | Bytes | SHA-256 |
| --- | ---: | --- |
| `laz-perf.js` | 89055 | `76712726788a65f0da0fb4476d9c46409a3b63af505e7d606c6eafdbb89a8372` |
| `laz-perf.wasm` | 214351 | `9c1802bc31b567dd4aa1ce9ab010e7e51e095dc42af8379b474ae6f221a9327f` |
| `LICENSE.txt` | 11347 | `959f77033ba56a3b146faf5c02f9162071f2d0bff4b8b6f1c2193a4b41127d39` |
| `package.json` | 23 | `d18f120db5a99372058690f86a1332484253e3cdbdbddbc0f5effa6a392b26e0` (local, one line: marks the directory as CommonJS so that Node can `require` the UMD glue) |

Not vendored because nothing loads them: the `lib/node`, `lib/web` and `lib/worker` variants (the same decoder
bound to one environment each), the `index.js` wrappers and the type definitions.

Licence notes. Apache License 2.0; the licence text accompanies the files as section 4 requires. The files are
unmodified and the npm package carries no NOTICE file.

The decoder is compiled without C++ exception support: when it meets a damaged stream it aborts, and the
adapter then discards the instance and creates a new one for the next file.

## Updating a runtime

Work outside the repository so that no `node_modules` or tarballs are left behind:

```
mkdir -p /tmp/vendor-work && cd /tmp/vendor-work
npm pack <package>@<version>
tar -xzf <package>-<version>.tgz
```

then copy exactly the files listed in the tables above:

| Runtime | From `package/` | To |
| --- | --- | --- |
| occt | `dist/occt-import-js.js`, `dist/occt-import-js.wasm`, `dist/license.occt-import-js.txt`, `dist/license.occt.txt` | `vendor/occt/` |
| h5wasm | `dist/esm/hdf5_hl.js`, `dist/esm/hdf5_util.js`, `LICENSE.txt` | `vendor/h5wasm/` |
| web-ifc | `web-ifc-api.js`, `web-ifc.wasm`, `LICENSE.md` | `vendor/web-ifc/` |
| web-ifc | `web-ifc-api-node.js`, `web-ifc-node.wasm` | `vendor/web-ifc/node/` |
| laz-perf | `lib/laz-perf.js`, `lib/laz-perf.wasm` (and `COPYING` from the repository as `LICENSE.txt`) | `vendor/laz-perf/` |

Afterwards run `sha256sum` on the copied files and update the version, sizes and digests above together with
`library.version` and `sizeBytes` in the adapter, run `node --test tests/mesh/` and regenerate `precache.json`
(`npm run precache`) — the runtimes are part of the offline cache.

## Adding another runtime

Vendor it the same way under a new directory, record it in the table at the top with its file table, and write
an adapter in `src/geometry/adapters/` exporting `id`, `title`, `formats`, `extensions`, `library`, `sizeBytes`,
`available()` and `parse()`. The runtime must be loaded lazily and relative to `import.meta.url`, and has to
work in the browser, in a module worker and in Node. An adapter whose runtime cannot be shipped returns `false`
from `available()` and rejects in `parse()` with a message that names a format to convert to.
