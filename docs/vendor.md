# Vendored third-party runtimes

Everything under `vendor/` is third-party code shipped unmodified so that the application works offline. Nothing
is fetched from a CDN at run time. Each runtime is loaded lazily, relative to the module that uses it, the first
time a file of the corresponding format is opened.

| Directory | Library | Version | Licence | Used by | Size |
| --- | --- | --- | --- | --- | --- |
| `vendor/occt/` | occt-import-js (OpenCASCADE Technology compiled to WebAssembly) | 0.0.23 | LGPL-2.1 (library and OCCT, with the Open CASCADE exception) | `src/geometry/adapters/occt.js` — STEP, IGES, BREP | 7 700 902 bytes |

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
| `package.json` | 23 | (local, one line: marks the directory as CommonJS so that Node can `require` the UMD glue) |

Licence notes. Both the wrapper and OpenCASCADE are under the GNU Lesser General Public License 2.1; the two
licence texts are shipped beside the binaries. The library is used as a separate, dynamically loaded file and is
not modified or linked into application code, so a user can replace it with another build of the same interface
(drop-in files of the same names). The corresponding source is available from the upstream repository at the
tag of the version above.

`tests/mesh/adapters.test.js` checks the two runtime files against the sizes and digests in this table.

### Updating

Work outside the repository so that no `node_modules` or tarballs are left behind:

```
mkdir -p /tmp/vendor-work && cd /tmp/vendor-work
npm pack occt-import-js@<version>
tar -xzf occt-import-js-<version>.tgz
cp package/dist/occt-import-js.js package/dist/occt-import-js.wasm \
   package/dist/license.occt-import-js.txt package/dist/license.occt.txt  <repo>/vendor/occt/
sha256sum <repo>/vendor/occt/*
```

Then update the version, sizes and digests above, `library.version` and `sizeBytes` in
`src/geometry/adapters/occt.js`, run `node --test tests/mesh/` and regenerate `precache.json`
(`npm run precache`) — the WebAssembly file is part of the offline cache.

## Runtimes that are not vendored

These adapters exist as stubs (`available()` returns `false`, `parse` rejects with a message naming a format to
convert to). To enable one, vendor its runtime the same way under a new directory, record it in the table above
and replace the stub with a full adapter exporting the same names (`id`, `title`, `formats`, `extensions`,
`library`, `sizeBytes`, `available`, `parse`).

| Adapter | Runtime it would need | Formats |
| --- | --- | --- |
| `src/geometry/adapters/hdf5.js` | h5wasm | HDF5: NetCDF-4 CF grids, MED, CGNS |
| `src/geometry/adapters/web-ifc.js` | web-ifc | IFC2x3 / IFC4 |
| `src/geometry/adapters/laz.js` | laz-perf | LAZ |
| `src/geometry/adapters/sqlite.js` | sql.js | GeoPackage |
