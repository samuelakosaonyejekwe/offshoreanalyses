# Geometry import and export

Reference for `src/geometry/registry.js`, `detect.js`, `units.js`, `zip.js`, `xml.js`, `exporters.js` and
`src/geometry/importers/*`. The WebAssembly readers in `src/geometry/adapters/*`, the mesh tools and the viewer are
described in `docs/mesh.md`.

Everything here is plain ES modules with no DOM and no `fetch`: the same code runs in Node, on the page and in a
Web Worker. Tests: `node --test tests/geometry/`.

## 1. Entry points (`registry.js`)

```js
import { FORMATS, importFile, formatFor, acceptString, adapterStatus, GeometryFormatError } from './geometry/registry.js';

const { format, models, warnings, meta } = await importFile({ name: file.name, bytes }, { units: 'mm' });
```

| Export | Purpose |
| --- | --- |
| `FORMATS` | The compatibility matrix (§5). Each entry: `{ id, name, extensions, kinds, category, status, via, suites, notes, load }`, plus `nativeSubset` / `missing` (adapter formats) or `convert` (formats that are not read). |
| `importFile(input, opts) -> Promise<{ format, models, warnings, meta }>` | Detects the format, opens `.gz` and `.zip` containers, assembles file sets, runs the importer or adapter. |
| `formatFor(name, bytes?) -> entry \| null` | Registry entry for a file from its name and first bytes (synchronous; ZIPs are identified from their listing). |
| `acceptString()` | Value for `<input type="file" accept>`: every known extension plus `.zip`, `.gz` and shapefile side-cars. |
| `adapterStatus() -> Promise<{ [name]: { available, title, sizeBytes, library } }>` | Which WebAssembly readers this build really contains. Use it to grey out adapter formats in the matrix. |
| `GeometryFormatError` | Thrown with `reason` = `'convert'` (format not read: the message says what to export instead), `'adapter-missing'` or `'unrecognised'`. |

`input` is `{ name, text?, bytes?: Uint8Array, files?: Map<string, Uint8Array> }`. Give `files` for formats made of
several files (shapefile `.shp/.dbf/.prj`, glTF with `.bin`, OpenFOAM `polyMesh`, ESRI grid with `.prj`, node and
member CSVs); a ZIP of the same files works too.

`status` values:

* `native` — read here in JavaScript; `via` is the importer module and `load()` resolves to `{ parse }`.
* `adapter` — needs the reader named in `via` (`adapter:occt`, `adapter:web-ifc`, `adapter:hdf5`, `adapter:laz`,
  `adapter:sqlite`). If the adapter module is absent, or exports `available()` returning `false`, `importFile`
  throws `GeometryFormatError('adapter-missing')` whose message ends with the conversion to make (`missing`).
  STEP and IGES carry a `nativeSubset`: that part is read without the adapter, the adapter is called only when
  the file contains exact B-rep geometry, and if the adapter is missing or fails on the file, whatever the
  built-in reader found is returned with a warning.
* `convert` — never read (Parasolid, ACIS, DWG, binary DXF, E57, Fluent, HydroStar, native CAD files).
  `importFile` throws `GeometryFormatError('convert')` with the text of `convert`.

Options understood by `importFile` itself: `format` (registry id, skips detection — archives are still opened),
`fullGeometry` (IFC → `ifc-full`), `adapters` (`{ occt: module }` to supply an adapter, `{ occt: null }` to
declare it unavailable).

## 2. Importer contract

```js
parse(input, opts) -> { models, warnings, meta }     // may return a Promise
```

Importers never throw on records they do not understand — those are counted in `meta.skipped` and summarised in
`warnings`. They throw only when the file is structurally unreadable, with a message that names the file and says
what was expected. Truncated files and header counts larger than the file are detected rather than trusted.
Text formats are read line by line over the bytes (`LineCursor`), without building one string for the file.

### 2.1 `meta`

Always present:

| Field | Content |
| --- | --- |
| `format` | registry id |
| `counts` | `{ models, byKind, vertices, triangles, cells, nodes, members, points, loops }` (those that apply) |
| `bbox` | `{ min, max }` in metres, world coordinates |
| `units` | `{ unit, scale, basis }` — `basis` is `'file'`, `'option'` or `'assumed'`; `scale` is metres per file unit |
| `crs` | `null`, or `{ source, kind, epsg, name, basis, projection?, origin?, zone? }` |
| `elementTypes` | histogram of what was read (canonical element names, or the format's own entity names) |
| `skipped` | what was not read, by record type |

Added by `importFile`: `detected` (`{ id, reason, containers }`), `fileName`, and `via` when an adapter produced the
result. Format-specific fields (schema, mapping, materials, layers, symmetry, decimation, needsAdapter …) are listed
with each format in §5.

### 2.2 Common options

| Option | Meaning |
| --- | --- |
| `units` | Length unit of the file (`'m'`, `'mm'`, `'in'`, `'ft'`, `'us-ft'` …). Overrides what the file declares. |
| `scale` | Metres per file unit, when no named unit fits. |
| `upAxis` | `'y'` rotates y-up meshes to z up (STL, OBJ, PLY, OFF; glTF is y-up by definition and always rotated). |
| `vertical` | `'elevation'` or `'depth'` — sign convention of the vertical values. Depths are negated: z is always up. |
| `zUnits` | Unit of raster elevations when it differs from the horizontal unit. |
| `crs` | `'EPSG:xxxx'` or `'geographic'` when the file does not say (or says it wrongly). |
| `projection`, `origin`, `utmZone` | How degrees become metres: `'local'` (default; tangent plane about `origin: { lon, lat }` or the data centre), `'utm'`, `'none'`. Give every import of a project the same `origin` so that the models line up. |
| `maxCells` | Rasters larger than this are block-averaged (default 4 000 000). |
| `maxPoints` | LAS clouds larger than this are thinned by stride (default 5 000 000). |
| `mapping`, `schema`, `sectionUnits`, `kinds` | Tables, JSON and XML: see §4. |

Format-specific: `variable` (NetCDF), `expandSymmetry` (GDF, NEMOH), `dxfLines` (`'auto' | 'frame' | 'section'`),
`polygons: 'polyline'` (vector formats), `axisOrder: 'xy' | 'yx'` (GML), `tessellation` (IGES surface segments per
knot span, 2–16), `chord` (airfoil chord in metres), `grid: false` (keep lattice soundings as points).

## 3. Conventions

**Units.** Models leave the importers in metres. Formats without a unit (STL, OBJ, PLY, OFF, Nastran, Abaqus,
most tables) are taken as metres, `meta.units.basis` is `'assumed'`, and a warning says so — with a hint when
the size is implausible for a structure (an extent above 5000 suggests millimetres). Nothing is rescaled
automatically. Material constants in solver decks (`meta.materials`) are reported as written, not converted.

**Axes.** z up, right-handed. Surface triangles keep the winding of the source.

**Coordinate reference systems.** No datum transformation is performed. Projected coordinates are kept as they
are (scaled to metres when the CRS unit is feet). Geographic coordinates are projected with
`units.geoProjector` — by default an equirectangular tangent plane about the origin (accurate to about 0.01 % at
10 km, 1 % at 100 km from the origin), optionally UTM. A regular longitude/latitude raster stays a regular grid on
the tangent plane; it cannot in UTM, so grids ignore `projection: 'utm'`. Files with no CRS are read as projected
metres; when the extent would also fit degrees, a warning asks for `crs`. Vertical datums are reported
(`meta.crs.verticalDatum`, `meta.verticalCrs`) and never converted.

**Large coordinates.** `surface` and `points` positions are `Float32Array`. When any coordinate exceeds 10 km the
positions are stored relative to an integer-metre `origin` on the model (world = positions + origin), which keeps
sub-millimetre resolution at UTM northings. `frame`, `polyline`, `section` and `volume` coordinates are doubles.

**Grids.** `z[iy·nx + ix]` at `x = x0 + ix·dx`, `y = y0 + iy·dy`, `dx, dy > 0`: row 0 is the southern row.
Cell-registered rasters are stored at cell centres. No-data is `NaN`. `crs` is the label of the source CRS.

**Element node order** (`ELEM` in `importers/common.js`) — the Gmsh / Abaqus corner order used by `src/mesh`:
`tet10` mid-sides on edges 01, 12, 20, 03, 13, 23; `hex20` bottom edges, top edges, vertical edges;
`wedge6` base triangle with its normal pointing to the top triangle; `line3` end, end, middle. Every reader
permutes its native order into this one (tests check mid-side nodes against corner mid-points for each format).
`hex27`, `wedge15/18` and `pyr13/14` are reduced to `hex20`, `wedge6` and `pyr5` with a warning.

**Sections.** `{ type: 'tube', D, t }` for circular hollow sections, otherwise `{ A, Iy, Iz, J, shape, … }` with
y the strong axis (Eurocode 3 naming). Nastran I1, Abaqus I11 and the IFC / catalogue strong-axis inertia map to
`Iy`. Designations are resolved by `sections.resolveSection`: pipes / CHS / tubulars, RHS / SHS / box, plates and
bars by dimensions; IPE 80–600, HE A / HE B 100–1000; a table of common AISC W shapes. Anything else keeps its
name in `member.sectionName` and is listed in a warning.

**Extra fields on models** (beyond CONTRACT §4):

| Kind | Fields |
| --- | --- |
| all FE models | `nodeIds` (source node numbers), `sets` (`{ nodes, members \| elements \| cells }` → index arrays) |
| `surface` | `origin`, `elements` (untriangulated shells `{ type, conn, ids, tags }`), `quads` (panel meshes), `stations` (hull offsets), `attrs` (PLY vertex properties), `nodeData` / `nodeDataComponents` (VTK point data); groups may carry `thickness`, `material`, `type` |
| `volume` | `cells[i].ids`, `cells[i].tags`, `tagNames`, `nodeData` |
| `frame` | `nodeNames`; members: `id`, `sectionName`, `orient`, `orientNode`, `offsets`, `releases`, `mark`, `cardinalPoint`, `rotationDeg`, `globalId`, `name`, `tag`, `container` |
| `polyline` | `attrs.kp` (metres along the line), other attribute columns, `attrUnits`, `properties`, `closed`, `layer` |
| `points` | `origin`, `attrs` (`names`, `intensity`, `classification`, `rgb`, `normals`, feature properties) |
| `section` | `layers` (`[{ name, loops: [index] }]`), `properties`, `elevation`, `chord` |

## 4. Tables, JSON and XML (schema mapping)

CSV, TSV, XYZ and whitespace tables, arrays of records in JSON, and repeated record elements in XML have no
standard geometry schema. The reader maps columns to **roles** by header name (case and punctuation ignored;
`Depth (m)`, `x[mm]`, `KP_km` carry units) and decides what the table is. Both decisions are returned —
`meta.schema` (`{ table, type, basis }`) and `meta.mapping` (`{ role: { column, index, unit } }`) — and both can be
forced: `opts.schema` and `opts.mapping = { role: 'column name' | columnIndex }`.

| `schema` | Recognised when | Result |
| --- | --- | --- |
| `points` | x / y (/ z, depth) or lon / lat columns; header-less rows `[id] x y z` | `points`; a `grid` when the points form a regular lattice (missing nodes → NaN) |
| `route` | a KP / chainage / station column, a `line` / `route` / `branch` column, a file or table name such as *route, riser, cable, pipe, mooring*, or `opts.kinds = ['polyline']` | one `polyline` per line value, with `attrs.kp` (metres; computed from the 3-D length when the file has none) and every other column as an attribute. Chainage + depth without x / y gives a profile in the x–z plane. |
| `frame` | a table with two node-reference columns (`node_a / node_b`, `n1 / n2`, `start / end`, `from / to`, `i / j`) plus a node table (`id, x, y, z`); or one table with `x1 … z2` | `frame`; sections from `D` + `t` (or inner diameter), `area / Iy / Iz / J`, or a `section` name |
| `hull` | `station`, `waterline`, `half_breadth` columns, or a matrix whose first header cell is `station` (rows) or `waterline` and whose other headers are numbers | lofted `surface`, mirrored about y = 0, groups `starboard` / `port`; open at deck and ends |
| `blade` | `r` / `span` / `radius` and `chord` | `polyline` along +z with `attrs.r`, `chord` (m), `twistDeg`, `relThickness`, `airfoil` |
| `layout` | a text id / name column with coordinates | `points` with `attrs.names` and the other columns |
| `section` | an outer-diameter column without coordinates (or forced) | `section` of circular components: `layers[i] = { name, x, y, od, id, loops }` |
| `panels` | `x1 y1 z1 … x4 y4 z4` | `surface` with `quads` |

Units: a unit in a header wins for that column; otherwise `opts.units`, a `# units: mm` comment or a `units` key in
the JSON / XML root; otherwise metres. `KP` is kilometres by definition, `chainage / station / distance` follow the
file unit. `D` / `t` without a unit follow the file unit — a warning is raised above 20 m; use `sectionUnits`.
A `depth` / `wd` / `sounding` column is negated; `z` / `elevation` is not.

JSON additionally understands `{ nodes: { "N1": [x, y, z] } }`, nested `position: [x, y, z]` / `start` / `end`,
columnar `{ x: [], y: [], z: [] }`, and hull `stations: [{ x, offsets: [[y, z]] }]`. A file with
`"format": "offshore-analyses-geometry"` is the application's own lossless format (§6).

## 5. Format matrix

`status / via`: *native*, *adapter:name* or *convert*. "+ subset" means part of the format is read natively.
The Suites column lists the suites whose specification names the format ("—": an extra that no suite asks for).

| id | Format | Extensions | status / via | Model kinds | Suites |
| --- | --- | --- | --- | --- | --- |
| `stl` | STL, ASCII and binary | .stl | native | surface | all |
| `obj` | Wavefront OBJ | .obj | native | surface, polyline, points | all |
| `ply` | Stanford PLY, ASCII and binary | .ply | native | surface, points | all |
| `off` | Geomview OFF | .off | native | surface, points | — |
| `gltf`, `glb` | glTF 2.0 | .gltf .glb | native | surface, polyline, points | 1 2 4 6 7 8 10–16 |
| `3mf` | 3MF | .3mf | native | surface | 1 11 15 16 |
| `table` | Delimited table | .csv .tsv .tab .xyz .txt .dat | native | points, grid, polyline, frame, surface, section | all |
| `route-table` | Route / centreline table | (as table) | native | polyline | 8 9 10 12–16 |
| `frame-table` | Node and member tables | (as table) | native | frame, points | 1 2 3 4 6 8 15 |
| `hull-offsets` | Hull offset table | (as table) | native | surface | 11 15 |
| `blade-table` | Blade table | (as table) | native | polyline | 16 |
| `panel-csv` | Panel list | .csv .txt .dat | native | surface | 11 15 16 |
| `table-json` | JSON tables | .json | native | as table | 1–4 6 8–16 |
| `xml-model` | XML tables | .xml | native | as table | 1–4 6 8–16 |
| `oa-json` | Application geometry JSON | .json | native | all kinds | — |
| `msh` | Gmsh MSH 2.2 / 4.1 ASCII | .msh | native | frame, surface, volume, points | all |
| `nastran` | Nastran bulk data | .bdf .nas | native | frame, surface, volume, points | all |
| `abaqus` | Abaqus input deck | .inp | native | frame, surface, volume, points | all |
| `ansys-cdb` | ANSYS CDB archive | .cdb | native | frame, surface, volume, points | 1–6 8 10–16 |
| `unv` | I-DEAS universal file | .unv | native | frame, surface, volume, points | 1–6 8–16 |
| `vtk` | VTK legacy, ASCII and binary | .vtk | native | frame, surface, volume, points, grid | all |
| `vtu` | VTK XML | .vtu .vtp | native | frame, surface, volume, points | all |
| `su2` | SU2 mesh | .su2 | native | frame, surface, volume | 7 8 11 15 16 |
| `openfoam` | OpenFOAM polyMesh | directory / zip | native (boundary patches only) | surface | 7 8 11 15 16 |
| `exodus` | Exodus II on classic NetCDF | .e .exo .ex2 .g | native | frame, surface, volume, points | all |
| `exodus-hdf5` | Exodus II on NetCDF-4 | .e .exo .ex2 .g | adapter:hdf5 | frame, surface, volume | all |
| `med` | MED (SALOME) | .med .rmed | adapter:hdf5 | frame, surface, volume | all |
| `cgns` | CGNS | .cgns | adapter:hdf5 | surface, volume | 3 7 8 11 13 15 16 |
| `fluent` | Fluent mesh / case | .msh .cas | convert | — | 7 |
| `esri-ascii` | ESRI ASCII grid | .asc .agr | native | grid | 3 4 5 7–10 12 13 15 |
| `surfer-grd` | Surfer grid DSAA / DSBB / DSRR | .grd | native | grid | — |
| `geotiff` | GeoTIFF | .tif .tiff | native | grid | 3–16 |
| `netcdf` | NetCDF classic: CF grid, GMT grid, UGRID mesh | .nc .cdf | native | grid, points, surface | 7 9 10 16 |
| `netcdf4` | NetCDF-4 | .nc .nc4 | adapter:hdf5 | grid, surface | 7 9 10 16 |
| `hdf5` | HDF5 with a documented layout | .h5 .hdf5 | adapter:hdf5 | grid, points, surface, volume | 3 |
| `geojson` | GeoJSON | .geojson .json | native | points, polyline, section | 7–10 12–16 |
| `wkt` | WKT / EWKT | .wkt | native | points, polyline, section, surface | — |
| `shapefile` | ESRI Shapefile | .shp + .dbf + .prj | native | points, polyline, section | 7 8 13 15 16 |
| `landxml` | LandXML | .xml | native | surface, polyline, points, frame | 3 4 5 7 8 10 13 15 16 |
| `gml` | GML / CityGML / InfraGML geometry | .gml .xml | native | points, polyline, section, surface | 3 5 7–10 12–16 |
| `kml` | KML / KMZ | .kml .kmz | native | points, polyline, section | — |
| `gpx` | GPX | .gpx | native | points, polyline | — |
| `geopackage` | GeoPackage | .gpkg | adapter:sqlite | points, polyline, section, grid | 5 7 8 13 15 16 |
| `las` | LAS 1.0–1.4, uncompressed | .las | native | points | all |
| `laz` | LAZ | .laz | adapter:laz | points | all |
| `pts`, `ptx` | Leica PTS / PTX | .pts .ptx | native | points | 1 2 4 5 7 8 11 14 15 16 |
| `xyzrgb` | XYZ + colour / intensity / normals | .xyzrgb .xyzn .txt | native | points | 1 2 4 5 7 8 11 14 15 16 |
| `e57` | ASTM E57 | .e57 | convert | — | all |
| `dxf` | AutoCAD DXF, ASCII | .dxf | native | frame, surface, polyline, section, points | all |
| `dxf-binary` | AutoCAD binary DXF | .dxf | convert | — | — |
| `dwg` | AutoCAD DWG | .dwg | convert | — | all |
| `sdnf` | SDNF 3.0 | .sdnf | native | frame, surface | 1–6 16 |
| `cis2` | CIS/2 STEP file | .stp .step | native | frame | 1 2 5 6 |
| `ifc` | IFC 2x3 / 4 / 4.3, linear elements and alignments | .ifc | native | frame, polyline | 1–7 10 15 16 |
| `ifc-full` | IFC, geometry of all products | .ifc | adapter:web-ifc | surface | 1–7 10 15 16 |
| `step` | STEP AP203 / AP214 / AP242 | .step .stp | adapter:occt + subset | surface, polyline, points | all |
| `iges` | IGES | .igs .iges | adapter:occt + subset | frame, polyline, surface, points, section | all |
| `brep` | OpenCascade BREP | .brep .brp | adapter:occt | surface | — |
| `parasolid` | Parasolid | .x_t .x_b | convert | — | all |
| `acis` | ACIS SAT / SAB | .sat .sab | convert | — | all |
| `native-cad` | SolidWorks, Inventor, CATIA, NX, Creo, Rhino, Revit, JT, FBX, SketchUp, DGN … | many | convert | — | 2 9 11 |
| `svg` | SVG vector profile | .svg | native | section, polyline | 12 13 14 |
| `airfoil` | Airfoil coordinates, Selig / Lednicer | .dat | native | section | 16 |
| `gdf` | WAMIT GDF, low-order | .gdf | native | surface | 11 16 |
| `nemoh` | NEMOH mesh | .dat .nemoh | native | surface | 11 16 |
| `hydrostar` | HydroStar mesh | .hst | convert | — | 11 |

The table is a summary; `FORMATS[i].notes` in the registry is the text shown to users and a test
(`tests/geometry/registry.test.js`) checks that every native entry loads, every adapter entry says what to do
without its adapter, and every convert entry names the export to make.

### 5.1 What each native reader keeps and leaves out

**Surfaces** (`importers/surfaces.js`, `gltf.js`)

* STL — facets welded into an indexed mesh (bit-identical corners), degenerate facets dropped, ASCII solids →
  groups. `UNITS=` and `ORIGIN=` in the binary header are honoured (the exporter writes them).
* OBJ — `v`, `f` (n-gons triangulated, negative indices), `l`, `p`; groups by `g`, else `o`, else `usemtl`.
  Not kept: normals, texture coordinates, materials, free-form curves and surfaces.
* PLY — vertex and face elements, ASCII / little / big endian; every other scalar vertex property → `attrs`
  (`rgb` assembled from red / green / blue). Other elements are read past.
* glTF / GLB — default scene flattened through node matrices or TRS, triangles / strips / fans, lines, points,
  sparse and normalised accessors, embedded base64 or side-car buffers. Not decoded: Draco, meshopt; not kept:
  materials, skins, animations.
* 3MF — objects, components, build transforms, unit; production-extension `path` references.

**Solver decks and meshes** (`gmsh.js`, `nastran.js`, `abaqus.js`, `ansys.js`, `unv.js`, `vtk.js`, `cfd.js`,
`netcdf.js`). Line elements → `frame`, shells → `surface` (with the untriangulated `elements`), solids →
`volume`. Property / physical-group / block tags become surface groups, member groups and `tagNames`; named sets
become `sets`. `meta.materials` lists elastic constants where the deck has them.

* Gmsh — MSH 2.x and 4.1 ASCII; physical names; MSH 4.1 entity → physical mapping. Binary and 4.0: rejected.
* Nastran — small / large / free field, continuations, `INCLUDE` of supplied files; GRID in CORD2R/C/S systems
  defined in the basic system; PBARL / PBEAML types ROD, TUBE, BAR, BOX and symmetric I. Skipped: loads, SPC/MPC,
  rigid elements, sets, PCOMP lay-ups (name only), CORD1x and chained systems.
* Abaqus — flat decks and `*PART` / `*INSTANCE` (translation + rotation), `GENERATE` sets, sets of sets,
  `*BEAM SECTION` PIPE / BOX / I / RECT / CIRC, `*BEAM GENERAL SECTION`, shell and solid sections, `*MATERIAL`.
  Skipped: steps, loads, contact, springs, masses, connectors, `*NODE, INPUT=`.
* ANSYS CDB — blocked `NBLOCK` / `EBLOCK,…,SOLID`, element routines mapped through `ET` (beams 188/189, pipes
  288/289, shells 181/281/63, planes 182/183, solids 185/186/187/285 and their legacy equivalents), degenerate
  solids recognised, `CMBLOCK` ranges, `SECTYPE` / `SECDATA`, `R`, `MP`.
* UNV — 2411, 2412, 2467 / 2477 / 2435 / 2452, 164.
* VTK — legacy ASCII and binary (incl. the 5.1 offsets / connectivity layout) for unstructured, poly, structured,
  rectilinear and image data; XML `.vtu` / `.vtp` with ascii, inline base64, appended base64 / raw, zlib,
  UInt32 / UInt64 headers. Point data kept; cell data counted as skipped.
* SU2 — volume or surface elements plus markers. OpenFOAM — boundary patches only (see the matrix note).
* Exodus II (classic NetCDF) — blocks, block names and ids, node sets; side sets skipped.

**Terrain and geospatial** (`grids.js`, `geotiff.js`, `netcdf.js`, `vector.js`, `shapefile.js`, `gisxml.js`,
`pointcloud.js`)

* ESRI ASCII — corner or centre registration, `cellsize` or `dx` / `dy`, NODATA, `.prj` side-car.
* GeoTIFF — classic and BigTIFF, both byte orders, strips or tiles, none / LZW / Deflate / PackBits, predictors
  2 and 3, 8/16/32-bit integers and 32/64-bit floats, first band, tie point + pixel scale or unrotated matrix,
  pixel-is-area / pixel-is-point, GeoKeys → EPSG, `GDAL_NODATA`. Decoding block-averages straight into the output
  lattice, so memory is bounded by `maxCells`. Rejected with the `gdal` command to run: JPEG / LERC / ZSTD / WebP,
  rotated rasters, palette and 1-bit images.
* NetCDF classic — CDF-1/2/5; variable choice by `opts.variable`, name, then position; packing, fill and valid
  range; descending latitude; `positive: down`; `grid_mapping` WKT or EPSG; legacy GMT grids; UGRID 2-D meshes.
* GeoJSON / WKT / Shapefile / GML / KML / GPX — all go through one feature collector: points → one `points`
  model with the properties as attributes; lines → `polyline` each; flat polygons → `section` (outer loop first,
  holes after); 3-D polygons → closed polylines; TIN / polyhedral / CityGML surfaces → `surface`.
  GML honours `srsName`, `srsDimension` and the latitude-first order of URN-style geographic CRS names.
* LandXML — TIN surfaces (faces flagged invisible dropped, normals turned upward), breaklines, CgPoints,
  alignments with stationing in `attrs.kp` and profile elevations (PVIs and parabolic curves; vertices are
  inserted where the profile bends), plan features, pipe networks → `frame` with tube sections. Spirals are
  drawn as a quadratic curve tangent at both ends, with the stated length for stationing.
* LAS — formats 0–10; `attrs.classification`, `intensity`, `rgb`; CRS from GeoKey or WKT VLRs; stride thinning.

**Drawings, structural and CAD text** (`dxf.js`, `sdnf.js`, `step.js`, `ifc.js`, `iges.js`, `profiles.js`,
`panels.js`)

* DXF — see the header of `dxf.js` for the entity → kind mapping. LINE entities become a frame unless the drawing
  is flat and *every* line / arc closes into loops (a profile drawing), in which case they are assembled into
  section loops; `dxfLines` overrides. Blocks are expanded through INSERT (scale, rotation, arrays, nesting);
  entities on layer 0 inside a block take the INSERT's layer.
* SDNF — Packet 10 members and Packet 20 plates, read as a token stream so that record wrapping does not matter.
* Part 21 (`P21` in `step.js`) — statements are indexed once and parsed on demand: strings with `''` and
  `\X2\…\X0\` escapes, typed values, nested lists, complex instances, comments.
* STEP — native subset as listed in the matrix; `meta.brep` counts the exact entities and
  `meta.needsAdapter = 'occt'` when they exist. Assembly placements are not applied by the subset reader.
* CIS/2 — elements and located parts; attributes are located by the type of entity they reference, so the reader
  does not depend on attribute positions of a particular LPM revision.
* IFC — units, placement chains, `Axis` or extruded `Body` representation (through mapped items and boolean
  clippings), profiles (circle, hollow circle, rectangle, hollow rectangle, I; others by name), material and
  material profile sets (including those attached to the type object), spatial container, alignments
  (polylines, indexed polycurves, composite curves; 4.3 curve segments on lines and circles, other transition
  curves as chords). `meta.otherProducts` counts what only the full-geometry adapter can show.
* IGES — fixed-format ASCII; Hollerith strings; transformation matrices; entities as listed in the matrix;
  `meta.kernelEntities` counts what needs the kernel.
* SVG — paths (all commands, arcs included), basic shapes, nested transforms, layers from named groups.
* Airfoils — Selig and Lednicer layouts detected from the data; `meta.airfoil` gives thickness and camber.
* GDF / NEMOH — panels kept as `quads`; symmetry expanded with reversed winding in mirror images.

## 6. Exporters (`exporters.js`)

| id | Function | Kinds | Notes |
| --- | --- | --- | --- |
| `stl` | `exportStl(model) -> Uint8Array` | surface, grid | binary; normals computed; header `UNITS=m ORIGIN=x y z` |
| `obj` | `exportObj(model \| models) -> string` | surface, grid, frame, polyline, points | groups kept; world coordinates |
| `ply` | `exportPly(model, { binary }) -> Uint8Array \| string` | surface, grid, points | double coordinates, `rgb` |
| `vtk` | `exportVtk(model, { scalars, title }) -> string` | surface, volume, frame, polyline, points, grid | legacy ASCII; nodal scalars (n values) or vectors (3n); a model's `nodeData` is written too |
| `csv` | `exportCsv(model) -> string` | points, polyline, frame, grid, section | frames as `[nodes]` + `[members]`, read back by the table importer |
| `geojson` | `exportGeoJson(models, { crs, origin }) -> string` | points, polyline, frame, section | with `origin` → WGS 84 longitude / latitude (RFC 7946); otherwise projected with a named `crs` member |
| `esri-ascii` | `exportEsriAscii(grid, { nodata }) -> string` | grid | `cellsize`, or `dx` / `dy` for non-square cells |
| `json` | `exportJson(models, { crs, source }) -> string` | all | application format: `{ format: 'offshore-analyses-geometry', version: 1, units: 'm', models }`, typed arrays as `{ $typed, base64 }` (little-endian bytes); lossless, NaN included |

`exportModels(id, models, { fileName }) -> { data, fileName, mime }` checks kinds and builds the file name;
`EXPORTERS` lists the formats for menus. Every exporter adds a model's `origin` back, except STL (see the table).

## 7. Detection (`detect.js`)

Order: container magic (gzip, ZIP) → binary signatures (GLB, HDF5, NetCDF, TIFF, LAS, E57, SQLite, Surfer,
shapefile, binary DXF, DWG, binary STL by exact length) → content sniffing of the first 64 kB (Part 21 schema →
`ifc` / `cis2` / `step`; XML root element; JSON markers; keyword lines of the solver decks; IGES column 73;
DXF group codes; airfoil, NEMOH, panel and coloured-point row patterns) → file extension. Extensions that are
ambiguous in practice (`.dat`, `.txt`, `.xml`, `.json`, `.msh`, `.grd`, `.asc`, `.stp`, `.nc`) are decided by
content. `detectFileSet(names)` recognises polyMesh directories, shapefile sets, 3MF packages and glTF with
side-cars. `meta.detected.reason` records why a format was chosen.

## 8. Not supported

* Exact B-rep without the `occt` adapter (STEP solids, IGES trimmed surfaces); Parasolid, ACIS, DWG, binary DXF
  and native CAD files with or without it.
* ACIS solids embedded in DXF (3DSOLID, REGION, SURFACE); DXF hatches, text, dimensions.
* IFC geometry other than linear-element axes and alignments without the `web-ifc` adapter; ifcXML.
* HDF5-based files without the `hdf5` adapter: NetCDF-4, MED, CGNS, Exodus on NetCDF-4; ADF-based CGNS always.
* LAZ without the `laz` adapter; E57; GeoPackage without the `sqlite` adapter.
* Binary Gmsh and MSH 4.0; Fluent meshes; OpenFOAM cells (boundary patches are read) and binary polyMesh;
  ANSYS `EBLOCK` without the SOLID key; UNV datasets 15 / 781.
* GeoTIFF with JPEG / LERC / ZSTD compression or a rotated transform; TopoJSON; VTK XML image / structured /
  parallel files; Draco / meshopt-compressed glTF.
* Datum transformations and vertical datum conversion; projections other than the local tangent plane and UTM.
* Assembly placements in STEP read by the built-in subset reader; IFC cardinal-point offsets of member axes.
* Lofting a blade surface from a blade table; closing a hull lofted from offsets.

## 9. Adding a format

1. Write `parse(input, opts)` in `importers/` using the helpers of `common.js` (`lineCursor`, `Grow`, `Welder`,
   `FeMesh`, `surfaceModel`, `pointsModel`, `finish`) and `units.js` (`resolveUnits`, `geoProjector`,
   `verticalConvention`). Return through `finish()` so that counts and bounds are filled.
2. Add the entry to `FORMATS` with honest `kinds`, `notes` and `suites`, and the id to `detect.js`
   (signature or sniffing rule, and `EXTENSION_MAP`).
3. Add tests that build a minimal file in the test, check known geometry, the report fields, a malformed-input
   message, and `assertTruncationSafe` for binary layouts.
