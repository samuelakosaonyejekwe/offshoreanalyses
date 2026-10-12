# Mesh tools, 3-D viewer and CAD adapters

Reference for `src/mesh/*`, `src/ui/viewer3d.js` and `src/geometry/adapters/*`. Models are the kinds of
`docs/CONTRACT.md` §4. Everything in `src/mesh` is pure, DOM-free and runs under Node; lengths are metres, z is up.

Conventions shared by all files

* Triangle surfaces: `positions` (x, y, z per vertex), `indices` (three ids per face). Outward orientation is
  counter-clockwise seen from outside, so the divergence-theorem volume is positive.
* Surface `groups` are `{ name, start, count }` in **index units** (as WebGL draw ranges): face range
  `start / 3 … (start + count) / 3`. Every operation that changes faces carries the groups through.
* Volume cells use the Gmsh / Abaqus / CGNS corner order (`CELL_FACES`, `CELL_EDGES` in `topology.js`). Quadratic
  cells (tet10, hex20) are treated on their corner nodes.
* Grid: `z[iy·nx + ix]` at `(x0 + ix·dx, y0 + iy·dy)`; NaN marks no data.
* Functions never modify their input; invalid arguments throw `RangeError` naming the argument.
* Tests: `node --test tests/mesh/` (benchmarks are written out in the test names).

## topology.js

| Export | Purpose |
| --- | --- |
| `bounds(positions, dim = 3) -> { min, max, diagonal }` | axis-aligned bounds |
| `weldVertices(positions, tol = 0) -> { positions, map, count, merged, rep }` | spatial-hash merge; `map` old → new |
| `weld(surface, tol?) -> { model, merged, map }` | weld and remap a surface (default tol 1e-6 × diagonal) |
| `buildEdges(indices, nVerts) -> { edgeCount, edgeVerts, faceEdges, edgeFaceStart, edgeHalfEdges, edgeFaceCount(e) }` | edge–face adjacency |
| `boundaryLoops(indices, nVerts, edges?) -> { loops, closed, edgeCount }` | chains of one-face edges |
| `nonManifoldEdges(indices, nVerts, edges?) -> [[i, j]]` | edges with more than two faces |
| `nonManifoldVertices(indices, nVerts, edges?) -> number[]` | bow-tie vertices |
| `connectedComponents(indices, nVerts, by = 'edge' \| 'vertex', edges?) -> { count, faceComponent, sizes }` | |
| `duplicateFaces(indices, nVerts) -> number[]`, `degenerateFaces(positions, indices, relTol = 1e-7) -> number[]` | face ids |
| `orientConsistently(positions, indices, nVerts?) -> { indices, flipped, orientable, components, componentVolume, componentClosed }` | consistent winding, closed parts turned outward |
| `eulerCharacteristic(indices, nVerts, edges?) -> { V, E, F, chi, components, boundaryLoops, genus }` | χ = V − E + F = 2C − 2g − b |
| `boundaryFaces(volume) -> { tris, quads, triCell, quadCell, interiorFaces, overSharedFaces }` | faces owned by one cell (typed hash table, suited to millions of cells) |
| `cellVolume(positions, conn, offset, type) -> number` | signed cell volume |
| `toSurface(model) -> surface \| null` | surface as is, volume boundary (with `faceCell`), grid height-field |
| `selfIntersections(positions, indices, { maxPairs = 1000 }) -> { pairs, truncated }` | crossing face pairs (touching and coplanar overlap are not reported) |
| `inspect(model, { tolerance?, weld?, selfIntersections? }) -> { ok, watertight, reasons?, issues, stats, kind, name }` | any kind, see below |
| `watertight(model, opts?) -> { watertight, reasons }` | verdict for a surface or a volume boundary |
| `faceGroupIds(surface) -> { ids, names }`, `groupsFromIds(indices, ids, names) -> { indices, groups, order }` | groups ↔ per-face ids |
| `CELL_INFO`, `CELL_FACES`, `CELL_EDGES` | cell tables |

`inspect` issues are `{ code, count, severity: 'error' | 'warn' | 'info', message, edges?, elements?, vertices? }`;
`ok` is true when nothing has severity `error`. `edges`, `elements` and `vertices` (first 5000) index the model as
given and can be passed straight to `viewer.highlight`.

| Kind | Codes |
| --- | --- |
| surface | `non-finite-coordinates`, `index-out-of-range`, `empty`, `unwelded-vertices`, `degenerate-faces`, `duplicate-faces`, `open-boundary`, `non-manifold-edges`, `non-manifold-vertices`, `unreferenced-vertices`, `not-orientable`, `inconsistent-orientation`, `multiple-components`, `self-intersections` |
| volume | `unsupported-cell-type`, `bad-connectivity-length`, `index-out-of-range`, `inverted-cells`, `zero-volume-cells`, `unreferenced-nodes`, `empty`, `over-shared-faces`, `non-conforming-boundary` |
| frame | `index-out-of-range`, `zero-length-members`, `duplicate-members`, `coincident-nodes`, `unreferenced-nodes`, `free-ends`, `disconnected-parts`, `members-without-section` |
| polyline / points / grid / section | `non-finite-coordinates`, `repeated-points`, `too-few-points`, `duplicate-points`, `bad-grid-size`, `zero-spacing`, `no-data-cells`, `degenerate-loops`, `repeated-closing-point`, `self-intersecting-loops`, `empty`, `unknown-kind` |

Surface `stats`: `vertices, weldedVertices, faces, edges, eulerCharacteristic, genus, components, boundaryLoops,
boundaryEdges, closed, manifold, orientable, outward, area, volume (NaN when open), bbox,
selfIntersectionsChecked`. A surface is inspected after welding within the tolerance (STL-style soups are
reported as `unwelded-vertices`, an `info`); pass `weld: false` to judge the connectivity exactly as stored.
The self-intersection search is on by default up to 60 000 faces.

## measure.js

| Export | Returns |
| --- | --- |
| `surfaceArea(positions, indices)` | m² |
| `signedVolume(positions, indices)` | m³, positive for closed outward surfaces |
| `inertiaSolid(positions, indices, { density = 1 })` | `{ volume, mass, centroid, inertia, inertiaOrigin }` — products carry the minus sign (Ixy = −∫xy dm) |
| `inertiaShell(positions, indices, { arealDensity \| thickness, density })` | `{ area, mass, centroid, inertia, inertiaOrigin }` thin-shell |
| `principalAxes(inertia, mass?)` | `{ moments (ascending), axes, radiiOfGyration? }` |
| `aabb(positions)`, `obb(positions)` | `{ min, max, size, center, diagonal }`, `{ center, axes, halfSizes, size, volume }` (PCA box) |
| `projectedArea(positions, indices, direction, { resolution = 1024 })` | `{ area, cellSize, width, height, centroid, errorBound }` silhouette (wind / current area) |
| `sectionProperties(loops)` | `{ area, centroid, Ixx, Iyy, Ixy, I1, I2, thetaDeg, perimeter, rx, ry }` — holes by nesting |
| `planeCut(positions, indices, { point, normal })` | `{ loops, loops2d, closed, basis, origin, area, centroid, centroid2d, Iuu, Ivv, Iuv, I1, I2, perimeter }` — waterplanes, frame sections |
| `volumeBelowPlane(positions, indices, { point, normal })` | `{ volume, centroid, wettedArea }` — displacement, centre of buoyancy |
| `wallThickness(positions, indices, { samples = 256, seed = 1 })` | `{ min, max, mean, median, p05, p95, samples, misses, values, faces }` |
| `rayTriangle(org, dir, P, ia, ib, ic)` | ray parameter or Infinity |
| `polylineMeasure(points)` | `{ length, horizontalLength, chainage, curvature, inclinationDeg, headingDeg, minRadius, maxInclinationDeg, zMin, zMax }` |
| `gridSample(grid, x, y)` | bilinear z or NaN |
| `gridStats(grid)` | `{ nodes, valid, noData, zMin, zMax, zMean, zStd, depthMin, depthMax, slopeMeanDeg, slopeMaxDeg, slopeP95Deg, roughness, extentX, extentY, area }` |
| `gridProfile(grid, points, { step?, stride = 3 })` | `{ chainage, x, y, z, lineZ, gap, slopeDeg, outside }` seabed under a route |
| `planeBasis(normal)` | `{ u, v, n }` |
| `measure(model)` | headline numbers for any kind |

## repair.js

`repair(model, opts) -> { model, actions: [{ code, count, message }], flags? }`

| Kind | Sequence | Action codes |
| --- | --- | --- |
| surface | weld → remove degenerate → remove duplicate → drop tiny components → split bow-tie vertices → fill small holes → consistent outward orientation → drop unused vertices | `welded-vertices`, `removed-degenerate-faces`, `removed-duplicate-faces`, `dropped-small-components`, `split-non-manifold-vertices`, `filled-holes`, `holes-left-open`, `flipped-faces`, `not-orientable`, `removed-unreferenced-vertices`, `non-manifold-edges-remain` |
| volume | drop zero-volume cells, reorder inverted cells, drop unused nodes, remap sets | `reordered-inverted-cells`, `removed-zero-volume-cells`, `inverted-cells-remain`, `removed-unreferenced-nodes` |
| frame | `healFrame` | `merged-nodes`, `split-members`, `removed-zero-length`, `removed-duplicates`, `removed-unused-nodes`; `flags: { freeEndNodes, overlappingMembers, shortMembers, isolatedParts }`; model gets `memberOrigin` |
| polyline, points | drop invalid and repeated points | `removed-points` |
| grid | fill no-data gaps up to `maxGap` (2) nodes | `filled-no-data` |
| section | clean loops, outer counter-clockwise, holes clockwise | `cleaned-loop-points`, `reversed-loops`, `removed-degenerate-loops` |

Surface options: `{ tolerance, weld, fillHoles, maxHoleEdges = 32, fixOrientation, splitNonManifoldVertices,
dropSmallComponents, minComponentFaces = 4, minComponentAreaFraction = 1e-4 }` (booleans default true).
Frame options: `{ tolerance (default max(1 mm, 1e-4 × diagonal)), split = true, removeUnusedNodes = true }`.

Single steps: `removeDegenerateFaces`, `removeDuplicateFaces`, `removeUnreferencedVertices`,
`fillHoles(model, { maxEdges })`, `fixOrientation`, `splitNonManifoldVertices`,
`dropSmallComponents(model, { minFaces, minAreaFraction })`, `healFrame(frame, opts)`.

Not attempted: resolving self-intersections, closing large gaps by surface reconstruction, and separating edges
shared by more than two sheets. These are reported and left to the engineer.

## quality.js

`assess(model, { limits?, bins = 10, panel? }) -> { perElement, histogram, histograms, summary, grade, diagnosis, panel? }`

* `perElement`: `aspectRatio, minAngle, maxAngle, skewness, scaledJacobian, size` (Float32Array, element order)
  and `rating` (Uint8Array: 0 good, 1 marginal, 2 poor) — ready for a colour plot.
* `histogram`: `{ metric: 'scaledJacobian', edges, counts }`; `histograms` also has `skewness`, `minAngle`,
  `aspectRatio`.
* `summary`: `{ elements, byType, <metric>: { min, max, mean, std }, size: { …, ratio }, marginal, poor, inverted, worst: { element, metric, value } }`.
* `grade`: A all good · B no poor and ≤ 5 % marginal · C ≤ 1 % poor, none inverted · D otherwise · `-` no elements.
* `diagnosis`: plain-language sentences (what is wrong, why it matters, what to do).
* `panel` (surfaces): result of `panelCheck`.

Metrics follow the Verdict library: triangle aspect ratio L_max·ΣL / (4√3·A); quadrilateral / hexahedron /
wedge / pyramid aspect ratio longest ÷ shortest edge; tetrahedron aspect ratio R / 3r; equiangle skewness
(tetrahedron: equivolume); scaled Jacobian normalised to 1 for the ideal element. Default limits
(`QUALITY_LIMITS`, `[warn, fail]`) are the usual solver warning levels and can be overridden per family
(`tri`, `quad`, `tet`, `hex`, `wedge`, `pyr`).

Other exports: `triQuality(P, a, b, c)`, `quadQuality(P, a, b, c, d)`, `cellQuality(P, conn, offset, type)`,
`panelFrequencyLimit(panelSize, { panelsPerWavelength = 7, depth = ∞, g }) -> { minWavelength, maxOmega, maxFrequencyHz, minPeriod }`,
`panelCheck(surface, { panelsPerWavelength, depth, g, waterline = 0, minPeriod? }) -> { panels, panelSizeMax,
panelSizeP95, panelSizeMean, minWavelength, maxOmega, maxFrequencyHz, minPeriod, normalsConsistent,
normalsOutward, reversedPanels, closed, panelsAboveWaterline, wettedPanels, requiredPanelSize?, adequate?, notes }`.
The frequency limit uses λ = n·h with the 95th-percentile wetted panel size and ω² = g·k·tanh(k·d).

## refine.js

| Export | Notes |
| --- | --- |
| `refine(model, { project? })` | one uniform level: surface 1→4, volume 1→8, frame members halved, polyline and grid spacing halved |
| `refineTriangles(positions, indices, project?)`, `refineQuads(positions, quads, project?)` | element-level 1→4 |
| `refineVolume(model, { project? })` | tet 1→8 (Bey), hex 1→8, wedge 1→8, pyramid 1→6 pyramids + 4 tetrahedra; conforming for mixed meshes; blocks carry `parent`; `sets.cells` mapped to children |
| `refineTo(model, targetSize, { project? })` | surface: longest-edge bisection until every edge ≤ target; volume: uniform levels; frame: `subdivideMembers`; polyline / grid: resampling |
| `subdivideMembers(frame, { parts \| maxLength })` | adds `memberParent` |
| `resamplePolyline(polyline, { spacing \| count })` | arc-length resampling, numeric attributes interpolated, `attrs.kp` added |
| `resampleGrid(grid, { dx, dy? } \| { factor })` | bilinear when refining, cell mean when coarsening |
| `decimate(surface, targetCount, { borderWeight = 100, maxNormalTurnDeg = 80 })` | quadric-error edge collapse; keeps manifoldness, boundaries and group borders; adds `decimation: { from, to, maxError }` |
| `coarsenTo(model, targetCount, opts?)` | surface → `decimate`; points → every k-th; polyline, grid → resampled. Volume and frame models throw |
| `smooth(surface, { iterations = 10, lambda = 0.5, taubin = true, mu = −lambda, fixBoundary, fixGroupBorders, fixed })` | Laplacian / Taubin; boundaries and group borders held |

`project(x, y, z) -> [x, y, z]` snaps new nodes to the true geometry (needed for curved bodies in a convergence
study: without it refinement only subdivides the facets).

## generate.js

All bodies are closed, outward-oriented `surface` models (Float32 positions) with named groups.

| Export | Groups |
| --- | --- |
| `box({ size, center, divisions \| elementSize })` | `x-`, `x+`, `y-`, `y+`, `z-`, `z+` |
| `cylinder({ radius, height, base, segments, elementSize, heightDivisions, capRings })` | `bottom`, `side`, `top` |
| `cone({ radiusBottom, radiusTop, height, … })` | `bottom`, `side`, `top` (a zero radius gives an apex) |
| `sphere({ radius, center, subdivisions = 3 \| elementSize })` | `surface` (icosphere, 20·4ˢ faces) |
| `tube({ a, b, radius, radiusB?, … })` | `bottom` (at a), `side`, `top` (at b) |
| `barge({ length, beam, depth, draft, rakeBow, rakeStern, rakeHeight, elementSize \| divisions })` | `stern`, `bow`, `starboard`, `port`, `bottom`, `deck` |
| `wigleyHull({ length, beam, draft, nx = 40, nz = 12, closed = true })` | `hull-starboard`, `hull-port`, `deck` |
| `semiSubmersible({ pontoonLength, pontoonWidth, pontoonHeight, pontoonSpacing, columnDiameter, columnHeight, columnsPerPontoon = 2, columnSpacing, draft, segments = 16, elementSize })` | `pontoons`, `columns`, `column-tops` — two closed hulls (χ = 4), column feet imprinted in the pontoon decks |
| `spar({ diameter, draft, freeboard = 10, upperDiameter, taperTop, taperBottom, segments, elementSize })` | `keel`, `lower`, `taper`, `upper`, `top` |
| `monopile({ diameter, zTop, zToe, topDiameter, taperTop, taperBottom, wallThickness, segments, elementSize })` | `toe`, `shaft`, `taper`, `upper`, `top`, `inner` (hollow: genus 1) |

Floating bodies have the waterline at z = 0 (keel at −draft); bow at +x, port at +y.

* `triangulateSection(section, { size?, minAngle = 25, maxVertices = 200000 })` — Ruppert Delaunay refinement of
  a section with holes (even–odd rule). Returns a planar `surface` (z = 0, normal +z) with `xy` (Float64Array),
  `boundary` (vertex pairs) and `stats: { vertices, triangles, minAngleDeg, maxEdge, area, truncated }`. The
  boundary is kept exactly (segments are subdivided); the result is Delaunay. Loops must be simple and must not
  touch. Input angles below about 30° cannot be improved and are left as they are.
* `delaunay(xy) -> Uint32Array` — unconstrained Delaunay triangulation of a point set.
* `pointsToGrid(points, { dx, dy?, power = 2, neighbours = 12, radius?, bounds?, name? }) -> grid` — inverse
  distance weighting; nodes with no point within the radius (default 5 × spacing) are NaN.
* `gridToSurface(grid, { stride = 1, name? }) -> surface` with `gridNode`.

## convergence.js

Index convention of Celik et al. (2008): 1 = finest grid; r21 = h2 / h1 > 1. Errors and GCI are fractions.

| Export | Returns |
| --- | --- |
| `classify(f1, f2, f3, { r21, r32, tol })` | `{ type: 'monotone' \| 'oscillatory' \| 'divergent' \| 'oscillatory-divergent' \| 'converged' \| 'indeterminate', R }` |
| `observedOrder({ f1, f2, f3, r21, r32 })` | `{ p, s, iterations, converged }` — three-grid formula solved iteratively |
| `richardson({ f1, f2, r21, p })` | extrapolated value |
| `gci({ f1, f2, r21, p, Fs = 1.25 })` | `{ gci, ea, eext, extrapolated }` — Fs 1.25 for three grids, 3.0 for two |
| `asymptoticCheck({ f1, f2, f3, r21, r32, p, tol = 0.05 })` | `{ ratio, inRange, gci21, gci32 }` — GCI32 / (r21^p · GCI21) ≈ 1 |
| `threeGrid({ f1, f2, f3, r21, r32 \| h }, { formalOrder?, Fs? })` | `{ type, R, p, pUsed, extrapolated, gciFine, gciMedium, ea21, eext21, asymptoticRatio, inAsymptoticRange, Fs, notes }` |
| `fitOrder(h, f, { pMin, pMax, weights })` | `{ f0, C, p, residualStd, atBound }` least squares for more than three grids |
| `selectLevel({ levels: [{ h, value, cost? }], extrapolated, tolerance, Fs, costExponent })` | `{ index, level, errors, relativeCost, saving, met }` coarsest admissible level |
| `runStudy({ levels, evaluate, tolerance = 0.01, onProgress, isSize, dimension, formalOrder, signal })` | Promise of `{ levels: [{ level, h, value, cost, error }], type, R, p, pUsed, extrapolated, gciFine, gciMedium, ea21, asymptoticRatio, inAsymptoticRange, Fs, fit, selected, converged, tolerance, notes }` |
| `errorNorms(numerical, exact, weights?)` | `{ l1, l2, linf, relL2, n }` |
| `orderFromErrors(h, errors)` | `{ p, pairwise, coefficient }` |

`runStudy` matches the `mesh` ladder of CONTRACT §2.3: a level is a count (h ∝ level^(−1/dimension)) unless
`isSize` is true. Benchmark in the tests: Celik et al. (2008) Table 1, φ = 6.063, 5.972, 5.863 with r21 = 1.5,
r32 = 1.333 gives p = 1.53, φ_ext = 6.168, GCI_fine = 2.2 %.

## viewer3d.js

```js
import { createViewer } from './ui/viewer3d.js';
const viewer = createViewer(container, { edges, view: 'iso', ortho: false, label, controls: true, renderer: 'auto' });
```

| Method | Notes |
| --- | --- |
| `setScene({ models, styles? }) -> ids` | replaces everything and fits; a scene without height is shown from above |
| `addModel(model, style?) -> id`, `removeModel(id)`, `setStyle(id, style)` | style: `{ color, opacity, edges, edgeColor, lineWidth, pointSize, zScale, maxPoints = 300000, maxGridNodes = 400, colourBy: 'z' \| 'none', visible, range, smooth }` |
| `setNodeValues(id, Float32Array \| null, { label, unit, range: [min, max] })` | one value per model node; NaN is drawn grey |
| `setDeformed(id, positions \| null, scale = 1 \| 'auto') -> scale` | displayed shape = base + scale·(deformed − base); the undeformed outline stays as a ghost |
| `highlight(id, { edges, elements, vertices })` | ids as in the model (and as returned by `inspect`); omit the second argument to clear |
| `clear()`, `fit()`, `setView('iso' \| 'top' \| 'front' \| 'side')`, `setEdges(on)`, `setProjection('perspective' \| 'orthographic')`, `setVerticalScale(k)` | |
| `onPick(cb) -> unsubscribe`, `pick(x, y)` | `cb({ id, kind, name, point, element, node, value } \| null)` |
| `snapshot() -> PNG data URL`, `resize()`, `dispose()`, `getCamera()`, `setCamera(cam)`, `summary()`, `canvas`, `mode` | `mode` is `'webgl2'`, `'webgl'` or `'2d'` |

What is drawn: surfaces shaded with optional edges; volumes by their boundary faces (face outlines as edges);
frames, polylines and section outlines as lines of constant screen width; point clouds decimated to `maxPoints`
and coloured by elevation; grids as a shaded height-field coloured by elevation with automatic vertical
exaggeration (shown beside the dimensions; `setVerticalScale` applies one factor to the whole scene, e.g. a
route over bathymetry); `model3d` plot objects with `nodes`, `lines`, `tris`, `nodeValues`, `valueLabel`,
`valueUnit`, `deformed`.

Interaction: drag orbits (arcball), right / middle drag or Shift + drag pans, wheel zooms about the cursor; one
finger orbits, two fingers pan and pinch; keys: arrows orbit (Shift + arrows pan), `+` / `−` zoom, `F` fit,
`E` edges, `1`–`4` views. The canvas is focusable, has an `aria-label` and is described by a text summary of the
scene (`summary()`).

Behaviour: renders only when something changes; device pixel ratio capped at 2; recovers from WebGL context
loss; falls back to a 2-D canvas wireframe when WebGL is unavailable; colours come from the CSS custom
properties `--surface-2`, `--ink`, `--accent`, `--grid` and are re-read on the window `themechange` event; nodal
values use a perceptually uniform, colour-blind-safe sequential map with a colour bar (label, unit, min, max).
Buffers are stored relative to each model's centre, so UTM-scale coordinates keep their precision.

Pure functions exported for reuse and tested under Node: matrices (`mat4Identity`, `mat4Multiply`,
`mat4Perspective`, `mat4Ortho`, `mat4FromQuat`, `projectPoint`), quaternions (`quatNormalize`, `quatMultiply`,
`quatFromAxisAngle`, `quatRotate`, `quatFromBasis`), arcball (`arcballVector`, `arcballRotation`,
`orbitArcball`, `orbitTurntable`), camera (`defaultCamera`, `viewOrientation`, `cameraAxes`, `cameraMatrices`,
`fitDistance`, `fitCamera`, `panTarget`, `zoomCamera`, `screenRay`), colour (`colourMap`, `colourMapTexture`,
`colourMapStops`, `luminance`, `valueRange`, `normaliseValues`, `parseColour`, `formatNumber`), buffers
(`vertexNormals`, `triEdgePairs`, `uniqueEdges`, `extractBoundary`, `decimateIndices`, `autoExaggeration`,
`gridMesh`, `buildItem`, `applyDeformed`, `itemEdges`, `itemBounds`, `unionBounds`, `thickLineBuffers`,
`expandIndexed`, `highlightGeometry`), picking (`pickTriangles`, `pickScreen`), text (`dimensionText`,
`describeScene`) and `shaderSource`.

Limits: lines and mesh edges are not depth-sorted against transparent faces; transparent models are drawn
unsorted; the 2-D fallback draws wireframe, lines and points only (faces are not shaded or coloured by value).

## Adapters (`src/geometry/adapters`)

Every adapter exports `id`, `title`, `formats`, `extensions`, `library` (`{ name, version, licence }`),
`sizeBytes` (what a browser downloads on first use), `available()` and
`async parse(input, opts) -> { models, warnings, meta }` with `input = { name, bytes | text }`, the same shape
as the importers, and `sniff(input)`, which says whether the content is theirs. All five are part of this
build. The four that need a runtime also export `load()` (starts the runtime once; `parse` calls it); the
runtime is fetched relative to the module on first use and works on the page, in a module worker and under
Node. `docs/vendor.md` records every vendored file.

| Adapter | Runtime | Formats | Result |
| --- | --- | --- | --- |
| `occt.js` | occt-import-js 0.0.23, 7.7 MB | STEP (AP203 / 214 / 242), IGES, BREP | one `surface` per solid or shell, one group per B-rep face |
| `hdf5.js` | h5wasm 0.10.3, 4.2 MB | NetCDF-4 (CF grids, UGRID), MED, CGNS, Exodus II on NetCDF-4, generic HDF5 | `grid`, `points`, `surface`, `frame`, `volume` |
| `web-ifc.js` | web-ifc 0.0.78, 7.5 MB | IFC2x3 / IFC4 / IFC4x3, all products | one `surface` per IFC class, one group per element |
| `laz.js` | laz-perf 0.0.7, 0.3 MB | LAZ | one `points` model |
| `sqlite.js` | none (reads the SQLite file format itself) | GeoPackage vector layers | `points`, `polyline`, `section` |

`available()` is true when WebAssembly exists and, under Node, the vendored files are on disk. When a runtime
cannot start, `parse` rejects with a message that names the missing piece and what to convert the file to
(`convertTo`).

### occt.js

Options `{ format, quality: 'coarse' | 'normal' | 'fine', linearDeflection [m], angularDeflection [rad],
units (BREP only), upAxis }`. Metres; STEP and IGES declare their unit, BREP is taken as millimetres.

### hdf5.js

The convention is recognised from the object tree, never from the file name; `meta.flavour` says which one was
found and `meta.tree = { entries: [{ path, type, shape?, dtype? }], truncated }` lists the first 2000 objects.

| `meta.flavour` | Recognised by | Models |
| --- | --- | --- |
| `med` | root group `ENS_MAA` | frame, surface, volume |
| `cgns` | root label "Root Node of HDF5 File", a `CGNSLibraryVersion` or a `CGNSBase_t` node | frame, surface, volume |
| `exodus-netcdf4` | datasets `coordx` or `coord` and `connect<n>` | frame, surface, volume |
| `netcdf4-ugrid` | a variable with `cf_role = "mesh_topology"` and `face_node_connectivity` | surface |
| `netcdf4-cf-grid` | a numeric variable of two or more dimensions on x/y or lon/lat axes | grid; points when an axis is uneven |
| `generic` | anything else | none, unless `opts.datasets` names them |

Options: those of the grid, vector and mesh importers (`units`, `scale`, `crs`, `vertical`, `zUnits`,
`projection`, `origin`) and `maxCells`, `variable`, `mesh`, `datasets`.

**CF grids.** The variable is `opts.variable`, else the first of `elevation` (GEBCO, EMODnet), `z` (ETOPO, GMT),
`Band1`, `bathymetry`, `topo`, `depth`, `depth_average` … found on recognisable axes, else one with a height or
depth standard name, else the first gridded variable (a warning lists the others). Axes are the attached
dimension scales; without them, one-dimensional variables of matching length named or attributed as x / y
(`lon`, `lat`, `x`, `y`, `easting`, `northing`, `COLUMNS`, `LINES`, `axis`, `standard_name`, degree units).
Leading dimensions (time, level) are taken at index 0. `scale_factor`, `add_offset`, `_FillValue`,
`missing_value`, `valid_range`, `valid_min` and `valid_max` are applied; `positive = "down"` or a depth name
flips the sign so that z is up. The variable is read in strips of rows and block-averaged while it is read to at
most `opts.maxCells` nodes (default 1 500 000; 0 keeps every node), no-data cells being left out of each
average; `meta.decimation = { factor, source: [nx, ny] }`. Row 0 of the result is the southern row whatever the
storage order. Longitude/latitude grids are placed on the local tangent plane about the grid centre
(`projection: 'none'` keeps degrees). Deflate and shuffle are decoded; plug-in filters are not.

**UGRID.** Node coordinates, a node elevation variable when one exists (depth made negative) and
`face_node_connectivity` with its `start_index` and fill value; quadrilaterals and polygons are triangulated.

**MED.** The first mesh (or `opts.mesh`) and its first computation step. Coordinates and connectivities are
stored component by component and are re-interlaced. Corner nodes are mapped from the MED order to the order of
this application (canonical node k = MED node `map[k]`):

| MED type | Read as | `map` |
| --- | --- | --- |
| `SE2`, `TR3`, `QU4` | line2, tri3, quad4 | identity |
| `TE4` | tet4 | `[0, 2, 1, 3]` |
| `HE8` | hex8 | `[0, 3, 2, 1, 4, 7, 6, 5]` |
| `PE6` | wedge6 | `[0, 2, 1, 3, 5, 4]` |
| `PY5` | pyr5 | `[0, 3, 2, 1, 4]` |

MED numbers the base face of a solid cell with the opposite sense of rotation, so without the mapping every
cell would have a negative volume. The tables agree with files written by the MED 4.1 library through Gmsh 4.15.
Quadratic cells (`SE3`, `TR6`, `TR7`, `QU8`, `QU9`, `T10`, `H20`, `H27`, `P15`, `P18`, `P13`) are read on their
corner nodes with a warning. Families: each group named by a family of cells becomes `sets.cells` /
`sets.elements` / `sets.members`, each group of a family of nodes becomes `sets.nodes`; surface groups and cell
tags are named after the groups (`"a + b"` for a family in several groups, the family name when it has none).
The axis unit (`UNI`) is used when it is a length unit. Not read: polygons, polyhedra, structured MED grids,
fields, further computation steps.

**CGNS.** Every zone of every base, in name order. Unstructured zones: all `Elements_t` sections — `BAR_2/3`,
`TRI_3/6`, `QUAD_4/8/9`, `TETRA_4/10`, `PYRA_5/13/14`, `PENTA_6/15/18`, `HEXA_8/20/27` and `MIXED`; the CGNS
node orders are those of this application except for `HEXA_20` / `HEXA_27`, whose top-face and vertical mid-edge
nodes are exchanged. Structured zones become hex8 (3-D), quad4 (2-D) or line2 (1-D) cells on the vertex
lattice, i running fastest. Each section (or structured zone) is an element set and a surface group; boundary
conditions of unstructured zones given as `PointList` / `PointRange` / `ElementList` / `ElementRange` become
element sets (grid location `FaceCenter`, `CellCenter` …) or node sets (`Vertex`), under the name of the
condition, of its family and of the names listed in that family. `DimensionalUnits` of the base gives the
length unit. Zones are placed side by side: node and element numbers of later zones are offset, and nodes on
zone interfaces are not merged. Not read: `NGON_n` / `NFACE_n` sections, rind planes, cylindrical or spherical
coordinates, boundary conditions on structured index ranges, solutions. ADF files are refused with the
`cgnsconvert -h` advice.

**Exodus II on NetCDF-4.** Element blocks (`BAR`, `TRI`, `QUAD` / `SHELL`, `TET`, `HEX`, `WEDGE`, `PYRAMID`,
linear and quadratic) with their ids and names as sets and tags, and node sets. Side sets are counted, not read.

**Generic HDF5.** Nothing is guessed. Without `opts.datasets` the result has no models and one warning that
lists the datasets. `opts.datasets` takes dataset paths:

| `opts.datasets` | Result |
| --- | --- |
| `{ grid, x?, y?, x0?, y0?, dx?, dy?, nodata? }` | `grid` from a `[ny, nx]` dataset; evenly spaced axis datasets or origin and spacing |
| `{ points }` or `{ x, y, z }` | `points` from an `[n, 3]` dataset or three equally long ones |
| `{ positions, elements: { tri3 \| quad4 \| tet4 \| hex8 \| wedge6 \| pyr5 \| line2: path }, indexBase? }` | frame / surface / volume; connectivity in the order of this application; 1-based indices are detected when they run to n |

Limits: the whole file has to be in memory (its bytes are mounted in the library's in-memory file system
without a copy, and the library reads the datasets it is asked for on top of that); a failed read reports the
innermost line of the HDF5 error stack.

### web-ifc.js

Options `{ units, circleSegments = 12 (4–128), classes: ['IfcPlate', …], openings = false }`. One `surface` per
IFC class with `ifcClass` and `source`; groups carry `name`, `expressId`, `globalId`, `tag`, `color`;
`meta = { schema, units, engine, circleSegments, classes: { IfcBeam: n }, elements: [{ class, name, globalId,
expressId, model, group, triangles }] }`. Placements, the project length unit (SI prefixes and conversion-based
units) and opening subtraction are evaluated by the engine; results are in metres, z up. Opening and space
volumes are left out and counted in a warning unless `openings` is true or their class is asked for. Faces are
tessellated separately, so vertices on shared edges are duplicated (`repair` welds them). Reached through
`importFile(input, { fullGeometry: true })`; without that option IFC files go to `importers/ifc.js`, which gives
frame models with section properties for linear elements.

### laz.js

Options: those of the LAS importer (`units`, `vertical`, `projection`, `origin`) and `maxPoints` (default
2 000 000; 0 keeps all). `readLazHeader(bytes, fileName)` returns the public header and the "laszip encoded"
record (VLR 22204: compressor, version, chunk size, items) without decompressing anything. `parse` decodes the
records in file order — the stream cannot be entered in the middle — keeps every n-th one so that at most
`maxPoints` remain, and hands them to the LAS importer as an uncompressed image with the original header and
VLRs. Scale and offset, `attrs.classification`, `attrs.intensity` (when any is non-zero), `attrs.rgb`, the
coordinate reference system and unit handling are therefore exactly those of a LAS file; `meta` adds
`compression = { compressor, version, chunkSize, items }`, `decoder`, `pointCount` and `thinning`.

Read: point formats 0–3 (pointwise chunked compressor) and 6–8 (layered compressor), fixed and variable chunk
sizes, extra bytes. Refused with the reason: files without the laszip record, the unchunked stream of early
LASzip versions, waveform formats 4, 5, 9 and 10, files of 2 GB or more (the decoder needs the whole file in its
32-bit memory; a file takes about twice its size while it is read). LAZ carries no checksum: damaged bytes
either stop the decoder ("could not be decompressed") or decode into wrong coordinates, in which case points
outside the bounding box of the header are counted in a warning. An uncompressed LAS file named `.laz` is read
directly.

### sqlite.js

GeoPackage 1.0–1.4 feature tables through a read-only reader of the SQLite 3 file format written for this
application (table B-trees, overflow chains, all text encodings; `openSqlite(bytes)`, `gpkgHeader`, `parseWkb`
and `parseColumns` are exported). Options: those of the vector importer (`projection`, `origin`, `utmZone`,
`crs`, `units`, `vertical`, `polygons`) and `layers` (table names), `maxFeatures`. Points become one `points`
model per layer with the attribute columns as per-point attributes, lines one `polyline` each, flat polygons one
`section` each (polygons with varying z: closed polylines); longitude/latitude layers are projected about one
common origin. `meta = { application, sqlite: { pageSize, encoding, pages, journal }, layers: [...], tiles:
[...], crs, units }`. Tile pyramids (imagery, tiled gridded coverage) are listed in `meta.tiles`, not decoded;
`WITHOUT ROWID` tables are reported, not read; a database in WAL mode is read as of its last checkpoint, with a
warning.
