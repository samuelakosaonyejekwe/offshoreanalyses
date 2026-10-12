# Solver interoperability and the panel-method added mass

Two things live here:

* `src/interop/` — writers that turn the app's own models into ready-to-run input for established open-source
  solvers, and readers that bring the results back;
* `src/engine/hydro/bem.js` — a boundary-element solver for the added mass of a hull in the infinite-frequency
  and zero-frequency limits.

Everything is plain ES modules without DOM or network access, so it runs in the browser and under Node 20.

## 1. What was run against real solvers, and what was not

Checked on the development machine with `which` and the usual install locations:

| Solver | Found | Status of the hand-off |
| --- | --- | --- |
| OpenFOAM v2406 (openfoam.com) | yes, `/usr/lib/openfoam/openfoam2406` | **Executed end to end** (wave tank and current case, see §5.3) |
| CalculiX (`ccx`) | no | **Not executed.** The deck follows the documented input format and is re-read by the app's Abaqus importer; it has not been run through `ccx` |
| NEMOH (`preProc`, `solver`, `postProc`) | no | **Not executed.** Input files follow the NEMOH v3.0 manual; the result readers are tested on files laid out as the manual describes, not on real NEMOH output |
| Capytaine (Python package) | no (`import capytaine` fails; `xarray` missing too) | **Not executed.** The generated script passes `python3 -m py_compile` (syntax only); its API calls are written from the Capytaine 2 documentation |
| WAMIT | no (commercial) | **Not executed.** Only the GDF geometry is written; it is re-read by the app's GDF importer |

Consequences to keep in mind until someone runs these cases:

* a keyword or column that a particular solver build rejects would not have been caught for NEMOH, Capytaine,
  WAMIT or CalculiX;
* the conversion of NEMOH / Capytaine excitation phases to the app's convention (complex conjugate, see §3.2)
  follows the documented time conventions and has not been confirmed against a solved benchmark;
* the unit of the phase column of NEMOH's `ExcitationForce.tec` (radians or degrees) is detected from the
  values, with an explicit override;
* the meaning of the per-block number in NEMOH's `CA.dat` / `CM.dat` (period or frequency) is read from the
  header word, with an explicit override.

The registry (`SOLVER_EXPORTS`) carries the same information in the `executed` flag of each entry.

## 2. Common conventions

* Every writer returns `{ files: [{ path, text | bytes }], notes: string[] }`. A path ending in `/` is an empty
  folder. `notes` are meant to be shown to the user: they report what was changed (clipped triangles, reversed
  normals, replaced sections), what is not modelled and how to run the case.
* `packCase(files, { root })` packs a case into one ZIP (stored, no compression; local headers, central
  directory and CRC-32 written by the module; `Allrun`, `Allclean`, `.sh` and `.py` get mode 0755).
* SI units, z up, still-water level at z = 0. Node and element numbers written to solvers are the app's
  0-based indices plus one.
* Invalid input throws `RangeError` naming the argument.

```js
import { SOLVER_EXPORTS, exportCase, packCase } from './src/interop/index.js';
const { files, notes } = exportCase('nemoh', hullSurface, { omegas, headingsDeg, cog });
const zip = packCase(files, { root: 'hull_nemoh' });      // Uint8Array
```

`SOLVER_EXPORTS` entries: `{ id, solver, purpose, kinds, extension, executed, write(model, options),
read(input, options), reads, notes }`; ids `nemoh`, `capytaine`, `wamit-gdf`, `calculix`,
`openfoam-wave-tank`, `openfoam-current`.

## 3. Potential-flow solvers (`src/interop/nemoh.js`)

The wetted mesh is the hull clipped at z = 0 (`wettedPanels` of the BEM module): flat triangles with normals
into the fluid, written as quadrilaterals with a repeated node; faces lying in the waterplane are dropped; an
inward-oriented surface is detected (negative displaced volume) and reversed. No symmetry plane is used.

### 3.1 NEMOH v3.0

`writeNemohCase({ surface, rho, g, depth, omegas, headingsDeg, cog, name })` writes

| File | Content |
| --- | --- |
| `Nemoh.cal` | environment, one body, 6 degrees of freedom and 6 generalised forces (rotations and moments about `cog`), frequency line with the v3 type flag (`1` = rad/s), heading range; IRF, Kochin, free-surface, RAO and QTF output switched off |
| `input_solver.txt` | v3 solver settings (Gauss order 2, LU) |
| `mesh/<name>.dat` | solver-format mesh: `2 0`, node table, `0 0. 0. 0.`, panel table, `0 0 0 0` |
| `results/` | empty folder |
| `ID.dat`, `input.txt` | only read by NEMOH v2; harmless for v3 |

NEMOH takes uniform ranges (count, min, max). Lists that are not uniformly spaced are replaced by the uniform
range with the same count and limits, and a note says so. Deep water is written as depth `0.`.

Not written: `Mesh.cal` (the mesh is already in solver format, so NEMOH's `mesh` pre-processor is not needed),
the hydrostatics / mass files of the `Mechanics` folder (the app keeps mass and hydrostatics itself; RAO
post-processing inside NEMOH is therefore off), QTF input.

### 3.2 Reading NEMOH results

`readNemohResults({ files, rho?, g?, depth?, phaseUnit?, conjugate?, abscissa?, Ainf?, A0? })` returns the
container of `src/engine/hydro/hydrodb.js` (`createHydroDb`), so `raoSolve`, `hydroAt`, `checkHydroDb` and the
Cummins model work on it directly.

* radiation: `RadiationCoefficients.tec` (one zone per radiating mode, columns ω, then (A, B) per force mode),
  or `CM.dat` (added mass) with `CA.dat` (damping);
* excitation: `ExcitationForce.tec` (one zone per heading, columns ω, then (modulus, phase) per mode);
* `rho`, `g`, `depth` from `Nemoh.cal` when it is among the files;
* NEMOH and Capytaine write complex amplitudes for exp(−iωt); the app uses exp(+iωt), so the excitation is
  conjugated (`conjugate: false` switches this off);
* single body with six modes only; multi-body and generalised-mode results are rejected;
* NEMOH's files do not contain the infinite-frequency added mass: pass `Ainf` (for example from
  `addedMassBem`, §6) when the time-domain model needs it.

### 3.3 WAMIT GDF

`writeWamitGdf({ surface, ulen, g, name })` writes a low-order GDF (`ULEN GRAV`, `ISX ISY` = 0 0, `NPAN`, four
vertices per panel, counter-clockwise seen from the fluid). The WAMIT control files are not written. WAMIT
output is read with the existing `hydroDbFromWamit` (`.1`, `.3`) using the same `ULEN`.

### 3.4 Capytaine

`writeCapytaineScript({ … same arguments … })` writes `run_capytaine.py` and `mesh/<name>.dat`. The script
loads the mesh, solves all radiation and diffraction problems with `BEMSolver.fill_dataset`, saves the complete
dataset as NetCDF and writes `<name>_hydro.csv` in the app's CSV schema (`kind,omega,headingDeg,i,j,re,im`),
converted to exp(+iωt), with `# rho=`, `# g=`, `# depth=` comment lines and, in deep water, the
infinite-frequency added mass (`AINF` rows). `readCapytaineCsv(text)` reads that CSV. The NetCDF file itself is
not read by this module.

## 4. CalculiX (`src/interop/calculix.js`)

Target: CalculiX CrunchiX 2.20 / 2.21, Abaqus-style `.inp`.

`writeCalculixCase(model, { material, supports, loads, analysis, steps, elementOrder, thickness,
genericSection, name })`

| Model | Elements | Section |
| --- | --- | --- |
| `frame` | B31, or B32 with generated mid-side nodes (`elementOrder: 2`) | `*BEAM SECTION` PIPE for tubes, CIRC / RECT / BOX when the imported section carries that shape; orientation line per section |
| `surface` | S3 | `*SHELL SECTION` with the group's thickness or `thickness` |
| `volume` | C3D4, C3D10, C3D8, C3D20, C3D6 | `*SOLID SECTION`; 5-node pyramids are split into two C3D4 |

Material: `*ELASTIC`, `*DENSITY`, `*PLASTIC` (perfectly plastic, bilinear from a tangent modulus, or a table).
Supports: `*NSET` + `*BOUNDARY` (`'fixed'`, `'pinned'`, or a list of degrees of freedom; non-zero prescribed
values allowed). Loads: nodal forces and moments (`*CLOAD`), gravity and pressure (`*DLOAD`), frame line loads
lumped to the member ends. Step: `*STATIC`, `*STATIC` with `NLGEOM`, `*BUCKLE`, `*FREQUENCY`, `*DYNAMIC, DIRECT`.
Output: `*NODE FILE` (U, RF), `*EL FILE` (S, E, PEEQ), `*NODE PRINT`, `*EL PRINT`; `OUTPUT=2D` for beams and
shells so that results return on the original nodes.

Limits, stated plainly:

* one material for the whole model, one step;
* generic beam sections `{ A, Iy, Iz, J }` have no CalculiX equivalent. By default they are written as the
  rectangle with the same `Iy` and `Iz` (area, torsion and mass then differ; the note gives the area ratio).
  `genericSection: 'general'` writes `*BEAM GENERAL SECTION` instead, which Abaqus reads and CalculiX does not;
* line loads are lumped (q·L/2 per end node), without the fixed-end moments of a consistent load vector;
* a dynamic step applies the loads at full value from t = 0; no `*AMPLITUDE` is written;
* contact, springs, point masses, member releases and offsets, multi-point constraints are not written.

Readers:

* `readCalculixFrd(text, { nodeCount, positions, step, scale })` — ASCII `.frd`: every nodal block (DISP,
  STRESS, …) of every step, parsed by fixed columns, mapped to node id − 1; returns `displacement`, `stress`,
  `vonMises`, and for the 3-D viewer `nodeValues` (von Mises stress when available, else displacement
  magnitude), `valueLabel`, `valueUnit`, `deformed`. Binary `.frd` is not read.
* `readCalculixDat(text, { nodeCount, positions, scale })` — printed tables: displacements, reaction forces and
  their total, integration-point stresses averaged to one von Mises value per element.

Both readers are tested on files written in the test in the documented layout, not on files produced by `ccx`.

## 5. OpenFOAM (`src/interop/openfoam.js`)

Target: **OpenFOAM v2406** from openfoam.com (v2312 uses the same dictionaries). The wave boundary conditions
(`waveVelocity`, `waveAlpha`, `constant/waveProperties` with `waveModel`) exist only in that line. For the
Foundation releases (v11, v12) `blockMeshDict`, `snappyHexMeshDict`, `decomposeParDict`, `setFieldsDict`,
`constant/g` and the `forces` function object carry over; `transportProperties` / `turbulenceProperties`
become `physicalProperties` / `momentumTransport`, the wave set-up uses that line's own wave library, and the
solver is started with `foamRun`. No Foundation-specific files are written.

### 5.1 Wave tank

`writeOpenfoamWaveTank({ surface, domain: { length, width, height, waterDepth, xMin?, yMin? }, wave: { H, T,
theory, rampTime? }, current, cellSize, refinement: { structure, box, freeSurface, boxMargin }, endTime,
turbulence, processors, stl, cofr, writeInterval, fluid, name })`

Files: `constant/triSurface/<name>.stl` (binary or ASCII, one solid named `structure`), `system/blockMeshDict`,
`snappyHexMeshDict` (refinement on the structure, in a box around it and in the band |z| ≤ H), `controlDict`
(with the `forces` function object on patch `structure`), `fvSchemes`, `fvSolution`, `decomposeParDict`,
`setFieldsDict`, `constant/transportProperties`, `turbulenceProperties` plus the other variant as
`turbulenceProperties.kOmegaSST` or `.laminar`, `g`, `waveProperties`, `0/alpha.water` (+ `.orig`), `U`,
`p_rgh`, `k`, `omega`, `nut`, `Allrun`, `Allclean`.

Patches: `inlet` (wave generation), `outlet` (wave absorption), `sides` (slip), `bottom` (wall), `atmosphere`,
`structure` (wall). Wave theories: `airy` → StokesI, `stokes2` → StokesII, `stokes5` → StokesV, `cnoidal`.

Not modelled: a current combined with the waves (the wave boundary conditions of this release do not superpose
one; `current` is reported in the notes and ignored), irregular or oblique waves, moving or floating bodies
(no dynamic mesh), boundary-layer cells (`addLayers` is off). The outlet absorbs with shallow-water theory.

### 5.2 Current past a submerged body

`writeOpenfoamCurrentCase({ surface, domain: { length, width, height }, current, cellSize, refinement,
iterations, turbulence, turbulenceIntensity, fluid, processors, stl, cofr, name })` — `simpleFoam`, k-ω SST
with wall functions (or laminar), velocity inlet, pressure outlet, slip on the four lateral faces, no free
surface. Inlet k = 1.5 (I·U)², ω = √k / (0.09^¼ · 0.07 · D) with D the cross-flow size of the body.

### 5.3 Runs performed (OpenFOAM v2406, build `_1653fa08-20260127`)

| Case | Set-up | Result |
| --- | --- | --- |
| Wave tank, serial | vertical cylinder D = 0.4 m through bottom and surface; tank 6 × 1.6 × 0.7 m, depth 0.4 m; StokesII, H = 0.06 m, T = 1.2 s; 0.1 m cells, one refinement level; `endTime` 1.2 s | `blockMesh`, `snappyHexMesh` (24 472 cells, `checkMesh`: Mesh OK), `setFields`, `interFoam` to t = 1.2 s in 57 steps (34 s CPU). The wave model reported reference depth 0.4 m and wave length 1.9358 m; linear dispersion for T = 1.2 s, d = 0.4 m gives 1.9358 m. Water fraction 0.5714 = 0.4 / 0.7 |
| Wave tank, from the ZIP, 2 processes | same case packed with `packCase`, extracted with `unzip`, started with `./Allrun` | all steps including `decomposePar`, parallel `interFoam`, `reconstructPar` completed; F_x at t = 1.2 s: −0.69870 N (serial run: −0.69870 N) |
| Wave tank, k-ω SST variant | same mesh, `turbulenceProperties.kOmegaSST` copied over, `endTime` 0.1 s | `interFoam` selected `kOmegaSST` and ran 11 steps to the end |
| Current case | sphere D = 1 m in a 12 × 6 × 6 m box, U = 1 m/s, 0.5 m cells, two levels on the body; `./Allrun`, 60 iterations | 11 779 cells, Mesh OK, `simpleFoam` completed; F_x = 147.4 N at iteration 60 (C_d ≈ 0.37 with ρ = 1025 kg/m³) |

`force.dat` and `moment.dat` of these runs were read with `readOpenfoamForces`; lines of the first run are the
fixture of the reader's test.

What these runs show and what they do not: the dictionaries are accepted by the mesher and the solvers, the
boundary conditions are consistent, and the scripts work after a ZIP round trip. They are **smoke tests on very
coarse meshes** (1.2 cells per wave height; a 60-iteration steady run that is not converged). The forces above
are not validated numbers. Mesh and time-step convergence, y+ and reflection checks remain with the user.

### 5.4 Forces reader

`readOpenfoamForces(text, { momentText })` returns `{ format, cofr, t, Fx, Fy, Fz, Mx, My, Mz, pressure,
viscous, n }` from either layout:

* columns (`force.dat` and `moment.dat`, OpenFOAM v2106 and later) — tested on real v2406 output;
* bracketed (`forces.dat`, OpenFOAM up to v2012 and the Foundation releases; with or without the porous
  group) — tested on a file written by hand in that layout, since no such release is installed.

## 6. Panel-method added mass (`src/engine/hydro/bem.js`)

`addedMassBem({ surface, rho, limit, cog, depth?, formulation?, solver?, imagePairs? })` →
`{ A (6×6, SI, symmetric), panels, conditionNumber, symmetryError, area, volume, limit, formulation, solver }`
and `wettedPanels(surface)`.

| `limit` | Free-surface condition | Image of the hull about z = 0 |
| --- | --- | --- |
| `'infinite'` (ω → ∞) | φ = 0 | opposite sign |
| `'zero'` (ω → 0) | ∂φ/∂z = 0 (rigid lid) | equal sign |
| `'unbounded'` | none | none |

A∞ is the constant added mass of the Cummins time-domain equation.

**This is not a radiation/diffraction solver.** It gives no frequency dependence, no wave damping and no
excitation forces; for those use the hand-off of §3.

Two discretisations on the flat triangles of the wetted hull, collocation at the centroids, panel integrals in
closed form:

* `formulation: 'potential'` (default) — Green's identity with piecewise-constant potential and the known
  source strength (Morino);
* `formulation: 'source'` — constant-strength source panels (Hess & Smith), kept as an independent cross-check.

The default is the potential formulation because it is markedly more accurate on the same mesh:

| Deeply submerged sphere, A₁₁ / (½ρV) | 80 panels | 320 | 1280 | observed order in panel size |
| --- | --- | --- | --- | --- |
| potential | 0.9216 | 0.9797 | 0.9949 | 1.95, 2.00 |
| source | 1.0116 | 1.0374 | 1.0267 | not monotone against the true sphere; 1.1, 1.05 against ½ρV of the faceted body |

Other benchmarks in `tests/engine/hydro.bem.test.js`:

* floating hemisphere, 656 panels: A₃₃(∞) = ρπa³/3 exactly (the double body with the opposite-sign image is a
  rigid sphere translating vertically) — computed 0.9949 of it; A₁₁(0) = ρπa³/3 — computed 0.9949; A₃₃(0) from
  the Legendre series derived in the test, 0.8310 · ρ(2/3)πa³ (Hulme 1982) — computed 0.8255 · ρ(2/3)πa³;
  with 2592 panels A₃₃(∞) is 0.9987 of the exact value;
* horizontal cylinder, deeply submerged: sway and heave 0.966 ρπr²L at L/D = 20 and lower at L/D = 5 (end
  effects), approaching the strip value from below;
* the returned matrix is exactly symmetric and positive definite; the reference point transfers moments
  correctly.

Statements about accuracy that should not be over-read:

* the raw collocation matrix is not exactly symmetric. The returned `A` is its symmetric part and
  `symmetryError` reports the asymmetry that was removed (0.6 % for a 448-panel raked barge about an offset
  reference point; 0.16 % for a 2398-panel box barge);
* bodies with sharp convex edges converge more slowly than the sphere: for a 20 × 6 m barge at 2 m draught the
  heave added mass changes by 2 % and surge by 5 % between 448 and 1792 panels. Run two meshes;
* finite depth is implemented for `limit: 'infinite'` only (alternating image series about free surface and
  seabed, truncated after `imagePairs` pairs with averaging of the last two partial sums; 12 pairs changed the
  barge result at 1 m under-keel clearance by less than 0.1 % against 24). It has no closed-form benchmark in
  the tests: it is checked for the deep-water limit, for monotonic growth as the clearance closes, and against
  the source formulation. `limit: 'zero'` with a finite depth throws, because the rigid-lid heave added mass
  of a surface-piercing body is unbounded between two parallel walls;
* thin plates and appendages thinner than a panel are outside the range of a constant-panel method.

Cost: dense storage (16 n² bytes), GMRES with an LU fallback (`solver: 'lu'` forces LU and returns the
condition estimate). Measured under Node 20 on a laptop CPU (Core i5-7300U): 656 panels 0.25 s, 1280 panels
1.6 s, 2592 panels 3.1 s. The limit is 6000 wetted panels.

## 7. Tests

```
node --test tests/interop/ tests/engine/hydro.bem.test.js
```

* NEMOH mesh and GDF: panel and node counts, wetted area and volume, normals into the fluid for a box and a
  cylinder; both files are read back by the app's panel importers;
* NEMOH and Capytaine result readers: fixtures in the documented layouts → database → `raoSolve`;
* CalculiX: a cantilever deck is re-read by the app's Abaqus importer to the same nodes, elements, section,
  material and sets; shells and all solid types likewise; `.frd` and `.dat` readers on fixed-format fixtures;
* OpenFOAM: folder contents, patch consistency across fields, balanced dictionaries, STL read back by the STL
  importer, both `forces` layouts;
* ZIP: read back by `src/geometry/zip.js`, CRC-32 check value, executable bits.
