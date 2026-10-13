# Solver interoperability and the panel-method added mass

Two things live here:

* `src/interop/` — writers that turn the app's own models into ready-to-run input for established open-source
  solvers, and readers that bring the results back;
* `src/engine/hydro/bem.js` — a boundary-element solver for the added mass of a hull in the infinite-frequency
  and zero-frequency limits.

Everything is plain ES modules without DOM or network access, so it runs in the browser and under Node 20.

## 1. What was run against real solvers, and what was not

State of 2026-10-13. "Executed" means: a case written by the application's writer was unpacked, run by the solver
without editing, and its results were read back by the application's reader. The open solvers were installed
outside the repository (virtual environment, conda-forge package, build from source); nothing of them is kept
here except the small result files under `tests/interop/fixtures/` (provenance in its `README.md`).

| Solver | Version run | Status of the hand-off |
| --- | --- | --- |
| OpenFOAM | v2406 (openfoam.com) | **Executed** (wave tank and current case, §5.3) |
| Capytaine | 3.0.0 (Python 3.12) | **Executed**: barge 20 × 8 × 2 m (768 panels, 12 frequencies × 3 headings, and 5 × 3), vertical cylinder D = 8 m, 40 m draught (1128 panels, 12 frequencies), hemisphere a = 5 m (2592 panels); deep water; CSV read back |
| NEMOH | v3, master `dcc8113` of 2026-06-24, built from source with gfortran | **Executed**: the barge, 5 frequencies × 3 headings, `preProc`, `solver`, `postProc`; `RadiationCoefficients.tec`, `ExcitationForce.tec`, `CM.dat`, `CA.dat` read back |
| CalculiX (`ccx`) | 2.23 (conda-forge) | **Executed**: static, `*FREQUENCY` and `*BUCKLE` steps of a cantilever and a portal frame (B32R with PIPE; B31 and B32 with RECT), static step of a plate (S6, S3); `.dat` and `.frd` read back |
| WAMIT | — (proprietary) | **Not executed.** The `.gdf` files of the three hulls were loaded with an independent GDF reader (Capytaine `load_mesh`): panel count, vertex count after welding, wetted area, displaced volume and normal direction equal those of the source mesh (§3.3) |

The registry (`SOLVER_EXPORTS`) carries the same information in the `executed` flag and the `notes` of each entry.

### 1.1 Defects the runs exposed, and their fixes

| Where | What the solver showed | Fix |
| --- | --- | --- |
| NEMOH writer | `preProc` stopped with "Cannot open file './barge.dat'": the mesh named in `Nemoh.cal` is opened relative to the case folder, and `mesh/` is where preProc writes its own files | the mesh is written to the case folder as `<name>.dat` |
| NEMOH reader, `CM.dat` / `CA.dat` | each line of a block belongs to one *radiating* degree of freedom, so the block was read transposed (A15 and A51 exchanged: 229 t·m against 256 t·m for the barge); the frequency is printed with four decimals only, which made the consistency check with `ExcitationForce.tec` fail for frequencies such as 2π/7 | block transposed; within the rounding the frequencies of `ExcitationForce.tec` are used |
| Capytaine script | with Capytaine 3.0.0 `separate_complex_values(ds).to_netcdf(...)` raised a `TypeError` (categorical coordinate) after all problems had been solved, and no CSV was written | the CSV is written first; NetCDF goes through `cpt.export_dataset` when it exists, and a failure there is reported, not fatal |
| Capytaine script | only the infinite-frequency added mass was requested | `A0` rows (omega = 0) added; `readCapytaineCsv` returns `db.A0` |
| CalculiX writer, tubes | "*BEAM SECTION of type PIPE can only be used for B32R elements": every deck with a tubular member was rejected | a model with a PIPE or BOX section is written with B32R elements and generated mid-side nodes |
| CalculiX writer, beam orientation | without an orientation of its own, a member got the local 1-axis (0, 0, −1): a rectangle b × h lay on its side (tip deflection 59.4 mm instead of 14.9 mm) and `Iy` of a generic section acted in the horizontal plane | default 1-axis = local y of `frame3d.js` (horizontal normal to the member; (0, 1, 0) for vertical members) |
| CalculiX writer, shells | S3 shells gave 0.12 of the plate-theory deflection (they are expanded to wedges that lock in bending) | S6 shells with generated mid-side nodes are the default for surfaces; `elementOrder: 1` still writes S3, with a warning note |
| CalculiX writer, shell pressure | a positive `P` moved the plate along its normal, opposite to the documented meaning of a positive pressure | the pressure of shells is written with the opposite sign |
| CalculiX writer, supports | a supported edge of quadratic elements was free at its generated mid-side node | mid-side nodes whose two end nodes are in a support set are added to it (beams, and shell edges that are boundary edges or interior lines) |
| CalculiX `.dat` reader | `*FREQUENCY` and `*BUCKLE` steps print no nodal table: the reader threw; shell stress rows end with the name of the expanded element and the whole table was dropped | eigenvalues, natural frequencies and buckling factors are returned (`eigenvalues`, `frequencies`, `bucklingFactors`); the trailing name is ignored |

### 1.2 The phase convention, settled

NEMOH and Capytaine write complex amplitudes for exp(−iωt). The application uses exp(+iωt) with the heading as
the propagation direction. Conjugating the excitation (done by `readNemohResults` and inside the Capytaine script)
is correct; it was confirmed on the solved barge in the longest waves, where the answer is known:

* heave: F₃ → ρ g A_wp e^(−kT) − ω² A₃₃ + iω B₃₃. At ω = 0.1 rad/s Capytaine gives 1 596 404 + 129.0 i N/m after
  conversion; ρ g A_wp e^(−kT) − ω² A₃₃ = 1 596 421 and ω B₃₃ = 129.4. The imaginary part is positive — the force
  leads the elevation by the radiation-damping term; without the conjugation it would lag. NEMOH at 0.25 rad/s:
  Im F₃ = 4651, ω B₃₃ = 4653;
* surge, wave towards +x: F₁ → + iω² (ρV + A₁₁): 3911 i against 3918; the sign reverses for 180°;
* the application's own long-wave model (`floaterRAO`) puts the surge force of the cylinder at +90.0°, converted
  Capytaine at +89.99° (25 s) to +69° (4 s, diffraction);
* NEMOH and Capytaine, read through the two independent readers, agree within 0.05 % (added mass), 0.1 %
  (damping) and 0.3 % (real and imaginary parts of the excitation) on the same mesh.

NEMOH v3 writes the phase column of `ExcitationForce.tec` in radians and rad/s in `CM.dat` / `CA.dat`.

### 1.3 What still was not executed

* WAMIT: neither the solver nor a file written by it (`hydroDbFromWamit` is tested on files laid out as the manual
  describes);
* NEMOH and Capytaine in finite depth and with more than one body; NEMOH's QTF module; Capytaine 2.x (the
  script keeps the 2.x NetCDF call as a fallback, which was not run);
* CalculiX: solid (volume) models, `*PLASTIC`, `NLGEOM`, `*DYNAMIC`, gravity and line loads, BOX and CIRC
  sections, and generic sections written as equivalent rectangles;
* OpenFOAM Foundation releases (§5).

Two things observed and left as they are: the total of the reaction forces that CalculiX prints for the supported
nodes of an expanded shell did not equal the applied load (23.3 kN of 30 kN for the plate), and with `OUTPUT=2D`
the nodal stresses of shells are mid-surface values. Use the displacements, or request section forces by hand.
The application's Abaqus importer does not know the element type S6; decks of surface models are read back by it
only when written with `elementOrder: 1`.

### 1.4 The built-in models against these solvers

`src/spec/cross-solver.json` holds every number; the application displays it. Summary (relative difference of the
built-in value from the reference):

| Built-in model | Reference | Result |
| --- | --- | --- |
| `addedMassBem`, ω → ∞ and ω → 0, barge / cylinder / hemisphere, six diagonal terms | Capytaine at ω = inf and 0, same meshes | translations within 3.5 % (largest: cylinder surge), rotations within 5.7 % (barge roll). Hemisphere heave at ω → ∞: built-in −0.13 %, Capytaine +1.5 % from the exact ρπa³/3 |
| `boxBargeRAO` (method `'panel'`), barges 20 × 8 × 2 m, 100 × 20 × 5 m and 80 × 24 × 3 m, head and beam seas | Capytaine + `raoSolve`, centre of gravity on the waterline, dry gyradii 0.35 B and 0.25 L, the same 5 % of critical added roll damping | heave and pitch within 5 % at every compared period (4–25 s, 5–25 s, 4–20 s) except heave in beam seas of the third barge at 4 s (−7.7 %, wavelength = 1.04 B); within 1.5 % from 1.15 × the heave natural period upward. Roll within 5 % from 5 s, 7 s and 4 s respectively (elsewhere: +6.3 % at 4 s for the first barge, +3.1 % against its 1728-panel mesh; +6 % at 6 s and +11 % at 5 s for the second, where the roll is 0.02 and 0.006 rad/m) |
| `boxBargeRAO` (method `'closed-form'`), barge 20 × 8 × 2 m | same | heave within 5 % for T ≥ 8 s (twice the heave natural period of the model), −16 % at 6 s, −31 % at 5 s in head seas and −15 %, −33 %, −42 % at 6, 5, 4 s in beam seas; pitch within 4 % for T ≥ 6 s, −19 % at 5 s, −63 % at 4 s; roll within 5 % for T ≥ 7 s, +8.6 % at 6 s, −8.6 % at 5 s, −15 % at 4 s |
| `floaterRAO`, spar D = 8 m, draught 40 m, hemisphere end mass, free-end loss, Haskind damping | same, no added damping in either | surge force within 5 % for T ≥ 6 s (+6.4 % at 5 s, +47 % at 4 s); heave force and heave RAO within 3.1 % for T ≥ 5 s, the heave resonance (13.1 s) included; surge RAO within 5 % at 5–12 s and 25 s, pitch RAO at 5–12 s and 20–25 s; around the pitch resonance (14.80 s against 14.63 s, radiation damping only) −17 … +24 % |
| `analyzeFrame` | CalculiX `*STATIC` | cantilevers within 1.4 %; portal sway +7.5 % (CalculiX joints are stiff blocks of the tube size) |
| `modalAnalysis` | CalculiX `*FREQUENCY` | cantilever within 1.9 %, portal within 3.5 % (six modes each) |
| `bucklingAnalysis` | CalculiX `*BUCKLE` | cantilever within 1.5 %, portal within 5.6 % (four modes each) |
| `navierPlate` | CalculiX S6 shells | centre deflection within 0.7 % |

The regression tests keep the built-in models within these recorded differences.

What the comparison changed in the floater models (`src/engine/hydro/floater.js`):

| Difference found | Cause | Now |
| --- | --- | --- |
| barge heave in beam seas tended to 1 + kT | the published closed form approximates the depth decay of the bottom pressure by exp(−\|k cos β\| T), which is 1 in beam seas | exp(−kT) and the breadth factor sin(½k_yB)/(½k_yB): the exact integral of the undisturbed pressure over the bottom at any heading; unchanged in head seas |
| barge roll 3–12 % high from 6 s, two to three times too high at 4–5 s | Froude–Krylov moment alone (the diffraction moment lowers it by 1.5 % at 25 s, 9 % at 10 s, 38 % at 5 s), added inertia only through the gyradius given, no radiation damping; the comparison also gave the model 0.4 B including added inertia against 0.35 B dry plus the solver's added inertia | relative-motion diffraction moment −ω² A44 × (wave slope at the keel); strip added inertia A44 from DNV-RP-C205 Table D-1 (8 % above Capytaine's); radiation damping from the Haskind relation with the model's own moment (within 12 % of Capytaine's B44 from 5 s); `kxxDry` input |
| barge heave 15 % low at 1.5 × the heave natural period and 30–40 % low near it, pitch 20–60 % low below 1.25 ×, roll within 15 % below 1.3 × the roll natural period: for ship-size draughts the whole range of wave periods of 5–15 s | the closed form takes the heave added mass equal to the displacement and a two-dimensional sectional damping; the roll model had a strip added inertia and no sway coupling | `boxBargeRAO` solves the radiation problem of the box with a three-dimensional panel method and the free-surface Green function (`src/engine/hydro/radiation.js`): frequency-dependent added mass, wave damping from the radiated energy, exciting forces from the Haskind relation, sway–roll and surge–pitch coupling. Differences from Capytaine as in the table above; the closed form remains as `method: 'closed-form'` (and with forward speed) |
| spar heave force +7 % at 10 s, +24 % at 7 s | wave kinematics sampled at the strip centres: the heave force of a deep column is a small difference, and the sampling error was of its size (10 points of the 24 % at 7 s); the rest is the half-disc end added mass, 65 % of the solver's | kinematics integrated exactly along each strip; end option `'hemisphere'` (2/3 πρa³, 2 % above Capytaine) |
| spar surge and heave RAOs ±50 % near resonance | 2 % of critical damping assumed where the solver has radiation damping only; heave natural period 12.96 s against 13.10 s (end added mass); pitch added inertia 14 % high (no end effect on the strip added mass), pitch natural period 15.10 s against 14.63 s | radiation damping from the Haskind relation with the model's own excitation (within 8 % of Capytaine's B11, B33, B55 for T ≥ 10 s, whose own coefficients satisfy the relation to 2 % in surge and pitch and 6 % in heave); `free1`/`free2` take the transverse added mass lost around a free end (DNV-RP-C205 Table D-2) off at that end: pitch added inertia 2 % high, natural period 14.80 s |

One property of the reference itself: on the 768-panel barge the sway–roll added-mass terms of Capytaine are not
symmetric (A₂₄ = −10 700, A₄₂ = −27 200 kg·m at 25 s). Solving with A₂₄ on both sides raises the reference roll by
3–5 % at 7–25 s, which is the size of the remaining difference there.

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
| `<name>.dat` | solver-format mesh in the case folder (preProc opens it relative to that folder and writes its own files to `mesh/`): `2 0`, node table, `0 0. 0. 0.`, panel table, `0 0 0 0` |
| `results/` | empty folder |
| `ID.dat`, `input.txt` | only read by NEMOH v2; harmless for v3 |

Run: `preProc <case>`, `solver <case>`, `postProc <case>`. NEMOH takes uniform ranges (count, min, max). Lists that are not uniformly spaced are replaced by the uniform
range with the same count and limits, and a note says so. Deep water is written as depth `0.`.

Not written: `Mesh.cal` (the mesh is already in solver format, so NEMOH's `mesh` pre-processor is not needed),
the hydrostatics / mass files of the `Mechanics` folder (the app keeps mass and hydrostatics itself; RAO
post-processing inside NEMOH is therefore off), QTF input.

### 3.2 Reading NEMOH results

`readNemohResults({ files, rho?, g?, depth?, phaseUnit?, conjugate?, abscissa?, Ainf?, A0? })` returns the
container of `src/engine/hydro/hydrodb.js` (`createHydroDb`), so `raoSolve`, `hydroAt`, `checkHydroDb` and the
Cummins model work on it directly.

* radiation: `RadiationCoefficients.tec` (one zone per radiating mode, columns ω, then (A, B) per force mode),
  or `CM.dat` (added mass) with `CA.dat` (damping) — blocks with one line per radiating mode, frequency in rad/s
  with four decimals (the seven-digit frequencies of `ExcitationForce.tec` replace them when that file is there);
* excitation: `ExcitationForce.tec` (one zone per heading, columns ω, then (modulus, phase in radians) per mode);
* `rho`, `g`, `depth` from `Nemoh.cal` when it is among the files;
* NEMOH and Capytaine write complex amplitudes for exp(−iωt); the app uses exp(+iωt), so the excitation is
  conjugated (`conjugate: false` switches this off) — confirmed on solved cases, §1.2;
* single body with six modes only; multi-body and generalised-mode results are rejected;
* NEMOH's files do not contain the infinite-frequency added mass: pass `Ainf` (for example from
  `addedMassBem`, §6) when the time-domain model needs it.

### 3.3 WAMIT GDF

`writeWamitGdf({ surface, ulen, g, name })` writes a low-order GDF (`ULEN GRAV`, `ISX ISY` = 0 0, `NPAN`, four
vertices per panel, counter-clockwise seen from the fluid). The WAMIT control files are not written. WAMIT
output is read with the existing `hydroDbFromWamit` (`.1`, `.3`) using the same `ULEN`.

Checked with an independent reader (Capytaine 3.0.0 `load_mesh(..., file_format="gdf")`), since WAMIT is not
available:

| Hull | Panels (app / reader) | Vertices | Wetted area, m² (app / reader) | Volume, m³ (app / reader) |
| --- | --- | --- | --- | --- |
| barge 20 × 8 × 2 m | 768 / 768 | 413 / 413 | 272 / 272.0 | 320 / 320.0 |
| cylinder D = 8 m, 40 m draught | 1128 / 1128 | 577 / 577 | 1052.1344 / 1052.1344 | 1987.7302 / 1987.7302 |
| hemisphere a = 5 m | 2592 / 2592 | 1345 / 1345 | 156.89192 / 156.89192 | 261.23368 / 261.23368 |

The volume is ∫ z n_z dS with the reader's own normals (right-hand rule on the vertex order); it is positive and
equals ∫ x n_x dS, so the normals point into the fluid on every face family.

### 3.4 Capytaine

`writeCapytaineScript({ … same arguments … })` writes `run_capytaine.py` and `mesh/<name>.dat`. The script
loads the mesh, solves all radiation and diffraction problems with `BEMSolver.fill_dataset` and writes
`<name>_hydro.csv` in the app's CSV schema (`kind,omega,headingDeg,i,j,re,im`), converted to exp(+iωt), with
`# solver=`, `# rho=`, `# g=`, `# depth=` comment lines and, in deep water, the added mass at ω = ∞ and ω = 0
(`AINF` and `A0` rows). Afterwards it saves the complete dataset as NetCDF (`cpt.export_dataset`, or the 2.x
call when that function is missing); if that step fails the CSV is already complete. `readCapytaineCsv(text)`
reads the CSV. The NetCDF file itself is not read by this module. Run with Capytaine 3.0.0 (§1).

## 4. CalculiX (`src/interop/calculix.js`)

Target: CalculiX CrunchiX 2.20 – 2.23, Abaqus-style `.inp`; run with 2.23 (§1).

`writeCalculixCase(model, { material, supports, loads, analysis, steps, elementOrder, thickness,
genericSection, name })`

| Model | Elements | Section |
| --- | --- | --- |
| `frame` | B31, or B32 with generated mid-side nodes (`elementOrder: 2`); B32R for the whole model when it holds a PIPE or BOX section, because CalculiX has these sections for B32R only | `*BEAM SECTION` PIPE for tubes, CIRC / RECT / BOX when the section carries that shape; orientation line per section: the member's own, else the horizontal normal to the member (local y of `frame3d.js`), so that h and `Iy` act in the vertical plane |
| `surface` | S6 with generated mid-side nodes (default), S3 with `elementOrder: 1` (far too stiff in bending in CalculiX, §1.1) | `*SHELL SECTION` with the group's thickness or `thickness` |
| `volume` | C3D4, C3D10, C3D8, C3D20, C3D6 | `*SOLID SECTION`; 5-node pyramids are split into two C3D4 |

Material: `*ELASTIC`, `*DENSITY`, `*PLASTIC` (perfectly plastic, bilinear from a tangent modulus, or a table).
Supports: `*NSET` + `*BOUNDARY` (`'fixed'`, `'pinned'`, or a list of degrees of freedom; non-zero prescribed
values allowed); the generated mid-side node of an edge whose two ends are supported is supported with them. Loads: nodal forces and moments (`*CLOAD`), gravity and pressure (`*DLOAD`; a positive pressure on a shell acts against the element normal and is written
with the opposite sign, because CalculiX applies a positive `P` along the normal), frame line loads
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
  their total, integration-point stresses averaged to one von Mises value per element; from a `*FREQUENCY` step
  `eigenvalues` (rad²/s²) and `frequencies` (Hz), from a `*BUCKLE` step `bucklingFactors`.

In the `.frd` file of a frequency step every mode is one result set whose `time` is the frequency in Hz; in a
buckling step the first set is the pre-load state (`time` 0) and the following ones carry the buckling factor.
Both readers are tested on files produced by `ccx` 2.23 (`tests/interop/fixtures/calculix/`) as well as on files
written in the tests in the documented layout.

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

Both limits were compared with Capytaine 3.0.0 (radiation problems at ω = inf and ω = 0) on identical meshes of a
barge, a vertical cylinder and a hemisphere: translational terms within 3.5 %, rotational terms within 5.7 %
(§1.4, `src/spec/cross-solver.json`).

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
* CalculiX: a cantilever deck is re-read by the app's Abaqus importer to the same elements, section, material and
  sets; S3 shells and all solid types likewise; `.frd` and `.dat` readers on fixed-format fixtures;
* OpenFOAM: folder contents, patch consistency across fields, balanced dictionaries, STL read back by the STL
  importer, both `forces` layouts;
* ZIP: read back by `src/geometry/zip.js`, CRC-32 check value, executable bits;
* `cross-solver.test.js`, on the stored solver output of §1 (no solver needed): the readers return the recorded
  values from real Capytaine, NEMOH and CalculiX files; the phase convention checks of §1.2; NEMOH against
  Capytaine; every built-in model stays within the difference recorded in `src/spec/cross-solver.json` (a change
  that moves a result further from the reference fails); the writer properties that the runs showed to be
  necessary (B32R for tubes, beam orientation, S6, pressure sign, supported mid-side nodes, mesh location for
  NEMOH, CSV before NetCDF); the GDF measurements of the independent reader.
