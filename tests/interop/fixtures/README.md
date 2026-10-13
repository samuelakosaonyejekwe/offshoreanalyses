# Solver output used as reference

Every file in this folder was written by an external solver from a case that the application exported, on
2026-10-13 (Linux x86-64). Nothing here was edited by hand except where "excerpt" says so. The tests in
`../cross-solver.test.js` read these files with the application's readers; the solvers are not needed to run
them. The numbers derived from them are in `src/spec/cross-solver.json`.

Solvers

| Solver | Version | How it was obtained |
| --- | --- | --- |
| Capytaine | 3.0.0 (Python 3.12, numpy 2.5.4, xarray 2026.9.0) | `pip install capytaine` in a virtual environment |
| NEMOH | v3, master branch, commit `dcc8113` of 2026-06-24 (after tag v3.0.3) | built from https://gitlab.com/lheea/Nemoh with CMake and gfortran, LAPACK/BLAS from conda-forge |
| CalculiX CrunchiX (`ccx`) | 2.23 | conda-forge package `calculix` |

Common data: sea water 1025 kg/m³, g = 9.80665 m/s², deep water, rotations and moments about the origin on the
waterline, z up. Steel E = 210 GPa, nu = 0.3, 7850 kg/m³.

## `capytaine/` — CSV written by the exported `run_capytaine.py` (`writeCapytaineScript`)

Hulls from `src/mesh/generate.js`, clipped at the waterline by the writer:

| File | Hull | Panels | Frequencies and headings |
| --- | --- | --- | --- |
| `barge_hydro.csv` | `barge({ length: 20, beam: 8, depth: 4, draft: 2, divisions: [20, 8, 8] })` | 768 | 0.1 rad/s and periods 25, 20, 16, 14, 12, 10, 8, 7, 6, 5, 4 s; 0°, 90°, 180° |
| `cylinder_hydro.csv` | `cylinder({ radius: 4, height: 44, base: [0, 0, -40], segments: 24, heightDivisions: 22, capRings: 4 })` | 1128 | same frequencies; 0° |
| `hemisphere_hydro.csv` | `sphere({ radius: 5, subdivisions: 4 })` | 2592 | 0.5 and 1.0 rad/s; 0° |
| `barge5_hydro.csv` | the barge again | 768 | 0.25, 0.5, 0.75, 1.0, 1.25 rad/s; 0°, 90°, 180° (the NEMOH case) |
| `bargefine_hydro.csv` | `barge({ length: 20, beam: 8, depth: 4, draft: 2, divisions: [30, 12, 12] })` | 1728 | periods 25, 8, 6, 5, 4 s; 90°, 135°, 180° |
| `barge100_hydro.csv` | `barge({ length: 100, beam: 20, depth: 10, draft: 5, divisions: [40, 10, 8] })` | 1600 | periods 25, 20, 16, 14, 12, 10, 9, 8, 7, 6, 5 s; 90°, 135°, 180° |
| `barge80_hydro.csv` | `barge({ length: 80, beam: 24, depth: 6, draft: 3, divisions: [40, 12, 6] })` | 1584 | periods 20, 16, 14, 12, 10, 9, 8, 7, 6, 5, 4 s; 90°, 135°, 180° |

The last three files are the references of the box motion model (`boxBargeRAO`, `src/engine/hydro/radiation.js`): the
same barge on a finer mesh, a barge of ship proportions (L/B = 5, B/T = 4) and a wide shallow one (B/T = 8). They
were written by the same exported script, run with Capytaine 3.0.0 as the others.

Rows `AINF` and `A0` are radiation problems solved at omega = inf and omega = 0. The excitation is already
conjugated to exp(+i omega t) by the script.

## `nemoh/` — NEMOH results of the exported barge case (`writeNemohCase`)

Barge as above, 768 panels, 0.25 … 1.25 rad/s (5 values), headings 0°, 90°, 180°. Run: `preProc .`, `solver .`,
`postProc .` in the unpacked case folder. `Nemoh.cal` is the input as exported; `RadiationCoefficients.tec`,
`ExcitationForce.tec`, `CM.dat` (added mass) and `CA.dat` (damping) are from `results/`.

## `calculix/` — printed tables (`.dat`) and one ASCII results file (`.frd`) of exported decks (`writeCalculixCase`)

| File | Case |
| --- | --- |
| `cantilever_static.dat`, `.frd` | tube 500 × 20 mm, 10 m along x, 10 members, fixed at node 1, 10 kN downwards at the tip; B32R with PIPE section |
| `cantilever_frequency.dat` | same beam, `*FREQUENCY`, 6 modes |
| `cantilever_buckle.dat` | same beam, 1 MN axial compression at the tip, `*BUCKLE`, 4 modes |
| `rect_static_b31.dat`, `rect_static_b32.dat` | the cantilever with a rectangle b = 0.2 m, h = 0.4 m (depth vertical), B31 and B32 |
| `portal_static.dat` | portal frame in the x–z plane, columns 5 m, beam 6 m, 5 members each, same tube, fixed bases, 100 kN in x at the top of the first column |
| `portal_frequency.dat` | same frame, 6 modes |
| `portal_buckle.dat` | same frame, 1 MN downwards on each column top, 4 modes |
| `plate_static.displacements.dat` | plate 2 × 1.5 m, t = 20 mm, 12 × 10 cells of two triangles, edges held in z, 10 kPa; S6 shells. **Excerpt**: the displacement table of `plate_static.dat`; the reaction and stress tables that follow it in the original file (292 kB in total) were cut off |
| `plate_static_s3.displacements.dat` | same plate with S3 shells (`elementOrder: 1`). **Excerpt** as above |

## `gdf/capytaine_gdf_reader.json`

Not solver output: what Capytaine's mesh reader (`cpt.load_mesh(..., file_format="gdf")`) reported for the `.gdf`
files that `writeWamitGdf` wrote for the three hulls above — panel and vertex counts, wetted area, the volume
integrals of z n_z and x n_x over the panels with the reader's own normals, vertical extent and the area-weighted
sum of the normals. WAMIT itself was not available.
