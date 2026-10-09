# Building a suite

You are writing the analysis modules of one suite of Offshore Analyses. The numerical engines exist already
under `src/engine/` — your job is to turn them into complete, self-explanatory engineering analyses.

## Read before writing anything

1. `docs/CONTRACT.md` — all of it. §2 is the module contract, §2.2a (explanations and recommendations) is
   mandatory, §3 is the project bus.
2. `src/spec/suite-NN.json` for your suite — all of it: `analyses` (what must be covered), `models` (equations and
   models that must be applied), `inputs` and `outputs` (the industrial data the suite must accept and produce),
   `ic` / `bc` (the conditions each analysis needs), `cvv` (the verification the suite is held to).
3. `docs/engine/*.md` — one line per engine export. Open the engine source for exact argument names and units
   before calling anything; do not guess a signature.

## What to deliver

`src/suites/sNN/index.js` plus one file per module (shared helpers go in files whose name starts with `_`).
You own that folder only. Do not edit engines, tests, core or UI. If an engine function is wrong or missing,
work around it inside your suite and say so in your report (function, symptom, what you did).

### Coverage

* Every entry of `analyses` must be covered by at least one module (`covers`). Give each analysis its own module
  wherever the physics, the inputs or the acceptance criteria differ; merge only genuinely inseparable ones.
  Expect roughly as many modules as analyses, give or take a few.
* Across the suite, every item in the spec's `outputs` list must be produced by some module as a KPI, check,
  chart or table, and the `inputs` categories must be accepted by the modules that need them.
* Every model family named in `models` that an engine implements should be used by some module. High-fidelity
  methods the engines do not provide (3-D CFD, explicit nonlinear FE, full potential-flow BEM and the like) are
  stated in `limits` as outside the module — never implied.

### Inputs

* Industrial defaults that describe a realistic asset, so each module runs meaningfully out of the box.
* Tag every input with `role`; initial conditions are `ic`, boundary conditions (supports, restraints, applied
  environmental and interface loads, prescribed motions, hydraulic/electrical end conditions) are `bc`.
  Reject incompatible combinations (over- or under-constrained) with a clear `Error`.
* Bind site and interface values to the bus (`bind`) so pulled site data and upstream results flow in; publish
  what downstream suites need (`publish` + static `publishes`); report material and marine-operation quantities
  in `quantities` (`steel_t`, `concrete_m3`, `rock_t`, `cable_km`, `pipe_km`, `copper_t`, `aluminium_t`,
  `polymer_t`, `vessel_days: { kind: days }`).
* Use the materials library for materials (`type: 'select', source: 'materials:<category>'`, `ctx.material(id)`).
* Where a real model would come from a file, accept it: `type: 'geometry'` with the right `kinds` (frame for
  jackets/topsides/frames, surface for hulls and solids, polyline for routes/centrelines, grid for bathymetry,
  section for cross-sections, points for layouts), always with a parametric fallback built from the form fields
  when the value is `null`. Tables (`type: 'table'`) for soil layers, load cases, scatter diagrams, layer stacks,
  line segments, cash-flow lines; `series` for time histories and curves.
* Mark the variables a designer would actually change with `design: true` and give the uncertain ones an
  `uncertain` default distribution.

### Results

* KPIs: the few numbers an engineer reports, each with `explain`. Acceptance limits via `limit` + `sense`.
* Checks: every code or criterion check with demand, capacity, `uc`, `ref` (standard and clause) and `explain`.
* Charts: the diagrams an engineer expects for that analysis (profiles along depth/arc length/chainage, curves vs
  the governing variable, time histories, spectra, RAOs, interaction diagrams with the operating point,
  envelopes, scatter/heat maps, polar plots, mode shapes, 3-D `model3d` views with nodal values and deformed
  shape for frames/meshes) — each with `explain`. Keep series ≤ ~1500 points.
* Tables: the tabulated results a report needs — each with `explain`.
* `notes`: what governed and why, in plain words.
* `advice`: per §2.2a — at least one `engineering` and one `economic` item on every run, pass or fail, each with
  `text`, `action`, `basis`; at least one `sustainability`. Compute the numbers you quote (required thickness,
  mass saved, tonnes CO₂e using 1.9 t CO₂e/t steel, 0.3 t/m³ concrete, 3.206 t CO₂ per t marine gas oil as
  screening factors and say they are screening factors). Do not state costs as facts: express economic
  consequences through quantities (tonnes, metres, vessel days, downtime hours, energy) and, if you give money,
  through a clearly-labelled indicative unit rate exposed as an input of role `cost`.

### Numerics

* Any discretisation is an input named in `mesh` (and `meshTime` for time steps) with at least three levels at
  ratio 2 and a `qoi` that converges; the default level must be fine enough for a ≈ 2 % answer and fast enough
  for < 1 s. Heavy options (long time-domain runs, fine meshes) are opt-in through inputs.
* `compute` is deterministic: fixed seeds for anything random.
* No NaN: throw `Error` with a message a user can act on.

### Verification

* At least one `verify` case per module against an independent closed-form, hand-calculated or published value
  (state the source). Derive the expected number in a comment. Where no independent benchmark exists, use
  limiting-case physics (e.g. zero wave height → hydrostatic; rigid limit; single-phase limit) rather than
  `pending`. `pending` is the last resort and must be rare.

## Gate

`SUITE=NN node --test tests/suites.test.js` must pass completely: contract shape, full coverage of the specified
analyses, default run, explanations on every output, engineering + economic + sustainability advice,
reproducibility, structured-cloneability, verification cases, and the discretisation study for meshed modules.
Run it often. Also run a few modules by hand with non-default inputs (extremes of the allowed ranges) and make
sure they neither throw unhelpfully nor return nonsense.

## House rules

No attribution lines, generator credits or tool mentions anywhere. Comments explain engineering intent (the
reference, the assumption, the unit trap). Do not commit.

## Report back (short)

Module list (id → analyses covered → fidelity), gate result (tests passed/failed), the bus keys you publish and
bind, engine problems you worked around, and anything in the spec you could not cover and why.
