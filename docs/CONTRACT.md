# Offshore Analyses — code contract

Static, offline-first web application. No build step, no framework, no CDN: plain ES modules served as files.
Everything under `src/engine` and `src/suites` must also run under Node 20 (`node --test`), so it must not touch
`window`, `document`, `fetch`, or any browser API.

```
src/engine/<area>/*.js     pure numerical engines (shared by all suites)
src/suites/sNN/*.js        analysis modules for suite NN (declarative, see §2) + index.js
src/spec/*.json            the 17-suite specification (analyses, models, I/O, IC/BC, CVV, geometry formats)
src/geometry, src/mesh     importers, mesh tools
src/core, src/ui, src/data app shell, generic renderers and studies, live-data connectors
tests/**                   node:test suites
```

## 1. Engine conventions

* ES modules, named exports, no default export, no dependencies outside `src/engine`.
* SI units everywhere: m, kg, s, N, Pa, J, W, K (temperatures in °C only where the argument name ends in `C`),
  radians for angles (argument names ending in `Deg` are degrees). State units in JSDoc for every argument and
  every returned field.
* Functions take a single options object when there are more than three arguments and return plain objects with
  named fields (never positional arrays of mixed meaning).
* Deterministic: any randomness uses the seeded generator in `src/engine/math/rng.js`.
* Every public function has JSDoc stating: what it computes, the reference (standard + clause, or paper), and the
  validity range. Invalid input throws `RangeError` with a message naming the argument.
* Empirical correlations must carry their published coefficients exactly. If you are not certain of a coefficient,
  do not guess: derive the quantity from first principles, or expose the coefficient as an argument with a
  documented typical range.
* Every engine file has a test file `tests/engine/<area>.<file>.test.js` (`node:test` + `node:assert/strict`).
  Tests are benchmarks: closed-form solutions, hand calculations written out in the test, conservation or
  equilibrium checks, convergence-order checks, or published worked examples (cite them in the test name).
  Run with `node --test tests/`.
* Performance matters (phones): typed arrays for anything large, no O(n³) on n > ~1500 without a banded/sparse
  path, no allocations inside inner loops.

## 2. Analysis module contract

One file per analysis module in `src/suites/sNN/`, default-exporting a plain object:

```js
export default {
  id: 's08-on-bottom-stability',      // unique, 'sNN-' prefix, kebab-case
  suite: 8,
  title: 'On-bottom stability',
  summary: 'One sentence: what it answers.',
  covers: [1, 9],                      // 1-based indices into analyses[] of src/spec/suite-NN.json
  fidelity: 'analytical',              // 'analytical' | 'numerical' | 'reduced-order' | 'screening'
  standards: ['DNV-RP-F109 (2021) §4.4'],
  models: ['Morison equation', 'JONSWAP spectrum', 'Coulomb friction'],
  limits: 'Assumptions and validity range, 1–3 sentences. Say plainly what is not modelled.',
  publishes: ['pipeline.submergedWeight'],  // bus keys this module writes (static list, drives the wiring map; must match `publish` in compute)
  searchTerms: 'subsea pipeline on-bottom stability',  // optional: query for the live literature feed
  inputs: [ /* see 2.1 */ ],
  mesh: { key: 'nElem', label: 'Elements', levels: [10, 20, 40, 80], qoi: 'maxDisp' },   // optional, see 2.3
  compute(inp, ctx) { /* see 2.2 — may be async */ },
  verify: [ /* see 2.4 — at least one */ ],
};
```

`src/suites/sNN/index.js`:

```js
import a from './on-bottom-stability.js';
export default { n: 8, modules: [a /* , ... */] };
```

### 2.1 Inputs

```js
{ key: 'D', label: 'Steel outer diameter', unit: 'm', type: 'number', default: 0.3239, min: 0.05, max: 2,
  group: 'Pipe', role: 'geometry', help: 'Short hint shown under the field.',
  bind: 'pipeline.D',          // optional: default comes from the project bus when that key exists (see §3)
  design: true,                // optional: the design advisor may vary it to satisfy failed checks
  uncertain: { dist: 'normal', cov: 0.03 } }   // optional: default distribution for Monte Carlo studies
```

* `type`: `number` | `int` | `select` (`options: [{ value, label }]`) | `bool` | `text` |
  `table` (`columns: [{ key, label, unit, type }]`, `default: [ {..row..} ]`) |
  `series` (`default: { x: [], y: [] }`, with `xLabel`, `yLabel`, `xUnit`, `yUnit`) |
  `geometry` (`kinds: ['frame' | 'surface' | 'volume' | 'polyline' | 'grid' | 'points' | 'section']`; value is an
  imported model, see §4, or `null`).
* `unit` is the unit the user types in. It may be an engineering unit (`mm`, `MPa`, `kN`, `deg`, `°C`, `bar`,
  `t`, `MW`, `USD` …). **`compute` receives the values exactly as displayed** — convert to SI yourself at the top
  of `compute` and say so in a comment. Results are likewise reported in the engineering unit you declare.
* `role`: `geometry` | `material` | `load` | `environment` | `soil` | `ic` | `bc` | `operation` | `control` |
  `criteria` | `cost`. Every module must tag its initial conditions (`ic`) and boundary conditions (`bc`)
  explicitly — the app reports the exact IC/BC set used for each run.
* `group` is the form section heading. Order inputs as an engineer would fill them in.
* Give every input a realistic industrial default so the module runs out of the box, and `min`/`max` bounds that
  reject non-physical values.
* `showIf: { key: 'soil', in: ['clay'] }` hides an input unless another input has one of the listed values.
* Materials: use `type: 'select'` with `source: 'materials:steel'` (or `concrete`, `chain`, `wire`, `rope`,
  `polymer`, `soil`, `rock`, `conductor`, `insulation`, `fluid`) — the app fills the options from the materials
  library, and `ctx.material(id)` returns the property record.

### 2.2 `compute(inp, ctx)` result

```js
return {
  kpis: [   // headline numbers; key must be unique within the module
    { key: 'sfLat', label: 'Lateral stability safety factor', value: 1.32, unit: '-', limit: 1.0, sense: 'min',
      explain: 'Resisting over driving lateral force in the design storm; below 1.0 the pipe is predicted to move.' },
    //  sense 'min': ok when value >= limit; 'max': ok when value <= limit; omit limit for plain information
  ],
  checks: [ // code / acceptance checks; uc = demand / capacity
    { name: 'Lateral stability (absolute)', demand: 812, capacity: 1070, unit: 'N/m', uc: 0.76,
      ref: 'DNV-RP-F109 Eq. 3.38' },
  ],
  plots: [  // rendered by the shared plot kit
    { type: 'line', title: 'Stability vs wave height',
      x: { label: 'Hs', unit: 'm' }, y: { label: 'Safety factor', unit: '-' },
      series: [{ name: 'Lateral', x: [..], y: [..] }, { name: 'Limit', x: [..], y: [..], dash: true }] },
    // other types:
    //  'scatter' (same shape as line), 'bar' ({ categories: [], series: [{ name, y: [] }] }),
    //  'area', 'hist' ({ values: [], bins }), 'polar' ({ series: [{ name, thetaDeg: [], r: [] }] }),
    //  'heatmap' / 'contour' ({ xs: [], ys: [], z: [[row per y]], zLabel, zUnit }),
    //  'tornado' ({ base, items: [{ name, low, high }] }),
    //  'model3d' ({ nodes: [[x,y,z]], lines: [[i,j]], tris: [[i,j,k]], nodeValues: [], valueLabel, valueUnit,
    //              deformed: [[x,y,z]] })   // z is up
  ],
  tables: [ { title: 'Load cases', columns: [{ key: 'case', label: 'Case' }, { key: 'F', label: 'Force', unit: 'kN' }],
              rows: [{ case: '100-yr', F: 12.3 }] } ],
  notes: ['Plain-language remarks about what governed and why.'],
  advice: [ { severity: 'fail', text: 'Lateral stability is not met.', action: 'Increase concrete coating to ≥ 55 mm or trench.' } ],
  publish: { 'pipeline.submergedWeight': 1234.5 },   // SI values on the project bus (§3)
  quantities: { steel_t: 120, concrete_m3: 40, vessel_days: { pipelay: 12 } },   // feeds Economics + sustainability
};
```

* Arrays in plots/tables should stay small enough to render on a phone (≤ ~2000 points per series; decimate).
* `ctx` provides: `ctx.site` (project site record, §3, may be partly empty), `ctx.bus` (object of published
  values), `ctx.material(id)`, `ctx.progress(fraction, message)`, `ctx.signal` (AbortSignal). Never assume
  `ctx.site` fields exist — fall back to the module input.
* `compute` must be pure with respect to its inputs: same inputs → same outputs (the generic studies call it
  hundreds of times for sensitivity, Monte Carlo, calibration, optimisation and mesh convergence). Keep a single
  evaluation under ~1 s on a laptop at default settings.
* Throw `Error` with a user-readable message for inconsistent input (e.g. over- or under-constrained boundary
  conditions, non-converged solve). Do not return NaN silently.

### 2.2a Explanations and recommendations are mandatory

Every piece of generated output must explain itself to an engineer who has not read the standard that day:

* every **KPI** has `explain`: one or two sentences — what the quantity is, how to read it, and what a good or bad
  value means for the asset (e.g. "Ratio of resisting to driving lateral force on the pipe in the 100-year storm.
  Below 1.0 the pipe is predicted to move; DNV-RP-F109 asks for ≥ 1.0 with the safety-class factor included.").
* every **check** has `explain`: what failure mode it guards against and what happens physically if it is exceeded.
* every **plot** has `explain`: what is plotted and what to look for (the feature that matters: a peak, a crossing
  of the limit line, a knee).
* every **table** has `explain`: what each row represents and how to use it.
* `advice` is never empty and always speaks to BOTH engineering and economics. Each item is
  `{ kind: 'engineering' | 'economic', severity, text, action, basis, sustainability? }`. Every run returns at
  least one `engineering` item and at least one `economic` item. The economic item states the cost, schedule or
  value consequence of the result or of the recommended action in the engineer's terms (tonnes, vessel days,
  downtime, energy yield, inspection frequency → money), with `basis` naming the economic practice (e.g. AACE
  estimate classes, whole-life costing per ISO 15686-5, risk-based inspection per DNV-RP-G101 / API 580,
  LCOE conventions, value-of-information). In the Economics suite the `engineering` item points back to the
  technical driver behind the number (the quantity, yield, availability or design life that moves it most).
  Item fields:
  `text` states the finding with its number; `action` says what to do, quantified where the result allows it;
  `basis` names the current practice it rests on (standard + edition/clause, recommended practice, or established
  industry guidance — use the latest edition you know and say which); at least one item per run carries
  `sustainability`: the material, energy, emissions, seabed-footprint, durability or end-of-life consequence of the
  result or of the recommended action (e.g. "Each extra millimetre of wall adds ~X t of steel ≈ Y t CO₂e; a
  higher grade reaches the same utilisation with less mass."). Give advice for passing results too: where the
  reserve is, what could be optimised, what to confirm before relying on the result.
  For the Economics suite `basis` cites economic/financial practice (e.g. AACE RP 18R-97 estimate classes,
  IEA/IRENA LCOE conventions, project-finance covenants) and `sustainability` covers carbon cost, taxonomy
  alignment, circularity or social/environmental externalities as relevant.
* Wording: plain, specific, no filler. Write for a practising engineer on a phone.

### 2.3 Discretisation and convergence (`mesh`)

Any module that discretises (elements, segments, panels, grid cells, time steps, frequency bins, samples) exposes
the control as an input and names it in `mesh`. The app then runs the refinement ladder, computes observed order,
Richardson-extrapolated value and Grid Convergence Index, and selects the coarsest level meeting the tolerance.
`qoi` is the `key` of the KPI to converge on. A second ladder may be given as `meshTime: { key, levels, qoi }` for
time-step convergence. The control is assumed to be a count (elements, panels, steps): representative size h ∝ 1/value.
If the control is itself a size (element length, time step) add `isSize: true`. Levels must be at least three,
with a constant refinement ratio (normally 2).

### 2.4 Verification cases (`verify`)

```js
verify: [{
  name: 'Simply supported span, f1 = (π/2)·sqrt(EI/(m L⁴))',
  source: 'Blevins, Formulas for Natural Frequency and Mode Shape, Table 8-1',
  inputs: { L: 40, /* overrides of defaults */ },
  expect: [{ kpi: 'f1', value: 0.6183, tol: 0.01 }],     // tol is relative; use abs: for absolute
}]
```

These run in the app's Verification page and in `node --test`. Use independent closed-form or published values —
not numbers copied from the module's own output. A module with no independent benchmark available must say so:
`verify: [{ name: 'No independent benchmark available', pending: true }]`.

The app never claims more than the evidence: modules whose `verify` cases pass are shown as *Numerically verified*;
everything else is *Implemented*. Nothing is labelled experimentally or field validated — users validate against
their own measurements with the built-in calibration/validation tools.

## 3. Project bus (cross-suite wiring)

Modules exchange values through named keys. `publish` writes; `bind` on an input reads (as the default).
All bus values are SI. Canonical keys (add new ones under your suite's prefix, `dot.separated`):

| key | meaning (SI) |
| --- | --- |
| `site.lat`, `site.lon` | decimal degrees |
| `site.depth` | still-water depth, m (positive) |
| `site.Hs1`, `site.Tp1`, `site.Hs100`, `site.Tp100`, `site.Hmax100` | 1-yr / 100-yr sea state, m, s |
| `site.Uc1`, `site.Uc100`, `site.UcBed100` | surface / near-bed current, m/s |
| `site.Uw1`, `site.Uw100` | 1-hour mean wind at 10 m, m/s |
| `site.tideRange`, `site.surge100` | m |
| `site.rhoW`, `site.Tsea`, `site.Tbed`, `site.Tair` | kg/m³, °C |
| `site.pga475` | peak ground acceleration, m/s² |
| `site.su0`, `site.suGrad`, `site.phiDeg`, `site.gammaSub`, `site.d50` | soil: Pa, Pa/m, deg, N/m³, m |
| `site.windWeibullA`, `site.windWeibullK`, `site.hubWindMean` | m/s, –, m/s |
| `site.ghiAnnual` | kWh/m²/yr |
| `topsides.weight`, `topsides.cogZ`, `topsides.Fx`, `topsides.My` | N, m, N, N·m at jacket interface |
| `jacket.baseShear`, `jacket.otm`, `jacket.pileAxialMax`, `jacket.pileShearMax`, `jacket.f1` | N, N·m, N, N, Hz |
| `foundation.Kh`, `foundation.Kr`, `foundation.Khr`, `foundation.Kv` | pile-head stiffness: N/m, N·m/rad, N, N/m |
| `owt.thrustRated`, `owt.rna`, `owt.f1`, `owt.mudlineMoment`, `owt.mudlineShear` | N, kg, Hz, N·m, N |
| `vessel.disp`, `vessel.GM`, `vessel.offsetMax`, `vessel.heaveSig`, `vessel.rollSigDeg`, `vessel.Tn3`, `vessel.Tn4` | kg, m, m, m, deg, s, s |
| `vessel.Fmean` | mean environmental force on the floater, N |
| `mooring.Tmax`, `mooring.Kxx`, `mooring.anchorH`, `mooring.anchorV`, `mooring.pretension` | N, N/m, N, N, N |
| `riser.topTension`, `riser.hangoffAngleDeg`, `riser.tdpTension`, `riser.reactionH` | N, deg, N, N |
| `pipeline.D`, `pipeline.t`, `pipeline.length`, `pipeline.submergedWeight`, `pipeline.endForce`, `pipeline.endExpansion`, `pipeline.dP`, `pipeline.Tout` | m, m, m, N/m, N, m, Pa, °C |
| `cable.length`, `cable.ampacity`, `cable.lossesFrac`, `cable.topTension` | m, A, –, N |
| `umbilical.topTension`, `flexible.topTension` | N |
| `install.liftWeight`, `install.hookLoad`, `install.operabilityFrac`, `install.vesselDays` | N, N, –, days |
| `energy.ratedPower`, `energy.aep`, `energy.capacityFactor`, `energy.wakeLossFrac`, `energy.availability` | W, Wh/yr, –, –, – |
| `quantities.steel_t`, `quantities.concrete_m3`, `quantities.rock_t`, `quantities.cable_km`, `quantities.pipe_km` | t, m³, t, km, km |
| `life.fatigueYears`, `life.designYears` | years |
| `econ.capex`, `econ.opexAnnual`, `econ.npv`, `econ.lcoe` | project currency, /yr, –, per MWh |

## 4. Imported geometry model

Importers return one or more *models* of these kinds (arrays are typed arrays; lengths in metres after unit scaling):

```js
{ kind: 'surface', name, positions: Float32Array(3n), indices: Uint32Array(3m), groups?: [{ name, start, count }] }
{ kind: 'volume',  name, positions, cells: [{ type: 'tet4'|'tet10'|'hex8'|'hex20'|'wedge6'|'pyr5', conn: Uint32Array }], sets?: {..} }
{ kind: 'frame',   name, nodes: Float64Array(3n), members: [{ a, b, section?: { type: 'tube', D, t } | { A, Iy, Iz, J }, material?, group? }] }
{ kind: 'polyline',name, points: Float64Array(3n), attrs?: { kp?: Float64Array, ... } }        // routes, centrelines
{ kind: 'grid',    name, nx, ny, x0, y0, dx, dy, z: Float32Array(nx*ny), crs?: string }      // bathymetry / DEM (z up, m)
{ kind: 'points',  name, positions: Float32Array(3n), attrs?: {..} }                             // point clouds, soundings
{ kind: 'section', name, loops: [Float64Array(2n)], layers?: [...] }                              // 2-D cross-sections
```

## 5. House rules

* No attribution lines, generator credits or tool mentions in code, comments, docs or commit messages.
* Comments explain engineering intent (reference, assumption, unit trap), not what the next line does.
* Do not commit. Do not edit files outside the paths you were assigned.
