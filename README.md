# Offshore Analyses

Engineering analysis suite for fixed offshore platforms and jackets, topsides, offshore wind turbine foundations,
subsea support structures, wharves, dolphins and other marine structures, subsea pipelines, risers, mooring
systems, floating structures, umbilicals, cables, flexible pipes, installation operations, offshore
renewable-energy systems — and the economics of all of them.

**Open the app:** https://samuelakosaonyejekwe.github.io/offshoreanalyses/

Mirror (served from this repository through a CDN, independent of GitHub Pages):
https://raw.githack.com/samuelakosaonyejekwe/offshoreanalyses/main/index.html

## What it does

* **17 connected suites, 292 analysis modules.** Every analysis in the 17-suite specification (261 of 261) is
  covered by a working module. What one analysis computes (loads, motions, reactions, stiffnesses, quantities,
  energy) becomes the default input of the analyses that depend on it.
* **Any site on Earth.** Waves, wind, currents, tide, water depth, seismicity and solar resource are pulled live
  by the browser from public sources and reduced to design values; measured data overrides them.
* **Your geometry and data.** 68 CAD, mesh, structural, centreline, bathymetry, point-cloud and tabular formats
  are recognised; geometry is checked, repaired, refined or coarsened, and viewed in 3-D. The compatibility table
  in the app states exactly what is read natively, what loads an adapter on first use, and the few proprietary
  formats that must be exported to a neutral format first.
* **Every result explains itself.** Each number, check, chart and table carries a plain-language explanation;
  each run ends with engineering and economic recommendations, their basis, and the sustainability consequence.
* **Studies on every module.** Sensitivity ranking, parameter sweeps, uncertainty propagation, discretisation
  convergence with automatic mesh selection (observed order, Richardson extrapolation, Grid Convergence Index),
  calibration against measurements and validation against independent data.
* **Managed assurance.** A code-factor register lists every standards-derived factor (1,242 of them) with a
  sign-off against the governing edition; an evidence register records calibration and validation against
  measurements with the range covered; the modelled sea level can be checked against tide-gauge measurements.
* **Hand-off to high-fidelity solvers.** Ready-to-run cases for OpenFOAM, NEMOH, Capytaine, WAMIT and CalculiX are
  written from the app's own models and their results read back; a built-in boundary-element solver gives the
  added-mass matrix of any imported hull. `docs/interop.md` states which exports were run with the solver.
* **Works offline.** Installable on desktop, Android and iOS; all calculations run on the device; projects stay
  on the device.

## Run it locally

No build step. Any static file server works:

```
npm run serve        # http://localhost:8080
npm test             # ~3,500 tests: engine benchmarks, importers, mesh tools and the suite contract gate
npm run build        # pack each suite into one file (dist/), refresh the catalogue and the offline file list
```

## Host it anywhere

The repository root is the website. Publish the folder as-is on any static host (GitHub Pages, Cloudflare Pages,
Netlify, GitLab Pages, an S3 bucket, an intranet server). Add the address to `mirrors.json` so installed copies
can fall back to it when another host is unreachable.

## Layout

```
src/engine    numerical engines (pure ES modules, tested under Node)
src/suites    analysis modules, one folder per suite, plus the generated catalogue
src/spec      the 17-suite specification shown inside the app
src/geometry  file importers, exporters and adapters      src/mesh   mesh tools and convergence toolkit
src/interop   case writers and result readers for external solvers
src/core      project store, analysis worker, generic studies
dist          packed suites and worker (built by `npm run build`; a test fails if they are stale)
src/data      live-data connectors and site characterisation
src/ui        application shell and views
vendor        third-party readers loaded on first use (versions and licences in docs/vendor.md)
docs          module contract and engine reference
```

## Engineering responsibility

The software performs calculations; it does not replace a qualified engineer or an independent design
verification. Each module states its method, reference standard, limits and verification status inside the app.
Modules whose benchmark cases pass are numerically verified — that shows the equations are solved correctly, not
that a model has been validated against a physical asset; the built-in calibration and validation tools are
there for that, with your own measurements. Open global datasets are for screening and early design; site
investigation and measured data take precedence.
