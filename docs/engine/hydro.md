# `src/engine/hydro` — hydrodynamics and floating-body engines

SI units throughout; angles in radians unless the name ends in `Deg`. z is up; for wave problems z = 0 is the
still-water level and the seabed is at z = −depth. Wave heading = direction the waves travel towards, measured from
+x towards +y (0° following, 180° head sea). Complex amplitudes use exp(+iωt) per unit wave amplitude; a positive
phase is a lead on the crest at the origin. Floating position for mesh hydrostatics:
`p_earth = Ry(trim)·Rx(heel)·p_body − (0, 0, draft)` (heel positive = port side up, trim positive = bow down,
draft = depth of the body origin below the waterline). Tests: `tests/engine/hydro.<file>.test.js`.

## morison.js
- `morisonForcePerLength({ rho, D, Cd, Cm, Ca?, axis, u, a, vBody?, aBody? }) -> { f, fDrag, fInertia }` — N/m vectors; normal-component Morison with relative velocity. DNV-RP-C205 §6.2.
- `memberWaveForce({ p1, p2, D, Cd, Cm, rho, g, wave: { H, T, depth, headingDeg, phaseDeg }, current?, t, stretching: 'none'|'wheeler'|'vertical', nGauss, nSeg, about }) -> { F, M }` — N, N·m; Gauss integration over the wetted length in Airy kinematics. DNV-RP-C205 §6.2; Wheeler (1970).
- `memberForceHistory({ …memberWaveForce, nSteps }) -> { t, Fx, Fy, Fz, Mx, My, Mz, Fh, Mh, FhMax, tFhMax, MhMax, tMhMax, FzMax }` — s, N, N·m over one wave cycle with maxima.
- `pileForceClosedForm({ H, T, depth, D, Cd, Cm, rho, g }) -> { FD, FI, Fmax, MD, MI, Mmax, k, phaseFmaxDeg }` — N, N·m about the mudline; closed-form vertical pile. Morison et al. (1950).
- `dragCoefficientSteady(roughnessRatio) -> C_DS` — 0.65 / (29 + 4 log10 Δ)/20 / 1.05. DNV-RP-C205 §6.7.1.5.
- `wakeAmplification(KC, CDS) -> psi` — KC < 12 only (throws above). DNV-RP-C205 §6.7.2.2.
- `inertiaCoefficient(KC, CDS) -> { Ca, Cm }` — DNV-RP-C205 §6.9.1.
- `keuleganCarpenter(Um, T, D) -> KC` — m/s, s, m.
- `marineGrowth({ D, thickness, rhoGrowth = 1325, rhoWater, g }) -> { De, massPerLength, submergedWeightPerLength }` — m, kg/m, N/m; geometry, density is an input.
- `maccamyFuchs({ H, T, depth, D, rho, g }) -> { F, M, CmEff, phaseLagDeg, ka, A, k, forcePerLength(z) }` — N, N·m (mudline), –, deg; linear diffraction on a vertical cylinder. MacCamy & Fuchs (1954).
- `rectangleAddedMassCoefficient(aOverB) -> Ca` — DNV-RP-C205 App. D Table D-1 (log-linear interpolation between tabulated ratios).
- `rectangleRollInertiaCoefficient(ratio) -> beta` — rectangle rotating about its centre, I_a = β ρ π a⁴ per unit length on the longer half-side; ratio of the sides 1, 2, 5, 10, ∞ → 0.234, 0.15, 0.15, 0.147, 0.125. DNV-RP-C205 App. D Table D-1. Half of it is the roll added inertia of a floating rectangular section about its waterline centre (high-frequency limit).
- `finiteCylinderAddedMassCoefficient(lengthOverDiameter) -> Ca` — circular cylinder of finite length moving normal to its axis, on its own volume; l/D 1.2, 2.5, 5, 9, ∞ → 0.62, 0.78, 0.90, 0.96, 1.0, linear between, 1 − 0.36/(l/D) beyond 9. DNV-RP-C205 App. D Table D-2.
- `addedMass2D({ shape: 'circle'|'ellipse'|'plate'|'square'|'rectangle', a, b?, rho }) -> { Ca, referenceArea, addedMassPerLength }` — –, m², kg/m. DNV-RP-C205 App. D Table D-1.
- `addedMass3D({ shape: 'sphere'|'disc', a, rho }) -> { addedMass, Ca, referenceVolume }` — kg; sphere ½ρ(4/3)πa³, disc (8/3)ρa³. DNV-RP-C205 App. D Table D-2.
- `memberBuoyancy({ p1, p2, D, t?, flooded, waterline, rho, g }) -> { volume, force, centre, submergedLength }` — m³, N, m. Archimedes.
- `slamForcePerLength({ D, v, Cs, rho }) -> N/m` — ½ρCsDv²; Cs is an input (π von Kármán … 5.15 Campbell & Weynberg). DNV-RP-C205 §8.
- `vortexSheddingFrequency({ U, D, St = 0.2 }) -> Hz` — St·U/D. DNV-RP-C205 §9.

## hydrostatics.js
- `bodyToEarth(p, draft, heel, trim) -> [x, y, z]`, `earthToBody(p, draft, heel, trim) -> [x, y, z]` — m; frame transformation defined above.
- `clipIntegrals(mesh, draft, heel, trim) -> { V, Mx, My, Mz, Aw, Sx, Sy, Jxx, Jyy, Jxy, wet }` — m³, m⁴, m², m³, m⁴, m²; exact integrals of the faceted closed surface below the flotation plane (divergence theorem).
- `bodyIntegrals(mesh, lost, draft, heel, trim) -> same` — hull minus permeability-weighted lost-buoyancy compartments.
- `meshHydrostatics({ mesh, draft, heel, trim, rho, g, cog?, mass?, length?, lost? }) -> { volume, displacement, buoyancy, cob, cobBody, wettedArea, waterplaneArea, cof, cofBody, Ixx, Iyy, Ixy, I1, I2, principalAngle, KB, BMt, BMl, KMt, KMl, GMt, GMl, TPC, massPerDraft, MCT, MCTcm, stiffness }` — SI; TPC t/cm, MCT N·m per m of trim, MCTcm t·m/cm, stiffness 6×6 about the waterline origin. Newman (1977) ch. 6; Biran & López-Pulido (2014).
- `boxMesh({ xMin, xMax, yMin, yMax, zMin, zMax })`, `cylinderMesh({ D, zMin, zMax, x0, y0, nTheta })`, `mergeMeshes(meshes) -> surface model` — closed outward-oriented meshes (contract §4), m.
- `boxHydrostatics({ L, B, T, KG?, rho, g }) -> { volume, displacement, waterplaneArea, KB, BMt, BMl, KMt, KMl, GMt, GMl, TPC, MCT, C33, C44, C55 }` — closed form (KB = T/2, BM = B²/12T).
- `verticalCylinderHydrostatics({ D, T, KG?, rho, g }) -> { volume, displacement, waterplaneArea, KB, BM, KM, GM, C33, C44 }` — closed form (BM = D²/16T).
- `multiColumnHydrostatics({ columns: [{ x, y, D | a, b, zBottom }], submerged: [{ volume, x, y, z }], T, KG?, rho, g }) -> { volume, displacement, waterplaneArea, cof, cob, KB, Ixx, Iyy, BMt, BMl, KMt, KMl, GMt, GMl, C33, C44, C55 }` — parallel-axis closed form.
- `draftForVolume({ mesh, volume, heel, trim, lost? }) -> draft` — m; bisection at fixed attitude.
- `solveEquilibrium({ mesh, mass, cog, rho, guess?, fixHeel, fixTrim, lost?, tol, maxIter }) -> { draft, heel, trim, iterations, converged, residual }` — m, rad; Newton on volume and B–G alignment.
- `freeSurfaceCorrection({ tanks: [{ l, b, rho } | { mesh, fillLevel, rho }], volume, rho }) -> { GGt, GGl, fsmT, fsmL, perTank }` — m, kg·m; GG′ = Σρ_t·i/(ρ∇). IMO 2008 IS Code Part B ch. 3.

## stability.js
- `gzCurve({ mesh, mass, cog, rho, heelsDeg, side, freeTrim, trim, lost?, openings? }) -> { heelDeg, gz, draft, trimDeg, openingHeight, downfloodDeg }` — deg, m; re-clipped righting levers at constant displacement.
- `areaUnderGz(heelDeg, gz, aDeg, bDeg) -> m·rad`.
- `gzProperties({ heelDeg, gz, downfloodDeg? }) -> { gzMax, heelAtGzMaxDeg, vanishingDeg, equilibriumDeg, rangeDeg, area0to30, area0to40, area30to40, gzAt30, gm0, limitDeg }` — m, deg, m·rad.
- `imoIntactCriteria({ heelDeg, gz, gm0?, downfloodDeg? }) -> { pass, checks: [{ name, value, limit, unit, ok, ref }] }` — IMO 2008 IS Code Part A §2.2 (0.055 / 0.09 / 0.03 m·rad, 0.20 m, 25°, 0.15 m).
- `windHeelingLever({ V, areas: [{ A, z, Cs, Ch }], zLateralResistance, displacement, rhoAir = 1.222, g }) -> { force, moment, lever }` — N, N·m, m. IMO MODU Code 2009 §3.2.3 (Cs, Ch are inputs).
- `moduWindCriterion({ heelDeg, gz, lever0, variation: 'constant'|'cos'|'cos2', unitType, downfloodDeg? }) -> { firstInterceptDeg, secondInterceptDeg, limitDeg, areaRighting, areaHeeling, ratio, required, positiveRange, pass }` — area ratio 1.4 (surface, self-elevating) or 1.3 (column-stabilised). IMO MODU Code 2009 §3.3.
- `damageStability({ mesh, mass, cog, rho, compartments: [{ mesh | box, permeability }], heelsDeg?, marginLine?, openings? }) -> { equilibrium, side, curve, properties, marginLine, lostVolume }` — lost-buoyancy method, final stage.
- `progressiveFlooding({ mesh, mass, cog, compartment, opening: { point, area, Cd }, dt, tMax, rho, g }) -> { t, floodMass, draft, heelDeg, trimDeg, level, head, flow, equalised, sunk, timeToEqualise }` — s, kg, m, deg, m³/s; Torricelli orifice + quasi-static re-equilibration (explicit Euler, one compartment).
- `ballastPlan({ tanks: [{ name, x, y, capacity, volume, rho }], hydro: { displacement, waterplaneArea, lcf, GMt, GMl, rho }, current, target, mode: 'transfer'|'fill-discharge', pumpRate?, nPumps }) -> { tanks, achieved, error, volumeMoved, pumpTime, feasible }` — m³, kg, m, deg, s; bounded least squares on small-angle hydrostatics.
- `incliningTest({ displacement, shifts: [{ mass, distance, deflection, pendulumLength }], KM? }) -> { GM, KG, perShift, spread }` — m. IMO 2008 IS Code Part B Annex 1.

## floater.js
- `buildFloater({ members, mass?, cog, gyradii, rho, g, depth, Kext?, nStrips }) -> model { strips, waterplanes, M, A, C, Kext, mass, volume, cob, waterplaneArea, maxDimension, extent, validOmegaMax, … }` — member types `cylinder`, `rect`, `plate`; 6×6 matrices about the still-water origin. Cylinder ends: `end1`/`end2` = `'none'`, `'half-disc'` (4/3 ρa³), `'hemisphere'` (2/3 πρa³), `'disc'` or kg for the axial added mass; `free1`/`free2` = true takes the transverse added mass lost around a free end, (1 − C_A(2l/D)) ρ Ca πD²/4 · l (DNV-RP-C205 Table D-2), off at that end. Faltinsen (1990) ch. 3; DNV-RP-C205 App. D. Valid for D/λ < 0.2, no diffraction.
- `naturalPeriods(model) -> { uncoupled: [6], coupled }` — s.
- `waveSpectrum(omega, { Hs, Tp, gamma }) -> Float64Array` — JONSWAP, m²·s/rad. DNV-RP-C205 §3.5.5.
- `radiationDamping({ model, omega, nHeadings? }) -> number[][][]` — 6×6 per frequency, N·s/m …; Haskind relation B_ij = k/(8πρg c_g) ∫ Re{X_i X_j*} dβ with the model's own excitation in waves from all directions (deep water: ω³/(4πρg³) ∫ …). Newman (1962).
- `floaterRAO({ model, omega, headingDeg, seaState?, radiation = 'haskind', extraDampingRatio = 0, extraDamping?, viscous, maxIter, tol }) -> { omega, period, headingDeg, depth, g, re, im, amp, phase, force: { re, im }, radiationDamping, viscousDamping, iterations, converged, validOmegaMax }` — m/m, rad/m, rad, N/m; long-wave Froude–Krylov + inertia excitation with the wave kinematics integrated exactly along each strip, radiation damping from the Haskind relation (`radiation: 'none'` to leave it out), `extraDampingRatio` (former name `radiationDampingRatio`) as an added fraction of critical for damping that is not modelled, stochastically linearised drag (Borgman 1967).
- `pointMotion({ rao, point }) -> { omega, re, im, amp, velAmp, accAmp }` — [x, y, z] arrays, m/m, (m/s)/m, (m/s²)/m.
- `responseStatistics({ omega, amp, S, duration }) -> { m0, m2, sigma, significant, Tz, nCycles, mpm, spectrum }` — response units; MPM = √(2 m0 ln N).
- `relativeMotion({ rao, point, waveFactor }) -> { omega, re, im, amp }` — m/m, wave minus vertical point motion.
- `airGap({ rao, point, staticAirGap, S, duration, waveFactor }) -> { sigma, mpmRelative, minAirGap, exceedancePerCycle, Tz }` — m. DNV-RP-C205 §8.2 (undisturbed wave).
- `operability({ seaStates: [{ Hs, Tp, prob }], evaluate(Hs, Tp), limits }) -> { operability, downtime, perState, governingCounts }` — fractions.
- `limitingHs({ omega, amp, limit, Tp, gamma, statistic, duration }) -> { Tp, HsLimit }` — s, m; linear responses only.
- `boxBargeRAO({ L, B, T, omega, headingDeg, speed, Cb, GMt?, kxx?, kxxDry?, kyy?, rollDampingRatio, rollRadiation, method?, rho, g }) -> { omega, heave, pitch, roll, rollPhase, heavePhase, pitchPhase, heaveNaturalPeriod, rollNaturalPeriod, rollAddedInertia, rollRadiationDamping, validPeriodMin: { heave, roll }, method }` — m/m, rad/m, rad.
  - `method: 'panel'` (default at zero speed): added mass, wave damping and exciting forces of the box from the panel solution of `radiation.js` (`boxHydrodynamics`: 14 frequencies up to the wavelength equal to the beam, interpolated, cached per geometry), equations of motion about the waterline centre in three groups: heave; surge–pitch (mass gyradius `kyy`, default 0.25 L); sway–roll (C44 = Δ g GMt, `rollDampingRatio` of critical added to the wave damping). Centre of gravity from GMt for the box, zG = B²/(12T) − T/2 − GMt (waterline when GMt is absent). `kxx` is the roll gyradius including added inertia, so the roll natural period is 2π kxx/√(g GMt); `kxxDry` is the mass gyradius and the computed added inertia is added. Cb < 1: heave and pitch from the box of breadth B·Cb, roll from the box of breadth B with the mass ρLBT·Cb. `validPeriodMin` = period of the wave whose length is the beam.
  - `method: 'closed-form'` (default with forward speed, the panel solution being for zero speed): heave and pitch by Jensen, Mansour & Olsen (2004) with the bottom-pressure factors exp(−kT) and sin(½k_yB)/(½k_yB) at every heading (the paper's exp(−|k cos β|T) is exact in head seas only); roll by a 1-DOF model (not from the paper): Froude–Krylov moment, relative-motion diffraction −ω²A44·(wave slope at the keel), strip added inertia A44 from DNV-RP-C205 Table D-1, Haskind radiation damping. `validPeriodMin` = twice 2π√(2T/g) (heave) and 1.3 × the roll natural period.

Validity of the floater models, from the comparison with Capytaine 3.0.0 (`src/spec/cross-solver.json`, `docs/interop.md` §1.4). Three boxes were compared, each with the centre of gravity on the waterline, dry gyradii 0.35 B and 0.25 L, 5 % of critical added roll damping in both solutions: 20 × 8 × 2 m (L/B 2.5, B/T 4; Capytaine meshes of 768 and 1728 panels), 100 × 20 × 5 m (L/B 5, B/T 4; 1600 panels) and 80 × 24 × 3 m (L/B 3.3, B/T 8; 1584 panels).

| Model | Within 5 % of the panel solution | Outside |
| --- | --- | --- |
| `boxBargeRAO` `'panel'`, heave and pitch, head, quartering and beam seas | head and beam seas: every compared period from 0.72 × the heave natural period upward (20 × 8 × 2: 4–25 s, natural period 4.4 s; 100 × 20 × 5: 5–25 s, 7.0 s; 80 × 24 × 3: 5–20 s, 6.6 s). The largest differences, −3 … −5 %, lie at 0.75–0.9 × the heave natural period, where the response falls steeply and follows the added mass (heave added mass of the model 1–2 % above Capytaine's); within 1.5 % from 1.15 × upward | at 135° −5.3 % at one period (100 × 20 × 5, 6 s); 80 × 24 × 3 at 4 s (wavelength = 1.04 B, 0.6 × the heave natural period): heave in beam seas −7.7 %, pitch at 135° −7.7 % |
| `boxBargeRAO` `'panel'`, roll, beam and quartering seas | from 0.75–0.85 × the roll natural period upward (20 × 8 × 2: 5–25 s, +6.3 % at 4 s against the 768-panel mesh and +3.1 % against the 1728-panel mesh; 100 × 20 × 5: 7–25 s; 80 × 24 × 3: 4–20 s). The reference is uncertain by 1–3 % in roll (its sway–roll added-mass terms are not symmetric); with their mean the model is within 3.1 % from 0.84 × the roll natural period | 100 × 20 × 5 below 0.75 × the roll natural period: +6 % at 6 s, +11 % at 5 s (beam), −7 % at 5 s (135°), where the roll is below 0.02 rad/m. The resonant peak follows the gyradius and the added damping that are entered |
| `boxBargeRAO` `'closed-form'` heave, head and beam seas | wave periods ≥ 2 × the heave natural period 2π√(2T/g) (`validPeriodMin.heave`) | low: −15 % at 1.5 ×, −30 … −40 % near it (added mass taken equal to the displacement, two-dimensional sectional damping) |
| `boxBargeRAO` `'closed-form'` pitch, head seas | periods ≥ 1.5 × the heave natural period | −19 % at 1.25 ×, −63 % at 1.0 × (wavelength ≈ L) |
| `boxBargeRAO` `'closed-form'` roll, beam seas | periods ≥ 1.3 × the roll natural period (`validPeriodMin.roll`) | within 15 % below; no sway–roll coupling, strip added inertia |
| `floaterRAO` forces, single column | D/λ < 0.2 (`validOmegaMax`): surge and heave force from 5–6 s for D = 8 m | surge force +6 % at D/λ = 0.2, +47 % at 0.32 (diffraction) |
| `floaterRAO` motions, single column | heave at all periods above the limit, resonance included, with `end1: 'hemisphere'`; surge and pitch outside 0.9–1.4 × the pitch natural period | within 25 % near the pitch resonance: the peak has radiation damping only and follows the natural period (1.2 % high with `free1`, 3 % without) |

Not covered by `boxBargeRAO` `'panel'`: hull forms other than a box (bilge radius, rake, ship lines), finite water depth, forward speed, wavelengths shorter than the beam (the coefficients of the last solved frequency are held), viscous roll damping other than the fraction entered. The column compared for `floaterRAO`: D = 8 m with 40 m draught.

## radiation.js
Three-dimensional radiation problem of a surface-piercing hull in deep water at zero speed: constant-potential quadrilateral panels, Green's identity with G = 1/r + 1/r′ + 2K[F(X, Y) − iπ e^Y J0(X)] (K = ω²/g, X = KR, Y = K(z + ζ)). The Rankine parts are integrated in closed form (routines of `bem.js`), the wave part is tabulated: F(X, Y) = e^Y[F(X, 0) − ∫₀^|Y| e^t/√(X² + t²) dt] with F(X, 0) = −(π/2)[H0(X) + Y0(X)] (Newman 1985), the logarithmic singularity split off, cubic interpolation on a grid uniform in √X and √|Y|. Planes of symmetry x = 0 and y = 0 reduce the unknowns to a quarter of the hull. Exciting forces from the Haskind relation (no diffraction solution). Irregular frequencies are removed with null-field equations at interior points (Schenck 1968).
- `waveGreenFunction(X, Y, out?) -> { F, FX, FY, J0, J1 }` — wave part and its derivatives; 0 ≤ X ≤ 100, Y ≤ 0; within 3·10⁻⁵ of direct quadrature of the defining principal-value integral.
- `hullPanels(quads) -> panel set` — 12 numbers per panel (four vertices, counter-clockwise seen from the fluid; a triangle repeats its last vertex), m.
- `boxQuarterPanels({ L, B, T, nx, ny, nz, grading = 1 }) -> Float64Array` — quarter of a box, cosine spacing towards the edges.
- `radiationSetup({ panels, symmetry: { x?, y? }, ref = [0,0,0], interiorPoints? }) -> setup` — frequency-independent influence matrices.
- `radiationSolve({ setup, omega, rho, g, modes?, nullField = true, split = true }) -> { omega, A, B, phi, symmetryError }` — 6×6 added mass and damping about `ref` (kg … N·m·s), potentials per mode on the panels.
- `haskindExcitation({ setup, phi, omega, headingDeg, rho, g }) -> { re, im }` — N/m, N·m/m; X_i = iωρ∫(φ0 n_i − φ_i ∂φ0/∂n) dS, 2 × 2 Gauss rule per panel. Newman (1962).
- `haskindDamping({ setup, phi, omega, rho, g, nHeadings? }) -> number[6][6]` — B_ij = ω³/(4πρg³)∫Re{X_i X_j*}dβ, symmetric and positive semi-definite, the heading integral in closed form (Fourier–Bessel series) when the hull has both planes of symmetry; agrees with the pressure integral of `radiationSolve` within 1–2 % except for roll about the waterline, where the pressure integral converges slowly.
- `boxHydrodynamics({ L, B, T, rho, g, mesh?, nOmega = 14, lambdaMin?, modes? }) -> model` — box solved on a frequency grid (quarter mesh of 5 panels over the half beam, up to 14 over the half length, 2–6 over the draught: 80–130 unknowns), in non-dimensional form (unit beam, density and gravity) and kept by L/B and T/B, so a change of scale, density or gravity is not solved again; damping from `haskindDamping`, null-field equations from 0.9 × the first irregular frequency (K = κ coth κT, κ = π√(1/L² + 1/B²)); `model.unit` holds the non-dimensional solution.
- `boxCoefficients(model, omega) -> { A, B }`, `boxExcitation(model, omega, headingDeg) -> { re, im }` — cubic interpolation between the solved frequencies (potentials interpolated, the incident wave exact); end values held outside the range.

Checks (`tests/engine/hydro.radiation.test.js`): floating hemisphere, zero-frequency added mass A33/ρV = 0.826 and A11/ρV = 0.498 on 81 panels per quarter against 0.830951 and 0.5 (Hulme 1982, J. Fluid Mech. 121, as quoted by Storti & D'Elía 2001, arXiv:math/0110302, eq. 38), the faceted volume being 1.3 % below that of the hemisphere; the same hemisphere against Capytaine at Ka = 0.13 and 0.51: added mass −1.4 … −2.0 %, damping −0.7 … −3.5 %, exciting forces within 1.4 %; energy balance B33 = ω³|X3|²/(2ρg³) and B11 = ω³|X1|²/(4ρg³) within 2 %; box 20 × 8 × 2 m against Capytaine at 4–25 s: diagonal added mass within 3 % (roll 4 %), damping within 5 % of its largest value, exciting forces within 3 % (roll moment 6 %) and phases within 0.03 rad; convergence with the panel count; removal of the irregular frequency at 2.45 rad/s.

### Performance of the box model

Measured with `process.cpuUsage()` on a shared machine under load (load average 15–25 on 4 cores; a loop of 10⁶ multiply–add operations took 5–7 ms of processor time there, against about 1 ms on an unloaded desktop processor, which was not available for measurement). Hull 100 × 20 × 5 m, 127 panels on the quarter, 14 frequencies, three classes of modes:

| | Before (complex elimination at every frequency) | Now |
| --- | --- | --- |
| new hull, code warm (from about the sixth hull of a session) | 1.07–1.27 s | 0.50–0.56 s |
| first hull of a session (partly in code the engine has not yet optimised) | 3.1–3.9 s | 2.7 s |
| second to fifth hull | 1.2–1.9 s | 0.6–1.3 s |
| hull 20 × 8 × 2 m (111 panels), code warm | 0.73–0.81 s | 0.41–0.43 s |
| 200 hulls with L, B, T drawn independently within ±20 %, heave, pitch and roll at 40 frequencies each | 272 CPU-s (1.35 s per hull) | 113 CPU-s (0.55 s per hull) |
| same proportions at another scale, density or gravity | as a new hull | not solved again (Froude similarity) |
| same hull again, 40 frequencies, one heading | 30–80 ms | 30–60 ms |

Where the time went, and what was done (shares of a warm new hull before → after):
- complex elimination, 35 % → real elimination 15 % plus complex 7 %: the system is solved as K_r + iUVᵀ, with the outgoing-wave part e^{K(z+ζ)}J0(KR) written as a series of 5–15 products per class (Graf's addition theorem); the three or so frequencies that carry null-field rows keep the complex elimination;
- assembly of the wave part, 25 % → about 20 %: square roots, logarithms and reciprocals of the pair distances are kept with the hull, so that a frequency costs table look-ups and arithmetic only; numbers are passed between the hot functions in typed arrays (as arguments they are allocated on the heap, which had the garbage collector at 10–15 %);
- Rankine integrals, 15 % → 5 %: closed form within 2.5 panel diameters, 2 × 2 Gauss beyond (the coefficients change by less than 3·10⁻⁵);
- Haskind damping, 12 % → 1 %: heading integral in closed form from the same series;
- the table of the Green function is built once per session and lazily (about 40 ms for the range a barge needs); it was 2 % of a new hull and is not shipped as an asset.

The results differ from those of the complex elimination by less than 3·10⁻⁵ of the coefficients (the split solution sums the Bessel series where the other path interpolates J0 and J1 in the table); the comparison with Capytaine above is unchanged to the 0.1 % shown. The results are deterministic: the same hull gives the same numbers in any order of calls.

Not adopted: interpolating between solved neighbouring hulls in (L/B, B/T). Tried on 100 × 20 × 5 m with the four neighbours at a spacing h in ln(L/B) and ln(B/T), same mesh topology, potentials interpolated and the Haskind integral taken on the actual hull: h = 0.1 gave heave and pitch within 0.5 % for wave periods above twice the shortest solved one and 1–4 % below, roll −0.5 … −1 %; h = 0.25 gave up to 2 % in heave at the head-sea cancellation, 3–6 % in roll and 5–30 % near the wavelength equal to the beam. At the spacing that keeps the error below 1 % a 200-sample study touches about as many lattice points as half its samples, and a sensitivity study more points than samples, so the direct solution is used for every hull.

A new hull therefore costs about half a second of processor time on that machine once the code is warm, and several times that for the first hulls of a session; `tests/engine/hydro.radiation.test.js` guards the warm figure against a reference loop timed alongside.


## hydrodb.js
- `createHydroDb({ omega, A, B, headingsDeg, Fre, Fim, Ainf?, A0?, rho, g, depth }) -> db` — SI container, sorted.
- `parseWamit1(text, { rho, L, periodIs }) -> { omega, A, B, Ainf, A0 }` — A = ρL^k Ā, B = ρL^k ωB̄ (k = 3, 4, 5). WAMIT v7 User Manual.
- `parseWamit3(text, { rho, g, L, periodIs }) -> { omega, headingsDeg, Fre, Fim }` — X = ρgL^m X̄ (m = 2, 3); also reads `.2`.
- `parseWamit4(text, { L, periodIs }) -> { omega, headingsDeg, re, im }` — RAOs, m/m and rad/m.
- `hydroDbFromWamit({ text1, text3?, rho, g, L, depth, periodIs }) -> db`.
- `parseHydroJson(json) -> db`, `parseHydroCsv(text, { rho, g, depth }) -> db` — schemas documented at the top of the file.
- `hydroAt(db, omega, headingDeg?) -> { A, B, F: { re, im } }` — linear interpolation in ω and heading.
- `raoSolve({ db, M, C, Kext?, Bext?, omega?, headingDeg }) -> { omega, period, headingDeg, depth, g, re, im, amp, phase }` — m/m, rad/m.
- `checkHydroDb(db, { tol }) -> { ok, negativeDamping, asymmetryA, asymmetryB, haskind }` — B ≥ 0, symmetry, B_ii = k/(8πρgC_g)∫|X_i|²dβ. Newman (1977) ch. 6.

## timedomain.js
- `retardationKernel({ omega, B, dt, tMax }) -> { dt, t, n, K }` — K(t) = (2/π)∫B cos ωt dω, N/m. Cummins (1962).
- `infiniteFrequencyAddedMass({ omega, A, kernel }) -> { Ainf, spread }` — kg. Ogilvie (1964).
- `waveComponents({ omega, S, seed }) -> { omega, amp, phase }` — m, rad; seeded.
- `waveExcitationForce({ components, Fre, Fim }) -> f(t) -> number[]` — N.
- `simulateCummins({ M, Ainf?, Blin?, C, kernel?, force?(t, x, v), x0?, v0?, dt, nSteps }) -> { t, x, v }` — fixed-step predictor–corrector with trapezoidal convolution; constant-coefficient model when no kernel.
- `freeDecay({ mass, stiffness, b1, b2, kernel?, x0, dt, duration }) -> { t, x, v }`.
- `decayAnalysis({ t, x, removeMean, inertia? }) -> { period, peaks, logDecrement, dampingRatio, p1, p2, b1, b2, r2 }` — s, 1/s, 1/unit; p–q fit of linear + quadratic damping. Faltinsen (1990) ch. 7.

## drift.js
- `maruoDrift({ a, R, T?, depth, rho, g }) -> N/m` — ½ρga²R²(1 + 2kh/sinh 2kh). Maruo (1960); Longuet-Higgins (1977).
- `cylinderDriftShortWave({ a, D, rho, g }) -> { force, coefficient }` — N, N/m²; (2/3)ρga²r, short-wave upper bound.
- `meanDriftForce({ omega, S, D }) -> N` — 2∫S·D dω.
- `newmanDriftSpectrum({ omega, S, D, mu }) -> Float64Array` — N²·s/rad; 8∫S(ω)S(ω+μ)D(ω+μ/2)² dω. Newman (1974).
- `slowDriftResponse({ mass, stiffness, damping | dampingRatio, mu, SF, meanForce, duration }) -> { naturalPeriod, meanOffset, sigma, sigmaWhiteNoise, extreme, nCycles, damping }` — s, m.
- `ocimfForce({ kind: 'wind'|'current', V, angleDeg, Cx, Cy, Cxy, AT, AL, L, T, rho }) -> { Fx, Fy, Mz }` — N, N·m; coefficient curves are inputs. OCIMF (1994).
- `projectedAreaForce({ V, directionDeg, Afront, Aside, Cfront, Cside, xSide, yFront, rho }) -> { Fx, Fy, Mz }` — cross-flow principle.
- `meanEnvironmentalForce({ headingsDeg, loads: [{ name, directionDeg, force(relDeg) }] }) -> { headingDeg, Fx, Fy, Mz, F, worst, best }` — N, N·m.

## impact.js
- `wedgeWaterEntry({ deadriseDeg, V, t, halfBeam?, rho }) -> { vonKarman: { c, force, addedMass }, wagner: { c, force, addedMass, peakPressure, Cp }, ratio, chineWetted }` — m, N/m, kg/m, Pa. von Kármán (1929); Wagner (1932).
- `slammingPressureCoefficient({ shape: 'cylinder'|'flat'|'wedge', deadriseDeg? }) -> C_Pa` — 5.15 / 2π / 2.5/(tan β)^1.1. DNV-RP-C205 §8.
- `cylinderSlamCoefficient({ s, D }) -> Cs` — 5.15[D/(D + 19s) + 0.107s/D]. Campbell & Weynberg (1980).
- `slammingProbability({ m0r, m0v, d, L?, vThreshold?, g }) -> { probability, ratePerHour, vThreshold, emergenceProbability }` — Ochi (1964), threshold 0.093√(gL) by default.
- `deckWetnessProbability({ m0r, m2r?, freeboard }) -> { probability, ratePerHour }`.
- `damBreakDeckFlow({ h0, x, t, rho, g }) -> { h, u, frontSpeed, frontPosition, arrived, stagnationPressure, momentumFlux, hydrostaticPressure }` — m, m/s, Pa, N/m. Ritter (1892).
- `sloshingPeriods({ shape: 'rectangular'|'cylindrical', L | R, h, g }) -> { periods, omegas }` — s, rad/s; linear theory.
- `sloshingResonanceCheck({ sloshingPeriod, motionPeriods, band = 0.2 }) -> { ratios, resonant, any }` — band is a screening input.
- `waveInDeckForce({ crestElevation, deckBottom, deckTop?, width, u, U, Cd, waveKinematicsFactor, currentBlockageFactor, rho }) -> { force, area, elevation, inundation }` — N, m², m. API RP 2A-WSD §C17.6.2 (Cd is an input; table in the JSDoc).
- `vimScreening({ U, Tn, D, lockInRange = [4, 10], AoverD?, St }) -> { Ur, lockIn, sheddingPeriod, UrResonance, amplitude, criticalSpeeds }` — band and A/D are inputs.

## dp.js
- `allocateThrust({ thrusters: [{ name, x, y, type: 'azimuth'|'tunnel'|'fixed', Tmax, Tmin?, angleDeg?, forbidden? }], demand: { Fx, Fy, Mz }, method: 'bounded'|'pinv', disabled? }) -> { thrusters: [{ name, T, angleDeg, Tx, Ty, utilisation }], achieved, residual, feasible, maxUtilisation, demandFraction, iterations }` — N, N·m, deg. Fossen (2011) §12.3.
- `dpCapability({ thrusters, envForce(windSpeed, directionDeg), directionsDeg?, vMax, failures: 'single'|'none'|number[][], maxUtilisation }) -> { directionDeg, intact, failure, failureCase }` — limiting wind speed m/s per direction. IMCA M 140; DNV-ST-0111.
- `lqr({ A, B, Q, R }) -> { K, P, residual, iterations }` — continuous algebraic Riccati equation by the matrix sign function. Roberts (1980).
- `dpLqrGains({ M, D, Q, R }) -> { Kp, Kd, K, residual }` — 3×3 gains for the low-speed DP model.
- `dpPidGains({ M, D, omegaN, zeta }) -> { Kp, Ki, Kd }` — pole placement. Fossen (2011) §12.2.
- `simulateDp({ M, D, thrusters, setpoint, initial, gains, envForce?(t), waveMotion?, noiseStd, observer?, thrustRate, azimuthRateDeg, dt, duration, statsStart?, seed }) -> { t, x, y, psiDeg, footprint: { maxRadius, r95, stdX, stdY, stdPsiDeg, maxPsiDevDeg, meanX, meanY }, utilisation, saturatedFraction }` — m, deg; PID + steady-state Kalman wave filter + rate-limited thrusters.

## _wavehelp.js (private helpers)
- `G`, `dispersion(omega, depth, g)`, `depthRatios`, `prepWave`, `airyEta`, `airyPoint`, `jonswap(omega, { Hs, Tp, gamma })`, `spectrumMoments(omega, S)`, `csolve(Ar, Ai, br, bi)`, `boundedLeastSquares(A, b, lo, hi, opts)`, `needPos`, `needVec3` — not part of the public API; use `src/engine/metocean` for wave theory and spectra in new code.
