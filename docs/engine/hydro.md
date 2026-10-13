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
- `boxBargeRAO({ L, B, T, omega, headingDeg, speed, Cb, GMt?, kxx?, kxxDry?, rollDampingRatio, rollRadiation, rho, g }) -> { omega, heave, pitch, roll, rollPhase, heaveNaturalPeriod, rollNaturalPeriod, rollAddedInertia, rollRadiationDamping, validPeriodMin: { heave, roll } }` — m/m, rad/m; heave and pitch by Jensen, Mansour & Olsen (2004) with the bottom-pressure factors exp(−kT) and sin(½k_yB)/(½k_yB) at every heading (the paper's exp(−|k cos β|T) is exact in head seas only); roll by a 1-DOF model (not from the paper): Froude–Krylov moment, relative-motion diffraction −ω²A44·(wave slope at the keel), strip added inertia A44 from DNV-RP-C205 Table D-1 (used for the inertia when `kxxDry` is given; `kxx` includes added inertia), Haskind radiation damping, `rollDampingRatio` as damping in addition to it.

Validity of the closed-form floater models, from the comparison with Capytaine (`src/spec/cross-solver.json`, `docs/interop.md` §1.4):

| Model | Within 5 % of the panel solution | Outside |
| --- | --- | --- |
| `boxBargeRAO` heave, head and beam seas | wave periods ≥ 2 × the heave natural period 2π√(2T/g) (`validPeriodMin.heave`) | low: −15 % at 1.5 ×, −30 … −40 % near it (added mass taken equal to the displacement, two-dimensional sectional damping) |
| `boxBargeRAO` pitch, head seas | periods ≥ 1.5 × the heave natural period | −19 % at 1.25 ×, −63 % at 1.0 × (wavelength ≈ L) |
| `boxBargeRAO` roll, beam seas | periods ≥ 1.3 × the roll natural period (`validPeriodMin.roll`) | within 15 % below; near resonance the response follows the natural period (strip added inertia 6–10 % above the panel value for L/B = 2.5: +1.3 % in period) and the added damping; no sway–roll coupling |
| `floaterRAO` forces, single column | D/λ < 0.2 (`validOmegaMax`): surge and heave force from 5–6 s for D = 8 m | surge force +6 % at D/λ = 0.2, +47 % at 0.32 (diffraction) |
| `floaterRAO` motions, single column | heave at all periods above the limit, resonance included, with `end1: 'hemisphere'`; surge and pitch outside 0.9–1.4 × the pitch natural period | within 25 % near the pitch resonance: the peak has radiation damping only and follows the natural period (1.2 % high with `free1`, 3 % without) |

The one hull of each kind that was compared: barge 20 × 8 × 2 m (L/B = 2.5, B/T = 4), column D = 8 m with 40 m draught.


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
