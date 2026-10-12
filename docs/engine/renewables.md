# `src/engine/aero` — wake/energy-yield, tidal, wave-energy, solar and storage engines

SI units throughout (m, kg, s, N, Pa, W); angles in degrees where the name ends in `Deg`, temperatures in °C where
the name ends in `C`. Energies are in Wh where the field name ends in `Wh`. Wind direction `wdDeg` is the direction
the wind comes from, clockwise from north (+y north, +x east); tidal `flowDirDeg` is the direction the current flows
toward. Wave spectra are one-sided in angular frequency, S(ω) in m²·s/rad. Solar azimuths are clockwise from north,
longitude east positive, times in ms since epoch (UTC). Tests: `tests/engine/aero.<file>.test.js`.

## wake.js
- `inductionFromCt(Ct) -> a` — axial induction ½(1 − √(1 − Ct)), Ct limited to 0.999. Actuator-disc momentum theory.
- `jensenDeficit({ Ct, x, D, k = 0.05 }) -> { deficit, wakeDiameter }` — (-), m; top-hat (1 − √(1 − Ct))/(1 + kx/R)². Jensen (1983); Katić et al. (1986).
- `frandsenDeficit({ Ct, x, D, alpha, n = 2 }) -> { deficit, wakeDiameter, beta }` — momentum-conserving top-hat with Dw = D(β^(n/2) + αx/D)^(1/n); α has no default. Frandsen et al. (2006) Wind Energy 9:39.
- `wakeGrowthFromTI(TI) -> kStar` — k* = 0.3837 TI + 0.003678. Niayifar & Porté-Agel (2016).
- `gaussianDeficit({ Ct, x, D, kStar, r = 0 }) -> { deficit, centreDeficit, sigma, epsilon, nearWake }` — (-), m; near-wake radicand clamped and flagged. Bastankhah & Porté-Agel (2014).
- `yawDeflection({ Ct, yawDeg, x, D, TI, ky, kz, alpha = 2.32, beta = 0.154 }) -> { deflection, theta0, x0, sigmaY, sigmaZ }` — m, rad; wake-centre deflection behind a yawed rotor. Bastankhah & Porté-Agel (2016).
- `meanderedGaussian({ centreDeficit, sigma, sigmaMeanderY, sigmaMeanderZ }) -> { centreDeficit, sigmaY, sigmaZ }` — fixed-frame mean of a meandering Gaussian wake (convolution). After Larsen et al. (2008).
- `addedTurbulenceFrandsen({ Ct, xOverD }) -> I_add` — 1/(1.5 + 0.8(x/D)/√Ct). Frandsen (2007); IEC 61400-1 Annex E.
- `addedTurbulenceCrespo({ Ct, xOverD, I0 }) -> { added, inRange }` — 0.73 a^0.8325 I0^0.0325 (x/D)^−0.32. Crespo & Hernández (1996).
- `effectiveTurbulence({ I0, neighbours: [{ distanceD, Ct }], m, pw = 0.06 }) -> { Ieff, wakeTI }` — fatigue-effective TI with Wöhler exponent m. IEC 61400-1:2019 Annex E.
- `superpose(deficits, method = 'rss') -> deficit` — 'linear' | 'rss' | 'energy' | 'max'. Lissaman (1979); Katić et al. (1986).
- `circleOverlapArea(R1, R2, d) -> area` — m²; lens formula.
- `tophatRotorAverage({ deficit, wakeRadius, rotorRadius, offset }) -> { deficit, overlapFraction }` — exact partial-wake average of a top-hat.
- `rotorPoints(nR = 3, nTheta = 8) -> { y, z, n }` — equal-area ring quadrature points on the unit disc.
- `gaussianRotorAverage({ centreDeficit, sigma, dy, dz, rotorRadius, points }) -> deficit` — rotor-disc average of a Gaussian wake by quadrature.
- `rotorEquivalentWindSpeed({ profile, hubHeight, D, nSeg = 40 }) -> { rews, hubSpeed, ratio }` — m/s; (Σ v³ A_i/A)^(1/3). IEC 61400-12-1:2017.
- `inductionOnAxis({ a, x, R }) -> u/U∞` — 1 − a(1 + ξ/√(1 + ξ²)), ξ = x/R (negative upstream). Vortex cylinder; Branlard & Gaunaa (2015).
- `inductionUpstream({ a, x, r, R }) -> Δu/U∞` — off-axis far-field induction (equivalent point source, capped by the on-axis value); zero for x > 0. Derived; cf. Gribben & Hawkes (2019).
- `deepArrayFrandsen({ Ct, sx, sy, hubHeight, z0, H, kappa = 0.4 }) -> { z00, ct, speedRatio, powerRatio }` — infinite-farm effective roughness (m) and hub-speed ratio for a reference height H. Frandsen (1992).
- `farmFlow({ turbines: [{ x, y, D, hubHeight, ct, power }], U, wdDeg, TI = 0.06, model = 'gaussian'|'jensen', k, kStar, superposition = 'rss', deficitReference = 'freestream'|'local', turbulenceModel = 'crespo'|'frandsen'|'none', reflect, blockage, maxIter, tol, points }) -> { U, TI, Ct, power, deficit, blockage, order, totalPower, grossPower, efficiency, iterations }` — m/s, -, -, W per turbine; steady farm solver, upstream → downstream sweep, optional image turbines and summed-induction blockage.

## aep.js
- `curvePower({ u, p }, U) -> P` — W; linear interpolation, zero outside the table.
- `cubicPowerCurve({ ratedPower, D, Cp = 0.45, efficiency = 0.94, rho = 1.225, cutIn = 3, cutOut = 25, smoothing = Infinity, step = 0.5 }) -> { u, p, cp, ct, ratedSpeed, ratedPower, D, rho }` — generic pitch-regulated curve; Ct from actuator-disc momentum theory.
- `powerCurveFromRotor({ rotor: U => { power, thrust }, speeds, D, rho, ratedPower, efficiency }) -> { u, p, cp, ct, thrust, D, rho }` — curve from a rotor/BEM callback.
- `parseTurbineCsv(text, { D?, rho? }) -> { u, p, cp, ct, thrust }` — m/s, W, -, -, N; NREL turbine-models CSV layout (wind speed, power, Cp, thrust, Ct with bracketed units).
- `airDensityCorrection(curve, { rho, rho0 = 1.225, regulation = 'pitch'|'stall' }) -> curve` — wind-speed scaling (ρ/ρ0)^(1/3) or power scaling. IEC 61400-12-1:2017.
- `airDensity({ pressure, temperatureC, relativeHumidity = 0 }) -> rho` — kg/m³; moist air. IEC 61400-12-1.
- `weibull(A, k) -> { pdf, cdf, mean, std, cubeMean, moment }` — m/s; mean = A Γ(1 + 1/k).
- `weibullScaleFromMean(mean, k = 2) -> A` — m/s.
- `extrapolateWind({ U, z1, z2, alpha | z0 }) -> U2` — m/s; power law or neutral log law. IEC 61400-1 / 61400-3-1 wind profile.
- `aepWeibull({ curve, A, k = 2, uMin, uMax, n = 4000, hoursPerYear = 8760 }) -> { aepWh, meanPower }` — Simpson integration of P(u) f(u).
- `aepBins({ curve, A, k = 2, hoursPerYear }) -> { aepWh, meanPower }` — method of bins, N_h Σ [F(V_i) − F(V_(i−1))](P_(i−1) + P_i)/2. IEC 61400-12-1:2017.
- `normaliseRose(sectors: [{ dirDeg, frequency, A, k }]) -> sectors` — frequencies scaled to sum to 1.
- `farmAep({ turbines, rose, uMin = 0.5, uMax = 30, du = 1, subSectors = 1, wake = {}, hoursPerYear }) -> { grossWh, aepWh, arrayEfficiency, wakeLoss, perTurbine, perSector, ratedPower, capacityFactor }` — sector × speed-bin loop over `farmFlow`.
- `lossCascade({ grossWh, losses, ratedPower, hoursPerYear }) -> { netWh, totalLoss, steps, capacityFactor, fullLoadHours }` — multiplicative loss chain (wake, blockage, electrical, availability, environmental, curtailment).
- `normalQuantile(p) -> z` — standard normal quantile.
- `aepUncertainty({ p50Wh, windSpeed, sensitivity, energy, interAnnual, years, levels = [75, 90, 99] }) -> { sigma, sigmaWh, components, p50Wh, p75Wh, p90Wh, p99Wh, exceedance }` — root-sum-square of independent components; Pxx = P50(1 − z_xx σ), Pxx exceeded with probability xx %.
- `mcp({ site, reference, longTermReference, method = 'linear'|'varianceRatio' }) -> { slope, intercept, r, r2, predict, longTermMean, concurrentMean, correctionFactor, predicted }` — measure–correlate–predict. Rogers et al. (2005).
- `timeSeriesProduction({ speeds, curve, dtHours = 1, times, lossFactor, ratedPower }) -> { power, energyWh, meanPower, capacityFactor, monthly, daily, diurnal }` — W, Wh; UTC monthly/daily/diurnal aggregation.
- `capacityValue({ generation, demand, ratedPower, topFraction = 0.1 }) -> { capacityValue, steps }` — peak-period capacity factor (capacity-credit proxy).
- `turbulenceSummary({ mean, std, binWidth = 1, minSpeed = 3 }) -> { meanTI, bins, ti15 }` — binned TI; representative = mean + 1.28 σ. IEC 61400-1.
- `shearSummary({ u1, u2, z1, z2, minSpeed }) -> { alpha, alphaMeanOfSamples, z0, n }` — power-law exponent and log-law roughness (m).
- `extremeWind({ annualMaxima, returnPeriod = 50, gustFactor = 1, method = 'mom' }) -> { value, loc, scale, returnPeriod }` — m/s; Gumbel fit (reliability/extreme.js), V_T = μ − β ln(−ln(1 − 1/T)).
- `generationProfile({ power, dtHours = 1, prices, exportLimit }) -> { energyMWh, curtailedMWh, totalMWh, revenue, capturePrice }` — revenue-ready profile.

## layout.js
- `polygonArea(polygon) -> area` — signed m², shoelace formula.
- `pointInPolygon(pt, polygon) -> boolean` — ray casting.
- `distanceToBoundary(pt, polygon) -> d` — m.
- `checkLayout({ points, polygon, exclusions, minSpacing, setback }) -> { feasible, outside, excluded, spacingViolations, minDistance, penalty }` — boundary, exclusion and spacing constraints.
- `gridLayout({ polygon, D, sx, sy, rowAngleDeg = 0, stagger, exclusions, setback, origin, maxCount }) -> { points, count }` — regular/staggered grid clipped to the site; spacings in rotor diameters.
- `rowAngleForWind(wdDeg) -> rowAngleDeg` — row direction perpendicular to a wind direction.
- `boundaryLayout({ polygon, nBoundary, setback, infill, minSpacing }) -> { points, boundaryCount, count }` — perimeter positions with optional grid infill.
- `makeAepObjective({ turbine, rose, aep }) -> points => Wh` — annual-energy objective from `farmAep`.
- `optimiseLayout({ initial, objective, polygon, minSpacing, exclusions, setback, method = 'pattern'|'random'|'pso', maxEvals, step, minStep, seed, particles }) -> { points, value, initialValue, improvement, history, evals, feasible }` — seeded constrained layout optimisation with improvement history.
- `spacingSweep({ grid, sxValues, syValues, objective }) -> [{ sx, sy, count, value, valuePerTurbine }]` — spacing sensitivity.
- `orientationSweep({ grid, anglesDeg, objective }) -> { results, best }` — row-orientation sensitivity.
- `capacityDensity({ points, ratedPower, polygon }) -> { installed, hullArea, hullDensity, siteArea, siteDensity }` — W, m², W/m² (≡ MW/km²).
- `minimumSpanningTree(points, root = 0) -> { edges, length }` — m; Prim's algorithm.
- `cableLayout({ turbines, substation, maxPerString, routeFactor }) -> { edges, length, strings, mstLength }` — m; capacitated strings by the Esau–Williams (1966) heuristic, MST when unconstrained.

## tidal.js
- `TIDAL_SPEEDS_DEG_PER_HOUR` — M2, S2, N2, K2, K1, O1, P1, Q1, M4 angular speeds (deg/h).
- `constituentPeriod(name) -> { periodHours, period, omega }` — h, s, rad/s.
- `synthesizeTide({ constituents: [{ name, amplitude, phaseDeg }], t | duration, dt, mean }) -> { t, u }` — s, m/s; harmonic synthesis without nodal corrections.
- `springNeap({ springPeak, neapPeak }) -> { constituents, springNeapPeriod, envelope }` — M2 + S2 estimator; period in s.
- `velocityDistribution(u, { binWidth = 0.1 }) -> { speed, exceedance, bins, mean, max, cubeRootMeanCube }` — exceedance curve and histogram of |u|.
- `currentProfile({ model = 'power'|'log', uRef, zRef, n = 7, z0 }) -> z => u` — m/s at height z above the bed. Soulsby (1997).
- `rotorAveragedVelocity({ profile, hubHeight, D, nSeg }) -> { mean, cubeMean, powerVelocity, hub }` — disc averages of u and u³.
- `actuatorRotor({ a = 1/3, efficiency = 0.75 }) -> U => { Cp, Ct }` — ideal actuator-disc fallback rotor.
- `rotorFromSolver({ solve, rotor, D, tsr, pitchDeg, rho }) -> U => { Cp, Ct }` — adapter for an injected BEM solver `solve({ rotor, U, omega, pitchDeg, rho })`.
- `tidalPowerCurve({ D, ratedPower, rotor, rho = 1025, cutIn = 0.7, uMax = 5, step = 0.1, efficiency, parkedCt }) -> { u, p, cp, ct, thrust, ratedSpeed, D, rho }` — power and thrust curve with rated cap.
- `cavitationCheck({ U, tipSpeed, D, hubDepth, cpMin, rho, patm, pv = 1700, a = 0, margin = 1 }) -> { sigma, sigmaRequired, cavitates, relativeSpeed, tipDepth, maxTipSpeed, utilisation }` — tip cavitation number at top dead centre vs −Cp,min. Batten et al. (2008).
- `channelBlockage({ B, alpha4, Fr = 0 }) -> { alpha2, alpha4, beta4, Ct, Cp, h4OverH, dhOverH, CtLocal, CpLocal }` — actuator disc in a blocked open channel (1-D momentum). Garrett & Cummins (2007); Houlsby et al. (2008); Whelan et al. (2009).
- `channelBlockageLimit({ B, Fr = 0 }) -> { CpMax, alpha4, Ct, alpha2, beta4, garrettCummins, dhOverH }` — maximum Cp; (16/27)(1 − B)^−2 for Fr = 0.
- `blockageCorrection({ Ct, Cp, B, Fr = 0 }) -> { speedRatio, CtFree, CpFree, alpha2, alpha4 }` — channel → open-water coefficients (equal thrust and disc velocity). Barnsley & Wellicome (1990); Bahaj et al. (2007).
- `tidalArrayFlow({ turbines, U, flowDirDeg, TI = 0.1, wake }) -> farmFlow result` — Gaussian wakes with seabed reflection.
- `bidirectionalYield({ floodDirDeg, ebbDirDeg, headingDeg, floodEnergy, ebbEnergy, exponent = 3 }) -> { headingDeg, floodMisalignmentDeg, ebbMisalignmentDeg, floodFactor, ebbFactor, energy, loss, asymmetryDeg }` — fixed-yaw cos^p loss with flow asymmetry.
- `tidalEnergy({ curve, speeds | distribution, availability, velocityFactor, hoursPerYear, ratedPower }) -> { aepWh, meanPower, capacityFactor, ratedFraction }` — annual energy from a series or distribution.
- `extremeThrust({ D, Ct, springPeak, surgeCurrent, wave: { H, T, depth }, hubHeight, rho, duration }) -> { designVelocity, waveVelocity, thrust, meanThrust, thrustRange, cycles }` — N; spring + surge + linear-wave orbital velocity; wave-induced thrust range and cycle count.
- `supportLoads({ thrust, hubHeight, columnDiameter, profile, columnTop, Cd = 1, Cm = 2, wave, rho, nSeg }) -> { baseShear, baseMoment, columnDrag, columnDragMoment, inertiaShear, inertiaMoment }` — N, N·m about the seabed; thrust + Morison on the column.
- `kineticPowerDensity(u, rho = 1025) -> { mean, max, meanSpeed, maxSpeed, cubeRootMeanCube }` — W/m²; ½ρ|u|³ statistics.
- `fluxResource({ meanPowerDensity, area, sif, hoursPerYear }) -> { fluxPower, extractablePower, extractableWh }` — flux method with significant impact factor. Black & Veatch (2005).
- `tidalRangeEbb({ basinArea, range, turbineArea, sluiceArea, period, cdTurbine, cdSluice, startHead, minHead, efficiency, ratedPower, rho, dt, cycles }) -> { energyPerCycleWh, potentialPerCycleWh, efficiency, meanPower, peakPower, annualWh, generatingFraction, t, sea, basin, power }` — 0-D ebb-generation barrage/lagoon model; potential ½ρgAR². Prandle (1984).

## wec.js
- `dispersion(omega, depth = Infinity) -> { k, wavelength, c, cg }` — linear dispersion.
- `wavePower({ Hm0, Te, depth, rho = 1025 }) -> J` — W/m; ρg²Hm0²Te/(64π) deep, ρg(Hm0²/16)c_g(Te) finite depth.
- `regularWavePower({ amplitude, omega, depth, rho }) -> J` — W/m; ½ρga²c_g.
- `waveSpectrum({ Hs, Tp, gamma = 1, omega }) -> S` — m²·s/rad; Bretschneider (γ = 1) / JONSWAP. DNV-RP-C205 §3.5.5.
- `spectralPeriods({ omega, S }) -> { Hm0, Te, Tz, m0 }` — Te = m−1/m0.
- `scatterResource({ hs, te, prob, depth, rho }) -> { meanJ, annualWhPerM, J, contribution }` — mean annual wave power from a scatter diagram prob[iHs][iTe].
- `monthlyResource({ hs, te, monthly, depth, rho }) -> J[]` — W/m per month.
- `cylinderHeave({ radius, draft, omega, depth, rho, addedMassCoeff = 1/6, diffractionFactor = 1 }) -> { A, B, K, Fe, k, cg, displacedMass }` — truncated-cylinder heave approximations: K = ρgA_w, A = ρD³/6, Froude–Krylov excitation, Haskind damping B = k|Fe|²/(4ρg c_g). Newman (1962); Falnes (2002).
- `pointAbsorber({ mass, A, B, K, Fe, omega, amplitude, control = 'passive'|'reactive'|'fixed', Bpto, Kpto, Bvisc, limits: { stroke, velocity, force } }) -> { power, stroke, velocity, ptoForce, Bpto, Kpto, constrained, feasible, powerUnconstrained, powerMax }` — 1-DOF frequency-domain absorber; reactive optimum |Feζ|²/(8B). Falnes (2002).
- `captureWidth({ power, J, width }) -> { captureWidth, ratio }` — m, -.
- `heaveLimits({ omega, H, volume, depth, rho }) -> { maxCaptureWidth, radiationLimit, budalLimit, limit }` — λ/2π limit and Budal bound (π/4)ρgVH/T. Budal & Falnes (1975, 1980).
- `cylinderBody({ radius, draft, mass, Bvisc, depth, rho, ... }) -> omega => { M, B, K, Fe }` — body callback for the spectral functions.
- `spectralResponse({ omega, S, body, Bpto, Kpto }) -> { power, sigmaX, sigmaV, sigmaF, Tz }` — irregular-wave response by superposition; P̄ = B_pto σ_v².
- `powerMatrix({ hs, tp, body, Bpto = 'optimal', Kpto, gamma, ratedPower, efficiency, omega }) -> { hs, tp, power, Bpto }` — W, power[iHs][iTp] with rating cap.
- `matrixEnergy({ power, prob, ratedPower, availability, arrayFactor, devices, hoursPerYear }) -> { aepWh, meanPower, capacityFactor }` — annual energy from power matrix × scatter diagram.
- `radiationKernel({ omega, B, dt, tMax }) -> { dt, K }` — N/m; K(t) = (2/π)∫B cos ωt dω. Cummins (1962); Ogilvie (1964).
- `infiniteAddedMass({ A, omega, kernel }) -> Ainf` — kg; Ogilvie relation.
- `simulateCummins({ mass, Ainf, stiffness, kernel, excitation, dt, duration, x0, v0, pto: { type: 'linear'|'coulomb'|'hydraulic', ... }, mooringStiffness, dragCoeff, endStop, latching, declutch, extraForce, kernelLength }) -> { t, x, v, ptoForce, pressure, absorbedEnergy, meanAbsorbedPower, meanGeneratedPower, maxStroke, maxVelocity, maxPtoForce }` — RK4 time-domain model with radiation convolution, nonlinear PTO, end stops, latching and declutching.
- `twoBodyAbsorber({ omega, amplitude, body1, body2, A12, B12, Bpto, Kpto, optimise }) -> { power, Bpto, x1, x2, relative, ptoForce }` — float + reaction body, PTO on relative heave. Falnes (1999).
- `flapHydrodynamics2D({ omega, depth, rho }) -> { excitation, damping, waveRatio, k, cg }` — per metre width; bottom-hinged thin flap from wavemaker theory + 2-D Haskind. Dean & Dalrymple (1991).
- `flapStiffness({ volume, buoyancyArm, mass, gravityArm, rho }) -> C` — N·m/rad; g(ρV z_b − m z_g).
- `flapConverter({ width, depth, omega, amplitude, inertia, addedInertia, stiffness, Bpto = 'optimal', Bvisc, widthFactor, maxAngle, rho }) -> { power, angle, Bpto, ptoTorque, captureWidth, captureWidthRatio, excitationMoment, radiationDamping, constrained }` — 1-DOF pitch model of an oscillating wave surge converter.
- `captureWidthDevice({ J, captureWidthRatio, length, efficiency, ratedPower }) -> P` — W; parameterised attenuator/terminator (user-supplied CWR).
- `overtoppingPower({ q, crestFreeboard, width, efficiency, rho }) -> P` — W; ηρg q R_c w with user-supplied overtopping rate.
- `arrayPower({ devices, power, q = 1 }) -> P` — W; interaction factor q as an argument.
- `survivalResponse({ omega, S, body, Bpto, Kpto, duration = 10800, strokeLimit }) -> { maxStroke, maxVelocity, maxPtoForce, sigmaX, cycles, strokeUtilisation }` — most probable maxima σ√(2 ln N) in an extreme sea state.
- `slammingForce({ Cs, area, velocity, rho }) -> F` — N; ½ρC_sAv² hook. DNV-RP-C205 §8.6.

## owc.js
- `columnResonance({ draft, addedLength = 0 }) -> { period, omega }` — T = 2π√((d + l_a)/g).
- `chamberCompliance({ volume, patm = 101325, gamma = 1.4 }) -> C` — m³/Pa; V0/(γ p_atm). Sarmento & Falcão (1985).
- `columnBody({ radius, draft, addedLength, depth, rho, excitationFactor, Bvisc }) -> { area, column }` — rigid-piston water column coefficients.
- `owcRegular({ column, area, omega, amplitude, compliance, turbine = 'wells'|'impulse', conductance, ki, efficiency }) -> { pneumaticPower, electricalPower, pressure, flow, stroke, conductance, efficiency }` — W, Pa, m³/s, m; piston + linearised isentropic chamber + turbine (impulse by harmonic linearisation).
- `optimalTurbineDamping({ column, area, omega, amplitude, compliance, depth, rho }) -> { conductance, damping, pneumaticPower, captureWidth }` — optimum Wells turbine damping (Pa·s/m³).
- `owcIrregular({ column, area, omega, S, compliance, turbine, conductance, ki, efficiency }) -> { pneumaticPower, electricalPower, sigmaP, sigmaQ, sigmaX, conductance, efficiency }` — spectral solution; impulse turbine by statistical linearisation.
- `owcPowerMatrix({ hs, tp, column, area, compliance, turbine, conductance, ki, efficiency, gamma, ratedPower, prob, availability, omega }) -> { power, aepWh, meanPower, capacityFactor }` — power matrix and annual energy.

## pv.js
- `solarPosition({ time, lat, lon }) -> { zenithDeg, elevationDeg, apparentElevationDeg, azimuthDeg, declinationDeg, equationOfTimeMin, hourAngleDeg, distanceAU }` — NOAA Solar Calculator algorithm (Meeus).
- `sunriseSunset({ time, lat, lon, horizonZenithDeg = 90.833 }) -> { sunrise, sunset, solarNoon, dayLengthHours, polar }` — ms epoch, h.
- `extraterrestrialIrradiance(time, solarConstant = 1361) -> E0n` — W/m²; G_sc/R².
- `airMass(zenithDeg, pressure = 101325) -> { relative, absolute }` — Kasten & Young (1989).
- `clearSkyHaurwitz(zenithDeg) -> GHI` — W/m²; 1098 cos z exp(−0.059/cos z). Haurwitz (1945) per Reno et al. (2012).
- `clearSkyIneichen({ zenithDeg, linkeTurbidity, altitude, dniExtra, pressure }) -> { ghi, dni, dhi }` — Ineichen & Perez (2002).
- `erbsDecomposition({ ghi, zenithDeg, dniExtra }) -> { kt, kd, dhi, dni }` — Erbs, Klein & Duffie (1982).
- `erbsDiffuseFraction(kt) -> kd` — the Erbs correlation.
- `angleOfIncidence({ zenithDeg, azimuthDeg, tiltDeg, surfaceAzimuthDeg }) -> { cosAoi, aoiDeg }`.
- `planeOfArray({ dni, dhi, ghi, zenithDeg, azimuthDeg, tiltDeg, surfaceAzimuthDeg, albedo = 0.08, model = 'haydavies'|'isotropic', dniExtra }) -> { total, beam, skyDiffuse, groundReflected, aoiDeg, circumsolar }` — W/m². Liu & Jordan (1963); Hay & Davies (1980).
- `rearIrradiance({ ghi, dhi, tiltDeg, albedo, groundShadeFactor }) -> G_rear` — W/m²; isotropic view-factor bifacial model.
- `iamAshrae(aoiDeg, b0 = 0.05) -> IAM` — 1 − b0(1/cos θ − 1).
- `iamPhysical(aoiDeg, { n = 1.526, K = 4, L = 0.002 }) -> IAM` — Fresnel + Bouguer. De Soto et al. (2006).
- `cellTempNoct({ poa, tempAirC, noctC = 45 }) -> Tc` — °C.
- `cellTempFaiman({ poa, tempAirC, windSpeed, u0 = 25, u1 = 6.84, waterCoolingOffsetC = 0 }) -> Tm` — °C. Faiman (2008); floating-PV coefficients as arguments.
- `cellTempSandia({ poa, tempAirC, windSpeed, a = -3.47, b = -0.0594, deltaT = 3 }) -> Tc` — °C. King et al. (2004).
- `lambertW(x) -> W0` — principal branch, x ≥ −1/e.
- `singleDiodeCurrent(params, V) -> I` / `singleDiodeVoltage(params, I) -> V` — explicit Lambert-W solutions of the five-parameter model; params { IL, I0, a, Rs, Rsh }. Jain & Kapoor (2004).
- `singleDiodeCurve(params, nPoints = 0) -> { isc, voc, vmp, imp, pmp, fillFactor, v, i }` — characteristic points and sampled curve.
- `singleDiodeAt(ref, { G, tempCellC }) -> params` — irradiance/temperature translation. De Soto et al. (2006).
- `extractSingleDiode({ Isc, Voc, Imp, Vmp, Ns, alphaIsc, betaVoc | n }) -> { IL, I0, a, Rs, Rsh, n, alphaIsc, Ns, residualBetaVoc }` — five parameters from datasheet values. De Soto et al. (2006).
- `pvwattsDc({ pdc0, poa, tempCellC, gammaPdc = -0.0035 }) -> Pdc` — W. Dobos (2014).
- `seriesStringCurve({ groups: [{ params, fraction, G, tempCellC }], diodeDrop = 0.5, nPoints }) -> { i, v, pmp, vmp, imp, isc, voc, sumOfGroupPmp, mismatchLoss }` — partial shading with bypass diodes (I–V summation).
- `parallelStrings(strings, nPoints) -> { v, i, pmp, vmp, imp, sumOfStringPmp, mismatchLoss }` — parallel aggregation with blocking diodes.
- `inverterPvwatts({ pdc, pac0, etaNom = 0.96 }) -> { pac, efficiency, clipped }` — Dobos (2014).
- `inverterSandia({ pdc, vdc, inverter: { Paco, Pdco, Vdco, Pso, C0, C1, C2, C3, Pnt } }) -> { pac, clipped }` — King et al. (2007).
- `inverterQuadratic({ pdc, pacRated, p0, k1, k2 }) -> { pac, efficiency, clipped }` — quadratic loss model. Schmidt & Sauer (1996).
- `europeanEfficiency(efficiencyAt) -> etaEU` — 0.03 η5 + 0.06 η10 + 0.13 η20 + 0.10 η30 + 0.48 η50 + 0.20 η100.
- `waveTiltEffect({ amplitudeDeg, distribution = 'harmonic'|'gaussian', ...planeOfArray inputs }) -> { mean, still, factor }` — POA averaged over a wave-induced tilt oscillation.
- `simulatePv({ times, ghi, tempAirC, windSpeed, lat, lon, tiltDeg, surfaceAzimuthDeg, pdc0, dtHours, albedo, transposition, iam, soiling, bifaciality, groundShadeFactor, thermal, thermalParams, gammaPdc, dcLoss, dcAcRatio, inverterEfficiency, acLoss, availability, degradationPerYear, year, tiltAmplitudeDeg, solarConstant }) -> { pac, pdc, poa, tempCellC, energyWh, dcEnergyWh, poaWhPerM2, ghiWhPerM2, clippedWh, specificYield, performanceRatio, capacityFactor, monthly, daily }` — full energy chain from a GHI/temperature/wind series.
- `syntheticYear({ monthlyGhi, monthlyTempC, diurnalRangeC, lat, lon, year, solarConstant }) -> { times, ghi, tempAirC, clearness }` — hourly year from monthly climatology (constant clearness index per month). Liu & Jordan (1960).
- `optimalTilt({ system, tiltMin, tiltMax, tiltStep, azimuthsDeg }) -> { tiltDeg, surfaceAzimuthDeg, energyWh, table }` — tilt/azimuth search on annual energy.

## floatingpv.js
- `floatBuoyancy({ floatVolume, floatHeight, deadMass, liveMass, snowMass, rho = 1000, minFreeboard }) -> { cases, waterplaneArea }` — draft, freeboard and reserve buoyancy for dead / +live / +snow cases. Archimedes; DNV-RP-0584.
- `islandStability({ mass, kg, draft, length, width, fillRatio | floats, windForce, windLever, personMass, personOffset, rho }) -> { waterplaneInertia, BM, KB, GM, heelingMoment, heelDeg, stable }` — initial stability and heel of a float-island segment.
- `panelWindLoad({ windSpeed, panelArea, tiltDeg, cN, rows, shelterFactor, rho = 1.225 }) -> { dynamicPressure, normalForcePerRow, drag, lift, dragFirstRow, liftFirstRow }` — N; c_N and shelter factor are mandatory code/wind-tunnel inputs (EN 1991-1-4 §7.3; ASCE 7).
- `floatWaveLoad({ H, T, depth, submergedVolume, projectedArea, Cm = 2, Cd = 1, rho }) -> { inertia, drag, max, velocity, acceleration, wavelength }` — N; long-wave Morison force on a small float. DNV-RP-C205 §6.
- `meanDriftForce({ Hs, width, reflection, regular, rho }) -> F` — N; ρgR²Hs²/16 per metre (reflection coefficient as argument). Maruo (1960).
- `currentDrag({ speed, area, Cd, rho }) -> F` — N.
- `combineLoads(loads: [{ force, dirDeg }]) -> { Fx, Fy, magnitude, dirDeg }` — vector sum of environmental loads.
- `mooringLoadSharing({ Fx, Fy, lines: [{ azimuthDeg, stiffness, pretension }] }) -> { offsetX, offsetY, offset, tensions, maxTension, slackLines, loadedLines, sharing }` — rigid-island elastic load sharing with slack-line removal.
- `mooringLevelRange({ anchorDistance, depthMin, depthMax, load, allowableOffset }) -> { length, slackAtLowWater, excursionAtLowWater, elasticStrain, requiredStiffness }` — geometric slack for a water-level range and stiffness requirement.
- `anchorCount({ load, anchorCapacity, loadFactor, sharing }) -> { count, loadPerAnchor, utilisation }`.
- `hingedRaftResponse({ H, wavelength, moduleLength, moduleWidth, fillRatio, rho }) -> { waveSlope, modulePitch, relativeRotation, relativeRotationDeg, mismatch, connectorShear, lengthRatio, regime }` — rad, m, N; quasi-static articulated-raft screening (derived in the JSDoc).
- `evaporationSaving({ evaporationRate, coveredArea, reduction }) -> m³/yr` — user-supplied evaporation rate and reduction.
- `coolingYieldGain({ deltaTempC, gammaPdc }) -> gain` — γ ΔT_c.
- `floatingCableLength({ horizontalDistance, offset, levelRange, slackFactor = 1.1 }) -> { length, slack }` — m.
- `exposureClass({ Hs, thresholds = [1, 3], labels }) -> { class, index }` — site screening by Hs with editable thresholds.
- `FLOATING_PV_ELECTRICAL_NOTES` — qualitative cable/earthing reminders.

## storage.js
- `batteryStep({ request, soc, dtHours, battery: { capacityWh, cRate | powerW, roundTripEfficiency | etaCharge/etaDischarge, socMin, socMax } }) -> { power, soc, storedDeltaWh, lossWh }` — W (+ discharge / − charge), Wh.
- `runBattery({ request, battery, dtHours }) -> { power, soc, chargedWh, dischargedWh, lossWh, storedDeltaWh, equivalentFullCycles }` — battery over a request series; E_in η_c = E_out/η_d + ΔE.
- `batteryDegradation({ years, cyclesPerYear, calendarFadePerYear, cycleFadePerCycle, endOfLife = 0.8 }) -> { capacityFraction, lifeYears, lifeCycles }` — linear calendar + cycle fade (rates as arguments).
- `electrolyser({ power, ratedPower, specificConsumption, consumptionFactor, minLoad = 0.1, waterMargin = 1.5, dtHours }) -> { hydrogenKg, waterLitres, energyUsedWh, energyUnusedWh, capacityFactor, meanSpecificConsumption, efficiencyLhv, onFraction, hydrogenRate }` — kg, L; stoichiometric water 8.94 kg/kg × margin.
- `dispatchFirm({ generation, target, battery, dtHours, exportSurplus, tolerance }) -> { delivered, soc, battery, firmness, firmnessWithoutStorage, energyShiftedWh, shortfallWh, surplusWh, equivalentFullCycles, deliveredWh }` — firm-power dispatch.
- `dispatchCurtailment({ generation, exportLimit, battery, dtHours }) -> { exported, soc, curtailedWithoutWh, curtailedWithWh, curtailmentAvoidedWh, energyShiftedWh, exportedWh, equivalentFullCycles, cableUtilisation }` — curtailment absorption at an export limit.
- `dispatchArbitrage({ prices, battery, dtHours, method = 'threshold'|'lp', lowQuantile, highQuantile, lpMaxSteps = 96 }) -> { power, soc, revenue, chargedWh, dischargedWh, equivalentFullCycles, method }` — price arbitrage; LP via math/optim `linprog` with greedy fallback.
- `complementarity({ series, rated, exportLimit }) -> { capacityFactors, combinedCapacityFactor, correlation, cov, combinedCov, cableUtilisation, curtailedFraction, combined }` — hybrid wind/solar/wave statistics on a shared export cable.
