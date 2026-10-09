# `src/engine/lines` — rod, riser, VIV and installation engines

SI units throughout (m, kg, s, N, Pa, rad; `...Deg` arguments in degrees). x horizontal, z up. All rod-based
analyses are planar (2-D): out-of-plane current, torsion, cross offsets and lateral seabed friction are not modelled.
One line per export: signature → returned fields — units — reference.

## rod.js

- `solveRodStatic(model) -> { converged, iterations, steps, outerIterations, residual, s, x, z, angle, Te, Tw, shear, moment, curvature, strain, reaction, friction, endA, endB, state }` — model: `sections[{ length, EA, EI, w, nElem?, D?, Cdn?, Cdt?, Ai?, Ae?, m?, ma?, maT?, contactOffset? }]`, `endA/endB { x?, z?, angle?, kRot?, Fx?, Fz?, M? }`, `pointLoads`, `seabed { z, kn, mu?, kt?, profile? }`, `supports[{ xc, zc, R, k, side?, angMin?, angMax? }]`, `current`, `pressure`, `initial ('straight' | 'cable' | { x, z, angle? })`, `nSteps`, `rampWeight`, `rampLoads`, `rampSupports`, `tol`, `maxIter`; lengths m, forces N, moments N·m, curvature 1/m, reactions N/m, angles rad — corotational Euler–Bernoulli beam (Crisfield 1991, Vol. 1 Ch. 7), effective tension (Sparks 2007 Ch. 2–3), Morison drag (DNV-RP-C205 §6.2).
- `rodModes(result, { nModes = 6 }) -> { s, frequenciesHz, modes[{ fHz, omega, periodS, phi, axial, dx, dz, curvature }], converged }` — Hz, rad/s, s; phi normalised to 1; modal curvature 1/m² — tangent stiffness + consistent mass with added mass, subspace iteration (Bathe §11.6).
- `rodHarmonicResponse(result, { omega, end, xAmp, zAmp, phaseZ, rotAmp, rho, zeta, tol, maxIter }) -> { s, xAmp, zAmp, normalAmp, tensionAmp, momentAmp, curvatureAmp, iterations, converged }` — m, N, N·m, 1/m — frequency-domain response with harmonically linearised drag, c = (8/3π)·½ρ·Cd·D·|v̂| (DNV-RP-C205 §6.2; Faltinsen §7).
- `cableShape(model) -> { s, x, z, angle }` — m, rad — lumped-chain catenary between two pinned ends, optional flat seabed (Irvine §1.2).

## riser.js

- `pipeProperties({ D, t, E, rhoSteel, rhoInt, rhoExt, coatings, buoyancy, Ca, g }) -> { As, Ai, Ae, I, EA, EI, Dhydro, massSteel, massCoating, massBuoyancy, massContents, m, ma, wDry, w, displacedArea }` — m², m⁴, N, N·m², kg/m, N/m — apparent weight bookkeeping (Sparks 2007 Ch. 2).
- `effectiveTension({ Tw, pi, pe, Ai, Ae }) -> number` — N — Te = Tw − pi·Ai + pe·Ae (Sparks 2007 Eq. 2.5).
- `trueWallTension({ Te, pi, pe, Ai, Ae }) -> number` — N — inverse of the above.
- `wallStresses({ Te, M, pi, pe, D, t }) -> { Tw, axial, bendingOuter, hoopInner, hoopOuter, radialInner, radialOuter, vonMisesInner, vonMisesOuter, vonMises }` — N, Pa — Lamé thick-wall stresses + von Mises (Timoshenko & Goodier §28).
- `catenaryRiser({ h, w, hangoffAngleDeg | topTension | H }) -> { H, a, topTension, V, suspendedLength, horizontalSpan, hangoffAngleDeg, tdpCurvature, topCurvature }` — N, m, deg, 1/m — inextensible catenary (Irvine §1.2).
- `boundaryLayerTdp({ EI, H, w, xi }) -> { lambda, catenaryCurvature, pointReaction, tdpShear, curvature[] }` — m, 1/m, N — κ = (w/H)(1 − e^(−ξ/λ)), λ = √(EI/H) (Plunkett 1967; Sparks 2007 Ch. 13).
- `boundaryLayerHangoff({ EI, T, deltaAngle, kRot? }) -> { lambda, endMoment, endCurvature, jointRotation, rotationalStiffness }` — m, N·m, 1/m, rad, N·m/rad — M0 = Δθ/(1/k + 1/√(EI·T)) (Sparks 2007 Ch. 13).
- `trenchProfile({ zSeabed, xStart, length, depth }) -> (x) => z` — m — idealised raised-cosine trench (shape is an assumption).
- `solveRiser(cfg) -> rod result + { topTension, bottomTension, maxTension, minTension, hangoffAngleDeg, bottomAngleDeg, maxCurvature, sMaxCurvature, minBendRadius, maxMoment, tdp, sag, hog, config }` — cfg: `sections` (top → seabed), `top`, `bottom`, `topCondition`, `bottomCondition`, `seabed`, `current`, `pressure`, `pointLoads`, `nElem` — N, deg, 1/m, m, N·m — generic riser/cable/umbilical statics on rod.js.
- `solveScr({ pipe, h, hangoffAngleDeg | topTension | H, groundedLength, seabed, topCondition, topSections, nElem, current, pressure }) -> solveRiser result + { catenary }` — steel catenary riser with flex/stress joint and touchdown zone (API RP 2RD §5).
- `solveWaveRiser({ pipe, buoyant, upperLength, buoyancyLength, lowerLength, top, bottom, seabed, bottomCondition, pointLoads, nElem }) -> solveRiser result + { sagClearance, hogElevation, netBuoyancy }` — m, N — lazy / steep / pliant wave (API RP 17B §7).
- `solveBuoyancyLength({ ...solveWaveRiser args, targetSagClearance, minLength, maxLength, tol }) -> { buoyancyLength, lowerLength, result }` — m — Brent search on buoyancy-section length.
- `riserOffsetCases(base, offsets, { heave }) -> [{ offset, topTension, bottomTension, hangoffAngleDeg, maxCurvature, minTension, tdpX, tdpS, tdpExcursion, result }]` — m, N, deg, 1/m — near/far in-plane vessel offsets (API RP 2RD §5.4).
- `quasiStaticHeave(base, { heaveAmp | raoHeave, waveAmp }) -> { heaveAmp, topTensionMax, topTensionMin, tensionAmplitude, tdpRange, maxCurvature, minTension }` — m, N, 1/m — quasi-static hang-off heave.
- `solveTtr({ length, pipe, topTension | tensionFactor, offset, current, kFlexLower, kFlexUpper, tensionerStiffness, nElem, pressure }) -> { topTension, bottomTension, tensionFactor, lowerAngleDeg, upperAngleDeg, stroke, maxMoment, maxCurvature, maxDeflection, lowerJointMoment, upperJointMoment, rod }` — N, deg, m, N·m, 1/m — top-tensioned riser as tensioned beam (API RP 16Q §3; Sparks 2007 Ch. 7–8).
- `drillingRiserOffsets({ ...solveTtr args, offsets, meanLimitDeg = 2, maxLimitDeg = 4 }) -> { rows[{ offset, lowerAngleDeg, upperAngleDeg, topTension, stroke, ucMean, ucMax }], offsetAtMeanLimit, offsetAtMaxLimit }` — m, deg, N — flex-joint angles vs operating limits (API RP 16Q; limits are arguments).
- `driftOff({ mass, force, dragCoefficient, limitOffset, x0, edsTime, reactionTime, nPoints }) -> { terminalVelocity, timeConstant, timeToLimit, redOffset, yellowOffset, t, x, v }` — m/s, s, m — closed-form drift-off x = v∞·τ·ln cosh(t/τ) (first principles).
- `solveFreeStandingRiser({ length, pipe, canNetUplift, jumper { length, w, vesselX, vesselZ }, kFlexBase, current, canDragForce, nElem }) -> { topOffset, topSetDown, baseAngleDeg, topTension, baseTension, jumperH, jumperTowerV, jumperVesselV, jumperVesselTension, maxMoment, iterations, rod }` — m, deg, N, N·m — tensioned beam + jumper catenary (documented simple model).
- `stressAlong(result, { D, t, pIntRef, zIntRef, rhoInt, rhoExt, zSurface, g }) -> { s, vonMises, axial, hoop, max, sMax }` — m, Pa — wallStresses at every node.

## risercheck.js

- `vonMisesApi2rd({ T, M, pi, pe, D, t, sigmaY, Cf = 1.0, Ca = 2/3 }) -> { sigmaR, sigmaTheta, sigmaZ, vonMises, allowable, uc }` — Pa — API RP 2RD (1998) §5.2.3, Cf 1.0 / 1.2 / 1.5 operating / extreme / survival.
- `burstPressureDnv({ D, t, fy, fu }) -> number` — Pa — pb = (2/√3)(2t/(D − t))·min(fy, fu/1.15) (DNV-OS-F201 Eq. 5.15).
- `burstCheckDnv({ pli, pe, D, t, fy, fu, gammaSC, gammaM }) -> { pb, demand, capacity, uc }` — Pa — DNV-OS-F201 Eq. 5.14.
- `burstPressureApi1111({ D, t, S, U, thin }) -> number` — Pa — 0.45(S + U)·ln(D/Di) or 0.90(S + U)·t/(D − t) (API RP 1111).
- `hoopStress({ pi, pe, D, t }) -> { barlow, meanDiameter }` — Pa — thin-wall hoop stress.
- `collapsePressureDnv({ D, t, fy, E, nu, f0, alphaFab }) -> { pc, pel, pp }` — Pa — collapse cubic with ovality (DNV-OS-F201 Eq. 5.18–5.21; DNV-ST-F101 Sec. 13 closed form).
- `collapseCheckDnv({ ...collapse args, pe, pmin, gammaSC, gammaM }) -> { pc, demand, capacity, uc }` — Pa — DNV-OS-F201 Eq. 5.17.
- `collapsePressureApi({ D, t, S, E, nu }) -> { Pc, Py, Pe }` — Pa — Pc = Py·Pe/√(Py² + Pe²) (API RP 1111 / API RP 2RD).
- `propagatingBuckling({ D, t, fy, alphaFab, pe, pmin, gammaC, gammaM, gammaSC }) -> { pprDnv, ppApi, capacity, uc }` — Pa — 35·fy·αfab·(t/D)^2.5 (DNV) and 24·S·(t/D)^2.4 (API RP 1111).
- `combinedLoadingDnv({ D, t, fy, fu, MF, ME, MA, TeF, TeE, TeA, pld, pe, pmin, gammaF, gammaE, gammaA, gammaSC, gammaM, E, nu, f0, alphaFab }) -> { Md, Ted, Mk, Tk, alphaC, pb, pc, mode, uc }` — N·m, N, Pa — LRFD combined loading, internal (Eq. 5.28) and external (Eq. 5.30) overpressure (DNV-OS-F201 / DNV-ST-F201 §5).
- `localBucklingStrain({ strain, D, t, alphaH, alphaGw, pmin, pe, pb, gammaEps }) -> { epsC, capacity, uc }` — εc = 0.78(t/D − 0.01)(1 + 5.75(pmin − pe)/pb)·αh^−1.5·αgw (DNV-ST-F101 Sec. 5).
- `flexJointCheck({ angleDeg, allowableDeg, kRot }) -> { uc, moment }` — N·m — rotation utilisation.
- `stressJointTaper({ deltaAngle, T, D, t, E, sigmaBendAllow, nPoints }) -> { length, curvature, lambdaPipe, rootOD, rootEI, rootMoment, xi, EI, OD }` — m, 1/m, N·m², N·m — constant-curvature taper EI(ξ) = EIp(1 + ξ/λ + ½(ξ/λ)²) (first-principles tensioned-beam statics).
- `connectorUtilisation({ T, M, envelope[{ T, M }] }) -> { uc, capacityT, capacityM }` — N, N·m — ray utilisation against a rated capacity polygon.
- `fatigueHotSpot({ dT, dM, D, t, scf, radius, cycles, snCyclesToFailure }) -> { stressRange, axialRange, bendingRange, damage }` — Pa — Δσ = SCF(ΔT/A + ΔM·r/I), Miner damage with caller S–N (DNV-RP-C203 §2).

## viv.js

- `strouhalNumber(Re, { stHighRe = 0.2 }) -> number` — St = 0.212(1 − 21.2/Re) (50–150), 0.198(1 − 19.7/Re) (250–2×10⁵), blend between, argument above (Roshko 1954; Blevins §3.2).
- `sheddingFrequency({ U, D, St }) -> number` — Hz — fs = St·U/D.
- `reducedVelocity({ U, fn, D }) -> number` — VR = U/(fn·D).
- `stabilityParameter({ me, zeta, rho, D, St }) -> { Ks, SG, delta }` — KS = 4π·me·ζ/(ρD²), SG = 2π·St²·KS (DNV-RP-C205 §9.1).
- `lockIn({ Vr, Ks, inline, crossflow, crossflowPeak, ksInlineMax }) -> { inline, crossflow, crossflowPeak }` — DNV-RP-C205 §9 windows as arguments.
- `inlineOnsetF105({ Ks, gammaOn, gammaK }) -> number` — in-line onset reduced velocity (DNV-RP-F105 §4.3.5).
- `keuleganCarpenter({ Um, T, D }) -> number` — KC = Um·T/D (wave-induced VIV note).
- `responseAmplitude({ SG, gamma, model, coefficients, override }) -> number` — A/D; Griffin 1.29γ/(1 + 0.43·SG)^3.35 or Sarpkaya 0.32γ/√(0.06 + SG²) (Blevins Table 3-2).
- `dragAmplification({ Cd0, aRmsOverD, aOverD, model }) -> { Cd, factor }` — 1 + 1.043(2·Arms/D)^0.65 (Vandiver 1983) or 1 + 2.1·A/D (Blevins 1990).
- `tensionedBeamModes({ L, EI, m, T | T0, T1, nModes, nPoints, nBasis, forceRitz }) -> { s, modes[{ n, fHz, wkbHz, phi, curvature }] }` — m, Hz, 1/m² — closed form (uniform T) or Rayleigh–Ritz (linearly varying T).
- `hangingStringFrequency({ n, L, m, T0, T1, besselJ0, besselY0 }) -> number` — Hz — exact Bessel solution for linearly varying tension.
- `hydrodynamicDamping({ zone, omega, D, V, aOverD, rho, nu, Csw, Cvl, Cvh }) -> number` — N·s/m² — still-water / low- / high-reduced-velocity damping (Venugopal 1996 form; coefficients are arguments).
- `vivModalScreening({ s, modes, U, D, mass, zeta, rho, St, vrWindow, Cl, responseModel, aOverDOverride, suppression, damping }) -> { modes[{ index, fHz, powerInLength, powerInFraction, vrMin, vrMax, me, zetaEff, Ks, SG, gamma, aOverD, aRmsOverD, aOverDEnergy, power, dragFactor, Dref, excited }], dominant, excited }` — m, kg/m, W — modal power-in screening in sheared current (Vandiver 1993; DNV-RP-C205 §9; DNV-RP-F204 §4).
- `vivFatigue({ s, modes, screening, E, stressDiameter, scf, snCyclesToFailure, modesUsed }) -> { s, stressAmplitude, damagePerYear, maxDamagePerYear, sMax, lifeYears }` — Pa, 1/yr, m, yr — modal curvature → stress → damage (DNV-RP-F204 §4).
- `vivFatigueLife({ ...screening and fatigue args, cases[{ probability, U }] }) -> { s, damagePerYear, lifeYears, sMax, cases }` — current-profile scatter weighting (DNV-RP-F204 §4.2).
- `suppressionCoverage({ ...screening args, efficiency, targetAOverD, from, tol }) -> { coveredLength, coverageFraction, zone, aOverD, achievable }` — m — strake / fairing coverage solver.
- `spanVivScreening({ frequenciesHz, U, D, St, Ks, ranges }) -> [{ fHz, Vr, fsOverFn, inline, crossflow, crossflowPeak }]` — free-span / jumper hook with caller frequencies.
- `wakeOscillator({ massRatio, zeta, Ur, St, A = 12, eps = 0.3, CL0 = 0.3, CD = 2.0, CM = 1, cycles, stepsPerCycle }) -> { Ur, aOverD, liftMagnification, frequencyRatio, peakUr, peakAOverD, lockInAmplitude, mu, gamma, M, SG }` — van der Pol wake oscillator with acceleration coupling (Facchinetti, de Langre & Biolley 2004).

## pipelay.js

- `sLay({ pipe, D, hingeHeight, stingerRadius, stingerLength, tension, groundedLength, seabed, rollerStiffness, nElem }) -> { topTension, bottomTension, overbendCurvature, overbendStrain, overbendStrainAnalytic, sagbendCurvature, sagbendMoment, sagbendStrain, sSagbend, liftOffS, departureAngleDeg, bearsOnTip, stingerTipAngleDeg, tdpX, layBack, layBackFromHinge, suspendedLength, maxRollerReaction, rod }` — N, 1/m, N·m, m, deg, N/m — S-lay over a circular stinger on rod.js (DNV-ST-F101 Sec. 10).
- `minimumLayTension({ ...sLay args, sagbendStrainLimit, tensionMin, tensionMax, tol }) -> { tension, result }` — N — root search on sLay.
- `jLay({ pipe, D, depth, towerAngleDeg, seabed, nElem }) -> { topTension, bottomTension, sagbendCurvature, sagbendStrain, sagbendMoment, layBack, suspendedLength, towerMoment, rod }` — N, 1/m, N·m, m.
- `stiffenedCatenary({ EI, H, w, h, D, EA, topFree, nPoints }) -> { lambda0, lambdaTop, suspendedLength, topTension, sagbendCurvature, sSagbend, sagbendMoment, sagbendStrain, sigma, curvature }` — boundary-layer stiffened catenary (Plunkett 1967; Dixon & Rutledge 1968).
- `reelLay({ D, t, reelRadius, alignerRadius, straightenerRadius, fy, E, nu, f0, pe, pc }) -> { reelStrain, alignerStrain, straightenerStrain, yieldCurvature, reelCurvature, reelMoment, residualCurvature, residualRadius, ovalityBrazier, ovalityF101 }` — 1/m, N·m, m — reel bending strain, elastic–perfectly-plastic spring-back, ovalisation (Brazier 1927; DNV-ST-F101 Sec. 5).
- `cableLay({ w, h, H | topTension | departureAngleDeg, chuteRadius, mbr, swpLimit, payoutSpeed, vesselSpeed }) -> { H, topTension, departureAngleDeg, layBack, suspendedLength, tdpBendRadius, sidewallPressure, ucChuteRadius, ucTdpRadius, ucSidewall, slack }` — N, deg, m, N/m — catenary cable lay, SWP = T/R (DNV-RP-0360 §4).
- `cableSlackLay({ w, d, V, h, Cd, rho }) -> { angleDeg, topTension, layBack, suspendedLength }` — straight-line slack lay at the critical angle (Zajac 1957).
- `abandonmentRecovery({ H, wPipe, wWire, depth, headHeights, headWeight }) -> [{ headHeight, headTension, vesselTension, wireAngleDeg, pipeSuspendedLength }]` — m, N, deg — A&R wire tension vs head height.
- `capstan({ T1, mu, theta }) -> number` — N — T2 = T1·e^(μθ).
- `shorePullIn({ tailTension, startFactor, segments }) -> { winchTension, tensions }` — N — friction + capstan accumulation along a pull route.
- `ploughTowForce({ Cw, Wsub, Cs, gammaSub, depth, Cd, speed }) -> { force, friction, staticSoil, rate }` — N — three-term tow-force hook (Cathie & Wintgens 2001 form; coefficients are arguments).
- `layDynamicTension({ staticTension, w, EA, m, suspendedLength, topAngleDeg, heaveAmp, surgeAmp, period }) -> { axialAmp, lower, upper, dafLower, dafUpper, minTensionUpper }` — m, N — quasi-static + inertia lower bound and drag-locked elastic upper bound.
- `limitingSeaState({ staticTension, maxTension, minTension, hsGrid, dafOfHs }) -> { limitingDaf, allowableHs }` — limiting dynamic amplification hook.
- `layCriteriaF101({ grade, criterion, strainLimit, overbendStrain, D, t, fy, sagbendTe, sagbendMoment, pe, stressFactor }) -> { overbendLimit, ucOverbend, sagbendStress, sagbendAllowable, ucSagbend }` — Pa — simplified laying criteria: overbend strain by grade, sagbend 0.87·fy (DNV-ST-F101 Sec. 10).

## lowering.js

- `staticTension({ mass, volume, cableWeight, depths, rho, g }) -> { objectWeightInWater, rows[{ depth, tension }] }` — N, m — DNV-RP-N103 §5.2.1.
- `axialResonance({ mass, addedMass, EA, length, cableMass, kc, theta }) -> { periodExact, periodLumped, omegaExact, omegaLumped, stiffness }` — s, rad/s, N/m — (ωL/c)·tan(ωL/c) = mL/M′ and the one-third-cable-mass lumped form (DNV-RP-N103 §5.3.3).
- `heaveResponse({ mass, addedMass, EA, length, cableMass, kc, heaveAmp, period, Cd, area, rho, linearDamping, staticTension, theta }) -> { motionAmp, motionRatio, dynamicForce, daf, minTension, maxTension, slack, naturalPeriod, damping, iterations }` — m, N, s, N·s/m — SDOF with linearised quadratic drag, slack-wire check (DNV-RP-N103 §5.3).
- `heaveCompensator({ mass, k, c, period, activeResidual }) -> { transmissibility, efficiency, naturalPeriod, dampingRatio }` — spring–damper stage transmissibility (DNV-RP-N103 §5.5).
- `pendulumPeriod({ length, amplitudeDeg, g }) -> number` — s — 2π√(L/g)(1 + θ0²/16).
- `currentOffset({ length, objectWeight, cableWeight, cableDiameter, cableCd, objectDragArea, current, rho, nSteps }) -> { offset, objectDrag, topAngleDeg, bottomAngleDeg, offsetProfile }` — m, N, deg — horizontal offset in current (DNV-RP-N103 §5.2.2).
- `landingVelocity({ loweringSpeed, objectMotionAmp, period }) -> { velocity, motionVelocity }` — m/s.
- `craneTipMotion({ heaveAmp, rollAmpDeg, pitchAmpDeg, b, l, T3, T4, T5 }) -> { eta, velocity, acceleration }` — m, m/s, m/s² — DNV-RP-N103 §4.3.3.
- `waveKinematicsN103({ Hs, Tz, d, amplitudeFactor = 0.9, g }) -> { zetaA, velocity, acceleration, tzInRange }` — m, m/s, m/s² — ζa = 0.9·Hs and exponential depth decay (DNV-RP-N103 §4.3.2–4.3.4).
- `splashZoneN103({ mass, volume, addedMass, slamArea, Cs, dragArea, CD, waterplaneArea, vc, etaCt, vCt, aCt, zetaA, vW, aW, staticMin, staticMax, rho, g }) -> { slamVelocity, Fslam, Frho, FM, FD, Fhyd, staticMin, staticMax, slackUc, slackOk, Ftotal, dafConv }` — m/s, N — Simplified Method forces and slack-sling criterion Fhyd ≤ 0.9·Fstatic-min (DNV-RP-N103 §4.3–4.4).
- `snapLoadN103({ mass, addedMass, K, vff, vr }) -> { Fsnap, vSnap, C }` — N, m/s — Fsnap = vsnap·√(K(M + A33)) (DNV-RP-N103 §4.7.2).
