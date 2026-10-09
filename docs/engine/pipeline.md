# `src/engine/pipeline` — subsea pipeline mechanical design engines

SI units throughout (m, kg, s, N, Pa, J, K); `…C` arguments in °C, `…Deg` in degrees, other angles in radians.
Every function takes one options object unless shown otherwise and returns a plain object with named fields; the
JSDoc of each function lists every argument and returned field with its unit, the reference and the validity range.
`[..]` marks optional arguments.

Sign conventions: effective axial force is tension-positive in `expansion.js` and `span.js`; forces in `buckling.js`
(P, P0) are compressive-positive; beam loads in `crossing.js` are downward-positive, deflections upward-positive.

## pipeprops.js — section, weights, materials, standard sizes

| export | signature → result (units) | reference |
| --- | --- | --- |
| `pipeSection` | `({ D, t })` m → `{ Di, As, Ai, Ae, I, J, Z, rGyr, DtRatio }` (m, m², m⁴, m³) | annulus geometry |
| `pipeWeights` | `({ D, t, [rhoSteel], [coatings: [{ thickness, density, absorption }]], [rhoContents], [rhoWater], [rhoFlood], [marineGrowth], [g] })` → masses kg/m, `weightDry`, `buoyancy`, `submergedEmpty/Flooded/Operating` N/m, `sg*` | Archimedes; layer geometry |
| `steelGrade` | `(grade)` → `{ smys, smts, E, nu, alpha, rho }` Pa, 1/K, kg/m³ | API 5L PSL 2 Table 7; DNV-ST-F101 Table 7-5 |
| `listSteelGrades` | `()` → `string[]` | — |
| `temperatureDerating` | `(tempC, [{ curve }])` → de-rating Pa | DNV-ST-F101 Fig. 5-2 (C-Mn / 13Cr), override by `curve` |
| `fabricationFactor` | `(process)` → α_fab | DNV-ST-F101 Table 5-5 |
| `materialStrengthFactor` | `([supplementaryU])` → α_U | DNV-ST-F101 Table 5-3 |
| `characteristicStrength` | `({ smys, smts, [tempC], [alphaU], [fyTemp], [fuTemp], [curve] })` → `{ fy, fu }` Pa | DNV-ST-F101 Eq. 5.4, 5.5 |
| `npsOuterDiameter` | `(nps)` → D m | ASME B36.10M |
| `standardPipeTable` | `([{ sizes, schedules }])` → `[{ nps, D, walls: [{ schedule, t, Di, massSteel }] }]` | ASME B36.10M (STD, XS, 40, 80, 160) |
| `scheduleWallThickness` | `(nps, schedule)` → t m | ASME B36.10M |

## strength.js — wall thickness and limit states

| export | signature → result (units) | reference |
| --- | --- | --- |
| `dnvResistanceFactors` | `(safetyClass, [limitState])` → `{ gammaM, gammaScPc, gammaSc, gammaEps, alphaMpt, alphaSpt }` | DNV-ST-F101 Tables 5-1, 5-2, 5-8, 5-9 |
| `characteristicWallThickness` | `({ t, [tFab], [tFabFrac], [tCorr], [condition] })` → `{ t1, t2 }` m | DNV-ST-F101 Table 5-6 |
| `burstPressure` | `({ D, t, fy, fu })` → `{ pb, fcb }` Pa | DNV-ST-F101 Eq. 5.8, 5.9 |
| `localIncidentalPressure` | `({ pRef, [gammaInc], [rhoCont], [hRef], [hLocal] })` → Pa | DNV-ST-F101 Eq. 4.1 |
| `millTestPressure` | `({ D, tMin, smys, smts, [k] })` → Pa | DNV-ST-F101 Eq. 7.3 |
| `pressureContainment` | `({ D, t, [tFab], [tFabFrac], [tCorr], fy, fu, pLocalInc, [pe], [safetyClass], [pLocalTest], [pMill], … })` → `{ t1, pb, pbDesign, demand, capacity, uc, governing, testCheck }` | DNV-ST-F101 Eq. 5.6, 5.7 |
| `requiredWallThicknessBurst` | `({ D, fy, fu, pLocalInc, [pe], [safetyClass], [tCorr], [tFab], [tFabFrac] })` → `{ t1, tNominal }` m | Eq. 5.6 / 5.8 inverted |
| `collapsePressure` | `({ D, t, fy, [E], [nu], [alphaFab], [f0] })` → `{ pc, pel, pp, residual }` Pa | DNV-ST-F101 Eq. 5.10–5.12 |
| `collapseCheck` | `({ D, t1, fy, pe, [pMin], [safetyClass], … })` → `{ pc, demand, capacity, uc }` | DNV-ST-F101 Eq. 5.10 |
| `requiredWallThicknessCollapse` | `({ D, fy, pe, [pMin], … })` → t1 m | bisection on `collapseCheck` |
| `propagationPressure` | `({ D, t2, fy, [alphaFab], [pe], [pMin], [safetyClass] })` → `{ ppr, capacity, uc, arrestorsRequired }` | DNV-ST-F101 Eq. 5.15, 5.16 |
| `buckleArrestorCapacity` | `({ D, t2, tBA, LBA, fy, [fyBA], [pe], … })` → `{ ppr, pprBA, pX, capacity, uc }` | DNV-ST-F101 Eq. 5.17, 5.18 (coefficients to be confirmed) |
| `plasticCapacities` | `({ D, t, fy })` → `{ Sp N, Mp N·m }` | DNV-ST-F101 Eq. 5.20, 5.21 |
| `localBucklingLoadControlled` | `({ D, t2, fy, fu, MSd, SSd, [pi], [pe], [pMin], [safetyClass], … })` → `{ mode, uc, alphaC, alphaP, beta, allowableMoment, warnings }` | DNV-ST-F101 Eq. 5.19–5.24, 5.28 |
| `localBucklingDisplacementControlled` | `({ D, t2, fy, fu, epsSd, [pMin], [pe], [alphaH], [alphaGw], [safetyClass] })` → `{ mode, epsC, epsRd, uc }` | DNV-ST-F101 Eq. 5.30–5.32 |
| `equivalentStressCheck` | `({ D, t, fy, eta, [pi], [pe], [N], [M], [T] })` → stresses Pa, `uc` | von Mises; ASD format |
| `b318DesignFactor` | `(locationClass)` → F | ASME B31.8 Table 841.1.6-1, A842.2.2-1 |
| `b318TemperatureFactor` | `(tempC)` → T | ASME B31.8 Table 841.1.8-1 |
| `asmeStressCheck` | `({ D, t, smys, pi, [pe], [F1], [F2], [F3], [restrained], [deltaT], [axialForce], [M], [T], [criterion] })` → `{ SH, SLmax, SLmin, Scombined, ucHoop, ucLong, ucCombined, tRequiredHoop }` | ASME B31.8 §841, §833, §A842; ASME B31.4 §402–403 |
| `lameStresses` | `({ D, t, [pi], [pe], [r] })` → `{ sigmaH, sigmaR, sigmaAxialCapped, … }` Pa | Lamé |
| `hydrotestPressure` | `({ pDesign, [factor], [D], [t], [smys], [pe] })` → `{ pTest, hoopStress, hoopFractionSmys }` | ASME B31.4 §437.4.1; B31.8 §A847.2 |
| `bendingStrainRadius` | `({ D, [R], [epsAllow], [sigmaAllow], [E] })` → `{ strain, minRadiusStrain, minRadiusElastic }` | beam kinematics |
| `brazierOvalisation` | `({ D, t, curvature, [E], [nu] })` → `{ zeta, ovality, moment, brazierMoment, criticalCurvature }` | Brazier (1927) |
| `ovalisationBendingPressure` | `({ D, t, strain, [f0], [pe], [pc], [c1], [c2], [limit] })` → `{ ovality, uc }` | DNV-ST-F101 Eq. 5.34 (coefficients exposed) |
| `elbowFactors` | `({ t, R, r })` → `{ h, k, i }` | ASME B31.4 Table 402.1-1; B31.8 App. E |
| `bendWallThickness` | `({ D, t, R })` → `{ tExtrados, tIntrados, thinningFraction, lorenzIntrados, lorenzExtrados }` | kinematics; ASME B31.3 §304.2.1 |

## stability.js — on-bottom stability (pipes, cables, umbilicals)

Active load model: Morison with worst-phase search (`method: 'morison'`, default) or caller-supplied peak coefficients
(`method: 'coefficients'`). The DNV-RP-F109 peak-coefficient tables are not built in; each result names the method used.

| export | signature → result (units) | reference |
| --- | --- | --- |
| `waveNumber` | `(omega, depth, [g])` → k rad/m | linear dispersion |
| `jonswapGamma` | `(Hs, Tp)` → γ | DNV-RP-F109 §3.4.3 |
| `jonswapSpectrum` | `(omega, { Hs, Tp, gamma }, [g])` → m²·s/rad | DNV-RP-F109 Eq. 3.5–3.7 |
| `nearBedSpectralVelocity` | `({ Hs, Tp, depth, [gamma], [nOmega] })` → `{ Us m/s, Tu s, M0, M2, gamma, Tn }` | DNV-RP-F109 Eq. 3.8–3.12 |
| `designOscillation` | `({ Us, Tu, [durationS], [RD], [gamma], [Tn], [kT] })` → `{ Ustar m/s, Tstar s, tau, kT }` | DNV-RP-F109 Eq. 3.15–3.17 |
| `currentOverPipe` | `({ Vref, zRef, z0, D, [angleDeg] })` → m/s | DNV-RP-F109 Eq. 3.3 |
| `loadReductionFactors` | `({ [zp], D, [trenchAngleDeg], [permeable] })` → `{ rPenY, rPenZ, rTrY, rTrZ, rPermZ, rTotY, rTotZ }` | DNV-RP-F109 Eq. 3.19–3.24 |
| `passiveResistance` | `({ soil, FC, D, zp, [gammaSoilSub], [su], [gammaSoilDry] })` → `{ FR N/m, kappa, Gc }` | DNV-RP-F109 §3.4.6 |
| `initialPenetration` | `({ soil, FC, D, … })` → `{ zpi m, kappa, Gc }` | DNV-RP-F109 §3.4.6 |
| `verticalStability` | `({ ws, b, [gammaW] })` → `{ uc, sg }` | DNV-RP-F109 Eq. 3.1 |
| `morisonCycle` | `({ D, Ustar, Vstar, Tstar, [rhoW], [CD], [CM], [CL], [nPhase] })` → `{ phaseDeg, FY, FZ N/m, FYmax, FZmax, K, M }` | Morison; DNV RP E305 §5.3 |
| `absoluteLateralStability` | `({ ws, D, mu, Ustar, Vstar, [Tstar], gammaSC, [method], [CY], [CZ], [Fw], [zp], [soil], … })` → `{ method, ucLateral, ucLift, uc, FY, FZ, FR, phaseDeg, K, M }` | DNV-RP-F109 Eq. 3.38, 3.39 |
| `requiredSubmergedWeight` | `({ …absoluteLateralStability args without ws, [targetUc] })` → `{ ws N/m, check }` | bisection |
| `requiredConcreteThickness` | `({ pipe, [concreteDensity], [concreteAbsorption], stability, [gammaW], [targetUc], [tMax], [tStep] })` → `{ thickness m, ws, Dtotal, sg, check }` | DNV-RP-F109 Eq. 3.1, 3.38, 3.39 |

## pipesoil.js — pipe–soil interaction

| export | signature → result (units) | reference |
| --- | --- | --- |
| `contactGeometry` | `({ D, z })` → `{ B m, beta rad, Abm m², zeta, arc m }` | geometry; wedging per DNV-RP-F114 §4.3 |
| `verticalBearingUndrained` | `({ D, z, su0, [suGrad], [gammaSub] })` → `{ Qv N/m, suInvert, bearingFactor, buoyancyTerm }` | DNV-RP-F114 §4.2 |
| `verticalBearingDrained` | `({ D, z, phiDeg, gammaSub, [nGamma], [surchargeDepth] })` → `{ Qv N/m, Nq, Ngamma, dq, B, z0 }` | DNV-RP-F114 §4.2; Brinch Hansen / Vesić |
| `embedment` | `({ D, W, model, [kLay], [kDyn], … })` → `{ z, zStatic m, zOverD, V N/m, geometry }` | DNV-RP-F114 §4.2 |
| `axialResistance` | `({ model, D, z, [W], [mu], [alpha], [su], [residualRatio] })` → `{ breakout, residual N/m, zeta, equivalentFriction }` | DNV-RP-F114 §4.3 |
| `lateralResistance` | `({ [model], soil, D, z, W, mu, [muResidual], [su], [gammaSub], [phiDeg], [bermFactor], [bermFactorResidual] })` → `{ breakout, residual, passive N/m, frictionBreakout, frictionResidual }` | Coulomb + Rankine; DNV-RP-F109 §3.4.6 (`'verley'`) |
| `mobilisationCurve` | `({ breakout, [residual], xBreakout, [xResidual] })` → `{ x, f, force(u), stiffness }` | tri-linear idealisation |
| `frictionRange` | `({ LE, HE, [pLow], [pHigh] })` → `{ mu, sigma, BE, mean, std, cov, quantile(p) }` | lognormal quantile matching |
| `soilStiffness` | `({ D, massRatio, [soilType], [CV], [CL], [nu], [KVS] })` → `{ KV, KL, KVS N/m² }` | DNV-RP-F105 (2006) §7.4, Tables 7-5/7-6 |
| `upliftResistanceSand` | `({ D, H, gammaSub, f })` → N/m | Schaminée et al. (1990) |
| `upliftResistanceClay` | `({ D, H, su, [gammaSub], [Nc] })` → `{ R, shallow, deep N/m, governing }` | vertical slip / local failure |
| `requiredCover` | `({ soil, R, D, [gammaSub], [f], [su], [Nc] })` → H m | inversion of the uplift models |
| `coverQuantities` | `({ D, H, crestWidth, sideSlope, length, rhoSolid, [porosity] })` → `{ area m², volume m³, mass kg, baseWidth m, gammaSub N/m³, prismWeight N/m }` | geometry |

## expansion.js — end expansion and walking

| export | signature → result (units) | reference |
| --- | --- | --- |
| `restrainedAxialForce` | `({ [H], deltaPi, Ai, As, deltaT, [E], [alpha], [nu] })` → `{ S0 N, pressureTerm, thermalTerm, freeStrain }` | DNV-RP-F110 Eq. 6 |
| `temperatureProfile` | `({ x, TinC, TambC, [lambda] \| massFlow, cp, U, Dref })` → °C | steady energy balance |
| `pressureProfile` | `({ x, pIn, [frictionGradient], [rho], [z] })` → Pa | steady momentum balance |
| `endExpansion` | `({ L, As, Ai, friction, deltaT, [deltaPi], [H], [Rhot], [Rcold], [kHot], [kCold], [n] })` → `{ x, S0, S N, regime, anchorHot, anchorCold m, expansionHot, expansionCold m, Rhot, Rcold N, maxCompression }` | friction build-up; DNV-RP-F110 §6 |
| `endExpansionUniform` | `({ P, f, EA, [L], [R] })` → `{ anchorLength m, expansion m, regime }` | closed form |
| `spoolGuidedCantilever` | `({ delta, D, I, [Lleg], [sigmaAllow], [E], [sif] })` → `{ force N, moment N·m, stress Pa, stiffness N/m, requiredLeg m }` | guided cantilever method |
| `walkingSlope` | `({ deltaP, W, L, mu, slopeRad, EA })` → `{ walkPerCycle m, anchorSeparation m, fullyMobilised }` | Carr, Sinclair & Bruton (2006), OTC 17945 |
| `walkingEndTension` | `({ deltaP, SR, W, L, mu, EA })` → same | Carr, Sinclair & Bruton (2006) |
| `walkingArrestForce` | `({ W, L, mu, [slopeRad], [SR] })` → `{ holdBackTension, rigidAnchorTension, rigidAnchorCompression }` N | equilibrium of the mobilised friction profiles |
| `thermalTransientSteps` | `({ L, nElem, alphaDeltaT, [lambda], [pressureStrain], [nHeat], [nCool] })` → `Float64Array[]` free strain per element | idealised start-up transient |
| `uniformCycleSteps` | `({ nElem, freeStrain, [nUp], [nDown] })` → `Float64Array[]` | — |
| `axialWalkingModel` | `({ L, EA, W, mu, steps, [cycles], [slopeRad], [tensionAtEnd], [tensionAtStart], [uMob], [anchor], [tol] })` → `{ x, u m, force N, walkPerCycle m, walkRate m/cycle, hotForce, coldForce, anchorForceMin/Max N }` | 1-D elastic bar with Coulomb friction (numerical) |

## buckling.js — global buckling

| export | signature → result (units) | reference |
| --- | --- | --- |
| `HOBBS_CONSTANTS` | `{ 1, 2, 3, 4, inf: { k1…k5 } }` | Hobbs (1984) Table 1 |
| `hobbsLateral` | `({ mode, L, EI, EA, w, muLateral, [muAxial], [D] })` → `{ P, P0 N, amplitude m, moment N·m, strain, slipLength m }` | Hobbs (1984) |
| `hobbsVertical` | `({ L, EI, EA, w, muAxial, [D] })` → `{ P, P0, amplitude, moment, slope, strain }` | Hobbs (1984) |
| `hobbsMinimum` | `({ mode \| 'vertical', EI, EA, w, …, [alpha], [deltaPi], [Ai], [H] })` → `{ L, P0min, P, amplitude, moment, strain, deltaT K }` | Hobbs (1984); numerical minimum |
| `hobbsCurve` | `({ mode, Lmin, Lmax, [n], … })` → `{ L, P0, P, amplitude }` | Hobbs (1984) |
| `infiniteModeForce` | `({ EA, EI, fL, [S0] })` → `{ Sinf N, L m, constant, susceptible, uc }` | Hobbs mode ∞ minimised; DNV-RP-F110 screening form |
| `palmerPhiW` | `(phiL)` → Φ_w | Palmer et al. (1990), OTC 6335 |
| `propImperfectionLength` | `({ EI, delta, w })` → L m | beam statics |
| `upheavalDownload` | `({ P, EI, delta, [L], [wInstall], wPipe, [safetyFactor], [cover] })` → `{ phiL, phiW, L, wRequired, coverResistance N/m, coverHeight m }` | Palmer et al. (1990); Schaminée et al. (1990) |
| `buckleSpacing` | `({ P0, Ppost, f, EA, [spacing], [PcritUnplanned], [PcritTrigger] })` → `{ slipLength m, feedInIsolated, feedIn m, forceAtAnchor N, maxSpacing, minSpacing m }` | axial friction equilibrium |

## span.js — free spans

| export | signature → result (units) | reference |
| --- | --- | --- |
| `boundaryCoefficients` | `(type, [leffOverL])` → `{ C1…C6, C4mid, C5mid }` | DNV-RP-F105 (2006) Table 6-1 |
| `concreteStiffnessFactor` | `({ EIconc, EIsteel, [kc] })` → CSF | DNV-RP-F105 Eq. 6.5 |
| `effectiveSpanLength` | `({ L, K, EI, [CSF] })` → `{ ratio, Leff m, beta }` | DNV-RP-F105 Eq. 6.13 |
| `addedMassCoefficient` | `(gapRatio)` → C_a | DNV-RP-F105 Eq. 6.21 |
| `spanResponse` | `({ L, D, Ds, t, EI, me, [E], [Seff], [q], [CSF], [boundary], [KV], [KL], [KVS] })` → `{ LeffCF, LeffIL, LeffStatic m, PcrCF, PcrIL N, deflection m, momentStatic N·m, stressStatic Pa, f1IL, f1CF Hz, AIL, ACF Pa, warnings }` | DNV-RP-F105 §6.7 |
| `stabilityParameter` | `({ me, zetaT, D, [rhoW], [gammaK] })` → `{ Ks, Ksd }` | DNV-RP-F105 Eq. 4.6 |
| `inlineTurbulenceReduction` | `({ Ic, [thetaRel] })` → `{ R1, R2 }` | DNV-RP-F105 Eq. 4.9–4.10 (to be confirmed) |
| `inlineResponse` | `({ VR, Ksd, [gammaOn], [R1], [R2], [alpha] })` → `{ AyOverD, psiAlpha, VRonset, VR1, VR2, VRend, AY1, AY2 }` | DNV-RP-F105 §4.3 |
| `crossflowResponse` | `({ VR, [Ksd], [alpha], [KC], [fRatio], [gapRatio], [trenchRatio], [massRatio], [gammaOn] })` → `{ AzOverD, Rk, VRonset, VR1, VR2, VRend, AZ1 }` | DNV-RP-F105 §4.4 |
| `vivScreening` | `({ f1IL, f1CF, Uc100, [Uw1], D, L, VRonsetIL, VRonsetCF, [gammaIL], [gammaCF] })` → `{ ilOk, cfOk, ucIL, ucCF, alphaBar }` | DNV-RP-F105 §2.3 |
| `vivFatigue` | `({ f1IL, f1CF, AIL, ACF, D, Ksd, current: { scale, shape, loc }, [seaStates], sn, … })` → `{ damagePerYearIL, damagePerYearCF, lifeILYears, lifeCFYears, lifeYears, maxStressRangeIL/CF Pa }` | DNV-RP-F105 §2.4, §4.3–4.4; Miner |
| `spanUls` | `({ response, D, Uextreme, [alpha], [Ksd], [KC], [sigmaDrag], [gammaS], [sigmaAllow] })` → `{ sigmaStatic, sigmaEnv, sigmaTotal Pa, uc }` | DNV-RP-F105 §2.5 |
| `maxSpanLength` | `({ criterion, span, [sigmaAllow], [U], [VRonset], [gammaF], [Lmin], [Lmax] })` → `{ L m, response }` | bisection on `spanResponse` |

## integrity.js — defect assessment

| export | signature → result (units) | reference |
| --- | --- | --- |
| `b31gOriginal` | `({ D, t, d, L, smys, [safetyFactor] })` → `{ z, M, flowStress, failureStress, failurePressure, safePressure Pa }` | ASME B31G (original) |
| `b31gModified` | same → same | ASME B31G Modified (0.85 dL); Kiefner & Vieth (1989) |
| `effectiveArea` | `({ D, t, x, depth, smys, [safetyFactor] })` → `{ failureStress, failurePressure, safePressure, effectiveLength, effectiveArea, iStart, iEnd }` | ASME B31G Level 2 (RSTRENG) |
| `dnvF101Single` | `({ D, t, d, L, fu, [F2], [F1], [gammaM], [gammaD], [epsD], [stdDepth] })` → `{ Q, failurePressure, safePressure, allowablePartA Pa }` | DNV-RP-F101 Part B; Part A format with caller factors |
| `remainingLife` | `({ assess, defect, maop, depthRate, [lengthRate], [depthLimit], [horizon], [reinspectFraction], [interval] })` → `{ remainingLife yr, governing, reinspectionInterval yr, safePressureNow, deratedPressure Pa }` | ASME B31G §1.6 loop |
| `dentStrain` | `({ t, R0, R1, R2, d, L, [limit] })` → `{ eps1, eps2, eps3, epsInside, epsOutside, epsMax, uc }` | ASME B31.8 App. R |
| `dentDepthCheck` | `({ d, D, [onWeld], [withStressRaiser] })` → `{ depthRatio, limit, acceptable, note }` | ASME B31.8 §851.4.1 |
| `wallLossExtreme` | `({ blockMaxima, returnPeriod, [t] })` → `{ loc, scale, predictedMax m, remainingWall m, exceedanceOfWall }` | Gumbel (1958) |
| `crackAssessment` | `({ flaw, assess })` → `{ acceptable, uc, detail }` | dispatcher to the fracture engine |
| `fitnessForService` | `({ D, t, d, L, smys, maop, [smts], [profile], [safetyFactor], [F2] })` → `{ levels, acceptable, governingLevel, depthRatio, action }` | ASME B31G level hierarchy |

## crossing.js — crossings, supports, impact

| export | signature → result (units) | reference |
| --- | --- | --- |
| `beamOnSupports` | `({ L, EI, [w], [supports: [{ x, height, stiffness, unilateral }]], [seabed: { stiffness, elevation }], [pointLoads], [ends], [nElem], [D], [E] })` → `{ x, deflection m, rotation, moment N·m, shear N, reactions: [{ x, R N, active }], seabedReaction N/m, spans, maxMoment, maxStress Pa, totalReaction, totalLoad }` | beam FEM with unilateral contact |
| `crossingSeparation` | `({ axisElevation, Dcross, topOfExisting, [minSeparation], [settlement] })` → `{ separation m, ok, margin }` | DNV-ST-F101 crossing separation (0.3 m default) |
| `crossingQuantities` | `({ [mattresses], [rockBerms], [rhoWater] })` → mattress count / volume m³ / mass kg / submerged weight N, rock volume m³ / mass kg | geometry |
| `droppedObjectImpact` | `({ mass, volume, CD, area, [Ca], [rhoW] })` → `{ terminalVelocity m/s, kineticEnergy, impactEnergy J, addedMass kg }` | DNV-RP-F107 §5 |
| `impactEnergy` | `({ mass, [addedMass], velocity })` → J | mechanics |
| `steelPipeDentEnergy` | `({ D, t, fy, [dentRatio], [energy] })` → `{ energy J, dentRatio, mp N }` | DNV-RP-F107 §4.3 |
| `concreteCoatingEnergy` | `({ [mode], Y, b, [h], [D], x0 })` → J | DNV-RP-F107 §4.4 (to be confirmed) |
| `impactProtection` | `({ impactEnergy, [absorbed], D, t, fy, [dentLimit] })` → `{ energyToSteel J, dentRatio, ok, additionalProtection J }` | energy balance |
| `trawlPullover` | `({ [Fp] \| CF, V, mass, kWarp, L, EI, [q], [ends], [D], [E] })` → `{ Fp N, momentMax N·m, deflection m, stress Pa, soilHoldsLoad }` | DNV-RP-F111 force format (coefficients supplied); beam statics |
| `rockSizeShields` | `({ thetaCr, [tau] \| U, f, [rhoRock], [rhoW] })` → `{ D50 m, tau Pa }` | Shields (1936) |

## Not built in (supply as arguments, or use the stated alternative)

* DNV-RP-F109 peak load coefficient tables (C_Y*, C_Z*), γ_SC tables by region, and the generalised (displacement) method.
* DNV-RP-F101 Part A partial-factor tables, interacting-defect rules and combined-loading equations.
* DNV-RP-F111 pull-over coefficient and warp-stiffness expressions; DNV-RP-F107 coefficient tables.
* Duplex / super-duplex temperature de-rating curves; girth-weld factor curve α_gw of DNV-ST-F101.
* Analytical thermal-transient walking expression (the numerical walking model covers it).
* Finite-element lateral buckling with imperfections, multi-span and force-model wave loading on spans,
  EPRG dent–gouge model.
