# Pile-foundation geotechnical engines (`src/engine/geotech`)

SI units throughout (m, N, Pa, kg, s, rad; `…Deg` arguments in degrees). Depth `z` is measured downwards from the
original mudline. A *curve object* is `{ type, pu, p(x), k(x), meta }` with `k = dp/dx`; it is used for p–y
(N/m vs m), t–z (Pa vs m), Q–z (N vs m) and user-supplied base / moment springs.

Lateral sign convention (`pile.js`, `monopile.js`, `pilegroup.js`): `y` positive with the head load `H`;
`theta = −dy/dz`, positive when the head tilts with the load; positive head moment `M = H·e` for a load acting
above the head. `[H; M] = [[KHH, KHM], [KHM, KMM]]·[y; theta]` with `KHM < 0`.

## soilprofile.js

| Export | Signature → result (units) | Reference |
| --- | --- | --- |
| `validateLayers` | `(layers) → { valid, errors[] }`; layer fields `type`, `top` m, `thickness` m, `gammaSub` N/m³, `suTop` Pa, `suGrad` Pa/m, `phiDeg`, `Dr` 0–1, `eps50`, `OCR`, `G0` Pa or `Vs` m/s, `k` N/m³, `PI` %, `St` | data model |
| `createProfile` | `(layers, { gammaW = 10052 N/m³ }) → { layers, depth m, gammaW, scour }` | data model |
| `layerIndexAt` | `(profile, z m) → index` (interface belongs to the lower layer) | — |
| `verticalEffectiveStress` | `(profile, z m) → σ′v Pa`, including scour overburden reduction | ∫γ′dz; API RP 2A-WSD C6.3 scour commentary / API RP 2GEO / ISO 19901-4 |
| `propertiesAt` | `(profile, z m) → { type, sigmaV Pa, su Pa, phiDeg, Dr, eps50, OCR, G0 Pa, k N/m³, zBelowSurface m, layer, … }` | G0 = ρVs² |
| `applyScour` | `(profile, { general m, local m = 1.5D, D m, reductionDepth m = 6D }) → profile` | API RP 2A-WSD C6.3 (1.5 D local scour, 6 D overburden reduction depth) |
| `scourBaseDepth` | `(profile) → m` | — |
| `sampleProfile` | `(profile, { dz = 0.5 m }) → { z m, sigmaV Pa, su Pa, phiDeg, type }` | — |
| `EXAMPLE_PROFILES`, `exampleProfile` | `(name) → profile`; `softClayNC`, `stiffClayOC`, `mediumDenseSand`, `denseSand`, `layeredSandClay` | generic demonstration values, not site data |

## py.js

| Export | Signature → result (units) | Reference |
| --- | --- | --- |
| `matlockPu` | `({ su Pa, sigmaV Pa, z m, D m, J = 0.5 }) → { pu N/m, Np, shallow, zOverZR }` | Matlock (1970); API RP 2GEO §8.5 / RP 2A-WSD §6.8.2 |
| `matlockTransitionDepth` | `({ su Pa, gammaSub N/m³, D m, J }) → XR m` | API RP 2A-WSD §6.8.2 |
| `pySoftClay` | `({ su, sigmaV, z, D, eps50, J, loading 'static'｜'cyclic', zOverZR, yLinearRatio = 0.1 }) → curve N/m` | Matlock (1970); API RP 2A-WSD §6.8.3 |
| `pyStiffClayNoFreeWater` | `({ su, sigmaV, z, D, eps50, J, loading, nCycles, yLinearRatio }) → curve N/m` | Reese & Welch (1975) |
| `sandCoefficients` | `(phiDeg 20–45) → { C1, C2, C3 }` | Reese, Cox & Koop (1974) wedge/flow expressions behind the API figure |
| `sandInitialModulus` | `(phiDeg 29–45) → k N/m³` = (0.008085 φ^2.45 − 26.09)·10⁶ | closed-form approximation of the API k–φ′ figure (below water table) |
| `sandPu`, `sandTransitionDepth` | `({ phiDeg, sigmaV Pa, z m, D m }) → { pu, pus, pud N/m, shallow }`; `(phiDeg, D) → m` | API RP 2A-WSD §6.8.6 |
| `pySand` | `({ phiDeg, sigmaV, z, D, k N/m³, loading }) → curve N/m` (A = max(3 − 0.8 z/D, 0.9) static, 0.9 cyclic) | O'Neill & Murchison (1983); API RP 2A-WSD §6.8.7 |
| `pyWeakRock` | `({ qur Pa, Eir Pa, z m, D m, alphaR = 1, krm }) → curve N/m` | Reese (1997) |
| `pyJeanjean` | `({ su Pa, su0 Pa, su1 Pa/m, Gmax Pa, z, D, yLinear m }) → curve N/m` | Jeanjean (2009) OTC 20158 |
| `pyUser`, `pyLinear`, `zeroCurve` | `({ y[], p[] })`, `(stiffness)`, `()` → curve | generic |
| `scalePy` | `(curve, { pMult = 1, yMult = 1 }) → curve` (p-multiplier, y-multiplier, cyclic degradation) | Brown et al. (1988); Focht & Koch (1973) |
| `samplePy`, `initialStiffness` | `(curve, { yMax m, n }) → { y, p }`; `(curve, yRef m) → N/m²` | — |
| `pyCurveAt` | `(profile, z m, { D m, loading, J, stiffClaySu = 96e3 Pa, stiffClayCyclicFactor = 1, pMult, yMult }) → curve N/m` | API RP 2GEO §8.5 |

## tz.js

| Export | Signature → result (units) | Reference |
| --- | --- | --- |
| `API_TZ_CLAY`, `API_QZ` | tabulated z/D and mobilised fractions | API RP 2GEO §8.4 / RP 2A-WSD §6.7.2–6.7.3 |
| `tzClay` | `({ tmax Pa, D m, residualRatio = 0.9 }) → curve Pa` | API RP 2A-WSD §6.7.2 |
| `tzSand` | `({ tmax Pa, zPeak = 0.0025 m }) → curve Pa` | API RP 2A-WSD §6.7.2 |
| `qz` | `({ Qp N, D m }) → curve N` (compression only) | API RP 2A-WSD §6.7.3 |
| `hyperbolicCurve` | `({ pmax, ki, Rf = 1, odd = true }) → curve` | Kondner (1963) hyperbola |
| `scaleCurve`, `sampleCurve` | `(curve, factor)`; `(curve, { zMax m, n, zMin })` | — |

## axialcapacity.js

| Export | Signature → result (units) | Reference |
| --- | --- | --- |
| `API_SAND_TABLE` | β, fLim Pa, Nq, qLim Pa by soil description | API RP 2GEO Table 1 |
| `densityDescription`, `apiSandParameters` | `(Dr) → string`; `(layer) → { beta, fLim Pa, Nq, qLim Pa, category }` | API RP 2GEO Table 1 |
| `clayAlpha` | `(psi = su/σ′v) → α` | API RP 2A-WSD §6.4.2 |
| `unitResistanceAt` | `(profile, z m, { betaFactor = 1 }) → { f Pa, q Pa, type, sigmaV Pa }` | API RP 2A-WSD §6.4.2–6.4.3 |
| `axialCapacity` | `({ profile, L m, D m, t m, closedEnded, plugMode, internalFrictionFactor, tensionFactor, setupFactor, betaFactor, dz }) → { Qfo, Qfi, QpPlugged, QpUnplugged, Qp, Qc, Qt N, plugged, qTip Pa, z, f, QfoCum }` | API RP 2A-WSD §6.4–6.5 |
| `capacityVsPenetration` | `({ …, Lmin, Lmax, step m }) → { L m, Qc, Qt, Qfo, Qp N, plugged }` | as above |
| `requiredPenetration` | `({ …, load N, mode, safetyFactor, resistanceFactor, step, tol }) → { L m, capacity N, designResistance N, found }` | API RP 2A-WSD §6.3.4 (1.5 / 2.0 as arguments) |
| `cptShaftFriction`, `CPT_SHAFT_COEFFICIENTS` | `({ qc Pa, sigmaV Pa, z, L, D, t m, deltaCvDeg, coeff, pa = 100e3 }) → f Pa`; sets for simplified ICP-05 and offshore UWA-05, compression / tension | API RP 2GEO Annex C unified CPT expression |
| `cptEndBearingUWA` | `({ qcAvg Pa, D, t }) → { q Pa, Qp N, Ar }` | offshore UWA-05 (Lehane et al. 2005) |
| `groupBlockCapacityClay` | `({ B, W, L m, suShaft, suBase Pa, Nc = 9, nPiles, Qsingle N }) → { Qblock, Qsum, Qgroup N, efficiency, governs }` | Terzaghi & Peck block failure |

## pile.js

| Export | Signature → result (units) | Reference |
| --- | --- | --- |
| `solveLateralPile` | `({ L, stickUp, nElem, D, t ｜ sections ｜ EI, E, timoshenko, shearCoeff, axialLoad N, profile, loading, pyOptions, pMult, yMult ｜ pyAt(z, D), base { shear, moment }, mThetaAt(z), head { H N, M N·m, y m, theta rad, kTrans N/m, kRot N·m/rad, fixedRotation }, tol, maxIter, throwOnFail }) → { z m, y m, theta rad, moment N·m, shear N, p N/m, pu N/m, mobilisation, sigmaBending Pa, head, toe, maxMoment, depthMaxMoment, depthFirstZeroDeflection, convergence }` | beam on nonlinear Winkler foundation, Hermite elements + P–Δ; Reese & Van Impe (2011) |
| `pileHeadCurves` | `(o, { factors }) → { factor, H N, M N·m, y m, theta rad, maxMoment }` | as above |
| `pileHeadStiffness` | `(o) → { tangent, secant: { KHH N/m, KHM N, KMM N·m/rad, fHH, fHM, fMM, asymmetry }, working }` | static condensation |
| `equivalentCantilever` | `({ KHH, KHM, KMM, EI N·m² }) → { lateral, rotation, coupling, fixedHead } m` | beam theory |
| `ultimateLateralCapacity` | `(o, { eccentricity m, yLimit m = 0.1D, nSteps, tol }) → { Hult N, criterion, yAtHult m, maxMoment N·m, curve }` | load stepping |
| `criticalPileLength` | `(o, { Lmin, Lmax m, nLengths, tol }) → { Lcrit m, L, y, theta, toeY }` | definition by converged head response |
| `steelUtilisation` | `({ N N, M N·m, D, t m, fy Pa, usageFactor }) → { sigmaAxial, sigmaBending, sigmaMax Pa, utilisation }` | elastic N/A + M·D/2I |
| `solveAxialPile` | `({ profile, L, D, t, E ｜ EA, nElem, mode, headDisplacement m ｜ headLoad N, residualRatio, zPeakSand, tensionFactor, setupFactor, plugMode, tzAt(z), tipCurve }) → { z, w m, force N, tMob Pa, headLoad N, headDisplacement m, tipLoad N, QsUltimate, QpUltimate N, convergence }` | t–z / Q–z load transfer, API RP 2A-WSD §6.6–6.7 |
| `axialLoadSettlement` | `(o, { wMax m = 0.15D, steps }) → { w m, Q N, tipLoad N, capacity N, wAtCapacity m, headStiffness N/m, forceAtCapacity N }` | as above |
| `hetenyiSemiInfinite` | `({ H, M, k N/m², EI, z[] }) → { beta 1/m, y0 m, theta0 rad, KHH, KHM, KMM, y, theta, moment }` | Hetényi (1946) |
| `bromsClay`, `bromsSand` | `({ su Pa ｜ gammaSub, phiDeg, D, L, e m, My N·m }) → { HuShort, HuLong, Hu N, mode, … }` | Broms (1964a, b) |
| `randolphLateral` | `({ H, M, EI, D, G0 Pa, Ggrad Pa/m, nu }) → { y0 m, theta0 rad, lc m, Gc Pa, rhoc, KHH, KHM, KMM }` | Randolph (1981) |
| `randolphWrothAxial` | `({ L, D, EA N, GL Pa, rho, xi, nu }) → { stiffness N/m, zeta, muL }` | Randolph & Wroth (1978) |

## pilegroup.js (global z up; loads at the cap reference point)

| Export | Signature → result (units) | Reference |
| --- | --- | --- |
| `pileLocalStiffness` | `({ KHH, KHM, KMM, Kaxial N/m, Ktorsion N·m/rad, head 'fixed'｜'pinned' }) → 6×6` | stiffness method |
| `analysePileGroup` | `({ piles: [{ x, y, z m, batterDeg, azimuthDeg, stiffness ｜ K }], load: { Fx, Fy, Fz N, Mx, My, Mz N·m } }) → { displacement { ux, uy, uz m, rx, ry, rz rad }, K, piles: [{ axial N (compression +), shear1, shear2, shear N, moment1, moment2, moment N·m, torsion N·m, local }], equilibriumError }` | rigid-cap stiffness method, Poulos & Davis (1980) |
| `P_MULTIPLIER_TABLE_AASHTO`, `pMultiplier` | `({ spacingOverD, row, table }) → fm` | AASHTO LRFD Table 10.7.2.4-1 (default table; any table may be passed) |
| `axialInteractionFactor` | `({ s, D, L m, nu, rho }) → α` | Randolph & Wroth (1979) |
| `lateralInteractionFactor` | `({ s, D m, EI, Gc Pa, rhoc, psiDeg, coefficient = 0.6 }) → α` | Randolph (1981) |
| `groupStiffnessFromInteraction` | `({ piles: [{ x, y }], K1 N/m, alphaFn(s, psiDeg) }) → { Kgroup N/m, efficiency, loadShare, settlementRatio }` | Poulos (1968) superposition |
| `converseLabarre` | `({ m, n, D, s m }) → Eg` | Converse–Labarre (Bowles) |
| `groupTorsionalCapacity` | `({ piles: [{ x, y m, Hu N, Tu N·m }], centre }) → { Tult N·m, centre }` | statics |
| `pileUtilisationTable` | `({ piles, capacity { compression, tension, lateral N, moment N·m } }) → { rows, maxUtilisation, governing }` | demand / capacity |

## monopile.js

| Export | Signature → result (units) | Reference |
| --- | --- | --- |
| `rigidPileStiffness` | `({ L m, kRef N/m², n, zRef m }) → { KLL N/m, KLR N, KRR N·m/rad }` | integration of Winkler springs; Poulos & Davis (1980) rigid-pile forms |
| `flexiblePileStiffness` | `({ EI, modulus 'constant'｜'linear', k N/m², nh N/m³ }) → { KLL, KLR, KRR }` | Hetényi (1946); Poulos & Davis (1980) after Barber (1953) |
| `flexiblePileStiffnessContinuum` | `({ D m, EI, Es Pa, profile 'homogeneous'｜'linear'｜'parabolic' }) → { KLL, KLR, KRR, Ep Pa }` | EN 1998-5 Annex C Table C.1 (Gazetas) |
| `monopileStiffness` | `(solveLateralPile options) → { tangent, secant: { KLL, KLR, KRR }, working }` | numerical (pile.js) |
| `towerNaturalFrequency` | `({ segments: [{ L, Dbot, Dtop, t m, E Pa, rho kg/m³, addedMass kg/m ｜ EI, massPerLength }], topMass kg, topInertia kg·m², KLL, KLR, KRR, nModes, elementsPerSegment }) → { frequencies Hz, omega rad/s, modeShapes, x m }` | beam FE eigenproblem on a three-spring foundation |
| `mudlineRotationCheck` | `({ thetaTotal rad, thetaAccumulated rad, limitTotalDeg = 0.5, limitAccumulatedDeg = 0.25 }) → { totalDeg, accumulatedDeg, ucTotal, ucAccumulated, ok }` | configurable SLS limits (DNV-ST-0126 practice) |
| `cyclicRotationAccumulation` | `({ thetaStatic rad, N, Tb, Tc, exponent = 0.31 }) → { deltaTheta rad, ratio }` | LeBlanc, Houlsby & Byrne (2010) |
| `embeddedLengthCriteria` | `(o, { Lmin, Lmax m, nLengths, tol }) → { L, headY, headTheta, toeY, toeTheta, zeroToeKick m, verticalTangent m, critical m }` | zero toe-kick / vertical tangent / critical length |
| `apiMinimumWallThickness` | `(D m) → t m` = 0.00635 + D/100 | API RP 2A-WSD §6.10.6 |
| `monopileSizingScreen` | `({ D, t m, dOverTLimit = 60 }) → { tMin m, wallOk, dOverT, bucklingCheckRequired, area m², inertia m⁴, massPerLength kg/m }` | screening hook |

## driving.js

| Export | Signature → result (units) | Reference |
| --- | --- | --- |
| `GENERIC_HAMMERS` | `{ ratedEnergy J, ramMass kg, efficiency }` for six generic hydraulic hammers | generic demonstration values, not manufacturer data |
| `impactVelocity` | `({ ratedEnergy J, ramMass kg, efficiency }) → m/s` | energy balance |
| `waveEquation` | `({ hammer, cushionStiffness N/m, helmetMass kg, pile { length m, area m² ｜ D, t, E, rho, nSegments }, shaftRu N[], toeRu N, shaftQuake, toeQuake m, shaftDamping, toeDamping s/m, dtFactor = 0.9, pileDamping, duration s, history }) → { set m, blowsPerMetre 1/m, maxCompressiveStress Pa, maxTensileStress Pa, maxToeDisplacement m, enthru J, impactEnergy J, energyDrift, waveSpeed m/s, dt s, impedance N·s/m, history }` | Smith (1960) |
| `rodImpactForce` | `({ ramMass kg, impactVelocity m/s, area m², E, rho, t s }) → { force N, stress Pa, peakStress Pa, impedance }` | analytic rigid mass on elastic rod |
| `clayPileCapacityFactor` | `(OCR) → Fp` = 0.5 OCR^0.3 | Semple & Gemeinhardt (1981) |
| `srdProfile` | `({ profile, D, t m, Lmax m, dz m, upperPlugged, frictionFactorAt }) → { z m, fUnit Pa, qUnit Pa, isClay, cases: { coringLower, coringUpper, pluggedLower, pluggedUpper: { total, shaft, toe N, shaftPerLength N/m } } }` | Stevens, Wiltsie & Turton (1982) |
| `bearingGraph` | `({ …waveEquation, srd N[], shaftFraction, penetration m, distribution }) → { srd N, blowsPerMetre, blowsPer025m, set m, maxCompressiveStress, maxTensileStress Pa, enthru J }` | Smith (1960) |
| `drivabilityAnalysis` | `({ profile, hammer, pile, srdCase, depths m[], quake, damping, refusalBlowsPer025m = 250, … }) → { depth m, srd N, blowsPerMetre, blowsPer025m, stresses Pa, enthru J, refusal, refusalDepth m, totalBlows, fatigueBlocks [{ depthFrom, depthTo m, stressRange Pa, blows }] }` | Smith (1960) + Stevens et al. (1982) |
| `selfWeightPenetration` | `({ z m[], srd N[], pileWeight N, hammerWeight N }) → { swpPile m, swpPileHammer m, runIntervals, runRisk }` | vertical equilibrium |
| `pileSubmergedWeight` | `({ D, t, length m, rho, rhoWater kg/m³ }) → N` | Archimedes |

## cpt.js

| Export | Signature → result (units) | Reference |
| --- | --- | --- |
| `correctedConeResistance` | `({ qc Pa, u2 Pa, a }) → qt Pa` | Lunne, Robertson & Powell (1997) |
| `cptRatios` | `({ qt, fs, sigmaV0, sigmaV0Eff, u2, u0 Pa }) → { Rf %, Bq, Qt, Fr %, qnet Pa }` | Robertson (1990) |
| `soilBehaviourType`, `sbtZone` | `({ qt, fs, sigmaV0, sigmaV0Eff, pa }) → { Ic, Qtn, Fr %, n, zone, description, group }`; `(Ic) → { zone, description, group }` | Robertson (2009); Robertson & Wride (1998) |
| `relativeDensity` | `({ qc Pa, sigmaV0Eff Pa, method 'baldi'｜'jamiolkowski', C0 = 157, C1 = 0.55, C2 = 2.41, a = 0.268, b = 0.675 }) → Dr (fraction)` | Baldi et al. (1986); Jamiolkowski et al. (2003) |
| `frictionAngleDeg` | `({ qc, sigmaV0Eff, method 'kulhawy-mayne'｜'robertson-campanella' }) → deg` | Kulhawy & Mayne (1990); Robertson & Campanella (1983) |
| `undrainedShearStrength` | `({ qt, sigmaV0 Pa, Nkt }) → su Pa` | Lunne et al. (1997) |
| `overconsolidationRatio` | `({ qt, sigmaV0, sigmaV0Eff, k = 0.33 }) → OCR` | Kulhawy & Mayne (1990) |
| `smallStrainShearModulus` | `({ method 'rix-stokoe'｜'robertson', qc ｜ qt, sigmaV0, sigmaV0Eff, Ic, rho kg/m³, A = 1634, nq, ns, c1 = 0.55, c0 = 1.68 }) → { Gmax Pa, Vs m/s }` | Rix & Stokoe (1991); Robertson (2009) |
| `unitWeightEstimate` | `({ qt Pa, Rf %, gammaW = 9810 N/m³ }) → γ N/m³` | Robertson & Cabal (2010) |
| `cptToLayers` | `({ z m[], qc Pa[], fs Pa[], u2 Pa[], a, Nkt, gammaW, minThickness m, sandThresholdIc = 2.6, clayThresholdIc = 2.95, eps50 }) → { layers (for createProfile), points }` | Ic-threshold layering |
| `sptCorrection` | `({ N, energyRatioPct, CB, CS, CR, sigmaV0Eff Pa, CNmax = 1.7 }) → { N60, CN, N160 }` | Skempton (1986); Liao & Whitman (1986) |

## siteinvestigation.js

From site-investigation records to the soil values the suites share. Reading and reduction only: the correlations
are those of `cpt.js`.

| Export | Signature → result (units) | Reference |
| --- | --- | --- |
| `parseCptTable` | `(text, { depth, qc, fs, u2 } fallback units) → { z m[], qc Pa[], fs Pa[], u2 Pa[]｜null, qcIsQt, format, units, skipped, warnings[] }`; reads GEF files (`#COLUMNINFO` quantity numbers 1/11 depth, 2 cone resistance, 3 local friction, 6 pore pressure u2; `#COLUMNVOID`), the `SCPT` group of AGS 3/4 files (`SCPT_DPTH`, `SCPT_RES`, `SCPT_FRES`, `SCPT_PWP2`), and delimited text (comma, semicolon, tab, blanks) with named columns — or none, taken as depth, qc, fs, u2. Units in the header win over the fallback. Rows without depth > 0, cone resistance > 0 or a sleeve reading are skipped and counted; rows are sorted and repeated depths dropped | file formats: GEF-CPT-Report (CUR/Deltares); AGS data format |
| `PRESSURE_UNITS`, `LENGTH_UNITS` | factors to Pa (`Pa`, `kPa`, `MPa`, `bar`, `kgf/cm²`, `tsf`, `psi`) and to m (`m`, `cm`, `mm`, `ft`) | — |
| `fitStrengthLine` | `(z m[], su Pa[]) → { su0 Pa, suGrad Pa/m, n, r2, throughOrigin }`; least squares; a negative mudline value is replaced by the best line through the origin, a negative slope by the mean | — |
| `cptDesignValues` | `({ z, qc, fs, u2, qcIsQt, a = 0.8, Nkt = 15, gammaW = 10052, phiMethod, drMethod }) → { su0, suGrad, suRange (Nkt 10 and 20), suFit, phiDeg, phiRange [min, median, max], Dr, gammaSub, gammaSubClay, gammaSubSand, share { clay, silt, sand }, depth, layers, points }`: strength line through su = (qt − σv0)/Nkt of every clay or silt reading (at least 3); φ′ = lower quartile of the sand readings (at least 3), each limited to 20–45°; γ′ of the soil that makes up most of the sounding | Lunne, Robertson & Powell (1997); Robertson (2009); Kulhawy & Mayne (1990); Jamiolkowski et al. (2003); Robertson & Cabal (2010) |
| `boreholeDesignValues` | `(layers [{ top, bottom m, soil 'clay'｜'sand', suTop, suBottom Pa, phiDeg, gammaSub N/m³ }]) → { su0, suGrad, suFit, phiDeg, gammaSub, share, depth }`: strength line through the clay layers sampled every 0.25 m; lowest sand friction angle; thickness-weighted γ′ of the main soil | — |

Notes: the lower quartile as the characteristic friction angle and the least-squares strength line are reductions
chosen for screening a single sounding; a design profile is a cautious estimate over all soundings of a site and
is the geotechnical engineer's to draw. Median grain size is not derived from a cone test.

## Not implemented, and items to confirm against the printed sources

* Reese, Cox & Koop (1975) stiff clay with free water (coefficients As, Ac exist only as figures). Stiff clay uses
  the Matlock shape with an explicit cyclic factor; the Reese & Welch no-free-water curve is offered separately.
* Fugro-05 and NGI-05 CPT methods; CPT end bearing other than offshore UWA-05; Alm & Hamre (2001) SRD;
  vibratory driving; underwater noise; residual-stress wave-equation analysis.
* Shadlou & Bhattacharya stiffness expressions and the Arany et al. closed-form frequency coefficients: the
  stiffness is computed numerically (`monopileStiffness`) and the frequency by the beam eigenproblem.
* Poulos & Davis chart-based elastic influence factors: the Randolph (1981) and Randolph & Wroth (1978) closed
  forms are provided as the elastic-continuum references.
* LeBlanc et al. Tb and Tc are chart functions and are arguments.
* Layering corrections to p–y curves (Georgiadis equivalent depth) are not applied.
* Values held from the published sources that should be confirmed against the printed documents before
  certification use: the k(φ′) fit to the API figure; the simplified ICP-05 / offshore UWA-05 coefficient sets;
  the AASHTO p-multiplier table; the 1.5 D / 6 D default scour depths; the free-head coefficient (≈ 0.4) of the
  Randolph lateral interaction factor (argument); the Stevens et al. upper-bound plugged multipliers
  (arguments; sand 1.3 / 1.5 by default, clay 1.0); clause numbers of API RP 2GEO (cited at section level).
* Wave equation: below the stability limit a cushionless rigid impact overshoots the peak segment stress by up
  to about 20 % (analytic peak ρ c v0); head force, ENTHRU and blow count are unaffected (see `waveEquation`).
* Observed convergence of the BNWF solver with element count: second order (lumped nodal springs), verified
  against Hetényi in `tests/engine/geotech.pile.test.js`.
