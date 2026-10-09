# `src/engine/xsection` — cross-section mechanics of flexible pipes, umbilicals and cables

SI units throughout (m, N, Pa, kg, s, K; angles in rad unless the name ends in `Deg`; temperatures in °C only
where the name ends in `C`). Lay angle `alpha` is measured from the product axis; `lay` is +1/−1; `psi` is the
angular position from the extrados. Every function throws `RangeError` naming the offending argument.
One line per export: signature → returned fields (units) — reference.

## helix.js — helically wound layers

- `helixGeometry({ R, alpha, wireWidth, nWires | fillFraction })` → `pitch`, `wireLengthPerPitch` (m), `wireLengthPerMetre` (–), `nWires`, `fillFraction`, `gap` (m), `curvature` (1/m), `torsion` (rad/m) — geometry of the circular helix; F = n b/(2πR cos α).
- `rectWireSection({ width, thickness })` → `A` (m²), `In`, `It`, `J` (m⁴) — Roark Table 10.1 rectangle torsion constant.
- `roundWireSection({ diameter })` → `A` (m²), `In`, `It`, `J` (m⁴).
- `wireStrain({ alpha, R, lay, epsAxial, twist, dR })` → `strain` (–) — ε_w = cos²α ε + lay R sinα cosα τ + sin²α ΔR/R (derived; Féret & Bournazel 1987, Witz & Tan 1992).
- `helixLayerStiffness({ n, E, A, R, alpha, lay })` → `kAxial` = nEA cos³α (N), `kTorsion` = nEA R² sin²α cosα (N m²), `kCoupling` = lay nEA R sinα cos²α (N m), `kAxialRadial` (N), `kTorsionRadial` (N m), `kRadial` (N) — Witz & Tan (1992); thin-layer radial coupling derived in the file.
- `helixContactPressure({ n, wireForce, R, alpha })` → `pressure` = n F sin²α/(2πR² cosα) (Pa), `lineLoadPerWire` (N/m), `radialForcePerLength` (N/m) — derived from wire curvature; Féret & Bournazel (1987).
- `wireCurvatureChange({ alpha, kappa, psi, curve })` → `dKn`, `dKt` (1/m), `dTw` (rad/m) — loxodromic (cos⁴α, (1+sin²α)cosα, sinα cos³α) and geodesic (3cos²α, 0, cotα(1−3sin²α)) paths, derived on the torus; Sævik (1992/2011), Féret & Bournazel (1987).
- `wireBendingStress({ E, thickness, width, alpha, kappa, psi, curve })` → `sigmaNormal`, `sigmaTransverse`, `sigmaCorner` (Pa), `dKn`, `dKt`, `dTw` — Sævik (2011).
- `helixBendingStiffness({ n, E, A, R, alpha, In, It, J, G, curve })` → `stick` = nEA R² cos³α/2, `slip` (N m²) — strain-energy derivation; Witz & Tan (1992).
- `criticalCurvature({ E, A, R, alpha, n, frictionPerLength | mu, qOuter, qInner })` → `frictionPerLength` (N/m), `kappaCr` = f/(EA cos²α sinα), `kappaFullSlip` = (π/2)κcr (1/m), `sigmaStickAtOnset`, `sigmaFrictionMax` = πRf/(2 sinα A) (Pa), `momentAtOnset`, `frictionMoment` (N m), `slipAmplitudePerCurvature` (m²) — shear-lag balance; Féret & Bournazel (1987), Sævik (2011).
- `wireFrictionStress({ E, A, R, alpha, kappa, psi, frictionPerLength })` → `sigma` (Pa), `state`, `psi0` (rad), `slipped` — stick / partial-slip / full-slip distribution (derived).
- `wireStressAtPosition({ E, A, R, alpha, thickness, width, kappa, psi, frictionPerLength, sigmaAxial, curve })` → `sigmaAxial`, `sigmaFriction`, `sigmaNormal`, `sigmaTransverse`, `sigmaMaxCorner`, `sigmaMinCorner` (Pa), `state`.
- `momentCurvatureLoop({ EIstick, EIslip, frictionMoment, kappaAmp, nPoints })` → `kappaY` (1/m), `Mmax` (N m), `energyPerCycle` (J/m), `dampingRatio` (–), `secantStiffness` (N m²), `kappa[]`, `moment[]` — bilinear hysteresis, equivalent damping E_d/(4πE_s) (Chopra §3.9).

## flexible.js — unbonded flexible pipe (API 17J / 17B)

- `API17J_FACTORS` (constant) — permissible utilisation: tensile armour 0.67 / 0.85 / 0.91, pressure armour 0.55 / 0.85 / 0.91 (normal / abnormal / hydrotest) — API Spec 17J Table 6 as recalled; overridable.
- `SHEATH_STRAIN_LIMITS` (constant) — PE, PA 7.7 %; PVDF 7.0 % static, 3.5 % dynamic — API Spec 17J §5.3.1 as recalled.
- `carcassBucklingFactor(depth)` → `factor` (–): 0.67 (≤ 300 m) to 0.85 (≥ 900 m), linear between — API Spec 17J Table 6 as recalled.
- `buildFlexibleStack({ innerDiameter, layers })` → `innerRadius`, `outerRadius` (m), `layers[]` (with `rIn`, `rOut`, `R`, `kind`, `nA` m², `n`, `alpha`), `sealInner`, `sealOuter` — API RP 17B layer terminology.
- `axisymmetricResponse(stack, { tension, tensionType, pInt, pExt, annulus, pAnnulus, torque, axialStrain, twist, deltaT, initialGaps })` → `axialStrain` (–), `twist` (rad/m), `wallForce`, `endCapForce`, `tensionEffective` (N), `torque` (N m), `layers[]` (forces N, N m, N/m; displacements m; pressures and stresses Pa), `interfaces[]` (`contactPressure` Pa, `gap` m, `closed`), `stiffness` {K11 N, K12 = K21 N m, K22 N m²}, `residuals`, `iterations` — layer-compatibility system with active-set contact, derived in the file header (Féret & Bournazel 1987; Witz & Tan 1992).
- `layerUtilisations(stack, response, { condition, factors })` → `layers[]` (`stress`, `capacity` Pa, `factor`, `utilisation`), `governing` — API Spec 17J Table 6.
- `burstPressure(stack, { pExt, tension, basis })` → `pressureArmour`, `tensileArmour`, `burst` (Pa), `governing` — API RP 17B §8 failure modes; elastic scaling of the axisymmetric solution.
- `axialFailureTension(stack, { pInt, pExt, basis, twist })` → `tension` (N), `layer`.
- `torsionLimits(stack, { loads, factor })` → `tightening`, `opening` {`stressLimit`, `separation`, `governing` (N m), `mode`} — wire-stress and gap-opening (bird-caging onset) criteria.
- `bendingProperties(stack, { response, mu, curve })` → `EIstick`, `EIslip` (N m²), `frictionMoment` (N m), `kappaCrMin` (1/m), `layers[]` — helix.js formulas with solved interface pressures; Sævik (2011).
- `minimumBendRadius(stack, { storageFactor, staticFactor, dynamicFactor })` → `lockingRadius`, `sheathStaticRadius`, `sheathDynamicRadius`, `storage`, `static`, `dynamic` (m), `lockingLayer` — gap-closure geometry r_out/lockStrain; API Spec 17J §5.3.1 factors 1.1 / 1.0 / 1.5 as recalled.
- `massAndWeight(stack, { rhoContent, rhoSeawater, g, annulusFlooded })` → `massEmpty`, `massSeawaterFilled`, `massOperating`, `annulusWaterMass`, `displacedMass` (kg/m), `contentArea`, `annulusVoidArea` (m²), `weightAirEmpty`, `submergedEmpty`, `submergedFlooded`, `submergedOperating` (N/m).
- `equivalentBeamProperties(stack, { loads, mu, weight })` → `EA`, `EAfixed` (N), `GJ`, `GJfixed` (N m²), `coupling` (N m), `EIstick`, `EIslip` (N m²), `kappaCr` (1/m), `frictionMoment` (N m), `outerDiameter`, `innerDiameter` (m), `mass`.
- `thermalProfile(stack, { TinC, ToutC, hIn, hOut })` → `heatLoss` (W/m), `UvalueInner`, `UvalueOuter` (W/m²K), `resistance` (m K/W), `interfaceTempsC[]`, `layerMeanTempsC[]` (°C) — composite-cylinder conduction (Incropera §3.3).
- `pressureThermalElongation(stack, { pInt, pExt, deltaT, annulus })` → `strainPressure`, `strainThermal`, `strainTotal` (–), `twistPressure`, `twistThermal` (rad/m).
- `exampleFlexibleDesign(id)` → generic demonstration stacks `'riser-6in'`, `'flowline-10in'` (flagged `generic: true`; not manufacturer designs).

## flexcollapse.js — collapse and armour instability

- `interlockedRingStiffness({ E, Imin, pitch, nTendons, K })` → `EI` = K n E I_min/L_p (N m per m), `tEq` (m) — API RP 17B §8.4 equivalent ring stiffness (K is profile specific).
- `solidRingStiffness({ E, thickness, fill, nu })` → `EI` (N m) — E t³/12, optional plane-strain factor.
- `ringElasticCollapse({ EI, R })` → `pCr` = 3EI/R³ (Pa) — Levy (1884); Timoshenko & Gere Ch. 7.
- `combinedElasticCollapse({ rings })` → `pCr` = Σ 3EI_i/R_i³ (Pa), `contributions[]` — API RP 17B §8.4.
- `ovalisedRingCollapse({ pCr, yieldStress, R, thickness, sectionModulus, ovality })` → `pCollapse`, `pYield`, `pCr` (Pa), `w0` (m) — Timoshenko first-yield ovalised ring, derived in the file.
- `collapseAssessment({ carcass, pressureArmour, annulus, ovality, curvatureFactor, pDesign, factor, rho, g })` → `pElastic`, `pCollapse` (Pa), `collapseDepth`, `allowableDepth` (m), `utilisation`, `mode` — wet (carcass alone) / dry (carcass + pressure armour).
- `birdcagingCheck({ compression, layers, tape, sheath, pExt, pAnnulus, annulusFlooded })` → `demandPressure`, `capacityPressure`, `tapeCapacity`, `sheathCapacity`, `pressureSupport` (Pa), `utilisation`, `criticalCompression` (N) — outward pressure C tan²α/(2πR²), derived.
- `wireFoundationBuckling({ E, In, k | restraintE, restraintThickness, restraintR, R, alpha, n })` → `Pcr` = 2√(kEI) (N), `halfWavelength` (m), `k` (N/m²), `layerCompression` (N) — Hetényi, Beams on Elastic Foundation.
- `lateralBucklingScreen({ E, It, R, alpha, n, lengthFactor, compression })` → `Pcr` (N), `bucklingLength` (m), `layerCompression` (N), `utilisation` — Euler screening formula (derived; limitations stated in the JSDoc).
- `overbendCheck({ curvature, lockingRadius })` → `utilisation`, `locked`, `marginRadius` (m).

## flexfatigue.js — armour fatigue, wear, end fitting

- `armourStressRange({ E, A, R, alpha, thickness, width, frictionPerLength, kappaRange, curve })` → `frictionRange`, `stickRange`, `normalBendingRange`, `transverseBendingRange`, `rangeExtrados`, `rangeNeutralAxis`, `governing` (Pa), `slipped` — Sævik (2011).
- `meanStressCorrection({ stressRange, meanStress, uts, method })` → `equivalentRange` (Pa), `factor` — Goodman / Gerber.
- `snCycles({ stressRange, curve })` → `cycles` — N = a S⁻ᵐ (S in Pa), optional second slope and cut-off.
- `genericArmourSN({ uts, environment, environmentFactor, f, marin })` → `{ a, m, limit?, generic: true, label }` — textbook estimate (Shigley §6-7/6-8), unqualified; corrosive environments need a caller-supplied knock-down.
- `minerDamage({ bins, curve, uts, meanStressMethod, scf, corrosionFactor, dff, exposureYears })` → `damage`, `damagePerYear` (1/yr), `life`, `lifeWithDff` (yr), `bins[]` — Palmgren–Miner; DFF 10 per API Spec 17J.
- `armourFatigueCoefficients(stack, layerIndex, { tensionLow, tensionHigh, mu, loads })` → `meanStress0` (Pa), `meanStressPerTension` (Pa/N), `friction0` (N/m), `frictionPerTension` (1/m), `wire`.
- `pressureCycleStressRange(stack, layerIndex, { pLow, pHigh, loads })` → `stressRange`, `meanStress`, `stressLow`, `stressHigh` (Pa) — pressure-armour cycling.
- `armourFatigueDamage({ wire, bins, meanStress0, meanStressPerTension, friction0, frictionPerTension, curve, uts, meanStressMethod, scf, corrosionFactor, dff, exposureYears, wireCurve })` → `minerDamage` fields + `stress[]` (Pa, N/m) — bins of (tension N, curvature range 1/m, cycles).
- `wireSlipPerCycle({ R, alpha, kappaRange })` → `slipRange` = R² cos²α Δκ/sinα, `slidingPerCycle` (m) — Féret & Bournazel (1987), derived.
- `archardWear({ wearCoefficient, contactPressure, slidingPerCycle, cycles, thickness })` → `wearDepth`, `remainingThickness` (m), `stressFactor`, `wornThrough` — Archard (1953), h = k p s.
- `sheathCreepIntoGaps({ gapWidth, intrusionDepth, pitch, sheathThickness, fillShape, allowableReduction })` → `thicknessLoss` (m), `reduction`, `utilisation` — volume balance; 30 % limit of API Spec 17J as recalled.
- `sheathStrainCheck({ rOut, kappaMax, strainLimit, kappaRange, strainRangeLimit })` → `strain`, `utilisation`, `strainRange`, `rangeUtilisation`.
- `endFittingAnchorage({ bondStress, perimeter, embeddedLength, wireArea, uts, nWires, coneHalfAngle, muCone, muWire })` → `bondCapacity`, `capacity`, `wireStrength` (N), `wedgeAmplification`, `wedgeCapacityFraction`, `efficiency`, `selfReleasing` — bond × embedded length and wedge statics.
- `multiaxialCriterion({ criterion, amplitude, mean, k, limit, nGrid })` → `equivalent`, `shearTerm`, `normalTerm` (Pa), `utilisation` — Sines, Crossland, Dang Van, Findley (Socie & Marquis 2000).

## annulus.js — annulus environment

- `R_GAS` (constant) — molar gas constant 8.314462618 J/(mol K) (CODATA 2018).
- `H2S_SOUR_LIMIT` (constant) — 300 Pa, ISO 15156-2 sour-service threshold 0.3 kPa (0.05 psi) H2S partial pressure.
- `permeabilityFromCm3(value)` → `Pe` (mol m/(m² s Pa)) from cm³(STP) cm/(cm² s bar).
- `arrheniusPermeability({ Pe0, Ea, T })` → `Pe` — Pe0 exp(−Ea/RT).
- `annulusPermeation({ components, ri, ro, length, T, areaFactor })` → `rates[]` (`Pe`, `molarRate` mol/s), `totalMolarRate` (mol/s), `moleFractions` — n = Pe Δf 2πL/ln(ro/ri) (Crank; API RP 17B).
- `annulusFreeVolume({ layers, length, extraAreaVoid })` → `voidArea`, `annulusArea` (m²), `volume` (m³), `voidFraction`.
- `annulusPressureBuildUp({ molarRate, volume, T, p0, pVent, Tstd, pStd })` → `pressureRate` (Pa/s), `timeToVent`, `ventInterval` (s), `ventRateStd` (m³/s std), `molesPerVent` (mol) — ideal gas.
- `waterSaturationPressure(TC)` → `pSat` (Pa) — Antoine equation, 1–100 °C.
- `annulusClassification({ moleFractions, pAnnulus, TC, sourLimit })` → `partialPressures` (Pa), `classification` (sour / sweet / benign), `h2sRatio`, `waterCondenses` — ISO 15156-2.
- `annulusPH({ pCO2, pK1, pKH, confinementShift, pHOverride })` → `pH`, `pHPureWater` — carbonic-acid equilibrium (Stumm & Morgan); confinement effect is an argument.
- `annulusFlooding({ freeVolume, floodedLength, breachDepth, floodedFraction, rho, g, pAtm })` → `waterVolume` (m³), `addedMassPerLength` (kg/m), `addedWeightPerLength` (N/m), `annulusPressure` (Pa abs), `consequences[]`.
- `annulusVolumeTest({ method, p1, p2, T, injectedStdVolume, refVolume, refPressure, designVolume, pressureDrift, Tstd, pStd })` → `freeVolume` (m³), `floodedFraction`, `leakMolarRate` (mol/s), `leakStdRate` (m³/s std) — ideal gas / Boyle.
- `rgdScreening({ pressureDrop, duration, qualifiedRate, qualifiedPressure, initialPressure })` → `rate` (Pa/s), `rateUtilisation`, `pressureUtilisation`, `risk` — against the qualified blistering-test envelope (API Spec 17J).

## umbilical.js — umbilical / dynamic cable assembler

- `UMBILICAL_UTILISATION` (constant) — { normal 0.67, abnormal 0.8, accidental 1.0 } placeholder utilisation factors in the ISO 13628-5 / API 17E pattern; replace with the governing edition's values.
- `umbilicalComponent(c)` → component with `A` (m²), `I` (m⁴), `EA` (N), `EI`, `GJ` (N m²), `mass` (kg/m), `boreArea` (m²).
- `assembleUmbilical(def, { rhoSeawater, g })` → `outerDiameter` (m), `massEmpty`, `massFilled` (kg/m), `weightAir`, `submergedEmpty`, `submergedFilled` (N/m), `K11` (N), `K12` (N m), `K22` (N m²), `EA`, `EAfixed`, `GJ`, `GJfixed`, `EIstick`, `EIslip`, `torquePerTension` (m), `twistPerTension`, `torqueBalanced`, `elements[]` — helical-tendon model (helix.js).
- `umbilicalLoadSharing(section, { tension, torque, twistRestrained })` → `axialStrain`, `twist` (rad/m), `torque` (N m), `elements[]` (`force` N, `share`, `stress` Pa, `lineLoad`, `contactInner`, `contactOuter` N/m).
- `umbilicalUtilisation(section, { tension, curvature, mu, eta, pExt, twistRestrained })` → `utilisation`, `governing`, `elements[]` (`vonMises`, `axialStress`, `bendingStress`, `frictionStress` Pa) — von Mises of tubes incl. pressure (Lamé), tension share, E d/2 κ and friction stress; ISO 13628-5 / API 17E.
- `umbilicalCapacityCurve(section, { eta, mu, pExt, nPoints, kappaMax, twistRestrained })` → `curvature[]` (1/m), `allowableTension[]` (N), `governing[]`, `kappaLimit`, `tensionLimit`.
- `umbilicalMbr(section, { tension, eta, mu, pExt, twistRestrained })` → `mbr` (m), `curvature` (1/m).
- `ringRadialLoadMoment({ P, R, nLoads })` → `momentAtLoad`, `momentBetween` (N m/m), `coefficient` — M = (PR/2)(1/θ − cot θ), θ = π/n; 0.3183 PR for two opposing loads (Roark Table 9.2).
- `umbilicalCrushCheck({ tubeOd, tubeWall, allowableStress, padLoad, nTracks, elementsUnderPad, alpha, tension, chuteRadius, mu, safetyFactor })` → `tubeLineLoadCapacity`, `allowablePadLoad`, `tubeLineLoadTensioner`, `chuteLineLoad` (N/m), `utilisationTensioner`, `utilisationChute`, `requiredTensionerLength`, `minimumChuteRadius` (m).
- `umbilicalCoreTemperature({ cores, ambientC, bundleResistance, sheath, hExternal, outerDiameter })` → `coreTempsC[]`, `sheathInnerTempC`, `surfaceTempC` (°C), `totalLoss` (W/m), `outerResistance` (K m/W) — series thermal resistances (IEC 60287 method, simplest form).
- `umbilicalTubeFatigue(section, { element, bins, curve, mu, scf, dff, exposureYears, twistRestrained })` → `minerDamage` fields + `ranges[]` (Pa).
- `tubePressureDrop({ flow, d, length, rho, viscosity, roughness })` → `dp` (Pa), `velocity` (m/s), `reynolds`, `frictionFactor`, `regime` — Darcy–Weisbach, 64/Re or Haaland (1983).
- `hydraulicResponseTime({ length, d, viscosity, bulkModulus, volumetricExpansion, fraction })` → `tau1`, `time` (s), `resistance` (Pa s/m³), `capacitance` (m³/Pa) — first mode of the laminar line diffusion equation, τ₁ = (4/π²) r c L².

## tubes.js — small-bore tubes, hoses, cores, fibres

- `vonMises({ s1, s2, s3, tau })` → `vonMises` (Pa).
- `lameStresses({ pInt, pExt, ri, ro, r, ends, axialForce })` → `radial`, `hoop`, `axial`, `vonMises` (Pa) — Lamé; Timoshenko & Goodier §28.
- `thinWallStresses({ pInt, pExt, D, t, barlow })` → `hoop`, `axial` (Pa).
- `thickWallLimitPressures({ ri, ro, yieldStress, uts })` → `pYield`, `pPlastic`, `pBurstEstimate` (Pa) — Lamé + von Mises; (2/√3) σ ln(ro/ri) (Hill).
- `b313PressureDesign({ D, S, t | P, E, W, Y, c })` → `P` (Pa), `t` (m), `Y`, `thickWall` — ASME B31.3 §304.1.2 Eq. (3a).
- `dnvBurstPressure({ D, t, fy, fu })` → `pb`, `fcb`, `yieldLimit`, `burstLimit` (Pa) — DNV-ST-F101 Eq. (5.8).
- `tubeCollapsePressure({ D, t, E, fy, nu, ovality, alphaFab, enforceMinimum })` → `pc`, `pel`, `pp` (Pa), `f0` — DNV-ST-F101 Eq. (5.10) cubic.
- `tubeCombinedUtilisation({ D, t, pInt, pExt, tension, curvature | moment, E, yieldStress, eta, ends })` → `vonMisesMax`, `bendingStressOuter`, `axialMembrane`, `hoopBore` (Pa), `utilisation`, `location`.
- `proofTestPressure({ designPressure, proofFactor, ri, ro, yieldStress })` → `testPressure`, `vonMisesBore` (Pa), `utilisation` — 1.5 × design pressure default (argument).
- `tubeBendThinning({ t, D, bendRadius })` → `tExtrados`, `tIntrados` (m), `thinningFraction`, `bendingStrain` — t R/(R ± r), derived.
- `hoseRating({ burstPressure, burstRatio, volumetricExpansion, length, pressureChange })` → `designPressure` (Pa), `expansionVolume` (m³).
- `insulationResistance({ resistivity, conductorDiameter, insulationDiameter, length })` → `resistance` (Ω), `resistanceKm` (Ω km) — ρ ln(D/d)/(2πL).
- `looseTubeStrainWindow({ layLength, pitchRadius, clearance, efl, cableStrain, proofStrain, allowableFraction })` → `tensileWindow`, `contractionWindow`, `fibreStrain`, `allowableStrain` (–), `utilisation` — helix geometry, derived.

## cablemech.js — submarine power and telecom cables

- `GENERIC_CABLE_MATERIALS` (constant) — indicative handbook E (Pa), rho (kg/m³), alphaT (1/K), flagged generic.
- `GENERIC_TELECOM_CABLES` (constant) — LW, LWP, SA, DA indicative diameter (m), mass (kg/m), weight in water (N/m), NTTS/NOTS/NPTS (N), flagged generic.
- `buildCable(def, { rhoSeawater, g })` → `outerDiameter`, `coreDiameter` (m), `mass` (kg/m), `weightAir`, `submergedWeight` (N/m), `specificGravity`, `EA`, `EAfixed` (N), `GJ`, `GJfixed` (N m²), `coupling` (N m), `torquePerTension` (m), `torqueBalanced`, `EImin`, `EImax` (N m²), `conductorArea` (m²), `layers[]` — CIGRE TB 623 in intent.
- `cablePullingTension(cable, { conductor, conductorStress, armourUtilisation })` → `conductorLimit`, `armourLimit`, `armourStrength` (N) — 50 N/mm² Cu / 30 N/mm² Al defaults (arguments).
- `cableHandlingLimits({ outerDiameter, installFactor, staticFactor, tension, bendRadius, sidewallLimit, squeeze, squeezeLimit, nTracks, mu, safetyFactor })` → `mbrInstall`, `mbrStatic`, `minRadiusForSidewall`, `gripLength` (m), `sidewallPressure` (N/m), utilisations.
- `cableImpactClass({ energy, bands })` → `label`, `band`, `exceedsAll`.
- `cableCoiling(cable, { coilRadius })` → `coilable`, `twistRate` = 1/R (rad/m), `armour[]` (strain), `note`.
- `cableLayTension({ submergedWeight, depth, bottomTension, dynamicFactor })` → `staticTopTension` = H + w d, `designTension` = 1.3 w d + H (N), `departureAngle` (rad), `suspendedLength`, `layback`, `touchdownRadius` (m) — catenary; Electra 171 / CIGRE TB 623 form as recalled.
- `cableOnBottomStability({ submergedWeight, outerDiameter, velocity, mu, Cd, Cl, rho, safetyFactor, mass })` → `drag`, `lift`, `resistance`, `requiredWeight` (N/m), `utilisation`, `criticalVelocity` (m/s), `specificGravity`, `stable`.
- `conductorThermalThrust({ area, Eeff, alphaT, deltaT, snakeWavelength, initialAmplitude, residualThrust })` → `thrust` = E A α ΔT (N), `freeStrain`, `snakeAmplitude` (m), `snakeCurvature` (1/m).
- `telecomTensionUtilisation({ transientTension, operatingTension, permanentTension, ntts, nots, npts })` → `transient`, `operating`, `permanent`, `governing` — ITU-T G.972 definitions.
- `telecomSinkAngle({ weightWater, outerDiameter, shipSpeed, Cd, rho, depth })` → `hydrodynamicConstant` (rad m/s), `hydrodynamicConstantDegKnots`, `sinkAngle`, `sinkAngleSmall` (rad), `sinkAngleDeg`, `layback` (m), `sinkTime` (s) — Zajac (1957), α ≈ H/V.
- `telecomSlack({ sinkAngle, seabedSlope, shipSpeed, slack | payoutSpeed })` → `payoutSpeed` (m/s), `slack`, `slackPercent`, `speedRatio`, `zeroSlackRatio` — Vc/Vs = [(1+e) sinα + sinβ]/sin(α+β), Zajac (1957).

## burial.js — protection by burial and rock

- `BPI_CLASSES` (constant) — BPI 1/2/3 threat definitions (Mole et al. 1997; Allan 1998).
- `burialProtectionIndex({ referenceDepth, burialDepth, targetBpi })` → `bpi`, `requiredDepth` (m), `protectedAgainst`.
- `anchorPenetration({ flukeLength | anchorMass, flukeCoefficient, flukeAngleDeg, soilFactor })` → `flukeLength`, `penetration` (m) — CBRA structure; all factors arguments.
- `cbraStrikeProbability({ classes, burialDepth, pTraffic, pWd })` → `pStrike` (1/yr), `returnPeriod` (yr), `classes[]` — Carbon Trust CBRA (2015): P = P_traffic P_wd Σ D/(V·8760 h) P_incident.
- `depthOfLowering({ classes, targetProbability, fishingPenetration, margin, pTraffic, pWd })` → `recommendedDepth` (m), `governing`, `pStrike`, `returnPeriod`, `pStrikeUnburied`.
- `GENERIC_FISHING_PENETRATION` (constant) — editable indicative gear × soil penetration table (m), flagged generic.
- `fishingGearPenetration({ gear, soil, safetyFactor, table })` → `penetration`, `design` (m), `generic`.
- `trenchingRate({ power, efficiency, specificEnergy | strength, strengthFactor, trenchWidth, trenchDepth | trenchArea })` → `rate` (m/s), `ratePerHour` (m/h), `specificEnergy` (J/m³), `trenchArea` (m²) — specific-energy model (Teale 1965 concept).
- `ploughTowForce({ soil, depth, ploughWeight, gammaSub, speed, Cw, Cs, Cd, Cc, su, phiDeg, shareWidth, mu })` → `empirical`, `staticPart`, `ratePart`, `rankine` (N) — Cathie & Wintgens (2001) form with caller coefficients; Rankine passive wedge.
- `naturalInfillTime({ trenchArea, transportRate, trapEfficiency })` → `time` (s), `timeDays`.
- `rockBerm({ productDiameter, cover, crestWidth, slope, length, bulkDensity | rockDensity, porosity, overdump, requiredAbsorptionCover })` → `height`, `baseWidth` (m), `area` (m²), `volume` (m³), `mass` (kg), `tonnes`, `coverUtilisation`.
- `rockSizeShields({ bedShear | shearVelocity, shieldsCritical, rockDensity, rhoW, g, d50 })` → `requiredD50` (m), `bedShear` (Pa), `utilisation` — Shields (1936); CIRIA C683.
- `mattressCount({ length, mattressLength, overlap, endOverrun, coverWidth, mattressWidth, mattressMass })` → `perRow`, `rows`, `count`, `totalMass` (kg).
- `burialThermalRecord({ depthToTop, outerDiameter, soilResistivity, backfillResistivity })` → `burialDepthAxis` (m), `soilResistivity`, `backfillResistivity` (K m/W), `u`, `T4` (K m/W) — IEC 60287-2-1 single buried cable.
- `sandwaveExposure({ height, wavelength, migrationRate, initialCover, targetCover, phase })` → `time` (unit of the migration rate), `minimumCover` (m), `exposed`, `period` — kinematic sinusoidal bedform.

## ancillaries.js — bend stiffeners, buoyancy, pull-in, clamps

- `bendStiffener({ tension, angle, length, rootDiameter, tipDiameter, innerDiameter, Epu, EIpipe, linear, nSteps, mbr, strainLimit })` → `s[]` (m), `curvature[]` (1/m), `angle[]` (rad), `moment[]` (N m), `rootCurvature`, `tipCurvaturePipe`, `maxCurvature` (1/m), `maxCurvatureAt` (m), `rootMoment`, `rootMomentStiffener` (N m), `rootShear`, `rootAxial` (N), `maxPuStrain`, `curvaturePerAngle`, `mbrUtilisation`, `strainUtilisation` — elastica with varying EI, RK4 shooting; Boef & Out (1990).
- `sizeBendStiffener({ loadCases, lengths, rootDiameterMax, tipDiameter, innerDiameter, Epu, EIpipe, mbr, strainLimit, nDiameters })` → `feasible`, `length`, `rootDiameter` (m), `volume` (m³), `maxCurvature` (1/m), `maxPuStrain`, `utilisation`.
- `bendRestrictor({ pitch, lockAngle, nElements, mbr, tension, angle, momentCapacity, shearCapacity })` → `lockingRadius` = p/(2 sin(δ/2)), `lockedLength` (m), `totalAngle` (rad), `rootMoment` = T R (1 − cos θ) (N m), `rootShear` (N), utilisations, `coversAngle`.
- `bellmouth({ radius, maxAngle, throatDiameter, mbr, tension })` → `length`, `flare`, `exitDiameter`, `arcLength` (m), `mbrUtilisation`, `lineLoad` (N/m).
- `GENERIC_FOAM_DENSITY` (constant) — editable indicative syntactic-foam density (kg/m³) vs depth rating (m), flagged generic.
- `buoyancyModules({ outerDiameter, innerDiameter, moduleLength, foamDensity | depthRating, table, hardwareWeight, lossAllowance, requiredBuoyancy, sectionLength, clampForce, mu, axialLoad, slipSafetyFactor, rhoW, g })` → `foamDensity` (kg/m³), `foamVolume` (m³), `netBuoyancy` (N), `pitch` (m), `count`, `providedBuoyancy` (N/m), `slipCapacity` (N), `slipUtilisation`.
- `midWaterArch({ netBuoyancy, T1, angle1, T2, angle2, tetherCapacity, nTethers })` → `tetherVertical`, `tetherHorizontal`, `tetherTension`, `perTether` (N), `tetherAngle` (rad), `slack`, `utilisation`.
- `capstan({ tensionIn, mu, wrapAngle })` → `tensionOut` (N), `ratio` — T₂ = T₁ e^{μθ}.
- `pullInTension({ backTension, segments, capacity })` → `pullTension` (N), `tensions[]` (N), `maxSidewall` (N/m), `utilisation` — cable-pulling equations (Rifenburg 1953).
- `clampSpacing({ load, momentAllowable, clampCapacity, EI, deflectionLimit, length })` → `spacing`, `spacingMoment`, `spacingClamp`, `spacingDeflection` (m), `governing`, `count`.
- `impactProtection({ impactEnergy | mass, dropHeight | velocity, ratedEnergy, length, shellLength, g })` → `impactEnergy` (J), `utilisation`, `shells`.
- `hangOffClamp({ tension, nBolts, boltPreload, mu, diameter, clampLength, faces, safetyFactor, allowablePressure })` → `clampingForce`, `slipCapacity` (N), `slipUtilisation`, `contactPressure` (Pa), `pressureUtilisation`.
- `terminationLandingLoad({ mass, addedMass, submergedWeight, velocity, stiffness, capacity })` → `dynamicForce`, `peakLoad` = W + √(W² + k m v²) (N), `daf`, `kineticEnergy` (J), `utilisation`.

## Items to confirm before design use

Values recalled from standards and exposed as overridable arguments or constants: API Spec 17J Table 6 utilisation
factors and the depth-dependent carcass factor; 17J sheath strain limits and the 1.1 / 1.0 / 1.5 MBR factors; 17J
30 % creep limit; API RP 17B default ovality 0.2 %; ISO 13628-5 / API 17E utilisation factors, proof factor 1.5
and hose burst ratio 4; the 1.3 w d + H laying test tension; conductor pulling stresses. Clause numbers quoted
for API RP 17B and API Spec 17J are indicative. Generic tables (`exampleFlexibleDesign`, `genericArmourSN`,
`GENERIC_*`) are demonstration data only.
