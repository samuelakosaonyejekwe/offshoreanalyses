# `src/engine/coastal` — coastal hydraulics and sediment engines

SI units throughout (m, kg, s, N, Pa); angles in radians unless the name ends in `Deg`; `g` defaults to
9.80665 m/s². Forces and moments on walls are per metre run. Sediment transport rates are solid volume per unit
width [m²/s]; overtopping discharges are [m³/s per m]. Invalid input throws `RangeError`.

## transform.js

- `linearWave({ T, h, g? }) -> { k, L, C, n, Cg, L0, C0, omega }` — Airy dispersion and kinematics [rad/m, m, m/s, -, m/s, m, m/s, rad/s]. Dean & Dalrymple (1991).
- `shoalingCoefficient({ T, h, g? }) -> Ks` — linear shoaling from deep water, Ks = √(Cg0/Cg) [-]. Dean & Dalrymple (1991) §4.8.
- `refractionSnell({ T, h, theta0Deg, g? }) -> { theta, thetaDeg, Kr }` — Snell's law over straight parallel contours [rad, deg, -]. Dean & Dalrymple (1991) §4.8.
- `transformProfile({ H0, T, theta0Deg?, depths, gammaB?, g? }) -> { Ks, Kr, H, thetaDeg, L, broken }` — Ks·Kr along a depth profile with optional depth-limited cap [-, -, m, deg, m, flag].
- `breakingMiche({ T, h, g? }) -> { Hmax, L, steepnessMax }` — H/L = 0.142 tanh(kh) [m, m, -]. Miche (1944).
- `breakingMcCowan({ h, index? }) -> Hb` — Hb = 0.78 hb [m]. McCowan (1894).
- `breakingWeggel({ h, T, slope, g? }) -> { Hb, index, a, b }` — Hb/hb = b − a Hb/(gT²) [m, -, -, -]. Weggel (1972); SPM (1984).
- `breakingGoda({ h, T, slope, A?, slopeFactor?, g? }) -> { Hb, index, L0 }` — Goda breaker index, A = 0.17, K = 15 [m, -, m]. Goda (1975, 2000).
- `godaSurfZoneHeights({ H0, T, h, slope, Ks?, g? }) -> { Hs, Hmax, Ks, beta0, beta1, betaMax, surfZone }` — random-wave H1/3 and Hmax in the surf zone [m, m, -…]. Goda (1975, 2000).
- `battjesStiveGamma(s0) -> γ` — γ = 0.5 + 0.4 tanh(33 s0) [-]. Battjes & Stive (1985).
- `battjesJanssenQb(Hrms/Hm) -> Qb` — fraction of broken waves [-]. Battjes & Janssen (1978).
- `energyFluxBalance({ x, depth, Hrms0, Tp, theta0Deg?, gamma?, alpha?, fw?, rho?, eta0?, setup?, hMin?, g? }) -> { x, Hrms, eta, h, thetaDeg, Qb, Db, Df, energyFlux, Sxx, setupMax, nWet }` — cross-shore energy-flux balance with bore dissipation, bottom friction and radiation-stress set-up [m, m, m, m, deg, -, W/m², W/m², W/m, N/m, m, -]. Battjes & Janssen (1978); Longuet-Higgins & Stewart (1964); Thornton & Guza (1983).
- `fresnel(x) -> { C, S }` — Fresnel integrals [-]. Abramowitz & Stegun (1964) §7.3.
- `diffractionSemiInfinite({ x, y, L, theta0Deg? }) -> { Kd, re, im, phase }` — Sommerfeld solution behind a semi-infinite breakwater [-, -, -, rad]. Penney & Price (1952).
- `diffractionGap({ x, y, L, B, theta0Deg? }) -> { Kd, re, im }` — breakwater gap by superposition of two tips, B ≳ L [-]. Penney & Price (1952).
- `reflectionSeeligAhrens({ xi, a?, b? }) -> Kr` — Kr = a ξ²/(b + ξ²), defaults a = 1.0, b = 5.5 [-]. Seelig & Ahrens (1981).
- `waveCurrent({ T, h, U, H0?, g? }) -> { blocked, k, L, sigma, Tr, Cr, Cgr, Ca, Cga, heightRatio, H, Ublock }` — Doppler-shifted dispersion, wave-action height change and blocking [rad/m, m, rad/s, s, m/s…, -, m, m/s]. Longuet-Higgins & Stewart (1961); Peregrine (1976).

## runup.js

- `iribarren({ tanAlpha, H, T, g? }) -> { xi, L0, steepness, breakerType }` — surf-similarity number [-, m, -]. Battjes (1974).
- `huntRunup({ H, T, tanAlpha, g? }) -> { R, xi, valid }` — R/H = ξ for breaking waves [m, -]. Hunt (1959).
- `ROUGHNESS_FACTORS` — γf by armour type [-]. EurOtop (2018) ch. 6.
- `obliqueFactor({ betaDeg, kind? }) -> γβ` — 1 − 0.0022|β| (run-up), 1 − 0.0033|β| (overtopping), 1 − 0.0063|β| (rubble) [-]. EurOtop (2018).
- `bermFactor({ B, Lberm, db, Hm0, Ru2 }) -> γb` — berm influence factor, 0.6–1 [-]. EurOtop (2018) ch. 5.
- `roughnessFactorSurging({ gammaF, xi }) -> γf,surging` — linear increase to 1 between ξ = 1.8 and 10 [-]. EurOtop (2018) §6.2.
- `runupEurOtop({ Hm0, Tm10, tanAlpha, gammaB?, gammaF?, gammaBeta?, design?, surgingRoughness?, maxRatio?, g? }) -> { Ru2, ratio, xi, regime }` — Ru2% mean (1.65; 1.0(4 − 1.5/√ξ)) or design (1.75; 1.07) [m, -, -]. EurOtop (2018) Eq. 5.1, 5.2.
- `overtoppingSlope({ Hm0, Tm10, tanAlpha, Rc, gammaB?, gammaF?, gammaBeta?, gammaV?, gammaStar?, design?, g? }) -> { q, qBreaking, qMax, qDimensionless, xi, regime }` — mean discharge on slopes [m³/s/m]. EurOtop (2018) Eq. 5.10–5.13.
- `overtoppingSteepSlope({ Hm0, cotAlpha, Rc, gammaF?, gammaBeta?, g? }) -> { q, a, b }` — very steep slopes, 0 ≤ cot α ≤ 2 [m³/s/m, -, -]. EurOtop (2018) ch. 5.
- `overtoppingVertical({ Hm0, Rc, Tm10?, h?, foreshore?, design?, g? }) -> { q, qDimensionless, regime, impulsiveness }` — plain vertical walls: no foreshore, non-impulsive, impulsive [m³/s/m, -, -, -]. EurOtop (2018) ch. 7.
- `overtoppingRubbleMound({ Hm0, Tm10, tanAlpha, Rc, gammaF, gammaBeta?, Gc?, design?, g? }) -> { q, qNoCrest, Cr, xi, regime }` — armoured slope with crest reduction Cr = 3.06 exp(−1.5 Gc/Hm0) [m³/s/m, m³/s/m, -, -]. EurOtop (2018) §6.3.
- `TOLERABLE_OVERTOPPING` — tolerable q [m³/s/m] and Vmax [m³/m] for people, vehicles, structures, property. EurOtop (2018) ch. 3.
- `classifyOvertopping({ q, Vmax?, limits? }) -> [{ id, group, description, qLimit, VmaxLimit, utilisationQ, utilisationV, ok }]` — utilisation against each limit [-].
- `overtoppingVolumes({ q, Hm0, Tm10, Tm, duration, Rc?, Ru2?, Pov?, b?, g? }) -> { Pov, Nw, Now, a, b, Vmax, Vmean, exceedance(V), volumeAt(p) }` — Weibull distribution of individual volumes and maximum volume [-, -, -, m³/m, -, m³/m, m³/m]. EurOtop (2018) ch. 5.
- `seaLevelRiseSensitivity({ rises, Rc, Hm0, h?, depthLimited?, evaluate }) -> [{ rise, Rc, Hm0, h, q, factor, crestRaise }]` — overtopping versus sea-level rise and the crest raise that restores present q [m, m, m, m, m³/s/m, -, m].

## armour.js

- `relativeDensity({ rhoS, rhoW? }) -> Δ` — ρs/ρw − 1 [-].
- `HUDSON_KD` — no-damage KD (breaking / non-breaking, trunk) [-]. SPM (1984) Table 7-8.
- `hudson({ H, KD, cotAlpha, rhoS, rhoW?, g? }) -> { M50, W50, Dn50, Ns, delta }` — M50 = ρs H³/(KD Δ³ cot α) [kg, N, m, -, -]. Hudson (1959); SPM (1984).
- `vanDerMeerRock({ Hs, Tm, cotAlpha, P, S, N, delta, rhoS?, g? }) -> { Ns, Dn50, M50, xi, xiCr, regime }` — plunging/surging rock formulae [-, m, kg, -, -]. Van der Meer (1988); Rock Manual (2007) §5.2.2.2.
- `vanDerMeerDamage({ Hs, Tm, cotAlpha, P, N, delta, Dn50, g? }) -> { S, Ns, regime }` — damage level of a given armour [-]. Van der Meer (1988).
- `vanDerMeerShallow({ Hs, H2, Tm10, cotAlpha, P, S, N, delta, cpl?, cs?, rhoS?, g? }) -> { Ns, Dn50, M50, xi, xiCr, regime }` — shallow-water modification, cpl = 8.4, cs = 1.3 [-, m, kg, -, -]. Van Gent et al. (2003).
- `concreteUnitVanDerMeer({ unit, Hs, Tm, Nod, N, delta, rhoS?, g? }) -> { Ns, Dn, M, som }` — cubes and tetrapods [-, m, kg, -]. Van der Meer (1988).
- `stabilityNumberSize({ Hs, Ns, delta, rhoS? }) -> { Dn, M, volume }` — unit size from a design stability number [m, kg, m³].
- `toeStabilityVanDerMeer({ Hs, ht, h, delta, Nod?, rhoS? }) -> { Ns, Dn50, M50, valid }` — Hs/(Δ Dn50) = (2 + 6.2 (ht/h)^2.7) Nod^0.15 [-, m, kg]. Van der Meer (1998).
- `massToDn({ M50, rhoS, shape? }) -> { Dn50, D50 }` — nominal and sieve diameter [m]. Rock Manual (2007) §3.4.
- `dnToMass({ Dn50, rhoS }) -> M50` — [kg].
- `armourLayer({ Dn50, n?, kt?, porosity?, rhoS? }) -> { thickness, unitsPerArea, massPerArea }` — t = n kt Dn50, N = n kt (1 − nv)/Dn50² [m, 1/m², kg/m²]. SPM (1984) Eq. 7-121, 7-122.
- `filterRules({ D15f, D50f?, D85b, D15b, D50b?, retentionLimit?, permeabilityLimit?, uniformityLimit? }) -> { retention, permeability, uniformity, ok }` — geometrically closed filter criteria [-]. Terzaghi & Peck (1948).
- `underlayerMass({ M50armour, ratio? }) -> M50` — underlayer mass, ratio 10–15 [kg]. SPM (1984); Rock Manual (2007).
- `rubbleMoundQuantities({ crestWidth, height, cotFront, cotRear, layers, core? }) -> { totalArea, baseWidth, layers: [{ name, area, solidVolume, mass, crestWidth, height }] }` — per-layer volumes of a trapezoidal section [m³/m, m, m³/m, kg/m]. Plane geometry.

## caisson.js

- `takahashiImpulsive({ H, d, h, BM, L }) -> { alphaI, alphaI0, alphaI1 }` — impulsive pressure coefficient [-]. Takahashi, Tanimoto & Shimosako (1994).
- `godaPressure({ Hmax, T, h, d, hPrime, hc, hb, B, betaDeg?, BM?, lambda1?, lambda2?, lambda3?, rho?, g? }) -> { etaStar, L, alpha1, alpha2, alpha3, alphaI, alphaStar, p1, p2, p3, p4, pu, P, Mp, armP, U, Mu, hcStar }` — Goda wave pressures, horizontal force, uplift and moments [m, m, -…, Pa…, N/m, N·m/m, m, N/m, N·m/m, m]. Goda (1974, 2000) §4.2.
- `caissonStability({ P, Mp, U, Mu, W, t, B, mu? }) -> { sfSliding, sfOverturning, We, te, eccentricity, pHeel, pToe, contactWidth, fullContact }` — sliding, overturning and heel bearing pressure [-, -, N/m, m, m, Pa, Pa, m]. Goda (2000) §4.3.
- `sainflou({ H, T, h, hc?, rho?, g? }) -> { L, h0, p1, p2, crest: { force, moment, arm, pTop }, trough: { force, moment, arm, pSurface } }` — standing-wave pressures [m, m, Pa, Pa, N/m, N·m/m, m, Pa]. Sainflou (1928); SPM (1984).
- `gravityWallStability({ W, xW, B, hHigh, hLow, mu?, upliftFactor?, Hextra?, Mextra?, rho?, g? }) -> { Fh, Fl, U, armU, N, H, sfSliding, sfOverturning, sfFlotation, resultantFromToe, eccentricity, pMax, pMin, contactWidth, middleThird }` — hydrostatic thrust, uplift, sliding, overturning and base pressure [N/m…, m, -, Pa, m]. Statics.
- `flotationCheck({ weight, displacedVolume, resisting?, upliftExtra?, rho?, g? }) -> { buoyancy, sf, netDownward, floats, ballastForSf(sf) }` — buoyancy safety factor [N, -, N]. Archimedes.

## sediment.js

- `dimensionlessGrainSize({ d, rhoS?, rho?, nu?, g? }) -> D*` — [-]. van Rijn (1984); Soulsby (1997) Eq. 75.
- `settlingVelocitySoulsby({ d, rhoS?, rho?, nu?, g? }) -> ws` — [m/s]. Soulsby (1997) Eq. 102.
- `settlingVelocityVanRijn({ d, rhoS?, rho?, nu?, g? }) -> ws` — three size ranges [m/s]. van Rijn (1984b).
- `criticalShields(Dstar) -> θcr` — 0.30/(1 + 1.2 D*) + 0.055 [1 − exp(−0.020 D*)] [-]. Soulsby & Whitehouse (1997).
- `criticalShearStress({ d, rhoS?, rho?, nu?, g? }) -> { Dstar, thetaCr, tauCr, uStarCr }` — [-, -, Pa, m/s].
- `shieldsParameter({ tau, d, rhoS?, rho?, g? }) -> θ` — [-]. Shields (1936).
- `mobilityCheck({ tau, d, rhoS?, rho?, nu?, g? }) -> { theta, thetaCr, ratio, mobile, Dstar }` — threshold-of-motion check [-].
- `grainRoughness(d50) -> { ks, z0 }` — ks = 2.5 d50, z0 = d50/12 [m]. Soulsby (1997).
- `bedformRoughness({ height, length, ar? }) -> { z0, ks }` — z0 = ar Δ²/λ [m]. Soulsby (1997).
- `currentRipples(d50) -> { length, height, valid }` — λ = 1000 d50, Δ = λ/7 [m]. Soulsby (1997).
- `bedShearCurrent({ U, h, method?, z0?, n?, ks?, rho?, g? }) -> { tau, uStar, Cd, chezy }` — log law, Manning or Colebrook–White [Pa, m/s, -, m^½/s]. Soulsby (1997) §3.3.
- `waveOrbital({ H, T, h, g? }) -> { Uw, A, k }` — near-bed orbital velocity and excursion [m/s, m, rad/m]. Linear theory.
- `waveFrictionFactor({ A, method?, z0?, ks? }) -> fw` — Soulsby fw = 1.39 (A/z0)^−0.52 or Swart (1974) [-].
- `bedShearWave({ Uw, T, method?, z0?, ks?, rho? }) -> { tau, fw, A }` — τw = ½ ρ fw Uw² [Pa, -, m].
- `combinedShearSoulsby({ tauC, tauW, phiDeg? }) -> { tauMean, tauMax }` — DATA2 wave–current stresses [Pa]. Soulsby (1995; 1997 Eq. 69–70).
- `meyerPeterMuller({ theta, d, rhoS?, rho?, thetaCr?, coef?, g? }) -> { qb, phi }` — Φ = 8 (θ − 0.047)^1.5 [m²/s, -]. Meyer-Peter & Müller (1948).
- `engelundHansen({ U, tau, d50, rhoS?, rho?, g? }) -> { qt, phi, theta }` — total load, f Φ = 0.1 θ^2.5 [m²/s, -, -]. Engelund & Hansen (1967).
- `vanRijnBedLoad({ U, h, d50, d90?, rhoS?, rho?, nu?, g? }) -> { qb, T, Dstar, uStarGrain, uStarCr }` — 0.053 T^2.1/D*^0.3 form [m²/s, -, -, m/s, m/s]. van Rijn (1984a).
- `vanRijnSuspendedLoad({ U, h, d50, d90?, ks?, a?, sigmaS?, rhoS?, rho?, nu?, g? }) -> { qs, ca, a, F, Z, Zprime, ds, ws, uStar, T }` — reference-concentration method [m²/s, -, m, -, -, -, m, m/s, m/s, -]. van Rijn (1984b).
- `vanRijnSimplified({ U, h, d50, d90?, rhoS?, rho?, nu?, g? }) -> { qb, qs, qt, Ucr }` — simplified current-alone formulae [m²/s, m/s]. van Rijn (1984); Soulsby (1997) Eq. 133.
- `soulsbyVanRijn({ U, Urms, h, d50, d90?, bedSlope?, z0?, rhoS?, rho?, nu?, g? }) -> { qt, qb, qs, Ucr, Cd }` — total load under waves plus current [m²/s, m/s, -]. Soulsby (1997) Eq. 136.
- `rouseProfile({ z, h, a, ca, ws, uStar }) -> { c, rouse }` — suspended concentration profile [unit of ca, -]. Rouse (1937).
- `isbashStoneSize({ U, delta, C?, g? }) -> D` — D = U²/(2 g Δ C²), C = 0.86 or 1.20 [m]. Isbash (1936).
- `shieldsRockSize({ U, h, delta, psiCr?, ksFactor?, slopeDeg?, reposeDeg?, g? }) -> { d, chezy, tau, ksl }` — Shields/Chézy rock size under a current [m, m^½/s, m²/s², -]. Rock Manual (2007) §5.2.

## swe.js

- `solveSwe1d({ dx, zb, h0 | eta0, u0?, tEnd, x0?, cfl?, manning?, limiter?, hDry?, left?, right?, outputTimes?, landSide?, hShore?, g?, maxSteps? }) -> { x, zb, t, h, u, eta, snapshots, shoreline: { t, x, z, maxRunup, maxRunupTime }, envelope: { hMax, etaMax, uMax, inundationLimit }, volume: { initial, final, boundaryNet, clipped, error, relativeError }, steps, dtMin, dtMax }` — 1-D nonlinear shallow-water finite-volume solver (HLL, MUSCL, hydrostatic reconstruction, wetting/drying, Manning, SSP-RK2); boundaries `wall | transmissive | periodic | stage | discharge | wave` [m, m/s, s, m² per unit width]. Toro (2001); Audusse et al. (2004).
- `damBreakAnalytic({ hL, hR?, x, t, g? }) -> { h, u, hm, um, shockSpeed }` — Stoker (wet bed) and Ritter (dry bed) dam-break solutions [m, m/s, m, m/s, m/s]. Stoker (1957); Ritter (1892).

## morpho.js

- `exnerEvolve({ dx, zb0, transport(zb, t), tEnd, porosity?, scheme?, limiter?, cfl?, dt?, dtMax?, left?, right?, maxSteps? }) -> { zb, t, steps, courantMax, mass: { initial, final, boundaryNet, error, relativeError } }` — Exner bed evolution, upwind or flux-limited Lax–Wendroff, with sediment mass diagnostic [m, s, m²]. Exner (1925); LeVeque (2002).
- `deanParameter(ws) -> A` — A = 0.067 ws^0.44 with ws in cm/s (argument in m/s) [m^⅓]. Dean (1987).
- `deanProfile({ y, A }) -> { h, widthToDepth(h) }` — h = A y^⅔ [m]. Dean (1977).
- `bruunRule({ S, hClosure, B, L?, A? }) -> { retreat, L, slope, ratio }` — R = S L/(B + h*) [m, m, -, -]. Bruun (1962).
- `depthOfClosure({ He, Te, method?, g? }) -> hc` — Hallermeier 2.28 He − 68.5 He²/(g Te²) or Birkemeier 1.75/57.9 [m]. Hallermeier (1978); Birkemeier (1985).
- `depthOfClosureOuter({ HsMean, sigmaHs, TsMean, d50, g? }) -> hi` — (Hs − 0.3 σ) Ts √(g/(5000 d50)) [m]. Hallermeier (1981).
