# `src/engine/metocean` — exports

SI units throughout; `z` is elevation above still-water level (m, up), `d` depth (m), angles in rad, spectra
`S(ω)` in m²·s/rad, wind spectra in (m/s)²/Hz. Accelerations `ax`, `az` are local (∂/∂t); `axTotal`, `azTotal`
are particle accelerations. Tests: `tests/engine/metocean.<file>.test.js`.

## waves.js

| export | signature → result (units) | reference |
| --- | --- | --- |
| `G`, `RHO_SEA` | 9.80665 m/s², 1025 kg/m³ | — |
| `wavenumber` | `({ T s, d m, U = 0 m/s, g }) → k rad/m` — exact linear dispersion, Doppler shift, throws when blocked | Dean & Dalrymple (1991) §3.4; DNV-RP-C205 §3.2.2 |
| `waveProperties` | `({ T, d, H?, U, rho, g }) → { k, L m, omega, omegaRel rad/s, c, cRel, cg, cgRel m/s, n, kd, L0 m, E J/m², P W/m, steepness, ursell, regime }` | Dean & Dalrymple ch. 3–4 |
| `wavelength`, `celerity`, `groupVelocity` | `(T s, d m, g?) → m, m/s, m/s` | linear theory |
| `energyDensity` | `(H m, rho?, g?) → J/m²` = ρgH²/8 | linear theory |
| `energyFlux` | `({ H, T, d, rho, g }) → W/m` = E·cg | linear theory |
| `depthRegime` | `(d/L) → 'deep' \| 'intermediate' \| 'shallow'` (limits 1/2, 1/20) | SPM (1984) Table 2-1 |
| `ursellNumber`, `steepness` | `(H, L, d) → HL²/d³`; `(H, L) → H/L` | Ursell (1953) |
| `airyKinematics` | `({ H, T, d, x, z, t, phase, U, stretching: 'none'\|'wheeler'\|'linear'\|'vertical', clipAboveSurface, rho, g, k? }) → { eta m, u, w m/s, ax, az m/s², pDyn Pa, k, omega, theta, wet }` | DNV-RP-C205 §3.2.2, §3.5.2; Wheeler (1970) |
| `stokes2Kinematics` | `({ H, T, d, x, z, t, rho, g }) → { eta, eta1, eta2 m, u, w, ax, az, pDyn, k, omega, stokesDrift m/s, secondaryCrest, wet }` | Dean & Dalrymple §11.2; SPM Eqs. 2-45…2-52 |
| `stokes5Coefficients` | `(kd) → { S, C0, C2, C4, A[i][j], B[i][j] }` (-) | Fenton (1985) Table 1, C4 per Fenton (1990) |
| `stokes5Wave` | `({ H, d, T \| L, currentEulerian = 0 \| currentStokes, g }) → { k, L, T, c, omega, eps, uBar, currentEulerian, currentStokes m/s, Q m²/s, R m²/s², etaAmp[] m, velAmp[] m/s, crest, trough m, ursell, valid }` | Fenton (1985) |
| `stokes5Kinematics` | `(sol, { x, z, t, rho }) → { eta, u, w, ax, az, axTotal, azTotal, pDyn, wet }` | Fenton (1985) Eqs. 12–16 |
| `ellipticKE` | `(m, m1 = 1 − m) → { K, E }` | Abramowitz & Stegun §17.6 (AGM) |
| `jacobiElliptic` | `(u, m, m1?) → { sn, cn, dn }` | Abramowitz & Stegun §16.4 |
| `solitaryWave` | `({ H, d, x, z, t, rho, g }) → { eta, u, w, ax, az, pDyn, c m/s, beta 1/m, volume m², energy J/m, wet }` | Boussinesq (1872); Munk (1949) |
| `cnoidalWave` | `({ H, d, T \| L, g }) → { m, m1, K, E, L, T, c, crest, trough, ursell }` (first order, KdV) | Korteweg & de Vries (1895); Wiegel (1960) |
| `cnoidalKinematics` | `(sol, { x, z, t, rho }) → { eta, u, w, ax, az, pDyn, wet }` | as above; mass conservation for u |
| `breakingHeightMiche` | `(L m, d m) → Hb m` = 0.142 L tanh(2πd/L) | Miche (1944) |
| `breakingHeightMcCowan` | `(d m, kappa = 0.78) → Hb m` | McCowan (1894) |
| `breakingHeightGoda` | `({ T, d, slope, A = 0.17, K = 15, g }) → Hb m` | Goda (1974); K = 11 in Goda (2010) |
| `breakingLimits` | `({ T, d, slope, g }) → { miche, mccowan, goda, governing, L }` m | as above |
| `selectWaveTheory` | `({ H, T, d, g }) → { theory, streamFunctionOrder, HgT2, dgT2, ursell, steepness, dOverL, HoverHb, regime, approximate: true }` — region boundaries approximated | API RP 2A-WSD Fig. 5.3 / ISO 19901-1 Fig. A.4 |
| `shoalingRefraction` | `({ T, d, d0 = ∞, theta0 rad, H0, g }) → { Ks, Kr, theta rad, H m, k, L, c, cg }` | Snell; Dean & Dalrymple §4.8 |

## streamfunction.js

| export | signature → result (units) | reference |
| --- | --- | --- |
| `streamFunctionWave` | `({ H, d, T \| L, order = 10 (1…30), currentEulerian = 0 \| currentStokes, heightSteps?, tol = 1e-12, g }) → { k, L, T, c, omega, currentEulerian, currentStokes m/s, B0 m/s, B[] m²/s, Q m²/s, R m²/s², etaNodes[] m, etaCoef[] m, crest, trough m, report: { converged, iterations, heightSteps, residual, kinematicResidual, dynamicResidual } }` | Rienecker & Fenton (1981); Fenton (1988); Dean (1965) |
| `streamFunctionEta` | `(sol, x m, t s) → eta m` (on the streamline ψ = −Q) | Fenton (1988) |
| `streamFunctionProfile` | `(sol, n = 73) → { x m, eta m }` over one wavelength | — |
| `streamFunctionKinematics` | `(sol, { x, z, t, rho }) → { eta, u, w, ax, az, axTotal, azTotal, p Pa, pDyn Pa, wet }` up to the instantaneous surface | Rienecker & Fenton (1981) Eqs. 6–9 |

Checks in the tests: Airy limit; both free-surface conditions < 1e-6 between collocation points; agreement with
Stokes 5th order; deep-water c²k/g against the Stokes series at ak = 0.2; near-limiting deep-water phase speed
≈ 1.093 √(g/k) near ak = 0.436 (Longuet-Higgins 1975; Cokelet 1977).

## spectra.js

| export | signature → result (units) | reference |
| --- | --- | --- |
| `piersonMoskowitz` | `(omega rad/s, Hs m, Tp s) → S m²·s/rad` | DNV-RP-C205 §3.5.5.1 |
| `jonswapGamma` | `(Hs m, Tp s) → γ` from Tp/√Hs | DNV-RP-C205 §3.5.5.5 |
| `jonswap` | `(omega, { Hs, Tp, gamma?, sigmaA = 0.07, sigmaB = 0.09 }) → S`, Aγ = 1 − 0.287 ln γ | DNV-RP-C205 §3.5.5.2 |
| `jonswapPeriodRatios` | `(gamma) → { Tm01OverTp, TzOverTp }` (1 ≤ γ ≤ 7) | DNV-RP-C205 §3.5.5.4 |
| `tmaFactor` | `(omega, d m, g?) → φ` = tanh²kd / (1 + 2kd/sinh 2kd) | Kitaigorodskii et al. (1975); Bouws et al. (1985) |
| `tma` | `(omega, { Hs, Tp, gamma?, d, g? }) → S` | Bouws et al. (1985); DNV-RP-C205 §3.5.6 |
| `bretschneider` | `(omega, { Hs, T, periodType: 'Tp'\|'Tz'\|'T1' }) → S` | Bretschneider (1959); ITTC (1978); ISSC (1964) |
| `ochiHubble` | `(omega, { Hs1, Tp1, lambda1, Hs2, Tp2, lambda2 }) → S` | Ochi & Hubble (1976) |
| `spreadingCosN` | `(theta, thetaP, n) → D 1/rad` | DNV-RP-C205 §3.5.8 |
| `spreadingCos2s` | `(theta, thetaP, s) → D 1/rad` | Longuet-Higgins et al. (1963) |
| `tabulatedMoments` | `(omega[], S[]) → { mMinus1, m0, m1, m2, m4 }` m²·(rad/s)ⁿ | DNV-RP-C205 §3.5.3 |
| `spectrumMoments` | `(fn, { wMin, wMax, n = 4000, tail = false }) → moments` (Simpson, optional ω⁻⁵ tail) | DNV-RP-C205 §3.5.3 |
| `spectralParameters` | `(moments) → { Hm0 m, Tz, Tm01, Te, Tc s, epsilon, nu }` | Cartwright & Longuet-Higgins (1956) |
| `rayleighExceedance` | `(h m, Hs m) → P` = exp(−2h²/Hs²) | Longuet-Higgins (1952) |
| `mostProbableMaxHeight` | `(Hs, N) → m` = Hs √(½ ln N) | Longuet-Higgins (1952) |
| `expectedMaxHeight` | `(Hs, N, 'asymptotic' \| 'exact') → m` | Longuet-Higgins (1952) |
| `maxWaveHeight` | `({ Hs, Tz, duration s }) → { N, mostProbable, expected m, mostProbableOverHs, expectedOverHs }` | DNV-RP-C205 §3.5.11 |
| `forristallCrest` | `({ Hs, T1, d, type: '2D'\|'3D', N?, beta3Ur2 = 0.284, g }) → { alpha, beta, S1, Ur, k1, scale m, mostProbableMax, expectedMax m }` | Forristall (2000); DNV-RP-C205 §3.5.10.4 |
| `forristallExceedance` | `(x m, params) → P(crest > x)` | Forristall (2000) |
| `nearBedVelocity` | `(fn, { d, wMin, wMax, n, g }) → { Us m/s, Tu s, M0, M2 }` | DNV-RP-F109 §3.4.3 |
| `synthesizeSeaState` | `(fn, { wMin, wMax, duration s, dt s, nFreq = 200, seed = 1, x, d, g }) → { t, eta, components: { omega, amp, phase, k } }` | DNV-RP-C205 §3.3.2 |
| `regularWaveEquivalents` | `({ Hs, Tp, Tz?, Te?, duration = 10800 s, cLow = 6.5, cHigh = 11 }) → { energy, significant, design }` (H m, T s) | energy equivalence; NORSOK N-003 period range |

## wind.js

| export | signature → result (units) | reference |
| --- | --- | --- |
| `RHO_AIR`, `KARMAN` | 1.225 kg/m³, 0.4 | — |
| `powerLawProfile` | `(z m, Uref m/s, zRef = 10, alpha = 0.12) → m/s` | DNV-RP-C205 §2.3.2 |
| `logProfile` | `(z, Uref, zRef = 10, z0 = 0.0002 m) → m/s` | DNV-RP-C205 §2.3.2 |
| `frictionVelocity` | `(U, z, z0) → u* m/s` | DNV-RP-C205 §2.3.2 |
| `equivalentShearExponent` | `(z1, z2, z0) → α` | DNV-RP-C205 §2.3.2 |
| `froyaWind` | `({ U0 m/s (1 h, 10 m), z = 10, t = 3600 s }) → { u, Uz m/s, Iu, C, gustFactor, sigma m/s }` | ISO 19901-1 §A.7.3; NORSOK N-003 |
| `froyaReferenceSpeed` | `({ u, z, t }) → U0 m/s` (inverse) | as above |
| `convertWindSpeed` | `({ u, zFrom, tFrom, zTo, tTo }) → { u, U0 }` m/s | as above |
| `gustFactor` | `({ U0, z = 10, t = 3 }) → u(z,t)/U(z)` | ISO 19901-1 Eq. A.6 |
| `froyaSpectrum` | `(f Hz, { U0, z }) → S (m/s)²/Hz` | ISO 19901-1 §A.7.4 |
| `IEC_IREF`, `IEC_VREF` | `{ 'A+': 0.18, A: 0.16, B: 0.14, C: 0.12 }`, `{ I: 50, II: 42.5, III: 37.5 }` m/s | IEC 61400-1 Table 1 |
| `iecTurbulenceScale` | `(z m) → Λ1 m` | IEC 61400-1 §6.3.1 |
| `iecNtm` | `(Vhub, Iref) → { sigma1 m/s, TI }` | IEC 61400-1 §6.3.2.3 |
| `iecEtm` | `({ Vhub, Iref, Vave }) → { sigma1, TI }` | IEC 61400-1 §6.3.3.4 |
| `kaimalSpectrum` | `(f Hz, { Vhub, sigma1, z, component: 'u'\|'v'\|'w' }) → S (m/s)²/Hz` | IEC 61400-1 Annex C |
| `iecExtremeWind` | `({ Vref, zhub, z? }) → { Ve50, Ve1, V50, V1, sigma1 }` m/s | IEC 61400-1 §6.3.3.2 |
| `iecEog` | `({ Vhub, Ve1, sigma1, D m, zhub, t?, T = 10.5 s }) → { Vgust, T, V, Vmax, Vmin }` | IEC 61400-1 §6.3.3.3 |
| `API_SHAPE_COEFFICIENTS` | beams 1.5, building sides 1.5, cylindrical sections 0.5, overall platform projected area 1.0 | API RP 2A-WSD wind force shape coefficients |
| `airDensity` | `({ temperatureC = 15, pressure = 101325 Pa, relativeHumidity = 0 }) → kg/m³` | ideal gas; Tetens (1930) |
| `windForce` | `({ U, A m², Cs = 1, Ch = 1, rho, alpha = π/2 }) → { F N, q Pa }` | API RP 2A-WSD; DNV-RP-C205 §5.3 |
| `heightCoefficient` | `(z, zRef = 10, alpha = 0.12) → Ch` = (z/zRef)^(2α) | derived from the profile |
| `weibullWind` | `({ A m/s, k, U?, rho }) → { pdf, cdf, mean, std, meanCube, powerDensity W/m², mode, median }` | IEC 61400-1 §6.3.2.1 |
| `weibullFromMeanStd` | `(mean, std) → { A, k }` | method of moments |
| `weibullFitMle` | `(speeds[]) → { A, k, n, nExcluded, logLik }` | Cohen (1965) |

## current.js

| export | signature → result (units) | reference |
| --- | --- | --- |
| `SEABED_ROUGHNESS` | z0 (m) by soil: silt 5e-6 … boulder 4e-2 | DNV-RP-F109 Table 3-1 |
| `currentPowerLaw` | `(z, { v0 m/s, d, alpha = 1/7 }) → m/s` | DNV-RP-C205 §4.1.4.2 |
| `currentLinear` | `(z, { v0, d0 = 50 m, d }) → m/s` (wind-generated) | DNV-RP-C205 §4.1.4.2 |
| `currentBilinear` | `(z, { vSurface, vBreak, vBed, zBreak, d }) → m/s` | piecewise-linear design profile |
| `currentLogNearBed` | `(zBed m, { vRef, zRef, z0 }) → m/s` | DNV-RP-F109 Eq. 3.2 |
| `meanCurrentOverHeight` | `({ vRef, zRef, z0, D, theta = π/2 }) → m/s` | DNV-RP-F109 Eq. 3.3 |
| `currentTabulated` | `(z, { zs[], vs[] }) → m/s` | site data |
| `windDrivenSurfaceCurrent` | `(U1h10 m/s, factor = 0.03) → m/s` (percentage rule, 0.015–0.03) | DNV-RP-C205 §4.1.4.4 |
| `currentTideAndWind` | `(z, { vTide, vWind, d, alpha, d0 }) → m/s` | DNV-RP-C205 §4.1.4 |
| `stretchedElevation` | `({ z, eta, d, method: 'linear'\|'nonlinear'\|'none', k? }) → z' m` | API RP 2A-WSD current stretching; ISO 19901-1 §A.9 |
| `stretchedCurrent` | `(profileFn, { z, eta, d, method, k }) → m/s` | as above |
| `depthAveragedCurrent` | `(profileFn, { d, n = 400 }) → m/s` | definition |
| `nearBedWaveCurrent` | `({ H, T, d, Uc, angle rad, zBed, g }) → { Uw m/s, Aw m, k, Tr s, Umax, Umin m/s, ratio }` | linear theory with Doppler shift |

## windows.js

Series are objects of equal-length arrays at a constant interval `dtHours`; limits are
`[{ key, max?, min?, alpha? }]` (alpha scales `max`). Durations in hours.

| export | signature → result (units) | reference |
| --- | --- | --- |
| `operationReferencePeriod` | `({ plannedHours, contingencyHours, maxReferenceHours = 96, maxPlannedHours = 72 }) → { TR, TPOP, TC h, weatherRestricted }` | DNV-ST-N001 |
| `forecastLimit` | `(designLimit, alpha) → alpha·designLimit` (alpha from the DNV table, supplied by the caller) | DNV-ST-N001 alpha factor |
| `workabilityMask` | `(series, limits) → Uint8Array` | DNV-ST-N001 |
| `findWindows` | `(mask, dtHours, minHours = 0) → [{ startIndex, endIndex, startHours, durationHours }]` | — |
| `operability` | `({ series, limits, dtHours, requiredHours?, timeMs? }) → { workableFraction, startableFraction, meanWaitingHours, maxWaitingHours, p90WaitingHours, windows, monthly[12], seasonal[4], mask, … }` | DNV-ST-N001; DNV-RP-N103 |
| `persistenceTable` | `({ values, dtHours, thresholds[], durations[] h, sense: 'below'\|'above' }) → [{ threshold, occurrence, nEvents, meanDurationHours, maxDurationHours, eventExceedance[], timeFraction[] }]` | Graham (1982) |
| `markovStationary` | `(P, { tol, maxIter }) → π[]` | Norris (1997) |
| `markovFit` | `(values[], edges[]) → { edges, nStates, counts, P, occupancy, stationary }` | Anastasiou & Tsekos (1996) |
| `markovSimulate` | `(model, { nSteps, seed = 1, startState?, sample: 'uniform'\|'mid' }) → { states, values }` | as above |
| `campaignMonteCarlo` | `({ series, dtHours, operations: [{ name, limits, netHours, contingencyHours }], nTrials = 1000, seed = 1, timeMs?, startMonths? }) → { completionHours: { p10, p50, p90, mean, min, max }, downtimeHours: { p10, p50, p90, mean }, operations[], samples }` | DNV-ST-N001; DNV-RP-N103 |

## Not implemented

* Torsethaugen double-peak spectrum (the simplified DNV-RP-C205 form has many fitted constants; omitted rather
  than reproduced from memory). Use `ochiHubble` for two-peak seas.
* Tabulated height coefficients of mobile-unit codes (use `heightCoefficient` or a profile).
* DNV-ST-N001 alpha-factor table values (pass the tabulated value as `alpha`).
* Higher-order cnoidal theory and limiting (breaking) stream-function waves above about 97 % of the maximum height.
