# `src/engine/lines` — mooring lines, cable mechanics, interference and berth mooring

SI units throughout (m, kg, s, N, rad; arguments ending in `Deg` are degrees). One line per export:
signature — what it returns (units) — reference. Property estimators and built-in coefficient curves are
indicative and every constant can be overridden.

## catenary.js — single uniform line in a vertical plane

Anchor at (0, 0), fairlead at (X, Z); H, V = horizontal / vertical tension components at the fairlead.

- `inextensibleCatenary({ H, w, z }) -> { s, x, T, V, angle, a }` — suspended length [m], horizontal projection [m], tension [N], angle [rad] at height z above the horizontal-tangent point — Faltinsen (1990) Eqs. 8.14–8.20.
- `inextensibleCatenaryShape({ H, w, x }) -> { z, s, T, angle }` — shape functions z = a(cosh(x/a) − 1) etc. [m, m, N, rad] — Faltinsen (1990) Eq. 8.16–8.18.
- `inextensibleHorizontalStiffness({ H, w, z }) -> number` — dH/dX of a grounded inextensible line [N/m] — Faltinsen (1990) Eq. 8.26.
- `catenaryInitialGuess({ X, Z, L, w, EA }) -> { H, V }` — starting values for the Newton solve [N] — Peyrot & Goulois (1979) via Jonkman (2007) Eq. 2-39/2-40.
- `solveShooting(fwd, { X, Z, H0, V0, scaleL, scaleF, Hmax, checkSlack, tol }) -> { H, V, slack, iterations, residual }` — damped Newton in (ln H, V) with nested-bisection fallback and slack detection for any forward function [N, m] — numerical method (shared with multiseg.js).
- `solveCatenary({ X, Z, L, w, EA, CB = 0, seabed = true, slope = 0, nProfile, guess }) -> { H, V, T, angle, anchorH, anchorV, anchorT, anchorAngle, upliftAngle, groundedLength, suspendedLength, touchdownX, touchdownZ, stretchedLength, stiffness[2][2], slack, profile }` — elastic catenary with frictional horizontal/sloping seabed, taut, vertical, slack, buoyant (w < 0) and weightless cases [N, rad, m, N/m] — Irvine (1981); Jonkman (2007) NREL/TP-500-41958 Eqs. 2-35…2-38.
- `catenaryProfile(sol, n) -> { s, x, z, T }` — coordinates and tension along the unstretched arc length [m, m, m, N] — Jonkman (2007) Eqs. 2-41…2-44.
- `lengthForFairleadTension({ X, Z, w, EA, CB, seabed, slope, T }) -> { L, solution }` — unstretched length giving a target fairlead tension [m] — inverse of solveCatenary (Brent).
- `lengthForFairleadAngle({ X, Z, w, EA, …, angle }) -> { L, solution }` — length giving a target fairlead angle from horizontal [m] — inverse of solveCatenary.
- `spanForFairleadTension({ Z, L, w, EA, …, T }) -> { X, solution }` — anchor radius giving a target pretension [m] — inverse of solveCatenary.
- `minLengthNoUplift({ Z, w, T, EA = Infinity }) -> { L, H, V, X }` — minimum length for zero anchor uplift at a design tension [m, N, N, m] — Faltinsen (1990) Eq. 8.20; elastic form from Jonkman (2007) Eq. 2-38.

## multiseg.js — multi-segment line with clump weights / buoys

- `compileLine({ segments: [{ L, w, EA, CB, name }], points: [{ junction, weight }], seabed = true, CB }) -> line` — validated, pre-processed line (anchor → fairlead) for repeated solves [m, N/m, N] — data structure.
- `solveMultiSegment({ X, Z, segments, points, seabed, CB, line, guess, nProfile, stiffness }) -> { H, V, T, angle, anchorH, anchorV, anchorT, upliftAngle, groundedLength, suspendedLength, touchdownX, stretchedLength, sagClearance, stiffness[2][2], segments[], junctions[], slack, profile }` — shooting on (H, V) with segment-wise elastic catenaries, touchdown in any segment, buoyant sections (lazy wave), point loads at junctions, friction on the grounded part [N, rad, m, N/m] — Irvine (1981); Jonkman (2007) section 2.5.

## mooringsystem.js — quasi-static mooring system

System = `{ lines: [{ anchor: [x,y,z], fairlead: [x,y,z] (body), segments, points, seabed, CB, MBL, turret }], turretCenter }`; pose = `[x, y, z, roll, pitch, yaw]`.

- `rotationMatrix(roll, pitch, yaw) -> number[3][3]` — R = Rz·Ry·Rx [-] — rigid-body kinematics.
- `lineForces(system, pose, { exclude, nProfile, warm }) -> { F[6], lines[], maxTension, maxTensionLine }` — mooring load on the floater and per-line results [N, N·m] — API RP 2SK quasi-static analysis.
- `equilibrium(system, { Fext, dofs = [0,1,5], pose0, hydroK, hydroRef, exclude, tol, maxIter }) -> { pose, offset, offsetDirection, residual, F, lines, maxTension, … }` — static equilibrium by damped Newton [m, rad, N] — API RP 2SK (mean offset).
- `restoringCurve(system, { direction, offsets, pose0, exclude }) -> { offset, restoring, transverse, yawMoment, maxTension, lineTensions, maxUplift, minGrounded }` — force–offset curve [m, N, N·m, rad] — Faltinsen (1990) section 8.2.
- `stiffnessMatrix(system, { pose, dx, dtheta, exclude }) -> number[6][6]` — linearised mooring stiffness by central differences [N/m, N/rad, N·m/rad] — Faltinsen (1990) section 8.2.
- `setPretension(system, { target, mode = 'length' | 'anchorRadius', segment, pose }) -> { system, tensions, adjustments }` — adjusts a segment length or the anchor radius to reach the target pretension [N, m] — inverse solve (Brent).
- `brokenLineAnalysis(system, { Fext, dofs, …, overshootReduction = 1, lines }) -> { intact, damaged[] }` — intact and one-line-broken equilibria with transient overshoot by the work–energy method [m, N] — API RP 2SK damaged/transient condition; energy method from first principles.
- `API_2SK_SAFETY_FACTORS` — { intact: { dynamic: 1.67, quasi-static: 2.0 }, damaged: { dynamic: 1.25, quasi-static: 1.43 } } — API RP 2SK (3rd ed.) tension limits.
- `checkApi2sk({ tension, MBL, condition, analysis, requiredSF }) -> { safetyFactor, required, uc, ok }` — tension safety-factor check [-] — API RP 2SK (3rd ed.).
- `DNV_E301_FACTORS` — γmean/γdyn by limit state (ULS/ALS), consequence class (1/2) and analysis type — DNV-OS-E301 (2010–2018) Ch.2 Sec.2.
- `checkDnvE301({ Tmean, Tdyn, MBL, limitState, consequenceClass, analysis, strengthFactor = 0.95, gammaMean, gammaDyn }) -> { Td, Sc, uc, ok, … }` — partial-factor utilisation u = (Tmean·γmean + Tdyn·γdyn)/Sc [N, -] — DNV-OS-E301 Ch.2 Sec.2.
- `checkAnchorUplift({ upliftAngle, allowable = 0, anchorV }) -> { ok, … }` — anchor uplift angle check [rad] — API RP 2SK anchor requirements.
- `checkGroundedLength({ groundedLength, required }) -> { margin, ok }` — grounded length check [m] — design practice.
- `checkOffset({ offset, depth, limitFraction }) -> { offsetFraction, uc, ok }` — offset as a fraction of water depth [-] — project limit (riser/umbilical).
- `responseStatistics({ sigma, N }) -> { sig, max }` — significant 2σ and maximum σ√(2 ln N) — API RP 2SK response statistics.
- `combineApi2sk({ mean, lfSig, lfMax, wfSig, wfMax }) -> { max, case1, case2, governing }` — max(mean + LFmax + WFsig, mean + LFsig + WFmax) — API RP 2SK (3rd ed.).
- `quasiDynamicTension(system, { meanPose, direction, lf: { sig, max }, wf: { sig, max }, daf = 1, exclude }) -> { maxOffset, lines[], maxTension }` — line tensions at the two API motion combinations with DAF on the WF increment [N] — API RP 2SK.
- `tendonStiffness(tendons) -> { axial[], surge, sway, heave, roll, pitch, yaw }` — TLP tendon-set stiffness [N/m, N·m/rad] — Faltinsen (1990) ch. 8; Chakrabarti (2005) ch. 8.
- `tendonOffset({ tendons, offset, waterplaneStiffness }) -> { setDown, restoring, tensions, angles, maxTension, minTension }` — set-down and restoring force vs offset [m, N, rad] — TLP statics.
- `tendonTensions({ tendons, Fz, Mx, My, hydroHeave, hydroRoll, hydroPitch }) -> { heave, roll, pitch, tensions, maxTension, minTension, slack }` — tendon tensions under vertical load and overturning moments [m, rad, N] — rigid-body statics; API RP 2T tension components.
- `naturalPeriod({ mass, addedMass, stiffness }) -> number` — 2π√((M + A)/K) [s] — SDOF dynamics.
- `spreadLayout({ nLines, anchorRadius, fairleadRadius, fairleadZ, depth, heading, segments, points, MBL, CB }) -> system` — symmetric spread layout.
- `clusterLayout({ nGroups, linesPerGroup, groupSpread, … }) -> system` — grouped spread layout.
- `turretLayout({ turretCenter, … }) -> system` — turret layout (lines flagged `turret`, no yaw coupling through the bearing).
- `calmHawser({ Fmean, hawserStiffness, buoyStiffness, mass, addedMass, surgeVelocity, MBL }) -> { Tmean, Tmax, dT, hawserExtension, buoyOffset, seriesStiffness, surgePeriod, utilisation }` — simple CALM hawser series-spring / energy model [N, m, s] — first principles (screening).
- `segmentDistance([p1, q1], [p2, q2]) -> { distance, s, t }` — closest distance of two 3-D segments [m] — Ericson (2005) section 5.1.9.
- `polylineMinDistance(A, B) -> { distance, iA, iB, sA, sB }` — minimum distance between polylines [m] — geometry.
- `clearances(system, { pose, nProfile, assets: [{ name, points, radius }], exclude }) -> { lineLine[], lineAsset[], minLineLine, minLineAsset }` — line-to-line and line-to-asset clearances [m] — API RP 2SK clearance requirements.

## lineprops.js — indicative component properties

- `CHAIN_GRADE_C` — break-load constants c: R3 0.0223, R3S 0.0249, R4 0.0274, R4S 0.0304, R5 0.0320 — IACS UR W22 / DNV-OS-E302.
- `CHAIN_COEFFS` — mass 19.9e3·d² (studless) / 21.9e3·d² (stud-link) kg/m, EA 0.854e11·d² / 1.01e11·d² N, d in m — OrcaFlex chain estimators.
- `chainMBL({ d, grade, c }) -> number` — c·d²(44 − 0.08d) kN, d in mm, returned in N — IACS UR W22.
- `chainProperties({ d, grade, type, massCoeff, eaCoeff, c, rhoSteel, rhoWater, g }) -> { d, mass, w, EA, MBL, dVol }` — [m, kg/m, N/m, N, N, m] — as above; submerged weight from displaced steel volume.
- `WIRE_COEFFS` — mass / EA / MBL coefficients for six-strand (steel or fibre core) and spiral strand — OrcaFlex wire estimators (≡ EN 12385-4 K·d²·Rr, grade 1770); spiral strand indicative.
- `wireRopeProperties({ d, construction, massCoeff, eaCoeff, mblCoeff, dOuter, … }) -> { d, mass, w, EA, MBL, dVol }` — as above.
- `FIBRE_COEFFS` — polyester / nylon (8-strand, OrcaFlex estimators) and HMPE (indicative) coefficients and fibre densities.
- `fibreRopeProperties({ d, material, massCoeff, eaCoeff, mblCoeff, rhoFibre, … }) -> { d, mass, w, EA, MBL, dVol }` — as above.
- `FIBRE_STIFFNESS_PRESETS` — polyester-dynamic (α = 18.5, β = 33), polyester-static (α = 13, β = 0) — François & Davies (2008); BV NI 432.
- `fibreRopeStiffness({ MBL, Tmean, alpha, beta, preset }) -> { EA, Kr, alpha, beta }` — EA = (α + β·Tmean/MBL)·MBL [N] — François & Davies (2008).
- `corrodedDiameter({ d, rate, years }) -> number` — diameter after corrosion/wear allowance [m] — DNV-OS-E301 Ch.2 Sec.2 (rate supplied by the user).
- `HYDRO_COEFFS` — drag coefficients: stud chain 2.6/1.4, studless chain 2.4/1.15, six-strand rope 1.8, spiral rope 1.6 (unsheathed) / 1.2 (sheathed), fibre rope 1.6 — DNV-OS-E301 Ch.2 Sec.2.
- `hydroCoefficients({ type, Cdn, Cdt, Can = 1, Cat, ropeCdt = 0 }) -> { Cdn, Cdt, Can, Cat, source }` — drag per DNV-OS-E301; added mass and rope tangential drag are documented assumptions.
- `marineGrowth({ d, thickness, rhoGrowth = 1325, component, mu, Cdn, … }) -> { mass, w, CdnGrowth, dGrowth }` — marine growth mass, weight and drag increase [kg/m, N/m, -, m] — DNV-OS-E301 Ch.2 Sec.2.
- `toSegment(props, L, name) -> { L, w, EA, MBL, mass, d, name }` — segment record for multiseg.js.

## lumpedmass.js — 3-D lumped-mass line dynamics

- `airyWaveNumber(omega, depth, g) -> number` — k from ω² = g·k·tanh(k·h) [rad/m] — DNV-RP-C205 section 3.
- `buildLumpedLine({ anchor, fairlead, sections: [{ L, EA, mass, d, nSeg, Cdn, Cdt, Can, Cat, BA | zeta, EI }], points, rho, g, seabed: { z, k, c | zeta, mu, vc }, current, wave, waterLevel, allowCompression }) -> model` — node/segment model initialised on the catenary shape — Hall & Goupee (2015).
- `relaxLine(model, { tol, maxIter, cfl }) -> { converged, iterations, residual, fairleadTension, anchorTension, segmentTension, H, groundedNodes, touchdownS }` — static state by dynamic relaxation with kinetic damping [N, m] — Day (1965); Barnes (1999).
- `simulateLine(model, { motion, duration, dt, integrator = 'rk4' | 'semi-implicit', cfl, dtInternal, segments, quasiStatic, statsFrom, relax }) -> { t, fairleadTension, anchorTension, segmentTension, quasiStaticTension, touchdownS, energy, envelopes, static, stats: { …, daf }, steps, dtInternal }` — time-domain response to prescribed fairlead motion with sub-stepping on the axial wave speed [s, N, m, J] — Hall & Goupee (2015); Morison equation per DNV-RP-C205 section 6.
- `detectSnapLoads({ t, T, Tref, slackFraction = 0.1, spikeFactor = 1 }) -> { slack, snap, events[], minTension, maxTension, threshold }` — slack → re-tension screening [s, N] — slack criterion in the style of DNV-RP-N103 / DNV-OS-E301 (thresholds are arguments).
- `lineBreakTransient({ system, brokenLine, Fext, mass: [mx, my, Izz], damping, t0, duration, dt }) -> { t, x, y, yaw, offset, maxTension, intact, damaged, maxOffset, maxExcursion, staticExcursion, overshootFactor, peakTension, peakTensionLine, peakTensionTime, lineMaxTension }` — surge/sway/yaw floater transient after a line release with quasi-static mooring restoring [s, m, rad, N] — API RP 2SK transient analysis; first-principles equations of motion.

## interference.js — current deflection, wake shielding, clearance

- `currentDeflection({ L, w, d, Cd, current, mode = 'top-tensioned' | 'free-hanging', topTension, endWeight, endForce, rho, n, dragFactor }) -> { s, y, slope, Te, f, yMax, sAtMax, maxAngle, topAngle, bottomAngle, totalDrag }` — tensioned-string deflection under steady current [m, -, N, N/m, rad] — Sparks (2007) ch. 6–7.
- `huseWake({ x, y, D, Cd, U0, k1 = 0.25, k2 = 1.0 }) -> { deficit, velocity, halfWidth, xs, dragFactor }` — wake velocity deficit and shielding factor [m/s, m, -] — Huse (1993) OTC 7237; DNV-RP-F203.
- `blevinsWakeDrag({ x, y, D, Cd, Cd0, a1 = 1.0, a2 = 4.5 }) -> { Cd, dragFactor }` — downstream drag coefficient in a wake [-] — Blevins (2005) JOMAE 127.
- `clearanceUnderCurrent({ upstream, downstream, current, spacing, transverseOffset, wakeModel, wake, n, rho, maxIter, tol }) -> { s, clearance, distance, dragFactor, upstream, downstream, minClearance, sAtMin, collision, iterations, converged }` — clearance of two lines along their length with wake shielding [m] — DNV-RP-F203 approach.
- `impactVelocityEstimate({ force, gap, mass, stiffness }) -> { velocity, energy, reaches }` — work–energy impact velocity bound [m/s, J/m] — first principles (screening).

## berthmooring.js — moored ship at a berth

- `KNOT` — 1852/3600 m/s.
- `GENERIC_COEFFICIENTS` — generic coefficient magnitudes by vessel type (tanker, lng, container, bulk); placeholders, not OCIMF/SIGTTO data.
- `genericCoefficientCurves({ type, kind, Cx0, Cy0, Cxy0, shallowFactor }) -> { angleDeg, Cx, Cy, Cxy, generic: true }` — generic smooth curves vs angle of attack [-].
- `environmentalForce({ kind, speed, thetaDeg, vessel, curves, type, rho, shallowFactor }) -> { Fx, Fy, Mz, Cx, Cy, Cxy, generic }` — wind/current force in ship axes [N, N, N·m] — OCIMF MEG force definitions.
- `MEG_LINE_LIMITS` — line load limits: steel 0.55, synthetic 0.50, polyamide 0.45 of MBL — OCIMF MEG3 (MEG4 WLL equivalent).
- `ocimfStandardCases({ windStepDeg, windKn }) -> cases[]` — 60 kn wind with 3 kn 0°/180°, 2 kn 10°/170°, 0.75 kn 90° current [m/s, deg] — OCIMF MEG standard environmental criteria.
- `berthEquilibrium({ vessel, pose0, lines, fenders, wind, current, extra, reference, condition, limits, manifold, armEnvelope, bollards, tol, maxIter }) -> { pose, excursion, residual, environment, lines[], bollards[], fenders[], checks, converged, iterations }` — 3-DOF static equilibrium with tension-only lines and compression-only fenders; line, bollard, fender and loading-arm checks [m, deg, N, rad] — OCIMF MEG static mooring analysis; BS 6349-4.
- `runBerthCases(model, cases) -> { cases[], governing, maxLineUc, maxFenderUc, maxSurge, maxSway }` — envelope over environmental cases.
