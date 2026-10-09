import test from 'node:test';
import assert from 'node:assert/strict';
import {
  bendStiffener, sizeBendStiffener, bendRestrictor, bellmouth, GENERIC_FOAM_DENSITY, buoyancyModules,
  midWaterArch, capstan, pullInTension, clampSpacing, impactProtection, hangOffClamp, terminationLandingLoad,
} from '../../src/engine/xsection/ancillaries.js';

const close = (a, b, rel = 1e-9, abs = 0) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b) + abs, `${a} vs ${b}`);

test('bend stiffener reduces to the constant-EI tensioned beam: kappa0 = theta sqrt(T/EI), exponential decay', () => {
  const T = 200e3, EI = 50e3, theta = 0.05, lam = Math.sqrt(T / EI);
  // no stiffener material: the "stiffener" is just a length of bare pipe
  const r = bendStiffener({ tension: T, angle: theta, length: 1.5, rootDiameter: 0.2, tipDiameter: 0.2, innerDiameter: 0.2, Epu: 0, EIpipe: EI, linear: true, nSteps: 400 });
  close(r.rootCurvature, theta * lam, 1e-8);
  for (const i of [50, 200, 400]) close(r.curvature[i], theta * lam * Math.exp(-lam * r.s[i]), 1e-7);
  close(r.rootMoment, EI * theta * lam, 1e-8);
  // nonlinear elastica: exact root curvature 2 lambda sin(theta/2); small angle agrees with the linear value
  const nl = bendStiffener({ tension: T, angle: theta, length: 1.5, rootDiameter: 0.2, tipDiameter: 0.2, innerDiameter: 0.2, Epu: 0, EIpipe: EI });
  close(nl.rootCurvature, 2 * lam * Math.sin(theta / 2), 1e-8);
  close(nl.rootCurvature, theta * lam, 1e-3);
  const big = bendStiffener({ tension: T, angle: 1.0, length: 1.5, rootDiameter: 0.2, tipDiameter: 0.2, innerDiameter: 0.2, Epu: 0, EIpipe: EI });
  close(big.rootCurvature, 2 * lam * Math.sin(0.5), 1e-8);
  close(big.rootShear, T * Math.sin(1.0)); close(big.rootAxial, T * Math.cos(1.0));
});

test('uniform sleeve (two-zone tensioned beam): closed-form root curvature of the linear problem', () => {
  const T = 150e3, EIp = 40e3, theta = 0.1, L = 1.2, Epu = 60e6, Do = 0.5, di = 0.2;
  const EI1 = EIp + Epu * Math.PI * (Do ** 4 - di ** 4) / 64, l1 = Math.sqrt(T / EI1), c = Math.sqrt(T * EIp);
  // psi = theta - phi = theta cosh(l1 s) + B sinh(l1 s); tip condition -EI1 psi'(L) = c psi(L)
  const B = -theta * (EI1 * l1 * Math.sinh(l1 * L) + c * Math.cosh(l1 * L)) / (EI1 * l1 * Math.cosh(l1 * L) + c * Math.sinh(l1 * L));
  const r = bendStiffener({ tension: T, angle: theta, length: L, rootDiameter: Do, tipDiameter: Do, innerDiameter: di, Epu, EIpipe: EIp, linear: true, nSteps: 400 });
  close(r.rootCurvature, -B * l1, 1e-7);
  const psiL = theta * Math.cosh(l1 * L) + B * Math.sinh(l1 * L);
  close(r.tipCurvaturePipe, c * psiL / EIp, 1e-7);
  close(r.rootMomentStiffener, r.rootMoment * (1 - EIp / EI1), 1e-12);
  // the bare pipe just beyond a stiff sleeve sees the largest curvature
  assert.ok(r.tipCurvaturePipe > r.rootCurvature);
  close(r.maxCurvature, r.tipCurvaturePipe, 1e-12);
});

test('tapered stiffener: equilibrium along the beam, tip condition, RK4 convergence and benefit over a bare pipe', () => {
  const o = { tension: 300e3, angle: 0.35, length: 2.5, rootDiameter: 0.75, tipDiameter: 0.26, innerDiameter: 0.25, Epu: 70e6, EIpipe: 60e3, mbr: 2.5, strainLimit: 0.1 };
  const r = bendStiffener({ ...o, nSteps: 800 });
  // M' = -T sin(theta - phi) by central differences on the returned arrays
  for (const i of [100, 400, 700]) {
    const dM = (r.moment[i + 1] - r.moment[i - 1]) / (r.s[i + 1] - r.s[i - 1]);
    close(dM, -o.tension * Math.sin(o.angle - r.angle[i]), 1e-4);
  }
  const n = r.s.length - 1;
  close(r.moment[n], 2 * Math.sqrt(o.tension * o.EIpipe) * Math.sin((o.angle - r.angle[n]) / 2), 1e-8);
  // global moment equilibrium: root moment = T x lever arm of the tension line about the root. Beyond the tip the
  // pipe carries M = T x (distance to the asymptote), so M0 = T [x_L sin(theta) - y_L cos(theta)] + M(L)
  let x = 0, y = 0;
  for (let i = 0; i < n; i++) { const h = r.s[i + 1] - r.s[i], pm = 0.5 * (r.angle[i] + r.angle[i + 1]); x += h * Math.cos(pm); y += h * Math.sin(pm); }
  close(r.rootMoment, o.tension * (x * Math.sin(o.angle) - y * Math.cos(o.angle)) + r.moment[n], 1e-5);
  const coarse = bendStiffener({ ...o, nSteps: 100 });
  close(coarse.rootCurvature, r.rootCurvature, 1e-6);
  const bare = 2 * Math.sqrt(o.tension / o.EIpipe) * Math.sin(o.angle / 2);
  assert.ok(r.maxCurvature < 0.6 * bare, 'stiffener spreads the curvature');
  close(r.mbrUtilisation, r.maxCurvature * 2.5, 1e-12);
  assert.ok(r.maxPuStrain > 0 && r.strainUtilisation === r.maxPuStrain / 0.1);
  assert.throws(() => bendStiffener({ ...o, tipDiameter: 0.2 }), /tipDiameter/);
});

test('sizing loop returns a feasible least-volume design that meets the MBR', () => {
  const base = { tipDiameter: 0.26, innerDiameter: 0.25, Epu: 70e6, EIpipe: 60e3 };
  const z = sizeBendStiffener({ ...base, loadCases: [{ tension: 300e3, angle: 0.3 }, { tension: 150e3, angle: 0.45 }], lengths: [1.5, 2.0, 2.5, 3.0], rootDiameterMax: 1.0, mbr: 2.5, strainLimit: 0.1 });
  assert.equal(z.feasible, true);
  assert.ok(z.utilisation <= 1);
  for (const lc of [{ tension: 300e3, angle: 0.3 }, { tension: 150e3, angle: 0.45 }]) {
    const r = bendStiffener({ ...base, ...lc, length: z.length, rootDiameter: z.rootDiameter, nSteps: 200 });
    assert.ok(r.maxCurvature * 2.5 <= 1 + 1e-9);
  }
  close(z.volume, Math.PI * z.length * (z.rootDiameter ** 2 + z.rootDiameter * 0.26 + 0.26 ** 2) / 12 - Math.PI * 0.25 ** 2 * z.length / 4, 1e-12);
  assert.equal(sizeBendStiffener({ ...base, loadCases: [{ tension: 3e6, angle: 1.2 }], lengths: [0.5], rootDiameterMax: 0.3, mbr: 5 }).feasible, false);
});

test('bend restrictor: joints of the locked polygon lie on the locking radius; root moment T R (1 - cos theta)', () => {
  const p = 0.3, delta = 0.12, n = 8;
  const r = bendRestrictor({ pitch: p, lockAngle: delta, nElements: n, mbr: 2.0, tension: 50e3, angle: 0.6, momentCapacity: 60e3 });
  // build the polygon and compute the circumradius from three consecutive joints
  const pts = [[0, 0]]; let a = delta / 2;
  for (let i = 0; i < 3; i++) { const q = pts[pts.length - 1]; pts.push([q[0] + p * Math.cos(a), q[1] + p * Math.sin(a)]); a += delta; }
  const [A, B, C] = pts; const ab = Math.hypot(B[0] - A[0], B[1] - A[1]), bc = Math.hypot(C[0] - B[0], C[1] - B[1]), ca = Math.hypot(A[0] - C[0], A[1] - C[1]);
  const area = Math.abs((B[0] - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (B[1] - A[1])) / 2;
  close(r.lockingRadius, ab * bc * ca / (4 * area), 1e-10);
  close(r.lockingRadius, p / delta, 1e-3);
  close(r.totalAngle, 0.96); assert.equal(r.coversAngle, true);
  close(r.rootMoment, 50e3 * r.lockingRadius * (1 - Math.cos(0.6)));
  close(r.mbrUtilisation, 2.0 / r.lockingRadius);
  assert.equal(bendRestrictor({ pitch: p, lockAngle: delta, nElements: 3, tension: 1, angle: 0.6 }).coversAngle, false);
});

test('bellmouth constant-curvature profile', () => {
  const b = bellmouth({ radius: 3, maxAngle: Math.PI / 6, throatDiameter: 0.4, mbr: 2.4, tension: 90e3 });
  close(b.length, 1.5); close(b.flare, 3 * (1 - Math.cos(Math.PI / 6))); close(b.exitDiameter, 0.4 + 2 * b.flare);
  close(b.mbrUtilisation, 0.8); close(b.lineLoad, 30e3);
});

test('buoyancy modules: net buoyancy, count and pitch, clamp slip by hand', () => {
  const V = Math.PI * (1.2 ** 2 - 0.3 ** 2) / 4 * 1.5;
  const b = buoyancyModules({ outerDiameter: 1.2, innerDiameter: 0.3, moduleLength: 1.5, foamDensity: 450, hardwareWeight: 300, lossAllowance: 0.05, requiredBuoyancy: 800, sectionLength: 120, clampForce: 200e3, mu: 0.1, slipSafetyFactor: 1.5 });
  const B = (1025 - 450) * V * 9.80665 * 0.95 - 300;
  close(b.netBuoyancy, B, 1e-12);
  close(b.pitch, B / 800, 1e-12);
  assert.equal(b.count, Math.ceil(800 * 120 / B));
  assert.ok(b.providedBuoyancy >= 800);
  close(b.slipUtilisation, 1.5 * B / 20e3, 1e-12);
  // table lookup (generic, editable): 1500 m between 1000 m (400) and 2000 m (480) -> 440
  const t = buoyancyModules({ outerDiameter: 1.2, innerDiameter: 0.3, moduleLength: 1.5, depthRating: 1500 });
  close(t.foamDensity, 440); assert.equal(t.genericDensity, true); assert.equal(GENERIC_FOAM_DENSITY.generic, true);
  close(buoyancyModules({ outerDiameter: 1.2, innerDiameter: 0.3, moduleLength: 1.5, depthRating: 150, table: { depth: [100, 200], density: [300, 320] } }).foamDensity, 310);
});

test('mid-water arch: force polygon closes', () => {
  const B = 400e3, T1 = 120e3, a1 = 1.2, T2 = 90e3, a2 = 1.0;
  const r = midWaterArch({ netBuoyancy: B, T1, angle1: a1, T2, angle2: a2, tetherCapacity: 500e3, nTethers: 2 });
  // forces on the arch: buoyancy up, legs down-and-out, tether (-Th, -Tv)
  const fx = -T1 * Math.cos(a1) + T2 * Math.cos(a2) + r.tetherHorizontal;
  const fy = B - T1 * Math.sin(a1) - T2 * Math.sin(a2) - r.tetherVertical;
  close(fx, 0, 0, 1e-6); close(fy, 0, 0, 1e-6);
  close(r.tetherTension, Math.hypot(r.tetherVertical, r.tetherHorizontal));
  close(r.utilisation, r.tetherTension / 2 / 500e3);
  assert.equal(r.slack, false);
  assert.equal(midWaterArch({ netBuoyancy: 100e3, T1, angle1: a1, T2, angle2: a2 }).slack, true);
});

test('capstan T2 = T1 exp(mu theta) and pull-in march through a J-tube by hand', () => {
  close(capstan({ tensionIn: 10e3, mu: 0.3, wrapAngle: Math.PI / 2 }).tensionOut, 10e3 * Math.exp(0.15 * Math.PI));
  close(capstan({ tensionIn: 10e3, mu: 0.25, wrapAngle: Math.PI }).ratio, 2.19328, 1e-5);
  // capstan relation as the integral of dT = mu T dtheta
  let T = 1; const n = 100000; for (let i = 0; i < n; i++) T *= 1 + 0.3 * (2 / n);
  close(capstan({ tensionIn: 1, mu: 0.3, wrapAngle: 2 }).tensionOut, T, 1e-4);
  const r = pullInTension({ backTension: 5e3, capacity: 100e3, segments: [
    { type: 'straight', length: 50, inclination: 0, mu: 0.3, weight: 400 },
    { type: 'bend', angle: Math.PI / 2, mu: 0.3, radius: 3 },
    { type: 'straight', length: 30, inclination: Math.PI / 2, mu: 0.3, weight: 400 },
  ] });
  const t1 = 5e3 + 400 * 50 * 0.3, t2 = t1 * Math.exp(0.3 * Math.PI / 2), t3 = t2 + 400 * 30;
  close(r.tensions[0], t1); close(r.tensions[1], t2); close(r.pullTension, t3, 1e-12);
  close(r.maxSidewall, t2 / 3); close(r.utilisation, t3 / 100e3);
});

test('clamp spacing, impact protection, hang-off clamp and landing load by hand', () => {
  const c = clampSpacing({ load: 200, momentAllowable: 600, clampCapacity: 2000, EI: 5e3, deflectionLimit: 0.02, length: 100 });
  close(c.spacingMoment, 6); close(c.spacingClamp, 10); close(c.spacingDeflection, (384 * 5e3 * 0.02 / 1000) ** 0.25);
  assert.equal(c.governing, 'deflection'); assert.equal(c.count, Math.ceil(100 / c.spacing) + 1);
  const i = impactProtection({ mass: 500, dropHeight: 2, ratedEnergy: 15e3, length: 30, shellLength: 2 });
  close(i.impactEnergy, 500 * 9.80665 * 2); close(i.utilisation, i.impactEnergy / 15e3); assert.equal(i.shells, 15);
  close(impactProtection({ mass: 500, velocity: 3, ratedEnergy: 15e3 }).impactEnergy, 2250);
  const h = hangOffClamp({ tension: 150e3, nBolts: 8, boltPreload: 100e3, mu: 0.25, diameter: 0.2, clampLength: 0.8, allowablePressure: 10e6 });
  close(h.slipCapacity, 0.25 * 800e3 * 2 / 2); close(h.slipUtilisation, 0.75); close(h.contactPressure, 5e6); close(h.pressureUtilisation, 0.5);
  const l = terminationLandingLoad({ mass: 20e3, addedMass: 10e3, submergedWeight: 170e3, velocity: 0.5, stiffness: 5e6 });
  close(l.peakLoad, 170e3 + Math.sqrt(170e3 ** 2 + 5e6 * 30e3 * 0.25));
  // energy balance at maximum compression x = F/k
  const x = l.peakLoad / 5e6;
  close(0.5 * 30e3 * 0.25 + 170e3 * x, 0.5 * 5e6 * x * x, 1e-10);
  close(terminationLandingLoad({ mass: 20e3, submergedWeight: 170e3, velocity: 0, stiffness: 5e6 }).daf, 2);
});
