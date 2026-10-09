import test from 'node:test';
import assert from 'node:assert/strict';
import {
  createPolar, parseAeroDynPolar, parseCsvPolar, parseSeligCoordinates, viternaCdMax, extrapolateViterna,
  polarLinearFit, stallDelayDuSelig, stallDelaySnel, blendPolars, polarByThickness, liftSlopeThickness,
  genericPolar, cylinderPolar, cpMinEstimate, cavitationCheck, normalTangential, reynoldsNumber,
} from '../../src/engine/aero/polar.js';

const D2R = Math.PI / 180;
const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} got ${a}, expected ${b} +- ${tol}`);

// a simple attached-flow table: Cl = 2 pi (alpha + 2 deg), rounded stall at 14 deg
function baseTable() {
  const alphaDeg = []; const cl = []; const cd = []; const cm = [];
  for (let a = -10; a <= 14; a += 1) {
    alphaDeg.push(a);
    const lin = 2 * Math.PI * (a + 2) * D2R;
    cl.push(a <= 10 ? lin : lin - 0.012 * (a - 10) ** 2 * 2 * Math.PI);
    cd.push(0.008 + 0.0001 * a * a);
    cm.push(-0.05);
  }
  return createPolar({ name: 'test', alphaDeg, cl, cd, cm, thickness: 0.18 });
}

test('createPolar: linear interpolation, end clamping and wrapping of alpha', () => {
  const p = createPolar({ alphaDeg: [-180, 0, 90, 180], cl: [0, 0, 0.9, 0], cd: [0.02, 0.01, 1.3, 0.02] });
  close(p.clAt(45 * D2R), 0.45, 1e-12);
  close(p.cdAt(45 * D2R), 0.655, 1e-12);
  close(p.clAt((45 + 360) * D2R), 0.45, 1e-9, 'wrapped');
  close(p.clAt((45 - 720) * D2R), 0.45, 1e-9, 'wrapped');
  assert.equal(p.cmAt(0.3), 0);
  assert.throws(() => createPolar({ alphaDeg: [0, 0], cl: [0, 1], cd: [0, 0] }), RangeError);
  assert.throws(() => createPolar({ alphaDeg: [0, 1], cl: [0], cd: [0, 0] }), RangeError);
});

test('reynoldsNumber: 80 m/s, 2 m chord, air nu = 1.46e-5 -> 1.096e7', () => {
  close(reynoldsNumber(80, 2, 1.46e-5), 160 / 1.46e-5, 1e-3);
});

test('parseAeroDynPolar: header skipped, Reynolds number in millions, table read', () => {
  const txt = [
    'AeroDyn airfoil file.  Compatible with AeroDyn v13.0.',
    'Sample section',
    'more header',
    '     1         Number of airfoil tables in this file',
    '     6.0       Table ID parameter (Reynolds number in millions)',
    '    9.00       Stall angle (deg)',
    '   -4.0   -0.20   0.0090  -0.08',
    '    0.0    0.25   0.0070  -0.08',
    '    4.0    0.70   0.0080  -0.08',
    '    8.0    1.10   0.0120  -0.07',
    '   12.0    1.30   0.0300  -0.06',
  ].join('\n');
  const p = parseAeroDynPolar(txt, { thickness: 0.21 });
  assert.deepEqual(p.alphaDeg, [-4, 0, 4, 8, 12]);
  assert.equal(p.Re, 6e6);
  close(p.clAt(2 * D2R), 0.475, 1e-12);
  close(p.cmAt(10 * D2R), -0.065, 1e-12);
  assert.throws(() => parseAeroDynPolar('no numbers here'), RangeError);
});

test('parseCsvPolar: header row, separators and unsorted rows', () => {
  const p = parseCsvPolar('alpha,cl,cd\n5; 0.8; 0.01\n0, 0.3, 0.008\n-5\t-0.2\t0.009\n');
  assert.deepEqual(p.alphaDeg, [-5, 0, 5]);
  close(p.clAt(2.5 * D2R), 0.55, 1e-12);
});

function naca4(t, n, mCam = 0, pCam = 0.4) {
  // NACA 4-digit thickness distribution (closed trailing edge variant not used: standard -0.1015 coefficient)
  const yt = (x) => 5 * t * (0.2969 * Math.sqrt(x) - 0.1260 * x - 0.3516 * x * x + 0.2843 * x ** 3 - 0.1015 * x ** 4);
  const yc = (x) => (mCam === 0 ? 0 : (x < pCam ? mCam / pCam ** 2 * (2 * pCam * x - x * x) : mCam / (1 - pCam) ** 2 * (1 - 2 * pCam + 2 * pCam * x - x * x)));
  const xs = [];
  for (let i = 0; i <= n; i++) xs.push(0.5 * (1 - Math.cos(Math.PI * i / n)));
  const up = xs.map((x) => [x, yc(x) + yt(x)]);
  const lo = xs.map((x) => [x, yc(x) - yt(x)]);
  return { up, lo };
}

test('parseSeligCoordinates: NACA 0012 thickness 12 % at 30 % chord, zero camber, LE radius 1.1019 t^2', () => {
  const { up, lo } = naca4(0.12, 200);
  const lines = ['NACA 0012'];
  for (let i = up.length - 1; i >= 0; i--) lines.push(`${up[i][0].toFixed(7)} ${up[i][1].toFixed(7)}`);
  for (let i = 1; i < lo.length; i++) lines.push(`${lo[i][0].toFixed(7)} ${lo[i][1].toFixed(7)}`);
  const g = parseSeligCoordinates(lines.join('\n'));
  assert.equal(g.name, 'NACA 0012');
  close(g.thickness, 0.12, 2e-4);
  close(g.thicknessPos, 0.30, 0.01);
  close(g.camber, 0, 1e-6);
  close(g.leRadius, 1.1019 * 0.12 * 0.12, 0.10 * 1.1019 * 0.0144, 'LE radius within 10 %');
  // area of the 4-digit section: 0.685 t c^2 (integral of the thickness polynomial, open trailing edge)
  close(g.area, 0.6851 * 0.12, 2e-4);
});

test('parseSeligCoordinates: Lednicer format and cambered NACA 4412 (camber 4 % at 40 %)', () => {
  const { up, lo } = naca4(0.12, 120, 0.04, 0.4);
  const lines = ['NACA 4412 (Lednicer)', `${up.length}. ${lo.length}.`, ''];
  for (const p of up) lines.push(`${p[0].toFixed(7)} ${p[1].toFixed(7)}`);
  lines.push('');
  for (const p of lo) lines.push(`${p[0].toFixed(7)} ${p[1].toFixed(7)}`);
  const g = parseSeligCoordinates(lines.join('\n'));
  close(g.camber, 0.04, 3e-4);
  close(g.camberPos, 0.40, 0.01);
  close(g.thickness, 0.12, 5e-4);
});

test('viternaCdMax: 1.11 + 0.018 AR up to AR 50, 2.01 beyond (Viterna and Corrigan 1982)', () => {
  close(viternaCdMax(10), 1.29, 1e-12);
  close(viternaCdMax(17), 1.416, 1e-12);
  assert.equal(viternaCdMax(60), 2.01);
});

test('extrapolateViterna: published coefficients A1, A2, B1, B2; continuity at stall; flat-plate limits', () => {
  const base = baseTable();
  const AR = 10;
  const p = extrapolateViterna(base, { aspectRatio: AR });
  const as = 14 * D2R;
  const cls = base.cl[base.cl.length - 1];
  const cds = base.cd[base.cd.length - 1];
  const cdMax = 1.11 + 0.018 * AR;
  const A2 = (cls - cdMax * Math.sin(as) * Math.cos(as)) * Math.sin(as) / Math.cos(as) ** 2;
  const B2 = (cds - cdMax * Math.sin(as) ** 2) / Math.cos(as);
  close(p.viterna.A1, cdMax / 2, 1e-12); close(p.viterna.A2, A2, 1e-12);
  close(p.viterna.B1, cdMax, 1e-12); close(p.viterna.B2, B2, 1e-12);
  // continuity at the matching point
  close(p.clAt(as), cls, 1e-9); close(p.cdAt(as), cds, 1e-9);
  // hand values at 44 deg (a generated grid node)
  const a = 44 * D2R;
  close(p.clAt(a), 0.5 * cdMax * Math.sin(2 * a) + A2 * Math.cos(a) ** 2 / Math.sin(a), 1e-9);
  close(p.cdAt(a), cdMax * Math.sin(a) ** 2 + B2 * Math.cos(a), 1e-9);
  // plate normal to the flow
  close(p.clAt(Math.PI / 2), 0, 1e-9); close(p.cdAt(Math.PI / 2), cdMax, 1e-9);
  close(p.cdAt(-Math.PI / 2), cdMax, 1e-9);
  // zero lift at +-180 deg, full coverage, reflected branch reduced by 0.7
  close(p.clAt(Math.PI), 0, 1e-9); close(p.clAt(-Math.PI), 0, 1e-9);
  assert.equal(p.alphaDeg[0], -180); assert.equal(p.alphaDeg[p.alphaDeg.length - 1], 180);
  close(p.clAt(136 * D2R), -0.7 * p.clAt(44 * D2R), 1e-9);
  for (let i = 1; i < p.alphaDeg.length; i++) assert.ok(p.alphaDeg[i] > p.alphaDeg[i - 1]);
  // the original table is untouched inside its range
  close(p.clAt(3 * D2R), base.clAt(3 * D2R), 1e-12);
  // centre of pressure at mid chord for the plate at 90 deg: Cm(c/4) = -Cn/4
  close(p.cmAt(Math.PI / 2), -cdMax * 0.25, 1e-9);
});

test('polarLinearFit: zero-lift angle and slope of a linear table', () => {
  const f = polarLinearFit(baseTable());
  close(f.alpha0, -2 * D2R, 1e-9);
  close(f.slope, 2 * Math.PI, 1e-6);
});

test('stallDelayDuSelig: f_L by hand for c/r = 0.3, r/R = 0.3, tip-speed ratio 7 (Du and Selig 1998)', () => {
  const p = extrapolateViterna(baseTable(), { aspectRatio: 10 });
  const q = stallDelayDuSelig(p, { cOverR: 0.3, rOverR: 0.3, tsr: 7 });
  // Lambda = 7/sqrt(50) = 0.98995; exponent = 1/(0.98995*0.3) = 3.3672; 0.3^3.3672 = 0.017353
  // f_L = (1/2pi) (1.6*0.3/0.1267 * (1 - 0.017353)/(1 + 0.017353) - 1) = 0.4232
  close(q.stallDelayFactor, 0.4232, 2e-4);
  // attached flow unchanged, post-stall lift raised towards the potential line, drag changed per Eggers
  close(q.clAt(4 * D2R), p.clAt(4 * D2R), 1e-9);
  const a = 20 * D2R;
  const clp = 2 * Math.PI * (a + 2 * D2R);
  const dCl = 0.42323 * (clp - p.clAt(a));
  close(q.clAt(a) - p.clAt(a), dCl, 2e-3);
  close(q.cdAt(a) - p.cdAt(a), dCl * (Math.sin(a) - 0.12 * Math.cos(a)) / (Math.cos(a) + 0.12 * Math.sin(a)), 2e-3);
  // negative angles are not modified
  close(q.clAt(-30 * D2R), p.clAt(-30 * D2R), 1e-12);
});

test('stallDelaySnel: Cl3D = Cl2D + 3 (c/r)^2 (Cl,pot - Cl2D)', () => {
  const p = extrapolateViterna(baseTable(), { aspectRatio: 10 });
  const q = stallDelaySnel(p, { cOverR: 0.25 });
  const a = 20 * D2R;
  close(q.clAt(a), p.clAt(a) + 3 * 0.0625 * (2 * Math.PI * (a + 2 * D2R) - p.clAt(a)), 1e-6);
  close(q.cdAt(a), p.cdAt(a), 1e-12);
});

test('blendPolars / polarByThickness: linear in thickness between family members', () => {
  const a = genericPolar({ thickness: 0.21 });
  const b = genericPolar({ thickness: 0.25 });
  const mid = polarByThickness([b, a], 0.23);
  close(mid.thickness, 0.23, 1e-12);
  for (const deg of [-3, 2, 7, 20, 90]) {
    close(mid.clAt(deg * D2R), 0.5 * (a.clAt(deg * D2R) + b.clAt(deg * D2R)), 1e-9);
    close(mid.cdAt(deg * D2R), 0.5 * (a.cdAt(deg * D2R) + b.cdAt(deg * D2R)), 1e-9);
  }
  assert.equal(polarByThickness([a, b], 0.1), a);
  assert.equal(polarByThickness([a, b], 0.4), b);
  assert.throws(() => blendPolars(a, b, 1.5), RangeError);
});

test('liftSlopeThickness: 2 pi (1 + 0.77 t/c) eta', () => {
  close(liftSlopeThickness(0, 1), 2 * Math.PI, 1e-12);
  close(liftSlopeThickness(0.18, 0.9), 2 * Math.PI * 1.1386 * 0.9, 1e-12);
});

test('genericPolar: prescribed slope, zero-lift angle, stall point, drag bucket and Viterna beyond; labelled generic', () => {
  const p = genericPolar({ thickness: 0.18, alpha0Deg: -4, alphaStallDeg: 13, clMax: 1.45, cdMin: 0.006, kBucket: 0.005, aspectRatio: 17 });
  assert.match(p.label, /generic parametric polar/);
  assert.equal(p.generic, true);
  const m = 2 * Math.PI * (1 + 0.77 * 0.18) * 0.9;
  close(p.clAt(-4 * D2R), 0, 1e-9);
  close((p.clAt(0) - p.clAt(-2 * D2R)) / (2 * D2R), m, 1e-6);
  close(p.clAt(13 * D2R), 1.45, 1e-9);
  // the stall point is the lift maximum on the positive side
  let mx = -Infinity;
  for (let a = 0; a <= 90; a += 0.25) mx = Math.max(mx, p.clAt(a * D2R));
  close(mx, 1.45, 1e-9);
  // minimum drag at the lift of zero incidence, drag rises quadratically with the linear lift
  close(p.cdAt(0), 0.006, 1e-9);
  close(p.cdAt(4 * D2R), 0.006 + 0.005 * (m * 4 * D2R) ** 2, 1e-9);
  // flat plate at 90 deg with Cd,max = 1.11 + 0.018 * 17
  close(p.cdAt(Math.PI / 2), 1.416, 1e-9);
  // thin-airfoil moment of the equivalent parabolic camber line
  close(p.cmAt(2 * D2R), Math.PI * (-4 * D2R) / 2, 1e-9);
  // mirror-image negative stall
  close(p.clAt(-21 * D2R), -1.45, 1e-9);
  // inconsistent stall point is rejected
  assert.throws(() => genericPolar({ thickness: 0.18, alpha0Deg: -4, alphaStallDeg: 13, clMax: 2.5 }), RangeError);
  assert.throws(() => genericPolar({ thickness: 0.9 }), RangeError);
});

test('genericPolar defaults exist for the whole 18-40 % thickness range', () => {
  for (const t of [0.18, 0.21, 0.25, 0.30, 0.35, 0.40, 0.33]) {
    const p = genericPolar({ thickness: t });
    assert.ok(p.params.clMax > 1 && p.params.clMax < 1.6);
    assert.ok(Number.isFinite(p.clAt(0.1)) && p.cdAt(0.1) > 0);
  }
});

test('cylinderPolar: no lift, constant drag', () => {
  const p = cylinderPolar({ cd: 0.35 });
  for (const a of [-3, -1, 0, 0.5, 2.5]) { assert.equal(p.clAt(a), 0); assert.equal(p.cdAt(a), 0.35); }
});

test('cpMinEstimate: ellipse at zero lift 1 - (1 + t/c)^2; lift adds Cl/4 to the velocity ratio', () => {
  close(cpMinEstimate({ cl: 0, thickness: 0.2 }), 1 - 1.2 ** 2, 1e-12);
  close(cpMinEstimate({ cl: 0.8, thickness: 0.2 }), 1 - 1.4 ** 2, 1e-12);
  close(cpMinEstimate({ cl: -0.8, thickness: 0.2 }), 1 - 1.4 ** 2, 1e-12);
});

test('cavitationCheck: sigma by hand at 10 m submergence and 12 m/s', () => {
  // (101325 + 1025*9.80665*10 - 1700) / (0.5*1025*144) = 200143.2/73800 = 2.7120
  const c = cavitationCheck({ cpMin: -2.0, W: 12, depth: 10 });
  close(c.sigma, (101325 + 1025 * 9.80665 * 10 - 1700) / 73800, 1e-9);
  close(c.sigma, 2.712, 1e-3);
  assert.equal(c.cavitates, false);
  assert.equal(cavitationCheck({ cpMin: -3.0, W: 12, depth: 10 }).cavitates, true);
});

test('normalTangential: resolution into the rotor plane', () => {
  const r = normalTangential(1.0, 0.01, 10 * D2R);
  close(r.cn, Math.cos(10 * D2R) + 0.01 * Math.sin(10 * D2R), 1e-12);
  close(r.ct, Math.sin(10 * D2R) - 0.01 * Math.cos(10 * D2R), 1e-12);
});
