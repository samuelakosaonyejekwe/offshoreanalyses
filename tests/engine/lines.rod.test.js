import test from 'node:test';
import assert from 'node:assert/strict';
import { solveRodStatic, rodModes, rodHarmonicResponse, cableShape } from '../../src/engine/lines/rod.js';

const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

function cantilever(nElem, PL2overEI = 1) {
  const L = 10, EI = 1e6, P = PL2overEI * EI / (L * L);
  return solveRodStatic({
    sections: [{ length: L, EA: 1e12, EI, w: 0 }], nElem,
    endA: { x: 0, z: 0, angle: 0 }, endB: { Fz: -P },
  });
}

test('Euler elastica: tip-loaded cantilever PL²/EI = 1 (Bisshopp & Drucker 1945; Mattiasson 1981: δv/L = 0.3017, δh/L = 0.0564)', () => {
  const r = cantilever(40);
  const n = r.x.length - 1;
  assert.ok(Math.abs(-r.z[n] / 10 - 0.3017) < 3e-4, `vertical ${-r.z[n] / 10}`);
  assert.ok(Math.abs((10 - r.x[n]) / 10 - 0.0564) < 3e-4, `horizontal ${(10 - r.x[n]) / 10}`);
  // root moment = P × horizontal lever arm (statics of the deformed shape)
  assert.ok(rel(Math.abs(r.moment[0]), 1e4 * r.x[n]) < 1e-3);
  assert.ok(r.residual < 1e-5 * 1e4);
});

test('Euler elastica PL²/EI = 2 and 5 (Mattiasson 1981 table: δv/L = 0.49346, 0.71379; δh/L = 0.16064, 0.38763)', () => {
  const a = cantilever(60, 2), b = cantilever(60, 5);
  const n = a.x.length - 1;
  assert.ok(Math.abs(-a.z[n] / 10 - 0.49346) < 4e-4);
  assert.ok(Math.abs((10 - a.x[n]) / 10 - 0.16064) < 4e-4);
  assert.ok(Math.abs(-b.z[n] / 10 - 0.71379) < 6e-4);
  assert.ok(Math.abs((10 - b.x[n]) / 10 - 0.38763) < 6e-4);
});

test('mesh convergence of the corotational rod is second order (elastica tip deflection, n = 4, 8, 16)', () => {
  const tip = (n) => { const r = cantilever(n); return -r.z[n]; };
  const f1 = tip(4), f2 = tip(8), f3 = tip(16);
  const order = Math.log2(Math.abs(f1 - f2) / Math.abs(f2 - f3));
  assert.ok(order > 1.7 && order < 2.4, `observed order ${order}`);
});

test('EI → 0 reproduces the elastic-free catenary between two supports (closed form y = a cosh(x/a))', () => {
  const L = 120, X = 100, w = 50;
  // sqrt(L² − h²) = 2a sinh(X/2a) with h = 0, solved by bisection here
  let lo = 1, hi = 1e4;
  for (let k = 0; k < 200; k++) { const a = 0.5 * (lo + hi); if (2 * a * Math.sinh(X / (2 * a)) > L) lo = a; else hi = a; }
  const a = 0.5 * (lo + hi);
  const sag = a * (Math.cosh(X / (2 * a)) - 1), H = w * a, Tend = w * (a + sag);
  const r = solveRodStatic({
    sections: [{ length: L, EA: 1e11, EI: 1e-2, w }], nElem: 240,
    endA: { x: 0, z: 0 }, endB: { x: X, z: 0 }, initial: 'cable',
  });
  let zmin = 0; for (const v of r.z) zmin = Math.min(zmin, v);
  assert.ok(rel(-zmin, sag) < 2e-4, `sag ${-zmin} vs ${sag}`);
  assert.ok(rel(-r.endA.Fx, H) < 2e-4, `H ${-r.endA.Fx} vs ${H}`);
  assert.ok(rel(r.endA.Fz + r.endB.Fz, w * L) < 1e-6);       // vertical equilibrium
  assert.ok(rel(r.Te[0], Tend) < 1e-3);
  assert.ok(rel(r.Te[120], H) < 1e-3);
  // the starting chain itself is already the catenary
  const c = cableShape({ sections: [{ length: L, EA: 1e11, EI: 1, w }], nElem: 240, endA: { x: 0, z: 0 }, endB: { x: X, z: 0 } });
  assert.ok(rel(-Math.min(...c.z), sag) < 1e-3);
});

test('catenary with touchdown on a stiff seabed: lay-back and top tension vs closed form (T = H + w h)', () => {
  const w = 800, h = 300, H = 60e3, a = H / w;
  const S = Math.sqrt(h * (h + 2 * a)), Xs = a * Math.asinh(S / a), Lb = 150, L = S + Lb;
  const r = solveRodStatic({
    sections: [{ length: L, EA: 5e11, EI: 10, w }], nElem: 400,
    endA: { x: 0, z: 0 }, endB: { x: Lb + Xs, z: h }, seabed: { z: 0, kn: 1e7 }, initial: 'cable',
  });
  const n = r.x.length - 1;
  assert.ok(rel(r.Te[n], H + w * h) < 1e-3, `top tension ${r.Te[n]}`);
  assert.ok(rel(-r.endA.Fx, H) < 2e-3);
  assert.ok(rel(Math.tan(r.angle[n]), S / a) < 5e-3);
  // seabed carries exactly the grounded weight
  let Rv = 0; for (let k = 0; k <= n; k++) Rv += r.reaction[k] * (k === 0 || k === n ? 0.5 : 1) * (L / n);
  assert.ok(rel(Rv + r.endB.Fz + r.endA.Fz, w * L) < 1e-6);
  assert.ok(rel(Rv, w * Lb) < 0.01);
});

test('tensioned beam under uniform load: midspan deflection vs closed form with k = sqrt(T/EI)', () => {
  const L = 50, EI = 2e7, T = 4e5, q = 100, k = Math.sqrt(T / EI);
  // y(L/2) = q L²/(8T) − (q/(T k²)) (1 − 1/cosh(kL/2))   (Timoshenko, Strength of Materials II, tie-rod)
  const exact = q * L * L / (8 * T) - (q / (T * k * k)) * (1 - 1 / Math.cosh(k * L / 2));
  const Mexact = (q / (k * k)) * (1 - 1 / Math.cosh(k * L / 2));
  const r = solveRodStatic({
    sections: [{ length: L, EA: 1e13, EI, w: q }], nElem: 100,
    endA: { x: 0, z: 0 }, endB: { z: 0, Fx: T },
  });
  assert.ok(rel(-r.z[50], exact) < 2e-3, `${-r.z[50]} vs ${exact}`);
  assert.ok(rel(Math.abs(r.moment[50]), Mexact) < 3e-3, `${r.moment[50]} vs ${Mexact}`);
  assert.ok(rel(r.Te[50], T) < 1e-3);
});

test('rod pulled over a rigid circular support takes curvature 1/R where fully seated', () => {
  const R = 20, EI = 1e5;
  const r = solveRodStatic({
    sections: [{ length: 30, EA: 1e10, EI, w: 0 }], nElem: 150,
    endA: { x: 0, z: 0, angle: 0 }, endB: { Fz: -3000 },
    supports: [{ xc: 0, zc: -R, R, k: 1e8 }], nSteps: 20,
  });
  for (let k = 10; k <= 50; k++) {          // s = 2 … 10 m, well inside the wrapped length
    assert.ok(rel(-r.curvature[k], 1 / R) < 0.01, `s=${r.s[k]} curvature ${r.curvature[k]}`);
    assert.ok(Math.abs(Math.hypot(r.x[k], r.z[k] + R) - R) < 1e-3);
  }
  assert.ok(rel(Math.abs(r.moment[20]), EI / R) < 0.01);
  // beyond lift-off the rod is free: no contact reaction at the tip
  assert.equal(r.reaction[150], 0);
});

test('seabed Coulomb friction: tension of a pulled on-bottom line decays at μ·w per metre', () => {
  const w = 1000, mu = 0.5, L = 200, F = 40e3;
  const r = solveRodStatic({
    sections: [{ length: L, EA: 1e9, EI: 1e4, w }], nElem: 200,
    endA: { x: 0, z: 0 }, endB: { Fx: F }, seabed: { z: 0, kn: 1e6, mu, kt: 1e7 }, nSteps: 10,
  });
  // mobilised length F/(μ w) = 80 m: tension at 40 m from the pulled end = F − μ w 40
  assert.ok(rel(r.Te[160], F - mu * w * 40) < 0.02, `${r.Te[160]}`);
  assert.ok(Math.abs(r.Te[60]) < 0.02 * F);
  assert.ok(rel(r.friction[190], -mu * w) < 0.01);
});

test('natural frequencies of a tensioned pinned beam: fn = (n/2L) sqrt(T/m) sqrt(1 + n²π²EI/(T L²))', () => {
  const L = 100, EI = 5e6, T = 2e5, m = 120, ma = 80;
  const r = solveRodStatic({
    sections: [{ length: L, EA: 1e11, EI, w: 0, m, ma }], nElem: 80,
    endA: { x: 0, z: 0 }, endB: { z: 0, Fx: T },
  });
  const md = rodModes(r, { nModes: 4 });
  for (let n = 1; n <= 4; n++) {
    const fn = (n / (2 * L)) * Math.sqrt(T / (m + ma)) * Math.sqrt(1 + n * n * Math.PI ** 2 * EI / (T * L * L));
    assert.ok(rel(md.modes[n - 1].fHz, fn) < 1e-3, `mode ${n}: ${md.modes[n - 1].fHz} vs ${fn}`);
  }
  // mode 1 is a half sine with curvature (π/L)² at midspan per unit amplitude
  assert.ok(rel(Math.abs(md.modes[0].curvature[40]), (Math.PI / L) ** 2) < 5e-3);
  assert.ok(Math.abs(Math.abs(md.modes[0].phi[40]) - 1) < 1e-6);
});

test('harmonic end excitation of a taut string below resonance: y(x) = a sin(kx)/sin(kL)', () => {
  const L = 100, T = 1e5, m = 50, a = 0.1;
  const r = solveRodStatic({
    sections: [{ length: L, EA: 1e10, EI: 1, w: 0, m }], nElem: 100,
    endA: { x: 0, z: 0 }, endB: { x: L * (1 + T / 1e10), z: 0 },
  });
  const f1 = Math.sqrt(T / m) / (2 * L);
  const omega = 2 * Math.PI * 0.6 * f1, k = omega * Math.sqrt(m / T);
  const h = rodHarmonicResponse(r, { omega, end: 'B', zAmp: a });
  assert.ok(h.converged);
  for (const i of [25, 50, 75]) assert.ok(rel(h.zAmp[i], Math.abs(a * Math.sin(k * i) / Math.sin(k * L))) < 2e-3);
  // with drag the response at resonance stays finite and below the undamped off-resonance trend
  const rd = solveRodStatic({
    sections: [{ length: L, EA: 1e10, EI: 1, w: 0, m, D: 0.3, Cdn: 1.2 }], nElem: 100,
    endA: { x: 0, z: 0 }, endB: { x: L * (1 + T / 1e10), z: 0 },
  });
  const hd = rodHarmonicResponse(rd, { omega: 2 * Math.PI * f1, end: 'B', zAmp: a });
  assert.ok(hd.converged && Math.max(...hd.zAmp) < 3 && Math.max(...hd.zAmp) > a);
});

test('current drag on a vertical tensioned string: midspan offset = q L²/(8T) with q = ½ρ Cd D U²', () => {
  const L = 200, T = 2e6, D = 0.5, Cd = 1.0, U = 1.0, rho = 1025;
  const q = 0.5 * rho * Cd * D * U * U;
  const r = solveRodStatic({
    sections: [{ length: L, EA: 1e11, EI: 1e3, w: 0, D, Cdn: Cd }], nElem: 100,
    endA: { x: 0, z: 0 }, endB: { x: 0, Fz: T }, initialAngle: Math.PI / 2, current: { rho, U },
  });
  assert.ok(rel(r.x[50], q * L * L / (8 * T)) < 5e-3, `${r.x[50]}`);
});

test('effective vs true wall tension: Tw = Te + pi·Ai − pe·Ae, and pressure-induced axial strain', () => {
  const L = 100, Ai = 0.03, Ae = 0.04, EA = 2e9, T = 5e5, pi = 2e7;
  const r = solveRodStatic({
    sections: [{ length: L, EA, EI: 1e6, w: 0, Ai, Ae }], nElem: 10,
    endA: { x: 0, z: 0 }, endB: { z: 0, Fx: T },
    pressure: { rhoInt: 0, pIntRef: pi, zIntRef: 0, rhoExt: 0, zSurface: 0, nu: 0.3 },
  });
  assert.ok(rel(r.Te[5], T) < 1e-9);
  assert.ok(rel(r.Tw[5], T + pi * Ai) < 1e-9);
  // elongation = (Te + (1 − 2ν) pi Ai) L / EA
  assert.ok(rel(r.x[10] - L, (T + 0.4 * pi * Ai) * L / EA) < 1e-6);
});

test('invalid input throws RangeError naming the argument', () => {
  assert.throws(() => solveRodStatic({ sections: [{ length: -1, EA: 1, EI: 1, w: 0 }] }), /sections\[0\]\.length/);
  assert.throws(() => solveRodStatic({ sections: [] }), RangeError);
});
