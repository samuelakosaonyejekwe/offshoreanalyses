import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  gauss1d, gaussTriangle, invertSmall, rectangleMesh, annulusMesh, tubeWallMesh, plateWithHoleMesh, mappedMesh,
  refineUniform, meshArea, boundaryEdges, nodesWhere, edgesWhere, solvePlane, zzError, solveHeat2d, solvePlate,
  bucklePlate, modalPlate, solveShell, stiffenedPanelMesh,
} from '../../src/engine/structural/fem2d.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b} (rel ${rel})`);

// Smooth interior distortion that leaves the outline of the rectangle unchanged.
const wobble = (lx, ly, amp) => (x, y) => {
  const s = Math.sin((Math.PI * x) / lx) * Math.sin((Math.PI * y) / ly);
  return [x + amp * lx * s * (0.6 + y / ly), y - amp * ly * s * (0.4 + x / lx)];
};

test('quadrature rules integrate polynomials exactly', () => {
  for (let n = 1; n <= 4; n++) {
    const g = gauss1d(n);
    const deg = 2 * n - 2;
    let s = 0;
    for (let i = 0; i < n; i++) s += g.w[i] * g.x[i] ** deg;
    close(s, 2 / (deg + 1), 1e-13, `gauss1d(${n})`);
  }
  // ∫ ξ^p η^q over the unit triangle = p! q! / (p + q + 2)!
  const t7 = gaussTriangle(7);
  const fact = [1, 1, 2, 6, 24, 120, 720, 5040];
  for (const [p, q] of [[0, 0], [2, 1], [3, 2], [5, 0], [1, 4]]) {
    let s = 0;
    for (let i = 0; i < 7; i++) s += t7.w[i] * t7.xi[i] ** p * t7.eta[i] ** q;
    close(s, (fact[p] * fact[q]) / fact[p + q + 2], 1e-13, `triangle ${p},${q}`);
  }
  const a = Float64Array.from([4, 1, 2, 1, 3, 0, 2, 0, 5]);
  const inv = invertSmall(Float64Array.from(a), 3);
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
    let s = 0;
    for (let k = 0; k < 3; k++) s += a[i * 3 + k] * inv[k * 3 + j];
    assert.ok(Math.abs(s - (i === j ? 1 : 0)) < 1e-13);
  }
});

test('mesh generators: areas, boundary sets and refinement', () => {
  for (const type of ['T3', 'Q4', 'T6', 'Q8']) {
    const m = rectangleMesh({ lx: 2, ly: 1, nx: 4, ny: 3, type, distort: wobble(2, 1, 0.08) });
    close(meshArea(m).area, 2, 1e-12, `rectangle ${type}`);
    assert.equal(m.edges.left.length, 3);
    assert.equal(m.edges.top.length, 4);
    assert.equal(boundaryEdges(m).length, 14);
  }
  // Quadratic elements follow the arc: area of a quarter annulus to high accuracy.
  const ann = annulusMesh({ ri: 1, ro: 2, nr: 3, ntheta: 8, type: 'Q8' });
  close(meshArea(ann).area, (Math.PI * (4 - 1)) / 4, 2e-5, 'annulus Q8');
  const hole = plateWithHoleMesh({ a: 1, w: 5, h: 4, nr: 6, ntheta: 4, type: 'Q8' });
  close(meshArea(hole).area, 20 - Math.PI / 4, 2e-4, 'plate with hole');
  assert.ok(hole.edges.hole.length === 8 && hole.edges.right.length > 0 && hole.edges.symX.length === 6);
  const r = refineUniform(rectangleMesh({ lx: 1, ly: 1, nx: 2, ny: 2, type: 'T3' }));
  assert.equal(r.ne, 32);
  close(meshArea(r).area, 1, 1e-13);
  assert.throws(() => rectangleMesh({ lx: -1, ly: 1, nx: 1, ny: 1 }), RangeError);
});

test('patch test: every plane element reproduces a constant strain state exactly on distorted meshes', () => {
  const E = 210e9; const nu = 0.3;
  // u = a0 + a1 x + a2 y, v = b0 + b1 x + b2 y
  const ux = (x, y) => 1e-4 + 2e-4 * x - 1.5e-4 * y;
  const uy = (x, y) => -2e-4 + 0.5e-4 * x + 3e-4 * y;
  const exx = 2e-4; const eyy = 3e-4; const gxy = -1.5e-4 + 0.5e-4;
  const c = E / (1 - nu * nu);
  const sxx = c * (exx + nu * eyy); const syy = c * (eyy + nu * exx); const sxy = (E / (2 * (1 + nu))) * gxy;
  for (const [type, q4] of [['T3'], ['Q4', 'incompatible'], ['Q4', 'full'], ['T6'], ['Q8']]) {
    const mesh = rectangleMesh({ lx: 2, ly: 1.2, nx: 4, ny: 3, type, distort: wobble(2, 1.2, 0.07) });
    const bn = new Set();
    for (const ed of boundaryEdges(mesh)) for (const n of ed) bn.add(n);
    const fix = [...bn].map((n) => ({ nodes: [n], ux: ux(mesh.nodes[2 * n], mesh.nodes[2 * n + 1]), uy: uy(mesh.nodes[2 * n], mesh.nodes[2 * n + 1]) }));
    const r = solvePlane({ mesh, E, nu, mode: 'stress', fix, q4: q4 || 'incompatible' });
    let err = 0;
    for (let n = 0; n < mesh.nn; n++) {
      err = Math.max(err, Math.abs(r.u[2 * n] - ux(mesh.nodes[2 * n], mesh.nodes[2 * n + 1])), Math.abs(r.u[2 * n + 1] - uy(mesh.nodes[2 * n], mesh.nodes[2 * n + 1])));
    }
    assert.ok(err < 1e-12, `${type} ${q4 || ''} displacement error ${err}`);
    for (let n = 0; n < mesh.nn; n++) {
      close(r.stress.sxx[n], sxx, 1e-8, `${type} sxx`);
      close(r.stress.syy[n], syy, 1e-8, `${type} syy`);
      close(r.stress.sxy[n], sxy, 1e-8, `${type} sxy`);
    }
    // Interior nodes carry no reaction; the whole body is in equilibrium.
    assert.ok(Math.abs(r.reactionSum.fx) < 1e-6 * sxx && Math.abs(r.reactionSum.fy) < 1e-6 * sxx);
  }
});

test('cantilever in plane stress: incompatible-mode Q4 and Q8 match Timoshenko beam theory', () => {
  // Tip shear load P on a cantilever L × h × t: δ = P L³/(3 E I) + P L/(κ G A), κ = 5/6 (rectangular section).
  const E = 200e9; const nu = 0.3; const L = 2; const h = 0.2; const t = 0.05; const P = 1000;
  const I = (t * h ** 3) / 12;
  const G = E / (2 * (1 + nu));
  const exact = (P * L ** 3) / (3 * E * I) + (P * L) / ((5 / 6) * G * h * t);
  const run = (type, nx, ny, q4) => {
    const mesh = rectangleMesh({ lx: L, ly: h, nx, ny, type });
    const r = solvePlane({
      mesh, E, nu, thickness: t, fix: [{ nodes: mesh.nodeSets.left, ux: 0, uy: 0 }],
      edgeLoads: [{ edges: mesh.edges.right, ty: -P / (h * t) }], q4,
    });
    close(r.reactionSum.fy, P, 1e-9, 'vertical reaction');
    close(r.appliedSum.fy, -P, 1e-12, 'applied load');
    const tip = nodesWhere(mesh, (x, y) => Math.abs(x - L) < 1e-9 && Math.abs(y - h / 2) < 1e-9)[0];
    return -r.u[2 * tip + 1];
  };
  close(run('Q4', 20, 2, 'incompatible'), exact, 0.015, 'Q4 incompatible');
  close(run('Q8', 10, 2), exact, 0.01, 'Q8');
  close(run('T6', 10, 2), exact, 0.015, 'T6');
  // The plain bilinear element is far too stiff on the same mesh (shear locking) — the reason for the default.
  const locked = run('Q4', 20, 2, 'full');
  assert.ok(locked < 0.95 * exact, `fully integrated Q4 should lock, got ratio ${locked / exact}`);
});

test('T3 convergence: strain-energy error falls at second order (Lamé cylinder, smooth solution)', () => {
  const E = 200e9; const nu = 0.3; const ri = 1; const ro = 2; const p = 1e6;
  // Exact work of the pressure per unit thickness on the quarter ring: U = ½ p u(ri) (π ri / 2), plane strain.
  const uRi = (((1 + nu) / E) * (p * ri * ri)) / (ro * ro - ri * ri) * ((1 - 2 * nu) * ri + (ro * ro) / ri);
  const Uex = 0.5 * p * uRi * ((Math.PI * ri) / 2);
  const run = (n) => {
    const mesh = annulusMesh({ ri, ro, nr: n, ntheta: 2 * n, type: 'T3' });
    const r = solvePlane({ mesh, E, nu, mode: 'strain', fix: [{ nodes: mesh.nodeSets.start, uy: 0 }, { nodes: mesh.nodeSets.end, ux: 0 }], edgeLoads: [{ edges: mesh.edges.inner, pressure: p }] });
    return Math.abs(Uex - r.strainEnergy) / Uex;
  };
  const e1 = run(4); const e2 = run(8); const e3 = run(16);
  const o1 = Math.log2(e1 / e2); const o2 = Math.log2(e2 / e3);
  assert.ok(e3 < 0.005, `fine-mesh error ${e3}`);
  assert.ok(o1 > 1.7 && o1 < 2.3 && o2 > 1.8 && o2 < 2.2, `observed orders ${o1}, ${o2}`);
});

test('Kirsch problem: stress concentration at a circular hole tends to 3 (plate 20 hole radii wide)', () => {
  const E = 210e9; const nu = 0.3; const s0 = 100e6;
  const mesh = plateWithHoleMesh({ a: 1, w: 20, h: 20, nr: 16, ntheta: 8, grading: 40, type: 'Q8' });
  const r = solvePlane({
    mesh, E, nu, mode: 'stress',
    fix: [{ nodes: mesh.nodeSets.symX, uy: 0 }, { nodes: mesh.nodeSets.symY, ux: 0 }],
    edgeLoads: [{ edges: mesh.edges.right, tx: s0 }],
  });
  // Peak hoop stress at (0, a); Howland's finite-width correction for d/W = 0.05 is below 1 %.
  const top = nodesWhere(mesh, (x, y) => Math.abs(x) < 1e-9 && Math.abs(y - 1) < 1e-9)[0];
  close(r.stress.sxx[top] / s0, 3.0, 0.02, 'Kt');
  // Compressive hoop stress −σ0 at (a, 0).
  const side = nodesWhere(mesh, (x, y) => Math.abs(y) < 1e-9 && Math.abs(x - 1) < 1e-9)[0];
  close(r.stress.syy[side] / s0, -1.0, 0.05, 'hoop at 0°');
  close(r.reactionSum.fx, -s0 * 20, 1e-8, 'equilibrium');
});

test('Lamé thick cylinder: plane-strain annulus and axisymmetric tube wall', () => {
  const E = 200e9; const nu = 0.3; const ri = 0.1; const ro = 0.2; const p = 50e6;
  const hoop = (r) => ((p * ri * ri) / (ro * ro - ri * ri)) * (1 + (ro * ro) / (r * r));
  const radial = (r) => ((p * ri * ri) / (ro * ro - ri * ri)) * (1 - (ro * ro) / (r * r));
  // Plane strain radial displacement: u = (1+ν)/E · p ri²/(ro²−ri²) · [(1−2ν) r + ro²/r]
  const uR = (r) => (((1 + nu) / E) * (p * ri * ri)) / (ro * ro - ri * ri) * ((1 - 2 * nu) * r + (ro * ro) / r);
  const ann = annulusMesh({ ri, ro, nr: 6, ntheta: 12, type: 'Q8' });
  const r1 = solvePlane({
    mesh: ann, E, nu, mode: 'strain',
    fix: [{ nodes: ann.nodeSets.start, uy: 0 }, { nodes: ann.nodeSets.end, ux: 0 }],
    edgeLoads: [{ edges: ann.edges.inner, pressure: p }],
  });
  const ni = nodesWhere(ann, (x, y) => Math.abs(y) < 1e-12 && Math.abs(x - ri) < 1e-12)[0];
  const no = nodesWhere(ann, (x, y) => Math.abs(y) < 1e-12 && Math.abs(x - ro) < 1e-12)[0];
  close(r1.stress.syy[ni], hoop(ri), 0.01, 'hoop inner (plane strain)');
  close(r1.stress.syy[no], hoop(ro), 0.01, 'hoop outer');
  close(r1.stress.sxx[ni], radial(ri), 0.03, 'radial inner = -p');
  close(r1.u[2 * ni], uR(ri), 1e-3, 'u(ri)');
  close(r1.stress.szz[ni], nu * (hoop(ri) + radial(ri)), 0.02, 'axial plane-strain stress');

  // Axisymmetric: long tube with axial strain suppressed (plane strain along z).
  for (const type of ['Q4', 'Q8', 'T6']) {
    const tube = tubeWallMesh({ ri, t: ro - ri, length: 0.05, nr: type === 'Q4' ? 24 : 8, nz: 2, type });
    const r2 = solvePlane({
      mesh: tube, E, nu, mode: 'axisym',
      fix: [{ nodes: tube.nodeSets.bottom, uy: 0 }, { nodes: tube.nodeSets.top, uy: 0 }],
      edgeLoads: [{ edges: tube.edges.inner, pressure: p }],
    });
    const a = nodesWhere(tube, (x, y) => Math.abs(y) < 1e-12 && Math.abs(x - ri) < 1e-12)[0];
    close(r2.u[2 * a], uR(ri), type === 'Q4' ? 3e-3 : 1e-3, `axisym u(ri) ${type}`);
    close(r2.stress.szz[a], hoop(ri), type === 'Q4' ? 0.03 : 0.01, `axisym hoop ${type}`);
    // Total radial load = p · 2π ri · L
    close(r2.appliedSum.fx, p * 2 * Math.PI * ri * 0.05, 1e-9, 'pressure resultant');
  }
});

test('Zienkiewicz–Zhu indicator decreases with refinement and tracks the true energy error', () => {
  const E = 200e9; const nu = 0.3;
  const run = (n) => {
    const mesh = plateWithHoleMesh({ a: 1, w: 4, h: 4, nr: n, ntheta: n, grading: 4, type: 'Q4' });
    const r = solvePlane({ mesh, E, nu, fix: [{ nodes: mesh.nodeSets.symX, uy: 0 }, { nodes: mesh.nodeSets.symY, ux: 0 }], edgeLoads: [{ edges: mesh.edges.right, tx: 1e6 }] });
    return { zz: zzError({ mesh, result: r, E, nu }), U: r.strainEnergy };
  };
  const c = run(4); const m = run(8); const f = run(16);
  assert.ok(m.zz.eta < 0.62 * c.zz.eta && f.zz.eta < 0.62 * m.zz.eta, `eta ${c.zz.eta}, ${m.zz.eta}, ${f.zz.eta}`);
  // Displacement method: energy increases to the exact value; the true error energy on the middle mesh,
  // estimated by Richardson extrapolation, agrees with the ZZ estimate within a factor of about 2.
  const Uex = f.U + (f.U - m.U) / 3;
  const trueErr = Math.sqrt(2 * (Uex - m.U));
  assert.ok(m.zz.errorNorm > 0.4 * trueErr && m.zz.errorNorm < 2.5 * trueErr, `zz ${m.zz.errorNorm} vs true ${trueErr}`);
  assert.equal(m.zz.elemError.length, 2 * 8 * 8);
});

test('steady heat conduction: hollow cylinder (logarithmic profile), convective slab, internal generation', () => {
  const k = 45; const ri = 0.1; const ro = 0.16; const Ti = 150; const To = 20; const L = 0.2;
  const tube = tubeWallMesh({ ri, t: ro - ri, length: L, nr: 12, nz: 2, type: 'Q8' });
  const r = solveHeat2d({ mesh: tube, k, axisym: true, fixT: [{ nodes: tube.nodeSets.inner, T: Ti }, { nodes: tube.nodeSets.outer, T: To }] });
  // T(r) = Ti + (To − Ti) ln(r/ri)/ln(ro/ri);  Q = 2π k L (Ti − To)/ln(ro/ri)
  const mid = nodesWhere(tube, (x, y) => Math.abs(x - 0.13) < 1e-9 && Math.abs(y) < 1e-9)[0];
  close(r.T[mid], Ti + ((To - Ti) * Math.log(0.13 / ri)) / Math.log(ro / ri), 1e-5, 'T(mid)');
  const Q = (2 * Math.PI * k * L * (Ti - To)) / Math.log(ro / ri);
  // Net heat through all prescribed nodes vanishes (what enters at the bore leaves at the outer surface).
  assert.ok(Math.abs(r.heatIn) < 1e-6 * Q);
  // Radial flow from the centroid flux of the innermost element: q_r · 2π r L at its mid-radius.
  const qin = r.flux.qx[0] * 2 * Math.PI * L * (ri + (ro - ri) / 24);
  close(qin, Q, 2e-3, 'radial heat flow');

  // Slab with convection on both faces: q = (T1 − T2)/(1/h1 + L/k + 1/h2)
  const slab = rectangleMesh({ lx: 0.3, ly: 0.1, nx: 6, ny: 1, type: 'Q4' });
  const h1 = 25; const h2 = 400; const kk = 1.4;
  const s = solveHeat2d({ mesh: slab, k: kk, convection: [{ edges: slab.edges.left, h: h1, Tinf: 800 }, { edges: slab.edges.right, h: h2, Tinf: 20 }] });
  const q = (800 - 20) / (1 / h1 + 0.3 / kk + 1 / h2);
  close(s.flux.qx[0], q, 1e-9, 'slab flux');
  const l = nodesWhere(slab, (x) => x < 1e-12)[0];
  close(s.T[l], 800 - q / h1, 1e-9, 'hot-face temperature');
  // Uniform internal generation with both faces held at 0: T_max = Q L²/(8k)
  const gm = rectangleMesh({ lx: 0.3, ly: 0.1, nx: 8, ny: 1, type: 'Q8' });
  const g = solveHeat2d({ mesh: gm, k: kk, source: 5000, fixT: [{ nodes: gm.nodeSets.left, T: 0 }, { nodes: gm.nodeSets.right, T: 0 }] });
  close(g.maxT, (5000 * 0.3 * 0.3) / (8 * kk), 1e-9, 'parabolic profile');
  close(-g.heatIn, 5000 * 0.3 * 0.1, 1e-9, 'generated heat leaves through the faces');
});

const plateEdgesNodes = (mesh) => {
  const s = new Set();
  for (const ed of boundaryEdges(mesh)) for (const n of ed) s.add(n);
  return Int32Array.from(s);
};

test('MITC4 plate patch test: constant curvature reproduced on a distorted mesh', () => {
  const E = 70e9; const nu = 0.25; const t = 0.02;
  // w = (c1 x² + c2 y² + c3 x y)/2 with zero transverse shear: ry = −w,x, rx = w,y
  const c1 = 0.02; const c2 = -0.03; const c3 = 0.015;
  const w = (x, y) => 0.5 * (c1 * x * x + c2 * y * y + c3 * x * y);
  const rx = (x, y) => c2 * y + 0.5 * c3 * x;
  const ry = (x, y) => -(c1 * x + 0.5 * c3 * y);
  const mesh = rectangleMesh({ lx: 1, ly: 0.8, nx: 4, ny: 4, type: 'Q4', distort: wobble(1, 0.8, 0.06) });
  const bn = plateEdgesNodes(mesh);
  const fix = [...bn].map((n) => { const x = mesh.nodes[2 * n]; const y = mesh.nodes[2 * n + 1]; return { nodes: [n], w: w(x, y), rx: rx(x, y), ry: ry(x, y) }; });
  const r = solvePlate({ mesh, E, nu, t, fix });
  const D = (E * t ** 3) / (12 * (1 - nu * nu));
  let err = 0;
  for (let n = 0; n < mesh.nn; n++) err = Math.max(err, Math.abs(r.w[n] - w(mesh.nodes[2 * n], mesh.nodes[2 * n + 1])));
  assert.ok(err < 1e-9, `w error ${err}`);
  // κx = ∂ry/∂x = −c1, κy = −∂rx/∂y = −c2
  for (let n = 0; n < mesh.nn; n++) {
    close(r.moments.mx[n], D * (-c1 - nu * c2), 1e-6, 'mx');
    close(r.moments.my[n], D * (-c2 - nu * c1), 1e-6, 'my');
  }
});

test('Mindlin plate vs Navier/Timoshenko: w = 0.00406 q a⁴/D (simply supported), 0.00126 q a⁴/D (clamped)', () => {
  const E = 210e9; const nu = 0.3; const a = 1; const t = 0.01; const q = 1000;
  const D = (E * t ** 3) / (12 * (1 - nu * nu));
  const mesh = rectangleMesh({ lx: a, ly: a, nx: 24, ny: 24, type: 'Q4' });
  const edge = plateEdgesNodes(mesh);
  const centre = nodesWhere(mesh, (x, y) => Math.abs(x - 0.5) < 1e-9 && Math.abs(y - 0.5) < 1e-9)[0];
  // Hard simple support: w = 0 and the rotation about the edge normal suppressed.
  const ss = solvePlate({
    mesh, E, nu, t, pressure: q,
    fix: [{ nodes: edge, w: 0 }, { nodes: mesh.nodeSets.left, rx: 0 }, { nodes: mesh.nodeSets.right, rx: 0 }, { nodes: mesh.nodeSets.bottom, ry: 0 }, { nodes: mesh.nodeSets.top, ry: 0 }],
  });
  close((ss.w[centre] * D) / (q * a ** 4), 0.00406, 0.01, 'alpha SS');
  close(ss.moments.mx[centre] / (q * a * a), 0.0479, 0.015, 'beta SS');
  close(ss.reactionSum, -q * a * a, 1e-8, 'vertical equilibrium');
  const cl = solvePlate({ mesh, E, nu, t, pressure: q, fix: [{ nodes: edge, w: 0, rx: 0, ry: 0 }] });
  close((cl.w[centre] * D) / (q * a ** 4), 0.00126, 0.015, 'alpha clamped');
  // No shear locking: the same mesh at t/a = 1/1000 keeps the thin-plate answer.
  const thin = solvePlate({ mesh, E, nu, t: 0.001, pressure: q, fix: [{ nodes: edge, w: 0, rx: 0, ry: 0 }] });
  close((thin.w[centre] * (D / 1000)) / (q * a ** 4), 0.00126, 0.015, 'alpha clamped, very thin');
});

test('plate buckling k = 4 and fundamental frequency of a simply supported square plate', () => {
  const E = 210e9; const nu = 0.3; const a = 1; const t = 0.01; const rho = 7850;
  const D = (E * t ** 3) / (12 * (1 - nu * nu));
  const mesh = rectangleMesh({ lx: a, ly: a, nx: 16, ny: 16, type: 'Q4' });
  const fix = [{ nodes: plateEdgesNodes(mesh), w: 0 }, { nodes: mesh.nodeSets.left, rx: 0 }, { nodes: mesh.nodeSets.right, rx: 0 }, { nodes: mesh.nodeSets.bottom, ry: 0 }, { nodes: mesh.nodeSets.top, ry: 0 }];
  const b = bucklePlate({ mesh, E, nu, t, Nx: 1, fix });
  close((b.factors[0] * a * a) / (Math.PI ** 2 * D), 4.0, 0.015, 'k uniaxial');
  // Pure shear, simply supported square: k = 9.34 (Timoshenko & Gere Table 9-10)
  const s = bucklePlate({ mesh, E, nu, t, Nxy: 1, fix });
  close((s.factors[0] * a * a) / (Math.PI ** 2 * D), 9.34, 0.03, 'k shear');
  // ω11 = π² (1/a² + 1/b²) √(D/(ρ t))
  const m = modalPlate({ mesh, E, nu, t, rho, fix, nModes: 3 });
  close(m.omega[0], 2 * Math.PI ** 2 * Math.sqrt(D / (rho * t)), 0.01, 'omega11');
  close(m.omega[1], 5 * Math.PI ** 2 * Math.sqrt(D / (rho * t)), 0.02, 'omega12');
  close(m.omega[2], m.omega[1], 1e-3, 'degenerate pair');
});

test('facet shell: an inclined flat plate reproduces the plate solution; membrane patch state is exact', () => {
  const E = 210e9; const nu = 0.3; const t = 0.012; const q = 5000;
  const m2 = rectangleMesh({ lx: 1.2, ly: 0.8, nx: 12, ny: 8, type: 'Q4' });
  const edge = plateEdgesNodes(m2);
  const ref = solvePlate({ mesh: m2, E, nu, t, pressure: q, fix: [{ nodes: edge, w: 0, rx: 0, ry: 0 }] });
  // Rotate the plate into a general orientation.
  const R = (() => {
    const a = 0.7; const b = -0.4;
    const ca = Math.cos(a); const sa = Math.sin(a); const cb = Math.cos(b); const sb = Math.sin(b);
    return [[cb, 0, sb], [sa * sb, ca, -sa * cb], [-ca * sb, sa, ca * cb]];
  })();
  const nodes = new Float64Array(3 * m2.nn);
  for (let n = 0; n < m2.nn; n++) for (let k = 0; k < 3; k++) nodes[3 * n + k] = R[k][0] * m2.nodes[2 * n] + R[k][1] * m2.nodes[2 * n + 1];
  const r = solveShell({ nodes, quads: m2.conn, E, nu, t, pressure: q, fix: [{ nodes: edge, ux: 0, uy: 0, uz: 0, rx: 0, ry: 0, rz: 0 }] });
  const c = nodesWhere(m2, (x, y) => Math.abs(x - 0.6) < 1e-9 && Math.abs(y - 0.4) < 1e-9)[0];
  const normal = [R[0][2], R[1][2], R[2][2]];
  const wn = r.u[6 * c] * normal[0] + r.u[6 * c + 1] * normal[1] + r.u[6 * c + 2] * normal[2];
  close(wn, ref.w[c], 1e-6, 'normal deflection');
  const F = q * 1.2 * 0.8;
  close(r.reactionSum.fx, -F * normal[0], 1e-7, 'reaction x');
  close(r.reactionSum.fz, -F * normal[2], 1e-7, 'reaction z');

  // In-plane uniform tension on the inclined plate: N = σ t everywhere.
  const s0 = 80e6;
  const dir = [R[0][0], R[1][0], R[2][0]];
  const right = m2.nodeSets.right; const left = m2.nodeSets.left;
  const loads = [];
  const dy = 0.8 / 8;
  for (const n of right) {
    const y = m2.nodes[2 * n + 1];
    const trib = y < 1e-9 || y > 0.8 - 1e-9 ? dy / 2 : dy;
    loads.push({ node: n, fx: s0 * t * trib * dir[0], fy: s0 * t * trib * dir[1], fz: s0 * t * trib * dir[2] });
  }
  // Restrain the left edge along the plate x-direction through a statically determinate set plus the normal.
  const fixM = [];
  for (const n of left) fixM.push({ nodes: [n], ux: 0, uy: 0, uz: 0 });
  const rm = solveShell({ nodes, quads: m2.conn, E, nu: 0, t, nodalLoads: loads, fix: fixM });
  for (let e = 0; e < rm.elem.nx.length; e++) close(rm.elem.nx[e], s0 * t, 1e-6, 'membrane force');
  close(rm.maxVonMises, s0, 1e-6, 'von Mises');
});

test('facet shell: Scordelis–Lo roof, mid-side vertical deflection 0.3024 (MacNeal & Harder 1985)', () => {
  const R = 25; const L = 50; const t = 0.25; const E = 4.32e8; const nu = 0; const phi = (40 * Math.PI) / 180;
  const n = 16;
  // Quarter model: x from the diaphragm (0) to mid-span (L/2), θ from the crown (0) to the free edge (40°).
  const m2 = mappedMesh({ ns: n, nt: n, type: 'Q4', map: (s, tt) => [s * (L / 2), tt * phi] });
  const nodes = new Float64Array(3 * m2.nn);
  for (let i = 0; i < m2.nn; i++) {
    const th = m2.nodes[2 * i + 1];
    nodes[3 * i] = m2.nodes[2 * i]; nodes[3 * i + 1] = R * Math.sin(th); nodes[3 * i + 2] = R * Math.cos(th);
  }
  const r = solveShell({
    nodes, quads: m2.conn, E, nu, t, surfaceLoad: [0, 0, -90],
    fix: [
      { nodes: m2.nodeSets.s0, uy: 0, uz: 0 }, // rigid diaphragm
      { nodes: m2.nodeSets.s1, ux: 0, ry: 0, rz: 0 }, // mid-span symmetry
      { nodes: m2.nodeSets.t0, uy: 0, rx: 0, rz: 0 }, // crown symmetry
    ],
  });
  const p = nodesWhere(m2, (x, th) => Math.abs(x - L / 2) < 1e-9 && Math.abs(th - phi) < 1e-9)[0];
  close(-r.u[6 * p + 2], 0.3024, 0.04, 'free-edge mid-span deflection');
  // Total load = 90 × (quarter surface area)
  close(r.appliedSum.fz, -90 * (L / 2) * R * phi, 2e-3, 'load resultant (faceted area)');
  close(r.reactionSum.fz, -r.appliedSum.fz, 1e-8, 'equilibrium');
});

test('stiffened panel shell model: equilibrium and stiffening effect', () => {
  const E = 210e9; const nu = 0.3;
  const geo = { a: 2.4, b: 2.4, hw: 0.2, tp: 0.012, tw: 0.01, nx: 8, nyPerBay: 2, nWeb: 2 };
  const run = (nStiff, bf = 0, tf = 0) => {
    const m = stiffenedPanelMesh({ ...geo, nStiff, bf, tf });
    const all = [...m.nodeSets.x0, ...m.nodeSets.x1, ...m.nodeSets.y0, ...m.nodeSets.y1];
    const r = solveShell({ nodes: m.nodes, quads: m.quads, E, nu, t: m.thickness, fix: [{ nodes: all, ux: 0, uy: 0, uz: 0 }], pressure: Array.from({ length: m.ne }, (_, e) => (e < m.elemSets.plate.length ? 2e4 : 0)) });
    close(r.reactionSum.fz, -2e4 * 2.4 * 2.4, 1e-7, 'reaction');
    return r.maxDisp;
  };
  const d0 = run(0); const d2 = run(2); const dT = run(2, 0.1, 0.012);
  assert.ok(d2 < 0.25 * d0, `stiffeners reduce deflection: ${d0} -> ${d2}`);
  assert.ok(dT < d2, 'flanges add stiffness');
});

test('input validation', () => {
  const mesh = rectangleMesh({ lx: 1, ly: 1, nx: 2, ny: 2 });
  assert.throws(() => solvePlane({ mesh, E: -1, nu: 0.3, fix: [] }), RangeError);
  assert.throws(() => solvePlane({ mesh, E: 1, nu: 0.5, fix: [] }), RangeError);
  assert.throws(() => solvePlane({ mesh, E: 1, nu: 0.3, fix: [] }), /no supports/);
  assert.throws(() => solvePlate({ mesh: rectangleMesh({ lx: 1, ly: 1, nx: 2, ny: 2, type: 'T3' }), E: 1, nu: 0.3, t: 0.1, fix: [] }), RangeError);
  assert.throws(() => edgesWhere({}, () => true), RangeError);
});
