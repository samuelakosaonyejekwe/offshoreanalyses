import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  plateRigidity, navierPlate, ritzPlate, clampedPlateCoefficients, ritzPlateBuckling, ritzPlateFrequencies,
  largeDeflectionPlate, plateBucklingCoefficient, plateSlenderness, effectiveWidthRatio, columnReduction,
  plateUltimateStrength, stiffenerWithPlate, stiffenerColumnBuckling, stiffenerTripping, orthotropicPlate,
  grillageBeamGrid, hullGirderSection, hullGirderUltimate, deckPatchLoad,
} from '../../src/engine/structural/plates.js';
import { rectangleMesh, solvePlate, boundaryEdges, nodesWhere } from '../../src/engine/structural/fem2d.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.abs(b), `${msg || ''} got ${a}, expected ${b} (rel ${rel})`);
const steel = { E: 210e9, nu: 0.3 };

test('Navier: square simply supported plate, w = 0.00406 q a⁴/D, M = 0.0479 q a² (Timoshenko Table 8)', () => {
  const r = navierPlate({ a: 2, b: 2, t: 0.02, ...steel, load: { type: 'uniform', q: 10e3 } });
  close(r.alpha, 0.00406, 2e-3, 'alpha');
  close(r.betaX, 0.0479, 2e-3, 'beta');
  close(r.betaY, r.betaX, 1e-12);
  close(r.D, (210e9 * 0.02 ** 3) / (12 * 0.91), 1e-12);
  assert.ok(Math.abs(r.mxy) < 1e-9 * r.mx && Math.abs(r.qx) < 1e-6 * 10e3);
  // b/a = 2: α = 0.01013, βx (short span) = 0.1017, βy = 0.0464 (same table)
  const r2 = navierPlate({ a: 1, b: 2, t: 0.01, ...steel, load: { type: 'uniform', q: 1 } });
  close(r2.alpha, 0.01013, 2e-3); close(r2.betaX, 0.1017, 2e-3); close(r2.betaY, 0.0464, 3e-3);
  // Long plate tends to the strip solution 5 q a⁴ / (384 D)
  const strip = navierPlate({ a: 1, b: 12, t: 0.01, ...steel, load: { type: 'uniform', q: 1 } });
  close(strip.alpha, 5 / 384, 1e-3);
});

test('Navier: central point load on a square plate, w = 0.01160 P a²/D (Timoshenko §34); patch → point limit', () => {
  const a = 1.5; const P = 5000;
  const D = plateRigidity({ ...steel, t: 0.012 });
  const r = navierPlate({ a, b: a, t: 0.012, ...steel, load: { type: 'point', P, x0: a / 2, y0: a / 2 }, nTerms: 199 });
  close((r.w * D) / (P * a * a), 0.0116, 2e-3, 'point-load coefficient');
  const p = navierPlate({ a, b: a, t: 0.012, ...steel, load: { type: 'patch', P, x0: a / 2, y0: a / 2, u: 0.01, v: 0.01 }, nTerms: 199 });
  close(p.w, r.w, 2e-3, 'small patch ≈ point load');
  // A patch covering the whole plate is the uniform load.
  const full = navierPlate({ a, b: a, t: 0.012, ...steel, load: { type: 'patch', P: 1000 * a * a, x0: a / 2, y0: a / 2, u: a, v: a } });
  const uni = navierPlate({ a, b: a, t: 0.012, ...steel, load: { type: 'uniform', q: 1000 } });
  close(full.w, uni.w, 1e-10); close(full.mx, uni.mx, 1e-10);
  // Reciprocity (Maxwell–Betti): deflection at B due to load at A equals deflection at A due to load at B.
  const wAB = navierPlate({ a, b: a, t: 0.012, ...steel, load: { type: 'point', P, x0: 0.3, y0: 0.4 }, x: 1.1, y: 0.9 }).w;
  const wBA = navierPlate({ a, b: a, t: 0.012, ...steel, load: { type: 'point', P, x0: 1.1, y0: 0.9 }, x: 0.3, y: 0.4 }).w;
  close(wAB, wBA, 1e-12);
});

test('Ritz plate: reproduces Navier for simple supports and the classical clamped-plate coefficients', () => {
  const base = { a: 1.2, b: 0.9, t: 0.01, ...steel, load: { type: 'uniform', q: 2000 } };
  const n = navierPlate(base);
  const r = ritzPlate(base);
  close(r.w, n.w, 1e-5, 'w'); close(r.mx, n.mx, 1e-3, 'mx'); close(r.my, n.my, 1e-3, 'my');
  const pn = navierPlate({ ...base, load: { type: 'patch', P: 500, x0: 0.4, y0: 0.3, u: 0.3, v: 0.2 }, x: 0.5, y: 0.4, nTerms: 199 });
  const pr = ritzPlate({ ...base, load: { type: 'patch', P: 500, x0: 0.4, y0: 0.3, u: 0.3, v: 0.2 }, x: 0.5, y: 0.4, nBasis: 12 });
  close(pr.w, pn.w, 2e-3, 'patch load');
  // All edges clamped, square: w = 0.00126 q a⁴/D, edge moment −0.0513 q a², centre moment 0.0231 q a²
  // (Timoshenko & Woinowsky-Krieger Table 35).
  const c = clampedPlateCoefficients({ aspect: 1 });
  close(c.alpha, 0.00126, 5e-3, 'alpha');
  close(c.betaEdgeLong, -0.0513, 0.01, 'edge moment');
  close(c.betaCentreShortSpan, 0.0231, 0.01, 'centre moment');
  // b/a = 2: α = 0.00254, edge −0.0829, centre 0.0412 / 0.0158
  const c2 = clampedPlateCoefficients({ aspect: 2 });
  close(c2.alpha, 0.00254, 5e-3); close(c2.betaEdgeLong, -0.0829, 0.01); close(c2.betaCentreShortSpan, 0.0412, 0.01);
  close(c2.betaCentreLongSpan, 0.0158, 0.02); close(c2.betaEdgeShort, -0.0571, 0.01);
  // Same coefficients whichever side is the long one.
  const c05 = clampedPlateCoefficients({ aspect: 0.5 });
  close(c05.alpha, c2.alpha, 1e-9); close(c05.betaEdgeLong, c2.betaEdgeLong, 1e-9);
});

test('Ritz plate agrees with the MITC4 finite-element plate for mixed edges (clamped–free–simply supported)', () => {
  const a = 1; const b = 0.8; const t = 0.008; const q = 3000;
  const edges = { x0: 'C', x1: 'F', y0: 'S', y1: 'S' };
  const r = ritzPlate({ a, b, t, ...steel, edges, load: { type: 'uniform', q }, x: a, y: b / 2, nBasis: 10 });
  const mesh = rectangleMesh({ lx: a, ly: b, nx: 30, ny: 24 });
  const fe = solvePlate({
    mesh, ...steel, t, pressure: q,
    fix: [{ nodes: mesh.nodeSets.left, w: 0, rx: 0, ry: 0 }, { nodes: mesh.nodeSets.bottom, w: 0, ry: 0 }, { nodes: mesh.nodeSets.top, w: 0, ry: 0 }],
  });
  const n = nodesWhere(mesh, (x, y) => Math.abs(x - a) < 1e-9 && Math.abs(y - b / 2) < 1e-9)[0];
  close(r.w, fe.w[n], 0.01, 'free-edge mid-point deflection');
  assert.ok(boundaryEdges(mesh).length > 0);
  assert.throws(() => ritzPlate({ a, b, t, ...steel, edges: { x0: 'F', x1: 'F', y0: 'F', y1: 'F' }, load: { type: 'uniform', q } }), /supported/);
});

test('plate buckling coefficients: k = 4 (uniaxial), biaxial k = 2, shear 9.34, clamped 10.07 / 14.71', () => {
  close(plateBucklingCoefficient({ aspect: 1 }).k, 4, 1e-12);
  close(plateBucklingCoefficient({ aspect: 3 }).k, 4, 1e-12);
  assert.equal(plateBucklingCoefficient({ aspect: 3 }).m, 3);
  // a/b = √2 is the transition between one and two half-waves: k = 4.5
  close(plateBucklingCoefficient({ aspect: Math.SQRT2 }).k, 4.5, 1e-12);
  close(plateBucklingCoefficient({ aspect: 0.5 }).k, (0.5 + 2) ** 2, 1e-12);
  close(plateBucklingCoefficient({ aspect: 1, loading: 'biaxial', ratio: 1 }).k, 2, 1e-12);
  // Transverse tension equal to the compression: k = 8 for a square plate (m = 2, n = 1: 25/3.. minimised)
  const kt = plateBucklingCoefficient({ aspect: 1, loading: 'biaxial', ratio: -1 }).k;
  close(kt, Math.min(...[1, 2, 3, 4].map((m) => (m * m + 1) ** 2 / (m * m - 1) || Infinity).filter(Number.isFinite)), 1e-12);
  close(plateBucklingCoefficient({ aspect: 1, loading: 'shear' }).k, 9.34, 1e-12);
  close(plateBucklingCoefficient({ aspect: 4, loading: 'shear' }).k, 5.34 + 4 / 16, 1e-12);
  // Ritz buckling reproduces the exact simply supported values …
  const b = { a: 1, b: 1, t: 0.01, ...steel };
  close(ritzPlateBuckling({ ...b, sigmaX: 1 }).k, 4, 1e-4, 'Ritz k = 4');
  close(ritzPlateBuckling({ ...b, sigmaX: 1, sigmaY: 1 }).k, 2, 1e-4, 'Ritz biaxial');
  close(ritzPlateBuckling({ ...b, tau: 1 }).k, 9.34, 3e-3, 'Ritz shear (Timoshenko & Gere Table 9-10: 9.34)');
  close(ritzPlateBuckling({ ...b, a: 2, sigmaX: 1 }).k, 4, 1e-3, 'two half-waves');
  // … and the classical clamped results: uniaxial square 10.07 (Timoshenko & Gere Table 9-15), shear 14.71
  // (Budiansky & Connor 1948).
  close(plateBucklingCoefficient({ aspect: 1, edges: 'clamped' }).k, 10.07, 5e-3, 'clamped uniaxial');
  close(plateBucklingCoefficient({ aspect: 1, loading: 'shear', edges: 'clamped' }).k, 14.71, 5e-3, 'clamped shear');
  // Loaded edges simply supported, unloaded edges clamped, long plate: k → 6.97
  close(ritzPlateBuckling({ a: 1.98, b: 1, t: 0.01, ...steel, edges: { x0: 'S', x1: 'S', y0: 'C', y1: 'C' }, sigmaX: 1, nBasis: 12 }).k, 6.97, 5e-3, 'SS loaded / clamped sides');
  assert.throws(() => plateBucklingCoefficient({ aspect: 50 }), RangeError);
});

test('Ritz plate frequencies: simply supported ω_mn = π² (m²/a² + n²/b²) √(D/ρt); clamped square λ = 35.99', () => {
  const a = 1.2; const b = 0.8; const t = 0.01; const rho = 7850;
  const D = plateRigidity({ ...steel, t });
  const f = ritzPlateFrequencies({ a, b, t, ...steel, rho, nModes: 3 });
  const w = (m, n) => Math.PI ** 2 * (m * m / (a * a) + n * n / (b * b)) * Math.sqrt(D / (rho * t));
  close(f.omega[0], w(1, 1), 1e-6); close(f.omega[1], w(2, 1), 1e-4); close(f.omega[2], w(1, 2), 1e-3);
  // Fully clamped square plate: ω a² √(ρ t / D) = 35.99 (Leissa, "Vibration of Plates", Table 4.19)
  const c = ritzPlateFrequencies({ a: 1, b: 1, t, ...steel, rho, edges: { x0: 'C', x1: 'C', y0: 'C', y1: 'C' }, nModes: 1 });
  close(c.omega[0] * Math.sqrt((rho * t) / D), 35.99, 1e-3);
  // Added mass lowers the frequency by √(m/(m + ma)).
  const wet = ritzPlateFrequencies({ a, b, t, ...steel, rho, addedMass: rho * t * 3, nModes: 1 });
  close(wet.omega[0], f.omega[0] / 2, 1e-9);
});

test('large deflection: small-load limit equals the one-term linear solution; membrane stiffening; cubic residual', () => {
  const a = 1; const t = 0.005; const E = 210e9; const nu = 0.3;
  const D = plateRigidity({ E, nu, t });
  const small = largeDeflectionPlate({ a, b: a, t, E, nu, q: 1 });
  close(small.w0, (16 * 1) / (Math.PI ** 6 * D * 4), 1e-6, 'linear one-term limit 16q/(π⁶ D (1/a²+1/b²)²)');
  close(small.w0Linear / ((0.00406 * 1) / D), 1, 0.03, 'one-term linear vs exact Navier (2.5 % high)');
  const q = 50e3;
  const im = largeDeflectionPlate({ a, b: a, t, E, nu, q });
  const mv = largeDeflectionPlate({ a, b: a, t, E, nu, q, edges: 'movable' });
  assert.ok(im.w0 < mv.w0 && mv.w0 < im.w0Linear, 'restraint stiffens');
  // Hand calculation of the governing cubic for the square plate with immovable edges (ν = 0.3):
  //   4 D π⁴ w0/a⁴ + E t π⁴ w0³/a⁴ · [2.6/(8·0.91) + 1/8] = 16 q/π²
  const res = (4 * D * Math.PI ** 4 * im.w0) / a ** 4 + ((E * t * Math.PI ** 4 * im.w0 ** 3) / a ** 4) * (2.6 / 7.28 + 0.125) - (16 * q) / Math.PI ** 2;
  assert.ok(Math.abs(res) < 1e-9 * ((16 * q) / Math.PI ** 2), `cubic residual ${res}`);
  assert.ok(im.ratio > 1 && im.membraneX > 0 && im.bendingX > 0);
  // Membrane-dominated limit: w0 ∝ q^(1/3)
  const big1 = largeDeflectionPlate({ a, b: a, t: 0.001, E, nu, q: 1e5 });
  const big2 = largeDeflectionPlate({ a, b: a, t: 0.001, E, nu, q: 8e5 });
  close(big2.w0 / big1.w0, 2, 0.02, 'cube-root law');
});

test('slenderness, effective width and column curves (EN 1993-1-5 §4.4, EN 1993-1-1 §6.3.1)', () => {
  // b/t = 100, S355: σE = 0.9038 E (t/b)² = 18.98 MPa; σcr = 75.92 MPa; λp = √(355/75.92) = 2.162
  const s = plateSlenderness({ b: 1, t: 0.01, E: 210e9, nu: 0.3, fy: 355e6 });
  close(s.sigmaE, 18.98e6, 1e-3); close(s.lambdaP, 2.1624, 1e-3);
  // EN 1993-1-5 form: λp = (b/t) / (28.4 ε √k), ε = √(235/fy)
  close(s.lambdaP, 100 / (28.4 * Math.sqrt(235 / 355) * 2), 2e-3, 'EN closed form');
  close(s.beta, 100 * Math.sqrt(355e6 / 210e9), 1e-12);
  close(effectiveWidthRatio({ method: 'en1993-internal', lambdaP: 2.1624 }), (2.1624 - 0.22) / 2.1624 ** 2, 1e-12);
  assert.equal(effectiveWidthRatio({ method: 'en1993-internal', lambdaP: 0.673 }), 1);
  close(effectiveWidthRatio({ method: 'winter', lambdaP: 1.5 }), (1 - 0.22 / 1.5) / 1.5, 1e-12);
  // ψ = 1: EN internal formula is Winter's
  close(effectiveWidthRatio({ method: 'en1993-internal', lambdaP: 1.5, psi: 1 }), effectiveWidthRatio({ method: 'winter', lambdaP: 1.5 }), 1e-12);
  close(effectiveWidthRatio({ method: 'en1993-outstand', lambdaP: 1.2 }), (1.2 - 0.188) / 1.44, 1e-12);
  close(effectiveWidthRatio({ method: 'vonKarman', lambdaP: 2 }), 0.5, 1e-12);
  close(effectiveWidthRatio({ method: 'faulkner', beta: 2 }), 0.75, 1e-12);
  assert.equal(effectiveWidthRatio({ method: 'faulkner', beta: 0.8 }), 1);
  // Column curve b at λ = 1: Φ = 0.5 (1 + 0.34·0.8 + 1) = 1.136, χ = 1/(1.136 + √(1.136² − 1)) = 0.5970
  close(columnReduction({ lambda: 1, curve: 'b' }), 0.597, 1e-3);
  assert.equal(columnReduction({ lambda: 0.15 }), 1);
  // Perry–Robertson identity: (σE − σ)(fy − σ) = η σE σ with η = α(λ − 0.2)
  const lam = 1.3; const chi = columnReduction({ lambda: lam, curve: 'c' });
  const fy = 1; const sE = fy / (lam * lam); const sig = chi * fy;
  close((sE - sig) * (fy - sig), 0.49 * (lam - 0.2) * sE * sig, 1e-10, 'Perry–Robertson');
});

test('ultimate strength of unstiffened plates', () => {
  const r = plateUltimateStrength({ a: 3, b: 1, t: 0.01, E: 210e9, fy: 355e6 });
  close(r.k, 4, 1e-12);
  close(r.sigmaU, r.rho * 355e6, 1e-3, 'long plate: plate-type behaviour governs');
  // Wide short panel (a/b = 0.25): column-like behaviour reduces the strength below the plate value.
  const w = plateUltimateStrength({ a: 0.5, b: 2, t: 0.012, E: 210e9, fy: 355e6 });
  assert.ok(w.rhoC < w.rho && w.rhoC > 0);
  // Shear: stocky web reaches yield, τu = fy/√3
  const st = plateUltimateStrength({ a: 1, b: 1, t: 0.03, E: 210e9, fy: 355e6, tau: 100e6 });
  close(st.tauU, 355e6 / Math.sqrt(3), 1e-12);
  close(st.utilisation, 100e6 / st.tauU, 1e-12);
  // Slender web, rigid end post, λw ≥ 1.08: χw = 1.37/(0.7 + λw)
  const sl = plateUltimateStrength({ a: 2, b: 2, t: 0.008, E: 210e9, fy: 355e6, rigidEndPost: true });
  assert.ok(sl.lambdaW > 1.08);
  close(sl.chiW, 1.37 / (0.7 + sl.lambdaW), 1e-12);
  const fa = plateUltimateStrength({ a: 3, b: 1, t: 0.01, E: 210e9, fy: 355e6, method: 'faulkner' });
  const beta = 100 * Math.sqrt(355e6 / 210e9);
  close(fa.rho, 2 / beta - 1 / beta ** 2, 1e-12);
});

test('stiffened panel: section properties (hand calculation), column buckling and tripping', () => {
  // Plate 600 × 10 with flat bar 150 × 10: A = 6000 + 1500 = 7500 mm²; z̄ = 1500·80/7500 = 16 mm above plate
  // mid-plane; I = 600·10³/12 + 6000·16² + 10·150³/12 + 1500·64² = 50 000 + 1 536 000 + 2 812 500 + 6 144 000
  //   = 10 542 500 mm⁴
  const s = stiffenerWithPlate({ be: 0.6, tp: 0.01, hw: 0.15, tw: 0.01 });
  close(s.A, 7500e-6, 1e-12); close(s.zNA, 0.016, 1e-12); close(s.I, 10542500e-12, 1e-12);
  close(s.Zflange, 10542500e-12 / (0.155 - 0.016), 1e-12);
  const c = stiffenerColumnBuckling({ L: 3, s: 0.6, tp: 0.01, hw: 0.15, tw: 0.01, E: 210e9, fy: 355e6 });
  assert.ok(c.rho < 1 && c.be < 0.6);
  close(c.sigmaE, (Math.PI ** 2 * 210e9 * c.section.I) / (c.section.A * 9), 1e-12, 'Euler stress of effective column');
  close(c.chi, columnReduction({ lambda: c.lambda, curve: 'b' }), 1e-12);
  close(c.Pu, c.sigmaU * c.sectionFull.A, 1e-12);
  assert.ok(c.sigmaU < 355e6);
  // Stocky plating: fully effective, Euler stress from the gross section
  const full = stiffenerColumnBuckling({ L: 3, s: 0.3, tp: 0.02, hw: 0.15, tw: 0.01, E: 210e9, fy: 235e6 });
  assert.equal(full.rho, 1);
  // Tripping of a flat bar: for long spans σT → G It/Ip = G (tw/hw)²; here 80.77 GPa · (10/150)² = 359 MPa
  const t = stiffenerTripping({ L: 100, hw: 0.15, tw: 0.01, E: 210e9, nu: 0.3 });
  close(t.sigmaT, (210e9 / 2.6) * (0.01 / 0.15) ** 2, 1e-3);
  const tee = stiffenerTripping({ L: 3, hw: 0.3, tw: 0.01, bf: 0.1, tf: 0.015, E: 210e9, fy: 355e6 });
  const G = 210e9 / 2.6;
  const It = (0.3 * 0.01 ** 3 + 0.1 * 0.015 ** 3) / 3; const Ip = 0.3 ** 3 * 0.01 / 3 + 0.09 * 0.1 * 0.015 + 0.1 ** 3 * 0.015 / 12;
  const Iw = 0.015 * 0.1 ** 3 * 0.09 / 12;
  close(tee.sigmaT, (G * It + Math.PI ** 2 * 210e9 * Iw / 9) / Ip, 1e-12);
  assert.ok(tee.sigmaTinelastic <= 355e6);
});

test('grillage: orthotropic plate reduces to the isotropic Navier solution; beam grid vs beam theory', () => {
  const a = 2; const b = 1.5; const t = 0.01; const q = 1000;
  const D = plateRigidity({ ...steel, t });
  const iso = navierPlate({ a, b, t, ...steel, load: { type: 'uniform', q } });
  const o = orthotropicPlate({ a, b, Dx: D, Dy: D, H: D, q });
  close(o.wCentre, iso.w, 1e-4);
  // Isotropic buckling check: Nx,cr = k π² D / b² with k from the closed form
  close(o.NxCr, (plateBucklingCoefficient({ aspect: a / b }).k * Math.PI ** 2 * D) / (b * b), 1e-12);
  // Single girder crossing a very flexible transverse: simply supported beam with central load, δ = P L³/(48 EI)
  const g = grillageBeamGrid({ a: 6, b: 4, nGirders: 1, nBeams: 1, EIx: 2e7, EIy: 1e-3, pointLoads: [{ i: 1, j: 1, P: 1e4 }] });
  close(g.wMax, (1e4 * 6 ** 3) / (48 * 2e7), 1e-6, 'girder alone');
  close(g.girderMomentMax, (1e4 * 6) / 4, 1e-6, 'M = PL/4');
  close(g.reactionSum, -1e4, 1e-9);
  // Two crossing members sharing a central load: P splits in proportion to 48EI/L³ (compatibility at the node)
  const kx = (48 * 2e7) / 6 ** 3; const ky = (48 * 1e7) / 4 ** 3;
  const x = grillageBeamGrid({ a: 6, b: 4, nGirders: 1, nBeams: 1, EIx: 2e7, EIy: 1e7, pointLoads: [{ i: 1, j: 1, P: 1e4 }] });
  close(x.wMax, 1e4 / (kx + ky), 1e-9, 'load sharing');
  close(x.beamMomentMax, ((1e4 * ky) / (kx + ky)) * 4 / 4, 1e-9);
  // Clamped ends: δ = P L³/(192 EI)
  const c = grillageBeamGrid({ a: 6, b: 4, nGirders: 1, nBeams: 1, EIx: 2e7, EIy: 1e-3, pointLoads: [{ i: 1, j: 1, P: 1e4 }], support: 'clamped' });
  close(c.wMax, (1e4 * 6 ** 3) / (192 * 2e7), 1e-6);
  // Dense grid under uniform pressure approaches the orthotropic plate with H = 0 (torsionless members)
  const n = 11;
  const dense = grillageBeamGrid({ a: 6, b: 6, nGirders: n, nBeams: n, EIx: 1e6, EIy: 1e6, q: 5000 });
  const sp = 6 / (n + 1);
  const smeared = orthotropicPlate({ a: 6, b: 6, Dx: 1e6 / sp, Dy: 1e6 / sp, H: 0, q: 5000 });
  close(dense.wMax, smeared.wCentre, 0.03, 'dense grillage vs smeared plate');
});

test('hull girder: box-section modulus and progressive collapse between first yield and plastic moment', () => {
  // Idealised box 20 m wide × 10 m deep, 20 mm plating, lumped: deck and bottom 0.4 m² each, sides 2 × 0.2 m²
  // split into 10 strips. I = 2 · 0.4 · 5² + 2 · (0.02 · 10³/12) = 20 + 3.333 = 23.33 m⁴; Z = 4.667 m³
  const el = [{ area: 0.4, z: 10, fy: 355e6 }, { area: 0.4, z: 0, fy: 355e6 }];
  for (let i = 0; i < 10; i++) el.push({ area: 0.04, z: 0.5 + i, fy: 355e6, i0: (0.04 * 1) / 12 });
  const s = hullGirderSection({ elements: el });
  close(s.zNA, 5, 1e-12); close(s.I, 20 + 2 * (0.02 * 1000) / 12, 1e-9); close(s.Ztop, s.I / 5, 1e-12);
  // No buckling: ultimate = plastic moment = fy (0.4·5·2 + 2·0.2·2.5·... ) = fy · (4 + 0.4·2.5·... ) lumped
  const u = hullGirderUltimate({ elements: el, E: 210e9, maxStrainFactor: 30, nSteps: 300 });
  let Mp = 0;
  for (const e of el) Mp += e.area * 355e6 * Math.abs(e.z - 5);
  close(u.Mp, Mp, 1e-12);
  close(u.sagging.Mu, Mp, 5e-3, 'sagging reaches the plastic moment without buckling');
  close(u.hogging.Mu, u.sagging.Mu, 1e-9, 'symmetric section');
  close(u.MyieldFirst, 355e6 * s.Ztop, 1e-12);
  // Initial slope of the M–κ curve is E I (lumped strips: without own inertia of the strips)
  const EIlumped = 210e9 * (s.I - 10 * (0.04 / 12));
  close(u.sagging.moment[1] / u.sagging.curvature[1], EIlumped, 1e-6, 'elastic stiffness');
  // Deck buckling at 60 % of yield: sagging capacity drops, hogging (deck in tension) is hardly affected
  const el2 = el.map((e) => (e.z === 10 ? { ...e, sigmaU: 0.6 * 355e6 } : e));
  const u2 = hullGirderUltimate({ elements: el2, E: 210e9, maxStrainFactor: 30, nSteps: 300 });
  assert.ok(u2.sagging.Mu < 0.85 * u.sagging.Mu && u2.sagging.Mu > 0.6 * u.sagging.Mu);
  close(u2.hogging.Mu, u.hogging.Mu, 1e-6);
  // Exact value for the weakened section: force balance 0.6·0.4 fy (deck) + side compression = bottom + side tension
  assert.throws(() => hullGirderUltimate({ elements: [{ area: 1, z: 0, fy: 1, sigmaU: 2 }, { area: 1, z: 1, fy: 1 }], E: 1 }), RangeError);
});

test('deck patch load: consistent with the plate solvers; clamped panel governed by the edge moment', () => {
  const o = { a: 2.4, b: 0.7, t: 0.012, E: 210e9, P: 50e3, u: 0.3, v: 0.25, fy: 355e6 };
  const ss = deckPatchLoad(o);
  const n = navierPlate({ a: 2.4, b: 0.7, t: 0.012, ...steel, load: { type: 'patch', P: 50e3, x0: 1.2, y0: 0.35, u: 0.3, v: 0.25 }, nTerms: 199 });
  close(ss.w, n.w, 1e-12);
  close(ss.sigmaMax, (6 * Math.max(n.mx, n.my)) / 0.012 ** 2, 1e-12);
  const cl = deckPatchLoad({ ...o, edges: 'clamped' });
  assert.ok(cl.w < ss.w && cl.edgeMomentMax > 0 && cl.utilisation > 0);
  assert.throws(() => deckPatchLoad({ ...o, u: 5 }), RangeError);
});
