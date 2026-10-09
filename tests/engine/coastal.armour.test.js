import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  relativeDensity, HUDSON_KD, hudson, vanDerMeerRock, vanDerMeerDamage, vanDerMeerShallow,
  concreteUnitVanDerMeer, stabilityNumberSize, toeStabilityVanDerMeer, massToDn, dnToMass, armourLayer,
  filterRules, underlayerMass, rubbleMoundQuantities,
} from '../../src/engine/coastal/armour.js';

const g = 9.80665;
const rel = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg || ''} got ${a}, expected ${b}`);

test('Hudson formula: hand calculation and equivalence with Ns = (KD cot α)^(1/3)', () => {
  // H = 3 m, rough angular quarrystone non-breaking KD = 4, cot α = 2, ρs = 2650, ρw = 1025:
  // Δ = 1.5854, M50 = 2650·27/(4·1.5854³·2) = 2245 kg
  const delta = relativeDensity({ rhoS: 2650, rhoW: 1025 });
  rel(delta, 2650 / 1025 - 1, 1e-14);
  const r = hudson({ H: 3, KD: HUDSON_KD['quarrystone-rough-angular'].nonBreaking, cotAlpha: 2, rhoS: 2650 });
  rel(r.M50, (2650 * 27) / (4 * delta ** 3 * 2), 1e-12);
  assert.ok(Math.abs(r.M50 - 2245) < 3, `M50 = ${r.M50}`);
  rel(r.Ns, 2, 1e-12);
  rel(3 / (delta * r.Dn50), r.Ns, 1e-12);
  rel(r.W50, r.M50 * g, 1e-12);
  // Mass scales with H³ and inversely with KD.
  rel(hudson({ H: 6, KD: 4, cotAlpha: 2, rhoS: 2650 }).M50, 8 * r.M50, 1e-12);
  rel(hudson({ H: 3, KD: 8, cotAlpha: 2, rhoS: 2650 }).M50, r.M50 / 2, 1e-12);
  assert.equal(HUDSON_KD.tetrapod.breaking, 7.0);
  assert.throws(() => hudson({ H: 3, KD: 4, cotAlpha: 2, rhoS: 900 }), RangeError);
});

test('Van der Meer (1988): hand calculation of the plunging formula', () => {
  // Hs = 4 m, Tm = 8 s, cot α = 3, P = 0.4, S = 2, N = 3000, Δ = 1.6
  // som = 2π·4/(g·64) = 0.04004, ξm = (1/3)/√som = 1.666
  // ξcr = [6.2·0.4^0.31·√(1/3)]^(1/0.9) = 3.09 → plunging
  // Ns = 6.2·0.4^0.18·(2/√3000)^0.2·1.666^−0.5 = 2.10 → Dn50 = 4/(1.6·2.10) = 1.19 m
  const r = vanDerMeerRock({ Hs: 4, Tm: 8, cotAlpha: 3, P: 0.4, S: 2, N: 3000, delta: 1.6, rhoS: 2650 });
  const xi = 1 / 3 / Math.sqrt((2 * Math.PI * 4) / (g * 64));
  rel(r.xi, xi, 1e-12);
  rel(r.xiCr, (6.2 * 0.4 ** 0.31 * Math.sqrt(1 / 3)) ** (1 / 0.9), 1e-12);
  assert.equal(r.regime, 'plunging');
  const Ns = 6.2 * 0.4 ** 0.18 * (2 / Math.sqrt(3000)) ** 0.2 / Math.sqrt(xi);
  rel(r.Ns, Ns, 1e-12);
  assert.ok(Math.abs(r.Ns - 2.10) < 0.01 && Math.abs(r.Dn50 - 1.19) < 0.01, `Ns ${r.Ns}, Dn50 ${r.Dn50}`);
  rel(r.M50, 2650 * r.Dn50 ** 3, 1e-12);
});

test('Van der Meer: plunging and surging formulae intersect at ξcr; surging branch by hand; cot α >= 4 is plunging only', () => {
  const o = { Hs: 3, cotAlpha: 2, P: 0.5, S: 3, N: 2000, delta: 1.6 };
  const xiCr = (6.2 * 0.5 ** 0.31 * Math.sqrt(0.5)) ** (1 / 1.0);
  const pl = 6.2 * 0.5 ** 0.18 * (3 / Math.sqrt(2000)) ** 0.2 / Math.sqrt(xiCr);
  const su = 1.0 * 0.5 ** -0.13 * (3 / Math.sqrt(2000)) ** 0.2 * Math.sqrt(2) * xiCr ** 0.5;
  rel(pl, su, 1e-12, 'intersection');
  // Period giving ξm just below / above ξcr
  const Tat = (xi) => Math.sqrt((2 * Math.PI * o.Hs) / g) * xi * o.cotAlpha;
  const below = vanDerMeerRock({ ...o, Tm: Tat(xiCr * 0.999) });
  const above = vanDerMeerRock({ ...o, Tm: Tat(xiCr * 1.001) });
  assert.equal(below.regime, 'plunging');
  assert.equal(above.regime, 'surging');
  rel(below.Ns, above.Ns, 2e-3);
  const s = vanDerMeerRock({ ...o, Tm: Tat(5) });
  rel(s.Ns, 0.5 ** -0.13 * (3 / Math.sqrt(2000)) ** 0.2 * Math.sqrt(2) * 5 ** 0.5, 1e-10);
  // Stability is lowest at the transition.
  assert.ok(below.Ns < vanDerMeerRock({ ...o, Tm: Tat(2) }).Ns && above.Ns < s.Ns);
  const gentle = vanDerMeerRock({ ...o, cotAlpha: 4, Tm: 30 });
  assert.equal(gentle.regime, 'plunging');
  assert.throws(() => vanDerMeerRock({ ...o, Tm: 8, P: 0.05 }), RangeError);
});

test('Van der Meer damage is the inverse of the design formula and S grows with √N up to 7500 waves', () => {
  const o = { Hs: 4, Tm: 9, cotAlpha: 2.5, P: 0.4, N: 3000, delta: 1.62 };
  const d = vanDerMeerRock({ ...o, S: 5 });
  rel(vanDerMeerDamage({ ...o, Dn50: d.Dn50 }).S, 5, 1e-10);
  const s1 = vanDerMeerDamage({ ...o, N: 1000, Dn50: d.Dn50 }).S;
  const s4 = vanDerMeerDamage({ ...o, N: 4000, Dn50: d.Dn50 }).S;
  rel(s4 / s1, 2, 1e-12);
  rel(vanDerMeerDamage({ ...o, N: 20000, Dn50: d.Dn50 }).S, vanDerMeerDamage({ ...o, N: 7500, Dn50: d.Dn50 }).S, 1e-12);
  // Damage ∝ Hs^5 at constant ξ: 10 % more wave height at the same steepness gives 61 % more damage.
  const k = Math.sqrt(1.1);
  rel(vanDerMeerDamage({ ...o, Hs: 4.4, Tm: 9 * k, Dn50: d.Dn50 }).S / 5, 1.1 ** 5, 1e-10);
});

test('Van der Meer shallow-water modification (Van Gent et al. 2003): cpl = 8.4, cs = 1.3, Hs/H2% factor', () => {
  const o = { Hs: 2.5, H2: 3.2, Tm10: 9, cotAlpha: 3, P: 0.4, S: 2, N: 1000, delta: 1.6 };
  const r = vanDerMeerShallow(o);
  const xi = 1 / 3 / Math.sqrt((2 * Math.PI * 2.5) / (g * 81));
  rel(r.xiCr, ((8.4 / 1.3) * 0.4 ** 0.31 * Math.sqrt(1 / 3)) ** (1 / 0.9), 1e-12);
  assert.equal(r.regime, 'plunging');
  rel(r.Ns, 8.4 * 0.4 ** 0.18 * (2 / Math.sqrt(1000)) ** 0.2 * (2.5 / 3.2) / Math.sqrt(xi), 1e-12);
  // With the Rayleigh ratio H2%/Hs = 1.4 the shallow formula is the deep one with 8.4/1.4 = 6.0 ≈ 6.2.
  const deepLike = vanDerMeerShallow({ ...o, H2: 1.4 * 2.5 });
  const deep = vanDerMeerRock({ Hs: 2.5, Tm: 9, cotAlpha: 3, P: 0.4, S: 2, N: 1000, delta: 1.6 });
  rel(deepLike.Ns, deep.Ns * (8.4 / 1.4 / 6.2), 1e-10);
});

test('concrete armour units (Van der Meer 1988): no-damage limits and hand values', () => {
  const som = (2 * Math.PI * 5) / (g * 100);
  const c0 = concreteUnitVanDerMeer({ unit: 'cube', Hs: 5, Tm: 10, Nod: 0, N: 3000, delta: 1.34 });
  rel(c0.Ns, som ** -0.1, 1e-12);
  const c = concreteUnitVanDerMeer({ unit: 'cube', Hs: 5, Tm: 10, Nod: 0.5, N: 3000, delta: 1.34, rhoS: 2400 });
  rel(c.Ns, ((6.7 * 0.5 ** 0.4) / 3000 ** 0.3 + 1.0) * som ** -0.1, 1e-12);
  rel(c.M, 2400 * (5 / (1.34 * c.Ns)) ** 3, 1e-12);
  const t = concreteUnitVanDerMeer({ unit: 'tetrapod', Hs: 5, Tm: 10, Nod: 0.5, N: 3000, delta: 1.34 });
  rel(t.Ns, ((3.75 * Math.sqrt(0.5)) / 3000 ** 0.25 + 0.85) * som ** -0.2, 1e-12);
  assert.ok(t.Ns > c.Ns);
  assert.throws(() => concreteUnitVanDerMeer({ unit: 'dolos', Hs: 5, Tm: 10, Nod: 0, N: 1000, delta: 1.3 }), RangeError);
  const s = stabilityNumberSize({ Hs: 5, Ns: 2.7, delta: 1.34, rhoS: 2400 });
  rel(s.Dn, 5 / (1.34 * 2.7), 1e-12);
  rel(s.M, 2400 * s.volume, 1e-12);
});

test('toe stability (Van der Meer 1998): Hs/(Δ Dn50) = (2 + 6.2 (ht/h)^2.7) Nod^0.15', () => {
  // ht/h = 0.6, Nod = 2: Ns = (2 + 6.2·0.6^2.7)·2^0.15 = (2 + 1.561)·1.110 = 3.95
  const r = toeStabilityVanDerMeer({ Hs: 3, ht: 4.8, h: 8, delta: 1.6, Nod: 2, rhoS: 2650 });
  rel(r.Ns, (2 + 6.2 * 0.6 ** 2.7) * 2 ** 0.15, 1e-12);
  assert.ok(Math.abs(r.Ns - 3.95) < 0.01, `Ns = ${r.Ns}`);
  rel(r.Dn50, 3 / (1.6 * r.Ns), 1e-12);
  assert.equal(r.valid, true);
  // A deeper toe needs smaller rock; outside 0.4 < ht/h < 0.9 the result is flagged.
  assert.ok(toeStabilityVanDerMeer({ Hs: 3, ht: 6.5, h: 8, delta: 1.6 }).Dn50 < r.Dn50);
  assert.equal(toeStabilityVanDerMeer({ Hs: 3, ht: 2, h: 8, delta: 1.6 }).valid, false);
  assert.throws(() => toeStabilityVanDerMeer({ Hs: 3, ht: 9, h: 8, delta: 1.6 }), RangeError);
});

test('grading helpers: M50 <-> Dn50, layer thickness, placing density', () => {
  const { Dn50, D50 } = massToDn({ M50: 5000, rhoS: 2650 });
  rel(Dn50, (5000 / 2650) ** (1 / 3), 1e-12);
  rel(D50, Dn50 / 0.84, 1e-12);
  rel(dnToMass({ Dn50, rhoS: 2650 }), 5000, 1e-12);
  // Two layers, kt = 1.0, porosity 0.37: t = 2 Dn50, N = 2·0.63/Dn50², mass = N·M50
  const l = armourLayer({ Dn50, rhoS: 2650 });
  rel(l.thickness, 2 * Dn50, 1e-12);
  rel(l.unitsPerArea, 1.26 / Dn50 ** 2, 1e-12);
  rel(l.massPerArea, l.unitsPerArea * 5000, 1e-12);
  // Consistency: placed mass per area = thickness × (1 − porosity) × ρs
  rel(l.massPerArea, l.thickness * 0.63 * 2650, 1e-12);
  rel(underlayerMass({ M50armour: 5000 }), 500, 1e-12);
});

test('filter rules (Terzaghi retention and permeability, uniformity)', () => {
  const ok = filterRules({ D15f: 0.02, D50f: 0.04, D85b: 0.006, D15b: 0.002, D50b: 0.004 });
  rel(ok.retention.ratio, 0.02 / 0.006, 1e-12);
  rel(ok.permeability.ratio, 10, 1e-12);
  rel(ok.uniformity.ratio, 10, 1e-12);
  assert.equal(ok.ok, true);
  const washout = filterRules({ D15f: 0.05, D85b: 0.006, D15b: 0.002 });
  assert.equal(washout.retention.ok, false);
  assert.equal(washout.ok, false);
  const tight = filterRules({ D15f: 0.006, D85b: 0.006, D15b: 0.002 });
  assert.equal(tight.permeability.ok, false);
});

// Shoelace area of a polygon given as [x, z] pairs.
function polyArea(p) {
  let a = 0;
  for (let i = 0; i < p.length; i++) {
    const [x1, z1] = p[i];
    const [x2, z2] = p[(i + 1) % p.length];
    a += x1 * z2 - x2 * z1;
  }
  return Math.abs(a) / 2;
}

test('cross-section quantities agree with an independent polygon (shoelace) construction of the offset layers', () => {
  const B = 6;
  const H = 10;
  const cf = 2;
  const cr = 1.5;
  const t1 = 2.4;
  const t2 = 1.1;
  const q = rubbleMoundQuantities({
    crestWidth: B, height: H, cotFront: cf, cotRear: cr,
    layers: [{ name: 'armour', thickness: t1, porosity: 0.37, rhoS: 2650 }, { name: 'under', thickness: t2, porosity: 0.35, rhoS: 2650 }],
    core: { porosity: 0.3, rhoS: 2600 },
  });
  // Profile at inward offset d: slope lines shifted by d along their normals, crest lowered by d.
  const profile = (d) => {
    const zc = H - d;
    // front face (seaward, x decreasing to the left): outer line through (0, H) with dx/dz = −cf
    const xf = (z) => -cf * (H - z) + d * Math.sqrt(1 + cf * cf);
    const xr = (z) => B + cr * (H - z) - d * Math.sqrt(1 + cr * cr);
    return [[xf(0), 0], [xr(0), 0], [xr(zc), zc], [xf(zc), zc]];
  };
  const a0 = polyArea(profile(0));
  const a1 = polyArea(profile(t1));
  const a2 = polyArea(profile(t1 + t2));
  rel(q.totalArea, a0, 1e-12);
  rel(q.totalArea, B * H + 0.5 * H * H * (cf + cr), 1e-12);
  rel(q.layers[0].area, a0 - a1, 1e-12, 'armour');
  rel(q.layers[1].area, a1 - a2, 1e-12, 'underlayer');
  rel(q.layers[2].area, a2, 1e-12, 'core');
  rel(q.layers.reduce((s, l) => s + l.area, 0), q.totalArea, 1e-12);
  rel(q.layers[0].mass, (a0 - a1) * 0.63 * 2650, 1e-12);
  rel(q.baseWidth, B + H * (cf + cr), 1e-12);
  // Vertical-sided block: a layer of thickness t removes t from each side and from the top.
  const v = rubbleMoundQuantities({ crestWidth: 5, height: 4, cotFront: 0, cotRear: 0, layers: [{ name: 'skin', thickness: 0.5 }] });
  rel(v.layers[1].area, 4 * 3.5, 1e-12);
  assert.throws(() => rubbleMoundQuantities({ crestWidth: 1, height: 4, cotFront: 2, cotRear: 2, layers: [{ name: 'a', thickness: 3 }] }), RangeError);
});
