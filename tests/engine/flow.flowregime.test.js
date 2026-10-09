import test from 'node:test';
import assert from 'node:assert/strict';
import {
  stratifiedGeometry, taitelDuklerResidual, equilibriumLevel, taitelDuklerParameters, taitelDukler,
  taitelBarneaDuklerVertical, beggsBrillPattern, flowRegime, regimeMap, levelFromHoldup, stratifiedForces,
} from '../../src/engine/flow/flowregime.js';
import { G } from '../../src/engine/flow/driftflux.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-300), `${msg ?? ''} got ${a}, expected ${b}`);
// Air–water at atmospheric conditions in a 25 mm pipe: the system of Taitel & Dukler (1976) Fig. 4 / Mandhane data.
const aw25 = { D: 0.025, rhoL: 1000, rhoG: 1.2, muL: 1e-3, muG: 1.8e-5, sigma: 0.072 };
const DEG = Math.PI / 180;

test('stratified geometry: half-full pipe by hand and area closure', () => {
  const g = stratifiedGeometry(0.5);
  close(g.AL, Math.PI / 8, 1e-14); close(g.AG, Math.PI / 8, 1e-14);
  close(g.SL, Math.PI / 2, 1e-14); close(g.SG, Math.PI / 2, 1e-14);
  close(g.Si, 1, 1e-14); close(g.holdup, 0.5, 1e-14);
  close(g.DL, 4 * (Math.PI / 8) / (Math.PI / 2), 1e-14);
  close(g.DG, 4 * (Math.PI / 8) / (Math.PI / 2 + 1), 1e-14);
  for (const h of [0.05, 0.3, 0.77, 0.99]) {
    const q = stratifiedGeometry(h);
    close(q.AL + q.AG, Math.PI / 4, 1e-13);
    close(q.SL + q.SG, Math.PI, 1e-13);
    // dÃL/dh̃ by central difference
    const d = (stratifiedGeometry(h + 1e-6).AL - stratifiedGeometry(h - 1e-6).AL) / 2e-6;
    close(q.dAL, d, 1e-6);
    close(levelFromHoldup(q.holdup), h, 1e-9);
  }
});

test('Taitel–Dukler equilibrium level satisfies the dimensionless momentum balance (X, Y formulation)', () => {
  for (const X2 of [1e-4, 1e-2, 1, 100, 1e4]) {
    for (const Y of [-5, 0, 10, 1000]) {
      const h = equilibriumLevel({ X2, Y });
      assert.ok(h > 0 && h < 1);
      const L = X2 * (() => { const g = stratifiedGeometry(h); return (g.uL * g.DL) ** -0.2 * g.uL ** 2 * g.SL / g.AL; })();
      const r = taitelDuklerResidual(h, X2, Y);
      assert.ok(Math.abs(r) < 1e-8 * Math.max(L, Math.abs(4 * Y), 1), `X2=${X2} Y=${Y} residual ${r}`);
    }
  }
  // Level rises monotonically with X (more liquid) and falls with Y (downward slope drains the pipe).
  let prev = 0;
  for (const X2 of [1e-4, 1e-2, 1, 100, 1e4]) { const h = equilibriumLevel({ X2 }); assert.ok(h > prev); prev = h; }
  assert.ok(equilibriumLevel({ X2: 1, Y: 100 }) < equilibriumLevel({ X2: 1, Y: 0 }));
});

test('Taitel–Dukler dimensionless groups by hand (turbulent–turbulent, 25 mm air–water)', () => {
  const p = taitelDuklerParameters({ ...aw25, vsl: 0.1, vsg: 10, theta: -1 * DEG });
  const ReL = 1000 * 0.1 * 0.025 / 1e-3;
  const ReG = 1.2 * 10 * 0.025 / 1.8e-5;
  const gL = 4 * 0.046 * ReL ** -0.2 / 0.025 * 1000 * 0.1 ** 2 / 2;
  const gG = 4 * 0.046 * ReG ** -0.2 / 0.025 * 1.2 * 10 ** 2 / 2;
  close(p.X2, gL / gG, 1e-13);
  close(p.Y, 998.8 * G * Math.sin(1 * DEG) / gG, 1e-12);
  close(p.F, Math.sqrt(1.2 / 998.8) * 10 / Math.sqrt(0.025 * G * Math.cos(1 * DEG)), 1e-13);
  close(p.K, p.F * Math.sqrt(ReL), 1e-13);
  close(p.T, Math.sqrt(gL / (998.8 * G * Math.cos(1 * DEG))), 1e-13);
});

test('Taitel & Dukler (1976) map, air–water 25 mm horizontal: the four reference points', () => {
  const at = (vsg, vsl) => taitelDukler({ ...aw25, vsl, vsg });
  const strat = at(1, 0.01);
  assert.match(strat.regime, /^stratified/);
  assert.equal(strat.regime, 'stratified-smooth');
  assert.ok(Math.abs(strat.residual) < 1e-8);
  assert.equal(at(1, 1).regime, 'intermittent');
  assert.equal(at(30, 0.1).regime, 'annular');
  assert.equal(at(0.5, 5).regime, 'dispersed-bubble');
});

test('Taitel & Dukler (1976) map: qualitative boundaries and inclination effects', () => {
  const at = (vsg, vsl, deg = 0) => taitelDukler({ ...aw25, vsl, vsg, theta: deg * DEG }).regime;
  // stratified → intermittent with rising liquid rate at low gas rate (boundary near vsl ≈ 0.1–0.2 m/s)
  assert.match(at(1, 0.05), /^stratified/);
  assert.equal(at(1, 0.5), 'intermittent');
  // smooth → wavy with rising gas rate, then annular
  assert.equal(at(2, 0.01), 'stratified-smooth');
  assert.equal(at(10, 0.01), 'stratified-wavy');
  assert.equal(at(40, 0.01), 'annular');
  // intermittent → dispersed bubble with rising liquid rate
  assert.equal(at(0.5, 1), 'intermittent');
  assert.equal(at(0.5, 8), 'dispersed-bubble');
  // a slight upward slope removes low-rate stratified flow; a downward slope extends it (paper Figs. 6–7)
  assert.equal(at(1, 0.01, 1), 'intermittent');
  assert.match(at(1, 0.3, -1), /^stratified/);
  // Barnea's level criterion (0.35 instead of 0.5) keeps unstable flow with 0.35 < hL/D < 0.5 intermittent
  const p = { ...aw25, vsl: 0.2, vsg: 6 };
  const base = taitelDukler(p);
  assert.ok(base.hL > 0.35 && base.hL < 0.5);
  assert.equal(base.regime, 'annular');
  assert.equal(taitelDukler({ ...p, annularLevel: 0.35 }).regime, 'intermittent');
});

test('Taitel, Barnea & Dukler (1980) vertical map: hand-calculated transition velocities for air–water', () => {
  // annular: vsg = 3.1 (σ g Δρ)^¼ / √ρG = 3.1 × (0.072 × 9.80665 × 998.8)^¼ / √1.2 = 3.1 × 5.153 / 1.0954 = 14.58 m/s
  const a = taitelBarneaDuklerVertical({ ...aw25, D: 0.05, vsl: 0.1, vsg: 20 });
  assert.equal(a.regime, 'annular');
  close(20 / a.annularRatio, 14.58, 2e-3);
  assert.notEqual(taitelBarneaDuklerVertical({ ...aw25, D: 0.05, vsl: 0.1, vsg: 12 }).regime, 'annular');
  // dispersed bubble in 50 mm: vm = 4 × 0.05^0.429 × (7.2e-5)^0.089 / (1e-6)^0.072 × (9.80665 × 0.9988)^0.446 = 3.54 m/s
  const d = taitelBarneaDuklerVertical({ ...aw25, D: 0.05, vsl: 4, vsg: 0.5 });
  assert.equal(d.regime, 'dispersed-bubble');
  close(4.5 / d.dispersedRatio, 3.54, 3e-3);
  // bubble flow needs [ρL² g D²/(Δρ σ)]^¼ ≥ 4.36: not at 25 mm, yes at 100 mm
  assert.equal(taitelBarneaDuklerVertical({ ...aw25, vsl: 0.5, vsg: 0.05 }).bubblePossible, false);
  const b = taitelBarneaDuklerVertical({ ...aw25, D: 0.1, vsl: 0.5, vsg: 0.05 });
  assert.equal(b.regime, 'bubble');
  // bubble → slug at vsl = 3 vsg − 1.15 (g Δρ σ/ρL²)^¼ = 3 × 0.05 − 0.1874
  close(b.bubbleMargin, 0.5 - (0.15 - 1.15 * (G * 998.8 * 0.072 / 1e6) ** 0.25), 1e-12);
  assert.equal(taitelBarneaDuklerVertical({ ...aw25, D: 0.1, vsl: 0.1, vsg: 1 }).regime, 'slug');
  // churn is undeveloped slug flow: lE/D = 40.6 (vm/√(gD) + 0.22)
  const c = taitelBarneaDuklerVertical({ ...aw25, D: 0.05, vsl: 0.5, vsg: 1.5, LoverD: 50 });
  close(c.entranceLengthOverD, 40.6 * (2 / Math.sqrt(G * 0.05) + 0.22), 1e-13);
  assert.equal(c.regime, 'churn');
  assert.equal(taitelBarneaDuklerVertical({ ...aw25, D: 0.05, vsl: 0.5, vsg: 1.5, LoverD: 500 }).regime, 'slug');
});

test('Beggs & Brill (1973) pattern boundaries L1–L4 by hand at λ = 0.5', () => {
  // vm chosen so that λ = 0.5
  const r = beggsBrillPattern({ vsl: 1, vsg: 1, D: 0.1 });
  close(r.L1, 256.3, 1e-3);      // 316 × 0.5^0.302
  close(r.L2, 0.005120, 1e-3);   // 0.0009252 × 2^2.4684
  close(r.L3, 0.2735, 1e-3);     // 0.10 × 2^1.4516
  close(r.L4, 53.37, 1e-3);      // 0.5 × 2^6.738
  close(r.Fr, 4 / (G * 0.1), 1e-14);
  assert.equal(r.pattern, 'intermittent'); // λ ≥ 0.4 and L3 < Fr ≤ L4
  assert.equal(beggsBrillPattern({ vsl: 0.01, vsg: 0.01, D: 0.1 }).pattern, 'segregated'); // Fr = 4.1e-4 < L2
  assert.equal(beggsBrillPattern({ vsl: 0.03, vsg: 0.03, D: 0.1 }).pattern, 'segregated');
  assert.equal(beggsBrillPattern({ vsl: 0.1, vsg: 0.1, D: 0.1 }).pattern, 'transition');   // Fr = 0.041 in [L2, L3]
  assert.equal(beggsBrillPattern({ vsl: 5, vsg: 5, D: 0.1 }).pattern, 'distributed');      // Fr = 102 > L4
  assert.equal(beggsBrillPattern({ vsl: 0.005, vsg: 2, D: 0.1 }).pattern, 'segregated');   // λ < 0.01, Fr < L1
});

test('flowRegime: model selection by inclination with the documented blend', () => {
  const p = { ...aw25, D: 0.05, vsl: 0.3, vsg: 1 };
  const h = flowRegime({ ...p, theta: 2 * DEG });
  assert.equal(h.model, 'taitel-dukler'); assert.equal(h.interpolated, false);
  const v = flowRegime({ ...p, theta: 85 * DEG });
  assert.equal(v.model, 'taitel-barnea-dukler'); assert.equal(v.regime, 'slug');
  const m = flowRegime({ ...p, theta: 30 * DEG });
  assert.equal(m.interpolated, true); assert.equal(m.model, 'taitel-dukler'); assert.ok(m.vertical && m.horizontal);
  assert.equal(flowRegime({ ...p, theta: 50 * DEG }).model, 'taitel-barnea-dukler');
  assert.equal(flowRegime({ ...p, theta: -30 * DEG }).extrapolated, true);
});

test('regimeMap: grid classification is consistent with point evaluation and contains the classic regions', () => {
  const m = regimeMap({ ...aw25, model: 'taitel-dukler', nx: 24, ny: 20, point: { vsg: 1, vsl: 1 } });
  assert.equal(m.xs.length, 24); assert.equal(m.ys.length, 20); assert.equal(m.z.length, 20); assert.equal(m.z[0].length, 24);
  for (const name of ['stratified-smooth', 'stratified-wavy', 'intermittent', 'annular', 'dispersed-bubble']) assert.ok(m.legend.includes(name), name);
  assert.equal(m.point.regime, 'intermittent');
  for (const [iy, ix] of [[0, 0], [7, 13], [19, 23], [12, 3]]) {
    assert.equal(m.regimes[iy][ix], taitelDukler({ ...aw25, vsg: m.xs[ix], vsl: m.ys[iy] }).regime);
    assert.equal(m.legend[m.z[iy][ix]], m.regimes[iy][ix]);
  }
  const bb = regimeMap({ D: 0.1, model: 'beggs-brill', nx: 10, ny: 10 });
  assert.ok(bb.legend.includes('segregated') && bb.legend.includes('distributed'));
  assert.throws(() => regimeMap({ ...aw25, model: 'nope' }), RangeError);
});

test('stratified shear balance: static fluid leaves only the buoyancy term; smooth interface uses fG', () => {
  const p = { D: 0.1, rhoL: 850, rhoG: 20, muL: 2e-3, muG: 1.5e-5 };
  const s = stratifiedForces({ ...p, h: 0.4, vL: 0, vG: 0, theta: 0.1 });
  close(s.residual, -(850 - 20) * G * Math.sin(0.1), 1e-13);
  close(s.gradGravity, (s.holdup * 850 + (1 - s.holdup) * 20) * G * Math.sin(0.1), 1e-13);
  const f = stratifiedForces({ ...p, h: 0.4, vL: 0.5, vG: 5 });
  assert.equal(f.fI, f.fG);
  close(f.tauWL, f.fL * 850 * 0.25 / 2, 1e-13);
  close(f.tauI, f.fG * 20 * 4.5 ** 2 / 2, 1e-13);
  close(f.gradFriction, (f.tauWL * f.SL + f.tauWG * f.SG) / (Math.PI * 0.01 / 4), 1e-13);
  assert.equal(stratifiedForces({ ...p, h: 0.4, vL: 0.5, vG: 5, interfacial: 'constant' }).fI, 0.0142);
  // Andritsos–Hanratty: no enhancement below vsg,t = 5 √(1.2/ρG), enhancement above it
  const lo = stratifiedForces({ ...p, h: 0.3, vL: 0.5, vG: 1, interfacial: 'andritsos' });
  const hi = stratifiedForces({ ...p, h: 0.3, vL: 0.5, vG: 8, interfacial: 'andritsos' });
  assert.equal(lo.fI, lo.fG);
  const vsg = 8 * (1 - hi.holdup);
  close(hi.fI / hi.fG, 1 + 15 * Math.sqrt(0.3) * (vsg / (5 * Math.sqrt(1.2 / 20)) - 1), 1e-12);
});

test('invalid arguments throw RangeError', () => {
  assert.throws(() => taitelDukler({ ...aw25, vsl: 0, vsg: 1 }), /vsl/);
  assert.throws(() => taitelDukler({ ...aw25, vsl: 1, vsg: 1, rhoG: 2000 }), /rhoG/);
  assert.throws(() => stratifiedGeometry(1.2), /h/);
  assert.throws(() => beggsBrillPattern({ vsl: 0, vsg: 0, D: 0.1 }), RangeError);
});
