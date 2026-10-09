import test from 'node:test';
import assert from 'node:assert/strict';
import {
  bendForce, slugBendLoad, spoolForceHistory, dynamicLoadFactor, slugFatigueCycles, allowableSlugVelocity, momentumCheck,
} from '../../src/engine/flow/slugloads.js';

const close = (a, b, rel, msg) => assert.ok(Math.abs(a - b) <= rel * Math.max(Math.abs(b), 1e-300), `${msg ?? ''} got ${a}, expected ${b}`);
const D = 0.2;
const A = Math.PI * D * D / 4;

test('bend momentum force: ρ A v² √(2(1 − cos θ)) — √2 for 90°, 2 for 180°, 0 for a straight pipe', () => {
  const F = 800 * A * 25;
  const b90 = bendForce({ rho: 800, v: 5, D });
  close(b90.axialForce, F, 1e-14);
  close(b90.resultant, Math.SQRT2 * F, 1e-14);
  close(Math.hypot(b90.Fx, b90.Fy), b90.resultant, 1e-14);
  close(bendForce({ rho: 800, v: 5, D, bendAngle: Math.PI }).resultant, 2 * F, 1e-14);
  assert.equal(bendForce({ rho: 800, v: 5, D, bendAngle: 0 }).resultant, 0);
  close(bendForce({ rho: 800, v: 5, D, bendAngle: Math.PI / 4 }).resultant, F * Math.sqrt(2 - Math.SQRT2), 1e-13);
});

test('slug load at a bend: force range from the slug/film density difference', () => {
  const l = slugBendLoad({ rhoSlug: 700, rhoFilm: 100, v: 6, D });
  close(l.forceRange, Math.SQRT2 * 600 * A * 36, 1e-13);
  close(l.slugForce - l.filmForce, l.forceRange, 1e-13);
  assert.throws(() => slugBendLoad({ rhoSlug: 50, rhoFilm: 100, v: 6, D }), /rhoSlug/);
});

test('dynamic load factor, undamped rectangular pulse: 2 sin(π td/Tn) for td < Tn/2, else 2 (Biggs 1964)', () => {
  const fn = 5;
  for (const ratio of [0.1, 0.25, 0.4]) {
    close(dynamicLoadFactor({ naturalFrequency: fn, duration: ratio / fn }).dlf, 2 * Math.sin(Math.PI * ratio), 1e-5, `td/Tn=${ratio}`);
  }
  close(dynamicLoadFactor({ naturalFrequency: fn, duration: 2 / fn }).dlf, 2, 1e-6);
  close(dynamicLoadFactor({ naturalFrequency: fn, duration: 0.5 / fn }).dlf, 2, 1e-5);
});

test('dynamic load factor: slow ramp is quasi-static, ramp of one natural period gives 1, damping lowers the peak', () => {
  const fn = 5;
  close(dynamicLoadFactor({ naturalFrequency: fn, duration: 1, riseTime: 4 }).dlf, 1, 0.02);
  // undamped response to a ramp-then-hold load: DLF = 1 + |sin(π tr/Tn)|/(π tr/Tn); tr = Tn → exactly 1
  close(dynamicLoadFactor({ naturalFrequency: fn, duration: 3, riseTime: 0.2 }).dlf, 1, 1e-5);
  close(dynamicLoadFactor({ naturalFrequency: fn, duration: 3, riseTime: 0.1 }).dlf, 1 + 1 / (Math.PI * 0.5), 1e-5);
  // damped step: overshoot 1 + exp(−ζ π/√(1 − ζ²))
  const z = 0.05;
  close(dynamicLoadFactor({ naturalFrequency: fn, duration: 3, damping: z }).dlf, 1 + Math.exp(-z * Math.PI / Math.sqrt(1 - z * z)), 1e-5);
});

test('spool force history: pulse timing, plateau, rise time and overlap of two bends', () => {
  const v = 5;
  const h = spoolForceHistory({ D, rhoSlug: 700, rhoFilm: 100, v, bends: [{ s: 0, radius: 0.5 }, { s: 20, bendAngle: Math.PI / 4, radius: 0.5 }],
    slugs: [{ tFront: 1, length: 10 }, { tFront: 8, length: 0.2 }], duration: 14, dt: 0.001 });
  const load = slugBendLoad({ rhoSlug: 700, rhoFilm: 100, v, D });
  close(h.range[0], load.forceRange, 1e-13);
  close(h.riseTime[0], 0.5 * (Math.PI / 2) / v, 1e-13);
  const at = (b, t) => h.forces[b][Math.round(t / 0.001)];
  close(at(0, 0.5), load.filmForce, 1e-13);              // before the slug
  close(at(0, 2.0), load.slugForce, 1e-13);              // slug body in the bend (arrives 1 s, leaves 3 s + rise)
  close(at(0, 1 + h.riseTime[0] / 2), load.filmForce + 0.5 * load.forceRange, 1e-2);
  close(at(0, 3.5), load.filmForce, 1e-13);
  close(h.peak[0], load.slugForce, 1e-13);
  // second bend 20 m downstream: the front arrives 4 s later
  const l2 = slugBendLoad({ rhoSlug: 700, rhoFilm: 100, v, D, bendAngle: Math.PI / 4 });
  close(at(1, 4.9), l2.filmForce, 1e-13);
  close(at(1, 6.0), l2.slugForce, 1e-13);
  // a slug shorter than the bend arc (0.2 m < 0.785 m) never loads the bend fully
  let pk = 0;
  for (let i = Math.round(7.5 / 0.001); i < Math.round(9 / 0.001); i++) pk = Math.max(pk, h.forces[0][i]);
  assert.ok(pk < load.slugForce && pk > load.filmForce);
  close(pk, load.filmForce + (0.2 / (0.5 * Math.PI / 2)) * load.forceRange, 2e-2);
});

test('fatigue-cycle list: one cycle per slug, stress range = DLF × ΔF × transfer factor', () => {
  const c = slugFatigueCycles({ forceRange: 2000, stressPerForce: 5000, slugFrequency: 0.1, slugFlowFraction: 0.5, dlf: 1.6 });
  close(c.cyclesPerYear, 0.1 * 365.25 * 86400 * 0.5, 1e-14);
  close(c.cycles[0].stressRange, 1.6 * 2000 * 5000, 1e-14);
  const b = slugFatigueCycles({ forceRange: 2000, stressPerForce: 5000, slugFrequency: 0.1, bins: [{ probability: 0.7, dlf: 1.2 }, { probability: 0.3, dlf: 2 }] });
  close(b.cycles[0].cyclesPerYear + b.cycles[1].cyclesPerYear, b.cyclesPerYear, 1e-14);
  close(b.maxStressRange, 2 * 2000 * 5000, 1e-14);
  assert.throws(() => slugFatigueCycles({ forceRange: 1, stressPerForce: 1, slugFrequency: 1, bins: [{ probability: 0.5, dlf: 1 }] }), /sum/);
});

test('allowable slug velocity inverts the bend-load equation; ρv² check', () => {
  const v = allowableSlugVelocity({ allowableForce: 5000, rhoSlug: 700, rhoFilm: 100, D, dlf: 2 });
  close(2 * slugBendLoad({ rhoSlug: 700, rhoFilm: 100, v, D }).forceRange, 5000, 1e-12);
  const m = momentumCheck({ rho: 300, v: 4, limit: 6000 });
  close(m.momentum, 4800, 1e-14); close(m.utilisation, 0.8, 1e-14); assert.equal(m.ok, true);
  assert.equal(momentumCheck({ rho: 300, v: 4 }).ok, null);
  assert.equal(momentumCheck({ rho: 300, v: 5, limit: 6000 }).ok, false);
});
