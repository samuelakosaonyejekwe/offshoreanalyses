import { test } from 'node:test';
import assert from 'node:assert/strict';
import { form, sorm, importanceSampling } from '../../src/engine/reliability/form.js';
import { monteCarlo } from '../../src/engine/reliability/sampling.js';
import { normCdf, normPdf } from '../../src/engine/math/stats.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ''} expected ${b}, got ${a} (tol ${tol})`);

test('form: linear limit state with normal variables — beta = (muR - muS)/sqrt(sR² + sS²) exactly', () => {
  const vars = { R: { dist: 'normal', mean: 10, std: 1.5 }, S: { dist: 'normal', mean: 6, std: 1.2 } };
  const sd = Math.sqrt(1.5 ** 2 + 1.2 ** 2);
  const beta = 4 / sd;
  const r = form((x) => x.R - x.S, vars);
  assert.ok(r.converged);
  close(r.beta, beta, 1e-9);
  close(r.pf, normCdf(-beta), 1e-12);
  // direction cosines: resistance negative, load positive (u* = beta alpha)
  close(r.alpha.R, -1.5 / sd, 1e-7);
  close(r.alpha.S, 1.2 / sd, 1e-7);
  close(r.importance.R + r.importance.S, 1, 1e-12);
  // design point: R* = muR - beta sR²/sd = S*
  close(r.designPoint.R, 10 - (beta * 1.5 ** 2) / sd, 1e-7);
  close(r.designPoint.S, 6 + (beta * 1.2 ** 2) / sd, 1e-7);
  close(r.designPoint.R, r.designPoint.S, 1e-7);
  assert.ok(r.iterations <= 3);
});

test('form: three-variable linear case with a deterministic term and {mean, cov} inputs', () => {
  // g = R - D - L - c ; muR 300 (cov 0.1), D 100 (cov 0.05), L 80 (cov 0.25), c = 20
  const vars = {
    R: { dist: 'normal', mean: 300, cov: 0.1 },
    D: { dist: 'normal', mean: 100, cov: 0.05 },
    L: { dist: 'normal', mean: 80, cov: 0.25 },
    c: { dist: 'deterministic', value: 20 },
  };
  const beta = (300 - 100 - 80 - 20) / Math.sqrt(30 ** 2 + 5 ** 2 + 20 ** 2);
  const r = form((x) => x.R - x.D - x.L - x.c, vars);
  close(r.beta, beta, 1e-9);
  assert.equal(r.designPoint.c, 20);
  close(r.alpha.c, 0, 1e-12);
});

test('form: lognormal R and S, g = R - S — closed form beta = ln[(muR/muS) sqrt((1+VS²)/(1+VR²))] / sqrt(ln((1+VR²)(1+VS²)))', () => {
  const muR = 150;
  const VR = 0.12;
  const muS = 80;
  const VS = 0.25;
  const vars = { R: { dist: 'lognormal', mean: muR, cov: VR }, S: { dist: 'lognormal', mean: muS, cov: VS } };
  const beta = Math.log((muR / muS) * Math.sqrt((1 + VS * VS) / (1 + VR * VR))) / Math.sqrt(Math.log((1 + VR * VR) * (1 + VS * VS)));
  const r = form((x) => x.R - x.S, vars);
  assert.ok(r.converged);
  close(r.beta, beta, 1e-6);
  // the same surface written as ln R - ln S or R/S - 1 gives the same index (FORM invariance)
  close(form((x) => Math.log(x.R) - Math.log(x.S), vars).beta, beta, 1e-8);
  close(form((x) => x.R / x.S - 1, vars).beta, beta, 1e-6);
  close(r.designPoint.R, r.designPoint.S, 1e-4);
  // Monte Carlo agrees with the (here exact) FORM probability
  const mc = monteCarlo((x) => x.R - x.S, vars, { n: 400000, seed: 8 });
  const pf = normCdf(-beta);
  close(mc.pf, pf, 4 * Math.sqrt(pf / 400000), `pf exact ${pf}`);
  assert.ok(mc.ciPf[0] < r.pf && r.pf < mc.ciPf[1]);
});

test('form: correlated normal variables — beta = (muR - muS)/sqrt(sR² + sS² - 2 rho sR sS)', () => {
  const vars = { R: { dist: 'normal', mean: 10, std: 1.5 }, S: { dist: 'normal', mean: 6, std: 1.2 } };
  for (const rho of [0.5, -0.4]) {
    const beta = 4 / Math.sqrt(1.5 ** 2 + 1.2 ** 2 - 2 * rho * 1.5 * 1.2);
    close(form((x) => x.R - x.S, vars, { corr: { 'R,S': rho } }).beta, beta, 1e-8, `rho=${rho}`);
  }
  // correlated lognormals: ln R, ln S jointly normal with rho_ln = ln(1 + rho VR VS)/(zR zS)
  const lv = { R: { dist: 'lognormal', mean: 150, cov: 0.12 }, S: { dist: 'lognormal', mean: 80, cov: 0.25 } };
  const zR = Math.sqrt(Math.log(1 + 0.12 ** 2));
  const zS = Math.sqrt(Math.log(1 + 0.25 ** 2));
  const rhoLn = Math.log(1 + 0.3 * 0.12 * 0.25) / (zR * zS);
  const num = Math.log(150) - 0.5 * zR * zR - (Math.log(80) - 0.5 * zS * zS);
  const beta = num / Math.sqrt(zR * zR + zS * zS - 2 * rhoLn * zR * zS);
  close(form((x) => x.R - x.S, lv, { corr: [['R', 'S', 0.3]] }).beta, beta, 1e-6);
});

test('form: non-normal marginals (Weibull resistance, Gumbel load) — FORM, importance sampling and Monte Carlo agree', () => {
  const vars = {
    R: { dist: 'weibull', scale: 12, shape: 8 },
    S: { dist: 'gumbel', loc: 4, scale: 0.8 },
    k: { dist: 'uniform', min: 0.9, max: 1.1 },
  };
  const g = (x) => x.k * x.R - x.S;
  const r = form(g, vars);
  assert.ok(r.converged);
  // the design point lies on the limit-state surface
  close(g(r.designPoint), 0, 1e-4);
  const mc = monteCarlo(g, vars, { n: 400000, seed: 15 });
  assert.ok(mc.nFail > 300);
  // FORM is a first-order approximation here (negative curvature): within 25 % of the simulated probability
  assert.ok(Math.abs(r.pf / mc.pf - 1) < 0.25, `FORM ${r.pf} vs MC ${mc.pf}`);
  // the same design point is found from different starting points
  close(form(g, vars, { u0: [0, 3, 0] }).beta, r.beta, 1e-7);
  const is = importanceSampling(g, vars, { designPoint: r.designPoint, n: 20000, seed: 3 });
  assert.ok(is.cov < 0.05, `cov ${is.cov}`);
  close(is.pf, mc.pf, 4 * (is.cov * is.pf + mc.covPf * mc.pf), `IS ${is.pf} vs MC ${mc.pf}`);
  // SORM improves on FORM
  const s2 = sorm(g, vars);
  assert.ok(s2.valid);
  assert.ok(Math.abs(s2.pfBreitung / mc.pf - 1) < 0.08, `SORM ${s2.pfBreitung}, FORM ${r.pf}, MC ${mc.pf}`);
  assert.ok(Math.abs(s2.pfHohenbichler / mc.pf - 1) < 0.06, `SORM-HR ${s2.pfHohenbichler} vs MC ${mc.pf}`);
});

test('form: non-linear limit state g = x1³ + x2³ - 18 with N(10, 5) variables — symmetric design point x1 = x2 = 9^(1/3)', () => {
  const vars = { x1: { dist: 'normal', mean: 10, std: 5 }, x2: { dist: 'normal', mean: 10, std: 5 } };
  const g = (x) => x.x1 ** 3 + x.x2 ** 3 - 18;
  const r = form(g, vars);
  assert.ok(r.converged);
  // by symmetry the design point has x1 = x2 = 9^(1/3), so beta = sqrt(2) (10 - 9^(1/3))/5
  const beta = (Math.SQRT2 * (10 - Math.cbrt(9))) / 5;
  close(r.beta, beta, 1e-6);
  close(r.designPoint.x1, Math.cbrt(9), 1e-5);
  close(r.alpha.x1, -Math.SQRT1_2, 1e-5);
});

test('form: origin inside the failure domain gives a negative beta and pf > 0.5', () => {
  const vars = { R: { dist: 'normal', mean: 5, std: 1 }, S: { dist: 'normal', mean: 7, std: 1 } };
  const r = form((x) => x.R - x.S, vars);
  close(r.beta, -Math.SQRT2, 1e-9);
  close(r.pf, normCdf(Math.SQRT2), 1e-10);
});

test('form: a limit state that ignores the variables is rejected', () => {
  assert.throws(() => form(() => 1, { a: { dist: 'normal', mean: 0, std: 1 } }), /gradient/);
  assert.throws(() => form('nope', {}), RangeError);
});

test('sorm: linear limit state has zero curvature and pfBreitung = pfForm', () => {
  const vars = { a: { dist: 'normal', mean: 0, std: 1 }, b: { dist: 'normal', mean: 0, std: 1 }, c: { dist: 'normal', mean: 0, std: 1 } };
  const r = sorm((x) => 3 - (x.a + 2 * x.b - x.c) / Math.sqrt(6), vars);
  close(r.beta, 3, 1e-8);
  for (const k of r.curvatures) close(k, 0, 1e-6);
  close(r.pfBreitung, normCdf(-3), 1e-9);
  close(r.betaGeneralised, 3, 1e-6);
});

test('sorm: paraboloid g = beta + kappa (u1² + u2²)/2 - u3 — curvatures recovered, Breitung vs exact integral', () => {
  const sn = { dist: 'normal', mean: 0, std: 1 };
  const vars = { u1: sn, u2: sn, u3: sn };
  const beta = 3;
  for (const kappa of [0.15, -0.1]) {
    const r = sorm((x) => beta + 0.5 * kappa * (x.u1 ** 2 + x.u2 ** 2) - x.u3, vars);
    close(r.beta, beta, 1e-7);
    assert.equal(r.curvatures.length, 2);
    close(r.curvatures[0], kappa, 1e-5);
    close(r.curvatures[1], kappa, 1e-5);
    close(r.pfBreitung, normCdf(-beta) / (1 + beta * kappa), 1e-9);
    // exact: t = (u1² + u2²)/2 is Exp(1) distributed, pf = ∫ Phi(-(beta + kappa t)) e^-t dt (Simpson, t in [0, 40])
    const n = 40000;
    const h = 40 / n;
    let s = 0;
    for (let i = 0; i <= n; i++) {
      const t = i * h;
      const w = i === 0 || i === n ? 1 : i % 2 ? 4 : 2;
      s += w * normCdf(-(beta + kappa * t)) * Math.exp(-t);
    }
    const exact = (s * h) / 3;
    // Breitung is asymptotic: within 8 % at beta = 3; the Hohenbichler–Rackwitz form within 3 %
    assert.ok(Math.abs(r.pfBreitung / exact - 1) < 0.08, `kappa ${kappa}: Breitung ${r.pfBreitung} vs exact ${exact}`);
    assert.ok(Math.abs(r.pfHohenbichler / exact - 1) < 0.03, `kappa ${kappa}: HR ${r.pfHohenbichler} vs exact ${exact}`);
    // and both are much closer than FORM
    assert.ok(Math.abs(r.pfBreitung / exact - 1) < Math.abs(r.pfForm / exact - 1));
    close(normCdf(-r.betaGeneralised), r.pfBreitung, 1e-12);
  }
  assert.ok(normPdf(0) > 0);
});

test('importanceSampling: small probability of a linear limit state (beta = 4, pf = 3.17e-5) with 10^4 samples', () => {
  const vars = { R: { dist: 'normal', mean: 0, std: 1 }, S: { dist: 'normal', mean: 0, std: 1 } };
  const g = (x) => 4 - (x.R + x.S) / Math.SQRT2;
  const pf = normCdf(-4);
  const withDp = importanceSampling(g, vars, { designPoint: { R: 4 / Math.SQRT2, S: 4 / Math.SQRT2 }, n: 10000, seed: 5 });
  assert.ok(withDp.cov < 0.03);
  close(withDp.pf, pf, 4 * withDp.cov * pf);
  // design point found automatically by FORM, reproducible
  const auto = importanceSampling(g, vars, { n: 10000, seed: 5 });
  close(auto.pf, withDp.pf, 1e-9 * pf + 1e-12);
  assert.throws(() => importanceSampling(g, vars, { designPoint: { R: 1 } }), /designPoint\.S/);
});
