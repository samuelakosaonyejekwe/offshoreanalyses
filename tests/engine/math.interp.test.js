import test from 'node:test';
import assert from 'node:assert/strict';
import {
  interp1, makeLinear, makePchip, makeSpline, interp2, linspace, logspace, resample, decimate,
} from '../../src/engine/math/interp.js';

function close(a, b, tol, msg = '') {
  assert.ok(Math.abs(a - b) <= tol, `${msg} got ${a}, expected ${b} (tol ${tol})`);
}

test('interp1: interior, nodes, and the three extrapolation modes', () => {
  const xs = [0, 1, 3];
  const ys = [0, 10, 30];
  assert.equal(interp1(xs, ys, 0.5), 5);
  assert.equal(interp1(xs, ys, 2), 20);
  assert.equal(interp1(xs, ys, 1), 10);
  assert.equal(interp1(xs, ys, 3), 30);
  assert.equal(interp1(xs, ys, -1), 0, 'clamp is the default');
  assert.equal(interp1(xs, ys, 9), 30);
  assert.equal(interp1(xs, ys, 4, { extrapolate: 'linear' }), 40);
  assert.equal(interp1(xs, ys, -1, { extrapolate: 'linear' }), -10);
  assert.throws(() => interp1(xs, ys, 4, { extrapolate: 'error' }), RangeError);
  assert.equal(interp1(new Float64Array(xs), new Float64Array(ys), 0.25), 2.5);
  assert.ok(Number.isNaN(interp1(xs, ys, NaN)));
});

test('interp1: rejects unsorted tables, mismatched lengths and unknown modes', () => {
  assert.throws(() => interp1([0, 2, 1], [0, 1, 2], 0.5), RangeError);
  assert.throws(() => interp1([0, 1, 1], [0, 1, 2], 0.5), RangeError);
  assert.throws(() => interp1([0, 1], [0, 1, 2], 0.5), RangeError);
  assert.throws(() => interp1([0, 1], [0, 1], 0.5, { extrapolate: 'nearest' }), RangeError);
});

test('makeLinear: matches interp1 on a long table (binary search) and copies its inputs', () => {
  const xs = [];
  const ys = [];
  for (let i = 0; i <= 1000; i++) {
    xs.push(i * i * 1e-3);
    ys.push(Math.sin(i * 0.01));
  }
  const f = makeLinear(xs, ys);
  for (const x of [0, 0.00037, 12.5, 499.9, 1000, 1e9, -5]) assert.equal(f(x), interp1(xs, ys, x));
  ys[0] = 99;
  assert.notEqual(f(0), 99);
  const g = makeLinear([0, 1], [0, 2], { extrapolate: 'linear' });
  assert.equal(g(2), 4);
});

test('makeSpline (not-a-knot): reproduces a cubic exactly, including first and second derivatives', () => {
  const p = (x) => 2 * x * x * x - 3 * x * x + 0.5 * x - 7;
  const d1 = (x) => 6 * x * x - 6 * x + 0.5;
  const d2 = (x) => 12 * x - 6;
  const xs = [-2, -1.2, -0.1, 0.4, 1.5, 1.9, 3];
  const s = makeSpline(xs, xs.map(p), { bc: 'not-a-knot' });
  for (const x of [-2, -1.7, -0.5, 0, 0.123, 1, 1.7, 2.5, 3]) {
    close(s(x), p(x), 1e-12, `value at ${x}`);
    close(s.d1(x), d1(x), 1e-11, `d1 at ${x}`);
    close(s.d2(x), d2(x), 1e-10, `d2 at ${x}`);
  }
  // four points: a single cubic through all of them
  const s4 = makeSpline([0, 1, 2, 4], [0, 1, 2, 4].map(p), { bc: 'not-a-knot' });
  close(s4(3), p(3), 1e-12);
  // three points: the interpolating parabola
  const s3 = makeSpline([0, 1, 3], [1, 2, 10], { bc: 'not-a-knot' });
  close(s3(2), 5, 1e-12, 'parabola x² + 1 through (0,1), (1,2), (3,10)');
  close(s3.d2(0.3), 2, 1e-12);
});

test('makeSpline (clamped): reproduces a cubic exactly when the exact end slopes are supplied', () => {
  const p = (x) => x * x * x - x;
  const xs = [0, 0.7, 1.1, 2];
  const s = makeSpline(xs, xs.map(p), { bc: 'clamped', d0: -1, dn: 11 });
  for (const x of [0, 0.2, 0.9, 1.5, 2]) {
    close(s(x), p(x), 1e-13);
    close(s.d1(x), 3 * x * x - 1, 1e-12);
  }
  close(s.d1(0), -1, 1e-13);
  close(s.d1(2), 11, 1e-12);
  // two points with clamped ends: the cubic Hermite segment
  const h = makeSpline([0, 1], [0, 0], { bc: 'clamped', d0: 1, dn: 1 });
  close(h(0.5), 0, 1e-15);
  close(h.d1(0), 1, 1e-14);
  close(h.d1(1), 1, 1e-14);
  close(h.d1(0.5), -0.5, 1e-14);
});

test('makeSpline (natural): interpolates, zero end curvature, C² continuity at the knots', () => {
  const xs = [0, 1, 2.5, 3, 4.2, 6];
  const ys = [1, 3, 2, 5, 4, 0];
  const s = makeSpline(xs, ys);
  for (let i = 0; i < xs.length; i++) close(s(xs[i]), ys[i], 1e-13);
  close(s.d2(xs[0]), 0, 1e-12);
  close(s.d2(xs[5]), 0, 1e-12);
  for (let i = 1; i < 5; i++) {
    const e = 1e-7;
    close(s.d1(xs[i] - e), s.d1(xs[i] + e), 1e-5, 'd1 continuous');
    close(s.d2(xs[i] - e), s.d2(xs[i] + e), 1e-5, 'd2 continuous');
  }
  // a natural spline through collinear points is that line
  const l = makeSpline([0, 1, 3, 7], [1, 3, 7, 15]);
  close(l(5), 11, 1e-13);
  close(l.d1(2), 2, 1e-13);
});

test('makeSpline: fourth-order convergence for sin x with not-a-knot ends (error ratio 16 per halving)', () => {
  const errs = [32, 64, 128, 256].map((n) => {
    const xs = linspace(0, Math.PI, n + 1);
    const s = makeSpline(xs, Array.from(xs, Math.sin), { bc: 'not-a-knot' });
    let e = 0;
    for (let i = 0; i <= 2000; i++) {
      const x = (Math.PI * i) / 2000;
      e = Math.max(e, Math.abs(s(x) - Math.sin(x)));
    }
    return e;
  });
  for (let k = 1; k < errs.length; k++) {
    const order = Math.log2(errs[k - 1] / errs[k]);
    assert.ok(order > 3.6 && order < 4.6, `observed order ${order}`);
  }
  assert.ok(errs[3] < 1e-9);
});

test('makeSpline: extrapolation modes and argument checks', () => {
  const xs = [0, 1, 2, 3];
  const ys = [0, 1, 8, 27];
  const lin = makeSpline(xs, ys, { bc: 'not-a-knot' });
  // default: straight continuation with the end slope (27 at x = 3 for x³)
  close(lin(4), 27 + 27, 1e-11);
  close(lin.d1(4), 27, 1e-11);
  assert.equal(lin.d2(4), 0);
  const ext = makeSpline(xs, ys, { bc: 'not-a-knot', extrapolate: 'extend' });
  close(ext(4), 64, 1e-11);
  const cl = makeSpline(xs, ys, { extrapolate: 'clamp' });
  assert.equal(cl(-5), 0);
  assert.equal(cl(10), 27);
  assert.throws(() => makeSpline(xs, ys, { extrapolate: 'error' })(3.5), RangeError);
  assert.throws(() => makeSpline([0], [1]), RangeError);
  assert.throws(() => makeSpline(xs, ys, { bc: 'periodic' }), RangeError);
  // two points: straight line whatever the end condition
  close(makeSpline([0, 2], [1, 5])(1), 3, 1e-15);
});

test('makePchip: interpolates the nodes and reproduces straight lines exactly', () => {
  const xs = [0, 1, 2.5, 4];
  const ys = [1, 3, 6, 9];
  const f = makePchip(xs, ys);
  for (let i = 0; i < 4; i++) close(f(xs[i]), ys[i], 1e-14);
  for (const x of [0.3, 1.7, 3.3]) close(f(x), 1 + 2 * x, 1e-13);
  close(f.d1(2), 2, 1e-13);
  close(makePchip([0, 2], [0, 4])(1.5), 3, 1e-15);
  assert.equal(makePchip([1], [7])(3), 7);
});

test('makePchip: monotone step-like data gives no overshoot, where the cubic spline overshoots', () => {
  const xs = [0, 1, 2, 3, 4, 5, 6];
  const ys = [0, 0, 0, 1, 1, 1, 1];
  const p = makePchip(xs, ys);
  const s = makeSpline(xs, ys);
  let prev = -Infinity;
  let splineMax = -Infinity;
  let splineMin = Infinity;
  for (let i = 0; i <= 600; i++) {
    const x = i / 100;
    const v = p(x);
    assert.ok(v >= -1e-15 && v <= 1 + 1e-15, `pchip overshoot ${v} at ${x}`);
    assert.ok(v >= prev - 1e-15, 'pchip is monotone');
    prev = v;
    splineMax = Math.max(splineMax, s(x));
    splineMin = Math.min(splineMin, s(x));
  }
  assert.ok(splineMax > 1.02 && splineMin < -0.02, 'the C² spline oscillates on the same data');
  // flat segments stay exactly flat
  assert.equal(p(0.5), 0);
  assert.equal(p(4.5), 1);
});

test('makePchip: local extremum in the data is preserved (zero slope at the peak) and d1 matches a finite difference', () => {
  const xs = [0, 1, 2, 3, 4];
  const ys = [0, 2, 5, 2, 0];
  const p = makePchip(xs, ys);
  assert.equal(p.d1(2), 0);
  for (let i = 0; i <= 400; i++) assert.ok(p(i / 100) <= 5 + 1e-14);
  const h = 1e-6;
  close(p.d1(1.3), (p(1.3 + h) - p(1.3 - h)) / (2 * h), 1e-7);
  // extrapolation: clamp by default, linear / extend on request
  assert.equal(p(9), 0);
  const pl = makePchip(xs, ys, { extrapolate: 'linear' });
  close(pl(5), 0 + pl.d1(4) * 1, 1e-14);
  assert.throws(() => makePchip(xs, ys, { extrapolate: 'error' })(-1), RangeError);
});

test('interp2: bilinear surface a + bx + cy + dxy is reproduced exactly; clamped outside the grid', () => {
  const f = (x, y) => 2 + 3 * x - y + 0.5 * x * y;
  const xs = [0, 1, 4];
  const ys = [-1, 2];
  const z = ys.map((y) => xs.map((x) => f(x, y)));
  for (const [x, y] of [[0, -1], [4, 2], [0.3, 0.5], [2.5, 1.9], [1, -0.2]]) close(interp2(xs, ys, z, x, y), f(x, y), 1e-13);
  assert.equal(interp2(xs, ys, z, -10, -10), f(0, -1));
  assert.equal(interp2(xs, ys, z, 10, 10), f(4, 2));
  close(interp2(xs, ys, z, 10, 0.5), f(4, 0.5), 1e-13);
  // z is indexed z[iy][ix]
  assert.equal(interp2([0, 1], [0, 1], [[0, 1], [10, 11]], 1, 0), 1);
  assert.equal(interp2([0, 1], [0, 1], [[0, 1], [10, 11]], 0, 1), 10);
  assert.equal(interp2([0, 1], [5], [[3, 5]], 0.5, 99), 4, 'single row');
  assert.throws(() => interp2(xs, ys, [[1, 2, 3]], 0, 0), RangeError);
});

test('linspace and logspace: end points exact, spacing uniform', () => {
  const a = linspace(0, 1, 11);
  assert.ok(a instanceof Float64Array);
  assert.equal(a.length, 11);
  assert.equal(a[0], 0);
  assert.equal(a[10], 1);
  close(a[3], 0.3, 1e-16);
  assert.deepEqual(Array.from(linspace(5, 5, 1)), [5]);
  assert.deepEqual(Array.from(linspace(2, -2, 3)), [2, 0, -2]);
  const l = logspace(-2, 2, 5);
  const expected = [0.01, 0.1, 1, 10, 100];
  for (let i = 0; i < 5; i++) close(l[i] / expected[i], 1, 1e-15);
  assert.throws(() => linspace(0, 1, 0), RangeError);
});

test('resample: linear data is preserved, end points kept, uniform output grid', () => {
  const xs = [0, 0.1, 0.5, 2, 2.1, 10];
  const ys = xs.map((x) => 3 * x - 2);
  const { x, y } = resample(xs, ys, 41);
  assert.equal(x.length, 41);
  assert.equal(x[0], 0);
  assert.equal(x[40], 10);
  for (let i = 0; i < 41; i++) close(y[i], 3 * x[i] - 2, 1e-12);
  // piecewise-linear data: values agree with interp1
  const zs = [1, 4, -2, 0, 7, 3];
  const r = resample(xs, zs, 101);
  for (let i = 0; i < 101; i++) close(r.y[i], interp1(xs, zs, r.x[i]), 1e-12);
  assert.throws(() => resample(xs, zs, 1), RangeError);
});

test('decimate: keeps end points, global extremes and order; output no longer than maxPoints', () => {
  const n = 100000;
  const xs = new Float64Array(n);
  const ys = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    xs[i] = i * 0.01;
    ys[i] = Math.sin(i * 0.003) + 0.2 * Math.sin(i * 0.71);
  }
  ys[31415] = 9; // isolated spike
  ys[77777] = -8;
  const { x, y } = decimate(xs, ys, 2000);
  assert.ok(x.length <= 2000 && x.length > 1900);
  assert.equal(x[0], xs[0]);
  assert.equal(x[x.length - 1], xs[n - 1]);
  assert.equal(y[0], ys[0]);
  let mx = -Infinity;
  let mn = Infinity;
  for (let i = 0; i < y.length; i++) {
    mx = Math.max(mx, y[i]);
    mn = Math.min(mn, y[i]);
    if (i > 0) assert.ok(x[i] > x[i - 1], 'abscissae stay ascending');
  }
  assert.equal(mx, 9);
  assert.equal(mn, -8);
});

test('decimate: short series returned unchanged as copies; argument checks', () => {
  const r = decimate([0, 1, 2], [5, 6, 7], 10);
  assert.deepEqual(Array.from(r.x), [0, 1, 2]);
  assert.deepEqual(Array.from(r.y), [5, 6, 7]);
  const d = decimate([0, 1, 2, 3, 4, 5, 6, 7], [0, 5, -5, 1, 2, 9, -9, 0], 4);
  assert.ok(d.x.length <= 4);
  assert.deepEqual(Array.from(d.y), [0, 9, -9, 0]);
  assert.equal(Math.max(...d.y), 9);
  assert.equal(Math.min(...d.y), -9);
  assert.throws(() => decimate([0, 1], [0], 10), RangeError);
  assert.throws(() => decimate([0, 1], [0, 1], 2), RangeError);
});
