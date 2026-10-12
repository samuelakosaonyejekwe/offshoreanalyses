import assert from 'node:assert/strict';

export const close = (a, b, tol, msg = '') => assert.ok(Math.abs(a - b) <= tol, `${msg} expected ${b} ± ${tol}, got ${a}`);
export const rel = (a, b, tol, msg = '') => assert.ok(Math.abs(a - b) <= tol * Math.abs(b), `${msg} expected ${b} within ${tol * 100} %, got ${a}`);
export const closeV = (a, b, tol, msg = '') => { assert.equal(a.length, b.length); for (let i = 0; i < a.length; i++) close(a[i], b[i], tol, `${msg}[${i}]`); };

/** Regular polygon loop (counter-clockwise) as packed x, y. */
export const circle = (r, n, cx = 0, cy = 0) => Float64Array.from({ length: 2 * n }, (_, i) => (i % 2 ? cy + r * Math.sin(2 * Math.PI * ((i - 1) / 2) / n) : cx + r * Math.cos(2 * Math.PI * (i / 2) / n)));

/** Unit-cube hexahedral mesh with n cells per side. */
export function hexBlock(n, size = 1) {
  const m = n + 1;
  const positions = new Float64Array(3 * m * m * m);
  const id = (i, j, k) => (k * m + j) * m + i;
  for (let k = 0; k < m; k++) for (let j = 0; j < m; j++) for (let i = 0; i < m; i++) positions.set([size * i / n, size * j / n, size * k / n], 3 * id(i, j, k));
  const conn = [];
  for (let k = 0; k < n; k++) for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) conn.push(id(i, j, k), id(i + 1, j, k), id(i + 1, j + 1, k), id(i, j + 1, k), id(i, j, k + 1), id(i + 1, j, k + 1), id(i + 1, j + 1, k + 1), id(i, j + 1, k + 1));
  return { kind: 'volume', name: 'block', positions, cells: [{ type: 'hex8', conn: Uint32Array.from(conn) }] };
}
