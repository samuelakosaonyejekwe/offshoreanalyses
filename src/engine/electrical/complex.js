/**
 * Minimal complex arithmetic for phasor calculations (line models, sequence impedances, load flow).
 * A complex number is a plain object { re, im }. All functions are pure and return new objects.
 * Reference: any text on complex analysis; principal branches are used for arg, sqrt and log.
 */

/**
 * Build a complex number.
 * @param {number} re real part (any unit)
 * @param {number} [im=0] imaginary part (same unit)
 * @returns {{re:number, im:number}}
 */
export function cx(re, im = 0) {
  return { re, im };
}

/** Coerce a real number or complex object into a complex object. @param {number|{re:number,im:number}} z @returns {{re:number,im:number}} */
export function toCx(z) {
  return typeof z === 'number' ? { re: z, im: 0 } : z;
}

/** Sum a + b. @param {{re:number,im:number}} a @param {{re:number,im:number}} b @returns {{re:number,im:number}} */
export function cadd(a, b) {
  return { re: a.re + b.re, im: a.im + b.im };
}

/** Difference a − b. @param {{re:number,im:number}} a @param {{re:number,im:number}} b @returns {{re:number,im:number}} */
export function csub(a, b) {
  return { re: a.re - b.re, im: a.im - b.im };
}

/** Product a·b. @param {{re:number,im:number}} a @param {{re:number,im:number}} b @returns {{re:number,im:number}} */
export function cmul(a, b) {
  return { re: a.re * b.re - a.im * b.im, im: a.re * b.im + a.im * b.re };
}

/**
 * Quotient a / b (Smith's algorithm, avoids overflow for badly scaled operands).
 * @param {{re:number,im:number}} a @param {{re:number,im:number}} b @returns {{re:number,im:number}}
 */
export function cdiv(a, b) {
  if (b.re === 0 && b.im === 0) throw new RangeError('cdiv: divisor b is zero');
  if (Math.abs(b.re) >= Math.abs(b.im)) {
    const r = b.im / b.re;
    const d = b.re + b.im * r;
    return { re: (a.re + a.im * r) / d, im: (a.im - a.re * r) / d };
  }
  const r = b.re / b.im;
  const d = b.re * r + b.im;
  return { re: (a.re * r + a.im) / d, im: (a.im * r - a.re) / d };
}

/** Product of a complex number and a real scalar. @param {{re:number,im:number}} a @param {number} s @returns {{re:number,im:number}} */
export function cscale(a, s) {
  return { re: a.re * s, im: a.im * s };
}

/** Negation −a. @param {{re:number,im:number}} a @returns {{re:number,im:number}} */
export function cneg(a) {
  return { re: -a.re, im: -a.im };
}

/** Complex conjugate. @param {{re:number,im:number}} a @returns {{re:number,im:number}} */
export function cconj(a) {
  return { re: a.re, im: -a.im };
}

/** Reciprocal 1/a. @param {{re:number,im:number}} a @returns {{re:number,im:number}} */
export function cinv(a) {
  return cdiv({ re: 1, im: 0 }, a);
}

/** Modulus |a|. @param {{re:number,im:number}} a @returns {number} */
export function cabs(a) {
  return Math.hypot(a.re, a.im);
}

/** Argument of a, rad in (−π, π]. @param {{re:number,im:number}} a @returns {number} */
export function carg(a) {
  return Math.atan2(a.im, a.re);
}

/**
 * Complex number from polar form.
 * @param {number} mag modulus (any unit)
 * @param {number} angle argument, rad
 * @returns {{re:number,im:number}}
 */
export function cpolar(mag, angle) {
  return { re: mag * Math.cos(angle), im: mag * Math.sin(angle) };
}

/** Principal square root (real part ≥ 0). @param {{re:number,im:number}} a @returns {{re:number,im:number}} */
export function csqrt(a) {
  const m = Math.hypot(a.re, a.im);
  if (m === 0) return { re: 0, im: 0 };
  const re = Math.sqrt(0.5 * (m + Math.abs(a.re)));
  const im = a.im / (2 * re);
  if (a.re >= 0) return { re, im };
  return { re: Math.abs(im), im: a.im >= 0 ? re : -re };
}

/** Exponential e^a. @param {{re:number,im:number}} a @returns {{re:number,im:number}} */
export function cexp(a) {
  const e = Math.exp(a.re);
  return { re: e * Math.cos(a.im), im: e * Math.sin(a.im) };
}

/** Principal natural logarithm. @param {{re:number,im:number}} a @returns {{re:number,im:number}} */
export function clog(a) {
  return { re: Math.log(Math.hypot(a.re, a.im)), im: Math.atan2(a.im, a.re) };
}

/** Hyperbolic sine. @param {{re:number,im:number}} a @returns {{re:number,im:number}} */
export function csinh(a) {
  return { re: Math.sinh(a.re) * Math.cos(a.im), im: Math.cosh(a.re) * Math.sin(a.im) };
}

/** Hyperbolic cosine. @param {{re:number,im:number}} a @returns {{re:number,im:number}} */
export function ccosh(a) {
  return { re: Math.cosh(a.re) * Math.cos(a.im), im: Math.sinh(a.re) * Math.sin(a.im) };
}

/** Hyperbolic tangent. @param {{re:number,im:number}} a @returns {{re:number,im:number}} */
export function ctanh(a) {
  return cdiv(csinh(a), ccosh(a));
}

/**
 * Parallel combination a·b/(a + b) of two impedances.
 * @param {{re:number,im:number}} a @param {{re:number,im:number}} b @returns {{re:number,im:number}}
 */
export function cparallel(a, b) {
  return cdiv(cmul(a, b), cadd(a, b));
}

/**
 * Solve the complex linear system A x = b by Gaussian elimination with partial pivoting.
 * Intended for small dense systems (n up to a few tens: mesh equations, sequence networks).
 * @param {{re:number,im:number}[][]} A n×n matrix (not modified)
 * @param {{re:number,im:number}[]} b right-hand side, length n (not modified)
 * @returns {{re:number,im:number}[]} solution vector
 */
export function csolve(A, b) {
  const n = A.length;
  if (!Array.isArray(b) || b.length !== n) throw new RangeError('csolve: b must have the same length as A');
  const M = A.map((row, i) => {
    if (row.length !== n) throw new RangeError('csolve: A must be square');
    return [...row.map(toCx), toCx(b[i])];
  });
  for (let k = 0; k < n; k++) {
    let p = k;
    let best = cabs(M[k][k]);
    for (let i = k + 1; i < n; i++) {
      const v = cabs(M[i][k]);
      if (v > best) { best = v; p = i; }
    }
    if (best === 0) throw new RangeError('csolve: A is singular');
    if (p !== k) { const t = M[p]; M[p] = M[k]; M[k] = t; }
    for (let i = k + 1; i < n; i++) {
      const f = cdiv(M[i][k], M[k][k]);
      if (f.re === 0 && f.im === 0) continue;
      for (let j = k; j <= n; j++) M[i][j] = csub(M[i][j], cmul(f, M[k][j]));
    }
  }
  const x = new Array(n);
  for (let i = n - 1; i >= 0; i--) {
    let s = M[i][n];
    for (let j = i + 1; j < n; j++) s = csub(s, cmul(M[i][j], x[j]));
    x[i] = cdiv(s, M[i][i]);
  }
  return x;
}

/**
 * Product of two 2×2 complex transmission (ABCD) matrices, M = M1·M2.
 * @param {{A,B,C,D}} m1 first two-port (nearest the sending end)
 * @param {{A,B,C,D}} m2 second two-port
 * @returns {{A:{re,im},B:{re,im},C:{re,im},D:{re,im}}}
 */
export function abcdMul(m1, m2) {
  return {
    A: cadd(cmul(m1.A, m2.A), cmul(m1.B, m2.C)),
    B: cadd(cmul(m1.A, m2.B), cmul(m1.B, m2.D)),
    C: cadd(cmul(m1.C, m2.A), cmul(m1.D, m2.C)),
    D: cadd(cmul(m1.C, m2.B), cmul(m1.D, m2.D)),
  };
}
