// Boundary-element (panel-method) added mass of a hull in the two frequency limits of the linear free-surface
// problem, and of a body in unbounded fluid.
//
// Method: flat triangular panels on the wetted hull, collocation at the panel centroids. The free surface is
// replaced by an image of the hull about z = 0:
//
//   limit 'infinite' (omega -> infinity):  free-surface condition  phi = 0        -> image of OPPOSITE sign
//   limit 'zero'     (omega -> 0):         free-surface condition  d phi/dz = 0   -> image of EQUAL sign (rigid lid)
//   limit 'unbounded':                      no free surface, no image (deeply submerged body)
//
// With the Green function G = 1/r + s/r' (s = -1, +1, 0), the body condition d phi_k / dn = n_k for radiation
// mode k (n: unit normal pointing into the fluid, n_4..6 = (r - cog) × n) and A_jk = -rho ∫ phi_k n_j dS, two
// classical discretisations are offered:
//
//   formulation 'potential' (default): Green's third identity on the hull with piecewise-constant potential
//     and the known constant source strength n_k on every panel (Morino 1974),
//       2 pi phi_i - sum_{j != i} phi_j ∫_j dG/dn dS = - sum_j n_kj ∫_j G dS.
//     Observed convergence on a sphere: second order in the panel size (error ∝ 1/panels).
//   formulation 'source': constant-strength source panels (Hess & Smith 1964),
//       phi_k(x) = sum_j sigma_kj ∫_j G dS,  -2 pi sigma_i + sum_{j != i} sigma_j n_i · grad ∫_j G dS = n_ki.
//     Observed convergence: first order in the panel size (error ∝ 1/sqrt(panels)), several times less
//     accurate than 'potential' on the same mesh; kept as an independent cross-check.
//
// Panel integrals of 1/r, of its gradient and of the normal dipole (the solid angle of the triangle) are
// evaluated in closed form for every panel pair (edge logarithms, van Oosterom–Strackee solid angle),
// including the self-influence (jump term ±2 pi, vanishing principal value on a flat panel).
//
// References: Hess & Smith, "Calculation of nonlifting potential flow about arbitrary three-dimensional bodies",
// J. Ship Research 8 (1964); Morino & Kuo, "Subsonic potential aerodynamics for complex configurations: a
// general theory", AIAA J. 12 (1974); Newman, "Distributions of sources and normal dipoles over a quadrilateral panel",
// J. Eng. Math. 20 (1986); van Oosterom & Strackee, IEEE Trans. Biomed. Eng. 30 (1983) (solid angle of a
// triangle); Newman, Marine Hydrodynamics (1977) §4.13 (added mass), §6.18 (frequency limits).
//
// NOT a radiation/diffraction solver: there is no frequency dependence, no wave damping and no wave excitation.
// For A(omega), B(omega) and exciting forces export the hull to NEMOH, Capytaine or WAMIT (src/interop/nemoh.js)
// and read the results back into the container of hydrodb.js.

const MAX_PANELS = 6000;

/**
 * Wetted part of a triangulated hull: every triangle is clipped at the still-water plane z = 0, the part with
 * z < 0 is kept, and faces lying in the plane z = 0 (a lid, or the deck of a body floating flush) are dropped.
 * Normals are oriented into the fluid: the sign is fixed so that the displaced volume ∫ z n_z dS is positive,
 * which works for closed surfaces and for hulls that are already open at the waterline.
 * Reference: Sutherland & Hodgman, "Reentrant polygon clipping", Comm. ACM 17 (1974).
 * Validity: consistently oriented triangle surface, z up, waterline at z = 0; at least one triangle below z = 0.
 * @param {{positions:ArrayLike<number>, indices:ArrayLike<number>, origin?:number[]}} surface surface model
 *   (CONTRACT §4), coordinates in m; `origin` (m) is added to the positions when present
 * @param {{tol?:number}} [opts] tol: snapping distance to the waterplane, m (default 1e-6 × bounding-box diagonal)
 * @returns {{n:number, vertices:Float64Array, centroids:Float64Array, normals:Float64Array, areas:Float64Array,
 *   area:number, volume:number, clipped:number, flipped:boolean}} n panels; vertices 9 per panel (m), centroids
 *   and unit normals (into the fluid) 3 per panel, areas m², total wetted area m², displaced volume m³,
 *   number of source triangles cut by the waterplane, and whether the input orientation was reversed
 */
export function wettedPanels(surface, { tol } = {}) {
  if (!surface || !surface.positions || !surface.indices) throw new RangeError('surface must be a triangle surface with positions and indices');
  const P = surface.positions; const I = surface.indices;
  if (I.length < 3 || I.length % 3 !== 0) throw new RangeError('surface.indices must hold three vertex indices per triangle');
  const o = surface.origin || [0, 0, 0];
  let lo = [Infinity, Infinity, Infinity]; let hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < P.length; i += 3) for (let k = 0; k < 3; k++) { const v = P[i + k] + o[k]; if (v < lo[k]) lo[k] = v; if (v > hi[k]) hi[k] = v; }
  const diag = Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
  if (!(diag > 0) || !Number.isFinite(diag)) throw new RangeError('surface.positions must be finite and span a volume');
  const eps = tol === undefined ? 1e-6 * diag : tol;
  if (!(eps >= 0)) throw new RangeError('tol must be >= 0');
  const out = [];
  let clipped = 0;
  const poly = new Float64Array(12); const tri = new Float64Array(9);
  const push = (ax, ay, az, bx, by, bz, cx, cy, cz) => {
    const ux = bx - ax; const uy = by - ay; const uz = bz - az; const vx = cx - ax; const vy = cy - ay; const vz = cz - az;
    const a2 = Math.hypot(uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx);
    if (a2 > 1e-14 * diag * diag) out.push(ax, ay, az, bx, by, bz, cx, cy, cz);
  };
  for (let f = 0; f < I.length; f += 3) {
    let below = 0; let above = 0;
    for (let k = 0; k < 3; k++) {
      const q = 3 * I[f + k];
      tri[3 * k] = P[q] + o[0]; tri[3 * k + 1] = P[q + 1] + o[1];
      let z = P[q + 2] + o[2];
      if (Math.abs(z) <= eps) z = 0;
      tri[3 * k + 2] = z;
      if (z < 0) below++; else if (z > 0) above++;
    }
    if (below === 0) continue; // dry, or lying in the waterplane
    if (above === 0) { push(tri[0], tri[1], tri[2], tri[3], tri[4], tri[5], tri[6], tri[7], tri[8]); continue; }
    clipped++;
    let m = 0;
    for (let k = 0; k < 3; k++) {
      const a = 3 * k; const b = 3 * ((k + 1) % 3);
      const za = tri[a + 2]; const zb = tri[b + 2];
      if (za <= 0) { poly[3 * m] = tri[a]; poly[3 * m + 1] = tri[a + 1]; poly[3 * m + 2] = za; m++; }
      if ((za < 0 && zb > 0) || (za > 0 && zb < 0)) {
        const s = za / (za - zb);
        poly[3 * m] = tri[a] + s * (tri[b] - tri[a]); poly[3 * m + 1] = tri[a + 1] + s * (tri[b + 1] - tri[a + 1]); poly[3 * m + 2] = 0; m++;
      }
    }
    for (let k = 1; k + 1 < m; k++) push(poly[0], poly[1], poly[2], poly[3 * k], poly[3 * k + 1], poly[3 * k + 2], poly[3 * k + 3], poly[3 * k + 4], poly[3 * k + 5]);
  }
  const n = out.length / 9;
  if (n === 0) throw new RangeError('surface has no wetted panels (nothing lies below z = 0)');
  const vertices = Float64Array.from(out);
  const centroids = new Float64Array(3 * n); const normals = new Float64Array(3 * n); const areas = new Float64Array(n);
  let volume = 0; let area = 0;
  for (let p = 0; p < n; p++) {
    const q = 9 * p;
    const ux = vertices[q + 3] - vertices[q]; const uy = vertices[q + 4] - vertices[q + 1]; const uz = vertices[q + 5] - vertices[q + 2];
    const vx = vertices[q + 6] - vertices[q]; const vy = vertices[q + 7] - vertices[q + 1]; const vz = vertices[q + 8] - vertices[q + 2];
    const nx = uy * vz - uz * vy; const ny = uz * vx - ux * vz; const nz = ux * vy - uy * vx;
    const l = Math.hypot(nx, ny, nz);
    normals[3 * p] = nx / l; normals[3 * p + 1] = ny / l; normals[3 * p + 2] = nz / l;
    areas[p] = l / 2; area += l / 2;
    for (let k = 0; k < 3; k++) centroids[3 * p + k] = (vertices[q + k] + vertices[q + 3 + k] + vertices[q + 6 + k]) / 3;
    volume += centroids[3 * p + 2] * (nz / 2); // ∫ z n_z dS, exact on a flat triangle
  }
  let flipped = false;
  if (volume < 0) {
    flipped = true; volume = -volume;
    for (let p = 0; p < n; p++) {
      const q = 9 * p;
      for (let k = 0; k < 3; k++) { const t = vertices[q + 3 + k]; vertices[q + 3 + k] = vertices[q + 6 + k]; vertices[q + 6 + k] = t; normals[3 * p + k] = -normals[3 * p + k]; }
    }
  }
  return { n, vertices, centroids, normals, areas, area, volume, clipped, flipped };
}

/** Per-panel edge data used by the closed-form integrals. */
function panelGeometry(pan) {
  const { n, vertices: V, normals: N } = pan;
  const T = new Float64Array(9 * n); const M = new Float64Array(9 * n); const L = new Float64Array(3 * n); const diam = new Float64Array(n);
  for (let p = 0; p < n; p++) {
    for (let e = 0; e < 3; e++) {
      const a = 9 * p + 3 * e; const b = 9 * p + 3 * ((e + 1) % 3);
      const tx = V[b] - V[a]; const ty = V[b + 1] - V[a + 1]; const tz = V[b + 2] - V[a + 2];
      const l = Math.hypot(tx, ty, tz);
      const q = 9 * p + 3 * e;
      T[q] = tx / l; T[q + 1] = ty / l; T[q + 2] = tz / l;
      // in-plane outward normal of the edge, m = t × n
      M[q] = T[q + 1] * N[3 * p + 2] - T[q + 2] * N[3 * p + 1];
      M[q + 1] = T[q + 2] * N[3 * p] - T[q] * N[3 * p + 2];
      M[q + 2] = T[q] * N[3 * p + 1] - T[q + 1] * N[3 * p];
      L[3 * p + e] = l;
      if (l > diam[p]) diam[p] = l;
    }
  }
  return { T, M, L, diam };
}

/**
 * ∫ 1/r dS over panel j and its gradient with respect to the field point, closed form.
 * out[0] = potential, out[1..3] = gradient, out[4] = signed solid angle (negative when the field point is on the
 * side the normal points to; minus the normal-dipole integral). `self` takes the principal value on the panel.
 */
function panelIntegral(pan, geo, j, px, py, pz, self, out) {
  const V = pan.vertices; const N = pan.normals; const M = geo.M; const L = geo.L;
  const q = 9 * j;
  const ax = V[q] - px; const ay = V[q + 1] - py; const az = V[q + 2] - pz;
  const bx = V[q + 3] - px; const by = V[q + 4] - py; const bz = V[q + 5] - pz;
  const cx = V[q + 6] - px; const cy = V[q + 7] - py; const cz = V[q + 8] - pz;
  const ra = Math.sqrt(ax * ax + ay * ay + az * az); const rb = Math.sqrt(bx * bx + by * by + bz * bz); const rc = Math.sqrt(cx * cx + cy * cy + cz * cz);
  let pot = 0; let gx = 0; let gy = 0; let gz = 0;
  // edge a-b
  let l = L[3 * j]; let den = ra + rb - l; if (den < 1e-14 * l) den = 1e-14 * l;
  let Q = Math.log((ra + rb + l) / den);
  let d = ax * M[q] + ay * M[q + 1] + az * M[q + 2];
  pot += d * Q; gx -= M[q] * Q; gy -= M[q + 1] * Q; gz -= M[q + 2] * Q;
  // edge b-c
  l = L[3 * j + 1]; den = rb + rc - l; if (den < 1e-14 * l) den = 1e-14 * l;
  Q = Math.log((rb + rc + l) / den);
  d = bx * M[q + 3] + by * M[q + 4] + bz * M[q + 5];
  pot += d * Q; gx -= M[q + 3] * Q; gy -= M[q + 4] * Q; gz -= M[q + 5] * Q;
  // edge c-a
  l = L[3 * j + 2]; den = rc + ra - l; if (den < 1e-14 * l) den = 1e-14 * l;
  Q = Math.log((rc + ra + l) / den);
  d = cx * M[q + 6] + cy * M[q + 7] + cz * M[q + 8];
  pot += d * Q; gx -= M[q + 6] * Q; gy -= M[q + 7] * Q; gz -= M[q + 8] * Q;
  let om = 0;
  if (!self) {
    const nx = N[3 * j]; const ny = N[3 * j + 1]; const nz = N[3 * j + 2];
    // signed solid angle: negative when the field point is on the side the normal points to
    const triple = ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + az * (bx * cy - by * cx);
    const dd = ra * rb * rc + (ax * bx + ay * by + az * bz) * rc + (ax * cx + ay * cy + az * cz) * rb + (bx * cx + by * cy + bz * cz) * ra;
    om = 2 * Math.atan2(triple, dd);
    const h = -(ax * nx + ay * ny + az * nz);
    pot += h * om; gx += om * nx; gy += om * ny; gz += om * nz;
  }
  out[0] = pot; out[1] = gx; out[2] = gy; out[3] = gz; out[4] = om;
}

/** In-place LU with partial pivoting on a row-major n×n Float64Array. Returns the pivot vector. */
function luInPlace(a, n) {
  const piv = new Int32Array(n);
  for (let k = 0; k < n; k++) {
    let p = k; let big = Math.abs(a[k * n + k]);
    for (let i = k + 1; i < n; i++) { const v = Math.abs(a[i * n + k]); if (v > big) { big = v; p = i; } }
    if (!(big > 0)) throw new Error('addedMassBem: the influence matrix is singular (duplicate or degenerate panels?)');
    piv[k] = p;
    const rk = k * n;
    if (p !== k) { const rp = p * n; for (let j = 0; j < n; j++) { const t = a[rk + j]; a[rk + j] = a[rp + j]; a[rp + j] = t; } }
    const inv = 1 / a[rk + k];
    for (let i = k + 1; i < n; i++) {
      const ri = i * n; const f = a[ri + k] * inv;
      a[ri + k] = f;
      if (f !== 0) for (let j = k + 1; j < n; j++) a[ri + j] -= f * a[rk + j];
    }
  }
  return piv;
}

function luSolveInPlace(a, piv, n, x) {
  for (let k = 0; k < n; k++) { const p = piv[k]; if (p !== k) { const t = x[k]; x[k] = x[p]; x[p] = t; } }
  for (let i = 1; i < n; i++) { let s = x[i]; const ri = i * n; for (let j = 0; j < i; j++) s -= a[ri + j] * x[j]; x[i] = s; }
  for (let i = n - 1; i >= 0; i--) { let s = x[i]; const ri = i * n; for (let j = i + 1; j < n; j++) s -= a[ri + j] * x[j]; x[i] = s / a[ri + i]; }
}

/** Solves Aᵀ x = b from the factors of P A = L U. */
function luSolveTransposeInPlace(a, piv, n, x) {
  for (let i = 0; i < n; i++) { let s = x[i]; for (let j = 0; j < i; j++) s -= a[j * n + i] * x[j]; x[i] = s / a[i * n + i]; }
  for (let i = n - 1; i >= 0; i--) { let s = x[i]; for (let j = i + 1; j < n; j++) s -= a[j * n + i] * x[j]; x[i] = s; }
  for (let k = n - 1; k >= 0; k--) { const p = piv[k]; if (p !== k) { const t = x[k]; x[k] = x[p]; x[p] = t; } }
}

/** 1-norm condition estimate (Hager 1984; Higham, ACM TOMS 14 (1988)) from the LU factors. */
function conditionEstimate(norm1, lu, piv, n) {
  const x = new Float64Array(n).fill(1 / n); const z = new Float64Array(n);
  let est = 0;
  for (let it = 0; it < 5; it++) {
    luSolveInPlace(lu, piv, n, x);
    let s = 0;
    for (let i = 0; i < n; i++) { s += Math.abs(x[i]); z[i] = x[i] >= 0 ? 1 : -1; }
    if (s <= est && it > 0) break;
    est = s;
    luSolveTransposeInPlace(lu, piv, n, z);
    let jm = 0; let zm = 0; let zx = 0;
    for (let i = 0; i < n; i++) { if (Math.abs(z[i]) > zm) { zm = Math.abs(z[i]); jm = i; } zx += z[i] * x[i]; }
    if (zm <= zx && it > 0) break;
    x.fill(0); x[jm] = 1;
  }
  return norm1 * est;
}

function matVec(a, n, x, y) {
  for (let i = 0; i < n; i++) { const row = i * n; let s = 0; for (let j = 0; j < n; j++) s += a[row + j] * x[j]; y[i] = s; }
}

/**
 * Restarted GMRES (Saad & Schultz 1986) with modified Gram–Schmidt and Givens rotations. The panel-method
 * matrices are second-kind (±2 pi I + compact), so a few tens of iterations reach round-off.
 * Returns true when the relative residual tol was reached; x holds the solution.
 */
function gmres(a, n, b, x, { tol = 1e-11, m = 60, restarts = 4 } = {}) {
  let bn = 0;
  for (let i = 0; i < n; i++) bn += b[i] * b[i];
  bn = Math.sqrt(bn);
  x.fill(0);
  if (bn === 0) return true;
  const Vb = []; const H = new Float64Array((m + 1) * m); const cs = new Float64Array(m); const sn = new Float64Array(m); const g = new Float64Array(m + 1);
  const w = new Float64Array(n); const y = new Float64Array(m);
  for (let rs = 0; rs <= restarts; rs++) {
    matVec(a, n, x, w);
    let beta = 0;
    for (let i = 0; i < n; i++) { w[i] = b[i] - w[i]; beta += w[i] * w[i]; }
    beta = Math.sqrt(beta);
    if (beta <= tol * bn) return true;
    if (!Vb[0]) Vb[0] = new Float64Array(n);
    for (let i = 0; i < n; i++) Vb[0][i] = w[i] / beta;
    g.fill(0); g[0] = beta;
    let k = 0; let done = false;
    for (; k < m && !done; k++) {
      matVec(a, n, Vb[k], w);
      for (let q = 0; q <= k; q++) {
        const v = Vb[q]; let h = 0;
        for (let i = 0; i < n; i++) h += w[i] * v[i];
        H[q * m + k] = h;
        for (let i = 0; i < n; i++) w[i] -= h * v[i];
      }
      let hn = 0;
      for (let i = 0; i < n; i++) hn += w[i] * w[i];
      hn = Math.sqrt(hn);
      if (!Vb[k + 1]) Vb[k + 1] = new Float64Array(n);
      if (hn > 0) for (let i = 0; i < n; i++) Vb[k + 1][i] = w[i] / hn;
      for (let q = 0; q < k; q++) { const t = cs[q] * H[q * m + k] + sn[q] * H[(q + 1) * m + k]; H[(q + 1) * m + k] = -sn[q] * H[q * m + k] + cs[q] * H[(q + 1) * m + k]; H[q * m + k] = t; }
      const hk = H[k * m + k]; const den = Math.hypot(hk, hn);
      cs[k] = hk / den; sn[k] = hn / den;
      H[k * m + k] = den;
      g[k + 1] = -sn[k] * g[k]; g[k] = cs[k] * g[k];
      if (Math.abs(g[k + 1]) <= tol * bn || hn === 0) done = true;
    }
    for (let i = k - 1; i >= 0; i--) { let s = g[i]; for (let j = i + 1; j < k; j++) s -= H[i * m + j] * y[j]; y[i] = s / H[i * m + i]; }
    for (let q = 0; q < k; q++) { const v = Vb[q]; const c = y[q]; for (let i = 0; i < n; i++) x[i] += c * v[i]; }
    if (done) return true;
  }
  return false;
}

/**
 * Added-mass matrix of a hull by a panel method, in the infinite-frequency limit (phi = 0 on z = 0,
 * opposite-sign image), the zero-frequency limit (rigid lid, equal-sign image) or in unbounded fluid.
 *
 * A∞ (limit 'infinite') is the constant added mass of the Cummins time-domain equation; A(0) (limit 'zero') is
 * the low-frequency asymptote used for slow-drift and manoeuvring motions.
 *
 * This is NOT a frequency-dependent radiation/diffraction solver: it returns no wave damping, no A(omega) and
 * no excitation forces. Use the NEMOH / Capytaine / WAMIT hand-off in src/interop/nemoh.js for those.
 *
 * Reference: Morino & Kuo (1974) (formulation 'potential'); Hess & Smith (1964) (formulation 'source');
 * Newman (1977) §4.13 and §6.18; closed-form panel integrals after Newman (1986).
 * Validity: consistently oriented, non-self-intersecting triangle surface, z up, waterline at z = 0, between 4
 * and 6000 wetted panels (dense O(n²) storage: 2 × 8 n² bytes). Flat panels with centroid collocation: on a
 * faceted sphere the 'potential' formulation converges at second order in the panel size and the 'source'
 * formulation at first order. Sharp convex edges (boxes, barges) converge more slowly; refine towards the
 * edges and check with a second mesh. Thin plates (two faces closer than a panel size) are outside the range.
 * Finite depth is supported for limit 'infinite' only (alternating image series about the free surface and the
 * flat seabed): with a rigid lid in finite depth the heave added mass of a surface-piercing body is unbounded,
 * so 'zero' with a finite depth throws. The seabed must lie below the keel.
 *
 * The raw collocation matrix is not exactly symmetric (a discretisation error that vanishes with refinement);
 * the returned A is the symmetric part and `symmetryError` reports the asymmetry that was removed.
 * @param {object} o
 * @param {{positions:ArrayLike<number>, indices:ArrayLike<number>, origin?:number[]}} o.surface hull surface, m
 * @param {number} [o.rho=1025] water density, kg/m³
 * @param {'infinite'|'zero'|'unbounded'} [o.limit='infinite'] frequency limit (see the file header)
 * @param {number[]} [o.cog=[0,0,0]] reference point of the rotations and moments, m
 * @param {number} [o.depth] water depth, m (omit or Infinity for deep water)
 * @param {'potential'|'source'} [o.formulation='potential'] discretisation (see the file header)
 * @param {'gmres'|'lu'} [o.solver='gmres'] linear solver; 'gmres' falls back to LU when it does not converge,
 *   'lu' additionally returns the condition estimate
 * @param {number} [o.imagePairs=12] seabed image pairs kept in finite depth (alternating series, Euler-averaged)
 * @returns {{A:number[][], panels:number, conditionNumber:number|null, symmetryError:number, area:number,
 *   volume:number, limit:string, formulation:string, solver:string}} A: 6×6 added mass about cog in SI (kg for
 *   i, j ≤ 3; kg·m for mixed terms; kg·m² for i, j > 3), symmetric; panels: number of wetted panels;
 *   conditionNumber: 1-norm estimate of the influence matrix (-), null unless the LU solver ran;
 *   symmetryError: max |A_ij - A_ji| / sqrt(A_ii A_jj) of the raw matrix (-); wetted area m² and displaced
 *   volume m³ of the panel model; the solver actually used
 */
export function addedMassBem({ surface, rho = 1025, limit = 'infinite', cog = [0, 0, 0], depth, formulation = 'potential', solver = 'gmres', imagePairs = 12 } = {}) {
  if (!(rho > 0) || !Number.isFinite(rho)) throw new RangeError('rho must be > 0');
  if (!['infinite', 'zero', 'unbounded'].includes(limit)) throw new RangeError("limit must be 'infinite', 'zero' or 'unbounded'");
  if (!['potential', 'source'].includes(formulation)) throw new RangeError("formulation must be 'potential' or 'source'");
  if (!['gmres', 'lu'].includes(solver)) throw new RangeError("solver must be 'gmres' or 'lu'");
  if (!cog || cog.length !== 3 || !(Number.isFinite(cog[0]) && Number.isFinite(cog[1]) && Number.isFinite(cog[2]))) throw new RangeError('cog must be [x, y, z] in metres');
  const finiteDepth = depth !== undefined && depth !== null && Number.isFinite(depth);
  if (depth !== undefined && depth !== null && !(depth > 0)) throw new RangeError('depth must be > 0');
  if (finiteDepth && limit !== 'infinite') throw new RangeError("depth: finite depth is only supported with limit 'infinite' (the rigid-lid heave added mass is unbounded in finite depth)");
  if (!Number.isInteger(imagePairs) || imagePairs < 1) throw new RangeError('imagePairs must be a positive integer');
  const pan = wettedPanels(surface);
  const n = pan.n;
  if (n > MAX_PANELS) throw new RangeError(`surface has ${n} wetted panels; the dense solver accepts at most ${MAX_PANELS} (coarsen the mesh)`);
  if (n < 4) throw new RangeError('surface must have at least 4 wetted panels');
  const C = pan.centroids; const N = pan.normals; const Ar = pan.areas;
  if (finiteDepth) {
    let zmin = 0;
    for (let i = 2; i < pan.vertices.length; i += 3) if (pan.vertices[i] < zmin) zmin = pan.vertices[i];
    if (!(depth > -zmin)) throw new RangeError('depth must exceed the draught of the hull');
  }
  const geo = panelGeometry(pan);
  const s = limit === 'infinite' ? -1 : limit === 'zero' ? 1 : 0;
  const src = formulation === 'source';
  // K: system matrix; S: single-layer (potential of unit source density) matrix
  const K = new Float64Array(n * n); const S = new Float64Array(n * n);
  const r = new Float64Array(5);
  // seabed image copies beyond the first: z-shift and sign; the last pair carries weight 1/2 (Euler averaging of
  // the last two partial sums of the alternating series)
  const farShift = []; const farSign = [];
  if (finiteDepth) {
    for (let m = 1; m <= imagePairs + 1; m++) {
      const w = m <= imagePairs ? 1 : 0.5;
      const sg = m % 2 === 0 ? 1 : -1;
      farShift.push(2 * m * depth, -2 * m * depth); farSign.push(sg * w, sg * w);
    }
  }
  for (let i = 0; i < n; i++) {
    const px = C[3 * i]; const py = C[3 * i + 1]; const pz = C[3 * i + 2];
    const nx = N[3 * i]; const ny = N[3 * i + 1]; const nz = N[3 * i + 2];
    const row = i * n;
    for (let j = 0; j < n; j++) {
      panelIntegral(pan, geo, j, px, py, pz, i === j, r);
      let pot = r[0];
      // kernel: normal velocity of a unit source panel ('source') or signed solid angle ('potential')
      let ker = src ? (i === j ? -2 * Math.PI : r[1] * nx + r[2] * ny + r[3] * nz) : (i === j ? 2 * Math.PI : r[4]);
      if (s !== 0) {
        // image of panel j about z = 0, evaluated by reflecting the field point instead
        panelIntegral(pan, geo, j, px, py, -pz, false, r);
        pot += s * r[0];
        ker += s * (src ? r[1] * nx + r[2] * ny - r[3] * nz : r[4]);
      }
      if (finiteDepth) {
        const cx = C[3 * j]; const cy = C[3 * j + 1]; const cz = C[3 * j + 2];
        const mx = N[3 * j]; const my = N[3 * j + 1]; const mz = N[3 * j + 2];
        const near2 = 64 * geo.diam[j] * geo.diam[j];
        for (let f = 0; f < farShift.length; f++) {
          const sh = farShift[f]; const sg = farSign[f];
          // translated copy (source points at zeta + shift, sign sg), then mirrored copy (-zeta + shift, sign -sg);
          // both are evaluated on the original panel at a moved field point, as a point source/dipole when far
          for (let mir = 0; mir < 2; mir++) {
            const fz = mir === 0 ? pz - sh : sh - pz;
            const zs = mir === 0 ? 1 : -1;
            const sign = mir === 0 ? sg : -sg;
            const dx = px - cx; const dy = py - cy; const dz = fz - cz;
            const d2 = dx * dx + dy * dy + dz * dz;
            if (d2 > near2) {
              const d1 = Math.sqrt(d2); const g3 = Ar[j] / (d2 * d1);
              pot += sign * Ar[j] / d1;
              ker += sign * (src ? -g3 * (dx * nx + dy * ny + zs * dz * nz) : -g3 * (dx * mx + dy * my + dz * mz));
            } else {
              panelIntegral(pan, geo, j, px, py, fz, false, r);
              pot += sign * r[0];
              ker += sign * (src ? r[1] * nx + r[2] * ny + zs * r[3] * nz : r[4]);
            }
          }
        }
      }
      K[row + j] = ker; S[row + j] = pot;
    }
  }
  // generalised normals of the six rigid-body modes
  const gn = [];
  for (let k = 0; k < 6; k++) gn.push(new Float64Array(n));
  for (let i = 0; i < n; i++) {
    const rx = C[3 * i] - cog[0]; const ry = C[3 * i + 1] - cog[1]; const rz = C[3 * i + 2] - cog[2];
    const nx = N[3 * i]; const ny = N[3 * i + 1]; const nz = N[3 * i + 2];
    gn[0][i] = nx; gn[1][i] = ny; gn[2][i] = nz;
    gn[3][i] = ry * nz - rz * ny; gn[4][i] = rz * nx - rx * nz; gn[5][i] = rx * ny - ry * nx;
  }
  const rhs = []; const sol = [];
  for (let k = 0; k < 6; k++) {
    const b = new Float64Array(n);
    if (src) b.set(gn[k]); else { matVec(S, n, gn[k], b); for (let i = 0; i < n; i++) b[i] = -b[i]; }
    rhs.push(b); sol.push(new Float64Array(n));
  }
  let used = solver; let conditionNumber = null;
  if (solver === 'gmres') {
    for (let k = 0; k < 6 && used === 'gmres'; k++) if (!gmres(K, n, rhs[k], sol[k])) used = 'lu';
  }
  if (used === 'lu') {
    let norm1 = 0;
    const colSum = new Float64Array(n);
    for (let i = 0; i < n; i++) { const row = i * n; for (let j = 0; j < n; j++) colSum[j] += Math.abs(K[row + j]); }
    for (let j = 0; j < n; j++) if (colSum[j] > norm1) norm1 = colSum[j];
    const piv = luInPlace(K, n);
    for (let k = 0; k < 6; k++) { sol[k].set(rhs[k]); luSolveInPlace(K, piv, n, sol[k]); }
    conditionNumber = conditionEstimate(norm1, K, piv, n);
  }
  const raw = Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]);
  const phi = new Float64Array(n);
  for (let k = 0; k < 6; k++) {
    if (src) matVec(S, n, sol[k], phi); else phi.set(sol[k]);
    for (let i = 0; i < n; i++) {
      const w = -rho * phi[i] * Ar[i];
      for (let m = 0; m < 6; m++) raw[m][k] += w * gn[m][i];
    }
  }
  let symmetryError = 0;
  const A = Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]);
  for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) {
    A[i][j] = 0.5 * (raw[i][j] + raw[j][i]);
    if (!Number.isFinite(A[i][j])) throw new Error('addedMassBem: the solution is not finite (check the mesh for coincident or degenerate panels)');
    if (j > i) { const sc = Math.sqrt(Math.abs(raw[i][i] * raw[j][j])); if (sc > 0) symmetryError = Math.max(symmetryError, Math.abs(raw[i][j] - raw[j][i]) / sc); }
  }
  return { A, panels: n, conditionNumber, symmetryError, area: pan.area, volume: pan.volume, limit, formulation, solver: used };
}
