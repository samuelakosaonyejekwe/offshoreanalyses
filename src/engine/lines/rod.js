/**
 * Planar (2-D) geometrically nonlinear slender-rod statics and linearised dynamics.
 *
 * Formulation: two-node corotational Euler–Bernoulli beam elements (Crisfield, "Non-linear Finite Element
 * Analysis of Solids and Structures", Vol. 1, Ch. 7), exact for arbitrarily large rigid-body rotation of each
 * element with small strains in the corotated frame. The rod is formulated in *effective* tension
 * Te = Tw − pi·Ai + pe·Ae (Sparks, "Fundamentals of Marine Riser Mechanics", Ch. 2–3): the distributed load
 * is the apparent (submerged) weight and the axial strain is (Te + (1 − 2ν)(pi·Ai − pe·Ae))/EA.
 *
 * Coordinates: x horizontal, z vertical upward, angles in radians measured from +x towards +z. Arc length s
 * runs from end A (s = 0) to end B. The model is strictly two-dimensional: out-of-plane current, torsion and
 * lateral seabed friction are not represented (3-D is not implemented).
 */

import { subspaceIteration } from '../math/eigen.js';

const G0 = 9.80665;

/* ------------------------------------------------------------------------------------------------------- */
/* Banded matrix with room for partial-pivoting fill: row-major, entry (i, j) at a[i*W + j - i + kl].       */

function bandCreate(n, hb) {
  const W = 3 * hb + 1;
  return { n, hb, W, a: new Float64Array(n * W), piv: new Int32Array(n) };
}
function bandAdd(B, i, j, v) { B.a[i * B.W + j - i + B.hb] += v; }
function bandGet(B, i, j) { return Math.abs(i - j) > B.hb ? 0 : B.a[i * B.W + j - i + B.hb]; }
function bandSet(B, i, j, v) { B.a[i * B.W + j - i + B.hb] = v; }

function bandFactor(B) {
  const { n, hb, W, a, piv } = B;
  for (let k = 0; k < n; k++) {
    const rmax = Math.min(n - 1, k + hb);
    const cmax = Math.min(n - 1, k + 2 * hb);
    let p = k;
    let best = Math.abs(a[k * W + hb]);
    for (let r = k + 1; r <= rmax; r++) {
      const v = Math.abs(a[r * W + k - r + hb]);
      if (v > best) { best = v; p = r; }
    }
    if (!(best > 1e-300)) throw new Error(`rod: singular stiffness at degree of freedom ${k} (unrestrained mechanism?)`);
    piv[k] = p;
    if (p !== k) {
      for (let j = k; j <= cmax; j++) {
        const ik = k * W + j - k + hb;
        const ip = p * W + j - p + hb;
        const t = a[ik]; a[ik] = a[ip]; a[ip] = t;
      }
    }
    const d = a[k * W + hb];
    for (let r = k + 1; r <= rmax; r++) {
      const ir = r * W - r + hb;
      const f = a[ir + k] / d;
      a[ir + k] = f;
      if (f !== 0) {
        const ikr = k * W - k + hb;
        for (let j = k + 1; j <= cmax; j++) a[ir + j] -= f * a[ikr + j];
      }
    }
  }
}
function bandSolve(B, b) {
  const { n, hb, W, a, piv } = B;
  for (let k = 0; k < n; k++) {
    const p = piv[k];
    if (p !== k) { const t = b[k]; b[k] = b[p]; b[p] = t; }
    const bk = b[k];
    if (bk !== 0) {
      const rmax = Math.min(n - 1, k + hb);
      for (let r = k + 1; r <= rmax; r++) b[r] -= a[r * W + k - r + hb] * bk;
    }
  }
  for (let k = n - 1; k >= 0; k--) {
    const cmax = Math.min(n - 1, k + 2 * hb);
    const ik = k * W - k + hb;
    let s = b[k];
    for (let j = k + 1; j <= cmax; j++) s -= a[ik + j] * b[j];
    b[k] = s / a[ik + k];
  }
  return b;
}

/* ------------------------------------------------------------------------------------------------------- */

function num(v, name, { min = -Infinity, max = Infinity, strictMin = false } = {}) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < min || v > max || (strictMin && v <= min)) {
    throw new RangeError(`rod: argument ${name} = ${v} is outside its valid range`);
  }
  return v;
}

function wrapPi(a) { return Math.atan2(Math.sin(a), Math.cos(a)); }

function buildMesh(model) {
  const { sections } = model;
  if (!Array.isArray(sections) || sections.length === 0) throw new RangeError('rod: argument sections must be a non-empty array');
  let Ltot = 0;
  sections.forEach((sec, k) => {
    num(sec.length, `sections[${k}].length`, { min: 0, strictMin: true });
    num(sec.EA, `sections[${k}].EA`, { min: 0, strictMin: true });
    num(sec.EI, `sections[${k}].EI`, { min: 0, strictMin: true });
    num(sec.w, `sections[${k}].w`);
    Ltot += sec.length;
  });
  const nTot = Math.round(num(model.nElem ?? 100, 'nElem', { min: 1, max: 20000 }));
  const counts = sections.map((sec) => Math.max(1, Math.round(sec.nElem ?? (nTot * sec.length) / Ltot)));
  const ne = counts.reduce((a, b) => a + b, 0);
  const mesh = {
    ne, nn: ne + 1, ndof: 3 * (ne + 1), Ltot,
    L0: new Float64Array(ne), EA: new Float64Array(ne), EI: new Float64Array(ne), w: new Float64Array(ne),
    D: new Float64Array(ne), Cdn: new Float64Array(ne), Cdt: new Float64Array(ne),
    Ai: new Float64Array(ne), Ae: new Float64Array(ne), m: new Float64Array(ne), ma: new Float64Array(ne),
    maT: new Float64Array(ne), sec: new Int32Array(ne), s: new Float64Array(ne + 1),
    trib: new Float64Array(ne + 1), P0: new Float64Array(ne),
    dragX: new Float64Array(ne + 1), dragZ: new Float64Array(ne + 1),
    offset: new Float64Array(ne + 1),
  };
  let e = 0;
  sections.forEach((sec, k) => {
    const h = sec.length / counts[k];
    for (let i = 0; i < counts[k]; i++, e++) {
      mesh.L0[e] = h; mesh.EA[e] = sec.EA; mesh.EI[e] = sec.EI; mesh.w[e] = sec.w;
      mesh.D[e] = sec.D ?? 0; mesh.Cdn[e] = sec.Cdn ?? 0; mesh.Cdt[e] = sec.Cdt ?? 0;
      mesh.Ai[e] = sec.Ai ?? 0; mesh.Ae[e] = sec.Ae ?? 0;
      mesh.m[e] = sec.m ?? 0; mesh.ma[e] = sec.ma ?? 0; mesh.maT[e] = sec.maT ?? 0;
      mesh.sec[e] = k;
      mesh.s[e + 1] = mesh.s[e] + h;
      mesh.trib[e] += h / 2; mesh.trib[e + 1] += h / 2;
      const off = sec.contactOffset ?? 0;
      mesh.offset[e] = Math.max(mesh.offset[e], off); mesh.offset[e + 1] = off;
    }
  });
  mesh.hRef = Ltot / ne;
  return mesh;
}

function nearestNode(mesh, s) {
  let lo = 0, hi = mesh.nn - 1;
  while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (mesh.s[mid] <= s) lo = mid; else hi = mid; }
  return (s - mesh.s[lo] <= mesh.s[hi] - s) ? lo : hi;
}

/* Lumped (bar-chain) cable shape between two pinned ends, used only as a starting configuration. */
function cableGuess(mesh, model, xA, zA, xB, zB) {
  const ne = mesh.ne;
  const zs = model.seabed ? model.seabed.z : null;
  const tolZ = 1e-6 * mesh.Ltot + 1e-9;
  const groundA = zs !== null && Math.abs(zA - zs) <= tolZ;
  const groundB = zs !== null && Math.abs(zB - zs) <= tolZ && !groundA;
  const rev = groundB;                      // march from the grounded end
  const idx = (i) => (rev ? ne - 1 - i : i);
  const x0 = rev ? xB : xA, z0 = rev ? zB : zA, x1 = rev ? xA : xB, z1 = rev ? zA : zB;
  const dirx = x1 >= x0 ? 1 : -1;
  const X = Math.abs(x1 - x0), h = z1 - z0;
  const grounded = groundA || groundB;
  let Wabs = 0;
  for (let e = 0; e < ne; e++) Wabs += Math.abs(mesh.w[e]) * mesh.L0[e];
  Wabs = Math.max(Wabs, 1e-9);
  const px = new Float64Array(ne + 1), pz = new Float64Array(ne + 1), ang = new Float64Array(ne);
  const wt = mesh.w;
  const march = (H, V0) => {
    let V = V0, x = 0, z = 0, lifted = !grounded;
    px[0] = 0; pz[0] = 0;
    for (let i = 0; i < ne; i++) {
      const e = idx(i);
      const Wl = wt[e] * mesh.L0[e];
      if (i > 0) V += 0.5 * (Wl + wt[idx(i - 1)] * mesh.L0[idx(i - 1)]); else V += 0.5 * Wl;
      let Ve = V;
      if (!lifted) { if (Ve <= 0) Ve = 0; else lifted = true; }
      const T = Math.hypot(H, Ve);
      const Ls = mesh.L0[e] * (1 + T / mesh.EA[e]);
      x += Ls * H / T; z += Ls * Ve / T;
      ang[i] = Math.atan2(Ve, H);
      px[i + 1] = x; pz[i + 1] = z;
    }
    return [x - X, z - h];
  };
  // Nested bisection: for a given H the end elevation rises monotonically with the end vertical force V0, and
  // with V0(H) so determined the horizontal span grows monotonically with H.
  const solveV = (H) => {
    let lo = -Wabs - 1, hi = 2 * Wabs + 1e4 * H;
    for (let k = 0; k < 70; k++) {
      const mid = 0.5 * (lo + hi);
      if (march(H, mid)[1] > 0) hi = mid; else lo = mid;
    }
    return 0.5 * (lo + hi);
  };
  let lo = Math.log(1e-8 * Wabs), hi = Math.log(1e5 * Wabs);
  for (let k = 0; k < 60; k++) {
    const mid = 0.5 * (lo + hi), H = Math.exp(mid);
    if (march(H, solveV(H))[0] > 0) hi = mid; else lo = mid;
  }
  const Hs = Math.exp(0.5 * (lo + hi));
  march(Hs, solveV(Hs));
  const u = new Float64Array(mesh.ndof);
  for (let i = 0; i <= ne; i++) {
    const node = rev ? ne - i : i;
    u[3 * node] = x0 + dirx * px[i];
    u[3 * node + 1] = z0 + pz[i];
  }
  fillAngles(mesh, u);
  return u;
}

function fillAngles(mesh, u) {
  const ne = mesh.ne;
  let prev = 0;
  for (let i = 0; i <= ne; i++) {
    let b0 = null, b1 = null;
    if (i > 0) b0 = Math.atan2(u[3 * i + 1] - u[3 * i - 2], u[3 * i] - u[3 * i - 3]);
    if (i < ne) b1 = Math.atan2(u[3 * i + 4] - u[3 * i + 1], u[3 * i + 3] - u[3 * i]);
    let th;
    if (b0 === null) th = b1; else if (b1 === null) th = b0; else th = b0 + 0.5 * wrapPi(b1 - b0);
    if (i > 0) th = prev + wrapPi(th - prev);
    u[3 * i + 2] = th; prev = th;
  }
}

function currentAt(current, z) {
  if (!current) return 0;
  if (typeof current.profile === 'function') return current.profile(z);
  if (current.z && current.U) {
    const zs = current.z, U = current.U, n = zs.length;
    if (z <= zs[0]) return U[0];
    if (z >= zs[n - 1]) return U[n - 1];
    let k = 0; while (k < n - 2 && zs[k + 1] < z) k++;
    return U[k] + (U[k + 1] - U[k]) * (z - zs[k]) / (zs[k + 1] - zs[k]);
  }
  return current.U ?? 0;
}

function updateLaggedLoads(mesh, model, u) {
  const { ne } = mesh;
  const cur = model.current, pr = model.pressure;
  mesh.dragX.fill(0); mesh.dragZ.fill(0);
  for (let e = 0; e < ne; e++) {
    const i = 3 * e;
    const dx = u[i + 3] - u[i], dz = u[i + 4] - u[i + 1];
    const Ln = Math.hypot(dx, dz), c = dx / Ln, s = dz / Ln;
    const zm = 0.5 * (u[i + 1] + u[i + 4]);
    if (cur && mesh.D[e] > 0) {
      const rho = cur.rho ?? 1025;
      const U = currentAt(cur, zm);
      const Un = -s * U, Ut = c * U;
      const fn = 0.5 * rho * mesh.Cdn[e] * mesh.D[e] * Math.abs(Un) * Un;
      const ft = 0.5 * rho * mesh.Cdt[e] * Math.PI * mesh.D[e] * Math.abs(Ut) * Ut;
      const fx = (-s * fn + c * ft) * Ln * 0.5, fz = (c * fn + s * ft) * Ln * 0.5;
      mesh.dragX[e] += fx; mesh.dragX[e + 1] += fx; mesh.dragZ[e] += fz; mesh.dragZ[e + 1] += fz;
    }
    if (pr) {
      const g = pr.g ?? G0;
      const pi = (pr.pIntRef ?? 0) + (pr.rhoInt ?? 0) * g * ((pr.zIntRef ?? 0) - zm);
      const pe = Math.max(0, (pr.rhoExt ?? 1025) * g * ((pr.zSurface ?? 0) - zm));
      mesh.P0[e] = (1 - 2 * (pr.nu ?? 0.3)) * (Math.max(pi, 0) * mesh.Ai[e] - pe * mesh.Ae[e]);
    }
  }
}

// element scratch vectors (r, q and the two rotation rows of the corotational B matrix)
const SCR_R = new Float64Array(6), SCR_Q = new Float64Array(6);
const SCR_B1 = Float64Array.of(0, 0, 1, 0, 0, 0), SCR_B2 = Float64Array.of(0, 0, 0, 0, 0, 1);

function seabedZ(sb, x) { return typeof sb.profile === 'function' ? sb.profile(x) : sb.z; }
// Axial shear stiffness of the seabed before sliding [Pa]. Friction needs an elastic stage; when a friction
// coefficient is given without kt the normal stiffness is used (mobilisation over a slip of mu × penetration).
function seabedKt(sb) { return sb.kt !== undefined ? sb.kt : sb.mu > 0 ? sb.kn : 0; }

/**
 * Residual (internal − external) and tangent stiffness. `out` (optional) receives element forces and contact
 * reactions for post-processing.
 */
function assemble(mesh, model, st, u, fac, R, B, out) {
  const { ne, nn } = mesh;
  R.fill(0);
  if (B) B.a.fill(0);
  for (let e = 0; e < ne; e++) {
    const i = 3 * e;
    const dx = u[i + 3] - u[i], dz = u[i + 4] - u[i + 1];
    const Ln = Math.hypot(dx, dz), c = dx / Ln, s = dz / Ln;
    const beta = Math.atan2(dz, dx);
    const t1 = wrapPi(u[i + 2] - beta), t2 = wrapPi(u[i + 5] - beta);
    const L0 = mesh.L0[e], ka = mesh.EA[e] / L0, kb = mesh.EI[e] / L0;
    const N = ka * (Ln - L0) - mesh.P0[e];
    const M1 = kb * (4 * t1 + 2 * t2), M2 = kb * (2 * t1 + 4 * t2);
    const Q = (M1 + M2) / Ln;
    if (Math.abs(N) > fac.Nmax) fac.Nmax = Math.abs(N);
    // r = [-c,-s,0,c,s,0], q = [s,-c,0,-s,c,0]
    R[i] += -N * c - Q * s; R[i + 1] += -N * s + Q * c; R[i + 2] += M1;
    R[i + 3] += N * c + Q * s; R[i + 4] += N * s - Q * c; R[i + 5] += M2;
    if (out) { out.N[e] = N; out.M1[e] = M1; out.M2[e] = M2; out.Ln[e] = Ln; }
    if (B) {
      const r = SCR_R, q = SCR_Q, b1 = SCR_B1, b2 = SCR_B2;
      r[0] = -c; r[1] = -s; r[3] = c; r[4] = s;
      q[0] = s; q[1] = -c; q[3] = -s; q[4] = c;
      b1[0] = -s / Ln; b1[1] = c / Ln; b1[3] = s / Ln; b1[4] = -c / Ln;
      b2[0] = b1[0]; b2[1] = b1[1]; b2[3] = b1[3]; b2[4] = b1[4];
      const g1 = N / Ln, g2 = (M1 + M2) / (Ln * Ln);
      for (let a = 0; a < 6; a++) {
        for (let b = 0; b < 6; b++) {
          const v = ka * r[a] * r[b]
            + kb * (4 * b1[a] * b1[b] + 2 * b1[a] * b2[b] + 2 * b2[a] * b1[b] + 4 * b2[a] * b2[b])
            + g1 * q[a] * q[b] + g2 * (r[a] * q[b] + q[a] * r[b]);
          bandAdd(B, i + a, i + b, v);
        }
      }
    }
    // apparent weight (dead load, lumped)
    const Wh = 0.5 * fac.weight * (st.w0 ? st.w0[e] + fac.lam * (mesh.w[e] - st.w0[e]) : mesh.w[e]) * L0;
    R[i + 1] += Wh; R[i + 4] += Wh;
  }
  for (let k = 0; k < nn; k++) {
    R[3 * k] -= fac.load * (st.Fx[k] + mesh.dragX[k]);
    R[3 * k + 1] -= fac.load * (st.Fz[k] + mesh.dragZ[k]);
    R[3 * k + 2] -= fac.load * st.Mn[k];
  }
  if (out) out.q.fill(0);
  const sb = model.seabed;
  if (sb) {
    const kn = sb.kn, mu = sb.mu ?? 0, kt = seabedKt(sb);
    for (let k = 0; k < nn; k++) {
      const x = u[3 * k], z = u[3 * k + 1];
      const zs = seabedZ(sb, x) + mesh.offset[k];
      const pen = zs - z;
      if (pen < 0) continue;
      const kl = kn * mesh.trib[k];
      const Rn = kl * pen;
      R[3 * k + 1] -= Rn;
      if (out) out.q[k] += kn * pen;
      let slope = 0;
      if (typeof sb.profile === 'function') { const hx = 1e-4 * mesh.hRef + 1e-6; slope = (sb.profile(x + hx) - sb.profile(x - hx)) / (2 * hx); }
      if (B) { bandAdd(B, 3 * k + 1, 3 * k + 1, kl); if (slope !== 0) bandAdd(B, 3 * k + 1, 3 * k, -kl * slope); }
      if (mu > 0 && kt > 0) {
        if (Number.isNaN(st.anchor[k])) st.anchor[k] = x;       // first contact: friction spring starts unloaded
        const ks = kt * mesh.trib[k];
        const ft = ks * (x - st.anchor[k]);
        if (Math.abs(ft) <= mu * Rn) {
          R[3 * k] += ft;
          if (B) bandAdd(B, 3 * k, 3 * k, ks);
          if (out) out.fr[k] = -ft / mesh.trib[k];
        } else {
          const sg = Math.sign(ft);
          R[3 * k] += sg * mu * Rn;
          if (B && !fac.symmetric) bandAdd(B, 3 * k, 3 * k + 1, -sg * mu * kl);
          if (out) out.fr[k] = -sg * mu * Rn / mesh.trib[k];
        }
      }
    }
  }
  if (model.supports) {
    for (const sp of model.supports) {
      const inside = sp.side === 'inside';
      for (let k = 0; k < nn; k++) {
        const px = u[3 * k] - sp.xc, pz = u[3 * k + 1] - sp.zc;
        const d = Math.hypot(px, pz);
        if (d === 0) continue;
        const Rc = sp.R + (inside ? -1 : 1) * ((sp.offset ?? 0) + mesh.offset[k]);
        const pen = inside ? d - Rc : Rc - d;
        if (pen <= 0) continue;
        // Angular window: no contact outside it, and the stiffness rises smoothly from zero to its full value over
        // the arc length `edge` inside each end (ph = taper factor with zero slope at both ends of the taper, dph
        // and ddph its first and second derivatives with respect to the polar angle). Without the taper a node
        // sliding past the end of the support would lose its whole reaction at once, and the discrete problem
        // would have no equilibrium when the rod bears on the end. The contact force derives from the energy
        // ½·k·ph·pen², so the tangent is symmetric.
        let ph = 1, dph = 0, ddph = 0;
        if (sp.angMin !== undefined || sp.angMax !== undefined) {
          const a = Math.atan2(pz, px);
          const dLo = sp.angMin !== undefined ? (a - sp.angMin) * sp.R : Infinity;
          const dHi = sp.angMax !== undefined ? (sp.angMax - a) * sp.R : Infinity;
          const dIn = Math.min(dLo, dHi);
          if (dIn <= 0) continue;
          const edge = sp.edge ?? mesh.trib[k];
          if (dIn < edge) {
            const ta = dIn / edge, da = (dLo < dHi ? 1 : -1) * sp.R / edge;
            ph = ta * ta * (3 - 2 * ta); dph = 6 * ta * (1 - ta) * da; ddph = (6 - 12 * ta) * da * da;
          }
        }
        const kl = fac.support * sp.k * mesh.trib[k];
        const sgn = inside ? -1 : 1;
        const erx = px / d, erz = pz / d;
        const nx = sgn * erx, nz = sgn * erz;
        R[3 * k] -= kl * ph * pen * nx; R[3 * k + 1] -= kl * ph * pen * nz;
        if (out) out.q[k] += fac.support * sp.k * ph * pen;
        if (B) {
          const g = -sgn * pen / d;
          bandAdd(B, 3 * k, 3 * k, kl * ph * (nx * nx + g * (1 - nx * nx)));
          bandAdd(B, 3 * k, 3 * k + 1, kl * ph * (nx * nz - g * nx * nz));
          bandAdd(B, 3 * k + 1, 3 * k, kl * ph * (nx * nz - g * nx * nz));
          bandAdd(B, 3 * k + 1, 3 * k + 1, kl * ph * (nz * nz + g * (1 - nz * nz)));
        }
        if (dph !== 0) {
          // gradient of the polar angle: ea/d with ea = (−erz, erx); its Hessian is −(er·eaᵀ + ea·erᵀ)/d²
          const ax = -erz / d, az = erx / d, c1 = kl * pen * dph, c2 = 0.5 * kl * pen * pen * dph, c3 = 0.5 * kl * pen * pen * ddph;
          R[3 * k] += c2 * ax; R[3 * k + 1] += c2 * az;
          if (B) {
            bandAdd(B, 3 * k, 3 * k, -2 * c1 * nx * ax - 2 * c2 * erx * ax + c3 * ax * ax);
            bandAdd(B, 3 * k, 3 * k + 1, -c1 * (nx * az + ax * nz) - c2 * (erx * az + ax * erz) + c3 * ax * az);
            bandAdd(B, 3 * k + 1, 3 * k, -c1 * (nx * az + ax * nz) - c2 * (erx * az + ax * erz) + c3 * ax * az);
            bandAdd(B, 3 * k + 1, 3 * k + 1, -2 * c1 * nz * az - 2 * c2 * erz * az + c3 * az * az);
          }
        }
      }
    }
  }
  for (const rs of st.rotSprings) {
    R[rs.dof] += rs.k * (u[rs.dof] - rs.th0);
    if (B) bandAdd(B, rs.dof, rs.dof, rs.k);
  }
}

/**
 * Total potential energy [J] whose gradient is the residual of `assemble`, for a frictionless rod on a level seabed
 * (the only case in which the model is conservative): element axial and bending strain energy, weight and load
 * potentials, and the contact-spring energies. `acc[0]` receives the sum of the absolute values of the terms, a
 * measure of the rounding noise of the result.
 */
function potentialEnergy(mesh, model, st, u, fac, acc) {
  const { ne, nn } = mesh;
  let E = 0, Ea = 0;
  for (let e = 0; e < ne; e++) {
    const i = 3 * e;
    const dx = u[i + 3] - u[i], dz = u[i + 4] - u[i + 1];
    const Ln = Math.hypot(dx, dz), beta = Math.atan2(dz, dx);
    const t1 = wrapPi(u[i + 2] - beta), t2 = wrapPi(u[i + 5] - beta);
    const L0 = mesh.L0[e], dL = Ln - L0;
    const Es = 0.5 * (mesh.EA[e] / L0) * dL * dL + 2 * (mesh.EI[e] / L0) * (t1 * t1 + t1 * t2 + t2 * t2);
    const Wh = 0.5 * fac.weight * (st.w0 ? st.w0[e] + fac.lam * (mesh.w[e] - st.w0[e]) : mesh.w[e]) * L0;
    const Ew = Wh * (u[i + 1] + u[i + 4]) - mesh.P0[e] * dL;
    E += Es + Ew; Ea += Es + Math.abs(Ew);
  }
  for (let k = 0; k < nn; k++) {
    const El = -fac.load * ((st.Fx[k] + mesh.dragX[k]) * u[3 * k] + (st.Fz[k] + mesh.dragZ[k]) * u[3 * k + 1] + st.Mn[k] * u[3 * k + 2]);
    E += El; Ea += Math.abs(El);
  }
  let Ec = 0;
  const sb = model.seabed;
  if (sb) {
    for (let k = 0; k < nn; k++) {
      const pen = sb.z + mesh.offset[k] - u[3 * k + 1];
      if (pen > 0) Ec += 0.5 * sb.kn * mesh.trib[k] * pen * pen;
    }
  }
  for (const sp of model.supports ?? []) {
    const inside = sp.side === 'inside';
    for (let k = 0; k < nn; k++) {
      const px = u[3 * k] - sp.xc, pz = u[3 * k + 1] - sp.zc;
      const d = Math.hypot(px, pz);
      if (d === 0) continue;
      const Rc = sp.R + (inside ? -1 : 1) * ((sp.offset ?? 0) + mesh.offset[k]);
      const pen = inside ? d - Rc : Rc - d;
      if (pen <= 0) continue;
      let ph = 1;
      if (sp.angMin !== undefined || sp.angMax !== undefined) {
        const a = Math.atan2(pz, px);
        const dIn = Math.min(sp.angMin !== undefined ? (a - sp.angMin) * sp.R : Infinity, sp.angMax !== undefined ? (sp.angMax - a) * sp.R : Infinity);
        if (dIn <= 0) continue;
        const ta = Math.min(1, dIn / (sp.edge ?? mesh.trib[k]));
        ph = ta * ta * (3 - 2 * ta);
      }
      Ec += 0.5 * fac.support * sp.k * mesh.trib[k] * ph * pen * pen;
    }
  }
  for (const rs of st.rotSprings) { const dth = u[rs.dof] - rs.th0; Ec += 0.5 * rs.k * dth * dth; }
  if (acc) acc[0] = Ea + Ec;
  return E + Ec;
}

function applyConstraint(B, rhs, c, d) {
  const lo = Math.max(0, c - B.hb), hi = Math.min(B.n - 1, c + B.hb);
  for (let i = lo; i <= hi; i++) {
    if (i === c) continue;
    if (d !== 0) rhs[i] -= bandGet(B, i, c) * d;
    bandSet(B, i, c, 0); bandSet(B, c, i, 0);
  }
  bandSet(B, c, c, 1);
  rhs[c] = d;
}

function endSpec(end, name) {
  const e = end ?? {};
  for (const k of ['x', 'z', 'angle', 'kRot', 'Fx', 'Fz', 'M']) if (e[k] !== undefined) num(e[k], `${name}.${k}`);
  if (e.kRot !== undefined && e.angle === undefined) throw new RangeError(`rod: argument ${name}.kRot needs ${name}.angle (spring reference angle)`);
  return e;
}

/**
 * Static equilibrium of a planar slender rod with bending stiffness under apparent weight, point loads, steady
 * current drag, seabed and rigid circular support contact.
 *
 * Method: corotational Euler–Bernoulli beam elements, Newton–Raphson with step limiting and adaptive
 * load/displacement stepping (restart from the last converged state, with a damped line-search iteration as the
 * fallback when plain Newton steps do not settle); Morison drag and the pressure end-cap strain term are lagged
 * in an outer fixed-point loop. Reference: Crisfield (1991) Vol. 1 Ch. 7; Sparks (2007) Ch. 2–3 (effective tension);
 * DNV-RP-C205 §6.2 (Morison normal/tangential drag).
 * Validity: planar, small strain, linear elastic; rod stress-free shape is straight; seabed reaction is a
 * linear unilateral (lift-off) Winkler spring acting vertically; supports are frictionless.
 *
 * @param {object} model
 * @param {Array<object>} model.sections ordered from end A to end B; each
 *   `{ length [m], EA [N], EI [N·m²], w [N/m apparent (submerged) weight per unstretched length, + down, negative
 *   for buoyant sections], nElem? [-], D? [m drag diameter], Cdn? [-], Cdt? [- on wetted area πD], Ai?, Ae? [m²
 *   internal / external cross-section areas], m? [kg/m mass incl. contents], ma? [kg/m normal added mass],
 *   maT? [kg/m tangential added mass], contactOffset? [m, axis-to-contact-surface distance] }`
 * @param {number} [model.nElem=100] total element count distributed by length where a section gives no nElem
 * @param {object} [model.endA] `{ x?, z? [m prescribed position; omit to leave free], angle? [rad, clamped tangent
 *   direction A→B; with kRot it is the spring reference angle], kRot? [N·m/rad flex-joint stiffness], Fx?, Fz? [N],
 *   M? [N·m] }`
 * @param {object} [model.endB] same fields as endA
 * @param {Array<object>} [model.pointLoads] `{ s [m], Fx?, Fz? [N, +z up: buoy positive, clump negative], M? [N·m] }`
 *   applied at the nearest node
 * @param {object} [model.seabed] `{ z [m elevation], kn [N/m per m = Pa, normal stiffness], mu? [-], kt? [Pa axial
 *   shear stiffness before sliding], profile? (x) => z [m] (trench / slope; overrides z) }`. Axial friction is an
 *   elastic–sliding (Coulomb) spring anchored where each node first touches down: it acts when mu > 0, with
 *   kt defaulting to kn when it is not given (full friction mobilised after a slip of mu × penetration);
 *   kt = 0 switches friction off explicitly. Friction develops only through movement along the seabed after
 *   first contact, so it depends on the load path (ramping) as in the real installation.
 * @param {Array<object>} [model.supports] rigid circles `{ xc, zc, R [m], k [Pa penalty stiffness per unit length],
 *   side? 'outside' | 'inside', angMin?, angMax? [rad polar angle window about the centre], edge? [m arc length inside
 *   each end of the window over which the stiffness tapers smoothly to zero; default one element length], offset? [m] }`
 * @param {object} [model.current] `{ rho? [kg/m³], U? [m/s uniform, +x] | profile? (z) => U | z: [], U: [] }`
 * @param {object} [model.pressure] `{ rhoInt [kg/m³], pIntRef [Pa at zIntRef], zIntRef [m], rhoExt [kg/m³],
 *   zSurface [m], nu? = 0.3, g? }`
 * @param {'straight'|'cable'|object} [model.initial='straight'] starting shape: stress-free straight line from end A
 *   at `initialAngle`; a lumped-catenary between the two prescribed end positions; or `{ x, z, angle? }` node
 *   arrays (e.g. a previous result on the same mesh)
 * @param {number} [model.initialAngle=0] rad, direction of the straight starting shape
 * @param {number} [model.nSteps] initial number of load/displacement steps (default 10 for 'straight', else 4)
 * @param {boolean} [model.rampWeight] ramp the weight with the step factor (default true for 'straight')
 * @param {boolean} [model.rampLoads=true] ramp point/end loads and drag with the step factor (set false when the
 *   starting shape is already close to equilibrium under the full loads)
 * @param {boolean} [model.rampSupports=true] ramp the support penalty stiffness with the step factor
 * @param {number} [model.tol=1e-8] force residual tolerance relative to the reference force
 * @param {number} [model.maxIter=60] Newton iterations per step
 * @returns {{ converged: boolean, iterations: number, steps: number, residual: number, outerIterations: number,
 *   s: Float64Array, x: Float64Array, z: Float64Array, angle: Float64Array, Te: Float64Array, Tw: Float64Array,
 *   shear: Float64Array, moment: Float64Array, curvature: Float64Array, reaction: Float64Array,
 *   friction: Float64Array, strain: Float64Array, endA: object, endB: object, state: object }}
 *   s, x, z [m]; angle [rad, nodal tangent]; Te, Tw (effective / true wall tension) and shear [N]; moment [N·m];
 *   curvature [1/m]; reaction (normal contact) and friction [N/m]; strain [-] axial; endA/endB
 *   `{ Fx, Fz [N], M [N·m] }` are the forces exerted by the support on the rod end; residual [N] is the final
 *   max-norm force imbalance; `state` is opaque and feeds rodModes / rodHarmonicResponse.
 */
export function solveRodStatic(model) {
  if (!model || typeof model !== 'object') throw new RangeError('rod: argument model must be an object');
  const mesh = buildMesh(model);
  const { ne, nn, ndof } = mesh;
  const A = endSpec(model.endA, 'endA'), Bn = endSpec(model.endB, 'endB');
  if (model.seabed) { num(model.seabed.kn, 'seabed.kn', { min: 0, strictMin: true });
    if (model.seabed.mu !== undefined) num(model.seabed.mu, 'seabed.mu', { min: 0 });
    if (model.seabed.kt !== undefined) num(model.seabed.kt, 'seabed.kt', { min: 0 }); if (typeof model.seabed.profile !== 'function') num(model.seabed.z, 'seabed.z'); }
  (model.supports ?? []).forEach((sp, k) => { num(sp.R, `supports[${k}].R`, { min: 0, strictMin: true }); num(sp.k, `supports[${k}].k`, { min: 0, strictMin: true }); num(sp.xc, `supports[${k}].xc`); num(sp.zc, `supports[${k}].zc`); });

  const init = model.initial ?? 'straight';
  let u = new Float64Array(ndof);
  let w0 = null;
  if (init === 'straight') {
    const a0 = model.initialAngle ?? 0;
    const xa = A.x ?? 0, za = A.z ?? 0;
    for (let k = 0; k < nn; k++) { u[3 * k] = xa + mesh.s[k] * Math.cos(a0); u[3 * k + 1] = za + mesh.s[k] * Math.sin(a0); u[3 * k + 2] = a0; }
  } else if (init === 'cable') {
    if (A.x === undefined || A.z === undefined || Bn.x === undefined || Bn.z === undefined) throw new RangeError("rod: argument initial = 'cable' needs both end positions prescribed");
    u = cableGuess(mesh, model, A.x, A.z, Bn.x, Bn.z);
  } else if (init && init.x && init.z) {
    if (init.x.length !== nn || init.z.length !== nn) throw new RangeError(`rod: argument initial.x/z must have ${nn} nodes`);
    for (let k = 0; k < nn; k++) { u[3 * k] = init.x[k]; u[3 * k + 1] = init.z[k]; }
    if (init.angle && init.angle.length === nn) for (let k = 0; k < nn; k++) u[3 * k + 2] = init.angle[k]; else fillAngles(mesh, u);
    if (init.state && init.state.mesh && init.state.mesh.ne === ne && model.rampWeight) w0 = Float64Array.from(init.state.mesh.w);
  } else throw new RangeError("rod: argument initial must be 'straight', 'cable' or { x, z }");

  const st = {
    Fx: new Float64Array(nn), Fz: new Float64Array(nn), Mn: new Float64Array(nn),
    anchor: new Float64Array(nn).fill(NaN), rotSprings: [], cons: [], w0,
  };
  const addLoad = (k, o) => { st.Fx[k] += o.Fx ?? 0; st.Fz[k] += o.Fz ?? 0; st.Mn[k] += o.M ?? 0; };
  for (const pl of model.pointLoads ?? []) { num(pl.s, 'pointLoads[].s', { min: 0, max: mesh.Ltot * (1 + 1e-9) }); addLoad(nearestNode(mesh, pl.s), pl); }
  [[A, 0], [Bn, ne]].forEach(([e, k]) => {
    addLoad(k, e);
    if (e.x !== undefined) st.cons.push({ dof: 3 * k, v0: u[3 * k], vt: e.x });
    if (e.z !== undefined) st.cons.push({ dof: 3 * k + 1, v0: u[3 * k + 1], vt: e.z });
    if (e.angle !== undefined) {
      if (e.kRot !== undefined) st.rotSprings.push({ dof: 3 * k + 2, k: e.kRot, th0: e.angle });
      else { const v0 = u[3 * k + 2]; st.cons.push({ dof: 3 * k + 2, v0, vt: v0 + wrapPi(e.angle - v0) }); }
    }
  });
  const isCon = new Uint8Array(ndof);
  for (const c of st.cons) isCon[c.dof] = 1;

  const rampWeight = model.rampWeight ?? (init === 'straight');
  const tol = model.tol ?? 1e-8, maxIter = model.maxIter ?? 60;
  let Fref = 1;
  for (let e = 0; e < ne; e++) Fref += Math.abs(mesh.w[e]) * mesh.L0[e];
  for (let k = 0; k < nn; k++) Fref = Math.max(Fref, Math.abs(st.Fx[k]), Math.abs(st.Fz[k]), Math.abs(st.Mn[k]) / mesh.hRef);

  const R = new Float64Array(ndof), du = new Float64Array(ndof);
  const B = bandCreate(ndof, 5);
  const sb = model.seabed;
  const setAnchors = (first) => {
    if (!sb || !(sb.mu > 0) || !(seabedKt(sb) > 0)) return;
    for (let k = 0; k < nn; k++) {
      const x = u[3 * k], pen = seabedZ(sb, x) + mesh.offset[k] - u[3 * k + 1];
      if (pen < 0) { st.anchor[k] = NaN; continue; }
      if (Number.isNaN(st.anchor[k]) || first) { st.anchor[k] = x; continue; }
      const lim = sb.mu * sb.kn * pen / seabedKt(sb);
      const dxa = x - st.anchor[k];
      if (Math.abs(dxa) > lim) st.anchor[k] = x - Math.sign(dxa) * lim;
    }
  };
  setAnchors(true);
  const hasLagged = !!(model.current || model.pressure);
  let totalIter = 0, lastRes = 0;
  const maxTrans = model.maxStep ?? 0.05 * mesh.Ltot, maxRot = 0.4;

  const newton = (lam, growMax = 12) => {
    const fac = { lam, weight: rampWeight && !w0 ? lam : 1, load: (model.rampLoads ?? true) ? lam : 1, support: (model.rampSupports ?? true) ? Math.max(lam, 1e-3) : 1 };
    let prev = Infinity, grow = 0, lastStep = Infinity;
    for (let it = 0; it < maxIter; it++) {
      fac.Nmax = 0;
      assemble(mesh, model, st, u, fac, R, B, null);
      const Fr = Fref + fac.Nmax;
      let rn = 0;
      for (let i = 0; i < ndof; i++) {
        if (isCon[i]) continue;
        const v = Math.abs(i % 3 === 2 ? R[i] / mesh.hRef : R[i]);
        if (v > rn) rn = v;
      }
      let cn = 0;
      for (let i = 0; i < ndof; i++) du[i] = -R[i];
      for (const c of st.cons) {
        const d = c.v0 + lam * (c.vt - c.v0) - u[c.dof];
        cn = Math.max(cn, Math.abs(d));
        applyConstraint(B, du, c.dof, d);
      }
      if (!Number.isFinite(rn)) return false;
      lastRes = rn;
      // second test: the Newton correction has reached round-off (very stiff EA limits the attainable residual)
      if (cn <= 1e-11 * (mesh.Ltot + 1) && (rn <= tol * Fr || (lastStep <= 1e-12 * mesh.Ltot && rn <= 1e-4 * Fr))) return true;
      if (rn > prev) { grow++; if (grow > growMax) return false; }
      prev = rn;
      bandFactor(B); bandSolve(B, du);
      let mt = 0, mr = 0;
      for (let k = 0; k < nn; k++) {
        mt = Math.max(mt, Math.abs(du[3 * k]), Math.abs(du[3 * k + 1]));
        mr = Math.max(mr, Math.abs(du[3 * k + 2]));
      }
      const sc = Math.min(1, maxTrans / (mt + 1e-300), maxRot / (mr + 1e-300));
      for (let i = 0; i < ndof; i++) u[i] += sc * du[i];
      lastStep = sc < 1 ? Infinity : Math.max(mt, mr * mesh.hRef);
      totalIter++;
    }
    return false;
  };

  // Globalised equilibrium iteration, used when the plain Newton iteration above does not settle (typically a
  // starting shape whose contact state is far from the final one, so that full Newton steps make nodes chatter in
  // and out of contact). Each correction solves (K + mu·diag|K|)·du = −R and is accepted only if it makes
  // progress; otherwise it is halved, and when no fraction of it helps the damping mu is raised
  // (Levenberg–Marquardt / pseudo-transient continuation; Nocedal & Wright, "Numerical Optimization", Ch. 3 and
  // §11.2). Progress is measured on the total potential energy (sufficient-decrease test along a descent
  // direction), which cannot stall at the slope changes of the contact forces and leads to a stable equilibrium.
  // Where the model has no potential (seabed friction, sloping seabed), or the energy change is below its rounding
  // noise, the weighted imbalance |R| is used instead. mu decays to zero as full steps are accepted, so the final
  // iterations are plain Newton and the convergence test is the same as above. Because every accepted step
  // improves a merit function, the outcome does not hinge on the rounding of any single iteration.
  const trial = new Float64Array(ndof), ancTrial = new Float64Array(nn), diag = new Float64Array(ndof), acc = [0];
  const conservative = !(sb && ((sb.mu > 0 && seabedKt(sb) > 0) || typeof sb.profile === 'function'));
  const robust = (lam) => {
    const fac = { lam, weight: rampWeight && !w0 ? lam : 1, load: (model.rampLoads ?? true) ? lam : 1, support: (model.rampSupports ?? true) ? Math.max(lam, 1e-3) : 1 };
    for (const c of st.cons) u[c.dof] = c.v0 + lam * (c.vt - c.v0);
    const norms = () => {
      let rn = 0, phi = 0;
      for (let i = 0; i < ndof; i++) {
        if (isCon[i]) continue;
        const v = Math.abs(i % 3 === 2 ? R[i] / mesh.hRef : R[i]);
        if (v > rn) rn = v;
        phi += v * v;
      }
      return [rn, phi];
    };
    let mu = 0, lastStep = Infinity;
    const budget = 12 * maxIter;
    for (let it = 0; it < budget; it++) {
      fac.Nmax = 0;
      assemble(mesh, model, st, u, fac, R, B, null);
      const [rn, phi] = norms();
      if (!Number.isFinite(rn)) return false;
      lastRes = rn;
      if (rn <= tol * (Fref + fac.Nmax) || (lastStep <= 1e-12 * mesh.Ltot && rn <= 1e-4 * (Fref + fac.Nmax))) return true;
      const E0 = conservative ? potentialEnergy(mesh, model, st, u, fac, acc) : 0, noise = 1e-13 * acc[0];
      for (let i = 0; i < ndof; i++) { du[i] = -R[i]; diag[i] = Math.abs(bandGet(B, i, i)); }
      if (mu > 0) for (let i = 0; i < ndof; i++) bandAdd(B, i, i, mu * diag[i]);
      for (const c of st.cons) applyConstraint(B, du, c.dof, 0);
      let solved = true;
      try { bandFactor(B); bandSolve(B, du); } catch (err) { solved = false; }
      let accepted = false, al = 0, full = 0;
      if (solved) {
        let mt = 0, mr = 0, slope = 0;
        for (let k = 0; k < nn; k++) {
          mt = Math.max(mt, Math.abs(du[3 * k]), Math.abs(du[3 * k + 1]));
          mr = Math.max(mr, Math.abs(du[3 * k + 2]));
        }
        for (let i = 0; i < ndof; i++) if (!isCon[i]) slope += R[i] * du[i];
        // with an indefinite tangent the correction can point uphill in energy: only more damping helps then
        if (Number.isFinite(mt) && Number.isFinite(mr) && !(conservative && slope >= 0)) {
          full = Math.min(1, maxTrans / (mt + 1e-300), maxRot / (mr + 1e-300));
          ancTrial.set(st.anchor);
          al = full;
          for (let ls = 0; ls < 7; ls++, al *= 0.5) {
            for (let i = 0; i < ndof; i++) trial[i] = u[i] + al * du[i];
            if (conservative && -al * slope > noise) {
              const E1 = potentialEnergy(mesh, model, st, trial, fac, null);
              if (Number.isFinite(E1) && E1 <= E0 + 1e-4 * al * slope) { accepted = true; break; }
            } else {
              assemble(mesh, model, st, trial, fac, R, null, null);
              const ph = norms()[1];
              if (Number.isFinite(ph) && ph <= (1 - 1e-4 * al) * phi) { accepted = true; break; }
              st.anchor.set(ancTrial);
            }
          }
          if (accepted) lastStep = al < 1 ? Infinity : Math.max(mt, mr * mesh.hRef);
        }
      }
      totalIter++;
      if (accepted) {
        u.set(trial);
        if (al === full) mu = mu < 1e-9 ? 0 : 0.2 * mu; else if (al < 0.2 * full) mu = Math.max(4 * mu, 1e-6);
      } else {
        mu = Math.max(10 * mu, 1e-6);
        if (mu > 1e12) return false;                         // no correction makes progress: no equilibrium nearby
      }
    }
    return false;
  };

  // Does the load factor change the problem materially? If weight, loads and supports are applied in full from
  // the start and the prescribed ends barely move, cutting the step poses the same problem again.
  let travel = 0, anyLoad = hasLagged;
  for (const c of st.cons) travel = Math.max(travel, Math.abs(c.vt - c.v0) * (c.dof % 3 === 2 ? mesh.hRef : 1));
  for (let k = 0; k < nn; k++) if (st.Fx[k] !== 0 || st.Fz[k] !== 0 || st.Mn[k] !== 0) anyLoad = true;
  const stepped = rampWeight || !!w0 || ((model.rampLoads ?? true) && anyLoad) || ((model.rampSupports ?? true) && (model.supports ?? []).length > 0)
    || travel > 1e-3 * mesh.Ltot;

  const saved = Float64Array.from(u), savedAnchor = Float64Array.from(st.anchor);
  let lam = 0, steps = 0, fails = 0;
  let dl = 1 / Math.max(1, Math.round(model.nSteps ?? (init === 'straight' ? 10 : 4)));
  while (lam < 1 - 1e-12) {
    const lt = Math.min(1, lam + dl);
    if (hasLagged) updateLaggedLoads(mesh, model, u);
    const before = totalIter;
    let ok;
    // Started close to equilibrium with nothing being stepped, a plain Newton iteration whose imbalance keeps
    // rising is wandering between contact states; where it would end up is then a matter of rounding, so it is
    // stopped early in favour of the globalised iteration.
    try { ok = newton(lt, stepped ? 12 : 3); } catch (err) { ok = false; }
    if (!ok) {
      // restart from the last converged state; switch to the globalised iteration at once when a smaller step
      // would pose the same problem again, otherwise after repeated failures at shrinking steps
      u.set(saved); st.anchor.set(savedAnchor);
      fails++;
      if (!stepped || fails % 3 === 0) {
        try { ok = robust(lt); } catch (err) { ok = false; }
        if (!ok) { u.set(saved); st.anchor.set(savedAnchor); }
        if (!ok && !stepped) dl = 0;
      }
    } else fails = 0;
    if (ok) {
      lam = lt; steps++;
      setAnchors(false);
      saved.set(u); savedAnchor.set(st.anchor);
      if (fails === 0 && totalIter - before < 8) dl *= 1.6;
      fails = 0;
    } else {
      dl *= 0.5;
      if (dl < 1e-7) {
        throw new Error(`rod: static solution did not converge (load factor ${lam.toFixed(4)}, residual ${lastRes.toExponential(2)} N, `
          + `${(lastRes / Fref).toExponential(1)} of the reference force): no stable equilibrium was found from the starting shape — `
          + 'check that the end conditions and tension can carry the weight, that the rod can reach its supports, and refine the mesh or start closer to the expected shape');
      }
    }
  }
  let outer = 0;
  if (hasLagged) {
    const maxOuter = model.maxOuter ?? 60;
    let change = Infinity;
    const prevU = new Float64Array(ndof);
    for (; outer < maxOuter && change > 1e-9 * mesh.Ltot; outer++) {
      prevU.set(u);
      updateLaggedLoads(mesh, model, u);
      let ok;
      try { ok = newton(1); } catch (err) { ok = false; }
      if (!ok) { u.set(prevU); try { ok = robust(1); } catch (err) { ok = false; } }
      if (!ok) throw new Error('rod: static solution did not converge in the drag / pressure iteration (the current or pressure load changes the shape too abruptly: check the current profile and drag coefficients)');
      change = 0;
      for (let k = 0; k < nn; k++) change = Math.max(change, Math.abs(u[3 * k] - prevU[3 * k]), Math.abs(u[3 * k + 1] - prevU[3 * k + 1]));
    }
    if (change > 1e-6 * mesh.Ltot) throw new Error('rod: drag / pressure fixed-point iteration did not converge');
  }
  st.w0 = null;
  return postProcess(mesh, model, st, u, { iterations: totalIter, steps, outer });
}

function postProcess(mesh, model, st, u, info) {
  const { ne, nn, ndof } = mesh;
  const R = new Float64Array(ndof);
  const out = { N: new Float64Array(ne), M1: new Float64Array(ne), M2: new Float64Array(ne), Ln: new Float64Array(ne), q: new Float64Array(nn), fr: new Float64Array(nn) };
  const fac = { lam: 1, weight: 1, load: 1, support: 1 };
  assemble(mesh, model, st, u, fac, R, null, out);
  const x = new Float64Array(nn), z = new Float64Array(nn), angle = new Float64Array(nn);
  const Te = new Float64Array(nn), Tw = new Float64Array(nn), shear = new Float64Array(nn);
  const moment = new Float64Array(nn), curvature = new Float64Array(nn), strain = new Float64Array(nn);
  const pr = model.pressure;
  let res = 0;
  const isCon = new Uint8Array(ndof);
  for (const c of st.cons) isCon[c.dof] = 1;
  for (let i = 0; i < ndof; i++) if (!isCon[i]) res = Math.max(res, Math.abs(i % 3 === 2 ? R[i] / mesh.hRef : R[i]));
  const Qe = (e) => (out.M1[e] + out.M2[e]) / out.Ln[e];
  for (let k = 0; k < nn; k++) {
    x[k] = u[3 * k]; z[k] = u[3 * k + 1]; angle[k] = u[3 * k + 2];
    const eL = k - 1, eR = k;
    let T, V, M, kap, eps;
    if (k === 0 || k === ne) {
      const e0 = k === 0 ? 0 : ne - 1, e1 = ne > 1 ? (k === 0 ? 1 : ne - 2) : e0;
      // element values are mid-element values: extrapolate linearly to the end node
      T = 1.5 * out.N[e0] - 0.5 * out.N[e1];
      V = 1.5 * Qe(e0) - 0.5 * Qe(e1);
      M = k === 0 ? -out.M1[0] : out.M2[ne - 1];
      kap = M / mesh.EI[e0];
      eps = (T + mesh.P0[e0]) / mesh.EA[e0];
    } else {
      const a = mesh.L0[eR] / (mesh.L0[eL] + mesh.L0[eR]);
      T = a * out.N[eL] + (1 - a) * out.N[eR];
      V = a * Qe(eL) + (1 - a) * Qe(eR);
      M = 0.5 * (out.M2[eL] - out.M1[eR]);
      kap = 0.5 * (out.M2[eL] / mesh.EI[eL] - out.M1[eR] / mesh.EI[eR]);
      eps = a * (out.N[eL] + mesh.P0[eL]) / mesh.EA[eL] + (1 - a) * (out.N[eR] + mesh.P0[eR]) / mesh.EA[eR];
    }
    Te[k] = T; shear[k] = V; moment[k] = M; curvature[k] = kap; strain[k] = eps;
    if (pr) {
      const g = pr.g ?? G0, e = Math.min(k, ne - 1);
      const pi = Math.max(0, (pr.pIntRef ?? 0) + (pr.rhoInt ?? 0) * g * ((pr.zIntRef ?? 0) - z[k]));
      const pe = Math.max(0, (pr.rhoExt ?? 1025) * g * ((pr.zSurface ?? 0) - z[k]));
      Tw[k] = T + pi * mesh.Ai[e] - pe * mesh.Ae[e];
    } else Tw[k] = T;
  }
  const endF = (k) => ({ Fx: R[3 * k], Fz: R[3 * k + 1], M: R[3 * k + 2] });
  return {
    converged: true, iterations: info.iterations, steps: info.steps, outerIterations: info.outer, residual: res,
    s: Float64Array.from(mesh.s), x, z, angle, Te, Tw, shear, moment, curvature, strain,
    reaction: Float64Array.from(out.q), friction: Float64Array.from(out.fr),
    endA: endF(0), endB: endF(ne),
    state: { mesh, model, st, u: Float64Array.from(u), elemN: out.N, elemLn: out.Ln },
  };
}

/* ------------------------------------------------------------------------------------------------------- */
/* Linearised dynamics about the static shape                                                               */

function assembleKM(state) {
  const { mesh, model, st, u } = state;
  const { ne, ndof } = mesh;
  const K = bandCreate(ndof, 5), M = bandCreate(ndof, 5);
  const R = new Float64Array(ndof);
  assemble(mesh, model, st, u, { lam: 1, weight: 1, load: 1, support: 1, symmetric: true }, R, K, null);
  for (let e = 0; e < ne; e++) {
    const i = 3 * e;
    const dx = u[i + 3] - u[i], dz = u[i + 4] - u[i + 1];
    const Ln = Math.hypot(dx, dz), c = dx / Ln, s = dz / Ln, L = mesh.L0[e];
    const mA = (mesh.m[e] + mesh.maT[e]) * L / 6, mT = (mesh.m[e] + mesh.ma[e]) * L / 420;
    const ml = [
      [2 * mA, 0, 0, mA, 0, 0],
      [0, 156 * mT, 22 * L * mT, 0, 54 * mT, -13 * L * mT],
      [0, 22 * L * mT, 4 * L * L * mT, 0, 13 * L * mT, -3 * L * L * mT],
      [mA, 0, 0, 2 * mA, 0, 0],
      [0, 54 * mT, 13 * L * mT, 0, 156 * mT, -22 * L * mT],
      [0, -13 * L * mT, -3 * L * L * mT, 0, -22 * L * mT, 4 * L * L * mT],
    ];
    // T: local = T global, block [[c, s, 0], [-s, c, 0], [0, 0, 1]]
    const T = [[c, s, 0], [-s, c, 0], [0, 0, 1]];
    for (let a = 0; a < 6; a++) {
      for (let b = 0; b < 6; b++) {
        let v = 0;
        const ab = a < 3 ? 0 : 3, bb = b < 3 ? 0 : 3;
        for (let p = 0; p < 3; p++) {
          const ta = T[p][a - ab];
          if (ta === 0) continue;
          for (let q = 0; q < 3; q++) v += ta * ml[ab + p][bb + q] * T[q][b - bb];
        }
        if (v !== 0) bandAdd(M, i + a, i + b, v);
      }
    }
  }
  return { K, M };
}

function bandMul(B, x, y) {
  const { n, hb, W, a } = B;
  for (let i = 0; i < n; i++) {
    let s = 0;
    const lo = Math.max(0, i - hb), hi = Math.min(n - 1, i + hb), base = i * W - i + hb;
    for (let j = lo; j <= hi; j++) s += a[base + j] * x[j];
    y[i] = s;
  }
  return y;
}

function modalCurvature(mesh, u, dU, kap) {
  const { ne } = mesh;
  kap.fill(0);
  for (let e = 0; e < ne; e++) {
    const i = 3 * e;
    const dx = u[i + 3] - u[i], dz = u[i + 4] - u[i + 1];
    const Ln = Math.hypot(dx, dz), c = dx / Ln, s = dz / Ln;
    const dbeta = (-s * (dU[i + 3] - dU[i]) + c * (dU[i + 4] - dU[i + 1])) / Ln;
    const t1 = dU[i + 2] - dbeta, t2 = dU[i + 5] - dbeta;
    const k1 = -(4 * t1 + 2 * t2) / mesh.L0[e], k2 = (2 * t1 + 4 * t2) / mesh.L0[e];
    const f1 = e === 0 ? 1 : 0.5, f2 = e === ne - 1 ? 1 : 0.5;
    kap[e] += f1 * k1; kap[e + 1] += f2 * k2;
  }
}

/**
 * Natural frequencies and mode shapes of small in-plane vibrations about a static rod solution.
 *
 * Method: tangent stiffness (elastic + geometric + contact springs) with a consistent mass matrix (structural
 * mass in all directions, added mass normal/tangential to the local axis), lowest modes by subspace iteration.
 * Reference: Bathe, "Finite Element Procedures", §11.6; added mass per DNV-RP-C205 §6.9.
 * Validity: linear vibration about a stable equilibrium with every rigid-body motion restrained; contact
 * status frozen; sections must carry `m` (and `ma` for submerged parts).
 *
 * @param {object} result output of solveRodStatic
 * @param {object} [opts]
 * @param {number} [opts.nModes=6] number of lowest modes [-]
 * @returns {{ s: Float64Array, frequenciesHz: number[], modes: Array<{ fHz: number, omega: number, periodS: number,
 *   phi: Float64Array, axial: Float64Array, dx: Float64Array, dz: Float64Array, curvature: Float64Array }>,
 *   converged: boolean }} s [m]; fHz [Hz]; omega [rad/s]; phi / axial are the displacement components normal /
 *   tangential to the static axis, scaled so max|phi| = 1 [-]; dx, dz global components [-]; curvature is the
 *   modal curvature per unit phi amplitude [1/m per m = 1/m²]
 */
export function rodModes(result, { nModes = 6 } = {}) {
  if (!result || !result.state) throw new RangeError('rod: argument result must come from solveRodStatic');
  const { mesh, st, u } = result.state;
  const { ndof, nn } = mesh;
  num(nModes, 'nModes', { min: 1, max: ndof });
  let mtot = 0;
  for (let e = 0; e < mesh.ne; e++) mtot += mesh.m[e] * mesh.L0[e];
  if (!(mtot > 0)) throw new RangeError('rod: argument sections[].m (mass per length) must be positive for dynamics');
  const { K, M } = assembleKM(result.state);
  const dummy = new Float64Array(ndof);
  for (const c of st.cons) {
    applyConstraint(K, dummy, c.dof, 0);
    applyConstraint(M, dummy, c.dof, 0);
  }
  let kmax = 0, mmax = 0;
  for (let i = 0; i < ndof; i++) { kmax = Math.max(kmax, Math.abs(bandGet(K, i, i))); mmax = Math.max(mmax, bandGet(M, i, i)); }
  // restrained freedoms are parked at a frequency far above the band of interest
  for (const c of st.cons) { bandSet(K, c.dof, c.dof, kmax); bandSet(M, c.dof, c.dof, 1e-9 * mmax); }
  bandFactor(K);
  const eig = subspaceIteration({
    n: ndof, nModes: Math.round(nModes),
    solveK: (b) => bandSolve(K, Float64Array.from(b)),
    mulM: (xv) => bandMul(M, xv, new Float64Array(ndof)),
    tol: 1e-9, maxIter: 80, seed: 7,
  });
  const modes = [];
  const kap = new Float64Array(nn);
  for (let k = 0; k < eig.values.length; k++) {
    const v = eig.vectors[k];
    const lam = eig.values[k];
    if (!(lam > 0)) throw new Error('rod: non-positive eigenvalue — the static configuration is not stable');
    const phi = new Float64Array(nn), axial = new Float64Array(nn), dx = new Float64Array(nn), dz = new Float64Array(nn);
    let big = 0;
    for (let i = 0; i < nn; i++) {
      const th = u[3 * i + 2];
      dx[i] = v[3 * i]; dz[i] = v[3 * i + 1];
      phi[i] = -Math.sin(th) * dx[i] + Math.cos(th) * dz[i];
      axial[i] = Math.cos(th) * dx[i] + Math.sin(th) * dz[i];
      if (Math.abs(phi[i]) > Math.abs(big)) big = phi[i];
    }
    const sc = big !== 0 ? 1 / big : 1;
    const dU = Float64Array.from(v, (q) => q * sc);
    modalCurvature(mesh, u, dU, kap);
    for (let i = 0; i < nn; i++) { phi[i] *= sc; axial[i] *= sc; dx[i] *= sc; dz[i] *= sc; }
    const omega = Math.sqrt(lam);
    modes.push({ fHz: omega / (2 * Math.PI), omega, periodS: 2 * Math.PI / omega, phi, axial, dx, dz, curvature: Float64Array.from(kap) });
  }
  return { s: Float64Array.from(mesh.s), frequenciesHz: modes.map((m) => m.fHz), modes, converged: eig.converged };
}

/**
 * Frequency-domain response of the rod to harmonic motion of one end about the static shape, with
 * harmonically linearised Morison drag damping (iterated on the velocity amplitude).
 *
 * Method: (K − ω²M + iωC) U = −Z·U_prescribed with the tangent stiffness and consistent mass of rodModes and
 * C from equivalent linearisation of |v|v for a sinusoid, c = (8/3π)·½ρ·Cd·D·|v̂| per unit length (normal) and the
 * same form on the wetted area for the tangential direction. Reference: DNV-RP-C205 §6.2; Faltinsen, "Sea Loads
 * on Ships and Offshore Structures", §7 (equivalent linearisation); Sparks (2007) Ch. 15.
 * Validity: small harmonic motion about the static shape, contact status frozen, steady current not included
 * in the drag linearisation, planar motion.
 *
 * LIMITATION AT TOUCHDOWN — READ BEFORE USING momentAmp FOR FATIGUE. With the contact status frozen, the node
 * where the static shape lands is a fixed support for the dynamic motion: the suspended span rotates about it and
 * all the dynamic bending of the touchdown zone is concentrated at that one node. In reality the touchdown point
 * travels along the pipe and each material point only sees the static curvature profile slide past it. The
 * frozen-contact momentAmp in the touchdown zone is therefore not the bending range a material point sees: on a
 * stiff seabed with a large touchdown excursion it has been found to over-predict the touchdown stress range of a
 * steel catenary riser about threefold, while for small excursions it can equally fall below the travelling
 * estimate. Do not use it for touchdown fatigue without the correction below. Away from the touchdown zone the
 * result is unaffected.
 * `touchdownSmoothing: { travel }` replaces the bending amplitude inside the touchdown zone by that kinematic
 * picture: for each node within 8λ + travel on the suspended side and 4λ + travel on the grounded side of a
 * static touchdown point (λ = √(EI/Te), the bending boundary-layer length), the curvature amplitude becomes half
 * the range of the static curvature |κ(s)| over s ± travel (nine sample points), and momentAmp = EI × that.
 * `travel` is the amplitude of the touchdown-point excursion along the rod [m] for the motion analysed: obtain it
 * from two static solutions with the end displaced by ± the motion amplitude (difference of the touchdown arc
 * length), or from the inextensible catenary relations (Irvine, "Cable Structures", §1.2). The frozen-contact
 * values are kept in momentAmpFrozen / curvatureAmpFrozen. Reference for the travelling-touchdown model:
 * Pesce, Aranha & Martins (1998), "The soil rigidity effect in the touchdown boundary-layer of a catenary riser:
 * static problem", ISOPE; Sparks (2007) Ch. 13 (boundary layer at touchdown).
 *
 * @param {object} result output of solveRodStatic (both translations of the excited end must be prescribed)
 * @param {object} opts
 * @param {number} opts.omega circular frequency [rad/s]
 * @param {'A'|'B'} [opts.end='A'] excited end
 * @param {number} [opts.xAmp=0] horizontal motion amplitude [m]
 * @param {number} [opts.zAmp=0] vertical motion amplitude [m]
 * @param {number} [opts.phaseZ=0] phase of z relative to x [rad]
 * @param {number} [opts.rotAmp=0] rotation amplitude of a clamped end [rad]
 * @param {number} [opts.rho=1025] fluid density for drag linearisation [kg/m³]
 * @param {number} [opts.zeta=0] structural damping ratio applied as stiffness-proportional damping at ω [-]
 * @param {number} [opts.tol=1e-4] relative convergence tolerance on the velocity amplitudes [-]
 * @param {number} [opts.maxIter=60]
 * @param {{travel: number}} [opts.touchdownSmoothing] post-process the touchdown zone with the travelling
 *   touchdown model (see above); travel = touchdown excursion amplitude along the rod [m], >= 0
 * @returns {{ s: Float64Array, xAmp: Float64Array, zAmp: Float64Array, normalAmp: Float64Array,
 *   tensionAmp: Float64Array, momentAmp: Float64Array, curvatureAmp: Float64Array, iterations: number,
 *   converged: boolean, momentAmpFrozen?: Float64Array, curvatureAmpFrozen?: Float64Array,
 *   touchdown?: {s: number, node: number, lambda: number, from: number, to: number}[] }} displacement amplitudes
 *   [m], dynamic effective tension amplitude [N], bending moment amplitude [N·m], curvature amplitude [1/m], all at
 *   nodes; with touchdownSmoothing also the frozen-contact arrays and, per static touchdown point, its arc
 *   length [m], node, boundary-layer length λ [m] and the smoothed interval of arc length [m]
 */
export function rodHarmonicResponse(result, opts) {
  if (!result || !result.state) throw new RangeError('rod: argument result must come from solveRodStatic');
  const { omega, end = 'A', xAmp = 0, zAmp = 0, phaseZ = 0, rotAmp = 0, rho = 1025, zeta = 0, tol = 1e-4, maxIter = 60, touchdownSmoothing } = opts ?? {};
  num(omega, 'omega', { min: 0, strictMin: true });
  if (touchdownSmoothing !== undefined) {
    if (!touchdownSmoothing || typeof touchdownSmoothing !== 'object') throw new RangeError('rod: argument touchdownSmoothing must be { travel }');
    num(touchdownSmoothing.travel, 'touchdownSmoothing.travel', { min: 0 });
    if (!result.reaction || !result.curvature) throw new RangeError('rod: argument result must carry reaction and curvature for touchdownSmoothing');
  }
  num(zeta, 'zeta', { min: 0 });
  const { mesh, st, u } = result.state;
  const { ndof, nn, ne } = mesh;
  const { K, M } = assembleKM(result.state);
  const kNode = end === 'B' ? ne : 0;
  const isCon = new Uint8Array(ndof);
  for (const c of st.cons) isCon[c.dof] = 1;
  if (!isCon[3 * kNode] || !isCon[3 * kNode + 1]) throw new RangeError('rod: the excited end must have prescribed x and z');
  const Ure = new Float64Array(ndof), Uim = new Float64Array(ndof);
  const pre = new Float64Array(2 * ndof);
  pre[2 * (3 * kNode)] = xAmp;
  pre[2 * (3 * kNode + 1)] = zAmp * Math.cos(phaseZ); pre[2 * (3 * kNode + 1) + 1] = zAmp * Math.sin(phaseZ);
  if (isCon[3 * kNode + 2]) pre[2 * (3 * kNode + 2)] = rotAmp;
  const cn = new Float64Array(nn), ct = new Float64Array(nn), dn = new Float64Array(nn), dt = new Float64Array(nn);
  for (let e = 0; e < ne; e++) {
    const a = 0.5 * (8 / (3 * Math.PI)) * 0.5 * rho * mesh.D[e] * mesh.L0[e];
    dn[e] += a * mesh.Cdn[e]; dn[e + 1] += a * mesh.Cdn[e];
    dt[e] += a * mesh.Cdt[e] * Math.PI; dt[e + 1] += a * mesh.Cdt[e] * Math.PI;
  }
  const v0 = omega * Math.hypot(xAmp, zAmp) * 0.5;
  const vn = new Float64Array(nn).fill(v0), vt = new Float64Array(nn).fill(v0);
  const Z = bandCreate(2 * ndof, 11);
  const rhs = new Float64Array(2 * ndof);
  let it = 0, converged = false;
  for (; it < maxIter && !converged; it++) {
    Z.a.fill(0); rhs.fill(0);
    for (let k = 0; k < nn; k++) { cn[k] = dn[k] * vn[k]; ct[k] = dt[k] * vt[k]; }
    for (let i = 0; i < ndof; i++) {
      const lo = Math.max(0, i - 5), hi = Math.min(ndof - 1, i + 5);
      for (let j = lo; j <= hi; j++) {
        const kij = bandGet(K, i, j);
        const re = kij - omega * omega * bandGet(M, i, j);
        let im = 2 * zeta * kij;
        const ni = (i / 3) | 0;
        if (ni === ((j / 3) | 0) && i % 3 < 2 && j % 3 < 2) {
          const sn = Math.sin(u[3 * ni + 2]), cs = Math.cos(u[3 * ni + 2]);
          const ni_ = i % 3 === 0 ? -sn : cs, nj_ = j % 3 === 0 ? -sn : cs;
          const ti_ = i % 3 === 0 ? cs : sn, tj_ = j % 3 === 0 ? cs : sn;
          im += omega * (cn[ni] * ni_ * nj_ + ct[ni] * ti_ * tj_);
        }
        if (re !== 0) { bandAdd(Z, 2 * i, 2 * j, re); bandAdd(Z, 2 * i + 1, 2 * j + 1, re); }
        if (im !== 0) { bandAdd(Z, 2 * i, 2 * j + 1, -im); bandAdd(Z, 2 * i + 1, 2 * j, im); }
      }
    }
    for (const c of st.cons) {
      applyConstraint(Z, rhs, 2 * c.dof, pre[2 * c.dof]);
      applyConstraint(Z, rhs, 2 * c.dof + 1, pre[2 * c.dof + 1]);
    }
    bandFactor(Z); bandSolve(Z, rhs);
    for (let i = 0; i < ndof; i++) { Ure[i] = rhs[2 * i]; Uim[i] = rhs[2 * i + 1]; }
    let err = 0, ref = 1e-12;
    for (let k = 0; k < nn; k++) {
      const th = u[3 * k + 2], sn = Math.sin(th), cs = Math.cos(th);
      const a = omega * Math.hypot(-sn * Ure[3 * k] + cs * Ure[3 * k + 1], -sn * Uim[3 * k] + cs * Uim[3 * k + 1]);
      const b = omega * Math.hypot(cs * Ure[3 * k] + sn * Ure[3 * k + 1], cs * Uim[3 * k] + sn * Uim[3 * k + 1]);
      if (dn[k] > 0) { err = Math.max(err, Math.abs(a - vn[k])); ref = Math.max(ref, a); }
      if (dt[k] > 0) { err = Math.max(err, Math.abs(b - vt[k])); ref = Math.max(ref, b); }
      vn[k] = 0.5 * (vn[k] + a); vt[k] = 0.5 * (vt[k] + b);
    }
    if (err <= tol * ref) converged = true;
  }
  const xA = new Float64Array(nn), zA = new Float64Array(nn), nA = new Float64Array(nn);
  for (let k = 0; k < nn; k++) {
    const th = u[3 * k + 2], sn = Math.sin(th), cs = Math.cos(th);
    xA[k] = Math.hypot(Ure[3 * k], Uim[3 * k]); zA[k] = Math.hypot(Ure[3 * k + 1], Uim[3 * k + 1]);
    nA[k] = Math.hypot(-sn * Ure[3 * k] + cs * Ure[3 * k + 1], -sn * Uim[3 * k] + cs * Uim[3 * k + 1]);
  }
  const kr = new Float64Array(nn), ki = new Float64Array(nn);
  modalCurvature(mesh, u, Ure, kr); modalCurvature(mesh, u, Uim, ki);
  const curvatureAmp = new Float64Array(nn), momentAmp = new Float64Array(nn), tensionAmp = new Float64Array(nn);
  const dN = new Float64Array(ne);
  for (let e = 0; e < ne; e++) {
    const i = 3 * e;
    const dx = u[i + 3] - u[i], dz = u[i + 4] - u[i + 1];
    const Ln = Math.hypot(dx, dz), c = dx / Ln, s = dz / Ln, ka = mesh.EA[e] / mesh.L0[e];
    const er = c * (Ure[i + 3] - Ure[i]) + s * (Ure[i + 4] - Ure[i + 1]);
    const ei = c * (Uim[i + 3] - Uim[i]) + s * (Uim[i + 4] - Uim[i + 1]);
    dN[e] = ka * Math.hypot(er, ei);
  }
  for (let k = 0; k < nn; k++) {
    curvatureAmp[k] = Math.hypot(kr[k], ki[k]);
    momentAmp[k] = curvatureAmp[k] * mesh.EI[Math.min(k, ne - 1)];
    tensionAmp[k] = k === 0 ? dN[0] : k === ne ? dN[ne - 1] : 0.5 * (dN[k - 1] + dN[k]);
  }
  const out = { s: Float64Array.from(mesh.s), xAmp: xA, zAmp: zA, normalAmp: nA, tensionAmp, momentAmp, curvatureAmp, iterations: it, converged };
  if (touchdownSmoothing) {
    out.momentAmpFrozen = Float64Array.from(momentAmp);
    out.curvatureAmpFrozen = Float64Array.from(curvatureAmp);
    out.touchdown = smoothTouchdown(mesh, result, touchdownSmoothing.travel, curvatureAmp, momentAmp);
  }
  return out;
}

// Travelling-touchdown post-process of rodHarmonicResponse: inside each touchdown zone the curvature amplitude is
// half the range of the static |curvature| swept over ± travel. Returns the touchdown points treated.
function smoothTouchdown(mesh, result, travel, curvatureAmp, momentAmp) {
  const { nn, ne, s } = mesh;
  const kap = new Float64Array(nn);
  for (let k = 0; k < nn; k++) kap[k] = Math.abs(result.curvature[k]);
  const at = (x) => {
    if (x <= s[0]) return kap[0];
    if (x >= s[nn - 1]) return kap[nn - 1];
    let lo = 0, hi = nn - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (s[m] <= x) lo = m; else hi = m; }
    return kap[lo] + (kap[hi] - kap[lo]) * (x - s[lo]) / (s[hi] - s[lo]);
  };
  const touching = (k) => result.reaction[k] > 0;
  const points = [];
  for (let k = 0; k < nn - 1; k++) {
    if (touching(k) === touching(k + 1)) continue;
    // a supported rod end resting on the seabed carries its load through the support: not a touchdown point
    if (k < 2 || k > nn - 4) continue;
    // node kT is the first grounded node; dir points from the grounded side to the suspended side
    const kT = touching(k) ? k : k + 1, dir = touching(k) ? 1 : -1;
    const lambda = Math.sqrt(mesh.EI[Math.min(kT, ne - 1)] / Math.max(Math.abs(result.Te[kT]), 1));
    const a = s[kT] + dir * (8 * lambda + travel), b = s[kT] - dir * (4 * lambda + travel);
    points.push({ s: s[kT], node: kT, lambda, from: Math.max(Math.min(a, b), s[0]), to: Math.min(Math.max(a, b), s[nn - 1]) });
  }
  const done = new Uint8Array(nn);
  for (const p of points) {
    for (let k = 0; k < nn; k++) {
      if (done[k] || s[k] < p.from || s[k] > p.to) continue;
      let mn = Infinity, mx = -Infinity;
      for (let q = -4; q <= 4; q++) { const v = at(s[k] + travel * q / 4); if (v < mn) mn = v; if (v > mx) mx = v; }
      curvatureAmp[k] = 0.5 * (mx - mn);
      momentAmp[k] = curvatureAmp[k] * mesh.EI[Math.min(k, ne - 1)];
      done[k] = 1;
    }
  }
  return points;
}

/**
 * Lumped-chain catenary through two pinned points (no bending stiffness), exposed for cross-checks of the rod
 * solver and as a fast configuration estimate.
 *
 * Method: bar chain with nodal weights, shooting on the end force (H, V); with a seabed at the elevation of one
 * end the chain may rest on it before lifting off. Reference: Irvine, "Cable Structures", §1.2–1.3.
 * Validity: no current, frictionless flat seabed, both ends pinned.
 *
 * @param {object} model same `sections`, `nElem`, `endA`, `endB`, `seabed` as solveRodStatic
 * @returns {{ s: Float64Array, x: Float64Array, z: Float64Array, angle: Float64Array }} node arrays [m, m, m, rad]
 */
export function cableShape(model) {
  const mesh = buildMesh(model);
  const A = endSpec(model.endA, 'endA'), Bn = endSpec(model.endB, 'endB');
  for (const [e, n] of [[A, 'endA'], [Bn, 'endB']]) { num(e.x, `${n}.x`); num(e.z, `${n}.z`); }
  const u = cableGuess(mesh, model, A.x, A.z, Bn.x, Bn.z);
  const nn = mesh.nn, x = new Float64Array(nn), z = new Float64Array(nn), angle = new Float64Array(nn);
  for (let k = 0; k < nn; k++) { x[k] = u[3 * k]; z[k] = u[3 * k + 1]; angle[k] = u[3 * k + 2]; }
  return { s: Float64Array.from(mesh.s), x, z, angle };
}
