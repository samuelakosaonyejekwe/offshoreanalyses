/**
 * Single uniform mooring line in a vertical plane: inextensible closed forms and the elastic catenary with a
 * frictional (horizontal or sloping) seabed, as used in quasi-static mooring codes.
 *
 * Frame: anchor at the origin, fairlead at (X, Z), X >= 0 horizontal, Z positive upward. Unstretched arc length s
 * runs from the anchor (0) to the fairlead (L). H and V are the horizontal and vertical components of the
 * line tension at the fairlead (V > 0 when the line pulls the fairlead downward).
 *
 * References:
 *  - Irvine, H.M. (1981) Cable Structures, MIT Press, ch. 1 (elastic catenary).
 *  - Jonkman, J.M. (2007) Dynamics Modeling and Loads Analysis of an Offshore Floating Wind Turbine,
 *    NREL/TP-500-41958, section 2.5 (Eqs. 2-35 to 2-41: suspended line, line resting on a seabed with friction
 *    coefficient CB, initial guess after Peyrot & Goulois 1979).
 *  - Faltinsen, O.M. (1990) Sea Loads on Ships and Offshore Structures, ch. 8 (inextensible catenary, Eq. 8.26).
 */

import { brent } from '../math/roots.js';

const asinh = Math.asinh;
const hypot = Math.hypot;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}

function finite(v, name) {
  if (typeof v !== 'number' || !Number.isFinite(v)) throw new RangeError(`${name} must be a finite number`);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Closed-form inextensible catenary                                                                            */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Inextensible catenary hanging from a horizontal-tangent point (touchdown or sag bend) up to a point z above it.
 * T = H + w z, s = sqrt(z (z + 2H/w)), x = (H/w) acosh(1 + w z / H).
 * Reference: Faltinsen (1990) Eqs. 8.14-8.20. Valid for w > 0, H > 0, z >= 0, no elasticity.
 * @param {object} o
 * @param {number} o.H horizontal tension [N]
 * @param {number} o.w submerged weight per unit length [N/m]
 * @param {number} o.z height above the horizontal-tangent point [m]
 * @returns {{s:number,x:number,T:number,V:number,angle:number,a:number}} suspended length s [m], horizontal
 *   projection x [m], tension T [N], vertical tension component V [N], angle from horizontal [rad],
 *   catenary parameter a = H/w [m]
 */
export function inextensibleCatenary({ H, w, z }) {
  finite(H, 'H'); finite(w, 'w'); finite(z, 'z');
  need(H > 0, 'H must be > 0');
  need(w > 0, 'w must be > 0');
  need(z >= 0, 'z must be >= 0');
  const a = H / w;
  const s = Math.sqrt(z * (z + 2 * a));
  const x = a * Math.acosh(1 + z / a);
  const V = w * s;
  return { s, x, T: H + w * z, V, angle: Math.atan2(V, H), a };
}

/**
 * Shape of the inextensible catenary measured from its horizontal-tangent point:
 * z(x) = a (cosh(x/a) - 1), s(x) = a sinh(x/a), T(x) = H cosh(x/a), a = H/w.
 * Reference: Faltinsen (1990) Eq. 8.16-8.18. Valid for w > 0, H > 0.
 * @param {object} o
 * @param {number} o.H horizontal tension [N]
 * @param {number} o.w submerged weight per unit length [N/m]
 * @param {number} o.x horizontal distance from the horizontal-tangent point [m]
 * @returns {{z:number,s:number,T:number,angle:number}} elevation [m], arc length [m], tension [N], angle [rad]
 */
export function inextensibleCatenaryShape({ H, w, x }) {
  finite(H, 'H'); finite(w, 'w'); finite(x, 'x');
  need(H > 0, 'H must be > 0');
  need(w > 0, 'w must be > 0');
  const a = H / w;
  return { z: a * (Math.cosh(x / a) - 1), s: a * Math.sinh(x / a), T: H * Math.cosh(x / a), angle: Math.atan(Math.sinh(x / a)) };
}

/**
 * Horizontal stiffness dH/dX of an inextensible catenary line with part of its (constant) length resting on a
 * frictionless horizontal seabed: k = w / (acosh(1 + w h / H) - 2 / sqrt(1 + 2 H / (w h))).
 * Reference: Faltinsen (1990) Eq. 8.26. Valid for w > 0, H > 0, h > 0, anchor never lifted.
 * @param {object} o
 * @param {number} o.H horizontal tension [N]
 * @param {number} o.w submerged weight per unit length [N/m]
 * @param {number} o.z height of the fairlead above the seabed [m]
 * @returns {number} horizontal stiffness [N/m]
 */
export function inextensibleHorizontalStiffness({ H, w, z }) {
  finite(H, 'H'); finite(w, 'w'); finite(z, 'z');
  need(H > 0, 'H must be > 0');
  need(w > 0, 'w must be > 0');
  need(z > 0, 'z must be > 0');
  const a = H / w;
  return w / (Math.acosh(1 + z / a) - 2 / Math.sqrt(1 + 2 * a / z));
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Shooting solver shared with the multi-segment line                                                           */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Starting values for (H, V) of a uniform line. Reference: Peyrot & Goulois (1979) as quoted by Jonkman (2007)
 * Eq. 2-39/2-40; the taut branch is replaced by the straight elastic-bar estimate when EA is given.
 * @param {object} o
 * @param {number} o.X horizontal span [m]
 * @param {number} o.Z vertical span [m]
 * @param {number} o.L unstretched length [m]
 * @param {number} o.w weight per unit length [N/m] (either sign, non-zero)
 * @param {number} [o.EA] axial stiffness [N]
 * @returns {{H:number,V:number}} starting tension components at the fairlead [N]
 */
export function catenaryInitialGuess({ X, Z, L, w, EA }) {
  const chord = hypot(X, Z);
  let lam;
  if (X <= 1e-12 * L) lam = 1e6;
  else if (L <= chord) lam = 0.2;
  else lam = Math.sqrt(3 * ((L * L - Z * Z) / (X * X) - 1));
  lam = Math.max(lam, 1e-3);
  let H = Math.abs(w * X / (2 * lam));
  const zt = Z / Math.tanh(lam);
  let V = w > 0 ? 0.5 * w * (zt + L) : 0.5 * w * (L - zt);
  if (EA && L < chord) {
    const T = EA * (chord - L) / L;
    H = Math.max(H, T * X / chord);
    V = T * Z / chord + 0.5 * w * L;
  }
  if (!(H > 0)) H = Math.abs(w) * L * 1e-3;
  return { H, V };
}

/**
 * Solve X(H, V) = X, Z(H, V) = Z for the fairlead tension components of a line whose end position is given by a
 * forward (shooting) function. Damped Newton in (ln H, V) with a central-difference Jacobian. A step can carry
 * the iterate into the region where the whole line lies on the seabed (V below the lift-off value): there the
 * end position does not depend on V, the Jacobian is singular and Newton has no direction. The iterate is then
 * put back on the curve Z(H, V) = Z at its current H by a one-dimensional bisection in V (Z is monotonic in V)
 * and Newton continues from that point, where the line reaches the fairlead and the Jacobian is regular. Only if
 * that also fails is the nested bisection used (X is monotonic in H along Z = const), followed by
 * a Newton polish. A slack line (X not larger than the span reachable with zero horizontal tension) is detected
 * and returned with H = Hmin and slack = true.
 * @param {(H:number, V:number, out:{X:number,Z:number}) => void} fwd forward function filling out.X, out.Z [m]
 * @param {object} o
 * @param {number} o.X target horizontal span [m]
 * @param {number} o.Z target vertical span [m]
 * @param {number} o.H0 starting H [N] (> 0)
 * @param {number} o.V0 starting V [N]
 * @param {number} o.scaleL length scale for the residual [m]
 * @param {number} o.scaleF force scale (total line weight) [N]
 * @param {number} o.Hmax upper bound on H for the bisection fallback [N]
 * @param {boolean} [o.checkSlack=false] test for the slack-line state before iterating
 * @param {number} [o.tol=1e-11] residual tolerance relative to scaleL [-]
 * @returns {{H:number,V:number,slack:boolean,iterations:number,residual:number}} H, V [N]; residual [m]
 */
export function solveShooting(fwd, o) {
  const { X, Z, scaleL, scaleF, Hmax } = o;
  const tol = (o.tol ?? 1e-11) * scaleL;
  const out = { X: 0, Z: 0 };
  const Hmin = 1e-9 * scaleF;
  let evals = 0;

  const innerV = (H, Vg) => {
    // Z(H, V) is non-decreasing in V: bracket then bisect.
    let step = scaleF;
    let lo = Vg, hi = Vg;
    fwd(H, lo, out); evals++;
    let k = 0;
    while (out.Z > Z && k++ < 80) { lo -= step; step *= 2; fwd(H, lo, out); evals++; }
    step = scaleF; k = 0;
    fwd(H, hi, out); evals++;
    while (out.Z < Z && k++ < 80) { hi += step; step *= 2; fwd(H, hi, out); evals++; }
    for (let it = 0; it < 200; it++) {
      const mid = 0.5 * (lo + hi);
      if (mid === lo || mid === hi) break;
      fwd(H, mid, out); evals++;
      if (out.Z < Z) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
  };

  if (o.checkSlack) {
    const Vs = innerV(Hmin, Math.max(o.V0, 0));
    fwd(Hmin, Vs, out);
    if (out.X >= X - tol && Math.abs(out.Z - Z) < 1e-6 * scaleL) {
      return { H: Hmin, V: Vs, slack: true, iterations: evals, residual: Math.abs(out.Z - Z) };
    }
  }

  const newton = (H0, V0, maxIter) => {
    let lh = Math.log(Math.min(Math.max(H0, Hmin), Hmax));
    let V = V0;
    fwd(Math.exp(lh), V, out); evals++;
    let r1 = out.X - X, r2 = out.Z - Z;
    let f = hypot(r1, r2);
    let restores = 0;
    for (let it = 0; it < maxIter && f > tol; it++) {
      const eh = 1e-6;
      const dv = 1e-6 * Math.max(Math.abs(V), 1e-3 * scaleF);
      fwd(Math.exp(lh + eh), V, out); const xp = out.X, zp = out.Z;
      fwd(Math.exp(lh - eh), V, out);
      const a11 = (xp - out.X) / (2 * eh), a21 = (zp - out.Z) / (2 * eh);
      fwd(Math.exp(lh), V + dv, out); const xq = out.X, zq = out.Z;
      fwd(Math.exp(lh), V - dv, out);
      const a12 = (xq - out.X) / (2 * dv), a22 = (zq - out.Z) / (2 * dv);
      evals += 4;
      const det = a11 * a22 - a12 * a21;
      // singular (or, relative to its entries, numerically singular) Jacobian: restore Z at this H and go on
      const degenerate = !(Math.abs(det) > 1e-12 * (Math.abs(a11 * a22) + Math.abs(a12 * a21))) || !Number.isFinite(det);
      if (degenerate) {
        if (restores++ >= 4) break;
        V = innerV(Math.exp(lh), Math.max(V, 0));
        fwd(Math.exp(lh), V, out); evals++;
        r1 = out.X - X; r2 = out.Z - Z; f = hypot(r1, r2);
        continue;
      }
      let dl = -(a22 * r1 - a12 * r2) / det;
      let dV = -(-a21 * r1 + a11 * r2) / det;
      if (Math.abs(dl) > 2.5) { const c = 2.5 / Math.abs(dl); dl *= c; dV *= c; }
      let alpha = 1, ok = false;
      for (let ls = 0; ls < 16; ls++) {
        const lhn = lh + alpha * dl, Vn = V + alpha * dV;
        fwd(Math.exp(lhn), Vn, out); evals++;
        const fn = hypot(out.X - X, out.Z - Z);
        if (fn < f) { lh = lhn; V = Vn; r1 = out.X - X; r2 = out.Z - Z; f = fn; ok = true; break; }
        alpha *= 0.5;
      }
      if (!ok) break;
    }
    return { H: Math.exp(lh), V, f };
  };

  let s = newton(o.H0, o.V0, 60);
  if (!(s.f <= tol * 1e3)) {
    // Nested bisection: g(ln H) = X(H, V(H)) - X is increasing.
    let Vg = o.V0;
    const g = (lh) => { const H = Math.exp(lh); Vg = innerV(H, Vg); fwd(H, Vg, out); return out.X - X; };
    let lo = Math.log(Hmin), hi = Math.log(Hmax);
    const glo = g(lo);
    if (glo >= 0) {
      const V = Vg;
      fwd(Hmin, V, out);
      return { H: Hmin, V, slack: true, iterations: evals, residual: Math.abs(out.Z - Z) };
    }
    if (g(hi) < 0) throw new Error('Line solve failed: required span exceeds what the line can reach within its elastic range');
    for (let it = 0; it < 60; it++) {
      const mid = 0.5 * (lo + hi);
      if (g(mid) < 0) lo = mid; else hi = mid;
    }
    const Hb = Math.exp(0.5 * (lo + hi));
    const Vb = innerV(Hb, Vg);
    s = newton(Hb, Vb, 30);
    if (!(s.f <= 1e-6 * scaleL)) throw new Error('Line solve did not converge');
  }
  return { H: s.H, V: s.V, slack: false, iterations: evals, residual: s.f };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Elastic catenary with seabed                                                                                 */
/* ------------------------------------------------------------------------------------------------------------ */

// Forward evaluation for w > 0. State st receives X, Z and the internal partition of the line.
function forward(H, V, p, st) {
  const { L, w, EA, seabed, tanB, cosB, sinB, q } = p;
  let Ls = L, V0 = V - w * L, LB = 0, G = 0;
  if (seabed && V0 < H * tanB) {
    Ls = (V - H * tanB) / w;
    if (Ls < 0) Ls = 0;
    LB = L - Ls;
    V0 = V - w * Ls;
    const T0 = H / cosB;
    G = (q > 0 && T0 < q * LB) ? T0 * T0 / (2 * q) : T0 * LB - 0.5 * q * LB * LB;
  }
  const V1 = V0 + w * Ls;
  const Lbs = LB + G / EA;
  const xs = (H / w) * (asinh(V1 / H) - asinh(V0 / H)) + H * Ls / EA;
  const zs = (hypot(H, V1) - hypot(H, V0)) / w + (V0 * Ls + 0.5 * w * Ls * Ls) / EA;
  st.X = Lbs * cosB + xs;
  st.Z = Lbs * sinB + zs;
  st.Ls = Ls; st.LB = LB; st.V0 = V0; st.Lbs = Lbs;
}

function validate(o) {
  for (const k of ['X', 'Z', 'L', 'w', 'EA']) finite(o[k], k);
  need(o.X >= 0, 'X must be >= 0');
  need(o.L > 0, 'L must be > 0');
  need(o.EA > 0, 'EA must be > 0');
}

/**
 * Elastic catenary between an anchor at (0, 0) and a fairlead at (X, Z): solves for the fairlead tension
 * components (H, V) by damped Newton with a robust starting guess.
 *
 * Cases handled: fully suspended line; line partly resting on the seabed with friction coefficient CB (the anchor
 * tension is reduced by the friction on the grounded length, down to zero); horizontal or uniformly sloping
 * seabed through the anchor (slope positive when the seabed rises towards the fairlead; friction is assumed
 * fully mobilised against sliding towards the fairlead); taut and vertical lines; slack lines (X smaller than
 * the span covered by the grounded length with the remainder hanging vertically: H = 0, slack = true); buoyant
 * lines (w < 0, always fully suspended, solved by mirroring); weightless lines (straight elastic bar).
 *
 * Reference: Jonkman (2007) NREL/TP-500-41958 Eqs. 2-35 to 2-38 (horizontal seabed). The sloping-seabed form is
 * the same derivation with the touchdown tangent at the slope angle and the along-slope tension gradient
 * w (sin b + CB cos b) on the grounded part. Validity: linear-elastic line, small strain, no bending stiffness,
 * static, anchor on the seabed when seabed = true; a buoyant line reaching the free surface is not modelled.
 *
 * @param {object} o
 * @param {number} o.X horizontal anchor-fairlead distance [m] (>= 0)
 * @param {number} o.Z fairlead elevation above the anchor [m]
 * @param {number} o.L unstretched line length [m]
 * @param {number} o.w submerged weight per unit length [N/m] (negative for a buoyant line)
 * @param {number} o.EA axial stiffness [N]
 * @param {number} [o.CB=0] seabed friction coefficient [-]
 * @param {boolean} [o.seabed=true] seabed contact considered (anchor on the seabed)
 * @param {number} [o.slope=0] seabed slope angle [rad], positive rising from anchor towards fairlead
 * @param {number} [o.nProfile=0] number of profile points to return (0 = none)
 * @param {{H:number,V:number}} [o.guess] warm-start values [N]
 * @returns {object} H, V, T [N] and angle [rad from horizontal] at the fairlead; anchorH, anchorV, anchorT [N],
 *   anchorAngle [rad from horizontal] and upliftAngle [rad from the seabed]; groundedLength, suspendedLength
 *   (unstretched) [m]; touchdownX, touchdownZ [m]; stretchedLength [m]; stiffness = [[dH/dX, dH/dZ],
 *   [dV/dX, dV/dZ]] [N/m]; slack, converged flags; evaluations = forward (shooting) evaluations used by the
 *   solve [-]; profile {s, x, z, T} when nProfile > 0.
 */
export function solveCatenary(o) {
  validate(o);
  const { X, Z, L, EA } = o;
  const w = o.w;
  const CB = o.CB ?? 0;
  const slope = o.slope ?? 0;
  need(CB >= 0, 'CB must be >= 0');
  need(Math.abs(slope) < Math.PI / 2 - 1e-6, 'slope must be within (-pi/2, pi/2)');
  const chord = hypot(X, Z);

  // Weightless line: straight linear-elastic bar, zero tension when slack.
  if (Math.abs(w) * L <= 1e-12 * EA) return weightless(o, chord);

  const mirror = w < 0;
  const seabed = mirror ? false : (o.seabed ?? true);
  const wz = Math.abs(w);
  const Zs = mirror ? -Z : Z;
  if (seabed) need(Z >= X * Math.tan(slope) - 1e-9 * L, 'Z places the fairlead below the seabed');
  const p = {
    L, w: wz, EA, CB, seabed,
    tanB: Math.tan(slope), cosB: Math.cos(slope), sinB: Math.sin(slope),
    q: wz * (Math.sin(slope) + CB * Math.cos(slope)), slope, X, Z: Zs, mirror,
  };
  const st = { X: 0, Z: 0, Ls: 0, LB: 0, V0: 0, Lbs: 0 };
  const fwd = (H, V, out) => { forward(H, V, p, st); out.X = st.X; out.Z = st.Z; };
  const scaleF = wz * L;

  let H, V, slack = false, residual = 0, evaluations = 0;
  if (X <= 1e-10 * L && !(seabed && L > Zs)) {
    // Vertical line (taut, or hanging in a bight when there is no seabed): H = 0 and V follows from the
    // Z equation alone, which is monotonic in V.
    H = 1e-9 * scaleF;
    const out = { X: 0, Z: 0 };
    const big = scaleF + EA * (Math.abs(Zs) / L + 2);
    let lo = -big, hi = big;
    for (let i = 0; i < 200; i++) {
      const mid = 0.5 * (lo + hi);
      if (mid === lo || mid === hi) break;
      fwd(H, mid, out);
      if (out.Z < Zs) lo = mid; else hi = mid;
    }
    V = 0.5 * (lo + hi);
    slack = L > Math.abs(Zs);
  } else {
    const g = o.guess && o.guess.H > 0 ? { H: o.guess.H, V: mirror ? -o.guess.V : o.guess.V }
      : catenaryInitialGuess({ X, Z: Zs, L, w: wz, EA });
    const stretchMax = 1 + 2 * scaleF / EA;
    const checkSlack = seabed && L * stretchMax >= 0.9 * (X + Math.abs(Zs));
    const r = solveShooting(fwd, {
      X, Z: Zs, H0: g.H, V0: g.V, scaleL: L, scaleF, Hmax: 1e3 * (scaleF + EA), checkSlack,
    });
    H = r.H; V = r.V; slack = r.slack; residual = r.residual; evaluations = r.iterations;
  }
  forward(H, V, p, st);

  // Tangent stiffness: inverse of the flexibility d(X, Z)/d(H, V) by central differences.
  let k11 = 0, k12 = 0, k21 = 0, k22 = 0;
  if (!slack) {
    const t = { X: 0, Z: 0, Ls: 0, LB: 0, V0: 0, Lbs: 0 };
    const dh = 1e-6 * H, dv = 1e-6 * Math.max(Math.abs(V), 1e-3 * scaleF);
    forward(H + dh, V, p, t); const xp = t.X, zp = t.Z;
    forward(H - dh, V, p, t);
    const a11 = (xp - t.X) / (2 * dh), a21 = (zp - t.Z) / (2 * dh);
    forward(H, V + dv, p, t); const xq = t.X, zq = t.Z;
    forward(H, V - dv, p, t);
    const a12 = (xq - t.X) / (2 * dv), a22 = (zq - t.Z) / (2 * dv);
    const det = a11 * a22 - a12 * a21;
    k11 = a22 / det; k12 = -a12 / det; k21 = -a21 / det; k22 = a11 / det;
  }

  const sg = mirror ? -1 : 1;
  const grounded = st.LB > 0;
  let aH, aV;
  if (grounded) {
    const TA = Math.max(H / p.cosB - p.q * st.LB, 0);
    aH = TA * p.cosB; aV = TA * p.sinB;
  } else { aH = H; aV = V - wz * L; }
  const Hout = slack ? 0 : H;
  if (slack) { aH = 0; if (grounded) aV = 0; }
  const Vf = sg * V, aVs = sg * aV;
  // Stretched length: L + (1/EA) * integral of T ds.
  const V1 = st.V0 + wz * st.Ls;
  const Hs = Math.max(H, 1e-300);
  const intT = (V1 * hypot(H, V1) - st.V0 * hypot(H, st.V0) + H * H * (asinh(V1 / Hs) - asinh(st.V0 / Hs))) / (2 * wz);
  const sol = {
    H: Hout, V: Vf, T: hypot(Hout, Vf), angle: Math.atan2(Vf, Hout),
    anchorH: aH, anchorV: aVs, anchorT: hypot(aH, aVs), anchorAngle: Math.atan2(aVs, aH),
    upliftAngle: (aH > 0 || aVs !== 0) ? Math.atan2(aVs, aH) - (seabed ? slope : 0) : 0,
    groundedLength: st.LB, suspendedLength: st.Ls,
    touchdownX: grounded ? Math.min(st.Lbs * p.cosB, X) : 0, touchdownZ: grounded ? Math.min(st.Lbs * p.cosB, X) * p.tanB : 0,
    stretchedLength: st.Lbs + st.Ls + intT / EA,
    stiffness: [[k11, sg * k12], [sg * k21, k22]],
    slack, converged: true, residual, evaluations,
    params: { ...p, H, V, Ls: st.Ls, LB: st.LB, V0: st.V0, Lbs: st.Lbs },
  };
  if (o.nProfile > 0) sol.profile = catenaryProfile(sol, o.nProfile);
  return sol;
}

function weightless(o, chord) {
  const { X, Z, L, EA } = o;
  const taut = chord > L;
  const T = taut ? EA * (chord - L) / L : 0;
  const cx = chord > 0 ? X / chord : 1, cz = chord > 0 ? Z / chord : 0;
  const H = T * cx, V = T * cz;
  let K = [[0, 0], [0, 0]];
  if (taut) {
    // K = (EA/L) e e^T + (T/chord)(I - e e^T)
    const ka = EA / L, kg = T / chord;
    K = [[ka * cx * cx + kg * (1 - cx * cx), (ka - kg) * cx * cz], [(ka - kg) * cx * cz, ka * cz * cz + kg * (1 - cz * cz)]];
  }
  const sol = {
    H, V, T, angle: Math.atan2(Z, X), anchorH: H, anchorV: V, anchorT: T, anchorAngle: Math.atan2(Z, X),
    upliftAngle: Math.atan2(Z, X), groundedLength: 0, suspendedLength: L, touchdownX: 0, touchdownZ: 0,
    stretchedLength: taut ? chord : L, stiffness: K, slack: !taut, converged: true, residual: 0,
    params: { weightless: true, X, Z, L, EA, T },
  };
  if (o.nProfile > 0) sol.profile = catenaryProfile(sol, o.nProfile);
  return sol;
}

/**
 * Coordinates and tension along a solved line at n points equally spaced in unstretched arc length.
 * Reference: Jonkman (2007) Eqs. 2-41 to 2-44 (profile of the suspended and grounded parts).
 * @param {object} sol result of solveCatenary
 * @param {number} n number of points (>= 2)
 * @returns {{s:Float64Array,x:Float64Array,z:Float64Array,T:Float64Array}} unstretched arc length from the anchor
 *   [m], horizontal and vertical coordinates relative to the anchor [m], tension [N]
 */
export function catenaryProfile(sol, n) {
  need(Number.isInteger(n) && n >= 2, 'n must be an integer >= 2');
  const p = sol.params;
  const s = new Float64Array(n), x = new Float64Array(n), z = new Float64Array(n), T = new Float64Array(n);
  if (p.weightless) {
    for (let i = 0; i < n; i++) {
      const f = i / (n - 1);
      s[i] = f * p.L; x[i] = f * p.X; z[i] = f * p.Z; T[i] = p.T;
    }
    return { s, x, z, T };
  }
  const { L, w, EA, H, V0, LB, Lbs, q, cosB, sinB, mirror } = p;
  const T0 = H / cosB;
  const sZero = (q > 0) ? Math.max(LB - T0 / q, 0) : 0; // zero-tension length next to the anchor
  const Hs = Math.max(H, 1e-300);
  const a0 = asinh(V0 / Hs), h0 = hypot(H, V0);
  const xg = Lbs * cosB, zg = Lbs * sinB;
  for (let i = 0; i < n; i++) {
    const si = (i / (n - 1)) * L;
    s[i] = si;
    let xi, zi, Ti;
    if (si <= LB) {
      let e = 0;
      Ti = 0;
      if (si > sZero) {
        e = (T0 - q * LB) * (si - sZero) + 0.5 * q * (si * si - sZero * sZero);
        Ti = T0 - q * (LB - si);
      }
      const pos = si + e / EA;
      xi = pos * cosB; zi = pos * sinB;
    } else {
      const sg = si - LB;
      const Vs = V0 + w * sg;
      xi = xg + (H / w) * (asinh(Vs / Hs) - a0) + H * sg / EA;
      zi = zg + (hypot(H, Vs) - h0) / w + (V0 * sg + 0.5 * w * sg * sg) / EA;
      Ti = hypot(H, Vs);
    }
    if (sol.slack && xi > p.X) xi = p.X; // surplus line lies in a heap; keep the plotted profile within the span
    x[i] = xi; z[i] = mirror ? -zi : zi; T[i] = Ti;
  }
  return { s, x, z, T };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Inverse problems                                                                                             */
/* ------------------------------------------------------------------------------------------------------------ */

function scanBracket(f, x0, factor, maxSteps, name) {
  // March geometrically from x0 until f changes sign; returns [a, b].
  let a = x0, fa = f(a);
  for (let i = 0; i < maxSteps; i++) {
    const b = a * factor, fb = f(b);
    if (fa === 0) return [a, a];
    if (fa * fb <= 0) return [a, b];
    a = b; fa = fb;
  }
  throw new RangeError(`${name}: target cannot be reached for this geometry`);
}

/**
 * Unstretched line length that gives a target fairlead tension for a fixed anchor-fairlead geometry.
 * Tension falls monotonically with length while part of the line can rest on the seabed; the search starts
 * from the taut side and stops at the first crossing. Reference: elastic catenary, Jonkman (2007) section 2.5.
 * @param {object} o X, Z [m], w [N/m], EA [N], CB, seabed, slope as in solveCatenary
 * @param {number} o.T target fairlead tension [N]
 * @returns {{L:number,solution:object}} length [m] and the catenary solution at that length
 */
export function lengthForFairleadTension(o) {
  finite(o.T, 'T');
  need(o.T > 0, 'T must be > 0');
  const chord = hypot(o.X, o.Z);
  need(chord > 0, 'X and Z must not both be zero');
  const f = (L) => solveCatenary({ ...o, L, nProfile: 0, guess: undefined }).T - o.T;
  const L0 = chord / (1 + (2 * o.T + Math.abs(o.w) * chord) / o.EA);
  const [a, b] = scanBracket(f, L0, 1.02, 400, 'T');
  const L = a === b ? a : brent(f, a, b, { tol: 1e-12 * chord });
  return { L, solution: solveCatenary({ ...o, L }) };
}

/**
 * Unstretched line length that gives a target fairlead angle (from the horizontal) for a fixed geometry.
 * @param {object} o X, Z [m], w [N/m], EA [N], CB, seabed, slope as in solveCatenary
 * @param {number} o.angle target fairlead angle from horizontal [rad]
 * @returns {{L:number,solution:object}} length [m] and the catenary solution
 */
export function lengthForFairleadAngle(o) {
  finite(o.angle, 'angle');
  const chord = hypot(o.X, o.Z);
  need(chord > 0, 'X and Z must not both be zero');
  need(o.angle > Math.atan2(o.Z, o.X) && o.angle < Math.PI / 2, 'angle must lie between the chord angle and pi/2');
  const f = (L) => solveCatenary({ ...o, L, nProfile: 0, guess: undefined }).angle - o.angle;
  const [a, b] = scanBracket(f, chord * 0.999, 1.02, 400, 'angle');
  const L = a === b ? a : brent(f, a, b, { tol: 1e-12 * chord });
  return { L, solution: solveCatenary({ ...o, L }) };
}

/**
 * Horizontal anchor-fairlead distance (anchor radius) that gives a target fairlead tension (pretension) for a
 * line of fixed length. Tension rises monotonically with the span.
 * @param {object} o Z [m], L [m], w [N/m], EA [N], CB, seabed, slope as in solveCatenary
 * @param {number} o.T target fairlead tension [N]
 * @returns {{X:number,solution:object}} horizontal span [m] and the catenary solution
 */
export function spanForFairleadTension(o) {
  finite(o.T, 'T');
  need(o.T > 0, 'T must be > 0');
  const { Z, L, w, EA } = o;
  const Lmax = L * (1 + (2 * o.T + Math.abs(w) * L) / EA);
  need(Lmax > Math.abs(Z), 'L is too short to reach the fairlead elevation at this tension');
  const hi = Math.sqrt(Lmax * Lmax - Z * Z);
  const lo = Math.max(L - Math.abs(Z), 0) * 0.999 + 1e-9 * L;
  const f = (X) => solveCatenary({ ...o, X, nProfile: 0, guess: undefined }).T - o.T;
  need(f(lo) < 0, 'T is below the tension of the slack (vertically hanging) line');
  const X = brent(f, lo, hi, { tol: 1e-12 * L });
  return { X, solution: solveCatenary({ ...o, X }) };
}

/**
 * Minimum line length for zero anchor uplift at a design fairlead tension (horizontal seabed, no friction
 * needed): the line is tangent to the seabed exactly at the anchor.
 * Inextensible: L = z sqrt(2 T / (w z) - 1) (Faltinsen 1990 Eq. 8.20 with H = T - w z). With EA the unstretched
 * suspended length V/w is found from Z = (T - H)/w + V^2/(2 EA w), V^2 = T^2 - H^2 (Jonkman 2007 Eq. 2-38).
 * @param {object} o
 * @param {number} o.Z fairlead height above the seabed [m]
 * @param {number} o.w submerged weight per unit length [N/m] (> 0)
 * @param {number} o.T design fairlead tension [N] (must exceed w Z)
 * @param {number} [o.EA=Infinity] axial stiffness [N]
 * @returns {{L:number,H:number,V:number,X:number}} minimum unstretched length [m], horizontal and vertical
 *   tension at the fairlead [N], horizontal anchor-fairlead distance at that state [m]
 */
export function minLengthNoUplift({ Z, w, T, EA = Infinity }) {
  finite(Z, 'Z'); finite(w, 'w'); finite(T, 'T');
  need(Z > 0, 'Z must be > 0');
  need(w > 0, 'w must be > 0');
  need(EA > 0, 'EA must be > 0');
  need(T > w * Z, 'T must exceed w*Z (weight of the vertical hanging length)');
  let H;
  if (EA === Infinity) H = T - w * Z;
  else {
    const g = (h) => (T - h) / w + (T * T - h * h) / (2 * EA * w) - Z;
    need(g(0) > 0, 'T must exceed the tension of the vertically hanging stretched line');
    H = brent(g, 0, T, { tol: 1e-13 * T });
  }
  const V = Math.sqrt(T * T - H * H);
  const L = V / w;
  const X = (H / w) * asinh(V / H) + (EA === Infinity ? 0 : H * L / EA);
  return { L, H, V, X };
}
