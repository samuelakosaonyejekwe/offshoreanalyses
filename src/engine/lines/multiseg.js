/**
 * Multi-segment mooring line / cable in a vertical plane (chain-wire-chain, chain-polyester-chain, lazy-wave
 * buoyancy sections) with point clump weights and buoys at the segment junctions.
 *
 * Each segment is an elastic catenary (Irvine 1981; Jonkman 2007, NREL/TP-500-41958 section 2.5). The line is
 * solved by shooting on the fairlead tension components (H, V): starting from the fairlead the vertical tension
 * component is carried down segment by segment (dV/ds = w, jump of P at a point load), giving the anchor-fairlead
 * spans X(H, V) and Z(H, V) in closed form.
 *
 * Seabed (horizontal, at the anchor elevation): the line touches down where the vertical tension component
 * reaches zero, provided that everything between that point and the anchor is non-buoyant (w >= 0 and point
 * loads >= 0). Otherwise the zero crossing is a suspended sag bend (lazy-wave) and the line continues. On the
 * grounded part the tension falls at CB w per unit length towards the anchor, never below zero. A clump weight
 * lying on the seabed is carried by the seabed. A sag bend that dips below the seabed is not prevented; it is
 * reported through sagClearance < 0.
 *
 * Frame and sign conventions as in catenary.js: anchor at (0, 0), fairlead at (X, Z); segments are listed from
 * the anchor to the fairlead.
 */

import { solveShooting, catenaryInitialGuess } from './catenary.js';

const asinh = Math.asinh;
const hypot = Math.hypot;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}

/**
 * Validate and pre-process a line definition for repeated solves.
 * @param {object} o
 * @param {Array<{L:number,w:number,EA:number,CB?:number,name?:string}>} o.segments anchor -> fairlead; L unstretched
 *   length [m], w submerged weight per length [N/m] (negative = buoyant), EA [N], CB seabed friction [-]
 * @param {Array<{junction:number,weight:number,name?:string}>} [o.points] point loads; junction k (1..n-1) is the
 *   node between segments k-1 and k; weight is the net submerged weight [N], negative for a buoy
 * @param {boolean} [o.seabed=true] horizontal seabed at the anchor elevation
 * @param {number} [o.CB=0] default seabed friction coefficient for segments without their own [-]
 * @returns {object} compiled line (typed arrays; pass as `line` to solveMultiSegment)
 */
export function compileLine(o) {
  need(Array.isArray(o.segments) && o.segments.length > 0, 'segments must be a non-empty array');
  const n = o.segments.length;
  const L = new Float64Array(n), w = new Float64Array(n), EA = new Float64Array(n), CB = new Float64Array(n);
  const P = new Float64Array(n + 1);
  let Ltot = 0, scaleF = 0, EAmin = Infinity, Wnet = 0;
  o.segments.forEach((sg, i) => {
    need(Number.isFinite(sg.L) && sg.L > 0, `segments[${i}].L must be > 0`);
    need(Number.isFinite(sg.w), `segments[${i}].w must be finite`);
    need(Number.isFinite(sg.EA) && sg.EA > 0, `segments[${i}].EA must be > 0`);
    L[i] = sg.L; w[i] = sg.w; EA[i] = sg.EA; CB[i] = sg.CB ?? o.CB ?? 0;
    need(CB[i] >= 0, `segments[${i}].CB must be >= 0`);
    Ltot += sg.L; scaleF += Math.abs(sg.w) * sg.L; Wnet += sg.w * sg.L; EAmin = Math.min(EAmin, sg.EA);
  });
  for (const pt of o.points ?? []) {
    need(Number.isInteger(pt.junction) && pt.junction >= 1 && pt.junction <= n - 1,
      'points[].junction must be an integer between 1 and segments.length - 1');
    need(Number.isFinite(pt.weight), 'points[].weight must be finite');
    P[pt.junction] += pt.weight; scaleF += Math.abs(pt.weight); Wnet += pt.weight;
  }
  const seabed = o.seabed ?? true;
  const canGround = new Uint8Array(n);
  let ok = seabed;
  for (let i = 0; i < n; i++) {
    if (i > 0 && P[i] < 0) ok = false;
    if (w[i] < 0) ok = false;
    canGround[i] = ok ? 1 : 0;
  }
  if (scaleF === 0) scaleF = 1e-6 * EAmin;
  return {
    n, L, w, EA, CB, P, canGround, seabed, Ltot, scaleF, EAmin, Wnet,
    names: o.segments.map((sg, i) => sg.name ?? `segment ${i + 1}`),
    // work arrays filled by the forward pass
    Vtop: new Float64Array(n), Vbot: new Float64Array(n), Lsus: new Float64Array(n),
    TgTop: new Float64Array(n), TgBot: new Float64Array(n), dxS: new Float64Array(n), dzS: new Float64Array(n),
    dxG: new Float64Array(n), st: { X: 0, Z: 0, anchorH: 0, anchorV: 0, stretch: 0 },
  };
}

// Forward (shooting) pass from the fairlead down to the anchor.
function forward(H, V, c) {
  const { n, L, w, EA, CB, P, canGround } = c;
  let Vt = V, X = 0, Z = 0, grounded = false, Tg = 0, stretch = 0;
  for (let i = n - 1; i >= 0; i--) {
    const Li = L[i], wi = w[i], EAi = EA[i];
    let Lg = 0;
    c.dxS[i] = 0; c.dzS[i] = 0; c.dxG[i] = 0; c.TgTop[i] = 0; c.TgBot[i] = 0;
    if (!grounded && canGround[i] && Vt <= 0) { grounded = true; Tg = H; }
    if (!grounded) {
      let Ls = Li, V0 = Vt - wi * Li;
      if (canGround[i] && wi > 0 && V0 < 0) { Ls = Vt / wi; V0 = 0; }
      const T1 = hypot(H, Vt);
      let dx, dz, sInt;
      if (Math.abs(wi) * Ls <= 1e-12 * T1) {
        // weightless (or vanishing) piece: straight elastic bar
        dx = H * Ls * (1 / T1 + 1 / EAi); dz = Vt * Ls * (1 / T1 + 1 / EAi); sInt = T1 * Ls;
      } else {
        const T0 = hypot(H, V0);
        const a1 = asinh(Vt / H), a0 = asinh(V0 / H);
        dx = (H / wi) * (a1 - a0) + H * Ls / EAi;
        dz = (T1 - T0) / wi + (V0 * Ls + 0.5 * wi * Ls * Ls) / EAi;
        sInt = (Vt * T1 - V0 * T0 + H * H * (a1 - a0)) / (2 * wi);
      }
      X += dx; Z += dz; stretch += sInt / EAi;
      c.dxS[i] = dx; c.dzS[i] = dz;
      c.Vtop[i] = Vt; c.Vbot[i] = V0; c.Lsus[i] = Ls;
      if (Ls < Li) { grounded = true; Tg = H; Lg = Li - Ls; } else Vt = V0 - P[i];
    } else {
      Lg = Li; c.Vtop[i] = 0; c.Vbot[i] = 0; c.Lsus[i] = 0;
    }
    if (Lg > 0) {
      const q = CB[i] * wi;
      c.TgTop[i] = Tg;
      let G;
      if (q > 0 && Tg < q * Lg) { G = Tg * Tg / (2 * q); Tg = 0; } else { G = Tg * Lg - 0.5 * q * Lg * Lg; Tg -= q * Lg; }
      c.TgBot[i] = Tg;
      c.dxG[i] = Lg + G / EAi;
      X += c.dxG[i]; stretch += G / EAi;
    }
  }
  const st = c.st;
  st.X = X; st.Z = Z; st.stretch = stretch;
  st.anchorH = grounded ? Tg : H;
  st.anchorV = grounded ? 0 : c.Vbot[0];
  st.grounded = grounded;
}

/**
 * Static configuration of a multi-segment line between an anchor at (0, 0) and a fairlead at (X, Z).
 *
 * Reference: segment-wise elastic catenary (Irvine 1981 ch. 1; Jonkman 2007 Eqs. 2-35 to 2-38) with point-load
 * jump conditions; shooting on (H, V) by damped Newton with a bisection fallback. Validity: static, linear-elastic
 * segments, no bending stiffness, horizontal seabed at the anchor elevation, grounded clump weights carry no
 * friction, sag bends are not prevented from penetrating the seabed (check sagClearance).
 *
 * @param {object} o
 * @param {number} o.X horizontal anchor-fairlead distance [m] (>= 0)
 * @param {number} o.Z fairlead elevation above the anchor [m]
 * @param {Array} [o.segments] as in compileLine (ignored when o.line is given)
 * @param {Array} [o.points] as in compileLine
 * @param {boolean} [o.seabed=true] seabed contact considered
 * @param {number} [o.CB=0] default seabed friction coefficient [-]
 * @param {object} [o.line] line pre-compiled with compileLine (fast path for repeated solves)
 * @param {{H:number,V:number}} [o.guess] warm start [N]
 * @param {number} [o.nProfile=0] profile points per segment (0 = no profile)
 * @param {boolean} [o.stiffness=true] compute the tangent stiffness
 * @returns {object} H, V, T [N], angle [rad] at the fairlead; anchorH, anchorV, anchorT [N], upliftAngle [rad];
 *   groundedLength, suspendedLength (unstretched) [m]; touchdownX [m]; stretchedLength [m]; sagClearance [m]
 *   (elevation of the lowest suspended sag bend above the anchor level, null when there is none; negative
 *   means the sag bend would penetrate the seabed); stiffness [[dH/dX, dH/dZ], [dV/dX, dV/dZ]] [N/m];
 *   segments[] with name, Ttop, Tbottom [N], angleTop, angleBottom [rad from horizontal], Vtop, Vbottom [N],
 *   groundedLength, suspendedLength [m]; junctions[] with index, weight, Vabove, Vbelow [N]; slack flag;
 *   evaluations = number of forward (shooting) evaluations used by the solve [-];
 *   profile {s, x, z, T, seg} when nProfile > 0.
 */
export function solveMultiSegment(o) {
  const c = o.line ?? compileLine(o);
  const { X, Z } = o;
  need(Number.isFinite(X) && X >= 0, 'X must be finite and >= 0');
  need(Number.isFinite(Z), 'Z must be finite');
  if (c.seabed) need(Z >= -1e-9 * c.Ltot, 'Z places the fairlead below the seabed');
  const st = c.st;
  const fwd = (H, V, out) => { forward(H, V, c); out.X = st.X; out.Z = st.Z; };

  let g = o.guess;
  if (!(g && g.H > 0 && Number.isFinite(g.V))) {
    let wbar = c.Wnet / c.Ltot;
    const wref = c.scaleF / c.Ltot;
    if (Math.abs(wbar) < 0.05 * wref) wbar = 0.05 * wref;
    g = catenaryInitialGuess({ X, Z, L: c.Ltot, w: wbar, EA: c.EAmin });
  }
  const checkSlack = c.seabed && c.Ltot * (1 + 2 * c.scaleF / c.EAmin) >= 0.9 * (X + Math.abs(Z));
  const r = solveShooting(fwd, {
    X, Z, H0: g.H, V0: g.V, scaleL: c.Ltot, scaleF: c.scaleF, Hmax: 1e3 * (c.scaleF + c.EAmin), checkSlack,
  });
  const H = r.H, V = r.V;

  let K = [[0, 0], [0, 0]];
  if (!r.slack && (o.stiffness ?? true)) {
    const dh = 1e-6 * H, dv = 1e-6 * Math.max(Math.abs(V), 1e-3 * c.scaleF);
    forward(H + dh, V, c); const xp = st.X, zp = st.Z;
    forward(H - dh, V, c);
    const a11 = (xp - st.X) / (2 * dh), a21 = (zp - st.Z) / (2 * dh);
    forward(H, V + dv, c); const xq = st.X, zq = st.Z;
    forward(H, V - dv, c);
    const a12 = (xq - st.X) / (2 * dv), a22 = (zq - st.Z) / (2 * dv);
    const det = a11 * a22 - a12 * a21;
    K = [[a22 / det, -a12 / det], [-a21 / det, a11 / det]];
  }
  forward(H, V, c);

  const n = c.n;
  const segments = [];
  let grounded = 0, suspended = 0, tdx = 0;
  // Elevation of suspended sag bends (vertical tension component changing sign from negative below to positive
  // above), walking down from the fairlead. The touchdown point itself is not a sag bend.
  let zTop = Z, minClr = Infinity;
  for (let i = n - 1; i >= 0; i--) {
    const Ls = c.Lsus[i];
    if (!(Ls > 0)) break;
    const zb = zTop - c.dzS[i];
    const V0 = c.Vbot[i], wi = c.w[i];
    if (wi > 0 && V0 < 0 && c.Vtop[i] > 0) {
      minClr = Math.min(minClr, zb + (H - hypot(H, V0)) / wi - V0 * V0 / (2 * wi * c.EA[i]));
    }
    if (i > 0 && Ls >= c.L[i] && c.Lsus[i - 1] > 0 && V0 > 0 && V0 - c.P[i] < 0) minClr = Math.min(minClr, zb);
    zTop = zb;
  }
  for (let i = 0; i < n; i++) {
    const Ls = c.Lsus[i], Lg = c.L[i] - Ls;
    grounded += Lg; suspended += Ls;
    tdx += c.dxG[i];
    const Vt = c.Vtop[i], Vb = c.Vbot[i];
    const onBed = Lg > 0;
    segments.push({
      index: i, name: c.names[i],
      Ttop: Ls > 0 ? hypot(H, Vt) : c.TgTop[i],
      Tbottom: onBed ? c.TgBot[i] : hypot(H, Vb),
      angleTop: Ls > 0 ? Math.atan2(Vt, H) : 0,
      angleBottom: onBed ? 0 : Math.atan2(Vb, H),
      Vtop: Vt, Vbottom: onBed ? 0 : Vb,
      groundedLength: Lg, suspendedLength: Ls,
    });
  }
  const junctions = [];
  for (let k = 1; k < n; k++) {
    junctions.push({
      index: k, weight: c.P[k], suspended: c.Lsus[k - 1] >= c.L[k - 1] * (1 - 1e-14) && c.Lsus[k] > 0,
      Vabove: c.Vbot[k], Vbelow: c.Vtop[k - 1],
    });
  }
  const aH = r.slack && st.grounded ? 0 : st.anchorH, aV = st.anchorV;
  const Hout = r.slack ? 0 : H;
  const sol = {
    H: Hout, V, T: hypot(Hout, V), angle: Math.atan2(V, Hout),
    anchorH: aH, anchorV: aV, anchorT: hypot(aH, aV), upliftAngle: (aH > 0 || aV !== 0) ? Math.atan2(aV, aH) : 0,
    groundedLength: grounded, suspendedLength: suspended, touchdownX: Math.min(tdx, X),
    stretchedLength: c.Ltot + st.stretch, sagClearance: minClr === Infinity ? null : minClr,
    stiffness: K, segments, junctions, slack: r.slack, converged: true, residual: r.residual, evaluations: r.iterations,
  };
  if (o.nProfile > 0) sol.profile = profile(c, H, X, Z, o.nProfile, r.slack);
  return sol;
}

function profile(c, H, X, Z, nPer, slack) {
  const n = c.n;
  const m = n * (nPer + 1);
  const s = new Float64Array(m), x = new Float64Array(m), z = new Float64Array(m), T = new Float64Array(m);
  const seg = new Int32Array(m);
  // positions of the top of each segment, walking down from the fairlead
  let xTop = X, zTop = Z, sTop = c.Ltot;
  for (let i = n - 1; i >= 0; i--) {
    const Li = c.L[i], wi = c.w[i], EAi = c.EA[i], Ls = c.Lsus[i], Lg = Li - Ls;
    const V0 = c.Vbot[i];
    const xb = xTop - c.dxS[i], zb = zTop - c.dzS[i]; // bottom of the suspended piece
    const T0 = hypot(H, V0), a0 = asinh(V0 / H);
    const q = c.CB[i] * wi, Tg = c.TgTop[i];
    for (let j = 0; j <= nPer; j++) {
      const k = i * (nPer + 1) + j;
      const sl = (j / nPer) * Li; // local arc length from the segment bottom
      seg[k] = i; s[k] = sTop - Li + sl;
      if (sl >= Lg && Ls > 0) {
        const sg = sl - Lg;
        const Vs = V0 + wi * sg;
        const Ts = hypot(H, Vs);
        if (Math.abs(wi) * Ls <= 1e-12 * hypot(H, c.Vtop[i])) {
          x[k] = xb + H * sg * (1 / Ts + 1 / EAi); z[k] = zb + Vs * sg * (1 / Ts + 1 / EAi);
        } else {
          x[k] = xb + (H / wi) * (asinh(Vs / H) - a0) + H * sg / EAi;
          z[k] = zb + (Ts - T0) / wi + (V0 * sg + 0.5 * wi * sg * sg) / EAi;
        }
        T[k] = Ts;
      } else {
        // grounded: distance down from the top of the grounded piece
        const d = Lg - sl;
        let e, Tk;
        if (q > 0 && d > Tg / q) { e = Tg * Tg / (2 * q); Tk = 0; } else { e = Tg * d - 0.5 * q * d * d; Tk = Tg - q * d; }
        x[k] = xb - (d + e / EAi); z[k] = zb; T[k] = Tk;
        if (slack) x[k] = Math.min(Math.max(x[k], 0), X);
      }
    }
    xTop = xb - c.dxG[i]; zTop = zb; sTop -= Li;
  }
  return { s, x, z, T, seg };
}
