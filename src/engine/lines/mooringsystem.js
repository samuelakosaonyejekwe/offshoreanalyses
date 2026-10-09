/**
 * Quasi-static mooring system analysis: spread / turret / taut systems of N multi-segment lines attached to a
 * rigid floater, tendon (TLP) sets, single-point hawser model, layout generators, code checks and clearances.
 *
 * Global frame: z up, origin at the still-water level above the mooring centre unless the user chooses otherwise.
 * Pose q = [x, y, z, roll, pitch, yaw] (m, rad) locates the body origin and its orientation,
 * R = Rz(yaw) Ry(pitch) Rx(roll). Generalised mooring load F = [Fx, Fy, Fz, Mx, My, Mz] is expressed in the
 * global axes with moments about the (displaced) body origin.
 *
 * System definition:
 *   { lines: [{ anchor: [x, y, z] (global, m), fairlead: [x, y, z] (body frame, m),
 *               segments: [...], points: [...], seabed?: true, CB?: 0, MBL?: N, name?, turret?: false }],
 *     turretCenter?: [x, y, z] (body frame, m) }
 * with `segments` / `points` as in multiseg.js (anchor -> fairlead). Lines flagged `turret` hang from a turret
 * at turretCenter that keeps its earth-fixed heading: their fairlead offsets from the turret centre do not rotate
 * with the hull and they transmit no yaw moment through the bearing other than (R c) x F.
 *
 * References: API RP 2SK (3rd ed., 2005) sections 5-7; DNV-OS-E301 (2010-2018 editions) Ch.2 Sec.2;
 * Faltinsen (1990) Sea Loads ch. 8; Jonkman (2007) NREL/TP-500-41958 section 2.5.
 */

import { solveMultiSegment, compileLine } from './multiseg.js';
import { solveDense } from '../math/linalg.js';
import { brent } from '../math/roots.js';

const hypot = Math.hypot;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}

function vec3(v, name) {
  need(Array.isArray(v) && v.length === 3 && v.every(Number.isFinite), `${name} must be [x, y, z] of finite numbers`);
}

function checkSystem(system) {
  need(system && Array.isArray(system.lines) && system.lines.length > 0, 'system.lines must be a non-empty array');
  system.lines.forEach((ln, i) => { vec3(ln.anchor, `lines[${i}].anchor`); vec3(ln.fairlead, `lines[${i}].fairlead`); });
}

function pose6(pose) {
  if (!pose) return [0, 0, 0, 0, 0, 0];
  need(pose.length === 6 || pose.length === 3, 'pose must have 6 entries [x, y, z, roll, pitch, yaw] or 3 [x, y, yaw]');
  const q = pose.length === 3 ? [pose[0], pose[1], 0, 0, 0, pose[2]] : Array.from(pose);
  need(q.every(Number.isFinite), 'pose entries must be finite');
  return q;
}

/**
 * Rotation matrix R = Rz(yaw) Ry(pitch) Rx(roll) taking body-frame vectors to the global frame.
 * @param {number} roll [rad] @param {number} pitch [rad] @param {number} yaw [rad]
 * @returns {number[][]} 3x3 rotation matrix [-]
 */
export function rotationMatrix(roll, pitch, yaw) {
  const cr = Math.cos(roll), sr = Math.sin(roll), cp = Math.cos(pitch), sp = Math.sin(pitch);
  const cy = Math.cos(yaw), sy = Math.sin(yaw);
  return [
    [cy * cp, cy * sp * sr - sy * cr, cy * sp * cr + sy * sr],
    [sy * cp, sy * sp * sr + cy * cr, sy * sp * cr - cy * sr],
    [-sp, cp * sr, cp * cr],
  ];
}

const rot = (R, v) => [
  R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2],
  R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
  R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2],
];

/**
 * Mooring line forces on the floater for a given pose (quasi-static: each line is the static multi-segment
 * catenary between its anchor and the displaced fairlead, in its own vertical plane).
 * Reference: API RP 2SK section 5.3 (quasi-static analysis); line model per multiseg.js.
 * Validity: no current/wave load on the lines, no line dynamics, no line-line or line-hull contact.
 * @param {object} system system definition (see file header)
 * @param {number[]} [pose] [x, y, z, roll, pitch, yaw] (m, rad) or [x, y, yaw]; default all zero
 * @param {object} [opts]
 * @param {number[]} [opts.exclude] indices of lines treated as broken (zero force)
 * @param {number} [opts.nProfile=0] profile points per segment returned for each line (3-D coordinates)
 * @param {Array} [opts.warm] caller-owned array used to warm-start the line solves between calls
 * @returns {{F:number[],lines:object[],maxTension:number,maxTensionLine:number}} F = total mooring load
 *   [N, N, N, N m, N m, N m]; per line: index, name, broken, H, V, T (fairlead tension) [N], anchorH, anchorV,
 *   anchorT [N], upliftAngle [rad], groundedLength [m], span X, Z [m], heading [rad] (direction from fairlead to
 *   anchor), fairlead [x, y, z] global [m], force [Fx, Fy, Fz] on the body [N], stiffness 2x2 [N/m],
 *   utilisation T/MBL [-] when MBL is given, solution (multiseg result), points [[x, y, z]] when nProfile > 0
 */
export function lineForces(system, pose, opts = {}) {
  checkSystem(system);
  const q = pose6(pose);
  const R = rotationMatrix(q[3], q[4], q[5]);
  const F = [0, 0, 0, 0, 0, 0];
  const out = [];
  const excl = opts.exclude ?? [];
  const tc = system.turretCenter;
  let tcg = null;
  if (tc) { vec3(tc, 'turretCenter'); tcg = rot(R, tc); }
  let Tmax = 0, iMax = -1;
  system.lines.forEach((ln, i) => {
    if (excl.includes(i)) {
      out.push({ index: i, name: ln.name ?? `line ${i + 1}`, broken: true, H: 0, V: 0, T: 0, force: [0, 0, 0] });
      return;
    }
    let arm, armYaw;
    if (ln.turret && tcg) {
      arm = [tcg[0] + ln.fairlead[0] - tc[0], tcg[1] + ln.fairlead[1] - tc[1], tcg[2] + ln.fairlead[2] - tc[2]];
      armYaw = tcg;
    } else { arm = rot(R, ln.fairlead); armYaw = arm; }
    const fg = [q[0] + arm[0], q[1] + arm[1], q[2] + arm[2]];
    const dx = ln.anchor[0] - fg[0], dy = ln.anchor[1] - fg[1];
    const X = hypot(dx, dy), Z = fg[2] - ln.anchor[2];
    const warm = opts.warm ? opts.warm[i] : undefined;
    // compiled line definitions are kept only in the caller-owned warm-start store, never across calls
    let cl;
    if (opts.warm) {
      if (!opts.warm.compiled) opts.warm.compiled = [];
      cl = opts.warm.compiled[i] ?? (opts.warm.compiled[i] = compileLine(ln));
    } else cl = compileLine(ln);
    const sol = solveMultiSegment({ X, Z, line: cl, guess: warm, nProfile: opts.nProfile ?? 0 });
    if (opts.warm && !sol.slack) opts.warm[i] = { H: sol.H, V: sol.V };
    const ux = X > 0 ? dx / X : 1, uy = X > 0 ? dy / X : 0;
    const f = [sol.H * ux, sol.H * uy, -sol.V];
    F[0] += f[0]; F[1] += f[1]; F[2] += f[2];
    F[3] += arm[1] * f[2] - arm[2] * f[1];
    F[4] += arm[2] * f[0] - arm[0] * f[2];
    F[5] += armYaw[0] * f[1] - armYaw[1] * f[0];
    const rec = {
      index: i, name: ln.name ?? `line ${i + 1}`, broken: false,
      H: sol.H, V: sol.V, T: sol.T, anchorH: sol.anchorH, anchorV: sol.anchorV, anchorT: sol.anchorT,
      upliftAngle: sol.upliftAngle, groundedLength: sol.groundedLength, X, Z, heading: Math.atan2(dy, dx),
      fairlead: fg, force: f, stiffness: sol.stiffness, slack: sol.slack, solution: sol,
    };
    if (ln.MBL > 0) rec.utilisation = sol.T / ln.MBL;
    if (sol.profile) {
      const pr = sol.profile, pts = new Array(pr.x.length);
      // profile x runs from the anchor towards the fairlead
      for (let k = 0; k < pts.length; k++) pts[k] = [ln.anchor[0] - ux * pr.x[k], ln.anchor[1] - uy * pr.x[k], ln.anchor[2] + pr.z[k]];
      rec.points = pts;
    }
    if (sol.T > Tmax) { Tmax = sol.T; iMax = i; }
    out.push(rec);
  });
  return { F, lines: out, maxTension: Tmax, maxTensionLine: iMax };
}

function mat6(K, name) {
  need(Array.isArray(K) && K.length === 6 && K.every((r) => Array.isArray(r) && r.length === 6), `${name} must be a 6x6 matrix`);
}

/**
 * Static equilibrium pose of the moored floater under a mean external load, by damped Newton iteration with a
 * finite-difference Jacobian on the selected degrees of freedom.
 * Residual: F_mooring(q) + Fext - hydroK (q - hydroRef) = 0 for the active DOFs. For 6-DOF equilibrium the
 * caller supplies the hydrostatic stiffness (and includes weight/buoyancy imbalance in Fext, or chooses hydroRef
 * as the free-floating pose under the mooring pretension).
 * Reference: API RP 2SK section 5.3 (mean offset under steady wind, current and mean wave drift).
 * @param {object} system system definition
 * @param {object} [o]
 * @param {number[]} [o.Fext] external load [Fx, Fy, Fz, Mx, My, Mz] in global axes about the body origin [N, N m]
 * @param {number[]} [o.dofs=[0,1,5]] active DOF indices (0 surge ... 5 yaw); others stay at pose0
 * @param {number[]} [o.pose0] starting pose (m, rad)
 * @param {number[][]} [o.hydroK] 6x6 hydrostatic (or other linear) stiffness [N/m, N/rad, N m/rad]
 * @param {number[]} [o.hydroRef] pose at which the hydroK restoring load is zero (default pose0)
 * @param {number[]} [o.exclude] broken line indices
 * @param {number} [o.tol=1e-8] residual tolerance relative to the force scale [-]
 * @param {number} [o.maxIter=60]
 * @returns {{pose:number[],offset:number,offsetDirection:number,residual:number[],iterations:number,
 *   converged:boolean,F:number[],lines:object[],maxTension:number,maxTensionLine:number}} equilibrium pose
 *   (m, rad), horizontal offset from pose0 [m] and its direction [rad], residual load on active DOFs, line results
 */
export function equilibrium(system, o = {}) {
  checkSystem(system);
  const Fext = o.Fext ?? [0, 0, 0, 0, 0, 0];
  need(Array.isArray(Fext) && Fext.length === 6 && Fext.every(Number.isFinite), 'Fext must have 6 finite entries');
  const dofs = o.dofs ?? [0, 1, 5];
  need(dofs.length > 0 && dofs.every((d) => Number.isInteger(d) && d >= 0 && d < 6), 'dofs must be indices 0..5');
  const q0 = pose6(o.pose0);
  const ref = o.hydroRef ? pose6(o.hydroRef) : q0;
  if (o.hydroK) mat6(o.hydroK, 'hydroK');
  const tol = o.tol ?? 1e-8;
  const warm = [];
  const exclude = o.exclude;
  let last = null;

  const resid = (q) => {
    const r = lineForces(system, q, { exclude, warm });
    last = r;
    const res = new Array(dofs.length);
    for (let a = 0; a < dofs.length; a++) {
      const i = dofs[a];
      let v = r.F[i] + Fext[i];
      if (o.hydroK) for (let j = 0; j < 6; j++) v -= o.hydroK[i][j] * (q[j] - ref[j]);
      res[a] = v;
    }
    return res;
  };

  // scales: forces by the sum of line tensions, moments by that times the largest fairlead radius
  const r0 = lineForces(system, q0, { exclude, warm });
  let Fs = 0, Lr = 1, span = Infinity;
  r0.lines.forEach((l) => { if (!l.broken) { Fs += l.T; span = Math.min(span, l.X + Math.abs(l.Z)); } });
  system.lines.forEach((l) => { Lr = Math.max(Lr, hypot(l.fairlead[0], l.fairlead[1], l.fairlead[2])); });
  Fs = Math.max(Fs, hypot(Fext[0], Fext[1], Fext[2]), 1);
  const scale = dofs.map((d) => (d < 3 ? Fs : Fs * Lr));
  const merit = (res) => { let m = 0; for (let a = 0; a < res.length; a++) m += (res[a] / scale[a]) ** 2; return Math.sqrt(m); };
  const maxStep = Math.max(0.1 * (Number.isFinite(span) ? span : 100), 1);

  const q = q0.slice();
  let res = resid(q), f = merit(res), it = 0, converged = f < tol;
  const maxIter = o.maxIter ?? 60;
  const n = dofs.length;
  for (; it < maxIter && !converged; it++) {
    const J = [];
    for (let a = 0; a < n; a++) J.push(new Array(n).fill(0));
    for (let b = 0; b < n; b++) {
      const h = dofs[b] < 3 ? 1e-3 * Math.max(1, 0.01 * maxStep) : 1e-5;
      const qp = q.slice(), qm = q.slice();
      qp[dofs[b]] += h; qm[dofs[b]] -= h;
      const rp = resid(qp), rm = resid(qm);
      for (let a = 0; a < n; a++) J[a][b] = (rp[a] - rm[a]) / (2 * h);
    }
    let dq;
    try { dq = solveDense(J, res.map((v) => -v)); } catch (e) {
      throw new Error('Mooring equilibrium failed: the mooring stiffness is singular for the selected degrees of freedom');
    }
    need(dq.every(Number.isFinite), 'mooring stiffness is singular for the selected degrees of freedom');
    let lim = 1;
    for (let b = 0; b < n; b++) {
      const cap = dofs[b] < 3 ? maxStep : 0.2;
      if (Math.abs(dq[b]) * lim > cap) lim = cap / Math.abs(dq[b]);
    }
    let alpha = lim, ok = false;
    for (let ls = 0; ls < 20; ls++) {
      const qn = q.slice();
      for (let b = 0; b < n; b++) qn[dofs[b]] += alpha * dq[b];
      let rn = null, fn = Infinity;
      try { rn = resid(qn); fn = merit(rn); } catch (e) { fn = Infinity; }
      if (fn < f) { for (let j = 0; j < 6; j++) q[j] = qn[j]; res = rn; f = fn; ok = true; break; }
      alpha *= 0.5;
    }
    if (!ok) break;
    converged = f < tol;
  }
  if (!converged && f < 1e3 * tol) converged = true;
  if (!converged) throw new Error('Mooring equilibrium did not converge: the external load may exceed the restoring capacity');
  const r = lineForces(system, q, { exclude, nProfile: o.nProfile ?? 0 });
  return {
    pose: q, offset: hypot(q[0] - q0[0], q[1] - q0[1]), offsetDirection: Math.atan2(q[1] - q0[1], q[0] - q0[0]),
    residual: res, iterations: it, converged, F: r.F, lines: r.lines, maxTension: r.maxTension, maxTensionLine: r.maxTensionLine,
  };
}

/**
 * Restoring force - offset curve in a given horizontal direction (other DOFs held at pose0).
 * Reference: API RP 2SK section 5.3; Faltinsen (1990) section 8.2.
 * @param {object} system system definition
 * @param {object} o
 * @param {number} o.direction offset direction in the horizontal plane [rad from global x]
 * @param {number[]} o.offsets offsets along that direction [m]
 * @param {number[]} [o.pose0] base pose (m, rad)
 * @param {number[]} [o.exclude] broken line indices
 * @returns {{offset:number[],restoring:number[],transverse:number[],yawMoment:number[],maxTension:number[],
 *   lineTensions:number[][],maxUplift:number[],minGrounded:number[]}} restoring force opposing the offset [N],
 *   transverse force [N], yaw moment [N m], largest fairlead tension [N], fairlead tension per line per offset [N],
 *   largest anchor uplift angle [rad], smallest grounded length [m]
 */
export function restoringCurve(system, o) {
  need(Number.isFinite(o.direction), 'direction must be finite');
  need(Array.isArray(o.offsets) && o.offsets.length > 0, 'offsets must be a non-empty array');
  const q0 = pose6(o.pose0);
  const c = Math.cos(o.direction), s = Math.sin(o.direction);
  const res = { offset: [], restoring: [], transverse: [], yawMoment: [], maxTension: [], lineTensions: [], maxUplift: [], minGrounded: [] };
  const warm = [];
  for (const d of o.offsets) {
    const q = q0.slice();
    q[0] += d * c; q[1] += d * s;
    const r = lineForces(system, q, { exclude: o.exclude, warm });
    res.offset.push(d);
    res.restoring.push(-(r.F[0] * c + r.F[1] * s));
    res.transverse.push(-r.F[0] * s + r.F[1] * c);
    res.yawMoment.push(r.F[5]);
    res.maxTension.push(r.maxTension);
    res.lineTensions.push(r.lines.map((l) => l.T));
    let up = -Infinity, gr = Infinity;
    for (const l of r.lines) if (!l.broken) { up = Math.max(up, l.upliftAngle); gr = Math.min(gr, l.groundedLength); }
    res.maxUplift.push(up); res.minGrounded.push(gr);
  }
  return res;
}

/**
 * Linearised 6x6 mooring stiffness K_ij = -dF_i/dq_j at a pose, by central finite differences of the line forces.
 * K is symmetric where the net mooring load vanishes and the lines are frictionless (conservative); with a net
 * mean load the rotational rows carry the usual load-dependent non-symmetric terms.
 * Reference: Faltinsen (1990) section 8.2; API RP 2SK section 5.3.3.
 * @param {object} system system definition
 * @param {object} [o]
 * @param {number[]} [o.pose] pose (m, rad); default zero
 * @param {number} [o.dx=0.01] translation step [m]
 * @param {number} [o.dtheta=1e-5] rotation step [rad]
 * @param {number[]} [o.exclude] broken line indices
 * @returns {number[][]} 6x6 stiffness [N/m, N/rad, N m/m, N m/rad]
 */
export function stiffnessMatrix(system, o = {}) {
  const q = pose6(o.pose);
  const dx = o.dx ?? 0.01, dth = o.dtheta ?? 1e-5;
  need(dx > 0, 'dx must be > 0');
  need(dth > 0, 'dtheta must be > 0');
  const K = [];
  for (let i = 0; i < 6; i++) K.push(new Array(6).fill(0));
  for (let j = 0; j < 6; j++) {
    const h = j < 3 ? dx : dth;
    const qp = q.slice(), qm = q.slice();
    qp[j] += h; qm[j] -= h;
    const Fp = lineForces(system, qp, { exclude: o.exclude }).F, Fm = lineForces(system, qm, { exclude: o.exclude }).F;
    for (let i = 0; i < 6; i++) K[i][j] = -(Fp[i] - Fm[i]) / (2 * h);
  }
  return K;
}

function cloneSystem(system) {
  return {
    ...system,
    lines: system.lines.map((l) => ({
      ...l, anchor: l.anchor.slice(), fairlead: l.fairlead.slice(),
      segments: l.segments.map((s) => ({ ...s })), points: l.points ? l.points.map((p) => ({ ...p })) : undefined,
    })),
  };
}

/**
 * Pretension set-up: adjust either one segment length of each line or the anchor radius so that every line
 * reaches its target fairlead tension at the given pose.
 * @param {object} system system definition (not modified)
 * @param {object} o
 * @param {number|number[]} o.target target fairlead tension, one value or one per line [N]
 * @param {'length'|'anchorRadius'} [o.mode='length'] quantity adjusted
 * @param {number} [o.segment] index of the segment whose length is adjusted (default: the longest)
 * @param {number[]} [o.pose] pose at which the pretension applies (m, rad)
 * @returns {{system:object,tensions:number[],adjustments:number[]}} adjusted copy of the system, achieved
 *   fairlead tensions [N], and per line the new segment length or the new horizontal anchor-fairlead span [m]
 */
export function setPretension(system, o) {
  checkSystem(system);
  const mode = o.mode ?? 'length';
  need(mode === 'length' || mode === 'anchorRadius', "mode must be 'length' or 'anchorRadius'");
  const sys = cloneSystem(system);
  const n = sys.lines.length;
  const targets = Array.isArray(o.target) ? o.target : new Array(n).fill(o.target);
  need(targets.length === n && targets.every((t) => Number.isFinite(t) && t > 0), 'target must be > 0 (one value or one per line)');
  const pose = pose6(o.pose);
  const base = lineForces(system, pose);
  const tensions = [], adjustments = [];
  sys.lines.forEach((ln, i) => {
    const b = base.lines[i];
    const ux = Math.cos(b.heading), uy = Math.sin(b.heading);
    let iseg = o.segment;
    if (iseg === undefined) { iseg = 0; ln.segments.forEach((s, k) => { if (s.L > ln.segments[iseg].L) iseg = k; }); }
    need(Number.isInteger(iseg) && iseg >= 0 && iseg < ln.segments.length, 'segment index out of range');
    const L0 = ln.segments[iseg].L, X0 = b.X;
    const tension = (v) => {
      const segs = mode === 'length' ? ln.segments.map((s, k) => (k === iseg ? { ...s, L: v } : s)) : ln.segments;
      const X = mode === 'length' ? X0 : v;
      return solveMultiSegment({ X, Z: b.Z, segments: segs, points: ln.points, seabed: ln.seabed, CB: ln.CB, stiffness: false }).T;
    };
    const f = (v) => tension(v) - targets[i];
    const v0 = mode === 'length' ? L0 : X0;
    // tension falls with length and rises with span: march from the current value until the target is bracketed
    let a = v0, fa = f(a);
    const up = mode === 'length' ? fa > 0 : fa < 0;
    let bnd = a, fb = fa, found = fa === 0;
    for (let k = 0; k < 400 && !found; k++) {
      bnd = a * (up ? 1.01 : 1 / 1.01);
      try { fb = f(bnd); } catch (e) { break; }
      if (fa * fb <= 0) { found = true; break; }
      a = bnd; fa = fb;
    }
    need(found, `target pretension of line ${i + 1} cannot be reached by adjusting the ${mode === 'length' ? 'segment length' : 'anchor radius'}`);
    const v = fa === 0 ? a : brent(f, Math.min(a, bnd), Math.max(a, bnd), { tol: 1e-12 * v0 });
    if (mode === 'length') ln.segments[iseg].L = v;
    else { ln.anchor[0] = b.fairlead[0] + ux * v; ln.anchor[1] = b.fairlead[1] + uy * v; }
    tensions.push(tension(v)); adjustments.push(v);
  });
  return { system: sys, tensions, adjustments };
}

/**
 * Intact and one-line-broken equilibria (each line removed in turn) with a transient overshoot estimate.
 *
 * Transient (energy method): immediately after the break the floater is at rest at the intact equilibrium q0
 * while the damaged system has its equilibrium at q1. Moving along the straight path q(s) = q0 + s (q1 - q0)
 * under the constant mean load, the undamped floater comes to rest again where the net work of the unbalanced
 * load vanishes: W(s) = integral_0^s [Fext + F_damaged(q)] . (q1 - q0) ds = 0. For a linear system this gives
 * s = 2 (100 % overshoot); a hardening mooring gives less. Damping and added-mass memory are neglected, so the
 * estimate is an upper bound for the path-constrained motion; `overshootReduction` (0..1) scales the excursion
 * beyond the damaged equilibrium to allow for damping (1 = undamped).
 * Reference: API RP 2SK section 5.5 / 7.2 (damaged condition and transient analysis); energy method from first
 * principles (work-energy theorem).
 * @param {object} system system definition
 * @param {object} [o] Fext, dofs, pose0, hydroK, hydroRef as in equilibrium()
 * @param {number} [o.overshootReduction=1] fraction of the undamped overshoot retained [-]
 * @param {number[]} [o.lines] indices of the lines to break in turn (default all)
 * @returns {{intact:object,damaged:object[]}} intact = equilibrium() result; damaged[k] = { brokenLine,
 *   equilibrium (result), transientPose (m, rad), overshootFactor s [-], transientOffset [m] (from pose0),
 *   transientTensions [N] per line, maxTransientTension [N], maxTransientLine }
 */
export function brokenLineAnalysis(system, o = {}) {
  const intact = equilibrium(system, o);
  const dofs = o.dofs ?? [0, 1, 5];
  const Fext = o.Fext ?? [0, 0, 0, 0, 0, 0];
  const q0 = intact.pose;
  const base = pose6(o.pose0);
  const ref = o.hydroRef ? pose6(o.hydroRef) : base;
  const red = o.overshootReduction ?? 1;
  need(red >= 0 && red <= 1, 'overshootReduction must be within [0, 1]');
  const list = o.lines ?? system.lines.map((_, i) => i);
  const damaged = [];
  for (const k of list) {
    const eq = equilibrium(system, { ...o, exclude: [k], pose0: q0, hydroRef: ref });
    const q1 = eq.pose;
    const d = q1.map((v, j) => v - q0[j]);
    const poseAt = (s) => q0.map((v, j) => v + s * d[j]);
    const warm = [];
    const power = (s) => {
      const q = poseAt(s);
      const r = lineForces(system, q, { exclude: [k], warm });
      let p = 0;
      for (const i of dofs) {
        let v = r.F[i] + Fext[i];
        if (o.hydroK) for (let j = 0; j < 6; j++) v -= o.hydroK[i][j] * (q[j] - ref[j]);
        p += v * d[i];
      }
      return p;
    };
    let sStop = 1;
    const dn = hypot(...d);
    if (dn > 1e-9) {
      // Simpson integration of the work in steps of ds until it changes sign
      const ds = 0.05;
      let s = 0, W = 0, p0 = power(0), found = false;
      for (let it = 0; it < 200; it++) {
        let pm, p1;
        try { pm = power(s + ds / 2); p1 = power(s + ds); } catch (e) { break; }
        const Wn = W + ds * (p0 + 4 * pm + p1) / 6;
        if (s + ds > 1 && Wn <= 0) {
          // W is locally quadratic in s: interpolate with the end slopes (Hermite -> solve by bisection)
          const Wof = (t) => { const pt = p0 + (p1 - p0) * t; return W + ds * t * (p0 + pt) / 2; };
          let lo = 0, hi = 1;
          for (let b = 0; b < 60; b++) { const mid = 0.5 * (lo + hi); if (Wof(mid) > 0) lo = mid; else hi = mid; }
          sStop = s + ds * 0.5 * (lo + hi); found = true; break;
        }
        s += ds; W = Wn; p0 = p1;
      }
      if (!found) sStop = 2;
    }
    const sT = 1 + red * (sStop - 1);
    const qT = poseAt(sT);
    const rT = lineForces(system, qT, { exclude: [k] });
    damaged.push({
      brokenLine: k, equilibrium: eq, transientPose: qT, overshootFactor: sT,
      transientOffset: hypot(qT[0] - base[0], qT[1] - base[1]),
      transientTensions: rT.lines.map((l) => l.T), maxTransientTension: rT.maxTension, maxTransientLine: rT.maxTensionLine,
    });
  }
  return { intact, damaged };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Code checks                                                                                                  */
/* ------------------------------------------------------------------------------------------------------------ */

/** API RP 2SK (3rd ed. 2005) tension limits (60 / 50 / 80 / 70 % MBL) expressed as safety factors on MBL. */
export const API_2SK_SAFETY_FACTORS = Object.freeze({
  intact: Object.freeze({ dynamic: 1.67, 'quasi-static': 2.0 }),
  damaged: Object.freeze({ dynamic: 1.25, 'quasi-static': 1.43 }),
});

/**
 * Line tension check against the API RP 2SK (3rd ed.) tension-limit table: required safety factor on the minimum
 * breaking load of 1.67 (intact, dynamic analysis), 2.0 (intact, quasi-static), 1.25 (damaged, dynamic),
 * 1.43 (damaged, quasi-static) - i.e. tension limits of 60 / 50 / 80 / 70 % MBL.
 * @param {object} o
 * @param {number} o.tension maximum line tension [N]
 * @param {number} o.MBL minimum breaking load (after any corrosion/wear allowance) [N]
 * @param {'intact'|'damaged'} [o.condition='intact']
 * @param {'dynamic'|'quasi-static'} [o.analysis='dynamic']
 * @param {number} [o.requiredSF] override of the required safety factor [-]
 * @returns {{safetyFactor:number,required:number,uc:number,ok:boolean,ref:string}} achieved and required safety
 *   factor [-], utilisation = required / achieved [-]
 */
export function checkApi2sk({ tension, MBL, condition = 'intact', analysis = 'dynamic', requiredSF }) {
  need(Number.isFinite(tension) && tension >= 0, 'tension must be >= 0');
  need(Number.isFinite(MBL) && MBL > 0, 'MBL must be > 0');
  need(condition in API_2SK_SAFETY_FACTORS, "condition must be 'intact' or 'damaged'");
  need(analysis in API_2SK_SAFETY_FACTORS.intact, "analysis must be 'dynamic' or 'quasi-static'");
  const required = requiredSF ?? API_2SK_SAFETY_FACTORS[condition][analysis];
  need(required > 0, 'requiredSF must be > 0');
  const sf = tension > 0 ? MBL / tension : Infinity;
  const uc = tension * required / MBL;
  return { safetyFactor: sf, required, uc, ok: uc <= 1, ref: 'API RP 2SK (3rd ed.) tension limits' };
}

/**
 * DNV-OS-E301 partial safety factors (editions 2010-2018, Ch.2 Sec.2 Tables for ULS and ALS).
 * dynamic: [gammaMean, gammaDyn]; quasi-static: single factor applied to mean and dynamic tension alike.
 */
export const DNV_E301_FACTORS = Object.freeze({
  ULS: Object.freeze({
    1: Object.freeze({ dynamic: [1.10, 1.50], 'quasi-static': [1.70, 1.70] }),
    2: Object.freeze({ dynamic: [1.40, 2.10], 'quasi-static': [2.50, 2.50] }),
  }),
  ALS: Object.freeze({
    1: Object.freeze({ dynamic: [1.00, 1.10], 'quasi-static': [1.10, 1.10] }),
    2: Object.freeze({ dynamic: [1.00, 1.25], 'quasi-static': [1.35, 1.35] }),
  }),
});

/**
 * Line tension check in the DNV-OS-E301 partial-factor format: u = (Tmean gammaMean + Tdyn gammaDyn) / Sc <= 1,
 * with characteristic strength Sc = 0.95 Smbs when strength statistics are not available.
 * Reference: DNV-OS-E301 (2010-2018) Ch.2 Sec.2 (design equation, Sc, partial safety factors). The factors
 * differ between editions; pass gammaMean / gammaDyn to apply the edition required by the project.
 * @param {object} o
 * @param {number} o.Tmean characteristic mean tension [N]
 * @param {number} o.Tdyn characteristic dynamic tension (most probable max minus mean) [N]
 * @param {number} o.MBL minimum breaking strength Smbs [N]
 * @param {'ULS'|'ALS'} [o.limitState='ULS']
 * @param {1|2} [o.consequenceClass=1]
 * @param {'dynamic'|'quasi-static'} [o.analysis='dynamic']
 * @param {number} [o.strengthFactor=0.95] Sc / Smbs [-]
 * @param {number} [o.gammaMean] override [-]
 * @param {number} [o.gammaDyn] override [-]
 * @returns {{Td:number,Sc:number,uc:number,ok:boolean,gammaMean:number,gammaDyn:number,ref:string}} design
 *   tension [N], characteristic strength [N], utilisation [-]
 */
export function checkDnvE301({ Tmean, Tdyn, MBL, limitState = 'ULS', consequenceClass = 1, analysis = 'dynamic', strengthFactor = 0.95, gammaMean, gammaDyn }) {
  need(Number.isFinite(Tmean) && Tmean >= 0, 'Tmean must be >= 0');
  need(Number.isFinite(Tdyn) && Tdyn >= 0, 'Tdyn must be >= 0');
  need(Number.isFinite(MBL) && MBL > 0, 'MBL must be > 0');
  need(limitState in DNV_E301_FACTORS, "limitState must be 'ULS' or 'ALS'");
  need(consequenceClass === 1 || consequenceClass === 2, 'consequenceClass must be 1 or 2');
  need(analysis === 'dynamic' || analysis === 'quasi-static', "analysis must be 'dynamic' or 'quasi-static'");
  need(strengthFactor > 0 && strengthFactor <= 1, 'strengthFactor must be within (0, 1]');
  const g = DNV_E301_FACTORS[limitState][consequenceClass][analysis];
  const gm = gammaMean ?? g[0], gd = gammaDyn ?? g[1];
  const Td = Tmean * gm + Tdyn * gd;
  const Sc = strengthFactor * MBL;
  return { Td, Sc, uc: Td / Sc, ok: Td <= Sc, gammaMean: gm, gammaDyn: gd, ref: 'DNV-OS-E301 Ch.2 Sec.2 (partial safety factors)' };
}

/**
 * Anchor uplift check: uplift angle at the anchor against an allowable angle (zero for drag-embedment anchors
 * in the intact condition per API RP 2SK section 6.5 unless the anchor is qualified for uplift).
 * @param {object} o
 * @param {number} o.upliftAngle line angle above the seabed at the anchor [rad]
 * @param {number} [o.allowable=0] allowable uplift angle [rad]
 * @param {number} [o.anchorV] vertical anchor load [N] (reported back)
 * @returns {{upliftAngle:number,allowable:number,ok:boolean,anchorV:number}}
 */
export function checkAnchorUplift({ upliftAngle, allowable = 0, anchorV = 0 }) {
  need(Number.isFinite(upliftAngle), 'upliftAngle must be finite');
  need(Number.isFinite(allowable) && allowable >= 0, 'allowable must be >= 0');
  return { upliftAngle, allowable, ok: upliftAngle <= allowable + 1e-9, anchorV };
}

/**
 * Grounded length check: line length resting on the seabed against a required minimum.
 * @param {object} o
 * @param {number} o.groundedLength [m]
 * @param {number} [o.required=0] required minimum grounded length [m]
 * @returns {{groundedLength:number,required:number,margin:number,ok:boolean}} margin = grounded - required [m]
 */
export function checkGroundedLength({ groundedLength, required = 0 }) {
  need(Number.isFinite(groundedLength) && groundedLength >= 0, 'groundedLength must be >= 0');
  need(Number.isFinite(required) && required >= 0, 'required must be >= 0');
  return { groundedLength, required, margin: groundedLength - required, ok: groundedLength >= required };
}

/**
 * Offset check as a percentage of water depth. The limit is project specific (riser / umbilical / gangway
 * limits) and must be supplied.
 * @param {object} o
 * @param {number} o.offset horizontal offset [m]
 * @param {number} o.depth water depth [m]
 * @param {number} o.limitFraction allowable offset / depth [-] (e.g. 0.08 for 8 % of depth)
 * @returns {{offsetFraction:number,limitFraction:number,uc:number,ok:boolean}}
 */
export function checkOffset({ offset, depth, limitFraction }) {
  need(Number.isFinite(offset) && offset >= 0, 'offset must be >= 0');
  need(Number.isFinite(depth) && depth > 0, 'depth must be > 0');
  need(Number.isFinite(limitFraction) && limitFraction > 0, 'limitFraction must be > 0');
  const fr = offset / depth;
  return { offsetFraction: fr, limitFraction, uc: fr / limitFraction, ok: fr <= limitFraction };
}

/**
 * Significant and maximum values of a narrow-banded Gaussian response from its standard deviation:
 * significant = 2 sigma, maximum = sigma sqrt(2 ln N), N = number of cycles in the storm.
 * Reference: API RP 2SK (3rd ed.) section 5.4 / Appendix (statistics of low- and wave-frequency responses).
 * @param {object} o
 * @param {number} o.sigma standard deviation of the response [response unit]
 * @param {number} o.N number of response cycles in the storm duration [-] (> 1)
 * @returns {{sig:number,max:number}} significant and expected maximum single amplitude [response unit]
 */
export function responseStatistics({ sigma, N }) {
  need(Number.isFinite(sigma) && sigma >= 0, 'sigma must be >= 0');
  need(Number.isFinite(N) && N > 1, 'N must be > 1');
  return { sig: 2 * sigma, max: sigma * Math.sqrt(2 * Math.log(N)) };
}

/**
 * API RP 2SK combination of mean, low-frequency (LF) and wave-frequency (WF) components of a response
 * (offset or tension): max of (mean + LFmax + WFsig) and (mean + LFsig + WFmax).
 * Reference: API RP 2SK (3rd ed.) section 5.4.4 (Eqs. for maximum offset and maximum tension).
 * @param {object} o
 * @param {number} o.mean mean value
 * @param {number} o.lfSig significant LF single amplitude @param {number} o.lfMax maximum LF single amplitude
 * @param {number} o.wfSig significant WF single amplitude @param {number} o.wfMax maximum WF single amplitude
 * @returns {{max:number,case1:number,case2:number,governing:string}} combined maximum (same unit as the inputs)
 */
export function combineApi2sk({ mean, lfSig, lfMax, wfSig, wfMax }) {
  for (const [k, v] of Object.entries({ mean, lfSig, lfMax, wfSig, wfMax })) need(Number.isFinite(v), `${k} must be finite`);
  const case1 = mean + lfMax + wfSig, case2 = mean + lfSig + wfMax;
  return { max: Math.max(case1, case2), case1, case2, governing: case1 >= case2 ? 'mean + LFmax + WFsig' : 'mean + LFsig + WFmax' };
}

/**
 * Quasi-dynamic line tensions: the floater is displaced from the mean pose by the two API RP 2SK motion
 * combinations (LFmax + WFsig and LFsig + WFmax) along a given direction and the line tensions are read from the
 * quasi-static system. The wave-frequency tension increment is multiplied by a dynamic amplification factor
 * `daf` (ratio of dynamic to quasi-static WF tension amplitude, e.g. from lumpedmass.js) to allow for line
 * dynamics. Reference: API RP 2SK (3rd ed.) section 5.4.4 and 7.3 (quasi-static vs dynamic analysis).
 * @param {object} system system definition
 * @param {object} o
 * @param {number[]} o.meanPose mean equilibrium pose (m, rad)
 * @param {number} o.direction direction of the dynamic excursion [rad]
 * @param {{sig:number,max:number}} o.lf low-frequency motion single amplitudes [m]
 * @param {{sig:number,max:number}} o.wf wave-frequency motion single amplitudes [m]
 * @param {number} [o.daf=1] dynamic amplification of the WF tension increment [-] (>= 0)
 * @param {number[]} [o.exclude] broken line indices
 * @returns {{maxOffset:number,lines:Array<{index:number,Tmean:number,case1:number,case2:number,Tmax:number,
 *   Tdyn:number,governing:string}>,maxTension:number,maxTensionLine:number}} dynamic excursion [m]; per line the
 *   mean tension, the two combinations, their maximum and Tdyn = Tmax - Tmean [N]
 */
export function quasiDynamicTension(system, o) {
  const q = pose6(o.meanPose);
  need(Number.isFinite(o.direction), 'direction must be finite');
  for (const k of ['lf', 'wf']) need(o[k] && o[k].sig >= 0 && o[k].max >= 0, `${k} must be { sig, max } with non-negative values`);
  const daf = o.daf ?? 1;
  need(daf >= 0, 'daf must be >= 0');
  const c = Math.cos(o.direction), s = Math.sin(o.direction);
  const at = (d) => { const p = q.slice(); p[0] += d * c; p[1] += d * s; return lineForces(system, p, { exclude: o.exclude }).lines; };
  const mean = at(0);
  const a1 = at(o.lf.max), b1 = at(o.lf.max + o.wf.sig);
  const a2 = at(o.lf.sig), b2 = at(o.lf.sig + o.wf.max);
  const lines = [];
  let Tmax = 0, iMax = -1;
  mean.forEach((m, i) => {
    if (m.broken) return;
    const case1 = a1[i].T + daf * (b1[i].T - a1[i].T), case2 = a2[i].T + daf * (b2[i].T - a2[i].T);
    const T = Math.max(case1, case2, m.T);
    lines.push({ index: i, Tmean: m.T, case1, case2, Tmax: T, Tdyn: T - m.T, governing: case1 >= case2 ? 'mean + LFmax + WFsig' : 'mean + LFsig + WFmax' });
    if (T > Tmax) { Tmax = T; iMax = i; }
  });
  return { maxOffset: Math.max(o.lf.max + o.wf.sig, o.lf.sig + o.wf.max), lines, maxTension: Tmax, maxTensionLine: iMax };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Tendons (TLP)                                                                                                */
/* ------------------------------------------------------------------------------------------------------------ */

function checkTendons(tendons) {
  need(Array.isArray(tendons) && tendons.length > 0, 'tendons must be a non-empty array');
  tendons.forEach((t, i) => {
    for (const k of ['x', 'y', 'L', 'EA', 'T0']) need(Number.isFinite(t[k]), `tendons[${i}].${k} must be finite`);
    need(t.L > 0 && t.EA > 0 && t.T0 >= 0, `tendons[${i}] needs L > 0, EA > 0, T0 >= 0`);
  });
}

/**
 * Linear stiffness of a tendon set about the platform origin (tendons vertical at zero offset).
 * Heave / roll / pitch from the axial stiffness EA/L; surge / sway / yaw from the pretension (inverted-pendulum
 * geometric stiffness T/L). Reference: Faltinsen (1990) section 8 (TLP); Chakrabarti (2005) Handbook of Offshore
 * Engineering ch. 8.
 * @param {Array<{x:number,y:number,L:number,EA:number,T0:number}>} tendons top positions x, y [m], length [m],
 *   axial stiffness EA [N], pretension T0 [N]
 * @returns {{axial:number[],surge:number,sway:number,heave:number,roll:number,pitch:number,yaw:number}} axial
 *   stiffness per tendon [N/m]; surge, sway, heave [N/m]; roll, pitch, yaw [N m/rad]
 */
export function tendonStiffness(tendons) {
  checkTendons(tendons);
  const r = { axial: [], surge: 0, sway: 0, heave: 0, roll: 0, pitch: 0, yaw: 0 };
  for (const t of tendons) {
    const k = t.EA / t.L, kg = t.T0 / t.L;
    r.axial.push(k);
    r.surge += kg; r.sway += kg; r.heave += k;
    r.roll += k * t.y * t.y; r.pitch += k * t.x * t.x; r.yaw += kg * (t.x * t.x + t.y * t.y);
  }
  return r;
}

/**
 * Set-down, restoring force and tendon tensions of a TLP at a horizontal offset (platform kept level).
 * Each tendon of pretensioned length L spans sqrt((L - d)^2 + x^2) at offset x and set-down d; its tension is
 * T0 + (EA/L)(span - L). Vertical equilibrium sum(T cos a) = sum(T0) + rho g Aw d fixes d. In the inextensible
 * limit d = L - sqrt(L^2 - x^2) and F = sum(T) x / L.
 * Reference: Faltinsen (1990) section 8; Chakrabarti (2005) ch. 8 (TLP offset / set-down).
 * @param {object} o
 * @param {Array} o.tendons as in tendonStiffness
 * @param {number} o.offset horizontal offset [m]
 * @param {number} [o.waterplaneStiffness=0] rho g Aw [N/m] (extra buoyancy per metre of set-down)
 * @returns {{setDown:number,restoring:number,tensions:number[],angles:number[],maxTension:number,minTension:number}}
 *   set-down [m], horizontal restoring force [N], tendon tensions [N], tendon angles from vertical [rad]
 */
export function tendonOffset({ tendons, offset, waterplaneStiffness = 0 }) {
  checkTendons(tendons);
  need(Number.isFinite(offset), 'offset must be finite');
  need(waterplaneStiffness >= 0, 'waterplaneStiffness must be >= 0');
  const x = Math.abs(offset);
  const Lmin = Math.min(...tendons.map((t) => t.L));
  let T0sum = 0;
  for (const t of tendons) T0sum += t.T0;
  const state = (d) => {
    let Fv = 0, Fh = 0;
    const T = [], ang = [];
    for (const t of tendons) {
      const span = hypot(t.L - d, x);
      const Ti = Math.max(t.T0 + (t.EA / t.L) * (span - t.L), 0);
      T.push(Ti); ang.push(Math.atan2(x, t.L - d));
      Fv += Ti * (t.L - d) / span; Fh += Ti * x / span;
    }
    return { g: Fv - T0sum - waterplaneStiffness * d, Fh, T, ang };
  };
  let d = 0;
  if (x > 0) {
    need(x < Lmin, 'offset must be smaller than the tendon length');
    const hi = Lmin - Math.sqrt(Lmin * Lmin - x * x); // inextensible set-down bounds the elastic one
    d = state(hi).g >= 0 ? hi : brent((v) => state(v).g, 0, hi, { tol: 1e-13 * Lmin });
  }
  const s = state(d);
  return { setDown: d, restoring: s.Fh, tensions: s.T, angles: s.ang, maxTension: Math.max(...s.T), minTension: Math.min(...s.T) };
}

/**
 * Tendon tensions of a level TLP under a vertical load change and overturning moments (rigid platform, linear
 * tendons): heave z, roll rx, pitch ry follow from the 3x3 stiffness; dT_i = k_i (z + y_i rx - x_i ry).
 * Reference: rigid-body statics; API RP 2T section 7 (minimum / maximum tendon tension components).
 * @param {object} o
 * @param {Array} o.tendons as in tendonStiffness
 * @param {number} [o.Fz=0] additional upward load on the platform [N]
 * @param {number} [o.Mx=0] moment about x [N m] @param {number} [o.My=0] moment about y [N m]
 * @param {number} [o.hydroHeave=0] rho g Aw [N/m]
 * @param {number} [o.hydroRoll=0] hydrostatic roll stiffness [N m/rad] @param {number} [o.hydroPitch=0] [N m/rad]
 * @returns {{heave:number,roll:number,pitch:number,tensions:number[],maxTension:number,minTension:number,slack:boolean}}
 *   heave [m], roll, pitch [rad], tendon tensions [N]; slack = true when any tension would be negative
 */
export function tendonTensions({ tendons, Fz = 0, Mx = 0, My = 0, hydroHeave = 0, hydroRoll = 0, hydroPitch = 0 }) {
  checkTendons(tendons);
  const K = [[hydroHeave, 0, 0], [0, hydroRoll, 0], [0, 0, hydroPitch]];
  for (const t of tendons) {
    const k = t.EA / t.L, a = [1, t.y, -t.x];
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) K[i][j] += k * a[i] * a[j];
  }
  const u = solveDense(K, [Fz, Mx, My]);
  const T = tendons.map((t) => t.T0 + (t.EA / t.L) * (u[0] + t.y * u[1] - t.x * u[2]));
  const mn = Math.min(...T);
  return { heave: u[0], roll: u[1], pitch: u[2], tensions: T, maxTension: Math.max(...T), minTension: mn, slack: mn < 0 };
}

/**
 * Natural period of a single degree of freedom, T = 2 pi sqrt((M + A) / K): the hook for TLP (or any moored
 * floater) natural periods given the tendon / mooring stiffness and the mass plus added mass.
 * @param {object} o
 * @param {number} o.mass structural mass or inertia [kg or kg m^2]
 * @param {number} [o.addedMass=0] added mass or inertia [kg or kg m^2]
 * @param {number} o.stiffness restoring stiffness [N/m or N m/rad]
 * @returns {number} natural period [s]
 */
export function naturalPeriod({ mass, addedMass = 0, stiffness }) {
  need(Number.isFinite(mass) && mass > 0, 'mass must be > 0');
  need(Number.isFinite(addedMass) && addedMass >= 0, 'addedMass must be >= 0');
  need(Number.isFinite(stiffness) && stiffness > 0, 'stiffness must be > 0');
  return 2 * Math.PI * Math.sqrt((mass + addedMass) / stiffness);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Layout generators                                                                                            */
/* ------------------------------------------------------------------------------------------------------------ */

function makeLine(az, o, extra) {
  const c = Math.cos(az), s = Math.sin(az);
  const cx = o.center ? o.center[0] : 0, cy = o.center ? o.center[1] : 0;
  return {
    name: extra.name,
    anchor: [cx + o.anchorRadius * c, cy + o.anchorRadius * s, -o.depth],
    fairlead: [cx + o.fairleadRadius * c, cy + o.fairleadRadius * s, o.fairleadZ ?? 0],
    segments: o.segments.map((sg) => ({ ...sg })),
    points: o.points ? o.points.map((p) => ({ ...p })) : undefined,
    seabed: o.seabed ?? true, CB: o.CB ?? 0, MBL: o.MBL, ...extra,
  };
}

function checkLayout(o) {
  for (const k of ['anchorRadius', 'fairleadRadius', 'depth']) need(Number.isFinite(o[k]) && o[k] >= 0, `${k} must be >= 0`);
  need(o.anchorRadius > o.fairleadRadius, 'anchorRadius must exceed fairleadRadius');
  need(o.depth > 0, 'depth must be > 0');
  need(Array.isArray(o.segments) && o.segments.length > 0, 'segments must be a non-empty array');
}

/**
 * Symmetric spread mooring: n identical lines equally spaced in azimuth.
 * @param {object} o
 * @param {number} o.nLines number of lines (>= 1)
 * @param {number} o.anchorRadius horizontal radius of the anchors from the centre [m]
 * @param {number} o.fairleadRadius horizontal radius of the fairleads [m]
 * @param {number} [o.fairleadZ=0] fairlead elevation in the body frame [m] (negative below the origin)
 * @param {number} o.depth water depth; anchors at z = -depth [m]
 * @param {number} [o.heading=0] azimuth of line 1 [rad]
 * @param {Array} o.segments line make-up (anchor -> fairlead) @param {Array} [o.points] point loads
 * @param {number} [o.MBL] minimum breaking load [N] @param {number} [o.CB=0] seabed friction [-]
 * @returns {object} system definition
 */
export function spreadLayout(o) {
  checkLayout(o);
  need(Number.isInteger(o.nLines) && o.nLines >= 1, 'nLines must be an integer >= 1');
  const lines = [];
  for (let i = 0; i < o.nLines; i++) lines.push(makeLine((o.heading ?? 0) + (2 * Math.PI * i) / o.nLines, o, { name: `line ${i + 1}` }));
  return { depth: o.depth, lines };
}

/**
 * Grouped (clustered) spread mooring: nGroups clusters equally spaced in azimuth, each with linesPerGroup lines
 * fanned symmetrically about the cluster axis.
 * @param {object} o as spreadLayout, with
 * @param {number} o.nGroups number of clusters @param {number} o.linesPerGroup lines per cluster
 * @param {number} o.groupSpread angle between adjacent lines within a cluster [rad]
 * @returns {object} system definition
 */
export function clusterLayout(o) {
  checkLayout(o);
  need(Number.isInteger(o.nGroups) && o.nGroups >= 1, 'nGroups must be an integer >= 1');
  need(Number.isInteger(o.linesPerGroup) && o.linesPerGroup >= 1, 'linesPerGroup must be an integer >= 1');
  need(Number.isFinite(o.groupSpread) && o.groupSpread >= 0, 'groupSpread must be >= 0');
  const lines = [];
  for (let g = 0; g < o.nGroups; g++) {
    const axis = (o.heading ?? 0) + (2 * Math.PI * g) / o.nGroups;
    for (let k = 0; k < o.linesPerGroup; k++) {
      const az = axis + (k - (o.linesPerGroup - 1) / 2) * o.groupSpread;
      lines.push(makeLine(az, o, { name: `group ${g + 1} line ${k + 1}`, group: g }));
    }
  }
  return { depth: o.depth, lines };
}

/**
 * Turret mooring: lines attached to a turret of radius fairleadRadius centred at turretCenter (body frame).
 * The lines are flagged `turret` so that the hull weathervanes about the turret without twisting the mooring.
 * @param {object} o as spreadLayout (or with nGroups / linesPerGroup / groupSpread for clusters), with
 * @param {number[]} o.turretCenter [x, y, z] of the turret in the body frame [m]
 * @returns {object} system definition
 */
export function turretLayout(o) {
  vec3(o.turretCenter, 'turretCenter');
  const base = { ...o, center: [o.turretCenter[0], o.turretCenter[1]], fairleadZ: o.turretCenter[2] + (o.fairleadZ ?? 0) };
  const sys = o.nGroups ? clusterLayout(base) : spreadLayout(base);
  sys.lines.forEach((l) => { l.turret = true; });
  sys.turretCenter = o.turretCenter.slice();
  return sys;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Single-point mooring (CALM hawser), simple model                                                             */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Simple CALM-buoy hawser model: tanker - hawser - buoy in series along the load direction.
 * Static: hawser tension = mean environmental force; hawser extension T/kh, buoy offset T/kb.
 * Dynamic peak (energy method): a tanker surging away at speed v is arrested by the series spring
 * ks = 1/(1/kh + 1/kb): 0.5 (M + A) v^2 = 0.5 ks d^2, so dT = v sqrt((M + A) ks).
 * Linear springs, no damping, no fishtailing or wave-frequency buoy motion: a screening model only (not a
 * substitute for OCIMF / time-domain analysis). Derived from first principles (work-energy).
 * @param {object} o
 * @param {number} o.Fmean mean environmental force on the tanker [N] (>= 0)
 * @param {number} o.hawserStiffness hawser axial stiffness kh [N/m]
 * @param {number} o.buoyStiffness horizontal stiffness of the buoy mooring kb [N/m]
 * @param {number} [o.mass=0] tanker mass [kg] @param {number} [o.addedMass=0] surge added mass [kg]
 * @param {number} [o.surgeVelocity=0] tanker surge speed away from the buoy [m/s]
 * @param {number} [o.MBL] hawser breaking load [N]
 * @returns {{Tmean:number,Tmax:number,dT:number,hawserExtension:number,buoyOffset:number,seriesStiffness:number,
 *   surgePeriod:number|null,utilisation:number|null}} tensions [N], extensions [m], series stiffness [N/m],
 *   surge natural period [s], Tmax / MBL [-]
 */
export function calmHawser({ Fmean, hawserStiffness, buoyStiffness, mass = 0, addedMass = 0, surgeVelocity = 0, MBL }) {
  need(Number.isFinite(Fmean) && Fmean >= 0, 'Fmean must be >= 0');
  need(Number.isFinite(hawserStiffness) && hawserStiffness > 0, 'hawserStiffness must be > 0');
  need(Number.isFinite(buoyStiffness) && buoyStiffness > 0, 'buoyStiffness must be > 0');
  need(mass >= 0 && addedMass >= 0, 'mass and addedMass must be >= 0');
  need(Number.isFinite(surgeVelocity), 'surgeVelocity must be finite');
  const ks = 1 / (1 / hawserStiffness + 1 / buoyStiffness);
  const M = mass + addedMass;
  const dT = Math.abs(surgeVelocity) * Math.sqrt(M * ks);
  const Tmax = Fmean + dT;
  return {
    Tmean: Fmean, Tmax, dT, hawserExtension: Tmax / hawserStiffness, buoyOffset: Tmax / buoyStiffness, seriesStiffness: ks,
    surgePeriod: M > 0 ? 2 * Math.PI * Math.sqrt(M / ks) : null, utilisation: MBL > 0 ? Tmax / MBL : null,
  };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Clearances                                                                                                   */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * Minimum distance between two 3-D line segments A = [P1, Q1] and B = [P2, Q2].
 * Reference: Ericson (2005) Real-Time Collision Detection, section 5.1.9 (closest points of two segments).
 * @param {number[][]} segA [P1, Q1], end points [x, y, z] [m]
 * @param {number[][]} segB [P2, Q2], end points [x, y, z] [m]
 * @returns {{distance:number,s:number,t:number}} distance [m] and the closest-point parameters (0..1) on each
 */
export function segmentDistance(segA, segB) {
  need(Array.isArray(segA) && segA.length === 2 && Array.isArray(segB) && segB.length === 2, 'segments must be [P, Q] pairs');
  const p1 = segA[0], q1 = segA[1], p2 = segB[0], q2 = segB[1];
  const d1 = [q1[0] - p1[0], q1[1] - p1[1], q1[2] - p1[2]];
  const d2 = [q2[0] - p2[0], q2[1] - p2[1], q2[2] - p2[2]];
  const r = [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]];
  const a = d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2];
  const e = d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2];
  const f = d2[0] * r[0] + d2[1] * r[1] + d2[2] * r[2];
  const clamp = (v) => Math.min(Math.max(v, 0), 1);
  let s, t;
  if (a <= 1e-300 && e <= 1e-300) { s = 0; t = 0; } else if (a <= 1e-300) { s = 0; t = clamp(f / e); } else {
    const c = d1[0] * r[0] + d1[1] * r[1] + d1[2] * r[2];
    if (e <= 1e-300) { t = 0; s = clamp(-c / a); } else {
      const b = d1[0] * d2[0] + d1[1] * d2[1] + d1[2] * d2[2];
      const den = a * e - b * b;
      s = den > 1e-14 * a * e ? clamp((b * f - c * e) / den) : 0;
      t = (b * s + f) / e;
      if (t < 0) { t = 0; s = clamp(-c / a); } else if (t > 1) { t = 1; s = clamp((b - c) / a); }
    }
  }
  const dx = r[0] + s * d1[0] - t * d2[0], dy = r[1] + s * d1[1] - t * d2[1], dz = r[2] + s * d1[2] - t * d2[2];
  return { distance: hypot(dx, dy, dz), s, t };
}

/**
 * Minimum distance between two polylines (arrays of [x, y, z]); a single point is allowed for either.
 * @param {number[][]} A @param {number[][]} B vertices [m]
 * @returns {{distance:number,iA:number,iB:number,sA:number,sB:number}} minimum distance [m], indices of the
 *   closest segments and the parameters along them [-]
 */
export function polylineMinDistance(A, B) {
  need(Array.isArray(A) && A.length >= 1 && Array.isArray(B) && B.length >= 1, 'polylines must have at least one vertex');
  let best = { distance: Infinity, iA: 0, iB: 0, sA: 0, sB: 0 };
  const na = Math.max(A.length - 1, 1), nb = Math.max(B.length - 1, 1);
  for (let i = 0; i < na; i++) {
    const a0 = A[i], a1 = A[Math.min(i + 1, A.length - 1)];
    for (let j = 0; j < nb; j++) {
      const r = segmentDistance([a0, a1], [B[j], B[Math.min(j + 1, B.length - 1)]]);
      if (r.distance < best.distance) best = { distance: r.distance, iA: i, iB: j, sA: r.s, sB: r.t };
    }
  }
  return best;
}

/**
 * Clearances between all pairs of mooring lines and between each line and subsea assets (pipelines, risers,
 * umbilicals, structures given as polylines with an exclusion radius), at a given pose.
 * Reference: API RP 2SK section 6.6 (clearance requirements); geometry by polylineMinDistance.
 * @param {object} system system definition
 * @param {object} [o]
 * @param {number[]} [o.pose] pose (m, rad)
 * @param {number} [o.nProfile=40] profile points per segment used for each line
 * @param {Array<{name:string,points:number[][],radius?:number}>} [o.assets] asset polylines [m], radius [m]
 * @param {number[]} [o.exclude] broken line indices
 * @returns {{lineLine:Array<{a:number,b:number,distance:number}>,lineAsset:Array<{line:number,asset:string,
 *   distance:number,clearance:number}>,minLineLine:number,minLineAsset:number}} distances [m]; clearance =
 *   distance - asset radius [m]
 */
export function clearances(system, o = {}) {
  const r = lineForces(system, o.pose, { nProfile: o.nProfile ?? 40, exclude: o.exclude });
  const live = r.lines.filter((l) => !l.broken);
  const lineLine = [], lineAsset = [];
  let mLL = Infinity, mLA = Infinity;
  for (let i = 0; i < live.length; i++) {
    for (let j = i + 1; j < live.length; j++) {
      const d = polylineMinDistance(live[i].points, live[j].points).distance;
      lineLine.push({ a: live[i].index, b: live[j].index, distance: d });
      mLL = Math.min(mLL, d);
    }
    for (const as of o.assets ?? []) {
      need(Array.isArray(as.points) && as.points.length >= 1, 'assets[].points must be a non-empty array of [x, y, z]');
      const d = polylineMinDistance(live[i].points, as.points).distance;
      const cl = d - (as.radius ?? 0);
      lineAsset.push({ line: live[i].index, asset: as.name, distance: d, clearance: cl });
      mLA = Math.min(mLA, cl);
    }
  }
  return { lineLine, lineAsset, minLineLine: mLL, minLineAsset: mLA };
}
