/**
 * Moored ship at a berth: wind and current forces from coefficient curves, elastic mooring lines (with tails
 * and pretension) to bollards / hooks, nonlinear fenders, three-degree-of-freedom static equilibrium, and checks
 * of line, bollard and fender loads and of the vessel excursion at the manifold.
 *
 * Frames. Global (berth) frame: x along the berth face, y pointing away from the berth into the water, z up
 * from chart datum. Ship frame: x forward, y to port, origin at midships on the centreline; fairlead heights are
 * given above the keel, so their global elevation is waterLevel - draft + z.
 * Angle of attack theta: direction the wind or current comes FROM, measured from the bow, positive towards port
 * (0 = head-on, 90 = on the port beam, 180 = from astern), in degrees. (OCIMF tabulates against an angle measured
 * from the stern; convert when entering published curves.)
 *
 * Force definitions (OCIMF form):
 *   wind:    Fx = 0.5 rho_a Cx V^2 A_T,   Fy = 0.5 rho_a Cy V^2 A_L,   Mz = 0.5 rho_a Cxy V^2 A_L LBP
 *   current: Fx = 0.5 rho_w Cx V^2 LBP T, Fy = 0.5 rho_w Cy V^2 LBP T, Mz = 0.5 rho_w Cxy V^2 LBP^2 T
 * with Fx along ship x, Fy along ship y and Mz positive bow-to-port.
 *
 * References: OCIMF Mooring Equipment Guidelines (MEG3 2008 / MEG4 2018) - force definitions, standard
 * environmental criteria and line load limits; PIANC WG 33 / WG 211 (fender reaction curves); BS 6349-4.
 * The built-in coefficient curves are generic smooth shapes for preliminary work only and are NOT the OCIMF or
 * SIGTTO data sets: supply ship-specific curves for design.
 */

import { solveDense } from '../math/linalg.js';

/** One international knot in m/s. */
export const KNOT = 1852 / 3600;

const DEG = Math.PI / 180;
const hypot = Math.hypot;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}

function interp(xs, ys, x) {
  const n = xs.length;
  if (x <= xs[0]) return ys[0];
  if (x >= xs[n - 1]) return ys[n - 1];
  let k = 0;
  while (x > xs[k + 1]) k++;
  return ys[k] + (ys[k + 1] - ys[k]) * (x - xs[k]) / (xs[k + 1] - xs[k]);
}

/**
 * Generic coefficient magnitudes by vessel type: Cx = -Cx0 cos(theta), Cy = -Cy0 sin(theta),
 * Cxy = -Cxy0 sin(2 theta). Placeholders of realistic order of magnitude (deep-water current values), not
 * published data; every value can be overridden or replaced by tabulated curves.
 */
export const GENERIC_COEFFICIENTS = Object.freeze({
  tanker: Object.freeze({ wind: { Cx0: 0.9, Cy0: 0.9, Cxy0: 0.12 }, current: { Cx0: 0.04, Cy0: 0.6, Cxy0: 0.15 } }),
  lng: Object.freeze({ wind: { Cx0: 1.0, Cy0: 1.0, Cxy0: 0.14 }, current: { Cx0: 0.04, Cy0: 0.6, Cxy0: 0.15 } }),
  container: Object.freeze({ wind: { Cx0: 0.8, Cy0: 0.95, Cxy0: 0.14 }, current: { Cx0: 0.04, Cy0: 0.6, Cxy0: 0.15 } }),
  bulk: Object.freeze({ wind: { Cx0: 0.9, Cy0: 0.85, Cxy0: 0.12 }, current: { Cx0: 0.04, Cy0: 0.6, Cxy0: 0.15 } }),
});

/**
 * Generic coefficient curves (clearly not OCIMF data) tabulated every 10 degrees of the angle of attack.
 * @param {object} o
 * @param {'tanker'|'lng'|'container'|'bulk'} [o.type='tanker']
 * @param {'wind'|'current'} o.kind
 * @param {number} [o.Cx0] @param {number} [o.Cy0] @param {number} [o.Cxy0] magnitude overrides [-]
 * @param {number} [o.shallowFactor=1] multiplier on the lateral and yaw current coefficients for small
 *   under-keel clearance [-] (deep water 1; lateral current loads grow several-fold as depth/draft -> 1.1)
 * @returns {{angleDeg:number[],Cx:number[],Cy:number[],Cxy:number[],generic:true}} curves over 0..360 deg
 */
export function genericCoefficientCurves({ type = 'tanker', kind, Cx0, Cy0, Cxy0, shallowFactor = 1 }) {
  need(type in GENERIC_COEFFICIENTS, "type must be 'tanker', 'lng', 'container' or 'bulk'");
  need(kind === 'wind' || kind === 'current', "kind must be 'wind' or 'current'");
  need(shallowFactor > 0, 'shallowFactor must be > 0');
  const g = GENERIC_COEFFICIENTS[type][kind];
  const sf = kind === 'current' ? shallowFactor : 1;
  const a = Cx0 ?? g.Cx0, b = (Cy0 ?? g.Cy0) * sf, c = (Cxy0 ?? g.Cxy0) * sf;
  const out = { angleDeg: [], Cx: [], Cy: [], Cxy: [], generic: true };
  for (let d = 0; d <= 360; d += 10) {
    out.angleDeg.push(d);
    out.Cx.push(-a * Math.cos(d * DEG)); out.Cy.push(-b * Math.sin(d * DEG)); out.Cxy.push(-c * Math.sin(2 * d * DEG));
  }
  return out;
}

function checkCurves(c, name) {
  need(c && Array.isArray(c.angleDeg) && c.angleDeg.length >= 2, `${name}.angleDeg must be an array`);
  for (const k of ['Cx', 'Cy', 'Cxy']) need(Array.isArray(c[k]) && c[k].length === c.angleDeg.length, `${name}.${k} must match angleDeg in length`);
  for (let i = 1; i < c.angleDeg.length; i++) need(c.angleDeg[i] > c.angleDeg[i - 1], `${name}.angleDeg must be ascending`);
  const last = c.angleDeg[c.angleDeg.length - 1];
  need(c.angleDeg[0] === 0 && (last === 180 || last === 360), `${name}.angleDeg must span 0..180 or 0..360 degrees`);
}

function coeffAt(c, thetaDeg) {
  let t = ((thetaDeg % 360) + 360) % 360;
  const half = c.angleDeg[c.angleDeg.length - 1] === 180;
  let sgn = 1;
  if (half && t > 180) { t = 360 - t; sgn = -1; } // port/starboard symmetry: Cx even, Cy and Cxy odd
  return { Cx: interp(c.angleDeg, c.Cx, t), Cy: sgn * interp(c.angleDeg, c.Cy, t), Cxy: sgn * interp(c.angleDeg, c.Cxy, t) };
}

/**
 * Wind or current force on the vessel in ship axes from coefficient curves (OCIMF force definitions).
 * Reference: OCIMF MEG (wind and current drag force equations); see the file header for the reference areas.
 * @param {object} o
 * @param {'wind'|'current'} o.kind
 * @param {number} o.speed wind speed (at the reference height of the coefficients) or current speed [m/s]
 * @param {number} o.thetaDeg angle of attack (direction the flow comes from, from the bow towards port) [deg]
 * @param {object} o.vessel { LBP [m], draft [m] (current), AT transverse (head-on) and AL lateral windage [m^2] }
 * @param {object} [o.curves] { angleDeg, Cx, Cy, Cxy } over 0..180 (symmetric ship) or 0..360; default generic
 * @param {string} [o.type='tanker'] vessel type for the generic curves
 * @param {number} [o.rho] fluid density [kg/m^3]; default 1.28 (air, value used with the OCIMF wind
 *   coefficients) or 1025 (sea water)
 * @param {number} [o.shallowFactor=1] see genericCoefficientCurves (generic current curves only)
 * @returns {{Fx:number,Fy:number,Mz:number,Cx:number,Cy:number,Cxy:number,generic:boolean}} surge and sway
 *   force [N], yaw moment [N m] in ship axes; the coefficients used
 */
export function environmentalForce({ kind, speed, thetaDeg, vessel, curves, type = 'tanker', rho, shallowFactor = 1 }) {
  need(kind === 'wind' || kind === 'current', "kind must be 'wind' or 'current'");
  need(Number.isFinite(speed) && speed >= 0, 'speed must be >= 0');
  need(Number.isFinite(thetaDeg), 'thetaDeg must be finite');
  need(vessel && vessel.LBP > 0, 'vessel.LBP must be > 0');
  let c = curves;
  if (c) checkCurves(c, 'curves'); else c = genericCoefficientCurves({ type, kind, shallowFactor });
  const k = coeffAt(c, thetaDeg);
  let q, Ax, Ay;
  if (kind === 'wind') {
    need(vessel.AT > 0 && vessel.AL > 0, 'vessel.AT and vessel.AL must be > 0');
    q = 0.5 * (rho ?? 1.28) * speed * speed; Ax = vessel.AT; Ay = vessel.AL;
  } else {
    need(vessel.draft > 0, 'vessel.draft must be > 0');
    q = 0.5 * (rho ?? 1025) * speed * speed; Ax = vessel.LBP * vessel.draft; Ay = Ax;
  }
  return { Fx: q * k.Cx * Ax, Fy: q * k.Cy * Ay, Mz: q * k.Cxy * Ay * vessel.LBP, ...k, generic: !curves };
}

/**
 * Line load limits as a fraction of the line MBL: 55 % steel wire, 50 % synthetic (other than polyamide),
 * 45 % polyamide (nylon) - OCIMF MEG3 section 1; MEG4 expresses the same idea as a working load limit of 50 % of
 * the ship design MBL for synthetic lines and 55 % for steel wire. Override per line with `limitFraction`.
 */
export const MEG_LINE_LIMITS = Object.freeze({ steel: 0.55, synthetic: 0.50, polyamide: 0.45 });

/**
 * OCIMF standard environmental criteria for tanker berth moorings as a list of load cases: 60 kn wind from any
 * direction, simultaneously with 3 kn current at 0 / 180 deg, or 2 kn at 10 / 170 deg, or 0.75 kn from the
 * direction of maximum beam current loading (90 deg), each on either side.
 * Reference: OCIMF MEG3 section 2 / MEG4 section 3 (standard environmental criteria; wind is the 30-second
 * average at 10 m elevation).
 * @param {object} [o]
 * @param {number} [o.windStepDeg=30] spacing of the wind directions around the compass [deg]
 * @param {number} [o.windKn=60] wind speed [kn]
 * @returns {Array<{name:string,wind:{speed:number,thetaDeg:number},current:{speed:number,thetaDeg:number}}>}
 *   speeds in m/s, angles of attack relative to the bow [deg]
 */
export function ocimfStandardCases({ windStepDeg = 30, windKn = 60 } = {}) {
  need(windStepDeg > 0 && windStepDeg <= 180, 'windStepDeg must be within (0, 180]');
  need(windKn >= 0, 'windKn must be >= 0');
  const currents = [[3, 0], [3, 180], [2, 10], [2, 350], [2, 170], [2, 190], [0.75, 90], [0.75, 270]];
  const cases = [];
  for (let w = 0; w < 360 - 1e-9; w += windStepDeg) {
    for (const [kn, th] of currents) {
      cases.push({
        name: `wind ${windKn} kn from ${w} deg, current ${kn} kn from ${th} deg`,
        wind: { speed: windKn * KNOT, thetaDeg: w }, current: { speed: kn * KNOT, thetaDeg: th },
      });
    }
  }
  return cases;
}

function prepare(o) {
  const v = o.vessel;
  need(v && v.LBP > 0 && v.beam > 0, 'vessel.LBP and vessel.beam must be > 0');
  need(o.pose0 && Number.isFinite(o.pose0.x) && Number.isFinite(o.pose0.y), 'pose0 must give { x, y, headingDeg }');
  const side = v.side ?? 'starboard';
  need(side === 'port' || side === 'starboard', "vessel.side must be 'port' or 'starboard'");
  const ref = { waterLevel: o.reference?.waterLevel ?? 0, draft: o.reference?.draft ?? v.draft ?? 0 };
  const cond = { waterLevel: o.condition?.waterLevel ?? ref.waterLevel, draft: o.condition?.draft ?? ref.draft };
  const pose0 = [o.pose0.x, o.pose0.y, (o.pose0.headingDeg ?? 0) * DEG];
  const lines = (o.lines ?? []).map((ln, i) => {
    const nm = ln.name ?? `line ${i + 1}`;
    need(Array.isArray(ln.fairlead) && ln.fairlead.length === 3 && ln.fairlead.every(Number.isFinite), `lines[${i}].fairlead must be [x, y, z above keel]`);
    need(Array.isArray(ln.bollard) && ln.bollard.length === 3 && ln.bollard.every(Number.isFinite), `lines[${i}].bollard must be [x, y, z]`);
    need(ln.MBL > 0, `lines[${i}].MBL must be > 0`);
    const pre = ln.pretension ?? 0;
    need(pre >= 0 && pre < ln.MBL, `lines[${i}].pretension must be within [0, MBL)`);
    const material = ln.material ?? 'synthetic';
    const limit = ln.limitFraction ?? o.limits?.[material] ?? MEG_LINE_LIMITS[material];
    need(limit > 0 && limit <= 1, `lines[${i}]: no load limit for material '${material}' (give limitFraction)`);
    const rec = { name: nm, type: ln.type ?? 'line', def: ln, material, limit, MBL: ln.MBL, pre, bollardId: ln.bollardId ?? ln.bollard.join(',') };
    const lRef = lineGeom(rec, pose0, ref).len;
    if (ln.curve) {
      // nonlinear load-elongation: load fraction of MBL against strain of the whole line
      const c = ln.curve;
      need(Array.isArray(c.strain) && Array.isArray(c.loadFraction) && c.strain.length === c.loadFraction.length && c.strain.length >= 2,
        `lines[${i}].curve needs strain[] and loadFraction[] of equal length`);
      for (let k = 1; k < c.strain.length; k++) need(c.strain[k] > c.strain[k - 1] && c.loadFraction[k] > c.loadFraction[k - 1], `lines[${i}].curve must be strictly increasing`);
      need(c.strain[0] === 0 && c.loadFraction[0] === 0, `lines[${i}].curve must start at (0, 0)`);
      const ep = interp(c.loadFraction, c.strain, pre / ln.MBL);
      rec.L0 = lRef / (1 + ep);
      rec.tension = (len) => (len > rec.L0 ? ln.MBL * extrap(c.strain, c.loadFraction, (len - rec.L0) / rec.L0) : 0);
      rec.k = null;
    } else {
      let k = ln.stiffness;
      if (k === undefined) {
        need(ln.EA > 0, `lines[${i}] needs EA (or stiffness, or curve)`);
        const total = ln.length ?? lRef + (ln.inboardLength ?? 0);
        let flex = 0, main = total;
        if (ln.tail) {
          need(ln.tail.EA > 0 && ln.tail.L > 0 && ln.tail.L < total, `lines[${i}].tail needs EA > 0 and 0 < L < line length`);
          flex += ln.tail.L / ln.tail.EA; main -= ln.tail.L;
        }
        flex += main / ln.EA;
        k = 1 / flex;
      }
      need(k > 0, `lines[${i}].stiffness must be > 0`);
      rec.k = k;
      rec.L0 = lRef - pre / k;
      rec.tension = (len) => (len > rec.L0 ? k * (len - rec.L0) : 0);
    }
    return rec;
  });
  const ys = side === 'starboard' ? -v.beam / 2 : v.beam / 2;
  const pb = v.parallelBody ?? [-0.35 * v.LBP, 0.35 * v.LBP];
  need(pb.length === 2 && pb[1] > pb[0], 'vessel.parallelBody must be [xAft, xFwd]');
  const fenders = (o.fenders ?? []).map((f, i) => {
    need(Array.isArray(f.position) && f.position.length === 2 && f.position.every(Number.isFinite), `fenders[${i}].position must be [x, y]`);
    need(f.height > 0, `fenders[${i}].height must be > 0`);
    const c = f.curve;
    need(c && Array.isArray(c.deflection) && Array.isArray(c.reaction) && c.deflection.length === c.reaction.length && c.deflection.length >= 2,
      `fenders[${i}].curve needs deflection[] and reaction[] of equal length`);
    need(c.deflection[0] === 0 && c.reaction[0] === 0, `fenders[${i}].curve must start at (0, 0)`);
    for (let k = 1; k < c.deflection.length; k++) need(c.deflection[k] > c.deflection[k - 1], `fenders[${i}].curve.deflection must be ascending`);
    const n = f.normal ?? [0, 1];
    const nl = hypot(n[0], n[1]);
    need(nl > 0, `fenders[${i}].normal must be non-zero`);
    return { name: f.name ?? `fender ${i + 1}`, p: f.position, n: [n[0] / nl, n[1] / nl], h: f.height, curve: c, rated: f.ratedReaction ?? c.reaction[c.reaction.length - 1] };
  });
  return { v, side, ys, pb, ref, cond, pose0, lines, fenders };
}

function extrap(xs, ys, x) {
  const n = xs.length;
  if (x >= xs[n - 1]) return ys[n - 1] + (ys[n - 1] - ys[n - 2]) * (x - xs[n - 1]) / (xs[n - 1] - xs[n - 2]);
  return interp(xs, ys, x);
}

function lineGeom(rec, pose, cond) {
  const f = rec.def.fairlead, b = rec.def.bollard;
  const c = Math.cos(pose[2]), s = Math.sin(pose[2]);
  const rx = c * f[0] - s * f[1], ry = s * f[0] + c * f[1];
  const fg = [pose[0] + rx, pose[1] + ry, cond.waterLevel - cond.draft + f[2]];
  const d = [b[0] - fg[0], b[1] - fg[1], b[2] - fg[2]];
  return { len: hypot(d[0], d[1], d[2]), d, rx, ry, fg };
}

function evaluate(m, pose, env, detail) {
  const c = Math.cos(pose[2]), s = Math.sin(pose[2]);
  let Fx = 0, Fy = 0, Mz = 0;
  // environment (ship axes -> global); the angle of attack follows the ship heading
  const envOut = [];
  for (const e of env) {
    const theta = e.fromDeg !== undefined ? e.fromDeg - pose[2] / DEG : e.thetaDeg + (m.pose0[2] - pose[2]) / DEG;
    const f = e.fixed ?? environmentalForce({ kind: e.kind, speed: e.speed, thetaDeg: theta, vessel: m.v, curves: e.curves, type: m.v.type, rho: e.rho, shallowFactor: e.shallowFactor });
    Fx += c * f.Fx - s * f.Fy; Fy += s * f.Fx + c * f.Fy; Mz += f.Mz;
    if (detail) envOut.push({ kind: e.kind, thetaDeg: theta, Fx: f.Fx, Fy: f.Fy, Mz: f.Mz, generic: !!f.generic });
  }
  const lines = [];
  for (const ln of m.lines) {
    const g = lineGeom(ln, pose, m.cond);
    const T = ln.tension(g.len);
    if (T > 0) {
      const fx = T * g.d[0] / g.len, fy = T * g.d[1] / g.len;
      Fx += fx; Fy += fy; Mz += g.rx * fy - g.ry * fx;
    }
    if (detail) {
      const hor = hypot(g.d[0], g.d[1]);
      lines.push({
        name: ln.name, type: ln.type, tension: T, slack: !(T > 0), length: g.len, unstretchedLength: ln.L0, stiffness: ln.k,
        utilisation: T / ln.MBL, limitFraction: ln.limit, uc: T / (ln.limit * ln.MBL), ok: T <= ln.limit * ln.MBL,
        verticalAngle: Math.atan2(g.d[2], hor), horizontalAngle: Math.atan2(g.d[1], g.d[0]),
        horizontalTension: g.len > 0 ? T * hor / g.len : 0, verticalForce: g.len > 0 ? T * g.d[2] / g.len : 0,
        bollardId: ln.bollardId, forceOnBollard: g.len > 0 ? [-T * g.d[0] / g.len, -T * g.d[1] / g.len, -T * g.d[2] / g.len] : [0, 0, 0],
      });
    }
  }
  const fenders = [];
  for (const f of m.fenders) {
    // hull side: points c + s e1 + ys e2 with e1 = (cos, sin), e2 = (-sin, cos); intersect with the fender axis
    const px = f.p[0] - pose[0], py = f.p[1] - pose[1];
    const ne2 = -f.n[0] * s + f.n[1] * c;
    let defl = 0, R = 0, sAlong = NaN;
    if (Math.abs(ne2) > 1e-9) {
      const t = (m.ys - (-px * s + py * c)) / ne2;
      sAlong = px * c + py * s + t * (f.n[0] * c + f.n[1] * s);
      if (sAlong >= m.pb[0] && sAlong <= m.pb[1]) {
        defl = f.h - t;
        if (defl > 0) {
          R = extrap(f.curve.deflection, f.curve.reaction, defl);
          if (R < 0) R = 0;
          const cx = px + t * f.n[0], cy = py + t * f.n[1]; // contact point relative to the ship origin
          Fx += R * f.n[0]; Fy += R * f.n[1]; Mz += cx * R * f.n[1] - cy * R * f.n[0];
        } else defl = 0;
      }
    }
    if (detail) {
      const dmax = f.curve.deflection[f.curve.deflection.length - 1];
      fenders.push({ name: f.name, deflection: defl, reaction: R, contactX: sAlong, ratedReaction: f.rated, uc: R / f.rated, overCompressed: defl > dmax, ok: R <= f.rated && defl <= dmax });
    }
  }
  return { F: [Fx, Fy, Mz], lines, fenders, env: envOut };
}

/**
 * Three-degree-of-freedom (surge, sway, yaw) static equilibrium of a ship moored alongside a berth under wind,
 * current and optional additional (e.g. passing-ship) loads, with tension-only mooring lines and
 * compression-only fenders. Damped Newton iteration with a finite-difference Jacobian.
 *
 * Lines: { name, type, fairlead: [x, y, z above keel] (ship frame), bollard: [x, y, z] (global), MBL [N],
 *   material: 'steel' | 'synthetic' | 'polyamide', pretension [N],
 *   EA [N] with optional length [m] (total line length incl. the part on board; default: geometric
 *   fairlead-bollard distance + inboardLength) and tail: { EA, L } (series tail), or stiffness [N/m] directly,
 *   or curve: { strain: [], loadFraction: [] } (load / MBL against strain, through the origin, increasing),
 *   limitFraction (override of the load limit), bollardId }.
 * The unstretched length is set so that the line carries its pretension at pose0 in the reference condition
 * (water level and draft); a different `condition` (tide, draft) changes the fairlead elevation and thereby the
 * line lengths, tensions and vertical angles.
 * Fenders: { name, position: [x, y] of the fender base on the berth (global), normal: [nx, ny] (default
 *   [0, 1]), height (undeflected projection) [m], curve: { deflection: [] [m], reaction: [] [N] },
 *   ratedReaction [N] }. The hull is the straight side y_ship = -beam/2 (starboard side to) or +beam/2 over the
 *   parallel body; fender friction is neglected.
 * Environment: wind / current as { speed [m/s], fromDeg (global direction the flow comes from, counter-clockwise
 *   from the global x axis) or thetaDeg (angle of attack at pose0), curves, rho, shallowFactor }.
 * Passing-ship or other loads: extra = { surge, sway, yaw } in ship axes [N, N, N m] (e.g. peak values from a
 *   Flory (2002)-type estimate or a dedicated passing-ship calculation; no such regression is built in).
 *
 * Reference: OCIMF MEG (static mooring analysis method, standard criteria, line limits); BS 6349-4 (fenders).
 * Validity: static; rigid ship, level trim; no wave-induced motion or dynamic amplification; fenders act along
 * their own axis without friction.
 *
 * @param {object} o
 * @param {object} o.vessel { LBP, beam, draft [m], AT, AL [m^2], type, side: 'starboard' | 'port' (side
 *   alongside), parallelBody: [xAft, xFwd] [m] (default +-0.35 LBP) }
 * @param {{x:number,y:number,headingDeg:number}} o.pose0 initial position of the ship origin [m] and heading
 *   (bow direction, counter-clockwise from global x) [deg]
 * @param {Array} [o.lines] @param {Array} [o.fenders]
 * @param {object} [o.wind] @param {object} [o.current] @param {object} [o.extra]
 * @param {{waterLevel:number,draft:number}} [o.reference] condition at which the pretension is set
 * @param {{waterLevel:number,draft:number}} [o.condition] condition analysed (default = reference)
 * @param {object} [o.limits] line limit fractions by material (default MEG_LINE_LIMITS)
 * @param {number[]} [o.manifold] [x, y] of the manifold in the ship frame [m] (default origin)
 * @param {{surge:number,sway:number}} [o.armEnvelope] allowable +- excursion of the manifold along / normal to
 *   the berth [m]
 * @param {Array<{id:string,SWL:number}>} [o.bollards] safe working loads by bollardId [N]
 * @param {number} [o.tol=1e-9] @param {number} [o.maxIter=100]
 * @returns {object} pose { x, y, headingDeg }; excursion { surge, sway (ship origin, along / normal to the
 *   initial heading) [m], yawDeg, manifoldX, manifoldY (global) [m] }; residual [Fx, Fy, Mz] (global, [N, N,
 *   N m]); environment (forces in ship axes); lines[] (tension [N], utilisation, limitFraction, uc, ok,
 *   verticalAngle, horizontalAngle [rad], horizontalTension, verticalForce [N], forceOnBollard [N]);
 *   bollards[] (id, load [N], horizontal [N], vertical (uplift) [N], SWL, uc); fenders[] (deflection [m],
 *   reaction [N], uc, overCompressed); checks { maxLineUc, maxFenderUc, maxBollardUc, armOk };
 *   converged, iterations
 */
export function berthEquilibrium(o) {
  const m = prepare(o);
  const env = [];
  for (const kind of ['wind', 'current']) {
    const e = o[kind];
    if (!e) continue;
    need(Number.isFinite(e.speed) && e.speed >= 0, `${kind}.speed must be >= 0`);
    need(Number.isFinite(e.fromDeg) || Number.isFinite(e.thetaDeg), `${kind} needs fromDeg (global) or thetaDeg (relative to the bow)`);
    env.push({ kind, ...e });
  }
  if (o.extra) {
    const x = o.extra;
    env.push({ kind: 'extra', thetaDeg: 0, fixed: { Fx: x.surge ?? 0, Fy: x.sway ?? 0, Mz: x.yaw ?? 0, generic: false } });
  }
  const tol = o.tol ?? 1e-9, maxIter = o.maxIter ?? 100;
  const q = m.pose0.slice();
  // force scale: applied loads at the start, pretensions, or a nominal 1 kN
  let Fs = 1e3;
  for (const ln of m.lines) Fs = Math.max(Fs, ln.pre);
  {
    const lone = { ...m, lines: [], fenders: [] };
    const e0 = evaluate(lone, q, env, false).F;
    Fs = Math.max(Fs, hypot(e0[0], e0[1]), Math.abs(e0[2]) / m.v.LBP);
  }
  const sc = [Fs, Fs, Fs * m.v.LBP];
  const merit = (F) => hypot(F[0] / sc[0], F[1] / sc[1], F[2] / sc[2]);
  let F = evaluate(m, q, env, false).F, f = merit(F), it = 0, lam = 0;
  const h = [1e-5, 1e-5, 1e-7];
  let drift = 0;
  for (; it < maxIter && f > tol; it++) {
    const J = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let b = 0; b < 3; b++) {
      const qp = q.slice(), qm = q.slice();
      qp[b] += h[b]; qm[b] -= h[b];
      const Fp = evaluate(m, qp, env, false).F, Fm = evaluate(m, qm, env, false).F;
      for (let a = 0; a < 3; a++) J[a][b] = (Fp[a] - Fm[a]) / (2 * h[b]);
    }
    let accepted = false;
    for (let tries = 0; tries < 12 && !accepted; tries++) {
      // Levenberg-type regularisation keeps the step defined when lines are slack or a fender is out of contact
      const A = J.map((row, a) => row.map((val, b) => (a === b ? val - lam * (a < 2 ? Fs : Fs * m.v.LBP * m.v.LBP) : val)));
      let dq;
      try { dq = solveDense(A, F.map((x) => -x)); } catch (e) { dq = null; }
      if (dq && dq.every(Number.isFinite)) {
        let cap = 1;
        const lim = [0.5, 0.5, 0.02];
        for (let b = 0; b < 3; b++) if (Math.abs(dq[b]) * cap > lim[b]) cap = lim[b] / Math.abs(dq[b]);
        let alpha = cap;
        for (let ls = 0; ls < 12; ls++) {
          const qn = [q[0] + alpha * dq[0], q[1] + alpha * dq[1], q[2] + alpha * dq[2]];
          const Fn = evaluate(m, qn, env, false).F, fn = merit(Fn);
          if (fn < f) { q[0] = qn[0]; q[1] = qn[1]; q[2] = qn[2]; F = Fn; f = fn; accepted = true; break; }
          alpha *= 0.5;
        }
      }
      if (accepted) lam *= 0.3; else lam = lam === 0 ? 1e-3 : lam * 10;
    }
    if (!accepted) {
      // No Newton step reduces the unbalanced load: the ship is drifting free of every line and fender (zero
      // stiffness). Move it a short distance along the unbalanced load until something engages.
      let step = 0.25;
      for (let k = 0; k < 8 && !accepted; k++, step *= 0.25) {
        const qn = [q[0] + step * F[0] / (sc[0] * f), q[1] + step * F[1] / (sc[1] * f), q[2] + (step / m.v.LBP) * F[2] / (sc[2] * f)];
        const Fn = evaluate(m, qn, env, false).F, fn = merit(Fn);
        if (fn <= f * (1 + 1e-12)) { q[0] = qn[0]; q[1] = qn[1]; q[2] = qn[2]; F = Fn; f = fn; accepted = true; }
      }
      if (!accepted || ++drift > 400) break;
      lam = 0;
    }
  }
  const converged = f <= Math.max(tol, 1e-7);
  if (!converged) {
    throw new Error('Berth mooring equilibrium did not converge: the vessel is not restrained in every direction (check for slack lines, missing springs or fenders not in contact)');
  }
  const r = evaluate(m, q, env, true);
  const c0 = Math.cos(m.pose0[2]), s0 = Math.sin(m.pose0[2]);
  const dx = q[0] - m.pose0[0], dy = q[1] - m.pose0[1];
  const man = o.manifold ?? [0, 0];
  const mp = (p) => [p[0] + Math.cos(p[2]) * man[0] - Math.sin(p[2]) * man[1], p[1] + Math.sin(p[2]) * man[0] + Math.cos(p[2]) * man[1]];
  const m1 = mp(q), m0 = mp(m.pose0);
  const excursion = { surge: dx * c0 + dy * s0, sway: -dx * s0 + dy * c0, yawDeg: (q[2] - m.pose0[2]) / DEG, manifoldX: m1[0] - m0[0], manifoldY: m1[1] - m0[1] };
  // bollard / hook loads: vector sum of the lines on the same bollard
  const byId = new Map();
  for (const ln of r.lines) {
    const b = byId.get(ln.bollardId) ?? [0, 0, 0];
    b[0] += ln.forceOnBollard[0]; b[1] += ln.forceOnBollard[1]; b[2] += ln.forceOnBollard[2];
    byId.set(ln.bollardId, b);
  }
  const swl = new Map((o.bollards ?? []).map((b) => [b.id, b.SWL]));
  const bollards = [];
  let maxB = 0;
  for (const [id, b] of byId) {
    const load = hypot(b[0], b[1], b[2]);
    const S = swl.get(id);
    const uc = S > 0 ? load / S : null;
    if (uc !== null) maxB = Math.max(maxB, uc);
    bollards.push({ id, load, horizontal: hypot(b[0], b[1]), vertical: b[2], force: b, SWL: S ?? null, uc });
  }
  let maxL = 0, maxF = 0;
  for (const ln of r.lines) maxL = Math.max(maxL, ln.uc);
  for (const fd of r.fenders) maxF = Math.max(maxF, fd.uc);
  const arm = o.armEnvelope;
  const armOk = arm ? Math.abs(excursion.manifoldX) <= arm.surge && Math.abs(excursion.manifoldY) <= arm.sway : null;
  return {
    pose: { x: q[0], y: q[1], headingDeg: q[2] / DEG }, excursion, residual: r.F, environment: r.env,
    lines: r.lines, bollards, fenders: r.fenders,
    checks: { maxLineUc: maxL, maxFenderUc: maxF, maxBollardUc: swl.size ? maxB : null, armOk },
    converged, iterations: it,
  };
}

/**
 * Run a list of environmental cases (e.g. ocimfStandardCases) through berthEquilibrium and collect the envelope.
 * @param {object} model arguments of berthEquilibrium without wind / current
 * @param {Array<{name:string,wind?:object,current?:object,extra?:object}>} cases
 * @returns {{cases:Array<{name:string,maxLineUc:number,maxLine:string,maxFenderUc:number,surge:number,sway:number,
 *   converged:boolean,error?:string}>,governing:{line:object|null,fender:object|null},maxLineUc:number,
 *   maxFenderUc:number,maxSurge:number,maxSway:number}} per-case summary and envelope (excursions in m)
 */
export function runBerthCases(model, cases) {
  need(Array.isArray(cases) && cases.length > 0, 'cases must be a non-empty array');
  const out = [];
  let gL = null, gF = null, mL = 0, mF = 0, mS = 0, mW = 0;
  for (const cs of cases) {
    try {
      const r = berthEquilibrium({ ...model, wind: cs.wind, current: cs.current, extra: cs.extra ?? model.extra });
      let worst = null;
      for (const ln of r.lines) if (!worst || ln.uc > worst.uc) worst = ln;
      const rec = {
        name: cs.name, maxLineUc: r.checks.maxLineUc, maxLine: worst ? worst.name : '', maxFenderUc: r.checks.maxFenderUc,
        surge: r.excursion.manifoldX, sway: r.excursion.manifoldY, converged: true,
      };
      out.push(rec);
      if (rec.maxLineUc >= mL) { mL = rec.maxLineUc; gL = rec; }
      if (rec.maxFenderUc >= mF) { mF = rec.maxFenderUc; gF = rec; }
      mS = Math.max(mS, Math.abs(rec.surge)); mW = Math.max(mW, Math.abs(rec.sway));
    } catch (e) {
      if (e instanceof RangeError) throw e;
      out.push({ name: cs.name, maxLineUc: NaN, maxLine: '', maxFenderUc: NaN, surge: NaN, sway: NaN, converged: false, error: e.message });
    }
  }
  return { cases: out, governing: { line: gL, fender: gF }, maxLineUc: mL, maxFenderUc: mF, maxSurge: mS, maxSway: mW };
}
