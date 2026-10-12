/**
 * Three-dimensional lumped-mass dynamic model of a mooring line / cable (the formulation used by MoorDyn-type
 * codes), plus a simple coupled floater model for the line-break transient.
 *
 * The line is divided into N segments; nodes 0 (anchor, fixed) ... N (fairlead, prescribed motion). Each segment
 * is an axial spring-damper (EA, internal damping BA); mass, net weight, added mass, Morison drag and seabed
 * contact are lumped at the nodes (half of each adjacent segment). Optional bending stiffness acts through the
 * discrete bending energy (EI / l)(1 - u_a . u_b) at interior nodes.
 *
 * References:
 *  - Hall, M. & Goupee, A. (2015) "Validation of a lumped-mass mooring line model with DeepCwind
 *    semisubmersible model test data", Ocean Engineering 104, 590-603 (model equations: node mass and added-mass
 *    matrix, segment tension and damping, normal / tangential drag, seabed spring-damper).
 *  - Morison et al. (1950); DNV-RP-C205 section 6 (relative-velocity Morison load on slender members).
 *  - Airy wave kinematics: DNV-RP-C205 section 3.
 *
 * Global frame: z up. When waves or a current profile are used, z = waterLevel (default 0) is the still-water
 * level. All quantities SI.
 */

import { solveMultiSegment } from './multiseg.js';
import { equilibrium, lineForces } from './mooringsystem.js';

const hypot = Math.hypot;

function need(cond, msg) {
  if (!cond) throw new RangeError(msg);
}

function vec3(v, name) {
  need(Array.isArray(v) && v.length === 3 && v.every(Number.isFinite), `${name} must be [x, y, z] of finite numbers`);
}

/**
 * Airy (linear) wave number from the dispersion relation w^2 = g k tanh(k h), by Newton iteration.
 * Reference: DNV-RP-C205 section 3.2. Valid for h > 0, T > 0.
 * @param {number} omega angular frequency [rad/s]
 * @param {number} depth water depth [m]
 * @param {number} [g=9.80665] [m/s^2]
 * @returns {number} wave number k [rad/m]
 */
export function airyWaveNumber(omega, depth, g = 9.80665) {
  need(omega > 0, 'omega must be > 0');
  need(depth > 0, 'depth must be > 0');
  const k0 = omega * omega / g;
  let k = k0 / Math.sqrt(Math.tanh(k0 * depth));
  for (let i = 0; i < 50; i++) {
    const th = Math.tanh(k * depth);
    const f = g * k * th - omega * omega;
    const df = g * th + g * k * depth * (1 - th * th);
    const dk = f / df;
    k -= dk;
    if (Math.abs(dk) < 1e-14 * k) break;
  }
  return k;
}

/**
 * Build a lumped-mass line model.
 *
 * @param {object} o
 * @param {number[]} o.anchor [x, y, z] of the fixed lower end [m]
 * @param {number[]} o.fairlead [x, y, z] of the upper end at rest [m]
 * @param {Array<object>} o.sections line make-up from anchor to fairlead; each section:
 *   L unstretched length [m]; EA [N]; mass in air [kg/m]; d volume-equivalent diameter [m] (buoyancy, added mass
 *   and drag area); nSeg number of segments; Cdn, Cdt drag coefficients [-] (normal on d l, tangential on
 *   pi d l; default 1.2, 0.05); Can, Cat added-mass coefficients [-] (default 1.0, 0.0); BA internal damping
 *   [N s] (segment damping force BA (dl/dt)/l0) or zeta, the damping ratio of a segment's own axial mode [-]
 *   (BA = zeta l0 sqrt(EA mass); default zeta = 0); EI bending stiffness [N m^2] (default 0)
 * @param {Array<{junction:number,mass?:number,volume?:number,CdA?:number,Ca?:number}>} [o.points] clump weights /
 *   buoys at section junctions (junction k between sections k-1 and k): mass [kg], displaced volume [m^3],
 *   drag area CdA [m^2], added-mass coefficient Ca [-] on the displaced volume
 * @param {number} [o.rho=1025] water density [kg/m^3] @param {number} [o.g=9.80665] [m/s^2]
 * @param {object|false} [o.seabed] { z (default anchor z) [m], k normal stiffness per unit contact area d l
 *   [Pa/m] (default 3e6), c normal damping per unit area [Pa s/m] or zeta, the damping ratio of a node on the
 *   contact spring [-] (default zeta = 0.5), mu Coulomb friction coefficient [-] (default 0), vc friction
 *   regularisation speed [m/s] (default 0.01; friction = mu N v/max(|v|, vc), limited so that it cannot
 *   reverse a node within one step) }; false disables seabed contact
 * @param {object} [o.current] steady current: { speed [m/s], direction [rad] } uniform, or
 *   { z: number[] ascending [m], speed: number[] [m/s], direction [rad] } linear in z, clamped at the ends
 * @param {object} [o.wave] Airy wave: { H height [m], T period [s], direction [rad], depth [m], phase [rad] }
 *   (kinematics evaluated at min(z, waterLevel), no stretching)
 * @param {number} [o.waterLevel=0] still-water level z [m]
 * @param {boolean} [o.allowCompression=false] let segments carry compression (for bending-stiff cables)
 * @returns {object} model (opaque; pass to relaxLine / simulateLine). Fields of interest: N, s (node arc
 *   lengths [m]), r (node positions, Float64Array 3(N+1)) [m], dtCritical [s]
 */
export function buildLumpedLine(o) {
  vec3(o.anchor, 'anchor'); vec3(o.fairlead, 'fairlead');
  need(Array.isArray(o.sections) && o.sections.length > 0, 'sections must be a non-empty array');
  const rho = o.rho ?? 1025, g = o.g ?? 9.80665;
  need(rho >= 0, 'rho must be >= 0');
  let N = 0;
  o.sections.forEach((sc, i) => {
    for (const k of ['L', 'EA', 'mass', 'd']) need(Number.isFinite(sc[k]) && sc[k] > 0, `sections[${i}].${k} must be > 0`);
    need(Number.isInteger(sc.nSeg) && sc.nSeg >= 1, `sections[${i}].nSeg must be an integer >= 1`);
    N += sc.nSeg;
  });
  need(N >= 2, 'the line needs at least 2 segments in total');
  const nn = N + 1;
  const m = {
    N, rho, g, anchor: o.anchor.slice(), fairlead0: o.fairlead.slice(),
    l0: new Float64Array(N), EA: new Float64Array(N), BA: new Float64Array(N), EI: new Float64Array(N),
    s: new Float64Array(nn), sec: new Int32Array(N),
    mass: new Float64Array(nn), W: new Float64Array(nn), an: new Float64Array(nn), at: new Float64Array(nn),
    kdn: new Float64Array(nn), kdt: new Float64Array(nn), kpt: new Float64Array(nn),
    fin: new Float64Array(nn), fit: new Float64Array(nn), Ab: new Float64Array(nn), EIb: new Float64Array(nn),
    r: new Float64Array(3 * nn), v: new Float64Array(3 * nn), a: new Float64Array(3 * nn), F: new Float64Array(3 * nn),
    T: new Float64Array(N), len: new Float64Array(N), ux: new Float64Array(N), uy: new Float64Array(N), uz: new Float64Array(N),
    Fbed: new Float64Array(nn),
    allowCompression: !!o.allowCompression, waterLevel: o.waterLevel ?? 0,
    frictionOn: true, wavesOn: true, dragOnMotion: true, time: 0, hNow: 1,
  };
  let k = 0, sAcc = 0, cmax = 0, dtc = Infinity;
  const junctionNode = [0];
  const segDefs = [];
  o.sections.forEach((sc, i) => {
    const l0 = sc.L / sc.nSeg;
    const Cdn = sc.Cdn ?? 1.2, Cdt = sc.Cdt ?? 0.05, Can = sc.Can ?? 1.0, Cat = sc.Cat ?? 0.0;
    need(Cdn >= 0 && Cdt >= 0 && Can >= 0 && Cat >= 0, `sections[${i}] hydrodynamic coefficients must be >= 0`);
    const BA = sc.BA ?? (sc.zeta ?? 0) * l0 * Math.sqrt(sc.EA * sc.mass);
    need(BA >= 0, `sections[${i}].BA must be >= 0`);
    const area = Math.PI * sc.d * sc.d / 4;
    const c = Math.sqrt(sc.EA / sc.mass);
    cmax = Math.max(cmax, c);
    // highest axial mode of a segment chain: w = 2 c / l0; damping ratio of that mode: BA w / (2 EA)
    const zh = BA * c / (sc.EA * l0);
    // bending waves: highest discrete mode w = 4 sqrt(EI / m) / l0^2
    const dta = (l0 / c) / (1 + 2 * zh);
    const dtb = sc.EI > 0 ? l0 * l0 / (2 * Math.sqrt(sc.EI / sc.mass)) : Infinity;
    dtc = Math.min(dtc, 1 / Math.sqrt(1 / (dta * dta) + 1 / (dtb * dtb)));
    for (let j = 0; j < sc.nSeg; j++, k++) {
      m.l0[k] = l0; m.EA[k] = sc.EA; m.BA[k] = BA; m.EI[k] = sc.EI ?? 0; m.sec[k] = i;
      sAcc += l0; m.s[k + 1] = sAcc;
      for (const nd of [k, k + 1]) {
        m.mass[nd] += 0.5 * sc.mass * l0;
        m.W[nd] += 0.5 * (sc.mass - rho * area) * g * l0;
        m.an[nd] += 0.5 * rho * area * l0 * Can;
        m.at[nd] += 0.5 * rho * area * l0 * Cat;
        m.fin[nd] += 0.5 * rho * area * l0 * (1 + Can);
        m.fit[nd] += 0.5 * rho * area * l0 * (1 + Cat);
        m.kdn[nd] += 0.25 * rho * Cdn * sc.d * l0;
        m.kdt[nd] += 0.25 * rho * Cdt * Math.PI * sc.d * l0;
        m.Ab[nd] += 0.5 * sc.d * l0;
      }
    }
    junctionNode.push(k);
    segDefs.push({ L: sc.L, w: (sc.mass - rho * area) * g, EA: sc.EA });
  });
  for (let i = 1; i < N; i++) {
    const EIi = 0.5 * (m.EI[i - 1] + m.EI[i]);
    m.EIb[i] = EIi / (0.5 * (m.l0[i - 1] + m.l0[i]));
  }
  m.hasBending = m.EIb.some((x) => x > 0);
  const pts = [];
  for (const p of o.points ?? []) {
    need(Number.isInteger(p.junction) && p.junction >= 1 && p.junction <= o.sections.length - 1,
      'points[].junction must be an integer between 1 and sections.length - 1');
    const nd = junctionNode[p.junction];
    const pm = p.mass ?? 0, pv = p.volume ?? 0;
    need(pm >= 0 && pv >= 0, 'points[] mass and volume must be >= 0');
    m.mass[nd] += pm; m.W[nd] += (pm - rho * pv) * g;
    const ca = rho * pv * (p.Ca ?? 0);
    m.an[nd] += ca; m.at[nd] += ca; m.fin[nd] += rho * pv + ca; m.fit[nd] += rho * pv + ca;
    m.kpt[nd] += 0.5 * rho * (p.CdA ?? 0);
    pts.push({ junction: p.junction, weight: (pm - rho * pv) * g });
  }
  m.dtCritical = dtc;
  m.Ltot = sAcc;

  // seabed
  if (o.seabed === false) m.bed = null;
  else {
    const b = o.seabed ?? {};
    m.bed = { z: b.z ?? o.anchor[2], k: b.k ?? 3e6, zeta: b.zeta ?? 0.5, mu: b.mu ?? 0, vc: b.vc ?? 0.01 };
    need(m.bed.k > 0 && m.bed.zeta >= 0 && m.bed.mu >= 0 && m.bed.vc > 0, 'seabed needs k > 0, zeta >= 0, mu >= 0, vc > 0');
    if (b.c !== undefined) need(b.c >= 0, 'seabed.c must be >= 0');
    // node contact damper [N s/m] and the explicit stability limit of the contact spring-damper
    m.cbed = new Float64Array(nn);
    for (let i = 0; i < nn; i++) {
      const M = m.mass[i] + m.an[i], kk = m.bed.k * m.Ab[i];
      const om = Math.sqrt(kk / M);
      m.cbed[i] = b.c !== undefined ? b.c * m.Ab[i] : 2 * m.bed.zeta * Math.sqrt(kk * M);
      const z = m.cbed[i] / (2 * M * om);
      // stable step of the damped contact oscillator is (2 / om)(sqrt(1 + z^2) - z); dtCritical is half the limit
      m.dtCritical = Math.min(m.dtCritical, (Math.sqrt(1 + z * z) - z) / om);
    }
  }
  // current
  if (o.current) {
    const c = o.current;
    need(Number.isFinite(c.direction ?? 0), 'current.direction must be finite');
    if (Array.isArray(c.z)) {
      need(Array.isArray(c.speed) && c.speed.length === c.z.length && c.z.length >= 1, 'current.z and current.speed must have equal length');
      for (let i = 1; i < c.z.length; i++) need(c.z[i] > c.z[i - 1], 'current.z must be ascending');
      m.cur = { z: Float64Array.from(c.z), u: Float64Array.from(c.speed), cx: Math.cos(c.direction ?? 0), cy: Math.sin(c.direction ?? 0) };
    } else {
      need(Number.isFinite(c.speed), 'current.speed must be finite');
      m.cur = { z: Float64Array.from([0]), u: Float64Array.from([c.speed]), cx: Math.cos(c.direction ?? 0), cy: Math.sin(c.direction ?? 0) };
    }
  } else m.cur = null;
  // wave
  if (o.wave) {
    const w = o.wave;
    need(w.H >= 0 && w.T > 0 && w.depth > 0, 'wave needs H >= 0, T > 0, depth > 0');
    const om = 2 * Math.PI / w.T;
    const kk = airyWaveNumber(om, w.depth, g);
    m.wave = { A: w.H / 2, om, k: kk, h: w.depth, cx: Math.cos(w.direction ?? 0), cy: Math.sin(w.direction ?? 0), ph: w.phase ?? 0 };
  } else m.wave = null;

  // initial shape: static catenary in the vertical plane through the two ends
  const dx = o.fairlead[0] - o.anchor[0], dy = o.fairlead[1] - o.anchor[1];
  const X = hypot(dx, dy), Z = o.fairlead[2] - o.anchor[2];
  const ex = X > 0 ? dx / X : 1, ey = X > 0 ? dy / X : 0;
  const onBed = !!m.bed && Math.abs(m.bed.z - o.anchor[2]) < 1e-6 * Math.max(1, sAcc);
  m.static = { segments: segDefs, points: pts, seabed: onBed, ex, ey };
  let prof = null;
  try { prof = solveMultiSegment({ X, Z, segments: segDefs, points: pts, seabed: onBed && Z >= 0, nProfile: 200, stiffness: false }).profile; } catch (e) { prof = null; }
  let j = 0;
  for (let i = 0; i <= N; i++) {
    let px, pz;
    if (prof) {
      const si = Math.min(m.s[i], prof.s[prof.s.length - 1]);
      while (j < prof.s.length - 2 && prof.s[j + 1] < si) j++;
      const ds = prof.s[j + 1] - prof.s[j];
      const f = ds > 0 ? Math.min(Math.max((si - prof.s[j]) / ds, 0), 1) : 0;
      px = prof.x[j] + f * (prof.x[j + 1] - prof.x[j]); pz = prof.z[j] + f * (prof.z[j + 1] - prof.z[j]);
    } else { px = X * m.s[i] / sAcc; pz = Z * m.s[i] / sAcc; }
    m.r[3 * i] = o.anchor[0] + ex * px; m.r[3 * i + 1] = o.anchor[1] + ey * px; m.r[3 * i + 2] = o.anchor[2] + pz;
  }
  m.r[0] = o.anchor[0]; m.r[1] = o.anchor[1]; m.r[2] = o.anchor[2];
  m.r[3 * N] = o.fairlead[0]; m.r[3 * N + 1] = o.fairlead[1]; m.r[3 * N + 2] = o.fairlead[2];
  return m;
}

// Node accelerations for the state (r, v) at time t. Fills m.T (segment tensions), m.F (net node forces, all
// nodes) and a (interior nodes).
function accel(m, t, r, v, a) {
  const { N, l0, EA, BA, T, len, ux, uy, uz, F } = m;
  const comp = m.allowCompression;
  for (let s = 0; s < N; s++) {
    const i = 3 * s;
    const dx = r[i + 3] - r[i], dy = r[i + 4] - r[i + 1], dz = r[i + 5] - r[i + 2];
    const l = Math.sqrt(dx * dx + dy * dy + dz * dz);
    const x = dx / l, y = dy / l, z = dz / l;
    let Ts = EA[s] * (l - l0[s]) / l0[s];
    if (Ts > 0 || comp) {
      Ts += BA[s] * ((v[i + 3] - v[i]) * x + (v[i + 4] - v[i + 1]) * y + (v[i + 5] - v[i + 2]) * z) / l0[s];
      if (Ts < 0 && !comp) Ts = 0;
    } else Ts = 0;
    T[s] = Ts; len[s] = l; ux[s] = x; uy[s] = y; uz[s] = z;
  }
  F.fill(0);
  for (let s = 0; s < N; s++) {
    const i = 3 * s, fx = T[s] * ux[s], fy = T[s] * uy[s], fz = T[s] * uz[s];
    F[i] += fx; F[i + 1] += fy; F[i + 2] += fz;
    F[i + 3] -= fx; F[i + 4] -= fy; F[i + 5] -= fz;
  }
  if (m.hasBending) {
    for (let i = 1; i < N; i++) {
      const kb = m.EIb[i];
      if (kb === 0) continue;
      const sa = i - 1, sb = i;
      const c = ux[sa] * ux[sb] + uy[sa] * uy[sb] + uz[sa] * uz[sb];
      const ga = kb / len[sa], gb = kb / len[sb];
      const gax = ga * (ux[sb] - c * ux[sa]), gay = ga * (uy[sb] - c * uy[sa]), gaz = ga * (uz[sb] - c * uz[sa]);
      const gbx = gb * (ux[sa] - c * ux[sb]), gby = gb * (uy[sa] - c * uy[sb]), gbz = gb * (uz[sa] - c * uz[sb]);
      const j = 3 * i;
      F[j + 3] += gbx; F[j + 4] += gby; F[j + 5] += gbz;
      F[j - 3] -= gax; F[j - 2] -= gay; F[j - 1] -= gaz;
      F[j] += gax - gbx; F[j + 1] += gay - gby; F[j + 2] += gaz - gbz;
    }
  }
  const cur = m.cur, wave = m.wavesOn ? m.wave : null, bed = m.bed;
  const wl = m.waterLevel;
  for (let nd = 0; nd <= N; nd++) {
    const j = 3 * nd;
    // unit tangent from the neighbouring nodes
    const ja = nd > 0 ? j - 3 : j, jb = nd < N ? j + 3 : j;
    let qx = r[jb] - r[ja], qy = r[jb + 1] - r[ja + 1], qz = r[jb + 2] - r[ja + 2];
    const ql = Math.sqrt(qx * qx + qy * qy + qz * qz);
    qx /= ql; qy /= ql; qz /= ql;
    let fx = F[j], fy = F[j + 1], fz = F[j + 2] - m.W[nd];
    // fluid kinematics
    let Ux = 0, Uy = 0, Uz = 0, Ax = 0, Ay = 0, Az = 0;
    const zn = r[j + 2];
    if (cur) {
      const cz = cur.z, cu = cur.u, n = cz.length;
      let u;
      if (n === 1 || zn <= cz[0]) u = cu[0];
      else if (zn >= cz[n - 1]) u = cu[n - 1];
      else {
        let k = 0;
        while (zn > cz[k + 1]) k++;
        u = cu[k] + (cu[k + 1] - cu[k]) * (zn - cz[k]) / (cz[k + 1] - cz[k]);
      }
      Ux = u * cur.cx; Uy = u * cur.cy;
    }
    if (wave) {
      const zz = Math.min(zn, wl) - wl;
      const th = wave.k * (r[j] * wave.cx + r[j + 1] * wave.cy) - wave.om * t + wave.ph;
      const kh = wave.k * wave.h, kz = wave.k * (zz + wave.h);
      // ratios cosh(k(z+h))/sinh(kh), sinh(k(z+h))/sinh(kh) written to stay finite in deep water
      const e1 = Math.exp(kz - kh), e2 = Math.exp(-kz - kh), den = 1 - Math.exp(-2 * kh);
      const ch = (e1 + e2) / den, sh = (e1 - e2) / den;
      const ct = Math.cos(th), st = Math.sin(th);
      const uh = wave.A * wave.om * ch * ct, ah = wave.A * wave.om * wave.om * ch * st;
      Ux += uh * wave.cx; Uy += uh * wave.cy; Uz += wave.A * wave.om * sh * st;
      Ax = ah * wave.cx; Ay = ah * wave.cy; Az = -wave.A * wave.om * wave.om * sh * ct;
    }
    // Morison drag on the relative velocity, split into tangential and normal parts
    // (during static relaxation the drag is that of the fluid on the line at rest: see relaxLine)
    const vs = m.dragOnMotion ? 1 : 0;
    const rx = Ux - vs * v[j], ry = Uy - vs * v[j + 1], rz = Uz - vs * v[j + 2];
    const rt = rx * qx + ry * qy + rz * qz;
    const nx = rx - rt * qx, ny = ry - rt * qy, nz = rz - rt * qz;
    const nmag = Math.sqrt(nx * nx + ny * ny + nz * nz);
    const kt = m.kdt[nd] * Math.abs(rt) * rt, kn = m.kdn[nd] * nmag;
    fx += kn * nx + kt * qx; fy += kn * ny + kt * qy; fz += kn * nz + kt * qz;
    if (m.kpt[nd] > 0) {
      const rm = m.kpt[nd] * Math.sqrt(rx * rx + ry * ry + rz * rz);
      fx += rm * rx; fy += rm * ry; fz += rm * rz;
    }
    if (wave) {
      // Froude-Krylov plus added-mass force of the undisturbed fluid acceleration
      const atn = Ax * qx + Ay * qy + Az * qz;
      const fn = m.fin[nd], ft = m.fit[nd];
      fx += fn * (Ax - atn * qx) + ft * atn * qx;
      fy += fn * (Ay - atn * qy) + ft * atn * qy;
      fz += fn * (Az - atn * qz) + ft * atn * qz;
    }
    m.Fbed[nd] = 0;
    if (bed && zn < bed.z) {
      let fn = m.Ab[nd] * bed.k * (bed.z - zn) - m.cbed[nd] * v[j + 2];
      if (fn < 0) fn = 0;
      fz += fn; m.Fbed[nd] = fn;
      if (bed.mu > 0 && m.frictionOn) {
        const vh = Math.sqrt(v[j] * v[j] + v[j + 1] * v[j + 1]);
        if (vh > 0) {
          let ff = bed.mu * fn / Math.max(vh, bed.vc);
          const cap = (m.mass[nd] + m.an[nd]) / m.hNow;
          if (ff > cap) ff = cap;
          fx -= ff * v[j]; fy -= ff * v[j + 1];
        }
      }
    }
    F[j] = fx; F[j + 1] = fy; F[j + 2] = fz;
    if (nd > 0 && nd < N) {
      // (m + an) a + (at - an)(q . a) q = F
      const mn = m.mass[nd] + m.an[nd], mt = m.mass[nd] + m.at[nd];
      const qa = (m.at[nd] - m.an[nd]) * (fx * qx + fy * qy + fz * qz) / mt;
      a[j] = (fx - qa * qx) / mn; a[j + 1] = (fy - qa * qy) / mn; a[j + 2] = (fz - qa * qz) / mn;
    }
  }
}

function endTension(m, nd) {
  const j = 3 * nd;
  return hypot(m.F[j], m.F[j + 1], m.F[j + 2]);
}

/**
 * Static configuration by dynamic relaxation with kinetic damping, starting from the catenary shape: the
 * equations of motion are integrated with the fairlead held still (waves and seabed friction switched off, the
 * steady current kept) and all velocities are reset to zero whenever the kinetic energy passes a peak, until the
 * largest unbalanced node force is below tol times the force scale.
 * The pseudo-motion is a numerical device, so the fluid drag is evaluated for the line at rest in the current
 * (the static drag load) and not on the relative velocity: with the relative velocity a line that has to swing
 * a long way out of its starting plane — a cross current — can only creep there at about the current speed,
 * through tens of seconds of pseudo-time, while kinetic damping alone removes the energy in a few oscillations.
 * The equilibrium found is the same (all velocities are zero there).
 * Reference: Day (1965) / Barnes (1999) dynamic relaxation with kinetic damping.
 * @param {object} model model from buildLumpedLine (updated in place: positions set to the static state)
 * @param {object} [o]
 * @param {number} [o.tol=1e-7] unbalanced-force tolerance relative to the larger of the node weight and the
 *   mean segment tension [-]
 * @param {number} [o.maxIter=2000000] maximum number of pseudo-time steps
 * @param {number} [o.cfl=0.5] pseudo-time step as a fraction of the critical step [-]
 * @returns {{converged:boolean,iterations:number,residual:number,fairleadTension:number,anchorTension:number,
 *   segmentTension:Float64Array,H:number,groundedNodes:number,touchdownS:number}} residual = largest
 *   unbalanced force / force scale [-]; tensions [N]; fairleadTension / anchorTension are the end reactions
 *   (segment tension combined with the end node's own half-segment loads); H = horizontal component of the
 *   fairlead reaction [N]; touchdownS arc length of the touchdown point [m]
 */
export function relaxLine(model, o = {}) {
  const m = model;
  const tol = o.tol ?? 1e-7, maxIter = o.maxIter ?? 2000000;
  const dt = (o.cfl ?? 0.5) * m.dtCritical;
  const { N, r, v, a } = m;
  const wOn = m.wavesOn, fOn = m.frictionOn;
  m.wavesOn = false; m.frictionOn = false;
  const dOn = m.dragOnMotion;
  m.dragOnMotion = false;
  v.fill(0);
  let kePrev = 0, it = 0, resid = Infinity, converged = false, fresh = true;
  let wmax = 0;
  for (let i = 0; i <= N; i++) wmax = Math.max(wmax, Math.abs(m.W[i]));
  for (; it < maxIter; it++) {
    accel(m, 0, r, v, a);
    if (fresh) {
      // velocities are zero: m.F is the static unbalanced force
      let fm = 0, tm = 0;
      for (let i = 1; i < N; i++) fm = Math.max(fm, hypot(m.F[3 * i], m.F[3 * i + 1], m.F[3 * i + 2]));
      for (let s = 0; s < N; s++) tm += m.T[s] / N;
      resid = fm / Math.max(wmax, tm, 1e-9 * m.EA[0]);
      if (resid < tol) { converged = true; break; }
      fresh = false;
    }
    let ke = 0;
    for (let i = 1; i < N; i++) {
      const j = 3 * i;
      v[j] += dt * a[j]; v[j + 1] += dt * a[j + 1]; v[j + 2] += dt * a[j + 2];
      // flush vanishing velocities of grounded nodes (subnormal arithmetic is very slow)
      if (Math.abs(v[j]) < 1e-100) v[j] = 0;
      if (Math.abs(v[j + 1]) < 1e-100) v[j + 1] = 0;
      if (Math.abs(v[j + 2]) < 1e-100) v[j + 2] = 0;
      ke += (m.mass[i] + m.an[i]) * (v[j] * v[j] + v[j + 1] * v[j + 1] + v[j + 2] * v[j + 2]);
    }
    if (ke < kePrev) {
      // kinetic-energy peak passed: step back half a step towards the peak and restart from rest
      for (let i = 1; i < N; i++) { const j = 3 * i; r[j] -= 0.5 * dt * v[j]; r[j + 1] -= 0.5 * dt * v[j + 1]; r[j + 2] -= 0.5 * dt * v[j + 2]; }
      v.fill(0); kePrev = 0; fresh = true;
      continue;
    }
    kePrev = ke;
    for (let i = 1; i < N; i++) { const j = 3 * i; r[j] += dt * v[j]; r[j + 1] += dt * v[j + 1]; r[j + 2] += dt * v[j + 2]; }
  }
  v.fill(0);
  accel(m, 0, r, v, a);
  m.wavesOn = wOn; m.frictionOn = fOn; m.dragOnMotion = dOn;
  const td = touchdown(m);
  m.Tstatic = Float64Array.from(m.T);
  m.staticFairlead = endTension(m, N);
  m.staticAnchor = endTension(m, 0);
  return {
    converged, iterations: it, residual: resid, fairleadTension: m.staticFairlead, anchorTension: m.staticAnchor,
    segmentTension: Float64Array.from(m.T), H: hypot(m.F[3 * N], m.F[3 * N + 1]), groundedNodes: td.count, touchdownS: td.s,
  };
}

function touchdown(m) {
  if (!m.bed) return { count: 0, s: 0 };
  const { N, r, s, l0 } = m;
  // a node counts as grounded while it is within its own static spring deflection (plus 1 mm) of the seabed
  let i = 0, count = 0;
  for (let k = 0; k <= N; k++) if (r[3 * k + 2] <= m.bed.z + 1e-3) count++;
  while (i < N && r[3 * (i + 1) + 2] <= m.bed.z + 1e-3) i++;
  if (i >= N) return { count, s: s[N] };
  if (r[3 * i + 2] > m.bed.z + 1e-3) return { count, s: 0 };
  const z0 = r[3 * i + 2], z1 = r[3 * (i + 1) + 2];
  const f = Math.min(Math.max((m.bed.z + 1e-3 - z0) / (z1 - z0), 0), 1);
  return { count, s: s[i] + f * l0[i] };
}

function makeMotion(spec) {
  if (!spec) return (t, d, w) => { d[0] = d[1] = d[2] = 0; w[0] = w[1] = w[2] = 0; };
  if (typeof spec === 'function') {
    return (t, d, w) => {
      const q = spec(t);
      need(q && Array.isArray(q.r) && Array.isArray(q.v), 'motion function must return { r: [3], v: [3] }');
      d[0] = q.r[0]; d[1] = q.r[1]; d[2] = q.r[2]; w[0] = q.v[0]; w[1] = q.v[1]; w[2] = q.v[2];
    };
  }
  if (spec.type === 'harmonic') {
    vec3(spec.amplitude, 'motion.amplitude');
    const om = spec.omega ?? (2 * Math.PI / spec.period);
    need(Number.isFinite(om) && om > 0, 'motion needs omega or period > 0');
    const ph = Array.isArray(spec.phase) ? spec.phase : [spec.phase ?? 0, spec.phase ?? 0, spec.phase ?? 0];
    const ramp = spec.ramp ?? 2 * Math.PI / om;
    need(ramp >= 0, 'motion.ramp must be >= 0');
    return (t, d, w) => {
      let R = 1, dR = 0;
      if (t < ramp) { const x = t / ramp; R = x * x * (3 - 2 * x); dR = 6 * x * (1 - x) / ramp; }
      for (let k = 0; k < 3; k++) {
        const sn = Math.sin(om * t + ph[k]), cs = Math.cos(om * t + ph[k]);
        d[k] = spec.amplitude[k] * R * sn;
        w[k] = spec.amplitude[k] * (R * om * cs + dR * sn);
      }
    };
  }
  if (spec.type === 'table') {
    const tt = spec.t;
    need(Array.isArray(tt) && tt.length >= 2, 'motion.t must have at least 2 entries');
    for (let i = 1; i < tt.length; i++) need(tt[i] > tt[i - 1], 'motion.t must be ascending');
    const cols = [spec.x, spec.y, spec.z].map((c, k) => {
      const arr = c ?? new Array(tt.length).fill(0);
      need(arr.length === tt.length, `motion.${'xyz'[k]} must have the same length as motion.t`);
      return arr;
    });
    let i = 0;
    return (t, d, w) => {
      if (t <= tt[0]) { for (let k = 0; k < 3; k++) { d[k] = cols[k][0]; w[k] = 0; } return; }
      const n = tt.length;
      if (t >= tt[n - 1]) { for (let k = 0; k < 3; k++) { d[k] = cols[k][n - 1]; w[k] = 0; } return; }
      if (t < tt[i]) i = 0;
      while (t > tt[i + 1]) i++;
      const h = tt[i + 1] - tt[i], f = (t - tt[i]) / h;
      for (let k = 0; k < 3; k++) { d[k] = cols[k][i] + f * (cols[k][i + 1] - cols[k][i]); w[k] = (cols[k][i + 1] - cols[k][i]) / h; }
    };
  }
  throw new RangeError("motion must be a function, { type: 'harmonic', ... } or { type: 'table', ... }");
}

function energies(m, r, v) {
  let ke = 0, se = 0, pe = 0;
  const { N } = m;
  for (let i = 1; i < N; i++) {
    const j = 3 * i;
    ke += 0.5 * m.mass[i] * (v[j] * v[j] + v[j + 1] * v[j + 1] + v[j + 2] * v[j + 2]);
    pe += m.W[i] * r[j + 2];
  }
  for (let s = 0; s < N; s++) {
    const e = m.len[s] - m.l0[s];
    if (e > 0 || m.allowCompression) se += 0.5 * m.EA[s] * e * e / m.l0[s];
  }
  return { ke, se, pe };
}

/**
 * Time-domain simulation of the line with the fairlead driven by a prescribed 3-D motion and the anchor fixed.
 * The static state is found first (relaxLine) unless the model has already been relaxed.
 *
 * Integration: classical RK4 or semi-implicit (symplectic) Euler, sub-stepped so that the internal step does
 * not exceed cfl x the critical step l0 / sqrt(EA / m) of the stiffest segment (reduced for internal damping).
 *
 * Fairlead motion (displacement relative to the rest position):
 *   { type: 'harmonic', amplitude: [ax, ay, az] [m], period [s] | omega [rad/s], phase [rad] | [3], ramp [s] }
 *     d = R(t) a sin(w t + phase), R a smooth 0 -> 1 ramp over `ramp` (default one period);
 *   { type: 'table', t: [], x: [], y: [], z: [] } linear interpolation;
 *   or a function t -> { r: [dx, dy, dz], v: [vx, vy, vz] }.
 *
 * @param {object} model model from buildLumpedLine (state is advanced in place)
 * @param {object} o
 * @param {object|function} [o.motion] fairlead motion (default: fixed)
 * @param {number} o.duration simulated time [s]
 * @param {number} o.dt output interval [s] (also the largest internal step)
 * @param {'rk4'|'semi-implicit'} [o.integrator='rk4']
 * @param {number} [o.cfl=0.5] internal step / critical step [-]
 * @param {number} [o.dtInternal] fixed internal step [s] (overrides cfl; for time-step convergence studies)
 * @param {number[]} [o.segments] indices of segments whose tension history is stored
 * @param {boolean} [o.quasiStatic=true] also evaluate the quasi-static (catenary) fairlead tension at every
 *   output time, for the dynamic amplification factor
 * @param {number} [o.statsFrom=0] start time of the window used for envelopes and statistics [s]
 * @param {object} [o.relax] options for relaxLine
 * @returns {object} t [s]; fairleadTension, anchorTension (end reactions) [N]; segmentTension (history per
 *   requested segment) [N]; quasiStaticTension [N]; touchdownS [m]; energy { kinetic, strain, potential,
 *   fairleadWork } [J] (structural kinetic energy of interior nodes, elastic strain energy, net-weight potential
 *   energy, cumulative work done on the line by the fairlead); envelopes { Tmax, Tmin per segment [N]; xMin,
 *   xMax, yMin, yMax, zMin, zMax per node [m] }; static { fairleadTension, anchorTension, segmentTension };
 *   stats { fairleadMax, fairleadMin, fairleadMean, anchorMax, anchorMin, touchdownMin, touchdownMax,
 *   quasiStaticMax, daf }; steps (internal steps), dtInternal [s]
 */
export function simulateLine(model, o) {
  const m = model;
  need(Number.isFinite(o.duration) && o.duration > 0, 'duration must be > 0');
  need(Number.isFinite(o.dt) && o.dt > 0 && o.dt <= o.duration, 'dt must be > 0 and <= duration');
  const integ = o.integrator ?? 'rk4';
  need(integ === 'rk4' || integ === 'semi-implicit', "integrator must be 'rk4' or 'semi-implicit'");
  if (!m.Tstatic) {
    const rs = relaxLine(m, o.relax);
    if (!rs.converged) throw new Error('Static relaxation of the line did not converge; raise relax.maxIter or loosen relax.tol');
  }
  const { N, r, v, a } = m;
  const nn = N + 1;
  const hMax = o.dtInternal ?? (o.cfl ?? 0.5) * m.dtCritical;
  need(hMax > 0, 'dtInternal must be > 0');
  const nSub = Math.max(1, Math.ceil(o.dt / hMax - 1e-9));
  const h = o.dt / nSub;
  m.hNow = h;
  const nOut = Math.floor(o.duration / o.dt + 1e-9) + 1;
  const motion = makeMotion(o.motion);
  const f0 = m.fairlead0, d = [0, 0, 0], w = [0, 0, 0];
  const setFair = (t, rr, vv) => {
    motion(t, d, w);
    const j = 3 * N;
    rr[j] = f0[0] + d[0]; rr[j + 1] = f0[1] + d[1]; rr[j + 2] = f0[2] + d[2];
    vv[j] = w[0]; vv[j + 1] = w[1]; vv[j + 2] = w[2];
  };
  const segs = o.segments ?? [];
  for (const s of segs) need(Number.isInteger(s) && s >= 0 && s < N, 'segments[] index out of range');
  const res = {
    t: new Float64Array(nOut), fairleadTension: new Float64Array(nOut), anchorTension: new Float64Array(nOut),
    segmentTension: segs.map(() => new Float64Array(nOut)), touchdownS: new Float64Array(nOut),
    quasiStaticTension: (o.quasiStatic ?? true) ? new Float64Array(nOut) : null,
    energy: { kinetic: new Float64Array(nOut), strain: new Float64Array(nOut), potential: new Float64Array(nOut), fairleadWork: new Float64Array(nOut) },
    envelopes: {
      Tmax: new Float64Array(N).fill(-Infinity), Tmin: new Float64Array(N).fill(Infinity),
      xMin: new Float64Array(nn).fill(Infinity), xMax: new Float64Array(nn).fill(-Infinity),
      yMin: new Float64Array(nn).fill(Infinity), yMax: new Float64Array(nn).fill(-Infinity),
      zMin: new Float64Array(nn).fill(Infinity), zMax: new Float64Array(nn).fill(-Infinity),
    },
    static: { fairleadTension: m.staticFairlead, anchorTension: m.staticAnchor, segmentTension: Float64Array.from(m.Tstatic) },
  };
  const env = res.envelopes;
  const t0s = o.statsFrom ?? 0;
  const rt = new Float64Array(3 * nn), vt = new Float64Array(3 * nn);
  const k1 = new Float64Array(3 * nn), k2 = new Float64Array(3 * nn), k3 = new Float64Array(3 * nn), k4 = new Float64Array(3 * nn);
  const n3 = 3 * N;
  let work = 0, pPrev = 0, first = true;
  let t = m.time;
  const tStart = t;
  setFair(t, r, v);
  let iOut = 0;
  const totalSteps = (nOut - 1) * nSub;
  for (let step = 0; step <= totalSteps; step++) {
    t = tStart + step * h;
    accel(m, t, r, v, integ === 'rk4' ? k1 : a);
    // diagnostics at the start of the step (state and segment tensions are consistent here)
    const jf = 3 * N, sl = N - 1;
    const p = m.T[sl] * (m.ux[sl] * v[jf] + m.uy[sl] * v[jf + 1] + m.uz[sl] * v[jf + 2]);
    if (!first) work += 0.5 * h * (p + pPrev);
    pPrev = p; first = false;
    if (t - tStart >= t0s - 1e-12) {
      for (let s = 0; s < N; s++) { const Ts = m.T[s]; if (Ts > env.Tmax[s]) env.Tmax[s] = Ts; if (Ts < env.Tmin[s]) env.Tmin[s] = Ts; }
      for (let i = 0; i < nn; i++) {
        const j = 3 * i;
        if (r[j] < env.xMin[i]) env.xMin[i] = r[j]; if (r[j] > env.xMax[i]) env.xMax[i] = r[j];
        if (r[j + 1] < env.yMin[i]) env.yMin[i] = r[j + 1]; if (r[j + 1] > env.yMax[i]) env.yMax[i] = r[j + 1];
        if (r[j + 2] < env.zMin[i]) env.zMin[i] = r[j + 2]; if (r[j + 2] > env.zMax[i]) env.zMax[i] = r[j + 2];
      }
    }
    if (step % nSub === 0) {
      res.t[iOut] = t - tStart;
      res.fairleadTension[iOut] = endTension(m, N);
      res.anchorTension[iOut] = endTension(m, 0);
      segs.forEach((s, k) => { res.segmentTension[k][iOut] = m.T[s]; });
      res.touchdownS[iOut] = touchdown(m).s;
      const e = energies(m, r, v);
      res.energy.kinetic[iOut] = e.ke; res.energy.strain[iOut] = e.se; res.energy.potential[iOut] = e.pe; res.energy.fairleadWork[iOut] = work;
      if (res.quasiStaticTension) {
        const X = hypot(r[jf] - r[0], r[jf + 1] - r[1]), Z = r[jf + 2] - r[2];
        let Tq = NaN;
        try { Tq = solveMultiSegment({ X, Z, segments: m.static.segments, points: m.static.points, seabed: m.static.seabed && Z >= 0, stiffness: false }).T; } catch (e2) { Tq = NaN; }
        res.quasiStaticTension[iOut] = Tq;
      }
      iOut++;
    }
    if (step === totalSteps) break;
    if (integ === 'rk4') {
      const hh = 0.5 * h;
      for (let j = 3; j < n3; j++) { rt[j] = r[j] + hh * v[j]; vt[j] = v[j] + hh * k1[j]; }
      rt[0] = r[0]; rt[1] = r[1]; rt[2] = r[2];
      setFair(t + hh, rt, vt);
      accel(m, t + hh, rt, vt, k2);
      for (let j = 3; j < n3; j++) { rt[j] = r[j] + hh * v[j] + hh * hh * k1[j]; vt[j] = v[j] + hh * k2[j]; }
      accel(m, t + hh, rt, vt, k3);
      for (let j = 3; j < n3; j++) { rt[j] = r[j] + h * v[j] + h * hh * k2[j]; vt[j] = v[j] + h * k3[j]; }
      setFair(t + h, rt, vt);
      accel(m, t + h, rt, vt, k4);
      for (let j = 3; j < n3; j++) {
        r[j] += h * v[j] + (h * h / 6) * (k1[j] + k2[j] + k3[j]);
        const vn = v[j] + (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
        v[j] = (vn < 1e-100 && vn > -1e-100) ? 0 : vn; // avoid subnormal arithmetic on resting nodes
      }
    } else {
      for (let j = 3; j < n3; j++) {
        const vn = v[j] + h * a[j];
        v[j] = (vn < 1e-100 && vn > -1e-100) ? 0 : vn;
        r[j] += h * v[j];
      }
    }
    setFair(t + h, r, v);
  }
  m.time = tStart + totalSteps * h;

  // statistics over the window t >= statsFrom
  let fMax = -Infinity, fMin = Infinity, fSum = 0, aMax = -Infinity, aMin = Infinity, tdMin = Infinity, tdMax = -Infinity, qMax = -Infinity, cnt = 0;
  for (let i = 0; i < nOut; i++) {
    if (res.t[i] < t0s - 1e-12) continue;
    cnt++;
    const T = res.fairleadTension[i];
    fMax = Math.max(fMax, T); fMin = Math.min(fMin, T); fSum += T;
    aMax = Math.max(aMax, res.anchorTension[i]); aMin = Math.min(aMin, res.anchorTension[i]);
    tdMin = Math.min(tdMin, res.touchdownS[i]); tdMax = Math.max(tdMax, res.touchdownS[i]);
    if (res.quasiStaticTension && Number.isFinite(res.quasiStaticTension[i])) qMax = Math.max(qMax, res.quasiStaticTension[i]);
  }
  const Ts = m.staticFairlead;
  res.stats = {
    fairleadMax: fMax, fairleadMin: fMin, fairleadMean: fSum / Math.max(cnt, 1), anchorMax: aMax, anchorMin: aMin,
    touchdownMin: tdMin, touchdownMax: tdMax, quasiStaticMax: qMax > -Infinity ? qMax : null,
    daf: qMax > -Infinity && qMax - Ts > 1e-9 * Ts ? (fMax - Ts) / (qMax - Ts) : null,
  };
  res.steps = totalSteps; res.dtInternal = h;
  return res;
}

/**
 * Snap-load screening of a tension history. A slack event starts when the tension falls to or below
 * max(0, slackFraction x reference tension) - the usual "slack line" criterion (minimum tension <= 0, or below
 * about 10 % of the static pretension, as applied in snap-load assessments in the style of DNV-RP-N103 /
 * DNV-OS-E301; the fraction is an argument) - and ends when it rises above that level again. The largest
 * tension reached before the next slack event is reported as the re-tensioning peak; the event is flagged as a
 * snap load when that peak exceeds spikeFactor x the reference tension.
 * @param {object} o
 * @param {ArrayLike<number>} o.t time [s] @param {ArrayLike<number>} o.T tension [N]
 * @param {number} o.Tref reference (static pre-) tension [N]
 * @param {number} [o.slackFraction=0.1] slack threshold as a fraction of Tref [-]
 * @param {number} [o.spikeFactor=1.0] peak ratio above which a re-tensioning is called a snap [-]
 * @returns {{slack:boolean,snap:boolean,events:Array<{tSlack:number,tRetension:number,slackDuration:number,
 *   tPeak:number,peak:number,ratio:number,snap:boolean}>,minTension:number,maxTension:number,threshold:number}}
 *   ratio = peak / Tref [-]
 */
export function detectSnapLoads({ t, T, Tref, slackFraction = 0.1, spikeFactor = 1.0 }) {
  need(t && T && t.length === T.length && t.length >= 2, 't and T must have the same length (>= 2)');
  need(Number.isFinite(Tref) && Tref > 0, 'Tref must be > 0');
  need(slackFraction >= 0 && slackFraction < 1, 'slackFraction must be within [0, 1)');
  need(spikeFactor > 0, 'spikeFactor must be > 0');
  const thr = Math.max(0, slackFraction * Tref);
  const events = [];
  let mn = Infinity, mx = -Infinity, cur = null;
  let state = 'pre'; // 'pre' | 'slack' | 'rise' (tracking the largest tension until the next slack event)
  for (let i = 0; i < T.length; i++) {
    const Ti = T[i];
    mn = Math.min(mn, Ti); mx = Math.max(mx, Ti);
    if (state === 'slack') {
      if (Ti > thr) {
        cur.tRetension = t[i]; cur.slackDuration = t[i] - cur.tSlack; cur.peak = Ti; cur.tPeak = t[i];
        state = 'rise';
      }
      continue;
    }
    if (state === 'rise' && Ti > cur.peak) { cur.peak = Ti; cur.tPeak = t[i]; }
    if (Ti <= thr) {
      cur = { tSlack: t[i], tRetension: NaN, slackDuration: NaN, tPeak: NaN, peak: NaN, ratio: NaN, snap: false };
      events.push(cur);
      state = 'slack';
    }
  }
  for (const e of events) {
    if (Number.isFinite(e.peak)) { e.ratio = e.peak / Tref; e.snap = e.peak > spikeFactor * Tref; }
  }
  return { slack: events.length > 0, snap: events.some((e) => e.snap), events, minTension: mn, maxTension: mx, threshold: thr };
}

/**
 * Line-break transient of a moored floater in surge / sway / yaw: the floater (constant mass incl. added mass,
 * linear damping) starts at rest at its intact equilibrium under the mean load; at t0 one line is released and
 * the motion is integrated with RK4, the mooring load being the quasi-static restoring force of the remaining
 * lines at the instantaneous position (mooringsystem.lineForces). This is a simple one-way quasi-static
 * coupling: no line dynamics, frequency-independent added mass, no wave-frequency motion.
 * Reference: API RP 2SK section on transient analysis after a line failure; equations of motion from first
 * principles.
 * @param {object} o
 * @param {object} o.system mooring system definition (mooringsystem.js)
 * @param {number} o.brokenLine index of the released line
 * @param {number[]} [o.Fext] mean external load [Fx, Fy, Fz, Mx, My, Mz] [N, N m]
 * @param {number[]} o.mass [surge, sway] mass incl. added mass [kg] and yaw inertia incl. added inertia [kg m^2]
 * @param {number[]} [o.damping=[0,0,0]] linear damping [N s/m, N s/m, N m s/rad]
 * @param {number} [o.t0=0] time of the break [s]
 * @param {number} o.duration [s] @param {number} o.dt time step [s]
 * @returns {object} t [s]; x, y [m], yaw [rad]; offset [m] from the intact equilibrium; maxTension history [N];
 *   intact, damaged (equilibrium() results); maxOffset [m] (from the zero pose), maxExcursion [m] (from the
 *   intact equilibrium), staticExcursion [m] (intact -> damaged equilibrium), overshootFactor = maxExcursion /
 *   staticExcursion [-]; peakTension [N], peakTensionLine, peakTensionTime [s]; lineMaxTension per line [N]
 */
export function lineBreakTransient(o) {
  need(o.system && Array.isArray(o.system.lines), 'system must be a mooring system definition');
  const nL = o.system.lines.length;
  need(Number.isInteger(o.brokenLine) && o.brokenLine >= 0 && o.brokenLine < nL, 'brokenLine index out of range');
  need(Array.isArray(o.mass) && o.mass.length === 3 && o.mass.every((x) => x > 0), 'mass must be [mx, my, Izz] > 0');
  const C = o.damping ?? [0, 0, 0];
  need(C.length === 3 && C.every((x) => x >= 0), 'damping must be three values >= 0');
  need(o.duration > 0 && o.dt > 0 && o.dt <= o.duration, 'duration and dt must be > 0 with dt <= duration');
  const Fext = o.Fext ?? [0, 0, 0, 0, 0, 0];
  const t0 = o.t0 ?? 0;
  const intact = equilibrium(o.system, { Fext });
  const damaged = equilibrium(o.system, { Fext, exclude: [o.brokenLine], pose0: intact.pose });
  const base = intact.pose;
  const warmI = [], warmD = [];
  let lastLines = null;
  const acc = (t, y, out) => {
    const pose = [y[0], y[1], base[2], base[3], base[4], y[2]];
    const broken = t >= t0;
    const r = lineForces(o.system, pose, broken ? { exclude: [o.brokenLine], warm: warmD } : { warm: warmI });
    lastLines = r;
    out[0] = y[3]; out[1] = y[4]; out[2] = y[5];
    out[3] = (r.F[0] + Fext[0] - C[0] * y[3]) / o.mass[0];
    out[4] = (r.F[1] + Fext[1] - C[1] * y[4]) / o.mass[1];
    out[5] = (r.F[5] + Fext[5] - C[2] * y[5]) / o.mass[2];
  };
  const n = Math.floor(o.duration / o.dt + 1e-9) + 1;
  const res = { t: new Float64Array(n), x: new Float64Array(n), y: new Float64Array(n), yaw: new Float64Array(n), offset: new Float64Array(n), maxTension: new Float64Array(n) };
  const y = [base[0], base[1], base[5], 0, 0, 0];
  const k1 = new Array(6), k2 = new Array(6), k3 = new Array(6), k4 = new Array(6), yt = new Array(6);
  const lineMax = new Array(nL).fill(0);
  let peak = 0, peakLine = -1, peakTime = 0, maxExc = 0, maxOff = 0;
  const h = o.dt;
  for (let i = 0; i < n; i++) {
    const t = i * h;
    acc(t, y, k1);
    res.t[i] = t; res.x[i] = y[0]; res.y[i] = y[1]; res.yaw[i] = y[2];
    res.offset[i] = hypot(y[0] - base[0], y[1] - base[1]);
    res.maxTension[i] = lastLines.maxTension;
    maxExc = Math.max(maxExc, res.offset[i]); maxOff = Math.max(maxOff, hypot(y[0], y[1]));
    lastLines.lines.forEach((l, k) => { if (l.T > lineMax[k]) lineMax[k] = l.T; });
    if (lastLines.maxTension > peak) { peak = lastLines.maxTension; peakLine = lastLines.maxTensionLine; peakTime = t; }
    if (i === n - 1) break;
    for (let j = 0; j < 6; j++) yt[j] = y[j] + 0.5 * h * k1[j];
    acc(t + 0.5 * h, yt, k2);
    for (let j = 0; j < 6; j++) yt[j] = y[j] + 0.5 * h * k2[j];
    acc(t + 0.5 * h, yt, k3);
    for (let j = 0; j < 6; j++) yt[j] = y[j] + h * k3[j];
    acc(t + h, yt, k4);
    for (let j = 0; j < 6; j++) y[j] += (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
  }
  const stat = hypot(damaged.pose[0] - base[0], damaged.pose[1] - base[1]);
  return {
    ...res, intact, damaged, maxOffset: maxOff, maxExcursion: maxExc, staticExcursion: stat,
    overshootFactor: stat > 0 ? maxExc / stat : null, peakTension: peak, peakTensionLine: peakLine, peakTensionTime: peakTime,
    lineMaxTension: lineMax,
  };
}
