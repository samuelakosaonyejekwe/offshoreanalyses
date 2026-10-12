// Environmental loads and static response of the monopile + tower cantilever, shared by the strength, extreme-load
// and foundation modules: Morison wave and current load on the pile (Airy kinematics with Wheeler stretching,
// MacCamy–Fuchs inertia coefficient for large diameters), wind drag on the tower, section forces by direct
// integration from the top and deflection by integration of the curvature on the foundation springs.
import { airyKinematics, wavenumber } from '../../engine/metocean/waves.js';
import { maccamyFuchs } from '../../engine/hydro/morison.js';
import { need, sectionAt, G, clamp } from './_common.js';

/**
 * Distributed wave + current load on the vertical pile at the phase that maximises the mudline moment.
 * z = 0 at still water; stations from the mudline to the instantaneous crest.
 * @returns {{z:number[], f:number[] (N/m), shear:number (N), moment:number (N·m about the mudline), phaseDeg:number,
 *   eta:number (m), cm:number, k:number, ka:number, uMax:number (m/s)}}
 */
export function wavePileLoad({ H, T, depth, D, Cd = 1.0, Cm = 2.0, Uc = 0, curExp = 1 / 7, rho = 1025, n = 60, stretching = 'wheeler', diffraction = true }) {
  if (!(H > 0)) {
    // current only
    const z = [], f = []; let shear = 0, moment = 0;
    for (let i = 0; i < n; i++) { const zz = -depth + (i + 0.5) * depth / n, u = Uc * ((zz + depth) / depth) ** curExp, q = 0.5 * rho * Cd * D * u * Math.abs(u); z.push(zz); f.push(q); shear += q * depth / n; moment += q * depth / n * (zz + depth); }
    return { z, f, shear, moment, phaseDeg: 0, eta: 0, cm: Cm, k: 0, ka: 0, uMax: Uc, dz: depth / n };
  }
  need(H < 0.9 * depth, `The wave height (${H.toFixed(1)} m) exceeds the depth-limited breaking height in ${depth} m of water: reduce it to about 0.78 × depth.`);
  const k = wavenumber({ T, d: depth });
  const ka = k * D / 2;
  // MacCamy–Fuchs effective inertia coefficient where diffraction matters (D/λ > 0.2 ⇔ ka > 0.63); never above the Morison value
  const cm = diffraction ? Math.min(Cm, maccamyFuchs({ H, T, depth, D, rho }).CmEff) : Cm;
  let best = null;
  for (let ip = 0; ip < 72; ip++) {
    const t = ip / 72 * T;
    const eta = airyKinematics({ H, T, d: depth, z: 0, t, k, stretching: 'none', clipAboveSurface: false }).eta;
    const top = stretching === 'none' ? 0 : eta, L = top + depth, dz = L / n;
    let shear = 0, moment = 0, uMax = 0; const z = new Array(n), f = new Array(n);
    for (let i = 0; i < n; i++) {
      const zz = -depth + (i + 0.5) * dz;
      const kin = airyKinematics({ H, T, d: depth, z: zz, t, k, stretching, clipAboveSurface: false });
      const u = kin.u + Uc * clamp((zz + depth) / depth, 0, 1.2) ** curExp;
      const q = 0.5 * rho * Cd * D * u * Math.abs(u) + rho * cm * Math.PI * D * D / 4 * kin.ax;
      z[i] = zz; f[i] = q; shear += q * dz; moment += q * dz * (zz + depth); uMax = Math.max(uMax, Math.abs(u));
    }
    if (!best || moment > best.moment) best = { z, f, shear, moment, phaseDeg: ip * 5, eta, cm, k, ka, uMax, dz };
  }
  return best;
}

/**
 * Static response of the cantilevered support structure.
 * @param {object} o segs (support segments), nStations, point loads [{ z, F (N), M (N·m) }], distributed loads
 *   [{ z[], f[] (N/m), dz }], q(z) callback (N/m) for wind, vertical load at the top (N), foundation ('fixed' or
 *   { Kh, Kr, Khr }), E (Pa), rhoW, depth
 * @returns stations with z, D, t, N (compression +), V, M, EI, deflection u and slope
 */
export function cantilever({ segs, n = 80, points = [], distributed = [], wind = null, topWeight = 0, foundation = 'fixed', E = 210e9, depth }) {
  const zMud = segs[0].zBase, zTop = segs[segs.length - 1].zTop, dz = (zTop - zMud) / n;
  const z = Array.from({ length: n + 1 }, (_, i) => zMud + i * dz);
  // lateral load lumped at the stations' tributary lengths (mid-point rule on a 4× finer grid for the wind callback)
  const V = new Array(n + 1).fill(0), M = new Array(n + 1).fill(0), Nax = new Array(n + 1).fill(0);
  const loads = [];   // { z, F, M }
  for (const p of points) loads.push({ z: p.z, F: p.F || 0, M: p.M || 0 });
  for (const d of distributed) for (let i = 0; i < d.z.length; i++) if (d.z[i] <= zTop) loads.push({ z: d.z[i], F: d.f[i] * d.dz, M: 0 });
  if (wind) { const m = 4 * n, h = (zTop - Math.max(zMud, 0)) / m; for (let i = 0; i < m; i++) { const zz = Math.max(zMud, 0) + (i + 0.5) * h; loads.push({ z: zz, F: wind(zz) * h, M: 0 }); } }
  for (let i = 0; i <= n; i++) {
    let v = 0, m = 0;
    for (const L of loads) if (L.z >= z[i] - 1e-9) { v += L.F; m += L.F * (L.z - z[i]) + L.M; }
    V[i] = v; M[i] = m;
  }
  // axial force: RNA weight plus steel weight above, minus buoyancy of the submerged shell steel
  const sec = z.map((zz) => sectionAt(segs, zz));
  const dens = (zz) => (segs.find((s) => zz >= s.zBase - 1e-9 && zz <= s.zTop + 1e-9) || segs[segs.length - 1]).density;
  let w = topWeight; Nax[n] = w;
  for (let i = n - 1; i >= 0; i--) {
    const zm = 0.5 * (z[i] + z[i + 1]), s = sectionAt(segs, zm), A = Math.PI * (s.D - s.t) * s.t;
    w += A * dz * (dens(zm) - (zm < 0 ? 1025 : 0)) * G; Nax[i] = w;
  }
  const EI = sec.map((s) => E * Math.PI / 64 * (s.D ** 4 - (s.D - 2 * s.t) ** 4));
  // mudline displacement and rotation from the foundation springs: [H; M] = [[Kh, Khr], [Khr, Kr]]·[u; θ]
  let u0 = 0, th0 = 0;
  if (foundation !== 'fixed') { const det = foundation.Kh * foundation.Kr - foundation.Khr ** 2; u0 = (foundation.Kr * V[0] - foundation.Khr * M[0]) / det; th0 = (foundation.Kh * M[0] - foundation.Khr * V[0]) / det; }
  const u = new Array(n + 1).fill(0), th = new Array(n + 1).fill(0); u[0] = u0; th[0] = th0;
  for (let i = 1; i <= n; i++) {
    const c0 = M[i - 1] / EI[i - 1], c1 = M[i] / EI[i];
    th[i] = th[i - 1] + 0.5 * (c0 + c1) * dz;
    u[i] = u[i - 1] + th[i - 1] * dz + (2 * c0 + c1) * dz * dz / 6;
  }
  return { z, D: sec.map((s) => s.D), t: sec.map((s) => s.t), N: Nax, V, M, EI, u, slope: th, dz, zMud, zTop, depth };
}

/** Wind drag per unit height on the tower [N/m] for a power-law profile referred to hub height. */
export const towerWind = ({ segs, Uhub, zHub, alpha = 0.14, Cd = 0.7, rhoAir = 1.225 }) => (z) => {
  if (z <= 0) return 0;
  const U = Uhub * (Math.max(z, 1) / zHub) ** alpha, D = sectionAt(segs, z).D;
  return 0.5 * rhoAir * Cd * D * U * U;
};
