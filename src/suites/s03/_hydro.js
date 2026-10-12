// Shared hydrodynamic loading on a vertical surface-piercing pile: Morison drag and inertia integrated over depth
// with wave and current velocities added before squaring (nonlinear wave–current interaction in the drag term),
// stepped through the wave phase. Kinematics from Airy theory (integration to the still-water level) or from the
// Stokes fifth-order or stream-function solution (integration to the instantaneous free surface).
import { airyKinematics, wavenumber, selectWaveTheory, breakingLimits, stokes5Wave, stokes5Kinematics } from '../../engine/metocean/waves.js';
import { streamFunctionWave, streamFunctionKinematics } from '../../engine/metocean/streamfunction.js';
import { currentPowerLaw } from '../../engine/metocean/current.js';
import { maccamyFuchs } from '../../engine/hydro/morison.js';
import { need, run } from './_common.js';

// Stream-function solutions are deterministic functions of (H, d, T): they are kept for the session so that
// depth-refinement and sensitivity studies do not repeat the nonlinear solve.
const SFCACHE = new Map();

/**
 * Wind drag on the tower and on the dry part of the monopile in a power-law sheared wind, about the mudline.
 * @returns {{F:number, M:number}} force N and moment N·m
 */
export function towerWindLoad({ V, hubHeight, shear, cd, rhoAir, D, zInt, towerDb, towerDt, zTop, depth, level = 0, n = 60 }) {
  const zLow = Math.max(level, 0.1);
  let F = 0, M = 0;
  for (let i = 0; i < n; i++) {
    const dz = (zTop - zLow) / n, z = zLow + (i + 0.5) * dz;
    const Dz = z < zInt ? D : towerDb + ((towerDt - towerDb) * (z - zInt)) / (zTop - zInt);
    const U = V * (z / hubHeight) ** shear;
    const f = 0.5 * rhoAir * cd * Dz * U * U * dz;
    F += f; M += f * (z + depth);
  }
  return { F, M };
}

/**
 * Wave and current load on a vertical pile.
 * @param {object} o H m, T s, d m (water depth incl. tide and surge), D m (incl. marine growth), Cd, Cm, mcf bool,
 *   rho kg/m³, uc0 m/s surface current, alpha current-profile exponent, nZ integration stations,
 *   theory 'airy' | 'nonlinear' (Stokes 5th order or stream function to the free surface), nPhase
 * @returns {object} { shearMax N, momentMax N·m (about the mudline, same phase search), phaseDeg, crest m, k rad/m,
 *   L m, cmEff, FD, FI (drag and inertia shear at the governing phase), z[], f[] (force per length N/m at the phase
 *   of maximum moment), theory, recommended, breaking }
 */
export function waveCurrentLoad({ H, T, d, D, Cd, Cm, mcf, rho = 1025, uc0 = 0, alpha = 1 / 7, nZ = 40, theory = 'airy', nPhase = 72 }) {
  need(d > 0 && D > 0, 'Water depth and pile diameter must be positive.');
  const k = wavenumber({ T, d }), L = (2 * Math.PI) / k, omega = (2 * Math.PI) / T;
  const lim = breakingLimits({ T, d });
  const rec = H > 0 ? selectWaveTheory({ H, T, d }) : { theory: 'airy' };
  const breaking = H >= lim.governing;
  const cmEff = mcf && H > 0 ? maccamyFuchs({ H: Math.max(H, 1e-3), T, depth: d, D, rho }).CmEff : Cm;
  const area = (Math.PI * D * D) / 4;
  const uc = (z) => (uc0 > 0 ? currentPowerLaw(Math.min(z, 0), { v0: uc0, d, alpha }) : 0);
  let kin, used = 'airy';
  if (H <= 0) kin = () => ({ u: 0, ax: 0, eta: 0 });
  else if (theory === 'nonlinear' && !breaking) {
    // Stokes fifth order where it is the recommended theory (deep to intermediate water, moderate steepness);
    // the stream function elsewhere and as the fallback.
    const st5 = rec.theory === 'stokes5';
    const key = `${st5 ? 's5' : 'sf'}|${H}|${d}|${T}`;
    let sol = SFCACHE.get(key);
    if (sol === undefined) {
      try { sol = st5 ? stokes5Wave({ H, d, T }) : streamFunctionWave({ H, d, T, order: Math.max(rec.streamFunctionOrder || 8, 8) }); } catch { sol = null; }
      if (SFCACHE.size > 200) SFCACHE.clear();
      SFCACHE.set(key, sol);
    }
    let useSt5 = st5;
    if (!sol && st5) { try { sol = streamFunctionWave({ H, d, T, order: 10 }); useSt5 = false; } catch { sol = null; } }
    if (sol) {
      used = useSt5 ? 'Stokes fifth order' : 'stream function';
      const kf = useSt5 ? stokes5Kinematics : streamFunctionKinematics;
      kin = (z, ph) => { const q = kf(sol, { x: 0, z, t: ph / omega, rho }); return { u: q.u, ax: q.axTotal, eta: q.eta }; };
    }
  }
  if (!kin) {
    // Airy; with the nonlinear request (breaking-limited or solver failure) Wheeler stretching carries kinematics to the crest
    const stretching = theory === 'nonlinear' ? 'wheeler' : 'none';
    if (theory === 'nonlinear') used = 'Airy with Wheeler stretching';
    kin = (z, ph) => { const q = airyKinematics({ H, T, d, z, t: ph / omega, k, stretching, clipAboveSurface: false, rho }); return { u: q.u, ax: q.ax, eta: q.eta }; };
  }
  let best = null;
  for (let i = 0; i < nPhase; i++) {
    const ph = (2 * Math.PI * i) / nPhase;
    const eta = theory === 'nonlinear' && H > 0 ? kin(0, ph).eta : 0;
    const top = Math.max(eta, -0.9 * d), dz = (top + d) / nZ;
    let V = 0, Mo = 0, FD = 0, FI = 0;
    const zs = [], fs = [];
    for (let j = 0; j < nZ; j++) {
      const z = -d + (j + 0.5) * dz;
      const q = kin(Math.min(z, top), ph);
      const ur = q.u + uc(z);
      const fd = 0.5 * rho * Cd * D * ur * Math.abs(ur), fi = rho * cmEff * area * q.ax;
      FD += fd * dz; FI += fi * dz; V += (fd + fi) * dz; Mo += (fd + fi) * dz * (z + d);
      zs.push(z); fs.push(fd + fi);
    }
    if (!best || Math.abs(Mo) > Math.abs(best.momentMax)) best = { shearAtM: V, momentMax: Mo, phaseDeg: (ph * 180) / Math.PI, crest: eta, FD, FI, z: zs, f: fs, shearMax: best ? best.shearMax : 0 };
    if (Math.abs(V) > Math.abs(best.shearMax)) best.shearMax = V;
  }
  let crest = 0;
  if (H > 0) { crest = theory === 'nonlinear' ? run(() => kin(0, 0).eta, 'Wave crest') : H / 2; }
  return { ...best, shearMax: Math.abs(best.shearMax), momentMax: Math.abs(best.momentMax), crestElevation: crest, k, L, cmEff, theory: used, recommended: rec.theory, breaking, Hb: lim.governing, ka: 0.5 * k * D };
}
