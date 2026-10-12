// Shared vessel and environment inputs and the mean environmental force model (wind, current, wave drift)
// used by the extreme-environment, coupled-analysis and reliability modules.
import { projectedAreaForce, meanDriftForce } from '../../engine/hydro/drift.js';
import { jonswap } from '../../engine/metocean/spectra.js';
import { DEG, G } from './_common.js';

export const RHO_AIR = 1.226;
export const LIMIT_ENV = 'Environmental forces use projected-area drag for wind and current and a constant wave-reflection coefficient for mean drift (energy above the draught only); coefficients are user inputs, not vessel-specific wind-tunnel or diffraction data. Full potential-flow diffraction/radiation and full quadratic transfer functions are outside this module.';

export function envInputs({ Hs = 12, Tp = 14.5, Uw = 36, Uc = 1.0, suffix = '100', group = 'Environment (design storm)' } = {}) {
  return [
    { key: 'Hs', label: 'Significant wave height', unit: 'm', type: 'number', default: Hs, min: 0, max: 25, group, role: 'environment', bind: `site.Hs${suffix}`, uncertain: { dist: 'gumbel', cov: 0.1 },
      help: 'Design sea state (100-year for permanent moorings; 5–10 year for mobile units away from other structures).' },
    { key: 'Tp', label: 'Spectral peak period', unit: 's', type: 'number', default: Tp, min: 3, max: 30, group, role: 'environment', bind: `site.Tp${suffix}` },
    { key: 'gammaJ', label: 'JONSWAP peak factor', unit: '-', type: 'number', default: 2.5, min: 1, max: 7, group, role: 'environment', help: '1 = Pierson–Moskowitz; 2–3.3 for storm seas.' },
    { key: 'Uw', label: 'Wind speed (1-hour mean at 10 m)', unit: 'm/s', type: 'number', default: Uw, min: 0, max: 80, group, role: 'environment', bind: `site.Uw${suffix}`, uncertain: { dist: 'gumbel', cov: 0.1 } },
    { key: 'Uc', label: 'Surface current speed', unit: 'm/s', type: 'number', default: Uc, min: 0, max: 4, group, role: 'environment', bind: `site.Uc${suffix}`, uncertain: { dist: 'gumbel', cov: 0.15 } },
  ];
}

export function hullInputs({ group = 'Vessel' } = {}) {
  return [
    { key: 'Lv', label: 'Hull length at waterline', unit: 'm', type: 'number', default: 100, min: 10, max: 450, group, role: 'geometry' },
    { key: 'Bv', label: 'Hull beam at waterline', unit: 'm', type: 'number', default: 80, min: 5, max: 150, group, role: 'geometry' },
    { key: 'draught', label: 'Draught', unit: 'm', type: 'number', default: 20, min: 1, max: 250, group, role: 'geometry' },
    { key: 'Afront', label: 'Wind area seen head-on', unit: 'm²', type: 'number', default: 2000, min: 0, max: 30000, group, role: 'geometry', help: 'Projected above-water area facing the bow, including topsides and derrick.' },
    { key: 'Aside', label: 'Wind area seen beam-on', unit: 'm²', type: 'number', default: 2400, min: 0, max: 40000, group, role: 'geometry' },
    { key: 'CdWind', label: 'Wind drag coefficient', unit: '-', type: 'number', default: 1.0, min: 0.3, max: 2, group, role: 'load', uncertain: { dist: 'normal', cov: 0.1 }, help: 'On projected area; 0.9–1.2 for bluff offshore units (use wind-tunnel shape and height coefficients where available).' },
    { key: 'AcFront', label: 'Underwater area seen head-on', unit: 'm²', type: 'number', default: 1000, min: 0, max: 20000, group, role: 'geometry', help: 'Projected submerged area: columns and pontoons for a semi-submersible, beam × draught for a ship.' },
    { key: 'AcSide', label: 'Underwater area seen beam-on', unit: 'm²', type: 'number', default: 1200, min: 0, max: 40000, group, role: 'geometry' },
    { key: 'CdCur', label: 'Current drag coefficient', unit: '-', type: 'number', default: 1.0, min: 0.1, max: 2.5, group, role: 'load', uncertain: { dist: 'normal', cov: 0.15 } },
    { key: 'R2', label: 'Wave-drift reflection coefficient R²', unit: '-', type: 'number', default: 0.1, min: 0, max: 1, group, role: 'load', uncertain: { dist: 'lognormal', cov: 0.3 },
      help: 'Share of incident wave energy reflected by the hull (Maruo). About 0.05–0.2 for semi-submersibles, 0.3–0.6 for ship-shaped hulls in beam seas. Replace with diffraction-analysis drift coefficients for design.' },
  ];
}

const DRIFT_GRID = 80;
/** Frequency grid and JONSWAP ordinates covering the energetic band. */
export function seaGrid(inp, n = DRIFT_GRID) {
  const wp = 2 * Math.PI / inp.Tp;
  const omega = Array.from({ length: n }, (_, i) => wp * (0.5 + 3.0 * i / (n - 1)));
  const S = omega.map((w) => (inp.Hs > 0 ? jonswap(w, { Hs: inp.Hs, Tp: inp.Tp, gamma: inp.gammaJ }) : 0));
  return { omega, S };
}

/**
 * Mean drift coefficient D(ω) = F/a² [N/m²] for a projected waterline width W: Maruo's ½·ρ·g·R²·W, with only the
 * wave energy flux above the draught d being reflected (fraction 1 − exp(−2kd), deep-water kinematics).
 */
export function driftCoeff(omega, { R2, width, draught, rho = 1025 }) {
  return omega.map((w) => 0.5 * rho * G * R2 * width * (1 - Math.exp(-2 * w * w / G * draught)));
}

/**
 * Mean environmental force in global axes [N, N, N·m] for wind, waves and current travelling towards the given
 * directions [deg from global x], with the vessel at heading `headingDeg`.
 */
export function envForce(inp, { windDir, waveDir, curDir, headingDeg = 0, rhoW = 1025 }) {
  const comp = (F, relDeg) => {
    const a = (relDeg + headingDeg) * DEG, c = Math.cos(headingDeg * DEG), s = Math.sin(headingDeg * DEG);
    return { Fx: F.Fx * c - F.Fy * s, Fy: F.Fx * s + F.Fy * c, Mz: F.Mz || 0, a };
  };
  const wind = comp(projectedAreaForce({ V: inp.Uw, directionDeg: windDir - headingDeg, Afront: inp.Afront, Aside: inp.Aside, Cfront: inp.CdWind, Cside: inp.CdWind, rho: RHO_AIR }), windDir - headingDeg);
  const cur = comp(projectedAreaForce({ V: inp.Uc, directionDeg: curDir - headingDeg, Afront: inp.AcFront, Aside: inp.AcSide, Cfront: inp.CdCur, Cside: inp.CdCur, rho: rhoW }), curDir - headingDeg);
  const rel = (waveDir - headingDeg) * DEG;
  const width = inp.Bv * Math.abs(Math.cos(rel)) + inp.Lv * Math.abs(Math.sin(rel));
  let Fd = 0;
  if (inp.Hs > 0) { const { omega, S } = seaGrid(inp); Fd = meanDriftForce({ omega, S, D: driftCoeff(omega, { R2: inp.R2, width, draught: inp.draught, rho: rhoW }) }); }
  const drift = { Fx: Fd * Math.cos(waveDir * DEG), Fy: Fd * Math.sin(waveDir * DEG), Mz: 0 };
  const Fx = wind.Fx + cur.Fx + drift.Fx, Fy = wind.Fy + cur.Fy + drift.Fy;
  return { Fx, Fy, Mz: wind.Mz + cur.Mz, F: Math.hypot(Fx, Fy), wind: Math.hypot(wind.Fx, wind.Fy), current: Math.hypot(cur.Fx, cur.Fy), drift: Fd, width };
}
