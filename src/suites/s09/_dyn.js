// Frequency-domain riser dynamics shared by the vessel-motion, coupled-dynamic and wave-fatigue modules:
// harmonic response of the solved riser to hang-off motion, with Morison drag linearised at a stated amplitude.
import { rodHarmonicResponse } from '../../engine/lines/rod.js';
import { jonswap } from '../../engine/metocean/spectra.js';
import { catenaryRiser } from '../../engine/lines/riser.js';
import { G, run } from './_common.js';
import { hangoffRao } from './_vessel.js';

export const dynInputs = () => [
  { key: 'zeta', label: 'Structural damping ratio', unit: '-', type: 'number', default: 0.005, min: 0, max: 0.05, group: 'Dynamic model', role: 'material',
    help: 'Stiffness-proportional damping of the steel pipe (0.3–1 %). Hydrodynamic drag damping is added by the solver.' },
  { key: 'tdzModel', label: 'Touchdown bending model', type: 'select', default: 'travelling', group: 'Dynamic model', role: 'control',
    options: [{ value: 'travelling', label: 'Travelling touchdown point (static profile swept along the seabed)' }, { value: 'frozen', label: 'Frozen contact (linearised, conservative)' }],
    help: 'The linearised solver freezes seabed contact, which piles all touchdown bending onto one node. The travelling model instead sweeps the static curvature profile over the touchdown travel caused by the hang-off motion.' },
  { key: 'nFreq', label: 'Frequencies in the response sweep', type: 'int', default: 8, min: 5, max: 40, group: 'Dynamic model', role: 'control',
    help: 'Transfer functions are solved at this many wave periods and interpolated in between. More points resolve resonant peaks at the cost of run time.' },
  { key: 'Tmin', label: 'Shortest wave period in the sweep', unit: 's', type: 'number', default: 4, min: 2, max: 10, group: 'Dynamic model', role: 'control' },
  { key: 'Tmax', label: 'Longest wave period in the sweep', unit: 's', type: 'number', default: 25, min: 12, max: 40, group: 'Dynamic model', role: 'control' },
];

/** Harmonic response for hang-off surge / heave amplitudes [m] (heave 90° out of phase with surge). */
export function harmonic(base, { omega, xAmp, zAmp, rho, zeta }) {
  const n = base.s.length;
  if (!(Math.abs(xAmp) > 1e-9 || Math.abs(zAmp) > 1e-9)) {
    const z = () => new Float64Array(n);
    return { tensionAmp: z(), momentAmp: z(), curvatureAmp: z(), normalAmp: z(), xAmp: z(), zAmp: z(), converged: true, iterations: 0 };
  }
  // the riser runs towards +x, so vessel surge towards the riser shortens the span; sign is immaterial for a harmonic
  return run(() => rodHarmonicResponse(base, { omega, end: 'A', xAmp, zAmp, phaseZ: Math.PI / 2, rho, zeta, tol: 1e-3 }), 'The harmonic riser response did not solve');
}

const geomCache = new WeakMap();
/**
 * Touchdown travel along the riser [m of arc length] per metre of horizontal (far) and vertical hang-off motion,
 * from the inextensible catenary with fixed total length: X − S is fixed by the anchor position, so
 * d(X − S) = dx for an offset and ∂(X − S)/∂h·dz + ∂(X − S)/∂a·da = 0 for heave (Irvine §1.2 relations).
 */
function tdzGeometry(base, pp) {
  let g = geomCache.get(base);
  if (g) return g;
  const w = pp.section.w, h = base.config.top.z, a0 = base.catenary.a;
  const f = (a, hh) => { const c = catenaryRiser({ h: hh, w, H: w * a }); return { d: c.horizontalSpan - c.suspendedLength, S: c.suspendedLength }; };
  const ea = 1e-4 * a0, eh = 1e-4 * h, f0 = f(a0, h), fa = f(a0 + ea, h), fh = f(a0, h + eh);
  const dda = (fa.d - f0.d) / ea, ddh = (fh.d - f0.d) / eh, dSa = (fa.S - f0.S) / ea, dSh = (fh.S - f0.S) / eh;
  const n = base.s.length, kap = new Float64Array(n);
  for (let k = 0; k < n; k++) kap[k] = Math.abs(base.curvature[k]);
  let kT = n - 1; if (base.tdp) kT = base.tdp.node;
  g = { dSdx: dSa / dda, dSdz: dSh - dSa * ddh / dda, kap, sT: base.s[kT], lambda: Math.sqrt(pp.section.EI / Math.max(base.Te[kT], 1)) };
  geomCache.set(base, g);
  return g;
}
function interpS(s, v, x) {
  const n = s.length;
  if (x <= s[0]) return v[0];
  if (x >= s[n - 1]) return v[n - 1];
  let lo = 0, hi = n - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (s[m] <= x) lo = m; else hi = m; }
  return v[lo] + (v[hi] - v[lo]) * (x - s[lo]) / (s[hi] - s[lo]);
}

/**
 * Stress amplitude along the riser [Pa]: axial from dynamic tension plus outer-fibre bending, added in phase.
 * With `mv = { base, xAmp, zAmp }` the bending part inside the touchdown zone is taken from the static curvature
 * profile travelling ± the touchdown excursion (half the range seen by each material point) instead of the
 * frozen-contact moment amplitude.
 */
export function stressAmp(resp, pp, scf = 1, mv) {
  const As = pp.props.As, c = pp.D / 2 / pp.props.I, out = new Float64Array(resp.tensionAmp.length);
  for (let k = 0; k < out.length; k++) out[k] = scf * (resp.tensionAmp[k] / As + resp.momentAmp[k] * c);
  if (mv && mv.base.catenary) {
    const g = tdzGeometry(mv.base, pp), s = mv.base.s;
    const delta = Math.hypot(g.dSdx * mv.xAmp, g.dSdz * mv.zAmp);
    const lo = g.sT - 8 * g.lambda - delta, hi = g.sT + 4 * g.lambda + delta, Er = pp.mat.E * pp.D / 2;
    for (let k = 0; k < out.length; k++) {
      if (s[k] < lo || s[k] > hi) continue;
      let mn = Infinity, mx = -Infinity;
      for (let q = -4; q <= 4; q++) { const v = interpS(s, g.kap, s[k] + delta * q / 4); if (v < mn) mn = v; if (v > mx) mx = v; }
      out[k] = scf * (resp.tensionAmp[k] / As + 0.5 * (mx - mn) * Er);
    }
  }
  return out;
}
/** Touchdown travel amplitude [m] for given hang-off motion amplitudes. */
export function tdpTravel(base, pp, xAmp, zAmp) { const g = tdzGeometry(base, pp); return Math.hypot(g.dSdx * xAmp, g.dSdz * zAmp); }
export const mvOf = (inp, base, xAmp, zAmp) => (inp.tdzModel === 'frozen' ? undefined : { base, xAmp, zAmp });

/**
 * Transfer functions per metre of wave amplitude at nFreq periods, drag linearised at wave amplitude `aLin` [m]
 * (limited to 1/7 steepness at short periods). Returns ω ascending [rad/s] and, per ω, the response arrays.
 */
export function sweep(base, pp, inp, aLin) {
  const nf = inp.nFreq, om = [];
  if (!(inp.Tmax > inp.Tmin)) throw new Error('The longest sweep period must exceed the shortest.');
  for (let i = 0; i < nf; i++) om.push(2 * Math.PI / (inp.Tmax + (inp.Tmin - inp.Tmax) * i / (nf - 1)));
  const rao = hangoffRao(inp, om);
  const out = [];
  for (let i = 0; i < nf; i++) {
    const T = 2 * Math.PI / om[i];
    const a = Math.max(1e-3, Math.min(aLin, 0.5 * (G * T * T / (2 * Math.PI)) / 7));
    const r = harmonic(base, { omega: om[i], xAmp: rao.surge[i] * a, zAmp: rao.heave[i] * a, rho: pp.rhoW, zeta: inp.zeta });
    out.push({ omega: om[i], T, a, surge: rao.surge[i], heave: rao.heave[i], resp: r, mv: mvOf(inp, base, rao.surge[i] * a, rao.heave[i] * a) });
  }
  return { omega: om, rows: out, source: rao.source };
}

/** Linear interpolation of a per-frequency scalar onto ω. */
export function interpAt(omega, vals, w) {
  const n = omega.length;
  if (w <= omega[0]) return vals[0];
  if (w >= omega[n - 1]) return vals[n - 1];
  let k = 0; while (omega[k + 1] < w) k++;
  return vals[k] + (vals[k + 1] - vals[k]) * (w - omega[k]) / (omega[k + 1] - omega[k]);
}

/** JONSWAP wave spectrum on a uniform ω grid spanning the sweep, with the interpolation weights onto the sweep. */
export function seaGrid(sw, { Hs, Tp, gamma }, m = 120) {
  const so = sw.omega, w0 = so[0], w1 = so[so.length - 1], dw = (w1 - w0) / m, om = [], S = [];
  const lo = new Int32Array(m), fr = new Float64Array(m);
  for (let j = 0, k = 0; j < m; j++) {
    const w = w0 + (j + 0.5) * dw; om.push(w); S.push(jonswap(w, { Hs, Tp, gamma }));
    while (k < so.length - 2 && so[k + 1] < w) k++;
    lo[j] = k; fr[j] = (w - so[k]) / (so[k + 1] - so[k]);
  }
  return { omega: om, S, dw, lo, fr };
}

/** Response variance and zero-crossing rate from a transfer amplitude per metre of wave amplitude. */
export function responseStats(sw, grid, vals) {
  let m0 = 0, m2 = 0;
  for (let j = 0; j < grid.omega.length; j++) { const k = grid.lo[j], h = vals[k] + (vals[k + 1] - vals[k]) * grid.fr[j]; const e = grid.S[j] * h * h * grid.dw; m0 += e; m2 += e * grid.omega[j] ** 2; }
  const sigma = Math.sqrt(m0), Tz = m0 > 0 && m2 > 0 ? 2 * Math.PI * Math.sqrt(m0 / m2) : 0;
  return { m0, m2, sigma, Tz };
}
export const mpm = (sigma, Tz, duration = 10800) => (sigma > 0 && Tz > 0 ? sigma * Math.sqrt(2 * Math.log(Math.max(2, duration / Tz))) : 0);
