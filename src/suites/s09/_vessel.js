// Host-vessel motion at the riser hang-off. Heave and pitch come from the panel solution of a box of the vessel's
// length, beam and draft (engine/hydro: boxBargeRAO), which stands in for the radiation–diffraction database of the
// actual hull; a user RAO table overrides them.
import { boxBargeRAO } from '../../engine/hydro/floater.js';
import { makeLinear } from '../../engine/math/interp.js';
import { G, run, adv, fmt } from './_common.js';

/** Limit of the box-vessel model, appended to the `limits` of the modules that use it. */
export const LIMIT_VESSEL = ' Vessel heave and pitch are those of a rectangular box with the length, beam and draft entered (three-dimensional panel solution, deep water, zero speed, pitch gyradius 0.25 L, centre of gravity on the waterline): within 5 % of an independent panel code for a box, for waves longer than the beam. The shape of the actual hull, forward speed and shallow water are not represented; a run with wave periods below the wave whose length equals the beam carries a warning unless a hang-off RAO table is given.';

/** Shortest wave period of the checked range of the box-vessel model: deep-water wavelength equal to the beam, s. */
export const vesselPeriodMin = (inp) => Math.sqrt(2 * Math.PI * inp.vesB / G);

/**
 * Warning (advice items of severity 'warn', or none) when the wave period T [s] that carries the response lies
 * below the range in which the box-vessel model was checked against an independent panel solution.
 */
export function vesselRangeAdvice(inp, T) {
  if ((inp.raoTable || []).filter((r) => Number.isFinite(r.T) && r.T > 0).length >= 2 || !(T > 0)) return [];
  const Tmin = vesselPeriodMin(inp);
  if (T >= Tmin) return [];
  return [adv('engineering', 'warn',
    `Wave period ${fmt(T)} s is below the range of the box-vessel motion model: its panel solution is carried down to the wave whose length equals the beam (${fmt(inp.vesB)} m, period ${fmt(Tmin)} s). Shorter waves use the coefficients of that last period, and the hang-off motion there, although small, is not verified.`,
    `Enter the hang-off RAOs of the vessel motion analysis (table of surge and heave at the porch) for periods below ${fmt(Tmin)} s.`,
    'DNV-RP-C205 (2019) §7.3 (wave-frequency motions by radiation–diffraction analysis); comparison of the box panel model with Capytaine 3.0.0 on barges of 20 × 8 × 2 m, 100 × 20 × 5 m and 80 × 24 × 3 m (cross-solver verification of this application).')];
}

export function vesselInputs({ group = 'Host vessel motion' } = {}) {
  return [
    { key: 'vesL', label: 'Vessel length', unit: 'm', type: 'number', default: 300, min: 40, max: 500, group, role: 'bc',
      help: 'Waterline length of the host (FPSO, semi pontoon length or barge). Used by the box-vessel heave / pitch transfer functions.' },
    { key: 'vesB', label: 'Vessel beam', unit: 'm', type: 'number', default: 55, min: 10, max: 120, group, role: 'bc' },
    { key: 'vesT', label: 'Vessel draft', unit: 'm', type: 'number', default: 18, min: 2, max: 40, group, role: 'bc' },
    { key: 'vesHeading', label: 'Wave heading relative to the vessel', unit: 'deg', type: 'number', default: 180, min: 0, max: 360, group, role: 'bc',
      help: '180° = head sea, 90° = beam sea. The riser plane is taken along the vessel axis (the rod model is in-plane only).' },
    { key: 'hangArm', label: 'Hang-off distance from midship', unit: 'm', type: 'number', default: 40, min: 0, max: 250, group, role: 'bc',
      help: 'Longitudinal lever arm that turns pitch into vertical motion at the porch.' },
    { key: 'surgeScale', label: 'Surge transfer-function multiplier', unit: '-', type: 'number', default: 1, min: 0, max: 3, group, role: 'bc',
      help: 'Scales the long-wave surge estimate (orbital motion averaged over the hull). Calibrate against the vessel motion analysis.' },
    { key: 'raoTable', label: 'Hang-off motion RAOs (optional, overrides the vessel model)', type: 'table', group, role: 'bc',
      columns: [{ key: 'T', label: 'Wave period', unit: 's', type: 'number' }, { key: 'surge', label: 'Surge at hang-off', unit: 'm/m', type: 'number' }, { key: 'heave', label: 'Heave at hang-off', unit: 'm/m', type: 'number' }],
      default: [], help: 'Leave empty to use the box-vessel model. Give at least two periods; values are amplitudes per metre of wave amplitude at the porch.' },
  ];
}

/**
 * Surge and vertical motion amplitude at the hang-off per metre of wave amplitude for each ω [rad/s].
 * Heave and pitch: panel solution of the box hull (boxBargeRAO); heave and arm × pitch are added in phase
 * (conservative). Surge: long-wave Froude–Krylov estimate |cos β|·|sin(kL/2)/(kL/2)|·e^(−kT/2).
 */
export function hangoffRao(inp, omega) {
  const rows = (inp.raoTable || []).filter((r) => Number.isFinite(r.T) && r.T > 0).sort((a, b) => a.T - b.T);
  if (rows.length === 1) throw new Error('The hang-off RAO table needs at least two periods (or leave it empty to use the vessel model).');
  if (rows.length >= 2) {
    const fs = makeLinear(rows.map((r) => r.T), rows.map((r) => Math.max(0, r.surge || 0)));
    const fh = makeLinear(rows.map((r) => r.T), rows.map((r) => Math.max(0, r.heave || 0)));
    return { surge: omega.map((w) => fs(2 * Math.PI / w)), heave: omega.map((w) => fh(2 * Math.PI / w)), source: 'user RAO table' };
  }
  const r = run(() => boxBargeRAO({ L: inp.vesL, B: inp.vesB, T: inp.vesT, omega, headingDeg: inp.vesHeading }), 'Vessel transfer functions');
  const cb = Math.abs(Math.cos(inp.vesHeading * Math.PI / 180));
  const surge = omega.map((w) => { const k = w * w / G, q = 0.5 * k * inp.vesL * cb; return inp.surgeScale * cb * Math.abs(q > 1e-6 ? Math.sin(q) / q : 1) * Math.exp(-0.5 * k * inp.vesT); });
  const heave = omega.map((w, i) => Math.abs(r.heave[i]) + inp.hangArm * Math.abs(r.pitch[i]));
  return { surge, heave, source: 'box-vessel panel model' };
}
