// Host-vessel motion at the riser hang-off. The closed-form box-vessel transfer functions stand in for a
// potential-flow (diffraction/radiation) database, which is outside this suite; a user RAO table overrides them.
import { boxBargeRAO } from '../../engine/hydro/floater.js';
import { makeLinear } from '../../engine/math/interp.js';
import { G, run } from './_common.js';

export function vesselInputs({ group = 'Host vessel motion' } = {}) {
  return [
    { key: 'vesL', label: 'Vessel length', unit: 'm', type: 'number', default: 300, min: 40, max: 500, group, role: 'bc',
      help: 'Waterline length of the host (FPSO, semi pontoon length or barge). Used by the closed-form heave / pitch transfer functions.' },
    { key: 'vesB', label: 'Vessel beam', unit: 'm', type: 'number', default: 55, min: 10, max: 120, group, role: 'bc' },
    { key: 'vesT', label: 'Vessel draft', unit: 'm', type: 'number', default: 18, min: 2, max: 40, group, role: 'bc' },
    { key: 'vesHeading', label: 'Wave heading relative to the vessel', unit: 'deg', type: 'number', default: 180, min: 0, max: 360, group, role: 'bc',
      help: '180° = head sea, 90° = beam sea. The riser plane is taken along the vessel axis (the rod model is in-plane only).' },
    { key: 'hangArm', label: 'Hang-off distance from midship', unit: 'm', type: 'number', default: 60, min: 0, max: 250, group, role: 'bc',
      help: 'Longitudinal lever arm that turns pitch into vertical motion at the porch.' },
    { key: 'surgeScale', label: 'Surge transfer-function multiplier', unit: '-', type: 'number', default: 1, min: 0, max: 3, group, role: 'bc',
      help: 'Scales the long-wave surge estimate (orbital motion averaged over the hull). Calibrate against the vessel motion analysis.' },
    { key: 'raoTable', label: 'Hang-off motion RAOs (optional, overrides the vessel model)', type: 'table', group, role: 'bc',
      columns: [{ key: 'T', label: 'Wave period', unit: 's', type: 'number' }, { key: 'surge', label: 'Surge at hang-off', unit: 'm/m', type: 'number' }, { key: 'heave', label: 'Heave at hang-off', unit: 'm/m', type: 'number' }],
      default: [], help: 'Leave empty to use the closed-form vessel model. Give at least two periods; values are amplitudes per metre of wave amplitude at the porch.' },
  ];
}

/**
 * Surge and vertical motion amplitude at the hang-off per metre of wave amplitude for each ω [rad/s].
 * Heave and pitch: Jensen, Mansour & Olsen (2004) closed forms (amplitudes only, so heave and arm × pitch are added
 * in phase — conservative). Surge: long-wave Froude–Krylov estimate |cos β|·|sin(kL/2)/(kL/2)|·e^(−kT/2).
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
  // the closed forms carry the sign of the sinc term; only amplitudes are meaningful here
  const heave = omega.map((w, i) => Math.abs(r.heave[i]) + inp.hangArm * Math.abs(r.pitch[i]));
  return { surge, heave, source: 'closed-form box-vessel model' };
}
