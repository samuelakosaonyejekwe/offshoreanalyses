// Power and thrust curve of the generating unit, shared by the energy-yield and array modules: taken from the rotor
// model (BEM schedule), from a pasted table (published / supplier curve) or from a generic cubic curve.
import { cubicPowerCurve } from '../../engine/aero/aep.js';
import { turbineInputs, getTurbine, need, run, SUPPLIER } from './_common.js';

// Indicative 5-MW-class curve at 1.225 kg/m³ (close to the NREL 5-MW reference definition): paste the supplier curve over it.
export const DEFAULT_CURVE = [
  [3, 40, 0.92], [4, 178, 0.88], [5, 404, 0.84], [6, 738, 0.81], [7, 1187, 0.79], [8, 1771, 0.78], [9, 2519, 0.77], [10, 3448, 0.76], [11, 4563, 0.74],
  [12, 5000, 0.53], [13, 5000, 0.39], [14, 5000, 0.30], [15, 5000, 0.24], [16, 5000, 0.19], [17, 5000, 0.16], [18, 5000, 0.13], [19, 5000, 0.11],
  [20, 5000, 0.095], [21, 5000, 0.082], [22, 5000, 0.071], [23, 5000, 0.063], [24, 5000, 0.055], [25, 5000, 0.049],
].map(([u, p, ct]) => ({ u, p, ct }));

export function curveInputs({ group = 'Generating unit' } = {}) {
  const tab = { key: 'curveSource', in: ['table'] }, cub = { key: 'curveSource', in: ['cubic'] }, notRotor = { key: 'curveSource', in: ['table', 'cubic'] };
  return [
    { key: 'curveSource', label: 'Power-curve source', type: 'select', default: 'rotor', group, role: 'geometry',
      options: [{ value: 'rotor', label: 'Rotor model (BEM schedule)' }, { value: 'table', label: 'Pasted power-curve table' }, { value: 'cubic', label: 'Generic cubic curve (Cp, rating, diameter)' }],
      help: `A published or supplier curve can be pasted as a table. ${SUPPLIER}` },
    ...turbineInputs({ group, hub: false }).map((i) => ({ ...i, ...(i.key === 'ratedMW' ? {} : { showIf: { key: 'curveSource', in: ['rotor'] } }) })),
    { key: 'curve', label: 'Power and thrust curve', type: 'table', group, role: 'geometry', showIf: tab,
      columns: [{ key: 'u', label: 'Wind speed', unit: 'm/s', type: 'number' }, { key: 'p', label: 'Electrical power', unit: 'kW', type: 'number' }, { key: 'ct', label: 'Thrust coefficient', unit: '-', type: 'number' }],
      default: DEFAULT_CURVE.map((r) => ({ ...r })), help: 'Ascending wind speeds from cut-in to cut-out at the reference air density of 1.225 kg/m³. The Ct column drives the wake model; leave it at 0 to derive Ct from momentum theory.' },
    { key: 'D', label: 'Rotor diameter', unit: 'm', type: 'number', default: 126, min: 10, max: 320, group, role: 'geometry', showIf: notRotor, design: true },
    { key: 'ratedCubicMW', label: 'Rated power', unit: 'MW', type: 'number', default: 5, min: 0.1, max: 30, group, role: 'geometry', showIf: cub },
    { key: 'cpCubic', label: 'Electrical power coefficient below rated', unit: '-', type: 'number', default: 0.44, min: 0.2, max: 0.52, group, role: 'geometry', showIf: cub, uncertain: { dist: 'normal', cov: 0.04 } },
  ];
}

/** Curve record { u [m/s], p [W], ct [-], D [m], ratedPower [W], cutIn, cutOut, label } at the reference density 1.225 kg/m³. */
export function getCurve(inp) {
  if (inp.curveSource === 'table') {
    const rows = (inp.curve || []).filter((r) => Number.isFinite(r.u) && Number.isFinite(r.p));
    need(rows.length >= 3, 'The power-curve table needs at least three rows (wind speed, power).');
    for (let i = 1; i < rows.length; i++) need(rows[i].u > rows[i - 1].u, `Power-curve table: wind speeds must ascend (row ${i + 1}).`);
    need(rows.every((r) => r.p >= 0 && r.u >= 0), 'Power-curve table: wind speed and power cannot be negative.');
    const p = rows.map((r) => r.p * 1e3), rated = Math.max(...p);
    need(rated > 0, 'The power-curve table has no positive power.');
    const D = inp.D, A = Math.PI * D * D / 4;
    // Ct from the table where given; otherwise ideal actuator-disc momentum theory from the electrical Cp at 94 % efficiency
    const ct = rows.map((r, i) => {
      if (r.ct > 0) { need(r.ct < 1.2, `Power-curve table: thrust coefficient ${r.ct} in row ${i + 1} is not physical.`); return r.ct; }
      const cp = Math.min(p[i] / 0.94 / (0.5 * 1.225 * A * Math.max(r.u, 0.1) ** 3), 0.59);
      let lo = 0, hi = 1 / 3;
      for (let k = 0; k < 50; k++) { const a = 0.5 * (lo + hi); if (4 * a * (1 - a) ** 2 < cp) lo = a; else hi = a; }
      return 4 * lo * (1 - lo);
    });
    return { u: rows.map((r) => r.u), p, ct, D, ratedPower: rated, cutIn: rows[0].u, cutOut: rows[rows.length - 1].u, label: 'pasted table' };
  }
  if (inp.curveSource === 'cubic') {
    const c = run(() => cubicPowerCurve({ ratedPower: inp.ratedCubicMW * 1e6, D: inp.D, Cp: inp.cpCubic, smoothing: 10 }), 'Generic cubic power curve');
    return { u: c.u, p: c.p, ct: c.ct, D: inp.D, ratedPower: c.ratedPower, cutIn: c.u[0], cutOut: c.u[c.u.length - 1], label: 'generic cubic curve' };
  }
  const tb = getTurbine(inp), s = tb.schedule;
  return { u: [...s.U], p: [...s.power], ct: s.ct.map((c) => Math.min(Math.max(c, 0), 0.999)), D: tb.D, ratedPower: tb.ratedPower, cutIn: s.U[0], cutOut: s.U[s.U.length - 1], label: tb.rotor.name };
}
