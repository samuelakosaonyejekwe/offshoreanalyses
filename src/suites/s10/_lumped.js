// Shared set-up of the single-line lumped-mass dynamic model (line-dynamics and snap-load modules):
// inputs, model construction, static start state and the driven simulation.
import { buildLumpedLine, relaxLine, simulateLine } from '../../engine/lines/lumpedmass.js';
import { solveCatenary } from '../../engine/lines/catenary.js';
import { DEG, need, run, lineInputs, lineFrom, checkSpan } from './_common.js';

export const LIMIT_LM = 'Single line in still water or a steady in-plane current, driven by a prescribed harmonic fairlead motion; the vessel is not coupled back. Lumped-mass model with Morison drag and added mass, linear axial stiffness, no bending or torsion; seabed is a linear normal spring with regularised Coulomb friction (friction force ramps linearly below a small slip speed, so true stick is approximated). Consistent-mass finite elements, Cosserat rods and three-dimensional seabed terrain are outside this module.';

export function dynInputs({ X = 930, ax = 2.0, az = 1.0, period = 10, nCycles = 4, nSeg = 20 } = {}) {
  return [
    ...lineInputs({ X }),
    { key: 'Cdn', label: 'Normal drag coefficient', unit: '-', type: 'number', default: 2.4, min: 0, max: 4, group: 'Hydrodynamics', role: 'load', uncertain: { dist: 'normal', cov: 0.15 },
      help: 'On the nominal diameter. DNV-OS-E301: studless chain 2.4, stud-link 2.6, six-strand rope 1.8, spiral strand 1.6 (1.2 sheathed), fibre rope 1.6. Indicative — calibrate against submerged-line tests where available.' },
    { key: 'Cdt', label: 'Tangential drag coefficient', unit: '-', type: 'number', default: 1.15, min: 0, max: 2, group: 'Hydrodynamics', role: 'load', help: 'Studless chain 1.15, stud-link 1.4, ropes ≈ 0–0.1.' },
    { key: 'Can', label: 'Normal added-mass coefficient', unit: '-', type: 'number', default: 1.0, min: 0, max: 3, group: 'Hydrodynamics', role: 'load' },
    { key: 'zeta', label: 'Axial material damping ratio', unit: '-', type: 'number', default: 0.02, min: 0, max: 0.5, group: 'Hydrodynamics', role: 'material',
      help: 'Damping ratio of a segment\'s own axial mode; a few per cent suppresses numerical ringing without affecting wave-frequency tension. Fibre ropes have real hysteresis of this order.' },
    { key: 'Uc', label: 'Steady current (in the line plane)', unit: 'm/s', type: 'number', default: 0, min: -3, max: 3, group: 'Hydrodynamics', role: 'bc', bind: 'site.Uc100',
      help: 'Uniform current along the anchor → fairlead direction (negative: towards the anchor). The system-level solver ignores current on lines; this module includes it.' },
    { key: 'mu', label: 'Seabed friction coefficient', unit: '-', type: 'number', default: 0.6, min: 0, max: 2, group: 'Seabed', role: 'bc',
      help: 'Regularised Coulomb friction on grounded nodes, in all horizontal directions.' },
    { key: 'kBed', label: 'Seabed normal stiffness', unit: 'kPa/m', type: 'number', default: 3000, min: 50, max: 50000, group: 'Seabed', role: 'bc',
      help: 'Contact pressure per metre of penetration (soft clay ≈ 100–500, sand ≈ 3000–10000). Stiffer seabeds shorten the stable time step.' },
    { key: 'ax', label: 'Fairlead surge amplitude', unit: 'm', type: 'number', default: ax, min: 0, max: 30, group: 'Fairlead motion', role: 'bc',
      help: 'Prescribed horizontal motion amplitude in the line plane (boundary condition at the vessel end).' },
    { key: 'az', label: 'Fairlead heave amplitude', unit: 'm', type: 'number', default: az, min: 0, max: 30, group: 'Fairlead motion', role: 'bc' },
    { key: 'phaseDeg', label: 'Heave phase lead over surge', unit: 'deg', type: 'number', default: 90, min: -180, max: 180, group: 'Fairlead motion', role: 'bc',
      help: '90° gives the elliptical fairlead path typical of wave-frequency motion.' },
    { key: 'period', label: 'Motion period', unit: 's', type: 'number', default: period, min: 2, max: 400, group: 'Fairlead motion', role: 'bc', bind: 'site.Tp100' },
    { key: 'rampCycles', label: 'Start-up ramp', unit: 'cycles', type: 'number', default: 1, min: 0.25, max: 5, group: 'Initial condition', role: 'ic',
      help: 'The line starts at rest in its static equilibrium (dynamic relaxation of the lumped model); the motion is ramped in smoothly over this many cycles and statistics start after it.' },
    { key: 'nCycles', label: 'Simulated cycles', unit: '-', type: 'number', default: nCycles, min: 2, max: 60, group: 'Numerical', role: 'control',
      help: 'Total motion cycles including the ramp. More cycles cost proportionally more time.' },
    { key: 'nSeg', label: 'Line segments', type: 'int', default: nSeg, min: 6, max: 160, group: 'Numerical', role: 'control',
      help: 'Lumped-mass segments along the line. 20 resolves wave-frequency tension to about 1 %; run the refinement study to confirm.' },
    { key: 'cfl', label: 'Time step / critical step', unit: '-', type: 'number', default: 0.5, min: 0.05, max: 0.9, group: 'Numerical', role: 'control',
      help: 'Internal RK4 step as a fraction of the axial-wave stability limit of the stiffest segment.' },
  ];
}

/** Build, relax and drive the line. All values SI. */
export function runLine(inp, ctx) {
  const p = lineFrom(inp, ctx);
  const { L, X, h } = inp;
  checkSpan(L, X, h, p.EA, p.MBL);
  need(Number.isInteger(inp.nSeg) && inp.nSeg >= 2, 'The number of segments must be an integer of at least 2.');
  const cat = run(() => solveCatenary({ X, Z: h, L, w: p.w, EA: p.EA, CB: 0 }), 'Static catenary');
  need(!cat.slack, 'The line is completely slack in this geometry (fairlead too close to the anchor): increase the horizontal distance or shorten the line.');
  // Drag coefficients are quoted on the nominal diameter (DNV-OS-E301: force per length = ½ρ·C·d·u²); the engine applies the
  // normal coefficient to the volume-equivalent diameter and the tangential one to its circumference π·d, hence the rescaling.
  const model = run(() => buildLumpedLine({
    anchor: [0, 0, -h], fairlead: [X, 0, 0], rho: inp.rhoW ?? 1025,
    sections: [{ L, EA: p.EA, mass: p.mass, d: p.dVol, nSeg: inp.nSeg, Cdn: inp.Cdn * p.d / p.dVol, Cdt: inp.Cdt * p.d / (Math.PI * p.dVol), Can: inp.Can, Cat: 0, zeta: inp.zeta }],
    seabed: { z: -h, k: inp.kBed * 1e3, mu: inp.mu },             // kPa/m → Pa/m
    current: inp.Uc ? { speed: Math.abs(inp.Uc), direction: inp.Uc >= 0 ? 0 : Math.PI } : undefined,
  }), 'Lumped-mass model');
  const relax = run(() => relaxLine(model, { tol: 1e-6, maxIter: 400000 }), 'Static relaxation');
  need(relax.converged, 'The static start state did not converge. Use fewer segments, a softer seabed or check that the line can reach the fairlead.');
  const period = inp.period, duration = inp.nCycles * period;
  const moving = inp.ax > 0 || inp.az > 0;
  const steps = duration / (inp.cfl * model.dtCritical);
  need(!(steps > 150000), `This run needs about ${Math.round(steps / 1000)} thousand internal time steps (${inp.nCycles} cycles of ${period} s at a stable step of ${(inp.cfl * model.dtCritical * 1e3).toPrecision(2)} ms). Reduce the number of cycles, the motion period, the number of segments or the seabed stiffness.`);
  const nOut = Math.min(1500, Math.max(200, Math.round(inp.nCycles * 60)));
  const res = run(() => simulateLine(model, {
    motion: moving ? { type: 'harmonic', amplitude: [inp.ax, 0, inp.az], period, phase: [0, 0, inp.phaseDeg * DEG], ramp: inp.rampCycles * period } : undefined,
    duration, dt: duration / nOut, cfl: inp.cfl, statsFrom: Math.min(inp.rampCycles * period, 0.6 * duration), quasiStatic: true,
  }), 'Time-domain simulation');
  return { p, cat, model, relax, res, duration, period };
}
