// Ring-stiffness inputs and builders shared by the dry-collapse and wet-collapse modules.
import { interlockedRingStiffness, solidRingStiffness } from '../../engine/xsection/flexcollapse.js';
import { carcassBucklingFactor } from '../../engine/xsection/flexible.js';
import { need, run, CONFIRM } from './_common.js';

export function ringInputs() {
  const g = 'Carcass and pressure-armour ring stiffness';
  return [
    { key: 'ringMethod', label: 'Carcass ring stiffness from', type: 'select', default: 'profile', options: [{ value: 'profile', label: 'Interlocked profile (I_min, pitch, K)' }, { value: 'equivalent', label: 'Equivalent solid-ring thickness' }], group: g, role: 'geometry',
      help: 'The interlocked strip is replaced by a ring of the same bending stiffness: EI = K·n·E·I_min/pitch, or E·t_eq³/12. The profile itself (contact, unlocking) is not modelled.' },
    { key: 'Imin', label: 'Carcass profile minimum second moment of area', unit: 'mm⁴', type: 'number', default: 290, min: 1, max: 1e5, group: g, role: 'geometry', showIf: { key: 'ringMethod', in: ['profile'] }, design: true,
      uncertain: { dist: 'normal', cov: 0.08 }, help: 'Smallest second moment of the folded strip cross-section, from the manufacturer\'s profile drawing.' },
    { key: 'pitch', label: 'Carcass pitch', unit: 'mm', type: 'number', default: 16, min: 2, max: 100, group: g, role: 'geometry', showIf: { key: 'ringMethod', in: ['profile'] }, help: 'Axial advance per turn of the carcass strip.' },
    { key: 'nTendons', label: 'Number of carcass strips', type: 'int', default: 1, min: 1, max: 4, group: g, role: 'geometry', showIf: { key: 'ringMethod', in: ['profile'] } },
    { key: 'Kring', label: 'Profile factor K', unit: '-', type: 'number', default: 1, min: 0.3, max: 1.5, group: g, role: 'geometry', showIf: { key: 'ringMethod', in: ['profile'] },
      help: 'Accounts for lay angle and overlap of the profile; it is profile specific and should come from collapse tests (API RP 17B).' },
    { key: 'tEq', label: 'Carcass equivalent ring thickness', unit: 'mm', type: 'number', default: 6, min: 0.5, max: 40, group: g, role: 'geometry', showIf: { key: 'ringMethod', in: ['equivalent'] }, design: true,
      help: 'Thickness of the solid ring with the same bending stiffness as the carcass.' },
    { key: 'paTeq', label: 'Pressure-armour equivalent ring thickness', unit: 'mm', type: 'number', default: 6.4, min: 0, max: 40, group: g, role: 'geometry',
      help: 'Solid-ring thickness equivalent to the Z/C/T profile in ring bending (about 0.7–0.9 of the layer thickness). 0 = ignore the pressure armour.' },
    { key: 'ovality', label: 'Initial ovality (Dmax − Dmin)/(Dmax + Dmin)', unit: '%', type: 'number', default: 0.2, min: 0, max: 5, group: 'Imperfections', role: 'ic',
      uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Manufacturing ovality of the carcass. API RP 17B suggests 0.2 % as a minimum when no measurements exist. ' + CONFIRM },
    { key: 'curvFactor', label: 'Curvature knock-down factor', unit: '-', type: 'number', default: 1, min: 0.3, max: 1, group: 'Imperfections', role: 'ic',
      help: 'Collapse pressure of the bent pipe over that of the straight pipe (1 = straight). Take it from collapse tests at the minimum bend radius; typical values 0.85–0.95.' },
    { key: 'factorRule', label: 'Permissible collapse utilisation', type: 'select', default: 'depth', options: [{ value: 'depth', label: 'API 17J depth rule (0.67 at ≤ 300 m to 0.85 at ≥ 900 m)' }, { value: 'manual', label: 'Manual value' }], group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'factorManual', label: 'Manual permissible utilisation', unit: '-', type: 'number', default: 0.85, min: 0.3, max: 1, group: 'Acceptance criteria', role: 'criteria', showIf: { key: 'factorRule', in: ['manual'] } },
  ];
}

/** Carcass and pressure-armour ring records (SI) for the collapse engines. */
export function ringsFrom(inp, S) {
  need(S.idx.carcass >= 0, 'The stack has no carcass: a smooth-bore pipe relies on the pressure armour and bore pressure against collapse and is outside this module.');
  const C = S.stack.layers[S.idx.carcass];
  need(C.yield > 0, 'Enter the carcass yield strength in the layer table.');
  const ring = inp.ringMethod === 'profile'
    ? run(() => interlockedRingStiffness({ E: C.E, Imin: inp.Imin * 1e-12, pitch: inp.pitch / 1e3, nTendons: inp.nTendons, K: inp.Kring }), 'Carcass ring stiffness')
    : { ...solidRingStiffness({ E: C.E, thickness: inp.tEq / 1e3 }), tEq: inp.tEq / 1e3 };
  // membrane area per unit length = metal fraction × layer thickness; extreme fibre at half the layer thickness
  const a = C.fill * C.t;
  const carcass = { EI: ring.EI, R: C.R, thickness: a, sectionModulus: 2 * (ring.EI / C.E) / C.t, yieldStress: C.yield };
  let pa = null;
  if (S.idx.pa >= 0 && inp.paTeq > 0) {
    const P = S.stack.layers[S.idx.pa];
    pa = { EI: solidRingStiffness({ E: P.E, thickness: inp.paTeq / 1e3 }).EI, R: P.R };
  }
  const factor = inp.factorRule === 'depth' ? carcassBucklingFactor(inp.depth).factor : inp.factorManual;
  return { carcass, pa, factor, tEq: ring.tEq, C };
}
