// Shared inputs and helpers of the anchor modules: design loads at the anchor, clay strength profile,
// consolidation set-up and cyclic degradation factors, and safety-factor defaults.
import { ANCHOR_SAFETY_FACTORS, anchorDesignCheck } from '../../engine/geotech/anchors.js';
import { DEG, need, CONFIRM } from './_common.js';

export const LIMIT_ANCHOR = 'Capacity methods are limit-equilibrium / plasticity solutions with a linear undrained strength profile; layered soils, soil–structure finite-element analysis, installation disturbance beyond the set-up factor and detailed cyclic contour-diagram procedures are outside this module. Safety factors are overridable defaults.';

/** Anchor design loads (boundary condition at the anchor end of the line). */
export function anchorLoadInputs({ H = 2750, V = 0 } = {}) {
  const g = 'Design load at the anchor';
  return [
    { key: 'anchorH', label: 'Horizontal line load', unit: 'kN', type: 'number', default: H, min: 0, max: 60000, group: g, role: 'bc', bind: 'mooring.anchorH', uncertain: { dist: 'gumbel', cov: 0.15 },
      help: 'Maximum horizontal load at the anchor (mudline) from the tension-and-offset module, before safety factors.' },
    { key: 'anchorV', label: 'Vertical (uplift) line load', unit: 'kN', type: 'number', default: V, min: 0, max: 60000, group: g, role: 'bc', bind: 'mooring.anchorV',
      help: 'Uplift component at the anchor. Zero for catenary systems with grounded line; taut-leg systems load the anchor at 30–45°.' },
    { key: 'condition', label: 'Mooring condition', type: 'select', default: 'intact', options: [{ value: 'intact', label: 'Intact' }, { value: 'damaged', label: 'One line broken' }], group: g, role: 'criteria' },
    { key: 'service', label: 'Service', type: 'select', default: 'permanent', options: [{ value: 'permanent', label: 'Permanent mooring' }, { value: 'mobile', label: 'Mobile (MODU) mooring' }], group: g, role: 'criteria' },
  ];
}

/** Clay strength profile su = su0 + k·z. */
export function clayInputs({ su0 = 5, k = 1.5 } = {}) {
  return [
    { key: 'su0', label: 'Undrained shear strength at mudline', unit: 'kPa', type: 'number', default: su0, min: 0, max: 500, group: 'Soil', role: 'soil', bind: 'site.su0', uncertain: { dist: 'lognormal', cov: 0.2 },
      help: 'From CPT / laboratory tests (DSS strength is the usual reference for anchors).' },
    { key: 'suGrad', label: 'Strength gradient with depth', unit: 'kPa/m', type: 'number', default: k, min: 0, max: 20, group: 'Soil', role: 'soil', bind: 'site.suGrad', uncertain: { dist: 'lognormal', cov: 0.15 },
      help: 'Normally consolidated clay: typically 1.2–2.0 kPa/m.' },
  ];
}

/** Set-up, cyclic and safety-factor inputs. `kind` picks the default API factor row. */
export function capacityFactorInputs({ sf = 1.5 } = {}) {
  return [
    { key: 'setUp', label: 'Consolidation set-up factor', unit: '-', type: 'number', default: 1.0, min: 0.3, max: 2, group: 'Capacity factors', role: 'ic',
      help: 'Ratio of capacity at the time of loading to the reference static capacity. Initial condition of the soil after installation: below 1 soon after installation in sensitive clay (remoulding), rising towards 1.0–1.5 with consolidation.' },
    { key: 'cyclic', label: 'Cyclic loading factor', unit: '-', type: 'number', default: 0.9, min: 0.4, max: 1.5, group: 'Capacity factors', role: 'soil', uncertain: { dist: 'normal', cov: 0.1 },
      help: 'Cyclic degradation combined with loading-rate gain in the design storm (commonly 0.8–1.2 in clay; site-specific from cyclic DSS testing).' },
    { key: 'sfRequired', label: 'Required safety factor (0 = API RP 2SK table)', unit: '-', type: 'number', default: 0, min: 0, max: 5, group: 'Capacity factors', role: 'criteria',
      help: `Factor of safety on holding capacity. 0 uses the API RP 2SK (3rd ed.) table for the anchor type, service and condition (e.g. ${sf} here for permanent intact). ${CONFIRM}` },
  ];
}

/** Required safety factor from the API table unless overridden. */
export function requiredSf(inp, kind) {
  if (inp.sfRequired > 0) return inp.sfRequired;
  const row = ANCHOR_SAFETY_FACTORS['API RP 2SK'][kind];
  const v = row?.[inp.service]?.[inp.condition];
  need(v != null, `API RP 2SK gives no ${inp.condition} factor for this anchor type in ${inp.service} service: enter the required safety factor explicitly.`);
  return v;
}

/** Resultant and inclination of the anchor load (SI). */
export function anchorLoad(inp) {
  const H = inp.anchorH * 1e3, V = inp.anchorV * 1e3;
  const T = Math.hypot(H, V);
  need(T > 0, 'The anchor load is zero: enter a horizontal or vertical design load.');
  return { H, V, T, angle: Math.atan2(V, H), angleDeg: Math.atan2(V, H) / DEG };
}

export const designCheck = (o) => anchorDesignCheck(o);
