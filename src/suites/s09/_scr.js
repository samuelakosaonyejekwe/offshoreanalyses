// Steel-catenary-riser model shared by the modules that start from a converged static SCR shape
// (static configuration, vessel-motion response, wave fatigue, touchdown zone, hang-off loads).
import { solveScr, stressAlong } from '../../engine/lines/riser.js';
import { pipeInputs, waterInputs, hydroInputs, pipeFrom, run, need, DEG, CONFIRM } from './_common.js';

export function scrInputs({ angle = 12, grounded = 400, flex = 'flex' } = {}) {
  return [
    ...pipeInputs(),
    ...waterInputs(),
    ...hydroInputs(),
    { key: 'hangoffDepth', label: 'Hang-off depth below still water', unit: 'm', type: 'number', default: 15, min: 0, max: 200, group: 'Top end (hang-off)', role: 'bc',
      help: 'Elevation of the hang-off porch / riser balcony below the surface.' },
    { key: 'hangoffAngle', label: 'Nominal hang-off angle from vertical', unit: 'deg', type: 'number', default: angle, min: 3, max: 30, group: 'Top end (hang-off)', role: 'bc', design: true,
      help: 'Sets the catenary: smaller angle = shorter, steeper riser with lower top tension but a tighter touchdown bend.' },
    { key: 'topType', label: 'Top-end rotational restraint', type: 'select', default: flex, group: 'Top end (hang-off)', role: 'bc',
      options: [{ value: 'pinned', label: 'Pinned (ideal free rotation)' }, { value: 'flex', label: 'Flex joint (rotational spring)' }, { value: 'clamped', label: 'Clamped (stress joint)' }] },
    { key: 'kFlex', label: 'Flex-joint rotational stiffness', unit: 'kN·m/deg', type: 'number', default: 20, min: 0.5, max: 2000, group: 'Top end (hang-off)', role: 'bc', showIf: { key: 'topType', in: ['flex'] },
      uncertain: { dist: 'lognormal', cov: 0.25 }, help: 'Elastomeric flex joints: 5–60 kN·m/deg depending on size, pressure and temperature (supplier data).' },
    { key: 'groundedLength', label: 'Pipe modelled on the seabed', unit: 'm', type: 'number', default: grounded, min: 50, max: 3000, group: 'Bottom end (seabed)', role: 'bc',
      help: 'Length kept beyond the nominal touchdown point, pinned at its far end (anchor / PLET).' },
    { key: 'kn', label: 'Seabed vertical stiffness', unit: 'kPa', type: 'number', default: 1000, min: 5, max: 1e5, group: 'Bottom end (seabed)', role: 'soil',
      uncertain: { dist: 'lognormal', cov: 0.5 }, help: 'Force per metre of pipe per metre of penetration (kN/m/m). Soft clay 20–200, stiff clay / sand 1000+.' },
    { key: 'mu', label: 'Axial seabed friction coefficient', unit: '-', type: 'number', default: 0.4, min: 0, max: 1.2, group: 'Bottom end (seabed)', role: 'soil' },
    { key: 'pInt', label: 'Internal pressure at the hang-off', unit: 'bar', type: 'number', default: 200, min: 0, max: 1500, group: 'Internal fluid', role: 'bc',
      help: 'Operating pressure at the top; the contents head is added with depth.' },
  ];
}

export const criteriaInputs = () => [
  { key: 'Cf', label: 'API RP 2RD design case factor Cf', unit: '-', type: 'number', default: 1.0, min: 1, max: 1.5, group: 'Acceptance criteria', role: 'criteria',
    help: `1.0 operating, 1.2 extreme / temporary, 1.35 hydrotest, 1.5 survival (API RP 2RD Table 2). ${CONFIRM}` },
  { key: 'Ca', label: 'API RP 2RD allowable stress factor Ca', unit: '-', type: 'number', default: 0.667, min: 0.4, max: 1, group: 'Acceptance criteria', role: 'criteria',
    help: `Basic allowable = Ca × yield, Ca = 2/3 in API RP 2RD (1998) §5.2.3. ${CONFIRM}` },
];

/** Build and solve the SCR. Seabed is z = 0, hang-off at (0, h). Returns SI results. */
export function buildScr(inp, ctx, { nElem = 160, current, pinned = false } = {}) {
  const pp = pipeFrom(inp, ctx);
  need(pp.section.w > 0, `The riser is buoyant (apparent weight ${Math.round(pp.section.w)} N/m): a catenary riser needs positive submerged weight. Increase wall thickness or contents density.`);
  const h = inp.depth - inp.hangoffDepth;
  need(h > 20, 'Hang-off depth must leave at least 20 m between the hang-off and the seabed.');
  const topCondition = pinned || inp.topType === 'pinned' ? undefined
    : inp.topType === 'flex' ? { type: 'flex', kRot: inp.kFlex * 1e3 / DEG } : { type: 'clamped' };
  const pressure = { rhoInt: inp.rhoInt, pIntRef: (inp.pInt ?? 0) * 1e5, zIntRef: h, rhoExt: pp.rhoW, zSurface: inp.depth, nu: pp.mat.nu };
  const base = run(() => solveScr({
    pipe: pp.section, h, hangoffAngleDeg: inp.hangoffAngle, groundedLength: inp.groundedLength,
    seabed: { kn: inp.kn * 1e3, mu: inp.mu }, topCondition, nElem, current, pressure,
  }), 'The catenary riser did not solve');
  need(base.converged, 'The static riser solution did not converge. Try a larger hang-off angle, a longer grounded length or a stiffer seabed.');
  return { pp, h, base, pressure };
}

/** von Mises stress along a solved riser with the hydrostatic internal / external pressure profile [Pa]. */
export const stressOf = (res, pp, pressure) => stressAlong(res, { D: pp.D, t: pp.t, pIntRef: pressure.pIntRef, zIntRef: pressure.zIntRef, rhoInt: pressure.rhoInt, rhoExt: pressure.rhoExt, zSurface: pressure.zSurface });
