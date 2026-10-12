// Torsional response: twist under torque, torsional stiffness with free and restrained extension, torque induced
// by tension and pressure (torque balance) and the torque limits in the lay-tightening and lay-opening directions.
import { axisymmetricResponse, torsionLimits } from '../../engine/xsection/flexible.js';
import { stackInputs, stackFrom, loadInputs, loadsFrom, capacityOf, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, linspace, sev,
  LIMIT_AXI, GENERIC, STD_17J, STD_17B, CONFIRM, DEG, CO2_MGO, MGO_T_PER_DAY } from './_common.js';

export default {
  id: 's14-torsion',
  suite: 14,
  title: 'Torsional response',
  summary: 'Gives the twist and wire stress under torque, the torque a tensioned pipe generates when its ends cannot rotate, and the allowable torque in each direction.',
  covers: [4],
  fidelity: 'analytical',
  standards: ['API RP 17B (5th ed.) §8 (torsional failure modes)', 'API Spec 17J (4th ed.) Table 6'],
  models: ['Torsional elasticity of helical layers', 'Coupled axial–torsional stiffness matrix', 'Unilateral interlayer contact (gap opening as bird-caging onset)'],
  limits: LIMIT_AXI + ' The lay-opening limit is the torque at which an interface next to the tensile armour first opens — an onset criterion, not a post-buckling capacity; wire buckling under torsion is not modelled. ' + GENERIC,
  publishes: ['flexible.GJ', 'flexible.torqueLimit'],
  searchTerms: 'flexible pipe torsion tensile armour torque balance birdcaging',
  inputs: [
    ...stackInputs(),
    ...loadInputs({ torque: true }).map((i) => (i.key === 'torque' ? { ...i, default: 20 } : i)),
    { key: 'ends', label: 'Axial end condition', type: 'select', default: 'free', options: [{ value: 'free', label: 'Free to extend (tension prescribed)' }, { value: 'fixed', label: 'Fixed length (axial strain = 0)' }], group: 'Loads and pressure boundaries', role: 'bc',
      help: 'Free: the effective tension is applied and the pipe may lengthen as it twists. Fixed: the length is held and the tension input is ignored.' },
    { key: 'length', label: 'Length twisted', unit: 'm', type: 'number', default: 100, min: 1, max: 10000, group: 'Pipe construction', role: 'geometry', help: 'Length over which the torque acts, for the total end rotation.' },
    { key: 'fTa', label: 'Tensile armour utilisation factor', unit: '-', type: 'number', default: 0.67, min: 0.3, max: 1, group: 'Acceptance criteria', role: 'criteria', help: 'API 17J Table 6 normal operation. ' + CONFIRM },
  ],
  compute(inp, ctx) {
    const S = stackFrom(inp, ctx);
    const base = loadsFrom(inp);
    const tq = base.torque; delete base.torque;
    const withEnds = (o) => (inp.ends === 'fixed' ? { ...o, axialStrain: 0 } : o);
    const r = run(() => axisymmetricResponse(S.stack, withEnds({ ...base, torque: tq })), 'Axisymmetric solution');
    const K = r.stiffness;
    const GJfree = K.K22 - K.K12 * K.K21 / K.K11, GJfixed = K.K22;
    const held = run(() => axisymmetricResponse(S.stack, { ...base, twist: 0 }), 'Restrained-twist solution'); // torsionally restrained ends
    const free = run(() => axisymmetricResponse(S.stack, base), 'Torque-free solution');
    const hasTa = S.idx.ta.length > 0;
    const kpis = [
      kpi('twist', 'Twist', r.twist / DEG, 'deg/m', 'Rotation per metre of pipe under the applied torque, tension and pressure.'),
      kpi('rotation', 'End rotation over the length', r.twist * inp.length / DEG, 'deg', 'Total relative rotation of the two ends. Large rotations must be released by a swivel or allowed for at the end fittings.'),
      kpi('GJ', 'Torsional stiffness (free to extend)', GJfree / 1e3, 'kN·m²', 'Torque per unit twist when the pipe may change length; use this as GJ in global analysis.'),
      kpi('GJfixed', 'Torsional stiffness (fixed length)', GJfixed / 1e3, 'kN·m²', 'Torque per unit twist with the length held; the upper bound of the torsional stiffness.'),
      kpi('coupling', 'Tension–twist coupling', 0.5 * (K.K12 + K.K21) / 1e3, 'kN·m', 'Off-diagonal stiffness term: torque per unit axial strain. Zero for a perfectly torque-balanced armour pair.'),
      kpi('heldTorque', 'Torque induced with ends held', held.torque / 1e3, 'kN·m', 'Torque the pipe generates under tension and pressure when its ends cannot rotate; it must be carried by the end fittings.'),
      kpi('freeTwist', 'Twist of the free pipe', free.twist / DEG, 'deg/m', 'Twist under tension and pressure alone, with no torque — the torque-balance quality of the armour.'),
    ];
    const checks = [], plots = [], tables = [];
    let ucT = 0, lim = null, preOpen = false;
    if (hasTa) {
      const limits = run(() => torsionLimits(S.stack, { loads: base, factor: inp.fTa }), 'Torsion limits');
      // An interface next to the armour that is already open without torque (e.g. no external pressure clamping the outer
      // sheath) makes the gap-onset criterion uninformative: fall back to the wire-stress limit and say so.
      const eff = (l) => (l.separation === 0 && l.stressLimit ? { ...l, governing: l.stressLimit, mode: 'stress', preOpen: true } : { ...l, governing: l.governing ?? 0 });
      const tight = eff(limits.tightening), open = eff(limits.opening);
      preOpen = !!(tight.preOpen || open.preOpen);
      const dir = tq !== 0 && Math.sign(tq) === Math.sign(open.governing || -Math.sign(tight.governing || 1)) ? open : tight;
      lim = Math.min(Math.abs(tight.governing), Math.abs(open.governing));
      ucT = Math.abs(dir.governing) > 0 ? Math.abs(tq) / Math.abs(dir.governing) : 99;
      checks.push(chk(`Applied torque vs ${dir === open ? 'lay-opening' : 'lay-tightening'} limit (${dir.mode})`, tq / 1e3, Math.abs(dir.governing) / 1e3, 'kN·m', `${STD_17B}; wire stress at ${inp.fTa} × structural capacity (${STD_17J})`,
        dir.mode === 'separation' ? 'Guards against the armour layers separating: in the lay-opening direction the outer armour unwinds, loses radial support and can bird-cage.' : 'Guards against overstressing the armour wires in torsion; beyond the limit the wires exceed their permissible stress.'));
      kpis.push(
        kpi('torqueUc', 'Torque utilisation', ucT, '-', 'Applied torque over the allowable torque in its direction. Above 1.0 the torsion criterion is not met.', 1, 'max'),
        kpi('torqueTight', 'Allowable torque, lay-tightening', Math.abs(tight.governing) / 1e3, 'kN·m', `Torque that tightens the outer armour, limited by ${tight.mode === 'stress' ? 'wire stress' : 'layer separation'}.`),
        kpi('torqueOpen', 'Allowable torque, lay-opening', Math.abs(open.governing) / 1e3, 'kN·m', `Torque that unwinds the outer armour, limited by ${open.mode === 'stress' ? 'wire stress' : 'layer separation (bird-caging onset)'}.`),
      );
      const span = 1.2 * Math.max(Math.abs(tight.governing), Math.abs(open.governing), Math.abs(tq), 1e3);
      const ts = linspace(-span, span, 13);
      const sols = ts.map((t) => axisymmetricResponse(S.stack, withEnds({ ...base, torque: t })));
      plots.push(plot('line', 'Armour wire stress vs torque', 'Torque', 'kN·m', 'Wire stress', 'MPa', [
        ...S.idx.ta.map((i) => ({ name: S.stack.layers[i].name, x: ts.map((t) => round(t / 1e3)), y: sols.map((q) => round(q.layers[i].wireStress / 1e6)) })),
        { name: 'Permissible', x: [round(-span / 1e3), round(span / 1e3)], y: Array(2).fill(round(inp.fTa * capacityOf(S.stack.layers[S.idx.ta[0]]) / 1e6)), dash: true },
      ], 'Stress in each tensile armour layer as torque is applied. One layer tightens and picks up load while the other unloads; a kink marks an interface opening.'));
      plots.push(plot('line', 'Torque–twist response', 'Twist', 'deg/m', 'Torque', 'kN·m', [
        { name: 'Response', x: sols.map((q) => round(q.twist / DEG)), y: ts.map((t) => round(t / 1e3)) },
        { name: 'Applied', x: [round(r.twist / DEG)], y: [round(tq / 1e3)] },
      ], 'Torque against twist. The slope is the torsional stiffness; it is lower in the lay-opening direction once layers separate.'));
      tables.push(tbl('Torque limits', [col('dir', 'Direction'), col('stress', 'Wire-stress limit', 'kN·m'), col('sep', 'Separation limit', 'kN·m'), col('gov', 'Governing', 'kN·m'), col('mode', 'Mode')],
        [['Lay-tightening', tight], ['Lay-opening', open]].map(([d, l]) => ({ dir: d, stress: l.stressLimit === null ? null : round(l.stressLimit / 1e3), sep: l.separation === null ? null : round(l.separation / 1e3), gov: round(l.governing / 1e3), mode: l.mode ?? 'none' })),
        'Signed torque limits in each direction. The wire-stress limit is the torque at which a wire reaches its permissible stress; the separation limit is where an armour interface first opens.'));
    } else {
      const span = Math.max(Math.abs(tq), 1e3);
      plots.push(plot('line', 'Torque–twist response', 'Twist', 'deg/m', 'Torque', 'kN·m', [{ name: 'Response', x: [round(-span / GJfree / DEG), round(span / GJfree / DEG)], y: [round(-span / 1e3), round(span / 1e3)] }],
        'Linear torque–twist line of a section without helical armour; the slope is the torsional stiffness.'));
      tables.push(tbl('Stiffness matrix', [col('term', 'Term'), col('value', 'Value')], [{ term: 'K11 [N]', value: round(K.K11) }, { term: 'K12 [N·m]', value: round(K.K12) }, { term: 'K22 [N·m²]', value: round(K.K22) }],
        'Tangent stiffness relating axial force and torque to axial strain and twist.'));
    }
    return {
      kpis, checks, plots, tables,
      notes: [hasTa ? `Torque utilisation ${fmt(ucT)}; the lower of the two directional limits is ${fmt(lim / 1e3)} kN·m.` : 'No tensile armour in the stack: only stiffness and twist are reported.',
        `With the ends held, tension and pressure induce ${fmt(held.torque / 1e3)} kN·m.`,
        ...(preOpen ? ['An interface next to the tensile armour is already open without torque (nothing clamps the outer layers in this load case), so the separation criterion is not informative and the wire-stress limit is used.'] : []),
        ...(hasTa && lim === 0 ? ['The tensile armour already exceeds its permissible stress without torque: no torque can be allowed in this load case.'] : []), GENERIC],
      advice: [
        adv('engineering', hasTa ? sev(ucT) : 'info', hasTa ? `Applied torque ${fmt(tq / 1e3)} kN·m is ${fmt(ucT * 100)} % of the allowable in its direction; end rotation ${fmt(r.twist * inp.length / DEG)}° over ${fmt(inp.length)} m.` : `Twist ${fmt(r.twist / DEG)} deg/m under ${fmt(tq / 1e3)} kN·m.`,
          ucT > 1 ? 'Release the torque with a swivel at the pull-in head or re-orient the lay so that residual torque acts in the lay-tightening direction.' : 'Keep installation-induced twist (reel misalignment, tensioner tracks) in the lay-tightening direction and record the residual twist at tie-in.',
          `${STD_17B}; ${STD_17J} Table 6 (factor is a recalled default — confirm)`),
        adv('economic', 'info', `Torque induced with rotation prevented is ${fmt(Math.abs(held.torque) / 1e3)} kN·m at this tension; untwisting a pipe offshore typically costs of the order of half a vessel day per event.`,
          'Specify the allowable twist per metre in the installation procedure and monitor the twist line during lay — it avoids a recovery and re-lay.',
          'Installation risk allowance per DNV-ST-N001 operation planning; schedule contingency at AACE Class 3 level',
          `Half a vessel day avoided ≈ ${fmt(0.5 * MGO_T_PER_DAY)} t of marine gas oil ≈ ${fmt(0.5 * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening factors: ${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂ per t fuel).`),
      ],
      publish: { 'flexible.GJ': GJfree, ...(lim !== null ? { 'flexible.torqueLimit': lim } : {}) },
    };
  },
  verify: [{
    // Isotropic thin tube in torsion: GJ = G·π(ro⁴ − ri⁴)/2. E = 1 GPa, ν = 0.4 (PA11 library value) → G = 357.14 MPa.
    // ri = 0.100, ro = 0.110: ro⁴ − ri⁴ = 1.4641e-4 − 1.0e-4 = 4.641e-5 → GJ = 357.14e6·π·4.641e-5/2 = 26 036 N·m².
    // Twist under 5 kN·m = 5000/26036 = 0.19204 rad/m = 11.003 deg/m.
    name: 'Polymer tube in torsion: GJ = G·π(ro⁴ − ri⁴)/2',
    source: 'Timoshenko & Goodier, Theory of Elasticity, §104 (circular shaft); hand calculation',
    inputs: { boreD: 200, Te: 0, pInt: 0, depth: 0, pAnn: 0, torque: 5, matSheath: 'pa11', layers: [{ name: 'Tube', type: 'pressureSheath', t: 10, alphaDeg: 0, lay: 0, fill: 1, wireW: 0, E: 1, yield: 0, uts: 0, lock: 0 }] },
    expect: [{ kpi: 'GJ', value: 26.036, tol: 0.002 }, { kpi: 'twist', value: 11.003, tol: 0.002 }],
  }],
};
