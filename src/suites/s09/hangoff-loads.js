// Hang-off and connection loads of a catenary riser: porch reactions, flex-joint rotation and moment,
// stress-joint sizing, connector capacity utilisation and pipe stress at the top connection.
import { riserOffsetCases, boundaryLayerHangoff } from '../../engine/lines/riser.js';
import { flexJointCheck, stressJointTaper, connectorUtilisation, vonMisesApi2rd } from '../../engine/lines/risercheck.js';
import {
  G, DEG, need, run, fmt, round, adv, check, kpi, plot, meshInput, finish, handWeight, currentInputs, LIMIT_2D, CONFIRM, CO2_STEEL,
} from './_common.js';
import { scrInputs, criteriaInputs, buildScr, currentOf } from './_scr.js';

const A12 = 985 * Math.sin(12 * DEG) / (1 - Math.sin(12 * DEG));

export default {
  id: 's09-hangoff-loads', suite: 9, title: 'Hang-off and connection loads',
  summary: 'Reactions, rotation and bending moment at the riser hang-off, with flex-joint, stress-joint and connector utilisation for near, mean and far vessel positions.',
  covers: [10], fidelity: 'numerical',
  standards: ['API RP 2RD (1998, R2013) §5', 'DNV-ST-F201 (2018) Sec. 5–6', 'API Spec 17G / ISO 13628-7 (connector capacity charts)'],
  models: ['Corotational beam riser model', 'Tensioned-beam end boundary layer M = Δθ/(1/k + 1/√(EI·T))', 'Rotational spring flex joint', 'Constant-curvature tapered stress joint', 'Connector capacity envelope (tension–moment)', 'von Mises working stress'],
  limits: `${LIMIT_2D} Vessel rotation is added as a quasi-static angle at the porch; torsion at the connection and out-of-plane moments are not produced. The connector envelope is whatever chart the user supplies — the default polygon is illustrative.`,
  publishes: ['riser.hangoffMoment', 'riser.hangoffRotationDeg', 'riser.hangoffTensionMax'],
  searchTerms: 'riser hang-off flex joint stress joint connector load',
  inputs: [
    ...scrInputs(),
    ...currentInputs(),
    { key: 'offset', label: 'Vessel offset (near / far)', unit: 'm', type: 'number', default: 40, min: 0, max: 300, group: 'Top end (hang-off)', role: 'bc', bind: 'vessel.offsetMax', uncertain: { dist: 'normal', cov: 0.2 } },
    { key: 'pitchAmp', label: 'Vessel rotation amplitude at the porch', unit: 'deg', type: 'number', default: 4, min: 0, max: 20, group: 'Top end (hang-off)', role: 'bc', uncertain: { dist: 'normal', cov: 0.15 },
      help: 'Extreme single-amplitude pitch / roll in the riser plane. It rotates the porch relative to the riser and is taken by the flex or stress joint.' },
    { key: 'daf', label: 'Dynamic amplification on top tension', unit: '-', type: 'number', default: 1.15, min: 1, max: 2.5, group: 'Top end (hang-off)', role: 'load',
      help: 'Ratio of extreme dynamic to static top tension, from the vessel-motion or coupled dynamic analysis.' },
    { key: 'envelope', label: 'Connector capacity envelope', type: 'table', group: 'Connection rating', role: 'criteria',
      columns: [{ key: 'T', label: 'Tension', unit: 'kN', type: 'number' }, { key: 'M', label: 'Bending moment', unit: 'kN·m', type: 'number' }],
      default: [{ T: 6000, M: 0 }, { T: 4000, M: 900 }, { T: 0, M: 1200 }, { T: -2500, M: 0 }, { T: 0, M: -1200 }, { T: 4000, M: -900 }],
      help: 'Vertices of the rated tension–moment envelope in order around the origin (manufacturer chart at the design pressure). At least three points enclosing the origin.' },
    { key: 'sjAllow', label: 'Allowable bending stress in a stress joint', unit: 'MPa', type: 'number', default: 250, min: 50, max: 700, group: 'Connection rating', role: 'criteria',
      help: 'Used to size an alternative tapered stress joint for the same rotation (steel ≈ 0.5–0.6 σy; titanium higher).' },
    ...criteriaInputs(),
    meshInput(160),
  ],
  mesh: { key: 'nElem', label: 'Beam elements', levels: [80, 160, 320, 640], qoi: 'tensionMax' },

  compute(inp, ctx) {
    // engineering units → SI: kN, kN·m, MPa, bar, deg handled below
    const hh = inp.depth - inp.hangoffDepth;
    need(inp.offset <= 0.25 * hh, `A vessel offset of ${inp.offset} m exceeds 25 % of the hang-off height (${fmt(hh)} m). Reduce the offset.`);
    const { pp, base, pressure } = buildScr(inp, ctx, { nElem: inp.nElem, current: currentOf(inp, 1) });
    const cases = [{ name: 'Mean', res: base }];
    if (inp.offset > 0) {
      const oc = run(() => riserOffsetCases(base, [-inp.offset, inp.offset]), 'Vessel-offset cases did not solve');
      oc.forEach((o) => need(o.result.converged, `The offset case ${o.offset} m did not converge; reduce the offset.`));
      cases.unshift({ name: `Near −${inp.offset} m`, res: oc[0].result }); cases.push({ name: `Far +${inp.offset} m`, res: oc[1].result });
    }
    const porch = base.catenary.hangoffAngleDeg, EI = pp.section.EI;
    const kRot = inp.topType === 'flex' ? inp.kFlex * 1e3 / DEG : undefined;
    const pitch = inp.pitchAmp * DEG;
    const piTop = inp.pInt * 1e5, peTop = pp.rhoW * G * inp.hangoffDepth;
    for (const c of cases) {
      const r = c.res;
      c.T = r.topTension; c.Fx = Math.abs(r.endA.Fx); c.Fz = Math.abs(r.endA.Fz); c.angle = r.hangoffAngleDeg;
      c.rotStatic = Math.abs(r.hangoffAngleDeg - porch);
      c.Mstatic = Math.abs(r.moment[0]);
      // extra end moment when the porch pitches by ±pitch relative to the riser (Sparks 2007 Ch. 13)
      c.Mpitch = inp.topType === 'pinned' || pitch === 0 ? 0 : boundaryLayerHangoff({ EI, T: c.T * inp.daf, deltaAngle: pitch, kRot }).endMoment;
      c.Td = c.T * inp.daf; c.Md = c.Mstatic + c.Mpitch;
      // joint rotation: a flex joint takes the full relative angle of a pinned-like end; a clamped end takes none
      c.rot = inp.topType === 'clamped' ? 0 : c.rotStatic + inp.pitchAmp;
      const Tw = c.Td + piTop * pp.props.Ai - peTop * pp.props.Ae;
      c.vm = vonMisesApi2rd({ T: Tw, M: c.Md, pi: piTop, pe: peTop, D: pp.D, t: pp.t, sigmaY: pp.mat.fy, Cf: inp.Cf, Ca: inp.CaAllow });
    }
    const env = (inp.envelope || []).filter((v) => Number.isFinite(v.T) && Number.isFinite(v.M)).map((v) => ({ T: v.T * 1e3, M: v.M * 1e3 }));
    need(env.length >= 3, 'The connector capacity envelope needs at least three (tension, moment) vertices enclosing the origin.');
    for (const c of cases) c.conn = run(() => connectorUtilisation({ T: c.Td, M: c.Md, envelope: env }), 'Connector envelope could not be evaluated (list the vertices in order around the origin, enclosing it)');
    const by = (f) => cases.reduce((a, c) => (f(c) > f(a) ? c : a), cases[0]);
    const gT = by((c) => c.Td), gM = by((c) => c.Md), gC = by((c) => c.conn.uc), gV = by((c) => c.vm.uc), gR = by((c) => c.rot);
    const lam = Math.sqrt(EI / gT.Td);
    const rotTotal = Math.max(...cases.map((c) => c.rotStatic)) + inp.pitchAmp;
    const sj = run(() => stressJointTaper({ deltaAngle: rotTotal * DEG, T: gT.Td, D: pp.D, t: pp.t, E: pp.mat.E, sigmaBendAllow: inp.sjAllow * 1e6 }), 'Stress-joint sizing failed');
    // steel mass of the tapered joint: integrate the wall area along the taper (constant bore)
    const Di = pp.D - 2 * pp.t; let sjMass = 0;
    for (let i = 1; i < sj.xi.length; i++) { const a = (d) => Math.PI / 4 * (d * d - Di * Di); sjMass += 0.5 * (a(sj.OD[i]) + a(sj.OD[i - 1])) * (sj.xi[i] - sj.xi[i - 1]) * pp.mat.rho; }

    const kpis = [
      kpi('tensionMax', 'Design top tension (incl. dynamic factor)', gT.Td / 1e3, 'kN', `Largest hang-off effective tension over the vessel positions multiplied by the dynamic factor (${gT.name} case). It is the axial design load of the porch, flex joint and first connector.`),
      kpi('reactionV', 'Vertical porch reaction (mean, static)', cases.find((c) => c.name === 'Mean').Fz / 1e3, 'kN', 'Static vertical load on the porch in the mean position: the suspended apparent weight of the riser.'),
      kpi('reactionH', 'Largest horizontal porch reaction (static)', Math.max(...cases.map((c) => c.Fx)) / 1e3, 'kN', 'Horizontal pull on the hull, largest in the far position. It loads the porch laterally and enters the mooring balance.'),
      kpi('angleMean', 'Riser angle at the porch (mean)', base.hangoffAngleDeg, 'deg', 'Riser axis angle from vertical at the hang-off in the mean position, to be compared with the fabricated porch angle.'),
      kpi('rotation', 'Joint rotation demand', gR.rot, 'deg', 'Largest relative rotation between riser and porch: static change with offset and current plus the vessel rotation amplitude. A clamped (stress-joint) top takes it as bending instead.', inp.topType === 'flex' ? inp.flexAllow : null, 'max'),
      kpi('moment', 'Design bending moment at the connection', gM.Md / 1e3, 'kN·m', 'Bending moment transmitted to the porch and the first connector: static moment from the solved riser plus the moment caused by vessel rotation.'),
      kpi('momentPitch', 'Moment from vessel rotation alone', gM.Mpitch / 1e3, 'kN·m', 'End moment Δθ/(1/k + 1/√(EI·T)) generated when the porch rotates by the pitch amplitude. A softer flex joint reduces it; a clamped end gives Δθ·√(EI·T).'),
      kpi('lambda', 'Bending boundary-layer length at the top', lam, 'm', 'Length √(EI/T) over which the end moment decays down the riser. The first weld below the joint usually sits inside it and sees most of the bending.'),
      kpi('connUc', 'Connector utilisation', gC.conn.uc, '-', 'Design tension and moment as a fraction of the rated envelope along the same load ray. Above 1.0 the connector is outside its rated capacity.', 1, 'max'),
      kpi('vmUc', 'Pipe stress utilisation at the top', gV.vm.uc, '-', 'von Mises stress in the plain pipe just below the joint over the API RP 2RD allowable, with the design tension, moment and pressure.', 1, 'max'),
      kpi('sjLength', 'Equivalent tapered stress-joint length', sj.length, 'm', 'Length of a constant-curvature taper that would absorb the same rotation within the allowable bending stress. Zero means plain pipe already accommodates the rotation.'),
      kpi('sjRootOD', 'Stress-joint root diameter', sj.rootOD * 1e3, 'mm', 'Outer diameter needed at the clamped root of the tapered joint; with the length it sizes the forging.'),
    ];
    const checks = [
      check('Connector capacity envelope', gC.conn.uc, 1, '-', 'API Spec 17G / ISO 13628-7 capacity chart (manufacturer rating)', 'Guards the hang-off connector against combined tension and bending beyond its rated envelope, where leakage, preload loss or structural failure of the connection can occur.'),
      check('Pipe stress at the connection (von Mises)', gV.vm.vonMises / 1e6, gV.vm.allowable / 1e6, 'MPa', 'API RP 2RD (1998) §5.2.3, Table 2', 'Guards against yielding of the riser pipe immediately below the joint under design tension, end moment and pressure.'),
    ];
    if (inp.topType === 'flex') {
      const fj = flexJointCheck({ angleDeg: gR.rot, allowableDeg: inp.flexAllow, kRot });
      checks.push(check('Flex-joint rotation', gR.rot, inp.flexAllow, 'deg', 'Flex-joint supplier rating (project data)', 'Limits the cone angle of the elastomeric joint. Beyond the rating the elastomer is over-strained and metal parts can contact, which quickly damages the joint.'));
      kpis.push(kpi('flexMoment', 'Flex-joint restoring moment at full rotation', fj.moment / 1e3, 'kN·m', 'Moment k·θ the joint would return if the whole rotation demand were taken in the elastomer: an upper bound on the moment into the porch.'));
    }

    const env2 = [...env, env[0]];
    const plots = [
      plot('Connector capacity envelope', 'Tension', 'kN', 'Bending moment', 'kN·m', [{ name: 'Rated envelope', x: env2.map((v) => v.T / 1e3), y: env2.map((v) => v.M / 1e3) }, ...cases.map((c) => ({ name: c.name, x: [0, c.Td / 1e3], y: [0, c.Md / 1e3], dash: true }))],
        'Rated tension–moment envelope with the design load point of each vessel position (end of each dashed ray). A point outside the polygon exceeds the connector rating.'),
      plot('Hang-off loads by vessel position', 'Case', '', 'Load', 'kN | kN·m', [], 'Design tension, static horizontal reaction and design moment for each vessel position; the far case governs tension, the near case usually governs rotation.', 'bar'),
      plot('Tapered stress-joint profile', 'Distance from tip towards root', 'm', 'Outer diameter', 'mm', [{ name: 'Outer diameter', x: sj.xi, y: sj.OD.map((d) => d * 1e3) }], 'Outer diameter of a constant-curvature stress joint that absorbs the rotation demand. The profile grows from the pipe diameter at the tip to the root diameter at the clamped end.'),
      plot('Bending moment below the hang-off', 'Arc length from hang-off', 'm', 'Bending moment', 'kN·m', cases.map((c) => { const xs = [], ys = []; for (let k = 0; k < c.res.s.length && c.res.s[k] <= 12 * lam; k++) { xs.push(c.res.s[k]); ys.push(c.res.moment[k] / 1e3); } return { name: c.name, x: xs, y: ys }; }),
        'Static bending moment in the first metres below the joint. It decays over the boundary-layer length; welds inside this zone carry the hang-off bending.'),
    ];
    plots[1] = { type: 'bar', title: 'Hang-off loads by vessel position', categories: cases.map((c) => c.name), x: { label: 'Vessel position', unit: '' }, y: { label: 'Load', unit: 'kN or kN·m' },
      series: [{ name: 'Design tension (kN)', y: cases.map((c) => c.Td / 1e3) }, { name: 'Horizontal reaction (kN)', y: cases.map((c) => c.Fx / 1e3) }, { name: 'Design moment (kN·m)', y: cases.map((c) => c.Md / 1e3) }], explain: plots[1].explain };
    const tables = [{ title: 'Hang-off load cases', columns: [{ key: 'case', label: 'Position' }, { key: 'T', label: 'Static tension', unit: 'kN' }, { key: 'Td', label: 'Design tension', unit: 'kN' }, { key: 'Fx', label: 'Horizontal reaction', unit: 'kN' }, { key: 'Fz', label: 'Vertical reaction', unit: 'kN' }, { key: 'ang', label: 'Riser angle', unit: 'deg' }, { key: 'rot', label: 'Joint rotation', unit: 'deg' }, { key: 'M', label: 'Design moment', unit: 'kN·m' }, { key: 'uc', label: 'Connector utilisation', unit: '-' }],
      rows: cases.map((c) => ({ case: c.name, T: round(c.T / 1e3, 5), Td: round(c.Td / 1e3, 5), Fx: round(c.Fx / 1e3, 4), Fz: round(c.Fz / 1e3, 5), ang: round(c.angle, 4), rot: round(c.rot, 3), M: round(c.Md / 1e3, 4), uc: round(c.conn.uc, 3) })),
      explain: 'Connection forces, moments and utilisation for each vessel position; these are the interface loads to give to the hull and connector designers.' }];
    const ucMax = Math.max(...checks.map((c) => c.uc)), gov = checks.find((c) => c.uc === ucMax);
    const advice = [
      adv('engineering', ucMax > 1 ? 'fail' : ucMax > 0.85 ? 'warn' : 'ok',
        `${gov.name} governs at ${fmt(ucMax)}: design tension ${fmt(gT.Td / 1e3)} kN, design moment ${fmt(gM.Md / 1e3)} kN·m, rotation demand ${fmt(gR.rot)}°.`,
        ucMax > 1 ? (gov.name.startsWith('Flex') ? `Select a joint rated for at least ±${fmt(gR.rot * 1.1)}°, or realign the porch to the mean riser angle of ${fmt(base.hangoffAngleDeg)}° to remove the static part.`
          : gov.name.startsWith('Connector') ? `Select a connector whose envelope contains ${fmt(gC.Td / 1e3)} kN with ${fmt(gC.Md / 1e3)} kN·m, or soften the top restraint to cut the moment.`
            : `Reduce the end moment (softer flex joint than ${inp.kFlex} kN·m/deg) or use an upset first joint of about ${fmt(pp.t * 1e3 * gV.vm.uc / 0.9)} mm wall.`)
          : `Set the porch angle to the mean riser angle (${fmt(base.hangoffAngleDeg)}°) so the static rotation is split evenly between near and far, and pass the moment range to the hang-off fatigue assessment.`,
        `API RP 2RD (1998, R2013) §5; connector rating per API Spec 17G charts. ${CONFIRM}`),
      adv('economic', 'info',
        `The alternative to the flex joint is a tapered stress joint ${fmt(sj.length)} m long with a ${fmt(sj.rootOD * 1e3)} mm root, about ${fmt(sjMass / 1e3)} t of forged steel. A flex joint is a long-lead elastomeric item with a finite service life; a stress joint has no elastomer but loads the porch with ${fmt(sj.rootMoment / 1e3)} kN·m.`,
        'Compare the two on whole-life terms: flex-joint replacement needs the riser to be depressurised and a vessel on site, while the stress joint adds forging mass and porch steel up front.',
        'Whole-life costing (ISO 15686-5); long-lead procurement risk per AACE RP 18R-97 estimate classes.',
        `Stress-joint forging ≈ ${fmt(sjMass / 1e3)} t of steel ≈ ${fmt(sjMass / 1e3 * CO2_STEEL)} t CO₂e at ${CO2_STEEL} t CO₂e/t (screening factor). Avoiding a mid-life flex-joint change-out avoids an offshore campaign and its fuel.`),
    ];
    return finish({
      kpis, checks, plots, tables,
      notes: [`Top restraint: ${inp.topType}${inp.topType === 'flex' ? ` (${inp.kFlex} kN·m/deg)` : ''}; porch aligned with the nominal catenary angle ${fmt(porch)}°.`, `${gov.name} governs (utilisation ${fmt(ucMax)}).`],
      advice,
      publish: { 'riser.hangoffMoment': gM.Md, 'riser.hangoffRotationDeg': gR.rot, 'riser.hangoffTensionMax': gT.Td },
      quantities: { steel_t: sjMass / 1e3 },
    });
  },

  verify: [
    {
      // Pinned top, no current, no offset: vertical reaction = suspended apparent weight w·S,
      // S = sqrt(h(h + 2a)), a = h·sinθ/(1 − sinθ); top tension T = w(h + a); DAF = 1.
      name: 'Catenary hang-off reactions: V = w·S, T = w·(h + a)', source: 'Irvine, Cable Structures (1981) §1.2',
      inputs: { topType: 'pinned', Us: 0, Ub: 0, offset: 0, pInt: 0, daf: 1, nElem: 320 },
      expect: [{ kpi: 'reactionV', value: handWeight() * Math.sqrt(985 * (985 + 2 * A12)) / 1e3, tol: 0.01 }, { kpi: 'tensionMax', value: handWeight() * (985 + A12) / 1e3, tol: 0.005 }],
    },
    {
      // Clamped top rotated by Δθ = 4°: M = Δθ·√(EI·T). EI = 207e9·π/64·(0.2731⁴ − 0.2319⁴) = 2.714e7 N·m²,
      // T = w(h + a) ≈ 1.1887e6 N → √(EI·T) = 5.680e6 N·m/rad → M = 0.06981·5.680e6 ≈ 396.5 kN·m.
      name: 'Clamped-end moment from porch rotation M = Δθ·√(EI·T)', source: 'Sparks, Fundamentals of Marine Riser Mechanics (2007) Ch. 13',
      inputs: { topType: 'clamped', Us: 0, Ub: 0, offset: 0, pInt: 0, daf: 1, pitchAmp: 4, nElem: 320 },
      expect: [{ kpi: 'momentPitch', value: 4 * DEG * Math.sqrt(207e9 * Math.PI / 64 * (0.2731 ** 4 - 0.2319 ** 4) * handWeight() * (985 + A12)) / 1e3, tol: 0.01 }],
    },
  ],
};
