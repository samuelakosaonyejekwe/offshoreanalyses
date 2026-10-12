// Structural strength of a jacket or tripod foundation for an offshore wind turbine: 3-D frame analysis under
// turbine interface loads, self-weight and buoyancy, and Morison wave and current loading, with tubular member
// code checks, pile-head reactions and deflections.
import { jacketModel, fromImportedFrame } from '../../engine/structural/builders.js';
import { analyzeFrame, refineModel, sectionStress } from '../../engine/structural/frame3d.js';
import { resolveSection } from '../../engine/structural/sections.js';
import { tubularMemberCheck } from '../../engine/structural/tubular.js';
import { airyKinematics, wavenumber } from '../../engine/metocean/waves.js';
import { currentPowerLaw } from '../../engine/metocean/current.js';
import { need, run, fmt, round, adv, check, worst, sev, steelCo2Text, steelInputs, steelOf, LIMIT_CODE, CONFIRM } from './_common.js';

/** Morison load per unit length (global vector, N/m) on a member at its mid-point for one wave phase. */
function morison(pa, pb, D, kin, rho, Cd, Cm) {
  const e = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
  const L = Math.hypot(e[0], e[1], e[2]);
  const ex = e[0] / L, ey = e[1] / L, ez = e[2] / L;
  const proj = (v) => { const d = v[0] * ex + v[1] * ey + v[2] * ez; return [v[0] - d * ex, v[1] - d * ey, v[2] - d * ez]; };
  const un = proj(kin.u), an = proj(kin.a);
  const um = Math.hypot(un[0], un[1], un[2]);
  const kd = 0.5 * rho * Cd * D * um, ki = rho * Cm * (Math.PI * D * D) / 4;
  return [kd * un[0] + ki * an[0], kd * un[1] + ki * an[1], kd * un[2] + ki * an[2]];
}

export default {
  id: 's03-structural-strength', suite: 3, title: 'Jacket and tripod structural strength',
  summary: 'Member forces, stresses and code utilisations of a jacket or tripod substructure under turbine interface loads, weight, buoyancy, waves and current, with pile-head reactions, deflections and steel mass.',
  covers: [2], fidelity: 'numerical',
  standards: ['ISO 19902:2020 §13 (tubular members)', 'DNV-ST-0126 (2021) §4 (structural design), §5', 'IEC 61400-3-1:2019 §7.6 (partial safety factors)', 'API RP 2A-WSD (21st ed.) §3.2–3.3 (comparison)', 'DNV-RP-C205 (2019) §6 (Morison loads on slender members)', 'NORSOK N-004 (2013) §6'],
  models: ['3-D frame analysis by the direct stiffness method (Timoshenko beam elements)', 'Force and moment equilibrium with rigid distribution of interface loads to the leg tops', 'Morison equation on inclined members with normal-velocity decomposition', 'Airy wave theory and power-law current', 'Buoyant self-weight (marine method)', 'ISO 19902 tubular member checks: tension, compression (column buckling), bending, shear, combined', 'von Mises equivalent stress'],
  limits: `Linear elastic frame with rigid joints and supports at the mudline (fixed, pinned or imported); joint punching-shear strength, joint flexibility, fatigue, pile–soil nonlinearity and the transition piece are not modelled. Wave kinematics are Airy without stretching, evaluated at member mid-points; the phase is chosen for maximum base shear. ${LIMIT_CODE} Shell finite-element analysis of nodes and nonlinear pushover are outside the module.`,
  searchTerms: 'offshore wind jacket tripod foundation structural analysis tubular member check ISO 19902',
  inputs: [
    { key: 'type', label: 'Substructure', type: 'select', default: 'jacket4', group: 'Geometry', role: 'geometry', options: [{ value: 'jacket4', label: 'Four-legged jacket' }, { value: 'jacket3', label: 'Three-legged jacket / tripod frame' }, { value: 'import', label: 'Imported frame model' }] },
    { key: 'frame', label: 'Imported frame (SDNF, IFC structural, beam centre-lines)', type: 'geometry', kinds: ['frame'], default: null, group: 'Geometry', role: 'geometry', showIf: { key: 'type', in: ['import'] }, help: 'Members need tube sections; the lowest nodes are supported and the highest nodes receive the interface loads. Without a file the parametric jacket is used.' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 40, min: 10, max: 90, group: 'Geometry', role: 'environment', bind: 'site.depth' },
    { key: 'zTop', label: 'Jacket top (interface) elevation', unit: 'm', type: 'number', default: 18, min: 5, max: 35, group: 'Geometry', role: 'geometry' },
    { key: 'topW', label: 'Leg spacing at the top', unit: 'm', type: 'number', default: 14, min: 4, max: 40, group: 'Geometry', role: 'geometry', design: true },
    { key: 'batter', label: 'Leg batter (1 horizontal : n vertical)', unit: '-', type: 'number', default: 10, min: 3, max: 100, group: 'Geometry', role: 'geometry', design: true },
    { key: 'bays', label: 'Bracing bays', type: 'int', default: 4, min: 1, max: 8, group: 'Geometry', role: 'geometry' },
    { key: 'bracing', label: 'Bracing pattern', type: 'select', default: 'X', group: 'Geometry', role: 'geometry', options: [{ value: 'X', label: 'X' }, { value: 'K', label: 'K' }, { value: 'diagonal', label: 'Single diagonal' }] },
    { key: 'legD', label: 'Leg diameter', unit: 'm', type: 'number', default: 1.2, min: 0.4, max: 5, group: 'Members', role: 'geometry', design: true },
    { key: 'legT', label: 'Leg wall thickness', unit: 'mm', type: 'number', default: 40, min: 10, max: 120, group: 'Members', role: 'geometry', design: true },
    { key: 'brD', label: 'Brace diameter', unit: 'm', type: 'number', default: 0.6, min: 0.2, max: 3, group: 'Members', role: 'geometry', design: true },
    { key: 'brT', label: 'Brace wall thickness', unit: 'mm', type: 'number', default: 20, min: 6, max: 80, group: 'Members', role: 'geometry', design: true },
    ...steelInputs(),
    { key: 'support', label: 'Support at the mudline', type: 'select', default: 'fixed', group: 'Boundary conditions', role: 'bc', options: [{ value: 'fixed', label: 'Fixed (piles grouted, stiff soil)' }, { value: 'pinned', label: 'Pinned' }], help: 'Idealisation of the pile or bucket heads. Pinned supports give statically clear leg forces; fixed supports add leg bending.' },
    { key: 'Fx', label: 'Horizontal interface load (characteristic)', unit: 'MN', type: 'number', default: 2.2, min: 0, max: 50, group: 'Interface loads', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Tower-base shear from the turbine load calculation (rotor thrust, tower drag, inertia).' },
    { key: 'My', label: 'Overturning interface moment (characteristic)', unit: 'MN·m', type: 'number', default: 170, min: 0, max: 3000, group: 'Interface loads', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Tower-base bending moment at the jacket top, acting with the horizontal load.' },
    { key: 'Mz', label: 'Torsion about the vertical axis', unit: 'MN·m', type: 'number', default: 12, min: 0, max: 500, group: 'Interface loads', role: 'bc', help: 'Yaw moment from the rotor; it loads the braces of all faces.' },
    { key: 'Fz', label: 'Vertical interface load', unit: 'MN', type: 'number', default: 9, min: 0, max: 100, group: 'Interface loads', role: 'bc', help: 'Weight of rotor–nacelle assembly, tower and transition piece.' },
    { key: 'dirDeg', label: 'Load direction from the x-axis', unit: 'deg', type: 'number', default: 0, min: 0, max: 90, group: 'Interface loads', role: 'bc', help: '0° loads along a face; the diagonal direction (45° for four legs) maximises the force in one leg.' },
    { key: 'Hw', label: 'Design wave height', unit: 'm', type: 'number', default: 14, min: 0, max: 35, group: 'Waves and current', role: 'bc', bind: 'site.Hmax100', uncertain: { dist: 'gumbel', cov: 0.1 } },
    { key: 'Tw', label: 'Wave period', unit: 's', type: 'number', default: 11, min: 4, max: 20, group: 'Waves and current', role: 'bc', bind: 'site.Tp100' },
    { key: 'Uc', label: 'Surface current', unit: 'm/s', type: 'number', default: 1.0, min: 0, max: 3, group: 'Waves and current', role: 'bc', bind: 'site.Uc100' },
    { key: 'Cd', label: 'Drag coefficient', unit: '-', type: 'number', default: 1.05, min: 0.5, max: 1.6, group: 'Waves and current', role: 'environment', help: 'ISO 19902: 0.65 smooth, 1.05 rough (marine growth).' },
    { key: 'Cm', label: 'Inertia coefficient', unit: '-', type: 'number', default: 1.2, min: 1, max: 2.2, group: 'Waves and current', role: 'environment', help: 'ISO 19902: 1.6 smooth, 1.2 rough.' },
    { key: 'mg', label: 'Marine growth thickness', unit: 'mm', type: 'number', default: 50, min: 0, max: 200, group: 'Waves and current', role: 'environment' },
    { key: 'selfW', label: 'Include self-weight and buoyancy', type: 'bool', default: true, group: 'Boundary conditions', role: 'load' },
    { key: 'gammaF', label: 'Load factor on environmental and interface loads', unit: '-', type: 'number', default: 1.35, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: `IEC 61400-3-1 normal design situations. ${CONFIRM}` },
    { key: 'gammaG', label: 'Load factor on weight', unit: '-', type: 'number', default: 1.1, min: 0.9, max: 1.5, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'nSub', label: 'Elements per member', type: 'int', default: 2, min: 1, max: 8, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'baseShear' },
  compute(inp, ctx) {
    const steel = steelOf(inp, ctx);
    const mat = { E: steel.E, nu: steel.nu, rho: steel.rho, fy: steel.fy };
    let base;
    if (inp.type === 'import' && inp.frame) {
      base = run(() => fromImportedFrame(inp.frame, { defaultSection: { type: 'tube', D: inp.brD, t: inp.brT / 1e3 }, defaultMaterial: mat, supports: 'lowest', supportType: inp.support }), 'Imported frame');
      base.loadCases = [{ name: 'selfWeight', selfWeight: { buoyancy: { waterline: 0 } } }];
    } else {
      base = run(() => jacketModel({ legs: inp.type === 'jacket3' ? 3 : 4, waterDepth: inp.depth, topElevation: inp.zTop, topWidthX: inp.topW, batter: inp.batter, bays: inp.bays, bracing: inp.bracing,
        legSection: { D: inp.legD, t: inp.legT / 1e3 }, braceSection: { D: inp.brD, t: inp.brT / 1e3 }, piles: { mode: inp.support }, material: mat }), 'Jacket model');
    }
    const ref = run(() => refineModel(base, inp.nSub), 'Mesh refinement');
    const model = ref.model;
    const nodes = model.nodes.map((p) => Array.from(p));
    const zMin = Math.min(...nodes.map((p) => p[2])), zMax = Math.max(...nodes.map((p) => p[2]));
    const tol = 1e-3 * Math.max(zMax - zMin, 1);
    // Interface nodes are taken from the unrefined model (refinement appends its nodes, so indices are unchanged):
    // only structural joints at the top level receive load, not the intermediate nodes of top horizontals.
    const tops = base.meta?.topNodes ?? base.nodes.map((p, i) => [p, i]).filter(([p]) => p[2] > zMax - tol).map(([, i]) => i);
    const mud = base.meta?.mudlineNodes ?? (model.supports || []).map((s) => s.node);
    need(tops.length >= 1 && mud.length >= 1, 'The frame needs at least one top node for the interface loads and one supported node.');
    const d = -zMin;
    // Interface loads distributed rigidly to the top nodes about their centroid.
    const cx = tops.reduce((a, i) => a + nodes[i][0], 0) / tops.length, cy = tops.reduce((a, i) => a + nodes[i][1], 0) / tops.length;
    const th = (inp.dirDeg * Math.PI) / 180, c = Math.cos(th), s = Math.sin(th);
    const gF = inp.gammaF, Fh = gF * inp.Fx * 1e6, Mo = gF * inp.My * 1e6, Mt = gF * inp.Mz * 1e6, Fv = inp.gammaG * inp.Fz * 1e6;
    // Moment vector of a horizontal force along (c, s) at height: about the horizontal axis (−s, c).
    const Mx = -Mo * s, Myy = Mo * c;
    let Sxx = 0, Syy = 0, Srr = 0;
    for (const i of tops) { const x = nodes[i][0] - cx, y = nodes[i][1] - cy; Sxx += x * x; Syy += y * y; Srr += x * x + y * y; }
    const nodal = tops.map((i) => {
      const x = nodes[i][0] - cx, y = nodes[i][1] - cy;
      const fz = -Fv / tops.length + (Syy > 0 ? (Mx * y) / Syy : 0) - (Sxx > 0 ? (Myy * x) / Sxx : 0);
      return { node: i, F: [(Fh * c) / tops.length - (Srr > 0 ? (Mt * y) / Srr : 0), (Fh * s) / tops.length + (Srr > 0 ? (Mt * x) / Srr : 0), fz, 0, 0, 0] };
    });
    // Morison loads on wet members: phase scan for the largest base shear along the load direction.
    const H = inp.Hw, T = inp.Tw, k = H > 0 ? wavenumber({ T, d }) : 0, om = (2 * Math.PI) / T;
    const mems = model.members.map((m, i) => {
      const pa = nodes[m.a], pb = nodes[m.b];
      const sec = m.section && m.section.D ? m.section : null;
      return { i, pa, pb, mid: [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])], L: Math.hypot(pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]), D: sec ? sec.D : 0, t: sec ? sec.t : 0, group: m.group || 'member' };
    });
    // Wet part of each member: clipped at the still-water level so that refinement does not move the loaded length.
    const wet = [];
    for (const m of mems) {
      if (!(m.D > 0)) continue;
      const za = m.pa[2], zb = m.pb[2];
      if (za >= 0 && zb >= 0) continue;
      let x1 = 0, x2 = 1;
      if (za < 0 && zb > 0) x2 = -za / (zb - za); else if (zb < 0 && za > 0) x1 = za / (za - zb);
      const xm = 0.5 * (x1 + x2);
      wet.push({ ...m, x1, x2, Lw: (x2 - x1) * m.L, mid: [m.pa[0] + xm * (m.pb[0] - m.pa[0]), m.pa[1] + xm * (m.pb[1] - m.pa[1]), Math.min(m.pa[2] + xm * (m.pb[2] - m.pa[2]), 0)] });
    }
    const loadsAt = (ph) => wet.map((m) => {
      const xs = m.mid[0] * c + m.mid[1] * s; // coordinate along the wave direction
      const q = H > 0 ? airyKinematics({ H, T, d, x: xs, z: m.mid[2], t: ph / om, k, clipAboveSurface: false }) : { u: 0, w: 0, ax: 0, az: 0 };
      const uc = inp.Uc > 0 ? currentPowerLaw(m.mid[2], { v0: inp.Uc, d }) : 0;
      const uh = q.u + uc;
      return { m, f: morison(m.pa, m.pb, m.D + (2 * inp.mg) / 1e3, { u: [uh * c, uh * s, q.w], a: [q.ax * c, q.ax * s, q.az] }, 1025, inp.Cd, inp.Cm) };
    });
    let best = { shear: -Infinity, ph: 0, loads: [] };
    const nPh = H > 0 ? 24 : 1;
    for (let j = 0; j < nPh; j++) {
      const ph = (2 * Math.PI * j) / nPh;
      const ls = loadsAt(ph);
      const sh = ls.reduce((a, l) => a + (l.f[0] * c + l.f[1] * s) * l.m.Lw, 0);
      if (sh > best.shear) best = { shear: sh, ph, loads: ls };
    }
    const waveShear = Math.max(best.shear, 0);
    const waveOtm = best.loads.reduce((a, l) => a + (l.f[0] * c + l.f[1] * s) * l.m.Lw * (l.m.mid[2] - zMin), 0);
    model.loadCases = [
      ...(inp.selfW ? [{ name: 'weight', selfWeight: { buoyancy: { waterline: 0 } } }] : []),
      { name: 'interface', nodal },
      { name: 'wave', member: best.loads.map((l) => ({ member: l.m.i, type: 'uniform', q: l.f.map((v) => v * gF), axes: 'global', xi1: l.m.x1, xi2: l.m.x2 })) },
    ];
    model.combinations = [{ name: 'ULS', factors: { ...(inp.selfW ? { weight: inp.gammaG } : {}), interface: 1, wave: 1 } }];
    const an = run(() => analyzeFrame(model, { cases: ['ULS'] }), 'Frame analysis');
    const res = an.results.ULS;
    need(res.equilibrium.relative < 1e-6, 'The frame solution is not in equilibrium (mechanism or ill-conditioned model): check supports and connectivity.');
    // Member checks (ISO 19902) at both ends of every element, with the unbraced length of the parent member.
    const parentLen = new Map();
    ref.memberMap.forEach((subs, p) => { const m0 = base.members[p]; const a = base.nodes[m0.a], b = base.nodes[m0.b]; const L = Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]); subs.forEach((q) => parentLen.set(q, { L, p })); });
    const byParent = new Map();
    let steelKg = 0, vmMax = 0;
    for (const mf of res.memberForces) {
      const m = mems[mf.member];
      if (!m.D) continue;
      steelKg += Math.PI * (m.D - m.t) * m.t * m.L * steel.rho;
      const pl = parentLen.get(mf.member) ?? { L: m.L, p: mf.member };
      for (const end of [mf.endA, mf.endB]) {
        const ck = tubularMemberCheck({ code: 'ISO', section: { D: m.D, t: m.t }, fy: steel.fy, forces: end, L: pl.L, options: { E: steel.E } });
        const vm = sectionStress({ section: resolveSection({ type: 'tube', D: m.D, t: m.t }), ...end }).vonMises;
        vmMax = Math.max(vmMax, vm);
        const cur = byParent.get(pl.p);
        // utilisations beyond the elastic critical load come back unbounded: report them capped
        const ucm = Number.isFinite(ck.governing.uc) ? Math.min(ck.governing.uc, 99) : 99;
        if (!cur || ucm > cur.uc) byParent.set(pl.p, { p: pl.p, uc: ucm, mode: ck.governing.name, ref: ck.governing.ref, N: end.N, M: Math.hypot(end.My, end.Mz), V: Math.hypot(end.Vy, end.Vz), group: m.group, D: m.D, t: m.t, L: pl.L, vm });
      }
    }
    const list = [...byParent.values()].sort((a, b) => b.uc - a.uc);
    need(list.length > 0, 'No tubular members were found in the model: member sections must be tubes with diameter and thickness.');
    const gov = list[0];
    const grp = (g) => list.filter((q) => q.group === g)[0];
    const legW = grp('leg') ?? gov, brW = list.filter((q) => q.group !== 'leg')[0] ?? gov;
    // Reactions: pile axial (compression positive downward on the soil) and shear.
    const reac = res.reactions.filter((r) => mud.includes(r.node));
    const axial = reac.map((r) => r.R[2]), shearR = reac.map((r) => Math.hypot(r.R[0], r.R[1]));
    const compMax = Math.max(...axial, 0), tensMax = Math.max(...axial.map((v) => -v), 0);
    const u = res.displacements;
    const disp = (i) => Math.hypot(u[6 * i], u[6 * i + 1], u[6 * i + 2]);
    const topDisp = Math.max(...tops.map(disp));
    const baseShear = Math.hypot(res.reactionSum[0], res.reactionSum[1]);
    const otm = Fh * (zMax - zMin) + Mo + gF * waveOtm;
    const fyd = steel.fy / 1.1;
    const checks = [
      check(`Most utilised leg member (${legW.mode})`, legW.uc, 1, 'utilisation', `${legW.ref}; ISO 19902:2020 §13`, 'Guards the legs, which carry the overturning moment as axial push–pull, against yielding or column buckling between bracing levels. Failure of a leg is not redundant.'),
      check(`Most utilised brace or horizontal (${brW.mode})`, brW.uc, 1, 'utilisation', `${brW.ref}; ISO 19902:2020 §13`, 'Guards the bracing, which carries shear and torsion, against buckling. A buckled brace sheds load to its neighbours and can start progressive collapse.'),
      check('Largest von Mises stress', vmMax / 1e6, fyd / 1e6, 'MPa', 'DNV-ST-0126 (2021) §4.4 (γM = 1.1)', 'Equivalent stress from axial force, bending, shear and torsion in the most stressed section. Exceeding design yield means local plastic deformation at that section.'),
    ];
    const g2 = worst(checks), ok = g2.uc <= 1;
    // Thickness that brings the governing member to unit utilisation (stress ∝ 1/t for a thin tube; buckling strength changes little).
    const tReq = gov.t * Math.min(gov.uc, 4), groupLen = mems.filter((m) => m.group === gov.group).reduce((a, m) => a + m.L, 0);
    const dMass = Math.PI * gov.D * (tReq - gov.t) * groupLen * steel.rho / 1e3;
    const scale = topDisp > 0 ? (0.05 * (zMax - zMin)) / topDisp : 0;
    const lines = model.members.map((m) => [m.a, m.b]);
    // nodal value: utilisation of the most utilised member framing into the node
    const nodeUc = new Array(nodes.length).fill(0);
    ref.memberMap.forEach((subs, p) => { const q = byParent.get(p); if (!q) return; subs.forEach((mi) => { const m = model.members[mi]; nodeUc[m.a] = Math.max(nodeUc[m.a], q.uc); nodeUc[m.b] = Math.max(nodeUc[m.b], q.uc); }); });
    return {
      kpis: [
        { key: 'ucMax', label: 'Largest member utilisation', unit: '-', value: gov.uc, limit: 1, sense: 'max', explain: `Governing ISO 19902 check over all members (${gov.group}, ${gov.mode.toLowerCase()}). Must not exceed 1 under factored loads.` },
        { key: 'ucLeg', label: 'Largest leg utilisation', unit: '-', value: legW.uc, explain: 'Most utilised leg segment: axial compression with bending, including column buckling over the bay height.' },
        { key: 'ucBrace', label: 'Largest brace utilisation', unit: '-', value: brW.uc, explain: 'Most utilised brace or horizontal. Braces usually govern in torsion-dominated or wave-dominated cases.' },
        { key: 'baseShear', label: 'Factored base shear', unit: 'MN', value: baseShear / 1e6, explain: 'Sum of the horizontal support reactions: interface shear plus wave and current load on the substructure, with load factors.' },
        { key: 'otm', label: 'Factored overturning moment at the mudline', unit: 'MN·m', value: otm / 1e6, explain: 'Interface moment plus horizontal loads times their height above the mudline. The legs resist it as a push–pull couple over the base width.' },
        { key: 'waveShear', label: 'Wave and current force on the substructure (characteristic)', unit: 'MN', value: waveShear / 1e6, explain: `Morison force on all wet members at the wave phase of maximum base shear (${fmt(best.ph * 180 / Math.PI)}°). Jackets are transparent: this is far smaller than on a monopile in the same wave.` },
        { key: 'pileCompMax', label: 'Largest pile-head compression', unit: 'MN', value: compMax / 1e6, explain: 'Largest vertical support reaction pushing a pile or bucket into the soil: the design load for axial pile capacity.' },
        { key: 'pileTensMax', label: 'Largest pile-head tension', unit: 'MN', value: tensMax / 1e6, explain: 'Largest uplift on a support. Sustained tension is critical for suction buckets in sand and reduces pile fatigue life.' },
        { key: 'pileShearMax', label: 'Largest pile-head shear', unit: 'MN', value: Math.max(...shearR, 0) / 1e6, explain: 'Largest horizontal support reaction, for the lateral pile or bucket check.' },
        { key: 'topDisp', label: 'Displacement at the interface', unit: 'mm', value: topDisp * 1e3, explain: 'Largest displacement of the jacket top under factored loads; it adds to the tower-top motion and indicates the substructure stiffness.' },
        { key: 'vmMax', label: 'Largest von Mises stress', unit: 'MPa', value: vmMax / 1e6, limit: fyd / 1e6, sense: 'max', explain: 'Peak equivalent stress in any member section under factored loads.' },
        { key: 'steelT', label: 'Primary steel mass', unit: 't', value: steelKg / 1e3, explain: 'Mass of all tubular members in the model, without nodes reinforcement, transition piece, piles and secondary steel (add 25–40 % for a complete jacket).' },
      ],
      checks,
      plots: [
        { type: 'model3d', title: 'Member utilisation and deformed shape', nodes, lines, nodeValues: nodeUc, valueLabel: 'Member utilisation', valueUnit: '-', deformed: nodes.map((p, i) => [p[0] + scale * u[6 * i], p[1] + scale * u[6 * i + 1], p[2] + scale * u[6 * i + 2]]),
          explain: 'Frame model coloured by the utilisation of the members at each node, with the deformed shape scaled to 5 % of the height. Look for a concentration of high utilisation in one leg or bay.' },
        { type: 'bar', title: 'Utilisation of the fifteen most loaded members', categories: list.slice(0, 15).map((q) => `${q.group} ${q.p}`), series: [{ name: 'Utilisation', y: list.slice(0, 15).map((q) => q.uc) }], x: { label: 'Member' }, y: { label: 'Utilisation', unit: '-' },
          explain: 'Ranking of members by code utilisation. A steep fall after the first few bars means local strengthening is efficient; a flat ranking means the whole frame is evenly used.' },
        { type: 'bar', title: 'Vertical reactions at the supports', categories: reac.map((r, i) => `Support ${i + 1}`), series: [{ name: 'Compression (+) / tension (−), MN', y: axial.map((v) => v / 1e6) }], x: { label: 'Support' }, y: { label: 'Vertical reaction', unit: 'MN' },
          explain: 'Push–pull of the foundation under the overturning moment, including weight. Negative bars are supports in tension.' },
      ],
      tables: [
        { title: 'Most utilised members', columns: [{ key: 'm', label: 'Member' }, { key: 'g', label: 'Group' }, { key: 'D', label: 'D', unit: 'm' }, { key: 't', label: 't', unit: 'mm' }, { key: 'L', label: 'Length', unit: 'm' }, { key: 'N', label: 'Axial (tension +)', unit: 'MN' }, { key: 'M', label: 'Moment', unit: 'MN·m' }, { key: 'uc', label: 'Utilisation' }, { key: 'mode', label: 'Governing check' }],
          rows: list.slice(0, 20).map((q) => ({ m: q.p, g: q.group, D: round(q.D), t: round(q.t * 1e3), L: round(q.L), N: round(q.N / 1e6), M: round(q.M / 1e6), uc: round(q.uc), mode: q.mode })),
          explain: 'Twenty most utilised members with their governing section forces and failure mode: the list to work through when resizing.' },
        { title: 'Support reactions (factored)', columns: [{ key: 's', label: 'Support' }, { key: 'x', label: 'x', unit: 'm' }, { key: 'y', label: 'y', unit: 'm' }, { key: 'Rz', label: 'Vertical', unit: 'MN' }, { key: 'Rh', label: 'Horizontal', unit: 'MN' }, { key: 'Mr', label: 'Moment', unit: 'MN·m' }],
          rows: reac.map((r, i) => ({ s: i + 1, x: round(nodes[r.node][0]), y: round(nodes[r.node][1]), Rz: round(r.R[2] / 1e6), Rh: round(shearR[i] / 1e6), Mr: round(Math.hypot(r.R[3], r.R[4]) / 1e6) })),
          explain: 'Forces delivered to each pile or bucket head. Pass the largest compression, tension and shear to the foundation capacity checks.' },
      ],
      notes: [`${model.members.length} elements, ${nodes.length} nodes; ${inp.type === 'import' && inp.frame ? 'imported frame' : `${inp.type === 'jacket3' ? 'three' : 'four'}-legged, ${inp.bays} bays of ${inp.bracing} bracing, base width ${fmt(base.meta.baseWidth[0])} m`}.`,
        `Governing: ${gov.group} member ${gov.p}, ${gov.mode.toLowerCase()} at ${fmt(gov.uc * 100)} %.`, `Wave and current add ${fmt(gF * waveShear / Math.max(baseShear, 1) * 100)} % of the base shear; the interface loads the rest.`,
        'Tubular joints (punching shear, chord plastification) and fatigue need separate checks.'],
      advice: [
        ok ? adv('engineering', sev(g2.uc), `All members pass: legs ${fmt(legW.uc * 100)} %, braces ${fmt(brW.uc * 100)} %, peak stress ${fmt(vmMax / 1e6)} MPa. Pile-head loads: ${fmt(compMax / 1e6)} MN compression, ${fmt(tensMax / 1e6)} MN tension, ${fmt(Math.max(...shearR) / 1e6)} MN shear.`,
          `Run the diagonal load direction and the other design situations, then check joints and fatigue; the ${gov.group} group could be thinned towards ${fmt(Math.max(tReq, 0.008) * 1e3)} mm (from ${fmt(gov.t * 1e3)} mm) on strength alone.`,
          'ISO 19902:2020 §13–14; DNV-ST-0126 (2021) §4; IEC 61400-3-1:2019 §7.6',
          `Resizing the ${gov.group} group to unit utilisation changes steel by ${fmt(dMass)} t — ${steelCo2Text(dMass)} Jacket mass is usually set by fatigue at the joints, so confirm that before removing steel.`)
          : adv('engineering', 'fail', `${g2.name} fails at ${fmt(g2.uc * 100)} % (member ${gov.p}, ${gov.group}, D = ${fmt(gov.D)} m, t = ${fmt(gov.t * 1e3)} mm).`,
            gov.mode.toLowerCase().includes('compression') || gov.mode.toLowerCase().includes('buckling') ? `Increase the ${gov.group} wall to about ${fmt(tReq * 1.03 * 1e3)} mm or its diameter (buckling strength grows with D), or add a bracing level to shorten the ${fmt(gov.L)} m unbraced length.` : `Increase the ${gov.group} wall to about ${fmt(tReq * 1.03 * 1e3)} mm, widen the base (leg force ∝ 1/base width) or use a higher steel grade.`,
            'ISO 19902:2020 §13.2–13.4; DNV-ST-0126 (2021) §4.4',
            `Thickening the ${gov.group} group adds ${fmt(dMass)} t — ${steelCo2Text(dMass)} A wider footprint reduces leg forces with less added steel but increases bracing length.`),
        adv('economic', 'info', `The frame holds ${fmt(steelKg / 1e3)} t of primary steel for a factored overturning moment of ${fmt(otm / 1e6)} MN·m (${fmt(steelKg / 1e3 / Math.max(otm / 1e6, 1e-9))} t per MN·m).`,
          'Jacket cost follows the number of joints and welding hours more than tonnage: fewer bays, standardised brace sizes and automated-weld node geometry reduce fabrication cost even where they add some steel. Compare with a monopile on installed cost including piling or buckets.',
          'Fabrication cost models for jackets (man-hours per joint and per tonne); AACE RP 18R-97 Class 4 estimate'),
      ],
      quantities: { steel_t: round(steelKg / 1e3) },
    };
  },
  verify: [
    { name: 'Pinned four-legged jacket under a horizontal top load: leg reaction = F·H / (2·B)', source: 'Moment equilibrium about the mudline with antisymmetric vertical reactions (statics): F = 1 MN at the top, H = 40 + 18 = 58 m, base width B = 14 + 2 × 58/10 = 25.6 m → R = 1 × 58 / (2 × 25.6) = 1.1328 MN',
      inputs: { type: 'jacket4', support: 'pinned', selfW: false, Hw: 0, Uc: 0, Fx: 1, My: 0, Mz: 0, Fz: 0, gammaF: 1, depth: 40, zTop: 18, topW: 14, batter: 10 }, expect: [{ kpi: 'pileCompMax', value: 1.1328, tol: 0.001 }, { kpi: 'baseShear', value: 1.0, tol: 1e-6 }] },
    { name: 'Interface moment only: leg reaction = M / (2·B)', source: 'Statics of the leg couple: M = 100 MN·m about the horizontal axis, base width 25.6 m, two legs per side → R = 100 / (2 × 25.6) = 1.9531 MN', inputs: { type: 'jacket4', support: 'pinned', selfW: false, Hw: 0, Uc: 0, Fx: 0, My: 100, Mz: 0, Fz: 0, gammaF: 1 }, expect: [{ kpi: 'pileCompMax', value: 1.9531, tol: 0.001 }] },
  ],
};
