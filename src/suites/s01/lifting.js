// Lifting and installation: four-point (or as-imported) lift of the upright jacket with slings to a single hook,
// sling and padeye design loads with dynamic and skew factors, member utilisation in the lift, hook load during
// lowering through the water line, hoist-line resonance and the on-bottom stability of the unpiled jacket.
import { G, DEG, need, run, fmt, round, adv, check, fsCheck, linspace, steelCo2, vesselCo2, CODE_CAVEAT, MORISON_CAVEAT, CONFIRM, BASIS_MEMBER, BASIS_COST } from './_common.js';
import { geomInputs, envInputs, criteriaInputs, meshInputs, buildJacket, designWave, solve, memberRows, byGroup, displacements, nodeValuesFromMembers, plot3d, memberCols, memberTableRows, groupChecks, testFrame } from './_model.js';
import { padeyeCheck } from '../../engine/structural/steelcheck.js';
import { axialResonance, heaveResponse } from '../../engine/lines/lowering.js';

const crit = criteriaInputs({ gG: 1.3 }).filter((i) => i.key !== 'gE');
const envI = envInputs({ H: 3.5, T: 7.5, Uc: 0.5, Uw: 15, surge: 0.5, Hlabel: 'Installation wave height, unpiled phase (maximum individual wave)', Hbind: null, Tbind: 'site.Tp1', Ucbind: 'site.Uc1', Uwbind: 'site.Uw1', group: 'Installation environment' });

export default {
  id: 's01-lifting', suite: 1, title: 'Lifting and installation',
  summary: 'Hook load, sling forces and padeye utilisation for the offshore lift of the jacket, member utilisation while hanging, hook load and line dynamics during lowering, and sliding, overturning and mudmat bearing of the jacket standing unpiled on the seabed.',
  covers: [17], fidelity: 'numerical',
  standards: ['DNV-ST-N001 (2021) §16 (lifting operations: DAF, skew load, sling and lift-point factors)', 'ISO 19901-6:2009 §18–§19 (lifting, installation)', 'DNV-RP-N103 (2019) §5 (deep-water lowering)', 'API RP 2A-WSD 22nd ed. §5.4.3, §12', 'ISO 19902:2020 §8.3 and §22 (installation, on-bottom stability)'],
  models: ['Static equilibrium of the lifted body on sling (bar) elements with a single hook', 'Dynamic amplification, skew load and weight contingency factors', 'Padeye checks by section statics (bearing, tear-out, net tension, base)', 'Buoyancy of members during lowering (marine method)', 'Hoist-line axial resonance and heave response (single degree of freedom)', 'Morison loading and Coulomb sliding / overturning / bearing of the unpiled jacket'],
  limits: `Upright single-hook lift with slings attached at the leg tops; spreader bars, two-crane lifts, launching from a barge and upending are not simulated (upending is a sequence of static positions with changing buoyancy that needs the model in each position). Lowering dynamics are a single-degree-of-freedom estimate below the wave zone, and mudmat bearing is a uniform-strength undrained check. Lift factors are criteria to be confirmed against the governing edition. ${MORISON_CAVEAT} ${CODE_CAVEAT}`,
  searchTerms: 'jacket lift installation sling padeye hook load on-bottom stability mudmat',
  inputs: [
    ...geomInputs(),
    { key: 'swCont', label: 'Weight contingency', unit: '-', type: 'number', default: 1.1, min: 1, max: 1.5, group: 'Lift weight and factors', role: 'load', help: 'On modelled steel for cans, anodes, mudmats, installation aids and weighing tolerance.' },
    { key: 'rigW', label: 'Rigging weight', unit: 't', type: 'number', default: 45, min: 0, max: 2000, group: 'Lift weight and factors', role: 'load', help: 'Slings, shackles and spreaders below the hook; about 3–5 % of the lift weight.' },
    { key: 'DAF', label: 'Dynamic amplification factor', unit: '-', type: 'number', default: 1.15, min: 1, max: 2, group: 'Lift weight and factors', role: 'criteria', help: `Offshore lift in air by DNV-ST-N001: about 1.25 for 100–300 t, 1.20 for 300–1000 t, 1.15 for 1000–2500 t, 1.10 above. ${CONFIRM}` },
    { key: 'SKL', label: 'Skew load factor', unit: '-', type: 'number', default: 1.25, min: 1, max: 2, group: 'Lift weight and factors', role: 'criteria', help: '1.25 for a statically indeterminate four-sling lift with matched slings; 1.0 for a determinate arrangement.' },
    { key: 'gC', label: 'Consequence factor for lift points and members framing into them', unit: '-', type: 'number', default: 1.3, min: 1, max: 1.6, group: 'Lift weight and factors', role: 'criteria', help: `1.3 for lift points and attachments, 1.15 for main members, 1.0 otherwise (DNV-ST-N001). Applied to the padeye; members use the action factor. ${CONFIRM}` },
    { key: 'slingAng', label: 'Sling angle from horizontal', unit: 'deg', type: 'number', default: 60, min: 30, max: 89, group: 'Rigging', role: 'geometry', design: true, help: 'Sets the hook height above the leg tops; 60° is the usual minimum.' },
    { key: 'slingD', label: 'Sling diameter', unit: 'mm', type: 'number', default: 180, min: 10, max: 500, group: 'Rigging', role: 'geometry' },
    { key: 'slingMBL', label: 'Sling minimum breaking load', unit: 't', type: 'number', default: 2200, min: 1, max: 20000, group: 'Rigging', role: 'material', design: true },
    { key: 'slingSF', label: 'Nominal sling safety factor', unit: '-', type: 'number', default: 3.0, min: 2, max: 6, group: 'Rigging', role: 'criteria', help: `Product of load, consequence, reduction, wear and material factors; not less than about 3.0 for steel slings. ${CONFIRM}` },
    { key: 'craneCap', label: 'Crane capacity at the lift radius', unit: 't', type: 'number', default: 2500, min: 10, max: 30000, group: 'Rigging', role: 'operation' },
    { key: 'peT', label: 'Padeye main plate thickness', unit: 'mm', type: 'number', default: 100, min: 10, max: 300, group: 'Padeye', role: 'geometry', design: true },
    { key: 'peR', label: 'Padeye main plate radius', unit: 'mm', type: 'number', default: 330, min: 30, max: 1000, group: 'Padeye', role: 'geometry' },
    { key: 'peHole', label: 'Pin hole diameter', unit: 'mm', type: 'number', default: 186, min: 10, max: 600, group: 'Padeye', role: 'geometry' },
    { key: 'pePin', label: 'Shackle pin diameter', unit: 'mm', type: 'number', default: 180, min: 10, max: 600, group: 'Padeye', role: 'geometry' },
    { key: 'peCheekT', label: 'Cheek plate thickness (each side)', unit: 'mm', type: 'number', default: 50, min: 0, max: 150, group: 'Padeye', role: 'geometry' },
    { key: 'peCheekR', label: 'Cheek plate radius', unit: 'mm', type: 'number', default: 270, min: 0, max: 900, group: 'Padeye', role: 'geometry' },
    { key: 'peBase', label: 'Padeye base width', unit: 'mm', type: 'number', default: 1100, min: 50, max: 3000, group: 'Padeye', role: 'geometry' },
    { key: 'peH', label: 'Hole centre above the base', unit: 'mm', type: 'number', default: 400, min: 30, max: 1500, group: 'Padeye', role: 'geometry' },
    { key: 'wireEA', label: 'Crane hoist line axial stiffness EA (all falls)', unit: 'MN', type: 'number', default: 6000, min: 10, max: 1e5, group: 'Lowering', role: 'bc' },
    { key: 'tipHeave', label: 'Crane tip heave amplitude', unit: 'm', type: 'number', default: 0.5, min: 0, max: 5, group: 'Lowering', role: 'bc', help: 'Single amplitude at the crane tip in the installation sea state.' },
    { key: 'legsFlooded', label: 'Legs flood freely during lowering', type: 'bool', default: true, group: 'Lowering', role: 'ic', help: 'Sealed legs give buoyancy that unloads the hook and must be ballasted to set the jacket down.' },
    ...envI.filter((i) => !['theory'].includes(i.key)), envI.find((i) => i.key === 'theory'),
    { key: 'mudA', label: 'Mudmat area per leg', unit: 'm²', type: 'number', default: 80, min: 1, max: 1000, group: 'On-bottom stability (unpiled)', role: 'geometry', design: true },
    { key: 'suMud', label: 'Undrained shear strength at the mudline', unit: 'kPa', type: 'number', default: 20, min: 1, max: 500, group: 'On-bottom stability (unpiled)', role: 'soil', bind: 'site.su0' },
    { key: 'muSoil', label: 'Mudmat sliding coefficient', unit: '-', type: 'number', default: 0.3, min: 0.05, max: 1, group: 'On-bottom stability (unpiled)', role: 'soil', help: 'Friction or adhesion ratio between mudmat and seabed; with skirts the soil strength governs instead.' },
    { key: 'fsStab', label: 'Required factor of safety, unpiled', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'On-bottom stability (unpiled)', role: 'criteria', help: `1.5 against sliding and overturning and 2.0 on bearing (factor × 1.33) are usual installation values (API RP 2A / ISO 19902 §22). ${CONFIRM}` },
    { key: 'liftDays', label: 'Heavy-lift vessel days on site', unit: 'days', type: 'number', default: 3, min: 0.2, max: 60, group: 'Rigging', role: 'operation' },
    ...crit, ...meshInputs({ nSub: 2, nPhase: 16 }),
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'maxDisp' },
  compute(inp, ctx) {
    // SI conversion: t → kg → N, mm → m, MN → N, kPa → Pa.
    const J0 = buildJacket({ ...inp, found: 'fixed' }, ctx, { found: 'fixed' }), M = J0.model, P = M.nodes, nT = J0.topNodes.length;
    const dry = () => M.members.forEach((m) => { m.flooded = false; });
    // hook above the centre of gravity (plan), slings from the leg tops
    const mass = J0.orig.map((m) => m.sec.A * m.L * J0.steel.rho), mt = mass.reduce((s, v) => s + v, 0), cg = [0, 1].map((k) => J0.orig.reduce((s, m, i) => s + mass[i] * 0.5 * (J0.base.nodes[m.a][k] + J0.base.nodes[m.b][k]), 0) / mt);
    const rad = Math.max(...J0.topNodes.map((n) => Math.hypot(P[n][0] - cg[0], P[n][1] - cg[1]))), hHook = Math.max(rad * Math.tan(inp.slingAng * DEG), 5), hook = P.length;
    P.push([cg[0], cg[1], J0.zTop + hHook]);
    // slings as bar elements: metallic area about 0.45 of the gross circle, E 100 GPa (wire rope / cable-laid sling)
    const As = 0.45 * Math.PI * (inp.slingD / 1e3) ** 2 / 4;
    M.materials.sling = { E: 1.0e11, nu: 0.3, rho: 0, fy: 1.5e9 };
    const s0 = M.members.length;
    J0.topNodes.forEach((n) => M.members.push({ a: n, b: hook, section: { A: As, Iy: 1e-8, Iz: 1e-8, J: 1e-8 }, material: 'sling', type: 'truss', group: 'sling' }));
    // weak horizontal springs at the leg bases and tops stop the free pendulum and yaw modes; the hook is above the
    // centre of gravity, so they take no load
    const kS = Math.max(1e4, mt * G / 20), soft = [kS, kS, 0, 0, 0, kS];
    const lift = () => { M.supports = [{ node: hook, fix: 'pinned' }, ...J0.mudNodes.map((n) => ({ node: n, spring: soft })), ...J0.topNodes.filter((n) => !J0.mudNodes.includes(n)).map((n) => ({ node: n, spring: [kS, kS, 0, 0, 0, 0] }))]; };
    lift(); dry();
    const fStat = inp.swCont, fDyn = inp.swCont * inp.DAF, fMem = inp.swCont * inp.DAF * inp.gG;
    const cases = [{ name: 'static', nodal: [], member: [], selfWeight: { factor: fStat } }, { name: 'member', nodal: [], member: [], selfWeight: { factor: fMem } }];
    // lowering: water line rising up the jacket (z measured in the in-place datum)
    const wl = linspace(J0.zMud, J0.zTop, 7);
    wl.forEach((z, i) => cases.push({ name: `wl${i}`, nodal: [], member: [], selfWeight: { factor: fStat, buoyancy: { waterline: z, rhoW: inp.rhoW } } }));
    const flood = (on) => M.members.forEach((m, e) => { m.flooded = on && e < s0 && J0.orig[J0.elems[e]?.o]?.flooded; });
    const res = run(() => solve(J0, cases.slice(0, 2), 'Lift analysis'), 'Lift analysis');
    const slingN = (r) => J0.topNodes.map((n, k) => r.memberForces[s0 + k].endA.N), sStat = slingN(res.static), W = sStat.reduce((s, v, k) => s + v * (P[hook][2] - P[J0.topNodes[k]][2]) / Math.hypot(...[0, 1, 2].map((q) => P[hook][q] - P[J0.topNodes[k]][q])), 0);
    need(W > 0 && sStat.every((v) => Number.isFinite(v)), 'The lift model did not carry the weight: check that the top nodes are connected to the structure.');
    flood(inp.legsFlooded);
    const resW = solve(J0, cases.slice(2), 'Lowering analysis'), hookW = wl.map((z, i) => { const r = resW[`wl${i}`]; return slingN(r).reduce((s, v, k) => s + v * (P[hook][2] - P[J0.topNodes[k]][2]) / Math.hypot(...[0, 1, 2].map((q) => P[hook][q] - P[J0.topNodes[k]][q])), 0); });
    dry();
    const rig = inp.rigW * 1e3 * G, hookStat = W + rig, hookDyn = hookStat * inp.DAF, sMax = Math.max(...sStat), slingDes = sMax * inp.DAF * inp.SKL, slingCap = inp.slingMBL * 1e3 * G / inp.slingSF;
    // member check in the lift (sling elements excluded)
    const rows = memberRows(J0, res.member, inp, { name: 'lift' }), groups = byGroup(rows), gov = rows.reduce((a, r) => (r.uc > a.uc ? r : a)), dp = displacements(J0, res.member);
    // padeye at the most loaded sling; the sling leaves the padeye plane by the given out-of-plane share
    const pe = run(() => padeyeCheck({ F: slingDes * inp.gC * inp.gG, fy: J0.steel.fy, t: inp.peT / 1e3, R: inp.peR / 1e3, dHole: inp.peHole / 1e3, dPin: inp.pePin / 1e3, tCheek: inp.peCheekT / 1e3, rCheek: inp.peCheekR / 1e3, baseWidth: inp.peBase / 1e3, holeHeight: inp.peH / 1e3, angleDeg: 90 - inp.slingAng, outOfPlaneFraction: 0.05, allow: { bearing: 0.9, shear: 0.58, tension: 0.9, combined: 0.9, weld: 0.58 } }), 'Padeye check (check hole, pin and plate dimensions)');
    need(inp.pePin < inp.peHole && inp.peHole < 2 * inp.peR, 'Padeye: the pin must fit the hole and the hole must be smaller than the plate.');
    // lowering dynamics near the seabed: full submerged weight on the hoist line
    const depth = inp.depth, Wsub = Math.max(hookW[hookW.length - 1], 0), hor = (m) => Math.sqrt(Math.max(0, 1 - ((J0.base.nodes[m.b][2] - J0.base.nodes[m.a][2]) / m.L) ** 2));
    // heave added mass: displaced mass of each tube times the sine of its angle to the vertical
    const mAdd = J0.orig.reduce((s, m) => s + (m.tube ? inp.rhoW * Math.PI * m.D ** 2 / 4 * m.L * hor(m) : 0), 0);
    const lw = { mass: mt * inp.swCont, addedMass: mAdd, EA: inp.wireEA * 1e6, length: Math.max(depth * 0.3 + hHook + 30, 20), cableMass: 60 };
    const ar = axialResonance(lw), hv = Wsub > 0 && inp.tipHeave > 0 ? heaveResponse({ ...lw, heaveAmp: inp.tipHeave, period: inp.T, Cd: 1.0, area: J0.orig.filter((m) => m.group === 'horizontal' || m.group === 'plan').reduce((s, m) => s + m.D * m.L, 0), rho: inp.rhoW, staticTension: Wsub }) : { daf: 1, minTension: Wsub, slack: false, dynamicForce: 0 };
    // on-bottom stability before piling: legs flooded, mudmats on the seabed, installation wave
    M.members.length = s0; P.length = hook; M.supports = J0.mudNodes.map((n) => ({ node: n, fix: 'fixed' }));
    flood(true);
    const fi = { ...inp, block: inp.block ?? 1, mg: 0 }, dw = designWave(J0, { ...fi, Wtop: 0, FxTop: 0, MyTop: 0, windA: 0 }, inp.heading ?? 0, { Cd: Math.min(inp.Cd, 0.7) });
    const rb = solve(J0, [{ name: 'onbottom', nodal: [], member: [], selfWeight: { factor: inp.swCont, buoyancy: { waterline: inp.surge, rhoW: inp.rhoW } } }]).onbottom, Wb = rb.reactions.reduce((s, x) => s + x.R[2], 0);
    const Fh = Math.abs(dw.sc.FhMax), Mo = Math.abs(dw.sc.MhMax), Fz = Math.min(...dw.sc.Fz), cb = Math.cos((inp.heading ?? 0) * DEG), sb = Math.sin((inp.heading ?? 0) * DEG);
    const xs = J0.mudNodes.map((n) => (P[n][0] - J0.centre[0]) * cb + (P[n][1] - J0.centre[1]) * sb), bHalf = Math.max(...xs), Weff = Wb + Math.min(Fz, 0);
    const fsSlide = Fh > 0 ? inp.muSoil * Weff / Fh : 999, fsOT = Mo > 0 ? Weff * bHalf / Mo : 999, sx2 = xs.reduce((s, v) => s + v * v, 0), qMax = (Weff / J0.mudNodes.length + (sx2 > 0 ? Mo * bHalf / sx2 : 0)) / inp.mudA, qUlt = 6.05 * inp.suMud * 1e3, fsBear = qMax > 0 ? qUlt / qMax : 999;
    const steel = J0.steelT * inp.swCont, hookT = hookDyn / G / 1e3;
    const checks = [
      check('Crane capacity at the lift radius', hookT, inp.craneCap, 't', 'DNV-ST-N001 (2021) §16.2 (dynamic hook load ≤ crane curve)', 'Dynamic hook load (lift weight with contingency, rigging and DAF) against the crane capacity at the radius. Exceeding the curve overloads the crane or forces a shorter radius the vessel may not reach.'),
      check('Sling load against breaking load / safety factor', slingDes / 1e6, slingCap / 1e6, 'MN', `DNV-ST-N001 (2021) §16.4 (MBL/γsf with γsf = ${inp.slingSF}); skew load factor ${inp.SKL}`, 'Most loaded sling with dynamic amplification and skew load against its minimum breaking load divided by the nominal safety factor. A sling failure drops the jacket.'),
      check(`Padeye — ${pe.governing.name}`, pe.governing.demand / 1e6, pe.governing.capacity / 1e6, 'MPa', `${pe.governing.ref}; design load with consequence factor ${inp.gC}`, 'Governing stress in the padeye (pin bearing, tear-out, net section, cheek weld or base) under the factored sling load at the sling angle with 5 % out-of-plane. Lift points get the highest factors because they have no redundancy.'),
      ...groupChecks(groups, inp, 'lift'),
      fsCheck('Unpiled jacket: sliding', fsSlide, inp.fsStab, 'factor', 'ISO 19902:2020 §22.5; API RP 2A-WSD §12.4 (on-bottom stability before piling)', 'Friction under the mudmats against the wave and current force of the installation sea state. If it slides before the first piles are driven the jacket is off location.'),
      fsCheck('Unpiled jacket: overturning', fsOT, inp.fsStab, 'factor', 'ISO 19902:2020 §22.5; rigid-body moment balance about the leeward mudmats', 'Stabilising moment of the submerged weight against the wave overturning moment. Below the requirement the weather limit for the unpiled phase must be lowered or ballast added.'),
      fsCheck('Unpiled jacket: mudmat bearing', fsBear, inp.fsStab * 4 / 3, 'factor', `qult = 6.05·su (square footing, undrained) against the peak mudmat pressure; factor ${fmt(inp.fsStab * 4 / 3)}`, 'Peak pressure under the most loaded mudmat against the undrained bearing capacity of the seabed. If exceeded the jacket settles unevenly and cannot be levelled for piling.'),
      check('Hoist line stays taut during lowering', hv.dynamicForce / 1e6, Math.max(Wsub, 1) / 1e6, 'MN', 'DNV-RP-N103 (2019) §5.3 (dynamic force below static tension: no slack)', 'Dynamic line force from crane-tip heave against the submerged weight. If the dynamic force exceeds the static tension the line goes slack and snaps tight.'),
    ];
    const ok = checks.every((c) => c.uc <= 1), bad = checks.filter((c) => c.uc > 1), gW = checks.reduce((a, c) => (c.uc > a.uc ? c : a));
    return {
      kpis: [
        { key: 'liftWeight', label: 'Lift weight (with contingency)', unit: 't', value: W / G / 1e3, explain: 'Jacket weight in air with the weight contingency, as carried by the slings. The basis for vessel selection.' },
        { key: 'hookStatic', label: 'Static hook load', unit: 'MN', value: hookStat / 1e6, explain: 'Lift weight plus rigging weight hanging from the hook in air.' },
        { key: 'hookDyn', label: 'Dynamic hook load', unit: 'MN', value: hookDyn / 1e6, explain: 'Static hook load times the dynamic amplification factor: the load to compare with the crane curve.' },
        { key: 'hookT', label: 'Dynamic hook load', unit: 't', value: hookT, limit: inp.craneCap, sense: 'max', explain: 'The same in tonnes against the crane capacity at the lift radius.' },
        { key: 'slingMax', label: 'Design load of the most loaded sling', unit: 'MN', value: slingDes / 1e6, limit: slingCap / 1e6, sense: 'max', explain: 'Largest static sling force times DAF and skew load factor. It also sizes shackles and the padeye.' },
        { key: 'padeyeUc', label: 'Padeye utilisation', unit: '-', value: pe.governing.uc, limit: 1, sense: 'max', explain: 'Governing stress ratio of the padeye under the factored sling load.' },
        { key: 'ucLift', label: 'Highest member utilisation in the lift', unit: '-', value: gov.uc, limit: 1, sense: 'max', explain: 'Largest member utilisation while hanging: the top bay is in compression between the sling attachments, which no in-place case produces.' },
        { key: 'hookSub', label: 'Hook load fully submerged', unit: 'MN', value: Wsub / 1e6, explain: 'Static hook load when the jacket is below the surface with the chosen flooding. A small or negative value means the jacket floats and must be ballasted.' },
        { key: 'Tres', label: 'Hoist-line resonance period near the seabed', unit: 's', value: ar.periodLumped, explain: 'Axial natural period of the jacket on the hoist line. It should be well below the wave periods of the installation sea state to avoid resonant heave.' },
        { key: 'fsSliding', label: 'Unpiled: factor of safety against sliding', unit: '-', value: Math.min(fsSlide, 999), limit: inp.fsStab, sense: 'min', explain: 'Mudmat friction over the horizontal wave and current load in the installation sea state.' },
        { key: 'fsOverturn', label: 'Unpiled: factor of safety against overturning', unit: '-', value: Math.min(fsOT, 999), limit: inp.fsStab, sense: 'min', explain: 'Restoring over overturning moment about the leeward mudmats.' },
        { key: 'fsBearing', label: 'Unpiled: factor of safety on mudmat bearing', unit: '-', value: Math.min(fsBear, 999), limit: inp.fsStab * 4 / 3, sense: 'min', explain: 'Undrained bearing capacity over the peak mudmat pressure from weight and overturning.' },
        { key: 'maxDisp', label: 'Largest displacement in the lift', unit: 'mm', value: dp.max * 1e3, explain: 'Largest nodal displacement while hanging, relative to the hook, including sling stretch.' },
      ],
      checks,
      plots: [
        plot3d(J0, { title: 'Member utilisation while hanging from the hook', values: nodeValuesFromMembers(J0, rows.map((r) => r.uc)), label: 'Utilisation', unit: '-', u: res.member.displacements, explain: 'Utilisation of the jacket members in the lift with the deflected shape. The inclined slings squeeze the top frame together: look at the top horizontals and plan braces.' }),
        { type: 'line', title: 'Static hook load while lowering through the water line', x: { label: 'Water line elevation on the jacket (in-place datum)', unit: 'm' }, y: { label: 'Hook load', unit: 'MN' }, series: [{ name: inp.legsFlooded ? 'Legs flooding' : 'Legs sealed', x: wl.map((z) => round(z)), y: hookW.map((v) => round((v + rig) / 1e6)) }, { name: 'In air', x: [round(wl[0]), round(wl[wl.length - 1])], y: [round(hookStat / 1e6), round(hookStat / 1e6)], dash: true }], explain: 'Hook load as the jacket is lowered and more of it is submerged. The drop is the buoyancy; a curve reaching zero means the jacket floats before it reaches the seabed.' },
        { type: 'bar', title: 'Sling forces', categories: sStat.map((v, k) => `Sling ${k + 1}`), series: [{ name: 'Static (MN)', y: sStat.map((v) => round(v / 1e6)) }, { name: 'With DAF and skew (MN)', y: sStat.map((v) => round(v * inp.DAF * inp.SKL / 1e6)) }], y: { label: 'Sling force', unit: 'MN' }, explain: 'Force in each sling from the frame analysis, and with dynamic and skew load factors. Unequal static forces show an eccentric centre of gravity.' },
        { type: 'bar', title: 'Padeye utilisation by check', categories: pe.checks.map((c) => c.name), series: [{ name: 'Utilisation', y: pe.checks.map((c) => round(c.uc)) }], y: { label: 'Utilisation', unit: '-' }, explain: 'Stress ratio of each padeye check. Bearing and tear-out depend on plate thickness and edge distance; the base check on the padeye width and the sling angle.' },
      ],
      tables: [
        { title: 'Lift weight and hook load', columns: [{ key: 'c', label: 'Item' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }], rows: [{ c: 'Modelled steel', v: round(J0.steelT), u: 't' }, { c: 'Lift weight with contingency', v: round(W / G / 1e3), u: 't' }, { c: 'Rigging', v: inp.rigW, u: 't' }, { c: 'Static hook load', v: round(hookStat / G / 1e3), u: 't' }, { c: 'Dynamic hook load', v: round(hookT), u: 't' }, { c: 'Hook height above the leg tops', v: round(hHook), u: 'm' }, { c: 'Sling length', v: round(Math.hypot(rad, hHook)), u: 'm' }, { c: 'Hook load fully submerged', v: round(Wsub / G / 1e3), u: 't' }], explain: 'Build-up of the hook load from modelled steel to the dynamic value, with the rigging geometry that follows from the sling angle.' },
        { title: 'Most utilised members in the lift', columns: memberCols, rows: memberTableRows(rows, 12), explain: 'Members with the highest utilisation while hanging, with the action factor and DAF applied.' },
        { title: 'On-bottom stability before piling', columns: [{ key: 'c', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }], rows: [{ c: 'Submerged weight on the mudmats', v: round(Wb / 1e6), u: 'MN' }, { c: 'Wave and current base shear', v: round(Fh / 1e6), u: 'MN' }, { c: 'Overturning moment', v: round(Mo / 1e6), u: 'MN·m' }, { c: 'Largest wave uplift', v: round(Math.min(Fz, 0) / 1e6), u: 'MN' }, { c: 'Peak mudmat pressure', v: round(qMax / 1e3), u: 'kPa' }, { c: 'Bearing capacity 6.05·su', v: round(qUlt / 1e3), u: 'kPa' }, { c: 'Dynamic line force from tip heave', v: round(hv.dynamicForce / 1e6), u: 'MN' }], explain: 'Loads and resistances of the jacket standing on its mudmats in the installation sea state, and the line dynamics during lowering.' },
      ],
      notes: [`Lift weight ${fmt(W / G / 1e3)} t, dynamic hook load ${fmt(hookT)} t (DAF ${fmt(inp.DAF)}) against ${fmt(inp.craneCap)} t; most loaded sling ${fmt(slingDes / 1e6)} MN with skew factor ${fmt(inp.SKL)} against ${fmt(slingCap / 1e6)} MN.`,
        `Padeye governed by ${pe.governing.name.toLowerCase()} at ${fmt(pe.governing.uc)}; highest member utilisation ${fmt(gov.uc)} in ${gov.label.toLowerCase()} member ${gov.i}.`,
        `Lowering: hook load falls to ${fmt(Wsub / 1e6)} MN submerged (${inp.legsFlooded ? 'legs flooding' : 'legs sealed'}); line resonance ${fmt(ar.periodLumped)} s against a ${fmt(inp.T)} s wave${hv.slack ? ' — the line goes slack' : ''}.`,
        `Unpiled in H = ${fmt(inp.H)} m (${dw.wf.label}): sliding ${fmt(Math.min(fsSlide, 999))}, overturning ${fmt(Math.min(fsOT, 999))}, bearing ${fmt(Math.min(fsBear, 999))}.`],
      advice: [
        ok ? adv('engineering', 'ok', `The lift and set-down are feasible: hook load ${fmt(hookT)} t of ${fmt(inp.craneCap)} t, sling ${fmt(slingDes / slingCap)}, padeye ${fmt(pe.governing.uc)}, members ${fmt(gov.uc)}, unpiled factors ${fmt(Math.min(fsSlide, fsOT, 999))} or better.`, 'Confirm the weight and centre of gravity by weighing, check the padeye for the sling angle tolerance (±5°) and the leg locally under it, and set the weather limit of the unpiled phase from the on-bottom factors with the forecast uncertainty (alpha factor).', 'DNV-ST-N001 (2021) §16; ISO 19901-6:2009 §18–§19; ISO 19902:2020 §22', `Lift weight drives vessel class and fuel: ${vesselCo2(inp.liftDays)} Each 100 t of jacket steel avoided is also about ${fmt(100 * 1.9)} t CO₂e of material.`)
          : adv('engineering', 'fail', `Installation criteria are not met: ${bad.map((c) => `${c.name} (${fmt(c.uc)})`).join('; ')}.`, gW.name.startsWith('Unpiled') ? `Lower the installation wave limit, enlarge the mudmats (about ${fmt(inp.mudA * (gW.name.includes('bearing') ? gW.uc : 1))} m² per leg for bearing) or add ballast and skirts before release.` : gW.name.startsWith('Padeye') ? `Thicken the padeye main plate to about ${fmt(inp.peT * gW.uc / 0.9)} mm or add cheek plates.` : gW.name.startsWith('Sling') ? `Use slings of at least ${fmt(slingDes * inp.slingSF / G / 1e3)} t breaking load or a steeper sling angle.` : gW.name.startsWith('Crane') ? 'Use a larger crane vessel, reduce lift weight or split the lift.' : 'Add temporary lift bracing in the top bay or use a spreader frame so the slings pull vertically.', BASIS_MEMBER + '; DNV-ST-N001 (2021) §16', `Temporary lift steel is removed after installation; a spreader frame of about ${steelCo2(steel * 0.02)} is reusable on the next lift.`),
        adv('economic', 'info', `A ${fmt(hookT)} t dynamic hook load needs a heavy-lift vessel for about ${fmt(inp.liftDays)} day(s) on site.`, 'Vessel capacity steps are the cost drivers: staying under the next crane class is worth more than the steel it takes. Weather waiting in the unpiled phase is the main schedule risk, so the on-bottom factors translate directly into workable days.', `${BASIS_COST}; DNV-ST-N001 (2021) §2 (operation planning and weather windows)`),
      ],
      quantities: { steel_t: round(steel), vessel_days: { heavyLift: inp.liftDays } },
    };
  },
  verify: [
    // Single vertical tube lifted by one vertical sling: m = 7850·0.20570·95 = 153 398 kg → W = 1.5043 MN; with a weight contingency
    // of 1.0 and no rigging the static hook load is W and the dynamic hook load W·1.15 = 1.7300 MN.
    { name: 'Hook load of a single-sling lift: W = ρ·A·L·g, dynamic = W·DAF', source: 'Statics; DNV-ST-N001 §16.2 definition of dynamic hook load',
      inputs: { geom: testFrame('column'), swCont: 1, rigW: 0, DAF: 1.15, H: 0, Uc: 0, Uw: 0, nSub: 1 },
      expect: [{ kpi: 'hookStatic', value: 1.5043, tol: 0.001 }, { kpi: 'hookDyn', value: 1.7300, tol: 0.001 }] },
  ],
};
