// Offshore lift of the topsides: lift weight and hook load with the code factors, sling forces from a frame
// analysis of the deck hanging in its four-sling rigging, sling and padeye checks, member utilisation in the lift
// condition, tilt, crane capacity and the dynamic hook load from crane-tip motion.
import { padeyeCheck } from '../../engine/structural/steelcheck.js';
import { axialResonance, heaveResponse, pendulumPeriod } from '../../engine/lines/lowering.js';
import { deckInputs, buildDeck, loadCase, weights, solveDeck, memberForces, memberChecks, groupSummary, framePlot, nodeValuesFromMembers, equipmentOf, need, run, fmt, round, linspace, adv, check, fsCheck, steelCo2Text, vesselCo2Text, CONFIRM, DECK_LIMITS, G, DEG } from './_common.js';

export default {
  id: 's02-lift-installation', suite: 2, title: 'Lifting and installation',
  summary: 'Lift weight, hook load, sling forces, sling and padeye utilisation, member utilisation of the deck hanging in its rigging, tilt, crane capacity and the dynamic amplification from crane-tip motion for a four-point offshore lift.',
  covers: [13], fidelity: 'numerical',
  standards: ['DNV-ST-N001 (2021) §16 (lifting operations: weight, DAF, skew load, sling and lift-point design)', 'ISO 19901-6:2009 §18 (lifting operations)', 'DNV-RP-N103 (2019) §4.3, §5.3 (crane-tip motion, dynamic hoist-line force)', 'API RP 2A-WSD (2014) §5.4 (lifting forces, padeye factors)', 'ISO 19902:2020 §13 / EN 1993-1-1 §6 (member checks in the lift condition)'],
  models: ['Quasi-static lifting analysis', 'Dynamic amplification approach', 'Matrix stiffness method with truss slings (multi-point support stiffness)', 'Static equilibrium', 'Single-degree-of-freedom hoist-line dynamics', 'Padeye section statics (bearing, tear-out, net section)'],
  limits: `${DECK_LIMITS} Four slings to a single hook above the centre of gravity, linear elastic and always in tension; spreader bars, multi-hook lifts, sling length tolerances (covered only by the skew load factor), shackles, tugger lines and set-down impact are not modelled. Dynamic hoist-line response is a single-degree-of-freedom model in air driven by a harmonic crane-tip motion — vessel motion analysis supplies that motion.`,
  publishes: ['install.liftWeight', 'install.hookLoad', 'install.vesselDays'],
  searchTerms: 'offshore heavy lift topsides sling forces padeye DAF skew load factor installation',
  inputs: [
    ...deckInputs({}, ['legs', 'support']),
    { key: 'rigging', label: 'Rigging weight', unit: 't', type: 'number', default: 25, min: 0, max: 2000, group: 'Lift weight and factors', role: 'load', help: 'Slings, shackles and spreaders below the hook.' },
    { key: 'cont', label: 'Weight contingency factor', unit: '-', type: 'number', default: 1.05, min: 1, max: 1.3, group: 'Lift weight and factors', role: 'criteria', help: `1.03 for a weighed structure, 1.05–1.10 for a calculated weight. ${CONFIRM}` },
    { key: 'daf', label: 'Dynamic amplification factor', unit: '-', type: 'number', default: 1.15, min: 1, max: 2, group: 'Lift weight and factors', role: 'bc', uncertain: { dist: 'normal', cov: 0.05 }, help: `DNV-ST-N001 Table 16-1 offshore: 1.30 below 100 t, 1.20 to 1000 t, 1.15 to 2500 t, 1.10 above. ${CONFIRM}` },
    { key: 'skl', label: 'Skew load factor', unit: '-', type: 'number', default: 1.25, min: 1, max: 1.75, group: 'Lift weight and factors', role: 'criteria', help: `Sling length tolerance in a statically indeterminate four-sling lift: 1.25 (DNV-ST-N001 §16.5). 1.0 with a load-equalising arrangement. ${CONFIRM}` },
    { key: 'cogF', label: 'Centre-of-gravity inaccuracy factor', unit: '-', type: 'number', default: 1.05, min: 1, max: 1.3, group: 'Lift weight and factors', role: 'criteria', help: CONFIRM },
    { key: 'gC', label: 'Consequence factor for lift points and their members', unit: '-', type: 'number', default: 1.3, min: 1, max: 1.5, group: 'Lift weight and factors', role: 'criteria', help: `DNV-ST-N001 Table 16-7: 1.30 lift points and attachments, 1.15 members framing into them, 1.0 other members (the higher value is applied to all members here). ${CONFIRM}` },
    { key: 'alpha', label: 'Smallest sling angle to the horizontal', unit: 'deg', type: 'number', default: 60, min: 30, max: 85, group: 'Rigging', role: 'geometry', design: true, help: 'The hook is placed above the centre of gravity at the height that gives the flattest sling this angle. Below 45° sling forces and deck compression rise steeply.' },
    { key: 'sling', label: 'Sling construction', type: 'select', source: 'materials:wire', default: 'wire-sixstrand', group: 'Rigging', role: 'material' },
    { key: 'dSling', label: 'Sling diameter', unit: 'mm', type: 'number', default: 170, min: 10, max: 500, group: 'Rigging', role: 'geometry', design: true },
    { key: 'sfReq', label: 'Required sling safety factor on the breaking load', unit: '-', type: 'number', default: 3, min: 1.5, max: 6, group: 'Rigging', role: 'criteria', help: `Nominal factor γsf on the dynamic sling force including bending and termination reductions (DNV-ST-N001 §16.7: not less than about 2.3–3.0). ${CONFIRM}` },
    { key: 'pt', label: 'Padeye main plate thickness', unit: 'mm', type: 'number', default: 120, min: 10, max: 300, group: 'Padeyes', role: 'geometry', design: true },
    { key: 'pR', label: 'Padeye main plate radius', unit: 'mm', type: 'number', default: 400, min: 40, max: 1500, group: 'Padeyes', role: 'geometry' },
    { key: 'pHole', label: 'Pin hole diameter', unit: 'mm', type: 'number', default: 186, min: 15, max: 800, group: 'Padeyes', role: 'geometry' },
    { key: 'pPin', label: 'Shackle pin diameter', unit: 'mm', type: 'number', default: 180, min: 12, max: 780, group: 'Padeyes', role: 'geometry' },
    { key: 'pCheek', label: 'Cheek plate thickness (each side)', unit: 'mm', type: 'number', default: 60, min: 0, max: 200, group: 'Padeyes', role: 'geometry' },
    { key: 'pRc', label: 'Cheek plate radius', unit: 'mm', type: 'number', default: 320, min: 0, max: 1400, group: 'Padeyes', role: 'geometry' },
    { key: 'pBase', label: 'Padeye base width', unit: 'mm', type: 'number', default: 1200, min: 80, max: 4000, group: 'Padeyes', role: 'geometry' },
    { key: 'pH', label: 'Hole centre above the base', unit: 'mm', type: 'number', default: 380, min: 40, max: 2000, group: 'Padeyes', role: 'geometry' },
    { key: 'oop', label: 'Out-of-plane sling force', unit: '%', type: 'number', default: 5, min: 0, max: 20, group: 'Padeyes', role: 'bc', help: 'Lateral component from sling misalignment: 3 % plus the actual angle is usual, 5 % minimum.' },
    { key: 'craneCap', label: 'Crane capacity at the lift radius', unit: 't', type: 'number', default: 2500, min: 10, max: 30000, group: 'Crane vessel', role: 'criteria' },
    { key: 'tipAmp', label: 'Crane-tip vertical motion amplitude', unit: 'm', type: 'number', default: 0.25, min: 0, max: 5, group: 'Crane vessel', role: 'bc', bind: 'vessel.heaveSig', help: 'Single amplitude of the crane tip in the limiting sea state, from the vessel motion analysis.' },
    { key: 'tipT', label: 'Crane-tip motion period', unit: 's', type: 'number', default: 8, min: 2, max: 30, group: 'Crane vessel', role: 'bc' },
    { key: 'hoistEA', label: 'Hoist system axial stiffness EA (all falls)', unit: 'MN', type: 'number', default: 12000, min: 10, max: 1e6, group: 'Crane vessel', role: 'bc', help: 'Sum over the reeved falls of the main hoist wire.' },
    { key: 'hoistL', label: 'Hoist line length, crane tip to hook', unit: 'm', type: 'number', default: 40, min: 2, max: 500, group: 'Crane vessel', role: 'ic' },
    { key: 'days', label: 'Heavy-lift vessel days for the operation', unit: 'days', type: 'number', default: 3, min: 0.2, max: 60, group: 'Crane vessel', role: 'operation' },
    { key: 'tiltLim', label: 'Tilt limit of the lifted deck', unit: 'deg', type: 'number', default: 2, min: 0.2, max: 10, group: 'Acceptance criteria', role: 'criteria', help: 'Stab-in of the legs normally needs the deck within about 2° of level.' },
  ],
  compute(inp, ctx) {
    // SI: t → kg, mm → m, MN → N. Legs hang from the deck in the lift, so the in-place supports are not used.
    const base = { ...inp, legs: inp.legs ?? 'corners', support: 'fixed' };
    const deck = buildDeck(base, ctx), eqs = equipmentOf(deck, base), w = weights(deck, base, { eqs });
    const wire = ctx.material(inp.sling); need(wire && wire.mblCoeff, `Sling construction '${inp.sling}' is not in the materials library.`);
    const d = inp.dSling / 1e3, MBL = wire.mblCoeff * d * d, EAs = wire.eaCoeff * d * d;
    // Lift points: the top-level nodes nearest the four plan corners.
    const top = deck.levelNodes[deck.levelNodes.length - 1], P = deck.model.nodes, zTop = deck.levelsZ[deck.levelsZ.length - 1];
    const corners = [[deck.x0, deck.y0], [deck.x0 + deck.Lx, deck.y0], [deck.x0 + deck.Lx, deck.y0 + deck.Ly], [deck.x0, deck.y0 + deck.Ly]];
    const lp = [...new Set(corners.map(([cx, cy]) => top.reduce((b, n) => (Math.hypot(P[n][0] - cx, P[n][1] - cy) < Math.hypot(P[b][0] - cx, P[b][1] - cy) ? n : b), top[0])))];
    need(lp.length === 4, 'Four distinct lift points are needed at the corners of the top level: the top deck has too few nodes.');
    const dist = lp.map((n) => Math.hypot(P[n][0] - w.cog[0], P[n][1] - w.cog[1]));
    need(Math.min(...dist) > 0.1, 'A lift point coincides with the centre of gravity in plan: the rigging cannot be arranged.');
    const Hh = Math.max(...dist) * Math.tan(inp.alpha * DEG), hook = [w.cog[0], w.cog[1], zTop + Hh];
    const Ls = lp.map((n) => Math.hypot(P[n][0] - hook[0], P[n][1] - hook[1], Hh)), ang = lp.map((n, i) => Math.asin(Hh / Ls[i]) / DEG);
    // Frame hanging in its slings: truss slings to a pinned hook, weak horizontal springs at two lift points hold the
    // pendulum and yaw rigid-body modes (they carry no load when the hook is above the centre of gravity).
    const hookNode = P.length, f = inp.cont * inp.daf;
    const extra = {
      nodes: [...P, hook], materials: { ...deck.model.materials, sling: { E: 1e11, nu: 0.3, rho: 0 } }, supports: [{ node: hookNode, fix: 'pinned' }, { node: lp[0], spring: [1e5, 1e5, 0, 0, 0, 0] }, { node: lp[2], spring: [1e5, 1e5, 0, 0, 0, 0] }],
      members: [...deck.model.members, ...lp.map((n) => ({ a: n, b: hookNode, type: 'truss', material: 'sling', section: { A: EAs / 1e11, Iy: 1e-8, Iz: 1e-8, J: 1e-8 }, group: 'Sling' }))],
    };
    const { ref, sol } = solveDeck(deck, [loadCase(deck, base, { gDead: f, gEquip: f, gLive: 0, eqs })], 1, extra);
    const r = sol.results.design, nM = deck.model.members.length;
    const T = lp.map((_, i) => { const mf = r.memberForces[ref.memberMap[nM + i][0]]; return Math.max(mf.endA.N, mf.endB.N, 0); });
    need(T.every((v) => Number.isFinite(v)), 'The lift analysis did not return sling forces: check the frame connectivity.');
    const vert = T.reduce((s, v, i) => s + v * Math.sin(ang[i] * DEG), 0), Wlift = w.total * inp.cont, hookStatic = Wlift + inp.rigging * 1e3 * G, hookDyn = hookStatic * inp.daf;
    need(Math.abs(vert - Wlift * inp.daf) < 5e-3 * Wlift * inp.daf, 'Equilibrium check failed: the vertical components of the sling forces do not balance the lift weight (a sling may have gone slack).');
    const Tdes = T.map((v) => v * inp.skl * inp.cogF), iS = Tdes.indexOf(Math.max(...Tdes)), Tmax = Tdes[iS];
    // Statically determinate estimate for comparison: weight shared by the inverse plan distance along each diagonal pair.
    const share = lp.map((n, i) => { const j = (i + 2) % 4; return dist[j] / (dist[i] + dist[j]) / 2; }), Tstat = share.map((s, i) => s * Wlift * inp.daf / Math.sin(ang[i] * DEG));
    // Members in the lift condition, with the consequence factor on the forces.
    const fm = memberForces(deck, ref, r).map((x) => ({ ...x, N: x.N * inp.gC, My: x.My * inp.gC, Mz: x.Mz * inp.gC, Vy: x.Vy * inp.gC, Vz: x.Vz * inp.gC, T: x.T * inp.gC }));
    const rows = memberChecks(deck, fm), gm = rows.reduce((a, c) => (c.uc > a.uc ? c : a)), groups = groupSummary(rows);
    const u = r.displacements, zl = lp.map((n) => u[6 * n + 2]);
    const tilt = Math.max(Math.abs(zl[0] - zl[2]) / Math.hypot(P[lp[0]][0] - P[lp[2]][0], P[lp[0]][1] - P[lp[2]][1]), Math.abs(zl[1] - zl[3]) / Math.hypot(P[lp[1]][0] - P[lp[3]][0], P[lp[1]][1] - P[lp[3]][1])) / DEG;
    const pe = run(() => padeyeCheck({ F: Tmax * inp.gC, fy: deck.steel.fy, t: inp.pt / 1e3, R: inp.pR / 1e3, dHole: inp.pHole / 1e3, dPin: inp.pPin / 1e3, tCheek: inp.pCheek / 1e3, rCheek: inp.pRc / 1e3, aCheekWeld: Math.max(0.5 * inp.pCheek, 6) / 1e3, baseWidth: inp.pBase / 1e3, holeHeight: inp.pH / 1e3, angleDeg: 90 - ang[iS], outOfPlaneFraction: inp.oop / 100 }), 'Padeye check');
    // Hoist-line dynamics in air: SDOF on the hoist stiffness driven by the crane-tip motion.
    const mLift = hookStatic / G, res = axialResonance({ mass: mLift, EA: inp.hoistEA * 1e6, length: inp.hoistL, cableMass: 0 });
    const dyn = (Tp) => heaveResponse({ mass: mLift, EA: inp.hoistEA * 1e6, length: inp.hoistL, cableMass: 0, heaveAmp: Math.max(inp.tipAmp, 1e-6), period: Tp, rho: 1.225, staticTension: hookStatic });
    const hr = run(() => dyn(inp.tipT), 'Hoist-line dynamics'), dafCalc = 1 + hr.dynamicForce / hookStatic;
    const pend = pendulumPeriod({ length: inp.hoistL + Hh + (zTop - w.cog[2]) });
    const checks = [
      fsCheck(`Sling ${iS + 1}: safety factor on the breaking load`, MBL / Tmax, inp.sfReq, 'DNV-ST-N001 §16.7 (nominal sling safety factor on the dynamic sling force with skew load and centre-of-gravity factors)', 'Minimum breaking load of the sling over its design force. Below the required factor the sling may fail by overload, with the wear, bending and termination losses a sling sees in service.'),
      check('Padeye', pe.governing.uc, 1, 'stress ratio', `${pe.governing.name} — ${pe.governing.ref || 'allowable stress on the governing section'}; design force includes the consequence factor ${inp.gC}`, 'Bearing, tear-out, net-section, cheek weld and base stresses of the most loaded lift point. A padeye failure drops the load: it is designed with the highest consequence factor and full NDT.'),
      check(`Member strength in the lift condition (${gm.group.toLowerCase()} ${gm.i + 1})`, gm.uc, 1, 'interaction', `${gm.ref}; forces × consequence factor ${inp.gC}`, 'The deck hangs from four points instead of standing on its legs: girders at the lift points go into compression from the inclined slings and legs hang in tension. Above 1.0 a member needs temporary strengthening for the lift.'),
      check('Hook load against crane capacity', hookStatic / G / 1e3, inp.craneCap, 't', 'DNV-ST-N001 §16.3 (static hook load ≤ rated capacity at the radius; crane charts include their own dynamic allowance)', 'Static hook load (lift weight with contingency plus rigging) against the crane curve at the lift radius. Exceeding it means a bigger vessel, a shorter radius or a split lift.'),
      check('Dynamic factor from crane-tip motion against the factor used', dafCalc, inp.daf, 'DAF', 'DNV-RP-N103 §5.3 (dynamic hoist-line force); DNV-ST-N001 §16.2.5', 'Dynamic factor computed from the hoist-line stiffness and the crane-tip motion in the limiting sea state. If it exceeds the code factor used, the sea-state limit must be lowered or the higher value used in design.'),
      check('Tilt of the lifted deck', tilt, inp.tiltLim, 'deg', 'DNV-ST-N001 §16.4 (hook above the centre of gravity; tilt within stab-in tolerance)', 'Out-of-level of the deck hanging in its slings from elastic sling stretch and frame deflection. Excess tilt prevents the legs stabbing into the jacket together.'),
    ];
    const gov = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = gov.uc <= 1;
    const nN = ref.model.nodes.length, vals = nodeValuesFromMembers({ ...ref, nodeMap: ref.nodeMap.slice(0, nM) }, rows.map((x) => x.uc)); while (vals.length < nN) vals.push(0);
    const Tps = linspace(3, 16, 40), dReq = inp.dSling * Math.sqrt(inp.sfReq / (MBL / Tmax));
    return {
      kpis: [
        { key: 'liftWeight', label: 'Lift weight (with contingency)', unit: 't', value: Wlift / G / 1e3, explain: 'Structure and equipment weight times the contingency factor, without rigging: the weight reported to the lift contractor.' },
        { key: 'hookStatic', label: 'Static hook load', unit: 'kN', value: hookStatic / 1e3, explain: 'Lift weight plus rigging weight hanging from the crane hook.' },
        { key: 'hookLoad', label: 'Dynamic hook load', unit: 't', value: hookDyn / G / 1e3, explain: 'Static hook load times the dynamic amplification factor: the design load for crane, hook and vessel.' },
        { key: 'slingFE', label: 'Largest sling force from the frame analysis', unit: 'kN', value: Math.max(...T) / 1e3, explain: 'Largest sling tension with weight contingency and dynamic factor, shared according to the stiffness of deck and slings (before skew load and centre-of-gravity factors).' },
        { key: 'slingDesign', label: 'Design sling force', unit: 'kN', value: Tmax / 1e3, explain: 'Largest sling force including the skew load and centre-of-gravity factors: the force for sling, shackle and padeye selection.' },
        { key: 'slingSF', label: 'Sling safety factor', unit: '-', value: MBL / Tmax, limit: inp.sfReq, sense: 'min', explain: `Breaking load of the ${inp.dSling} mm sling (${fmt(MBL / 1e3)} kN) over the design sling force.` },
        { key: 'hookHeight', label: 'Hook height above the lift points', unit: 'm', value: Hh, explain: 'Height of the hook over the top deck needed for the chosen sling angle: check it against the crane hook height at the lift radius.' },
        { key: 'ucPadeye', label: 'Padeye utilisation', unit: '-', value: round(pe.governing.uc, 5), limit: 1, sense: 'max', explain: `Governing padeye stress ratio (${pe.governing.name}).` },
        { key: 'ucMember', label: 'Highest member utilisation in the lift', unit: '-', value: round(gm.uc, 5), limit: 1, sense: 'max', explain: 'Governing member utilisation with the deck hanging in its slings and the consequence factor applied.' },
        { key: 'tilt', label: 'Tilt in the slings', unit: 'deg', value: tilt, limit: inp.tiltLim, sense: 'max', explain: 'Out-of-level between diagonally opposite lift points.' },
        { key: 'dafCalc', label: 'Dynamic factor from crane-tip motion', unit: '-', value: dafCalc, limit: inp.daf, sense: 'max', explain: `1 + dynamic hoist-line force / static hook load for a crane-tip amplitude of ${inp.tipAmp} m at ${inp.tipT} s.` },
        { key: 'Thoist', label: 'Natural period of the load on the hoist line', unit: 's', value: res.periodLumped, explain: 'Vertical bounce period of the lifted weight on the hoist wires. It should be well below the wave periods; resonance multiplies the hook load.' },
        { key: 'Tpend', label: 'Pendulum period of the hanging deck', unit: 's', value: pend, explain: 'Swing period of the load below the crane tip; vessel roll or pitch near this period makes the load difficult to control.' },
      ],
      checks,
      plots: [
        framePlot({ title: 'Deck in its rigging: utilisation and deformed shape', ref, u, values: vals, valueLabel: 'Utilisation (with consequence factor)', valueUnit: '-', explain: 'The deck hanging from the hook with its four slings. Members near the lift points carry the horizontal sling components in compression; the legs hang in tension.' }),
        { type: 'bar', title: 'Sling forces', categories: lp.map((n, i) => `Sling ${i + 1} (${fmt(P[n][0], 3)}, ${fmt(P[n][1], 3)})`), series: [{ name: 'Frame analysis', y: T.map((v) => round(v / 1e3)) }, { name: 'Statics (diagonal pairs)', y: Tstat.map((v) => round(v / 1e3)) }, { name: 'Design (skew and CoG factors)', y: Tdes.map((v) => round(v / 1e3)) }, { name: 'Breaking load / required factor', y: lp.map(() => round(MBL / inp.sfReq / 1e3)) }], y: { label: 'Sling force', unit: 'kN' },
          explain: 'Force in each sling from the frame analysis and from simple statics, the design value with the skew load and centre-of-gravity factors, and the allowable. Unequal bars follow the centre-of-gravity offset.' },
        { type: 'line', title: 'Dynamic factor against crane-tip motion period', x: { label: 'Motion period', unit: 's' }, y: { label: 'Dynamic factor', unit: '-' }, series: [{ name: 'Computed', x: Tps.map((v) => round(v)), y: Tps.map((v) => round(Math.min(1 + dyn(v).dynamicForce / hookStatic, 5))) }, { name: 'Factor used in design', x: [3, 16], y: [inp.daf, inp.daf], dash: true }],
          explain: 'Dynamic hook-load factor for the given crane-tip amplitude over a range of motion periods. Where the curve rises above the dashed design factor the lift must not be performed, which sets the limiting wave period.' },
        { type: 'bar', title: 'Member utilisation in the lift by group', categories: groups.map((g) => g.group), series: [{ name: 'Utilisation', y: groups.map((g) => round(g.uc)) }], y: { label: 'Utilisation', unit: '-' }, explain: 'Governing utilisation of each member group while hanging in the slings; compare with the in-place values to find members sized by the lift.' },
      ],
      tables: [
        { title: 'Lift weight report', columns: [{ key: 'item', label: 'Item' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }], rows: [{ item: 'Structure and equipment (net)', v: round(w.total / G / 1e3), u: 't' }, { item: 'With contingency', v: round(Wlift / G / 1e3), u: 't' }, { item: 'Rigging', v: inp.rigging, u: 't' }, { item: 'Static hook load', v: round(hookStatic / G / 1e3), u: 't' }, { item: 'Dynamic hook load', v: round(hookDyn / G / 1e3), u: 't' }, { item: 'Centre of gravity x, y, z', v: `${fmt(w.cog[0])}, ${fmt(w.cog[1])}, ${fmt(w.cog[2])}`, u: 'm' }, { item: 'Hook position x, y, z', v: `${fmt(hook[0])}, ${fmt(hook[1])}, ${fmt(hook[2])}`, u: 'm' }],
          explain: 'Build-up of the hook load from the net weight, and the hook position above the centre of gravity. This is the weight statement the crane vessel is selected on.' },
        { title: 'Slings', columns: [{ key: 'n', label: 'Sling' }, { key: 'L', label: 'Length', unit: 'm' }, { key: 'ang', label: 'Angle to horizontal', unit: 'deg' }, { key: 'T', label: 'Force (frame)', unit: 'kN' }, { key: 'Td', label: 'Design force', unit: 'kN' }, { key: 'sf', label: 'Safety factor' }],
          rows: lp.map((_, i) => ({ n: i + 1, L: round(Ls[i]), ang: round(ang[i]), T: round(T[i] / 1e3), Td: round(Tdes[i] / 1e3), sf: round(MBL / Tdes[i]) })), explain: 'Length, angle and force of each sling. Lengths are hook to lift point (add shackles and terminations when ordering).' },
      ],
      notes: [`${gov.name} governs at ${fmt(gov.uc)}. Lift weight ${fmt(Wlift / G / 1e3)} t, dynamic hook load ${fmt(hookDyn / G / 1e3)} t; ${gm.group.toLowerCase()} ${gm.i + 1} governs the members (${gm.ref}).`, `Sling vertical components sum to ${fmt(vert / 1e3, 6)} kN against ${fmt(Wlift * inp.daf / 1e3, 6)} kN (equilibrium).`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `The lift is feasible: hook load ${fmt(hookDyn / G / 1e3)} t dynamic, sling safety factor ${fmt(MBL / Tmax)} (required ${inp.sfReq}), padeye ${fmt(pe.governing.uc)}, members ${fmt(gm.uc)}, tilt ${fmt(tilt)}°.`,
          `Order slings as matched pairs within ±0.25 % length, locate the padeyes in line with the sling directions (${fmt(ang[iS])}° at the governing point), and limit the lift to sea states where the crane-tip motion keeps the dynamic factor below ${inp.daf} (computed ${fmt(dafCalc)}).`, 'DNV-ST-N001 (2021) §16; ISO 19901-6:2009 §18; DNV-RP-N103 (2019) §5.3',
          `One lift with ${fmt(inp.days)} vessel days: ${vesselCo2Text(inp.days)} A single integrated lift avoids the offshore hook-up emissions and hours of a modular installation.`)
          : adv('engineering', 'fail', `${gov.name} is not met (utilisation ${fmt(gov.uc)}): design sling force ${fmt(Tmax / 1e3)} kN, sling factor ${fmt(MBL / Tmax)}, padeye ${fmt(pe.governing.uc)}, members ${fmt(gm.uc)}, hook load ${fmt(hookDyn / G / 1e3)} t.`,
            gov.name.startsWith('Sling') ? `Use slings of at least ${fmt(Math.ceil(dReq / 5) * 5, 3)} mm, or doubled slings, or steepen the sling angle above ${inp.alpha}°.` : gov.name.startsWith('Padeye') ? `Thicken the padeye to about ${fmt(Math.ceil(inp.pt * pe.governing.uc / 0.9 / 5) * 5, 3)} mm total with cheeks, or enlarge its radius; governing mode: ${pe.governing.name}.` : gov.name.startsWith('Member') ? `Strengthen the ${gm.group.toLowerCase()} members framing into the lift points by about ${fmt((gm.uc - 1) * 100, 2)} %, steepen the slings, or add a spreader frame to remove the horizontal sling components from the deck.` : gov.name.startsWith('Hook') ? `Reduce the lift weight by ${fmt(hookStatic / G / 1e3 - inp.craneCap)} t (later installation of equipment), or select a vessel with more capacity at the radius.` : gov.name.startsWith('Dynamic') ? `Lower the limiting sea state so the crane-tip amplitude stays below ${fmt(inp.tipAmp * (inp.daf - 1) / Math.max(dafCalc - 1, 1e-6))} m, or design for a DAF of ${fmt(dafCalc)}.` : 'Shorten or lengthen sling pairs to level the deck, or move the hook over the measured centre of gravity after weighing.',
            'DNV-ST-N001 (2021) §16.2–16.9; ISO 19901-6:2009 §18', `Temporary lift strengthening of about ${steelCo2Text(Math.max(gm.uc - 1, 0.02) * deck.frameSteelT * 0.2)} is small beside a second lift campaign: ${vesselCo2Text(inp.days)}`),
        adv('economic', 'info', `The lift needs a crane of at least ${fmt(hookStatic / G / 1e3)} t at the radius for ${fmt(inp.days)} vessel day(s); the lift weight of ${fmt(Wlift / G / 1e3)} t decides the vessel class.`,
          'Hold a lift-weight budget with the contingency from the first weight report: crossing a crane-capacity step late means a larger vessel at a much higher day rate, or splitting the deck with offshore hook-up.', 'DNV-ST-N001 (2021) §5 (weight control for marine operations); ISO 19901-5:2021; AACE RP 18R-97'),
      ],
      publish: { 'install.liftWeight': Wlift, 'install.hookLoad': hookDyn, 'install.vesselDays': inp.days },
      quantities: { vessel_days: { heavyLift: inp.days } },
    };
  },
  verify: [{
    // Symmetric deck with a 2000 t block at the plan centre, weight 2059.63 t (hand take-off of the in-place module)
    // = 20 198.1 kN, lifted at the four corners of a 30 × 20 m deck with slings at 60° to the horizontal, all factors 1:
    //   T = W/(4 sin 60°) = 20 198.1/3.4641 = 5830.7 kN per sling; hook height = √(15² + 10²)·tan 60° = 31.22 m.
    name: 'Symmetric four-point lift: T = W / (4 sin α)', source: 'Statics of a four-sling lift (hand calculation); DNV-ST-N001 §16.5',
    inputs: { cont: 1, daf: 1, skl: 1, cogF: 1, rigging: 0, alpha: 60, live: 0, tPlate: 0, secSteel: 0, nx: 1, ny: 1, nLevels: 1, bracing: 'none', columns: 'corners', equip: [{ name: 'Block', W: 2000, x: 15, y: 10, level: 1, a: 0, b: 0, zc: 0 }] },
    expect: [{ kpi: 'slingFE', value: 5830.7, tol: 2e-3 }, { kpi: 'hookStatic', value: 20198.1, tol: 1e-3 }, { kpi: 'hookHeight', value: 31.22, tol: 1e-3 }],
  }],
};
