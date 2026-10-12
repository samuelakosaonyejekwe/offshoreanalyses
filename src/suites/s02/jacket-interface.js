// Integrated topside–jacket load transfer: the deck frame on rigid or elastic jacket supports with imposed
// differential leg displacement, operating and storm load cases, leg reactions and their envelope, the resultant
// at the interface for the jacket analysis, and the strength of the deck-leg stab-in section.
import { tubularMemberCheck } from '../../engine/structural/tubular.js';
import { deckInputs, solverInput, buildDeck, loadCase, weights, solveDeck, framePlot, nodeDisp, interfaceResultant, equipmentOf, need, run, fmt, round, adv, check, steelCo2Text, CONFIRM, DECK_LIMITS, G } from './_common.js';

export default {
  id: 's02-jacket-interface', suite: 2, title: 'Topside-to-jacket load transfer',
  summary: 'Leg reactions, forces and moments passed from the topsides to the jacket for operating and storm cases, with jacket interface flexibility, differential leg displacement, uplift and the strength of the deck-leg connection section.',
  covers: [12], fidelity: 'numerical',
  standards: ['ISO 19901-3:2014 §9 (interface with the support structure)', 'ISO 19902:2020 §9.10 (action combinations), §13 (tubular members), §12 (modelling of interfaces)', 'API RP 2A-WSD (2014) §6', 'NORSOK N-004 (2022) §6.3'],
  models: ['Static equilibrium', 'Matrix stiffness method (3-D direct stiffness)', 'Elastic support (interface spring) model', 'Imposed support displacement', 'Elastic beam-column analysis'],
  limits: `${DECK_LIMITS} The jacket is represented by uncoupled translational and rotational springs at each leg (or rigid supports): coupling between legs through the jacket, pile–soil nonlinearity and wave loading on the jacket itself belong to the jacket analysis. Stab-in guides, shim plates and the weld detail of the connection are not checked.`,
  publishes: ['topsides.weight', 'topsides.cogZ', 'topsides.Fx', 'topsides.My', 'topsides.legMax', 'quantities.steel_t'],
  searchTerms: 'topsides jacket interface reactions load transfer deck leg stab-in support stiffness',
  inputs: [
    ...deckInputs(),
    { key: 'iface', label: 'Jacket interface model', type: 'select', default: 'spring', group: 'Jacket interface', role: 'bc', options: [{ value: 'rigid', label: 'Rigid (restraint as selected above)' }, { value: 'spring', label: 'Elastic springs at each leg' }], help: 'Springs replace the fixed or pinned restraint at every leg; at least the translational stiffnesses must be positive.' },
    { key: 'Kh', label: 'Horizontal stiffness per leg', unit: 'MN/m', type: 'number', default: 150, min: 0.1, max: 1e5, group: 'Jacket interface', role: 'bc', showIf: { key: 'iface', in: ['spring'] }, uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Jacket top stiffness seen by one leg, from a unit-load case in the jacket model (pile-head stiffness in series with the jacket frame).' },
    { key: 'Kv', label: 'Vertical stiffness per leg', unit: 'MN/m', type: 'number', default: 1500, min: 1, max: 1e6, group: 'Jacket interface', role: 'bc', showIf: { key: 'iface', in: ['spring'] }, uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'Kr', label: 'Rotational stiffness per leg', unit: 'MN·m/rad', type: 'number', default: 2000, min: 0, max: 1e7, group: 'Jacket interface', role: 'bc', showIf: { key: 'iface', in: ['spring'] }, help: '0 for a pinned stab-in; a welded leg on a stiff jacket leg is of the order of 4EI/L of the jacket leg bay.' },
    { key: 'dSettle', label: 'Imposed vertical displacement of one corner leg', unit: 'mm', type: 'number', default: 10, min: -200, max: 200, group: 'Jacket interface', role: 'bc', help: 'Differential settlement, fabrication tolerance or jacket deformation at the leg nearest the origin (downward positive). Applied in its own load case.' },
    { key: 'Hx', label: 'Storm lateral load along x (characteristic)', unit: 'kN', type: 'number', default: 650, min: 0, max: 5e4, group: 'Environmental and interface loads', role: 'bc', bind: 'topsides.windFx', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Wind resultant from the wind module, shared between levels in proportion to elevation.' },
    { key: 'Hy', label: 'Storm lateral load along y (characteristic)', unit: 'kN', type: 'number', default: 950, min: 0, max: 5e4, group: 'Environmental and interface loads', role: 'bc', bind: 'topsides.windFy' },
    { key: 'ax', label: 'Horizontal acceleration of the interface', unit: 'g', type: 'number', default: 0.03, min: 0, max: 1, group: 'Environmental and interface loads', role: 'bc', help: 'Wave-induced acceleration of the jacket top (or vessel motion for a floating support) acting on all topsides mass with the storm.' },
    { key: 'gGop', label: 'Operating: factor on permanent actions', unit: '-', type: 'number', default: 1.3, min: 0.8, max: 1.6, group: 'Acceptance criteria', role: 'criteria', help: `ISO 19902 gravity combination 1.3 G + 1.5 Q. ${CONFIRM}` },
    { key: 'gQop', label: 'Operating: factor on variable actions', unit: '-', type: 'number', default: 1.5, min: 0.8, max: 1.8, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'gGst', label: 'Storm: factor on permanent and variable actions', unit: '-', type: 'number', default: 1.1, min: 0.8, max: 1.6, group: 'Acceptance criteria', role: 'criteria', help: `ISO 19902 extreme combination 1.1 G + 1.1 Q + 1.35 E. ${CONFIRM}` },
    { key: 'gE', label: 'Storm: factor on environmental actions', unit: '-', type: 'number', default: 1.35, min: 0, max: 1.8, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'gMin', label: 'Storm with minimum weight: factor on permanent actions', unit: '-', type: 'number', default: 0.9, min: 0.5, max: 1.1, group: 'Acceptance criteria', role: 'criteria', help: `Favourable gravity for the uplift check (0.9 G, no variable load). ${CONFIRM}` },
    { key: 'redist', label: 'Allowed leg reaction change from leg mismatch', unit: '%', type: 'number', default: 15, min: 1, max: 50, group: 'Acceptance criteria', role: 'criteria', help: 'Project criterion on the redistribution caused by the imposed displacement, as a share of the governing leg reaction.' },
    solverInput(1),
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'legMax' },
  compute(inp, ctx) {
    // SI: MN/m → N/m, mm → m, kN → N, g → fraction of gravity in the load direction vector.
    const deck = buildDeck(inp, ctx), eqs = equipmentOf(deck, inp), zs = deck.levelsZ, zsum = zs.reduce((s, z) => s + z, 0);
    const spring = inp.iface === 'spring';
    if (spring) need(inp.Kh > 0 && inp.Kv > 0, 'Elastic interface: the horizontal and vertical stiffnesses must be positive, otherwise the topsides is unrestrained.');
    need(!(spring && inp.Kr === 0 && deck.supNodes.length < 3), 'With pinned springs the topsides needs at least three legs to be stable.');
    const lat = (hx, hy) => { const e = []; deck.levelNodes.forEach((ns, k) => { for (const n of ns) e.push({ node: n, F: [(hx * 1e3 * zs[k]) / zsum / ns.length, (hy * 1e3 * zs[k]) / zsum / ns.length, 0, 0, 0, 0] }); }); return e; };
    const a = inp.ax * inp.gE;
    const defs = [
      { name: 'Operating', lc: { gDead: inp.gGop, gEquip: inp.gGop, gLive: inp.gQop } },
      { name: 'Storm +x', lc: { gDead: inp.gGst, gEquip: inp.gGst, gLive: inp.gGst, dir: [a, 0, -1], extra: lat(inp.gE * inp.Hx, 0) } },
      { name: 'Storm −x', lc: { gDead: inp.gGst, gEquip: inp.gGst, gLive: inp.gGst, dir: [-a, 0, -1], extra: lat(-inp.gE * inp.Hx, 0) } },
      { name: 'Storm +y', lc: { gDead: inp.gGst, gEquip: inp.gGst, gLive: inp.gGst, dir: [0, a, -1], extra: lat(0, inp.gE * inp.Hy) } },
      { name: 'Storm −y', lc: { gDead: inp.gGst, gEquip: inp.gGst, gLive: inp.gGst, dir: [0, -a, -1], extra: lat(0, -inp.gE * inp.Hy) } },
      { name: 'Storm +x, minimum weight', lc: { gDead: inp.gMin, gEquip: inp.gMin, gLive: 0, dir: [a / inp.gMin, 0, -1], extra: lat(inp.gE * inp.Hx, 0) } },
      { name: 'Storm +y, minimum weight', lc: { gDead: inp.gMin, gEquip: inp.gMin, gLive: 0, dir: [0, a / inp.gMin, -1], extra: lat(0, inp.gE * inp.Hy) } },
      { name: 'Operating + leg displacement', lc: { gDead: inp.gGop, gEquip: inp.gGop, gLive: inp.gQop }, settle: true },
      { name: 'Weight only (unfactored)', lc: { gDead: 1, gEquip: 1, gLive: 0 } },
    ];
    // The settling leg: the support nearest the plan origin. With springs the displacement is imposed as the force
    // K·δ on the spring (equivalent to moving its far end); with rigid supports as a prescribed settlement.
    const sN = deck.supNodes.reduce((b, n) => { const p = deck.model.nodes[n], q = deck.model.nodes[b]; return Math.hypot(p[0] - deck.x0, p[1] - deck.y0) < Math.hypot(q[0] - deck.x0, q[1] - deck.y0) ? n : b; }, deck.supNodes[0]);
    const kh = inp.Kh * 1e6, kv = inp.Kv * 1e6, kr = inp.Kr * 1e6, dS = inp.dSettle / 1e3;
    const cases = defs.map((d, i) => {
      const l = loadCase(deck, inp, { name: `c${i}`, ...d.lc, eqs });
      if (d.settle && dS !== 0) { if (spring) l.nodal.push({ node: sN, F: [0, 0, -kv * dS, 0, 0, 0] }); else l.settlements = [{ node: sN, d: [0, 0, -dS, 0, 0, 0] }]; }
      return l;
    });
    const supports = spring ? deck.supNodes.map((n) => ({ node: n, spring: [kh, kh, kv, kr, kr, Math.max(kr, 1e3)] })) : deck.model.supports;
    const { ref, sol } = solveDeck(deck, cases, inp.nSub, { supports });
    // Reactions: at rigid supports from the solver; at springs R = −K u (plus the imposed-displacement force).
    const reactionsOf = (r, d) => deck.supNodes.map((n) => {
      const p = deck.model.nodes[n];
      if (!spring) { const s = r.reactions.find((x) => x.node === n); return { node: n, x: p[0], y: p[1], R: s ? Array.from(s.R) : [0, 0, 0, 0, 0, 0] }; }
      const u = r.displacements, k = [kh, kh, kv, kr, kr, Math.max(kr, 1e3)], R = k.map((kk, j) => -kk * u[6 * n + j]);
      if (d.settle && n === sN) R[2] -= kv * dS;
      return { node: n, x: p[0], y: p[1], R };
    });
    const res = defs.map((d, i) => {
      const r = sol.results[`c${i}`], legs = reactionsOf(r, d), ir = interfaceResultant(deck, { reactions: legs.map((l) => ({ node: l.node, R: l.R })) });
      return { name: d.name, r, legs, R: ir.R, maxV: Math.max(...legs.map((l) => l.R[2])), minV: Math.min(...legs.map((l) => l.R[2])), maxH: Math.max(...legs.map((l) => Math.hypot(l.R[0], l.R[1]))), maxM: Math.max(...legs.map((l) => Math.hypot(l.R[3], l.R[4]))) };
    });
    const w = weights(deck, inp, { eqs }), wRes = res[res.length - 1], env = res.slice(0, 8);
    need(Math.abs(wRes.R[2] - w.total) < 2e-3 * w.total + 1, 'Equilibrium check failed: the vertical reactions do not balance the topsides weight. Check the interface stiffnesses.');
    const legMax = Math.max(...env.map((x) => x.maxV)), legMin = Math.min(...env.map((x) => x.minV)), hMax = Math.max(...env.map((x) => x.maxH)), mMax = Math.max(...env.map((x) => x.maxM));
    const gC = env.reduce((p, c) => (c.maxV > p.maxV ? c : p)), storm = env.slice(1, 5), FxMax = Math.max(...storm.map((x) => Math.abs(x.R[0]))), FyMax = Math.max(...storm.map((x) => Math.abs(x.R[1]))), MyMax = Math.max(...storm.map((x) => Math.abs(x.R[4]))), MxMax = Math.max(...storm.map((x) => Math.abs(x.R[3])));
    // Deck-leg section at the interface under each leg's reaction set (ISO 19902 §13).
    const col = deck.sections.column; let stab = { uc: 0, name: '', c: '' };
    for (const c of env) for (const l of c.legs) {
      const k = run(() => tubularMemberCheck({ code: 'ISO', section: col, fy: deck.steel.fy, L: zs[0], forces: { N: -l.R[2], Vy: Math.hypot(l.R[0], l.R[1]), Vz: 0, T: l.R[5], My: Math.hypot(l.R[3], l.R[4]) + Math.hypot(l.R[0], l.R[1]) * 0, Mz: 0 }, options: { E: deck.steel.E } }), 'Deck-leg check');
      if (k.governing.uc > stab.uc) stab = { uc: k.governing.uc, name: k.governing.name, c: c.name, leg: l };
    }
    const settle = res[7], dLeg = settle.legs.map((l, i) => l.R[2] - res[0].legs[i].R[2]), dMax = Math.max(...dLeg.map(Math.abs));
    const mean = wRes.R[2] / deck.supNodes.length, uneven = Math.max(...wRes.legs.map((l) => l.R[2])) / mean;
    const checks = [
      check(`Deck-leg section at the interface (${stab.c})`, stab.uc, 1, 'interaction', `ISO 19902:2020 §13 ${stab.name}`, 'Axial force, shear and moment at the foot of the most loaded deck leg against the resistance of the leg tube. Above 1.0 the leg yields or buckles at the stab-in and the connection cannot transfer the topsides load.'),
      check('Leg reaction change from the imposed displacement', dMax / 1e3, inp.redist / 100 * legMax / 1e3, 'kN', `Project criterion: redistribution within ${inp.redist} % of the governing leg reaction (ISO 19901-3 §9, ISO 19902 §12 fabrication and installation tolerances)`, `A ${inp.dSettle} mm vertical mismatch at one leg redistributes load between legs through the stiffness of the deck. A stiff deck on a stiff jacket is sensitive: large changes mean shimming and weigh-off tolerances must be tight.`),
    ];
    if (legMin < 0) checks.push(check('Leg uplift', -legMin / 1e3, legMax / 1e3, 'kN', 'ISO 19901-3 §9; ISO 19902 §9.10 (0.9 G favourable)', 'A leg goes into tension in the storm with minimum weight. The stab-in weld or bolted connection must carry it; compare the uplift with that capacity (shown here against the largest compression for scale).'));
    const gov = checks.reduce((p, c) => (c.uc > p.uc ? c : p)), ok = gov.uc <= 1;
    const u = gC.r.displacements, nN = ref.model.nodes.length, legT = deck.orig.filter((o) => o.group === 'Deck leg').reduce((s, o) => s + o.sec.A * o.L * deck.steel.rho, 0) / 1e3;
    return {
      kpis: [
        { key: 'weight', label: 'Topsides operating weight', unit: 't', value: w.total / G / 1e3, explain: 'Unfactored structure, plate, secondary steel and equipment carried by the jacket.' },
        { key: 'cogZ', label: 'Centre of gravity above the interface', unit: 'm', value: w.cog[2], explain: 'Height of the operating centre of gravity: the lever arm of inertia loads in storm, earthquake and transport.' },
        { key: 'cogX', label: 'Centre of gravity offset in x from the plan centre', unit: 'm', value: w.cog[0] - (deck.x0 + deck.Lx / 2), explain: 'Eccentricity of the weight about the centre of the leg group; it makes the leg reactions unequal under gravity alone.' },
        { key: 'legMax', label: 'Largest leg compression', unit: 'kN', value: legMax / 1e3, explain: `Largest factored vertical reaction of one leg over all cases (${gC.name}): the design load for the jacket leg and pile below.` },
        { key: 'legMin', label: 'Smallest leg reaction', unit: 'kN', value: legMin / 1e3, explain: 'Lowest vertical reaction over all cases including storm with minimum weight; negative is uplift.' },
        { key: 'legShear', label: 'Largest leg shear', unit: 'kN', value: hMax / 1e3, explain: 'Largest horizontal reaction of one leg.' },
        { key: 'legMoment', label: 'Largest leg moment', unit: 'kN·m', value: mMax / 1e3, explain: 'Largest bending moment passed through one leg connection (zero for a pinned interface).' },
        { key: 'Fx', label: 'Storm shear at the interface along x', unit: 'kN', value: FxMax / 1e3, explain: 'Factored horizontal resultant passed to the jacket in the x storm cases (wind plus inertia).' },
        { key: 'My', label: 'Storm overturning moment about y', unit: 'kN·m', value: MyMax / 1e3, explain: 'Factored moment of the topsides loads about the plan centre of the interface in the x storm cases, including the weight eccentricity.' },
        { key: 'reactV', label: 'Sum of vertical reactions (weight case)', unit: 'kN', value: wRes.R[2] / 1e3, explain: 'Equilibrium check: equals the unfactored operating weight.' },
        { key: 'uneven', label: 'Largest leg share over the mean (weight only)', unit: '-', value: uneven, explain: 'How unevenly gravity is shared between legs: 1.0 is perfectly even. Values well above 1 come from eccentric equipment or stiffness differences.' },
        { key: 'dSettleLoad', label: 'Leg reaction change from the imposed displacement', unit: 'kN', value: dMax / 1e3, explain: 'Largest change of a leg reaction caused by the imposed vertical displacement of one leg.' },
      ],
      checks,
      plots: [
        framePlot({ title: `Displacement and deformed shape — ${gC.name}`, ref, u, values: nodeDisp(u, nN).map((v) => v * 1e3), valueLabel: 'Displacement', valueUnit: 'mm', explain: 'Frame displacement for the case giving the largest leg reaction. With an elastic interface part of the movement is the jacket springs compressing and rotating.' }),
        { type: 'bar', title: 'Vertical leg reactions by load case', categories: res[0].legs.map((l) => `(${fmt(l.x, 3)}, ${fmt(l.y, 3)})`), series: env.map((c) => ({ name: c.name, y: c.legs.map((l) => round(l.R[2] / 1e3)) })), y: { label: 'Vertical reaction', unit: 'kN' }, explain: 'Vertical reaction of each leg, labelled by plan position, for every load case. The spread between storm directions is the push–pull from overturning; the settlement case shows the redistribution through the deck.' },
        { type: 'bar', title: 'Interface resultants by load case', categories: env.map((c) => c.name), series: [{ name: 'Fx (kN)', y: env.map((c) => round(c.R[0] / 1e3)) }, { name: 'Fy (kN)', y: env.map((c) => round(c.R[1] / 1e3)) }, { name: 'Mx (10 kN·m)', y: env.map((c) => round(c.R[3] / 1e4)) }, { name: 'My (10 kN·m)', y: env.map((c) => round(c.R[4] / 1e4)) }], y: { label: 'Resultant', unit: 'kN, 10 kN·m' }, explain: 'Horizontal force and overturning moment of the support reactions about the plan centre of the interface for each case: these are the loads the jacket analysis applies at its top (with the opposite sign).' },
      ],
      tables: [
        { title: 'Interface load cases', columns: [{ key: 'name', label: 'Case' }, { key: 'Fx', label: 'Fx', unit: 'kN' }, { key: 'Fy', label: 'Fy', unit: 'kN' }, { key: 'Fz', label: 'Fz', unit: 'kN' }, { key: 'Mx', label: 'Mx', unit: 'kN·m' }, { key: 'My', label: 'My', unit: 'kN·m' }, { key: 'Mz', label: 'Mz', unit: 'kN·m' }, { key: 'maxV', label: 'Max leg', unit: 'kN' }, { key: 'minV', label: 'Min leg', unit: 'kN' }],
          rows: res.map((c) => ({ name: c.name, Fx: round(c.R[0] / 1e3), Fy: round(c.R[1] / 1e3), Fz: round(c.R[2] / 1e3), Mx: round(c.R[3] / 1e3), My: round(c.R[4] / 1e3), Mz: round(c.R[5] / 1e3), maxV: round(c.maxV / 1e3), minV: round(c.minV / 1e3) })),
          explain: 'Resultant of the support reactions about the plan centre of the interface (forces on the topsides; reverse the sign for loads on the jacket) and the extreme leg reactions of each case.' },
        { title: `Leg reactions — ${gC.name}`, columns: [{ key: 'x', label: 'x', unit: 'm' }, { key: 'y', label: 'y', unit: 'm' }, { key: 'Rx', label: 'Rx', unit: 'kN' }, { key: 'Ry', label: 'Ry', unit: 'kN' }, { key: 'Rz', label: 'Rz', unit: 'kN' }, { key: 'Mx', label: 'Mx', unit: 'kN·m' }, { key: 'My', label: 'My', unit: 'kN·m' }],
          rows: gC.legs.map((l) => ({ x: round(l.x), y: round(l.y), Rx: round(l.R[0] / 1e3), Ry: round(l.R[1] / 1e3), Rz: round(l.R[2] / 1e3), Mx: round(l.R[3] / 1e3), My: round(l.R[4] / 1e3) })), explain: 'Reactions of every leg in the governing case: the six-component set to apply to the corresponding jacket leg.' },
      ],
      notes: [`${gC.name} gives the largest leg reaction of ${fmt(legMax / 1e3)} kN; the weight case shares gravity ${fmt(uneven)} : 1 between the most loaded leg and the mean.`, `Interface: ${spring ? `springs Kh ${inp.Kh} MN/m, Kv ${inp.Kv} MN/m, Kr ${inp.Kr} MN·m/rad per leg` : `rigid (${inp.support})`}. Vertical equilibrium: ${fmt(wRes.R[2] / 1e3, 6)} kN against a weight of ${fmt(w.total / 1e3, 6)} kN.`,
        `Storm resultants: Fx ${fmt(FxMax / 1e3)} kN, Fy ${fmt(FyMax / 1e3)} kN, Mx ${fmt(MxMax / 1e3)} kN·m, My ${fmt(MyMax / 1e3)} kN·m (factored).`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `The interface transfers the loads: largest leg reaction ${fmt(legMax / 1e3)} kN, deck-leg utilisation ${fmt(stab.uc)} (${stab.c}), ${legMin < 0 ? `uplift ${fmt(-legMin / 1e3)} kN` : 'no uplift'}.`,
          `Issue the reaction table to the jacket analysis and iterate once with the jacket-top stiffness it returns; set the weigh-off and shimming tolerance from the ${fmt(dMax / 1e3)} kN change per ${inp.dSettle} mm of leg mismatch.`, 'ISO 19901-3:2014 §9; ISO 19902:2020 §9.10, §12–§13',
          `Every tonne on the topsides is carried again by jacket legs, piles and foundation: the ${fmt(w.total / G / 1e3)} t reported here drives steel below the interface as well (deck legs alone ${steelCo2Text(legT)}).`)
          : adv('engineering', 'fail', `${gov.name} is not met (utilisation ${fmt(gov.uc)}): largest leg reaction ${fmt(legMax / 1e3)} kN, deck-leg utilisation ${fmt(stab.uc)}, reaction change from leg mismatch ${fmt(dMax / 1e3)} kN.`,
            gov.name.startsWith('Deck-leg') ? `Increase the leg wall to about ${fmt(inp.ct * stab.uc / 0.9, 3)} mm at the stab-in (a thickened can of one diameter length is enough), or add legs to share the load.` : gov.name.startsWith('Leg uplift') ? `Design the stab-in weld or holding-down for ${fmt(-legMin / 1e3)} kN of tension, or lower the centre of gravity and wind area.` : 'Tighten the leg elevation tolerance, shim at installation after weigh-off, or soften the load path with a bearing detail; a stiffer deck on a stiff jacket cannot absorb the mismatch.',
            'ISO 19901-3:2014 §9; ISO 19902:2020 §13', `A thickened leg can at the four stab-ins adds about ${steelCo2Text(Math.max(stab.uc / 0.9 - 1, 0.05) * legT)} — small against the consequence of an overstressed interface.`),
        adv('economic', 'info', `The jacket and its piles are sized by the ${fmt(legMax / 1e3)} kN leg reaction and the storm resultant of ${fmt(Math.max(FxMax, FyMax) / 1e3)} kN reported here; the weight is ${fmt(w.total / G / 1e3)} t with its centre of gravity ${fmt(w.cog[2])} m above the interface.`,
          'Freeze the interface loads with an agreed contingency (typically 5–10 % on weight, more on centre-of-gravity position early on) before jacket fabrication: late topsides weight growth is the most common cause of jacket and pile re-work.', 'ISO 19901-5:2021 (weight control and reporting); AACE RP 18R-97'),
      ],
      publish: { 'topsides.weight': w.total, 'topsides.cogZ': w.cog[2], 'topsides.Fx': FxMax, 'topsides.My': MyMax, 'topsides.legMax': legMax, 'quantities.steel_t': deck.steelT },
      quantities: { steel_t: deck.steelT },
    };
  },
  verify: [{
    // One bay, one level at z = 4 m on four rigid legs, a 2000 t block at the plan centre, no plate or live load.
    // Weight from the hand take-off of the in-place module: 2059.63 t → ΣRz = 20 198 kN, shared equally: the largest
    // share over the mean is 1.000. Storm +x: lateral load 1000 kN (factor 1.0) applied at the deck, no inertia:
    //   Fx = 1000 kN; My about the interface = H·z = 1000·4 = 4000 kN·m.
    name: 'Statics of a symmetric deck: ΣRz = W, Fx = H, My = H·z', source: 'Equilibrium (hand calculation)',
    inputs: { iface: 'rigid', live: 0, Hx: 1000, Hy: 0, ax: 0, gE: 1, gGst: 1, tPlate: 0, secSteel: 0, nx: 1, ny: 1, nLevels: 1, bracing: 'none', columns: 'corners', dSettle: 0, equip: [{ name: 'Block', W: 2000, x: 15, y: 10, level: 1, a: 0, b: 0, zc: 0 }] },
    expect: [{ kpi: 'weight', value: 2059.63, tol: 0.001 }, { kpi: 'reactV', value: 20198.1, tol: 0.001 }, { kpi: 'Fx', value: 1000, tol: 0.001 }, { kpi: 'My', value: 4000, tol: 0.002 }, { kpi: 'uneven', value: 1, tol: 0.001 }],
  }],
};
