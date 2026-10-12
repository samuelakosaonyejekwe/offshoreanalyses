// In-place strength of the topsides frame: factored gravity (dead, equipment, live) plus a lateral environmental
// load, three-dimensional frame analysis and member code checks.
import { deckInputs, solverInput, buildDeck, loadCase, weights, solveDeck, memberForces, memberChecks, groupSummary, nodeValuesFromMembers, framePlot, nodeDisp, interfaceResultant,
  equipmentOf, need, fmt, round, adv, check, steelCo2Text, CONFIRM, DECK_LIMITS, G } from './_common.js';

export default {
  id: 's02-inplace-strength', suite: 2, title: 'In-place structural strength',
  summary: 'Member forces, moments, shears, torsion and utilisation of girders, beams, columns, legs and braces under factored dead, equipment, live and lateral loads, with the weight report.',
  covers: [1], fidelity: 'numerical',
  standards: ['ISO 19901-3:2014 §7–§9 (topsides structure, actions and strength)', 'ISO 19902:2020 §9.10 (partial action factors), §13 (tubular members)', 'EN 1993-1-1:2005 §6.2, §6.3 (I-section members)', 'NORSOK N-004 (2022) §6'],
  models: ['Static equilibrium', 'Matrix stiffness method (3-D direct stiffness)', 'Euler–Bernoulli and Timoshenko beam theories', 'Saint-Venant torsion', 'Generalised Hooke’s law', 'Elastic beam-column analysis', 'von Mises yielding'],
  limits: `${DECK_LIMITS} Loads are one factored combination; envelope the operating and extreme combinations by re-running. Class 4 sections are checked with gross properties (flagged by the engine), and plate-girder shear buckling is not checked.`,
  publishes: ['topsides.weight', 'topsides.cogZ', 'quantities.steel_t', 'topsides.maxUc'],
  searchTerms: 'offshore topsides in-place structural analysis deck frame member utilisation',
  inputs: [
    ...deckInputs(),
    { key: 'Hx', label: 'Lateral environmental load along x', unit: 'kN', type: 'number', default: 900, min: -5e4, max: 5e4, group: 'Environmental and interface loads', role: 'bc', bind: 'topsides.windFx', uncertain: { dist: 'lognormal', cov: 0.2 },
      help: 'Wind (and wave-in-deck where applicable) resultant, shared between deck levels in proportion to their elevation. The wind module publishes its value.' },
    { key: 'Hy', label: 'Lateral environmental load along y', unit: 'kN', type: 'number', default: 300, min: -5e4, max: 5e4, group: 'Environmental and interface loads', role: 'bc' },
    { key: 'gG', label: 'Partial factor on permanent actions', unit: '-', type: 'number', default: 1.1, min: 0.8, max: 1.6, group: 'Acceptance criteria', role: 'criteria', help: `ISO 19902 extreme-environment combination uses 1.1 G; gravity-only operating combination uses 1.3 G. ${CONFIRM}` },
    { key: 'gQ', label: 'Partial factor on variable actions', unit: '-', type: 'number', default: 1.1, min: 0.5, max: 1.8, group: 'Acceptance criteria', role: 'criteria', help: `1.1 Q with extreme environment; 1.5 Q for the gravity-only combination. ${CONFIRM}` },
    { key: 'gE', label: 'Partial factor on environmental actions', unit: '-', type: 'number', default: 1.35, min: 0, max: 1.8, group: 'Acceptance criteria', role: 'criteria', help: `ISO 19902 gives 1.35 on the extreme environmental action. ${CONFIRM}` },
    { key: 'gM0', label: 'γM0 (cross-section resistance)', unit: '-', type: 'number', default: 1.0, min: 1, max: 1.3, group: 'Acceptance criteria', role: 'criteria', help: `EN 1993-1-1 recommends 1.0; NORSOK N-004 uses 1.15. ${CONFIRM}` },
    { key: 'gM1', label: 'γM1 (member buckling resistance)', unit: '-', type: 'number', default: 1.0, min: 1, max: 1.3, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'Kbuck', label: 'Effective length factor of columns and braces', unit: '-', type: 'number', default: 1.0, min: 0.5, max: 2.5, group: 'Acceptance criteria', role: 'criteria', help: 'Take it from the buckling module (effective lengths of the sway mode) for unbraced frames.' },
    { key: 'latFrac', label: 'Laterally unrestrained fraction of a beam bay', unit: '-', type: 'number', default: 0.25, min: 0.02, max: 1, group: 'Acceptance criteria', role: 'criteria', help: 'Stringers and welded plate restrain the compression flange: 0.25 means restraints at quarter points; 1.0 means a bare beam.' },
    solverInput(2),
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'maxDisp' },
  compute(inp, ctx) {
    // SI: mm → m, kN → N, kPa → Pa, t → kg inside the shared builder.
    const deck = buildDeck(inp, ctx), eqs = equipmentOf(deck, inp);
    const zs = deck.levelsZ, zsum = zs.reduce((s, z) => s + z, 0);
    const extra = [];
    deck.levelNodes.forEach((ns, k) => { for (const n of ns) extra.push({ node: n, F: [(inp.Hx * 1e3 * inp.gE * zs[k]) / zsum / ns.length, (inp.Hy * 1e3 * inp.gE * zs[k]) / zsum / ns.length, 0, 0, 0, 0] }); });
    const lc = loadCase(deck, inp, { gDead: inp.gG, gEquip: inp.gG, gLive: inp.gQ, extra, eqs });
    const { ref, sol } = solveDeck(deck, [lc], inp.nSub);
    const r = sol.results.design, u = r.displacements, nN = ref.model.nodes.length;
    const forces = memberForces(deck, ref, r);
    const rows = memberChecks(deck, forces, { K: inp.Kbuck, gammaM0: inp.gM0, gammaM1: inp.gM1, latFrac: inp.latFrac });
    const groups = groupSummary(rows), gov = rows.reduce((a, c) => (c.uc > a.uc ? c : a));
    const disp = nodeDisp(u, nN), maxDisp = Math.max(...disp);
    const w = weights(deck, inp, { eqs }), ir = interfaceResultant(deck, r);
    const applied = inp.gG * (w.frame + w.plate + w.equip) + inp.gQ * (inp.live * 1e3 * deck.planArea * zs.length);
    const reactV = ir.R[2], baseShear = Math.hypot(ir.R[0], ir.R[1]), maxLeg = Math.max(...ir.rows.map((x) => x.Rz)), minLeg = Math.min(...ir.rows.map((x) => x.Rz));
    need(Math.abs(reactV - applied) < 1e-3 * applied + 1, 'Equilibrium check failed: the vertical reactions do not balance the applied loads.');
    const maxVM = Math.max(...rows.map((x) => x.vm));
    const checks = groups.map((g) => check(`Member strength — ${g.group}`, g.uc, 1, 'interaction', g.ref,
      `Governing axial, bending, shear and beam-column interaction among the ${g.n} ${g.group.toLowerCase()} member(s) (member ${g.member + 1}). Above 1.0 the member yields or buckles under the factored loads.`));
    if (minLeg < 0) checks.push(check('Leg uplift at the interface', -minLeg / 1e3, maxLeg / 1e3, 'kN', 'ISO 19901-3 §9 (interface design)', 'A deck leg in net tension must be held by its weld or stab-in connection; compare the uplift with the connection tension capacity.'));
    const maxUc = gov.uc, ok = maxUc <= 1;
    const govSteel = deck.orig.filter((o) => o.group === gov.group).reduce((s, o) => s + o.sec.A * o.L * deck.steel.rho, 0) / 1e3;
    const target = 0.9, dMass = govSteel * (maxUc / target - 1);
    const ucNode = nodeValuesFromMembers(ref, rows.map((x) => x.uc));
    const advice = [
      ok ? adv('engineering', maxUc > 0.9 ? 'warn' : 'ok', `Highest member utilisation ${fmt(maxUc)} in the ${gov.group.toLowerCase()} group (${gov.ref}); maximum displacement ${fmt(maxDisp * 1e3)} mm.`,
        maxUc < 0.6 ? `The ${gov.group.toLowerCase()} group has reserve: its section area could drop by about ${fmt((1 - maxUc / target) * 100, 2)} % before reaching ${target}. Run the gravity-only combination (1.3 G + 1.5 Q, no environment) and the lift, transport and blast cases before trimming — they often govern topsides steel.`
          : 'Keep the sizes; confirm the governing members under the gravity-only combination (1.3 G + 1.5 Q) and check joints, lateral restraint of compression flanges and the weight contingency.',
        'ISO 19901-3:2014 §7–§9; ISO 19902:2020 §9.10, §13; EN 1993-1-1:2005 §6.3',
        maxUc < 0.6 ? `Trimming the ${gov.group.toLowerCase()} group to ${target} utilisation would save about ${steelCo2Text(-dMass)}` : `Primary and secondary steel analysed: ${steelCo2Text(deck.steelT)}`)
        : adv('engineering', 'fail', `The ${gov.group.toLowerCase()} group is overstressed: utilisation ${fmt(maxUc)} (${gov.ref}), member ${gov.i + 1}, N = ${fmt(gov.N / 1e3)} kN, My = ${fmt(gov.My / 1e3)} kN·m.`,
          `Increase the section modulus or area of the ${gov.group.toLowerCase()} members by about ${fmt((maxUc / target - 1) * 100, 2)} % (to reach ${target}), add a column or leg under the loaded bay, or move heavy equipment toward the column lines.`,
          'ISO 19902:2020 §13; EN 1993-1-1:2005 §6.2–6.3',
          `Strengthening the ${gov.group.toLowerCase()} group adds about ${steelCo2Text(Math.max(dMass, 0.05 * govSteel))} A deeper section at the same mass is the lower-carbon route to bending capacity.`),
      adv('economic', 'info', `Topsides operating weight ${fmt(w.total / G / 1e3)} t, of which structural steel ${fmt(deck.steelT)} t (${fmt((deck.steelT / (w.total / G / 1e3)) * 100, 2)} %): frame ${fmt(deck.frameSteelT)} t, plate and secondary ${fmt(deck.plateSteelT)} t.`,
        'Track this weight against the lift-vessel capacity and the jacket design load with a weight contingency (typically 5–10 % at detailed design). Every tonne removed from the topsides also removes jacket and pile steel below it.',
        'ISO 19901-5:2021 (weight management); AACE RP 18R-97 estimate classes'),
    ];
    return {
      kpis: [
        { key: 'maxUc', label: 'Highest member utilisation', unit: '-', value: round(maxUc, 5), limit: 1, sense: 'max', explain: 'Largest strength or beam-column utilisation of any member under the factored loads. Above 1.0 the member is overstressed.' },
        { key: 'maxDisp', label: 'Maximum nodal displacement', unit: 'mm', value: maxDisp * 1e3, explain: 'Largest translation of any node under the factored loads (use the serviceability module for unfactored deflection limits).' },
        { key: 'maxVM', label: 'Highest equivalent stress (envelope)', unit: 'MPa', value: maxVM / 1e6, explain: 'Largest von Mises stress in any member, combining the extreme axial, bending, shear and torsion values along the member as if they coincided (conservative).' },
        { key: 'reactV', label: 'Total vertical reaction', unit: 'kN', value: reactV / 1e3, explain: 'Sum of vertical reactions at the jacket interface; it equals the factored applied vertical load (equilibrium check).' },
        { key: 'baseShear', label: 'Total horizontal reaction', unit: 'kN', value: baseShear / 1e3, explain: 'Horizontal force passed to the jacket at the interface under the factored lateral loads.' },
        { key: 'maxLeg', label: 'Largest leg reaction', unit: 'kN', value: maxLeg / 1e3, explain: 'Largest vertical reaction of a single deck leg: the design load for the leg-to-jacket connection and the jacket leg below.' },
        { key: 'weight', label: 'Operating weight (unfactored)', unit: 't', value: w.total / G / 1e3, explain: 'Structure, plate, secondary steel and equipment, without live load: the weight carried by the jacket in operation.' },
        { key: 'cogZ', label: 'Centre of gravity above the interface', unit: 'm', value: w.cog[2], explain: 'Height of the operating centre of gravity above the jacket interface: it drives overturning in storms, earthquakes and sea transport.' },
        { key: 'steelT', label: 'Structural steel', unit: 't', value: deck.steelT, explain: 'Frame members plus deck plate and the secondary steel allowance.' },
      ],
      checks,
      plots: [
        framePlot({ title: 'Member utilisation and deformed shape', ref, u, values: ucNode, valueLabel: 'Utilisation', valueUnit: '-', explain: 'Frame coloured by the highest utilisation of the members meeting at each node. Look for values near or above 1.0 and for the load path from heavy equipment to the legs.' }),
        { type: 'bar', title: 'Governing utilisation by member group', categories: groups.map((g) => g.group), series: [{ name: 'Utilisation', y: groups.map((g) => round(g.uc)) }], y: { label: 'Utilisation', unit: '-' }, explain: 'Highest utilisation in each member group; anything above 1.0 fails, and a very low bar marks a group that could be lightened.' },
        { type: 'bar', title: 'Leg reactions at the jacket interface', categories: ir.rows.map((x) => `(${fmt(x.x, 3)}, ${fmt(x.y, 3)})`), series: [{ name: 'Vertical', y: ir.rows.map((x) => round(x.Rz / 1e3)) }, { name: 'Horizontal', y: ir.rows.map((x) => round(Math.hypot(x.Rx, x.Ry) / 1e3)) }], y: { label: 'Reaction', unit: 'kN' }, explain: 'Vertical and horizontal reaction of each leg, labelled by plan position. Uneven bars show eccentric equipment or overturning from the lateral load.' },
      ],
      tables: [
        { title: 'Governing member forces by group', columns: [{ key: 'group', label: 'Group' }, { key: 'n', label: 'Members' }, { key: 'member', label: 'Governing member' }, { key: 'N', label: 'N', unit: 'kN' }, { key: 'My', label: 'My', unit: 'kN·m' }, { key: 'Mz', label: 'Mz', unit: 'kN·m' }, { key: 'V', label: 'V', unit: 'kN' }, { key: 'T', label: 'T', unit: 'kN·m' }, { key: 'vm', label: 'von Mises', unit: 'MPa' }, { key: 'uc', label: 'UC' }, { key: 'ref', label: 'Governing check' }],
          rows: groups.map((g) => ({ group: g.group, n: g.n, member: g.member + 1, N: round(g.N / 1e3), My: round(g.My / 1e3), Mz: round(g.Mz / 1e3), V: round(g.V / 1e3), T: round(g.T / 1e3), vm: round(g.vm / 1e6), uc: round(g.uc), ref: g.ref })),
          explain: 'One row per member group with the forces of its most utilised member (N positive in tension; extreme values along the member). Use it to size each group.' },
        { title: 'Weight report', columns: [{ key: 'item', label: 'Item' }, { key: 'W', label: 'Weight', unit: 't' }], rows: [{ item: 'Frame steel', W: round(w.frame / G / 1e3) }, { item: 'Deck plate and secondary steel', W: round(w.plate / G / 1e3) }, { item: 'Equipment (operating)', W: round(w.equip / G / 1e3) }, { item: 'Total operating weight', W: round(w.total / G / 1e3) }, { item: 'Live load (unfactored, all levels)', W: round((inp.live * 1e3 * deck.planArea * zs.length) / G / 1e3) }],
          explain: 'Unfactored weight breakdown. The total and its centre of gravity are published to the jacket, lift and transport analyses.' },
      ],
      notes: [`${gov.group} member ${gov.i + 1} governs at ${fmt(maxUc)} (${gov.ref}).`, `Combination: ${inp.gG} G + ${inp.gQ} Q + ${inp.gE} E. Centre of gravity at (${fmt(w.cog[0])}, ${fmt(w.cog[1])}, ${fmt(w.cog[2])}) m.`,
        `Vertical equilibrium: reactions ${fmt(reactV / 1e3, 6)} kN against applied ${fmt(applied / 1e3, 6)} kN.`],
      advice,
      publish: { 'topsides.weight': w.total, 'topsides.cogZ': w.cog[2], 'quantities.steel_t': deck.steelT, 'topsides.maxUc': maxUc },
      quantities: { steel_t: deck.steelT },
    };
  },
  verify: [{
    // One bay, one level, a 2000 t block at the plan centre, no plate, no live or lateral load: the total reaction
    // must equal γG × g × (equipment + frame steel), and the weight report must reproduce the hand take-off below.
    name: 'Vertical equilibrium: total reaction equals γG × operating weight (no live, no lateral load)',
    source: 'Statics: ΣFz = 0 (hand check, reactV = 1.3 · weight · g)',
    inputs: { live: 0, Hx: 0, Hy: 0, gG: 1.3, tPlate: 0, secSteel: 0, nx: 1, ny: 1, nLevels: 1, bracing: 'none', columns: 'corners', equip: [{ name: 'Block', W: 2000, x: 15, y: 10, level: 1, a: 0, b: 0, zc: 0 }] },
    // Frame steel: girders 2 × 30 m × A_g, beams 2 × 20 m × A_b, legs 4 × 4 m × A_c.
    // A_g = 2·0.5·0.04 + 1.42·0.02 = 0.0684 m²; A_b = 2·0.45·0.035 + 1.13·0.018 = 0.05184 m²; A_c = π(0.914−0.032)·0.032 = 0.088668 m².
    // Mass = 7850·(60·0.0684 + 40·0.05184 + 16·0.088668) = 7850·7.5963 = 59 631 kg → 59.63 t.
    // weight = 2059.63 t; reactV = 1.3 · 2059.63 · 9.80665 = 26 257.6 kN.
        // Girders and beams (48.49 t) and the block sit at z = 4 m, the legs (11.14 t) are centred at z = 2 m:
    // cogZ = (2048.49·4 + 11.14·2)/2059.63 = 3.989 m.
    expect: [{ kpi: 'weight', value: 2059.63, tol: 0.001 }, { kpi: 'reactV', value: 26257.6, tol: 0.001 }, { kpi: 'cogZ', value: 3.989, tol: 0.002 }],
  }],
};
