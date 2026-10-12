// Serviceability of the topsides frame: vertical deflections, member sag, lateral drift and joint rotations under
// the unfactored (characteristic) combination.
import { memberDiagram } from '../../engine/structural/frame3d.js';
import { deckInputs, solverInput, buildDeck, loadCase, solveDeck, framePlot, equipmentOf, fmt, round, adv, check, steelCo2Text, CONFIRM, DECK_LIMITS } from './_common.js';

export default {
  id: 's02-deflection-serviceability', suite: 2, title: 'Deflection and serviceability',
  summary: 'Deck vertical deflections, member sag, lateral displacement, storey drift and joint rotations under characteristic loads against span-ratio limits.',
  covers: [4], fidelity: 'numerical',
  standards: ['ISO 19901-3:2014 §7.5 (serviceability)', 'ISO 19902:2020 §9 (serviceability limit states)', 'EN 1990:2002 A1.4 (deflection criteria)', 'NORSOK N-001 (2020) §7.2.5'],
  models: ['Static equilibrium', 'Matrix stiffness method (3-D direct stiffness)', 'Euler–Bernoulli and Timoshenko beam theories', 'Shear deformation relationships', 'Moment–curvature integration along members'],
  limits: `${DECK_LIMITS} Deflections are elastic and first-order; creep, bolt slip, fabrication camber and jacket deformation are not included. Limits are project criteria.`,
  searchTerms: 'offshore deck deflection serviceability limit span ratio lateral drift topsides',
  inputs: [
    ...deckInputs(),
    { key: 'swFac', label: 'Self-weight factor', unit: '-', type: 'number', default: 1, min: 0, max: 1.5, group: 'Permanent and equipment loads', role: 'load', help: '1 includes the structural self-weight; 0 isolates the deflection from equipment and live load (useful when the deck is cambered for dead load).' },
    { key: 'Hx', label: 'Lateral load along x (characteristic)', unit: 'kN', type: 'number', default: 650, min: -5e4, max: 5e4, group: 'Environmental and interface loads', role: 'bc', bind: 'topsides.windFx', help: 'Wind resultant for drift, shared between levels in proportion to elevation.' },
    { key: 'Hy', label: 'Lateral load along y (characteristic)', unit: 'kN', type: 'number', default: 0, min: -5e4, max: 5e4, group: 'Environmental and interface loads', role: 'bc' },
    { key: 'defRatio', label: 'Vertical deflection limit (span / value)', unit: '-', type: 'number', default: 250, min: 100, max: 1000, group: 'Acceptance criteria', role: 'criteria', help: `L/250 total is common; L/360 where brittle finishes or sensitive equipment are supported. ${CONFIRM}` },
    { key: 'driftRatio', label: 'Lateral drift limit (height / value)', unit: '-', type: 'number', default: 300, min: 100, max: 1000, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'rotLimit', label: 'Joint rotation limit', unit: 'mrad', type: 'number', default: 5, min: 0.5, max: 50, group: 'Acceptance criteria', role: 'criteria', help: 'Slope limit at equipment supports (alignment of rotating machinery, drainage, piping flanges). Vendor data governs.' },
    solverInput(2),
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'dAbs' },
  compute(inp, ctx) {
    const deck = buildDeck(inp, ctx), eqs = equipmentOf(deck, inp);
    const zs = deck.levelsZ, zsum = zs.reduce((s, z) => s + z, 0), extra = [];
    deck.levelNodes.forEach((ns, k) => { for (const n of ns) extra.push({ node: n, F: [(inp.Hx * 1e3 * zs[k]) / zsum / ns.length, (inp.Hy * 1e3 * zs[k]) / zsum / ns.length, 0, 0, 0, 0] }); });
    const { ref, sol } = solveDeck(deck, [loadCase(deck, inp, { gDead: inp.swFac, extra, eqs })], inp.nSub);
    const r = sol.results.design, u = r.displacements, P = ref.model.nodes, nN = P.length;
    // Vertical deflection relative to the mean settlement of the leg heads (the level-1 nodes above the supports).
    const heads = deck.supNodes.map((s) => { const p = deck.model.nodes[s]; return deck.levelNodes[0].find((n) => Math.hypot(deck.model.nodes[n][0] - p[0], deck.model.nodes[n][1] - p[1]) < 1e-6); }).filter((n) => n !== undefined);
    const headZ = heads.length ? heads.reduce((s, n) => s + u[6 * n + 2], 0) / heads.length : 0;
    let dAbs = 0, dRel = 0, rot = 0, lat = 0, latNode = 0;
    const vert = new Array(nN);
    for (let i = 0; i < nN; i++) {
      vert[i] = -u[6 * i + 2]; dAbs = Math.max(dAbs, Math.abs(u[6 * i + 2])); if (P[i][2] >= zs[0] - 1e-6) dRel = Math.max(dRel, Math.abs(u[6 * i + 2] - headZ));   // leg nodes below the cellar deck are not deck deflection
      rot = Math.max(rot, Math.hypot(u[6 * i + 3], u[6 * i + 4]));
      const h = Math.hypot(u[6 * i], u[6 * i + 1]); if (h > lat) { lat = h; latNode = i; }
    }
    // Storey drift between consecutive levels (mean horizontal displacement of each level).
    const levU = deck.levelNodes.map((ns) => [ns.reduce((s, n) => s + u[6 * n], 0) / ns.length, ns.reduce((s, n) => s + u[6 * n + 1], 0) / ns.length]);
    const storeys = zs.map((z, k) => { const below = k === 0 ? [0, 0] : levU[k - 1], h = k === 0 ? z : z - zs[k - 1], d = Math.hypot(levU[k][0] - below[0], levU[k][1] - below[1]); return { level: k + 1, z, h, ux: levU[k][0], uy: levU[k][1], drift: d, ratio: d > 0 ? h / d : Infinity }; });
    const worstDrift = storeys.reduce((a, c) => (c.drift / c.h > a.drift / a.h ? c : a));
    // Member sag relative to the chord through its ends.
    let sag = { v: 0, L: 1, i: 0, group: '' };
    deck.orig.forEach((o, i) => {
      if (!o.horizontal) return;
      for (const k of ref.memberMap[i]) {
        const d = memberDiagram(r, k, 9), n = d.s.length, L = d.s[n - 1];
        for (let j = 0; j < n; j++) { const v = Math.abs(d.uz[j] - d.uz[0] - (d.uz[n - 1] - d.uz[0]) * d.s[j] / L); if (v / o.L > sag.v / sag.L) sag = { v, L: o.L, i, group: o.group }; }
      }
    });
    // With sub-elements the chord is that of the sub-element: add the node offsets from the member chord.
    deck.orig.forEach((o, i) => {
      if (!o.horizontal) return;
      const ns = ref.nodeMap[i], a = ns[0], b = ns[ns.length - 1];
      ns.forEach((n, j) => { const v = Math.abs(u[6 * n + 2] - u[6 * a + 2] - (u[6 * b + 2] - u[6 * a + 2]) * j / (ns.length - 1)); if (v / o.L > sag.v / sag.L) sag = { v, L: o.L, i, group: o.group }; });
    });
    const spanG = Math.min(deck.Lx, deck.Ly) / (inp.legs === 'all' && deck.grid ? Math.max(inp.nx, inp.ny) : 1);
    const spanGlobal = inp.legs === 'all' && deck.grid ? deck.span : spanG;
    const limG = spanGlobal / inp.defRatio, limS = sag.L / inp.defRatio, H = zs[zs.length - 1], limL = H / inp.driftRatio;
    const checks = [
      check('Deck deflection relative to the leg heads', dRel * 1e3, limG * 1e3, 'mm', `Span / ${inp.defRatio} on the ${fmt(spanGlobal)} m span between leg lines (ISO 19901-3 §7.5; EN 1990 A1.4)`, 'Limits the overall sag of the deck between its supports. Excess sag misaligns piping and rotating equipment, ponds water and is visible.'),
      check(`Member sag — ${sag.group || 'beam'} ${sag.i + 1}`, sag.v * 1e3, limS * 1e3, 'mm', `Member span / ${inp.defRatio}`, 'Deflection of a single girder or beam relative to the chord through its ends: the local curvature that equipment feet and plating must follow.'),
      check('Lateral displacement at the top level', lat * 1e3, limL * 1e3, 'mm', `Height / ${inp.driftRatio}`, 'Sidesway of the topsides under wind. Large drift strains bridge landings, piping and conductors crossing the interface and feels uncomfortable.'),
      check(`Storey drift — level ${worstDrift.level}`, worstDrift.drift * 1e3, (worstDrift.h / inp.driftRatio) * 1e3, 'mm', `Storey height / ${inp.driftRatio}`, 'Relative sway between consecutive deck levels, which racks cladding, firewalls and vertical piping.'),
      check('Joint rotation', rot * 1e3, inp.rotLimit, 'mrad', 'Project / vendor alignment criterion', 'Largest slope of a joint. Excess rotation misaligns machinery trains and flanged connections seated on the deck.'),
    ];
    const gov = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = gov.uc <= 1;
    // Deflection ∝ 1/I ∝ about 1/h² at constant area: depth increase needed to bring the governing vertical check to 1.0.
    const depthFac = Math.sqrt(Math.max(checks[0].uc, checks[1].uc));
    const beamT = deck.orig.filter((o) => o.horizontal).reduce((s, o) => s + o.sec.A * o.L * deck.steel.rho, 0) / 1e3;
    return {
      kpis: [
        { key: 'dAbs', label: 'Largest vertical displacement', unit: 'mm', value: dAbs * 1e3, explain: 'Largest absolute vertical movement of any node, including the elastic shortening of legs and columns.' },
        { key: 'dRel', label: 'Deck deflection relative to the leg heads', unit: 'mm', value: dRel * 1e3, limit: limG * 1e3, sense: 'max', explain: 'Vertical deflection measured from the mean level of the leg heads: the sag of the deck between its supports.' },
        { key: 'spanRatio', label: 'Span over deflection', unit: '-', value: dRel > 0 ? Math.min(spanGlobal / dRel, 1e5) : 1e5, limit: inp.defRatio, sense: 'min', explain: 'Span between leg lines divided by the relative deflection; compare with L/250 to L/360.' },
        { key: 'sag', label: 'Largest member sag', unit: 'mm', value: sag.v * 1e3, limit: limS * 1e3, sense: 'max', explain: 'Largest deflection of one girder or beam relative to the chord through its ends.' },
        { key: 'lateral', label: 'Largest lateral displacement', unit: 'mm', value: lat * 1e3, limit: limL * 1e3, sense: 'max', explain: 'Largest horizontal movement of any node under the characteristic lateral load.' },
        { key: 'drift', label: 'Worst storey drift ratio (h / Δ)', unit: '-', value: Math.min(worstDrift.ratio, 1e5), limit: inp.driftRatio, sense: 'min', explain: 'Storey height over relative sway of the worst storey; higher is stiffer.' },
        { key: 'rot', label: 'Largest joint rotation', unit: 'mrad', value: rot * 1e3, limit: inp.rotLimit, sense: 'max', explain: 'Largest slope at a joint (rotation about a horizontal axis).' },
      ],
      checks,
      plots: [
        framePlot({ title: 'Vertical deflection and deformed shape', ref, u, values: vert.map((v) => v * 1e3), valueLabel: 'Downward deflection', valueUnit: 'mm', explain: 'Frame coloured by downward deflection under characteristic loads. The largest values should sit mid-way between leg lines or under the heaviest equipment.' }),
        { type: 'line', title: 'Lateral displacement over height', x: { label: 'Mean lateral displacement', unit: 'mm' }, y: { label: 'Elevation above the interface', unit: 'm' }, series: [{ name: 'Level sway', x: [0, ...storeys.map((s) => round(Math.hypot(s.ux, s.uy) * 1e3))], y: [0, ...zs] }, { name: `Limit H/${inp.driftRatio}`, x: [0, round(limL * 1e3)], y: [0, H], dash: true }],
          explain: 'Sway profile of the deck levels against the drift limit line. A kink marks a soft storey; a braced frame shows a nearly straight, steep line.' },
        { type: 'bar', title: 'Serviceability utilisation', categories: checks.map((c) => c.name), series: [{ name: 'Utilisation', y: checks.map((c) => round(c.uc)) }], y: { label: 'Demand / limit', unit: '-' }, explain: 'Each deflection criterion as a fraction of its limit; anything above 1.0 needs stiffening or a camber.' },
      ],
      tables: [{ title: 'Storey displacements', columns: [{ key: 'level', label: 'Level' }, { key: 'z', label: 'Elevation', unit: 'm' }, { key: 'ux', label: 'ux', unit: 'mm' }, { key: 'uy', label: 'uy', unit: 'mm' }, { key: 'drift', label: 'Storey drift', unit: 'mm' }, { key: 'ratio', label: 'h / Δ' }],
        rows: storeys.map((s) => ({ level: s.level, z: round(s.z), ux: round(s.ux * 1e3), uy: round(s.uy * 1e3), drift: round(s.drift * 1e3), ratio: Number.isFinite(s.ratio) ? round(s.ratio) : '∞' })), explain: 'Mean horizontal displacement of each deck level and the drift of the storey below it. Use the h/Δ column against the drift criterion.' }],
      notes: [`${gov.name} governs at ${fmt(gov.uc)} of its limit.`, `Largest lateral displacement at node (${fmt(P[latNode][0])}, ${fmt(P[latNode][1])}, ${fmt(P[latNode][2])}) m.`, 'Loads are characteristic (all factors 1.0).'],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `All serviceability criteria are met; ${gov.name.toLowerCase()} uses ${fmt(gov.uc * 100)} % of its limit (deck sag ${fmt(dRel * 1e3)} mm = span/${fmt(spanGlobal / Math.max(dRel, 1e-9), 3)}).`,
          'Specify a fabrication camber equal to the dead-load deflection for the long girders, and give the equipment vendor the predicted slope at its supports.', 'ISO 19901-3:2014 §7.5; EN 1990:2002 A1.4', `Stiffness, not strength, often sizes long deck girders: the beams and girders hold ${steelCo2Text(beamT)}`)
          : adv('engineering', 'fail', `${gov.name} exceeds its limit: ${fmt(gov.demand)} mm (or mrad) against ${fmt(gov.capacity)} (utilisation ${fmt(gov.uc)}).`,
            gov.name.startsWith('Lateral') || gov.name.startsWith('Storey') ? 'Add vertical bracing in the perimeter faces or increase the leg and column diameter; bracing is far more effective than heavier moment-frame members.' : `Increase the girder and beam depth by about ${fmt((depthFac - 1) * 100, 2)} % at the same area, add an intermediate leg or column line, or pre-camber for the dead-load part.`,
            'ISO 19901-3:2014 §7.5; EN 1990:2002 A1.4', `Deeper sections at constant area cost no steel; a ${fmt((depthFac ** 2 - 1) * 100, 2)} % stiffness gain by added area instead would need about ${steelCo2Text(beamT * (depthFac ** 2 - 1))}`),
        adv('economic', 'info', `Deflection control is achieved here with ${fmt(beamT)} t of girders and beams; the governing criterion is at ${fmt(gov.uc * 100)} % of its limit.`,
          'Where deflection governs, camber and depth are cheaper than tonnes: a cambered girder is a fabrication setting, while added flange area is paid for in steel, welding and lift weight.', 'ISO 19901-5:2021 weight control; whole-life costing per ISO 15686-5'),
      ],
    };
  },
  verify: [{
    // One bay, one level, a 2000 t block at the plan centre spread equally to the four corner nodes above the legs,
    // self-weight off: the deck does not bend and every node settles by the elastic shortening of the legs,
    // δ = P h/(E A) with P = 2000e3·9.80665/4 = 4.9033 MN, h = 4 m, A = π(0.914 − 0.032)·0.032 = 0.088668 m²,
    // E = 210 GPa → δ = 4.9033e6·4/(210e9·0.088668) = 1.0533 mm; relative deflection and sag are zero.
    name: 'Axial shortening of the deck legs: δ = P h / (E A)', source: 'Hooke’s law for a bar (hand calculation)',
    inputs: { swFac: 0, live: 0, Hx: 0, Hy: 0, tPlate: 0, secSteel: 0, nx: 1, ny: 1, nLevels: 1, bracing: 'none', columns: 'corners', equip: [{ name: 'Block', W: 2000, x: 15, y: 10, level: 1, a: 0, b: 0, zc: 0 }] },
    expect: [{ kpi: 'dAbs', value: 1.0533, tol: 0.002 }, { kpi: 'dRel', value: 0, abs: 1e-6 }, { kpi: 'sag', value: 0, abs: 1e-6 }],
  }],
};
