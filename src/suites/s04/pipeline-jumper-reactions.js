// Pipeline and jumper reaction loads: connector loads resolved to the foundation for operating, test and shut-in conditions.
import { mudmatClassicalCapacity } from '../../engine/geotech/mudmat.js';
import { fmt, round, need, adv, check, steelCo2Text, vesselCo2Text, matInputs, soilInputs, soilFrom, CONFIRM } from './_common.js';

const DEFAULT_CONN = [
  { name: 'Flowline hub A (10 in)', x: 5.0, y: 2.5, z: 1.5, Fx: -60, Fy: 15, Fz: 25, Mx: 10, My: 40, Mz: 20, th: 1 },
  { name: 'Flowline hub B (10 in)', x: 5.0, y: -2.5, z: 1.5, Fx: -40, Fy: -10, Fz: 25, Mx: -8, My: 30, Mz: -15, th: 1 },
  { name: 'Well jumper 1 (6 in)', x: -5.0, y: 3.5, z: 2.0, Fx: 20, Fy: 8, Fz: 15, Mx: 12, My: -18, Mz: 10, th: 0 },
  { name: 'Well jumper 2 (6 in)', x: -5.0, y: -3.5, z: 2.0, Fx: 20, Fy: -8, Fz: 15, Mx: -12, My: -18, Mz: -10, th: 0 },
  { name: 'Umbilical termination', x: 0, y: 7.0, z: 1.2, Fx: 0, Fy: -12, Fz: 8, Mx: 5, My: 0, Mz: 0, th: 0 },
];

export default {
  id: 's04-pipeline-jumper-reactions', suite: 4, title: 'Pipeline and jumper reaction loads',
  summary: 'Resolves pipeline, jumper and umbilical connector loads, including thermal expansion reactions, into foundation loads for operating, hydrotest and shut-in conditions and checks connectors and foundation.',
  covers: [11], fidelity: 'analytical',
  standards: ['ISO 13628-1:2005 / API RP 17A §5.5 (interface loads)', 'ISO 13628-4 / API 17D (connector load capacity)', 'ISO 13628-15 / API 17P (manifold structures)', 'DNV-ST-F101 (2021) §5 (pipeline end expansion)', 'ISO 19901-4:2016 §7 (foundation check)'],
  models: ['Force and moment equilibrium (load transfer to mudline)', 'Linear jumper / spool stiffness for thermal expansion reactions', 'Meyerhof effective-area method', 'Torsion-equivalent horizontal load', 'Brinch Hansen bearing capacity'],
  limits: 'Connector loads are user values per condition scaled by factors; the expansion reaction is linear (spool stiffness × end displacement) and no pipe–soil or spool flexibility analysis is performed here. The structure is rigid between hubs and foundation. Connector capacity is a linear force–moment interaction, not a manufacturer capacity chart.',
  publishes: ['subsea.designV', 'subsea.designH', 'subsea.designM', 'subsea.designT'],
  searchTerms: 'jumper spool connector loads manifold foundation thermal expansion reaction',
  inputs: [
    { key: 'conns', label: 'Connection loads in operation', type: 'table', group: 'Interface loads', role: 'bc',
      columns: [{ key: 'name', label: 'Connection', type: 'text' }, { key: 'x', label: 'x', unit: 'm', type: 'number' }, { key: 'y', label: 'y', unit: 'm', type: 'number' }, { key: 'z', label: 'z above mudline', unit: 'm', type: 'number' },
        { key: 'Fx', label: 'Fx', unit: 'kN', type: 'number' }, { key: 'Fy', label: 'Fy', unit: 'kN', type: 'number' }, { key: 'Fz', label: 'Fz (down +)', unit: 'kN', type: 'number' },
        { key: 'Mx', label: 'Mx', unit: 'kN·m', type: 'number' }, { key: 'My', label: 'My', unit: 'kN·m', type: 'number' }, { key: 'Mz', label: 'Mz', unit: 'kN·m', type: 'number' }, { key: 'th', label: 'Takes expansion (1/0)', type: 'number' }],
      default: DEFAULT_CONN, help: 'Loads from the connected line on the structure at each hub, in structure axes (x along B, y along L, origin at the mat centre on the mudline). Rows flagged 1 also receive the expansion reaction along −x.' },
    { key: 'Fexp', label: 'Pipeline end expansion force', unit: 'kN', type: 'number', default: 150, min: 0, max: 1e5, group: 'Thermal expansion', role: 'bc', bind: 'pipeline.endForce', uncertain: { dist: 'lognormal', cov: 0.25 },
      help: 'Upper bound of the axial force the line can deliver (fully restrained end). The reaction used is the lesser of this and stiffness × expansion, shared by the flagged connections.' },
    { key: 'dExp', label: 'Pipeline end expansion', unit: 'm', type: 'number', default: 0.25, min: 0, max: 5, group: 'Thermal expansion', role: 'bc', bind: 'pipeline.endExpansion' },
    { key: 'kSpool', label: 'Spool / jumper stiffness in the expansion direction', unit: 'kN/m', type: 'number', default: 400, min: 0, max: 1e6, group: 'Thermal expansion', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Force per unit end displacement of the tie-in spool; flexible M- or Z-spools are typically a few hundred kN/m.' },
    { key: 'fTest', label: 'Hydrotest factor on connection loads', unit: '-', type: 'number', default: 1.25, min: 0, max: 3, group: 'Conditions', role: 'bc', help: 'Pressure-proportional scaling of the tabulated loads at system pressure test, without thermal expansion.' },
    { key: 'fShut', label: 'Shut-in factor on connection loads', unit: '-', type: 'number', default: 0.6, min: 0, max: 3, group: 'Conditions', role: 'bc' },
    { key: 'fCold', label: 'Shut-in factor on the expansion reaction', unit: '-', type: 'number', default: -0.5, min: -1.5, max: 1.5, group: 'Conditions', role: 'bc', help: 'Negative when the line contracts on cool-down and pulls the structure the other way (residual lock-in).' },
    { key: 'W', label: 'Submerged weight of structure and equipment', unit: 'kN', type: 'number', default: 1000, min: 1, max: 2e5, group: 'Structure', role: 'ic', bind: 'subsea.submergedWeight' },
    { key: 'cx', label: 'Centre of gravity offset x', unit: 'm', type: 'number', default: 0.2, min: -30, max: 30, group: 'Structure', role: 'ic', bind: 'subsea.cogX' },
    { key: 'cy', label: 'Centre of gravity offset y', unit: 'm', type: 'number', default: 0.3, min: -40, max: 40, group: 'Structure', role: 'ic', bind: 'subsea.cogY' },
    { key: 'Henv', label: 'Coexisting environmental load (x)', unit: 'kN', type: 'number', default: 60, min: 0, max: 1e5, group: 'Structure', role: 'bc', bind: 'subsea.hydroFH' },
    { key: 'zEnv', label: 'Lever arm of the environmental load', unit: 'm', type: 'number', default: 1.8, min: 0, max: 40, group: 'Structure', role: 'bc', bind: 'subsea.hydroArm' },
    ...matInputs(),
    { key: 'D', label: 'Skirt depth', unit: 'm', type: 'number', default: 0.5, min: 0, max: 10, group: 'Mudmat geometry', role: 'geometry' },
    ...soilInputs(),
    { key: 'Fallow', label: 'Connector allowable force', unit: 'kN', type: 'number', default: 300, min: 1, max: 1e5, group: 'Acceptance criteria', role: 'criteria', help: 'Resultant external force capacity of the connector at design pressure, from the manufacturer capacity chart.' },
    { key: 'Mallow', label: 'Connector allowable bending moment', unit: 'kN·m', type: 'number', default: 250, min: 1, max: 1e5, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'gammaF', label: 'Load factor on interface and environmental loads', unit: '-', type: 'number', default: 1.3, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'gammaM', label: 'Material factor on soil strength', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
  ],
  compute(inp) {
    // SI: kN → N and kN·m → N·m are applied at the foundation check; sums are carried in kN and kN·m.
    const conns = (inp.conns || []).filter((c) => c.name !== undefined);
    need(conns.length > 0, 'Give at least one connection.');
    const nTh = conns.filter((c) => c.th >= 1).length;
    const Fth = Math.min(inp.Fexp, inp.kSpool * inp.dExp);
    const s = soilFrom(inp);
    const { B, L, D } = inp;
    const su = s.suAt(D + Math.min(B, L) / 4);
    const conds = [{ name: 'Operating', fT: 1, fE: 1 }, { name: 'Hydrotest', fT: inp.fTest, fE: 0 }, { name: 'Shut-in / cool-down', fT: inp.fShut, fE: inp.fCold }];
    const results = conds.map((c) => {
      let Fx = inp.Henv, Fy = 0, Fz = 0, Mx = 0, My = inp.Henv * inp.zEnv, Mz = 0;
      const per = conns.map((k) => {
        const fx = c.fT * (k.Fx || 0) - (k.th >= 1 && nTh ? c.fE * Fth / nTh : 0), fy = c.fT * (k.Fy || 0), fzD = c.fT * (k.Fz || 0);
        const mx = c.fT * (k.Mx || 0), my = c.fT * (k.My || 0), mz = c.fT * (k.Mz || 0);
        // Moment about the mat centre at mudline: r × F + M, with F = (fx, fy, −fzD).
        Fx += fx; Fy += fy; Fz += fzD;
        Mx += k.y * -fzD - k.z * fy + mx; My += k.z * fx - k.x * -fzD + my; Mz += k.x * fy - k.y * fx + mz;
        const F = Math.hypot(fx, fy, fzD), M = Math.hypot(mx, my, mz);
        return { name: k.name, F, M, uc: F / inp.Fallow + M / inp.Mallow };
      });
      const g = inp.gammaF;
      // Weight (unfactored) at its centre of gravity adds vertical load and moment.
      const V = inp.W + g * Fz;
      const MB = g * My + inp.W * inp.cx, ML = g * Mx - inp.W * inp.cy;
      const H = g * Math.hypot(Fx, Fy), T = g * Mz;
      let ucB = 99, ucS = 99, err = null;
      if (V > 0) {
        try {
          const r = mudmatClassicalCapacity({ B, L, V: V * 1e3, HB: Math.abs(g * Fx) * 1e3, HL: Math.abs(g * Fy) * 1e3, MB: Math.abs(MB) * 1e3, ML: Math.abs(ML) * 1e3, T: Math.abs(T) * 1e3,
            su: s.clay ? su : undefined, phi: s.phi, c: s.c, gamma: s.gamma, D, materialFactor: inp.gammaM });
          ucB = r.ucBearing; ucS = r.ucSliding;
        } catch (e) { err = e.message; }
      } else err = 'net uplift on the foundation';
      const connMax = per.reduce((a, p) => (p.uc > a.uc ? p : a));
      return { ...c, Fx, Fy, Fz, V, H, MB, ML, M: Math.hypot(MB, ML), T, ucB: Math.min(ucB, 99), ucS: Math.min(ucS, 99), err, per, connMax };
    });
    const score = (r) => Math.max(r.ucB, r.ucS);
    const gov = results.reduce((a, r) => (score(r) > score(a) ? r : a));
    const govC = results.reduce((a, r) => (r.connMax.uc > a.connMax.uc ? r : a));
    const checks = [];
    for (const r of results) {
      checks.push({ name: `Foundation bearing — ${r.name}`, demand: round(r.ucB), capacity: 1, unit: 'utilisation', uc: r.ucB, ref: 'ISO 19901-4:2016 §7.3 (effective area, Brinch Hansen factors)', explain: 'Guards against bearing failure of the mat under the vertical load, eccentricity and inclination produced by the interface loads in this condition.' });
      checks.push({ name: `Foundation sliding — ${r.name}`, demand: round(r.ucS), capacity: 1, unit: 'utilisation', uc: r.ucS, ref: 'ISO 19901-4:2016 §7.4; DNV-RP-C212 §5.4 (torsion-equivalent load)', explain: 'Guards against the structure being dragged or twisted along the seabed by pipeline and jumper pulls.' });
    }
    checks.push(check(`Connector interaction — ${govC.connMax.name} (${govC.name})`, govC.connMax.uc, 1, 'F/Fallow + M/Mallow', 'ISO 13628-4 / API 17D connector external load capacity (linear interaction, screening)', 'Guards against overloading the connector, which first shows as loss of preload and seal leakage.'));
    const ok = checks.every((c) => c.uc <= 1);
    const hubs = conns.map((c) => c.name);
    const dSpool = inp.kSpool * inp.dExp > inp.Fexp;
    const advice = [
      ok ? adv('engineering', 'ok', `Governing condition for the foundation is ${gov.name}: V = ${fmt(gov.V)} kN, H = ${fmt(gov.H)} kN, M = ${fmt(gov.M)} kN·m, T = ${fmt(gov.T)} kN·m (bearing ${fmt(gov.ucB)}, sliding ${fmt(gov.ucS)}). Expansion reaction ${fmt(Fth)} kN.`,
        'Replace the tabulated hub loads with the results of the spool flexibility analysis for each condition, and confirm connector capacity against the manufacturer chart at the coexisting pressure.',
        'ISO 13628-1 / API RP 17A §5.5; ISO 13628-4 / API 17D; ISO 19901-4:2016 §7',
        `A more flexible spool reduces the expansion reaction in proportion; halving the stiffness would cut ${fmt(Fth / 2)} kN of horizontal load and the mat steel that resists it.`)
        : adv('engineering', 'fail', `${checks.reduce((a, c) => (c.uc > a.uc ? c : a)).name} is exceeded. Governing foundation condition ${gov.name}: bearing ${fmt(gov.ucB)}, sliding ${fmt(gov.ucS)}${gov.err ? ` (${gov.err})` : ''}; highest connector interaction ${fmt(govC.connMax.uc)}.`,
          govC.connMax.uc > 1 ? `Reduce the load at "${govC.connMax.name}": add flexibility to the spool (the expansion reaction is ${fmt(Fth)} kN at ${fmt(inp.kSpool)} kN/m), or select a connector rated for ${fmt(govC.connMax.F)} kN and ${fmt(govC.connMax.M)} kN·m.` : `Lower the spool stiffness to about ${fmt(inp.kSpool / Math.max(score(gov), 1.01))} kN/m, re-balance opposing hubs to cancel torsion (${fmt(gov.T)} kN·m), or enlarge the mat / add skirts.`,
          'ISO 13628-1 / API RP 17A §5.5; ISO 13628-4 / API 17D; ISO 19901-4:2016 §7',
          `Spool flexibility is the low-material fix: a few metres of additional pipe (≈ 0.1 t/m for 10 in) against a mat enlargement of several tonnes — ${steelCo2Text(5)}`),
      adv('economic', ok ? 'info' : 'warn', `The expansion reaction (${fmt(Fth)} kN, ${dSpool ? 'capped by the line end force' : 'set by spool stiffness'}) accounts for ${fmt(Fth / (Math.abs(results[0].Fx) + 1e-9) * 100)} % of the operating horizontal load.`,
        'Fix spool geometry and structure together: metrology and spool fabrication are on the installation critical path, and a connector overload found late means a new spool (weeks) or an extra tie-in campaign (about 2 vessel days).',
        'AACE RP 18R-97 estimate classes; schedule-risk practice for subsea tie-ins',
        vesselCo2Text(2)),
    ];
    return {
      kpis: [
        { key: 'V', label: 'Design vertical load', unit: 'kN', value: gov.V, explain: 'Submerged weight plus factored vertical connector loads in the governing condition, compression positive.' },
        { key: 'H', label: 'Design horizontal load', unit: 'kN', value: gov.H, explain: 'Resultant factored horizontal load at mudline from connectors, expansion reaction and coexisting environmental load.' },
        { key: 'M', label: 'Design overturning moment', unit: 'kN·m', value: gov.M, explain: 'Resultant factored moment about the mat centre at mudline, including forces times their height and the weight eccentricity.' },
        { key: 'T', label: 'Design torsion', unit: 'kN·m', value: Math.abs(gov.T), explain: 'Factored moment about the vertical axis from horizontal pulls offset from the mat centre; it twists the mat on the seabed.' },
        { key: 'MB', label: 'Moment giving eccentricity along B', unit: 'kN·m', value: Math.abs(gov.MB), explain: 'Component of the design moment that shifts the reaction along the width B.' },
        { key: 'Fth', label: 'Thermal expansion reaction', unit: 'kN', value: Fth, explain: 'Force delivered by the expanding line through the tie-in spools: the lesser of spool stiffness × expansion and the line end force.' },
        { key: 'ucBearing', label: 'Foundation bearing utilisation', unit: '-', value: gov.ucB, limit: 1, sense: 'max', explain: 'Highest bearing utilisation over the three conditions (99 when the loads cannot be carried at all).' },
        { key: 'ucSliding', label: 'Foundation sliding utilisation', unit: '-', value: Math.max(...results.map((r) => r.ucS)), limit: 1, sense: 'max', explain: 'Highest sliding utilisation over the three conditions, torsion included.' },
        { key: 'ucConn', label: 'Highest connector interaction', unit: '-', value: govC.connMax.uc, limit: 1, sense: 'max', explain: 'Largest F/Fallow + M/Mallow over all connections and conditions; above 1 the connector is outside its external-load envelope.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Foundation loads by condition', categories: results.map((r) => r.name),
          series: [{ name: 'V (kN)', y: results.map((r) => round(r.V)) }, { name: 'H (kN)', y: results.map((r) => round(r.H)) }, { name: 'M (kN·m)', y: results.map((r) => round(r.M)) }, { name: '|T| (kN·m)', y: results.map((r) => round(Math.abs(r.T))) }], y: { label: 'Load', unit: 'kN, kN·m' },
          explain: 'Resultant loads at the mat centre for each operating condition. Compare the bars to see which condition drives sliding (H, T) and which drives bearing (V, M).' },
        { type: 'bar', title: 'Connector interaction by connection', categories: hubs, series: results.map((r) => ({ name: r.name, y: r.per.map((p) => round(p.uc)) })), y: { label: 'F/Fallow + M/Mallow', unit: '-' },
          explain: 'Utilisation of each connector in each condition. Bars above 1.0 identify the hub whose spool needs more flexibility.' },
        { type: 'line', title: 'Expansion reaction against spool stiffness', x: { label: 'Spool stiffness', unit: 'kN/m' }, y: { label: 'Reaction', unit: 'kN' },
          series: [{ name: 'Reaction', x: Array.from({ length: 21 }, (_, i) => round(inp.kSpool * i / 10)), y: Array.from({ length: 21 }, (_, i) => Math.min(inp.Fexp, inp.kSpool * i / 10 * inp.dExp)) }, { name: 'Present design', x: [inp.kSpool], y: [Fth] }],
          explain: 'Force delivered to the structure as the tie-in spool is made stiffer or softer. The plateau is the fully restrained line end force.' },
      ],
      tables: [
        { title: 'Foundation loads and utilisations', columns: [{ key: 'name', label: 'Condition' }, { key: 'V', label: 'V', unit: 'kN' }, { key: 'Fx', label: 'ΣFx (unfactored)', unit: 'kN' }, { key: 'Fy', label: 'ΣFy (unfactored)', unit: 'kN' }, { key: 'H', label: 'H', unit: 'kN' }, { key: 'MB', label: 'M (along B)', unit: 'kN·m' }, { key: 'ML', label: 'M (along L)', unit: 'kN·m' }, { key: 'T', label: 'T', unit: 'kN·m' }, { key: 'ucB', label: 'Bearing' }, { key: 'ucS', label: 'Sliding' }],
          rows: results.map((r) => ({ name: r.name, V: round(r.V), Fx: round(r.Fx), Fy: round(r.Fy), H: round(r.H), MB: round(r.MB), ML: round(r.ML), T: round(r.T), ucB: round(r.ucB), ucS: round(r.ucS) })),
          explain: 'One row per condition: the factored resultants at the mat centre on the mudline and the bearing and sliding utilisations they produce.' },
        { title: 'Connector loads in the governing condition', columns: [{ key: 'name', label: 'Connection' }, { key: 'F', label: 'Resultant force', unit: 'kN' }, { key: 'M', label: 'Resultant moment', unit: 'kN·m' }, { key: 'uc', label: 'Interaction' }],
          rows: govC.per.map((p) => ({ name: p.name, F: round(p.F), M: round(p.M), uc: round(p.uc) })),
          explain: `Unfactored resultant force and moment at each hub in the ${govC.name} condition, for comparison with the connector capacity chart.` },
      ],
      notes: [`${gov.name} governs the foundation; ${govC.name} governs the connectors.`, `Expansion reaction ${fmt(Fth)} kN shared by ${nTh} connection(s), acting along −x.`, gov.err ? `Foundation check could not be completed in the governing condition: ${gov.err}.` : 'All conditions could be evaluated on the effective area.'],
      advice,
      publish: { 'subsea.designV': gov.V * 1e3, 'subsea.designH': gov.H * 1e3, 'subsea.designM': gov.M * 1e3, 'subsea.designT': Math.abs(gov.T) * 1e3 },
      quantities: {},
    };
  },
  verify: [
    { name: 'Single horizontal pull offset from the centre: M = F·z, T = F·y', source: 'Statics: Fx = 100 kN at (0, 3, 2) m → My = 200 kN·m, Mz = −300 kN·m, H = 100 kN, V = W = 1000 kN',
      inputs: { conns: [{ name: 'hub', x: 0, y: 3, z: 2, Fx: 100, Fy: 0, Fz: 0, Mx: 0, My: 0, Mz: 0, th: 0 }], Fexp: 0, W: 1000, cx: 0, cy: 0, Henv: 0, gammaF: 1, fTest: 0.5, fShut: 0.5, fCold: 0 },
      expect: [{ kpi: 'M', value: 200, tol: 1e-9 }, { kpi: 'T', value: 300, tol: 1e-9 }, { kpi: 'H', value: 100, tol: 1e-9 }, { kpi: 'V', value: 1000, tol: 1e-9 }] },
    { name: 'Expansion reaction: F = min(k·δ, F_end)', source: 'Linear spool: 400 kN/m × 0.25 m = 100 kN < 150 kN',
      inputs: { kSpool: 400, dExp: 0.25, Fexp: 150 }, expect: [{ kpi: 'Fth', value: 100, tol: 1e-9 }] },
  ],
};
