// Decommissioning and removal lifts: probabilistic lift weight of an aged structure (corrosion loss, marine growth,
// grout, trapped water, debris), hook load against crane capacity, integrity of the corroded lifting region and
// the loading of each lift in the removal sequence.
import { propagate } from '../../engine/reliability/sampling.js';
import { G, DEG, need, run, round, fmt, kpi, check, adv, bar, tbl, worst, tube, dafN001, LIFT_SITES, vesselCo2Text, CO2_STEEL, CONFIRM, N001, LIMIT_N001 } from './_common.js';

const DEFAULT_SEQ = [
  { name: 'Topside module M1', weight: 2200, aged: 0, capacity: 3300 }, { name: 'Topside module M2', weight: 1800, aged: 0, capacity: 3300 }, { name: 'Module support frame', weight: 1100, aged: 0, capacity: 2500 },
  { name: 'Jacket upper section', weight: 2600, aged: 1, capacity: 4200 }, { name: 'Jacket lower section', weight: 2900, aged: 1, capacity: 5000 },
];

export default {
  id: 's15-decommissioning', suite: 15, title: 'Decommissioning and removal lifts',
  summary: 'Estimates the removal weight of an aged structure with its uncertainty, the dynamic hook load against crane capacity, the strength of the corroded structure at the lift points, and the loading of every lift in the removal sequence.',
  covers: [18], fidelity: 'analytical',
  standards: [`${N001} Sec. 16 (lifting) and Sec. 17 (removal operations: weight and CoG uncertainty, structural condition)`, 'ISO 19901-6 (2009) §19 (decommissioning and removal)', 'ISO 19902 (2020) §23–24 (assessment of existing structures, corroded members)', 'OSPAR Decision 98/3 (disposal of disused offshore installations)', 'Latin hypercube sampling: McKay, Beckman & Conover (1979)'],
  models: ['Probabilistic weight model (Latin hypercube Monte Carlo)', 'Corrosion-informed residual section (uniform wall loss)', 'Temporary load redistribution to the lift points (rigid body, skew factor)', 'Dynamic lifting (DAF by weight class)', 'Axial plus bending stress of the corroded member at the lift point', 'Operation sequencing of the removal lifts'],
  limits: `${LIMIT_N001} Weight components are independent random variables with the distributions entered; corrosion is uniform wall loss without pitting, cracks or dented members, and the lift-point check is a nominal section check, not a damaged-structure finite-element analysis. Cutting, soil breakout of pile stubs and grout bond are represented by allowances. The seeded Monte Carlo is reproducible but only as good as the survey data behind the inputs.`,
  publishes: ['install.removalWeightP90', 'install.removalHookLoad'],
  searchTerms: 'decommissioning removal lift weight uncertainty marine growth corrosion residual strength single lift',
  inputs: [
    { key: 'Wdry', label: 'As-built steel weight of the piece (weight report)', unit: 't', type: 'number', default: 2600, min: 5, max: 60000, group: 'Removal weight', role: 'ic', help: 'Initial condition of the assessment: the documented weight before ageing effects.' },
    { key: 'covW', label: 'Uncertainty of the as-built weight (CoV)', unit: '%', type: 'number', default: 5, min: 0, max: 30, group: 'Removal weight', role: 'load', help: 'Old weight reports and undocumented modifications: 3–10 %.' },
    { key: 'corr', label: 'Mean steel loss by corrosion', unit: '%', type: 'number', default: 4, min: 0, max: 40, group: 'Removal weight', role: 'ic', help: 'Mass lost from the as-built steel; from thickness surveys.' },
    { key: 'covCorr', label: 'Uncertainty of the corrosion loss (CoV)', unit: '%', type: 'number', default: 40, min: 0, max: 150, group: 'Removal weight', role: 'load' },
    { key: 'Amg', label: 'Surface carrying marine growth', unit: 'm²', type: 'number', default: 5000, min: 0, max: 5e5, group: 'Removal weight', role: 'geometry' },
    { key: 'tmg', label: 'Mean marine growth thickness', unit: 'mm', type: 'number', default: 40, min: 0, max: 400, group: 'Removal weight', role: 'ic' },
    { key: 'covMg', label: 'Uncertainty of marine growth (CoV)', unit: '%', type: 'number', default: 35, min: 0, max: 150, group: 'Removal weight', role: 'load' },
    { key: 'rhoMg', label: 'Marine growth density in air (drained)', unit: 'kg/m³', type: 'number', default: 1325, min: 800, max: 1600, group: 'Removal weight', role: 'load', help: 'Wet growth lifted through the surface keeps most of its water; hard growth is denser.' },
    { key: 'grout', label: 'Grout, pile stubs and other attached mass', unit: 't', type: 'number', default: 180, min: 0, max: 20000, group: 'Removal weight', role: 'load' },
    { key: 'covGrout', label: 'Uncertainty of attached mass (CoV)', unit: '%', type: 'number', default: 20, min: 0, max: 100, group: 'Removal weight', role: 'load' },
    { key: 'Vflood', label: 'Volume of members that may hold water', unit: 'm³', type: 'number', default: 250, min: 0, max: 50000, group: 'Removal weight', role: 'geometry', help: 'Flooded or possibly flooded members and compartments that are not confirmed to drain during the lift.' },
    { key: 'drainMin', label: 'Share of that water retained: least', unit: '%', type: 'number', default: 10, min: 0, max: 100, group: 'Removal weight', role: 'load' },
    { key: 'drainMax', label: 'Share of that water retained: most', unit: '%', type: 'number', default: 70, min: 0, max: 100, group: 'Removal weight', role: 'load' },
    { key: 'debris', label: 'Debris, sediment and drill cuttings: most likely', unit: 't', type: 'number', default: 40, min: 0, max: 10000, group: 'Removal weight', role: 'load' },
    { key: 'debrisMax', label: 'Debris, sediment and drill cuttings: upper bound', unit: 't', type: 'number', default: 160, min: 0, max: 20000, group: 'Removal weight', role: 'load' },
    { key: 'fCont', label: 'Contingency factor on the P90 weight', unit: '-', type: 'number', default: 1.05, min: 1, max: 1.5, group: 'Criteria', role: 'criteria', help: `Covers what the survey cannot see. ${CONFIRM}` },
    { key: 'Wrig', label: 'Rigging and lifting-tool weight', unit: 't', type: 'number', default: 150, min: 0, max: 3000, group: 'Lift', role: 'load' },
    { key: 'breakout', label: 'Breakout allowance at lift-off', unit: '%', type: 'number', default: 5, min: 0, max: 50, group: 'Lift', role: 'bc', help: 'Soil suction, friction in guides and incomplete cuts at the moment of lift-off, as a share of the weight.' },
    { key: 'where', label: 'Lift location', type: 'select', default: 'offshore', options: LIFT_SITES, group: 'Lift', role: 'operation' },
    { key: 'dafUser', label: 'DAF override (0 = table by weight class)', unit: '-', type: 'number', default: 0, min: 0, max: 3, group: 'Criteria', role: 'criteria', help: CONFIRM },
    { key: 'craneCap', label: 'Crane capacity at the removal radius', unit: 't', type: 'number', default: 4200, min: 5, max: 60000, group: 'Lift', role: 'bc', design: true },
    { key: 'nPts', label: 'Number of lift points', unit: '-', type: 'int', default: 4, min: 1, max: 16, group: 'Lift points on the aged structure', role: 'bc' },
    { key: 'skl', label: 'Skew / load-sharing factor on lift-point loads', unit: '-', type: 'number', default: 1.25, min: 1, max: 2, group: 'Criteria', role: 'criteria', help: `1.25 for four slings to one hook without a load-equalising system; higher if the CoG is uncertain. ${CONFIRM}` },
    { key: 'slingDeg', label: 'Sling angle to the horizontal', unit: 'deg', type: 'number', default: 65, min: 30, max: 90, group: 'Lift points on the aged structure', role: 'geometry' },
    { key: 'Dleg', label: 'Leg / member diameter at the lift point', unit: 'mm', type: 'number', default: 1500, min: 100, max: 6000, group: 'Lift points on the aged structure', role: 'geometry' },
    { key: 'tLeg', label: 'As-built wall thickness', unit: 'mm', type: 'number', default: 45, min: 5, max: 150, group: 'Lift points on the aged structure', role: 'geometry' },
    { key: 'tLoss', label: 'Measured wall loss at the lift point', unit: 'mm', type: 'number', default: 6, min: 0, max: 100, group: 'Lift points on the aged structure', role: 'ic', uncertain: { dist: 'lognormal', cov: 0.4 }, help: 'Actual structural condition from the thickness survey; splash-zone members often lose 0.3–0.5 mm/year where coating has failed.' },
    { key: 'ecc', label: 'Eccentricity of the lifting attachment from the member axis', unit: 'm', type: 'number', default: 0.6, min: 0, max: 10, group: 'Lift points on the aged structure', role: 'geometry' },
    { key: 'mat', label: 'Steel grade of the existing structure', type: 'select', source: 'materials:steel', default: 's355', group: 'Lift points on the aged structure', role: 'material' },
    { key: 'gM', label: 'Material factor for the aged steel', unit: '-', type: 'number', default: 1.25, min: 1, max: 2, group: 'Criteria', role: 'criteria', help: 'Higher than for new steel when material certificates are missing or the condition is uncertain.' },
    { key: 'fLoad', label: 'Load × consequence factor at lift points', unit: '-', type: 'number', default: 1.69, min: 1, max: 3, group: 'Criteria', role: 'criteria', help: `1.3 × 1.3 for lift points and their supporting members. ${CONFIRM}` },
    { key: 'seq', label: 'Removal sequence', type: 'table', group: 'Removal sequence', role: 'operation',
      columns: [{ key: 'name', label: 'Lift', type: 'text' }, { key: 'weight', label: 'As-built weight', unit: 't', type: 'number' }, { key: 'aged', label: 'Submerged piece (1) or dry (0)', type: 'number' }, { key: 'capacity', label: 'Crane capacity at its radius', unit: 't', type: 'number' }],
      default: DEFAULT_SEQ, help: 'Every lift of the campaign in order. Submerged pieces (1) take the full weight growth found for the assessed piece; dry pieces (0) take only the weight-report uncertainty and the contingency.' },
    { key: 'hoursPerLift', label: 'Vessel time per lift incl. cutting and rigging', unit: 'h', type: 'number', default: 36, min: 1, max: 500, group: 'Removal sequence', role: 'operation' },
    { key: 'nSamples', label: 'Monte Carlo samples', unit: '-', type: 'int', default: 2000, min: 100, max: 50000, group: 'Numerics', role: 'control' },
    { key: 'seed', label: 'Random seed', unit: '-', type: 'int', default: 11, min: 1, max: 99999, group: 'Numerics', role: 'control' },
  ],
  compute(inp, ctx) {
    // % → fraction, mm → m, t kept as tonnes in the weight model.
    const steel = ctx.material(inp.mat); need(steel, `Steel '${inp.mat}' is not in the materials library.`);
    need(inp.drainMax >= inp.drainMin, 'The largest retained-water share must not be below the least.');
    need(inp.debrisMax >= inp.debris, 'The upper bound of debris must not be below the most likely value.');
    const c = (pct) => pct / 100;
    // Random variables only where an uncertainty is entered; the rest are constants.
    const vars = {}, fixed = {};
    const put = (k, mean, cov, dist = 'normal') => { if (cov > 0 && mean > 0) vars[k] = { dist, mean, std: mean * cov }; else fixed[k] = mean; };
    put('W', inp.Wdry, c(inp.covW)); put('corr', c(inp.corr), c(inp.covCorr), 'lognormal'); put('tmg', inp.tmg / 1e3, c(inp.covMg), 'lognormal'); put('grout', inp.grout, c(inp.covGrout));
    if (inp.drainMax > inp.drainMin) vars.ret = { dist: 'uniform', min: c(inp.drainMin), max: c(inp.drainMax) }; else fixed.ret = c(inp.drainMin);
    if (inp.debrisMax > inp.debris) vars.debris = { dist: 'triangular', min: 0, mode: inp.debris, max: inp.debrisMax }; else fixed.debris = inp.debris;
    const parts = (x) => { const v = { ...fixed, ...x }; const steelW = Math.max(v.W, 0) * (1 - Math.min(Math.max(v.corr, 0), 0.9)); return { steel: steelW, mg: inp.Amg * Math.max(v.tmg, 0) * inp.rhoMg / 1e3, grout: Math.max(v.grout, 0), water: inp.Vflood * 1.025 * Math.min(Math.max(v.ret, 0), 1), debris: Math.max(v.debris, 0) }; };
    const total = (x) => { const p = parts(x); return p.steel + p.mg + p.grout + p.water + p.debris; };
    const meanX = Object.fromEntries(Object.entries(vars).map(([k, d]) => [k, d.dist === 'uniform' ? (d.min + d.max) / 2 : d.dist === 'triangular' ? (d.min + d.mode + d.max) / 3 : d.mean]));
    const pm = parts(meanX), Wnom = pm.steel + pm.mg + pm.grout + pm.water + pm.debris;
    let st, samples = [];
    if (Object.keys(vars).length) {
      const mc = run(() => propagate((x) => ({ W: total(x) }), vars, { n: inp.nSamples, seed: inp.seed }), 'Weight Monte Carlo');
      st = mc.stats.W; samples = Array.from(mc.outputs.W);
    } else st = { mean: Wnom, std: 0, p10: Wnom, p50: Wnom, p90: Wnom, min: Wnom, max: Wnom };
    const Wdes = st.p90 * inp.fCont, SHL = Wdes * (1 + c(inp.breakout)) + inp.Wrig;
    const daf = inp.dafUser > 0 ? inp.dafUser : dafN001(SHL, inp.where), DHL = SHL * daf;
    // Probability that the static hook load (no contingency) exceeds the crane capacity.
    const capW = (inp.craneCap - inp.Wrig) / (1 + c(inp.breakout));
    const pEx = samples.length ? samples.filter((w) => w > capW).length / samples.length : (Wnom > capW ? 1 : 0);
    // Sensitivity: swing of the total weight when each variable moves from its 5th to 95th percentile (others at their means).
    const tor = Object.entries(vars).map(([k, d]) => {
      const lo = d.dist === 'uniform' ? d.min + 0.05 * (d.max - d.min) : d.dist === 'triangular' ? d.min + 0.2 * (d.mode - d.min) : Math.max(d.mean - 1.645 * d.std, 0);
      const hi = d.dist === 'uniform' ? d.min + 0.95 * (d.max - d.min) : d.dist === 'triangular' ? d.max - 0.2 * (d.max - d.mode) : d.mean + 1.645 * d.std;
      const nm = { W: 'As-built weight', corr: 'Corrosion loss', tmg: 'Marine growth', grout: 'Grout and attached mass', ret: 'Retained water', debris: 'Debris and sediment' }[k];
      return { name: nm, low: round(total({ ...meanX, [k]: lo })), high: round(total({ ...meanX, [k]: hi })) };
    }).sort((a, b) => Math.abs(b.high - b.low) - Math.abs(a.high - a.low));
    // Corroded member at the lift point: tension from the vertical sling component, bending from the horizontal one.
    const D = inp.Dleg / 1e3, tRes = (inp.tLeg - inp.tLoss) / 1e3;
    need(tRes > 0.002, 'The measured wall loss leaves no usable wall at the lift point: choose another lift point or reinforce.');
    const sec = tube(D, tRes), sec0 = tube(D, inp.tLeg / 1e3);
    const Fv = DHL * 1e3 * G / inp.nPts * inp.skl, th = inp.slingDeg * DEG, Tsl = Fv / Math.sin(th);
    const sig = inp.fLoad * (Fv / sec.A + (Tsl * Math.cos(th) * inp.ecc) / sec.W), fd = steel.fy / inp.gM;
    const sig0 = inp.fLoad * (Fv / sec0.A + (Tsl * Math.cos(th) * inp.ecc) / sec0.W);
    const dt = D / tRes;
    const growth = Wdes / inp.Wdry; // design removal weight over the as-built weight
    const growthDry = (1 + 1.2816 * c(inp.covW)) * inp.fCont; // P90 of the weight-report uncertainty alone
    const seq = (inp.seq || []).map((r) => ({ name: String(r.name || 'Lift'), w: Number(r.weight), cap: Number(r.capacity), aged: Number(r.aged) !== 0 })).filter((r) => r.w > 0 && r.cap > 0).map((r) => {
      const shl = r.w * (r.aged ? growth : growthDry) * (1 + c(inp.breakout)) + inp.Wrig * Math.min(r.w / inp.Wdry, 1.5);
      const d = inp.dafUser > 0 ? inp.dafUser : dafN001(shl, inp.where);
      return { ...r, shl, daf: d, dhl: shl * d, uc: shl / r.cap };
    });
    const checks = [
      check('Static hook load (P90 weight + contingency) against crane capacity', SHL, inp.craneCap, 't', `${N001} Sec. 16–17 (removal: weight contingency; crane capacity)`, 'An aged structure is heavier than its weight report. If the real weight exceeds the chart the lift must be aborted with the piece already cut free.'),
      check('Corroded member at the lift point: axial + bending', sig / 1e6, fd / 1e6, 'MPa', `${N001} Sec. 17 (structural condition for removal); ISO 19902 (2020) §23 (assessment with measured thickness)`, 'Nominal stress in the corroded section from the lift-point load. Old structure was not designed to hang from these points; thin walls can tear or buckle locally.'),
      check('Probability of exceeding crane capacity (limit 1 %)', pEx * 100, 1, '%', 'Reliability target for the lift weight (project risk criterion; DNV-ST-N001 Sec. 17 requires documented weight uncertainty)', 'Chance that the static hook load from the weight distribution alone is above the crane capacity. It should be negligible before committing to a single lift.'),
      ...seq.filter((r) => r.uc > 1).map((r) => check(`Sequence: ${r.name}`, r.shl, r.cap, 't', `${N001} Sec. 16 (crane capacity)`, 'A lift of the removal sequence exceeds its crane capacity once the weight growth is applied: split the piece or change the crane configuration.')),
    ];
    const gov = worst(checks), ok = gov.uc <= 1;
    const days = Math.max(seq.length, 1) * inp.hoursPerLift / 24;
    const steelRec = seq.length ? seq.reduce((s, r) => s + r.w, 0) * (1 - c(inp.corr)) : pm.steel;
    return {
      kpis: [
        kpi('Wnom', 'Removal weight at mean values', Wnom, 't', 'Corroded steel plus marine growth, grout, retained water and debris with every variable at its mean.'),
        kpi('W50', 'Removal weight P50', st.p50, 't', 'Median of the weight distribution.'),
        kpi('W90', 'Removal weight P90', st.p90, 't', 'Weight not exceeded with 90 % probability: the basis of the design hook load.'),
        kpi('Wstd', 'Standard deviation of the removal weight', st.std, 't', 'Spread of the weight estimate. A large value means a survey (thickness, growth, flooded-member detection) would pay for itself.'),
        kpi('growth', 'Design removal weight / as-built weight', growth, '-', 'How much heavier the piece is taken to be than its weight report, including P90 margin and contingency.'),
        kpi('SHL', 'Static hook load', SHL, 't', 'P90 weight × contingency × breakout allowance, plus rigging.', inp.craneCap, 'max'),
        kpi('daf', 'Dynamic amplification factor', daf, '-', 'By weight class for the location, or the override.'),
        kpi('DHL', 'Dynamic hook load', DHL, 't', 'Static hook load × DAF: the design load of rigging and lift points.'),
        kpi('pEx', 'Probability of exceeding crane capacity', pEx * 100, '%', 'Share of weight samples whose static hook load is above the chart value.', 1, 'max'),
        kpi('sig', 'Stress in the corroded member at the lift point', sig / 1e6, 'MPa', 'Factored axial plus bending stress on the measured residual wall.', fd / 1e6, 'max'),
        kpi('resA', 'Residual section area', 100 * sec.A / sec0.A, '%', 'Corroded over as-built cross-section area at the lift point.'),
        kpi('sigRatio', 'Stress increase from corrosion', sig / sig0, '-', 'Stress on the corroded wall over the stress on the as-built wall for the same load.'),
        kpi('nLifts', 'Lifts in the removal sequence', seq.length, '-', 'Number of crane lifts planned.'),
        kpi('days', 'Vessel time for the removal lifts', days, 'days', 'Lifts × time per lift including cutting and rigging, without weather.'),
      ],
      checks,
      plots: [
        ...(samples.length ? [{ type: 'hist', title: 'Removal weight distribution', values: samples.map((v) => round(v)), bins: 30, x: { label: 'Removal weight', unit: 't' }, y: { label: 'Samples', unit: '-' },
          explain: 'Spread of the lift weight from the uncertain components. Compare its right tail with the crane capacity less rigging and breakout allowance.' }] : []),
        ...(tor.length ? [{ type: 'tornado', title: 'What drives the weight uncertainty', base: round(Wnom), items: tor, x: { label: 'Removal weight', unit: 't' },
          explain: 'Swing of the total weight when one component moves across its likely range with the others at their means. The top bars show which survey would narrow the estimate most.' }] : []),
        bar('Removal sequence: hook load against crane capacity', seq.length ? seq.map((r) => r.name) : ['This lift'], 'Load', 't', seq.length ? [{ name: 'As-built weight', y: seq.map((r) => round(r.w)) }, { name: 'Static hook load', y: seq.map((r) => round(r.shl)) }, { name: 'Crane capacity', y: seq.map((r) => round(r.cap)) }] : [{ name: 'Static hook load', y: [round(SHL)] }, { name: 'Crane capacity', y: [inp.craneCap] }],
          'Each lift of the campaign with the weight growth applied. Any lift whose hook load approaches its capacity bar is a candidate for splitting or for a weight-verification survey.'),
      ],
      tables: [
        tbl('Weight build-up at mean values', [['k', 'Component'], ['v', 'Weight', 't'], ['s', 'Share', '%']], [['Corroded steel', pm.steel], ['Marine growth', pm.mg], ['Grout and attached mass', pm.grout], ['Retained water', pm.water], ['Debris and sediment', pm.debris], ['Total', Wnom]].map(([k, v]) => ({ k, v: round(v), s: round(100 * v / Wnom) })),
          'Where the removal weight comes from. Items other than steel are those a weight report does not contain.'),
        tbl('Removal sequence loading', [['n', 'Lift'], ['w', 'As-built weight', 't'], ['shl', 'Static hook load', 't'], ['daf', 'DAF'], ['dhl', 'Dynamic hook load', 't'], ['cap', 'Crane capacity', 't'], ['uc', 'Utilisation']], seq.map((r) => ({ n: r.name, w: round(r.w), shl: round(r.shl), daf: round(r.daf, 3), dhl: round(r.dhl), cap: round(r.cap), uc: round(r.uc, 3) })),
          'One row per lift in removal order with the design weight growth, breakout allowance and dynamic factor applied; the utilisation is static hook load over crane capacity.'),
      ],
      notes: [`Removal weight P50 ${fmt(st.p50)} t, P90 ${fmt(st.p90)} t against ${inp.Wdry} t as built (design factor ${fmt(growth)}).`, `Lift-point member ${inp.Dleg} × ${fmt(inp.tLeg - inp.tLoss)} mm residual (D/t = ${fmt(dt)}): ${fmt(sig / 1e6)} MPa.`, `"${gov.name}" governs at ${fmt(gov.uc)}.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `The piece can be lifted: static hook load ${fmt(SHL)} t against ${inp.craneCap} t, ${fmt(pEx * 100)} % chance of exceeding the chart, ${fmt(sig / 1e6)} MPa in the corroded lift-point member against ${fmt(fd / 1e6)} MPa.`,
          `${tor.length ? `Narrow the largest uncertainty first ("${tor[0].name}", swing ${fmt(Math.abs(tor[0].high - tor[0].low))} t) by survey before the lift, ` : 'Confirm the inputs by survey before the lift, '}fit load cells for a controlled lift-off with a hold point at ${fmt(0.9 * SHL)} t, and verify the lift points by inspection and, where D/t exceeds about 60, by a local buckling check.`, `${N001} Sec. 17; ISO 19901-6 §19`)
          : adv('engineering', 'fail', `"${gov.name}" is not met (${fmt(gov.uc)}): static hook load ${fmt(SHL)} t against ${inp.craneCap} t, exceedance probability ${fmt(pEx * 100)} %, lift-point stress ${fmt(sig / 1e6)} MPa.`,
            gov.name.startsWith('Corroded') ? `Reinforce the lift point (clamp or doubler sized for ${fmt(Fv * inp.fLoad / 1e3)} kN), move it to sounder steel, or add lift points to share the load.` : gov.name.startsWith('Sequence') ? 'Split the overloaded piece with an additional cut, or re-plan it at a shorter crane radius.' : `Reduce the weight before the lift (remove marine growth: ${fmt(pm.mg)} t; drain or vent members: ${fmt(pm.water)} t; clear debris), split the piece, or use a crane with at least ${fmt(Math.ceil(SHL / 100) * 100)} t at this radius.`,
            `${N001} Sec. 16–17; ISO 19901-6 §19`),
        adv('economic', 'info', `The removal sequence is ${seq.length || 1} lift(s), about ${fmt(days)} vessel days; the weight estimate carries ±${fmt(st.std)} t (one standard deviation).`,
          `A pre-removal survey that halves the weight uncertainty lowers the P90 hook load by roughly ${fmt(0.64 * st.std)} t: compare that with the cost of an extra cut and lift (${fmt(inp.hoursPerLift)} vessel hours) if a piece turns out too heavy offshore.`,
          'Value of information for survey scope; removal priced as vessel days per lift (AACE RP 18R-97); OSPAR 98/3 for the disposal obligation',
          `About ${fmt(steelRec)} t of steel returns to shore for recycling, avoiding up to ${fmt(steelRec * CO2_STEEL)} t CO₂e if it displaces primary steel (screening factor ${CO2_STEEL} t CO₂e/t); the removal campaign itself burns about ${vesselCo2Text(days, 60)}`),
      ],
      publish: { 'install.removalWeightP90': st.p90 * 1e3 * G, 'install.removalHookLoad': DHL * 1e3 * G },
      quantities: { vessel_days: { removalLifts: round(days) } },
    };
  },
  verify: [
    { name: 'Deterministic weight build-up with all uncertainties set to zero', source: 'Hand calculation: steel 2600·(1 − 0.04) = 2496 t; marine growth 9000 m²·0.06 m·1.325 t/m³ = 715.5 t; grout 180 t; water 250·1.025·0.4 = 102.5 t; debris 40 t → 3534.0 t; P50 = P90 = mean',
      inputs: { Amg: 9000, tmg: 60, covW: 0, covCorr: 0, covMg: 0, covGrout: 0, drainMin: 40, drainMax: 40, debrisMax: 40 }, expect: [{ kpi: 'Wnom', value: 2496 + 715.5 + 180 + 102.5 + 40, tol: 1e-9 }, { kpi: 'W90', value: 3534, tol: 1e-9 }, { kpi: 'Wstd', value: 0, abs: 1e-9 }] },
    { name: 'Standard deviation of a sum of independent variables: σ² = Σσᵢ²', source: 'Variance addition with only the as-built weight (5 % of 2600 t = 130 t) and grout (20 % of 180 t = 36 t) uncertain: σ = √(130²·0.96² + 36²) = 129.9 t (Latin hypercube, 5000 samples, within 2 %)',
      inputs: { Amg: 9000, tmg: 60, covCorr: 0, covMg: 0, drainMin: 40, drainMax: 40, debrisMax: 40, nSamples: 5000 }, expect: [{ kpi: 'Wstd', value: Math.hypot(130 * 0.96, 36), tol: 0.02 }, { kpi: 'W50', value: 3534, tol: 0.003 }] },
    { name: 'Corroded tube section: area ratio and stress increase', source: 'Geometry: D = 1500 mm, wall 45 → 39 mm: A = π·(D − t)·t → 178 997 mm² / 205 696 mm² = 87.02 %',
      inputs: {}, expect: [{ kpi: 'resA', value: 100 * ((1.5 - 0.039) * 0.039) / ((1.5 - 0.045) * 0.045), tol: 1e-9 }] },
  ],
};
