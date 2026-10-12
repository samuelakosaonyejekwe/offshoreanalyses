// Structural integrity and life extension of an existing topsides: corrosion wastage applied to the frame and
// re-check of strength over time, consumed and remaining fatigue life, engineering critical assessment and growth
// of a recorded crack, and inspection updating of the failure probability with a risk-based inspection plan.
import { corrosionWastage } from '../../engine/fatigue/degradation.js';
import { snCurve, snConstants, SN_CURVE_NAMES } from '../../engine/fatigue/sn.js';
import { minerWeibull, fatigueLife } from '../../engine/fatigue/damage.js';
import { parisLawBs7910, surfaceCrackGrowth, assessSurfaceCrack, criticalCrackSize } from '../../engine/fatigue/fracture.js';
import { podPreset, crackPopulation, updateAfterInspection, planInspections, populationFailureProbability } from '../../engine/fatigue/inspection.js';
import { gamma } from '../../engine/math/special.js';
import { deckInputs, solverInput, buildDeck, loadCase, solveDeck, memberForces, memberChecks, groupSummary, framePlot, nodeValuesFromMembers, equipmentOf, need, run, fmt, round, adv, check, fsCheck, steelCo2Text, CONFIRM, DECK_LIMITS } from './_common.js';

const CAP = 1e4;   // years reported when a life is effectively unlimited

export default {
  id: 's02-life-extension', suite: 2, title: 'Structural integrity and life extension',
  summary: 'Re-assessment of an ageing topsides: corroded-section strength now and at the end of the extension, consumed and remaining fatigue life, fitness-for-service of a recorded crack, updated failure probability after inspection and the inspection plan and strengthening needed.',
  covers: [14], fidelity: 'numerical',
  standards: ['ISO 19901-9:2019 (structural integrity management)', 'NORSOK N-006 (2015) §6–§9 (assessment of existing structures)', 'DNV-RP-C203 (2021) §2.4, §5 (S–N fatigue)', 'BS 7910:2019 clauses 7–8 (fracture assessment, fatigue crack growth)', 'DNV-RP-C210 (2019) (probabilistic inspection planning)', 'API RP 2SIM (2014)'],
  models: ['Corrosion wastage model (linear after coating breakdown)', 'Matrix stiffness method with reduced sections', 'Residual strength assessment', 'S–N fatigue with Weibull long-term stress ranges (Miner sum)', 'Fracture mechanics: Paris law, Newman–Raju stress intensity, failure assessment diagram', 'Bayesian inspection updating with probability-of-detection curves', 'Monte Carlo crack population'],
  limits: `${DECK_LIMITS} Corrosion is uniform thinning of every member by the same wastage (pitting and local section loss need the measured thickness map); fatigue and fracture are assessed for one governing detail with a Weibull stress-range distribution, a semi-elliptical surface crack in a plate and a constant geometry factor in the probabilistic model. Nonlinear load redistribution after member failure is not modelled.`,
  publishes: ['life.fatigueYears', 'topsides.maxUc'],
  searchTerms: 'life extension structural integrity management ageing offshore topsides corrosion fatigue crack inspection planning',
  inputs: [
    ...deckInputs(),
    { key: 'age', label: 'Age of the structure', unit: 'yr', type: 'number', default: 25, min: 0, max: 80, group: 'Service history', role: 'ic' },
    { key: 'ext', label: 'Life extension sought', unit: 'yr', type: 'number', default: 15, min: 1, max: 50, group: 'Service history', role: 'operation', design: true },
    { key: 'coatLife', label: 'Effective coating life', unit: 'yr', type: 'number', default: 15, min: 0, max: 60, group: 'Corrosion', role: 'ic', help: 'Years before the coating broke down and steel started to waste; maintenance painting restarts it.' },
    { key: 'cRate', label: 'Corrosion rate after coating breakdown', unit: 'mm/yr', type: 'number', default: 0.1, min: 0, max: 1, group: 'Corrosion', role: 'environment', uncertain: { dist: 'lognormal', cov: 0.4 }, help: 'Per exposed surface, from thickness surveys. Marine atmosphere 0.05–0.15 mm/yr for bare carbon steel; calibrate against measured wall loss.' },
    { key: 'Hx', label: 'Lateral environmental load along x', unit: 'kN', type: 'number', default: 900, min: -5e4, max: 5e4, group: 'Loads for the strength re-check', role: 'bc', bind: 'topsides.windFx' },
    { key: 'gG', label: 'Partial factor on permanent and variable actions', unit: '-', type: 'number', default: 1.1, min: 0.8, max: 1.6, group: 'Loads for the strength re-check', role: 'criteria', help: `Assessment of existing structures may use reduced factors where loads are measured (NORSOK N-006). ${CONFIRM}` },
    { key: 'gE', label: 'Partial factor on environmental actions', unit: '-', type: 'number', default: 1.35, min: 0, max: 1.8, group: 'Loads for the strength re-check', role: 'criteria', help: CONFIRM },
    { key: 'snName', label: 'S–N curve of the governing detail', type: 'select', default: 'F', group: 'Fatigue of the governing detail', role: 'criteria', options: SN_CURVE_NAMES.map((v) => ({ value: v, label: `DNV-RP-C203 ${v}` })) },
    { key: 'snEnv', label: 'Environment', type: 'select', default: 'air', group: 'Fatigue of the governing detail', role: 'environment', options: [{ value: 'air', label: 'Air (coated)' }, { value: 'seawater-cp', label: 'Seawater with cathodic protection' }, { value: 'free-corrosion', label: 'Free corrosion' }] },
    { key: 'dS0', label: 'Largest nominal stress range in a year', unit: 'MPa', type: 'number', default: 80, min: 1, max: 800, group: 'Fatigue of the governing detail', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Stress range exceeded once per year at the detail (wave-induced support motion, crane and machinery cycles), from the cyclic loading history or strain measurements.' },
    { key: 'hW', label: 'Weibull shape of the stress ranges', unit: '-', type: 'number', default: 1, min: 0.4, max: 2, group: 'Fatigue of the governing detail', role: 'bc', help: '0.8–1.1 for wave-dominated loading; higher for narrow-banded machinery cycles.' },
    { key: 'nu', label: 'Stress cycles per year', unit: '-', type: 'number', default: 5e6, min: 1e3, max: 1e9, group: 'Fatigue of the governing detail', role: 'bc' },
    { key: 'scf', label: 'Stress concentration factor', unit: '-', type: 'number', default: 1.3, min: 1, max: 10, group: 'Fatigue of the governing detail', role: 'geometry' },
    { key: 'tDet', label: 'Plate thickness at the detail', unit: 'mm', type: 'number', default: 25, min: 5, max: 150, group: 'Fatigue of the governing detail', role: 'geometry' },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 2, min: 1, max: 10, group: 'Fatigue of the governing detail', role: 'criteria', help: CONFIRM },
    { key: 'aCrack', label: 'Recorded crack depth', unit: 'mm', type: 'number', default: 1.5, min: 0.05, max: 100, group: 'Recorded crack', role: 'ic', help: 'Depth of the deepest surface-breaking indication at the detail from the last inspection (use the sizing value plus tolerance).' },
    { key: 'cCrack', label: 'Recorded crack length (2c)', unit: 'mm', type: 'number', default: 15, min: 0.2, max: 2000, group: 'Recorded crack', role: 'ic' },
    { key: 'sigMax', label: 'Peak membrane stress at the detail', unit: 'MPa', type: 'number', default: 180, min: 1, max: 700, group: 'Recorded crack', role: 'bc', help: 'Extreme tensile stress across the crack (primary), including the stress concentration.' },
    { key: 'resFrac', label: 'Residual stress as a fraction of yield', unit: '-', type: 'number', default: 0.5, min: 0, max: 1, group: 'Recorded crack', role: 'ic', help: 'BS 7910: 1.0 for as-welded joints transverse to the weld, lower after post-weld heat treatment or proof loading.' },
    { key: 'Kmat', label: 'Fracture toughness K_mat', unit: 'MPa√m', type: 'number', default: 120, min: 20, max: 400, group: 'Recorded crack', role: 'material', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'From CTOD or Charpy correlation at the minimum design temperature; older steels can be well below modern grades.' },
    { key: 'inspect', label: 'Last inspection of the detail', type: 'select', default: 'none', group: 'Inspection', role: 'ic', options: [{ value: 'none', label: 'Not inspected' }, { value: 'clear', label: 'Inspected at the present age, nothing found' }, { value: 'found', label: 'Inspected at the present age, crack sized as recorded' }] },
    { key: 'method', label: 'Inspection method', type: 'select', default: 'ec-normal', group: 'Inspection', role: 'operation', options: [{ value: 'ec-good', label: 'Eddy current / MPI / ACFM, good conditions' }, { value: 'ec-normal', label: 'Eddy current / MPI / ACFM, normal conditions' }, { value: 'ut', label: 'Ultrasonic testing' }, { value: 'visual-easy', label: 'Close visual, easy access' }] },
    { key: 'pfTarget', label: 'Target annual failure probability', unit: '-', type: 'number', default: 1e-3, min: 1e-5, max: 0.1, group: 'Inspection', role: 'criteria', help: `1e-3 for details with redundancy and low consequence, 1e-4 for high consequence. ${CONFIRM}` },
    { key: 'nMC', label: 'Samples in the crack population', type: 'int', default: 4000, min: 500, max: 100000, group: 'Solver', role: 'control', help: 'Monte Carlo samples (fixed seed). The smallest probability resolved is about 10 / samples.' },
    solverInput(1),
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'ucEnd' },
  compute(inp, ctx) {
    // SI: mm → m, MPa → Pa, MPa√m → Pa√m.
    const wast = (t) => (inp.cRate > 0 ? corrosionWastage({ time: t, model: 'linear', coatingLife: inp.coatLife, rate: inp.cRate / 1e3 }).depth * 1e3 : 0);   // mm per surface
    // Corroded inputs: tubes lose wall from the outside, I-section webs and flanges and the plate from both faces.
    const thin = (d) => { const o = { ...inp }; for (const k of ['gtw', 'gtf', 'btw', 'btf']) o[k] = inp[k] - 2 * d; for (const k of ['ct', 'brt']) if (inp[k] != null) o[k] = inp[k] - d; if (inp.tPlate > 0) o.tPlate = Math.max(inp.tPlate - 2 * d, 0.1); for (const k of ['gtw', 'gtf', 'btw', 'btf', 'ct']) need(o[k] > 1, `After ${fmt(d)} mm of wastage per surface the ${k} thickness is consumed: the structure cannot be assessed as a frame at that age.`); return o; };
    const tEnd = inp.age + inp.ext, times = [0, inp.age / 2, inp.age, inp.age + inp.ext / 2, tEnd];
    const states = times.map((t) => {
      const d = wast(t), ci = thin(d), deck = buildDeck(ci, ctx), eqs = equipmentOf(deck, ci), zs = deck.levelsZ, zsum = zs.reduce((s, z) => s + z, 0), extra = [];
      deck.levelNodes.forEach((ns, k) => { for (const n of ns) extra.push({ node: n, F: [(inp.Hx * 1e3 * inp.gE * zs[k]) / zsum / ns.length, 0, 0, 0, 0, 0] }); });
      const { ref, sol } = solveDeck(deck, [loadCase(deck, ci, { gDead: inp.gG, gEquip: inp.gG, gLive: inp.gG, extra, eqs })], inp.nSub);
      const rows = memberChecks(deck, memberForces(deck, ref, sol.results.design));
      return { t, d, deck, ref, r: sol.results.design, rows, gov: rows.reduce((a, c) => (c.uc > a.uc ? c : a)), groups: groupSummary(rows) };
    });
    const [s0, , sNow, , sEnd] = states;
    // Fatigue: Weibull scale from the yearly maximum; stress grows as the section thins (ratio of governing-group area).
    const curve = run(() => snCurve(inp.snName, inp.snEnv), 'S–N curve'), sn = snConstants(curve, { thickness: inp.tDet / 1e3, scf: inp.scf });
    const dRate = (fac) => minerWeibull({ sn, shape: inp.hW, cycles: inp.nu, maxStressRange: inp.dS0 * 1e6 * fac }).damage;
    const tRef = Math.min(inp.gtf, inp.btf, inp.ct), stressFac = (d) => tRef / Math.max(tRef - 2 * d, 0.2 * tRef);
    // Consumed damage: yearly rates integrated with the stress factor of each year.
    let D0 = 0; for (let y = 0; y < Math.ceil(inp.age); y++) D0 += dRate(stressFac(wast(y + 0.5))) * Math.min(1, inp.age - y);
    const rateNow = dRate(stressFac(sNow.d)), rateEnd = dRate(stressFac(sEnd.d)), rateFut = 0.5 * (rateNow + rateEnd);
    const fl = fatigueLife({ damageRate: rateFut, designLife: inp.ext, dff: inp.dff, accumulatedDamage: D0 });
    const remFat = Math.max(Math.min(fl.remainingLife, CAP), 0), ucFat = Math.min(fl.utilisation, 99);
    // Crack: fitness for service now, critical size, and growth under the equivalent constant-amplitude range.
    const st = s0.deck.steel, t = inp.tDet / 1e3, a0 = inp.aCrack / 1e3, c0 = inp.cCrack / 2e3;
    need(a0 < 0.9 * t, 'The recorded crack depth must be less than 90 % of the plate thickness: a near-through crack needs a through-thickness assessment and immediate repair.');
    need(c0 >= a0 / 2, 'The recorded crack length 2c must be at least the crack depth (a/c ≤ 2 for the surface-crack solution).');
    const ecaArgs = { t, Pm: inp.sigMax * 1e6, Qm: inp.resFrac * st.fy, yieldStrength: st.fy, tensileStrength: st.fu, E: st.E, toughness: inp.Kmat * 1e6 };
    const eca = run(() => assessSurfaceCrack({ ...ecaArgs, a: a0, c: c0 }), 'Fracture assessment');
    const crit = run(() => criticalCrackSize({ ...ecaArgs, aspectRatio: Math.min(a0 / c0, 2) }), 'Critical crack size');
    const aCrit = crit.limitedBy === 'wall' ? t : Math.max(crit.a, 0);
    // The growth threshold is set to zero: with a variable-amplitude spectrum represented by one equivalent range the
    // large cycles keep the crack growing even when the equivalent ΔK is below the constant-amplitude threshold.
    const law = { ...parisLawBs7910({ environment: inp.snEnv === 'air' ? 'air' : 'marine', curve: 'simplified' }), threshold: 0 }, m = law.stages[0].m, A = law.stages[0].A;
    const q = (inp.dS0 * 1e6) / Math.log(inp.nu) ** (1 / inp.hW), Seq = q * gamma(1 + m / inp.hW) ** (1 / m) * inp.scf * stressFac(sNow.d);
    let crackYears = 0, grow = { history: { cycles: [0], a: [a0], c: [c0] }, status: 'unstable now' };
    if (aCrit > a0 * 1.001) {
      grow = run(() => surfaceCrackGrowth({ a0, c0, t, membraneRange: Seq, law, maxCycles: CAP * inp.nu, throughFraction: Math.min(aCrit / t, 1), stepFraction: 0.02 }), 'Crack growth');
      crackYears = Math.min(Number.isFinite(grow.cycles) ? grow.cycles / inp.nu : CAP, CAP);
    }
    // Inspection updating: prior crack population (DNV-RP-C210 initial sizes), outcome of the last inspection, plan.
    const podPar = podPreset(inp.method);
    const prior = run(() => crackPopulation({ a0: { dist: 'exponential', mean: 0.11e-3 }, C: { dist: 'lognormal', mean: A, cov: 0.5 }, m, stressRange: Seq, stressFactor: { dist: 'lognormal', mean: 1, cov: 0.2 }, Y: 1.12, cyclesPerYear: inp.nu, aCrit: Math.max(aCrit, 1.05 * a0), n: Math.round(inp.nMC), seed: 11 }), 'Crack population');
    let pop = prior;
    if (inp.inspect !== 'none' && inp.age > 0) pop = run(() => updateAfterInspection(prior, inp.inspect === 'clear' ? { time: inp.age, pod: podPar, detected: false } : { time: inp.age, pod: podPar, detected: true, measuredSize: a0, sizingStd: Math.max(0.3 * a0, 0.3e-3) }), 'Inspection updating');
    const yrs = []; for (let y = Math.floor(inp.age) + 1; y <= Math.ceil(tEnd); y++) yrs.push(y);
    const pfPrior = populationFailureProbability(prior, yrs, { survivedTo: inp.age }), pfPost = populationFailureProbability(pop, yrs, { survivedTo: inp.age });
    const plan = run(() => planInspections({ population: pop, pod: podPar, targetAnnualPf: inp.pfTarget, serviceLife: tEnd, maxInspections: 30 }), 'Inspection planning');
    const future = plan.inspectionTimes.filter((y) => y > inp.age), fin = (v) => (Number.isFinite(v) ? v : 0);
    const pfEnd = fin(pfPost.annualPf[yrs.length - 1]), pfMaxPlan = Math.max(0, ...plan.table.filter((r) => r.year > inp.age).map((r) => fin(r.annualPf)));
    const checks = [
      check(`Member strength at the end of the extension (${sEnd.gov.group.toLowerCase()} ${sEnd.gov.i + 1})`, sEnd.gov.uc, 1, 'interaction', sEnd.gov.ref, `Governing utilisation with every member thinned by ${fmt(sEnd.d)} mm per exposed surface at year ${fmt(tEnd)}. Above 1.0 the member needs strengthening, load reduction or coating repair before that age.`),
      check('Fatigue damage at the end of the extension', ucFat, 1, 'Miner sum × DFF', `DNV-RP-C203 §2.4, curve ${inp.snName} (${inp.snEnv}), DFF ${inp.dff}; NORSOK N-006 §8`, 'Consumed plus future Miner damage of the governing detail times the design fatigue factor. Above 1.0 the calculated fatigue life is used up within the extension and continued service must rest on inspection or improvement of the detail.'),
      check('Recorded crack: fracture assessment (FAD)', 1 / Math.max(eca.margin, 1e-6), 1, 'Kr / allowable Kr', 'BS 7910:2019 clause 7, Option 1 failure assessment diagram', 'Position of the recorded crack on the failure assessment diagram under the peak stress with welding residual stress. Above 1.0 the crack is predicted to be unstable: repair before the next storm.'),
      fsCheck('Recorded crack: years of growth to the critical size', Math.max(crackYears, 1e-3), inp.ext, 'BS 7910:2019 clause 8 (Paris law, Newman–Raju)', `Time for the recorded crack to grow from ${inp.aCrack} mm to the critical depth of ${fmt(aCrit * 1e3)} mm under the equivalent stress range. It must exceed the extension, or the inspection interval if the crack is monitored.`),
      check('Annual failure probability with the inspection plan', Math.max(pfMaxPlan, 1e-9), inp.pfTarget, '1/yr', 'DNV-RP-C210 (target annual probability of fatigue failure); ISO 19901-9 §9', 'Largest annual probability of fatigue failure of the detail during the extension when the planned inspections are carried out and find nothing. Above the target more frequent or more sensitive inspection is needed.'),
    ];
    const gov = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = gov.uc <= 1;
    // Strengthening: area the failing groups need to return to 0.9 utilisation at the end of the extension.
    const addT = sEnd.groups.filter((g) => g.uc > 1).reduce((s, g) => s + (g.uc / 0.9 - 1) * sEnd.deck.orig.filter((o) => o.group === g.group).reduce((x, o) => x + o.sec.A * o.L * st.rho, 0) / 1e3, 0);
    const lost = Math.max(s0.deck.frameSteelT - sEnd.deck.frameSteelT, 0);
    const hist = grow.history, step = Math.max(1, Math.floor(hist.a.length / 200)), gi = hist.a.map((_, i) => i).filter((i) => i % step === 0 || i === hist.a.length - 1);
    return {
      kpis: [
        { key: 'wallLoss', label: 'Wastage per surface today', unit: 'mm', value: sNow.d, explain: 'Thickness lost from each exposed steel surface since the coating broke down. Compare with thickness survey readings and adjust the rate.' },
        { key: 'wallLossEnd', label: 'Wastage per surface at the end of the extension', unit: 'mm', value: sEnd.d, explain: 'Projected loss if no coating maintenance is done.' },
        { key: 'ucNew', label: 'Member utilisation as built', unit: '-', value: round(s0.gov.uc, 5), explain: 'Governing utilisation with the original thicknesses, for reference.' },
        { key: 'ucNow', label: 'Member utilisation today', unit: '-', value: round(sNow.gov.uc, 5), limit: 1, sense: 'max', explain: 'Governing utilisation with today’s corroded sections.' },
        { key: 'ucEnd', label: 'Member utilisation at the end of the extension', unit: '-', value: round(sEnd.gov.uc, 5), limit: 1, sense: 'max', explain: 'Governing utilisation with the sections projected to the end of the extension.' },
        { key: 'damage', label: 'Fatigue damage consumed', unit: '-', value: D0, explain: 'Characteristic Miner sum accumulated to date at the governing detail (without the design fatigue factor). 1.0 is the calculated end of life.' },
        { key: 'damageRate', label: 'Fatigue damage rate today', unit: '1/yr', value: rateNow, explain: 'Miner damage accumulated per year at today’s section thickness.' },
        { key: 'fatigueRemaining', label: 'Remaining fatigue life', unit: 'yr', value: remFat, limit: inp.ext, sense: 'min', explain: `Years from today until the consumed plus future damage times the design fatigue factor reaches 1.0 (reported as ${CAP} when effectively unlimited).` },
        { key: 'aCrit', label: 'Critical crack depth', unit: 'mm', value: aCrit * 1e3, explain: 'Depth at which a crack of the recorded shape becomes unstable under the peak stress (equal to the thickness when the wall is breached first).' },
        { key: 'crackYears', label: 'Crack growth life', unit: 'yr', value: crackYears, limit: inp.ext, sense: 'min', explain: 'Years for the recorded crack to reach the critical depth; sets the longest inspection interval if the crack is left in place.' },
        { key: 'pfEnd', label: 'Annual failure probability in the last year (no further inspection)', unit: '1/yr', value: pfEnd, explain: 'Probability that the detail fails by fatigue in the final year of the extension, given survival to today and the last inspection result.' },
        { key: 'nInsp', label: 'Inspections needed during the extension', unit: '-', value: future.length, explain: 'Number of inspections of the detail that keep the annual failure probability below the target, assuming each finds nothing.' },
        { key: 'strengthen', label: 'Strengthening steel', unit: 't', value: addT, explain: 'Steel needed to bring the member groups that fail at the end of the extension back to a utilisation of 0.9.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Utilisation and wastage over the life', x: { label: 'Age', unit: 'yr' }, y: { label: 'Utilisation (-) / wastage (mm)', unit: '-' }, series: [{ name: 'Governing member utilisation', x: times.map((v) => round(v)), y: states.map((s) => round(s.gov.uc)) }, { name: 'Wastage per surface (mm)', x: times.map((v) => round(v)), y: states.map((s) => round(s.d)) }, { name: 'Limit', x: [0, tEnd], y: [1, 1], dash: true }, { name: 'Today', x: [inp.age, inp.age], y: [0, 1], dash: true }],
          explain: 'Governing member utilisation as the sections thin, with the wastage itself. The age at which the utilisation crosses 1.0 is the latest date for coating repair or strengthening.' },
        framePlot({ title: 'Utilisation at the end of the extension', ref: sEnd.ref, u: sEnd.r.displacements, values: nodeValuesFromMembers(sEnd.ref, sEnd.rows.map((x) => x.uc)), valueLabel: 'Utilisation', valueUnit: '-', explain: 'Frame coloured by utilisation with the projected corroded sections: the members to strengthen or re-coat first are the ones nearest 1.0.' }),
        { type: 'line', title: 'Growth of the recorded crack', x: { label: 'Years from today', unit: 'yr' }, y: { label: 'Crack size', unit: 'mm' }, series: [{ name: 'Depth a', x: gi.map((i) => round(Math.min(hist.cycles[i] / inp.nu, CAP))), y: gi.map((i) => round(hist.a[i] * 1e3)) }, { name: 'Half-length c', x: gi.map((i) => round(Math.min(hist.cycles[i] / inp.nu, CAP))), y: gi.map((i) => round(hist.c[i] * 1e3)) }, { name: 'Critical depth', x: [0, Math.max(round(crackYears), 1)], y: [round(aCrit * 1e3), round(aCrit * 1e3)], dash: true }],
          explain: 'Depth and half-length of the recorded crack against time under the equivalent stress range. Growth accelerates sharply near the end: most of the life is spent while the crack is small, which is why early detection matters.' },
        { type: 'line', title: 'Annual failure probability of the detail', x: { label: 'Age', unit: 'yr' }, y: { label: 'Annual failure probability', unit: '1/yr' }, series: [{ name: 'Before the last inspection result', x: yrs, y: Array.from(pfPrior.annualPf, (v) => fin(v)) }, { name: 'After the last inspection result', x: yrs, y: Array.from(pfPost.annualPf, (v) => fin(v)) }, { name: 'With the inspection plan', x: plan.table.filter((r) => r.year > inp.age).map((r) => r.year), y: plan.table.filter((r) => r.year > inp.age).map((r) => fin(r.annualPf)) }, { name: 'Target', x: [yrs[0], yrs[yrs.length - 1]], y: [inp.pfTarget, inp.pfTarget], dash: true }],
          explain: 'Probability of fatigue failure per year through the extension: prior, updated for the last inspection outcome, and with the planned inspections (each drop is an inspection that finds nothing). The planned curve should stay below the target line.' },
      ],
      tables: [
        { title: 'Utilisation by member group', columns: [{ key: 'group', label: 'Group' }, { key: 'u0', label: 'As built' }, { key: 'u1', label: 'Today' }, { key: 'u2', label: 'End of extension' }, { key: 'ref', label: 'Governing check' }],
          rows: sEnd.groups.map((g) => ({ group: g.group, u0: round(s0.groups.find((x) => x.group === g.group)?.uc ?? 0), u1: round(sNow.groups.find((x) => x.group === g.group)?.uc ?? 0), u2: round(g.uc), ref: g.ref })), explain: 'Governing utilisation of each member group at three ages. Thin-walled groups lose the largest share of their section and rise fastest.' },
        { title: 'Inspection plan for the governing detail', columns: [{ key: 'n', label: 'Inspection' }, { key: 'year', label: 'At age', unit: 'yr' }, { key: 'from', label: 'Years from today', unit: 'yr' }], rows: future.length ? future.map((y, i) => ({ n: i + 1, year: round(y), from: round(y - inp.age) })) : [{ n: '-', year: '-', from: 'none needed' }],
          explain: `Inspection times with the selected method (${inp.method}) that keep the annual failure probability below the target, assuming no crack is found. A finding resets the plan.` },
      ],
      notes: [...(s0.deck.imported ? ['Imported members keep their as-modelled sections: only the default tube and plate take the computed wastage, so enter measured thicknesses in the imported model.'] : []), `${gov.name} governs at ${fmt(gov.uc)}.`, `Wastage ${fmt(sNow.d)} mm today and ${fmt(sEnd.d)} mm at year ${fmt(tEnd)}; frame steel lost by then ${fmt(lost)} t.`,
        `Fatigue: consumed damage ${fmt(D0)}, remaining life ${fmt(remFat)} yr at DFF ${inp.dff}. Crack: FAD margin ${fmt(eca.margin)} (${eca.governingPoint} point), critical depth ${fmt(aCrit * 1e3)} mm (${crit.limitedBy}), growth ${grow.status}.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.85 ? 'warn' : 'ok', `The structure supports a ${inp.ext}-year extension: utilisation ${fmt(sEnd.gov.uc)} at year ${fmt(tEnd)}, remaining fatigue life ${fmt(remFat)} yr, crack growth life ${fmt(crackYears)} yr, ${future.length} inspection(s) of the governing detail.`,
          `Confirm the corrosion rate with a thickness survey of the ${sEnd.gov.group.toLowerCase()} group, inspect the governing detail ${future.length ? `at age ${future.map((y) => fmt(y, 3)).join(', ')}` : 'at the regular interval'} with ${inp.method}, and re-run after each survey.`, 'ISO 19901-9:2019 §8–§10; NORSOK N-006 (2015) §6–§9; DNV-RP-C210 (2019)',
          `Keeping the existing ${fmt(s0.deck.steelT)} t of steel in service for ${inp.ext} more years avoids about ${fmt(s0.deck.steelT * 1.9)} t CO₂e of replacement steel (screening factor 1.9 t CO₂e/t); coating maintenance that stops the projected ${fmt(lost)} t of wastage is the lowest-carbon measure.`)
          : adv('engineering', 'fail', `${gov.name} is not met (utilisation ${fmt(gov.uc)}): member utilisation ${fmt(sEnd.gov.uc)} at year ${fmt(tEnd)}, remaining fatigue life ${fmt(remFat)} yr, crack growth life ${fmt(crackYears)} yr.`,
            gov.name.startsWith('Member') ? `Re-coat before the utilisation reaches 1.0 and strengthen the failing groups with about ${fmt(addT)} t of steel (doublers, added flange plates or a new brace), or reduce deck loads.`
              : gov.name.startsWith('Fatigue') ? `Improve the detail (toe grinding or hammer peening gives a factor of about 2–4 on life), reduce the stress concentration from ${inp.scf}, or justify continued service by inspection at intervals no longer than ${fmt(Math.max(crackYears / 2, 1), 2)} yr.`
                : gov.name.includes('FAD') ? 'Repair the crack now: grind out and re-weld, or bolt a clamp or doubler across it, and lower the load on the member until it is done.'
                  : gov.name.includes('growth') ? `Grind out or repair the crack, or monitor it at intervals no longer than ${fmt(Math.max(crackYears / 2, 0.5), 2)} yr with crack-arrest holes or a bolted repair ready.` : `Inspect more often or with a more sensitive method than ${inp.method}; the plan needs ${future.length} inspection(s) and still exceeds the target.`,
            'ISO 19901-9:2019 §9–§11; NORSOK N-006 (2015) §8; BS 7910:2019 clauses 7–8', `Strengthening needs about ${steelCo2Text(Math.max(addT, 0.5))} Local repair and improved coating extend the life of ${fmt(s0.deck.steelT)} t of existing steel for a small fraction of the emissions of replacement.`),
        adv('economic', 'info', `The extension rests on ${future.length} planned inspection(s) of the governing detail and on stopping ${fmt(sEnd.d - sNow.d)} mm of further wastage; strengthening identified: ${fmt(addT)} t.`,
          'Compare the cost of the inspection and coating campaign with the value of the deferred replacement or decommissioning; schedule inspections with planned shutdowns and use the updated failure probability to drop inspections that add no information.', 'Risk-based inspection per DNV-RP-G101 / API RP 580; value of information; whole-life costing per ISO 15686-5'),
      ],
      publish: { 'life.fatigueYears': remFat, 'topsides.maxUc': sEnd.gov.uc },
      quantities: { steel_t: addT },
    };
  },
  verify: [
    { // Linear wastage after coating breakdown: d = r (t − Tc) = 0.1·(25 − 15) = 1.0 mm today, 0.1·(40 − 15) = 2.5 mm at year 40.
      name: 'Linear corrosion wastage after coating breakdown: d = r (t − Tc)', source: 'Hand calculation (constant-rate model; Paik & Thayamballi 2003 linear form)',
      inputs: { age: 25, ext: 15, coatLife: 15, cRate: 0.1 }, expect: [{ kpi: 'wallLoss', value: 1.0, tol: 1e-6 }, { kpi: 'wallLossEnd', value: 2.5, tol: 1e-6 }] },
    { // Single-slope curve (D, free corrosion: log a = 11.687, m = 3), Weibull h = 1, no corrosion, t = tRef, SCF = 1:
      // q = Δσ0/ln(n0) = 154.249/ln(5e6) = 10.000 MPa; D per year = n0 q³ Γ(1 + 3)/a = 5e6·1000·6/10^11.687 = 0.06168;
      // 25 years → 1.542.
      name: 'Weibull closed-form Miner damage on a single-slope S–N curve: D = n0 q^m Γ(1 + m/h)/a', source: 'DNV-RP-C203 (2021) Eq. (5.1.1) with Table 2-4 curve D (hand calculation)',
      inputs: { snName: 'D', snEnv: 'free-corrosion', dS0: 154.249, hW: 1, nu: 5e6, scf: 1, tDet: 25, cRate: 0, age: 25 }, expect: [{ kpi: 'damageRate', value: 0.06168, tol: 0.005 }, { kpi: 'damage', value: 1.542, tol: 0.005 }] },
  ],
};
