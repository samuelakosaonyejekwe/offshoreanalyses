// Pipeline integrity, remaining life and reliability: remaining strength of corrosion defects from inspection data
// (ASME B31G, Modified B31G, effective area, DNV-RP-F101), growth to the next inspection and remaining life, dent
// screening, extreme-value wall loss and the time-dependent failure probability of the governing defect by FORM.
import { b31gOriginal, b31gModified, effectiveArea, dnvF101Single, remainingLife, dentStrain, dentDepthCheck, wallLossExtreme, fitnessForService } from '../../engine/pipeline/integrity.js';
import { form } from '../../engine/reliability/form.js';
import { pipeInputs, pipeFrom, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, sev, worst, linspace,
  steelPerMetre, vesselCo2, CONFIRM, RBI_BASIS, CO2_STEEL, BAR, MPA } from './_common.js';

const G_D = 'Inspection data';
const G_G = 'Defect growth';
const G_R = 'Reliability model';
const G_C = 'Acceptance criteria';
const DEFECTS = [
  { kp: 1.2, depth: 22, length: 120 }, { kp: 3.85, depth: 38, length: 260 }, { kp: 4.1, depth: 15, length: 60 },
  { kp: 6.7, depth: 45, length: 90 }, { kp: 8.3, depth: 30, length: 400 },
];
const TARGET = { low: 1e-3, medium: 1e-4, high: 1e-5 };

export default {
  id: 's08-integrity',
  suite: 8,
  title: 'Pipeline integrity, remaining life and reliability',
  summary: 'Assesses corrosion defects and dents from inspection, ranks the critical segments, and gives the remaining life, re-inspection date and failure probability of the governing defect.',
  covers: [18],
  fidelity: 'analytical',
  standards: ['ASME B31G-2012 (R2017) Level 1 and Level 2 (effective area)', 'DNV-RP-F101 (2019) Part B (allowable stress) for single defects', 'ASME B31.8 (2022) §851.4 and Appendix R (dents)', 'DNV-RP-G101 / API RP 580 (risk-based inspection)', 'DNV-ST-F101 (2021) Sec. 2 (target failure probabilities)'],
  models: ['Corrosion defect remaining-strength models (B31G, Modified B31G, RSTRENG effective area, DNV-RP-F101)', 'Linear defect growth and remaining life', 'Dent strain from curvatures', 'Gumbel extreme-value extrapolation of wall loss', 'First-order reliability method (FORM) on the burst limit state'],
  limits: 'Single, isolated, blunt corrosion defects under internal pressure only: interacting defects, axial compression or bending on the defect, cracks, gouges and dent–gouge combinations need the specific methods of the codes. DNV-RP-F101 Part A partial factors are not tabulated. The reliability model uses assumed distributions (typed below) for sizing error, growth rate, strength, pressure and model uncertainty — replace them with tool-specific and site data. Target failure probabilities are recalled values to be confirmed.',
  publishes: ['pipeline.remainingLife', 'pipeline.safePressure'],
  searchTerms: 'pipeline corrosion defect assessment B31G DNV-RP-F101 remaining life reliability inspection',
  inputs: [
    ...pipeInputs({ omit: ['tCoat', 'rhoCoat', 'tConc', 'concrete', 'rhoCont', 'wsGiven', 'rhoW'] }),
    num('maop', 'Maximum allowable operating pressure', 'bar', 150, 1, 1400, 'Pipe', 'bc', 'Pressure the defects must sustain (differential pressure at the defect).', { uncertain: { dist: 'normal', cov: 0.03 } }),
    { key: 'defects', label: 'Corrosion defects from inspection', type: 'table', group: G_D, role: 'ic',
      columns: [{ key: 'kp', label: 'Chainage', unit: 'km', type: 'number' }, { key: 'depth', label: 'Depth', unit: '% wall', type: 'number' }, { key: 'length', label: 'Axial length', unit: 'mm', type: 'number' }],
      default: DEFECTS, help: 'Metal-loss features reported by in-line inspection: peak depth as a percentage of the wall and axial length.' },
    { key: 'profile', label: 'River-bottom profile of the deepest defect (optional)', type: 'series', default: { x: [], y: [] }, xLabel: 'Axial position', xUnit: 'mm', yLabel: 'Depth', yUnit: 'mm', group: G_D, role: 'ic',
      help: 'At least three depth readings along the defect enable the Level 2 effective-area assessment.' },
    { key: 'blockMax', label: 'Deepest wall loss per inspected block (optional)', type: 'series', default: { x: [], y: [] }, xLabel: 'Block', xUnit: '-', yLabel: 'Maximum wall loss', yUnit: 'mm', group: G_D, role: 'ic',
      help: 'At least five block maxima from a partial inspection give a Gumbel estimate of the deepest pit over the whole line.' },
    num('blocks', 'Total blocks in the line / blocks inspected', '-', 20, 1.1, 1e5, G_D, 'operation', 'Return period for the extreme-value extrapolation: total length over inspected block length.'),
    num('dentDepth', 'Dent depth (0 = none)', 'mm', 6, 0, 200, G_D, 'ic', 'Deepest plain dent reported.'),
    num('dentLen', 'Dent length', 'mm', 250, 10, 3000, G_D, 'ic', 'Axial length of the dent.'),
    num('dentR1', 'Dent transverse radius of curvature (0 = not measured)', 'mm', 400, -5000, 5000, G_D, 'ic', 'Negative for a re-entrant dent. With the longitudinal radius it gives the dent strain.'),
    num('dentR2', 'Dent longitudinal radius of curvature', 'mm', -3000, -20000, 20000, G_D, 'ic', 'Negative for a re-entrant profile.'),
    { key: 'dentOnWeld', label: 'Dent affects a weld', type: 'bool', default: false, group: G_D, role: 'ic', help: 'Dents on girth or seam welds are limited to 2 % of the diameter.' },
    sel('method', 'Remaining-strength method', 'modified', [['modified', 'Modified B31G (0.85 dL)'], ['original', 'ASME B31G (original)'], ['dnv', 'DNV-RP-F101 Part B']], G_C, 'control', 'Method used for acceptance, ranking and remaining life; all are reported for the governing defect.'),
    num('sf', 'Safety factor on failure pressure (B31G)', '-', 1.39, 1.1, 2.5, G_C, 'criteria', `1/0.72 = 1.39 is customary; not below 1.25. ${CONFIRM}`),
    num('F2', 'Usage factor F2 (DNV-RP-F101 Part B)', '-', 0.72, 0.4, 0.9, G_C, 'criteria', `Operational usage factor; total factor F = 0.9·F2. ${CONFIRM}`),
    num('depthRate', 'Depth growth rate', 'mm/yr', 0.2, 0, 3, G_G, 'operation', 'Corrosion growth of the defect depth, from run-to-run comparison or the corrosion model.', { bind: 'pipeline.corrosionRate', uncertain: { dist: 'lognormal', cov: 0.4 } }),
    num('lengthRate', 'Length growth rate', 'mm/yr', 2, 0, 100, G_G, 'operation', 'Axial growth of the defect.'),
    num('interval', 'Planned re-inspection interval', 'years', 5, 0.5, 30, G_G, 'operation', 'Time to the next in-line inspection.'),
    num('life', 'Remaining service life required', 'years', 20, 1, 100, G_G, 'operation', 'Period over which the line is to remain in service.', { bind: 'life.designYears' }),
    sel('safetyClass', 'Safety class for the target reliability', 'medium', [['low', 'Low'], ['medium', 'Medium'], ['high', 'High']], G_R, 'criteria', `Annual target failure probability 1e-3 / 1e-4 / 1e-5 (DNV-ST-F101 Sec. 2, ultimate limit states). ${CONFIRM}`),
    num('sdDepth', 'Depth sizing standard deviation', '% wall', 8, 0.5, 30, G_R, 'criteria', 'In-line inspection tool accuracy (±10 % wall at 80 % confidence ≈ 7.8 % standard deviation).'),
    num('covRate', 'Coefficient of variation of the growth rate', '-', 0.4, 0, 1.5, G_R, 'criteria', 'Uncertainty of the corrosion growth rate.'),
    num('covFu', 'Coefficient of variation of tensile strength', '-', 0.04, 0.01, 0.15, G_R, 'criteria', 'Line-pipe tensile strength scatter; the mean is taken as 1.08 × SMTS.'),
    num('covModel', 'Coefficient of variation of the burst model', '-', 0.1, 0.01, 0.4, G_R, 'criteria', 'Scatter of test burst pressure over the DNV-RP-F101 capacity equation (mean bias 1.05).'),
    num('covP', 'Coefficient of variation of the annual maximum pressure', '-', 0.03, 0.005, 0.2, G_R, 'criteria', 'Annual maximum operating pressure (Gumbel) relative to MAOP.'),
  ],
  compute(inp, ctx) {
    const p = pipeFrom(inp, ctx);
    const D = p.D, t = p.t, maop = inp.maop * BAR;
    const rows = (inp.defects || []).map((r) => ({ kp: Number(r.kp), dPct: Number(r.depth), L: Number(r.length) / 1e3 })).filter((r) => r.dPct > 0 && r.L > 0);
    need(rows.length > 0, 'No corrosion defect is listed: enter at least one feature with a depth and a length.');
    for (const r of rows) need(r.dPct < 100, `The defect at ${fmt(r.kp)} km has a depth of ${fmt(r.dPct)} % of the wall: a through-wall defect cannot be assessed, repair it.`);
    const assess = inp.method === 'original' ? (a) => b31gOriginal({ ...a, safetyFactor: inp.sf }) : inp.method === 'dnv' ? (a) => dnvF101Single({ D: a.D, t: a.t, d: a.d, L: a.L, fu: p.smts, F2: inp.F2 }) : (a) => b31gModified({ ...a, safetyFactor: inp.sf });
    const dRate = inp.depthRate / 1e3, lRate = inp.lengthRate / 1e3;
    const res = rows.map((r) => {
      const d = Math.min(r.dPct / 100, 0.95) * t;
      const def = { D, t, d, L: r.L, smys: p.smys };
      const a = run(() => assess(def), `Defect at ${fmt(r.kp)} km`);
      const life = run(() => remainingLife({ assess, defect: def, maop, depthRate: dRate, lengthRate: lRate, depthLimit: 0.8, horizon: 100, reinspectFraction: 0.5, interval: inp.interval }), 'Remaining life');
      return { ...r, d, a, life, erf: maop / Math.max(a.safePressure, 1), pAtInsp: life.deratedPressure ?? a.safePressure };
    });
    const g = res.reduce((a, r) => (r.life.remainingLife < a.life.remainingLife || (r.life.remainingLife === a.life.remainingLife && r.erf > a.erf) ? r : a), res[0]);
    const gErf = res.reduce((a, r) => (r.erf > a.erf ? r : a), res[0]);
    const gd = { D, t, d: g.d, L: g.L, smys: p.smys };
    const o1 = b31gOriginal({ ...gd, safetyFactor: inp.sf }), o2 = b31gModified({ ...gd, safetyFactor: inp.sf }), o3 = dnvF101Single({ D, t, d: g.d, L: g.L, fu: p.smts, F2: inp.F2 });
    const prof = inp.profile && inp.profile.x && inp.profile.x.length >= 3 ? inp.profile : null;
    let lvl2 = null;
    if (prof) { try { lvl2 = effectiveArea({ D, t, x: prof.x.map((v) => v / 1e3), depth: prof.y.map((v) => Math.min(v / 1e3, 0.95 * t)), smys: p.smys, safetyFactor: inp.sf }); } catch { lvl2 = null; } }
    const ffs = fitnessForService({ D, t, d: g.d, L: g.L, smys: p.smys, maop, smts: p.smts, safetyFactor: inp.sf, F2: inp.F2, ...(prof ? { profile: { x: prof.x.map((v) => v / 1e3), depth: prof.y.map((v) => Math.min(v / 1e3, 0.95 * t)) } } : {}) });
    // dent
    const dent = inp.dentDepth > 0 ? dentDepthCheck({ d: inp.dentDepth / 1e3, D, onWeld: !!inp.dentOnWeld }) : null;
    const dStr = inp.dentDepth > 0 && inp.dentR1 !== 0 && inp.dentR2 !== 0 ? dentStrain({ t, R0: D / 2, R1: inp.dentR1 / 1e3, R2: inp.dentR2 / 1e3, d: inp.dentDepth / 1e3, L: inp.dentLen / 1e3, limit: inp.dentOnWeld ? 0.04 : 0.06 }) : null;
    // extreme-value wall loss
    const bm = inp.blockMax && inp.blockMax.y && inp.blockMax.y.length >= 5 ? inp.blockMax.y.map((v) => v / 1e3) : null;
    const ext = bm ? run(() => wallLossExtreme({ blockMaxima: bm, returnPeriod: inp.blocks, t }), 'Extreme-value wall loss') : null;
    // time-dependent reliability of the governing defect (DNV-RP-F101 capacity with model uncertainty)
    const target = TARGET[inp.safetyClass];
    const vars = (rateMean) => ({
      dRel: { dist: 'normal', mean: g.d / t, std: inp.sdDepth / 100 }, fu: { dist: 'lognormal', mean: 1.08 * p.smts, cov: inp.covFu }, Xm: { dist: 'normal', mean: 1.05, cov: inp.covModel },
      P: { dist: 'gumbel', mean: maop, cov: inp.covP }, ...(rateMean > 0 && inp.covRate > 0 ? { rate: { dist: 'lognormal', mean: rateMean, cov: inp.covRate } } : {}),
    });
    const limit = (yr) => (x) => {
      const d = Math.max(x.dRel, 0) * t + (x.rate ?? dRate) * yr, L = g.L + lRate * yr;
      if (d >= 0.999 * t) return -x.P;
      const Q = Math.sqrt(1 + 0.31 * L * L / (D * t));
      return x.Xm * 2 * t * x.fu / (D - t) * (1 - d / t) / (1 - d / (t * Q)) - x.P;
    };
    const years = linspace(0, Math.max(inp.life, inp.interval), 9);
    const rel = years.map((yr) => { try { const f = form(limit(yr), vars(dRate), { maxIter: 60 }); return { yr, beta: f.beta, pf: f.pf, alpha: f.importance }; } catch { return { yr, beta: -8, pf: 1, alpha: {} }; } });
    const pfNow = rel[0].pf;
    const atInsp = (() => { try { return form(limit(inp.interval), vars(dRate), { maxIter: 60 }); } catch { return { pf: 1, beta: -8, importance: {} }; } })();
    let lifeRel = years[years.length - 1];
    for (let i = 0; i < rel.length; i++) if (rel[i].pf > target) { lifeRel = i === 0 ? 0 : rel[i - 1].yr + (rel[i].yr - rel[i - 1].yr) * (Math.log(target) - Math.log(Math.max(rel[i - 1].pf, 1e-300))) / (Math.log(rel[i].pf) - Math.log(Math.max(rel[i - 1].pf, 1e-300))); break; }
    const imp = Object.entries(atInsp.importance || {}).sort((a, b) => b[1] - a[1]);
    const checks = [
      chk('Safe pressure of the worst defect now', maop, gErf.a.safePressure, 'Pa', inp.method === 'dnv' ? 'DNV-RP-F101 (2019) Part B' : 'ASME B31G-2012 §2.2',
        'Guards against the most severe corrosion defect failing by plastic collapse of the remaining ligament at the operating pressure. Above 1.0 the pressure must be reduced or the defect repaired.'),
      chk('Deepest defect against the 80 % wall limit', Math.max(...res.map((r) => r.dPct)), 80, '% wall', 'ASME B31G-2012 §1.2; DNV-RP-F101 (85 %)',
        'Defects deeper than 80 % of the wall are outside the methods and may leak; they are repaired regardless of the calculated pressure.'),
      chk('Remaining life against the re-inspection interval', inp.interval, Math.max(g.life.remainingLife, 1e-6), 'years', 'ASME B31G §1.6; ASME B31.8S §7; DNV-RP-G101',
        'Guards against a defect growing to failure before the next inspection finds it; the defect must stay acceptable until then with margin.'),
      chk('Annual failure probability at the next inspection', atInsp.pf, target, '-', 'DNV-ST-F101 (2021) Sec. 2 target failure probabilities (recalled); FORM',
        'Compares the probability that the governing defect bursts, with inspection sizing and growth uncertainty, against the target for the safety class.'),
    ];
    if (dent) checks.push(chk('Plain dent depth', dent.depthRatio, dent.limit, '-', 'ASME B31.8 (2022) §851.4.1',
      'Guards against a dent deep enough to restrict pigs and to crack under pressure cycling; dents with gouges or on welds are far more severe.'));
    if (dStr) checks.push(chk('Dent strain', dStr.epsMax, dStr.limit, '-', 'ASME B31.8 (2022) Appendix R',
      'Guards against excessive local strain in the dent, which exhausts ductility and promotes cracking.'));
    const gov = worst(checks);
    // acceptance curve: allowable depth against length at MAOP for the selected method
    const Ls = linspace(0.02, Math.max(2.5 * Math.max(...res.map((r) => r.L)), 0.5), 30);
    const curve = Ls.map((L) => { let lo = 0.01, hi = 0.8; if (assess({ D, t, d: hi * t, L, smys: p.smys }).safePressure >= maop) return 80; if (assess({ D, t, d: lo * t, L, smys: p.smys }).safePressure < maop) return 0; for (let i = 0; i < 30; i++) { const m = 0.5 * (lo + hi); if (assess({ D, t, d: m * t, L, smys: p.smys }).safePressure >= maop) lo = m; else hi = m; } return lo * 100; });
    const ty = linspace(0, Math.min(Math.max(g.life.remainingLife * 1.2, inp.interval * 1.5), 100), 25);
    const pT = ty.map((y) => { const d = Math.min(g.d + dRate * y, 0.95 * t); return assess({ D, t, d, L: g.L + lRate * y, smys: p.smys }).safePressure; });
    const nRepair = res.filter((r) => r.erf > 1 || r.dPct >= 80 || r.life.remainingLife < inp.interval).length;
    const table = res.map((r) => ({ kp: round(r.kp, 4), depth: round(r.dPct), length: round(r.L * 1e3), pf: round(r.a.failurePressure / BAR), ps: round(r.a.safePressure / BAR), erf: round(r.erf, 3), life: round(Math.min(r.life.remainingLife, 100)), next: round(Math.min(r.life.reinspectionInterval, 100)),
      action: r.erf > 1 || r.dPct >= 80 ? 'repair or de-rate now' : r.life.remainingLife < inp.interval ? 'repair before the next inspection' : r.life.remainingLife < inp.life ? 'monitor' : 'acceptable' }))
      .sort((a, b) => a.life - b.life || b.erf - a.erf);
    return {
      kpis: [
        kpi('erf', 'Estimated repair factor of the worst defect', gErf.erf, '-', 'Operating pressure over the safe pressure of the most severe defect. Above 1.0 the defect is not acceptable at MAOP.', 1, 'max'),
        kpi('pSafe', 'Safe operating pressure (governing defect)', gErf.a.safePressure / BAR, 'bar', 'Lowest safe pressure among the listed defects by the selected method: the pressure the line may operate at today.', inp.maop, 'min'),
        kpi('remLife', 'Remaining life of the governing defect', Math.min(g.life.remainingLife, 100), 'years', `Time until the governing defect (${fmt(g.kp)} km) no longer sustains MAOP or reaches 80 % depth, at the stated growth rates (limited by ${g.life.governing}).`, inp.interval, 'min'),
        kpi('reinspect', 'Recommended re-inspection interval', Math.min(g.life.reinspectionInterval, 100), 'years', 'Half the remaining life of the governing defect: inspect again no later than this.'),
        kpi('pDerated', 'Safe pressure at the next planned inspection', Math.min(...res.map((r) => r.pAtInsp)) / BAR, 'bar', 'Lowest safe pressure among the defects after growth over the planned inspection interval.'),
        kpi('pfB31G', 'Failure pressure, original B31G', o1.failurePressure / BAR, 'bar', 'Predicted burst pressure of the governing defect by original B31G (parabolic area, flow stress 1.1 SMYS).'),
        kpi('pfMod', 'Failure pressure, Modified B31G', o2.failurePressure / BAR, 'bar', 'Predicted burst pressure by Modified B31G (0.85 dL area, flow stress SMYS + 69 MPa).'),
        kpi('pfDnv', 'Failure pressure, DNV-RP-F101', o3.failurePressure / BAR, 'bar', 'Predicted burst pressure by the DNV-RP-F101 capacity equation based on tensile strength.'),
        kpi('pfLevel2', 'Failure pressure, effective area', lvl2 ? lvl2.failurePressure / BAR : 0, 'bar', 'Level 2 (RSTRENG) result from the river-bottom profile; zero when no profile is supplied.'),
        kpi('pfNow', 'Annual failure probability now', pfNow, '-', 'Probability that the governing defect fails at the annual maximum pressure today, with sizing, strength and model uncertainty.', target, 'max'),
        kpi('pfInsp', 'Annual failure probability at the next inspection', atInsp.pf, '-', 'Same probability after growth over the planned inspection interval, including growth-rate uncertainty.', target, 'max'),
        kpi('betaInsp', 'Reliability index at the next inspection', atInsp.beta, '-', 'FORM reliability index; 3.7 corresponds to 1e-4 and 4.3 to 1e-5.'),
        kpi('lifeRel', 'Time to reach the target failure probability', lifeRel, 'years', 'Reliability-based remaining life: when the annual failure probability of the governing defect exceeds the target for the safety class.'),
        kpi('nRepair', 'Defects needing action before the next inspection', nRepair, '-', 'Number of listed defects that are unacceptable now or would become so before the planned inspection.'),
        kpi('dentRatio', 'Dent depth / diameter', dent ? dent.depthRatio * 100 : 0, '%', 'Plain dents are limited to 6 % of the diameter (2 % on welds).', dent ? dent.limit * 100 : undefined, 'max'),
        kpi('wallMax', 'Extreme-value wall loss over the whole line', ext ? ext.predictedMax * 1e3 : 0, 'mm', 'Deepest pit expected in the uninspected line from the Gumbel fit of block maxima; zero when no block data are supplied.'),
      ],
      checks,
      plots: [
        plot('scatter', 'Defect acceptance at MAOP', 'Defect length', 'mm', 'Defect depth', '% wall', [
          { name: 'Acceptance limit', x: Ls.map((v) => round(v * 1e3)), y: curve.map((v) => round(v)) }, { name: 'Reported defects', x: res.map((r) => round(r.L * 1e3)), y: res.map((r) => round(r.dPct)) },
          { name: `After ${fmt(inp.interval)} years of growth`, x: res.map((r) => round((r.L + lRate * inp.interval) * 1e3)), y: res.map((r) => round(Math.min(r.dPct + dRate * inp.interval / t * 100, 100))) },
        ], 'Defects below the curve are acceptable at MAOP with the safety factor. The second set of points shows where each defect will be at the next inspection: any that cross the curve need action first.'),
        plot('line', 'Safe pressure of the governing defect vs time', 'Time', 'years', 'Safe pressure', 'bar', [
          { name: 'Safe pressure', x: ty.map((v) => round(v)), y: pT.map((v) => round(v / BAR)) }, { name: 'MAOP', x: [0, round(ty[ty.length - 1])], y: [inp.maop, inp.maop], dash: true },
        ], 'Decline of the safe pressure as the governing defect grows. The crossing with MAOP is the remaining life; re-inspect at about half that time.'),
        plot('line', 'Failure probability of the governing defect', 'Time', 'years', 'log10 annual failure probability', '-', [
          { name: 'FORM', x: rel.map((r) => round(r.yr)), y: rel.map((r) => round(Math.log10(Math.max(r.pf, 1e-30)))) }, { name: 'Target', x: [0, round(years[years.length - 1])], y: [Math.log10(target), Math.log10(target)], dash: true },
        ], 'Annual failure probability with inspection sizing and growth uncertainty. Where it rises through the target line the defect must have been repaired or re-inspected.'),
      ],
      tables: [
        tbl('Critical segments: defects ranked by remaining life', [col('kp', 'Chainage', 'km'), col('depth', 'Depth', '% wall'), col('length', 'Length', 'mm'), col('pf', 'Failure pressure', 'bar'), col('ps', 'Safe pressure', 'bar'), col('erf', 'ERF', '-'), col('life', 'Remaining life', 'years'), col('next', 'Re-inspect within', 'years'), col('action', 'Action')], table.slice(0, 200),
          'Each row is a reported defect, most critical first: its failure and safe pressure, repair factor (above 1 = not acceptable), remaining life at the stated growth and the action it calls for.'),
        tbl('Assessment levels for the governing defect', [col('name', 'Method'), col('ps', 'Safe pressure', 'bar'), col('uc', 'MAOP / safe pressure', '-'), col('ok', 'Acceptable')], ffs.levels.map((l) => ({ name: l.name, ps: round(l.safePressure / BAR), uc: round(l.uc, 3), ok: l.acceptable ? 'yes' : 'no' })),
          'The same defect assessed by each available method, from the most to the least conservative; a higher level may accept a defect that Level 1 rejects.'),
        tbl('Sensitivity of the failure probability', [col('v', 'Variable'), col('imp', 'Importance (α²)', '-')], imp.map(([k, v]) => ({ v: { dRel: 'Defect depth (sizing error)', fu: 'Tensile strength', Xm: 'Burst model uncertainty', P: 'Operating pressure', rate: 'Corrosion growth rate' }[k] || k, imp: round(v, 3) })),
          'Share of each uncertain variable in the failure probability at the next inspection. The largest entry shows what better data would improve most.'),
      ],
      notes: [`Governing defect at ${fmt(g.kp)} km: ${fmt(g.dPct)} % deep, ${fmt(g.L * 1e3)} mm long; remaining life ${fmt(Math.min(g.life.remainingLife, 100))} years (limited by ${g.life.governing}).`,
        `Fitness-for-service summary: ${ffs.action}`,
        dent ? `Dent: ${dent.note}${dStr ? ` Maximum strain ${fmt(dStr.epsMax * 100)} %.` : ''}` : 'No dent reported.',
        ext ? `Gumbel extrapolation: deepest pit over the whole line ${fmt(ext.predictedMax * 1e3)} mm, probability of exceeding the wall ${fmt(ext.exceedanceOfWall)}.` : 'No block-maxima data: extreme-value extrapolation not performed.'],
      advice: [
        adv('engineering', sev(gov.uc), `Worst defect ERF ${fmt(gErf.erf)} (safe pressure ${fmt(gErf.a.safePressure / BAR)} bar against MAOP ${fmt(inp.maop)} bar); governing remaining life ${fmt(Math.min(g.life.remainingLife, 100))} years; failure probability at the next inspection ${fmt(atInsp.pf)} (target ${fmt(target)}).`,
          gov.uc > 1 ? (/[Dd]ent/.test(gov.name) ? 'The dent exceeds the screening limit: inspect it for cracks and gouges, assess it by a strain-based or fatigue method, and repair it (sleeve or cut-out) if it cannot be justified.' : gErf.erf > 1 ? `De-rate to ${fmt(gErf.a.safePressure / BAR)} bar or repair the defect at ${fmt(gErf.kp)} km (clamp or cut-out); ${nRepair} defect(s) need action before the next inspection.` : `Shorten the inspection interval to ${fmt(Math.min(g.life.reinspectionInterval, lifeRel))} years or repair the defect at ${fmt(g.kp)} km; ${nRepair} defect(s) need action.`)
            : `All defects are acceptable to the next inspection in ${fmt(inp.interval)} years. Re-inspect within ${fmt(Math.min(g.life.reinspectionInterval, 100))} years and confirm the growth rate by run-to-run comparison.`,
          'ASME B31G-2012; DNV-RP-F101 (2019); ASME B31.8S §7; DNV-RP-G101 (risk-based inspection)'),
        adv('economic', 'info', `${nRepair} of ${res.length} defects need intervention before the planned inspection; the governing defect allows an inspection interval of up to ${fmt(Math.min(g.life.reinspectionInterval, 100))} years.`,
          nRepair > 0 ? 'A repair clamp is typically several vessel days per location; a temporary pressure de-rating costs deferred production instead — rank the two by total cost over the period to the next shutdown.' : `Extending the inspection interval from ${fmt(inp.interval)} towards ${fmt(Math.min(g.life.reinspectionInterval, lifeRel, 100))} years saves an in-line inspection campaign; the largest contributor to the failure probability (${imp[0] ? ({ dRel: 'depth sizing error', fu: 'tensile strength', Xm: 'model uncertainty', P: 'pressure', rate: 'growth rate' }[imp[0][0]] || imp[0][0]) : 'n/a'}) shows which data would justify it.`,
          RBI_BASIS, `${vesselCo2(Math.max(nRepair, 1) * 4)} Each avoided cut-out saves about ${fmt(steelPerMetre(p, p.t) * 12 / 1e3)} t of replacement pipe (${fmt(steelPerMetre(p, p.t) * 12 / 1e3 * CO2_STEEL)} t CO₂e, screening factor) and the hydrocarbon release risk of an unplanned failure.`),
      ],
      publish: { 'pipeline.remainingLife': Math.min(g.life.remainingLife, 100), 'pipeline.safePressure': gErf.a.safePressure },
      quantities: { vessel_days: { repair: nRepair * 4 } },
    };
  },
  verify: [{
    // D = 500 mm, t = 10 mm, d = 5 mm, L = 200 mm, X52 (SMYS 360 MPa): z = L²/(D t) = 8.
    // Original B31G: M = √(1 + 0.8 z) = 2.7203; S_F = 1.1·360·(1 − ⅔·0.5)/(1 − ⅔·0.5/2.7203) = 300.9 MPa; P_F = 2 S_F t/D = 12.03 MPa = 120.3 bar.
    // Modified B31G: M = √(1 + 0.6275 z − 0.003375 z²) = 2.4091; S_F = (360 + 69)·(1 − 0.85·0.5)/(1 − 0.85·0.5/2.4091) = 299.5 MPa; P_F = 119.8 bar.
    name: 'ASME B31G original and Modified B31G failure pressure of a 50 % deep defect',
    source: 'ASME B31G-2012 §2.2; Kiefner & Vieth (1989); hand calculation',
    inputs: { D: 500, t: 10, steel: 'api5l-x52', maop: 60, defects: [{ kp: 1, depth: 50, length: 200 }], depthRate: 0.1, lengthRate: 0, dentDepth: 0 },
    expect: [{ kpi: 'pfB31G', value: 120.3, tol: 0.003 }, { kpi: 'pfMod', value: 119.8, tol: 0.003 }],
  }],
};
