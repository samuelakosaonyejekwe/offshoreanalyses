// Station-keeping reliability: probability of line failure in the annual extreme storm by FORM, SORM, Monte Carlo
// and importance sampling on a response surface of the mooring system, and the probability of losing station
// through a second failure in the damaged condition.
import { form, sorm, importanceSampling } from '../../engine/reliability/form.js';
import { monteCarlo } from '../../engine/reliability/sampling.js';
import { normInv } from '../../engine/math/stats.js';
import { need, run, round, fmt, adv, check, sev, linspace, interp1, systemInputs, buildSystem, fastSystem, DEG, LIMIT_QS, LIMIT_PROPS } from './_common.js';

const sysMode = { key: 'loadModel', in: ['system'] }, dirMode = { key: 'loadModel', in: ['direct'] };

export default {
  id: 's10-reliability',
  suite: 10,
  title: 'Station-keeping reliability',
  summary: 'Estimates the annual probability that a mooring line breaks in the extreme storm and that a second line follows, and shows which uncertainties drive it.',
  covers: [15],
  fidelity: 'reduced-order',
  standards: ['DNV-OS-E301 (2018) Ch.2 Sec.2 (target annual failure probabilities behind the partial factors)', 'ISO 2394 (2015) General principles on reliability for structures', 'DNV-RP-C210 (2019) Probabilistic methods for inspection planning (reliability methods)', 'ISO 19901-7 (2013) Annex A (reliability basis)'],
  models: ['First-order reliability method (Hasofer–Lind–Rackwitz–Fiessler)', 'Second-order reliability method (Breitung)', 'Crude Monte Carlo simulation and design-point importance sampling', 'Gumbel annual-maximum environmental force', 'Quasi-static response surface of line tension vs mean force (intact and one line broken)', 'Series system of lines and conditional second-failure (progressive failure) probability'],
  limits: `The limit state is line breakage g = strength − tension, with tension from a quasi-static response surface along one load direction and a random multiplier for low- and wave-frequency dynamics; anchors, connectors, fatigue and human or operational failures are not in the model, so the result is a notional reliability for comparing designs, not a predicted failure frequency. Directional and seasonal statistics, correlation between lines and time-dependent degradation are not modelled. ${LIMIT_QS} ${LIMIT_PROPS}`,
  searchTerms: 'mooring system reliability FORM Monte Carlo annual failure probability',
  inputs: [
    { key: 'loadModel', label: 'Tension model', type: 'select', default: 'system', options: [{ value: 'system', label: 'Mooring system response surface' }, { value: 'direct', label: 'Direct: strength and tension distributions' }], group: 'Model', role: 'control',
      help: 'System: line tension is computed from the mooring layout for a random annual-maximum force. Direct: enter the statistics of line strength and annual-maximum tension yourself.' },
    ...systemInputs().map((i) => ({ ...i, showIf: i.showIf ?? sysMode })),
    { key: 'loadDirDeg', label: 'Direction the load acts towards', unit: 'deg', type: 'number', default: 0, min: -180, max: 360, group: 'Environmental load', role: 'bc', showIf: sysMode, help: 'Governing direction from the extreme-environment module (boundary condition: load applied to the vessel).' },
    { key: 'Fchar', label: 'Characteristic (100-year) mean force', unit: 'kN', type: 'number', default: 3000, min: 10, max: 60000, group: 'Environmental load', role: 'bc', bind: 'vessel.Fmean', showIf: sysMode,
      help: 'Mean environmental force with an annual exceedance probability of 1 %. The annual maximum is modelled as Gumbel with this 100-year value.' },
    { key: 'Fcov', label: 'Coefficient of variation of the annual-maximum force', unit: '-', type: 'number', default: 0.25, min: 0.05, max: 0.8, group: 'Environmental load', role: 'environment', showIf: sysMode,
      help: 'Year-to-year variability of the largest storm load: about 0.15–0.25 for North Sea-type climates, 0.4 or more for hurricane-dominated areas.' },
    { key: 'dynMean', label: 'Dynamic factor on load-induced tension, mean', unit: '-', type: 'number', default: 1.4, min: 1, max: 5, group: 'Environmental load', role: 'load', showIf: sysMode,
      help: 'Ratio of maximum to mean tension increase above pretension (low- and wave-frequency dynamics): Tmax = T0 + factor × (Tmean − T0). From the tension-and-offset module.' },
    { key: 'dynCov', label: 'Dynamic factor, coefficient of variation', unit: '-', type: 'number', default: 0.15, min: 0, max: 0.6, group: 'Environmental load', role: 'load', showIf: sysMode },
    { key: 'Smean', label: 'Annual-maximum line tension, mean', unit: 'kN', type: 'number', default: 5000, min: 1, max: 60000, group: 'Direct model', role: 'load', showIf: dirMode },
    { key: 'Scov', label: 'Annual-maximum line tension, coefficient of variation', unit: '-', type: 'number', default: 0.3, min: 0.01, max: 1, group: 'Direct model', role: 'load', showIf: dirMode },
    { key: 'Rmean', label: 'Line strength, mean', unit: 'kN', type: 'number', default: 9000, min: 1, max: 80000, group: 'Direct model', role: 'material', showIf: dirMode },
    { key: 'distType', label: 'Distribution type (direct model)', type: 'select', default: 'normal', options: [{ value: 'normal', label: 'Normal strength and tension' }, { value: 'lognormal', label: 'Lognormal strength and tension' }, { value: 'gumbel', label: 'Lognormal strength, Gumbel tension' }], group: 'Direct model', role: 'control', showIf: dirMode },
    { key: 'strengthBias', label: 'Mean strength / catalogue MBL', unit: '-', type: 'number', default: 1.1, min: 0.8, max: 1.5, group: 'Resistance', role: 'material', showIf: sysMode,
      help: 'Actual break strength of new chain and rope typically exceeds the catalogue minimum by 5–15 %; reduce it for corroded or aged lines (initial condition of the components).' },
    { key: 'Rcov', label: 'Line strength, coefficient of variation', unit: '-', type: 'number', default: 0.07, min: 0.01, max: 0.4, group: 'Resistance', role: 'material', help: 'About 0.05 for new chain, 0.05–0.10 for wire and fibre rope, higher for degraded components.' },
    { key: 'modelCov', label: 'Model uncertainty on tension, coefficient of variation', unit: '-', type: 'number', default: 0.1, min: 0, max: 0.5, group: 'Resistance', role: 'criteria', showIf: sysMode, help: 'Lognormal multiplier (median 1) on the calculated tension for analysis-method uncertainty.' },
    { key: 'pfTarget', label: 'Target annual failure probability (single line)', unit: '-', type: 'number', default: 1e-4, min: 1e-7, max: 0.1, group: 'Acceptance criteria', role: 'criteria',
      help: 'DNV-OS-E301 calibrated its ULS factors to about 10⁻⁴ per year for consequence class 1 and 10⁻⁵ for class 2. Confirm the target with the governing standard and the operator\'s risk criteria.' },
    { key: 'pfLossTarget', label: 'Target annual probability of loss of station', unit: '-', type: 'number', default: 1e-5, min: 1e-8, max: 0.1, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'nMC', label: 'Monte Carlo samples', type: 'int', default: 20000, min: 1000, max: 400000, group: 'Numerical', role: 'control', help: 'Crude Monte Carlo resolves probabilities down to about 10 / samples; smaller values come from FORM/SORM and importance sampling.' },
    { key: 'seed', label: 'Random seed', type: 'int', default: 1, min: 1, max: 99999, group: 'Numerical', role: 'control' },
  ],
  compute(inp) {
    let g, vars, labels, b = null, nLines = 1, surf = null, gDam = null, critName = 'line', MBL = 0;
    if (inp.loadModel === 'direct') {
      const R = { mean: inp.Rmean * 1e3, cov: inp.Rcov }, S = { mean: inp.Smean * 1e3, cov: inp.Scov };
      vars = { R: { dist: inp.distType === 'normal' ? 'normal' : 'lognormal', ...R }, S: { dist: inp.distType === 'normal' ? 'normal' : inp.distType === 'gumbel' ? 'gumbel' : 'lognormal', ...S } };
      g = (x) => x.R - x.S;
      labels = { R: 'Line strength', S: 'Annual-maximum tension' };
      MBL = inp.Rmean * 1e3;
    } else {
      b = buildSystem(inp); nLines = b.nLines; MBL = b.MBL;
      const th = inp.loadDirDeg * DEG, cd = Math.cos(th), sd = Math.sin(th);
      const fs = fastSystem(b.system, { maxOffset: Math.max(0.2 * b.depth, 40) });
      const calm = fs.solve([0, 0, 0, 0, 0, 0]);
      // Gumbel annual maximum with the 100-year value Fchar: x100 = loc + scale·y, y = −ln(−ln 0.99); mean = loc + 0.5772·scale; std = π·scale/√6
      const y100 = -Math.log(-Math.log(0.99)), cs = Math.PI / Math.sqrt(6);
      const Fc = inp.Fchar * 1e3, mean = Fc / (1 + inp.Fcov * (y100 - 0.5772156649) / cs), std = inp.Fcov * mean;
      // response surface: largest mean line tension vs force, intact and with the critical line removed
      const Fgrid = linspace(0, Math.max(3.0 * Fc, mean + 12 * std), 28);
      const build = (exclude) => {
        const F = [], T = []; let start = [calm.x, calm.y];
        for (const f of Fgrid) {
          let r; try { r = fs.solve([f * cd, f * sd, 0, 0, 0, 0], start, exclude); } catch (e) { break; }
          if (r.maxTension > 3 * MBL) { F.push(f); T.push(r.maxTension); break; }
          F.push(f); T.push(r.maxTension); start = [r.x, r.y];
        }
        need(F.length >= 3, 'The mooring cannot balance even small loads in this configuration: check the layout before a reliability analysis.');
        return { F, T, T0: T[0] };
      };
      const intact = build([]);
      const crit = fs.solve([Fc * cd, Fc * sd, 0, 0, 0, 0], [calm.x, calm.y]).maxTensionLine;
      critName = b.system.lines[crit].name;
      const dam = build([crit]);
      surf = { intact, dam, Fc, mean, std };
      const tension = (s, F) => {
        // beyond the last solved point the line is at its elastic limit: extrapolate with the end slope (steep)
        const n = s.F.length;
        if (F <= s.F[n - 1]) return interp1(s.F, s.T, Math.max(F, 0));
        return s.T[n - 1] + (s.T[n - 1] - s.T[n - 2]) / (s.F[n - 1] - s.F[n - 2]) * (F - s.F[n - 1]);
      };
      vars = { R: { dist: 'lognormal', mean: inp.strengthBias * MBL, cov: inp.Rcov }, F: { dist: 'gumbel', mean, std } };
      if (inp.dynCov > 0) vars.D = { dist: 'lognormal', mean: inp.dynMean, cov: inp.dynCov };
      if (inp.modelCov > 0) vars.X = { dist: 'lognormal', mean: Math.sqrt(1 + inp.modelCov ** 2), cov: inp.modelCov };   // median 1
      const gOf = (s) => (x) => { const Tm = tension(s, x.F); return x.R - (x.X ?? 1) * (s.T0 + (x.D ?? inp.dynMean) * Math.max(Tm - s.T0, 0)); };
      g = gOf(intact); gDam = gOf(dam);
      labels = { R: 'Line strength', F: 'Annual-maximum environmental force', D: 'Dynamic tension factor', X: 'Model uncertainty' };
    }
    const f1 = run(() => form(g, vars), 'FORM');
    need(Number.isFinite(f1.beta), 'FORM did not return a finite reliability index: check the distributions.');
    let s2 = null; try { s2 = sorm(g, vars); } catch (e) { s2 = null; }
    const mc = run(() => monteCarlo(g, vars, { n: inp.nMC, seed: inp.seed, keepInputs: false }), 'Monte Carlo');
    let is = null; try { is = importanceSampling(g, vars, { designPoint: f1.designPoint, n: 4000, seed: inp.seed }); } catch (e) { is = null; }
    const pfForm = f1.pf, pfSorm = s2 && s2.valid ? s2.pfBreitung : pfForm;
    const pfBest = is && is.cov < 0.2 && is.pf > 0 ? is.pf : pfSorm;
    // series system of lines in the same storm: bounded by the governing line (fully correlated) and n × pf (independent)
    const pfAnyUpper = Math.min(1, nLines * pfBest);
    // one-line-broken exposure: the same limit state on the response surface of the damaged system. (A storm that overloads the
    // intact system would overload the damaged one too, so the informative number is the overload probability while a line is
    // already missing for another reason.)
    let p2 = 0;
    if (gDam) p2 = run(() => form(gDam, vars), 'FORM (damaged)').pf;
    const betaT = -normInv(inp.pfTarget);
    const imp = Object.keys(vars).map((k) => ({ k, name: labels[k] || k, a2: f1.importance[k] ?? 0, alpha: f1.alpha[k] ?? 0, x: f1.designPoint[k] }));
    imp.sort((p, q) => q.a2 - p.a2);
    const betaBest = -normInv(Math.min(Math.max(pfBest, 1e-300), 1 - 1e-16));
    const checks = [
      check('Annual probability of line failure', pfBest, inp.pfTarget, '-', 'DNV-OS-E301 (2018) Ch.2 Sec.2 calibration target; ISO 2394 (2015) §8',
        'Compares the notional annual probability that the critical line breaks in the extreme storm with the target. Exceeding it means the design is less reliable than the code calibration assumes.'),
    ];
    if (gDam) checks.push(check('Annual probability of overload with one line already broken', p2, Math.max(inp.pfLossTarget / 1e-2, inp.pfTarget), '-', 'ISO 19901-7 (2013) §8.9 (redundancy); target = loss-of-station target ÷ assumed 10⁻² annual line-loss rate from other causes',
      'Probability that, with the critical line already missing for a year, a remaining line is overloaded. Multiplied by the rate of line loss from fatigue, wear or damage it gives the chance of losing station by a second failure.'));
    const gov = checks.reduce((a, c) => (c.uc > a.uc ? c : a), checks[0]);
    const lossViaDamage = gDam ? 1e-2 * p2 : 0;
    const pfStation = Math.min(1, pfBest + lossViaDamage);
    // strength change needed to meet the target (lognormal-ish scaling: shift in ln R by Δβ·sqrt(σlnR² + σlnS²) ≈ Δβ·σg/|∇|)
    const sigLn = Math.sqrt(Math.log(1 + inp.Rcov ** 2) + (inp.loadModel === 'direct' ? Math.log(1 + inp.Scov ** 2) : Math.log(1 + inp.Fcov ** 2) + Math.log(1 + inp.dynCov ** 2) + Math.log(1 + inp.modelCov ** 2)));
    const strengthFactor = Math.exp((betaT - betaBest) * sigLn);
    return {
      kpis: [
        { key: 'beta', label: 'Reliability index β (FORM), annual', value: round(f1.beta, 5), unit: '-', limit: round(betaT, 4), sense: 'min', explain: 'Distance from the mean point to the failure surface in standard-normal space. β = 3.1 corresponds to 10⁻³ per year, 3.7 to 10⁻⁴, 4.3 to 10⁻⁵.' },
        { key: 'pfForm', label: 'Failure probability (FORM)', value: Number(pfForm.toPrecision(4)), unit: '1/yr', explain: 'Φ(−β): first-order estimate of the annual probability that the critical line\'s tension exceeds its strength.' },
        { key: 'pfSorm', label: 'Failure probability (SORM)', value: Number(pfSorm.toPrecision(4)), unit: '1/yr', explain: 'Second-order correction for the curvature of the failure surface; close to FORM when the limit state is nearly linear.' },
        { key: 'pfMC', label: 'Failure probability (Monte Carlo)', value: Number(mc.pf.toPrecision(4)), unit: '1/yr', explain: `Fraction of ${inp.nMC} random years in which the line fails. Zero means the probability is below what this sample size can resolve (about ${fmt(10 / inp.nMC)}).` },
        { key: 'pfIS', label: 'Failure probability (importance sampling)', value: Number((is ? is.pf : pfSorm).toPrecision(4)), unit: '1/yr', explain: 'Sampling concentrated around the design point; an independent check on FORM/SORM that works at small probabilities.' },
        { key: 'pf', label: 'Annual probability of line failure (best estimate)', value: Number(pfBest.toPrecision(4)), unit: '1/yr', limit: inp.pfTarget, sense: 'max', explain: 'Importance-sampling value when it has converged, otherwise SORM. This is the single-line reliability indicator to compare with the target.' },
        { key: 'pfAnyLine', label: 'Upper bound for failure of any line', value: Number(pfAnyUpper.toPrecision(4)), unit: '1/yr', explain: 'Series-system bound n × pf (lines treated as independent). The truth lies between the single-line value and this bound because all lines see the same storm.' },
        { key: 'pfDamaged', label: 'Annual probability of overload with one line broken', value: Number(p2.toPrecision(4)), unit: '1/yr', explain: 'Same limit state evaluated on the damaged system (critical line removed): how exposed the unit is while a line is missing.' },
        { key: 'pfStation', label: 'Station-keeping failure indicator', value: Number(pfStation.toPrecision(4)), unit: '1/yr', limit: Math.max(inp.pfTarget, inp.pfLossTarget), sense: 'max', explain: 'Overload of an intact line plus an assumed 10⁻² per year line loss from other causes followed by overload in the damaged condition. A comparative indicator of station-keeping performance.' },
        { key: 'strengthNeeded', label: 'Strength multiplier to reach the target', value: round(strengthFactor, 4), unit: '-', explain: 'Approximate factor on line strength (MBL) that would bring the single-line reliability exactly to the target; below 1 means reserve.' },
        { key: 'mcMeanMargin', label: 'Mean safety margin / MBL', value: round(mc.mean / MBL, 4), unit: '-', explain: 'Average of strength minus tension over the simulated years, as a share of the break load. The spread around it, not the mean, decides the reliability.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Importance of the uncertain variables (FORM α²)', categories: imp.map((v) => v.name), series: [{ name: 'α² (%)', y: imp.map((v) => round(100 * v.a2, 4)) }],
          explain: 'Share of the failure probability attributable to each uncertainty. The tallest bar is where better data (metocean statistics, break tests, dynamic analysis) improves the estimate most.' },
        { type: 'hist', title: 'Distribution of the safety margin (Monte Carlo)', values: Array.from(mc.samples.subarray(0, Math.min(mc.samples.length, 4000)), (v) => round(v / MBL, 4)), bins: 30, x: { label: '(Strength − tension) / MBL', unit: '-' },
          explain: 'Histogram of strength minus tension in each simulated year, normalised by the break load. Failure is the tail to the left of zero; a long left tail signals rare, very large storm loads.' },
      ],
      tables: [
        { title: 'Design point (most likely failure combination)', columns: [{ key: 'name', label: 'Variable' }, { key: 'x', label: 'Value at design point' }, { key: 'alpha', label: 'Direction cosine α' }, { key: 'a2', label: 'Importance α²', unit: '%' }],
          rows: imp.map((v) => ({ name: v.name, x: v.k === 'D' || v.k === 'X' ? round(v.x, 4) : `${fmt(v.x / 1e3, 4)} kN`, alpha: round(v.alpha, 3), a2: round(100 * v.a2, 3) })),
          explain: 'The most probable combination of values at which the line fails. Positive α marks load-type variables (failure when high), negative α resistance-type variables (failure when low).' },
        { title: 'Reliability methods compared', columns: [{ key: 'm', label: 'Method' }, { key: 'pf', label: 'Annual pf' }, { key: 'beta', label: 'β' }, { key: 'note', label: 'Note' }],
          rows: [{ m: 'FORM', pf: Number(pfForm.toPrecision(3)), beta: round(f1.beta, 4), note: `${f1.iterations} iterations${f1.converged ? '' : ' (not converged)'}` },
            { m: 'SORM (Breitung)', pf: Number(pfSorm.toPrecision(3)), beta: round(-normInv(Math.min(Math.max(pfSorm, 1e-300), 1 - 1e-16)), 4), note: s2 && s2.valid ? 'curvature correction applied' : 'not applicable: FORM value used' },
            { m: 'Monte Carlo', pf: Number(mc.pf.toPrecision(3)), beta: round(Number.isFinite(mc.beta) ? mc.beta : 99, 4), note: `${mc.nFail} failures in ${mc.n}; 95 % interval ${Number(mc.ciPf[0].toPrecision(2))} – ${Number(mc.ciPf[1].toPrecision(2))}` },
            { m: 'Importance sampling', pf: Number((is ? is.pf : 0).toPrecision(3)), beta: round(is && Number.isFinite(is.beta) ? is.beta : 99, 4), note: is ? `estimator c.o.v. ${fmt(is.cov)}` : 'not available' }],
          explain: 'Four independent estimates of the same probability. Agreement between them (within the Monte Carlo interval where it has failures) is the check that the result is numerically sound.' },
      ],
      notes: [`β = ${fmt(f1.beta)} (annual), pf ≈ ${Number(pfBest.toPrecision(3))} per year against a target of ${inp.pfTarget}; most important variable: ${imp[0].name} (${fmt(100 * imp[0].a2)} %).`,
        inp.loadModel === 'system' ? `Critical line ${critName}; annual-maximum force Gumbel with mean ${fmt(surf.mean / 1e3)} kN, c.o.v. ${inp.Fcov}, 100-year value ${inp.Fchar} kN. With that line removed the annual overload probability is ${Number(p2.toPrecision(3))}.` : 'Direct model: strength and tension distributions as entered.',
        'Notional reliability for comparison of design options; it excludes anchors, connectors, fatigue and operational causes.'],
      advice: [
        adv('engineering', sev(gov.uc), `Annual probability of line failure ≈ ${Number(pfBest.toPrecision(3))} (β = ${fmt(betaBest)}) against a target of ${inp.pfTarget} (β = ${fmt(betaT)}); ${imp[0].name.toLowerCase()} carries ${fmt(100 * imp[0].a2)} % of the uncertainty.`,
          pfBest > inp.pfTarget ? `Raise line strength by about ${fmt(100 * (strengthFactor - 1))} % (diameter +${fmt(100 * (Math.sqrt(strengthFactor) - 1))} %), or reduce the dominant uncertainty: ${imp[0].k === 'F' || imp[0].k === 'S' ? 'longer metocean records or a site-specific hindcast narrow the extreme-load distribution' : imp[0].k === 'D' ? 'a coupled time-domain analysis narrows the dynamic factor' : imp[0].k === 'R' ? 'break-test certificates narrow the strength distribution' : 'benchmark the analysis against model tests'}.`
            : `The design meets the target with a strength reserve of about ${fmt(100 * (1 - strengthFactor))} %. Keep the assumptions valid in service: confirm corrosion and wear do not erode the mean strength (bias ${inp.loadModel === 'system' ? inp.strengthBias : 1}).`,
          'ISO 2394 (2015) §8 and Annex; DNV-OS-E301 (2018) Ch.2 Sec.2 (reliability basis of the partial factors); DNV-RP-C210 (2019) for the methods.'),
        adv('economic', 'info', gDam ? `With one line missing the annual overload probability rises from ${Number(pfBest.toPrecision(2))} to ${Number(p2.toPrecision(2))} (${fmt(p2 / Math.max(pfBest, 1e-300))} ×).` : `Reliability changes by about one order of magnitude for each ${fmt(100 * (Math.exp(0.55 * sigLn) - 1))} % of strength at this level of uncertainty.`,
          gDam ? `Use the ratio to set the repair urgency: every month spent in the damaged condition carries about ${Number((p2 / 12).toPrecision(2))} probability of a second failure. Compare the cost of a held spare line and a call-off vessel with that exposure multiplied by the consequence of losing station.` : 'Use the sensitivity to trade line strength against data: a metocean study or break-test programme that halves the dominant uncertainty can be worth more than a size step in chain.',
          'Risk-based decision making per ISO 31000 (2018) and ISO 2394 (2015) Annex G (economic optimisation of target reliability); value-of-information analysis.',
          `${strengthFactor >= 1 ? 'Reaching the target by strength alone adds' : 'The reserve over the target corresponds to'} about ${fmt(100 * Math.abs(strengthFactor - 1))} % of line strength${b ? ` (≈ ${fmt(Math.abs(b.steel_t * (strengthFactor - 1)))} t of line steel ≈ ${fmt(Math.abs(b.steel_t * (strengthFactor - 1)) * 1.9)} t CO₂e at a screening factor of 1.9 t CO₂e/t)` : ''}; reducing uncertainty with data achieves the same reliability without the material, and a loss of station carries spill and seabed-damage consequences that dwarf either.`),
      ],
      quantities: b ? { steel_t: round(b.steel_t) } : {},
    };
  },
  verify: [
    { // Normal strength R ~ N(9000, 900) kN and normal tension S ~ N(5000, 1500) kN: β = (9000 − 5000)/sqrt(900² + 1500²)
      // = 4000/1749.3 = 2.2866, pf = Φ(−2.2866) = 0.01111. Monte Carlo with 100 000 samples must agree within sampling error.
      name: 'Normal strength and load: β = (μR − μS)/√(σR² + σS²)',
      source: 'Cornell (1969) reliability index; Melchers & Beck, Structural Reliability Analysis and Prediction, §1.4',
      inputs: { loadModel: 'direct', distType: 'normal', Rmean: 9000, Rcov: 0.1, Smean: 5000, Scov: 0.3, nMC: 100000 },
      expect: [{ kpi: 'beta', value: 2.2866, tol: 0.001 }, { kpi: 'pfForm', value: 0.01111, tol: 0.005 }, { kpi: 'pfMC', value: 0.01111, tol: 0.08 }, { kpi: 'pfIS', value: 0.01111, tol: 0.1 }] },
    { // Lognormal R and S: β = ln[(μR/μS)·sqrt((1 + VS²)/(1 + VR²))] / sqrt(ln[(1 + VR²)(1 + VS²)]) with μR/μS = 2, VR = 0.1, VS = 0.3:
      // numerator ln(2·sqrt(1.09/1.01)) = ln(2.0777) = 0.73126; denominator sqrt(ln(1.1009)) = 0.31006 → β = 2.3585.
      name: 'Lognormal strength and load (closed form)',
      source: 'Melchers & Beck, Structural Reliability Analysis and Prediction, §1.4.3 (lognormal safety index)',
      inputs: { loadModel: 'direct', distType: 'lognormal', Rmean: 10000, Rcov: 0.1, Smean: 5000, Scov: 0.3 },
      expect: [{ kpi: 'beta', value: 2.3585, tol: 0.002 }] },
  ],
};
