// Time-dependent structural reliability and service life of a coastal structure: annual failure probability by
// FORM for a chosen limit state (armour stability, overtopping, caisson sliding or a generic R − S), tracked
// through the life as the sea rises and the structure deteriorates, with a Monte Carlo cross-check.
import { form } from '../../engine/reliability/form.js';
import { monteCarlo } from '../../engine/reliability/sampling.js';
import { vanDerMeerRock, relativeDensity } from '../../engine/coastal/armour.js';
import { overtoppingRubbleMound } from '../../engine/coastal/runup.js';
import { godaPressure } from '../../engine/coastal/caisson.js';
import { normCdf, normInv } from '../../engine/math/stats.js';
import { G, CO2_ROCK, CO2_CONC, EDITION_NOTE, need, run, finite, r, f, linspace, check, eng, eco, slrInputs, slrRows, SLR_BASIS } from './_common.js';

const LS = { armour: 'Armour stability (Van der Meer)', overtopping: 'Overtopping discharge (EurOtop)', sliding: 'Caisson sliding (Goda)', generic: 'Generic resistance − load' };

export default {
  id: 's07-reliability-service-life', suite: 7, title: 'Reliability and service life',
  summary: 'Annual and lifetime probability of failure of a coastal structure for a chosen limit state, how it grows with sea-level rise and deterioration, and the year the target reliability is lost.',
  covers: [12], fidelity: 'reduced-order',
  standards: ['ISO 2394:2015 General principles on reliability for structures', 'EN 1990:2002 Annex B and C (reliability index targets)', 'PIANC WG 12 (1992) / PIANC WG 28 (2003) partial-coefficient and reliability methods for breakwaters', 'CIRIA C683 Rock Manual (2007) §2.3 and §5.2.2.2 (coefficient scatter)', 'ISO 13822:2010 Assessment of existing structures', 'EurOtop (2018) §5.3 (model uncertainty)'],
  models: ['First-order reliability method (HL-RF)', 'Crude Monte Carlo simulation', 'Gumbel annual-maximum wave and surge climate', 'Van der Meer rock stability', 'EurOtop rubble-mound overtopping', 'Goda wave pressure with Coulomb sliding', 'Time-dependent deterioration and sea-level rise', 'Series system in time (independent years)'],
  limits: `Component reliability of one limit state with independent annual maxima; correlation between years (resistance common to all years), system effects between failure modes, statistical uncertainty of the extreme-value fits and human error are not modelled, so the lifetime probability is an upper estimate for resistance-dominated cases. Distributions and model-uncertainty factors are inputs to be fitted to site data. ${EDITION_NOTE}`,
  publishes: ['marine.betaAnnual', 'marine.pfLife', 'life.designYears'],
  searchTerms: 'structural reliability breakwater probabilistic design FORM service life sea level rise',
  inputs: [
    { key: 'ls', label: 'Limit state', type: 'select', default: 'armour', options: Object.entries(LS).map(([value, label]) => ({ value, label })), group: 'Limit state', role: 'criteria' },
    { key: 'life', label: 'Required service life', unit: 'years', type: 'number', default: 50, min: 5, max: 200, group: 'Limit state', role: 'criteria', bind: 'life.designYears' },
    { key: 'pfTarget', label: 'Tolerable annual probability of failure', unit: '-', type: 'number', default: 0.01, min: 1e-7, max: 0.3, group: 'Limit state', role: 'criteria', help: 'By consequence class: about 1e-2 for repairable damage to a breakwater, 1e-3 to 1e-4 for flood defence of occupied land, 1e-5 to 1e-6 where life is at risk (ISO 2394; EN 1990 Annex B). Confirm with the owner and regulator.' },
    { key: 'HsMean', label: 'Mean annual-maximum Hs offshore of the toe', unit: 'm', type: 'number', default: 3.6, min: 0.2, max: 15, group: 'Load climate', role: 'bc', help: 'Mean of the yearly largest significant wave height; the 100-year value follows from the Gumbel tail.' },
    { key: 'HsCov', label: 'Coefficient of variation of annual-maximum Hs', unit: '-', type: 'number', default: 0.18, min: 0.03, max: 0.6, group: 'Load climate', role: 'bc' },
    { key: 'Tm', label: 'Mean wave period in the storm', unit: 's', type: 'number', default: 9, min: 3, max: 20, group: 'Load climate', role: 'bc' },
    { key: 'surgeMean', label: 'Mean annual-maximum surge', unit: 'm', type: 'number', default: 0.6, min: 0, max: 5, group: 'Load climate', role: 'bc', bind: 'site.surge100' },
    { key: 'surgeStd', label: 'Standard deviation of annual-maximum surge', unit: 'm', type: 'number', default: 0.2, min: 0.01, max: 2, group: 'Load climate', role: 'bc' },
    { key: 'hToe', label: 'Toe depth at design high water (no surge)', unit: 'm', type: 'number', default: 8, min: 0.5, max: 60, group: 'Load climate', role: 'bc', bind: 'marine.hToe' },
    { key: 'gammaBr', label: 'Depth-limited Hs / depth at the toe', unit: '-', type: 'number', default: 0.55, min: 0.3, max: 0.8, group: 'Load climate', role: 'criteria', help: 'The toe wave is the smaller of the offshore value and this ratio times the local depth.' },
    ...slrInputs({ group: 'Climate change', slr: 0.6 }),
    { key: 'slrCov', label: 'Coefficient of variation of the sea-level-rise allowance', unit: '-', type: 'number', default: 0.4, min: 0.05, max: 1, group: 'Climate change', role: 'bc', help: 'The allowance above is the mean rise at the end of the service life; it is applied as a quadratic (accelerating) path in time.' },
    { key: 'Dn50', label: 'Armour Dn50 as placed', unit: 'm', type: 'number', default: 1.75, min: 0.1, max: 4, group: 'Armour', role: 'ic', bind: 'marine.armourDn50', design: true, showIf: { key: 'ls', in: ['armour'] } },
    { key: 'DnCov', label: 'Coefficient of variation of Dn50', unit: '-', type: 'number', default: 0.05, min: 0.01, max: 0.2, group: 'Armour', role: 'material', showIf: { key: 'ls', in: ['armour'] } },
    { key: 'cotA', label: 'Front slope cot α', unit: '-', type: 'number', default: 2, min: 1.5, max: 6, group: 'Armour', role: 'geometry', showIf: { key: 'ls', in: ['armour', 'overtopping'] } },
    { key: 'Sd', label: 'Damage level defining failure', unit: '-', type: 'number', default: 8, min: 2, max: 17, group: 'Armour', role: 'criteria', showIf: { key: 'ls', in: ['armour'] }, help: '2 for start of damage (serviceability); 8 for failure of a 1:2 slope (Rock Manual Table 5.23 — confirm).' },
    { key: 'P', label: 'Notional permeability', unit: '-', type: 'number', default: 0.4, min: 0.1, max: 0.6, group: 'Armour', role: 'material', showIf: { key: 'ls', in: ['armour'] } },
    { key: 'rhoS', label: 'Armour density', unit: 'kg/m³', type: 'number', default: 2650, min: 2000, max: 3200, group: 'Armour', role: 'material', showIf: { key: 'ls', in: ['armour'] } },
    { key: 'wear', label: 'Loss of armour size per year (abrasion, breakage)', unit: '%/yr', type: 'number', default: 0.1, min: 0, max: 2, group: 'Deterioration', role: 'ic', help: 'Deterioration rate from inspection records; applied to Dn50, to the caisson weight or to the generic resistance.' },
    { key: 'crest', label: 'Crest freeboard at design high water today', unit: 'm', type: 'number', default: 6, min: 0.3, max: 30, group: 'Crest', role: 'geometry', design: true, showIf: { key: 'ls', in: ['overtopping'] } },
    { key: 'gammaF', label: 'Roughness factor γf', unit: '-', type: 'number', default: 0.4, min: 0.3, max: 1, group: 'Crest', role: 'material', showIf: { key: 'ls', in: ['overtopping'] }, help: 'Recalled from EurOtop (2018) Table 6.2; confirm against the current edition.' },
    { key: 'Gc', label: 'Armour crest width', unit: 'm', type: 'number', default: 4, min: 0, max: 30, group: 'Crest', role: 'geometry', showIf: { key: 'ls', in: ['overtopping'] } },
    { key: 'qTol', label: 'Discharge defining failure', unit: 'l/s/m', type: 'number', default: 10, min: 0.01, max: 500, group: 'Crest', role: 'criteria', showIf: { key: 'ls', in: ['overtopping'] } },
    { key: 'settle', label: 'Crest settlement per year', unit: 'mm/yr', type: 'number', default: 3, min: 0, max: 50, group: 'Deterioration', role: 'ic', showIf: { key: 'ls', in: ['overtopping'] } },
    { key: 'B', label: 'Caisson width', unit: 'm', type: 'number', default: 15, min: 2, max: 60, group: 'Caisson', role: 'geometry', design: true, bind: 'marine.caissonWidth', showIf: { key: 'ls', in: ['sliding'] } },
    { key: 'hc', label: 'Caisson crest above design high water', unit: 'm', type: 'number', default: 4, min: 0, max: 20, group: 'Caisson', role: 'geometry', showIf: { key: 'ls', in: ['sliding'] } },
    { key: 'rhoC', label: 'Average density of the filled caisson', unit: 'kg/m³', type: 'number', default: 2100, min: 1400, max: 2600, group: 'Caisson', role: 'material', showIf: { key: 'ls', in: ['sliding'] } },
    { key: 'muMean', label: 'Mean base friction coefficient', unit: '-', type: 'number', default: 0.64, min: 0.3, max: 0.9, group: 'Caisson', role: 'bc', showIf: { key: 'ls', in: ['sliding'] }, help: 'Takayama & Ikeda (1992): mean 0.64, cov 0.16 for concrete on rubble.' },
    { key: 'muCov', label: 'Coefficient of variation of friction', unit: '-', type: 'number', default: 0.16, min: 0.03, max: 0.4, group: 'Caisson', role: 'bc', showIf: { key: 'ls', in: ['sliding'] } },
    { key: 'Rmean', label: 'Mean resistance', unit: '-', type: 'number', default: 10, min: 0.001, max: 1e9, group: 'Generic', role: 'criteria', showIf: { key: 'ls', in: ['generic'] } },
    { key: 'Rcov', label: 'Coefficient of variation of resistance', unit: '-', type: 'number', default: 0.1, min: 0.01, max: 0.6, group: 'Generic', role: 'criteria', showIf: { key: 'ls', in: ['generic'] } },
    { key: 'Smean', label: 'Mean annual-maximum load', unit: '-', type: 'number', default: 6, min: 0.001, max: 1e9, group: 'Generic', role: 'bc', showIf: { key: 'ls', in: ['generic'] } },
    { key: 'Scov', label: 'Coefficient of variation of load', unit: '-', type: 'number', default: 0.25, min: 0.01, max: 1, group: 'Generic', role: 'bc', showIf: { key: 'ls', in: ['generic'] } },
    { key: 'Sdist', label: 'Load distribution', type: 'select', default: 'gumbel', options: [{ value: 'gumbel', label: 'Gumbel (annual maxima)' }, { value: 'normal', label: 'Normal' }], group: 'Generic', role: 'bc', showIf: { key: 'ls', in: ['generic'] } },
    { key: 'mCov', label: 'Model uncertainty (coefficient of variation)', unit: '-', type: 'number', default: 0.065, min: 0, max: 0.6, group: 'Limit state', role: 'criteria', help: 'Scatter of the formula itself: about 6.5 % on the Van der Meer coefficient, a factor on EurOtop discharge exponents, 10–20 % on Goda forces.' },
    { key: 'nMC', label: 'Monte Carlo samples for the cross-check', type: 'int', default: 10000, min: 1000, max: 400000, group: 'Solver', role: 'control' },
    { key: 'unitCost', label: 'Indicative cost of a failure (repair + consequence)', unit: 'USD', type: 'number', default: 5e6, min: 0, max: 1e11, group: 'Risk', role: 'cost', help: 'Indicative figure used only to express annual risk; replace with a project estimate.' },
  ],
  compute(inp) {
    const ls = inp.ls, N = Math.round(inp.life);
    const delta = relativeDensity({ rhoS: inp.rhoS, rhoW: 1025 });
    const slrAt = (t) => inp.slr * (t / inp.life) ** 2; // accelerating path reaching the allowance at the end of life
    /** Random variables and limit-state function for structure age t [years]. */
    const problem = (t, slrEnd = inp.slr) => {
      const det = Math.max(0.05, 1 - inp.wear / 100 * t), slrMean = slrEnd * (t / inp.life) ** 2;
      const vars = { m: inp.mCov > 0 ? { dist: 'lognormal', mean: 1, cov: inp.mCov } : { dist: 'deterministic', value: 1 } };
      if (ls === 'generic') {
        vars.R = { dist: 'normal', mean: inp.Rmean * det, cov: inp.Rcov };
        vars.S = { dist: inp.Sdist, mean: inp.Smean, cov: inp.Scov };
        return { vars, g: (x) => x.m * x.R - x.S };
      }
      vars.Hs = { dist: 'gumbel', mean: inp.HsMean, cov: inp.HsCov };
      vars.surge = { dist: 'gumbel', mean: Math.max(inp.surgeMean, 1e-3), std: inp.surgeStd };
      if (slrMean > 1e-6) vars.slr = { dist: 'lognormal', mean: slrMean, cov: inp.slrCov };
      const wl = (x) => Math.max(x.surge, -0.5 * inp.hToe) + (x.slr || 0);
      // smooth minimum of the offshore and the depth-limited wave height (p-norm, p = 8): the transition into the surf
      // zone is gradual in nature and a smooth limit state keeps the FORM gradient search well-behaved
      const toeHs = (x) => { const a = Math.max(x.Hs, 0.05), b = Math.max(inp.gammaBr * (inp.hToe + wl(x)), 0.05); return Math.pow(a ** -8 + b ** -8, -1 / 8); };
      if (ls === 'armour') {
        vars.Dn = { dist: 'normal', mean: inp.Dn50 * det, cov: inp.DnCov };
        return { vars, g: (x) => { const H = toeHs(x); return x.m * vanDerMeerRock({ Hs: H, Tm: inp.Tm, cotAlpha: inp.cotA, P: inp.P, S: inp.Sd, N: 3000, delta, rhoS: inp.rhoS }).Ns * delta * Math.max(x.Dn, 0.01) - H; } };
      }
      if (ls === 'overtopping') {
        const settled = inp.settle / 1000 * t;
        // model factor acts on the discharge through a log-normal multiplier with exponent scaled to the formula scatter
        return { vars, g: (x) => { const H = toeHs(x), Rc = Math.max(0.02, inp.crest - settled - wl(x)); const q = overtoppingRubbleMound({ Hm0: H, Tm10: inp.Tm * 1.1, tanAlpha: 1 / inp.cotA, Rc, gammaF: inp.gammaF, Gc: inp.Gc }).q; return Math.log(inp.qTol / 1000) - Math.log(Math.max(q, 1e-30)) - 10 * Math.log(x.m); } };
      }
      vars.mu = { dist: 'lognormal', mean: inp.muMean, cov: inp.muCov };
      return { vars, g: (x) => {
        const w = wl(x), h = inp.hToe + w, H = toeHs(x), hc = Math.max(inp.hc - w, 0);
        const gp = godaPressure({ Hmax: Math.pow((1.8 * H) ** -8 + (0.78 * h) ** -8, -1 / 8), T: inp.Tm * 1.15, h, d: 0.8 * h, hPrime: 0.9 * h, hc, hb: 1.05 * h, B: inp.B });
        const W = (inp.rhoC * det * G * inp.B * (0.9 * inp.hToe + inp.hc) - 1030 * G * inp.B * 0.9 * h);
        return x.mu * (W - gp.U) - x.m * gp.P;
      } };
    };
    const solveAt = (t, slrEnd) => {
      const p = problem(t, slrEnd);
      let res; try { res = form(p.g, p.vars, { maxIter: 200 }); } catch (e) { throw new Error(`Reliability analysis at year ${t}: ${e.message}`); }
      const beta = Math.max(-8, Math.min(8, res.beta));
      return { beta, pf: normCdf(-beta), alpha: res.alpha, design: res.designPoint, importance: res.importance, converged: res.converged, p };
    };
    const years = Array.from(new Set(linspace(0, N, Math.min(N + 1, 6)).map((v) => Math.round(v))));
    const pts = years.map((t) => ({ t, ...solveAt(t) }));
    for (const q of pts) finite(q.beta, `Reliability index at year ${q.t}`);
    // lifetime: independent years, annual pf interpolated between the analysed ages (log-linear)
    const pfAt = (t) => { let i = 1; while (i < pts.length - 1 && pts[i].t < t) i++; const a = pts[i - 1], b = pts[i]; const w = b.t > a.t ? (t - a.t) / (b.t - a.t) : 0; return Math.exp(Math.log(Math.max(a.pf, 1e-300)) * (1 - w) + Math.log(Math.max(b.pf, 1e-300)) * w); };
    let surv = 1; const cum = [];
    for (let y = 1; y <= N; y++) { surv *= 1 - pfAt(y); cum.push(1 - surv); }
    const pfLife = 1 - surv, betaLife = -normInv(Math.min(Math.max(pfLife, 1e-15), 1 - 1e-15));
    const first = pts[0], last = pts[pts.length - 1];
    // service life: first year the annual pf exceeds the target
    let lifeYears = N + 0; let found = false;
    for (let y = 0; y <= N; y++) { if (pfAt(Math.max(y, 0.0001)) > inp.pfTarget) { lifeYears = y; found = true; break; } }
    const betaTarget = -normInv(inp.pfTarget);
    // Monte Carlo cross-check at the end of life
    const mc = run(() => monteCarlo(last.p.g, last.p.vars, { n: Math.round(inp.nMC), seed: 7 }), 'Monte Carlo cross-check');
    const imp = Object.entries(last.importance).sort((a, b) => b[1] - a[1]);
    const names = { Hs: 'Annual-maximum wave height', surge: 'Storm surge', slr: 'Sea-level rise', Dn: 'Armour size', mu: 'Base friction', m: 'Model uncertainty', R: 'Resistance', S: 'Load' };
    const scen = slrRows(inp).map((s) => { const q = solveAt(N, s.rise); return { name: s.name, rise: s.rise, beta: r(q.beta, 3), pf: r(q.pf, 3), ok: q.pf <= inp.pfTarget ? 'yes' : 'no' }; });
    const risk = last.pf * inp.unitCost, riskNow = first.pf * inp.unitCost;
    const ok = !found;
    // strengthening target: resistance scale that restores the target at end of life (FORM sensitivity, first order)
    const top = imp[0] ? imp[0][0] : 'm';
    const dBeta = Math.max(0, betaTarget - last.beta);
    const armourUp = ls === 'armour' ? inp.Dn50 * (Math.exp(dBeta * Math.hypot(inp.DnCov, inp.mCov, 0.12)) - 1) : 0;

    const advice = [
      ok
        ? eng('ok', `${LS[ls]}: annual failure probability rises from ${f(first.pf, 2)} today (β = ${f(first.beta)}) to ${f(last.pf, 2)} at year ${N} (β = ${f(last.beta)}), within the tolerable ${f(inp.pfTarget, 2)}; lifetime probability ${f(pfLife, 2)}.`,
          `The target is held for the whole ${N}-year life. ${names[top] || top} carries ${f((imp[0]?.[1] ?? 0) * 100, 2)} % of the uncertainty at the design point — spend investigation effort there before spending on the structure.`,
          'ISO 2394:2015 §8 and Annex G; EN 1990:2002 Annex C; PIANC WG 28 (2003)', 'A structure shown to hold its reliability target needs no precautionary strengthening — the material is saved, and the monitoring that confirms it is far lighter than rock or concrete.')
        : eng('fail', `${LS[ls]}: the tolerable annual failure probability of ${f(inp.pfTarget, 2)} (β = ${f(betaTarget)}) is exceeded from year ${lifeYears} — today ${f(first.pf, 2)}, at year ${N} ${f(last.pf, 2)} (β = ${f(last.beta)}); lifetime probability ${f(pfLife, 2)}.`,
          `${lifeYears === 0 ? 'The structure does not meet the target now.' : `Plan intervention before year ${lifeYears}.`} The index must rise by ${f(dBeta, 2)}: ${ls === 'armour' ? `about ${f(armourUp, 2)} m more Dn50 (${f(((1 + armourUp / inp.Dn50) ** 3 - 1) * 100, 2)} % more mass), a flatter slope or a berm` : ls === 'overtopping' ? 'a higher or wider crest, a rougher armour or a wave wall' : ls === 'sliding' ? 'a wider or heavier caisson, a friction mat or a rear rubble backing' : 'a higher resistance or a lower load'}; ${names[top] || top} dominates (${f((imp[0]?.[1] ?? 0) * 100, 2)} %).`,
          'ISO 2394:2015 §8; ISO 13822:2010 §7 (existing structures); PIANC WG 28 (2003)',
          ls === 'armour' ? `Upsizing the armour by ${f(((1 + armourUp / inp.Dn50) ** 3 - 1) * 100, 2)} % in mass thickens the layer by ${f(armourUp / inp.Dn50 * 100, 2)} %: per 1000 t of armour that is about ${f(armourUp / inp.Dn50 * 1000, 3)} t more rock (≈ ${f(armourUp / inp.Dn50 * 1000 * CO2_ROCK, 2)} t CO₂e, screening factor). Timing the upgrade to the trigger year avoids placing it decades early.` : `Strengthening adds concrete or rock (screening factors ${f(CO2_CONC)} t CO₂e/m³ and ${f(CO2_ROCK)} t CO₂e/t); deferring it to the trigger year identified here avoids building for a sea level that may arrive later than assumed.`),
      eng(scen[scen.length - 1].ok === 'no' ? 'warn' : 'info', `Sea-level rise drives the trend: at the end of life the annual failure probability ranges from ${f(scen[0].pf, 2)} (${scen[0].name}) to ${f(scen[scen.length - 1].pf, 2)} (${scen[scen.length - 1].name}).`,
        'Adopt an adaptive pathway: design for the central scenario, verify the high scenario can be met by a later upgrade, and review at each new projection or after any storm that causes damage.', SLR_BASIS),
      eng(Math.abs(mc.pf - last.pf) > 3 * Math.sqrt(Math.max(last.pf, 1e-9) / inp.nMC) + 0.3 * last.pf ? 'warn' : 'info', `Cross-check at year ${N}: FORM gives ${f(last.pf, 3)}, Monte Carlo ${f(mc.pf, 3)} (95 % interval ${f(mc.ciPf[0], 2)}–${f(mc.ciPf[1], 2)}, ${mc.nFail} failures in ${mc.n}).`,
        mc.nFail < 10 ? 'Too few failures were sampled to confirm the FORM result: increase the sample count or rely on FORM with a SORM or importance-sampling check.' : 'Where the two differ markedly the limit state is strongly curved (depth-limited waves cap the load): prefer the Monte Carlo value.', 'ISO 2394:2015 Annex E; Rackwitz & Fiessler (1978); Melchers & Beck (2018) ch. 3–4'),
      eco(ok ? 'info' : 'warn', `With the indicative failure cost of ${f(inp.unitCost, 3)} USD the annual risk grows from ${f(riskNow, 3)} to ${f(risk, 3)} USD per year over the life.`,
        'Justify strengthening when its annualised cost is below the risk it removes, and inspection when it can change that decision; the importance factors show which measurement has the highest value of information.', 'Risk-based optimisation per ISO 2394:2015 Annex G; value-of-information analysis; whole-life costing per ISO 15686-5',
        'Risk-based timing of upgrades spreads rock and concrete use over the life instead of front-loading it for the worst-case scenario.'),
    ];

    return {
      kpis: [
        { key: 'betaAnnual', label: 'Annual reliability index today', value: r(first.beta), unit: '-', explain: 'Distance of the limit state from the mean in standard-normal space for a one-year reference period at the present age. Higher is safer; 2.3 ≈ 1 in 100 per year, 3.1 ≈ 1 in 1000.' },
        { key: 'pfAnnual', label: 'Annual failure probability today', value: r(first.pf, 4), unit: '-', limit: inp.pfTarget, sense: 'max', explain: 'Probability that the limit state is exceeded in a year at the present age, from FORM.' },
        { key: 'betaEnd', label: 'Annual reliability index at the end of life', value: r(last.beta), unit: '-', explain: 'The same index after sea-level rise and deterioration over the service life. The drop from today is the price of ageing and climate change.' },
        { key: 'pfEnd', label: 'Annual failure probability at the end of life', value: r(last.pf, 4), unit: '-', limit: inp.pfTarget, sense: 'max', explain: 'Probability of failure in the final year of the service life; compared with the tolerable annual value.' },
        { key: 'pfLife', label: 'Lifetime failure probability', value: r(pfLife, 4), unit: '-', explain: 'Probability of at least one failure during the service life, treating years as independent. An upper estimate when resistance uncertainty dominates.' },
        { key: 'betaLife', label: 'Lifetime reliability index', value: r(betaLife), unit: '-', explain: 'Reliability index equivalent to the lifetime failure probability, for comparison with code targets stated for the design life (EN 1990: 3.8 for 50 years, ultimate limit states of buildings).' },
        { key: 'serviceLife', label: 'Years until the target is exceeded', value: lifeYears, unit: 'years', limit: N, sense: 'min', explain: 'Age at which the annual failure probability first exceeds the tolerable value — the reliability-based service life. Equal to the required life when never exceeded.' },
        { key: 'pfMC', label: 'Monte Carlo failure probability at the end of life', value: r(mc.pf, 4), unit: '-', explain: 'Independent estimate by direct sampling of the same limit state. Agreement with FORM confirms that the limit state is not strongly nonlinear.' },
        { key: 'betaTarget', label: 'Target annual reliability index', value: r(betaTarget), unit: '-', explain: 'Reliability index equivalent to the tolerable annual failure probability.' },
        { key: 'risk', label: 'Annual risk at the end of life (indicative)', value: r(risk), unit: 'USD/yr', explain: 'Failure probability times the indicative failure cost entered. A planning figure for comparing strengthening and monitoring options, not an estimate.' },
      ],
      checks: [
        check('Annual failure probability today', first.pf, inp.pfTarget, '-', 'ISO 2394:2015 Annex G; EN 1990 Annex B', 'Guards against operating a structure whose yearly chance of failing exceeds what the owner and society accept for its consequence class.'),
        check('Annual failure probability at the end of life', last.pf, inp.pfTarget, '-', 'ISO 2394:2015 §8; ISO 13822:2010', 'The same criterion with sea-level rise and deterioration included: it decides whether the structure reaches its intended life without upgrading.'),
      ],
      plots: [
        { type: 'line', title: 'Reliability index through the service life', x: { label: 'Age', unit: 'years' }, y: { label: 'Annual reliability index β', unit: '-' },
          series: [{ name: 'Annual β', x: pts.map((q) => q.t), y: pts.map((q) => r(q.beta, 3)) }, { name: 'Target', x: [0, N], y: [r(betaTarget, 3), r(betaTarget, 3)], dash: true }],
          explain: 'Annual reliability as the sea rises and the structure deteriorates. Where the curve falls through the target line is the reliability-based end of service life.' },
        { type: 'line', title: 'Cumulative probability of failure', x: { label: 'Age', unit: 'years' }, y: { label: 'Probability of at least one failure', unit: '-' },
          series: [{ name: 'Cumulative', x: cum.map((_, i) => i + 1), y: cum.map((v) => r(v, 4)) }],
          explain: 'Chance that the limit state has been exceeded at least once by a given age. Its slope steepens as the annual probability grows late in life.' },
        { type: 'bar', title: 'What drives the failure probability (end of life)', categories: imp.map(([k]) => names[k] || k), series: [{ name: 'Importance α²', y: imp.map(([, v]) => r(v, 3)) }],
          explain: 'Share of the uncertainty at the design point due to each variable (squared FORM direction cosines, summing to 1). The tallest bar is where better data or tighter control buys most reliability.' },
      ],
      tables: [
        { title: 'Design point at the end of life', columns: [{ key: 'v', label: 'Variable' }, { key: 'x', label: 'Value at the design point' }, { key: 'a', label: 'Sensitivity α', unit: '-' }, { key: 'i', label: 'Importance', unit: '%' }],
          rows: Object.keys(last.design).map((k) => ({ v: names[k] || k, x: r(last.design[k]), a: r(last.alpha[k], 3), i: r(last.importance[k] * 100, 3) })),
          explain: 'The most likely combination of values at failure. Positive α marks load-type variables, negative α resistance-type; the design-point values are what a partial-factor design should cover.' },
        { title: 'Sea-level-rise scenarios at the end of life', columns: [{ key: 'name', label: 'Scenario' }, { key: 'rise', label: 'Rise', unit: 'm' }, { key: 'beta', label: 'Annual β', unit: '-' }, { key: 'pf', label: 'Annual pf', unit: '-' }, { key: 'ok', label: 'Meets target' }],
          rows: scen, explain: 'End-of-life reliability with each scenario taken as the mean sea-level rise. The first scenario that fails the target marks the limit of the present design.' },
        { title: 'Reliability by age', columns: [{ key: 't', label: 'Age', unit: 'years' }, { key: 'slr', label: 'Mean rise', unit: 'm' }, { key: 'beta', label: 'β', unit: '-' }, { key: 'pf', label: 'Annual pf', unit: '-' }],
          rows: pts.map((q) => ({ t: q.t, slr: r(slrAt(q.t), 3), beta: r(q.beta, 3), pf: r(q.pf, 3) })), explain: 'Annual reliability at each analysed age with the mean sea-level rise reached by then; use it to set inspection and review dates.' },
      ],
      notes: [`${LS[ls]}; ${names[top] || top} dominates the uncertainty at the end of life. FORM ${pts.every((q) => q.converged) ? 'converged at every age' : 'did not fully converge at some ages — check with Monte Carlo'}.`,
        `Sea-level rise follows a quadratic path to ${f(inp.slr)} m at year ${N}; deterioration ${f(inp.wear)} % per year.`],
      advice,
      publish: { 'marine.betaAnnual': first.beta, 'marine.pfLife': pfLife, 'life.designYears': N },
    };
  },
  verify: [
    { name: 'Normal R − S: β = (μR − μS)/√(σR² + σS²)', source: 'Cornell (1969); Melchers & Beck (2018) Eq. 1.22. μR = 10, σR = 1, μS = 6, σS = 1.5 → β = 4/√3.25 = 2.2188, pf = Φ(−2.2188) = 0.01325; 50 independent years → 1 − (1 − 0.01325)^50 = 0.4867',
      inputs: { ls: 'generic', Rmean: 10, Rcov: 0.1, Smean: 6, Scov: 0.25, Sdist: 'normal', mCov: 0, wear: 0, slr: 0, life: 50 },
      expect: [{ kpi: 'betaAnnual', value: 2.2188, tol: 0.002 }, { kpi: 'pfAnnual', value: 0.01325, tol: 0.01 }, { kpi: 'pfLife', value: 0.4867, tol: 0.01 }, { kpi: 'pfMC', value: 0.01325, tol: 0.15 }] },
  ],
};
