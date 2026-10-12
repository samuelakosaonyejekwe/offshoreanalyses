// Mooring line fatigue: long-term Miner damage from a scatter of sea states with wave- and low-frequency tension
// components, using DNV-OS-E301 S–N curves or API RP 2SK T–N curves, with a rainflow cross-check.
import { mooringSnConstants, mooringNominalArea, mooringTnApi, chainBreakingLoad } from '../../engine/fatigue/sn.js';
import { bimodalDamageRate } from '../../engine/fatigue/spectral.js';
import { rainflowFourPoint, rangeHistogram, damageEquivalentLoad } from '../../engine/fatigue/rainflow.js';
import { minerSumCycles, fatigueLife } from '../../engine/fatigue/damage.js';
import { makeRng } from '../../engine/math/rng.js';
import { YEAR_S, need, run, round, fmt, adv, check, sev, steelCo2Text, vesselCo2Text, CONFIRM } from './_common.js';

const COMPONENTS = [
  { value: 'chain-studless', label: 'Studless chain', dnv: 'chain-studless', api: 'chain-studless' },
  { value: 'chain-stud', label: 'Stud-link chain', dnv: 'chain-stud', api: 'chain-stud' },
  { value: 'connecting-link', label: 'Connecting link (Kenter / D-shackle)', dnv: 'chain-studless', api: 'connecting-link' },
  { value: 'rope-stranded', label: 'Six/multi-strand wire rope', dnv: 'rope-stranded', api: 'rope-six-strand' },
  { value: 'rope-spiral', label: 'Spiral-strand wire rope', dnv: 'rope-spiral', api: 'rope-spiral-strand' },
];
const SCATTER = [
  { Hs: 1, Tp: 7, prob: 30, Tmean: 1320, sdWF: 3, TzWF: 6, sdLF: 2, TzLF: 100 },
  { Hs: 2, Tp: 8, prob: 28, Tmean: 1350, sdWF: 7, TzWF: 6.5, sdLF: 5, TzLF: 100 },
  { Hs: 3, Tp: 9, prob: 18, Tmean: 1400, sdWF: 11, TzWF: 7, sdLF: 8, TzLF: 100 },
  { Hs: 4, Tp: 10, prob: 11, Tmean: 1470, sdWF: 17, TzWF: 7.5, sdLF: 13, TzLF: 100 },
  { Hs: 5, Tp: 11, prob: 6, Tmean: 1560, sdWF: 23, TzWF: 8, sdLF: 19, TzLF: 100 },
  { Hs: 6, Tp: 12, prob: 4, Tmean: 1660, sdWF: 30, TzWF: 8.5, sdLF: 25, TzLF: 100 },
  { Hs: 8, Tp: 13, prob: 2.5, Tmean: 1880, sdWF: 44, TzWF: 9.5, sdLF: 38, TzLF: 100 },
  { Hs: 10, Tp: 14.5, prob: 0.5, Tmean: 2100, sdWF: 61, TzWF: 10.5, sdLF: 53, TzLF: 100 },
];

/** Narrow-band Gaussian record with standard deviation sd and mean period Tz (sum of sinusoids, seeded phases). */
function narrowBand(sd, Tz, t, rng, n = 24, band = 0.12) {
  const out = new Float64Array(t.length);
  if (!(sd > 0)) return out;
  const a = sd * Math.sqrt(2 / n), w0 = 2 * Math.PI / Tz;
  for (let k = 0; k < n; k++) {
    const w = w0 * (1 + band * (2 * (k + 0.5) / n - 1)), ph = 2 * Math.PI * rng.next();
    for (let i = 0; i < t.length; i++) out[i] += a * Math.cos(w * t[i] + ph);
  }
  return out;
}

export default {
  id: 's10-line-fatigue',
  suite: 10,
  title: 'Mooring line fatigue',
  summary: 'Sums fatigue damage of a chain or wire-rope segment over the long-term sea-state scatter and returns the fatigue life against the design life and safety factor.',
  covers: [6],
  fidelity: 'analytical',
  standards: ['DNV-OS-E301 (2018) Ch.2 Sec.2 [6] (fatigue limit state)', 'API RP 2SK (3rd ed., 2005) §6.4 and T–N curves', 'ISO 19901-7 (2013) §9', 'ASTM E1049-85 (2017) rainflow counting'],
  models: ['Chain and wire S–N curves (DNV-OS-E301) and T–N curves (API RP 2SK)', 'Palmgren–Miner linear damage', 'Narrow-band Rayleigh damage per sea state', 'Dual narrow-band combination of low- and wave-frequency tension (Jiao & Moan 1990)', 'Rainflow cycle counting of a synthesised bimodal tension history', 'Corrosion / wear reduction of the nominal cross-section'],
  limits: 'Tension-tension fatigue of the straight line only: out-of-plane bending of top chain at fairleads and chain stoppers, bending over sheaves, interlink wear and fibre-rope fatigue are not covered. Tension standard deviations per sea state are inputs from the coupled or tension analysis; the S–N and T–N constants are those recalled from the standards and the design fatigue factor is an overridable default to confirm against the governing edition.',
  publishes: ['life.fatigueYears'],
  searchTerms: 'mooring chain fatigue T-N curve S-N Miner dual narrow band',
  inputs: [
    { key: 'component', label: 'Component', type: 'select', default: 'chain-studless', options: COMPONENTS.map(({ value, label }) => ({ value, label })), group: 'Component', role: 'material',
      help: 'The fatigue-critical segment: usually the top chain or the chain near touchdown; connecting links have the lowest curve.' },
    { key: 'd', label: 'Nominal diameter', unit: 'mm', type: 'number', default: 95, min: 20, max: 250, group: 'Component', role: 'geometry', design: true },
    { key: 'MBL', label: 'Minimum breaking load', unit: 'kN', type: 'number', default: 9000, min: 100, max: 60000, group: 'Component', role: 'material', bind: 'mooring.MBL', help: 'Used as reference strength for wire rope (API) and for the mean-load ratio.' },
    { key: 'standard', label: 'Fatigue curve', type: 'select', default: 'dnv', options: [{ value: 'dnv', label: 'DNV-OS-E301 S–N (nominal stress)' }, { value: 'api', label: 'API RP 2SK T–N (tension / reference break strength)' }], group: 'Fatigue resistance', role: 'material',
      help: 'DNV: nc = aD·s^−m with s the nominal stress range (two bar areas for chain). API: N·R^M = K with R the tension range over the ORQ break strength of the same diameter (chain) or the rope MBL.' },
    { key: 'corrRate', label: 'Corrosion / wear rate on diameter', unit: 'mm/yr', type: 'number', default: 0.3, min: 0, max: 1.5, group: 'Fatigue resistance', role: 'material', uncertain: { dist: 'lognormal', cov: 0.3 },
      help: 'DNV-OS-E301 guidance: 0.2–0.4 mm/yr in the splash zone and at the seabed, 0.1–0.3 mid-catenary (higher in warm water). Site-specific.' },
    { key: 'corrFrac', label: 'Share of the life-time allowance used for fatigue', unit: '-', type: 'number', default: 0.5, min: 0, max: 1, group: 'Fatigue resistance', role: 'criteria', help: 'DNV-OS-E301: half of the corrosion allowance is deducted from the diameter in the fatigue calculation.' },
    { key: 'scatter', label: 'Sea-state scatter with tension response', type: 'table', group: 'Fatigue loading', role: 'load',
      columns: [{ key: 'Hs', label: 'Hs', unit: 'm', type: 'number' }, { key: 'Tp', label: 'Tp', unit: 's', type: 'number' }, { key: 'prob', label: 'Probability', unit: '%', type: 'number' }, { key: 'Tmean', label: 'Mean tension', unit: 'kN', type: 'number' },
        { key: 'sdWF', label: 'WF tension std. dev.', unit: 'kN', type: 'number' }, { key: 'TzWF', label: 'WF period', unit: 's', type: 'number' }, { key: 'sdLF', label: 'LF tension std. dev.', unit: 'kN', type: 'number' }, { key: 'TzLF', label: 'LF period', unit: 's', type: 'number' }],
      default: SCATTER,
      help: 'One row per sea-state (and direction) bin of the wave scatter diagram with the tension response of the line in that bin: mean, wave-frequency and low-frequency standard deviations and their zero-crossing periods, from the coupled analysis. Probabilities are normalised to 100 %.' },
    { key: 'method', label: 'Combination of LF and WF damage', type: 'select', default: 'dual-narrow-band', options: [{ value: 'dual-narrow-band', label: 'Dual narrow-band (Jiao–Moan)' }, { value: 'combined-spectrum', label: 'Combined spectrum (conservative)' }, { value: 'simple-summation', label: 'Simple summation (can be non-conservative)' }], group: 'Fatigue loading', role: 'control',
      help: 'DNV-OS-E301 allows simple summation only when one component clearly dominates; otherwise combined spectrum or dual narrow-band.' },
    { key: 'designLife', label: 'Design service life', unit: 'years', type: 'number', default: 25, min: 1, max: 60, group: 'Acceptance criteria', role: 'criteria', bind: 'life.designYears' },
    { key: 'dff', label: 'Design fatigue (safety) factor', unit: '-', type: 'number', default: 5, min: 1, max: 15, group: 'Acceptance criteria', role: 'criteria',
      help: `DNV-OS-E301 (2018): γF = 5 for lines not regularly inspected ashore, rising to 8 when adjacent lines' damages are close; API RP 2SK: 3. ${CONFIRM}` },
    { key: 'D0', label: 'Fatigue damage already accumulated', unit: '-', type: 'number', default: 0, min: 0, max: 5, group: 'Initial condition', role: 'ic',
      help: 'Characteristic Miner sum consumed so far (previous service, tow, or a re-used line). Initial condition for remaining-life assessment.' },
    { key: 'yearsInService', label: 'Years already in service', unit: 'years', type: 'number', default: 0, min: 0, max: 60, group: 'Initial condition', role: 'ic' },
    { key: 'rfSeconds', label: 'Rainflow check: record length', unit: 's', type: 'number', default: 3600, min: 600, max: 21600, group: 'Numerical', role: 'control', help: 'Length of the synthesised tension record of the most damaging bin used to cross-check the spectral damage by rainflow counting.' },
    { key: 'seed', label: 'Rainflow check: random seed', type: 'int', default: 7, min: 1, max: 99999, group: 'Numerical', role: 'control' },
  ],
  compute(inp) {
    const comp = COMPONENTS.find((c) => c.value === inp.component);
    const d = inp.d / 1e3;
    const allowance = inp.corrFrac * inp.corrRate / 1e3 * inp.designLife;          // m on diameter
    need(allowance < 0.5 * d, `The corrosion allowance used for fatigue (${fmt(allowance * 1e3)} mm) removes more than half the diameter: reduce the rate or increase the diameter.`);
    const rows = (inp.scatter || []).filter((r) => r.prob > 0);
    need(rows.length > 0, 'The scatter table has no row with a positive probability.');
    const pSum = rows.reduce((a, r) => a + r.prob, 0);
    // S–N constants in a neutral form; "stress" unit is Pa for DNV and N (tension) for API
    let sn, toS, label, refStrength;
    if (inp.standard === 'dnv') {
      sn = mooringSnConstants(comp.dnv);
      const area = run(() => mooringNominalArea({ type: comp.dnv, diameter: d, corrosionAllowance: allowance }), 'Nominal area').area;
      toS = (T) => T / area; label = `DNV-OS-E301 ${comp.dnv}, nominal area ${fmt(area * 1e6)} mm²`; refStrength = inp.MBL * 1e3;
    } else {
      const isChain = comp.api.startsWith('chain') || comp.api === 'connecting-link';
      refStrength = isChain ? chainBreakingLoad({ diameter: d - allowance, grade: 'ORQ' }) : inp.MBL * 1e3 * ((d - allowance) / d) ** 2;
      const meanAll = rows.reduce((a, r) => a + r.prob * (r.Tmean || 0), 0) / pSum * 1e3;
      const Lm = Math.min(Math.max(meanAll / refStrength, 0), 0.9);
      const { M, K } = mooringTnApi(comp.api, Lm);
      sn = { m1: M, a1: K * refStrength ** M, m2: null, a2: null, sKnee: 0 };
      toS = (T) => T; label = `API RP 2SK ${comp.api}, reference strength ${fmt(refStrength / 1e3)} kN${isChain ? ' (ORQ)' : ''}, mean-load ratio ${fmt(Lm)}`;
    }
    const bins = rows.map((r, i) => {
      need(r.TzWF > 0 && r.TzLF > 0, `Scatter row ${i + 1}: periods must be positive.`);
      need(r.sdWF >= 0 && r.sdLF >= 0, `Scatter row ${i + 1}: standard deviations cannot be negative.`);
      const bm = bimodalDamageRate({ sigmaLF: toS(r.sdLF * 1e3), nuLF: 1 / r.TzLF, sigmaWF: toS(r.sdWF * 1e3), nuWF: 1 / r.TzWF, sn, method: inp.method });
      const p = r.prob / pSum;
      return { ...r, p, rate: bm.damageRate, annual: p * bm.damageRate * YEAR_S, wfOnly: p * bm.rateWF * YEAR_S, lfOnly: p * bm.rateLF * YEAR_S, rho: bm.rho };
    });
    const annual = bins.reduce((a, r) => a + r.annual, 0);
    need(annual > 0, 'All tension standard deviations are zero: there is no cyclic loading to accumulate damage.');
    const remainingDesign = Math.max(inp.designLife - inp.yearsInService, 0);
    const life = fatigueLife({ damageRate: annual, designLife: remainingDesign, dff: inp.dff, accumulatedDamage: inp.D0 });
    const dom = bins.reduce((a, r) => (r.annual > a.annual ? r : a), bins[0]);
    // --- rainflow cross-check on the most damaging bin
    const dtR = Math.min(dom.TzWF / 16, 0.6), nR = Math.min(Math.round(inp.rfSeconds / dtR), 60000);
    const t = Float64Array.from({ length: nR }, (_, i) => i * dtR);
    const rng = makeRng(inp.seed);
    const wf = narrowBand(dom.sdWF * 1e3, dom.TzWF, t, rng), lf = narrowBand(dom.sdLF * 1e3, dom.TzLF, t, rng, 12, 0.1);
    const hist = Float64Array.from(wf, (v, i) => dom.Tmean * 1e3 + v + lf[i]);
    const rf = rainflowFourPoint(Array.from(hist, toS), { residue: 'repeat' });
    const rfRate = minerSumCycles(rf.cycles, sn) / (nR * dtR);
    const rfRatio = dom.rate > 0 ? rfRate / dom.rate : 1;
    const hg = rangeHistogram(rf.cycles.map((c) => ({ ...c, range: c.range / toS(1e3) })), { bins: 24 });   // ranges back to kN
    const nCyc = rf.cycles.reduce((a, c) => a + c.count, 0);
    const del = damageEquivalentLoad({ cycles: rf.cycles.map((c) => ({ ...c, range: c.range / toS(1e3) })), m: sn.m1, nRef: nCyc });
    const uc = life.utilisation ?? (inp.dff * inp.D0);
    // diameter that would bring the utilisation to 1 (damage ∝ d^(−2m) for a given tension spectrum)
    const dReq = inp.d * uc ** (1 / (2 * sn.m1));
    const massPerM = 19.9e3 * d * d;                                              // kg/m, studless chain estimator for the quantity note
    const checks = [
      check('Fatigue damage over the design life (with safety factor)', uc, 1, '-', inp.standard === 'dnv' ? 'DNV-OS-E301 Ch.2 Sec.2 [6.1]: 1 − dc·γF ≥ 0' : 'API RP 2SK (3rd ed.) §6.4: fatigue life ≥ 3 × service life',
        'Guards against fatigue cracking and parting of the line within its service life. The Miner sum over the remaining design life, times the safety factor, plus damage already consumed must stay below 1.'),
    ];
    const order = bins.map((r, i) => i);
    return {
      kpis: [
        { key: 'annualDamage', label: 'Annual fatigue damage', value: round(annual, 5), unit: '1/yr', explain: 'Characteristic Miner damage accumulated per year, without the safety factor. Its reciprocal is the characteristic fatigue life.' },
        { key: 'life', label: 'Characteristic fatigue life', value: round(Math.min(life.characteristicLife, 1e6), 5), unit: 'years', explain: 'Years to a Miner sum of 1 from new, without safety factor.' },
        { key: 'lifeFactored', label: 'Fatigue life with safety factor', value: round(Math.min(life.designLife, 1e6), 5), unit: 'years', limit: inp.designLife, sense: 'min', explain: 'Characteristic life divided by the design fatigue factor; must be at least the design service life.' },
        { key: 'remaining', label: 'Predicted remaining service life', value: round(Math.min(life.remainingLife, 1e6), 5), unit: 'years', explain: 'Years the line can still serve from now, allowing for the safety factor and the damage already accumulated.' },
        { key: 'uc', label: 'Fatigue utilisation', value: round(uc, 4), unit: '-', limit: 1, sense: 'max', explain: 'Safety factor × (accumulated + future damage over the remaining design life). Above 1 the fatigue limit state is not met.' },
        { key: 'domShare', label: 'Share of damage from the worst bin', value: round(100 * dom.annual / annual, 4), unit: '%', explain: 'How concentrated the damage is in one sea state. Mooring fatigue is usually spread over moderate, frequent seas rather than the extreme storm.' },
        { key: 'lfShare', label: 'Low-frequency share of damage (if acting alone)', value: round(100 * bins.reduce((a, r) => a + r.lfOnly, 0) / Math.max(bins.reduce((a, r) => a + r.lfOnly + r.wfOnly, 0), 1e-300), 4), unit: '%', explain: 'Relative weight of slow-drift cycling against wave-frequency cycling in the damage; decides whether stiffness/damping or line dynamics is the lever.' },
        { key: 'rfRatio', label: 'Rainflow / spectral damage, worst bin', value: round(rfRatio, 4), unit: '-', explain: 'Miner damage from rainflow counting of a synthesised record divided by the spectral estimate for the same bin. Values near 1 (0.7–1.3 with sampling scatter) support the spectral method.' },
        { key: 'del', label: 'Damage-equivalent tension range, worst bin', value: round(del, 5), unit: 'kN', explain: 'Constant tension range that would do the same damage per cycle as the counted spectrum in the most damaging sea state.' },
        { key: 'dRequired', label: 'Diameter for utilisation 1.0', value: round(dReq, 4), unit: 'mm', limit: inp.d, sense: 'max', explain: 'Diameter at which the same tension spectrum gives a utilisation of exactly 1 (damage scales with d^−2m). Indicative, since a different diameter also changes the tensions slightly.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Annual damage by sea-state bin', categories: bins.map((r) => `Hs ${r.Hs} / Tp ${r.Tp}`), series: [{ name: 'Annual damage', y: bins.map((r) => round(r.annual, 4)) }],
          explain: 'Contribution of each sea state to the yearly Miner sum (probability × damage rate). The tallest bars are where reducing tension variation pays off.' },
        { type: 'bar', title: 'Rainflow tension-range histogram, worst bin', categories: hg.centres ? hg.centres.map((v) => String(round(v, 3))) : hg.edges.slice(0, -1).map((v, i) => String(round(0.5 * (v + hg.edges[i + 1]), 3))), series: [{ name: 'Cycles', y: Array.from(hg.counts, (v) => round(v, 4)) }],
          explain: 'Counted tension ranges (kN) of the synthesised record for the most damaging bin. Many small wave-frequency cycles ride on fewer, larger low-frequency ones; the large-range tail carries most of the damage.' },
        { type: 'line', title: 'Synthesised tension record, worst bin (first part)', x: { label: 'Time', unit: 's' }, y: { label: 'Tension', unit: 'kN' },
          series: [{ name: 'Tension', x: Array.from(t.subarray(0, Math.min(nR, 1200)), (v) => round(v)), y: Array.from(hist.subarray(0, Math.min(nR, 1200)), (v) => round(v / 1e3)) }],
          explain: 'Bimodal tension history built from the bin\'s standard deviations and periods: the slow undulation is low-frequency drift, the ripple is wave-frequency motion.' },
        { type: 'line', title: 'Damage accumulation over service life', x: { label: 'Years from now', unit: 'years' }, y: { label: 'Miner sum × safety factor', unit: '-' },
          series: [{ name: 'Factored damage', x: [0, Math.max(remainingDesign, 1) * 1.5], y: [inp.dff * inp.D0, inp.dff * (inp.D0 + annual * Math.max(remainingDesign, 1) * 1.5)] }, { name: 'Limit', x: [0, Math.max(remainingDesign, 1) * 1.5], y: [1, 1], dash: true },
            { name: 'End of design life', x: [remainingDesign, remainingDesign], y: [0, 1], dash: true }],
          explain: 'Factored Miner sum growing linearly with time. The line must reach the end of design life below 1; where it crosses 1 is the predicted remaining life.' },
      ],
      tables: [
        { title: 'Damage by sea-state bin', columns: [{ key: 'Hs', label: 'Hs', unit: 'm' }, { key: 'Tp', label: 'Tp', unit: 's' }, { key: 'p', label: 'Probability', unit: '%' }, { key: 'sdWF', label: 'WF std. dev.', unit: 'kN' }, { key: 'sdLF', label: 'LF std. dev.', unit: 'kN' }, { key: 'rho', label: 'Combination factor ρ' }, { key: 'annual', label: 'Annual damage' }, { key: 'share', label: 'Share', unit: '%' }],
          rows: order.map((i) => bins[i]).map((r) => ({ Hs: r.Hs, Tp: r.Tp, p: round(100 * r.p, 3), sdWF: r.sdWF, sdLF: r.sdLF, rho: round(r.rho, 3), annual: Number(r.annual.toPrecision(3)), share: round(100 * r.annual / annual, 3) })),
          explain: 'Per-bin damage with the dual narrow-band factor ρ applied to the combined-spectrum damage. Use the shares to decide which sea states need a refined dynamic analysis.' },
      ],
      notes: [`${label}; slope m = ${sn.m1}. Diameter reduced by ${fmt(allowance * 1e3)} mm for fatigue (${inp.corrFrac} × ${inp.corrRate} mm/yr × ${inp.designLife} yr).`,
        `Worst bin Hs ${dom.Hs} m carries ${fmt(100 * dom.annual / annual)} % of the damage; rainflow counting of a ${fmt(nR * dtR)} s record gives ${fmt(rfRatio)} × the spectral damage.`,
        'Out-of-plane bending at the fairlead is not included and often governs top-chain fatigue.'],
      advice: [
        adv('engineering', sev(uc), `Characteristic life ${fmt(life.characteristicLife)} years; with a safety factor of ${inp.dff} the usable life is ${fmt(life.designLife)} years against ${inp.designLife} required (utilisation ${fmt(uc)}).`,
          uc > 1 ? `Increase the diameter to about ${fmt(dReq)} mm (damage ∝ d^−${2 * sn.m1}), reduce the tension variation (softer system, more damping), or plan replacement of this segment after ${fmt(life.remainingLife)} years.`
            : `Remaining life is ${fmt(life.remainingLife)} years. Check out-of-plane bending of the top chain separately and confirm the corrosion rate at the first in-service inspection.`,
          'DNV-OS-E301 (2018) Ch.2 Sec.2 [6]; API RP 2SK (3rd ed., 2005) §6.4; ISO 19901-7 (2013) §9. Confirm the fatigue safety factor against the governing edition.'),
        adv('economic', uc > 1 ? 'warn' : 'info', uc > 1 ? `The segment would need replacing after about ${fmt(life.remainingLife)} years of a ${inp.designLife}-year life.` : `Fatigue leaves a reserve factor of ${fmt(1 / Math.max(uc, 1e-6))} on the required safety factor.`,
          uc > 1 ? `Compare one mid-life change-out (about 2–4 anchor-handling vessel days per line plus the segment) with ${fmt(100 * ((dReq / inp.d) ** 2 - 1))} % more chain mass up front; for top chain a change-out from the platform is often the cheaper option.` : uc < 0.5 ? `Inspection of this segment can follow the standard 5-yearly class cycle; a diameter of about ${fmt(dReq)} mm would still meet fatigue (strength and corrosion allowance permitting).` : 'Keep the diameter; schedule in-service measurement of the worst-loaded links at the class renewal survey.',
          'Risk-based inspection planning per DNV-RP-G101 / API RP 2I (in-service inspection of mooring hardware); whole-life costing per ISO 15686-5.',
          `${steelCo2Text(massPerM * 100 / 1e3)} (per 100 m of this chain). ${uc > 1 ? vesselCo2Text(3) + ' An early change-out repeats the segment\'s footprint; a slightly larger diameter avoids it.' : 'Meeting the life without replacement avoids a second set of steel and an offshore campaign.'}`),
      ],
      publish: { 'life.fatigueYears': Math.min(life.designLife, 1e6) },
      quantities: { steel_t: round(massPerM * 100 / 1e3) },
    };
  },
  verify: [
    { // One sea state, wave-frequency only: σT = 100 kN on 100 mm studless chain, no corrosion. Area = 2·π·0.1²/4 = 0.015708 m²,
      // σ = 6.366 MPa. Narrow-band: D/yr = ν·T·(2√2σ)^m·Γ(1 + m/2)/aD = (31 557 600/10)·(18.006)³·1.3293/6e10 = 0.4082.
      name: 'Narrow-band damage, DNV studless chain S–N curve',
      source: 'DNV-OS-E301 Ch.2 Sec.2 [6.3] closed-form narrow-band damage with aD = 6.0e10, m = 3',
      inputs: { component: 'chain-studless', d: 100, standard: 'dnv', corrRate: 0, method: 'dual-narrow-band', scatter: [{ Hs: 3, Tp: 12, prob: 100, Tmean: 1000, sdWF: 100, TzWF: 10, sdLF: 0, TzLF: 100 }] },
      expect: [{ kpi: 'annualDamage', value: 0.40819, tol: 0.002 }, { kpi: 'life', value: 2.4498, tol: 0.002 }] },
    { // Same loading with the API T–N curve: N·R³ = 316, R = ΔT / ORQ break strength; ORQ(100 mm) = 0.0211·100²·(44 − 8) = 7596 kN.
      // D/yr = 3 155 760·(2√2·100/7596)³·1.3293/316 = 0.6854.
      name: 'Narrow-band damage, API studless chain T–N curve',
      source: 'API RP 2SK (3rd ed.) T–N curve for studless chain (M = 3.0, K = 316), hand calculation',
      inputs: { component: 'chain-studless', d: 100, standard: 'api', corrRate: 0, method: 'dual-narrow-band', scatter: [{ Hs: 3, Tp: 12, prob: 100, Tmean: 1000, sdWF: 100, TzWF: 10, sdLF: 0, TzLF: 100 }] },
      expect: [{ kpi: 'annualDamage', value: 0.68538, tol: 0.002 }] },
  ],
};
