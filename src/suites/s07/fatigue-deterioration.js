// Fatigue and deterioration of a marine structure over its service life: S–N fatigue of a welded steel detail
// under the long-term wave stress distribution, crack growth, steel corrosion wastage, chloride-induced
// corrosion initiation and abrasion of concrete, and a condition-state projection.
import { snCurve, snConstants, SN_CURVE_NAMES, weibullClosedFormDamage, weibullScaleFromMax } from '../../engine/fatigue/sn.js';
import { fatigueLife } from '../../engine/fatigue/damage.js';
import { corrosionWastage, remainingWall, markovDeterioration } from '../../engine/fatigue/degradation.js';
import { parisLawBs7910, parisLifeClosedForm } from '../../engine/fatigue/fracture.js';
import { erfc, erfinv, gamma as gammaFn } from '../../engine/math/special.js';
import { YEAR_S, CO2_STEEL, CO2_CONC, need, run, finite, r, f, linspace, check, eng, eco } from './_common.js';

export default {
  id: 's07-fatigue-deterioration', suite: 7, title: 'Fatigue, corrosion and remaining life',
  summary: 'Fatigue life of a welded steel detail, crack growth, steel wastage, chloride and abrasion attack on concrete — and which of them ends the service life first.',
  covers: [9], fidelity: 'analytical',
  standards: ['DNV-RP-C203 (2021) §2.4, §5', 'BS 7910:2019 clause 8', 'fib Model Code for Service Life Design (Bulletin 34, 2006) / fib MC 2010 §7.8', 'EN 1992-1-1 §4 and EN 206 exposure classes XS1–XS3', 'BS 6349-1-4:2021 (materials, corrosion rates)', 'ISO 12944 / NORSOK M-501 (coating life)'],
  models: ['S–N fatigue with Miner summation', 'Weibull long-term stress-range distribution', 'Paris-law fracture mechanics', 'Linear and Guedes Soares–Garbatov corrosion wastage', 'Fick\'s second law chloride ingress', 'Linear abrasion', 'Markov condition-state deterioration'],
  limits: 'Nominal-stress S–N fatigue of one detail with a two-parameter Weibull stress spectrum; the stress range itself must come from a structural analysis. Chloride ingress uses a constant apparent diffusion coefficient with an optional ageing exponent and gives initiation only (no propagation or cracking). Corrosion rates and transition probabilities are project inputs, not code values.',
  publishes: ['life.fatigueYears', 'marine.remainingLife'],
  searchTerms: 'fatigue S-N marine structure corrosion chloride ingress service life remaining life',
  inputs: [
    { key: 'age', label: 'Present age of the structure', unit: 'years', type: 'number', default: 15, min: 0, max: 150, group: 'Initial condition', role: 'ic', help: 'Damage, wastage and chloride penetration accumulated to date are computed from this age.' },
    { key: 'life', label: 'Required service life from new', unit: 'years', type: 'number', default: 50, min: 5, max: 200, group: 'Initial condition', role: 'criteria', bind: 'life.designYears' },
    { key: 'sMax', label: 'Largest stress range in the service life', unit: 'MPa', type: 'number', default: 100, min: 1, max: 800, group: 'Fatigue loading', role: 'bc', design: true, uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Nominal stress range at the detail exceeded once in the service life, from the wave-load and structural analysis.' },
    { key: 'hW', label: 'Weibull shape of the stress-range distribution', unit: '-', type: 'number', default: 0.9, min: 0.4, max: 2, group: 'Fatigue loading', role: 'bc', uncertain: { dist: 'normal', cov: 0.08 }, help: 'About 0.8–1.0 for wave loading on fixed structures; higher for dynamically sensitive ones.' },
    { key: 'Tz', label: 'Mean stress-cycle period', unit: 's', type: 'number', default: 6, min: 0.5, max: 30, group: 'Fatigue loading', role: 'bc', help: 'Mean zero-crossing period of the waves (or traffic cycle interval).' },
    { key: 'curve', label: 'S–N curve', type: 'select', default: 'F', options: SN_CURVE_NAMES.map((n) => ({ value: n, label: n })), group: 'Steel detail', role: 'material', help: 'DNV-RP-C203 detail category: D butt weld, F fillet/attachment, W3 weld root.' },
    { key: 'env', label: 'Environment', type: 'select', default: 'seawater-cp', options: [{ value: 'air', label: 'Air' }, { value: 'seawater-cp', label: 'Sea water with cathodic protection' }, { value: 'free-corrosion', label: 'Free corrosion' }], group: 'Steel detail', role: 'bc' },
    { key: 'steel', label: 'Steel grade', type: 'select', source: 'materials:steel', default: 's355', group: 'Steel detail', role: 'material' },
    { key: 't', label: 'Plate / wall thickness', unit: 'mm', type: 'number', default: 30, min: 4, max: 150, group: 'Steel detail', role: 'geometry', design: true },
    { key: 'scf', label: 'Stress concentration factor', unit: '-', type: 'number', default: 1.2, min: 1, max: 10, group: 'Steel detail', role: 'geometry' },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 3, min: 1, max: 10, group: 'Criteria', role: 'criteria', help: '1 to 10 depending on consequence and access for inspection (DNV-ST-0126 / ISO 19902).' },
    { key: 'a0', label: 'Initial flaw depth', unit: 'mm', type: 'number', default: 0.5, min: 0.05, max: 20, group: 'Steel detail', role: 'ic', help: 'Detectable or assumed weld-toe defect at the last inspection.' },
    { key: 'Y', label: 'Crack geometry factor Y', unit: '-', type: 'number', default: 1.12, min: 0.5, max: 3, group: 'Steel detail', role: 'geometry' },
    { key: 'tMin', label: 'Minimum acceptable steel thickness', unit: 'mm', type: 'number', default: 20, min: 1, max: 150, group: 'Corrosion of steel', role: 'criteria' },
    { key: 'rate', label: 'Corrosion rate after coating breakdown', unit: 'mm/yr', type: 'number', default: 0.3, min: 0.005, max: 2, group: 'Corrosion of steel', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.4 }, help: 'Indicative uncoated carbon steel: 0.1–0.2 immersed, 0.3–0.5 splash zone, more with accelerated low-water corrosion. Use site measurements.' },
    { key: 'coat', label: 'Coating / protection life', unit: 'years', type: 'number', default: 20, min: 0, max: 60, group: 'Corrosion of steel', role: 'ic' },
    { key: 'conc', label: 'Concrete', type: 'select', source: 'materials:concrete', default: 'c40', group: 'Concrete durability', role: 'material' },
    { key: 'cover', label: 'Cover to reinforcement', unit: 'mm', type: 'number', default: 60, min: 15, max: 150, group: 'Concrete durability', role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.12 } },
    { key: 'Da', label: 'Apparent chloride diffusion coefficient at 28 days', unit: '1e-12 m²/s', type: 'number', default: 3, min: 0.05, max: 50, group: 'Concrete durability', role: 'material', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'From migration or ponding tests: about 1–3 for slag or fly-ash blends, 5–15 for plain Portland cement.' },
    { key: 'ageing', label: 'Ageing exponent', unit: '-', type: 'number', default: 0.4, min: 0, max: 0.7, group: 'Concrete durability', role: 'material', help: 'D(t) = D28 (28 d / t)^α: about 0.3 Portland, 0.45 slag, 0.6 fly ash (fib Bulletin 34). 0 keeps D constant.' },
    { key: 'Cs', label: 'Surface chloride content', unit: '% binder', type: 'number', default: 3.5, min: 0.3, max: 8, group: 'Concrete durability', role: 'bc', help: 'Exposure boundary condition: about 3–5 % in the splash and tidal zone (XS3), 1.5–3 % in spray (XS1).' },
    { key: 'Ccr', label: 'Critical chloride content', unit: '% binder', type: 'number', default: 0.6, min: 0.1, max: 2, group: 'Concrete durability', role: 'criteria' },
    { key: 'abr', label: 'Abrasion rate', unit: 'mm/yr', type: 'number', default: 0.2, min: 0, max: 10, group: 'Concrete durability', role: 'bc', help: 'Loss of concrete surface from sand, shingle or ice; shingle beaches can exceed 1 mm/yr.' },
    { key: 'pStay', label: 'Annual probability of staying in a condition state', unit: '-', type: 'number', default: 0.96, min: 0.5, max: 0.999, group: 'Condition model', role: 'ic', help: 'From inspection records: the share of elements that keep their condition grade from one year to the next.' },
  ],
  compute(inp, ctx) {
    const steel = ctx.material(inp.steel), conc = ctx.material(inp.conc);
    need(steel && conc, 'The selected steel or concrete is not in the materials library.');
    need(inp.tMin < inp.t, 'The minimum acceptable thickness must be less than the as-built thickness.');
    need(inp.Ccr < inp.Cs, 'The critical chloride content must be below the surface content, otherwise corrosion never initiates.');
    need(inp.a0 < inp.t / 2, 'The initial flaw must be less than half the plate thickness.');
    // --- S–N fatigue: Weibull closed form with n0 cycles in the service life
    const n0 = inp.life * YEAR_S / inp.Tz, nYear = YEAR_S / inp.Tz;
    const sn = run(() => snConstants(snCurve(inp.curve, inp.env), { thickness: inp.t / 1000, scf: inp.scf }), 'S–N curve');
    const q = weibullScaleFromMax(inp.sMax * 1e6, inp.hW, n0);
    const Dlife = finite(weibullClosedFormDamage({ sn, scale: q, shape: inp.hW, cycles: n0 }), 'Fatigue damage');
    const rateD = Dlife / inp.life;
    const fl = fatigueLife({ damageRate: Math.max(rateD, 1e-30), designLife: inp.life, dff: inp.dff, accumulatedDamage: 0 });
    const fatLife = Math.min(fl.designLife, 9999), fatRemain = Math.max(0, fatLife - inp.age), ucFat = inp.dff * Dlife;
    // --- crack growth from the initial flaw to half thickness under the damage-equivalent constant range
    const law = parisLawBs7910({ environment: inp.env === 'air' ? 'air' : 'marine', curve: 'simplified' }).stages[0];
    const sEq = q * inp.scf * Math.pow(gammaFn(1 + law.m / inp.hW), 1 / law.m);
    const Ncrack = parisLifeClosedForm({ a0: inp.a0 / 1000, af: inp.t / 2000, stressRange: sEq, Y: inp.Y, A: law.A, m: law.m });
    const crackLife = Math.min(Ncrack / nYear, 9999);
    // --- steel wastage
    const corr = { model: 'linear', coatingLife: inp.coat, rate: inp.rate / 1000 };
    const wall = remainingWall({ thickness: inp.t / 1000, minThickness: inp.tMin / 1000, age: inp.age, surfaces: 1, corrosion: corr, horizon: 500 });
    const steelEnd = Number.isFinite(wall.endOfLifeAge) ? wall.endOfLifeAge : 9999;
    // --- chloride initiation: C(x,t) = Cs erfc[x / (2 √(D(t) t))], D(t) = D28 (t28/t)^α  →  solve for t
    const D28 = inp.Da * 1e-12, t28 = 28 * 86400, x = inp.cover / 1000, zc = erfinv(1 - inp.Ccr / inp.Cs);
    const Dt = x * x / (4 * zc * zc); // required D(t)·t
    const tInit = inp.ageing > 0 ? Math.pow(Dt / (D28 * Math.pow(t28, inp.ageing)), 1 / (1 - inp.ageing)) / YEAR_S : Dt / D28 / YEAR_S;
    const clAt = (xx, yrs) => { const tt = Math.max(yrs, 0.01) * YEAR_S, D = D28 * Math.pow(t28 / tt, inp.ageing); return inp.Cs * erfc(xx / (2 * Math.sqrt(D * tt))); };
    // abrasion consumes cover: time for half the cover to be lost
    const abrLife = inp.abr > 0 ? 0.5 * inp.cover / inp.abr : 9999;
    const concLife = Math.min(tInit, abrLife, 9999);
    // --- condition states 1 (as new) … 5 (failed), one-step deterioration with the stated staying probability
    const p = inp.pStay, P = [[p, 1 - p, 0, 0, 0], [0, p, 1 - p, 0, 0], [0, 0, p, 1 - p, 0], [0, 0, 0, p, 1 - p], [0, 0, 0, 0, 1]];
    const mk = markovDeterioration({ P, initial: [1, 0, 0, 0, 0], years: Math.round(inp.life) });
    const lives = { 'S–N fatigue (with DFF)': fatLife, 'Steel wastage': steelEnd, 'Chloride initiation': Math.min(tInit, 9999), 'Abrasion of cover': abrLife };
    const [govName, govLife] = Object.entries(lives).sort((a, b) => a[1] - b[1])[0];
    const remaining = Math.max(0, govLife - inp.age);
    const ok = govLife >= inp.life;
    // strengthening quantities: thickness for the fatigue life (damage ∝ S^m ∝ t^-m for a bending-governed detail, screening), cover for chlorides
    const tReq = ucFat > 1 ? inp.t * Math.pow(ucFat, 1 / sn.m1) : inp.t;
    const coverReq = tInit < inp.life ? inp.cover * Math.pow(inp.life / tInit, (1 - inp.ageing) / 2) : inp.cover;
    const yrs = linspace(0, Math.max(inp.life, 1), 51);
    const depths = linspace(0, Math.max(120, 1.5 * inp.cover), 40);

    const advice = [
      ok
        ? eng('ok', `Service life is governed by ${govName.toLowerCase()} at ${f(govLife)} years — beyond the ${f(inp.life)} years required; ${f(remaining)} years remain at the present age of ${f(inp.age)}.`,
          `Fatigue utilisation with DFF ${f(inp.dff)} is ${f(ucFat, 2)}; chloride initiation at ${f(Math.min(tInit, 9999))} years; steel reaches minimum thickness at ${f(steelEnd)} years. Confirm with inspection at about half the remaining life of the governing mechanism.`,
          'DNV-RP-C203 (2021) §5; fib Bulletin 34 (2006); BS 6349-1-4:2021', 'Demonstrated reserve life supports life extension rather than replacement — the largest carbon saving available to an existing structure.')
        : eng('fail', `Service life is cut short by ${govName.toLowerCase()} at ${f(govLife)} years against ${f(inp.life)} required (${f(remaining)} years left at age ${f(inp.age)}).`,
          `${ucFat > 1 ? `Fatigue: reduce the stress range by ${f((1 - Math.pow(1 / ucFat, 1 / sn.m1)) * 100, 2)} % (about ${f(tReq)} mm plate), grind or peen the weld toe, or improve the detail category. ` : ''}${tInit < inp.life ? `Chlorides: ${f(coverReq)} mm cover or a lower-diffusion binder would reach ${f(inp.life)} years; for an existing structure apply cathodic protection or a surface treatment before initiation. ` : ''}${steelEnd < inp.life ? `Steel: recoat or fit anodes before year ${f(inp.coat)} or add ${f((inp.life - steelEnd) * inp.rate, 2)} mm corrosion allowance. ` : ''}${abrLife < inp.life ? 'Abrasion: add a sacrificial layer or steel/granite facing in the abrasion zone.' : ''}`,
          'DNV-RP-C203 (2021) §5 and Appendix D (improvement); fib Bulletin 34 (2006); EN 1992-1-1 §4.4; BS 6349-1-4:2021',
          `Extra plate from ${f(inp.t)} to ${f(tReq)} mm adds ${f((tReq - inp.t) * 7.85, 3)} kg/m² of steel (≈ ${f((tReq - inp.t) * 7.85 * CO2_STEEL, 3)} kg CO₂e/m²); extra cover adds ${f(Math.max(0, coverReq - inp.cover) / 1000, 3)} m³/m² of concrete (≈ ${f(Math.max(0, coverReq - inp.cover) * CO2_CONC, 3)} kg CO₂e/m²) — weld improvement and low-diffusion binders achieve the life with far less material (screening factors).`),
      eng('info', `Crack growth from a ${f(inp.a0)} mm flaw to half thickness takes ${f(crackLife)} years at the damage-equivalent range of ${f(sEq / 1e6)} MPa; chloride at the bar is ${f(clAt(x, inp.age), 2)} % of binder now (critical ${f(inp.Ccr)} %).`,
        `Set the inspection interval at no more than half the crack-growth life (${f(crackLife / 2, 2)} years) with a method able to find ${f(inp.a0)} mm flaws, and take chloride profiles at the next inspection to calibrate the diffusion coefficient.`, 'BS 7910:2019 clause 8; DNV-RP-C210 (2021) probabilistic inspection planning; fib Bulletin 34'),
      eco(ok ? 'info' : 'warn', `Expected condition after ${f(inp.life)} years is state ${f(mk.expectedState[mk.expectedState.length - 1], 2)} of 5, with ${f(mk.pWorst[mk.pWorst.length - 1] * 100, 2)} % of elements in the worst state; mean time to the worst state is ${f(mk.meanTimeToWorst[0])} years.`,
        'Plan intervention before initiation: preventive measures (recoating, cathodic protection, surface treatment) applied before the governing mechanism starts cost a fraction of repairing section loss afterwards. Use the condition projection to phase the maintenance budget.',
        'Risk-based inspection per DNV-RP-G101 / API RP 580; whole-life costing per ISO 15686-5; fib Bulletin 44 (condition control)',
        'Preventive maintenance avoids break-out and replacement of concrete and steel: repair after corrosion has started typically consumes several times the material of early protection.'),
    ];

    return {
      kpis: [
        { key: 'life', label: 'Governing service life from new', value: r(govLife), unit: 'years', limit: inp.life, sense: 'min', explain: 'Shortest of the fatigue, steel-wastage, chloride and abrasion lives (crack growth sets the inspection interval instead). It is the age at which the first mechanism reaches its limit.' },
        { key: 'remaining', label: 'Remaining life at the present age', value: r(remaining), unit: 'years', explain: 'Governing life less the present age: the time left before intervention is needed on the governing mechanism.' },
        { key: 'damage', label: 'Miner fatigue damage over the required life', value: r(Dlife, 4), unit: '-', explain: 'Sum of cycle ratios n/N over the service life without the design fatigue factor. Failure is predicted at 1.0.' },
        { key: 'ucFat', label: 'Fatigue utilisation (DFF × damage)', value: r(ucFat, 4), unit: '-', limit: 1, sense: 'max', explain: 'Damage multiplied by the design fatigue factor. Above 1.0 the detail does not reach the required life with the required safety.' },
        { key: 'fatLife', label: 'Design fatigue life (with DFF)', value: r(fatLife), unit: 'years', explain: 'Years to reach a Miner sum of 1/DFF. Capped at 9999 for negligible damage.' },
        { key: 'crackLife', label: 'Crack-growth life from the present flaw', value: r(crackLife), unit: 'years', explain: 'Time for the assumed flaw to grow to half thickness by Paris-law growth. It sets the inspection interval, not the design life.' },
        { key: 'steelLife', label: 'Age at which steel reaches minimum thickness', value: r(steelEnd), unit: 'years', explain: 'Coating life plus the time to corrode the allowance at the stated rate. 9999 means not reached within 500 years.' },
        { key: 'tNow', label: 'Remaining steel thickness now', value: r(Math.max(0, wall.remainingThickness) * 1000), unit: 'mm', explain: 'As-built thickness less the wastage accumulated to the present age. Compare with thickness measurements from inspection.' },
        { key: 'tInit', label: 'Time to chloride-induced corrosion initiation', value: r(Math.min(tInit, 9999)), unit: 'years', explain: 'Age at which the chloride content at the reinforcement reaches the critical value. After it, bars corrode and cover cracks within typically 5–15 years.' },
        { key: 'clNow', label: 'Chloride content at the bar now', value: r(clAt(x, inp.age), 3), unit: '% binder', limit: inp.Ccr, sense: 'max', explain: 'Predicted chloride at cover depth at the present age. Compare with dust samples; above the critical content corrosion has started.' },
        { key: 'abrLife', label: 'Time to abrade half the cover', value: r(abrLife), unit: 'years', explain: 'Age at which abrasion has removed half the concrete cover, after which chloride protection is lost quickly. 9999 means no abrasion.' },
        { key: 'state', label: 'Expected condition state at end of life (1–5)', value: r(mk.expectedState[mk.expectedState.length - 1], 3), unit: '-', explain: 'Mean condition grade of the element population at the end of the required life: 1 as new, 5 failed. Drives the maintenance budget.' },
      ],
      checks: [
        check('Fatigue (Miner sum × DFF)', ucFat, 1, '-', 'DNV-RP-C203 (2021) §5.1; DNV-ST-0126 design fatigue factors', 'Guards against fatigue cracking of the welded detail from wave cycling. Exceeded, a through-thickness crack is expected within the life.'),
        check('Steel thickness over the required life', inp.life, steelEnd, 'years', 'BS 6349-1-4:2021 §4; ISO 12944', 'Guards against section loss below the minimum thickness needed for strength. Exceeded, the member is overstressed or perforates.'),
        check('Chloride initiation over the required life', inp.life, Math.min(tInit, 9999), 'years', 'fib Bulletin 34 (2006) §B2; EN 206 XS classes', 'Guards against reinforcement corrosion starting within the service life. Once initiated, spalling and loss of bar section follow.'),
        check('Abrasion of cover over the required life', inp.life, abrLife, 'years', 'BS 6349-1-4:2021; CIRIA C683 §3.6', 'Guards against sand, shingle or ice wearing away the concrete cover, which removes the protection of the reinforcement.'),
      ],
      plots: [
        { type: 'line', title: 'Deterioration with age', x: { label: 'Age', unit: 'years' }, y: { label: 'Utilisation of each limit', unit: '-' },
          series: [{ name: 'Fatigue (DFF × damage)', x: yrs.map((y) => r(y, 3)), y: yrs.map((y) => r(inp.dff * rateD * y, 3)) }, { name: 'Steel wastage / allowance', x: yrs.map((y) => r(y, 3)), y: yrs.map((y) => r(corrosionWastage({ time: y, ...corr }).depth * 1000 / (inp.t - inp.tMin), 3)) }, { name: 'Chloride at bar / critical', x: yrs.map((y) => r(y, 3)), y: yrs.map((y) => r(clAt(x, y) / inp.Ccr, 3)) }, { name: 'Limit', x: [0, yrs[50]], y: [1, 1], dash: true }, { name: 'Present age', x: [inp.age, inp.age], y: [0, 1], dash: true }],
          explain: 'Each mechanism as a share of its limit over the life of the structure. The first curve to cross 1.0 ends the service life; the vertical line marks today.' },
        { type: 'line', title: 'Chloride profile in the cover', x: { label: 'Depth from surface', unit: 'mm' }, y: { label: 'Chloride', unit: '% binder' },
          series: [{ name: `Now (${f(inp.age)} yr)`, x: depths.map((d) => r(d, 3)), y: depths.map((d) => r(clAt(d / 1000, inp.age), 3)) }, { name: `End of life (${f(inp.life)} yr)`, x: depths.map((d) => r(d, 3)), y: depths.map((d) => r(clAt(d / 1000, inp.life), 3)) }, { name: 'Critical content', x: [0, depths[39]], y: [inp.Ccr, inp.Ccr], dash: true }, { name: 'Bar', x: [inp.cover, inp.cover], y: [0, inp.Cs], dash: true }],
          explain: 'Chloride content with depth now and at the end of the required life. Corrosion initiates when the profile reaches the critical line at the bar position.' },
        { type: 'area', title: 'Condition-state distribution over time', x: { label: 'Age', unit: 'years' }, y: { label: 'Share of elements', unit: '-' },
          series: [0, 1, 2, 3, 4].map((sIdx) => ({ name: `State ${sIdx + 1}`, x: mk.distribution.map((_, y) => y), y: mk.distribution.map((d) => r(d[sIdx], 3)) })),
          explain: 'Share of elements in each condition grade through the life, from the Markov model. The growth of states 4 and 5 is the repair workload to plan for.' },
      ],
      tables: [
        { title: 'Life by mechanism', columns: [{ key: 'm', label: 'Mechanism' }, { key: 'life', label: 'Life from new', unit: 'years' }, { key: 'rem', label: 'Remaining', unit: 'years' }, { key: 'ok', label: 'Meets required life' }],
          rows: Object.entries(lives).map(([m, l]) => ({ m, life: r(l), rem: r(Math.max(0, l - inp.age)), ok: l >= inp.life ? 'yes' : 'no' })), explain: 'Predicted life of each deterioration mechanism and what is left at the present age. Sort by remaining life to prioritise inspection and preventive work.' },
        { title: 'Fatigue parameters', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['Cycles in the service life', n0, '-'], ['Weibull scale q', q / 1e6, 'MPa'], ['S–N slope m1', sn.m1, '-'], ['Knee stress (nominal)', sn.sKnee / 1e6, 'MPa'], ['Damage-equivalent range at the hot spot', sEq / 1e6, 'MPa'], ['Paris exponent', law.m, '-'], ['Thickness for the fatigue life', tReq, 'mm'], ['Cover for the chloride life', coverReq, 'mm']].map(([qq, v, u]) => ({ q: qq, v: r(v), u })),
          explain: 'Intermediate fatigue quantities and the thickness and cover that would meet the required life, for checking against a hand calculation.' },
      ],
      notes: [`${govName} governs. S–N curve ${inp.curve} (${inp.env}), thickness and SCF corrections included.`, `Chloride model: D28 = ${f(inp.Da)}e-12 m²/s, ageing ${f(inp.ageing)}; a constant-D calculation would give ${f(Dt / D28 / YEAR_S)} years.`],
      advice,
      publish: { 'life.fatigueYears': fatLife, 'marine.remainingLife': remaining },
    };
  },
  verify: [
    { name: 'Linear corrosion: age at minimum thickness = coating life + allowance / rate', source: 'Hand calculation: 20 mm plate, 12 mm minimum, 0.3 mm/yr after a 10-year coating life → 10 + 8/0.3 = 36.67 years',
      inputs: { t: 20, tMin: 12, rate: 0.3, coat: 10, a0: 0.5 }, expect: [{ kpi: 'steelLife', value: 36.67, tol: 0.01 }] },
    { name: 'Fick\'s second law: initiation time t = [x / (2 erfc⁻¹(Ccr/Cs))]² / D', source: 'Crank (1975) Eq. 2.45; fib Bulletin 34 Eq. B2.1-1 with constant D. Cover 50 mm, D = 1e-12 m²/s, Cs = 0.6, Ccr = 0.1: erfc⁻¹(0.16667) = 0.9779 (erf 0.9779 = 0.8333), t = (0.05/1.9558)²/1e-12 = 6.536e8 s = 20.71 years',
      inputs: { cover: 50, Da: 1, ageing: 0, Cs: 0.6, Ccr: 0.1 }, expect: [{ kpi: 'tInit', value: 20.71, tol: 0.01 }] },
    { name: 'Single-slope Weibull closed form D = n0 q^m Γ(1 + m/h)/a', source: 'DNV-RP-C203 (2021) Eq. 5.1.1 for a free-corrosion curve (single slope). Curve F free corrosion: log a = 11.378, m = 3; t = 25 mm (no thickness correction), SCF 1; Δσ0 = 100 MPa, h = 1, Tz = 6 s, 20 years: n0 = 1.0519e8, q = 100/ln(n0) = 5.413 MPa, D = 1.0519e8 × 5.413³ × Γ(4)/10^11.378 = 0.4192',
      inputs: { curve: 'F', env: 'free-corrosion', t: 25, tMin: 10, scf: 1, sMax: 100, hW: 1, Tz: 6, life: 20, age: 0 }, expect: [{ kpi: 'damage', value: 0.4192, tol: 0.02 }] },
  ],
};
