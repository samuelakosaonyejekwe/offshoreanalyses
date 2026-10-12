// Tensile armour stress and fatigue: wire stress ranges from tension and curvature load bins (friction stick–slip
// plus local wire bending), mean-stress correction, S–N life and Palmgren–Miner damage for every armour layer.
import { armourFatigueCoefficients, armourFatigueDamage, genericArmourSN, snCycles, multiaxialCriterion } from '../../engine/xsection/flexfatigue.js';
import { axisymmetricResponse } from '../../engine/xsection/flexible.js';
import { stackInputs, stackFrom, loadInputs, loadsFrom, muInput, kpi, chk, adv, plot, bar, tbl, col, run, need, fmt, round, linspace, sev, LIMIT_AXI, GENERIC, STD_17J, STD_17B, CONFIRM, CO2_MGO, MGO_T_PER_DAY } from './_common.js';

const ENV_HELP = 'There are no built-in sweet or sour armour-wire S–N curves. For a corrosive annulus the knock-down on stress must come from qualified fatigue tests on the actual wire in the actual annulus environment.';

/** S–N inputs shared with the service-life module. */
export function snInputs({ source = 'generic' } = {}) {
  const g = 'Fatigue resistance';
  return [
    { key: 'snSource', label: 'Armour-wire S–N curve', type: 'select', default: source, options: [{ value: 'generic', label: 'Generic textbook estimate from tensile strength (unqualified)' }, { value: 'custom', label: 'Curve N = a·S⁻ᵐ from the wire qualification programme (enter a and m)' }], group: g, role: 'material',
      help: 'The generic curve is a textbook estimate for polished steel in air and is for demonstration only. Use the qualified curve from the wire test programme for design.' },
    { key: 'snLogA', label: 'log₁₀ a (S in MPa)', unit: '-', type: 'number', default: 21.6, min: 8, max: 40, group: g, role: 'material', showIf: { key: 'snSource', in: ['custom'] }, help: 'Intercept of the S–N curve for the actual wire and annulus environment, stress range in MPa. The default (with m = 6: about 10⁶ cycles at 400 MPa) is a placeholder of realistic magnitude for high-strength wire in dry air — it is not a qualified curve.' },
    { key: 'snM', label: 'Slope m', unit: '-', type: 'number', default: 6, min: 2, max: 20, group: g, role: 'material', showIf: { key: 'snSource', in: ['custom'] } },
    { key: 'env', label: 'Annulus environment', type: 'select', default: 'air', options: [{ value: 'air', label: 'Dry annulus (in-air curve)' }, { value: 'sweet', label: 'Sweet: water with CO₂' }, { value: 'sour', label: 'Sour: water with H₂S' }], group: g, role: 'ic', showIf: { key: 'snSource', in: ['generic'] },
      help: ENV_HELP },
    { key: 'envFactor', label: 'Environment knock-down on stress (from qualified tests)', unit: '-', type: 'number', default: 1, min: 1, max: 5, group: g, role: 'material', showIf: { key: 'snSource', in: ['generic'] },
      uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Factor ≥ 1 dividing the fatigue strength for a sweet or sour annulus. It has no generic value: enter the result of the qualification tests. Must exceed 1 when the environment is not dry air.' },
    { key: 'marin', label: 'Surface / size / reliability factor on the endurance limit', unit: '-', type: 'number', default: 0.3, min: 0.1, max: 1, group: g, role: 'material', showIf: { key: 'snSource', in: ['generic'] },
      uncertain: { dist: 'normal', cov: 0.15 }, help: 'Combined Marin factor applied to the polished-specimen endurance limit of the generic estimate (as-rolled surface, wire size, reliability). 1.0 reproduces the textbook polished-specimen value.' },
    { key: 'meanMethod', label: 'Mean-stress correction', type: 'select', default: 'goodman', options: [{ value: 'goodman', label: 'Goodman' }, { value: 'gerber', label: 'Gerber' }, { value: 'none', label: 'None (curve already at the service R-ratio)' }], group: g, role: 'control' },
    { key: 'scf', label: 'Stress concentration factor', unit: '-', type: 'number', default: 1.3, min: 1, max: 3, group: g, role: 'material', help: 'On the stress range, e.g. at the end-fitting entrance or a wire weld.' },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 10, min: 1, max: 20, group: 'Acceptance criteria', role: 'criteria', help: 'API 17J requires a factor of 10 on the armour fatigue life. ' + CONFIRM },
    { key: 'designLife', label: 'Design service life', unit: 'years', type: 'number', default: 25, min: 1, max: 60, group: 'Acceptance criteria', role: 'criteria', bind: 'life.designYears' },
  ];
}
/** S–N curve record for the engine (S in Pa). */
export function snFrom(inp, uts) {
  if (inp.snSource === 'custom') return { curve: { a: 10 ** (inp.snLogA + 6 * inp.snM), m: inp.snM }, label: `user-entered curve log a = ${inp.snLogA}, m = ${inp.snM} (to be traceable to qualified wire tests)`, generic: false };
  need(inp.env === 'air' || inp.envFactor > 1, `A ${inp.env} annulus needs an environment knock-down greater than 1 taken from qualified wire fatigue tests — there is no generic value. Enter it, or supply the qualified S–N curve directly.`);
  need(uts > 0, 'The generic S–N estimate needs the wire tensile strength in the layer table.');
  const curve = run(() => genericArmourSN({ uts, environment: inp.env, environmentFactor: inp.env === 'air' ? undefined : inp.envFactor, marin: inp.marin ?? 1 }), 'S–N curve');
  return { curve, label: `generic unqualified ${inp.env} estimate`, generic: true };
}
export const cap = (v, hi = 1e5) => (Number.isFinite(v) ? Math.min(v, hi) : hi);

export default {
  id: 's14-armour-fatigue',
  suite: 14,
  title: 'Tensile armour stress and fatigue',
  summary: 'Turns tension and curvature load bins into armour-wire stress ranges and accumulates fatigue damage for each tensile armour layer.',
  covers: [8],
  fidelity: 'analytical',
  standards: ['API Spec 17J (4th ed.) §5.3.4 (service life, design fatigue factor 10)', 'API RP 17B (5th ed.) §8 (fatigue methodology)', 'DNV-RP-C203 §2.4 (S–N format)'],
  models: ['Armour-wire friction stick–slip stress (Sævik 2011)', 'Local wire bending on a bent helix', 'S–N curves with Goodman / Gerber mean-stress correction', 'Palmgren–Miner damage', 'Multiaxial fatigue criteria (Sines, Crossland, Dang Van, Findley critical plane)'],
  limits: LIMIT_AXI + ' Contact pressure is constant through each cycle. No built-in sweet/sour wire S–N curves exist: the generic curve is an unqualified in-air textbook estimate and any environment knock-down must come from qualified test data. Fretting and wire-weld fatigue are not modelled. ' + GENERIC,
  publishes: ['flexible.armourFatigueLife'],
  searchTerms: 'flexible riser tensile armour fatigue S-N curve corrosion annulus',
  inputs: [
    ...stackInputs(),
    ...loadInputs({ depth: 30 }).filter((i) => i.key !== 'Te').map((i) => (i.key === 'depth' ? { ...i, bind: undefined, help: 'Depth of the fatigue hot spot (hang-off and bend-stiffener region ≈ 0–50 m; touchdown = full water depth). External pressure raises contact pressure and with it the friction stress.' } : i)),
    muInput(),
    { key: 'bins', label: 'Load bins (per year)', type: 'table', group: 'Global loading', role: 'bc',
      columns: [{ key: 'label', label: 'Sea state / case', type: 'text' }, { key: 'Te', label: 'Mean effective tension', unit: 'kN', type: 'number' }, { key: 'dKappa', label: 'Curvature range', unit: '1/m', type: 'number' }, { key: 'cycles', label: 'Cycles per year', type: 'number' }],
      default: [
        { label: 'Hs 0–2 m', Te: 450, dKappa: 0.003, cycles: 3.2e6 },
        { label: 'Hs 2–4 m', Te: 520, dKappa: 0.008, cycles: 1.3e6 },
        { label: 'Hs 4–6 m', Te: 600, dKappa: 0.02, cycles: 2.0e5 },
        { label: 'Hs 6–8 m', Te: 700, dKappa: 0.045, cycles: 1.5e4 },
        { label: 'Hs > 8 m', Te: 800, dKappa: 0.08, cycles: 600 },
      ], help: 'Mean tension, peak-to-peak curvature and yearly cycle count at the hot spot (hang-off / bend-stiffener region or touchdown) from global dynamic analysis.' },
    { key: 'wirePath', label: 'Wire path in bending', type: 'select', default: 'loxodromic', options: [{ value: 'loxodromic', label: 'Loxodromic' }, { value: 'geodesic', label: 'Geodesic' }], group: 'Interlayer contact', role: 'control' },
    ...snInputs(),
    { key: 'multiaxial', label: 'Multiaxial infinite-life screen', type: 'select', default: 'findley', options: [{ value: 'none', label: 'None' }, { value: 'sines', label: 'Sines' }, { value: 'crossland', label: 'Crossland' }, { value: 'dangvan', label: 'Dang Van' }, { value: 'findley', label: 'Findley (critical plane)' }], group: 'Multiaxial screen', role: 'control',
      help: 'Screens the largest bin (axial stress cycle plus contact pressure) against an infinite-life criterion. The constants are material data from tests; the defaults are indicative only.' },
    { key: 'maK', label: 'Normal-stress sensitivity k', unit: '-', type: 'number', default: 0.3, min: 0, max: 1.5, group: 'Multiaxial screen', role: 'material' },
    { key: 'maLimit', label: 'Criterion limit', unit: 'MPa', type: 'number', default: 350, min: 20, max: 1500, group: 'Multiaxial screen', role: 'material', help: 'Right-hand side of the criterion (e.g. torsion endurance limit). Indicative default — replace with test data.' },
  ],
  compute(inp, ctx) {
    const S = stackFrom(inp, ctx);
    need(S.idx.ta.length > 0, 'The stack has no tensile armour layer.');
    const rowsIn = (inp.bins || []).filter((b) => Number(b.cycles) > 0);
    need(rowsIn.length > 0, 'Enter at least one load bin with a positive number of cycles.');
    for (const b of rowsIn) need(Number.isFinite(Number(b.Te)) && Number(b.dKappa) >= 0, `Load bin '${b.label ?? ''}': tension and curvature range must be numbers (curvature range ≥ 0).`);
    const loads = loadsFrom({ ...inp, Te: 0 }); delete loads.tension;
    const bins = rowsIn.map((b) => ({ tension: b.Te * 1e3, curvatureRange: Number(b.dKappa), cycles: Number(b.cycles) }));
    const Ts = bins.map((b) => b.tension), Tlo = Math.min(...Ts), Thi = Math.max(...Ts);
    const per = S.idx.ta.map((i) => {
      const L = S.stack.layers[i];
      need(L.wireWidth > 0, `${L.name}: enter the wire width in the layer table — local wire bending stress needs the wire section.`);
      const sn = snFrom(inp, L.uts);
      const co = run(() => armourFatigueCoefficients(S.stack, i, { tensionLow: Tlo === Thi ? 0.9 * Tlo - 1e3 : Tlo, tensionHigh: Thi, mu: inp.mu, loads }), 'Armour stress coefficients');
      const meanTop = co.meanStress0 + co.meanStressPerTension * Thi;
      need(inp.meanMethod === 'none' || meanTop < 0.98 * L.uts, `${L.name}: the mean wire stress at ${fmt(Thi / 1e3)} kN (${fmt(meanTop / 1e6)} MPa) reaches the wire tensile strength — the armour fails statically under this tension and pressure, so a fatigue life is meaningless. Check the loads and strengths.`);
      const d = run(() => armourFatigueDamage({ wire: co.wire, bins, meanStress0: co.meanStress0, meanStressPerTension: co.meanStressPerTension, friction0: co.friction0, frictionPerTension: co.frictionPerTension,
        curve: sn.curve, uts: L.uts, meanStressMethod: inp.meanMethod, scf: inp.scf, dff: inp.dff, exposureYears: 1, wireCurve: inp.wirePath }), `${L.name} fatigue damage`);
      return { L, i, sn, co, d };
    });
    const gov = per.reduce((a, c) => (c.d.damagePerYear > a.d.damagePerYear ? c : a), per[0]);
    const life = cap(gov.d.life), lifeDff = cap(gov.d.lifeWithDff);
    const dmgDesign = gov.d.damagePerYear * inp.dff * inp.designLife;
    const rngMax = Math.max(...gov.d.stress.map((s) => s.range)) * inp.scf;
    const meanMax = Math.max(...gov.d.stress.map((s) => s.meanStress));
    const iWorst = gov.d.bins.reduce((a, b, k) => (b.damage > gov.d.bins[a].damage ? k : a), 0);
    const checks = per.map((p) => chk(`${p.L.name}: fatigue damage over the design life × DFF`, p.d.damagePerYear * inp.dff * inp.designLife, 1, '-', `${STD_17J} §5.3.4 (DFF ${inp.dff}); Palmgren–Miner`,
      'Guards against fatigue cracking of the armour wires. Above 1.0 the calculated life divided by the design fatigue factor is shorter than the design life.'));
    let ma = null;
    if (inp.multiaxial !== 'none') {
      // largest bin: axial stress amplitude and mean along the wire, contact pressure as a compressive transverse mean stress
      const kBig = gov.d.stress.reduce((a, s, k) => (s.range > gov.d.stress[a].range ? k : a), 0);
      const resp = axisymmetricResponse(S.stack, { ...loads, tension: bins[kBig].tension });
      const q = Math.max(resp.layers[gov.i].contactIn, resp.layers[gov.i].contactOut, 0) / gov.L.fill;
      const amp = gov.d.stress[kBig].range * inp.scf / 2;
      ma = run(() => multiaxialCriterion({ criterion: inp.multiaxial, amplitude: [amp, 0, 0, 0, 0, 0], mean: [gov.d.stress[kBig].meanStress, -q, 0, 0, 0, 0], k: inp.maK, limit: inp.maLimit * 1e6, nGrid: 24 }), 'Multiaxial criterion');
      checks.push(chk(`${gov.L.name}: ${inp.multiaxial} infinite-life screen, largest bin`, ma.equivalent / 1e6, inp.maLimit, 'MPa', 'Socie & Marquis, Multiaxial Fatigue (2000) Ch. 5 — constants from material tests',
        'Screens whether the largest stress cycle, combined with contact pressure, stays below the multiaxial endurance limit. Above 1.0 that bin is damaging and the Miner sum governs.'));
    }
    const lab = rowsIn.map((b, k) => String(b.label || `Bin ${k + 1}`));
    const sGrid = linspace(Math.log10(Math.max(rngMax / 30, 1e6)), Math.log10(Math.max(rngMax * 3, 2e8)), 25);
    const snPts = sGrid.map((ls) => ({ ls, n: snCycles({ stressRange: 10 ** ls, curve: gov.sn.curve }).cycles })).filter((p) => Number.isFinite(p.n) && p.n > 0);
    const ok = lifeDff >= inp.designLife, noDamage = !(gov.d.damagePerYear > 0);
    const kAllow = dmgDesign > 0 ? (1 / dmgDesign) ** (1 / gov.sn.curve.m) : Infinity; // factor on all stress ranges that would bring the damage to 1
    return {
      kpis: [
        kpi('lifeDff', 'Fatigue life with design fatigue factor', lifeDff, 'years', `Calculated life of ${gov.L.name} divided by the design fatigue factor of ${inp.dff}. It must exceed the design life. Values of 100 000 mean no damaging cycles.`, inp.designLife, 'min'),
        kpi('life', 'Calculated fatigue life (unfactored)', life, 'years', 'Years to a Miner sum of 1.0 for the governing layer, before any safety factor.'),
        kpi('damageDesign', 'Damage over the design life × DFF', dmgDesign, '-', 'Miner sum for the design life multiplied by the design fatigue factor. Must not exceed 1.0.', 1, 'max'),
        kpi('damagePerYear', 'Damage per year', gov.d.damagePerYear, '1/yr', 'Unfactored Miner damage accumulated each year by the governing layer.'),
        kpi('stressRangeMax', 'Largest wire stress range', rngMax / 1e6, 'MPa', 'Largest hot-spot stress range of any bin, including the stress concentration factor. Friction stress saturates at full slip; local wire bending keeps growing with curvature.'),
        kpi('meanStressMax', 'Largest wire mean stress', meanMax / 1e6, 'MPa', 'Highest mean wire stress from tension and pressure; it raises the damage through the mean-stress correction.'),
        kpi('govLayer', 'Governing layer (position in the stack)', gov.i + 1, '-', `Row number of the governing tensile armour layer in the layer table (${gov.L.name}).`),
        kpi('multiaxialUc', 'Multiaxial screen utilisation', ma ? ma.utilisation : 0, '-', 'Equivalent multiaxial stress of the largest bin over the criterion limit; below 1.0 that cycle is non-damaging by the chosen criterion. Zero when the screen is switched off.'),
      ],
      checks,
      plots: [
        bar('Fatigue damage by load bin', lab, 'Damage per year', '1/yr', per.map((p) => ({ name: p.L.name, y: p.d.bins.map((b) => round(b.damage, 4)) })),
          'Share of the yearly damage from each bin. The bin that dominates is usually a moderate sea state with many cycles, not the storm.'),
        bar('Wire stress range by load bin', lab, 'Stress range', 'MPa', per.map((p) => ({ name: p.L.name, y: p.d.stress.map((s) => round(s.range * inp.scf / 1e6)) })),
          'Hot-spot stress range per bin. Once the wires slip the friction part stops growing, which flattens the increase with curvature.'),
        plot('line', 'S–N curve and load bins', 'log₁₀ cycles', '-', 'log₁₀ stress range (MPa)', '-', [
          { name: 'S–N curve', x: snPts.map((p) => round(Math.log10(p.n))), y: snPts.map((p) => round(p.ls - 6)) },
          { name: 'Load bins (equivalent range, applied cycles)', x: gov.d.bins.map((b, k) => round(Math.log10(Math.max(bins[k].cycles, 1)))), y: gov.d.bins.map((b) => round(Math.log10(Math.max(b.equivalentRange, 1)) - 6)) },
        ], 'Load bins against the S–N curve on log axes. Points close to or right of the line are the damaging ones; points far below the curve contribute little.'),
      ],
      tables: [tbl(`Load bins — ${gov.L.name}`, [col('bin', 'Bin'), col('Te', 'Tension', 'kN'), col('dk', 'Curvature range', '1/m'), col('n', 'Cycles/yr'), col('mean', 'Mean stress', 'MPa'), col('range', 'Stress range', 'MPa'), col('eq', 'Equivalent range', 'MPa'), col('N', 'Cycles to failure'), col('d', 'Damage/yr'), col('slip', 'Slip')],
        gov.d.bins.map((b, k) => ({ bin: lab[k], Te: round(bins[k].tension / 1e3), dk: bins[k].curvatureRange, n: bins[k].cycles, mean: round(gov.d.stress[k].meanStress / 1e6), range: round(gov.d.stress[k].range * inp.scf / 1e6), eq: round(b.equivalentRange / 1e6),
          N: Number.isFinite(b.cyclesToFailure) ? round(b.cyclesToFailure) : 'no damage', d: round(b.damage, 4), slip: gov.d.stress[k].slipped ? 'full slip' : 'stick / partial' })),
        'One row per load bin for the governing layer: wire stresses, mean-stress-corrected range, allowable cycles and yearly damage.')],
      notes: [`${gov.L.name} governs; bin '${lab[iWorst]}' contributes most of the damage.`, `S–N basis: ${gov.sn.label}.${gov.sn.generic ? ' This curve is unqualified — results are indicative only.' : ''}`,
        inp.snSource === 'generic' && inp.env !== 'air' ? `Environment knock-down ${inp.envFactor} on stress as entered — it must be traceable to qualified tests.` : ENV_HELP, GENERIC],
      advice: [
        adv('engineering', ok ? (lifeDff < 1.5 * inp.designLife ? 'warn' : 'ok') : 'fail', `Factored fatigue life of ${gov.L.name} is ${noDamage ? 'unlimited (no damaging cycles)' : fmt(lifeDff) + ' years'} against a design life of ${inp.designLife} years (largest stress range ${fmt(rngMax / 1e6)} MPa).`,
          ok ? `${noDamage ? 'All stress ranges are below the cut-off of the curve' : 'Reserve factor ' + fmt(lifeDff / inp.designLife) + ' on life'}. Replace the generic S–N curve by the qualified wire curve for the annulus environment and re-run before relying on it.` : `Stress ranges must fall by about ${fmt((1 - kAllow) * 100)} % (a factor ${fmt(kAllow)}): lengthen or stiffen the bend stiffener to cut curvature range, lower the friction (fresh anti-wear tape), or keep the annulus dry.`,
          `${STD_17J} §5.3.4; ${STD_17B}; Sævik (2011)`),
        adv('economic', ok ? 'info' : 'warn', noDamage ? 'No load bin exceeds the fatigue limit of the S–N curve used: armour fatigue does not limit the service life on this basis.' : `At ${fmt(gov.d.damagePerYear * 100)} % unfactored damage per year the riser reaches its factored limit after ${fmt(lifeDff)} years.`,
          ok ? `A life-extension case can claim ${noDamage ? 'the full extension period' : 'up to ' + fmt(Math.max(0, lifeDff - inp.designLife)) + ' extra years'} on fatigue, subject to annulus testing; that defers a riser replacement of the order of 10–20 vessel days.` : `Plan a replacement or re-termination at year ${fmt(lifeDff)} unless the curvature is reduced; each replacement is of the order of 10–20 vessel days plus deferred production.`,
          'Life-extension assessment per NORSOK Y-002 / ISO 15686-5 whole-life costing; risk-based inspection per DNV-RP-G101',
          `Extending service instead of replacing saves the whole pipe (steel and polymer that cannot be readily recycled as a composite) and about ${fmt(15 * MGO_T_PER_DAY * CO2_MGO)} t CO₂ for a 15-day campaign (screening: ${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂ per t).`),
      ],
      publish: { 'flexible.armourFatigueLife': lifeDff },
    };
  },
  verify: [{
    // Friction switched off (μ = 0): the range is local wire bending only. Wire 12 × 4 mm at α = 35°, E = 207 GPa,
    // Δκ = 0.02 1/m, loxodromic path. Normal bending at the extrados: E·(t/2)·cos⁴α·Δκ = 207e9·0.002·0.4503·0.02 = 3.73 MPa.
    // Transverse bending at the neutral axis: E·(b/2)·(1 + sin²α)·cosα·Δκ = 207e9·0.006·1.3290·0.81915·0.02 = 27.04 MPa
    // (governs). S–N: log a = 12 (MPa), m = 3 → N = 1e12/27.04³ = 5.058e7; 1e6 cycles/yr → damage 0.01977 per year.
    name: 'Local wire bending stress range and single-block Miner damage (no friction)',
    source: 'Sævik (2011) Marine Structures 24, 97–111 (wire curvature on a bent helix); Basquin/Miner hand calculation',
    inputs: { mu: 0, snSource: 'custom', snLogA: 12, snM: 3, meanMethod: 'none', scf: 1, multiaxial: 'none', bins: [{ label: 'single', Te: 500, dKappa: 0.02, cycles: 1e6 }] },
    expect: [{ kpi: 'stressRangeMax', value: 27.04, tol: 0.005 }, { kpi: 'damagePerYear', value: 0.01977, tol: 0.015 }],
  }],
};
