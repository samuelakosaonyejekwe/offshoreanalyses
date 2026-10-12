// Fatigue and cyclic loading: stress ranges in tubes and armour wires from tension and curvature ranges (tension
// share, own-axis bending, friction), S–N curves and Palmgren–Miner summation; optional rainflow counting of a
// measured or simulated tension history.
import { umbilicalTubeFatigue } from '../../engine/xsection/umbilical.js';
import { snCurve, SN_CURVE_NAMES } from '../../engine/fatigue/sn.js';
import { snCycles } from '../../engine/xsection/flexfatigue.js';
import { rainflowFourPoint, rangeHistogram } from '../../engine/fatigue/rainflow.js';
import { fatigueLife } from '../../engine/fatigue/damage.js';
import { sectionInputs, buildSection, sectionMasses, adv, check, worst, sev, fmt, round, run, need, linspace, YEAR_S, CO2_STEEL, LIMIT_RIGID, STD_UMB, vesselCo2Text } from './_common.js';

const ENV = [{ value: 'air', label: 'Air (dry interstices)' }, { value: 'seawater-cp', label: 'Sea water with cathodic protection' }, { value: 'free-corrosion', label: 'Free corrosion' }];
// DNV-RP-C203 curve record → { a [Pa^m], m, m2, nKnee } for the cross-section fatigue engine (stress in Pa)
const toCurve = (c) => ({ a: 10 ** (c.logA1 + 6 * c.m1), m: c.m1, ...(c.m2 ? { m2: c.m2, nKnee: c.kneeCycles } : {}) });

export default {
  id: 's12-fatigue', suite: 12, title: 'Fatigue and cyclic loading',
  summary: 'Miner damage and fatigue life of every steel tube and armour wire from the long-term tension and curvature spectrum, with the critical element and load case.',
  covers: [6], fidelity: 'analytical',
  standards: [`${STD_UMB} §7.4 (fatigue analysis, design fatigue factor 10)`, 'DNV-RP-C203 (2021) §2 (S–N curves, Miner summation)', 'ASTM E1049-85 (2017) (rainflow counting)'],
  models: ['S–N relationships', 'Palmgren–Miner damage accumulation', 'Rainflow counting', 'Stress concentration factors', 'Armour-wire and tube fatigue: tension share, own-axis bending and stick–slip friction stress'],
  limits: `${LIMIT_RIGID} S–N curves are DNV-RP-C203 classes used as stand-ins: seam-welded super-duplex tube and galvanised wire need the manufacturer's qualified curves. Fretting, wear and local-strain (low-cycle) fatigue from reeling are not modelled here.`,
  publishes: ['life.fatigueYears', 'umbilical.fatigueDamageRate'],
  searchTerms: 'umbilical steel tube fatigue S-N Miner friction stress',
  inputs: [
    ...sectionInputs(),
    { key: 'spectrum', label: 'Long-term load spectrum at the critical section', type: 'table', group: 'Cyclic loading', role: 'load',
      columns: [{ key: 'name', label: 'Load case', type: 'text' }, { key: 'meanT', label: 'Mean tension', unit: 'kN', type: 'number' }, { key: 'dT', label: 'Tension range', unit: 'kN', type: 'number' }, { key: 'dK', label: 'Curvature range', unit: '1/m', type: 'number' }, { key: 'cycles', label: 'Cycles per year', type: 'number' }],
      default: [
        { name: 'Hs 0–1.5 m', meanT: 80, dT: 2, dK: 0.0008, cycles: 2.6e6 },
        { name: 'Hs 1.5–2.5 m', meanT: 80, dT: 4, dK: 0.0016, cycles: 1.7e6 },
        { name: 'Hs 2.5–4 m', meanT: 82, dT: 6.5, dK: 0.003, cycles: 6.0e5 },
        { name: 'Hs 4–6 m', meanT: 84, dT: 10, dK: 0.0045, cycles: 1.2e5 },
        { name: 'Hs 6–8 m', meanT: 86, dT: 14, dK: 0.007, cycles: 1.5e4 },
        { name: 'Hs > 8 m', meanT: 90, dT: 21, dK: 0.0105, cycles: 9e2 },
      ],
      help: 'Tension and curvature ranges at the hang-off (or touchdown) per sea-state class with their yearly cycle counts, from the dynamic analyses.' },
    { key: 'history', label: 'Tension history for rainflow counting (optional)', type: 'series', default: { x: [], y: [] }, xLabel: 'Time', xUnit: 's', yLabel: 'Tension', yUnit: 'kN', group: 'Cyclic loading', role: 'load',
      help: 'Measured or simulated tension record. Its rainflow cycles, scaled to a year, are added to the table spectrum.' },
    { key: 'histFraction', label: 'Fraction of the year represented by the history', unit: '-', type: 'number', default: 0.1, min: 0, max: 1, group: 'Cyclic loading', role: 'load' },
    { key: 'kPerT', label: 'Curvature range per unit tension range (history)', unit: '1/m per kN', type: 'number', default: 0.0007, min: 0, max: 0.05, group: 'Cyclic loading', role: 'load', help: 'Converts rainflow tension ranges to accompanying curvature ranges at the section (from the vessel-motion module).' },
    { key: 'snTube', label: 'S–N curve — steel tubes', type: 'select', default: 'C1', options: SN_CURVE_NAMES.map((n) => ({ value: n, label: `DNV-RP-C203 ${n}` })), group: 'Fatigue resistance', role: 'material',
      help: 'Stand-in class for seam-welded tube; use the manufacturer-qualified curve where available.' },
    { key: 'snWire', label: 'S–N curve — armour wires / strength members', type: 'select', default: 'B1', options: SN_CURVE_NAMES.map((n) => ({ value: n, label: `DNV-RP-C203 ${n}` })), group: 'Fatigue resistance', role: 'material' },
    { key: 'snEnv', label: 'Environment', type: 'select', default: 'seawater-cp', options: ENV, group: 'Fatigue resistance', role: 'material' },
    { key: 'scf', label: 'Stress concentration factor (tube welds)', unit: '-', type: 'number', default: 1.15, min: 1, max: 4, group: 'Fatigue resistance', role: 'material', uncertain: { dist: 'lognormal', cov: 0.1 }, help: 'Girth / seam weld misalignment of tubes; armour wires use 1.0.' },
    { key: 'mu', label: 'Friction coefficient between elements', unit: '-', type: 'number', default: 0.2, min: 0, max: 0.6, group: 'Fatigue resistance', role: 'material', uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 10, min: 1, max: 20, group: 'Acceptance criteria', role: 'criteria', help: '10 for non-inspectable umbilical components (ISO 13628-5 / API 17E pattern) — confirm against the governing edition.' },
    { key: 'life', label: 'Design service life', unit: 'years', type: 'number', default: 25, min: 1, max: 60, group: 'Acceptance criteria', role: 'criteria', bind: 'life.designYears' },
    { key: 'D0', label: 'Damage already accumulated', unit: '-', type: 'number', default: 0.005, min: 0, max: 1, group: 'Initial condition', role: 'ic', help: 'Characteristic Miner sum consumed before service: reeling, load-out, installation, earlier service.' },
    { key: 'lengthKm', label: 'Dynamic section length', unit: 'km', type: 'number', default: 0.74, min: 0.02, max: 50, group: 'Quantities', role: 'operation', bind: 'cable.length' },
  ],
  compute(inp, ctx) {
    const b = buildSection(inp, ctx), s = b.sec;
    const targets = s.elements.filter((e) => e.component.yieldStress && e.component.A > 0);
    need(targets.length > 0, 'No steel tube or wire (component with a yield stress) to assess for fatigue.');
    const rows = (inp.spectrum || []).filter((r) => r.cycles > 0);
    for (const r of rows) need(r.dT >= 0 && r.dK >= 0 && r.meanT >= 0, `Load case '${r.name}': ranges and mean tension cannot be negative.`);
    const bins = rows.map((r) => ({ name: r.name, tensionRange: r.dT * 1e3, curvatureRange: r.dK, cycles: r.cycles, meanTension: r.meanT * 1e3 }));
    // rainflow of the optional history, scaled to one year
    let rf = null;
    const hy = inp.history?.y || [];
    if (hy.length >= 3) {
      need(inp.histFraction > 0, 'Give the fraction of the year that the tension history represents.');
      const hx = inp.history.x || [], dur = hx.length === hy.length ? hx[hx.length - 1] - hx[0] : 0;
      need(dur > 0, 'The tension history needs increasing time values of the same length as the tension values.');
      rf = run(() => rainflowFourPoint(hy), 'Rainflow counting');
      // 24 range classes, each assigned its upper edge (conservative binning)
      const hist = rangeHistogram(rf.cycles, { bins: 24 });
      const mean = hy.reduce((a, v) => a + v, 0) / hy.length, scale = inp.histFraction * YEAR_S / dur;
      for (let i = 0; i < hist.counts.length; i++) {
        const rg = hist.edges[i + 1];
        if (hist.counts[i] > 0) bins.push({ name: `Rainflow ≤ ${fmt(rg)} kN`, tensionRange: rg * 1e3, curvatureRange: rg * inp.kPerT, cycles: hist.counts[i] * scale, meanTension: Math.max(mean, 0) * 1e3 });
      }
    }
    need(bins.length > 0, 'The load spectrum is empty: enter at least one load case with cycles, or a tension history.');
    const cTube = toCurve(run(() => snCurve(inp.snTube, inp.snEnv), 'Tube S–N curve')), cWire = toCurve(run(() => snCurve(inp.snWire, inp.snEnv), 'Wire S–N curve'));
    const res = targets.map((e) => {
      const tube = e.component.id > 0, curve = tube ? cTube : cWire, scf = tube ? inp.scf : 1;
      const f = run(() => umbilicalTubeFatigue(s, { element: e.name, bins, curve, mu: inp.mu, scf, dff: inp.dff, exposureYears: 1 }), `Fatigue of ${e.name}`);
      const fl = fatigueLife({ damageRate: f.damagePerYear, designLife: inp.life, dff: inp.dff, accumulatedDamage: inp.D0 });
      return { name: e.name, tube, curve, scf, f, fl, rate: f.damagePerYear, maxRange: Math.max(...f.ranges.map((r) => r.total)) * scf };
    });
    const gov = res.reduce((a, r) => (r.rate > a.rate ? r : a), res[0]);
    const checks = res.map((r) => check(`${r.name}: fatigue damage over the design life`, inp.dff * (inp.D0 + r.rate * inp.life), 1, '-', `${STD_UMB} §7.4; DNV-RP-C203 (2021) §2.2 (Miner sum × DFF ≤ 1)`,
      r.tube ? 'Accumulated Miner damage of the tube including the design fatigue factor. Above 1.0 a through-wall fatigue crack — loss of a hydraulic or chemical line — is predicted within the design life at the required safety level.'
        : 'Accumulated Miner damage of the wire including the design fatigue factor. Above 1.0 wire breaks are predicted, shedding load to the remaining wires and the tubes.'));
    const gc = worst(checks);
    const gb = gov.f.bins.map((q, i) => ({ ...bins[i], N: q.cyclesToFailure, d: q.damage, ...gov.f.ranges[i] }));
    const dTot = gb.reduce((a, q) => a + q.d, 0) || 1e-300;
    const worstBin = gb.reduce((a, q) => (q.d > a.d ? q : a), gb[0]);
    const lifeDff = Number.isFinite(gov.fl.designLife) ? gov.fl.designLife : 1e6, remaining = Number.isFinite(gov.fl.remainingLife) ? gov.fl.remainingLife : 1e6;
    const yrs = linspace(0, Math.max(inp.life * 1.5, 1), 31);
    const sGrid = linspace(Math.log10(Math.max(1e6, gov.maxRange / 20)), Math.log10(Math.max(gov.maxRange * 2, 5e7)), 40).map((l) => 10 ** l);
    const masses = sectionMasses(b), steelT = masses.steel * inp.lengthKm;
    return {
      kpis: [
        { key: 'damageRate', label: 'Damage per year, governing element', value: round(gov.rate, 5), unit: '1/yr', explain: `Characteristic Miner damage accumulated each year in "${gov.name}" without the design fatigue factor. Its inverse is the calculated life.` },
        { key: 'lifeChar', label: 'Calculated fatigue life (no DFF)', value: round(Math.min(gov.fl.characteristicLife, 1e6)), unit: 'years', explain: 'Years to a Miner sum of 1.0 on the design S–N curve. This is not the usable life: divide by the design fatigue factor.' },
        { key: 'lifeDesign', label: 'Fatigue design life (with DFF)', value: round(lifeDff), unit: 'years', limit: inp.life, sense: 'min', explain: `Calculated life divided by the design fatigue factor ${inp.dff}. It must exceed the service life of ${inp.life} years.` },
        { key: 'remaining', label: 'Remaining life from today', value: round(Math.max(remaining, 0)), unit: 'years', explain: 'Design life left after deducting the damage already accumulated (installation, earlier service).' },
        { key: 'utilisation', label: 'Fatigue utilisation at the design life', value: round(gc.uc), unit: '-', limit: 1, sense: 'max', explain: 'DFF × (initial damage + yearly damage × service life). Above 1.0 the fatigue criterion is not met.' },
        { key: 'maxRange', label: 'Largest stress range (incl. SCF)', value: round(gov.maxRange / 1e6), unit: 'MPa', explain: 'Stress range of the most severe load case in the governing element; compare with the knee of the S–N curve to see whether high-cycle or storm cycles dominate.' },
        { key: 'worstShare', label: 'Share of damage from the worst load case', value: round(worstBin.d / dTot * 100), unit: '%', explain: `Fraction of the governing element's damage caused by "${worstBin.name}". A dominant moderate sea state means the result is sensitive to its cycle count and curvature range.` },
        { key: 'nCycles', label: 'Cycles per year in the spectrum', value: round(bins.reduce((a, q) => a + q.cycles, 0)), unit: '-', explain: 'Total number of counted cycles per year, table plus rainflow history. About 5 million a year is typical for wave-frequency loading.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: `Damage by load case — ${gov.name}`, categories: gb.map((q) => q.name), series: [{ name: 'Damage per year [-]', y: gb.map((q) => q.d) }],
          explain: 'Which part of the long-term spectrum consumes the fatigue life of the governing element. The tallest bar is the load case to analyse most carefully (and the one whose reduction buys the most life).' },
        { type: 'bar', title: `Stress-range components — ${gov.name}`, categories: gb.map((q) => q.name), series: [{ name: 'Tension share [MPa]', y: gb.map((q) => round(q.tensionPart / 1e6)) }, { name: 'Own-axis bending [MPa]', y: gb.map((q) => round(q.bendingPart / 1e6)) }, { name: 'Friction [MPa]', y: gb.map((q) => round(q.frictionPart / 1e6)) }],
          explain: 'Make-up of the stress range per load case: tension variation, bending of the element about its own axis and stick–slip friction. Bending-dominated ranges are cured by a longer bend stiffener, friction-dominated ones by lower contact pressure.' },
        { type: 'line', title: 'S–N curve and load cases (log–log)', x: { label: 'log10 cycles', unit: '-' }, y: { label: 'log10 stress range', unit: 'log10 MPa' },
          series: [{ name: `S–N (${gov.tube ? inp.snTube : inp.snWire}, ${inp.snEnv})`, x: sGrid.map((v) => Math.log10(snCycles({ stressRange: v, curve: gov.curve }).cycles)), y: sGrid.map((v) => Math.log10(v / 1e6)) },
            { name: 'Load cases (cycles per year)', x: gb.map((q) => Math.log10(Math.max(q.cycles, 1))), y: gb.map((q) => Math.log10(Math.max(q.total * gov.scf, 1) / 1e6)) }],
          explain: 'Design S–N curve of the governing element with each load case plotted at its yearly cycle count. Points close to the curve consume life quickly; the horizontal gap is the number of years that case alone would take to fail.' },
        { type: 'line', title: 'Damage accumulation', x: { label: 'Service time', unit: 'years' }, y: { label: 'Miner sum', unit: '-' },
          series: [{ name: gov.name, x: yrs, y: yrs.map((y) => inp.D0 + gov.rate * y) }, { name: `Allowable 1/DFF = ${fmt(1 / inp.dff)}`, x: [0, yrs[yrs.length - 1]], y: [1 / inp.dff, 1 / inp.dff], dash: true }, { name: 'Design life', x: [inp.life, inp.life], y: [0, Math.max(1 / inp.dff, inp.D0 + gov.rate * inp.life)], dash: true }],
          explain: 'Linear growth of the Miner sum from the initial damage. The crossing with the allowable line is the end of the fatigue design life; it must lie to the right of the design-life marker.' },
      ],
      tables: [
        { title: 'Fatigue by element (critical locations)', explain: 'One row per steel tube or wire type: yearly damage, calculated life, design life with DFF and utilisation at the service life. The first-ranked element is the critical fatigue location in the cross-section.',
          columns: [{ key: 'name', label: 'Element' }, { key: 'sn', label: 'S–N curve' }, { key: 'maxRange', label: 'Max range', unit: 'MPa' }, { key: 'rate', label: 'Damage / year', unit: '1/yr' }, { key: 'life', label: 'Life (no DFF)', unit: 'yr' }, { key: 'lifeD', label: 'Design life', unit: 'yr' }, { key: 'uc', label: 'Utilisation', unit: '-' }],
          rows: res.slice().sort((a, c) => c.rate - a.rate).map((r) => ({ name: r.name, sn: r.tube ? inp.snTube : inp.snWire, maxRange: round(r.maxRange / 1e6), rate: round(r.rate, 3), life: round(Math.min(r.fl.characteristicLife, 1e6)), lifeD: round(Math.min(r.fl.designLife, 1e6)), uc: round(inp.dff * (inp.D0 + r.rate * inp.life)) })) },
        { title: `Load cases — ${gov.name}`, explain: 'One row per load case for the governing element: stress-range components, cycles per year, cycles to failure on the S–N curve and the resulting yearly damage.',
          columns: [{ key: 'name', label: 'Load case' }, { key: 'dS', label: 'Stress range', unit: 'MPa' }, { key: 'n', label: 'Cycles / yr' }, { key: 'N', label: 'Cycles to failure' }, { key: 'd', label: 'Damage / yr' }],
          rows: gb.map((q) => ({ name: q.name, dS: round(q.total * gov.scf / 1e6), n: round(q.cycles), N: Number.isFinite(q.N) ? round(q.N) : null, d: round(q.d, 3) })) },
      ],
      notes: [
        `Critical element: ${gov.name}; critical load case: ${worstBin.name} (${fmt(worstBin.d / dTot * 100)} % of its damage).`,
        `Fatigue design life ${fmt(lifeDff)} years with DFF ${inp.dff} against a service life of ${inp.life} years.`,
        rf ? `Rainflow counting of the history found ${rf.cycles.length} cycles; they were scaled to ${fmt(inp.histFraction * 100)} % of a year and added to the table spectrum.` : 'No tension history supplied: the table spectrum alone was used.',
      ],
      advice: [
        adv('engineering', sev(gc.uc), `${gov.name} governs with a fatigue design life of ${fmt(lifeDff)} years (utilisation ${fmt(gc.uc)}); ${worstBin.name} contributes ${fmt(worstBin.d / dTot * 100)} % of the damage, with ${fmt(worstBin.bendingPart / Math.max(worstBin.total, 1) * 100)} % of its stress range from bending.`,
          gc.uc > 1 ? `Reduce the stress range by about ${fmt((1 - gc.uc ** (-1 / gov.curve.m)) * 100)} % (S–N slope ${gov.curve.m}): lengthen the bend stiffener to cut curvature ranges, re-tune the configuration to lower hang-off angle variation, or select a tube with a better qualified weld class.`
            : 'Replace the stand-in S–N classes with the manufacturer-qualified tube and wire curves and repeat; keep the load spectrum under revision control because life scales with the 3rd–5th power of the ranges.',
          `${STD_UMB} §7.4; DNV-RP-C203 (2021) §2.2–2.4`),
        adv('economic', gc.uc > 1 ? 'warn' : 'info', `At the calculated rate the dynamic section would need replacing after ${fmt(lifeDff)} years on the DFF basis; the service life is ${inp.life} years.`,
          gc.uc > 1 ? `Budget ${Math.max(1, Math.ceil(inp.life / Math.max(lifeDff, 0.5)) - 1)} replacement(s) of the ${fmt(inp.lengthKm * 1e3)} m dynamic section, or invest in the curvature reduction above: a 20 % cut in stress range roughly doubles the life.` : `The reserve (${fmt(lifeDff / inp.life)} × the service life) supports a life-extension case later; record the as-installed spectrum so the remaining life can be re-assessed rather than replacing on age.`,
          'Whole-life costing ISO 15686-5; life-extension assessment NORSOK U-009 / DNV-RP-G101 risk-based approach',
          `One replacement of the dynamic section re-manufactures ≈ ${fmt(steelT)} t of steel ≈ ${fmt(steelT * CO2_STEEL)} t CO₂e (screening ${CO2_STEEL} t/t) and needs a lay campaign: ${vesselCo2Text(6)}`),
      ],
      publish: { 'life.fatigueYears': lifeDff, 'umbilical.fatigueDamageRate': gov.rate },
      quantities: { steel_t: steelT },
    };
  },
  verify: [
    // Δσ = 40 kN / 238.51 mm² = 167.71 MPa; DNV-RP-C203 B1 in air: log N = 15.117 − 4·log(167.71) → N = 1.6549e6;
    // 1e5 cycles/yr → d = 0.060425 /yr, life = 16.549 yr, with DFF 10 → 1.6549 yr.
    { name: 'Single tube, constant-amplitude tension: N = a·Δσ⁻ᵐ, d = n/N',
      source: 'DNV-RP-C203 (2021) Table 2-1 curve B1 (m = 4, log a = 15.117); hand calculation',
      inputs: { sheathT: 0, layers: [], snTube: 'B1', snEnv: 'air', scf: 1, D0: 0, spectrum: [{ name: 'CA', meanT: 50, dT: 40, dK: 0, cycles: 1e5 }],
        components: [{ name: 'Tube', kind: 'tube', layer: 0, count: 1, od: 31.8, wall: 2.6, E: 200, rho: 7800, fy: 550, p: 0, rhoF: 0 }] },
      expect: [{ kpi: 'damageRate', value: 0.060425, tol: 0.002 }, { kpi: 'lifeDesign', value: 1.6549, tol: 0.002 }, { kpi: 'maxRange', value: 167.71, tol: 0.002 }] },
    // Reproducible stress history: 0, 50, 0, 50, … kN (20 reversals pairs → 10 full cycles of range 50 kN in 100 s).
    // Scaled to a year at fraction 1: n = 10 × 31 557 600 / 100 = 3.15576e6. Δσ = 50e3/238.51e-6 = 209.63 MPa →
    // N = 10^15.117/209.63⁴ = 6.7787e5 → d = 4.6554 /yr.
    { name: 'Rainflow of a constant-amplitude history reproduces n/N',
      source: 'ASTM E1049 rainflow (constant amplitude → one cycle per period); DNV-RP-C203 B1 air',
      inputs: { sheathT: 0, layers: [], snTube: 'B1', snEnv: 'air', scf: 1, D0: 0, spectrum: [], histFraction: 1, kPerT: 0,
        history: { x: Array.from({ length: 21 }, (_, i) => i * 5), y: Array.from({ length: 21 }, (_, i) => (i % 2 ? 50 : 0)) },
        components: [{ name: 'Tube', kind: 'tube', layer: 0, count: 1, od: 31.8, wall: 2.6, E: 200, rho: 7800, fy: 550, p: 0, rhoF: 0 }] },
      expect: [{ kpi: 'damageRate', value: 4.6554, tol: 0.003 }] },
  ],
};
