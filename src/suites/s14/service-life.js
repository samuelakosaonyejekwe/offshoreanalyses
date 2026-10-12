// Flexible pipe fatigue and service-life assessment: rainflow counting of a curvature history, tensile-armour
// fatigue, pressure-armour pressure cycling, wear and broken-wire effects, accumulated damage and remaining life,
// with the annulus volume test as the inspection input.
import { rainflowFourPoint } from '../../engine/fatigue/rainflow.js';
import { axisymmetricResponse, bendingProperties } from '../../engine/xsection/flexible.js';
import { armourFatigueCoefficients, armourFatigueDamage, pressureCycleStressRange, minerDamage, wireSlipPerCycle, archardWear } from '../../engine/xsection/flexfatigue.js';
import { annulusFreeVolume, annulusVolumeTest } from '../../engine/xsection/annulus.js';
import { stackInputs, stackFrom, loadInputs, loadsFrom, muInput, kpi, chk, adv, plot, bar, tbl, col, run, need, fmt, round, thin, LIMIT_AXI, GENERIC, STD_17J, STD_17B, CO2_MGO, MGO_T_PER_DAY, BAR } from './_common.js';
import { snInputs, snFrom, cap } from './armour-fatigue.js';

const YEAR = 31557600;
// Deterministic irregular curvature record used as the default history (three superposed wave components).
const DEF_T = Array.from({ length: 601 }, (_, i) => i);
const DEF_K = DEF_T.map((t) => Number((0.0018 * Math.sin(2 * Math.PI * t / 9) + 0.001 * Math.sin(2 * Math.PI * t / 13.7 + 1) + 0.0005 * Math.sin(2 * Math.PI * t / 5.3 + 2)).toFixed(7)));

export default {
  id: 's14-service-life',
  suite: 14,
  title: 'Fatigue and service-life assessment',
  summary: 'Combines armour fatigue from a counted curvature history, pressure-armour pressure cycling, wear, broken wires and inspection findings into accumulated damage and remaining service life.',
  covers: [15],
  fidelity: 'analytical',
  standards: ['API Spec 17J (4th ed.) §5.3.4 (service life)', 'API RP 17B (5th ed.) §8, §13 (integrity management)', 'ASTM E1049-85 (rainflow counting)', 'DNV-RP-G101 / API RP 580 (risk-based inspection)'],
  models: ['Rainflow counting (four-point rule, repeating history)', 'Palmgren–Miner damage', 'Armour-wire fatigue: friction stick–slip and local wire bending', 'Pressure-armour pressure-cycle fatigue', 'Archard wear and net-section stress increase', 'Inspection-based condition update (annulus volume test, broken wires)'],
  limits: LIMIT_AXI + ' The curvature history is taken as stationary and repeated all year at one mean tension; long-term scatter needs one run per sea state (or the load-bin module). No built-in sweet/sour S–N curves: the generic curve is unqualified and the environment knock-down must come from qualified tests. Polymer ageing (e.g. PA11 hydrolysis) and corrosion wastage are not modelled. ' + GENERIC,
  publishes: ['life.fatigueYears', 'flexible.remainingLife'],
  searchTerms: 'flexible riser service life remaining life fatigue integrity management annulus test',
  inputs: [
    ...stackInputs(),
    ...loadInputs({ depth: 30 }).map((i) => (i.key === 'depth' ? { ...i, bind: undefined, help: 'Depth of the critical section (hang-off region ≈ 0–50 m).' } : i)),
    muInput(),
    { key: 'history', label: 'Curvature time history at the critical section', type: 'series', default: { x: DEF_T, y: DEF_K }, xLabel: 'Time', xUnit: 's', yLabel: 'Curvature', yUnit: '1/m', group: 'Global loading', role: 'bc',
      help: 'Representative record from global dynamic analysis or monitoring. It is rainflow counted and repeated over the year.' },
    { key: 'occurrence', label: 'Fraction of the year this record represents', unit: '%', type: 'number', default: 100, min: 0.1, max: 100, group: 'Global loading', role: 'bc', help: 'Use 100 % for a record that represents the long-term average sea state; run one case per sea state and add the damage rates for a scatter diagram.' },
    { key: 'pLow', label: 'Pressure cycle: low', unit: 'bar', type: 'number', default: 20, min: 0, max: 1500, group: 'Pressure history', role: 'operation' },
    { key: 'pHigh', label: 'Pressure cycle: high', unit: 'bar', type: 'number', default: 200, min: 0, max: 1500, group: 'Pressure history', role: 'operation' },
    { key: 'nPress', label: 'Pressure cycles per year', type: 'number', default: 24, min: 0, max: 1e6, group: 'Pressure history', role: 'operation', help: 'Shutdown / start-up cycles between the two pressures.' },
    ...snInputs({ source: 'custom' }),
    { key: 'yearsInService', label: 'Years in service', unit: 'years', type: 'number', default: 10, min: 0, max: 60, group: 'Condition and inspection', role: 'ic' },
    { key: 'priorDamage', label: 'Known accumulated Miner damage (0 = compute from years in service)', unit: '-', type: 'number', default: 0, min: 0, max: 1, group: 'Condition and inspection', role: 'ic', help: 'Unfactored damage already consumed, e.g. from a previous assessment or a different earlier service.' },
    { key: 'brokenWires', label: 'Broken wires in the outer tensile armour', type: 'int', default: 0, min: 0, max: 40, group: 'Condition and inspection', role: 'ic', help: 'From inspection; the remaining wires carry proportionally more stress.' },
    { key: 'wearK', label: 'Wire wear coefficient (0 = intact anti-wear tape)', unit: 'm²/N', type: 'number', default: 0, min: 0, max: 1e-12, group: 'Condition and inspection', role: 'material', help: 'Archard coefficient for wire wear once the tape is lost; from wear tests.' },
    { key: 'testDone', label: 'Annulus volume test available', type: 'bool', default: true, group: 'Condition and inspection', role: 'ic' },
    { key: 'testLength', label: 'Annulus length tested', unit: 'm', type: 'number', default: 1500, min: 1, max: 50000, group: 'Condition and inspection', role: 'geometry', showIf: { key: 'testDone', in: [true] } },
    { key: 'testP1', label: 'Test: initial annulus pressure (absolute)', unit: 'kPa', type: 'number', default: 101.3, min: 1, max: 2000, group: 'Condition and inspection', role: 'ic', showIf: { key: 'testDone', in: [true] } },
    { key: 'testP2', label: 'Test: final annulus pressure (absolute)', unit: 'kPa', type: 'number', default: 150, min: 1, max: 3000, group: 'Condition and inspection', role: 'ic', showIf: { key: 'testDone', in: [true] } },
    { key: 'testVol', label: 'Test: nitrogen injected (standard litres)', unit: 'L', type: 'number', default: 720, min: 0.1, max: 1e6, group: 'Condition and inspection', role: 'ic', showIf: { key: 'testDone', in: [true] }, help: 'Gas volume at standard conditions that raised the annulus pressure from the initial to the final value.' },
    { key: 'reqRemaining', label: 'Remaining life required', unit: 'years', type: 'number', default: 15, min: 0, max: 60, group: 'Acceptance criteria', role: 'criteria' },
  ],
  compute(inp, ctx) {
    const S = stackFrom(inp, ctx);
    need(S.idx.ta.length > 0, 'The stack has no tensile armour layer.');
    const hx = inp.history?.x || [], hy = inp.history?.y || [];
    need(hx.length >= 4 && hx.length === hy.length, 'The curvature history needs at least four points with matching time and curvature columns.');
    need(hy.every(Number.isFinite) && hx.every(Number.isFinite), 'The curvature history contains non-numeric values.');
    const duration = hx[hx.length - 1] - hx[0];
    need(duration > 0, 'The curvature history must span a positive time.');
    need(inp.pHigh >= inp.pLow, 'Pressure cycle: the high pressure must not be below the low pressure.');
    const rf = rainflowFourPoint(hy, { residue: 'repeat' });
    let cyc = rf.cycles.filter((c) => c.range > 0);
    need(cyc.length > 0, 'The curvature history contains no cycles (constant signal).');
    const scale = YEAR * inp.occurrence / 100 / duration;
    if (cyc.length > 400) { // group into 60 range classes at the upper edge (conservative)
      const rMax = Math.max(...cyc.map((c) => c.range)), nb = 60, acc = new Array(nb).fill(0);
      for (const c of cyc) acc[Math.min(nb - 1, Math.floor(c.range / rMax * nb))] += c.count;
      cyc = acc.map((n, k) => ({ range: (k + 1) / nb * rMax, count: n })).filter((c) => c.count > 0);
    }
    const loads = loadsFrom(inp); const Tm = loads.tension; delete loads.tension;
    const bins = cyc.map((c) => ({ tension: Tm, curvatureRange: c.range, cycles: c.count * scale }));
    const cyclesPerYear = bins.reduce((a, b) => a + b.cycles, 0);
    const resp = run(() => axisymmetricResponse(S.stack, { ...loads, tension: Tm }), 'Axisymmetric solution');
    const bend = bendingProperties(S.stack, { response: resp, mu: inp.mu });
    const outer = S.idx.ta[S.idx.ta.length - 1];
    const per = S.idx.ta.map((i, k) => {
      const L = S.stack.layers[i];
      need(L.wireWidth > 0, `${L.name}: enter the wire width in the layer table.`);
      const n = L.n ?? 1;
      need(i !== outer || inp.brokenWires < n, `Broken wires (${inp.brokenWires}) must be fewer than the ${n} wires of the outer armour.`);
      const broken = i === outer ? n / (n - inp.brokenWires) : 1;
      // wire wear: sliding accumulated by the cycles that exceed the slip-onset curvature
      const kcr = bend.layers[k].kappaCr ?? 0;
      let sliding = 0;
      for (const b of bins) if (b.curvatureRange / 2 > kcr) sliding += b.cycles * wireSlipPerCycle({ R: L.R, alpha: L.alpha, kappaRange: b.curvatureRange }).slidingPerCycle;
      const q = Math.max(resp.layers[i].contactIn, resp.layers[i].contactOut, 0) / L.fill;
      const wear = archardWear({ wearCoefficient: inp.wearK, contactPressure: q, slidingPerCycle: sliding, cycles: Math.max(inp.yearsInService, 0), thickness: L.wireThickness ?? L.t });
      need(!wear.wornThrough, `${L.name}: the wear coefficient entered wears the wire through within the years in service — check the value.`);
      const sn = snFrom(inp, L.uts);
      const co = run(() => armourFatigueCoefficients(S.stack, i, { tensionLow: 0.8 * Tm - 1e3, tensionHigh: Tm + 1e3, mu: inp.mu, loads }), 'Armour stress coefficients');
      const factor = inp.scf * broken * wear.stressFactor;
      const meanTop = (co.meanStress0 + co.meanStressPerTension * Tm) * broken * wear.stressFactor;
      need(inp.meanMethod === 'none' || meanTop < 0.98 * L.uts, `${L.name}: the mean wire stress (${fmt(meanTop / 1e6)} MPa, including broken-wire and wear factors) reaches the wire tensile strength — the layer fails statically, so a fatigue life is meaningless. Check tension, pressure and the number of broken wires.`);
      const d = run(() => armourFatigueDamage({ wire: co.wire, bins, meanStress0: co.meanStress0 * broken * wear.stressFactor, meanStressPerTension: co.meanStressPerTension * broken * wear.stressFactor, friction0: co.friction0, frictionPerTension: co.frictionPerTension,
        curve: sn.curve, uts: L.uts, meanStressMethod: inp.meanMethod, scf: factor, dff: inp.dff, exposureYears: 1 }), `${L.name} fatigue`);
      return { L, sn, d, factor, wear, broken };
    });
    const ta = per.reduce((a, c) => (c.d.damagePerYear > a.d.damagePerYear ? c : a), per[0]);
    let paRate = 0, paRange = 0;
    if (S.idx.pa >= 0 && inp.nPress > 0 && inp.pHigh > inp.pLow) {
      const P = S.stack.layers[S.idx.pa];
      const pc = run(() => pressureCycleStressRange(S.stack, S.idx.pa, { pLow: inp.pLow * BAR, pHigh: inp.pHigh * BAR, loads: { ...loads, tension: Tm } }), 'Pressure-cycle stress range');
      paRange = pc.stressRange;
      paRate = minerDamage({ bins: [{ stressRange: pc.stressRange, cycles: inp.nPress, meanStress: pc.meanStress }], curve: snFrom(inp, P.uts).curve, uts: P.uts, meanStressMethod: inp.meanMethod, scf: 1, dff: inp.dff, exposureYears: 1 }).damagePerYear;
    }
    const rate = Math.max(ta.d.damagePerYear, paRate);
    const govTa = ta.d.damagePerYear >= paRate;
    const done = inp.priorDamage > 0 ? inp.priorDamage : rate * inp.yearsInService;
    const allow = 1 / inp.dff;
    const remaining = rate > 0 ? cap(Math.max(0, (allow - done) / rate)) : (done < allow ? 1e5 : 0);
    const total = cap(inp.yearsInService + remaining);
    let flood = 0, vFree = 0, vDesign = 0;
    if (inp.testDone) {
      need(inp.testP2 > inp.testP1, 'Annulus test: the final pressure must exceed the initial pressure.');
      const between = S.stack.layers.filter((L) => L.index > S.stack.sealInner && L.index < S.stack.sealOuter);
      need(between.length > 0, 'There is no annulus between two sheaths to test.');
      vDesign = annulusFreeVolume({ layers: between.map((L) => ({ rIn: L.rIn, rOut: L.rOut, fill: L.kind === 'helical' ? Math.min(1, L.nA / Math.cos(L.alpha) / (Math.PI * (L.rOut ** 2 - L.rIn ** 2))) : 1 })), length: inp.testLength }).volume;
      const t = run(() => annulusVolumeTest({ method: 'pressure', p1: inp.testP1 * 1e3, p2: inp.testP2 * 1e3, injectedStdVolume: inp.testVol / 1e3, designVolume: vDesign }), 'Annulus volume test');
      vFree = t.freeVolume; flood = Math.min(1, Math.max(0, t.floodedFraction ?? 0));
    }
    const dryCurve = inp.snSource === 'generic' && inp.env === 'air';
    const envMismatch = flood > 0.05 && dryCurve;
    const interval = Math.max(0.5, Math.min(5, remaining / 2));
    const checks = [
      chk('Accumulated damage × DFF at the required remaining life', (done + rate * inp.reqRemaining) * inp.dff, 1, '-', `${STD_17J} §5.3.4 (DFF ${inp.dff}); Palmgren–Miner`,
        'Guards against fatigue failure of the governing metallic layer before the end of the required remaining service. Above 1.0 the factored life is used up early.'),
      chk('Tensile armour: damage per year × DFF × total required life', ta.d.damagePerYear * inp.dff * (inp.yearsInService + inp.reqRemaining), 1, '-', `${STD_17J} §5.3.4; ${STD_17B} §8`,
        'Armour-wire fatigue from bending cycles at the critical section, including stress concentration, broken-wire and wear factors.'),
      chk('Pressure armour: damage per year × DFF × total required life', paRate * inp.dff * (inp.yearsInService + inp.reqRemaining), 1, '-', `${STD_17B} §8 (pressure cycling of the pressure armour)`,
        'Fatigue of the pressure-armour profile from shutdown and start-up pressure cycles.'),
    ];
    const order = ta.d.bins.map((b, k) => ({ k, d: b.damage })).sort((a, b) => b.d - a.d).slice(0, 12);
    const [tx, ty] = thin([hx, hy], 600);
    const yrs = [0, inp.yearsInService, inp.yearsInService + remaining > 1e4 ? inp.yearsInService + 2 * inp.reqRemaining : inp.yearsInService + Math.max(remaining, inp.reqRemaining)];
    const dmgAt = (y) => (y <= inp.yearsInService ? done * (inp.yearsInService > 0 ? y / inp.yearsInService : 0) : done + rate * (y - inp.yearsInService)) * inp.dff;
    const ok = remaining >= inp.reqRemaining;
    const rangeHist = cyc.map((c) => c.range);
    return {
      kpis: [
        kpi('remainingLife', 'Remaining service life (with DFF)', remaining, 'years', `Years until accumulated damage times the design fatigue factor of ${inp.dff} reaches 1.0 for the governing layer. 100 000 means no damaging cycles.`, inp.reqRemaining, 'min'),
        kpi('totalLife', 'Total factored service life', total, 'years', 'Years already in service plus the remaining life.'),
        kpi('damageRate', 'Damage per year (governing layer)', rate, '1/yr', `Unfactored Miner damage per year of the ${govTa ? ta.L.name : 'pressure armour'}.`),
        kpi('damageToDate', 'Accumulated damage to date', done, '-', 'Unfactored Miner sum consumed so far (entered, or damage rate × years in service). The factored limit is 1/DFF.'),
        kpi('lifeTA', 'Tensile armour factored life', cap(ta.d.lifeWithDff), 'years', 'Fatigue life of the governing tensile armour layer divided by the design fatigue factor.'),
        kpi('lifePA', 'Pressure armour factored life', paRate > 0 ? cap(1 / (paRate * inp.dff)) : 1e5, 'years', 'Fatigue life of the pressure armour under the pressure cycles, divided by the design fatigue factor.'),
        kpi('cyclesPerYear', 'Bending cycles per year', cyclesPerYear, '-', 'Rainflow-counted cycles of the record scaled to one year.'),
        kpi('maxRange', 'Largest counted curvature range', Math.max(...rangeHist), '1/m', 'Largest peak-to-peak curvature cycle found by rainflow counting.'),
        kpi('stressFactor', 'Stress factor from wear and broken wires', ta.factor / inp.scf, '-', 'Increase of wire stress from lost section (wear) and broken wires in the governing layer; 1.0 for an intact layer.'),
        kpi('paRange', 'Pressure-armour stress range per pressure cycle', paRange / 1e6, 'MPa', 'Change of pressure-armour wire stress between the low and high bore pressure.'),
        kpi('floodedFraction', 'Annulus flooded fraction from the volume test', flood * 100, '%', 'Share of the design free volume that the test gas could not occupy, i.e. taken up by water. Above a few per cent the annulus is wet.'),
        kpi('criticalMode', 'Critical damage location (1 = tensile armour, 2 = pressure armour)', govTa ? 1 : 2, '-', 'Which metallic layer consumes fatigue life fastest and therefore sets the remaining life.'),
        kpi('inspectionInterval', 'Suggested inspection interval', interval, 'years', 'Half the remaining life, capped at 5 years and not below 6 months — a risk-based starting point for annulus testing and visual inspection.'),
      ],
      checks,
      plots: [
        plot('line', 'Curvature history', 'Time', 's', 'Curvature', '1/m', [{ name: 'Curvature', x: tx.map((v) => round(v)), y: ty.map((v) => round(v)) }],
          'The record that is rainflow counted. Check that it is long enough to contain the large, rare cycles as well as the many small ones.'),
        plot('hist', 'Rainflow range distribution', 'Curvature range', '1/m', 'Cycles', '-', [], 'Histogram of the counted curvature ranges. The few large ranges on the right usually produce most of the damage.', { values: rangeHist.map((v) => round(v)), bins: 20 }),
        bar('Damage by range class (largest contributors)', order.map((o) => fmt(bins[o.k].curvatureRange)), 'Damage per year', '1/yr', [{ name: ta.L.name, y: order.map((o) => round(o.d, 4)) }],
          'The curvature ranges (1/m) that contribute most to yearly damage of the governing armour layer, largest first.'),
        plot('line', 'Factored damage over time', 'Years in service', 'years', 'Damage × DFF', '-', [
          { name: 'Damage × DFF', x: yrs.map((y) => round(y)), y: yrs.map((y) => round(dmgAt(y))) }, { name: 'Limit', x: [0, round(yrs[2])], y: [1, 1], dash: true }, { name: 'Today', x: [inp.yearsInService], y: [round(done * inp.dff)] }],
        'Accumulated factored damage from first service. The crossing with the limit line marks the end of the factored fatigue life.'),
      ],
      tables: [tbl('Failure-mode utilisation and remaining life', [col('mode', 'Layer / mode'), col('rate', 'Damage per year', '1/yr'), col('life', 'Factored life', 'years'), col('factor', 'Stress factor', '-'), col('sn', 'S–N basis')], [
        ...per.map((p) => ({ mode: `${p.L.name}: bending fatigue`, rate: round(p.d.damagePerYear, 4), life: round(cap(p.d.lifeWithDff)), factor: round(p.factor, 4), sn: p.sn.label })),
        { mode: 'Pressure armour: pressure cycling', rate: round(paRate, 4), life: round(paRate > 0 ? cap(1 / (paRate * inp.dff)) : 1e5), factor: 1, sn: S.idx.pa >= 0 ? snFrom(inp, S.stack.layers[S.idx.pa].uts).label : 'no pressure armour' },
      ], 'One row per metallic layer and damage mechanism. The shortest factored life is the critical location and sets the remaining service life.')],
      notes: [`${govTa ? ta.L.name : 'Pressure armour'} is critical; ${fmt(cyclesPerYear)} cycles per year counted from a ${fmt(duration)} s record.`,
        inp.testDone ? `Annulus test: free volume ${fmt(vFree * 1e3)} L against ${fmt(vDesign * 1e3)} L design → ${fmt(flood * 100)} % flooded.` : 'No annulus volume test entered: the annulus condition is assumed, not verified.',
        envMismatch ? 'The annulus test indicates water in the annulus, but the fatigue curve is for dry air — the result is unconservative until a qualified wet-annulus curve or knock-down is entered.' : `S–N basis: ${ta.sn.label}.`, GENERIC],
      advice: [
        adv('engineering', envMismatch ? 'warn' : (ok ? 'ok' : 'fail'), `Remaining factored life ${remaining >= 1e5 ? 'is not limited by fatigue' : fmt(remaining) + ' years'} against ${inp.reqRemaining} years required; ${fmt(done * inp.dff * 100)} % of the factored life is used after ${fmt(inp.yearsInService)} years.`,
          envMismatch ? `The volume test shows ${fmt(flood * 100)} % of the annulus is flooded: re-run with the environment set to sweet or sour and the knock-down from qualified wire tests, and plan a sheath repair.` : (ok ? 'Keep the basis valid: repeat the annulus volume test at the suggested interval and update the curvature history from monitoring data.' : (done >= allow ? 'The factored fatigue life is already used up on this basis: verify the load history and S–N data, inspect now, and plan replacement or re-termination.' : `Reduce the damage rate by a factor ${fmt(rate * inp.reqRemaining / (allow - done))} (smaller curvature at the stiffener, dry annulus) or plan replacement in year ${fmt(total)}.`)),
          `${STD_17J} §5.3.4; ${STD_17B} §13 (integrity management); ASTM E1049-85 (2017)`),
        adv('economic', ok ? 'info' : 'warn', `On this basis the riser can serve ${remaining >= 1e5 ? 'beyond any planned life' : 'until year ' + fmt(total)}; suggested inspection interval ${fmt(interval)} years.`,
          ok ? `Life extension of ${remaining >= 1e5 ? 'the field life' : fmt(remaining - inp.reqRemaining) + ' years beyond the requirement'} is supportable: an annulus test every ${fmt(interval)} years (a few hours of platform time) protects a replacement campaign of 10–20 vessel days.` : `Budget a replacement or re-termination in year ${fmt(total)} — roughly 10–20 vessel days plus the pipe — or invest in curvature reduction now; compare on a net-present-cost basis.`,
          'Risk-based inspection per DNV-RP-G101 (2021) / API RP 580; life-extension economics per ISO 15686-5',
          `Each year of extended service avoids manufacturing a replacement pipe and a campaign of about 15 vessel days ≈ ${fmt(15 * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening: ${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂ per t fuel).`),
      ],
      publish: { 'life.fatigueYears': total, 'flexible.remainingLife': remaining },
    };
  },
  verify: [{
    // 50 sine cycles of curvature range 0.02 1/m in 500 s → 50·31 557 600/500 = 3.1558e6 cycles per year. With friction
    // off the wire stress range is local bending only: E·(b/2)·(1 + sin²α)·cosα·Δκ = 27.04 MPa (wire 12 × 4 mm, α = 35°).
    // S–N log a = 12 (MPa), m = 3 → N = 1e12/27.04³ = 5.058e7; damage rate = 3.1558e6/5.058e7 = 0.06239 per year.
    // New pipe, DFF 10: remaining life = (1/10)/0.06239 = 1.603 years.
    name: 'Constant-amplitude history: rainflow count, Miner rate and remaining life',
    source: 'ASTM E1049 (a sine wave of k periods counts as k cycles of its range); Basquin/Miner hand calculation; Sævik (2011) wire bending',
    inputs: { mu: 0, snSource: 'custom', snLogA: 12, snM: 3, meanMethod: 'none', scf: 1, nPress: 0, yearsInService: 0, priorDamage: 0, brokenWires: 0, wearK: 0, occurrence: 100, dff: 10, testDone: false,
      history: { x: Array.from({ length: 1001 }, (_, i) => i * 0.5), y: Array.from({ length: 1001 }, (_, i) => 0.01 * Math.sin(2 * Math.PI * (i * 0.5) / 10 + Math.PI / 2)) } },
    expect: [{ kpi: 'cyclesPerYear', value: 3.1558e6, tol: 0.015 }, { kpi: 'damageRate', value: 0.06239, tol: 0.02 }, { kpi: 'remainingLife', value: 1.603, tol: 0.02 }],
  }],
};
