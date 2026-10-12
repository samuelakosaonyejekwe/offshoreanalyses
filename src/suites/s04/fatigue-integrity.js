// Structural fatigue and integrity: S–N damage from environmental and operational cycles with corrosion and prior damage.
import { snCurve, snConstants, cyclesFromConstants, SN_CURVE_NAMES } from '../../engine/fatigue/sn.js';
import { minerWeibull, minerSum, fatigueLife } from '../../engine/fatigue/damage.js';
import { corrosionWastage } from '../../engine/fatigue/degradation.js';
import { need, run, fmt, round, adv, check, vesselCo2Text, CONFIRM } from './_common.js';

const DEFAULT_BLOCKS = [
  { name: 'Start-up / shut-down thermal cycles', range: 110, perYear: 24 },
  { name: 'Pressure fluctuation and slugging', range: 18, perYear: 200000 },
  { name: 'Flow-induced vibration of small-bore piping supports', range: 9, perYear: 3000000 },
];

export default {
  id: 's04-fatigue-integrity', suite: 4, title: 'Structural fatigue and integrity',
  summary: 'Miner damage at a welded detail from wave-induced and operational stress cycles, with corrosion wastage, existing damage and design fatigue factor, giving fatigue life and remaining service life.',
  covers: [14], fidelity: 'analytical',
  standards: ['DNV-RP-C203 (2021) §2.4, §5.1', 'DNV-ST-0126 / DNV-OS-C101 (design fatigue factors)', 'ISO 19902:2020 §16', 'DNV-RP-B401 (2021) (cathodic protection context)'],
  models: ['S–N fatigue assessment (DNV-RP-C203 curves)', 'Palmgren–Miner cumulative damage', 'Weibull long-term stress-range distribution (closed form)', 'Linear corrosion wastage after coating breakdown', 'Thickness and stress-concentration correction'],
  limits: 'Nominal-stress S–N approach at one detail; the long-term wave-induced stress ranges follow a two-parameter Weibull distribution given by its largest range, and operational cycles are constant-amplitude blocks. Corrosion raises stress in inverse proportion to remaining thickness. Crack growth by fracture mechanics, variable-amplitude sequence effects and hot-spot finite-element stresses are outside the module.',
  publishes: ['life.fatigueYears', 'subsea.fatigueUc'],
  searchTerms: 'subsea structure fatigue S-N Miner DNV-RP-C203 remaining life',
  inputs: [
    { key: 'snName', label: 'S–N curve (detail class)', type: 'select', default: 'F', options: SN_CURVE_NAMES.map((n) => ({ value: n, label: `DNV-RP-C203 ${n}` })), group: 'Detail', role: 'material',
      help: 'Fillet-welded attachments F–F3; full-penetration butt welds D–E; cruciform and load-carrying fillets F1–W. Take the class from the welding and inspection specification.' },
    { key: 'snEnv', label: 'Environment', type: 'select', default: 'seawater-cp', group: 'Detail', role: 'material', options: [{ value: 'air', label: 'Air' }, { value: 'seawater-cp', label: 'Sea water with cathodic protection' }, { value: 'free-corrosion', label: 'Free corrosion' }] },
    { key: 'scf', label: 'Stress concentration factor', unit: '-', type: 'number', default: 1.3, min: 1, max: 10, group: 'Detail', role: 'geometry', uncertain: { dist: 'lognormal', cov: 0.1 } },
    { key: 'thk', label: 'Plate / wall thickness', unit: 'mm', type: 'number', default: 20, min: 4, max: 150, group: 'Detail', role: 'geometry', design: true },
    { key: 'dS0', label: 'Largest wave-induced stress range', unit: 'MPa', type: 'number', default: 55, min: 0, max: 800, group: 'Environmental cycles', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.15 },
      help: 'Nominal stress range exceeded once in the reference number of cycles, from the in-place analysis under cyclic hydrodynamic load.' },
    { key: 'n0', label: 'Reference number of cycles', unit: '-', type: 'number', default: 1e8, min: 1e3, max: 1e10, group: 'Environmental cycles', role: 'bc', help: 'Number of cycles in which the largest range occurs once; 10⁸ corresponds to about 20 years of wave action.' },
    { key: 'hW', label: 'Weibull shape parameter', unit: '-', type: 'number', default: 0.9, min: 0.4, max: 2, group: 'Environmental cycles', role: 'bc', uncertain: { dist: 'normal', cov: 0.08 }, help: 'About 0.8–1.0 for drag-dominated near-bed loading; 1.0–1.2 for inertia-dominated.' },
    { key: 'cycYr', label: 'Wave cycles per year', unit: '1/yr', type: 'number', default: 5e6, min: 1e3, max: 5e7, group: 'Environmental cycles', role: 'bc' },
    { key: 'blocks', label: 'Operational stress cycles', type: 'table', group: 'Operational cycles', role: 'bc',
      columns: [{ key: 'name', label: 'Source', type: 'text' }, { key: 'range', label: 'Stress range', unit: 'MPa', type: 'number' }, { key: 'perYear', label: 'Cycles per year', type: 'number' }],
      default: DEFAULT_BLOCKS, help: 'Constant-amplitude blocks from thermal expansion reactions, pressure cycles and vibration.' },
    { key: 'age', label: 'Present age of the structure', unit: 'yr', type: 'number', default: 0, min: 0, max: 100, group: 'Initial state', role: 'ic', help: '0 for a new design; the service years already elapsed for a life-extension assessment.' },
    { key: 'D0', label: 'Fatigue damage already accumulated', unit: '-', type: 'number', default: 0, min: 0, max: 5, group: 'Initial state', role: 'ic', help: 'Characteristic Miner sum consumed so far (fabrication, transport, installation and past service), without the design fatigue factor.' },
    { key: 'corr0', label: 'Wall loss already measured', unit: 'mm', type: 'number', default: 0, min: 0, max: 50, group: 'Initial state', role: 'ic' },
    { key: 'life', label: 'Required service life from now', unit: 'yr', type: 'number', default: 25, min: 1, max: 100, group: 'Service', role: 'operation', bind: 'life.designYears' },
    { key: 'corrRate', label: 'Corrosion rate after coating / CP depletion', unit: 'mm/yr', type: 'number', default: 0.1, min: 0, max: 2, group: 'Service', role: 'environment', uncertain: { dist: 'lognormal', cov: 0.4 }, help: 'Order 0.1–0.2 mm/yr for immersed uncoated carbon steel; zero while cathodic protection is effective.' },
    { key: 'protLife', label: 'Remaining life of coating and cathodic protection', unit: 'yr', type: 'number', default: 20, min: 0, max: 100, group: 'Service', role: 'environment', help: 'Corrosion starts when the protection is exhausted; from the anode design (DNV-RP-B401).' },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 3, min: 1, max: 10, group: 'Acceptance criteria', role: 'criteria', help: `1 for inspectable non-critical details up to 10 for uninspectable critical details; 3 is common for subsea structures without planned inspection. ${CONFIRM}` },
  ],
  compute(inp) {
    // SI: MPa → Pa, mm → m.
    const curve = run(() => snCurve(inp.snName, inp.snEnv), 'S–N curve');
    const t0 = (inp.thk - inp.corr0) / 1e3;
    need(t0 > 0.002, 'The measured wall loss leaves less than 2 mm of plate.');
    const blocks = (inp.blocks || []).filter((b) => b.range > 0 && b.perYear > 0);
    need(blocks.length > 0 || inp.dS0 > 0, 'Give a wave-induced stress range or at least one operational cycle block.');
    // Annual characteristic damage for a given remaining thickness (stress ∝ 1/t for membrane action).
    const rate = (t) => {
      const sn = snConstants(curve, { thickness: t, scf: inp.scf });
      const f = t0 / t;
      const dW = inp.dS0 > 0 ? minerWeibull({ sn, shape: inp.hW, cycles: inp.n0, maxStressRange: inp.dS0 * 1e6 * f }).damage * inp.cycYr / inp.n0 : 0;
      const ms = blocks.length ? minerSum({ ranges: blocks.map((b) => b.range * 1e6 * f), counts: blocks.map((b) => b.perYear), sn }) : { damage: 0, contributions: [] };
      return { dW, dB: ms.damage, parts: Array.from(ms.contributions), total: dW + ms.damage, sn };
    };
    const thkAt = (yr) => t0 - corrosionWastage({ time: yr, model: 'linear', coatingLife: inp.protLife, rate: inp.corrRate / 1e3 }).depth;
    const r0 = rate(t0);
    need(r0.total > 0, 'All stress ranges are below the fatigue limit of the curve: no damage accumulates. Check the stress ranges.');
    // March year by year from now: cumulative damage with thickness loss; stop when DFF·D reaches 1.
    const yrs = [0], Dc = [inp.D0], tk = [t0 * 1e3];
    let D = inp.D0, remaining = Infinity, tMin = t0;
    for (let y = 1; y <= 300; y++) {
      const t = thkAt(y - 0.5);
      if (t <= 0.002) { if (!Number.isFinite(remaining)) remaining = y - 1; break; }
      const dr = rate(t).total;
      const Dn = D + dr;
      if (!Number.isFinite(remaining) && inp.dff * Dn >= 1) remaining = y - 1 + (1 / inp.dff - D) / dr;
      D = Dn;
      if (y <= Math.max(inp.life * 1.5, 30)) { yrs.push(y); Dc.push(D); tk.push(t * 1e3); }
      if (y === Math.round(inp.life)) tMin = t;
      if (Number.isFinite(remaining) && y > inp.life * 1.5) break;
    }
    if (!Number.isFinite(remaining)) remaining = 300;
    remaining = Math.max(remaining, 0);
    const iLife = Math.min(Math.round(inp.life), Dc.length - 1);
    const Dlife = Dc[iLife];
    const uc = inp.dff * Dlife;
    const fl = fatigueLife({ damageRate: r0.total, designLife: inp.life, dff: inp.dff, accumulatedDamage: inp.D0 });
    const src = [{ name: 'Wave-induced (Weibull)', d: r0.dW }, ...blocks.map((b, i) => ({ name: b.name, d: r0.parts[i] || 0 }))];
    const top = src.reduce((a, s) => (s.d > a.d ? s : a));
    // S–N curve for the chart.
    const sr = [], nn = [];
    for (let i = 0; i <= 40; i++) { const s = 5 * 10 ** (i * 0.05); const N = cyclesFromConstants(r0.sn, s * 1e6); if (Number.isFinite(N) && N > 0) { sr.push(round(s)); nn.push(N); } }
    const m = curve.m1;
    // Stress reduction (or thickness increase) that would bring utilisation to 1: damage ∝ S^m on the upper branch.
    const sFac = uc > 0 ? uc ** (-1 / m) : 1;
    const checks = [
      { name: 'Fatigue damage at end of required life (with DFF)', demand: round(uc), capacity: 1, unit: 'DFF × Miner sum', uc: Math.min(uc, 99), ref: 'DNV-RP-C203 (2021) §2.2 Eq. (2.2.1) with design fatigue factor', explain: 'Guards against fatigue crack initiation at the detail within the required life. Above 1 a through-thickness crack becomes likely before end of life.' },
      check('Required life against remaining fatigue life', inp.life, remaining, 'yr', 'DNV-RP-C203 (2021) §5; ISO 19902 §16', 'Compares the service life still required with the life left after the damage already consumed and the corrosion expected.'),
      check('Remaining thickness at end of life', 0.5 * inp.thk, tMin * 1e3, 'mm', 'Integrity screening: at least half of the nominal thickness to remain', 'Guards against corrosion consuming the section once coating and anodes are depleted, which raises all stresses and accelerates fatigue.'),
    ];
    const ok = checks.every((c) => c.uc <= 1);
    const inspInt = Math.max(Math.min(remaining / 3, 10), 1);
    const advice = [
      ok ? adv('engineering', 'ok', `Damage at ${fmt(inp.life)} years is ${fmt(Dlife)} (${fmt(uc)} with DFF ${inp.dff}); remaining fatigue life ${fmt(remaining)} years. "${top.name}" contributes ${fmt(top.d / r0.total * 100)} % of the annual damage.`,
        uc < 0.3 ? `There is reserve: the stress range could rise by a factor of ${fmt(sFac)} (or the detail class drop one step) before the limit is reached. Confirm the SCF by local analysis and the cycle counts from operating data.` : 'Confirm the SCF by local finite-element analysis and keep cathodic protection effective to the end of life; the result is sensitive to both.',
        'DNV-RP-C203 (2021) §2.4, §5.1; ISO 19902:2020 §16',
        `A structure that reaches its design life without repair avoids an intervention campaign: ${vesselCo2Text(3)}`)
        : adv('engineering', 'fail', `Fatigue utilisation ${fmt(uc)} with DFF ${inp.dff}; remaining life ${fmt(remaining)} years against ${fmt(inp.life)} required. "${top.name}" dominates (${fmt(top.d / r0.total * 100)} % of annual damage).`,
          `Reduce the stress range at the detail by a factor of ${fmt(sFac)} — for example thickness ${fmt(inp.thk / sFac)} mm instead of ${fmt(inp.thk)} mm, a lower SCF by weld profiling or grinding (DNV-RP-C203 §7), or a better detail class — or reduce the "${top.name}" cycles at source.`,
          'DNV-RP-C203 (2021) §2.4, §5.1, §7 (improvement methods); ISO 19902:2020 §16',
          `Weld toe grinding roughly doubles to triples life with no added steel; the thickness route adds about ${fmt((1 / sFac - 1) * 100)} % to the plate mass at the detail, at about 1.9 t CO₂e per tonne of added steel (screening factor).`),
      adv('economic', ok ? 'info' : 'warn', `With ${fmt(remaining)} years of remaining life and DFF ${inp.dff}, a first close visual / CP inspection at about ${fmt(inspInt)} years is consistent with the predicted damage growth.`,
        ok ? 'Use the reserve to justify a lower inspection frequency or a life extension rather than adding steel; each avoided ROV inspection campaign saves about 2 vessel days.' : 'Weigh improvement at fabrication (hours of grinding) against in-service inspection at a higher frequency or a repair clamp later (several vessel days plus production deferral).',
        'Risk-based inspection per DNV-RP-G101 / DNV-RP-C210 (probabilistic inspection planning); whole-life costing per ISO 15686-5'),
    ];
    return {
      kpis: [
        { key: 'dYear', label: 'Annual fatigue damage (present)', unit: '1/yr', value: r0.total, explain: 'Characteristic Miner damage accumulated per year at the present thickness, without the design fatigue factor.' },
        { key: 'lifeChar', label: 'Characteristic fatigue life', unit: 'yr', value: Math.min(fl.characteristicLife, 1e6), explain: 'Years to a Miner sum of 1 at the present damage rate, from new and without the design fatigue factor (capped at 10⁶).' },
        { key: 'lifeDesign', label: 'Design fatigue life', unit: 'yr', value: Math.min(fl.designLife, 1e6), explain: 'Characteristic life divided by the design fatigue factor: the life that may be claimed in design.' },
        { key: 'remaining', label: 'Remaining service life', unit: 'yr', value: remaining, limit: inp.life, sense: 'min', explain: 'Years from now until DFF × damage reaches 1, including damage already consumed and the stress increase from corrosion (capped at 300).' },
        { key: 'Dlife', label: 'Damage at end of required life', unit: '-', value: Dlife, explain: 'Characteristic Miner sum at the end of the required life, including the starting damage.' },
        { key: 'uc', label: 'Fatigue utilisation', unit: '-', value: uc, limit: 1, sense: 'max', explain: 'Design fatigue factor times the Miner sum at end of life. Above 1.0 the detail does not meet its fatigue requirement.' },
        { key: 'tEnd', label: 'Thickness at end of required life', unit: 'mm', value: tMin * 1e3, explain: 'Remaining plate thickness after corrosion once the protection is depleted.' },
        { key: 'shareWave', label: 'Share of damage from waves', unit: '%', value: r0.dW / r0.total * 100, explain: 'Part of the annual damage caused by wave-induced cycles; the remainder comes from operational cycles.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Damage accumulation', x: { label: 'Years from now', unit: 'yr' }, y: { label: 'Miner sum', unit: '-' },
          series: [{ name: 'Characteristic damage', x: yrs, y: Dc }, { name: `Allowable (1/DFF = ${fmt(1 / inp.dff)})`, x: [0, yrs[yrs.length - 1]], y: [1 / inp.dff, 1 / inp.dff], dash: true }, { name: 'Required life', x: [inp.life, inp.life], y: [0, Math.max(Dlife, 1 / inp.dff)], dash: true }],
          explain: 'Cumulative damage from the present state. The curve steepens once corrosion starts; it must stay below the allowable line up to the required life.' },
        { type: 'scatter', title: 'S–N curve and operational blocks', x: { label: 'Cycles to failure / cycles in the required life', unit: '-', log: true }, y: { label: 'Stress range (nominal)', unit: 'MPa', log: true },
          series: [{ name: `${inp.snName} (${inp.snEnv}), SCF and thickness included`, x: nn, y: sr }, { name: 'Blocks over the required life', x: blocks.map((b) => b.perYear * inp.life), y: blocks.map((b) => b.range) }, { name: 'Largest wave range', x: [1], y: [Math.max(inp.dS0, 1e-3)] }],
          explain: 'The design S–N curve against the applied blocks. Blocks close to or right of the curve consume life quickly; blocks below the knee stress contribute little.' },
        { type: 'bar', title: 'Annual damage by source', categories: src.map((s) => s.name), series: [{ name: 'Damage per year', y: src.map((s) => s.d) }], y: { label: 'Damage', unit: '1/yr' },
          explain: 'Which cycle source consumes the fatigue life. Mitigation should target the tallest bar.' },
        { type: 'line', title: 'Remaining thickness', x: { label: 'Years from now', unit: 'yr' }, y: { label: 'Thickness', unit: 'mm' }, series: [{ name: 'Thickness', x: yrs, y: tk }],
          explain: 'Plate thickness against time: constant while coating and anodes protect, then falling at the corrosion rate.' },
      ],
      tables: [{ title: 'Damage by source', columns: [{ key: 'name', label: 'Source' }, { key: 'd', label: 'Damage per year' }, { key: 'share', label: 'Share', unit: '%' }, { key: 'life', label: 'Life if alone', unit: 'yr' }],
        rows: src.map((s) => ({ name: s.name, d: round(s.d, 3), share: round(s.d / r0.total * 100), life: s.d > 0 ? round(1 / s.d) : '∞' })),
        explain: 'Annual characteristic damage of each source at the present thickness, its share and the life it alone would give.' }],
      notes: [`S–N curve ${inp.snName} (${inp.snEnv}): m = ${curve.m1}, knee at ${fmt(curve.kneeStress / 1e6)} MPa.`, `"${top.name}" governs the damage.`, inp.corrRate > 0 && inp.protLife < inp.life ? `Corrosion starts after ${inp.protLife} years and removes ${fmt((t0 - tMin) * 1e3)} mm by the end of the required life.` : 'No corrosion loss within the required life.'],
      advice,
      publish: { 'life.fatigueYears': remaining, 'subsea.fatigueUc': Math.min(uc, 99) },
      quantities: { vessel_days: { inspection: round(2 * Math.min(Math.max(Math.floor(inp.life / inspInt), 1), 10)) } },
    };
  },
  verify: [
    { name: 'Single block on curve D in air: N = 10^12.164 / Δσ³', source: 'DNV-RP-C203 Table 2-1 (D: m = 3, log a = 12.164): Δσ = 100 MPa → N = 1.4588e6; 1e5 cycles/yr → damage 0.06855 per year, life 14.59 yr',
      inputs: { snName: 'D', snEnv: 'air', scf: 1, thk: 25, dS0: 0, blocks: [{ name: 'b', range: 100, perYear: 1e5 }], corrRate: 0, dff: 1, D0: 0, age: 0, corr0: 0 },
      expect: [{ kpi: 'dYear', value: 1e5 / (10 ** 12.164 / 1e6), tol: 1e-3 }, { kpi: 'lifeChar', value: 14.588, tol: 1e-3 }, { kpi: 'remaining', value: 14.588, tol: 2e-3 }] },
    { name: 'Thickness correction: N scales with (t/t_ref)^(−k·m)', source: 'DNV-RP-C203 §2.4.3: curve D, k = 0.2, t = 50 mm → stress factor 2^0.2 = 1.1487, damage × 1.1487³ = 1.5157',
      inputs: { snName: 'D', snEnv: 'air', scf: 1, thk: 50, dS0: 0, blocks: [{ name: 'b', range: 100, perYear: 1e5 }], corrRate: 0, dff: 1 },
      expect: [{ kpi: 'dYear', value: 1e5 / (10 ** 12.164 / 1e6) * 2 ** 0.6, tol: 2e-3 }] },
  ],
};
