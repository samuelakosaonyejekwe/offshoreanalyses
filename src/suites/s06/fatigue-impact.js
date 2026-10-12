// Fatigue of the dolphin piles under repeated berthing impacts and cyclic mooring / wave loads:
// S–N curves and Palmgren–Miner summation, with optional rainflow counting of a measured load history.
import { snCurve, cyclesToFailure, SN_CURVE_NAMES } from '../../engine/fatigue/sn.js';
import { fatigueLife } from '../../engine/fatigue/damage.js';
import { rainflowFourPoint } from '../../engine/fatigue/rainflow.js';
import { need, run, fmt, round, cap, adv, check, CONFIRM, NOT_HIFI, ECON_BASIS, vesselCo2Text } from './_common.js';

const BERTH_DEFAULT = [{ share: 50, R: 500 }, { share: 30, R: 900 }, { share: 15, R: 1300 }, { share: 5, R: 1575 }];
const CYCLIC_DEFAULT = [{ range: 80, n: 500000 }, { range: 200, n: 40000 }, { range: 450, n: 1500 }, { range: 800, n: 40 }];

export default {
  id: 's06-fatigue-impact', suite: 6, title: 'Fatigue under repeated berthing and cyclic loads',
  summary: 'Cumulative fatigue damage and remaining life of the critical pile weld from the berthing-impact spectrum and cyclic mooring or wave loads.',
  covers: [7], fidelity: 'analytical',
  standards: ['DNV-RP-C203 (2021) §2 (S–N curves, Miner sum)', 'ISO 19902:2020 §16', 'BS 6349-2:2019 §6 (fatigue of maritime structures)', 'ASTM E1049-85 (rainflow counting)'],
  models: ['S–N fatigue relationships (DNV-RP-C203)', 'Palmgren–Miner cumulative damage', 'Impact-cycle fatigue (one stress cycle per berthing)', 'Rainflow counting (four-point rule)', 'Thickness correction', 'Design fatigue factor and accumulated-damage start state'],
  limits: 'Nominal stress scaled linearly from load (stress per unit load taken at the design load level, so soil nonlinearity makes small loads slightly conservative). One full stress cycle 0 → peak → 0 per berthing; dynamic ringing after impact is not counted. Crack-growth (fracture mechanics) assessment of found defects is outside this module. ' + NOT_HIFI,
  publishes: ['life.fatigueYears', 'dolphin.fatigueDamageRate'],
  searchTerms: 'fatigue repeated berthing impact pile weld S-N Miner dolphin',
  inputs: [
    { key: 'sigPer', label: 'Nominal stress per unit lateral load', unit: 'MPa/MN', type: 'number', default: 150, min: 0.1, max: 5000, group: 'Stress transfer', role: 'bc', bind: 'dolphin.stressPerLoad', uncertain: { dist: 'normal', cov: 0.1 },
      help: 'Bending stress at the critical weld per MN of horizontal load on the dolphin, from the global strength or flexible-dolphin module.' },
    { key: 'scf', label: 'Stress concentration factor', unit: '-', type: 'number', default: 1.3, min: 1, max: 6, group: 'Stress transfer', role: 'geometry', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Girth-weld misalignment or pile-to-cap connection geometry (DNV-RP-C203 §3.3).' },
    { key: 'curve', label: 'S–N curve', type: 'select', default: 'D', group: 'Detail', role: 'material', options: SN_CURVE_NAMES.map((v) => ({ value: v, label: `DNV-RP-C203 curve ${v}` })), help: 'D for a two-sided girth weld, E or F for single-sided or attachment welds.' },
    { key: 'env', label: 'Environment', type: 'select', default: 'seawater-cp', group: 'Detail', role: 'environment', options: [{ value: 'air', label: 'Air (above splash zone, coated)' }, { value: 'seawater-cp', label: 'Sea water with cathodic protection' }, { value: 'free-corrosion', label: 'Free corrosion (splash zone, no protection)' }] },
    { key: 't', label: 'Wall thickness at the weld', unit: 'mm', type: 'number', default: 28, min: 6, max: 150, group: 'Detail', role: 'geometry', bind: 'dolphin.pileT', help: 'Above the reference thickness the S–N strength is reduced (size effect).' },
    { key: 'fy', label: 'Yield strength at the weld', unit: 'MPa', type: 'number', default: 355, min: 200, max: 1000, group: 'Detail', role: 'material' },
    { key: 'nBerth', label: 'Berthings per year on this dolphin', unit: '1/yr', type: 'number', default: 300, min: 0, max: 20000, group: 'Loading history', role: 'operation', uncertain: { dist: 'normal', cov: 0.15 } },
    { key: 'berthing', label: 'Berthing spectrum', type: 'table', group: 'Loading history', role: 'bc', columns: [{ key: 'share', label: 'Share of berthings', unit: '%', type: 'number' }, { key: 'R', label: 'Dolphin reaction', unit: 'kN', type: 'number' }],
      default: BERTH_DEFAULT, help: 'Each berthing is one stress cycle from zero to the reaction and back. Buckling fenders reach near-rated reaction even in light berthings.' },
    { key: 'cyclic', label: 'Cyclic mooring / wave load spectrum', type: 'table', group: 'Loading history', role: 'bc', columns: [{ key: 'range', label: 'Load range', unit: 'kN', type: 'number' }, { key: 'n', label: 'Cycles per year', unit: '1/yr', type: 'number' }],
      default: CYCLIC_DEFAULT, help: 'Storm and moored-vessel load ranges on the dolphin (from the mooring and environmental modules or measurements).' },
    { key: 'history', label: 'Measured load history (optional)', type: 'series', default: { x: [], y: [] }, xLabel: 'Time', xUnit: 's', yLabel: 'Dolphin load', yUnit: 'kN', group: 'Loading history', role: 'bc', help: 'Rainflow-counted and scaled to one year by the recorded duration below; added to the tabulated spectra.' },
    { key: 'histDays', label: 'Duration represented by the measured history', unit: 'days', type: 'number', default: 30, min: 0.01, max: 3650, group: 'Loading history', role: 'bc' },
    { key: 'D0', label: 'Fatigue damage already accumulated', unit: '-', type: 'number', default: 0, min: 0, max: 5, group: 'Service state', role: 'ic', help: 'Characteristic Miner sum consumed to date (0 for a new structure; age × damage rate for an existing one).' },
    { key: 'life', label: 'Required further service life', unit: 'years', type: 'number', default: 30, min: 1, max: 120, group: 'Service state', role: 'criteria', bind: 'life.designYears' },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 3, min: 1, max: 10, group: 'Service state', role: 'criteria', help: `1 inspectable and redundant, 2–3 submerged or not inspectable, up to 10 for critical non-inspectable details (DNV-ST-0126 / ISO 19902). ${CONFIRM}` },
  ],
  compute(inp) {
    const curve = run(() => snCurve(inp.curve, inp.env), 'S–N curve');
    const k = inp.sigPer * inp.scf; // hot-spot MPa per MN  ≡ Pa per N
    const thickness = inp.t / 1e3;
    const blocks = [];
    const shareSum = (inp.berthing || []).reduce((s, r) => s + (r.share > 0 ? r.share : 0), 0);
    need(shareSum > 0 || (inp.cyclic || []).length > 0 || inp.history?.y?.length > 2, 'Give at least one berthing class, cyclic load block or measured history.');
    for (const r of inp.berthing || []) if (r.share > 0 && r.R > 0) blocks.push({ src: 'Berthing', load: r.R, n: inp.nBerth * r.share / shareSum });
    for (const r of inp.cyclic || []) if (r.n > 0 && r.range > 0) blocks.push({ src: 'Cyclic', load: r.range, n: r.n });
    let nRain = 0;
    if (inp.history?.y?.length > 2) {
      const cyc = rainflowFourPoint(inp.history.y, { residue: 'half' }).cycles;
      const perYear = 365.25 / inp.histDays;
      // bin the counted cycles into 12 classes of load range to keep the table readable
      const rmax = Math.max(...cyc.map((c) => c.range), 0);
      if (rmax > 0) {
        const bins = new Array(12).fill(0);
        for (const c of cyc) bins[Math.min(11, Math.floor(c.range / rmax * 12))] += c.count;
        bins.forEach((n, i) => { if (n > 0) { blocks.push({ src: 'Measured (rainflow)', load: rmax * (i + 1) / 12, n: n * perYear }); nRain += n; } });
      }
    }
    need(blocks.length > 0, 'The loading history holds no cycles: check the tables.');
    let rate = 0;
    for (const b of blocks) {
      b.range = k * b.load * 1e3;                       // Pa (kN → N)
      b.N = cyclesToFailure({ curve, stressRange: b.range / inp.scf, thickness, scf: inp.scf });
      b.d = Number.isFinite(b.N) && b.N > 0 ? b.n / b.N : 0;
      rate += b.d;
    }
    const rateB = blocks.filter((b) => b.src === 'Berthing').reduce((s, b) => s + b.d, 0);
    const maxRange = Math.max(...blocks.map((b) => b.range));
    const safeRate = Math.max(rate, 1e-12);
    const fl = fatigueLife({ damageRate: safeRate, designLife: inp.life, dff: inp.dff, accumulatedDamage: inp.D0 });
    const util = fl.utilisation ?? inp.dff * (inp.D0 + safeRate * inp.life);
    const lifeChar = cap(fl.characteristicLife, 1e5), remaining = cap(Math.max(fl.remainingLife, 0), 1e5);
    // allowable berthings per year keeping the utilisation at 1 with the other sources unchanged
    const other = rate - rateB;
    const nAllow = rateB > 0 ? inp.nBerth * Math.max((1 / inp.dff - inp.D0) / inp.life - other, 0) / rateB : 1e6;
    // stress reduction needed (m ≈ 3 on the governing branch)
    const sFac = util > 1 ? util ** (-1 / 3) : 1;
    const ok = util <= 1;
    const checks = [
      check('Fatigue damage at the end of the required life (with design fatigue factor)', util, 1, '-', 'DNV-RP-C203 (2021) §2.2 Eq. (2.2.1); ISO 19902:2020 §16.12', 'Miner sum including damage to date, times the design fatigue factor. Above 1 a through-thickness crack at the weld is predicted within the required life.'),
      check('Largest hot-spot stress range against twice yield (low-cycle limit)', maxRange / 1e6, 2 * inp.fy, 'MPa', 'DNV-RP-C203 (2021) §2.4 validity of the high-cycle S–N approach', 'The S–N curves apply to nominally elastic cycling. Ranges approaching twice yield need a strain-based low-cycle assessment.'),
    ];
    const order = blocks.slice().sort((a, b) => b.d - a.d);
    const sLine = []; const nLine = [];
    for (let i = 0; i <= 24; i++) { const s = 10 ** (1 + i * 2 / 24); const N = cyclesToFailure({ curve, stressRange: s * 1e6 / inp.scf, thickness, scf: inp.scf }); if (Number.isFinite(N) && N > 0) { sLine.push(round(s)); nLine.push(round(Math.log10(N))); } }
    const yrs = Array.from({ length: 13 }, (_, i) => round(inp.life * i / 12));
    const advice = [
      adv('engineering', ok ? 'ok' : 'fail', `Damage rate ${rate.toExponential(2)} per year (${fmt(rate > 0 ? rateB / rate * 100 : 0)} % from berthing impacts); characteristic life ${fmt(lifeChar)} years, ${fmt(remaining)} years remaining with DFF ${fmt(inp.dff)}.`,
        ok ? `Utilisation ${fmt(util)} at ${fmt(inp.life)} years. Schedule the first weld inspection of the critical section at about ${fmt(Math.min(remaining, inp.life) / 3)} years and confirm the berthing reaction spectrum with fender-panel load records.`
          : `Reduce the hot-spot stress range by ${fmt((1 - sFac) * 100)} % (thicker can at the weld, ground weld toes to move up one S–N class, or a softer fender), or limit berthings to ${fmt(nAllow)} per year, or inspect at intervals that justify a lower design fatigue factor.`,
        'DNV-RP-C203 (2021) §2; ISO 19902:2020 §16; DNV-RP-C210 for inspection planning',
        `Extending the fatigue life by weld improvement avoids pile replacement: one underwater inspection or repair campaign is about ${vesselCo2Text(2)}`),
      adv('economic', ok ? 'info' : 'warn', `The detail supports about ${fmt(nAllow)} berthings per year at the required ${fmt(inp.life)}-year life; the present plan is ${fmt(inp.nBerth)}.`,
        ok ? 'Use the reserve to plan inspection by risk rather than by calendar: a fatigue utilisation below about 0.3 supports extending underwater weld inspection intervals.' : 'Compare the cost of weld-toe grinding at fabrication (small, one-off) with in-service underwater inspection every few years and the downtime of a berth out of service for repair.',
        `Risk-based inspection per DNV-RP-G101 / API RP 580; ${ECON_BASIS}`),
    ];
    return {
      kpis: [
        { key: 'damageRate', label: 'Fatigue damage per year', unit: '1/yr', value: rate, explain: 'Characteristic Miner sum accumulated in one year of operation (without the design fatigue factor). Its inverse is the characteristic fatigue life.' },
        { key: 'life', label: 'Characteristic fatigue life', unit: 'years', value: lifeChar, explain: 'Years to a Miner sum of 1.0 from new, without the design fatigue factor. It should exceed the service life times the design fatigue factor.' },
        { key: 'remaining', label: 'Remaining fatigue life', unit: 'years', value: remaining, limit: inp.life, sense: 'min', explain: 'Years from now until the allowable damage is used up, including the design fatigue factor and the damage already accumulated.' },
        { key: 'util', label: 'Fatigue utilisation at the required life', unit: '-', value: util, limit: 1, sense: 'max', explain: 'Design fatigue factor times the damage at the end of the required life. Below 1 the detail is acceptable; below about 0.3 it needs no in-service inspection for fatigue.' },
        { key: 'maxRange', label: 'Largest hot-spot stress range', unit: 'MPa', value: maxRange / 1e6, explain: 'Stress range of the heaviest block including the stress concentration factor. Fatigue damage grows with its cube.' },
        { key: 'berthShare', label: 'Share of damage from berthing impacts', unit: '%', value: rate > 0 ? rateB / rate * 100 : 0, explain: 'Part of the yearly damage caused by berthing; the rest comes from cyclic mooring, wave or measured loads. It shows which load to reduce.' },
        { key: 'nAllow', label: 'Allowable berthings per year', unit: '1/yr', value: Math.min(nAllow, 1e6), explain: 'Berthing frequency at which the utilisation reaches 1.0 at the required life with the other load sources unchanged: an operating limit for the berth.' },
        { key: 'cycles', label: 'Stress cycles per year', unit: '1/yr', value: blocks.reduce((s, b) => s + b.n, 0), explain: 'Total number of counted cycles in a year from all sources.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Fatigue damage by load block', categories: order.slice(0, 12).map((b) => `${b.src} ${fmt(b.load)} kN`), series: [{ name: 'Damage per year', y: order.slice(0, 12).map((b) => Number(b.d.toPrecision(4))) }], y: { label: 'Damage per year', unit: '1/yr' },
          explain: 'Contribution of each load class to the yearly damage, largest first. Typically a few heavy classes dominate: that is where a stress reduction pays.' },
        { type: 'line', title: 'S–N curve and load blocks', x: { label: 'log10 cycles', unit: '-' }, y: { label: 'Hot-spot stress range', unit: 'MPa' },
          series: [{ name: `Curve ${inp.curve} (${inp.env})`, x: nLine, y: sLine }, { name: 'Applied cycles in the required life', x: blocks.map((b) => round(Math.log10(Math.max(b.n * inp.life, 1)))), y: blocks.map((b) => round(b.range / 1e6)) }],
          explain: 'Each point is a load block: its stress range against the number of cycles in the required life. Points close to or right of the curve consume most of the life.' },
        { type: 'line', title: 'Accumulated damage over the required life', x: { label: 'Years from now', unit: 'yr' }, y: { label: 'Miner sum × design fatigue factor', unit: '-' },
          series: [{ name: 'Design damage', x: yrs, y: yrs.map((y) => round(inp.dff * (inp.D0 + rate * y))) }, { name: 'Limit', x: [0, inp.life], y: [1, 1], dash: true }],
          explain: 'Growth of factored fatigue damage with time from the present state. The crossing with the limit line is the end of the fatigue life.' },
      ],
      tables: [{ title: 'Damage calculation', columns: [{ key: 'src', label: 'Source' }, { key: 'load', label: 'Load (range)', unit: 'kN' }, { key: 'n', label: 'Cycles per year' }, { key: 'range', label: 'Hot-spot range', unit: 'MPa' }, { key: 'N', label: 'Cycles to failure' }, { key: 'd', label: 'Damage per year' }],
        rows: blocks.map((b) => ({ src: b.src, load: round(b.load), n: round(b.n), range: round(b.range / 1e6), N: Number.isFinite(b.N) ? Number(b.N.toPrecision(4)) : 'no damage', d: Number(b.d.toPrecision(4)) })),
        explain: 'One row per load block: stress range, cycles to failure on the S–N curve and the Miner contribution nᵢ/Nᵢ per year.' }],
      notes: [`S–N curve ${inp.curve}, ${inp.env}; thickness ${fmt(inp.t)} mm (reference ${fmt(curve.tRef * 1e3)} mm, exponent ${curve.k}).`, ...(nRain ? [`${fmt(nRain)} rainflow cycles counted in the measured history and scaled by ${fmt(365.25 / inp.histDays)} to one year.`] : []), ok ? 'Fatigue life meets the requirement.' : 'Fatigue life does not meet the requirement.'],
      advice,
      publish: { 'life.fatigueYears': cap(fl.designLife, 1e5), 'dolphin.fatigueDamageRate': rate },
      quantities: {},
    };
  },
  verify: [
    // Curve D in air: log N = 12.164 − 3·log Δσ; Δσ = 100 MPa → N = 10^6.164 = 1.4588e6 cycles; 1000 cycles/yr → life 1458.8 yr.
    { name: 'Single block on DNV curve D in air: N = 10^12.164 / Δσ³', source: 'DNV-RP-C203 Table 2-1, curve D (m = 3, log a = 12.164); Miner sum of one block',
      inputs: { sigPer: 100, scf: 1, curve: 'D', env: 'air', t: 25, nBerth: 1000, berthing: [{ share: 100, R: 1000 }], cyclic: [], D0: 0, dff: 1, life: 30 },
      expect: [{ kpi: 'life', value: 1458.8, tol: 0.002 }, { kpi: 'damageRate', value: 6.8549e-4, tol: 0.002 }, { kpi: 'maxRange', value: 100, tol: 1e-9 }] },
    // Two blocks: 100 MPa × 1000/yr and 200 MPa × 100/yr → N2 = 1.4588e6/8 = 1.8235e5; D = 6.8549e-4 + 5.4839e-4 = 1.23388e-3 per year.
    { name: 'Miner sum of two blocks', source: 'Palmgren–Miner rule with the DNV-RP-C203 curve D in air (hand calculation)',
      inputs: { sigPer: 100, scf: 1, curve: 'D', env: 'air', t: 25, nBerth: 1000, berthing: [{ share: 100, R: 1000 }], cyclic: [{ range: 2000, n: 100 }], D0: 0, dff: 1, life: 30 },
      expect: [{ kpi: 'damageRate', value: 1.23388e-3, tol: 0.002 }] },
  ],
};
