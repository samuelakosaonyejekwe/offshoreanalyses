// Fatigue and vibration of a piled wharf bent: natural frequencies with the entrained water, steady-state response
// to wave, crane and vortex-shedding excitation, and S–N fatigue of the pile-head connection under the yearly
// spectrum of berthing, wave, mooring and vessel-induced load cycles (with rainflow counting of a measured record).
import { modalAnalysis } from '../../engine/structural/framemodal.js';
import { harmonicResponse, dynamicAmplification } from '../../engine/structural/framedynamics.js';
import { snCurve, cyclesToFailure, stressRangeAtCycles, SN_CURVE_NAMES } from '../../engine/fatigue/sn.js';
import { fatigueLife } from '../../engine/fatigue/damage.js';
import { rainflowFourPoint } from '../../engine/fatigue/rainflow.js';
import { vortexSheddingFrequency } from '../../engine/hydro/morison.js';
import { G, need, run, r, f, cap, check, eng, eco, linspace, bentInputs, pileInputs, bentModel, solveBent, subdivideBent, framePlot, steelCo2, CONFIRM, NOT_HIFI, ECON_BASIS } from './_common.js';

const SPECTRUM_DEFAULT = [
  { name: 'Berthing, normal', range: 600, n: 150 }, { name: 'Berthing, hard (5 % of calls)', range: 900, n: 8 },
  { name: 'Mooring line load cycles', range: 150, n: 20000 }, { name: 'Wave loading, operating sea', range: 40, n: 2000000 },
  { name: 'Wave loading, storm', range: 120, n: 20000 }, { name: 'Passing ships and wash', range: 30, n: 5000 }, { name: 'Crane travel and braking', range: 70, n: 60000 },
];
const EXC_DEFAULT = [
  { name: 'Waves (peak period 7 s)', freq: 0.143, F: 60 }, { name: 'Crane trolley travel', freq: 0.35, F: 25 }, { name: 'Crane hoist / slewing', freq: 1.6, F: 15 }, { name: 'Vehicle traffic', freq: 2.5, F: 8 },
];

export default {
  id: 's05-fatigue-vibration', suite: 5, title: 'Fatigue and vibration',
  summary: 'Natural frequencies of the bent, resonance margins and steady-state vibration under wave, crane and traffic excitation, vortex shedding of the piles, and the fatigue life of the pile-head connection under the yearly load spectrum.',
  covers: [11], fidelity: 'numerical',
  standards: ['DNV-RP-C203 (2021) Fatigue design of offshore steel structures', 'EN 1993-1-9:2005 (fatigue)', 'DNV-RP-C205 (2019) §9 (vortex shedding)', 'ISO 10137:2007 (serviceability of structures against vibration)', 'ISO 19902:2020 §16 (fatigue)', 'ASTM E1049-85 (rainflow counting)'],
  models: ['Modal analysis (consistent mass, added water mass)', 'Steady-state harmonic response by modal superposition', 'Dynamic amplification factor of a damped oscillator', 'Vortex-shedding frequency (Strouhal)', 'S–N fatigue curves with thickness correction', 'Rainflow cycle counting (four-point rule)', 'Palmgren–Miner damage'],
  limits: `Linear elastic bent on equivalent-fixity piles; the stress at the pile-head weld is the nominal bending plus axial stress times the stress concentration factor entered (no local shell or weld-toe finite-element analysis). Load ranges are applied as horizontal force ranges on the bent at deck level. Fatigue of reinforcement in the concrete deck and crack-growth (fracture mechanics) assessment are not included. ${NOT_HIFI}`,
  publishes: ['wharf.f1', 'life.fatigueYears', 'wharf.fatigueDamageRate'],
  searchTerms: 'wharf pile fatigue S-N vibration natural frequency crane vortex shedding',
  inputs: [
    ...bentInputs(), ...pileInputs(),
    { key: 'gSdl', label: 'Superimposed dead load', unit: 'kPa', type: 'number', default: 2.5, min: 0, max: 30, group: 'Mass', role: 'ic' },
    { key: 'qMass', label: 'Cargo load present as mass', unit: 'kPa', type: 'number', default: 15, min: 0, max: 200, group: 'Mass', role: 'ic', help: 'Average cargo on the deck during operations; vibration is checked for the loaded deck (lower frequency) — run again with zero for the empty deck.' },
    { key: 'mCrane', label: 'Crane mass carried by this bent', unit: 't', type: 'number', default: 60, min: 0, max: 2000, group: 'Mass', role: 'ic' },
    { key: 'zeta', label: 'Damping ratio', unit: '%', type: 'number', default: 2, min: 0.3, max: 15, group: 'Vibration', role: 'criteria', help: 'About 1–2 % for welded steel piles with a concrete deck at small amplitude, 3–5 % with soil and hydrodynamic damping at larger amplitude.' },
    { key: 'exc', label: 'Harmonic excitation at deck level', type: 'table', group: 'Vibration', role: 'bc', columns: [{ key: 'name', label: 'Source', type: 'text' }, { key: 'freq', label: 'Frequency', unit: 'Hz', type: 'number' }, { key: 'F', label: 'Force amplitude on the bent', unit: 'kN', type: 'number' }],
      default: EXC_DEFAULT, help: 'Horizontal force amplitude and frequency of each periodic source acting on one bent: wave loading at the peak period, crane travel, hoisting and slewing harmonics from the supplier, traffic.' },
    { key: 'vLim', label: 'Deck vibration velocity limit', unit: 'mm/s', type: 'number', default: 15, min: 1, max: 100, group: 'Vibration', role: 'criteria', help: `Peak velocity for equipment and operator tolerance; crane suppliers and ISO 10137 give project values. ${CONFIRM}` },
    { key: 'Uc', label: 'Current speed at the piles', unit: 'm/s', type: 'number', default: 1.0, min: 0, max: 4, group: 'Vibration', role: 'environment', bind: 'site.Uc100' },
    { key: 'St', label: 'Strouhal number', unit: '-', type: 'number', default: 0.2, min: 0.15, max: 0.3, group: 'Vibration', role: 'criteria' },
    { key: 'spectrum', label: 'Yearly spectrum of horizontal load ranges on the bent', type: 'table', group: 'Fatigue loading', role: 'bc', columns: [{ key: 'name', label: 'Source', type: 'text' }, { key: 'range', label: 'Load range', unit: 'kN', type: 'number' }, { key: 'n', label: 'Cycles per year', unit: '1/yr', type: 'number' }],
      default: SPECTRUM_DEFAULT, help: 'Range (peak to trough) of the horizontal force on one bent and its yearly count, from the berthing, mooring, environmental and vessel-induced loading modules and the port call statistics.' },
    { key: 'history', label: 'Measured load history (optional, rainflow counted)', type: 'series', default: { x: [], y: [] }, xLabel: 'Time', xUnit: 's', yLabel: 'Horizontal load on the bent', yUnit: 'kN', group: 'Fatigue loading', role: 'bc', help: 'A record from strain gauges or load pins, converted to bent load. Its cycles are counted by rainflow and scaled to a year.' },
    { key: 'histDays', label: 'Duration represented by the measured history', unit: 'days', type: 'number', default: 1, min: 0.001, max: 3650, group: 'Fatigue loading', role: 'bc' },
    { key: 'sigPer', label: 'Nominal stress per kN of bent load (0 = from the frame analysis)', unit: 'MPa/kN', type: 'number', default: 0, min: 0, max: 10, group: 'Fatigue detail', role: 'bc', help: 'Stress at the pile-head weld per unit horizontal load on the bent. Enter a value from strain-gauge calibration to override the frame analysis.' },
    { key: 'scf', label: 'Stress concentration factor at the weld', unit: '-', type: 'number', default: 1.5, min: 1, max: 6, group: 'Fatigue detail', role: 'criteria', help: 'Pile-to-deck connection: 1.0–1.3 for a pile cast into the deck with shear keys, 1.5–2.5 for a welded cap plate or stiffened collar. Use detailed analysis for unusual details.' },
    { key: 'curve', label: 'S–N curve', type: 'select', default: 'D', group: 'Fatigue detail', role: 'criteria', options: SN_CURVE_NAMES.map((v) => ({ value: v, label: `DNV-RP-C203 curve ${v}` })), help: 'D for a full-penetration butt weld ground or with good profile on hot-spot stress; E–F for attachments on nominal stress.' },
    { key: 'env', label: 'Environment', type: 'select', default: 'seawater-cp', group: 'Fatigue detail', role: 'bc', options: [{ value: 'air', label: 'Air (above the splash zone, coated)' }, { value: 'seawater-cp', label: 'Seawater with cathodic protection' }, { value: 'free-corrosion', label: 'Free corrosion (splash zone, coating lost)' }] },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 3, min: 1, max: 10, group: 'Fatigue detail', role: 'criteria', help: `1–2 inspectable and non-critical, 3 not inspectable or splash zone, up to 10 for critical uninspectable details (DNV-ST-0126 / ISO 19902). ${CONFIRM}` },
    { key: 'life', label: 'Design life', unit: 'years', type: 'number', default: 50, min: 5, max: 120, group: 'Fatigue detail', role: 'criteria', bind: 'life.designYears' },
    { key: 'D0', label: 'Fatigue damage already accumulated', unit: '-', type: 'number', default: 0, min: 0, max: 1, group: 'Fatigue detail', role: 'ic', help: 'For an existing wharf: Miner sum consumed to date (without the design fatigue factor).' },
    { key: 'nModes', label: 'Modes', unit: '-', type: 'int', default: 8, min: 1, max: 24, group: 'Numerics', role: 'control' },
    { key: 'nSub', label: 'Elements per pile', unit: '-', type: 'int', default: 4, min: 1, max: 32, group: 'Numerics', role: 'control', help: 'Subdivision of each pile for the mass distribution; the local pile frequency needs at least four.' },
  ],
  mesh: { key: 'nSub', label: 'Elements per pile', levels: [1, 2, 4, 8], qoi: 'fPile' },
  compute(inp, ctx) {
    const B = bentModel({ ...inp, crane: [], Ffender: 1, Fmoor: 0, Fenv: 0, qLive: 0 }, ctx);
    const mDeck = ((inp.slabT * 25e3 + inp.gSdl * 1e3 + inp.qMass * 1e3) * B.sy * Math.max(B.deckLen, 0)) / G + inp.mCrane * 1e3;
    const heads = B.meta.headNodes, model = { ...subdivideBent(B.model, B.meta, inp.nSub).model, nodeMasses: heads.map((n) => ({ node: n, m: mDeck / heads.length })) };
    const nFree = (model.nodes.length - B.meta.toeNodes.length) * 2;
    const modal = run(() => modalAnalysis(model, { nModes: Math.max(1, Math.min(Math.round(inp.nModes), nFree)), waterline: 0, Ca: 1 }), 'Modal analysis');
    const f1 = modal.frequencies[0], zeta = inp.zeta / 100;
    // first mode with significant motion of the piles relative to the deck (local pile bending), for vortex shedding
    const swayIdx = 0, localIdx = modal.frequencies.findIndex((_, k) => k > 0 && modal.massFraction[k].x < 0.2);
    const fPile = localIdx > 0 ? modal.frequencies[localIdx] : modal.frequencies[Math.min(1, modal.frequencies.length - 1)];
    // harmonic response to each source and the frequency sweep for a unit force
    const exc = (inp.exc || []).map((e, i) => { need(e.freq > 0 && e.F >= 0, `Excitation ${i + 1}: frequency must be positive and the force not negative.`); return { name: String(e.name || `Source ${i + 1}`), freq: Number(e.freq), F: Number(e.F) * 1e3 }; });
    const out = [{ node: B.headSea, dof: 0 }], force = (F) => [{ node: B.headSea, dof: 0, amplitude: F }];
    const stat = harmonicResponse(modal, { frequencies: [1e-4], forces: force(1), outputs: out, zeta }).amplitude[0][0]; // static flexibility, m/N
    const resp = exc.map((e) => { const a = harmonicResponse(modal, { frequencies: [e.freq], forces: force(Math.max(e.F, 1e-9)), outputs: out, zeta }).amplitude[0][0]; return { ...e, amp: e.F > 0 ? a : 0, vel: e.F > 0 ? 2 * Math.PI * e.freq * a : 0, daf: dynamicAmplification({ ratio: e.freq / f1, zeta }).daf, r: e.freq / f1 }; });
    const worstV = resp.length ? resp.reduce((m, e) => (e.vel > m.vel ? e : m), resp[0]) : { vel: 0, name: '–', amp: 0, daf: 1, r: 0, freq: 0 };
    const near = resp.length ? resp.reduce((m, e) => (Math.abs(e.r - 1) < Math.abs(m.r - 1) ? e : m), resp[0]) : { r: 0, name: '–' };
    const sep = Math.abs(near.r - 1);
    // vortex shedding of a pile in the current: reduced velocity on the local pile frequency
    const fs = vortexSheddingFrequency({ U: inp.Uc, D: B.P.D, St: inp.St }), Vr = fPile > 0 ? inp.Uc / (fPile * B.P.D) : 0;
    // stress at the pile head per unit horizontal load on the bent
    const S = solveBent(B, inp), unit = S.combo({ B: 1 }); // case B holds 1 kN at the seaward pile head
    const perPile = unit.piles.map((p) => ({ id: p.id, s: p.Mhead / B.P.Z + Math.abs(p.N) / B.P.A, sMax: p.M / B.P.Z + Math.abs(p.N) / B.P.A }));
    const hot = perPile.reduce((m, p) => (p.s > m.s ? p : m), perPile[0]);
    const k = inp.sigPer > 0 ? inp.sigPer * 1e6 : hot.s; // Pa per kN, nominal
    const curve = run(() => snCurve(inp.curve, inp.env), 'S–N curve'), thickness = B.P.tNom;
    const blocks = [];
    for (const [i, q] of (inp.spectrum || []).entries()) { need(Number(q.range) >= 0 && Number(q.n) >= 0, `Spectrum row ${i + 1}: range and cycles must not be negative.`); if (q.range > 0 && q.n > 0) blocks.push({ src: String(q.name || `Block ${i + 1}`), load: Number(q.range), n: Number(q.n) }); }
    let nRain = 0;
    if (inp.history?.y?.length > 2) {
      const cyc = rainflowFourPoint(inp.history.y, { residue: 'half' }).cycles, perYear = 365.25 / inp.histDays, rmax = Math.max(...cyc.map((c) => c.range), 0);
      if (rmax > 0) { const bins = new Array(12).fill(0); for (const c of cyc) bins[Math.min(11, Math.floor(c.range / rmax * 12))] += c.count; bins.forEach((n, i) => { if (n > 0) { blocks.push({ src: 'Measured (rainflow)', load: rmax * (i + 1) / 12, n: n * perYear }); nRain += n; } }); }
    }
    need(blocks.length > 0, 'The fatigue loading holds no cycles: give at least one spectrum row or a measured history.');
    let rate = 0;
    for (const b of blocks) { b.range = k * b.load * inp.scf; b.N = cyclesToFailure({ curve, stressRange: k * b.load, thickness, scf: inp.scf }); b.d = Number.isFinite(b.N) && b.N > 0 ? b.n / b.N : 0; rate += b.d; }
    const safe = Math.max(rate, 1e-12), fl = fatigueLife({ damageRate: safe, designLife: inp.life, dff: inp.dff, accumulatedDamage: inp.D0 });
    const util = fl.utilisation ?? inp.dff * (inp.D0 + safe * inp.life), lifeChar = cap(fl.characteristicLife, 1e5), remaining = cap(Math.max(fl.remainingLife, 0), 1e5);
    const gov = blocks.reduce((m, b) => (b.d > m.d ? b : m), blocks[0]), maxRange = Math.max(...blocks.map((b) => b.range));
    const nTot = blocks.reduce((s, b) => s + b.n, 0);
    // constant-amplitude range that the detail tolerates for the same yearly count over the life (allowable equivalent)
    const sAllow = stressRangeAtCycles({ curve, cycles: Math.max(nTot * inp.life * inp.dff, 1), thickness, scf: 1 });
    const sFac = util > 1 ? util ** (-1 / 3) : 1;
    const checks = [
      check('Fatigue utilisation of the pile-head connection', util, 1, '-', 'DNV-RP-C203 (2021) §2.2; ISO 19902:2020 §16', 'Design fatigue factor times the Miner sum at the end of the design life. Above 1 a crack is expected to start at the weld toe within the life: it lets water in, grows under every berthing and ends in loss of the pile-head moment connection.'),
      check(`Vibration velocity of the deck (${worstV.name})`, worstV.vel * 1e3, inp.vLim, 'mm/s', 'ISO 10137:2007 Annex C; crane supplier limits', 'Peak velocity of the steady deck vibration under the strongest periodic source. Excess vibration disturbs crane positioning and operators and adds fatigue cycles not in the spectrum.'),
      check(`Separation from resonance (${near.name})`, 0.2, Math.max(sep, 1e-3), '-', 'DNV-RP-C205 (2019) §9; common practice of ±20 % separation', 'Distance of the nearest excitation frequency from the sway frequency against a 20 % band. Inside the band the response is amplified by up to 1/(2ζ).'),
      check('Vortex-shedding reduced velocity of the piles', Vr, 4.0, '-', 'DNV-RP-C205 (2019) §9.1.4 (cross-flow lock-in from about Vr = 3–5 in water)', 'Current speed over pile frequency times diameter. Above about 4 cross-flow lock-in can start: sustained pile vibration with a very large number of stress cycles.'),
    ];
    const w = checks.reduce((m, c) => (c.uc > m.uc ? c : m)), ok = w.uc <= 1;
    const fr = linspace(0.02, Math.max(3 * f1, 3), 240), sweep = harmonicResponse(modal, { frequencies: fr, forces: force(1), outputs: out, zeta }).amplitude[0];
    const mode = modal.modes[swayIdx]; let mx = 0; for (let i = 0; i < model.nodes.length; i++) mx = Math.max(mx, Math.hypot(mode[6 * i], mode[6 * i + 2]));
    const steelBent = B.pileLen * B.P.massPerM / 1e3;
    return {
      kpis: [
        { key: 'f1', label: 'Sway natural frequency', value: r(f1), unit: 'Hz', explain: `Lowest natural frequency of the bent with the mass entered and the water moving with the piles (period ${f(1 / f1)} s). Long-period waves and crane motions near it are amplified.` },
        { key: 'fPile', label: 'Local pile frequency', value: r(fPile, 6), unit: 'Hz', explain: 'First mode in which the piles bend between deck and fixity with little deck motion: the frequency that vortex shedding can excite.' },
        { key: 'kStat', label: 'Static sway stiffness of the bent', value: r(1 / stat / 1e6), unit: 'MN/m', explain: 'Horizontal force per unit deck movement from the modal model. It links load ranges to movement and should agree with the static frame analysis.' },
        { key: 'amp', label: 'Largest vibration amplitude', value: r(worstV.amp * 1e3), unit: 'mm', explain: `Steady-state deck displacement under ${worstV.name} (dynamic amplification ${f(worstV.daf)}).` },
        { key: 'vel', label: 'Largest vibration velocity', value: r(worstV.vel * 1e3), unit: 'mm/s', limit: inp.vLim, sense: 'max', explain: 'Peak velocity 2πf × amplitude of the strongest source, against the limit for equipment and people.' },
        { key: 'dafMax', label: 'Largest dynamic amplification factor', value: r(resp.length ? Math.max(...resp.map((e) => e.daf)) : 1), unit: '-', explain: 'Dynamic over static response of the sway mode for the excitation closest to resonance. Values well above 1 call for stiffening, detuning or damping.' },
        { key: 'fShed', label: 'Vortex-shedding frequency', value: r(fs), unit: 'Hz', explain: `Shedding frequency St·U/D of a pile in the current; the reduced velocity on the pile frequency is ${f(Vr)}.` },
        { key: 'sigPerKN', label: 'Hot-spot stress per kN of bent load', value: r(k * inp.scf / 1e6), unit: 'MPa/kN', explain: `Stress at the pile-head weld of pile ${hot.id} per unit horizontal load, including the stress concentration factor.` },
        { key: 'maxRange', label: 'Largest hot-spot stress range', value: r(maxRange / 1e6), unit: 'MPa', explain: 'Stress range of the largest load block. Above about twice the yield stress the detail is in low-cycle fatigue and the S–N method no longer applies.' },
        { key: 'rate', label: 'Fatigue damage per year', value: r(rate, 3), unit: '1/yr', explain: `Miner sum accumulated in one year of operation; ${gov.src} contributes ${f(gov.d / safe * 100)} %.` },
        { key: 'life', label: 'Characteristic fatigue life', value: r(lifeChar), unit: 'years', explain: 'Years to a Miner sum of 1 on the design S–N curve (mean minus two standard deviations), without the design fatigue factor.' },
        { key: 'util', label: 'Fatigue utilisation at the end of the design life', value: r(util), unit: '-', limit: 1, sense: 'max', explain: 'Design fatigue factor × (accumulated + future damage). It must not exceed 1.' },
        { key: 'remaining', label: 'Remaining design fatigue life', value: r(remaining), unit: 'years', explain: 'Years from now until the factored damage reaches 1, allowing for the damage already accumulated.' },
        { key: 'sAllow', label: 'Allowable constant stress range for this cycle count', value: r(sAllow / 1e6), unit: 'MPa', explain: 'Hot-spot stress range that the detail could carry for the same total number of cycles over the life with the design fatigue factor: a yardstick for the spectrum.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Frequency response of the deck to a unit force', x: { label: 'Frequency', unit: 'Hz' }, y: { label: 'Dynamic over static displacement', unit: '-' },
          series: [{ name: 'Amplification', x: Array.from(fr, (v) => r(v)), y: Array.from(sweep, (v) => r(v / stat)) }, ...resp.map((e) => ({ name: e.name, x: [r(e.freq), r(e.freq)], y: [0, r(e.daf)], dash: true }))],
          explain: 'Dynamic amplification of deck sway against forcing frequency, with each excitation source marked. A source standing on a peak is in resonance; the peak height is set by damping.' },
        framePlot(model, { disp: mode, scale: mx > 0 ? r(2 / mx, 3) : 1, title: `Sway mode (${f(f1)} Hz)`, explain: 'Shape of the lowest mode. All piles bend in double curvature with the deck translating: the stress peaks are at the pile heads and at fixity.' }),
        { type: 'bar', title: 'Fatigue damage by load source', categories: blocks.map((b) => b.src), series: [{ name: 'Share of yearly damage (%)', y: blocks.map((b) => r(b.d / safe * 100, 3)) }], x: { label: 'Source', unit: '' }, y: { label: 'Share of damage', unit: '%' },
          explain: 'Which loads consume the fatigue life. A few large berthing cycles and millions of small wave cycles can contribute equally; the tallest bar is where mitigation pays.' },
      ],
      tables: [
        { title: 'Fatigue load blocks', columns: [{ key: 'src', label: 'Source' }, { key: 'load', label: 'Load range', unit: 'kN' }, { key: 'range', label: 'Hot-spot range', unit: 'MPa' }, { key: 'n', label: 'Cycles per year' }, { key: 'N', label: 'Cycles to failure' }, { key: 'd', label: 'Damage per year' }],
          rows: blocks.map((b) => ({ src: b.src, load: r(b.load), range: r(b.range / 1e6), n: r(b.n), N: Number.isFinite(b.N) ? r(b.N, 3) : 'no damage', d: r(b.d, 3) })), explain: 'Each block of the spectrum with its stress range, the endurance on the S–N curve and the damage it adds each year (Palmgren–Miner).' },
        { title: 'Vibration response by source', columns: [{ key: 'name', label: 'Source' }, { key: 'freq', label: 'Frequency', unit: 'Hz' }, { key: 'r', label: 'f / f1' }, { key: 'daf', label: 'Amplification' }, { key: 'amp', label: 'Amplitude', unit: 'mm' }, { key: 'vel', label: 'Velocity', unit: 'mm/s' }],
          rows: resp.map((e) => ({ name: e.name, freq: r(e.freq), r: r(e.r, 3), daf: r(e.daf, 3), amp: r(e.amp * 1e3, 3), vel: r(e.vel * 1e3, 3) })), explain: 'Steady-state response of the deck to each periodic source: frequency ratio to the sway mode, amplification, displacement and velocity amplitude.' },
      ],
      notes: [`${curve.name} curve, ${inp.env}, wall ${f(thickness * 1e3)} mm, stress concentration ${f(inp.scf)}; ${gov.src} governs the damage.${nRain ? ` ${f(nRain)} cycles were counted from the measured history.` : ''}`,
        `Vibrating mass ${f(modal.totalMass.x / 1e3)} t; modes used: ${modal.frequencies.slice(0, 4).map((v) => f(v)).join(', ')} Hz. Nearest excitation at ${f(near.r)} of the sway frequency.`],
      advice: [
        ok ? eng('ok', `Fatigue utilisation ${f(util)} after ${f(inp.life)} years (characteristic life ${f(lifeChar)} years); sway frequency ${f(f1)} Hz, largest vibration ${f(worstV.vel * 1e3)} mm/s from ${worstV.name}.`,
          `Inspect the pile-head welds of pile ${hot.id} first: it carries the highest stress per unit load. Keep the cathodic protection working, because in free corrosion the life shortens several-fold; confirm the stress concentration factor for the detail as built.`, 'DNV-RP-C203 (2021); ISO 19902:2020 §16; ISO 10137:2007',
          `The fatigue reserve means no extra steel is needed: the present ${f(steelBent)} t of pile steel per bent serves the full life if protection is maintained, avoiding replacement (${steelCo2(steelBent)} per bent).`)
          : eng('fail', `${w.name}: ${f(w.uc)}. Fatigue utilisation ${f(util)}, characteristic life ${f(lifeChar)} years; largest vibration ${f(worstV.vel * 1e3)} mm/s.`,
            w.name.startsWith('Fatigue') ? `Reduce the hot-spot stress range by ${f((1 - sFac) * 100)} %: grind or profile the weld (a better S–N class), lower the stress concentration with a thicker pile-head section, or cut the governing cycles (${gov.src}).` : w.name.startsWith('Vortex') ? 'Raise the pile frequency (shorter free length, bracing) or fit strakes or fairings over the upper third of the wetted length.' : `Detune the bent from ${near.name}: change stiffness (rakers) or mass, or add damping; a 20 % separation is the usual target.`,
            'DNV-RP-C203 (2021) §2, §7 (improvement); DNV-RP-C205 (2019) §9',
            `Weld improvement adds no steel. Thickening the pile heads by 20 % over the top 3 m instead adds about ${f(0.2 * B.P.massPerM * 3 * B.nPileMembers / 1e3)} t per bent: ${steelCo2(0.2 * B.P.massPerM * 3 * B.nPileMembers / 1e3)}.`),
        eco(util > 0.7 ? 'warn' : 'info', `The connection uses ${f(util * 100)} % of its design fatigue budget over ${f(inp.life)} years; ${gov.src} causes ${f(gov.d / safe * 100)} % of the damage.`,
          'Plan weld inspection by risk: set the first inspection at about a third of the design fatigue life and the interval from the findings. An inspectable, repairable detail can justify a lower design fatigue factor than an inaccessible one.', `Risk-based inspection per DNV-RP-G101 and DNV-RP-C210; ${ECON_BASIS}`),
      ],
      publish: { 'wharf.f1': f1, 'life.fatigueYears': Math.min(lifeChar / inp.dff, 1e5), 'wharf.fatigueDamageRate': rate },
    };
  },
  verify: [
    // DNV-RP-C203 curve D in air: log N = 12.164 − 3 log Δσ. Δσ = 100 MPa (1 MPa/kN × 100 kN, SCF 1, wall 25 mm = reference thickness):
    // N = 10^12.164/100³ = 1.4588e6; n = 14 588 cycles/yr → damage 0.0100 per year → characteristic life 100 years.
    { name: 'Miner sum for one block on the D curve in air', source: 'DNV-RP-C203 (2021) Table 2-1, Eq. 2.4.1 (hand calculation)',
      inputs: { t: 25, corr: 0, sigPer: 1, scf: 1, curve: 'D', env: 'air', dff: 1, D0: 0, spectrum: [{ name: 'Block', range: 100, n: 14588 }] }, expect: [{ kpi: 'rate', value: 0.01, tol: 0.002 }, { kpi: 'life', value: 100, tol: 0.002 }] },
    // One pile as a cantilever with a lumped deck mass (see the seismic module): k = 4.873 MN/m, m = 500 t → f = (1/2π)√(k/m) = 0.4968 Hz
    // (pile and entrained water lower it by under 1 %). Forcing at 0.1 Hz: r = 0.2013, DAF = 1/(1 − r²) = 1.042; amplitude = 10 kN/k·DAF = 2.14 mm.
    { name: 'One-degree-of-freedom bent: f = (1/2π)√(k/m) and off-resonance amplitude', source: 'Chopra, Dynamics of Structures §2–3 (hand calculation)',
      inputs: { nPiles: 1, overhang: 0, D: 1, t: 25, corr: 0, deckLevel: 0.5, depth: 5, fixDepth: 5, mCrane: 500, qMass: 0, gSdl: 0, zeta: 2, exc: [{ name: 'Test', freq: 0.1, F: 10 }], Uc: 0 },
      expect: [{ kpi: 'f1', value: 0.4968, tol: 0.015 }, { kpi: 'amp', value: 2.14, tol: 0.025 }, { kpi: 'kStat', value: 4.873, tol: 0.02 }] },
  ],
};
