// Deck vibration and vibration-induced fatigue: natural frequencies and mode shapes of the topsides frame with
// its equipment mass, resonance margins against the machine orders, steady-state response to the unbalanced
// force by modal superposition and the fatigue life of the most stressed member end under the cyclic stress.
import { modalAnalysis } from '../../engine/structural/framemodal.js';
import { harmonicResponse } from '../../engine/structural/framedynamics.js';
import { memberForcesFromDisplacements } from '../../engine/structural/frame3d.js';
import { unbalanceForce, separationMargin, campbellData } from '../../engine/structural/machinery.js';
import { snCurve, snConstants, cyclesFromConstants, SN_CURVE_NAMES } from '../../engine/fatigue/sn.js';
import { fatigueLife } from '../../engine/fatigue/damage.js';
import { deckInputs, solverInput, buildDeck, equipmentOf, dynamicModel, stressOf, framePlot, need, run, fmt, round, adv, check, fsCheck, steelCo2Text, CONFIRM, DECK_LIMITS } from './_common.js';

const LIFE_CAP = 1e4;   // years reported when the stress range gives an effectively unlimited life

export default {
  id: 's02-fatigue-vibration', suite: 2, title: 'Deck vibration and vibration fatigue',
  summary: 'Natural frequencies and mode shapes of the deck with its equipment, resonance margins against machine orders, floor vibration amplitudes under the unbalanced force and the fatigue life of the cyclically stressed framing.',
  covers: [6], fidelity: 'numerical',
  standards: ['ISO 19901-3:2014 §7.6 (vibration)', 'ISO 20816-1:2016 (vibration severity, zone boundaries)', 'ISO 21940-11:2016 (balance quality grades)', 'DNV-RP-C203 (2021) §2.4, §5 (S–N curves, Miner sum)', 'NORSOK S-002 (2018) §5.5 (whole-body vibration at work areas)'],
  models: ['Eigenvalue analysis (consistent mass, subspace iteration)', 'Modal superposition', 'Forced harmonic response', 'Frequency response functions', 'Rotor unbalance forcing model', 'Campbell diagram', 'Mass–stiffness–damping formulation', 'Fatigue vibration model (S–N curve, Miner sum)'],
  limits: `${DECK_LIMITS} The response uses the lowest global modes of the beam model with one modal damping ratio: local plate and stringer modes above the retained modes, machine skid flexibility and isolator dynamics are not represented, so results for excitation above the highest retained mode are indicative only. Damping is an input, not identified from measurements.`,
  publishes: ['topsides.f1', 'life.fatigueYears'],
  searchTerms: 'offshore deck vibration machinery resonance natural frequency fatigue topsides',
  inputs: [
    ...deckInputs(),
    { key: 'machRow', label: 'Equipment row carrying the machine', type: 'int', default: 2, min: 1, max: 40, group: 'Machine', role: 'load', help: 'Row number in the equipment table; the dynamic force acts at the nearest frame node of its deck level.' },
    { key: 'speed', label: 'Operating speed', unit: 'rpm', type: 'number', default: 440, min: 30, max: 20000, group: 'Machine', role: 'bc', design: true, help: 'Running speed of the rotor or crankshaft (the 1× excitation frequency is speed / 60).' },
    { key: 'speedMin', label: 'Lowest continuous speed', unit: 'rpm', type: 'number', default: 400, min: 0, max: 20000, group: 'Machine', role: 'bc', help: 'Lower end of the operating speed range for the separation-margin check (start-up passes through lower speeds briefly).' },
    { key: 'rotorMass', label: 'Rotor mass', unit: 't', type: 'number', default: 6, min: 0.01, max: 500, group: 'Machine', role: 'load' },
    { key: 'grade', label: 'Balance quality grade G', unit: 'mm/s', type: 'number', default: 6.3, min: 0.4, max: 630, group: 'Machine', role: 'load', help: 'ISO 21940-11: G 2.5 turbines and compressors, G 6.3 general machinery, G 16–40 reciprocating engines.' },
    { key: 'service', label: 'Service factor on the unbalance', unit: '-', type: 'number', default: 2, min: 1, max: 10, group: 'Machine', role: 'load', help: 'Allowance for deterioration of the balance in service (project value).' },
    { key: 'Fextra', label: 'Additional harmonic force (vendor data)', unit: 'kN', type: 'number', default: 12, min: 0, max: 2000, group: 'Machine', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.3 },
      help: 'Primary unbalanced force of a reciprocating machine or other certified 1× force amplitude, added to the rotor unbalance.' },
    { key: 'zeta', label: 'Modal damping ratio', unit: '%', type: 'number', default: 2, min: 0.2, max: 10, group: 'Dynamic properties', role: 'ic', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Welded steel decks: 1–2 % bare, 2–3 % with outfitting. Calibrate against a modal test where available.' },
    { key: 'liveFrac', label: 'Fraction of the live load present as mass', unit: '-', type: 'number', default: 0.25, min: 0, max: 1, group: 'Dynamic properties', role: 'ic' },
    { key: 'nModes', label: 'Modes retained', type: 'int', default: 24, min: 3, max: 40, group: 'Dynamic properties', role: 'control' },
    { key: 'marginReq', label: 'Required separation margin', unit: '%', type: 'number', default: 20, min: 5, max: 50, group: 'Acceptance criteria', role: 'criteria', help: `±20 % between natural and excitation frequencies is common practice for machinery supports. ${CONFIRM}` },
    { key: 'vLimit', label: 'Vibration velocity limit (rms)', unit: 'mm/s', type: 'number', default: 4.5, min: 0.2, max: 50, group: 'Acceptance criteria', role: 'criteria', help: `Zone B/C boundary of the applicable part of ISO 20816 for the machine class, or the work-area limit. ${CONFIRM}` },
    { key: 'snName', label: 'S–N curve of the governing detail', type: 'select', default: 'F', group: 'Fatigue', role: 'criteria', options: SN_CURVE_NAMES.map((v) => ({ value: v, label: `DNV-RP-C203 ${v}` })), help: 'F for fillet-welded attachments and beam end connections; D for full-penetration butt welds.' },
    { key: 'scf', label: 'Stress concentration factor', unit: '-', type: 'number', default: 1.5, min: 1, max: 10, group: 'Fatigue', role: 'criteria', help: 'Geometric SCF at the connection not contained in the S–N class (from the local stress module).' },
    { key: 'hours', label: 'Running hours per year', unit: 'h', type: 'number', default: 8000, min: 1, max: 8766, group: 'Fatigue', role: 'operation' },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 2, min: 1, max: 10, group: 'Fatigue', role: 'criteria', help: `DNV / NORSOK N-004: 1–3 for accessible topsides details by consequence, up to 10 where not inspectable. ${CONFIRM}` },
    { key: 'life', label: 'Design life', unit: 'yr', type: 'number', default: 25, min: 1, max: 100, group: 'Fatigue', role: 'criteria', bind: 'life.designYears' },
    solverInput(1),
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'f1' },
  compute(inp, ctx) {
    // SI: rpm → Hz, t → kg, kN → N, % → fraction.
    const deck = buildDeck(inp, ctx), eqs = equipmentOf(deck, inp), zeta = inp.zeta / 100;
    need(inp.speedMin <= inp.speed, 'The lowest continuous speed must not exceed the operating speed.');
    need(eqs.length >= 1, 'Enter at least one equipment item: the machine sits on an equipment row.');
    const mq = eqs[Math.min(Math.round(inp.machRow), eqs.length) - 1];
    let node = -1, bd = Infinity;
    for (const n of deck.levelNodes[mq.k]) { const c = deck.model.nodes[n], d = Math.hypot(c[0] - mq.x, c[1] - mq.y); if (d < bd) { bd = d; node = n; } }
    const ref = dynamicModel(deck, inp, { liveFrac: inp.liveFrac, nSub: inp.nSub, eqs });
    const modal = run(() => modalAnalysis(ref.model, { nModes: Math.round(inp.nModes), tol: 1e-6 }), 'Modal analysis (check supports and bracing: a mechanism has no frequencies)');
    const fn = Array.from(modal.frequencies);
    need(fn.length && fn.every((f) => Number.isFinite(f) && f > 0), 'The modal analysis returned no positive natural frequency: the frame is a mechanism or has no mass.');
    const fe = inp.speed / 60;
    const ub = unbalanceForce({ mass: inp.rotorMass * 1e3, gradeMmS: inp.grade, speedRpm: inp.speed, serviceFactor: inp.service });
    const F0 = ub.force + inp.Fextra * 1e3;
    need(F0 > 0, 'The harmonic force is zero: enter a rotor mass and balance grade or a vendor force.');
        // Separation from the 1× and 2× lines over the operating speed range, for the modes that matter: those whose
    // resonant velocity at the machine support, F(f) φ²/(2ζω)/√2 with F ∝ f², would exceed the velocity limit.
    const margins = [];
    fn.forEach((f, k) => {
      const ph = modal.modes[k], wk = 2 * Math.PI * f, vRes = F0 * (f / fe) ** 2 * (ph[6 * node + 2] ** 2 + ph[6 * node + 1] ** 2) / (2 * zeta * wk) / Math.SQRT2;
      for (const ord of [1, 2]) { const lo = ord * inp.speedMin / 60, hi = ord * fe, m = f >= lo && f <= hi ? 0 : Math.min(Math.abs(lo - f), Math.abs(hi - f)) / f; margins.push({ mode: k + 1, f, ord, margin: m, vRes, live: vRes * 1e3 > inp.vLimit }); }
    });
    const liveModes = margins.filter((m) => m.live);
    const mgov = liveModes.length ? liveModes.reduce((a, c) => (c.margin < a.margin ? c : a)) : { mode: 0, f: 0, ord: 1, margin: 1 };
    const kNear = fn.reduce((b, f, k) => (Math.abs(f - fe) < Math.abs(fn[b] - fe) ? k : b), 0);
    const sep = separationMargin({ excitationHz: fe, naturalHz: fn[kNear], required: inp.marginReq / 100 });
    // Response at the operating speed: the unbalance rotates, so the vertical and horizontal components are 90° apart
    // and their responses are combined by square root of the sum of squares (re, im of both).
    const nN = ref.model.nodes.length, outs = [];
    for (let i = 0; i < nN; i++) for (let d = 0; d < 6; d++) outs.push({ node: i, dof: d });
    const fields = [];
    for (const dof of [2, 1]) {
      const h = run(() => harmonicResponse(modal, { frequencies: [fe], forces: [{ node, dof, amplitude: F0 }], outputs: outs, zeta }), 'Harmonic response');
      const re = new Float64Array(6 * nN), im = new Float64Array(6 * nN);
      for (let j = 0; j < outs.length; j++) { re[j] = h.re[j][0]; im[j] = h.im[j][0]; }
      fields.push(re, im);
    }
    const amp = new Array(nN).fill(0), U = new Float64Array(6 * nN);
    for (let i = 0; i < nN; i++) { let s = 0; for (const f of fields) s += f[6 * i] ** 2 + f[6 * i + 1] ** 2 + f[6 * i + 2] ** 2; amp[i] = Math.sqrt(s); }
    for (let j = 0; j < 6 * nN; j++) U[j] = Math.hypot(fields[0][j], fields[1][j]);
    let dMax = 0, nMax = 0; amp.forEach((v, i) => { if (v > dMax) { dMax = v; nMax = i; } });
    const w = 2 * Math.PI * fe, vRms = dMax * w / Math.SQRT2, acc = dMax * w * w;
    // Cyclic stress: member end forces of each response field, stress amplitudes combined by SRSS (upper bound).
    const owner = new Array(ref.model.members.length); ref.memberMap.forEach((ks, i) => { for (const k of ks) owner[k] = i; });
    const s2 = new Float64Array(ref.model.members.length);
    for (const f of fields) {
      const mf = memberForcesFromDisplacements(modal.system, f);
      mf.forEach((q, k) => { let s = 0; for (const o of [0, 6]) { const ss = stressOf(deck.orig[owner[k]].sec, { N: q[o], Vy: q[o + 1], Vz: q[o + 2], T: q[o + 3], My: q[o + 4], Mz: q[o + 5] }); s = Math.max(s, Math.abs(ss.sigmaMax), Math.abs(ss.sigmaMin)); } s2[k] += s * s; });
    }
    let sAmp = 0, kS = 0; s2.forEach((v, k) => { if (v > sAmp * sAmp) { sAmp = Math.sqrt(v); kS = k; } });
    const range = 2 * sAmp * inp.scf, og = deck.orig[owner[kS]];
    const sn = snConstants(run(() => snCurve(inp.snName, 'air'), 'S–N curve'), { thickness: og.sec.tf || og.sec.t || 0.025 });
    const Nf = cyclesFromConstants(sn, range), nYear = fe * 3600 * inp.hours, dRate = Number.isFinite(Nf) && Nf > 0 ? nYear / Nf : 0;
    const fl = dRate > 0 ? fatigueLife({ damageRate: dRate, designLife: inp.life, dff: inp.dff }) : null;
    const lifeYears = fl ? Math.min(fl.designLife, LIFE_CAP) : LIFE_CAP, ucFat = fl ? Math.min(fl.utilisation, 99) : 0;
    // Frequency response at the machine node (vertical force, vertical response) over the speed range.
    const fTop = Math.max(1.5 * fe, 1.2 * fn[Math.min(fn.length - 1, 5)]), fs = [];
    for (let i = 0; i < 240; i++) fs.push(fTop * (i + 1) / 240);
    for (const f of fn) if (f < fTop) fs.push(f);
    fs.push(fe); fs.sort((a, b) => a - b);
    const frf = run(() => harmonicResponse(modal, { frequencies: fs, forces: [{ node, dof: 2, amplitude: 1 }], outputs: [{ node, dof: 2 }], zeta }), 'Frequency response');
    // Velocity amplitude under the speed-dependent force F(f) = F0 (f/fe)² (unbalance grows with speed squared).
    const vSweep = fs.map((f, i) => frf.amplitude[0][i] * F0 * (f / fe) ** 2 * 2 * Math.PI * f / Math.SQRT2 * 1e3);
    const camp = run(() => campbellData({ modes: fn.slice(0, 8), orders: [1, 2], speedMinRpm: 0, speedMaxRpm: 1.3 * inp.speed, operatingRangeRpm: [inp.speedMin, inp.speed], required: inp.marginReq / 100, nPoints: 41 }), 'Campbell diagram');
    const mReq = inp.marginReq / 100;
    const checks = [
      fsCheck('Separation of natural frequencies from the 1× and 2× orders', Math.max(mgov.margin, 1e-4), mReq, 'ISO 19901-3 §7.6; project separation margin', `Smallest distance between the machine order lines over the operating speed range and any retained mode able to exceed the velocity limit at resonance (${liveModes.length ? `mode ${mgov.mode} at ${fmt(mgov.f)} Hz against ${mgov.ord}×` : 'no retained mode is able to'}). Inside the margin the deck runs at or near resonance and the amplitude is controlled only by damping.`),
      check('Deck vibration velocity (rms)', vRms * 1e3, inp.vLimit, 'mm/s', 'ISO 20816-1 zone B/C boundary (machine class value) / NORSOK S-002 work-area limit', 'Largest steady-state vibration velocity of the deck at the operating speed. Above the limit the machine trips on vibration, bearings and piping connections wear and the work area becomes uncomfortable.'),
      check('Vibration fatigue of the framing', ucFat, 1, 'Miner sum × DFF', `DNV-RP-C203 §2.4, curve ${inp.snName} in air, SCF ${inp.scf}, DFF ${inp.dff}`, 'Accumulated fatigue damage of the most stressed member end over the design life under the steady cyclic stress. Machinery cycles accumulate at millions per day, so the stress range must stay near the endurance range of the detail.'),
    ];
    const ok = checks.every((c) => c.uc <= 1), gov = checks.reduce((a, c) => (c.uc > a.uc ? c : a));
    const pm = ref.model.nodes[nMax];
    // Stiffening needed: f ∝ √(k/m) and k ∝ I ∝ depth² at constant area → depth factor for the governing mode to leave the band.
    const fTarget = fn[kNear] < fe ? fn[kNear] : fe * (1 + mReq), depthFac = fn[kNear] >= fe ? fTarget / fn[kNear] : 1;
    const beamT = deck.orig.filter((o) => o.horizontal).reduce((s, o) => s + o.sec.A * o.L * deck.steel.rho, 0) / 1e3;
    return {
      kpis: [
        { key: 'f1', label: 'Fundamental natural frequency', unit: 'Hz', value: fn[0], explain: 'Lowest natural frequency of the topsides with its equipment mass. Global sway modes are usually 0.5–3 Hz; vertical deck modes 3–15 Hz.' },
        { key: 'fNear', label: 'Natural frequency nearest to the running speed', unit: 'Hz', value: fn[kNear], explain: `The retained mode closest to the 1× excitation at ${fmt(fe)} Hz: the mode most likely to be excited in service.` },
        { key: 'margin', label: 'Smallest separation margin', unit: '%', value: mgov.margin * 100, limit: inp.marginReq, sense: 'min', explain: 'Smallest relative distance between the 1× or 2× order over the operating speed range and any retained mode that could exceed the velocity limit at resonance. Zero means such a mode lies inside the range; 100 % means no retained mode responds enough to matter.' },
        { key: 'F0', label: 'Harmonic force amplitude', unit: 'kN', value: F0 / 1e3, explain: 'Rotating unbalance m·G·Ω with the service factor plus the vendor force: the 1× force applied at the machine support.' },
        { key: 'disp', label: 'Largest displacement amplitude', unit: 'µm', value: dMax * 1e6, explain: 'Zero-to-peak steady-state displacement of the most responsive node at the operating speed.' },
        { key: 'vRms', label: 'Largest vibration velocity (rms)', unit: 'mm/s', value: vRms * 1e3, limit: inp.vLimit, sense: 'max', explain: 'Root-mean-square velocity of the most responsive node: the quantity that vibration severity standards limit.' },
        { key: 'acc', label: 'Largest acceleration amplitude', unit: 'm/s²', value: acc, explain: 'Peak acceleration at the same node; compare with equipment and comfort limits (about 0.1–0.5 m/s² for manned areas).' },
        { key: 'dSigma', label: 'Hot-spot stress range', unit: 'MPa', value: range / 1e6, explain: 'Twice the cyclic stress amplitude at the most stressed member end, times the stress concentration factor: the range entered in the S–N curve.' },
        { key: 'fatigueLife', label: 'Design fatigue life', unit: 'yr', value: lifeYears, limit: inp.life, sense: 'min', explain: `Life at the design fatigue factor under the steady machine vibration (reported as ${LIFE_CAP} years when the damage is negligible).` },
      ],
      checks,
      plots: [
        framePlot({ title: `Mode ${kNear + 1} (${fmt(fn[kNear])} Hz) — nearest to the running speed`, ref, u: modal.modes[kNear], values: Array.from({ length: nN }, (_, i) => Math.hypot(modal.modes[kNear][6 * i], modal.modes[kNear][6 * i + 1], modal.modes[kNear][6 * i + 2])), valueLabel: 'Mode amplitude (mass-normalised)', valueUnit: '1/√kg', explain: 'Shape of the mode closest to the 1× excitation. A machine placed at an antinode of this mode excites it most; moving it toward a column line reduces the response.' }),
        framePlot({ title: 'Vibration amplitude at the operating speed', ref, u: U, values: amp.map((v) => v * 1e6), valueLabel: 'Displacement amplitude', valueUnit: 'µm', explain: 'Steady-state displacement amplitude of every node under the unbalanced force. The response should decay away from the machine; a remote hot area indicates a mode being driven.' }),
        { type: 'line', title: 'Velocity response at the machine support against speed', x: { label: 'Excitation frequency', unit: 'Hz' }, y: { label: 'Vertical velocity (rms)', unit: 'mm/s' },
          series: [{ name: 'Response to the unbalance (force ∝ speed²)', x: fs.map((f) => round(f)), y: vSweep.map((v) => round(v)) }, { name: 'Limit', x: [fs[0], fs[fs.length - 1]], y: [inp.vLimit, inp.vLimit], dash: true }, { name: 'Operating speed', x: [fe, fe], y: [0, Math.max(inp.vLimit, ...vSweep)], dash: true }],
          explain: 'Vertical velocity of the machine support as the speed runs up. Peaks are resonances with deck modes; the operating line should sit in a valley, clear of the peaks by the separation margin.' },
        { type: 'line', title: 'Campbell diagram', x: { label: 'Speed', unit: 'rpm' }, y: { label: 'Frequency', unit: 'Hz' },
          series: [...camp.modeHz.map((m, k) => ({ name: `Mode ${k + 1}`, x: Array.from(camp.speedRpm, (v) => round(v)), y: Array.from(m, (v) => round(v)) })), ...camp.orderHz.map((o, k) => ({ name: `${k + 1}× order`, x: Array.from(camp.speedRpm, (v) => round(v)), y: Array.from(o, (v) => round(v)), dash: true }))],
          explain: 'Natural frequencies (horizontal lines) against the 1× and 2× excitation lines. Every crossing is a critical speed; none should fall between the lowest continuous and the operating speed.' },
      ],
      tables: [
        { title: 'Natural frequencies', columns: [{ key: 'mode', label: 'Mode' }, { key: 'f', label: 'Frequency', unit: 'Hz' }, { key: 'T', label: 'Period', unit: 's' }, { key: 'mx', label: 'Mass x', unit: '%' }, { key: 'my', label: 'Mass y', unit: '%' }, { key: 'mz', label: 'Mass z', unit: '%' }, { key: 'm1', label: 'Margin to 1×', unit: '%' }, { key: 'vRes', label: 'Velocity if resonant', unit: 'mm/s rms' }],
          rows: fn.map((f, k) => ({ mode: k + 1, f: round(f), T: round(1 / f), mx: round(modal.massFraction[k].x * 100, 3), my: round(modal.massFraction[k].y * 100, 3), mz: round(modal.massFraction[k].z * 100, 3), m1: round(Math.abs(fe - f) / f * 100, 3), vRes: round(margins[2 * k].vRes * 1e3, 3) })),
          explain: 'One row per retained mode with its effective mass fraction per direction (which way it moves), its distance from the running speed and the velocity it would reach at the machine support if the unbalance ran exactly at its frequency. Modes whose resonant velocity is below the limit need no separation.' },
        { title: 'Critical speeds', columns: [{ key: 'mode', label: 'Mode' }, { key: 'order', label: 'Order' }, { key: 'rpm', label: 'Speed', unit: 'rpm' }, { key: 'hz', label: 'Frequency', unit: 'Hz' }, { key: 'inRange', label: 'In operating range' }],
          rows: camp.crossings.slice(0, 40).map((c) => ({ mode: c.mode + 1, order: `${c.order}×`, rpm: round(c.speedRpm), hz: round(c.hz), inRange: c.inOperatingRange ? 'yes' : 'no' })),
          explain: 'Speeds at which an excitation order meets a natural frequency. Those marked as in the operating range are resonances in continuous service; the others are passed during start-up and shut-down.' },
      ],
      notes: [`${gov.name} governs at ${fmt(gov.uc)}. The 1× excitation at ${fmt(fe)} Hz is ${sep.tuning === 'low' ? 'above' : 'below'} mode ${kNear + 1} (${fmt(fn[kNear])} Hz): the support is ${sep.tuning}-tuned with ${fmt(sep.margin * 100)} % separation.`,
        `Largest response at node (${fmt(pm[0])}, ${fmt(pm[1])}, ${fmt(pm[2])}) m; most stressed member is a ${og.group.toLowerCase()} with a nominal stress amplitude of ${fmt(sAmp / 1e6)} MPa.`,
        `${fn.length} modes up to ${fmt(fn[fn.length - 1])} Hz retained${fe > fn[fn.length - 1] ? ' — the running speed lies above them, so the response is indicative: retain more modes or model the local deck' : ''}.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.85 ? 'warn' : 'ok', `Vibration velocity ${fmt(vRms * 1e3)} mm/s rms against ${inp.vLimit} mm/s; separation margin ${fmt(mgov.margin * 100)} % (required ${inp.marginReq} %); vibration fatigue life ${fmt(lifeYears)} years.`,
          'Keep the framing. Confirm the natural frequencies by a modal (hammer) test before commissioning and re-run with the measured damping; ask the vendor for certified unbalanced forces and moments.', 'ISO 19901-3:2014 §7.6; ISO 20816-1:2016; DNV-RP-C203 (2021) §2.4',
          `Stiffness, not strength, is what the girders and beams (${fmt(beamT)} t) provide here; avoiding resonance by placement costs no steel, whereas stiffening by mass would add ${steelCo2Text(0.1 * beamT)}`)
          : adv('engineering', 'fail', `${gov.name} is not met (utilisation ${fmt(gov.uc)}): velocity ${fmt(vRms * 1e3)} mm/s rms, margin ${fmt(mgov.margin * 100)} %, fatigue life ${fmt(lifeYears)} years.`,
            mgov.margin < mReq && liveModes.length ? `Move mode ${mgov.mode} (${fmt(mgov.f)} Hz) out of the ${mgov.ord}× band: ${depthFac > 1 ? `raise it to ${fmt(fTarget)} Hz with about ${fmt((depthFac - 1) * 100, 2)} % more girder and beam depth, or ` : ''}add a column or knee brace under the machine, relocate the machine onto a column line, or mount it on isolators tuned below ${fmt(fe / 2.5)} Hz.` : `Reduce the force (balance to a finer grade than G ${inp.grade}), add isolators, or stiffen the support beams under the machine; for fatigue, improve the detail class or reduce the SCF of ${inp.scf}.`,
            'ISO 19901-3:2014 §7.6; ISO 21940-11:2016; DNV-RP-C203 (2021) §2.4',
            `A column or brace under the machine is a few tonnes; raising the frequency by depth alone needs no extra steel at constant area, while ${fmt((depthFac ** 2 - 1) * 100, 2)} % more stiffness by area would cost about ${steelCo2Text(beamT * Math.max(depthFac ** 2 - 1, 0.05))}`),
        adv('economic', ok ? 'info' : 'warn', `A deck resonance found at commissioning costs downtime: each day of a ${fmt(inp.hours / 8766 * 100, 2)} % duty machine out of service is ${fmt(inp.hours / 365, 3)} running hours lost, and retrofit stiffening is hot work on a live deck.`,
          'Settle machine position, support stiffness and isolator selection at the design stage and budget one modal test per machine skid; it is far cheaper than vibration troubleshooting offshore.', 'Whole-life costing per ISO 15686-5; reliability-centred maintenance practice (ISO 14224 failure data for rotating equipment)'),
      ],
      publish: { 'topsides.f1': fn[0], 'life.fatigueYears': lifeYears },
    };
  },
  verify: [{
    // One bay, one level, four fixed-base legs CHS 300 × 10 of h = 20 m under a deep deck grid, 100 t at the plan
    // centre (shared equally by the four corner nodes). The stiff girders hold the leg heads against rotation, so
    // the sway stiffness is k = 4 · 12EI/h³ with I = π(0.3⁴ − 0.28⁴)/64 = 9.5896e-5 m⁴:
    //   k = 48 · 210e9 · 9.5896e-5 / 20³ = 120 829 N/m.
    // Mass: block 100 000 kg + girders and beams 7850·(60·0.0684 + 40·0.05184) = 48 492 kg + 13/35 of the legs
    // (fixed–guided sway shape) 0.3714 · 4 · 20 · π·0.29·0.01·7850 = 2 125 kg → M = 150 617 kg.
    //   f1 = √(k/M)/(2π) = √(120 829/150 617)/(2π) = 0.1425 Hz.
    name: 'Sway frequency of a stiff deck on four fixed-base legs: f = √(48EI/(h³ M))/(2π)', source: 'Chopra, Dynamics of Structures, §1.3 (shear frame); 13/35 column mass per Blevins Table 8-1 consistent mass',
    inputs: { live: 0, tPlate: 0, secSteel: 0, nx: 1, ny: 1, nLevels: 1, zLow: 20, cD: 300, ct: 10, bracing: 'none', columns: 'corners', nSub: 4, machRow: 1, equip: [{ name: 'Block', W: 100, x: 15, y: 10, level: 1, a: 0, b: 0, zc: 0 }] },
    expect: [{ kpi: 'f1', value: 0.1425, tol: 0.02 }],
  }],
};
