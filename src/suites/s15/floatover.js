// Float-over deck mating: two-mass spring–damper model of barge and deck with compression-only deck support
// units and leg mating units, ballast-driven load transfer and wave-induced barge heave.
import { rk4 } from '../../engine/math/ode.js';
import { makeRng } from '../../engine/math/rng.js';
import { waveSpectrum } from '../../engine/hydro/floater.js';
import { G, need, round, fmt, kpi, check, adv, line, tbl, thin, worst, linspace, vesselCo2Text, seaInputs, CONFIRM, LIMIT_ROM } from './_common.js';

export default {
  id: 's15-floatover', suite: 15, title: 'Float-over and topside mating',
  summary: 'Simulates the transfer of a topside from the float-over barge to the jacket legs as the barge is ballasted down in waves, and reports mating loads, impact velocity, the number of re-contacts, clearance margins and the ballast and time needed.',
  covers: [14], fidelity: 'reduced-order',
  standards: ['DNV-ST-N001 (2021) Sec. 14 (float-over: clearances, load transfer, ballast capacity and redundancy)', 'ISO 19901-6 (2009) §17 (float-over installation)', 'DNV-RP-C205 (2019) §3.5 (JONSWAP spectrum)'],
  models: ['Two-degree-of-freedom vertical dynamics of barge and deck (Newton)', 'Leg mating unit and deck support unit stiffness and damping (compression-only contact)', 'Impact response at first contact', 'Ballast-induced displacement variation at a prescribed rate', 'Hydrostatic restoring with added mass and linear damping', 'Froude–Krylov heave excitation of a box barge in a JONSWAP sea (seeded phases)', 'Mechanical load transfer between barge and jacket'],
  limits: `${LIMIT_ROM} Only heave is simulated: the barge is a box with constant added mass and damping, roll and pitch do not redistribute load between legs (all legs share equally), the jacket is rigid, and horizontal mating is a clearance comparison rather than a fender and surge simulation. One seeded wave realisation is run — repeat with other seeds for statistics. A coupled multi-body time-domain analysis with the real LMU load–deflection curves is required for design.`,
  publishes: ['install.floatoverLmuLoad', 'install.floatoverBallast'],
  searchTerms: 'float-over topside mating leg mating unit load transfer barge ballast simulation',
  inputs: [
    { key: 'Wd', label: 'Deck (topside) weight', unit: 't', type: 'number', default: 12000, min: 100, max: 60000, group: 'Deck', role: 'load', bind: 'install.liftWeight', uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'nLeg', label: 'Number of mating legs', unit: '-', type: 'int', default: 8, min: 4, max: 16, group: 'Deck', role: 'geometry' },
    { key: 'L', label: 'Barge length', unit: 'm', type: 'number', default: 180, min: 50, max: 400, group: 'Float-over barge', role: 'geometry' },
    { key: 'B', label: 'Barge beam (between jacket legs)', unit: 'm', type: 'number', default: 42, min: 15, max: 100, group: 'Float-over barge', role: 'geometry' },
    { key: 'T', label: 'Barge draft at the start of mating', unit: 'm', type: 'number', default: 7, min: 2, max: 25, group: 'Float-over barge', role: 'ic', help: 'Initial condition: draft with the deck on board, docked in the jacket slot, before ballasting down.' },
    { key: 'Ca', label: 'Heave added mass / displaced mass', unit: '-', type: 'number', default: 1.0, min: 0.2, max: 4, group: 'Float-over barge', role: 'load', uncertain: { dist: 'normal', cov: 0.2 }, help: 'About 1 for a wide shallow box; from the vessel hydrodynamic data when available.' },
    { key: 'zeta', label: 'Heave damping ratio', unit: '-', type: 'number', default: 0.15, min: 0.02, max: 1, group: 'Float-over barge', role: 'load', uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'gap0', label: 'Initial vertical gap, stabbing cone to LMU', unit: 'm', type: 'number', default: 0.6, min: 0.05, max: 5, group: 'Mating units', role: 'ic', design: true, help: 'Mating clearance when ballasting starts. It must exceed the wave-induced relative motion or the cones strike early.' },
    { key: 'kLmu', label: 'LMU vertical stiffness per leg', unit: 'MN/m', type: 'number', default: 60, min: 1, max: 5000, group: 'Mating units', role: 'bc', design: true, uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Elastomer stack stiffness over its working stroke (secant value).' },
    { key: 'zLmu', label: 'LMU damping ratio', unit: '-', type: 'number', default: 0.1, min: 0, max: 1, group: 'Mating units', role: 'bc' },
    { key: 'capLmu', label: 'LMU capacity per leg', unit: 't', type: 'number', default: 2600, min: 10, max: 20000, group: 'Mating units', role: 'bc', design: true },
    { key: 'kDsu', label: 'Total DSU (deck support unit) stiffness', unit: 'MN/m', type: 'number', default: 1200, min: 10, max: 50000, group: 'Mating units', role: 'bc', help: 'All deck support units on the barge together.' },
    { key: 'capDsu', label: 'Total DSU capacity', unit: 't', type: 'number', default: 20000, min: 100, max: 2e5, group: 'Mating units', role: 'bc' },
    { key: 'rate', label: 'Ballast rate as draft increase', unit: 'm/min', type: 'number', default: 0.15, min: 0.005, max: 2, group: 'Ballast', role: 'bc', design: true, help: 'Prescribed boundary condition: how fast the free-floating barge would sink. Rapid-ballast systems reach 0.1–0.3 m/min.' },
    { key: 't0', label: 'Time before ballasting starts', unit: 's', type: 'number', default: 120, min: 0, max: 1200, group: 'Ballast', role: 'operation', help: 'Barge held in position; the relative motion in this period shows the clearance margin.' },
    ...seaInputs({ Hs: 1.0, Tp: 8, hsMax: 5 }),
    { key: 'seed', label: 'Wave phase seed', unit: '-', type: 'int', default: 3, min: 1, max: 99999, group: 'Sea state', role: 'control' },
    { key: 'xH', label: 'Horizontal relative motion per metre of Hs (significant amplitude)', unit: 'm/m', type: 'number', default: 0.15, min: 0, max: 2, group: 'Sea state', role: 'bc', help: 'Surge/sway of the stabbing cone relative to the leg with mooring and fenders engaged, from the mating analysis.' },
    { key: 'rCone', label: 'Capture radius of the LMU receptor', unit: 'm', type: 'number', default: 0.5, min: 0.05, max: 3, group: 'Mating units', role: 'geometry' },
    { key: 'clearReq', label: 'Required vertical clearance before mating', unit: 'm', type: 'number', default: 0.25, min: 0, max: 3, group: 'Criteria', role: 'criteria', help: `Margin between the initial gap and the largest relative heave while holding. ${CONFIRM}` },
    { key: 'exitGap', label: 'Barge-to-deck gap required to withdraw', unit: 'm', type: 'number', default: 0.5, min: 0.05, max: 5, group: 'Criteria', role: 'criteria' },
    { key: 'tEnd', label: 'Simulated time', unit: 's', type: 'number', default: 1500, min: 200, max: 14400, group: 'Numerics', role: 'control' },
    { key: 'stepRate', label: 'Time steps per second', unit: '1/s', type: 'int', default: 10, min: 10, max: 400, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'stepRate', label: 'Time steps per second', levels: [10, 20, 40, 80], qoi: 'tTransfer' },
  compute(inp) {
    // t → kg, MN/m → N/m, m/min → m/s.
    const rho = 1025, md = inp.Wd * 1e3, Wd = md * G;
    const Aw = inp.L * inp.B, Kh = rho * G * Aw, Mb = rho * Aw * inp.T - md;
    need(Mb > 0, 'The deck is heavier than the barge displacement at this draft: increase the draft or the barge size.');
    const Meff = Mb * 1 + inp.Ca * rho * Aw * inp.T, Cb = 2 * inp.zeta * Math.sqrt(Kh * (Meff + md));
    const kL = inp.kLmu * 1e6 * inp.nLeg, cL = 2 * inp.zLmu * Math.sqrt(kL * md), kD = inp.kDsu * 1e6, cD = 2 * 0.1 * Math.sqrt(kD * md);
    const r = inp.rate / 60, g0 = inp.gap0, t0 = inp.t0;
    // Effective wave elevation for heave excitation of a box in head seas: Smith factor × length-averaging factor.
    const nW = 60, om = linspace(0.25, 2.0, nW), S = inp.Hs > 0 ? waveSpectrum(om, { Hs: inp.Hs, Tp: inp.Tp, gamma: 3.3 }) : new Float64Array(nW);
    const rng = makeRng(inp.seed), dw = om[1] - om[0];
    const amp = om.map((o, i) => { const k = o * o / G, x = k * inp.L / 2; return Math.sqrt(2 * S[i] * dw) * Math.exp(-k * inp.T) * Math.abs(Math.sin(x) / x); });
    const ph = om.map(() => 2 * Math.PI * rng.next());
    const nSteps = Math.round(inp.tEnd * inp.stepRate), h = inp.tEnd / nSteps;
    const eta = new Float64Array(2 * nSteps + 1);
    for (let k = 0; k <= 2 * nSteps; k++) { const t = 0.5 * k * h; let s = 0; for (let i = 0; i < nW; i++) s += amp[i] * Math.cos(om[i] * t + ph[i]); eta[k] = s; }
    const fDsu = (zb, zd, vb, vd) => Math.max(0, Wd + kD * (zb - zd) + cD * (vb - vd));
    const fLmu = (zd, vd) => (zd < -g0 ? Math.max(0, kL * (-zd - g0) - cL * vd) : 0);
    const f = (t, y) => {
      const zb = y[0], zd = y[1], vb = y[2], vd = y[3];
      const Fd = fDsu(zb, zd, vb, vd), Fl = fLmu(zd, vd);
      const zeq = -r * Math.max(t - t0, 0), e = eta[Math.min(Math.round(2 * t / h), 2 * nSteps)];
      return [vb, vd, (-Kh * (zb - zeq) - Cb * vb + Kh * e - (Fd - Wd)) / Meff, (Fd + Fl - Wd) / md];
    };
    const sol = rk4(f, 0, [0, 0, 0, 0], inp.tEnd, nSteps);
    const n = sol.t.length, tt = sol.t;
    const FL = new Float64Array(n), FD = new Float64Array(n), ZD = new Float64Array(n), ZB = new Float64Array(n), GAP = new Float64Array(n);
    let tContact = NaN, tTransfer = NaN, FLmax = 0, FDmax = 0, vImp = 0, nImpact = 0, relMax = 0, tExit = NaN, prevL = 0, lastD = 0;
    for (let k = 0; k < n; k++) {
      const y = sol.y[k];
      need(Number.isFinite(y[0]) && Number.isFinite(y[1]), 'The mating simulation diverged. Increase the number of time steps per second.');
      const Fd = fDsu(y[0], y[1], y[2], y[3]), Fl = fLmu(y[1], y[3]);
      FL[k] = Fl; FD[k] = Fd; ZD[k] = y[1]; ZB[k] = y[0]; GAP[k] = (y[1] - y[0]) - Wd / kD;
      if (Fl > 0 && prevL === 0) { nImpact++; if (Number.isNaN(tContact)) { tContact = tt[k]; vImp = Math.abs(y[3]); } }
      if (Fd > 0) lastD = tt[k];
      if (tt[k] <= t0) relMax = Math.max(relMax, Math.abs(y[1]));
      FLmax = Math.max(FLmax, Fl); FDmax = Math.max(FDmax, Fd); prevL = Fl;
    }
    const transferred = lastD < inp.tEnd - 2 * h && Number.isFinite(tContact);
    if (transferred) tTransfer = lastD;
    if (transferred) { // exit: last time the barge-to-deck gap is below the withdrawal clearance
      let lastSmall = tTransfer; for (let k = 0; k < n; k++) if (tt[k] > tTransfer && GAP[k] < inp.exitGap) lastSmall = tt[k];
      if (lastSmall < inp.tEnd - 2 * h) tExit = lastSmall;
    }
    need(Number.isFinite(tContact), `The stabbing cones never reach the leg mating units within ${inp.tEnd} s: at ${inp.rate} m/min the ${g0} m gap needs ${fmt(t0 + g0 / r)} s. Extend the simulated time or increase the ballast rate.`);
    // Quasi-static references: contact after g0/r; full transfer after a further W·(1/kLMU + 1/kDSU + 1/Kh)/r of ballasting.
    const tTransRef = Wd * (1 / kL + 1 / kD + 1 / Kh) / r;
    const hMpm = 1.86 * inp.xH * inp.Hs; // most probable maximum ≈ 1.86 × significant amplitude for about 1000 cycles
    const ballast = r * (Number.isFinite(tExit) ? tExit - t0 : inp.tEnd - t0) * Aw;
    const checks = [
      check('LMU load per leg', FLmax / inp.nLeg / G / 1e3, inp.capLmu, 't', 'DNV-ST-N001 (2021) Sec. 14 (float-over: mating loads on LMUs and legs)', 'Largest load on a leg mating unit during transfer, including impact and wave-induced variation. Overload crushes the elastomer or overstresses the leg stabbing point.'),
      check('DSU load', FDmax / G / 1e3, inp.capDsu, 't', 'DNV-ST-N001 (2021) Sec. 14 (deck support units and barge deck loads)', 'Largest load on the deck support units: the deck weight amplified by barge heave before and during transfer.'),
      check('Vertical clearance while holding before mating', relMax + inp.clearReq, g0, 'm (motion + margin against gap)', 'DNV-ST-N001 (2021) Sec. 14 (mating clearances)', 'The wave-induced vertical motion of the cones must stay inside the initial gap with margin, otherwise they strike the LMUs before the controlled ballast-down starts.'),
      check('Horizontal motion against LMU capture radius', hMpm, inp.rCone, 'm', 'DNV-ST-N001 (2021) Sec. 14 (horizontal mating tolerances)', 'The stabbing cone must stay within the receptor as it enters. Larger horizontal motion means the cone lands on the rim: side loads and possible jamming.'),
    ];
    if (!transferred) checks.push(check('Load transfer completed within the simulated time', 1.5, 1, '- (not completed)', 'DNV-ST-N001 (2021) Sec. 14 (ballast capacity for timely load transfer)', 'The deck is still partly on the barge at the end: the hammering stage is not passed. Extend the time or raise the ballast rate.'));
    const gov = worst(checks), ok = gov.uc <= 1;
    const [tp, a1, a2, a3, a4, a5] = thin([tt, FL, FD, ZD, ZB, GAP], 500);
    const dafL = FLmax / Wd, minutes = (Number.isFinite(tExit) ? tExit : inp.tEnd) / 60;
    return {
      kpis: [
        kpi('tContact', 'Time of first contact on the LMUs', tContact, 's', 'When the stabbing cones first touch. In calm water it equals the start time plus gap / ballast rate; waves bring it forward.'),
        kpi('tTransfer', 'Time when the deck is fully on the jacket', transferred ? tTransfer : inp.tEnd, 's', 'Last instant with load on the deck support units. Equal to the simulated time when the transfer is not complete.'),
        kpi('dtTransfer', 'Duration of load transfer (first contact to full transfer)', (transferred ? tTransfer : inp.tEnd) - tContact, 's', 'The critical stage with load shared between barge and jacket. Quasi-static value: W·(1/kLMU + 1/kDSU + 1/Kh)/rate.'),
        kpi('dtRef', 'Quasi-static transfer duration', tTransRef, 's', 'Closed-form duration without waves, for comparison with the simulated value.'),
        kpi('FlmuLeg', 'Maximum LMU load per leg', FLmax / inp.nLeg / G / 1e3, 't', 'Peak mating load on one leg, all legs sharing equally.', inp.capLmu, 'max'),
        kpi('dafL', 'Maximum LMU load / deck weight', dafL, '-', 'Dynamic factor on the final static leg load. Values well above 1 mean wave hammering after the weight is across.'),
        kpi('FdsuMax', 'Maximum DSU load', FDmax / G / 1e3, 't', 'Peak load on the deck support units, deck weight plus barge heave inertia.', inp.capDsu, 'max'),
        kpi('vImp', 'Impact velocity at first contact', vImp, 'm/s', 'Relative vertical velocity of the deck and the jacket when the cones touch; it sets the impact energy the LMU elastomers must absorb.'),
        kpi('nImpact', 'Number of separate contacts', nImpact, '-', 'How many times the cones touch, separate and touch again before staying in contact: a measure of hammering.'),
        kpi('relMax', 'Largest relative heave while holding', relMax, 'm', 'Vertical motion of the deck relative to the jacket before ballasting starts.', g0 - inp.clearReq, 'max'),
        kpi('hMpm', 'Horizontal relative motion, most probable maximum', hMpm, 'm', '1.86 × the significant amplitude entered per metre of Hs.', inp.rCone, 'max'),
        kpi('tExit', 'Time when the barge is clear to withdraw', Number.isFinite(tExit) ? tExit : inp.tEnd, 's', 'After this the gap between the deck underside and the barge supports stays above the required withdrawal clearance.'),
        kpi('ballast', 'Ballast taken on for mating and separation', ballast, 'm³', 'Ballast rate × waterplane area × duration from the start of ballasting to barge exit (or to the end of the simulation).'),
        kpi('Theave', 'Barge heave natural period with the deck on', 2 * Math.PI * Math.sqrt((Meff + md) / Kh), 's', 'If the wave peak period is close to it, barge heave and the mating loads are amplified.'),
      ],
      checks,
      plots: [
        line('Load transfer history', 'Time', 's', 'Load', 'MN', [{ name: 'LMUs (all legs)', x: tp, y: a1.map((v) => round(v / 1e6)) }, { name: 'DSUs', x: tp, y: a2.map((v) => round(v / 1e6)) }, { name: 'Deck weight', x: [0, inp.tEnd], y: [Wd / 1e6, Wd / 1e6], dash: true }],
          'The deck weight moves from the barge supports to the legs as the barge is ballasted down. Spikes on the LMU curve are wave impacts; the stage between first contact and zero DSU load is the one to keep short.'),
        line('Vertical motions and barge-to-deck gap', 'Time', 's', 'Displacement', 'm', [{ name: 'Deck', x: tp, y: a3.map((v) => round(v, 4)) }, { name: 'Barge', x: tp, y: a4.map((v) => round(v, 4)) }, { name: 'Barge-to-deck gap', x: tp, y: a5.map((v) => round(Math.max(v, 0), 4)) }, { name: 'LMU contact level', x: [0, inp.tEnd], y: [-g0, -g0], dash: true }],
          'The deck follows the barge down until it lands on the LMUs, then stays while the barge keeps sinking and the gap opens. The gap curve must rise above the withdrawal clearance and stay there.'),
      ],
      tables: [tbl('Mating stages', [['k', 'Stage'], ['t', 'Time', 's'], ['v', 'Value'], ['u', 'Unit']], [
        { k: 'Ballasting starts', t: round(t0), v: round(relMax, 3), u: 'm largest relative heave while holding' },
        { k: 'First contact on LMUs', t: round(tContact), v: round(vImp, 3), u: 'm/s impact velocity' },
        { k: 'Deck fully on jacket', t: round(transferred ? tTransfer : inp.tEnd), v: round(FLmax / 1e6), u: 'MN peak LMU load (all legs)' },
        { k: 'Barge clear to withdraw', t: round(Number.isFinite(tExit) ? tExit : inp.tEnd), v: round(ballast), u: 'm³ ballast taken on' },
      ], 'The sequence of the mating with the time and governing value of each stage, for the installation procedure and the ballast plan.')],
      notes: [`First contact at ${fmt(tContact)} s (calm-water value ${fmt(t0 + g0 / r)} s), ${transferred ? `full transfer at ${fmt(tTransfer)} s` : 'transfer not completed'}; ${nImpact} separate contact(s).`, `Peak LMU load ${fmt(dafL)} × deck weight; barge heave natural period ${fmt(2 * Math.PI * Math.sqrt((Meff + md) / Kh))} s.`, `"${gov.name}" governs at ${fmt(gov.uc)}.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 || nImpact > 3 ? 'warn' : 'ok', `Mating completes in ${fmt(minutes)} min: peak LMU load ${fmt(FLmax / inp.nLeg / G / 1e3)} t per leg (${fmt(dafL)} × static), impact velocity ${fmt(vImp)} m/s, ${nImpact} contact(s).`,
          `Run several wave seeds and the Hs/Tp range of the installation criteria; keep the transfer stage (${fmt((transferred ? tTransfer : inp.tEnd) - tContact)} s) as short as the ballast system allows and verify LMU stroke for ${fmt(FLmax / kL)} m of compression.`, 'DNV-ST-N001 (2021) Sec. 14; ISO 19901-6 §17')
          : adv('engineering', 'fail', `"${gov.name}" is not met (${fmt(gov.uc)}): LMU ${fmt(FLmax / inp.nLeg / G / 1e3)} t per leg, relative heave ${fmt(relMax)} m in a ${g0} m gap, horizontal motion ${fmt(hMpm)} m.`,
            gov.name.startsWith('LMU') ? 'Soften the LMUs (longer elastomer stroke), increase the ballast rate to shorten the hammering stage, or lower the limiting sea state.' : gov.name.startsWith('Vertical') ? `Increase the initial gap to at least ${fmt(relMax + inp.clearReq)} m or reduce the limiting Hs to about ${fmt(inp.Hs * (g0 - inp.clearReq) / Math.max(relMax, 1e-6))} m.` : gov.name.startsWith('Horizontal') ? `Enlarge the receptor to ${fmt(hMpm)} m capture radius or stiffen the surge/sway fendering and mooring.` : gov.name.startsWith('Load transfer') ? 'Extend the simulated time or raise the ballast rate.' : 'Add deck support units or reduce barge heave with a lower limiting sea state.',
            'DNV-ST-N001 (2021) Sec. 14; ISO 19901-6 §17'),
        adv('economic', 'info', `Float-over installs the ${fmt(inp.Wd)} t deck in one piece without a heavy-lift vessel; the mating itself takes about ${fmt(minutes)} min and ${fmt(ballast)} m³ of ballast at ${inp.rate} m/min.`,
          'The limiting sea state of the mating (typically Hs 1.0–1.5 m with swell limits) decides the waiting time: use the weather-window module with this limit, and compare the rapid-ballast system and LMU cost against offshore hook-up of a lifted, modularised deck.',
          'Integrated-deck float-over versus modular lift: installed-cost comparison incl. offshore hook-up hours (AACE RP 18R-97); weather risk per DNV-ST-N001 Sec. 2',
          `An integrated deck completed onshore avoids months of offshore hook-up; each flotel/construction-vessel day avoided is about ${vesselCo2Text(1, 30)}`),
      ],
      publish: { 'install.floatoverLmuLoad': FLmax / inp.nLeg, 'install.floatoverBallast': ballast },
      quantities: { vessel_days: { floatoverBarge: round(minutes / 1440 + 1, 3) } },
    };
  },
  verify: [
    { name: 'Calm-water mating: contact after gap/rate, transfer after W·(1/kLMU + 1/kDSU + 1/Kh)/rate', source: 'Quasi-static springs in series (statics): gap 0.6 m at 0.15 m/min → 240 s after the 120 s hold = 360 s; W = 117.68 MN, kLMU = 480 MN/m, kDSU = 1200 MN/m, Kh = 1025·9.80665·180·42 = 75.99 MN/m → 1.892 m of ballasting → 756.7 s',
      inputs: { Hs: 0 }, expect: [{ kpi: 'tContact', value: 360, tol: 0.01 }, { kpi: 'dtTransfer', value: 12e6 * 9.80665 * (1 / 480e6 + 1 / 1200e6 + 1 / (1025 * 9.80665 * 180 * 42)) / 0.0025, tol: 0.02 }, { kpi: 'dafL', value: 1, tol: 0.01 }] },
  ],
};
