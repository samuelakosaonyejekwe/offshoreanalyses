// Vortex-induced motion of a deep-draft floater in current: transverse natural period on the mooring, reduced
// velocity and lock-in band, motion amplitude from a user amplitude ratio, drag amplification, the resulting offset
// and restoring demand, and the yearly exposure to lock-in from a Weibull current climate.
import { vimScreening } from '../../engine/hydro/impact.js';
import { LIMIT_CFD, need, run, round, rnd, fmt, adv, check, sev, worst, linspace, steelCo2Text } from './_common.js';

export default {
  id: 's11-vim',
  suite: 11,
  title: 'Vortex-induced motion',
  summary: 'Whether the design current locks vortex shedding onto the sway natural period of a spar or deep-draft column hull, how large the resulting motion and offset are, and how many hours and cycles per year the unit spends in lock-in.',
  covers: [14],
  fidelity: 'screening',
  standards: ['DNV-RP-C205 (2019) §9 (vortex-induced oscillations, reduced velocity, drag amplification)', 'DNV-OS-E301 (2021) Ch.2 Sec.2 [2.8] (VIM in mooring design)', 'API RP 2SK (3rd ed., 2005, R2015) §5.6 and Appendix H (VIM design guidance)', 'ISO 19901-7 (2013) §8.5'],
  models: ['Single-degree-of-freedom sway natural period on the linearised mooring stiffness', 'Reduced velocity Ur = U·Tn/D and Strouhal shedding period', 'Lock-in band with user amplitude ratio A/D (model-test data)', 'Drag amplification Cd = Cd0·(1 + k·A/D)', 'Weibull current-speed distribution for lock-in exposure'],
  limits: `Screening only: the lock-in band and the amplitude ratio A/D are user inputs to be taken from model tests of the actual hull with its appendages and strakes; no amplitude is predicted. One uniform current over the draft, one transverse mode, linear mooring stiffness; multi-column wake interaction, heading dependence and sheared currents are not modelled. ${LIMIT_CFD}`,
  publishes: ['vessel.vimAmplitude'],
  searchTerms: 'vortex induced motion spar semi-submersible reduced velocity strakes',
  inputs: [
    { key: 'D', label: 'Hull / column diameter', unit: 'm', type: 'number', default: 30, min: 3, max: 60, group: 'Hull geometry', role: 'geometry', design: true, help: 'Diameter of the spar hull or of one column of a deep-draft semi-submersible (the characteristic shedding dimension).' },
    { key: 'draft', label: 'Draft', unit: 'm', type: 'number', default: 90, min: 5, max: 250, group: 'Hull geometry', role: 'ic', help: 'Submerged length of the column exposed to current.' },
    { key: 'nCol', label: 'Number of columns', unit: '-', type: 'int', default: 1, min: 1, max: 8, group: 'Hull geometry', role: 'geometry', help: '1 for a spar; 4 for a deep-draft semi-submersible. Projected area and column added mass scale with it.' },
    { key: 'massT', label: 'Displacement', unit: 't', type: 'number', default: 65000, min: 500, max: 2e6, group: 'Mass properties', role: 'ic', bind: 'vessel.disp', uncertain: { dist: 'normal', cov: 0.03 }, help: 'Mass of the unit in this loading condition.' },
    { key: 'Ca', label: 'Sway added-mass coefficient of the columns', unit: '-', type: 'number', default: 1.0, min: 0.3, max: 2, group: 'Mass properties', role: 'bc', help: 'Added mass = Ca·ρ·(πD²/4)·draft per column; 1.0 for a long circular cylinder.' },
    { key: 'maExtraT', label: 'Additional sway added mass (pontoons, plates)', unit: 't', type: 'number', default: 0, min: 0, max: 2e6, group: 'Mass properties', role: 'bc', help: 'Added mass of pontoons, heave plates or a truss section not covered by the column term.' },
    { key: 'Kmoor', label: 'Horizontal mooring stiffness', unit: 'kN/m', type: 'number', default: 250, min: 5, max: 50000, group: 'Station keeping', role: 'bc', bind: 'mooring.Kxx', uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Linearised restoring stiffness at the mean offset. It sets the sway natural period and the offset.' },
    { key: 'Hcap', label: 'Horizontal restoring capacity of the mooring', unit: 'kN', type: 'number', default: 12000, min: 100, max: 500000, group: 'Station keeping', role: 'criteria', help: 'Horizontal force the intact mooring can resist at its allowable line tension (from the mooring suite).' },
    { key: 'Uc', label: 'Design current speed', unit: 'm/s', type: 'number', default: 1.2, min: 0, max: 4, group: 'Current', role: 'environment', bind: 'site.Uc100', uncertain: { dist: 'normal', cov: 0.15 }, help: 'Extreme (loop/eddy or 100-year) current, taken uniform over the draft.' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1040, group: 'Current', role: 'environment', bind: 'site.rhoW', help: 'Density of the water.' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 1500, min: 50, max: 4000, group: 'Current', role: 'environment', bind: 'site.depth', help: 'Used for the offset limit as a percentage of depth.' },
    { key: 'weibA', label: 'Weibull scale of the current speed', unit: 'm/s', type: 'number', default: 0.4, min: 0.02, max: 3, group: 'Current', role: 'environment', help: 'Scale parameter of the long-term distribution of current speed at the site.' },
    { key: 'weibK', label: 'Weibull shape of the current speed', unit: '-', type: 'number', default: 1.5, min: 0.5, max: 5, group: 'Current', role: 'environment', help: 'Shape parameter (1.2–2 is typical).' },
    { key: 'UrLow', label: 'Lock-in onset reduced velocity', unit: '-', type: 'number', default: 4, min: 2, max: 8, group: 'VIM response (model-test data)', role: 'bc', help: 'Lower bound of the transverse lock-in band; 4–5 for bare circular hulls. User input.' },
    { key: 'UrHigh', label: 'Lock-in end reduced velocity', unit: '-', type: 'number', default: 10, min: 5, max: 20, group: 'VIM response (model-test data)', role: 'bc', help: 'Upper bound of the lock-in band. User input.' },
    { key: 'AoverD', label: 'Transverse amplitude ratio A/D in lock-in (bare hull)', unit: '-', type: 'number', default: 0.5, min: 0, max: 1.5, group: 'VIM response (model-test data)', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Single amplitude over diameter from model tests: 0.3–0.9 for bare spars. No amplitude is predicted by the module.' },
    { key: 'strakeEff', label: 'Strake effectiveness (reduction of A/D)', unit: '-', type: 'number', default: 0.6, min: 0, max: 0.95, group: 'VIM response (model-test data)', role: 'bc', design: true, help: 'Fraction by which helical strakes reduce the amplitude ratio (0 = no strakes; 0.5–0.8 for well-designed strakes).' },
    { key: 'St', label: 'Strouhal number', unit: '-', type: 'number', default: 0.2, min: 0.1, max: 0.3, group: 'VIM response (model-test data)', role: 'bc', help: 'Shedding frequency = St·U/D; 0.2 in the subcritical range, higher beyond the drag crisis.' },
    { key: 'Cd0', label: 'Drag coefficient without VIM', unit: '-', type: 'number', default: 0.7, min: 0.3, max: 2, group: 'VIM response (model-test data)', role: 'load', help: 'In-line drag coefficient on the projected area D × draft (0.6–0.8 smooth at high Reynolds number, 1.0–1.4 with strakes).' },
    { key: 'kAmp', label: 'Drag amplification per unit A/D', unit: '-', type: 'number', default: 2.1, min: 0, max: 4, group: 'VIM response (model-test data)', role: 'load', help: 'Cd = Cd0·(1 + k·A/D); DNV-RP-C205 quotes 2.1 for vibrating cylinders.' },
    { key: 'offsetMaxPct', label: 'Allowable offset', unit: '% of depth', type: 'number', default: 5, min: 0.5, max: 20, group: 'Acceptance criteria', role: 'criteria', help: 'Riser or umbilical limit on the total excursion in the extreme current condition.' },
  ],
  compute(inp) {
    need(inp.UrHigh > inp.UrLow, 'The lock-in band is inverted: the end reduced velocity must exceed the onset value.');
    const rho = inp.rhoW; const D = inp.D; const K = inp.Kmoor * 1e3; const m = inp.massT * 1e3;
    const ma = inp.Ca * rho * Math.PI * D * D / 4 * inp.draft * inp.nCol + inp.maExtraT * 1e3;
    const Tn = 2 * Math.PI * Math.sqrt((m + ma) / K);
    const aRatio = inp.AoverD * (1 - inp.strakeEff);
    const band = [inp.UrLow, inp.UrHigh];
    const at = (U) => {
      const v = run(() => vimScreening({ U, Tn, D, lockInRange: band, AoverD: aRatio, St: inp.St }), 'VIM screening');
      const A = v.amplitude || 0; const Cd = inp.Cd0 * (1 + inp.kAmp * A / D);
      const F = 0.5 * rho * Cd * D * inp.draft * inp.nCol * U * U; const mean = F / K;
      return { ...v, A, Cd, F, mean, total: Math.hypot(mean, A) };
    };
    const d = at(inp.Uc);
    // the worst current is not always the largest: scan up to the design current
    const us = linspace(0, Math.max(inp.Uc, d.criticalSpeeds[1] * 1.15, 0.1), 80); const sweep = us.map(at);
    const env = sweep.filter((s, i) => us[i] <= inp.Uc + 1e-12).reduce((a, s) => (s.total > a.total ? s : a), d);
    const offAllow = inp.offsetMaxPct / 100 * inp.depth;
    // long-term exposure: P(U1 < U < U2) of the Weibull current climate
    const cdf = (u) => 1 - Math.exp(-((u / inp.weibA) ** inp.weibK));
    const pLock = cdf(d.criticalSpeeds[1]) - cdf(d.criticalSpeeds[0]);
    const hours = pLock * 8766; const cycles = hours * 3600 / Tn;
    const checks = [
      check('Total offset in the design current (mean drift + VIM amplitude)', env.total, offAllow, 'm', 'DNV-OS-E301 (2021) Ch.2 Sec.2 [2.8]; API RP 2SK (2005) §5.6 — offset limits with VIM included', 'Guards risers and umbilicals, whose top angle and touchdown zone tolerate only a limited excursion. VIM adds a slow transverse oscillation to the mean current offset.'),
      check('Mooring restoring demand at the VIM extreme', K * env.total / 1e3, inp.Hcap, 'kN', 'DNV-OS-E301 (2021) Ch.2 Sec.2 (line tension with VIM-amplified drag)', 'Guards the mooring lines: the amplified drag and the transverse excursion both add line tension, and the many VIM cycles also consume chain fatigue life.'),
    ];
    const gov = worst(checks);
    // three-start strakes, height 0.13·D, 16 mm plate, 10 % allowance for the helix length
    const strakeSteel = 3 * (0.13 * D) * (inp.draft * 1.1) * 0.016 * 7.85 * inp.nCol;
    return {
      kpis: [
        { key: 'Tn', label: 'Sway natural period on the mooring', value: round(Tn, 5), unit: 's', explain: 'Period of the transverse oscillation of the hull on its mooring, including added mass. VIM locks onto this period.' },
        { key: 'addedMass', label: 'Sway added mass', value: round(ma / 1e3, 5), unit: 't', explain: 'Mass of water accelerated with the hull in sway; for a deep cylinder it is about equal to the displacement.' },
        { key: 'Ur', label: 'Reduced velocity at the design current', value: round(d.Ur, 5), unit: '-', explain: `U·Tn/D. Lock-in is assumed between ${fmt(inp.UrLow)} and ${fmt(inp.UrHigh)} (user band); shedding and natural frequency coincide at 1/St = ${fmt(1 / inp.St)}.` },
        { key: 'lockIn', label: 'Lock-in at the design current (1 = yes)', value: d.lockIn ? 1 : 0, unit: '-', explain: 'Whether the design current lies inside the lock-in band. Lock-in is common for spars and is then designed for rather than avoided.' },
        { key: 'sheddingT', label: 'Vortex-shedding period at the design current', value: round(Number.isFinite(d.sheddingPeriod) ? d.sheddingPeriod : 0, 5), unit: 's', explain: 'D/(St·U): compare with the sway natural period (zero when there is no current).' },
        { key: 'Uonset', label: 'Current speed at lock-in onset', value: round(d.criticalSpeeds[0], 4), unit: 'm/s', explain: 'Lowest current speed at which VIM is expected; currents below it give no lock-in.' },
        { key: 'Uend', label: 'Current speed at the end of lock-in', value: round(d.criticalSpeeds[1], 4), unit: 'm/s', explain: 'Current speed above which shedding is too fast to lock onto the sway mode.' },
        { key: 'amplitude', label: 'VIM amplitude', value: round(d.A, 4), unit: 'm', explain: 'Transverse single amplitude = effective A/D × diameter when locked in (zero outside the band).' },
        { key: 'CdAmp', label: 'Amplified drag coefficient', value: round(d.Cd, 4), unit: '-', explain: 'In-line drag coefficient including the increase caused by the transverse motion.' },
        { key: 'Fdrag', label: 'Current drag force', value: round(d.F / 1e3, 4), unit: 'kN', explain: 'Mean in-line current force on the submerged columns with the amplified drag coefficient.' },
        { key: 'meanOffset', label: 'Mean current offset', value: round(d.mean, 4), unit: 'm', explain: 'Steady in-line offset = drag force / mooring stiffness.' },
        { key: 'totalOffset', label: 'Largest total offset up to the design current', value: round(env.total, 4), unit: 'm', limit: round(offAllow, 4), sense: 'max', explain: 'Vector sum of the mean in-line offset and the transverse VIM amplitude, maximised over all currents up to the design value.' },
        { key: 'lockHours', label: 'Hours per year in lock-in', value: round(hours, 4), unit: 'h/yr', explain: 'Time the long-term current distribution spends inside the lock-in band: the exposure for mooring and riser fatigue.' },
        { key: 'vimCycles', label: 'VIM cycles per year', value: round(cycles, 4), unit: '1/yr', explain: 'Number of slow transverse oscillations per year, the cycle count for the VIM fatigue contribution of the mooring lines.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Reduced velocity and VIM amplitude versus current', x: { label: 'Current speed', unit: 'm/s' }, y: { label: 'Ur (-) / amplitude (m)', unit: '-' },
          series: [{ name: 'Reduced velocity', x: rnd(us), y: rnd(sweep.map((s) => s.Ur)) }, { name: 'VIM amplitude (m)', x: rnd(us), y: rnd(sweep.map((s) => s.A)) }, { name: 'Lock-in onset', x: [0, round(us[us.length - 1])], y: [inp.UrLow, inp.UrLow], dash: true }, { name: 'Lock-in end', x: [0, round(us[us.length - 1])], y: [inp.UrHigh, inp.UrHigh], dash: true }],
          explain: 'The reduced velocity rises linearly with current. Between the two dashed lines the hull is taken as locked in and the amplitude jumps to the model-test value; look at where the design current sits in that band.' },
        { type: 'line', title: 'Offset versus current', x: { label: 'Current speed', unit: 'm/s' }, y: { label: 'Offset', unit: 'm' },
          series: [{ name: 'Mean in-line offset', x: rnd(us), y: rnd(sweep.map((s) => s.mean)) }, { name: 'Total (mean + VIM)', x: rnd(us), y: rnd(sweep.map((s) => s.total)) }, { name: 'Allowable', x: [0, round(us[us.length - 1])], y: [round(offAllow), round(offAllow)], dash: true }],
          explain: 'The mean offset grows with the square of current and steps up inside lock-in because the drag is amplified. The total must stay under the allowable line for every current up to the design value.' },
        { type: 'line', title: 'Long-term current exceedance', x: { label: 'Current speed', unit: 'm/s' }, y: { label: 'Probability of exceedance', unit: '-' },
          series: [{ name: 'Weibull climate', x: rnd(us), y: rnd(us.map((u) => 1 - cdf(u))) }, { name: 'Lock-in onset', x: [round(d.criticalSpeeds[0]), round(d.criticalSpeeds[0])], y: [0, 1], dash: true }],
          explain: 'Share of time the current exceeds each speed. The value at the lock-in onset line is the fraction of the year with VIM (less the small share beyond the end of the band).' },
      ],
      tables: [
        { title: 'VIM screening summary', columns: [{ key: 'U', label: 'Current', unit: 'm/s' }, { key: 'Ur', label: 'Ur' }, { key: 'lock', label: 'Lock-in' }, { key: 'A', label: 'Amplitude', unit: 'm' }, { key: 'Cd', label: 'Cd' }, { key: 'mean', label: 'Mean offset', unit: 'm' }, { key: 'tot', label: 'Total offset', unit: 'm' }],
          rows: sweep.filter((s, i) => i % 8 === 0 || i === sweep.length - 1).map((s) => ({ U: round(s.Ur * D / Tn, 3), Ur: round(s.Ur, 4), lock: s.lockIn ? 'yes' : 'no', A: round(s.A, 4), Cd: round(s.Cd, 4), mean: round(s.mean, 4), tot: round(s.total, 4) })),
          explain: 'One row per current speed: reduced velocity, whether the hull is inside the lock-in band, and the resulting amplitude, drag coefficient and offsets.' },
      ],
      notes: [
        `Sway natural period ${fmt(Tn)} s (mass ${fmt(m / 1e3)} t + added mass ${fmt(ma / 1e3)} t on ${fmt(inp.Kmoor)} kN/m). Lock-in between ${fmt(d.criticalSpeeds[0])} and ${fmt(d.criticalSpeeds[1])} m/s.`,
        `At ${fmt(inp.Uc)} m/s: Ur = ${fmt(d.Ur)}, ${d.lockIn ? `locked in, amplitude ${fmt(d.A)} m (A/D ${fmt(aRatio)} after ${fmt(100 * inp.strakeEff)} % strake reduction), drag coefficient ${fmt(d.Cd)}` : 'outside the lock-in band'}.`,
        `The lock-in band and A/D are inputs from model tests, not predictions. Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at utilisation ${fmt(gov.uc)}; total offset ${fmt(env.total)} m against ${fmt(offAllow)} m, ${fmt(hours)} h/yr in lock-in.`,
          gov.uc > 1 ? `Reduce the amplitude ratio to about ${fmt(Math.max(Math.sqrt(Math.max(offAllow ** 2 - env.mean ** 2, 0)) / D, 0))} (better strake coverage or height), or stiffen the mooring: a stiffness of ${fmt(inp.Kmoor * gov.uc)} kN/m cuts the mean offset in proportion and raises the onset current to ${fmt(d.criticalSpeeds[0] * Math.sqrt(gov.uc))} m/s.`
            : `Margin on offset is ${fmt(offAllow - env.total)} m. Confirm A/D and the band by towing-tank VIM tests of the appended hull at the governing headings, and add ${fmt(cycles)} cycles per year of ${fmt(2 * d.A || 2 * aRatio * D)} m range to the mooring fatigue analysis.`,
          'DNV-OS-E301 (2021) Ch.2 Sec.2 [2.8]; API RP 2SK (2005) Appendix H; DNV-RP-C205 (2019) §9.'),
        adv('economic', 'info', `Strakes are the usual VIM mitigation: three-start strakes 0.13·D high over the draft are about ${fmt(strakeSteel)} t of plate, and they raise the in-line drag (and mooring load) by 30–60 %.`,
          'Trade strake steel and added mooring load against the mooring fatigue life recovered; a VIM model-test campaign is cheap compared with one extra chain size on every line.',
          'Value-of-information for model tests; whole-life costing per ISO 15686-5 (mooring replacement interval).',
          `${steelCo2Text(strakeSteel)} Fewer VIM cycles extend chain life and defer a mooring replacement campaign with its anchor-handler fuel.`),
      ],
      publish: { 'vessel.vimAmplitude': d.A },
    };
  },
  verify: [
    { // m = 65 000 t; ma = 1.0·1025·(π·30²/4)·90 = 65 207 t; Tn = 2π·√(130.207e6/250e3) = 2π·22.822 = 143.39 s.
      // U = 1 m/s: Ur = 143.39/30 = 4.780 (inside 4–10); shedding period = 30/(0.2·1) = 150 s; A = 0.5·30 = 15 m (no strakes).
      // Cd = 0.7·(1 + 2.1·0.5) = 1.435; F = ½·1025·1.435·(30·90)·1² = 1985.7 kN; mean offset = 7.943 m; total = √(7.943² + 15²) = 16.97 m.
      name: 'Spar: natural period, reduced velocity, amplified drag and offset',
      source: 'DNV-RP-C205 (2019) §9.1 (Ur = U·T/D, Cd amplification) — hand calculation',
      inputs: { massT: 65000, Uc: 1, strakeEff: 0, Kmoor: 250 },
      expect: [{ kpi: 'Tn', value: 143.39, tol: 1e-3 }, { kpi: 'Ur', value: 4.780, tol: 1e-3 }, { kpi: 'sheddingT', value: 150, tol: 1e-9 }, { kpi: 'amplitude', value: 15, tol: 1e-9 }, { kpi: 'Fdrag', value: 1985.7, tol: 1e-3 }, { kpi: 'totalOffset', value: 16.97, tol: 2e-3 }] },
    { // Lock-in exposure: P = exp(−(U1/A)^k) − exp(−(U2/A)^k); with Tn = 143.39 s, U1 = 4·30/143.39 = 0.8369, U2 = 2.0922 m/s,
      // A = 0.4, k = 1.5: (2.0922)^1.5 = 3.0263 → e^−3.0263 = 0.04849; (5.2305)^1.5 = 11.962 → 6.4e-6; P = 0.04849 → 425.0 h/yr.
      name: 'Hours per year inside the lock-in band from the Weibull current climate',
      source: 'Two-parameter Weibull distribution — hand calculation',
      inputs: { massT: 65000, Kmoor: 250 },
      expect: [{ kpi: 'lockHours', value: 425.0, tol: 5e-3 }] },
    { name: 'No current: no shedding, no offset', source: 'Limiting case U = 0', inputs: { Uc: 0 },
      expect: [{ kpi: 'totalOffset', value: 0, abs: 1e-12 }, { kpi: 'amplitude', value: 0, abs: 1e-12 }] },
  ],
};
