// Mooring redundancy: every single-line failure in turn (damaged equilibrium), the transient after the worst break
// (time-domain surge/sway and work–energy overshoot) and a progressive-failure check on the remaining lines.
import { brokenLineAnalysis, checkApi2sk, checkDnvE301, checkOffset } from '../../engine/lines/mooringsystem.js';
import { DEG, need, run, round, fmt, adv, check, worst, sev, thin, systemInputs, loadInputs, criteriaInputs, meanState, steelAdvice, LIMIT_QS, LIMIT_PROPS } from './_common.js';

export default {
  id: 's10-line-break',
  suite: 10,
  title: 'Redundancy and line-break analysis',
  summary: 'Removes each line in turn, finds the new vessel position and tensions, follows the transient after the worst break and tests whether a second line would fail.',
  covers: [11],
  fidelity: 'reduced-order',
  standards: ['API RP 2SK (3rd ed., 2005) §5.5, §6.2–6.3 (damaged and transient conditions)', 'DNV-OS-E301 (2018) Ch.2 Sec.2 [4.3] (ALS)', 'ISO 19901-7 (2013) §8.9 (redundancy check)'],
  models: ['Quasi-static catenary system with one line removed', 'Newtonian surge–sway transient of the vessel on the remaining lines (RK4)', 'Work–energy overshoot estimate (first principles)', 'API safety factors and DNV ALS partial factors for the damaged condition', 'Progressive (cascading) failure check: second line removed when its transient tension exceeds its break load'],
  limits: `${LIMIT_QS} The transient uses constant mass and added mass, linear damping and quasi-static lines with the heading held; wave-frequency motion is not superposed on the transient, and the overshoot from the energy method is an undamped upper bound along a straight path. Dynamic tension is added as an input ratio on the mean tension. ${LIMIT_PROPS}`,
  publishes: ['mooring.TmaxDamaged', 'mooring.offsetDamaged'],
  searchTerms: 'mooring line failure transient redundancy damaged condition',
  inputs: [
    ...systemInputs(),
    ...loadInputs(),
    { key: 'breakLine', label: 'Line to break (0 = find the worst)', type: 'int', default: 0, min: 0, max: 24, group: 'Damage case', role: 'bc',
      help: 'Boundary condition of the damage case: this line\'s load-carrying capability is removed at the event time. 0 breaks each line in turn and follows the one that leaves the highest tension.' },
    { key: 'tBreak', label: 'Time of the break', unit: 's', type: 'number', default: 20, min: 0, max: 600, group: 'Damage case', role: 'ic',
      help: 'The vessel is at rest at its intact equilibrium under the mean load until this time (initial condition), then the line is released.' },
    { key: 'dispT', label: 'Displacement', unit: 't', type: 'number', default: 50000, min: 500, max: 600000, group: 'Vessel', role: 'geometry', bind: 'vessel.disp' },
    { key: 'Ca', label: 'Horizontal added-mass coefficient', unit: '-', type: 'number', default: 0.15, min: 0, max: 2, group: 'Vessel', role: 'load' },
    { key: 'zeta', label: 'Low-frequency damping', unit: '% critical', type: 'number', default: 15, min: 0, max: 100, group: 'Vessel', role: 'load', uncertain: { dist: 'lognormal', cov: 0.3 },
      help: 'Damping of the slow surge/sway motion on the damaged mooring. 0 reproduces the undamped energy-method overshoot.' },
    { key: 'dynRatio', label: 'Dynamic / mean tension increment ratio', unit: '-', type: 'number', default: 0.4, min: 0, max: 3, group: 'Vessel', role: 'load', bind: 'mooring.dynRatio',
      help: 'Low- and wave-frequency tension added to the damaged-condition mean tension, as a fraction of the mean tension (from the intact tension-and-offset analysis: Tdyn / Tmean).' },
    { key: 'duration', label: 'Transient duration', unit: 'natural periods', type: 'number', default: 4, min: 1.5, max: 12, group: 'Numerical', role: 'control' },
    { key: 'stepsPerPeriod', label: 'Time steps per natural period', type: 'int', default: 80, min: 20, max: 640, group: 'Numerical', role: 'control', help: 'RK4 steps per surge natural period of the damaged system.' },
    ...criteriaInputs({ damaged: true, offsetDefault: 12 }),
  ],
  meshTime: { key: 'stepsPerPeriod', label: 'Steps per natural period', levels: [40, 80, 160], qoi: 'transientOffset' },
  compute(inp) {
    const { b, sys, fs, Fext, calm, eq } = meanState(inp);
    const n = b.nLines;
    need(inp.breakLine <= n, `Line ${inp.breakLine} does not exist: the system has ${n} lines.`);
    // --- every single-line failure: damaged static equilibrium
    const cases = [];
    for (let k = 0; k < n; k++) {
      try {
        const d = fs.solve(Fext, [eq.x, eq.y], [k]);
        cases.push({ k, d, offset: Math.hypot(d.x - calm.x, d.y - calm.y), shift: Math.hypot(d.x - eq.x, d.y - eq.y), Tmax: d.maxTension, line: d.maxTensionLine, ok: true });
      } catch (e) { cases.push({ k, ok: false, offset: Infinity, shift: Infinity, Tmax: Infinity, line: -1 }); }
    }
    const lost = cases.filter((c) => !c.ok);
    need(lost.length < n, 'No single-line-broken case reaches equilibrium: the system has no redundancy at this load.');
    const pick = inp.breakLine > 0 ? cases[inp.breakLine - 1] : cases.filter((c) => c.ok).reduce((a, c) => (c.Tmax > a.Tmax ? c : a));
    need(pick.ok, `With ${sys.lines[pick.k].name} broken the remaining lines cannot balance the load: station keeping is lost in this damage case.`);
    const kB = pick.k;
    // --- transient: M q'' + C q' = F_lines(q) + Fext in surge and sway, heading held
    const M = inp.dispT * 1e3 * (1 + inp.Ca);
    const h0 = 0.05;
    const fx0 = fs.forces(pick.d.x, pick.d.y, [kB]), fx1 = fs.forces(pick.d.x + h0, pick.d.y, [kB]), fy1 = fs.forces(pick.d.x, pick.d.y + h0, [kB]);
    const dir = pick.shift > 1e-6 ? Math.atan2(pick.d.y - eq.y, pick.d.x - eq.x) : inp.loadDirDeg * DEG;
    const cd = Math.cos(dir), sd = Math.sin(dir);
    const Kxx = -(fx1.F[0] - fx0.F[0]) / h0, Kyy = -(fy1.F[1] - fx0.F[1]) / h0, Kxy = -(fy1.F[0] - fx0.F[0]) / h0;
    const Kd = Math.max(Kxx * cd * cd + 2 * Kxy * cd * sd + Kyy * sd * sd, 1);          // stiffness along the excursion
    const Tn = 2 * Math.PI * Math.sqrt(M / Kd);
    const C = 2 * (inp.zeta / 100) * Math.sqrt(Kd * M);
    const dt = Tn / inp.stepsPerPeriod, tEnd = inp.tBreak + inp.duration * Tn, N = Math.ceil(tEnd / dt);
    const acc = (t, y) => {
      const r = t >= inp.tBreak ? fs.forces(y[0], y[1], [kB]) : fs.forces(y[0], y[1]);
      return { a: [y[2], y[3], (r.F[0] + Fext[0] - C * y[2]) / M, (r.F[1] + Fext[1] - C * y[3]) / M], r };
    };
    const y = [eq.x, eq.y, 0, 0];
    const hist = { t: [], off: [], exc: [], T: [] };
    let peakT = 0, peakLine = -1, peakTime = 0, maxExc = 0, maxOff = 0;
    const lineMax = new Array(n).fill(0);
    let slackMin = Infinity;
    for (let i = 0; i <= N; i++) {
      const t = i * dt;
      let k1;
      try { k1 = acc(t, y); } catch (e) { throw new Error(`The transient after the break of ${sys.lines[kB].name} overloads a remaining line beyond its elastic range at t = ${fmt(t)} s: the system cannot arrest the vessel. Add damping, lines or strength.`); }
      const exc = Math.hypot(y[0] - eq.x, y[1] - eq.y), off = Math.hypot(y[0] - calm.x, y[1] - calm.y);
      hist.t.push(t); hist.off.push(off); hist.exc.push(exc); hist.T.push(k1.r.maxTension);
      maxExc = Math.max(maxExc, exc); maxOff = Math.max(maxOff, off);
      k1.r.lines.forEach((l, j) => { if (!l.broken) { if (l.T > lineMax[j]) lineMax[j] = l.T; if (t >= inp.tBreak) slackMin = Math.min(slackMin, l.T); } });
      if (t >= inp.tBreak && k1.r.maxTension > peakT) { peakT = k1.r.maxTension; peakLine = k1.r.maxTensionLine; peakTime = t; }
      if (i === N) break;
      const st = (kk, f) => y.map((v, j) => v + f * dt * kk.a[j]);
      const k2 = acc(t + dt / 2, st(k1, 0.5)), k3 = acc(t + dt / 2, st(k2, 0.5)), k4 = acc(t + dt, st(k3, 1));
      for (let j = 0; j < 4; j++) y[j] += dt / 6 * (k1.a[j] + 2 * k2.a[j] + 2 * k3.a[j] + k4.a[j]);
    }
    const overshoot = pick.shift > 1e-9 ? maxExc / pick.shift : 1;
    // --- engine work–energy overshoot for the same line (undamped bound)
    let energy = null;
    try { energy = brokenLineAnalysis(sys, { Fext, dofs: [0, 1], lines: [kB] }).damaged[0]; } catch (e) { energy = null; }
    // --- design tensions in the damaged condition (mean + dynamic allowance) and transient
    const TmeanD = pick.Tmax, TdynD = inp.dynRatio * TmeanD, TmaxD = TmeanD + TdynD;
    const apiD = checkApi2sk({ tension: TmaxD, MBL: b.MBL, condition: 'damaged', requiredSF: inp.sfDamaged });
    const als = checkDnvE301({ Tmean: TmeanD, Tdyn: TdynD, MBL: b.MBL, limitState: 'ALS', strengthFactor: inp.strengthFactor, gammaMean: inp.gammaMean, gammaDyn: inp.gammaDyn });
    const off = checkOffset({ offset: Math.max(maxOff, pick.offset), depth: b.depth, limitFraction: inp.offsetLimitPct / 100 });
    // --- progressive failure: if the transient peak exceeds the break load, remove that line too
    let second = null, cascade = 0;
    if (peakT * (1 + inp.dynRatio) > b.MBL) {
      cascade = 1;
      try { second = fs.solve(Fext, [pick.d.x, pick.d.y], [kB, peakLine]); } catch (e) { second = null; cascade = 2; }
    }
    let two = null;
    if (!cascade) { try { two = fs.solve(Fext, [pick.d.x, pick.d.y], [kB, pick.line]); } catch (e) { two = null; } }
    const sfTwo = two ? b.MBL / two.maxTension : 0;
    const checks = [
      check('Damaged-condition tension, API safety factor', TmaxD * apiD.required / 1e3, b.MBL / 1e3, 'kN', 'API RP 2SK (3rd ed.) §6.3 (damaged, dynamic)',
        'Guards against a second line parting once the vessel has settled at the one-line-broken equilibrium in the design storm.'),
      check('Damaged-condition tension, DNV ALS', als.Td / 1e3, als.Sc / 1e3, 'kN', 'DNV-OS-E301 Ch.2 Sec.2 [4.3]',
        'Accidental limit state with one line missing, in partial-factor format.'),
      check('Transient peak tension', peakT / 1e3, b.MBL / 1e3, 'kN', 'API RP 2SK (3rd ed.) §5.5 (transient analysis): no break in the transient',
        'Guards against the overshoot while the vessel swings to its new position: the peak may exceed the settled damaged tension and trigger a second failure before the vessel comes to rest.'),
      check('Damaged / transient offset', Math.max(maxOff, pick.offset), inp.offsetLimitPct / 100 * b.depth, 'm', 'API RP 2SK (3rd ed.) §6.2 (damaged-condition offset limit)',
        'Guards against riser, umbilical or neighbouring-structure damage during and after the excursion that follows a line failure.'),
    ];
    const g = worst(checks);
    const [pt, pe, po, pT] = thin([hist.t, hist.exc, hist.off, hist.T], 600);
    const sfDam = b.MBL / TmaxD;
    return {
      kpis: [
        { key: 'intactOffset', label: 'Intact mean offset', value: round(eq.offset, 5), unit: 'm', explain: 'Vessel offset under the mean load before the break; the starting point of the transient.' },
        { key: 'intactTmax', label: 'Intact largest mean tension', value: round(eq.maxTension / 1e3, 5), unit: 'kN', explain: 'Mean tension of the most loaded line with all lines present.' },
        { key: 'staticShift', label: 'Shift to the damaged equilibrium', value: round(pick.shift, 5), unit: 'm', explain: 'Distance between the intact and the one-line-broken equilibrium positions for the governing break.' },
        { key: 'damagedOffset', label: 'Post-failure mean offset', value: round(pick.offset, 5), unit: 'm', explain: 'Offset from the calm-water position at which the vessel settles with one line broken.' },
        { key: 'transientOffset', label: 'Maximum transient excursion', value: round(maxExc, 5), unit: 'm', explain: 'Largest distance from the intact position reached during the swing after the break; larger than the static shift because the vessel overshoots.' },
        { key: 'overshoot', label: 'Overshoot factor (time domain)', value: round(overshoot, 4), unit: '-', explain: 'Maximum transient excursion over the static shift. 2.0 is the undamped linear limit; damping and a stiffening mooring reduce it.' },
        { key: 'overshootEnergy', label: 'Overshoot factor (energy method, undamped)', value: round(energy ? energy.overshootFactor : overshoot, 4), unit: '-', explain: 'Independent estimate from the work–energy balance along a straight path without damping; an upper bound to compare with the time-domain value.' },
        { key: 'TmeanDamaged', label: 'Damaged mean tension', value: round(TmeanD / 1e3, 5), unit: 'kN', explain: 'Mean tension of the most loaded remaining line at the damaged equilibrium.' },
        { key: 'TmaxDamaged', label: 'Damaged maximum tension (with dynamics)', value: round(TmaxD / 1e3, 5), unit: 'kN', explain: 'Damaged mean tension plus the dynamic allowance; the design tension for the damaged-condition checks.' },
        { key: 'sfDamaged', label: 'Damaged-condition safety factor', value: round(sfDam, 4), unit: '-', limit: inp.sfDamaged, sense: 'min', explain: 'Break load over damaged maximum tension. API RP 2SK asks for 1.25 with dynamic analysis.' },
        { key: 'Tpeak', label: 'Transient peak tension', value: round(peakT / 1e3, 5), unit: 'kN', limit: round(b.MBL / 1e3, 5), sense: 'max', explain: 'Largest mean-level line tension during the swing; it must stay below the break load for the failure not to cascade.' },
        { key: 'Tn', label: 'Natural period, damaged system', value: round(Tn, 4), unit: 's', explain: 'Period of the swing on the remaining lines; the transient peak arrives about half a period after the break.' },
        { key: 'cascade', label: 'Progressive failure indicator', value: cascade, unit: '-', limit: 0, sense: 'max', explain: '0: remaining lines survive the transient. 1: a second line is predicted to part but the system finds a new equilibrium. 2: station keeping is lost.' },
        { key: 'sfTwoBroken', label: 'Safety factor with two lines broken', value: round(second ? b.MBL / second.maxTension : sfTwo, 4), unit: '-', explain: 'Break load over mean tension if the most loaded remaining line also fails (0 when no equilibrium exists). Indicates whether the system would survive a double failure.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Vessel excursion after the line break', x: { label: 'Time', unit: 's' }, y: { label: 'Distance', unit: 'm' },
          series: [{ name: 'Excursion from intact position', x: pt.map((v) => round(v)), y: pe.map((v) => round(v)) }, { name: 'Offset from calm-water position', x: pt.map((v) => round(v)), y: po.map((v) => round(v)) },
            { name: 'Static shift', x: [0, hist.t[hist.t.length - 1]], y: [pick.shift, pick.shift], dash: true }, { name: 'Offset limit', x: [0, hist.t[hist.t.length - 1]], y: [inp.offsetLimitPct / 100 * b.depth, inp.offsetLimitPct / 100 * b.depth], dash: true }],
          explain: 'The vessel is released at the break time, swings past its new equilibrium (dashed static shift) and decays onto it. The first peak is the transient excursion to hold against the offset limit.' },
        { type: 'line', title: 'Largest line tension during the transient', x: { label: 'Time', unit: 's' }, y: { label: 'Tension', unit: 'kN' },
          series: [{ name: 'Largest remaining-line tension', x: pt.map((v) => round(v)), y: pT.map((v) => round(v / 1e3)) }, { name: 'Break load', x: [0, hist.t[hist.t.length - 1]], y: [b.MBL / 1e3, b.MBL / 1e3], dash: true }, { name: 'Damaged equilibrium', x: [0, hist.t[hist.t.length - 1]], y: [TmeanD / 1e3, TmeanD / 1e3], dash: true }],
          explain: 'Tension jumps at the break as neighbours take the load, peaks at maximum excursion and settles at the damaged-equilibrium level. A peak touching the break-load line means a cascading failure.' },
        { type: 'bar', title: 'Single-line-failure cases', categories: cases.map((c) => `${sys.lines[c.k].name} broken`),
          series: [{ name: 'Largest remaining tension (kN)', y: cases.map((c) => (c.ok ? round(c.Tmax / 1e3) : 0)) }, { name: 'Post-failure offset (m × 100)', y: cases.map((c) => (c.ok ? round(c.offset * 100) : 0)) }],
          explain: 'Consequence of losing each line. Breaking a windward line is worst for tension; a zero bar marks a case with no equilibrium. Offsets are scaled by 100 to share the axis.' },
      ],
      tables: [
        { title: 'Single-line-failure response', columns: [{ key: 'broken', label: 'Broken line' }, { key: 'off', label: 'Post-failure offset', unit: 'm' }, { key: 'shift', label: 'Shift from intact', unit: 'm' }, { key: 'T', label: 'Largest remaining mean tension', unit: 'kN' }, { key: 'line', label: 'In line' }, { key: 'sf', label: 'MBL / (mean + dynamic)' }, { key: 'status', label: 'Status' }],
          rows: cases.map((c) => (c.ok ? { broken: sys.lines[c.k].name, off: round(c.offset), shift: round(c.shift), T: round(c.Tmax / 1e3), line: sys.lines[c.line].name, sf: round(b.MBL / (c.Tmax * (1 + inp.dynRatio)), 3), status: b.MBL / (c.Tmax * (1 + inp.dynRatio)) >= inp.sfDamaged ? 'meets factor' : 'below factor' }
            : { broken: sys.lines[c.k].name, off: 0, shift: 0, T: 0, line: '–', sf: 0, status: 'no equilibrium: station keeping lost' })),
          explain: 'Each row removes one line under the same mean load. Use it to identify critical lines and failure scenarios: the lowest safety factor is the governing damaged case.' },
        { title: 'Peak tension per line during the transient', columns: [{ key: 'line', label: 'Line' }, { key: 'T', label: 'Peak mean-level tension', unit: 'kN' }, { key: 'u', label: 'Peak / MBL', unit: '%' }],
          rows: sys.lines.map((l, j) => ({ line: l.name + (j === kB ? ' (broken)' : ''), T: round(lineMax[j] / 1e3), u: round(100 * lineMax[j] / b.MBL, 3) })),
          explain: 'Largest tension each line sees from before the break to the end of the simulation. The neighbours of the broken line take the largest increase.' },
      ],
      notes: [`Governing break: ${sys.lines[kB].name}. Vessel shifts ${fmt(pick.shift)} m statically and ${fmt(maxExc)} m at the transient peak (overshoot ${fmt(overshoot)}; undamped energy method ${fmt(energy ? energy.overshootFactor : NaN)}).`,
        `${sys.lines[peakLine >= 0 ? peakLine : pick.line].name} peaks at ${fmt(peakT / 1e3)} kN at t = ${fmt(peakTime)} s; damaged natural period ${fmt(Tn)} s.`,
        cascade === 0 ? `No cascade: with ${sys.lines[pick.line].name} also removed the safety factor on mean tension would be ${fmt(sfTwo)}.` : cascade === 1 ? 'A second line is predicted to part in the transient; the system then finds a two-lines-broken equilibrium.' : 'A second line is predicted to part and no equilibrium exists afterwards: loss of station keeping.'],
      advice: [
        adv('engineering', cascade ? 'fail' : sev(g.uc), `Losing ${sys.lines[kB].name} moves the vessel ${fmt(pick.shift)} m (transient ${fmt(maxExc)} m) and raises the largest tension from ${fmt(eq.maxTension / 1e3)} to ${fmt(TmeanD / 1e3)} kN mean, ${fmt(peakT / 1e3)} kN transient; damaged safety factor ${fmt(sfDam)} against ${inp.sfDamaged}. Governing: ${g.name} (${fmt(g.uc)}).`,
          g.uc > 1 || cascade ? (g.name.includes('offset') ? `Limit the excursion: closer line spacing around the dominant storm direction or higher pretension; ${fmt(Math.max(maxOff, pick.offset) - inp.offsetLimitPct / 100 * b.depth)} m must be removed, or the risers given a damaged-condition allowance.` : `Add redundancy where it matters: an extra line next to ${sys.lines[kB].name}, or ${fmt(100 * (Math.sqrt(Math.max(g.uc, 1.01)) - 1))} % more diameter on the windward lines; thruster assistance is the alternative for the transient.`)
            : `The system is single-failure tolerant with ${fmt(100 * (1 - g.uc))} % reserve. Keep a spare line segment and connectors available so the damaged condition is short-lived.`,
          'API RP 2SK (3rd ed., 2005) §5.5 and §6.3; DNV-OS-E301 (2018) Ch.2 Sec.2 [4.3]; ISO 19901-7 (2013) §8.9. Confirm factors against the governing edition.'),
        adv('economic', 'info', `Redundancy is bought with line count: this system has ${n} lines carrying ${fmt(b.steel_t)} t; one line is ${fmt(b.steel_t / n)} t plus an anchor and a day of installation.`,
          g.uc > 1 || cascade ? `One additional line (≈ ${fmt(b.steel_t / n)} t) is usually cheaper than up-sizing all ${n} lines by ${fmt(100 * (Math.max(g.uc, 1.01) - 1))} % (≈ ${fmt(b.steel_t * (Math.max(g.uc, 1.01) - 1))} t); price both.` : `A line failure would stop production until repaired: weigh the cost of holding a spare line and a call-off vessel contract against the downtime of a damaged-condition weather restriction.`,
          'Value-of-redundancy assessment with whole-life costing per ISO 15686-5; spare-holding strategy per ISO 19901-7 (2013) §14 (in-service integrity management).',
          `An extra line adds about ${fmt(b.steel_t / n)} t of steel ≈ ${fmt(b.steel_t / n * 1.9)} t CO₂e (screening factor 1.9 t CO₂e/t) and one more anchor footprint, against the environmental consequence of a drifting unit and damaged risers.`),
        steelAdvice(b, Math.max(checks[0].uc, checks[1].uc), 'damaged-condition tension'),
      ],
      publish: { 'mooring.TmaxDamaged': TmaxD, 'mooring.offsetDamaged': Math.max(maxOff, pick.offset) },
      quantities: { steel_t: round(b.steel_t) },
    };
  },
  verify: [
    { // Eight taut, practically weightless spring lines (EA = 100 MN, L = 990 m, span 1000 m, fairlead 1 m above the anchor level):
      // stretch s0 = 1000.0005 − 990 = 10.0005 m, k = EA/L = 101.01 kN/m, pretension T0 = k·s0 = 1010.15 kN. No external load.
      // Break the line at 0°. Remaining stiffness along x: in-line line k, four lines at ±45°/±135°: 4·(½k + ½T0/ℓ), two transverse
      // lines 2·T0/ℓ → K = 3k + 4T0/ℓ = 303.03 + 4.04 = 307.07 kN/m. Unbalanced force T0 → static shift δ = 1010.15/307.07 = 3.290 m.
      // Linear and undamped → the transient peak is 2δ = 6.579 m (overshoot factor 2).
      name: 'Taut linear-spring spread: static shift and undamped overshoot after one line fails',
      source: 'Hand calculation (linear springs with geometric stiffness; step release of a linear undamped oscillator overshoots by 100 %)',
      inputs: { nLines: 8, headingDeg: 0, depth: 30, fairleadDraft: 29, fairleadRadius: 1, anchorRadius: 1001, CB: 0, Fmean: 0, breakLine: 1, zeta: 0, stepsPerPeriod: 160, duration: 2, dynRatio: 0,
        segments: [{ name: 'spring', type: 'custom', d: 0, grade: '', L: 990, w: 0.01, EA: 100, MBL: 20000 }] },
      expect: [{ kpi: 'staticShift', value: 3.290, tol: 0.005 }, { kpi: 'overshoot', value: 2.0, tol: 0.01 }, { kpi: 'overshootEnergy', value: 2.0, tol: 0.02 }, { kpi: 'intactTmax', value: 1010.15, tol: 0.002 }] },
  ],
};
