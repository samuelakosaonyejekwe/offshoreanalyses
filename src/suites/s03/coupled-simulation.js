// Integrated aero-hydro-servo-elastic time-domain simulation of a bottom-fixed offshore wind turbine: rotor and
// controller, modal support structure on foundation springs, turbulent wind and irregular waves, with optional
// shutdown and gust events.
import { simulateCoupled } from '../../engine/aero/coupled.js';
import { rotorEffectiveTurbulence, ntmSigma, etmSigma, ewm, eogGust, eogSpeed, turbulenceScale } from '../../engine/aero/turbulence.js';
import { need, run, fmt, round, adv, check, worst, sev, steelCo2Text, thin, tube, turbineInputs, supportInputs, foundationInputs, steelInputs, steelOf, turbine, operatingPoint, LIMIT_AERO, LIMIT_COUPLED, LIMIT_CODE, LIMIT_HIFI, CONFIRM } from './_common.js';
import { dynamicsInputs, supportStructure } from './_dynamics.js';

const DT_WIND = 0.1; // the wind record is synthesised on a fixed grid so that time-step refinement sees the same inflow

export default {
  id: 's03-coupled-simulation', suite: 3, title: 'Integrated aero-hydro-servo-elastic simulation',
  summary: 'Time-domain response of rotor, controller, tower and monopile to turbulent wind, irregular waves and current, including emergency stop, grid loss and gust events, with the resulting loads at the tower base and the mudline.',
  covers: [1], fidelity: 'reduced-order',
  standards: ['IEC 61400-3-1:2019 §7 (design situations and load cases)', 'IEC 61400-1:2019 §6.3 (wind conditions), §7.4', 'DNV-ST-0437 (2016) §4 (load calculations)', 'DNV-ST-0126 (2021) §4', 'DNV-RP-C205 (2019) §3, §6', 'Jonkman et al. (2009) NREL/TP-500-38060'],
  models: ['Blade element momentum theory with Prandtl tip/hub loss and Glauert (Buhl) high-thrust correction (coefficient surfaces)', 'Actuator-disc momentum theory and Øye dynamic inflow (optional)', 'Variable-speed torque law and gain-scheduled PI blade-pitch control', 'Kaimal turbulence with rotor averaging; IEC NTM, ETM and extreme operating gust', 'Linear irregular waves (JONSWAP), Morison loads with relative velocity, MacCamy–Fuchs diffraction, current drag', 'Modal support structure on three-spring foundation (coupled equations of motion, RK4)', 'Rainflow counting and damage-equivalent loads'],
  limits: `${LIMIT_COUPLED} ${LIMIT_AERO} Actuator-line aerodynamics and Mann uniform-shear turbulence are not available: blade element momentum surfaces and the Kaimal model are used. Short records with one seed give indicative extremes only: IEC 61400-1 asks for at least six 10-minute realisations per condition. ${LIMIT_CODE} ${LIMIT_HIFI}`,
  searchTerms: 'aero-hydro-servo-elastic simulation offshore wind turbine monopile reduced order',
  inputs: [
    ...turbineInputs(),
    ...supportInputs({ embed: false }),
    ...steelInputs(),
    ...foundationInputs(),
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'windModel', label: 'Wind condition', type: 'select', default: 'ntm', group: 'Wind', role: 'bc',
      options: [{ value: 'ntm', label: 'Normal turbulence (NTM)' }, { value: 'etm', label: 'Extreme turbulence (ETM, class I mean wind)' }, { value: 'steady', label: 'Steady wind' }, { value: 'eog', label: 'Extreme operating gust (EOG) in steady wind' }] },
    { key: 'Vhub', label: 'Mean wind speed at hub height', unit: 'm/s', type: 'number', default: 12, min: 2, max: 70, group: 'Wind', role: 'bc', uncertain: { dist: 'normal', cov: 0.05 } },
    { key: 'Iref', label: 'Reference turbulence intensity', unit: '-', type: 'number', default: 0.14, min: 0.05, max: 0.25, group: 'Wind', role: 'bc', help: 'IEC class A 0.16, B 0.14, C 0.12.' },
    { key: 'seed', label: 'Random seed (wind and waves)', type: 'int', default: 1, min: 1, max: 99999, group: 'Wind', role: 'ic', help: 'Fixes the turbulence and wave realisation so that runs are reproducible.' },
    { key: 'Hs', label: 'Significant wave height', unit: 'm', type: 'number', default: 2.5, min: 0, max: 15, group: 'Waves and current', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.1 } },
    { key: 'Tp', label: 'Peak wave period', unit: 's', type: 'number', default: 7, min: 2, max: 22, group: 'Waves and current', role: 'bc' },
    { key: 'gamma', label: 'JONSWAP peak enhancement', unit: '-', type: 'number', default: 3.3, min: 1, max: 7, group: 'Waves and current', role: 'bc' },
    { key: 'Uc', label: 'Current speed (uniform)', unit: 'm/s', type: 'number', default: 0.5, min: 0, max: 3, group: 'Waves and current', role: 'bc', bind: 'site.Uc1' },
    ...dynamicsInputs(),
    { key: 'state', label: 'Turbine state at the start', type: 'select', default: 'production', group: 'Operation', role: 'ic',
      options: [{ value: 'production', label: 'Power production at the steady operating point' }, { value: 'parked', label: 'Parked, blades feathered' }],
      help: 'Initial rotor speed and pitch. Outside the cut-in to cut-out range the turbine is always parked.' },
    { key: 'event', label: 'Event during the run', type: 'select', default: 'none', group: 'Operation', role: 'operation',
      options: [{ value: 'none', label: 'None' }, { value: 'stop', label: 'Normal shutdown' }, { value: 'estop', label: 'Emergency stop' }, { value: 'gridloss', label: 'Loss of grid (emergency stop without generator torque)' }] },
    { key: 'tEvent', label: 'Time of the event (or gust start)', unit: 's', type: 'number', default: 80, min: 1, max: 3000, group: 'Operation', role: 'operation' },
    { key: 'inflow', label: 'Dynamic inflow', type: 'bool', default: false, group: 'Operation', role: 'control', help: 'Lags the rotor induction after fast pitch or wind changes (Øye time constant).' },
    { key: 'towerDrag', label: 'Wind drag on the tower', type: 'bool', default: true, group: 'Operation', role: 'load' },
    { key: 'gammaF', label: 'Load factor on the simulated extreme', unit: '-', type: 'number', default: 1.35, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: `IEC 61400-3-1 normal design situations (1.1 abnormal). ${CONFIRM}` },
    { key: 'gammaM', label: 'Material factor on yield', unit: '-', type: 'number', default: 1.1, min: 1, max: 1.5, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'duration', label: 'Simulated time after the start-up transient', unit: 's', type: 'number', default: 120, min: 30, max: 3600, group: 'Numerics', role: 'control', help: 'Longer records (600 s and more) give better statistics but take proportionally longer.' },
    { key: 'dt', label: 'Time step', unit: 's', type: 'number', default: 0.05, min: 0.005, max: 0.1, group: 'Numerics', role: 'control' },
    { key: 'nElem', label: 'Beam elements of the support structure', type: 'int', default: 30, min: 8, max: 160, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Beam elements', levels: [15, 30, 60], qoi: 'f1' },
  meshTime: { key: 'dt', label: 'Time step', levels: [0.1, 0.05, 0.025], isSize: true, qoi: 'mudMomentMax' },
  compute(inp, ctx) {
    const steel = steelOf(inp, ctx);
    const tb = turbine(inp);
    const str = supportStructure(inp, tb, steel, { nElem: inp.nElem, nModes: 2 });
    const f1 = str.modes[0].frequency, f2 = str.modes[1].frequency;
    need(2 * Math.PI * f2 * inp.dt < 1.9, `The time step ${inp.dt} s is too large for the second structural mode (${fmt(f2)} Hz): use ${fmt(1.5 / (2 * Math.PI * f2), 2)} s or less.`);
    const op = operatingPoint(tb, inp.Vhub);
    const parked = inp.state === 'parked' || op.parked;
    const tT = 20, T = inp.duration + tT;
    const tEv = tT + inp.tEvent;
    // Wind record (rotor-effective speed).
    let sigma1 = 0, vGust = 0, wind;
    const R = tb.R;
    if (inp.windModel === 'ntm' || inp.windModel === 'etm') {
      sigma1 = inp.windModel === 'ntm' ? ntmSigma({ Vhub: inp.Vhub, Iref: inp.Iref }) : etmSigma({ Vhub: inp.Vhub, Iref: inp.Iref, Vave: 10 });
      if (parked && inp.Vhub > tb.rotor.cutOut) sigma1 = 0.11 * inp.Vhub; // IEC extreme wind model turbulence
      const u = run(() => rotorEffectiveTurbulence({ Vhub: inp.Vhub, sigma1, hubHeight: inp.hubHeight, R, duration: T + 1, dt: DT_WIND, seed: inp.seed }), 'Turbulence synthesis').u;
      wind = (t) => { const x = Math.min(Math.max(t / DT_WIND, 0), u.length - 1.000001), i = Math.floor(x); return Math.max(0.1, u[i] + (x - i) * (u[i + 1] - u[i])); };
    } else if (inp.windModel === 'eog') {
      const s1 = ntmSigma({ Vhub: inp.Vhub, Iref: inp.Iref });
      vGust = eogGust({ Vhub: inp.Vhub, sigma1: s1, D: 2 * R, lambda1: turbulenceScale(inp.hubHeight), Ve1: ewm({ Vref: 50, z: inp.hubHeight, hubHeight: inp.hubHeight }).Ve1 });
      wind = (t) => eogSpeed(t - tEv, { Vhub: inp.Vhub, Vgust: vGust });
    } else wind = inp.Vhub;
    const cmd = { stop: 'stop', estop: 'estop', gridloss: 'estop' }[inp.event];
    const events = !parked && cmd ? (t) => (t >= tEv ? { command: cmd, gridLoss: inp.event === 'gridloss' } : null) : undefined;
    const args = {
      surface: tb.surface, turbine: { inertia: tb.rotor.drivetrainInertia, hubHeight: inp.hubHeight }, structure: str,
      wind, waves: inp.Hs > 0 || inp.Uc > 0 ? { Hs: Math.max(inp.Hs, 1e-4), Tp: inp.Tp, gamma: inp.gamma, seed: inp.seed + 500, depth: inp.depth, macCamyFuchs: inp.mcf, current: inp.Uc } : null,
      dt: inp.dt, duration: T, tTransient: tT, rho: tb.rhoAir, rhoWater: inp.rhoW, shearExp: 0.14, towerWindDrag: inp.towerDrag, dynamicInflow: inp.inflow, events,
      ...(parked ? { controller: null, initial: { omega: 0, pitchDeg: 90 }, aeroOverride: { ct: tb.parkedCt, cq: 0 } } : { controller: tb.controller, initial: { omega: op.omega, pitchDeg: op.pitchDeg } }),
    };
    const r = run(() => simulateCoupled(args), 'Coupled simulation');
    const st = r.stats;
    for (const k of ['mudlineMoment', 'towerBaseMoment', 'thrust', 'towerTopDisp']) need(Number.isFinite(st[k].absMax), 'The simulation diverged (non-finite response): reduce the time step.');
    const pile = tube(inp.D, inp.t / 1e3), tw = tube(inp.towerDb, inp.towerTb / 1e3);
    const fyd = steel.fy / inp.gammaM;
    const sigMud = inp.gammaF * st.mudlineMoment.absMax / pile.W, sigTB = inp.gammaF * st.towerBaseMoment.absMax / tw.W;
    const rpm = (v) => (v * 30) / Math.PI;
    const overspeed = rpm(st.rotorSpeed.max) / tb.rpmRated;
    const checks = [
      check('Monopile bending stress at the mudline (factored extreme)', sigMud / 1e6, fyd / 1e6, 'MPa', 'IEC 61400-3-1:2019 §7.6; DNV-ST-0126 (2021) §4.4', 'Bending stress from the largest simulated mudline moment times the load factor. Exceedance means the pile wall yields in this design situation.'),
      check('Tower-base bending stress (factored extreme)', sigTB / 1e6, fyd / 1e6, 'MPa', 'IEC 61400-3-1:2019 §7.6; DNV-ST-0126 (2021) §4.4', 'Bending stress at the tower bottom flange level from the largest simulated moment times the load factor. Shell buckling is checked separately.'),
      check('Rotor overspeed', overspeed, 1.2, '× rated', 'IEC 61400-1:2019 §8 (protection functions); typical supplier trip level 1.15–1.25 × rated', 'Guards against the rotor running away during gusts or grid loss. Above the trip level the safety system must stop the turbine and blade and drive-train loads rise sharply.'),
    ];
    const gov = worst(checks), ok = gov.uc <= 1;
    const i0 = Math.round(tT / inp.dt);
    const cut = (a) => a.subarray(i0);
    const [tt, w, th, mm, tbm, xd, om, pi, pw, eta] = thin([cut(r.t), cut(r.wind), cut(r.thrust), cut(r.mudlineMoment), cut(r.towerBaseMoment), cut(r.towerTopDisp), cut(r.rotorSpeed), cut(r.pitchDeg), cut(r.power), cut(r.waveElevation)], 900);
    const t0 = tt.map((v) => round(v - tT, 5));
    const psd = r.psd ? thin([r.psd.f, r.psd.towerBaseMoment], 500) : [[], []];
    const kMax = psd[0].findIndex((f) => f > Math.max(1.2, 1.5 * f2));
    const pf = kMax > 0 ? psd[0].slice(1, kMax) : psd[0].slice(1), pp = (kMax > 0 ? psd[1].slice(1, kMax) : psd[1].slice(1)).map((v) => v / 1e12);
    const lever = st.thrust.mean !== 0 ? st.mudlineMoment.mean / st.thrust.mean : 0;
    const tReq = inp.t * sigMud / fyd;
    const dT = (tube(inp.D, Math.min(Math.max(tReq, 15), 0.4 * inp.D * 1e3) / 1e3).A - pile.A) * steel.rho * (inp.depth + inp.zInt) / 1e3;
    const chans = [['Wind speed', 'wind', 1, 'm/s'], ['Wave elevation', 'waveElevation', 1, 'm'], ['Rotor speed', 'rotorSpeed', 30 / Math.PI, 'rpm'], ['Blade pitch', 'pitchDeg', 1, 'deg'], ['Electrical power', 'power', 1e-6, 'MW'], ['Rotor thrust', 'thrust', 1e-3, 'kN'],
      ['Tower-top displacement', 'towerTopDisp', 1, 'm'], ['Tower-top acceleration', 'towerTopAcc', 1, 'm/s²'], ['Tower-base moment', 'towerBaseMoment', 1e-6, 'MN·m'], ['Tower-base shear', 'towerBaseShear', 1e-6, 'MN'], ['Mudline moment', 'mudlineMoment', 1e-6, 'MN·m'], ['Mudline shear', 'mudlineShear', 1e-6, 'MN']];
    return {
      kpis: [
        { key: 'f1', label: 'First fore-aft natural frequency', unit: 'Hz', value: f1, explain: 'Fundamental frequency of the support-structure model used in the simulation; compare with the natural-frequency module.' },
        { key: 'powerMean', label: 'Mean electrical power', unit: 'MW', value: st.power.mean / 1e6, explain: 'Average power over the record. At or above rated wind it should equal the rated power; a shortfall shows the controller is saturating or the turbine stopped.' },
        { key: 'thrustMean', label: 'Mean rotor thrust', unit: 'kN', value: st.thrust.mean / 1e3, explain: 'Average aerodynamic thrust on the rotor, the main source of overturning moment in production.' },
        { key: 'thrustMax', label: 'Maximum rotor thrust', unit: 'kN', value: st.thrust.max / 1e3, explain: 'Largest thrust in the record. Around rated wind it exceeds the steady rated thrust because pitch control lags gusts.' },
        { key: 'mudMomentMax', label: 'Maximum mudline moment', unit: 'MN·m', value: st.mudlineMoment.absMax / 1e6, explain: 'Largest overturning moment at the mudline in the record (characteristic, unfactored): wind thrust, tower drag, wave and current load and dynamic amplification together.' },
        { key: 'mudMomentMean', label: 'Mean mudline moment', unit: 'MN·m', value: st.mudlineMoment.mean / 1e6, explain: 'Average overturning moment. The ratio of maximum to mean shows how dynamic the load case is.' },
        { key: 'mudMomentStd', label: 'Std. deviation of mudline moment', unit: 'MN·m', value: st.mudlineMoment.std / 1e6, explain: 'Fluctuating part of the mudline moment, the driver of fatigue damage.' },
        { key: 'mudShearMax', label: 'Maximum mudline shear', unit: 'MN', value: st.mudlineShear.absMax / 1e6, explain: 'Largest horizontal force at the mudline; with the moment it loads the pile–soil system.' },
        { key: 'mudShearMean', label: 'Mean mudline shear', unit: 'MN', value: st.mudlineShear.mean / 1e6, explain: 'Average horizontal force at the mudline: mean thrust, tower drag and current drag. Waves contribute to the fluctuation, not to the mean.' },
        { key: 'tbMomentMax', label: 'Maximum tower-base moment', unit: 'MN·m', value: st.towerBaseMoment.absMax / 1e6, explain: 'Largest bending moment at the tower base (interface level), which sizes the bottom flange and the transition piece.' },
        { key: 'dispMax', label: 'Maximum tower-top displacement', unit: 'm', value: st.towerTopDisp.absMax, explain: 'Largest fore-aft displacement at hub height, to compare with blade-to-tower clearance.' },
        { key: 'accMax', label: 'Maximum tower-top acceleration', unit: 'm/s²', value: st.towerTopAcc.absMax, explain: 'Largest fore-aft nacelle acceleration; turbines alarm at roughly 0.15–0.3 g.' },
        { key: 'rpmMax', label: 'Maximum rotor speed', unit: 'rpm', value: rpm(st.rotorSpeed.max), limit: 1.2 * tb.rpmRated, sense: 'max', explain: `Largest rotor speed in the record against 120 % of rated (${fmt(tb.rpmRated)} rpm).` },
        { key: 'delMud', label: 'Damage-equivalent mudline moment (1 Hz, m = 4)', unit: 'MN·m', value: r.del.mudlineMoment / 1e6, explain: 'Constant-amplitude moment range that, applied once per second, does the same fatigue damage as the record for an S–N slope of 4. Used to compare load cases, not as a stress.' },
        { key: 'leverArm', label: 'Mean mudline moment / mean thrust', unit: 'm', value: lever, explain: 'Effective height of the mean horizontal load above the mudline. With thrust only it equals hub height plus water depth; tower drag and current add moment without adding rotor thrust and so raise the ratio.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Wind speed and rotor thrust', x: { label: 'Time', unit: 's' }, y: { label: 'Wind (m/s) / thrust (×100 kN)', unit: '' }, series: [{ name: 'Rotor-effective wind, m/s', x: t0, y: w }, { name: 'Thrust, ×100 kN', x: t0, y: th.map((v) => v / 1e5) }],
          explain: 'Inflow and the resulting thrust. Below rated wind the thrust follows the wind; above rated the pitch controller reduces thrust as the wind rises.' },
        { type: 'line', title: 'Mudline and tower-base bending moment', x: { label: 'Time', unit: 's' }, y: { label: 'Moment', unit: 'MN·m' }, series: [{ name: 'Mudline', x: t0, y: mm.map((v) => v / 1e6) }, { name: 'Tower base', x: t0, y: tbm.map((v) => v / 1e6) }],
          explain: 'Overturning moment histories. Look for the largest peak, for wave-frequency oscillation and for ringing at the natural frequency after an event.' },
        { type: 'line', title: 'Rotor speed, pitch and power', x: { label: 'Time', unit: 's' }, y: { label: 'rpm / deg / MW', unit: '' }, series: [{ name: 'Rotor speed, rpm', x: t0, y: om.map(rpm) }, { name: 'Blade pitch, deg', x: t0, y: pi }, { name: 'Power, MW', x: t0, y: pw.map((v) => v / 1e6) }],
          explain: 'Controller behaviour: speed should stay near rated above rated wind while pitch moves. After a stop the pitch runs to feather and speed decays.' },
        { type: 'line', title: 'Tower-top displacement and wave elevation', x: { label: 'Time', unit: 's' }, y: { label: 'm', unit: 'm' }, series: [{ name: 'Tower-top displacement', x: t0, y: xd }, { name: 'Wave elevation', x: t0, y: eta }],
          explain: 'Fore-aft motion of the nacelle beside the sea surface elevation at the pile. Lightly damped oscillation at the natural period shows resonance.' },
        { type: 'line', title: 'Spectrum of the tower-base moment', x: { label: 'Frequency', unit: 'Hz' }, y: { label: 'Spectral density', unit: '(MN·m)²/Hz' }, series: [{ name: 'Tower-base moment', x: pf, y: pp }, { name: 'f1', x: [f1, f1], y: [0, Math.max(...pp, 0)], dash: true }, { name: 'Wave peak', x: [1 / inp.Tp, 1 / inp.Tp], y: [0, Math.max(...pp, 0)], dash: true }],
          explain: 'Frequency content of the response: low-frequency turbulence, the wave band and the structural resonance. The relative heights show which mechanism drives fatigue.' },
      ],
      tables: [
        { title: 'Channel statistics', columns: [{ key: 'c', label: 'Channel' }, { key: 'mean', label: 'Mean' }, { key: 'std', label: 'Std. dev.' }, { key: 'min', label: 'Min' }, { key: 'max', label: 'Max' }, { key: 'u', label: 'Unit' }],
          rows: chans.map(([c, k, f, u]) => ({ c, mean: round(st[k].mean * f), std: round(st[k].std * f), min: round(st[k].min * f), max: round(st[k].max * f), u })),
          explain: 'Statistics of each output channel after the start-up transient. Use the maxima as characteristic loads of this realisation and the standard deviations for fatigue comparison.' },
        ...(r.psd ? [{ title: 'Spectral peaks of the tower-base moment', columns: [{ key: 'n', label: 'Source' }, { key: 'f', label: 'Expected', unit: 'Hz' }, { key: 'fp', label: 'Found', unit: 'Hz' }, { key: 'p', label: 'Peak / median level', unit: '-' }],
          rows: r.psd.peaks.filter((p) => Number.isFinite(p.fPeak)).map((p) => ({ n: p.name, f: round(p.f), fp: round(p.fPeak), p: round(p.prominence) })),
          explain: 'Spectral peaks near the rotor, structural and wave frequencies with their prominence. The model has no 1P/3P excitation, so energy near those lines comes from turbulence and resonance only.' }] : []),
      ],
      notes: [`${tb.name}; ${parked ? 'rotor parked and feathered' : `production at ${fmt(inp.Vhub)} m/s (steady thrust ${fmt(op.thrust / 1e3)} kN)`}; ${inp.windModel.toUpperCase()} wind${sigma1 > 0 ? `, σ1 = ${fmt(sigma1)} m/s` : ''}${vGust ? `, gust ${fmt(vGust)} m/s` : ''}.`,
        `One realisation of ${inp.duration} s: extremes scatter between seeds by 5–15 %.`, `${gov.name} governs (utilisation ${fmt(gov.uc)}).`,
        'Rotor loads come from generic polars: use them for preliminary foundation loading and replace them with the turbine supplier\'s load documents.'],
      advice: [
        ok ? adv('engineering', sev(gov.uc), `The simulated extreme mudline moment is ${fmt(st.mudlineMoment.absMax / 1e6)} MN·m (factored stress ${fmt(sigMud / 1e6)} MPa, ${fmt(sigMud / fyd * 100)} % of design yield); rotor speed peaks at ${fmt(overspeed * 100)} % of rated.`,
          'Run the other design situations (extreme turbulence near rated, emergency stop, gust with grid loss, parked in the 50-year storm) and at least six seeds each before taking the envelope into the foundation checks; the rotor model uses generic polars and rigid blades, so the final loads must come from the turbine supplier\'s load documents.',
          'IEC 61400-3-1:2019 §7.4 (Table 2); IEC 61400-1:2019 §7.5 (six 10-minute realisations)',
          `At this load level a monopile wall of about ${fmt(Math.max(tReq, 15))} mm would reach design yield in bending above the mudline, ${fmt(dT)} t against the present wall — ${steelCo2Text(dT)} Fatigue and buckling normally govern before yield.`)
          : adv('engineering', 'fail', `${gov.name}: ${fmt(gov.demand)} against ${fmt(gov.capacity)} ${gov.unit}.`,
            gov.name.includes('overspeed') ? 'Review the pitch rate and the emergency-stop logic with the turbine supplier; a faster pitch-to-feather or earlier trip limits overspeed.' : `Increase the wall thickness at the governing section to about ${fmt((gov.name.includes('Tower') ? inp.towerTb * sigTB : inp.t * sigMud) / fyd * 1.02)} mm or the diameter, or use a higher steel grade.`,
            'IEC 61400-3-1:2019 §7.6; DNV-ST-0126 (2021) §4.4',
            `Meeting the stress check by wall thickness changes the monopile above the mudline by ${fmt(dT)} t — ${steelCo2Text(dT)}`),
        adv('economic', 'info', `Mean power in this condition is ${fmt(st.power.mean / 1e6)} MW; the damage-equivalent mudline moment is ${fmt(r.del.mudlineMoment / 1e6)} MN·m at 1 Hz.`,
          'Use the reduced-order model to rank design situations and to size the foundation concept, then spend the detailed load-simulation budget (supplier aeroelastic code, six seeds per bin) on the governing cases only.',
          'Front-end loading practice (AACE RP 18R-97 Class 4 → 3); IEC 61400-3-1 integrated load analysis as the contractual load basis'),
      ],
    };
  },
  verify: [
    { name: 'Steady wind, calm sea, no tower drag: mean mudline moment / mean thrust = hub height + water depth', source: 'Moment equilibrium of the support structure about the mudline (statics): 90 m + 25 m = 115 m', inputs: { windModel: 'steady', Vhub: 9, Hs: 0, Uc: 0, towerDrag: false, duration: 60 }, expect: [{ kpi: 'leverArm', value: 115, tol: 0.002 }] },
    { name: 'Parked in light wind with a steady current only: mean mudline shear = ½ρ·Cd·D·U²·d', source: 'Morison drag on a vertical cylinder in a uniform current (DNV-RP-C205 §6.2): 0.5 × 1025 × 0.9 × 6 × 1.0² × 25 = 69.2 kN; rotor and tower wind loads are below 0.2 kN at 2 m/s',
      inputs: { windModel: 'steady', Vhub: 2, Hs: 0, Uc: 1, Cd: 0.9, D: 6, depth: 25, towerDrag: false, duration: 60 }, expect: [{ kpi: 'mudShearMean', value: 0.0692, tol: 0.02 }] },
  ],
};
