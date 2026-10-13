// Coupled aero-hydro-servo-elastic dynamic analysis: reduced-order time-domain simulation of a bottom-fixed
// (monopile) or floating offshore wind turbine in turbulent wind and irregular waves with the baseline torque and
// pitch controller — rotor, drivetrain, tower and platform response, section loads and damage-equivalent loads.
import { simulateCoupled, buildSupportStructure } from '../../engine/aero/coupled.js';
import { surfaceLookup } from '../../engine/aero/bem.js';
import { ntmSigma, turbulenceScale, ewm, eogGust, eogSpeed, rotorEffectiveTurbulence } from '../../engine/aero/turbulence.js';
import { turbineInputs, getTurbine, tunedController, steadyAt, supportInputs, supportSegments, foundationInputs, foundationOf, adv, check, worst, sev, fmt, round, run, need, thin, linspace, G, RPM, DEG,
  LIMIT_POLAR, LIMIT_COUPLED, LIMIT_CODE, LIMIT_CONTROL, CFD_OUT, CONFIRM, SUPPLIER, steelCo2Text, segmentsMass } from './_common.js';

const fixed = { key: 'concept', in: ['monopile'] }, floating = { key: 'concept', in: ['floating'] };
const on = (list, showIf) => list.map((i) => (i.showIf ? i : { ...i, showIf }));

export default {
  id: 's16-coupled-dynamics', suite: 16, title: 'Coupled aero-hydro-servo-elastic dynamics',
  summary: 'Time-domain response of the turbine, controller, support structure and (for floaters) platform to turbulent wind and irregular waves: rotor speed, pitch, power, thrust, tower and mudline loads and their damage-equivalent values.',
  covers: [1], fidelity: 'reduced-order',
  standards: ['IEC 61400-1:2019 §6.3 (normal turbulence and extreme operating gust), §7.4 (design load cases)', 'IEC 61400-3-1:2019 (fixed offshore wind turbines)', 'IEC TS 61400-3-2:2019 (floating offshore wind turbines)', 'DNV-ST-0437 (2016) Loads and site conditions for wind turbines', 'NREL/TP-500-38060 (5-MW reference turbine and baseline controller)'],
  models: ['Fully coupled aero-hydro-servo-elastic simulation (reduced order)', 'Rotor Cp/Ct/Cq surfaces from blade element momentum theory', 'Dynamic inflow (Øye time constant)', 'Variable-speed generator torque control and gain-scheduled PI pitch control', 'Supervisory start-up / shutdown / emergency-stop logic', 'Euler–Bernoulli beam finite elements and modal reduction of the support structure', 'Coupled lateral / rocking foundation springs', 'Three-degree-of-freedom floating-body dynamics (surge, heave, pitch)', 'Kaimal turbulence with rotor averaging', 'JONSWAP irregular waves, Airy kinematics, Morison equation, MacCamy–Fuchs diffraction', 'Rainflow counting and damage-equivalent loads', 'Fourth-order Runge–Kutta integration'],
  limits: `${LIMIT_COUPLED} ${LIMIT_POLAR} Floaters use constant added mass and damping with a long-wave Froude–Krylov estimate of the wave forces (import RAO-based forces from a diffraction analysis for design); side-side and torsional response are not computed. ${CFD_OUT} A 10-minute run with at least six seeds per wind speed is required by IEC 61400-1 for design loads; the default run is shorter. ${LIMIT_CONTROL} ${LIMIT_CODE}`,
  publishes: ['owt.delMudline', 'owt.towerTopAccMax', 'owt.mudlineMomentOp', 'vessel.pitchMaxDeg'],
  searchTerms: 'aero-hydro-servo-elastic offshore wind turbine simulation monopile floating',
  inputs: [
    ...turbineInputs(),
    { key: 'concept', label: 'Support concept', type: 'select', default: 'monopile', group: 'Support structure', role: 'geometry', options: [{ value: 'monopile', label: 'Bottom-fixed: monopile and tower (flexible)' }, { value: 'floating', label: 'Floating platform (surge, heave, pitch; rigid tower)' }] },
    ...on(supportInputs({ steel: false }).filter((i) => i.key !== 'depth'), fixed),
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 30, min: 5, max: 1500, group: 'Site', role: 'environment', bind: 'site.depth', help: 'Monopile concept: 5–70 m. Floating concept: any depth (wave forces use it).' },
    { key: 'zeta', label: 'Structural damping ratio', unit: '-', type: 'number', default: 0.01, min: 0.002, max: 0.05, group: 'Support structure', role: 'material', showIf: fixed, uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Steel, soil and hydrodynamic damping lumped per mode (about 1 % of critical in total for monopiles; aerodynamic damping comes from the model).' },
    ...on(foundationInputs(), fixed),
    { key: 'pMass', label: 'Platform mass including tower, RNA and ballast', unit: 't', type: 'number', default: 14070, min: 500, max: 80000, group: 'Floating platform', role: 'geometry', showIf: floating, bind: 'vessel.disp', help: 'Semi-submersible of the OC4 DeepCwind size by default.' },
    { key: 'pInertia', label: 'Pitch inertia about the waterline', unit: '10⁶ t·m²', type: 'number', default: 12, min: 0.1, max: 500, group: 'Floating platform', role: 'geometry', showIf: floating },
    { key: 'aSurge', label: 'Surge added mass', unit: 't', type: 'number', default: 8500, min: 0, max: 1e5, group: 'Floating platform', role: 'bc', showIf: floating },
    { key: 'aHeave', label: 'Heave added mass', unit: 't', type: 'number', default: 14700, min: 0, max: 2e5, group: 'Floating platform', role: 'bc', showIf: floating },
    { key: 'aPitch', label: 'Pitch added inertia', unit: '10⁶ t·m²', type: 'number', default: 7.2, min: 0, max: 500, group: 'Floating platform', role: 'bc', showIf: floating },
    { key: 'kSurge', label: 'Surge stiffness (mooring)', unit: 'kN/m', type: 'number', default: 71, min: 1, max: 50000, group: 'Floating platform', role: 'bc', showIf: floating, bind: 'mooring.Kxx', help: 'Station-keeping restraint from the mooring module.' },
    { key: 'kHeave', label: 'Heave stiffness (hydrostatic)', unit: 'kN/m', type: 'number', default: 3836, min: 10, max: 1e6, group: 'Floating platform', role: 'bc', showIf: floating },
    { key: 'kPitch', label: 'Pitch stiffness (hydrostatic and mooring)', unit: 'MN·m/rad', type: 'number', default: 1300, min: 10, max: 1e6, group: 'Floating platform', role: 'bc', showIf: floating, help: 'ρg∇·GM plus the mooring contribution.' },
    { key: 'pDamp', label: 'Linear damping ratio of the platform modes', unit: '-', type: 'number', default: 0.05, min: 0, max: 0.5, group: 'Floating platform', role: 'bc', showIf: floating, uncertain: { dist: 'lognormal', cov: 0.4 } },
    { key: 'pDraft', label: 'Platform draft', unit: 'm', type: 'number', default: 20, min: 3, max: 150, group: 'Floating platform', role: 'geometry', showIf: floating },
    { key: 'omegaN', label: 'Pitch-controller natural frequency', unit: 'rad/s', type: 'number', default: 0.2, min: 0.05, max: 1, group: 'Control system', role: 'control', showIf: floating, help: 'Detuned below the platform pitch frequency to avoid negative damping (Larsen & Hanson 2007; 0.6 rad/s is the bottom-fixed setting).' },
    { key: 'windModel', label: 'Wind condition', type: 'select', default: 'ntm', group: 'Wind', role: 'bc', options: [{ value: 'ntm', label: 'Normal turbulence model (rotor-effective Kaimal)' }, { value: 'steady', label: 'Steady wind' }, { value: 'eog', label: 'Extreme operating gust on steady wind' }], help: 'Inflow boundary condition at hub height.' },
    { key: 'Vhub', label: 'Mean hub-height wind speed', unit: 'm/s', type: 'number', default: 12, min: 3, max: 30, group: 'Wind', role: 'bc', design: true, uncertain: { dist: 'normal', cov: 0.05 } },
    { key: 'Iref', label: 'Reference turbulence intensity Iref', unit: '-', type: 'number', default: 0.12, min: 0.05, max: 0.2, group: 'Wind', role: 'bc', help: `IEC 61400-1 Table 1: A 0.16, B 0.14, C 0.12 (offshore sites are usually C or below). ${CONFIRM}` },
    { key: 'Vref', label: 'Reference wind speed Vref (turbine class)', unit: 'm/s', type: 'number', default: 50, min: 30, max: 60, group: 'Wind', role: 'bc', showIf: { key: 'windModel', in: ['eog'] }, help: 'Class I 50 m/s, II 42.5 m/s, III 37.5 m/s; sets the gust magnitude.' },
    { key: 'windSeed', label: 'Turbulence seed', type: 'int', default: 3, min: 1, max: 99999, group: 'Wind', role: 'ic' },
    { key: 'yawErr', label: 'Yaw misalignment', unit: 'deg', type: 'number', default: 0, min: -30, max: 30, group: 'Wind', role: 'bc' },
    { key: 'Hs', label: 'Significant wave height', unit: 'm', type: 'number', default: 2, min: 0, max: 14, group: 'Waves and current', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Sea state associated with the wind speed (normal sea state). 0 = still water.' },
    { key: 'Tp', label: 'Peak period', unit: 's', type: 'number', default: 7, min: 3, max: 20, group: 'Waves and current', role: 'bc' },
    { key: 'gamma', label: 'JONSWAP peak-enhancement factor', unit: '-', type: 'number', default: 3.3, min: 1, max: 7, group: 'Waves and current', role: 'bc' },
    { key: 'Uc', label: 'Current speed (uniform)', unit: 'm/s', type: 'number', default: 0.3, min: 0, max: 3, group: 'Waves and current', role: 'bc', showIf: fixed },
    { key: 'macCamy', label: 'MacCamy–Fuchs diffraction correction of the inertia force', type: 'bool', default: true, group: 'Waves and current', role: 'bc', showIf: fixed },
    { key: 'waveSeed', label: 'Wave seed', type: 'int', default: 1, min: 1, max: 99999, group: 'Waves and current', role: 'ic' },
    { key: 'initState', label: 'Initial turbine state', type: 'select', default: 'steady', group: 'Initial state', role: 'ic', options: [{ value: 'steady', label: 'Steady operating point of the mean wind speed' }, { value: 'parked', label: 'Parked, blades feathered (start-up from standstill)' }],
      help: 'Rotor speed, pitch and controller mode at t = 0. Structural deflections start from zero and settle during the discarded transient.' },
    { key: 'event', label: 'Supervisory event', type: 'select', default: 'none', group: 'Control system', role: 'control', options: [{ value: 'none', label: 'None (continuous production)' }, { value: 'stop', label: 'Normal shutdown' }, { value: 'estop', label: 'Emergency stop' }, { value: 'gridloss', label: 'Grid loss with emergency stop' }] },
    { key: 'tEvent', label: 'Event time', unit: 's', type: 'number', default: 150, min: 1, max: 3500, group: 'Control system', role: 'control', showIf: { key: 'event', in: ['stop', 'estop', 'gridloss'] } },
    { key: 'pitchRate', label: 'Pitch-rate limit', unit: 'deg/s', type: 'number', default: 8, min: 1, max: 20, group: 'Control system', role: 'control', help: 'NREL 5-MW baseline: 8 deg/s.' },
    { key: 'dynInflow', label: 'Dynamic inflow (Øye lag of the induced velocity)', type: 'bool', default: false, group: 'Control system', role: 'control', help: 'Off = quasi-steady rotor aerodynamics. The lag is integrated to first order in the time step: halve the step when it is switched on and check the time-step study.' },
    { key: 'overspeed', label: 'Over-speed limit (fraction of rated rotor speed)', unit: '-', type: 'number', default: 1.2, min: 1.05, max: 1.5, group: 'Acceptance', role: 'criteria', help: `Protection trip level; the rotor must stay below it in every load case. ${CONFIRM}` },
    { key: 'accLimit', label: 'Allowable nacelle fore-aft acceleration', unit: 'm/s²', type: 'number', default: 3, min: 0.5, max: 8, group: 'Acceptance', role: 'criteria', help: 'Supplier limit for nacelle components (commonly 0.2–0.4 g).' },
    { key: 'pitchLimit', label: 'Allowable platform pitch', unit: 'deg', type: 'number', default: 10, min: 2, max: 20, group: 'Acceptance', role: 'criteria', showIf: floating },
    { key: 'wohler', label: 'Wöhler exponent for damage-equivalent loads', unit: '-', type: 'number', default: 4, min: 3, max: 12, group: 'Acceptance', role: 'criteria' },
    { key: 'duration', label: 'Simulated duration', unit: 's', type: 'number', default: 300, min: 60, max: 3600, group: 'Solver', role: 'control', help: 'Includes a 60 s start-up transient that is discarded. Use 660 s for a 10-minute design run.' },
    { key: 'dt', label: 'Time step', unit: 's', type: 'number', default: 0.05, min: 0.005, max: 0.1, group: 'Solver', role: 'control' },
    { key: 'nElem', label: 'Beam elements of the support structure', type: 'int', default: 40, min: 6, max: 160, group: 'Solver', role: 'control', showIf: fixed },
  ],
  mesh: { key: 'nElem', label: 'Beam elements', levels: [10, 20, 40, 80], qoi: 'f1' },
  meshTime: { key: 'dt', label: 'Time step (s)', levels: [0.1, 0.05, 0.025], qoi: 'mudMomStd', isSize: true },
  compute(inp) {
    const tb = getTurbine(inp), rotor = tb.rotor, R = rotor.tipRadius, fix = inp.concept === 'monopile';
    need(inp.hubHeight > R + 5, `The hub height (${inp.hubHeight} m) leaves no clearance under the blade tip (rotor radius ${fmt(R)} m).`);
    need(inp.duration / inp.dt <= 2.4e5, 'Too many time steps: increase the time step or shorten the run.');
    const tTrans = Math.min(60, inp.duration / 3);
    // support structure or floater
    let structure = null, floater = null, segs = null;
    if (fix) {
      need(inp.depth <= 70, 'The monopile concept is limited to 70 m water depth in this module: select the floating concept.');
      segs = supportSegments(inp, { zHub: inp.hubHeight, hubOffset: 2.4 * tb.scale });
      structure = run(() => buildSupportStructure({ segments: segs, rnaMass: rotor.rnaMass, zHub: inp.hubHeight, zTowerBase: inp.transZ, foundation: foundationOf(inp), nElem: inp.nElem, nModes: 2, zeta: inp.zeta }), 'Support-structure model');
    } else {
      const M = inp.pMass * 1e3, I = inp.pInertia * 1e9, A = [inp.aSurge * 1e3, inp.aHeave * 1e3, inp.aPitch * 1e9], K = [inp.kSurge * 1e3, inp.kHeave * 1e3, inp.kPitch * 1e6];
      const Mt = [M, M, I];
      const damping = Mt.map((m, i) => 2 * inp.pDamp * Math.sqrt(K[i] * (m + A[i])));
      // long-wave Froude–Krylov estimate of the first-order wave forces per unit wave amplitude (deep-water decay over the draft)
      const f = linspace(0.03, 0.4, 38), vol = M / 1025, Aw = K[1] / (1025 * G), zc = inp.pDraft / 2;
      const surge = f.map((fi) => { const w = 2 * Math.PI * fi, k = w * w / G; return (1025 * vol + A[0]) * w * w * Math.exp(-k * zc); });
      const heave = f.map((fi) => { const w = 2 * Math.PI * fi, k = w * w / G; return Math.max(0, 1025 * G * Aw - w * w * A[1]) * Math.exp(-k * inp.pDraft); });
      floater = { mass: Mt, addedMass: A, stiffness: K, damping, quadDamping: [0.5 * 1025 * 1.0 * inp.pDraft * 50, 0, 0], refZ: 0, towerBaseZ: 10, rnaMass: rotor.rnaMass,
        forceRao: { f, surge: { amp: surge, phaseDeg: f.map(() => 90) }, heave: { amp: heave, phaseDeg: f.map(() => 0) }, pitch: { amp: surge.map((s) => s * zc), phaseDeg: f.map(() => -90) } } };
    }
    const controller = { ...(fix ? tb.controller : tunedController({ rotor, surface: tb.surface, schedule: tb.schedule, omegaN: inp.omegaN })), pitchRate: inp.pitchRate * DEG, overspeedFactor: inp.overspeed + 0.1 };
    // wind series
    const sigma1 = ntmSigma({ Vhub: inp.Vhub, Iref: inp.Iref });
    let wind, gust = 0;
    const tGust = Math.max(tTrans + 30, 0.5 * inp.duration);
    if (inp.windModel === 'ntm') {
      // synthesised on a fixed 0.1 s grid so that the realisation does not change with the solver time step
      // (a power-of-two number of samples keeps the spectral synthesis on the fast transform)
      const dtw = 0.1, nw = 2 ** Math.ceil(Math.log2((inp.duration + 2) / dtw)), ws = run(() => rotorEffectiveTurbulence({ Vhub: inp.Vhub, sigma1, hubHeight: inp.hubHeight, R, duration: nw * dtw, dt: dtw, seed: inp.windSeed }), 'Turbulent wind').u;
      wind = (t) => { const x = Math.max(0, t) / dtw, i = Math.min(ws.length - 2, Math.floor(x)); return ws[i] + (ws[i + 1] - ws[i]) * (x - i); };
    }
    else if (inp.windModel === 'eog') {
      const Ve1 = ewm({ Vref: inp.Vref, z: inp.hubHeight, hubHeight: inp.hubHeight }).Ve1;
      gust = eogGust({ Vhub: inp.Vhub, sigma1, D: 2 * R, lambda1: turbulenceScale(inp.hubHeight), Ve1 });
      wind = (t) => eogSpeed(t - tGust, { Vhub: inp.Vhub, Vgust: gust });
    } else wind = inp.Vhub;
    const st0 = steadyAt(tb, inp.Vhub);
    const parked = inp.initState === 'parked' || st0.omega === 0;
    const initial = parked ? { omega: 0.02, pitchDeg: 90, mode: 'startup' } : { omega: st0.omega, pitchDeg: st0.pitchDeg, mode: 'production' };
    const events = inp.event === 'none' ? undefined : (t) => (t >= inp.tEvent && t < inp.tEvent + 2 * inp.dt ? (inp.event === 'stop' ? { command: 'stop' } : { command: 'estop', gridLoss: inp.event === 'gridloss' }) : (inp.event === 'gridloss' && t >= inp.tEvent ? { gridLoss: true } : {}));
    const sim = run(() => simulateCoupled({ surface: tb.surface, turbine: { inertia: rotor.drivetrainInertia, hubHeight: inp.hubHeight }, controller, ...(fix ? { structure } : { floater }), wind,
      ...(inp.Hs > 0 ? { waves: { Hs: inp.Hs, Tp: inp.Tp, gamma: inp.gamma, seed: inp.waveSeed, depth: inp.depth, macCamyFuchs: fix && inp.macCamy, current: fix ? inp.Uc : 0 } } : {}),
      dt: inp.dt, duration: inp.duration, tTransient: tTrans, yawErrorDeg: inp.yawErr, shearExp: 0.14, towerWindDrag: fix, dynamicInflow: inp.dynInflow, initial, events, del: { m: inp.wohler } }), 'Coupled simulation');
    const s = sim.stats;
    for (const k of ['thrust', 'mudlineMoment', 'rotorSpeed']) need(Number.isFinite(s[k].mean) && Number.isFinite(s[k].absMax), 'The simulation diverged (non-finite response): reduce the time step, or check the stiffness and damping inputs.');
    const f1 = fix ? structure.modes[0].frequency : 0, f2 = fix ? structure.modes[1].frequency : 0;
    // aerodynamic damping of the first fore-aft mode from the thrust slope at frozen rotor speed and pitch (Salzmann & van der Tempel 2005)
    let zetaAero = 0;
    if (fix && st0.omega > 0) {
      const lk = surfaceLookup(tb.surface, st0.omega * tb.surface.R / inp.Vhub, st0.pitchDeg), A = tb.surface.area;
      const dTdV = 1.225 * A * inp.Vhub * lk.ct - 0.5 * 1.225 * A * st0.omega * tb.surface.R * lk.dCtDTsr;
      const m = structure.modes[0];
      zetaAero = Math.max(0, m.phiHub ** 2 * dTdV / (2 * Math.sqrt(m.stiffness * m.mass)));
    }
    const rpmMax = s.rotorSpeed.max / RPM, rpmRated = rotor.ratedRotorSpeed / RPM;
    const pitchMaxDeg = floater ? Math.max(Math.abs(s['platform.pitch'].max), Math.abs(s['platform.pitch'].min)) / DEG : 0;
    const checks = [
      check('Maximum rotor speed against the over-speed limit', rpmMax, inp.overspeed * rpmRated, 'rpm', 'IEC 61400-1:2019 §8 (control and protection functions)', 'Highest rotor speed in the run against the protection trip level. If it is reached the safety system shuts the turbine down; repeated trips mean the pitch controller is too slow for the gusts.'),
      check('Nacelle fore-aft acceleration', s.towerTopAcc.absMax, inp.accLimit, 'm/s²', 'Supplier nacelle acceleration limit; DNV-ST-0437 §4 (load effects)', 'Peak tower-top acceleration against the allowable for nacelle equipment. Exceeding it calls for more damping (tower damper, controller tuning) or a stiffer support.'),
      ...(floater ? [check('Maximum platform pitch', pitchMaxDeg, inp.pitchLimit, 'deg', 'IEC TS 61400-3-2:2019; DNV-ST-0119 (2021) floating wind turbine structures (operational tilt limits)', 'Largest platform pitch against the operating limit. Large pitch loses energy (cos³ of the tilt), raises tower-base bending from the RNA weight and can trip the turbine.')] : []),
    ];
    const gov = worst(checks);
    const steel = fix ? segmentsMass(segs) / 1e3 : 0;
    const th = (arrs) => thin(arrs, 700);
    const [t, w, eta, rs, pd, pw, thr, ttd, mm, tbm] = th([sim.t, sim.wind, sim.waveElevation, sim.rotorSpeed, sim.pitchDeg, sim.power, sim.thrust, sim.towerTopDisp, sim.mudlineMoment, sim.towerBaseMoment]);
    const T = t.map((v) => round(v, 3));
    const psd = sim.psd;
    const ip = psd ? thin([Array.from(psd.f), Array.from(psd.towerBaseMoment)], 300) : null;
    const chan = [['Wind speed', 'wind', 1, 'm/s'], ['Wave elevation', 'waveElevation', 1, 'm'], ['Rotor speed', 'rotorSpeed', 1 / RPM, 'rpm'], ['Blade pitch', 'pitchDeg', 1, 'deg'], ['Electrical power', 'power', 1e-6, 'MW'], ['Rotor thrust', 'thrust', 1e-3, 'kN'],
      ['Tower-top displacement', 'towerTopDisp', 1, 'm'], ['Tower-top acceleration', 'towerTopAcc', 1, 'm/s²'], ['Tower-base moment', 'towerBaseMoment', 1e-6, 'MN·m'], ['Tower-base shear', 'towerBaseShear', 1e-3, 'kN'], ['Mudline moment', 'mudlineMoment', 1e-6, 'MN·m'], ['Mudline shear', 'mudlineShear', 1e-3, 'kN'],
      ...(floater ? [['Platform surge', 'platform.surge', 1, 'm'], ['Platform heave', 'platform.heave', 1, 'm'], ['Platform pitch', 'platform.pitch', 1 / DEG, 'deg']] : [])];
    const stopped = sim.mode[sim.mode.length - 1];
    const lever = fix ? inp.hubHeight + inp.depth : inp.hubHeight - 10;
    return {
      kpis: [
        { key: 'mudMomMax', label: fix ? 'Maximum mudline bending moment' : 'Maximum tower-base bending moment', value: round(s.mudlineMoment.absMax / 1e6, 5), unit: 'MN·m', explain: 'Largest fore-aft overturning moment after the start-up transient (characteristic, no load factor). It combines thrust, wave load and the dynamic response of the first mode.' },
        { key: 'mudMomMean', label: 'Mean bending moment at the same section', value: round(s.mudlineMoment.mean / 1e6, 5), unit: 'MN·m', explain: 'Time-average moment, essentially mean thrust times its lever arm. The difference to the maximum is the dynamic part.' },
        { key: 'mudMomStd', label: 'Standard deviation of that moment', value: round(s.mudlineMoment.std / 1e6, 5), unit: 'MN·m', explain: 'Fluctuation of the overturning moment: the fatigue driver. It grows with turbulence, wave height and proximity of the first natural frequency to the wave or rotor frequencies.' },
        { key: 'delMud', label: `Damage-equivalent moment (m = ${inp.wohler}, 1 Hz)`, value: round(sim.del.mudlineMoment / 1e6, 5), unit: 'MN·m', explain: 'Constant moment range that, applied once per second, causes the same fatigue damage as the simulated history for the chosen Wöhler exponent. It is the input to the fatigue module.' },
        { key: 'thrustMean', label: 'Mean rotor thrust', value: round(s.thrust.mean / 1e3, 5), unit: 'kN', explain: 'Average aerodynamic thrust in the run; it should match the steady thrust curve at the mean wind speed.' },
        { key: 'thrustMax', label: 'Maximum rotor thrust', value: round(s.thrust.max / 1e3, 5), unit: 'kN', explain: 'Peak thrust, typically during a gust just below rated wind speed before the pitch controller reacts.' },
        { key: 'powerMean', label: 'Mean electrical power', value: round(s.power.mean / 1e6, 5), unit: 'MW', explain: 'Average generator output over the analysed window.' },
        { key: 'rpmMax', label: 'Maximum rotor speed', value: round(rpmMax, 4), unit: 'rpm', limit: round(inp.overspeed * rpmRated, 4), sense: 'max', explain: 'Highest rotor speed; the margin to the over-speed trip shows how well the pitch loop rejects gusts.' },
        { key: 'ttDispMax', label: 'Maximum tower-top displacement', value: round(s.towerTopDisp.absMax, 4), unit: 'm', explain: 'Peak fore-aft deflection of the tower top (elastic part; for floaters it includes platform motion).' },
        { key: 'ttAccMax', label: 'Maximum nacelle fore-aft acceleration', value: round(s.towerTopAcc.absMax, 4), unit: 'm/s²', limit: inp.accLimit, sense: 'max', explain: 'Peak acceleration at the tower top, compared with the limit for nacelle components.' },
        { key: 'f1', label: 'First fore-aft natural frequency', value: round(f1, 5), unit: 'Hz', explain: 'First bending frequency of the monopile–tower–RNA system with the foundation springs and hydrodynamic added mass (0 for the rigid-tower floating model). It should lie between the 1P and 3P bands.' },
        { key: 'zetaAero', label: 'Aerodynamic damping of the first mode', value: round(zetaAero * 100, 4), unit: '% of critical', explain: 'Damping from the thrust slope at the operating point (frozen rotor speed and pitch). Typically 2–6 % in production and near zero when parked, which is why idling fatigue matters for monopiles.' },
        { key: 'momentArm', label: 'Mean moment / mean thrust', value: round(s.thrust.mean > 1 ? s.mudlineMoment.mean / s.thrust.mean : 0, 5), unit: 'm', explain: 'Effective lever arm of the mean load. In steady wind without waves or tower drag it must equal the height of the hub above the section — a load-transfer consistency check.' },
        { key: 'pitchMax', label: 'Maximum platform pitch', value: round(pitchMaxDeg, 4), unit: 'deg', explain: 'Largest platform pitch angle (0 for the bottom-fixed concept).' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Wind speed and wave elevation', x: { label: 'Time', unit: 's' }, y: { label: 'Wind (m/s) · wave (m)', unit: '-' }, series: [{ name: 'Rotor-effective wind speed (m/s)', x: T, y: w.map((v) => round(v, 4)) }, { name: 'Wave elevation (m)', x: T, y: eta.map((v) => round(v, 4)) }],
          explain: 'The two environmental inputs of the run. Compare the gusts and wave groups with the response plots below to see what drives each load peak.' },
        { type: 'line', title: 'Rotor speed, pitch and power', x: { label: 'Time', unit: 's' }, y: { label: 'rpm · deg · MW', unit: '-' }, series: [{ name: 'Rotor speed (rpm)', x: T, y: rs.map((v) => round(v / RPM, 4)) }, { name: 'Blade pitch (deg)', x: T, y: pd.map((v) => round(v, 4)) }, { name: 'Electrical power (MW)', x: T, y: pw.map((v) => round(v / 1e6, 4)) }],
          explain: 'Controller response. Below rated the pitch stays at fine and rotor speed follows the wind; above rated the pitch moves to hold speed and power. After a stop event the pitch runs to feather and power falls to zero.' },
        { type: 'line', title: 'Rotor thrust and tower-top displacement', x: { label: 'Time', unit: 's' }, y: { label: 'Thrust (MN) · displacement (m)', unit: '-' }, series: [{ name: 'Thrust (MN)', x: T, y: thr.map((v) => round(v / 1e6, 4)) }, { name: 'Tower-top displacement (m)', x: T, y: ttd.map((v) => round(v, 4)) }],
          explain: 'Thrust and the resulting fore-aft motion of the tower top. Oscillation at the first natural frequency superposed on the slow thrust variation shows the dynamic amplification.' },
        { type: 'line', title: fix ? 'Mudline and tower-base bending moment' : 'Tower-base bending moment', x: { label: 'Time', unit: 's' }, y: { label: 'Moment', unit: 'MN·m' }, series: [{ name: fix ? 'Mudline' : 'Tower base', x: T, y: mm.map((v) => round(v / 1e6, 4)) }, ...(fix ? [{ name: 'Tower base', x: T, y: tbm.map((v) => round(v / 1e6, 4)) }] : [])],
          explain: 'Section bending moments. The mudline moment adds the wave load on the monopile to the rotor loads; its cycles are counted for the damage-equivalent load.' },
        ...(ip ? [{ type: 'line', title: 'Spectrum of the tower-base moment', x: { label: 'Frequency', unit: 'Hz' }, y: { label: 'log₁₀ PSD', unit: '(N·m)²/Hz' }, series: [{ name: 'Tower-base moment', x: ip[0].slice(1).map((v) => round(v, 4)), y: ip[1].slice(1).map((v) => round(Math.log10(Math.max(v, 1)), 4)) }, ...(fix ? [{ name: 'First mode', x: [f1, f1], y: [Math.log10(Math.max(Math.min(...ip[1].slice(1)), 1)), Math.log10(Math.max(...ip[1]))], dash: true }] : [])],
          explain: 'Power spectral density of the tower-base moment (log scale). Peaks at the wave frequency, at the first natural frequency and, for floaters, at the platform periods show where the fatigue energy sits.' }] : []),
        ...(floater ? [{ type: 'line', title: 'Platform motions', x: { label: 'Time', unit: 's' }, y: { label: 'Surge, heave (m) · pitch (deg)', unit: '-' }, series: (() => { const [a, b, c] = th([sim.platform.surge, sim.platform.heave, sim.platform.pitch]); return [{ name: 'Surge (m)', x: T, y: a.map((v) => round(v, 4)) }, { name: 'Heave (m)', x: T, y: b.map((v) => round(v, 4)) }, { name: 'Pitch (deg)', x: T, y: c.map((v) => round(v / DEG, 4)) }]; })(),
          explain: 'Surge, heave and pitch of the platform. A slowly growing pitch oscillation above rated wind speed is the signature of controller-induced negative damping.' }] : []),
        ...(fix ? [{ type: 'line', title: 'Fore-aft mode shapes', x: { label: 'Normalised lateral displacement', unit: '-' }, y: { label: 'Elevation', unit: 'm' }, series: structure.modes.map((m, i) => ({ name: `Mode ${i + 1} (${fmt(m.frequency)} Hz)`, x: Array.from(m.shape, (v) => round(v, 4)), y: Array.from(structure.nodes, (v) => round(v, 3)) })),
          explain: 'First two bending mode shapes of the support structure, scaled to unit hub displacement. Displacement at the mudline shows the flexibility of the foundation.' }] : []),
      ],
      tables: [
        { title: 'Response statistics', columns: [{ key: 'ch', label: 'Channel' }, { key: 'mean', label: 'Mean' }, { key: 'std', label: 'Std. dev.' }, { key: 'min', label: 'Min' }, { key: 'max', label: 'Max' }, { key: 'u', label: 'Unit' }],
          rows: chan.filter(([, k]) => s[k]).map(([ch, k, f, u]) => ({ ch, mean: round(s[k].mean * f, 5), std: round(s[k].std * f, 5), min: round(s[k].min * f, 5), max: round(s[k].max * f, 5), u })),
          explain: 'Mean, standard deviation and extremes of every channel after the start-up transient. Characteristic values: apply the IEC partial load factors (1.35 normal, 1.1 abnormal) for design.' },
        ...(psd ? [{ title: 'Spectral peaks of the tower-base moment', columns: [{ key: 'name', label: 'Source' }, { key: 'f', label: 'Expected', unit: 'Hz' }, { key: 'fp', label: 'Peak found', unit: 'Hz' }, { key: 'prom', label: 'Prominence (× median)', unit: '-' }],
          rows: psd.peaks.map((p) => ({ name: p.name, f: round(p.f, 4), fp: round(p.fPeak, 4), prom: round(p.prominence, 4) })),
          explain: 'Spectral level near each expected excitation or natural frequency relative to the median level. A high prominence at a natural frequency marks resonant response.' }] : []),
      ],
      notes: [
        `${rotor.name} on a ${fix ? `monopile in ${inp.depth} m of water (first mode ${fmt(f1)} Hz, second ${fmt(f2)} Hz)` : 'floating platform (rigid tower)'}; ${inp.windModel === 'ntm' ? `turbulent wind ${inp.Vhub} m/s, σ₁ = ${fmt(sigma1)} m/s` : inp.windModel === 'eog' ? `extreme operating gust of ${fmt(gust)} m/s on ${inp.Vhub} m/s at t = ${fmt(tGust)} s` : `steady wind ${inp.Vhub} m/s`}; sea state Hs ${inp.Hs} m, Tp ${inp.Tp} s.`,
        `Analysed window ${fmt(inp.duration - tTrans, 4)} s after a ${fmt(tTrans)} s transient; controller mode at the end: ${stopped}.`,
        `Governing check: ${gov.name} (utilisation ${fmt(gov.uc)}). Mean lever arm ${fmt(s.thrust.mean > 1 ? s.mudlineMoment.mean / s.thrust.mean : 0, 4)} m against ${fmt(lever, 4)} m geometric.`,
        'Characteristic loads from one short realisation: IEC 61400-1 asks for six 10-minute seeds per wind speed and extrapolation for ultimate loads.',
      ],
      advice: [
        adv('engineering', sev(gov.uc), `Peak ${fix ? 'mudline' : 'tower-base'} moment ${fmt(s.mudlineMoment.absMax / 1e6)} MN·m (mean ${fmt(s.mudlineMoment.mean / 1e6)}, damage-equivalent ${fmt(sim.del.mudlineMoment / 1e6)} MN·m); rotor speed peaks at ${fmt(rpmMax)} rpm and nacelle acceleration at ${fmt(s.towerTopAcc.absMax)} m/s².`,
          gov.uc > 1 ? `'${gov.name}' is exceeded (utilisation ${fmt(gov.uc)}): ${gov === checks[0] ? 'raise the pitch-rate limit or the controller bandwidth' : gov === checks[1] ? 'add tower damping or move the first frequency away from the wave peak' : `lower the controller natural frequency below the platform pitch frequency or add pitch stiffness`}, then re-run.` : `Response is within the limits. Run the full wind-speed range with six seeds each for design loads, and ${fix ? `keep the first frequency (${fmt(f1)} Hz) clear of the wave peak (${fmt(1 / inp.Tp)} Hz) and of 1P/3P` : 'import diffraction-based wave forces'}; ${SUPPLIER}`,
          'IEC 61400-1:2019 §7.4 and IEC 61400-3-1:2019 (design load cases, six seeds per bin); reduced-order simulation for concept design, to be confirmed by a full aeroelastic code before certification.',
          fix ? `Support structure of this run: ${fmt(steel)} t of steel. ${steelCo2Text(steel)} Lower dynamic loads (damping, frequency placement) translate directly into wall thickness and tonnes.` : 'For floaters, controller tuning that limits platform pitch reduces tower and mooring steel at no material cost — one of the cheapest carbon savings in the design.'),
        adv('economic', 'info', `The damage-equivalent moment of ${fmt(sim.del.mudlineMoment / 1e6)} MN·m and the aerodynamic damping of ${fmt(zetaAero * 100)} % are what size the fatigue-governed wall thickness${fix ? ` of a ${fmt(steel)} t support structure` : ''}.`,
          `Fatigue damage scales with the load range to the power ${inp.wohler}: a 10 % lower damage-equivalent load (tower damper, better frequency placement) allows about ${fmt((1 - 0.9) * 100)} % less section modulus, i.e. roughly ${fmt(steel * 0.1)} t of steel, or ${fmt((1 / 0.9 ** inp.wohler - 1) * 100)} % more fatigue life at the same weight.`,
          'Whole-life costing per ISO 15686-5 and load-mitigation value studies (IEA Wind Task 37): control features are weighed against steel mass and life extension, not against energy alone.'),
      ],
      publish: { 'owt.delMudline': sim.del.mudlineMoment, 'owt.towerTopAccMax': s.towerTopAcc.absMax, 'owt.mudlineMomentOp': s.mudlineMoment.absMax, 'vessel.pitchMaxDeg': pitchMaxDeg },
      quantities: { steel_t: round(steel, 4) },
    };
  },
  verify: [
    { name: 'Uniform cantilever with tip mass: f₁ = (1/2π)·√(3EI / ((M + 33/140·mL)·L³))', source: 'Rayleigh quotient with the static deflection shape; Blevins, Formulas for Natural Frequency and Mode Shape, Table 8-8',
      // D = 6 m, t = 50 mm: I = π/64·(6⁴ − 5.9⁴) = 4.1363 m⁴, A = π·5.95·0.05 = 0.93462 m², m = 7336.8 kg/m; L = 5 + 87.6 = 92.6 m → mL = 679.4 t
      // k = 3·210e9·4.1363/92.6³ = 3.2819 MN/m; M_eff = 350 000 + 0.2357·679 388 = 510 130 kg → f = √(3.2819e6/510 130)/2π = 0.4037 Hz
      // (5 m of water: added mass and enclosed water act next to the clamp where the mode shape vanishes)
      inputs: { foundation: 'fixed', depth: 5, monoD: 6, monoT: 50, towerBaseD: 6, towerTopD: 6, towerBaseT: 50, towerTopT: 50, towerMassFactor: 1, hubHeight: 90, transZ: 10, Hs: 0, windModel: 'steady', Vhub: 8, duration: 90, nElem: 40 },
      expect: [{ kpi: 'f1', value: 0.4037, tol: 0.02 }] },
    { name: 'Load transfer in steady wind without waves: mean mudline moment = mean thrust × (hub height + depth)', source: 'Statics of the cantilevered support structure (equilibrium check)',
      // tower wind drag adds a small moment: lever arm within 3 % of 90 + 30 = 120 m
      inputs: { windModel: 'steady', Vhub: 10, Hs: 0, Uc: 0, yawErr: 0, duration: 200 }, expect: [{ kpi: 'momentArm', value: 120, tol: 0.03 }] },
  ],
};
