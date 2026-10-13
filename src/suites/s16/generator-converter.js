// Drivetrain, generator and power-conversion chain of a wind (or tidal) turbine: gearbox losses, permanent-magnet
// synchronous generator in the d–q frame, back-to-back converter losses, transformer, the resulting efficiency
// curve and annual conversion efficiency, and the torsional response of the drivetrain to a loss of generator torque.
import { pmsmSteadyState, converterLoss, gearboxEfficiency, optimalTorqueGain } from '../../engine/electrical/machines.js';
import { transformer } from '../../engine/electrical/circuit.js';
import { weibull } from '../../engine/aero/aep.js';
import { turbineInputs, getTurbine, optimumCp, site, adv, check, worst, sev, fmt, round, run, need, thin, RPM, LIMIT_POLAR, LIMIT_CODE, LIMIT_CONTROL, CONFIRM, energyCo2Text, CO2_COPPER } from './_common.js';

export default {
  id: 's16-generator-converter', suite: 16, title: 'Drivetrain, generator and power converter',
  summary: 'Efficiency of gearbox, permanent-magnet generator, converter and transformer over the operating range, annual conversion efficiency, electrical operating limits and the drivetrain torsional response to a sudden loss of generator torque.',
  covers: [12, 13], fidelity: 'analytical',
  standards: ['IEC 61400-1:2019 §9–10 (mechanical and electrical systems)', 'IEC 61400-4:2012 (gearboxes)', 'IEC 60034-1 / IEC 60034-2-1 (rotating machines: rating, losses and efficiency)', 'IEC 61800-4 / IEC 62477-1 (power converters)', 'IEC 60076-1 (power transformers) and IEC 60076-16 (wind-turbine transformers)', 'IEC 61400-21-1:2019 (electrical characteristics, grid connection)'],
  models: ['Electrical machine d–q axis model of the permanent-magnet synchronous generator', 'Electromagnetic torque relationship Te = 1.5·p·ψm·iq (maximum torque per ampere)', 'Faraday induction: speed-proportional back-EMF', 'Copper, iron and friction loss models', 'Averaged converter model: conduction and switching losses of a two-level IGBT bridge', 'Transformer equivalent circuit from test data', 'Gearbox stage and no-load losses', 'Generator torque control: optimal torque gain k·Ω² (maximum power point tracking)', 'Two-inertia torsional drivetrain with Runge–Kutta integration', 'Weibull-weighted annual conversion efficiency'],
  limits: `Fundamental-frequency, steady-state electrical models with averaged converter losses: harmonics, switching transients, field weakening limits of a real design, thermal behaviour and grid-fault ride-through control are not modelled. Machine parameters are derived from per-unit values and are generic; device loss data are placeholders. The torsional model has two inertias (rotor and generator) on one shaft stiffness. ${LIMIT_POLAR} Maxwell field (finite-element) machine design and switching-level converter simulation are outside the module. ${LIMIT_CONTROL} ${LIMIT_CODE}`,
  publishes: ['drivetrain.efficiencyRated', 'drivetrain.efficiencyAnnual', 'drivetrain.shaftTorqueMin'],
  searchTerms: 'wind turbine PMSG full converter efficiency drivetrain torsional grid loss',
  inputs: [
    ...turbineInputs({ hub: false }),
    { key: 'concept', label: 'Drivetrain concept', type: 'select', default: 'geared', group: 'Drivetrain', role: 'geometry', options: [{ value: 'geared', label: 'Three-stage gearbox, high-speed generator' }, { value: 'medium', label: 'Two-stage gearbox, medium-speed generator' }, { value: 'direct', label: 'Direct drive' }] },
    { key: 'ratio', label: 'Gearbox ratio', unit: '-', type: 'number', default: 97, min: 1, max: 200, group: 'Drivetrain', role: 'geometry', showIf: { key: 'concept', in: ['geared', 'medium'] }, help: 'NREL 5-MW reference: 97. Medium-speed designs use 30–40.' },
    { key: 'p0', label: 'Gearbox no-load loss / rated power', unit: '%', type: 'number', default: 0.5, min: 0, max: 3, group: 'Drivetrain', role: 'operation', showIf: { key: 'concept', in: ['geared', 'medium'] } },
    { key: 'kShaft', label: 'Shaft torsional stiffness (low-speed side)', unit: 'MN·m/rad', type: 'number', default: 867.6, min: 10, max: 1e5, group: 'Drivetrain', role: 'material', help: 'NREL 5-MW: 867.6 MN·m/rad. Scale with rating to the power 1.67 for a scaled rotor.' },
    { key: 'cShaft', label: 'Shaft torsional damping', unit: 'MN·m·s/rad', type: 'number', default: 6.2, min: 0, max: 1000, group: 'Drivetrain', role: 'material', help: 'NREL 5-MW: 6.215 MN·m·s/rad (about 5 % of critical).' },
    { key: 'Vll', label: 'Generator rated line voltage', unit: 'kV', type: 'number', default: 0.69, min: 0.4, max: 13.8, group: 'Generator', role: 'geometry', help: '690 V up to about 6 MW; 3.3 kV for larger machines.' },
    { key: 'polePairs', label: 'Pole pairs', type: 'int', default: 3, min: 1, max: 300, group: 'Generator', role: 'geometry', help: 'Sets the electrical frequency: f = p·n/60. Direct-drive machines have 80–200.' },
    { key: 'xd', label: 'Synchronous reactance', unit: 'pu', type: 'number', default: 0.4, min: 0.1, max: 1.5, group: 'Generator', role: 'material', help: 'Per unit of the rated impedance V²/S (non-salient machine: Ld = Lq).' },
    { key: 'rs', label: 'Stator resistance', unit: 'pu', type: 'number', default: 0.01, min: 0.001, max: 0.08, group: 'Generator', role: 'material', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Copper loss at rated current as a fraction of rated power: about 1 % for high-speed, 3–5 % for direct drive.' },
    { key: 'iron', label: 'Iron loss at rated speed / rated power', unit: '%', type: 'number', default: 0.8, min: 0, max: 4, group: 'Generator', role: 'material' },
    { key: 'fric', label: 'Friction and windage at rated speed / rated power', unit: '%', type: 'number', default: 0.3, min: 0, max: 3, group: 'Generator', role: 'material' },
    { key: 'emf', label: 'Back-EMF at rated speed / rated voltage', unit: '-', type: 'number', default: 0.92, min: 0.6, max: 1.1, group: 'Generator', role: 'material', help: 'Sets the permanent-magnet flux linkage (Faraday: e = p·ω·ψm).' },
    { key: 'Vdc', label: 'DC-link voltage', unit: 'kV', type: 'number', default: 1.1, min: 0.6, max: 25, group: 'Converter', role: 'bc', help: 'Must exceed the peak line voltage: about 1.1 kV for 690 V, 5.4 kV for 3.3 kV.' },
    { key: 'fsw', label: 'Switching frequency', unit: 'Hz', type: 'number', default: 2500, min: 300, max: 10000, group: 'Converter', role: 'control' },
    { key: 'nPar', label: 'Converter modules in parallel', type: 'int', default: 6, min: 1, max: 40, group: 'Converter', role: 'geometry', help: 'Generic 1700 V / 1000 A-class IGBT modules: enough in parallel to carry the rated current with margin.' },
    { key: 'iModule', label: 'Rated rms current per module', unit: 'A', type: 'number', default: 900, min: 100, max: 4000, group: 'Converter', role: 'criteria' },
    { key: 'uk', label: 'Transformer short-circuit voltage', unit: '%', type: 'number', default: 6, min: 2, max: 16, group: 'Transformer', role: 'material' },
    { key: 'pk', label: 'Transformer load loss at rated current / rated power', unit: '%', type: 'number', default: 0.9, min: 0.1, max: 3, group: 'Transformer', role: 'material' },
    { key: 'pfe', label: 'Transformer no-load loss / rated power', unit: '%', type: 'number', default: 0.12, min: 0.02, max: 1, group: 'Transformer', role: 'material' },
    site('A', { role: 'bc' }), site('k', { role: 'bc' }),
    { key: 'event', label: 'Torsional event', type: 'select', default: 'gridloss', group: 'Dynamic response', role: 'ic', options: [{ value: 'gridloss', label: 'Loss of generator torque from rated (grid loss)' }, { value: 'step', label: 'Generator torque step of 20 % at rated' }], help: 'Initial state: steady rated operation; the event occurs at t = 0.5 s with constant aerodynamic torque (pitch has not yet reacted).' },
    { key: 'revLimit', label: 'Allowable shaft torque reversal / rated torque', unit: '-', type: 'number', default: 0.75, min: 0, max: 1.5, group: 'Acceptance', role: 'criteria', help: 'Agreed with the gearbox and coupling supplier; reversals above about half the rated torque need a managed torque ramp-down.' },
    { key: 'etaTarget', label: 'Required efficiency at rated power (shaft to medium voltage)', unit: '%', type: 'number', default: 90, min: 80, max: 99, group: 'Acceptance', role: 'criteria', help: `From rotor shaft to the medium-voltage side of the turbine transformer. The NREL 5-MW definition lumps shaft-to-generator-terminal losses into 94.4 %. ${CONFIRM}` },
    { key: 'tSim', label: 'Simulated time of the torsional event', unit: 's', type: 'number', default: 4, min: 1, max: 20, group: 'Solver', role: 'control' },
    { key: 'dt', label: 'Time step', unit: 's', type: 'number', default: 0.005, min: 0.0005, max: 0.05, group: 'Solver', role: 'control' },
  ],
  meshTime: { key: 'dt', label: 'Time step (s)', levels: [0.02, 0.01, 0.005], qoi: 'shaftMin', isSize: true },
  compute(inp) {
    const tb = getTurbine(inp), rotor = tb.rotor, s = tb.schedule;
    const N = inp.concept === 'direct' ? 1 : inp.ratio;
    const stages = inp.concept === 'geared' ? [0.01, 0.01, 0.015] : inp.concept === 'medium' ? [0.01, 0.01] : [];
    const Prated = rotor.ratedPower, etaLump = rotor.drivetrainEfficiency, PshaftR = Prated / etaLump;   // rotor shaft power at rated
    const wR = rotor.ratedRotorSpeed, wg = wR * N;
    // generator parameters from per-unit data on the rated apparent power
    const V = inp.Vll * 1e3, S = Prated / 0.95, Zb = V * V / S, we = inp.polePairs * wg;
    need(we / (2 * Math.PI) > 2 && we / (2 * Math.PI) < 400, `The electrical frequency at rated speed is ${fmt(we / (2 * Math.PI))} Hz: adjust the pole pairs (${inp.polePairs}) to the generator speed of ${fmt(wg / RPM)} rpm (10–120 Hz is usual).`);
    const psiM = inp.emf * Math.SQRT2 * V / Math.sqrt(3) / we, L = inp.xd * Zb / we, Rs = inp.rs * Zb;
    const ke = inp.iron / 100 * Prated / (we * we) * 0.6, kh = inp.iron / 100 * Prated / we * 0.4, kw = inp.fric / 100 * Prated / wg ** 3;
    const machine = { polePairs: inp.polePairs, psiM, Ld: L, Lq: L, Rs, kh, ke, windage: kw };
    need(inp.Vdc * 1e3 > V * Math.SQRT2 * 1.02, `The DC-link voltage (${inp.Vdc} kV) must exceed the peak line voltage of the generator (${fmt(V * Math.SQRT2 / 1e3)} kV).`);
    const chain = (Paero, omegaR) => {
      const gb = stages.length || inp.p0 > 0 ? gearboxEfficiency({ stageLosses: stages, noLoadLossFraction: inp.concept === 'direct' ? 0 : inp.p0 / 100, power: Paero, ratedPower: PshaftR }) : { efficiency: 1, outputPower: Paero };
      const Pg = Math.max(gb.outputPower, 0), wgen = omegaR * N;
      if (Pg <= 0 || wgen <= 0) return null;
      // shaft power into the generator = |Te·ω| + iron + mechanical losses: iterate the electromagnetic torque
      let Te = -Pg / wgen, g = null;
      for (let it = 0; it < 6; it++) { g = pmsmSteadyState({ ...machine, speed: wgen, torque: Te, strategy: 'mtpa' }); Te = -(Pg - g.ironLoss - g.mechanicalLoss) / wgen; }
      const Pel = -g.electricalPower;
      if (!(Pel > 0)) return null;
      const m = Math.min(1.15, g.Vll * Math.SQRT2 / (inp.Vdc * 1e3) * 2 / Math.sqrt(3));
      const cg = converterLoss({ Irms: g.Irms, Vdc: inp.Vdc * 1e3, modulationIndex: Math.min(m, 1), cosPhi: -Math.abs(g.powerFactor), fsw: inp.fsw, parallel: inp.nPar });
      const Igrid = (Pel - cg.total) / (Math.sqrt(3) * V * 0.98);
      const cn = converterLoss({ Irms: Math.max(Igrid, 0), Vdc: inp.Vdc * 1e3, modulationIndex: 0.9, cosPhi: 0.98, fsw: inp.fsw, parallel: inp.nPar });
      const Pconv = Pel - cg.total - cn.total;
      const tr = transformer({ Sr: S, uk: inp.uk / 100, Pk: inp.pk / 100 * Prated, P0: inp.pfe / 100 * Prated, loading: Math.max(Pconv, 0) / 0.98 / S, powerFactor: 0.98 });
      const Pout = Pconv - tr.losses;
      return { Paero, gb, g, cg, cn, tr, Pg, Pel, Pconv, Pout, eta: Pout / Paero, m };
    };
    const rated = run(() => chain(PshaftR, wR), 'Rated operating point');
    need(rated && rated.Pout > 0, 'The conversion chain delivers no power at rated: check the generator and converter data.');
    // efficiency over the operating schedule and Weibull-weighted annual efficiency
    const pts = s.U.map((U, i) => { const c = chain(s.aeroPower[i], s.omega[i]); return { U, c }; }).filter((q) => q.c && q.c.Pout > 0);
    const wb = weibull(inp.weibA, inp.weibK);
    let eIn = 0, eOut = 0;
    for (let i = 0; i < pts.length; i++) { const lo = i ? 0.5 * (pts[i].U + pts[i - 1].U) : pts[i].U - 0.25, hi = i < pts.length - 1 ? 0.5 * (pts[i].U + pts[i + 1].U) : pts[i].U + 0.25, w = wb.cdf(hi) - wb.cdf(lo); eIn += w * pts[i].c.Paero; eOut += w * pts[i].c.Pout; }
    const etaAnnual = eIn > 0 ? eOut / eIn : 0;
    // control: optimal torque gain on the generator shaft
    const opt = optimumCp();
    const kOpt = optimalTorqueGain({ rho: 1.225, R: rotor.rotorRadius, CpMax: opt.cpMax, lambdaOpt: opt.tsrOpt, gearRatio: N }).kGenerator;
    // torsional event: two inertias on the low-speed shaft, constant aerodynamic torque
    const Jr = rotor.rotorInertia, Jg = rotor.generatorInertia * 97 ** 2 * (inp.concept === 'direct' ? 4 : 1);   // generator inertia referred to the low-speed shaft (direct-drive ring generator taken four times larger)
    const k = inp.kShaft * 1e6, cs = inp.cShaft * 1e6, T0 = PshaftR / wR, Tg1 = inp.event === 'gridloss' ? 0 : 0.8 * T0, tEv = 0.5;
    const n = Math.round(inp.tSim / inp.dt);
    need(n <= 2e5 && n >= 20, 'Choose a time step and duration that give between 20 and 200 000 steps.');
    const wn = Math.sqrt(k * (1 / Jr + 1 / Jg));
    need(wn * inp.dt < 1.5, `The time step is too large for the torsional mode (${fmt(wn / (2 * Math.PI))} Hz): use ${fmt(0.5 / wn, 2)} s or less.`);
    // state: [twist φ, twist rate, rotor speed]
    let y = [T0 / k, 0, wR];
    const f = (t, q) => { const Tg = t < tEv ? T0 : Tg1, Ts = k * q[0] + cs * q[1]; return [q[1], (T0 - Ts) / Jr - (Ts - Tg) / Jg, (T0 - Ts) / Jr]; };
    const tt = new Float64Array(n + 1), ts = new Float64Array(n + 1), wgS = new Float64Array(n + 1);
    let TsMin = T0, TsMax = T0;
    for (let i = 0; i <= n; i++) {
      const t = i * inp.dt, Ts = k * y[0] + cs * y[1];
      tt[i] = t; ts[i] = Ts; wgS[i] = (y[2] - y[1]) * N; TsMin = Math.min(TsMin, Ts); TsMax = Math.max(TsMax, Ts);
      const k1 = f(t, y), k2 = f(t + inp.dt / 2, y.map((v, j) => v + inp.dt / 2 * k1[j])), k3 = f(t + inp.dt / 2, y.map((v, j) => v + inp.dt / 2 * k2[j])), k4 = f(t + inp.dt, y.map((v, j) => v + inp.dt * k3[j]));
      y = y.map((v, j) => v + inp.dt / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]));
    }
    const overspeed1s = (() => { const i = Math.min(n, Math.round((tEv + 1) / inp.dt)); return wgS[i] / (wR * N); })();
    const zeta = cs * wn / (2 * k);
    const Icap = inp.nPar * inp.iModule;
    const checks = [
      check('Conversion efficiency at rated power', inp.etaTarget, rated.eta * 100, '%', 'IEC 60034-2-1 / IEC 61400-4 loss accounting; project requirement', 'Demand is the required efficiency from rotor shaft to the medium-voltage terminals, capacity the calculated value. A shortfall is energy lost as heat in gearbox, generator, converter and transformer — and cooling capacity to match.'),
      check('Generator current against the converter rating', rated.g.Irms, Icap, 'A', 'IEC 62477-1 / converter module rating', 'Rated stator current against the continuous current capability of the parallel converter modules. Over-current trips the converter or overheats the semiconductors.'),
      check('Generator terminal voltage against the DC-link limit', rated.g.Vll * Math.SQRT2 / 1e3, inp.Vdc * Math.sqrt(3) / 2 * 1.155, 'kV', 'Sinusoidal PWM with third-harmonic injection: peak line voltage ≤ Vdc', 'Peak line voltage of the generator at rated load against what the bridge can synthesise from the DC link. If exceeded the converter loses current control (field weakening would be needed).'),
      check('Shaft torque reversal after the torsional event', Math.max(0, -TsMin) / T0, Math.max(inp.revLimit, 1e-6), '-', 'IEC 61400-4:2012 §6 (torque reversals); drivetrain design practice', 'Magnitude of negative shaft torque after the event relative to rated torque, against the agreed allowable. Torque reversals unload and reload gear teeth and bearings: a main cause of gearbox damage.'),
    ];
    const gov = worst(checks);
    const lossR = [['Gearbox', PshaftR - rated.Pg], ['Generator copper', rated.g.copperLoss], ['Generator iron', rated.g.ironLoss], ['Generator friction', rated.g.mechanicalLoss], ['Generator-side converter', rated.cg.total], ['Grid-side converter', rated.cn.total], ['Transformer', rated.tr.losses]];
    const hours = 8760 * 0.45, dE = (v) => v * hours * Prated / 1e6;   // MWh/yr at a 45 % capacity factor per unit of efficiency
    const [tp, tsp, wgp] = thin([tt, ts, wgS], 800);
    return {
      kpis: [
        { key: 'etaRated', label: 'Efficiency at rated power (shaft to medium voltage)', value: round(rated.eta * 100, 4), unit: '%', limit: inp.etaTarget, sense: 'min', explain: 'Electrical power after the transformer over rotor shaft power at rated operation. Geared PMSG chains reach 93–95 %, direct drive 94–96 %.' },
        { key: 'etaAnnual', label: 'Annual energy-weighted conversion efficiency', value: round(etaAnnual * 100, 4), unit: '%', explain: 'Energy delivered over shaft energy across the wind-speed distribution. Part-load efficiency matters more than the rated value because the turbine spends most hours below rated.' },
        { key: 'etaGear', label: 'Gearbox efficiency at rated', value: round(rated.gb.efficiency * 100, 4), unit: '%', explain: 'Product of the stage efficiencies less the no-load loss (100 % for direct drive).' },
        { key: 'etaGen', label: 'Generator efficiency at rated', value: round(rated.Pel / rated.Pg * 100, 4), unit: '%', explain: 'Electrical output over shaft input of the generator: copper, iron and friction losses.' },
        { key: 'etaConv', label: 'Converter efficiency at rated (both bridges)', value: round(rated.Pconv / rated.Pel * 100, 4), unit: '%', explain: 'Output of the grid-side bridge over generator output: conduction and switching losses of the back-to-back converter.' },
        { key: 'etaTrafo', label: 'Transformer efficiency at rated', value: round(rated.tr.efficiency * 100, 4), unit: '%', explain: 'From the load and no-load losses of the turbine transformer at this loading.' },
        { key: 'iRated', label: 'Generator phase current at rated', value: round(rated.g.Irms, 5), unit: 'A', limit: Icap, sense: 'max', explain: 'Stator rms current at rated torque; it sizes the generator winding, cables in the tower and the converter.' },
        { key: 'vRated', label: 'Generator line voltage at rated', value: round(rated.g.Vll, 5), unit: 'V', explain: 'Terminal voltage from back-EMF plus the reactance and resistance drops at rated current.' },
        { key: 'freq', label: 'Electrical frequency at rated speed', value: round(rated.g.frequency, 4), unit: 'Hz', explain: 'Stator frequency p·n/60; the full converter decouples it from the grid frequency.' },
        { key: 'torqueGen', label: 'Generator electromagnetic torque at rated', value: round(Math.abs(rated.g.torque) / 1e3, 5), unit: 'kN·m', explain: 'Air-gap torque 1.5·p·ψm·iq: the reaction torque that the controller commands.' },
        { key: 'kOpt', label: 'Optimal torque gain on the generator shaft', value: round(kOpt, 5), unit: 'N·m/(rad/s)²', explain: 'Gain of the torque law T = k·Ω² that keeps the rotor at its optimum tip-speed ratio below rated (maximum power point tracking).' },
        { key: 'fTors', label: 'Drivetrain torsional natural frequency', value: round(wn / (2 * Math.PI), 4), unit: 'Hz', explain: `Free–free frequency of rotor against generator on the shaft (damping ratio ${fmt(zeta * 100)} %).` },
        { key: 'shaftMin', label: 'Minimum shaft torque after the event / rated torque', value: round(TsMin / T0, 5), unit: '-', explain: 'Lowest shaft torque in the oscillation that follows the event. Negative values are torque reversals; an undamped drivetrain swings to (J_g − J_r)/(J_g + J_r) after a full loss of generator torque.' },
        { key: 'overspeed', label: 'Generator speed one second after the event / rated', value: round(overspeed1s, 5), unit: '-', explain: 'Speed rise while the aerodynamic torque is unbalanced and before the pitch system acts: it sets the required pitch rate and the over-speed protection margin.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Efficiency of the conversion chain versus wind speed', x: { label: 'Hub-height wind speed', unit: 'm/s' }, y: { label: 'Efficiency', unit: '%' }, series: [{ name: 'Total (shaft to medium voltage)', x: pts.map((q) => q.U), y: pts.map((q) => round(q.c.eta * 100, 3)) }, { name: 'Gearbox', x: pts.map((q) => q.U), y: pts.map((q) => round(q.c.gb.efficiency * 100, 3)) }, { name: 'Generator', x: pts.map((q) => q.U), y: pts.map((q) => round(q.c.Pel / q.c.Pg * 100, 3)) }, { name: 'Converter', x: pts.map((q) => q.U), y: pts.map((q) => round(q.c.Pconv / q.c.Pel * 100, 3)) }, { name: 'Transformer', x: pts.map((q) => q.U), y: pts.map((q) => round(q.c.tr.efficiency * 100, 3)) }],
          explain: 'Efficiency of each stage and of the whole chain over the operating wind speeds. Constant (no-load) losses pull the curve down at low wind; where it flattens is where the chain works best.' },
        { type: 'bar', title: 'Loss breakdown at rated power', categories: lossR.map((q) => q[0]), series: [{ name: 'Loss (kW)', y: lossR.map((q) => round(q[1] / 1e3, 4)) }], y: { label: 'Loss', unit: 'kW' },
          explain: 'Where the power is lost at rated operation. The largest bars are the targets for efficiency improvement and set the cooling demand of each component.' },
        { type: 'line', title: 'Generator torque–speed schedule', x: { label: 'Generator speed', unit: 'rpm' }, y: { label: 'Electromagnetic torque', unit: 'kN·m' }, series: [{ name: 'Operating schedule', x: pts.map((q) => round(q.c.g.frequency * 60 / inp.polePairs, 3)), y: pts.map((q) => round(Math.abs(q.c.g.torque) / 1e3, 4)) }, { name: 'Optimal torque law k·Ω²', x: pts.map((q) => round(q.c.g.frequency * 60 / inp.polePairs, 3)), y: pts.map((q) => round(kOpt * (q.c.g.frequency * 2 * Math.PI / inp.polePairs) ** 2 / 1e3, 4)), dash: true }],
          explain: 'Generator torque against speed along the operating schedule with the quadratic optimal-torque law. The schedule follows the law in the variable-speed region and rises steeply to rated torque at rated speed.' },
        { type: 'line', title: 'Shaft torque and generator speed after the torsional event', x: { label: 'Time', unit: 's' }, y: { label: 'Per unit of rated', unit: '-' }, series: [{ name: 'Shaft torque / rated', x: tp.map((v) => round(v, 4)), y: tsp.map((v) => round(v / T0, 4)) }, { name: 'Generator speed / rated', x: tp.map((v) => round(v, 4)), y: wgp.map((v) => round(v / (wR * N), 4)) }],
          explain: 'Low-speed shaft torque and generator speed after the event at 0.5 s. The shaft rings at the torsional frequency; with a full loss of generator torque the torque swings negative and the speed climbs until the pitch system intervenes.' },
      ],
      tables: [
        { title: 'Operating points', columns: [{ key: 'U', label: 'Wind speed', unit: 'm/s' }, { key: 'Ps', label: 'Shaft power', unit: 'kW' }, { key: 'n', label: 'Generator speed', unit: 'rpm' }, { key: 'T', label: 'Torque', unit: 'kN·m' }, { key: 'I', label: 'Current', unit: 'A' }, { key: 'V', label: 'Line voltage', unit: 'V' }, { key: 'Po', label: 'Output', unit: 'kW' }, { key: 'eta', label: 'Efficiency', unit: '%' }],
          rows: pts.filter((_, i) => i % 2 === 0).map((q) => ({ U: q.U, Ps: round(q.c.Paero / 1e3, 5), n: round(q.c.g.frequency * 60 / inp.polePairs, 4), T: round(Math.abs(q.c.g.torque) / 1e3, 4), I: round(q.c.g.Irms, 4), V: round(q.c.g.Vll, 4), Po: round(q.c.Pout / 1e3, 5), eta: round(q.c.eta * 100, 3) })),
          explain: 'Shaft power, generator speed, torque, current, voltage and delivered power at each wind speed. The output column is an electrical power curve at the medium-voltage terminals.' },
        { title: 'Generator parameters (derived from per-unit data)', columns: [{ key: 'q', label: 'Parameter' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [{ q: 'Rated apparent power', v: round(S / 1e6, 4), u: 'MVA' }, { q: 'Permanent-magnet flux linkage ψm', v: round(psiM, 5), u: 'Wb' }, { q: 'Synchronous inductance Ld = Lq', v: round(L * 1e3, 5), u: 'mH' }, { q: 'Stator resistance', v: round(Rs * 1e3, 5), u: 'mΩ' }, { q: 'Electrical time constant L/R', v: round(L / Rs * 1e3, 4), u: 'ms' }, { q: 'Modulation index at rated', v: round(rated.m, 4), u: '-' }, { q: 'Generator inertia on the low-speed shaft', v: round(Jg / 1e6, 5), u: '10⁶ kg·m²' }, { q: 'Rotor inertia', v: round(Jr / 1e6, 5), u: '10⁶ kg·m²' }],
          explain: 'Equivalent-circuit parameters of the generator and the inertias used in the torsional model. The electrical time constant bounds how fast torque can be changed by the converter.' },
      ],
      notes: [
        `${rotor.name}, ${inp.concept === 'direct' ? 'direct drive' : `gearbox 1:${N}`}: rated shaft power ${fmt(PshaftR / 1e6, 4)} MW at ${fmt(wR / RPM)} rpm → ${fmt(rated.Pout / 1e6, 4)} MW at medium voltage (${fmt(rated.eta * 100, 4)} %; the rotor definition assumes ${fmt(etaLump * 100)} %).`,
        `Generator: ${fmt(rated.g.Irms)} A, ${fmt(rated.g.Vll)} V, ${fmt(rated.g.frequency)} Hz, power factor ${fmt(Math.abs(rated.g.powerFactor))}; converter ${inp.nPar} modules (${fmt(Icap)} A).`,
        `Torsional mode ${fmt(wn / (2 * Math.PI))} Hz with ${fmt(zeta * 100)} % damping; after the event the shaft torque swings to ${fmt(TsMin / T0)} of rated and the generator speed is ${fmt(overspeed1s)} of rated after one second.`,
        `Governing check: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `The chain converts ${fmt(rated.eta * 100, 4)} % of shaft power at rated and ${fmt(etaAnnual * 100, 4)} % over the year; after a ${inp.event === 'gridloss' ? 'loss of generator torque' : 'torque step'} the shaft torque reaches ${fmt(TsMin / T0)} of rated.`,
          gov === checks[3] && gov.uc > 1 ? `Limit torque reversals: ramp the generator torque down over about ${fmt(2 * Math.PI / wn * 1.5, 2)} s (one and a half torsional periods) through the DC-link chopper instead of dropping it, or add an active drivetrain damper in the torque loop.` : gov.uc > 1 ? (gov === checks[0] ? `Reduce the largest loss (${lossR.reduce((a, b) => (b[1] > a[1] ? b : a))[0]}, ${fmt(lossR.reduce((a, b) => (b[1] > a[1] ? b : a))[1] / 1e3)} kW) — or accept ${fmt(rated.eta * 100, 4)} % in the yield assessment.` : gov === checks[1] ? `Add converter modules: ${Math.ceil(rated.g.Irms / inp.iModule * 1.1)} are needed for ${fmt(rated.g.Irms)} A.` : `Raise the DC-link voltage to about ${fmt(rated.g.Vll * Math.SQRT2 / 1e3 * 1.1)} kV or lower the back-EMF ratio.`) : 'The electrical chain is within its ratings. Replace the per-unit generator data and the generic device losses by supplier data, and check fault ride-through and harmonic performance per IEC 61400-21-1.',
          'IEC 61400-1:2019 §9–10 with IEC 60034 loss definitions and IEC 61400-4 for the gearbox; fault ride-through per grid code and IEC 61400-21-1.',
          `Losses at rated are ${fmt((1 - rated.eta) * PshaftR / 1e3)} kW of heat to be removed by the cooling system. A permanent-magnet generator of this size uses rare-earth magnets (order of 100 kg/MW for geared, 600 kg/MW for direct drive) and several tonnes of copper (≈ ${CO2_COPPER} t CO₂e/t): efficiency gains should be weighed against material intensity.`),
        adv('economic', 'info', `One percentage point of conversion efficiency is worth about ${fmt(dE(0.01))} MWh a year per turbine at a 45 % capacity factor; the annual energy-weighted efficiency is ${fmt(etaAnnual * 100, 4)} % against ${fmt(etaLump * 100)} % assumed in the rotor model.`,
          `Use the annual value in the yield assessment${Math.abs(etaAnnual - etaLump) > 0.005 ? ` (the difference changes the energy by ${fmt((etaAnnual / etaLump - 1) * 100)} %)` : ''}. Compare drivetrain concepts on whole-life cost: gearbox replacements (about one per turbine life for three-stage designs) against the heavier direct-drive generator.`,
          'LCOE-based drivetrain comparison (IEA Wind Task 37 practice) with whole-life costing per ISO 15686-5; efficiency valued at the energy price over the design life.',
          energyCo2Text(dE(0.01))),
      ],
      publish: { 'drivetrain.efficiencyRated': rated.eta, 'drivetrain.efficiencyAnnual': etaAnnual, 'drivetrain.shaftTorqueMin': TsMin },
    };
  },
  verify: [
    { name: 'Undamped two-inertia drivetrain after loss of generator torque: minimum shaft torque / rated = (J_g − J_r)/(J_g + J_r)', source: 'Free vibration of a two-mass torsional oscillator about its new equilibrium (hand calculation with NREL 5-MW inertias)',
      // J_r = 38 759 228 kg·m², J_g = 534.116 × 97² = 5 025 497 kg·m²: equilibrium torque T·J_g/(J_r + J_g) = 0.11478 T, swing amplitude T·J_r/(J_r + J_g) = 0.88522 T → minimum −0.7704 T
      inputs: { event: 'gridloss', cShaft: 0, tSim: 3, dt: 0.002 }, expect: [{ kpi: 'shaftMin', value: -0.7704, tol: 0.003 }] },
    { name: 'Three-stage gearbox at rated: η = 0.99 × 0.99 × 0.985 − 0.005', source: 'Stage-loss product less the no-load loss (IEC 61400-4 loss accounting; hand calculation)',
      // 0.99·0.99·0.985 = 0.96540; − 0.005 = 0.96040 → 96.04 %
      inputs: { concept: 'geared', p0: 0.5 }, expect: [{ kpi: 'etaGear', value: 96.04, tol: 0.0005 }] },
  ],
};
