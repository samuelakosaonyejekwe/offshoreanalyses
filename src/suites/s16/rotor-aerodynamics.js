// Offshore wind turbine aerodynamic loading: blade-element-momentum solution of the rotor at one operating point
// (thrust, torque, mechanical power, blade-root moments, spanwise loading) and the steady thrust / power schedule
// over the operating wind-speed range.
import { bemSolve, maxCp, surfaceLookup, BETZ_LIMIT, glauertIdealCp } from '../../engine/aero/bem.js';
import { nrel5mw, createRotor } from '../../engine/aero/rotors.js';
import { genericPolar, cylinderPolar } from '../../engine/aero/polar.js';
import { unsteadyBem, harmonicAmplitudes } from '../../engine/aero/unsteady.js';
import { region2Gain } from '../../engine/aero/control.js';
import { tipSpeedCheck } from '../../engine/aero/blade.js';
import { turbineInputs, getTurbine, buildTurbine, optimumCp, steadyAt, adv, check, sev, fmt, round, run, need, linspace, RPM, DEG, LIMIT_POLAR, LIMIT_CODE, CFD_OUT, SUPPLIER, CO2_GFRP, energyCo2Text } from './_common.js';

// Blade planform of the NREL 5-MW reference at its 17 analysis nodes (Jonkman et al. 2009, Table 3-1); t/c = 1 marks the cylindrical root.
const REF_BLADE = nrel5mw().stations.map((q) => ({ r: q.r, chord: q.chord, twist: q.twistDeg, tc: q.thickness }));
const tab = { key: 'bladeMode', in: ['table'] };

/** Rotor from a user blade table: generic parametric polars by relative thickness, gross properties scaled from the reference. */
function customRotor(inp) {
  const rows = (inp.blade || []).filter((q) => q.r > 0 && q.chord > 0);
  need(rows.length >= 5, 'The blade table needs at least five stations (radius, chord, twist, relative thickness).');
  for (let i = 1; i < rows.length; i++) need(rows[i].r > rows[i - 1].r, `Blade table: radii must ascend (row ${i + 1}).`);
  need(rows[0].r > inp.hubR && rows[rows.length - 1].r < inp.tipR, `Blade table: stations must lie between the hub radius (${inp.hubR} m) and the tip radius (${inp.tipR} m).`);
  need(rows[rows.length - 1].r >= 0.9 * inp.tipR && rows[0].r <= inp.hubR + 0.3 * (inp.tipR - inp.hubR), `Blade table: the stations must cover the blade from near the root to at least 90 % of the tip radius (${fmt(0.9 * inp.tipR)} m); they run from ${fmt(rows[0].r)} m to ${fmt(rows[rows.length - 1].r)} m.`);
  need(rows.every((q) => q.tc > 0.05 && q.tc <= 1 && Math.abs(q.twist) <= 45), 'Blade table: relative thickness must be between 0.05 and 1 (1 = cylinder) and twist within ±45°.');
  const ref = nrel5mw(), sc = inp.tipR / ref.tipRadius, { stations, polars, ...base } = ref;
  const wR = inp.rpmRated * RPM, P = inp.ratedCustomMW * 1e6;
  return run(() => createRotor({ ...base, name: 'User blade geometry', label: 'user planform with generic parametric polars by thickness', blades: inp.nBlades, hubRadius: inp.hubR, tipRadius: inp.tipR, hubHeight: inp.hubHeight,
    ratedPower: P, ratedRotorSpeed: wR, maxRotorSpeed: wR, minRotorSpeed: 0.57 * wR, bladeLength: inp.tipR - inp.hubR, bladeMass: ref.bladeMass * sc ** 3 * inp.nBlades / inp.nBlades, rnaMass: ref.rnaMass * sc ** 3,
    rotorInertia: ref.rotorInertia * sc ** 5, drivetrainInertia: ref.drivetrainInertia * sc ** 5,
    stations: rows.map((q) => ({ r: q.r, chord: q.chord, twistDeg: q.twist, thickness: q.tc, airfoil: q.tc >= 0.99 ? 'cylinder' : `generic ${Math.round(q.tc * 100)} %`, polar: q.tc >= 0.99 ? cylinderPolar({ cd: 0.5 }) : genericPolar({ thickness: Math.min(Math.max(q.tc, 0.18), 0.4) }) })) }), 'Rotor definition');
}

export default {
  id: 's16-rotor-aerodynamics', suite: 16, title: 'Rotor aerodynamic loading (BEM)',
  summary: 'Rotor thrust, torque, mechanical power and blade loading at an operating point from blade-element-momentum theory, with the thrust and power schedule over the wind-speed range.',
  covers: [4], fidelity: 'numerical',
  standards: ['IEC 61400-1:2019 §6.3 (wind conditions) and §7.4 (design load cases)', 'IEC 61400-3-1:2019 (fixed offshore wind turbines)', 'NREL/TP-500-38060 (Jonkman et al. 2009, 5-MW reference turbine)', 'NREL/TP-500-36881 (Moriarty & Hansen 2005, AeroDyn theory)'],
  models: ['Blade element momentum theory (Ning 2014 residual)', 'Unsteady blade-element momentum with Øye dynamic inflow', 'Beddoes–Leishman dynamic stall', 'Potential-flow tower shadow', 'Actuator-disc momentum theory and Betz limit', 'Prandtl tip- and hub-loss corrections', 'Glauert / Buhl high-thrust correction', 'Glauert ideal rotor with wake rotation', 'Skewed-wake (Pitt–Peters) correction for yaw and tilt', 'Power-law wind shear', 'Tip-speed ratio and power coefficient formulations', 'Viterna post-stall extrapolation of generic polars'],
  limits: `${LIMIT_POLAR} Rigid straight blades: aeroelastic deflection and turbulence are not in this module (see the coupled dynamics module); the unsteady analysis is deterministic (shear, tilt, tower shadow) with independent blades. Lifting-line and vortex-lattice methods are not provided. ${CFD_OUT} ${LIMIT_CODE}`,
  publishes: ['owt.thrustRated', 'owt.rna', 'owt.ratedPower', 'owt.rotorDiameter', 'owt.hubHeight'],
  searchTerms: 'blade element momentum wind turbine rotor thrust NREL 5-MW reference',
  inputs: [
    ...turbineInputs(),
    { key: 'bladeMode', label: 'Blade geometry', type: 'select', default: 'reference', group: 'Blade geometry', role: 'geometry', options: [{ value: 'reference', label: 'Rotor model selected above' }, { value: 'table', label: 'Entered blade table' }],
      help: 'A blade description (span stations, chord, twist, relative thickness of the airfoil sections) can be pasted. Sections use generic parametric polars by thickness.' },
    { key: 'blade', label: 'Blade stations', type: 'table', group: 'Blade geometry', role: 'geometry', showIf: tab, columns: [{ key: 'r', label: 'Radius', unit: 'm', type: 'number' }, { key: 'chord', label: 'Chord', unit: 'm', type: 'number' }, { key: 'twist', label: 'Twist', unit: 'deg', type: 'number' }, { key: 'tc', label: 'Relative thickness t/c (1 = cylinder)', unit: '-', type: 'number' }],
      default: REF_BLADE.map((q) => ({ ...q })), help: 'Default: the NREL 5-MW planform. Radii measured along the blade from the rotor axis, ascending, strictly between hub and tip radius.' },
    { key: 'nBlades', label: 'Number of blades', type: 'int', default: 3, min: 2, max: 4, group: 'Blade geometry', role: 'geometry', showIf: tab },
    { key: 'hubR', label: 'Hub radius', unit: 'm', type: 'number', default: 1.5, min: 0.2, max: 8, group: 'Blade geometry', role: 'geometry', showIf: tab },
    { key: 'tipR', label: 'Tip radius', unit: 'm', type: 'number', default: 63, min: 5, max: 160, group: 'Blade geometry', role: 'geometry', showIf: tab, design: true },
    { key: 'ratedCustomMW', label: 'Rated electrical power', unit: 'MW', type: 'number', default: 5, min: 0.1, max: 30, group: 'Blade geometry', role: 'geometry', showIf: tab },
    { key: 'rpmRated', label: 'Rated rotor speed', unit: 'rpm', type: 'number', default: 12.1, min: 3, max: 40, group: 'Blade geometry', role: 'control', showIf: tab },
    { key: 'U', label: 'Hub-height wind speed', unit: 'm/s', type: 'number', default: 11.4, min: 1, max: 40, group: 'Inflow', role: 'bc', design: true, uncertain: { dist: 'normal', cov: 0.05 }, help: 'Undisturbed 10-minute mean wind speed at hub height (far-field inflow boundary). 11.4 m/s is the rated speed of the NREL 5-MW reference.' },
    { key: 'shearExp', label: 'Wind-shear exponent', unit: '-', type: 'number', default: 0.14, min: 0, max: 0.4, group: 'Inflow', role: 'bc', help: `Power-law exponent: IEC 61400-3-1 uses 0.14 offshore for normal conditions (0.11 for extreme).` },
    { key: 'yawDeg', label: 'Yaw misalignment', unit: 'deg', type: 'number', default: 0, min: -30, max: 30, group: 'Inflow', role: 'bc', help: 'Angle between wind and rotor axis; the skewed-wake correction is reliable to about 30°.' },
    { key: 'rhoAir', label: 'Air density', unit: 'kg/m³', type: 'number', default: 1.225, min: 1.0, max: 1.35, group: 'Inflow', role: 'environment', uncertain: { dist: 'normal', cov: 0.02 } },
    { key: 'opMode', label: 'Operating state', type: 'select', default: 'schedule', group: 'Operating state', role: 'ic',
      options: [{ value: 'schedule', label: 'From the variable-speed pitch-regulated schedule' }, { value: 'manual', label: 'Prescribed rotor speed and pitch' }],
      help: 'Initial turbine state. The schedule follows optimum tip-speed ratio below rated and pitches to feather above.' },
    { key: 'rpm', label: 'Rotor speed', unit: 'rpm', type: 'number', default: 12.1, min: 0.5, max: 30, group: 'Operating state', role: 'ic', showIf: { key: 'opMode', in: ['manual'] } },
    { key: 'pitchDeg', label: 'Collective blade pitch', unit: 'deg', type: 'number', default: 0, min: -5, max: 90, group: 'Operating state', role: 'ic', showIf: { key: 'opMode', in: ['manual'] }, help: 'Positive to feather.' },
    { key: 'tipLoss', label: 'Prandtl tip and hub loss', type: 'bool', default: true, group: 'Model options', role: 'control' },
    { key: 'unsteady', label: 'Unsteady blade loads over a revolution (dynamic inflow, Beddoes–Leishman dynamic stall, shear, tower shadow)', type: 'bool', default: true, group: 'Model options', role: 'control' },
    { key: 'nAz', label: 'Time steps per revolution', type: 'int', default: 48, min: 12, max: 384, group: 'Model options', role: 'control', showIf: { key: 'unsteady', in: [true] }, help: 'Resolution of the unsteady blade-element solution; the tower passage needs about 40 or more.' },
    { key: 'tipLimit', label: 'Allowable blade-tip speed', unit: 'm/s', type: 'number', default: 90, min: 60, max: 120, group: 'Acceptance', role: 'criteria', help: 'Project limit for leading-edge erosion and noise; offshore rotors are commonly limited to 80–100 m/s.' },
  ],
  mesh: { key: 'nAz', label: 'Time steps per revolution', levels: [24, 48, 96], qoi: 'flap1P' },
  compute(inp) {
    const custom = inp.bladeMode === 'table';
    const tb = custom ? buildTurbine(customRotor(inp), { rho: inp.rhoAir }) : getTurbine(inp, { rho: inp.rhoAir }), rotor = tb.rotor, R = rotor.tipRadius;
    need(inp.hubHeight > R + 5, `The hub height (${inp.hubHeight} m) leaves no clearance under the blade tip (rotor radius ${fmt(R)} m).`);
    let omega, pitch;
    if (inp.opMode === 'manual') { omega = inp.rpm * RPM; pitch = inp.pitchDeg; } else {
      const st = steadyAt(tb, inp.U);
      need(st.omega > 0, `At ${inp.U} m/s the turbine is outside its operating range (cut-in ${rotor.cutIn} m/s, cut-out ${rotor.cutOut} m/s): choose a prescribed rotor speed and pitch to analyse a parked or idling rotor.`);
      omega = st.omega; pitch = st.pitchDeg;
    }
    const bem = run(() => bemSolve({ rotor, U: inp.U, omega, pitchDeg: pitch, rho: inp.rhoAir, yawDeg: inp.yawDeg, shearExp: inp.shearExp, hubHeight: inp.hubHeight, tipLoss: inp.tipLoss, hubLoss: inp.tipLoss }), 'BEM solution');
    const opt = custom ? run(() => maxCp({ rotor, rho: inp.rhoAir }), 'Optimum power coefficient') : optimumCp();
    // unsteady blade loads: three revolutions in sheared inflow with tower shadow, harmonics of the last one
    let un = null;
    if (inp.unsteady && omega > 0.05) {
      const nPer = inp.nAz, Trev = 2 * Math.PI / omega, nSt = 3 * nPer;
      const u = run(() => unsteadyBem({ rotor, dt: Trev / nPer, nSteps: nSt, wind: inp.U, omega, pitchDeg: pitch, rho: inp.rhoAir, shearExp: inp.shearExp, hubHeight: inp.hubHeight, towerShadow: !!rotor.tower, dynamicInflow: 'oye', dynamicStall: 'bl' }), 'Unsteady BEM');
      const last = Array.from(u.flapMoment[0].subarray(nSt - nPer)), hf = harmonicAmplitudes(last, 3), ht = harmonicAmplitudes(Array.from(u.thrust.subarray(nSt - nPer)), 2 * rotor.blades);
      un = { az: Array.from(u.azimuth.subarray(nSt - nPer), (a) => ((a / DEG) % 360 + 360) % 360), flap: last, mean: hf.mean, p1: hf.harmonics[0].amp, p2: hf.harmonics[1].amp, thrustNP: ht.harmonics[rotor.blades - 1].amp, max: Math.max(...last), min: Math.min(...last) };
    }
    const k2 = region2Gain({ R: rotor.rotorRadius, cpMax: opt.cpMax, tsrOpt: opt.tsrOpt, rho: inp.rhoAir, gearboxRatio: rotor.gearboxRatio });
    const tip = tipSpeedCheck({ omega: Math.max(omega, rotor.ratedRotorSpeed), radius: R, limit: inp.tipLimit });
    const elec = bem.power * rotor.drivetrainEfficiency;
    const sp = bem.spanwise, s = tb.schedule;
    const alphaMax = Math.max(...sp.filter((q) => q.r > 0.6 * R).map((q) => q.alpha / DEG));   // outer 40 % of the span (thin sections that carry the load)
    const checks = [
      check('Blade-tip speed at rated rotor speed', tip.tipSpeed, inp.tipLimit, 'm/s', 'Project criterion (leading-edge erosion, cf. DNV-RP-0573)', 'Tip speed against the project limit. Rain-erosion damage of the leading edge grows with a high power of impact speed, so exceeding the limit shortens blade repair intervals.'),
      check('Power coefficient against the Betz limit', bem.Cp, BETZ_LIMIT, '-', 'Betz (1920); Burton et al., Wind Energy Handbook §3.2', 'No rotor can extract more than 16/27 of the kinetic power through its swept area. A value above the limit means the operating point is outside the validity of momentum theory.'),
      check('Outboard angle of attack against stall', Math.max(alphaMax, 0), 13, 'deg', 'Generic polar stall angle of the 18–21 % sections (12.5–13°)', 'Largest angle of attack on the outer 40 % of the blade against the stall angle of the generic polars. Beyond it the section stalls: lift collapses, drag rises and loads become unsteady.'),
    ];
    const tsrs = linspace(3, 13, 41);
    const cpFine = tsrs.map((t) => surfaceLookup(tb.surface, t, tb.schedule.pitchFineDeg).cp);
    const cp4 = tsrs.map((t) => surfaceLookup(tb.surface, t, 4).cp);
    const steel = 0, gfrp = rotor.bladeMass * rotor.blades / 1e3;
    const aepLossPerPct = 0.01 * rotor.ratedPower * 0.45 * 8760 / 1e6;   // MWh/yr per 1 % of Cp at a 45 % capacity factor
    const gov = checks.reduce((a, c) => (c.uc > a.uc ? c : a));
    return {
      kpis: [
        { key: 'thrust', label: 'Rotor thrust', value: round(bem.thrust / 1e3, 5), unit: 'kN', explain: 'Axial aerodynamic force on the rotor at this operating point. It is the main overturning load on tower and foundation in power production and peaks close to rated wind speed.' },
        { key: 'torque', label: 'Aerodynamic torque (low-speed shaft)', value: round(bem.torque / 1e3, 5), unit: 'kN·m', explain: 'Torque delivered by the blades to the main shaft. It sizes the drivetrain and, multiplied by rotor speed, gives the mechanical power.' },
        { key: 'aeroPower', label: 'Mechanical (aerodynamic) power', value: round(bem.power / 1e6, 5), unit: 'MW', explain: 'Torque times rotor speed before drivetrain and electrical losses. Compare with the rated electrical power divided by the drivetrain efficiency.' },
        { key: 'elecPower', label: 'Electrical power', value: round(elec / 1e6, 5), unit: 'MW', explain: 'Mechanical power times the lumped drivetrain efficiency. Above rated wind speed the pitch controller holds it at the rating.' },
        { key: 'cp', label: 'Power coefficient Cp', value: round(bem.Cp, 4), unit: '-', limit: round(BETZ_LIMIT, 4), sense: 'max', explain: 'Fraction of the kinetic power through the swept area that reaches the shaft. Modern rotors reach 0.45–0.50 below rated; the Betz limit is 0.593.' },
        { key: 'ct', label: 'Thrust coefficient Ct', value: round(bem.Ct, 4), unit: '-', explain: 'Thrust over ½ρAU². It drives the wake deficit behind the turbine: about 0.75–0.85 below rated, falling quickly once the blades pitch.' },
        { key: 'tsr', label: 'Tip-speed ratio', value: round(bem.tsr, 4), unit: '-', explain: 'Blade-tip speed over wind speed. The controller holds it near the optimum below rated; it falls above rated as rotor speed is capped.' },
        { key: 'cpMax', label: 'Maximum power coefficient of the rotor', value: round(opt.cpMax, 4), unit: '-', explain: `Best Cp over tip-speed ratio and pitch (at λ = ${fmt(opt.tsrOpt)} and pitch ${fmt(opt.pitchOptDeg, 2)}°). The published NREL 5-MW value is 0.482 at λ = 7.55 with the measured polars; the difference shows the uncertainty of the generic polars.` },
        { key: 'thrustRated', label: 'Peak steady thrust over the schedule', value: round(tb.thrustRated / 1e3, 5), unit: 'kN', explain: 'Largest steady thrust between cut-in and cut-out, reached just below rated wind speed. It is the reference thrust passed to the support-structure and foundation modules.' },
        { key: 'flapMoment', label: 'Blade-root flapwise moment (azimuth mean)', value: round(bem.flapMoment / 1e3, 5), unit: 'kN·m', explain: 'Out-of-plane aerodynamic bending moment at the blade root, averaged over a revolution. It governs spar-cap design and tip deflection towards the tower.' },
        { key: 'edgeMoment', label: 'Blade-root edgewise moment (aerodynamic)', value: round(bem.edgeMoment / 1e3, 5), unit: 'kN·m', explain: 'In-plane aerodynamic moment at the blade root that produces the torque. The gravity moment of the blade weight, reversing every revolution, is added in the fatigue module.' },
        { key: 'flap1P', label: 'Blade-root flapwise moment, 1P amplitude (unsteady)', value: round(un ? un.p1 / 1e3 : 0, 5), unit: 'kN·m', explain: 'Once-per-revolution variation of the flapwise root moment from wind shear, shaft tilt and tower shadow, with dynamic inflow and dynamic stall (0 when the unsteady analysis is switched off). It is the main deterministic fatigue load of the blade and pitch bearing.' },
        { key: 'thrustNP', label: 'Rotor thrust variation at blade-passing frequency', value: round(un ? un.thrustNP / 1e3 : 0, 5), unit: 'kN', explain: 'Amplitude of the thrust oscillation at the blade-passing frequency caused by tower shadow and shear: the 3P excitation of the tower.' },
        { key: 'tipSpeed', label: 'Blade-tip speed at rated rotor speed', value: round(tip.tipSpeed, 4), unit: 'm/s', limit: inp.tipLimit, sense: 'max', explain: 'Rated rotor speed times tip radius. Higher tip speed lowers torque and drivetrain mass but raises leading-edge erosion and noise.' },
        { key: 'ratedWind', label: 'Rated wind speed of the schedule', value: round(tb.ratedWind, 4), unit: 'm/s', explain: 'Wind speed at which rated electrical power is first reached. The NREL 5-MW reference quotes 11.4 m/s.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Spanwise blade loading', x: { label: 'Radius', unit: 'm' }, y: { label: 'Force per unit length', unit: 'kN/m' },
          series: [{ name: 'Normal (thrust-wise)', x: sp.map((q) => round(q.r, 4)), y: sp.map((q) => round(q.Np / 1e3, 4)) }, { name: 'Tangential (torque-wise)', x: sp.map((q) => round(q.r, 4)), y: sp.map((q) => round(q.Tp / 1e3, 4)) }],
          explain: 'Aerodynamic force per metre of blade, azimuth-averaged. The normal load grows towards the tip and drops at the very tip because of the Prandtl loss; the area under the tangential curve times radius is the torque.' },
        { type: 'line', title: 'Induction and angle of attack along the blade', x: { label: 'Radius', unit: 'm' }, y: { label: 'a (-) · α (deg/10)', unit: '-' },
          series: [{ name: 'Axial induction a', x: sp.map((q) => round(q.r, 4)), y: sp.map((q) => round(q.a, 4)) }, { name: 'Angle of attack / 10', x: sp.map((q) => round(q.r, 4)), y: sp.map((q) => round(q.alpha / DEG / 10, 4)) }, { name: 'Betz optimum a = 1/3', x: [sp[0].r, sp[sp.length - 1].r], y: [1 / 3, 1 / 3], dash: true }],
          explain: 'Axial induction factor and angle of attack (divided by 10) at each station. Induction near 1/3 over the outer blade means the rotor works close to the momentum-theory optimum; angles above about 1.2 on this scale are stalled.' },
        { type: 'line', title: 'Power coefficient versus tip-speed ratio', x: { label: 'Tip-speed ratio', unit: '-' }, y: { label: 'Cp', unit: '-' },
          series: [{ name: 'Fine pitch', x: tsrs.map((v) => round(v, 3)), y: cpFine.map((v) => round(v, 4)) }, { name: 'Pitch 4°', x: tsrs.map((v) => round(v, 3)), y: cp4.map((v) => round(v, 4)) },
            { name: 'Glauert ideal rotor', x: tsrs.map((v) => round(v, 3)), y: tsrs.map((t) => round(glauertIdealCp(t), 4)) }, { name: 'Operating point', x: [round(bem.tsr, 3)], y: [round(bem.Cp, 4)] }],
          explain: 'Cp curves of the rotor at fine pitch and with 4° of pitch, the ideal-rotor bound with wake rotation and the operating point. The peak of the fine-pitch curve defines the tip-speed ratio the torque controller tracks below rated.' },
        { type: 'line', title: 'Steady thrust and power over the operating range', x: { label: 'Hub-height wind speed', unit: 'm/s' }, y: { label: 'Thrust (kN) · power (kW/10)', unit: '-' },
          series: [{ name: 'Thrust (kN)', x: s.U, y: s.thrust.map((v) => round(v / 1e3, 4)) }, { name: 'Electrical power (kW/10)', x: s.U, y: s.power.map((v) => round(v / 1e4, 4)) }],
          explain: 'Steady thrust and electrical power from cut-in to cut-out. Thrust peaks at rated wind speed and then falls as the blades pitch: that peak sets the operating overturning moment on the foundation.' },
        { type: 'line', title: 'Rotor speed and pitch schedule', x: { label: 'Hub-height wind speed', unit: 'm/s' }, y: { label: 'Rotor speed (rpm) · pitch (deg)', unit: '-' },
          series: [{ name: 'Rotor speed (rpm)', x: s.U, y: s.omega.map((v) => round(v / RPM, 4)) }, { name: 'Pitch (deg)', x: s.U, y: s.pitchDeg.map((v) => round(v, 3)) }],
          explain: 'Controller set-points: rotor speed rises with wind speed to its rated value, then pitch increases to hold rated power. The speed range defines the 1P and 3P excitation bands used in the resonance check.' },
        ...(un ? [{ type: 'line', title: 'Blade-root flapwise moment over a revolution (unsteady)', x: { label: 'Azimuth of the blade (0° = up)', unit: 'deg' }, y: { label: 'Flapwise root moment', unit: 'kN·m' },
          series: (() => { const o = un.az.map((a, i) => [a, un.flap[i]]).sort((p, q) => p[0] - q[0]); return [{ name: 'Unsteady BEM', x: o.map((q) => round(q[0], 2)), y: o.map((q) => round(q[1] / 1e3, 4)) }, { name: 'Azimuth mean (steady BEM)', x: [0, 360], y: [bem.flapMoment / 1e3, bem.flapMoment / 1e3], dash: true }]; })(),
          explain: 'Flapwise root moment of one blade against azimuth with dynamic inflow and dynamic stall. The slow 1P swing comes from shear and tilt; the sharp dip near 180° is the tower passage.' }] : []),
        { type: 'polar', title: 'Single-blade thrust over a revolution', series: [{ name: 'Blade thrust (kN)', thetaDeg: bem.azimuthal.psi.map((p) => round(p / DEG, 2)), r: bem.azimuthal.bladeThrust.map((v) => round(v / 1e3, 4)) }],
          explain: 'Thrust on one blade against azimuth (0° = blade up). Shear, tilt and yaw make it vary once per revolution: that 1P variation is a main source of blade and shaft fatigue.' },
      ],
      tables: [
        { title: 'Blade stations', columns: [{ key: 'r', label: 'Radius', unit: 'm' }, { key: 'chord', label: 'Chord', unit: 'm' }, { key: 'twist', label: 'Twist', unit: 'deg' }, { key: 'tc', label: 't/c', unit: '-' }, { key: 'a', label: 'a', unit: '-' }, { key: 'ap', label: "a'", unit: '-' }, { key: 'alpha', label: 'α', unit: 'deg' }, { key: 'cl', label: 'Cl', unit: '-' }, { key: 'Np', label: 'Normal load', unit: 'kN/m' }, { key: 'Tp', label: 'Tangential load', unit: 'kN/m' }, { key: 'Re', label: 'Reynolds number', unit: '10⁶' }],
          rows: sp.map((q, i) => ({ r: round(q.r, 4), chord: round(rotor.stations[i].chord, 4), twist: round(rotor.stations[i].twistDeg, 3), tc: round(rotor.stations[i].thickness ?? 1, 3), a: round(q.a, 3), ap: round(q.ap, 4), alpha: round(q.alpha / DEG, 3), cl: round(q.cl, 3), Np: round(q.Np / 1e3, 4), Tp: round(q.Tp / 1e3, 4), Re: round(q.Re / 1e6, 3) })),
          explain: 'One row per blade station: planform, induction factors, angle of attack, lift coefficient and loads. Use it to see where the blade works away from its design incidence or outside the Reynolds range of the polars.' },
        { title: 'Operating schedule', columns: [{ key: 'U', label: 'Wind speed', unit: 'm/s' }, { key: 'region', label: 'Region' }, { key: 'rpm', label: 'Rotor speed', unit: 'rpm' }, { key: 'pitch', label: 'Pitch', unit: 'deg' }, { key: 'P', label: 'Electrical power', unit: 'kW' }, { key: 'T', label: 'Thrust', unit: 'kN' }, { key: 'cp', label: 'Cp', unit: '-' }, { key: 'ct', label: 'Ct', unit: '-' }],
          rows: s.U.map((u, i) => ({ U: u, region: s.region[i], rpm: round(s.omega[i] / RPM, 4), pitch: round(s.pitchDeg[i], 3), P: round(s.power[i] / 1e3, 5), T: round(s.thrust[i] / 1e3, 5), cp: round(s.cp[i], 4), ct: round(s.ct[i], 4) })),
          explain: 'Steady power and thrust curve of the rotor model with the control region at each wind speed. The power and Ct columns can be pasted into the energy-yield and wake modules as a curve table.' },
      ],
      notes: [
        `Operating point: ${fmt(inp.U)} m/s, ${fmt(omega / RPM, 4)} rpm, pitch ${fmt(pitch, 3)}° (${inp.opMode === 'manual' ? 'prescribed' : 'from the schedule'}); BEM ${bem.converged ? 'converged at every station' : 'did NOT converge at every station — treat the result with caution'}.`,
        `Governing check: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        `Region-2 torque gain k = ${fmt(k2.kLss, 4)} N·m/(rad/s)² on the rotor shaft (${fmt(k2.kHssRpm, 4)} N·m/rpm² on the generator shaft).`,
        `Rotor: ${rotor.name}, diameter ${fmt(2 * R, 4)} m, ${rotor.blades} blades, RNA mass ${fmt(rotor.rnaMass / 1e3, 4)} t. Generic polars: the maximum Cp of ${fmt(opt.cpMax)} is ${fmt(100 * (opt.cpMax / 0.482 - 1), 2)} % from the published 0.482 of the reference turbine.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `Rotor thrust is ${fmt(bem.thrust / 1e3)} kN and Cp ${fmt(bem.Cp)} at ${fmt(inp.U)} m/s; the peak steady thrust over the schedule is ${fmt(tb.thrustRated / 1e3)} kN at ${fmt(tb.ratedWind)} m/s.`,
          gov.uc > 1 ? `Bring '${gov.name}' within its limit (utilisation ${fmt(gov.uc)}): reduce rotor speed or change pitch, then re-run.` : `Carry the peak thrust of ${fmt(tb.thrustRated / 1e3)} kN to the support-structure and foundation modules; ${SUPPLIER}`,
          'IEC 61400-1:2019 §7.4 (DLC 1.x power production) with steady BEM per Moriarty & Hansen (2005); peak-shaving of thrust near rated is current practice for large rotors.',
          `Blades of this rotor: about ${fmt(gfrp)} t of glass-fibre composite ≈ ${fmt(gfrp * CO2_GFRP)} t CO₂e (screening factor ${CO2_GFRP} t CO₂e/t) with no established recycling route at scale — blade life extension and leading-edge protection are the main sustainability levers.`),
        adv('economic', 'info', `Each 1 % of power coefficient lost below rated (soiling, erosion, pitch misalignment) costs about ${fmt(aepLossPerPct)} MWh a year per turbine at a 45 % capacity factor.`,
          'Check the supplier power curve against this estimate at contract stage and include a leading-edge inspection in the O&M plan; a thrust limit (peak shaving) can be traded against about 1–2 % energy if the foundation governs.',
          'LCOE convention (IEA / IRENA): annual energy is the denominator, so a percentage of energy is worth about the same percentage of lifetime cost; power-curve warranty practice per IEC 61400-12-1:2017.',
          energyCo2Text(aepLossPerPct)),
      ],
      publish: { 'owt.thrustRated': tb.thrustRated, 'owt.rna': rotor.rnaMass, 'owt.ratedPower': rotor.ratedPower, 'owt.rotorDiameter': 2 * R, 'owt.hubHeight': inp.hubHeight },
      quantities: { steel_t: steel, polymer_t: round(gfrp, 4) },
    };
  },
  verify: [
    { name: 'NREL 5-MW reference: mechanical power 5.297 MW at 11.4 m/s, 12.1 rpm, zero pitch (generic polars within 5 %)',
      source: 'Jonkman, Butterfield, Musial & Scott (2009) NREL/TP-500-38060 Table 1-1: rated electrical 5 MW at 94.4 % efficiency → 5.297 MW mechanical',
      // 5.0 MW / 0.944 = 5.2966 MW; tip-speed ratio = (12.1·π/30 rad/s)·R/11.4 with the engine's coned radius 63·cos(2.5°) = 62.94 m → 6.996
      inputs: { opMode: 'manual', U: 11.4, rpm: 12.1, pitchDeg: 0, shearExp: 0, yawDeg: 0 }, expect: [{ kpi: 'aeroPower', value: 5.297, tol: 0.05 }, { kpi: 'tsr', value: 6.996, tol: 0.005 }] },
    { name: 'NREL 5-MW reference: peak power coefficient 0.482 (generic polars within 5 %)',
      source: 'Jonkman et al. (2009) NREL/TP-500-38060 §7.2: Cp,max = 0.482 at tip-speed ratio 7.55',
      inputs: {}, expect: [{ kpi: 'cpMax', value: 0.482, tol: 0.05 }] },
  ],
};
