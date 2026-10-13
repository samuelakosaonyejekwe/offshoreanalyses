// Design load cases of a bottom-fixed offshore wind turbine on the reduced-order coupled simulator: a subset of
// IEC 61400-3-1 Table 2 run as time-domain simulations and tabulated as ultimate (characteristic and design) and
// fatigue (damage-equivalent) loads at the tower base and the mudline, with the section stress and Miner-sum checks
// that follow from them. Engineering units (mm, MN·m, MPa) are converted to SI at the top of compute.
import { runDesignLoadCases, LOAD_SAFETY_FACTORS } from '../../engine/aero/dlc.js';
import { buildSupportStructure } from '../../engine/aero/coupled.js';
import { ewm } from '../../engine/aero/turbulence.js';
import { snCurve, cyclesToFailure, SN_CURVE_NAMES } from '../../engine/fatigue/sn.js';
import { turbineInputs, getTurbine, supportInputs, supportSegments, foundationInputs, foundationOf, segmentsMass, sectionAt, site, adv, check, worst, sev, fmt, round, run, need, G, DEG,
  LIMIT_COUPLED, LIMIT_POLAR, LIMIT_CODE, LIMIT_CONTROL, CONFIRM, SUPPLIER, CO2_STEEL, steelCo2Text, costInput, money } from './_common.js';

const curves = SN_CURVE_NAMES.filter((c) => ['B1', 'B2', 'C', 'C1', 'C2', 'D', 'E', 'F', 'F1', 'F3', 'G'].includes(c));
// case id → [input key, on by default, short description, IEC type]
const CASES = [
  ['1.2', 'run12', true, 'Power production, normal turbulence, normal sea (fatigue)', 'F'],
  ['1.3', 'run13', true, 'Power production, extreme turbulence, normal sea', 'N'],
  ['1.4', 'run14', false, 'Power production, extreme coherent gust with direction change', 'N'],
  ['1.6', 'run16', false, 'Power production, normal turbulence, severe sea', 'N'],
  ['2.3', 'run23', false, 'Power production, extreme operating gust with grid loss', 'A'],
  ['5.1', 'run51', false, 'Emergency stop', 'N'],
  ['6.1', 'run61', true, 'Parked, 50-year wind and extreme sea', 'N'],
  ['6.2', 'run62', false, 'Parked, 50-year wind, grid loss (yaw errors to ±180°)', 'A'],
  ['7.1', 'run71', false, 'Parked with fault, 1-year wind and sea', 'A'],
];
const sectionW = (D, t) => Math.PI * (D ** 4 - (D - 2 * t) ** 4) / (32 * D);
const sectionA = (D, t) => Math.PI * (D - t) * t;

export default {
  id: 's16-design-load-cases', suite: 16, title: 'Design load cases (IEC 61400-3-1 subset)',
  summary: 'Runs a subset of the IEC design load cases in the time domain and tabulates the ultimate and fatigue loads at the tower base and the mudline, with the partial factors, the governing case and the resulting section utilisations.',
  covers: [1, 6, 10], fidelity: 'reduced-order',
  standards: ['IEC 61400-3-1:2019 Table 2 (design load cases) and §7', 'IEC 61400-1:2019 §6.3 (wind conditions), §7.6 (ultimate and fatigue limit states), Table 3 (partial safety factors for loads)', 'DNV-ST-0437 (2021) (loads and site conditions for wind turbines)', 'DNV-ST-0126 (2021) §4 (support structures: material and design fatigue factors)', 'DNV-RP-C203 (2021) (S–N curves, Miner sum)'],
  models: ['Fully coupled aero-hydro-servo-elastic simulation (reduced order)', 'Normal and extreme turbulence models (NTM, ETM)', 'Extreme operating gust and extreme coherent gust with direction change (EOG, ECD)', 'Extreme wind speed model, turbulent (EWM)', 'Irregular linear waves with Morison loading on the monopile', 'Variable-speed torque and gain-scheduled PI pitch control with emergency-stop and grid-loss logic', 'Joint wind–wave bins with Weibull probabilities', 'Rainflow counting, damage-equivalent load and Miner sum', 'Euler–Bernoulli beam finite elements and modal reduction of the support structure'],
  limits: `${LIMIT_COUPLED} A deliberately simplified subset of the load-case table: one wave direction aligned with the wind, no wind–wave misalignment or water-level variation, direction changes entering as the cosine of the misalignment, no statistical extrapolation (DLC 1.1), no fault cases other than grid loss, and by default one short seed per condition where the standard requires at least six 10-minute realisations. The parked-rotor thrust coefficients are inputs. Section checks use bending plus self-weight at the mudline and tower base only; shell buckling and the embedded pile are in the strength and foundation modules. Results are for concept screening of support-structure loads, not for certification. ${LIMIT_POLAR} ${LIMIT_CONTROL} ${LIMIT_CODE}`,
  publishes: ['owt.mudlineMoment', 'owt.mudlineShear', 'owt.f1', 'owt.towerBaseMoment'],
  searchTerms: 'IEC 61400-3-1 design load cases offshore wind turbine ultimate fatigue loads monopile tower base',
  inputs: [
    ...turbineInputs(),
    // Walls a little heavier than the shared reference support structure, so that the built-in example also carries
    // the fatigue case with the design fatigue factors (the reference walls are sized on the ultimate cases).
    ...supportInputs().map((i) => (i.key === 'monoT' ? { ...i, default: 90 } : i.key === 'towerBaseT' ? { ...i, default: 40 } : i)),
    ...foundationInputs(),
    { key: 'zeta', label: 'Structural damping ratio', unit: '-', type: 'number', default: 0.01, min: 0.002, max: 0.05, group: 'Support structure', role: 'material', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Steel, soil and hydrodynamic damping lumped per mode; aerodynamic damping comes from the model.' },
    { key: 'Vref', label: 'Reference wind speed Vref (turbine class)', unit: 'm/s', type: 'number', default: 50, min: 30, max: 60, group: 'Wind conditions', role: 'bc', help: `Class I 50 m/s, II 42.5 m/s, III 37.5 m/s: the 50-year 10-minute mean at hub height. ${CONFIRM}` },
    { key: 'Iref', label: 'Reference turbulence intensity Iref', unit: '-', type: 'number', default: 0.12, min: 0.05, max: 0.2, group: 'Wind conditions', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 }, help: `IEC 61400-1 Table 1: A 0.16, B 0.14, C 0.12. ${CONFIRM}` },
    site('A'), site('k'),
    { key: 'seaTable', label: 'Normal sea state by wind speed (empty = generic fetch-limited wind sea)', type: 'table', group: 'Sea state', role: 'bc', default: [],
      columns: [{ key: 'U', label: 'Hub wind speed', unit: 'm/s', type: 'number' }, { key: 'Hs', label: 'Hs', unit: 'm', type: 'number' }, { key: 'Tp', label: 'Tp', unit: 's', type: 'number' }],
      help: 'Expected Hs and Tp conditional on the hub wind speed, from the site scatter data (at least two rows).' },
    site('Hs100', { label: '50-year significant wave height (severe and extreme sea state)' }), site('Tp100'), site('Hs1'), site('Tp1'),
    { key: 'ctParked', label: 'Thrust coefficient of the parked, feathered rotor', unit: '-', type: 'number', default: 0.014, min: 0.002, max: 0.1, group: 'Parked rotor', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.3 }, help: `Referred to the swept area, small yaw error, including nacelle drag: about 0.01–0.02. ${SUPPLIER}` },
    { key: 'ctBroad', label: 'Thrust coefficient with the wind broadside to the feathered blades', unit: '-', type: 'number', default: 0.045, min: 0.005, max: 0.2, group: 'Parked rotor', role: 'bc', help: 'Used for DLC 6.2 between the aligned value and this one in proportion to |sin(yaw error)|: blade planform area × drag coefficient over swept area, about 0.03–0.06.' },
    { key: 'ctFault', label: 'Thrust coefficient of the parked rotor with a fault (DLC 7.1)', unit: '-', type: 'number', default: 0.03, min: 0.002, max: 0.3, group: 'Parked rotor', role: 'bc', help: 'For example one blade seized away from feather.' },
    ...CASES.map(([id, key, on, text]) => ({ key, label: `DLC ${id}: ${text}`, type: 'bool', default: on, group: 'Load cases to run', role: 'control', help: on ? 'Part of the default set.' : 'Optional: adds simulations and run time.' })),
    { key: 'duration', label: 'Simulated time per stochastic run', unit: 's', type: 'number', default: 60, min: 30, max: 660, group: 'Simulation', role: 'control', help: 'After the transient. 600 s is the design standard: short runs under-sample the extremes and the largest fatigue cycles.' },
    { key: 'nSeeds', label: 'Seeds per condition', type: 'int', default: 1, min: 1, max: 6, group: 'Simulation', role: 'control', help: 'IEC 61400-1 asks for at least six 10-minute seeds per stochastic condition.' },
    { key: 'binWidth', label: 'Wind-speed bin width for DLC 1.2', unit: 'm/s', type: 'number', default: 5.5, min: 1, max: 8, group: 'Simulation', role: 'control', help: 'IEC uses 2 m/s bins; wider bins are faster and coarser.' },
    { key: 'seed', label: 'Random seed', type: 'int', default: 1, min: 1, max: 99999, group: 'Simulation', role: 'ic' },
    { key: 'dt', label: 'Time step', unit: 's', type: 'number', default: 0.1, min: 0.0125, max: 0.2, group: 'Simulation', role: 'control' },
    { key: 'nElem', label: 'Beam elements of the support structure', type: 'int', default: 20, min: 6, max: 160, group: 'Simulation', role: 'control' },
    { key: 'designLife', label: 'Design life', unit: 'yr', type: 'number', default: 25, min: 5, max: 50, group: 'Acceptance', role: 'criteria', bind: 'life.designYears' },
    { key: 'avail', label: 'Time in power production within the operating wind range', unit: '%', type: 'number', default: 95, min: 50, max: 100, group: 'Acceptance', role: 'ic' },
    { key: 'gammaM', label: 'Material factor for yield', unit: '-', type: 'number', default: 1.1, min: 1, max: 1.3, group: 'Acceptance', role: 'criteria', help: `DNV-ST-0126: 1.10. ${CONFIRM}` },
    { key: 'snMono', label: 'S–N curve of the monopile girth weld (seawater with cathodic protection)', type: 'select', default: 'D', group: 'Acceptance', role: 'material', options: curves.map((c) => ({ value: c, label: `DNV-RP-C203 curve ${c}` })) },
    { key: 'snTower', label: 'S–N curve of the tower-base weld (in air)', type: 'select', default: 'D', group: 'Acceptance', role: 'material', options: curves.map((c) => ({ value: c, label: `DNV-RP-C203 curve ${c}` })) },
    { key: 'scf', label: 'Stress concentration factor at the welds', unit: '-', type: 'number', default: 1.15, min: 1, max: 3, group: 'Acceptance', role: 'material', uncertain: { dist: 'lognormal', cov: 0.1 } },
    { key: 'dffMono', label: 'Design fatigue factor, monopile', unit: '-', type: 'number', default: 3, min: 1, max: 10, group: 'Acceptance', role: 'criteria', help: `DNV-ST-0126: 3 where inspection and repair are not planned. ${CONFIRM}` },
    { key: 'dffTower', label: 'Design fatigue factor, tower', unit: '-', type: 'number', default: 2, min: 1, max: 10, group: 'Acceptance', role: 'criteria', help: CONFIRM },
    costInput('steelRate', 'Fabricated and installed support-structure steel', 'USD/t', 4000, 50000, 'Used to express a change of wall thickness as money.'),
  ],
  mesh: { key: 'nElem', label: 'Beam elements', levels: [8, 16, 32, 64], qoi: 'f1' },
  meshTime: { key: 'dt', label: 'Time step (s)', levels: [0.2, 0.1, 0.05], qoi: 'delMud', isSize: true },
  compute(inp, ctx) {
    const tb = getTurbine(inp), rotor = tb.rotor, R = rotor.tipRadius;
    need(inp.hubHeight > R + 5, `The hub height (${inp.hubHeight} m) leaves no clearance under the blade tip (rotor radius ${fmt(R)} m).`);
    const st = ctx.material(inp.steel);
    need(st && st.fy > 0 && st.E > 0, `Steel '${inp.steel}' is not in the materials library.`);
    const ids = CASES.filter(([, key]) => inp[key]).map(([id]) => id);
    need(ids.length > 0, 'Select at least one design load case.');
    const nRuns = ids.reduce((s, id) => s + ({ '1.2': Math.max(1, Math.round((rotor.cutOut - rotor.cutIn) / inp.binWidth)) * inp.nSeeds, '1.3': 4 * inp.nSeeds, '1.4': 6, '1.6': 4 * inp.nSeeds, '2.3': 12, '5.1': 4 * inp.nSeeds, '6.1': 3 * inp.nSeeds, '6.2': 6 * inp.nSeeds, '7.1': inp.nSeeds })[id], 0);
    need(nRuns * (inp.duration + 20) / inp.dt <= 3e6, `The selection needs ${nRuns} simulations of ${fmt(inp.duration + 20)} s at ${inp.dt} s steps: reduce the cases, seeds or duration, or raise the time step.`);
    const segs = supportSegments(inp, { zHub: inp.hubHeight, hubOffset: 2.4 * tb.scale, E: st.E, rhoSteel: st.rho ?? 7850 });
    const structure = run(() => buildSupportStructure({ segments: segs, rnaMass: rotor.rnaMass, zHub: inp.hubHeight, zTowerBase: inp.transZ, foundation: foundationOf(inp), nElem: Math.round(inp.nElem), nModes: 2, zeta: inp.zeta }), 'Support-structure model');
    const f1 = structure.modes[0].frequency;
    const model = { surface: tb.surface, turbine: { inertia: rotor.drivetrainInertia, hubHeight: inp.hubHeight }, controller: tb.controller, structure, dt: inp.dt, shearExp: 0.14 };
    const table = (inp.seaTable || []).filter((r) => r.U > 0 && r.Hs >= 0 && r.Tp > 0).sort((a, b) => a.U - b.U);
    const siteRec = { hubHeight: inp.hubHeight, rotorDiameter: 2 * R, waterDepth: inp.depth, Vref: inp.Vref, Iref: inp.Iref, ratedWind: tb.ratedWind, cutIn: rotor.cutIn, cutOut: rotor.cutOut,
      weibullA: inp.weibA, weibullK: inp.weibK, Hs50: inp.Hs100, Tp50: inp.Tp100, Hs1: inp.Hs1, Tp1: inp.Tp1, ...(table.length >= 2 ? { seaState: table } : {}) };
    // sections and weld S–N curves (stress in Pa from the section moment in N·m)
    const secM = sectionAt(segs, -inp.depth), secT = sectionAt(segs, inp.transZ + 1e-6);
    const Wm = sectionW(secM.D, secM.t), Wt = sectionW(secT.D, secT.t);
    const cvM = snCurve(inp.snMono, 'seawater-cp'), cvT = snCurve(inp.snTower, 'air');
    const snDamage = ({ channel, ranges, counts }) => {
      const mud = channel === 'mudlineMoment', W = mud ? Wm : Wt, curve = mud ? cvM : cvT, t = mud ? secM.t : secT.t;
      let d = 0;
      for (let i = 0; i < ranges.length; i++) { const sr = ranges[i] / W; if (sr > 1e4) d += counts[i] / cyclesToFailure({ curve, stressRange: sr, thickness: t, scf: inp.scf }); }
      return d;
    };
    const ctYaw = (yaw) => inp.ctParked + (inp.ctBroad - inp.ctParked) * Math.abs(Math.sin(yaw * DEG));
    ctx.progress?.(0.05, `Running ${nRuns} simulations`);
    const res = run(() => runDesignLoadCases({ model, site: siteRec, cases: ids, nSeeds: inp.nSeeds, duration: inp.duration, tTransient: 20, seed: inp.seed, binWidth: inp.binWidth,
      lifetimeYears: inp.designLife, availability: inp.avail / 100, wohlerM: 4, nEq: 1e7, snDamage, parkedCt: ctYaw, parkedFaultCt: inp.ctFault }), 'Design load cases');
    const ult = res.cases, fat = res.fatigue;
    const e = ewm({ Vref: inp.Vref, z: inp.hubHeight, hubHeight: inp.hubHeight });
    // governing ultimate case per quantity (design values)
    const top = (k) => ult.reduce((a, c) => (c.design[k] > (a ? a.design[k] : -1) ? c : a), null);
    const gMud = top('mudlineMoment'), gTow = top('towerBaseMoment'), gShear = top('mudlineShear');
    const Md = gMud ? gMud.design.mudlineMoment : 0, Mt = gTow ? gTow.design.towerBaseMoment : 0, Vd = gShear ? gShear.design.mudlineShear : 0;
    // self-weight above each section (favourable effects ignored: load factor 1.0 on gravity for the stress sum)
    const mTower = segmentsMass([segs[1]]), mAll = segmentsMass(segs), Nm = (mAll + rotor.rnaMass) * G, Nt = (mTower + rotor.rnaMass) * G;
    const sigM = Md / Wm + Nm / sectionA(secM.D, secM.t), sigT = Mt / Wt + Nt / sectionA(secT.D, secT.t), fyd = st.fy / inp.gammaM;
    const dMud = fat ? fat.damage.mudlineMoment : 0, dTow = fat ? fat.damage.towerBaseMoment : 0;
    const opFrac = fat ? fat.bins.reduce((s, b) => s + b.prob, 0) : Math.exp(-((rotor.cutIn / inp.weibA) ** inp.weibK)) - Math.exp(-((rotor.cutOut / inp.weibA) ** inp.weibK));

    const checks = [];
    if (gMud) checks.push(check(`Yield at the mudline (DLC ${gMud.dlc})`, sigM / 1e6, fyd / 1e6, 'MPa', `IEC 61400-3-1:2019 §7.6.2 with γf = ${gMud.psf}; DNV-ST-0126 (2021) §4.4 (γM = ${inp.gammaM})`,
      'Extreme-fibre stress from the governing design mudline moment plus self-weight against design yield. Exceeding it means first yield of the monopile shell under the factored extreme load.'));
    if (gTow) checks.push(check(`Yield at the tower base (DLC ${gTow.dlc})`, sigT / 1e6, fyd / 1e6, 'MPa', `IEC 61400-3-1:2019 §7.6.2 with γf = ${gTow.psf}; DNV-ST-0126 (2021) §4.4`,
      'Extreme-fibre stress at the tower-base section, where the wall steps down and the flange sits, against design yield.'));
    if (fat) {
      checks.push(check('Fatigue damage at the mudline weld over the design life (DLC 1.2)', inp.dffMono * dMud, 1, '-', `DNV-RP-C203 (2021) §2.2 curve ${inp.snMono}, seawater with cathodic protection; DFF = ${inp.dffMono}`,
        'Miner sum of the rainflow-counted stress ranges of the power-production bins, scaled to the design life, times the design fatigue factor. Above 1 a fatigue crack is expected at the girth weld within the life.'));
      checks.push(check('Fatigue damage at the tower base over the design life (DLC 1.2)', inp.dffTower * dTow, 1, '-', `DNV-RP-C203 (2021) curve ${inp.snTower} in air; DFF = ${inp.dffTower}`,
        'The same Miner sum at the tower-base weld. Parked and transient cycles are not included, so the true damage is higher.'));
    }
    const gov = worst(checks);
    const mn = (v) => round(v / 1e6, 5);
    const byId = Object.fromEntries(ult.map((c) => [c.dlc, c]));
    const meta = Object.fromEntries(CASES.map(([id, , , text, type]) => [id, { text, type }]));
    const kpis = [
      { key: 'f1', label: 'First fore-aft natural frequency', value: round(f1, 5), unit: 'Hz', explain: `First bending mode of the support structure with the rotor–nacelle mass. It must stay clear of the rotor frequency 1P (${fmt(rotor.minRotorSpeed / (2 * Math.PI))}–${fmt(rotor.ratedRotorSpeed / (2 * Math.PI))} Hz) and the blade-passing frequency 3P.` },
      { key: 'v50', label: '50-year 10-minute mean wind speed at hub height', value: round(e.V50, 5), unit: 'm/s', explain: 'Extreme wind speed model (turbulent): V50 at hub height equals Vref; the parked cases 6.1 and 6.2 are simulated at this mean with 11 % turbulence.' },
      { key: 'opFrac', label: 'Time between cut-in and cut-out', value: round(opFrac * 100, 6), unit: '%', explain: 'Weibull probability of the operating wind range: the share of the life in which the fatigue of DLC 1.2 accumulates (before availability).' },
      { key: 'nRuns', label: 'Simulations run', value: nRuns, unit: '-', explain: `Time-domain runs behind the table, each ${fmt(inp.duration)} s after a 20 s transient. The standard asks for six 10-minute seeds per stochastic condition.` },
    ];
    if (gMud) kpis.push(
      { key: 'mUlt', label: 'Design mudline moment (governing case)', value: mn(Md), unit: 'MN·m', explain: `Largest characteristic mudline overturning moment times its partial factor, from DLC ${gMud.dlc} (${gMud.governing.mudlineMoment}). This is the moment the monopile and its foundation are sized for.` },
      { key: 'vUlt', label: 'Design mudline shear (governing case)', value: mn(Vd), unit: 'MN', explain: `Largest factored base shear, from DLC ${gShear.dlc}. With the moment it is the load pair for the lateral pile analysis.` },
      { key: 'mUltTower', label: 'Design tower-base moment (governing case)', value: mn(Mt), unit: 'MN·m', explain: `Largest factored fore-aft moment at the tower base / transition, from DLC ${gTow.dlc}: the load on the flange or grouted connection.` },
      { key: 'psfGov', label: 'Partial safety factor of the governing case', value: gMud.psf, unit: '-', explain: 'IEC 61400-1 Table 3: 1.35 for normal (N) and 1.1 for abnormal (A) design situations, applied to the characteristic extreme.' },
      { key: 'sigMud', label: 'Design stress at the mudline', value: round(sigM / 1e6, 4), unit: 'MPa', limit: round(fyd / 1e6, 5), sense: 'max', explain: `Bending stress M/W (section modulus ${fmt(Wm, 4)} m³) plus self-weight stress, against design yield fy/γM.` },
      { key: 'sigTower', label: 'Design stress at the tower base', value: round(sigT / 1e6, 4), unit: 'MPa', limit: round(fyd / 1e6, 5), sense: 'max', explain: `Bending plus self-weight stress at the tower base (section modulus ${fmt(Wt, 4)} m³) against design yield.` });
    for (const c of ult) kpis.push({ key: `m${c.dlc.replace('.', '')}`, label: `DLC ${c.dlc}: characteristic mudline moment`, value: mn(c.characteristic.mudlineMoment), unit: 'MN·m', explain: `${meta[c.dlc].text}. Largest absolute mudline moment over ${c.runs.length} run(s) before the partial factor of ${c.psf}; governing run: ${c.governing.mudlineMoment}.` });
    if (byId['6.1']) kpis.push({ key: 'psf61', label: 'DLC 6.1: design / characteristic moment', value: round(byId['6.1'].design.mudlineMoment / byId['6.1'].characteristic.mudlineMoment, 6), unit: '-', explain: 'Ratio of the design to the characteristic mudline moment of the parked 50-year case: the normal-situation load factor.' });
    if (fat) kpis.push(
      { key: 'delMud', label: 'Lifetime damage-equivalent mudline moment (m = 4, 10⁷ cycles)', value: mn(fat.del.mudlineMoment), unit: 'MN·m', explain: 'Constant moment range that causes the lifetime damage of power production in ten million cycles: the single number for comparing designs, sites and controllers.' },
      { key: 'delTower', label: 'Lifetime damage-equivalent tower-base moment', value: mn(fat.del.towerBaseMoment), unit: 'MN·m', explain: 'The same measure at the tower base.' },
      { key: 'dmgMud', label: 'Design fatigue damage at the mudline weld', value: round(inp.dffMono * dMud, 4), unit: '-', limit: 1, sense: 'max', explain: `Miner sum over ${fmt(inp.designLife)} years times the design fatigue factor ${inp.dffMono}. Short single-seed runs under-sample the largest cycles: treat it as a lower-bound screening value.` },
      { key: 'dmgTower', label: 'Design fatigue damage at the tower base', value: round(inp.dffTower * dTow, 4), unit: '-', limit: 1, sense: 'max', explain: `Miner sum times the design fatigue factor ${inp.dffTower} at the tower-base weld.` });

    const steel = (mAll) / 1e3, dT = 0.001, dSteel = Math.PI * (secM.D - secM.t) * dT * (inp.depth + inp.transZ) * (st.rho ?? 7850) / 1e3;
    const ucY = gMud ? sigM / fyd : 0, ucF = fat ? inp.dffMono * dMud : 0;
    const tNeed = Math.max(ucY, ucF ** (1 / 4)) * inp.monoT;   // stress ∝ 1/t; damage ∝ stress^m with m ≈ 4 between the S–N slopes

    return {
      kpis, checks,
      plots: [
        ...(ult.length ? [{ type: 'bar', title: 'Mudline moment by design load case', categories: ult.map((c) => `DLC ${c.dlc}`), x: { label: 'Design load case', unit: '' }, y: { label: 'Mudline moment', unit: 'MN·m' },
          series: [{ name: 'Characteristic', y: ult.map((c) => mn(c.characteristic.mudlineMoment)) }, { name: 'Design (× γf)', y: ult.map((c) => mn(c.design.mudlineMoment)) }],
          explain: 'Characteristic and factored mudline moment of every ultimate case. The tallest design bar governs; abnormal cases carry the smaller factor, so a case can lead on characteristic load and still not govern.' },
        { type: 'bar', title: 'Tower-base moment and mudline shear by design load case', categories: ult.map((c) => `DLC ${c.dlc}`), x: { label: 'Design load case', unit: '' }, y: { label: 'Design value (MN·m, MN)', unit: '' },
          series: [{ name: 'Tower-base moment, design (MN·m)', y: ult.map((c) => mn(c.design.towerBaseMoment)) }, { name: 'Mudline shear, design × 10 (MN)', y: ult.map((c) => mn(10 * c.design.mudlineShear)) }],
          explain: 'Wind-driven cases load the tower base almost as much as the mudline; wave-driven cases show up in the shear and add moment mainly below the water line.' }] : []),
        ...(fat ? [{ type: 'bar', title: 'Damage-equivalent mudline moment by wind bin (DLC 1.2)', categories: fat.bins.map((b) => `${fmt(b.U)} m/s`), x: { label: 'Wind-speed bin', unit: '' }, y: { label: '1-Hz damage-equivalent moment / share', unit: '' },
          series: [{ name: '1-Hz DEL (MN·m)', y: fat.bins.map((b) => mn(b.del.mudlineMoment)) }, { name: 'Probability (%)', y: fat.bins.map((b) => round(b.prob * 100, 4)) }],
          explain: 'Damage-equivalent moment of each wind bin with its probability. The bin that combines a high value with a high probability contributes most to the lifetime damage, usually just above rated wind speed.' }] : []),
      ],
      tables: [
        { title: 'Design load case table', columns: [{ key: 'dlc', label: 'DLC' }, { key: 'desc', label: 'Design situation' }, { key: 'type', label: 'Type' }, { key: 'psf', label: 'γf' }, { key: 'mChar', label: 'Mudline moment, char.', unit: 'MN·m' }, { key: 'mDes', label: 'Mudline moment, design', unit: 'MN·m' }, { key: 'tChar', label: 'Tower-base moment, char.', unit: 'MN·m' }, { key: 'tDes', label: 'Tower-base moment, design', unit: 'MN·m' }, { key: 'vDes', label: 'Mudline shear, design', unit: 'MN' }, { key: 'gov', label: 'Governing run' }],
          rows: [...ult.map((c) => ({ dlc: c.dlc, desc: c.description, type: c.type === 'N' ? 'Ultimate, normal' : 'Ultimate, abnormal', psf: c.psf, mChar: mn(c.characteristic.mudlineMoment), mDes: mn(c.design.mudlineMoment), tChar: mn(c.characteristic.towerBaseMoment), tDes: mn(c.design.towerBaseMoment), vDes: mn(c.design.mudlineShear), gov: c.governing.mudlineMoment })),
            ...(fat ? [{ dlc: '1.2', desc: 'Power production, NTM, NSS', type: 'Fatigue', psf: LOAD_SAFETY_FACTORS.F, mChar: mn(fat.del.mudlineMoment), mDes: mn(fat.del.mudlineMoment), tChar: mn(fat.del.towerBaseMoment), tDes: mn(fat.del.towerBaseMoment), vDes: 0, gov: 'damage-equivalent range, m = 4, 10⁷ cycles' }] : [])],
          explain: 'One row per load case: the characteristic extreme over its runs, the design value with the partial factor of its type, and the run that produced it. The fatigue row lists damage-equivalent moment ranges instead of extremes.' },
        ...(fat ? [{ title: 'Wind–wave bins of DLC 1.2', columns: [{ key: 'U', label: 'Wind speed', unit: 'm/s' }, { key: 'p', label: 'Probability', unit: '%' }, { key: 'Hs', label: 'Hs', unit: 'm' }, { key: 'Tp', label: 'Tp', unit: 's' }, { key: 'delM', label: '1-Hz DEL mudline', unit: 'MN·m' }, { key: 'delT', label: '1-Hz DEL tower base', unit: 'MN·m' }],
          rows: fat.bins.map((b) => ({ U: round(b.U), p: round(b.prob * 100), Hs: round(b.Hs), Tp: round(b.Tp), delM: mn(b.del.mudlineMoment), delT: mn(b.del.towerBaseMoment) })),
          explain: 'The joint wind and sea-state description behind the fatigue case, with the damage-equivalent moments of each bin. Replace the generic sea states with the site scatter table.' }] : []),
      ],
      notes: [...res.notes, `${ids.length} load case(s), ${nRuns} simulation(s) of ${fmt(inp.duration)} s at ${inp.dt} s steps; support-structure frequency ${fmt(f1)} Hz.`,
        `Not run: ${CASES.filter(([id]) => !ids.includes(id)).map(([id]) => id).join(', ') || 'none'}. Switch the remaining cases on before quoting a governing load.`],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name}: ${fmt(gov.demand)} against ${fmt(gov.capacity)}${gov.unit === '-' ? '' : ` ${gov.unit}`} (utilisation ${fmt(gov.uc)}).${gMud ? ` Governing ultimate case DLC ${gMud.dlc}: design mudline moment ${fmt(Md / 1e6)} MN·m, shear ${fmt(Vd / 1e6)} MN.` : ''}${fat ? ` Lifetime damage-equivalent mudline moment ${fmt(fat.del.mudlineMoment / 1e6)} MN·m.` : ''}`,
          gov.uc > 1 ? `Raise the monopile wall at the mudline from ${inp.monoT} mm to about ${fmt(Math.ceil(tNeed))} mm (or the diameter, which also raises the natural frequency), then re-run with every load case, six seeds and 600 s records.` : `The sections carry the cases run. Before relying on the loads: switch on the remaining cases (${CASES.filter(([id]) => !ids.includes(id)).map(([id]) => id).join(', ') || 'none left'}), use six 600 s seeds per condition and the site sea-state table, and pass the design mudline moment and shear to the foundation and strength modules.`,
          'IEC 61400-3-1:2019 Table 2 and §7.6; IEC 61400-1:2019 Table 3 (partial factors); DNV-ST-0126 (2021) §4; DNV-RP-C203 (2021).',
          `${steelCo2Text(steel)} Each millimetre of monopile wall is ${fmt(dSteel)} t of steel ≈ ${fmt(dSteel * CO2_STEEL)} t CO₂e: load reduction by control (peak shaving near rated, damping of the tower mode) saves that steel on every position.`),
        adv('economic', 'info', `The support structure holds ${fmt(steel)} t of steel; at the indicative ${fmt(inp.steelRate)} USD/t one millimetre of monopile wall is ${money(dSteel * inp.steelRate)} per position.`,
          gov.uc > 1 ? `Meeting the governing check by wall thickness alone costs about ${money(Math.max(tNeed - inp.monoT, 0) * dSteel * inp.steelRate)} per position; compare with a larger diameter or a load-reducing controller setting before adding steel.` : `The reserve of ${fmt((1 - gov.uc) * 100)} % on the governing check is worth up to ${money(Math.max(inp.monoT - tNeed, 0) * dSteel * inp.steelRate)} of steel per position if it survives the full load-case set; spend on the full simulation campaign before optimising it away.`,
          'AACE RP 18R-97 (Class 4 quantities from concept-stage loads); whole-life costing per ISO 15686-5; value of information of the full IEC load-case campaign.'),
      ],
      publish: { ...(gMud ? { 'owt.mudlineMoment': Md, 'owt.mudlineShear': Vd, 'owt.towerBaseMoment': Mt } : {}), 'owt.f1': f1 },
      quantities: { steel_t: round(steel) },
    };
  },
  verify: [
    { // IEC 61400-1 Table 3: the design value of a normal (N) situation is 1.35 × the characteristic extreme.
      // Extreme wind speed model: the 50-year 10-minute mean at hub height equals Vref (42.5 m/s for class II).
      name: 'DLC 6.1 parked: load factor 1.35 on the characteristic extreme and V50(hub) = Vref', source: 'IEC 61400-1:2019 Table 3 (normal design situation) and §6.3.3.2 (extreme wind speed model)',
      inputs: { run12: false, run13: false, run61: true, Vref: 42.5, duration: 30 },
      expect: [{ kpi: 'psf61', value: 1.35, tol: 1e-9 }, { kpi: 'psfGov', value: 1.35, tol: 1e-12 }, { kpi: 'v50', value: 42.5, tol: 1e-9 }] },
    { // Weibull A = 10 m/s, k = 2 (Rayleigh), NREL 5-MW cut-in 3 and cut-out 25 m/s:
      // P(3 < V < 25) = exp(−0.3²) − exp(−2.5²) = 0.913931 − 0.001930 = 0.912001.
      name: 'DLC 1.2 bins: Weibull probability of the operating wind range', source: 'Hand calculation: F(V) = 1 − exp(−(V/A)^k)',
      inputs: { run12: true, run13: false, run61: false, weibA: 10, weibK: 2, duration: 30, binWidth: 8 },
      expect: [{ kpi: 'opFrac', value: 91.2001, tol: 1e-5 }] },
  ],
};
