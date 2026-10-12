// Extreme environmental loading on a monopile-supported turbine: 50-year wind with the rotor parked, extreme and
// reduced wave heights with nonlinear kinematics and current, extreme operating gust near rated wind, breaking-wave
// slamming, air gap, and the governing factored shear and overturning moment at the mudline.
import { ewm, eogGust, ntmSigma, turbulenceScale } from '../../engine/aero/turbulence.js';
import { surfaceLookup } from '../../engine/aero/bem.js';
import { slamForcePerLength } from '../../engine/hydro/morison.js';
import { dynamicAmplification } from '../../engine/structural/framedynamics.js';
import { runDesignLoadCases } from '../../engine/aero/dlc.js';
import { need, run, fmt, round, adv, check, worst, sev, steelCo2Text, tube, turbineInputs, supportInputs, foundationInputs, steelInputs, steelOf, turbine, operatingPoint, LIMIT_AERO, LIMIT_COUPLED, LIMIT_CODE, LIMIT_HIFI, CONFIRM } from './_common.js';
import { waveCurrentLoad, towerWindLoad } from './_hydro.js';
import { dynamicsInputs, supportStructure } from './_dynamics.js';

export default {
  id: 's03-extreme-loading', suite: 3, title: 'Extreme environmental loading',
  summary: 'Governing extreme shear and overturning moment at the mudline from the 50-year storm with the rotor parked, the extreme operating gust in a severe sea state and breaking-wave slamming, with load factors and the platform air gap.',
  covers: [13], fidelity: 'analytical',
  standards: ['IEC 61400-3-1:2019 §7.4 Table 2 (DLC 1.6, 2.3, 6.1, 6.2), §7.6.2 (partial safety factors)', 'IEC 61400-1:2019 §6.3.3 (extreme wind and gust models)', 'DNV-ST-0437 (2016) §2.4 and §4', 'DNV-RP-C205 (2019) §3.3 (stream function), §6, §8.4 (wave slamming)', 'DNV-ST-0126 (2021) §4.3 (air gap)', 'Wienke & Oumeraci (2005) Coastal Engineering 52'],
  models: ['IEC extreme wind speed model (EWM) and extreme operating gust (EOG)', 'Stokes fifth-order and stream-function wave theory with integration to the free surface; Airy with Wheeler stretching as fallback', 'Morison equation with MacCamy–Fuchs inertia coefficient and nonlinear wave–current interaction', 'Parked-rotor thrust from blade element momentum theory at standstill', 'Frozen-pitch gust thrust from the coefficient surface', 'Breaking-wave slamming (Wienke–Oumeraci, curling factor)', 'Single-mode dynamic amplification of the wave load', 'Simplified IEC design-load-case subset on the coupled simulator (opt-in)'],
  limits: `Quasi-static combination of co-directional wind, wave and current extremes with a single-mode amplification of the wave load; the design-load-case driver is a simplified subset of IEC 61400-3-1 with one short seed per condition. ${LIMIT_AERO} ${LIMIT_COUPLED} ${LIMIT_CODE} ${LIMIT_HIFI}`,
  publishes: ['owt.mudlineMoment', 'owt.mudlineShear'],
  searchTerms: 'offshore wind 50-year storm extreme wave monopile load DLC 6.1 breaking wave slamming',
  inputs: [
    ...turbineInputs(),
    ...supportInputs({ embed: false }),
    ...steelInputs(),
    ...foundationInputs(),
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'Vref', label: 'Reference wind speed (50-year, 10-minute, hub height)', unit: 'm/s', type: 'number', default: 50, min: 25, max: 70, group: 'Extreme wind', role: 'bc', uncertain: { dist: 'gumbel', cov: 0.08 }, help: 'IEC class I 50, II 42.5, III 37.5 m/s; tropical-cyclone sites need a site-specific value.' },
    { key: 'Iref', label: 'Reference turbulence intensity', unit: '-', type: 'number', default: 0.14, min: 0.05, max: 0.25, group: 'Extreme wind', role: 'bc' },
    { key: 'cdTower', label: 'Tower drag coefficient', unit: '-', type: 'number', default: 0.6, min: 0.3, max: 1.2, group: 'Extreme wind', role: 'load' },
    { key: 'Hs50', label: 'Extreme significant wave height', unit: 'm', type: 'number', default: 6, min: 1, max: 18, group: 'Extreme sea', role: 'bc', bind: 'site.Hs100', uncertain: { dist: 'gumbel', cov: 0.1 }, help: '50-year (IEC) or 100-year value; depth-limited breaking is checked.' },
    { key: 'Tp50', label: 'Associated wave period', unit: 's', type: 'number', default: 10.5, min: 5, max: 22, group: 'Extreme sea', role: 'bc', bind: 'site.Tp100' },
    { key: 'Hmax', label: 'Extreme individual wave height (0 = 1.86 Hs)', unit: 'm', type: 'number', default: 0, min: 0, max: 35, group: 'Extreme sea', role: 'bc', bind: 'site.Hmax100' },
    { key: 'Uc50', label: 'Extreme surface current', unit: 'm/s', type: 'number', default: 1.2, min: 0, max: 4, group: 'Extreme sea', role: 'bc', bind: 'site.Uc100', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'surge', label: 'Storm surge', unit: 'm', type: 'number', default: 1.0, min: 0, max: 8, group: 'Extreme sea', role: 'bc', bind: 'site.surge100' },
    { key: 'tide', label: 'Tidal range', unit: 'm', type: 'number', default: 2, min: 0, max: 12, group: 'Extreme sea', role: 'bc', bind: 'site.tideRange', help: 'High water is taken as half the range above mean sea level.' },
    { key: 'Hs1', label: 'Severe sea state with the operating gust (Hs)', unit: 'm', type: 'number', default: 4, min: 0, max: 14, group: 'Operating extreme', role: 'bc', bind: 'site.Hs1', help: '1-year significant wave height combined with the extreme operating gust near rated wind.' },
    { key: 'Tp1', label: 'Its wave period', unit: 's', type: 'number', default: 8.5, min: 4, max: 20, group: 'Operating extreme', role: 'bc', bind: 'site.Tp1' },
    { key: 'Uc1', label: 'Current with the operating gust', unit: 'm/s', type: 'number', default: 0.8, min: 0, max: 3, group: 'Operating extreme', role: 'bc', bind: 'site.Uc1' },
    ...dynamicsInputs(),
    { key: 'mg', label: 'Marine growth thickness', unit: 'mm', type: 'number', default: 50, min: 0, max: 200, group: 'Dynamics', role: 'environment' },
    { key: 'slam', label: 'Add breaking-wave slamming to the storm cases', type: 'bool', default: false, group: 'Breaking waves', role: 'load', help: 'Required where plunging breakers can occur at the foundation (steep seabed slope or depth-limited waves).' },
    { key: 'curl', label: 'Curling factor λ', unit: '-', type: 'number', default: 0.4, min: 0.1, max: 0.6, group: 'Breaking waves', role: 'load', showIf: { key: 'slam', in: [true] }, help: 'Share of the breaker crest height that impacts the pile; about 0.4–0.5 for plunging breakers (Wienke & Oumeraci 2005).' },
    { key: 'gammaN', label: 'Load factor, normal cases', unit: '-', type: 'number', default: 1.35, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: `IEC 61400-3-1 Table 3, normal design situations. ${CONFIRM}` },
    { key: 'gammaA', label: 'Load factor, abnormal cases', unit: '-', type: 'number', default: 1.1, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: `IEC 61400-3-1 Table 3, abnormal design situations. ${CONFIRM}` },
    { key: 'gammaM', label: 'Material factor on yield', unit: '-', type: 'number', default: 1.1, min: 1, max: 1.5, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'runDlc', label: 'Also run the time-domain load-case subset (slow)', type: 'bool', default: false, group: 'Numerics', role: 'control', help: 'Runs simplified DLC 1.3, 2.3 and 6.1 on the coupled simulator with one short seed each and lists them beside the quasi-static cases.' },
    { key: 'nZ', label: 'Integration stations over depth', type: 'int', default: 40, min: 8, max: 320, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nZ', label: 'Depth stations', levels: [10, 20, 40, 80], qoi: 'waveMoment50' },
  compute(inp, ctx) {
    const steel = steelOf(inp, ctx);
    const tb = turbine(inp);
    const lever = inp.hubHeight + inp.depth;
    const level = inp.surge + 0.5 * inp.tide;
    const d = inp.depth + level;
    const Dh = inp.D + (2 * inp.mg) / 1e3;
    const q = (ct, V) => ct * 0.5 * tb.rhoAir * tb.area * V * V;
    const e = ewm({ Vref: inp.Vref, z: inp.hubHeight, hubHeight: inp.hubHeight });
    const tw = (V, shear) => towerWindLoad({ V, hubHeight: inp.hubHeight, shear, cd: inp.cdTower, rhoAir: tb.rhoAir, D: inp.D, zInt: inp.zInt, towerDb: inp.towerDb, towerDt: inp.towerDt, zTop: tb.zTowerTop, depth: inp.depth, level });
    // First-mode amplification of the wave load with the rotor parked (no aerodynamic damping).
    let f1 = NaN;
    try { f1 = supportStructure(inp, tb, steel, { nElem: 24, nModes: 1 }).modes[0].frequency; } catch { f1 = NaN; }
    const zeta = inp.zetaS / 100;
    const dafOf = (T) => (Number.isFinite(f1) ? dynamicAmplification({ ratio: 1 / T / f1, zeta }).daf : 1);
    const wave = (H, T, uc) => run(() => waveCurrentLoad({ H, T, d, D: Dh, Cd: inp.Cd, Cm: inp.Cm, mcf: inp.mcf, rho: inp.rhoW, uc0: uc, nZ: inp.nZ, theory: 'nonlinear', nPhase: 36 }), 'Extreme wave loading');
    const Hmax = inp.Hmax > 0 ? inp.Hmax : 1.86 * inp.Hs50;
    const w50 = wave(Hmax, inp.Tp50, inp.Uc50);
    const Hlim = Math.min(Hmax, 0.999 * w50.Hb), limited = Hmax >= w50.Hb;
    const wUse = limited ? wave(Hlim, inp.Tp50, inp.Uc50) : w50;
    const Hred = Math.min(1.3 * inp.Hs50, Hlim);
    const wRed = wave(Hred, inp.Tp50, inp.Uc50);
    const w1 = wave(Math.min(1.86 * inp.Hs1, Hlim), inp.Tp1, inp.Uc1);
    // Slamming on the storm cases (Wienke & Oumeraci): F = λ·ηb·ρ·R·C²·2π, applied at the crest.
    let slam = { F: 0, M: 0, peak: 0, tau: 0, dlf: 0 };
    if (inp.slam) {
      const C = wUse.L / inp.Tp50;
      const Fpk = slamForcePerLength({ D: Dh, v: C, Cs: 2 * Math.PI, rho: inp.rhoW }) * inp.curl * Math.max(wUse.crestElevation, 0);
      // The impact lasts τ = 13D/(64C) (Wienke & Oumeraci), far shorter than the natural period: the equivalent static
      // force is the peak times the impulse response factor π·f1·τ of a triangular pulse (≤ 1).
      const tau = (13 * Dh) / (64 * C);
      const dlf = Number.isFinite(f1) ? Math.min(1, Math.PI * f1 * tau) : 1;
      slam = { F: Fpk * dlf, M: Fpk * dlf * (d + wUse.crestElevation * (1 - 0.5 * inp.curl)), peak: Fpk, tau, dlf };
    }
    // Parked rotor: feathered and aligned (6.1), broadside after loss of yaw (6.2): flat-plate drag on the blade planform.
    // Two of three blades are taken as fully exposed (drag coefficient 1.3 on their planform).
    const ctAligned = tb.parkedCt, ctBroad = Math.max(tb.rotor.solidity * 1.3 * (2 / 3), ctAligned);
    // Extreme operating gust near rated wind with pitch frozen during the gust rise.
    const Vr = tb.windAtMaxThrust, op = operatingPoint(tb, Vr);
    const sig = ntmSigma({ Vhub: Vr, Iref: inp.Iref });
    const vg = eogGust({ Vhub: Vr, sigma1: sig, D: 2 * tb.R, lambda1: turbulenceScale(inp.hubHeight), Ve1: e.Ve1 });
    const Vpk = Vr + 0.74 * vg;
    const Tgust = q(surfaceLookup(tb.surface, (op.omega * tb.R) / Vpk, op.pitchDeg).ct, Vpk);
    const mk = (id, name, V, T, twL, wv, type, slamOn) => {
      const dW = dafOf(wv.T ?? inp.Tp50);
      const F = T + twL.F + wv.shearMax * dW + (slamOn ? slam.F : 0), M = T * lever + twL.M + wv.momentMax * dW + (slamOn ? slam.M : 0);
      const g = type === 'N' ? inp.gammaN : inp.gammaA;
      return { id, name, V, thrust: T, tower: twL.F, wave: wv.shearMax * dW, waveM: wv.momentMax * dW, F, M, type, g, Fd: g * F, Md: g * M, daf: dW };
    };
    const cases = [
      mk('6.1b', 'Parked, 50-year gust wind Ve50 + reduced wave height', e.Ve50, q(ctAligned, e.Ve50), tw(e.Ve50, 0.11), { ...wRed, T: inp.Tp50 }, 'N', inp.slam),
      mk('6.1c', 'Parked, reduced wind 1.1·Vref + extreme wave height', 1.1 * inp.Vref, q(ctAligned, 1.1 * inp.Vref), tw(1.1 * inp.Vref, 0.11), { ...wUse, T: inp.Tp50 }, 'N', inp.slam),
      mk('6.2', 'Parked, grid loss: rotor broadside, Ve50 + reduced wave height', e.Ve50, q(ctBroad, e.Ve50), tw(e.Ve50, 0.11), { ...wRed, T: inp.Tp50 }, 'A', inp.slam),
      mk('1.6 / 1.4', 'Production near rated, operating gust + severe sea state', Vpk, Tgust, tw(Vpk, 0.14), { ...w1, T: inp.Tp1 }, 'N', false),
      mk('2.3', 'Operating gust with loss of grid + severe sea state', Vpk, Tgust, tw(Vpk, 0.14), { ...w1, T: inp.Tp1 }, 'A', false),
    ];
    const gov = cases.reduce((a, c) => (c.Md > a.Md ? c : a), cases[0]);
    const govF = cases.reduce((a, c) => (c.Fd > a.Fd ? c : a), cases[0]);
    // Optional time-domain subset through the design-load-case driver.
    let dlc = null;
    if (inp.runDlc) {
      const str = supportStructure(inp, tb, steel, { nElem: 24, nModes: 2 });
      dlc = run(() => runDesignLoadCases({
        model: { surface: tb.surface, turbine: { inertia: tb.rotor.drivetrainInertia, hubHeight: inp.hubHeight }, controller: tb.controller, structure: str, rho: tb.rhoAir, rhoWater: inp.rhoW, dt: 0.05 },
        site: { hubHeight: inp.hubHeight, rotorDiameter: 2 * tb.R, waterDepth: d, Vref: inp.Vref, Iref: inp.Iref, ratedWind: tb.windAtMaxThrust, cutIn: tb.rotor.cutIn, cutOut: tb.rotor.cutOut,
          seaState: () => ({ Hs: inp.Hs1 * 0.6, Tp: inp.Tp1 }), Hs50: inp.Hs50, Tp50: inp.Tp50, Hs1: inp.Hs1, Tp1: inp.Tp1, macCamyFuchs: inp.mcf },
        cases: ['1.3', '2.3', '6.1'], nSeeds: 1, duration: 100, tTransient: 20, windSpeeds: [tb.windAtMaxThrust], parkedCt: ctAligned, yawErrors61Deg: [0],
      }), 'Time-domain design load cases');
    }
    const pile = tube(inp.D, inp.t / 1e3);
    const sigD = gov.Md / pile.W, fyd = steel.fy / inp.gammaM;
    const crestAbs = level + wUse.crestElevation;
    const gapReq = Math.max(0.2 * inp.Hs50, 1), gap = inp.zInt - crestAbs;
    const checks = [
      check('Monopile bending stress at the mudline (governing factored moment)', sigD / 1e6, fyd / 1e6, 'MPa', 'IEC 61400-3-1:2019 §7.6; DNV-ST-0126 (2021) §4.4', 'Elastic bending stress in the pile wall at the mudline under the governing factored extreme. Exceedance means yielding in the design storm; shell buckling is checked separately.'),
      check('Air gap between the extreme crest and the interface platform', gapReq, Math.max(gap, 1e-9), 'm (required / available)', 'DNV-ST-0126 (2021) §4.3.2: air gap ≥ 20 % of Hs,50, minimum 1 m', 'Guards the platform, boat landing top and tower door against wave impact. A crest reaching the platform adds large slamming loads that are not in the design loads.'),
      check('Extreme wave height against the breaking limit at the site', Hmax, w50.Hb, 'm', 'DNV-RP-C205 (2019) §3.4.6', 'Above the limit the extreme wave is depth-limited: the height is capped at the breaking height and plunging-breaker slamming must be considered.'),
    ];
    checks[1].uc = gap > 0 ? Math.min(gapReq / gap, 99) : 99;
    const gc = worst(checks.slice(0, 2)), ok = gc.uc <= 1;
    const tReq = inp.t * sigD / fyd;
    const dMass = (tube(inp.D, Math.min(Math.max(tReq, 15), 0.4 * inp.D * 1e3) / 1e3).A - pile.A) * steel.rho * (inp.depth + inp.zInt) / 1e3;
    return {
      kpis: [
        { key: 'designMoment', label: 'Governing design moment at the mudline', unit: 'MN·m', value: gov.Md / 1e6, explain: `Largest factored overturning moment over the load cases (case ${gov.id}, load factor ${gov.g}). This is the ULS action for the pile, the soil and the grouted or bolted connection.` },
        { key: 'designShear', label: 'Governing design shear at the mudline', unit: 'MN', value: govF.Fd / 1e6, explain: `Largest factored horizontal force over the load cases (case ${govF.id}).` },
        { key: 'charMoment', label: 'Characteristic moment of the governing case', unit: 'MN·m', value: gov.M / 1e6, explain: 'Unfactored moment of the case that governs after factoring; published for the geotechnical checks, which apply their own factors.' },
        { key: 'charShear', label: 'Characteristic shear of the governing case', unit: 'MN', value: gov.F / 1e6, explain: 'Unfactored horizontal force acting with the characteristic moment.' },
        { key: 've50', label: '50-year 3-second gust at hub height', unit: 'm/s', value: e.Ve50, explain: 'IEC steady extreme wind model: 1.4 times the reference wind speed at hub height.' },
        { key: 'thrustParked', label: 'Parked rotor thrust in the 50-year gust', unit: 'kN', value: q(ctAligned, e.Ve50) / 1e3, explain: 'Drag of the feathered, aligned rotor. Small compared with operating thrust — provided yaw control keeps the rotor aligned.' },
        { key: 'thrustGust', label: 'Rotor thrust in the extreme operating gust', unit: 'kN', value: Tgust / 1e3, explain: `Thrust at the gust peak (${fmt(Vpk)} m/s) with pitch and rotor speed frozen at the rated-wind operating point: an upper estimate of the controller's lag.` },
        { key: 'waveShear50', label: 'Extreme wave + current force', unit: 'MN', value: wUse.shearMax / 1e6, explain: 'Largest horizontal force of the extreme wave with current over the wave cycle, from nonlinear kinematics integrated to the free surface, before dynamic amplification.' },
        { key: 'waveMoment50', label: 'Extreme wave + current moment', unit: 'MN·m', value: wUse.momentMax / 1e6, explain: 'Largest overturning moment of the extreme wave with current about the mudline, before dynamic amplification.' },
        { key: 'crest', label: 'Extreme crest elevation above mean sea level', unit: 'm', value: crestAbs, limit: inp.zInt - gapReq, sense: 'max', explain: 'Crest of the extreme wave on top of surge and high tide. It must stay below the interface platform by the required air gap.' },
        { key: 'dafWave', label: 'Dynamic amplification of the storm wave load', unit: '-', value: dafOf(inp.Tp50), explain: 'Amplification of wave loading at the storm wave period by the first structural mode with the rotor parked. Grows as the natural period approaches the wave period.' },
        { key: 'slamForce', label: 'Breaking-wave slamming force (peak)', unit: 'MN', value: slam.peak / 1e6, explain: `Peak impulsive force of a plunging breaker on the pile (zero unless slamming is switched on). It lasts about ${fmt(slam.tau * 1e3)} ms, so the load cases carry the equivalent static value, the peak times the impulse response factor ${fmt(slam.dlf)}.` },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Design moment at the mudline by load case', categories: cases.map((c) => c.id), series: [{ name: 'Rotor', y: cases.map((c) => c.g * c.thrust * lever / 1e6) }, { name: 'Tower drag', y: cases.map((c) => c.g * (c.M - c.thrust * lever - c.waveM - (c.type && inp.slam && c.id.startsWith('6') ? slam.M : 0)) / 1e6) }, { name: 'Waves and current', y: cases.map((c) => c.g * c.waveM / 1e6) }], x: { label: 'Load case' }, y: { label: 'Factored moment', unit: 'MN·m' },
          explain: 'Factored overturning moment of each design load case split by source. Storm cases are wave-dominated, operating cases thrust-dominated; the tallest total governs.' },
        { type: 'line', title: 'Extreme wave force over depth', x: { label: 'Force per length', unit: 'kN/m' }, y: { label: 'Elevation above storm water level', unit: 'm' },
          series: [{ name: `Extreme wave H = ${fmt(Hlim)} m (${wUse.theory})`, x: wUse.f.map((v) => v / 1e3), y: wUse.z }, { name: `Reduced wave H = ${fmt(Hred)} m`, x: wRed.f.map((v) => v / 1e3), y: wRed.z, dash: true }],
          explain: 'Hydrodynamic force per metre at the phase of maximum moment, up to the instantaneous free surface. The crest zone carries the highest load and the longest lever arm.' },
      ],
      tables: [
        { title: 'Design load cases (quasi-static)', columns: [{ key: 'id', label: 'DLC' }, { key: 'name', label: 'Situation' }, { key: 'V', label: 'Hub wind', unit: 'm/s' }, { key: 'T', label: 'Rotor thrust', unit: 'kN' }, { key: 'W', label: 'Wave + current', unit: 'kN' }, { key: 'F', label: 'Shear, char.', unit: 'MN' }, { key: 'M', label: 'Moment, char.', unit: 'MN·m' }, { key: 'g', label: 'γf' }, { key: 'Md', label: 'Moment, design', unit: 'MN·m' }],
          rows: cases.map((c) => ({ id: c.id, name: c.name, V: round(c.V), T: round(c.thrust / 1e3), W: round(c.wave / 1e3), F: round(c.F / 1e6), M: round(c.M / 1e6), g: c.g, Md: round(c.Md / 1e6) })),
          explain: 'One row per design situation with characteristic and factored loads at the mudline. Use the governing row for the ULS checks and keep the others for the load envelope.' },
        ...(dlc ? [{ title: 'Time-domain load cases on the coupled simulator (one short seed)', columns: [{ key: 'dlc', label: 'DLC' }, { key: 'd', label: 'Situation' }, { key: 'Mc', label: 'Mudline moment, char.', unit: 'MN·m' }, { key: 'Md', label: 'Mudline moment, design', unit: 'MN·m' }, { key: 'Fc', label: 'Mudline shear, char.', unit: 'MN' }, { key: 'gv', label: 'Governing run' }],
          rows: dlc.table.map((r) => ({ dlc: r.dlc, d: r.description, Mc: round(r.mudlineMomentChar / 1e6), Md: round(r.mudlineMomentDesign / 1e6), Fc: round(r.mudlineShearChar / 1e6), gv: r.governing })),
          explain: 'Extremes from the simplified time-domain load cases (linear waves, rigid blades, one seed of 100 s) for comparison with the quasi-static rows: production cases include controller dynamics, storm cases include resonance but not nonlinear crest loading.' }] : []),
        { title: 'Extreme wave parameters', columns: [{ key: 'k', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['Storm water depth (surge + high tide)', d, 'm'], ['Extreme wave height used', Hlim, 'm'], ['Breaking height at the site', w50.Hb, 'm'], ['Wavelength', wUse.L, 'm'], ['Crest above storm water level', wUse.crestElevation, 'm'], ['Wave theory used', wUse.theory, ''], ['Recommended theory', wUse.recommended, ''], ['Effective inertia coefficient', wUse.cmEff, '-'], ['First natural frequency', Number.isFinite(f1) ? f1 : 0, 'Hz'], ['Gust magnitude of the EOG', vg, 'm/s'], ['Broadside rotor thrust coefficient (6.2)', ctBroad, '-'], ['Slamming: equivalent static force', slam.F / 1e6, 'MN']].map(([k, v, u]) => ({ k, v: typeof v === 'number' ? round(v) : v, u })),
          explain: 'Wave, water-level and gust parameters behind the load cases, for checking against the metocean design basis.' },
      ],
      notes: [`Case ${gov.id} governs the moment (${gov.name.toLowerCase()}): ${fmt(gov.M / 1e6)} MN·m characteristic, ${fmt(gov.Md / 1e6)} MN·m factored.`,
        limited ? `The extreme wave is depth-limited: H = ${fmt(Hmax)} m exceeds the breaking height ${fmt(w50.Hb)} m and was capped; consider slamming.` : `The extreme wave (${fmt(Hmax)} m) is below the breaking height (${fmt(w50.Hb)} m).`,
        `Wave kinematics: ${wUse.theory}; storm wave load amplified by ${fmt(dafOf(inp.Tp50))} for the first mode.`, ...(dlc ? dlc.notes : []),
        'Rotor loads come from generic polars: replace them with the supplier\'s extreme load tables before design.'],
      advice: [
        ok ? adv('engineering', sev(gc.uc), `Governing design moment ${fmt(gov.Md / 1e6)} MN·m (case ${gov.id}); pile stress ${fmt(sigD / 1e6)} MPa is ${fmt(sigD / fyd * 100)} % of design yield and the air gap is ${fmt(gap)} m against ${fmt(gapReq)} m required.`,
          `Carry M = ${fmt(gov.M / 1e6)} MN·m and H = ${fmt(gov.F / 1e6)} MN (characteristic) into the pile–soil, capacity and buckling checks, as preliminary loads to be replaced by the supplier's extreme load tables; confirm the load factors against the governing edition of IEC 61400-3-1. ${limited || checks[2].uc > 0.9 ? 'Waves are near the breaking limit: run the slamming option and a dynamic check of the impulsive load.' : 'Confirm the storm wave period range: the shortest associated period gives the highest amplification.'}`,
          'IEC 61400-3-1:2019 §7.4 and §7.6.2; DNV-ST-0437 (2016) §4; DNV-RP-C205 (2019) §8.4',
          `At the governing moment a wall of about ${fmt(Math.max(tReq, 15))} mm satisfies yield at the mudline, ${fmt(dMass)} t relative to the present wall above the mudline — ${steelCo2Text(dMass)}`)
          : adv('engineering', 'fail', `${gc.name}: ${fmt(gc.demand)} against ${fmt(gc.capacity)} ${gc.unit}.`,
            gc.name.includes('Air gap') ? `Raise the interface level to at least ${fmt(crestAbs + gapReq)} m above mean sea level, or design the platform and attachments for wave impact.` : `Increase the monopile wall at the mudline to about ${fmt(tReq * 1.02)} mm or the diameter (stress ∝ 1/D²t), or use a higher steel grade.`,
            'DNV-ST-0126 (2021) §4.3.2 and §4.4; IEC 61400-3-1:2019 §7.6',
            `Raising the interface by ${fmt(Math.max(crestAbs + gapReq - inp.zInt, 0))} m adds ${fmt(Math.max(crestAbs + gapReq - inp.zInt, 0) * pile.A * steel.rho / 1e3)} t of pile; the wall change is ${fmt(dMass)} t — ${steelCo2Text(dMass)}`),
        adv('economic', 'info', `Waves and current give ${fmt(gov.waveM / gov.M * 100)} % of the governing characteristic moment; the factored storm moment is ${fmt(gov.Md / Math.max(cases[3].Md, 1) * 100)} % of the factored operating-gust moment.`,
          gov.id === '6.2' && gov.waveM < 0.5 * gov.M ? 'Loss of yaw control in the storm governs through drag on the misaligned parked rotor: back-up power for the yaw system (so that the rotor stays aligned, case 6.1) removes this case for the cost of a battery or diesel set, far less than the steel it drives. Confirm the misaligned-rotor load with the supplier.' : gov.id.startsWith('6') ? 'The storm governs: metocean extremes (Hs, surge, current) deserve site-specific hindcast and extreme-value analysis, since a 10 % lower extreme wave height reduces the governing moment and the steel roughly in proportion.' : 'Operation governs: negotiate thrust-limiting control (peak shaving) with the turbine supplier — a few per cent of thrust removed near rated wind reduces foundation steel for well under 1 % of energy yield.',
          'Value-of-information for metocean studies; LCOE trade-off between energy yield and foundation CAPEX (AACE RP 18R-97 Class 4 quantities)'),
      ],
      publish: { 'owt.mudlineMoment': gov.M, 'owt.mudlineShear': gov.F },
    };
  },
  verify: [
    { name: 'IEC extreme wind model: Ve50 = 1.4·Vref at hub height', source: 'IEC 61400-1:2019 Eq. (15): 1.4 × 50 = 70 m/s', inputs: { Vref: 50 }, expect: [{ kpi: 've50', value: 70, tol: 1e-9 }] },
    { name: 'Small-amplitude limit of the nonlinear wave load: F = ρ·Cm·(πD²/4)·g·(H/2)·tanh(kd)', source: 'Morison inertia force in a linear wave (Dean & Dalrymple 1991): H = 0.2 m, T = 10 s, d = 25 m, D = 6 m, Cm = 2 → k = 0.04820 1/m, F = 1025 × 2 × 28.274 × 9.80665 × 0.1 × 0.8352 = 47.47 kN',
      inputs: { Hs50: 1, Hmax: 0.2, Tp50: 10, Uc50: 0, surge: 0, tide: 0, mg: 0, Cd: 0, Cm: 2, mcf: false, D: 6, depth: 25, nZ: 80 }, expect: [{ kpi: 'waveShear50', value: 0.04747, tol: 0.01 }] },
  ],
};
