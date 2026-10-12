// Extreme environmental load analysis: IEC extreme wind, gust and sea-state models combined into deterministic
// design load cases for a bottom-fixed offshore wind turbine (power production with extreme gust, parked in the
// 50-year storm with reduced or extreme wave height, parked with loss of yaw control, optional breaking-wave
// slam), with the mudline and tower-base loads, the crest elevation and air gap, and a survival check of the pile.
import { ewm, eogGust, eogSpeed, ntmSigma, etmSigma, turbulenceScale } from '../../engine/aero/turbulence.js';
import { LOAD_SAFETY_FACTORS } from '../../engine/aero/dlc.js';
import { bemSolve, surfaceLookup } from '../../engine/aero/bem.js';
import { tubularSectionStress } from '../../engine/aero/blade.js';
import { wavenumber, breakingLimits } from '../../engine/metocean/waves.js';
import { slamForcePerLength } from '../../engine/hydro/morison.js';
import { wavePileLoad, cantilever, towerWind } from './_loads.js';
import { turbineInputs, getTurbine, steadyAt, supportInputs, supportSegments, site, adv, check, worst, sev, fmt, round, run, need, linspace, G,
  LIMIT_POLAR, LIMIT_CODE, CONFIRM, SUPPLIER, steelCo2Text } from './_common.js';

export default {
  id: 's16-extreme-loads', suite: 16, title: 'Extreme environmental loads and survival',
  summary: 'Characteristic and design mudline loads of a bottom-fixed turbine for the governing extreme load cases — operating gust, 50-year storm parked, loss of yaw control and breaking-wave slam — with crest elevation, air gap and a survival check.',
  covers: [10], fidelity: 'analytical',
  standards: ['IEC 61400-1:2019 §6.3.2 (extreme wind speed model, extreme operating gust, extreme turbulence) and Table 3 (load factors)', 'IEC 61400-3-1:2019 §6.4 and Table 1 (marine conditions and design load cases 1.3, 1.6, 2.3, 6.1b, 6.1c, 6.2)', 'DNV-ST-0437 (2016) §4 (load cases)', 'DNV-RP-C205 (2021) §3 and §8 (extreme waves, slamming)', 'Wienke & Oumeraci (2005) Coastal Engineering 52 (breaking-wave impact on a slender pile)'],
  models: ['IEC extreme wind speed model (steady and reduced)', 'IEC extreme operating gust', 'Normal and extreme turbulence models', 'Rayleigh most-probable-maximum and reduced wave heights', 'Stokes second-order crest elevation', 'Depth-limited breaking (Miche, McCowan)', 'Morison equation with Airy kinematics and Wheeler stretching, MacCamy–Fuchs inertia coefficient', 'Slamming force ½ρCsDv² with curling factor', 'Rotor thrust from Ct surfaces (frozen controller during the gust)', 'Partial load factors'],
  limits: `Deterministic, quasi-static load combinations with co-directional wind, waves and current: turbulence-driven extremes, dynamic response and statistical extrapolation require the coupled simulator with at least six seeds per case. The parked-rotor thrust coefficient is an input (steady BEM at standstill with yaw error is outside its validity); slam loads use a constant slamming coefficient and curling factor without the dynamic response of the structure. ${LIMIT_POLAR} ${LIMIT_CODE}`,
  publishes: ['owt.mudlineMoment', 'owt.mudlineShear', 'owt.mudlineMomentDesign', 'owt.crestElevation'],
  searchTerms: 'IEC 61400-3 design load cases extreme wind speed model 50-year wave monopile',
  inputs: [
    ...turbineInputs(),
    ...supportInputs(),
    { key: 'Vref', label: 'Reference wind speed Vref (10-min, 50-year, hub height)', unit: 'm/s', type: 'number', default: 50, min: 25, max: 70, group: 'Extreme wind', role: 'bc', uncertain: { dist: 'normal', cov: 0.07 }, help: `IEC wind turbine class: I 50, II 42.5, III 37.5 m/s; class T (typhoon) 57 m/s. Use the site 50-year value when known. ${CONFIRM}` },
    { key: 'Iref', label: 'Reference turbulence intensity', unit: '-', type: 'number', default: 0.12, min: 0.05, max: 0.2, group: 'Extreme wind', role: 'bc', help: 'IEC category A 0.16, B 0.14, C 0.12.' },
    { key: 'dafOp', label: 'Dynamic amplification of the operating thrust', unit: '-', type: 'number', default: 1.2, min: 1, max: 2, group: 'Extreme wind', role: 'bc', help: 'Applied to the gust thrust for tower dynamics (the gust itself comes from the extreme operating gust model).' },
    { key: 'ctParked', label: 'Thrust coefficient of the parked, feathered rotor', unit: '-', type: 'number', default: 0.014, min: 0.002, max: 0.1, group: 'Extreme wind', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.3 }, help: `Referred to the swept area, including ±8° yaw error and nacelle drag: about 0.01–0.02 for a three-bladed feathered rotor. ${SUPPLIER}` },
    site('Hs100', { group: 'Extreme sea state' }), site('Tp100', { group: 'Extreme sea state' }), site('Hs1', { group: 'Extreme sea state' }), site('Tp1', { group: 'Extreme sea state' }), site('Uc100', { group: 'Extreme sea state' }), site('Uc1', { group: 'Extreme sea state' }),
    { key: 'surge', label: 'Storm surge plus highest tide above mean sea level', unit: 'm', type: 'number', default: 2.5, min: 0, max: 10, group: 'Extreme sea state', role: 'bc', bind: 'site.surge100', help: 'Raises the still-water level for the storm cases and the air-gap check.' },
    { key: 'cdBroad', label: 'Drag coefficient of the blades broadside to the wind', unit: '-', type: 'number', default: 1.3, min: 0.8, max: 2, group: 'Extreme wind', role: 'bc', help: 'Loss of yaw control (DLC 6.2): flat-plate drag on two thirds of the blade planform. A screening value — the supplier load document governs.' },
    { key: 'Cd', label: 'Drag coefficient of the pile', unit: '-', type: 'number', default: 1.05, min: 0.5, max: 1.6, group: 'Hydrodynamic coefficients', role: 'bc', help: 'Rough (marine growth): 1.05 per DNV-RP-C205.' },
    { key: 'Cm', label: 'Inertia coefficient of the pile', unit: '-', type: 'number', default: 2.0, min: 1.2, max: 2.2, group: 'Hydrodynamic coefficients', role: 'bc' },
    { key: 'slam', label: 'Include breaking-wave slam in the extreme-wave case', type: 'bool', default: false, group: 'Hydrodynamic coefficients', role: 'bc', help: 'Relevant where the extreme wave is close to the depth-limited breaking height or the seabed slopes.' },
    { key: 'curl', label: 'Curling factor λ', unit: '-', type: 'number', default: 0.4, min: 0.1, max: 0.6, group: 'Hydrodynamic coefficients', role: 'bc', showIf: { key: 'slam', in: [true] }, help: 'Fraction of the crest height that hits the pile: 0.4–0.5 for plunging breakers (Wienke & Oumeraci 2005).' },
    { key: 'gammaN', label: 'Load factor, normal cases', unit: '-', type: 'number', default: LOAD_SAFETY_FACTORS.N, min: 1, max: 1.7, group: 'Acceptance', role: 'criteria', help: `IEC 61400-1 Table 3: 1.35. ${CONFIRM}` },
    { key: 'gammaA', label: 'Load factor, abnormal cases', unit: '-', type: 'number', default: LOAD_SAFETY_FACTORS.A, min: 1, max: 1.5, group: 'Acceptance', role: 'criteria', help: `IEC 61400-1 Table 3: 1.1. ${CONFIRM}` },
    { key: 'gammaM', label: 'Material factor', unit: '-', type: 'number', default: 1.1, min: 1, max: 1.3, group: 'Acceptance', role: 'criteria' },
    { key: 'gapMin', label: 'Required air gap below the transition platform', unit: 'm', type: 'number', default: 1.5, min: 0, max: 5, group: 'Acceptance', role: 'criteria', help: 'DNV-ST-0126: 20 % of Hs,50 with a minimum of 1 m is common practice.' },
  ],
  compute(inp, ctx) {
    const tb = getTurbine(inp), rotor = tb.rotor, R = rotor.tipRadius, A = tb.surface.area;
    need(inp.hubHeight > R + 5, 'The hub height leaves no clearance under the blade tip.');
    const st = ctx.material(inp.steel);
    need(st && st.fy > 0, `Steel '${inp.steel}' is not in the materials library.`);
    const segs = supportSegments(inp, { zHub: inp.hubHeight, hubOffset: 2.4 * tb.scale, E: st.E });
    const x = ewm({ Vref: inp.Vref, z: inp.hubHeight, hubHeight: inp.hubHeight });
    // --- operating gust (DLC 2.3 / 1.3 style): rated wind, rotor speed and pitch frozen during the 10.5 s gust ---
    const Vr = tb.ratedWind, op = steadyAt(tb, Vr);
    const sigma1 = ntmSigma({ Vhub: Vr, Iref: inp.Iref });
    const vg = eogGust({ Vhub: Vr, sigma1, D: 2 * R, lambda1: turbulenceScale(inp.hubHeight), Ve1: x.Ve1 });
    const Vpk = Vr + 0.74 * vg;
    const ctG = surfaceLookup(tb.surface, op.omega * tb.surface.R / Vpk, op.pitchDeg).ct;
    const Tgust = 0.5 * 1.225 * A * Vpk * Vpk * ctG * inp.dafOp;
    // extreme turbulence (DLC 1.3): peak of the turbulent thrust with a 3σ gust on the frozen operating point
    const sEtm = etmSigma({ Vhub: Vr, Iref: inp.Iref, Vave: 0.2 * inp.Vref });
    const Vetm = Vr + 2 * sEtm, Tetm = 0.5 * 1.225 * A * Vetm * Vetm * surfaceLookup(tb.surface, op.omega * tb.surface.R / Vetm, op.pitchDeg + 0.6 * (Vetm - Vr)).ct * inp.dafOp;
    // --- parked rotor ---
    const Tpark = (V) => 0.5 * 1.225 * A * V * V * inp.ctParked;
    const edge = run(() => bemSolve({ rotor, U: x.Ve50, omega: 0, pitchDeg: 90 }), 'Parked-rotor BEM');
    const planform = rotor.solidity * rotor.sweptArea;
    const Tside = 0.5 * 1.225 * inp.cdBroad * planform * (2 / 3) * x.Ve50 ** 2;    // loss of yaw: blades broadside to the wind, flat-plate drag on two thirds of the planform
    // --- waves ---
    const d50 = inp.depth + inp.surge;
    const brk = breakingLimits({ T: inp.Tp100, d: d50 });
    const cap = (H, dep = d50) => Math.min(H, 0.98 * brk.governing, 0.85 * dep);   // no wave higher than the depth-limited breaking height
    const Hext = cap(1.86 * inp.Hs100), Hred = cap(1.3 * inp.Hs100), H1 = cap(1.86 * inp.Hs1);
    const k50 = wavenumber({ T: inp.Tp100, d: d50 }), L50 = 2 * Math.PI / k50;
    const crest = Hext / 2 + k50 * Hext * Hext / 16 * Math.cosh(k50 * d50) * (2 + Math.cosh(2 * k50 * d50)) / Math.sinh(k50 * d50) ** 3;
    const crestElev = inp.surge + crest, gap = inp.transZ - crestElev;
    const cases = [
      { id: 'DLC 1.3', name: 'Power production, extreme turbulence, normal sea', T: Tetm, Uw: Vetm, alpha: 0.14, H: cap(1.86 * inp.Hs1 * 0.6, inp.depth), Tp: inp.Tp1, Uc: inp.Uc1, depth: inp.depth, g: inp.gammaN, type: 'N' },
      { id: 'DLC 1.6', name: 'Power production at rated, severe sea state', T: tb.thrustRated * inp.dafOp, Uw: Vr, alpha: 0.14, H: cap(Hext, inp.depth), Tp: inp.Tp100, Uc: inp.Uc1, depth: inp.depth, g: inp.gammaN, type: 'N' },
      { id: 'DLC 2.3', name: 'Extreme operating gust with loss of grid', T: Tgust, Uw: Vpk, alpha: 0.14, H: cap(H1 * 0.6, inp.depth), Tp: inp.Tp1, Uc: inp.Uc1, depth: inp.depth, g: inp.gammaA, type: 'A' },
      { id: 'DLC 6.1b', name: 'Parked, steady extreme wind Ve50, reduced wave height', T: Tpark(x.Ve50), Uw: x.Ve50, alpha: 0.11, H: Hred, Tp: inp.Tp100, Uc: inp.Uc100, depth: d50, g: inp.gammaN, type: 'N' },
      { id: 'DLC 6.1c', name: 'Parked, reduced wind 1.1·Vref, extreme wave height', T: Tpark(1.1 * inp.Vref), Uw: 1.1 * inp.Vref, alpha: 0.11, H: Hext, Tp: inp.Tp100, Uc: inp.Uc100, depth: d50, g: inp.gammaN, type: 'N', slam: inp.slam },
      { id: 'DLC 6.2', name: 'Parked, loss of yaw control, blades broadside', T: Tside, Uw: x.Ve50, alpha: 0.11, H: Hred, Tp: inp.Tp100, Uc: inp.Uc100, depth: d50, g: inp.gammaA, type: 'A' },
    ];
    let slamRec = null;
    const res = cases.map((c) => {
      const sg = c.depth === inp.depth ? segs : supportSegments({ ...inp, depth: c.depth, transZ: inp.transZ - inp.surge }, { zHub: inp.hubHeight - inp.surge, hubOffset: 2.4 * tb.scale, E: st.E });
      const zHub = c.depth === inp.depth ? inp.hubHeight : inp.hubHeight - inp.surge;
      const wv = run(() => wavePileLoad({ H: c.H, T: c.Tp, depth: c.depth, D: inp.monoD, Cd: inp.Cd, Cm: inp.Cm, Uc: c.Uc, n: 50 }), `${c.id} wave load`);
      const pts = [{ z: zHub, F: c.T, M: 0 }];
      let Fslam = 0;
      if (c.slam) {
        // Wienke & Oumeraci (2005): F = λ·η_b·ρ·R·C²·2π at first contact, i.e. slam force per length with Cs = 2π and v = wave celerity over the curled height
        const C = (2 * Math.PI / c.Tp) / wv.k, etaB = crest;
        Fslam = slamForcePerLength({ D: inp.monoD, v: C, Cs: 2 * Math.PI }) * inp.curl * etaB;
        pts.push({ z: etaB * (1 - inp.curl / 2), F: Fslam, M: 0 });
        slamRec = { F: Fslam, C, etaB };
      }
      const sol = cantilever({ segs: sg, n: 60, points: pts, distributed: [wv], wind: towerWind({ segs: sg, Uhub: c.Uw, zHub, alpha: c.alpha }), topWeight: rotor.rnaMass * G, foundation: 'fixed', E: st.E, depth: c.depth });
      const lever = zHub + c.depth, Mthr = c.T * lever, Mwave = wv.moment + (Fslam ? Fslam * (crest * (1 - inp.curl / 2) + c.depth) : 0), Mwind = sol.M[0] - Mthr - Mwave;
      const iTB = sol.z.reduce((a, z, i) => (Math.abs(z - (inp.transZ - (c.depth - inp.depth))) < Math.abs(sol.z[a] - (inp.transZ - (c.depth - inp.depth))) ? i : a), 0);
      return { ...c, V: sol.V[0], M: sol.M[0], Mthr, Mwave, Mwind, Vwave: wv.shear + Fslam, Md: sol.M[0] * c.g, Vd: sol.V[0] * c.g, Mtb: sol.M[iTB], N: sol.N[0], wave: wv };
    });
    const gov = res.reduce((a, r) => (r.Md > a.Md ? r : a), res[0]), govChar = res.reduce((a, r) => (r.M > a.M ? r : a), res[0]);
    const sec = tubularSectionStress({ D: inp.monoD, t: inp.monoT / 1e3, N: gov.N * 1.1, M: gov.Md, V: gov.Vd });
    const checks = [
      check(`Monopile yield at the mudline, ${gov.id}`, sec.vonMises / 1e6, st.fy / inp.gammaM / 1e6, 'MPa', 'DNV-ST-0126 (2021) §4.4 (ULS); IEC 61400-3-1:2019 §7.6', 'Factored von Mises stress at the mudline section in the governing extreme case against design yield. Above 1 the pile would yield in the design storm: survival is not demonstrated.'),
      check('Platform elevation against the extreme crest plus air gap', crestElev + inp.gapMin, inp.transZ, 'm', 'DNV-ST-0126 (2021) §4.3 (air gap); IEC 61400-3-1:2019 §6.4', 'Demand is the extreme crest elevation (surge and tide included) plus the required air gap, capacity the elevation of the transition platform. If the crest reaches the platform, wave-in-deck slamming loads the platform, boat landing and J-tubes.'),
      check('Extreme wave height against the breaking limit', 1.86 * inp.Hs100, brk.governing, 'm', 'DNV-RP-C205 (2021) §3.4.6 (Miche / depth-limited breaking)', 'Rayleigh maximum wave height against the height at which a wave of this period breaks in this depth. Above 1 the wave is depth-limited (capped in the loads) and breaking-wave slam must be considered.'),
    ];
    const gc = worst(checks);
    const tg = linspace(-2, 12.5, 59), eog = tg.map((t) => eogSpeed(t, { Vhub: Vr, Vgust: vg }));
    const zz = linspace(10, inp.hubHeight + R, 30);
    const tReq = inp.monoT * checks[0].uc;
    const pileMass = Math.PI * (inp.monoD - inp.monoT / 1e3) * inp.monoT / 1e3 * (inp.depth + inp.transZ) * 7850 / 1e3;
    const mn = (v) => round(v / 1e6, 5);
    return {
      kpis: [
        { key: 'mudMoment', label: 'Governing characteristic mudline moment', value: mn(govChar.M), unit: 'MN·m', explain: `Largest unfactored overturning moment over the load cases (${govChar.id}). This is the characteristic load passed to the foundation and geotechnical modules.` },
        { key: 'mudShear', label: 'Mudline shear in that case', value: mn(govChar.V), unit: 'MN', explain: 'Unfactored horizontal force at the mudline in the case that governs the moment.' },
        { key: 'mudMomentDesign', label: 'Governing design mudline moment', value: mn(gov.Md), unit: 'MN·m', explain: `Largest factored moment (${gov.id}, γf = ${gov.g}). Abnormal cases carry 1.1 and may not govern after factoring even when their characteristic load is highest.` },
        { key: 'towerBaseMoment', label: 'Characteristic tower-base moment (governing case)', value: mn(Math.max(...res.map((r) => r.Mtb))), unit: 'MN·m', explain: 'Largest unfactored moment at the tower base: the interface load for the transition piece and flange.' },
        { key: 've50', label: '50-year extreme wind speed Ve50 (3-s gust, hub)', value: round(x.Ve50, 4), unit: 'm/s', explain: 'IEC steady extreme wind model: 1.4 × Vref at hub height. The turbine is parked and feathered at this speed.' },
        { key: 'vgust', label: 'Extreme operating gust magnitude at rated wind', value: round(vg, 4), unit: 'm/s', explain: 'IEC gust amplitude; the wind peaks at Vhub + 0.74·Vgust half-way through the 10.5 s gust while the controller has not yet reacted.' },
        { key: 'thrustGust', label: 'Rotor thrust in the operating gust', value: round(Tgust / 1e3, 5), unit: 'kN', explain: 'Thrust at the gust peak with rotor speed and pitch frozen at their rated-wind values, times the dynamic amplification. Often the largest rotor load of all.' },
        { key: 'thrustParked', label: 'Parked-rotor thrust at Ve50', value: round(Tpark(x.Ve50) / 1e3, 5), unit: 'kN', explain: 'Drag of the feathered rotor and nacelle in the 50-year gust from the entered thrust coefficient. Much smaller than operating thrust: storms load the pile through waves, not the rotor.' },
        { key: 'hmax', label: 'Extreme wave height used', value: round(Hext, 4), unit: 'm', explain: `Most probable maximum 1.86·Hs in a 3-hour storm, capped at the breaking height of ${fmt(brk.governing)} m for this depth and period.` },
        { key: 'crestElev', label: 'Extreme crest elevation above mean sea level', value: round(crestElev, 4), unit: 'm', explain: 'Stokes second-order crest of the extreme wave on top of surge and tide. The transition platform and all secondary steel should sit above it with an air gap.' },
        { key: 'airGap', label: 'Air gap under the transition platform', value: round(gap, 4), unit: 'm', limit: inp.gapMin, sense: 'min', explain: 'Platform elevation minus extreme crest elevation. Negative means the wave reaches the platform.' },
        { key: 'ucSurvival', label: 'Survival utilisation of the pile at the mudline', value: round(checks[0].uc, 4), unit: '-', limit: 1, sense: 'max', explain: 'Factored stress over design yield in the governing extreme case: below 1 the primary structure survives the design storm (soil capacity is checked in the foundation module).' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Characteristic mudline moment by load case', categories: res.map((r) => r.id), series: [{ name: 'Rotor thrust', y: res.map((r) => mn(r.Mthr)) }, { name: 'Wind on tower', y: res.map((r) => mn(r.Mwind)) }, { name: 'Wave and current', y: res.map((r) => mn(r.Mwave)) }], y: { label: 'Moment', unit: 'MN·m' },
          explain: 'Contribution of rotor thrust, tower wind and hydrodynamic load to the mudline moment in each case. Operating cases are thrust-driven, storm cases wave-driven: which one governs decides whether rotor or site sets the foundation.' },
        { type: 'bar', title: 'Design (factored) mudline moment and shear', categories: res.map((r) => r.id), series: [{ name: 'Moment (MN·m)', y: res.map((r) => mn(r.Md)) }, { name: 'Shear × 10 (MN)', y: res.map((r) => round(r.Vd / 1e5, 5)) }], y: { label: 'MN·m · MN×10', unit: '-' },
          explain: 'Factored loads per case (1.35 normal, 1.1 abnormal). The tallest moment bar is the governing ULS case for the pile section.' },
        { type: 'line', title: 'Extreme operating gust at rated wind speed', x: { label: 'Time', unit: 's' }, y: { label: 'Hub-height wind speed', unit: 'm/s' }, series: [{ name: 'Wind speed', x: tg.map((v) => round(v, 3)), y: eog.map((v) => round(v, 4)) }],
          explain: 'IEC "Mexican hat" gust: a dip, a sharp rise to the peak at 5.25 s and a second dip. The rise is faster than the pitch system, so thrust overshoots.' },
        { type: 'line', title: 'Extreme wind profile and rotor span', x: { label: 'Wind speed', unit: 'm/s' }, y: { label: 'Height above sea level', unit: 'm' }, series: [{ name: 'Ve50 (3-s gust)', x: zz.map((z) => round(ewm({ Vref: inp.Vref, z, hubHeight: inp.hubHeight }).Ve50, 4)), y: zz.map((z) => round(z, 3)) }, { name: 'V50 (10-min mean)', x: zz.map((z) => round(ewm({ Vref: inp.Vref, z, hubHeight: inp.hubHeight }).V50, 4)), y: zz.map((z) => round(z, 3)) }],
          explain: 'Height profile of the 50-year gust and 10-minute mean wind speeds (power-law exponent 0.11). The tower drag integrates this profile.' },
        { type: 'line', title: `Hydrodynamic load on the pile, ${gov.id}`, x: { label: 'Load per unit height', unit: 'kN/m' }, y: { label: 'Elevation above storm still water', unit: 'm' }, series: [{ name: 'Wave and current', x: gov.wave.f.map((v) => round(v / 1e3, 4)), y: gov.wave.z.map((v) => round(v, 3)) }],
          explain: 'Distributed wave and current force at the phase of maximum overturning moment in the governing case. Load concentrates near the surface, where particle velocity and acceleration are highest.' },
      ],
      tables: [
        { title: 'Load cases', columns: [{ key: 'id', label: 'Case' }, { key: 'name', label: 'Description' }, { key: 'T', label: 'Rotor thrust', unit: 'kN' }, { key: 'Uw', label: 'Hub wind', unit: 'm/s' }, { key: 'H', label: 'Wave height', unit: 'm' }, { key: 'V', label: 'Shear (char.)', unit: 'MN' }, { key: 'M', label: 'Moment (char.)', unit: 'MN·m' }, { key: 'g', label: 'γf', unit: '-' }, { key: 'Md', label: 'Moment (design)', unit: 'MN·m' }, { key: 'Mtb', label: 'Tower-base moment (char.)', unit: 'MN·m' }],
          rows: res.map((r) => ({ id: r.id, name: r.name, T: round(r.T / 1e3, 5), Uw: round(r.Uw, 4), H: round(r.H, 4), V: mn(r.V), M: mn(r.M), g: r.g, Md: mn(r.Md), Mtb: mn(r.Mtb) })),
          explain: 'One row per deterministic load case with its environment, characteristic and factored mudline loads and the tower-base moment. Case numbers follow IEC 61400-3-1 Table 1 in simplified, quasi-static form.' },
        { title: 'Extreme environment', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [{ q: 'Ve50 / Ve1 (3-s gust, hub)', v: `${fmt(x.Ve50)} / ${fmt(x.Ve1)}`, u: 'm/s' }, { q: 'V50 / V1 (10-min mean, hub)', v: `${fmt(x.V50)} / ${fmt(x.V1)}`, u: 'm/s' }, { q: 'Rated wind, NTM σ₁, ETM σ₁', v: `${fmt(Vr)} / ${fmt(sigma1)} / ${fmt(sEtm)}`, u: 'm/s' }, { q: 'EOG magnitude / peak speed', v: `${fmt(vg)} / ${fmt(Vpk)}`, u: 'm/s' },
            { q: 'Extreme / reduced wave height', v: `${fmt(Hext)} / ${fmt(Hred)}`, u: 'm' }, { q: 'Breaking height (Miche / McCowan)', v: `${fmt(brk.miche)} / ${fmt(brk.mccowan)}`, u: 'm' }, { q: 'Wavelength at the extreme period', v: fmt(L50, 4), u: 'm' }, { q: 'Crest above storm still water', v: fmt(crest), u: 'm' }, ...(slamRec ? [{ q: 'Slam force / wave celerity', v: `${fmt(slamRec.F / 1e3)} kN / ${fmt(slamRec.C)}`, u: 'kN / m/s' }] : [])],
          explain: 'Extreme wind and wave parameters behind the load cases, for comparison with the site metocean design basis.' },
      ],
      notes: [
        `Governing design case: ${gov.id} (${gov.name}) with ${fmt(gov.Md / 1e6)} MN·m factored; governing characteristic case: ${govChar.id} with ${fmt(govChar.M / 1e6)} MN·m and ${fmt(govChar.V / 1e6)} MN.`,
        `Edge-on BEM thrust coefficient of the feathered rotor at standstill: ${fmt(edge.Ct, 2)} (lower bound without yaw error or nacelle); ${inp.ctParked} is used.`,
        `Extreme wave ${fmt(Hext)} m (1.86·Hs = ${fmt(1.86 * inp.Hs100)} m, breaking limit ${fmt(brk.governing)} m) in ${fmt(d50)} m of storm water depth; crest ${fmt(crestElev)} m above mean sea level, transition at ${inp.transZ} m.`,
        `Governing check: ${gc.name} (utilisation ${fmt(gc.uc)}). Quasi-static cases: confirm with time-domain simulations and the turbine supplier's load document.`,
      ],
      advice: [
        adv('engineering', sev(gc.uc), `The governing extreme case is ${gov.id}: design mudline moment ${fmt(gov.Md / 1e6)} MN·m, survival utilisation of the pile ${fmt(checks[0].uc)}; air gap ${fmt(gap)} m against ${inp.gapMin} m required.`,
          checks[0].uc > 1 ? `Increase the monopile wall at the mudline to about ${fmt(tReq, 3)} mm or the diameter by ${fmt((Math.sqrt(checks[0].uc) - 1) * 100)} %.` : gap < inp.gapMin ? `Raise the transition platform to at least ${fmt(crestElev + inp.gapMin)} m above mean sea level.` : checks[2].uc > 1 ? 'The extreme wave is depth-limited: switch on the breaking-wave slam case and check the seabed slope for plunging breakers.' : `The structure survives the extreme cases with margin. Pass ${fmt(govChar.M / 1e6)} MN·m and ${fmt(govChar.V / 1e6)} MN to the foundation module; ${SUPPLIER}`,
          'IEC 61400-3-1:2019 Table 1 (DLC 1.3, 1.6, 2.3, 6.1, 6.2) with the load factors of IEC 61400-1:2019 Table 3; DNV-ST-0126 (2021) for the structural check and air gap.',
          `Monopile above the mudline: ${fmt(pileMass)} t. ${steelCo2Text(pileMass)} Designing for the true site extremes instead of the generic class avoids steel that never works.`),
        adv('economic', 'info', `${gov.Mwave > gov.Mthr ? 'Waves dominate' : 'Rotor thrust dominates'} the governing case (${fmt(Math.max(gov.Mwave, gov.Mthr) / gov.M * 100)} % of the moment): ${gov.Mwave > gov.Mthr ? 'the site metocean extremes, not the turbine, size the pile' : 'the turbine load level, not the storm, sizes the pile'}.`,
          gov.Mwave > gov.Mthr ? `A site-specific wave study (directional extremes, depth-limited breaking, joint probability with wind) usually reduces the design wave by 5–15 %; with moment ∝ H to H², that is worth roughly ${fmt(pileMass * 0.07)} t of pile steel per position.` : `Thrust peak shaving in the controller (a 5–10 % lower gust thrust for about 1 % of energy) would lower the design moment by about ${fmt(gov.Mthr / gov.M * 7.5)} % — weigh ${fmt(pileMass * gov.Mthr / gov.M * 0.075)} t of steel per position against the energy.`,
          'Value-of-information practice for metocean campaigns and LCOE-based controller trade-offs (IEA Wind Task 37); Class 4 estimate per AACE RP 18R-97 on steel tonnage.'),
      ],
      publish: { 'owt.mudlineMoment': govChar.M, 'owt.mudlineShear': govChar.V, 'owt.mudlineMomentDesign': gov.Md, 'owt.crestElevation': crestElev },
      quantities: { steel_t: round(pileMass, 4) },
    };
  },
  verify: [
    { name: 'IEC extreme wind and gust models: Ve50 = 1.4·Vref and Vgust = min{1.35(Ve1 − Vhub), 3.3σ₁/(1 + 0.1D/Λ₁)}', source: 'IEC 61400-1:2019 Eqs. (11), (14)–(17): hand calculation for class I, category C, D = 126 m, hub height 90 m',
      // Ve50 = 1.4 × 50 = 70 m/s; Ve1 = 0.8 × 70 = 56 m/s. The schedule's rated wind speed is 11.50 m/s: σ₁ = 0.12 × (0.75 × 11.50 + 5.6) = 1.707 m/s;
      // Λ₁ = 42 m → 3.3 × 1.707 / (1 + 0.1 × 126/42) = 4.333 m/s (the other bound, 1.35 × (56 − 11.5) = 60 m/s, does not govern)
      inputs: { Vref: 50, Iref: 0.12, hubHeight: 90 }, expect: [{ kpi: 've50', value: 70, tol: 1e-6 }, { kpi: 'vgust', value: 4.333, tol: 0.005 }] },
    { name: 'Parked-rotor thrust: T = ½ρ·A·Ct·Ve50² with the coned swept area', source: 'Hand calculation: A = π·(63·cos 2.5°)² = 12 445 m², ½·1.225·12 445·0.014·70² = 522.9 kN',
      inputs: { Vref: 50, ctParked: 0.014 }, expect: [{ kpi: 'thrustParked', value: 522.9, tol: 0.003 }] },
  ],
};
