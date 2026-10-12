// Nearshore wave transformation to the toe of a coastal structure: shoaling, refraction, depth-induced breaking,
// set-up, diffraction behind a breakwater, reflection and wave–current interaction.
import { linearWave, transformProfile, energyFluxBalance, godaSurfZoneHeights, breakingGoda, breakingWeggel, breakingMiche, breakingMcCowan,
  diffractionSemiInfinite, diffractionGap, reflectionSeeligAhrens, waveCurrent, shoalingCoefficient } from '../../engine/coastal/transform.js';
import { iribarren } from '../../engine/coastal/runup.js';
import { selectWaveTheory } from '../../engine/metocean/waves.js';
import { G, DEG, need, run, finite, r, f, interp, linspace, thin, check, eng, eco, slrInputs, slrRows, SLR_BASIS } from './_common.js';

/** Cross-shore bed profile (x shoreward, z bed elevation to datum) from an imported model or the table. */
function bedProfile(inp) {
  const g = inp.bathy;
  let pts = [];
  if (g && g.kind === 'grid') {
    // central row of the grid, taken along +x
    const j = Math.floor(g.ny / 2);
    for (let i = 0; i < g.nx; i++) pts.push([g.x0 + i * g.dx, g.z[j * g.nx + i]]);
  } else if (g && (g.kind === 'polyline' || g.kind === 'points')) {
    const p = g.points || g.positions;
    let s = 0;
    for (let i = 0; i < p.length / 3; i++) {
      if (i > 0) s += Math.hypot(p[3 * i] - p[3 * i - 3], p[3 * i + 1] - p[3 * i - 2]);
      pts.push([s, p[3 * i + 2]]);
    }
  } else {
    pts = (inp.profile || []).map((q) => [Number(q.x), Number(q.z)]);
  }
  pts = pts.filter((q) => Number.isFinite(q[0]) && Number.isFinite(q[1])).sort((a, b) => a[0] - b[0]);
  need(pts.length >= 2, 'The bed profile needs at least two points (distance, bed level).');
  for (let i = 1; i < pts.length; i++) need(pts[i][0] > pts[i - 1][0], 'Bed-profile distances must increase strictly from the offshore boundary to the structure.');
  return { xs: pts.map((q) => q[0]), zs: pts.map((q) => q[1]) };
}

function march(inp, bed, wl, n) {
  const x = linspace(bed.xs[0], bed.xs[bed.xs.length - 1], n);
  const depth = x.map((xi) => wl - interp(bed.xs, bed.zs, xi));
  need(depth[0] > 1, `The offshore boundary is in only ${f(depth[0])} m of water at this water level: move the boundary seaward or check the datum.`);
  need(depth[n - 1] > 0.3, `The structure toe is dry or nearly dry (depth ${f(depth[n - 1])} m) at this water level: there is no wave to transform. Raise the water level or lower the toe.`);
  need(Math.min(...depth) > 0.05, 'The bed profile rises above the water level seaward of the structure (a bar or shoal dries). Use the flood-wave module for emerged profiles.');
  const e = run(() => energyFluxBalance({ x, depth, Hrms0: inp.Hs0 / Math.SQRT2, Tp: inp.Tp, theta0Deg: inp.thetaDeg, gamma: inp.gammaB, fw: inp.fw, rho: inp.rhoW }), 'Energy-flux balance');
  need(e.nWet === n, 'The wave march stopped before the structure because the water became too shallow. Check the profile and water level.');
  return { x, depth, e };
}

export default {
  id: 's07-wave-transformation', suite: 7, title: 'Nearshore wave transformation',
  summary: 'Carries the offshore sea state to the structure toe — shoaling, refraction, breaking, set-up, diffraction, reflection and wave–current interaction.',
  covers: [3], fidelity: 'reduced-order',
  standards: ['CIRIA C683 Rock Manual (2007) §4.2', 'USACE Coastal Engineering Manual EM 1110-2-1100 Part II', 'Goda, Random Seas and Design of Maritime Structures (2000) ch. 3'],
  models: ['Airy linear wave theory', 'Wave shoaling', 'Snell refraction', 'Battjes–Janssen breaking dissipation', 'Radiation-stress set-up', 'Goda surf-zone wave heights', 'Penney–Price diffraction', 'Seelig–Ahrens reflection', 'Wave-action conservation on a current', 'Wave-theory selection (Stokes / cnoidal / stream function)'],
  limits: 'One cross-shore line with straight parallel contours and a narrow-banded sea; diffraction is the constant-depth Sommerfeld solution and the gap solution is a two-tip superposition valid for gaps wider than about one wavelength. Boussinesq, mild-slope, spectral wave-action and 3-D CFD models are outside the module, as are harbour resonance and 2-D refraction over irregular bathymetry.',
  publishes: ['marine.HsToe', 'marine.HmaxToe', 'marine.hToe', 'marine.TpToe', 'marine.setupToe', 'marine.thetaToeDeg'],
  searchTerms: 'nearshore wave transformation breaking set-up diffraction breakwater',
  inputs: [
    { key: 'Hs0', label: 'Significant wave height at the offshore boundary', unit: 'm', type: 'number', default: 5.5, min: 0.05, max: 20, group: 'Offshore wave boundary', role: 'bc', bind: 'site.Hs100', uncertain: { dist: 'lognormal', cov: 0.1 },
      help: 'Hm0 of the design storm at the first point of the bed profile.' },
    { key: 'Tp', label: 'Peak period', unit: 's', type: 'number', default: 11, min: 2, max: 25, group: 'Offshore wave boundary', role: 'bc', bind: 'site.Tp100', uncertain: { dist: 'normal', cov: 0.08 } },
    { key: 'thetaDeg', label: 'Wave angle from shore-normal at the boundary', unit: 'deg', type: 'number', default: 25, min: 0, max: 75, group: 'Offshore wave boundary', role: 'bc' },
    { key: 'bathy', label: 'Bathymetry / shore-normal profile (optional import)', type: 'geometry', kinds: ['grid', 'polyline', 'points'], default: null, group: 'Bathymetry', role: 'geometry',
      help: 'A grid uses its central row along +x; a polyline or sorted sounding line uses chainage and z. Leave empty to use the table.' },
    { key: 'profile', label: 'Bed profile (offshore boundary → structure toe)', type: 'table', group: 'Bathymetry', role: 'ic',
      columns: [{ key: 'x', label: 'Distance', unit: 'm', type: 'number' }, { key: 'z', label: 'Bed level to datum', unit: 'm', type: 'number' }],
      default: [{ x: 0, z: -22 }, { x: 700, z: -14 }, { x: 1200, z: -10 }, { x: 1600, z: -8 }, { x: 1900, z: -7 }],
      help: 'Initial seabed. The last row is the toe of the structure.' },
    { key: 'tide', label: 'Tide level above datum', unit: 'm', type: 'number', default: 1.2, min: -5, max: 10, group: 'Water level', role: 'bc', help: 'Design high water (e.g. MHWS) above the profile datum.' },
    { key: 'surge', label: 'Storm surge', unit: 'm', type: 'number', default: 0.8, min: 0, max: 8, group: 'Water level', role: 'bc', bind: 'site.surge100', uncertain: { dist: 'lognormal', cov: 0.25 } },
    ...slrInputs(),
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Water level', role: 'environment', bind: 'site.rhoW' },
    { key: 'gammaB', label: 'Breaker parameter γ', unit: '-', type: 'number', default: 0.8, min: 0.4, max: 1.2, group: 'Breaking and friction', role: 'criteria', uncertain: { dist: 'normal', cov: 0.1 },
      help: 'Ratio of maximum wave height to depth in the dissipation model (Battjes & Janssen 0.8). A calibration parameter against measured nearshore wave heights.' },
    { key: 'fw', label: 'Wave friction factor', unit: '-', type: 'number', default: 0.01, min: 0, max: 0.3, group: 'Breaking and friction', role: 'bc', help: 'Bed boundary condition: 0 for no bottom friction, 0.01–0.02 for sand, up to 0.1–0.2 over reefs.' },
    { key: 'layout', label: 'Breakwater plan layout', type: 'select', default: 'semi', options: [{ value: 'none', label: 'No sheltering' }, { value: 'semi', label: 'Single breakwater arm (semi-infinite)' }, { value: 'gap', label: 'Entrance gap between two arms' }], group: 'Diffraction', role: 'geometry' },
    { key: 'gapB', label: 'Entrance gap width', unit: 'm', type: 'number', default: 200, min: 10, max: 2000, group: 'Diffraction', role: 'geometry', showIf: { key: 'layout', in: ['gap'] }, design: true,
      help: 'The two-tip superposition needs a gap wider than about one local wavelength.' },
    { key: 'incDeg', label: 'Wave travel direction from the breakwater line', unit: 'deg', type: 'number', default: 90, min: 20, max: 160, group: 'Diffraction', role: 'bc', help: '90° is normal incidence.' },
    { key: 'xd', label: 'Point of interest: distance along the breakwater from the tip (or gap centre)', unit: 'm', type: 'number', default: 120, min: -2000, max: 2000, group: 'Diffraction', role: 'geometry' },
    { key: 'yd', label: 'Point of interest: distance behind the breakwater', unit: 'm', type: 'number', default: 150, min: 1, max: 3000, group: 'Diffraction', role: 'geometry' },
    { key: 'cotA', label: 'Structure front slope cot α', unit: '-', type: 'number', default: 2, min: 0.2, max: 10, group: 'Reflection', role: 'geometry' },
    { key: 'reflA', label: 'Reflection coefficient a', unit: '-', type: 'number', default: 0.6, min: 0.1, max: 1, group: 'Reflection', role: 'material', help: 'Seelig & Ahrens Kr = a ξ²/(b + ξ²): rubble mound a = 0.6, b = 6.6; smooth slope a = 1.0, b = 5.5.' },
    { key: 'reflB', label: 'Reflection coefficient b', unit: '-', type: 'number', default: 6.6, min: 1, max: 20, group: 'Reflection', role: 'material' },
    { key: 'Uc', label: 'Current along the wave direction at the toe (+ following, − opposing)', unit: 'm/s', type: 'number', default: 0.4, min: -3, max: 3, group: 'Current', role: 'bc', help: 'Tidal or river-mouth current; an opposing current steepens the waves and can block them.' },
    { key: 'kLong', label: 'Longshore-current coefficient', unit: '-', type: 'number', default: 1.17, min: 0.5, max: 2, group: 'Current', role: 'criteria', help: 'Mid-surf velocity V = k √(g Hb,rms) sin αb cos αb (Komar & Inman 1970, k = 1.17). Calibrate against measured currents.' },
    { key: 'nx', label: 'Cross-shore grid points', type: 'int', default: 200, min: 20, max: 3200, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nx', label: 'Cross-shore grid points', levels: [50, 100, 200, 400], qoi: 'HsToe' },
  compute(inp) {
    const bed = bedProfile(inp);
    const steep = inp.Hs0 / (G * inp.Tp * inp.Tp / (2 * Math.PI));
    need(steep < 0.08, `The boundary sea state is steeper than waves can be (Hs/L0 = ${f(steep, 2)}; wind seas rarely exceed 0.05–0.06): lengthen the period or reduce the wave height.`);
    const wl = inp.tide + inp.surge + inp.slr;
    const n = Math.round(inp.nx);
    const { x, depth, e } = march(inp, bed, wl, n);
    const last = n - 1;
    const hToe = depth[last];
    const HsToe = finite(e.Hrms[last] * Math.SQRT2, 'Wave height at the toe');
    const thetaToe = e.thetaDeg[last];
    const setupToe = e.eta[last];

    // linear (non-dissipative) reference: deep-water equivalent of the boundary wave, then Ks·Kr along the profile
    const lw0 = linearWave({ T: inp.Tp, h: depth[0] });
    const sinDeep = Math.sin(inp.thetaDeg * DEG) * lw0.C0 / lw0.C;
    need(sinDeep < 0.999, 'The boundary wave angle is too oblique to be consistent with Snell refraction from deep water. Reduce the angle.');
    const thetaDeep = Math.asin(sinDeep) / DEG;
    const lin0 = transformProfile({ H0: 1, T: inp.Tp, theta0Deg: thetaDeep, depths: [depth[0]] });
    const Hdeep = inp.Hs0 / lin0.H[0];
    const lin = transformProfile({ H0: Hdeep, T: inp.Tp, theta0Deg: thetaDeep, depths: depth, gammaB: inp.gammaB });
    const KsToe = shoalingCoefficient({ T: inp.Tp, h: hToe }) / shoalingCoefficient({ T: inp.Tp, h: depth[0] });
    const KrToe = lin.Kr[last] / lin.Kr[0];

    // breaking limits at the toe; local bed slope over the last tenth of the profile
    const iS = Math.max(0, Math.floor(0.9 * last));
    const slope = Math.min(0.2, Math.max(0.002, (depth[iS] - hToe) / (x[last] - x[iS])));
    const Ts = inp.Tp / 1.05; // T1/3 ≈ Tp/1.05 for single-peaked wind seas (Goda 2000)
    const hbGoda = breakingGoda({ h: hToe, T: Ts, slope }).Hb;
    const hbWeggel = run(() => breakingWeggel({ h: hToe, T: Ts, slope }).Hb, 'Weggel breaker height');
    const hbMiche = breakingMiche({ T: Ts, h: hToe }).Hmax;
    const hbMcCowan = breakingMcCowan({ h: hToe });
    const goda = run(() => godaSurfZoneHeights({ H0: Hdeep * lin.Kr[last], T: Ts, h: hToe, slope }), 'Goda surf-zone heights');
    const HmaxToe = Math.min(goda.Hmax, 1.8 * Math.max(HsToe, goda.Hs));
    const lwToe = linearWave({ T: inp.Tp, h: hToe });
    const theory = selectWaveTheory({ H: Math.min(HmaxToe, 0.99 * Math.min(hbMiche, hbMcCowan)), T: Ts, d: hToe });

    // first node where the linear wave is depth-limited: breaker line for the longshore current
    let ib = Array.from(lin.broken).findIndex((b) => b > 0);
    const brokenBeforeToe = ib >= 0;
    if (ib < 0) ib = last;
    const HbRms = lin.H[ib] / Math.SQRT2, thB = lin.thetaDeg[ib] * DEG;
    const Vlong = inp.kLong * Math.sqrt(G * HbRms) * Math.sin(thB) * Math.cos(thB);

    // diffraction at the point of interest and over the lee area
    const L = lwToe.L;
    const kd = (px, py) => (inp.layout === 'semi' ? diffractionSemiInfinite({ x: px, y: py, L, theta0Deg: inp.incDeg }).Kd
      : inp.layout === 'gap' ? diffractionGap({ x: px, y: py, L, B: inp.gapB, theta0Deg: inp.incDeg }).Kd : 1);
    const Kd = finite(kd(inp.xd, inp.yd), 'Diffraction coefficient');
    const xi = iribarren({ tanAlpha: 1 / inp.cotA, H: HsToe, T: inp.Tp }).xi;
    const Krefl = reflectionSeeligAhrens({ xi, a: inp.reflA, b: inp.reflB });
    const wc = run(() => waveCurrent({ T: inp.Tp, h: hToe, U: inp.Uc, H0: HsToe }), 'Wave–current interaction');
    const HsCurrent = wc.blocked ? 0 : wc.H;

    // sea-level-rise pathway: the toe wave of a depth-limited site grows with the water level
    const scen = slrRows(inp).map((s) => {
      const m = march(inp, bed, inp.tide + inp.surge + s.rise, Math.min(n, 200));
      return { name: s.name, rise: s.rise, h: r(m.depth[m.depth.length - 1]), Hs: r(m.e.Hrms[m.depth.length - 1] * Math.SQRT2), setup: r(m.e.eta[m.depth.length - 1], 3) };
    });
    const dHdSlr = scen.length > 1 && scen[scen.length - 1].rise > scen[0].rise ? (scen[scen.length - 1].Hs - scen[0].Hs) / (scen[scen.length - 1].rise - scen[0].rise) : 0;
    const depthLimited = HsToe / hToe > 0.45;

    const [px, pH, pLin, pEta, pQb, pZ] = thin([x, Array.from(e.Hrms, (v) => v * Math.SQRT2), Array.from(lin.H), Array.from(e.eta), Array.from(e.Qb), x.map((xi2) => interp(bed.xs, bed.zs, xi2))], 300);
    const plots = [
      { type: 'line', title: 'Significant wave height across the profile', x: { label: 'Distance from offshore boundary', unit: 'm' }, y: { label: 'Hs', unit: 'm' },
        series: [{ name: 'Energy-flux balance (breaking, friction)', x: px, y: pH.map((v) => r(v)) }, { name: 'Linear Ks·Kr with depth cap', x: px, y: pLin.map((v) => r(v)), dash: true }],
        explain: 'Wave height from the offshore boundary to the structure toe. Where the solid line falls below the linear line the sea is losing energy to breaking and friction: the toe wave is then depth-limited and grows with every rise in water level.' },
      { type: 'line', title: 'Mean water level and breaking fraction', x: { label: 'Distance from offshore boundary', unit: 'm' }, y: { label: 'Set-up (m) / fraction broken (-)', unit: '' },
        series: [{ name: 'Wave set-up η', x: px, y: pEta.map((v) => r(v, 3)) }, { name: 'Fraction of broken waves Qb', x: px, y: pQb.map((v) => r(v, 3)) }],
        explain: 'Set-down seaward of the breaker line turns into set-up through the surf zone; the set-up at the toe adds to the still-water level used for overtopping and wall loads. Qb near 1 means a saturated surf zone.' },
      { type: 'area', title: 'Bed profile and water level', x: { label: 'Distance from offshore boundary', unit: 'm' }, y: { label: 'Level to datum', unit: 'm' },
        series: [{ name: 'Bed', x: px, y: pZ.map((v) => r(v)) }, { name: 'Still-water level incl. surge and sea-level rise', x: [px[0], px[px.length - 1]], y: [wl, wl], dash: true }],
        explain: 'The initial seabed used for the march and the design still-water level. Check that it matches the surveyed shore-normal section at the structure.' },
      { type: 'line', title: 'Toe wave height versus sea-level rise', x: { label: 'Sea-level rise', unit: 'm' }, y: { label: 'Hs at toe', unit: 'm' },
        series: [{ name: 'Hs at toe', x: scen.map((s) => s.rise), y: scen.map((s) => s.Hs) }],
        explain: 'How the design wave at the structure grows as the sea rises. A steep line marks a depth-limited site where armour size, overtopping and wall loads all increase with climate change.' },
    ];
    if (inp.layout !== 'none') {
      const span = Math.max(3 * L, Math.abs(inp.xd) * 1.3, inp.layout === 'gap' ? inp.gapB : 0), reach = Math.max(4 * L, inp.yd * 1.3);
      const xs = linspace(-span, span, 41), ys = linspace(reach / 24, reach, 24);
      plots.push({ type: 'heatmap', title: 'Diffraction coefficient in the lee', xs: xs.map((v) => r(v)), ys: ys.map((v) => r(v)), z: ys.map((yy) => xs.map((xx) => r(kd(xx, yy), 3))), zLabel: 'Kd', zUnit: '-',
        x: { label: inp.layout === 'gap' ? 'Along breakwater from gap centre' : 'Along breakwater from tip', unit: 'm' }, y: { label: 'Distance behind breakwater', unit: 'm' },
        explain: 'Ratio of local to incident wave height behind the breakwater. Values below about 0.3 are well sheltered; lobes above 1 near the shadow line are where berths and small craft see amplified waves.' });
    }

    const checks = [
      check('Toe wave against the Goda breaking limit', HmaxToe, hbGoda, 'm', 'Goda (2000) §3.5; Goda (1975) breaker index', 'Compares the highest wave reaching the toe with the largest wave the local depth can carry. At 1.0 the design wave is breaking on the structure, which brings impulsive loads and makes the design wave a function of water level.'),
    ];
    if (inp.layout === 'gap') checks.push(check('Entrance gap at least one wavelength wide', L, inp.gapB, 'm', 'Penney & Price (1952) §6', 'The gap solution superposes two independent breakwater tips. If the gap is narrower than about one wavelength the tips interact and the coefficient is unreliable; a numerical harbour model is then needed.'));
    if (wc.blocked) checks.push(check('Wave blocking by the opposing current', Math.abs(inp.Uc), Math.abs(wc.Ublock), 'm/s', 'Peregrine (1976)', 'When the opposing current exceeds the blocking speed the waves cannot propagate against it and break at the blocking point, concentrating energy and scour there.'));

    const advice = [];
    advice.push(depthLimited
      ? eng('warn', `The toe wave (Hs = ${f(HsToe)} m in ${f(hToe)} m of water, Hs/h = ${f(HsToe / hToe, 2)}) is depth-limited: it rises by about ${f(dHdSlr, 2)} m for each metre of sea-level rise.`,
        `Design armour, overtopping and wall loads for the end-of-life water level (rise ${f(inp.slr)} m here) and keep a raising or toe-berm option for the higher scenarios.`, SLR_BASIS,
        'A depth-limited structure sized for today\'s water level needs heavier armour or a higher crest later; allowing for the rise now avoids a second mobilisation and its rock and fuel.')
      : eng('ok', `The toe wave (Hs = ${f(HsToe)} m, Hs/h = ${f(HsToe / hToe, 2)}) is not depth-limited; it changes by ${f(dHdSlr, 2)} m per metre of sea-level rise.`,
        'Design on the offshore extreme statistics; confirm the offshore boundary sea state with a hindcast extreme-value analysis and re-check if the foreshore erodes.', 'CIRIA C683 Rock Manual (2007) §4.2.4; Goda (2000) ch. 3',
        'Offshore-governed design waves are insensitive to sea-level rise, so the armour quantity fixed today remains adequate through the scenarios.'));
    advice.push(eng('info', `Recommended wave theory at the toe: ${theory.theory}${theory.theory === 'streamFunction' ? ` (order ${theory.streamFunctionOrder})` : ''} (Ursell number ${f(theory.ursell)}). Set-up at the toe is ${f(setupToe, 2)} m and the longshore current about ${f(Vlong, 2)} m/s.`,
      'Add the set-up to the still-water level in the overtopping and stability modules, and calibrate γ and the friction factor against measured nearshore wave heights before final design.', 'USACE Coastal Engineering Manual Part II ch. 3–4; Battjes & Janssen (1978)'));
    if (inp.layout === 'gap' && inp.gapB < L) advice.push(eng('warn', `The gap (${f(inp.gapB)} m) is narrower than the local wavelength (${f(L)} m): the two-tip diffraction solution is outside its validity.`, 'Use a mild-slope or Boussinesq harbour model for the entrance, or treat the result as indicative only.', 'Penney & Price (1952); PIANC WG 24 harbour-tranquillity practice'));
    advice.push(eco('info', `Every 0.1 m of extra toe wave height raises rock armour mass by about ${f(((1 + 0.1 / HsToe) ** 3 - 1) * 100, 2)} % (M50 ∝ H³).`,
      'Spend on a nearshore wave measurement or a calibrated 2-D wave model where the toe wave is uncertain by more than about 10 %: the saving in armour tonnage normally exceeds the survey cost.', 'Value-of-information reasoning; AACE RP 18R-97 estimate classes for the design stage',
      'Armour mass scales with the cube of the design wave, so a 10 % over-estimate of Hs commits about 33 % more quarried rock and its haulage emissions.'));

    return {
      kpis: [
        { key: 'HsToe', label: 'Significant wave height at the toe', value: r(HsToe), unit: 'm', explain: 'Hm0 reaching the structure after shoaling, refraction, friction and breaking. This is the design wave for armour, overtopping and wall loads.' },
        { key: 'HmaxToe', label: 'Maximum wave height at the toe (Goda)', value: r(HmaxToe), unit: 'm', explain: 'Highest individual wave in the design sea at the toe from Goda\'s surf-zone model, capped at 1.8 Hs. Used for wall pressures and pier forces.' },
        { key: 'hToe', label: 'Water depth at the toe', value: r(hToe), unit: 'm', explain: 'Still-water depth at the structure including tide, surge and the sea-level-rise allowance (set-up not included).' },
        { key: 'HsOverH', label: 'Hs / depth at the toe', value: r(HsToe / hToe, 3), unit: '-', explain: 'Above about 0.45–0.5 the sea state is depth-limited: larger offshore waves do not increase the toe wave, but a higher water level does.' },
        { key: 'KsToe', label: 'Shoaling coefficient boundary → toe', value: r(KsToe), unit: '-', explain: 'Linear change of wave height from the change of group velocity alone. Above 1 the wave grows as it slows in shallow water.' },
        { key: 'KrToe', label: 'Refraction coefficient boundary → toe', value: r(KrToe), unit: '-', explain: 'Reduction of wave height as crests turn towards the depth contours and spread along a longer crest (Snell\'s law). 1.0 for normal incidence.' },
        { key: 'thetaToe', label: 'Wave angle at the toe', value: r(thetaToe), unit: 'deg', explain: 'Angle between the wave direction and the shore-normal at the structure. It enters the obliquity factors for run-up, overtopping and wall pressure.' },
        { key: 'setupToe', label: 'Wave set-up at the toe', value: r(setupToe, 3), unit: 'm', explain: 'Rise of the mean water level driven by breaking (radiation stress). Positive values add to the still-water level at the structure; negative is set-down outside the surf zone.' },
        { key: 'HbToe', label: 'Breaking wave height at the toe (Goda)', value: r(hbGoda), unit: 'm', explain: 'Largest individual wave the toe depth and foreshore slope can sustain. Waves above it have broken seaward.' },
        { key: 'LToe', label: 'Wavelength at the toe', value: r(L), unit: 'm', explain: 'Linear-theory wavelength at the peak period in the toe depth; the length scale for diffraction, berm widths and wall-pressure distribution.' },
        { key: 'Kd', label: 'Diffraction coefficient at the point of interest', value: r(Kd, 3), unit: '-', explain: 'Local over incident wave height behind the breakwater. Multiply the toe wave by Kd to get the wave at the berth or lee structure.' },
        { key: 'Krefl', label: 'Reflection coefficient of the structure', value: r(Krefl, 3), unit: '-', explain: 'Reflected over incident wave height (Seelig & Ahrens). High reflection raises standing-wave heights and toe scour in front of the structure and disturbs navigation.' },
        { key: 'Vlong', label: 'Mid-surf longshore current', value: r(Vlong, 3), unit: 'm/s', explain: 'Wave-driven current along the shore from oblique breaking. It drives littoral drift past the structure; zero for normal incidence.' },
        { key: 'HsCurrent', label: 'Toe wave height on the current', value: r(HsCurrent), unit: 'm', explain: 'Wave height after the Doppler shift and wave-action change on the stated current. An opposing current steepens the waves; zero means the waves are blocked.' },
      ],
      checks, plots,
      tables: [
        { title: 'Breaking limits at the toe', columns: [{ key: 'm', label: 'Criterion' }, { key: 'Hb', label: 'Breaking height', unit: 'm' }, { key: 'idx', label: 'Hb / h', unit: '-' }],
          rows: [{ m: 'Goda (1975)', Hb: r(hbGoda), idx: r(hbGoda / hToe, 3) }, { m: 'Weggel (1972)', Hb: r(hbWeggel), idx: r(hbWeggel / hToe, 3) }, { m: 'Miche (1944)', Hb: r(hbMiche), idx: r(hbMiche / hToe, 3) }, { m: 'McCowan (1894)', Hb: r(hbMcCowan), idx: 0.78 }],
          explain: 'Four independent estimates of the largest wave that can exist in the toe depth. Their spread shows the uncertainty of the breaker index; design on the Goda or Weggel value for sloping foreshores.' },
        { title: 'Sea-level-rise scenarios', columns: [{ key: 'name', label: 'Scenario' }, { key: 'rise', label: 'Rise', unit: 'm' }, { key: 'h', label: 'Toe depth', unit: 'm' }, { key: 'Hs', label: 'Hs at toe', unit: 'm' }, { key: 'setup', label: 'Set-up', unit: 'm' }],
          rows: scen, explain: 'The same storm re-run at each sea-level-rise scenario. Use the row matching the end of the design life for the structure checks and the higher rows to size the adaptation allowance.' },
        { title: 'Wave regime at the toe', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }],
          rows: [{ q: 'Recommended wave theory', v: theory.theory }, { q: 'Ursell number', v: f(theory.ursell) }, { q: 'Relative depth d/L', v: f(theory.dOverL) }, { q: 'Local foreshore slope', v: `1:${f(1 / slope)}` }, { q: 'Iribarren number on the structure', v: f(xi) }, { q: 'Breaking before the toe', v: brokenBeforeToe ? 'yes' : 'no' }, { q: 'Deep-water equivalent Hs', v: `${f(Hdeep)} m` }],
          explain: 'Dimensionless numbers that decide which wave theory and which breaking, run-up and stability formulae apply at the structure.' },
      ],
      notes: [
        `${depthLimited ? 'Depth-limited' : 'Offshore-governed'} design wave: Hs ${f(inp.Hs0)} m at the boundary becomes ${f(HsToe)} m at the toe (linear theory alone would give ${f(lin.H[last])} m).`,
        'Hmax and the breaking limits use T1/3 = Tp/1.05. Nonlinear (Stokes, cnoidal, solitary) kinematics are selected here but applied in the wave-load modules.',
      ],
      advice,
      publish: { 'marine.HsToe': HsToe, 'marine.HmaxToe': HmaxToe, 'marine.hToe': hToe, 'marine.TpToe': inp.Tp, 'marine.setupToe': setupToe, 'marine.thetaToeDeg': thetaToe },
    };
  },
  verify: [
    { name: 'Linear shoaling 50 m → 10 m, T = 10 s, normal incidence, no breaking or friction',
      source: 'Dean & Dalrymple (1991) §4.8: H2/H1 = √(Cg1/Cg2). Hand calculation: Cg(50 m) = 8.550 m/s (L = 151.25 m, n = 0.5653), Cg(10 m) = 8.068 m/s (L = 92.36 m, n = 0.8736) → ratio 1.0294',
      inputs: { Hs0: 0.1, Tp: 10, thetaDeg: 0, profile: [{ x: 0, z: -50 }, { x: 4000, z: -10 }], tide: 0, surge: 0, slr: 0, fw: 0, Uc: 0, layout: 'none' },
      expect: [{ kpi: 'HsToe', value: 0.10294, tol: 0.005 }, { kpi: 'KsToe', value: 1.0294, tol: 0.002 }] },
    { name: 'Snell refraction 30° at 50 m → 10 m, T = 10 s',
      source: 'Snell\'s law sin θ/C = const: C(50) = 15.125 m/s, C(10) = 9.2356 m/s → θ = asin(0.5 × 9.2356/15.125) = 17.78°',
      inputs: { Hs0: 0.1, Tp: 10, thetaDeg: 30, profile: [{ x: 0, z: -50 }, { x: 4000, z: -10 }], tide: 0, surge: 0, slr: 0, fw: 0, Uc: 0, layout: 'none' },
      expect: [{ kpi: 'thetaToe', value: 17.78, tol: 0.005 }] },
  ],
};
