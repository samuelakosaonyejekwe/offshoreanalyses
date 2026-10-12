// Soil–pile interaction of a single dolphin pile: beam on nonlinear Winkler foundation with p–y, t–z and Q–z
// springs; profiles along the pile, pile-head stiffness matrix, static against cyclic response.
import { solveLateralPile, pileHeadStiffness, axialLoadSettlement } from '../../engine/geotech/pile.js';
import { pyCurveAt } from '../../engine/geotech/py.js';
import { DEG, need, run, fmt, round, cap, adv, check, arr, linspace, pileInputs, soilInputs, pileProps, profileFrom, latOpts, steelCo2Text, NOT_HIFI, ECON_BASIS } from './_common.js';

export default {
  id: 's06-soil-pile-interaction', suite: 6, title: 'Soil–pile interaction',
  summary: 'Deflection, rotation, bending moment, shear and soil reaction along a laterally and axially loaded dolphin pile, with p–y curves, pile-head stiffness and cyclic degradation.',
  covers: [5], fidelity: 'numerical',
  standards: ['API RP 2GEO (2011, add. 2014) §8.4–8.5', 'ISO 19901-4:2016 §8', 'DNV-RP-C212 (2019) §5', 'Reese & Van Impe (2011) Single Piles and Pile Groups under Lateral Loading'],
  models: ['Winkler beam on elastic foundation', 'p–y curves: Matlock soft clay, API sand (static and cyclic)', 't–z shaft and Q–z toe load transfer', 'Group p-multiplier', 'Timoshenko shear deformation (optional)', 'P–Δ geometric stiffness under axial load', 'Hetényi linear solution (linear-spring option)'],
  limits: 'Uncoupled Winkler springs: soil reactions at different depths do not interact, and lateral and axial responses are solved separately. Elastic pile. Cyclic response uses the API cyclic curves, not an explicit cycle-by-cycle degradation model. Layering corrections to the p–y curves are not applied. ' + NOT_HIFI,
  publishes: ['foundation.Kh', 'foundation.Kr', 'foundation.Khr', 'foundation.Kv', 'dolphin.pileHeadDisp', 'dolphin.pileHeadRot'],
  searchTerms: 'p-y curves laterally loaded pile beam on nonlinear Winkler foundation',
  inputs: [
    ...pileInputs(),
    { key: 'headFix', label: 'Pile head', type: 'select', default: 'fixed', group: 'Pile', role: 'bc', options: [{ value: 'fixed', label: 'Rotation restrained (cap)' }, { value: 'free', label: 'Free, with applied moment' }] },
    { key: 'timo', label: 'Include shear deformation (Timoshenko)', type: 'bool', default: false, group: 'Pile', role: 'control', help: 'Matters only for short, large-diameter piles (length to diameter below about 10).' },
    ...soilInputs(),
    { key: 'H', label: 'Lateral load at the pile head', unit: 'kN', type: 'number', default: 480, min: 0.1, max: 1e5, group: 'Head loads', role: 'bc', bind: 'dolphin.pileShearMax', uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'M', label: 'Moment at the pile head (free head)', unit: 'kN·m', type: 'number', default: 0, min: -1e6, max: 1e6, group: 'Head loads', role: 'bc', showIf: { key: 'headFix', in: ['free'] }, help: 'Positive for a load acting above the head (adds to the deflection).' },
    { key: 'N', label: 'Axial compression', unit: 'kN', type: 'number', default: 3200, min: 0, max: 2e5, group: 'Head loads', role: 'bc', bind: 'dolphin.pileAxialMax', help: 'Enters the lateral solution through the P–Δ term and loads the t–z model.' },
    { key: 'yLimit', label: 'Allowable deflection at the mudline', unit: 'mm', type: 'number', default: 120, min: 1, max: 2000, group: 'Acceptance', role: 'criteria', help: 'Serviceability limit on the pile at the seabed (about 10 % of the diameter is a common ceiling before the soil is regarded as failed).' },
    { key: 'nElem', label: 'Elements over the embedded length', unit: '-', type: 'int', default: 60, min: 10, max: 640, group: 'Numerical', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Pile elements (embedded length)', levels: [20, 40, 80, 160], qoi: 'yHead' },
  compute(inp, ctx) {
    const P = pileProps(inp, ctx);
    const profile = profileFrom(inp, P.D);
    const fixed = inp.headFix === 'fixed';
    // SI: kN → N, kN·m → N·m
    const head = fixed ? { H: inp.H * 1e3, fixedRotation: true } : { H: inp.H * 1e3, M: inp.M * 1e3 };
    const base = latOpts(inp, P, profile, { axialLoad: inp.N * 1e3, timoshenko: !!inp.timo });
    const s = run(() => solveLateralPile({ ...base, head }), 'Lateral pile analysis');
    const K = run(() => pileHeadStiffness({ ...base, head }), 'Pile-head stiffness');
    const py = (inp.soilModel ?? 'py') === 'py';
    // The other loading type, for the degradation ratio
    let other = null;
    if (py) { try { other = solveLateralPile({ ...base, loading: inp.loading === 'cyclic' ? 'static' : 'cyclic', head }); } catch { other = null; } }
    const ax = run(() => axialLoadSettlement({ profile, L: inp.Lemb, D: P.D, t: P.t, E: P.E, nElem: 40 }, { steps: 40 }), 'Axial pile response');
    const Kv = 1 / (1 / ax.headStiffness + inp.Lfree / (P.E * P.A));
    const iMud = s.z.findIndex((z) => z >= -1e-9), yMud = s.y[Math.max(iMud, 0)];
    const iM = Array.from(s.moment, Math.abs).reduce((b, m, i, a) => (m > a[b] ? i : b), 0);
    const sig = inp.N * 1e3 / P.A + Math.abs(s.maxMoment) / P.Z;
    const wAx = inp.N * 1e3 / Kv;
    const checks = [
      check('Deflection at the mudline', Math.abs(yMud) * 1e3, inp.yLimit, 'mm', 'API RP 2GEO (2011) §8.5 (serviceability, project limit)', 'Lateral movement of the pile at the seabed. Large values mean the upper soil has yielded, a gap opens behind the pile and stiffness degrades with each cycle.'),
      check('Combined stress in the pile', sig / 1e6, P.fy / 1e6, 'MPa', 'ISO 19902:2020 §13.3 (elastic indicator)', 'Axial plus bending stress at the section of maximum moment against yield. The full member check with buckling is in the global strength module.'),
      check('Axial load against t–z capacity', inp.N, ax.capacity / 1e3, 'kN', 'API RP 2GEO (2011) §8.4', 'Head compression against the peak of the load–settlement curve. Close to 1 the pile plunges.'),
    ];
    const w = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = w.uc <= 1;
    const depth = Array.from(s.z, (z) => -z);
    const plots = [
      { type: 'line', title: 'Lateral deflection along the pile', x: { label: 'Deflection', unit: 'mm' }, y: { label: 'Elevation relative to mudline', unit: 'm' },
        series: [{ name: py ? `${inp.loading} p–y` : 'Deflection', x: arr(Array.from(s.y, (v) => v * 1e3)), y: arr(depth) }, ...(other ? [{ name: `${inp.loading === 'cyclic' ? 'static' : 'cyclic'} p–y`, x: arr(Array.from(other.y, (v) => v * 1e3)), y: arr(depth), dash: true }] : [])],
        explain: 'Deflected shape from pile head to toe. Look for the depth of the first zero crossing (active length) and whether the toe still moves — a moving toe means the pile is too short to be fixed.' },
      { type: 'line', title: 'Bending moment and shear along the pile', x: { label: 'Moment (kN·m) / shear (kN)', unit: '' }, y: { label: 'Elevation relative to mudline', unit: 'm' },
        series: [{ name: 'Moment (kN·m)', x: arr(Array.from(s.moment, (v) => v / 1e3)), y: arr(depth) }, { name: 'Shear (kN)', x: arr(Array.from(s.shear, (v) => v / 1e3)), y: arr(depth) }],
        explain: 'Internal forces along the pile. The moment peak (at the cap for a fixed head, a few diameters below the mudline for a free head) sizes the wall; shear changes sign there.' },
      { type: 'line', title: 'Soil reaction and ultimate resistance', x: { label: 'Soil reaction', unit: 'kN/m' }, y: { label: 'Elevation relative to mudline', unit: 'm' },
        series: [{ name: 'Mobilised p', x: arr(Array.from(s.p, (v) => v / 1e3)), y: arr(depth) }, ...(py ? [{ name: 'Ultimate pu', x: arr(Array.from(s.pu, (v) => (Number.isFinite(v) ? v / 1e3 : 0))), y: arr(depth), dash: true }, { name: '−pu', x: arr(Array.from(s.pu, (v) => (Number.isFinite(v) ? -v / 1e3 : 0))), y: arr(depth), dash: true }] : [])],
        explain: 'Soil reaction per metre of pile against its ultimate value. Where the solid line touches the dashed envelope the soil has yielded; the reaction reverses below the pivot depth.' },
      { type: 'line', title: 'Axial load–settlement of the pile (t–z, Q–z)', x: { label: 'Head settlement', unit: 'mm' }, y: { label: 'Head load', unit: 'kN' },
        series: [{ name: 'Pile head', x: ax.w.map((v) => round(v * 1e3)), y: ax.Q.map((v) => round(v / 1e3)) }, { name: 'Tip', x: ax.w.map((v) => round(v * 1e3)), y: ax.tipLoad.map((v) => round(v / 1e3)) }, { name: 'Applied', x: [0, round(ax.w[ax.w.length - 1] * 1e3)], y: [inp.N, inp.N], dash: true }],
        explain: 'Load transfer in the axial direction: shaft friction mobilises within a few millimetres, end bearing needs settlement of several per cent of the diameter. The applied load should sit on the stiff initial branch.' },
    ];
    if (py) {
      const sc = inp.scour ?? 0, zs = [1, 2, 4, 8].map((k) => sc + k * P.D).filter((z) => z < inp.Lemb);
      const ys = linspace(0, 0.25 * P.D, 26);
      plots.push({ type: 'line', title: 'p–y curves at selected depths', x: { label: 'Deflection y', unit: 'mm' }, y: { label: 'Soil resistance p', unit: 'kN/m' },
        series: zs.map((z) => { const c = pyCurveAt(profile, z, { D: P.D, loading: inp.loading, pMult: inp.pMult }); return { name: `z = ${fmt(z)} m`, x: ys.map((y) => round(y * 1e3)), y: ys.map((y) => round(c.p(y) / 1e3)) }; }),
        explain: 'Lateral soil resistance curves used as springs at one, two, four and eight diameters below the scoured surface. Resistance and stiffness grow with depth; cyclic curves plateau lower near the surface.' });
    }
    const ratio = other ? (inp.loading === 'cyclic' ? s.head.y / other.head.y : other.head.y / s.head.y) : 1;
    const advice = [
      adv('engineering', ok ? 'ok' : 'fail', `Head deflection ${fmt(s.head.y * 1e3)} mm, mudline deflection ${fmt(yMud * 1e3)} mm, maximum moment ${fmt(Math.abs(s.maxMoment) / 1e3)} kN·m at ${fmt(s.z[iM])} m; ${w.name.toLowerCase()} governs at ${fmt(w.uc)}.`,
        ok ? `Use the secant head stiffness (KHH ${fmt(K.secant.KHH / 1e6)} MN/m) in the structural model at this load level and the tangent stiffness for small dynamic increments.${other ? (ratio > 1.01 ? ` Cyclic curves give ${fmt(ratio)} times the static head deflection: design for the cyclic case.` : ' Static and cyclic curves coincide at this load level (deflections below the cyclic degradation threshold).') : ''}`
          : (w.name.startsWith('Deflection') ? 'Increase the diameter (stiffness grows with D⁴) or protect against scour; lengthening the pile does not help once the toe is fixed.' : w.name.startsWith('Combined') ? 'Thicken the wall over the zone of maximum moment (a heavier can of 6–8 diameters) rather than over the whole pile.' : 'Lengthen the pile or reduce the axial load.'),
        'API RP 2GEO (2011, add. 2014) §8.4–8.5; ISO 19901-4:2016 §8',
        `A thickened can only around the moment peak avoids carrying the heavy wall over the full length: over ${fmt(inp.Lfree + inp.Lemb)} m, each millimetre of wall is ${steelCo2Text(Math.PI * P.D * 0.001 * (inp.Lfree + inp.Lemb) * 7.85)}`),
      adv('economic', 'info', `The pile works laterally over about ${fmt(s.depthFirstZeroDeflection ?? inp.Lemb)} m of embedment; soil below adds axial capacity only.`,
        'Calibrate the p–y stiffness against a lateral load test or instrumented berthing (use the linear-spring option and the calibration tool): a measured stiffness narrows the deflection prediction and can remove a diameter increment.', `${ECON_BASIS}; value-of-information for pile load testing`),
    ];
    return {
      kpis: [
        { key: 'yHead', label: 'Pile-head deflection', unit: 'mm', value: s.head.y * 1e3, explain: 'Lateral movement of the pile head under the applied load, including the free length above the seabed and soil nonlinearity.' },
        { key: 'thetaHead', label: 'Pile-head rotation', unit: 'deg', value: s.head.theta / DEG, explain: 'Rotation of the pile head (zero for a restrained head). It tilts a mono-pile dolphin and its fender panel.' },
        { key: 'yMud', label: 'Deflection at the mudline', unit: 'mm', value: yMud * 1e3, limit: inp.yLimit, sense: 'max', explain: 'Movement of the pile where it enters the soil; the measure of how hard the upper soil is working.' },
        { key: 'Mmax', label: 'Maximum bending moment', unit: 'kN·m', value: Math.abs(s.maxMoment) / 1e3, explain: 'Largest bending moment along the pile, the quantity that sizes the wall thickness.' },
        { key: 'zMmax', label: 'Depth of maximum moment', unit: 'm', value: s.z[iM], explain: 'Position of the moment peak relative to the mudline (negative above it). A thickened section should be centred here.' },
        { key: 'Vmax', label: 'Maximum shear force', unit: 'kN', value: Math.max(...Array.from(s.shear, Math.abs)) / 1e3, explain: 'Largest shear along the pile; below the mudline it can exceed the head load because the soil reaction reverses.' },
        { key: 'sigma', label: 'Combined axial and bending stress', unit: 'MPa', value: sig / 1e6, limit: P.fy / 1e6, sense: 'max', explain: 'N/A + M/Z at the moment peak against yield: an elastic indicator of the section utilisation.' },
        { key: 'mobil', label: 'Peak soil mobilisation p/pu', unit: '-', value: Number.isFinite(s.maxMobilisation) ? cap(s.maxMobilisation, 1.5) : 0, explain: 'Largest ratio of soil reaction to ultimate resistance. 1.0 means the soil has reached its limit at that depth (zero for linear springs).' },
        { key: 'zZero', label: 'Depth of first zero deflection', unit: 'm', value: s.depthFirstZeroDeflection ?? inp.Lemb, explain: 'Depth of the pivot point. Embedment beyond roughly twice this depth adds no lateral stiffness.' },
        { key: 'KHH', label: 'Head stiffness KHH (secant)', unit: 'MN/m', value: K.secant.KHH / 1e6, explain: 'Lateral force per unit head translation with rotation held, at the working load. Feed to the structural model as a foundation spring.' },
        { key: 'KMM', label: 'Head stiffness KMM (secant)', unit: 'MN·m/rad', value: K.secant.KMM / 1e6, explain: 'Moment per unit head rotation with translation held.' },
        { key: 'KHM', label: 'Head coupling stiffness KHM (secant)', unit: 'MN', value: K.secant.KHM / 1e6, explain: 'Cross-coupling between translation and rotation; negative in the sign convention used, and never negligible for a pile.' },
        { key: 'Kv', label: 'Axial head stiffness', unit: 'MN/m', value: Kv / 1e6, explain: 'Axial load per unit head settlement from the t–z model in series with the free length. It governs how push–pull shares between piles.' },
        { key: 'wAx', label: 'Axial settlement under the applied load', unit: 'mm', value: wAx * 1e3, explain: 'Elastic shortening plus soil settlement of the pile head under the axial load.' },
        { key: 'cycRatio', label: 'Cyclic / static head deflection', unit: '-', value: ratio, explain: 'How much more the pile deflects with the degraded cyclic p–y curves. Values well above 1 show that repeated loading governs serviceability.' },
      ],
      checks, plots,
      tables: [{ title: 'Pile-head stiffness matrices', columns: [{ key: 'm', label: 'Matrix' }, { key: 'KHH', label: 'KHH', unit: 'MN/m' }, { key: 'KHM', label: 'KHM', unit: 'MN' }, { key: 'KMM', label: 'KMM', unit: 'MN·m/rad' }],
        rows: [['Secant at the working load', K.secant], ['Tangent at the working load', K.tangent]].map(([m, k]) => ({ m, KHH: round(k.KHH / 1e6), KHM: round(k.KHM / 1e6), KMM: round(k.KMM / 1e6) })),
        explain: 'Head stiffness [H; M] = K·[y; θ] for coupling to a structural model: secant for static analysis at this load, tangent for vibration and small increments about it.' }],
      notes: [`${w.name} governs (utilisation ${fmt(w.uc)}).`, `Converged in ${s.convergence.iterations} iterations (residual ${s.convergence.residual.toExponential(1)}).`, s.depthFirstZeroDeflection == null ? 'The deflection does not change sign along the pile: it behaves as a short rigid pile.' : `Active length about ${fmt(s.depthFirstZeroDeflection)} m.`],
      advice,
      publish: { 'foundation.Kh': K.secant.KHH, 'foundation.Kr': K.secant.KMM, 'foundation.Khr': K.secant.KHM, 'foundation.Kv': Kv, 'dolphin.pileHeadDisp': s.head.y, 'dolphin.pileHeadRot': s.head.theta },
      quantities: {},
    };
  },
  verify: [
    // Long free-head pile on uniform linear springs, load at the mudline (Hetényi 1946):
    // EI = 210e9·π/64·(1 − 0.96⁴) = 1.5530e9 N·m²; k = 20 MN/m²; β = (k/4EI)^¼ = 0.23821 1/m; βL = 9.5 (long);
    // y0 = 2Hβ/k = 2·5e5·0.23821/2e7 = 11.91 mm; θ0 = 2Hβ²/k = 2.837e-3 rad = 0.16256°; Mmax = 0.3224·H/β = 676.7 kN·m.
    { name: 'Hetényi semi-infinite beam on elastic foundation', source: 'Hetényi (1946) Beams on Elastic Foundation, semi-infinite beam with end load',
      inputs: { D: 1, t: 20, corr: 0, steel: 's355', Lfree: 0, Lemb: 40, scour: 0, soilModel: 'linear', kLin: 20000, pMult: 1, headFix: 'free', H: 500, M: 0, N: 0, nElem: 160, timo: false },
      expect: [{ kpi: 'yHead', value: 11.91, tol: 0.005 }, { kpi: 'thetaHead', value: 0.16256, tol: 0.005 }, { kpi: 'Mmax', value: 676.7, tol: 0.01 }] },
  ],
};
