// Soil–structure interaction of a wharf pile: beam on nonlinear Winkler (p–y) springs under the pile-head loads,
// pile-head stiffness matrix for the deck model, the depth of apparent fixity used by the bent frame, the axial
// load–settlement (t–z, Q–z) response, and the head shear that keeps the seabed movement within its limit.
import { solveLateralPile, pileHeadStiffness, equivalentCantilever, axialLoadSettlement } from '../../engine/geotech/pile.js';
import { pyCurveAt } from '../../engine/geotech/py.js';
import { pMultiplier } from '../../engine/geotech/pilegroup.js';
import { need, run, r, f, arr, check, eng, eco, linspace, maxWhere, thin, pileInputs, pileProps, soilInputs, profileFrom, pileSoilModelInputs, latOpts, steelCo2, CONFIRM, ECON_BASIS, NOT_HIFI } from './_common.js';

export default {
  id: 's05-soil-structure', suite: 5, title: 'Soil–structure interaction of the piles',
  summary: 'Deflection, moment and soil reaction along a wharf pile on nonlinear soil springs, the pile-head stiffness and depth of fixity for the deck frame, the axial load–settlement curve and the permissible head shear.',
  covers: [6], fidelity: 'numerical',
  standards: ['API RP 2GEO (2011, add. 2014) §8.4–8.5', 'ISO 19901-4:2016 §8', 'BS 6349-2:2019 §7', 'ASCE 61-14 §4 (soil–structure interaction of piers and wharves)', 'AASHTO LRFD (2020) §10.7.2.4 (p-multipliers)'],
  models: ['Beam on nonlinear Winkler foundation (finite elements)', 'Lateral soil resistance curves: Matlock soft clay, API sand', 'Axial load transfer: t–z shaft and Q–z toe curves', 'Pile-group interaction by p-multipliers', 'Equivalent cantilever (apparent fixity)', 'P–Δ effect of the axial load'],
  limits: `Single pile with soil reactions uncoupled between depths (Winkler); group shadowing enters through the p-multiplier only. The pile is elastic; the sloping under-wharf bed and the stiffness of rock dykes are not represented beyond the layer table. ${NOT_HIFI} Calibrate against lateral load tests or instrumented berthings where available.`,
  publishes: ['wharf.fixityDepth', 'foundation.Kh', 'foundation.Kr', 'foundation.Khr', 'foundation.Kv', 'wharf.pileHeadDeflection'],
  searchTerms: 'p-y curves pile wharf soil structure interaction apparent fixity',
  inputs: [
    ...pileInputs(), ...soilInputs(), ...pileSoilModelInputs(),
    { key: 'soilModel', label: 'Lateral soil model', type: 'select', default: 'py', group: 'Pile model', role: 'soil', options: [{ value: 'py', label: 'Nonlinear p–y curves from the layer table' }, { value: 'linear', label: 'Linear Winkler springs (constant modulus)' }] },
    { key: 'kLin', label: 'Winkler modulus (force per length per deflection)', unit: 'kPa', type: 'number', default: 20000, min: 100, max: 2e6, group: 'Pile model', role: 'soil', showIf: { key: 'soilModel', in: ['linear'] }, help: 'Use to match a measured stiffness with the calibration tool.' },
    { key: 'rowSpacing', label: 'Row spacing in the load direction (0 = use the p-multiplier entered)', unit: 'm', type: 'number', default: 0, min: 0, max: 20, group: 'Pile model', role: 'geometry', help: 'When given, the p-multiplier is taken from the AASHTO table for the row below.' },
    { key: 'rowNo', label: 'Row position (1 = leading)', unit: '-', type: 'int', default: 1, min: 1, max: 6, group: 'Pile model', role: 'geometry', showIf: { key: 'soilModel', in: ['py'] } },
    { key: 'H', label: 'Shear at the pile head', unit: 'kN', type: 'number', default: 400, min: 1, max: 20000, group: 'Pile-head loads', role: 'bc', bind: 'wharf.pileShearMax', uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'M', label: 'Moment at the pile head (free head only)', unit: 'kN·m', type: 'number', default: 0, min: -1e5, max: 1e5, group: 'Pile-head loads', role: 'bc', showIf: { key: 'headFix', in: ['free'] } },
    { key: 'N', label: 'Axial compression', unit: 'kN', type: 'number', default: 3000, min: 0, max: 60000, group: 'Pile-head loads', role: 'bc', bind: 'wharf.pileAxialMax' },
    { key: 'timo', label: 'Include shear deformation of the pile', type: 'bool', default: true, group: 'Numerics', role: 'control' },
    { key: 'yLimit', label: 'Allowable deflection at the seabed', unit: 'mm', type: 'number', default: 40, min: 2, max: 500, group: 'Acceptance', role: 'criteria', help: `Project serviceability limit; 25–50 mm is common for wharf piles under berthing. ${CONFIRM}` },
  ],
  mesh: { key: 'nElem', label: 'Pile elements (embedded length)', levels: [20, 40, 80, 160], qoi: 'yHead' },
  compute(inp, ctx) {
    const P = pileProps(inp, ctx), profile = profileFrom(inp, P.D), fixed = inp.headFix === 'fixed', py = inp.soilModel !== 'linear';
    let pm = inp.pMult;
    if (py && inp.rowSpacing > 0) { const sD = inp.rowSpacing / P.D; pm = sD >= 5 ? pMultiplier({ spacingOverD: 5, row: Math.round(inp.rowNo) }) : pMultiplier({ spacingOverD: Math.max(sD, 3), row: Math.round(inp.rowNo) }); need(sD >= 2.5, `A row spacing of ${f(sD)} diameters is below the range of the p-multiplier table (3 D): treat the rows as a block.`); }
    const sc = inp.scour ?? 0, extra = { axialLoad: inp.N * 1e3, timoshenko: !!inp.timo, pMult: pm };
    if (!py) { const k = inp.kLin * 1e3; extra.pyAt = (z) => (z >= sc ? k : 0); }
    const base = latOpts(inp, P, profile, extra);
    const headOf = (H) => (fixed ? { H, fixedRotation: true } : { H, M: inp.M * 1e3 });
    const hint = 'Lateral pile analysis';
    const s = run(() => solveLateralPile({ ...base, head: headOf(inp.H * 1e3) }), `${hint} (the head shear may exceed the lateral capacity of the pile in this soil)`);
    const stiffHint = 'The pile-head stiffness is not positive at this load: the axial load is at or above the buckling load of the pile in this soil, or the pile is too slender or too short to be restrained. Reduce the axial load or the free length, or increase the diameter or penetration.';
    let K; try { K = pileHeadStiffness({ ...base, head: headOf(inp.H * 1e3) }); } catch { throw new Error(stiffHint); }
    need(K.secant.KHH > 0 && K.secant.KMM > 0 && K.secant.KHH * K.secant.KMM - K.secant.KHM ** 2 > 0, stiffHint);
    let other = null;
    if (py) { try { other = solveLateralPile({ ...base, loading: inp.loading === 'cyclic' ? 'static' : 'cyclic', head: headOf(inp.H * 1e3) }); } catch { other = null; } }
    const eq = equivalentCantilever({ KHH: K.secant.KHH, KHM: K.secant.KHM, KMM: K.secant.KMM, EI: P.EI });
    // equivalent cantilever measured from the pile head; the depth of fixity is what remains below the scoured bed
    const Le = fixed ? eq.fixedHead : eq.lateral, fixDepth = Le - inp.Lfree;
    need(Number.isFinite(fixDepth) && fixDepth > -0.5, 'The equivalent cantilever is shorter than the free length: the head stiffness is not usable at this load. Reduce the head shear or check the soil data.');
    const ax = run(() => axialLoadSettlement({ profile, L: inp.Lemb, D: P.D, t: P.t, E: P.E, nElem: 40 }, { steps: 40 }), 'Axial pile response');
    const Kv = 1 / (1 / ax.headStiffness + inp.Lfree / (P.E * P.A));
    const wAx = (() => { const N = inp.N * 1e3; for (let i = 1; i < ax.Q.length; i++) if (ax.Q[i] >= N) return ax.w[i - 1] + (ax.w[i] - ax.w[i - 1]) * (N - ax.Q[i - 1]) / (ax.Q[i] - ax.Q[i - 1]); return ax.w[ax.w.length - 1]; })() + inp.N * 1e3 * inp.Lfree / (P.E * P.A);
    const iMud = Math.max(s.z.findIndex((z) => z >= sc - 1e-9), 0), yMud = Math.abs(s.y[iMud]);
    const sig = inp.N * 1e3 / P.A + Math.abs(s.maxMoment) / P.Z;
    // head shear at which the seabed deflection reaches its limit
    const yAt = (H) => { try { const q = solveLateralPile({ ...base, head: headOf(H) }); return Math.abs(q.y[iMud]); } catch { return Infinity; } };
    const Hperm = maxWhere((H) => yAt(H) <= inp.yLimit / 1e3, 0.05 * inp.H * 1e3, 6 * inp.H * 1e3, 14);
    const checks = [
      check('Deflection at the seabed', yMud * 1e3, inp.yLimit, 'mm', 'API RP 2GEO (2011) §8.5 (serviceability, project limit)', 'Lateral movement of the pile at the scoured bed level. Large values mean the upper soil has yielded: a gap opens behind the pile and stiffness degrades with each berthing.'),
      check('Combined stress in the pile', sig / 1e6, P.fy / 1e6, 'MPa', 'ISO 19902:2020 §13.3 (elastic indicator)', 'Axial plus bending stress at the section of maximum moment against yield. The full member check with buckling is in the strength module.'),
      check('Axial load against the t–z capacity', inp.N, ax.capacity / 1e3, 'kN', 'API RP 2GEO (2011) §8.4', 'Head compression against the peak of the load–settlement curve. Near 1 the pile plunges; the factored check is in the capacity module.'),
      check('Soil mobilisation (largest p / pu)', Number.isFinite(s.maxMobilisation) ? s.maxMobilisation : 0, 1, '-', 'API RP 2GEO (2011) §8.5', 'Share of the ultimate lateral soil resistance in use at the most loaded depth. At 1 the soil there has yielded and load sheds deeper down the pile.'),
    ];
    const w = checks.slice(0, 3).reduce((a, c) => (c.uc > a.uc ? c : a)), ok = w.uc <= 1;
    const depth = Array.from(s.z, (z) => -z);
    const [dz, yy, mm, vv, pp, pu] = thin([depth, Array.from(s.y, (v) => v * 1e3), Array.from(s.moment, (v) => v / 1e3), Array.from(s.shear, (v) => v / 1e3), Array.from(s.p, (v) => v / 1e3), Array.from(s.pu, (v) => (Number.isFinite(v) ? v / 1e3 : 0))], 200);
    const plots = [
      { type: 'line', title: 'Lateral deflection along the pile', x: { label: 'Deflection', unit: 'mm' }, y: { label: 'Elevation relative to the seabed', unit: 'm' },
        series: [{ name: py ? `${inp.loading} p–y` : 'Linear springs', x: arr(yy), y: arr(dz) }, ...(other ? [{ name: `${inp.loading === 'cyclic' ? 'static' : 'cyclic'} p–y`, x: arr(thin([Array.from(other.y, (v) => v * 1e3)], 200)[0]), y: arr(dz), dash: true }] : [])],
        explain: 'Deflected shape from pile head to toe. Look for the depth of the first zero crossing (the active length) and whether the toe still moves: a moving toe means the pile is too short to be fixed.' },
      { type: 'line', title: 'Bending moment and shear along the pile', x: { label: 'Moment (kN·m) / shear (kN)', unit: '' }, y: { label: 'Elevation relative to the seabed', unit: 'm' },
        series: [{ name: 'Moment (kN·m)', x: arr(mm), y: arr(dz) }, { name: 'Shear (kN)', x: arr(vv), y: arr(dz) }],
        explain: 'Internal forces along the pile. A moment-connected pile has its peak at the deck and a second, opposite peak a few diameters below the bed; that in-ground peak cannot be inspected and sizes the wall.' },
      { type: 'line', title: 'Soil reaction and ultimate resistance', x: { label: 'Soil reaction', unit: 'kN/m' }, y: { label: 'Elevation relative to the seabed', unit: 'm' },
        series: [{ name: 'Mobilised p', x: arr(pp), y: arr(dz) }, ...(py ? [{ name: 'Ultimate pu', x: arr(pu), y: arr(dz), dash: true }, { name: '−pu', x: arr(pu.map((v) => -v)), y: arr(dz), dash: true }] : [])],
        explain: 'Soil reaction per metre of pile against its ultimate value. Where the solid line touches the dashed envelope the soil has yielded; the reaction reverses below the pivot depth.' },
      { type: 'line', title: 'Axial load–settlement of the pile (t–z, Q–z)', x: { label: 'Head settlement', unit: 'mm' }, y: { label: 'Head load', unit: 'kN' },
        series: [{ name: 'Pile head', x: ax.w.map((v) => r(v * 1e3)), y: ax.Q.map((v) => r(v / 1e3)) }, { name: 'Toe', x: ax.w.map((v) => r(v * 1e3)), y: ax.tipLoad.map((v) => r(v / 1e3)) }, { name: 'Applied', x: [0, r(ax.w[ax.w.length - 1] * 1e3)], y: [inp.N, inp.N], dash: true }],
        explain: 'Load transfer in the axial direction: shaft friction mobilises within a few millimetres, end bearing needs settlement of several per cent of the diameter. The applied load should sit on the stiff initial branch.' },
    ];
    if (py) {
      const zs = [1, 2, 4, 8].map((k) => sc + k * P.D).filter((z) => z < inp.Lemb), ys = linspace(0, 0.25 * P.D, 26);
      plots.push({ type: 'line', title: 'p–y curves at selected depths', x: { label: 'Deflection y', unit: 'mm' }, y: { label: 'Soil resistance p', unit: 'kN/m' },
        series: zs.map((z) => { const c = pyCurveAt(profile, z, { D: P.D, loading: inp.loading, pMult: pm }); return { name: `z = ${f(z)} m`, x: ys.map((y) => r(y * 1e3)), y: ys.map((y) => r(c.p(y) / 1e3)) }; }),
        explain: 'Lateral soil springs at one, two, four and eight diameters below the scoured bed. Resistance and stiffness grow with depth; cyclic curves plateau lower near the surface.' });
    }
    const ratio = other ? (inp.loading === 'cyclic' ? s.head.y / other.head.y : other.head.y / s.head.y) : 1;
    return {
      kpis: [
        { key: 'yHead', label: 'Pile-head deflection', value: r(s.head.y * 1e3, 7), unit: 'mm', explain: 'Lateral movement of the deck at this pile under the head shear. It is the single-pile estimate of deck sway and grows faster than the load once the soil yields.' },
        { key: 'thetaHead', label: 'Pile-head rotation', value: r(s.head.theta * 180 / Math.PI, 5), unit: 'deg', explain: 'Rotation at the head; zero for a moment-connected head by construction.' },
        { key: 'yMud', label: 'Deflection at the seabed', value: r(yMud * 1e3), unit: 'mm', limit: inp.yLimit, sense: 'max', explain: 'Movement of the pile where it enters the soil, against the serviceability limit. It controls gapping and cyclic softening.' },
        { key: 'Mmax', label: 'Largest bending moment', value: r(Math.abs(s.maxMoment) / 1e3), unit: 'kN·m', explain: `Peak moment along the pile, at ${f(s.depthMaxMoment)} m relative to the seabed (negative is above the bed).` },
        { key: 'Mhead', label: 'Moment at the pile head', value: r(Math.abs(s.head.M) / 1e3), unit: 'kN·m', explain: 'Moment passed between pile and deck beam: the demand on the pile-head connection.' },
        { key: 'KHH', label: 'Lateral head stiffness (secant)', value: r(K.secant.KHH / 1e6), unit: 'MN/m', explain: 'Force per unit translation of the pile head without rotation, at this load level. Use it as a spring in the deck model in place of a fixity depth.' },
        { key: 'KMM', label: 'Rotational head stiffness (secant)', value: r(K.secant.KMM / 1e6), unit: 'MN·m/rad', explain: 'Moment per unit rotation of the pile head without translation.' },
        { key: 'Kv', label: 'Axial head stiffness', value: r(Kv / 1e6), unit: 'MN/m', explain: 'Vertical stiffness of the pile from the t–z analysis in series with the free length; governs how the deck shares vertical load between piles.' },
        { key: 'fixDepth', label: 'Depth of apparent fixity below the seabed', value: r(fixDepth), unit: 'm', explain: `Depth at which a rigidly clamped pile has the same head stiffness (${f(fixDepth / P.D)} diameters). It is the value the bent frame analysis needs, and it lengthens as the load rises.` },
        { key: 'active', label: 'Active length below the seabed', value: r(s.depthFirstZeroDeflection ?? inp.Lemb), unit: 'm', explain: 'Depth of the first zero crossing of the deflection. Soil below adds axial capacity but almost no lateral stiffness.' },
        { key: 'wAx', label: 'Pile-head settlement under the axial load', value: r(wAx * 1e3), unit: 'mm', explain: 'Vertical movement of the pile head from the load-transfer curve plus elastic shortening of the free length.' },
        { key: 'Hperm', label: 'Maximum permissible head shear for the seabed limit', value: r(Hperm / 1e3), unit: 'kN', explain: 'Pile-head shear at which the seabed deflection reaches its limit: the serviceability capacity of the pile for berthing and mooring loads (capped at six times the applied shear).' },
      ],
      checks, plots,
      tables: [{ title: 'Pile-head stiffness matrix', columns: [{ key: 'term', label: 'Term' }, { key: 'sec', label: 'Secant' }, { key: 'tan', label: 'Tangent' }, { key: 'unit', label: 'Unit' }],
        rows: [{ term: 'KHH', sec: r(K.secant.KHH / 1e6), tan: r(K.tangent.KHH / 1e6), unit: 'MN/m' }, { term: 'KHM', sec: r(K.secant.KHM / 1e6), tan: r(K.tangent.KHM / 1e6), unit: 'MN' }, { term: 'KMM', sec: r(K.secant.KMM / 1e6), tan: r(K.tangent.KMM / 1e6), unit: 'MN·m/rad' }, { term: 'Kv', sec: r(Kv / 1e6), tan: r(Kv / 1e6), unit: 'MN/m' }],
        explain: 'Stiffness of the pile and soil below the head, [H; M] = K·[y; θ]. Use the secant values for the design load level and the tangent values for small dynamic increments about it.' }],
      notes: [`${py ? `p–y curves (${inp.loading}) with p-multiplier ${f(pm)}` : `Linear springs of ${f(inp.kLin)} kPa`}; equivalent cantilever ${f(Le)} m from the head.${other ? ` Cyclic over static head deflection ${f(ratio)}.` : ''}`,
        `${w.name} governs at ${f(w.uc)}; largest soil mobilisation ${f(Number.isFinite(s.maxMobilisation) ? s.maxMobilisation : 0)}.`],
      advice: [
        eng(ok ? 'ok' : 'fail', `Head deflection ${f(s.head.y * 1e3)} mm, seabed deflection ${f(yMud * 1e3)} mm against ${f(inp.yLimit)} mm, largest moment ${f(Math.abs(s.maxMoment) / 1e3)} kN·m; apparent fixity ${f(fixDepth)} m (${f(fixDepth / P.D)} D) below the bed.`,
          ok ? `Enter ${f(fixDepth)} m as the fixity depth of the bent frame (or the secant stiffness ${f(K.secant.KHH / 1e6)} MN/m as a spring) and re-run the strength module; the pile can take up to ${f(Hperm / 1e3)} kN of head shear within the seabed limit.`
            : w.name.startsWith('Deflection') ? `Limit the head shear to ${f(Hperm / 1e3)} kN, increase the diameter (stiffness grows with D⁴), add raker piles, or protect the bed against scour.` : w.name.startsWith('Combined') ? 'Thicken the wall over the zone of maximum moment (a heavier section of 6–8 diameters) rather than over the whole pile.' : 'Lengthen the pile or reduce the axial load.',
          'API RP 2GEO (2011, add. 2014) §8.4–8.5; ASCE 61-14 §4',
          `A thickened section only around the in-ground moment peak avoids carrying the heavy wall over the full length: over ${f(inp.Lfree + inp.Lemb)} m each millimetre of wall is ${steelCo2(Math.PI * P.D * 0.001 * (inp.Lfree + inp.Lemb) * 7.85)}.`),
        eco('info', `The pile works laterally over about ${f(s.depthFirstZeroDeflection ?? inp.Lemb)} m of its ${f(inp.Lemb)} m embedment; the fixity depth changes the deck sway and pile moments roughly with its cube.`,
          'Calibrate the soil stiffness against a lateral load test or instrumented berthings (linear-spring option with the calibration tool): a measured stiffness narrows the sway prediction and can remove a diameter increment across every bent.', `${ECON_BASIS}; value-of-information for pile load testing`),
      ],
      publish: { 'wharf.fixityDepth': Math.max(fixDepth, 0), 'foundation.Kh': K.secant.KHH, 'foundation.Kr': K.secant.KMM, 'foundation.Khr': K.secant.KHM, 'foundation.Kv': Kv, 'wharf.pileHeadDeflection': s.head.y },
    };
  },
  verify: [
    // Long free-head pile on uniform linear springs, load at the seabed (Hetényi 1946):
    // EI = 210e9·π/64·(1 − 0.96⁴) = 1.5530e9 N·m²; k = 20 MN/m²; β = (k/4EI)^¼ = 0.23821 1/m; βL = 9.5 (long);
    // y0 = 2Hβ/k = 2·5e5·0.23821/2e7 = 11.91 mm; θ0 = 2Hβ²/k = 2.837e-3 rad = 0.16256°; Mmax = 0.3224·H/β = 676.7 kN·m.
    { name: 'Hetényi semi-infinite beam on elastic foundation', source: 'Hetényi (1946) Beams on Elastic Foundation, semi-infinite beam with end load',
      inputs: { D: 1, t: 20, corr: 0, steel: 's355', Lfree: 0, Lemb: 40, scour: 0, soilModel: 'linear', kLin: 20000, pMult: 1, headFix: 'free', H: 500, M: 0, N: 0, nElem: 160, timo: false },
      expect: [{ kpi: 'yHead', value: 11.91, tol: 0.005 }, { kpi: 'thetaHead', value: 0.16256, tol: 0.005 }, { kpi: 'Mmax', value: 676.7, tol: 0.01 }] },
    // Same pile, head stiffness of the semi-infinite beam: KHH = 4EIβ³ = 4·1.5530e9·0.23821³ = 83.96 MN/m; fixed-head
    // equivalent cantilever Le = (12EI/KHH)^⅓ = (12·1.5530e9/83.96e6)^⅓ = 6.055 m = depth of fixity with no free length.
    { name: 'Head stiffness KHH = 4EIβ³ and the equivalent fixed-head cantilever', source: 'Hetényi (1946); cantilever stiffness 12EI/L³ (first principles)',
      inputs: { D: 1, t: 20, corr: 0, steel: 's355', Lfree: 0, Lemb: 40, scour: 0, soilModel: 'linear', kLin: 20000, pMult: 1, headFix: 'fixed', H: 500, N: 0, nElem: 160, timo: false },
      expect: [{ kpi: 'KHH', value: 83.96, tol: 0.01 }, { kpi: 'fixDepth', value: 6.055, tol: 0.01 }] },
  ],
};
