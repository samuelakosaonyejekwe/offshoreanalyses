// Pile–soil interaction of a monopile: beam on nonlinear Winkler foundation under mudline loads, mudline
// displacement and rotation, three-spring foundation stiffness, cyclic degradation and rotation accumulation.
import { solveLateralPile, pileHeadCurves, axialLoadSettlement } from '../../engine/geotech/pile.js';
import { monopileStiffness, mudlineRotationCheck, cyclicRotationAccumulation, flexiblePileStiffnessContinuum } from '../../engine/geotech/monopile.js';
import { pyCurveAt, samplePy } from '../../engine/geotech/py.js';
import { DEG, need, run, fmt, round, adv, check, worst, sev, steelCo2Text, tube, supportInputs, soilInputs, steelInputs, steelOf, soilProfile, soilLabel, LIMIT_PY, LIMIT_CODE, LIMIT_HIFI, CONFIRM } from './_common.js';
import { pyInputs, pileOptions } from './_pile.js';

const geo = supportInputs().filter((i) => ['D', 't', 'Lp'].includes(i.key));

export default {
  id: 's03-pile-soil-interaction', suite: 3, title: 'Pile–soil interaction and foundation stiffness',
  summary: 'Deflection, rotation, bending moment, shear and soil reaction along a laterally loaded monopile, the mudline stiffness matrix for the dynamic models, and the cyclic degradation and accumulated rotation against serviceability limits.',
  covers: [7], fidelity: 'numerical',
  standards: ['DNV-ST-0126 (2021) §7.6 (pile foundations, permanent rotation)', 'EN 1998-5:2004 Annex C (pile-head static stiffness)', 'API RP 2GEO (2011, add. 2014) §8.5 (p–y), §8.4 (t–z, Q–z)', 'ISO 19901-4:2016 §8', 'IEC 61400-3-1:2019 §11 (foundation design)', 'LeBlanc, Houlsby & Byrne (2010) Géotechnique 60(2)'],
  models: ['Euler–Bernoulli beam on nonlinear Winkler foundation (finite elements, P–Δ)', 'API p–y curves: Matlock soft clay, O\'Neill–Murchison sand', 'User-supplied p–y curves (monopile-specific soil reactions)', 'Cyclic p–y degradation and p-/y-multipliers', 'Three-spring foundation stiffness by static condensation', 'Elastic continuum pile-head stiffness (Gazetas; EN 1998-5 Annex C)', 't–z / Q–z axial load transfer', 'LeBlanc strain-accumulation law for cyclic rotation', 'Hetényi elastic foundation beam (verification)'],
  limits: `${LIMIT_PY} Winkler springs are uncoupled between depths; the pile is linear elastic. LeBlanc coefficients Tb and Tc are chart values entered by the user. ${LIMIT_CODE} ${LIMIT_HIFI}`,
  publishes: ['foundation.Kh', 'foundation.Kr', 'foundation.Khr', 'foundation.Kv'],
  searchTerms: 'monopile p-y curves foundation stiffness cyclic rotation accumulation PISA',
  inputs: [
    ...geo,
    ...steelInputs(),
    ...soilInputs(),
    ...pyInputs(),
    { key: 'Es', label: 'Soil Young\'s modulus at one diameter depth (continuum check)', unit: 'MPa', type: 'number', default: 120, min: 2, max: 2000, group: 'Soil reaction model', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.4 },
      help: 'Small-to-working-strain modulus for the elastic-continuum stiffness of EN 1998-5 Annex C, listed beside the p–y result as an independent estimate (modulus increasing with the square root of depth in sand, constant in clay).' },
    { key: 'H', label: 'Horizontal load at the mudline (characteristic)', unit: 'MN', type: 'number', default: 4.9, min: 0.01, max: 100, group: 'Loads at the mudline', role: 'bc', bind: 'owt.mudlineShear', uncertain: { dist: 'lognormal', cov: 0.15 },
      help: 'Wind, wave and current shear transferred by the tower and substructure (boundary condition at the pile head).' },
    { key: 'M', label: 'Overturning moment at the mudline (characteristic)', unit: 'MN·m', type: 'number', default: 280, min: 0, max: 5000, group: 'Loads at the mudline', role: 'bc', bind: 'owt.mudlineMoment', uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'V', label: 'Vertical load at the mudline', unit: 'MN', type: 'number', default: 12, min: 0, max: 300, group: 'Loads at the mudline', role: 'bc', help: 'Weight of turbine, tower and substructure; enters the P–Δ term and the steel check.' },
    { key: 'gammaF', label: 'Load factor for the steel check', unit: '-', type: 'number', default: 1.35, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: `IEC 61400-3-1 normal load cases. ${CONFIRM}` },
    { key: 'gammaM', label: 'Material factor on yield', unit: '-', type: 'number', default: 1.1, min: 1, max: 1.5, group: 'Acceptance criteria', role: 'criteria', help: `DNV-ST-0126 tubular sections. ${CONFIRM}` },
    { key: 'rotLimit', label: 'Total mudline rotation limit', unit: 'deg', type: 'number', default: 0.5, min: 0.05, max: 2, group: 'Acceptance criteria', role: 'criteria', help: 'Commonly 0.5° including installation tolerance; set by the turbine supplier.' },
    { key: 'accLimit', label: 'Accumulated rotation limit', unit: 'deg', type: 'number', default: 0.25, min: 0.02, max: 1, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'slsFrac', label: 'Serviceability load / characteristic extreme load', unit: '-', type: 'number', default: 0.6, min: 0.1, max: 1, group: 'Cyclic loading', role: 'load', help: 'Level of the repeated load used for accumulation (for example the load exceeded 100 times in the life).' },
    { key: 'nCycles', label: 'Equivalent number of load cycles at that level', unit: '-', type: 'number', default: 1e4, min: 1, max: 1e9, group: 'Cyclic loading', role: 'load', help: 'Lifetime storm history converted to an equivalent number of cycles at the serviceability load (strain-superposition rule of LeBlanc et al.); 10³–10⁵ is typical.' },
    { key: 'Tb', label: 'LeBlanc load-magnitude factor Tb', unit: '-', type: 'number', default: 0.06, min: 0, max: 0.5, group: 'Cyclic loading', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.4 },
      help: 'Chart value from LeBlanc et al. (2010) Fig. 9: grows with Mmax/MR and falls with relative density (about 0.03–0.15).' },
    { key: 'Tc', label: 'LeBlanc load-direction factor Tc', unit: '-', type: 'number', default: 1, min: 0, max: 4.5, group: 'Cyclic loading', role: 'soil', help: '1 for one-way loading, up to about 4 for partial two-way loading (Mmin/Mmax ≈ −0.6), 0 for static load.' },
    { key: 'nElem', label: 'Pile elements', type: 'int', default: 40, min: 10, max: 320, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Pile elements', levels: [20, 40, 80, 160], qoi: 'yHead' },
  compute(inp, ctx) {
    // SI: mm → m, MN → N, MPa → Pa.
    const steel = steelOf(inp, ctx);
    const H = inp.H * 1e6, M = inp.M * 1e6, V = inp.V * 1e6;
    need(H > 0 || M > 0, 'Enter a horizontal load or a moment: with no lateral action the pile does not deflect.');
    const profile = inp.pyMode === 'api' ? soilProfile(inp) : null;
    const opt = pileOptions(inp, profile, steel, { axialLoad: V });
    const sec = opt.sec;
    const r = run(() => solveLateralPile({ ...opt, head: { H, M } }), 'Lateral pile analysis at the characteristic load');
    const k = run(() => monopileStiffness({ ...opt, head: { H: H * inp.slsFrac, M: M * inp.slsFrac } }), 'Foundation stiffness');
    const K = k.secant;
    // Cyclic curves only exist for the API mode; elsewhere the static response is the reference.
    let rCyc = null, Kcyc = null, cycNote = null;
    if (inp.pyMode === 'api') {
      // The degraded curves may not carry the full extreme load: the static result then stands alone and a note is given.
      try {
        const rc = solveLateralPile({ ...opt, loading: 'cyclic', head: { H, M }, throwOnFail: false });
        if (rc.convergence.converged) rCyc = rc; else cycNote = 'With cyclic (degraded) p–y curves the pile does not reach equilibrium at the characteristic extreme load: the lateral capacity under cyclic loading is exceeded.';
        Kcyc = monopileStiffness({ ...opt, loading: 'cyclic', head: { H: H * inp.slsFrac, M: M * inp.slsFrac } }).secant;
      } catch { Kcyc = null; cycNote = cycNote || 'The cyclic-curve analysis did not converge at the serviceability load.'; }
    }
    const sls = run(() => solveLateralPile({ ...opt, head: { H: H * inp.slsFrac, M: M * inp.slsFrac } }), 'Lateral pile analysis at the serviceability load');
    const acc = cyclicRotationAccumulation({ thetaStatic: Math.abs(sls.head.theta), N: inp.nCycles, Tb: inp.Tb, Tc: inp.Tc });
    const rot = mudlineRotationCheck({ thetaTotal: Math.abs(r.head.theta) + acc.deltaTheta, thetaAccumulated: acc.deltaTheta, limitTotalDeg: inp.rotLimit, limitAccumulatedDeg: inp.accLimit });
    // Axial stiffness from the t–z model (API mode only); EA/L of the embedded shaft otherwise.
    let Kv = (steel.E * sec.A) / inp.Lp;
    if (profile) { try { Kv = axialLoadSettlement({ profile, L: inp.Lp, D: inp.D, t: inp.t / 1e3, E: steel.E, nElem: Math.min(inp.nElem, 60) }, { steps: 24 }).headStiffness; } catch { /* keep EA/L */ } }
    // Elastic-continuum estimate (flexible pile, Gazetas / EN 1998-5 Annex C) for comparison with the Winkler result.
    let Kc = null;
    try { Kc = flexiblePileStiffnessContinuum({ D: inp.D, EI: steel.E * sec.I, Es: inp.Es * 1e6, profile: inp.soil === 'clay' ? 'homogeneous' : 'parabolic' }); } catch { Kc = null; }
    const Md = inp.gammaF * r.maxMoment, Nd = inp.gammaF * V;
    const sig = Nd / sec.A + Md / sec.W, fyd = steel.fy / inp.gammaM;
    const checks = [
      check('Total mudline rotation (extreme + accumulated)', rot.totalDeg, inp.rotLimit, 'deg', 'DNV-ST-0126 (2021) §7.6.2; turbine supplier tilt limit', 'Guards against excessive tilt of the tower: beyond the limit the turbine may have to be shut down and bearings, yaw system and lift are outside their design envelope.'),
      check('Accumulated permanent rotation', rot.accumulatedDeg, inp.accLimit, 'deg', 'DNV-ST-0126 (2021) §7.6.2; LeBlanc et al. (2010)', 'Guards against ratcheting of the pile under repeated one-way loading. Exceedance means the tilt grows over the years and cannot be corrected after installation.'),
      check('Pile wall stress at maximum moment (factored)', sig / 1e6, fyd / 1e6, 'MPa', 'DNV-ST-0126 (2021) §4.4; EN 1993-1-1 §6.2 (elastic section)', 'Axial plus bending stress at the most stressed pile section below the mudline. Above the design yield the wall yields and local buckling of the can becomes possible.'),
      check('Soil resistance mobilised', r.maxMobilisation, 1, '-', 'API RP 2GEO (2011) §8.5', 'Largest ratio of soil reaction to ultimate resistance along the pile. Near 1.0 over a long length the soil is yielding and deflections grow quickly with load.'),
    ];
    const gov = worst(checks), ok = gov.uc <= 1;
    const mm = (a) => Array.from(a, (v) => v * 1e3), MN = (a) => Array.from(a, (v) => v / 1e6);
    const z = Array.from(r.z);
    // p–y curves for display at 1, 2 and 4 diameters.
    const pyPlot = [];
    for (const zd of [0.5, 1, 2, 4]) {
      const zz = zd * inp.D;
      if (zz >= inp.Lp) continue;
      const cv = opt.pyAt ? opt.pyAt(zz, inp.D) : pyCurveAt(profile, zz, { D: inp.D, loading: inp.loading, pMult: inp.pMult, yMult: inp.yMult });
      const s = typeof cv === 'number' ? { y: [0, 0.1 * inp.D], p: [0, cv * 0.1 * inp.D] } : samplePy(cv, { yMax: 0.1 * inp.D, n: 31 });
      pyPlot.push({ name: `z = ${fmt(zz)} m`, x: mm(s.y), y: Array.from(s.p, (v) => v / 1e3) });
    }
    const hc = run(() => pileHeadCurves({ ...opt, head: { H: 1.5 * H, M: 1.5 * M } }, { factors: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.667, 0.8, 0.9, 1.0] }), 'Head load–displacement curve');
    // Wall needed for the stress check at the same diameter (thin-wall W ≈ π D² t / 4).
    const tReq = (sig / fyd) * inp.t;
    const dMass = (tube(inp.D, Math.min(Math.max(tReq, 20), 0.45 * inp.D * 1e3) / 1e3).A - sec.A) * steel.rho * inp.Lp / 1e3;
    const degr = Kcyc ? Kcyc.KLL / K.KLL : 1;
    const rows = [];
    const step = Math.max(1, Math.floor(z.length / 40));
    for (let i = 0; i < z.length; i += step) rows.push({ z: round(z[i]), y: round(r.y[i] * 1e3), M: round(r.moment[i] / 1e6), Q: round(r.shear[i] / 1e6), p: round(r.p[i] / 1e3), mob: round(r.mobilisation[i]) });
    return {
      kpis: [
        { key: 'yHead', label: 'Mudline displacement', unit: 'mm', value: r.head.y * 1e3, explain: 'Lateral displacement of the pile at the mudline under the characteristic extreme load. Compare with the value assumed in the tower and secondary-steel design; tens of millimetres is typical.' },
        { key: 'thetaHead', label: 'Mudline rotation', unit: 'deg', value: Math.abs(r.head.theta) / DEG, limit: inp.rotLimit, sense: 'max', explain: 'Rotation of the pile at the mudline under the characteristic extreme load, before accumulation. It drives the tilt of the nacelle.' },
        { key: 'thetaAcc', label: 'Accumulated rotation after cycling', unit: 'deg', value: rot.accumulatedDeg, limit: inp.accLimit, sense: 'max', explain: 'Permanent rotation predicted by the LeBlanc power law after the stated number of cycles at the serviceability load. Usually limited to 0.25°.' },
        { key: 'maxMoment', label: 'Maximum bending moment in the pile', unit: 'MN·m', value: r.maxMoment / 1e6, explain: 'Largest bending moment along the pile (characteristic). It occurs a few metres below the mudline where the shear changes sign and sizes the thickest can.' },
        { key: 'zMaxMoment', label: 'Depth of maximum moment', unit: 'm', value: r.depthMaxMoment, explain: 'Depth below the mudline of the peak moment: the thickest can and any circumferential weld to avoid should be placed with this in mind.' },
        { key: 'toeKick', label: 'Toe displacement', unit: 'mm', value: r.toe.y * 1e3, explain: 'Displacement of the pile toe. A negative value ("toe kick") shows the pile rotating as a rigid body; values near zero show the embedded length is fully effective.' },
        { key: 'mobMax', label: 'Peak soil mobilisation', unit: '-', value: r.maxMobilisation, limit: 1, sense: 'max', explain: 'Soil reaction over ultimate soil resistance at the most loaded depth. Values near 1 near the surface are normal; a long yielded zone signals a pile that is too short.' },
        { key: 'Kh', label: 'Lateral stiffness K_LL', unit: 'GN/m', value: K.KLL / 1e9, explain: 'Secant mudline stiffness: horizontal force per unit displacement with rotation held. Feeds the natural-frequency and coupled models.' },
        { key: 'Kr', label: 'Rocking stiffness K_RR', unit: 'GN·m/rad', value: K.KRR / 1e9, explain: 'Secant mudline stiffness: moment per unit rotation with displacement held. The tower frequency is most sensitive to this term.' },
        { key: 'Khr', label: 'Coupling stiffness K_LR', unit: 'MN', value: K.KLR / 1e6, explain: 'Cross-coupling between lateral force and rotation (negative in this sign convention). Omitting it overestimates the foundation stiffness.' },
        { key: 'Kv', label: 'Vertical stiffness', unit: 'GN/m', value: Kv / 1e9, explain: 'Axial pile-head stiffness from the t–z and Q–z load-transfer curves (or EA/L when no profile is used).' },
        { key: 'cycDegr', label: 'Cyclic / static lateral stiffness', unit: '-', value: degr, explain: 'Ratio of the lateral stiffness with cyclic p–y curves to the static value. Below 1 the foundation softens under storm loading and the natural frequency drops.' },
        { key: 'ucSteel', label: 'Pile steel utilisation', unit: '-', value: sig / fyd, limit: 1, sense: 'max', explain: 'Factored axial plus bending stress over design yield at the section of maximum moment.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Pile deflection with depth', x: { label: 'Deflection', unit: 'mm' }, y: { label: 'Depth below mudline', unit: 'm' },
          series: [{ name: 'Static curves', x: mm(r.y), y: z.map((v) => -v) }, ...(rCyc ? [{ name: 'Cyclic curves', x: mm(rCyc.y), y: Array.from(rCyc.z, (v) => -v), dash: true }] : [])],
          explain: 'Deflected shape of the pile. Look for the depth of zero deflection (rotation point) and for the sign of the toe movement: a clear reverse movement at the toe means rigid-body rotation.' },
        { type: 'line', title: 'Bending moment and shear with depth', x: { label: 'Moment (MN·m) / shear (MN)', unit: '' }, y: { label: 'Depth below mudline', unit: 'm' },
          series: [{ name: 'Bending moment, MN·m', x: MN(r.moment), y: z.map((v) => -v) }, { name: 'Shear force ×10, MN', x: Array.from(r.shear, (v) => v / 1e5), y: z.map((v) => -v) }],
          explain: 'Moment and shear distributions below the mudline (shear scaled by ten for legibility). The moment peaks where the shear crosses zero; wall thickness can be reduced below that depth.' },
        { type: 'line', title: 'Soil reaction with depth', x: { label: 'Soil reaction', unit: 'kN/m' }, y: { label: 'Depth below mudline', unit: 'm' },
          series: [{ name: 'Mobilised reaction p', x: Array.from(r.p, (v) => v / 1e3), y: z.map((v) => -v) }, { name: 'Ultimate resistance ±pu', x: Array.from(r.pu, (v, i) => (Number.isFinite(v) ? Math.sign(r.p[i] || 1) * v / 1e3 : null)).filter((v) => v != null).slice(0, z.length), y: z.filter((_, i) => Number.isFinite(r.pu[i])).map((v) => -v), dash: true }],
          explain: 'Soil reaction per metre of pile against the ultimate resistance. Where the curves touch the soil has yielded; the reaction reverses below the rotation point.' },
        { type: 'line', title: 'p–y curves used', x: { label: 'Deflection y', unit: 'mm' }, y: { label: 'Soil reaction p', unit: 'kN/m' }, series: pyPlot,
          explain: 'Soil reaction curves at several depths. Check the initial slope (stiffness) and the plateau (ultimate resistance) against site-calibrated curves.' },
        { type: 'line', title: 'Mudline load–displacement curve', x: { label: 'Mudline displacement', unit: 'mm' }, y: { label: 'Horizontal load', unit: 'MN' },
          series: [{ name: 'Response', x: mm(hc.y), y: MN(hc.H) }, { name: 'Characteristic load', x: [r.head.y * 1e3], y: [inp.H] }],
          explain: 'Push-over response at the mudline with the moment-to-shear ratio kept constant, up to 1.5 times the characteristic load. Curvature shows the onset of soil yielding.' },
      ],
      tables: [
        { title: 'Response along the pile', columns: [{ key: 'z', label: 'Depth', unit: 'm' }, { key: 'y', label: 'Deflection', unit: 'mm' }, { key: 'M', label: 'Moment', unit: 'MN·m' }, { key: 'Q', label: 'Shear', unit: 'MN' }, { key: 'p', label: 'Soil reaction', unit: 'kN/m' }, { key: 'mob', label: 'p / pu', unit: '-' }], rows,
          explain: 'One row per depth station under the characteristic extreme load: use it to place can thickness changes and to compare with strain-gauge data from instrumented piles.' },
        { title: 'Mudline stiffness matrix', columns: [{ key: 'k', label: 'Term' }, { key: 'sec', label: 'Secant' }, { key: 'tan', label: 'Tangent' }, { key: 'cyc', label: 'Secant, cyclic curves' }, { key: 'con', label: 'Elastic continuum' }, { key: 'u', label: 'Unit' }],
          rows: [['K_LL', 'KLL', 1e9, 'GN/m'], ['K_RR', 'KRR', 1e9, 'GN·m/rad'], ['K_LR', 'KLR', 1e6, 'MN']].map(([n, kk, f, u]) => ({ k: n, sec: round(K[kk] / f), tan: round(k.tangent[kk] / f), cyc: Kcyc ? round(Kcyc[kk] / f) : '–', con: Kc ? round(Kc[kk] / f) : '–', u })),
          explain: 'Terms of [H; M] = K·[y; θ] at the serviceability load level. Use the secant values in the modal and coupled models; the cyclic column bounds the softening in a storm, and the elastic-continuum column (EN 1998-5 Annex C, flexible pile) is an independent small-strain estimate that is normally stiffer.' },
      ],
      notes: [`Soil: ${inp.pyMode === 'api' ? soilLabel(inp) + ', API-type curves' : inp.pyMode === 'user' ? 'user-supplied p–y curves' : `linear Winkler modulus ${inp.kLin} MPa`}.`,
        `${gov.name} governs (utilisation ${fmt(gov.uc)}).`, `Peak moment ${fmt(r.maxMoment / 1e6)} MN·m at ${fmt(r.depthMaxMoment)} m below the mudline; rotation point at ${r.depthFirstZeroDeflection ? fmt(r.depthFirstZeroDeflection) + ' m' : 'below the toe'}.`,
        'API-type curves are known to misjudge the stiffness of large-diameter monopiles; calibrate with p-/y-multipliers or user curves before using the stiffness for frequency design.', ...(cycNote ? [cycNote] : [])],
      advice: [
        ok ? adv('engineering', sev(gov.uc), `All pile–soil checks pass; ${gov.name.toLowerCase()} is the closest at ${fmt(gov.uc * 100)} %. Mudline rotation ${fmt(Math.abs(r.head.theta) / DEG)}° plus ${fmt(rot.accumulatedDeg)}° accumulated against ${inp.rotLimit}°.`,
          `Pass K_LL = ${fmt(K.KLL / 1e9)} GN/m, K_RR = ${fmt(K.KRR / 1e9)} GN·m/rad, K_LR = ${fmt(K.KLR / 1e6)} MN to the frequency check, and repeat with the cyclic stiffness (${fmt(degr * 100)} % of static) as the lower bound.`,
          'DNV-ST-0126 (2021) §7.6; API RP 2GEO (2011) §8.5; LeBlanc et al. (2010)',
          `The steel check uses ${fmt(sig / fyd * 100)} % of design yield: a wall of about ${fmt(Math.max(tReq, 20))} mm would reach 100 %, a change of ${fmt(dMass)} t over the embedded length — ${steelCo2Text(dMass)}`)
          : adv('engineering', 'fail', `${gov.name} is exceeded (${fmt(gov.demand)} against ${fmt(gov.capacity)} ${gov.unit}).`,
            gov.name.includes('stress') ? `Increase the wall to about ${fmt(tReq * 1.02)} mm at the depth of maximum moment (${fmt(r.depthMaxMoment)} m), or use a higher steel grade.` : `Increase the embedded length or the diameter: lateral stiffness grows roughly with D²–D³; re-run until rotation is below ${inp.rotLimit}° total and ${inp.accLimit}° accumulated.`,
            'DNV-ST-0126 (2021) §7.6 and §4.4; API RP 2GEO (2011) §8.5',
            `Each extra metre of this pile adds ${fmt(sec.A * steel.rho / 1e3)} t of steel ≈ ${fmt(sec.A * steel.rho / 1e3 * 1.9)} t CO₂e (screening factor 1.9 t CO₂e/t); a diameter increase is more material-efficient for stiffness than added length once the pile is flexible.`),
        adv('economic', ok ? 'info' : 'warn', `The embedded pile weighs ${fmt(sec.A * steel.rho * inp.Lp / 1e3)} t (${fmt(sec.A * steel.rho / 1e3)} t per metre). ${ok ? `Rotation reserve is ${fmt((1 - checks[0].uc) * 100)} %.` : 'The design must grow to meet the checks.'}`,
          ok ? 'Before trimming length, confirm the stiffness against measured frequencies of neighbouring turbines: a soil-reaction calibration (site-specific curves) typically saves more steel than it costs in testing.' : 'Compare a longer pile with a larger diameter on installed tonnes and on hammer capacity; request a site-specific soil-reaction calibration, which often recovers 10–20 % of pile mass.',
          'Value-of-information reasoning for site investigation (ISO 19901-8); quantity-based estimating per AACE RP 18R-97 Class 4'),
      ],
      publish: { 'foundation.Kh': K.KLL, 'foundation.Kr': K.KRR, 'foundation.Khr': K.KLR, 'foundation.Kv': Kv },
      quantities: { steel_t: round(sec.A * steel.rho * inp.Lp / 1e3) },
    };
  },
  verify: [
    { name: 'Long pile on a linear Winkler foundation, shear only: y0 = 2Hβ/k', source: 'Hetényi (1946) Beams on Elastic Foundation, semi-infinite beam: β = (k/4EI)^¼ with EI = 210e9 × π(6⁴ − 5.88⁴)/64 = 1.0367e12 N·m², k = 50 MPa → β = 0.05892 1/m; y0 = 2 × 1e6 × 0.05892 / 50e6 = 2.357 mm; θ0 = 2Hβ²/k = 1.389e-4 rad = 0.007957°',
      inputs: { pyMode: 'linear', kLin: 50, Lp: 80, H: 1, M: 0, V: 0, nElem: 80, D: 6, t: 60, steel: 's355' }, expect: [{ kpi: 'yHead', value: 2.357, tol: 0.01 }, { kpi: 'thetaHead', value: 0.007957, tol: 0.01 }] },
    { name: 'Lateral stiffness of the same pile: K_LL = 4EIβ³ (rotation held)', source: 'Hetényi (1946): head stiffness matrix of a semi-infinite beam on elastic foundation, KHH = 4EIβ³ = 4 × 1.0367e12 × 0.05892³ = 0.8483 GN/m; KMM = 2EIβ = 122.2 GN·m/rad',
      inputs: { pyMode: 'linear', kLin: 50, Lp: 80, H: 1, M: 0, V: 0, nElem: 80, D: 6, t: 60, steel: 's355' }, expect: [{ kpi: 'Kh', value: 0.8483, tol: 0.015 }, { kpi: 'Kr', value: 122.2, tol: 0.015 }] },
  ],
};
