// Internal tube pressure and integrity: thick-wall stresses, pressure-containment capacity (ASME B31.3 and
// DNV-ST-F101 forms), proof test, hose rating, pressure transients (Joukowsky surge), pressure drop, hydraulic
// response time and pressure-cycle fatigue of a steel umbilical tube.
import { lameStresses, thickWallLimitPressures, b313PressureDesign, dnvBurstPressure, proofTestPressure, hoseRating } from '../../engine/xsection/tubes.js';
import { tubePressureDrop, hydraulicResponseTime } from '../../engine/xsection/umbilical.js';
import { snCurve, cyclesToFailure, SN_CURVE_NAMES } from '../../engine/fatigue/sn.js';
import { adv, check, worst, sev, fmt, round, run, need, linspace, G, CO2_STEEL, STD_UMB, CONFIRM } from './_common.js';

export default {
  id: 's12-tube-pressure-integrity', suite: 12, title: 'Internal tube pressure and integrity',
  summary: 'Pressure-containment utilisation of a hydraulic or chemical tube at design, test and surge pressure, with flow pressure drop, response time and pressure-cycle fatigue.',
  covers: [13], fidelity: 'analytical',
  standards: ['ASME B31.3 (2022) §304.1.2 Eq. (3a) and §302.3.2 (allowable stress basis)', 'DNV-ST-F101 (2021) §5.4.2 Eq. 5.8 (pressure containment)', `${STD_UMB} §7 (tube design, proof test, hose burst ratio)`, 'Timoshenko & Goodier §28 (Lamé); Joukowsky (1898) surge'],
  models: ['Thin- and thick-wall pressure vessel relationships', 'Lamé elasticity with von Mises criterion', 'Pressure containment (burst) models', 'Joukowsky pressure transient with Korteweg wave speed', 'Darcy–Weisbach pressure drop and laminar line-diffusion response time', 'S–N fatigue under internal pressure cycles'],
  limits: 'The ISO 13628-5 / API 17E tube burst formula itself is not implemented: containment is checked with the ASME B31.3 and DNV-ST-F101 forms and with Lamé stresses against overridable utilisation factors. Proof factor 1.5, hose burst ratio 4 and utilisation factors are defaults to be confirmed against the governing edition. Surge is the instantaneous-closure bound (no line packing, friction attenuation or valve characteristic); a method-of-characteristics transient analysis is outside this module.',
  publishes: ['umbilical.tubeUtilisation', 'umbilical.tubeBurstPressure', 'umbilical.hydraulicResponseTime'],
  searchTerms: 'umbilical steel tube pressure containment burst proof test response time',
  inputs: [
    { key: 'od', label: 'Tube outside diameter', unit: 'mm', type: 'number', default: 16.5, min: 3, max: 120, group: 'Tube', role: 'geometry' },
    { key: 'wall', label: 'Nominal wall thickness', unit: 'mm', type: 'number', default: 1.9, min: 0.3, max: 20, group: 'Tube', role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'tol', label: 'Wall under-tolerance', unit: '%', type: 'number', default: 10, min: 0, max: 20, group: 'Tube', role: 'geometry', help: 'Manufacturing tolerance deducted for the code checks (typically ±10 % for seam-welded tube).' },
    { key: 'mat', label: 'Tube material', type: 'select', source: 'materials:steel', default: 'sdss-25cr', group: 'Tube', role: 'material' },
    { key: 'weld', label: 'Weld joint quality factor E', unit: '-', type: 'number', default: 1.0, min: 0.6, max: 1, group: 'Tube', role: 'material', help: 'ASME B31.3 Table A-1B: 1.0 for seamless or fully radiographed longitudinal welds, 0.8–0.85 otherwise.' },
    { key: 'length', label: 'Tube length (umbilical length)', unit: 'm', type: 'number', default: 5000, min: 10, max: 3e5, group: 'Tube', role: 'geometry', bind: 'cable.length' },
    { key: 'pd', label: 'Design pressure at the topside end', unit: 'bar', type: 'number', default: 517, min: 1, max: 2000, group: 'Pressure conditions', role: 'bc', uncertain: { dist: 'normal', cov: 0.03 }, help: 'Hydraulic boundary condition: maximum operating pressure at the supply (517 bar = 7500 psi).' },
    { key: 'depth', label: 'Water depth at the subsea end', unit: 'm', type: 'number', default: 350, min: 0, max: 4000, group: 'Pressure conditions', role: 'bc', bind: 'site.depth' },
    { key: 'rhoF', label: 'Fluid density', unit: 'kg/m³', type: 'number', default: 1060, min: 500, max: 2000, group: 'Fluid', role: 'material', help: 'Water-glycol control fluid ≈ 1060, methanol 792, MEG 1113.' },
    { key: 'visc', label: 'Fluid viscosity', unit: 'mPa·s', type: 'number', default: 3, min: 0.2, max: 500, group: 'Fluid', role: 'material', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'At seabed temperature (viscosity rises steeply in cold water).' },
    { key: 'bulk', label: 'Fluid bulk modulus', unit: 'GPa', type: 'number', default: 2.2, min: 0.5, max: 4, group: 'Fluid', role: 'material' },
    { key: 'flow', label: 'Flow rate', unit: 'L/min', type: 'number', default: 10, min: 0, max: 2000, group: 'Fluid', role: 'bc' },
    { key: 'rough', label: 'Bore roughness', unit: 'mm', type: 'number', default: 0.005, min: 0, max: 0.5, group: 'Fluid', role: 'material' },
    { key: 'cyclesYr', label: 'Full pressure cycles per year', unit: '-', type: 'number', default: 300, min: 0, max: 1e6, group: 'Service conditions', role: 'load', help: 'Depressurise / repressurise cycles of the line (valve operations, shutdowns).' },
    { key: 'cycleFrac', label: 'Cycle range / design pressure', unit: '-', type: 'number', default: 1, min: 0, max: 1, group: 'Service conditions', role: 'load' },
    { key: 'sn', label: 'S–N curve (tube seam weld)', type: 'select', default: 'C1', options: SN_CURVE_NAMES.map((n) => ({ value: n, label: `DNV-RP-C203 ${n}` })), group: 'Service conditions', role: 'material' },
    { key: 'life', label: 'Design service life', unit: 'years', type: 'number', default: 25, min: 1, max: 60, group: 'Service conditions', role: 'criteria', bind: 'life.designYears' },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 10, min: 1, max: 20, group: 'Service conditions', role: 'criteria' },
    { key: 'hoseBurst', label: 'Thermoplastic hose burst pressure (0 = no hose)', unit: 'bar', type: 'number', default: 1500, min: 0, max: 6000, group: 'Hose (if used)', role: 'material' },
    { key: 'hosePd', label: 'Hose design pressure', unit: 'bar', type: 'number', default: 345, min: 0, max: 2000, group: 'Hose (if used)', role: 'bc' },
    { key: 'hoseVe', label: 'Hose volumetric expansion', unit: 'mL/(m·bar)', type: 'number', default: 0.012, min: 0, max: 1, group: 'Hose (if used)', role: 'material', help: 'From hose qualification data; it dominates the response time of hose lines.' },
    { key: 'eta', label: 'Utilisation factor at design pressure', unit: '-', type: 'number', default: 0.67, min: 0.3, max: 1, group: 'Acceptance criteria', role: 'criteria', help: `Allowable von Mises stress / yield in normal operation. ${CONFIRM}` },
    { key: 'etaTest', label: 'Utilisation factor at proof test', unit: '-', type: 'number', default: 0.96, min: 0.5, max: 1, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'proof', label: 'Proof test factor', unit: '-', type: 'number', default: 1.5, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: `Test pressure / design pressure. ${CONFIRM}` },
    { key: 'burstRatio', label: 'Hose burst / design ratio', unit: '-', type: 'number', default: 4, min: 2, max: 6, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'gamma', label: 'DNV resistance factor γm·γSC', unit: '-', type: 'number', default: 1.31, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: 'DNV-ST-F101: γm = 1.15 × γSC,pressure containment = 1.138 (medium safety class).' },
    { key: 'surgeAllow', label: 'Allowable surge pressure / design pressure', unit: '-', type: 'number', default: 1.1, min: 1, max: 1.5, group: 'Acceptance criteria', role: 'criteria', help: 'Project criterion for incidental pressure (10 % over design is the DNV incidental-pressure convention).' },
  ],
  compute(inp, ctx) {
    const m = ctx.material(inp.mat); need(m?.fy && m?.fu, `Material '${inp.mat}' has no yield / tensile strength in the library.`);
    const D = inp.od / 1e3, t = inp.wall / 1e3, tMin = t * (1 - inp.tol / 100);
    need(t < D / 2, `Wall ${inp.wall} mm is not smaller than the tube radius ${inp.od / 2} mm.`);
    const ro = D / 2, ri = ro - t, riMin = ro - tMin, d = 2 * ri;
    const pd = inp.pd * 1e5, pe = 1025 * G * inp.depth, head = inp.rhoF * G * inp.depth, pSub = pd + head;
    const E = m.E, fy = m.fy, fu = m.fu;
    // operating stress at the subsea end (highest internal, with external hydrostatic), minimum wall
    const op = lameStresses({ pInt: pSub, pExt: pe, ri: riMin, ro });
    const top = lameStresses({ pInt: pd, pExt: 0, ri: riMin, ro });
    const vm = Math.max(op.vonMises, top.vonMises), govEnd = top.vonMises >= op.vonMises ? 'topside' : 'subsea';
    const lim = thickWallLimitPressures({ ri: riMin, ro, yieldStress: fy, uts: fu });
    // ASME B31.3 §302.3.2 basis: S = min(UTS/3, 2/3·yield)
    const S = Math.min(fu / 3, 2 * fy / 3);
    const b31 = run(() => b313PressureDesign({ D, S, t: tMin, E: inp.weld }), 'ASME B31.3 pressure design');
    const tReq = run(() => b313PressureDesign({ D, S, P: pd, E: inp.weld }), 'ASME B31.3 required wall').t / (1 - inp.tol / 100);
    const dnv = dnvBurstPressure({ D, t: tMin, fy, fu });
    const pt = proofTestPressure({ designPressure: pd, proofFactor: inp.proof, ri: riMin, ro, yieldStress: fy });
    // transient: Korteweg wave speed c = sqrt(K/ρ / (1 + K·d/(E·t))), Joukowsky Δp = ρ·c·Δv for instantaneous closure
    const K = inp.bulk * 1e9, c = Math.sqrt(K / inp.rhoF / (1 + K * d / (E * t)));
    const Q = inp.flow / 60000, v = Q / (Math.PI * d * d / 4), dpSurge = inp.rhoF * c * v;
    const dp = run(() => tubePressureDrop({ flow: Q, d, length: inp.length, rho: inp.rhoF, viscosity: inp.visc / 1e3, roughness: inp.rough / 1e3 }), 'Pressure drop');
    const rt = run(() => hydraulicResponseTime({ length: inp.length, d, viscosity: inp.visc / 1e3, bulkModulus: K, fraction: 0.9 }), 'Response time');
    const hose = inp.hoseBurst > 0 ? hoseRating({ burstPressure: inp.hoseBurst * 1e5, burstRatio: inp.burstRatio, volumetricExpansion: inp.hoseVe * 1e-6 / 1e5, length: inp.length, pressureChange: inp.hosePd * 1e5 }) : null;
    const rtHose = hose ? hydraulicResponseTime({ length: inp.length, d, viscosity: inp.visc / 1e3, bulkModulus: K, volumetricExpansion: inp.hoseVe * 1e-6 / 1e5, fraction: 0.9 }) : null;
    // pressure-cycle fatigue: hoop stress range at the bore for the pressure range (Lamé), weld SCF in the curve class
    const dP = inp.cycleFrac * pd, dS = dP * (ro * ro + riMin * riMin) / (ro * ro - riMin * riMin);
    const N = dS > 0 ? cyclesToFailure({ curve: snCurve(inp.sn, 'air'), stressRange: dS }) : Infinity;
    const dmgYr = inp.cyclesYr > 0 && Number.isFinite(N) ? inp.cyclesYr / N : 0;
    const checks = [
      check(`von Mises stress at design pressure (${govEnd} end)`, vm / 1e6, inp.eta * fy / 1e6, 'MPa', `${STD_UMB} §7.2 (utilisation ${inp.eta}); Lamé thick-wall stresses`, 'Bore von Mises stress of the closed-end tube at design pressure with minimum wall. Above the allowable the tube yields at the bore in normal operation.'),
      check('ASME B31.3 internal design pressure', pd / 1e5, b31.P / 1e5, 'bar', 'ASME B31.3 (2022) §304.1.2 Eq. (3a), S per §302.3.2', 'Design pressure against the code allowable for the minimum wall, P = 2·S·E·t/(D − 2·Y·t). It is the pressure-containment requirement for process-piping practice.'),
      check('DNV pressure containment (burst)', Math.max(pSub - pe, pd) / 1e5, dnv.pb / inp.gamma / 1e5, 'bar', 'DNV-ST-F101 (2021) §5.4.2 Eq. 5.8', 'Local differential pressure against the factored burst resistance (2/√3)·(2t/(D−t))·min(fy, fu/1.15). Exceeding it risks rupture of the tube wall.'),
      check('Proof test stress', pt.vonMisesBore / 1e6, inp.etaTest * fy / 1e6, 'MPa', `${STD_UMB} §11 (proof test at ${inp.proof} × design pressure)`, 'Bore von Mises stress during the factory or post-lay pressure test. The test must prove the tube without yielding it.'),
      check('Surge pressure on rapid valve closure', (pd + dpSurge) / 1e5, inp.surgeAllow * pd / 1e5, 'bar', 'Joukowsky Δp = ρ·c·Δv; incidental-pressure limit as a project criterion', 'Design pressure plus the water-hammer rise from stopping the flow instantly against the allowable incidental pressure. Repeated surges fatigue the tube and can lift relief valves.'),
      check('Pressure-cycle fatigue over the design life', inp.dff * dmgYr * inp.life, 1, '-', 'DNV-RP-C203 (2021) §2.2; design fatigue factor per ISO 13628-5 / API 17E pattern', 'Miner damage from pressurise / depressurise cycles of the bore hoop stress with the design fatigue factor. Above 1.0 a seam-weld fatigue crack is predicted within the service life.'),
    ];
    if (hose) checks.push(check('Hose design pressure (burst ratio)', inp.hosePd, hose.designPressure / 1e5, 'bar', `${STD_UMB} §7 (hose burst-to-design ratio ${inp.burstRatio}:1)`, 'Hose working pressure against burst pressure divided by the required ratio. The large ratio covers the creep rupture and ageing of the reinforcement.'));
    const gov = worst(checks);
    const rr = linspace(riMin, ro, 15).map((r) => ({ r, ...lameStresses({ pInt: pSub, pExt: pe, ri: riMin, ro, r }) }));
    const walls = linspace(Math.max(0.3, inp.wall * 0.5), Math.min(inp.od / 2 * 0.9, inp.wall * 2), 25);
    const tt = linspace(0, 5 * rt.tau1, 60);
    const xs = linspace(0, inp.length, 30);
    const steelT = Math.PI * (ro * ro - ri * ri) * (m.rho ?? 7800) * inp.length / 1e3;
    const steelReq = Math.PI * (ro * ro - (ro - Math.min(tReq, ro * 0.95)) ** 2) * (m.rho ?? 7800) * inp.length / 1e3;
    return {
      kpis: [
        { key: 'vmBore', label: 'von Mises stress at the bore (design pressure)', value: round(vm / 1e6), unit: 'MPa', limit: round(inp.eta * fy / 1e6), sense: 'max', explain: 'Combined hoop, radial and end-cap axial stress at the inner wall — the highest-stressed point of a pressurised tube — using the minimum wall.' },
        { key: 'utilisation', label: 'Pressure-containment utilisation (governing)', value: round(gov.uc), unit: '-', limit: 1, sense: 'max', explain: 'Largest demand / capacity ratio over the stress, code, burst, test, surge and fatigue checks.' },
        { key: 'pB313', label: 'ASME B31.3 allowable pressure', value: round(b31.P / 1e5), unit: 'bar', explain: 'Internal design pressure permitted by B31.3 for the minimum wall and the code allowable stress min(UTS/3, 2/3·yield).' },
        { key: 'pBurstDnv', label: 'DNV burst resistance', value: round(dnv.pb / 1e5), unit: 'bar', explain: 'Characteristic pressure-containment resistance in the DNV-ST-F101 form, before resistance factors.' },
        { key: 'pYield', label: 'First-yield pressure (Lamé)', value: round(lim.pYield / 1e5), unit: 'bar', explain: 'Internal pressure at which the bore of the closed-end tube first yields.' },
        { key: 'pPlastic', label: 'Ductile burst estimate', value: round((lim.pBurstEstimate ?? lim.pPlastic) / 1e5), unit: 'bar', explain: 'Fully plastic pressure (2/√3)·σu·ln(ro/ri): an estimate of the actual burst pressure for comparison with qualification burst tests.' },
        { key: 'burstMargin', label: 'Burst estimate / design pressure', value: round((lim.pBurstEstimate ?? lim.pPlastic) / pd), unit: '-', explain: 'Ratio of estimated burst to design pressure. Steel umbilical tubes typically show ratios of 2.5–4.' },
        { key: 'pTest', label: 'Proof test pressure', value: round(pt.testPressure / 1e5), unit: 'bar', explain: `Design pressure × ${inp.proof}. The test stress utilisation is ${fmt(pt.utilisation)} of yield.` },
        { key: 'tReq', label: 'Required nominal wall (B31.3)', value: round(tReq * 1e3), unit: 'mm', explain: 'Nominal wall that satisfies the B31.3 pressure design at the design pressure, including the under-tolerance.' },
        { key: 'waveSpeed', label: 'Pressure-wave speed', value: round(c), unit: 'm/s', explain: 'Speed of a pressure pulse in the fluid-filled elastic tube (Korteweg). It sets the travel time of a pressure change and the size of the surge.' },
        { key: 'surge', label: 'Joukowsky surge pressure rise', value: round(dpSurge / 1e5), unit: 'bar', explain: 'Pressure rise ρ·c·v if the flow is stopped faster than the wave return time 2L/c — the upper bound of the valve-closure transient.' },
        { key: 'dpFlow', label: 'Friction pressure drop along the tube', value: round(dp.dp / 1e5), unit: 'bar', explain: `Darcy–Weisbach loss at the stated flow (${dp.regime}, Re = ${fmt(dp.reynolds)}). It must leave enough pressure at the subsea end for the function served.` },
        { key: 'response', label: 'Hydraulic response time (90 %)', value: round(rt.time), unit: 's', explain: 'Time for the far-end pressure to reach 90 % of a step applied at the topside end (laminar diffusion of the pressure signal). It governs valve actuation and emergency shutdown timing.' },
        { key: 'cycleLife', label: 'Pressure-cycle fatigue life (no DFF)', value: round(dmgYr > 0 ? Math.min(1 / dmgYr, 1e6) : 1e6), unit: 'years', limit: inp.dff * inp.life, sense: 'min', explain: 'Years of pressure cycling to a Miner sum of 1.0; it must exceed the service life times the design fatigue factor.' },
        ...(hose ? [{ key: 'hoseDesign', label: 'Hose allowable design pressure', value: round(hose.designPressure / 1e5), unit: 'bar', limit: inp.hosePd, sense: 'min', explain: 'Hose burst pressure divided by the required burst ratio.' },
          { key: 'hoseResponse', label: 'Response time with hose expansion (90 %)', value: round(rtHose.time), unit: 's', explain: 'Response time of the same bore when the wall expands like a thermoplastic hose: the reason long control lines use steel tube.' }] : []),
      ],
      checks,
      plots: [
        { type: 'line', title: 'Stresses through the wall (subsea end, design pressure)', x: { label: 'Radius', unit: 'mm' }, y: { label: 'Stress', unit: 'MPa' },
          series: [{ name: 'Hoop', x: rr.map((q) => q.r * 1e3), y: rr.map((q) => q.hoop / 1e6) }, { name: 'Radial', x: rr.map((q) => q.r * 1e3), y: rr.map((q) => q.radial / 1e6) }, { name: 'Axial', x: rr.map((q) => q.r * 1e3), y: rr.map((q) => q.axial / 1e6) }, { name: 'von Mises', x: rr.map((q) => q.r * 1e3), y: rr.map((q) => q.vonMises / 1e6) }, { name: 'Allowable', x: [riMin * 1e3, ro * 1e3], y: [inp.eta * fy / 1e6, inp.eta * fy / 1e6], dash: true }],
          explain: 'Lamé distribution from bore to outer surface. In a thick small-bore tube the hoop stress falls markedly through the wall; the von Mises curve must stay below the dashed allowable at the bore.' },
        { type: 'line', title: 'Allowable pressure vs wall thickness', x: { label: 'Nominal wall', unit: 'mm' }, y: { label: 'Pressure', unit: 'bar' },
          series: [{ name: 'ASME B31.3', x: walls, y: walls.map((w) => b313PressureDesign({ D, S, t: w / 1e3 * (1 - inp.tol / 100), E: inp.weld }).P / 1e5) }, { name: 'DNV burst / γ', x: walls, y: walls.map((w) => dnvBurstPressure({ D, t: w / 1e3 * (1 - inp.tol / 100), fy, fu }).pb / inp.gamma / 1e5) }, { name: `Lamé at η = ${inp.eta}`, x: walls, y: walls.map((w) => inp.eta * thickWallLimitPressures({ ri: ro - w / 1e3 * (1 - inp.tol / 100), ro, yieldStress: fy }).pYield / 1e5) }, { name: 'Design pressure', x: [walls[0], walls[walls.length - 1]], y: [inp.pd, inp.pd], dash: true }, { name: 'Selected wall', x: [inp.wall, inp.wall], y: [0, inp.pd], dash: true }],
          explain: 'Pressure permitted by each criterion as the wall grows. The wall at which the lowest curve crosses the design-pressure line is the minimum acceptable wall; the selected wall should lie to its right.' },
        { type: 'line', title: 'Far-end pressure response to a step', x: { label: 'Time', unit: 's' }, y: { label: 'Far-end pressure / step', unit: '-' },
          series: [{ name: 'Steel tube', x: tt, y: tt.map((q) => Math.max(0, 1 - 4 / Math.PI * Math.exp(-q / rt.tau1))) }, ...(rtHose ? [{ name: 'Hose (with volumetric expansion)', x: tt, y: tt.map((q) => Math.max(0, 1 - 4 / Math.PI * Math.exp(-q / rtHose.tau1))), dash: true }] : [])],
          explain: 'First-mode solution of the laminar line-diffusion equation for a closed far end. The curve reaching 0.9 is the response time; a compliant hose stretches it many times over.' },
        { type: 'line', title: 'Pressure along the tube (flowing)', x: { label: 'Distance from topside', unit: 'm' }, y: { label: 'Internal pressure', unit: 'bar' },
          series: [{ name: 'Internal (with static head, linear depth)', x: xs, y: xs.map((x) => (pd + head * x / inp.length - dp.dp * x / inp.length) / 1e5) }, { name: 'External hydrostatic', x: xs, y: xs.map((x) => pe * x / inp.length / 1e5), dash: true }],
          explain: 'Internal pressure from the topside design value, increased by the static head of the fluid column and reduced by friction, against the external sea pressure. The gap is the differential pressure the wall must contain.' },
      ],
      tables: [{ title: 'Pressure-containment summary', explain: 'Each row is one criterion with the pressure or stress demand, its capacity and the utilisation. All use the minimum wall after the under-tolerance.',
        columns: [{ key: 'name', label: 'Criterion' }, { key: 'demand', label: 'Demand' }, { key: 'capacity', label: 'Capacity' }, { key: 'unit', label: 'Unit' }, { key: 'uc', label: 'Utilisation', unit: '-' }],
        rows: checks.map((q) => ({ name: q.name, demand: q.demand, capacity: q.capacity, unit: q.unit, uc: round(q.uc) })) }],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}). Minimum wall ${fmt(tMin * 1e3)} mm, D/t = ${fmt(D / tMin)}.`,
        `Pressure at the subsea end ${fmt(pSub / 1e5)} bar internal against ${fmt(pe / 1e5)} bar external; test pressure ${fmt(pt.testPressure / 1e5)} bar.`,
        `The ISO 13628-5 tube burst formula is not implemented here; B31.3, DNV and Lamé checks are used with overridable factors. ${CONFIRM}`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}. B31.3 allows ${fmt(b31.P / 1e5)} bar, DNV factored burst ${fmt(dnv.pb / inp.gamma / 1e5)} bar, design ${inp.pd} bar; surge adds ${fmt(dpSurge / 1e5)} bar.`,
          gov.uc > 1 ? (gov.name.includes('Surge') ? `Slow the valve closure to longer than the wave return time 2L/c = ${fmt(2 * inp.length / c)} s, or add an accumulator at the subsea end; alternatively reduce the flow velocity (${fmt(v)} m/s).` : `Increase the nominal wall to at least ${fmt(Math.max(tReq, t * gov.uc) * 1e3)} mm or select a higher-strength grade, then re-check bending and fatigue with the heavier tube.`)
            : 'Confirm the factors against the governing edition and the tube manufacturer’s qualified burst and proof tests; carry the tube stress at pressure into the tension-and-bending capacity curve.',
          `ASME B31.3 (2022) §304.1.2; DNV-ST-F101 (2021) §5.4.2; ${STD_UMB} §7 and §11`),
        adv('economic', 'info', `The tube contains ${fmt(steelT)} t of ${m.name} over ${fmt(inp.length / 1e3)} km; the code-minimum wall of ${fmt(tReq * 1e3)} mm would need ${fmt(steelReq)} t.`,
          tReq < t ? `The selected wall exceeds the B31.3 minimum by ${fmt((t / tReq - 1) * 100)} %: up to ${fmt(steelT - steelReq)} t of alloy per tube could be saved if bending, fatigue and collapse do not need it — or keep it as corrosion and fatigue margin.` : `The wall must grow by ${fmt((tReq / t - 1) * 100)} % (${fmt(steelReq - steelT)} t more alloy per tube); a higher-strength grade may reach the pressure rating with less mass and a smaller umbilical diameter.`,
          'Quantity-based comparison (AACE RP 18R-97 Class 3); alloy tube priced per tonne by the supplier',
          `${fmt(Math.abs(steelT - steelReq))} t of tube ≈ ${fmt(Math.abs(steelT - steelReq) * CO2_STEEL)} t CO₂e at the screening factor of ${CO2_STEEL} t CO₂e/t for steel; super-duplex stainless is typically 3–4 times that per tonne, so wall optimisation matters more than for carbon steel.`),
      ],
      publish: { 'umbilical.tubeUtilisation': gov.uc, 'umbilical.tubeBurstPressure': dnv.pb, 'umbilical.hydraulicResponseTime': rt.time },
      quantities: { steel_t: steelT },
    };
  },
  verify: [
    // ro = 8.25, ri = 6.35 mm, p = 69 MPa, no external pressure or head: closed-end Lamé + von Mises at the bore
    // σvm = √3·p·ro²/(ro² − ri²) = 1.7321 × 69 × 68.0625/27.74 = 293.2 MPa.
    // B31.3: S = min(800/3, 2·550/3) = 266.67 MPa; P = 2·S·t/(D − 2·0.4·t) = 2 × 266.67 × 1.9/(16.5 − 1.52) = 67.65 MPa = 676.5 bar.
    // DNV: pb = (2t/(D − t))·min(fy, fu/1.15)·2/√3 = (3.8/14.6) × 550 × 1.1547 = 165.30 MPa = 1653 bar.
    { name: 'Thick-wall tube: Lamé bore stress, ASME B31.3 Eq. (3a) and DNV-ST-F101 Eq. 5.8',
      source: 'Timoshenko & Goodier §28; ASME B31.3 §304.1.2; DNV-ST-F101 (2021) §5.4.2 — hand calculation',
      inputs: { pd: 690, depth: 0, tol: 0, flow: 0, cyclesYr: 0, hoseBurst: 0 },
      expect: [{ kpi: 'vmBore', value: 293.2, tol: 0.003 }, { kpi: 'pB313', value: 676.5, tol: 0.002 }, { kpi: 'pBurstDnv', value: 1653, tol: 0.002 }] },
    // Korteweg: c = √(K/ρ/(1 + K·d/(E·t))) = √(2.2e9/1060/1.07353) = 1390.4 m/s; v = Q/A = 1.6667e-4/1.26677e-4 = 1.3157 m/s;
    // Joukowsky Δp = ρ·c·v = 1060 × 1390.4 × 1.3157 = 19.39 bar.
    { name: 'Water hammer: Korteweg wave speed and Joukowsky pressure rise',
      source: 'Joukowsky (1898); Wylie & Streeter, Fluid Transients (1978) ch. 1 — hand calculation',
      inputs: { flow: 10 },
      expect: [{ kpi: 'waveSpeed', value: 1390.4, tol: 0.002 }, { kpi: 'surge', value: 19.39, tol: 0.003 }] },
  ],
};
