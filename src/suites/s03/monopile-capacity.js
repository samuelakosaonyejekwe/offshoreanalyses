// Monopile lateral and axial capacity: ultimate lateral resistance by load stepping on the p–y model, Broms
// rigid-pile cross-check, API axial shaft and end-bearing capacity, embedded-length criteria and steel quantity.
import { ultimateLateralCapacity, bromsSand, bromsClay } from '../../engine/geotech/pile.js';
import { embeddedLengthCriteria, monopileSizingScreen } from '../../engine/geotech/monopile.js';
import { axialCapacity, capacityVsPenetration } from '../../engine/geotech/axialcapacity.js';
import { need, run, fmt, round, adv, check, worst, sev, steelCo2Text, supportInputs, soilInputs, steelInputs, steelOf, soilProfile, soilLabel, LIMIT_PY, LIMIT_CODE, LIMIT_HIFI, CONFIRM } from './_common.js';
import { pyInputs, pileOptions } from './_pile.js';

const geo = supportInputs().filter((i) => ['depth', 'D', 't', 'Lp', 'zInt'].includes(i.key));
const py = pyInputs().filter((i) => i.key !== 'kLin').map((i) => (i.key === 'pyMode' ? { ...i, options: i.options.filter((o) => o.value !== 'linear') } : i));

export default {
  id: 's03-monopile-capacity', suite: 3, title: 'Monopile lateral and axial capacity',
  summary: 'Ultimate lateral capacity under combined shear and moment, axial compression and tension capacity, embedded-length criteria and the pile steel quantity, with load and resistance factors.',
  covers: [6], fidelity: 'numerical',
  standards: ['DNV-ST-0126 (2021) §7.6 (pile foundations)', 'API RP 2GEO (2011, add. 2014) §8.1–8.5', 'API RP 2A-WSD (21st ed.) §6.4–6.5, §6.10.6', 'ISO 19901-4:2016 §8', 'IEC 61400-3-1:2019 §7.6 (partial safety factors)', 'Broms (1964a, b)'],
  models: ['Nonlinear lateral soil resistance (p–y) with load stepping to the ultimate state', 'Broms rigid-pile ultimate lateral resistance', 'API axial pile resistance: α-method (clay), β-method with limits (sand), plugged / unplugged end bearing', 'Embedded-length criteria: zero toe-kick, vertical tangent, critical length', 'Force and moment equilibrium'],
  limits: `${LIMIT_PY} Lateral failure is defined by a head-displacement criterion (default 10 % of the diameter) or by loss of equilibrium, whichever comes first. Group effects, set-up and ageing are not included. ${LIMIT_CODE} ${LIMIT_HIFI}`,
  publishes: ['quantities.steel_t'],
  searchTerms: 'monopile lateral capacity axial capacity embedded length criteria API',
  inputs: [
    ...geo,
    ...steelInputs(),
    ...soilInputs(),
    ...py,
    { key: 'H', label: 'Horizontal load at the mudline (characteristic)', unit: 'MN', type: 'number', default: 4.9, min: 0.01, max: 100, group: 'Loads at the mudline', role: 'bc', bind: 'owt.mudlineShear', uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'M', label: 'Overturning moment at the mudline (characteristic)', unit: 'MN·m', type: 'number', default: 280, min: 0, max: 5000, group: 'Loads at the mudline', role: 'bc', bind: 'owt.mudlineMoment', uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'V', label: 'Vertical compression at the mudline', unit: 'MN', type: 'number', default: 12, min: 0, max: 300, group: 'Loads at the mudline', role: 'bc' },
    { key: 'gammaF', label: 'Load factor', unit: '-', type: 'number', default: 1.35, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: `IEC 61400-3-1 normal design situations. ${CONFIRM}` },
    { key: 'gammaRlat', label: 'Resistance factor, lateral', unit: '-', type: 'number', default: 1.15, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: `Applied to the ultimate lateral capacity. ${CONFIRM}` },
    { key: 'gammaRax', label: 'Resistance factor, axial', unit: '-', type: 'number', default: 1.25, min: 1, max: 2.5, group: 'Acceptance criteria', role: 'criteria', help: `Applied to the axial capacity (DNV-ST-0126 quotes 1.25 for piles in the ULS). ${CONFIRM}` },
    { key: 'yLimD', label: 'Lateral failure criterion: head displacement / D', unit: '-', type: 'number', default: 0.1, min: 0.02, max: 0.3, group: 'Acceptance criteria', role: 'criteria', help: 'Displacement at which the pile is considered to have failed laterally if equilibrium is still found.' },
    { key: 'nElem', label: 'Pile elements', type: 'int', default: 40, min: 10, max: 160, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Pile elements', levels: [20, 40, 80], qoi: 'Hult' },
  compute(inp, ctx) {
    const steel = steelOf(inp, ctx);
    const H = inp.H * 1e6, M = inp.M * 1e6, V = inp.V * 1e6;
    const profile = soilProfile(inp);
    const opt = pileOptions(inp, profile, steel);
    const sec = opt.sec, e = M / H, tW = inp.t / 1e3;
    const ult = run(() => ultimateLateralCapacity(opt, { eccentricity: e, yLimit: inp.yLimD * inp.D, nSteps: 30, tol: 1e-3 }), 'Ultimate lateral capacity');
    const Hd = inp.gammaF * H, HRd = ult.Hult / inp.gammaRlat;
    // Broms closed form for uniform profiles only (rigid short pile; long-pile branch needs the yield moment).
    const My = steel.fy * sec.W;
    let broms = null;
    if (inp.soil === 'sand') broms = bromsSand({ gammaSub: inp.gammaSub * 1e3, phiDeg: inp.phiDeg, D: inp.D, L: inp.Lp, e, My });
    else if (inp.soil === 'clay') broms = bromsClay({ su: (inp.su0 + 0.5 * inp.suGrad * inp.Lp) * 1e3, D: inp.D, L: inp.Lp, e, My });
    const ax = run(() => axialCapacity({ profile, L: inp.Lp, D: inp.D, t: tW }), 'Axial capacity');
    const Wsub = sec.A * (steel.rho - 1025) * 9.80665 * inp.Lp; // submerged weight of the embedded pile helps in tension
    const QcRd = ax.Qc / inp.gammaRax, QtRd = (ax.Qt + Wsub) / inp.gammaRax;
    const crit = run(() => embeddedLengthCriteria({ ...opt, L: undefined, head: { H: Hd, M: inp.gammaF * M } }, { Lmin: 2.5 * inp.D, Lmax: Math.min(Math.max(1.6 * inp.Lp, 9 * inp.D), profile.depth), nLengths: 14 }), 'Embedded-length criteria');
    const cvp = run(() => capacityVsPenetration({ profile, D: inp.D, t: tW, Lmin: 4, Lmax: Math.min(1.5 * inp.Lp, profile.depth), step: Math.max(1, Math.round(inp.Lp / 25)) }), 'Axial capacity against penetration');
    const size = monopileSizingScreen({ D: inp.D, t: tW });
    const checks = [
      check('Lateral capacity (shear with moment)', Hd / 1e6, HRd / 1e6, 'MN', 'DNV-ST-0126 (2021) §7.6.2; API RP 2GEO (2011) §8.5', 'Guards against the soil failing around the pile so that it rotates without limit. The load is applied at the eccentricity M/H; exceedance means collapse of the foundation in the design storm.'),
      check('Axial compression capacity', inp.gammaF * V / 1e6, QcRd / 1e6, 'MN', 'DNV-ST-0126 (2021) §7.6.3; API RP 2GEO (2011) §8.1', 'Shaft friction plus end bearing against the weight of turbine, tower and substructure. Rarely governs a monopile but must be shown.'),
      check('Pile first-yield moment against factored mudline moment', inp.gammaF * M / 1e6, My / 1.1 / 1e6, 'MN·m', 'DNV-ST-0126 (2021) §4.4 (elastic section, γM = 1.1)', 'Guards against the steel yielding before the soil resistance is mobilised: a plastic hinge in the pile would then govern the lateral capacity instead of the soil.'),
      check('Minimum wall thickness for driving', size.tMin * 1e3, inp.t, 'mm', 'API RP 2A-WSD §6.10.6: t ≥ 6.35 mm + D/100', 'Guards against local buckling and tip damage during driving. Thinner walls need a drivability and tip-integrity assessment.'),
    ];
    const critL = crit.critical ?? crit.verticalTangent ?? crit.zeroToeKick;
    const lenOk = critL == null ? null : inp.Lp >= critL;
    const gov = worst(checks), ok = gov.uc <= 1;
    const mPile = sec.A * steel.rho * (inp.Lp + inp.depth + inp.zInt) / 1e3;
    const perM = sec.A * steel.rho / 1e3;
    const LpOpt = critL != null ? Math.max(critL, 3 * inp.D) : inp.Lp;
    const dL = LpOpt - inp.Lp;
    const fin = (a) => a.map((v) => (Number.isFinite(v) ? v : null));
    const Lok = crit.L.filter((_, i) => Number.isFinite(crit.headY[i]));
    return {
      kpis: [
        { key: 'Hult', label: 'Ultimate lateral capacity', unit: 'MN', value: ult.Hult / 1e6, explain: `Horizontal load at the mudline, applied with the moment arm M/H = ${fmt(e)} m, at which the pile fails laterally (${ult.criterion}). Divide by the resistance factor before comparing with the factored load.` },
        { key: 'Mult', label: 'Moment at lateral failure', unit: 'MN·m', value: ult.Hult * e / 1e6, explain: 'Mudline moment accompanying the ultimate horizontal load at the same M/H ratio.' },
        { key: 'ucLat', label: 'Lateral utilisation', unit: '-', value: Hd / HRd, limit: 1, sense: 'max', explain: 'Factored horizontal load over factored lateral capacity. Above 1 the pile is too short or too slender for the design storm.' },
        { key: 'Qc', label: 'Axial compression capacity', unit: 'MN', value: ax.Qc / 1e6, explain: `Outer shaft friction plus the smaller of plugged and unplugged end bearing (${ax.plugged ? 'plugged' : 'unplugged / coring'} governs).` },
        { key: 'Qshaft', label: 'Outer shaft friction', unit: 'MN', value: ax.Qfo / 1e6, explain: 'Integrated unit friction on the outer pile wall (α-method in clay, β-method with limiting values in sand).' },
        { key: 'Qt', label: 'Axial tension capacity', unit: 'MN', value: ax.Qt / 1e6, explain: 'Shaft friction available in uplift, without pile weight. Relevant for jacket piles; a monopile is not loaded in tension.' },
        { key: 'ucAx', label: 'Axial utilisation', unit: '-', value: inp.gammaF * V / QcRd, limit: 1, sense: 'max', explain: 'Factored vertical load over factored compression capacity.' },
        { key: 'Lcrit', label: 'Critical embedded length', unit: 'm', value: crit.critical ?? 0, explain: 'Shortest length beyond which extra penetration no longer reduces mudline displacement and rotation (within 2 %). Zero when not reached in the scanned range.' },
        { key: 'LvertTan', label: 'Vertical-tangent length', unit: 'm', value: crit.verticalTangent ?? 0, explain: 'Shortest length at which the deflected pile is vertical at the toe; a traditional (conservative) monopile length criterion. Zero when not reached.' },
        { key: 'LzeroToe', label: 'Zero toe-kick length', unit: 'm', value: crit.zeroToeKick ?? 0, explain: 'Shortest length at which the toe no longer moves backwards. Zero when not reached in the scanned range.' },
        { key: 'dOverT', label: 'Diameter / thickness', unit: '-', value: size.dOverT, explain: 'Slenderness of the pile wall. Above about 60 a local buckling check is required; above about 120 driving damage and fabrication tolerances become critical.' },
        { key: 'tMinApi', label: 'API minimum wall for driving', unit: 'mm', value: size.tMin * 1e3, explain: 'Minimum wall thickness for a driven pile, 6.35 mm + D/100.' },
        { key: 'steelT', label: 'Monopile steel mass (toe to interface)', unit: 't', value: mPile, explain: 'Mass of the primary pile steel at constant wall thickness, from the toe to the tower interface. Add about 8–12 % for secondary steel and flanges.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Lateral load–displacement to failure', x: { label: 'Mudline displacement', unit: 'mm' }, y: { label: 'Horizontal load', unit: 'MN' },
          series: [{ name: 'Response', x: ult.curve.y.map((v) => v * 1e3), y: ult.curve.H.map((v) => v / 1e6) }, { name: 'Factored load', x: [0, ult.yAtHult * 1e3], y: [Hd / 1e6, Hd / 1e6], dash: true }, { name: 'Design capacity', x: [0, ult.yAtHult * 1e3], y: [HRd / 1e6, HRd / 1e6], dash: true }],
          explain: 'Load stepping to the ultimate state. The factored load line must lie below the design capacity line; a flat end of the curve is soil failure, a cut-off is the displacement criterion.' },
        { type: 'line', title: 'Mudline response against embedded length', x: { label: 'Embedded length', unit: 'm' }, y: { label: 'Mudline displacement', unit: 'mm' },
          series: [{ name: 'Mudline displacement at factored load', x: Lok, y: fin(crit.headY).filter((v) => v != null).map((v) => v * 1e3) }, { name: 'Present length', x: [inp.Lp, inp.Lp], y: [0, Math.max(...fin(crit.headY).filter((v) => v != null)) * 1e3], dash: true }],
          explain: 'Displacement falls steeply with length until the pile behaves as flexible, then levels off: the knee is the economic embedded length.' },
        { type: 'line', title: 'Axial capacity against penetration', x: { label: 'Capacity', unit: 'MN' }, y: { label: 'Penetration', unit: 'm' },
          series: [{ name: 'Compression', x: cvp.Qc.map((v) => v / 1e6), y: cvp.L.map((v) => -v) }, { name: 'Tension', x: cvp.Qt.map((v) => v / 1e6), y: cvp.L.map((v) => -v) }, { name: 'Factored load × γR', x: [inp.gammaF * V * inp.gammaRax / 1e6, inp.gammaF * V * inp.gammaRax / 1e6], y: [-cvp.L[0], -cvp.L[cvp.L.length - 1]], dash: true }],
          explain: 'Axial capacity as the pile is lengthened. The crossing with the dashed line is the penetration needed for the vertical load alone.' },
      ],
      tables: [
        { title: 'Capacity summary', columns: [{ key: 'k', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['Ultimate lateral capacity (p–y, load stepping)', ult.Hult / 1e6, 'MN'], ['Head displacement at failure', ult.yAtHult * 1e3, 'mm'], ['Maximum pile moment at failure', ult.maxMoment / 1e6, 'MN·m'], ['Pile first-yield moment fy·W', My / 1e6, 'MN·m'],
            ...(broms ? [[`Broms ultimate lateral load (${broms.mode} pile)`, broms.Hu / 1e6, 'MN']] : []),
            ['Outer shaft friction', ax.Qfo / 1e6, 'MN'], ['Inner shaft friction', ax.Qfi / 1e6, 'MN'], ['End bearing used', ax.Qp / 1e6, 'MN'], ['Unit end bearing', ax.qTip / 1e6, 'MPa'], ['Design tension resistance incl. pile weight', QtRd / 1e6, 'MN']].map(([k, v, u]) => ({ k, v: round(v), u })),
          explain: 'Components of the lateral and axial capacity. The Broms value is an independent rigid-pile estimate for uniform soil: agreement within a factor of about 1.5 is expected because the failure criteria differ.' },
      ],
      notes: [`Soil: ${inp.pyMode === 'user' ? 'user p–y curves for lateral response; ' : ''}${soilLabel(inp)}.`, `${gov.name} governs (utilisation ${fmt(gov.uc)}).`,
        critL == null ? 'No embedded-length criterion is met within the scanned range: the pile behaves as rigid; consider a longer pile.' : `Embedded length ${inp.Lp} m is ${lenOk ? 'at or above' : 'below'} the governing length criterion of ${fmt(critL)} m.`,
        `Lateral failure criterion reached: ${ult.criterion}.`],
      advice: [
        ok ? adv('engineering', lenOk === false ? 'warn' : sev(gov.uc), `Capacity checks pass: lateral ${fmt(Hd / HRd * 100)} %, axial ${fmt(inp.gammaF * V / QcRd * 100)} %. ${lenOk === false ? `The length is below the stiffness-based criterion (${fmt(critL)} m), so displacement, not capacity, will govern.` : 'Capacity is not the governing design driver for this monopile — stiffness, fatigue and natural frequency normally are.'}`,
          lenOk === false ? `Check rotation and natural frequency with the pile–soil interaction module before accepting ${inp.Lp} m; ${fmt(critL)} m would make the pile fully effective.` : `Confirm the API-type soil-reaction curves with site-specific (PISA-type) data and the resistance factors against the governing edition of DNV-ST-0126; with the present reserve the length could be reviewed down towards ${fmt(Math.max(LpOpt, 3 * inp.D))} m if rotation and frequency allow.`,
          'DNV-ST-0126 (2021) §7.6; API RP 2GEO (2011) §8',
          `Changing the embedded length from ${inp.Lp} m to ${fmt(LpOpt)} m changes the steel by ${fmt(dL * perM)} t — ${steelCo2Text(dL * perM)}`)
          : adv('engineering', 'fail', `${gov.name} fails: ${fmt(gov.demand)} against ${fmt(gov.capacity)} ${gov.unit} (utilisation ${fmt(gov.uc)}).`,
            gov.name.includes('wall') ? `Increase the wall to at least ${fmt(size.tMin * 1e3)} mm or demonstrate drivability and tip integrity by analysis.` : gov.name.includes('Lateral') ? `Lengthen the pile (lateral capacity grows roughly with L² to L³ in sand) or increase the diameter; target a design capacity of ${fmt(Hd / 1e6)} MN.` : 'Increase penetration or diameter until the factored capacity exceeds the factored load.',
            'DNV-ST-0126 (2021) §7.6; API RP 2GEO (2011) §8; API RP 2A-WSD §6.10.6',
            `Each extra metre of pile is ${fmt(perM)} t of steel ≈ ${fmt(perM * 1.9)} t CO₂e (screening factor 1.9 t CO₂e/t) and deeper driving raises hammer energy and underwater noise.`),
        adv('economic', 'info', `The monopile uses ${fmt(mPile)} t of primary steel (${fmt(perM)} t/m). Capacity utilisation is ${fmt(Math.max(Hd / HRd, inp.gammaF * V / QcRd) * 100)} %.`,
          `Treat embedded length as the main cost lever: every metre removed saves ${fmt(perM)} t of plate, rolling and welding and shortens driving; verify against rotation, frequency and fatigue before committing.`,
          'Quantity-based estimating per AACE RP 18R-97 (Class 4); foundation supply priced per tonne of fabricated primary steel'),
      ],
      publish: { 'quantities.steel_t': mPile },
      quantities: { steel_t: round(mPile) },
    };
  },
  verify: [
    { name: 'Outer shaft friction in overconsolidated clay by the API α-method', source: 'API RP 2A-WSD §6.4.2 hand integration: su = 100 kPa, γ′ = 8 kN/m³, ψ = su/σ′v; z < 12.5 m: f = 0.5·su^0.75·(γ′z)^0.25 → ∫ = 500.0 kN/m; 12.5–40 m: f = 0.5·√(su·γ′·z) → ∫ = 1968.5 kN/m; total 2468.5 kN/m × π × 2 m = 15.51 MN',
      inputs: { soil: 'clay', su0: 100, suGrad: 0, gammaSub: 8, D: 2, t: 40, Lp: 40, H: 0.2, M: 2, V: 1 }, expect: [{ kpi: 'Qshaft', value: 15.51, tol: 0.01 }] },
    { name: 'API minimum wall thickness for driving: t = 6.35 + D/100 mm', source: 'API RP 2A-WSD (21st ed.) Eq. 6.10.6-1: D = 6000 mm → 66.35 mm', inputs: { D: 6 }, expect: [{ kpi: 'tMinApi', value: 66.35, tol: 1e-6 }] },
  ],
};
