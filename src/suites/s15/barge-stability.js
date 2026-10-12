// Transport / installation vessel intact stability and ballast: loading condition, equilibrium, GZ curve,
// IMO and wind-overturning criteria, free-surface correction and the ballast transfer that levels the barge.
import { boxMesh, boxHydrostatics, solveEquilibrium, freeSurfaceCorrection } from '../../engine/hydro/hydrostatics.js';
import { gzCurve, gzProperties, areaUnderGz, imoIntactCriteria, windHeelingLever, moduWindCriterion, ballastPlan } from '../../engine/hydro/stability.js';
import { G, DEG, need, run, round, fmt, kpi, check, fsCheck, adv, line, tbl, worst, vesselCo2Text, CONFIRM, N001 } from './_common.js';

const DEFAULT_TANKS = [{ name: 'Wing tanks aft (4)', n: 4, l: 28, b: 8.5 }, { name: 'Wing tanks forward (4)', n: 4, l: 28, b: 8.5 }];

export default {
  id: 's15-barge-stability', suite: 15, title: 'Vessel stability and ballast',
  summary: 'Builds the loading condition of a cargo barge, finds draft, trim and heel, the righting-lever curve with free-surface correction, and checks intact stability, wind overturning and the ballast needed to bring the barge level.',
  covers: [6], fidelity: 'numerical',
  standards: [`${N001} Sec. 11 (intact stability of transport vessels: GM, range, wind area ratio)`, 'IMO 2008 IS Code (MSC.267(85)) Part A §2.2 and Part B ch. 3 (free surface)', 'IMO MODU Code (2009) §3.2–3.3 (wind heeling, area ratio)'],
  models: ['Archimedes buoyancy and hydrostatic equilibrium (mesh re-clipping)', 'Metacentric stability', 'Righting-arm (GZ) curve at constant displacement', 'Free-surface correction', 'Wind heeling lever from projected areas', 'Bounded least-squares ballast plan'],
  limits: 'Box-shaped hull (or the form fields as a box) with calm-water hydrostatics: no wave, dynamic or damaged-stability assessment, and the cargo is assumed not to immerse. The stability limits (minimum GM, range, wind speed and area ratio) are recalled DNV-ST-N001 and IMO values exposed as inputs to be confirmed; damage stability and longitudinal strength are outside the module.',
  publishes: ['install.bargeGM', 'install.bargeDraft', 'install.bargeDisp'],
  searchTerms: 'cargo barge intact stability GZ curve wind overturning ballast plan transport',
  inputs: [
    { key: 'hull', label: 'Hull geometry (optional import)', type: 'geometry', kinds: ['surface'], default: null, group: 'Barge', role: 'geometry', help: 'A closed hull surface with the keel at z = 0 and midship at x = 0. When empty, a box barge is built from the dimensions below.' },
    { key: 'Lb', label: 'Barge length', unit: 'm', type: 'number', default: 120, min: 20, max: 400, group: 'Barge', role: 'geometry' },
    { key: 'Bb', label: 'Barge beam', unit: 'm', type: 'number', default: 36, min: 6, max: 120, group: 'Barge', role: 'geometry', design: true },
    { key: 'Db', label: 'Barge depth', unit: 'm', type: 'number', default: 7.6, min: 2, max: 30, group: 'Barge', role: 'geometry' },
    { key: 'mLight', label: 'Lightship mass', unit: 't', type: 'number', default: 5000, min: 50, max: 2e5, group: 'Loading condition', role: 'load' },
    { key: 'kgLight', label: 'Lightship KG', unit: 'm', type: 'number', default: 4.2, min: 0.2, max: 30, group: 'Loading condition', role: 'geometry' },
    { key: 'mCargo', label: 'Cargo mass incl. seafastening and grillage', unit: 't', type: 'number', default: 6450, min: 0, max: 2e5, group: 'Loading condition', role: 'load', bind: 'install.liftWeight', uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'zCargo', label: 'Cargo CoG height above the deck', unit: 'm', type: 'number', default: 14, min: 0, max: 150, group: 'Loading condition', role: 'geometry', uncertain: { dist: 'normal', cov: 0.05 } },
    { key: 'xCargo', label: 'Cargo CoG from midship (forward positive)', unit: 'm', type: 'number', default: 4, min: -200, max: 200, group: 'Loading condition', role: 'ic', help: 'Initial condition of the stowage; drives trim.' },
    { key: 'yCargo', label: 'Cargo CoG from the centreline', unit: 'm', type: 'number', default: 0.3, min: -60, max: 60, group: 'Loading condition', role: 'ic', help: 'Drives heel.' },
    { key: 'mBallast', label: 'Ballast on board', unit: 't', type: 'number', default: 8500, min: 0, max: 5e5, group: 'Loading condition', role: 'ic', design: true },
    { key: 'kgBallast', label: 'Ballast KG', unit: 'm', type: 'number', default: 2.2, min: 0, max: 30, group: 'Loading condition', role: 'geometry' },
    { key: 'tanks', label: 'Slack tanks (free surface)', type: 'table', group: 'Loading condition', role: 'ic',
      columns: [{ key: 'name', label: 'Tank group', type: 'text' }, { key: 'n', label: 'Number slack', type: 'number' }, { key: 'l', label: 'Length', unit: 'm', type: 'number' }, { key: 'b', label: 'Breadth', unit: 'm', type: 'number' }],
      default: DEFAULT_TANKS, help: 'Tanks that are neither full nor empty. Pressed-up or empty tanks carry no free surface and should not be listed.' },
    { key: 'Awind', label: 'Windage area of cargo', unit: 'm²', type: 'number', default: 1300, min: 0, max: 50000, group: 'Wind', role: 'geometry' },
    { key: 'Cs', label: 'Shape coefficient of the cargo', unit: '-', type: 'number', default: 1.1, min: 0.3, max: 2, group: 'Wind', role: 'load', help: 'IMO MODU Code table: 1.0 large flat surfaces, 1.25–1.3 open truss work, 0.5 cylinders.' },
    { key: 'Vwind', label: 'Wind speed for the overturning check', unit: 'm/s', type: 'number', default: 36, min: 5, max: 80, group: 'Wind', role: 'bc', help: `Boundary condition of the check: 36 m/s is the customary intact value for an unrestricted tow; weather-restricted operations may use a lower documented value. ${CONFIRM}` },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'gmReq', label: 'Minimum GM (after free-surface correction)', unit: 'm', type: 'number', default: 1.0, min: 0, max: 20, group: 'Criteria', role: 'criteria', help: CONFIRM },
    { key: 'rangeReq', label: 'Minimum range of positive stability', unit: 'deg', type: 'number', default: 36, min: 10, max: 90, group: 'Criteria', role: 'criteria', help: `About 36° for large barges and 40° for small ones on unrestricted voyages; lower with a documented motion analysis. ${CONFIRM}` },
    { key: 'ratioReq', label: 'Minimum righting / wind-heeling area ratio', unit: '-', type: 'number', default: 1.4, min: 1, max: 3, group: 'Criteria', role: 'criteria', help: CONFIRM },
    { key: 'imoSet', label: 'IMO criteria set', type: 'select', default: 'pontoon', group: 'Criteria', role: 'criteria', options: [{ value: 'pontoon', label: 'Pontoons / deck-cargo barges (IS Code Part B §2.2)' }, { value: 'ship', label: 'General ship criteria (IS Code Part A §2.2)' }],
      help: `Wide, shallow barges reach maximum GZ well before 25° and cannot meet the ship criteria; the pontoon criteria (area to the angle of maximum GZ ≥ 0.08 m·rad, range ≥ 20° at L ≤ 100 m falling to 15° at L ≥ 150 m) apply instead. ${CONFIRM}` },
    { key: 'trimMax', label: 'Allowable static trim', unit: 'deg', type: 'number', default: 0.5, min: 0.05, max: 5, group: 'Criteria', role: 'criteria' },
    { key: 'heelMax', label: 'Allowable static heel', unit: 'deg', type: 'number', default: 0.5, min: 0.05, max: 5, group: 'Criteria', role: 'criteria' },
    { key: 'pumpCap', label: 'Ballast pump capacity', unit: 'm³/h', type: 'number', default: 2000, min: 10, max: 1e5, group: 'Barge', role: 'operation' },
  ],
  compute(inp) {
    // t → kg.
    const rho = inp.rhoW, L = inp.Lb, B = inp.Bb, D = inp.Db;
    const parts = [
      { n: 'Lightship', m: inp.mLight * 1e3, x: 0, y: 0, z: inp.kgLight },
      { n: 'Cargo', m: inp.mCargo * 1e3, x: inp.xCargo, y: inp.yCargo, z: D + inp.zCargo },
      { n: 'Ballast', m: inp.mBallast * 1e3, x: 0, y: 0, z: inp.kgBallast },
    ];
    const M = parts.reduce((s, p) => s + p.m, 0);
    const cog = ['x', 'y', 'z'].map((k) => parts.reduce((s, p) => s + p.m * p[k], 0) / M);
    const imported = inp.hull && inp.hull.kind === 'surface' && inp.hull.positions && inp.hull.indices;
    const mesh = imported ? inp.hull : boxMesh({ xMin: -L / 2, xMax: L / 2, yMin: -B / 2, yMax: B / 2, zMin: 0, zMax: D });
    const Tbox = M / (rho * L * B);
    need(imported || Tbox < D, `The barge sinks: the loading condition of ${fmt(M / 1e3)} t needs a draft of ${fmt(Tbox)} m, more than the depth of ${D} m. Remove ballast or cargo.`);
    const hs = run(() => boxHydrostatics({ L, B, T: Tbox, KG: cog[2], rho }), 'Hydrostatics');
    // Free surface: GG′ = Σ ρ·l·b³/12 / (ρ·∇) for each slack tank.
    const slack = [];
    for (const t of inp.tanks || []) { const n = Math.round(Number(t.n) || 0), l = Number(t.l), b = Number(t.b); if (n > 0 && l > 0 && b > 0) for (let i = 0; i < n; i++) slack.push({ l, b, rho }); }
    const fs = slack.length ? run(() => freeSurfaceCorrection({ tanks: slack, volume: M / rho, rho }), 'Free-surface correction') : { GGt: 0, GGl: 0 };
    const KGc = cog[2] + fs.GGt, GM = hs.KMt - KGc, GMl = hs.KMl - cog[2] - fs.GGl;
    need(GM > 0, `The barge is initially unstable: GM = ${fmt(GM)} m after a free-surface correction of ${fmt(fs.GGt)} m. Lower the centre of gravity, press up slack tanks or use a wider barge.`);
    const eq = run(() => solveEquilibrium({ mesh, mass: M, cog: [cog[0], cog[1], KGc], rho }), 'Floating equilibrium');
    const draft = eq.draft, trim = eq.trim / DEG, heel = eq.heel / DEG;
    const heels = []; for (let h = 0; h <= 70; h += 2.5) heels.push(h);
    const gz = run(() => gzCurve({ mesh, mass: M, cog: [cog[0], 0, KGc], rho, heelsDeg: heels }), 'Righting-lever curve');
    const pr = gzProperties({ heelDeg: gz.heelDeg, gz: gz.gz });
    const imo = imoIntactCriteria({ heelDeg: gz.heelDeg, gz: gz.gz, gm0: GM });
    const range = pr.vanishingDeg ?? 70;
    const wl = windHeelingLever({ V: inp.Vwind, areas: [{ A: inp.Awind, z: D - draft + inp.zCargo, Cs: inp.Cs, Ch: 1 }, { A: L * Math.max(D - draft, 0), z: (D - draft) / 2, Cs: 1, Ch: 1 }], zLateralResistance: -draft / 2, displacement: M });
    let wind = null;
    try { wind = moduWindCriterion({ heelDeg: gz.heelDeg, gz: gz.gz, lever0: wl.lever, variation: 'constant', unitType: 'surface' }); } catch { wind = null; }
    const ratio = wind && Number.isFinite(wind.ratio) ? wind.ratio : 0;
    const edge = Math.atan((D - draft) / (B / 2)) / DEG;
    const gz10 = gz.gz[4];
    // Ballast to bring the barge level: fore/aft tank groups for trim, port/starboard transfer for heel.
    const capG = 0.4 * L * B * D;
    const vol0 = Math.min(inp.mBallast * 1e3 / rho / 2, capG);
    const plan = run(() => ballastPlan({ tanks: [{ name: 'Aft', x: -L / 4, y: 0, capacity: capG, volume: vol0, rho }, { name: 'Forward', x: L / 4, y: 0, capacity: capG, volume: vol0, rho }],
      hydro: { displacement: M, waterplaneArea: hs.waterplaneArea, lcf: 0, GMt: GM, GMl: Math.max(GMl, 1), rho }, current: { draft, heelDeg: 0, trimDeg: Math.abs(trim) * Math.sign(cog[0] || 1) }, target: { draft, trimDeg: 0 }, mode: 'transfer' }), 'Ballast plan');
    const vTrim = plan.volumeMoved, vHeel = M * Math.abs(cog[1]) / (rho * B / 2), vBal = vTrim + vHeel;
    const checks = [
      fsCheck('Initial stability GM', GM, inp.gmReq, `${N001} Sec. 11 (intact stability); IMO 2008 IS Code Part A §2.2.4`, 'Metacentric height after free-surface correction. Too little GM gives large heel under wind or cargo offset and long, large rolls.'),
      fsCheck('Range of positive stability', range, inp.rangeReq, `${N001} Sec. 11 (range of stability)`, 'Heel angle at which the righting lever vanishes. A short range means a large roll or wind gust can capsize the tow.'),
      fsCheck('Wind overturning: area ratio', ratio, inp.ratioReq, `${N001} Sec. 11; IMO MODU Code §3.3 (area under the righting curve ≥ 1.4 × area under the wind heeling curve)`, 'Energy reserve against a steady wind: righting area to the second intercept over the wind-heeling area. Below the required ratio a gust can roll the barge past its stable range.'),
      check('Static trim', Math.abs(trim), inp.trimMax, 'deg', 'Project loading-condition limits; DNV-ST-N001 Sec. 11 (loading condition)', 'Trim at rest from the longitudinal position of cargo and ballast. Excessive trim reduces freeboard at one end and upsets towing behaviour.'),
      check('Static heel', Math.abs(heel), inp.heelMax, 'deg', 'Project loading-condition limits', 'Heel at rest from transverse offsets. It consumes part of the righting-lever curve before any weather acts.'),
    ];
    if (inp.imoSet === 'ship') checks.push(...imo.checks.map((c) => fsCheck(`IMO: ${c.name}`, c.value, c.limit, c.ref || 'IMO 2008 IS Code Part A §2.2', 'One of the general intact-stability criteria on the shape of the righting-lever curve; together they guard against capsize in a seaway.')));
    else {
      const aMax = areaUnderGz(gz.heelDeg, gz.gz, 0, pr.heelAtGzMaxDeg), rPont = L <= 100 ? 20 : L >= 150 ? 15 : 20 - 5 * (L - 100) / 50;
      checks.push(fsCheck('IMO pontoon: area under GZ to the angle of maximum GZ', aMax, 0.08, 'IMO 2008 IS Code Part B §2.2 (pontoons)', 'Righting energy available up to the peak of the GZ curve, in m·rad. It guards against capsize of a wide, shallow hull whose stability falls away once the deck edge immerses.'),
        fsCheck('IMO pontoon: range of stability', range, rPont, 'IMO 2008 IS Code Part B §2.2 (pontoons)', 'Minimum range of positive righting levers for a deck-cargo pontoon of this length.'));
    }
    const gov = worst(checks), ok = gov.uc <= 1;
    // Lowest cargo CoG height at which GM just meets the requirement (GM is linear in the cargo height).
    const zMaxCargo = inp.mCargo > 0 ? inp.zCargo + (GM - inp.gmReq) * M / (inp.mCargo * 1e3) : 0;
    return {
      kpis: [
        kpi('disp', 'Displacement', M / 1e3, 't', 'Total mass of lightship, cargo and ballast: the weight the hull must float.'),
        kpi('draft', 'Draft at the equilibrium position', draft, 'm', 'Depth of the keel origin below the waterline at the floating equilibrium.', D, 'max'),
        kpi('trim', 'Static trim', trim, 'deg', 'Positive bow down. Corrected with fore-and-aft ballast transfer.'),
        kpi('heel', 'Static heel', heel, 'deg', 'From the transverse cargo offset; corrected by port–starboard ballast transfer.'),
        kpi('KG', 'KG corrected for free surface', KGc, 'm', 'Height of the centre of gravity above the keel including the virtual rise from slack tanks.'),
        kpi('GM', 'Transverse metacentric height GM', GM, 'm', 'KM − KG with free-surface correction. Positive is stable; very large values give a stiff barge with short roll period and high cargo accelerations.', inp.gmReq, 'min'),
        kpi('fsc', 'Free-surface correction GG′', fs.GGt, 'm', 'Virtual rise of the centre of gravity from liquid moving in slack tanks. Pressing tanks up removes it.'),
        kpi('gzMax', 'Maximum righting lever', pr.gzMax, 'm', 'Peak of the GZ curve; with the displacement it gives the largest heeling moment the barge can resist statically.'),
        kpi('gz10', 'Righting lever at 10° heel', gz10, 'm', 'For a wall-sided hull GZ = sin φ·(GM + ½·BM·tan² φ); used here to cross-check the mesh hydrostatics.'),
        kpi('range', 'Range of positive stability', range, 'deg', 'Angle at which GZ returns to zero (reported as 70° when it is still positive at the end of the curve).', inp.rangeReq, 'min'),
        kpi('ratio', 'Righting / wind-heeling area ratio', ratio, '-', 'Area under the righting curve over area under the wind-heeling lever up to the second intercept.', inp.ratioReq, 'min'),
        kpi('windHeel', 'Steady wind heel angle', wind?.firstInterceptDeg ?? 0, 'deg', 'First intercept of the wind-heeling lever with the GZ curve: the steady list under the design wind.'),
        kpi('edge', 'Deck-edge immersion angle', edge, 'deg', 'Heel at which the deck edge reaches the water; beyond it the waterplane and the righting lever fall off quickly.'),
        kpi('Troll', 'Natural roll period (kxx = 0.4·B)', 2 * Math.PI * 0.4 * B / Math.sqrt(G * GM), 's', 'Estimate from GM with a radius of gyration of 0.4·B including added inertia. Short periods (stiff barge) mean high transverse accelerations on the cargo.'),
        kpi('vBal', 'Ballast to transfer to level the barge', vBal, 'm³', 'Volume moved between fore/aft tank groups to remove trim plus the port–starboard transfer to remove heel.'),
        kpi('zMaxCargo', 'Highest cargo CoG above deck meeting the GM requirement', zMaxCargo, 'm', 'Linear extrapolation of GM with cargo height: the margin available for weight-report changes in the vertical centre of gravity.'),
      ],
      checks,
      plots: [
        line('Righting-lever curve and wind heeling lever', 'Heel angle', 'deg', 'Lever', 'm', [{ name: 'Righting lever GZ', x: Array.from(gz.heelDeg), y: Array.from(gz.gz).map((v) => round(v)) }, { name: `Wind heeling lever (${inp.Vwind} m/s)`, x: [0, 70], y: [round(wl.lever), round(wl.lever)], dash: true }, { name: 'Initial slope GM', x: [0, 57.3 * Math.min(1, pr.gzMax / GM)], y: [0, Math.min(GM, pr.gzMax)], dash: true }],
          'The GZ curve at constant displacement. Read the steady wind heel at the first crossing with the wind lever, the peak near deck-edge immersion and the angle where GZ returns to zero.'),
      ],
      tables: [
        tbl('Loading condition', [['n', 'Item'], ['m', 'Mass', 't'], ['x', 'x from midship', 'm'], ['y', 'y from centreline', 'm'], ['z', 'Height above keel', 'm']], [...parts.map((p) => ({ n: p.n, m: round(p.m / 1e3), x: round(p.x), y: round(p.y), z: round(p.z) })), { n: 'Total', m: round(M / 1e3), x: round(cog[0]), y: round(cog[1]), z: round(cog[2]) }],
          'Mass and centre of gravity of each item and of the whole barge. Update it with the final weight report and ballast plan.'),
        tbl('Ballast plan to level the barge', [['t', 'Tank group'], ['dv', 'Change of volume', 'm³'], ['h', 'Pumping time', 'h']], [...plan.tanks.map((t) => ({ t: t.name, dv: round(t.deltaVolume), h: round(Math.abs(t.deltaVolume) / inp.pumpCap, 3) })), { t: 'Port ↔ starboard', dv: round(vHeel), h: round(vHeel / inp.pumpCap, 3) }],
          'Transfers that remove the static trim and heel at constant displacement, from the small-angle hydrostatic model. Re-check the final condition after pumping.'),
      ],
      notes: [`Displacement ${fmt(M / 1e3)} t at ${fmt(draft)} m draft; GM ${fmt(GM)} m (free surface ${fmt(fs.GGt)} m).`, `Wind lever ${fmt(wl.lever)} m at ${inp.Vwind} m/s; area ratio ${fmt(ratio)}.`, `"${gov.name}" governs at ${fmt(gov.uc)}.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `Intact stability is met: GM ${fmt(GM)} m, range ${fmt(range)}°, wind area ratio ${fmt(ratio)}; trim ${fmt(trim)}° and heel ${fmt(heel)}° at rest.`,
          `Transfer about ${fmt(vBal)} m³ of ballast to level the barge, keep the listed tanks either pressed up or empty to hold the free-surface correction at ${fmt(fs.GGt)} m, and check the one-compartment damage case separately. The cargo CoG may rise to ${fmt(zMaxCargo)} m above deck before GM falls to the limit.`, `${N001} Sec. 11; IMO 2008 IS Code`)
          : adv('engineering', 'fail', `"${gov.name}" is not met (${fmt(gov.uc)}): GM ${fmt(GM)} m, range ${fmt(range)}°, area ratio ${fmt(ratio)}, trim ${fmt(trim)}°, heel ${fmt(heel)}°.`,
            gov.name.startsWith('Static') ? `Transfer about ${fmt(vBal)} m³ of ballast (${fmt(vBal / inp.pumpCap)} h of pumping) to level the barge, or reposition the cargo.` : 'Add low ballast, press up slack tanks, lower the cargo on its grillage or use a wider barge: stability grows with the square of the beam.',
            `${N001} Sec. 11; IMO 2008 IS Code`),
        adv('economic', 'info', `The barge carries ${fmt(inp.mBallast)} t of ballast for ${fmt(inp.mCargo)} t of cargo; levelling takes about ${fmt(vBal / inp.pumpCap)} h of pumping.`,
          GM > 4 * inp.gmReq ? 'GM is several times the requirement: the barge is stiffer than needed, which raises cargo accelerations and seafastening steel. Less ballast or a smaller barge class would reduce both draft and tow resistance.' : 'Stability margin is moderate: lock the cargo weight and vertical centre of gravity early, because a late rise would force extra ballast, a larger barge or a restricted weather window.',
          'Barge class selection and ballast voyage planning; DNV-ST-N001 Sec. 11 for the governing criteria; AACE RP 18R-97 for the estimate class',
          `Each 1000 t of ballast carried adds about ${fmt(1000 / (rho * L * B / 1e3) * 100)} cm of draft and tow resistance; a two-week tow burns about ${vesselCo2Text(14, 12)}`),
      ],
      publish: { 'install.bargeGM': GM, 'install.bargeDraft': draft, 'install.bargeDisp': M },
      quantities: {},
    };
  },
  verify: [
    { name: 'Box barge GM = T/2 + B²/(12·T) − KG', source: 'Wall-sided box (Biran & López-Pulido ch. 2): M = 19 950 t, T = 19.95e6/(1025·120·36) = 4.5054 m; KG = (5000·4.2 + 6450·21.6 + 8500·2.2)/19950 = 8.9735 m; GM = 2.2527 + 1296/54.065 − 8.9735 = 17.250 m (no slack tanks)',
      inputs: { tanks: [], xCargo: 0, yCargo: 0 }, expect: [{ kpi: 'GM', value: 19.95e6 / (1025 * 120 * 36) / 2 + 36 * 36 / (12 * 19.95e6 / (1025 * 120 * 36)) - (5000 * 4.2 + 6450 * 21.6 + 8500 * 2.2) / 19950, tol: 1e-6 }, { kpi: 'draft', value: 19.95e6 / (1025 * 120 * 36), tol: 1e-4 }] },
    { name: 'Wall-sided righting lever: GZ = sin φ·(GM + ½·BM·tan² φ) at 10°', source: 'Wall-sided formula (Biran & López-Pulido ch. 5), valid before deck-edge immersion (9.75° here, so the formula over-predicts by less than 1 %): BM = 23.971 m, GM = 17.250 m → GZ(10°) = 0.17365·(17.250 + 0.5·23.971·0.031091) = 3.060 m',
      inputs: { tanks: [], xCargo: 0, yCargo: 0 }, expect: [{ kpi: 'gz10', value: 3.060, tol: 0.01 }] },
    { name: 'Free-surface correction GG′ = Σ l·b³/12 / ∇', source: 'IMO 2008 IS Code Part B ch. 3: eight slack tanks 28 × 8.5 m → 8·28·8.5³/12 = 11 464 m⁴; ∇ = 19 463 m³ → GG′ = 0.5890 m',
      inputs: {}, expect: [{ kpi: 'fsc', value: 8 * 28 * 8.5 ** 3 / 12 / (19.95e6 / 1025), tol: 1e-6 }] },
  ],
};
