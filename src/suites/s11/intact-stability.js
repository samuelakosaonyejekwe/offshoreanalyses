// Intact stability at large angles: righting-lever curve by re-clipping the hull, IMO general criteria for
// ship-shaped units and the MODU wind-heeling area-ratio criterion with user wind areas and coefficients.
import { gzCurve, gzProperties, imoIntactCriteria, windHeelingLever, moduWindCriterion, areaUnderGz } from '../../engine/hydro/stability.js';
import { hullInputs, hullFrom, hullEnclosedVolume, need, run, round, fmt, adv, checkMin, check, sev, worst, steelCo2Text } from './_common.js';

const interp = (xs, ys, x) => { for (let i = 1; i < xs.length; i++) if (x <= xs[i]) return ys[i - 1] + (ys[i] - ys[i - 1]) * (x - xs[i - 1]) / (xs[i] - xs[i - 1]); return ys[ys.length - 1]; };

export default {
  id: 's11-intact-stability',
  suite: 11,
  title: 'Intact stability and righting-arm curve',
  summary: 'Righting-lever (GZ) curve of the hull at constant displacement, wind heeling lever, and the IMO intact criteria and MODU Code wind-heeling area ratio with the margin on KG.',
  covers: [1],
  fidelity: 'numerical',
  standards: ['IMO MODU Code 2009 §3.2–3.3', 'IMO 2008 IS Code (MSC.267(85)) Part A §2.2', 'DNV-OS-C301 (2021) Ch.2 Sec.1', 'ABS MODU Rules Part 3 Ch.3'],
  models: ['Righting-arm calculation at constant displacement (re-clipped hull, free trim optional)', 'Nonlinear hydrostatic equilibrium', 'Wind heeling lever from projected areas with shape and height coefficients', 'Area-ratio (dynamic stability) criterion', 'Free-surface correction as a virtual rise of G'],
  limits: 'Calm-water righting levers of the watertight hull only: no deck box buoyancy, no waves, no dynamic roll. Wind shape and height coefficients are user inputs taken from the MODU Code tables or a wind-tunnel test; the heeling lever variation with heel is a chosen law. The IS Code weather criterion (§2.3) is not included. A tension-leg platform is assessed free-floating (tow or installation condition).',
  searchTerms: 'intact stability righting arm MODU wind heeling criterion',
  inputs: [
    ...hullInputs({ type: 'semi' }),
    { key: 'GGfs', label: 'Free-surface rise of G', unit: 'm', type: 'number', default: 0.25, min: 0, max: 5, group: 'Loading condition', role: 'ic', help: 'Virtual rise of the centre of gravity from slack tanks (from the hydrostatics module). Added to KG for the whole curve.' },
    { key: 'freeTrim', label: 'Free trim while heeling', type: 'bool', default: false, group: 'Loading condition', role: 'bc', help: 'Re-balance trim at every heel angle (required for asymmetric hulls; slower).' },
    { key: 'unitType', label: 'Unit type for the wind criterion', type: 'select', default: 'column-stabilised', group: 'Wind heeling', role: 'criteria',
      options: [{ value: 'column-stabilised', label: 'Column-stabilised (area ratio 1.3)' }, { value: 'surface', label: 'Surface unit / ship-shaped (area ratio 1.4)' }, { value: 'self-elevating', label: 'Self-elevating afloat (area ratio 1.4)' }], help: 'Sets the required ratio of righting to heeling energy.' },
    { key: 'Vwind', label: 'Wind speed for the heeling lever', unit: 'm/s', type: 'number', default: 51.5, min: 0, max: 80, group: 'Wind heeling', role: 'environment',
      help: 'MODU Code: 36 m/s (70 kn) normal operating and transit, 51.5 m/s (100 kn) severe storm.' },
    { key: 'areas', label: 'Wind areas above the waterline', type: 'table', group: 'Wind heeling', role: 'load',
      columns: [{ key: 'name', label: 'Block', type: 'text' }, { key: 'A', label: 'Projected area', unit: 'm²', type: 'number' }, { key: 'z', label: 'Centroid above waterline', unit: 'm', type: 'number' }, { key: 'Cs', label: 'Shape coeff.', type: 'number' }, { key: 'Ch', label: 'Height coeff.', type: 'number' }],
      default: [{ name: 'Columns above water', A: 540, z: 9, Cs: 0.5, Ch: 1.0 }, { name: 'Deck box', A: 720, z: 22, Cs: 1.0, Ch: 1.1 }, { name: 'Topsides and derrick', A: 1500, z: 42, Cs: 1.1, Ch: 1.2 }],
      help: 'Beam-on projected areas with MODU Code §3.2.3 shape (0.4 sphere, 0.5 cylinder, 1.0 flat, 1.25 derrick, 1.1 clustered deckhouses) and height coefficients (1.0 to 15.3 m, 1.1 to 30.5 m, 1.2 to 46 m, 1.3 to 61 m). They are user inputs.' },
    { key: 'zLat', label: 'Centre of lateral resistance below waterline', unit: 'm', type: 'number', default: 11, min: 0, max: 150, group: 'Wind heeling', role: 'load', help: 'Depth of the centroid of the underwater lateral area (about half the draft), the lower end of the heeling arm.' },
    { key: 'variation', label: 'Heeling lever variation with heel', type: 'select', default: 'cos2', group: 'Wind heeling', role: 'criteria',
      options: [{ value: 'constant', label: 'Constant (conservative)' }, { value: 'cos', label: 'cos φ' }, { value: 'cos2', label: 'cos² φ (ship-shaped, MODU §3.2.5)' }], help: 'How the wind lever reduces as the unit heels.' },
    { key: 'dfX', label: 'Downflooding opening x', unit: 'm', type: 'number', default: 0, min: -250, max: 250, group: 'Downflooding', role: 'bc', help: 'Longitudinal position of the lowest opening that cannot be closed weathertight.' },
    { key: 'dfY', label: 'Downflooding opening y (off centreline)', unit: 'm', type: 'number', default: 30, min: 0, max: 80, group: 'Downflooding', role: 'bc', help: 'Transverse position; the curve is heeled towards this side.' },
    { key: 'dfZ', label: 'Downflooding opening height above keel', unit: 'm', type: 'number', default: 48, min: 1, max: 250, group: 'Downflooding', role: 'bc', help: 'Height of the opening. Its immersion ends the usable range of stability.' },
    { key: 'heelMaxStatic', label: 'Allowable static wind heel', unit: 'deg', type: 'number', default: 15, min: 2, max: 30, group: 'Acceptance criteria', role: 'criteria', help: 'Project limit on the steady heel at the first intercept (not a MODU Code number).' },
    { key: 'heelRange', label: 'Largest heel angle of the curve', unit: 'deg', type: 'number', default: 60, min: 30, max: 85, group: 'Numerical controls', role: 'control', help: 'Upper end of the righting-lever curve.' },
    { key: 'dHeel', label: 'Heel step', unit: 'deg', type: 'number', default: 2.5, min: 0.5, max: 5, group: 'Numerical controls', role: 'control', help: 'Spacing of the computed points; 2.5° or finer for criteria based on areas.' },
  ],
  mesh: { key: 'dHeel', label: 'Heel step', levels: [5, 2.5, 1.25], qoi: 'area40', isSize: true },
  compute(inp) {
    const hull = hullFrom(inp, { nTheta: 24 });
    const cog = [hull.cog[0], hull.cog[1], inp.KG + inp.GGfs];
    const gm0 = hull.hs.KMt - inp.KG - inp.GGfs;
    need(gm0 > 0, `Initially unstable: GM = ${fmt(gm0)} m after the free-surface rise. KG must be below ${fmt(hull.hs.KMt - inp.GGfs)} m before a righting-lever curve is meaningful.`);
    const heels = []; for (let h = 0; h <= inp.heelRange + 1e-9; h += inp.dHeel) heels.push(round(h, 6));
    // heel towards the side of the opening (+y port opening immerses when port goes down → side −1 = starboard up)
    const curve = run(() => gzCurve({ mesh: hull.mesh, mass: hull.mass, cog, rho: hull.rho, heelsDeg: heels, side: -1, freeTrim: inp.freeTrim, openings: [{ name: 'downflooding', point: [inp.dfX, inp.dfY, inp.dfZ] }] }), 'Righting-lever curve');
    need(curve.gz.every(Number.isFinite), 'The righting-lever curve could not be evaluated over the whole heel range: reduce the largest heel angle.');
    const df = curve.downfloodDeg;
    const prop = gzProperties({ heelDeg: curve.heelDeg, gz: curve.gz, downfloodDeg: df });
    const areas = (inp.areas || []).filter((a) => a.A > 0).map((a) => ({ A: a.A, z: a.z, Cs: a.Cs > 0 ? a.Cs : 1, Ch: a.Ch > 0 ? a.Ch : 1 }));
    const wind = areas.length && inp.Vwind > 0 ? windHeelingLever({ V: inp.Vwind, areas, zLateralResistance: -inp.zLat, displacement: hull.mass }) : { force: 0, moment: 0, lever: 0 };
    const modu = run(() => moduWindCriterion({ heelDeg: curve.heelDeg, gz: curve.gz, lever0: Math.max(wind.lever, 1e-6), variation: inp.variation, unitType: inp.unitType, downfloodDeg: df }), 'Wind-heeling criterion');
    const imo = imoIntactCriteria({ heelDeg: curve.heelDeg, gz: curve.gz, gm0, downfloodDeg: df });
    const shipLike = inp.unitType === 'surface';
    const first = modu.firstInterceptDeg ?? 0;
    const ratio = Number.isFinite(modu.ratio) ? Math.min(modu.ratio, 999) : 999;
    const checks = [
      checkMin(`MODU wind-heeling area ratio (to ${fmt(modu.limitDeg)}°)`, ratio, modu.required, '-', 'IMO MODU Code 2009 §3.3.1',
        'Guards against capsize in a gust: the righting energy up to the second intercept or downflooding must exceed the wind heeling energy by the required factor, otherwise a dynamic roll can carry the unit past its range of stability.'),
      check('Static heel under the design wind', first, inp.heelMaxStatic, 'deg', 'Project criterion; MODU Code 2009 §3.3 first intercept',
        'Guards against an excessive steady list in the storm, which immerses openings and overloads equipment foundations even when the energy criterion is met.'),
      checkMin('Positive range of stability', prop.rangeDeg ?? (inp.heelRange - prop.equilibriumDeg), Math.max(modu.secondInterceptDeg ?? 0, first) > 0 ? Math.min(modu.limitDeg, inp.heelRange) : 10, 'deg', 'IMO MODU Code 2009 §3.3.2',
        'The righting lever must stay positive from upright to the second intercept; a curve that vanishes earlier means the unit capsizes before the wind lever is overcome.'),
    ];
    if (shipLike) for (const c of imo.checks) checks.push(checkMin(c.name, c.value, c.limit, c.unit, c.ref, 'IMO general intact criterion for ships: a minimum of righting energy, lever and initial stiffness established from casualty statistics; falling short means insufficient reserve against wave and wind heeling.'));
    const gov = worst(checks);
    const lever = curve.heelDeg.map((h) => wind.lever * (inp.variation === 'cos2' ? Math.cos(h * Math.PI / 180) ** 2 : inp.variation === 'cos' ? Math.cos(h * Math.PI / 180) : 1));
    // allowable KG: raising G by d lowers GZ by d·sinφ; find d at which the area ratio just meets the requirement
    let lo = -5; let hi = Math.min(gm0 - 0.02, 30);
    const ratioAt = (d) => { const gz = curve.gz.map((v, i) => v - d * Math.sin(curve.heelDeg[i] * Math.PI / 180)); try { const m = moduWindCriterion({ heelDeg: curve.heelDeg, gz, lever0: Math.max(wind.lever, 1e-6), variation: inp.variation, unitType: inp.unitType, downfloodDeg: df }); return m.pass ? 1 : -1; } catch { return -1; } };
    let kgMargin;
    if (ratioAt(hi) > 0) kgMargin = hi; else if (ratioAt(lo) < 0) kgMargin = lo; else { for (let i = 0; i < 30; i++) { const mid = 0.5 * (lo + hi); if (ratioAt(mid) > 0) lo = mid; else hi = mid; } kgMargin = lo; }
    const steel = 0.16 * hullEnclosedVolume(hull);
    return {
      kpis: [
        { key: 'gm0', label: 'Initial metacentric height (corrected)', value: round(gm0, 5), unit: 'm', explain: 'KM − KG − free-surface rise: the slope of the righting-lever curve at the origin.' },
        { key: 'gzMax', label: 'Maximum righting lever', value: round(prop.gzMax, 5), unit: 'm', explain: 'Peak of the GZ curve: the largest steady heeling lever the unit can resist.' },
        { key: 'heelGzMax', label: 'Heel at maximum GZ', value: round(prop.heelAtGzMaxDeg, 4), unit: 'deg', explain: 'Angle of the GZ peak; beyond it the righting moment falls. Ships must reach it at 25° or more.' },
        { key: 'gz10', label: 'Righting lever at 10°', value: round(interp(curve.heelDeg, curve.gz, 10), 5), unit: 'm', explain: 'Small-angle point of the curve, close to GM·sin 10° plus the form contribution; used to check the curve against the wall-sided formula.' },
        { key: 'range', label: 'Range of positive stability', value: round(prop.vanishingDeg ?? inp.heelRange, 4), unit: 'deg', explain: 'Angle at which GZ returns to zero (the computed range if it never does). Past it the unit capsizes.' },
        { key: 'area30', label: 'Area under GZ, 0–30°', value: round(prop.area0to30, 5), unit: 'm·rad', explain: 'Righting energy per unit weight up to 30°; the IS Code asks 0.055 m·rad for ships.' },
        { key: 'area40', label: 'Area under GZ, 0–40° or downflooding', value: round(prop.area0to40, 5), unit: 'm·rad', explain: 'Righting energy up to 40° or the downflooding angle if smaller; the IS Code asks 0.090 m·rad for ships.' },
        { key: 'windLever', label: 'Upright wind heeling lever', value: round(wind.lever, 5), unit: 'm', explain: 'Wind overturning moment divided by the weight of the unit, at zero heel.' },
        { key: 'windForce', label: 'Wind force', value: round(wind.force / 1e3, 5), unit: 'kN', explain: 'Total wind force on the listed blocks, ½·ρair·Cs·Ch·A·V².' },
        { key: 'heelWind', label: 'Static wind heel (first intercept)', value: round(first, 4), unit: 'deg', limit: inp.heelMaxStatic, sense: 'max', explain: 'Steady heel where the righting lever equals the wind lever.' },
        { key: 'heel2', label: 'Second intercept', value: round(modu.secondInterceptDeg ?? inp.heelRange, 4), unit: 'deg', explain: 'Angle where the falling GZ curve meets the wind lever again (the end of the computed range when it does not): the dynamic capsize boundary.' },
        { key: 'areaRatio', label: 'Righting / heeling area ratio', value: round(ratio, 4), unit: '-', limit: modu.required, sense: 'min', explain: 'Area under GZ divided by area under the wind lever up to the second intercept or downflooding. The MODU Code requires 1.3 for column-stabilised and 1.4 for surface units.' },
        { key: 'downflood', label: 'Downflooding angle', value: round(df ?? inp.heelRange, 4), unit: 'deg', explain: 'Heel at which the opening reaches the waterline (the end of the computed range if it stays dry). It truncates the usable righting energy.' },
        { key: 'kgMargin', label: 'KG margin on the wind criterion', value: round(kgMargin, 4), unit: 'm', explain: 'How far the centre of gravity can rise (negative: must come down) before the area-ratio criterion is just met. Allowable KG = KG + this margin.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Righting and wind heeling levers', x: { label: 'Heel', unit: 'deg' }, y: { label: 'Lever', unit: 'm' },
          series: [{ name: 'Righting lever GZ', x: curve.heelDeg, y: curve.gz.map((v) => round(v, 5)) }, { name: 'Wind heeling lever', x: curve.heelDeg, y: lever.map((v) => round(v, 5)), dash: true },
            ...(df != null ? [{ name: 'Downflooding', x: [round(df), round(df)], y: [0, round(prop.gzMax)], dash: true }] : [])],
          explain: 'The GZ curve against the wind lever. The first crossing is the steady heel, the second the capsize boundary; the area between the curves from first to second crossing is the reserve against a gust. A vertical line marks downflooding.' },
        { type: 'line', title: 'Draft and opening height while heeling', x: { label: 'Heel', unit: 'deg' }, y: { label: 'Height', unit: 'm' },
          series: [{ name: 'Draft at origin', x: curve.heelDeg, y: curve.draft.map((v) => round(v, 5)) }, { name: 'Opening above waterline', x: curve.heelDeg, y: curve.openingHeight.map((v) => round(v, 5)) }],
          explain: 'Constant-displacement draft and the height of the downflooding opening above the inclined waterline. Where the opening height crosses zero the unit starts to flood.' },
      ],
      tables: [
        { title: 'Righting-lever curve', columns: [{ key: 'h', label: 'Heel', unit: 'deg' }, { key: 'gz', label: 'GZ', unit: 'm' }, { key: 'lev', label: 'Wind lever', unit: 'm' }, { key: 'area', label: 'Area from 0', unit: 'm·rad' }, { key: 'open', label: 'Opening height', unit: 'm' }],
          rows: curve.heelDeg.map((h, i) => ({ h, gz: round(curve.gz[i], 5), lev: round(lever[i], 5), area: round(areaUnderGz(curve.heelDeg, curve.gz, 0, h), 5), open: round(curve.openingHeight[i], 4) })).filter((_, i) => i % Math.ceil(curve.heelDeg.length / 60) === 0),
          explain: 'Cross-curve data at the loaded displacement: one row per heel angle with the dynamic stability (area) accumulated from upright.' },
        { title: 'IMO IS Code general criteria', columns: [{ key: 'name', label: 'Criterion' }, { key: 'value', label: 'Value' }, { key: 'limit', label: 'Required' }, { key: 'unit', label: 'Unit' }, { key: 'ok', label: 'Met' }],
          rows: imo.checks.map((c) => ({ name: c.name, value: round(c.value, 5), limit: c.limit, unit: c.unit, ok: c.ok ? 'yes' : 'no' })),
          explain: 'The six general intact criteria for ships. They are enforced as checks for surface (ship-shaped) units and shown for information otherwise.' },
      ],
      notes: [
        `GZ peaks at ${fmt(prop.gzMax)} m at ${fmt(prop.heelAtGzMaxDeg)}°; ${df != null ? `downflooding at ${fmt(df)}° limits the usable range` : 'the opening stays dry over the computed range'}.`,
        `Wind ${fmt(inp.Vwind)} m/s gives ${fmt(wind.force / 1e3)} kN and a lever of ${fmt(wind.lever)} m; area ratio ${fmt(ratio)} against ${modu.required} required.`,
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at utilisation ${fmt(gov.uc)}; the wind criterion leaves a KG margin of ${fmt(kgMargin)} m.`,
          gov.uc > 1 ? `Lower KG by at least ${fmt(Math.max(-kgMargin, 0.05))} m, raise or make weathertight the downflooding opening (now immersing at ${fmt(df ?? inp.heelRange)}°), or reduce windage high up; then repeat for the damaged condition.`
            : `Publish allowable KG = ${fmt(inp.KG + kgMargin)} m for this draft in the operating manual, and repeat over the draft range (transit, operating, survival) and for the damaged condition, which often governs.`,
          'IMO MODU Code 2009 §3.3 (intact stability criteria) and §3.2.3 (wind force coefficients); IMO 2008 IS Code Part A §2.2.'),
        adv('economic', 'info', `A KG margin of ${fmt(kgMargin)} m corresponds to roughly ${fmt(Math.max(kgMargin, 0) * hull.mass / 1e3 / Math.max(hull.depth + 10 - inp.KG, 5))} t of additional variable deck load at 10 m above the hull top before the criterion is reached.`,
          'Use the margin deliberately: either carry more deck load per trip (fewer supply runs) or keep it as growth allowance for future topsides modules instead of re-qualifying the hull later.',
          'Weight and stability control per ISO 19901-5 (2021); whole-life costing per ISO 15686-5.',
          `Meeting a shortfall with larger columns adds steel: ${steelCo2Text(0.05 * steel)} for a 5 % hull-steel increase. Lowering windage or KG achieves the same margin without added material.`),
      ],
    };
  },
  verify: [
    { // Box 100 × 20 × 12 m, T = 5 m, KG = 4 m: GM = 5.1667, BM = 6.6667. Wall-sided formula (valid below bilge emergence at
      // 26.6°): GZ(10°) = sin10°·(GM + ½·BM·tan²10°) = 0.173648·(5.16667 + 0.5·6.66667·0.031091) = 0.9152 m.
      name: 'Wall-sided formula GZ = sin φ (GM + ½ BM tan² φ) for a box barge at 10°',
      source: 'Biran & López-Pulido, Ship Hydrostatics and Stability (2014) §5.3 — wall-sided formula',
      inputs: { hullType: 'barge', L: 100, B: 20, Dhull: 12, draft: 5, KG: 4, GGfs: 0, unitType: 'surface', dfY: 10, dfZ: 12, heelRange: 60 },
      expect: [{ kpi: 'gz10', value: 0.9152, tol: 2e-3 }, { kpi: 'gm0', value: 5.1667, tol: 1e-4 }, { kpi: 'downflood', value: 34.99, tol: 0.03 }] },
    { // Wind lever by hand for one block: F = ½·1.222·1.0·1.0·1000·40² = 977 600 N; arm = 10 + 4 = 14 m;
      // lever = F·arm/(Δ g) = 13 686 400/(10 250 000·9.80665) = 0.13616 m.
      name: 'Wind heeling lever of a single block',
      source: 'IMO MODU Code 2009 §3.2.3, F = ½ Cs Ch ρ V² A — hand calculation',
      inputs: { hullType: 'barge', L: 100, B: 20, Dhull: 12, draft: 5, KG: 4, GGfs: 0, Vwind: 40, zLat: 4, areas: [{ name: 'block', A: 1000, z: 10, Cs: 1, Ch: 1 }], dfY: 10, dfZ: 12 },
      expect: [{ kpi: 'windForce', value: 977.6, tol: 1e-3 }, { kpi: 'windLever', value: 0.13616, tol: 1e-3 }] },
  ],
};
