// Gravity-based foundation for an offshore wind turbine: overturning, sliding and bearing capacity of a circular
// base (with optional skirts) under factored turbine and environmental loads, contact pressure, settlement and tilt.
import { effectiveAreaCircle, overturningCheck, slidingResistance, isoUndrainedCapacity, apiDrainedCapacity, bearingCapacity, terzaghiBearing } from '../../engine/geotech/bearing.js';
import { elasticSettlement, footingRotation, consolidationSettlement, consolidationU } from '../../engine/geotech/settlement.js';
import { mccUndrainedStrengthRatio } from '../../engine/geotech/soilmodels.js';
import { G, DEG, need, run, fmt, round, cap, adv, check, worst, sev, concreteCo2Text, LIMIT_CODE, LIMIT_HIFI, CONFIRM } from './_common.js';

export default {
  id: 's03-gravity-base-stability', suite: 3, title: 'Gravity-base overturning, sliding and bearing capacity',
  summary: 'Stability of a circular gravity-based foundation: eccentricity and overturning, sliding at the base or skirt tip, bearing capacity on the effective area under inclined eccentric load, contact pressure, settlement and tilt.',
  covers: [8], fidelity: 'analytical',
  standards: ['DNV-ST-0126 (2021) §7.5 and Appendix G (gravity-based foundations, effective area)', 'ISO 19901-4:2016 §7 (shallow foundations)', 'API RP 2GEO (2011) §7', 'DNV-RP-C212 (2019) §5', 'IEC 61400-3-1:2019 §7.6', 'Terzaghi (1943); Meyerhof (1963); Brinch Hansen (1970); Vesic (1975)'],
  models: ['Force and moment equilibrium of a rigid base', 'Meyerhof effective area for an eccentrically loaded circular footing', 'Undrained bearing capacity with shape, depth and inclination factors (ISO 19901-4)', 'Drained bearing capacity with load-inclination factors (API RP 2A / Vesic)', 'Terzaghi and Meyerhof bearing-capacity formulations (comparison)', 'Mohr–Coulomb sliding with skirt passive resistance', 'Elastic settlement and rotation of a rigid circular footing', 'Terzaghi one-dimensional consolidation theory (settlement and time rate)', 'Modified Cam-Clay undrained strength ratio (cross-check of the strength profile)'],
  limits: `Rigid base on a uniform soil (undrained clay with strength increasing with depth, or drained sand); layered soils, punch-through, cyclic strength degradation and pore-pressure build-up are not included; consolidation is one-dimensional on a single clay layer. Hydrodynamic loads on the base itself must be included in the applied loads. ${LIMIT_CODE} ${LIMIT_HIFI}`,
  searchTerms: 'gravity based foundation offshore wind bearing capacity sliding overturning effective area',
  inputs: [
    { key: 'Db', label: 'Base diameter', unit: 'm', type: 'number', default: 32, min: 8, max: 60, group: 'Base geometry', role: 'geometry', design: true },
    { key: 'skirt', label: 'Skirt depth', unit: 'm', type: 'number', default: 0.8, min: 0, max: 8, group: 'Base geometry', role: 'geometry', design: true, help: 'Skirts move the sliding plane below weak surface soil and resist scour undermining.' },
    { key: 'concrete', label: 'Concrete', type: 'select', source: 'materials:concrete', default: 'c45', group: 'Base geometry', role: 'material' },
    { key: 'Vconc', label: 'Concrete volume', unit: 'm³', type: 'number', default: 2600, min: 100, max: 30000, group: 'Base geometry', role: 'geometry', design: true, help: 'Base slab, shaft and ice cone, all taken as submerged (conservative for stability).' },
    { key: 'ballastT', label: 'Ballast mass', unit: 't', type: 'number', default: 9000, min: 0, max: 80000, group: 'Base geometry', role: 'geometry', design: true, help: 'Sand or iron-ore ballast placed inside the base after set-down.' },
    { key: 'rhoBallast', label: 'Ballast saturated density', unit: 'kg/m³', type: 'number', default: 2000, min: 1300, max: 4500, group: 'Base geometry', role: 'material' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'Vtop', label: 'Vertical load from turbine and tower', unit: 'MN', type: 'number', default: 7, min: 0.5, max: 60, group: 'Loads at base level', role: 'bc', help: 'Weight of rotor–nacelle assembly, tower and interface steel transferred at the tower base.' },
    { key: 'H', label: 'Horizontal load (characteristic)', unit: 'MN', type: 'number', default: 6.5, min: 0, max: 200, group: 'Loads at base level', role: 'bc', bind: 'owt.mudlineShear', uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Wind, wave and current resultant including the wave load on the base itself.' },
    { key: 'M', label: 'Overturning moment at base level (characteristic)', unit: 'MN·m', type: 'number', default: 300, min: 0, max: 8000, group: 'Loads at base level', role: 'bc', bind: 'owt.mudlineMoment', uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'uplift', label: 'Wave uplift on the base (characteristic)', unit: 'MN', type: 'number', default: 4, min: 0, max: 200, group: 'Loads at base level', role: 'bc', help: 'Reduction of the vertical load by wave pressure under and over the base slab at the instant of maximum horizontal load.' },
    { key: 'soil', label: 'Soil response', type: 'select', default: 'sand', group: 'Soil', role: 'soil', options: [{ value: 'sand', label: 'Sand — drained (φ′)' }, { value: 'clay', label: 'Clay — undrained (su)' }] },
    { key: 'phiDeg', label: 'Effective friction angle φ′', unit: 'deg', type: 'number', default: 35, min: 20, max: 45, group: 'Soil', role: 'soil', bind: 'site.phiDeg', showIf: { key: 'soil', in: ['sand'] }, uncertain: { dist: 'normal', cov: 0.06 } },
    { key: 'su0', label: 'Undrained shear strength at mudline', unit: 'kPa', type: 'number', default: 80, min: 2, max: 600, group: 'Soil', role: 'soil', bind: 'site.su0', showIf: { key: 'soil', in: ['clay'] }, uncertain: { dist: 'lognormal', cov: 0.25 } },
    { key: 'suGrad', label: 'Strength gradient', unit: 'kPa/m', type: 'number', default: 2, min: 0, max: 30, group: 'Soil', role: 'soil', bind: 'site.suGrad', showIf: { key: 'soil', in: ['clay'] } },
    { key: 'gammaSub', label: 'Submerged unit weight γ′', unit: 'kN/m³', type: 'number', default: 10, min: 4, max: 13, group: 'Soil', role: 'soil', bind: 'site.gammaSub' },
    { key: 'Es', label: 'Soil Young\'s modulus', unit: 'MPa', type: 'number', default: 60, min: 1, max: 2000, group: 'Soil', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.4 }, help: 'Secant modulus at working stress for settlement and tilt (undrained Eu with ν = 0.5 for clay).' },
    { key: 'Hclay', label: 'Thickness of the compressible clay below the base', unit: 'm', type: 'number', default: 12, min: 0.5, max: 100, group: 'Consolidation (clay)', role: 'soil', showIf: { key: 'soil', in: ['clay'] } },
    { key: 'Cc', label: 'Compression index Cc', unit: '-', type: 'number', default: 0.25, min: 0.02, max: 2, group: 'Consolidation (clay)', role: 'soil', showIf: { key: 'soil', in: ['clay'] }, uncertain: { dist: 'lognormal', cov: 0.25 } },
    { key: 'e0', label: 'Initial void ratio', unit: '-', type: 'number', default: 0.9, min: 0.3, max: 3, group: 'Consolidation (clay)', role: 'ic', showIf: { key: 'soil', in: ['clay'] }, help: 'In-situ state of the clay before loading.' },
    { key: 'ocr', label: 'Overconsolidation ratio', unit: '-', type: 'number', default: 3, min: 1, max: 30, group: 'Consolidation (clay)', role: 'ic', showIf: { key: 'soil', in: ['clay'] }, help: 'Preconsolidation pressure over in-situ effective stress; recompression (Cr = Cc/6) applies below it.' },
    { key: 'cv', label: 'Coefficient of consolidation', unit: 'm²/yr', type: 'number', default: 4, min: 0.1, max: 200, group: 'Consolidation (clay)', role: 'soil', showIf: { key: 'soil', in: ['clay'] } },
    { key: 'phiCs', label: 'Critical-state friction angle of the clay', unit: 'deg', type: 'number', default: 28, min: 18, max: 40, group: 'Consolidation (clay)', role: 'soil', showIf: { key: 'soil', in: ['clay'] }, help: 'Used with Modified Cam-Clay to cross-check the undrained strength entered against the stress history.' },
    { key: 'rough', label: 'Base interface roughness tanδ / tanφ′', unit: '-', type: 'number', default: 0.8, min: 0.3, max: 1, group: 'Soil', role: 'soil', help: 'About 1.0 with skirts or a cast-in-place rough base, 0.6–0.8 for a smooth precast base on a gravel bed.' },
    { key: 'gammaF', label: 'Load factor on environmental loads', unit: '-', type: 'number', default: 1.35, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: `IEC 61400-3-1 normal design situations. ${CONFIRM}` },
    { key: 'gammaG', label: 'Load factor on favourable weight', unit: '-', type: 'number', default: 0.9, min: 0.7, max: 1.1, group: 'Acceptance criteria', role: 'criteria', help: `Permanent load acting favourably. ${CONFIRM}` },
    { key: 'gammaM', label: 'Material factor on soil strength', unit: '-', type: 'number', default: 1.25, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: `DNV-ST-0126: 1.25 on su (total stress), 1.15 on tanφ′ (effective stress). ${CONFIRM}` },
    { key: 'tiltLim', label: 'Allowable permanent tilt', unit: 'deg', type: 'number', default: 0.25, min: 0.05, max: 1, group: 'Acceptance criteria', role: 'criteria' },
  ],
  compute(inp, ctx) {
    const conc = ctx && ctx.material ? ctx.material(inp.concrete) : null;
    const rhoC = conc?.rhoReinforced ?? conc?.rho ?? 2500;
    const R = inp.Db / 2, A = Math.PI * R * R;
    need(inp.rhoBallast > inp.rhoW, 'Ballast density must exceed the water density.');
    // Submerged weights (N): concrete and ballast, both below water.
    const Wc = inp.Vconc * (rhoC - inp.rhoW) * G, Wb = inp.ballastT * 1e3 * (1 - inp.rhoW / inp.rhoBallast) * G;
    const Vk = inp.Vtop * 1e6 + Wc + Wb - inp.uplift * 1e6;
    need(Vk > 0, 'The wave uplift exceeds the submerged weight: the base would lift off. Add ballast or concrete.');
    const Hk = inp.H * 1e6, Mk = inp.M * 1e6 + Hk * inp.skirt; // moment transferred to skirt-tip level
    const Vd = inp.gammaG * (inp.Vtop * 1e6 + Wc + Wb) - inp.gammaF * inp.uplift * 1e6, Hd = inp.gammaF * Hk, Md = inp.gammaF * Mk;
    need(Vd > 0, 'The factored vertical load is not compressive: add ballast or reduce uplift.');
    const ot = overturningCheck({ V: Vk, M: Mk, B: inp.Db, shape: 'circle' });
    const ed = Md / Vd;
    const clay = inp.soil === 'clay';
    const phi = clay ? 0 : Math.atan(Math.tan(inp.phiDeg * DEG) / inp.gammaM), gam = inp.gammaSub * 1e3;
    const suBase = (inp.su0 + inp.suGrad * inp.skirt) * 1e3;
    let Qd = 0, Aeff = 0, Beff = 0, Leff = 0, Hres = 0, meyer = 0;
    const overturned = ed >= 0.98 * R;
    if (!overturned) {
      const ea = effectiveAreaCircle({ R, e: Math.max(ed, 0) });
      Aeff = ea.Aeff; Beff = ea.Beff; Leff = ea.Leff;
      if (clay) {
        Qd = run(() => isoUndrainedCapacity({ Beff, Leff, su0: suBase, kappa: inp.suGrad * 1e3, D: inp.skirt, H: Math.min(Hd, 0.999 * Aeff * suBase / inp.gammaM), materialFactor: inp.gammaM }), 'Undrained bearing capacity').Qd + gam * inp.skirt * Aeff;
        Hres = slidingResistance({ mode: 'undrained', A: Aeff, suBase, skirtDepth: inp.skirt, skirtWidth: inp.Db, suSkirt: (inp.su0 + 0.5 * inp.suGrad * inp.skirt) * 1e3, gamma: gam, materialFactor: inp.gammaM }).H;
        try { meyer = bearingCapacity({ method: 'meyerhof', B: Beff, L: Leff, D: inp.skirt, c: suBase / inp.gammaM, phi: 0, gamma: gam, V: Vd, HB: Math.min(Hd, 0.99 * Vd) }).Qult; } catch { meyer = 0; }
      } else {
        Qd = run(() => apiDrainedCapacity({ Beff, Leff, phi, gamma: gam, q: gam * inp.skirt, D: inp.skirt, V: Vd, H: Math.min(Hd, 0.95 * Vd * Math.tan(phi)) }), 'Drained bearing capacity').Q;
        const delta = Math.atan(inp.rough * Math.tan(phi));
        Hres = slidingResistance({ mode: 'drained', A: Aeff, V: Vd, delta: inp.skirt > 0 ? phi : delta, phi, gamma: gam, skirtDepth: inp.skirt, skirtWidth: inp.Db }).H;
        try { meyer = bearingCapacity({ method: 'meyerhof', B: Beff, L: Leff, D: inp.skirt, phi, gamma: gam, V: Vd, HB: Math.min(Hd, 0.9 * Vd) }).Qult; } catch { meyer = 0; }
      }
    }
    let terz = 0;
    try { terz = terzaghiBearing({ shape: 'circular', B: inp.Db, c: clay ? suBase : 0, phi: clay ? 0 : inp.phiDeg * DEG, gamma: gam, q: gam * inp.skirt }).qult * A; } catch { terz = 0; }
    const qMax = Aeff > 0 ? Vd / Aeff : Infinity;
    const E = inp.Es * 1e6, nu = clay ? 0.5 : 0.3;
    const settle = elasticSettlement({ q: Vk / A, B: inp.Db, E, nu, shape: 'circle', point: 'rigid' }).s;
    const tilt = footingRotation({ M: Mk, B: inp.Db, E, nu, shape: 'circle' }).theta / DEG;
    // Clay only: primary consolidation under the sustained weight (Terzaghi one-dimensional theory, two-way drainage)
    // and the undrained strength implied by Modified Cam-Clay for the stated stress history.
    let sCons = 0, t90 = 0, u5 = 0, suMcc = 0, suIn = 0;
    if (clay) {
      const qSus = (inp.Vtop * 1e6 + Wc + Wb) / A;
      sCons = run(() => consolidationSettlement({ q: qSus, B: inp.Db, shape: 'circle', sigmaV0Base: gam * inp.skirt, layers: [{ thickness: inp.Hclay, gamma: gam, Cc: inp.Cc, Cr: inp.Cc / 6, e0: inp.e0, ocr: inp.ocr }] }), 'Consolidation settlement').s;
      const Hd = inp.Hclay / 2;
      t90 = (0.848 * Hd * Hd) / inp.cv; u5 = consolidationU((inp.cv * 5) / (Hd * Hd));
      const Mcs = (6 * Math.sin(inp.phiCs * DEG)) / (3 - Math.sin(inp.phiCs * DEG));
      const zRef = Math.max(inp.skirt + 0.25 * inp.Db, 1);
      // mean effective stress with K0 = (1 − sinφ)·OCR^sinφ (Mayne & Kulhawy 1982)
      const K0 = (1 - Math.sin(inp.phiCs * DEG)) * inp.ocr ** Math.sin(inp.phiCs * DEG);
      suMcc = mccUndrainedStrengthRatio({ M: Mcs, lambda: 1, kappa: 0.2, ocr: inp.ocr }) * gam * zRef * (1 + 2 * K0) / 3;
      suIn = (inp.su0 + inp.suGrad * zRef) * 1e3;
    }
    const kern = ot.kernLimit;
    const checks = [
      check('Overturning: factored eccentricity against the base radius limit', ed, 0.6 * R, 'm', 'DNV-ST-0126 (2021) App. G; ISO 19901-4:2016 §7.2 (effective area, e ≤ 0.3 D)', 'Guards against the load resultant leaving the base. Beyond 0.3 diameters the compressed zone becomes a thin sliver, contact pressure rises sharply and the base rocks or overturns.'),
      check('Sliding at base / skirt-tip level', Hd / 1e6, Math.max(Hres, 1e-9) / 1e6, 'MN', 'DNV-ST-0126 (2021) §7.5.3; ISO 19901-4:2016 §7.3', 'Guards against the foundation being pushed sideways by wave and wind load. The resistance is friction or undrained strength on the effective area plus passive pressure on the skirts.'),
      check('Bearing capacity on the effective area', Vd / 1e6, Math.max(Qd, 1e-9) / 1e6, 'MN', clay ? 'ISO 19901-4:2016 §7.2.2; API RP 2GEO §7.3.1 (undrained)' : 'API RP 2A-WSD §6.13.2; ISO 19901-4:2016 §7.2.3 (drained)', 'Guards against a shear failure of the soil under the loaded edge of the base. Load inclination and eccentricity cut the capacity far below the central vertical value.'),
      check('No uplift under characteristic loads (resultant within the kern)', ot.e, kern, 'm', 'DNV-ST-0126 (2021) §7.5.5 (serviceability: full base contact)', 'Keeps the whole base in compression in the characteristic extreme. Repeated gapping pumps water and soil from under the base and softens the foundation.'),
      check('Permanent tilt from the extreme moment', tilt, inp.tiltLim, 'deg', 'DNV-ST-0126 (2021) §7.5.6; turbine supplier tilt tolerance', 'Elastic rotation of the base under the characteristic moment, an indicator of the tilt accumulated over storms. Beyond the tolerance the turbine cannot operate.'),
    ];
    if (overturned) { checks[1].uc = 99; checks[2].uc = 99; }
    const gov = worst(checks), ok = gov.uc <= 1;
    // Ballast needed to bring the governing stability utilisation to 1 (bisection on ballast mass; pure statics for eccentricity).
    const ballastFor = (limit) => { const Vneed = Md / limit; return Math.max((Vneed + inp.gammaF * inp.uplift * 1e6) / inp.gammaG - inp.Vtop * 1e6 - Wc, 0) / ((1 - inp.rhoW / inp.rhoBallast) * G) / 1e3; };
    const bEcc = ballastFor(0.6 * R), bKern = Math.max((Mk / kern + inp.uplift * 1e6 - inp.Vtop * 1e6 - Wc) / ((1 - inp.rhoW / inp.rhoBallast) * G) / 1e3, 0);
    const ds = [0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3].map((f) => f * inp.Db);
    const eOverR = ds.map((Dd) => { const Wb2 = Wb * (Dd / inp.Db) ** 2; const V2 = inp.gammaG * (inp.Vtop * 1e6 + Wc + Wb2) - inp.gammaF * inp.uplift * 1e6; return V2 > 0 ? Math.min(Md / V2 / (0.6 * Dd / 2), 5) : 5; });
    // Interaction: allowable H against V for sliding and bearing at the present eccentricity (display only).
    return {
      kpis: [
        { key: 'fsOT', label: 'Overturning safety factor (characteristic)', unit: '-', value: cap(ot.fsOverturning), explain: 'Stabilising moment V·R over the overturning moment at skirt-tip level, without factors. Values below about 2 normally fail the eccentricity limit once factors are applied.' },
        { key: 'ecc', label: 'Factored load eccentricity / radius', unit: '-', value: ed / R, limit: 0.6, sense: 'max', explain: 'Position of the factored load resultant relative to the base centre. 0.25 keeps full contact; the bearing check becomes very sensitive beyond about 0.5.' },
        { key: 'ucSlide', label: 'Sliding utilisation', unit: '-', value: checks[1].uc, limit: 1, sense: 'max', explain: 'Factored horizontal load over factored sliding resistance on the effective area including skirt passive resistance.' },
        { key: 'ucBear', label: 'Bearing utilisation', unit: '-', value: checks[2].uc, limit: 1, sense: 'max', explain: 'Factored vertical load over design bearing capacity of the effective area under the inclined, eccentric load.' },
        { key: 'QbearChar', label: 'Design bearing capacity', unit: 'MN', value: Qd / 1e6, explain: 'Vertical load the effective area can carry with the material factor applied, for the present load inclination and eccentricity.' },
        { key: 'Aeff', label: 'Effective bearing area', unit: 'm²', value: Aeff, explain: `Part of the base assumed to carry the load (Meyerhof), ${fmt(Aeff / A * 100)} % of the full base. It shrinks rapidly as the moment grows.` },
        { key: 'qMax', label: 'Design contact pressure on the effective area', unit: 'kPa', value: Number.isFinite(qMax) ? qMax / 1e3 : 0, explain: 'Factored vertical load divided by the effective area: the pressure the soil and the base slab must sustain at the loaded edge.' },
        { key: 'Wsub', label: 'Submerged weight of base and ballast', unit: 'MN', value: (Wc + Wb) / 1e6, explain: 'Stabilising weight after buoyancy. It is the cheapest source of stability: ballast adds weight without structural concrete.' },
        { key: 'settle', label: 'Immediate settlement', unit: 'mm', value: settle * 1e3, explain: 'Elastic settlement of the rigid base under the characteristic vertical load. Consolidation settlement of clay comes on top and is not computed here.' },
        { key: 'sCons', label: 'Consolidation settlement (clay)', unit: 'mm', value: sCons * 1e3, explain: 'Long-term primary consolidation of the clay layer under the sustained weight (zero for sand). It develops over years and must be allowed for in the interface level and the cable hang-off.' },
        { key: 't90', label: 'Time to 90 % consolidation (clay)', unit: 'years', value: t90, explain: 'Time for nine tenths of the consolidation settlement with drainage to both faces of the clay layer (zero for sand). Strength gain from consolidation follows the same clock.' },
        { key: 'tilt', label: 'Rotation under the characteristic moment', unit: 'deg', value: tilt, limit: inp.tiltLim, sense: 'max', explain: 'Elastic rotation of the base in the extreme event; repeated loading leaves a fraction of it as permanent tilt.' },
        { key: 'concrete', label: 'Concrete volume', unit: 'm³', value: inp.Vconc, explain: 'Structural concrete of the foundation, the main embodied-carbon item of a gravity base.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Utilisation by failure mode', categories: ['Overturning (eccentricity)', 'Sliding', 'Bearing', 'Uplift (kern)', 'Tilt'], series: [{ name: 'Utilisation', y: checks.map((c) => Math.min(c.uc, 3)) }], x: { label: 'Failure mode' }, y: { label: 'Utilisation', unit: '-' },
          explain: 'Demand over capacity for each stability mode (capped at 3 for display). The tallest bar governs the base diameter and ballast.' },
        { type: 'line', title: 'Eccentricity utilisation against base diameter', x: { label: 'Base diameter', unit: 'm' }, y: { label: 'Factored eccentricity / limit', unit: '-' },
          series: [{ name: 'Ballast scaled with base area', x: ds, y: eOverR }, { name: 'Limit', x: [ds[0], ds[ds.length - 1]], y: [1, 1], dash: true }, { name: 'Present design', x: [inp.Db], y: [ed / (0.6 * R)] }],
          explain: 'How the overturning utilisation falls as the base widens with ballast in proportion to its area. The crossing with the limit is the smallest stable diameter.' },
      ],
      tables: [
        { title: 'Loads and resistances', columns: [{ key: 'k', label: 'Quantity' }, { key: 'c', label: 'Characteristic' }, { key: 'd', label: 'Design' }, { key: 'u', label: 'Unit' }],
          rows: [['Vertical load at base', Vk / 1e6, Vd / 1e6, 'MN'], ['Horizontal load', Hk / 1e6, Hd / 1e6, 'MN'], ['Moment at skirt-tip level', Mk / 1e6, Md / 1e6, 'MN·m'], ['Eccentricity', ot.e, ed, 'm'], ['Load inclination H/V', Hk / Vk, Hd / Vd, '-'],
            ['Sliding resistance', '–', Hres / 1e6, 'MN'], ['Bearing capacity, code method', '–', Qd / 1e6, 'MN'], ['Bearing capacity, Meyerhof factors', '–', meyer / 1e6, 'MN'], ['Terzaghi capacity, central vertical load (unfactored strength)', terz / 1e6, '–', 'MN'],
            ['Effective width × length', '–', `${fmt(Beff)} × ${fmt(Leff)}`, 'm'], ['Ballast for the eccentricity limit', '–', bEcc, 't'], ['Ballast for full contact (kern)', bKern, '–', 't'],
            ...(clay ? [['Degree of consolidation after 5 years', u5 * 100, '–', '%'], ['Undrained strength at 0.25 D below the base, entered', suIn / 1e3, '–', 'kPa'], ['Undrained strength there from Modified Cam-Clay (κ/λ = 0.2)', suMcc / 1e3, '–', 'kPa']] : [])].map(([k, c, d2, u]) => ({ k, c: typeof c === 'number' ? round(c) : c, d: typeof d2 === 'number' ? round(d2) : d2, u })),
          explain: 'Characteristic and factored actions at the base with the resistances of each method. The Terzaghi value for a central vertical load shows how much capacity is lost to inclination and eccentricity.' },
      ],
      notes: [`${gov.name} governs (utilisation ${fmt(gov.uc)}).`, `${clay ? `Undrained clay, su = ${fmt(suBase / 1e3)} kPa at skirt-tip level` : `Drained sand, φ′ = ${inp.phiDeg}° (design ${fmt(phi / DEG)}°)`}; effective area ${fmt(Aeff)} m² of ${fmt(A)} m².`,
        overturned ? 'The factored resultant falls outside the base: sliding and bearing cannot be evaluated until weight or diameter is increased.' : `Load inclination H/V = ${fmt(Hd / Vd)} (design).`, 'Scour around the base edge reduces the effective area: see the scour module for protection sizing.'],
      advice: [
        ok ? adv('engineering', sev(gov.uc), `The base is stable: eccentricity ${fmt(ed / R)} R, sliding ${fmt(checks[1].uc * 100)} %, bearing ${fmt(checks[2].uc * 100)} %, tilt ${fmt(tilt)}°.`,
          `Confirm the cyclic soil strength and include the wave load on the base in the applied loads. Ballast could be reduced towards ${fmt(Math.max(bEcc, bKern))} t (from ${fmt(inp.ballastT)} t) before overturning or uplift governs; re-check sliding and bearing at that weight.`,
          'DNV-ST-0126 (2021) §7.5 and App. G; ISO 19901-4:2016 §7',
          `${concreteCo2Text(inp.Vconc)} Stability bought with ballast instead of concrete avoids most of that: ballast sand has a small fraction of the footprint per tonne.`)
          : adv('engineering', 'fail', `${gov.name}: ${fmt(gov.demand)} against ${fmt(gov.capacity)} ${gov.unit} (utilisation ${fmt(gov.uc)}).`,
            gov.name.includes('Sliding') ? `Add skirts (now ${inp.skirt} m) or ballast: sliding resistance grows with weight in sand and with base area and skirt depth in clay.` : gov.name.includes('Tilt') ? 'Widen the base (rotation ∝ 1/D³) or improve the ground (gravel bed, replacement of soft surface soil).' : `Increase ballast to about ${fmt(Math.max(bEcc, bKern))} t or widen the base; eccentricity falls in proportion to the vertical load.`,
            'DNV-ST-0126 (2021) §7.5; ISO 19901-4:2016 §7.2–7.3',
            `Reaching stability with ${fmt(Math.max(Math.max(bEcc, bKern) - inp.ballastT, 0))} t more ballast needs no extra concrete; widening the slab by 10 % instead adds roughly ${fmt(0.2 * inp.Vconc * 0.5)} m³ — ${concreteCo2Text(0.2 * inp.Vconc * 0.5)}`),
        adv('economic', 'info', `The foundation uses ${fmt(inp.Vconc)} m³ of concrete and ${fmt(inp.ballastT)} t of ballast for ${fmt((Wc + Wb) / 1e6)} MN of submerged weight.`,
          'Gravity bases trade cheap material against heavy logistics: check quay bearing capacity, float-out or heavy-lift method and seabed preparation (dredging and gravel bed) before comparing with a monopile on steel tonnage alone.',
          'AACE RP 18R-97 Class 4 quantities; installation method statement and marine-spread day rates as the cost drivers of gravity foundations'),
      ],
      quantities: { concrete_m3: round(inp.Vconc) },
    };
  },
  verify: [
    { name: 'Overturning safety factor of a circular base: FS = V·R / M', source: 'Statics: concrete 1000 m³ at 2500 − 1025 kg/m³ → 14.465 MN submerged, plus 5.535 MN from the turbine = 20 MN; R = 15 m, M = 100 MN·m → FS = 20 × 15 / 100 = 3.0',
      inputs: { Db: 30, skirt: 0, concrete: 'c45', Vconc: 1000, ballastT: 0, Vtop: 5.535, uplift: 0, H: 0, M: 100, rhoW: 1025 }, expect: [{ kpi: 'fsOT', value: 3.0, tol: 0.002 }] },
    { name: 'Undrained bearing capacity of a circular footing under central vertical load: Q = 5.14 × 1.18 × su × A', source: 'ISO 19901-4:2016 §7.2.2 with shape factor sc = 0.18 for uniform strength (exact solution 6.05·su, Eason & Shield 1960): 5.1416 × 1.18 × 50 kPa × π × 15² m² = 214.4 MN',
      inputs: { soil: 'clay', su0: 50, suGrad: 0, Db: 30, skirt: 0, H: 0, M: 0.0001, gammaM: 1, uplift: 0 }, expect: [{ kpi: 'QbearChar', value: 214.4, tol: 0.003 }] },
  ],
};
