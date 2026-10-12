// Foundation–soil interaction: monopile as a beam on nonlinear p–y springs (deflection, rotation, bending moment,
// soil mobilisation, pile-head stiffness matrix, embedded-length criteria, lateral capacity, accumulated rotation)
// or a gravity-base foundation (bearing capacity on the effective area, sliding, overturning and uplift).
import { solveLateralPile, ultimateLateralCapacity } from '../../engine/geotech/pile.js';
import { monopileStiffness, mudlineRotationCheck, cyclicRotationAccumulation, embeddedLengthCriteria } from '../../engine/geotech/monopile.js';
import { createProfile, applyScour, validateLayers } from '../../engine/geotech/soilprofile.js';
import { effectiveAreaCircle, bearingCapacity, slidingResistance, overturningCheck } from '../../engine/geotech/bearing.js';
import { tubularSectionStress } from '../../engine/aero/blade.js';
import { adv, check, worst, sev, fmt, round, run, need, thin, DEG, LIMIT_CODE, CONFIRM, steelCo2Text, CO2_STEEL, CO2_CONCRETE } from './_common.js';

const mono = { key: 'concept', in: ['monopile'] }, gbs = { key: 'concept', in: ['gravity'] };
const LAYERS = [{ type: 'sand', thickness: 8, gammaSub: 9.5, su: 0, suGrad: 0, phi: 33, eps50: 0 }, { type: 'clay', thickness: 6, gammaSub: 9, su: 90, suGrad: 3, phi: 0, eps50: 0.007 }, { type: 'sand', thickness: 60, gammaSub: 10, su: 0, suGrad: 0, phi: 38, eps50: 0 }];

function profileOf(inp, depthNeeded) {
  let layers;
  if (inp.soilSpec === 'sand') layers = [{ type: 'sand', top: 0, thickness: depthNeeded + 20, gammaSub: inp.gammaSub * 1e3, phiDeg: inp.phi }];
  else if (inp.soilSpec === 'clay') layers = [{ type: 'clay', top: 0, thickness: depthNeeded + 20, gammaSub: inp.gammaSub * 1e3, suTop: inp.su0 * 1e3, suGrad: inp.suGrad * 1e3, eps50: inp.eps50 }];
  else {
    const rows = (inp.layers || []).filter((r) => r.thickness > 0);
    need(rows.length >= 1, 'The soil-layer table needs at least one layer with a positive thickness.');
    let top = 0;
    layers = rows.map((r, i) => {
      const type = String(r.type || '').trim().toLowerCase();
      need(type === 'sand' || type === 'clay', `Soil layer ${i + 1}: type must be 'sand' or 'clay'.`);
      need(r.gammaSub > 2 && r.gammaSub < 15, `Soil layer ${i + 1}: submerged unit weight must be between 2 and 15 kN/m³.`);
      if (type === 'sand') need(r.phi >= 29 && r.phi <= 45, `Soil layer ${i + 1}: friction angle of sand must be between 29° and 45° (range of the API initial-modulus correlation).`); else need(r.su > 0, `Soil layer ${i + 1}: clay needs a positive undrained shear strength.`);
      const L = type === 'sand' ? { type, top, thickness: r.thickness, gammaSub: r.gammaSub * 1e3, phiDeg: r.phi } : { type, top, thickness: r.thickness, gammaSub: r.gammaSub * 1e3, suTop: r.su * 1e3, suGrad: (r.suGrad || 0) * 1e3, eps50: r.eps50 > 0 ? r.eps50 : 0.01 };
      top += r.thickness; return L;
    });
    need(top >= depthNeeded + 1, `The soil layers reach ${fmt(top)} m but the foundation needs ${fmt(depthNeeded + 1)} m: extend the last layer.`);
  }
  const v = validateLayers(layers);
  need(v.valid, `Soil profile: ${v.errors.join('; ')}`);
  return createProfile(layers);
}

export default {
  id: 's16-foundation-soil', suite: 16, title: 'Foundation–soil interaction (monopile p–y or gravity base)',
  summary: 'Lateral response, stiffness and capacity of a monopile in layered soil, or bearing, sliding and overturning of a gravity-base foundation, under the turbine mudline loads.',
  covers: [9], fidelity: 'numerical',
  standards: ['DNV-ST-0126 (2021) §7 (geotechnical design of monopiles and gravity-based foundations)', 'API RP 2GEO (2011, add. 2014) §8.5 (p–y curves) and §7 (shallow foundations)', 'ISO 19901-4:2016 (geotechnical and foundation design)', 'IEC 61400-3-1:2019 §7.6', 'LeBlanc, Houlsby & Byrne (2010) Géotechnique 60(2) (accumulated rotation)'],
  models: ['Nonlinear pile–soil resistance: API p–y curves for sand (O\'Neill & Murchison) and soft clay (Matlock), static or cyclic', 'Beam on nonlinear Winkler foundation (Hermite finite elements with P–Δ)', 'Pile-head stiffness matrix by static condensation (three-spring model)', 'Hetényi beam on elastic foundation (linear-spring option)', 'Embedded-length criteria (zero toe-kick, vertical tangent, critical length)', 'Foundation bearing capacity on the effective area (Vesic factors)', 'Sliding resistance and overturning statics of a gravity base', 'Scour as loss of soil support and overburden'],
  limits: `API p–y curves were derived for slender piles and are known to misrepresent the stiffness of large-diameter, rigid monopiles (PISA-type site-calibrated soil reactions including base shear and distributed moment are the current alternative): treat stiffness and rotation as preliminary. The pile is linear elastic; axial capacity, driveability, liquefaction and cyclic degradation of clay are not covered here. Gravity-base checks assume homogeneous soil under the base. ${LIMIT_CODE}`,
  publishes: ['foundation.Kh', 'foundation.Kr', 'foundation.Khr', 'foundation.headRotationDeg', 'quantities.steel_t', 'quantities.concrete_m3'],
  searchTerms: 'monopile p-y curves lateral stiffness accumulated rotation gravity base bearing capacity offshore wind',
  inputs: [
    { key: 'concept', label: 'Foundation concept', type: 'select', default: 'monopile', group: 'Foundation', role: 'geometry', options: [{ value: 'monopile', label: 'Monopile (laterally loaded pile)' }, { value: 'gravity', label: 'Gravity-base foundation' }] },
    { key: 'D', label: 'Pile outer diameter', unit: 'm', type: 'number', default: 6.5, min: 1, max: 14, group: 'Foundation', role: 'geometry', showIf: mono, design: true },
    { key: 't', label: 'Pile wall thickness', unit: 'mm', type: 'number', default: 80, min: 15, max: 160, group: 'Foundation', role: 'geometry', showIf: mono, design: true },
    { key: 'Lp', label: 'Embedded length', unit: 'm', type: 'number', default: 38, min: 5, max: 100, group: 'Foundation', role: 'geometry', showIf: mono, design: true, help: 'Penetration below the original mudline; typically 4–6 diameters for monopiles.' },
    { key: 'steel', label: 'Pile steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Foundation', role: 'material', showIf: mono },
    { key: 'gbsD', label: 'Base diameter', unit: 'm', type: 'number', default: 34, min: 8, max: 60, group: 'Foundation', role: 'geometry', showIf: gbs, design: true },
    { key: 'gbsW', label: 'Submerged weight of base and ballast', unit: 'MN', type: 'number', default: 85, min: 1, max: 600, group: 'Foundation', role: 'load', showIf: gbs, design: true, help: 'Concrete structure plus sand or rock ballast, buoyancy deducted.' },
    { key: 'gbsEmb', label: 'Base embedment (skirt or excavation depth)', unit: 'm', type: 'number', default: 1, min: 0, max: 10, group: 'Foundation', role: 'geometry', showIf: gbs },
    { key: 'gbsConc', label: 'Concrete volume of the base', unit: 'm³', type: 'number', default: 2500, min: 0, max: 30000, group: 'Foundation', role: 'geometry', showIf: gbs, help: 'For the quantity and carbon take-off only.' },
    { key: 'soilSpec', label: 'Soil description', type: 'select', default: 'sand', group: 'Soil', role: 'soil', options: [{ value: 'sand', label: 'Uniform sand' }, { value: 'clay', label: 'Uniform clay' }, { value: 'layers', label: 'Layer table' }, { value: 'linear', label: 'Linear springs (constant modulus of subgrade reaction)' }],
      help: 'Pile–soil boundary condition. The layer table and linear springs apply to the monopile; the gravity base uses the uniform sand or clay data.' },
    { key: 'phi', label: 'Friction angle of sand', unit: 'deg', type: 'number', default: 36, min: 29, max: 45, group: 'Soil', role: 'soil', bind: 'site.phiDeg', showIf: { key: 'soilSpec', in: ['sand'] }, uncertain: { dist: 'normal', cov: 0.05 } },
    { key: 'su0', label: 'Undrained shear strength at the mudline', unit: 'kPa', type: 'number', default: 60, min: 2, max: 600, group: 'Soil', role: 'soil', bind: 'site.su0', showIf: { key: 'soilSpec', in: ['clay'] }, uncertain: { dist: 'lognormal', cov: 0.25 } },
    { key: 'suGrad', label: 'Strength gradient with depth', unit: 'kPa/m', type: 'number', default: 2, min: 0, max: 20, group: 'Soil', role: 'soil', bind: 'site.suGrad', showIf: { key: 'soilSpec', in: ['clay'] } },
    { key: 'eps50', label: 'Strain at half the maximum stress ε50', unit: '-', type: 'number', default: 0.01, min: 0.002, max: 0.03, group: 'Soil', role: 'soil', showIf: { key: 'soilSpec', in: ['clay'] }, help: 'Matlock: 0.02 soft, 0.01 firm, 0.005 stiff clay.' },
    { key: 'gammaSub', label: 'Submerged unit weight', unit: 'kN/m³', type: 'number', default: 10, min: 4, max: 13, group: 'Soil', role: 'soil', bind: 'site.gammaSub', showIf: { key: 'soilSpec', in: ['sand', 'clay'] } },
    { key: 'kLin', label: 'Modulus of subgrade reaction', unit: 'MN/m²', type: 'number', default: 50, min: 0.5, max: 2000, group: 'Soil', role: 'soil', showIf: { key: 'soilSpec', in: ['linear'] }, help: 'Force per metre of pile per metre of deflection, constant with depth (Hetényi model).' },
    { key: 'layers', label: 'Soil layers (from the mudline down)', type: 'table', group: 'Soil', role: 'soil', showIf: { key: 'soilSpec', in: ['layers'] },
      columns: [{ key: 'type', label: 'Type (sand / clay)', type: 'text' }, { key: 'thickness', label: 'Thickness', unit: 'm', type: 'number' }, { key: 'gammaSub', label: "γ'", unit: 'kN/m³', type: 'number' }, { key: 'su', label: 'su at top', unit: 'kPa', type: 'number' }, { key: 'suGrad', label: 'su gradient', unit: 'kPa/m', type: 'number' }, { key: 'phi', label: 'φ', unit: 'deg', type: 'number' }, { key: 'eps50', label: 'ε50', unit: '-', type: 'number' }],
      default: LAYERS.map((r) => ({ ...r })), help: 'Strength columns that do not apply to a layer type are ignored. The last layer must extend below the pile toe.' },
    { key: 'scour', label: 'Local scour depth', unit: 'm', type: 'number', default: 0, min: 0, max: 20, group: 'Soil', role: 'ic', bind: 'scour.depth', showIf: mono, help: 'Initial seabed state: soil support removed around the pile, with reduced overburden to six diameters below (API RP 2A C6.3).' },
    { key: 'loading', label: 'p–y loading type', type: 'select', default: 'cyclic', group: 'Soil', role: 'soil', showIf: mono, options: [{ value: 'cyclic', label: 'Cyclic (storm) curves' }, { value: 'static', label: 'Static curves' }] },
    { key: 'H', label: 'Horizontal load at the mudline (characteristic)', unit: 'MN', type: 'number', default: 5.8, min: 0.01, max: 100, group: 'Loads at the mudline', role: 'bc', bind: 'owt.mudlineShear', uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'From the extreme-loads module (governing characteristic case).' },
    { key: 'M', label: 'Overturning moment at the mudline (characteristic)', unit: 'MN·m', type: 'number', default: 325, min: 0, max: 5000, group: 'Loads at the mudline', role: 'bc', bind: 'owt.mudlineMoment', uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'V', label: 'Vertical load from turbine, tower and substructure', unit: 'MN', type: 'number', default: 10, min: 0, max: 300, group: 'Loads at the mudline', role: 'bc', bind: 'owt.mudlineAxial' },
    { key: 'slsFrac', label: 'Serviceability load / characteristic extreme load', unit: '-', type: 'number', default: 0.5, min: 0.1, max: 1, group: 'Loads at the mudline', role: 'load', help: 'Load level at which the stiffness is linearised and cyclic rotation accumulates (operating loads are about half the extreme).' },
    { key: 'nCycles', label: 'Equivalent load cycles at that level', unit: '-', type: 'number', default: 1e5, min: 1, max: 1e9, group: 'Loads at the mudline', role: 'load', showIf: mono },
    { key: 'Tb', label: 'LeBlanc load-magnitude factor Tb', unit: '-', type: 'number', default: 0.05, min: 0, max: 0.5, group: 'Loads at the mudline', role: 'soil', showIf: mono, uncertain: { dist: 'lognormal', cov: 0.4 }, help: 'Chart value (LeBlanc et al. 2010): 0.03–0.15, growing with load level and falling with density.' },
    { key: 'gammaF', label: 'Load factor', unit: '-', type: 'number', default: 1.35, min: 1, max: 2, group: 'Acceptance', role: 'criteria', help: `IEC 61400-3-1 normal load cases. ${CONFIRM}` },
    { key: 'gammaM', label: 'Material factor (steel yield; soil resistance)', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.6, group: 'Acceptance', role: 'criteria', help: `DNV-ST-0126: 1.10 steel, 1.15 on tan φ (drained) and 1.25 on su (undrained). One value is applied here. ${CONFIRM}` },
    { key: 'rotLimit', label: 'Total mudline rotation limit', unit: 'deg', type: 'number', default: 0.5, min: 0.05, max: 2, group: 'Acceptance', role: 'criteria', showIf: mono, help: 'Commonly 0.5° including 0.25° installation tolerance; set by the turbine supplier.' },
    { key: 'nElem', label: 'Pile elements', type: 'int', default: 40, min: 10, max: 320, group: 'Solver', role: 'control', showIf: mono },
  ],
  mesh: { key: 'nElem', label: 'Pile elements', levels: [20, 40, 80, 160], qoi: 'yHead' },
  compute(inp, ctx) {
    const H = inp.H * 1e6, M = inp.M * 1e6, V = inp.V * 1e6;
    if (inp.concept === 'gravity') return gravityBase(inp, { H, M, V });
    const st = ctx.material(inp.steel);
    need(st && st.fy > 0, `Steel '${inp.steel}' is not in the materials library.`);
    const D = inp.D, t = inp.t / 1e3;
    need(t < D / 4, 'Wall thickness must be below a quarter of the diameter.');
    need(inp.Lp > inp.scour + 2 * D, `The scour depth of ${inp.scour} m leaves less than two diameters of embedment: lengthen the pile or protect against scour.`);
    const linear = inp.soilSpec === 'linear';
    let profile = null;
    if (!linear) { profile = profileOf(inp, inp.Lp); if (inp.scour > 0) profile = applyScour(profile, { general: 0, local: inp.scour, D }); }
    const base = { L: inp.Lp, D, t, E: st.E, nElem: inp.nElem, axialLoad: V, ...(linear ? { pyAt: () => inp.kLin * 1e6 } : { profile, loading: inp.loading }) };
    const r = run(() => solveLateralPile({ ...base, head: { H, M } }), 'Lateral pile analysis at the characteristic load (the soil cannot carry this load if it does not converge — lengthen or enlarge the pile)');
    const sl = { H: H * inp.slsFrac, M: M * inp.slsFrac };
    const sls = run(() => solveLateralPile({ ...base, head: sl }), 'Lateral pile analysis at the serviceability load');
    const K = run(() => monopileStiffness({ ...base, head: sl }), 'Pile-head stiffness').secant;
    const acc = cyclicRotationAccumulation({ thetaStatic: Math.abs(sls.head.theta), N: inp.nCycles, Tb: inp.Tb, Tc: 1 });
    const rot = mudlineRotationCheck({ thetaTotal: Math.abs(sls.head.theta) + acc.deltaTheta, thetaAccumulated: acc.deltaTheta, limitTotalDeg: inp.rotLimit, limitAccumulatedDeg: inp.rotLimit / 2 });
    // lateral capacity with the load eccentricity e = M/H (not defined for linear springs)
    let ult = null;
    if (!linear) { try { ult = ultimateLateralCapacity({ ...base, nElem: Math.min(inp.nElem, 40) }, { eccentricity: H > 0 ? M / H : 0, yLimit: 0.1 * D, nSteps: 24 }); } catch { ult = null; } }
    // embedded-length criteria on a coarse scan
    let crit = null;
    try { const { L, ...noL } = base; crit = embeddedLengthCriteria({ ...noL, nElem: Math.min(inp.nElem, 30), head: { H, M } }, { Lmin: Math.max(inp.scour + 2.5 * D, 0.5 * inp.Lp), Lmax: 1.5 * inp.Lp, nLengths: 9 }); } catch { crit = null; }
    const sec = tubularSectionStress({ D, t, N: V * inp.gammaF, M: r.maxMoment * inp.gammaF, V: 0 });
    const fyd = st.fy / 1.1;
    const checks = [
      check('Mudline rotation at the serviceability load, including accumulation', rot.totalDeg, inp.rotLimit, 'deg', 'DNV-ST-0126 (2021) §7.6 (permanent rotation); LeBlanc et al. (2010)', 'Rotation of the pile at the mudline under the repeated service load plus the rotation accumulated over the load cycles. Beyond the limit the turbine tilts permanently: yaw bearing and gearbox alignment suffer and the supplier warranty is at stake.'),
      check('Pile steel stress at the section of maximum moment', sec.vonMises / 1e6, fyd / 1e6, 'MPa', 'DNV-ST-0126 (2021) §4.4 (ULS, γM = 1.10)', `Factored bending and axial stress at the point of maximum moment, ${fmt(r.depthMaxMoment)} m below the mudline, against design yield. The moment peaks below the mudline because the upper soil is weak.`),
      ...(ult ? [check('Lateral soil capacity (pile-head deflection limit 0.1 D)', H * inp.gammaF / 1e6, ult.Hult / inp.gammaM / 1e6, 'MN', 'DNV-ST-0126 (2021) §7.6 (ULS lateral capacity); API RP 2GEO §8.5', 'Factored horizontal load (with its moment) against the load at which the soil reaction is exhausted or the head deflects by a tenth of the diameter, divided by the material factor. Exceeding it means a soil failure mechanism in the design storm.')] : []),
      ...(crit && crit.critical ? [check('Embedded length against the critical pile length', crit.critical, inp.Lp, 'm', 'Monopile design practice (critical length; DNV-ST-0126 §7.6 guidance)', 'Demand is the length beyond which more penetration no longer reduces head deflection, capacity the actual embedment. A shorter pile behaves as a rigid stub whose rotation is sensitive to soil variability.')] : []),
    ];
    const gov = worst(checks);
    const steelT = Math.PI * (D - t) * t * inp.Lp * (st.rho ?? 7850) / 1e3;
    const [z, y, mom, p, mob] = thin([r.z, r.y, r.moment, r.p, r.mobilisation], 160);
    const zs = z.map((v) => round(-v, 3));
    const dK = (a, b) => (b - a) / a * 100;
    return {
      kpis: [
        { key: 'yHead', label: 'Pile-head deflection at the characteristic load', value: round(r.head.y * 1e3, 5), unit: 'mm', explain: 'Lateral displacement of the pile at the mudline under the unfactored extreme load. A few tens of millimetres is normal for a monopile.' },
        { key: 'rotHead', label: 'Pile-head rotation at the characteristic load', value: round(Math.abs(r.head.theta) / DEG, 5), unit: 'deg', explain: 'Rotation at the mudline under the extreme load (elastic part; it recovers largely after the storm).' },
        { key: 'rotSls', label: 'Rotation at the serviceability load plus accumulation', value: round(rot.totalDeg, 5), unit: 'deg', limit: inp.rotLimit, sense: 'max', explain: `Rotation under the repeated service load (${fmt(Math.abs(sls.head.theta) / DEG)}°) plus the rotation accumulated over ${fmt(inp.nCycles)} cycles (${fmt(rot.accumulatedDeg)}°).` },
        { key: 'Mmax', label: 'Maximum pile bending moment (characteristic)', value: round(r.maxMoment / 1e6, 5), unit: 'MN·m', explain: `Largest moment in the pile, ${fmt(r.depthMaxMoment)} m below the mudline. It exceeds the mudline moment because the shear keeps adding moment until the soil reaction has balanced it.` },
        { key: 'Kh', label: 'Lateral stiffness K_LL (secant, service load)', value: round(K.KLL / 1e9, 5), unit: 'GN/m', explain: 'Force per unit head displacement with zero rotation. With K_LR and K_RR it is the foundation boundary condition of the natural-frequency and dynamic modules.' },
        { key: 'Kr', label: 'Rocking stiffness K_RR', value: round(K.KRR / 1e9, 5), unit: 'GN·m/rad', explain: 'Moment per unit head rotation with zero displacement. It controls the first natural frequency of the turbine more than the lateral term.' },
        { key: 'Khr', label: 'Coupling stiffness K_LR', value: round(K.KLR / 1e6, 5), unit: 'MN', explain: 'Off-diagonal term (negative): a horizontal head force also rotates the pile head.' },
        { key: 'mobilisation', label: 'Largest soil mobilisation p/pu', value: round(r.maxMobilisation, 4), unit: '-', explain: 'Highest ratio of soil reaction to ultimate resistance along the pile (0 for linear springs). Values near 1 over the upper diameters are usual; near 1 over the whole length means little reserve.' },
        { key: 'Hult', label: 'Lateral capacity at this load eccentricity', value: round(ult ? ult.Hult / 1e6 : 0, 5), unit: 'MN', explain: `Horizontal load at which the soil fails or the head deflects by 0.1 D (${ult ? ult.criterion : 'not evaluated for linear springs'}).` },
        { key: 'Lcrit', label: 'Critical embedded length', value: round(crit && crit.critical ? crit.critical : 0, 4), unit: 'm', explain: 'Shortest length whose head response is within 2 % of a much longer pile (0 if not reached in the scanned range). Longer piles add steel without stiffness.' },
        { key: 'ucSteel', label: 'Steel utilisation below the mudline', value: round(checks[1].uc, 4), unit: '-', limit: 1, sense: 'max', explain: 'Factored stress over design yield at the section of maximum moment.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Pile deflection', x: { label: 'Lateral deflection', unit: 'mm' }, y: { label: 'Elevation relative to the mudline', unit: 'm' }, series: [{ name: 'Characteristic load', x: y.map((v) => round(v * 1e3, 4)), y: zs }, { name: 'Zero', x: [0, 0], y: [0, -inp.Lp], dash: true }],
          explain: 'Deflected shape of the embedded pile. A flexible pile shows a point of zero deflection well above the toe; a toe that kicks back means the pile is short and rotates rigidly.' },
        { type: 'line', title: 'Pile bending moment', x: { label: 'Bending moment', unit: 'MN·m' }, y: { label: 'Elevation relative to the mudline', unit: 'm' }, series: [{ name: 'Characteristic load', x: mom.map((v) => round(v / 1e6, 4)), y: zs }],
          explain: 'Bending moment along the pile: it grows below the mudline to a maximum and then decays as the soil takes the load. The peak sets the governing wall thickness.' },
        { type: 'line', title: 'Soil reaction and mobilisation', x: { label: 'Reaction (MN/m) · mobilisation (-)', unit: '-' }, y: { label: 'Elevation relative to the mudline', unit: 'm' }, series: [{ name: 'Soil reaction p (MN/m)', x: p.map((v) => round(v / 1e6, 5)), y: zs }, { name: 'Mobilisation p/pu', x: mob.map((v) => round(Number.isFinite(v) ? v : 0, 4)), y: zs }],
          explain: 'Soil reaction per metre of pile and the fraction of the ultimate resistance mobilised. The reaction reverses sign below the rotation point; mobilisation of 1 marks soil that has yielded.' },
        ...(ult ? [{ type: 'line', title: 'Load–deflection curve at the pile head', x: { label: 'Head deflection', unit: 'mm' }, y: { label: 'Horizontal load', unit: 'MN' }, series: [{ name: 'H (with M = H·e)', x: ult.curve.y.map((v) => round(v * 1e3, 4)), y: ult.curve.H.map((v) => round(v / 1e6, 5)) }, { name: 'Factored load', x: [0, Math.max(...ult.curve.y) * 1e3], y: [H * inp.gammaF / 1e6, H * inp.gammaF / 1e6], dash: true }],
          explain: 'Horizontal load against head deflection with the moment growing in proportion. The curve softens as the upper soil yields; the factored load should sit on the stiff part.' }] : []),
        ...(crit ? [{ type: 'line', title: 'Head response versus embedded length', x: { label: 'Embedded length', unit: 'm' }, y: { label: 'Head deflection (mm) · rotation (0.01°)', unit: '-' }, series: [{ name: 'Head deflection (mm)', x: crit.L.map((v) => round(v, 3)), y: crit.headY.map((v) => round(v * 1e3, 4)) }, { name: 'Head rotation (0.01°)', x: crit.L.map((v) => round(v, 3)), y: crit.headTheta.map((v) => round(Math.abs(v) / DEG * 100, 4)) }, { name: 'Actual length', x: [inp.Lp, inp.Lp], y: [0, Math.max(...crit.headY) * 1e3], dash: true }],
          explain: 'Head deflection and rotation for a range of pile lengths under the characteristic load. Where the curves flatten, extra length buys nothing: that knee is the critical length.' }] : []),
      ],
      tables: [
        { title: 'Pile-head stiffness matrix (secant at the service load)', columns: [{ key: 'term', label: 'Term' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }, { key: 'note', label: 'Use' }],
          rows: [{ term: 'K_LL', v: round(K.KLL / 1e9, 5), u: 'GN/m', note: 'lateral' }, { term: 'K_LR', v: round(K.KLR / 1e6, 5), u: 'MN', note: 'coupling (negative)' }, { term: 'K_RR', v: round(K.KRR / 1e9, 5), u: 'GN·m/rad', note: 'rocking' }],
          explain: 'Three-spring foundation model [H; M] = K·[u; θ] at the mudline, linearised at the service load. Paste it into the natural-frequency, strength and coupled-dynamics modules (published on the project bus).' },
        { title: 'Pile response summary', columns: [{ key: 'q', label: 'Quantity' }, { key: 'c', label: 'Characteristic load' }, { key: 's', label: 'Service load' }, { key: 'u', label: 'Unit' }],
          rows: [{ q: 'Head deflection', c: round(r.head.y * 1e3, 4), s: round(sls.head.y * 1e3, 4), u: 'mm' }, { q: 'Head rotation', c: round(Math.abs(r.head.theta) / DEG, 5), s: round(Math.abs(sls.head.theta) / DEG, 5), u: 'deg' }, { q: 'Maximum moment', c: round(r.maxMoment / 1e6, 4), s: round(sls.maxMoment / 1e6, 4), u: 'MN·m' }, { q: 'Depth of maximum moment', c: round(r.depthMaxMoment, 3), s: round(sls.depthMaxMoment, 3), u: 'm' }, { q: 'Toe deflection', c: round(r.toe.y * 1e3, 4), s: round(sls.toe.y * 1e3, 4), u: 'mm' }, { q: 'Largest mobilisation', c: round(r.maxMobilisation, 4), s: round(sls.maxMobilisation, 4), u: '-' }],
          explain: 'Head and toe response and peak moment at the extreme and service load levels. The ratio of the two columns shows how nonlinear the soil response is.' },
      ],
      notes: [
        `Monopile ${fmt(D)} m × ${fmt(inp.t)} mm, ${fmt(inp.Lp)} m embedded (L/D = ${fmt(inp.Lp / D)}) in ${inp.soilSpec === 'layers' ? 'layered soil' : inp.soilSpec === 'linear' ? 'linear springs' : `uniform ${inp.soilSpec}`}${inp.scour > 0 ? `, ${inp.scour} m of local scour` : ''}; ${linear ? 'linear' : inp.loading} p–y curves; solver ${r.convergence.converged ? 'converged' : 'did not converge'} in ${r.convergence.iterations} iterations.`,
        `Governing check: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        crit ? `Length criteria: zero toe-kick ${crit.zeroToeKick ? `${fmt(crit.zeroToeKick)} m` : 'not reached'}, vertical tangent ${crit.verticalTangent ? `${fmt(crit.verticalTangent)} m` : 'not reached'}, critical length ${crit.critical ? `${fmt(crit.critical)} m` : 'not reached in the scanned range'}.` : 'Embedded-length scan did not converge for the shortest trial piles.',
        'API p–y curves under-predict the stiffness of large-diameter monopiles at small strain: calibrate against site data before fixing the natural frequency.',
      ],
      advice: [
        adv('engineering', sev(gov.uc), `Service-load rotation with accumulation is ${fmt(rot.totalDeg)}° against ${inp.rotLimit}°; steel utilisation below the mudline ${fmt(checks[1].uc)}${ult ? `, lateral capacity utilisation ${fmt(checks[2].uc)}` : ''}; stiffness K_LL = ${fmt(K.KLL / 1e9)} GN/m, K_RR = ${fmt(K.KRR / 1e9)} GN·m/rad.`,
          gov.uc > 1 ? (gov === checks[0] ? `Reduce rotation: lengthen the pile towards the critical length${crit && crit.critical ? ` (${fmt(crit.critical)} m)` : ''} or increase the diameter by about ${fmt((gov.uc ** 0.4 - 1) * 100)} %; scour protection removes the ${inp.scour} m loss of support.` : gov === checks[1] ? `Increase the wall in the zone ${fmt(Math.max(0, r.depthMaxMoment - D))}–${fmt(r.depthMaxMoment + D)} m below the mudline to about ${fmt(inp.t * gov.uc * 1.03, 3)} mm.` : `Increase embedment or diameter: the soil capacity or length criterion is exceeded by ${fmt((gov.uc - 1) * 100)} %.`) : `The pile is adequate. ${crit && crit.critical && inp.Lp > 1.15 * crit.critical ? `It is ${fmt(inp.Lp - crit.critical)} m longer than the critical length: test a shorter pile.` : 'Use the stiffness matrix in the natural-frequency module and confirm with site-calibrated soil reactions.'}`,
          'DNV-ST-0126 (2021) §7.6 with API RP 2GEO p–y curves; PISA-type soil-reaction models (Byrne et al. 2019) are current practice for large monopiles.',
          `Embedded pile: ${fmt(steelT)} t of steel ≈ ${fmt(steelT * CO2_STEEL)} t CO₂e (screening factor ${CO2_STEEL} t/t). Each metre of embedment is ${fmt(steelT / inp.Lp)} t; a shorter pile also means fewer hammer blows and less underwater noise.`),
        adv('economic', 'info', `The embedded length carries ${fmt(steelT)} t of steel (${fmt(steelT / inp.Lp)} t per metre) and its driving time.`,
          crit && crit.critical && inp.Lp > crit.critical ? `Shortening from ${fmt(inp.Lp)} m to the critical length of ${fmt(crit.critical)} m would save about ${fmt((inp.Lp - crit.critical) * steelT / inp.Lp)} t per position, provided rotation (${fmt(rot.totalDeg)}°) and capacity stay within limits.` : 'Additional site investigation (CPT at every position, cyclic laboratory tests) typically pays back through shorter piles: 1 m less embedment on every position is a material saving across the farm.',
          'Value-of-information reasoning for geotechnical campaigns; Class 4 estimate per AACE RP 18R-97 on pile tonnage and installation vessel time.'),
      ],
      publish: { 'foundation.Kh': K.KLL, 'foundation.Kr': K.KRR, 'foundation.Khr': K.KLR, 'foundation.headRotationDeg': rot.totalDeg, 'quantities.steel_t': steelT, 'quantities.concrete_m3': 0 },
      quantities: { steel_t: round(steelT, 4) },
    };

    function gravityBase(q, L) {
      need(q.soilSpec === 'sand' || q.soilSpec === 'clay', 'The gravity-base check needs uniform sand or uniform clay: choose one of the two soil descriptions.');
      const R = q.gbsD / 2, W = q.gbsW * 1e6, Vt = L.V + W;
      // moment at base level includes the shear acting over the embedment
      const Mb = L.M + L.H * q.gbsEmb, Hd = L.H * q.gammaF, Md = Mb * q.gammaF, Vd = Vt;   // favourable weight unfactored
      const e = Md / Vd;
      need(e < R * 0.98, `The load eccentricity (${fmt(e)} m) is outside the base (radius ${fmt(R)} m): the foundation overturns. Add ballast or enlarge the base.`);
      const ea = effectiveAreaCircle({ R, e });
      const sand = q.soilSpec === 'sand';
      const phiD = sand ? Math.atan(Math.tan(q.phi * DEG) / q.gammaM) : 0, suD = sand ? 0 : (q.su0 + q.suGrad * (q.gbsEmb + ea.Beff / 4)) * 1e3 / q.gammaM;
      const bc = run(() => bearingCapacity({ method: 'vesic', B: ea.Beff, L: ea.Leff, D: q.gbsEmb, c: suD, phi: phiD, gamma: q.gammaSub * 1e3, V: Vd, HB: Hd, useDepthFactors: false }), 'Bearing capacity (the horizontal load may exceed the sliding resistance of the effective area)');
      const sl = sand ? slidingResistance({ mode: 'drained', A: ea.Aeff, V: Vd, delta: Math.atan(Math.tan(q.phi * DEG) / q.gammaM) }) : slidingResistance({ mode: 'undrained', A: ea.Aeff, suBase: (q.su0 + q.suGrad * q.gbsEmb) * 1e3 / q.gammaM });
      const ot = overturningCheck({ V: Vt, M: Mb, B: q.gbsD, shape: 'circle' });
      const cks = [
        check('Bearing capacity on the effective area', Vd / 1e6, bc.Qult / 1e6, 'MN', 'DNV-ST-0126 (2021) §7.5 and App. G (effective area, Vesic factors)', 'Vertical load against the bearing resistance of the reduced contact area under eccentric, inclined loading with factored soil strength. Exceeding it means a bearing failure under the loaded edge in the design storm.'),
        check('Sliding at the base', Hd / 1e6, sl.H / 1e6, 'MN', 'DNV-ST-0126 (2021) §7.5 (sliding resistance)', 'Factored horizontal load against the base shear resistance on the effective area (friction V·tan δ in sand, su·A′ in clay). Above 1 the foundation slides: add weight, skirts or a larger base.'),
        check('Load eccentricity against the no-uplift limit D/8', Mb / Vt, ot.kernLimit, 'm', 'Statics of a rigid circular footing; DNV-ST-0126 §7.5 (contact area under characteristic loads)', 'Eccentricity of the characteristic load against the kern radius. Beyond it the windward edge lifts off each cycle, pumping water and soil from under the base.'),
      ];
      const g = worst(cks);
      const es = Array.from({ length: 30 }, (_, i) => (i + 0.5) / 30 * R * 0.95);
      const qCurve = es.map((ee) => { try { const a = effectiveAreaCircle({ R, e: ee }); return bearingCapacity({ method: 'vesic', B: a.Beff, L: a.Leff, D: q.gbsEmb, c: suD, phi: phiD, gamma: q.gammaSub * 1e3, V: Vd, HB: Hd, useDepthFactors: false }).Qult / 1e6; } catch { return 0; } });
      const co2 = q.gbsConc * CO2_CONCRETE;
      return {
        kpis: [
          { key: 'ucBearing', label: 'Bearing utilisation', value: round(cks[0].uc, 4), unit: '-', limit: 1, sense: 'max', explain: 'Vertical load over the bearing resistance of the effective area with factored soil strength.' },
          { key: 'ucSliding', label: 'Sliding utilisation', value: round(cks[1].uc, 4), unit: '-', limit: 1, sense: 'max', explain: 'Factored horizontal load over the sliding resistance of the base.' },
          { key: 'slideCap', label: 'Design sliding resistance', value: round(sl.H / 1e6, 5), unit: 'MN', explain: 'Base shear resistance on the effective area with factored soil strength.' },
          { key: 'ecc', label: 'Design load eccentricity', value: round(e, 4), unit: 'm', explain: 'Factored overturning moment over vertical load: the offset of the resultant from the centre of the base.' },
          { key: 'Aeff', label: 'Effective bearing area', value: round(ea.Aeff, 5), unit: 'm²', explain: 'Lens-shaped contact area centred on the load resultant; the bearing check uses its equivalent rectangle.' },
          { key: 'fsOverturning', label: 'Factor of safety against overturning (characteristic)', value: round(Math.min(ot.fsOverturning, 999), 4), unit: '-', explain: 'Stabilising over overturning moment about the leeward edge under unfactored loads. Values above about 1.5–2 are normally needed.' },
          { key: 'yHead', label: 'Mean contact pressure on the effective area', value: round(Vd / ea.Aeff / 1e3, 5), unit: 'kPa', explain: 'Vertical load over effective area under design loads (reported under this key so that the discretisation study, which applies to the pile model only, still finds its quantity).' },
        ],
        checks: cks,
        plots: [
          { type: 'line', title: 'Bearing resistance versus load eccentricity', x: { label: 'Eccentricity', unit: 'm' }, y: { label: 'Vertical resistance', unit: 'MN' }, series: [{ name: 'Bearing resistance', x: es.map((v) => round(v, 3)), y: qCurve.map((v) => round(v, 4)) }, { name: 'Vertical load', x: [0, R], y: [Vd / 1e6, Vd / 1e6], dash: true }, { name: 'Design eccentricity', x: [e, e], y: [0, Math.max(...qCurve)], dash: true }],
            explain: 'Bearing resistance falls quickly as the load moves off-centre because the effective area shrinks. The design point is where the two dashed lines cross: it must lie below the curve.' },
        ],
        tables: [{ title: 'Gravity-base load and resistance summary', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [{ q: 'Vertical load (turbine + base)', v: round(Vt / 1e6, 4), u: 'MN' }, { q: 'Design horizontal load', v: round(Hd / 1e6, 4), u: 'MN' }, { q: 'Design moment at base level', v: round(Md / 1e6, 4), u: 'MN·m' }, { q: 'Effective width × length', v: `${fmt(ea.Beff)} × ${fmt(ea.Leff)}`, u: 'm' }, { q: 'Bearing resistance', v: round(bc.Qult / 1e6, 4), u: 'MN' }, { q: 'Sliding resistance', v: round(sl.H / 1e6, 4), u: 'MN' }, { q: 'No-uplift eccentricity D/8', v: round(ot.kernLimit, 4), u: 'm' }],
          explain: 'Loads at the base and the resistances of the gravity foundation. Weight is the stabilising action, so it is not factored up.' }],
        notes: [`Gravity base ${fmt(q.gbsD)} m diameter, ${fmt(q.gbsW)} MN submerged weight on uniform ${q.soilSpec}; governing check: ${g.name} (utilisation ${fmt(g.uc)}).`, 'Cyclic degradation, settlement, scour undermining and installation tolerances need separate assessment.'],
        advice: [
          adv('engineering', sev(g.uc), `Bearing utilisation ${fmt(cks[0].uc)}, sliding ${fmt(cks[1].uc)}, eccentricity ${fmt(Mb / Vt)} m against a no-uplift limit of ${fmt(ot.kernLimit)} m.`,
            g.uc > 1 ? `Increase the stabilising weight to about ${fmt(q.gbsW * Math.max(g.uc, 1.05))} MN or the base diameter to about ${fmt(q.gbsD * Math.max(g.uc, 1.05) ** 0.4)} m; skirts raise sliding resistance on clay.` : 'The base is stable; check settlement, cyclic strength and scour protection around the perimeter.',
            'DNV-ST-0126 (2021) §7.5 and Appendix G; API RP 2GEO §7 (effective-area method).',
            `Concrete of the base: ${fmt(q.gbsConc)} m³ ≈ ${fmt(co2)} t CO₂e (screening factor ${CO2_CONCRETE} t/m³). Gravity bases avoid piling noise and can be removed at end of life.`),
          adv('economic', 'info', `A gravity base trades ${fmt(q.gbsConc)} m³ of concrete and ${fmt(q.gbsW)} MN of ballast for piling: no hammer spread, but seabed preparation and a heavy-lift or float-out operation.`,
            'Compare on installed cost per position with the monopile option for this soil: gravity bases are competitive on hard or rocky seabed where driving is not possible.', 'Concept selection by whole-life cost (ISO 15686-5) with Class 4 quantities (AACE RP 18R-97).'),
        ],
        publish: { 'quantities.concrete_m3': q.gbsConc, 'quantities.steel_t': 0 },
        quantities: { concrete_m3: q.gbsConc },
      };
    }
  },
  verify: [
    { name: 'Long pile on linear springs under a head force: y₀ = 2Hβ/k, β = (k/4EI)^¼', source: 'Hetényi (1946) Beams on Elastic Foundation, semi-infinite beam with a free end',
      // D = 6 m, t = 60 mm: I = π/64·(6⁴ − 5.88⁴) = 4.9387 m⁴, EI = 1.0371e12 N·m²; k = 50 MN/m² → β = (50e6/4.1485e12)^¼ = 0.058921 1/m (βL = 4.7 for L = 80 m)
      // y₀ = 2·5e6·0.058921/50e6 = 11.784 mm; K_LL of the long pile = 4EIβ³·... not used here
      inputs: { concept: 'monopile', soilSpec: 'linear', kLin: 50, D: 6, t: 60, Lp: 80, H: 5, M: 0, V: 0, scour: 0, steel: 's355', nElem: 80 }, expect: [{ kpi: 'yHead', value: 11.784, tol: 0.01 }] },
    { name: 'Gravity base on sand: design sliding resistance V·tan(δd) for a centric load', source: 'Coulomb friction on the full base; DNV-ST-0126 §7.5',
      // V = 10 + 90 = 100 MN, tan δd = tan 35° / 1.15 = 0.60887 → 60.89 MN (M = 0: effective area = full area)
      inputs: { concept: 'gravity', soilSpec: 'sand', phi: 35, gbsW: 90, V: 10, H: 5, M: 0, gbsEmb: 0, gammaM: 1.15 }, expect: [{ kpi: 'slideCap', value: 60.887, tol: 0.002 }] },
  ],
};
