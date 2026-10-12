// Touchdown and seabed interaction of a cable: catenary touchdown geometry, bending boundary layer and touchdown
// reaction, contact pressure and embedment, axial friction anchorage of the residual tension, and lateral
// on-bottom stability under near-bed current and wave-induced velocity.
import { catenaryRiser, boundaryLayerTdp } from '../../engine/lines/riser.js';
import { embedment, axialResistance, lateralResistance } from '../../engine/pipeline/pipesoil.js';
import { nearBedSpectralVelocity, designOscillation } from '../../engine/pipeline/stability.js';
import { cableOnBottomStability } from '../../engine/xsection/cablemech.js';
import { lineInputs, lineFrom, adv, check, worst, sev, fmt, round, run, need, linspace, G, STD_CABLE, CO2_ROCK, vesselCo2Text } from './_common.js';

export default {
  id: 's13-touchdown-seabed', suite: 13, title: 'Cable touchdown and seabed interaction',
  summary: 'Touchdown position and curvature, contact force and pressure, embedment, friction anchorage and on-bottom stability of a cable on the seabed.',
  covers: [5], fidelity: 'analytical',
  standards: ['DNV-RP-F114 (2021) §4 (pipe–soil interaction: embedment, axial and lateral resistance)', 'DNV-RP-F109 (2021) §3 (on-bottom stability, near-bed spectral velocity)', `${STD_CABLE}`, 'Sparks, Fundamentals of Marine Riser Mechanics (2007) ch. 13 (touchdown boundary layer)'],
  models: ['Nonlinear catenary equations', 'Tensioned Euler–Bernoulli beam boundary layer at touchdown', 'Beam on elastic foundation (Hetényi) for the peak contact reaction', 'Seabed contact mechanics: bearing capacity and embedment', 'Coulomb frictional sliding with passive soil resistance', 'Morison drag and lift on a cylinder on the seabed', 'JONSWAP near-bed spectral velocity'],
  limits: 'Flat uniform seabed, static touchdown (no trenching by cyclic motion, no liquefaction or scour); the lateral check is the absolute-stability form with user drag and lift coefficients, not the calibrated DNV-RP-F109 peak-load tables. Seabed mobility and free spans are treated in the burial and free-span modules.',
  publishes: ['cable.embedment', 'cable.lateralResistance', 'cable.tdpCurvature'],
  searchTerms: 'submarine cable touchdown seabed interaction embedment on-bottom stability',
  inputs: [
    ...lineInputs({ hydro: false }),
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 80, min: 3, max: 3000, group: 'Touchdown', role: 'environment', bind: 'site.depth' },
    { key: 'H', label: 'Residual bottom tension at touchdown', unit: 'kN', type: 'number', default: 8, min: 0.05, max: 2000, group: 'Touchdown', role: 'ic', bind: 'cable.tdpTension', design: true, uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Horizontal tension left in the cable at touchdown (laying or installed): the initial pretension of the seabed section.' },
    { key: 'Lbed', label: 'Cable length available on the seabed', unit: 'm', type: 'number', default: 500, min: 1, max: 5e5, group: 'Touchdown', role: 'bc', help: 'Length between touchdown and the next fixed point (burial start, clump weight, structure): it must anchor the bottom tension by friction.' },
    { key: 'soil', label: 'Seabed soil', type: 'select', default: 'sand', group: 'Seabed', role: 'soil', options: [{ value: 'sand', label: 'Sand (drained)' }, { value: 'clay', label: 'Clay (undrained)' }] },
    { key: 'phi', label: 'Friction angle', unit: 'deg', type: 'number', default: 32, min: 20, max: 45, group: 'Seabed', role: 'soil', bind: 'site.phiDeg', showIf: { key: 'soil', in: ['sand'] } },
    { key: 'su', label: 'Undrained shear strength at the mudline', unit: 'kPa', type: 'number', default: 5, min: 0.3, max: 300, group: 'Seabed', role: 'soil', bind: 'site.su0', showIf: { key: 'soil', in: ['clay'] }, uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'suGrad', label: 'Strength gradient', unit: 'kPa/m', type: 'number', default: 1.5, min: 0, max: 30, group: 'Seabed', role: 'soil', bind: 'site.suGrad', showIf: { key: 'soil', in: ['clay'] } },
    { key: 'gammaSub', label: 'Submerged unit weight of soil', unit: 'kN/m³', type: 'number', default: 9, min: 2, max: 13, group: 'Seabed', role: 'soil', bind: 'site.gammaSub' },
    { key: 'kSoil', label: 'Seabed vertical stiffness', unit: 'kPa', type: 'number', default: 500, min: 5, max: 1e5, group: 'Seabed', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.5 }, help: 'Line load per metre of penetration (kN/m/m); sets the spread of the touchdown reaction.' },
    { key: 'muAx', label: 'Axial friction coefficient', unit: '-', type: 'number', default: 0.5, min: 0.05, max: 1.5, group: 'Seabed', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Sand 0.4–0.7; for clay the adhesion factor on the undrained strength.' },
    { key: 'muLat', label: 'Lateral friction coefficient', unit: '-', type: 'number', default: 0.5, min: 0.05, max: 1.5, group: 'Seabed', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'DNV-RP-F109 uses 0.6 for sand and 0.2 for clay with a separate passive term.' },
    { key: 'kDyn', label: 'Dynamic embedment factor', unit: '-', type: 'number', default: 1.5, min: 1, max: 5, group: 'Seabed', role: 'soil', help: 'Extra embedment from lay-vessel motion working the cable into the seabed, typically 1–3 (DNV-RP-F114).' },
    { key: 'Uc', label: 'Near-bed current normal to the cable', unit: 'm/s', type: 'number', default: 0.4, min: 0, max: 4, group: 'Environment', role: 'bc', bind: 'site.UcBed100', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'Hs', label: 'Significant wave height', unit: 'm', type: 'number', default: 6, min: 0, max: 20, group: 'Environment', role: 'bc', bind: 'site.Hs100' },
    { key: 'Tp', label: 'Peak period', unit: 's', type: 'number', default: 11, min: 3, max: 25, group: 'Environment', role: 'bc', bind: 'site.Tp100' },
    { key: 'Cd', label: 'Drag coefficient on the seabed', unit: '-', type: 'number', default: 1.0, min: 0.3, max: 2.5, group: 'Environment', role: 'load', help: 'Cylinder resting on the bed in steady flow: about 0.7–1.2.' },
    { key: 'Cl', label: 'Lift coefficient on the seabed', unit: '-', type: 'number', default: 0.9, min: 0, max: 2, group: 'Environment', role: 'load' },
    { key: 'sf', label: 'Safety factor on lateral stability', unit: '-', type: 'number', default: 1.1, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'mbrAllow', label: 'Allowable bend radius at touchdown', unit: 'm', type: 'number', default: 3.2, min: 0.3, max: 60, group: 'Acceptance criteria', role: 'criteria', bind: 'cable.mbrDynamic' },
    { key: 'pAllow', label: 'Allowable contact line load on the cable', unit: 'kN/m', type: 'number', default: 20, min: 0.5, max: 500, group: 'Acceptance criteria', role: 'criteria', help: 'From the crush / sidewall-pressure test of the cable.' },
    { key: 'zMaxOverD', label: 'Maximum acceptable embedment / diameter', unit: '-', type: 'number', default: 1.0, min: 0.1, max: 5, group: 'Acceptance criteria', role: 'criteria', help: 'Self-burial beyond one diameter changes the thermal rating and makes recovery harder.' },
  ],
  compute(inp) {
    const line = lineFrom(inp), w = line.w, D = line.D, H = inp.H * 1e3;
    need(w > 0, `The cable is buoyant (apparent weight ${fmt(w)} N/m): it does not rest on the seabed.`);
    const cat = run(() => catenaryRiser({ h: inp.depth, w, H }), 'Catenary');
    const xi = linspace(0, 8 * Math.sqrt(line.EI / H), 60);
    const bl = run(() => boundaryLayerTdp({ EI: line.EI, H, w, xi }), 'Touchdown boundary layer');
    // peak contact line load: distributed weight plus the concentrated touchdown reaction spread by the elastic seabed
    // (Hetényi: point load P on an infinite beam on elastic foundation → p_max = P·β/2, β = (k/4EI)^¼)
    const k = inp.kSoil * 1e3, beta = (k / (4 * line.EI)) ** 0.25;
    const qPeak = w + bl.pointReaction * beta / 2, kLay = qPeak / w;
    const clay = inp.soil === 'clay';
    const soilArgs = clay ? { model: 'undrained', su0: inp.su * 1e3, suGrad: inp.suGrad * 1e3, gammaSub: inp.gammaSub * 1e3 } : { model: 'drained', phiDeg: inp.phi, gammaSub: inp.gammaSub * 1e3 };
    const emb = run(() => embedment({ D, W: w, kLay, kDyn: inp.kDyn, ...soilArgs }), 'Embedment (the soil may be too weak to carry the cable within the search depth)');
    const z = Math.max(emb.z, 1e-6), B = emb.geometry?.B ?? D;
    const suZ = (inp.su + inp.suGrad * z) * 1e3;
    const ax = clay ? axialResistance({ model: 'undrained', D, z, alpha: Math.min(inp.muAx, 1), su: suZ }) : axialResistance({ model: 'drained', D, z, W: w, mu: inp.muAx });
    // hydrodynamic demand: near-bed current plus design single oscillation of the storm
    const nb = inp.Hs > 0 ? nearBedSpectralVelocity({ Hs: inp.Hs, Tp: inp.Tp, depth: inp.depth }) : { Us: 0, Tu: inp.Tp, Tn: Math.sqrt(inp.depth / G) };
    const osc = nb.Us > 1e-6 ? designOscillation({ Us: nb.Us, Tu: nb.Tu, Tn: nb.Tn }) : { Ustar: 0 };
    const U = inp.Uc + osc.Ustar;
    const stab = run(() => cableOnBottomStability({ submergedWeight: w, outerDiameter: D, velocity: U, mu: inp.muLat, Cd: inp.Cd, Cl: inp.Cl, rho: line.rhoW, safetyFactor: inp.sf, mass: line.m }), 'On-bottom stability');
    const Wn = Math.max(w - stab.lift, 0);
    // when lift exceeds the weight the cable leaves the seabed: a vanishing contact force keeps the resistance model defined (utilisation capped below)
    const lat = lateralResistance({ soil: inp.soil, D, z, W: Math.max(Wn, 1e-3), mu: inp.muLat, su: suZ, gammaSub: inp.gammaSub * 1e3, phiDeg: inp.phi });
    const ucLat = lat.breakout > 0 ? Math.min(inp.sf * stab.drag / lat.breakout, 99) : 99;
    const Lanchor = H / Math.max(ax.breakout, 1e-9);
    const Rtdp = 1 / bl.catenaryCurvature;
    const checks = [
      check('Bend radius at touchdown', inp.mbrAllow, Rtdp, 'm', `${STD_CABLE} (minimum bend radius under tension)`, 'Demand is the allowable radius, capacity the catenary radius H/w at touchdown. Too little bottom tension lets the cable kink where it lands.'),
      check('Peak contact line load', qPeak / 1e3, inp.pAllow, 'kN/m', 'Hetényi, Beams on Elastic Foundation (1946); cable crush test', 'Seabed reaction per metre at touchdown (weight plus the concentrated boundary-layer reaction) against the tested crush load of the cable.'),
      check('Embedment', z / D, inp.zMaxOverD, '× D', 'DNV-RP-F114 §4.2 (embedment)', 'Predicted self-embedment relative to the diameter against the acceptable value. Deep embedment in soft clay raises the thermal resistance around the cable and the recovery load.'),
      check('Lateral stability (friction only, absolute)', Math.min(stab.utilisation, 99), 1, '-', 'DNV-RP-F109 §3.6 (absolute lateral static stability) in form; Morison drag and lift', 'Factored drag against Coulomb friction on the weight less lift. Above 1.0 the cable is predicted to slide sideways in the design storm, abrading the serving and dragging on crossings and joints.'),
      check('Lateral stability with passive soil resistance', ucLat, 1, '-', 'DNV-RP-F114 §4.4; DNV-RP-F109 §3.4.6 (passive resistance)', 'Same demand against friction plus the passive resistance of the soil ridge in front of the embedded cable. It shows how much the embedment helps.'),
      check('Friction anchorage of the bottom tension', Lanchor, inp.Lbed, 'm', 'Coulomb friction; DNV-RP-F114 §4.3 (axial resistance)', 'Length of cable on the seabed needed for axial friction to hold the residual tension, against the available length. If it is longer, the tension reaches the next fixed point (joint, burial transition or structure).'),
    ];
    const gov = worst(checks);
    const Us = linspace(0, Math.max(2 * stab.criticalVelocity, U * 1.3, 0.5), 40);
    const ucOfU = Us.map((u) => { const q = 0.5 * line.rhoW * D * u * u, r = inp.muLat * Math.max(w - inp.Cl * q, 0); return r > 0 ? Math.min(inp.sf * inp.Cd * q / r, 5) : 5; });
    const xs = linspace(0, cat.horizontalSpan, 50), a = cat.a;
    const reqExtra = Math.max(0, stab.requiredWeight - w);
    return {
      kpis: [
        { key: 'tdpRadius', label: 'Bend radius at touchdown', value: round(Rtdp), unit: 'm', limit: inp.mbrAllow, sense: 'min', explain: 'Catenary radius H/w where the cable meets the seabed; the bending stiffness makes the true curvature rise to this value over a few boundary-layer lengths.' },
        { key: 'lambda', label: 'Boundary-layer length', value: round(bl.lambda), unit: 'm', explain: 'Decay length √(EI/H) of the bending effect at touchdown. Within about three lengths the curvature differs from the catenary value.' },
        { key: 'layback', label: 'Touchdown distance from the surface point', value: round(cat.horizontalSpan), unit: 'm', explain: 'Horizontal distance from the top of the suspended span to the touchdown point (layback): the touchdown location for the given depth and bottom tension.' },
        { key: 'suspended', label: 'Suspended length', value: round(cat.suspendedLength), unit: 'm', explain: 'Length of cable hanging between the surface point and the touchdown.' },
        { key: 'topTension', label: 'Top tension of the suspended span', value: round(cat.topTension / 1e3), unit: 'kN', explain: 'H + w·d at the surface for this bottom tension and depth.' },
        { key: 'tdpReaction', label: 'Concentrated touchdown reaction', value: round(bl.pointReaction), unit: 'N', explain: 'Extra seabed reaction at the touchdown point caused by bending stiffness (w·λ): the contact force in excess of the distributed weight.' },
        { key: 'qPeak', label: 'Peak contact line load', value: round(qPeak), unit: 'N/m', explain: 'Largest seabed reaction per metre: cable weight plus the touchdown reaction spread over the elastic seabed.' },
        { key: 'contactPressure', label: 'Peak contact pressure', value: round(qPeak / Math.max(B, 1e-4) / 1e3), unit: 'kPa', explain: 'Peak line load divided by the contact width at the predicted embedment; compare with the soil bearing strength and the cable serving rating.' },
        { key: 'embedment', label: 'Embedment', value: round(z * 1e3), unit: 'mm', explain: 'Depth the cable sinks into the seabed under the touchdown load concentration and lay dynamics (bearing-capacity equilibrium).' },
        { key: 'axialResistance', label: 'Axial resistance', value: round(ax.breakout), unit: 'N/m', explain: 'Friction (sand) or adhesion (clay) per metre resisting axial sliding of the cable.' },
        { key: 'anchorLength', label: 'Friction anchor length', value: round(Lanchor), unit: 'm', limit: inp.Lbed, sense: 'max', explain: 'Seabed length needed to hold the bottom tension by axial friction, H / f_axial.' },
        { key: 'Udesign', label: 'Design near-bed velocity', value: round(U), unit: 'm/s', explain: 'Steady current plus the largest wave-induced oscillation of the storm at the seabed, normal to the cable.' },
        { key: 'Ucrit', label: 'Critical velocity for sliding', value: round(stab.criticalVelocity), unit: 'm/s', explain: 'Near-bed velocity at which factored drag equals friction on the weight less lift. Compare with the design velocity.' },
        { key: 'ucLateral', label: 'Lateral stability utilisation (friction only)', value: round(Math.min(stab.utilisation, 99)), unit: '-', limit: 1, sense: 'max', explain: 'Factored drag over frictional resistance; above 1.0 the unburied cable moves in the design storm.' },
        { key: 'wRequired', label: 'Submerged weight required for stability', value: round(stab.requiredWeight), unit: 'N/m', explain: 'Weight at which the unburied cable is just stable; the shortfall must come from burial, rock or mattresses.' },
        { key: 'sg', label: 'Specific gravity', value: round(stab.specificGravity), unit: '-', explain: 'Cable mass over displaced water. Power cables are typically 2–4; values near 1 float up in liquefied backfill.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Suspended catenary to touchdown', x: { label: 'Distance from touchdown', unit: 'm' }, y: { label: 'Height above seabed', unit: 'm' }, series: [{ name: 'Cable', x: xs, y: xs.map((x) => a * (Math.cosh(x / a) - 1)) }],
          explain: 'Static shape between the touchdown (origin) and the surface for the given bottom tension. A larger bottom tension flattens the shape and moves the touchdown away.' },
        { type: 'line', title: 'Curvature in the touchdown boundary layer', x: { label: 'Distance above touchdown along the cable', unit: 'm' }, y: { label: 'Curvature', unit: '1/m' },
          series: [{ name: 'With bending stiffness', x: xi, y: bl.curvature }, { name: 'Catenary value w/H', x: [0, xi[xi.length - 1]], y: [bl.catenaryCurvature, bl.catenaryCurvature], dash: true }, { name: 'Allowable 1/MBR', x: [0, xi[xi.length - 1]], y: [1 / inp.mbrAllow, 1 / inp.mbrAllow], dash: true }],
          explain: 'Curvature rises from zero on the seabed to the catenary value over a few boundary-layer lengths; it must stay below the allowable line.' },
        { type: 'line', title: 'Lateral stability against near-bed velocity', x: { label: 'Velocity normal to the cable', unit: 'm/s' }, y: { label: 'Utilisation', unit: '-' },
          series: [{ name: 'Factored drag / friction', x: Us, y: ucOfU }, { name: 'Limit', x: [0, Us[Us.length - 1]], y: [1, 1], dash: true }, { name: 'Design velocity', x: [U, U], y: [0, Math.min(Math.max(stab.utilisation, 1), 5)], dash: true }],
          explain: 'Utilisation grows with velocity squared and steepens as lift unloads the cable. The crossing of the limit is the critical velocity; the design velocity should lie to its left.' },
      ],
      tables: [{ title: 'Seabed resistances', explain: 'Each row is one resistance mechanism per metre of cable with its value, for use as spring limits in a global model.',
        columns: [{ key: 'item', label: 'Mechanism' }, { key: 'v', label: 'Resistance', unit: 'N/m' }, { key: 'note', label: 'Basis' }],
        rows: [
          { item: 'Axial breakout', v: round(ax.breakout), note: clay ? 'adhesion × strength × contact arc' : 'μ × weight × wedging factor' },
          { item: 'Axial residual', v: round(ax.residual), note: 'after sliding' },
          { item: 'Lateral friction (weight less lift)', v: round(lat.frictionBreakout), note: 'Coulomb' },
          { item: 'Lateral passive', v: round(lat.passive), note: 'soil ridge at the predicted embedment' },
          { item: 'Lateral breakout total', v: round(lat.breakout), note: 'friction + passive' },
          { item: 'Hydrodynamic drag (design)', v: round(stab.drag), note: `U = ${fmt(U)} m/s` },
          { item: 'Hydrodynamic lift (design)', v: round(stab.lift), note: `U = ${fmt(U)} m/s` },
        ] }],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        `Near-bed velocity: current ${inp.Uc} m/s + wave oscillation ${fmt(osc.Ustar)} m/s (significant ${fmt(nb.Us)} m/s, period ${fmt(nb.Tu)} s).`,
        `Touchdown load concentration factor ${fmt(kLay)}; embedment ${fmt(z / D)} diameters in ${inp.soil}.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}. Design velocity ${fmt(U)} m/s against a critical velocity of ${fmt(stab.criticalVelocity)} m/s; embedment ${fmt(z * 1e3)} mm; touchdown radius ${fmt(Rtdp)} m.`,
          gov.uc > 1 ? (gov.name.startsWith('Lateral') ? `The unburied cable needs ${fmt(reqExtra)} N/m more submerged weight, or burial / rock cover over the exposed length; otherwise accept movement only with a dedicated dynamic stability analysis.` : gov.name.startsWith('Bend') ? `Raise the bottom tension to at least ${fmt(w * inp.mbrAllow / 1e3)} kN during laying.` : gov.name.startsWith('Friction') ? `Provide a holdback (clump weight, burial transition or anchor clamp) rated for ${fmt(Math.max(0, H - ax.breakout * inp.Lbed) / 1e3)} kN.` : 'Reduce the touchdown load (lower bottom tension, slower lay) or confirm the cable crush capacity.')
            : 'Confirm soil parameters with site CPT / box-core data and repeat for the weakest and strongest soil units along the route; unburied lengths need this check for every exposure.',
          'DNV-RP-F109 (2021) §3; DNV-RP-F114 (2021) §4; DNV-RP-0360 (2016) §4'),
        adv('economic', reqExtra > 0 ? 'warn' : 'info', reqExtra > 0 ? `Stability by weight alone would need ${fmt(reqExtra / G)} kg/m more submerged mass — about ${fmt(reqExtra / G / 1e3 * 1000)} t per kilometre.` : `The cable is stable unburied with ${fmt((1 - stab.utilisation) * 100)} % reserve in the design storm.`,
          reqExtra > 0 ? 'Burial is normally cheaper than added armour: a jetting or plough spread covers 100–400 m/h, whereas rock placement needs of the order of 2–4 t/m. Compare vessel days rather than material cost.' : 'Where the reserve holds along the route, surface lay without rock avoids a fall-pipe vessel campaign; keep burial for protection against fishing and anchors only.',
          'AACE RP 18R-97 (quantity-based estimate); whole-life costing ISO 15686-5',
          `Rock cover at about 3 t/m ≈ ${fmt(3000 * CO2_ROCK)} t CO₂e per km of rock (screening ${CO2_ROCK} t CO₂e/t) plus the fall-pipe vessel: ${vesselCo2Text(1)} Burial by jetting leaves no imported material on the seabed.`),
      ],
      publish: { 'cable.embedment': z, 'cable.lateralResistance': lat.breakout, 'cable.tdpCurvature': bl.catenaryCurvature },
      quantities: {},
    };
  },
  verify: [
    // w = (63 − 1025·π·0.16²/4)·9.80665 = 415.73 N/m. H = 10 kN: touchdown radius H/w = 24.05 m; λ = sqrt(EI/H) = sqrt(12000/10000) = 1.0954 m;
    // point reaction w·λ = 455.4 N. Sliding: ½ρD(Cl + SF·Cd/μ) = 0.5·1025·0.16·(0.9 + 1.1·1.0/0.5) = 254.2 → Ucrit = sqrt(415.73/254.2) = 1.2788 m/s.
    { name: 'Touchdown radius H/w, boundary layer √(EI/H), and critical sliding velocity',
      source: 'Irvine (1981) §1.2; Sparks (2007) ch. 13; Morison drag/lift with Coulomb friction (hand calculation)',
      inputs: { H: 10, Hs: 0, Uc: 0.5 },
      expect: [{ kpi: 'tdpRadius', value: 24.05, tol: 0.003 }, { kpi: 'lambda', value: 1.0954, tol: 0.003 }, { kpi: 'tdpReaction', value: 455.4, tol: 0.005 }, { kpi: 'Ucrit', value: 1.2788, tol: 0.003 }] },
  ],
};
