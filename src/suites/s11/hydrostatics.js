// Hydrostatics and floating equilibrium of a closed hull surface: displacement, centre of buoyancy, waterplane
// properties, metacentric heights with free-surface correction, restoring stiffness and the hydrostatic curves.
import { meshHydrostatics, solveEquilibrium, freeSurfaceCorrection } from '../../engine/hydro/hydrostatics.js';
import { hullInputs, hullFrom, hullPlot, hullEnclosedVolume, G, DEG, need, run, round, fmt, adv, checkMin, check, sev, worst, linspace, steelCo2Text } from './_common.js';

export default {
  id: 's11-hydrostatics',
  suite: 11,
  title: 'Hydrostatics and floating equilibrium',
  summary: 'Displacement, buoyancy centre, waterplane properties, metacentric height with free-surface correction, equilibrium draft, heel and trim and the restoring stiffness of an imported or parametric hull.',
  covers: [1],
  fidelity: 'numerical',
  standards: ['IMO 2008 IS Code (MSC.267(85)) Part A §2.2, Part B ch. 3', 'IMO MODU Code 2009 §3.2', 'DNV-OS-C301 (2021) Ch.2 Sec.1', 'DNV-RP-C205 (2019) §7.1'],
  models: ['Archimedes\' principle', 'Exact volume and waterplane integrals of a closed triangulated surface (divergence theorem)', 'Metacentric stability theory (BM = I/∇)', 'Free-surface correction GG′ = Σρᵢi/(ρ∇)', 'Nonlinear floating equilibrium (Newton iteration on draft, heel and trim)', 'Linear hydrostatic restoring matrix'],
  limits: 'Calm-water hydrostatics of the faceted surface: exact for the mesh, so a coarse imported mesh is the main error. Free-surface correction is the small-angle wall-sided value. For a tension-leg platform the free-floating metacentric height is reported without the tendons, which provide the actual restoring. Hull deformation and trapped air are not modelled.',
  publishes: ['vessel.disp', 'vessel.GM', 'quantities.steel_t'],
  searchTerms: 'floating structure hydrostatics metacentric height equilibrium',
  inputs: [
    ...hullInputs({ type: 'semi' }),
    { key: 'massT', label: 'Mass of the unit (0 = from the draft)', unit: 't', type: 'number', default: 0, min: 0, max: 2e6, group: 'Loading condition', role: 'ic',
      help: 'Leave 0 to float the unit at the draft above. Enter the weight-report mass to solve the equilibrium draft instead.' },
    { key: 'lcgOff', label: 'Longitudinal offset of G from the centre of buoyancy', unit: 'm', type: 'number', default: 0, min: -30, max: 30, group: 'Loading condition', role: 'ic', help: 'Positive forward. A non-zero value trims the unit.' },
    { key: 'tcgOff', label: 'Transverse offset of G from the centreline', unit: 'm', type: 'number', default: 0, min: -15, max: 15, group: 'Loading condition', role: 'ic', help: 'Positive to port. A non-zero value heels the unit.' },
    { key: 'tanks', label: 'Slack tanks (free surfaces)', type: 'table', group: 'Loading condition', role: 'ic',
      columns: [{ key: 'name', label: 'Tank', type: 'text' }, { key: 'l', label: 'Surface length', unit: 'm', type: 'number' }, { key: 'b', label: 'Surface breadth', unit: 'm', type: 'number' }, { key: 'rho', label: 'Liquid density', unit: 'kg/m³', type: 'number' }],
      default: [{ name: 'Ballast P2', l: 18, b: 15, rho: 1025 }, { name: 'Ballast S2', l: 18, b: 15, rho: 1025 }, { name: 'Fuel oil', l: 10, b: 8, rho: 870 }],
      help: 'Partly filled tanks only: free-surface length (along x) and breadth (across). Full or empty tanks carry no correction.' },
    { key: 'gmMin', label: 'Required metacentric height', unit: 'm', type: 'number', default: 1.0, min: 0, max: 10, group: 'Acceptance criteria', role: 'criteria',
      help: 'Minimum GM after free-surface correction. IMO IS Code asks 0.15 m for ships; MODU operating manuals commonly require 1.0 m for column-stabilised units (0.3 m in temporary conditions).' },
    { key: 'freeboardMin', label: 'Required freeboard to the hull top', unit: 'm', type: 'number', default: 1.5, min: 0, max: 40, group: 'Acceptance criteria', role: 'criteria', help: 'Minimum height of the watertight hull top (deck edge or column top) above the still waterline.' },
    { key: 'heelMax', label: 'Allowable static heel / trim', unit: 'deg', type: 'number', default: 2, min: 0.1, max: 15, group: 'Acceptance criteria', role: 'criteria', help: 'Operating limit on the static inclination of the loaded unit.' },
    { key: 'nTheta', label: 'Facets around each circular column', unit: '-', type: 'int', default: 48, min: 8, max: 512, group: 'Numerical controls', role: 'control', help: 'Circumferential resolution of the parametric cylinders (area-preserving polygons).' },
  ],
  mesh: { key: 'nTheta', label: 'Facets per column', levels: [12, 24, 48, 96], qoi: 'KMt' },
  compute(inp) {
    const hull = hullFrom(inp, { nTheta: inp.nTheta });
    const rho = hull.rho; const mesh = hull.mesh;
    const mass = inp.massT > 0 ? inp.massT * 1e3 : hull.mass;
    const cog = [hull.cog[0] + inp.lcgOff, hull.cog[1] + inp.tcgOff, inp.KG];
    let draft = hull.T; let heel = 0; let trim = 0;
    const solved = inp.massT > 0 || inp.lcgOff !== 0 || inp.tcgOff !== 0;
    if (solved && hull.type !== 'tlp') {
      const full = hullEnclosedVolume(hull) * rho;
      need(mass < full, `The unit sinks: mass ${fmt(mass / 1e3)} t exceeds the buoyancy of the fully immersed hull (${fmt(full / 1e3)} t). Reduce the mass or enlarge the hull.`);
      const eq = run(() => solveEquilibrium({ mesh, mass, cog, rho, guess: { draft: hull.T } }), 'Floating equilibrium');
      need(eq.converged, 'No upright floating equilibrium was found near the starting attitude: the unit is probably unstable (negative GM) or the centre of gravity offset is too large. Lower KG or reduce the offsets.');
      draft = eq.draft; heel = eq.heel; trim = eq.trim;
    }
    const hs = meshHydrostatics({ mesh, draft, heel, trim, rho, cog, mass: hull.type === 'tlp' ? hull.mass : mass });
    const tanks = (inp.tanks || []).filter((t) => t.l > 0 && t.b > 0).map((t) => ({ l: t.l, b: t.b, rho: t.rho > 0 ? t.rho : 1025 }));
    const fsc = tanks.length ? freeSurfaceCorrection({ tanks, volume: hs.volume, rho }) : { GGt: 0, GGl: 0, perTank: [] };
    const GMt = hs.GMt - fsc.GGt; const GMl = hs.GMl - fsc.GGl;
    const freeboard = hull.depth - draft;
    const reserve = hullEnclosedVolume(hull) / hs.volume - 1;
    const incl = Math.max(Math.abs(heel), Math.abs(trim)) / DEG;
    const isTlp = hull.type === 'tlp';
    const checks = [
      isTlp ? check('Tendon pretension positive', 0.05 * hs.buoyancy / 1e6, hull.pretension / 1e6, 'MN', 'API RP 2T (3rd ed., 2010) §5 — minimum tendon tension',
        'A tension-leg platform is held upright by tendon tension, not by its waterplane: pretension must stay a healthy share of the buoyancy (5 % is used here as the floor) or tendons go slack in waves.')
        : checkMin('Metacentric height after free-surface correction', GMt, inp.gmMin, 'm', 'IMO 2008 IS Code Part A §2.2.4; IMO MODU Code 2009 §3.2',
          'Guards against a tender or unstable unit: with GM at or below zero the unit lolls or capsizes under the smallest heeling moment.'),
      checkMin('Freeboard to the top of the watertight hull', freeboard, inp.freeboardMin, 'm', 'International Convention on Load Lines 1966 / MODU Code 2009 §3.7 (project freeboard)',
        'Guards against loss of reserve buoyancy: with too little freeboard the deck edge or column top immerses at small heel and the righting lever collapses.'),
      check('Static heel / trim', incl, inp.heelMax, 'deg', 'Operating manual limit (project criterion)',
        'Guards against operating with a permanent list: equipment, drainage and process trains are designed for a small static inclination, and a list eats into the stability range.'),
    ];
    const gov = worst(checks);
    // hydrostatic curves over the practicable draft range
    const drafts = linspace(Math.max(0.15 * hull.depth, 0.3), 0.95 * hull.depth, 16);
    const curve = drafts.map((d) => { const h = meshHydrostatics({ mesh, draft: d, rho, cog }); return { d, disp: h.displacement / 1e3, KB: h.KB, KMt: h.KMt, KMl: h.KMl, Aw: h.waterplaneArea, TPC: h.TPC }; });
    const K = hs.stiffness;
    const steel = 0.16 * hullEnclosedVolume(hull); // screening structural density of floating hulls, t per m³ enclosed
    const kgMax = hs.KMt - fsc.GGt - inp.gmMin;
    const rows = [
      ['Displaced volume', hs.volume, 'm³'], ['Displacement', hs.displacement / 1e3, 't'], ['Buoyancy', hs.buoyancy / 1e6, 'MN'], ['Draft', draft, 'm'], ['Heel', heel / DEG, 'deg'], ['Trim', trim / DEG, 'deg'],
      ['LCB (from midship)', hs.cobBody[0], 'm'], ['TCB', hs.cobBody[1], 'm'], ['KB', hs.KB, 'm'], ['Waterplane area', hs.waterplaneArea, 'm²'], ['LCF', hs.cofBody[0], 'm'],
      ['Waterplane Ixx (heel)', hs.Ixx, 'm⁴'], ['Waterplane Iyy (trim)', hs.Iyy, 'm⁴'], ['BMt', hs.BMt, 'm'], ['BMl', hs.BMl, 'm'], ['KMt', hs.KMt, 'm'], ['KMl', hs.KMl, 'm'], ['KG', inp.KG, 'm'],
      ['Free-surface rise GG′ (transverse)', fsc.GGt, 'm'], ['GMt corrected', GMt, 'm'], ['GMl corrected', GMl, 'm'], ['Wetted surface', hs.wettedArea, 'm²'], ['TPC', hs.TPC, 't/cm'], ['MCT', hs.MCTcm, 't·m/cm'],
      ['Tendon pretension', hull.pretension / 1e6, 'MN'],
    ].map(([q, v, u]) => ({ q, v: round(v, 5), u }));
    return {
      kpis: [
        { key: 'disp', label: 'Displacement', value: round(hs.displacement / 1e3, 6), unit: 't', explain: 'Mass of water displaced at the floating position; equals the mass of the unit for a free floater (Archimedes), and exceeds it by the tendon pretension for a TLP.' },
        { key: 'volume', label: 'Displaced volume', value: round(hs.volume, 6), unit: 'm³', explain: 'Volume of the hull below the waterplane, integrated exactly over the closed surface mesh.' },
        { key: 'draft', label: 'Equilibrium draft', value: round(draft, 5), unit: 'm', explain: 'Depth of the keel origin below the still waterline at equilibrium. It is the input draft unless a mass or centre-of-gravity offset was given.' },
        { key: 'heelDeg', label: 'Heel (magnitude)', value: round(Math.abs(heel) / DEG, 4), unit: 'deg', limit: inp.heelMax, sense: 'max', explain: 'Static list about the longitudinal axis from a transverse offset of the centre of gravity. Zero for a symmetric loading.' },
        { key: 'trimDeg', label: 'Trim (magnitude)', value: round(Math.abs(trim) / DEG, 4), unit: 'deg', limit: inp.heelMax, sense: 'max', explain: 'Static inclination about the transverse axis from a longitudinal offset of the centre of gravity.' },
        { key: 'KB', label: 'Centre of buoyancy above keel (KB)', value: round(hs.KB, 5), unit: 'm', explain: 'Height of the centroid of the displaced volume; a deep, bottom-heavy displacement lowers it.' },
        { key: 'BMt', label: 'Transverse metacentric radius (BM)', value: round(hs.BMt, 5), unit: 'm', explain: 'Waterplane second moment divided by displaced volume. Widely spaced columns or a broad beam give a large BM; a spar has almost none and relies on KB above KG.' },
        { key: 'KMt', label: 'Transverse metacentre above keel (KM)', value: round(hs.KMt, 5), unit: 'm', explain: 'KB + BM: the highest position the centre of gravity may take before the unit becomes unstable.' },
        { key: 'GMt', label: 'Metacentric height, corrected (GM)', value: round(GMt, 5), unit: 'm', ...(isTlp ? {} : { limit: inp.gmMin, sense: 'min' }), explain: 'KM − KG − free-surface rise. The initial stability measure: larger is stiffer in roll (and shorter roll period); below the required value the unit is too tender.' },
        { key: 'GMl', label: 'Longitudinal metacentric height', value: round(GMl, 5), unit: 'm', explain: 'Same quantity for trim. Ship-shaped hulls have hundreds of metres; column-stabilised units have similar values in both directions.' },
        { key: 'GGt', label: 'Free-surface rise of G', value: round(fsc.GGt, 4), unit: 'm', explain: 'Virtual rise of the centre of gravity from liquid moving in slack tanks. Pressing tanks full or emptying them removes it.' },
        { key: 'Aw', label: 'Waterplane area', value: round(hs.waterplaneArea, 5), unit: 'm²', explain: 'Area cut by the waterline; sets the heave stiffness and tonnes per centimetre immersion.' },
        { key: 'TPC', label: 'Tonnes per centimetre immersion', value: round(hs.TPC, 4), unit: 't/cm', explain: 'Load that sinks the unit by one centimetre. Small for column-stabilised units: they are sensitive to deck-load changes.' },
        { key: 'C33', label: 'Heave restoring stiffness', value: round(K[2][2] / 1e3, 5), unit: 'kN/m', explain: 'ρ·g·Aw: vertical force per metre of sinkage, the spring of the heave natural period.' },
        { key: 'C44', label: 'Roll restoring stiffness', value: round((K[3][3] - fsc.GGt * mass * G) / 1e6, 5), unit: 'MN·m/rad', explain: 'Δ·g·GM including the free-surface loss: righting moment per radian of heel.' },
        { key: 'freeboard', label: 'Freeboard to hull top', value: round(freeboard, 4), unit: 'm', limit: inp.freeboardMin, sense: 'min', explain: 'Height of the top of the watertight hull above the waterline at the keel-origin draft.' },
        { key: 'reserve', label: 'Reserve buoyancy', value: round(100 * reserve, 4), unit: '%', explain: 'Enclosed hull volume above the waterline as a percentage of the displaced volume: what is left to absorb flooding or overload.' },
        { key: 'kgMax', label: 'Maximum allowable KG', value: round(kgMax, 5), unit: 'm', explain: 'Highest centre of gravity that still meets the required GM at this draft, with the present slack tanks.' },
      ],
      checks,
      plots: [
        hullPlot(mesh, draft, 'Hull surface and still-water pressure head'),
        { type: 'line', title: 'Hydrostatic curves', x: { label: 'Draft', unit: 'm' }, y: { label: 'KB, KM', unit: 'm' },
          series: [{ name: 'KMt', x: curve.map((c) => round(c.d)), y: curve.map((c) => round(Math.min(c.KMt, 20 * hull.depth))) }, { name: 'KB', x: curve.map((c) => round(c.d)), y: curve.map((c) => round(c.KB)) },
            { name: 'KG', x: [round(drafts[0]), round(drafts[drafts.length - 1])], y: [inp.KG, inp.KG], dash: true }],
          explain: 'Metacentre and buoyancy centre against draft, with the centre of gravity as a dashed line. The unit is initially stable wherever KM lies above KG; a semi-submersible shows a sharp drop of KM as the pontoons submerge (transit to operating draft).' },
        { type: 'line', title: 'Displacement and waterplane area versus draft', x: { label: 'Draft', unit: 'm' }, y: { label: 'Displacement', unit: 't' },
          series: [{ name: 'Displacement', x: curve.map((c) => round(c.d)), y: curve.map((c) => round(c.disp)) }, { name: 'Waterplane area × 10 (m²)', x: curve.map((c) => round(c.d)), y: curve.map((c) => round(c.Aw * 10)), dash: true }],
          explain: 'Deadweight scale: read the displacement at any draft. The slope is the waterplane area; a kink marks a change of waterplane (pontoon top, deck edge).' },
      ],
      tables: [
        { title: 'Hydrostatic particulars', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }], rows,
          explain: 'One row per hydrostatic quantity at the equilibrium floating position, in the form of a stability-booklet page.' },
        { title: 'Hydrostatic restoring matrix (heave, roll, pitch)', columns: [{ key: 'dof', label: 'DOF' }, { key: 'h', label: 'Heave', unit: 'N/m | N' }, { key: 'r', label: 'Roll', unit: 'N | N·m' }, { key: 'p', label: 'Pitch', unit: 'N | N·m' }],
          rows: [2, 3, 4].map((i) => ({ dof: ['Heave', 'Roll', 'Pitch'][i - 2], h: round(K[i][2], 5), r: round(K[i][3], 5), p: round(K[i][4], 5) })),
          explain: 'Linear restoring coefficients about the waterline origin, without free-surface correction; these are the C-matrix terms a motion analysis needs. Off-diagonal terms appear when the waterplane centroid is not at the origin.' },
        { title: 'Hydrostatic curves', columns: [{ key: 'd', label: 'Draft', unit: 'm' }, { key: 'disp', label: 'Displacement', unit: 't' }, { key: 'KB', label: 'KB', unit: 'm' }, { key: 'KMt', label: 'KMt', unit: 'm' }, { key: 'KMl', label: 'KMl', unit: 'm' }, { key: 'Aw', label: 'Aw', unit: 'm²' }, { key: 'TPC', label: 'TPC', unit: 't/cm' }],
          rows: curve.map((c) => Object.fromEntries(Object.entries(c).map(([k, v]) => [k, round(v, 5)]))),
          explain: 'Tabulated hydrostatics on even keel for a range of drafts, for loading-condition work away from the design draft.' },
      ],
      notes: [
        `${hull.imported ? 'Imported hull surface' : `Parametric ${hull.type} hull`}: ${mesh.indices.length / 3} triangles, displaced volume ${fmt(hs.volume, 5)} m³ at ${fmt(draft, 4)} m draft.`,
        solved ? `Equilibrium solved for mass ${fmt(mass / 1e3, 5)} t with G offset (${fmt(inp.lcgOff)}, ${fmt(inp.tcgOff)}) m: heel ${fmt(heel / DEG)}° (${heel >= 0 ? 'port' : 'starboard'} side up), trim ${fmt(trim / DEG)}° (${trim >= 0 ? 'bow' : 'stern'} down).` : 'The unit floats upright at the input draft; mass was derived from the displaced volume.',
        isTlp ? `Tension-leg platform: weight ${fmt(hull.mass / 1e3)} t is ${fmt(100 * (1 - inp.tendonPretensionFrac))} % of the buoyancy; the ${fmt(hull.pretension / 1e6)} MN difference is tendon pretension. The free-floating GM of ${fmt(GMt)} m is informative only.` : `Governing check: ${gov.name} at utilisation ${fmt(gov.uc)}.`,
      ],
      advice: [
        adv('engineering', isTlp ? 'info' : sev(gov.uc), isTlp ? `Buoyancy ${fmt(hs.buoyancy / 1e6)} MN, pretension ${fmt(hull.pretension / 1e6)} MN (${fmt(100 * inp.tendonPretensionFrac)} % of buoyancy).` : `GM is ${fmt(GMt)} m against ${fmt(inp.gmMin)} m required; free surfaces cost ${fmt(fsc.GGt)} m; freeboard ${fmt(freeboard)} m.`,
          isTlp ? 'Check minimum tendon tension in the extreme sea state and the free-floating stability for wet tow and installation, when the tendons are not connected.'
            : GMt < inp.gmMin ? `Lower KG to ≤ ${fmt(kgMax)} m (move ${fmt((inp.KG - kgMax) * mass / 1e3)} t·m of vertical moment down), press slack tanks full to recover ${fmt(fsc.GGt)} m, or widen the column spacing.`
              : `The centre of gravity can rise by ${fmt(kgMax - inp.KG)} m before the GM limit is reached; confirm KG with an inclining test and re-run the intact and damage stability modules, which normally govern the allowable KG.`,
          'IMO MODU Code 2009 ch. 3; IMO 2008 IS Code Part B ch. 3 and Annex 1 (inclining test); DNV-OS-C301 (2021).'),
        adv('economic', 'info', `Each metre of allowable-KG margin is worth about ${fmt(mass / 1e3 / Math.max(hull.depth - inp.KG + 10, 5))} t of extra deck payload at 10 m above the hull top; the present margin is ${fmt(kgMax - inp.KG)} m.`,
          'Treat KG margin as payload capacity: track it in the weight report from concept onwards, because recovering stability late means adding column diameter or ballast, both of which add steel and displacement.',
          'Weight control per ISO 19901-5 (2021); AACE RP 18R-97 class 4–3 estimates are driven by hull steel weight.',
          `Screening hull steel at 0.16 t per m³ of enclosed hull: ${steelCo2Text(steel)} A stability shortfall solved with fixed ballast instead of geometry carries no steel but permanently consumes payload.`),
      ],
      publish: { 'vessel.disp': hs.displacement, 'vessel.GM': GMt, 'quantities.steel_t': steel },
      quantities: { steel_t: round(steel, 4) },
    };
  },
  verify: [
    { // Box 100 × 20 m at 5 m draft, KG 4 m: ∇ = 10 000 m³, Δ = 10 250 t, KB = 2.5, BM = B²/(12T) = 400/60 = 6.6667,
      // KM = 9.1667, GM = 5.1667 m, TPC = 1.025·2000/100 = 20.5 t/cm.
      name: 'Box barge: KB = T/2, BM = B²/12T, Δ = ρLBT',
      source: 'Biran & López-Pulido, Ship Hydrostatics and Stability (2014) ch. 2 — closed form for a rectangular barge',
      inputs: { hullType: 'barge', L: 100, B: 20, Dhull: 12, draft: 5, KG: 4, tanks: [] },
      expect: [{ kpi: 'disp', value: 10250, tol: 1e-4 }, { kpi: 'KB', value: 2.5, tol: 1e-4 }, { kpi: 'BMt', value: 6.6667, tol: 1e-4 }, { kpi: 'GMt', value: 5.1667, tol: 1e-4 }, { kpi: 'TPC', value: 20.5, tol: 1e-4 }] },
    { // Vertical cylinder D = 30 m, T = 90 m: ∇ = π·15²·90 = 63 617 m³, KB = 45, BM = D²/(16T) = 900/1440 = 0.625 m.
      name: 'Spar: BM = D²/16T of a vertical circular cylinder',
      source: 'Closed form, I = πD⁴/64, ∇ = πD²T/4',
      inputs: { hullType: 'spar', Dspar: 30, Hspar: 105, draft: 90, KG: 38, tanks: [], nTheta: 96 },
      expect: [{ kpi: 'volume', value: 63617, tol: 1e-3 }, { kpi: 'BMt', value: 0.625, tol: 5e-3 }, { kpi: 'GMt', value: 7.625, tol: 2e-3 }] },
    { // Same box with G moved 0.10 m to port: wall-sided equilibrium tanφ·(GM + ½·BM·tan²φ) = 0.10 → tanφ = 0.019354,
      // φ = 1.1088°. Slack tank 20 × 10 m sea water: GG′ = l·b³/(12∇) = 20·1000/(12·10 000) = 0.1667 m → GM = 5.0 m
      // (the correction is reported, the mesh equilibrium itself is solved without it).
      name: 'Heel from a transverse shift of G and free-surface correction',
      source: 'Wall-sided formula and GG′ = ρᵢ·i/(ρ∇), Biran & López-Pulido (2014) ch. 5–6',
      inputs: { hullType: 'barge', L: 100, B: 20, Dhull: 12, draft: 5, KG: 4, tcgOff: 0.1, tanks: [{ name: 't', l: 20, b: 10, rho: 1025 }] },
      expect: [{ kpi: 'heelDeg', value: 1.1088, tol: 2e-3 }, { kpi: 'GGt', value: 0.16667, tol: 1e-3 }, { kpi: 'GMt', value: 5.0, tol: 1e-3 }, { kpi: 'draft', value: 5, tol: 1e-3 }] },
  ],
};
