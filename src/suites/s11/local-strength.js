// Local strength of a stiffened hull panel: plate and stiffener bending under lateral pressure, combined von Mises
// stress with the hull-girder stress, plate buckling, stiffener column buckling and tripping.
import { navierPlate, largeDeflectionPlate, plateSlenderness, plateUltimateStrength, stiffenerWithPlate, stiffenerColumnBuckling, stiffenerTripping } from '../../engine/structural/plates.js';
import { G, need, run, round, rnd, fmt, adv, check, sev, worst, linspace, steelCo2Text } from './_common.js';

function evaluate(inp, mat, tPlate) {
  const fy = mat.fy; const E = mat.E; const nu = mat.nu ?? 0.3;
  const tp = (tPlate - inp.tCorr) / 1e3; const s = inp.s / 1e3; const l = inp.l;
  need(tp > 0.003, 'Net plate thickness after the corrosion deduction is below 3 mm: increase the plate or reduce the deduction.');
  const st = { hw: inp.hw / 1e3, tw: Math.max(inp.tw - inp.tCorr, 3) / 1e3, bf: inp.bf / 1e3, tf: inp.bf > 0 ? Math.max(inp.tf - inp.tCorr, 3) / 1e3 : 0 };
  const p = inp.rhoW * G * inp.head + inp.pDyn * 1e3;
  const sx = inp.sigmaX * 1e6; const tau = inp.tau * 1e6;
  // plate between stiffeners as a clamped strip (long panel): edge moment p·s²/12 → σ = p·s²/(2·t²)
  const sigPlate = p * s * s / (2 * tp * tp);
  const nav = p > 0 ? navierPlate({ a: l, b: s, t: tp, E, nu, load: { type: 'uniform', q: p }, nTerms: 39 }) : { w: 0 };
  const ld = p > 0 ? largeDeflectionPlate({ a: l, b: s, t: tp, E, nu, q: p }) : { w0: 0, ratio: 0 };
  const sec = stiffenerWithPlate({ be: s, tp, ...st });
  const Mst = p * s * l * l / 12; // stiffener clamped at the web frames
  const sigStiff = Mst / sec.Zflange;
  const tauStiff = (p * s * l / 2) / (st.hw * st.tw);
  // flange at the frame: bending compression adds to the axial hull-girder stress
  const vm = Math.sqrt((sx + sigStiff) ** 2 + 3 * tau * tau);
  const sl = plateSlenderness({ b: s, t: tp, E, nu, fy, k: 4 });
  const pu = plateUltimateStrength({ a: l, b: s, t: tp, E, nu, fy, sigmaX: Math.max(sx, 0), tau: Math.abs(tau) });
  const col = stiffenerColumnBuckling({ L: l, s, tp, ...st, E, nu, fy });
  const trip = stiffenerTripping({ L: l, ...st, E, nu, fy });
  return { fy, E, tp, s, l, p, sx, tau, sigPlate, nav, ld, sec, sigStiff, tauStiff, vm, sl, pu, col, trip, st };
}

export default {
  id: 's11-local-strength',
  suite: 11,
  title: 'Local stress and buckling of stiffened panels',
  summary: 'Checks a stiffened hull panel under lateral sea or tank pressure combined with hull-girder stress: plate and stiffener bending, von Mises yield, plate buckling, stiffener column buckling and tripping, with the thickness needed to pass.',
  covers: [8],
  fidelity: 'analytical',
  standards: ['DNV-RP-C201 (2010) Part 1 (buckling strength of plated structures)', 'DNV-OS-C101 (2021) Ch.2 Sec.4 (LRFD, γM = 1.15)', 'EN 1993-1-5 (2006) §4.4, §4.5, §5, §10', 'IACS CSR-H (2022) Pt 1 Ch 6 and Ch 8'],
  models: ['Thin plate theory (Navier series, clamped strip)', 'Large-deflection plate response (membrane action)', 'Elastic plate buckling σcr = k·π²E/(12(1−ν²))·(t/b)²', 'Effective-width ultimate plate strength with shear interaction', 'Perry–Robertson stiffener/plate column buckling', 'Stiffener tripping (torsional buckling) with Johnson–Ostenfeld correction', 'von Mises yield criterion'],
  limits: 'One panel between web frames with identical longitudinals; stresses are nominal (no hot-spot or bracket-toe concentration). The combined-load buckling check uses the reduced-stress interaction of compression and shear; transverse compression and lateral-pressure interaction on buckling are not included. Nonlinear finite-element collapse analysis and welding-induced imperfections beyond the code curves are outside the module.',
  searchTerms: 'stiffened panel buckling ultimate strength hull plating lateral pressure',
  inputs: [
    { key: 's', label: 'Stiffener spacing', unit: 'mm', type: 'number', default: 800, min: 300, max: 1500, group: 'Panel geometry', role: 'geometry', design: true, help: 'Spacing of the longitudinals: the short edge of the plate field.' },
    { key: 'l', label: 'Stiffener span (web-frame spacing)', unit: 'm', type: 'number', default: 3.2, min: 1, max: 10, group: 'Panel geometry', role: 'geometry', design: true, help: 'Length of the plate field and span of the stiffeners.' },
    { key: 'tp', label: 'Plate thickness (as built)', unit: 'mm', type: 'number', default: 18, min: 6, max: 80, group: 'Panel geometry', role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.02 }, help: 'Gross plate thickness.' },
    { key: 'hw', label: 'Stiffener web height', unit: 'mm', type: 'number', default: 400, min: 80, max: 1000, group: 'Panel geometry', role: 'geometry', design: true, help: 'Web of the tee, angle or flat bar.' },
    { key: 'tw', label: 'Stiffener web thickness', unit: 'mm', type: 'number', default: 12, min: 5, max: 40, group: 'Panel geometry', role: 'geometry', help: 'Web thickness.' },
    { key: 'bf', label: 'Stiffener flange width (0 = flat bar)', unit: 'mm', type: 'number', default: 150, min: 0, max: 400, group: 'Panel geometry', role: 'geometry', help: 'Flange width.' },
    { key: 'tf', label: 'Stiffener flange thickness', unit: 'mm', type: 'number', default: 18, min: 0, max: 50, group: 'Panel geometry', role: 'geometry', help: 'Flange thickness.' },
    { key: 'tCorr', label: 'Corrosion deduction', unit: 'mm', type: 'number', default: 1.5, min: 0, max: 10, group: 'Panel geometry', role: 'ic', help: 'Thickness lost or reserved for corrosion, deducted from plate, web and flange for the checks (pre-existing wastage as an initial condition).' },
    { key: 'steel', label: 'Steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Material', role: 'material', help: 'Grade of plate and stiffeners.' },
    { key: 'head', label: 'Static pressure head', unit: 'm', type: 'number', default: 14, min: 0, max: 300, group: 'Loads', role: 'bc', help: 'Height of sea water (or tank liquid) above the panel in the design condition.' },
    { key: 'rhoW', label: 'Liquid density', unit: 'kg/m³', type: 'number', default: 1025, min: 500, max: 2000, group: 'Loads', role: 'load', bind: 'site.rhoW', help: 'Density of the liquid giving the static head.' },
    { key: 'pDyn', label: 'Dynamic pressure', unit: 'kPa', type: 'number', default: 35, min: 0, max: 2000, group: 'Loads', role: 'bc', uncertain: { dist: 'normal', cov: 0.15 }, help: 'Wave-induced, sloshing or slamming pressure added to the static head.' },
    { key: 'sigmaX', label: 'In-plane compressive stress along the stiffeners', unit: 'MPa', type: 'number', default: 110, min: -300, max: 600, group: 'Loads', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Hull-girder or pontoon bending stress at the panel (compression positive), from the global strength module.' },
    { key: 'tau', label: 'In-plane shear stress', unit: 'MPa', type: 'number', default: 30, min: 0, max: 300, group: 'Loads', role: 'bc', help: 'Hull-girder shear stress at the panel.' },
    { key: 'gammaM', label: 'Material factor γM', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.5, group: 'Acceptance criteria', role: 'criteria', help: 'Resistance factor on yield and buckling strength (DNV-OS-C101 LRFD: 1.15). The loads entered are taken as design (factored) values.' },
  ],
  compute(inp, ctx) {
    const mat = ctx.material(inp.steel);
    need(mat && mat.fy > 0, 'Select a steel from the materials library.');
    const r = run(() => evaluate(inp, mat, inp.tp), 'Panel evaluation');
    const fyd = r.fy / inp.gammaM;
    const mk = (e) => [
      check('Plate bending under lateral pressure (yield)', e.sigPlate / 1e6, fyd / 1e6, 'MPa', 'DNV-OS-C101 (2021) Ch.2 Sec.4 [5] (plating); clamped strip p·s²/(2t²)', 'Guards against permanent set of the plating between stiffeners; beyond first yield the plate dishes and relies on membrane action.'),
      check('Stiffener bending + axial, von Mises at the flange', e.vm / 1e6, fyd / 1e6, 'MPa', 'DNV-OS-C101 (2021) Ch.2 Sec.4 [6] (stiffeners); von Mises criterion', 'Guards against yielding of the stiffener flange at the web frame, where pressure bending adds to the hull-girder compression.'),
      check('Stiffener web shear at the support', e.tauStiff / 1e6, fyd / Math.sqrt(3) / 1e6, 'MPa', 'DNV-OS-C101 (2021) Ch.2 Sec.4 [6]; τ ≤ fy/(√3·γM)', 'Guards the stiffener web against shear yielding at the frame connection under the pressure load.'),
      check('Plate buckling (compression + shear)', e.pu.utilisation * inp.gammaM, 1, '-', 'EN 1993-1-5 (2006) §4.4, §5, §10 Eq. (10.5); DNV-RP-C201 Part 1 §6', 'Guards against buckling of the plate field between stiffeners; a buckled plate sheds load to the stiffeners and the panel loses stiffness.'),
      check('Stiffener column buckling', Math.max(e.sx, 0) / 1e6, e.col.sigmaU / inp.gammaM / 1e6, 'MPa', 'DNV-RP-C201 (2010) Part 1 §7; EN 1993-1-5 §4.5.3', 'Guards against the stiffener with its effective plate buckling as a column between web frames — the collapse mode of a compressed deck or bottom panel.'),
      check('Stiffener tripping', Math.max(e.sx, 0) / 1e6, (e.trip.sigmaTinelastic ?? e.trip.sigmaT) / inp.gammaM / 1e6, 'MPa', 'DNV-RP-C201 (2010) Part 1 §7.5 (torsional buckling)', 'Guards against sideways twisting of the stiffener about its attachment line, which removes its support to the plate suddenly.'),
    ];
    const checks = mk(r); const gov = worst(checks);
    // plate thickness needed for every plate-related check to pass (bisection on the gross thickness)
    const ucAt = (t) => { try { const c = mk(evaluate(inp, mat, t)); return Math.max(c[0].uc, c[3].uc, c[4].uc); } catch { return 99; } };
    let lo = inp.tCorr + 4; let hi = 120; let tReq = hi;
    if (ucAt(hi) <= 1) { for (let i = 0; i < 40; i++) { const m = 0.5 * (lo + hi); if (ucAt(m) <= 1) hi = m; else lo = m; } tReq = hi; }
    const As = r.st.hw * r.st.tw + r.st.bf * r.st.tf;
    const massM2 = (inp.tp / 1e3 + (inp.hw * inp.tw + inp.bf * inp.tf) * 1e-6 / r.s) * 7850; // kg/m², gross
    const dMass = (tReq - inp.tp) * 7.85; // kg/m² per panel area
    const ts = linspace(Math.max(inp.tCorr + 6, 8), Math.max(inp.tp * 1.8, 30), 25);
    const sweep = ts.map((t) => { try { return mk(evaluate(inp, mat, t)); } catch { return null; } });
    return {
      kpis: [
        { key: 'pressure', label: 'Design lateral pressure', value: round(r.p / 1e3, 5), unit: 'kPa', explain: 'Static head plus dynamic pressure acting on the panel.' },
        { key: 'sigPlate', label: 'Plate bending stress', value: round(r.sigPlate / 1e6, 5), unit: 'MPa', limit: round(fyd / 1e6, 4), sense: 'max', explain: 'Surface bending stress at the stiffener line of a clamped plate strip, p·s²/(2t²), on the net thickness.' },
        { key: 'sigStiff', label: 'Stiffener bending stress', value: round(r.sigStiff / 1e6, 5), unit: 'MPa', explain: 'Bending stress at the stiffener flange at the web frame, p·s·l²/(12·Z), with the full plate flange.' },
        { key: 'vonMises', label: 'von Mises stress at the stiffener flange', value: round(r.vm / 1e6, 5), unit: 'MPa', limit: round(fyd / 1e6, 4), sense: 'max', explain: 'Combined axial, bending and shear stress; must stay below the design yield stress fy/γM.' },
        { key: 'deflection', label: 'Plate deflection (large-deflection theory)', value: round(r.ld.w0 * 1e3, 4), unit: 'mm', explain: 'Centre deflection of the plate field including membrane stiffening (simply supported, immovable edges). Above about half the thickness, linear plate theory over-predicts.' },
        { key: 'sigmaCr', label: 'Elastic plate buckling stress', value: round(r.sl.sigmaCr / 1e6, 5), unit: 'MPa', explain: 'Critical stress of the simply supported plate field in compression, 4·π²E/(12(1−ν²))·(t/s)². Above yield the plate is stocky and yields before it buckles.' },
        { key: 'beta', label: 'Plate slenderness β', value: round(r.sl.beta, 4), unit: '-', explain: '(s/t)·√(fy/E). Below about 1.0 the plate is fully effective; hull plating is typically 1.5–2.5.' },
        { key: 'sigmaUplate', label: 'Plate ultimate compressive strength', value: round(r.pu.sigmaU / 1e6, 5), unit: 'MPa', explain: 'Average stress the plate field can carry after buckling, from the effective-width method.' },
        { key: 'sigmaUcol', label: 'Stiffener column ultimate stress', value: round(r.col.sigmaU / 1e6, 5), unit: 'MPa', explain: 'Average compressive stress at which the stiffener with its effective plate buckles between frames.' },
        { key: 'sigmaTrip', label: 'Tripping stress', value: round((r.trip.sigmaTinelastic ?? r.trip.sigmaT) / 1e6, 5), unit: 'MPa', explain: 'Stress at which the stiffener twists sideways, with plasticity correction. Wide flanges and short spans raise it.' },
        { key: 'ucMax', label: 'Governing utilisation', value: round(gov.uc, 4), unit: '-', limit: 1, sense: 'max', explain: 'Largest demand-to-capacity ratio among the local checks.' },
        { key: 'tReq', label: 'Plate thickness to pass the plate and column checks', value: round(tReq, 4), unit: 'mm', explain: 'Smallest gross plate thickness at which plate bending, plate buckling and column buckling all pass with the present stiffeners.' },
        { key: 'panelMass', label: 'Panel steel mass', value: round(massM2, 4), unit: 'kg/m²', explain: 'Plate plus stiffeners per square metre of panel: the number to multiply by panel area for weight and cost.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Utilisation against plate thickness', x: { label: 'Gross plate thickness', unit: 'mm' }, y: { label: 'Utilisation', unit: '-' },
          series: [[0, 'Plate bending'], [3, 'Plate buckling'], [4, 'Column buckling']].map(([j, name]) => ({ name, x: rnd(ts.filter((_, i) => sweep[i])), y: rnd(sweep.filter(Boolean).map((c) => Math.min(c[j].uc, 3))) }))
            .concat([{ name: 'Limit', x: [round(ts[0]), round(ts[ts.length - 1])], y: [1, 1], dash: true }, { name: 'As built', x: [inp.tp, inp.tp], y: [0, Math.min(gov.uc * 1.2 + 0.2, 3)], dash: true }]),
          explain: 'How the plate-related utilisations fall as the plate gets thicker. The thickness where the highest curve crosses 1.0 is the required thickness; the vertical line is the present plate.' },
        { type: 'bar', title: 'Demand and capacity of the local checks', categories: checks.map((c) => c.name.split(' (')[0]), y: { label: 'Utilisation', unit: '-' }, series: [{ name: 'Utilisation', y: checks.map((c) => round(c.uc, 4)) }],
          explain: 'All local checks side by side. The tallest bar shows the failure mode that governs the scantlings.' },
      ],
      tables: [
        { title: 'Section and strength data', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['Net plate thickness', r.tp * 1e3, 'mm'], ['Yield stress', r.fy / 1e6, 'MPa'], ['Stiffener + plate area', r.sec.A * 1e4, 'cm²'], ['Stiffener area', As * 1e4, 'cm²'], ['Second moment with plate', r.sec.I * 1e8, 'cm⁴'], ['Section modulus at flange', r.sec.Zflange * 1e6, 'cm³'],
            ['Plate effective width ratio', r.pu.rho, '-'], ['Column slenderness', r.col.lambda, '-'], ['Column reduction χ', r.col.chi, '-'], ['Shear buckling strength', r.pu.tauU / 1e6, 'MPa'], ['Linear plate deflection (simply supported)', r.nav.w * 1e3, 'mm'], ['Deflection / thickness', r.ld.ratio, '-']].map(([q, v, u]) => ({ q, v: round(v, 5), u })),
          explain: 'Intermediate section properties and strengths behind the checks, for hand verification.' },
      ],
      notes: [
        `Net plate ${fmt(r.tp * 1e3)} mm on ${fmt(inp.s)} mm spacing: β = ${fmt(r.sl.beta)}, elastic buckling ${fmt(r.sl.sigmaCr / 1e6)} MPa, ultimate ${fmt(r.pu.sigmaU / 1e6)} MPa.`,
        `Pressure ${fmt(r.p / 1e3)} kPa with ${fmt(inp.sigmaX)} MPa compression and ${fmt(inp.tau)} MPa shear. Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at utilisation ${fmt(gov.uc)}; the plate-related checks pass from ${fmt(tReq)} mm gross plate (now ${fmt(inp.tp)} mm).`,
          gov.uc > 1 ? (gov.name.startsWith('Stiffener bending') || gov.name.startsWith('Stiffener web') || gov.name.startsWith('Stiffener tripping') ? `Strengthen the stiffener rather than the plate: a deeper web or wider flange raises the section modulus (now ${fmt(r.sec.Zflange * 1e6)} cm³) and tripping stress, or halve the span with an extra web frame (bending falls with l²).` : `Increase the plate to ${fmt(Math.ceil(tReq * 2) / 2)} mm, or reduce the stiffener spacing: buckling strength rises with (t/s)².`)
            : `Reserve is ${fmt(100 * (1 - gov.uc))} %. Confirm the corrosion deduction against the coating and inspection plan, and check the bracket toes and cut-outs for fatigue, which nominal-stress checks do not cover.`,
          'DNV-RP-C201 (2010) Part 1; DNV-OS-C101 (2021) Ch.2 Sec.4; EN 1993-1-5 (2006).'),
        adv('economic', 'info', `The panel weighs ${fmt(massM2)} kg/m². Moving to the required plate changes it by ${fmt(dMass)} kg/m² (${fmt(100 * dMass / massM2)} %).`,
          'Compare thicker plate with closer stiffener spacing: closer spacing saves plate weight but adds weld length and fabrication hours, which usually cost more per tonne saved than the steel itself.',
          'Structural optimisation practice (weight versus fabrication cost); AACE RP 18R-97 unit-rate estimating.',
          `Per 1000 m² of panel the thickness change is ${fmt(dMass)} t of steel: ${steelCo2Text(dMass)} A higher grade lowers mass only where yield, not buckling, governs.`),
      ],
    };
  },
  verify: [
    { // s = 800 mm, t = 16 mm net, E = 210 GPa, ν = 0.3: σE = π²·210e9/(12·0.91)·(0.016/0.8)² = 75.92 MPa, σcr = 4·σE = 303.7 MPa.
      // Clamped strip under 100 kPa: σ = p·s²/(2t²) = 1e5·0.64/(2·2.56e-4) = 125.0 MPa.
      name: 'Elastic plate buckling stress and clamped-strip bending stress',
      source: 'Timoshenko & Gere, Theory of Elastic Stability §9.2 (k = 4); beam formula M = p·s²/12 — hand calculation',
      inputs: { s: 800, l: 4, tp: 16, tCorr: 0, steel: 's355', head: 0, pDyn: 100, sigmaX: 100, tau: 0 },
      expect: [{ kpi: 'sigmaCr', value: 303.7, tol: 0.002 }, { kpi: 'sigPlate', value: 125.0, tol: 1e-4 }, { kpi: 'pressure', value: 100, tol: 1e-6 }] },
    { // Flat-bar-free check of the stiffener: tee 400×12 + 120×18 on 800×16 plate, l = 4 m, p = 100 kPa, no axial stress:
      // M = p·s·l²/12 = 1e5·0.8·16/12 = 106 667 N·m. Section by hand: areas 12800, 4800, 2160 mm² at z = 0, 208, 417 mm from the
      // plate mid-plane → zNA = 96.11 mm, I = 12800·96.11² + (12·400³/12 + 4800·111.89²) + (120·18³/12 + 2160·320.89²) + 800·16³/12
      //   = 1.1824e8 + 6.4000e7 + 6.0093e7 + 5.8e4 + 2.2242e8 + 2.73e5 = 4.6508e8 mm⁴; flange tip at 426 − 96.11 = 329.89 mm →
      // Z = 1.4098e6 mm³, σ = 106 667e3/1.4098e6 = 75.66 MPa.
      name: 'Stiffener bending stress p·s·l²/(12·Z) with a hand-calculated section modulus',
      source: 'Parallel-axis theorem and beam theory — hand calculation',
      inputs: { s: 800, l: 4, tp: 16, tCorr: 0, hw: 400, tw: 12, bf: 120, tf: 18, steel: 's355', head: 0, pDyn: 100, sigmaX: 0, tau: 0 },
      expect: [{ kpi: 'sigStiff', value: 75.66, tol: 0.005 }] },
  ],
};
