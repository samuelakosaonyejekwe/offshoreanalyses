// Module support, skid and grillage: distribution of the module weight to its support points (rigid module on
// elastic supports with centre-of-gravity offset and overturning), beam-grid analysis of the supporting grillage
// under those reactions, skid runner bending between supports, and support adequacy.
import { grillageBeamGrid } from '../../engine/structural/plates.js';
import { iSection } from '../../engine/structural/sections.js';
import { steelOf, need, run, fmt, round, adv, check, steelCo2Text, CONFIRM, G } from './_common.js';

const SEC = (p, g, d) => [
  { key: `${p}h`, label: 'Depth', unit: 'mm', type: 'number', default: d[0], min: 80, max: 4000, group: g, role: 'geometry', design: true },
  { key: `${p}b`, label: 'Flange width', unit: 'mm', type: 'number', default: d[1], min: 40, max: 1500, group: g, role: 'geometry' },
  { key: `${p}tw`, label: 'Web thickness', unit: 'mm', type: 'number', default: d[2], min: 4, max: 80, group: g, role: 'geometry' },
  { key: `${p}tf`, label: 'Flange thickness', unit: 'mm', type: 'number', default: d[3], min: 6, max: 120, group: g, role: 'geometry', design: true },
];

/** 3 × 3 symmetric solve by Cramer's rule; returns null when the matrix is singular (supports collinear). */
function solve3(A, b) {
  const det = (m) => m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  const d = det(A), scale = Math.abs(A[0][0] * A[1][1] * A[2][2]);
  if (!(Math.abs(d) > 1e-10 * scale)) return null;
  return [0, 1, 2].map((c) => det(A.map((row, i) => row.map((v, j) => (j === c ? b[i] : v)))) / d);
}

export default {
  id: 's02-skid-grillage', suite: 2, title: 'Module support, skid and grillage',
  summary: 'Support reactions of a module or skid on several support points, load distribution and deflection of the supporting grillage of girders and beams, skid runner stresses and adequacy of each support.',
  covers: [9], fidelity: 'numerical',
  standards: ['ISO 19901-3:2014 §9.3 (module supports), §10 (grillages for load-out and transport)', 'DNV-ST-N001 (2021) §5.6, §11 (load distribution, grillage and seafastening)', 'EN 1993-1-1:2005 §6.2 (cross-section resistance)', 'NORSOK N-004 (2022) §6', 'ISO 19901-5:2021 (weighing and centre of gravity)'],
  models: ['Grillage analysis (beam-grid direct stiffness)', 'Multi-point support stiffness analysis (rigid body on elastic supports)', 'Static equilibrium', 'Euler–Bernoulli beam theory', 'Saint-Venant torsion', 'Three-moment (continuous beam) approximation for the skid runner'],
  limits: 'The module is rigid on linear support springs (no lift-off: a support in tension is reported, not released), and the grillage is a rectangular beam grid simply supported or clamped along its four edges with the support reactions applied at the nearest girder–beam intersections. Module flexibility, shimming sequence, friction and the local web crippling under each support are not modelled; shear deformation of deep grillage beams is neglected.',
  searchTerms: 'module support grillage skid multi-point support reactions load distribution offshore module',
  inputs: [
    { key: 'W', label: 'Module weight', unit: 't', type: 'number', default: 320, min: 1, max: 30000, group: 'Module', role: 'load', uncertain: { dist: 'normal', cov: 0.05 }, help: 'Weighed or estimated weight including contingency.' },
    { key: 'xg', label: 'Centre of gravity x', unit: 'm', type: 'number', default: 6.6, min: -100, max: 200, group: 'Module', role: 'geometry', uncertain: { dist: 'normal', cov: 0.03 }, help: 'In the grillage plan coordinates (origin at a grillage corner).' },
    { key: 'yg', label: 'Centre of gravity y', unit: 'm', type: 'number', default: 4.3, min: -100, max: 200, group: 'Module', role: 'geometry' },
    { key: 'hg', label: 'Centre of gravity above the supports', unit: 'm', type: 'number', default: 4, min: 0, max: 60, group: 'Module', role: 'geometry' },
    { key: 'ax', label: 'Horizontal acceleration along x', unit: 'g', type: 'number', default: 0.05, min: 0, max: 1.5, group: 'Module', role: 'bc', help: 'Platform motion, wind-equivalent or transport acceleration; it adds an overturning couple W·a·h to the support reactions.' },
    { key: 'ay', label: 'Horizontal acceleration along y', unit: 'g', type: 'number', default: 0.08, min: 0, max: 1.5, group: 'Module', role: 'bc' },
    { key: 'gF', label: 'Load factor on the module weight', unit: '-', type: 'number', default: 1.3, min: 1, max: 2, group: 'Module', role: 'criteria', help: `Combined partial factor and weight contingency. ${CONFIRM}` },
    { key: 'sup', label: 'Support points', type: 'table', group: 'Supports', role: 'bc', help: 'Plan position of each module support and its vertical stiffness (bearing, shim pack and local structure in series). At least three, not in one line.',
      columns: [{ key: 'name', label: 'Support', type: 'text' }, { key: 'x', label: 'x', unit: 'm', type: 'number' }, { key: 'y', label: 'y', unit: 'm', type: 'number' }, { key: 'k', label: 'Stiffness', unit: 'kN/mm', type: 'number' }],
      default: [{ name: 'A1', x: 3, y: 2, k: 800 }, { name: 'A2', x: 6, y: 2, k: 800 }, { name: 'A3', x: 9, y: 2, k: 800 }, { name: 'B1', x: 3, y: 6, k: 800 }, { name: 'B2', x: 6, y: 6, k: 800 }, { name: 'B3', x: 9, y: 6, k: 800 }] },
    { key: 'Rcap', label: 'Design capacity of one support', unit: 'kN', type: 'number', default: 1500, min: 10, max: 1e5, group: 'Supports', role: 'criteria', help: 'Bearing pad, support stool or local web capacity per support point.' },
    { key: 'a', label: 'Grillage length (x)', unit: 'm', type: 'number', default: 12, min: 1, max: 80, group: 'Grillage', role: 'geometry' },
    { key: 'b', label: 'Grillage width (y)', unit: 'm', type: 'number', default: 8, min: 1, max: 60, group: 'Grillage', role: 'geometry' },
    { key: 'nG', label: 'Girders along x (inside the edges)', type: 'int', default: 3, min: 1, max: 12, group: 'Grillage', role: 'geometry', design: true },
    { key: 'nB', label: 'Beams along y (inside the edges)', type: 'int', default: 3, min: 1, max: 16, group: 'Grillage', role: 'geometry', design: true },
    { key: 'edge', label: 'Grillage edge support', type: 'select', default: 'simple', group: 'Grillage', role: 'bc', options: [{ value: 'simple', label: 'Simply supported on the deck girders' }, { value: 'clamped', label: 'Clamped (continuous over stiff girders)' }] },
    { key: 'q', label: 'Area load on the grillage', unit: 'kPa', type: 'number', default: 2, min: 0, max: 50, group: 'Grillage', role: 'load', help: 'Self-weight of the grillage, plating and any live load around the module.' },
    ...SEC('g', 'Grillage girders (along x)', [900, 350, 16, 30]),
    ...SEC('t', 'Grillage beams (along y)', [900, 350, 16, 30]),
    ...SEC('r', 'Skid runners', [500, 300, 14, 25]),
    { key: 'nRun', label: 'Skid runners', type: 'int', default: 2, min: 1, max: 8, group: 'Skid runners', role: 'geometry', help: 'Longitudinal base beams of the module that span between the support points.' },
    { key: 'mat', label: 'Steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Grillage', role: 'material' },
    { key: 'gM', label: 'Material factor γM', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.5, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'defRatio', label: 'Grillage deflection limit (span / value)', unit: '-', type: 'number', default: 300, min: 100, max: 1500, group: 'Acceptance criteria', role: 'criteria' },
  ],
  compute(inp, ctx) {
    // SI: t → kg, kN/mm → N/m, mm → m, kPa → Pa, g → m/s².
    const st = steelOf(ctx, inp.mat), W = inp.W * 1e3 * G * inp.gF;
    const sup = (inp.sup || []).filter((s) => s && Number.isFinite(Number(s.x)) && Number.isFinite(Number(s.y))).map((s, i) => ({ name: String(s.name || `S${i + 1}`), x: Number(s.x), y: Number(s.y), k: Number(s.k) * 1e6 }));
    need(sup.length >= 3, 'Enter at least three support points: fewer cannot hold the module against tipping (under-constrained).');
    need(sup.every((s) => s.k > 0), 'Every support needs a positive stiffness.');
    for (const s of sup) need(s.x >= 0 && s.x <= inp.a && s.y >= 0 && s.y <= inp.b, `Support ${s.name} at (${s.x}, ${s.y}) m lies outside the ${inp.a} × ${inp.b} m grillage.`);
    // Rigid module on springs: w(x, y) = w0 + θx·(y − yg) − θy·(x − xg); ΣF = W, ΣM about the centre of gravity = overturning couples.
    const dx = sup.map((s) => s.x - inp.xg), dy = sup.map((s) => s.y - inp.yg), k = sup.map((s) => s.k), S = (f) => sup.reduce((a, _, i) => a + f(i), 0);
    const A = [[S((i) => k[i]), S((i) => k[i] * dy[i]), S((i) => -k[i] * dx[i])], [S((i) => k[i] * dy[i]), S((i) => k[i] * dy[i] * dy[i]), S((i) => -k[i] * dx[i] * dy[i])], [S((i) => -k[i] * dx[i]), S((i) => -k[i] * dx[i] * dy[i]), S((i) => k[i] * dx[i] * dx[i])]];
    const Mx = W * inp.ay * inp.hg, My = W * inp.ax * inp.hg;   // couples from horizontal inertia at the centre of gravity
    const reactions = (mx, my) => { const u = solve3(A, [W, mx, -my]); need(u, 'The support points are in one line (or too few): the module is free to tip. Add a support off the line.'); return sup.map((_, i) => k[i] * (u[0] + u[1] * dy[i] - u[2] * dx[i])); };
    const R0 = reactions(0, 0), cases = [[Mx, My], [Mx, -My], [-Mx, My], [-Mx, -My]].map(([p, q]) => reactions(p, q));
    const Rmax = sup.map((_, i) => Math.max(R0[i], ...cases.map((c) => c[i]))), Rmin = sup.map((_, i) => Math.min(R0[i], ...cases.map((c) => c[i])));
    const iMax = Rmax.indexOf(Math.max(...Rmax)), gc = cases.reduce((p, c) => (Math.max(...c) > Math.max(...p) ? c : p)), settle = Math.max(...gc.map((r, i) => r / k[i]));
    const mk = (p, what) => { need(2 * inp[`${p}tf`] < inp[`${p}h`] && inp[`${p}tw`] < inp[`${p}b`], `${what}: not a valid I-shape, check the thicknesses.`); return run(() => iSection({ h: inp[`${p}h`] / 1e3, b: inp[`${p}b`] / 1e3, tw: inp[`${p}tw`] / 1e3, tf: inp[`${p}tf`] / 1e3, welded: true }), what); };
    const sg = mk('g', 'Grillage girder'), sb = mk('t', 'Grillage beam'), sr = mk('r', 'Skid runner'), Gs = st.G;
    // Reactions of the governing case go to the nearest girder–beam intersections.
    const nG = Math.round(inp.nG), nB = Math.round(inp.nB), sx = inp.a / (nB + 1), sy = inp.b / (nG + 1), acc = new Map();
    let offMax = 0;
    sup.forEach((s, n) => { const i = Math.min(Math.max(Math.round(s.x / sx), 1), nB), j = Math.min(Math.max(Math.round(s.y / sy), 1), nG), key = `${i},${j}`; offMax = Math.max(offMax, Math.hypot(s.x - i * sx, s.y - j * sy)); acc.set(key, (acc.get(key) || 0) + gc[n]); });
    const pointLoads = [...acc].map(([key, P]) => { const [i, j] = key.split(',').map(Number); return { i, j, P }; });
    const gr = run(() => grillageBeamGrid({ a: inp.a, b: inp.b, nGirders: nG, nBeams: nB, EIx: st.E * sg.Iy, EIy: st.E * sb.Iy, GJx: Gs * sg.J, GJy: Gs * sb.J, q: inp.q * 1e3, pointLoads, support: inp.edge }), 'Grillage analysis');
    need(Math.abs(Math.abs(gr.reactionSum) - gr.appliedSum) < 1e-3 * gr.appliedSum + 1, 'Equilibrium check failed in the grillage analysis.');
    const fyd = st.fy / inp.gM, sigG = gr.girderMomentMax / sg.Wely, sigB = gr.beamMomentMax / sb.Wely, span = Math.min(inp.a, inp.b), wLim = span / inp.defRatio;
    // Skid runner: continuous beam over the supports of one line carrying its share of the weight uniformly; M ≈ w s²/10
    // between supports and w c²/2 on the overhangs is not known without the module length, so the span moment is checked.
    const xs = [...new Set(sup.map((s) => round(s.x, 6)))].sort((p, q) => p - q), gap = xs.length > 1 ? Math.max(...xs.slice(1).map((v, i) => v - xs[i])) : inp.a, Lrun = xs.length > 1 ? xs[xs.length - 1] - xs[0] + gap : inp.a;
    const wRun = W / Math.max(1, Math.round(inp.nRun)) / Lrun, Mrun = wRun * gap * gap / 10, sigR = Mrun / sr.Wely;
    const checks = [
      check(`Most loaded support (${sup[iMax].name})`, Rmax[iMax] / 1e3, inp.Rcap, 'kN', 'ISO 19901-3 §9.3; DNV-ST-N001 §5.6 (support reactions with centre-of-gravity inaccuracy and overturning)', 'Largest factored reaction of any support point over the gravity and overturning cases against the capacity of one support. Above 1.0 the bearing, stool or local structure under that support is overloaded.'),
      check('Grillage girder bending', sigG / 1e6, fyd / 1e6, 'MPa', `Elastic bending M/W ≤ fy/γM, γM = ${inp.gM} (EN 1993-1-1 §6.2.5; NORSOK N-004 §6)`, 'Largest bending stress in the girders along x under the support reactions and area load. Above the limit the girder yields and the grillage sheds load to its neighbours in a way this elastic model does not follow.'),
      check('Grillage beam bending', sigB / 1e6, fyd / 1e6, 'MPa', `Elastic bending M/W ≤ fy/γM, γM = ${inp.gM}`, 'Largest bending stress in the transverse beams along y.'),
      check('Grillage deflection', gr.wMax * 1e3, wLim * 1e3, 'mm', `Shorter span / ${inp.defRatio} (project criterion)`, 'Largest deflection of the grillage under the module. Deflection differences between support points redistribute the module reactions and strain piping and structural tie-ins.'),
      check('Skid runner bending between supports', sigR / 1e6, fyd / 1e6, 'MPa', 'Continuous-beam span moment w s²/10 against fy/γM', 'Bending of the module base beam between adjacent support points under its share of the weight. An overstressed runner means supports are too far apart for the skid.'),
    ];
    const minR = Math.min(...Rmin);
    if (minR < 0) checks.push(check('Support uplift', -minR / 1e3, Math.max(...R0) / 1e3, 'kN', 'DNV-ST-N001 §5.6 (no uplift without tie-down)', 'A support goes into tension in an overturning case: without a tie-down the module lifts off and the remaining supports take more than calculated here.'));
    const gov = checks.reduce((p, c) => (c.uc > p.uc ? c : p)), ok = gov.uc <= 1;
    const grT = (nG * inp.a * sg.A + nB * inp.b * sb.A) * st.rho / 1e3, mean = W / sup.length, uneven = Math.max(...R0) / mean;
    const sigGov = Math.max(sigG, sigB), depthFac = Math.sqrt(Math.max(sigGov / fyd, gr.wMax / wLim));
    return {
      kpis: [
        { key: 'Rmax', label: 'Largest support reaction', unit: 'kN', value: Rmax[iMax] / 1e3, limit: inp.Rcap, sense: 'max', explain: `Largest factored reaction over gravity and overturning cases (support ${sup[iMax].name}).` },
        { key: 'Rmin', label: 'Smallest support reaction', unit: 'kN', value: minR / 1e3, explain: 'Least loaded support over all cases; a negative value is uplift needing a tie-down.' },
        { key: 'Rgrav', label: 'Largest reaction under gravity alone', unit: 'kN', value: Math.max(...R0) / 1e3, explain: 'Largest reaction from the factored weight and centre-of-gravity position only: compare with the weighing report.' },
        { key: 'uneven', label: 'Largest gravity share over the mean', unit: '-', value: uneven, explain: 'How unevenly the supports share the weight: 1.0 is even. High values come from the centre-of-gravity offset and stiffness differences.' },
        { key: 'settle', label: 'Largest support compression', unit: 'mm', value: settle * 1e3, explain: 'Compression of the most loaded support spring: the shim adjustment needed to level the module.' },
        { key: 'wMax', label: 'Largest grillage deflection', unit: 'mm', value: gr.wMax * 1e3, limit: wLim * 1e3, sense: 'max', explain: 'Deflection of the grillage at its most deflected intersection under the governing case.' },
        { key: 'sigG', label: 'Girder bending stress', unit: 'MPa', value: sigG / 1e6, limit: fyd / 1e6, sense: 'max', explain: 'Largest elastic bending stress in the grillage girders (along x).' },
        { key: 'sigB', label: 'Beam bending stress', unit: 'MPa', value: sigB / 1e6, limit: fyd / 1e6, sense: 'max', explain: 'Largest elastic bending stress in the grillage beams (along y).' },
        { key: 'sigRun', label: 'Skid runner bending stress', unit: 'MPa', value: sigR / 1e6, limit: fyd / 1e6, sense: 'max', explain: 'Bending stress in a skid runner between the two supports that are furthest apart.' },
        { key: 'steelT', label: 'Grillage steel', unit: 't', value: grT, explain: 'Weight of the grillage girders and beams (edge members excluded).' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Support reactions', categories: sup.map((s) => s.name), series: [{ name: 'Gravity', y: R0.map((r) => round(r / 1e3)) }, { name: 'Largest', y: Rmax.map((r) => round(r / 1e3)) }, { name: 'Smallest', y: Rmin.map((r) => round(r / 1e3)) }, { name: 'Capacity', y: sup.map(() => inp.Rcap) }], y: { label: 'Reaction', unit: 'kN' },
          explain: 'Reaction of each support under gravity and its range over the overturning cases, against the support capacity. Supports nearest the centre of gravity carry most; a negative smallest value is uplift.' },
        { type: 'heatmap', title: 'Grillage deflection', xs: Array.from({ length: nB + 2 }, (_, i) => round(i * sx)), ys: Array.from({ length: nG + 2 }, (_, j) => round(j * sy)), z: gr.w.map((row) => row.map((v) => round(v * 1e3))), zLabel: 'Deflection', zUnit: 'mm', x: { label: 'x', unit: 'm' }, y: { label: 'y', unit: 'm' },
          explain: 'Deflection of the grillage at every girder–beam intersection (edges are supports). The dish should be centred under the module; a lopsided dish shows load concentrating on one girder.' },
        { type: 'scatter', title: 'Support layout and centre of gravity', x: { label: 'x', unit: 'm' }, y: { label: 'y', unit: 'm' }, series: [{ name: 'Supports', x: sup.map((s) => s.x), y: sup.map((s) => s.y), mode: 'markers' }, { name: 'Centre of gravity', x: [inp.xg], y: [inp.yg], mode: 'markers' }, { name: 'Grillage outline', x: [0, inp.a, inp.a, 0, 0], y: [0, 0, inp.b, inp.b, 0] }],
          explain: 'Plan of the support points, the module centre of gravity and the grillage outline. The centre of gravity must lie well inside the polygon of the supports for all of them to stay in compression.' },
      ],
      tables: [{ title: 'Support reactions', columns: [{ key: 'name', label: 'Support' }, { key: 'x', label: 'x', unit: 'm' }, { key: 'y', label: 'y', unit: 'm' }, { key: 'R0', label: 'Gravity', unit: 'kN' }, { key: 'Rmax', label: 'Largest', unit: 'kN' }, { key: 'Rmin', label: 'Smallest', unit: 'kN' }, { key: 'share', label: 'Gravity share', unit: '%' }, { key: 'uc', label: 'Utilisation' }],
        rows: sup.map((s, i) => ({ name: s.name, x: s.x, y: s.y, R0: round(R0[i] / 1e3), Rmax: round(Rmax[i] / 1e3), Rmin: round(Rmin[i] / 1e3), share: round(R0[i] / W * 100, 3), uc: round(Rmax[i] / 1e3 / inp.Rcap) })),
        explain: 'Factored reaction of each support under gravity and the envelope over the four overturning directions. Use the gravity column to check the weighing results and the envelope to design stools, shims and the structure below.' }],
      notes: [`${gov.name} governs at ${fmt(gov.uc)}. Factored weight ${fmt(W / 1e3)} kN on ${sup.length} supports; gravity sum ${fmt(R0.reduce((p, c) => p + c, 0) / 1e3)} kN.`,
        `Reactions are applied to the grillage at the nearest intersections (largest offset ${fmt(offMax)} m)${offMax > 0.3 * Math.min(sx, sy) ? ' — move grillage members under the supports or add local beams' : ''}.`, `Grillage ${nG} girder(s) × ${nB} beam(s), ${inp.edge} edges: applied ${fmt(gr.appliedSum / 1e3)} kN, reactions ${fmt(Math.abs(gr.reactionSum) / 1e3)} kN.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `The supports and grillage carry the module: largest reaction ${fmt(Rmax[iMax] / 1e3)} kN at ${sup[iMax].name} (${fmt(Rmax[iMax] / 1e3 / inp.Rcap * 100)} % of capacity), grillage stress ${fmt(sigGov / 1e6)} MPa, deflection ${fmt(gr.wMax * 1e3)} mm.`,
          `Shim the supports to the calculated compression (up to ${fmt(settle * 1e3)} mm) at set-down, confirm the reactions by weighing (gravity shares in the table), and check web stiffeners under each support for ${fmt(Rmax[iMax] / 1e3)} kN.`, 'ISO 19901-3:2014 §9.3; DNV-ST-N001 (2021) §5.6; ISO 19901-5:2021',
          `Grillage steel ${fmt(grT)} t (${steelCo2Text(grT)}). ${gov.uc < 0.6 ? `With utilisation ${fmt(gov.uc)} the members could be about ${fmt((1 - Math.sqrt(Math.max(gov.uc, 0.2) / 0.9)) * 100, 2)} % shallower.` : 'The members are reasonably used.'}`)
          : adv('engineering', 'fail', `${gov.name} is not met (utilisation ${fmt(gov.uc)}): largest reaction ${fmt(Rmax[iMax] / 1e3)} kN, grillage stress ${fmt(sigGov / 1e6)} MPa against ${fmt(fyd / 1e6)} MPa, deflection ${fmt(gr.wMax * 1e3)} mm.`,
            gov.name.startsWith('Most loaded') ? `Add a support near ${sup[iMax].name}, soften it with a bearing pad so its neighbours take more, or raise its capacity to ${fmt(Rmax[iMax] / 1e3)} kN.` : gov.name.startsWith('Support uplift') ? `Tie the module down at the lifting support for ${fmt(-minR / 1e3)} kN or widen the support base.` : gov.name.startsWith('Skid') ? `Add an intermediate support (gap now ${fmt(gap)} m) or deepen the runner to about ${fmt(inp.rh * Math.sqrt(gov.uc), 3)} mm.` : `Deepen the grillage members by about ${fmt((depthFac - 1) * 100, 2)} %, add a girder under the heavy support line, or clamp the edges onto the deck girders.`,
            'ISO 19901-3:2014 §9.3, §10; DNV-ST-N001 (2021) §5.6, §11; EN 1993-1-1:2005 §6.2', `Depth is the efficient route: ${fmt((depthFac - 1) * 100, 2)} % more depth at the same area, whereas the same gain by area would add about ${steelCo2Text(grT * Math.max(depthFac ** 2 - 1, 0.05))}`),
        adv('economic', 'info', `The grillage is ${fmt(grT)} t of steel for a ${inp.W} t module (${fmt(grT / inp.W * 100)} % of the module weight) and the reactions depend on a centre of gravity known to a tolerance.`,
          'Weigh the module before set-down and re-run with the measured centre of gravity; designing supports for a centre-of-gravity envelope is cheaper than re-shimming or reinforcing after a failed weigh-off.', 'ISO 19901-5:2021 weighing practice; AACE RP 18R-97 (contingency by definition level)'),
      ],
      quantities: { steel_t: grT },
    };
  },
  verify: [
    { // Rigid module on four equal supports at (±3, ±2) m about the support centre, centre of gravity 0.6 m off in x:
      // R = W/4 + W·ex·x/Σx² = 0.25W + W·0.6·3/36 = 0.30 W on the near side. W = 100 t → 980.665 kN → Rmax = 294.2 kN.
      name: 'Rigid body on four equal springs with an eccentric centre of gravity: R = W/n + W e x / Σx²', source: 'Statics of a rigid body on elastic supports (hand calculation)',
      inputs: { W: 100, gF: 1, ax: 0, ay: 0, xg: 6.6, yg: 4, sup: [{ name: 'A', x: 3, y: 2, k: 500 }, { name: 'B', x: 9, y: 2, k: 500 }, { name: 'C', x: 3, y: 6, k: 500 }, { name: 'D', x: 9, y: 6, k: 500 }] },
      expect: [{ kpi: 'Rmax', value: 294.2, tol: 1e-3 }, { kpi: 'uneven', value: 1.2, tol: 1e-3 }] },
    { // Two crossed simply supported beams 400 × 180 × 10 × 16 of 10 m span loaded at their crossing by 100 kN:
      // each carries half, w = P L³/(96 EI) with I = 2.5399e-4 m⁴ → 1e5·1000/(96·210e9·2.5399e-4) = 19.53 mm.
      name: 'Crossed beams with a central load: w = P L³/(96 EI)', source: 'Timoshenko, Strength of Materials Part I (two-beam grillage; hand calculation)',
      inputs: { W: 10.1972, gF: 1, ax: 0, ay: 0, q: 0, a: 10, b: 10, nG: 1, nB: 1, xg: 5, yg: 5, gh: 400, gb: 180, gtw: 10, gtf: 16, th: 400, tb: 180, ttw: 10, ttf: 16, sup: [{ name: 'A', x: 4.9, y: 4.9, k: 500 }, { name: 'B', x: 5.1, y: 4.9, k: 500 }, { name: 'C', x: 4.9, y: 5.1, k: 500 }, { name: 'D', x: 5.1, y: 5.1, k: 500 }] },
      expect: [{ kpi: 'wMax', value: 19.53, tol: 2e-3 }] },
  ],
};
