// Structural buckling and local stress: mudmat plating under soil pressure, stiffeners, panel buckling, braces and skirts.
import { navierPlate, ritzPlate, largeDeflectionPlate, plateUltimateStrength, stiffenerWithPlate, stiffenerColumnBuckling, columnReduction, plateBucklingCoefficient } from '../../engine/structural/plates.js';
import { tubeSection } from '../../engine/structural/sections.js';
import { tubularMemberCheck } from '../../engine/structural/tubular.js';
import { need, run, fmt, round, adv, check, worst, steelCo2Text, CONFIRM } from './_common.js';

export default {
  id: 's04-buckling-local-stress', suite: 4, title: 'Structural buckling and local stress',
  summary: 'Plate bending of the mudmat plating under soil pressure, stiffener bending, plate and stiffened-panel buckling, column buckling of braces and buckling of skirts during penetration.',
  covers: [7], fidelity: 'analytical',
  standards: ['DNV-RP-C201 (2010) buckling of plated structures', 'EN 1993-1-5:2006 §4–5, §10', 'ISO 19902:2020 §13 (tubular members)', 'DNV-OS-C101 (2019) Ch.2 Sec.4', 'EN 1993-1-1:2005 §6.3.1'],
  models: ['Kirchhoff plate bending (Navier series, Rayleigh–Ritz)', 'von Kármán large-deflection plate (one-term Galerkin)', 'Plate buckling (elastic critical stress and EN 1993-1-5 effective width)', 'Stiffened-panel column buckling (Perry–Robertson)', 'Euler buckling', 'von Mises yielding', 'Hooke’s law', 'ISO 19902 tubular beam-column interaction'],
  limits: 'Idealised single panels with uniform pressure and uniform in-plane stress; simply supported or fully clamped edges bracket the real restraint. Weld details, cut-outs, stiffener tripping and global shell finite-element effects are not modelled; nonlinear finite-element collapse analysis is outside the module.',
  publishes: ['subsea.plateUc'],
  searchTerms: 'stiffened plate buckling mudmat plate bending soil pressure',
  inputs: [
    { key: 'mat', label: 'Steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Material', role: 'material' },
    { key: 'gammaM', label: 'Material factor γM', unit: '-', type: 'number', default: 1.15, min: 1, max: 2, group: 'Material', role: 'criteria', help: `DNV-OS-C101 LRFD: 1.15. ${CONFIRM}` },
    { key: 'corr', label: 'Corrosion allowance (deducted from plate)', unit: 'mm', type: 'number', default: 1, min: 0, max: 10, group: 'Material', role: 'ic', help: 'Thickness lost over the design life on the exposed face.' },
    { key: 'a', label: 'Stiffener span (panel length)', unit: 'mm', type: 'number', default: 2400, min: 200, max: 12000, group: 'Mudmat plating', role: 'geometry', design: true },
    { key: 'b', label: 'Stiffener spacing (panel width)', unit: 'mm', type: 'number', default: 700, min: 100, max: 4000, group: 'Mudmat plating', role: 'geometry', design: true },
    { key: 'tp', label: 'Plate thickness', unit: 'mm', type: 'number', default: 12, min: 4, max: 80, group: 'Mudmat plating', role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'edges', label: 'Panel edge restraint', type: 'select', default: 'clamped', group: 'Mudmat plating', role: 'bc', options: [{ value: 'clamped', label: 'Clamped (continuous plating over stiffeners)' }, { value: 'ss', label: 'Simply supported (conservative for deflection)' }] },
    { key: 'hw', label: 'Stiffener web height', unit: 'mm', type: 'number', default: 200, min: 40, max: 1500, group: 'Stiffeners', role: 'geometry', design: true },
    { key: 'tw', label: 'Stiffener web thickness', unit: 'mm', type: 'number', default: 10, min: 4, max: 60, group: 'Stiffeners', role: 'geometry' },
    { key: 'bf', label: 'Stiffener flange width', unit: 'mm', type: 'number', default: 90, min: 0, max: 600, group: 'Stiffeners', role: 'geometry' },
    { key: 'tf', label: 'Stiffener flange thickness', unit: 'mm', type: 'number', default: 12, min: 0, max: 60, group: 'Stiffeners', role: 'geometry' },
    { key: 'q', label: 'Design soil pressure on the plate', unit: 'kPa', type: 'number', default: 60, min: 0, max: 5000, group: 'Loads', role: 'bc', bind: 'subsea.contactPmax', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Peak factored contact pressure under the mat, including landing amplification where relevant.' },
    { key: 'sigX', label: 'In-plane compressive stress in the plating', unit: 'MPa', type: 'number', default: 60, min: 0, max: 700, group: 'Loads', role: 'bc', help: 'From global bending of the mat or frame action, acting along the stiffeners.' },
    { key: 'tau', label: 'In-plane shear stress in the plating', unit: 'MPa', type: 'number', default: 20, min: 0, max: 400, group: 'Loads', role: 'bc' },
    { key: 'Dm', label: 'Brace / post outer diameter', unit: 'mm', type: 'number', default: 406.4, min: 50, max: 3000, group: 'Brace or post', role: 'geometry', design: true },
    { key: 'tm', label: 'Brace / post wall thickness', unit: 'mm', type: 'number', default: 12.7, min: 3, max: 100, group: 'Brace or post', role: 'geometry', design: true },
    { key: 'Lm', label: 'Unbraced length', unit: 'm', type: 'number', default: 5, min: 0.3, max: 60, group: 'Brace or post', role: 'geometry' },
    { key: 'Km', label: 'Effective length factor K', unit: '-', type: 'number', default: 1.0, min: 0.5, max: 2.5, group: 'Brace or post', role: 'bc', help: '1.0 pinned ends; 0.7–0.8 for braces welded into stiff chords; 2.0 for a cantilever post.' },
    { key: 'Nm', label: 'Axial compression in the member', unit: 'kN', type: 'number', default: 600, min: 0, max: 1e5, group: 'Brace or post', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'Mm', label: 'Bending moment in the member', unit: 'kN·m', type: 'number', default: 60, min: 0, max: 1e5, group: 'Brace or post', role: 'bc' },
    { key: 'hs', label: 'Skirt height', unit: 'mm', type: 'number', default: 500, min: 0, max: 10000, group: 'Skirt', role: 'geometry' },
    { key: 'ts', label: 'Skirt thickness', unit: 'mm', type: 'number', default: 15, min: 4, max: 80, group: 'Skirt', role: 'geometry', design: true },
    { key: 'qs', label: 'Skirt tip line load during penetration', unit: 'kN/m', type: 'number', default: 40, min: 0, max: 20000, group: 'Skirt', role: 'bc', help: 'Penetration resistance divided by the total skirt length.' },
    { key: 'Ks', label: 'Skirt effective length factor', unit: '-', type: 'number', default: 2.0, min: 0.7, max: 2.5, group: 'Skirt', role: 'bc', help: '2.0 for a skirt free at the tip before soil support develops; about 0.7–1.0 once embedded.' },
    { key: 'nBasis', label: 'Ritz trial functions per direction', type: 'int', default: 8, min: 3, max: 14, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nBasis', label: 'Ritz functions', levels: [3, 6, 12], qoi: 'sigPlate' },
  compute(inp, ctx) {
    // SI: mm → m, kPa/MPa → Pa, kN → N.
    const m = ctx.material(inp.mat);
    need(m, `Steel '${inp.mat}' is not in the materials library.`);
    const E = m.E, nu = m.nu ?? 0.3, fy = m.fy, fd = fy / inp.gammaM;
    const a = inp.a / 1e3, b = inp.b / 1e3, t = (inp.tp - inp.corr) / 1e3;
    need(t > 0.002, `Net plate thickness after the corrosion allowance is ${fmt(t * 1e3)} mm: increase the plate or reduce the allowance.`);
    need(a >= b, 'The stiffener span must be at least the stiffener spacing (panel length ≥ width).');
    const q = inp.q * 1e3, sigX = inp.sigX * 1e6, tau = inp.tau * 1e6;
    // --- plate bending under lateral pressure
    const load = { type: 'uniform', q: Math.max(q, 1e-6) };
    let mMax, wMax;
    if (inp.edges === 'ss') {
      const r = run(() => navierPlate({ a, b, t, E, nu, load }), 'Plate bending');
      mMax = Math.max(Math.abs(r.mx), Math.abs(r.my)); wMax = Math.abs(r.w);
    } else {
      const r = run(() => ritzPlate({ a, b, t, E, nu, edges: { x0: 'C', x1: 'C', y0: 'C', y1: 'C' }, load, nBasis: inp.nBasis }), 'Plate bending');
      mMax = Math.max(Math.abs(r.mx), Math.abs(r.my), ...Object.values(r.edgeMoments).map(Math.abs)); wMax = r.wMax;
    }
    const sigPlate = 6 * mMax / (t * t);
    // Tension face of the plate bending across the stiffeners together with compression along them (opposite signs).
    const vm = Math.sqrt(sigPlate * sigPlate + sigX * sigX + sigPlate * sigX + 3 * tau * tau);
    const ld = run(() => largeDeflectionPlate({ a, b, t, E, nu, q: Math.max(q, 1e-6) }), 'Large-deflection plate');
    // --- stiffener with attached plate as a continuous beam between frames: M = q·s·L²/12
    const hw = inp.hw / 1e3, tw = inp.tw / 1e3, bf = inp.bf / 1e3, tf = inp.tf / 1e3;
    const st = run(() => stiffenerWithPlate({ be: b, tp: t, hw, tw, bf, tf }), 'Stiffener section');
    const Mst = q * b * a * a / 12;
    const sigSt = Mst / Math.min(st.Zplate, st.Zflange) + sigX;
    // --- plate and panel buckling under in-plane stress
    const pu = run(() => plateUltimateStrength({ a, b, t, E, nu, fy, sigmaX: sigX, tau }), 'Plate buckling');
    const kb = plateBucklingCoefficient({ aspect: Math.min(a / b, 20) });
    const sc = run(() => stiffenerColumnBuckling({ L: a, s: b, tp: t, hw, tw, bf, tf, E, nu, fy }), 'Stiffener column buckling');
    // --- tubular brace or post
    const D = inp.Dm / 1e3, tm = inp.tm / 1e3, Lm = inp.Lm, N = inp.Nm * 1e3, Mm = inp.Mm * 1e3;
    need(tm < D / 2, 'The member wall thickness must be smaller than its radius.');
    const sec = tubeSection({ D, t: tm });
    const Pcr = Math.PI ** 2 * E * sec.Iy / (inp.Km * Lm) ** 2;
    const lam = Math.sqrt(sec.A * fy / Pcr);
    const chi = columnReduction({ lambda: lam, curve: 'a' });
    const tub = run(() => tubularMemberCheck({ code: 'ISO', section: sec, fy, L: Lm, forces: { N: -N, My: Mm }, options: { Ky: inp.Km, Kz: inp.Km, E } }), 'Tubular member check');
    const ucTub = Number.isFinite(tub.governing.uc) ? Math.min(tub.governing.uc, 99) : 99;
    // --- skirt as a wide plate column under the tip reaction
    const hs = inp.hs / 1e3, ts = inp.ts / 1e3, qs = inp.qs * 1e3;
    const sigSk = qs / ts;
    const sigEsk = hs > 0 ? Math.PI ** 2 * E * ts * ts / (12 * (1 - nu * nu) * (inp.Ks * hs) ** 2) : Infinity;
    const chiSk = hs > 0 ? columnReduction({ lambda: Math.sqrt(fy / sigEsk), curve: 'c' }) : 1;
    const checks = [
      check('Plate bending under soil pressure', sigPlate / 1e6, fd / 1e6, 'MPa', 'DNV-OS-C101 Ch.2 Sec.4 (yield, plate bending); Timoshenko & Woinowsky-Krieger', 'Guards against yielding of the mat plating between stiffeners. If exceeded the plate dishes permanently and sheds load to the stiffeners.'),
      check('Plate equivalent (von Mises) stress', vm / 1e6, fd / 1e6, 'MPa', 'DNV-OS-C101 Ch.2 Sec.4 [2] von Mises criterion', 'Combines bending, in-plane and shear stress in the plate. If exceeded yielding spreads through the thickness.'),
      check('Stiffener bending with in-plane stress', sigSt / 1e6, fd / 1e6, 'MPa', 'DNV-RP-C201 §7 (stiffener as beam with effective plate flange)', 'Guards against yielding of the stiffener–plate beam spanning between frames under soil pressure.'),
      { name: 'Plate panel buckling (compression and shear)', demand: round(pu.utilisation * inp.gammaM), capacity: 1, unit: 'interaction', uc: Math.min(pu.utilisation * inp.gammaM, 99), ref: 'EN 1993-1-5 §4.4, §5, §10 Eq. (10.5)', explain: 'Guards against buckling of the plate between stiffeners under in-plane compression and shear. If exceeded the plate buckles and load redistributes to the stiffeners.' },
      check('Stiffened-panel column buckling', sigX / 1e6, sc.sigmaU / inp.gammaM / 1e6, 'MPa', 'DNV-RP-C201 §7 / EN 1993-1-5 §4.5.3 (Perry–Robertson)', 'Guards against the stiffener with its plate flange buckling as a column between frames — a sudden, non-ductile panel collapse.'),
      { name: `Tubular member (${tub.governing.name})`, demand: round(ucTub), capacity: 1, unit: 'interaction', uc: ucTub, ref: tub.governing.ref || 'ISO 19902 §13', explain: 'Beam-column check of the brace or post under axial compression and bending. If exceeded the member buckles or yields and the frame loses a load path.' },
      check('Skirt buckling during penetration', sigSk / 1e6, chiSk * fy / inp.gammaM / 1e6, 'MPa', 'EN 1993-1-1 §6.3.1 column curve c on a unit plate strip; DNV-RP-C212 §6 (skirt structural integrity)', 'Guards against the skirt plate buckling under tip resistance before the soil supports it laterally; a buckled skirt cannot be penetrated further.'),
    ];
    const w = worst(checks), ok = w.uc <= 1;
    const tReq = (inp.tp - inp.corr) * Math.sqrt(Math.max(checks[0].uc, checks[1].uc)) + inp.corr;
    const tonnesPerMm = 7.85 * 1e-3; // t per m² per mm of plate
    const ucs = checks.map((c) => round(c.uc));
    const tsw = [], sgw = [];
    for (let i = 0; i <= 12; i++) { const tt = (inp.tp * (0.6 + i * 0.08) - inp.corr) / 1e3; if (tt > 0.002) { tsw.push(round(tt * 1e3 + inp.corr)); sgw.push(sigPlate * (t / tt) ** 2 / 1e6); } }
    const lams = [], chis = [];
    for (let i = 0; i <= 30; i++) { const x = i * 0.1; lams.push(round(x)); chis.push(columnReduction({ lambda: x, curve: 'a' })); }
    const advice = [
      ok ? adv('engineering', 'ok', `All local checks pass; "${w.name}" governs at ${fmt(w.uc)}. Plate bending stress is ${fmt(sigPlate / 1e6)} MPa against ${fmt(fd / 1e6)} MPa.`,
        tReq < inp.tp - 0.5 ? `The plate could be reduced to about ${fmt(Math.ceil(tReq))} mm for bending; check the buckling utilisation (${fmt(checks[3].uc)}) and minimum handling thickness before doing so.` : 'The plate thickness is close to its optimum; confirm weld details and cut-outs with a local finite-element model before fabrication.',
        'DNV-RP-C201 (2010); DNV-OS-C101 (2019) Ch.2 Sec.4; EN 1993-1-5:2006',
        `Each millimetre of plate on a 140 m² mat is ${steelCo2Text(140 * tonnesPerMm)}`)
        : adv('engineering', 'fail', `"${w.name}" is exceeded at ${fmt(w.uc)}.`,
          w.name.startsWith('Plate') ? `Increase the plate to about ${fmt(Math.ceil(tReq))} mm, or close the stiffener spacing to ${fmt(inp.b / Math.sqrt(w.uc))} mm (bending stress scales with spacing squared).`
            : w.name.startsWith('Tubular') ? (N > 0.5 * Pcr ? `The member is too slender: Euler load ${fmt(Pcr / 1e3)} kN against ${fmt(N / 1e3)} kN applied. Increase the diameter to about ${fmt(inp.Dm * (2.5 * N / Pcr) ** (1 / 3))} mm (stiffness grows with D³) or halve the unbraced length.` : `Increase the wall to about ${fmt(inp.tm * w.uc * 1.05)} mm or shorten the unbraced length; Euler load is ${fmt(Pcr / 1e3)} kN against ${fmt(N / 1e3)} kN applied.`)
              : w.name.startsWith('Skirt') ? `Thicken the skirt to about ${fmt(inp.ts * Math.sqrt(w.uc) * 1.05)} mm or add vertical skirt stiffeners.` : `Deepen the stiffener web (section modulus needed ≈ ${fmt(Math.min(st.Zplate, st.Zflange) * w.uc * 1e6)} cm³) or add an intermediate frame to halve the span.`,
          'DNV-RP-C201 (2010); ISO 19902:2020 §13; EN 1993-1-5:2006',
          w.name.startsWith('Plate') || w.name.startsWith('Stiff') ? `Closer stiffener spacing is usually lighter than thicker plate: 1 mm of plate on a 140 m² mat is ${steelCo2Text(140 * tonnesPerMm)}` : `The member weighs ${fmt(sec.A * Lm * 7.85)} t; a larger diameter at the same wall adds stiffness far faster than mass. Per tonne of added steel: 1.9 t CO₂e (screening factor).`),
      adv('economic', ok ? 'info' : 'warn', `Plate thickness drives mat weight: ${fmt(inp.tp)} mm plate is ${fmt(inp.tp * tonnesPerMm * 1e3)} kg/m²; the bending-required thickness is ${fmt(tReq)} mm.`,
        ok ? 'Trade plate thickness against stiffener count in fabrication man-hours: thinner plate with closer stiffeners saves steel but adds weld length. Fix the choice with the fabricator before issuing drawings.' : 'Resolve the overstress on paper: a change after fabrication means cutting and re-welding plate in the yard, typically weeks on the critical path.',
        'AACE RP 18R-97 estimate classes; design-for-fabrication practice (weld length against plate tonnage)'),
    ];
    return {
      kpis: [
        { key: 'sigPlate', label: 'Plate bending stress', unit: 'MPa', value: sigPlate / 1e6, limit: fd / 1e6, sense: 'max', explain: 'Largest extreme-fibre bending stress in the panel (edge mid-point when clamped, centre when simply supported) under the soil pressure.' },
        { key: 'vm', label: 'Plate equivalent stress', unit: 'MPa', value: vm / 1e6, limit: fd / 1e6, sense: 'max', explain: 'von Mises combination of plate bending, in-plane compression and shear. Above the design yield the plate yields locally.' },
        { key: 'wPlate', label: 'Plate deflection', unit: 'mm', value: wMax * 1e3, explain: 'Largest lateral deflection of the panel by small-deflection theory. Beyond about half the thickness membrane action stiffens the plate.' },
        { key: 'wLarge', label: 'Plate deflection with membrane action', unit: 'mm', value: ld.w0 * 1e3, explain: 'Centre deflection from the one-term von Kármán solution (simply supported, edges held): the realistic value when deflection approaches the thickness.' },
        { key: 'sigSt', label: 'Stiffener stress', unit: 'MPa', value: sigSt / 1e6, limit: fd / 1e6, sense: 'max', explain: 'Bending stress of the stiffener with its plate flange spanning between frames plus the in-plane stress.' },
        { key: 'sigCr', label: 'Plate elastic critical stress', unit: 'MPa', value: pu.sigmaCr / 1e6, explain: 'Elastic buckling stress of the panel in compression. Well above yield means the plate is stocky and fully effective.' },
        { key: 'ucPanel', label: 'Plate buckling utilisation', unit: '-', value: pu.utilisation * inp.gammaM, limit: 1, sense: 'max', explain: 'Reduced-stress interaction of compression and shear against the ultimate plate strength, material factor included.' },
        { key: 'sigU', label: 'Stiffened-panel ultimate stress', unit: 'MPa', value: sc.sigmaU / 1e6, explain: 'Average compressive stress at which the stiffener with effective plate buckles as a column between frames.' },
        { key: 'Pcr', label: 'Euler buckling load of the member', unit: 'kN', value: Pcr / 1e3, explain: 'Elastic critical load π²EI/(KL)² of the brace or post. The applied compression should be a small fraction of it.' },
        { key: 'chi', label: 'Column reduction factor of the member', unit: '-', value: chi, explain: 'Buckling reduction factor on the squash load (EN 1993-1-1 curve a) at the member slenderness.' },
        { key: 'ucMember', label: 'Tubular member utilisation', unit: '-', value: ucTub, limit: 1, sense: 'max', explain: 'Governing ISO 19902 beam-column utilisation of the brace or post (capped at 99 beyond the Euler load).' },
        { key: 'tReq', label: 'Plate thickness required for bending', unit: 'mm', value: tReq, explain: 'Gross thickness (corrosion allowance included) at which the plate bending and equivalent stress just meet the design yield; stress scales with 1/t².' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Utilisation by check', categories: checks.map((c) => c.name), series: [{ name: 'Utilisation', y: ucs }], y: { label: 'Utilisation', unit: '-' },
          explain: 'All local strength and buckling utilisations side by side. Bars above 1.0 need strengthening; very low bars show where material can be saved.' },
        { type: 'line', title: 'Plate bending stress against thickness', x: { label: 'Gross plate thickness', unit: 'mm' }, y: { label: 'Bending stress', unit: 'MPa' },
          series: [{ name: 'Bending stress', x: tsw, y: sgw }, { name: 'Design yield', x: [tsw[0], tsw[tsw.length - 1]], y: [fd / 1e6, fd / 1e6], dash: true }],
          explain: 'Bending stress falls with the square of thickness. The crossing with the design-yield line gives the minimum plate.' },
        { type: 'line', title: 'Column curve and member slenderness', x: { label: 'Relative slenderness λ', unit: '-' }, y: { label: 'Reduction factor χ', unit: '-' },
          series: [{ name: 'Curve a', x: lams, y: chis }, { name: 'Brace / post', x: [round(lam)], y: [chi] }, { name: 'Skirt strip (curve c)', x: [round(Math.sqrt(fy / sigEsk))], y: [chiSk] }],
          explain: 'Where the brace and the skirt sit on the buckling curve. Points on the steep part (λ above about 1) are sensitive to effective length and imperfections.' },
      ],
      tables: [{ title: 'Section and buckling data', columns: [{ key: 'k', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
        rows: [['Net plate thickness', t * 1e3, 'mm'], ['Panel aspect a/b', a / b, '-'], ['Plate buckling coefficient k', kb.k, '-'], ['Plate slenderness λp', pu.lambdaP, '-'], ['Plate effectiveness ρ', pu.rho, '-'], ['Shear buckling stress', pu.tauCr / 1e6, 'MPa'],
          ['Stiffener + plate area', st.A * 1e4, 'cm²'], ['Stiffener + plate inertia', st.I * 1e8, 'cm⁴'], ['Section modulus at flange', st.Zflange * 1e6, 'cm³'], ['Stiffener moment qsL²/12', Mst / 1e3, 'kN·m'], ['Stiffener column slenderness', sc.lambda, '-'],
          ['Member slenderness λ', lam, '-'], ['Member area', sec.A * 1e4, 'cm²'], ['Skirt Euler stress', Number.isFinite(sigEsk) ? sigEsk / 1e6 : 0, 'MPa'], ['Skirt compressive stress', sigSk / 1e6, 'MPa']].map(([k, v, u]) => ({ k, v: round(v), u })),
        explain: 'Intermediate section properties and slenderness parameters behind each check, for comparison with a hand calculation or a code-check package.' }],
      notes: [`"${w.name}" governs at a utilisation of ${fmt(w.uc)}.`, wMax > t / 2 ? 'Plate deflection exceeds half the thickness: small-deflection bending stresses are conservative, see the membrane-action value.' : 'Plate deflection is within the small-deflection range.'],
      advice,
      publish: { 'subsea.plateUc': Math.max(checks[0].uc, checks[1].uc) },
      quantities: {},
    };
  },
  verify: [
    { name: 'Simply supported square plate, uniform load: Mmax = 0.0479·q·a²', source: 'Timoshenko & Woinowsky-Krieger, Theory of Plates and Shells, Table 8 (ν = 0.3): σ = 6·0.0479·10e3·1²/0.01² = 28.74 MPa; w = 0.00406·q·a⁴/D',
      inputs: { a: 1000, b: 1000, tp: 10, corr: 0, q: 10, edges: 'ss', sigX: 0, tau: 0, mat: 's355' },
      expect: [{ kpi: 'sigPlate', value: 28.74, tol: 0.005 }, { kpi: 'wPlate', value: 0.00406 * 10e3 / (210e9 * 1e-6 / (12 * 0.91)) * 1e3, tol: 0.005 }] },
    { name: 'Clamped square plate, uniform load: edge moment 0.0513·q·a²', source: 'Timoshenko & Woinowsky-Krieger, Table 35: σ = 6·0.0513·10e3/0.01² = 30.78 MPa',
      inputs: { a: 1000, b: 1000, tp: 10, corr: 0, q: 10, edges: 'clamped', sigX: 0, tau: 0, mat: 's355', nBasis: 12 },
      expect: [{ kpi: 'sigPlate', value: 30.78, tol: 0.02 }] },
    { name: 'Euler load of a pinned tube: Pcr = π²EI/L²', source: 'Euler (1744): D = 400 mm, t = 16 mm, L = 6 m, E = 210 GPa',
      inputs: { Dm: 400, tm: 16, Lm: 6, Km: 1, mat: 's355' },
      expect: [{ kpi: 'Pcr', value: Math.PI ** 2 * 210e9 * (Math.PI / 64) * (0.4 ** 4 - 0.368 ** 4) / 36 / 1e3, tol: 1e-4 }] },
  ],
};
