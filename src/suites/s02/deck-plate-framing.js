// Component design of deck plating and framing by classical theory: plate panel (Kirchhoff series and Mindlin
// finite elements), wheel/patch load, stiffener with effective plate flange, supporting beam and column.
import { navierPlate, ritzPlate, deckPatchLoad, stiffenerWithPlate, largeDeflectionPlate, plateRigidity } from '../../engine/structural/plates.js';
import { rectangleMesh, solvePlate } from '../../engine/structural/fem2d.js';
import { iSection, tubeSection } from '../../engine/structural/sections.js';
import { ec3MemberCheck } from '../../engine/structural/steelcheck.js';
import { tubularMemberCheck } from '../../engine/structural/tubular.js';
import { steelOf, need, run, fmt, round, adv, check, meshPlot, steelCo2Text, CONFIRM, G } from './_common.js';

export default {
  id: 's02-deck-plate-framing', suite: 2, title: 'Deck plate, stiffener, beam and column design',
  summary: 'Stresses, deflections and utilisations of a deck plate panel, its stiffeners, the supporting beam and a column under area and patch loads.',
  covers: [3], fidelity: 'analytical',
  standards: ['ISO 19901-3:2014 §7.6, §9 (deck plating and framing)', 'DNV-RP-C201 (2010) §5–§7 (plates and stiffeners under lateral load)', 'EN 1993-1-1:2005 §6.2–6.3', 'EN 1993-1-5:2006 §3 (effective width)', 'ISO 19902:2020 §13 (tubular columns)'],
  models: ['Kirchhoff–Love plate theory (Navier series, Rayleigh–Ritz)', 'Mindlin–Reissner plate theory (MITC4 finite elements)', 'Large-deflection plate analysis (von Kármán, one-term Galerkin)', 'Euler–Bernoulli beam theory', 'Shear deformation relationships', 'Slope-deflection / fixed-end moments for continuous members', 'Elastic beam-column analysis'],
  limits: 'Each component is checked on its own with idealised end conditions (panel edges all simply supported or all clamped, stiffener and beam simply supported or fully continuous). Elastic stresses only: yield-line and membrane reserve of the plate is reported but not used in the checks. Frame action is covered by the in-place module.',
  searchTerms: 'deck plate stiffener design lateral pressure wheel load effective width offshore topsides',
  inputs: [
    { key: 's', label: 'Stiffener spacing (panel width)', unit: 'mm', type: 'number', default: 700, min: 150, max: 3000, group: 'Plate panel', role: 'geometry', design: true },
    { key: 'Ls', label: 'Stiffener span (panel length)', unit: 'm', type: 'number', default: 3.3, min: 0.3, max: 15, group: 'Plate panel', role: 'geometry', design: true, help: 'Distance between the deck beams that support the stiffeners.' },
    { key: 'tp', label: 'Plate thickness', unit: 'mm', type: 'number', default: 10, min: 3, max: 60, group: 'Plate panel', role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'edges', label: 'Panel edge condition', type: 'select', default: 'clamped', group: 'Plate panel', role: 'bc', options: [{ value: 'ss', label: 'Simply supported (single loaded panel, conservative)' }, { value: 'clamped', label: 'Clamped (continuous plating, all panels loaded)' }] },
    { key: 'mat', label: 'Steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Plate panel', role: 'material' },
    { key: 'q', label: 'Design area load (factored)', unit: 'kPa', type: 'number', default: 20, min: 0, max: 500, group: 'Loads', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Factored permanent plus variable deck load, e.g. 1.3 × 2 + 1.5 × 10 kPa for a laydown area.' },
    { key: 'P', label: 'Design patch load (factored)', unit: 'kN', type: 'number', default: 35, min: 0, max: 5000, group: 'Loads', role: 'bc', help: 'Wheel, skid shoe or equipment foot on one panel, at its centre.' },
    { key: 'pu', label: 'Patch length (along the stiffener)', unit: 'mm', type: 'number', default: 300, min: 20, max: 3000, group: 'Loads', role: 'load' },
    { key: 'pv', label: 'Patch width (across the stiffeners)', unit: 'mm', type: 'number', default: 300, min: 20, max: 3000, group: 'Loads', role: 'load' },
    { key: 'hw', label: 'Stiffener web height', unit: 'mm', type: 'number', default: 240, min: 40, max: 1200, group: 'Stiffener', role: 'geometry', design: true },
    { key: 'tw', label: 'Stiffener web thickness', unit: 'mm', type: 'number', default: 10, min: 4, max: 60, group: 'Stiffener', role: 'geometry' },
    { key: 'bf', label: 'Stiffener flange width (0 = flat bar)', unit: 'mm', type: 'number', default: 100, min: 0, max: 600, group: 'Stiffener', role: 'geometry' },
    { key: 'tf', label: 'Stiffener flange thickness', unit: 'mm', type: 'number', default: 12, min: 0, max: 60, group: 'Stiffener', role: 'geometry' },
    { key: 'cont', label: 'Stiffener and beam continuity', type: 'select', default: 'continuous', group: 'Stiffener', role: 'bc', options: [{ value: 'simple', label: 'Simply supported spans (M = wL²/8)' }, { value: 'continuous', label: 'Continuous over supports (fixed-end moments, M = wL²/12)' }] },
    { key: 'Lb', label: 'Beam span', unit: 'm', type: 'number', default: 10, min: 1, max: 40, group: 'Supporting beam', role: 'geometry' },
    { key: 'bh', label: 'Beam depth', unit: 'mm', type: 'number', default: 700, min: 150, max: 4000, group: 'Supporting beam', role: 'geometry', design: true },
    { key: 'bb', label: 'Beam flange width', unit: 'mm', type: 'number', default: 300, min: 80, max: 1500, group: 'Supporting beam', role: 'geometry' },
    { key: 'btw', label: 'Beam web thickness', unit: 'mm', type: 'number', default: 12, min: 4, max: 80, group: 'Supporting beam', role: 'geometry' },
    { key: 'btf', label: 'Beam flange thickness', unit: 'mm', type: 'number', default: 22, min: 6, max: 120, group: 'Supporting beam', role: 'geometry', design: true },
    { key: 'cD', label: 'Column diameter', unit: 'mm', type: 'number', default: 610, min: 100, max: 3000, group: 'Column', role: 'geometry', design: true },
    { key: 'ct', label: 'Column wall thickness', unit: 'mm', type: 'number', default: 25, min: 5, max: 120, group: 'Column', role: 'geometry', design: true },
    { key: 'cL', label: 'Column height', unit: 'm', type: 'number', default: 7, min: 1, max: 30, group: 'Column', role: 'geometry' },
    { key: 'cK', label: 'Column effective length factor', unit: '-', type: 'number', default: 1.0, min: 0.5, max: 2.5, group: 'Column', role: 'bc', help: '0.7–1.0 braced; above 1.0 for sway frames (see the buckling module).' },
    { key: 'cN', label: 'Column design axial compression', unit: 'kN', type: 'number', default: 4500, min: 0, max: 2e5, group: 'Column', role: 'bc' },
    { key: 'cM', label: 'Column design end moment', unit: 'kN·m', type: 'number', default: 600, min: 0, max: 2e5, group: 'Column', role: 'bc' },
    { key: 'gM', label: 'Material factor γM', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.5, group: 'Acceptance criteria', role: 'criteria', help: `1.15 in NORSOK N-004 / DNV-OS-C101; 1.0 in EN 1993-1-1. ${CONFIRM}` },
    { key: 'defRatio', label: 'Deflection limit (span / value)', unit: '-', type: 'number', default: 250, min: 100, max: 1000, group: 'Acceptance criteria', role: 'criteria', help: 'Checked under the factored load divided by 1.4 as an estimate of the characteristic load.' },
    { key: 'nEl', label: 'Plate elements across the panel width', type: 'int', default: 12, min: 4, max: 48, group: 'Solver', role: 'control', help: 'MITC4 plate elements across the stiffener spacing; near-square elements along the span.' },
  ],
  mesh: { key: 'nEl', label: 'Elements across the panel', levels: [6, 12, 24], qoi: 'wFE' },
  compute(inp, ctx) {
    // SI: mm → m, kPa → Pa, kN → N.
    const st = steelOf(ctx, inp.mat), fyd = st.fy / inp.gM;
    const b = inp.s / 1e3, a = inp.Ls, t = inp.tp / 1e3, q = inp.q * 1e3, P = inp.P * 1e3;
    need(a >= b, 'The stiffener span must not be shorter than the stiffener spacing (the panel is defined with its long side along the stiffeners).');
    need(q > 0 || P > 0, 'Enter an area load or a patch load.');
    need(inp.pu / 1e3 <= a && inp.pv / 1e3 <= b, 'The patch load footprint must fit inside one plate panel.');
    const clamped = inp.edges === 'clamped', E = st.E, nu = st.nu, D = plateRigidity({ E, nu, t });
    // --- plate, uniform load (series) -----------------------------------------------------------------------
    let wPl = 0, sPl = 0;
    if (q > 0) {
      if (clamped) { const r = run(() => ritzPlate({ a, b, t, E, nu, edges: { x0: 'C', x1: 'C', y0: 'C', y1: 'C' }, load: { type: 'uniform', q }, nBasis: 10 }), 'Plate (Ritz)'); wPl = r.wMax; sPl = Math.max(Math.abs(r.sigmaX), Math.abs(r.sigmaY), ...Object.values(r.edgeMoments).map((m) => 6 * Math.abs(m) / (t * t))); }
      else { const r = run(() => navierPlate({ a, b, t, E, nu, load: { type: 'uniform', q } }), 'Plate (Navier)'); wPl = Math.abs(r.w); sPl = Math.max(Math.abs(r.sigmaX), Math.abs(r.sigmaY)); }
    }
    // --- plate, uniform load (Mindlin FE) -------------------------------------------------------------------
    const nyE = inp.nEl, nxE = Math.max(inp.nEl, Math.min(5 * inp.nEl, Math.round(inp.nEl * a / b)));
    const mesh = rectangleMesh({ lx: a, ly: b, nx: nxE, ny: nyE });
    const edgeNodes = [...new Set([...mesh.nodeSets.left, ...mesh.nodeSets.right, ...mesh.nodeSets.bottom, ...mesh.nodeSets.top])];
    const fe = run(() => solvePlate({ mesh, E, nu, t, pressure: Math.max(q, 1), fix: [clamped ? { nodes: edgeNodes, w: 0, rx: 0, ry: 0 } : { nodes: edgeNodes, w: 0 }] }), 'Plate finite elements');
    const feScale = q > 0 ? 1 : 0, wFE = fe.maxW * feScale;
    let sFE = 0; for (let i = 0; i < mesh.nn; i++) sFE = Math.max(sFE, 6 * Math.abs(fe.moments.mx[i]) / (t * t), 6 * Math.abs(fe.moments.my[i]) / (t * t));
    sFE *= feScale;
    const ld = q > 0 ? largeDeflectionPlate({ a, b, t, E, nu, q, edges: 'immovable' }) : { w0: 0, ratio: 0, w0Linear: 0 };
    // --- plate, patch load -----------------------------------------------------------------------------------
    const patch = P > 0 ? run(() => deckPatchLoad({ a, b, t, E, nu, P, u: inp.pu / 1e3, v: inp.pv / 1e3, edges: clamped ? 'clamped' : 'ss', fy: st.fy }), 'Patch load') : { sigmaMax: 0, w: 0, vonMises: 0 };
    // --- stiffener with effective plate flange (shear lag: be ≤ L/4) ----------------------------------------
    const cM = inp.cont === 'continuous' ? 12 : 8, cD = inp.cont === 'continuous' ? 1 : 5, cP = inp.cont === 'continuous' ? 8 : 4, cPD = inp.cont === 'continuous' ? 192 : 48;
    const be = Math.min(b, a / 4);
    const sp = run(() => stiffenerWithPlate({ be, tp: t, hw: inp.hw / 1e3, tw: inp.tw / 1e3, bf: inp.bf / 1e3, tf: inp.tf / 1e3 }), 'Stiffener section');
    const wS = q * b, Ms = wS * a * a / cM + P * a / cP, Vs = wS * a / 2 + P, sigS = Ms / Math.min(sp.Zflange, sp.Zplate), tauS = Vs / ((inp.hw / 1e3) * (inp.tw / 1e3));
    const dS = (cD * wS * a ** 4 / (384 * E * sp.I) + P * a ** 3 / (cPD * E * sp.I)) / 1.4;
    // --- supporting beam: line load from the stiffeners over its tributary width Ls ---------------------------
    const bsec = run(() => iSection({ h: inp.bh / 1e3, b: inp.bb / 1e3, tw: inp.btw / 1e3, tf: inp.btf / 1e3, welded: true }), 'Beam section');
    const wB = q * a + 1.3 * bsec.A * st.rho * G, L = inp.Lb, Mb = wB * L * L / cM, Vb = wB * L / 2;
    const bc = run(() => ec3MemberCheck({ section: bsec, fy: st.fy, E, G: st.G, MyEd: Mb, VzEd: Vb, Lcry: L, Lcrz: a, LLT: a, gammaM0: inp.gM, gammaM1: inp.gM }), 'Beam check');
    const dB = (cD * wB * L ** 4 / (384 * E * bsec.Iy) + (wB * L * L / 8) / (st.G * bsec.Asz)) / 1.4;
    // --- column ------------------------------------------------------------------------------------------------
    need(inp.ct < inp.cD / 2, 'Column wall thickness must be smaller than its radius.');
    const csec = tubeSection({ D: inp.cD / 1e3, t: inp.ct / 1e3 });
    const cc = run(() => tubularMemberCheck({ code: 'ISO', section: csec, fy: st.fy, L: inp.cL, forces: { N: -inp.cN * 1e3, My: inp.cM * 1e3, Mz: 0, Vy: 0, T: 0 }, options: { Ky: inp.cK, Kz: inp.cK, E } }), 'Column check');
    const pEuler = Math.PI ** 2 * E * csec.Iy / (inp.cK * inp.cL) ** 2;
    const sPlate = Math.max(sPl, sFE);
    const checks = [
      check('Plate bending — area load', sPlate / 1e6, fyd / 1e6, 'MPa', 'DNV-RP-C201 §5 (lateral load); elastic plate theory, fy/γM', 'Extreme-fibre bending stress of the plate panel under the area load (larger of series and finite-element values). Exceeding it starts surface yielding; a plate still has yield-line reserve, but permanent dishing follows.'),
      check('Plate bending — patch load', patch.sigmaMax / 1e6, fyd / 1e6, 'MPa', 'ISO 19901-3 §7.6; elastic plate theory under a rectangular patch', 'Local bending under a wheel or equipment foot. Exceeding it dents the plate between stiffeners.'),
      check('Plate deflection', (Math.max(wPl, wFE) / 1.4) * 1e3, (b / inp.defRatio) * 1e3, 'mm', `Serviceability: spacing / ${inp.defRatio}`, 'Limits visible dishing and ponding between stiffeners under characteristic load.'),
      check('Stiffener bending', sigS / 1e6, fyd / 1e6, 'MPa', 'DNV-RP-C201 §7; EN 1993-1-5 §3 (effective width for shear lag)', 'Bending stress of the stiffener with its effective plate flange at the more stressed fibre. Exceeding it yields the stiffener flange and the panel loses its support.'),
      check('Stiffener web shear', tauS / 1e6, fyd / Math.sqrt(3) / 1e6, 'MPa', 'EN 1993-1-1 §6.2.6', 'Shear stress in the stiffener web at the support; exceeding it yields the web in shear near the beam connection.'),
      check('Stiffener deflection', dS * 1e3, (a / inp.defRatio) * 1e3, 'mm', `Serviceability: span / ${inp.defRatio}`, 'Limits sag of the stringers under characteristic load, which affects equipment alignment and drainage.'),
      check(`Beam strength — ${bc.governing.name}`, bc.governing.uc, 1, 'interaction', `EN 1993-1-1 ${bc.governing.ref || '§6.2–6.3'}`, 'Bending, shear and lateral-torsional buckling of the deck beam between stringer restraints. Above 1.0 the beam yields or buckles sideways.'),
      check('Beam deflection', dB * 1e3, (L / inp.defRatio) * 1e3, 'mm', `Serviceability: span / ${inp.defRatio} (bending + shear deflection)`, 'Limits sag of the deck beam under characteristic load.'),
      check(`Column — ${cc.governing.name}`, cc.governing.uc, 1, 'interaction', 'ISO 19902:2020 §13.2–13.3', 'Axial compression with bending of the tubular column, including column buckling. Above 1.0 the column buckles or yields.'),
    ];
    const gov = checks.reduce((p, c) => (c.uc > p.uc ? c : p)), ok = gov.uc <= 1;
    // Required plate thickness from σ ∝ 1/t² for the bending checks.
    const tReq = inp.tp * Math.sqrt(Math.max(sPlate, patch.sigmaMax) / fyd);
    const kgm2 = t * st.rho + (sp.As * st.rho) / b, kgm2Req = (tReq / 1e3) * st.rho + (sp.As * st.rho) / b;
    const sc = wFE > 0 ? 0.1 * b / (fe.maxW) : 0;
    return {
      kpis: [
        { key: 'wPlate', label: 'Plate deflection (series solution)', unit: 'mm', value: wPl * 1e3, explain: 'Largest deflection of the panel under the factored area load from Kirchhoff plate theory (Navier series or Ritz).' },
        { key: 'wFE', label: 'Plate deflection (Mindlin finite elements)', unit: 'mm', value: wFE * 1e3, explain: 'Same quantity from MITC4 plate elements; it should agree with the series value within a few per cent once the mesh has converged.' },
        { key: 'sigPlate', label: 'Plate bending stress — area load', unit: 'MPa', value: sPlate / 1e6, limit: fyd / 1e6, sense: 'max', explain: 'Extreme-fibre stress at the panel centre or its clamped edges. Below fy/γM the plate stays elastic.' },
        { key: 'sigPatch', label: 'Plate bending stress — patch load', unit: 'MPa', value: patch.sigmaMax / 1e6, limit: fyd / 1e6, sense: 'max', explain: 'Largest surface stress under the wheel or equipment foot and at the panel edges.' },
        { key: 'wOverT', label: 'Large-deflection ratio w/t', unit: '-', value: ld.ratio, explain: 'Centre deflection over thickness with membrane action. Above about 0.5 small-deflection theory overestimates the bending stress and the membrane reserve becomes significant.' },
        { key: 'sigStiff', label: 'Stiffener bending stress', unit: 'MPa', value: sigS / 1e6, limit: fyd / 1e6, sense: 'max', explain: 'Stress at the critical fibre of the stiffener acting with its effective plate flange.' },
        { key: 'ucBeam', label: 'Beam utilisation', unit: '-', value: round(bc.governing.uc, 5), limit: 1, sense: 'max', explain: 'Governing EN 1993-1-1 utilisation of the supporting deck beam (bending, shear, lateral-torsional buckling).' },
        { key: 'ucColumn', label: 'Column utilisation', unit: '-', value: round(cc.governing.uc, 5), limit: 1, sense: 'max', explain: 'Governing ISO 19902 utilisation of the tubular column under axial compression and bending.' },
        { key: 'pEuler', label: 'Column Euler load', unit: 'kN', value: pEuler / 1e3, explain: 'Elastic critical load π²EI/(KL)²: the column design load should be a small fraction of it.' },
        { key: 'kgm2', label: 'Plate and stiffener steel', unit: 'kg/m²', value: kgm2, explain: 'Mass of plating and stiffeners per unit deck area: the figure to minimise when choosing spacing against thickness.' },
      ],
      checks,
      plots: [
        meshPlot({ title: 'Plate deflection (finite elements)', mesh, values: Array.from(fe.w, (v) => Math.abs(v) * feScale * 1e3), valueLabel: 'Deflection', valueUnit: 'mm', z: Array.from(fe.w, (v) => -Math.abs(v) * sc),
          explain: `One panel between stiffeners (x along the stiffeners) under the area load, deflection exaggerated ${fmt(sc, 2)} times. The dish should be symmetric with its maximum at the centre.` }),
        { type: 'bar', title: 'Component utilisation', categories: checks.map((c) => c.name), series: [{ name: 'Utilisation', y: checks.map((c) => round(c.uc)) }], y: { label: 'Utilisation', unit: '-' }, explain: 'All component checks side by side: the tallest bar governs; an even profile near 0.8–0.9 is an efficient design.' },
        { type: 'line', title: 'Plate stress vs thickness', x: { label: 'Plate thickness', unit: 'mm' }, y: { label: 'Bending stress', unit: 'MPa' },
          series: [{ name: 'Area load', x: [6, 8, 10, 12, 15, 20, 25].map((v) => v), y: [6, 8, 10, 12, 15, 20, 25].map((v) => round(sPlate / 1e6 * (inp.tp / v) ** 2)) }, { name: 'Patch load', x: [6, 8, 10, 12, 15, 20, 25], y: [6, 8, 10, 12, 15, 20, 25].map((v) => round(patch.sigmaMax / 1e6 * (inp.tp / v) ** 2)) }, { name: 'fy/γM', x: [6, 25], y: [fyd / 1e6, fyd / 1e6], dash: true }],
          explain: 'Elastic bending stress scales with 1/t²: read off the thickness where each curve drops below the dashed design strength.' },
      ],
      tables: [{ title: 'Component results', columns: [{ key: 'c', label: 'Component' }, { key: 'M', label: 'Moment', unit: 'kN·m' }, { key: 'V', label: 'Shear / axial', unit: 'kN' }, { key: 's', label: 'Stress', unit: 'MPa' }, { key: 'd', label: 'Deflection', unit: 'mm' }],
        rows: [{ c: 'Plate (per metre width)', M: round(sPlate * t * t / 6 / 1e3), V: round(q * b / 2 / 1e3), s: round(sPlate / 1e6), d: round(Math.max(wPl, wFE) * 1e3) }, { c: 'Stiffener + effective plate', M: round(Ms / 1e3), V: round(Vs / 1e3), s: round(sigS / 1e6), d: round(dS * 1.4 * 1e3) },
          { c: 'Supporting beam', M: round(Mb / 1e3), V: round(Vb / 1e3), s: round(Mb / bsec.Wely / 1e6), d: round(dB * 1.4 * 1e3) }, { c: 'Column', M: inp.cM, V: inp.cN, s: round((inp.cN * 1e3 / csec.A + inp.cM * 1e3 / csec.Wely) / 1e6), d: 0 }],
        explain: 'Design moment, shear (axial force for the column), elastic stress and deflection under the factored loads for each component in the load path from plate to column.' }],
      notes: [`${gov.name} governs at ${fmt(gov.uc)}.`, `Plate rigidity D = ${fmt(D / 1e3)} kN·m; effective plate flange of the stiffener ${fmt(be * 1e3)} mm of ${inp.s} mm.`, `Series and finite-element plate deflections: ${fmt(wPl * 1e3)} and ${fmt(wFE * 1e3)} mm (${nxE} × ${nyE} elements).`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `All components pass; ${gov.name.toLowerCase()} governs at ${fmt(gov.uc)}. Plate needs ${fmt(tReq, 2)} mm for elastic bending against ${inp.tp} mm provided.`,
          tReq < 0.8 * inp.tp ? `The plate has reserve: ${Math.max(6, Math.ceil(tReq))} mm would satisfy bending — check the minimum thickness for corrosion, welding distortion and dropped objects before reducing.` : 'Keep the scantlings; confirm the in-plane (buckling) check of the plate when it acts as the flange of the deck girders.',
          'DNV-RP-C201 (2010) §5–§7; ISO 19901-3:2014 §7.6', `Plate and stiffeners weigh ${fmt(kgm2)} kg/m²; at the required thickness ${fmt(kgm2Req)} kg/m², i.e. ${steelCo2Text(Math.abs(kgm2 - kgm2Req) * 600 / 1e3)} per 600 m² deck.`)
          : adv('engineering', 'fail', `${gov.name} fails at utilisation ${fmt(gov.uc)} (${fmt(gov.demand)} against ${fmt(gov.capacity)} ${gov.unit}).`,
            gov.name.startsWith('Plate') ? `Increase the plate to at least ${Math.ceil(tReq)} mm or reduce the stiffener spacing to about ${fmt(inp.s / Math.sqrt(gov.uc), 2)} mm.` : gov.name.startsWith('Stiffener') ? `Increase the stiffener section modulus by ${fmt((gov.uc - 1) * 100, 2)} % (deeper web is most efficient) or shorten the stiffener span.` : gov.name.startsWith('Beam') ? `Increase the beam section modulus by about ${fmt((gov.uc / 0.9 - 1) * 100, 2)} % or add an intermediate support.` : `Increase the column wall thickness or diameter by about ${fmt((gov.uc / 0.9 - 1) * 100, 2)} % in area, or brace the column to shorten its effective length.`,
            'DNV-RP-C201 (2010); EN 1993-1-1:2005 §6; ISO 19902:2020 §13', `Closer stiffeners usually meet a plate check with less steel than a thicker plate: 1 mm of plate is ${fmt(st.rho / 1e3)} kg/m² (${fmt(st.rho / 1e3 * 1.9)} kg CO₂e/m² at the 1.9 t/t screening factor).`),
        adv('economic', 'info', `Plate and stiffeners use ${fmt(kgm2)} kg of steel per m² of deck, with ${fmt(1e3 / inp.s, 2)} stiffener weld runs per metre width.`,
          'Balance plate thickness against stiffener count: fewer, deeper stiffeners with a thicker plate cut welding and coating hours, which dominate fabricated deck cost per tonne; standardise on one plate thickness per deck.', 'Design-for-fabrication practice; AACE RP 18R-97 estimate classes'),
      ],
    };
  },
  verify: [
    { // D = 210e9·0.01³/(12·0.91) = 19 230.8 N·m; w = 0.00406 q a⁴/D = 0.00406·1e4/19230.8 = 2.111 mm; M = 0.0479 q a² → σ = 6M/t² = 28.74 MPa.
      name: 'Simply supported square plate, uniform load: w = 0.00406 q a⁴/D, M = 0.0479 q a²', source: 'Timoshenko & Woinowsky-Krieger, Theory of Plates and Shells, Table 8',
      inputs: { s: 1000, Ls: 1, tp: 10, edges: 'ss', q: 10, P: 0, nEl: 24 }, expect: [{ kpi: 'wPlate', value: 2.111, tol: 0.005 }, { kpi: 'wFE', value: 2.111, tol: 0.015 }, { kpi: 'sigPlate', value: 28.74, tol: 0.02 }] },
    { // Clamped square plate: w = 0.00126 q a⁴/D = 0.6552 mm; edge moment 0.0513 q a² → σ = 30.78 MPa.
      name: 'Clamped square plate, uniform load: w = 0.00126 q a⁴/D, edge M = 0.0513 q a²', source: 'Timoshenko & Woinowsky-Krieger, Theory of Plates and Shells, Table 35',
      inputs: { s: 1000, Ls: 1, tp: 10, edges: 'clamped', q: 10, P: 0, nEl: 24 }, expect: [{ kpi: 'wPlate', value: 0.6552, tol: 0.01 }, { kpi: 'wFE', value: 0.6552, tol: 0.03 }, { kpi: 'sigPlate', value: 30.78, tol: 0.04 }] },
  ],
};
