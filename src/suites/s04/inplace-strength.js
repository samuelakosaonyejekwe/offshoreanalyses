// In-place structural strength of a subsea template / protection frame by three-dimensional frame analysis.
import { analyzeFrame, memberDiagram, sectionStress, refineModel } from '../../engine/structural/frame3d.js';
import { bucklingAnalysis } from '../../engine/structural/framebuckling.js';
import { resolveSection } from '../../engine/structural/sections.js';
import { tubularMemberCheck } from '../../engine/structural/tubular.js';
import { fromImportedFrame } from '../../engine/structural/builders.js';
import { need, run, fmt, round, adv, steelCo2Text, CONFIRM } from './_common.js';

/** Parametric template: base grid of beams, four corner posts and a top perimeter frame. */
function templateModel(inp) {
  const { Lx, Ly, hT, nx, ny } = inp;
  const nodes = [], members = [], id = (i, j) => i * (ny + 1) + j;
  for (let i = 0; i <= nx; i++) for (let j = 0; j <= ny; j++) nodes.push([-Lx / 2 + Lx * i / nx, -Ly / 2 + Ly * j / ny, 0]);
  for (let i = 0; i < nx; i++) for (let j = 0; j <= ny; j++) members.push({ a: id(i, j), b: id(i + 1, j), section: 'base', group: 'Base beam' });
  for (let i = 0; i <= nx; i++) for (let j = 0; j < ny; j++) members.push({ a: id(i, j), b: id(i, j + 1), section: 'base', group: 'Base beam' });
  const corners = [id(0, 0), id(nx, 0), id(nx, ny), id(0, ny)];
  const top = corners.map((c) => { nodes.push([nodes[c][0], nodes[c][1], hT]); return nodes.length - 1; });
  corners.forEach((c, k) => members.push({ a: c, b: top[k], section: 'post', group: 'Post' }));
  for (let k = 0; k < 4; k++) members.push({ a: top[k], b: top[(k + 1) % 4], section: 'top', group: 'Top frame' });
  const base = []; for (let i = 0; i <= nx; i++) for (let j = 0; j <= ny; j++) base.push(id(i, j));
  let inner = []; for (let i = 1; i < nx; i++) for (let j = 1; j < ny; j++) inner.push(id(i, j));
  if (!inner.length) inner = base;
  return { nodes, members, base, corners, top, inner, hub: id(nx, Math.round(ny / 2)), trib: (Lx * Ly) / base.length };
}

export default {
  id: 's04-inplace-strength', suite: 4, title: 'In-place structural strength',
  summary: 'Member forces, stresses, deflections, support reactions and elastic buckling factor of the template frame under weight, equipment, environmental and interface loads.',
  covers: [2], fidelity: 'numerical',
  standards: ['ISO 19902:2020 §13 (tubular member strength)', 'ISO 13628-1:2005 / API RP 17A (subsea structures)', 'DNV-OS-C101 (2019) Ch.2 Sec.4', 'NORSOK U-001 (2015) §5 (subsea structural design)'],
  models: ['Force and moment equilibrium', 'Euler–Bernoulli and Timoshenko beam formulations (3-D direct stiffness)', 'Hooke’s law', 'von Mises yielding', 'Linear (eigenvalue) frame buckling', 'Winkler soil springs under the base', 'Archimedes buoyancy of flooded members'],
  limits: 'Linear elastic frame of tubular members with rigid joints; soil is a bed of linear springs or rigid supports, without lift-off. Joint punching capacity, plate and shell action of mudmat plating, weld details and nonlinear finite-element collapse are outside the module. Imported frames keep their own sections; non-tubular sections are checked by von Mises stress only.',
  publishes: ['subsea.frameUc', 'subsea.frameSteel', 'quantities.steel_t'],
  searchTerms: 'subsea template protection structure frame analysis in-place strength',
  inputs: [
    { key: 'geom', label: 'Imported frame model (optional)', type: 'geometry', kinds: ['frame'], default: null, group: 'Geometry', role: 'geometry', help: 'SDNF, structured JSON or FE beam model. When empty the parametric template below is analysed. Imported frames are supported at their lowest nodes; equipment weight is spread over the remaining lower nodes and the horizontal load over the top nodes.' },
    { key: 'Lx', label: 'Template length (x, load direction)', unit: 'm', type: 'number', default: 10, min: 1, max: 80, group: 'Geometry', role: 'geometry', design: true, bind: 'subsea.matB' },
    { key: 'Ly', label: 'Template width (y)', unit: 'm', type: 'number', default: 14, min: 1, max: 80, group: 'Geometry', role: 'geometry', design: true, bind: 'subsea.matL' },
    { key: 'hT', label: 'Height of the top frame', unit: 'm', type: 'number', default: 3.5, min: 0.5, max: 30, group: 'Geometry', role: 'geometry' },
    { key: 'nx', label: 'Base bays along x', type: 'int', default: 4, min: 1, max: 12, group: 'Geometry', role: 'geometry' },
    { key: 'ny', label: 'Base bays along y', type: 'int', default: 4, min: 1, max: 12, group: 'Geometry', role: 'geometry' },
    { key: 'Db', label: 'Base beam diameter', unit: 'mm', type: 'number', default: 508, min: 60, max: 3000, group: 'Members', role: 'geometry', design: true },
    { key: 'tb', label: 'Base beam wall thickness', unit: 'mm', type: 'number', default: 16, min: 4, max: 100, group: 'Members', role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'Dp', label: 'Post diameter', unit: 'mm', type: 'number', default: 406.4, min: 60, max: 3000, group: 'Members', role: 'geometry', design: true },
    { key: 'tpst', label: 'Post wall thickness', unit: 'mm', type: 'number', default: 16, min: 4, max: 100, group: 'Members', role: 'geometry', design: true },
    { key: 'Dt', label: 'Top frame diameter', unit: 'mm', type: 'number', default: 323.9, min: 60, max: 3000, group: 'Members', role: 'geometry', design: true },
    { key: 'tt', label: 'Top frame wall thickness', unit: 'mm', type: 'number', default: 12.7, min: 4, max: 100, group: 'Members', role: 'geometry', design: true },
    { key: 'mat', label: 'Steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Members', role: 'material' },
    { key: 'Kbuck', label: 'Effective length factor for member checks', unit: '-', type: 'number', default: 1.0, min: 0.5, max: 2.5, group: 'Members', role: 'criteria' },
    { key: 'support', label: 'Support condition', type: 'select', default: 'springs', group: 'Boundary conditions', role: 'bc',
      options: [{ value: 'springs', label: 'Mudmat on soil springs (all base nodes)' }, { value: 'pinned', label: 'Pinned at the four corners' }, { value: 'fixed', label: 'Fixed at the four corners (piles)' }] },
    { key: 'ks', label: 'Vertical subgrade modulus', unit: 'kPa/m', type: 'number', default: 400, min: 1, max: 1e6, group: 'Boundary conditions', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.4 }, showIf: { key: 'support', in: ['springs'] }, help: 'Soil pressure per unit settlement; horizontal springs are taken at half this value.' },
    { key: 'Weq', label: 'Equipment submerged weight on the base', unit: 'kN', type: 'number', default: 800, min: 0, max: 2e5, group: 'Loads', role: 'bc', uncertain: { dist: 'normal', cov: 0.05 }, help: 'Manifold, trees, valves and modules, spread over the interior base nodes.' },
    { key: 'Wroof', label: 'Vertical load on the top frame', unit: 'kN', type: 'number', default: 150, min: 0, max: 2e5, group: 'Loads', role: 'bc', help: 'Roof hatches, protection covers and trawl-gear or intervention loads, shared by the top corners.' },
    { key: 'Hx', label: 'Horizontal load on the top frame (x)', unit: 'kN', type: 'number', default: 120, min: 0, max: 1e5, group: 'Loads', role: 'bc', bind: 'subsea.hydroFH', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Hydrodynamic or trawl load, shared by the top corners.' },
    { key: 'Fhub', label: 'Pipeline / jumper axial reaction at the hub (x)', unit: 'kN', type: 'number', default: 150, min: -1e5, max: 1e5, group: 'Loads', role: 'bc', bind: 'pipeline.endForce', help: 'Applied at the mid-point of the +x base edge.' },
    { key: 'Vhub', label: 'Vertical reaction at the hub (down +)', unit: 'kN', type: 'number', default: 50, min: -1e5, max: 1e5, group: 'Loads', role: 'bc' },
    { key: 'swFac', label: 'Self-weight factor', unit: '-', type: 'number', default: 1.0, min: 0, max: 2, group: 'Loads', role: 'bc', help: 'Submerged self-weight of the flooded members; 0 switches it off.' },
    { key: 'gammaF', label: 'Load factor on all loads', unit: '-', type: 'number', default: 1.3, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: `ULS action factor. ${CONFIRM}` },
    { key: 'defLimit', label: 'Deflection limit (span / value)', unit: '-', type: 'number', default: 200, min: 50, max: 1000, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'nSub', label: 'Elements per member', type: 'int', default: 2, min: 1, max: 16, group: 'Solver', role: 'control', help: 'Beam elements are exact at their nodes for nodal loads: subdividing adds sampling points along the members (mid-span deflection) and refines the buckling mode.' },
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'maxDisp' },
  compute(inp, ctx) {
    // SI: mm → m, kN → N, kPa/m → Pa/m.
    const mat = ctx.material(inp.mat);
    need(mat, `Steel '${inp.mat}' is not in the materials library.`);
    for (const [d, t, n] of [[inp.Db, inp.tb, 'base beam'], [inp.Dp, inp.tpst, 'post'], [inp.Dt, inp.tt, 'top frame']]) need(t < d / 2, `Wall thickness of the ${n} must be smaller than its radius.`);
    const steel = { E: mat.E, nu: mat.nu ?? 0.3, rho: mat.rho ?? 7850, fy: mat.fy };
    const f = inp.gammaF, kN = 1e3 * f;
    let model, loadNodes, topNodes, hub, supNodes, trib, span;
    if (inp.geom) {
      need(inp.geom.kind === 'frame', 'The imported model must be a frame (nodes and members).');
      model = run(() => fromImportedFrame(inp.geom, { defaultSection: { type: 'tube', D: inp.Db / 1e3, t: inp.tb / 1e3 }, defaultMaterial: 'steel', supportType: inp.support === 'fixed' ? 'fixed' : 'pinned' }), 'Imported frame');
      model.members.forEach((m) => { m.material = 'steel'; m.flooded = true; });
      supNodes = model.supports.map((s) => s.node);
      const zs = model.nodes.map((n) => n[2]), zmax = Math.max(...zs), zmin = Math.min(...zs), xs = model.nodes.map((n) => n[0]), ys = model.nodes.map((n) => n[1]);
      need(zmax - zmin > 1e-6, 'The imported frame is flat: it needs members above the support level to be analysed here.');
      topNodes = model.nodes.map((n, i) => i).filter((i) => zs[i] >= zmax - 1e-3 * (zmax - zmin));
      loadNodes = model.nodes.map((n, i) => i).filter((i) => !supNodes.includes(i) && !topNodes.includes(i));
      if (!loadNodes.length) loadNodes = topNodes;
      hub = loadNodes[0]; span = Math.max(Math.max(...xs) - Math.min(...xs), Math.max(...ys) - Math.min(...ys));
      trib = (Math.max(...xs) - Math.min(...xs)) * (Math.max(...ys) - Math.min(...ys)) / supNodes.length;
      if (inp.support === 'springs') model.supports = supNodes.map((n) => ({ node: n, spring: [0.5 * inp.ks * 1e3 * trib, 0.5 * inp.ks * 1e3 * trib, inp.ks * 1e3 * trib, 0, 0, 0] }));
      need(supNodes.length >= 3 || inp.support === 'fixed', 'The imported frame has fewer than three support nodes at its lowest level: it is under-constrained. Use fixed supports or add base nodes.');
    } else {
      const g = templateModel(inp);
      model = { nodes: g.nodes, members: g.members.map((m) => ({ ...m, material: 'steel', flooded: true })) };
      loadNodes = g.inner; topNodes = g.top; hub = g.hub; trib = g.trib; span = Math.max(inp.Lx / inp.nx, inp.Ly / inp.ny);
      supNodes = inp.support === 'springs' ? g.base : g.corners;
      model.supports = inp.support === 'springs'
        ? g.base.map((n) => ({ node: n, spring: [0.5 * inp.ks * 1e3 * trib, 0.5 * inp.ks * 1e3 * trib, inp.ks * 1e3 * trib, 0, 0, 0] }))
        : g.corners.map((n) => ({ node: n, fix: inp.support }));
    }
    model.materials = { steel };
    model.sections = { base: { type: 'tube', D: inp.Db / 1e3, t: inp.tb / 1e3 }, post: { type: 'tube', D: inp.Dp / 1e3, t: inp.tpst / 1e3 }, top: { type: 'tube', D: inp.Dt / 1e3, t: inp.tt / 1e3 } };
    const nodal = [];
    for (const n of loadNodes) nodal.push({ node: n, F: [0, 0, -inp.Weq * kN / loadNodes.length, 0, 0, 0] });
    for (const n of topNodes) nodal.push({ node: n, F: [inp.Hx * kN / topNodes.length, 0, -inp.Wroof * kN / topNodes.length, 0, 0, 0] });
    nodal.push({ node: hub, F: [inp.Fhub * kN, 0, -inp.Vhub * kN, 0, 0, 0] });
    const lc = { name: 'design', nodal };
    if (inp.swFac > 0) lc.selfWeight = { factor: inp.swFac * f, buoyancy: { waterline: 1e5, rhoW: 1025 } };
    model.loadCases = [lc];
    const secOf = (m) => resolveSection(typeof m.section === 'string' ? model.sections[m.section] : m.section);
    const orig = model.members.map((m) => { const a = model.nodes[m.a], b = model.nodes[m.b]; return { sec: secOf(m), L: Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]), group: m.group || 'Member' }; });
    const steelT = orig.reduce((s, o) => s + o.sec.A * o.L * steel.rho, 0) / 1e3;
    const ref = run(() => refineModel(model, inp.nSub), 'Mesh refinement');
    const sol = run(() => analyzeFrame(ref.model, { cases: ['design'] }), 'Frame analysis');
    const r = sol.results.design;
    need(r && r.displacements, 'The frame analysis returned no result: check supports and connectivity.');
    const P0 = ref.model.nodes, nN = P0.length, u = r.displacements, disp = new Array(nN);
    let maxDisp = 0, maxVert = 0;
    for (let i = 0; i < nN; i++) { disp[i] = Math.hypot(u[6 * i], u[6 * i + 1], u[6 * i + 2]); maxDisp = Math.max(maxDisp, disp[i]); maxVert = Math.max(maxVert, Math.abs(u[6 * i + 2])); }
    need(Number.isFinite(maxDisp), 'The frame is unstable (mechanism): add supports or check that all members are connected.');
    // Flexural deflection: vertical displacement relative to the rigid-body plane through the support nodes.
    let flex = maxVert;
    if (supNodes.length >= 3) {
      const K = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], rhs = [0, 0, 0];
      for (const n of supNodes) { const v = [1, model.nodes[n][0], model.nodes[n][1]]; for (let i = 0; i < 3; i++) { rhs[i] += v[i] * u[6 * n + 2]; for (let j = 0; j < 3; j++) K[i][j] += v[i] * v[j]; } }
      const det3 = (q) => q[0][0] * (q[1][1] * q[2][2] - q[1][2] * q[2][1]) - q[0][1] * (q[1][0] * q[2][2] - q[1][2] * q[2][0]) + q[0][2] * (q[1][0] * q[2][1] - q[1][1] * q[2][0]);
      const det = det3(K);
      if (Math.abs(det) > 1e-9) {
        const pl = [0, 1, 2].map((c) => det3(K.map((row, i) => row.map((v, j) => (j === c ? rhs[i] : v)))) / det);
        flex = 0;
        for (let i = 0; i < nN; i++) flex = Math.max(flex, Math.abs(u[6 * i + 2] - pl[0] - pl[1] * P0[i][0] - pl[2] * P0[i][1]));
      }
    }
    // Member checks on the original members, using the extreme section forces over their sub-elements.
    const rows = orig.map((o, i) => {
      let N = 0, My = 0, Mz = 0, V = 0, T = 0;
      for (const k of ref.memberMap[i]) {
        const d = memberDiagram(r, k, 5);
        const n1 = d.N[0], n2 = d.N[d.N.length - 1];
        if (Math.abs(n1) > Math.abs(N)) N = n1; if (Math.abs(n2) > Math.abs(N)) N = n2;
        My = Math.max(My, d.maxAbs.My); Mz = Math.max(Mz, d.maxAbs.Mz); V = Math.max(V, Math.hypot(d.maxAbs.Vy, d.maxAbs.Vz)); T = Math.max(T, d.maxAbs.T);
      }
      const ss = sectionStress({ section: o.sec, N, Vy: V, T, My, Mz });
      let uc = ss.vonMises / (steel.fy / 1.15), ref2 = 'von Mises, γM = 1.15';
      if (o.sec.type === 'tube') { const c = tubularMemberCheck({ code: 'ISO', section: o.sec, fy: steel.fy, L: o.L, forces: { N, My, Mz, Vy: V, T }, options: { Ky: inp.Kbuck, Kz: inp.Kbuck, E: steel.E } }); uc = c.governing.uc; ref2 = c.governing.name; }
      uc = Number.isFinite(uc) ? Math.min(uc, 99) : 99;
      return { i, group: o.group, L: o.L, N, M: Math.hypot(My, Mz), V, vm: ss.vonMises, uc, ref: ref2 };
    });
    const gov = rows.reduce((a, c) => (c.uc > a.uc ? c : a));
    const maxVM = Math.max(...rows.map((x) => x.vm));
    const reacts = r.reactions.map((x) => ({ node: x.node, Rx: x.R[0], Ry: x.R[1], Rz: x.R[2] }));
    const reactV = reacts.reduce((s, x) => s + x.Rz, 0), baseShear = Math.hypot(reacts.reduce((s, x) => s + x.Rx, 0), reacts.reduce((s, x) => s + x.Ry, 0));
    const maxReact = Math.max(...reacts.map((x) => x.Rz));
    const pSoil = inp.support === 'springs' ? maxReact / trib : 0;
    let lambda = 999;
    // Eigenvalue buckling on the refined mesh for small frames; one element per member otherwise (within about 1 % for sway modes, much faster).
    try { const bk = bucklingAnalysis(ref.model.nodes.length <= 60 ? ref.model : model, { loadCase: 'design', nModes: 1 }); if (bk.factors[0] > 0 && Number.isFinite(bk.factors[0])) lambda = Math.min(bk.factors[0], 999); } catch { lambda = 999; }
    const defAllow = span / inp.defLimit;
    const byGroup = {};
    for (const x of rows) { const g = byGroup[x.group] || (byGroup[x.group] = { group: x.group, uc: 0, N: 0, M: 0, vm: 0, n: 0 }); g.n++; if (x.uc > g.uc) { g.uc = x.uc; g.N = x.N; g.M = x.M; g.vm = x.vm; g.ref = x.ref; } }
    const groups = Object.values(byGroup);
    const checks = [
      ...groups.map((g) => ({ name: `Member strength — ${g.group}`, demand: round(g.uc), capacity: 1, unit: 'interaction', uc: Math.min(g.uc, 99), ref: `ISO 19902:2020 §13 (${g.ref})`, explain: 'Governing axial, bending, shear and beam-column interaction of this member group. Above 1.0 the member yields or buckles and the frame must redistribute load.' })),
      { name: 'Flexural deflection', demand: round(flex * 1e3), capacity: round(defAllow * 1e3), unit: 'mm', uc: Math.min(flex / defAllow, 99), ref: `Serviceability: bay span / ${inp.defLimit}`, explain: 'Limits bending deflection (rigid-body settlement and tilt on the supports removed) so that equipment interfaces and connector alignment are maintained.' },
      { name: 'Elastic frame buckling', demand: 3, capacity: round(lambda), unit: 'load factor (required / achieved)', uc: Math.min(3 / lambda, 99), ref: 'Linear eigenvalue buckling; a factor above about 3 allows first-order member checks (cf. EN 1993-1-1 §5.2.1 uses 10 for sway frames)', explain: 'Factor on the design loads at which the elastic frame buckles. A low factor means second-order effects amplify the member forces.' },
    ];
    const maxUc = gov.uc, ok = checks.every((c) => c.uc <= 1);
    // 3-D view with the deformed shape scaled for visibility.
    const size = Math.max(...ref.model.nodes.map((n) => Math.max(Math.abs(n[0]), Math.abs(n[1]), Math.abs(n[2]))));
    const sc = maxDisp > 0 ? 0.05 * size / maxDisp : 0;
    const P = ref.model.nodes;
    const govGroupSteel = orig.filter((o) => o.group === gov.group).reduce((s, o) => s + o.sec.A * o.L * steel.rho, 0) / 1e3;
    const advice = [
      ok ? adv('engineering', 'ok', `Highest member utilisation ${fmt(maxUc)} in the ${gov.group.toLowerCase()} group (${gov.ref}); maximum displacement ${fmt(maxDisp * 1e3)} mm; elastic buckling factor ${fmt(lambda)}.`,
        maxUc < 0.6 ? `The ${gov.group.toLowerCase()} members have reserve: wall thickness could come down by roughly ${fmt((1 - maxUc / 0.85) * 100)} % before reaching 0.85. Check joints, lifting and transport load cases, which often govern over in-place.` : 'Utilisation is in an efficient range. Verify tubular joints (punching shear) and the lifting and transport cases before fixing sizes.',
        'ISO 19902:2020 §13; ISO 13628-1 / API RP 17A; NORSOK U-001 (2015)',
        maxUc < 0.6 ? `Trimming the ${gov.group.toLowerCase()} group to 0.85 utilisation would save about ${steelCo2Text(govGroupSteel * (1 - maxUc / 0.85))}` : `The frame holds ${steelCo2Text(steelT)}`)
        : adv('engineering', 'fail', `${checks.reduce((a, c) => (c.uc > a.uc ? c : a)).name} is exceeded (member utilisation ${fmt(maxUc)}, flexural deflection ${fmt(flex * 1e3)} mm against ${fmt(defAllow * 1e3)} mm).`,
          maxUc > 1 ? `Increase the wall thickness of the ${gov.group.toLowerCase()} members by about ${fmt((maxUc / 0.9 - 1) * 100)} % (to reach 0.9), add a diagonal brace, or add a base bay to shorten the spans.` : 'Stiffen the base grid (larger diameter rather than thicker wall — stiffness grows with D³) or add bays.',
          'ISO 19902:2020 §13; ISO 13628-1 / API RP 17A; NORSOK U-001 (2015)',
          `Strengthening the ${gov.group.toLowerCase()} group adds about ${steelCo2Text(govGroupSteel * Math.max(maxUc / 0.9 - 1, 0.05))} A larger diameter at the same mass is the lower-carbon route to stiffness.`),
      adv('economic', ok ? 'info' : 'warn', `Primary steel in the analysed frame is ${fmt(steelT)} t; the ${gov.group.toLowerCase()} group (${fmt(govGroupSteel)} t) sets the size.`,
        'Fabricated subsea steel is paid per tonne and lifted per tonne: keep the frame under the next crane-capacity or deck-area threshold of the installation vessel, and standardise on two or three tubular sizes to cut fabrication hours.',
        'AACE RP 18R-97 estimate classes; design-to-installation practice (DNV-ST-N001 lift weight contingencies)'),
    ];
    return {
      kpis: [
        { key: 'maxUc', label: 'Highest member utilisation', unit: '-', value: maxUc, limit: 1, sense: 'max', explain: 'Largest strength or beam-column utilisation of any member under the factored loads. Above 1.0 the member is overstressed.' },
        { key: 'maxDisp', label: 'Maximum nodal displacement', unit: 'mm', value: maxDisp * 1e3, explain: 'Largest translation of any node, including rigid-body settlement on the soil springs.' },
        { key: 'flex', label: 'Flexural deflection', unit: 'mm', value: flex * 1e3, limit: defAllow * 1e3, sense: 'max', explain: 'Largest vertical deflection relative to the rigid-body plane through the supports: the bending of the frame itself.' },
        { key: 'maxVM', label: 'Highest equivalent stress (envelope)', unit: 'MPa', value: maxVM / 1e6, explain: `Largest von Mises stress in any member, combining the extreme axial, bending, shear and torsion values along the member as if they coincided (conservative envelope). Compare with yield ${fmt(steel.fy / 1e6)} MPa; the member utilisation is the governing code check.` },
        { key: 'reactV', label: 'Total vertical reaction', unit: 'kN', value: reactV / 1e3, explain: 'Sum of vertical support or soil-spring reactions; it must equal the applied vertical load (equilibrium check).' },
        { key: 'baseShear', label: 'Total horizontal reaction', unit: 'kN', value: baseShear / 1e3, explain: 'Sum of horizontal reactions: the load passed to the foundation in sliding.' },
        { key: 'maxReact', label: 'Largest vertical support reaction', unit: 'kN', value: maxReact / 1e3, explain: 'Largest single support or spring reaction. For a mudmat divide by the tributary area to obtain the local contact pressure.' },
        { key: 'pSoil', label: 'Peak soil contact pressure', unit: 'kPa', value: pSoil / 1e3, explain: 'Largest spring reaction over its tributary area (soil-spring supports only): compare with the bearing capacity.' },
        { key: 'lambda', label: 'Elastic buckling load factor', unit: '-', value: lambda, limit: 3, sense: 'min', explain: 'Multiplier on the design loads at elastic frame buckling. 999 means no compressive buckling mode was found.' },
        { key: 'steelT', label: 'Frame steel mass', unit: 't', value: steelT, explain: 'Mass of the tubular members analysed (no plating, skirts or secondary steel).' },
      ],
      checks,
      plots: [
        { type: 'model3d', title: 'Frame with displacement and deformed shape', nodes: P.map((n) => [round(n[0]), round(n[1]), round(n[2])]), lines: ref.model.members.map((m) => [m.a, m.b]), nodeValues: disp.map((d) => round(d * 1e3)), valueLabel: 'Displacement', valueUnit: 'mm',
          deformed: P.map((n, i) => [round(n[0] + sc * u[6 * i]), round(n[1] + sc * u[6 * i + 1]), round(n[2] + sc * u[6 * i + 2])]),
          explain: `Nodal displacement magnitude on the frame, with the deformed shape exaggerated ${fmt(sc)} times. Look for sway of the posts and dishing of the base grid.` },
        { type: 'bar', title: 'Utilisation by member group', categories: groups.map((g) => g.group), series: [{ name: 'Utilisation', y: groups.map((g) => round(g.uc)) }], y: { label: 'Utilisation', unit: '-' },
          explain: 'Highest utilisation in each member group. An uneven picture shows where steel can be moved from one group to another.' },
        { type: 'bar', title: 'Vertical support reactions', categories: reacts.map((x) => `n${x.node}`), series: [{ name: 'Rz (kN)', y: reacts.map((x) => round(x.Rz / 1e3)) }], y: { label: 'Reaction', unit: 'kN' },
          explain: 'Vertical reaction at each support or soil spring. Negative values mean tension, which a mudmat on soil cannot provide: check lift-off in the soil-interaction module.' },
      ],
      tables: [
        { title: 'Member groups', columns: [{ key: 'group', label: 'Group' }, { key: 'n', label: 'Members' }, { key: 'N', label: 'Axial (governing member)', unit: 'kN' }, { key: 'M', label: 'Moment', unit: 'kN·m' }, { key: 'vm', label: 'von Mises', unit: 'MPa' }, { key: 'uc', label: 'Utilisation' }, { key: 'ref', label: 'Governing check' }],
          rows: groups.map((g) => ({ group: g.group, n: g.n, N: round(g.N / 1e3), M: round(g.M / 1e3), vm: round(g.vm / 1e6), uc: round(g.uc), ref: g.ref })),
          explain: 'One row per member group with the forces in its most utilised member (tension positive) and the code clause that governs.' },
        { title: 'Most utilised members', columns: [{ key: 'i', label: 'Member' }, { key: 'group', label: 'Group' }, { key: 'L', label: 'Length', unit: 'm' }, { key: 'N', label: 'Axial', unit: 'kN' }, { key: 'M', label: 'Moment', unit: 'kN·m' }, { key: 'V', label: 'Shear', unit: 'kN' }, { key: 'uc', label: 'Utilisation' }],
          rows: rows.slice().sort((a, b) => b.uc - a.uc).slice(0, 12).map((x) => ({ i: x.i, group: x.group, L: round(x.L), N: round(x.N / 1e3), M: round(x.M / 1e3), V: round(x.V / 1e3), uc: round(x.uc) })),
          explain: 'The twelve most utilised members with their extreme section forces, for targeted strengthening or a detailed joint check.' },
      ],
      notes: [`${gov.group} member ${gov.i} governs (${gov.ref}).`, `Vertical equilibrium: reactions ${fmt(reactV / 1e3)} kN.`, inp.support === 'springs' ? `Soil springs: ${fmt(inp.ks)} kPa/m on ${fmt(trib)} m² per node.` : `Supports: ${inp.support} at ${supNodes.length} nodes.`,
        reacts.some((x) => x.Rz < -1) && inp.support === 'springs' ? 'Some soil springs are in tension: the linear model overstates restraint there.' : 'No support is in tension.'],
      advice,
      publish: { 'subsea.frameUc': maxUc, 'subsea.frameSteel': steelT * 1e3, 'quantities.steel_t': steelT },
      quantities: { steel_t: round(steelT) },
    };
  },
  verify: [
    { name: 'Statics of a four-legged frame: Rmax = W/4 + H·h/(2·Lx), ΣRz = W, ΣRx = H', source: 'Global force and moment equilibrium with symmetry about the x–z plane: 800/4 + 100·3.5/(2·10) = 217.5 kN',
      inputs: { support: 'pinned', Lx: 10, Ly: 14, hT: 3.5, Weq: 800, Wroof: 0, Hx: 100, Fhub: 0, Vhub: 0, swFac: 0, gammaF: 1 },
      expect: [{ kpi: 'maxReact', value: 217.5, tol: 1e-5 }, { kpi: 'reactV', value: 800, tol: 1e-6 }, { kpi: 'baseShear', value: 100, tol: 1e-6 }] },
    { name: 'Uniform settlement on soil springs: δ = W / (ks·A_tributary·n)', source: 'Vertical equilibrium of four equal corner springs, each 400 kPa/m × 35 m² = 14 000 kN/m, loaded by 140 kN: δ = 10 mm, uniform, so the frame carries no bending',
      inputs: { support: 'springs', nx: 1, ny: 1, Lx: 10, Ly: 14, ks: 400, Weq: 560, Wroof: 0, Hx: 0, Fhub: 0, Vhub: 0, swFac: 0, gammaF: 1, nSub: 1 },
      expect: [{ kpi: 'maxDisp', value: 560 / (400 * 140) * 1e3, tol: 1e-4 }, { kpi: 'reactV', value: 560, tol: 1e-6 }] },
  ],
};
