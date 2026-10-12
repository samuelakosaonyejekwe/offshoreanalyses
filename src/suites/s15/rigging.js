// Lift rigging: four-point statically indeterminate sling distribution, skew load, sling / shackle / spreader /
// padeye checks and the lifted object as a beam between its lift points.
import { solveDense as solve } from '../../engine/math/linalg.js';
import { G, DEG, need, round, fmt, kpi, check, adv, bar, line, tbl, worst, dafN001, liftWeightInputs, chiEC3, tube, padeyeStress, steelCo2Text, CONFIRM, N001, LIMIT_N001 } from './_common.js';

export default {
  id: 's15-rigging', suite: 15, title: 'Lift rigging and sling loads',
  summary: 'Distributes the hook load to four slings with the centre of gravity offset and sling length misfit, then checks slings, shackles, spreader bar, padeyes and the lifted object between its lift points.',
  covers: [2], fidelity: 'analytical',
  standards: [`${N001} Sec. 16 (skew load, CoG envelope, sling and shackle factors, lift-point consequence factors)`, 'EN 1993-1-1 (2005) §6.3.1 (spreader buckling)', 'EN 13414 / ISO 2408 (wire rope slings, minimum breaking load)', 'EN 13889 (shackles, working load limit)'],
  models: ['Rigid-body equilibrium', 'Sling force distribution (four-point, one redundancy)', 'Nonlinear cable and sling stiffness (linearised axial stiffness EA/L)', 'Elastic spreader-beam analysis', 'Lifting-point stress assessment', 'Beam bending of the lifted object'],
  limits: `${LIMIT_N001} The lifted object and the hook are rigid, slings are linear axial springs and the compatibility condition is written for small displacements, so the skew result is a first-order estimate; padeye stresses are plain section statics without the weld, cheek-plate and out-of-plane checks of a detailed design, and local finite-element analysis of lift points is outside the module.`,
  publishes: ['install.slingLoadMax', 'install.liftPointLoad', 'install.skewFactor'],
  searchTerms: 'four point lift sling load distribution skew load factor padeye spreader bar DNV-ST-N001',
  inputs: [
    ...liftWeightInputs({ W: 1200, rig: 60 }),
    { key: 'a', label: 'Lift-point spacing along x', unit: 'm', type: 'number', default: 30, min: 1, max: 200, group: 'Lift points', role: 'geometry', design: true, help: 'Four lift points at the corners of an a × b rectangle; numbered 1 (−x, −y), 2 (+x, −y), 3 (+x, +y), 4 (−x, +y).' },
    { key: 'b', label: 'Lift-point spacing along y', unit: 'm', type: 'number', default: 20, min: 1, max: 200, group: 'Lift points', role: 'geometry', design: true },
    { key: 'ex', label: 'CoG offset from the lift-point centre, x', unit: 'm', type: 'number', default: 1.5, min: -50, max: 50, group: 'Lifted object', role: 'ic', uncertain: { dist: 'normal', cov: 0.3 }, help: 'Initial condition of the suspended object: the hook sits vertically above the centre of gravity.' },
    { key: 'ey', label: 'CoG offset from the lift-point centre, y', unit: 'm', type: 'number', default: 0.8, min: -50, max: 50, group: 'Lifted object', role: 'ic', uncertain: { dist: 'normal', cov: 0.3 } },
    { key: 'fCog', label: 'CoG inaccuracy factor', unit: '-', type: 'number', default: 1.05, min: 1, max: 1.3, group: 'Criteria', role: 'criteria', help: `Applied to sling and lift-point loads when no CoG envelope is analysed; 1.0 if the offsets above already are the envelope extremes. ${CONFIRM}` },
    { key: 'hHook', label: 'Hook (or spreader) height above the lift-point plane', unit: 'm', type: 'number', default: 36, min: 1, max: 200, group: 'Rigging', role: 'bc', design: true, help: 'Boundary condition: all four slings meet at the hook. Sets the sling angle; DNV-ST-N001 asks for at least about 60° to the horizontal unless justified.' },
    { key: 'dSling', label: 'Sling diameter', unit: 'mm', type: 'number', default: 220, min: 10, max: 600, group: 'Rigging', role: 'geometry', design: true },
    { key: 'mblCoef', label: 'Sling breaking-load coefficient (MBL = c·d²)', unit: 'MPa', type: 'number', default: 500, min: 100, max: 1200, group: 'Rigging', role: 'material', help: 'Minimum breaking load per unit of nominal diameter squared: about 450–630 MPa for cable-laid and six-strand steel slings, 550–600 MPa for HMPE. Replace with the sling certificate value.' },
    { key: 'eSling', label: 'Sling apparent modulus (on the nominal area πd²/4)', unit: 'GPa', type: 'number', default: 25, min: 2, max: 200, group: 'Rigging', role: 'material', uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'About 20–35 GPa for cable-laid steel slings, up to 60 GPa for single-part wire rope; sets the redundancy split and the hoisting stiffness.' },
    { key: 'sklMode', label: 'Skew load', type: 'select', default: 'factor', group: 'Criteria', role: 'criteria', options: [{ value: 'factor', label: 'Skew load factor on the statically distributed loads' }, { value: 'misfit', label: 'Computed from the sling length misfit' }],
      help: 'A four-sling single-hook lift is statically indeterminate: unequal sling lengths shift load onto one diagonal.' },
    { key: 'skl', label: 'Skew load factor SKL', unit: '-', type: 'number', default: 1.25, min: 1, max: 2, group: 'Criteria', role: 'criteria', showIf: { key: 'sklMode', in: ['factor'] }, help: `1.25 for an indeterminate four-sling lift with matched slings, 1.0 for a determinate arrangement, about 1.1 with a floating spreader. ${CONFIRM}` },
    { key: 'misfit', label: 'Length difference between the two diagonal sling pairs', unit: '%', type: 'number', default: 0.25, min: 0, max: 2, group: 'Criteria', role: 'criteria', showIf: { key: 'sklMode', in: ['misfit'] }, help: `Fabrication tolerance as a fraction of the sling length: one diagonal pair shorter than the other by this amount (about 0.25 % for matched slings). ${CONFIRM}` },
    { key: 'gSling', label: 'Sling nominal safety factor on MBL', unit: '-', type: 'number', default: 3.0, min: 1.5, max: 8, group: 'Criteria', role: 'criteria', help: `Product of the load, consequence, termination/bending reduction, wear and material factors of DNV-ST-N001 Sec. 16; about 3.0 for steel slings with good terminations. ${CONFIRM}` },
    { key: 'wllShackle', label: 'Shackle working load limit', unit: 't', type: 'number', default: 600, min: 1, max: 3000, group: 'Rigging', role: 'material', design: true },
    { key: 'spreader', label: 'Spreader bar between the long-side pairs', type: 'bool', default: false, group: 'Spreader', role: 'bc', help: 'When on, two spreader bars of length b take the transverse sling component so that the slings below are vertical in the y–z plane.' },
    { key: 'Dsp', label: 'Spreader tube diameter', unit: 'mm', type: 'number', default: 914, min: 100, max: 3000, group: 'Spreader', role: 'geometry', design: true, showIf: { key: 'spreader', in: [true] } },
    { key: 'tsp', label: 'Spreader wall thickness', unit: 'mm', type: 'number', default: 25, min: 4, max: 120, group: 'Spreader', role: 'geometry', design: true, showIf: { key: 'spreader', in: [true] } },
    { key: 'mat', label: 'Steel of padeyes and spreader', type: 'select', source: 'materials:steel', default: 's355', group: 'Lift points', role: 'material' },
    { key: 'tPad', label: 'Padeye total plate thickness at the hole', unit: 'mm', type: 'number', default: 200, min: 10, max: 600, group: 'Lift points', role: 'geometry', design: true, help: 'Main plate plus both cheek plates.' },
    { key: 'rPad', label: 'Padeye outer radius about the hole', unit: 'mm', type: 'number', default: 420, min: 30, max: 2000, group: 'Lift points', role: 'geometry', design: true },
    { key: 'dPin', label: 'Shackle pin diameter', unit: 'mm', type: 'number', default: 250, min: 10, max: 800, group: 'Lift points', role: 'geometry' },
    { key: 'holeGap', label: 'Hole diameter minus pin diameter', unit: 'mm', type: 'number', default: 6, min: 0, max: 60, group: 'Lift points', role: 'geometry', help: 'Keep within about 4 % of the pin diameter so that the bearing stress stays as assumed.' },
    { key: 'fLoad', label: 'Load factor for structural lift-point checks', unit: '-', type: 'number', default: 1.3, min: 1, max: 2, group: 'Criteria', role: 'criteria', help: CONFIRM },
    { key: 'fCons', label: 'Consequence factor for lift points and spreaders', unit: '-', type: 'number', default: 1.3, min: 1, max: 2, group: 'Criteria', role: 'criteria', help: `1.30 for lift points, attachments and spreaders, 1.15 for members framing directly into the lift point, 1.0 elsewhere. ${CONFIRM}` },
    { key: 'gM', label: 'Material factor', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.5, group: 'Criteria', role: 'criteria', help: CONFIRM },
    { key: 'Zobj', label: 'Section modulus of the object between lift points (x direction)', unit: 'm³', type: 'number', default: 0.9, min: 1e-4, max: 500, group: 'Lifted object', role: 'geometry', help: 'Both main girders together; used for the global bending stress of the object hanging on its lift points.' },
    { key: 'Iobj', label: 'Second moment of area of the object (x direction)', unit: 'm⁴', type: 'number', default: 1.4, min: 1e-5, max: 5000, group: 'Lifted object', role: 'geometry' },
    { key: 'Lobj', label: 'Object length (x direction)', unit: 'm', type: 'number', default: 44, min: 1, max: 400, group: 'Lifted object', role: 'geometry', help: 'Overall length; the parts beyond the lift points cantilever.' },
    { key: 'zCog', label: 'CoG height above the lift-point plane', unit: 'm', type: 'number', default: -9, min: -150, max: 50, group: 'Lifted object', role: 'ic', help: 'Negative when the lift points are above the centre of gravity (stable). Positive values need a stability check against the sling geometry.' },
  ],
  compute(inp, ctx) {
    // t → kg, mm → m, GPa/MPa → Pa, % → fraction.
    const steel = ctx.material(inp.mat); need(steel, `Steel '${inp.mat}' is not in the materials library.`);
    const fy = steel.fy;
    const Wc = inp.W * inp.fWeight, SHL = Wc + inp.Wrig;
    const daf = inp.dafUser > 0 ? inp.dafUser : dafN001(SHL, inp.where);
    const a = inp.a, b = inp.b, H = inp.hHook;
    need(Math.abs(inp.ex) < a / 2 && Math.abs(inp.ey) < b / 2, 'The centre of gravity lies outside the lift-point rectangle: the object would capsize on the hook. Move the lift points or add a counterweight.');
    const d = inp.dSling / 1e3, EA = inp.eSling * 1e9 * Math.PI * d * d / 4, MBL = inp.mblCoef * 1e6 * d * d;
    const F = Wc * 1e3 * G * daf; // dynamic object weight carried by the four slings, N
    // Lift points relative to the centre; hook above the CoG at (ex, ey, H).
    const P = [[-a / 2, -b / 2], [a / 2, -b / 2], [a / 2, b / 2], [-a / 2, b / 2]];
    const geo = P.map(([x, y]) => {
      // With spreader bars along y the slings above the bar carry the transverse component, those below are in the x–z plane.
      const dy = inp.spreader ? 0 : y - inp.ey, dx = x - inp.ex;
      const L = Math.hypot(dx, dy, H);
      return { x, y, L, sin: H / L, hor: Math.hypot(dx, dy) };
    });
    // Unknowns: vertical sling components V1..V4. Three equilibrium equations + one compatibility equation:
    // rigid object and rigid hook → the vertical lift-point displacements are coplanar: δ1 − δ2 + δ3 − δ4 = 0,
    // δi = Vi·Li/(EA·sin²θi) − mi/sinθi with mi the length misfit (short positive) of sling i.
    const f = geo.map((g) => g.L / (EA * g.sin * g.sin));
    const mis = inp.sklMode === 'misfit' ? inp.misfit / 100 : 0;
    const sgn = [1, -1, 1, -1];
    const solveV = (m) => {
      const A = [[1, 1, 1, 1], geo.map((g) => g.x), geo.map((g) => g.y), f.map((fi, i) => sgn[i] * fi)];
      const rhs = [F, F * inp.ex, F * inp.ey, geo.reduce((s, g, i) => s + sgn[i] * (sgn[i] * 0.5 * m * g.L) / g.sin, 0)];
      return Array.from(solve(A, rhs));
    };
    const V0 = solveV(0);
    // Misfit on either diagonal: take the worse sign for every sling.
    const Vp = mis > 0 ? solveV(mis) : V0, Vm = mis > 0 ? solveV(-mis) : V0;
    const skl = inp.sklMode === 'factor' ? inp.skl : 1;
    const T0 = V0.map((v, i) => v / geo[i].sin);
    const T = geo.map((g, i) => Math.max(Vp[i], Vm[i]) / g.sin * skl * inp.fCog);
    need(Math.min(...Vp, ...Vm) > 0, 'One sling goes slack (negative load): the centre of gravity is too close to an edge or the sling misfit is too large for this stiffness. Lengthen the slings, tighten the length tolerance or move the lift points.');
    const Tmax = Math.max(...T), iMax = T.indexOf(Tmax);
    const sklEff = Tmax / (T0[iMax] * inp.fCog);
    const angMin = Math.min(...geo.map((g) => Math.asin(g.sin) / DEG));
    const checks = [
      check('Sling: dynamic load against MBL / safety factor', Tmax / 1e3, MBL / inp.gSling / 1e3, 'kN', `${N001} Sec. 16 (sling design load and nominal safety factor)`, 'Guards against sling rupture. The breaking load is divided by the nominal safety factor that covers load, consequence, termination, bending and wear.'),
      check('Shackle: static sling load against WLL', Tmax / daf / 1e3, inp.wllShackle * G, 'kN', `${N001} Sec. 16 (shackles: WLL not less than the static sling load)`, 'Shackles are rated by working load limit, which already contains their own safety factor; the static sling load must not exceed it.'),
      check('Sling angle to the horizontal', 60, angMin, 'deg (required / achieved)', `${N001} Sec. 16 (minimum sling angle)`, 'Flat slings multiply the sling tension and put large compression into the object between the lift points.'),
    ];
    // Padeye, LRFD: design load = load factor × consequence factor × dynamic sling load, resistance fy/γM.
    const Fd = Tmax * inp.fLoad * inp.fCons;
    const ps = padeyeStress({ F: Fd, t: inp.tPad / 1e3, R: inp.rPad / 1e3, dHole: (inp.dPin + inp.holeGap) / 1e3, dPin: inp.dPin / 1e3 });
    need(2 * inp.rPad > inp.dPin + inp.holeGap, 'The padeye radius must exceed half the hole diameter.');
    const fd = fy / inp.gM;
    checks.push(
      check('Padeye: pin bearing', ps.bearing / 1e6, 1.5 * fd / 1e6, 'MPa', `${N001} Sec. 16 (lift points); EN 1993-1-8 §3.13 (pin bearing 1.5·fy/γM)`, 'Crushing of the plate under the pin. Exceedance ovalises the hole and overloads the cheek-plate welds.'),
      check('Padeye: tear-out shear', ps.tearOut / 1e6, fd / Math.sqrt(3) / 1e6, 'MPa', `${N001} Sec. 16; shear yield fy/(√3·γM) on two planes from hole to edge`, 'The pin tearing out through the plate edge on two shear planes: the classic padeye failure.'),
      check('Padeye: net-section tension', ps.netTension / 1e6, fd / 1e6, 'MPa', `${N001} Sec. 16; EN 1993-1-1 §6.2.3 (net section)`, 'Tension across the plate beside the hole. Exceedance splits the padeye at the hole.'));
    let Nsp = 0, steelSp = 0, chi = 1, NbRd = 0;
    if (inp.spreader) {
      // Each spreader bar (length b, pinned) carries the horizontal y-component of the two slings above it as compression.
      const tanY = (b / 2) / H;
      Nsp = Math.max(Vp[0] + Vp[3], Vp[1] + Vp[2], Vm[0] + Vm[3], Vm[1] + Vm[2]) / 2 * tanY * skl * inp.fCog * inp.fLoad * inp.fCons;
      const sec = tube(inp.Dsp / 1e3, inp.tsp / 1e3);
      const lam = (b / sec.r) / (Math.PI * Math.sqrt((steel.E || 2.1e11) / fy));
      chi = chiEC3(lam, 0.49); NbRd = chi * sec.A * fy / inp.gM; steelSp = 2 * sec.mass * b / 1e3;
      checks.push(check('Spreader bar: flexural buckling', Nsp / 1e3, NbRd / 1e3, 'kN', 'EN 1993-1-1 §6.3.1 (curve c for welded tubes); DNV-ST-N001 consequence factor for spreaders', 'Buckling of the spreader in compression. A spreader has no alternative load path: if it buckles the slings close in and the lift is lost.'));
    }
    // Object as a beam on two lift lines a apart with overhangs: uniformly distributed dynamic weight.
    const w = F / inp.Lobj, c = (inp.Lobj - a) / 2;
    need(c >= 0, 'The object length must be at least the lift-point spacing along x.');
    const Msup = w * c * c / 2, Mmid = w * a * a / 8 - Msup, Mobj = Math.max(Math.abs(Msup), Math.abs(Mmid));
    const sigObj = Mobj * inp.fLoad / inp.Zobj;
    const E = steel.E || 2.1e11;
    const defl = Math.abs(w * a * a / (384 * E * inp.Iobj) * (5 * a * a - 24 * c * c));
    // Horizontal sling components compress the object between the lift points along x.
    const Ncomp = Math.max(...geo.map((g, i) => T[i] * Math.abs(g.x - inp.ex) / g.L)) * 2;
    checks.push(check('Lifted object: global bending between lift points', sigObj / 1e6, fd / 1e6, 'MPa', `${N001} Sec. 16 (structural check of the lifted object, load factor on the dynamic weight)`, 'Bending of the object hanging on its lift points with its weight uniformly distributed. Exceedance yields the main girders during the lift.'));
    // Tilt if the hook is not adjusted: none (hook above CoG). Stability: CoG must lie below the sling apex.
    const stab = H - inp.zCog;
    checks.push(check('Suspended stability: CoG below the hook', Math.max(inp.zCog, 0), H, 'm', 'Rigid-body stability of a suspended load (metacentre at the sling apex)', 'A centre of gravity at or above the sling convergence point makes the load top-heavy: it can roll over on the hook.'));
    const gov = worst(checks), ok = gov.uc <= 1;
    // Sling diameter that just satisfies the sling check.
    const dReq = Math.sqrt(Tmax * inp.gSling / (inp.mblCoef * 1e6)) * 1e3;
    const tPadReq = inp.tPad * Math.max(ps.bearing / (1.5 * fd), ps.tearOut / (fd / Math.sqrt(3)), ps.netTension / fd);
    const padMass = (tPadReq - inp.tPad) / 1e3 * Math.PI * (inp.rPad / 1e3) ** 2 * 7.85 * 4; // t for four padeyes (disc approximation)
    const names = ['1 (−x, −y)', '2 (+x, −y)', '3 (+x, +y)', '4 (−x, +y)'];
    const hs = [0.6, 0.8, 1, 1.25, 1.5, 2].map((k) => k * H);
    return {
      kpis: [
        kpi('Tmax', 'Maximum dynamic sling load', Tmax / 1e3, 'kN', 'Largest sling tension including DAF, CoG factor and skew: the design load of sling, shackle and lift point.', MBL / inp.gSling / 1e3, 'max'),
        kpi('Tstat', 'Sling load from statics alone (no skew, no CoG factor)', T0[iMax] / 1e3, 'kN', 'Tension in the governing sling from equilibrium and elastic compatibility with perfectly matched slings.'),
        kpi('sklEff', 'Effective skew load factor', sklEff, '-', 'Ratio of the governing sling load with skew to the matched-sling value. With the misfit option it is computed from sling stiffness and tolerance; stiff short slings give high values.'),
        kpi('angle', 'Smallest sling angle to the horizontal', angMin, 'deg', 'Flatter slings carry more tension for the same vertical load. Below about 60° the horizontal component grows quickly.', 60, 'min'),
        kpi('MBLreq', 'Required sling minimum breaking load', Tmax * inp.gSling / 1e3, 'kN', 'Dynamic sling load times the nominal safety factor: the certificate value to order.'),
        kpi('dReq', 'Required sling diameter', dReq, 'mm', 'Diameter that just meets the breaking-load requirement with the coefficient entered.'),
        kpi('liftPoint', 'Lift-point design load', Fd / 1e3, 'kN', 'Dynamic sling load times load and consequence factors: the load for padeye, weld and local structure checks.'),
        kpi('ucPad', 'Padeye governing utilisation', Math.max(ps.bearing / (1.5 * fd), ps.tearOut / (fd / Math.sqrt(3)), ps.netTension / fd), '-', 'Largest of bearing, tear-out and net-section utilisations at the pin hole.', 1, 'max'),
        kpi('Nsp', 'Spreader bar design compression', Nsp / 1e3, 'kN', 'Factored axial force in each spreader bar from the transverse sling component; zero without a spreader.'),
        kpi('sigObj', 'Object global bending stress (factored)', sigObj / 1e6, 'MPa', 'Extreme-fibre stress of the object as a beam on its lift points under the dynamic weight.', fd / 1e6, 'max'),
        kpi('defl', 'Object mid-span deflection during the lift', defl * 1e3, 'mm', 'Elastic sag (or hog, with long overhangs) between the lift lines. Large values distort interfaces such as stabbing guides and piping.'),
        kpi('Ncomp', 'Compression induced in the object between lift lines', Ncomp / 1e3, 'kN', 'Horizontal sling components push the lift points together along x; the deck framing between them must carry this.'),
        kpi('stab', 'Hook height above the centre of gravity', stab, 'm', 'Distance from the sling apex down to the CoG. It must be positive, and comfortably so, for the load to hang stably.', 0, 'min'),
      ],
      checks,
      plots: [
        bar('Sling loads', names, 'Tension', 'kN', [{ name: 'Statics, matched slings', y: T0.map((v) => v / 1e3) }, { name: 'Design (skew + CoG factor)', y: T.map((v) => v / 1e3) }, { name: 'MBL / safety factor', y: T.map(() => MBL / inp.gSling / 1e3) }],
          'Tension in each of the four slings. The spread between slings comes from the CoG offset; the gap between the two bars is the skew and CoG allowance.'),
        line('Governing sling load against hook height', 'Hook height above lift points', 'm', 'Sling tension', 'kN', [{ name: 'Design sling load', x: hs, y: hs.map((h) => Tmax * (Math.hypot(geo[iMax].hor, h) / h) / (geo[iMax].L / H) / 1e3) }, { name: 'MBL / safety factor', x: [hs[0], hs[hs.length - 1]], y: [MBL / inp.gSling / 1e3, MBL / inp.gSling / 1e3], dash: true }],
          'Longer rigging steepens the slings and lowers the tension, at the price of hook height. Read off the shortest rigging that keeps the sling within its capacity.'),
      ],
      tables: [tbl('Sling and lift-point results', [['s', 'Sling'], ['L', 'Length', 'm'], ['ang', 'Angle to horizontal', 'deg'], ['V', 'Vertical component (statics)', 'kN'], ['T0', 'Tension (statics)', 'kN'], ['T', 'Design tension', 'kN'], ['uc', 'Sling utilisation']],
        geo.map((g, i) => ({ s: names[i], L: round(g.L), ang: round(Math.asin(g.sin) / DEG), V: round(V0[i] / 1e3), T0: round(T0[i] / 1e3), T: round(T[i] / 1e3), uc: round(T[i] * inp.gSling / MBL, 3) })),
        'One row per sling: geometry, the elastic share of the dynamic weight and the design tension with skew and CoG factor. Use the design tension for sling, shackle and padeye selection.')],
      notes: [`DAF ${fmt(daf)}; sling ${iMax + 1} governs with ${fmt(Tmax / 1e3)} kN.`, inp.sklMode === 'misfit' ? `Skew computed from a ${inp.misfit} % length difference between the diagonal pairs: effective factor ${fmt(sklEff)}.` : `Skew load factor ${inp.skl} applied to the statically distributed loads.`, `"${gov.name}" governs at ${fmt(gov.uc)}.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `All rigging and lift-point checks pass; "${gov.name}" governs at ${fmt(gov.uc)}. Required sling MBL ${fmt(Tmax * inp.gSling / 1e3)} kN (about ${fmt(dReq)} mm).`,
          `Order slings as matched pairs within the length tolerance assumed, and confirm the CoG envelope: an extra ${fmt(0.02 * a)} m of offset along x raises the governing sling by about ${fmt(100 * 0.02 * a * 2 / a)} %.`, `${N001} Sec. 16`)
          : adv('engineering', 'fail', `"${gov.name}" is exceeded (${fmt(gov.uc)}).`,
            gov.name.startsWith('Sling:') ? `Increase the sling to at least ${fmt(Math.ceil(dReq / 5) * 5)} mm, or raise the hook to steepen the slings.` : gov.name.startsWith('Padeye') ? `Increase the padeye thickness at the hole to about ${fmt(Math.ceil(tPadReq / 5) * 5)} mm (cheek plates) or enlarge the radius.` : gov.name.startsWith('Shackle') ? `Use shackles of at least ${fmt(Math.ceil(Tmax / daf / G / 1e3 / 50) * 50)} t WLL.` : gov.name.startsWith('Sling angle') ? `Raise the hook to at least ${fmt(Math.max(...geo.map((g) => g.hor)) * Math.tan(60 * DEG))} m above the lift points or add a spreader.` : 'Strengthen the governing component or reduce the dynamic load by a lower limiting sea state.',
            `${N001} Sec. 16`),
        adv('economic', 'info', `Rigging is sized by ${fmt(Tmax / 1e3)} kN per sling; the effective skew allowance is ${fmt((sklEff - 1) * 100)} % on top of statics.`,
          inp.spreader ? `The two spreader bars weigh about ${fmt(steelSp)} t and are reusable: compare their hire against the heavier slings and padeyes needed without them.` : 'A determinate arrangement (three-point lift or a floating spreader) removes most of the skew allowance and can drop the sling and shackle size by one step; compare that against the extra rigging height and hire.',
          'Lift rigging is normally hired: cost follows sling MBL class and shackle WLL class (contractor rate cards); AACE RP 18R-97 for estimate class',
          tPadReq > inp.tPad ? `Thicker padeyes add about ${steelCo2Text(Math.max(padMass, 0))}` : `Lift points are cut off and scrapped after installation: four padeyes of this size are about ${steelCo2Text(inp.tPad / 1e3 * Math.PI * (inp.rPad / 1e3) ** 2 * 7.85 * 4)} Design them for reuse or recycling.`),
      ],
      publish: { 'install.slingLoadMax': Tmax, 'install.liftPointLoad': Fd, 'install.skewFactor': sklEff },
      quantities: inp.spreader ? { steel_t: round(steelSp) } : {},
    };
  },
  verify: [
    { name: 'Symmetric four-point lift: T = W·DAF / (4·sin θ)', source: 'Statics: CoG central, a = 30 m, b = 20 m, H = 36 m → L = √(15² + 10² + 36²) = 40.26 m, sin θ = 0.8942; W = 1260 t, DAF 1.15 → T = 1260e3·9.80665·1.15/(4·0.8942) = 3972.8 kN',
      inputs: { ex: 0, ey: 0, fCog: 1, skl: 1 }, expect: [{ kpi: 'Tmax', value: 1260e3 * 9.80665 * 1.15 / (4 * 36 / Math.hypot(15, 10, 36)) / 1e3, tol: 1e-9 }, { kpi: 'sklEff', value: 1, tol: 1e-9 }] },
    { name: 'Object beam without overhang: M = wL²/8, σ = M·γ/Z', source: 'Simply supported beam, Lobj = a = 30 m: M = F·a/8 with F = 1260e3·9.80665·1.15 N → σ = 1.3·F·30/8/0.9 = 76.97 MPa',
      inputs: { ex: 0, ey: 0, Lobj: 30 }, expect: [{ kpi: 'sigObj', value: 1.3 * 1260e3 * 9.80665 * 1.15 * 30 / 8 / 0.9 / 1e6, tol: 1e-9 }] },
    { name: 'Spreader compression: N = (V/2-pair)·tan of the transverse sling angle', source: 'Statics, central CoG: each bar carries (F/2)/2·(b/2)/H·factors = F/4·(10/36)·1.3·1.3 with skl = fCog = 1',
      inputs: { ex: 0, ey: 0, fCog: 1, skl: 1, spreader: true }, expect: [{ kpi: 'Nsp', value: 1260e3 * 9.80665 * 1.15 / 4 * (10 / 36) * 1.69 / 1e3, tol: 1e-9 }] },
  ],
};
