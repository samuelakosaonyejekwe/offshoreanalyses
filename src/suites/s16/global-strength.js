// Global structural and foundation strength: section forces, stresses, shell-buckling utilisation and deflection of
// the monopile and tower under factored rotor thrust, tower wind drag, wave and current loading and gravity, with
// the foundation reactions at the mudline.
import { tubularSectionStress } from '../../engine/aero/blade.js';
import { c202ShellCheck } from '../../engine/structural/shells.js';
import { wavePileLoad, cantilever, towerWind } from './_loads.js';
import { turbineInputs, getTurbine, supportInputs, supportSegments, segmentsMass, foundationInputs, foundationOf, site, adv, check, worst, sev, fmt, round, run, need, G, DEG,
  LIMIT_POLAR, LIMIT_CODE, CONFIRM, SUPPLIER, steelCo2Text, CO2_STEEL } from './_common.js';

export default {
  id: 's16-global-strength', suite: 16, title: 'Global strength of tower and monopile',
  summary: 'Bending moments, shear, axial force, stresses, shell-buckling utilisation and deflection along the monopile and tower for a factored combination of rotor thrust, wind, wave, current and gravity, with the mudline reactions.',
  covers: [5], fidelity: 'analytical',
  standards: ['DNV-ST-0126 (2021) Support structures for wind turbines §4 (ULS, load and material factors)', 'DNV-RP-C202 (2019) Buckling strength of shells §3', 'IEC 61400-3-1:2019 §7 (structural design)', 'IEC 61400-1:2019 Table 3 (partial safety factors for loads)', 'EN 1993-1-6:2007 (shell buckling, alternative method)', 'DNV-RP-C205 (2021) §6 (wave and current loads)'],
  models: ['Euler–Bernoulli cantilever: section forces by direct integration, deflection by curvature integration', 'Coupled lateral / rocking foundation springs', 'Morison equation with Airy kinematics and Wheeler stretching', 'MacCamy–Fuchs diffraction (effective inertia coefficient)', 'Power-law wind profile and drag on the tower', 'Tubular section stresses and von Mises criterion', 'DNV-RP-C202 shell buckling under combined axial, bending and shear stress'],
  limits: `Quasi-static beam analysis of an axisymmetric tubular structure: dynamic amplification enters through the user factor on thrust (use the coupled dynamics module for it), loads are co-directional, and local details (flanges, door opening, grouted or bolted connection, attachments) need separate checks. Column buckling and P-delta effects are not included. ${LIMIT_POLAR} Nonlinear shell finite elements are outside the module. ${LIMIT_CODE}`,
  publishes: ['owt.mudlineMomentUls', 'owt.mudlineShearUls', 'owt.mudlineAxial', 'quantities.steel_t'],
  searchTerms: 'offshore wind monopile tower ultimate limit state shell buckling DNV-RP-C202',
  inputs: [
    ...turbineInputs(),
    ...supportInputs(),
    { key: 'corr', label: 'Corrosion allowance on the monopile wall', unit: 'mm', type: 'number', default: 4, min: 0, max: 15, group: 'Support structure', role: 'material', help: 'Deducted from the monopile wall for strength (splash zone and submerged, design life × corrosion rate; DNV-RP-0416).' },
    { key: 'shellL', label: 'Shell length between ring stiffeners or flanges', unit: 'm', type: 'number', default: 20, min: 2, max: 60, group: 'Support structure', role: 'geometry', help: 'Buckling length l of the DNV-RP-C202 curvature parameter (can length between flanges, bulkheads or stiffeners).' },
    ...foundationInputs(),
    { key: 'thrustSrc', label: 'Rotor thrust', type: 'select', default: 'rotor', group: 'Loads', role: 'bc', options: [{ value: 'rotor', label: 'Peak steady thrust of the rotor model' }, { value: 'input', label: 'Entered value' }], help: 'Interface load applied at hub height.' },
    { key: 'thrust', label: 'Characteristic rotor thrust', unit: 'kN', type: 'number', default: 750, min: 0, max: 10000, group: 'Loads', role: 'bc', bind: 'owt.thrustRated', showIf: { key: 'thrustSrc', in: ['input'] }, uncertain: { dist: 'normal', cov: 0.1 } },
    { key: 'daf', label: 'Gust and dynamic amplification of thrust', unit: '-', type: 'number', default: 1.4, min: 1, max: 2.5, group: 'Loads', role: 'bc', uncertain: { dist: 'normal', cov: 0.08 }, help: 'Characteristic extreme operating thrust over steady peak thrust (turbulence and tower dynamics): 1.3–1.6 from load simulations.' },
    { key: 'Uhub', label: 'Hub-height wind speed for tower drag', unit: 'm/s', type: 'number', default: 20, min: 0, max: 80, group: 'Loads', role: 'bc', help: 'Gust speed acting on the tower with the rotor load (0 = none).' },
    site('Hs1', { group: 'Loads', label: 'Significant wave height combined with the thrust' }), site('Tp1', { group: 'Loads', label: 'Associated peak period' }), site('Uc1', { group: 'Loads', label: 'Surface current' }),
    { key: 'hFactor', label: 'Design wave height / Hs', unit: '-', type: 'number', default: 1.86, min: 1, max: 2, group: 'Loads', role: 'bc', help: 'Most probable maximum in a 3-hour sea state ≈ 1.86 Hs (Rayleigh).' },
    { key: 'Cd', label: 'Drag coefficient of the pile', unit: '-', type: 'number', default: 1.0, min: 0.5, max: 1.6, group: 'Loads', role: 'bc', uncertain: { dist: 'normal', cov: 0.15 }, help: 'DNV-RP-C205: 0.65 smooth, 1.05 rough (marine growth).' },
    { key: 'Cm', label: 'Inertia coefficient of the pile', unit: '-', type: 'number', default: 2.0, min: 1.2, max: 2.2, group: 'Loads', role: 'bc' },
    { key: 'stretching', label: 'Kinematics above still water', type: 'select', default: 'wheeler', group: 'Loads', role: 'bc', options: [{ value: 'wheeler', label: 'Wheeler stretching to the crest' }, { value: 'none', label: 'None (integrate to still water)' }] },
    { key: 'diffraction', label: 'MacCamy–Fuchs reduction of the inertia coefficient', type: 'bool', default: true, group: 'Loads', role: 'bc' },
    { key: 'gammaF', label: 'Load factor on environmental loads', unit: '-', type: 'number', default: 1.35, min: 1, max: 1.7, group: 'Acceptance', role: 'criteria', help: `IEC 61400-1 Table 3 / DNV-ST-0126: 1.35 normal, 1.1 abnormal. ${CONFIRM}` },
    { key: 'gammaG', label: 'Load factor on gravity', unit: '-', type: 'number', default: 1.1, min: 0.9, max: 1.35, group: 'Acceptance', role: 'criteria', help: `1.1 unfavourable, 0.9 favourable. ${CONFIRM}` },
    { key: 'gammaM', label: 'Material factor for yield', unit: '-', type: 'number', default: 1.1, min: 1, max: 1.3, group: 'Acceptance', role: 'criteria', help: `DNV-ST-0126: 1.10. The shell-buckling check uses the slenderness-dependent factor of DNV-RP-C202 (1.15–1.45). ${CONFIRM}` },
    { key: 'dtLimit', label: 'Maximum D/t of the monopile', unit: '-', type: 'number', default: 120, min: 40, max: 250, group: 'Acceptance', role: 'criteria', help: 'Fabrication, driving and handling limit on the slenderness of the pile cans (80–120 is usual).' },
    { key: 'dtTower', label: 'Maximum D/t of the tower', unit: '-', type: 'number', default: 250, min: 60, max: 400, group: 'Acceptance', role: 'criteria', help: 'Tower cans are not driven and reach D/t of 150–250; shell buckling is checked separately.' },
    { key: 'deflLimit', label: 'Allowable tower-top deflection / height', unit: '-', type: 'number', default: 0.0125, min: 0.002, max: 0.05, group: 'Acceptance', role: 'criteria', help: 'Serviceability limit on the characteristic deflection (1/80 by default; project-specific).' },
    { key: 'nStations', label: 'Stations along the structure', type: 'int', default: 80, min: 10, max: 640, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nStations', label: 'Stations', levels: [20, 40, 80, 160], qoi: 'topDefl' },
  compute(inp, ctx) {
    const tb = getTurbine(inp), rotor = tb.rotor;
    need(inp.hubHeight > rotor.tipRadius + 5, 'The hub height leaves no clearance under the blade tip.');
    const st = ctx.material(inp.steel);
    need(st && st.fy > 0 && st.E > 0, `Steel '${inp.steel}' is not in the materials library.`);
    need(inp.monoT - inp.corr > 10, 'The corrosion allowance leaves less than 10 mm of monopile wall.');
    const gross = supportSegments(inp, { zHub: inp.hubHeight, hubOffset: 2.4 * tb.scale, E: st.E, rhoSteel: st.rho ?? 7850 });
    const segs = supportSegments({ ...inp, monoT: inp.monoT - inp.corr }, { zHub: inp.hubHeight, hubOffset: 2.4 * tb.scale, E: st.E, rhoSteel: st.rho ?? 7850 });
    const T = (inp.thrustSrc === 'rotor' ? tb.thrustRated : inp.thrust * 1e3) * inp.daf;
    const Hd = Math.min(inp.hFactor * inp.Hs1, 0.78 * inp.depth);   // depth-limited: a higher wave has broken before it reaches the structure
    const wave = run(() => wavePileLoad({ H: Hd, T: inp.Tp1, depth: inp.depth, D: inp.monoD, Cd: inp.Cd, Cm: inp.Cm, Uc: inp.Uc1, n: Math.max(30, Math.round(inp.nStations / 2)), stretching: inp.stretching, diffraction: inp.diffraction }), 'Wave and current load');
    const wind = inp.Uhub > 0 ? towerWind({ segs, Uhub: inp.Uhub, zHub: inp.hubHeight }) : null;
    const fnd = foundationOf(inp);
    const sol = cantilever({ segs, n: inp.nStations, points: [{ z: inp.hubHeight, F: T, M: 0 }], distributed: [wave], wind, topWeight: rotor.rnaMass * G, foundation: fnd, E: st.E, depth: inp.depth });
    // hub thrust acts above the tower top: the cantilever routine carries loads above the top station as force + moment at the top
    const n = sol.z.length;
    const rows = sol.z.map((z, i) => {
      const Md = sol.M[i] * inp.gammaF, Vd = sol.V[i] * inp.gammaF, Nd = sol.N[i] * inp.gammaG;
      const s = tubularSectionStress({ D: sol.D[i], t: sol.t[i], N: Nd, M: Md, V: Vd });
      const ucY = s.vonMises / (st.fy / inp.gammaM);
      const b = c202ShellCheck({ r: (sol.D[i] - sol.t[i]) / 2, t: sol.t[i], l: inp.shellL, E: st.E, fy: st.fy, sigmaA: -s.sigmaAxial, sigmaM: -s.sigmaBending, tau: 0 });   // compression-side extreme fibre, where the beam shear stress vanishes
      return { z, D: sol.D[i], t: sol.t[i], Md, Vd, Nd, vm: s.vonMises, ucY, ucB: b.utilisation, lam: b.lambdaS, fksd: b.fksd, sj: b.sigmaJ };
    });
    const iY = rows.reduce((a, r, i) => (r.ucY > rows[a].ucY ? i : a), 0), iB = rows.reduce((a, r, i) => (r.ucB > rows[a].ucB ? i : a), 0);
    const iTB = rows.reduce((a, r, i) => (Math.abs(r.z - inp.transZ) < Math.abs(rows[a].z - inp.transZ) ? i : a), 0);
    const H = sol.zTop - sol.zMud, topDefl = sol.u[n - 1], rotMud = sol.slope[0] / DEG;
    const dtMono = inp.monoD / ((inp.monoT - inp.corr) / 1e3), dtTow = Math.max(inp.towerBaseD / (inp.towerBaseT / 1e3), inp.towerTopD / (inp.towerTopT / 1e3));
    const checks = [
      check(`Yield at the governing section (z = ${fmt(rows[iY].z)} m)`, rows[iY].vm / 1e6, st.fy / inp.gammaM / 1e6, 'MPa', 'DNV-ST-0126 (2021) §4.4 (ULS, von Mises, γM = 1.10)', 'Factored von Mises stress at the extreme fibre against design yield. Exceeding it means first yield of the shell under the design load combination.'),
      check(`Shell buckling at the governing section (z = ${fmt(rows[iB].z)} m)`, rows[iB].sj / 1e6, rows[iB].fksd / 1e6, 'MPa', 'DNV-RP-C202 (2019) §3.2 (unstiffened cylindrical shell, combined loading)', 'Equivalent design stress against the design shell-buckling strength. Thin cans buckle locally on the compression side before yield; the reserve after buckling is small.'),
      check('Yield at the mudline', rows[0].vm / 1e6, st.fy / inp.gammaM / 1e6, 'MPa', 'DNV-ST-0126 (2021) §4.4', 'Stress at the mudline section. The maximum pile moment actually occurs a few metres below the mudline: add the embedded-pile check of the foundation module.'),
      check('Yield at the tower base', rows[iTB].vm / 1e6, st.fy / inp.gammaM / 1e6, 'MPa', 'DNV-ST-0126 (2021) §4.4', 'Stress at the tower-base / transition section, where the wall thickness steps down and the flange or grouted connection sits.'),
      check('Slenderness D/t of the monopile', dtMono, inp.dtLimit, '-', 'Fabrication and driving practice (API RP 2A: t ≥ 6.35 mm + D/100); DNV-ST-0126', 'Diameter-to-thickness ratio of the corroded pile wall against the project limit. Very slender cans are sensitive to buckling, ovalisation in handling and driving damage.'),
      check('Slenderness D/t of the tower', dtTow, inp.dtTower, '-', 'Tower fabrication practice; EN 1993-1-6 fabrication tolerance classes', 'Largest diameter-to-thickness ratio of the tower cans against the project limit. Slender cans need tighter fabrication tolerances for the assumed buckling strength.'),
      check('Tower-top deflection (characteristic loads)', Math.abs(topDefl), inp.deflLimit * H, 'm', 'Serviceability criterion (project); DNV-ST-0126 §4.6', 'Lateral deflection of the tower top under unfactored loads against the allowable fraction of the height. Excessive deflection adds P-delta moment and affects nacelle equipment.'),
    ];
    const gov = worst(checks);
    const massT = segmentsMass(gross) / 1e3;
    // wall thickness that brings the governing section to unity: section modulus ∝ t for thin tubes
    const ucGov = Math.max(rows[iY].ucY, rows[iB].ucB), sGov = rows[iY].ucY >= rows[iB].ucB ? rows[iY] : rows[iB];
    const tReq = sGov.t * 1e3 * (rows[iY].ucY >= rows[iB].ucB ? ucGov : ucGov ** 0.7);
    const monoShare = segmentsMass([gross[0]]) / 1e3;
    const dMass = monoShare * (Math.max(ucGov, 0.3) - 1) * 0.6;   // indicative: 60 % of the monopile mass follows the governing thickness
    const zr = (v) => round(v, 3);
    return {
      kpis: [
        { key: 'mudMoment', label: 'Design bending moment at the mudline', value: round(rows[0].Md / 1e6, 5), unit: 'MN·m', explain: 'Factored overturning moment where the pile enters the seabed: thrust times its lever arm plus wind on the tower and wave and current on the pile.' },
        { key: 'mudShear', label: 'Design shear at the mudline', value: round(rows[0].Vd / 1e6, 5), unit: 'MN', explain: 'Factored horizontal force transferred to the soil. Waves usually contribute most of the shear but, with their short lever arm, less of the moment.' },
        { key: 'mudAxial', label: 'Design axial force at the mudline', value: round(rows[0].Nd / 1e6, 5), unit: 'MN', explain: 'Factored weight of RNA, tower and pile above the mudline (buoyancy of the submerged steel deducted).' },
        { key: 'waveShear', label: 'Characteristic wave and current force on the pile', value: round(wave.shear / 1e3, 5), unit: 'kN', explain: `Unfactored hydrodynamic base shear at the phase of maximum mudline moment (design wave ${fmt(Hd)} m, inertia coefficient ${fmt(wave.cm)}).` },
        { key: 'ucYield', label: 'Highest yield utilisation', value: round(rows[iY].ucY, 4), unit: '-', limit: 1, sense: 'max', explain: `Factored von Mises stress over design yield at the worst section (z = ${fmt(rows[iY].z)} m).` },
        { key: 'ucBuckling', label: 'Highest shell-buckling utilisation', value: round(rows[iB].ucB, 4), unit: '-', limit: 1, sense: 'max', explain: `DNV-RP-C202 utilisation at the worst section (z = ${fmt(rows[iB].z)} m, reduced slenderness ${fmt(rows[iB].lam)}). Buckling governs over yield for slender tower cans.` },
        { key: 'sigmaMax', label: 'Maximum von Mises stress', value: round(rows[iY].vm / 1e6, 5), unit: 'MPa', explain: 'Largest factored equivalent stress in the shell, at the compression-side extreme fibre.' },
        { key: 'topDefl', label: 'Tower-top deflection (characteristic)', value: round(topDefl, 5), unit: 'm', limit: round(inp.deflLimit * H, 4), sense: 'max', explain: 'Lateral displacement of the tower top under unfactored loads, including foundation rotation and translation.' },
        { key: 'rotMud', label: 'Rotation at the mudline (characteristic)', value: round(rotMud, 5), unit: 'deg', explain: 'Tilt of the pile at the mudline from the foundation springs. Turbine suppliers limit the permanent tilt to about 0.25°–0.5° in total.' },
        { key: 'towerBaseM', label: 'Design bending moment at the tower base', value: round(rows[iTB].Md / 1e6, 5), unit: 'MN·m', explain: 'Factored moment at the transition: the interface load for the flange, grouted or slip-joint connection.' },
        { key: 'steel', label: 'Steel mass of monopile (to mudline) and tower', value: round(massT, 5), unit: 't', explain: 'Primary steel above the mudline with the tower mass factor. The embedded pile length is added in the foundation module.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Design bending moment and shear', x: { label: 'Moment (MN·m) · shear ×10 (MN)', unit: '-' }, y: { label: 'Elevation above still water', unit: 'm' }, series: [{ name: 'Bending moment (MN·m)', x: rows.map((r) => round(r.Md / 1e6, 4)), y: rows.map((r) => zr(r.z)) }, { name: 'Shear × 10 (MN)', x: rows.map((r) => round(r.Vd / 1e5, 4)), y: rows.map((r) => zr(r.z)) }],
          explain: 'Factored bending moment and shear along the height. The moment grows almost linearly from the hub down; the kink below still water is the wave load on the pile.' },
        { type: 'line', title: 'Utilisation along the structure', x: { label: 'Utilisation', unit: '-' }, y: { label: 'Elevation above still water', unit: 'm' }, series: [{ name: 'Yield', x: rows.map((r) => round(r.ucY, 4)), y: rows.map((r) => zr(r.z)) }, { name: 'Shell buckling', x: rows.map((r) => round(r.ucB, 4)), y: rows.map((r) => zr(r.z)) }, { name: 'Limit', x: [1, 1], y: [sol.zMud, sol.zTop], dash: true }],
          explain: 'Yield and shell-buckling utilisation against elevation. An efficient design runs close to a constant utilisation; peaks show where the wall should be thickened and troughs where steel can be saved.' },
        { type: 'line', title: 'Deflected shape (characteristic loads)', x: { label: 'Lateral deflection', unit: 'm' }, y: { label: 'Elevation above still water', unit: 'm' }, series: [{ name: 'Deflection', x: sol.u.map((v) => round(v, 5)), y: sol.z.map(zr) }],
          explain: 'Lateral displacement under unfactored loads. The offset and slope at the mudline come from the foundation flexibility; curvature is concentrated where the moment is high and the section small.' },
        { type: 'line', title: 'Lateral load distribution (characteristic)', x: { label: 'Load per unit height', unit: 'kN/m' }, y: { label: 'Elevation above still water', unit: 'm' }, series: [{ name: 'Wave and current on the pile', x: wave.f.map((v) => round(v / 1e3, 4)), y: wave.z.map(zr) }, ...(wind ? [{ name: 'Wind on the tower × 10', x: sol.z.filter((z) => z > 0).map((z) => round(wind(z) / 100, 4)), y: sol.z.filter((z) => z > 0).map(zr) }] : [])],
          explain: 'Hydrodynamic load per metre of pile at the governing wave phase and wind drag per metre of tower (scaled). The rotor thrust acts as a point load at the hub.' },
        { type: 'model3d', title: 'Structure with yield utilisation and deflected shape', nodes: sol.z.map((z) => [0, 0, zr(z)]), lines: sol.z.slice(1).map((_, i) => [i, i + 1]), tris: [], nodeValues: rows.map((r) => round(Math.max(r.ucY, r.ucB), 4)), valueLabel: 'Utilisation (yield or buckling)', valueUnit: '-', deformed: sol.z.map((z, i) => [round(sol.u[i] * 20, 4), 0, zr(z)]),
          explain: 'Stick model of monopile and tower coloured by the larger of the yield and buckling utilisation, with the deflected shape magnified 20 times.' },
      ],
      tables: [
        { title: 'Sections', columns: [{ key: 'z', label: 'Elevation', unit: 'm' }, { key: 'D', label: 'Diameter', unit: 'm' }, { key: 't', label: 'Wall (net)', unit: 'mm' }, { key: 'N', label: 'Axial', unit: 'MN' }, { key: 'V', label: 'Shear', unit: 'MN' }, { key: 'M', label: 'Moment', unit: 'MN·m' }, { key: 'vm', label: 'von Mises', unit: 'MPa' }, { key: 'ucY', label: 'UC yield', unit: '-' }, { key: 'ucB', label: 'UC buckling', unit: '-' }],
          rows: rows.filter((_, i) => i % Math.max(1, Math.round(rows.length / 40)) === 0 || i === rows.length - 1).map((r) => ({ z: zr(r.z), D: round(r.D, 4), t: round(r.t * 1e3, 4), N: round(r.Nd / 1e6, 4), V: round(r.Vd / 1e6, 4), M: round(r.Md / 1e6, 5), vm: round(r.vm / 1e6, 4), ucY: round(r.ucY, 4), ucB: round(r.ucB, 4) })),
          explain: 'Factored section forces, stress and utilisations at stations from the mudline to the tower top. Use it to grade the can thicknesses along the height.' },
        { title: 'Load summary (characteristic)', columns: [{ key: 'load', label: 'Load' }, { key: 'F', label: 'Force', unit: 'kN' }, { key: 'arm', label: 'Lever arm above mudline', unit: 'm' }, { key: 'M', label: 'Mudline moment', unit: 'MN·m' }],
          rows: (() => { const Fw = sol.V[0] - T - wave.shear, Mw = sol.M[0] - T * (inp.hubHeight + inp.depth) - wave.moment; return [{ load: 'Rotor thrust (incl. amplification)', F: round(T / 1e3, 5), arm: round(inp.hubHeight + inp.depth, 4), M: round(T * (inp.hubHeight + inp.depth) / 1e6, 5) }, { load: 'Wind on the tower', F: round(Fw / 1e3, 5), arm: round(Fw > 1 ? Mw / Fw : 0, 4), M: round(Mw / 1e6, 5) }, { load: 'Wave and current on the pile', F: round(wave.shear / 1e3, 5), arm: round(wave.shear > 1 ? wave.moment / wave.shear : 0, 4), M: round(wave.moment / 1e6, 5) }, { load: 'Total', F: round(sol.V[0] / 1e3, 5), arm: round(sol.M[0] / sol.V[0], 4), M: round(sol.M[0] / 1e6, 5) }]; })(),
          explain: 'Contribution of each load to the base shear and mudline moment before load factors. The thrust dominates the moment because of its lever arm; waves dominate the shear.' },
      ],
      notes: [
        `Loads: thrust ${fmt(T / 1e3)} kN at ${inp.hubHeight} m (steady ${fmt(T / inp.daf / 1e3)} kN × ${inp.daf}), design wave ${fmt(Hd)} m${Hd < inp.hFactor * inp.Hs1 - 1e-9 ? ' (depth-limited at 0.78 × depth)' : ''} / ${inp.Tp1} s with ${inp.Uc1} m/s current (ka = ${fmt(wave.ka)}, Cm = ${fmt(wave.cm)}), tower wind at ${inp.Uhub} m/s; γf = ${inp.gammaF}.`,
        `Governing check: ${gov.name}, utilisation ${fmt(gov.uc)}. Yield peaks at z = ${fmt(rows[iY].z)} m, shell buckling at z = ${fmt(rows[iB].z)} m.`,
        `Mudline (characteristic): displacement ${fmt(sol.u[0] * 1e3)} mm, rotation ${fmt(rotMud)}°; tower-top deflection ${fmt(topDefl)} m (${fmt(topDefl / H * 100)} % of the height).`,
        `Steel ${st.name}: fy = ${fmt(st.fy / 1e6)} MPa nominal — reduce for plate thickness above 16 mm per EN 10025.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `Highest utilisation ${fmt(ucGov)} (${rows[iY].ucY >= rows[iB].ucB ? 'yield' : 'shell buckling'} at z = ${fmt(sGov.z)} m, wall ${fmt(sGov.t * 1e3)} mm); design mudline moment ${fmt(rows[0].Md / 1e6)} MN·m and shear ${fmt(rows[0].Vd / 1e6)} MN.`,
          ucGov <= 1 && gov.uc > 1 ? `Strength is adequate but '${gov.name}' is exceeded (utilisation ${fmt(gov.uc)}): ${gov.name.includes('deflection') ? 'stiffen the structure or the foundation (larger diameter is the efficient lever)' : 'thicken the wall or reduce the diameter of that element'}.` : ucGov > 1 ? `Increase the wall at z = ${fmt(sGov.z)} m to about ${fmt(tReq, 3)} mm (or the diameter by ${fmt((Math.sqrt(ucGov) - 1) * 100)} %), then grade the neighbouring cans; a higher steel grade helps yield but not elastic buckling.` : `Strength is adequate. The governing wall could be reduced to about ${fmt(Math.max(tReq, sGov.D * 1e3 / (sGov.z > inp.transZ ? inp.dtTower : inp.dtLimit)), 3)} mm if fatigue, the frequency window and the D/t limit allow — check the fatigue and natural-frequency modules before removing steel. ${SUPPLIER}`,
          'DNV-ST-0126 (2021) §4 (ULS with γf = 1.35 and γM = 1.10) and DNV-RP-C202 (2019) §3.2 for shell buckling; latest editions known.',
          `Support structure above the mudline: ${fmt(massT)} t. ${steelCo2Text(massT)} Each millimetre on the monopile wall is about ${fmt(monoShare / inp.monoT)} t of steel ≈ ${fmt(monoShare / inp.monoT * CO2_STEEL)} t CO₂e.`),
        adv('economic', 'info', `The structure carries ${fmt(massT)} t of primary steel; bringing the governing section to a utilisation of 1.0 would change the monopile mass by roughly ${dMass >= 0 ? '+' : ''}${fmt(dMass)} t.`,
          'Optimise can by can: thickness transitions at fabrication-friendly steps (5 mm), fewer plate thicknesses and circumferential welds reduce fabrication hours; ultimate strength rarely governs alone — let fatigue and frequency confirm any saving.',
          'Class 4 estimate practice (AACE RP 18R-97): primary steel tonnage × fabricated unit rate drives the foundation CAPEX; whole-life costing per ISO 15686-5 for the corrosion allowance versus coating and cathodic protection.'),
      ],
      publish: { 'owt.mudlineMomentUls': rows[0].Md, 'owt.mudlineShearUls': rows[0].Vd, 'owt.mudlineAxial': rows[0].Nd, 'quantities.steel_t': massT },
      quantities: { steel_t: round(massT, 4) },
    };
  },
  verify: [
    { name: 'Clamped uniform cantilever under a hub force: tip deflection FL³/3EI + (F·e)L²/2EI and base moment F·(L + e)', source: 'Euler–Bernoulli beam theory (Roark, Formulas for Stress and Strain, Table 8.1 cases 1a and 3a)',
      // D = 6 m, t = 50 mm (no corrosion): EI = 210e9 × 4.13628 = 8.6862e11 N·m²; L = 5 + 87.6 = 92.6 m; e = 2.4 m (hub above the tower top)
      // deflection = 1e6·92.6³/(3·8.6862e11) + 1e6·2.4·92.6²/(2·8.6862e11) = 0.30471 + 0.01185 = 0.31656 m; base moment = 1 MN × 95 m = 95 MN·m (γf = 1)
      inputs: { foundation: 'fixed', depth: 5, monoD: 6, monoT: 50, corr: 0, towerBaseD: 6, towerTopD: 6, towerBaseT: 50, towerTopT: 50, hubHeight: 90, transZ: 10, thrustSrc: 'input', thrust: 1000, daf: 1, Uhub: 0, Hs1: 0.5, hFactor: 1, Uc1: 0, Cd: 0.5, Cm: 1.2, gammaF: 1, gammaG: 1, steel: 's355' },
      // (a 0.5 m wave on 5 m of water adds < 0.1 % to the moment)
      expect: [{ kpi: 'topDefl', value: 0.31656, tol: 0.005 }, { kpi: 'mudMoment', value: 95, tol: 0.005 }] },
    { name: 'Inertia-dominated Morison force on a pile: F = ρ·Cm·(πD²/4)·g·(H/2)·tanh(kd)', source: 'Morison et al. (1950); Dean & Dalrymple (1991) Water Wave Mechanics, §8.3 (closed-form pile force)',
      // D = 6 m, d = 30 m, H = 4 m, T = 8 s: k = 0.06543 rad/m (kd = 1.963, tanh = 0.9612); F_I = 1025·2·28.274·9.80665·2·0.9612 = 1092.7 kN.
      // Drag amplitude ≈ 70 kN is 90° out of phase and below F_I/2, so the maximum equals F_I; kinematics integrated to still water, no current.
      inputs: { depth: 30, monoD: 6, Hs1: 4, hFactor: 1, Tp1: 8, Uc1: 0, Cd: 1, Cm: 2, stretching: 'none', diffraction: false }, expect: [{ kpi: 'waveShear', value: 1092.7, tol: 0.01 }] },
  ],
};
