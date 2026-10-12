// Buckling and ultimate strength of the monopile and tower shells: DNV-RP-C202 shell buckling along the height,
// global (Euler) buckling of the cantilever by geometric-stiffness eigenvalue analysis, plastic section capacity and
// the Brazier ovalisation limit.
import { c202ShellCheck, brazierMoment } from '../../engine/structural/shells.js';
import { monopileTowerModel } from '../../engine/structural/builders.js';
import { bucklingAnalysis } from '../../engine/structural/framebuckling.js';
import { G, need, run, fmt, round, adv, check, worst, sev, steelCo2Text, tube, turbineInputs, supportInputs, foundationInputs, steelInputs, steelOf, turbine, foundationOf, steelMass, LIMIT_CODE, CONFIRM } from './_common.js';

export default {
  id: 's03-buckling-ultimate', suite: 3, title: 'Buckling and ultimate strength',
  summary: 'Shell buckling utilisation of every can of the monopile and tower under factored axial force, bending and shear, the global buckling load of the structure and the plastic and ovalisation limits of the sections.',
  covers: [9], fidelity: 'numerical',
  standards: ['DNV-RP-C202 (2019) §3.2–3.4 (unstiffened cylindrical shells)', 'DNV-ST-0126 (2021) §4.4–4.6', 'EN 1993-1-6:2007 (shell buckling, cross-reference)', 'EN 1993-1-1:2005 §5.5 and §6.2 (section classification, plastic resistance)', 'IEC 61400-6:2020 §6 (towers and foundations)'],
  models: ['Shell buckling of unstiffened cylinders under combined loading (DNV-RP-C202)', 'Geometric stiffness and linear buckling eigenvalue analysis of the beam model', 'Euler buckling', 'Plastic section modulus and von Mises equivalent stress', 'Brazier ovalisation limit moment', 'Force and moment equilibrium along the height', 'Timoshenko beam elements'],
  limits: `Unstiffened cans between effective supports (flanges, ring stiffeners); openings such as the tower door and J-tube holes, flange connections, the grouted or bolted transition and geometric imperfections beyond the code knock-down factors are not modelled. Section forces come from a two-force model (rotor thrust at the hub, the rest at an equivalent height). ${LIMIT_CODE} Nonlinear shell finite elements with imperfections are outside the module.`,
  searchTerms: 'monopile tower shell buckling DNV-RP-C202 global buckling plastic moment',
  inputs: [
    ...turbineInputs(),
    ...supportInputs({ embed: false }),
    ...steelInputs(),
    ...foundationInputs(),
    { key: 'M', label: 'Overturning moment at the mudline (characteristic)', unit: 'MN·m', type: 'number', default: 280, min: 0, max: 5000, group: 'Loads', role: 'bc', bind: 'owt.mudlineMoment', uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'H', label: 'Shear at the mudline (characteristic)', unit: 'MN', type: 'number', default: 4.9, min: 0.01, max: 100, group: 'Loads', role: 'bc', bind: 'owt.mudlineShear', uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'T', label: 'Horizontal load at hub height in that case', unit: 'kN', type: 'number', default: 1700, min: 0, max: 20000, group: 'Loads', role: 'bc',
      help: 'Rotor thrust or parked-rotor drag acting with the mudline loads. It sets the moment in the tower; the remainder of the shear is placed at the height that reproduces the mudline moment.' },
    { key: 'gammaF', label: 'Load factor on environmental loads', unit: '-', type: 'number', default: 1.35, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: `IEC 61400-3-1 normal cases 1.35, abnormal 1.1. ${CONFIRM}` },
    { key: 'gammaG', label: 'Load factor on weight', unit: '-', type: 'number', default: 1.1, min: 0.9, max: 1.5, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'lShell', label: 'Shell length between effective supports', unit: 'm', type: 'number', default: 20, min: 1, max: 60, group: 'Shell buckling', role: 'geometry', design: true, help: 'Distance between flanges or ring stiffeners; shorter panels buckle at higher stress.' },
    { key: 'gammaMglobal', label: 'Material factor for global buckling', unit: '-', type: 'number', default: 1.2, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: 'Applied to the elastic critical load; shell buckling factors follow DNV-RP-C202 §3.1.' },
    { key: 'nTower', label: 'Elements in the tower (and in the pile)', type: 'int', default: 12, min: 3, max: 96, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nTower', label: 'Tower elements', levels: [6, 12, 24, 48], qoi: 'Pcr' },
  compute(inp, ctx) {
    const steel = steelOf(inp, ctx);
    const tb = turbine(inp);
    const zTop = tb.zTowerTop, zMud = -inp.depth;
    need(zTop > inp.zInt + 5, 'The tower top must be above the interface level.');
    const M0 = inp.gammaF * inp.M * 1e6, H0 = inp.gammaF * inp.H * 1e6, T = Math.min(inp.gammaF * inp.T * 1e3, H0);
    const leverHub = inp.hubHeight - zMud;
    need(M0 >= T * leverHub * 0.999 || H0 - T < 1, `The mudline moment (${fmt(inp.M)} MN·m) is smaller than the hub load times its lever (${fmt(inp.T * leverHub / 1e3)} MN·m): reduce the hub load or correct the moment.`);
    // Remainder of the shear at the height above the mudline that reproduces the mudline moment.
    const Fw = Math.max(H0 - T, 0), hw = Fw > 1 ? Math.min(Math.max((M0 - T * leverHub) / Fw, 0), leverHub) : 0;
    const Mat = (z) => T * (inp.hubHeight - z) + Fw * Math.max(hw - (z - zMud), 0);
    const Vat = (z) => T + (z - zMud < hw ? Fw : 0);
    const secAt = (z) => (z < inp.zInt ? tube(inp.D, inp.t / 1e3) : tube(inp.towerDb + ((inp.towerDt - inp.towerDb) * (z - inp.zInt)) / (zTop - inp.zInt), (inp.towerTb + ((inp.towerTt - inp.towerTb) * (z - inp.zInt)) / (zTop - inp.zInt)) / 1e3));
    // Weight above each station (buoyancy of the flooded pile neglected: conservative for compression).
    const nS = 48, zs = [], rows = [];
    for (let i = 0; i <= nS; i++) zs.push(zMud + ((zTop - zMud) * i) / nS);
    const wAbove = new Array(nS + 1).fill(0);
    let acc = tb.rnaMass * G;
    for (let i = nS; i >= 0; i--) {
      if (i < nS) { const zm = 0.5 * (zs[i] + zs[i + 1]); acc += secAt(zm).A * steel.rho * (zm >= inp.zInt ? inp.towerMassFactor : 1) * G * (zs[i + 1] - zs[i]); }
      wAbove[i] = acc;
    }
    const eps2 = 235e6 / steel.fy;
    let wst = null;
    for (let i = 0; i <= nS; i++) {
      const z = zs[i] === inp.zInt ? zs[i] + 1e-6 : zs[i];
      const s = secAt(z), N = inp.gammaG * wAbove[i], Mz = Mat(z), V = Vat(z);
      const sa = -N / s.A, sm = -Mz / s.W, tau = (2 * V) / s.A;
      const c = run(() => c202ShellCheck({ r: (s.D - s.t) / 2, t: s.t, l: inp.lShell, E: steel.E, nu: steel.nu, fy: steel.fy, sigmaA: sa, sigmaM: sm, tau }), 'Shell buckling check');
      const vm = Math.sqrt((sa + sm) ** 2 + 3 * tau * tau);
      const dt = s.D / s.t, cls = dt <= 50 * eps2 ? 1 : dt <= 70 * eps2 ? 2 : dt <= 90 * eps2 ? 3 : 4;
      const row = { z, D: s.D, t: s.t, N, M: Mz, sig: Math.abs(sa + sm), vm, uc: c.utilisation, fksd: c.fksd, lam: c.lambdaS, fEa: c.elastic.fEa, cls, Mp: steel.fy * s.Zp, sec: s };
      rows.push(row);
      if (!wst || row.uc > wst.uc) wst = row;
    }
    // Global buckling: beam model on the foundation springs under a unit top load.
    const fnd = foundationOf(inp);
    const big = 1e13;
    const spring = fnd === 'fixed' ? null : [
      [fnd.Kh, 0, 0, 0, fnd.Khr, 0], [0, fnd.Kh, 0, -fnd.Khr, 0, 0], [0, 0, big, 0, 0, 0],
      [0, -fnd.Khr, 0, fnd.Kr, 0, 0], [fnd.Khr, 0, 0, 0, fnd.Kr, 0], [0, 0, 0, 0, 0, big]];
    const model = run(() => monopileTowerModel({ waterDepth: inp.depth, pileD: inp.D, pileT: inp.t / 1e3, interfaceZ: inp.zInt, hubHeight: zTop, towerBaseD: inp.towerDb, towerTopD: inp.towerDt, towerBaseT: inp.towerTb / 1e3, towerTopT: inp.towerTt / 1e3,
      rnaMass: tb.rnaMass, foundation: spring ? { mode: 'spring', spring } : { mode: 'fixed' }, nPile: Math.max(2, Math.round(inp.nTower * (inp.depth + inp.zInt) / (zTop - inp.zInt))), nTower: inp.nTower, material: { E: steel.E, nu: steel.nu, rho: steel.rho, fy: steel.fy } }), 'Beam model');
    const Pref = 1e6;
    model.loadCases.push({ name: 'top', nodal: [{ node: model.meta.topNode, F: [0, 0, -Pref, 0, 0, 0] }] });
    const bk = run(() => bucklingAnalysis(model, { loadCase: 'top', nModes: 2 }), 'Global buckling analysis');
    const Pcr = bk.factors[0] * Pref;
    need(Pcr > 0 && Number.isFinite(Pcr), 'The buckling analysis returned no positive critical load: check the foundation stiffness and section sizes.');
    const Ntot = inp.gammaG * wAbove[0];
    const base = rows[0];
    const braz = brazierMoment({ r: (base.D - base.t) / 2, t: base.t, E: steel.E, nu: steel.nu });
    const ucPl = M0 / (base.Mp / 1.1);
    const checks = [
      check(`Shell buckling, worst station at z = ${fmt(wst.z)} m`, wst.sig / 1e6 * (wst.uc / Math.max(wst.sig / wst.fksd, 1e-12)), wst.fksd / 1e6, 'MPa (equivalent stress)', 'DNV-RP-C202 (2019) §3.2.1: σj,Sd ≤ fksd', 'Guards against local buckling of the shell wall under axial compression, bending and shear. A buckle forms on the compression side, the section loses its bending capacity and the structure collapses without warning.'),
      check('Global buckling under factored weight', Ntot / 1e6, Pcr / inp.gammaMglobal / 1e6, 'MN', 'DNV-ST-0126 (2021) §4.5; EN 1993-1-1 §5.2.1 (αcr)', 'Guards against sway buckling of the whole cantilever under its own weight and the rotor–nacelle mass. A ratio near 1 means second-order (P–Δ) moments become large.'),
      check('Plastic bending resistance at the mudline', M0 / 1e6, base.Mp / 1.1 / 1e6, 'MN·m', 'EN 1993-1-1 §6.2.5 (applicable to class 1–2 sections only)', 'Upper bound of the bending capacity if the wall could yield fully before buckling. For slender shells the buckling check above governs; this value shows the reserve a thicker wall would unlock.'),
      check('Bending moment against the Brazier ovalisation limit at the mudline', M0 / 1e6, braz.momentLimit / 1e6, 'MN·m', 'Brazier (1927); Calladine, Theory of Shell Structures §16', 'Guards against flattening of a long thin tube in bending, which reduces stiffness until the section snaps. Relevant for very high D/t.'),
    ];
    checks[0].uc = wst.uc; checks[0].demand = round(wst.uc * wst.fksd / 1e6, 5);
    const gov = worst(checks), ok = gov.uc <= 1;
    // Thickness for unit shell-buckling utilisation at the worst station (stress ∝ 1/t and strength rises with t: iterate).
    let tReq = wst.t;
    for (let it = 0; it < 30; it++) {
      const s = tube(wst.D, tReq);
      const c = c202ShellCheck({ r: (s.D - s.t) / 2, t: s.t, l: inp.lShell, E: steel.E, nu: steel.nu, fy: steel.fy, sigmaA: -wst.N / s.A, sigmaM: -wst.M / s.W, tau: (2 * Vat(wst.z)) / s.A });
      const next = Math.min(Math.max(tReq * c.utilisation ** 0.7, 0.008), 0.45 * wst.D);
      if (Math.abs(next - tReq) < 1e-5) { tReq = next; break; }
      tReq = next;
    }
    const sm = steelMass({ ...inp, Lp: 0 }, tb, steel.rho);
    const dMassPerM = (tube(wst.D, tReq).A - wst.sec.A) * steel.rho / 1e3;
    return {
      kpis: [
        { key: 'ucShell', label: 'Shell buckling utilisation (worst can)', unit: '-', value: wst.uc, limit: 1, sense: 'max', explain: `Design equivalent stress over design shell-buckling strength at the most utilised station (z = ${fmt(wst.z)} m, D/t = ${fmt(wst.D / wst.t)}). Must not exceed 1.` },
        { key: 'zWorst', label: 'Elevation of the worst can', unit: 'm', value: wst.z, explain: 'Where shell buckling governs. Usually just above the mudline or at a thickness step in the tower.' },
        { key: 'fksd', label: 'Design buckling strength at the worst can', unit: 'MPa', value: wst.fksd / 1e6, explain: 'Characteristic shell-buckling strength divided by the slenderness-dependent material factor (1.15 to 1.45). Compare with yield: the gap is the cost of slenderness.' },
        { key: 'lambdaS', label: 'Reduced shell slenderness at the worst can', unit: '-', value: wst.lam, explain: 'Square root of yield over elastic buckling stress for the load mix. Below 0.5 the shell is stocky and yield governs; above 1 elastic buckling governs.' },
        { key: 'fEa', label: 'Elastic axial buckling stress at the mudline can', unit: 'MPa', value: base.fEa / 1e6, explain: 'Elastic critical stress of the shell in pure axial compression including the code knock-down for imperfections.' },
        { key: 'Pcr', label: 'Global elastic buckling load (top load)', unit: 'MN', value: Pcr / 1e6, explain: 'Vertical load at the tower top at which the whole structure on its foundation would buckle sideways (first eigenvalue of the geometric-stiffness problem).' },
        { key: 'alphaCr', label: 'Global buckling load factor', unit: '-', value: Pcr / Ntot, explain: 'Critical load over factored weight. Above 10 second-order effects are negligible; below about 4 a second-order analysis is required.' },
        { key: 'Mp', label: 'Plastic moment at the mudline', unit: 'MN·m', value: base.Mp / 1e6, explain: 'Bending moment at full plastification of the pile section, fy times the plastic modulus. Reachable only by compact (class 1–2) sections.' },
        { key: 'sectionClass', label: 'Section class at the mudline', unit: '-', value: base.cls, explain: 'EN 1993-1-1 classification by D/t: 1–2 can develop the plastic moment, 3 reaches yield only, 4 is governed by shell buckling.' },
        { key: 'vmMax', label: 'Largest von Mises stress', unit: 'MPa', value: Math.max(...rows.map((r) => r.vm)) / 1e6, limit: steel.fy / 1.1 / 1e6, sense: 'max', explain: 'Largest factored equivalent stress along the height from axial force, bending and shear, against design yield.' },
        { key: 'Ntot', label: 'Factored weight at the mudline', unit: 'MN', value: Ntot / 1e6, explain: 'Rotor–nacelle assembly, tower and monopile above the mudline with the weight factor (buoyancy neglected).' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Shell buckling utilisation along the height', x: { label: 'Utilisation', unit: '-' }, y: { label: 'Elevation above still water', unit: 'm' },
          series: [{ name: 'Shell buckling (DNV-RP-C202)', x: rows.map((r) => r.uc), y: rows.map((r) => r.z) }, { name: 'von Mises / design yield', x: rows.map((r) => r.vm / (steel.fy / 1.1)), y: rows.map((r) => r.z), dash: true }, { name: 'Limit', x: [1, 1], y: [zMud, zTop], dash: true }],
          explain: 'Utilisation of each can from the mudline to the tower top. Peaks above the limit line show where the wall must be thickened; a flat profile near the limit is an efficient thickness distribution.' },
        { type: 'line', title: 'Design stress and buckling strength', x: { label: 'Stress', unit: 'MPa' }, y: { label: 'Elevation above still water', unit: 'm' },
          series: [{ name: 'Axial + bending stress', x: rows.map((r) => r.sig / 1e6), y: rows.map((r) => r.z) }, { name: 'Design buckling strength fksd', x: rows.map((r) => r.fksd / 1e6), y: rows.map((r) => r.z) }, { name: 'Design yield', x: [steel.fy / 1.1 / 1e6, steel.fy / 1.1 / 1e6], y: [zMud, zTop], dash: true }],
          explain: 'Factored compressive stress against the buckling strength at each elevation. The distance between buckling strength and yield is the penalty paid for slender walls.' },
        { type: 'line', title: 'Factored bending moment and axial force', x: { label: 'Moment (MN·m) / axial force ×10 (MN)', unit: '' }, y: { label: 'Elevation above still water', unit: 'm' },
          series: [{ name: 'Bending moment, MN·m', x: rows.map((r) => r.M / 1e6), y: rows.map((r) => r.z) }, { name: 'Axial compression ×10, MN', x: rows.map((r) => r.N / 1e5), y: rows.map((r) => r.z) }],
          explain: 'Section forces used in the checks: moment from the hub load and the lower horizontal resultant, compression from the weight above each section.' },
      ],
      tables: [
        { title: 'Shell check by station', columns: [{ key: 'z', label: 'Elevation', unit: 'm' }, { key: 'D', label: 'D', unit: 'm' }, { key: 't', label: 't', unit: 'mm' }, { key: 'dt', label: 'D/t' }, { key: 'N', label: 'N', unit: 'MN' }, { key: 'M', label: 'M', unit: 'MN·m' }, { key: 's', label: 'σ', unit: 'MPa' }, { key: 'f', label: 'fksd', unit: 'MPa' }, { key: 'uc', label: 'Utilisation' }, { key: 'c', label: 'Class' }],
          rows: rows.filter((_, i) => i % 3 === 0 || i === nS).map((r) => ({ z: round(r.z), D: round(r.D), t: round(r.t * 1e3), dt: round(r.D / r.t), N: round(r.N / 1e6), M: round(r.M / 1e6), s: round(r.sig / 1e6), f: round(r.fksd / 1e6), uc: round(r.uc), c: r.cls })),
          explain: 'Section size, factored forces, stress, design buckling strength and utilisation at stations along the height: the working table for sizing can thicknesses.' },
      ],
      notes: [`${gov.name} governs (utilisation ${fmt(gov.uc)}).`, `Shell length between supports ${inp.lShell} m; reduced slenderness ${fmt(wst.lam)} at the worst can.`,
        `Horizontal loads: ${fmt(T / 1e6)} MN at the hub and ${fmt(Fw / 1e6)} MN at ${fmt(hw)} m above the mudline (factored).`, 'Door openings, flanges and the transition connection need separate detailed checks.'],
      advice: [
        ok ? adv('engineering', sev(gov.uc), `Buckling checks pass: shell ${fmt(wst.uc * 100)} % at z = ${fmt(wst.z)} m, global load factor ${fmt(Pcr / Ntot)}. Section class at the mudline is ${base.cls}.`,
          `The worst can would reach 100 % with a wall of about ${fmt(tReq * 1e3)} mm (now ${fmt(wst.t * 1e3)} mm); confirm fabrication tolerances meet DNV-RP-C202 §3.1 (out-of-roundness, misalignment) since the strength assumes them.`,
          'DNV-RP-C202 (2019) §3.1–3.4; DNV-ST-0126 (2021) §4.4',
          `Shell buckling leaves room to change the worst can by ${fmt(dMassPerM)} t per metre — ${steelCo2Text(dMassPerM)} Fatigue and frequency usually limit how much of this can be taken.`)
          : adv('engineering', 'fail', `${gov.name}: utilisation ${fmt(gov.uc)}.`,
            gov.name.includes('Shell') ? `Thicken the can at z = ${fmt(wst.z)} m from ${fmt(wst.t * 1e3)} to about ${fmt(tReq * 1e3 * 1.02)} mm, or add a ring stiffener to halve the shell length (${inp.lShell} m).` : gov.name.includes('Global') ? 'Stiffen the foundation or increase diameters: the critical load scales with EI and with the rocking stiffness.' : `Increase the wall or diameter at the mudline: plastic modulus ≈ D²t, required ${fmt(M0 * 1.1 / steel.fy)} m³.`,
            'DNV-RP-C202 (2019) §3.2; DNV-ST-0126 (2021) §4.5; EN 1993-1-1 §6.2.5',
            `The thicker can adds ${fmt(dMassPerM)} t per metre — ${steelCo2Text(dMassPerM)} A ring stiffener achieves a similar gain in buckling strength with a fraction of that steel.`),
        adv('economic', 'info', `Primary steel above the mudline is ${fmt(sm.total / 1e3)} t; the design buckling strength at the worst can is ${fmt(wst.fksd / (steel.fy / 1.1) * 100)} % of design yield.`,
          wst.fksd < 0.8 * steel.fy / 1.1 ? 'A higher steel grade will not pay here because elastic buckling, not yield, limits the wall: spend on thickness or stiffeners rather than grade.' : 'Yield largely governs, so a higher grade (S420/S460) converts directly into thinner cans; compare the plate premium with the saved tonnage and welding volume.',
          'Fabrication cost drivers for rolled and welded cans (plate tonnage, weld volume ∝ t²); AACE RP 18R-97 Class 4 quantity estimate'),
      ],
      quantities: { steel_t: round(sm.total / 1e3) },
    };
  },
  verify: [
    { name: 'Global buckling of a uniform clamped cantilever: Pcr = π²EI/(4L²)', source: 'Euler (Timoshenko & Gere, Theory of Elastic Stability, §2.1): EI = 1.0371e12 N·m² (D = 6 m, t = 60 mm, E = 210 GPa), L = 25 + 87.6 = 112.6 m → Pcr = 201.8 MN',
      inputs: { fnd: 'fixed', D: 6, t: 60, towerDb: 6, towerDt: 6, towerTb: 60, towerTt: 60, depth: 25, hubHeight: 90, ratedMW: 5, nTower: 24 }, expect: [{ kpi: 'Pcr', value: 201.8, tol: 0.01 }] },
    { name: 'Elastic axial shell buckling stress, DNV-RP-C202 Table 3-2', source: 'Hand calculation: r = 2.965 m, t = 70 mm, l = 20 m, ν = 0.3 → Zl = l²/(rt)·√(1−ν²) = 1838.5, ξ = 0.702 Zl = 1290.6, ρ = 0.5(1 + r/150t)^-0.5 = 0.4415, C = √(1 + (ρξ)²) = 569.85, fEa = C·π²E/(12(1−ν²))·(t/l)² = 1324.9 MPa',
      inputs: { D: 6, t: 70, lShell: 20 }, expect: [{ kpi: 'fEa', value: 1324.9, tol: 0.002 }] },
  ],
};
