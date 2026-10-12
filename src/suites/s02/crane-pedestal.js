// Crane pedestal and lifting support: crane hook load with dynamic factor, offlead, sidelead and wind taken to the
// pedestal base; tubular strength, shell buckling, slew-ring flange bolts, top rotation, natural frequency and
// girth-weld fatigue; and the allowable hook load against radius set by the pedestal.
import { tubeSection } from '../../engine/structural/sections.js';
import { tubularMemberCheck } from '../../engine/structural/tubular.js';
import { c202ShellCheck } from '../../engine/structural/shells.js';
import { boltGroupElastic, boltResistance, BOLT_STRESS_AREA, BOLT_GRADES } from '../../engine/structural/steelcheck.js';
import { scfGirthWeld } from '../../engine/structural/scf.js';
import { snCurve, snConstants, cyclesFromConstants, SN_CURVE_NAMES } from '../../engine/fatigue/sn.js';
import { steelOf, need, run, fmt, round, linspace, adv, check, fsCheck, steelCo2Text, CONFIRM, G, DEG } from './_common.js';

export default {
  id: 's02-crane-pedestal', suite: 2, title: 'Crane pedestal and lifting support',
  summary: 'Pedestal base forces from the crane load chart point, stresses and utilisation of the pedestal tube, shell buckling, slew-ring bolt loads, rotation at the slew bearing, natural frequency, fatigue of the base weld and the hook load the pedestal allows at each radius.',
  covers: [8], fidelity: 'analytical',
  standards: ['EN 13852-1:2013 (offshore cranes) §5, Annex B (dynamic factors, offlead and sidelead)', 'API Spec 2C (2020) §5 (pedestal design loads)', 'DNV-ST-0378 (2021) (offshore and platform lifting appliances)', 'ISO 19902:2020 §13 (tubular members)', 'DNV-RP-C202 (2019) §3.2 (shell buckling)', 'EN 1993-1-8:2005 Table 3.4 (bolts)', 'DNV-RP-C203 (2021) §3.3.7, §2.4 (girth-weld SCF, S–N curves)'],
  models: ['Static equilibrium (cantilever statics)', 'Elastic beam-column analysis', 'Classical shell theory (buckling of an unstiffened cylinder)', 'Elastic bolt-group analysis', 'Dynamic amplification approach (hoisting factor)', 'Structural fatigue model (S–N, constant equivalent range)', 'Euler–Bernoulli cantilever stiffness'],
  limits: 'Pedestal as a prismatic cantilever tube fixed at the deck with code and shell-buckling checks on nominal stresses: the stiffened deck connection, ring stiffeners, access openings and the slew bearing itself are not modelled (use the local finite-element stress module for the deck junction and openings). Bolt forces are elastic without prying or preload loss. The crane dynamic factor is an input from the crane standard or motion analysis, not derived here.',
  publishes: ['topsides.craneMoment'],
  searchTerms: 'offshore crane pedestal design overturning moment slew bearing bolts shell buckling fatigue',
  inputs: [
    { key: 'swl', label: 'Hook load (safe working load at this radius)', unit: 't', type: 'number', default: 50, min: 0.5, max: 5000, group: 'Crane loads', role: 'load', design: true },
    { key: 'radius', label: 'Working radius', unit: 'm', type: 'number', default: 30, min: 2, max: 150, group: 'Crane loads', role: 'geometry', help: 'Horizontal distance from the slew axis to the hook for this load chart point.' },
    { key: 'daf', label: 'Dynamic amplification factor on the hook load', unit: '-', type: 'number', default: 1.5, min: 1, max: 4, group: 'Crane loads', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 }, help: `EN 13852-1 / API 2C: 1.3 for platform lifts, 1.5–2.5 and more for lifts from a supply vessel depending on sea state and hoist speed. ${CONFIRM}` },
    { key: 'mHook', label: 'Hook block and rigging mass', unit: 't', type: 'number', default: 2, min: 0, max: 200, group: 'Crane loads', role: 'load' },
    { key: 'mCrane', label: 'Crane dead weight above the pedestal', unit: 't', type: 'number', default: 160, min: 0, max: 5000, group: 'Crane loads', role: 'load' },
    { key: 'eCrane', label: 'Crane centre of gravity from the slew axis (toward the hook)', unit: 'm', type: 'number', default: 3, min: -20, max: 60, group: 'Crane loads', role: 'geometry', help: 'Negative when the machinery house and counterweight outweigh the boom.' },
    { key: 'offlead', label: 'Offlead (radial) force as a share of the hook load', unit: '%', type: 'number', default: 6, min: 0, max: 30, group: 'Crane loads', role: 'bc', help: `From vessel motion and hook swing: 2–3 % platform lifts, 5–10 % sea lifts (EN 13852-1 Annex B). ${CONFIRM}` },
    { key: 'sidelead', label: 'Sidelead (tangential) force as a share of the hook load', unit: '%', type: 'number', default: 3, min: 0, max: 20, group: 'Crane loads', role: 'bc' },
    { key: 'hTip', label: 'Boom tip above the pedestal top', unit: 'm', type: 'number', default: 28, min: 1, max: 150, group: 'Crane loads', role: 'geometry' },
    { key: 'Fwind', label: 'Wind force on the crane', unit: 'kN', type: 'number', default: 60, min: 0, max: 5000, group: 'Crane loads', role: 'environment', help: 'In-service wind on boom and house acting toward the hook side.' },
    { key: 'hWind', label: 'Wind force above the pedestal top', unit: 'm', type: 'number', default: 15, min: 0, max: 100, group: 'Crane loads', role: 'geometry' },
    { key: 'D', label: 'Pedestal outer diameter', unit: 'mm', type: 'number', default: 2500, min: 300, max: 8000, group: 'Pedestal', role: 'geometry', design: true },
    { key: 't', label: 'Pedestal wall thickness', unit: 'mm', type: 'number', default: 45, min: 8, max: 150, group: 'Pedestal', role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'H', label: 'Pedestal height above the deck', unit: 'm', type: 'number', default: 10, min: 0.5, max: 60, group: 'Pedestal', role: 'geometry' },
    { key: 'base', label: 'Restraint at the deck', type: 'select', default: 'fixed', group: 'Pedestal', role: 'bc', options: [{ value: 'fixed', label: 'Fixed (through-deck pedestal with bulkheads)' }], help: 'The pedestal is a cantilever: it must be continuous through at least one deck with shear plates or bulkheads to be fixed.' },
    { key: 'mat', label: 'Steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Pedestal', role: 'material' },
    { key: 'nBolt', label: 'Slew-ring bolts', type: 'int', default: 84, min: 8, max: 300, group: 'Slew-ring flange', role: 'geometry' },
    { key: 'dBc', label: 'Bolt circle diameter', unit: 'mm', type: 'number', default: 2650, min: 300, max: 9000, group: 'Slew-ring flange', role: 'geometry' },
    { key: 'bolt', label: 'Bolt size', type: 'select', default: 'M36', group: 'Slew-ring flange', role: 'geometry', options: Object.keys(BOLT_STRESS_AREA).map((v) => ({ value: v, label: v })) },
    { key: 'boltGrade', label: 'Bolt grade', type: 'select', default: '10.9', group: 'Slew-ring flange', role: 'material', options: Object.keys(BOLT_GRADES).map((v) => ({ value: String(v), label: String(v) })) },
    { key: 'gF', label: 'Partial factor on crane loads', unit: '-', type: 'number', default: 1.3, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: `Applied to all crane actions (the dynamic factor is separate). API 2C pedestal factors and EN 13852-1 partial factors differ — ${CONFIRM}` },
    { key: 'rotLim', label: 'Rotation limit at the slew bearing', unit: 'deg', type: 'number', default: 0.5, min: 0.02, max: 3, group: 'Acceptance criteria', role: 'criteria', help: 'Bearing manufacturer limit on tilt of the supporting flange under working load (unfactored).' },
    { key: 'snName', label: 'S–N curve of the base girth weld', type: 'select', default: 'D', group: 'Fatigue', role: 'criteria', options: SN_CURVE_NAMES.map((v) => ({ value: v, label: `DNV-RP-C203 ${v}` })) },
    { key: 'misalign', label: 'Fabrication misalignment at the weld', unit: 'mm', type: 'number', default: 4, min: 0, max: 30, group: 'Fatigue', role: 'geometry' },
    { key: 'lifts', label: 'Lifts per year', unit: '-', type: 'number', default: 6000, min: 1, max: 1e6, group: 'Fatigue', role: 'operation' },
    { key: 'kLoad', label: 'Cubic-mean hook load over the rated load', unit: '-', type: 'number', default: 0.5, min: 0.05, max: 1, group: 'Fatigue', role: 'operation', help: 'Load spectrum factor: (mean of load³)^(1/3) / SWL. About 0.4–0.6 for general platform cranes, higher for dedicated heavy duty.' },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 3, min: 1, max: 10, group: 'Fatigue', role: 'criteria', help: CONFIRM },
    { key: 'life', label: 'Design life', unit: 'yr', type: 'number', default: 25, min: 1, max: 100, group: 'Fatigue', role: 'criteria', bind: 'life.designYears' },
  ],
  compute(inp, ctx) {
    // SI: t → kg, mm → m, kN → N, % → fraction.
    const st = steelOf(ctx, inp.mat), D = inp.D / 1e3, t = inp.t / 1e3, H = inp.H;
    need(t < D / 2, 'The wall thickness must be smaller than the pedestal radius.');
    need(inp.dBc / 1e3 > 0.5 * D, 'The bolt circle diameter must be comparable to the pedestal diameter.');
    const sec = tubeSection({ D, t }), EI = st.E * sec.Iy;
    // Pedestal forces for a hook load P [N] at radius R, with load factor gF (1 for serviceability and fatigue).
    const forces = (swl, R, gF, daf = inp.daf) => {
      const P = daf * (swl + inp.mHook) * 1e3 * G, Wc = inp.mCrane * 1e3 * G, Ho = inp.offlead / 100 * P, Hs = inp.sidelead / 100 * P, Hw = inp.Fwind * 1e3;
      const at = (z) => { const Mi = P * R + Wc * inp.eCrane + Ho * (inp.hTip + z) + Hw * (inp.hWind + z), Mo = Hs * (inp.hTip + z); return { M: gF * Math.hypot(Mi, Mo), Mi: gF * Mi, Mo: gF * Mo }; };
      return { N: gF * (P + Wc), V: gF * Math.hypot(Ho + Hw, Hs), T: gF * Hs * R, top: at(0), base: at(H), P };
    };
    const utilisation = (swl, R) => {
      const f = forces(swl, R, inp.gF), Nb = f.N + sec.A * H * st.rho * G * inp.gF;
      const tub = tubularMemberCheck({ code: 'ISO', section: sec, fy: st.fy, L: H, forces: { N: -Nb, Vy: f.V, Vz: 0, T: f.T, My: f.base.M, Mz: 0 }, options: { Ky: 2, Kz: 2, E: st.E } });
      const tau = f.T * (D / 2) / sec.J + 2 * f.V / sec.A;
      const sh = c202ShellCheck({ r: (D - t) / 2, t, l: H, E: st.E, nu: st.nu, fy: st.fy, sigmaA: -Nb / sec.A, sigmaM: -f.base.M / sec.Wely, tau });
      const n = Math.round(inp.nBolt), rb = inp.dBc / 2e3, bolts = Array.from({ length: n }, (_, i) => [rb * Math.cos(2 * Math.PI * i / n), rb * Math.sin(2 * Math.PI * i / n)]);
      const bg = boltGroupElastic({ bolts, Vx: f.V, Mz: f.T, N: -f.N, My: f.top.M });
      const br = boltResistance({ d: Number(inp.bolt.slice(1)) / 1e3, As: BOLT_STRESS_AREA[inp.bolt], grade: inp.boltGrade, FvEd: bg.maxShear, FtEd: Math.max(bg.maxTension, 0) });
      return { f, Nb, tub, sh, bg, br, tau, uc: Math.max(tub.governing.uc, sh.utilisation, br.ucCombined) };
    };
    const u = run(() => utilisation(inp.swl, inp.radius), 'Pedestal check'), f = u.f;
    const sigA = u.Nb / sec.A, sigB = f.base.M / sec.Wely, sigMax = sigA + sigB;
    // Serviceability: tilt and deflection of the pedestal top under unfactored working loads.
    const s = forces(inp.swl, inp.radius, 1), rot = (s.top.M * H / EI + s.V * H * H / (2 * EI)) / DEG, defl = s.top.M * H * H / (2 * EI) + s.V * H ** 3 / (3 * EI);
    const mTop = (inp.mCrane + inp.swl + inp.mHook) * 1e3 + 0.24 * sec.A * H * st.rho, f1 = Math.sqrt(3 * EI / (mTop * H ** 3)) / (2 * Math.PI);
    // Fatigue of the base girth weld: each lift cycles the hook-load moment from zero to its dynamic value at the cubic-mean load.
    const scf = scfGirthWeld({ D, t, deltaM: inp.misalign / 1e3, delta0: 0.1 * t }).scf;
    const Mfat = inp.kLoad * inp.daf * inp.swl * 1e3 * G * inp.radius * (1 + inp.offlead / 100 * (inp.hTip + H) / inp.radius), dS = Mfat / sec.Wely * scf;
    const sn = snConstants(run(() => snCurve(inp.snName, 'air'), 'S–N curve'), { thickness: t }), Nf = cyclesFromConstants(sn, dS), lifeYears = Math.min(Nf / inp.lifts / inp.dff, 1e4);
    // Allowable hook load at each radius: the load at which the governing pedestal utilisation reaches 1.0.
    const allow = (R) => { if (utilisation(0, R).uc >= 1) return 0; let lo = 0, hi = Math.max(inp.swl, 1); while (utilisation(hi, R).uc < 1 && hi < 1e5) hi *= 2; for (let i = 0; i < 30; i++) { const mid = 0.5 * (lo + hi); if (utilisation(mid, R).uc < 1) lo = mid; else hi = mid; } return lo; };
    const radii = linspace(Math.max(5, 0.25 * inp.radius), 1.5 * inp.radius, 16), curve = radii.map(allow), allowHere = allow(inp.radius);
    const checks = [
      check('Pedestal tube strength at the deck', u.tub.governing.uc, 1, 'interaction', `ISO 19902:2020 §13 ${u.tub.governing.name}, K = 2 (cantilever)`, 'Axial compression, bending, shear and torsion at the base of the pedestal against the tubular member resistance. Above 1.0 the pedestal wall yields or buckles locally on the compression side under the factored crane load.'),
      check('Shell buckling of the pedestal wall', u.sh.utilisation, 1, 'σj / fksd', 'DNV-RP-C202 §3.2 (unstiffened cylinder under combined axial, bending and shear stress)', 'Equivalent stress against the design shell buckling strength. Thin, large-diameter pedestals buckle before they yield; ring stiffeners or a thicker wall raise the strength.'),
      check('Slew-ring flange bolts', u.br.ucCombined, 1, 'interaction', `EN 1993-1-8 Table 3.4, ${Math.round(inp.nBolt)} × ${inp.bolt} grade ${inp.boltGrade}; elastic bolt-group distribution`, 'Tension in the most loaded bolt from the overturning moment (less the share of compression) with shear from horizontal load and slew torque. Bolt overload opens the bearing joint and is the classic crane-loss failure mode.'),
      check('Tilt at the slew bearing (working load)', rot, inp.rotLim, 'deg', 'Slew bearing manufacturer limit; EN 13852-1 §5 (supporting structure stiffness)', 'Rotation of the pedestal top under unfactored working loads. Excess tilt overloads the bearing raceway locally and makes the boom drift when slewing.'),
      fsCheck('Fatigue life of the base girth weld', Math.max(lifeYears, 1e-3), inp.life, `DNV-RP-C203 curve ${inp.snName} in air, SCF ${fmt(scf)} (§3.3.7), DFF ${inp.dff}`, 'Design fatigue life of the pedestal-to-deck weld under the lift cycles at the cubic-mean load. Pedestals see a full stress cycle on every lift, so the base weld is fatigue-sensitive.'),
    ];
    const gov = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = gov.uc <= 1, mPed = sec.A * H * st.rho / 1e3;
    // 3-D view: axial stress on the pedestal surface (compression negative) from N and the moment growing toward the base.
    const nT = 24, nZ = 8, nodes = [], vals = [], tris = [];
    for (let k = 0; k <= nZ; k++) for (let i = 0; i < nT; i++) { const th = 2 * Math.PI * i / nT, z = H * k / nZ, M = forces(inp.swl, inp.radius, inp.gF); const Mz = M.top.M + (M.base.M - M.top.M) * (1 - k / nZ); nodes.push([round(D / 2 * Math.cos(th), 5), round(D / 2 * Math.sin(th), 5), round(z, 5)]); vals.push(round((-u.Nb / sec.A - Mz / sec.Wely * Math.cos(th)) / 1e6, 4)); }
    for (let k = 0; k < nZ; k++) for (let i = 0; i < nT; i++) { const a = k * nT + i, b = k * nT + (i + 1) % nT, c = a + nT, d = b + nT; tris.push([a, b, d], [a, d, c]); }
    const tReq = inp.t * Math.max(u.tub.governing.uc, u.sh.utilisation) / 0.9;
    return {
      kpis: [
        { key: 'Mbase', label: 'Overturning moment at the deck', unit: 'kN·m', value: f.base.M / 1e3, explain: 'Factored resultant moment at the pedestal base from the hook load at its radius, crane weight eccentricity, offlead, sidelead and wind.' },
        { key: 'Nbase', label: 'Vertical force at the deck', unit: 'kN', value: u.Nb / 1e3, explain: 'Factored hook load with its dynamic factor plus crane and pedestal weight: the vertical reaction into the deck structure.' },
        { key: 'Vbase', label: 'Horizontal force', unit: 'kN', value: f.V / 1e3, explain: 'Factored shear from offlead, sidelead and wind.' },
        { key: 'Tbase', label: 'Slew torque', unit: 'kN·m', value: f.T / 1e3, explain: 'Factored torsion about the pedestal axis from the sidelead force at the working radius.' },
        { key: 'sigB', label: 'Bending stress at the base', unit: 'MPa', value: sigB / 1e6, explain: 'Extreme-fibre bending stress M/W in the pedestal wall at the deck.' },
        { key: 'sigMax', label: 'Largest compressive stress', unit: 'MPa', value: sigMax / 1e6, explain: 'Axial plus bending stress on the compression side at the base, to compare with the yield and shell buckling strengths.' },
        { key: 'ucTube', label: 'Pedestal utilisation', unit: '-', value: round(Math.max(u.tub.governing.uc, u.sh.utilisation), 5), limit: 1, sense: 'max', explain: 'Larger of the tubular member and shell buckling utilisations.' },
        { key: 'boltT', label: 'Largest bolt tension', unit: 'kN', value: Math.max(u.bg.maxTension, 0) / 1e3, explain: 'Tension in the most loaded slew-ring bolt from the elastic distribution of moment and axial force.' },
        { key: 'ucBolt', label: 'Bolt utilisation', unit: '-', value: round(u.br.ucCombined, 5), limit: 1, sense: 'max', explain: 'Combined shear and tension utilisation of the most loaded bolt.' },
        { key: 'rot', label: 'Tilt of the pedestal top', unit: 'deg', value: rot, limit: inp.rotLim, sense: 'max', explain: 'Rotation at the slew bearing under unfactored working loads.' },
        { key: 'defl', label: 'Deflection of the pedestal top', unit: 'mm', value: defl * 1e3, explain: 'Horizontal movement of the pedestal top under working loads.' },
        { key: 'f1', label: 'Natural frequency of crane on pedestal', unit: 'Hz', value: f1, explain: 'First bending frequency with the crane and hook load as a top mass. It should be well above the hoisting and slewing excitation and wave frequencies (typically > 1 Hz).' },
        { key: 'fatigueLife', label: 'Design fatigue life of the base weld', unit: 'yr', value: lifeYears, limit: inp.life, sense: 'min', explain: 'Calculated life divided by the design fatigue factor.' },
        { key: 'swlAllow', label: 'Allowable hook load at this radius', unit: 't', value: allowHere, limit: inp.swl, sense: 'min', explain: 'Hook load at which the pedestal, shell or bolts reach utilisation 1.0 at the working radius. The rated load must not exceed it.' },
        { key: 'capUse', label: 'Lifting capacity utilisation', unit: '-', value: allowHere > 0 ? inp.swl / allowHere : 99, limit: 1, sense: 'max', explain: 'Rated hook load over the load the pedestal allows at this radius.' },
      ],
      checks,
      plots: [
        { type: 'model3d', title: 'Axial stress in the pedestal wall', nodes, tris, lines: [], nodeValues: vals, valueLabel: 'Axial stress (compression negative)', valueUnit: 'MPa', explain: 'Longitudinal stress around and along the pedestal under the factored loads, hook side at +x. The largest compression is at the base on the hook side; openings and attachments should be kept away from that zone.' },
        { type: 'line', title: 'Allowable hook load against radius', x: { label: 'Working radius', unit: 'm' }, y: { label: 'Hook load', unit: 't' }, series: [{ name: 'Allowed by the pedestal', x: radii.map((r) => round(r)), y: curve.map((v) => round(v)) }, { name: 'Rated point', x: [inp.radius], y: [inp.swl], mode: 'markers' }],
          explain: 'Hook load at which the pedestal reaches its limit, for each radius, with the same dynamic factor and lead angles. The crane load chart must lie below this curve everywhere.' },
        { type: 'bar', title: 'Pedestal utilisation by check', categories: checks.map((c) => c.name), series: [{ name: 'Utilisation', y: checks.map((c) => round(c.uc)) }], y: { label: 'Utilisation', unit: '-' }, explain: 'All pedestal checks side by side; the tallest bar is what limits the crane on this pedestal.' },
      ],
      tables: [{ title: 'Pedestal design loads', columns: [{ key: 'q', label: 'Quantity' }, { key: 'top', label: 'At the slew ring' }, { key: 'base', label: 'At the deck' }, { key: 'u', label: 'Unit' }],
        rows: [{ q: 'Vertical force', top: round(f.N / 1e3), base: round(u.Nb / 1e3), u: 'kN' }, { q: 'Horizontal force', top: round(f.V / 1e3), base: round(f.V / 1e3), u: 'kN' }, { q: 'Moment in the boom plane', top: round(f.top.Mi / 1e3), base: round(f.base.Mi / 1e3), u: 'kN·m' }, { q: 'Moment across the boom plane', top: round(f.top.Mo / 1e3), base: round(f.base.Mo / 1e3), u: 'kN·m' }, { q: 'Resultant moment', top: round(f.top.M / 1e3), base: round(f.base.M / 1e3), u: 'kN·m' }, { q: 'Torque', top: round(f.T / 1e3), base: round(f.T / 1e3), u: 'kN·m' }],
        explain: 'Factored crane loads at the slew-ring flange and at the deck. The deck values are the support reactions to apply to the deck structure and its local reinforcement.' }],
      notes: [`${gov.name} governs at ${fmt(gov.uc)}. D/t = ${fmt(D / t)}; shell slenderness λs = ${fmt(u.sh.lambdaS)}.`, `Dynamic hook load ${fmt(f.P / 1e3)} kN at ${inp.radius} m; fatigue stress range ${fmt(dS / 1e6)} MPa per lift with SCF ${fmt(scf)}.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `The pedestal carries ${inp.swl} t at ${inp.radius} m with DAF ${inp.daf}: utilisation ${fmt(Math.max(u.tub.governing.uc, u.sh.utilisation))}, bolts ${fmt(u.br.ucCombined)}, tilt ${fmt(rot)}°, fatigue life ${fmt(lifeYears)} yr; it would allow ${fmt(allowHere)} t at this radius.`,
          `Carry the base reactions (M = ${fmt(f.base.M / 1e3)} kN·m, N = ${fmt(u.Nb / 1e3)} kN) into the deck with full-penetration welds and shear plates aligned with the deck girders, and check every load chart point against the allowable curve.`, 'EN 13852-1:2013 §5; API Spec 2C (2020) §5; ISO 19902:2020 §13; DNV-RP-C202 (2019) §3.2',
          `The pedestal shell weighs ${fmt(mPed)} t (${steelCo2Text(mPed)}). ${gov.uc < 0.6 ? `With ${fmt((1 - gov.uc) * 100, 2)} % reserve the wall could be thinned toward ${fmt(Math.max(tReq, inp.t * 0.6), 2)} mm if fatigue and tilt allow.` : 'The reserve is modest: keep the wall.'}`)
          : adv('engineering', 'fail', `${gov.name} is not met (utilisation ${fmt(gov.uc)}): the pedestal allows ${fmt(allowHere)} t at ${inp.radius} m against ${inp.swl} t rated.`,
            gov.name.startsWith('Slew-ring') ? `Use ${Math.ceil(inp.nBolt * u.br.ucCombined / 0.9)} bolts, a larger bolt or a larger bolt circle, and specify controlled preload.` : gov.name.startsWith('Tilt') ? `Increase the diameter to about ${fmt(inp.D * (rot / inp.rotLim) ** (1 / 3), 3)} mm (stiffness grows with D³) or shorten the pedestal.` : gov.name.startsWith('Fatigue') ? `Improve the base weld (grind flush to curve C), tighten the misalignment tolerance below ${inp.misalign} mm, or thicken the base can to about ${fmt(inp.t * (inp.life / Math.max(lifeYears, 1e-3)) ** (1 / 3), 3)} mm.` : `Thicken the wall to about ${fmt(Math.ceil(tReq), 3)} mm over the lower part, add ring stiffeners, or de-rate the crane to ${fmt(allowHere)} t at this radius.`,
            'EN 13852-1:2013 §5; ISO 19902:2020 §13; DNV-RP-C202 (2019) §3.2; DNV-RP-C203 (2021)', `A thicker base can of two diameters length at ${fmt(Math.ceil(tReq), 3)} mm adds about ${steelCo2Text(Math.max(tReq / inp.t - 1, 0.05) * mPed * Math.min(2 * D / H, 1))} De-rating the crane costs no steel but restricts operations.`),
        adv('economic', 'info', `Lifting capacity utilisation is ${fmt(allowHere > 0 ? inp.swl / allowHere : 99)}: the pedestal ${allowHere > inp.swl ? `has ${fmt(allowHere - inp.swl)} t of spare hook capacity at this radius, which protects a future crane upgrade` : 'limits the crane below its rating, which means more or smaller lifts and longer supply-vessel time alongside'}.`,
          'Size the pedestal for the heaviest crane the platform may carry in its life: pedestal steel is cheap at fabrication, while strengthening a pedestal offshore means a crane outage and hot work.', 'API Spec 2C (2020) design load philosophy; whole-life costing per ISO 15686-5'),
      ],
      publish: { 'topsides.craneMoment': f.base.M },
      quantities: { steel_t: mPed },
    };
  },
  verify: [{
    // Hook load only: M = DAF · m g · R = 1.3·50 000·9.80665·30 = 19 123 kN·m. Tube 2000 × 40 mm:
    // W = π(D⁴ − d⁴)/(32 D) = π(2⁴ − 1.92⁴)/64 = 0.118323 m³ → σb = 19 123 e3/0.118323 = 161.6 MPa.
    name: 'Cantilever pedestal under a hook load: M = DAF·P·R, σ = M/W', source: 'Statics and elementary beam theory (hand calculation)',
    inputs: { swl: 50, radius: 30, daf: 1.3, mHook: 0, mCrane: 0, eCrane: 0, offlead: 0, sidelead: 0, Fwind: 0, gF: 1, D: 2000, t: 40, dBc: 2150 },
    expect: [{ kpi: 'Mbase', value: 19123, tol: 1e-3 }, { kpi: 'sigB', value: 161.6, tol: 2e-3 }],
  }],
};
