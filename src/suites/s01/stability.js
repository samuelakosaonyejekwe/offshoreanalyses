// Structural stability: eigenvalue buckling of the frame under the factored storm load, system effective lengths,
// column and local buckling strength of every tubular member, and hydrostatic hoop buckling of unflooded members.
import { G, need, fmt, round, adv, check, fsCheck, linspace, steelCo2, CODE_CAVEAT, MORISON_CAVEAT, CONFIRM, BASIS_COST } from './_common.js';
import { stdInputs, buildJacket, stormCase, solve, memberForces, linearBuckling, plot3d, steelMass, groupSteel, foundNote, testFrame } from './_model.js';
import { isoTubularStrengths, hoopBucklingCoefficient } from '../../engine/structural/tubular.js';

export default {
  id: 's01-stability', suite: 1, title: 'Structural stability and member buckling',
  summary: 'Elastic critical load factors and buckling modes of the jacket, effective lengths implied by the frame, and column, local and hydrostatic buckling capacity of every member.',
  covers: [6], fidelity: 'numerical',
  standards: ['ISO 19902:2020 §13.2.3 (column and local buckling), §13.2.6 (hydrostatic pressure), Table 13.5-1', 'API RP 2A-WSD 22nd ed. §6.2.2, §6.2.5', 'EN 1993-1-1:2005 §5.2.1'],
  models: ['Euler column buckling theory', 'Eigenvalue (linear) buckling with consistent geometric stiffness', 'ISO 19902 column curve and local buckling strength', 'Elastic hoop buckling of ring-unstiffened tubes', 'Effective length from the elastic critical force'],
  limits: `Eigenvalue buckling of the perfect elastic frame: imperfections, plasticity and post-buckling are represented only through the code column curve (see the nonlinear and pushover modules). Joints are rigid and unstiffened; ring-stiffened or grouted members and shell buckling of large cans are not covered. ${MORISON_CAVEAT} ${CODE_CAVEAT}`,
  searchTerms: 'tubular member buckling jacket effective length elastic critical load',
  inputs: [
    ...stdInputs({ mesh: { nSub: 2, nPhase: 16 } }),
    { key: 'nModes', label: 'Buckling modes', type: 'int', default: 3, min: 1, max: 8, group: 'Solver', role: 'control' },
    { key: 'alphaMin', label: 'Required elastic critical load factor', unit: '-', type: 'number', default: 3, min: 1, max: 20, group: 'Design criteria', role: 'criteria', help: 'Project criterion on the eigenvalue of the factored load set. EN 1993-1-1 treats second-order effects as negligible above 10; member strength is covered by the column curve.' },
    { key: 'slMax', label: 'Slenderness limit KL/r for compression members', unit: '-', type: 'number', default: 120, min: 60, max: 250, group: 'Design criteria', role: 'criteria', help: `Recommended upper limit for primary bracing in compression. ${CONFIRM}` },
    { key: 'dtMax', label: 'D/t limit', unit: '-', type: 'number', default: 120, min: 40, max: 300, group: 'Design criteria', role: 'criteria', help: `Validity limit of the ISO 19902 tubular strength formulas (D/t ≤ 120). ${CONFIRM}` },
    { key: 'gRh', label: 'Resistance factor — hoop buckling γR,h', unit: '-', type: 'number', default: 1.25, min: 1, max: 2, group: 'Design criteria', role: 'criteria', help: `ISO 19902 only. ${CONFIRM}` },
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'lambdaCr' },
  compute(inp, ctx) {
    // SI conversion inside the shared model; forces reported in kN and MN.
    const J = buildJacket(inp, ctx), st = stormCase(J, inp, { name: 'Storm' }), res = solve(J, [st.lc]).Storm;
    const bk = linearBuckling(J, 'Storm', inp.nModes), l1 = bk.factors[0], F = memberForces(J, res), E = J.steel.E, fy = J.steel.fy, rho = inp.rhoW ?? 1025;
    const rows = J.orig.filter((m) => m.tube).map((m) => {
      const s = isoTubularStrengths({ D: m.sec.D, t: m.sec.t, fy, E, L: m.L, K: m.K }), N = F[m.i].N, comp = Math.max(-N, 0);
      const PE = Math.PI ** 2 * E * m.sec.Iy / (m.K * m.L) ** 2, Kimp = comp > 0 ? Math.PI * Math.sqrt(E * m.sec.Iy / (l1 * comp)) / m.L : null;
      // hydrostatic head on unflooded members below the water line (design head with the gravity-load factor)
      const head = !m.flooded && m.zmid < 0 ? -m.zmid + (inp.surge || 0) : 0, p = rho * G * head * inp.gG, sh = p * m.sec.D / (2 * m.sec.t);
      const fhe = 2 * hoopBucklingCoefficient({ D: m.sec.D, t: m.sec.t, L: m.L }).Ch * E * m.sec.t / m.sec.D;
      return { i: m.i, m, N, comp, Pc: s.Pc, PE, lam: s.lambda, fc: s.fc, fyc: s.fyc, slender: m.K * m.L / s.r, dt: m.sec.D / m.sec.t, ucCol: comp / (s.Pc / inp.gRc), Kimp, sh, fh: s.fh ?? Math.min(fhe, fy), ucHoop: sh / ((s.fh ?? Math.min(fhe, fy)) / inp.gRh), fxe: s.fxe };
    });
    need(rows.length, 'The model has no tubular members to check.');
    const comp = rows.filter((r) => r.comp > 0), col = rows.reduce((a, r) => (r.ucCol > a.ucCol ? r : a)), sl = comp.reduce((a, r) => (r.slender > a.slender ? r : a), comp[0] ?? rows[0]);
    const dt = rows.reduce((a, r) => (r.dt > a.dt ? r : a)), hoop = rows.reduce((a, r) => (r.ucHoop > a.ucHoop ? r : a));
    const groups = {}; for (const r of rows) { const g = groups[r.m.group] || (groups[r.m.group] = { label: r.m.label, group: r.m.group, r }); if (r.ucCol > g.r.ucCol) g.r = r; }
    const gl = Object.values(groups), mode = bk.modes[0], amp = Array.from({ length: J.model.nodes.length }, (_, i) => Math.hypot(mode[6 * i], mode[6 * i + 1], mode[6 * i + 2]));
    // member that moves most in the first mode
    const part = J.orig.map((m, i) => Math.max(...J.memberMap[i].flatMap((e) => { const q = J.model.members[e]; return [amp[q.a], amp[q.b]]; }))), iMode = part.indexOf(Math.max(...part));
    const checks = [
      fsCheck('Elastic critical load factor of the frame', l1, inp.alphaMin, 'load factor', 'Eigenvalue buckling; EN 1993-1-1:2005 §5.2.1; project criterion', 'Multiple of the factored load set at which the perfect elastic frame buckles. A low value means a slender frame that is sensitive to imperfections and second-order effects.'),
      ...gl.map((g) => check(`Column buckling — ${g.label}`, g.r.comp / 1e6, g.r.Pc / inp.gRc / 1e6, 'MN', `ISO 19902:2020 §13.2.3, Eq. 13.2-5 with K = ${g.r.m.K}, γR,c = ${inp.gRc}`, 'Factored axial compression against the column buckling resistance of the most loaded member of the group. If exceeded the member bows sideways and sheds its load to the rest of the frame.')),
      check('Slenderness KL/r of compression members', sl.slender, inp.slMax, '-', 'API RP 2A-WSD 22nd ed. §6.2.2 commentary; project criterion', 'Very slender compression members lose strength quickly, vibrate in waves and are easily damaged in handling.'),
      check('Diameter to thickness ratio D/t', dt.dt, inp.dtMax, '-', 'ISO 19902:2020 §13.1 (validity of the tubular formulas)', 'Above this ratio the wall buckles locally before the member develops its strength and the code formulas no longer apply.'),
      check('Hydrostatic hoop buckling of unflooded members', hoop.sh / 1e6, hoop.fh / inp.gRh / 1e6, 'MPa', `ISO 19902:2020 §13.2.6, Eq. 13.2-23 to 13.2-29, γR,h = ${inp.gRh}`, 'Hoop stress from external water pressure against the hoop buckling strength of the unstiffened tube. If exceeded the member collapses inward unless it is flooded or ring-stiffened.'),
    ];
    const steel = steelMass(J, inp), ok = checks.every((c) => c.uc <= 1), bad = checks.filter((c) => c.uc > 1), gSteel = groupSteel(J, col.m.group);
    const lamAx = linspace(0.05, 3, 40), curve = lamAx.map((l) => (l <= 1.34 ? 1 - 0.278 * l * l : 0.9 / (l * l)));
    return {
      kpis: [
        { key: 'lambdaCr', label: 'Elastic critical load factor αcr', unit: '-', value: l1, limit: inp.alphaMin, sense: 'min', explain: 'Factor on the factored gravity and storm loads at which the first buckling mode of the frame occurs. It measures the global stability reserve of the elastic structure.' },
        { key: 'lambda2', label: 'Second buckling factor', unit: '-', value: bk.factors[1] ?? l1, explain: 'A second factor close to the first means several braces are near-equally critical, as in a regular jacket face.' },
        { key: 'maxColUc', label: 'Highest column buckling utilisation', unit: '-', value: col.ucCol, limit: 1, sense: 'max', explain: 'Factored compression over column buckling resistance for the governing member (axial only; the in-place module adds bending interaction).' },
        { key: 'maxSlender', label: 'Largest slenderness KL/r in compression', unit: '-', value: sl.slender, limit: inp.slMax, sense: 'max', explain: 'Slenderness of the most slender member that carries compression in this load case.' },
        { key: 'lamBar', label: 'Reduced slenderness λ of the governing column', unit: '-', value: col.lam, explain: 'Column slenderness parameter √(fyc/fe). Below about 0.5 the member is governed by yield, above 1.34 by elastic buckling.' },
        { key: 'Kimplied', label: 'Effective length factor implied by the frame (mode member)', unit: '-', value: rows.find((r) => r.i === iMode)?.Kimp ?? 0, explain: 'K back-calculated from the eigenvalue for the member that dominates the first mode: π√(EI/(αcr·N))/L. Compare with the K assumed in the code check; a larger implied value means the assumed K is unconservative.' },
        { key: 'maxDt', label: 'Largest D/t', unit: '-', value: dt.dt, limit: inp.dtMax, sense: 'max', explain: 'Most thin-walled member. Local buckling reduces strength once fy/fxe exceeds 0.170, that is D/t above about 60 for 355 MPa steel.' },
        { key: 'maxHoopUc', label: 'Highest hydrostatic hoop utilisation', unit: '-', value: hoop.ucHoop, limit: 1, sense: 'max', explain: 'External pressure hoop stress over hoop buckling strength for unflooded members; zero when all submerged members are flooded.' },
      ],
      checks,
      plots: [
        plot3d(J, { title: `First buckling mode (αcr = ${fmt(l1)})`, values: amp, label: 'Mode amplitude', unit: '-', shape: mode, explain: `Shape of the first elastic buckling mode, normalised to 1. The members that bow are the stability-critical ones (here member ${iMode}, ${J.orig[iMode].label.toLowerCase()}); a sway of the whole frame instead points to the foundation or leg stiffness.` }),
        { type: 'scatter', title: 'Compression members on the column curve', x: { label: 'Reduced slenderness λ', unit: '-' }, y: { label: 'Axial stress / local buckling strength', unit: '-' },
          series: [{ name: 'ISO 19902 column curve fc/fyc', x: lamAx.map((v) => round(v)), y: curve.map((v) => round(v)) }, { name: 'Members (σc/fyc)', x: comp.map((r) => round(r.lam)), y: comp.map((r) => round(r.comp / r.m.sec.A / r.fyc)) }],
          explain: 'Each point is a compression member: its slenderness against its factored compressive stress ratio. Points must lie below the curve (with the resistance factor on top); points far to the right are slender members that would benefit from an extra node or a larger diameter.' },
        { type: 'bar', title: 'Buckling load factors', categories: bk.factors.map((f, i) => `Mode ${i + 1}`), series: [{ name: 'Load factor', y: bk.factors.map((f) => round(f)) }], y: { label: 'Load factor', unit: '-' }, explain: 'Elastic critical load factors of the lowest modes. Closely spaced values indicate similar members buckling in turn.' },
      ],
      tables: [
        { title: 'Member buckling capacities (most utilised)', columns: [{ key: 'i', label: 'Member' }, { key: 'g', label: 'Group' }, { key: 'L', label: 'Length', unit: 'm' }, { key: 'K', label: 'K used' }, { key: 'Kimp', label: 'K implied' }, { key: 'sl', label: 'KL/r' }, { key: 'lam', label: 'λ' }, { key: 'N', label: 'Axial (tension +)', unit: 'kN' }, { key: 'PE', label: 'Euler load', unit: 'kN' }, { key: 'Pc', label: 'Column strength', unit: 'kN' }, { key: 'uc', label: 'Utilisation' }, { key: 'dt', label: 'D/t' }],
          rows: rows.slice().sort((a, b) => b.ucCol - a.ucCol).slice(0, 15).map((r) => ({ i: r.i, g: r.m.label, L: round(r.m.L), K: r.m.K, Kimp: r.Kimp == null ? '–' : round(Math.min(r.Kimp, 99), 3), sl: round(r.slender, 3), lam: round(r.lam, 3), N: round(r.N / 1e3), PE: round(r.PE / 1e3), Pc: round(r.Pc / 1e3), uc: round(r.ucCol), dt: round(r.dt, 3) })),
          explain: 'Euler load, representative column strength and utilisation of the fifteen most loaded compression members. "K implied" is meaningful for the members that take part in the first mode; lightly loaded members show large values.' },
        { title: 'Buckling modes', columns: [{ key: 'k', label: 'Mode' }, { key: 'f', label: 'Load factor' }, { key: 'm', label: 'Member with largest amplitude' }], rows: bk.factors.map((f, k) => { const md = bk.modes[k], a = J.orig.map((m, i) => Math.max(...J.memberMap[i].flatMap((e) => { const q = J.model.members[e]; return [Math.hypot(md[6 * q.a], md[6 * q.a + 1], md[6 * q.a + 2]), Math.hypot(md[6 * q.b], md[6 * q.b + 1], md[6 * q.b + 2])]; }))), i = a.indexOf(Math.max(...a)); return { k: k + 1, f: round(f), m: `${i} (${J.orig[i].label})` }; }), explain: 'Load factor and the member that deforms most in each mode: the instability modes in order of criticality.' },
      ],
      notes: [`First mode at αcr = ${fmt(l1)} is dominated by member ${iMode} (${J.orig[iMode].label.toLowerCase()}); governing column check: member ${col.i} (${col.m.label.toLowerCase()}), compression ${fmt(col.comp / 1e6)} MN against ${fmt(col.Pc / inp.gRc / 1e6)} MN.`,
        `Local buckling: fy/fxe up to ${fmt(Math.max(...rows.map((r) => fy / r.fxe)))} (reduction starts at 0.170). ${hoop.ucHoop > 0 ? `Deepest unflooded member sees ${fmt(hoop.sh / 1e6)} MPa hoop stress.` : 'No unflooded member is submerged.'}`,
        `Wave by ${st.wf.label}, heading ${fmt(st.beta)}°; the buckling solution converged in ${bk.iterations} subspace iterations.`, foundNote(J)],
      advice: [
        ok ? adv('engineering', 'ok', `Stability is adequate: αcr = ${fmt(l1)} (≥ ${fmt(inp.alphaMin)}), column utilisation ${fmt(col.ucCol)}, KL/r ≤ ${fmt(sl.slender)}, D/t ≤ ${fmt(dt.dt)}.`, col.ucCol < 0.5 ? `The governing ${col.m.label.toLowerCase()} use only ${fmt(col.ucCol * 100)} % of their column strength: a smaller diameter with the same slenderness limit, or one bay less, is worth testing against fatigue and installation cases.` : 'Keep the assumed effective length factors under review against the implied values in the table, and confirm that braces below the water line are flooded or checked for hydrostatic collapse.', 'ISO 19902:2020 §13.2.3, §13.2.6 and Table 13.5-1', `Buckling strength grows with diameter, not thickness: for slender braces a larger, thinner tube gives the same capacity with less steel, but more wave load. The ${col.m.label.toLowerCase()} hold ${steelCo2(gSteel)}`)
          : adv('engineering', 'fail', `Stability criteria are not met: ${bad.map((c) => `${c.name} (${fmt(c.uc)})`).join('; ')}.`, col.ucCol > 1 ? `Increase the diameter of the ${col.m.label.toLowerCase()} (radius of gyration needs about ${fmt((Math.sqrt(col.ucCol / 0.9) - 1) * 100)} % more in the elastic range) or halve the buckling length with an extra node.` : hoop.ucHoop > 1 ? 'Flood the affected members, add ring stiffeners or increase their wall thickness.' : 'Shorten or enlarge the slender members, or stiffen the foundation if the first mode is a sway of the frame.', 'ISO 19902:2020 §13.2.3, §13.2.6', `Enlarging the ${col.m.label.toLowerCase()} is of the order of ${steelCo2(gSteel * 0.15)}`),
        adv('economic', 'info', `Column strength is ${fmt(col.fc / col.fyc * 100)} % of the local buckling strength for the governing member (λ = ${fmt(col.lam)}); the ${col.m.label.toLowerCase()} weigh ${fmt(gSteel)} t.`, 'Slender braces pay for steel that buckling does not let them use. Compare the cost of an additional framing level (more joints and welds, shorter members) with heavier braces: joints dominate fabrication hours, tonnes dominate material and lift cost.', BASIS_COST),
      ],
      quantities: { steel_t: round(steel) },
    };
  },
  verify: [
    // Fixed-free column under tip load: Pcr = π²EI/(4L²) = π²·2.1e11·0.054485/(4·95²) = 3.1282 MN. With W = 1 MN, αcr = 3.128;
    // effective length factor K = 2.0 by definition.
    { name: 'Euler cantilever column: Pcr = π²EI/(2L)²', source: 'Euler (1744); Timoshenko & Gere, Theory of Elastic Stability §2.1',
      inputs: { geom: testFrame('column'), H: 0, Uc: 0, Uw: 0, windA: 0, swCont: 0, Wtop: 1, found: 'fixed', gG: 1, gE: 1, nSub: 8, nModes: 1 },
      expect: [{ kpi: 'lambdaCr', value: 3.1282, tol: 0.002 }, { kpi: 'Kimplied', value: 2.0, tol: 0.002 }] },
  ],
};
