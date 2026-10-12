// Fatigue of the tubular joints: deterministic wave blocks are stepped past the structure to give the nominal brace
// stress ranges, Efthymiou SCFs turn them into hot-spot ranges at eight points on both sides of each weld, and
// Palmgren–Miner with a DNV-RP-C203 S–N curve gives the life. A spectral (Dirlik) estimate cross-checks the governing joint.
import { need, run, fmt, round, adv, check, fsCheck, steelCo2, CODE_CAVEAT, MORISON_CAVEAT, CONFIRM, BASIS_COST } from './_common.js';
import { geomInputs, envInputs, foundationInputs, meshInputs, buildJacket, waveField, envLoader, solve, jointData, jointLoads, plot3d, steelMass, foundNote, testFrame } from './_model.js';
import { scfParameters, scfTY, hotSpotStresses } from '../../engine/structural/scf.js';
import { snCurve, snConstants, cyclesFromConstants, SN_CURVE_NAMES } from '../../engine/fatigue/sn.js';
import { longTermSpectralDamage } from '../../engine/fatigue/spectral.js';
import { fatigueLife } from '../../engine/fatigue/damage.js';
import { jonswap } from '../../engine/metocean/spectra.js';

const YEAR = 31557600;
const env = envInputs({ theory: 'wheeler' });
const pick = (k) => env.find((i) => i.key === k);

export default {
  id: 's01-fatigue', suite: 1, title: 'Fatigue and structural life of tubular joints',
  summary: 'Hot-spot stress ranges, Miner damage and fatigue life of every brace-to-chord joint from a long-term wave height distribution, with the critical locations, the damage by wave class and a spectral cross-check.',
  covers: [8], fidelity: 'numerical',
  standards: ['DNV-RP-C203 (2021) §2.4 (S–N curves), §3.3 and Appendix B (SCFs), §5 (simplified and spectral methods)', 'ISO 19902:2020 §16 (fatigue)', 'API RP 2A-WSD 22nd ed. §8', 'NORSOK N-004 (2022) §8 (design fatigue factors)'],
  models: ['Palmgren–Miner cumulative damage', 'S–N curves with thickness correction (DNV-RP-C203)', 'Efthymiou stress concentration factors', 'Hot-spot stress method (eight points, chord and brace side)', 'Deterministic wave-block fatigue with Morison loading', 'Spectral fatigue with stress transfer function and Dirlik rainflow-range density', 'JONSWAP spectrum'],
  limits: `Braces are treated as T/Y joints with parametric SCFs (no influence functions for K/KT load sharing, no weld-profile or grinding benefit, no joint flexibility), one wave heading carries the stated share of all waves, and dynamic response enters only through the amplification factor. The stress transfer function of the spectral check is linearised from the wave blocks. Crack growth and inspection updating are in the reassessment module. ${CODE_CAVEAT} ${MORISON_CAVEAT}`,
  publishes: ['life.fatigueYears'],
  searchTerms: 'tubular joint fatigue jacket hot spot stress Miner S-N curve',
  inputs: [
    ...geomInputs(),
    { key: 'blocks', label: 'Long-term wave blocks (one year)', type: 'table', group: 'Fatigue environment', role: 'bc',
      columns: [{ key: 'H', label: 'Wave height', unit: 'm', type: 'number' }, { key: 'T', label: 'Period', unit: 's', type: 'number' }, { key: 'n', label: 'Waves per year', unit: '-', type: 'number' }],
      default: [{ H: 0.75, T: 4.5, n: 3.2e6 }, { H: 2, T: 6, n: 1.4e6 }, { H: 3.5, T: 7.5, n: 3.0e5 }, { H: 5.5, T: 9, n: 5.0e4 }, { H: 8, T: 10.5, n: 3000 }, { H: 12, T: 12, n: 60 }],
      help: 'Individual wave height classes with their associated period and yearly number of waves (about 5 million waves per year in total), from the wave height exceedance curve of the site. More, narrower classes refine the result at the cost of run time.' },
    { key: 'scatter', label: 'Sea-state scatter diagram (spectral check)', type: 'table', group: 'Fatigue environment', role: 'bc',
      columns: [{ key: 'Hs', label: 'Hs', unit: 'm', type: 'number' }, { key: 'Tp', label: 'Tp', unit: 's', type: 'number' }, { key: 'p', label: 'Occurrence', unit: '%', type: 'number' }],
      default: [{ Hs: 0.5, Tp: 5, p: 25 }, { Hs: 1.5, Tp: 7, p: 38 }, { Hs: 2.5, Tp: 8, p: 22 }, { Hs: 3.5, Tp: 9, p: 10 }, { Hs: 5, Tp: 10.5, p: 4 }, { Hs: 7, Tp: 12, p: 0.9 }, { Hs: 9, Tp: 13.5, p: 0.1 }],
      help: 'Used for the spectral cross-check of the governing joint only. Occurrences are normalised to 100 %.' },
    { ...pick('heading'), group: 'Fatigue environment' },
    { key: 'dirFrac', label: 'Share of waves from this heading', unit: '-', type: 'number', default: 0.4, min: 0.05, max: 1, group: 'Fatigue environment', role: 'bc', help: '1.0 puts all waves on the analysed heading (conservative). 0.3–0.5 is typical for the dominant sector; analyse the other headings with their shares and add the damages.' },
    { ...pick('theory'), group: 'Fatigue environment' },
    { ...pick('Cd'), default: 0.7, group: 'Hydrodynamic coefficients' }, { ...pick('Cm'), default: 2.0, help: 'Fatigue waves are small and inertia dominated: about 2.0 with Cd 0.65–0.8 is usual for fatigue analysis.' },
    pick('mg'), pick('mgDepth'), { ...pick('kinFac'), default: 1.0 },
    { ...pick('appF'), default: 1.0, help: 'Factor on the Morison load of the modelled members. Appurtenances raise the global load but not the local wave load that bends a brace, so 1.0 is the usual fatigue basis; use a higher value only where conductor loads dominate the joint.' }, pick('rhoW'),
    { key: 'daf', label: 'Dynamic amplification on stress ranges', unit: '-', type: 'number', default: 1.1, min: 1, max: 3, group: 'Hydrodynamic coefficients', role: 'criteria', help: 'From the dynamic module at the short fatigue wave periods (1.0 for a stiff, quasi-static platform).' },
    ...foundationInputs(),
    { key: 'canT', label: 'Leg joint-can wall thickness', unit: 'mm', type: 'number', default: 55, min: 0, max: 200, group: 'Joint details', role: 'geometry', design: true, help: 'Chord wall at leg joints; 0 uses the leg wall.' },
    { key: 'chordC', label: 'Chord end fixity parameter C', unit: '-', type: 'number', default: 0.7, min: 0.5, max: 1, group: 'Joint details', role: 'bc' },
    { key: 'curve', label: 'S–N curve', type: 'select', default: 'T', group: 'Joint details', role: 'criteria', options: SN_CURVE_NAMES.map((n) => ({ value: n, label: n === 'T' ? 'T (tubular joints)' : n })), help: 'DNV-RP-C203 curve of the weld class: T for tubular joints with hot-spot stress, D for girth welds.' },
    { key: 'envSN', label: 'Corrosion environment', type: 'select', default: 'seawater-cp', group: 'Joint details', role: 'environment', options: [{ value: 'air', label: 'Air (above the splash zone, coated)' }, { value: 'seawater-cp', label: 'Seawater with cathodic protection' }, { value: 'free-corrosion', label: 'Free corrosion (splash zone, failed CP)' }] },
    { key: 'life', label: 'Design service life', unit: 'years', type: 'number', default: 25, min: 1, max: 100, group: 'Design criteria', role: 'criteria', bind: 'life.designYears' },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 3, min: 1, max: 10, group: 'Design criteria', role: 'criteria', help: `1 to 10 by consequence of failure and access for inspection (NORSOK N-004 Table 8-1; ISO 19902 Table 16.12-1): 10 for uninspectable, substantial-consequence joints below water. ${CONFIRM}` },
    { key: 'damage0', label: 'Damage already accumulated', unit: '-', type: 'number', default: 0, min: 0, max: 5, group: 'Design criteria', role: 'ic', help: 'Miner sum consumed before the start of the period (0 for a new structure).' },
    ...meshInputs({ nSub: 1, nPhase: 12 }),
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'maxRange' },
  meshTime: { key: 'nPhase', label: 'Wave phase steps', levels: [6, 12, 24], qoi: 'maxRange' },
  compute(inp, ctx) {
    // SI conversion: mm → m; stresses in Pa internally and MPa in the results. No current, wind or gravity: only ranges matter.
    const blocks = (inp.blocks || []).map((b) => ({ H: Number(b.H), T: Number(b.T), n: Number(b.n) * inp.dirFrac })).filter((b) => b.H > 0 && b.n > 0);
    need(blocks.length, 'Give at least one wave block with a positive height and number of waves.');
    need(blocks.every((b) => b.T >= 2 && b.T <= 30 && Number.isFinite(b.n)), 'Wave block periods must lie between 2 and 30 s.');
    const J = buildJacket(inp, ctx), jd = jointData(J, { canT: inp.canT }), fi = { ...inp, Uc: 0, surge: 0, block: 1 };
    need(jd.length, 'The model has no brace-to-chord joint to assess for fatigue.');
    const cases = [], nP = inp.nPhase;
    blocks.forEach((b, ib) => {
      const wf = run(() => waveField(fi, { H: b.H, T: b.T }), `Wave block ${ib + 1} (H ${b.H} m, T ${b.T} s)`), ld = envLoader(J, fi, wf, inp.heading);
      b.theory = wf.label;
      for (let i = 0; i < nP; i++) cases.push({ name: `${ib}:${i}`, nodal: [], member: ld.at(i * b.T / nP, true).items });
    });
    ctx.progress?.(0.3, 'Solving wave phases');
    const res = solve(J, cases), curve = snCurve(inp.curve, inp.envSN);
    // SCFs once per joint. Hot-spot stresses are linear in the three nominal stress components, so the influence
    // coefficients of the 16 points (chord side 0–7, brace side 8–15) are stored and the phase loop is plain arithmetic.
    const sc = jd.map((j) => {
      const p = scfParameters({ D: j.D, T: j.T, d: j.d, t: j.t, L: j.Lc, thetaDeg: Math.min(j.thetaDeg, 90) }), s = scfTY({ ...p, chordEnds: 'general', C: inp.chordC }), c = new Float64Array(48);
      const sides = [{ scfAC: s.axial.chordCrown, scfAS: s.axial.chordSaddle, scfMIP: s.ipb.chordCrown, scfMOP: s.opb.chordSaddle }, { scfAC: s.axial.braceCrown, scfAS: s.axial.braceSaddle, scfMIP: s.ipb.braceCrown, scfMOP: s.opb.braceSaddle }];
      sides.forEach((sd, q) => ['sigmaX', 'sigmaMy', 'sigmaMz'].forEach((key, m) => { const h = hotSpotStresses({ [key]: 1, ...sd }).sigma; for (let pt = 0; pt < 8; pt++) c[3 * (8 * q + pt) + m] = h[pt]; }));
      return { p, s, c, warn: s.warnings, sn: [snConstants(curve, { thickness: j.T }), snConstants(curve, { thickness: j.t })] };
    });
    const nJ = jd.length, dmg = jd.map(() => [new Float64Array(8), new Float64Array(8)]), rng = blocks.map(() => jd.map(() => [new Float64Array(8), new Float64Array(8)]));
    const hist = new Float64Array(nJ * 16 * nP);
    blocks.forEach((b, ib) => {
      for (let i = 0; i < nP; i++) {
        const L = jointLoads(J, res[`${ib}:${i}`], jd);
        for (let k = 0; k < nJ; k++) { const c = sc[k].c, sx = L[k].sx, si = L[k].sip, so = L[k].sop, o = (k * 16) * nP + i; for (let q = 0; q < 16; q++) hist[o + q * nP] = c[3 * q] * sx + c[3 * q + 1] * si + c[3 * q + 2] * so; }
      }
      // range of each point over the cycle; the sampled extremes are refined with a parabola through three phases
      const ext = (o, sgn) => { let im = 0; for (let i = 1; i < nP; i++) if (sgn * hist[o + i] > sgn * hist[o + im]) im = i; const y0 = hist[o + (im + nP - 1) % nP], y1 = hist[o + im], y2 = hist[o + (im + 1) % nP], den = y0 - 2 * y1 + y2; return nP >= 4 && sgn * den < 0 ? y1 - (y2 - y0) ** 2 / (8 * den) : y1; };
      for (let k = 0; k < nJ; k++) for (let q = 0; q < 16; q++) { const o = (k * 16 + q) * nP, r = (ext(o, 1) - ext(o, -1)) * inp.daf, sd = q >> 3, pt = q & 7; rng[ib][k][sd][pt] = r; if (r > 0) dmg[k][sd][pt] += b.n / cyclesFromConstants(sc[k].sn[sd], r); }
    });
    const rows = jd.map((j, k) => { let d = 0, q0 = 0, p0 = 0; for (let q = 0; q < 2; q++) for (let p = 0; p < 8; p++) if (dmg[k][q][p] > d) { d = dmg[k][q][p]; q0 = q; p0 = p; } return { k, j, d, q: q0, p: p0, life: d > 0 ? 1 / d : Infinity, rMax: Math.max(...rng.map((r) => r[k][q0][p0])) }; });
    const gov = rows.reduce((a, r) => (r.d > a.d ? r : a));
    need(gov.d > 0, 'No fatigue damage was computed: the wave blocks produce no stress range at the joints (below the fatigue threshold or no wave load).');
    const fl = fatigueLife({ damageRate: gov.d, designLife: inp.life, dff: inp.dff, accumulatedDamage: inp.damage0 }), lifeReq = inp.life * inp.dff, cap = (v) => Math.min(v, 1e6);
    const gb = blocks.map((b, ib) => { const r = rng[ib][gov.k][gov.q][gov.p]; return { ...b, r, N: r > 0 ? cyclesFromConstants(sc[gov.k].sn[gov.q], r) : Infinity }; }), gd = gb.map((b) => (Number.isFinite(b.N) ? b.n / b.N : 0));
    // spectral cross-check of the governing hot spot: transfer function = range per unit wave height at each block period
    let lifeSp = 1e6, spNote = 'Spectral check not run (no scatter diagram).';
    const sea = (inp.scatter || []).map((s) => ({ Hs: Number(s.Hs), Tp: Number(s.Tp), p: Number(s.p) })).filter((s) => s.Hs > 0 && s.Tp > 0 && s.p > 0);
    if (sea.length) {
      const tf = gb.map((b) => ({ f: 1 / b.T, h: b.r / b.H })).sort((a, b) => a.f - b.f), H = (f) => { if (f <= tf[0].f) return tf[0].h * (f / tf[0].f) ** 2; if (f >= tf[tf.length - 1].f) return tf[tf.length - 1].h; let i = 0; while (tf[i + 1].f < f) i++; return tf[i].h + (tf[i + 1].h - tf[i].h) * (f - tf[i].f) / (tf[i + 1].f - tf[i].f); };
      const ps = sea.reduce((s, q) => s + q.p, 0), f = Array.from({ length: 48 }, (_, i) => 0.03 + i * 0.01);
      const lt = run(() => longTermSpectralDamage({ states: sea.map((s) => ({ probability: s.p / ps, f, waveSpectrum: f.map((x) => 2 * Math.PI * jonswap(2 * Math.PI * x, { Hs: s.Hs, Tp: s.Tp })), transfer: f.map(H), label: `Hs ${s.Hs}` })), sn: sc[gov.k].sn[gov.q], duration: YEAR, method: 'dirlik' }), 'Spectral fatigue');
      lifeSp = lt.damage > 0 ? cap(1 / (lt.damage * inp.dirFrac)) : 1e6;
      spNote = `Spectral (Dirlik) life of the governing hot spot ${fmt(lifeSp)} years against ${fmt(cap(gov.life))} years from the wave blocks; the two long-term descriptions should agree within a factor of about 2–3.`;
    }
    const nShort = rows.filter((r) => r.life < lifeReq).length, top = rows.slice().sort((a, b) => b.d - a.d), nodeU = new Array(J.model.nodes.length).fill(0); rows.forEach((r) => { nodeU[r.j.node] = Math.max(nodeU[r.j.node], Math.min(r.d * lifeReq, 99)); });
    const steel = steelMass(J, inp), ok = fl.utilisation <= 1, sn = sc[gov.k].sn[gov.q], side = gov.q === 0 ? 'chord' : 'brace', m = gb[gb.length - 1].r >= sn.sKnee ? sn.m1 : (sn.m2 ?? sn.m1);
    // stress reduction needed: damage ∝ range^m on the governing slope; SCF ∝ about 1/T^1.5 for the chord can
    const sRed = fl.utilisation > 0 ? fl.utilisation ** (-1 / m) : 1, tReq = gov.j.T * (1 / sRed) ** (1 / 1.5), canPer = Math.PI * gov.j.D * 2.5 * gov.j.D * J.steel.rho / 1e3;
    const sRange = Array.from({ length: 30 }, (_, i) => 10 ** (0.6 + i * 0.06) * 1e6);
    return {
      kpis: [
        { key: 'fatigueLife', label: 'Calculated fatigue life (governing joint)', unit: 'years', value: cap(gov.life), limit: lifeReq, sense: 'min', explain: `Years to a Miner sum of 1.0 at the worst hot spot, without the design fatigue factor. It must exceed the service life times the factor (${fmt(inp.life)} × ${fmt(inp.dff)} = ${fmt(lifeReq)} years).` },
        { key: 'designLife', label: 'Design fatigue life (with design fatigue factor)', unit: 'years', value: cap(fl.designLife), explain: 'Calculated life divided by the design fatigue factor: the service life the joint is good for under the code.' },
        { key: 'damage', label: 'Factored damage over the service life', unit: '-', value: Math.min(fl.utilisation, 1e6), limit: 1, sense: 'max', explain: 'Design fatigue factor times the Miner sum accumulated over the service life, including damage already consumed. Above 1.0 the joint does not meet its fatigue requirement.' },
        { key: 'maxRange', label: 'Hot-spot stress range in the largest wave block', unit: 'MPa', value: gb[gb.length - 1].r / 1e6, explain: 'Stress range at the governing hot spot when the highest wave of the table passes, including the dynamic amplification.' },
        { key: 'damageYear', label: 'Annual damage of the governing joint', unit: '1/yr', value: gov.d, explain: 'Miner sum added per year at the worst hot spot; its inverse is the calculated life.' },
        { key: 'lifeSpectral', label: 'Spectral (Dirlik) life of the governing hot spot', unit: 'years', value: lifeSp, explain: 'Independent estimate from the scatter diagram, a JONSWAP spectrum and the stress transfer function. A large difference from the wave-block life points to an inconsistent long-term wave description.' },
        { key: 'nShort', label: 'Joints below the required life', unit: '-', value: nShort, limit: 0, sense: 'max', explain: `Number of joints whose calculated life is shorter than ${fmt(lifeReq)} years.` },
        { key: 'scfGov', label: 'Governing SCF at the critical hot spot', unit: '-', value: Math.max(...(gov.q === 0 ? [sc[gov.k].s.axial.chordSaddle, sc[gov.k].s.axial.chordCrown, sc[gov.k].s.ipb.chordCrown, sc[gov.k].s.opb.chordSaddle] : [sc[gov.k].s.axial.braceSaddle, sc[gov.k].s.axial.braceCrown, sc[gov.k].s.ipb.braceCrown, sc[gov.k].s.opb.braceSaddle])), explain: 'Largest Efthymiou factor on the governing side of the weld. Fatigue damage grows with the third to fifth power of this number.' },
      ],
      checks: [
        check(`Fatigue damage — node ${gov.j.node}, brace ${gov.j.brace} (${side} side, point ${gov.p + 1})`, Math.min(fl.utilisation, 1e6), 1, 'Miner sum × DFF', `DNV-RP-C203 (2021) Eq. 2.2.1 with curve ${inp.curve} (${inp.envSN}), DFF = ${inp.dff}`, 'Accumulated fatigue damage at the weld toe over the service life with the design fatigue factor. Above 1.0 a through-thickness crack is predicted within the factored life; the brace then floods and sheds its load.'),
        fsCheck('Calculated life against service life × DFF', cap(gov.life), lifeReq, 'years', 'ISO 19902:2020 §16.12; NORSOK N-004 §8', 'The same criterion in years: the calculated life of the worst joint must exceed the service life multiplied by the design fatigue factor.'),
      ],
      plots: [
        { type: 'bar', title: 'Annual damage by wave block (governing joint)', categories: gb.map((b) => `${b.H} m`), series: [{ name: 'Damage per year', y: gd.map((v) => round(v, 3)) }], y: { label: 'Miner damage', unit: '1/yr' }, explain: 'Share of the yearly damage caused by each wave-height class. The peak normally sits at moderate waves of 2–5 m: frequent enough and large enough. Extreme waves contribute little.' },
        { type: 'line', title: 'S–N curve and wave-block stress ranges (governing hot spot)', x: { label: 'Cycles (log10)', unit: '-' }, y: { label: 'Hot-spot stress range (log10 MPa)', unit: '-' }, series: [{ name: `Curve ${inp.curve}`, x: sRange.map((s) => round(Math.log10(cyclesFromConstants(sn, s)))), y: sRange.map((s) => round(Math.log10(s / 1e6))) }, { name: 'Applied cycles per year', x: gb.filter((b) => b.r > 0).map((b) => round(Math.log10(b.n))), y: gb.filter((b) => b.r > 0).map((b) => round(Math.log10(b.r / 1e6))) }],
          explain: 'Design S–N curve of the detail (with thickness correction) and the yearly cycles of each wave block at their stress range, on log scales. The closer a point is to the curve, the more that block damages the joint.' },
        { type: 'line', title: 'Stress range around the weld in the largest wave block', x: { label: 'Point (1, 5 crown; 3, 7 saddle)', unit: '-' }, y: { label: 'Hot-spot stress range', unit: 'MPa' }, series: [{ name: 'Chord side', x: [1, 2, 3, 4, 5, 6, 7, 8], y: Array.from(rng[blocks.length - 1][gov.k][0], (v) => round(v / 1e6)) }, { name: 'Brace side', x: [1, 2, 3, 4, 5, 6, 7, 8], y: Array.from(rng[blocks.length - 1][gov.k][1], (v) => round(v / 1e6)) }],
          explain: 'Stress range at the eight points of the governing joint for the largest wave. Saddle peaks (3, 7) show axial and out-of-plane control, crown peaks (1, 5) in-plane bending.' },
        plot3d(J, { title: 'Fatigue utilisation of the joints', values: nodeU, label: 'DFF × damage over service life', unit: '-', explain: 'Factored Miner damage of the worst joint at each node (capped at 99). Critical joints are usually near the waterline, where wave loads on the braces are largest.' }),
      ],
      tables: [
        { title: 'Critical fatigue locations', columns: [{ key: 'n', label: 'Node' }, { key: 'z', label: 'Elevation', unit: 'm' }, { key: 'b', label: 'Brace' }, { key: 'side', label: 'Side' }, { key: 'pt', label: 'Point' }, { key: 'beta', label: 'β' }, { key: 'gamma', label: 'γ' }, { key: 'r', label: 'Largest range', unit: 'MPa' }, { key: 'd', label: 'Damage per year' }, { key: 'life', label: 'Life', unit: 'years' }, { key: 'u', label: 'DFF × damage' }],
          rows: top.slice(0, 20).map((r) => ({ n: r.j.node, z: round(r.j.z, 3), b: r.j.brace, side: r.q === 0 ? 'chord' : 'brace', pt: r.p + 1, beta: round(sc[r.k].p.beta, 3), gamma: round(sc[r.k].p.gamma, 3), r: round(r.rMax / 1e6), d: round(r.d, 3), life: round(cap(r.life), 3), u: round(Math.min((inp.damage0 + r.d * inp.life) * inp.dff, 1e6), 3) })),
          explain: 'The twenty joints with the shortest lives: location around the weld, largest stress range, yearly damage, life and factored damage over the service life. This is the inspection priority list.' },
        { title: 'Wave blocks at the governing hot spot', columns: [{ key: 'H', label: 'H', unit: 'm' }, { key: 'T', label: 'T', unit: 's' }, { key: 'n', label: 'Waves per year' }, { key: 'th', label: 'Wave theory' }, { key: 'r', label: 'Stress range', unit: 'MPa' }, { key: 'N', label: 'Cycles to failure' }, { key: 'd', label: 'Damage per year' }],
          rows: gb.map((b, i) => ({ H: b.H, T: b.T, n: round(b.n), th: b.theory, r: round(b.r / 1e6), N: Number.isFinite(b.N) ? round(b.N, 3) : '∞', d: round(gd[i], 3) })), explain: 'For each wave class: the hot-spot range it causes, the cycles the S–N curve allows at that range and the resulting yearly damage. The sum of the last column is the annual damage.' },
      ],
      notes: [`Governing: node ${gov.j.node} (z = ${fmt(gov.j.z)} m), brace ${gov.j.brace}, ${side} side point ${gov.p + 1}: life ${fmt(cap(gov.life))} years against ${fmt(lifeReq)} required (DFF ${fmt(inp.dff)} × ${fmt(inp.life)} years).`,
        `${nShort} of ${rows.length} joints are below the required life. Dynamic amplification ${fmt(inp.daf)}, S–N curve ${inp.curve} (${inp.envSN}), chord wall at leg joints ${fmt(gov.j.T * 1e3)} mm at the governing joint.`, spNote,
        sc.some((s) => s.warn.length) ? `${sc.filter((s) => s.warn.length).length} joint(s) are outside the Efthymiou validity range (e.g. ${sc.find((s) => s.warn.length).warn[0]}).` : 'All joints are inside the Efthymiou validity range.', foundNote(J)],
      advice: [
        ok ? adv('engineering', 'ok', `Fatigue requirement met: calculated life ${fmt(cap(gov.life))} years ≥ ${fmt(lifeReq)} years; factored damage ${fmt(fl.utilisation)}.`, fl.utilisation < 0.3 ? `There is margin: the governing can could go down to about ${fmt(tReq * 1e3)} mm, or the design fatigue factor could be raised to make the joint inspection-free. Confirm with all headings and their directional shares.` : 'Keep the can thickness. Add the other wave headings with their directional shares, and plan the first in-service inspection of the listed joints at about half the design fatigue life.', 'DNV-RP-C203 (2021) §2.2, §5; ISO 19902:2020 §16.12', `Fatigue life responds to stress to the power ${fmt(m)}: a 10 % lower SCF gives about ${fmt((1.1 ** m - 1) * 100)} % more life for almost no steel. One millimetre of can wall on one joint is about ${steelCo2(canPer * 1e-3)}`)
          : adv('engineering', 'fail', `Fatigue life is too short at node ${gov.j.node}: ${fmt(cap(gov.life))} years against ${fmt(lifeReq)} required (factored damage ${fmt(fl.utilisation)}); ${nShort} joint(s) affected.`, `Hot-spot stress must drop by ${fmt((1 - sRed) * 100)} %: thicken the chord can to about ${fmt(tReq * 1e3)} mm, raise β with a larger brace, or use weld-toe grinding (a factor of about 2 on life in DNV-RP-C203 §7) on the listed joints.`, 'DNV-RP-C203 (2021) §2.4, §7 and Appendix B', `Thickening the governing can by ${fmt(Math.max(tReq - gov.j.T, 0.001) * 1e3)} mm is about ${steelCo2(canPer * Math.max(tReq - gov.j.T, 0.001))} per joint — far less than an offshore repair clamp weighing tonnes.`),
        adv('economic', ok ? 'info' : 'warn', `${nShort} joint(s) below the required life; a design fatigue factor of ${fmt(inp.dff)} implies in-service inspection of the critical joints.`, 'Trade steel against inspection: a factor of 10 makes a joint inspection-free for life, a factor of 2–3 commits to periodic underwater inspection (diver or ROV days every campaign). For the few governing joints a thicker can at fabrication is cheaper than decades of inspection.', `${BASIS_COST}; inspection planning per DNV-RP-C210 (2019) and ISO 19902:2020 §23`),
      ],
      publish: { 'life.fatigueYears': cap(fl.designLife) },
      quantities: { steel_t: round(steel) },
    };
  },
  verify: [
    // Horizontal brace stub (0.5 m × 12.5 mm, 5 m cantilever) at z = −10 m in 20 m of water, wave travelling normal to it:
    // Airy, H = 6 m, T = 8 s → k = 0.07082 rad/m; horizontal acceleration amplitude a = (H/2)ω²·cosh(k·10)/sinh(k·20) = 1.2043 m/s².
    // Inertia load (Cd = 0, Cm = 2) q = 1025·2·0.19635·1.2043 = 484.7 N/m → out-of-plane moment q·L²/2 = 6059 N·m,
    // W = 2.2766e-3 m³ → σ = 2.661 MPa. Chord saddle out-of-plane SCF (Efthymiou Eq. 10, β = 0.5, γ = 20, τ = 0.5, θ = 90°)
    // = γτβ(1.7 − 1.05β³) = 7.844 → range = 2·7.844·2.661 = 41.75 MPa. T curve, seawater with CP, second slope, with the
    // thickness correction for the 25 mm chord (reference 16 mm, k = 0.25): 41.75·(25/16)^0.25 = 46.68 MPa,
    // N = 10^16.13/46.68⁵ = 6.09e7 → with 1e6 waves per year the life is 60.9 years.
    { name: 'Single wave block on a cantilever brace: Δσ = 2·SCF·M/W, N = a·Δσ⁻ᵐ', source: 'Morison inertia load in Airy kinematics; Efthymiou OPB SCF (DNV-RP-C203 Table B-1 Eq. 10); DNV-RP-C203 Table 2-2 curve T (hand calculation)',
      inputs: { geom: testFrame('tjoint', { thetaDeg: 90 }), depth: 20, heading: 90, dirFrac: 1, theory: 'airy', Cd: 0, Cm: 2, mg: 0, kinFac: 1, appF: 1, daf: 1, canT: 0, found: 'fixed', blocks: [{ H: 6, T: 8, n: 1e6 }], scatter: [], nPhase: 24, nSub: 1 },
      expect: [{ kpi: 'maxRange', value: 41.75, tol: 0.01 }, { kpi: 'fatigueLife', value: 60.9, tol: 0.05 }] },
  ],
};
