// Fatigue and fracture mechanics of a pipeline girth weld: S–N damage from a stress-range histogram, an engineering
// critical assessment of a surface flaw on the failure assessment diagram, and Paris-law growth of the flaw to failure.
import { snCurve, cyclesToFailure, stressRangeAtCycles } from '../../engine/fatigue/sn.js';
import { parisLawBs7910, surfaceCrackGrowth, assessSurfaceCrack, criticalCrackSize, ctodToK, fadOption1 } from '../../engine/fatigue/fracture.js';
import { pipeInputs, pipeFrom, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, linspace, sev, worst,
  steelPerMetre, vesselCo2, CONFIRM, RBI_BASIS, CO2_STEEL, MPA } from './_common.js';

const G_L = 'Fatigue loading';
const G_S = 'S–N assessment';
const G_F = 'Flaw and fracture';
const G_C = 'Acceptance criteria';
const HIST = [{ range: 80, cycles: 5 }, { range: 50, cycles: 60 }, { range: 30, cycles: 1000 }, { range: 15, cycles: 2e4 }, { range: 8, cycles: 2e5 }];

export default {
  id: 's08-fatigue-fracture',
  suite: 8,
  title: 'Fatigue and fracture of girth welds',
  summary: 'Gives the fatigue life of the girth weld under its stress history and tells whether a weld flaw of given size is acceptable and how long it takes to grow through the wall.',
  covers: [12],
  fidelity: 'analytical',
  standards: ['DNV-RP-C203 (2021) S–N curves and Palmgren–Miner rule', 'BS 7910:2019 clause 7 (Option 1 FAD), clause 8 (fatigue crack growth), Annex M, Annex P', 'DNV-ST-F101 (2021) §5.4.8 (fatigue), Appendix A (ECA)', 'DNV-RP-F108 (2021) (ECA of girth welds)'],
  models: ['S–N curves with thickness and stress-concentration correction', 'Palmgren–Miner damage summation', 'Paris–Erdogan crack propagation (BS 7910 laws)', 'Newman–Raju stress intensity for a semi-elliptical surface crack', 'Failure assessment diagram, Option 1', 'CTOD to fracture toughness conversion'],
  limits: 'The flaw is assessed as a semi-elliptical surface crack in a flat plate of the wall thickness (no curvature or pipe-specific reference stress, no weld-toe magnification unless folded into the SCF), under tension and through-wall bending. Crack growth uses one equivalent constant-amplitude stress range with the growth threshold set to zero (conservative for spectra with many small cycles). Strain-based ECA for reeling and tearing analysis are not included. The design fatigue factor and the residual stress assumption are recalled defaults to be confirmed; sour or hydrogen service needs environment-specific crack growth data.',
  publishes: ['life.fatigueYears'],
  searchTerms: 'pipeline girth weld fatigue S-N engineering critical assessment BS 7910 Paris law',
  inputs: [
    ...pipeInputs({ omit: ['tCoat', 'rhoCoat', 'tConc', 'concrete', 'rhoCont', 'rhoW', 'wsGiven'] }),
    { key: 'hist', label: 'Annual stress-range histogram', type: 'table', group: G_L, role: 'load',
      columns: [{ key: 'range', label: 'Nominal axial stress range', unit: 'MPa', type: 'number' }, { key: 'cycles', label: 'Cycles per year', type: 'number' }],
      default: HIST, help: 'Blocks of nominal longitudinal stress range at the weld from VIV, wave loading, slugging or pressure and temperature cycling (e.g. from the free-span or slug-load results).' },
    num('scale', 'Scale factor on the stress ranges', '-', 1, 0.1, 5, G_L, 'load', 'Multiplies every block; use it for sensitivity or to apply a load factor.', { uncertain: { dist: 'lognormal', cov: 0.15 } }),
    num('life', 'Design life', 'years', 25, 1, 100, G_L, 'operation', 'Service life over which damage accumulates.', { bind: 'life.designYears' }),
    sel('sn', 'S–N curve', 'F1', ['C1', 'D', 'E', 'F', 'F1', 'F3'], G_S, 'material', 'DNV-RP-C203 class: D for ground-flush or high-quality welds, E/F cap, F1/F3 single-sided root.'),
    sel('env', 'Environment', 'seawater-cp', [['air', 'Air / dry internal'], ['seawater-cp', 'Sea water with cathodic protection'], ['free-corrosion', 'Free corrosion']], G_S, 'environment', 'Environment at the crack initiation site; sets the S–N curve and the crack growth law.'),
    num('scf', 'Stress concentration factor', '-', 1.2, 1, 4, G_S, 'material', 'Misalignment and local geometry not already in the S–N class.', { uncertain: { dist: 'lognormal', cov: 0.1 } }),
    num('dff', 'Design fatigue factor', '-', 6, 1, 10, G_C, 'criteria', `Allowable damage = 1/DFF. DNV-OS-F101 (2013) Table 5-11, carried into DNV-ST-F101: 3 / 6 / 10 for low / medium / high safety class; the default is the medium class. ${CONFIRM}`),
    num('a0', 'Flaw depth', 'mm', 2, 0.1, 60, G_F, 'ic', 'Height of the surface-breaking planar flaw, e.g. the inspection acceptance or the sized indication.', { uncertain: { dist: 'lognormal', cov: 0.2 } }),
    num('c0', 'Flaw half-length', 'mm', 15, 0.2, 500, G_F, 'ic', 'Half of the surface length of the flaw.'),
    num('Pm', 'Maximum membrane stress', 'MPa', 250, 0, 700, G_F, 'load', 'Largest longitudinal membrane stress in the extreme load case (installation or operation).', { uncertain: { dist: 'normal', cov: 0.1 } }),
    num('Pb', 'Maximum bending stress through the wall', 'MPa', 20, 0, 400, G_F, 'load', 'Through-wall bending from misalignment under the extreme load.'),
    sel('residual', 'Welding residual stress', 'yield', [['yield', 'As-welded: yield magnitude'], ['pwht', 'Stress relieved: 30 % of yield'], ['none', 'None']], G_F, 'ic', `Secondary membrane stress acting on the flaw (BS 7910 clause 7.1.8). ${CONFIRM}`),
    num('ctod', 'Fracture toughness CTOD', 'mm', 0.25, 0.01, 3, G_F, 'material', 'Characteristic crack-tip opening displacement of the weld metal / HAZ at minimum design temperature.', { uncertain: { dist: 'lognormal', cov: 0.3 } }),
    num('rBend', 'Bending fraction of the stress range', '-', 0.1, 0, 1, G_F, 'load', 'Part of each stress range acting as through-wall bending in the crack growth analysis.'),
    num('sfLife', 'Safety factor on crack-growth life', '-', 3, 1, 20, G_C, 'criteria', `Required ratio of calculated growth life to design life. ${CONFIRM}`),
  ],
  compute(inp, ctx) {
    const p = pipeFrom(inp, ctx);
    const t = p.t;
    const rows = (inp.hist || []).map((r) => ({ S: Number(r.range) * inp.scale * MPA, n: Number(r.cycles) })).filter((r) => r.S > 0 && r.n > 0);
    need(rows.length > 0, 'The stress-range histogram is empty: enter at least one block with a positive range and cycle count.');
    const a0 = inp.a0 / 1e3, c0 = inp.c0 / 1e3;
    need(a0 < 0.9 * t, `The flaw depth (${fmt(inp.a0)} mm) must be less than 90 % of the wall (${fmt(t * 1e3)} mm).`);
    need(c0 >= a0 * 0.5, 'The flaw half-length must be at least half its depth (a/c ≤ 2) for the Newman–Raju solution.');
    // S–N
    const curve = snCurve(inp.sn, inp.env);
    const blocks = rows.map((r) => { const N = cyclesToFailure({ curve, stressRange: r.S, thickness: t, scf: inp.scf }); return { ...r, N, d: r.n / N }; });
    const dYear = blocks.reduce((s, b) => s + b.d, 0);
    const lifeSN = dYear > 0 ? 1 / dYear : 1e9;
    const ucSN = dYear * inp.life * inp.dff;
    const nTot = rows.reduce((s, r) => s + r.n, 0);
    // FAD
    const Qm = inp.residual === 'yield' ? p.smys : inp.residual === 'pwht' ? 0.3 * p.smys : 0;
    const Kmat = ctodToK({ ctod: inp.ctod / 1e3, yieldStrength: p.smys, tensileStrength: p.smts, E: p.E, nu: p.nu }).K;
    const fadArgs = { t, Pm: inp.Pm * MPA, Pb: inp.Pb * MPA, Qm, yieldStrength: p.smys, tensileStrength: p.smts, E: p.E, toughness: Kmat };
    const fad = run(() => assessSurfaceCrack({ ...fadArgs, a: a0, c: c0 }), 'Failure assessment');
    const critSize = run(() => criticalCrackSize({ ...fadArgs, aspectRatio: a0 / c0 }), 'Critical flaw size');
    // crack growth under the equivalent constant-amplitude range (m = 3 weighting)
    const env = inp.env === 'air' ? 'air' : inp.env === 'seawater-cp' ? 'marine-cp' : 'marine';
    // the histogram is condensed to one equivalent range, so the growth threshold is set to zero: otherwise the small equivalent
    // range would arrest a flaw that the larger blocks do in fact drive
    const law = parisLawBs7910({ environment: env, threshold: 0 });
    const m = law.stages[0].m;
    const Seq = (rows.reduce((s, r) => s + r.n * r.S ** m, 0) / nTot) ** (1 / m) * inp.scf;
    const gr = run(() => surfaceCrackGrowth({ a0, c0, t, membraneRange: Seq * (1 - inp.rBend), bendingRange: Seq * inp.rBend, law, toughness: Kmat, peakFactor: Math.max((inp.Pm + inp.Pb) * MPA / Math.max(Seq, 1), 1), maxCycles: nTot * 2000 }), 'Crack growth');
    const lifeCG = Number.isFinite(gr.cycles) ? gr.cycles / nTot : 1e9;
    const ucCG = Math.min(inp.life * inp.sfLife / lifeCG, 99);
    // trajectory of the assessment point as the flaw grows
    const hist = gr.history;
    const idx = hist.a.length > 12 ? Array.from({ length: 12 }, (_, i) => Math.round(i * (hist.a.length - 1) / 11)) : hist.a.map((_, i) => i);
    const traj = idx.map((i) => { try { const f = assessSurfaceCrack({ ...fadArgs, a: Math.min(hist.a[i], 0.98 * t), c: hist.c[i] }); return { Lr: f.Lr, Kr: f.Kr, yr: hist.cycles[i] / nTot, a: hist.a[i] }; } catch { return null; } }).filter(Boolean);
    const lr = linspace(0, fad.LrMax, 40);
    const fadLine = lr.map((x) => fadOption1(x, { yieldStrength: p.smys, tensileStrength: p.smts, E: p.E }).Kr);
    const checks = [
      chk('S–N fatigue damage with design fatigue factor', dYear * inp.life, 1 / inp.dff, '-', 'DNV-RP-C203 (2021) §2.2; DNV-ST-F101 (2021) §5.4.8',
        'Guards against fatigue crack initiation and growth at the girth weld over the design life. Exceeding it means the weld is expected to crack before the factored life.'),
      chk('Fracture of the flaw under the extreme load (FAD)', 1, Math.max(fad.margin, 1e-6), '-', 'BS 7910:2019 clause 7, Option 1',
        'Guards against unstable fracture or plastic collapse from the flaw under the maximum stress and residual stress. A point outside the diagram means the flaw is not acceptable.'),
      chk('Crack growth life with safety factor', inp.life * inp.sfLife, lifeCG, 'years', 'BS 7910:2019 clause 8; DNV-RP-F108',
        'Guards against the flaw growing by fatigue to through-wall or to the critical size within the factored design life, which would cause a leak or fracture.'),
    ];
    const gov = worst(checks);
    const sAllow = stressRangeAtCycles({ curve, cycles: Math.max(nTot * inp.life * inp.dff, 1e4), thickness: t, scf: inp.scf });
    const sn = linspace(4, 9, 26).map((e) => 10 ** e);
    const snS = sn.map((N) => stressRangeAtCycles({ curve, cycles: N, thickness: t, scf: 1 }) / MPA);
    const scaleOk = ucSN > 0 ? (1 / ucSN) ** (1 / 3) : 99; // damage ∝ S^3 on the first slope
    return {
      kpis: [
        kpi('lifeSN', 'S–N fatigue life (unfactored)', Math.min(lifeSN, 1e9), 'years', 'One over the annual Miner damage. It must exceed the design life times the design fatigue factor.', inp.life * inp.dff, 'min'),
        kpi('damage', 'Miner damage over the design life', dYear * inp.life, '-', 'Accumulated damage without the design fatigue factor; the allowable is 1/DFF.', 1 / inp.dff, 'max'),
        kpi('ucSN', 'S–N utilisation', ucSN, '-', 'Damage × DFF. Above 1.0 the weld does not meet the fatigue requirement.', 1, 'max'),
        kpi('Kr', 'Fracture ratio K_r', fad.Kr, '-', 'Applied stress intensity (primary plus secondary) over the fracture toughness, with the plasticity interaction.'),
        kpi('Lr', 'Load ratio L_r', fad.Lr, '-', 'Reference stress over yield: proximity to plastic collapse of the remaining ligament.', fad.LrMax, 'max'),
        kpi('fadMargin', 'FAD margin', Math.min(fad.margin, 99), '-', 'Allowable K_r at the assessed L_r over the actual K_r. At least 1.0 for the flaw to be acceptable.', 1, 'min'),
        kpi('reserve', 'Reserve factor on primary stress', Math.min(fad.reserveFactor, 99), '-', 'Factor by which the primary stresses can rise before the point reaches the assessment line.'),
        kpi('Kmat', 'Fracture toughness K_mat', Kmat / 1e6, 'MPa√m', 'Toughness converted from the CTOD value.'),
        kpi('aCrit', 'Critical flaw depth', critSize.a * 1e3, 'mm', 'Deepest flaw of the same aspect ratio that stays inside the diagram under the extreme load; the wall thickness when even a through-wall flaw is stable.'),
        kpi('lifeCG', 'Crack growth life', Math.min(lifeCG, 1e9), 'years', `Time for the flaw to grow from its initial size to ${gr.status} under the equivalent stress range.`, inp.life * inp.sfLife, 'min'),
        kpi('Seq', 'Equivalent stress range at the weld', Seq / MPA, 'MPa', 'Constant-amplitude range giving the same damage as the histogram for the Paris exponent, including the SCF.'),
        kpi('nYear', 'Cycles per year', nTot, '-', 'Total number of stress cycles per year in the histogram.'),
        kpi('sAllow', 'Allowable constant stress range', sAllow / MPA, 'MPa', 'Constant nominal stress range that would just use the allowable damage with this number of cycles.'),
      ],
      checks,
      plots: [
        plot('line', 'S–N curve and loading blocks', 'log10 cycles', '-', 'Stress range', 'MPa', [
          { name: `${inp.sn} (${inp.env}), no SCF`, x: sn.map((N) => round(Math.log10(N))), y: snS.map((v) => round(v)) },
          { name: 'Blocks × SCF over the design life', x: blocks.map((b) => round(Math.log10(b.n * inp.life))), y: blocks.map((b) => round(b.S * inp.scf / MPA)) },
        ], 'The S–N curve with the histogram blocks plotted at their local stress range and life-time cycle count. Blocks close to the curve dominate the damage.'),
        plot('line', 'Failure assessment diagram', 'L_r', '-', 'K_r', '-', [
          { name: 'Option 1 assessment line', x: lr.map((v) => round(v)), y: fadLine.map((v) => round(v)) }, { name: 'Flaw today', x: [round(fad.Lr)], y: [round(fad.Kr)] },
          { name: 'Path as the flaw grows', x: traj.map((q) => round(q.Lr)), y: traj.map((q) => round(q.Kr)) },
        ], 'The flaw is acceptable while its point lies inside the line. The path shows the point moving outwards as fatigue deepens the flaw; where it crosses the line the flaw becomes critical.'),
        plot('line', 'Flaw depth vs time', 'Time', 'years', 'Flaw depth', 'mm', [
          { name: 'Depth', x: Array.from(hist.cycles, (c) => round(c / nTot)), y: Array.from(hist.a, (a) => round(a * 1e3)) }, { name: 'Wall thickness', x: [0, round(Math.min(lifeCG, 1e6))], y: [round(t * 1e3), round(t * 1e3)], dash: true },
        ], 'Growth of the flaw depth by fatigue. Growth accelerates as the flaw deepens, so most of the life is spent while the flaw is small.'),
      ],
      tables: [tbl('Damage by stress block', [col('S', 'Nominal range', 'MPa'), col('n', 'Cycles per year'), col('N', 'Cycles to failure'), col('d', 'Damage per year'), col('share', 'Share', '%')],
        blocks.map((b) => ({ S: round(b.S / MPA), n: b.n, N: round(Math.min(b.N, 1e15), 4), d: round(b.d, 4), share: round(dYear > 0 ? 100 * b.d / dYear : 0, 3) })),
        'Each histogram block with its allowable cycles and its share of the annual damage; the block with the largest share is the one to reduce.')],
      notes: [`S–N: ${inp.sn} in ${inp.env}, thickness ${fmt(t * 1e3)} mm, SCF ${fmt(inp.scf)}; annual damage ${fmt(dYear)}.`,
        `ECA: residual stress ${fmt(Qm / MPA)} MPa, K_mat ${fmt(Kmat / 1e6)} MPa√m; governing point at the ${fad.governingPoint} of the flaw; critical depth limited by ${critSize.limitedBy}.`,
        `Crack growth with the ${law.label} law ends by ${gr.status} after ${fmt(lifeCG)} years. ${gov.name} governs (${fmt(gov.uc)}).`],
      advice: [
        adv('engineering', sev(gov.uc), `S–N life ${fmt(lifeSN)} years (needs ${fmt(inp.life * inp.dff)}), FAD margin ${fmt(fad.margin)}, crack growth life ${fmt(lifeCG)} years (needs ${fmt(inp.life * inp.sfLife)}); governing utilisation ${fmt(gov.uc)}.`,
          gov.uc > 1 ? (gov.name.startsWith('S–N') ? `Reduce the stress ranges by about ${fmt((1 - scaleOk) * 100)} % (span correction, slug mitigation) or improve the weld to a higher S–N class.` : gov.name.startsWith('Fracture') ? `Tighten the flaw acceptance to below ${fmt(critSize.a * 1e3)} mm depth or specify a higher CTOD.` : 'Reduce the initial flaw size accepted at inspection or the cyclic stress; growth life scales with the cube of the stress range.')
            : `The weld meets all three criteria. The largest flaw depth acceptable at the extreme load is ${fmt(critSize.a * 1e3)} mm; set the inspection acceptance criteria below it with sizing error.`,
          'DNV-RP-C203 (2021); BS 7910:2019 clauses 7–8; DNV-ST-F101 (2021) App. A — DFF and residual stress recalled, confirm'),
        adv('economic', 'info', `The fatigue design uses ${fmt(Math.min(ucSN, 99) * 100)} % of the allowable damage; the allowable constant stress range is ${fmt(sAllow / MPA)} MPa.`,
          ucSN > 1 ? 'Weld-quality upgrades (root grinding, tighter high–low limits) cost hours per joint on the lay vessel; compare with removing the load at source, which usually costs less across thousands of welds.' : `Automated ultrasonic acceptance set from the critical depth of ${fmt(critSize.a * 1e3)} mm avoids unnecessary repairs; each avoided cut-out saves roughly an hour of lay-vessel time.`,
          RBI_BASIS, `A girth-weld failure means a repair clamp or spool replacement: about ${fmt(steelPerMetre(p, p.t) * 12 / 1e3)} t of pipe (${fmt(steelPerMetre(p, p.t) * 12 / 1e3 * CO2_STEEL)} t CO₂e, screening factor) and several vessel days. ${vesselCo2(5)}`),
      ],
      publish: { 'life.fatigueYears': Math.min(lifeSN, lifeCG, 1e9) },
    };
  },
  verify: [{
    // DNV-RP-C203 curve D in air: log N = 12.164 − 3 log Δσ. At Δσ = 100 MPa, N = 10^(12.164 − 6) = 1.4588e6 cycles
    // (t = 14.3 mm < 25 mm: no thickness correction). 1e5 cycles per year → damage 0.06855 per year → life 14.59 years.
    name: 'Single stress block on the D curve in air',
    source: 'DNV-RP-C203 (2021) Table 2-1; hand calculation',
    inputs: { D: 323.9, t: 14.3, hist: [{ range: 100, cycles: 1e5 }], scale: 1, sn: 'D', env: 'air', scf: 1 },
    expect: [{ kpi: 'lifeSN', value: 14.59, tol: 0.003 }],
  }],
};
