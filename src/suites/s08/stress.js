// Global and local stress analysis of a pipe section: thick-wall (Lamé) and thin-wall pressure stresses, axial stress
// for restrained and unrestrained conditions, bending from route curvature or an applied moment, torsion, von Mises and
// Tresca equivalent stress, and the DNV allowable-stress and ASME B31.8 stress criteria.
import { equivalentStressCheck, asmeStressCheck, lameStresses, bendingStrainRadius } from '../../engine/pipeline/strength.js';
import { restrainedAxialForce } from '../../engine/pipeline/expansion.js';
import { characteristicStrength } from '../../engine/pipeline/pipeprops.js';
import { pipeInputs, pipeFrom, depthInput, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, linspace, sev, worst,
  steelPerMetre, steelCo2, CONFIRM, LIMIT_NO_FE, ECON_BASIS, G, BAR, MPA } from './_common.js';

const G_L = 'Loads and restraint';
const G_C = 'Acceptance criteria';

export default {
  id: 's08-stress',
  suite: 8,
  title: 'Global and local pipe stress',
  summary: 'Combines pressure, temperature, axial restraint, bending and torsion into wall stresses and checks them against the code allowables.',
  covers: [2],
  fidelity: 'analytical',
  standards: ['DNV-ST-F101 (2021) Sec. 5 (allowable stress design format)', 'ISO 13623 (2017) §6.4', 'ASME B31.8 (2022) §833, §A842.2', 'ASME B31.4 (2022) §402–403'],
  models: ['Lamé thick-wall cylinder equations', 'Thin-wall (Barlow) hoop stress', 'Euler–Bernoulli beam bending', 'Restrained thermal expansion and Poisson effect', 'von Mises yield criterion', 'Tresca yield criterion'],
  limits: `Linear-elastic section analysis at one cross-section: no plasticity, ovalisation or residual stress from manufacture. The usage factor η and the ASME design factors are recalled defaults to be confirmed. ${LIMIT_NO_FE}`,
  publishes: ['pipeline.vonMises'],
  searchTerms: 'pipeline stress analysis von Mises hoop longitudinal restrained ASME B31.8',
  inputs: [
    ...pipeInputs({ omit: ['tCoat', 'rhoCoat', 'tConc', 'concrete', 'rhoCont', 'wsGiven'] }),
    num('tCorr', 'Corrosion allowance deducted', 'mm', 3, 0, 15, 'Pipe', 'geometry', 'Stresses are evaluated on the corroded wall.'),
    depthInput(120, G_L),
    num('pi', 'Internal pressure', 'bar', 200, 0, 1500, G_L, 'bc', 'Internal pressure at the section (design or operating, as the check requires).', { uncertain: { dist: 'normal', cov: 0.03 } }),
    num('Tinst', 'Installation temperature', '°C', 5, -10, 40, G_L, 'ic', 'Temperature at which the pipe was laid and restrained; reference for the thermal strain.', { bind: 'site.Tbed' }),
    num('Top', 'Operating temperature', '°C', 60, -20, 200, G_L, 'bc', 'Wall temperature at the section in operation.', { bind: 'pipeline.Tout', uncertain: { dist: 'normal', cov: 0.05 } }),
    sel('restraint', 'Axial restraint', 'restrained', [['restrained', 'Fully restrained (buried or far from the ends)'], ['unrestrained', 'Unrestrained with end cap (near a free end)'], ['force', 'Effective axial force specified']], G_L, 'bc',
      'Restrained sections carry the thermal and Poisson forces; free sections carry only the end-cap force.'),
    num('Seff', 'Effective axial force (tension positive)', 'kN', -1500, -50000, 50000, G_L, 'bc', 'Used when the effective force is specified, e.g. from the expansion module.', { bind: 'pipeline.effectiveForce', showIf: { key: 'restraint', in: ['force'] } }),
    num('Hlay', 'Residual lay tension', 'kN', 100, 0, 5000, G_L, 'ic', 'Effective tension locked in at installation (restrained case).', { showIf: { key: 'restraint', in: ['restrained'] } }),
    num('M', 'Bending moment', 'kN·m', 20, 0, 20000, G_L, 'bc', 'Bending moment at the section from spans, crossings or supports.', { uncertain: { dist: 'normal', cov: 0.1 } }),
    num('Rbend', 'Radius of route curvature (0 = straight)', 'm', 3000, 0, 100000, G_L, 'geometry', 'Elastic bending to a route curve or seabed undulation adds M = EI/R.'),
    num('Tq', 'Torque', 'kN·m', 0, 0, 5000, G_L, 'bc', 'Torsion about the pipe axis.'),
    num('eta', 'Usage factor η on yield', '-', 0.87, 0.5, 1, G_C, 'criteria', `Allowable equivalent stress as a fraction of yield (code and location dependent, 0.72–1.0). ${CONFIRM}`),
    num('F1', 'ASME hoop design factor', '-', 0.72, 0.4, 0.8, G_C, 'criteria', `ASME B31.8 Table A842.2.2-1. ${CONFIRM}`),
    num('F2', 'ASME longitudinal design factor', '-', 0.8, 0.5, 1, G_C, 'criteria', `ASME B31.8 Table A842.2.2-1. ${CONFIRM}`),
    num('F3', 'ASME combined design factor', '-', 0.9, 0.5, 1, G_C, 'criteria', `ASME B31.8 Table A842.2.2-1. ${CONFIRM}`),
  ],
  compute(inp, ctx) {
    const p0 = pipeFrom(inp, ctx);
    const t = p0.t - inp.tCorr / 1e3;
    need(t > 0.002, 'The corrosion allowance consumes the whole wall.');
    const p = pipeFrom({ ...inp, t: t * 1e3 }, ctx);
    const pe = p.rhoW * G * inp.depth, pi = inp.pi * BAR;
    const dT = inp.Top - inp.Tinst;
    const fy = characteristicStrength({ smys: p.smys, smts: p.smts, tempC: Math.max(inp.Top, 20), alphaU: 1 }).fy;
    // effective force S and true wall force N = S + p_i A_i − p_e A_e
    let S;
    if (inp.restraint === 'restrained') S = restrainedAxialForce({ H: inp.Hlay * 1e3, deltaPi: pi, Ai: p.Ai, As: p.As, deltaT: dT, E: p.E, alpha: p.alpha, nu: p.nu }).S0;
    else if (inp.restraint === 'force') S = inp.Seff * 1e3;
    else S = 0;
    const N = S + pi * p.Ai - pe * p.Ae;
    const Mc = inp.Rbend > 0 ? p.EI / inp.Rbend : 0;
    const M = inp.M * 1e3 + Mc, T = inp.Tq * 1e3;
    const eq = run(() => equivalentStressCheck({ D: p.D, t, fy, eta: inp.eta, pi, pe, N, M, T }), 'Equivalent stress');
    const lam = lameStresses({ D: p.D, t, pi, pe });
    const sA = N / p.As, sB = M * p.D / (2 * p.I), tau = T * p.D / (2 * 2 * p.I);
    // three-dimensional von Mises and Tresca through the wall on the tension (+) and compression (−) bending fibres
    const ri = p.D / 2 - t, ro = p.D / 2;
    const vm3 = (sh, sl, sr, tq) => Math.sqrt(0.5 * ((sh - sl) ** 2 + (sl - sr) ** 2 + (sr - sh) ** 2) + 3 * tq * tq);
    const tresca = (sh, sl, sr, tq) => {
      const c = 0.5 * (sh + sl), r = Math.hypot(0.5 * (sh - sl), tq);
      const pr = [c + r, c - r, sr];
      return Math.max(...pr) - Math.min(...pr);
    };
    const rs = linspace(ri, ro, 21);
    const wall = rs.map((r) => {
      const L = lameStresses({ D: p.D, t, pi, pe, r });
      const fb = sB * r / ro, tq = tau * r / ro;
      return { r, sh: L.sigmaH, sr: L.sigmaR, vmT: vm3(L.sigmaH, sA + fb, L.sigmaR, tq), vmC: vm3(L.sigmaH, sA - fb, L.sigmaR, tq), trT: tresca(L.sigmaH, sA + fb, L.sigmaR, tq), trC: tresca(L.sigmaH, sA - fb, L.sigmaR, tq) };
    });
    const vmMax = Math.max(...wall.map((w) => Math.max(w.vmT, w.vmC)));
    const trMax = Math.max(...wall.map((w) => Math.max(w.trT, w.trC)));
    const asme = run(() => asmeStressCheck({ D: p.D, t, smys: p.smys, pi, pe, F1: inp.F1, F2: inp.F2, F3: inp.F3, restrained: inp.restraint !== 'unrestrained', deltaT: inp.restraint === 'unrestrained' ? 0 : dT, E: p.E, alpha: p.alpha, nu: p.nu,
      axialForce: inp.restraint === 'restrained' ? inp.Hlay * 1e3 : 0, M, T, criterion: 'tresca' }), 'ASME stress check');
    const allow = inp.eta * fy;
    const checks = [
      chk('Equivalent stress (thin wall, von Mises)', eq.sigmaE, eq.allowable, 'Pa', 'DNV-ST-F101 (2021) Sec. 5 ASD format; ISO 13623 §6.4.2',
        'Guards against first yield of the wall under the combined hoop, longitudinal and shear stress. Exceeding the allowable means permanent deformation under the stated loads.'),
      chk('Equivalent stress (thick wall, Lamé + von Mises)', vmMax, allow, 'Pa', 'Lamé (1852) with the von Mises criterion',
        'Same yield check with the exact through-wall pressure stresses; it governs thick pipe where the bore hoop stress exceeds the thin-wall value.'),
      chk('Hoop stress (ASME)', asme.SH, asme.allowHoop, 'Pa', 'ASME B31.8 §841.1.1 / §A842.2.2',
        'Guards against hoop yielding under the pressure difference with the design factor F1.'),
      chk('Longitudinal stress (ASME)', Math.max(Math.abs(asme.SLmax), Math.abs(asme.SLmin)), asme.allowLong, 'Pa', 'ASME B31.8 §833.3 / §A842.2.2',
        'Guards against the axial plus bending stress exceeding the longitudinal allowable; in restrained hot lines it is usually compressive.'),
      chk('Combined stress (ASME, Tresca)', asme.Scombined, asme.allowCombined, 'Pa', 'ASME B31.8 §833.4 / §A842.2.2',
        'Guards against yield under the combination of hoop and longitudinal stress with the design factor F3.'),
    ];
    const gov = worst(checks);
    // allowable additional bending moment by the thin-wall equivalent-stress criterion
    const ms = linspace(0, Math.max(2 * M, 0.9 * fy * p.Z), 25);
    const eqM = ms.map((m) => equivalentStressCheck({ D: p.D, t, fy, eta: inp.eta, pi, pe, N, M: m, T }).sigmaE);
    let Mallow = ms[ms.length - 1];
    for (let i = 1; i < ms.length; i++) if (eqM[i] > allow) { Mallow = eqM[i - 1] > allow ? 0 : ms[i - 1] + (ms[i] - ms[i - 1]) * (allow - eqM[i - 1]) / (eqM[i] - eqM[i - 1]); break; }
    const ang = linspace(0, 360, 37);
    const vmAng = ang.map((a) => vm3(lam.sigmaHOuter, sA + sB * Math.cos(a * Math.PI / 180), -pe, tau));
    const rElastic = bendingStrainRadius({ D: p.D, sigmaAllow: Math.max(allow - Math.abs(sA), 1e6), E: p.E }).minRadiusElastic;
    // wall thickness at which the governing thin-wall utilisation reaches 1 (stress ∝ 1/t to first order)
    const tNeed = t * Math.max(eq.uc, vmMax / allow) + inp.tCorr / 1e3;
    return {
      kpis: [
        kpi('vonMises', 'Maximum von Mises stress', vmMax / MPA, 'MPa', 'Largest equivalent stress anywhere through the wall on the bending tension or compression side. It must stay below η times yield.', allow / MPA, 'max'),
        kpi('ucEq', 'Equivalent stress utilisation', Math.max(eq.uc, vmMax / allow), '-', 'von Mises stress over the allowable. Above 1.0 the section yields under the stated combination.', 1, 'max'),
        kpi('tresca', 'Maximum Tresca stress', trMax / MPA, 'MPa', 'Twice the maximum shear stress; equal to or up to 15 % above the von Mises value.'),
        kpi('hoopThin', 'Hoop stress, thin wall', eq.sigmaH / MPA, 'MPa', 'Mid-wall hoop stress (p_i − p_e)(D − t)/(2t) used by the design codes.'),
        kpi('hoopLame', 'Hoop stress at the bore, Lamé', lam.sigmaHInner / MPA, 'MPa', 'Exact hoop stress at the inner surface, the highest in the wall for internal overpressure.'),
        kpi('sigmaAxial', 'Axial wall stress', sA / MPA, 'MPa', 'True wall force over the steel area: thermal, Poisson, end-cap and lay-tension contributions. Negative is compression.'),
        kpi('sigmaBend', 'Bending stress at the outer fibre', sB / MPA, 'MPa', 'Stress from the applied moment plus the elastic route curvature.'),
        kpi('Seff', 'Effective axial force', S / 1e3, 'kN', 'Effective force (true wall force corrected for pressure). Negative values drive global buckling.'),
        kpi('Ntrue', 'True wall force', N / 1e3, 'kN', 'Axial force actually carried by the steel, which produces the axial stress.'),
        kpi('Mallow', 'Allowable bending moment', Mallow / 1e3, 'kN·m', 'Moment at which the thin-wall equivalent stress reaches the allowable with the pressure and axial force unchanged.'),
        kpi('ucAsme', 'ASME governing utilisation', Math.max(asme.ucHoop, asme.ucLong, asme.ucCombined), '-', 'Largest of the ASME hoop, longitudinal and combined stress utilisations.', 1, 'max'),
        kpi('Rmin', 'Minimum elastic bend radius', rElastic, 'm', 'Tightest route curve the pipe can follow elastically with the bending stress left after the axial stress is deducted.'),
      ],
      checks,
      plots: [
        plot('line', 'Stress through the wall', 'Radius', 'mm', 'Stress', 'MPa', [
          { name: 'Hoop (Lamé)', x: rs.map((r) => round(r * 1e3)), y: wall.map((w) => round(w.sh / MPA)) }, { name: 'Radial (Lamé)', x: rs.map((r) => round(r * 1e3)), y: wall.map((w) => round(w.sr / MPA)) },
          { name: 'von Mises, tension side', x: rs.map((r) => round(r * 1e3)), y: wall.map((w) => round(w.vmT / MPA)) }, { name: 'von Mises, compression side', x: rs.map((r) => round(r * 1e3)), y: wall.map((w) => round(w.vmC / MPA)) },
          { name: 'Allowable η·f_y', x: [round(ri * 1e3), round(ro * 1e3)], y: [round(allow / MPA), round(allow / MPA)], dash: true },
        ], 'Stress components from the bore to the outer surface. The von Mises curves must stay below the dashed allowable; the bore usually governs under internal pressure.'),
        plot('line', 'von Mises stress around the circumference', 'Angle from the bending tension fibre', 'deg', 'von Mises stress', 'MPa', [
          { name: 'Outer surface', x: ang, y: vmAng.map((v) => round(v / MPA)) }, { name: 'Allowable', x: [0, 360], y: [round(allow / MPA), round(allow / MPA)], dash: true },
        ], 'Equivalent stress at the outer surface around the pipe. With compressive axial force the peak sits on the compression side of the bend (180°).'),
        plot('line', 'Equivalent stress vs bending moment', 'Bending moment', 'kN·m', 'Equivalent stress', 'MPa', [
          { name: 'Thin-wall von Mises', x: ms.map((m) => round(m / 1e3)), y: eqM.map((v) => round(v / MPA)) }, { name: 'Allowable', x: [0, round(ms[ms.length - 1] / 1e3)], y: [round(allow / MPA), round(allow / MPA)], dash: true },
          { name: 'Operating point', x: [round(M / 1e3)], y: [round(eq.sigmaE / MPA)] },
        ], 'How much bending the section can take on top of pressure and axial force. The crossing with the allowable is the allowable moment.'),
      ],
      tables: [tbl('Stress components', [col('c', 'Component'), col('v', 'Value', 'MPa'), col('src', 'Origin')], [
        { c: 'Hoop, thin wall', v: round(eq.sigmaH / MPA), src: 'pressure difference' }, { c: 'Hoop at bore (Lamé)', v: round(lam.sigmaHInner / MPA), src: 'pressure difference' },
        { c: 'Hoop at outer surface (Lamé)', v: round(lam.sigmaHOuter / MPA), src: 'pressure difference' }, { c: 'Axial', v: round(sA / MPA), src: 'restraint, temperature, end cap' },
        { c: 'Bending (outer fibre)', v: round(sB / MPA), src: `moment ${fmt(inp.M)} kN·m + curvature ${fmt(Mc / 1e3)} kN·m` }, { c: 'Torsional shear', v: round(tau / MPA), src: 'torque' },
        { c: 'Longitudinal, tension fibre', v: round(eq.sigmaLTension / MPA), src: 'axial + bending' }, { c: 'Longitudinal, compression fibre', v: round(eq.sigmaLCompression / MPA), src: 'axial − bending' },
      ], 'Each row is one stress component at the section with the load that causes it, so the driver of the equivalent stress can be read off.')],
      notes: [`${gov.name} governs at ${fmt(gov.uc)}. Temperature rise ${fmt(dT)} K; yield at temperature ${fmt(fy / MPA)} MPa.`,
        inp.restraint === 'restrained' ? `Fully restrained effective force ${fmt(S / 1e3)} kN (compression negative).` : 'Section treated as axially free or with the specified effective force.',
        'Stresses are on the corroded wall; usage and design factors are recalled defaults.'],
      advice: [
        adv('engineering', sev(gov.uc), `Governing stress utilisation ${fmt(gov.uc)} (${gov.name}); maximum von Mises ${fmt(vmMax / MPA)} MPa against ${fmt(allow / MPA)} MPa allowable.`,
          gov.uc > 1 ? `Reduce the bending moment below ${fmt(Mallow / 1e3)} kN·m (span correction, larger route radius), or thicken the wall to about ${fmt(tNeed * 1e3)} mm.` : `Up to ${fmt(Mallow / 1e3)} kN·m of bending is available at this section. Confirm η and the ASME factors for the location class and check local buckling under the same loads.`,
          'DNV-ST-F101 (2021) Sec. 5 ASD format; ASME B31.8 (2022) §A842.2 — factors recalled, confirm'),
        adv('economic', 'info', `The section uses ${fmt(Math.min(gov.uc, 9.99) * 100)} % of its allowable stress; one millimetre of wall is ${fmt(steelPerMetre(p0, 1e-3))} kg of steel per metre of pipe.`,
          gov.uc > 1 ? `Curing the overstress by wall thickness costs about ${fmt(steelPerMetre(p0, Math.max(tNeed - p0.t, 0)))} kg/m; relieving the moment by span rectification or route smoothing is usually cheaper for a local hot spot.` : 'Stress is not the cost driver here; do not add wall for stress without a governing limit state.',
          ECON_BASIS, steelCo2(steelPerMetre(p0, 1e-3)).replace(' t of line-pipe steel', ' t of steel per kilometre for each millimetre of wall')),
      ],
      publish: { 'pipeline.vonMises': vmMax },
    };
  },
  verify: [{
    // D = 323.9 mm, t = 14.3 mm (no corrosion allowance), p_i = 20 MPa, p_e = 0:
    // Lamé bore hoop = p_i (r_i² + r_o²)/(r_o² − r_i²) = 20·(147.65² + 161.95²)/(161.95² − 147.65²) = 20·48028.3/4427.3 = 216.96 MPa.
    // Thin wall: p (D − t)/(2t) = 20·309.6/28.6 = 216.50 MPa.
    name: 'Lamé and thin-wall hoop stress under internal pressure',
    source: 'Timoshenko, Strength of Materials Part II, thick cylinders; hand calculation',
    inputs: { D: 323.9, t: 14.3, tCorr: 0, depth: 0, pi: 200, restraint: 'unrestrained', M: 0, Rbend: 0, Tq: 0 },
    expect: [{ kpi: 'hoopLame', value: 216.96, tol: 0.002 }, { kpi: 'hoopThin', value: 216.5, tol: 0.002 }],
  }],
};
