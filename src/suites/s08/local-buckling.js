// Burst, local buckling and combined loading: DNV-ST-F101 load-controlled and displacement-controlled local buckling
// criteria under bending, effective axial force and internal or external overpressure, with the moment–pressure
// interaction envelope and bending ovalisation.
import { characteristicStrength, fabricationFactor } from '../../engine/pipeline/pipeprops.js';
import { dnvResistanceFactors, plasticCapacities, localBucklingLoadControlled, localBucklingDisplacementControlled, pressureContainment,
  ovalisationBendingPressure, collapsePressure, burstPressure } from '../../engine/pipeline/strength.js';
import { pipeInputs, pipeFrom, depthInput, safetyClassInput, num, sel, kpi, chk, chkUc, adv, plot, tbl, col, run, need, fmt, round, linspace, sev, worst,
  steelPerMetre, steelCo2, CONFIRM, LIMIT_F101, LIMIT_NO_FE, STD_F101, ECON_BASIS, G, BAR, MPA } from './_common.js';

const G_L = 'Load effects';
const G_C = 'Acceptance criteria';

export default {
  id: 's08-local-buckling',
  suite: 8,
  title: 'Burst, local buckling and combined loading',
  summary: 'Checks the pipe wall against local buckling under combined bending, axial force and pressure, and draws the moment–pressure capacity envelope.',
  covers: [11],
  fidelity: 'analytical',
  standards: [`${STD_F101} §5.4.6 (local buckling — combined loading), Eq. 5.19–5.32`, `${STD_F101} §5.4.2 (pressure containment)`, `${STD_F101} §5.4.12 (ovalisation)`],
  models: ['Combined loading collapse criteria (load controlled)', 'Strain-based local buckling (displacement controlled)', 'Plastic moment and axial capacities', 'Burst pressure model', 'External pressure collapse with ovality', 'Bending ovalisation (Murphey & Langner form)'],
  limits: `${LIMIT_F101} Valid for 15 ≤ D/t ≤ 45 and |S_Sd|/S_p < 0.4. The girth-weld factor α_gw is not built in (1.0 up to D/t = 20, lower above — read it from the standard). The ovalisation coefficients are recalled. ${LIMIT_NO_FE}`,
  publishes: ['pipeline.allowableMoment'],
  searchTerms: 'pipeline local buckling combined loading DNV-ST-F101 load controlled displacement controlled',
  inputs: [
    ...pipeInputs({ omit: ['tCoat', 'rhoCoat', 'tConc', 'concrete', 'rhoCont', 'wsGiven'] }),
    num('tCorr', 'Corrosion allowance', 'mm', 3, 0, 15, 'Pipe', 'geometry', 'Deducted to obtain t2 in operation.'),
    sel('process', 'Pipe manufacturing process', 'seamless', ['seamless', 'UO', 'TRB', 'ERW', 'HFW', 'UOE'], 'Pipe', 'material', 'Sets the fabrication factor α_fab used in the collapse term.'),
    num('f0', 'Ovality (Dmax − Dmin)/D', '%', 1, 0.5, 5, 'Pipe', 'geometry', 'Initial out-of-roundness; not below 0.5 %.'),
    depthInput(120, G_L),
    num('pi', 'Internal pressure', 'bar', 200, 0, 1500, G_L, 'bc', 'Local incidental (or design) pressure acting with the load effects.', { uncertain: { dist: 'normal', cov: 0.03 } }),
    num('pMin', 'Minimum internal pressure', 'bar', 0, 0, 1500, G_L, 'bc', 'Lowest pressure that can act with the bending; used for external overpressure and for the strain criterion.'),
    num('Tdes', 'Design temperature', '°C', 80, -20, 200, G_L, 'bc', 'For yield and tensile de-rating.'),
    num('MF', 'Functional bending moment', 'kN·m', 120, 0, 50000, G_L, 'load', 'Characteristic moment from weight, spans, supports and thermal bending.', { uncertain: { dist: 'normal', cov: 0.1 } }),
    num('ME', 'Environmental bending moment', 'kN·m', 40, 0, 50000, G_L, 'load', 'Characteristic moment from waves, current or vessel motion.', { uncertain: { dist: 'gumbel', cov: 0.2 } }),
    num('SF', 'Functional effective axial force (tension +)', 'kN', -800, -50000, 50000, G_L, 'load', 'Characteristic effective force; compression is negative.', { bind: 'pipeline.effectiveForce' }),
    num('SE', 'Environmental effective axial force', 'kN', 0, -50000, 50000, G_L, 'load', 'Characteristic environmental axial force.'),
    num('epsF', 'Functional compressive bending strain', '%', 0.15, 0, 5, G_L, 'load', 'Characteristic longitudinal compressive strain for displacement-controlled conditions (reeling, laying over a stinger, buckled sections).'),
    num('epsE', 'Environmental compressive strain', '%', 0.03, 0, 5, G_L, 'load', 'Characteristic environmental strain.'),
    safetyClassInput(G_C),
    num('gammaF', 'Functional load effect factor γ_F', '-', 1.1, 1, 1.3, G_C, 'criteria', `Load combination b: 1.1 (combination a: 1.2). ${CONFIRM}`),
    num('gammaE', 'Environmental load effect factor γ_E', '-', 1.3, 0.7, 1.5, G_C, 'criteria', `Load combination b: 1.3 (combination a: 0.7). ${CONFIRM}`),
    num('gammaC', 'Condition load effect factor γ_c', '-', 1, 0.8, 1.2, G_C, 'criteria', `1.07 pipe resting on uneven seabed, 0.82 continuously stiff supported, 0.93 system pressure test, 1.0 otherwise. ${CONFIRM}`),
    num('alphaU', 'Material strength factor α_U', '-', 0.96, 0.9, 1, G_C, 'criteria', '0.96 normally; 1.00 with supplementary requirement U.'),
    num('alphaH', 'Strain-hardening ratio α_h (yield/tensile)', '-', 0.93, 0.8, 0.98, G_C, 'criteria', 'Maximum allowed yield-to-tensile ratio of the line pipe (0.93 unless otherwise specified).'),
    num('alphaGw', 'Girth-weld factor α_gw', '-', 1, 0.5, 1, G_C, 'criteria', '1.0 for D/t ≤ 20, reducing above — read it from the girth-weld factor figure of the governing DNV-ST-F101 edition.'),
    num('gammaEps', 'Strain resistance factor γ_ε (0 = table)', '-', 0, 0, 4, G_C, 'criteria', `Zero applies the recalled table value 2.0 / 2.5 / 3.3 for low / medium / high. ${CONFIRM}`),
    num('ovalLimit', 'Allowable ovalisation', '%', 3, 1, 6, G_C, 'criteria', 'DNV-ST-F101 limits flattening to 3 % in general.'),
  ],
  compute(inp, ctx) {
    const p = pipeFrom(inp, ctx);
    const t2 = p.t - inp.tCorr / 1e3;
    need(t2 > 0.002, 'The corrosion allowance consumes the whole wall.');
    const sc = inp.safetyClass, tab = dnvResistanceFactors(sc);
    const str = run(() => characteristicStrength({ smys: p.smys, smts: p.smts, tempC: Math.max(inp.Tdes, 20), alphaU: inp.alphaU }), 'Characteristic strength');
    const pe = p.rhoW * G * inp.depth, pi = inp.pi * BAR, pMin = inp.pMin * BAR;
    const aFab = fabricationFactor(inp.process), f0 = inp.f0 / 100;
    const MSd = (inp.MF * inp.gammaF * inp.gammaC + inp.ME * inp.gammaE) * 1e3;
    const SSd = (inp.SF * inp.gammaF * inp.gammaC + inp.SE * inp.gammaE) * 1e3;
    const epsSd = (inp.epsF * inp.gammaF * inp.gammaC + inp.epsE * inp.gammaE) / 100;
    const base = { D: p.D, t2, fy: str.fy, fu: str.fu, safetyClass: sc, alphaFab: aFab, f0, E: p.E, nu: p.nu };
    const lcc = run(() => localBucklingLoadControlled({ ...base, MSd, SSd, pi, pe, pMin }), 'Load-controlled local buckling');
    // the empty / depressurised state is a separate load case when the line can lose pressure
    const lccEmpty = run(() => localBucklingLoadControlled({ ...base, MSd, SSd, pi: pMin, pe, pMin }), 'Load-controlled local buckling (minimum pressure)');
    const dccArgs = { ...base, epsSd, pMin, pe, alphaH: inp.alphaH, alphaGw: inp.alphaGw };
    if (inp.gammaEps > 0) dccArgs.gammaEps = inp.gammaEps;
    const dcc = epsSd > 0 ? run(() => localBucklingDisplacementControlled(dccArgs), 'Displacement-controlled local buckling') : null;
    const cap = plasticCapacities({ D: p.D, t: t2, fy: str.fy });
    const burst = run(() => pressureContainment({ D: p.D, t: p.t, tCorr: inp.tCorr / 1e3, fy: str.fy, fu: str.fu, pLocalInc: pi, pe, safetyClass: sc, alphaU: inp.alphaU }), 'Pressure containment');
    const pc = collapsePressure({ D: p.D, t: t2, fy: str.fy, E: p.E, nu: p.nu, alphaFab: aFab, f0 });
    const oval = run(() => ovalisationBendingPressure({ D: p.D, t: t2, strain: Math.max(epsSd, MSd * p.D / (2 * p.EI)), f0, pe: Math.max(pe - pMin, 0), pc: pc.pc, limit: inp.ovalLimit / 100 }), 'Ovalisation');
    const checks = [
      chkUc(`Local buckling, load controlled (${lcc.mode} overpressure)`, lcc.uc, `${STD_F101} Eq. 5.19 / 5.28`,
        'Guards against the compressed side of the wall wrinkling under design moment, axial force and pressure. Exceeding 1.0 means a local buckle forms and the section loses its moment capacity.'),
      chkUc('Local buckling, load controlled at minimum pressure', lccEmpty.uc, `${STD_F101} Eq. 5.28`,
        'Same criterion with the line depressurised, when external overpressure adds to the bending. Exceeding it means the pipe can buckle and collapse during a shutdown or before commissioning.'),
      chk('Pressure containment (burst)', burst.demand, burst.capacity, 'Pa', `${STD_F101} Eq. 5.6`,
        'Guards against bursting under the internal overpressure alone; combined loading cannot be accepted if containment is not.'),
      chk('Ovalisation under bending and external pressure', oval.ovality, oval.limit, '-', `${STD_F101} Eq. 5.34 (recalled coefficients)`,
        'Guards against the cross-section flattening so far that pigs cannot pass and collapse resistance is lost.'),
    ];
    if (dcc) checks.splice(2, 0, chkUc(`Local buckling, displacement controlled (${dcc.mode} overpressure)`, dcc.uc, `${STD_F101} Eq. 5.30 / 5.31`,
      'Guards against local buckling where the pipe is forced to a given curvature (reel, stinger, sleeper, lateral buckle): the design compressive strain must stay below the factored critical strain.'));
    const gov = worst(checks);
    // capacity envelope: allowable design moment against internal pressure, with the axial force held
    const pbD = burst.pbDesign ?? burstPressure({ D: p.D, t: t2, fy: str.fy, fu: str.fu }).pb / (tab.gammaM * tab.gammaScPc);
    const ps = linspace(0, Math.max(pe + 0.98 * pbD, pi * 1.05), 30);
    const env = ps.map((q) => { try { return localBucklingLoadControlled({ ...base, MSd, SSd, pi: q, pe, pMin: Math.min(q, pMin) }).allowableMoment; } catch { return 0; } });
    const ss = linspace(-0.39 * cap.Sp, 0.39 * cap.Sp, 21);
    const envS = ss.map((s) => { try { return localBucklingLoadControlled({ ...base, MSd, SSd: s, pi, pe, pMin }).allowableMoment; } catch { return 0; } });
    const Mallow = Math.min(lcc.allowableMoment, lccEmpty.allowableMoment);
    const ucL = Math.max(lcc.uc, lccEmpty.uc);
    // wall for UC = 1 by bisection on the nominal thickness (load-controlled criteria only)
    let tReq = p.t;
    const ucAt = (tn) => { const tt = tn - inp.tCorr / 1e3; if (tt < 0.003) return 99; try { return Math.max(localBucklingLoadControlled({ ...base, t2: tt, MSd, SSd, pi, pe, pMin }).uc, localBucklingLoadControlled({ ...base, t2: tt, MSd, SSd, pi: pMin, pe, pMin }).uc); } catch { return 99; } };
    { let lo = inp.tCorr / 1e3 + 0.003, hi = p.D / 4; if (ucAt(hi) < 1) { for (let i = 0; i < 40; i++) { const m = 0.5 * (lo + hi); if (ucAt(m) > 1) lo = m; else hi = m; } tReq = hi; } }
    return {
      kpis: [
        kpi('ucLCC', 'Load-controlled utilisation', ucL, '-', 'Left-hand side of the combined-loading interaction equation for the worse of the pressurised and depressurised states. It must not exceed 1.0.', 1, 'max'),
        kpi('ucDCC', 'Displacement-controlled utilisation', dcc ? dcc.uc : 0, '-', 'Design compressive strain over the factored critical strain (with the external pressure term where relevant). Zero when no strain is entered.', 1, 'max'),
        kpi('Mallow', 'Allowable design moment', Mallow / 1e3, 'kN·m', 'Design moment that brings the load-controlled interaction to 1.0 with the axial force and pressure unchanged.'),
        kpi('MSd', 'Design bending moment', MSd / 1e3, 'kN·m', 'Characteristic moments multiplied by the load effect factors.'),
        kpi('SSd', 'Design effective axial force', SSd / 1e3, 'kN', 'Characteristic effective forces multiplied by the load effect factors; compression negative.'),
        kpi('Mp', 'Plastic moment capacity', cap.Mp / 1e3, 'kN·m', 'f_y (D − t2)² t2: the moment at which the whole corroded section has yielded in bending.'),
        kpi('Sp', 'Plastic axial capacity', cap.Sp / 1e3, 'kN', 'f_y π (D − t2) t2: the axial force that yields the whole corroded section.'),
        kpi('alphaC', 'Flow stress parameter α_c', lcc.alphaC, '-', 'Accounts for strain hardening; rises with the tensile-to-yield ratio and falls with D/t.'),
        kpi('epsC', 'Critical compressive strain', dcc ? dcc.epsC * 100 : 0, '%', 'Characteristic bending strain at which the wall buckles locally, before the strain resistance factor.'),
        kpi('epsSd', 'Design compressive strain', epsSd * 100, '%', 'Factored compressive strain demand.'),
        kpi('ucBurst', 'Pressure containment utilisation', burst.uc, '-', 'Internal overpressure over the factored burst resistance.', 1, 'max'),
        kpi('ovality', 'Ovalisation under bending', oval.ovality * 100, '%', 'Predicted (Dmax − Dmin)/D under the design strain and external overpressure.', inp.ovalLimit, 'max'),
        kpi('tReq', 'Wall thickness for load-controlled UC = 1', tReq * 1e3, 'mm', 'Nominal wall at which the load-controlled criterion is just met for these load effects.'),
        kpi('DtRatio', 'D/t2', p.D / t2, '-', 'Diameter-to-thickness ratio of the corroded wall; the criteria are calibrated for 15 to 45.', 45, 'max'),
      ],
      checks,
      plots: [
        plot('line', 'Moment capacity vs internal pressure', 'Internal pressure', 'bar', 'Allowable design moment', 'kN·m', [
          { name: 'Capacity envelope', x: ps.map((q) => round(q / BAR)), y: env.map((m) => round(m / 1e3)) },
          { name: 'Design point (operating)', x: [round(pi / BAR)], y: [round(MSd / 1e3)] }, { name: 'Design point (minimum pressure)', x: [round(pMin / BAR)], y: [round(MSd / 1e3)] },
        ], 'Interaction diagram: the allowable design moment at each internal pressure with the axial force held. The design points must lie under the envelope; the envelope closes to zero at the burst limit and dips at low pressure where external overpressure acts.'),
        plot('line', 'Moment capacity vs effective axial force', 'Design effective axial force', 'kN', 'Allowable design moment', 'kN·m', [
          { name: 'Capacity envelope', x: ss.map((s) => round(s / 1e3)), y: envS.map((m) => round(m / 1e3)) }, { name: 'Design point', x: [round(SSd / 1e3)], y: [round(MSd / 1e3)] },
        ], 'Allowable moment against axial force at the operating pressure. The curve is symmetric in tension and compression in the code equation; the point must lie below it.'),
      ],
      tables: [tbl('Interaction terms', [col('term', 'Term'), col('op', 'Operating pressure', '-'), col('min', 'Minimum pressure', '-')], [
        { term: 'Moment term γ M_Sd/(α_c M_p)', op: round(lcc.momentTerm), min: round(lccEmpty.momentTerm) }, { term: 'Axial term [γ S_Sd/(α_c S_p)]²', op: round(lcc.axialTerm), min: round(lccEmpty.axialTerm) },
        { term: 'Pressure term', op: round(lcc.pressureTerm), min: round(lccEmpty.pressureTerm) }, { term: 'Utilisation', op: round(lcc.uc), min: round(lccEmpty.uc) },
      ], 'The three terms of the load-controlled interaction equation for the two pressure states, showing which load effect uses the capacity.')],
      notes: [`${gov.name} governs at ${fmt(gov.uc)}. f_y = ${fmt(str.fy / MPA)} MPa, f_u = ${fmt(str.fu / MPA)} MPa at ${fmt(inp.Tdes)} °C; γ_m·γ_SC = ${fmt(tab.gammaM * tab.gammaSc)}.`,
        ...[...lcc.warnings, ...lccEmpty.warnings].filter((w, i, a) => a.indexOf(w) === i).map((w) => `Validity: ${w}.`),
        'Load-controlled criteria apply where the load is independent of the deformation; use the strain criterion only where the curvature is imposed.'],
      advice: [
        adv('engineering', sev(gov.uc), `Governing combined-loading utilisation ${fmt(gov.uc)} (${gov.name}); allowable design moment ${fmt(Mallow / 1e3)} kN·m against ${fmt(MSd / 1e3)} kN·m applied.`,
          gov.uc > 1 ? `Reduce the design moment below ${fmt(Mallow / 1e3)} kN·m (shorter spans, more supports) or increase the wall to about ${fmt(tReq * 1e3)} mm.` : `${fmt((Mallow - MSd) / 1e3)} kN·m of design moment remains. Confirm γ_F, γ_E, γ_c and α_gw against the governing edition before relying on it.`,
          `${STD_F101} §5.4.6 with Tables 5-1, 5-2, 5-8 (recalled — confirm)`),
        adv('economic', 'info', `The load-controlled criterion needs a ${fmt(tReq * 1e3)} mm wall against ${fmt(inp.t)} mm specified (${fmt(steelPerMetre(p, tReq - p.t))} kg/m difference).`,
          tReq > p.t ? `Thickening only the affected joints, or supporting the span, avoids adding ${fmt(steelPerMetre(p, tReq - p.t))} t per kilometre along the whole line.` : `Combined loading leaves ${fmt((p.t - tReq) * 1e3)} mm of wall unused at this section; check whether pressure containment or installation sets the wall before optimising.`,
          ECON_BASIS, steelCo2(Math.abs(steelPerMetre(p, tReq - p.t))).replace(' t of line-pipe steel', ' t of line-pipe steel per kilometre')),
      ],
      publish: { 'pipeline.allowableMoment': Mallow },
    };
  },
  verify: [{
    // Pure bending, no pressure, no axial force. D = 323.9, t2 = 14.3 mm, X65 at 20 °C with α_U = 0.96: f_y = 432, f_u = 513.6 MPa.
    // M_p = f_y (D − t)² t = 432e6·0.3096²·0.0143 = 592.1 kN·m; S_p = f_y π (D − t) t = 6008 kN.
    // β = (60 − 22.65)/90 = 0.4150; α_c = (1 − β) + β f_u/f_y = 1 + 0.4150·0.18889 = 1.0784.
    // UC = [γm γSC M/(α_c M_p)]² = [1.15·1.14·300/(1.0784·592.1)]² = 0.6159² = 0.3794 (medium safety class).
    name: 'Pure bending: plastic moment and load-controlled interaction',
    source: 'DNV-ST-F101 (2021) Eq. 5.19–5.24; hand calculation',
    inputs: { D: 323.9, t: 14.3, tCorr: 0, depth: 0, pi: 0, pMin: 0, Tdes: 20, MF: 300, ME: 0, SF: 0, SE: 0, epsF: 0, epsE: 0, gammaF: 1, gammaC: 1, safetyClass: 'medium' },
    expect: [{ kpi: 'Mp', value: 592.1, tol: 0.002 }, { kpi: 'Sp', value: 6008, tol: 0.002 }, { kpi: 'ucLCC', value: 0.3794, tol: 0.005 }],
  }],
};
