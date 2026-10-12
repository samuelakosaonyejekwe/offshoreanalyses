// Wall-thickness design for internal pressure (pressure containment / burst) and external pressure (system collapse,
// propagating buckling, buckle arrestors) to DNV-ST-F101, with the ASME B31.8 hoop criterion alongside.
import { characteristicStrength, fabricationFactor } from '../../engine/pipeline/pipeprops.js';
import { dnvResistanceFactors, characteristicWallThickness, burstPressure, localIncidentalPressure, millTestPressure, pressureContainment,
  requiredWallThicknessBurst, collapsePressure, collapseCheck, requiredWallThicknessCollapse, propagationPressure, buckleArrestorCapacity,
  asmeStressCheck, hydrotestPressure } from '../../engine/pipeline/strength.js';
import { pipeInputs, pipeFrom, lengthInput, depthInput, safetyClassInput, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, linspace, sev, worst,
  steelPerMetre, steelCo2, CONFIRM, LIMIT_F101, STD_F101, ECON_BASIS, G, BAR, MPA } from './_common.js';

const G_P = 'Pressure and temperature';
const G_T = 'Thickness allowances';
const G_C = 'Acceptance criteria';

export default {
  id: 's08-wall-thickness',
  suite: 8,
  title: 'Wall thickness for internal and external pressure',
  summary: 'Checks burst, system collapse and propagating buckling and returns the wall thickness each limit state requires.',
  covers: [10],
  fidelity: 'analytical',
  standards: [`${STD_F101} §5.4.2 (pressure containment), §5.4.4 (collapse), §5.4.5 (propagation buckling)`, 'ASME B31.8 (2022) §841.1.1, §A842.2', 'ASME B31.4 (2022) §403.2'],
  models: ['Pressure containment (burst) resistance', 'External pressure collapse with ovality (cubic interaction of elastic and plastic collapse)', 'Buckle propagation pressure', 'Integral buckle arrestor crossover pressure', 'Barlow hoop stress (ASME)'],
  limits: `${LIMIT_F101} Collapse and propagation expressions are calibrated for 15 ≤ D/t ≤ 45. The buckle arrestor coefficients are recalled. Local buckling under combined bending, axial force and pressure is in the combined-loading module; finite-element collapse analysis with measured imperfections is outside this module.`,
  publishes: ['pipeline.D', 'pipeline.t', 'pipeline.length', 'quantities.steel_t', 'quantities.pipe_km'],
  searchTerms: 'pipeline wall thickness burst collapse propagation buckling DNV-ST-F101',
  inputs: [
    ...pipeInputs({ omit: ['tCoat', 'rhoCoat', 'tConc', 'concrete', 'wsGiven'] }),
    lengthInput(10, 'Pipe'),
    depthInput(120, G_P),
    num('pd', 'Design pressure at the reference elevation', 'bar', 200, 1, 1500, G_P, 'bc', 'Internal design pressure (gauge) at the reference point, usually the topside or wellhead.', { design: true, uncertain: { dist: 'normal', cov: 0.03 } }),
    num('hRef', 'Elevation of the pressure reference', 'm', 20, -3500, 100, G_P, 'bc', 'Positive above sea level. The contents column between this point and the pipe adds to the local pressure.'),
    num('pMin', 'Minimum internal pressure', 'bar', 0, 0, 1500, G_P, 'bc', 'Lowest internal pressure that can be guaranteed (zero during installation or depressurisation).'),
    num('Tdes', 'Design temperature', '°C', 80, -20, 200, G_P, 'bc', 'Maximum design temperature, used for yield de-rating.'),
    num('tCorr', 'Corrosion allowance', 'mm', 3, 0, 15, G_T, 'geometry', 'Internal plus external corrosion allowance.', { design: true }),
    num('tFab', 'Fabrication minus tolerance', 'mm', 1, 0, 10, G_T, 'geometry', 'Negative wall-thickness tolerance of the line pipe (welded pipe ≈ 0.5–1 mm; seamless often 12.5 % of t).'),
    sel('process', 'Pipe manufacturing process', 'seamless', ['seamless', 'UO', 'TRB', 'ERW', 'HFW', 'UOE'], G_T, 'material', 'Sets the fabrication factor α_fab for collapse and propagation (1.00 seamless, 0.93 UO/TRB/ERW, 0.85 UOE).'),
    num('f0', 'Ovality (Dmax − Dmin)/D', '%', 1, 0.5, 5, G_T, 'geometry', 'Out-of-roundness from fabrication and installation; not below 0.5 %.', { uncertain: { dist: 'lognormal', cov: 0.3 } }),
    num('tBA', 'Buckle arrestor wall thickness (0 = none)', 'mm', 0, 0, 120, G_T, 'geometry', 'Integral arrestor thickness, checked only when propagation is not prevented by the pipe itself.'),
    num('LBA', 'Buckle arrestor length', 'm', 2, 0.2, 12, G_T, 'geometry', 'Length of the integral arrestor.'),
    safetyClassInput(G_C),
    sel('factorSource', 'Partial factors', 'table', [['table', 'Recalled DNV-ST-F101 tables by safety class'], ['manual', 'Typed below']], G_C, 'criteria', CONFIRM),
    num('gammaM', 'Material resistance factor γ_m', '-', 1.15, 1, 1.5, G_C, 'criteria', `ULS value 1.15. ${CONFIRM}`, { showIf: { key: 'factorSource', in: ['manual'] } }),
    num('gammaScPc', 'Safety class factor, pressure containment', '-', 1.138, 1, 1.6, G_C, 'criteria', `1.046 / 1.138 / 1.308 for low / medium / high. ${CONFIRM}`, { showIf: { key: 'factorSource', in: ['manual'] } }),
    num('gammaSc', 'Safety class factor, other limit states', '-', 1.14, 1, 1.6, G_C, 'criteria', `1.04 / 1.14 / 1.26 for low / medium / high. ${CONFIRM}`, { showIf: { key: 'factorSource', in: ['manual'] } }),
    num('gammaInc', 'Incidental / design pressure ratio', '-', 1.1, 1, 1.2, G_C, 'criteria', '1.10 for a typical pressure control system; 1.05 where the safety system justifies it.'),
    num('alphaU', 'Material strength factor α_U', '-', 0.96, 0.9, 1, G_C, 'criteria', '0.96 normally; 1.00 with supplementary requirement U.'),
    num('alphaSpt', 'System test factor α_spt', '-', 1.05, 1, 1.3, G_C, 'criteria', `Ratio of local test pressure to local incidental pressure (1.03 low, 1.05 medium/high). ${CONFIRM}`),
    num('F1', 'ASME hoop design factor', '-', 0.72, 0.4, 0.8, G_C, 'criteria', 'ASME B31.8 Table A842.2.2-1: 0.72 offshore pipeline, 0.50 riser. Confirm against the edition in force.'),
  ],
  compute(inp, ctx) {
    const p = pipeFrom(inp, ctx);
    const sc = inp.safetyClass;
    const tab = dnvResistanceFactors(sc);
    const man = inp.factorSource === 'manual';
    const gM = man ? inp.gammaM : tab.gammaM, gPc = man ? inp.gammaScPc : tab.gammaScPc, gSc = man ? inp.gammaSc : tab.gammaSc;
    const str = run(() => characteristicStrength({ smys: p.smys, smts: p.smts, tempC: Math.max(inp.Tdes, 20), alphaU: inp.alphaU }), 'Characteristic strength');
    const tCorr = inp.tCorr / 1e3, tFab = inp.tFab / 1e3, f0 = inp.f0 / 100;
    need(p.t - tCorr - tFab > 0.002, 'Nothing is left of the wall after the corrosion allowance and fabrication tolerance: increase the wall thickness.');
    const tw = characteristicWallThickness({ t: p.t, tFab, tCorr });
    const twInst = characteristicWallThickness({ t: p.t, tFab, tCorr, condition: 'test' });
    const pe = p.rhoW * G * inp.depth;
    const pli = localIncidentalPressure({ pRef: inp.pd * BAR, gammaInc: inp.gammaInc, rhoCont: inp.rhoCont, hRef: inp.hRef, hLocal: -inp.depth });
    const plt = inp.alphaSpt * pli;
    const pMill = millTestPressure({ D: p.D, tMin: p.t - tFab, smys: p.smys, smts: p.smts });
    const pcn = run(() => pressureContainment({ D: p.D, t: p.t, tFab, tCorr, fy: str.fy, fu: str.fu, pLocalInc: pli, pe, safetyClass: sc, pMill, alphaU: inp.alphaU, gammaM: gM, gammaSc: gPc }), 'Pressure containment');
    // the system-test term of Eq. 5.6 is an identity when the test pressure is the minimum α_spt·p_li, so the test is checked on its own (Eq. 5.7)
    const tst = run(() => pressureContainment({ D: p.D, t: p.t, tFab, tCorr, fy: str.fy, fu: str.fu, pLocalInc: pli, pe, safetyClass: sc, pLocalTest: plt, pMill, alphaU: inp.alphaU, gammaM: gM, gammaSc: gPc, alphaSpt: inp.alphaSpt }), 'System pressure test').testCheck;
    const aFab = fabricationFactor(inp.process);
    const pMin = inp.pMin * BAR;
    const colArgs = { D: p.D, fy: str.fy, pe, pMin, safetyClass: sc, E: p.E, nu: p.nu, alphaFab: aFab, f0, gammaM: gM, gammaSc: gSc };
    const col1 = run(() => collapseCheck({ ...colArgs, t1: tw.t1 }), 'Collapse');
    const prop = run(() => propagationPressure({ D: p.D, t2: tw.t2, fy: str.fy, alphaFab: aFab, pe, pMin, safetyClass: sc, gammaM: gM, gammaSc: gSc }), 'Propagation buckling');
    const tBurst = requiredWallThicknessBurst({ D: p.D, fy: str.fy, fu: str.fu, pLocalInc: pli, pe, safetyClass: sc, tCorr, tFab, gammaM: gM, gammaSc: gPc });
    const tColl = pe > pMin ? run(() => requiredWallThicknessCollapse(colArgs), 'Required collapse thickness') + tCorr + tFab : 0;
    // propagation: t2 from p_e − p_min = 35 f_y α_fab (t2/D)^2.5/(γm γSC)
    const tProp = pe > pMin ? p.D * ((pe - pMin) * gM * gSc / (35 * str.fy * aFab)) ** 0.4 + tCorr : 0;
    const asme = asmeStressCheck({ D: p.D, t: p.t - tCorr - tFab, smys: p.smys, pi: pli / inp.gammaInc, pe, F1: inp.F1 });
    const hydro = hydrotestPressure({ pDesign: pli / inp.gammaInc, D: p.D, t: p.t - tFab, smys: p.smys, pe });
    const checks = [
      chk('Pressure containment (burst), operation', pcn.demand, pcn.capacity, 'Pa', `${STD_F101} Eq. 5.6 (governing: ${pcn.governing})`,
        'Guards against bursting under the local incidental pressure: the net internal overpressure must stay below the factored burst resistance of the thinnest corroded wall. Exceeding it risks rupture.'),
      chk('System collapse under external pressure', col1.demand, col1.capacity, 'Pa', `${STD_F101} Eq. 5.10`,
        'Guards against the cross-section flattening under external overpressure when the pipe is empty or depressurised. Exceeding it collapses the pipe locally.'),
      chk('Propagating buckling (no arrestors)', prop.demand, prop.capacity, 'Pa', `${STD_F101} Eq. 5.16`,
        'Guards against a local buckle, once started, running along the line under external pressure. If it exceeds 1.0 buckle arrestors are required; it is not a failure of the pipe by itself.'),
      chk('Hoop stress (ASME B31.8 offshore)', asme.SH, asme.allowHoop, 'Pa', 'ASME B31.8 §A842.2.2',
        'Guards against yielding of the wall in the hoop direction at design pressure with the ASME design factor on SMYS.'),
    ];
    if (tst) checks.push(chk('Pressure containment, system pressure test', tst.demand, tst.capacity, 'Pa', `${STD_F101} Eq. 5.7`,
      'Guards against bursting or yielding during the system pressure test, when the wall has no corrosion loss but the pressure is highest.'));
    let arrestor = null;
    if (inp.tBA > 0) {
      arrestor = run(() => buckleArrestorCapacity({ D: p.D, t2: tw.t2, tBA: inp.tBA / 1e3, LBA: inp.LBA, fy: str.fy, alphaFab: aFab, pe, safetyClass: sc, gammaM: gM, gammaSc: gSc }), 'Buckle arrestor');
      checks.push(chk('Buckle arrestor crossover', pe, arrestor.capacity, 'Pa', `${STD_F101} Eq. 5.17–5.18 (recalled coefficients)`,
        'Guards against a propagating buckle crossing the arrestor; if exceeded the buckle runs through and the arrestor is ineffective.'));
    }
    const strength = checks.filter((c) => !c.name.startsWith('Propagating'));
    const gov = worst(strength);
    const tReq = Math.max(tBurst.tNominal, tColl);
    const dt = tReq - p.t, L = inp.Lkm * 1e3;
    const dSteel = steelPerMetre(p, dt) * L / 1e3;
    const ts = linspace(Math.max(0.4 * p.t, tCorr + tFab + 0.003), 2 * p.t, 25);
    const cap = ts.map((t) => {
      const w = characteristicWallThickness({ t, tFab, tCorr });
      return { b: burstPressure({ D: p.D, t: w.t1, fy: str.fy, fu: str.fu }).pb / (gM * gPc), c: collapsePressure({ D: p.D, t: w.t1, fy: str.fy, E: p.E, nu: p.nu, alphaFab: aFab, f0 }).pc / (gM * gSc), pr: 35 * str.fy * aFab * (w.t2 / p.D) ** 2.5 / (gM * gSc) };
    });
    const pcChar = collapsePressure({ D: p.D, t: tw.t1, fy: str.fy, E: p.E, nu: p.nu, alphaFab: aFab, f0 });
    const depthMax = col1.capacity / (p.rhoW * G) + pMin / (p.rhoW * G);
    const tx = ts.map((t) => round(t * 1e3));
    return {
      kpis: [
        kpi('ucBurst', 'Pressure containment utilisation', pcn.uc, '-', 'Net internal overpressure over the factored burst resistance (and the test-pressure limits). Above 1.0 the wall is too thin for the design pressure.', 1, 'max'),
        kpi('ucCollapse', 'Collapse utilisation', col1.uc, '-', 'External overpressure over the factored collapse resistance of the corroded wall. Above 1.0 the empty pipe collapses at this depth.', 1, 'max'),
        kpi('ucProp', 'Propagation utilisation', prop.uc, '-', 'External overpressure over the factored propagation pressure. Above 1.0 a buckle would run and arrestors are needed.'),
        kpi('tReq', 'Required nominal wall thickness', tReq * 1e3, 'mm', 'Largest of the burst and collapse requirements including corrosion allowance and fabrication tolerance. Round up to the next standard wall.'),
        kpi('tReqBurst', 'Wall thickness for burst', tBurst.tNominal * 1e3, 'mm', 'Nominal thickness at which the pressure containment criterion is just met.'),
        kpi('tReqCollapse', 'Wall thickness for collapse', tColl * 1e3, 'mm', 'Nominal thickness at which the system collapse criterion is just met at this depth (zero when there is no external overpressure).'),
        kpi('tReqProp', 'Wall thickness to stop propagation', tProp * 1e3, 'mm', 'Nominal thickness above which a buckle cannot propagate, so no arrestors are needed.'),
        kpi('pb', 'Burst resistance p_b(t1)', pcn.pb / BAR, 'bar', 'Characteristic pressure containment resistance of the minimum corroded wall, before partial factors.'),
        kpi('pc', 'Collapse pressure p_c(t1)', pcChar.pc / BAR, 'bar', 'Characteristic collapse pressure of the minimum corroded wall with the stated ovality.'),
        kpi('pel', 'Elastic collapse pressure', pcChar.pel / BAR, 'bar', 'Collapse pressure of a perfectly round elastic tube; governs thin walls.'),
        kpi('ppr', 'Propagation pressure', prop.ppr / BAR, 'bar', 'External pressure needed to keep an existing buckle running; always well below the collapse pressure.'),
        kpi('pli', 'Local incidental pressure', pli / BAR, 'bar', 'Design pressure times the incidental factor plus the contents head down to the pipe.'),
        kpi('pe', 'External pressure', pe / BAR, 'bar', 'Hydrostatic sea-water pressure at the pipe.'),
        kpi('depthMax', 'Limiting water depth for collapse', depthMax, 'm', 'Depth at which the factored collapse resistance is reached with the minimum internal pressure.'),
        kpi('hoopFrac', 'Hoop stress / SMYS at design pressure', asme.SH / p.smys, '-', 'Barlow hoop stress of the corroded wall as a fraction of SMYS; the ASME offshore limit is the design factor F1.', inp.F1, 'max'),
        kpi('pTest', 'Hydrotest pressure (1.25 × design)', hydro.pTest / BAR, 'bar', 'ASME strength test pressure at the pipe elevation and the reference for the test hoop stress.'),
      ],
      checks,
      plots: [
        plot('line', 'Design resistance vs wall thickness', 'Nominal wall thickness', 'mm', 'Pressure', 'bar', [
          { name: 'Burst resistance / γ', x: tx, y: cap.map((c) => round(c.b / BAR)) }, { name: 'Collapse resistance / γ', x: tx, y: cap.map((c) => round(c.c / BAR)) },
          { name: 'Propagation resistance / γ', x: tx, y: cap.map((c) => round(c.pr / BAR)) },
          { name: 'Internal overpressure demand', x: [tx[0], tx[tx.length - 1]], y: [round((pli - pe) / BAR), round((pli - pe) / BAR)], dash: true },
          { name: 'External overpressure demand', x: [tx[0], tx[tx.length - 1]], y: [round((pe - pMin) / BAR), round((pe - pMin) / BAR)], dash: true },
        ], 'Factored resistances as the wall thickens. Each resistance curve must lie above its dashed demand line at the chosen thickness; the crossing gives the required wall.'),
      ],
      tables: [tbl('Wall thickness by limit state', [col('ls', 'Limit state'), col('t', 'Required nominal thickness', 'mm'), col('uc', 'Utilisation at the specified wall', '-')], [
        { ls: 'Pressure containment', t: round(tBurst.tNominal * 1e3), uc: round(pcn.uc, 3) }, { ls: 'System collapse', t: round(tColl * 1e3), uc: round(col1.uc, 3) },
        { ls: 'Propagation (to avoid arrestors)', t: round(tProp * 1e3), uc: round(prop.uc, 3) }, { ls: 'ASME B31.8 hoop', t: round((asme.tRequiredHoop + tCorr + tFab) * 1e3), uc: round(asme.ucHoop, 3) },
      ], 'One row per limit state: the nominal wall thickness that just satisfies it and the utilisation of the specified wall. The largest strength requirement sets the wall.')],
      notes: [
        `${gov.name} governs (utilisation ${fmt(gov.uc)}). t1 = ${fmt(tw.t1 * 1e3)} mm in operation, ${fmt(twInst.t1 * 1e3)} mm before operation; characteristic f_y = ${fmt(str.fy / MPA)} MPa at ${fmt(inp.Tdes)} °C.`,
        prop.arrestorsRequired ? 'A buckle would propagate at this depth: provide buckle arrestors or thicken the wall to the propagation requirement.' : 'A buckle cannot propagate at this depth with this wall: no arrestors are needed.',
        !prop.inValidityRange ? 'D/t is outside the 15–45 calibration range of the propagation expression.' : `Partial factors: γ_m ${fmt(gM)}, γ_SC ${fmt(gPc)} (containment) / ${fmt(gSc)} (other), from ${man ? 'manual entry' : 'the recalled tables'}.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} is the governing strength check at ${fmt(gov.uc)}; the required nominal wall is ${fmt(tReq * 1e3)} mm against ${fmt(inp.t)} mm specified.`,
          gov.uc > 1 ? `Increase the wall to at least ${fmt(tReq * 1e3)} mm (next standard size), raise the grade, or lower the design pressure.` : `The wall has ${fmt(-dt * 1e3)} mm of margin over the strength requirement. Confirm the partial factors against the governing edition and check combined loading and installation before reducing it.`,
          `${STD_F101} §5.4.2–5.4.5 with Tables 5-1, 5-2 (recalled — confirm)`),
        adv('economic', 'info', `Each millimetre of wall is ${fmt(steelPerMetre(p, 1e-3) * L / 1e3)} t of steel over ${fmt(inp.Lkm)} km; the specified wall differs from the requirement by ${fmt(-dt * 1e3)} mm (${fmt(-dSteel)} t).`,
          prop.arrestorsRequired ? `Arrestors at regular spacing are normally cheaper than thickening the whole line to ${fmt(tProp * 1e3)} mm (+${fmt(steelPerMetre(p, Math.max(tProp - p.t, 0)) * L / 1e3)} t).` : dt < -0.001 ? `Optimising to the required wall would save about ${fmt(-dSteel)} t of line pipe; weigh it against corrosion-allowance and installation needs.` : 'The wall is matched to the pressure; a higher grade is the remaining lever on tonnage.',
          ECON_BASIS, steelCo2(steelPerMetre(p, 1e-3) * L / 1e3).replace(' t of line-pipe steel', ' t of line-pipe steel per millimetre of wall')),
      ],
      publish: { 'pipeline.D': p.D, 'pipeline.t': p.t, 'pipeline.length': L, 'quantities.steel_t': p.steelKgM * L / 1e3, 'quantities.pipe_km': inp.Lkm },
      quantities: { steel_t: round(p.steelKgM * L / 1e3), pipe_km: inp.Lkm },
    };
  },
  verify: [{
    // D = 323.9 mm, t1 = t = 14.3 mm (no allowances), X65 at 20 °C: f_y = 0.96·450 = 432 MPa, f_u/1.15 = 0.96·535/1.15 = 446.6 MPa
    // → f_cb = 432 MPa. p_b = 2t/(D − t)·f_cb·2/√3 = 2·14.3/309.6·432·1.1547 = 46.08 MPa = 460.8 bar.
    // Elastic collapse p_el = 2E(t/D)³/(1 − ν²) = 2·207e3·(14.3/323.9)³/0.91 = 39.15 MPa = 391.5 bar.
    name: 'Burst resistance and elastic collapse pressure of a 12-inch X65 pipe',
    source: 'DNV-ST-F101 (2021) Eq. 5.8, 5.9, 5.11; hand calculation',
    inputs: { D: 323.9, t: 14.3, tCorr: 0, tFab: 0, Tdes: 20, steel: 'api5l-x65' },
    expect: [{ kpi: 'pb', value: 460.8, tol: 0.003 }, { kpi: 'pel', value: 391.5, tol: 0.003 }],
  }],
};
