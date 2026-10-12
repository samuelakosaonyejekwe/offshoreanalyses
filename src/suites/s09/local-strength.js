// Local cross-sectional strength of a riser pipe section: API RP 2RD von Mises working stress and DNV LRFD
// combined loading under tension, bending, torsion and pressure, with the interaction diagram and required wall.
import { wallStresses } from '../../engine/lines/riser.js';
import { vonMisesApi2rd, combinedLoadingDnv, burstCheckDnv, hoopStress } from '../../engine/lines/risercheck.js';
import { bisect } from '../../engine/math/roots.js';
import {
  G, need, run, fmt, round, pipeInputs, waterInputs, steelTonnes, adv, check, kpi, plot, finish, steelCo2Text, CONFIRM, CO2_STEEL,
} from './_common.js';
import { criteriaInputs } from './_scr.js';

// plastic moment of the default section for the verification case: Mk = fy·αc·(D − t)²·t, αc = 0.6 + 0.4·fu/fy (D/t < 15)
const MK = 450e6 * (0.6 + 0.4 * 535 / 450) * (0.2731 - 0.0206) ** 2 * 0.0206;

export default {
  id: 's09-local-strength', suite: 9, title: 'Local cross-sectional strength',
  summary: 'Wall stresses and code utilisation of a riser cross-section under tension, bending, torsion and internal / external pressure, by API RP 2RD working stress and DNV LRFD combined loading.',
  covers: [15], fidelity: 'analytical',
  standards: ['API RP 2RD (1998, R2013) §5.2.3 von Mises criterion, Table 2 design case factors', 'DNV-ST-F201 (2018) Sec. 5 D combined loading (LRFD), Sec. 5 C bursting', 'DNV-ST-F101 (2021) Sec. 5'],
  models: ['Lamé thick-wall hoop and radial stress', 'True wall tension from effective tension', 'von Mises equivalent stress with torsional shear', 'LRFD combined axial tension–bending–pressure interaction with strain-hardening factor αc', 'Burst resistance'],
  limits: 'Plain pipe away from connectors, welds and discontinuities, linear elastic stresses with plastic section capacities from the code formulas. API STD 2RD method 1 is not implemented: the API RP 2RD (1998) von Mises format and the DNV LRFD criteria are. Torsion enters the von Mises stress only (the DNV interaction equation has no torsion term). Nonlinear local cross-section analysis, ovalisation under combined load and global–local finite-element submodelling are outside the module. Design factors are overridable defaults.',
  publishes: ['riser.localUc', 'riser.requiredWall'],
  searchTerms: 'riser cross section strength von Mises combined loading LRFD API 2RD',
  inputs: [
    ...pipeInputs().filter((i) => ['D', 't', 'mat'].includes(i.key)),
    ...waterInputs().filter((i) => i.key === 'rhoW'),
    { key: 'tCorr', label: 'Corrosion allowance', unit: 'mm', type: 'number', default: 0, min: 0, max: 10, group: 'Riser pipe', role: 'geometry' },
    { key: 'TeF', label: 'Effective tension, functional', unit: 'kN', type: 'number', default: 1190, min: -5000, max: 60000, group: 'Loads at the section', role: 'bc', bind: 'riser.topTension', uncertain: { dist: 'normal', cov: 0.05 } },
    { key: 'TeE', label: 'Effective tension, environmental', unit: 'kN', type: 'number', default: 180, min: -5000, max: 60000, group: 'Loads at the section', role: 'bc', bind: 'riser.dynTensionAmp', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'MF', label: 'Bending moment, functional', unit: 'kN·m', type: 'number', default: 60, min: 0, max: 60000, group: 'Loads at the section', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 } },
    { key: 'ME', label: 'Bending moment, environmental', unit: 'kN·m', type: 'number', default: 70, min: 0, max: 60000, group: 'Loads at the section', role: 'bc', bind: 'riser.hangoffMoment', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'torque', label: 'Torsional moment', unit: 'kN·m', type: 'number', default: 0, min: 0, max: 20000, group: 'Loads at the section', role: 'bc', help: 'From a 3-D global analysis or connection make-up; the 2-D riser models of this suite do not produce torsion.' },
    { key: 'pInt', label: 'Local internal design pressure', unit: 'bar', type: 'number', default: 250, min: 0, max: 1500, group: 'Loads at the section', role: 'bc' },
    { key: 'pMin', label: 'Local minimum internal pressure', unit: 'bar', type: 'number', default: 0, min: 0, max: 1500, group: 'Loads at the section', role: 'bc' },
    { key: 'secDepth', label: 'Water depth of the section', unit: 'm', type: 'number', default: 15, min: 0, max: 3600, group: 'Loads at the section', role: 'bc', help: 'Sets the external pressure at the section.' },
    ...criteriaInputs().map((c) => (c.key === 'Cf' ? { ...c, default: 1.2 } : c)),
    { key: 'gammaF', label: 'Functional load factor γF', unit: '-', type: 'number', default: 1.1, min: 0.9, max: 1.3, group: 'Acceptance criteria', role: 'criteria', help: `ULS 1.1 (DNV-ST-F201 Table 5-2). ${CONFIRM}` },
    { key: 'gammaE', label: 'Environmental load factor γE', unit: '-', type: 'number', default: 1.3, min: 0.7, max: 1.5, group: 'Acceptance criteria', role: 'criteria', help: `ULS 1.3. ${CONFIRM}` },
    { key: 'gammaSC', label: 'Safety-class factor γSC', unit: '-', type: 'number', default: 1.14, min: 1, max: 1.4, group: 'Acceptance criteria', role: 'criteria', help: `1.04 / 1.14 / 1.26 low / normal / high. ${CONFIRM}` },
    { key: 'gammaM', label: 'Material factor γm', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.3, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'f0', label: 'Ovality (Dmax − Dmin)/D', unit: '%', type: 'number', default: 0.5, min: 0.5, max: 3, group: 'Riser pipe', role: 'geometry' },
    { key: 'ucTarget', label: 'Target utilisation for the required wall', unit: '-', type: 'number', default: 0.9, min: 0.5, max: 1, group: 'Acceptance criteria', role: 'criteria' },
  ],

  compute(inp, ctx) {
    // engineering units → SI: mm → m, kN → N, kN·m → N·m, bar → Pa
    const mat = ctx.material(inp.mat); need(mat, `Steel '${inp.mat}' is not in the materials library.`);
    const D = inp.D / 1e3, tN = inp.t / 1e3, t = tN - inp.tCorr / 1e3;
    need(t > 0.002 && tN < D / 2.2, `Corroded wall ${fmt(t * 1e3)} mm is not usable; check wall thickness and corrosion allowance.`);
    const pi = inp.pInt * 1e5, pmin = inp.pMin * 1e5, pe = inp.rhoW * G * inp.secDepth;
    const Te = (inp.TeF + inp.TeE) * 1e3, M = (inp.MF + inp.ME) * 1e3, Tq = inp.torque * 1e3;
    const geo = (tt) => { const Di = D - 2 * tt; return { Ai: Math.PI / 4 * Di * Di, Ae: Math.PI / 4 * D * D, As: Math.PI / 4 * (D * D - Di * Di), I: Math.PI / 64 * (D ** 4 - Di ** 4) }; };
    const evalAt = (tt, o = {}) => {
      const g = geo(tt), te = o.Te ?? Te, m = o.M ?? M;
      const Tw = te + pi * g.Ai - pe * g.Ae;
      const api = vonMisesApi2rd({ T: Tw, M: m, pi, pe, D, t: tt, sigmaY: mat.fy, Cf: inp.Cf, Ca: inp.CaAllow });
      const tau = Tq * (D / 2) / (2 * g.I);
      const vm = Math.sqrt(api.vonMises ** 2 + 3 * tau * tau);
      return { g, Tw, api, tau, vm, ucApi: vm / api.allowable };
    };
    const dnvAt = (tt, o = {}) => combinedLoadingDnv({ D, t: tt, fy: mat.fy, fu: mat.fu, MF: o.MF ?? inp.MF * 1e3, ME: o.ME ?? inp.ME * 1e3, TeF: o.TeF ?? inp.TeF * 1e3, TeE: o.TeE ?? inp.TeE * 1e3, pld: pi, pe, pmin, gammaF: inp.gammaF, gammaE: inp.gammaE, gammaSC: inp.gammaSC, gammaM: inp.gammaM, E: mat.E, nu: mat.nu, f0: inp.f0 / 100 });
    const a = evalAt(t), dnv = run(() => dnvAt(t), 'DNV combined-loading check failed (check D/t ≤ 60 and that the tension is below the plastic capacity)');
    need(Number.isFinite(dnv.uc), `The net internal pressure (${fmt((pi - pe) / 1e6)} MPa) exceeds the burst resistance of the section (${fmt(dnv.pb / 1e6)} MPa): the combined-loading criterion cannot be evaluated. Increase the wall or reduce the pressure.`);
    const ws = wallStresses({ Te, M, pi, pe, D, t });
    const burst = burstCheckDnv({ pli: 1.1 * pi, pe, D, t, fy: mat.fy, fu: mat.fu, gammaSC: inp.gammaSC, gammaM: inp.gammaM });
    const hs = hoopStress({ pi, pe, D, t });
    // required corroded wall for the target utilisation on both criteria
    const f = (tt) => { let d; try { d = dnvAt(tt).uc; } catch { d = 99; } return Math.max(evalAt(tt).ucApi, Number.isFinite(d) ? d : 99) - inp.ucTarget; };
    const lo = 0.003, hi = D / 2.2;
    const tReq = f(hi) > 0 ? hi : f(lo) <= 0 ? lo : bisect(f, lo, hi, { tol: 1e-7 });
    const tReqNom = tReq + inp.tCorr / 1e3;

    const kpis = [
      kpi('vm', 'von Mises stress (incl. torsion)', a.vm / 1e6, 'MPa', 'Equivalent stress from axial (true wall tension ± bending), hoop, radial and torsional shear stress at the governing fibre, against the API RP 2RD allowable Cf·Ca·σy.', a.api.allowable / 1e6, 'max'),
      kpi('ucApi', 'API RP 2RD utilisation', a.ucApi, '-', 'von Mises stress over Cf·Ca·σy. Working-stress format: below 1.0 the section stays elastic with the code margin for the selected design case.', 1, 'max'),
      kpi('ucDnv', 'DNV LRFD combined-loading utilisation', dnv.uc, '-', `Left-hand side of the interaction equation for net ${dnv.mode} overpressure with factored moment and effective tension. It measures closeness to plastic local failure of the cross-section.`, 1, 'max'),
      kpi('Tw', 'True wall tension', a.Tw / 1e3, 'kN', 'Axial force carried by the steel: effective tension plus pi·Ai − pe·Ae. It gives the axial wall stress; effective tension alone would understate it under internal pressure.'),
      kpi('axial', 'Axial stress from true wall tension', ws.axial / 1e6, 'MPa', 'True wall tension over the steel area, uniform across the wall.'),
      kpi('bending', 'Bending stress at the outer fibre', ws.bendingOuter / 1e6, 'MPa', 'M·(D/2)/I. It adds to the axial stress on the tension side and subtracts on the compression side.'),
      kpi('hoop', 'Hoop stress at the bore (Lamé)', ws.hoopInner / 1e6, 'MPa', 'Circumferential stress at the inner surface from internal and external pressure by the thick-wall solution; it is highest at the bore.'),
      kpi('shear', 'Torsional shear stress', a.tau / 1e6, 'MPa', 'Shear stress at the outer surface from the torsional moment, T·r/J. Zero unless a torque is entered.'),
      kpi('Mk', 'Plastic bending moment capacity Mk', dnv.Mk / 1e3, 'kN·m', 'Characteristic plastic moment fy·αc·(D − t)²·t of the section, including the strain-hardening factor αc.'),
      kpi('Tk', 'Plastic axial capacity Tk', dnv.Tk / 1e3, 'kN', 'Characteristic plastic tension capacity fy·αc·π·(D − t)·t of the section.'),
      kpi('burstUc', 'Burst utilisation', burst.uc, '-', 'Local incidental net pressure (1.1 × design) over the factored burst resistance: pressure containment must be satisfied before combined loading is considered.', 1, 'max'),
      kpi('tReq', 'Nominal wall for the target utilisation', tReqNom * 1e3, 'mm', `Smallest nominal wall at which both the API and the DNV utilisation fall to ${inp.ucTarget} under the same loads (the loads themselves change little with wall).`),
      kpi('massPerM', 'Steel mass per metre', steelTonnes(D, tN, 1, mat.rho) * 1e3, 'kg/m', 'Mass of the steel pipe per metre at the nominal wall; multiplied by the riser length it gives the tonnage and, with buoyancy, the top tension.'),
    ];
    const checks = [
      check('von Mises working stress', a.vm / 1e6, a.api.allowable / 1e6, 'MPa', 'API RP 2RD (1998) §5.2.3, Table 2', 'Guards against first yield of the pipe wall under combined tension, bending, pressure and torsion. Beyond it the section is no longer within the elastic design basis.'),
      check(`Combined loading, net ${dnv.mode} overpressure`, dnv.uc, 1, '-', 'DNV-ST-F201 (2018) Sec. 5 D (Eq. 5.28 / 5.30)', 'Guards against plastic collapse or local buckling of the cross-section under factored moment, effective tension and net pressure.'),
      check('Pressure containment (burst)', burst.demand / 1e6, burst.capacity / 1e6, 'MPa', 'DNV-ST-F201 (2018) Sec. 5 C200', 'Guards against rupture of the wall under the local incidental pressure.'),
    ];

    // interaction diagram in the effective tension – moment plane (characteristic loads; DNV curve uses factored loads)
    const Tmax = dnv.Tk / (inp.gammaSC * inp.gammaM), xs = [], mD = [], mA = [];
    for (let i = 0; i <= 24; i++) {
      const T = Tmax * i / 24; xs.push(T / 1e3);
      const gD = (m) => { try { const u = dnvAt(t, { MF: m, ME: 0, TeF: T, TeE: 0 }).uc; return Number.isFinite(u) ? u - 1 : 1; } catch { return 1; } };
      const gA = (m) => evalAt(t, { Te: T, M: m }).ucApi - 1;
      const top = 2 * dnv.Mk;
      mD.push(gD(0) >= 0 ? 0 : gD(top) < 0 ? top / 1e3 : bisect(gD, 0, top, { tol: 1e-6 }) / 1e3);
      mA.push(gA(0) >= 0 ? 0 : gA(top) < 0 ? top / 1e3 : bisect(gA, 0, top, { tol: 1e-6 }) / 1e3);
    }
    // through-wall stress distribution (Lamé) on the tension side of bending
    const ro = D / 2, ri = ro - t, rs = [], sh = [], sr = [], sa = [];
    const A = (pi * ri * ri - pe * ro * ro) / (ro * ro - ri * ri), B = (pi - pe) * ri * ri * ro * ro / (ro * ro - ri * ri);
    for (let i = 0; i <= 20; i++) { const r = ri + t * i / 20; rs.push(r * 1e3); sh.push((A + B / (r * r)) / 1e6); sr.push((A - B / (r * r)) / 1e6); sa.push((ws.axial + M * r / a.g.I) / 1e6); }
    const tw = [], uA = [], uD = []; for (let i = 0; i <= 20; i++) { const tA = Math.max(0.004, 0.5 * t), tB = Math.min(1.8 * t, D / 2.3), tt = tA + (tB - tA) * i / 20; tw.push((tt + inp.tCorr / 1e3) * 1e3); uA.push(evalAt(tt).ucApi); let d; try { d = dnvAt(tt).uc; } catch { d = NaN; } uD.push(Number.isFinite(d) ? Math.min(d, 5) : 5); }
    const plots = [
      plot('Tension–moment interaction diagram', 'Effective tension', 'kN', 'Bending moment capacity', 'kN·m', [{ name: 'DNV LRFD limit (functional loads)', x: xs, y: mD }, { name: 'API RP 2RD von Mises limit', x: xs, y: mA }, { name: 'Operating point (characteristic)', x: [Te / 1e3], y: [M / 1e3] }],
        'Bending moment the section can carry at each effective tension at the stated pressures, for both codes. The operating point must lie inside the curves; its distance to them is the reserve.', 'scatter'),
      plot('Stress through the wall', 'Radius', 'mm', 'Stress', 'MPa', [{ name: 'Hoop', x: rs, y: sh }, { name: 'Radial', x: rs, y: sr }, { name: 'Axial (tension side)', x: rs, y: sa }], 'Hoop, radial and axial stress from the bore to the outer surface on the tension side of bending. Hoop stress peaks at the bore, bending at the outer fibre, so the governing point depends on the load mix.'),
      plot('Utilisation against wall thickness', 'Nominal wall thickness', 'mm', 'Utilisation', '-', [{ name: 'API RP 2RD', x: tw, y: uA }, { name: 'DNV LRFD', x: tw, y: uD }, { name: 'Limit', x: [tw[0], tw[20]], y: [1, 1], dash: true }, { name: 'Selected wall', x: [inp.t], y: [Math.max(a.ucApi, dnv.uc)] }], 'How both utilisations change with wall thickness under the same loads. The crossing of the limit line gives the minimum wall for each code.'),
    ];
    const tables = [{ title: 'Section stresses and capacities', columns: [{ key: 'item', label: 'Quantity' }, { key: 'value', label: 'Value' }, { key: 'unit', label: 'Unit' }],
      rows: [['Effective tension (characteristic)', Te / 1e3, 'kN'], ['True wall tension', a.Tw / 1e3, 'kN'], ['Bending moment (characteristic)', M / 1e3, 'kN·m'], ['Internal pressure', pi / 1e6, 'MPa'], ['External pressure', pe / 1e6, 'MPa'], ['API radial stress', a.api.sigmaR / 1e6, 'MPa'], ['API hoop stress', a.api.sigmaTheta / 1e6, 'MPa'], ['API axial stress (governing side)', a.api.sigmaZ / 1e6, 'MPa'], ['Hoop stress, mean-diameter form', hs.meanDiameter / 1e6, 'MPa'], ['Lamé von Mises, bore', ws.vonMisesInner / 1e6, 'MPa'], ['Lamé von Mises, outer surface', ws.vonMisesOuter / 1e6, 'MPa'], ['Design moment Md', dnv.Md / 1e3, 'kN·m'], ['Design effective tension Ted', dnv.Ted / 1e3, 'kN'], ['Strain-hardening factor αc', dnv.alphaC, '-'], ['Burst resistance pb', dnv.pb / 1e6, 'MPa']].map(([item, value, unit]) => ({ item, value: round(value, 5), unit })),
      explain: 'Stress components and capacities behind the two code checks, in the form needed for a calculation note or for comparison with a local finite-element model.' }];
    const ucMax = Math.max(...checks.map((c) => c.uc)), gov = checks.find((c) => c.uc === ucMax);
    const dWall = (tReqNom - tN) * 1e3, perMm = steelTonnes(D, 0.001, 1000, mat.rho);
    const advice = [
      adv('engineering', ucMax > 1 ? 'fail' : ucMax > 0.9 ? 'warn' : 'ok',
        `${gov.name} governs at ${fmt(ucMax)} (API ${fmt(a.ucApi)}, DNV ${fmt(dnv.uc)}, burst ${fmt(burst.uc)}). von Mises ${fmt(a.vm / 1e6)} MPa: axial ${fmt(ws.axial / 1e6)}, bending ${fmt(ws.bendingOuter / 1e6)}, hoop ${fmt(ws.hoopInner / 1e6)} MPa.`,
        ucMax > 1 ? `Increase the nominal wall to about ${fmt(tReqNom * 1e3)} mm for a utilisation of ${inp.ucTarget}, or select a higher grade: at the same wall, yield of at least ${fmt(mat.fy * ucMax / inp.ucTarget / 1e6)} MPa is needed. If bending dominates (${fmt(ws.bendingOuter / 1e6)} MPa), reducing the moment is more effective than adding wall.`
          : `The section has reserve; a nominal wall of ${fmt(tReqNom * 1e3)} mm would reach ${inp.ucTarget}. Before reducing it, check collapse, fatigue and the connector, which often govern instead of plain-pipe strength.`,
        `API RP 2RD (1998, R2013) §5.2.3 (Cf = ${inp.Cf}, Ca = ${inp.CaAllow}) and DNV-ST-F201 (2018) Sec. 5 D. API STD 2RD method 1 is not implemented. ${CONFIRM}`),
      adv('economic', dWall > 0 ? 'warn' : 'info',
        dWall > 0 ? `Strength needs ${fmt(dWall)} mm more wall: ${fmt(perMm * dWall)} t of steel per kilometre of riser and a proportionally higher top tension.` : `The wall is ${fmt(-dWall)} mm above what strength needs at a utilisation of ${inp.ucTarget}: ${fmt(perMm * -dWall)} t of steel per kilometre that is justified only if another limit state uses it.`,
        'Compare grades before adding wall: a higher-strength line pipe reaches the same utilisation with less steel, lower hang-off load and faster welding, at a higher price per tonne and tighter weld procedure requirements.',
        'Whole-life costing (ISO 15686-5); AACE RP 18R-97 class 4 quantities until the global analysis is final.',
        `Each millimetre of wall is ${fmt(perMm)} t of steel per kilometre ≈ ${fmt(perMm * CO2_STEEL)} t CO₂e (screening factor ${CO2_STEEL} t CO₂e/t). ${dWall > 0 ? steelCo2Text(perMm * dWall) : 'Using the reserve, or a higher grade, reduces embodied carbon and hull payload together.'}`),
    ];
    return finish({
      kpis, checks, plots, tables,
      notes: [`Checks on the corroded wall ${fmt(t * 1e3)} mm; net ${dnv.mode} overpressure in the DNV criterion; ${gov.name} governs (${fmt(ucMax)}).`, 'API utilisation is on characteristic loads with allowable stress; DNV utilisation is on factored loads with factored resistance — they are not expected to agree numerically.'],
      advice,
      publish: { 'riser.localUc': ucMax, 'riser.requiredWall': tReqNom },
      quantities: { steel_t: Math.max(0, perMm * dWall) },
    });
  },

  verify: [
    {
      // Pure tension, no pressure: σ = T/A, A = π/4·(273.1² − 231.9²) = 16341 mm²; T = 1000 kN → 61.20 MPa;
      // API utilisation = σ/(Cf·Ca·σy) = 61.20/(1.0·(2/3)·450) = 0.2040.
      name: 'Pure tension: σ = T/A and utilisation σ/(Cf·Ca·σy)', source: 'API RP 2RD (1998) §5.2.3 (hand calculation)',
      inputs: { TeF: 1000, TeE: 0, MF: 0, ME: 0, pInt: 0, pMin: 0, secDepth: 0, Cf: 1, CaAllow: 2 / 3, tCorr: 0 },
      expect: [{ kpi: 'vm', value: 1000e3 / (Math.PI / 4 * (0.2731 ** 2 - 0.2319 ** 2)) / 1e6, tol: 1e-4 }, { kpi: 'ucApi', value: 1000e3 / (Math.PI / 4 * (0.2731 ** 2 - 0.2319 ** 2)) / (2 / 3 * 450e6), tol: 1e-4 }],
    },
    {
      // Pure functional bending, no pressure or tension (DNV Eq. 5.28 reduces to γSC·γm·γF·MF/Mk):
      // Mk = fy·αc·(D − t)²·t, αc = (1 − β) + β·fu/fy with β = 0.4 for D/t = 13.3 < 15 → αc = 1.0756,
      // Mk = 450e6·1.0756·0.2525²·0.0206 = 635.6 kN·m; MF = 200 kN·m → uc = 1.14·1.15·1.1·200/635.6 = 0.4538.
      name: 'Pure bending in the DNV LRFD criterion: uc = γSC·γm·γF·M/Mk', source: 'DNV-ST-F201 (2018) Sec. 5 D, Eq. 5.28 with Mk = fy·αc·(D − t)²·t',
      inputs: { TeF: 0, TeE: 0, MF: 200, ME: 0, pInt: 0, pMin: 0, secDepth: 0, tCorr: 0 },
      expect: [{ kpi: 'Mk', value: MK / 1e3, tol: 1e-3 }, { kpi: 'ucDnv', value: 1.14 * 1.15 * 1.1 * 200e3 / MK, tol: 1e-3 }],
    },
  ],
};
