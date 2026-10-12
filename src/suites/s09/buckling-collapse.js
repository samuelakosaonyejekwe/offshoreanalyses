// Riser buckling and collapse: external-pressure collapse, propagating buckling, combined bending–tension–external
// pressure, strain-based local buckling, Brazier ovalisation and Euler buckling under effective compression.
import { collapsePressureDnv, collapseCheckDnv, collapsePressureApi, propagatingBuckling, combinedLoadingDnv, localBucklingStrain } from '../../engine/lines/risercheck.js';
import { brazierOvalisation } from '../../engine/pipeline/strength.js';
import { bisect } from '../../engine/math/roots.js';
import {
  G, need, run, fmt, round, pipeInputs, waterInputs, steelTonnes, adv, check, kpi, plot, finish, steelCo2Text, CONFIRM, CO2_STEEL,
} from './_common.js';

export default {
  id: 's09-buckling-collapse', suite: 9, title: 'Riser buckling and collapse',
  summary: 'Collapse under external pressure, propagating buckling, combined bending–tension–external-pressure utilisation, local buckling strain, ovalisation and column buckling of a riser section.',
  covers: [8], fidelity: 'analytical',
  standards: ['DNV-ST-F201 (2018) Sec. 5 D (hoop buckling, propagating buckling, combined loading)', 'DNV-ST-F101 (2021) Sec. 5, Sec. 13', 'API RP 1111 (2015) §4.3.2 collapse and propagation', 'API RP 2RD (1998, R2013) §5.4'],
  models: ['External pressure collapse cubic with ovality (elastic–plastic interaction)', 'Shell elastic buckling pressure 2E(t/D)³/(1 − ν²)', 'Propagating buckle pressure', 'LRFD combined axial tension–bending–pressure interaction', 'Strain-based local buckling', 'Brazier ovalisation', 'Euler column buckling under effective compression'],
  limits: 'Closed-form code criteria for plain pipe: API STD 2RD method 1 is not implemented (the API RP 2RD / RP 1111 collapse forms and the DNV LRFD criteria are). Localised cross-section deformation (wrinkling, dent interaction, buckle arrestor design) needs nonlinear shell finite elements and is outside the module. Safety-class, material and fabrication factors are overridable defaults recalled from the standards.',
  publishes: ['riser.collapseUc', 'riser.collapseDepth'],
  searchTerms: 'riser collapse external pressure propagating buckle local buckling',
  inputs: [
    ...pipeInputs().filter((i) => ['D', 't', 'mat'].includes(i.key)),
    ...waterInputs(),
    { key: 'tCorr', label: 'Corrosion allowance', unit: 'mm', type: 'number', default: 1.5, min: 0, max: 10, group: 'Riser pipe', role: 'geometry', help: 'Deducted from the wall for all checks (end of life).' },
    { key: 'tFab', label: 'Fabrication under-tolerance', unit: '%', type: 'number', default: 8, min: 0, max: 15, group: 'Riser pipe', role: 'geometry', help: 'Negative wall tolerance; deducted (with corrosion) for the pure collapse check (t1). Seamless pipe up to 12.5 %, welded about 5–8 %.' },
    { key: 'f0', label: 'Initial ovality (Dmax − Dmin)/D', unit: '%', type: 'number', default: 0.5, min: 0.5, max: 3, group: 'Riser pipe', role: 'geometry', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Not to be taken below 0.5 %. Reeling or bending increases it.' },
    { key: 'alphaFab', label: 'Fabrication factor αfab', type: 'select', default: '1', group: 'Riser pipe', role: 'criteria', options: [{ value: '1', label: 'Seamless (1.00)' }, { value: '0.93', label: 'UO / TRB / ERW (0.93)' }, { value: '0.85', label: 'UOE (0.85)' }] },
    { key: 'pmin', label: 'Minimum internal pressure at the section', unit: 'bar', type: 'number', default: 0, min: 0, max: 1500, group: 'Loads at the section', role: 'bc', help: '0 for an empty or depressurised riser (installation, blow-down): the governing collapse case.' },
    { key: 'TeF', label: 'Effective tension, functional', unit: 'kN', type: 'number', default: 240, min: -5000, max: 50000, group: 'Loads at the section', role: 'bc', bind: 'riser.tdpTension', help: 'Static effective tension at the section (touchdown zone by default). Negative = effective compression.' },
    { key: 'TeE', label: 'Effective tension, environmental', unit: 'kN', type: 'number', default: 60, min: -5000, max: 50000, group: 'Loads at the section', role: 'bc' },
    { key: 'MF', label: 'Bending moment, functional', unit: 'kN·m', type: 'number', default: 105, min: 0, max: 50000, group: 'Loads at the section', role: 'bc', bind: 'riser.maxMoment', uncertain: { dist: 'normal', cov: 0.1 } },
    { key: 'ME', label: 'Bending moment, environmental', unit: 'kN·m', type: 'number', default: 40, min: 0, max: 50000, group: 'Loads at the section', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'comp', label: 'Dynamic effective compression', unit: 'kN', type: 'number', default: 0, min: 0, max: 20000, group: 'Loads at the section', role: 'bc', help: 'Peak effective compression predicted by the dynamic analysis, if any.' },
    { key: 'Le', label: 'Effective buckling length', unit: 'm', type: 'number', default: 40, min: 1, max: 500, group: 'Loads at the section', role: 'bc', help: 'Half-wave length available for column buckling (distance between points of contraflexure in the sag bend).' },
    { key: 'gammaSC', label: 'Safety-class factor γSC', unit: '-', type: 'number', default: 1.14, min: 1, max: 1.4, group: 'Acceptance criteria', role: 'criteria', help: `1.04 / 1.14 / 1.26 for low / normal / high safety class (DNV-ST-F201 Table 5-3). ${CONFIRM}` },
    { key: 'gammaM', label: 'Material factor γm', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.3, group: 'Acceptance criteria', role: 'criteria', help: `1.15 for ULS and ALS. ${CONFIRM}` },
    { key: 'gammaF', label: 'Functional load factor γF', unit: '-', type: 'number', default: 1.1, min: 0.9, max: 1.3, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'gammaE', label: 'Environmental load factor γE', unit: '-', type: 'number', default: 1.3, min: 0.7, max: 1.5, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'gammaC', label: 'Propagation condition factor γc', unit: '-', type: 'number', default: 1.0, min: 0.9, max: 1, group: 'Acceptance criteria', role: 'criteria', help: '1.0 = no propagation allowed; 0.9 where propagation over a short distance is accepted (buckle arrestors).' },
    { key: 'gammaEps', label: 'Strain resistance factor γε', unit: '-', type: 'number', default: 2.5, min: 1.5, max: 3.5, group: 'Acceptance criteria', role: 'criteria', help: `2.0 / 2.5 / 3.3 for low / medium / high safety class (DNV-ST-F101). ${CONFIRM}` },
    { key: 'foApi', label: 'API collapse factor fo', unit: '-', type: 'number', default: 0.7, min: 0.5, max: 0.8, group: 'Acceptance criteria', role: 'criteria', help: `0.7 seamless / ERW, 0.6 cold-expanded pipe (API RP 1111). ${CONFIRM}` },
  ],

  compute(inp, ctx) {
    // engineering units → SI: mm → m, bar → Pa, kN → N, % → fraction
    const mat = ctx.material(inp.mat); need(mat, `Steel '${inp.mat}' is not in the materials library.`);
    const D = inp.D / 1e3, t = inp.t / 1e3, t2 = t - inp.tCorr / 1e3, t1 = t * (1 - inp.tFab / 100) - inp.tCorr / 1e3;
    need(t < D / 2.2, `Wall thickness ${inp.t} mm is too large for a ${inp.D} mm pipe (the bore would close). Reduce the wall or increase the diameter.`);
    need(t1 > 0.002, `Wall after tolerance and corrosion allowance is ${fmt(t1 * 1e3)} mm: too thin to assess. Increase the wall.`);
    need(D / t2 <= 60, `D/t = ${fmt(D / t2)} exceeds 60, outside the range of the combined-loading criteria.`);
    const f0 = inp.f0 / 100, aF = Number(inp.alphaFab), E = mat.E, nu = mat.nu, fy = mat.fy, fu = mat.fu;
    const pe = inp.rhoW * G * inp.depth, pmin = inp.pmin * 1e5;
    const col = collapsePressureDnv({ D, t: t1, fy, E, nu, f0, alphaFab: aF });
    const cc = collapseCheckDnv({ D, t: t1, fy, E, nu, f0, alphaFab: aF, pe, pmin, gammaSC: inp.gammaSC, gammaM: inp.gammaM });
    const api = collapsePressureApi({ D, t: t2, S: fy, E, nu });
    const prop = propagatingBuckling({ D, t: t2, fy, alphaFab: aF, pe, pmin, gammaC: inp.gammaC, gammaM: inp.gammaM, gammaSC: inp.gammaSC });
    const loads = { MF: inp.MF * 1e3, ME: inp.ME * 1e3, TeF: inp.TeF * 1e3, TeE: inp.TeE * 1e3 };
    const cl = (tt, o = {}) => combinedLoadingDnv({ D, t: tt, fy, fu, ...loads, pld: pmin, pe, pmin, gammaF: inp.gammaF, gammaE: inp.gammaE, gammaSC: inp.gammaSC, gammaM: inp.gammaM, E, nu, f0, alphaFab: aF, ...o });
    const comb = run(() => cl(t2), 'Combined-loading check failed');
    need(Number.isFinite(comb.uc), `The minimum internal pressure (${inp.pmin} bar) exceeds the burst resistance of the section net of external pressure: the combined-loading criterion cannot be evaluated. This module assesses external-pressure cases; use the local strength module for internal overpressure.`);
    const I = Math.PI / 64 * (D ** 4 - (D - 2 * t2) ** 4), EI = E * I, Mtot = (inp.MF + inp.ME) * 1e3;
    const kappa = Mtot / EI, strain = kappa * D / 2;
    const lb = localBucklingStrain({ strain, D, t: t2, pmin, pe, pb: comb.pb, gammaEps: inp.gammaEps });
    const br = brazierOvalisation({ D, t: t2, curvature: kappa, E, nu });
    const Pcr = Math.PI ** 2 * EI / inp.Le ** 2, comp = inp.comp * 1e3 + Math.max(0, -(inp.TeF + inp.TeE) * 1e3);
    const gam = inp.gammaSC * inp.gammaM;
    const depthCol = (col.pc / gam + pmin) / (inp.rhoW * G);
    // wall thickness (nominal) at which the DNV collapse check is just met
    const ucT = (tn) => collapseCheckDnv({ D, t: tn * (1 - inp.tFab / 100) - inp.tCorr / 1e3, fy, E, nu, f0, alphaFab: aF, pe, pmin, gammaSC: inp.gammaSC, gammaM: inp.gammaM }).uc - 1;
    const tLo = (inp.tCorr / 1e3 + 0.0025) / (1 - inp.tFab / 100), tHi = D / 2.2;
    const tReq = pe <= pmin || ucT(tLo) <= 0 ? tLo : ucT(tHi) > 0 ? tHi : bisect(ucT, tLo, tHi, { tol: 1e-7 });

    const kpis = [
      kpi('pe', 'External pressure at the section', pe / 1e6, 'MPa', 'Hydrostatic pressure at the water depth considered. With the minimum internal pressure it is the load that tries to flatten the pipe.'),
      kpi('pc', 'Characteristic collapse pressure (DNV)', col.pc / 1e6, 'MPa', 'Net external pressure at which the as-built pipe (minimum wall, initial ovality) collapses. It lies below both the elastic and the plastic collapse pressures.'),
      kpi('pel', 'Elastic collapse pressure', col.pel / 1e6, 'MPa', 'Buckling pressure of a perfect elastic ring, 2E(t/D)³/(1 − ν²). It governs thin-walled pipe (high D/t).'),
      kpi('pp', 'Plastic collapse pressure', col.pp / 1e6, 'MPa', 'Pressure at which the wall yields in hoop compression, 2·fy·αfab·t/D. It governs thick-walled pipe (low D/t).'),
      kpi('pcApi', 'Collapse pressure (API RP 1111)', api.Pc / 1e6, 'MPa', 'Collapse pressure from the API interaction of yield and elastic collapse, without an explicit ovality term; shown for comparison with the DNV value.'),
      kpi('collapseUc', 'Collapse utilisation (DNV)', cc.uc, '-', 'Net external pressure over the factored collapse resistance pc/(γSC·γm). Above 1.0 the pipe does not have the code margin against collapse when empty.', 1, 'max'),
      kpi('collapseDepth', 'Water depth at the collapse limit', depthCol, 'm', 'Depth at which the factored collapse resistance equals the net external pressure: the deepest this wall can go empty.'),
      kpi('tReq', 'Nominal wall needed for collapse', tReq * 1e3, 'mm', 'Smallest nominal wall thickness that satisfies the DNV collapse check at this depth with the stated tolerance, corrosion allowance and ovality.'),
      kpi('ppr', 'Propagating buckle pressure (DNV)', prop.pprDnv / 1e6, 'MPa', 'Pressure above which a buckle, once started by a dent or bending, runs along the pipe. It is much lower than the collapse pressure, so deep risers may rely on local buckling being prevented.'),
      kpi('propUc', 'Propagating buckle utilisation', prop.uc, '-', 'Net external pressure over the factored propagation pressure. Above 1.0 a local buckle would propagate; buckle arrestors or the acceptance of that risk are then needed.'),
      kpi('combUc', 'Combined loading utilisation (bending, tension, pressure)', comb.uc, '-', `Left-hand side of the DNV LRFD interaction for net ${comb.mode} overpressure with factored moment and effective tension. It must not exceed 1.0.`, 1, 'max'),
      kpi('strainUc', 'Local buckling strain utilisation', lb.uc, '-', 'Bending strain over the factored critical compressive strain for displacement-controlled bending (e.g. sag bend on the seabed).', 1, 'max'),
      kpi('ovality', 'Bending-induced ovalisation (Brazier)', br.ovality * 100, '%', 'Flattening of the cross-section caused by the bending curvature. It adds to the initial ovality and reduces collapse resistance in the sag bend.'),
      kpi('Pcr', 'Euler buckling load of the section', Pcr / 1e3, 'kN', 'Effective compression at which a pinned strut of the effective length buckles sideways, π²EI/Le². Compare with any dynamic compression at the touchdown.'),
      kpi('dOverT', 'D/t (corroded wall)', D / t2, '-', 'Diameter-to-thickness ratio. Below about 15 collapse is yield-governed; above about 30 it is elastic-buckling governed and very sensitive to ovality.'),
    ];
    const checks = [
      check('System collapse (external pressure)', cc.demand / 1e6, cc.capacity / 1e6, 'MPa', 'DNV-ST-F201 (2018) Sec. 5 D200; DNV-ST-F101 (2021) Eq. 5.10', 'Guards against flattening of the pipe under external pressure when the internal pressure is at its minimum. Collapse floods or severs the riser.'),
      check('Collapse (API RP 1111)', Math.max(0, pe - pmin) / 1e6, inp.foApi * api.Pc / 1e6, 'MPa', 'API RP 1111 (2015) §4.3.2.1: (Po − Pi) ≤ fo·Pc', 'Independent collapse criterion in the API format with collapse factor fo; it should agree in outcome with the DNV check.'),
      check('Propagating buckling', Math.max(0, pe - pmin) / 1e6, prop.capacity / 1e6, 'MPa', 'DNV-ST-F201 (2018) Sec. 5 D300; DNV-ST-F101 (2021) Eq. 5.16', 'Checks whether a local buckle would run along the riser. If it fails, a single dent or wrinkle can flatten a long length unless arrestors stop it.'),
      check(`Combined loading, net ${comb.mode} overpressure`, comb.uc, 1, '-', 'DNV-ST-F201 (2018) Sec. 5 D400–D500 (Eq. 5.28 / 5.30)', 'Interaction of factored bending moment, effective tension and net pressure. Exceeding 1.0 means local buckling of the wall under the combined load.'),
      check('Local buckling strain (displacement-controlled)', strain * 100, lb.capacity * 100, '%', 'DNV-ST-F101 (2021) Sec. 5 Eq. 5.30–5.31', 'Limits compressive bending strain where curvature is imposed by geometry. Beyond the critical strain the compression side wrinkles.'),
    ];
    if (comp > 0) checks.push(check('Column buckling under effective compression', comp / 1e3, Pcr / gam / 1e3, 'kN', 'Euler; DNV-ST-F201 (2018) Sec. 5 (effective compression)', 'Compares the effective compression with the factored Euler load. Above 1.0 the pipe bows sideways, which sharply increases bending at the touchdown.'));

    // collapse resistance against wall thickness
    const ts = [], pcs = [], pas = []; for (let i = 0; i <= 30; i++) { const tn = Math.max(tLo, 0.5 * t) + (1.6 * t - Math.max(tLo, 0.5 * t)) * i / 30; ts.push(tn * 1e3); pcs.push(collapsePressureDnv({ D, t: tn * (1 - inp.tFab / 100) - inp.tCorr / 1e3, fy, E, nu, f0, alphaFab: aF }).pc / gam / 1e6); pas.push(inp.foApi * collapsePressureApi({ D, t: tn - inp.tCorr / 1e3, S: fy, E, nu }).Pc / 1e6); }
    const fs = [], pf = []; for (let i = 0; i <= 20; i++) { const f = 0.005 + 0.025 * i / 20; fs.push(f * 100); pf.push(collapsePressureDnv({ D, t: t1, fy, E, nu, f0: f, alphaFab: aF }).pc / 1e6); }
    // moment capacity against external pressure at the given effective tension (combined criterion = 1)
    const pcl = col.pc / gam, ps = [], ms = [];
    for (let i = 0; i <= 16; i++) {
      const p = pmin + 0.98 * pcl * i / 16; ps.push(p / 1e6);
      const g = (m) => combinedLoadingDnv({ D, t: t2, fy, fu, MF: m, ME: 0, TeF: loads.TeF, TeE: loads.TeE, pld: pmin, pe: p, pmin, gammaF: 1, gammaE: inp.gammaE, gammaSC: inp.gammaSC, gammaM: inp.gammaM, E, nu, f0, alphaFab: aF }).uc - 1;
      let m = 0; try { const hi = 3 * fy * D * D * t2; m = g(0) >= 0 ? 0 : g(hi) < 0 ? hi : bisect(g, 0, hi, { tol: 1e-6 }); } catch { m = 0; }
      ms.push(m / 1e3);
    }
    const Md = (inp.gammaF * inp.MF + inp.gammaE * inp.ME);
    const plots = [
      plot('Collapse resistance against wall thickness', 'Nominal wall thickness', 'mm', 'Allowable net external pressure', 'MPa', [{ name: 'DNV pc/(γSC·γm)', x: ts, y: pcs }, { name: 'API fo·Pc', x: ts, y: pas }, { name: 'Net external pressure', x: [ts[0], ts[30]], y: [Math.max(0, pe - pmin) / 1e6, Math.max(0, pe - pmin) / 1e6], dash: true }, { name: 'Selected wall', x: [inp.t], y: [cc.capacity / 1e6] }],
        'Factored collapse resistance as the wall changes, with the pressure to be resisted (dashed). The crossing gives the minimum wall; the marker is the selected wall.'),
      plot('Collapse pressure against ovality', 'Initial ovality', '%', 'Characteristic collapse pressure', 'MPa', [{ name: 'pc', x: fs, y: pf }, { name: 'External pressure', x: [0.5, 3], y: [pe / 1e6, pe / 1e6], dash: true }], 'Sensitivity of collapse pressure to out-of-roundness. A steep fall means fabrication ovality and bending-induced ovalisation must be tightly controlled.'),
      plot('Moment capacity against external pressure', 'External pressure', 'MPa', 'Design bending moment capacity', 'kN·m', [{ name: 'Combined-loading limit', x: ps, y: ms }, { name: 'Design point', x: [pe / 1e6], y: [Md] }], 'Interaction diagram at the given effective tension: the factored bending moment the section can carry at each external pressure. The design point must lie inside the curve.', 'scatter'),
    ];
    const tables = [{ title: 'Buckling and collapse summary', columns: [{ key: 'item', label: 'Quantity' }, { key: 'value', label: 'Value' }, { key: 'unit', label: 'Unit' }],
      rows: [['Wall for collapse t1 (tolerance + corrosion)', t1 * 1e3, 'mm'], ['Wall for combined loading t2 (corrosion)', t2 * 1e3, 'mm'], ['External pressure', pe / 1e6, 'MPa'], ['Minimum internal pressure', pmin / 1e6, 'MPa'], ['Elastic collapse pressure', col.pel / 1e6, 'MPa'], ['Plastic collapse pressure', col.pp / 1e6, 'MPa'], ['Collapse pressure (DNV)', col.pc / 1e6, 'MPa'], ['Collapse pressure (API)', api.Pc / 1e6, 'MPa'], ['Propagation pressure (DNV)', prop.pprDnv / 1e6, 'MPa'], ['Propagation pressure (API)', prop.ppApi / 1e6, 'MPa'], ['Plastic moment capacity Mk', comb.Mk / 1e3, 'kN·m'], ['Plastic tension capacity Tk', comb.Tk / 1e3, 'kN'], ['Design moment Md', comb.Md / 1e3, 'kN·m'], ['Design effective tension Ted', comb.Ted / 1e3, 'kN'], ['Bending strain', strain * 100, '%'], ['Critical strain / γε', lb.capacity * 100, '%'], ['Brazier limit moment', br.brazierMoment / 1e3, 'kN·m']].map(([item, value, unit]) => ({ item, value: round(value, 5), unit })),
      explain: 'Resistances and load effects behind each check, for the calculation note. t1 is used for pure collapse, t2 for propagation and combined loading.' }];
    const ucMax = Math.max(...checks.map((c) => c.uc)), gov = checks.find((c) => c.uc === ucMax);
    const dT = tReq - t, perMm = steelTonnes(D, 0.001, 1000, mat.rho);
    const advice = [
      adv('engineering', ucMax > 1 ? 'fail' : ucMax > 0.85 ? 'warn' : 'ok',
        `${gov.name} governs at ${fmt(ucMax)}. Collapse resistance ${fmt(cc.capacity / 1e6)} MPa against ${fmt(cc.demand / 1e6)} MPa at ${inp.depth} m (D/t = ${fmt(D / t2)}); the wall is good to ${fmt(depthCol)} m when empty.`,
        ucMax > 1 ? (gov.name.startsWith('Propagating') ? `Either fit buckle arrestors and accept γc = 0.9, or increase the wall: propagation needs about ${fmt(D * ((Math.max(0, pe - pmin) * inp.gammaC * gam) / (35 * fy * aF)) ** 0.4 * 1e3 + inp.tCorr)} mm.`
          : gov.name.startsWith('System') || gov.name.startsWith('Collapse') ? `Increase the nominal wall to at least ${fmt(tReq * 1e3)} mm, tighten ovality below ${inp.f0} %, or keep a minimum internal pressure of about ${fmt((pe - cc.capacity) / 1e5)} bar in all conditions.`
            : `Reduce bending at the section (now ${fmt(Mtot / 1e3)} kN·m) or increase the wall; local buckling capacity grows roughly with t².`)
          : `Reserve against collapse is ${fmt(1 / Math.max(cc.uc, 1e-9))}. Confirm that the installation (empty) and blow-down cases really set the minimum internal pressure, and that measured ovality after reeling or bending stays below ${inp.f0} %.`,
        `DNV-ST-F201 (2018) Sec. 5 D and DNV-ST-F101 (2021) Sec. 5; API RP 1111 (2015) §4.3.2 for comparison. API STD 2RD method 1 is not implemented. ${CONFIRM}`),
      adv('economic', dT > 0 ? 'warn' : 'info',
        dT > 0 ? `Collapse needs ${fmt(dT * 1e3)} mm more wall: ${fmt(perMm * dT * 1e3)} t of steel per kilometre of riser, plus the extra top tension and welding time that go with it.`
          : `The wall is ${fmt(-dT * 1e3)} mm above the collapse minimum (${fmt(tReq * 1e3)} mm). If no other limit state needs it, that margin is ${fmt(perMm * -dT * 1e3)} t of steel per kilometre and the hull payload to carry it.`,
        'Check whether a higher grade, seamless pipe (αfab = 1.0) or tighter ovality reaches the same collapse margin with less wall before adding steel; in deep water collapse often sets the wall and therefore the top tension.',
        'Whole-life costing (ISO 15686-5): wall thickness drives pipe tonnage, welding time, lay tension and hull payload together.',
        `Each millimetre of wall is ${fmt(perMm)} t of steel per kilometre ≈ ${fmt(perMm * CO2_STEEL)} t CO₂e (screening factor ${CO2_STEEL} t CO₂e/t). ${dT > 0 ? steelCo2Text(perMm * dT * 1e3) : 'A higher grade reaches the same plastic collapse pressure with less mass.'}`),
    ];
    return finish({
      kpis, checks, plots, tables,
      notes: [`Collapse assessed on t1 = ${fmt(t1 * 1e3)} mm, combined loading on t2 = ${fmt(t2 * 1e3)} mm; net ${comb.mode} overpressure governs the combined criterion.`, `${gov.name} is the governing check (${fmt(ucMax)}).`],
      advice,
      publish: { 'riser.collapseUc': cc.uc, 'riser.collapseDepth': depthCol },
      quantities: { steel_t: Math.max(0, perMm * dT * 1e3) },
    });
  },

  verify: [
    {
      // No tolerance or corrosion: t1 = 20.6 mm, D = 273.1 mm, X65 (E = 207 GPa, ν = 0.3, fy = 450 MPa).
      // pel = 2E(t/D)³/(1 − ν²) = 2·207e9·(0.0206/0.2731)³/0.91 = 195.3 MPa; pp = 2·fy·t/D = 67.89 MPa.
      name: 'Elastic and plastic collapse pressures', source: 'DNV-ST-F101 (2021) Eq. 5.11–5.12; Timoshenko & Gere, Theory of Elastic Stability §7.5',
      inputs: { tCorr: 0, tFab: 0 },
      expect: [{ kpi: 'pel', value: 2 * 207e9 * (0.0206 / 0.2731) ** 3 / (1 - 0.09) / 1e6, tol: 1e-4 }, { kpi: 'pp', value: 2 * 450 * 0.0206 / 0.2731, tol: 1e-4 },
        // API RP 1111: Pc = Py·Pe/sqrt(Py² + Pe²) with Py = pp and Pe = pel above
        { kpi: 'pcApi', value: (() => { const pe = 2 * 207e9 * (0.0206 / 0.2731) ** 3 / 0.91 / 1e6, py = 2 * 450 * 0.0206 / 0.2731; return py * pe / Math.hypot(py, pe); })(), tol: 1e-4 }],
    },
    {
      // Euler load of a pinned strut: Pcr = π²EI/L², I = π/64·(D⁴ − Di⁴) with D = 273.1, Di = 231.9 mm, L = 40 m.
      name: 'Euler buckling load π²EI/L²', source: 'Euler; Timoshenko & Gere, Theory of Elastic Stability §2.1',
      inputs: { tCorr: 0, tFab: 0 },
      expect: [{ kpi: 'Pcr', value: Math.PI ** 2 * 207e9 * Math.PI / 64 * (0.2731 ** 4 - 0.2319 ** 4) / 40 ** 2 / 1e3, tol: 1e-4 }],
    },
  ],
};
