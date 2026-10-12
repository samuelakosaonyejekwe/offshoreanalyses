// Seabed bearing capacity of a mudmat under combined vertical, horizontal and moment loading.
import { bearingCapacity, terzaghiBearing, skemptonNc, isoUndrainedCapacity, bearingNq, bearingNc, bearingNgamma } from '../../engine/geotech/bearing.js';
import { mudmatClassicalCapacity, requiredMatSize } from '../../engine/geotech/mudmat.js';
import { druckerPragerParams, mohrCoulombTriaxialParams } from '../../engine/geotech/soilmodels.js';
import { DEG, need, run, fmt, round, adv, check, steelCo2Text, matInputs, soilInputs, loadInputs, soilFrom, CONFIRM, NOT_FE } from './_common.js';

const STEEL_T_PER_M2 = 0.16; // typical stiffened mudmat plating and skirts per plan area (screening); the frame is counted by the in-place module

export default {
  id: 's04-bearing-capacity', suite: 4, title: 'Seabed bearing capacity',
  summary: 'Ultimate and design bearing capacity of the mudmat by Terzaghi, Meyerhof, Brinch Hansen, Vesic, Skempton and ISO 19901-4, with the effective-area method for eccentric and inclined load.',
  covers: [3], fidelity: 'analytical',
  standards: ['ISO 19901-4:2016 §7.2–7.3', 'API RP 2GEO (2011) §7', 'DNV-RP-C212 (2019) §5', 'API RP 2A-WSD (21st ed.) §6.13'],
  models: ['Terzaghi bearing capacity', 'Meyerhof bearing capacity', 'Brinch Hansen bearing capacity', 'Vesic bearing capacity', 'Skempton undrained bearing capacity', 'Meyerhof effective-area method', 'Davis & Booker strength-gradient correction (user factor)', 'Mohr–Coulomb and matched Drucker–Prager strength'],
  limits: `Homogeneous soil below the base (undrained strength taken at B′/4 below skirt tip for a linearly increasing profile); punch-through, squeezing and layered failure are not checked here. The Davis & Booker factor F is a chart value supplied by the user. ${NOT_FE}`,
  publishes: ['subsea.bearingUc', 'subsea.Vult', 'subsea.matB', 'subsea.matL'],
  searchTerms: 'mudmat bearing capacity undrained clay effective area inclined eccentric load',
  inputs: [
    ...matInputs(),
    { key: 'D', label: 'Skirt depth (base embedment)', unit: 'm', type: 'number', default: 0.5, min: 0, max: 10, group: 'Mudmat geometry', role: 'geometry', design: true, help: 'The bearing base is taken at skirt-tip level. 0 for an unskirted mat.' },
    ...soilInputs(),
    ...loadInputs(),
    { key: 'method', label: 'Design method', type: 'select', default: 'hansen', group: 'Method and criteria', role: 'control',
      options: [{ value: 'hansen', label: 'Brinch Hansen (1970)' }, { value: 'vesic', label: 'Vesic (1973)' }, { value: 'meyerhof', label: 'Meyerhof (1963)' }, { value: 'iso', label: 'ISO 19901-4 undrained (clay)' }] },
    { key: 'Fdb', label: 'Davis & Booker factor F (ISO method)', unit: '-', type: 'number', default: 1.0, min: 1, max: 2, group: 'Method and criteria', role: 'criteria', showIf: { key: 'method', in: ['iso'] },
      help: 'Chart value for strength increasing with depth (1.0 uniform; up to about 1.5 for a rough base at large κB/su0). Not built in — read it from the chart.' },
    { key: 'gammaF', label: 'Load factor γf', unit: '-', type: 'number', default: 1.35, min: 1, max: 2, group: 'Method and criteria', role: 'criteria', help: `Set γf = 1 and use γm as a global factor for working-stress design. ${CONFIRM}` },
    { key: 'gammaM', label: 'Material factor γm', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Method and criteria', role: 'criteria', help: `ISO 19901-4: 1.5 undrained, 1.25 on tan φ′ drained. ${CONFIRM}` },
  ],
  compute(inp) {
    // SI: kN → N, kN·m → N·m, kPa → Pa, kN/m³ → N/m³, deg → rad.
    const { B, L, D } = inp;
    const s = soilFrom(inp);
    const V = inp.V * 1e3 * inp.gammaF, H = inp.H * 1e3 * inp.gammaF, M = inp.M * 1e3 * inp.gammaF;
    const eB = M / V;
    need(eB < B / 2, `Load eccentricity M/V = ${fmt(eB)} m lies outside the mudmat half-width ${fmt(B / 2)} m: the mat overturns. Increase B or reduce the moment.`);
    const Beff = B - 2 * eB, Aeff = Math.min(Beff, L) * Math.max(Beff, L);
    const zSu = D + Math.min(Beff, L) / 4;
    const su = s.suAt(zSu);
    const coh = s.clay ? su : s.c;
    const base = { B, L, D, c: coh, phi: s.phi, gamma: s.gamma, V, HB: H, eB };
    const tryCap = (method) => { try { return bearingCapacity({ ...base, method }); } catch (e) { return { error: e.message }; } };
    const res = { hansen: tryCap('hansen'), vesic: tryCap('vesic'), meyerhof: tryCap('meyerhof') };
    if (s.clay) {
      try {
        const bE = Math.min(Beff, L), lE = Math.max(Beff, L);
        const r = isoUndrainedCapacity({ Beff: bE, Leff: lE, su0: s.suAt(D), kappa: s.k, D, H, F: inp.Fdb });
        res.iso = { qult: r.qd + s.gamma * D, Aeff: r.Aeff, Qult: (r.qd + s.gamma * D) * r.Aeff };
      } catch (e) { res.iso = { error: e.message }; }
    }
    need(inp.method !== 'iso' || s.clay, 'The ISO 19901-4 undrained expression applies to clay: choose Hansen, Vesic or Meyerhof for sand.');
    const sel = res[inp.method];
    need(sel && !sel.error, `Bearing capacity (${inp.method}) cannot be evaluated: ${sel ? sel.error : 'method unavailable'}. The horizontal load probably exceeds the sliding resistance — check sliding first or enlarge the mat.`);
    // Net design capacity: strength-derived part divided by γm, overburden at base level carried at face value.
    const q0 = s.gamma * D;
    const Vd = (r) => ((r.qult - q0) / inp.gammaM + q0) * r.Aeff;
    const Vcap = Vd(sel);
    const uc = Vcap > 0 ? Math.min(V / Vcap, 99) : 99;
    // Vertical-centric reference values for the classical comparisons.
    const strip = terzaghiBearing({ shape: 'strip', B, c: coh, phi: s.phi, gamma: s.gamma, q: q0 });
    const terz = terzaghiBearing({ shape: Math.abs(L / B - 1) < 0.05 ? 'square' : 'strip', B, c: coh, phi: s.phi, gamma: s.gamma, q: q0 });
    const skNc = skemptonNc({ B: Math.min(B, L), L: Math.max(B, L), D });
    const skQ = s.clay ? skNc * su + q0 : NaN;
    // Capacity as a function of horizontal load (interaction) with the operating point.
    const hs = [], vs = [];
    const Hlim = s.clay ? su * Aeff : V * Math.tan(s.phi) + s.c * Aeff;
    for (let i = 0; i <= 24; i++) {
      const h = Hlim * 0.98 * i / 24;
      try { const r = bearingCapacity({ ...base, HB: h, method: inp.method === 'iso' ? 'hansen' : inp.method }); hs.push(h / 1e3); vs.push(Vd(r) / 1e3); } catch { /* beyond the validity of the inclination factors */ }
    }
    // Capacity against width, and the smallest mat meeting the bearing and sliding checks.
    const bs = [], ucB = [];
    for (let i = 0; i <= 16; i++) {
      const b = B * (0.5 + i * 0.0625);
      try { const r = bearingCapacity({ ...base, B: b, L: L * b / B, method: inp.method === 'iso' ? 'hansen' : inp.method }); const vc = ((r.qult - q0) / inp.gammaM + q0) * r.Aeff; bs.push(round(b)); ucB.push(vc > 0 ? Math.min(V / vc, 99) : 99); } catch { /* too small */ }
    }
    let req = null;
    try { req = requiredMatSize({ V, HB: H, MB: M, su: s.clay ? su : undefined, phi: s.phi, c: s.c, gamma: s.gamma, D, method: inp.method === 'vesic' ? 'vesic' : 'hansen', materialFactor: inp.gammaM, aspect: L / B, targetUc: 1, Bmin: 0.5, Bmax: 80 }); } catch { req = null; }
    const classical = run(() => { try { return mudmatClassicalCapacity({ B, L, V, HB: H, MB: M, su: s.clay ? su : undefined, phi: s.phi, c: s.c, gamma: s.gamma, D, method: inp.method === 'vesic' ? 'vesic' : 'hansen', materialFactor: inp.gammaM }); } catch { return null; } }, 'Mudmat capacity');
    const dp = !s.clay ? druckerPragerParams({ phi: s.phi, c: s.c }) : null;
    const mc = !s.clay ? mohrCoulombTriaxialParams({ phi: s.phi, c: s.c }) : null;
    const checks = [check(`Bearing capacity (${inp.method})`, V / 1e3, Vcap / 1e3, 'kN', 'ISO 19901-4:2016 §7.2–7.3 / API RP 2GEO §7.3', 'Guards against general shear failure of the soil under the mat. If exceeded the mat plunges and rotates until the load is shed or the soil strength at depth arrests it.')];
    if (classical) checks.push(check('Sliding on the effective area', classical.Hequiv / 1e3, classical.Hult / 1e3, 'kN', 'ISO 19901-4:2016 §7.4 / DNV-RP-C212 §5.4', 'Horizontal load against base shear on the effective area. If exceeded the bearing inclination factors lose meaning and the mat translates.'));
    const ok = checks.every((c) => c.uc <= 1);
    const areaNow = B * L, areaReq = req ? req.B * req.L : NaN;
    const dSteel = req ? (areaNow - areaReq) * STEEL_T_PER_M2 : 0;
    const advice = [];
    if (ok) advice.push(adv('engineering', 'ok', `Bearing utilisation ${fmt(uc)} with effective area ${fmt(Aeff)} m² of ${fmt(areaNow)} m² (eccentricity ${fmt(eB)} m).`,
      req ? `The smallest mat of the same proportions that still passes bearing and sliding is ${fmt(req.B)} m × ${fmt(req.L)} m (${fmt(areaReq)} m²). Before trimming, confirm su from site-specific CPT / T-bar data and check settlement and installation-stage loads.` : 'Confirm su from site-specific CPT / T-bar data and check settlement before optimising the footprint.',
      'ISO 19901-4:2016 §7; API RP 2GEO (2011) §7.3; DNV-RP-C212 (2019) §5',
      dSteel > 0 ? `Trimming to the required footprint would save about ${steelCo2Text(dSteel)} It also reduces the disturbed seabed footprint by ${fmt(areaNow - areaReq)} m².` : 'The footprint is already close to the minimum; no material saving is available from bearing.'));
    else advice.push(adv('engineering', 'fail', `Bearing or sliding on the effective area is exceeded (bearing ${fmt(uc)}${classical ? `, sliding ${fmt(classical.ucSliding)}` : ''}).`,
      req ? `Enlarge the mat to at least ${fmt(req.B)} m × ${fmt(req.L)} m, or deepen the skirts to reach stronger soil (su rises ${fmt(inp.suGrad ?? 0)} kPa per metre in clay), or reduce the eccentricity M/V = ${fmt(eB)} m by re-centring equipment.` : 'No mat up to 80 m wide satisfies the loads: reduce the horizontal load and moment, or change to a skirted or piled foundation.',
      'ISO 19901-4:2016 §7; API RP 2GEO (2011) §7.3; DNV-RP-C212 (2019) §5',
      req ? `The enlarged mat adds about ${steelCo2Text(Math.abs(dSteel))} Deeper skirts usually reach the same capacity with less plate.` : 'A foundation concept change carries a large material and installation footprint: resolve at concept stage.'));
    advice.push(adv('economic', ok ? 'info' : 'warn', `Mat steel scales with plan area at roughly ${STEEL_T_PER_M2} t/m² (screening): the present ${fmt(areaNow)} m² is about ${fmt(areaNow * STEEL_T_PER_M2)} t${req ? `, the required footprint ${fmt(areaReq * STEEL_T_PER_M2)} t` : ''}.`,
      'A site-specific geotechnical survey that narrows the su range usually pays for itself: a 20 % higher characteristic su reduces the required area roughly in proportion. Check that the mat still fits the installation vessel deck and crane radius.',
      'Value-of-information reasoning; AACE RP 18R-97 estimate classes; ISO 19901-8 marine soil investigations'));
    const rows = [
      ['Terzaghi (strip, vertical central load)', strip.qult], [`Terzaghi (${Math.abs(L / B - 1) < 0.05 ? 'square' : 'strip'} shape)`, terz.qult],
      ['Meyerhof (1963)', res.meyerhof.qult], ['Brinch Hansen (1970)', res.hansen.qult], ['Vesic (1973)', res.vesic.qult],
      ...(s.clay ? [['Skempton (1951), vertical central load', skQ], ['ISO 19901-4 / API RP 2GEO undrained', res.iso?.qult]] : []),
    ].map(([m, q]) => ({ method: m, qult: Number.isFinite(q) ? round(q / 1e3) : 'n/a' }));
    return {
      kpis: [
        { key: 'ucBearing', label: 'Bearing utilisation', unit: '-', value: uc, limit: 1, sense: 'max', explain: 'Factored vertical load over the design bearing capacity on the effective area. Above 1.0 the soil is predicted to fail in general shear (capped at 99 when the inclined load leaves no capacity).' },
        { key: 'qult', label: 'Ultimate bearing pressure', unit: 'kPa', value: sel.qult / 1e3, explain: 'Gross ultimate pressure on the effective area by the selected method, before the material factor.' },
        { key: 'Vult', label: 'Design vertical capacity', unit: 'kN', value: Vcap / 1e3, explain: 'Capacity available to the factored vertical load: (qult − overburden)/γm plus overburden, times the effective area.' },
        { key: 'qApplied', label: 'Effective contact pressure', unit: 'kPa', value: V / Aeff / 1e3, explain: 'Factored vertical load spread over the effective (Meyerhof) area. It rises quickly as eccentricity shrinks the effective width.' },
        { key: 'Aeff', label: 'Effective bearing area', unit: 'm²', value: Aeff, explain: 'Area of the reduced rectangle on which the eccentric load is centred. A low ratio to the gross area signals a moment-dominated foundation.' },
        { key: 'ecc', label: 'Load eccentricity M/V', unit: 'm', value: eB, limit: B / 6, sense: 'max', explain: 'Offset of the vertical resultant from the mat centre. Beyond B/6 (the kern) the heel of an unskirted mat lifts off.' },
        { key: 'suDesign', label: 'Strength used', unit: 'kPa', value: (s.clay ? su : s.c) / 1e3, explain: 'Undrained strength at B′/4 below skirt-tip level (clay) or the effective cohesion (sand) entering the capacity equation.' },
        { key: 'Nc', label: 'Bearing factor Nc', unit: '-', value: bearingNc(s.phi), explain: 'Prandtl cohesion factor: 5.14 for undrained clay, rising with friction angle.' },
        { key: 'Nq', label: 'Bearing factor Nq', unit: '-', value: bearingNq(s.phi), explain: 'Prandtl–Reissner surcharge factor: 1 for φ = 0, about 18 at 30°, 33 at 35°.' },
        { key: 'Ngamma', label: 'Bearing factor Nγ (Vesic)', unit: '-', value: bearingNgamma(s.phi, 'vesic'), explain: 'Self-weight factor 2(Nq + 1)tan φ; it carries the benefit of mat width in sand and is zero in undrained clay.' },
        { key: 'Breq', label: 'Required mat width (same proportions)', unit: 'm', value: req ? req.B : 99, explain: 'Smallest width for which both bearing and sliding utilisations reach 1.0 at the present L/B. 99 means no solution within 80 m.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Vertical capacity against horizontal load', x: { label: 'Factored horizontal load', unit: 'kN' }, y: { label: 'Design vertical capacity', unit: 'kN' },
          series: [{ name: 'Capacity envelope', x: hs, y: vs }, { name: 'Design point', x: [H / 1e3], y: [V / 1e3], points: true }],
          explain: 'How load inclination erodes the vertical capacity at the present eccentricity. The design point must lie below the curve; its vertical distance to the curve is the reserve.' },
        { type: 'bar', title: 'Ultimate bearing pressure by method', categories: rows.map((r) => r.method), series: [{ name: 'qult (kPa)', y: rows.map((r) => (typeof r.qult === 'number' ? r.qult : 0)) }], y: { label: 'qult', unit: 'kPa' },
          explain: 'Gross ultimate pressure by each author. The spread between methods is a measure of model uncertainty; the Terzaghi and Skempton bars ignore inclination and eccentricity.' },
        { type: 'line', title: 'Bearing utilisation against mat width', x: { label: 'Width B (L scaled in proportion)', unit: 'm' }, y: { label: 'Utilisation', unit: '-' },
          series: [{ name: 'Utilisation', x: bs, y: ucB }, { name: 'Limit', x: [bs[0] ?? 0, bs[bs.length - 1] ?? 1], y: [1, 1], dash: true }],
          explain: 'Utilisation for geometrically similar mats. Read the required width where the curve crosses 1.0.' },
      ],
      tables: [
        { title: 'Method comparison', columns: [{ key: 'method', label: 'Method' }, { key: 'qult', label: 'Ultimate pressure', unit: 'kPa' }], rows,
          explain: 'Gross ultimate bearing pressure by each implemented method for the same soil and geometry; use it to judge the sensitivity of the result to the choice of method.' },
        { title: 'Factors of the selected method', columns: [{ key: 'k', label: 'Factor' }, { key: 'v', label: 'Value' }],
          rows: Object.entries(sel.factors || {}).filter(([, v]) => typeof v === 'number').map(([k, v]) => ({ k, v: round(v) })).concat(dp ? [{ k: 'Drucker–Prager α (compression match)', v: round(dp.alpha) }, { k: 'Drucker–Prager k (kPa)', v: round(dp.k / 1e3) }, { k: 'Mohr–Coulomb triaxial M', v: round(mc.M) }] : []),
          explain: 'Shape (s), depth (d), inclination (i), base-tilt (b) and ground-slope (g) factors applied to each term, plus the equivalent Drucker–Prager and triaxial strength parameters for continuum checks in sand.' },
      ],
      notes: [`Selected method: ${inp.method}; strength taken at ${fmt(zSu)} m below mudline.`, `Effective width B′ = ${fmt(Beff)} m from eccentricity ${fmt(eB)} m.`, uc > 1 ? 'Bearing governs: the factored load exceeds the design capacity.' : `Reserve on bearing: ${fmt((1 / uc - 1) * 100)} %.`],
      advice,
      publish: { 'subsea.bearingUc': uc, 'subsea.Vult': Vcap, 'subsea.matB': B, 'subsea.matL': L },
      quantities: { steel_t: round(areaNow * STEEL_T_PER_M2) },
    };
  },
  verify: [
    { name: 'Square mat on uniform clay, vertical load: qult = (2+π)·su·(1 + (B/L)/Nc) = (2+π)·su + su', source: 'Vesic (1973) shape factor sc = 1 + (B/L)(Nq/Nc) with Nq = 1: qult = 5.1416·10 + 10 = 61.416 kPa',
      inputs: { B: 10, L: 10, D: 0, soil: 'clay', su0: 10, suGrad: 0, V: 1000, H: 0, M: 0, method: 'vesic', gammaF: 1, gammaM: 1 },
      expect: [{ kpi: 'qult', value: 61.416, tol: 1e-4 }, { kpi: 'Vult', value: 6141.6, tol: 1e-4 }, { kpi: 'Nc', value: 5.14159, tol: 1e-5 }] },
    { name: 'Bearing factors at φ′ = 30°: Nq = e^(π tan φ)·tan²(45°+φ/2), Nc = (Nq−1)cot φ, Nγ = 2(Nq+1)tan φ', source: 'Prandtl–Reissner; Vesic (1973) Table: Nq = 18.40, Nc = 30.14, Nγ = 22.40',
      inputs: { soil: 'sand', phiDeg: 30, cEff: 0, H: 0, M: 0, method: 'vesic' },
      expect: [{ kpi: 'Nq', value: 18.40, tol: 0.001 }, { kpi: 'Nc', value: 30.14, tol: 0.001 }, { kpi: 'Ngamma', value: 22.40, tol: 0.001 }] },
  ],
};
