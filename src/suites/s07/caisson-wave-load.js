// Wave pressure on an upright (caisson) breakwater or vertical seawall by Goda's method with the Takahashi
// impulsive coefficient, and the resulting sliding, overturning and heel-pressure margins.
import { godaPressure, caissonStability, sainflou } from '../../engine/coastal/caisson.js';
import { G, CO2_CONC, need, run, finite, r, f, linspace, check, sfCheck, eng, eco, slrInputs, slrRows, SLR_BASIS } from './_common.js';

export default {
  id: 's07-caisson-wave-load', suite: 7, title: 'Wave pressure on a vertical breakwater',
  summary: 'Goda wave pressures, uplift and impulsive breaking load on a caisson, with sliding, overturning and heel-pressure margins and the caisson width they require.',
  covers: [2, 5], fidelity: 'analytical',
  standards: ['Goda, Random Seas and Design of Maritime Structures (2000) §4.2–4.3', 'Technical Standards and Commentaries for Port and Harbour Facilities in Japan (OCDI 2020) Part III ch. 2', 'BS 6349-1-2:2016 / BS 6349-7', 'PROVERBS (Oumeraci et al. 2001) for impulsive loads'],
  models: ['Goda wave-pressure formula', 'Takahashi impulsive pressure coefficient', 'Sainflou standing-wave pressure', 'Hydrostatic pressure and Archimedes buoyancy', 'Sliding resistance (Coulomb friction)', 'Overturning stability', 'Linear base-pressure distribution'],
  limits: 'Quasi-static 2-D pressures on a rigid upright section per metre run; the short-duration dynamics of impulsive breaking (caisson rocking, sliding distance), 3-D wave obliquity along the wall and the bearing capacity of the mound and subsoil under eccentric inclined load are not modelled here. Impulsive conditions flagged by the module need hydraulic model tests; CFD and wave-flume precision is not implied.',
  publishes: ['marine.wallForce', 'marine.wallMoment', 'marine.wallUplift', 'marine.caissonWidth'],
  searchTerms: 'Goda wave pressure caisson breakwater sliding impulsive breaking',
  inputs: [
    { key: 'Hmax', label: 'Design (maximum) wave height at the wall', unit: 'm', type: 'number', default: 7, min: 0.3, max: 25, group: 'Wave boundary', role: 'bc', bind: 'marine.HmaxToe', uncertain: { dist: 'lognormal', cov: 0.1 }, help: '1.8 H1/3 seaward of the surf zone, or the surf-zone Hmax 5 H1/3 seaward of the wall. Before the sea-level-rise allowance.' },
    { key: 'T', label: 'Design wave period T1/3', unit: 's', type: 'number', default: 10, min: 3, max: 22, group: 'Wave boundary', role: 'bc', uncertain: { dist: 'normal', cov: 0.08 } },
    { key: 'betaDeg', label: 'Wave angle from the wall normal', unit: 'deg', type: 'number', default: 0, min: 0, max: 75, group: 'Wave boundary', role: 'bc', help: 'Goda recommends rotating the principal direction up to 15° towards the normal.' },
    { key: 'h', label: 'Water depth in front of the breakwater', unit: 'm', type: 'number', default: 12, min: 1, max: 60, group: 'Water level', role: 'bc', bind: 'marine.hToe', help: 'At the design still-water level before the sea-level-rise allowance.' },
    { key: 'hb', label: 'Depth 5 H1/3 seaward of the wall', unit: 'm', type: 'number', default: 13, min: 1, max: 80, group: 'Water level', role: 'bc', help: 'Deeper than the wall depth on a sloping foreshore; drives the breaking-pressure coefficient α2.' },
    { key: 'depthLimited', label: 'Design wave is depth-limited', type: 'bool', default: true, group: 'Water level', role: 'bc' },
    ...slrInputs(),
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1030, min: 990, max: 1050, group: 'Water level', role: 'environment', bind: 'site.rhoW' },
    { key: 'd', label: 'Depth above the mound armour', unit: 'm', type: 'number', default: 9, min: 0.5, max: 60, group: 'Section', role: 'geometry', design: true, help: 'A high mound (small d) promotes impulsive breaking.' },
    { key: 'hP', label: 'Depth of the caisson base below still water h′', unit: 'm', type: 'number', default: 10, min: 0.5, max: 60, group: 'Section', role: 'geometry' },
    { key: 'hc', label: 'Crest level above still water', unit: 'm', type: 'number', default: 4, min: 0, max: 20, group: 'Section', role: 'geometry', design: true },
    { key: 'B', label: 'Caisson width', unit: 'm', type: 'number', default: 15, min: 1, max: 60, group: 'Section', role: 'geometry', design: true },
    { key: 'BM', label: 'Mound berm width in front of the caisson', unit: 'm', type: 'number', default: 6, min: 0, max: 60, group: 'Section', role: 'geometry', design: true, help: 'Used for the impulsive coefficient αI; 0 ignores it.' },
    { key: 'rhoC', label: 'Average density of the filled caisson', unit: 'kg/m³', type: 'number', default: 2100, min: 1300, max: 2600, group: 'Section', role: 'material', uncertain: { dist: 'normal', cov: 0.03 }, help: 'Concrete shell plus sand or water ballast, averaged over the gross volume.' },
    { key: 'concFrac', label: 'Concrete share of the caisson volume', unit: '-', type: 'number', default: 0.28, min: 0.1, max: 1, group: 'Section', role: 'material' },
    { key: 'mu', label: 'Friction coefficient caisson / rubble mound', unit: '-', type: 'number', default: 0.6, min: 0.3, max: 0.9, group: 'Foundation boundary', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 }, help: '0.6 for concrete on rubble (Goda); up to 0.7–0.8 with friction mats.' },
    { key: 'pAllow', label: 'Allowable heel pressure on the mound', unit: 'kPa', type: 'number', default: 500, min: 100, max: 1500, group: 'Foundation boundary', role: 'criteria', help: 'Goda quotes 400–500 kPa on a rubble mound.' },
    { key: 'sfS', label: 'Required sliding safety factor', unit: '-', type: 'number', default: 1.2, min: 1, max: 2, group: 'Criteria', role: 'criteria' },
    { key: 'sfO', label: 'Required overturning safety factor', unit: '-', type: 'number', default: 1.2, min: 1, max: 2.5, group: 'Criteria', role: 'criteria' },
    { key: 'length', label: 'Length of breakwater', unit: 'm', type: 'number', default: 500, min: 1, max: 20000, group: 'Section', role: 'geometry' },
  ],
  compute(inp) {
    need(inp.d <= inp.hP + 1e-9, 'The depth above the mound armour cannot exceed the depth of the caisson base (the armour sits on the mound in front of the base).');
    need(inp.hP <= inp.h + 1e-9, 'The caisson base cannot be deeper than the seabed in front of it: h′ must not exceed h.');
    need(inp.hb >= inp.h - 1e-9, 'The depth 5 H1/3 seaward of the wall should not be less than the depth at the wall.');

    /** Loads and margins for a sea-level rise and a caisson width. */
    const state = (rise, B) => {
      const k = inp.depthLimited ? (inp.h + rise) / inp.h : 1;
      const hc = Math.max(inp.hc - rise, 0), hP = inp.hP + rise;
      const g = godaPressure({ Hmax: inp.Hmax * k, T: inp.T, h: inp.h + rise, d: inp.d + rise, hPrime: hP, hc, hb: inp.hb + rise, B, betaDeg: inp.betaDeg, BM: inp.BM > 0 ? inp.BM : undefined, rho: inp.rhoW });
      // weight in still water: gross weight less buoyancy of the submerged part
      const W = inp.rhoC * G * B * (inp.hP + inp.hc) - inp.rhoW * G * B * hP;
      need(W > g.U, `The caisson (average density ${f(inp.rhoC)} kg/m³) is lifted by buoyancy and wave uplift at this water level: increase the fill density or the crest height.`);
      const st = caissonStability({ P: g.P, Mp: g.Mp, U: g.U, Mu: g.Mu, W, t: B / 2, B, mu: inp.mu });
      return { g, st, W, hc, hP, H: inp.Hmax * k };
    };
    const s = run(() => state(inp.slr, inp.B), 'Goda wave pressure');
    const { g, st } = s;
    finite(g.P, 'Wave force');
    const margin = (B) => { try { const q = state(inp.slr, B); return Math.min(q.st.sfSliding / inp.sfS, q.st.sfOverturning / inp.sfO); } catch { return 0; } };
    let lo = 1, hi = 80;
    for (let i = 0; i < 50; i++) { const mid = 0.5 * (lo + hi); if (margin(mid) < 1) lo = mid; else hi = mid; }
    const Breq = hi;
    let sain = null; try { sain = sainflou({ H: s.H / 1.8, T: inp.T, h: inp.h + inp.slr, hc: s.hc, rho: inp.rhoW }); } catch { sain = null; }
    const impulsive = g.alphaI > g.alpha2 && g.alphaStar > 0.3;
    const pHeel = Number.isFinite(st.pHeel) ? st.pHeel : 1e9;

    const scen = slrRows(inp).map((sc) => { try { const q = state(sc.rise, inp.B); return { name: sc.name, rise: sc.rise, H: r(q.H, 3), P: r(q.g.P / 1e3), U: r(q.g.U / 1e3), sfS: r(q.st.sfSliding, 3), sfO: r(q.st.sfOverturning, 3), pHeel: r(q.st.pHeel / 1e3) }; } catch { return { name: sc.name, rise: sc.rise, H: 0, P: 0, U: 0, sfS: 0, sfO: 0, pHeel: 0 }; } });
    const bs = linspace(Math.max(2, 0.5 * Math.min(inp.B, Breq)), 1.6 * Math.max(inp.B, Breq), 30);
    const vol = inp.B * (inp.hP + inp.hc), conc = vol * inp.concFrac * inp.length;
    const dB = Breq - inp.B, dConc = dB * (inp.hP + inp.hc) * inp.concFrac * inp.length;
    const ok = st.sfSliding >= inp.sfS && st.sfOverturning >= inp.sfO && pHeel / 1e3 <= inp.pAllow;
    const gov = st.sfSliding / inp.sfS < st.sfOverturning / inp.sfO ? 'sliding' : 'overturning';

    const advice = [];
    advice.push(ok
      ? eng('ok', `The ${f(inp.B)} m caisson resists P = ${f(g.P / 1e3)} kN/m and U = ${f(g.U / 1e3)} kN/m with sliding SF ${f(st.sfSliding)} and overturning SF ${f(st.sfOverturning)}; heel pressure ${f(pHeel / 1e3)} kPa.`,
        `The required width is ${f(Breq)} m (${gov} governs), so ${f(-dB, 2)} m is reserve. Confirm the mound and subsoil bearing under the eccentric inclined load and check the higher sea-level scenarios before narrowing.`,
        'Goda (2000) §4.3; OCDI (2020) Part III ch. 2 (confirm current edition)',
        `Narrowing to the required width would save about ${f(Math.max(0, -dConc), 3)} m³ of concrete (≈ ${f(Math.max(0, -dConc) * CO2_CONC, 3)} t CO₂e, screening factor) but removes the sea-level-rise reserve.`)
      : eng('fail', `The ${f(inp.B)} m caisson does not meet the margins: sliding SF ${f(st.sfSliding)} (need ${f(inp.sfS)}), overturning SF ${f(st.sfOverturning)} (need ${f(inp.sfO)}), heel pressure ${f(pHeel / 1e3)} kPa (allow ${f(inp.pAllow)}).`,
        `Widen the caisson to about ${f(Breq)} m, or raise the fill density, add a friction mat (μ ≈ 0.7–0.8), or lower the crest and accept more overtopping; widening also spreads the heel pressure.`,
        'Goda (2000) §4.3; OCDI (2020) Part III ch. 2 (confirm current edition)',
        `Widening by ${f(Math.max(dB, 0), 2)} m adds about ${f(Math.max(dConc, 0), 3)} m³ of concrete over ${f(inp.length)} m (≈ ${f(Math.max(dConc, 0) * CO2_CONC, 3)} t CO₂e, screening factor); heavier fill or a friction mat reaches the same margin with almost no extra concrete.`));
    if (impulsive) advice.push(eng('warn', `Impulsive breaking is likely: αI = ${f(g.alphaI, 2)} exceeds α2 = ${f(g.alpha2, 2)} (mound depth d/h = ${f((inp.d + inp.slr) / (inp.h + inp.slr), 2)}, berm ${f(inp.BM)} m).`,
      'Lower or narrow the mound berm, or use a sloping-top or perforated caisson; confirm the pressures in a hydraulic model, because the impulsive load can be several times the Goda value for a few hundredths of a second.', 'Takahashi, Tanimoto & Shimosako (1994); PROVERBS (2001) parameter map'));
    const worst = scen[scen.length - 1];
    advice.push(eng(worst.sfS < inp.sfS ? 'warn' : 'info', `Sea-level rise: at +${f(worst.rise)} m (${worst.name}) the wave force becomes ${f(worst.P)} kN/m and the sliding SF ${f(worst.sfS)} (overturning ${f(worst.sfO)}).`,
      worst.sfS < inp.sfS ? 'Size the width for the adopted end-of-life scenario, or keep space on the harbour side for a rubble backing berm that can be added when the trigger level is reached.' : 'The section keeps its margins through the scenarios; overtopping, not stability, will set the adaptation date.', SLR_BASIS));
    advice.push(eco('info', `The section holds ${f(vol, 3)} m³ gross per metre; over ${f(inp.length)} m that is about ${f(conc, 3)} m³ of structural concrete plus ${f(vol * (1 - inp.concFrac) * inp.length, 3)} m³ of fill.`,
      `Each metre of caisson width costs about ${f((inp.hP + inp.hc) * inp.concFrac * inp.length, 3)} m³ of concrete over the length; compare with a higher-friction base or denser fill, which buy sliding margin without changing the casting yard and float-out draught.`,
      'AACE RP 18R-97 quantity-based estimating; whole-life costing per ISO 15686-5',
      `The concrete in the caissons embodies about ${f(conc * CO2_CONC, 3)} t CO₂e (screening factor 0.3 t/m³); low-clinker cements suit the massive, low-strength-demand sections.`));

    const zs = [-s.hP, 0, s.g.hcStar], ps = [g.p3, g.p1, g.p4];
    return {
      kpis: [
        { key: 'P', label: 'Horizontal wave force', value: r(g.P / 1e3), unit: 'kN/m', explain: 'Resultant of the Goda pressure diagram on the front face per metre of wall, under the wave crest. The driving load for sliding and overturning.' },
        { key: 'U', label: 'Wave uplift on the base', value: r(g.U / 1e3), unit: 'kN/m', explain: 'Triangular wave-induced uplift under the caisson, additional to buoyancy. It reduces the weight available for friction.' },
        { key: 'p1', label: 'Peak pressure at still-water level', value: r(g.p1 / 1e3), unit: 'kPa', explain: 'Largest ordinate of the Goda diagram. It sizes the front wall panels and shows the breaking contribution when α* is high.' },
        { key: 'Mp', label: 'Overturning moment of the wave force', value: r(g.Mp / 1e3), unit: 'kN·m/m', explain: 'Moment of the horizontal wave force about the caisson base; with the uplift moment it tries to rotate the caisson about its harbour-side heel.' },
        { key: 'etaStar', label: 'Height to which pressure acts η*', value: r(g.etaStar), unit: 'm', explain: 'Elevation above still water where the wave pressure falls to zero. A crest below η* is overtopped and the diagram is truncated there.' },
        { key: 'alphaStar', label: 'Breaking-pressure coefficient α*', value: r(g.alphaStar, 3), unit: '-', explain: 'Larger of Goda\'s α2 and the Takahashi impulsive coefficient αI. Values above about 0.3 mean breaking-wave impact dominates and model tests are needed.' },
        { key: 'sfSliding', label: 'Sliding safety factor', value: r(Math.min(st.sfSliding, 999), 4), unit: '-', limit: inp.sfS, sense: 'min', explain: 'Friction resistance μ(W − U) over the horizontal wave force. Below 1.0 the caisson slides on the mound; Japanese practice asks for 1.2.' },
        { key: 'sfOverturning', label: 'Overturning safety factor', value: r(Math.min(st.sfOverturning, 999), 4), unit: '-', limit: inp.sfO, sense: 'min', explain: 'Restoring moment of the buoyant weight over the moments of wave force and uplift about the heel. Below 1.0 the caisson rotates harbourwards.' },
        { key: 'pHeel', label: 'Heel bearing pressure', value: r(Math.min(pHeel, 1e8) / 1e3), unit: 'kPa', limit: inp.pAllow, sense: 'max', explain: 'Largest contact pressure under the harbour-side edge from the eccentric resultant. High values crush or punch the mound and tilt the caisson.' },
        { key: 'ecc', label: 'Eccentricity of the base resultant', value: r(st.eccentricity, 3), unit: 'm', explain: 'Offset of the resultant from mid-base. Beyond B/6 the seaward edge lifts off and the contact width shrinks.' },
        { key: 'Breq', label: 'Caisson width for the required margins', value: r(Breq), unit: 'm', explain: 'Smallest width that meets both the sliding and overturning factors in the design case. Compare with the width provided.' },
        { key: 'Psain', label: 'Sainflou standing-wave force (H1/3)', value: r(sain ? sain.crest.force / 1e3 : 0), unit: 'kN/m', explain: 'Total crest-phase force (hydrostatic included) for non-breaking standing waves, given for comparison in deep, non-breaking conditions. Zero when outside the formula range.' },
      ],
      checks: [
        sfCheck('Sliding on the rubble mound', st.sfSliding, inp.sfS, 'Goda (2000) Eq. 4.17; OCDI (2020)', 'Guards against the caisson being pushed harbourwards along the mound surface. Sliding is progressive: each large wave moves the caisson a little further.'),
        sfCheck('Overturning about the heel', st.sfOverturning, inp.sfO, 'Goda (2000) Eq. 4.18; OCDI (2020)', 'Guards against rotation about the harbour-side edge. In practice the mound fails in bearing under the heel before pure overturning.'),
        check('Heel pressure on the mound', pHeel / 1e3, inp.pAllow, 'kPa', 'Goda (2000) §4.3 (400–500 kPa on rubble)', 'Guards against crushing and punching of the rubble mound under the heel, which tilts the caisson and opens the joints.'),
      ],
      plots: [
        { type: 'line', title: 'Wave pressure on the front face', x: { label: 'Pressure', unit: 'kPa' }, y: { label: 'Elevation to still water', unit: 'm' },
          series: [{ name: 'Goda', x: [0, ...ps.map((p) => r(p / 1e3, 3)), 0], y: [zs[0], ...zs.map((z) => r(z, 3)), zs[2]] }],
          explain: 'Goda pressure diagram from the caisson base to the crest: peak at still-water level, reducing with depth and truncated at the crest. The enclosed area is the horizontal force.' },
        { type: 'line', title: 'Uplift under the base', x: { label: 'Distance from seaward edge', unit: 'm' }, y: { label: 'Wave uplift pressure', unit: 'kPa' },
          series: [{ name: 'Uplift', x: [0, inp.B], y: [r(g.pu / 1e3, 3), 0] }],
          explain: 'Wave-induced uplift, largest at the seaward toe and zero at the heel, additional to buoyancy. A base drain or skirt changes this boundary condition.' },
        { type: 'line', title: 'Safety factors versus caisson width', x: { label: 'Width B', unit: 'm' }, y: { label: 'Safety factor', unit: '-' },
          series: [{ name: 'Sliding', x: bs.map((b) => r(b, 3)), y: bs.map((b) => { try { return r(Math.min(state(inp.slr, b).st.sfSliding, 10), 3); } catch { return 0; } }) }, { name: 'Overturning', x: bs.map((b) => r(b, 3)), y: bs.map((b) => { try { return r(Math.min(state(inp.slr, b).st.sfOverturning, 10), 3); } catch { return 0; } }) }, { name: 'Required (sliding)', x: [bs[0], bs[29]], y: [inp.sfS, inp.sfS], dash: true }],
          explain: 'How the margins grow with width. The width where the lower curve crosses the required line is the minimum caisson; the slope shows what each extra metre buys.' },
      ],
      tables: [
        { title: 'Goda pressures and coefficients', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['α1 (period)', g.alpha1, '-'], ['α2 (breaking)', g.alpha2, '-'], ['αI (impulsive)', g.alphaI, '-'], ['α3 (depth)', g.alpha3, '-'], ['p1 at still water', g.p1 / 1e3, 'kPa'], ['p2 at seabed', g.p2 / 1e3, 'kPa'], ['p3 at caisson base', g.p3 / 1e3, 'kPa'], ['p4 at crest', g.p4 / 1e3, 'kPa'], ['pu uplift at toe', g.pu / 1e3, 'kPa'], ['Wavelength L', g.L, 'm'], ['Buoyant weight W', s.W / 1e3, 'kN/m'], ['Contact width', st.contactWidth, 'm']].map(([q, v, u]) => ({ q, v: r(v), u })),
          explain: 'The intermediate values of Goda\'s method for the design case, for checking against a hand calculation or a model-test report.' },
        { title: 'Sea-level-rise scenarios', columns: [{ key: 'name', label: 'Scenario' }, { key: 'rise', label: 'Rise', unit: 'm' }, { key: 'H', label: 'Hmax', unit: 'm' }, { key: 'P', label: 'Force', unit: 'kN/m' }, { key: 'U', label: 'Uplift', unit: 'kN/m' }, { key: 'sfS', label: 'Sliding SF', unit: '-' }, { key: 'sfO', label: 'Overturning SF', unit: '-' }, { key: 'pHeel', label: 'Heel pressure', unit: 'kPa' }],
          rows: scen, explain: 'The present section under each sea-level-rise scenario: deeper water admits a larger wave, adds buoyancy and lowers the freeboard. Zero rows mean the section floats or is outside the formula range.' },
      ],
      notes: [
        `${gov[0].toUpperCase()}${gov.slice(1)} governs. Design wave ${f(s.H)} m with ${f(inp.slr)} m of sea-level rise; buoyant weight ${f(s.W / 1e3)} kN/m.`,
        impulsive ? 'Impulsive breaking flagged: treat the result as a lower bound pending model tests.' : 'Pressures are in the pulsating (non-impulsive) range of the method.',
      ],
      advice,
      publish: { 'marine.wallForce': g.P, 'marine.wallMoment': g.Mp, 'marine.wallUplift': g.U, 'marine.caissonWidth': inp.B },
      quantities: { concrete_m3: r(conc) },
    };
  },
  verify: [
    { name: 'Goda pressure formula, hand evaluation',
      source: 'Goda (2000) §4.2 Eq. 4.2–4.11. Hmax = 7 m, T = 10 s, h = 12, d = 9, h′ = 10, hc = 4, hb = 13 m, ρ = 1030: L = 99.71 m, α1 = 0.8455, α2 = 0.0620, α3 = 0.8078, η* = 10.5 m, p1 = 64.17, p3 = 51.83, p4 = 39.72 kPa → P = ½(p1+p3)h′ + ½(p1+p4)hc = 787.8 kN/m; pu = 48.29 kPa, U = ½ pu B = 362.2 kN/m; W = (2100 × 14 − 1030 × 10) × 9.80665 × 15 = 2809.6 kN/m → sliding SF = 0.6 (W − U)/P = 1.864',
      inputs: { Hmax: 7, T: 10, h: 12, d: 9, hP: 10, hc: 4, hb: 13, B: 15, BM: 0, betaDeg: 0, slr: 0, rhoW: 1030, rhoC: 2100, mu: 0.6 },
      expect: [{ kpi: 'P', value: 787.8, tol: 0.002 }, { kpi: 'U', value: 362.2, tol: 0.002 }, { kpi: 'sfSliding', value: 1.864, tol: 0.003 }] },
  ],
};
