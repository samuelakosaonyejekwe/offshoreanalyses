// Deflection and settlement of a wharf in service: deck sway, cracked and long-term (creep and shrinkage) deflection
// of the deck beam, load–settlement of the piles and the differential settlement between rows, and the
// consolidation settlement of the backland with the step that develops at the wharf–backland transition.
import { momentCapacity, crackingMoment, crackedSection, crackedDeflection, effectiveModulus } from '../../engine/structural/concrete.js';
import { axialLoadSettlement } from '../../engine/geotech/pile.js';
import { consolidationSettlement, timeSettlementCurve, consolidationU, consolidationFD, secondaryCompression } from '../../engine/geotech/settlement.js';
import { YEAR_S, need, run, r, f, check, eng, eco, linspace, interp, bentInputs, bentLoadInputs, pileInputs, rcInputs, rcMember, soilInputs, profileFrom, bentModel, solveBent, framePlot, deformScale, concCo2, CONFIRM, NOT_HIFI, RC_NOTE, ECON_BASIS } from './_common.js';

const CLAY_DEFAULT = [{ name: 'Soft marine clay', thick: 5, gamma: 6.5, Cc: 0.45, Cr: 0.06, e0: 1.4, ocr: 1.2 }, { name: 'Firm clay', thick: 6, gamma: 8.5, Cc: 0.25, Cr: 0.04, e0: 0.9, ocr: 2.0 }];

/**
 * Notional creep coefficient φ(∞, t0) of EN 1992-1-1 Annex B (Eq. B.1–B.5, B.8c):
 *   φ0 = φRH·β(fcm)·β(t0);  φRH = [1 + (1 − RH/100)/(0.1·h0^(1/3))·α1]·α2 (α1 = α2 = 1 for fcm ≤ 35 MPa);
 *   β(fcm) = 16.8/√fcm;  β(t0) = 1/(0.1 + t0^0.2);  α1 = (35/fcm)^0.7, α2 = (35/fcm)^0.2.
 * fcm [MPa], RH [%], h0 = 2Ac/u [mm], t0 [days].
 */
export function creepCoefficient({ fcm, RH, h0, t0 }) {
  const a1 = fcm > 35 ? (35 / fcm) ** 0.7 : 1, a2 = fcm > 35 ? (35 / fcm) ** 0.2 : 1;
  const phiRH = (1 + ((1 - RH / 100) / (0.1 * Math.cbrt(h0))) * a1) * a2;
  return phiRH * (16.8 / Math.sqrt(fcm)) * (1 / (0.1 + t0 ** 0.2));
}

export default {
  id: 's05-deflection-settlement', suite: 5, title: 'Structural deflection and settlement',
  summary: 'Deck sway and long-term deck beam deflection with cracking, creep and shrinkage, settlement of the piles and its difference between rows, and the consolidation settlement of the backland against the tolerances for cranes, joints and paving.',
  covers: [10], fidelity: 'numerical',
  standards: ['EN 1992-1-1:2004 §7.4, Annex B (creep and shrinkage)', 'EN 1997-1:2004 §2.4.8, Annex H (limiting movements)', 'API RP 2GEO (2011) §8.4 (t–z, Q–z)', 'BS 6349-2:2019 §5 (serviceability)', 'Terzaghi (1943) one-dimensional consolidation'],
  models: ['Elastic frame analysis (matrix stiffness method)', 'Cracked-section analysis with tension stiffening (ζ interpolation)', 'Time-dependent creep (effective modulus) and shrinkage curvature', 'Axial pile load transfer (t–z, Q–z)', 'Elastic settlement of the pile shaft', 'Primary consolidation (e–log σ′ model, Boussinesq stress)', 'Terzaghi consolidation in time (closed form and finite differences)', 'Secondary compression'],
  limits: `Deck beam deflection is the EN 1992 interpolation between uncracked and cracked curvature at the critical section with a deflection coefficient for the span; pile settlement is for single piles (no group settlement interaction); backland settlement is one-dimensional consolidation under a wide fill. Modified Cam-Clay or other coupled finite-element settlement analysis is outside the module. ${RC_NOTE} ${NOT_HIFI}`,
  publishes: ['wharf.deckDeflection', 'wharf.pileSettlement', 'wharf.backlandSettlement', 'wharf.differentialSettlement'],
  searchTerms: 'wharf deck deflection creep pile settlement consolidation backland differential settlement crane rail',
  inputs: [
    ...bentInputs(), ...pileInputs(), ...bentLoadInputs(), ...rcInputs(),
    { key: 'psi2', label: 'ψ2 quasi-permanent share of cargo and crane loads', unit: '-', type: 'number', default: 0.6, min: 0, max: 1, group: 'Serviceability', role: 'criteria', help: `Long-term share of the variable loads; 0.8 for storage areas in EN 1990, less where stacks turn over. ${CONFIRM}` },
    { key: 'kDefl', label: 'Deflection coefficient of the span', unit: '-', type: 'number', default: 0.083, min: 0.03, max: 0.35, group: 'Serviceability', role: 'criteria', help: 'δ = K·L²·curvature at the critical section: 0.104 simply supported with uniform load, about 0.083 for an end span and 0.063 for an interior span of a continuous beam, 0.25 for a cantilever.' },
    { key: 'RH', label: 'Relative humidity', unit: '%', type: 'number', default: 80, min: 40, max: 100, group: 'Serviceability', role: 'environment' },
    { key: 't0', label: 'Age of the concrete at loading', unit: 'days', type: 'number', default: 28, min: 3, max: 365, group: 'Serviceability', role: 'ic' },
    { key: 'phiUser', label: 'Creep coefficient (0 = from EN 1992-1-1 Annex B)', unit: '-', type: 'number', default: 0, min: 0, max: 5, group: 'Serviceability', role: 'material' },
    { key: 'epsCs', label: 'Free shrinkage strain', unit: '×10⁻⁶', type: 'number', default: 300, min: 0, max: 800, group: 'Serviceability', role: 'material', help: 'Total drying plus autogenous shrinkage at the end of the design life; 250–450 ×10⁻⁶ for marine structures.' },
    { key: 'spanLim', label: 'Deck beam deflection limit (span / value)', unit: '-', type: 'number', default: 250, min: 100, max: 1000, group: 'Serviceability', role: 'criteria', help: `EN 1992-1-1 §7.4.1: span/250 for appearance and function under quasi-permanent loads; crane rails may need span/500 or better. ${CONFIRM}` },
    { key: 'swayLim', label: 'Deck sway limit (height / value)', unit: '-', type: 'number', default: 300, min: 50, max: 1000, group: 'Serviceability', role: 'criteria', help: CONFIRM },
    ...soilInputs(),
    { key: 'distLim', label: 'Angular distortion limit between pile rows (1 / value)', unit: '-', type: 'number', default: 500, min: 100, max: 3000, group: 'Soil', role: 'criteria', help: `Differential settlement over row spacing. 1/500 is a common limit for framed structures (EN 1997-1 Annex H); crane rail gauge and level tolerances may be tighter. ${CONFIRM}` },
    { key: 'nElem', label: 'Pile elements for the t–z analysis', unit: '-', type: 'int', default: 40, min: 10, max: 320, group: 'Numerics', role: 'control' },
    { key: 'fillH', label: 'Height of new fill behind the wharf', unit: 'm', type: 'number', default: 3, min: 0, max: 20, group: 'Backland', role: 'bc' },
    { key: 'fillGamma', label: 'Unit weight of the fill', unit: 'kN/m³', type: 'number', default: 19, min: 12, max: 23, group: 'Backland', role: 'material' },
    { key: 'qYard', label: 'Sustained yard surcharge', unit: 'kPa', type: 'number', default: 30, min: 0, max: 200, group: 'Backland', role: 'bc', help: 'Long-term average container or cargo load on the backland.' },
    { key: 'clay', label: 'Compressible layers below the fill', type: 'table', group: 'Backland', role: 'soil',
      columns: [{ key: 'name', label: 'Layer', type: 'text' }, { key: 'thick', label: 'Thickness', unit: 'm', type: 'number' }, { key: 'gamma', label: 'γ′', unit: 'kN/m³', type: 'number' }, { key: 'Cc', label: 'Cc', type: 'number' }, { key: 'Cr', label: 'Cr', type: 'number' }, { key: 'e0', label: 'e0', type: 'number' }, { key: 'ocr', label: 'OCR', type: 'number' }],
      default: CLAY_DEFAULT, help: 'From oedometer tests on borehole samples: compression and recompression indices, initial void ratio and overconsolidation ratio.' },
    { key: 'cv', label: 'Coefficient of consolidation', unit: 'm²/yr', type: 'number', default: 20, min: 0.05, max: 500, group: 'Backland', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.5 }, help: 'Equivalent value for the drainage provided: 1–5 m²/yr for intact clay (field values are often several times the laboratory value because of silt and sand partings), 10–50 m²/yr with prefabricated vertical drains at 1–2 m spacing (default).' },
    { key: 'drain', label: 'Drainage of the clay', type: 'select', default: 'double', group: 'Backland', role: 'bc', options: [{ value: 'double', label: 'Top and bottom (sand below)' }, { value: 'single', label: 'Top only' }] },
    { key: 'Calpha', label: 'Secondary compression index Cα', unit: '-', type: 'number', default: 0.01, min: 0, max: 0.1, group: 'Backland', role: 'soil' },
    { key: 'tHand', label: 'Time from filling to handover of the paving', unit: 'months', type: 'number', default: 18, min: 0, max: 240, group: 'Backland', role: 'ic', design: true, help: 'Surcharging and waiting time before the yard surface and crane beams are finished.' },
    { key: 'life', label: 'Design life', unit: 'years', type: 'number', default: 50, min: 5, max: 120, group: 'Backland', role: 'criteria', bind: 'life.designYears' },
    { key: 'sLim', label: 'Tolerable settlement after handover', unit: 'mm', type: 'number', default: 150, min: 10, max: 1000, group: 'Backland', role: 'criteria', help: 'Step at the wharf–backland transition and regrading tolerance of the yard; 100–300 mm for flexible paving with transition slabs, far less for rail-mounted cranes on grade.' },
  ],
  mesh: { key: 'nElem', label: 'Pile elements (t–z)', levels: [10, 20, 40, 80], qoi: 'sPile' },
  compute(inp, ctx) {
    // --- deck: sway and beam deflection
    const B = bentModel(inp, ctx), S = solveBent(B, inp);
    const qp = S.combo({ G: 1, Q: inp.psi2, C: inp.psi2 }), ch = S.combo({ G: 1, Q: 1, C: 1 }), berth = S.combo({ G: 1, Q: inp.psi2, C: inp.psi2, B: 1 });
    const rc = rcMember(inp, ctx, inp.capB, inp.capH), Ecm = rc.props.Ecm, Es = rc.Es;
    const h0 = 2 * inp.capB * inp.capH / (2 * (inp.capB + inp.capH)) * 1e3;
    const phi = inp.phiUser > 0 ? inp.phiUser : creepCoefficient({ fcm: rc.props.fcm / 1e6, RH: inp.RH, h0, t0: inp.t0 });
    const Eeff = effectiveModulus({ Ecm, phi });
    const span = inp.sx, Mqp = Math.max(qp.Msag, 1), Mch = Math.max(ch.Msag, 1);
    const defl = (M, Ec, beta) => {
      const un = crackingMoment({ section: rc.sag, fctm: rc.props.fctm, Ec, Es }), cr = crackedSection({ section: rc.sag, Ec, Es, M });
      const d = crackedDeflection({ M, Mcr: un.Mcr, EIuncracked: Ec * un.inertia, EIcracked: Ec * cr.inertia, beta, L: span, K: inp.kDefl });
      return { ...d, Mcr: un.Mcr, un, cr, Ec };
    };
    const dShort = defl(Mch, Ecm, 1), dLong = defl(Mqp, Eeff, 0.5);
    // shrinkage curvature 1/r = εcs·αe·S/I (EN 1992-1-1 Eq. 7.21), interpolated between uncracked and cracked sections
    const ae = Es / Eeff, eps = inp.epsCs * 1e-6;
    const kUn = eps * ae * (rc.AsB * (rc.d - dLong.un.centroid) - rc.AsT * (dLong.un.centroid - (inp.capH - rc.dHog))) / dLong.un.inertia;
    const kCr = eps * ae * (rc.AsB * (rc.d - dLong.cr.x)) / dLong.cr.inertia;
    const dShr = inp.kDefl * 1.25 * span * span * Math.max(dLong.zeta * kCr + (1 - dLong.zeta) * kUn, 0); // uniform curvature: K = 1/8 for a simple span, scaled with the span coefficient
    const dTotal = dLong.deflection + dShr, dAllow = span / inp.spanLim;
    const hFix = Math.max(...qp.piles.map((p) => p.L)), swayAllow = hFix / inp.swayLim, sway = Math.max(berth.sway, ch.sway);
    // --- piles: settlement of the most and least loaded rows under the characteristic combination
    const profile = profileFrom(inp, B.P.D);
    const ax = run(() => axialLoadSettlement({ profile, L: inp.Lemb, D: B.P.D, t: B.P.t, E: B.P.E, nElem: Math.round(inp.nElem) }, { steps: 40 }), 'Axial pile response');
    const vert = ch.piles.filter((p) => p.group !== 'raker');
    const settleOf = (p) => { const N = Math.max(-p.N, 0); need(N < ax.capacity, `Pile ${p.id} carries ${f(N / 1e3)} kN under service loads, above its ultimate axial capacity of ${f(ax.capacity / 1e3)} kN: increase the penetration.`); return interp(ax.Q, ax.w, N) + N * p.L / (B.P.E * B.P.A); };
    const sRows = vert.map((p) => ({ id: p.id, N: Math.max(-p.N, 0), s: settleOf(p), x: p.x }));
    const sMax = Math.max(...sRows.map((q) => q.s));
    let dDiff = 0, dist = 0; for (let i = 1; i < sRows.length; i++) { const dd = Math.abs(sRows[i].s - sRows[i - 1].s); if (dd > dDiff) dDiff = dd; dist = Math.max(dist, dd / Math.max(Math.abs(sRows[i].x - sRows[i - 1].x), 0.1)); }
    // --- backland: primary consolidation, its development in time and secondary compression
    const q = inp.fillH * inp.fillGamma * 1e3 + inp.qYard * 1e3;
    const layers = (inp.clay || []).map((c, i) => { need(c.thick > 0 && c.gamma > 0 && c.Cc > 0 && c.Cr >= 0 && c.e0 > 0 && c.ocr >= 1, `Clay layer ${i + 1}: thickness, γ′, Cc and e0 must be positive, Cr not negative and OCR at least 1.`); return { thickness: Number(c.thick), gamma: c.gamma * 1e3, Cc: Number(c.Cc), Cr: Number(c.Cr), e0: Number(c.e0), ocr: Number(c.ocr) }; });
    const Hc = layers.reduce((s, l) => s + l.thickness, 0);
    let sc = 0, sSec = 0, tp = 0, Uh = 1, curve = { t: [0], s: [0] }, fdU = 1;
    const tH = inp.tHand / 12 * YEAR_S, tL = inp.life * YEAR_S, cv = inp.cv / YEAR_S, Hdr = inp.drain === 'double' ? Hc / 2 : Hc;
    if (layers.length && q > 0) {
      // wide fill: the stress increment equals the surcharge over the full clay thickness (B much larger than the depth)
      sc = run(() => consolidationSettlement({ q, B: 2000, L: 2000, sigmaV0Base: 1e3, layers, nSub: 20 }), 'Consolidation settlement').s;
      const e0m = layers.reduce((s, l) => s + l.e0 * l.thickness, 0) / Hc;
      const times = linspace(0, tL, 60).map((t) => t).concat([tH]).sort((a, b) => a - b);
      const ts = run(() => timeSettlementCurve({ times: times.map((t) => Math.max(t, 1)), sc, cv, Hdr, Calpha: inp.Calpha, e0: e0m, H: Hc }), 'Time–settlement curve');
      curve = { t: times, s: ts.s }; tp = ts.tp;
      Uh = consolidationU(cv * tH / (Hdr * Hdr));
      sSec = tL > tp && inp.Calpha > 0 ? secondaryCompression({ Calpha: inp.Calpha, e0: e0m, H: Hc, t1: Math.max(tp, 1), t2: tL }) : 0;
      // finite-difference cross-check of the degree of consolidation at handover (uniform equivalent layer)
      if (tH > 0) { const mv = sc / (q * Hc); const fd = consolidationFD({ layers: [{ thickness: Hc, k: cv * mv * 9810, mv }], load: () => q, tEnd: tH, nSteps: 200, nCells: 40, top: 'drained', bottom: inp.drain === 'double' ? 'drained' : 'undrained' }); fdU = fd.U[fd.U.length - 1]; } else fdU = 0;
    }
    const sHand = interp(curve.t, curve.s, tH), sLife = curve.s[curve.s.length - 1], sPost = Math.max(sLife - sHand, 0);
    // waiting time that would bring the post-handover settlement within tolerance (bisection on the curve)
    let tNeed = tH; if (sPost > inp.sLim / 1e3) { let a = tH, b = tL; for (let i = 0; i < 40; i++) { const m = 0.5 * (a + b); if (sLife - interp(curve.t, curve.s, m) > inp.sLim / 1e3) a = m; else b = m; } tNeed = b; }
    const step = sPost + 0; // the piled deck does not follow the backland: the post-handover settlement appears as a step at the transition
    const checks = [
      check('Deck beam long-term deflection', dTotal * 1e3, dAllow * 1e3, 'mm', 'EN 1992-1-1:2004 §7.4.1(4), §7.4.3', 'Sag of the deck beam between pile rows under the quasi-permanent loads with cracking, creep and shrinkage. Excess sag ponds water on the deck and takes crane rails out of level.'),
      check('Deck sway (serviceability)', sway * 1e3, swayAllow * 1e3, 'mm', 'BS 6349-2:2019 §5; crane supplier tolerances', 'Horizontal movement of the deck under service berthing and operating loads. It strains expansion joints, services and crane rail alignment.'),
      check('Angular distortion between pile rows', dist, 1 / inp.distLim, '-', 'EN 1997-1:2004 Annex H', 'Differential pile settlement over the row spacing. Too much distortion cracks the deck beam over the stiffer piles and changes the crane rail cross-level.'),
      check('Backland settlement after handover', sPost * 1e3, inp.sLim, 'mm', 'EN 1997-1:2004 §2.4.9, §6.6.2; BS 6349-2:2019', 'Settlement of the fill and yard still to come after paving, which shows as a step against the piled deck. It trips equipment, breaks drainage falls and needs transition slabs or regrading.'),
    ];
    const w = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = w.uc <= 1;
    const sc2 = deformScale(B.model, qp.disp), hNeed = inp.capH * Math.cbrt(Math.max(dTotal / dAllow, 0.05));
    const yrs = curve.t.map((t) => r(t / YEAR_S, 4));
    return {
      kpis: [
        { key: 'dShort', label: 'Deck beam deflection, short term', value: r(dShort.deflection * 1e3), unit: 'mm', explain: 'Instantaneous sag under the characteristic loads with the short-term modulus and tension stiffening between cracks.' },
        { key: 'dLong', label: 'Deck beam deflection, long term', value: r(dTotal * 1e3), unit: 'mm', limit: r(dAllow * 1e3), sense: 'max', explain: `Sag under the quasi-permanent loads with creep (effective modulus) plus ${f(dShr * 1e3)} mm from shrinkage, against span over ${f(inp.spanLim)}.` },
        { key: 'phi', label: 'Creep coefficient', value: r(phi), unit: '-', explain: 'Ratio of creep strain to elastic strain at the end of life for this member size, humidity and age at loading. Long-term curvature is about (1 + φ) times the elastic value.' },
        { key: 'zeta', label: 'Cracking distribution coefficient ζ', value: r(dLong.zeta), unit: '-', explain: '0 for an uncracked beam, approaching 1 when the moment is far above the cracking moment. It weights the cracked and uncracked curvatures.' },
        { key: 'Mcr', label: 'Cracking moment of the deck beam', value: r(dShort.Mcr / 1e3), unit: 'kN·m', explain: `Moment at which the bottom fibre reaches the tensile strength. The quasi-permanent span moment is ${f(Mqp / 1e3)} kN·m.` },
        { key: 'sway', label: 'Deck sway, serviceability', value: r(sway * 1e3), unit: 'mm', limit: r(swayAllow * 1e3), sense: 'max', explain: 'Largest horizontal deck movement under the service berthing and characteristic combinations.' },
        { key: 'sPile', label: 'Largest pile-head settlement', value: r(sMax * 1e3, 7), unit: 'mm', explain: 'Settlement of the most loaded vertical pile under characteristic loads: load-transfer response of shaft and toe plus elastic shortening above the bed.' },
        { key: 'dDiff', label: 'Differential settlement between adjacent rows', value: r(dDiff * 1e3), unit: 'mm', explain: 'Largest difference in pile settlement between neighbouring rows, caused by their different loads. It is what bends the deck beam, not the total settlement.' },
        { key: 'distortion', label: 'Angular distortion (1 / value)', value: r(dist > 0 ? Math.min(1 / dist, 99999) : 99999), unit: '-', limit: inp.distLim, sense: 'min', explain: 'Row spacing over differential settlement. Larger is better; it must stay above the limit entered.' },
        { key: 'sc', label: 'Primary consolidation settlement of the backland', value: r(sc * 1e3), unit: 'mm', explain: `Final settlement of the clay under ${f(q / 1e3)} kPa of fill and yard load, from the e–log σ′ model.` },
        { key: 'U', label: 'Degree of consolidation at handover', value: r(Uh), unit: '-', explain: `Share of the primary settlement complete when the paving is finished (${f(fdU)} by the finite-difference solution). The remainder appears in service.` },
        { key: 'sPost', label: 'Backland settlement after handover', value: r(sPost * 1e3), unit: 'mm', limit: inp.sLim, sense: 'max', explain: `Primary settlement still to come plus secondary compression (${f(sSec * 1e3)} mm) over the design life: the step that grows at the wharf–backland transition.` },
        { key: 'tNeed', label: 'Waiting time needed before handover', value: r(tNeed / YEAR_S * 12), unit: 'months', explain: 'Time after filling at which the remaining settlement falls to the tolerance. Vertical drains or a surcharge shorten it.' },
      ],
      checks,
      plots: [
        framePlot(B.model, { disp: qp.disp, scale: sc2, title: `Bent under quasi-permanent loads (× ${sc2})`, values: B.model.nodes.map((_, i) => Math.abs(qp.disp[6 * i + 2]) * 1e3), valueLabel: 'Vertical displacement', valueUnit: 'mm',
          explain: 'Elastic displaced shape of the bent under long-term loads with piles on rigid supports at fixity; pile settlement from the t–z analysis is added separately.' }),
        { type: 'line', title: 'Pile load–settlement curve', x: { label: 'Head settlement', unit: 'mm' }, y: { label: 'Head load', unit: 'kN' },
          series: [{ name: 'Pile (t–z, Q–z)', x: ax.w.map((v) => r(v * 1e3)), y: ax.Q.map((v) => r(v / 1e3)) }, { name: 'Row loads', x: sRows.map((q2) => r((q2.s - q2.N * hFix / (B.P.E * B.P.A)) * 1e3)), y: sRows.map((q2) => r(q2.N / 1e3)) }],
          explain: 'Axial response of one pile with the service load of each row marked on it. Rows on the steep initial branch settle a few millimetres; a row near the knee signals too little penetration.' },
        { type: 'line', title: 'Backland settlement with time', x: { label: 'Time since filling', unit: 'years' }, y: { label: 'Settlement', unit: 'mm' },
          series: [{ name: 'Settlement', x: yrs, y: curve.s.map((v) => r(-v * 1e3)) }, { name: 'Handover', x: [r(tH / YEAR_S, 4), r(tH / YEAR_S, 4)], y: [0, r(-sLife * 1e3)], dash: true }],
          explain: 'Consolidation of the clay under the fill, followed by secondary compression. Everything to the right of the handover line happens under the finished yard and shows at the transition to the deck.' },
      ],
      tables: [
        { title: 'Pile settlement by row', columns: [{ key: 'id', label: 'Row' }, { key: 'N', label: 'Service compression', unit: 'kN' }, { key: 's', label: 'Settlement', unit: 'mm' }],
          rows: sRows.map((q2) => ({ id: q2.id, N: r(q2.N / 1e3), s: r(q2.s * 1e3) })), explain: 'Service load and settlement of each vertical pile row. Differences between neighbouring rows give the angular distortion of the deck.' },
        { title: 'Deflection build-up of the deck beam', columns: [{ key: 'item', label: 'Component' }, { key: 'value', label: 'Value' }, { key: 'unit', label: 'Unit' }],
          rows: [{ item: 'Span moment, characteristic', value: r(Mch / 1e3), unit: 'kN·m' }, { item: 'Span moment, quasi-permanent', value: r(Mqp / 1e3), unit: 'kN·m' }, { item: 'Cracking moment', value: r(dShort.Mcr / 1e3), unit: 'kN·m' }, { item: 'Short-term deflection', value: r(dShort.deflection * 1e3), unit: 'mm' }, { item: 'Long-term deflection (creep)', value: r(dLong.deflection * 1e3), unit: 'mm' }, { item: 'Shrinkage deflection', value: r(dShr * 1e3), unit: 'mm' }, { item: 'Effective modulus', value: r(Eeff / 1e9), unit: 'GPa' }, { item: 'Moment resistance (reference)', value: r(momentCapacity({ section: rc.sag, material: rc.material }).MRd / 1e3), unit: 'kN·m' }],
          explain: 'How the long-term sag is built from load, cracking, creep and shrinkage, so that each assumption can be checked against site measurements.' },
      ],
      notes: [`${w.name} governs at ${f(w.uc)}. Creep coefficient ${f(phi)}${inp.phiUser > 0 ? ' (entered)' : ` (EN 1992-1-1 Annex B, notional size ${f(h0)} mm)`}; the beam is ${dLong.zeta > 0 ? 'cracked' : 'uncracked'} under long-term loads.`,
        `Backland: ${f(Hc)} m of clay, drainage path ${f(Hdr)} m, end of primary consolidation after about ${f(tp / YEAR_S)} years; ${f(Uh * 100)} % complete at handover.`],
      advice: [
        ok ? eng('ok', `Deck beam sag ${f(dTotal * 1e3)} mm against ${f(dAllow * 1e3)} mm, sway ${f(sway * 1e3)} mm against ${f(swayAllow * 1e3)} mm, pile settlement ${f(sMax * 1e3)} mm with ${f(dDiff * 1e3)} mm between rows; the backland settles ${f(sPost * 1e3)} mm after handover.`,
          'Detail a transition slab and adjustable rail fixings for the step at the backland, and record deck levels at handover as the baseline for monitoring. Calibrate the cracked stiffness with a measured deflection under a known load if one is available.', 'EN 1992-1-1:2004 §7.4; EN 1997-1:2004 Annex H; BS 6349-2:2019 §5',
          `${dTotal / dAllow < 0.5 ? 'Deflection leaves a large reserve, so strength and detailing set the beam depth' : 'Deflection is close to governing the beam depth'}: on stiffness alone about ${f(Math.max(hNeed, 0.4))} m would meet the limit, ${f(Math.abs(Math.max(hNeed, 0.4) - inp.capH) * inp.capB * B.deckLen)} m³ per bent away from the present section (${concCo2((Math.max(hNeed, 0.4) - inp.capH) * inp.capB * B.deckLen)}).`)
          : eng('fail', `${w.name}: ${f(w.demand)} against ${f(w.capacity)} ${w.unit}.`,
            w.name.startsWith('Backland') ? `Wait about ${f(tNeed / YEAR_S * 12)} months before paving, or accelerate consolidation with vertical drains and a temporary surcharge of at least the yard load; otherwise design the transition for ${f(sPost * 1e3)} mm.` : w.name.startsWith('Deck beam') ? `Deepen the deck beam to about ${f(hNeed)} m, add compression steel to reduce creep curvature, or pre-camber the formwork by the long-term value.` : w.name.startsWith('Angular') ? 'Lengthen the heavily loaded piles so that all rows work on the stiff part of the load–settlement curve, or rebalance the row loads.' : 'Stiffen the bent with larger piles or rakers, or reduce the free length.',
            'EN 1992-1-1:2004 §7.4; EN 1997-1:2004 §2.4.8–2.4.9',
            `Deepening the beam to ${f(hNeed)} m adds ${f(Math.max(hNeed - inp.capH, 0) * inp.capB * B.deckLen)} m³ per bent: ${concCo2(Math.max(hNeed - inp.capH, 0) * inp.capB * B.deckLen)}. Waiting for consolidation or pre-loading with re-usable fill costs time rather than material.`),
        eco(ok ? 'info' : 'warn', `Handover ${f(inp.tHand)} months after filling leaves ${f(sPost * 1e3)} mm of backland settlement to maintain over ${f(inp.life)} years; the tolerance would be met after ${f(tNeed / YEAR_S * 12)} months.`,
          'Compare the cost of the waiting time or of ground treatment (drains, surcharge) with the life-cycle cost of regrading the yard and re-levelling rails: each regrading closes stacking area and interrupts operations.', `${ECON_BASIS}; observational method per EN 1997-1 §2.7`),
      ],
      publish: { 'wharf.deckDeflection': dTotal, 'wharf.pileSettlement': sMax, 'wharf.backlandSettlement': sPost, 'wharf.differentialSettlement': Math.max(dDiff, step) },
    };
  },
  verify: [
    // EN 1992-1-1 Annex B, C40 (fcm = 48 MPa), beam 1.2 × 1.4 m (h0 = 2·1.68/5.2 = 646 mm), RH 80 %, t0 = 28 days:
    // α1 = (35/48)^0.7 = 0.8016, α2 = (35/48)^0.2 = 0.9388; φRH = [1 + 0.2/(0.1·646^⅓)·0.8016]·0.9388 = 1.1129;
    // β(fcm) = 16.8/√48 = 2.4249; β(t0) = 1/(0.1 + 28^0.2) = 0.4885 → φ0 = 1.318.
    { name: 'Creep coefficient of EN 1992-1-1 Annex B', source: 'EN 1992-1-1:2004 Eq. B.1–B.5, B.8c (hand calculation)', inputs: { capB: 1.2, capH: 1.4, conc: 'c40', RH: 80, t0: 28, phiUser: 0 }, expect: [{ kpi: 'phi', value: 1.318, tol: 0.005 }] },
    // Normally consolidated clay under a wide load: s = Cc/(1+e0)·∫ log10[(γ′z + q)/(γ′z)] dz over H.
    // H = 4 m, γ′ = 8 kN/m³, q = 40 kPa (no fill, yard load only), Cc = 0.3, e0 = 1; with a 1 kPa overburden at the top:
    // ∫ = [(H + a + b)·ln(γ′(H + a + b)) − (H + a)·ln(γ′(H + a)) − (a + b)·ln(γ′(a + b)) + a·ln(γ′a)]/ln10, a = 0.125 m, b = 5 m → 2.473 m; s = 0.15·2.473 = 0.371 m.
    { name: 'Primary consolidation of normally consolidated clay under a wide load', source: 'Terzaghi & Peck e–log σ′ relation integrated over depth (hand integration)',
      inputs: { fillH: 0, qYard: 40, clay: [{ name: 'NC clay', thick: 4, gamma: 8, Cc: 0.3, Cr: 0.05, e0: 1, ocr: 1 }] }, expect: [{ kpi: 'sc', value: 371, tol: 0.01 }] },
    // Terzaghi: U = 50 % at Tv = 0.197. Double drainage of 4 m (Hdr = 2 m), cv = 2 m²/yr → t = 0.197·4/2 = 0.394 yr = 4.728 months.
    { name: 'Degree of consolidation U = 50 % at Tv = 0.197', source: 'Terzaghi (1943) one-dimensional consolidation',
      inputs: { fillH: 0, qYard: 40, cv: 2, drain: 'double', tHand: 4.728, clay: [{ name: 'NC clay', thick: 4, gamma: 8, Cc: 0.3, Cr: 0.05, e0: 1, ocr: 1 }] }, expect: [{ kpi: 'U', value: 0.5, tol: 0.01 }] },
  ],
};
