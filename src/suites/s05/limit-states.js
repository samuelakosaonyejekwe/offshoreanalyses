// Ultimate and serviceability limit states of a reinforced-concrete wharf member (deck beam, slab strip, pile cap or
// circular pile): the N–M interaction diagram and moment–curvature response, bending, shear, ductility and minimum
// reinforcement at the ultimate limit state, stress limits and crack width in service, the permissible moments,
// and an overall compliance table for the load combinations entered.
import { rcSection, circularBars, concreteProperties, interactionDiagram, momentCapacity, momentCurvature, crackingMoment, crackedSection, crackWidth, shearResistanceNoLinks, shearResistanceLinks } from '../../engine/structural/concrete.js';
import { need, run, r, f, check, eng, eco, thin, CONFIRM, NOT_HIFI, RC_NOTE, ECON_BASIS, CO2_STEEL, CO2_CONC } from './_common.js';

const COMBOS_DEFAULT = [
  { name: 'ULS operating (cargo + crane)', ls: 'ULS', N: 0, M: 3300, V: 2700 },
  { name: 'ULS berthing', ls: 'ULS', N: 350, M: -3500, V: 1900 },
  { name: 'ULS mooring + environment', ls: 'ULS', N: -250, M: -2600, V: 1500 },
  { name: 'SLS characteristic', ls: 'SLS-char', N: 0, M: 2300, V: 0 },
  { name: 'SLS quasi-permanent', ls: 'SLS-qp', N: 0, M: 1500, V: 0 },
  { name: 'SLS quasi-permanent, support', ls: 'SLS-qp', N: 0, M: -1700, V: 0 },
];

/** Section pair (sagging: bottom in tension; hogging: the member turned over) and its design data, SI. */
function build(inp, ctx) {
  const cm = ctx?.material ? ctx.material(inp.conc) : null, sm = ctx?.material ? ctx.material(inp.rebar) : null;
  need(cm && cm.fck > 0, `Concrete '${inp.conc}' is not in the materials library.`);
  need(sm && sm.fy > 0, `Reinforcing steel '${inp.rebar}' is not in the materials library.`);
  const c = inp.cover / 1e3, dl = inp.dLink / 1e3, nS = Math.round(inp.nStrips);
  const material = { fck: cm.fck, fyk: sm.fy, gammaC: inp.gammaC, gammaS: inp.gammaS, alphaCC: inp.alphaCC }, props = concreteProperties({ fck: cm.fck });
  if (inp.shape === 'circle') {
    const db = inp.dBot / 1e3, n = Math.round(inp.nBot);
    need(n >= 4, 'A circular section needs at least four bars.');
    const bars = run(() => circularBars({ D: inp.D, cover: c + dl + db / 2, nBars: n, barDiameter: db }), 'Bar layout');
    const sec = rcSection({ shape: 'circle', D: inp.D, bars, nStrips: nS }), As = n * Math.PI * db * db / 4;
    return { circle: true, sag: sec, hog: sec, material, props, b: 0.8 * inp.D, h: inp.D, d: sec.dMax, dHog: sec.dMax, AsB: As / 2, AsT: As / 2, As, db, dt: db, cover: c + dl, fck: cm.fck, fyk: sm.fy, Es: sm.E ?? 200e9, area: Math.PI * inp.D ** 2 / 4, spacing: Math.PI * (inp.D - 2 * (c + dl) - db) / n };
  }
  const b = inp.b, h = inp.h, db = inp.dBot / 1e3, dt = inp.dTop / 1e3, AsB = inp.nBot * Math.PI * db * db / 4, AsT = inp.nTop * Math.PI * dt * dt / 4;
  const d = h - c - dl - db / 2, dTopDepth = c + dl + dt / 2;
  need(d > dTopDepth + 0.02, `The section depth of ${f(h * 1e3)} mm is too small for the cover and bars given.`);
  need(AsB > 0, 'Give bottom reinforcement.');
  const lay = (arr) => arr.filter((q) => q.area > 0);
  const sag = rcSection({ shape: 'rect', b, h, bars: lay([{ depth: dTopDepth, area: AsT, diameter: dt }, { depth: d, area: AsB, diameter: db }]), nStrips: nS });
  const hog = AsT > 0 ? rcSection({ shape: 'rect', b, h, bars: lay([{ depth: h - d, area: AsB, diameter: db }, { depth: h - dTopDepth, area: AsT, diameter: dt }]), nStrips: nS }) : null;
  return { circle: false, sag, hog, material, props, b, h, d, dHog: h - dTopDepth, AsB, AsT, As: AsB + AsT, db, dt, cover: c + dl, fck: cm.fck, fyk: sm.fy, Es: sm.E ?? 200e9, area: b * h, spacing: (b - 2 * (c + dl) - db) / Math.max(inp.nBot - 1, 1), spacingTop: (b - 2 * (c + dl) - dt) / Math.max(inp.nTop - 1, 1) };
}

export default {
  id: 's05-limit-states', suite: 5, title: 'Ultimate and serviceability limit states',
  summary: 'Checks a reinforced-concrete deck beam, slab strip or pile for every load combination: bending with axial force, shear, ductility and minimum steel at the ultimate limit state, stresses and crack width in service, with the permissible moments and an overall compliance verdict.',
  covers: [12], fidelity: 'analytical',
  standards: ['EN 1992-1-1:2004 §3.1.7, §6.1, §6.2, §7.2, §7.3, §9.2.1.1', 'EN 1990:2002+A1 §6.4, §6.5 (limit states and combinations)', 'BS 6349-2:2019 §5–6', 'EN 1992-2:2005 (αcc and crack limits for exposed structures)'],
  models: ['Reinforced-concrete sectional analysis (strain compatibility on horizontal strips)', 'Parabola–rectangle concrete compression stress–strain model', 'Bilinear steel reinforcement plasticity', 'N–M interaction diagram', 'Moment–curvature relationship', 'Cracked-section elastic analysis', 'Crack width (EN 1992 bond model)', 'Shear resistance: concrete alone and variable-angle truss'],
  limits: `Uniaxial bending with axial force on a rectangular or circular section; biaxial bending, torsion, slenderness (second-order) effects, prestress and fatigue of the reinforcement are not covered. For circular sections shear and crack width are reported on an equivalent web of 0.8 D as an indication only. ${RC_NOTE} ${NOT_HIFI}`,
  publishes: ['wharf.sectionMRd', 'wharf.sectionVRd', 'wharf.sectionUc'],
  searchTerms: 'reinforced concrete limit state design crack width interaction diagram marine deck beam',
  inputs: [
    { key: 'shape', label: 'Section shape', type: 'select', default: 'rect', group: 'Section', role: 'geometry', options: [{ value: 'rect', label: 'Rectangular (beam, slab strip, pile cap)' }, { value: 'circle', label: 'Circular (reinforced-concrete pile or column)' }] },
    { key: 'b', label: 'Width', unit: 'm', type: 'number', default: 1.2, min: 0.15, max: 6, group: 'Section', role: 'geometry', showIf: { key: 'shape', in: ['rect'] }, help: '1.0 for a slab strip per metre.' },
    { key: 'h', label: 'Depth', unit: 'm', type: 'number', default: 1.4, min: 0.15, max: 5, group: 'Section', role: 'geometry', design: true, showIf: { key: 'shape', in: ['rect'] } },
    { key: 'D', label: 'Diameter', unit: 'm', type: 'number', default: 1.2, min: 0.3, max: 4, group: 'Section', role: 'geometry', showIf: { key: 'shape', in: ['circle'] } },
    { key: 'conc', label: 'Concrete class', type: 'select', source: 'materials:concrete', default: 'c40', group: 'Section', role: 'material' },
    { key: 'nBot', label: 'Bottom bars (all bars for a circular section)', unit: '-', type: 'int', default: 12, min: 2, max: 80, group: 'Reinforcement', role: 'geometry', design: true },
    { key: 'dBot', label: 'Bottom bar diameter', unit: 'mm', type: 'number', default: 32, min: 8, max: 50, group: 'Reinforcement', role: 'geometry', design: true },
    { key: 'nTop', label: 'Top bars', unit: '-', type: 'int', default: 12, min: 0, max: 80, group: 'Reinforcement', role: 'geometry', showIf: { key: 'shape', in: ['rect'] } },
    { key: 'dTop', label: 'Top bar diameter', unit: 'mm', type: 'number', default: 32, min: 8, max: 50, group: 'Reinforcement', role: 'geometry', showIf: { key: 'shape', in: ['rect'] } },
    { key: 'cover', label: 'Nominal cover to the links', unit: 'mm', type: 'number', default: 60, min: 15, max: 120, group: 'Reinforcement', role: 'geometry', uncertain: { dist: 'normal', cov: 0.12 }, help: 'Marine exposure XS3 typically needs 50–75 mm.' },
    { key: 'dLink', label: 'Link diameter (0 = none)', unit: 'mm', type: 'number', default: 16, min: 0, max: 25, group: 'Reinforcement', role: 'geometry' },
    { key: 'linkLegs', label: 'Link legs', unit: '-', type: 'int', default: 4, min: 0, max: 20, group: 'Reinforcement', role: 'geometry' },
    { key: 'linkS', label: 'Link spacing', unit: 'mm', type: 'number', default: 200, min: 50, max: 600, group: 'Reinforcement', role: 'geometry' },
    { key: 'rebar', label: 'Reinforcing steel', type: 'select', source: 'materials:steel', default: 'reinf-b500', group: 'Reinforcement', role: 'material' },
    { key: 'combos', label: 'Load combinations at the section', type: 'table', group: 'Actions', role: 'bc',
      columns: [{ key: 'name', label: 'Combination', type: 'text' }, { key: 'ls', label: 'Limit state (ULS / SLS-char / SLS-qp)', type: 'text' }, { key: 'N', label: 'Axial force (compression +)', unit: 'kN', type: 'number' }, { key: 'M', label: 'Moment (sagging +)', unit: 'kN·m', type: 'number' }, { key: 'V', label: 'Shear', unit: 'kN', type: 'number' }],
      default: COMBOS_DEFAULT, help: 'Design values for ULS rows, unfactored values for the serviceability rows, from the strength and deck modules. Negative moments put the top in tension.' },
    { key: 'gammaC', label: 'γc concrete', unit: '-', type: 'number', default: 1.5, min: 1, max: 2, group: 'Resistance factors', role: 'criteria', help: `EN 1992-1-1 §2.4.2.4; 1.2 for accidental situations. ${CONFIRM}` },
    { key: 'gammaS', label: 'γs reinforcement', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.5, group: 'Resistance factors', role: 'criteria', help: CONFIRM },
    { key: 'alphaCC', label: 'αcc long-term strength coefficient', unit: '-', type: 'number', default: 0.85, min: 0.7, max: 1, group: 'Resistance factors', role: 'criteria', help: 'EN 1992-1-1 §3.1.6 recommends 1.0; several national annexes and EN 1992-2 use 0.85. Confirm for the governing code.' },
    { key: 'k1c', label: 'Concrete stress limit, characteristic (fraction of fck)', unit: '-', type: 'number', default: 0.6, min: 0.3, max: 1, group: 'Serviceability criteria', role: 'criteria', help: `EN 1992-1-1 §7.2(2) for exposure classes XD, XF, XS. ${CONFIRM}` },
    { key: 'k2c', label: 'Concrete stress limit, quasi-permanent (fraction of fck)', unit: '-', type: 'number', default: 0.45, min: 0.2, max: 1, group: 'Serviceability criteria', role: 'criteria', help: 'Above this creep becomes non-linear (EN 1992-1-1 §7.2(3)).' },
    { key: 'k3s', label: 'Steel stress limit, characteristic (fraction of fyk)', unit: '-', type: 'number', default: 0.8, min: 0.4, max: 1, group: 'Serviceability criteria', role: 'criteria', help: 'EN 1992-1-1 §7.2(5).' },
    { key: 'wLim', label: 'Crack width limit', unit: 'mm', type: 'number', default: 0.3, min: 0.05, max: 0.5, group: 'Serviceability criteria', role: 'criteria', help: `EN 1992-1-1 Table 7.1N, XS exposure, quasi-permanent combination; 0.2 mm is often specified in the splash zone. ${CONFIRM}` },
    { key: 'nStrips', label: 'Concrete strips over the section depth', unit: '-', type: 'int', default: 200, min: 20, max: 1600, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nStrips', label: 'Strips over the depth', levels: [25, 50, 100, 200, 400], qoi: 'MRd' },
  compute(inp, ctx) {
    const S = build(inp, ctx), Ecm = S.props.Ecm;
    const combos = (inp.combos || []).map((q, i) => { const ls = String(q.ls || '').trim().toUpperCase(); need(['ULS', 'SLS-CHAR', 'SLS-QP'].includes(ls), `Combination ${i + 1}: the limit state must be ULS, SLS-char or SLS-qp.`); const N = Number(q.N) || 0, M = Number(q.M) || 0, V = Math.abs(Number(q.V) || 0); return { name: String(q.name || `Combination ${i + 1}`), ls, N: N * 1e3, M: M * 1e3, V: V * 1e3, hog: M < 0 }; });
    need(combos.some((q) => q.ls === 'ULS'), 'Give at least one ULS combination.');
    const secOf = (q) => { if (q.hog && !S.hog) throw new Error(`Combination '${q.name}' puts the top in tension but the section has no top bars.`); return q.hog ? S.hog : S.sag; };
    const dia = run(() => interactionDiagram({ section: S.sag, material: S.material, nPoints: 61 }), 'Interaction diagram');
    const diaH = S.hog && !S.circle ? interactionDiagram({ section: S.hog, material: S.material, nPoints: 61 }) : dia;
    const MRd0 = momentCapacity({ section: S.sag, material: S.material }), MRd0H = S.hog ? momentCapacity({ section: S.hog, material: S.material }) : { MRd: 0, xOverD: 0 };
    const AswS = inp.dLink > 0 && inp.linkLegs > 0 ? inp.linkLegs * Math.PI * (inp.dLink / 1e3) ** 2 / 4 / (inp.linkS / 1e3) : 0;
    const shearOf = (q) => { const d = q.hog ? S.dHog : S.d, Asl = q.hog ? S.AsT : S.AsB; const c0 = shearResistanceNoLinks({ bw: S.b, d, fck: S.fck, Asl, NEd: Math.max(Math.min(q.N, 0.2 * S.fck / inp.gammaC * S.area), 0), Ac: S.area, gammaC: inp.gammaC, alphaCC: inp.alphaCC }).VRdc; if (!(AswS > 0)) return { VRd: c0, VRdc: c0 }; const l = shearResistanceLinks({ bw: S.b, d, fck: S.fck, fywk: S.fyk, AswOverS: AswS, VEd: q.V, gammaC: inp.gammaC, gammaS: inp.gammaS, alphaCC: inp.alphaCC }); return { VRd: Math.max(l.VRd, c0), VRdc: c0 }; };
    const rows = combos.map((q) => {
      const sec = secOf(q), M = Math.abs(q.M);
      if (q.ls === 'ULS') {
        need(q.N < dia.NRdMax && q.N > dia.NRdMin, `Combination '${q.name}': the axial force of ${f(q.N / 1e3)} kN is outside the axial resistance of the section (${f(dia.NRdMin / 1e3)} to ${f(dia.NRdMax / 1e3)} kN).`);
        const mc = momentCapacity({ section: sec, material: S.material, N: q.N }), sh = shearOf(q);
        return { ...q, MRd: mc.MRd, xd: mc.xOverD, ucM: M / Math.max(mc.MRd, 1), VRd: sh.VRd, ucV: q.V / Math.max(sh.VRd, 1) };
      }
      // serviceability: uncracked below the cracking moment, cracked elastic section above it
      const un = crackingMoment({ section: sec, fctm: S.props.fctm, Ec: Ecm, Es: S.Es, N: q.N });
      const cracked = M > un.Mcr, cs = crackedSection({ section: sec, Ec: Ecm, Es: S.Es, M: Math.max(M, 1) });
      const sigC = cracked ? cs.sigmaC : M * un.centroid / un.inertia + q.N / un.area, sigS = cracked ? cs.sigmaS : (S.Es / Ecm) * M * (sec.dMax - un.centroid) / un.inertia;
      let wk = 0;
      if (q.ls === 'SLS-QP' && cracked) wk = crackWidth({ sigmaS: cs.sigmaS, As: q.hog ? S.AsT : S.AsB, b: S.b, h: S.h, d: q.hog ? S.dHog : S.d, x: cs.x, cover: S.cover, barDiameter: q.hog ? S.dt : S.db, fctEff: S.props.fctm, Ecm, Es: S.Es, barSpacing: q.hog ? (S.spacingTop ?? S.spacing) : S.spacing }).wk;
      const limC = (q.ls === 'SLS-CHAR' ? inp.k1c : inp.k2c) * S.fck;
      return { ...q, Mcr: un.Mcr, cracked, sigC, sigS, wk, ucC: sigC / limC, ucS: q.ls === 'SLS-CHAR' ? sigS / (inp.k3s * S.fyk) : 0, ucW: q.ls === 'SLS-QP' ? wk * 1e3 / inp.wLim : 0 };
    });
    const uls = rows.filter((q) => q.ls === 'ULS'), sls = rows.filter((q) => q.ls !== 'ULS');
    const gM = uls.reduce((m, q) => (q.ucM > m.ucM ? q : m), uls[0]), gV = uls.reduce((m, q) => (q.ucV > m.ucV ? q : m), uls[0]);
    const gC = sls.length ? sls.reduce((m, q) => (q.ucC > m.ucC ? q : m), sls[0]) : null, gS = sls.length ? sls.reduce((m, q) => (q.ucS > m.ucS ? q : m), sls[0]) : null, gW = sls.length ? sls.reduce((m, q) => (q.ucW > m.ucW ? q : m), sls[0]) : null;
    // ductility and minimum reinforcement (EN 1992-1-1 §5.6.3: x/d ≤ 0.45 up to C50; §9.2.1.1: As,min = 0.26 fctm/fyk·b·d ≥ 0.0013 b d)
    const xdLim = S.fck <= 50e6 ? 0.45 : 0.35, AsMin = Math.max(0.26 * S.props.fctm / S.fyk, 0.0013) * S.b * S.d;
    const mk = run(() => momentCurvature({ section: S.sag, material: S.material, N: 0, nPoints: 60 }), 'Moment–curvature');
    const mkMean = momentCurvature({ section: S.sag, material: { fck: S.props.fcm, fyk: S.fyk, gammaC: 1, gammaS: 1, alphaCC: 1 }, N: 0, nPoints: 60 });
    // permissible serviceability moment (sagging, N = 0): crack width and stress limits by bisection on the cracked section
    const slsOk = (M) => { const un = crackingMoment({ section: S.sag, fctm: S.props.fctm, Ec: Ecm, Es: S.Es }); if (M <= un.Mcr) return true; const cs = crackedSection({ section: S.sag, Ec: Ecm, Es: S.Es, M }); const wk = crackWidth({ sigmaS: cs.sigmaS, As: S.AsB, b: S.b, h: S.h, d: S.d, x: cs.x, cover: S.cover, barDiameter: S.db, fctEff: S.props.fctm, Ecm, Es: S.Es, barSpacing: S.spacing }).wk; return wk * 1e3 <= inp.wLim && cs.sigmaC <= inp.k2c * S.fck; };
    let Mqp = 0; { let a = 0, b = MRd0.MRd; if (slsOk(b)) Mqp = b; else { for (let i = 0; i < 40; i++) { const m = 0.5 * (a + b); if (slsOk(m)) a = m; else b = m; } Mqp = a; } }
    const checks = [
      check(`Bending with axial force (${gM.name})`, Math.abs(gM.M) / 1e3, gM.MRd / 1e3, 'kN·m', 'EN 1992-1-1:2004 §6.1', 'Design moment against the moment resistance at the coexisting axial force. If exceeded the tension steel yields and the compression zone crushes: flexural failure of the member.'),
      check(`Shear (${gV.name})`, gV.V / 1e3, gV.VRd / 1e3, 'kN', 'EN 1992-1-1:2004 §6.2.2, §6.2.3', 'Design shear against the resistance of the links (variable-angle truss) or of the concrete alone. Shear failure is sudden, so this check should keep a clear margin.'),
      check('Ductility: neutral axis depth x/d', Math.max(MRd0.xOverD, S.hog && !S.circle ? MRd0H.xOverD : 0), xdLim, '-', 'EN 1992-1-1:2004 §5.6.3(2)', 'Depth of the compression zone at failure. A shallow zone means the steel yields well before the concrete crushes, so the member gives warning and can redistribute moment.'),
      check('Minimum tension reinforcement', AsMin * 1e6, Math.min(S.AsB, S.hog && !S.circle ? S.AsT : S.AsB) * 1e6, 'mm²', 'EN 1992-1-1:2004 §9.2.1.1', 'Steel needed to carry the cracking moment against the steel provided. Below it the member can fail suddenly when it first cracks.'),
      ...(gC ? [check(`Concrete compressive stress in service (${gC.name})`, gC.sigC / 1e6, (gC.ls === 'SLS-CHAR' ? inp.k1c : inp.k2c) * S.fck / 1e6, 'MPa', 'EN 1992-1-1:2004 §7.2(2), (3)', 'Extreme-fibre concrete stress under service loads. Excess stress causes longitudinal micro-cracking and non-linear creep, which reduce durability.')] : []),
      ...(gS && gS.ucS > 0 ? [check(`Reinforcement stress in service (${gS.name})`, gS.sigS / 1e6, inp.k3s * S.fyk / 1e6, 'MPa', 'EN 1992-1-1:2004 §7.2(5)', 'Steel stress under the characteristic combination. Above the limit the bars may yield in service and cracks stay open permanently.')] : []),
      ...(gW && rows.some((q) => q.ls === 'SLS-QP') ? [check(`Crack width (${gW.name})`, gW.wk * 1e3, inp.wLim, 'mm', 'EN 1992-1-1:2004 §7.3.4, Table 7.1N', 'Calculated crack width under long-term loads. Wider cracks shorten the time for chlorides to reach the bars: the main durability control in the splash zone.')] : []),
    ];
    const w = checks.reduce((m, c) => (c.uc > m.uc ? c : m)), nFail = checks.filter((c) => c.uc > 1).length, ok = nFail === 0;
    const ucU = Math.max(gM.ucM, gV.ucV), ucSls = Math.max(gC?.ucC ?? 0, gS?.ucS ?? 0, gW?.ucW ?? 0);
    const AsNeed = S.AsB * Math.max(gM.ucM, 0.05), dAs = (AsNeed - S.AsB) * 7850 / 1e3; // t per metre of member
    const [kx, ky] = thin([Array.from(mk.kappa, (v) => v * 1e3), Array.from(mk.M, (v) => v / 1e3)], 120), [kmx, kmy] = thin([Array.from(mkMean.kappa, (v) => v * 1e3), Array.from(mkMean.M, (v) => v / 1e3)], 120);
    return {
      kpis: [
        { key: 'MRd', label: 'Moment resistance, sagging, no axial force', value: r(MRd0.MRd / 1e3, 8), unit: 'kN·m', explain: 'Design bending resistance with the bottom steel in tension: the largest factored sagging moment the section can carry. It is the maximum permissible ULS moment.' },
        { key: 'MRdHog', label: 'Moment resistance, hogging', value: r((S.circle ? MRd0.MRd : MRd0H.MRd) / 1e3), unit: 'kN·m', explain: 'Design bending resistance with the top steel in tension, for sections over supports and pile heads.' },
        { key: 'VRd', label: 'Shear resistance', value: r(gV.VRd / 1e3), unit: 'kN', explain: `Design shear resistance in the governing combination; the concrete alone would give ${f(shearOf(gV).VRdc / 1e3)} kN.` },
        { key: 'NRd', label: 'Axial compression resistance', value: r(dia.NRdMax / 1e3), unit: 'kN', explain: 'Squash load of the section (concrete and steel in uniform compression): the top of the interaction diagram.' },
        { key: 'ucULS', label: 'Governing ultimate utilisation', value: r(ucU), unit: '-', limit: 1, sense: 'max', explain: 'Largest of bending and shear demand over resistance across the ULS combinations. Above 1 the member does not meet the code strength.' },
        { key: 'ucSLS', label: 'Governing serviceability utilisation', value: r(ucSls), unit: '-', limit: 1, sense: 'max', explain: 'Largest of the concrete stress, steel stress and crack width ratios across the serviceability combinations.' },
        { key: 'xd', label: 'Neutral axis depth at failure x/d', value: r(MRd0.xOverD), unit: '-', limit: xdLim, sense: 'max', explain: 'Compression zone depth over effective depth at the ultimate moment. Low values mean a ductile, under-reinforced section.' },
        { key: 'ductility', label: 'Curvature ductility', value: r(Number.isFinite(mk.ductility) ? mk.ductility : 1), unit: '-', explain: 'Ultimate over first-yield curvature with design strengths. Above about 4 the section can form a plastic hinge, which matters for accidental and seismic loading.' },
        { key: 'Mcr', label: 'Cracking moment', value: r(crackingMoment({ section: S.sag, fctm: S.props.fctm, Ec: Ecm, Es: S.Es }).Mcr / 1e3), unit: 'kN·m', explain: 'Moment at which the tension face reaches the mean tensile strength. Below it the member is uncracked and stiff.' },
        { key: 'wk', label: 'Largest crack width, quasi-permanent', value: r((gW?.wk ?? 0) * 1e3), unit: 'mm', limit: inp.wLim, sense: 'max', explain: 'Characteristic crack width under the long-term combinations, against the limit for marine exposure.' },
        { key: 'sigS', label: 'Largest steel stress in service', value: r(Math.max(...sls.map((q) => q.sigS), 0) / 1e6), unit: 'MPa', explain: 'Stress in the tension bars on the cracked section under service loads. It drives crack width and is limited to avoid yielding.' },
        { key: 'Mqp', label: 'Maximum permissible quasi-permanent moment', value: r(Mqp / 1e3), unit: 'kN·m', explain: 'Largest sustained sagging moment for which the crack width and the long-term concrete stress stay within their limits: the serviceability rating of the member.' },
        { key: 'nFail', label: 'Checks not satisfied', value: nFail, unit: '-', limit: 0, sense: 'max', explain: 'Number of limit-state checks above unity. Zero means the member complies at both ultimate and serviceability limit states for every combination entered.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'N–M interaction diagram with the design points', x: { label: 'Moment', unit: 'kN·m' }, y: { label: 'Axial force (compression +)', unit: 'kN' },
          series: [{ name: 'Resistance, sagging', x: Array.from(dia.M, (v) => r(v / 1e3)), y: Array.from(dia.N, (v) => r(v / 1e3)) }, ...(S.circle ? [] : [{ name: 'Resistance, hogging', x: Array.from(diaH.M, (v) => r(-v / 1e3)), y: Array.from(diaH.N, (v) => r(v / 1e3)) }]), { name: 'ULS combinations', x: uls.map((q) => r((S.circle ? Math.abs(q.M) : q.M) / 1e3)), y: uls.map((q) => r(q.N / 1e3)) }],
          explain: 'Envelope of axial force and moment the section can resist, with each ULS combination as a point. Points inside the envelope pass; moderate compression increases the moment resistance up to the balance point.' },
        { type: 'line', title: 'Moment–curvature response (sagging)', x: { label: 'Curvature', unit: '1/km' }, y: { label: 'Moment', unit: 'kN·m' },
          series: [{ name: 'Design strengths', x: kx.map((v) => r(v)), y: ky.map((v) => r(v)) }, { name: 'Mean strengths', x: kmx.map((v) => r(v)), y: kmy.map((v) => r(v)), dash: true }],
          explain: 'Bending response from first loading to crushing of the concrete. The knee is first yield of the steel; the length of the plateau after it is the ductility available for redistribution.' },
        { type: 'bar', title: 'Utilisation by combination', categories: rows.map((q) => q.name), series: [{ name: 'Utilisation', y: rows.map((q) => r(q.ls === 'ULS' ? Math.max(q.ucM, q.ucV) : Math.max(q.ucC, q.ucS, q.ucW), 3)) }], x: { label: 'Combination', unit: '' }, y: { label: 'Utilisation', unit: '-' },
          explain: 'Governing demand over capacity for each combination (bending or shear at ULS; stress or crack width in service). Bars above 1 identify the combinations that need action.' },
      ],
      tables: [
        { title: 'Limit-state compliance by combination', columns: [{ key: 'name', label: 'Combination' }, { key: 'ls', label: 'Limit state' }, { key: 'N', label: 'N', unit: 'kN' }, { key: 'M', label: 'M', unit: 'kN·m' }, { key: 'res', label: 'Resistance / stress' }, { key: 'uc', label: 'Utilisation' }, { key: 'ok', label: 'Result' }],
          rows: rows.map((q) => { const uc = q.ls === 'ULS' ? Math.max(q.ucM, q.ucV) : Math.max(q.ucC, q.ucS, q.ucW); return { name: q.name, ls: q.ls, N: r(q.N / 1e3), M: r(q.M / 1e3), res: q.ls === 'ULS' ? `MRd ${f(q.MRd / 1e3)} kN·m, VRd ${f(q.VRd / 1e3)} kN` : `σc ${f(q.sigC / 1e6)} MPa, σs ${f(q.sigS / 1e6)} MPa${q.ls === 'SLS-QP' ? `, wk ${f(q.wk * 1e3)} mm` : ''}${q.cracked ? '' : ' (uncracked)'}`, uc: r(uc, 3), ok: uc <= 1 ? 'complies' : 'does not comply' }; }),
          explain: 'One row per load combination with its resistance or service stresses and the verdict. This is the overall serviceability and strength compliance statement of the member.' },
        { title: 'Section data', columns: [{ key: 'item', label: 'Quantity' }, { key: 'value', label: 'Value' }, { key: 'unit', label: 'Unit' }],
          rows: [{ item: 'Effective depth', value: r(S.d * 1e3), unit: 'mm' }, { item: 'Tension steel (sagging)', value: r(S.AsB * 1e6), unit: 'mm²' }, { item: 'Tension steel (hogging)', value: r(S.AsT * 1e6), unit: 'mm²' }, { item: 'Reinforcement ratio', value: r(S.AsB / (S.b * S.d) * 100, 3), unit: '%' }, { item: 'Minimum steel', value: r(AsMin * 1e6), unit: 'mm²' }, { item: 'Design concrete strength fcd', value: r(inp.alphaCC * S.fck / inp.gammaC / 1e6), unit: 'MPa' }, { item: 'Design steel strength fyd', value: r(S.fyk / inp.gammaS / 1e6), unit: 'MPa' }, { item: 'Balance point N', value: r(dia.balanced.N / 1e3), unit: 'kN' }, { item: 'Balance point M', value: r(dia.balanced.M / 1e3), unit: 'kN·m' }],
          explain: 'Derived section properties and design strengths used in every check, for comparison with the drawings and with a hand calculation.' },
      ],
      notes: [`${w.name} governs at ${f(w.uc)}; ${nFail === 0 ? 'all checks are satisfied' : `${nFail} check(s) are not satisfied`}. Design strengths with γc = ${f(inp.gammaC)}, γs = ${f(inp.gammaS)}, αcc = ${f(inp.alphaCC)}.`,
        `Ultimate utilisation ${f(ucU)}, serviceability utilisation ${f(ucSls)}: ${ucSls > ucU ? 'serviceability (durability) governs the reinforcement' : 'strength governs the reinforcement'}.`],
      advice: [
        ok ? eng('ok', `The member complies: ultimate utilisation ${f(ucU)} (MRd ${f(MRd0.MRd / 1e3)} kN·m, VRd ${f(gV.VRd / 1e3)} kN), serviceability utilisation ${f(ucSls)} (crack width ${f((gW?.wk ?? 0) * 1e3)} mm).`,
          `The section can take a sustained moment of ${f(Mqp / 1e3)} kN·m within the crack limit and a factored moment of ${f(MRd0.MRd / 1e3)} kN·m. Confirm the partial factors and crack limit for the governing code and exposure class before relying on the reserve.`, 'EN 1992-1-1:2004 §6–7; EN 1990:2002+A1 §6',
          `Bending needs about ${f(AsNeed * 1e6)} mm² of tension steel against ${f(S.AsB * 1e6)} mm² provided, a difference of ${f(Math.abs(dAs) * 1e3)} kg per metre of member (${f(Math.abs(dAs) * CO2_STEEL * 1e3)} kg CO₂e per metre, screening factor); keep enough steel for the crack limit, which protects the service life.`)
          : eng('fail', `${nFail} check(s) fail; ${w.name} governs: ${f(w.demand)} against ${f(w.capacity)} ${w.unit} (${f(w.uc)}).`,
            w.name.startsWith('Bending') ? `Strengthen in flexure: about ${f(AsNeed * 1e6)} mm² of tension steel is needed (bonded plates, near-surface bars or an enlarged section for an existing member), or limit the load so that the moment stays below ${f(gM.MRd / 1e3)} kN·m.` : w.name.startsWith('Shear') ? `Strengthen in shear: links of about ${f(AswS * gV.ucV * 1e6)} mm²/m are needed (closer links, external stirrups or wraps), or limit the shear to ${f(gV.VRd / 1e3)} kN.` : w.name.startsWith('Crack') ? `Reduce the steel stress under long-term load: more bars of smaller diameter, or limit the sustained moment to ${f(Mqp / 1e3)} kN·m.` : w.name.startsWith('Minimum') ? `Provide at least ${f(AsMin * 1e6)} mm² of tension steel.` : 'Increase the section depth or the compression steel so that the failing stress or ductility limit is met.',
            'EN 1992-1-1:2004 §6.1, §6.2, §7.2, §7.3; fib Bulletin 90 for externally bonded strengthening',
            `Extra steel of ${f(Math.max(dAs, 0) * 1e3)} kg per metre adds ${f(Math.max(dAs, 0) * CO2_STEEL * 1e3)} kg CO₂e per metre; enlarging the section by 100 mm instead adds ${f(S.b * 0.1)} m³ of concrete per metre (${f(S.b * 0.1 * CO2_CONC * 1e3)} kg CO₂e, screening factors).`),
        eco(ok ? 'info' : 'warn', `Strength is ${f(ucU * 100)} % used and serviceability ${f(ucSls * 100)} % used; the ratings are ${f(MRd0.MRd / 1e3)} kN·m (factored) and ${f(Mqp / 1e3)} kN·m (sustained).`,
          'Where serviceability governs, tightening crack control costs a little steel now and defers the first concrete repair by years; where strength governs an existing member, compare local strengthening with an operating load limit using the permissible moments.', `${ECON_BASIS}; service-life design per fib Model Code 2010 §7.8`),
      ],
      publish: { 'wharf.sectionMRd': MRd0.MRd, 'wharf.sectionVRd': gV.VRd, 'wharf.sectionUc': Math.max(ucU, ucSls) },
    };
  },
  verify: [
    // Singly reinforced rectangle, C30 (γc = 1.5, αcc = 1 → fcd = 20 MPa), B500 (fyd = 434.8 MPa), b = 0.30 m, h = 0.55 m,
    // 3 bars of 20 mm (As = 942.5 mm²), cover 40 mm → d = 0.55 − 0.04 − 0.01 = 0.50 m.
    // Parabola–rectangle block: resultant 0.8095·fcd·b·x acting 0.416·x below the top.
    // x = As·fyd/(0.8095·b·fcd) = 409.8e3/4.857e6 = 0.08437 m; MRd = 409.8·(0.50 − 0.416·0.08437) = 190.5 kN·m; x/d = 0.169.
    { name: 'Singly reinforced rectangular section: MRd = As·fyd·(d − 0.416 x)', source: 'EN 1992-1-1 §3.1.7 parabola–rectangle block coefficients (0.8095, 0.416): hand calculation',
      inputs: { shape: 'rect', b: 0.3, h: 0.55, conc: 'c30', nBot: 3, dBot: 20, nTop: 0, cover: 40, dLink: 0, linkLegs: 0, gammaC: 1.5, gammaS: 1.15, alphaCC: 1, nStrips: 400, combos: [{ name: 'ULS', ls: 'ULS', N: 0, M: 100, V: 10 }] },
      expect: [{ kpi: 'MRd', value: 190.5, tol: 0.005 }, { kpi: 'xd', value: 0.1687, tol: 0.01 }] },
    // Squash load of the same section: NRd = fcd·Ac + As·σs at εc2 = 2 ‰ (σs = 400 MPa < fyd): 20·0.165 MN + 942.5e-6·400 MN = 3300 + 377 = 3677 kN.
    { name: 'Squash load: NRd = fcd·Ac + As·Es·εc2', source: 'EN 1992-1-1 §6.1(5) uniform compression at εc2 (hand calculation)',
      inputs: { shape: 'rect', b: 0.3, h: 0.55, conc: 'c30', nBot: 3, dBot: 20, nTop: 0, cover: 40, dLink: 0, linkLegs: 0, gammaC: 1.5, gammaS: 1.15, alphaCC: 1, nStrips: 400, combos: [{ name: 'ULS', ls: 'ULS', N: 0, M: 100, V: 10 }] },
      expect: [{ kpi: 'NRd', value: 3677, tol: 0.005 }] },
  ],
};
