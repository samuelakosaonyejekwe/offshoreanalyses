// Anchor pile under inclined mooring load: lateral capacity by rigid rotation about the padeye, axial pull-out,
// combined-load interaction, and the embedded length required for the safety factor.
import { anchorPileCapacity } from '../../engine/geotech/anchors.js';
import { DEG, G, need, run, round, fmt, adv, check, worst, sev, linspace, steelCo2Text, vesselCo2Text } from './_common.js';
import { anchorLoadInputs, clayInputs, capacityFactorInputs, anchorLoad, LIMIT_ANCHOR } from './_anchor.js';

const clay = { key: 'soil', in: ['clay'] }, sand = { key: 'soil', in: ['sand'] };

export default {
  id: 's10-anchor-pile',
  suite: 10,
  title: 'Anchor pile capacity',
  summary: 'Checks a driven or drilled anchor pile for the lateral and uplift components of the mooring load at the padeye and finds the embedded length that meets the safety factors.',
  covers: [10],
  fidelity: 'analytical',
  standards: ['API RP 2SK (3rd ed., 2005) §6.5 and pile-anchor safety factors', 'API RP 2GEO (2011, add. 2014) §8 (axial and lateral pile capacity)', 'ISO 19901-4 (2016) §8', 'DNV-OS-E301 (2018) Ch.2 Sec.4'],
  models: ['Undrained lateral bearing p = Np·su·D (Murff & Hamilton 1993)', 'Drained lateral resistance p = 3·Kp·γ′·z·D (Broms 1964; Coulomb/Rankine passive pressure)', 'Rigid-pile rotation about the load attachment point', 'Axial pull-out: α-method in clay, K·tan δ·σ′v in sand', 'Power-law H–V interaction'],
  limits: `${LIMIT_ANCHOR} The lateral solution treats the pile as rigid, which is appropriate for L/D up to about 6–10; longer flexible piles are governed by bending and need a p–y beam-column analysis (see the foundations suite). Pile structural stress, driveability and cyclic gapping are not checked here.`,
  publishes: ['mooring.anchorCapacity'],
  searchTerms: 'anchor pile mooring lateral uplift capacity padeye',
  inputs: [
    ...anchorLoadInputs({ H: 2750, V: 1300 }),
    { key: 'D', label: 'Pile outer diameter', unit: 'm', type: 'number', default: 2.4, min: 0.4, max: 8, group: 'Pile', role: 'geometry', design: true },
    { key: 'L', label: 'Embedded length', unit: 'm', type: 'number', default: 26, min: 3, max: 120, group: 'Pile', role: 'geometry', design: true },
    { key: 'tw', label: 'Wall thickness', unit: 'mm', type: 'number', default: 50, min: 10, max: 150, group: 'Pile', role: 'geometry' },
    { key: 'steel', label: 'Steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Pile', role: 'material' },
    { key: 'zPadeye', label: 'Padeye depth below mudline', unit: 'm', type: 'number', default: 16, min: 0, max: 120, group: 'Pile', role: 'bc',
      help: 'Where the line attaches. Boundary condition of the pile: load applied at this depth, soil reaction along the shaft. Near the centroid of resistance (≈ L/2 to 2L/3) the pile translates; at the pile head it rotates and capacity drops sharply.' },
    { key: 'soil', label: 'Soil type', type: 'select', default: 'clay', options: [{ value: 'clay', label: 'Clay (undrained)' }, { value: 'sand', label: 'Sand (drained)' }], group: 'Soil', role: 'soil' },
    ...clayInputs({ su0: 5, k: 1.5 }).map((i) => ({ ...i, showIf: clay })),
    { key: 'alpha', label: 'Shaft adhesion factor α', unit: '-', type: 'number', default: 0.65, min: 0.1, max: 1, group: 'Soil', role: 'soil', showIf: clay },
    { key: 'NpConst', label: 'Lateral bearing factor Np (0 = Murff–Hamilton profile)', unit: '-', type: 'number', default: 0, min: 0, max: 14, group: 'Soil', role: 'soil', showIf: clay },
    { key: 'phiDeg', label: 'Friction angle', unit: 'deg', type: 'number', default: 32, min: 20, max: 45, group: 'Soil', role: 'soil', bind: 'site.phiDeg', showIf: sand, uncertain: { dist: 'normal', cov: 0.07 } },
    { key: 'KtanDelta', label: 'Shaft friction factor K·tan δ', unit: '-', type: 'number', default: 0.4, min: 0.1, max: 1.2, group: 'Soil', role: 'soil', showIf: sand, help: 'Tension loading: K ≈ 0.5–0.8, δ ≈ φ − 5°. 0.3–0.5 is typical for medium dense sand.' },
    { key: 'gammaSub', label: 'Submerged unit weight of soil', unit: 'kN/m³', type: 'number', default: 9, min: 2, max: 12, group: 'Soil', role: 'soil', bind: 'site.gammaSub' },
    ...capacityFactorInputs({ sf: 1.6 }),
    { key: 'sfAxial', label: 'Required axial safety factor (0 = API table)', unit: '-', type: 'number', default: 0, min: 0, max: 5, group: 'Capacity factors', role: 'criteria',
      help: 'API RP 2SK (3rd ed.) permanent intact: lateral 1.6, axial 2.0. The lateral factor is the "required safety factor" above.' },
  ],
  compute(inp, ctx) {
    const ld = anchorLoad(inp);
    const mat = ctx.material(inp.steel);
    need(mat, `Steel '${inp.steel}' is not in the materials library.`);
    const { D, L } = inp, t = inp.tw / 1e3;
    need(inp.zPadeye <= L, `The padeye (${inp.zPadeye} m) is below the pile tip (${L} m): raise the padeye or lengthen the pile.`);
    const tab = { permanent: { intact: [1.6, 2.0], damaged: [1.2, 1.5] }, mobile: { intact: [1.2, 1.5], damaged: [1.0, 1.2] } }[inp.service][inp.condition];
    const sfL = inp.sfRequired > 0 ? inp.sfRequired : tab[0], sfA = inp.sfAxial > 0 ? inp.sfAxial : tab[1];
    const fac = inp.setUp * inp.cyclic;
    const massOf = (len) => Math.PI * (D - t) * t * len * mat.rho;
    const capOf = (len, zP, angle) => {
      const W = massOf(len) * G * (1 - 1025 / mat.rho);
      const o = { soil: inp.soil, D, L: len, zPadeye: Math.min(zP, len), W, angle };
      if (inp.soil === 'clay') Object.assign(o, { su0: inp.su0 * 1e3, k: inp.suGrad * 1e3, alpha: inp.alpha, NpConst: inp.NpConst > 0 ? inp.NpConst : undefined });
      else Object.assign(o, { phi: inp.phiDeg * DEG, gamma: inp.gammaSub * 1e3, KtanDelta: inp.KtanDelta });
      const r = anchorPileCapacity(o);
      return { ...r, W, Hd: r.Hult * fac, Vd: W + (r.Vult - W) * fac };
    };
    if (inp.soil === 'clay') need(inp.su0 > 0 || inp.suGrad > 0, 'The clay has no strength: enter a mudline strength or a gradient.');
    const c = run(() => capOf(L, inp.zPadeye, ld.angle), 'Anchor pile capacity');
    const ea = L / D + 0.5, eb = L / (3 * D) + 4.5;
    const inter = (ld.H / c.Hd) ** ea + (ld.V / c.Vd) ** eb;
    const interF = (ld.H * sfL / c.Hd) ** ea + (ld.V * sfA / c.Vd) ** eb;
    const checks = [
      check('Lateral capacity', ld.H * sfL / 1e3, c.Hd / 1e3, 'kN', `API RP 2SK (3rd ed.) pile anchors, lateral factor ${sfL}`,
        'Guards against the pile being pushed through the soil sideways or rotating about a point below the padeye; exceeding it gives large, permanent pile-head displacement.'),
      check('Axial pull-out capacity', ld.V * sfA / 1e3, c.Vd / 1e3, 'kN', `API RP 2SK (3rd ed.) pile anchors, axial factor ${sfA}; API RP 2GEO §8.1`,
        'Guards against the pile being pulled out by the vertical component of the line load. Shaft friction in tension degrades under cyclic loading, hence the higher factor.'),
      { name: 'Combined H–V interaction (factored)', demand: round(interF, 4), capacity: 1, unit: '-', uc: interF, ref: 'Supachawarote, Randolph & Gourvenec (2004) interaction form',
        explain: 'Guards against failure under simultaneous lateral and uplift load: the factored load point must lie inside the interaction envelope (value ≤ 1).' },
      { name: 'Rigid-pile assumption (L/D)', demand: round(L / D, 3), capacity: 10, unit: '-', uc: L / D / 10, ref: 'Validity of the rigid-rotation solution (Broms 1964; Murff & Hamilton 1993)',
        explain: 'Not a strength check: above L/D ≈ 10 the pile bends rather than rotates, the lower soil is not mobilised and this module over-predicts lateral capacity.' },
    ];
    const g = worst(checks.slice(0, 3));
    // required length: smallest L (padeye kept at the same fraction of L) satisfying the three strength checks
    const frac = inp.zPadeye / L;
    const okAt = (len) => { try { const r = capOf(len, frac * len, ld.angle); const a = len / D + 0.5, b = len / (3 * D) + 4.5; return (ld.H * sfL / r.Hd) ** a + (ld.V * sfA / r.Vd) ** b <= 1 && ld.H * sfL <= r.Hd && ld.V * sfA <= r.Vd; } catch (e) { return false; } };
    let lo = 1.5 * D, hi = L;
    for (let i = 0; i < 30 && !okAt(hi); i++) hi *= 1.3;
    let Lreq = hi;
    if (okAt(hi)) { for (let i = 0; i < 22; i++) { const mid = 0.5 * (lo + hi); if (okAt(mid)) hi = mid; else lo = mid; } Lreq = hi; } else Lreq = 999;
    const zs = linspace(0, L, 15);
    const hz = zs.map((z) => capOf(L, z, ld.angle).Hd);
    const Ls = linspace(Math.max(2 * D, 0.4 * L), 1.6 * L, 12);
    const capL = Ls.map((len) => capOf(len, frac * len, ld.angle));
    const steel_t = massOf(L) / 1e3;
    return {
      kpis: [
        { key: 'Hult', label: 'Lateral capacity at the padeye', value: round(c.Hult / 1e3, 5), unit: 'kN', explain: 'Ultimate horizontal load at the padeye depth from full soil reaction on the rotating pile, before set-up and cyclic factors.' },
        { key: 'Vfriction', label: 'Shaft friction in tension', value: round((c.Vult - c.W) / 1e3, 5), unit: 'kN', explain: 'Pull-out resistance from external shaft friction alone (no reverse end bearing counted).' },
        { key: 'Vult', label: 'Axial pull-out capacity', value: round(c.Vult / 1e3, 5), unit: 'kN', explain: 'Shaft friction plus submerged pile weight.' },
        { key: 'fsLat', label: 'Lateral factor of safety', value: round(c.Hd / Math.max(ld.H, 1e-9), 4), unit: '-', limit: sfL, sense: 'min', explain: 'Design lateral resistance over the horizontal load component.' },
        { key: 'fsAx', label: 'Axial factor of safety', value: round(Math.min(c.Vd / Math.max(ld.V, 1e-9), 999), 4), unit: '-', limit: sfA, sense: 'min', explain: 'Design pull-out resistance over the vertical load component (999 when there is no uplift).' },
        { key: 'interaction', label: 'H–V interaction, unfactored load', value: round(inter, 4), unit: '-', explain: 'Position of the unfactored load in the interaction envelope; 1.0 is failure.' },
        { key: 'Tult', label: 'Inclined capacity at the load angle', value: round((c.Tult ?? 0) / 1e3, 5), unit: 'kN', explain: 'Ultimate line tension for the actual inclination from the interaction envelope (without set-up and cyclic factors).' },
        { key: 'zOpt', label: 'Optimal padeye depth', value: round(c.zOptimal, 4), unit: 'm', explain: 'Centroid of the lateral soil resistance: attaching here gives pure translation and the largest lateral capacity.' },
        { key: 'z0', label: 'Rotation depth', value: round(Number.isFinite(c.z0) ? c.z0 : 0, 4), unit: 'm', explain: 'Depth of the point about which the rigid pile rotates at failure; 0 is reported for pure translation.' },
        { key: 'Lreq', label: 'Required embedded length', value: round(Lreq, 4), unit: 'm', limit: L, sense: 'max', explain: 'Shortest pile (padeye at the same relative depth) that satisfies the lateral, axial and combined checks. 999 means no practical length works at this diameter.' },
        { key: 'steel', label: 'Pile steel mass', value: round(steel_t, 4), unit: 't', explain: 'Embedded pile mass, the main cost driver with the installation spread.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Lateral capacity vs padeye depth', x: { label: 'Padeye depth below mudline', unit: 'm' }, y: { label: 'Design lateral capacity', unit: 'kN' },
          series: [{ name: 'Capacity', x: zs.map((v) => round(v)), y: hz.map((v) => round(v / 1e3)) }, { name: 'Factored demand', x: [0, L], y: [ld.H * sfL / 1e3, ld.H * sfL / 1e3], dash: true }, { name: 'Chosen padeye', x: [inp.zPadeye], y: [round(c.Hd / 1e3)] }],
          explain: 'How strongly the lateral capacity depends on where the line is attached. The peak is the optimal padeye depth; capacity at the pile head is typically a third of it.' },
        { type: 'line', title: 'Capacity vs embedded length', x: { label: 'Embedded length', unit: 'm' }, y: { label: 'Design capacity', unit: 'kN' },
          series: [{ name: 'Lateral', x: Ls.map((v) => round(v)), y: capL.map((r) => round(r.Hd / 1e3)) }, { name: 'Axial pull-out', x: Ls.map((v) => round(v)), y: capL.map((r) => round(r.Vd / 1e3)) },
            { name: 'Factored lateral demand', x: [Ls[0], Ls[Ls.length - 1]], y: [ld.H * sfL / 1e3, ld.H * sfL / 1e3], dash: true }, { name: 'Factored axial demand', x: [Ls[0], Ls[Ls.length - 1]], y: [ld.V * sfA / 1e3, ld.V * sfA / 1e3], dash: true }],
          explain: 'Lateral and pull-out capacity as the pile is lengthened, against the factored demands. The length where both solid curves clear their dashed lines is the required embedment.' },
        { type: 'line', title: 'H–V interaction diagram', x: { label: 'Horizontal load', unit: 'kN' }, y: { label: 'Vertical load', unit: 'kN' },
          series: [{ name: 'Capacity envelope', x: linspace(0, 1, 30).map((u) => round(c.Hd * u / 1e3)), y: linspace(0, 1, 30).map((u) => round(c.Vd * Math.max(0, 1 - u ** ea) ** (1 / eb) / 1e3)) }, { name: 'Design load', x: [round(ld.H / 1e3)], y: [round(ld.V / 1e3)] }, { name: 'Factored load', x: [round(ld.H * sfL / 1e3)], y: [round(ld.V * sfA / 1e3)] }],
          explain: 'Interaction of lateral and uplift capacity with the unfactored and factored load points. Both points must lie inside the envelope.' },
      ],
      tables: [
        { title: 'Capacity summary', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['Lateral capacity (static)', c.Hult / 1e3, 'kN'], ['Lateral capacity (design)', c.Hd / 1e3, 'kN'], ['Pull-out capacity (static)', c.Vult / 1e3, 'kN'], ['Pull-out capacity (design)', c.Vd / 1e3, 'kN'], ['Submerged pile weight', c.W / 1e3, 'kN'], ['Horizontal demand', ld.H / 1e3, 'kN'], ['Vertical demand', ld.V / 1e3, 'kN'], ['Load angle at padeye', ld.angleDeg, 'deg'], ['L/D', L / D, '-']].map(([q, v, u]) => ({ q, v: round(v, 5), u })),
          explain: 'Static and design capacities next to the load components. Design values include the set-up and cyclic factors; weight is not factored.' },
      ],
      notes: [`${inp.soil === 'clay' ? 'Undrained clay' : 'Drained sand'} response; pile L/D = ${fmt(L / D)}, padeye at ${fmt(100 * frac)} % of the embedded length (optimum ${fmt(c.zOptimal)} m).`,
        `Safety factors: lateral ${sfL}, axial ${sfA} (${inp.sfRequired > 0 || inp.sfAxial > 0 ? 'user values where entered' : 'API RP 2SK table'}).`,
        L / D > 10 ? 'L/D exceeds 10: treat the lateral capacity as an upper bound and run a p–y analysis.' : 'Rigid-pile assumption is within its usual range.'],
      advice: [
        adv('engineering', sev(g.uc), `Lateral factor ${fmt(c.Hd / Math.max(ld.H, 1e-9))} (need ${sfL}), axial factor ${fmt(Math.min(c.Vd / Math.max(ld.V, 1e-9), 999))} (need ${sfA}), factored interaction ${fmt(interF)}; governing: ${g.name} (${fmt(g.uc)}).`,
          g.uc > 1 ? `Lengthen the pile to about ${fmt(Lreq)} m at this diameter${Math.abs(inp.zPadeye - c.zOptimal) > 0.15 * L ? `, and move the padeye from ${inp.zPadeye} m towards the optimum at ${fmt(c.zOptimal)} m` : ''}.`
            : `The pile could be as short as ${fmt(Lreq)} m. Confirm driveability and pile stresses at the padeye, and check lateral response with p–y curves if L/D approaches 10.`,
          'API RP 2SK (3rd ed., 2005) §6.5; API RP 2GEO (2011) §8; ISO 19901-4 (2016) §8.'),
        adv('economic', 'info', `Pile steel ${fmt(steel_t)} t; required length ${fmt(Lreq)} m against ${L} m provided.`,
          Lreq < 0.95 * L ? `Shortening to ${fmt(Lreq)} m saves about ${fmt(steel_t * (1 - Lreq / L))} t per pile and driving time; keep a driving tolerance of 1–2 m.` : Lreq > L ? `The extra ${fmt(Lreq - L)} m adds ${fmt(massOf(Lreq) / 1e3 - steel_t)} t per pile; compare with a larger diameter, which raises lateral capacity linearly at lower driving risk.` : 'Length is matched to the requirement; standardise it across anchor locations.',
          'Quantity-based Class 4 estimate per AACE RP 18R-97; pile installation planned per DNV-ST-N001 (2021).',
          `${steelCo2Text(Lreq < 0.95 * L ? steel_t * (1 - Lreq / L) : steel_t)} ${vesselCo2Text(1.5)} Driven piles are normally cut below mudline and left in place at decommissioning; underwater piling noise needs mitigation where marine mammals are present.`),
      ],
      publish: { 'mooring.anchorCapacity': (c.Tult ?? Math.min(c.Hd, c.Vd)) * fac },
      quantities: { steel_t: round(steel_t), vessel_days: { anchor_pile_installation: 1.5 } },
    };
  },
  verify: [
    { // Uniform clay su = 50 kPa, constant Np = 9, load at the centroid (mid-depth): pure translation,
      // H = Np·su·D·L = 9·50·2·20 = 18 000 kN. Shaft friction α·su·π·D·L = 0.65·50·π·2·20 = 4084 kN.
      name: 'Clay: translational lateral capacity and α-method shaft friction',
      source: 'Randolph & Gourvenec (2011) Offshore Geotechnical Engineering §5.2, §9.3; API RP 2GEO §8.1 (α-method)',
      inputs: { soil: 'clay', D: 2, L: 20, zPadeye: 10, su0: 50, suGrad: 0, NpConst: 9, alpha: 0.65 },
      expect: [{ kpi: 'Hult', value: 18000, tol: 0.005 }, { kpi: 'Vfriction', value: 4084, tol: 0.003 }, { kpi: 'zOpt', value: 10, tol: 0.01 }] },
    { // Sand: tension shaft friction ½·K tan δ·γ′·L²·π·D = 0.5·0.4·10·30²·π·2 = 11 310 kN.
      name: 'Sand: shaft friction in tension',
      source: 'API RP 2GEO (2011) §8.1.4 form f = K·σ′v·tan δ integrated over a uniform sand',
      inputs: { soil: 'sand', D: 2, L: 30, zPadeye: 20, KtanDelta: 0.4, gammaSub: 10, phiDeg: 35 },
      expect: [{ kpi: 'Vfriction', value: 11310, tol: 0.003 }] },
  ],
};
