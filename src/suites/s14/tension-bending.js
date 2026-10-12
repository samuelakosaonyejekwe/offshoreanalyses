// Axial tension and bending of the cross-section: axial failure and allowable tension, stick/slip bending
// stiffness, wire stress round the circumference under combined tension and curvature, the tension–curvature
// capacity envelope and the friction hysteresis loop.
import { axisymmetricResponse, axialFailureTension, bendingProperties } from '../../engine/xsection/flexible.js';
import { wireStressAtPosition, momentCurvatureLoop } from '../../engine/xsection/helix.js';
import { stackInputs, stackFrom, loadInputs, loadsFrom, muInput, massSplit, capacityOf, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, linspace, sev,
  LIMIT_AXI, GENERIC, STD_17J, STD_17B, CONFIRM, CO2_STEEL, DEG } from './_common.js';

export default {
  id: 's14-tension-bending',
  suite: 14,
  title: 'Axial tension and bending',
  summary: 'Finds the tension capacity, the bending stiffness before and after armour slip and the wire stress under combined tension and curvature, with the tension–curvature envelope.',
  covers: [3],
  fidelity: 'analytical',
  standards: ['API Spec 17J (4th ed.) §5.3, Table 6', 'API RP 17B (5th ed.) §8'],
  models: ['Helical wire bending and torsion (loxodromic / geodesic wire paths)', 'Interlayer stick–slip friction model (Féret & Bournazel 1987; Sævik 2011)', 'Coupled axial–bending stiffness', 'Nonlinear moment–curvature hysteresis (bilinear)', 'Elastic scaling to the wire tensile strength'],
  limits: LIMIT_AXI + ' Failure tension is elastic scaling of the wire stress to the tensile strength; contact pressure is taken uniform round the circumference and constant through the bending cycle. ' + GENERIC,
  publishes: ['flexible.failureTension', 'flexible.allowableTension'],
  searchTerms: 'flexible pipe bending stiffness stick slip tensile armour moment curvature',
  inputs: [
    ...stackInputs(),
    ...loadInputs(),
    { key: 'bendRadius', label: 'Imposed bend radius', unit: 'm', type: 'number', default: 5, min: 0.3, max: 5000, group: 'Loads and pressure boundaries', role: 'bc', design: true,
      help: 'Radius of curvature applied together with the tension (prescribed curvature = 1/radius), e.g. the smallest radius from global analysis.' },
    muInput(),
    { key: 'wirePath', label: 'Wire path in bending', type: 'select', default: 'loxodromic', options: [{ value: 'loxodromic', label: 'Loxodromic (wire keeps its lay angle)' }, { value: 'geodesic', label: 'Geodesic (wire slides sideways)' }], group: 'Interlayer contact', role: 'control',
      help: 'Loxodromic is the usual assumption for service loads; geodesic is the limit for large curvature with free transverse slip.' },
    { key: 'fTa', label: 'Tensile armour utilisation factor', unit: '-', type: 'number', default: 0.67, min: 0.3, max: 1, group: 'Acceptance criteria', role: 'criteria', help: 'API 17J Table 6: 0.67 normal, 0.85 abnormal, 0.91 hydrotest. ' + CONFIRM },
  ],
  compute(inp, ctx) {
    const S = stackFrom(inp, ctx);
    need(S.idx.ta.length > 0, 'The stack has no tensile armour layer: axial tension capacity cannot be assessed.');
    const loads = loadsFrom(inp);
    const kappa = 1 / inp.bendRadius;
    const r = run(() => axisymmetricResponse(S.stack, loads), 'Axisymmetric solution');
    const r0 = run(() => axisymmetricResponse(S.stack, { ...loads, tension: loads.tension + 1e5 }), 'Axisymmetric solution');
    const bend = run(() => bendingProperties(S.stack, { response: r, mu: inp.mu, curve: inp.wirePath }), 'Bending properties');
    const fail = run(() => axialFailureTension(S.stack, { pInt: loads.pInt, pExt: loads.pExt, basis: 'uts' }), 'Axial failure tension');

    // wire stress is linear in tension: σ = s0 + s1·T per layer
    const layers = S.idx.ta.map((i, k) => {
      const L = S.stack.layers[i], b = bend.layers[k];
      const s1 = (r0.layers[i].wireStress - r.layers[i].wireStress) / 1e5;
      const sAx = r.layers[i].wireStress;
      const n = L.n ?? 1, A = L.nA / n, t = L.wireThickness ?? L.t, w = L.wireWidth ?? 3 * t;
      const f = b.frictionPerLength ?? 0;
      const at = (kap, psi) => wireStressAtPosition({ E: L.E, A, R: L.R, alpha: L.alpha, thickness: t, width: w, kappa: kap, psi, frictionPerLength: f, sigmaAxial: sAx, curve: inp.wirePath });
      const psis = linspace(0, 2 * Math.PI, 49);
      const prof = psis.map((p) => at(kappa, p));
      const sMax = Math.max(...prof.map((q) => Math.max(Math.abs(q.sigmaMaxCorner), Math.abs(q.sigmaMinCorner))));
      const extra = (kap) => Math.max(...[0, Math.PI / 2, Math.PI].map((p) => { const q = at(kap, p); return Math.max(q.sigmaMaxCorner, -q.sigmaMinCorner) - (q.sigmaAxial); })) ;
      return { L, b, s1, sAx, cap: inp.fTa * capacityOf(L), psis, prof, sMax, extra, state: at(kappa, 0).state };
    });
    const allowT = Math.min(...layers.map((y) => (y.s1 > 0 ? loads.tension + (y.cap - y.sAx) / y.s1 : Infinity)));
    need(Number.isFinite(allowT), 'Tensile armour stress does not increase with tension for this stack — check lay angles and strengths.');
    const gov = layers.reduce((a, c) => (c.sMax / c.cap > a.sMax / a.cap ? c : a), layers[0]);
    const ucT = Math.max(...layers.map((y) => Math.abs(y.sAx) / y.cap));
    const ucComb = gov.sMax / gov.cap;
    const kY = bend.frictionMoment / Math.max(bend.EIstick - bend.EIslip, 1e-9);
    const loop = kappa > kY * 1.0001 ? momentCurvatureLoop({ EIstick: bend.EIstick, EIslip: bend.EIslip, frictionMoment: bend.frictionMoment, kappaAmp: kappa, nPoints: 21 }) : null;
    const moment = Math.min(bend.EIstick * kappa, bend.frictionMoment + bend.EIslip * kappa);

    // capacity envelope: allowable tension at each curvature, keeping friction at the present contact state
    const kMax = Math.max(2.5 * kappa, 0.4);
    const ks = linspace(0, kMax, 41);
    const env = ks.map((k) => Math.max(0, Math.min(...layers.map((y) => (y.s1 > 0 ? loads.tension + (y.cap - y.sAx - y.extra(k)) / y.s1 : Infinity)))));

    const checks = [
      chk('Tensile armour: axial stress under tension', Math.max(...layers.map((y) => Math.abs(y.sAx))) / 1e6, layers[0].cap / 1e6, 'MPa', `${STD_17J} (factor ${inp.fTa})`,
        'Guards against tensile failure of the armour wires under straight-pipe tension and pressure. Exceeding it erodes the margin to wire rupture.'),
      chk(`${gov.L.name}: combined tension + bending stress`, gov.sMax / 1e6, gov.cap / 1e6, 'MPa', `${STD_17J}; local bending per Sævik (2011)`,
        'Guards against overstressing the wire corner where axial, friction and local bending stresses add at the imposed curvature. Persistent exceedance leads to wire yielding and early fatigue cracking.'),
      chk('Effective tension vs axial failure tension', loads.tension / 1e3, fail.tension / 1e3, 'kN', `${STD_17B} (tensile failure of the tensile armour)`,
        'Compares the applied tension with the tension at which the most loaded wire reaches its tensile strength; the ratio is the true margin to parting the pipe.'),
    ];
    const ms = massSplit(S.stack);
    const taSteel = S.idx.ta.reduce((a, i) => a + S.stack.layers[i].rho * S.stack.layers[i].nA / Math.cos(S.stack.layers[i].alpha), 0); // kg/m
    const uc = Math.max(ucT, ucComb);
    return {
      kpis: [
        kpi('failTension', 'Axial failure tension', fail.tension / 1e3, 'kN', `Effective tension at which a wire of ${fail.layer} reaches its tensile strength (elastic scaling). The pipe parts near this load.`),
        kpi('allowTension', 'Allowable tension (straight pipe)', allowT / 1e3, 'kN', `Tension at which the most loaded tensile armour wire reaches ${inp.fTa} × structural capacity with the present pressures.`),
        kpi('tensionUc', 'Tension utilisation', ucT, '-', 'Tensile armour stress over the permissible stress for the straight pipe. Above 1.0 the tension criterion is not met.', 1, 'max'),
        kpi('combinedUc', 'Combined tension + bending utilisation', ucComb, '-', 'Peak wire-corner stress at the imposed bend radius over the permissible stress. Above 1.0 reduce curvature or tension.', 1, 'max'),
        kpi('EIstick', 'Bending stiffness before slip', bend.EIstick / 1e3, 'kN·m²', 'Stiffness while friction locks the wires to the pipe; closed form Σ n·E·A·R²·cos³α/2 plus the polymer layers.'),
        kpi('EIslip', 'Bending stiffness after slip', bend.EIslip / 1e3, 'kN·m²', 'Stiffness once the wires slide: polymer layers plus the wires bending about their own axes. Typically 1–3 % of the stick value.'),
        kpi('kappaCr', 'Slip-onset curvature', bend.kappaCrMin ?? 0, '1/m', 'Curvature at which the first armour layer starts to slide. Beyond it the pipe softens and dissipates energy every cycle.'),
        kpi('frictionMoment', 'Friction moment', bend.frictionMoment / 1e3, 'kN·m', 'Moment held by interlayer friction at full slip; grows with tension and pressure because they raise contact pressure.'),
        kpi('moment', 'Bending moment at the imposed radius', moment / 1e3, 'kN·m', 'Moment needed to hold the pipe at the imposed bend radius (bilinear stick–slip law).'),
        kpi('axialStrain', 'Axial strain', r.axialStrain * 100, '%', 'Extension per unit length under the applied tension and pressures.'),
        kpi('damping', 'Equivalent hysteretic damping', loop ? loop.dampingRatio : 0, '-', 'Equivalent viscous damping ratio of a bending cycle of this amplitude; zero while the wires stick.'),
      ],
      checks,
      plots: [
        plot('line', 'Tension–curvature capacity envelope', 'Curvature', '1/m', 'Allowable effective tension', 'kN', [
          { name: 'Allowable tension', x: ks.map((k) => round(k)), y: env.map((v) => round(v / 1e3)) },
          { name: 'Operating point', x: [round(kappa)], y: [round(loads.tension / 1e3)] },
        ], 'Allowable tension falls as curvature adds bending and friction stress to the wires. The operating point must lie below the curve; the knee marks full slip.'),
        plot('line', 'Wire stress round the circumference', 'Angle from extrados', 'deg', 'Stress', 'MPa', [
          ...layers.map((y) => ({ name: `${y.L.name} max corner`, x: y.psis.map((p) => round(p / DEG)), y: y.prof.map((q) => round(q.sigmaMaxCorner / 1e6)) })),
          { name: 'Permissible', x: [0, 360], y: [round(gov.cap / 1e6), round(gov.cap / 1e6)], dash: true },
        ], 'Peak corner stress of a wire as it winds round the pipe at the imposed curvature. Friction stress peaks at the extrados and intrados (0° and 180°), transverse wire bending at the neutral axis (90°, 270°).'),
        plot('line', 'Moment–curvature hysteresis loop', 'Curvature', '1/m', 'Moment', 'kN·m', loop
          ? [{ name: 'Steady cycle', x: loop.kappa.map((k) => round(k)), y: loop.moment.map((v) => round(v / 1e3)) }]
          : [{ name: 'Stick (no slip)', x: [-round(kappa), round(kappa)], y: [-round(bend.EIstick * kappa / 1e3), round(bend.EIstick * kappa / 1e3)] }],
        'Bending moment against curvature for a steady cycle at the imposed amplitude. The enclosed area is the energy dissipated by friction per cycle and per metre of pipe.'),
      ],
      tables: [tbl('Tensile armour layers', [col('layer', 'Layer'), col('axial', 'Axial wire stress', 'MPa'), col('peak', 'Peak corner stress', 'MPa'), col('allow', 'Permissible', 'MPa'), col('kcr', 'Slip-onset curvature', '1/m'), col('fric', 'Max friction stress', 'MPa'), col('state', 'State at imposed curvature')],
        layers.map((y) => ({ layer: y.L.name, axial: round(y.sAx / 1e6), peak: round(y.sMax / 1e6), allow: round(y.cap / 1e6), kcr: round(y.b.kappaCr ?? 0), fric: round((y.b.sigmaFrictionMax ?? 0) / 1e6), state: y.state })),
        'One row per tensile armour layer: stress from tension and pressure alone, the peak with bending added, and the slip state at the imposed curvature.')],
      notes: [`${gov.L.name} governs the combined check (${fmt(ucComb)}); the wires are in ${gov.state} at a bend radius of ${fmt(inp.bendRadius)} m.`,
        `Slip starts at a radius of ${bend.kappaCrMin ? fmt(1 / bend.kappaCrMin) : '∞'} m; stiffness drops from ${fmt(bend.EIstick / 1e3)} to ${fmt(bend.EIslip / 1e3)} kN·m².`, GENERIC],
      advice: [
        adv('engineering', sev(uc), `Tension utilisation ${fmt(ucT)}, combined tension + bending ${fmt(ucComb)}; failure tension ${fmt(fail.tension / 1e3)} kN against ${fmt(loads.tension / 1e3)} kN applied.`,
          uc > 1 ? `Limit the tension to ${fmt(allowT / 1e3)} kN at this pressure, or increase the bend radius; a ${fmt((uc - 1) * 100)} % larger armour area would also close the gap.` : `Margin of ${fmt((1 - uc) * 100)} %. Use the slip stiffness (${fmt(bend.EIslip / 1e3)} kN·m²) in global analysis and confirm the friction coefficient from a bending test.`,
          `${STD_17J}; ${STD_17B}; Sævik (2011) Marine Structures 24. Utilisation factor is a recalled default — confirm.`),
        adv('economic', 'info', `The tensile armour is ${fmt(taSteel)} t of steel per kilometre (${fmt(taSteel / Math.max(ms.steel, 1e-9) * 100)} % of the steel in the pipe) and is ${fmt(ucT * 100)} % utilised in tension.`,
          ucT < 0.5 ? `About ${fmt(taSteel * (1 - ucT / 0.9))} t/km of armour is not needed for this tension: a lighter armour or a single pipe design covering deeper water would make better use of it.` : 'The armour is well used; the cheaper lever is top tension (buoyancy distribution, configuration) rather than more steel.',
          'Material take-off per AACE Class 4 estimate practice; whole-life costing per ISO 15686-5',
          `Tensile armour steel embodies about ${fmt(taSteel * CO2_STEEL)} t CO₂e per kilometre (screening factor ${CO2_STEEL} t CO₂e/t). A higher-strength wire reaches the same capacity with less mass, provided it remains qualified for the annulus environment.`),
      ],
      publish: { 'flexible.failureTension': fail.tension, 'flexible.allowableTension': allowT },
    };
  },
  verify: [{
    // Stick bending stiffness, closed form: polymer tube E·π(ro⁴ − ri⁴)/4 plus, per armour layer, (nA)·E·R²·cos³α/2
    // with nA = fill·2πR·t·cosα. Tube: bore 200, t = 5 mm, E = 0.35 GPa → 0.35e9·π/4·(0.105⁴ − 0.1⁴) = 5 924 N·m².
    // Armour 1: R = 0.107, t = 4 mm, fill 0.9, α = 35° → nA = 1.9828e-3 m², EI = 1.9828e-3·207e9·0.107²·cos³35°/2 = 1.2915e6.
    // Armour 2: R = 0.111 → nA = 2.0569e-3 m², EI = 1.4418e6. Total = 2.7392e6 N·m² = 2 739.2 kN·m² (wire width unknown,
    // so the wires' own bending stiffness is not added).
    name: 'Stick bending stiffness of two helical layers on a polymer tube',
    source: 'Witz & Tan (1992) closed form Σ n·E·A·R²·cos³α/2; hand calculation',
    inputs: { boreD: 200, layers: [
      { name: 'Tube', type: 'pressureSheath', t: 5, alphaDeg: 0, lay: 0, fill: 1, wireW: 0, E: 0.35, yield: 0, uts: 0, lock: 0 },
      { name: 'Armour 1', type: 'tensileArmour', t: 4, alphaDeg: 35, lay: 1, fill: 0.9, wireW: 0, E: 207, yield: 1100, uts: 1400, lock: 0 },
      { name: 'Armour 2', type: 'tensileArmour', t: 4, alphaDeg: 35, lay: -1, fill: 0.9, wireW: 0, E: 207, yield: 1100, uts: 1400, lock: 0 }] },
    expect: [{ kpi: 'EIstick', value: 2739.2, tol: 0.003 }],
  }],
};
