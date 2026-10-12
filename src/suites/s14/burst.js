// Pressure containment and burst: burst pressure by the pressure-armour (hoop) and tensile-armour (axial) modes,
// utilisation at design pressure and at factory hydrotest, and the effect of tension on burst.
import { axisymmetricResponse, layerUtilisations, burstPressure } from '../../engine/xsection/flexible.js';
import { lameStresses } from '../../engine/xsection/tubes.js';
import { stackInputs, stackFrom, loadInputs, loadsFrom, factorInputs, factorsFrom, massSplit, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, linspace, sev, worst,
  LIMIT_AXI, GENERIC, STD_17J, STD_17B, CONFIRM, CO2_STEEL, BAR } from './_common.js';

export default {
  id: 's14-burst',
  suite: 14,
  title: 'Pressure containment and burst',
  summary: 'Calculates the burst pressure of the armour layers and the utilisation at design and hydrotest pressure.',
  covers: [5],
  fidelity: 'analytical',
  standards: ['API Spec 17J (4th ed.) §5.3, Table 6, §9.5 (hydrostatic test)', 'API RP 17B (5th ed.) §8 (burst)'],
  models: ['Internal pressure burst criterion (elastic scaling to wire tensile strength)', 'Pressure armour hoop equilibrium', 'Tensile armour end-cap equilibrium', 'Lamé thick-cylinder stresses in the pressure sheath', 'Multilayer axisymmetric analysis'],
  limits: LIMIT_AXI + ' Burst is the elastic wire stress scaled to the tensile strength — no plastic redistribution between layers, so it is a first-failure estimate, not a predicted test burst. Utilisation and test factors are recalled defaults. ' + GENERIC,
  publishes: ['flexible.burstPressure'],
  searchTerms: 'flexible pipe burst pressure armour hydrostatic test',
  inputs: [
    ...stackInputs(),
    ...loadInputs({ Te: 300 }).map((i) => (i.key === 'pInt' ? { ...i, label: 'Design pressure', design: true, help: 'Maximum bore pressure in operation (absolute), acting on the pressure sheath.' } : i)),
    { key: 'fatFactor', label: 'Hydrotest pressure / design pressure', unit: '-', type: 'number', default: 1.5, min: 1, max: 2, group: 'Pressure history', role: 'criteria',
      help: 'Factory acceptance test factor: commonly 1.5 for static flowlines and 1.3 for dynamic risers. ' + CONFIRM },
    ...factorInputs().filter((i) => i.key !== 'condition' && !/A$/.test(i.key)),
  ],
  compute(inp, ctx) {
    const S = stackFrom(inp, ctx);
    need(S.idx.pa >= 0 || S.idx.ta.length > 0, 'The stack has no armour layer: burst cannot be assessed.');
    const loads = loadsFrom(inp);
    const factors = factorsFrom({ ...inp, fTaA: inp.fTaN, fPaA: inp.fPaN });
    const burst = run(() => burstPressure(S.stack, { pExt: loads.pExt, tension: loads.tension, basis: 'uts' }), 'Burst pressure');
    const burstY = run(() => burstPressure(S.stack, { pExt: loads.pExt, tension: loads.tension, basis: 'yield' }), 'First-yield pressure');
    const rOp = run(() => axisymmetricResponse(S.stack, loads), 'Operating solution');
    const uOp = layerUtilisations(S.stack, rOp, { condition: 'normal', factors });
    const pFat = inp.fatFactor * loads.pInt;
    const rFat = run(() => axisymmetricResponse(S.stack, { pInt: pFat, pExt: 0, tension: 0, annulus: 'dry', pAnnulus: 0 }), 'Hydrotest solution'); // test on the factory floor: no external pressure
    const uFat = layerUtilisations(S.stack, rFat, { condition: 'hydrotest', factors });
    const mk = (u, cond, why) => u.layers.filter((x) => x.utilisation !== null).map((x) => chk(`${x.name}: ${cond}`, x.stress / 1e6, x.factor * x.capacity / 1e6, 'MPa', `${STD_17J} (factor ${x.factor})`, why));
    const checks = [
      ...mk(uOp, 'design pressure', 'Guards against burst in service: the wire stress at design pressure must stay below the permissible fraction of the structural capacity.'),
      ...mk(uFat, 'hydrotest', 'Guards against yielding the armour during the factory pressure test, which would leave permanent set in the pressure armour.'),
    ];
    need(checks.length > 0, 'No armour layer carries strength data (yield or tensile strength).');
    const gov = worst(checks);
    // pressure at which the governing normal-operation utilisation reaches 1 (stress linear in pressure about the operating state)
    const rHi = axisymmetricResponse(S.stack, { ...loads, pInt: loads.pInt + 1e6 });
    let pAllow = Infinity;
    uOp.layers.forEach((x) => {
      if (x.utilisation === null) return;
      const i = S.stack.layers.findIndex((L) => L.name === x.name);
      const ds = (rHi.layers[i].wireStress - rOp.layers[i].wireStress) / 1e6;
      if (ds > 0) pAllow = Math.min(pAllow, loads.pInt + (x.factor * x.capacity - rOp.layers[i].wireStress) / ds);
    });
    need(Number.isFinite(pAllow), 'Armour stress does not rise with pressure for this stack — check the layer table.');
    // pressure sheath as a thick cylinder between bore pressure and the contact pressure of the armour (Lamé)
    let sheath = null;
    if (S.idx.sheath >= 0) {
      const L = S.stack.layers[S.idx.sheath];
      sheath = lameStresses({ pInt: rOp.layers[S.idx.sheath].pIn, pExt: rOp.layers[S.idx.sheath].pOut, ri: L.rIn, ro: L.rOut, r: L.rIn, ends: 'open' });
    }
    const Ts = linspace(0, Math.max(2 * Math.abs(loads.tension), 1e6), 9);
    const bT = Ts.map((T) => burstPressure(S.stack, { pExt: loads.pExt, tension: T, basis: 'uts' }));
    const ps = linspace(0, burst.burst, 9);
    const sols = ps.map((p) => axisymmetricResponse(S.stack, { ...loads, pInt: p }));
    const lay = [S.idx.pa, ...S.idx.ta].filter((i) => i >= 0);
    const ms = massSplit(S.stack);
    const ratio = burst.burst / Math.max(loads.pInt, 1);
    const paL = S.idx.pa >= 0 ? S.stack.layers[S.idx.pa] : null;
    const paSteel = paL ? paL.rho * paL.nA / Math.cos(paL.alpha) : 0;
    return {
      kpis: [
        kpi('pBurst', 'Burst pressure', burst.burst / BAR, 'bar', `Internal pressure at which the ${burst.governing === 'pressureArmour' ? 'pressure armour (hoop mode)' : 'tensile armour (axial mode)'} wire reaches its tensile strength, by elastic scaling.`),
        kpi('pBurstPA', 'Burst pressure, pressure-armour mode', (burst.pressureArmour ?? burst.burst) / BAR, 'bar', 'Pressure at which the hoop-carrying pressure armour reaches its tensile strength.'),
        kpi('pBurstTA', 'Burst pressure, tensile-armour mode', (burst.tensileArmour ?? burst.burst) / BAR, 'bar', 'Pressure at which the end-cap load brings the tensile armour to its tensile strength.'),
        kpi('pYield', 'First-yield pressure', burstY.burst / BAR, 'bar', 'Pressure at which the first armour wire reaches yield; permanent deformation starts here.'),
        kpi('burstRatio', 'Burst pressure / design pressure', ratio, '-', 'Ratio of calculated burst to design pressure. Qualified designs usually show well above 2; the acceptance itself is by the layer utilisations.'),
        kpi('ucDesign', 'Utilisation at design pressure', Math.max(...uOp.layers.map((x) => x.utilisation ?? 0)), '-', 'Highest armour utilisation at design pressure with the normal-operation factors. Above 1.0 the pipe is under-designed for this pressure.', 1, 'max'),
        kpi('ucFat', 'Utilisation at hydrotest', Math.max(...uFat.layers.map((x) => x.utilisation ?? 0)), '-', 'Highest armour utilisation during the factory test (no external pressure) with the hydrotest factors.', 1, 'max'),
        kpi('pAllow', 'Allowable design pressure', pAllow / BAR, 'bar', 'Bore pressure at which the governing armour layer just reaches its normal-operation permissible stress at this depth and tension.'),
        kpi('sheathHoop', 'Pressure sheath hoop stress at the bore', sheath ? sheath.hoop / 1e6 : 0, 'MPa', 'Lamé hoop stress in the sheath between the bore pressure and the supporting armour; compressive values show the sheath is fully supported.'),
      ],
      checks,
      plots: [
        plot('line', 'Armour stress vs internal pressure', 'Internal pressure', 'bar', 'Wire stress', 'MPa',
          lay.map((i) => ({ name: S.stack.layers[i].name, x: ps.map((p) => round(p / BAR)), y: sols.map((q) => round(q.layers[i].wireStress / 1e6)) })),
          'Wire stress in each armour layer as the bore pressure rises to burst. The layer whose line reaches its tensile strength first sets the burst mode.'),
        plot('line', 'Burst pressure vs effective tension', 'Effective tension', 'kN', 'Burst pressure', 'bar', [
          { name: 'Pressure-armour mode', x: Ts.map((T) => round(T / 1e3)), y: bT.map((b) => round((b.pressureArmour ?? b.burst) / BAR)) },
          { name: 'Tensile-armour mode', x: Ts.map((T) => round(T / 1e3)), y: bT.map((b) => round((b.tensileArmour ?? b.burst) / BAR)) },
          { name: 'Design pressure', x: [0, round(Ts[Ts.length - 1] / 1e3)], y: [inp.pInt, inp.pInt], dash: true },
        ], 'Interaction of pressure and tension. Tension uses up tensile-armour capacity, so the axial burst mode falls with tension; the lower line governs.'),
      ],
      tables: [tbl('Armour utilisation', [col('layer', 'Layer'), col('sOp', 'Stress at design pressure', 'MPa'), col('uOp', 'Utilisation, design', '-'), col('sFat', 'Stress at hydrotest', 'MPa'), col('uFat', 'Utilisation, hydrotest', '-')],
        uOp.layers.filter((x) => x.utilisation !== null).map((x, k) => { const f = uFat.layers.filter((y) => y.utilisation !== null)[k]; return { layer: x.name, sOp: round(x.stress / 1e6), uOp: round(x.utilisation, 3), sFat: round(f.stress / 1e6), uFat: round(f.utilisation, 3) }; }),
        'One row per armour layer: wire stress and utilisation at design pressure (with external pressure and tension) and at the factory hydrotest (pressure alone).')],
      notes: [`Burst is governed by the ${burst.governing === 'pressureArmour' ? 'pressure armour' : 'tensile armour'} at ${fmt(burst.burst / BAR)} bar (${fmt(ratio)} × design pressure).`,
        `${gov.name} is the governing check (${fmt(gov.uc)}). Hydrotest pressure ${fmt(pFat / BAR)} bar.`, 'Burst is elastic scaling to tensile strength; a real burst test also mobilises plastic redistribution.', GENERIC],
      advice: [
        adv('engineering', sev(gov.uc), `Design pressure ${fmt(inp.pInt)} bar gives a governing utilisation of ${fmt(gov.uc)} (${gov.name}); calculated burst ${fmt(burst.burst / BAR)} bar.`,
          gov.uc > 1 ? `Lower the design pressure to ${fmt(Math.min(pAllow, loads.pInt / gov.uc) / BAR)} bar or increase the governing armour section by about ${fmt((gov.uc - 1) * 100)} %.` : `The pipe could carry up to ${fmt(pAllow / BAR)} bar on the normal-operation criterion. Confirm with the manufacturer's burst test report and the governing 17J edition.`,
          `${STD_17J}; ${STD_17B}. Factors and the test ratio are recalled defaults — confirm.`),
        adv('economic', 'info', `The pressure armour is ${fmt(paSteel)} t of steel per kilometre; design pressure uses ${fmt(Math.max(...uOp.layers.map((x) => x.utilisation ?? 0)) * 100)} % of the permissible stress.`,
          gov.uc < 0.7 ? `A rating about ${fmt((1 / Math.max(gov.uc, 0.05) - 1) * 100)} % higher is available without change — or roughly ${fmt(paSteel * (1 - gov.uc))} t/km of pressure armour could be saved by re-rating the structure to the actual design pressure.` : 'The structure is matched to the pressure rating; a lower shut-in pressure (HIPPS) is the only lever that would let a lighter pipe be used.',
          'Whole-life costing per ISO 15686-5; quantity-based estimating per AACE RP 18R-97',
          `Pressure armour steel embodies about ${fmt(paSteel * CO2_STEEL)} t CO₂e per kilometre of pipe (screening factor ${CO2_STEEL} t CO₂e/t); total steel in the wall is ${fmt(ms.steel)} t/km.`),
      ],
      publish: { 'flexible.burstPressure': burst.burst },
    };
  },
  verify: [{
    // Hoop-mode burst from ring statics. A thin liner (bore radius 100 mm, 0.5 mm thick) passes the pressure to the
    // pressure armour at r ≈ 100.25 mm: p·r = σ_pa·f·t·sin²α_pa + ΣN_θ,TA. The two tensile layers carry the end cap
    // p·π·r² = Σ σ·f·2πR·t·cos²α and therefore squeeze inwards with ΣN_θ,TA = Σ σ·f·t·sin²α = p·r²·tan²α/(2·R̄),
    // R̄ = 0.1105 m (mean of the layer radii 0.1085 and 0.1125).
    // p_burst = σ_u·f·t·sin²(89.9°)/(r − r²·tan²30°/(2R̄)) = 1000e6·0.9·0.006/(0.10025 − 0.010050·0.33333/0.221)
    //         = 5.4e6/0.085092 = 63.46 MPa = 634.6 bar.
    name: 'Hoop-mode burst of a pressure armour relieved by the tensile armour squeeze',
    source: 'Thin-ring and end-cap statics (Féret & Bournazel 1987); hand calculation',
    inputs: { boreD: 200, Te: 0, depth: 0, pAnn: 0, pInt: 100, layers: [
      { name: 'Liner', type: 'pressureSheath', t: 0.5, alphaDeg: 0, lay: 0, fill: 1, wireW: 0, E: 0.35, yield: 0, uts: 0, lock: 0 },
      { name: 'Pressure armour', type: 'pressureArmour', t: 6, alphaDeg: 89.9, lay: 1, fill: 0.9, wireW: 0, E: 207, yield: 800, uts: 1000, lock: 0 },
      { name: 'Tensile armour 1', type: 'tensileArmour', t: 4, alphaDeg: 30, lay: 1, fill: 0.9, wireW: 0, E: 207, yield: 1100, uts: 1400, lock: 0 },
      { name: 'Tensile armour 2', type: 'tensileArmour', t: 4, alphaDeg: 30, lay: -1, fill: 0.9, wireW: 0, E: 207, yield: 1100, uts: 1400, lock: 0 }] },
    expect: [{ kpi: 'pBurstPA', value: 634.6, tol: 0.01 }],
  }],
};
