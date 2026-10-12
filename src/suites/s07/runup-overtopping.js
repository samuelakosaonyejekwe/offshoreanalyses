// Wave run-up and overtopping of slopes, rubble mounds and vertical walls, with tolerable-discharge checks,
// hinterland flood volume and the sea-level-rise adaptation pathway.
import { iribarren, huntRunup, ROUGHNESS_FACTORS, obliqueFactor, bermFactor, runupEurOtop, overtoppingSlope, overtoppingSteepSlope, overtoppingVertical,
  overtoppingRubbleMound, TOLERABLE_OVERTOPPING, classifyOvertopping, overtoppingVolumes, seaLevelRiseSensitivity } from '../../engine/coastal/runup.js';
import { G, CO2_CONC, CO2_ROCK, EDITION_NOTE, need, run, finite, r, f, linspace, check, eng, eco, slrInputs, slrRows, SLR_BASIS } from './_common.js';

const ARMOUR = Object.keys(ROUGHNESS_FACTORS);

export default {
  id: 's07-runup-overtopping', suite: 7, title: 'Wave run-up and overtopping',
  summary: 'Run-up level, mean overtopping discharge, individual volumes and the crest level needed to keep overtopping tolerable, now and under sea-level rise.',
  covers: [4], fidelity: 'analytical',
  standards: ['EurOtop Manual (2018) ch. 3, 5, 6, 7', 'CIRIA C683 Rock Manual (2007) §5.1.1', 'EN 1991-1-7 / ISO 21650:2007 (actions from waves on coastal structures)'],
  models: ['Iribarren surf-similarity parameter', 'Hunt run-up', 'EurOtop run-up (Eq. 5.1–5.2)', 'EurOtop overtopping: slopes, rubble mounds, steep slopes, vertical walls', 'Weibull overtopping volumes', 'Hinterland storage (continuity)'],
  limits: `Empirical mean-value or design formulae for simple cross-sections with head-on to oblique long-crested waves; composite profiles, wave walls with bullnoses, wind spray and 3-D effects are not modelled, and neither CFD (VOF) nor physical-model precision is implied. ${EDITION_NOTE}`,
  publishes: ['marine.q', 'marine.Ru2', 'marine.Rc', 'marine.crestRequired'],
  searchTerms: 'wave overtopping EurOtop run-up seawall crest level sea level rise',
  inputs: [
    { key: 'type', label: 'Structure type', type: 'select', default: 'rubble', group: 'Structure', role: 'geometry',
      options: [{ value: 'rubble', label: 'Rubble-mound / armoured slope' }, { value: 'smooth', label: 'Smooth or revetted dike slope' }, { value: 'steep', label: 'Steep or battered wall (cot α ≤ 2)' }, { value: 'vertical', label: 'Plain vertical wall' }] },
    { key: 'crest', label: 'Crest level above datum', unit: 'm', type: 'number', default: 8.5, min: 0.5, max: 40, group: 'Structure', role: 'geometry', design: true, help: 'Top of armour or wall — the flood protection elevation.' },
    { key: 'cotA', label: 'Front slope cot α', unit: '-', type: 'number', default: 2, min: 0, max: 10, group: 'Structure', role: 'geometry', design: true, showIf: { key: 'type', in: ['rubble', 'smooth', 'steep'] } },
    { key: 'armour', label: 'Armour / surface roughness class', type: 'select', default: 'rock-2layer-permeable', group: 'Structure', role: 'material', options: ARMOUR.map((k) => ({ value: k, label: `${k} (γf = ${ROUGHNESS_FACTORS[k]})` })),
      help: 'Roughness factors recalled from EurOtop (2018) Table 6.2: confirm against the current edition.' },
    { key: 'gammaF', label: 'Roughness factor γf override (0 = use class)', unit: '-', type: 'number', default: 0, min: 0, max: 1, group: 'Structure', role: 'material', uncertain: { dist: 'normal', cov: 0.08 }, help: 'Calibrate against physical-model tests where available.' },
    { key: 'Gc', label: 'Armour crest width', unit: 'm', type: 'number', default: 4, min: 0, max: 40, group: 'Structure', role: 'geometry', design: true, showIf: { key: 'type', in: ['rubble'] } },
    { key: 'bermB', label: 'Berm width (0 = none)', unit: 'm', type: 'number', default: 0, min: 0, max: 60, group: 'Structure', role: 'geometry', design: true, showIf: { key: 'type', in: ['rubble', 'smooth'] } },
    { key: 'bermDb', label: 'Berm depth below still water (negative above)', unit: 'm', type: 'number', default: 0, min: -10, max: 15, group: 'Structure', role: 'geometry', showIf: { key: 'type', in: ['rubble', 'smooth'] } },
    { key: 'Hm0', label: 'Significant wave height at the toe', unit: 'm', type: 'number', default: 4.2, min: 0.1, max: 15, group: 'Wave boundary', role: 'bc', bind: 'marine.HsToe', uncertain: { dist: 'lognormal', cov: 0.1 },
      help: 'At the design water level before the sea-level-rise allowance below. If the bound value already includes the rise, set the allowance to zero.' },
    { key: 'Tm10', label: 'Spectral period Tm−1,0', unit: 's', type: 'number', default: 10, min: 2, max: 25, group: 'Wave boundary', role: 'bc', uncertain: { dist: 'normal', cov: 0.08 }, help: '≈ Tp/1.1 for a single-peaked spectrum.' },
    { key: 'betaDeg', label: 'Wave obliquity at the structure', unit: 'deg', type: 'number', default: 15, min: 0, max: 80, group: 'Wave boundary', role: 'bc', bind: 'marine.thetaToeDeg' },
    { key: 'durH', label: 'Storm-peak duration', unit: 'h', type: 'number', default: 3, min: 0.5, max: 24, group: 'Wave boundary', role: 'bc' },
    { key: 'swl', label: 'Design still-water level above datum (tide + surge + set-up)', unit: 'm', type: 'number', default: 2.0, min: -5, max: 15, group: 'Water level', role: 'bc' },
    { key: 'hToe', label: 'Toe depth at the design still-water level', unit: 'm', type: 'number', default: 9, min: 0.3, max: 60, group: 'Water level', role: 'bc', bind: 'marine.hToe' },
    { key: 'foreshore', label: 'Foreshore influences the waves (vertical wall)', type: 'bool', default: true, group: 'Water level', role: 'bc', showIf: { key: 'type', in: ['vertical'] } },
    { key: 'depthLimited', label: 'Toe wave is depth-limited (grows with water level)', type: 'bool', default: true, group: 'Water level', role: 'bc', help: 'Scales Hm0 in proportion to toe depth as the sea rises (constant breaker index).' },
    ...slrInputs(),
    { key: 'approach', label: 'Formula set', type: 'select', default: 'design', options: [{ value: 'mean', label: 'Mean-value approach (prediction, probabilistic use)' }, { value: 'design', label: 'Design / assessment approach (+1 σ)' }], group: 'Criteria', role: 'criteria' },
    { key: 'qTol', label: 'Tolerable mean discharge', unit: 'l/s/m', type: 'number', default: 5, min: 0.001, max: 500, group: 'Criteria', role: 'criteria',
      help: 'EurOtop (2018) ch. 3 guidance recalled here: 0.3–10 l/s/m for people depending on Hm0, 5–75 for vehicles, 1–5 for rubble-mound crests and rear slopes. Confirm against the current edition.' },
    { key: 'VTol', label: 'Tolerable maximum individual volume', unit: 'l/m', type: 'number', default: 2000, min: 10, max: 50000, group: 'Criteria', role: 'criteria' },
    { key: 'front', label: 'Length of overtopped frontage', unit: 'm', type: 'number', default: 400, min: 1, max: 50000, group: 'Hinterland', role: 'geometry' },
    { key: 'area', label: 'Hinterland area that ponds', unit: 'm²', type: 'number', default: 60000, min: 10, max: 1e9, group: 'Hinterland', role: 'geometry' },
    { key: 'drain', label: 'Drainage / return-flow capacity', unit: 'm³/s', type: 'number', default: 0.5, min: 0, max: 1000, group: 'Hinterland', role: 'bc', help: 'Downstream boundary: outfalls, pumps and return flow over the crest lee.' },
    { key: 'pondTol', label: 'Tolerable ponding depth', unit: 'm', type: 'number', default: 0.3, min: 0.01, max: 5, group: 'Hinterland', role: 'criteria' },
  ],
  compute(inp) {
    const gF = inp.type === 'vertical' ? 1 : (inp.gammaF > 0 ? inp.gammaF : ROUGHNESS_FACTORS[inp.armour]);
    const design = inp.approach === 'design';
    const vertical = inp.type === 'vertical';
    if (inp.type === 'steep') need(inp.cotA <= 2, 'The steep-slope formula covers cot α from 0 to 2: choose the slope type for flatter fronts.');
    if (inp.type === 'rubble' || inp.type === 'smooth') need(inp.cotA >= 1, 'Sloping-structure formulae need cot α ≥ 1: choose the steep or vertical type.');
    const tanA = vertical ? Infinity : 1 / Math.max(inp.cotA, 1e-6);
    const gBovt = obliqueFactor({ betaDeg: inp.betaDeg, kind: inp.type === 'rubble' ? 'rubble' : 'overtopping' });
    const gBru = obliqueFactor({ betaDeg: inp.betaDeg, kind: inp.type === 'rubble' ? 'rubble' : 'runup' });
    const duration = inp.durH * 3600, Tm = inp.Tm10 / 1.1;

    // state = { Rc, Hm0, h } → run-up and overtopping; used for the design case, the crest search and the scenarios
    const evalState = ({ Rc, Hm0, h }) => {
      const RcE = Math.max(Rc, 0.01);
      let Ru2, xi = NaN, regime, q, gB = 1;
      if (vertical) {
        const o = overtoppingVertical({ Hm0, Rc: RcE, Tm10: inp.Tm10, h: inp.foreshore ? h : undefined, foreshore: inp.foreshore, design });
        q = o.q; regime = o.regime;
        // 2 % level consistent with Pov = exp[−1.21 (Rc/Hm0)²] for non-impulsive walls (EurOtop 2018 ch. 7; Franco et al. 1994)
        Ru2 = Math.sqrt(-Math.log(0.02) / 1.21) * Hm0;
      } else {
        const ru0 = runupEurOtop({ Hm0, Tm10: inp.Tm10, tanAlpha: tanA, gammaF: gF, gammaBeta: gBru, design, surgingRoughness: inp.type === 'rubble' });
        if (inp.bermB > 0 && inp.type !== 'steep') gB = bermFactor({ B: inp.bermB, Lberm: inp.bermB + 2 * Hm0 * inp.cotA, db: inp.bermDb, Hm0, Ru2: ru0.Ru2 });
        const ru = runupEurOtop({ Hm0, Tm10: inp.Tm10, tanAlpha: tanA, gammaB: gB, gammaF: gF, gammaBeta: gBru, design, surgingRoughness: inp.type === 'rubble' });
        Ru2 = ru.Ru2; xi = ru.xi; regime = ru.regime;
        if (inp.type === 'rubble') q = overtoppingRubbleMound({ Hm0, Tm10: inp.Tm10, tanAlpha: tanA, Rc: RcE, gammaF: gF, gammaBeta: gBovt, Gc: inp.Gc, design }).q * (gB < 1 ? overtoppingSlope({ Hm0, Tm10: inp.Tm10, tanAlpha: tanA, Rc: RcE, gammaF: gF, gammaBeta: gBovt, gammaB: gB, design }).q / overtoppingSlope({ Hm0, Tm10: inp.Tm10, tanAlpha: tanA, Rc: RcE, gammaF: gF, gammaBeta: gBovt, design }).q : 1);
        else if (inp.type === 'smooth') q = overtoppingSlope({ Hm0, Tm10: inp.Tm10, tanAlpha: tanA, Rc: RcE, gammaB: gB, gammaF: gF, gammaBeta: gBovt, design }).q;
        else q = overtoppingSteepSlope({ Hm0, cotAlpha: inp.cotA, Rc: RcE, gammaF: gF, gammaBeta: gBovt }).q;
      }
      return { q: finite(q, 'Overtopping discharge'), Ru2, xi, regime, gB };
    };

    const Rc = inp.crest - inp.swl - inp.slr;
    need(Rc > 0, `The crest (${f(inp.crest)} m) is at or below the design water level including sea-level rise (${f(inp.swl + inp.slr)} m): the structure is submerged and overtopping formulae do not apply. Raise the crest or use the flood-wave module.`);
    const h = inp.hToe + inp.slr;
    const Hm0 = inp.depthLimited ? inp.Hm0 * h / inp.hToe : inp.Hm0;
    const s = run(() => evalState({ Rc, Hm0, h }), 'Overtopping');
    const q = s.q, qTol = inp.qTol / 1000;
    const xi = vertical ? 0 : s.xi;
    const hunt = vertical ? null : huntRunup({ H: Hm0, T: inp.Tm10, tanAlpha: tanA });
    const vol = q > 1e-9 ? overtoppingVolumes(vertical ? { q, Hm0, Tm10: inp.Tm10, Tm, duration, Pov: Math.min(1, Math.exp(-1.21 * (Rc / Hm0) ** 2)), b: 0.75 }
      : { q, Hm0, Tm10: inp.Tm10, Tm, duration, Rc, Ru2: s.Ru2 }) : null;
    const Vmax = vol && Number.isFinite(vol.Vmax) ? vol.Vmax : 0;
    const Pov = vol ? vol.Pov : 0;

    // crest freeboard that just meets the tolerable discharge at the design state (bisection: q falls monotonically with Rc)
    let lo = 0.01, hi = 40;
    for (let i = 0; i < 60; i++) { const mid = 0.5 * (lo + hi); if (evalState({ Rc: mid, Hm0, h }).q > qTol) lo = mid; else hi = mid; }
    const RcReq = hi, crestReq = inp.swl + inp.slr + RcReq, raise = crestReq - inp.crest;

    // scenarios: same storm, higher sea; crest raise restores the tolerable discharge
    const scen = slrRows(inp).map((sc) => {
      const hs = inp.hToe + sc.rise, Hs = inp.depthLimited ? inp.Hm0 * hs / inp.hToe : inp.Hm0, rc = inp.crest - inp.swl - sc.rise;
      const qs = evalState({ Rc: rc, Hm0: Hs, h: hs }).q;
      let a = 0.01, b = 40;
      for (let i = 0; i < 50; i++) { const mid = 0.5 * (a + b); if (evalState({ Rc: mid, Hm0: Hs, h: hs }).q > qTol) a = mid; else b = mid; }
      return { name: sc.name, rise: sc.rise, Rc: r(rc, 3), Hm0: r(Hs, 3), q: r(qs * 1000, 3), ratio: r(qs / qTol, 3), crest: r(inp.swl + sc.rise + b, 3), raise: r(Math.max(0, inp.swl + sc.rise + b - inp.crest), 3) };
    });
    const sens = run(() => seaLevelRiseSensitivity({ rises: linspace(0, 2, 21), Rc: inp.crest - inp.swl, Hm0: inp.Hm0, h: inp.hToe, depthLimited: inp.depthLimited, evaluate: (st) => evalState(st).q }), 'Sea-level-rise sensitivity');

    // hinterland: stored volume = (inflow − drainage) × duration (continuity)
    const Qin = q * inp.front, volIn = Qin * duration, volOut = Math.min(volIn, inp.drain * duration), stored = volIn - volOut, pond = stored / inp.area;
    const classes = classifyOvertopping({ q, Vmax: Vmax || undefined });

    const rcs = linspace(Math.max(0.1, 0.3 * Math.min(Rc, RcReq)), 1.6 * Math.max(Rc, RcReq), 40);
    const pass = q <= qTol;
    // quantities per metre of crest raise: a concrete wave wall 0.6 m thick, or the armour/core widening of a mound
    const concPerM = 0.6 * Math.max(raise, 0), rockPerM = inp.type === 'rubble' ? Math.max(raise, 0) * (inp.Gc + Math.max(raise, 0) * inp.cotA + (inp.crest + inp.hToe - inp.swl) * 2 * inp.cotA * 0.5) * 0.62 * 2.65 : 0;

    const advice = [];
    advice.push(pass
      ? eng('ok', `Mean overtopping is ${f(q * 1000)} l/s/m against ${f(inp.qTol)} l/s/m tolerable with ${f(inp.slr)} m of sea-level rise; the crest has ${f(-raise, 2)} m of reserve.`,
        `Keep the crest at ${f(inp.crest)} m; the reserve is used up at about ${f(sens.find((x) => x.q > qTol)?.rise ?? 2, 2)} m of total rise${(sens.find((x) => x.q > qTol)) ? '' : ' or more'} — confirm with physical-model tests for the final section.`,
        'EurOtop (2018) §3.3 tolerable discharges and §5–7 prediction formulae (confirm current edition)',
        'A crest that already tolerates the end-of-life sea level avoids a later raising campaign and its concrete or rock.')
      : eng('fail', `Mean overtopping is ${f(q * 1000)} l/s/m, ${f(q / qTol, 2)} × the tolerable ${f(inp.qTol)} l/s/m (Rc = ${f(Rc)} m, ${s.regime}).`,
        `Raise the crest by ${f(raise, 2)} m to ${f(crestReq)} m, or keep the crest and add a berm, a rougher armour (lower γf) or a recurved wave wall and re-run.`,
        'EurOtop (2018) §3.3 and §5–7 (confirm current edition); CIRIA C683 §5.1.1.3',
        `Raising by ${f(Math.max(raise, 0), 2)} m takes roughly ${f(concPerM, 2)} m³ of concrete per metre of wall (≈ ${f(concPerM * CO2_CONC, 2)} t CO₂e/m, screening factor)${rockPerM ? ` or about ${f(rockPerM, 2)} t of rock per metre for a full mound raise (≈ ${f(rockPerM * CO2_ROCK, 2)} t CO₂e/m)` : ''}; a berm or rougher armour can reach the same discharge with less material.`));
    const last = scen[scen.length - 1];
    advice.push(eng(last.raise > 0 ? 'warn' : 'info', `Adaptation pathway: at +${f(last.rise)} m (${last.name}) the discharge becomes ${f(last.q)} l/s/m and the crest would need to be ${f(last.crest)} m (${f(last.raise, 2)} m above the present crest).`,
      'Build for the end-of-life scenario adopted for the project and make the higher scenarios possible: size the foundation and crest slab for a future wave wall, and reserve land for a berm.', SLR_BASIS,
      'Designing the base for a later wave wall costs a few per cent more concrete now and avoids demolishing and rebuilding the crest when the trigger level is reached.'));
    advice.push(eco(pond > inp.pondTol ? 'warn' : 'info', `Over a ${f(inp.durH)} h storm peak ${f(volIn, 3)} m³ comes over ${f(inp.front)} m of frontage; ${f(stored, 3)} m³ stays behind after drainage, ponding ${f(pond, 2)} m deep.`,
      pond > inp.pondTol ? `Compare the flood damage with either ${f(Math.max(raise, 0), 2)} m of crest raising over ${f(inp.front)} m or ${f(stored / duration, 2)} m³/s of added drainage; choose on whole-life cost.` : 'No flood-damage cost is implied at this discharge; if the crest is being optimised, each 0.1 m lower saves material but multiplies the discharge — check the pathway table before trimming.',
      'Whole-life costing per ISO 15686-5; expected-annual-damage appraisal as in national flood-risk appraisal guidance',
      'Flooded hinterland carries clean-up, disruption and contaminated-water impacts; pumps add an operating energy demand that a passive higher crest does not.'));

    return {
      kpis: [
        { key: 'q', label: 'Mean overtopping discharge', value: r(q * 1000), unit: 'l/s/m', limit: inp.qTol, sense: 'max', explain: 'Average flow over each metre of crest during the storm peak. Compare with the tolerable discharge for the people, vehicles and structures behind the crest; it rises steeply as freeboard is lost.' },
        { key: 'Ru2', label: '2 % wave run-up height', value: r(s.Ru2), unit: 'm', explain: 'Level above still water exceeded by 2 % of incoming waves. A crest below this level is overtopped by more than one wave in fifty.' },
        { key: 'Rc', label: 'Crest freeboard in the design case', value: r(Rc), unit: 'm', explain: 'Crest level minus the design still-water level including the sea-level-rise allowance. Overtopping depends on Rc/Hm0.' },
        { key: 'RcOverH', label: 'Relative freeboard Rc / Hm0', value: r(Rc / Hm0, 3), unit: '-', explain: 'Governing dimensionless parameter of every overtopping formula; below about 1 the structure is heavily overtopped, above about 2 it is nearly dry.' },
        { key: 'xi', label: 'Iribarren number ξm−1,0', value: r(xi, 3), unit: '-', explain: 'Slope over the square root of wave steepness. Below about 2 waves plunge on the slope; above they surge, and run-up no longer grows with ξ. Zero for a vertical wall.' },
        { key: 'RuHunt', label: 'Hunt run-up (breaking waves)', value: r(hunt ? hunt.R : s.Ru2), unit: 'm', explain: 'Classical R = ξ H estimate for regular breaking waves on a smooth slope, given for comparison with the EurOtop 2 % value; not valid for surging waves or walls.' },
        { key: 'Pov', label: 'Share of waves overtopping', value: r(Pov * 100, 3), unit: '%', explain: 'Percentage of incoming waves that pass over the crest. A few per cent means occasional green-water events; tens of per cent means near-continuous flow.' },
        { key: 'Vmax', label: 'Largest individual overtopping volume', value: r(Vmax * 1000), unit: 'l/m', limit: inp.VTol, sense: 'max', explain: 'Expected maximum volume delivered by a single wave during the storm. It governs the hazard to people and vehicles better than the mean discharge.' },
        { key: 'crestReq', label: 'Crest level for the tolerable discharge', value: r(crestReq), unit: 'm', explain: 'Crest elevation at which the mean discharge equals the tolerable value in the design case. Compare with the present crest to read the raise needed or the reserve available.' },
        { key: 'raise', label: 'Crest raise needed (negative = reserve)', value: r(raise, 3), unit: 'm', explain: 'Difference between the required and the present crest level. Positive values are the strengthening measure; negative values show how much sea-level rise the present crest can still absorb.' },
        { key: 'floodVol', label: 'Water stored behind the crest in the storm', value: r(stored), unit: 'm³', explain: 'Overtopped volume over the frontage and storm duration less what the drainage removes. It is the inflow boundary condition of a hinterland flood model.' },
        { key: 'pond', label: 'Ponding depth in the hinterland', value: r(pond, 3), unit: 'm', limit: inp.pondTol, sense: 'max', explain: 'Stored volume spread over the ponding area. A level-pool estimate of flooding exposure; real depths concentrate in low spots.' },
      ],
      checks: [
        check('Mean discharge against the tolerable value', q * 1000, inp.qTol, 'l/s/m', 'EurOtop (2018) §3.3 (confirm current edition)', 'Guards against flow over the crest that endangers people and vehicles or erodes the crest and rear slope. Above the limit the crest must be closed in storms or raised.'),
        check('Maximum individual volume', Vmax * 1000, inp.VTol, 'l/m', 'EurOtop (2018) §3.3, Table 3.3 (confirm current edition)', 'Guards against a single large wave knocking down a person or pushing a vehicle. It can govern even when the mean discharge is small.'),
        check('Hinterland ponding', pond, inp.pondTol, 'm', 'Level-pool continuity; project flood-protection criterion', 'Guards against flooding of the protected area when overtopping exceeds the drainage capacity over the storm.'),
      ],
      plots: [
        { type: 'line', title: 'Overtopping discharge versus crest freeboard', x: { label: 'Freeboard Rc', unit: 'm' }, y: { label: 'q', unit: 'l/s/m' },
          series: [{ name: 'q', x: rcs.map((v) => r(v, 3)), y: rcs.map((v) => r(evalState({ Rc: v, Hm0, h }).q * 1000, 3)) }, { name: 'Tolerable', x: [rcs[0], rcs[rcs.length - 1]], y: [inp.qTol, inp.qTol], dash: true }, { name: 'Present crest', x: [r(Rc, 3)], y: [r(q * 1000, 3)] }],
          explain: 'Discharge in the design case as the crest is raised or lowered. The crossing with the dashed line is the required freeboard; the steepness shows how much each 0.1 m of crest is worth.' },
        { type: 'line', title: 'Overtopping versus sea-level rise', x: { label: 'Sea-level rise', unit: 'm' }, y: { label: 'q', unit: 'l/s/m' },
          series: [{ name: 'Present crest', x: sens.map((x) => r(x.rise, 3)), y: sens.map((x) => r(x.q * 1000, 3)) }, { name: 'Tolerable', x: [0, 2], y: [inp.qTol, inp.qTol], dash: true }],
          explain: 'Discharge over the unchanged crest as the sea rises. The rise at which the curve crosses the tolerable line is the adaptation trigger for raising the crest.' },
        { type: 'line', title: 'Crest raise that keeps present performance', x: { label: 'Sea-level rise', unit: 'm' }, y: { label: 'Crest raise', unit: 'm' },
          series: [{ name: 'Raise restoring today\'s discharge', x: sens.map((x) => r(x.rise, 3)), y: sens.map((x) => r(Number.isFinite(x.crestRaise) ? x.crestRaise : 10, 3)) }, { name: '1 : 1', x: [0, 2], y: [0, 2], dash: true }],
          explain: 'Crest raise needed to hold the present-day discharge. Where the curve runs above the 1:1 line the wave is depth-limited and the crest must rise faster than the sea.' },
      ].concat(vol ? [{ type: 'line', title: 'Individual overtopping volumes', x: { label: 'Volume per wave', unit: 'l/m' }, y: { label: 'Probability of exceedance per overtopping wave', unit: '-' },
        series: [(() => { const vs = linspace(0, Math.max(Vmax, vol.a * 3), 40); return { name: 'Weibull', x: vs.map((v) => r(v * 1000, 3)), y: vs.map((v) => r(Math.exp(-Math.pow(v / vol.a, vol.b)), 3)) }; })()],
        explain: 'Distribution of the volumes thrown over by individual overtopping waves. The long tail is why the largest wave carries many times the mean volume.' }] : []),
      tables: [
        { title: 'Sea-level-rise adaptation pathway', columns: [{ key: 'name', label: 'Scenario' }, { key: 'rise', label: 'Rise', unit: 'm' }, { key: 'Rc', label: 'Freeboard', unit: 'm' }, { key: 'Hm0', label: 'Hm0', unit: 'm' }, { key: 'q', label: 'q', unit: 'l/s/m' }, { key: 'ratio', label: 'q / tolerable', unit: '-' }, { key: 'crest', label: 'Required crest', unit: 'm' }, { key: 'raise', label: 'Raise', unit: 'm' }],
          rows: scen, explain: 'The design storm at each sea-level-rise scenario with the present crest, and the crest level each scenario would need. Use it to fix trigger levels for raising.' },
        { title: 'Tolerable overtopping classification', columns: [{ key: 'd', label: 'Hazard' }, { key: 'ql', label: 'q limit', unit: 'l/s/m' }, { key: 'u', label: 'q / limit', unit: '-' }, { key: 'ok', label: 'Result' }],
          rows: classes.map((c) => ({ d: c.description, ql: r(c.qLimit * 1000, 3), u: r(c.utilisationQ, 3), ok: c.ok ? 'tolerable' : 'exceeded' })),
          explain: 'The computed discharge and maximum volume against each tolerable limit for people, vehicles, structures and property. Limits are recalled from EurOtop (2018) ch. 3 and must be confirmed against the current edition.' },
      ],
      notes: [
        `${inp.type} structure, ${design ? 'design/assessment' : 'mean-value'} formulae, γf = ${f(gF)}, γβ = ${f(gBovt)}, γb = ${f(s.gB)}; regime: ${s.regime}.`,
        `Design case: Hm0 = ${f(Hm0)} m and h = ${f(h)} m after ${f(inp.slr)} m of sea-level rise${inp.depthLimited ? ' (depth-limited scaling)' : ''}.`,
      ],
      advice,
      publish: { 'marine.q': q, 'marine.Ru2': s.Ru2, 'marine.Rc': Rc, 'marine.crestRequired': crestReq },
      quantities: { concrete_m3: r(inp.type === 'rubble' ? 0 : concPerM * inp.front, 3), rock_t: r(rockPerM * inp.front, 3) },
    };
  },
  verify: [
    { name: 'Plain vertical wall without foreshore, mean value: q = 0.047 √(g Hm0³) exp[−(2.35 Rc/Hm0)^1.3]',
      source: 'EurOtop (2018) Eq. 7.1. Hand calculation, Hm0 = 2 m, Rc = 2 m: (2.35)^1.3 = 3.0366, exp(−3.0366) = 0.04800, √(9.80665 × 8) = 8.8574 → q = 0.047 × 0.04800 × 8.8574 = 0.01998 m³/s/m',
      inputs: { type: 'vertical', foreshore: false, approach: 'mean', Hm0: 2, crest: 4, swl: 2, slr: 0, betaDeg: 0, depthLimited: false },
      expect: [{ kpi: 'q', value: 19.98, tol: 0.005 }] },
    { name: 'Smooth 1:4 slope, mean run-up Ru2% = 1.65 γ ξ Hm0',
      source: 'EurOtop (2018) Eq. 5.1. Hand calculation, Hm0 = 2 m, Tm−1,0 = 6 s: L0 = 56.21 m, s = 0.03558, ξ = 0.25/√0.03558 = 1.3253 → Ru2 = 1.65 × 1.3253 × 2 = 4.374 m (below the 5.39 m maximum); Hunt R = ξ H = 2.651 m',
      inputs: { type: 'smooth', armour: 'smooth', cotA: 4, approach: 'mean', Hm0: 2, Tm10: 6, betaDeg: 0, crest: 9, swl: 2, slr: 0, depthLimited: false },
      expect: [{ kpi: 'Ru2', value: 4.374, tol: 0.005 }, { kpi: 'xi', value: 1.3253, tol: 0.005 }, { kpi: 'RuHunt', value: 2.651, tol: 0.005 }] },
  ],
};
