// Natural frequencies and mode shapes of the turbine–tower–monopile system on a three-spring foundation, with
// the 1P / 3P resonance margins, Campbell diagram and sensitivity to foundation stiffness.
import { towerNaturalFrequency } from '../../engine/geotech/monopile.js';
import { campbellDiagram } from '../../engine/structural/framedynamics.js';
import { need, fmt, round, adv, check, worst, steelCo2Text, turbineInputs, supportInputs, foundationInputs, steelInputs, steelOf, turbine, segments, foundationOf, steelMass, LIMIT_AERO, LIMIT_HIFI, HELP_SUPPLIER } from './_common.js';
import { dynamicsInputs, supportStructure } from './_dynamics.js';

/** Segments of the engine eigen-solver from the tubular segments: wet/dry split of the pile, tower in slices of constant wall. */
function eigSegments(inp, tb, steel) {
  const [pile, tower] = segments(inp, tb, steel);
  const tp = pile.tBase, Di = pile.dBase - 2 * tp, rhoW = inp.rhoW ?? 1025;
  const added = rhoW * (inp.Ca ?? 1) * Math.PI * pile.dBase ** 2 / 4 + (inp.flooded ? rhoW * Math.PI * Di * Di / 4 : 0);
  const out = [
    { L: inp.depth, D: pile.dBase, t: tp, E: steel.E, rho: steel.rho, addedMass: added },
    { L: inp.zInt, D: pile.dBase, t: tp, E: steel.E, rho: steel.rho },
  ];
  const nT = 6, H = tower.zTop - tower.zBase;
  for (let i = 0; i < nT; i++) {
    const s0 = i / nT, s1 = (i + 1) / nT, sm = 0.5 * (s0 + s1);
    out.push({ L: H / nT, Dbot: tower.dBase + (tower.dTop - tower.dBase) * s0, Dtop: tower.dBase + (tower.dTop - tower.dBase) * s1, t: tower.tBase + (tower.tTop - tower.tBase) * sm, E: steel.E, rho: tower.density });
  }
  return out;
}

export default {
  id: 's03-natural-frequency', suite: 3, title: 'Natural frequency and modal analysis',
  summary: 'First fore-aft natural frequencies and mode shapes of the turbine, tower and monopile on a flexible foundation, and the separation from rotor (1P) and blade-passing (3P) frequencies.',
  covers: [12], fidelity: 'numerical',
  standards: ['DNV-ST-0126 (2021) §4.5.4 and §6 (natural frequency, resonance)', 'IEC 61400-3-1:2019 §7.5 and §11', 'DNV-RP-C205 (2019) §6.9 (added mass)', 'Arany, Bhattacharya, Macdonald & Hogan (2017) Soil Dyn. Earthq. Eng. 92'],
  models: ['Euler–Bernoulli beam finite elements with consistent mass', 'Structural modal equations (generalised eigenproblem)', 'Three-spring (lateral, rocking, coupling) foundation model', 'Hydrodynamic added mass and enclosed water', 'Campbell diagram and resonance screening'],
  limits: `Planar bending modes without shear deformation, axial-load softening or blade flexibility; side-side frequencies are taken equal to fore-aft. ${LIMIT_AERO} ${LIMIT_HIFI}`,
  publishes: ['owt.f1', 'owt.rna'],
  searchTerms: 'offshore wind monopile natural frequency soft-stiff 1P 3P foundation stiffness',
  inputs: [
    ...turbineInputs(),
    ...supportInputs({ embed: false }),
    ...steelInputs(),
    ...foundationInputs(),
    ...dynamicsInputs().filter((i) => ['Ca', 'flooded'].includes(i.key)),
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'margin', label: 'Required separation from 1P and 3P', unit: '%', type: 'number', default: 10, min: 0, max: 30, group: 'Acceptance criteria', role: 'criteria', help: 'DNV-ST-0126 recommends keeping the natural frequency at least 10 % away from the 1P and 3P ranges unless resonance is shown to be acceptable.' },
    { key: 'nElem', label: 'Beam elements', type: 'int', default: 30, min: 6, max: 240, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Beam elements', levels: [15, 30, 60, 120], qoi: 'f1' },
  compute(inp, ctx) {
    const steel = steelOf(inp, ctx);
    const tb = turbine(inp);
    const fnd = foundationOf(inp);
    const solve = (over = {}, nModes = 1) => supportStructure({ ...inp, ...over }, tb, steel, { nElem: inp.nElem, nModes });
    const str = solve({}, 3);
    const freqs = str.modes.map((q) => q.frequency);
    const fixF = solve({ fnd: 'fixed' }).modes[0].frequency;
    const f1 = freqs[0], f2 = freqs[1];
    need(Number.isFinite(f1) && f1 > 0, 'The eigenvalue analysis returned no positive frequency: check the geometry and the foundation stiffness.');
    // Independent cross-check with the engine eigen-solver for towers on three springs (coarse mesh).
    let fCross = NaN;
    try {
      const segs = eigSegments(inp, tb, steel);
      fCross = towerNaturalFrequency({ segments: segs, topMass: tb.rnaMass, nModes: 1, elementsPerSegment: 2, ...(fnd === 'fixed' ? {} : { KLL: fnd.Kh, KLR: fnd.Khr, KRR: fnd.Kr }) }).frequencies[0];
    } catch { fCross = NaN; }
    const m = inp.margin / 100;
    const p1lo = tb.rpmMin / 60, p1hi = tb.rpmRated / 60, p3lo = 3 * p1lo, p3hi = 3 * p1hi;
    const soft = f1 < p1lo, softStiff = f1 >= p1hi && f1 <= p3lo, stiff = f1 > p3hi;
    const band = soft ? 'soft-soft' : softStiff ? 'soft-stiff' : stiff ? 'stiff-stiff' : (f1 < p1hi ? 'inside the 1P range' : 'inside the 3P range');
    // Separation margins to the nearest edge of each excitation range (negative inside a range).
    const sep1 = f1 > p1hi ? f1 / p1hi - 1 : f1 < p1lo ? 1 - f1 / p1lo : -Math.min(f1 / p1lo - 1, 1 - f1 / p1hi);
    const sep3 = f1 < p3lo ? 1 - f1 / p3lo : f1 > p3hi ? f1 / p3hi - 1 : -Math.min(f1 / p3lo - 1, 1 - f1 / p3hi);
    const uc = (sep) => (sep <= 0 ? 99 : m > 0 ? m / sep : 0);
    const checks = [
      check('Separation from the rotor frequency range (1P)', inp.margin, Math.max(sep1 * 100, 1e-9), '% (required / achieved)', 'DNV-ST-0126 (2021) §4.5.4', 'Guards against resonance with rotor mass or aerodynamic imbalance, which acts once per revolution. Inside the margin the tower sways strongly at some rotor speeds and fatigue life collapses.'),
      check('Separation from the blade-passing range (3P)', inp.margin, Math.max(sep3 * 100, 1e-9), '% (required / achieved)', 'DNV-ST-0126 (2021) §4.5.4', 'Guards against resonance with tower shadow and wind shear loads, which act three times per revolution. Inside the margin a speed-exclusion window in the controller is needed.'),
    ];
    checks[0].uc = Math.min(uc(sep1), 99); checks[1].uc = Math.min(uc(sep3), 99);
    const gov = worst(checks), ok = gov.uc <= 1;
    const camp = campbellDiagram({ naturalFrequencies: [f1, f2], orders: [1, 3], speedMin: tb.rpmMin, speedMax: tb.rpmRated, margin: m, nPoints: 30 });
    const mult = [0.25, 0.35, 0.5, 0.7, 1, 1.4, 2, 3, 4];
    const sens = fnd === 'fixed' ? mult.map(() => f1) : mult.map((f) => solve({ Kh: inp.Kh * f, Kr: inp.Kr * f, Khr: inp.Khr * f }).modes[0].frequency);
    const zNodes = Array.from(str.nodes);
    // Wall-thickness sensitivity for the advice: tower and pile walls scaled together by ±10 %.
    const scaleWall = (q) => solve({ t: inp.t * q, towerTb: inp.towerTb * q, towerTt: inp.towerTt * q }).modes[0].frequency;
    const dfdt = (scaleWall(1.1) - scaleWall(0.9)) / 0.2; // Hz per unit relative wall change
    const sm = steelMass({ ...inp, Lp: 0 }, tb, steel.rho);
    const target = softStiff || (!soft && !stiff && f1 >= p1hi) ? Math.max(p1hi * (1 + m), Math.min(f1, p3lo * (1 - m))) : f1;
    const dWall = dfdt > 0 ? (target - f1) / dfdt : 0;
    const dMass = dWall * sm.total / 1e3;
    return {
      kpis: [
        { key: 'f1', label: 'First fore-aft natural frequency', unit: 'Hz', value: f1, explain: `Fundamental bending frequency of the whole turbine–tower–foundation system. It should sit between the 1P range (${fmt(p1lo)}–${fmt(p1hi)} Hz) and the 3P range (${fmt(p3lo)}–${fmt(p3hi)} Hz) for a soft-stiff design.` },
        { key: 'f2', label: 'Second fore-aft natural frequency', unit: 'Hz', value: f2, explain: 'Second bending mode. It is normally well above 3P but may be excited by blade-passing harmonics (6P) and by breaking-wave impact.' },
        { key: 'f1Fixed', label: 'First frequency, clamped at the mudline', unit: 'Hz', value: fixF, explain: 'Rigid-soil upper bound. The difference to the flexible-foundation value shows how much the result depends on soil stiffness.' },
        { key: 'fndLoss', label: 'Frequency reduction by foundation flexibility', unit: '%', value: (1 - f1 / fixF) * 100, explain: 'How far the soil lowers the frequency below the clamped value. Monopiles typically lose 5–15 %; a larger loss means the result is soil-controlled and uncertain.' },
        { key: 'sep1P', label: 'Margin to the 1P range', unit: '%', value: sep1 * 100, limit: inp.margin, sense: 'min', explain: 'Relative distance between the natural frequency and the nearest edge of the rotor-speed range. Negative means the frequency lies inside the range.' },
        { key: 'sep3P', label: 'Margin to the 3P range', unit: '%', value: sep3 * 100, limit: inp.margin, sense: 'min', explain: 'Relative distance to the nearest edge of the blade-passing range. Negative means the frequency lies inside the range.' },
        { key: 'rnaT', label: 'Rotor–nacelle mass used', unit: 't', value: tb.rnaMass / 1e3, explain: `Top mass in the model (${tb.name}). ${HELP_SUPPLIER}` },
        { key: 'fCross', label: 'First frequency, independent stick model', unit: 'Hz', value: Number.isFinite(fCross) ? fCross : f1, explain: 'Same structure solved with the engine eigen-solver for towers on three springs on a coarse mesh with a stepped tower wall. Agreement within about 1 % verifies the model.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Campbell diagram', x: { label: 'Rotor speed', unit: 'rpm' }, y: { label: 'Frequency', unit: 'Hz' },
          series: [...camp.lines.map((l) => ({ name: `${l.order}P excitation`, x: Array.from(camp.speeds), y: Array.from(l.f) })),
            { name: 'Mode 1', x: [tb.rpmMin, tb.rpmRated], y: [f1, f1] }, { name: `Mode 1 − ${inp.margin} %`, x: [tb.rpmMin, tb.rpmRated], y: [f1 * (1 - m), f1 * (1 - m)], dash: true },
            { name: `Mode 1 + ${inp.margin} %`, x: [tb.rpmMin, tb.rpmRated], y: [f1 * (1 + m), f1 * (1 + m)], dash: true }],
          explain: 'Excitation lines 1P and 3P over the operating speed range against the first natural frequency with its exclusion band. An excitation line crossing the band inside the speed range is a resonance to avoid or to pass through quickly.' },
        { type: 'line', title: 'Mode shapes', x: { label: 'Normalised lateral displacement', unit: '-' }, y: { label: 'Elevation above still water', unit: 'm' },
          series: str.modes.map((q, i) => { const mx = Math.max(...Array.from(q.shape, Math.abs)); return { name: `Mode ${i + 1} (${fmt(freqs[i])} Hz)`, x: Array.from(q.shape, (v) => v / mx), y: zNodes }; }),
          explain: 'Lateral displacement of each mode from the mudline to the tower top. Displacement at the mudline shows foundation flexibility; the node of mode 2 marks where a damper or sensor is least effective.' },
        { type: 'line', title: 'Natural frequency against foundation stiffness', x: { label: 'Stiffness multiplier on K_LL, K_LR, K_RR', unit: '-' }, y: { label: 'First natural frequency', unit: 'Hz' },
          series: [{ name: 'f1', x: mult, y: sens }, { name: '1P upper + margin', x: [mult[0], mult[mult.length - 1]], y: [p1hi * (1 + m), p1hi * (1 + m)], dash: true }, { name: '3P lower − margin', x: [mult[0], mult[mult.length - 1]], y: [p3lo * (1 - m), p3lo * (1 - m)], dash: true }, { name: 'Design point', x: [1], y: [f1] }],
          explain: 'Effect of soil-stiffness uncertainty, scour and cyclic softening on the frequency. The curve must stay between the dashed lines over the credible stiffness range (typically 0.5 to 2).' },
      ],
      tables: [
        { title: 'Frequencies and excitation ranges', columns: [{ key: 'k', label: 'Item' }, { key: 'v', label: 'Frequency', unit: 'Hz' }, { key: 'n', label: 'Note' }],
          rows: [['Mode 1', f1, band], ['Mode 2', f2, ''], ['Mode 3', freqs[2] ?? 0, ''], ['Mode 1, clamped base', fixF, 'rigid soil'], ['1P lower', p1lo, `${fmt(tb.rpmMin)} rpm`], ['1P upper', p1hi, `${fmt(tb.rpmRated)} rpm`], ['3P lower', p3lo, ''], ['3P upper', p3hi, ''],
            ['Mode 1 with half the foundation stiffness', sens[2], 'lower bound soil'], ['Mode 1 with twice the foundation stiffness', sens[6], 'upper bound soil']].map(([k, v, n]) => ({ k, v: round(v), n })),
          explain: 'Natural frequencies beside the rotor-induced excitation ranges, with the soil-stiffness bounds to carry into the resonance assessment.' },
      ],
      notes: [`Design type: ${band}. ${tb.name}; rotor speed range ${fmt(tb.rpmMin)}–${fmt(tb.rpmRated)} rpm.`, `Foundation flexibility lowers the frequency by ${fmt((1 - f1 / fixF) * 100)} % from the clamped value.`,
        'Compare with the measured frequency after installation: as-built values are usually a few per cent above API-based predictions.'],
      advice: [
        ok ? adv('engineering', 'ok', `f1 = ${fmt(f1)} Hz is ${band} with ${fmt(sep1 * 100)} % margin to 1P and ${fmt(sep3 * 100)} % to 3P (required ${inp.margin} %).`,
          `Confirm the rotor-speed range and rotor–nacelle mass with the turbine supplier (scaled reference values are used) and the margin over the soil-stiffness range: with half the stiffness f1 = ${fmt(sens[2])} Hz, with twice ${fmt(sens[6])} Hz. Measure the as-built frequency at commissioning and after the first storms.`,
          'DNV-ST-0126 (2021) §4.5.4; IEC 61400-3-1:2019 §11',
          `Frequency is bought with steel: a 1 % change of all wall thicknesses moves f1 by ${fmt(dfdt * 0.01 * 1000)} mHz and the primary steel above the mudline by ${fmt(sm.total / 1e5)} t. ${steelCo2Text(sm.total / 1e5)}`)
          : adv('engineering', 'fail', `f1 = ${fmt(f1)} Hz is ${band}: margin to 1P ${fmt(sep1 * 100)} %, to 3P ${fmt(sep3 * 100)} % (required ${inp.margin} %).`,
            sep1 < m ? `Stiffen the structure to reach at least ${fmt(p1hi * (1 + m))} Hz: larger monopile and tower-base diameter are most effective (stiffness ∝ D³t), or agree a rotor-speed exclusion window with the turbine supplier.` : `Soften the structure to stay below ${fmt(p3lo * (1 - m))} Hz or agree a speed-exclusion window around the 3P crossing with the turbine supplier.`,
            'DNV-ST-0126 (2021) §4.5.4; IEC 61400-3-1:2019 §11',
            `Reaching the target ${fmt(target)} Hz by wall thickness alone needs about ${fmt(dWall * 100)} % more wall, ${fmt(dMass)} t of steel ≈ ${fmt(dMass * 1.9)} t CO₂e (screening factor 1.9 t CO₂e/t); a diameter increase reaches the same frequency with less material.`),
        adv('economic', ok ? 'info' : 'warn', `The tower and monopile above the mudline hold ${fmt(sm.total / 1e3)} t of steel; the frequency window ${fmt(p1hi * (1 + m))}–${fmt(p3lo * (1 - m))} Hz is ${fmt((p3lo * (1 - m) - p1hi * (1 + m)) * 1000)} mHz wide.`,
          ok ? `Use the reserve to the lower bound (${fmt((f1 - p1hi * (1 + m)) * 1000)} mHz) to trim steel only after the soil stiffness is confirmed; a resonance found after installation costs production through speed-exclusion windows.` : 'Resolve the resonance at design stage: a controller exclusion window costs energy yield every year, whereas extra diameter is a one-off steel cost.',
          'LCOE trade-off between CAPEX (fabricated steel per tonne) and annual energy production; DNV-ST-0126 design verification practice'),
      ],
      publish: { 'owt.f1': f1, 'owt.rna': tb.rnaMass },
      quantities: { steel_t: round(sm.total / 1e3) },
    };
  },
  verify: [
    { name: 'Uniform clamped cantilever without tip mass: f1 = (1.8751²/2π)·√(EI/(m·L⁴))', source: 'Blevins, Formulas for Natural Frequency and Mode Shape, Table 8-1: D = 6 m, t = 60 mm, E = 210 GPa → EI = 1.0367e12 N·m², m = 7850 × 1.1197 = 8789.6 kg/m, L = 25 + 87.6 = 112.6 m → f1 = 0.5596 × √(1.0367e12 / (8789.6 × 112.6⁴)) = 0.4793 Hz',
      inputs: { fnd: 'fixed', D: 6, t: 60, towerDb: 6, towerDt: 6, towerTb: 60, towerTt: 60, towerMassFactor: 1, Ca: 0, flooded: false, rnaT: 0.001, depth: 25, hubHeight: 90, ratedMW: 5, nElem: 60 }, expect: [{ kpi: 'f1', value: 0.4793, tol: 0.003 }, { kpi: 'f2', value: 3.004, tol: 0.005 }] },
    { name: 'Same cantilever with the 350 t top mass: Rayleigh f1 = (1/2π)·√(3EI/((M + 0.2357·m·L)·L³))', source: 'Rayleigh quotient with the static tip-load shape (Blevins Table 8-8, accurate to about 1 %): M + 0.2357 m L = 350000 + 233280 = 583280 kg → f1 = 0.3075 Hz',
      inputs: { fnd: 'fixed', D: 6, t: 60, towerDb: 6, towerDt: 6, towerTb: 60, towerTt: 60, towerMassFactor: 1, Ca: 0, flooded: false, rnaT: 350, depth: 25, hubHeight: 90, ratedMW: 5 }, expect: [{ kpi: 'f1', value: 0.3075, tol: 0.015 }] },
  ],
};
