// Wax deposition with thermal feedback, pigging interval and gel restart pressure of a waxy oil line, and mineral
// scale tendency of the produced water (calcite saturation, sulphate scale on mixing with sea water, deposition rate).
import { makeSolubilityCurve, waxDepositionLine, piggingInterval, gelRestartPressure, gellingScreen, diffusivityHaydukMinhas, waxAppearanceIdeal } from '../../engine/flow/wax.js';
import { SEAWATER_MG_PER_KG, toMolal, calciteSaturation, seawaterMixingScale, scaleDepositionRate } from '../../engine/flow/scale.js';
import { overallU } from '../../engine/flow/thermalpipe.js';
import { pipeInputs, pipeFrom, lengthInput, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, sev, worst, arr,
  vesselCo2, LIMIT_NO_CFD, ECON_BASIS, BAR, KELVIN, YEAR } from './_common.js';
import { thermalInputs, thermalFrom } from './_flow.js';

const G_O = 'Oil and operation';
const G_W = 'Wax';
const G_S = 'Produced water and scale';
const G_C = 'Acceptance criteria';
const WATER = [{ ion: 'Na', mgL: 30000 }, { ion: 'K', mgL: 400 }, { ion: 'Ca', mgL: 2500 }, { ion: 'Mg', mgL: 400 }, { ion: 'Ba', mgL: 250 }, { ion: 'Sr', mgL: 500 }, { ion: 'Cl', mgL: 52500 }, { ion: 'SO4', mgL: 10 }, { ion: 'HCO3', mgL: 500 }];
const PARAFFIN = [{ M: 170, z: 88 }, { M: 282, z: 6 }, { M: 352, z: 3.5 }, { M: 422, z: 1.8 }, { M: 492, z: 0.7 }];
const DAY = 86400;

export default {
  id: 's08-wax-scale',
  suite: 8,
  title: 'Wax deposition and scale',
  summary: 'Predicts wax build-up along a waxy oil line, when it must be pigged and what it takes to restart a gelled line, and screens the produced water for carbonate and sulphate scale.',
  covers: [15],
  fidelity: 'numerical',
  standards: ['Burger, Perkins & Striegler (1981) molecular-diffusion wax deposition', 'Singh et al. (2000) thermal feedback and ageing', 'Matzain et al. (2002) deposition model', 'Plummer & Busenberg (1982) carbonate equilibria', 'NACE TM0374 / Oddo & Tomson scale-tendency practice'],
  models: ['Wax appearance temperature and solubility curve', 'Ideal-solution solid–liquid equilibrium (Won) as a cross-check of the WAT', 'Molecular diffusion (Fick) wax deposition with thermal feedback', 'Matzain semi-empirical deposition with shear stripping', 'Gel restart pressure from wall yield stress', 'Aqueous speciation with Davies activity coefficients and saturation indices', 'Seawater–formation-water mixing line', 'Surface-reaction scale growth'],
  limits: `Single-phase oil with quasi-steady hydraulics and heat transfer; the wax solubility curve, deposit wax fraction and conductivity are user data (the curve should come from DSC or filtration, not from the ideal-solution model). Scale is an equilibrium screening with a Davies activity model (no ion pairing, no Pitzer terms), independent minerals and a user rate constant. Multiphase wax deposition, deposit ageing kinetics and CFD–DEM deposition are outside this module. ${LIMIT_NO_CFD}`,
  publishes: ['pipeline.piggingInterval', 'pipeline.waxThickness'],
  searchTerms: 'wax deposition pipeline pigging frequency gel restart scale saturation index barite calcite',
  inputs: [
    ...pipeInputs({ omit: ['rhoCont', 'wsGiven', 'rhoW', 'steel'], defaults: { tCoat: 50, rhoCoat: 750, tConc: 0 } }),
    lengthInput(20, 'Pipe'),
    num('rough', 'Wall roughness', 'mm', 0.045, 0.001, 3, 'Pipe', 'geometry', 'Clean-pipe roughness.'),
    num('qOil', 'Oil rate', 'Sm³/d', 4000, 10, 300000, G_O, 'bc', 'Stock-tank oil rate through the line.', { design: true, uncertain: { dist: 'normal', cov: 0.1 } }),
    num('rhoOil', 'Oil density', 'kg/m³', 850, 600, 1050, G_O, 'material', 'Density at line conditions.'),
    num('mu40', 'Oil viscosity at 40 °C', 'cP', 6, 0.2, 5000, G_O, 'material', 'Dynamic viscosity at 40 °C.', { uncertain: { dist: 'lognormal', cov: 0.2 } }),
    num('muB', 'Viscosity–temperature constant', 'K', 2500, 0, 8000, G_O, 'material', 'μ(T) = μ40·exp(B·(1/T − 1/313.15)); about 1500–4000 K for crude oils.'),
    num('cpOil', 'Oil specific heat', 'J/kgK', 2000, 1200, 3000, G_O, 'material', 'Heat capacity of the oil.'),
    num('kOil', 'Oil thermal conductivity', 'W/mK', 0.13, 0.08, 0.2, G_O, 'material', 'Conductivity of the oil.'),
    ...thermalInputs({ Tin: 60 }).map((i) => (i.key === 'kCoat' ? { ...i, default: 0.17 } : i)),
    num('WAT', 'Wax appearance temperature', '°C', 32, -10, 80, G_W, 'material', 'Cloud point from DSC or cross-polar microscopy at line pressure.', { uncertain: { dist: 'normal', cov: 0.05 } }),
    num('waxC', 'Dissolved wax at the WAT', 'wt %', 5, 0.1, 40, G_W, 'material', 'Wax content in solution at the wax appearance temperature.'),
    num('waxB', 'Slope of the solubility curve', '1/K', 0.05, 0.005, 0.3, G_W, 'material', 'C = C_WAT·exp(b·(T − WAT)) below the WAT, from the precipitation curve.'),
    { key: 'solCurve', label: 'Measured solubility curve (optional)', type: 'series', default: { x: [], y: [] }, xLabel: 'Temperature', xUnit: '°C', yLabel: 'Dissolved wax', yUnit: 'wt %', group: G_W, role: 'material',
      help: 'At least three points with increasing temperature replace the exponential curve.' },
    { key: 'paraffin', label: 'n-Paraffin distribution for the ideal-solution cross-check (optional)', type: 'table', group: G_W, role: 'material', columns: [{ key: 'M', label: 'Molar mass', unit: 'g/mol', type: 'number' }, { key: 'z', label: 'Amount', unit: 'mol %', type: 'number' }], default: PARAFFIN,
      help: 'Pseudo-components of the oil (solvent first, then the wax-forming n-paraffins). Used only to report the ideal solid–liquid-equilibrium wax appearance temperature, which over-predicts the measured value by 5–20 K.' },
    sel('waxModel', 'Deposition model', 'fick', [['fick', 'Molecular diffusion (Fick)'], ['matzain', 'Matzain semi-empirical']], G_W, 'control', 'Fick with a deposit wax fraction, or Matzain with enhancement and shear stripping.'),
    num('waxFrac', 'Wax fraction of the deposit', '-', 0.3, 0.05, 1, G_W, 'material', 'Fresh deposits hold 10–40 % wax with trapped oil; it rises with age.', { showIf: { key: 'waxModel', in: ['fick'] } }),
    num('kDep', 'Deposit thermal conductivity', 'W/mK', 0.25, 0.1, 0.5, G_W, 'material', 'The deposit insulates the wall and slows further deposition.'),
    num('thick0', 'Deposit thickness at the start', 'mm', 0, 0, 20, G_W, 'ic', 'Uniform deposit left after the last pig run (the model starts from a clean or uniformly coated wall).'),
    num('days', 'Operating time simulated', 'days', 60, 2, 720, G_W, 'control', 'Time since the last pig run.'),
    num('pour', 'Pour point', '°C', -5, -30, 60, G_W, 'material', 'ASTM D97 / D5853 pour point for the gelling screen.'),
    num('tauY', 'Gel yield stress at seabed temperature', 'Pa', 50, 0, 2000, G_W, 'material', 'Yield stress of the gelled oil after the planned shut-in time.'),
    { key: 'water', label: 'Formation water analysis', type: 'table', group: G_S, role: 'material', columns: [{ key: 'ion', label: 'Ion (Na, K, Ca, Mg, Ba, Sr, Fe, Cl, SO4, HCO3)', type: 'text' }, { key: 'mgL', label: 'Concentration', unit: 'mg/L', type: 'number' }], default: WATER,
      help: 'Produced-water ion analysis. Barium and strontium with sulphate from injected sea water form sulphate scale; calcium and bicarbonate form carbonate scale.' },
    num('Tsc', 'Temperature for the scale assessment', '°C', 80, 5, 200, G_S, 'bc', 'Water temperature at the location assessed.'),
    num('Psc', 'Pressure for the scale assessment', 'bar', 100, 1, 1000, G_S, 'bc', 'Pressure at the location assessed.'),
    num('pCO2', 'CO₂ partial pressure', 'bar', 1, 0.001, 100, G_S, 'bc', 'Sets the pH of the water; loss of CO₂ on pressure drop drives carbonate scale.'),
    num('kScale', 'Scale surface-rate constant', 'kg/m²s', 2e-9, 0, 1e-4, G_S, 'material', 'Rate constant of the first-order growth law. Strongly system dependent — a user argument, not a default to rely on.'),
    num('dose', 'Scale inhibitor dose', 'mg/L', 10, 0, 500, G_S, 'operation', 'Inhibitor concentration in the water phase.'),
    num('mic', 'Minimum inhibitor concentration', 'mg/L', 5, 0, 500, G_S, 'criteria', 'Laboratory MIC for the worst mixing ratio.'),
    num('tMaxWax', 'Allowable wax thickness', 'mm', 2, 0.2, 20, G_C, 'criteria', 'Typical pigging trigger 2–4 mm.'),
    num('dpRatio', 'Allowable pressure-drop increase ratio', '-', 1.2, 1.02, 3, G_C, 'criteria', 'Pigging trigger on frictional pressure drop relative to the clean line.'),
    num('pigPlan', 'Planned pigging interval', 'days', 14, 1, 365, G_C, 'operation', 'Pigging frequency of the operating philosophy.'),
    num('pRestart', 'Restart pressure available', 'bar', 120, 1, 1400, G_C, 'criteria', 'Differential pressure the pumps or wells can apply to break the gel.'),
    { key: 'nSeg', label: 'Segments along the line', type: 'int', default: 40, min: 10, max: 400, group: 'Numerics', role: 'control', help: 'Axial segments of the deposition march.' },
  ],
  mesh: { key: 'nSeg', label: 'Segments along the line', levels: [20, 40, 80], qoi: 'thickMax' },
  compute(inp, ctx) {
    const p = pipeFrom({ ...inp, rhoCont: 0, steel: 'api5l-x65' }, ctx);
    const L = inp.Lkm * 1e3, D = p.Di;
    const th = thermalFrom(inp, p);
    const Uwall = inp.Ugiven > 0 ? inp.Ugiven * p.D / D : run(() => overallU({ Di: D, layers: th.layers, hi: Infinity, ho: th.hExt, Dref: D }), 'Wall conductance').U;
    const mDot = inp.qOil / DAY * inp.rhoOil;
    const WAT = inp.WAT + KELVIN;
    const sc = inp.solCurve && inp.solCurve.x && inp.solCurve.x.length >= 3 ? inp.solCurve : null;
    if (sc) for (let i = 1; i < sc.x.length; i++) need(sc.x[i] > sc.x[i - 1], 'The solubility curve temperatures must increase.');
    const sol = sc ? makeSolubilityCurve({ table: sc.x.map((T, i) => ({ T: T + KELVIN, C: sc.y[i] / 100 })), WAT }) : makeSolubilityCurve({ exponential: { a: inp.waxC / 100, b: inp.waxB, Tref: WAT }, WAT });
    const mu = (T) => inp.mu40 / 1e3 * Math.exp(inp.muB * (1 / T - 1 / 313.15));
    let watIdeal = null;
    try {
      const pr = (inp.paraffin || []).map((r) => ({ M: Number(r.M) / 1e3, z: Number(r.z) })).filter((r) => r.M > 0 && r.z > 0);
      const zt = pr.reduce((a, r) => a + r.z, 0);
      if (pr.length >= 2) watIdeal = waxAppearanceIdeal({ components: pr.map((r) => ({ M: r.M, z: r.z / zt })) });
    } catch { watIdeal = null; }
    const tEnd = inp.days * DAY;
    const args = { D, L, n: inp.nSeg, mDot, Tin: inp.Tin + KELVIN, Tamb: inp.Tamb + KELVIN, Uwall, rhoOil: inp.rhoOil, muOil: mu, cpOil: inp.cpOil, kOil: inp.kOil, solubility: sol, WAT,
      model: inp.waxModel, waxFraction: inp.waxFrac, kDeposit: inp.kDep, rough: inp.rough / 1e3, tEnd, dt: tEnd / 120, maxOutput: 121 };
    const w = run(() => waxDepositionLine(args), 'Wax deposition');
    const t0 = inp.thick0 / 1e3;
    const thick = w.thickness.map((v) => v + t0), tMaxH = w.thicknessMax.map((v) => v + t0);
    const thickMax = Math.max(...thick);
    const pig = run(() => piggingInterval({ t: w.t, thicknessMax: tMaxH, dP: w.dP, maxThickness: inp.tMaxWax / 1e3, dPratio: inp.dpRatio }), 'Pigging interval');
    const pigDays = Number.isFinite(pig.interval) ? pig.interval / DAY : 9999;
    const areaLoss = 1 - ((D - 2 * thickMax) / D) ** 2;
    const iWat = w.Tsurface.findIndex((T) => T < WAT);
    const xWat = iWat >= 0 ? w.x[iWat] : L;
    const dPclean = w.dP[0], dPend = w.dP[w.dP.length - 1];
    const Tarr = w.Tarrival[w.Tarrival.length - 1];
    const gel = gellingScreen({ Tmin: inp.Tamb + KELVIN, pourPoint: inp.pour + KELVIN });
    const dpGel = gelRestartPressure({ tauY: Math.max(inp.tauY, 1e-9), L: gel.gelRisk ? L : Math.max(L - xWat, 0) || L, D });
    // scale
    const mgL = {}; for (const r of inp.water || []) { const k = String(r.ion || '').trim(); if (k && Number(r.mgL) > 0) mgL[k] = Number(r.mgL); }
    need(mgL.Ca > 0 && mgL.HCO3 > 0, 'The water analysis needs at least calcium (Ca) and bicarbonate (HCO3) for the carbonate scale assessment.');
    const fw = run(() => toMolal({ mgL }), 'Formation water analysis (check the ion names)');
    const sw = toMolal({ mgL: Object.fromEntries(Object.entries(SEAWATER_MG_PER_KG).map(([k, v]) => [k, v * 1.025])) });
    const Tsc = inp.Tsc + KELVIN;
    const cal = run(() => calciteSaturation({ molal: fw.molal, T: Tsc, pCO2: inp.pCO2 * BAR }), 'Calcite saturation');
    const mix = run(() => seawaterMixingScale({ seawater: sw.molal, formationWater: fw.molal, T: Tsc, P: inp.Psc * BAR, n: 41 }), 'Seawater mixing');
    const bar = mix.results.barite, cel = mix.results.celestite, anh = mix.results.anhydrite;
    const worstSul = [['barite', bar, 4480], ['celestite', cel, 3960], ['anhydrite', anh, 2960]].filter((m) => m[1]).reduce((a, m) => (m[1].worstMass > a[1].worstMass ? m : a));
    const depCal = scaleDepositionRate({ SR: Math.max(cal.SR, 1e-9), k: inp.kScale, rhoScale: 2710, inhibitorDose: inp.dose, mic: inp.mic });
    const depCalRaw = scaleDepositionRate({ SR: Math.max(cal.SR, 1e-9), k: inp.kScale, rhoScale: 2710 });
    const scaling = cal.SI > 0 || (worstSul[1].maxSI > 0);
    const checks = [
      chk('Pigging interval against the plan', inp.pigPlan, Math.max(pigDays, 1e-6), 'days', 'Operating criterion: wax thickness or pressure-drop trigger',
        'Guards against wax building up beyond the thickness or pressure-drop limit between pig runs, after which the pig may stick or the line loses capacity.'),
      chk('Gel restart pressure', dpGel, inp.pRestart * BAR, 'Pa', 'Wall yield-stress force balance, ΔP = 4·τ_y·L/D',
        'Guards against a shut-in line gelling so firmly that the available pressure cannot restart flow; exceeding it can mean a permanently blocked line.'),
      chk('Scale inhibitor dose against the minimum inhibitor concentration', scaling ? inp.mic : 0, Math.max(inp.dose, 1e-9), 'mg/L', 'Scale management practice (laboratory MIC at the worst mixing ratio)',
        'Guards against mineral scale growing on the pipe wall, valves and screens when the water is supersaturated; below the MIC the inhibitor does not suppress crystal growth.'),
    ];
    const gov = worst(checks);
    const xs = w.x.map((v) => round(v / 1e3, 4)), days = w.t.map((v) => round(v / DAY, 3));
    const pigsPerYear = 365 / Math.max(Math.min(pigDays, 365), 1);
    const Dwo = diffusivityHaydukMinhas({ T: Math.min(WAT, inp.Tin + KELVIN), mu: mu(Math.min(WAT, inp.Tin + KELVIN)) });
    return {
      kpis: [
        kpi('thickMax', 'Maximum wax thickness', thickMax * 1e3, 'mm', 'Largest deposit thickness along the line at the end of the simulated period without pigging; compare with the allowable thickness that triggers a pig run.'),
        kpi('pigInterval', 'Pigging interval', pigDays, 'days', `Time for the deposit to reach the thickness or pressure-drop trigger (governed by ${pig.governing}); 9999 means no trigger is approached.`, inp.pigPlan, 'min'),
        kpi('waxVolume', 'Wax deposit volume', w.waxVolume[w.waxVolume.length - 1], 'm³', 'Volume of deposit in the line at the end of the period: the amount a pig run must remove.'),
        kpi('rate', 'Peak deposition rate', (thickMax - t0) * 1e3 / inp.days, 'mm/day', 'Average growth rate of the thickest part of the deposit over the simulated period.'),
        kpi('areaLoss', 'Flow-area reduction at the thickest point', areaLoss * 100, '%', 'Loss of bore cross-section where the deposit is thickest.'),
        kpi('dPratio', 'Pressure drop relative to the clean line', dPend / Math.max(dPclean, 1e-9), '-', 'Frictional pressure drop with the deposit over the clean-pipe value (smaller bore and rougher wall).', inp.dpRatio, 'max'),
        kpi('xWAT', 'Distance to the first deposition', xWat / 1e3, 'km', 'Chainage at which the deposit surface first falls below the wax appearance temperature. The line length means no deposition.'),
        kpi('Tarrival', 'Arrival temperature', Tarr - KELVIN, '°C', 'Oil temperature at the outlet at the end of the period; the insulating deposit raises it slightly.'),
        kpi('watIdeal', 'Wax appearance temperature, ideal-solution model', watIdeal ? watIdeal - KELVIN : inp.WAT, '°C', 'Solid–liquid equilibrium estimate from the n-paraffin distribution (Won fusion properties, ideal solutions). It typically over-predicts by 5–20 K; the measured value entered above is used for deposition.'),
        kpi('pRestart', 'Gel restart pressure', dpGel / BAR, 'bar', 'Upper-bound pressure difference to yield the gelled oil over the gelled length.', inp.pRestart, 'max'),
        kpi('siCalcite', 'Calcite saturation index', cal.SI, '-', 'log10 of the ion activity product over the solubility product. Positive means calcium carbonate scale can form; above about 0.5–1 it usually does.'),
        kpi('pH', 'Water pH at the stated CO₂ pressure', cal.pH, '-', 'pH of the formation water in equilibrium with the CO₂ partial pressure.'),
        kpi('siSulphate', `Maximum ${worstSul[0]} saturation index on mixing`, worstSul[1].maxSI, '-', 'Highest sulphate-scale saturation index over all sea-water mixing ratios.'),
        kpi('massSulphate', 'Sulphate scale potential at the worst mixing ratio', worstSul[1].worstMass, 'kg/m³', 'Mass of scale that must precipitate per cubic metre of mixed brine to return to saturation.'),
        kpi('mixWorst', 'Worst sea-water fraction', worstSul[1].worstFraction * 100, '%', 'Sea-water mass fraction in the produced water at which the sulphate scale potential peaks.'),
        kpi('scaleRate', 'Calcite growth rate without inhibitor', depCalRaw.thicknessRate * YEAR * 1e3, 'mm/year', 'Thickness growth from the first-order surface-reaction law with the user rate constant; zero when undersaturated.'),
      ],
      checks,
      plots: [
        plot('line', 'Wax thickness along the line', 'Chainage', 'km', 'Deposit thickness', 'mm', [
          { name: `After ${fmt(inp.days)} days`, x: xs, y: thick.map((v) => round(v * 1e3)) }, { name: 'Allowable', x: [0, round(L / 1e3)], y: [inp.tMaxWax, inp.tMaxWax], dash: true },
        ], 'Deposit profile: zero where the wall is above the wax appearance temperature, a peak just downstream of it where the temperature gradient is steepest, then thinning as the oil approaches ambient.'),
        plot('line', 'Temperatures along the line', 'Chainage', 'km', 'Temperature', '°C', [
          { name: 'Bulk oil', x: xs, y: w.Tbulk.map((v) => round(v - KELVIN)) }, { name: 'Deposit surface', x: xs, y: w.Tsurface.map((v) => round(v - KELVIN)) },
          { name: 'Wax appearance temperature', x: [0, round(L / 1e3)], y: [inp.WAT, inp.WAT], dash: true },
        ], 'Deposition occurs where the surface temperature is below the dashed wax appearance temperature while the bulk is still warmer than the surface.'),
        plot('line', 'Deposit growth and pressure drop', 'Time', 'days', 'Maximum thickness', 'mm', [
          { name: 'Maximum thickness', x: days, y: tMaxH.map((v) => round(v * 1e3)) }, { name: 'Pressure-drop ratio × 10', x: days, y: Array.from(w.dP, (v) => round(10 * v / Math.max(dPclean, 1e-9))) },
        ], 'Growth of the thickest deposit with time and the resulting rise in frictional pressure drop. Growth slows as the deposit insulates the wall.'),
        plot('line', 'Sulphate scale on mixing with sea water', 'Sea-water fraction', '%', 'Scale potential', 'kg/m³', [
          ...(bar ? [{ name: 'Barite', x: mix.fraction.map((f) => round(f * 100)), y: bar.mass.map((v) => round(v)) }] : []),
          ...(cel ? [{ name: 'Celestite', x: mix.fraction.map((f) => round(f * 100)), y: cel.mass.map((v) => round(v)) }] : []),
          ...(anh ? [{ name: 'Anhydrite', x: mix.fraction.map((f) => round(f * 100)), y: anh.mass.map((v) => round(v)) }] : []),
        ], 'Mass of sulphate scale that can precipitate at each mixing ratio of injected sea water and formation water. The peak is the breakthrough condition to inhibit against.'),
      ],
      tables: [tbl('Scale tendency', [col('mineral', 'Mineral'), col('si', 'Saturation index', '-'), col('mass', 'Scale potential', 'kg/m³'), col('at', 'At')], [
        { mineral: 'Calcite (CaCO₃)', si: round(cal.SI, 3), mass: 0, at: `formation water, ${fmt(inp.pCO2)} bar CO₂, pH ${fmt(cal.pH)}` },
        ...[['Barite (BaSO₄)', bar], ['Celestite (SrSO₄)', cel], ['Anhydrite (CaSO₄)', anh]].filter((m) => m[1]).map(([n, m]) => ({ mineral: n, si: round(m.maxSI, 3), mass: round(m.worstMass, 3), at: `${fmt(m.worstFraction * 100)} % sea water` })),
      ], 'One row per mineral: the saturation index (positive = scaling) and the mass of scale available per cubic metre at the condition where it is worst.')],
      notes: [`Wall conductance ${fmt(Uwall)} W/m²K on the bore; wax diffusivity ${fmt(Dwo * 1e10)}·10⁻¹⁰ m²/s (Hayduk–Minhas); deposit starts at ${fmt(xWat / 1e3)} km.`,
        gel.gelRisk ? `Seabed temperature is within 5 K of the pour point or below it: the line can gel on shutdown (margin ${fmt(gel.margin)} K).` : `Seabed temperature is ${fmt(gel.margin)} K above the pour point: gelling on shutdown is unlikely.`,
        `Formation water: TDS ${fmt(fw.tds / 1e3)} g/L, charge balance ${fmt(fw.chargeBalance * 100)} %, ionic strength ${fmt(cal.I)} mol/kg. Scale kinetics use a user rate constant.`],
      advice: [
        adv('engineering', sev(gov.uc), `Wax reaches ${fmt(thickMax * 1e3)} mm in ${fmt(inp.days)} days; the pigging trigger is reached after ${fmt(pigDays)} days (plan ${fmt(inp.pigPlan)}); gel restart needs ${fmt(dpGel / BAR)} bar; calcite SI ${fmt(cal.SI)}.`,
          gov.uc > 1 ? (gov.name.startsWith('Pigging') ? `Pig at least every ${fmt(pigDays)} days, or insulate so the wall stays above ${fmt(inp.WAT)} °C over more of the line.` : gov.name.startsWith('Gel') ? `Displace the line with dead oil or diesel before a planned shutdown, inject pour-point depressant, or provide ${fmt(dpGel / BAR)} bar of restart pressure.` : `Raise the scale inhibitor dose above the MIC of ${fmt(inp.mic)} mg/L and confirm it at the worst mixing ratio of ${fmt(worstSul[1].worstFraction * 100)} % sea water.`)
            : `The planned ${fmt(inp.pigPlan)}-day pigging holds the deposit below ${fmt(inp.tMaxWax)} mm. Confirm the solubility curve and deposit wax fraction with laboratory data; they control the growth rate.`,
          'Burger et al. (1981); Singh et al. (2000); flow-assurance practice per ISO 13628-1 / API RP 17A'),
        adv('economic', 'info', `The deposit growth implies about ${fmt(pigsPerYear)} pig runs per year and ${fmt(w.waxVolume[w.waxVolume.length - 1])} m³ of wax to handle per run after ${fmt(inp.days)} days.`,
          pigDays < inp.pigPlan ? 'More frequent pigging costs operator time and deferred production per run; compare with wax inhibitor dosing or added insulation over the field life.' : `The pigging interval could be extended towards ${fmt(Math.min(pigDays, 365))} days, cutting the number of runs and the associated production interruptions.`,
          ECON_BASIS, `Each pig run from a support vessel is of the order of a vessel day. ${vesselCo2(pigsPerYear)} Wax and scale inhibitors are continuous chemical discharges in the produced water; thermal management avoids them.`),
      ],
      publish: { 'pipeline.piggingInterval': Math.min(pigDays, 9999) * DAY, 'pipeline.waxThickness': thickMax },
    };
  },
  verify: [{
    // Gel restart: ΔP = 4 τ_y L/D with τ_y = 50 Pa, L = 10 km, D_i = 323.9 − 2·14.3 = 295.3 mm, whole line gelled (seabed below the pour point):
    // ΔP = 4·50·10 000/0.2953 = 6.773 MPa = 67.73 bar.
    name: 'Gel restart pressure, ΔP = 4·τ_y·L/D',
    source: 'Force balance on a gelled plug (Davidson et al. 2004); hand calculation',
    inputs: { D: 323.9, t: 14.3, Lkm: 10, tauY: 50, pour: 20, Tamb: 4 },
    expect: [{ kpi: 'pRestart', value: 67.73, tol: 0.002 }],
  }, {
    // Wall everywhere above the wax appearance temperature (WAT below ambient): the oil stays undersaturated and nothing deposits.
    name: 'No deposition when the wall stays above the wax appearance temperature',
    source: 'Limiting case of the molecular-diffusion mechanism (no supersaturation at the wall)',
    inputs: { WAT: 0, Tamb: 4, Tin: 60 },
    expect: [{ kpi: 'thickMax', value: 0, abs: 1e-9 }],
  }],
};
