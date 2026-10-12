// Corrosion, erosion and wall-thickness degradation: CO2 corrosion rate along the line by the de Waard family of
// models with inhibitor availability, sand erosion in bends and straight pipe, erosional and sand-transport velocity
// limits, external oxygen corrosion at coating damage from electrochemical kinetics, and the wall-thickness
// projection over the design life.
import { deWaardMilliams, deWaardLotzMilliams, deWaardLotzDugstad, phCO2Saturated, inhibitedRate, corrosionAllowance, sourServiceScreening,
  massTransferBergerHau, oxygenCorrosionRate, mixedPotential, polarisationCurve, tafel } from '../../engine/flow/corrosion.js';
import { anodeDesign, coatingBreakdownConstants } from '../../engine/fatigue/degradation.js';
import { erosionalVelocity, erosionStraightPipe, erosionBend, erosionSalama, maxVelocityForErosion, criticalVelocityOroskarTurian } from '../../engine/flow/erosion.js';
import { pipeInputs, pipeFrom, lengthInput, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, sev, worst, linspace,
  steelPerMetre, steelCo2, CONFIRM, RBI_BASIS, CO2_STEEL, BAR, KELVIN } from './_common.js';

const G_P = 'Process conditions';
const G_C = 'CO₂ corrosion';
const G_E = 'Sand and erosion';
const G_X = 'External corrosion (electrochemical model)';
const G_A = 'Acceptance criteria';
const G_CP = 'External cathodic protection';

export default {
  id: 's08-corrosion-erosion',
  suite: 8,
  title: 'Corrosion, erosion and wall-thickness degradation',
  summary: 'Estimates internal CO₂ corrosion and sand erosion rates along the line, the corrosion allowance they consume, and when the wall reaches its minimum thickness.',
  covers: [16],
  fidelity: 'screening',
  standards: ['de Waard & Milliams (1975); de Waard, Lotz & Milliams (1991); de Waard, Lotz & Dugstad (1995)', 'NORSOK M-001 (2014) (inhibitor availability, corrosion allowance)', 'ISO 15156-2 / NACE MR0175 (sour service limit)', 'DNV-RP-O501 (2015) sand erosion', 'API RP 14E (1991, R2019) §2.5 erosional velocity', 'DNV-RP-F103 / DNV-RP-B401 (external corrosion context)'],
  models: ['Carbon dioxide corrosion correlations (de Waard family) with scale, pH and glycol factors', 'Corrosion mass-transfer resistance model', 'Inhibitor availability model', 'Particle-impact erosion model for bends and straight pipe', 'Erosional velocity correlation (API RP 14E) and Salama correlation', 'Sand transport critical velocity (Oroskar–Turian)', 'Butler–Volmer / Tafel kinetics and mixed-potential theory with a diffusion-limited cathode', 'Galvanic anode cathodic protection design (current demand, anode mass and resistance)', 'Linear wall-thickness deterioration'],
  limits: 'Empirical CO₂ corrosion correlations for sweet service on bare carbon steel: not valid in the H₂S-dominated regime, with organic acids, or under protective films other than the built-in scaling factor; oil wetting and top-of-line condensation are not modelled. Erosion models are for ductile steel with dilute sand and carry a factor of about three in uncertainty. Electrochemical parameters are user inputs with illustrative defaults. Localised (pitting, mesa, weld) attack, microbiological corrosion and CFD particle tracking are outside the module.',
  publishes: ['pipeline.corrosionRate', 'pipeline.wallLife'],
  searchTerms: 'CO2 corrosion rate de Waard pipeline sand erosion DNV-RP-O501 corrosion allowance',
  inputs: [
    ...pipeInputs({ omit: ['tCoat', 'rhoCoat', 'tConc', 'concrete', 'rhoCont', 'wsGiven', 'rhoW', 'steel'] }),
    lengthInput(10, 'Pipe'),
    num('pIn', 'Inlet pressure', 'bar', 100, 1, 1400, G_P, 'bc', 'Total pressure at the inlet.', { bind: 'pipeline.pIn' }),
    num('pOut', 'Outlet pressure', 'bar', 60, 1, 1400, G_P, 'bc', 'Total pressure at the outlet.'),
    num('Tin', 'Inlet temperature', '°C', 70, 0, 150, G_P, 'bc', 'Fluid temperature at the inlet.', { uncertain: { dist: 'normal', cov: 0.05 } }),
    num('Tout', 'Outlet temperature', '°C', 30, 0, 150, G_P, 'bc', 'Fluid temperature at the outlet.', { bind: 'pipeline.Tout' }),
    num('vLiq', 'Liquid velocity', 'm/s', 2, 0.05, 20, G_P, 'operation', 'Actual liquid velocity (for mass transfer and sand transport).', { uncertain: { dist: 'normal', cov: 0.15 } }),
    num('vMix', 'Mixture velocity', 'm/s', 6, 0.1, 80, G_P, 'operation', 'Gas plus liquid superficial velocity at the worst (low-pressure) location, for erosion.', { bind: 'pipeline.vMix', design: true, uncertain: { dist: 'normal', cov: 0.15 } }),
    num('rhoMix', 'Mixture density', 'kg/m³', 250, 5, 1100, G_P, 'operation', 'No-slip mixture density at the erosion location.'),
    num('muMix', 'Mixture viscosity', 'cP', 0.3, 0.01, 100, G_P, 'operation', 'Mixture viscosity at the erosion location.'),
    num('rhoLiq', 'Liquid density', 'kg/m³', 850, 500, 1200, G_P, 'material', 'For the sand transport limit.'),
    num('muLiq', 'Liquid viscosity', 'cP', 2, 0.1, 1000, G_P, 'material', 'For the sand transport limit.'),
    num('yCO2', 'CO₂ in the gas', 'mol %', 0.5, 0, 50, G_C, 'material', 'Mole fraction of CO₂; partial pressure = fraction × total pressure.', { uncertain: { dist: 'normal', cov: 0.1 } }),
    num('yH2S', 'H₂S in the gas', 'ppm (mol)', 20, 0, 200000, G_C, 'material', 'For sour-service and corrosion-regime screening.'),
    sel('model', 'CO₂ corrosion model', 'dwld', [['dwld', 'de Waard, Lotz & Dugstad (1995), flow dependent'], ['dwlm', 'de Waard, Lotz & Milliams (1991), scale and pH factors'], ['dwm', 'de Waard & Milliams (1975), worst-case nomogram']], G_C, 'control', 'The 1995 model adds a mass-transfer resistance that depends on velocity; the 1975 nomogram is an upper bound.'),
    num('pH', 'Actual water pH (0 = CO₂-saturated condensed water)', '-', 5.5, 0, 8, G_C, 'material', 'Formation water with bicarbonate has a higher pH, which lowers the rate.'),
    num('glycolWater', 'Water content of the aqueous phase', 'wt %', 100, 5, 100, G_C, 'operation', '100 = no glycol; MEG-rich phases reduce the rate (1991 model factor).'),
    num('avail', 'Inhibitor availability', '-', 0.97, 0, 1, G_C, 'operation', 'Fraction of time the corrosion inhibitor is at the required dose (0 = no inhibition; 0.90–0.95 typical, 0.99 with redundancy).'),
    num('crInh', 'Inhibited corrosion rate', 'mm/year', 0.1, 0.01, 1, G_C, 'operation', 'Residual rate while the inhibitor is effective (typically 0.1 mm/year).'),
    num('sand', 'Sand production', 'kg/day', 5, 0, 5000, G_E, 'bc', 'Sand mass rate carried by the stream.', { uncertain: { dist: 'lognormal', cov: 0.5 } }),
    num('dp', 'Sand particle size', 'µm', 250, 20, 2000, G_E, 'material', 'Median particle diameter.'),
    num('Rc', 'Bend radius / pipe diameter', '-', 5, 0.5, 50, G_E, 'geometry', '1.5 for standard elbows, 5 for induction bends.'),
    num('GF', 'Geometry factor', '-', 1, 1, 4, G_E, 'geometry', '1 for an isolated bend, 2–3 for bends in series (DNV-RP-O501).'),
    num('erosC', 'API RP 14E constant C', '-', 100, 50, 300, G_E, 'criteria', '100 continuous solids-free service; lower with sand.'),
    num('cO2', 'Dissolved oxygen in sea water', 'ppm', 7, 0, 12, G_X, 'environment', 'Oxygen available at a coating holiday without cathodic protection.'),
    num('vSea', 'Sea-water velocity at the pipe', 'm/s', 0.3, 0.01, 3, G_X, 'environment', 'Sets the oxygen mass-transfer coefficient.'),
    num('Ea', 'Anodic equilibrium potential', 'V', -0.62, -1.2, 0.2, G_X, 'material', `Iron dissolution reference potential (illustrative). ${CONFIRM}`),
    num('i0a', 'Anodic exchange current density', 'A/m²', 1e-3, 1e-7, 10, G_X, 'material', 'Illustrative default; fit to polarisation data.'),
    num('ba', 'Anodic Tafel slope', 'V/decade', 0.06, 0.02, 0.2, G_X, 'material', 'Typically 0.04–0.08 V/decade for iron.'),
    num('Ec', 'Cathodic equilibrium potential', 'V', 0.6, -0.5, 1.3, G_X, 'material', 'Oxygen reduction reference potential at sea-water pH (illustrative).'),
    num('i0c', 'Cathodic exchange current density', 'A/m²', 1e-6, 1e-10, 1, G_X, 'material', 'Illustrative default; fit to polarisation data.'),
    num('bc', 'Cathodic Tafel slope', 'V/decade', 0.12, 0.03, 0.3, G_X, 'material', 'About 0.12 V/decade for oxygen reduction.'),
    num('damage', 'Coating damage fraction without effective cathodic protection', '%', 0, 0, 100, G_X, 'ic', 'Share of the external surface corroding freely; zero for an intact coating with working CP.'),
    sel('coatCat', 'External coating category', 'III', [['I', 'I — one-layer epoxy'], ['II', 'II — multi-layer marine coating'], ['III', 'III — high-build / 3-layer polyolefin'], ['bare', 'Bare steel']], G_CP, 'material', `Sets the DNV-RP-B401 coating breakdown constants; pipeline coatings with field joints are treated in DNV-RP-F103 / ISO 15589-2. ${CONFIRM}`),
    num('icMean', 'Mean design current density', 'A/m²', 0.07, 0.01, 0.5, G_CP, 'environment', `Bare-steel current demand averaged over the life; depends on depth, temperature and region (DNV-RP-B401 tables, about 0.06–0.22 A/m²). ${CONFIRM}`),
    num('icFinal', 'Initial / final design current density', 'A/m²', 0.12, 0.01, 0.6, G_CP, 'environment', `Used for both the polarisation (initial) and re-polarisation (final) checks. ${CONFIRM}`),
    num('anodeMass', 'Net mass of one bracelet anode', 'kg', 60, 2, 2000, G_CP, 'geometry', 'Aluminium–zinc–indium alloy mass per bracelet.'),
    num('anodeArea', 'Exposed area of one bracelet anode', 'm²', 0.5, 0.02, 10, G_CP, 'geometry', 'Outer surface of the bracelet, for the anode resistance.'),
    num('rhoSea', 'Sea-water resistivity', 'Ω·m', 0.3, 0.15, 5, G_CP, 'environment', 'About 0.2–0.35 Ω·m in the open sea; higher in cold or brackish water and much higher in sediment.'),
    num('anodeSpacingMax', 'Maximum anode spacing', 'm', 300, 12, 2000, G_CP, 'criteria', `Attenuation limit on the distance between bracelets (ISO 15589-2 gives 300 m as a default maximum). ${CONFIRM}`),
    num('tCorr', 'Corrosion allowance provided', 'mm', 4, 0, 15, G_A, 'geometry', 'Wall thickness reserved for metal loss.', { design: true }),
    num('lossNow', 'Wall loss to date', 'mm', 0, 0, 30, G_A, 'ic', 'Metal loss already measured by inspection.'),
    num('life', 'Remaining design life', 'years', 25, 1, 100, G_A, 'operation', 'Period over which the wall must stay above the minimum.', { bind: 'life.designYears' }),
    num('erAllow', 'Allowable erosion rate', 'mm/year', 0.1, 0.005, 2, G_A, 'criteria', 'Company limit for erosive wear (often 0.05–0.1 mm/year).'),
    num('maxCA', 'Largest practicable corrosion allowance', 'mm', 8, 3, 15, G_A, 'criteria', 'Above this a corrosion-resistant alloy or cladding is normally selected (6–10 mm).'),
  ],
  compute(inp, ctx) {
    const p = pipeFrom({ ...inp, rhoCont: 0, steel: 'api5l-x65' }, ctx);
    const D = p.Di, L = inp.Lkm * 1e3;
    const n = 21, xs = linspace(0, L, n);
    const pts = xs.map((x) => { const f = x / L; return { x, P: (inp.pIn + (inp.pOut - inp.pIn) * f) * BAR, T: inp.Tin + (inp.Tout - inp.Tin) * f + KELVIN }; });
    const rateAt = ({ P, T }, vel = inp.vLiq) => {
      const pCO2 = Math.max(inp.yCO2 / 100 * P, 1);
      const pH = inp.pH > 0 ? inp.pH : undefined;
      if (inp.model === 'dwm') return { rate: deWaardMilliams({ T, pCO2 }), pCO2 };
      if (inp.model === 'dwlm') { const r = deWaardLotzMilliams({ T, pCO2, P, pH, waterInGlycol: inp.glycolWater }); return { rate: r.rate, pCO2, Tscale: r.Tscale, Fscale: r.Fscale }; }
      const gly = inp.glycolWater < 100 ? 10 ** (1.6 * (Math.log10(inp.glycolWater) - 2)) : 1;
      const r = deWaardLotzDugstad({ T, pCO2, U: vel, d: D, P, pH, extraFactor: gly });
      return { rate: r.rate, pCO2, Vr: r.Vr, Vm: r.Vm };
    };
    need(inp.yCO2 > 0 || inp.sand > 0 || inp.damage > 0, 'With no CO₂, no sand and no external exposure there is no degradation mechanism to assess.');
    const prof = pts.map((q) => (inp.yCO2 > 0 ? run(() => rateAt(q), 'CO₂ corrosion rate') : { rate: 0, pCO2: 0 }));
    const iMax = prof.reduce((a, r, i) => (r.rate > prof[a].rate ? i : a), 0);
    const crUn = prof[iMax].rate;
    const crInt = inp.yCO2 > 0 ? inhibitedRate({ uninhibited: crUn, inhibited: Math.min(inp.crInh, crUn), availability: inp.avail }) : 0;
    const pHsat = inp.yCO2 > 0 ? phCO2Saturated({ T: pts[iMax].T, pCO2: prof[iMax].pCO2 }) : 7;
    const sour = sourServiceScreening({ pH2S: inp.yH2S * 1e-6 * inp.pIn * BAR, pCO2: Math.max(inp.yCO2 / 100 * inp.pIn * BAR, 1e-9) });
    // erosion
    const mSand = inp.sand / 86400, dp = inp.dp * 1e-6, muM = inp.muMix / 1e3;
    const eBend = mSand > 0 ? run(() => erosionBend({ Up: inp.vMix, D, mSand, dp, rhoM: inp.rhoMix, muM, Rc: inp.Rc, GF: inp.GF }), 'Bend erosion').rate : 0;
    const eStr = mSand > 0 ? erosionStraightPipe({ Up: inp.vMix, D, mSand }) : 0;
    const eSal = mSand > 0 ? erosionSalama({ mSand, Vm: inp.vMix, dp, D, rhoM: inp.rhoMix }) : 0;
    const vEros = erosionalVelocity({ rho: inp.rhoMix, C: inp.erosC });
    const vMaxE = mSand > 0 ? maxVelocityForErosion({ model: 'bend', allowable: inp.erAllow, sandConcentration: mSand / (inp.vMix * Math.PI * D * D / 4), D, dp, rhoM: inp.rhoMix, muM, Rc: inp.Rc, GF: inp.GF }).v : 300;
    const Cs = mSand > 0 ? Math.max(mSand / 2650 / (inp.vLiq * Math.PI * D * D / 4), 1e-7) : 0;
    const vCrit = mSand > 0 ? criticalVelocityOroskarTurian({ d: dp, D, C: Math.min(Math.max(Cs, 1e-6), 0.3), rhoL: inp.rhoLiq, muL: inp.muLiq / 1e3 }) : 0;
    // external corrosion at coating damage: oxygen reduction under mass-transfer control (Berger & Hau on the outer diameter)
    const nuSea = 1.4e-6, Dox = 1.6e-9;
    const mt = massTransferBergerHau({ Re: Math.max(inp.vSea * p.D / nuSea, 1e3), Sc: nuSea / Dox, Ddiff: Dox, d: p.D });
    const ox = oxygenCorrosionRate({ cO2: inp.cO2 * 1e-3, km: mt.km });
    const anodic = { E0: inp.Ea, i0: inp.i0a, b: inp.ba }, cathodic = { E0: inp.Ec, i0: inp.i0c, b: inp.bc, iLim: Math.max(ox.iLim, 1e-9) };
    const mp = run(() => mixedPotential({ anodic, cathodic }), 'Mixed-potential solution');
    const pc = polarisationCurve({ anodic, cathodic, n: 61 });
    const crExt = mp.rate * inp.damage / 100;
    // external cathodic protection by bracelet anodes over the whole length
    const cb = coatingBreakdownConstants(inp.coatCat, 100);
    const cp = run(() => anodeDesign({ zones: [{ area: Math.PI * p.D * L, initialDensity: inp.icFinal, meanDensity: inp.icMean, finalDensity: inp.icFinal, a: cb.a, b: cb.b }], designLife: inp.life,
      anode: { type: 'bracelet', netMass: inp.anodeMass, area: inp.anodeArea, length: 0.5 }, resistivity: inp.rhoSea }), 'Cathodic protection design');
    const nAnode = Math.max(cp.anodes, Math.ceil(L / inp.anodeSpacingMax));
    const spacing = L / nAnode;
    // wall projection at the worst location
    const total = crInt + eBend + crExt;
    const ca = corrosionAllowance({ rate: Math.max(crInt + crExt, 1e-9), life: inp.life, maxAllowance: inp.maxCA, available: Math.max(inp.tCorr - inp.lossNow, 0) });
    const need_ = total * inp.life + inp.lossNow;
    const lifeWall = total > 0 ? Math.max(inp.tCorr - inp.lossNow, 0) / total : 1e6;
    const yrs = linspace(0, Math.max(inp.life, Math.min(lifeWall * 1.2, 100)), 26);
    const tMin = inp.t - inp.tCorr;
    const checks = [
      chk('Metal loss over the design life against the corrosion allowance', need_, inp.tCorr, 'mm', 'NORSOK M-001 (2014) §4 corrosion allowance; DNV-ST-F101 (2021) §6',
        'Guards against internal corrosion, erosion and external attack consuming more wall than was reserved, after which the pressure and collapse checks no longer hold.'),
      chk('Sand erosion rate in bends', eBend, inp.erAllow, 'mm/year', 'DNV-RP-O501 (2015) pipe-bend model',
        'Guards against sand cutting through the outer wall of bends, the highest-wear location; the model is conservative by about a factor of three.'),
      chk('Mixture velocity against the API RP 14E erosional velocity', inp.vMix, vEros, 'm/s', 'API RP 14E (1991, R2019) §2.5.a',
        'Empirical velocity limit above which erosion-corrosion and noise become likely; it contains no sand rate and is a screening value only.'),
    ];
    if (mSand > 0) checks.push(chk('Sand transport: liquid velocity above the settling limit', vCrit, inp.vLiq, 'm/s', 'Oroskar & Turian (1980)',
      'Guards against sand settling into a stationary bed, which promotes under-deposit corrosion and pigging problems.'));
    const gov = worst(checks);
    const km = xs.map((x) => round(x / 1e3, 4));
    const us = linspace(0.3, Math.max(3 * inp.vLiq, 6), 20);
    const crU = inp.yCO2 > 0 && inp.model === 'dwld' ? us.map((u) => rateAt(pts[iMax], u).rate) : us.map(() => crUn);
    const dCA = Math.max(need_ - inp.tCorr, 0);
    const Lr = inp.Lkm;
    return {
      kpis: [
        kpi('crUninhibited', 'Uninhibited CO₂ corrosion rate (worst location)', crUn, 'mm/year', 'Predicted bare-steel corrosion rate at the worst point along the line with no inhibitor.'),
        kpi('crInternal', 'Time-averaged internal corrosion rate', crInt, 'mm/year', 'Rate with the inhibitor availability model: inhibited rate while the inhibitor works, uninhibited rate for the rest of the time.'),
        kpi('erosionBend', 'Sand erosion rate in bends', eBend, 'mm/year', 'Metal loss rate at the bend extrados from particle impact (DNV-RP-O501).', inp.erAllow, 'max'),
        kpi('erosionStraight', 'Sand erosion rate in straight pipe', eStr, 'mm/year', 'Erosion of straight runs; normally negligible beside bends.'),
        kpi('erosionSalama', 'Bend erosion by the Salama correlation', eSal, 'mm/year', 'Independent empirical estimate for multiphase flow, for comparison.'),
        kpi('crExternal', 'External corrosion contribution', crExt, 'mm/year', 'Free-corrosion rate at coating damage (oxygen diffusion limited) weighted by the unprotected area fraction.'),
        kpi('rateTotal', 'Total wall-loss rate at the worst location', total, 'mm/year', 'Sum of internal corrosion, bend erosion and external contribution.'),
        kpi('caRequired', 'Corrosion allowance required', need_, 'mm', 'Wall loss to date plus the predicted loss over the remaining design life.', inp.tCorr, 'max'),
        kpi('wallLife', 'Time to consume the corrosion allowance', Math.min(lifeWall, 1e6), 'years', 'Years until the remaining allowance is used up at the predicted rate.', inp.life, 'min'),
        kpi('tEnd', 'Wall thickness at the end of life', Math.max(inp.t - need_, 0), 'mm', 'Remaining wall at the worst location after the design life.', tMin, 'min'),
        kpi('vErosional', 'API RP 14E erosional velocity', vEros, 'm/s', 'C/√ρ in SI units for the stated mixture density.'),
        kpi('vMaxErosion', 'Velocity for the allowable erosion rate', Math.min(vMaxE, 300), 'm/s', 'Mixture velocity at which the bend erosion reaches the allowable rate at the same sand concentration (300 = not reached).'),
        kpi('vCritSand', 'Critical velocity for sand transport', vCrit, 'm/s', 'Liquid velocity below which sand settles into a bed.'),
        kpi('pCO2', 'CO₂ partial pressure at the inlet', inp.yCO2 / 100 * inp.pIn, 'bar', 'Driving quantity of sweet corrosion.'),
        kpi('pHsat', 'pH of CO₂-saturated condensed water', pHsat, '-', 'pH of pure water in equilibrium with the CO₂ at the worst location; formation water with bicarbonate is higher.'),
        kpi('anodes', 'Bracelet anodes required', nAnode, '-', 'Number of bracelets from the larger of the mass, current and maximum-spacing requirements over the design life.'),
        kpi('anodeSpacing', 'Anode spacing', spacing, 'm', 'Distance between bracelets along the line.', inp.anodeSpacingMax, 'max'),
        kpi('anodeMassTotal', 'Total net anode mass', nAnode * inp.anodeMass / 1e3, 't', 'Aluminium alloy installed for external cathodic protection.'),
        kpi('cpCurrent', 'Mean protection current demand', cp.currentMean, 'A', 'Life-averaged current the anodes must deliver to the coating defects.'),
        kpi('Ecorr', 'External free-corrosion potential', mp.Ecorr, 'V', 'Mixed potential where iron dissolution balances oxygen reduction at a coating holiday.'),
        kpi('crFree', 'External free-corrosion rate of bare steel', mp.rate, 'mm/year', 'Penetration rate of unprotected bare steel in the stated sea water from the mixed-potential solution.'),
      ],
      checks,
      plots: [
        plot('line', 'CO₂ corrosion rate along the line', 'Chainage', 'km', 'Corrosion rate', 'mm/year', [
          { name: 'Uninhibited', x: km, y: prof.map((r) => round(r.rate)) }, { name: 'Time-averaged with inhibitor', x: km, y: prof.map((r) => round(inp.yCO2 > 0 ? inhibitedRate({ uninhibited: r.rate, inhibited: Math.min(inp.crInh, r.rate), availability: inp.avail }) : 0)) },
          { name: 'Temperature (°C ÷ 20)', x: km, y: pts.map((q) => round((q.T - KELVIN) / 20)), dash: true },
        ], 'Predicted rate from inlet to outlet. It usually peaks where temperature is 60–80 °C and falls towards the cold end; the location of the peak is where to place corrosion monitoring.'),
        plot('line', 'Wall thickness projection', 'Time', 'years', 'Wall thickness at the worst location', 'mm', [
          { name: 'Predicted wall', x: yrs.map((y) => round(y)), y: yrs.map((y) => round(Math.max(inp.t - inp.lossNow - total * y, 0))) },
          { name: 'Minimum wall (nominal − allowance)', x: [0, round(yrs[yrs.length - 1])], y: [round(tMin), round(tMin)], dash: true },
        ], 'Linear projection of the wall at the worst location. Where the line crosses the dashed minimum the corrosion allowance is used up.'),
        plot('line', 'Corrosion rate vs liquid velocity', 'Liquid velocity', 'm/s', 'Uninhibited rate', 'mm/year', [
          { name: 'Worst location', x: us.map((u) => round(u)), y: crU.map((v) => round(v)) }, { name: 'Operating point', x: [inp.vLiq], y: [round(crUn)] },
        ], 'Flow dependence of the corrosion rate from the mass-transfer resistance (1995 model). A flat line means the selected model has no velocity term.'),
        plot('line', 'Polarisation diagram at a coating holiday', 'Potential', 'V', 'log10 current density (A/m²)', '-', [
          { name: 'Anodic (iron dissolution)', x: pc.E.map((v) => round(v)), y: pc.iAnodic.map((v) => round(Math.log10(Math.max(v, 1e-12)))) },
          { name: 'Cathodic (oxygen reduction)', x: pc.E.map((v) => round(v)), y: pc.iCathodic.map((v) => round(Math.log10(Math.max(v, 1e-12)))) },
          { name: 'Corrosion point', x: [round(mp.Ecorr)], y: [round(Math.log10(Math.max(mp.icorr, 1e-12)))] },
        ], 'Evans diagram: the crossing of the anodic and cathodic curves is the corrosion potential and current. A flat cathodic branch shows the rate is limited by oxygen supply, not by kinetics.'),
      ],
      tables: [tbl('Degradation summary', [col('mech', 'Mechanism'), col('rate', 'Rate', 'mm/year'), col('loss', 'Loss over the design life', 'mm'), col('basis', 'Basis')], [
        { mech: 'Internal CO₂ corrosion', rate: round(crInt), loss: round(crInt * inp.life), basis: `${inp.model} model, availability ${fmt(inp.avail)}` },
        { mech: 'Sand erosion (bends)', rate: round(eBend), loss: round(eBend * inp.life), basis: `DNV-RP-O501, ${fmt(inp.sand)} kg/day, ${fmt(inp.dp)} µm` },
        { mech: 'External corrosion at coating damage', rate: round(crExt), loss: round(crExt * inp.life), basis: `${fmt(inp.damage)} % unprotected area, i_lim ${fmt(ox.iLim)} A/m²` },
        { mech: 'Total', rate: round(total), loss: round(total * inp.life), basis: 'sum at the worst location' },
      ], 'Each degradation mechanism with its rate and the wall it consumes over the remaining design life, so the dominant one is clear.')],
      notes: [`Cathodic protection: ${nAnode} bracelets at ${fmt(spacing)} m, governed by ${Math.ceil(L / inp.anodeSpacingMax) > cp.anodes ? 'the maximum spacing' : cp.governing}; coating breakdown a = ${fmt(cb.a)}, b = ${fmt(cb.b)} per year.`,
        `Worst internal corrosion at ${fmt(xs[iMax] / 1e3)} km (${fmt(pts[iMax].T - KELVIN)} °C, p_CO2 ${fmt(prof[iMax].pCO2 / BAR)} bar).`,
        sour.sourService ? `H₂S partial pressure exceeds 0.3 kPa: sour service — materials must meet ISO 15156 and the corrosion regime is ${sour.regime}; the CO₂ models lose validity as H₂S dominates.` : `Sweet service (H₂S below the ISO 15156 limit); corrosion regime ${sour.regime}.`,
        `Cathodic Tafel slope at 25 °C for α·n = 0.5 would be ${fmt(tafel({}).b)} V/decade; the electrochemical inputs are illustrative.`],
      advice: [
        adv('engineering', sev(gov.uc), `Predicted wall loss ${fmt(total)} mm/year at the worst location (${fmt(crInt)} corrosion, ${fmt(eBend)} erosion, ${fmt(crExt)} external): ${fmt(need_)} mm over ${fmt(inp.life)} years against ${fmt(inp.tCorr)} mm allowance.`,
          gov.uc > 1 ? (gov.name.startsWith('Metal') ? (need_ > inp.maxCA ? `An allowance of ${fmt(need_)} mm is beyond the practicable ${fmt(inp.maxCA)} mm: raise inhibitor availability, or select a corrosion-resistant alloy or clad pipe.` : `Increase the corrosion allowance to ${fmt(need_)} mm or raise the inhibitor availability above ${fmt(inp.avail)}.`)
            : gov.name.startsWith('Sand transport') ? `Keep the liquid velocity above ${fmt(vCrit)} m/s or pig regularly to remove settled sand.` : `Limit the mixture velocity to ${fmt(Math.min(vMaxE, vEros))} m/s, use longer-radius bends, or reduce sand at source.`)
            : `The allowance lasts ${fmt(lifeWall)} years. Place corrosion monitoring near ${fmt(xs[iMax] / 1e3)} km and confirm the rate with coupons and the first in-line inspection.`,
          'NORSOK M-001 (2014); de Waard et al. (1991, 1995); DNV-RP-O501 (2015); ISO 15156'),
        adv('economic', 'info', `The line needs ${fmt(need_)} mm of allowance; each millimetre of wall is ${fmt(steelPerMetre(p, 1e-3) * Lr)} t of steel over ${fmt(Lr)} km.`,
          dCA > 0 ? `Adding ${fmt(dCA)} mm of allowance costs about ${fmt(steelPerMetre(p, dCA / 1e3) * Lr)} t of line pipe; raising inhibitor availability from ${fmt(inp.avail)} to 0.99 is the operating-cost alternative to compare over the field life.` : `${fmt(inp.tCorr - need_)} mm of the allowance is not needed at this rate; the reserve can be traded against a longer inspection interval rather than removed.`,
          RBI_BASIS, `${steelCo2(steelPerMetre(p, 1e-3) * Lr).replace(' t of line-pipe steel', ' t of steel per millimetre of allowance')} Reliable inhibition extends asset life with far less material than extra wall, but adds a continuous chemical consumption.`),
      ],
      publish: { 'pipeline.corrosionRate': total / 1e3, 'pipeline.wallLife': Math.min(lifeWall, 1e6) },
      quantities: { aluminium_t: round(nAnode * inp.anodeMass / 1e3), pipe_km: inp.Lkm },
    };
  },
  verify: [{
    // de Waard & Milliams (1975): log10 V = 5.8 − 1710/T + 0.67 log10 p_CO2. T = 313.15 K, p_CO2 = 1 bar (2 mol % at 50 bar):
    // log10 V = 5.8 − 5.4606 = 0.3394 → V = 2.185 mm/year.
    // API RP 14E: v_e = 1.2199·C/√ρ = 1.2199·100/√100 = 12.20 m/s.
    name: 'de Waard–Milliams nomogram rate and API RP 14E erosional velocity',
    source: 'de Waard & Milliams, Corrosion 31(5) (1975); API RP 14E Eq. 2.14; hand calculation',
    inputs: { model: 'dwm', pIn: 50, pOut: 50, Tin: 40, Tout: 40, yCO2: 2, pH: 0, avail: 0, rhoMix: 100, erosC: 100 },
    expect: [{ kpi: 'crUninhibited', value: 2.185, tol: 0.003 }, { kpi: 'vErosional', value: 12.20, tol: 0.002 }],
  }],
};
