// Thermal rating and heat dissipation: IEC 60287 continuous current rating for the installation, temperatures of
// conductor, sheath, armour and surface at the load current, route sections and thermal hotspot, IEC 60853 cyclic
// and emergency ratings, and a 2-D finite-volume seabed temperature field with trench backfill as a cross-check.
import { rateCable, temperaturesAtCurrent, ratingEquation, externalThermalModel } from '../../engine/electrical/ampacity.js';
import { transientModel, cyclicRatingFactor, emergencyRating, responseTime, solveTrench2D } from '../../engine/electrical/thermaltransient.js';
import { cableInputs, buildElec, massesPerMetre, cableQuantities, embodiedCo2, adv, check, worst, sev, fmt, round, run, need, linspace, logspace, LIMIT_GENERIC, LIMIT_FE, TB880, CONFIRM, HOURS_YEAR, lossCo2Text } from './_common.js';

const DEFAULT_SECTIONS = [
  { name: 'Offshore, standard burial', length: 22000, type: 'buried', depth: 1.5, rho: 0.7, amb: 12 },
  { name: 'Sandwave field, deep burial', length: 1500, type: 'buried', depth: 3.5, rho: 0.7, amb: 12 },
  { name: 'Unburied at crossing (under rock)', length: 300, type: 'buried', depth: 0.8, rho: 1.0, amb: 12 },
  { name: 'Landfall HDD duct', length: 900, type: 'buried', depth: 8, rho: 1.0, amb: 14 },
  { name: 'Exposed on the seabed', length: 300, type: 'seabed', depth: 0, rho: 0, amb: 12 },
];
// Typical wind-farm daily cycle as fractions of the peak current (24 hourly values), used for the cyclic rating.
const DEFAULT_CYCLE = [0.55, 0.5, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0, 1.0, 0.95, 0.9, 0.85, 0.8, 0.8, 0.85, 0.9, 1.0, 1.0, 0.9, 0.8, 0.7, 0.65, 0.6];

export default {
  id: 's13-thermal-rating', suite: 13, title: 'Thermal rating and heat dissipation',
  summary: 'Continuous, cyclic and emergency current ratings of the cable in its installation, the temperatures of every layer at the load current, and the thermal hotspot along the route.',
  covers: [8], fidelity: 'numerical',
  standards: ['IEC 60287-1-1:2023 (current rating equations and losses)', 'IEC 60287-2-1:2023 (thermal resistances)', 'IEC 60853-2 (cyclic and emergency ratings)', 'CIGRE TB 880 (2022) Power cable rating examples (armour loss of three-core submarine cables)', 'CIGRE TB 640 (2015) A guide for rating calculations of insulated cables'],
  models: ['IEC 60287 current-rating method (thermal resistance network T1–T4)', 'Joule heating with temperature-dependent conductor resistivity (coupled electro-thermal fixed point)', 'Skin- and proximity-effect models', 'Sheath and armour loss factors', 'Dielectric heating', 'Fourier heat conduction: Kennelly image method and 2-D finite volumes with backfill regions', 'Two-zone soil drying-out', 'Transient heat diffusion: Van Wormer ladder network and exponential-integral soil response (IEC 60853)', 'Cyclic and emergency rating'],
  limits: `${LIMIT_GENERIC} ${TB880} The armour-loss reduction factor is therefore an input (1.0 = IEC). ${LIMIT_FE} The seabed field is steady 2-D conduction without pore-water convection (conservative in permeable sand); sections are rated independently without axial heat flow and with the loss factors of the base installation, and J-tube sections need the dedicated J-tube balance of the engine with project air-gap data. IEC 60287 / IEC 60853 constants (skin and proximity coefficients, covering resistivities, soil diffusivity, temperature limits) are overridable defaults to confirm against the governing edition.`,
  publishes: ['cable.ampacity', 'cable.conductorTemp', 'cable.lossPerM'],
  searchTerms: 'submarine cable current rating IEC 60287 ampacity seabed thermal resistivity',
  inputs: [
    ...cableInputs(),
    { key: 'U', label: 'System voltage (line-to-line)', unit: 'kV', type: 'number', default: 66, min: 1, max: 550, group: 'Electrical operating condition', role: 'bc' },
    { key: 'f', label: 'Frequency', unit: 'Hz', type: 'number', default: 50, min: 0, max: 60, group: 'Electrical operating condition', role: 'bc', help: '0 = DC (no skin, proximity, sheath, armour or dielectric loss).' },
    { key: 'I', label: 'Load current', unit: 'A', type: 'number', default: 580, min: 0, max: 5000, group: 'Electrical operating condition', role: 'bc', design: true, uncertain: { dist: 'normal', cov: 0.05 }, help: 'Initial energised state: steady current demand per phase at the hottest point.' },
    { key: 'armourLossFactor', label: 'Armour-loss reduction factor', unit: '-', type: 'number', default: 1.0, min: 0.05, max: 1, group: 'Electrical operating condition', role: 'criteria',
      help: `Multiplier on the IEC 60287 armour loss factor λ2. ${TB880} Measurements and 3-D FE typically support 0.2–0.6; 1.0 keeps the IEC value. Project decision to be justified.` },
    { key: 'thetaMax', label: 'Maximum conductor temperature', unit: '°C', type: 'number', default: 90, min: 40, max: 110, group: 'Electrical operating condition', role: 'criteria', help: `90 °C for XLPE and EPR (IEC 60502-2 / IEC 60840); 55 °C for mass-impregnated HVDC. ${CONFIRM}` },
    { key: 'inst', label: 'Installation', type: 'select', default: 'buried', group: 'Installation', role: 'bc', options: [{ value: 'buried', label: 'Buried in the seabed' }, { value: 'seabed', label: 'Lying on the seabed in water' }] },
    { key: 'depth', label: 'Burial depth to top of cable', unit: 'm', type: 'number', default: 1.5, min: 0.1, max: 30, group: 'Installation', role: 'ic', bind: 'cable.burialDepth', design: true, showIf: { key: 'inst', in: ['buried'] } },
    { key: 'rhoSoil', label: 'Soil thermal resistivity', unit: 'K·m/W', type: 'number', default: 0.7, min: 0.2, max: 4, group: 'Installation', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.2 }, showIf: { key: 'inst', in: ['buried'] }, help: 'Saturated marine sand 0.4–0.7, clay 0.7–1.0, organic mud above 1.0. Measure in situ along the route.' },
    { key: 'rhoBackfill', label: 'Trench backfill thermal resistivity', unit: 'K·m/W', type: 'number', default: 0.9, min: 0.2, max: 4, group: 'Installation', role: 'soil', showIf: { key: 'inst', in: ['buried'] }, help: 'Used in the 2-D field for the trench above and around the cable.' },
    { key: 'trenchW', label: 'Trench width', unit: 'm', type: 'number', default: 0.5, min: 0.1, max: 5, group: 'Installation', role: 'geometry', showIf: { key: 'inst', in: ['buried'] } },
    { key: 'ambient', label: 'Ambient (seabed / sea-water) temperature', unit: '°C', type: 'number', default: 12, min: -2, max: 35, group: 'Installation', role: 'bc', bind: 'site.Tbed', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Thermal boundary condition at the seabed surface and far field.' },
    { key: 'dry', label: 'Soil drying-out', type: 'bool', default: false, group: 'Installation', role: 'soil', showIf: { key: 'inst', in: ['buried'] }, help: 'Two-zone model for landfall sections above the water table; not relevant in saturated seabed.' },
    { key: 'rhoDry', label: 'Dry-soil thermal resistivity', unit: 'K·m/W', type: 'number', default: 2.5, min: 0.5, max: 6, group: 'Installation', role: 'soil', showIf: { key: 'dry', in: [true] } },
    { key: 'dThetaX', label: 'Critical temperature rise for drying', unit: 'K', type: 'number', default: 35, min: 5, max: 70, group: 'Installation', role: 'soil', showIf: { key: 'dry', in: [true] } },
    { key: 'hSea', label: 'Surface heat-transfer coefficient in water', unit: 'W/(m²·K)', type: 'number', default: 200, min: 5, max: 5000, group: 'Installation', role: 'bc', showIf: { key: 'inst', in: ['seabed'] }, help: 'Convective boundary to sea water; 100–1000 depending on current.' },
    { key: 'neighbours', label: 'Neighbouring heat sources', type: 'table', group: 'Installation', role: 'bc',
      columns: [{ key: 'name', label: 'Source', type: 'text' }, { key: 'x', label: 'Horizontal distance', unit: 'm', type: 'number' }, { key: 'depth', label: 'Depth', unit: 'm', type: 'number' }, { key: 'W', label: 'Heat', unit: 'W/m', type: 'number' }],
      default: [], help: 'Parallel circuits or crossing pipelines buried nearby (mutual heating by the image method).' },
    { key: 'sections', label: 'Route sections', type: 'table', group: 'Route', role: 'geometry',
      columns: [{ key: 'name', label: 'Section', type: 'text' }, { key: 'length', label: 'Length', unit: 'm', type: 'number' }, { key: 'type', label: 'Type (buried / seabed)', type: 'text' }, { key: 'depth', label: 'Depth to top', unit: 'm', type: 'number' }, { key: 'rho', label: 'Soil resistivity', unit: 'K·m/W', type: 'number' }, { key: 'amb', label: 'Ambient', unit: '°C', type: 'number' }],
      default: DEFAULT_SECTIONS.map((r) => ({ ...r })), help: 'Each section is rated separately; the lowest rating is the route rating (thermal hotspot). An empty table uses the single installation above.' },
    { key: 'cycle', label: 'Daily load cycle (fraction of peak)', type: 'series', default: { x: DEFAULT_CYCLE.map((_, i) => i), y: [...DEFAULT_CYCLE] }, xLabel: 'Hour', xUnit: 'h', yLabel: 'Current / peak', yUnit: '-', group: 'Cyclic and emergency rating', role: 'load',
      help: 'Equal steps over 24 h. A constant series gives a cyclic factor of 1.' },
    { key: 'diffusivity', label: 'Soil thermal diffusivity', unit: '10⁻⁶ m²/s', type: 'number', default: 0.5, min: 0.1, max: 2, group: 'Cyclic and emergency rating', role: 'soil', help: 'IEC 60853-2 suggests 0.5 × 10⁻⁶ m²/s for moist soil when not measured.' },
    { key: 'preload', label: 'Preload before the emergency', unit: '% of rating', type: 'number', default: 60, min: 0, max: 100, group: 'Cyclic and emergency rating', role: 'ic', help: 'Initial thermal state: steady current carried before the overload.' },
    { key: 'emDur', label: 'Emergency duration', unit: 'h', type: 'number', default: 6, min: 0.1, max: 500, group: 'Cyclic and emergency rating', role: 'operation' },
    { key: 'L', label: 'Route length for quantities', unit: 'km', type: 'number', default: 25, min: 0.05, max: 800, group: 'Route', role: 'geometry', bind: 'cable.length' },
    { key: 'cellsPerD', label: 'Finite-volume cells per cable diameter', type: 'int', default: 4, min: 3, max: 32, group: 'Solver', role: 'control', help: 'Resolution of the 2-D seabed temperature field around the cable.' },
  ],
  mesh: { key: 'cellsPerD', label: 'Cells per cable diameter', levels: [4, 8, 16], qoi: 'T4fv' },
  compute(inp, ctx) {
    const b = buildElec(inp, ctx), g = b.geom, De = g.De;
    const dc = !(inp.f > 0);
    const system = { ...(dc ? { dc: true } : { f: inp.f }), U: inp.U * 1e3, thetaMaxC: inp.thetaMax, armourLossReduction: inp.armourLossFactor };
    need(inp.thetaMax > inp.ambient + 5, 'The maximum conductor temperature must exceed the ambient temperature by at least 5 K.');
    const instOf = (type, depthTop, rho, amb) => {
      if (type === 'seabed') return { type: 'seabed', h: inp.hSea, ambientC: amb };
      need(depthTop > 0 && rho > 0, 'A buried section needs a positive depth and soil resistivity.');
      const ext = (inp.neighbours || []).filter((n) => n.W > 0).map((n) => { need(n.depth > 0 && Math.hypot(n.x, n.depth - depthTop) > De, `Neighbouring source '${n.name}' must be buried and clear of the cable.`); return { x: n.x, depth: n.depth, W: n.W }; });
      return { type: 'buried', depth: depthTop + De / 2, rhoSoil: rho, ambientC: amb, formation: 'single', ...(ext.length ? { externalSources: ext } : {}), ...(inp.dry ? { dryZone: { rhoDry: Math.max(inp.rhoDry, rho), criticalRiseK: inp.dThetaX } } : {}) };
    };
    const installation = instOf(inp.inst, inp.depth, inp.rhoSoil, inp.ambient);
    const rate = run(() => rateCable({ cable: g, installation, system }), 'Current rating');
    need(rate.I > 0, 'The dielectric loss and mutual heating alone use up the permissible temperature rise: no current can be carried in this installation.');
    let load;
    try { load = temperaturesAtCurrent({ cable: g, installation, system, I: inp.I }); } catch (e) { throw new Error(`No thermal equilibrium at ${inp.I} A (thermal runaway): the load current is far above the rating of ${fmt(rate.I)} A.`); }
    // Other sections, sweeps and the finite-volume cross-check reuse the converged AC resistance, loss factors and dielectric
    // loss of the base rating in the IEC 60287-1-1 §1.4 rating equation (their dependence on the installation is second order).
    const n = g.cores, al = g.alpha20, Rref = rate.R / (1 + al * (inp.thetaMax - 20));
    const fastRate = (inst) => {
      const ext = externalThermalModel(inst, g);
      return { T4: ext.T4, I: ratingEquation({ dTheta: inp.thetaMax - inst.ambientC - (ext.extRise || 0), R: rate.R, T1: g.T1, T2: g.T2, T3: g.T3, T4: ext.T4, n, Wd: rate.Wd, lambda1: rate.lambda1, lambda2: rate.lambda2, v: ext.v, dThetaX: ext.dThetaX }), ext };
    };
    // θ = θamb + I²·R(θ)·ΣT + Wd·Td with R linear in θ → closed form; a non-positive denominator is thermal runaway
    const fastTemp = (inst, I, T4) => {
      const k = I * I * Rref * (g.T1 + n * (1 + rate.lambda1) * g.T2 + n * (1 + rate.lambda1 + rate.lambda2) * (g.T3 + T4));
      const den = 1 - k * al;
      return den > 0 ? (inst.ambientC + rate.Wd * (0.5 * g.T1 + n * (g.T2 + g.T3 + T4)) + k * (1 - 20 * al)) / den : Infinity;
    };
    const secRows = (inp.sections || []).filter((s) => s.length > 0);
    let route = null;
    if (secRows.length) {
      let kp = 0;
      const secs = secRows.map((s, i) => {
        const type = String(s.type || 'buried').trim().toLowerCase();
        need(type === 'buried' || type === 'seabed', `Route section ${i + 1}: type must be 'buried' or 'seabed'.`);
        const inst = instOf(type, s.depth, s.rho, s.amb);
        need(inp.thetaMax > inst.ambientC, `Route section ${i + 1}: ambient temperature must be below the maximum conductor temperature.`);
        const r = run(() => fastRate(inst), `Route section ${i + 1}`);
        const row = { name: s.name || `Section ${i + 1}`, length: s.length, kpStart: kp, I: r.I, conductorAtLoadC: fastTemp(inst, inp.I, r.T4) };
        kp += s.length; return row;
      });
      const lim = secs.reduce((a, q) => (q.I < a.I ? q : a), secs[0]);
      route = { sections: secs, I: lim.I, limitingName: lim.name };
    }
    const Iroute = route ? Math.min(route.I, rate.I) : rate.I;
    // transient model (IEC 60853): cyclic and emergency ratings
    const buried = inp.inst === 'buried';
    const model = run(() => transientModel({ geometry: g, lambda1: rate.lambda1, lambda2: rate.lambda2, external: buried ? { type: 'buried', rhoSoil: inp.rhoSoil, diffusivity: inp.diffusivity * 1e-6, depth: inp.depth + De / 2 } : { type: 'resistance', T4: rate.T4eff } }), 'Transient thermal model');
    const cyc = (inp.cycle?.y || []).map(Number);
    need(cyc.length >= 2 && cyc.every((v) => v >= 0) && Math.max(...cyc) > 0, 'The load cycle needs at least two non-negative values with a positive peak.');
    const peak = Math.max(...cyc), cycN = cyc.map((v) => v / peak);
    const cy = run(() => cyclicRatingFactor({ model, cycle: cycN, dt: 86400 / cycN.length, nCycles: 10 }), 'Cyclic rating factor');
    const em = run(() => emergencyRating({ model, R: rate.R, preloadCurrent: inp.preload / 100 * rate.I, duration: inp.emDur * 3600, thetaMaxC: inp.thetaMax, ambientC: inp.ambient, Wd: rate.Wd }), 'Emergency rating');
    const tau = responseTime({ model });
    // 2-D finite-volume field: native soil with a backfilled trench; cable heat = total loss at the load current
    let fv = null, T4fv = rate.T4eff, Tfield = null;
    if (buried) {
      const zc = inp.depth + De / 2, Wl = Math.max(load.losses.total, 1);
      fv = run(() => solveTrench2D({ cables: [{ x: 0, depth: zc, diameter: De, W: Wl }], rhoSoil: inp.rhoSoil, ambientC: inp.ambient, h: De / inp.cellsPerD, tol: 1e-8, ratio: 1.4,
        regions: [{ x0: -inp.trenchW / 2, x1: inp.trenchW / 2, z0: 0, z1: zc + De / 2 + 0.05, rhoT: inp.rhoBackfill }] }), 'Finite-volume seabed field');
      T4fv = fv.cables[0].T4;
      // decimated field around the trench for the heat map
      const xe = fv.xEdges, ze = fv.zEdges, ix = [], iz = [];
      for (let i = 0; i < fv.nx; i++) { const xm = 0.5 * (xe[i] + xe[i + 1]); if (Math.abs(xm) <= Math.max(3 * zc, 2)) ix.push(i); }
      for (let k = 0; k < fv.nz; k++) { const zm = 0.5 * (ze[k] + ze[k + 1]); if (zm <= Math.max(3 * zc, 2)) iz.push(k); }
      const pick = (a, n) => (a.length <= n ? a : Array.from({ length: n }, (_, i) => a[Math.round(i * (a.length - 1) / (n - 1))]));
      const px = pick(ix, 60), pz = pick(iz, 40);
      Tfield = { xs: px.map((i) => round(0.5 * (xe[i] + xe[i + 1]), 4)), ys: pz.map((k) => round(-0.5 * (ze[k] + ze[k + 1]), 4)), z: pz.map((k) => px.map((i) => round(fv.T[k * fv.nx + i], 4))) };
    }
    const rateFv = buried ? run(() => fastRate({ type: 'custom', T4: T4fv, ambientC: inp.ambient }), 'Rating with the finite-volume T4') : rate;
    const Igov = Math.min(Iroute, rateFv.I);
    const checks = [
      check('Load current against the continuous rating (this installation)', inp.I, rate.I, 'A', 'IEC 60287-1-1 §1.4.1 (permissible current rating)', 'Demand is the load current, capacity the current that brings the conductor to its maximum temperature here. Above it the insulation ages faster than its design basis (life roughly halves per 8–10 K).'),
      check('Conductor temperature at the load current', load.conductorC - inp.ambient, inp.thetaMax - inp.ambient, 'K', 'IEC 60287-1-1 §1.4; IEC 60502-2 / IEC 60840 (maximum conductor temperature)', 'Conductor temperature rise above ambient at the load current against the permissible rise.'),
      ...(route ? [check('Load current against the route rating (hotspot section)', inp.I, route.I, 'A', 'IEC 60287 applied per route section; CIGRE TB 640', `Lowest section rating along the route — here '${route.limitingName}'. A short hot section limits the whole circuit unless it is mitigated locally.`)] : []),
      ...(buried ? [check('Load current against the rating with the 2-D backfill field', inp.I, rateFv.I, 'A', 'IEC 60287-1-1 §1.4 with T4 from 2-D finite volumes (IEC TR 62095 approach)', 'Rating recomputed with the external thermal resistance from the 2-D field that includes the trench backfill. A poorer backfill than the native soil lowers the rating.')] : []),
    ];
    const gov = worst(checks);
    const depths = linspace(0.3, Math.max(4, inp.depth * 1.5), 14), rhos = linspace(0.4, 1.6, 13);
    const sd = buried ? { I: depths.map((d) => fastRate({ ...installation, depth: d + De / 2 }).I) } : null;
    const sr = buried ? { I: rhos.map((r) => fastRate({ ...installation, rhoSoil: r, ...(inp.dry ? { dryZone: { rhoDry: Math.max(inp.rhoDry, r), criticalRiseK: inp.dThetaX } } : {}) }).I) } : null;
    const tt = logspace(600, 3.15e8, 50);
    const stepC = tt.map((t) => inp.ambient + rate.I ** 2 * rate.R * model.step(t) + rate.Wd * model.dielectricSteady);
    const mass = massesPerMetre(b), Lm = inp.L * 1e3;
    const lossKw = load.losses.total * Lm / 1e3, lossMwh = lossKw * HOURS_YEAR / 1e3 * cy.lossLoadFactor;
    const nodes = [{ n: 'Conductor', T: load.conductorC }, { n: 'Metallic sheath', T: load.sheathC }, { n: 'Armour', T: load.armourC }, { n: 'Cable surface', T: load.surfaceC }, { n: 'Ambient', T: inp.ambient }];
    const lam2Iec = inp.armourLossFactor > 0 ? rate.lambda2 / inp.armourLossFactor : rate.lambda2;
    return {
      kpis: [
        { key: 'ampacity', label: 'Continuous current rating', value: round(rate.I, 5), unit: 'A', explain: 'Current that brings the conductor to its maximum temperature in this installation with 100 % load factor (IEC 60287). It is the thermal ampacity of the section.' },
        { key: 'ampacityRoute', label: 'Route rating (governing section)', value: round(Igov), unit: 'A', explain: 'Lowest of the section ratings and of the rating with the 2-D backfill field: the continuous rating of the circuit.' },
        { key: 'conductorTemp', label: 'Conductor temperature at the load current', value: round(load.conductorC), unit: '°C', limit: inp.thetaMax, sense: 'max', explain: 'Steady conductor temperature at the entered load current with resistances at their own temperatures.' },
        { key: 'sheathTemp', label: 'Sheath temperature', value: round(load.sheathC), unit: '°C', explain: 'Metallic sheath temperature at the load current; it sets the sheath resistance and the thermal ageing of the core jacket.' },
        { key: 'surfaceTemp', label: 'Cable surface temperature', value: round(load.surfaceC), unit: '°C', explain: 'Temperature of the serving against the soil or water. In unsaturated soil a high value can dry the soil out; in sediment it is the value regulators ask about (e.g. the 2 K criterion).' },
        { key: 'utilisation', label: 'Thermal rating utilisation', value: round(inp.I / Igov), unit: '-', limit: 1, sense: 'max', explain: 'Load current divided by the governing continuous rating.' },
        { key: 'lossTotal', label: 'Total losses at the load current', value: round(load.losses.total), unit: 'W/m', explain: 'Conductor, dielectric, sheath and armour losses per metre of cable: the heat the seabed must carry away.' },
        { key: 'Rac', label: 'Conductor AC resistance at operating temperature', value: round(load.R * 1e6, 5), unit: 'µΩ/m', explain: 'DC resistance at the conductor temperature increased by skin and proximity effects.' },
        { key: 'lambda1', label: 'Sheath loss factor λ1', value: round(rate.lambda1), unit: '-', explain: 'Sheath losses (circulating and eddy currents) as a fraction of conductor losses.' },
        { key: 'lambda2', label: 'Armour loss factor λ2 (as used)', value: round(rate.lambda2), unit: '-', explain: `Armour losses as a fraction of conductor losses after the reduction factor (${inp.armourLossFactor}). The IEC value is ${fmt(lam2Iec)} and is known to be high for three-core submarine cables.` },
        { key: 'T1', label: 'T1 (insulation)', value: round(g.T1), unit: 'K·m/W', explain: 'Thermal resistance between one conductor and its sheath.' },
        { key: 'T4', label: 'T4 (external, IEC image method)', value: round(rate.T4eff), unit: 'K·m/W', explain: 'External thermal resistance of the surroundings for this installation; for a buried cable it grows with the logarithm of depth and linearly with soil resistivity.' },
        { key: 'T4fv', label: 'T4 from the 2-D field with backfill', value: round(T4fv, 5), unit: 'K·m/W', explain: 'External thermal resistance from the finite-volume temperature field including the trench backfill; compare with the IEC value for uniform soil.' },
        { key: 'cyclicFactor', label: 'Cyclic rating factor M', value: round(cy.M), unit: '-', explain: 'Permissible peak of the daily cycle divided by the continuous rating (IEC 60853): thermal inertia lets the peak exceed the steady rating.' },
        { key: 'cyclicRating', label: 'Permissible cyclic peak current', value: round(cy.M * rate.I), unit: 'A', explain: 'M × continuous rating for the entered load cycle.' },
        { key: 'emergencyRating', label: 'Emergency rating', value: round(em.I), unit: 'A', explain: `Current that can be carried for ${inp.emDur} h from a ${inp.preload} % preload without exceeding the maximum conductor temperature.` },
        { key: 'tau', label: 'Thermal response time (63 %)', value: round(tau / 86400), unit: 'days', explain: 'Time for the conductor to reach 63 % of its final rise after a load step; buried cables take weeks because the soil stores heat.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Temperature by layer at the load current', categories: nodes.map((n) => n.n), series: [{ name: 'Temperature (°C)', y: nodes.map((n) => round(n.T)) }],
          explain: 'Temperature from conductor to ambient. The largest step shows which thermal resistance controls the rating: the insulation (T1) or the surroundings (T4).' },
        { type: 'bar', title: 'Loss breakdown at the load current', categories: ['Conductor', 'Dielectric', 'Sheath', 'Armour'], series: [{ name: 'Loss (W/m)', y: [load.losses.conductor, load.losses.dielectric, load.losses.sheath, load.losses.armour].map((v) => round(v)) }],
          explain: 'Where the heat is generated. A large armour share is the IEC 60287 over-estimate for three-core cables — the reduction factor acts on that bar only.' },
        ...(sd ? [{ type: 'line', title: 'Rating against burial depth', x: { label: 'Depth to top of cable', unit: 'm' }, y: { label: 'Continuous rating', unit: 'A' }, series: [{ name: 'Rating', x: depths, y: sd.I }, { name: 'Load current', x: [depths[0], depths[13]], y: [inp.I, inp.I], dash: true }],
          explain: 'Deeper burial derates the cable. The depth at which the curve crosses the load current is the deepest acceptable burial for this soil.' },
        { type: 'line', title: 'Rating against soil thermal resistivity', x: { label: 'Soil thermal resistivity', unit: 'K·m/W' }, y: { label: 'Continuous rating', unit: 'A' }, series: [{ name: 'Rating', x: rhos, y: sr.I }, { name: 'Load current', x: [rhos[0], rhos[12]], y: [inp.I, inp.I], dash: true }],
          explain: 'Sensitivity to the least certain input. Read off the resistivity at which the rating falls to the load current and compare with the survey range.' }] : []),
        { type: 'line', title: 'Conductor temperature after a step to rated current', x: { label: 'log₁₀ time', unit: 's' }, y: { label: 'Conductor temperature', unit: '°C' }, series: [{ name: 'Step response', x: tt.map((t) => Math.log10(t)), y: stepC }, { name: 'Limit', x: [Math.log10(tt[0]), Math.log10(tt[49])], y: [inp.thetaMax, inp.thetaMax], dash: true }],
          explain: 'Transient heating from ambient at rated current: the cable itself responds within hours, the soil over weeks to years. The gap to the limit at a given time is the emergency capacity.' },
        ...(route ? [{ type: 'line', title: 'Conductor temperature along the route at the load current', x: { label: 'Chainage', unit: 'km' }, y: { label: 'Conductor temperature', unit: '°C' },
          series: [{ name: 'Conductor', x: route.sections.flatMap((s) => [s.kpStart / 1e3, (s.kpStart + s.length) / 1e3]), y: route.sections.flatMap((s) => { const v = Math.min(s.conductorAtLoadC ?? 0, 300); return [v, v]; }) }, { name: 'Limit', x: [0, route.sections.reduce((a, s) => a + s.length, 0) / 1e3], y: [inp.thetaMax, inp.thetaMax], dash: true }],
          explain: 'Conductor temperature section by section (no axial heat flow). The highest plateau is the thermal hotspot that limits the circuit; values capped at 300 °C mean thermal runaway.' }] : []),
        ...(Tfield ? [{ type: 'heatmap', title: 'Seabed temperature field around the cable', x: { label: 'Horizontal distance', unit: 'm' }, y: { label: 'Elevation relative to seabed', unit: 'm' }, xs: Tfield.xs, ys: Tfield.ys, z: Tfield.z, zLabel: 'Temperature', zUnit: '°C',
          explain: 'Steady temperature of the sediment at the load current, with the backfilled trench. Look at how far the warm zone reaches the seabed surface and whether the backfill concentrates the gradient.' }] : []),
      ],
      tables: [
        { title: 'Thermal network', explain: 'Thermal resistances and loss factors of the IEC 60287 network used in the rating equation, for checking against the manufacturer calculation.',
          columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [{ q: 'T1 insulation', v: round(g.T1), u: 'K·m/W' }, { q: 'T2 bedding / fillers', v: round(g.T2), u: 'K·m/W' }, { q: 'T3 serving', v: round(g.T3), u: 'K·m/W' }, { q: 'T4 external', v: round(rate.T4eff), u: 'K·m/W' },
            { q: 'Conductor R at θmax', v: round(rate.R * 1e6), u: 'µΩ/m' }, { q: 'Skin effect ys', v: round(rate.ys), u: '-' }, { q: 'Proximity effect yp', v: round(rate.yp), u: '-' }, { q: 'λ1 sheath', v: round(rate.lambda1), u: '-' }, { q: 'λ2 armour (used)', v: round(rate.lambda2), u: '-' },
            { q: 'Dielectric loss per core', v: round(rate.Wd), u: 'W/m' }, { q: 'Overall diameter', v: round(De * 1e3), u: 'mm' }, { q: 'Capacitance', v: round(g.C * 1e9), u: 'µF/km' }] },
        ...(route ? [{ title: 'Route sections', explain: 'One row per section: its continuous rating and the conductor temperature at the load current. The lowest rating is the route rating.',
          columns: [{ key: 'name', label: 'Section' }, { key: 'kp', label: 'Start', unit: 'km' }, { key: 'len', label: 'Length', unit: 'm' }, { key: 'I', label: 'Rating', unit: 'A' }, { key: 'T', label: 'Conductor at load', unit: '°C' }],
          rows: route.sections.map((s) => ({ name: s.name, kp: round(s.kpStart / 1e3), len: s.length, I: round(s.I), T: Number.isFinite(s.conductorAtLoadC) ? round(s.conductorAtLoadC) : 'runaway' })) }] : []),
      ],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        route ? `Thermal hotspot: '${route.limitingName}' at ${fmt(route.I)} A.` : 'No route sections entered: single installation rated.',
        `Armour loss factor: IEC ${fmt(lam2Iec)}, used ${fmt(rate.lambda2)} (reduction factor ${inp.armourLossFactor}).${rate.dried ? ' The soil dries out at the rating (two-zone model active).' : ''}`,
        buried ? `2-D field heat balance error ${fmt(fv.heatBalanceError * 100)} %; T4 ${fmt(T4fv)} against ${fmt(rate.T4eff)} K·m/W from the image method.` : 'Cable in water: external resistance from the surface coefficient only.',
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name}: ${fmt(gov.uc)}. Continuous rating ${fmt(rate.I)} A here, ${fmt(Igov)} A for the route; conductor ${fmt(load.conductorC)} °C at ${inp.I} A; cyclic peak ${fmt(cy.M * rate.I)} A; ${inp.emDur} h emergency ${fmt(em.I)} A.`,
          gov.uc > 1 ? `Mitigate the limiting section: shallower burial or thermally stable backfill there, a larger conductor on that length (${fmt(inp.area * gov.uc ** 2)} mm² scales losses back to the limit), or — where justified by measurement or FE — an armour-loss reduction factor below ${inp.armourLossFactor}.`
            : 'Confirm soil thermal resistivity by in-situ measurement along the route, document the armour-loss basis (IEC value or a justified reduction per CIGRE TB 880), and consider distributed temperature sensing on the fibre to verify the hotspot in service.',
          'IEC 60287-1-1:2023 / -2-1:2023; IEC 60853-2; CIGRE TB 880 (2022); CIGRE TB 640 (2015)'),
        adv('economic', 'info', `Losses at the load current are ${fmt(load.losses.total)} W/m — ${fmt(lossKw)} kW over ${inp.L} km, about ${fmt(lossMwh)} MWh a year with the loss-load factor of the cycle (${fmt(cy.lossLoadFactor)}).`,
          `Rating margin is ${fmt((1 - inp.I / Igov) * 100)} %. ${inp.I / Igov < 0.8 ? 'The conductor could be one size smaller on thermal grounds; check that the capitalised value of the extra losses does not exceed the metal saved (IEC 60287-3-2).' : 'The section is thermally well used; dynamic rating with fibre temperature monitoring (cyclic factor ' + fmt(cy.M) + ') can release a further ' + fmt((cy.M - 1) * 100) + ' % at peak without more copper.'}`,
          'IEC 60287-3-2 (economic optimisation of cable size); whole-life costing ISO 15686-5',
          `${lossCo2Text(lossMwh)} Cable metal over ${inp.L} km: ${fmt(mass.copper * Lm / 1e3)} t copper, ${fmt(mass.lead * Lm / 1e3)} t lead, ${fmt(mass.steel * Lm / 1e3)} t steel ≈ ${fmt(embodiedCo2(mass, Lm))} t CO₂e embodied (screening factors).`),
      ],
      publish: { 'cable.ampacity': Igov, 'cable.conductorTemp': load.conductorC, 'cable.lossPerM': load.losses.total },
      quantities: cableQuantities(mass, Lm),
    };
  },
  verify: [
    // Single-core 630 mm² Cu / XLPE, no sheath, no armour, buried 1.0 m to top, ρsoil 0.7, ambient 12 °C, 66 kV 50 Hz.
    // R′(90 °C) = 0.0283e-3 × (1 + 0.00393 × 70) = 3.6085e-5 Ω/m; xs² = 8πf/R′·1e-7 = 3.4823 → ys = xs⁴/(192 + 0.8·xs⁴) = 0.06012;
    // isolated cable → yp = 0 → R = 3.8254e-5 Ω/m. dc = 29.854 mm; T1 = 3.5/(2π)·ln(1 + 22/29.854) = 0.30755;
    // coverings in series (T3): PE jacket 2 mm on 51.854 mm 0.04139 + bedding 2 mm (ρ 6) 0.06605 + serving 4 mm 0.06989 = 0.17733;
    // De = 67.854 mm, L = 1.03393 m, u = 30.475 → T4 = 0.7/(2π)·ln(u + √(u² − 1)) = 0.45787 K·m/W.
    // C = 2.5e-9/(18·ln(49.854/31.854)) = 3.1006e-10 F/m; Wd = ωC·U0²·tanδ = 314.16 × 3.1006e-10 × 38105² × 0.001 = 0.1414 W/m.
    // I = sqrt((78 − 0.1414·(0.5·0.30755 + 0.17733 + 0.45787)) / (3.8254e-5·(0.30755 + 0.17733 + 0.45787))) = 1469.6 A.
    { name: 'IEC 60287 rating equation, single isolated buried cable (hand calculation of R, ys, T1, T3, T4, Wd)',
      source: 'IEC 60287-1-1 §1.4.1.1, §2.1.2; IEC 60287-2-1 §4.1.2.1, §4.2.2 (hand calculation)',
      inputs: { cores: '1', tSheath: 0, dWire: 0, depth: 1.0, sections: [], cellsPerD: 4 },
      expect: [{ kpi: 'ampacity', value: 1469.6, tol: 0.003 }, { kpi: 'T4', value: 0.45787, tol: 0.002 }, { kpi: 'T1', value: 0.30755, tol: 0.002 }] },
    // Uniform soil (backfill = native): the finite-volume T4 must reproduce Kennelly's image solution
    // T4 = ρ/(2π)·ln(u + √(u² − 1)) for the default three-core cable: De = 150.1 mm, L = 1.575 m, u = 20.99 → 0.4163 K·m/W.
    { name: '2-D finite-volume field in uniform soil reproduces the Kennelly image solution',
      source: 'IEC 60287-2-1 §4.2.2 (single isolated buried cable); Kennelly (1893)',
      inputs: { rhoBackfill: 0.7, sections: [], cellsPerD: 8 },
      expect: [{ kpi: 'T4fv', value: 0.4163, tol: 0.03 }] },
  ],
};
