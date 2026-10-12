// Rubble-mound breakwater and revetment armour: required rock or concrete-unit size, damage of the provided armour,
// toe berm, layer build-up, filter rules and section quantities, with the sea-level-rise pathway.
import { relativeDensity, HUDSON_KD, hudson, vanDerMeerRock, vanDerMeerDamage, vanDerMeerShallow, concreteUnitVanDerMeer, stabilityNumberSize,
  toeStabilityVanDerMeer, massToDn, armourLayer, filterRules, underlayerMass, rubbleMoundQuantities } from '../../engine/coastal/armour.js';
import { criticalCircle } from '../../engine/geotech/slope.js';
import { CO2_ROCK, CO2_CONC, CO2_MGO, MGO_T_PER_DAY, EDITION_NOTE, need, run, finite, r, f, linspace, check, eng, eco, slrInputs, slrRows, SLR_BASIS } from './_common.js';

const KD_CLASSES = Object.keys(HUDSON_KD);

export default {
  id: 's07-armour-stability', suite: 7, title: 'Breakwater armour stability',
  summary: 'Sizes the armour, underlayer and toe of a rubble mound, predicts the damage of the armour provided and takes off the section quantities.',
  covers: [1], fidelity: 'analytical',
  standards: ['CIRIA C683 Rock Manual (2007) §5.2.2', 'USACE Shore Protection Manual (1984) ch. 7', 'Van der Meer (1988) rock slopes and gravel beaches under wave attack', 'EN 13383-1 armourstone gradings'],
  models: ['Hudson armour stability', 'Van der Meer rock stability (plunging / surging)', 'Van Gent shallow-water modification', 'Van der Meer cubes and tetrapods', 'Van der Meer toe stability', 'Terzaghi filter rules', 'Bishop simplified slope stability (limit equilibrium)', 'Rubble-mound section geometry'],
  limits: `Trunk sections under head-on irregular waves; roundheads, oblique attack, overtopping damage to the rear slope, settlement, seismic and wave-induced pore-pressure effects on the mound are not covered; the slip check is static with buoyant weights (Bishop circles — Morgenstern–Price non-circular surfaces are outside the module). No single-layer concrete-unit stability numbers are built in — the design Ns is an input taken from the licensor's guidance. ${EDITION_NOTE}`,
  publishes: ['marine.armourM50', 'marine.armourDn50', 'quantities.rock_t'],
  searchTerms: 'rubble mound breakwater armour stability Van der Meer Hudson damage',
  inputs: [
    { key: 'Hs', label: 'Significant wave height at the toe', unit: 'm', type: 'number', default: 4.2, min: 0.3, max: 15, group: 'Wave boundary', role: 'bc', bind: 'marine.HsToe', uncertain: { dist: 'lognormal', cov: 0.1 }, help: 'Before the sea-level-rise allowance below.' },
    { key: 'Tm', label: 'Mean wave period Tm', unit: 's', type: 'number', default: 9, min: 2, max: 22, group: 'Wave boundary', role: 'bc', uncertain: { dist: 'normal', cov: 0.08 }, help: '≈ Tp/1.2 for a JONSWAP sea.' },
    { key: 'Nw', label: 'Number of waves in the design storm', unit: '-', type: 'number', default: 3000, min: 500, max: 7500, group: 'Wave boundary', role: 'bc', help: 'Storm duration / Tm; the formula saturates near 7500.' },
    { key: 'H2r', label: 'H2% / Hs at the toe', unit: '-', type: 'number', default: 1.3, min: 1.0, max: 1.5, group: 'Wave boundary', role: 'bc', help: '1.4 for Rayleigh-distributed deep-water waves, falling to 1.1–1.2 in a saturated surf zone.' },
    { key: 'hToe', label: 'Water depth at the toe of the structure', unit: 'm', type: 'number', default: 9, min: 1, max: 60, group: 'Water level', role: 'bc', bind: 'marine.hToe' },
    { key: 'depthLimited', label: 'Toe wave is depth-limited', type: 'bool', default: true, group: 'Water level', role: 'bc', help: 'Scales Hs with toe depth as the sea rises.' },
    ...slrInputs(),
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Water level', role: 'environment', bind: 'site.rhoW' },
    { key: 'unit', label: 'Armour type', type: 'select', default: 'rock', group: 'Armour', role: 'material',
      options: [{ value: 'rock', label: 'Quarry rock, double layer' }, { value: 'cube', label: 'Concrete cubes, double layer' }, { value: 'tetrapod', label: 'Tetrapods, double layer' }, { value: 'single', label: 'Single-layer interlocking unit (Ns from licensor)' }] },
    { key: 'rock', label: 'Rock source', type: 'select', source: 'materials:rock', default: 'rock-granite', group: 'Armour', role: 'material', help: 'Gives the armour density; concrete units use the concrete-unit record.' },
    { key: 'cotA', label: 'Front slope cot α', unit: '-', type: 'number', default: 2, min: 1.2, max: 6, group: 'Section', role: 'geometry', design: true },
    { key: 'M50', label: 'Armour mass provided (M50 or unit mass)', unit: 't', type: 'number', default: 14, min: 0.02, max: 120, group: 'Armour', role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.1 }, help: 'As-placed armour: median rock mass of the grading, or the concrete unit mass.' },
    { key: 'P', label: 'Notional permeability P', unit: '-', type: 'number', default: 0.4, min: 0.1, max: 0.6, group: 'Armour', role: 'material', help: '0.1 impermeable core, 0.4 armour–filter–core, 0.5 armour on permeable core, 0.6 homogeneous mound.' },
    { key: 'Sd', label: 'Tolerable damage level Sd', unit: '-', type: 'number', default: 2, min: 1, max: 15, group: 'Criteria', role: 'criteria', help: 'Start of damage 2 (cot α 1.5–2) to 3 (cot α 4–6); failure 8–17. Rock Manual Table 5.23 — confirm.' },
    { key: 'Nod', label: 'Tolerable displaced units Nod (concrete, toe)', unit: '-', type: 'number', default: 0.5, min: 0.05, max: 4, group: 'Criteria', role: 'criteria', help: 'Units displaced per strip one Dn wide: 0–0.5 start of damage, 2 severe (cubes).' },
    { key: 'NsUnit', label: 'Design stability number of the single-layer unit', unit: '-', type: 'number', default: 2.7, min: 1, max: 4, group: 'Criteria', role: 'criteria', showIf: { key: 'unit', in: ['single'] }, help: 'Hs/(Δ Dn) from the unit licensor\'s current design guidance; no value is built in.' },
    { key: 'kdClass', label: 'Hudson KD class', type: 'select', default: 'quarrystone-rough-angular', group: 'Criteria', role: 'criteria', options: KD_CLASSES.map((k) => ({ value: k, label: k })), help: 'No-damage trunk values recalled from SPM (1984) Table 7-8.' },
    { key: 'KD', label: 'Hudson KD override (0 = class value)', unit: '-', type: 'number', default: 0, min: 0, max: 40, group: 'Criteria', role: 'criteria' },
    { key: 'breaking', label: 'Breaking waves on the structure (Hudson)', type: 'bool', default: true, group: 'Criteria', role: 'criteria' },
    { key: 'hudsonH', label: 'Hudson design wave / Hs', unit: '-', type: 'number', default: 1.0, min: 1, max: 1.8, group: 'Criteria', role: 'criteria', help: '1.0 applies KD to Hs (Rock Manual usage); SPM (1984) applies it to H1/10 = 1.27 Hs, which doubles the mass.' },
    { key: 'ht', label: 'Water depth above the toe berm', unit: 'm', type: 'number', default: 6.5, min: 0.5, max: 55, group: 'Toe and layers', role: 'geometry', design: true },
    { key: 'nLayers', label: 'Armour layers', type: 'int', default: 2, min: 1, max: 3, group: 'Toe and layers', role: 'geometry' },
    { key: 'kt', label: 'Layer coefficient kt', unit: '-', type: 'number', default: 1.0, min: 0.7, max: 1.6, group: 'Toe and layers', role: 'material', help: 'Rough quarrystone ≈ 1.0, cubes 1.1, tetrapods 1.04.' },
    { key: 'nv', label: 'Armour-layer porosity', unit: '-', type: 'number', default: 0.37, min: 0.25, max: 0.65, group: 'Toe and layers', role: 'material' },
    { key: 'D15f', label: 'Underlayer D15', unit: 'm', type: 'number', default: 0.35, min: 0.005, max: 2, group: 'Toe and layers', role: 'material' },
    { key: 'D85b', label: 'Core D85', unit: 'm', type: 'number', default: 0.25, min: 0.001, max: 2, group: 'Toe and layers', role: 'material' },
    { key: 'D15b', label: 'Core D15', unit: 'm', type: 'number', default: 0.03, min: 0.0001, max: 1, group: 'Toe and layers', role: 'material' },
    { key: 'crestW', label: 'Crest width', unit: 'm', type: 'number', default: 5, min: 1, max: 40, group: 'Section', role: 'geometry' },
    { key: 'height', label: 'Mound height (bed to crest)', unit: 'm', type: 'number', default: 15, min: 2, max: 60, group: 'Section', role: 'ic', help: 'As-built section used for the quantities.' },
    { key: 'cotR', label: 'Rear slope cot α', unit: '-', type: 'number', default: 1.5, min: 1.2, max: 5, group: 'Section', role: 'geometry' },
    { key: 'length', label: 'Breakwater length', unit: 'm', type: 'number', default: 600, min: 1, max: 20000, group: 'Section', role: 'geometry' },
    { key: 'slip', label: 'Check deep-seated slip of the seaward slope (Bishop)', type: 'bool', default: true, group: 'Global stability', role: 'control' },
    { key: 'phiMound', label: 'Friction angle of the rock mound', unit: 'deg', type: 'number', default: 40, min: 30, max: 50, group: 'Global stability', role: 'material' },
    { key: 'cFound', label: 'Foundation soil cohesion / undrained strength', unit: 'kPa', type: 'number', default: 5, min: 0, max: 500, group: 'Global stability', role: 'soil', bind: 'site.su0' },
    { key: 'phiFound', label: 'Foundation soil friction angle', unit: 'deg', type: 'number', default: 30, min: 0, max: 45, group: 'Global stability', role: 'soil', bind: 'site.phiDeg' },
    { key: 'gamFound', label: 'Foundation soil submerged unit weight', unit: 'kN/m³', type: 'number', default: 9.5, min: 4, max: 13, group: 'Global stability', role: 'soil', bind: 'site.gammaSub' },
    { key: 'sfSlip', label: 'Required factor against slip', unit: '-', type: 'number', default: 1.3, min: 1, max: 2, group: 'Criteria', role: 'criteria' },
    { key: 'placeRate', label: 'Rock placing rate of the marine spread', unit: 't/day', type: 'number', default: 3000, min: 100, max: 30000, group: 'Construction', role: 'operation' },
  ],
  compute(inp, ctx) {
    const concrete = inp.unit !== 'rock';
    const mat = ctx.material(concrete ? 'concrete-unit' : inp.rock) || ctx.material(inp.rock);
    need(mat && mat.rho > 0, 'The selected armour material has no density in the materials library.');
    const rhoS = mat.rho, rockRho = (ctx.material(inp.rock) || mat).rho;
    const delta = relativeDensity({ rhoS, rhoW: inp.rhoW });
    need(delta > 0.5, 'Armour density must exceed about 1.5 × water density to be stable.');
    need(inp.ht < inp.hToe + inp.slr, 'The toe berm must lie below the water surface and above the bed: depth above the toe berm must be less than the toe depth.');
    const M50prov = inp.M50 * 1000;
    const DnProv = Math.cbrt(M50prov / rhoS);

    /** Required unit for a given wave height and toe depth. */
    const required = (Hs, h) => {
      const shallow = h / Hs < 3;
      if (inp.unit === 'rock') {
        const o = shallow ? vanDerMeerShallow({ Hs, H2: inp.H2r * Hs, Tm10: inp.Tm * 1.1, cotAlpha: inp.cotA, P: inp.P, S: inp.Sd, N: inp.Nw, delta, rhoS })
          : vanDerMeerRock({ Hs, Tm: inp.Tm, cotAlpha: inp.cotA, P: inp.P, S: inp.Sd, N: inp.Nw, delta, rhoS });
        return { Ns: o.Ns, Dn: o.Dn50, M: o.M50, regime: `${o.regime}${shallow ? ', shallow-water form' : ''}`, xi: o.xi };
      }
      if (inp.unit === 'single') { const o = stabilityNumberSize({ Hs, Ns: inp.NsUnit, delta, rhoS }); return { Ns: inp.NsUnit, Dn: o.Dn, M: o.M, regime: 'licensor stability number', xi: NaN }; }
      const o = concreteUnitVanDerMeer({ unit: inp.unit, Hs, Tm: inp.Tm, Nod: inp.Nod, N: inp.Nw, delta, rhoS });
      return { Ns: o.Ns, Dn: o.Dn, M: o.M, regime: `${inp.unit} formula`, xi: NaN };
    };
    const hDes = inp.hToe + inp.slr;
    const Hs = inp.depthLimited ? inp.Hs * hDes / inp.hToe : inp.Hs;
    const req = run(() => required(Hs, hDes), 'Armour stability');
    finite(req.M, 'Required armour mass');
    // Hudson KD: class table for rock, cubes and tetrapods; for a single-layer unit the identity Ns = (KD cot α)^⅓ is used
    const kdKey = inp.unit === 'cube' ? 'modified-cube' : inp.unit === 'tetrapod' ? 'tetrapod' : inp.kdClass;
    const KD = inp.KD > 0 ? inp.KD : inp.unit === 'single' ? inp.NsUnit ** 3 / inp.cotA : HUDSON_KD[kdKey][inp.breaking ? 'breaking' : 'nonBreaking'];
    const hud = hudson({ H: inp.hudsonH * Hs, KD, cotAlpha: inp.cotA, rhoS, rhoW: inp.rhoW });
    // damage of the armour actually provided (rock only: the Van der Meer formula inverted for S)
    const dmg = inp.unit === 'rock' ? run(() => vanDerMeerDamage({ Hs, Tm: inp.Tm, cotAlpha: inp.cotA, P: inp.P, N: inp.Nw, delta, Dn50: DnProv }), 'Damage of the provided armour') : null;
    const Sprov = dmg ? dmg.S : inp.Sd * Math.pow(req.Dn / DnProv, 5);
    const NsProv = Hs / (delta * DnProv);
    const displaced = Sprov * (1 - inp.nv) / DnProv; // stones displaced per metre of breakwater: Ae = S Dn50², volume per stone Dn50³

    const toe = run(() => toeStabilityVanDerMeer({ Hs, ht: inp.ht + inp.slr, h: hDes, delta: relativeDensity({ rhoS: rockRho, rhoW: inp.rhoW }), Nod: inp.Nod, rhoS: rockRho }), 'Toe stability');
    const layer = armourLayer({ Dn50: DnProv, n: inp.unit === 'single' ? 1 : inp.nLayers, kt: inp.kt, porosity: inp.nv, rhoS });
    const Mu = underlayerMass({ M50armour: M50prov }), DnU = massToDn({ M50: Mu, rhoS: rockRho }).Dn50;
    const filt = filterRules({ D15f: inp.D15f, D85b: inp.D85b, D15b: inp.D15b });
    const tU = 2 * DnU;
    need(layer.thickness + tU < inp.height * 0.8, 'Armour and underlayer are thicker than the mound allows: increase the mound height or reduce the armour size.');
    const qty = run(() => rubbleMoundQuantities({ crestWidth: inp.crestW, height: inp.height, cotFront: inp.cotA, cotRear: inp.cotR,
      layers: [{ name: 'Armour', thickness: layer.thickness, porosity: inp.nv, rhoS }, { name: 'Underlayer', thickness: tU, porosity: 0.35, rhoS: rockRho }], core: { porosity: 0.3, rhoS: rockRho } }), 'Section quantities (widen the crest or reduce the armour size so that the layers fit)');
    const perM = Object.fromEntries(qty.layers.map((l) => [l.name, l.mass / 1000]));
    const armourT = perM.Armour * inp.length, rockT = (concrete ? 0 : armourT) + (perM.Underlayer + perM.core) * inp.length;
    const concM3 = concrete ? qty.layers[0].solidVolume * inp.length : 0;
    const days = (rockT + (concrete ? armourT : 0)) / inp.placeRate;

    // deep-seated slip circle through the seaward slope and the foundation (Bishop's simplified method, buoyant weights)
    let slipFs = 99;
    if (inp.slip) {
      const Hm = inp.height, xc0 = inp.cotA * Hm, xr = xc0 + inp.crestW + inp.cotR * Hm;
      const surface = [[-3 * Hm, 0], [0, 0], [xc0, Hm], [xc0 + inp.crestW, Hm], [xr, 0], [xr + 3 * Hm, 0]];
      const soil = [{ yTop: Hm, gamma: 0.65 * (rockRho - inp.rhoW) * 9.80665, c: 0, phi: inp.phiMound * Math.PI / 180 }, { yTop: 0, gamma: inp.gamFound * 1e3, c: inp.cFound * 1e3, phi: inp.phiFound * Math.PI / 180 }];
      try { const cc = criticalCircle({ surface, soil, grid: 5, nSlices: 20, minDepth: 0.12 * Hm, xcRange: [-xc0, xc0], ycRange: [1.05 * Hm, 4 * Hm], yRigid: -Hm }); if (cc && Number.isFinite(cc.fs)) slipFs = cc.fs; } catch { slipFs = 99; }
    }
    const scen = slrRows(inp).map((s) => {
      const h = inp.hToe + s.rise, H = inp.depthLimited ? inp.Hs * h / inp.hToe : inp.Hs, q = required(H, h);
      const S = inp.unit === 'rock' ? vanDerMeerDamage({ Hs: H, Tm: inp.Tm, cotAlpha: inp.cotA, P: inp.P, N: inp.Nw, delta, Dn50: DnProv }).S : inp.Sd * Math.pow(q.Dn / DnProv, 5);
      return { name: s.name, rise: s.rise, Hs: r(H, 3), M: r(q.M / 1000, 3), uc: r(q.M / M50prov, 3), S: r(S, 3) };
    });
    const hs = linspace(0.5 * Hs, 1.5 * Hs, 31);
    const ok = req.M <= M50prov;
    const dM = req.M - M50prov;
    const armourDelta = armourT * (Math.cbrt(req.M / M50prov) - 1); // layer thickness ∝ Dn ∝ M^(1/3)

    const advice = [];
    advice.push(ok
      ? eng('ok', `The ${f(inp.M50)} t armour meets the ${f(req.M / 1000)} t required for Hs = ${f(Hs)} m (${req.regime}); predicted damage S = ${f(Sprov, 2)} against ${f(inp.Sd)} tolerated.`,
        `The reserve is ${f((1 - req.M / M50prov) * 100, 2)} % on mass. Before trimming it, confirm the design wave at the end-of-life sea level and test the section in a flume.`,
        'CIRIA C683 Rock Manual (2007) §5.2.2.2; Van der Meer (1988) — confirm current edition',
        `Reducing to the required mass would thin the armour layer and save about ${f(Math.max(0, -armourDelta), 3)} t of ${concrete ? 'units' : 'armour rock'} over ${f(inp.length)} m, at the cost of all reserve against sea-level rise.`)
      : eng('fail', `The ${f(inp.M50)} t armour is undersized: ${f(req.M / 1000)} t is required for Hs = ${f(Hs)} m (${req.regime}); predicted damage S = ${f(Sprov, 2)} against ${f(inp.Sd)} tolerated.`,
        `Increase the armour to at least ${f(req.M / 1000)} t (Dn50 ${f(req.Dn)} m), or flatten the slope, or add a toe berm to break the waves; for rock above the quarry yield switch to concrete units.`,
        'CIRIA C683 Rock Manual (2007) §5.2.2.2; Van der Meer (1988) — confirm current edition',
        `Heavier armour thickens the layer: about ${f(Math.max(0, armourDelta), 3)} t more ${concrete ? 'concrete units' : 'rock'} over ${f(inp.length)} m (≈ ${f(Math.max(0, armourDelta) * (concrete ? CO2_CONC / 2.35 : CO2_ROCK), 3)} t CO₂e, screening factor).`));
    const worst = scen[scen.length - 1];
    advice.push(eng(worst.uc > 1 ? 'warn' : 'info', `Sea-level rise: at +${f(worst.rise)} m (${worst.name}) the wave is ${f(worst.Hs)} m, the required mass ${f(worst.M)} t and the damage of the present armour S = ${f(worst.S, 2)}.`,
      worst.uc > 1 ? 'Either size the armour for the adopted end-of-life scenario now, or keep the section able to take a later overlay: a wide crest and a toe that can carry an extra layer.' : 'The armour stays stable through the scenarios; monitor profile surveys after major storms to confirm.', SLR_BASIS));
    if (slipFs < inp.sfSlip) advice.push(eng('fail', `Deep-seated slip: Bishop factor ${f(slipFs)} against ${f(inp.sfSlip)} required (foundation c = ${f(inp.cFound)} kPa, φ = ${f(inp.phiFound)}°).`, 'Flatten the slopes, add toe berms as counterweights, stage the construction to let the foundation consolidate, or replace or improve the soft layer.', 'Bishop (1955); EN 1997-1:2004 §11; CIRIA C683 §5.4.3.2', 'Counterweight berms use extra rock but avoid dredging and replacing the soft foundation, which has the larger seabed footprint.'));
    if (!filt.ok) advice.push(eng('warn', `The underlayer–core interface is not a closed filter (retention D15f/D85b = ${f(filt.retention.ratio)}, permeability D15f/D15b = ${f(filt.permeability.ratio)}).`, 'Regrade the underlayer, add an intermediate layer or a geotextile so that core material cannot wash out through the armour.', 'Terzaghi & Peck (1948) filter rules; CIRIA C683 §5.4.3.6'));
    advice.push(eco('info', `The section takes ${f(perM.Armour + perM.Underlayer + perM.core, 3)} t per metre (${f(rockT, 3)} t of rock${concrete ? ` and ${f(concM3, 3)} m³ of concrete units` : ''} in total) — about ${f(days, 3)} placing days at ${f(inp.placeRate)} t/day.`,
      'Match the armour grading to what the quarry yields: oversize blocks drive quarry waste and haul cost. Compare a flatter rock slope against concrete units on total tonnage and vessel days.', 'AACE RP 18R-97 Class 4 quantity-based estimating; CIRIA C683 ch. 3 (quarry yield) and ch. 9 (construction)',
      `Placing burns roughly ${f(days * MGO_T_PER_DAY, 3)} t of marine gas oil (≈ ${f(days * MGO_T_PER_DAY * CO2_MGO, 3)} t CO₂); rock supply adds about ${f(rockT * CO2_ROCK, 3)} t CO₂e${concrete ? ` and the concrete units about ${f(concM3 * CO2_CONC, 3)} t CO₂e` : ''} (screening factors). Rock can be recovered and re-used at end of life.`));

    return {
      kpis: [
        { key: 'M50req', label: 'Required armour mass', value: r(req.M / 1000), unit: 't', explain: 'Median rock mass (or unit mass) that keeps damage at the tolerated level in the design sea state. The provided armour should be at least this heavy.' },
        { key: 'Dn50req', label: 'Required nominal diameter Dn50', value: r(req.Dn), unit: 'm', explain: 'Cube-equivalent size of the required armour, (M50/ρs)^⅓. It sets the layer thickness and the grading class.' },
        { key: 'Ns', label: 'Stability number Hs/(Δ Dn50) required', value: r(req.Ns), unit: '-', explain: 'Dimensionless wave height the armour can stand at the tolerated damage. Rock typically 1.5–2.5; interlocking units up to about 3.' },
        { key: 'ucMass', label: 'Required / provided armour mass', value: r(req.M / M50prov, 3), unit: '-', limit: 1, sense: 'max', explain: 'Above 1.0 the armour is too light for the design wave. Because mass varies with H³ this ratio moves fast with wave height.' },
        { key: 'S', label: 'Damage level of the provided armour', value: r(Sprov, 3), unit: '-', limit: inp.Sd, sense: 'max', explain: 'Eroded cross-section area divided by Dn50². Around 2 is start of damage, 8 and above means the underlayer is exposed (failure). For concrete units it is scaled from the size ratio.' },
        { key: 'displaced', label: 'Stones displaced per metre run', value: r(displaced, 3), unit: '1/m', explain: 'Estimate of armour stones moved out of the eroded zone per metre of breakwater in the design storm; the maintenance volume to restore.' },
        { key: 'M50Hudson', label: 'Hudson armour mass', value: r(hud.M50 / 1000), unit: 't', explain: 'Classical Hudson no-damage mass with the stated KD, given as an independent cross-check. It ignores period, permeability and storm duration.' },
        { key: 'NsProv', label: 'Stability number of the provided armour', value: r(NsProv, 3), unit: '-', explain: 'Hs/(Δ Dn50) with the armour as placed. Compare with the required stability number: higher means more loaded.' },
        { key: 'DnToe', label: 'Toe-berm rock Dn50', value: r(toe.Dn50), unit: 'm', explain: 'Rock size that keeps the toe berm in place. A lost toe lets the armour slide; shallower toes need larger rock.' },
        { key: 'tArmour', label: 'Armour layer thickness', value: r(layer.thickness), unit: 'm', explain: 'n kt Dn50 for the armour provided. Used for crest levels, quantities and the underlayer geometry.' },
        { key: 'fsSlip', label: 'Factor of safety against deep-seated slip', value: r(Math.min(slipFs, 99)), unit: '-', limit: inp.sfSlip, sense: 'min', explain: 'Lowest Bishop factor of safety found for circular slips through the seaward slope and the foundation soil (buoyant weights, static). 99 when the check is switched off or no valid circle exists.' },
        { key: 'rockPerM', label: 'Material per metre of breakwater', value: r(perM.Armour + perM.Underlayer + perM.core), unit: 't/m', explain: 'Mass of armour, underlayer and core in one metre of the trunk section. Multiply by length for supply and placing.' },
        { key: 'rockTotal', label: 'Rock for the whole length', value: r(rockT), unit: 't', explain: 'Quarried rock in armour (if rock), underlayer and core over the stated length; drives quarry output, haulage and vessel days.' },
      ],
      checks: [
        check('Armour mass', req.M / 1000, inp.M50, 't', inp.unit === 'rock' ? 'Van der Meer (1988); Rock Manual (2007) Eq. 5.136–5.137' : inp.unit === 'single' ? 'Stability number from licensor guidance' : 'Van der Meer (1988) cubes / tetrapods', 'Guards against armour being rolled or lifted out by the design waves. If exceeded, units are displaced, the underlayer is exposed and the mound unravels.'),
        check('Hudson cross-check', hud.M50 / 1000, inp.M50, 't', 'Hudson (1959); SPM (1984) Eq. 7-116', 'Independent no-damage estimate. A large disagreement with Van der Meer signals long periods, low permeability or a long storm that Hudson cannot represent.'),
        check('Damage level', Sprov, inp.Sd, '-', 'Rock Manual (2007) Table 5.23', 'Guards against progressive erosion of the armour profile. Beyond the failure level the filter layer is exposed and washes out.'),
        { name: 'Deep-seated slip of the seaward slope', demand: inp.sfSlip, capacity: r(Math.min(slipFs, 99)), unit: '-', uc: r(Math.min(inp.sfSlip / Math.max(slipFs, 1e-6), 99)), ref: 'Bishop (1955); EN 1997-1 §11; Rock Manual (2007) §5.4.3', explain: 'Guards against a rotational slide of the mound and the soil beneath it. Soft foundations fail this way during construction, before any storm.' },
        check('Filter retention D15f / D85b', filt.retention.ratio, 5, '-', 'Terzaghi & Peck (1948); Rock Manual (2007) §5.4.3.6', 'Guards against core material being drawn out through the underlayer by wave action, which settles the crest and loosens the armour.'),
      ],
      plots: [
        { type: 'line', title: 'Required armour mass versus wave height', x: { label: 'Hs at toe', unit: 'm' }, y: { label: 'M50', unit: 't' },
          series: [{ name: 'Required', x: hs.map((v) => r(v, 3)), y: hs.map((v) => r(required(v, hDes).M / 1000, 3)) }, { name: 'Hudson', x: hs.map((v) => r(v, 3)), y: hs.map((v) => r(hudson({ H: inp.hudsonH * v, KD, cotAlpha: inp.cotA, rhoS, rhoW: inp.rhoW }).M50 / 1000, 3)) }, { name: 'Provided', x: [hs[0], hs[30]], y: [inp.M50, inp.M50], dash: true }],
          explain: 'Armour mass needed as the design wave changes. The crossing with the dashed line is the largest wave the present armour can take at the tolerated damage; the cubic rise shows why the design wave must be well founded.' },
        ...(inp.unit === 'rock' ? [{ type: 'line', title: 'Damage of the provided armour', x: { label: 'Hs at toe', unit: 'm' }, y: { label: 'Damage S', unit: '-' },
          series: [{ name: 'S', x: hs.map((v) => r(v, 3)), y: hs.map((v) => r(Math.min(40, vanDerMeerDamage({ Hs: v, Tm: inp.Tm, cotAlpha: inp.cotA, P: inp.P, N: inp.Nw, delta, Dn50: DnProv }).S), 3)) }, { name: 'Tolerable', x: [hs[0], hs[30]], y: [inp.Sd, inp.Sd], dash: true }],
          explain: 'Damage grows with the fifth power of wave height, so a storm 20 % above design does about 2.5 times the damage. Read the overload capacity before failure (S ≈ 8–17).' }] : []),
        { type: 'bar', title: 'Section quantities per metre', categories: qty.layers.map((l) => l.name), series: [{ name: 'Mass (t/m)', y: qty.layers.map((l) => r(l.mass / 1000, 3)) }],
          explain: 'Mass of each layer in one metre of breakwater. The core dominates the tonnage; the armour dominates the unit cost and the placing time.' },
      ],
      tables: [
        { title: 'Sea-level-rise pathway', columns: [{ key: 'name', label: 'Scenario' }, { key: 'rise', label: 'Rise', unit: 'm' }, { key: 'Hs', label: 'Hs', unit: 'm' }, { key: 'M', label: 'Required mass', unit: 't' }, { key: 'uc', label: 'Required / provided', unit: '-' }, { key: 'S', label: 'Damage S', unit: '-' }],
          rows: scen, explain: 'Design wave, required armour mass and predicted damage of the present armour for each sea-level-rise scenario. A ratio above 1 marks the scenario at which strengthening is needed.' },
        { title: 'Layer build-up', columns: [{ key: 'l', label: 'Layer' }, { key: 'M', label: 'M50', unit: 'kg' }, { key: 'Dn', label: 'Dn50', unit: 'm' }, { key: 't', label: 'Thickness', unit: 'm' }, { key: 'm', label: 'Mass', unit: 't/m' }],
          rows: [{ l: 'Armour (provided)', M: r(M50prov), Dn: r(DnProv, 3), t: r(layer.thickness, 3), m: r(perM.Armour, 3) }, { l: 'Underlayer (M50/10)', M: r(Mu), Dn: r(DnU, 3), t: r(tU, 3), m: r(perM.Underlayer, 3) }, { l: 'Core', M: 0, Dn: 0, t: 0, m: r(perM.core, 3) }, { l: 'Toe berm (required)', M: r(toe.M50), Dn: r(toe.Dn50, 3), t: r(2 * toe.Dn50, 3), m: 0 }],
          explain: 'Rock sizes and thicknesses from the armour down to the core, with the mass of each per metre run. Use it to choose standard gradings (EN 13383) and to check constructability.' },
      ],
      notes: [
        `${req.regime}; relative density Δ = ${f(delta)}; design Hs = ${f(Hs)} m at h = ${f(hDes)} m including ${f(inp.slr)} m of sea-level rise.`,
        `Hudson uses KD = ${f(KD)} with H = ${f(inp.hudsonH)} Hs. ${ok ? 'The provided armour governs with reserve.' : `The armour mass governs: short by ${f(dM / 1000)} t.`}`,
      ],
      advice,
      publish: { 'marine.armourM50': M50prov, 'marine.armourDn50': DnProv, 'quantities.rock_t': rockT },
      quantities: { rock_t: r(rockT), concrete_m3: r(concM3), vessel_days: { 'rock placing': r(days, 3) } },
    };
  },
  verify: [
    { name: 'Hudson formula M50 = ρs H³/(KD Δ³ cot α)', source: 'Hudson (1959); SPM (1984) Eq. 7-116. Hand calculation: H = 3 m, KD = 4, cot α = 2, ρs = 2650, ρw = 1025 → Δ = 1.5854, Δ³ = 3.9846, M50 = 2650 × 27/(4 × 3.9846 × 2) = 2244.6 kg',
      inputs: { Hs: 3, slr: 0, depthLimited: false, hToe: 20, ht: 15, KD: 4, hudsonH: 1, cotA: 2, rock: 'rock-granite', unit: 'rock' }, expect: [{ kpi: 'M50Hudson', value: 2.2446, tol: 0.002 }] },
    { name: 'Van der Meer plunging: Hs/(Δ Dn50) = 6.2 P^0.18 (S/√N)^0.2 ξm^−0.5', source: 'Van der Meer (1988); Rock Manual (2007) Eq. 5.136. Hand calculation: Hs = 3 m, Tm = 8 s, cot α = 2 → ξm = 2.885; P = 0.4, S = 2, N = 3000 → Ns = 6.2 × 0.8480 × 0.5158 × 0.5887 = 1.5966, Dn50 = 3/(1.5854 × 1.5966) = 1.185 m',
      inputs: { Hs: 3, Tm: 8, Nw: 3000, P: 0.4, Sd: 2, slr: 0, depthLimited: false, hToe: 20, ht: 15, cotA: 2, rock: 'rock-granite', unit: 'rock' }, expect: [{ kpi: 'Ns', value: 1.5966, tol: 0.003 }, { kpi: 'Dn50req', value: 1.185, tol: 0.003 }] },
  ],
};
