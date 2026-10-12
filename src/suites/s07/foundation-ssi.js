// Soil–structure interaction of a gravity foundation (caisson, pier footing, tower base): contact pressure,
// immediate and consolidation settlement in time, tilt, foundation springs, period lengthening and the
// wave-induced pore pressure in the seabed at the toe.
import { elasticSettlement, footingRotation, stressBelowRect, consolidationFD } from '../../engine/geotech/settlement.js';
import { wavePorePressure } from '../../engine/geotech/seepage.js';
import { momentaryLiquefaction } from '../../engine/geotech/liquefaction.js';
import { YEAR_S, CO2_CONC, need, run, finite, r, f, linspace, check, eng, eco, slrInputs, slrRows, SLR_BASIS } from './_common.js';

export default {
  id: 's07-foundation-ssi', suite: 7, title: 'Foundation settlement and soil–structure interaction',
  summary: 'Contact pressures, settlement with time, tilt and foundation stiffness of a gravity base on layered soil, the period shift it gives the structure and the wave-induced pore pressure at the toe.',
  covers: [7], fidelity: 'numerical',
  standards: ['EN 1997-1:2004 §6.6 and Annex F, H', 'DNV-RP-C212 (2021) §5–6', 'ISO 19901-4:2016 §7', 'Gazetas (1991) foundation vibrations; Veletsos & Meek (1974)', 'CIRIA C683 Rock Manual (2007) §5.4.5'],
  models: ['Boussinesq stress distribution', 'Elastic half-space settlement and rotation', 'Terzaghi one-dimensional consolidation (finite volumes, layered)', 'Darcy seepage', 'Rigid-footing impedance springs', 'Veletsos–Meek period lengthening', 'Wave-induced pore pressure (Yamamoto / Putnam)', 'Momentary liquefaction criterion'],
  limits: 'Linear-elastic soil with oedometric consolidation under the foundation centre-line: no plasticity, creep beyond primary consolidation, cyclic degradation or three-dimensional drainage, and static (zero-frequency) springs. Fully coupled wave–seabed poromechanics and finite-element soil models are outside the module; the pore-pressure result is the rigid-skeleton limit.',
  publishes: ['marine.foundKv', 'marine.foundKh', 'marine.foundKr', 'marine.settlement'],
  searchTerms: 'gravity foundation settlement consolidation tilt soil structure interaction springs wave induced pore pressure',
  inputs: [
    { key: 'B', label: 'Foundation width (in the plane of the moment)', unit: 'm', type: 'number', default: 15, min: 1, max: 100, group: 'Foundation', role: 'geometry', design: true, bind: 'marine.caissonWidth' },
    { key: 'L', label: 'Foundation length', unit: 'm', type: 'number', default: 30, min: 1, max: 500, group: 'Foundation', role: 'geometry' },
    { key: 'V', label: 'Net vertical load (buoyancy deducted)', unit: 'MN', type: 'number', default: 80, min: 0.1, max: 5000, group: 'Loads at the base', role: 'bc', uncertain: { dist: 'normal', cov: 0.05 } },
    { key: 'M', label: 'Overturning moment', unit: 'MN·m', type: 'number', default: 120, min: 0, max: 50000, group: 'Loads at the base', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'Hh', label: 'Horizontal load', unit: 'MN', type: 'number', default: 20, min: 0, max: 2000, group: 'Loads at the base', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'layers', label: 'Soil layers below the base', type: 'table', group: 'Ground', role: 'soil',
      columns: [{ key: 'name', label: 'Layer', type: 'text' }, { key: 'h', label: 'Thickness', unit: 'm', type: 'number' }, { key: 'E', label: 'Young\'s modulus', unit: 'MPa', type: 'number' }, { key: 'mv', label: 'Compressibility mv', unit: '1/MPa', type: 'number' }, { key: 'k', label: 'Permeability', unit: 'm/s', type: 'number' }],
      default: [{ name: 'Rubble bedding', h: 2, E: 120, mv: 0.008, k: 1e-2 }, { name: 'Dense sand', h: 6, E: 60, mv: 0.015, k: 1e-5 }, { name: 'Firm clay', h: 8, E: 25, mv: 0.12, k: 1e-10 }],
      help: 'Top to bottom from the underside of the foundation. mv from oedometer tests over the relevant stress range; permeability sets the consolidation rate.' },
    { key: 'nu', label: 'Poisson\'s ratio of the ground', unit: '-', type: 'number', default: 0.3, min: 0.1, max: 0.49, group: 'Ground', role: 'soil' },
    { key: 'drain', label: 'Drainage at the base of the profile', type: 'select', default: 'undrained', options: [{ value: 'undrained', label: 'Impermeable (rock, stiff clay)' }, { value: 'drained', label: 'Free-draining (sand, gravel)' }], group: 'Ground', role: 'bc', help: 'The top of the profile drains into the bedding layer.' },
    { key: 'gamSub', label: 'Submerged unit weight of the seabed', unit: 'kN/m³', type: 'number', default: 9.5, min: 4, max: 13, group: 'Ground', role: 'soil', bind: 'site.gammaSub' },
    { key: 'tRef', label: 'Time of interest after loading', unit: 'years', type: 'number', default: 2, min: 0.001, max: 200, group: 'Initial state and time', role: 'ic', help: 'For example the date of deck or crown-wall installation, or the end of the defects period.' },
    { key: 'life', label: 'Design service life', unit: 'years', type: 'number', default: 50, min: 1, max: 200, group: 'Initial state and time', role: 'criteria', bind: 'life.designYears' },
    { key: 'sPrior', label: 'Settlement already taken place', unit: 'mm', type: 'number', default: 0, min: 0, max: 5000, group: 'Initial state and time', role: 'ic', help: 'Existing settlement of an existing structure, added to the prediction.' },
    { key: 'sAllow', label: 'Allowable total settlement', unit: 'mm', type: 'number', default: 300, min: 5, max: 3000, group: 'Criteria', role: 'criteria' },
    { key: 'tiltAllow', label: 'Allowable tilt (1 in …)', unit: '-', type: 'number', default: 200, min: 20, max: 2000, group: 'Criteria', role: 'criteria' },
    { key: 'Tfix', label: 'Fixed-base natural period of the structure', unit: 's', type: 'number', default: 0.6, min: 0.02, max: 10, group: 'Structure', role: 'geometry' },
    { key: 'mass', label: 'Effective modal mass', unit: 't', type: 'number', default: 6000, min: 1, max: 1e6, group: 'Structure', role: 'geometry' },
    { key: 'heff', label: 'Effective height of the mass above the base', unit: 'm', type: 'number', default: 14, min: 0.5, max: 300, group: 'Structure', role: 'geometry' },
    { key: 'Hw', label: 'Design wave height at the toe', unit: 'm', type: 'number', default: 5, min: 0, max: 30, group: 'Wave boundary', role: 'bc', bind: 'marine.HmaxToe' },
    { key: 'Tw', label: 'Wave period', unit: 's', type: 'number', default: 10, min: 2, max: 25, group: 'Wave boundary', role: 'bc' },
    { key: 'd', label: 'Water depth at the toe', unit: 'm', type: 'number', default: 12, min: 1, max: 100, group: 'Wave boundary', role: 'bc', bind: 'marine.hToe' },
    ...slrInputs({ group: 'Wave boundary' }),
    { key: 'nCells', label: 'Cells per soil layer', type: 'int', default: 20, min: 4, max: 320, group: 'Solver', role: 'control' },
    { key: 'nt', label: 'Time steps', type: 'int', default: 200, min: 20, max: 6400, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nCells', label: 'Cells per layer', levels: [5, 10, 20, 40], qoi: 'Uref' },
  meshTime: { key: 'nt', label: 'Time steps', levels: [50, 100, 200, 400], qoi: 'Uref' },
  compute(inp) {
    const layers = (inp.layers || []).map((l, i) => ({ name: String(l.name || `Layer ${i + 1}`), h: Number(l.h), E: Number(l.E) * 1e6, mv: Number(l.mv) / 1e6, k: Number(l.k) }));
    need(layers.length >= 1, 'Give at least one soil layer below the foundation.');
    for (const l of layers) need(l.h > 0 && l.E > 0 && l.mv > 0 && l.k > 0, `Layer '${l.name}': thickness, modulus, compressibility and permeability must all be positive.`);
    const B = inp.B, L = Math.max(inp.L, B), V = inp.V * 1e6, M = inp.M * 1e6, Hh = inp.Hh * 1e6;
    const Htot = layers.reduce((a, l) => a + l.h, 0);
    // contact pressure: linear distribution, no tension
    const q = V / (B * L), e = M / V;
    need(e < B / 2, `The load eccentricity (${f(e)} m) is outside the base (B/2 = ${f(B / 2)} m): the foundation overturns. Widen the base or reduce the moment.`);
    const qMax = e <= B / 6 ? q * (1 + 6 * e / B) : 2 * V / (3 * L * (B / 2 - e)), qMin = e <= B / 6 ? q * (1 - 6 * e / B) : 0;
    // thickness-weighted modulus for the elastic estimates (harmonic: layers in series under vertical load)
    const E = Htot / layers.reduce((a, l) => a + l.h / l.E, 0), Gs = E / (2 * (1 + inp.nu));
    const si = elasticSettlement({ q, B, L, E, nu: inp.nu, shape: 'rect', point: 'rigid', H: Htot }).s;
    // primary consolidation: Σ mv Δσ h with the Boussinesq stress under the centre
    let sc = 0, z0 = 0; const prof = [], perLayer = [];
    for (const l of layers) {
      let sl = 0; const nSub = 10;
      for (let j = 0; j < nSub; j++) { const z = z0 + (j + 0.5) * l.h / nSub; const ds = stressBelowRect({ q, B, L, z }); sl += l.mv * ds * l.h / nSub; prof.push([z, ds]); }
      perLayer.push({ name: l.name, s: sl, cv: l.k / (l.mv * 9810) }); sc += sl; z0 += l.h;
    }
    // degree of consolidation in time from the layered 1-D solution (unit load; linear, so U is load-independent)
    const tRef = inp.tRef * YEAR_S, tLife = inp.life * YEAR_S;
    const fd = (tEnd) => run(() => consolidationFD({ layers: layers.map((l) => ({ thickness: l.h, k: l.k, mv: l.mv })), load: () => 1000, tEnd, nSteps: Math.round(inp.nt), nCells: Math.round(inp.nCells), top: 'drained', bottom: inp.drain }), 'Consolidation');
    const a = fd(tRef), bL = fd(tLife);
    const Uref = finite(a.U[a.U.length - 1], 'Degree of consolidation'), Ulife = bL.U[bL.U.length - 1];
    const prior = inp.sPrior / 1000;
    const sRef = prior + si + Uref * sc, sLife = prior + si + Ulife * sc, sFinal = prior + si + sc;
    const rot = footingRotation({ M, B, L, E, nu: inp.nu, shape: 'rect' });
    const tilt = Math.tan(rot.theta);
    // static springs of a rigid footing on a half-space through equivalent radii (Richart, Hall & Woods 1970; Gazetas 1991)
    const R0 = Math.sqrt(B * L / Math.PI), Rr = Math.pow(B ** 3 * L / (3 * Math.PI), 0.25);
    const Kv = 4 * Gs * R0 / (1 - inp.nu), Kh = 8 * Gs * R0 / (2 - inp.nu), Kr = 8 * Gs * Rr ** 3 / (3 * (1 - inp.nu));
    const kStr = 4 * Math.PI ** 2 * inp.mass * 1000 / inp.Tfix ** 2;
    const Tssi = inp.Tfix * Math.sqrt(1 + kStr / Kh + kStr * inp.heff ** 2 / Kr);
    const uh = Hh / Kh;
    // wave-induced pore pressure at the toe and depth of momentary liquefaction under the trough
    const pore = (rise) => {
      if (!(inp.Hw > 0)) return { p0: 0, zLiq: 0, k: 0 };
      const dd = inp.d + rise, w = wavePorePressure({ H: inp.Hw, d: dd, T: inp.Tw, z: 0, thickness: Htot });
      let zLiq = 0;
      for (const z of linspace(0.05, Math.min(Htot, 10), 60)) { if (momentaryLiquefaction({ p0: w.p0, z, gammaSub: inp.gamSub * 1e3, k: w.k }).liquefied) zLiq = z; else break; }
      return { p0: w.p0, zLiq, k: w.k };
    };
    const pp = pore(inp.slr);
    const scen = slrRows(inp).map((s) => { const p = pore(s.rise); return { name: s.name, rise: s.rise, p0: r(p.p0 / 1e3, 3), zLiq: r(p.zLiq, 3) }; });
    const okS = sLife * 1000 <= inp.sAllow, okT = tilt <= 1 / inp.tiltAllow;
    const tiltDen = tilt > 0 ? 1 / tilt : 1e6;
    // widening needed for the tilt criterion: rotation ∝ M/(B² L) for a given moment
    const Breq = okT ? B : B * Math.sqrt(tilt * inp.tiltAllow);
    const dConc = (Breq - B) * L * 1.0; // 1 m thick base slab extension, screening

    const ts = Array.from(bL.t), sT = Array.from(bL.settlement);
    const idx = linspace(0, ts.length - 1, Math.min(120, ts.length)).map((v) => Math.round(v));
    const advice = [
      okS && okT
        ? eng('ok', `Settlement is ${f(sLife * 1000)} mm at ${f(inp.life)} years (${f(sFinal * 1000)} mm final; ${f(si * 1000)} mm immediate) against ${f(inp.sAllow)} mm, and tilt 1 in ${f(tiltDen)} against 1 in ${f(inp.tiltAllow)}.`,
          `${f(Uref * 100, 2)} % of consolidation is complete after ${f(inp.tRef)} years: allow ${f((sLife - sRef) * 1000)} mm of further settlement in crest levels and joints placed at that time. Verify with settlement plates during construction.`,
          'EN 1997-1:2004 §6.6.2 and Annex F; DNV-RP-C212 (2021) §5',
          'Building the crest high by the predicted post-construction settlement avoids a later raising campaign and its concrete.')
        : eng('fail', `${okS ? '' : `Settlement reaches ${f(sLife * 1000)} mm at ${f(inp.life)} years against ${f(inp.sAllow)} mm allowed. `}${okT ? '' : `Tilt is 1 in ${f(tiltDen)} against 1 in ${f(inp.tiltAllow)} allowed (differential ${f(rot.differential * 1000)} mm across the base).`}`,
          `${okS ? '' : `Preload or stage the construction, replace or improve the compressible layer (${perLayer.slice().sort((x, y) => y.s - x.s)[0].name} gives ${f(perLayer.slice().sort((x, y) => y.s - x.s)[0].s * 1000)} mm), or over-build the crest by the predicted settlement. `}${okT ? '' : `Widen the base to about ${f(Breq)} m or reduce the moment to limit rotation.`}`,
          'EN 1997-1:2004 §6.6.2 and Annex H (limiting deformations); DNV-RP-C212 (2021) §5',
          `Widening the base slab to ${f(Breq)} m adds about ${f(Math.max(0, dConc), 3)} m³ of concrete (≈ ${f(Math.max(0, dConc) * CO2_CONC, 3)} t CO₂e, screening factor); preloading with reusable fill achieves the settlement target with almost no permanent material.`),
      eng(Tssi / inp.Tfix > 1.15 ? 'warn' : 'info', `Foundation springs: Kv = ${f(Kv / 1e6)} MN/m, Kh = ${f(Kh / 1e6)} MN/m, Kr = ${f(Kr / 1e9)} GN·m/rad. Soil flexibility lengthens the natural period from ${f(inp.Tfix)} s to ${f(Tssi)} s (× ${f(Tssi / inp.Tfix, 3)}).`,
        'Use these springs under the structural model for wave and seismic response; if the lengthened period approaches the wave period or the site period, run a dynamic analysis with frequency-dependent impedance.', 'Gazetas (1991); Veletsos & Meek (1974); ISO 19901-4:2016 §7; EN 1998-5 §6'),
      eng(pp.zLiq > 0.3 ? 'warn' : 'info', `Wave-induced bed pressure at the toe is ±${f(pp.p0 / 1e3)} kPa; momentary liquefaction under the trough reaches ${f(pp.zLiq, 2)} m below the seabed${scen.length ? ` (${f(scen[scen.length - 1].zLiq, 2)} m at +${f(scen[scen.length - 1].rise)} m sea level)` : ''}.`,
        pp.zLiq > 0.3 ? 'Protect the toe with a filter-compatible rock berm at least as thick as the liquefied depth, and check partly saturated (gassy) beds with a poro-elastic analysis.' : 'No significant loss of effective stress at the toe in the rigid-skeleton limit; confirm for loose or gassy sediments.', `Yamamoto et al. (1978); Sumer & Fredsøe (2002) ch. 10; ${SLR_BASIS}`),
      eco('info', `Post-construction settlement of ${f((sLife - sRef) * 1000)} mm must be absorbed by crest allowance or later raising over the ${f(inp.life)}-year life; ${f(Ulife * 100, 2)} % of consolidation occurs within it.`,
        'Compare a construction-stage preload or a longer wait before crown-wall casting with the whole-life cost of re-levelling; instrument the base so that inspection is triggered by measured, not assumed, settlement.', 'Whole-life costing per ISO 15686-5; observational method per EN 1997-1 §2.7',
        'An observational approach uses settlement data to avoid over-designing the foundation, saving concrete and ground-improvement material.'),
    ];

    return {
      kpis: [
        { key: 'sLife', label: 'Settlement at the end of the service life', value: r(sLife * 1000), unit: 'mm', limit: inp.sAllow, sense: 'max', explain: 'Existing settlement plus immediate settlement plus the consolidation reached at the end of the design life. It lowers the crest and strains joints and services.' },
        { key: 'sRef', label: 'Settlement at the time of interest', value: r(sRef * 1000), unit: 'mm', explain: 'Settlement accumulated by the stated time after loading. Anything built at that date sees only the remainder.' },
        { key: 'sFinal', label: 'Final settlement (immediate + full consolidation)', value: r(sFinal * 1000), unit: 'mm', explain: 'Long-term settlement when all excess pore pressure has dissipated; the upper bound for primary consolidation.' },
        { key: 'si', label: 'Immediate (elastic) settlement', value: r(si * 1000), unit: 'mm', explain: 'Settlement that occurs as the load is applied, from elastic distortion of the ground. It is built out during construction.' },
        { key: 'Uref', label: 'Degree of consolidation at the time of interest', value: r(Uref * 100), unit: '%', explain: 'Share of primary consolidation completed by the stated time, from the layered diffusion solution. Low values mean most settlement is still to come.' },
        { key: 'tilt', label: 'Tilt (1 in …)', value: r(Math.min(tiltDen, 1e6)), unit: '-', limit: inp.tiltAllow, sense: 'min', explain: 'Rigid-body rotation of the base under the moment, expressed as 1 in N. Smaller N is worse: it opens caisson joints and is visible above about 1 in 250.' },
        { key: 'qMax', label: 'Maximum contact pressure', value: r(qMax / 1e3), unit: 'kPa', explain: 'Peak pressure under the loaded edge from the vertical load and moment (linear distribution, no tension). Compare with the bearing resistance.' },
        { key: 'qMin', label: 'Minimum contact pressure', value: r(qMin / 1e3), unit: 'kPa', explain: 'Pressure under the unloaded edge. Zero means part of the base has lifted off.' },
        { key: 'Kv', label: 'Vertical foundation stiffness', value: r(Kv / 1e6), unit: 'MN/m', explain: 'Static vertical spring of the rigid base on the ground, for use under a structural model.' },
        { key: 'Kh', label: 'Horizontal foundation stiffness', value: r(Kh / 1e6), unit: 'MN/m', explain: 'Static sliding spring of the base. The horizontal load divided by it is the elastic base displacement.' },
        { key: 'Kr', label: 'Rocking foundation stiffness', value: r(Kr / 1e6), unit: 'MN·m/rad', explain: 'Static rocking spring about the long axis. It usually controls the period lengthening of tall structures.' },
        { key: 'Tssi', label: 'Natural period with soil flexibility', value: r(Tssi), unit: 's', explain: 'Fixed-base period lengthened by base sliding and rocking. A longer period changes the seismic demand and can approach wave periods.' },
        { key: 'uh', label: 'Elastic horizontal displacement of the base', value: r(uh * 1000), unit: 'mm', explain: 'Recoverable sliding of the base under the horizontal load, from the horizontal spring.' },
        { key: 'zLiq', label: 'Depth of momentary liquefaction at the toe', value: r(pp.zLiq, 3), unit: 'm', explain: 'Depth to which the upward pore-pressure gradient under a wave trough exceeds the submerged weight of the soil. Soil in this zone loses strength for a moment in each wave.' },
      ],
      checks: [
        check('Total settlement', sLife * 1000, inp.sAllow, 'mm', 'EN 1997-1 §6.6.2, Annex H', 'Guards against loss of crest level and damage to joints, pavements and services. Excess settlement turns into overtopping and repair cost.'),
        check('Tilt', tilt, 1 / inp.tiltAllow, '-', 'EN 1997-1 Annex H; DNV-RP-C212 (2021) §5', 'Guards against visible lean, opening of caisson joints and misalignment of crane rails or bearings.'),
        check('Base in full contact (e ≤ B/6)', e, B / 6, 'm', 'EN 1997-1 §6.5.4', 'Keeps the whole base in compression. Beyond the middle third the heel lifts, contact pressure escalates and rocking stiffness falls.'),
      ],
      plots: [
        { type: 'line', title: 'Settlement with time', x: { label: 'Time', unit: 'years' }, y: { label: 'Settlement', unit: 'mm' },
          series: [{ name: 'Settlement', x: idx.map((i) => r(ts[i] / YEAR_S, 3)), y: idx.map((i) => r((prior + si + (sT[i] / bL.sFinal) * sc) * 1000, 3)) }, { name: 'Allowable', x: [0, inp.life], y: [inp.sAllow, inp.sAllow], dash: true }],
          explain: 'Total settlement from loading to the end of the service life. The knee shows when consolidation is essentially complete; settlement to the right of the construction date must be allowed for in levels.' },
        { type: 'line', title: 'Stress increase and excess pore pressure with depth', x: { label: 'Stress increase (kPa) / remaining excess pore pressure (% of load)', unit: '' }, y: { label: 'Depth below base', unit: 'm' },
          series: [{ name: 'Boussinesq Δσ under centre', x: prof.map((p) => r(p[1] / 1e3, 3)), y: prof.map((p) => -r(p[0], 3)) }, { name: `Excess pore pressure at ${f(inp.tRef)} yr`, x: Array.from(a.u, (u) => r(u / 10, 3)), y: Array.from(a.z, (z) => -r(z, 3)) }],
          explain: 'How the foundation pressure spreads with depth and how much excess pore pressure is left at the time of interest. Layers still carrying pore pressure are the ones still settling.' },
        { type: 'bar', title: 'Consolidation settlement by layer', categories: perLayer.map((l) => l.name), series: [{ name: 'mm', y: perLayer.map((l) => r(l.s * 1000, 3)) }],
          explain: 'Contribution of each layer to the final consolidation settlement. The largest bar is the layer to improve, replace or preload.' },
      ],
      tables: [
        { title: 'Layer response', columns: [{ key: 'name', label: 'Layer' }, { key: 's', label: 'Final consolidation', unit: 'mm' }, { key: 'cv', label: 'cv', unit: 'm²/yr' }],
          rows: perLayer.map((l) => ({ name: l.name, s: r(l.s * 1000, 3), cv: r(l.cv * YEAR_S, 3) })), explain: 'Settlement contribution and coefficient of consolidation of each layer; low cv layers control how long settlement continues.' },
        { title: 'Sea-level-rise scenarios: seabed pore pressure at the toe', columns: [{ key: 'name', label: 'Scenario' }, { key: 'rise', label: 'Rise', unit: 'm' }, { key: 'p0', label: 'Bed pressure amplitude', unit: 'kPa' }, { key: 'zLiq', label: 'Liquefied depth', unit: 'm' }],
          rows: scen, explain: 'Wave pressure on the seabed and the depth of momentary liquefaction for the same wave in deeper water. Deeper water attenuates the pressure unless the wave itself grows.' },
      ],
      notes: [`Equivalent modulus ${f(E / 1e6)} MPa over ${f(Htot)} m; eccentricity ${f(e)} m (${e <= B / 6 ? 'inside' : 'outside'} the middle third).`, `${perLayer.slice().sort((x, y) => y.s - x.s)[0].name} governs the consolidation settlement.`],
      advice,
      publish: { 'marine.foundKv': Kv, 'marine.foundKh': Kh, 'marine.foundKr': Kr, 'marine.settlement': sLife },
      quantities: { concrete_m3: r(Math.max(0, dConc)) },
    };
  },
  verify: [
    { name: 'Terzaghi consolidation: U = 50 % at Tv = 0.197', source: 'Terzaghi (1943); Craig, Soil Mechanics ch. 7. Single layer H = 5 m, single drainage, k = 1e-9 m/s, mv = 1 /MPa → cv = k/(mv γw) = 1.0194e-7 m²/s, t50 = 0.197 × 25/cv = 4.831e7 s = 1.531 years',
      inputs: { layers: [{ name: 'Clay', h: 5, E: 20, mv: 1, k: 1e-9 }], drain: 'undrained', tRef: 1.531, nCells: 40, nt: 400 }, expect: [{ kpi: 'Uref', value: 50, tol: 0.02 }] },
    { name: 'Rigid footing springs on a half-space, Kv = 4 G R/(1 − ν), Kh = 8 G R/(2 − ν)', source: 'Richart, Hall & Woods (1970); Gazetas (1991). E = 50 MPa, ν = 0.3 → G = 19.23 MPa; B = L = 10 m → R = √(100/π) = 5.642 m; Kv = 4 × 19.23 × 5.642/0.7 = 620.0 MN/m; Kh = 8 × 19.23 × 5.642/1.7 = 510.6 MN/m',
      inputs: { B: 10, L: 10, nu: 0.3, layers: [{ name: 'Sand', h: 20, E: 50, mv: 0.02, k: 1e-4 }] }, expect: [{ kpi: 'Kv', value: 620.0, tol: 0.002 }, { kpi: 'Kh', value: 510.6, tol: 0.002 }] },
  ],
};
