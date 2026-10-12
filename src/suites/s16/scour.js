// Scour and seabed erosion at a monopile: sediment mobility, equilibrium scour depth in current, waves and their
// combination, time development, extent of the scour hole, and the design of a rock scour protection (stone size
// by the Shields criterion, filter rules, layer thickness and quantities).
import { keuleganCarpenter, pileScour, scourBreusers, scourTimeScale, scourDevelopment, criticalShields, bedShearStress, rockSizeShields, filterCriteria, armourTakeoff } from '../../engine/geotech/scour.js';
import { wavenumber } from '../../engine/metocean/waves.js';
import { adv, check, worst, sev, fmt, round, run, need, linspace, G, DEG, LIMIT_CODE, CONFIRM, CO2_ROCK, vesselCo2Text } from './_common.js';

const prot = { key: 'protect', in: [true] };
const sea = (key, label, def, bind, extra = {}) => ({ key, label, type: 'number', default: def, min: 0, group: 'Hydrodynamic conditions', role: 'bc', bind, ...extra });

export default {
  id: 's16-scour', suite: 16, title: 'Scour and scour protection at the foundation',
  summary: 'Seabed mobility, equilibrium and time-dependent scour depth around a monopile in current and waves, and the stone size, filter, thickness and tonnage of a rock scour protection.',
  covers: [17], fidelity: 'analytical',
  standards: ['DNV-ST-0126 (2021) §7.2 and App. (scour and scour protection)', 'DNV-RP-0618 (2022) Rock scour protection for monopiles', 'Sumer & Fredsøe (2002) The Mechanics of Scour in the Marine Environment', 'Soulsby (1997) Dynamics of Marine Sands', 'CIRIA C683 (2007) The Rock Manual §5.2 and §5.4', 'Whitehouse (1998) Scour at Marine Structures'],
  models: ['Shields mobility criterion (Soulsby critical Shields curve)', 'Current, wave and combined bed shear stress (Soulsby 1997)', 'Airy near-bed orbital velocity', 'Equilibrium pile scour in waves and current (Sumer & Fredsøe 2001)', 'Breusers pier-scour formula', 'Exponential scour development with the Sumer time scale', 'Stone sizing by amplified bed shear stress', 'Terzaghi filter rules', 'Annular armour-layer quantity take-off'],
  limits: `Empirical formulas for an isolated, slender, vertical circular pile on a non-cohesive, live-bed seabed: cohesive or layered soils, global scour and seabed-level change from sand waves, edge scour around the protection and pile groups are not modelled. Stone sizing is a static-stability (Shields) estimate with a user amplification factor; dynamic-stability design and physical model tests are current practice for final design. ${LIMIT_CODE}`,
  publishes: ['scour.depth', 'scour.depthUnprotected', 'quantities.rock_t', 'scour.rockD50'],
  searchTerms: 'monopile scour depth Sumer Fredsoe scour protection rock armour design',
  inputs: [
    { key: 'D', label: 'Pile diameter', unit: 'm', type: 'number', default: 6.5, min: 0.3, max: 15, group: 'Structure', role: 'geometry', design: true },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 30, min: 3, max: 100, group: 'Site', role: 'environment', bind: 'site.depth' },
    { key: 'd50', label: 'Median grain size of the seabed d50', unit: 'mm', type: 'number', default: 0.25, min: 0.06, max: 20, group: 'Site', role: 'soil', bind: 'site.d50', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Fine sand 0.1–0.25 mm, medium sand 0.25–0.5 mm. Non-cohesive seabed assumed.' },
    { key: 'd85', label: 'Seabed d85', unit: 'mm', type: 'number', default: 0.4, min: 0.08, max: 40, group: 'Site', role: 'soil', help: 'For the filter (retention) rule.' },
    sea('Uc', 'Depth-averaged current, normal conditions', 0.6, 'site.Uc1', { unit: 'm/s', max: 4, uncertain: { dist: 'normal', cov: 0.15 }, help: 'Tidal current that acts most of the time: it drives the long-term (equilibrium) scour.' }),
    sea('Hs', 'Significant wave height, normal conditions', 1.5, undefined, { unit: 'm', max: 10 }),
    sea('Tp', 'Peak period, normal conditions', 6.5, undefined, { unit: 's', min: 2, max: 20 }),
    sea('UcX', 'Depth-averaged current, extreme', 1.2, 'site.Uc100', { unit: 'm/s', max: 5 }),
    sea('HsX', 'Significant wave height, extreme', 8.5, 'site.Hs100', { unit: 'm', max: 16, help: 'Design storm for the stability of the scour protection.' }),
    sea('TpX', 'Peak period, extreme', 12.5, 'site.Tp100', { unit: 's', min: 3, max: 22 }),
    { key: 'S0', label: 'Scour depth already present', unit: 'm', type: 'number', default: 0, min: 0, max: 20, group: 'Initial state', role: 'ic', help: 'Initial seabed state from the latest survey (0 at installation).' },
    { key: 'nSigma', label: 'Standard deviations added to the mean scour depth', unit: '-', type: 'number', default: 0, min: 0, max: 2, group: 'Design basis', role: 'criteria', help: 'S/D = 1.3 has a standard deviation of 0.7 in current; DNV-ST-0126 quotes 1.3 D as the default design value (0 σ).' },
    { key: 'allow', label: 'Scour depth allowed for in the foundation design', unit: 'm', type: 'number', default: 0, min: 0, max: 25, group: 'Design basis', role: 'criteria', help: `0 = the pile was designed for a protected seabed. An unprotected design normally allows 1.3 D. ${CONFIRM}` },
    { key: 'repose', label: 'Angle of repose of the seabed sand', unit: 'deg', type: 'number', default: 30, min: 20, max: 40, group: 'Design basis', role: 'soil' },
    { key: 'protect', label: 'Install rock scour protection', type: 'bool', default: true, group: 'Scour protection', role: 'operation' },
    { key: 'rock', label: 'Armour rock', type: 'select', source: 'materials:rock', default: 'rock-granite', group: 'Scour protection', role: 'material', showIf: prot },
    { key: 'rockD50', label: 'Armour stone size d50', unit: 'mm', type: 'number', default: 700, min: 20, max: 1500, group: 'Scour protection', role: 'geometry', showIf: prot, design: true, help: 'Median sieve size of the armour grading (e.g. 10–60 kg ≈ 200 mm, 60–300 kg ≈ 380 mm, 300–1000 kg ≈ 650 mm).' },
    { key: 'alpha', label: 'Bed shear-stress amplification at the pile', unit: '-', type: 'number', default: 2.0, min: 1, max: 6, group: 'Scour protection', role: 'criteria', showIf: prot, uncertain: { dist: 'normal', cov: 0.15 }, help: 'About 4 beside a circular pile in steady current; 1.5–2.5 in wave-dominated storm flow at small Keulegan–Carpenter number, as for large monopiles (Sumer & Fredsøe 2002; Whitehouse 1998). Governs the stone size.' },
    { key: 'extent', label: 'Outer radius of the protection / pile diameter', unit: '-', type: 'number', default: 2.5, min: 1, max: 6, group: 'Scour protection', role: 'geometry', showIf: prot, design: true, help: 'Armour typically extends 4–5 D in diameter (radius 2–2.5 D).' },
    { key: 'layers', label: 'Armour thickness in stone diameters', unit: '-', type: 'number', default: 3, min: 1, max: 6, group: 'Scour protection', role: 'geometry', showIf: prot, help: 'At least 2–3 d50.' },
    { key: 'filterD50', label: 'Filter-layer stone size d50 (0 = no filter layer)', unit: 'mm', type: 'number', default: 30, min: 0, max: 200, group: 'Scour protection', role: 'geometry', showIf: prot },
    { key: 'filterD15', label: 'Filter-layer d15', unit: 'mm', type: 'number', default: 1.5, min: 0.1, max: 100, group: 'Scour protection', role: 'geometry', showIf: prot, help: 'Wide-graded filters (e.g. 1–100 mm) are used on sand because a narrow stone grading cannot retain it.' },
    { key: 'filterD85', label: 'Filter-layer d85', unit: 'mm', type: 'number', default: 90, min: 1, max: 400, group: 'Scour protection', role: 'geometry', showIf: prot },
    { key: 'filterT', label: 'Filter-layer thickness', unit: 'm', type: 'number', default: 0.5, min: 0, max: 2, group: 'Scour protection', role: 'geometry', showIf: prot },
    { key: 'rockRate', label: 'Rock placement rate of the installation vessel', unit: 't/day', type: 'number', default: 6000, min: 200, max: 30000, group: 'Scour protection', role: 'operation', showIf: prot, help: 'Fall-pipe or side-stone-dumping vessel, net of loading and transit.' },
    { key: 'nPos', label: 'Foundation positions', type: 'int', default: 60, min: 1, max: 400, group: 'Scour protection', role: 'geometry' },
  ],
  compute(inp, ctx) {
    const D = inp.D, h = inp.depth, d = inp.d50 / 1e3, rho = 1025, s = 2.65;
    need(inp.d85 >= inp.d50, 'The seabed d85 cannot be smaller than d50.');
    const orbital = (Hs, Tp) => { if (!(Hs > 0)) return 0; const k = wavenumber({ T: Tp, d: h }); return Math.PI * Hs / (Tp * Math.sinh(k * h)); };
    const cr = criticalShields({ d, s });
    const cond = (Uc, Hs, Tp, name) => {
      const Um = orbital(Hs, Tp), KC = Um > 0 ? keuleganCarpenter(Um, Tp, D) : 0;
      const tau = bedShearStress({ d, h, U: Uc, Uw: Um, T: Tp, rho });
      const theta = tau.tauMax / (rho * (s - 1) * G * d);
      const Ucw = Uc + Um > 0 ? Uc / (Uc + Um) : 1;
      // waves present: Sumer & Fredsøe combined formula; pure current: S/D = 1.3
      const sc = Um > 1e-6 ? pileScour({ KC, Ucw, nSigma: 0 }) : pileScour({ KC: Infinity, Ucw: 1, nSigma: inp.nSigma });
      let SoverD = sc.SoverD;
      if (Um > 1e-6 && Ucw > 0.7) SoverD = Math.max(SoverD, 1.3 * Math.min(1, (Ucw - 0.7) / 0.3)) + inp.nSigma * 0.7 * Math.min(1, (Ucw - 0.7) / 0.3);   // blend towards the current value where the formula's KC threshold suppresses a current-dominated case
      const live = theta > cr.thetaCr;
      return { name, Uc, Hs, Tp, Um, KC, tau, theta, Ucw, SoverD: live ? SoverD : SoverD * Math.max(0, Math.min(1, 2 * theta / cr.thetaCr - 1)), live };
    };
    const normal = cond(inp.Uc, inp.Hs, inp.Tp, 'Normal (tide and operational waves)');
    const curOnly = cond(inp.Uc, 0, inp.Tp, 'Current only');
    const storm = cond(inp.UcX, inp.HsX, inp.TpX, 'Extreme storm');
    const tideMax = cond(Math.max(inp.Uc, inp.UcX), 0, inp.TpX, 'Extreme current only');
    const all = [normal, curOnly, storm, tideMax];
    const gov = all.reduce((a, c) => (c.SoverD > a.SoverD ? c : a), all[0]);
    const Seq = gov.SoverD * D, Sbr = scourBreusers({ D, h });
    // time scale in the governing condition
    const mode = gov.Um > 1e-6 && gov.Ucw < 0.5 ? 'waves' : 'current';
    let T = Infinity;
    if (gov.live && Seq > 0) { try { T = scourTimeScale({ mode, D, theta: Math.max(gov.theta, 1e-4), d, s, delta: h, KC: Math.max(gov.KC, 0.1) }).T; } catch { T = Infinity; } }
    const days = [0.25, 0.5, 1, 2, 4, 7, 14, 30, 90, 365], dev = days.map((dd) => (Number.isFinite(T) ? Math.max(inp.S0, scourDevelopment(Seq, dd * 86400, T)) : inp.S0));
    const holeR = D / 2 + Seq / Math.tan(inp.repose * DEG);
    // scour protection
    let rk = null, req = null, take = null, filt = null, takeF = null, stormRock = null;
    if (inp.protect) {
      rk = ctx.material(inp.rock);
      need(rk && rk.rho > 1500, `Rock '${inp.rock}' is not in the materials library.`);
      const UmX = orbital(inp.HsX, inp.TpX);
      req = run(() => rockSizeShields({ h, U: inp.UcX, Uw: UmX, T: inp.TpX, rho, amplification: inp.alpha, rhoS: rk.rho }), 'Stone sizing');
      const dR = inp.rockD50 / 1e3;
      stormRock = bedShearStress({ d: dR, h, U: inp.UcX, Uw: UmX, T: inp.TpX, rho });
      take = armourTakeoff({ D, extentFactor: inp.extent, thickness: inp.layers * dR, porosity: 0.4, rhoS: rk.rho });
      // filter rules between armour and the layer below (filter stone or seabed); narrow armour grading assumed, d15 ≈ 0.6 d50
      if (inp.filterD50 > 0) need(inp.filterD15 < inp.filterD50 && inp.filterD50 < inp.filterD85, 'Filter grading: d15 < d50 < d85 is required.');
      const below = inp.filterD50 > 0 ? { D15b: inp.filterD15 / 1e3, D50b: inp.filterD50 / 1e3, D85b: inp.filterD85 / 1e3 } : { D15b: 0.6 * d, D50b: d, D85b: inp.d85 / 1e3 };
      filt = { armour: filterCriteria({ D15f: 0.6 * dR, D50f: dR, ...below }), bed: inp.filterD50 > 0 ? filterCriteria({ D15f: inp.filterD15 / 1e3, D50f: inp.filterD50 / 1e3, D15b: 0.6 * d, D50b: d, D85b: inp.d85 / 1e3 }) : null };
      if (inp.filterD50 > 0 && inp.filterT > 0) takeF = armourTakeoff({ D, extentFactor: inp.extent + 0.3, thickness: inp.filterT, porosity: 0.4, rhoS: rk.rho });
    }
    const designScour = inp.protect ? inp.S0 : Math.max(Seq, inp.S0);
    const rockT = inp.protect ? (take.mass + (takeF ? takeF.mass : 0)) / 1e3 : 0, rockAll = rockT * inp.nPos;
    const vDays = inp.protect ? rockAll / inp.rockRate : 0;
    const checks = [
      check('Design scour depth against the allowance in the foundation design', designScour, Math.max(inp.allow, 0.05 * D), 'm', 'DNV-ST-0126 (2021) §7.2.3 (scour to be accounted for or prevented)', 'Scour that will develop (zero behind an intact protection) against the depth the pile design assumed. More scour than allowed lowers the natural frequency and the lateral capacity and raises the pile moment.'),
      ...(inp.protect ? [
        check('Armour stone size against the required size in the design storm', req.d50 * 1e3, inp.rockD50, 'mm', 'DNV-RP-0618 (2022) / CIRIA C683 §5.2 (Shields approach with amplified shear stress)', 'Demand is the smallest stone that stays put under the amplified bed shear stress of the design storm, capacity the stone size chosen. Smaller stones are winnowed away and the protection fails progressively.'),
        check('Filter rule armour over the layer below (retention D15f / D85b < 5)', filt.armour.retention.ratio, 5, '-', 'CIRIA C683 §5.4.3.6 (geometrically closed filter)', 'Ratio of the fine fraction of the armour to the coarse fraction of the layer beneath. Above 5 the finer material is sucked out through the voids and the armour sinks.'),
        ...(filt.bed ? [check('Filter rule filter layer over the seabed (retention D15f / D85b < 5)', filt.bed.retention.ratio, 5, '-', 'CIRIA C683 §5.4.3.6', 'Same retention rule between the filter stone and the seabed sand. If it fails, sand is pumped out through the filter by wave pressure and the whole protection settles.')] : []),
        check('Armour thickness against two stone diameters', 2, inp.layers, 'd50', 'DNV-RP-0618 (2022) §4 (minimum layer thickness)', 'A layer thinner than about two stone diameters leaves gaps after placement tolerances and the first stone loss exposes the filter.'),
      ] : []),
    ];
    const gc = worst(checks);
    const xs = linspace(0.2, 2.5, 24), sVsU = xs.map((u) => cond(u, inp.Hs, inp.Tp, '').SoverD * D);
    const m50 = inp.protect ? rk.rho * (0.84 * inp.rockD50 / 1e3) ** 3 : 0;
    return {
      kpis: [
        { key: 'S', label: 'Equilibrium scour depth (unprotected)', value: round(Seq, 4), unit: 'm', explain: `Depth of the scour hole that develops at the pile without protection in the governing condition (${gov.name}). The default for current-dominated sites is 1.3 pile diameters.` },
        { key: 'SoverD', label: 'Scour depth / pile diameter', value: round(gov.SoverD, 4), unit: '-', explain: 'Relative scour depth: 1.3 in steady current, much less in waves at low Keulegan–Carpenter number, as for large monopiles in storm waves.' },
        { key: 'designScour', label: 'Design scour depth', value: round(designScour, 4), unit: 'm', limit: round(Math.max(inp.allow, 0.05 * D), 4), sense: 'max', explain: 'Scour to carry into the foundation design: the surveyed depth behind an intact protection, or the equilibrium depth if the seabed is left unprotected.' },
        { key: 'KC', label: 'Keulegan–Carpenter number, normal waves', value: round(normal.KC, 4), unit: '-', explain: 'Orbital excursion relative to the pile diameter. Below about 6 waves alone do not scour a pile: the horseshoe vortex has no time to form.' },
        { key: 'KCstorm', label: 'Keulegan–Carpenter number, extreme storm', value: round(storm.KC, 4), unit: '-', explain: 'The same number in the design storm; it is still small for large monopiles, so storms tend to backfill a current-scoured hole.' },
        { key: 'theta', label: 'Shields parameter, normal conditions', value: round(normal.theta, 4), unit: '-', limit: undefined, explain: `Dimensionless bed shear stress of the undisturbed flow. Above the critical value of ${fmt(cr.thetaCr)} the seabed is mobile everywhere (live bed); below it only the amplified flow at the pile moves sand (clear water).` },
        { key: 'thetaCr', label: 'Critical Shields parameter of the seabed', value: round(cr.thetaCr, 4), unit: '-', explain: 'Threshold of sediment motion for this grain size (Soulsby 1997).' },
        { key: 'T90', label: 'Time to 90 % of the equilibrium depth', value: round(Number.isFinite(T) ? 2.303 * T / 86400 : 0, 4), unit: 'days', explain: 'Scour grows as 1 − exp(−t/T): 90 % is reached after 2.3 time scales (0 = no development in the governing condition). Often only days to weeks — protection must follow pile installation closely.' },
        { key: 'holeR', label: 'Radius of the scour hole from the pile axis', value: round(holeR, 4), unit: 'm', explain: 'Extent of the hole with side slopes at the angle of repose. Cables must be buried or protected across it.' },
        { key: 'rockReq', label: 'Required armour stone size d50', value: round(req ? req.d50 * 1e3 : 0, 4), unit: 'mm', explain: 'Smallest stone that is statically stable in the design storm with the amplified shear stress at the pile (0 when no protection is selected).' },
        { key: 'rockT', label: 'Rock per position (armour and filter)', value: round(rockT, 5), unit: 't', explain: 'Mass of armour and filter stone for one foundation from the layer geometry at 40 % porosity.' },
        { key: 'rockAll', label: 'Rock for all positions', value: round(rockAll, 5), unit: 't', explain: 'Total quarry tonnage for the farm: the quantity that sets the rock-installation campaign.' },
        { key: 'breusers', label: 'Scour depth by Breusers et al. (current)', value: round(Sbr, 4), unit: 'm', explain: 'Independent pier-scour estimate 1.5·D·tanh(h/D) for comparison with the Sumer & Fredsøe value.' },
      ].map((k) => { if (k.limit === undefined) delete k.limit; return k; }),
      checks,
      plots: [
        { type: 'line', title: 'Scour development with time', x: { label: 'Time', unit: 'days' }, y: { label: 'Scour depth', unit: 'm' }, series: [{ name: 'S(t)', x: days, y: dev.map((v) => round(v, 4)) }, { name: 'Equilibrium', x: [days[0], days[days.length - 1]], y: [Seq, Seq], dash: true }],
          explain: 'Growth of the scour hole in the governing condition. The steep start shows how soon after piling the protection must be in place.' },
        { type: 'line', title: 'Equilibrium scour depth versus current speed (normal waves)', x: { label: 'Depth-averaged current', unit: 'm/s' }, y: { label: 'Scour depth', unit: 'm' }, series: [{ name: 'S', x: xs.map((v) => round(v, 3)), y: sVsU.map((v) => round(v, 4)) }, { name: 'Normal current', x: [inp.Uc, inp.Uc], y: [0, Math.max(...sVsU, 0.1)], dash: true }],
          explain: 'Sensitivity of the scour depth to the tidal current with the normal wave climate. The curve rises from clear-water to live-bed conditions and saturates at 1.3 D.' },
        { type: 'line', title: 'Scour hole and protection cross-section', x: { label: 'Distance from the pile axis', unit: 'm' }, y: { label: 'Level relative to the seabed', unit: 'm' },
          series: [{ name: 'Unprotected scour hole', x: [D / 2, D / 2, round(holeR, 3), round(holeR + D, 3)], y: [0, round(-Seq, 3), 0, 0] }, ...(inp.protect ? [{ name: 'Rock armour', x: [D / 2, round(take.outerRadius, 3), round(take.outerRadius, 3)], y: [round(inp.layers * inp.rockD50 / 1e3 + (takeF ? inp.filterT : 0), 3), round(inp.layers * inp.rockD50 / 1e3 + (takeF ? inp.filterT : 0), 3), 0] }] : []), { name: 'Pile wall', x: [D / 2, D / 2], y: [-Seq - 1, 3], dash: true }],
          explain: 'Radial section: the scour hole that would form without protection and the rock layer that prevents it. The protection should cover the footprint of the potential hole.' },
        { type: 'bar', title: 'Scour depth by condition', categories: all.map((c) => c.name), series: [{ name: 'S (m)', y: all.map((c) => round(c.SoverD * D, 4)) }], y: { label: 'Scour depth', unit: 'm' },
          explain: 'Equilibrium scour for each hydrodynamic condition. For large monopiles the tidal current, not the storm, usually gives the deepest hole.' },
      ],
      tables: [
        { title: 'Hydrodynamic conditions and seabed mobility', columns: [{ key: 'name', label: 'Condition' }, { key: 'Uc', label: 'Current', unit: 'm/s' }, { key: 'Um', label: 'Orbital velocity', unit: 'm/s' }, { key: 'KC', label: 'KC', unit: '-' }, { key: 'Ucw', label: 'Uc/(Uc+Um)', unit: '-' }, { key: 'theta', label: 'Shields θ', unit: '-' }, { key: 'bed', label: 'Seabed' }, { key: 'S', label: 'Scour depth', unit: 'm' }],
          rows: all.map((c) => ({ name: c.name, Uc: round(c.Uc, 3), Um: round(c.Um, 4), KC: round(c.KC, 4), Ucw: round(c.Ucw, 3), theta: round(c.theta, 4), bed: c.live ? 'live bed' : 'clear water', S: round(c.SoverD * D, 4) })),
          explain: 'One row per condition: flow at the bed, Keulegan–Carpenter number, mobility of the undisturbed seabed and the resulting equilibrium scour depth.' },
        ...(inp.protect ? [{ title: 'Scour protection', columns: [{ key: 'q', label: 'Item' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [{ q: 'Armour d50 (chosen / required)', v: `${inp.rockD50} / ${fmt(req.d50 * 1e3)}`, u: 'mm' }, { q: 'Median stone mass M50', v: round(m50, 4), u: 'kg' }, { q: 'Armour thickness', v: round(inp.layers * inp.rockD50 / 1e3, 4), u: 'm' }, { q: 'Outer radius', v: round(take.outerRadius, 4), u: 'm' }, { q: 'Armour mass per position', v: round(take.mass / 1e3, 5), u: 't' }, { q: 'Filter mass per position', v: round(takeF ? takeF.mass / 1e3 : 0, 5), u: 't' },
            { q: 'Bed shear stress on the armour in the storm (unamplified)', v: round(stormRock.tauMax, 4), u: 'Pa' }, { q: 'Rock-installation vessel days (all positions)', v: round(vDays, 4), u: 'days' }],
          explain: 'Geometry, grading and quantities of the rock protection. The required stone size should be confirmed by a dynamic-stability method or model tests before ordering rock.' }] : []),
      ],
      notes: [
        `Governing condition: ${gov.name} — S/D = ${fmt(gov.SoverD)} (${fmt(Seq)} m); the undisturbed seabed is ${gov.live ? 'mobile (live bed)' : 'not mobile (clear-water scour)'}; Breusers gives ${fmt(Sbr)} m.`,
        `Storm: orbital velocity ${fmt(storm.Um)} m/s, KC = ${fmt(storm.KC)} — ${storm.KC < 6 ? 'too small for wave-induced scour at this diameter' : 'wave-induced scour is active'}.`,
        inp.protect ? `Protection: armour d50 ${inp.rockD50} mm (required ${fmt(req.d50 * 1e3)} mm), ${fmt(inp.layers * inp.rockD50 / 1e3)} m thick to a radius of ${fmt(take.outerRadius)} m; ${fmt(rockT)} t per position.` : 'No protection: the equilibrium scour depth must be included in the foundation design and monitored.',
        `Governing check: ${gc.name} (utilisation ${fmt(gc.uc)}).`,
      ],
      advice: [
        adv('engineering', sev(gc.uc), `Unprotected equilibrium scour is ${fmt(Seq)} m (${fmt(gov.SoverD)} D) and would be 90 % developed in ${Number.isFinite(T) ? fmt(2.303 * T / 86400) : '–'} days; design scour depth ${fmt(designScour)} m against ${fmt(inp.allow)} m allowed.`,
          gc.uc > 1 ? (gc === checks[0] ? `Either protect the seabed or redesign the pile for ${fmt(Seq)} m of scour (longer embedment, lower natural frequency).` : gc.name.startsWith('Armour stone') ? `Increase the armour to d50 ≥ ${fmt(req.d50 * 1e3 * 1.05)} mm (M50 ≈ ${fmt(rk.rho * (0.84 * req.d50 * 1.05) ** 3)} kg) or justify a lower amplification factor by model tests.` : `Adjust the gradings: '${gc.name}' is not met (ratio ${fmt(gc.demand)}); add or coarsen the filter layer.`) : inp.protect ? 'The protection is statically stable and the filter rules are met: confirm edge-scour allowance (a falling apron of extra rock at the perimeter) and survey after the first winter.' : 'The pile design covers the scour: set up a bathymetric monitoring plan and a trigger level for remedial rock.',
          'DNV-ST-0126 (2021) §7.2 and DNV-RP-0618 (2022); Sumer & Fredsøe (2002) for scour depth, CIRIA C683 for filters.',
          inp.protect ? `Rock for ${inp.nPos} positions: ${fmt(rockAll)} t ≈ ${fmt(rockAll * CO2_ROCK)} t CO₂e at the quay (screening factor ${CO2_ROCK} t/t), plus placement: ${vesselCo2Text(vDays, 35)} The rock covers about ${fmt(take.area * inp.nPos / 1e4)} ha of seabed and becomes hard-substrate habitat.` : `Without rock there is no quarrying or seabed cover, but each pile needs about ${fmt(Math.PI * (D - 0.07) * 0.07 * Seq * 7.85)} t of extra steel to carry ${fmt(Seq)} m of scour.`),
        adv('economic', 'info', inp.protect ? `The protection needs ${fmt(rockAll)} t of rock and about ${fmt(vDays)} vessel days for ${inp.nPos} positions.` : `Leaving the seabed unprotected saves the rock campaign but lengthens every pile for ${fmt(Seq)} m of scour.`,
          inp.protect ? `Compare with designing the piles for scour: ${fmt(Seq)} m of extra free length is roughly ${fmt(Math.PI * (D - 0.07) * 0.07 * Seq * 7.85 * inp.nPos)} t of pile steel for the farm against ${fmt(rockAll)} t of rock and ${fmt(vDays)} vessel days; a thinner armour (${inp.layers} → 2.5 d50) saves ${fmt(rockAll * (1 - 2.5 / inp.layers) * (take.mass / (take.mass + (takeF ? takeF.mass : 0))))} t.` : 'Budget for scour monitoring (annual multibeam survey) and a contingency rock campaign; the cost of late remedial rock is several times that of rock placed during installation.',
          'Whole-life costing per ISO 15686-5 (protection now versus steel plus monitoring); quantities for a Class 4 estimate per AACE RP 18R-97.'),
      ],
      publish: { 'scour.depth': designScour, 'scour.depthUnprotected': Seq, 'quantities.rock_t': rockAll, 'scour.rockD50': inp.protect ? inp.rockD50 / 1e3 : 0 },
      quantities: { rock_t: round(rockAll, 4), ...(vDays > 0 ? { vessel_days: { rock_installation: round(vDays, 3) } } : {}) },
    };
  },
  verify: [
    { name: 'Steady current on a live bed: S = 1.3·D', source: 'Sumer, Fredsøe & Christiansen (1992); DNV-ST-0126 (2021) default equilibrium scour depth',
      // D = 6.5 m → 8.45 m; d50 = 0.2 mm in a 1.2 m/s current is live bed (θ well above θcr)
      inputs: { D: 6.5, Uc: 1.2, UcX: 1.2, Hs: 0, HsX: 0, d50: 0.2, d85: 0.3, nSigma: 0, protect: false }, expect: [{ kpi: 'S', value: 8.45, tol: 0.002 }, { kpi: 'SoverD', value: 1.3, tol: 0.002 }] },
    { name: 'Waves only on a slender pile: S/D = 1.3·{1 − exp[−0.03(KC − 6)]}', source: 'Sumer, Fredsøe & Christiansen (1992) J. Waterway, Port, Coastal Ocean Eng. 118(1); hand calculation',
      // h = 20 m, Hs = 4 m, T = 10 s: k = 0.05183 rad/m, sinh(kh) = 1.2326 → Um = π·4/(10·1.2326) = 1.0195 m/s; D = 1 m → KC = 10.195
      // S/D = 1.3·(1 − exp(−0.12585)) = 0.1537 → S = 0.154 m
      inputs: { D: 1, depth: 20, Uc: 0, UcX: 0, Hs: 4, Tp: 10, HsX: 4, TpX: 10, d50: 0.2, d85: 0.3, nSigma: 0, protect: false }, expect: [{ kpi: 'KC', value: 10.195, tol: 0.01 }, { kpi: 'S', value: 0.1537, tol: 0.04 }] },
  ],
};
