// Seabed scour around a monopile: sediment mobility, equilibrium and time-dependent local scour, its effect on
// foundation stiffness, pile response and natural frequency, and the rock size and quantity of scour protection.
import { criticalShields, bedShearStress, pileScour, pileScourRudolphBos, scourBreusers, scourTimeScale, scourDevelopment, rockSizeShields, armourTakeoff, filterCriteria, keuleganCarpenter } from '../../engine/geotech/scour.js';
import { nearBedWaveCurrent } from '../../engine/metocean/current.js';
import { applyScour } from '../../engine/geotech/soilprofile.js';
import { solveLateralPile } from '../../engine/geotech/pile.js';
import { monopileStiffness } from '../../engine/geotech/monopile.js';
import { G, DEG, NU_SEA, need, run, fmt, round, adv, check, worst, sev, vesselCo2Text, turbineInputs, supportInputs, steelInputs, steelOf, soilInputs, soilProfile, turbine, LIMIT_PY } from './_common.js';
import { pileOptions } from './_pile.js';
import { supportStructure } from './_dynamics.js';

export default {
  id: 's03-scour', suite: 3, title: 'Seabed scour and scour protection',
  summary: 'Local scour depth at the monopile under current and waves and its development in time, the loss of foundation stiffness, the shift of natural frequency and pile response, and the armour rock size and tonnage for protection.',
  covers: [10], fidelity: 'screening',
  standards: ['DNV-ST-0126 (2021) §7.7 and Appendix D (scour, S = 1.3 D)', 'DNV-RP-0618 (2022) rock scour protection for monopiles', 'Sumer & Fredsøe (2002) The Mechanics of Scour in the Marine Environment', 'Soulsby (1997) Dynamics of Marine Sands', 'CIRIA C683 The Rock Manual §5.2–5.4', 'API RP 2A-WSD §6.3 commentary (overburden reduction)'],
  models: ['Shields critical mobility and wave–current bed shear stress (Soulsby)', 'Keulegan–Carpenter dependent equilibrium scour (Sumer & Fredsøe), Rudolph & Bos low-KC form, Breusers current scour', 'Scour time scale and exponential development', 'Overburden reduction below the scour hole in the p–y model', 'Three-spring foundation stiffness and beam eigenfrequency with scour', 'Critical-shear rock sizing and Terzaghi filter rules'],
  limits: `Non-cohesive seabed; clay erosion and layered beds are not predicted. Empirical scour formulae were derived for slender piles in live-bed conditions and carry a scatter of about ±0.7 D. Edge scour of the protection, winnowing and seabed lowering by sand waves are not modelled. ${LIMIT_PY} Morphological modelling and CFD are outside the module.`,
  publishes: ['quantities.rock_t'],
  searchTerms: 'monopile scour depth scour protection rock size natural frequency',
  inputs: [
    ...turbineInputs(),
    ...supportInputs(),
    ...steelInputs(),
    ...soilInputs(),
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'd50', label: 'Median grain size d50', unit: 'mm', type: 'number', default: 0.25, min: 0.03, max: 60, group: 'Seabed', role: 'soil', bind: 'site.d50', uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'd85', label: 'Grain size d85', unit: 'mm', type: 'number', default: 0.5, min: 0.05, max: 120, group: 'Seabed', role: 'soil', help: 'From the grain-size distribution; used for the filter retention rule.' },
    { key: 'rhoS', label: 'Grain density', unit: 'kg/m³', type: 'number', default: 2650, min: 1500, max: 3500, group: 'Seabed', role: 'soil' },
    { key: 'S0', label: 'Scour already present', unit: 'm', type: 'number', default: 0, min: 0, max: 30, group: 'Seabed', role: 'ic', help: 'Initial seabed profile at the pile from the latest bathymetric survey; scour develops from this depth.' },
    { key: 'Sg', label: 'General seabed lowering', unit: 'm', type: 'number', default: 0, min: 0, max: 15, group: 'Seabed', role: 'ic', help: 'Regional lowering (sand-wave migration, global scour) over the design life.' },
    { key: 'Uc', label: 'Near-bed current speed (1 m above bed)', unit: 'm/s', type: 'number', default: 0.8, min: 0, max: 3, group: 'Flow', role: 'bc', bind: 'site.UcBed100', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'Hs', label: 'Significant wave height with that current', unit: 'm', type: 'number', default: 1.5, min: 0, max: 18, group: 'Flow', role: 'bc', help: 'Scour-forming condition: persistent tidal current with the usual sea state. Storm waves at low Keulegan–Carpenter number backfill rather than deepen the hole.' },
    { key: 'Tp', label: 'Peak period', unit: 's', type: 'number', default: 6.5, min: 3, max: 22, group: 'Flow', role: 'bc' },
    { key: 'HsStorm', label: 'Storm significant wave height (protection design)', unit: 'm', type: 'number', default: 6, min: 0, max: 18, group: 'Flow', role: 'bc', bind: 'site.Hs100', help: 'The armour must stay in place in the extreme sea state.' },
    { key: 'TpStorm', label: 'Storm peak period', unit: 's', type: 'number', default: 10.5, min: 3, max: 22, group: 'Flow', role: 'bc', bind: 'site.Tp100' },
    { key: 'hours', label: 'Duration of the condition', unit: 'h', type: 'number', default: 720, min: 0.1, max: 1e6, group: 'Flow', role: 'bc', help: 'Storm duration, or a long time for persistent tidal currents.' },
    { key: 'nSigma', label: 'Standard deviations on the scour depth', unit: '-', type: 'number', default: 0, min: 0, max: 3, group: 'Scour criteria', role: 'criteria', help: 'S/D = 1.3 + n·0.7 in steady current (Sumer & Fredsøe); DNV-ST-0126 quotes 1.3 D as the design value unless site data justify otherwise.' },
    { key: 'ampl', label: 'Bed shear amplification at the pile', unit: '-', type: 'number', default: 4, min: 1, max: 12, group: 'Scour criteria', role: 'criteria', help: 'About 4 at a circular pile in steady current, 2–3 in wave-dominated flow at low KC.' },
    { key: 'amplRock', label: 'Bed shear amplification for armour stability', unit: '-', type: 'number', default: 2.5, min: 1, max: 8, group: 'Scour criteria', role: 'criteria', help: 'Amplification on top of the armour layer in the storm; 2–3 in wave-dominated flow (Whitehouse 1998; DNV-RP-0618).' },
    { key: 'dfLim', label: 'Allowable natural-frequency reduction', unit: '%', type: 'number', default: 3, min: 0.5, max: 20, group: 'Scour criteria', role: 'criteria', help: 'Margin kept for scour in the frequency design; beyond it the 1P separation or the fatigue assessment must be revisited.' },
    { key: 'H', label: 'Horizontal load at the mudline', unit: 'MN', type: 'number', default: 4.9, min: 0.01, max: 100, group: 'Loads at the mudline', role: 'bc', bind: 'owt.mudlineShear' },
    { key: 'M', label: 'Overturning moment at the mudline', unit: 'MN·m', type: 'number', default: 280, min: 0, max: 5000, group: 'Loads at the mudline', role: 'bc', bind: 'owt.mudlineMoment' },
    { key: 'prot', label: 'Scour protection', type: 'select', default: 'rock', group: 'Protection', role: 'ic', options: [{ value: 'rock', label: 'Rock armour installed around the pile' }, { value: 'none', label: 'None — the pile is designed for scour' }],
      help: 'Initial seabed condition for the in-place analysis: with protection only the surveyed residual scour and the general lowering act on the pile.' },
    { key: 'd50Rock', label: 'Armour size provided d50', unit: 'mm', type: 'number', default: 500, min: 20, max: 2500, group: 'Protection', role: 'geometry', design: true, showIf: { key: 'prot', in: ['rock'] } },
    { key: 'rock', label: 'Armour rock', type: 'select', source: 'materials:rock', default: 'rock-granite', group: 'Protection', role: 'material' },
    { key: 'extent', label: 'Protection radius / pile diameter', unit: '-', type: 'number', default: 2.5, min: 1, max: 6, group: 'Protection', role: 'geometry', design: true, help: 'Outer radius of the armour layer from the pile axis; 2–3 D is common (DNV-RP-0618).' },
    { key: 'tFilter', label: 'Filter layer thickness', unit: 'm', type: 'number', default: 0.5, min: 0, max: 2, group: 'Protection', role: 'geometry' },
    { key: 'nElem', label: 'Pile elements', type: 'int', default: 40, min: 10, max: 160, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Pile elements', levels: [20, 40, 80], qoi: 'yIntact' },
  compute(inp, ctx) {
    const steel = steelOf(inp, ctx);
    const rock = ctx && ctx.material ? ctx.material(inp.rock) : null;
    const rhoR = rock?.rho ?? 2650;
    const d = inp.d50 / 1e3, s = inp.rhoS / inp.rhoW, h = inp.depth, D = inp.D;
    need(s > 1.05, 'Grain density must exceed the water density.');
    const z0 = d / 12; // grain roughness (Soulsby 1997)
    const T = inp.Tp, Uw = inp.Hs > 0 ? run(() => nearBedWaveCurrent({ H: inp.Hs, T, d: h, zBed: 0 }), 'Wave kinematics at the bed').Uw : 0;
    // Depth-averaged current consistent with a log profile through the 1 m reference value.
    const Ubar = inp.Uc * (Math.log(h / z0) - 1) / Math.log(1 / z0);
    const bs = run(() => bedShearStress({ d, h, U: Ubar, Uw, T, rho: inp.rhoW }), 'Bed shear stress');
    const cr = criticalShields({ d, s, nu: NU_SEA });
    const den = (inp.rhoS - inp.rhoW) * G * d;
    const theta = bs.tauMax / den, mobFar = theta / cr.thetaCr, mobNear = inp.ampl * mobFar;
    const scours = mobNear > 1;
    const KC = Uw > 0 ? keuleganCarpenter(Uw, T, D) : Infinity;
    const Ucw = inp.Uc + Uw > 0 ? inp.Uc / (inp.Uc + Uw) : 0;
    // The combined wave–current expression is calibrated for KC ≥ 4; applied to a large monopile with KC well below
    // that it predicts almost no scour even when the current dominates. Sumer & Fredsøe (2002, §3.5) note that for
    // Uc/(Uc+Um) above about 0.7 the scour approaches the steady-current value, which is applied here.
    const currentDominated = Ucw >= 0.7;
    const sf = pileScour({ KC: currentDominated ? Infinity : KC, Ucw: currentDominated ? 1 : Ucw, nSigma: inp.nSigma });
    let rb = null;
    if (!currentDominated && KC > 1 && KC < 10) { try { rb = pileScourRudolphBos({ KC, Ucw }); } catch { rb = null; } }
    const Sbreusers = scourBreusers({ D, h });
    // Design equilibrium depth: Sumer & Fredsøe, not below the Rudolph & Bos low-KC value where that applies.
    const Seq = scours ? Math.max(sf.SoverD, rb ? rb.SoverD : 0) * D : 0;
    let Tsc = Infinity;
    if (scours && theta > 0) {
      try { Tsc = (currentDominated || Ucw >= 0.5 || !Number.isFinite(KC) ? scourTimeScale({ mode: 'current', D, theta, d, s, delta: h }) : scourTimeScale({ mode: 'waves', D, theta, d, s, KC })).T; } catch { Tsc = Infinity; }
    }
    const t = inp.hours * 3600;
    const dev = (tt) => (Seq <= inp.S0 ? inp.S0 : inp.S0 + (Number.isFinite(Tsc) ? scourDevelopment(Seq - inp.S0, tt, Tsc) : 0));
    const St = dev(t), Sdes = Math.max(Seq, inp.S0);
    // Structural effect: pile response and stiffness without and with the design scour.
    const profile = soilProfile(inp);
    const opt0 = pileOptions({ ...inp, pyMode: 'api', loading: 'static' }, profile, steel);
    const Hn = inp.H * 1e6, Mn = inp.M * 1e6;
    const optS = Sdes + inp.Sg > 0.01 ? { ...opt0, profile: run(() => applyScour(profile, { general: inp.Sg, local: Sdes, D }), 'Scoured soil profile') } : opt0;
    need(inp.Lp - Sdes - inp.Sg > 1.5 * D, `With ${fmt(Sdes + inp.Sg)} m of scour only ${fmt(inp.Lp - Sdes - inp.Sg)} m of the pile remains embedded: protection or a longer pile is required before the response can be evaluated.`);
    const head = { H: Hn, M: Mn };
    const r0 = run(() => solveLateralPile({ ...opt0, head }), 'Pile response without scour'), rS = run(() => solveLateralPile({ ...optS, head }), 'Pile response with scour');
    const half = { H: 0.5 * Hn, M: 0.5 * Mn };
    const K0 = monopileStiffness({ ...opt0, head: half }).secant, KS = monopileStiffness({ ...optS, head: half }).secant;
    const tb = turbine(inp);
    const fOf = (K) => supportStructure({ ...inp, fnd: 'springs', Kh: K.KLL / 1e9, Kr: K.KRR / 1e9, Khr: K.KLR / 1e6 }, tb, steel, { nElem: 30, nModes: 1 }).modes[0].frequency;
    const f0 = fOf(K0), fS = fOf(KS);
    const dF = (1 - fS / f0) * 100;
    // Protected bed: only the surveyed residual scour and the general lowering remain.
    const protectedBed = inp.prot === 'rock';
    const Sres = inp.S0 + inp.Sg;
    let rR = r0, fR = f0;
    if (Sres > 0.01) {
      const optR = { ...opt0, profile: applyScour(profile, { general: inp.Sg, local: Math.max(inp.S0, 0.01), D }) };
      rR = run(() => solveLateralPile({ ...optR, head }), 'Pile response with residual scour');
      fR = fOf(monopileStiffness({ ...optR, head: half }).secant);
    }
    const dFres = (1 - fR / f0) * 100;
    // Protection: stable armour size, layer thickness, quantities and filter check against the seabed.
    let d50r = 0, rockNote = null;
    const UwS = inp.HsStorm > 0 ? nearBedWaveCurrent({ H: inp.HsStorm, T: inp.TpStorm, d: h, zBed: 0 }).Uw : 0;
    try { d50r = rockSizeShields({ h, U: Ubar, Uw: UwS, T: inp.TpStorm, rho: inp.rhoW, amplification: inp.amplRock, rhoS: rhoR, nu: NU_SEA }).d50; } catch (e) { rockNote = `No statically stable rock size was found for this loading (${e.message}): use a dynamically stable design or larger units.`; }
    // Layer of three stones of the size placed (the provided size when protection is installed, else the required size).
    const dArm = inp.prot === 'rock' ? Math.max(inp.d50Rock / 1e3, d50r) : d50r;
    const tArm = Math.max(3 * dArm, 0.5);
    const arm = armourTakeoff({ D, extentFactor: inp.extent, thickness: tArm, rhoS: rhoR });
    const fil = inp.tFilter > 0 ? armourTakeoff({ D, extentFactor: inp.extent + 0.5, thickness: inp.tFilter, rhoS: rhoR }) : { mass: 0, volume: 0, area: arm.area };
    const rockT = (arm.mass + fil.mass) / 1e3;
    // Filter rule between armour and seabed (one-layer check): retention D15f/D85b ≤ 5.
    const d15f = 0.6 * Math.max(d50r, 1e-4);
    let filt = null;
    try { filt = filterCriteria({ D15f: d15f, D50f: Math.max(d50r, 1e-4), D15b: 0.5 * d, D50b: d, D85b: inp.d85 / 1e3 }); } catch { filt = null; }
    const needFilter = filt ? !filt.retention.ok : true;
    const checksUnprotected = [
      check('Natural-frequency reduction by scour', dF, inp.dfLim, '%', 'DNV-ST-0126 (2021) §7.7.1 (scour to be included in the dynamic analysis)', 'Guards the separation between the natural frequency and the rotor frequency: scour lengthens the free pile and softens the foundation, moving the frequency towards 1P and raising fatigue loads.'),
      check('Mudline rotation increase by scour', Math.abs(rS.head.theta) / DEG, 0.5, 'deg', 'DNV-ST-0126 (2021) §7.6.2 (serviceability rotation)', 'Rotation of the pile at the original mudline level under the design load with the scour hole present. Scour that takes the rotation beyond the tilt tolerance requires protection or backfilling.'),
      check('Pile bending moment with scour against first yield', rS.maxMoment * 1.35 / 1e6, steel.fy / 1.1 * opt0.sec.W / 1e6, 'MN·m', 'DNV-ST-0126 (2021) §4.4; IEC 61400-3-1 §7.6', 'Scour moves the point of fixity down and raises the maximum moment in the pile. Exceedance means the embedded cans are too thin for the unprotected condition.'),
    ];
    const checksProtected = [
      check('Armour stone size: required over provided d50', d50r * 1e3, inp.d50Rock, 'mm', 'DNV-RP-0618 (2022) §5; CIRIA C683 §5.2.3 (critical-shear stability)', 'Guards the armour layer against being moved by the amplified wave and current shear stress in the design storm. Undersized stone is carried away and the pile scours as if unprotected.'),
      check('Natural-frequency reduction with the protected bed', dFres, inp.dfLim, '%', 'DNV-ST-0126 (2021) §7.7.1', 'Frequency loss from the surveyed residual scour and general seabed lowering that remain with protection in place.'),
      check('Mudline rotation with the protected bed', Math.abs(rR.head.theta) / DEG, 0.5, 'deg', 'DNV-ST-0126 (2021) §7.6.2 (serviceability rotation)', 'Rotation of the pile at the original mudline level under the design load with the residual scour: the protection must keep it within the tilt tolerance.'),
    ];
    const checks = protectedBed ? checksProtected : checksUnprotected;
    const gov = worst(checks), ok = gov.uc <= 1;
    const ts = [], ss = [];
    const tEnd = Number.isFinite(Tsc) ? Math.max(3 * Tsc, t) : t;
    for (let i = 0; i <= 40; i++) { const tt = (tEnd * i) / 40; ts.push(round(tt / 3600)); ss.push(dev(tt)); }
    // Frequency against scour depth.
    const sd = [0, 0.5, 1, 1.5, 2].map((q) => q * D).filter((q) => inp.Lp - q - inp.Sg > 1.6 * D);
    const fs = sd.map((q) => (q + inp.Sg < 0.01 ? f0 : fOf(monopileStiffness({ ...opt0, profile: applyScour(profile, { general: inp.Sg, local: Math.max(q, 0.01), D }), head: half }).secant)));
    const vdays = 0.6;
    return {
      kpis: [
        { key: 'Seq', label: 'Equilibrium local scour depth', unit: 'm', value: Seq, explain: `Depth the scour hole approaches if the condition persists (${fmt(Seq / D)} D). In steady current the mean value is 1.3 D; in waves at low Keulegan–Carpenter number it is much smaller.` },
        { key: 'St', label: 'Scour depth after the stated duration', unit: 'm', value: St, explain: 'Scour reached at the end of the condition, starting from the surveyed depth. Compare with the equilibrium value to see whether a single storm is long enough to develop it.' },
        { key: 'Tscour', label: 'Scour time scale', unit: 'h', value: Number.isFinite(Tsc) ? Tsc / 3600 : 0, explain: 'Time to reach 63 % of equilibrium. Hours for fine sand in strong flow; zero is shown when the bed is immobile.' },
        { key: 'thetaCr', label: 'Critical Shields parameter', unit: '-', value: cr.thetaCr, explain: 'Threshold of grain motion for the seabed sediment (Soulsby 1997).' },
        { key: 'mobility', label: 'Mobility ratio at the pile', unit: '-', value: mobNear, explain: 'Amplified Shields parameter over critical. Above 1 scour develops: live-bed when the far-field ratio also exceeds 1, clear-water otherwise.' },
        { key: 'KC', label: 'Keulegan–Carpenter number', unit: '-', value: Number.isFinite(KC) ? KC : 0, explain: 'Wave orbital excursion relative to the pile diameter (zero is shown without waves). Below about 6 waves alone cause little scour at a large monopile.' },
        { key: 'yHead', label: 'Mudline displacement with scour', unit: 'mm', value: rS.head.y * 1e3, explain: 'Lateral displacement at the original mudline level under the design loads with the scour hole present.' },
        { key: 'yIntact', label: 'Mudline displacement, intact seabed', unit: 'mm', value: r0.head.y * 1e3, explain: 'Lateral displacement at the mudline under the same loads before any scour: the reference for the increase caused by the scour hole.' },
        { key: 'f1', label: 'Natural frequency without scour', unit: 'Hz', value: f0, explain: 'First fore-aft frequency with the foundation stiffness of the intact seabed.' },
        { key: 'f1Scoured', label: 'Natural frequency with unprotected scour', unit: 'Hz', value: fS, explain: 'First frequency with the equilibrium scour hole and general lowering in the soil model, as it would be without protection. The drop moves the structure towards the rotor frequency.' },
        { key: 'dF', label: 'Frequency reduction if unprotected', unit: '%', value: dF, explain: 'Relative loss of natural frequency caused by the full scour hole; typically 2–8 % for unprotected monopiles in sand. It is the case for protection.' },
        { key: 'dFdesign', label: 'Frequency reduction in the design condition', unit: '%', value: protectedBed ? dFres : dF, limit: inp.dfLim, sense: 'max', explain: 'Frequency loss in the condition assessed: residual scour under the protection, or the full scour hole when the pile is left unprotected.' },
        { key: 'kLoss', label: 'Loss of rocking stiffness', unit: '%', value: (1 - KS.KRR / K0.KRR) * 100, explain: 'Reduction of the mudline rocking stiffness by the scour hole, the main cause of the frequency drop.' },
        { key: 'momentRise', label: 'Increase of maximum pile moment', unit: '%', value: (rS.maxMoment / r0.maxMoment - 1) * 100, explain: 'How much the peak bending moment in the pile grows because the fixity point moves down with scour.' },
        { key: 'rockD50', label: 'Stable armour size d50', unit: 'mm', value: d50r * 1e3, explain: 'Smallest median stone size that stays in place under the amplified wave–current bed shear stress (critical-shear criterion). Zero if no static solution exists.' },
        { key: 'rockT', label: 'Rock for scour protection', unit: 't', value: rockT, explain: 'Armour layer three stones thick (minimum 0.5 m) to the stated radius plus the filter layer, at 40 % porosity. Reported whether or not protection is selected, as the quantity it would take.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Scour development in time', x: { label: 'Time', unit: 'h' }, y: { label: 'Scour depth', unit: 'm' },
          series: [{ name: 'Local scour', x: ts, y: ss }, { name: 'Equilibrium', x: [0, ts[ts.length - 1]], y: [Sdes, Sdes], dash: true }, { name: 'End of condition', x: [inp.hours], y: [St] }],
          explain: 'Exponential growth of the scour hole towards equilibrium. If the curve is still steep at the end of the condition, longer events will scour deeper.' },
        { type: 'line', title: 'Natural frequency against scour depth', x: { label: 'Local scour depth', unit: 'm' }, y: { label: 'First natural frequency', unit: 'Hz' },
          series: [{ name: 'f1', x: sd, y: fs }, { name: 'Design scour', x: [Sdes], y: [fS] }, { name: 'Allowable reduction', x: [0, sd[sd.length - 1]], y: [f0 * (1 - inp.dfLim / 100), f0 * (1 - inp.dfLim / 100)], dash: true }],
          explain: 'Loss of natural frequency as the scour hole deepens. The crossing with the dashed line is the scour depth at which protection or intervention is needed.' },
        { type: 'line', title: 'Pile deflection with and without scour', x: { label: 'Deflection', unit: 'mm' }, y: { label: 'Depth below original mudline', unit: 'm' },
          series: [{ name: 'Intact seabed', x: Array.from(r0.y, (v) => v * 1e3), y: Array.from(r0.z, (v) => -v) }, { name: 'With scour', x: Array.from(rS.y, (v) => v * 1e3), y: Array.from(rS.z, (v) => -v), dash: true }],
          explain: 'Deflected shape of the pile under the same mudline loads. With scour the upper soil no longer supports the pile and the deflection and rotation grow.' },
      ],
      tables: [
        { title: 'Flow, sediment and scour estimates', columns: [{ key: 'k', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['Depth-averaged current', Ubar, 'm/s'], ['Wave orbital velocity at the bed', Uw, 'm/s'], ['Storm orbital velocity at the bed', UwS, 'm/s'], ['Maximum bed shear stress', bs.tauMax, 'Pa'], ['Far-field mobility θ/θcr', mobFar, '-'], ['Current ratio Uc/(Uc+Uw)', Ucw, '-'],
            ['Scour, Sumer & Fredsøe', sf.SoverD * D, 'm'], ['Scour, Rudolph & Bos (1 < KC < 10)', rb ? rb.SoverD * D : '–', 'm'], ['Scour, Breusers (current only)', Sbreusers, 'm'], ['DNV default 1.3 D', 1.3 * D, 'm'],
            ['Mudline displacement, intact', r0.head.y * 1e3, 'mm'], ['Mudline displacement, scoured', rS.head.y * 1e3, 'mm'], ['Mudline rotation, intact', Math.abs(r0.head.theta) / DEG, 'deg'], ['Mudline rotation, scoured', Math.abs(rS.head.theta) / DEG, 'deg'],
            ['K_RR intact', K0.KRR / 1e9, 'GN·m/rad'], ['K_RR scoured', KS.KRR / 1e9, 'GN·m/rad'], ['Armour layer thickness', tArm, 'm'], ['Armour mass', arm.mass / 1e3, 't'], ['Filter mass', fil.mass / 1e3, 't'], ['Protected seabed area', fil.area ?? arm.area, 'm²'], ['Filter retention D15f/D85b', filt ? filt.retention.ratio : '–', '-']]
            .map(([k, v, u]) => ({ k, v: typeof v === 'number' ? round(v) : v, u })),
          explain: 'Hydraulic inputs, alternative scour estimates and the structural consequences side by side, with the protection quantities. Compare the scour formulae before fixing the design value.' },
      ],
      notes: [scours ? (mobFar > 1 ? 'Live-bed conditions: the whole seabed is mobile.' : 'Clear-water scour: the bed moves only where the pile amplifies the flow.') : 'The seabed is not mobile even at the pile in this condition: no scour is predicted.',
        `${currentDominated ? 'Current-dominated flow: steady-current scour applies. ' : 'Wave-dominated flow: scour is limited by the low Keulegan–Carpenter number. '}Unprotected scour of ${fmt(Sdes)} m (${fmt(Sdes / D)} D) would lower the natural frequency from ${fmt(f0)} to ${fmt(fS)} Hz. ${inp.prot === 'rock' ? 'The bed is protected: the checks use the residual scour only.' : 'No protection: the checks use the full scour hole.'}`, needFilter ? 'The armour does not retain the seabed sand directly: a filter layer or geotextile is required.' : 'The armour grading retains the seabed material without a separate filter.', ...(rockNote ? [rockNote] : [])],
      advice: [
        protectedBed
          ? (ok ? adv('engineering', sev(gov.uc), `With rock protection the foundation keeps its stiffness: frequency loss ${fmt(dFres)} % (limit ${inp.dfLim} %). Without it the scour would reach ${fmt(Sdes)} m and the frequency would drop ${fmt(dF)} %.`,
            `Specify armour d50 ≥ ${fmt(d50r * 1e3)} mm (provided ${inp.d50Rock} mm), ${fmt(tArm)} m thick to ${fmt(inp.extent)} D radius${needFilter ? ` on a ${fmt(inp.tFilter)} m filter layer` : ''}; survey after the first winter for edge scour and falling apron movement.`,
            'DNV-RP-0618 (2022); DNV-ST-0126 (2021) §7.7; CIRIA C683 §5.2.3 and §5.4.3',
            `Protection uses ${fmt(rockT)} t of quarried rock over ${fmt(arm.area)} m² of seabed (it becomes hard-substrate habitat) and ${vesselCo2Text(vdays)} The alternative of ${fmt(Sdes)} m more pile is ${fmt(opt0.sec.A * steel.rho * Sdes / 1e3)} t of steel ≈ ${fmt(opt0.sec.A * steel.rho * Sdes / 1e3 * 1.9)} t CO₂e (screening 1.9 t/t).`)
            : adv('engineering', 'fail', `${gov.name}: ${fmt(gov.demand)} against ${fmt(gov.capacity)} ${gov.unit}.`,
              gov.name.startsWith('Armour') ? `Increase the armour to d50 ≥ ${fmt(d50r * 1e3)} mm (layer ${fmt(tArm)} m thick), or use a wider grading with a dynamically stable design verified by model tests.` : `Backfill the existing ${fmt(inp.S0)} m scour hole before placing the armour, or include the residual scour in the frequency and fatigue design.`,
              'DNV-RP-0618 (2022) §5; CIRIA C683 §5.2.3',
              `Larger stone raises the protection to about ${fmt(rockT)} t of rock on ${fmt(arm.area)} m² of seabed; ${vesselCo2Text(vdays)}`))
          : (ok ? adv('engineering', sev(gov.uc), `Unprotected design scour of ${fmt(Sdes)} m is tolerable: frequency drops ${fmt(dF)} % (limit ${inp.dfLim} %), rotation ${fmt(Math.abs(rS.head.theta) / DEG)}°, pile moment +${fmt((rS.maxMoment / r0.maxMoment - 1) * 100)} %.`,
            `Include ${fmt(Sdes)} m of scour in the fatigue and frequency design, and survey the seabed after the first winter and after major storms; trigger protection if scour exceeds ${fmt(sd.find((q, i) => fs[i] < f0 * (1 - inp.dfLim / 100)) ?? Sdes)} m.`,
            'DNV-ST-0126 (2021) §7.7; DNV-RP-0618 (2022); Sumer & Fredsøe (2002)',
            `Designing for scour avoids ${fmt(rockT)} t of quarried rock over ${fmt(arm.area)} m² of seabed and its placement: ${vesselCo2Text(vdays)}`)
            : adv('engineering', 'fail', `${gov.name}: ${fmt(gov.demand)} against ${fmt(gov.capacity)} ${gov.unit} with ${fmt(Sdes)} m of scour.`,
              `Install scour protection before or directly after piling: armour d50 ≥ ${fmt(d50r * 1e3)} mm, ${fmt(tArm)} m thick to ${fmt(inp.extent)} D radius on a ${fmt(inp.tFilter)} m filter (${fmt(rockT)} t), or lengthen the pile by about ${fmt(Sdes)} m.`,
              'DNV-RP-0618 (2022); DNV-ST-0126 (2021) §7.7; CIRIA C683 §5.2.3 and §5.4.3',
              `Protection needs ${fmt(rockT)} t of rock on ${fmt(arm.area)} m² of seabed (new hard-substrate habitat) and ${vesselCo2Text(vdays)} The alternative, ${fmt(Sdes)} m more pile, is ${fmt(opt0.sec.A * steel.rho * Sdes / 1e3)} t of steel ≈ ${fmt(opt0.sec.A * steel.rho * Sdes / 1e3 * 1.9)} t CO₂e (screening 1.9 t/t).`)),
        adv('economic', 'info', `Scour protection for this foundation is ${fmt(rockT)} t of rock and about ${fmt(vdays)} vessel day; the unprotected alternative is ${fmt(Sdes)} m of extra pile, ${fmt(opt0.sec.A * steel.rho * Sdes / 1e3)} t of steel.`,
          'Compare the two per foundation across the array: rock placement has a high mobilisation cost that is shared between positions, whereas extra pile length adds steel, driving time and noise at every position. Monitor-and-react is cheapest where the bed is only mobile in rare storms.',
          'Whole-life costing per ISO 15686-5; risk-based inspection and intervention planning per DNV-RP-G101'),
      ],
      publish: { 'quantities.rock_t': protectedBed ? rockT : 0 },
      quantities: protectedBed ? { rock_t: round(rockT), vessel_days: { rock_placement: vdays } } : {},
    };
  },
  verify: [
    { name: 'Equilibrium scour at a pile in steady current: S = 1.3·D', source: 'Sumer, Fredsøe & Christiansen (1992); DNV-ST-0126 (2021) §7.7: 1.3 × 6 m = 7.8 m', inputs: { Hs: 0, Uc: 1.2, D: 6, nSigma: 0, d50: 0.25 }, expect: [{ kpi: 'Seq', value: 7.8, tol: 1e-6 }] },
    { name: 'Critical Shields parameter of 0.2 mm sand: θcr = 0.30/(1+1.2D*) + 0.055(1 − e^(−0.02D*))', source: 'Soulsby (1997) Eq. 77 with D* = d[g(s−1)/ν²]^(1/3) = 4.121 (s = 2.65, ν = 1.36e-6 m²/s): θcr = 0.05046 + 0.00435 = 0.0548', inputs: { d50: 0.2, rhoS: 2650, rhoW: 1000 }, expect: [{ kpi: 'thetaCr', value: 0.05481, tol: 0.002 }] },
  ],
};
