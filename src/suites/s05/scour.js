// Scour and seabed erosion at a berth: bed mobility under tidal current, waves and propeller wash, local scour at
// the piles by several methods with its development in time, under-keel bed lowering by sediment continuity
// (Exner), the exposure of the pile foundation, and the rock needed to protect the bed.
import { keuleganCarpenter, pileScour, scourHec18, scourSheppardMelville, scourBreusers, scourTimeScale, scourDevelopment, criticalShields, bedShearStress, rockSizeShields, armourTakeoff } from '../../engine/geotech/scour.js';
import { waveOrbital, vanRijnSimplified, isbashStoneSize } from '../../engine/coastal/sediment.js';
import { propellerJet } from './_vessel.js';
import { G, need, run, r, f, check, eng, eco, linspace, CONFIRM, ECON_BASIS, CO2_ROCK, CO2_MGO, MGO_T_PER_DAY } from './_common.js';

export default {
  id: 's05-scour', suite: 5, title: 'Scour and seabed erosion at the berth',
  summary: 'Whether the bed moves under current, waves and propeller wash, how deep scour develops at the piles and under the keel and how fast, how much pile embedment is exposed, and the rock size and tonnage that protects the bed.',
  covers: [9], fidelity: 'reduced-order',
  standards: ['PIANC WG180 Guidelines for Protecting Berthing Structures from Scour Caused by Ships (2015)', 'Sumer & Fredsøe, The Mechanics of Scour in the Marine Environment (2002)', 'FHWA HEC-18, 5th ed. (2012) §7', 'CIRIA C683 The Rock Manual (2007) §5.2', 'Soulsby, Dynamics of Marine Sands (1997)', 'BS 6349-1-4:2021 (materials) and BS 6349-2:2019 §8'],
  models: ['Shields sediment mobility', 'Bed shear stress under current and waves (Soulsby)', 'Empirical pile-scour correlations (Sumer–Fredsøe, HEC-18, Sheppard–Melville, Breusers)', 'Exponential scour development in time', 'Propeller jet: efflux, Gaussian decay and bed velocity', 'Sediment continuity (Exner, conservative upwind finite volumes) with bed-load and suspended-load transport (van Rijn)', 'Isbash and Shields rock sizing'],
  limits: 'Empirical equilibrium scour at an isolated circular pile in non-cohesive sediment; pile-group and quay-wall effects, cohesive beds and rudder or twin-screw jet interaction are not modelled. The propeller jet is a free round jet with editable coefficients. Under-keel erosion is one-dimensional with the flow adapting instantly to the bed. Three-dimensional flow and morphodynamic modelling is outside the module; calibrate against repeat bathymetric surveys.',
  publishes: ['wharf.scourDepth', 'wharf.propBedVelocity', 'quantities.rock_t'],
  searchTerms: 'propeller scour berth pile scour protection PIANC WG180',
  inputs: [
    { key: 'h', label: 'Water depth at the berth', unit: 'm', type: 'number', default: 13, min: 2, max: 40, group: 'Site', role: 'bc', bind: 'site.depth' },
    { key: 'U', label: 'Depth-averaged tidal current', unit: 'm/s', type: 'number', default: 0.8, min: 0, max: 4, group: 'Site', role: 'bc', bind: 'site.Uc100', uncertain: { dist: 'normal', cov: 0.15 } },
    { key: 'Hs', label: 'Wave height at the berth', unit: 'm', type: 'number', default: 1.0, min: 0, max: 8, group: 'Site', role: 'bc' },
    { key: 'Tp', label: 'Wave period', unit: 's', type: 'number', default: 6, min: 2, max: 20, group: 'Site', role: 'bc' },
    { key: 'd50', label: 'Median grain size of the bed', unit: 'mm', type: 'number', default: 0.25, min: 0.06, max: 50, group: 'Sediment', role: 'soil', bind: 'site.d50', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'From grab samples or the borehole log. Below about 0.06 mm the bed is cohesive and these methods do not apply.' },
    { key: 'rhoS', label: 'Grain density', unit: 'kg/m³', type: 'number', default: 2650, min: 2000, max: 3200, group: 'Sediment', role: 'soil' },
    { key: 'rhoW', label: 'Water density', unit: 'kg/m³', type: 'number', default: 1025, min: 995, max: 1035, group: 'Sediment', role: 'environment', bind: 'site.rhoW' },
    { key: 'D', label: 'Pile diameter', unit: 'm', type: 'number', default: 1.067, min: 0.2, max: 4, group: 'Piles', role: 'geometry', bind: 'wharf.pileD' },
    { key: 'Lemb', label: 'Pile penetration below the design bed', unit: 'm', type: 'number', default: 28, min: 3, max: 90, group: 'Piles', role: 'geometry', bind: 'wharf.pileLemb' },
    { key: 'nPiles', label: 'Piles to protect', unit: '-', type: 'int', default: 40, min: 1, max: 2000, group: 'Piles', role: 'geometry', help: 'Piles of the seaward rows along the berth, where propeller wash and current act.' },
    { key: 'method', label: 'Design scour method', type: 'select', default: 'envelope', group: 'Piles', role: 'criteria', options: [{ value: 'envelope', label: 'Largest of Sumer–Fredsøe and Sheppard–Melville (tide, waves, propeller)' }, { value: 'sumer', label: 'Sumer–Fredsøe (waves and current)' }, { value: 'sheppard', label: 'Sheppard–Melville (current, clear-water to live-bed)' }, { value: 'hec18', label: 'HEC-18 (CSU) pier equation' }, { value: 'breusers', label: 'Breusers: 1.5 D tanh(h/D)' }] },
    { key: 'nSigma', label: 'Standard deviations added to the mean scour', unit: '-', type: 'number', default: 1, min: 0, max: 3, group: 'Piles', role: 'criteria', help: `The steady-current mean is 1.3 D with a standard deviation of 0.7 D; DNV-ST-0126 quotes 1.3 D as a default. ${CONFIRM}` },
    { key: 'general', label: 'General bed lowering', unit: 'm', type: 'number', default: 0.5, min: 0, max: 5, group: 'Piles', role: 'ic', help: 'Over-dredging tolerance and long-term bed level change from surveys, added to the local scour.' },
    { key: 'allow', label: 'Scour allowance in the pile design', unit: 'm', type: 'number', default: 3.0, min: 0, max: 15, group: 'Piles', role: 'criteria', design: true, help: 'Loss of bed level assumed in the strength and capacity checks of the piles.' },
    { key: 'elapsed', label: 'Duration of the flow event', unit: 'h', type: 'number', default: 72, min: 0.1, max: 1e5, group: 'Piles', role: 'ic', help: 'Cumulative time under the design current for the development curve.' },
    { key: 'P', label: 'Installed main-engine power', unit: 'kW', type: 'number', default: 20000, min: 100, max: 90000, group: 'Propeller wash', role: 'operation' },
    { key: 'fP', label: 'Share of power used at the berth', unit: '-', type: 'number', default: 0.15, min: 0.01, max: 1, group: 'Propeller wash', role: 'operation', uncertain: { dist: 'normal', cov: 0.3 }, help: 'About 0.05–0.15 for main propellers during berthing and unberthing; 1.0 for bow thrusters.' },
    { key: 'Dp', label: 'Propeller diameter', unit: 'm', type: 'number', default: 7, min: 0.5, max: 11, group: 'Propeller wash', role: 'geometry' },
    { key: 'hp', label: 'Height of the propeller axis above the bed', unit: 'm', type: 'number', default: 5.5, min: 0.5, max: 25, group: 'Propeller wash', role: 'bc', help: 'Under-keel clearance plus half the propeller diameter (plus the tip clearance to the keel line).' },
    { key: 'ducted', label: 'Ducted propeller or thruster', type: 'bool', default: false, group: 'Propeller wash', role: 'geometry' },
    { key: 'cV0', label: 'Efflux velocity coefficient', unit: '-', type: 'number', default: 1.15, min: 0.8, max: 1.8, group: 'Propeller wash', role: 'criteria', help: `V0 = c (P/(ρ Dp²))^(1/3). About 1.15 open and 1.17 ducted in PIANC WG180; ideal momentum theory gives 1.72 η^(1/3). ${CONFIRM}` },
    { key: 'cAxis', label: 'Jet axis decay coefficient', unit: '-', type: 'number', default: 2.8, min: 1.5, max: 6.2, group: 'Propeller wash', role: 'criteria', help: 'V_axis = A V0 D0/x. 2.8 for a propeller jet (Blaauw & van de Kaa); 6.2 for a low-turbulence free jet; lower with a rudder splitting the jet.' },
    { key: 'isbashC', label: 'Isbash coefficient for rock under the jet', unit: '-', type: 'number', default: 0.6, min: 0.4, max: 1.2, group: 'Scour protection', role: 'criteria', help: `D = U²/(2 g Δ C²). 0.86 exposed stone in ordinary turbulence, 1.2 embedded stone; propeller-jet guidance corresponds to about 0.55–0.7. ${CONFIRM}` },
    { key: 'rhoRock', label: 'Rock density', unit: 'kg/m³', type: 'number', default: 2650, min: 2200, max: 3200, group: 'Scour protection', role: 'material' },
    { key: 'amp', label: 'Shear amplification at the pile', unit: '-', type: 'number', default: 4, min: 1, max: 10, group: 'Scour protection', role: 'criteria', help: 'Bed shear stress next to the pile over the undisturbed value: about 4 in current, 2–4 in waves.' },
    { key: 'extent', label: 'Protection radius (multiple of the pile diameter)', unit: '-', type: 'number', default: 2.5, min: 1, max: 6, group: 'Scour protection', role: 'geometry' },
    { key: 'draft', label: 'Draft of the vessel alongside', unit: 'm', type: 'number', default: 11.5, min: 0, max: 25, group: 'Under-keel erosion', role: 'bc', help: 'Zero switches the vessel off: the bed then stays level (limiting case).' },
    { key: 'beam', label: 'Beam of the vessel', unit: 'm', type: 'number', default: 37, min: 5, max: 70, group: 'Under-keel erosion', role: 'geometry' },
    { key: 'keelAmp', label: 'Velocity under the keel over the cross-berth current', unit: '-', type: 'number', default: 1.5, min: 0.5, max: 8, group: 'Under-keel erosion', role: 'bc', help: 'Initial amplification of the current under the hull; about 1.2–2 from measurements and flow models, with depth / under-keel clearance as the two-dimensional upper bound when no flow can pass round bow and stern.' },
    { key: 'Ucross', label: 'Current component across the berth', unit: 'm/s', type: 'number', default: 0.4, min: 0, max: 3, group: 'Under-keel erosion', role: 'bc' },
    { key: 'days', label: 'Cumulative time alongside in this current', unit: 'days', type: 'number', default: 30, min: 0.1, max: 3650, group: 'Under-keel erosion', role: 'ic' },
    { key: 'nCell', label: 'Cells across the berth pocket', unit: '-', type: 'int', default: 160, min: 20, max: 1280, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nCell', label: 'Cells across the berth', levels: [40, 80, 160, 320], qoi: 'keelMean' },
  compute(inp) {
    const d = inp.d50 / 1e3, rho = inp.rhoW, s = inp.rhoS / rho, h = inp.h, D = inp.D;
    need(inp.rhoS > rho * 1.2, 'The grain density must be well above the water density.');
    const cs = criticalShields({ d, s }), tauCr = cs.tauCr * rho;
    const Uw = inp.Hs > 0 ? waveOrbital({ H: inp.Hs, T: inp.Tp, h }).Uw : 0;
    const bs = run(() => bedShearStress({ d, h, U: inp.U, Uw, T: inp.Tp, rho }), 'Bed shear stress');
    const theta = bs.tauMax / ((inp.rhoS - rho) * G * d);
    const Vc = Math.sqrt(tauCr / (rho * bs.Cd)); // depth-averaged current at the threshold of motion
    // propeller jet
    const jet = propellerJet({ P: inp.P * 1e3 * inp.fP, Dp: inp.Dp, hp: inp.hp, rho, ducted: !!inp.ducted, cV0: inp.cV0, cAxis: inp.cAxis });
    need(inp.hp > 0.5 * inp.Dp, 'The propeller axis must be more than half a propeller diameter above the bed.');
    const Vb = jet.VbMax;
    // local scour at a pile by each method
    const KC = Uw > 0 ? keuleganCarpenter(Uw, inp.Tp, D) : 0;
    const Ucw = inp.U + Uw > 0 ? inp.U / (inp.U + Uw) : 1;
    const live = theta > cs.thetaCr;
    const sf = live ? pileScour({ KC: Uw > 0 ? Math.max(KC, 1e-6) : Infinity, Ucw, nSigma: inp.nSigma }).SoverD * D : 0;
    const sm = (V) => (V > 0 ? scourSheppardMelville({ a: D, y1: h, V1: V, Vc, D50: d }).ys : 0);
    const smTide = sm(inp.U), smJet = sm(Vb);
    const hec = inp.U > 0 ? scourHec18({ a: D, y1: h, V1: inp.U }).ys : 0;
    const br = scourBreusers({ D, h });
    const methods = [{ name: 'Sumer–Fredsøe (waves + current, live bed)', S: sf }, { name: 'Sheppard–Melville (tidal current)', S: smTide }, { name: 'Sheppard–Melville (propeller bed velocity)', S: smJet }, { name: 'HEC-18 (tidal current)', S: hec }, { name: 'Breusers 1.5 D tanh(h/D)', S: br }];
    const Sloc = { envelope: Math.max(sf, smTide, smJet), sumer: sf, sheppard: Math.max(smTide, smJet), hec18: hec, breusers: br }[inp.method];
    const Stot = Sloc + inp.general;
    // development in time under the tidal current (live bed only)
    let Tsc = Infinity, St = 0;
    if (live && Sloc > 0) { Tsc = scourTimeScale({ mode: Uw > inp.U ? 'waves' : 'current', D, theta: Math.max(theta, 1e-6), d, s, delta: h, KC: Math.max(KC, 1e-6) }).T; St = scourDevelopment(Sloc, inp.elapsed * 3600, Tsc); }
    // under-keel bed lowering: Exner with van Rijn transport; the unit discharge under the keel is the cross-berth
    // current times the amplification entered times the initial clearance, and stays constant as the bed lowers
    const n = Math.round(inp.nCell), Lr = 5 * inp.beam, dx = Lr / n, q = inp.draft > 0 ? inp.keelAmp * inp.Ucross * (h - inp.draft) : inp.Ucross * h; // under the flat keel
    const xs = Array.from({ length: n }, (_, i) => (i + 0.5) * dx);
    // flat-bottomed hull section (the bilge radius is small against the beam of a cargo vessel)
    const Tk = xs.map((x) => (Math.abs(x - Lr / 2) <= inp.beam / 2 ? inp.draft : 0));
    need(inp.draft < h - 0.3, `The draft of ${f(inp.draft)} m leaves less than 0.3 m under the keel in ${f(h)} m of water.`);
    // unit discharge at each station, fixed in time: ambient current outside the hull, amplified under it
    const qx = Tk.map((t) => inp.Ucross * (1 + (inp.keelAmp - 1) * (inp.draft > 0 ? t / inp.draft : 0)) * (h - t));
    // the van Rijn transport relation is calibrated for sand of 0.1–2 mm: outside it the under-keel model is not run
    const sandy = d >= 1e-4 && d <= 2e-3;
    const qsOf = (zb, i) => { if (!sandy) return 0; const clr = Math.max(-zb - Tk[i], 0.2); const Ux = qx[i] / clr; return Ux > 0 ? vanRijnSimplified({ U: Ux, h: clr, d50: d, rhoS: inp.rhoS, rho }).qt : 0; };
    const zb0 = new Float64Array(n).fill(-h);
    let ex = { zb: zb0, mass: { relativeError: 0, error: 0, initial: 1 }, steps: 0 };
    const Ukeel0 = inp.draft > 0 ? q / (h - inp.draft) : inp.Ucross, qKeel0 = qsOf(-h, Math.floor(n / 2));
    let tSim = 0, atEq = false;
    if (q > 0 && inp.draft > 0 && qKeel0 > 0) {
      // Exner: (1 − p) ∂zb/∂t + ∂qs/∂x = 0 in conservation form. Transport is in +x and falls as the bed lowers
      // (larger clearance, lower velocity), so the bed celerity a = (∂qs/∂zb)/(1 − p) is positive everywhere and the
      // first-order upwind flux is the transport of the cell upstream. The step follows the largest celerity; the
      // run stops at the end time or once the bed has stopped changing (equilibrium).
      const p1 = 0.6, tEnd = inp.days * 86400, zb = Float64Array.from(zb0), qs = new Float64Array(n);
      let steps = 0, fluxNet = 0, moved = 0;
      while (steps < 6000 && tSim < tEnd - 1e-6 && !atEq) {
        let cel = 0;
        for (let i = 0; i < n; i++) { qs[i] = qsOf(zb[i], i); cel = Math.max(cel, Math.abs(qsOf(zb[i] - 0.01, i) - qs[i]) / 0.01 / p1); }
        const dtk = Math.min(0.5 * dx / Math.max(cel, 1e-12), tEnd - tSim, tEnd / 20);
        // upstream boundary: zero gradient (inflow equals the transport of the first cell)
        for (let i = n - 1; i >= 1; i--) { const dz = -(qs[i] - qs[i - 1]) * dtk / (p1 * dx); zb[i] += dz; moved = Math.max(moved, Math.abs(dz) / dtk); }
        fluxNet += (qs[0] - qs[n - 1]) * dtk; tSim += dtk; steps++;
        // equilibrium: the fastest bed-level change of the last 50 steps would move the bed by under 1 mm in the time left
        if (steps % 50 === 0) { atEq = moved * (tEnd - tSim) < 1e-3; moved = 0; }
      }
      let vol = 0; for (let i = 0; i < n; i++) vol += p1 * (zb[i] + h) * dx;
      ex = { zb, mass: { error: vol - fluxNet, initial: p1 * h * Lr }, steps };
    }
    let keel = 0, depo = 0, eroded = 0; for (let i = 0; i < n; i++) { keel = Math.max(keel, -h - ex.zb[i]); depo = Math.max(depo, ex.zb[i] + h); eroded += Math.max(-h - ex.zb[i], 0) * dx; }
    const keelMean = inp.draft > 0 ? eroded / inp.beam : 0; // eroded section area spread over the beam
    // equilibrium clearance: transport stops when the under-keel velocity falls to the threshold
    const UcrKeel = sandy ? vanRijnSimplified({ U: 0, h: Math.max(h - inp.draft, 0.2), d50: d, rhoS: inp.rhoS, rho }).Ucr : Infinity;
    const keelEq = inp.draft > 0 && q > 0 ? Math.max(q / UcrKeel - (h - inp.draft), 0) : 0; // first estimate: threshold velocity of the initial clearance
    const massErr = Math.abs(ex.mass.error) / Math.max(Math.abs(ex.mass.initial), 1e-9);
    // protection: rock against tide + waves at the pile (Shields) and against the jet (Isbash)
    const delta = (inp.rhoRock - rho) / rho;
    let dShields = 0;
    if (inp.U > 0 || Uw > 0) { try { dShields = rockSizeShields({ h, U: inp.U, Uw, T: inp.Tp, rho, amplification: inp.amp, rhoS: inp.rhoRock }).d50; } catch { throw new Error(`The current and wave loading at the pile is too severe for loose rock by the Shields approach (the stable stone would exceed one tenth of the ${f(h)} m depth): use concrete mattresses or grouted protection, or check the design wave at the berth.`); } }
    const dIsbash = isbashStoneSize({ U: Vb, delta, C: inp.isbashC });
    const dRock = Math.max(dShields, dIsbash, 0.02), tLayer = Math.max(2 * dRock, 0.3);
    const take = armourTakeoff({ D, extentFactor: inp.extent, thickness: tLayer, rhoS: inp.rhoRock });
    const rockT = take.mass * Math.round(inp.nPiles) / 1e3, days = Math.max(1, rockT / 1500);
    const thetaJet = bedShearStress({ d, h: Math.max(inp.hp, 1), U: Vb, rho }).tauC / ((inp.rhoS - rho) * G * d);
    const VbCrit = Vc * Math.sqrt(bs.Cd / bedShearStress({ d, h: Math.max(inp.hp, 1), U: 1, rho }).Cd);
    // power share at which the bed under the jet is at the threshold of motion (bed velocity ∝ P^(1/3))
    const fPerm = Math.min(1, inp.fP * (VbCrit / Math.max(Vb, 1e-9)) ** 3);
    const expo = Stot / inp.Lemb;
    const checks = [
      check('Total scour against the design allowance', Stot, inp.allow, 'm', 'BS 6349-2:2019 §8; API RP 2A-WSD §6.2 (scour in pile design)', 'General bed lowering plus local scour against the bed-level loss assumed in the pile checks. If exceeded the piles are longer and weaker in bending than designed and have less embedment.'),
      check('Under-keel bed lowering against the allowance', keel + inp.general, inp.allow, 'm', 'Exner sediment continuity; PIANC WG180 (2015) §7', 'Erosion beneath the moored vessel from the current squeezed under the keel, with general lowering. It undermines the front piles and any toe protection along the berth line.'),
    ];
    // Bed mobility (θ / θcr) is a regime indicator, not an acceptance criterion: a live bed is the normal state of a sandy
    // berth and its consequence is the scour assessed by the two checks above. It is reported as KPIs and advice.
    const structural = checks, w = structural.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = w.uc <= 1;
    const needProt = thetaJet > cs.thetaCr || !ok;
    const ts = Number.isFinite(Tsc) ? linspace(0, Math.max(5 * Tsc, inp.elapsed * 3600), 40) : [0, inp.elapsed * 3600];
    const jx = linspace(0.5 * jet.xMax, 4 * jet.xMax, 40);
    const stride = Math.max(1, Math.floor(n / 160));
    return {
      kpis: [
        { key: 'thetaCr', label: 'Critical Shields number', value: r(cs.thetaCr), unit: '-', explain: `Threshold of motion for ${f(inp.d50)} mm grains (dimensionless grain size ${f(cs.Dstar)}). The bed moves when the Shields number exceeds it.` },
        { key: 'theta', label: 'Shields number under tide and waves', value: r(theta), unit: '-', explain: 'Peak bed shear stress over the submerged grain weight per unit area. Its ratio to the critical value says whether the bed is at rest (clear water) or live.' },
        { key: 'mobAmbient', label: 'Bed mobility under tide and waves (θ / θcr)', value: r(theta / cs.thetaCr), unit: '-', explain: 'Ambient Shields number over its critical value (Soulsby 1997 Eq. 77). Above 1 the bed is live: scour holes reach equilibrium quickly but also refill. It sets the scour regime; it is not a pass/fail criterion.' },
        { key: 'mobJet', label: 'Bed mobility under the propeller jet (θ / θcr)', value: r(thetaJet / cs.thetaCr), unit: '-', explain: 'Shields number under the jet over its critical value (PIANC WG180 §8). Above 1 every manoeuvre moves bed material and a scour hole grows with the number of ship calls unless the bed is protected or power is limited.' },
        { key: 'Vc', label: 'Threshold depth-averaged current', value: r(Vc), unit: 'm/s', explain: 'Current speed at which this sediment starts to move at the berth depth.' },
        { key: 'V0', label: 'Propeller efflux velocity', value: r(jet.V0), unit: 'm/s', explain: 'Mean velocity of the jet just behind the propeller at the power share used, from the power and propeller diameter.' },
        { key: 'Vb', label: 'Largest bed velocity under the jet', value: r(Vb), unit: 'm/s', explain: `Peak near-bed velocity, reached about ${f(jet.xMax)} m behind the propeller. Compare with the threshold bed velocity of about ${f(VbCrit)} m/s for this sediment.` },
        { key: 'Sloc', label: 'Local scour depth at a pile (design)', value: r(Sloc), unit: 'm', explain: `Equilibrium scour hole depth by the selected method (${f(Sloc / D)} pile diameters).` },
        { key: 'Stot', label: 'Total scour (general + local)', value: r(Stot), unit: 'm', limit: inp.allow, sense: 'max', explain: 'Bed-level loss at the pile to carry into the pile strength, stiffness and capacity checks, against the allowance made in design.' },
        { key: 'St', label: 'Local scour after the event duration', value: r(St), unit: 'm', explain: `Depth reached after ${f(inp.elapsed)} h under the tidal current with the exponential development law (time scale ${Number.isFinite(Tsc) ? f(Tsc / 3600) + ' h' : 'not applicable: clear-water bed'}).` },
        { key: 'keelScour', label: 'Under-keel bed lowering', value: r(keel, 5), unit: 'm', explain: `Deepest erosion below the moored vessel from sediment continuity, ${atEq ? `at equilibrium (reached within ${f(tSim / 86400)} days)` : `after ${f(tSim / 86400)} days alongside`}. The uniform-flow estimate of the equilibrium is ${f(keelEq)} m.` },
        { key: 'keelMean', label: 'Mean bed lowering under the hull', value: r(keelMean, 5), unit: 'm', explain: 'Eroded cross-section area divided by the beam: the average loss of bed level beneath the vessel, and the measure of the volume that reappears as a bar on the downstream side.' },
        { key: 'massErr', label: 'Sediment mass balance error', value: massErr < 1e-9 ? 0 : r(massErr, 3), unit: '-', explain: 'Relative difference between the change of sediment volume in the model and the net flux through its ends. Erosion under the keel reappears as deposition downstream, so it should vanish: an error below 1e-9 of the sediment volume is arithmetic round-off and is reported as 0.' },
        { key: 'exposure', label: 'Pile embedment lost to scour', value: r(expo * 100), unit: '%', explain: 'Total scour as a share of the pile penetration. Above about 10 % re-check the axial and lateral pile capacity with the scoured profile.' },
        { key: 'dRock', label: 'Required rock size D50', value: r(dRock * 1e3), unit: 'mm', explain: `Larger of the Shields size at the pile under tide and waves (${f(dShields * 1e3)} mm) and the Isbash size under the propeller jet (${f(dIsbash * 1e3)} mm).` },
        { key: 'rock', label: 'Rock for pile scour protection', value: r(rockT), unit: 't', explain: `Rings of ${f(inp.extent)} D radius and ${f(tLayer)} m thickness around ${Math.round(inp.nPiles)} piles. A filter layer or geotextile is needed in addition.` },
        { key: 'fPerm', label: 'Maximum permissible power share at the berth', value: r(fPerm), unit: '-', explain: 'Share of installed power at which the unprotected bed under the jet is at the threshold of motion: an operating limit for manoeuvring without bed protection.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Local pile scour by method', categories: methods.map((m) => m.name), series: [{ name: 'Scour depth (m)', y: methods.map((m) => r(m.S)) }], x: { label: 'Method', unit: '' }, y: { label: 'Scour depth', unit: 'm' },
          explain: 'Equilibrium scour depth from each correlation for the same pile. The spread is the model uncertainty; the HEC-18 and Breusers values are bridge-pier envelopes and usually the highest.' },
        { type: 'line', title: 'Scour development in time', x: { label: 'Time', unit: 'h' }, y: { label: 'Scour depth', unit: 'm' },
          series: [{ name: 'S(t)', x: ts.map((t) => r(t / 3600)), y: ts.map((t) => r(Number.isFinite(Tsc) ? scourDevelopment(Sloc, t, Tsc) : 0)) }, { name: 'Allowance less general lowering', x: [0, r(ts[ts.length - 1] / 3600)], y: [r(Math.max(inp.allow - inp.general, 0)), r(Math.max(inp.allow - inp.general, 0))], dash: true }],
          explain: 'Growth of the scour hole under the design current. A short time scale means the hole is fully formed within one tide; a long one means surveys can catch it early.' },
        { type: 'line', title: 'Propeller jet velocity at the bed', x: { label: 'Distance behind the propeller', unit: 'm' }, y: { label: 'Velocity', unit: 'm/s' },
          series: [{ name: 'Bed velocity', x: jx.map((x) => r(x)), y: jx.map((x) => r(jet.at(x, inp.hp))) }, { name: 'Jet axis', x: jx.map((x) => r(x)), y: jx.map((x) => r(jet.axis(x))), dash: true }, { name: 'Threshold for the bed', x: [r(jx[0]), r(jx[jx.length - 1])], y: [r(VbCrit), r(VbCrit)], dash: true }],
          explain: 'Velocity along the bed behind the propeller. The length over which the solid line exceeds the threshold is the length of bed that needs protection.' },
        { type: 'line', title: 'Bed profile under the moored vessel', x: { label: 'Distance across the berth', unit: 'm' }, y: { label: 'Level', unit: 'm' },
          series: [{ name: 'Bed after erosion', x: xs.filter((_, i) => i % stride === 0).map((x) => r(x)), y: Array.from(ex.zb).filter((_, i) => i % stride === 0).map((z) => r(z, 5)) }, { name: 'Initial bed', x: [0, r(Lr)], y: [-h, -h], dash: true }, { name: 'Keel', x: xs.filter((_, i) => i % stride === 0).map((x) => r(x)), y: Tk.filter((_, i) => i % stride === 0).map((t) => r(-t)), dash: true }],
          explain: 'Cross-berth section with the hull and the eroded bed. Erosion concentrates under the keel where the flow is squeezed and the material settles in the wider section downstream.' },
      ],
      tables: [
        { title: 'Scour estimates', columns: [{ key: 'name', label: 'Method' }, { key: 'S', label: 'Scour depth', unit: 'm' }, { key: 'SD', label: 'S / D' }],
          rows: methods.map((m) => ({ name: m.name, S: r(m.S), SD: r(m.S / D, 3) })), explain: 'Each correlation with its result. The design value follows the method selected; report the range alongside it.' },
        { title: 'Flow and mobility', columns: [{ key: 'item', label: 'Quantity' }, { key: 'value', label: 'Value' }, { key: 'unit', label: 'Unit' }],
          rows: [{ item: 'Wave orbital velocity at the bed', value: r(Uw), unit: 'm/s' }, { item: 'Keulegan–Carpenter number', value: r(KC), unit: '-' }, { item: 'Bed shear stress, current', value: r(bs.tauC), unit: 'Pa' }, { item: 'Bed shear stress, waves', value: r(bs.tauW), unit: 'Pa' }, { item: 'Peak combined', value: r(bs.tauMax), unit: 'Pa' }, { item: 'Critical', value: r(tauCr), unit: 'Pa' }, { item: 'Under-keel velocity (initial)', value: r(inp.draft > 0 ? Ukeel0 : 0), unit: 'm/s' }, { item: 'Jet thrust', value: r(jet.thrust / 1e3), unit: 'kN' }],
          explain: 'Intermediate hydraulic quantities, for checking by hand and against measurements at the berth.' },
      ],
      notes: [`The ambient bed is ${live ? 'live (sediment in motion)' : 'at rest (clear-water conditions)'}; under the propeller jet the Shields number is ${f(thetaJet / cs.thetaCr)} times critical.`,
        `${sandy ? '' : `The bed-load and suspended-load relation covers sand of 0.1–2 mm, so under-keel erosion was not computed for ${f(inp.d50)} mm material. `}Design local scour ${f(Sloc)} m by ${inp.method === 'envelope' ? 'the envelope of methods' : inp.method}; under-keel lowering ${f(keel)} m after ${f(tSim / 86400)} days${atEq ? ' (equilibrium)' : tSim < inp.days * 86400 - 1 ? ` (the run stopped there: the step budget was used before the ${f(inp.days)} days entered)` : ''}, ${ex.steps} time steps, deposition up to ${f(depo)} m downstream.`],
      advice: [
        ok ? eng(needProt ? 'warn' : 'ok', `Total scour ${f(Stot)} m against an allowance of ${f(inp.allow)} m (${f(expo * 100)} % of the embedment); bed velocity under the propeller ${f(Vb)} m/s against a threshold of ${f(VbCrit)} m/s (Shields number ${f(thetaJet / cs.thetaCr)} times critical under the jet, ${f(theta / cs.thetaCr)} times under tide and waves).`,
          needProt ? `The jet moves the bed: protect the berth pocket with rock of D50 ≥ ${f(dRock * 1e3)} mm on a filter (or a mattress), or limit manoeuvring power to ${f(fPerm * 100)} % of installed power. Survey the bed after the first season.` : 'No protection is needed for the conditions entered. Confirm with a bathymetric survey after the first year and after any change of vessel class.',
          'PIANC WG180 (2015); Sumer & Fredsøe (2002); CIRIA C683 (2007) §5.2',
          `Protection for ${Math.round(inp.nPiles)} piles is ${f(rockT)} t of rock ≈ ${f(rockT * CO2_ROCK)} t CO₂e (screening factor ${CO2_ROCK} t/t) and covers ${f(take.area * Math.round(inp.nPiles))} m² of seabed habitat.`)
          : eng('fail', `${w.name}: ${f(w.demand)} m against ${f(w.capacity)} m.`,
            `Either re-check the piles with ${f(Math.max(Stot, keel + inp.general))} m of scour (strength, sway and capacity modules), or place scour protection: rock of D50 ≥ ${f(dRock * 1e3)} mm, ${f(tLayer)} m thick, to ${f(inp.extent)} D around each pile and across the propeller track.`, 'PIANC WG180 (2015) §8–9; CIRIA C683 (2007) §5.2.3',
            `Rock protection: ${f(rockT)} t ≈ ${f(rockT * CO2_ROCK)} t CO₂e (screening factor) plus about ${f(days)} placement day(s) ≈ ${f(days * MGO_T_PER_DAY * CO2_MGO)} t CO₂ of fuel. Longer piles instead would add steel on every pile.`),
        eco('info', `Protecting the berth needs about ${f(rockT)} t of rock and ${f(days)} marine-spread day(s); an operating limit of ${f(fPerm * 100)} % power avoids it but slows manoeuvres.`,
          'Compare three options over the life of the berth: design the piles for the scour, protect the bed once, or monitor with repeat surveys and intervene when the hole reaches the allowance. Monitoring is cheapest where the bed is at rest between ship calls.', `${ECON_BASIS}; risk-based inspection planning (DNV-RP-G101 principles)`),
      ],
      publish: { 'wharf.scourDepth': Stot, 'wharf.propBedVelocity': Vb, 'quantities.rock_t': rockT },
      quantities: { rock_t: r(rockT), vessel_days: { rock_placement: r(days, 3) } },
    };
  },
  verify: [
    // Soulsby (1997) Eq. 77 for d = 0.2 mm in sea water (s = 2650/1025 = 2.585, ν = 1.36e-6 m²/s):
    // D* = d·[g(s − 1)/ν²]^(1/3) = 2e-4·(15.548/1.8496e-12)^(1/3) = 4.067; θcr = 0.30/(1 + 1.2·4.067) + 0.055·(1 − e^(−0.0813)) = 0.0553.
    { name: 'Critical Shields number of fine sand', source: 'Soulsby (1997) Dynamics of Marine Sands, Eq. 77 (hand calculation)', inputs: { d50: 0.2 }, expect: [{ kpi: 'thetaCr', value: 0.0553, tol: 0.005 }] },
    // Strong steady current on a live bed, no waves: equilibrium scour 1.3 D (mean, no added standard deviation) = 1.3·1.0 = 1.30 m.
    { name: 'Live-bed scour at a pile in steady current: S = 1.3 D', source: 'Sumer, Fredsøe & Christiansen (1992); Sumer & Fredsøe (2002) §3.2', inputs: { D: 1, U: 1.5, Hs: 0, nSigma: 0, method: 'sumer', general: 0 }, expect: [{ kpi: 'Sloc', value: 1.3, tol: 1e-6 }] },
    // Propeller jet: P = 5 MW on a 5 m open propeller, axis 6 m above the bed:
    // V0 = 1.15·(5e6/(1025·25))^(1/3) = 1.15·5.800 = 6.670 m/s; D0 = 5/√2 = 3.536 m;
    // Vb = 2.8·e^(−½)/√(2·15.4)·V0·D0/hp = 0.3060·6.670·3.536/6 = 1.203 m/s;
    // Isbash with Δ = (2650 − 1025)/1025 = 1.585, C = 0.6: D = 1.203²/(2·9.807·1.585·0.36) = 0.1293 m.
    { name: 'Propeller efflux, bed velocity and Isbash stone size', source: 'Blaauw & van de Kaa (1978) jet relations as used in PIANC WG180; Isbash (1936) — hand calculation',
      inputs: { P: 5000, fP: 1, Dp: 5, hp: 6, ducted: false, cV0: 1.15, cAxis: 2.8, isbashC: 0.6, U: 0, Hs: 0 }, expect: [{ kpi: 'V0', value: 6.670, tol: 0.002 }, { kpi: 'Vb', value: 1.203, tol: 0.003 }, { kpi: 'dRock', value: 129.3, tol: 0.01 }] },
    // No vessel alongside: uniform flow over a level bed transports sediment without gradients, so the bed does not change.
    { name: 'Limiting case: no vessel, no under-keel erosion', source: 'Exner equation with ∂qs/∂x = 0 (uniform flow)', inputs: { draft: 0 }, expect: [{ kpi: 'keelScour', value: 0, abs: 1e-9 }, { kpi: 'massErr', value: 0, abs: 1e-9 }] },
  ],
};
