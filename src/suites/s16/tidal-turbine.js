// Tidal-stream turbine hydrodynamics and performance: optimum-rotor BEM (or actuator disc) performance, power and
// thrust curves, harmonic tide and annual energy, channel blockage, blade-tip cavitation, extreme thrust and the
// loads on the support column.
import { designOptimalRotor } from '../../engine/aero/rotors.js';
import { genericPolar, cpMinEstimate } from '../../engine/aero/polar.js';
import { bemSolve, BETZ_LIMIT } from '../../engine/aero/bem.js';
import { springNeap, synthesizeTide, velocityDistribution, currentProfile, rotorAveragedVelocity, actuatorRotor, tidalPowerCurve, cavitationCheck, channelBlockage, channelBlockageLimit,
  tidalEnergy, extremeThrust, supportLoads, kineticPowerDensity } from '../../engine/aero/tidal.js';
import { tubularSectionStress } from '../../engine/aero/blade.js';
import { waterVapourPressure } from '../../engine/flow/water.js';
import { site, adv, check, worst, sev, fmt, round, run, need, linspace, thin, G, DEG, LIMIT_TIDAL, LIMIT_POLAR, LIMIT_CODE, CFD_OUT, CONFIRM, steelCo2Text, energyCo2Text, CO2_STEEL } from './_common.js';

const bem = { key: 'rotorModel', in: ['bem'] };

/** Channel state with the same disc-referred thrust coefficient as the open-water rotor (equal disc operating point, Barnsley & Wellicome argument inverted). */
function blockedState(ctOpen, B, Fr) {
  const a = 0.5 * (1 - Math.sqrt(Math.max(1 - Math.min(ctOpen, 0.96), 0))), target = ctOpen / (1 - a) ** 2;
  if (B < 1e-4) return { gainCp: 1, gainCt: 1, alpha2: 1 - a };
  let lo = 0.02, hi = 0.995, st = null;
  const f = (x) => { try { return channelBlockage({ B, alpha4: x, Fr }); } catch { return null; } };
  for (let i = 0; i < 60; i++) {
    const mid = 0.5 * (lo + hi), s = f(mid);
    if (!s) { lo = mid; continue; }     // bypass flow would be critical at this loading: move to a lighter-loaded disc
    st = s;
    if (s.CtLocal > target) lo = mid; else hi = mid;
  }
  if (!st) throw new Error('The channel is so blocked that the bypass flow becomes critical: reduce the rotor area or widen the channel section.');
  const r = st.alpha2 / (1 - a);
  return { gainCp: r ** 3, gainCt: r ** 2, alpha2: st.alpha2 };
}

export default {
  id: 's16-tidal-turbine', suite: 16, title: 'Tidal-stream turbine performance and loading',
  summary: 'Power and thrust of a horizontal-axis tidal turbine, annual energy from the harmonic tide, channel blockage, tip cavitation margin and support-column loads.',
  covers: [3], fidelity: 'numerical',
  standards: ['IEC TS 62600-200:2013 (tidal energy converter power performance)', 'IEC TS 62600-201:2015 (tidal resource assessment)', 'IEC TS 62600-2:2019 (marine energy design requirements)', 'DNV-ST-0164 (2015) Tidal turbines', 'DNV-RP-C205 (2021) §6 (Morison loads)'],
  models: ['Blade element momentum theory adapted to tidal turbines', 'Actuator-disc hydrodynamics and Betz limit', 'Glauert–Schmitz optimum rotor design', 'Prandtl tip loss and Buhl high-thrust correction', 'Harmonic tide (M2 + S2 spring–neap)', 'Power-law current profile', 'One-dimensional channel blockage (Garrett & Cummins 2007; Houlsby et al. 2008)', 'Cavitation number against minimum pressure coefficient', 'Morison equation with linear (Airy) wave kinematics'],
  limits: `${LIMIT_TIDAL} The rotor is an optimum planform with one generic parametric polar: ${LIMIT_POLAR} Wave loading is quasi-steady with linear waves (no dynamic inflow or added mass of the rotor). ${CFD_OUT} ${LIMIT_CODE}`,
  publishes: ['tidal.ratedPower', 'tidal.aep', 'tidal.thrustMax', 'tidal.baseMoment', 'tidal.capacityFactor'],
  searchTerms: 'tidal stream turbine blade element momentum cavitation blockage Garrett Cummins',
  inputs: [
    { key: 'D', label: 'Rotor diameter', unit: 'm', type: 'number', default: 20, min: 2, max: 40, group: 'Rotor', role: 'geometry', design: true },
    { key: 'blades', label: 'Blades', type: 'int', default: 3, min: 2, max: 5, group: 'Rotor', role: 'geometry' },
    { key: 'ratedMW', label: 'Rated electrical power', unit: 'MW', type: 'number', default: 1.5, min: 0.05, max: 5, group: 'Rotor', role: 'geometry', design: true },
    { key: 'rotorModel', label: 'Rotor model', type: 'select', default: 'bem', group: 'Rotor', role: 'geometry', options: [{ value: 'bem', label: 'Optimum rotor, blade element momentum' }, { value: 'disc', label: 'Ideal actuator disc' }] },
    { key: 'tsr', label: 'Design tip-speed ratio', unit: '-', type: 'number', default: 4.5, min: 2, max: 8, group: 'Rotor', role: 'ic', showIf: bem, design: true, help: 'Rotor speed follows Ω = λU/R below rated. Tidal rotors run at 3.5–5 to keep tip speed below the cavitation limit.' },
    { key: 'alphaDes', label: 'Design angle of attack', unit: 'deg', type: 'number', default: 5, min: 1, max: 9, group: 'Rotor', role: 'geometry', showIf: bem, help: 'Sets the design lift coefficient on the generic 21 % thick polar.' },
    { key: 'induction', label: 'Axial induction factor', unit: '-', type: 'number', default: 1 / 3, min: 0.05, max: 0.45, group: 'Rotor', role: 'geometry', showIf: { key: 'rotorModel', in: ['disc'] }, help: '1/3 is the Lanchester–Betz optimum.' },
    { key: 'eff', label: 'Water-to-wire efficiency (drivetrain and electrical)', unit: '-', type: 'number', default: 0.88, min: 0.5, max: 1, group: 'Rotor', role: 'operation', uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'cutIn', label: 'Cut-in current speed', unit: 'm/s', type: 'number', default: 0.8, min: 0.2, max: 2, group: 'Rotor', role: 'control' },
    { key: 'depth', label: 'Water depth (lowest tide)', unit: 'm', type: 'number', default: 40, min: 8, max: 120, group: 'Site', role: 'environment', bind: 'site.depth' },
    { key: 'hubZ', label: 'Hub height above the seabed', unit: 'm', type: 'number', default: 18, min: 3, max: 80, group: 'Site', role: 'geometry', design: true },
    { key: 'width', label: 'Channel width', unit: 'm', type: 'number', default: 500, min: 20, max: 50000, group: 'Site', role: 'bc', help: 'Lateral extent of the flow section per device (lateral boundary). Blockage ratio = rotor area / (width × depth); a wide section leaves the flow unconstrained. The blockage effect is evaluated with one-dimensional actuator-disc theory.' },
    { key: 'springPeak', label: 'Mean spring peak current (surface)', unit: 'm/s', type: 'number', default: 3.6, min: 0.5, max: 6, group: 'Tidal resource', role: 'bc', bind: 'site.Uc1', uncertain: { dist: 'normal', cov: 0.07 }, help: 'Inflow boundary: peak along-stream speed at mean springs, at the free surface.' },
    { key: 'neapRatio', label: 'Neap peak / spring peak', unit: '-', type: 'number', default: 0.5, min: 0.2, max: 1, group: 'Tidal resource', role: 'bc', help: '1 = no spring–neap modulation (M2 only).' },
    { key: 'profileN', label: 'Current-profile power-law denominator n', unit: '-', type: 'number', default: 7, min: 4, max: 12, group: 'Tidal resource', role: 'bc', help: 'u(z) = u_surface·(z/h)^(1/n); 1/7 is the usual open-channel profile.' },
    { key: 'surge', label: 'Extreme current (spring peak plus surge)', unit: 'm/s', type: 'number', default: 4.2, min: 0.5, max: 7, group: 'Extreme conditions', role: 'bc', bind: 'site.Uc100' },
    site('Hs100', { group: 'Extreme conditions', default: 6 }), site('Tp100', { group: 'Extreme conditions', default: 11 }),
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 995, max: 1035, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'Tsea', label: 'Sea temperature', unit: '°C', type: 'number', default: 12, min: 0, max: 32, group: 'Site', role: 'environment', bind: 'site.Tsea', help: 'Sets the vapour pressure for the cavitation check.' },
    { key: 'colD', label: 'Support-column diameter', unit: 'm', type: 'number', default: 3, min: 0.5, max: 8, group: 'Support column', role: 'geometry', design: true },
    { key: 'colT', label: 'Support-column wall thickness', unit: 'mm', type: 'number', default: 40, min: 10, max: 120, group: 'Support column', role: 'geometry', design: true },
    { key: 'steel', label: 'Column steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Support column', role: 'material' },
    { key: 'gammaF', label: 'Load factor on extreme environmental loads', unit: '-', type: 'number', default: 1.35, min: 1, max: 2, group: 'Acceptance', role: 'criteria', help: `IEC TS 62600-2 / DNV-ST-0164 normal design situation. ${CONFIRM}` },
    { key: 'gammaM', label: 'Material factor', unit: '-', type: 'number', default: 1.1, min: 1, max: 1.5, group: 'Acceptance', role: 'criteria', help: CONFIRM },
    { key: 'cavMargin', label: 'Required cavitation margin σ / (−Cp,min)', unit: '-', type: 'number', default: 1.2, min: 1, max: 2, group: 'Acceptance', role: 'criteria' },
    { key: 'nUnits', label: 'Number of turbines', type: 'int', default: 4, min: 1, max: 200, group: 'Array', role: 'geometry' },
    { key: 'avail', label: 'Technical availability', unit: '%', type: 'number', default: 90, min: 40, max: 100, group: 'Array', role: 'ic', help: 'Starting availability state for the energy estimate; subsea access makes 85–95 % typical of early arrays.' },
    { key: 'nStations', label: 'Blade stations', type: 'int', default: 20, min: 6, max: 80, group: 'Solver', role: 'control', showIf: bem },
  ],
  mesh: { key: 'nStations', label: 'Blade stations', levels: [10, 20, 40], qoi: 'cpRotor' },
  compute(inp, ctx) {
    const R = inp.D / 2, rho = inp.rhoW, A = Math.PI * R * R, P = inp.ratedMW * 1e6;
    need(inp.hubZ - R > 1, `The rotor (radius ${fmt(R)} m) hits the seabed at a hub height of ${inp.hubZ} m: leave at least 1 m (usually 25 % of the depth) below the blade tip.`);
    need(inp.depth - inp.hubZ - R > 1, `The blade tip breaks the surface: hub ${inp.hubZ} m above the bed with radius ${fmt(R)} m in ${inp.depth} m of water. Lower the hub or reduce the diameter.`);
    need(inp.surge >= inp.springPeak * 0.5, 'The extreme current is less than half the spring peak: check the two inputs.');
    // rotor performance at the design tip-speed ratio (constant below rated)
    let cpH, ct, sol = null, rotor = null, clTip = 1.0, tipThick = 0.21, aTip = inp.induction;
    if (inp.rotorModel === 'bem') {
      const polar = genericPolar({ thickness: 0.21, aspectRatio: 8 });
      const clDesign = polar.clAt(inp.alphaDes * DEG);
      need(clDesign > 0.2, 'The design angle of attack gives too little lift on the generic polar.');
      rotor = run(() => designOptimalRotor({ tipRadius: R, hubRadius: 0.12 * R, blades: inp.blades, tsr: inp.tsr, clDesign, alphaDesignDeg: inp.alphaDes, polar, nStations: inp.nStations, maxChord: 0.14 * R, name: 'Optimum tidal rotor' }), 'Rotor design');
      const Uref = 2;
      sol = run(() => bemSolve({ rotor, U: Uref, omega: inp.tsr * Uref / R, pitchDeg: 0, rho, mu: 1.4e-3 }), 'BEM solution');
      cpH = sol.Cp; ct = sol.Ct;
      const tipSt = sol.spanwise[sol.spanwise.length - 2] || sol.spanwise[sol.spanwise.length - 1];
      clTip = Math.max(tipSt.cl, 0.1); aTip = Math.min(Math.max(tipSt.a, 0), 0.5);
    } else {
      const d = actuatorRotor({ a: inp.induction, efficiency: 1 })(1);
      cpH = d.Cp; ct = d.Ct;
    }
    need(cpH > 0.05, 'The rotor produces no useful power at the design point: change the tip-speed ratio or design incidence.');
    // channel blockage: same disc operating point in a bounded section
    const B = A / (inp.width * inp.depth);
    need(B < 0.5, `Blockage ratio ${fmt(B)} is beyond the validity of the one-dimensional model (keep rotor area below half the section).`);
    const hubProfile = currentProfile({ model: 'power', uRef: 1, zRef: inp.depth, n: inp.profileN });
    const avg = rotorAveragedVelocity({ profile: hubProfile, hubHeight: inp.hubZ, D: inp.D });
    const Fr = inp.springPeak * avg.mean / Math.sqrt(G * inp.depth);
    const blk = run(() => blockedState(ct, B, Fr), 'Channel blockage');
    const lim = run(() => channelBlockageLimit({ B, Fr }), 'Blockage limit');
    const cpEl = cpH * blk.gainCp * inp.eff, ctEff = Math.min(ct * blk.gainCt, 1.6);
    const curve = run(() => tidalPowerCurve({ D: inp.D, ratedPower: P, rotor: () => ({ Cp: cpEl, Ct: ctEff }), rho, cutIn: inp.cutIn, uMax: Math.max(5, inp.surge * 1.15), step: 0.1, efficiency: inp.eff }), 'Power curve');
    const iPeak = curve.thrust.indexOf(Math.max(...curve.thrust)), thrustPeak = curve.thrust[iPeak];
    const ratedSpeed = Number.isFinite(curve.ratedSpeed) ? curve.ratedSpeed : 0;
    // tide: M2 + S2 at the surface, converted to the rotor-averaged (power-weighted) speed
    const sn = springNeap({ springPeak: inp.springPeak, neapPeak: inp.springPeak * inp.neapRatio });
    const tide = synthesizeTide({ constituents: sn.constituents, duration: 365 * 86400, dt: 600 });
    const dist = velocityDistribution(tide.u, { binWidth: 0.1 });
    const en = tidalEnergy({ curve, speeds: tide.u, availability: inp.avail / 100, velocityFactor: avg.powerVelocity, ratedPower: P });
    const hubU = new Float64Array(tide.u.length); for (let i = 0; i < hubU.length; i++) hubU[i] = tide.u[i] * avg.hub;
    const kpd = kineticPowerDensity(hubU, rho);
    // cavitation at top dead centre at rated speed (highest tip speed), instantaneous surface at the wave trough of the 1-in-n sea is not included
    const Ucav = ratedSpeed > 0 ? ratedSpeed : inp.springPeak * avg.powerVelocity;
    const tipSpeed = inp.rotorModel === 'bem' ? inp.tsr * Ucav : 4.5 * Ucav;
    const pv = waterVapourPressure({ T: inp.Tsea + 273.15, S: 35 });
    const cpMin = cpMinEstimate({ cl: clTip, thickness: tipThick, kT: 1.2 });
    const cav = cavitationCheck({ U: Ucav, tipSpeed, D: inp.D, hubDepth: inp.depth - inp.hubZ, cpMin, rho, pv, a: aTip, margin: inp.cavMargin });
    // extreme thrust: operating at the extreme current with the extreme wave orbital velocity, and the column loads
    const ctExt = curve.ct[Math.min(curve.ct.length - 1, Math.max(0, Math.round((inp.surge * avg.mean - curve.u[0]) / 0.1)))] || ctEff;
    const Hmax = 1.86 * inp.Hs100;
    const wave = { H: Math.min(Hmax, 0.78 * inp.depth), T: inp.Tp100, depth: inp.depth };
    const ext = extremeThrust({ D: inp.D, Ct: Math.max(ctExt, 0.05), springPeak: inp.surge * avg.mean, wave, hubHeight: inp.hubZ, rho });
    const extThrust = Math.max(ext.thrust, thrustPeak);
    const sup = supportLoads({ thrust: extThrust, hubHeight: inp.hubZ, columnDiameter: inp.colD, profile: currentProfile({ model: 'power', uRef: inp.surge, zRef: inp.depth, n: inp.profileN }), Cd: 1.0, Cm: 2.0, wave, rho });
    const st = ctx.material(inp.steel);
    need(st && st.fy > 0, `Steel '${inp.steel}' is not in the materials library.`);
    need(inp.colT / 1e3 < inp.colD / 4, 'Column wall thickness must be below a quarter of the diameter.');
    // drag and inertia are 90° out of phase: combine by square root of sum of squares
    const Mbase = Math.hypot(sup.baseMoment, sup.inertiaMoment), Vbase = Math.hypot(sup.baseShear, sup.inertiaShear);
    const sec = tubularSectionStress({ D: inp.colD, t: inp.colT / 1e3, M: Mbase * inp.gammaF, V: Vbase * inp.gammaF });
    const colMass = Math.PI * (inp.colD - inp.colT / 1e3) * inp.colT / 1e3 * inp.hubZ * 7850 / 1e3;   // t per column
    const checks = [
      check('Blade-tip cavitation at rated speed (top dead centre)', cav.sigmaRequired, cav.sigma, '-', 'Batten et al. (2008) Renewable Energy 33; DNV-ST-0164 §7 (cavitation)', 'Demand is the required cavitation number (−Cp,min × margin), capacity the available one at the shallowest blade position. If demand exceeds capacity, vapour bubbles form on the suction side: erosion, noise and loss of lift.'),
      check('Support-column bending stress at the seabed (extreme current and wave)', sec.vonMises / 1e6, st.fy / inp.gammaM / 1e6, 'MPa', 'IEC TS 62600-2:2019 §7 (ULS); DNV-ST-0164; elastic section check', 'Factored von Mises stress at the column base against design yield. Exceeding it means the column yields in the extreme event; local shell buckling and the foundation must be checked separately.'),
      check('Power coefficient against the blocked-channel limit', cpH * blk.gainCp, lim.CpMax, '-', 'Garrett & Cummins (2007) J. Fluid Mech. 588; Houlsby et al. (2008)', 'Hydrodynamic Cp of the rotor in the channel against the theoretical maximum for this blockage and Froude number. A value above the limit would indicate inputs outside the model.'),
    ];
    const gov = worst(checks);
    // plots
    const nShow = Math.round(15 * 86400 / 600);
    const [tt, uu] = thin([Array.from(tide.t.subarray(0, nShow), (v) => v / 86400), Array.from(tide.u.subarray(0, nShow))], 600);
    const pShow = uu.map((u) => { const us = Math.abs(u) * avg.powerVelocity; const i = Math.min(curve.u.length - 1, Math.max(0, Math.round((us - curve.u[0]) / 0.1))); return us < curve.u[0] ? 0 : curve.p[i] / 1e6; });
    const [ex, sp] = thin([Array.from(dist.exceedance, (v) => v * 100), Array.from(dist.speed)], 200);
    const a4 = linspace(0.2, 0.95, 40);
    const cpB = a4.map((x) => { try { return channelBlockage({ B, alpha4: x, Fr }).Cp; } catch { return null; } });
    const cp0 = a4.map((x) => channelBlockage({ B: 0, alpha4: x, Fr: 0 }).Cp);
    const zz = linspace(0.5, inp.depth, 40);
    const aep = en.aepWh * inp.nUnits, steelT = colMass * inp.nUnits;
    const dAvail = aep * 0.01 / (inp.avail / 100) / 1e6;
    return {
      kpis: [
        { key: 'cpRotor', label: 'Hydrodynamic power coefficient (open water)', value: round(cpH, 5), unit: '-', limit: round(BETZ_LIMIT, 4), sense: 'max', explain: 'Shaft power over ½ρAU³ at the design tip-speed ratio in unbounded flow. Well-designed tidal rotors reach 0.40–0.48; the Betz limit is 0.593.' },
        { key: 'ctRotor', label: 'Thrust coefficient (open water)', value: round(ct, 5), unit: '-', explain: 'Thrust over ½ρAU² at the design point: 0.8–0.9 for a rotor near optimum loading. Water is 800 times denser than air, so thrust per unit area is very large.' },
        { key: 'ratedSpeed', label: 'Rated current speed', value: round(ratedSpeed, 4), unit: 'm/s', explain: 'Rotor-averaged current speed at which rated power is first reached (0 if the rating is never reached). Matching it to roughly 70 % of the spring peak balances energy against structure.' },
        { key: 'thrustMax', label: 'Peak operating thrust', value: round(thrustPeak / 1e3, 5), unit: 'kN', explain: 'Largest rotor thrust over the operating range, at rated speed. Above rated the blades shed load, so this is the governing operating load on the support.' },
        { key: 'aep', label: 'Annual energy of the array', value: round(aep / 1e9, 5), unit: 'GWh/yr', explain: 'Energy of all turbines over a year of the synthesised tide, including availability. Tidal energy is predictable years ahead, unlike wind.' },
        { key: 'cf', label: 'Capacity factor', value: round(en.capacityFactor * 100, 4), unit: '%', explain: 'Mean power over rated power. 30–45 % is typical; it rises if the generator is down-rated relative to the spring peak.' },
        { key: 'powerDensity', label: 'Mean kinetic power density at hub height', value: round(kpd.mean / 1e3, 5), unit: 'kW/m²', explain: 'Time-mean of ½ρ|u|³ at hub height: the raw resource. Commercial sites are usually above 1.5–2 kW/m².' },
        { key: 'blockage', label: 'Blockage ratio', value: round(B, 5), unit: '-', explain: 'Rotor swept area over the flow section assigned to one device. Above a few percent the confined flow raises power and thrust beyond open-water values.' },
        { key: 'cpMaxBlocked', label: 'Maximum Cp in the blocked channel', value: round(lim.CpMax, 5), unit: '-', explain: 'Upper bound of the power coefficient for this blockage and Froude number from one-dimensional momentum theory; it tends to 16/27·(1 − B)⁻² for a rigid lid.' },
        { key: 'blockGain', label: 'Power gain from blockage', value: round((blk.gainCp - 1) * 100, 4), unit: '%', explain: 'Increase of power at the same disc operating point compared with open water. Thrust rises with it, so the structure pays for the gain.' },
        { key: 'sigma', label: 'Cavitation number at the blade tip', value: round(cav.sigma, 4), unit: '-', limit: round(cav.sigmaRequired, 4), sense: 'min', explain: 'Available pressure margin over the dynamic pressure of the tip at top dead centre. It must exceed the suction peak of the section (−Cp,min) times the required margin.' },
        { key: 'tipSpeed', label: 'Blade-tip speed at rated', value: round(tipSpeed, 4), unit: 'm/s', limit: round(cav.maxTipSpeed, 4), sense: 'max', explain: 'Tip speed at rated current against the largest tip speed that still avoids cavitation at this submergence.' },
        { key: 'thrustExt', label: 'Extreme rotor thrust (current, surge and wave)', value: round(extThrust / 1e3, 5), unit: 'kN', explain: 'Thrust at the extreme current with the orbital velocity of the extreme wave added at hub height, or the peak operating thrust if larger.' },
        { key: 'baseMoment', label: 'Extreme overturning moment at the seabed', value: round(Mbase / 1e6, 5), unit: 'MN·m', explain: 'Rotor thrust times hub height plus drag and inertia on the column (unfactored). It sizes the column, the foundation and the anchoring.' },
        { key: 'thrustRange', label: 'Wave-induced thrust range in the extreme sea', value: round(ext.thrustRange / 1e3, 5), unit: 'kN', explain: 'Thrust fluctuation per wave cycle from orbital velocity superposed on the current. It is the fatigue driver for blades and support in wave-exposed sites.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Power and thrust curves', x: { label: 'Rotor-averaged current speed', unit: 'm/s' }, y: { label: 'Power (kW) · thrust (kN)', unit: '-' },
          series: [{ name: 'Electrical power (kW)', x: curve.u.map((v) => round(v, 3)), y: curve.p.map((v) => round(v / 1e3, 4)) }, { name: 'Thrust (kN)', x: curve.u.map((v) => round(v, 3)), y: curve.thrust.map((v) => round(v / 1e3, 4)) }],
          explain: 'Electrical power and rotor thrust against current speed. Power follows the cube law to the rating; thrust peaks at rated speed and then falls as the blades pitch — the knee is the governing operating load.' },
        { type: 'line', title: 'Tide and power over a spring–neap cycle', x: { label: 'Time', unit: 'days' }, y: { label: 'Speed (m/s) · power (MW)', unit: '-' },
          series: [{ name: 'Surface current (m/s)', x: tt.map((v) => round(v, 4)), y: uu.map((v) => round(v, 4)) }, { name: 'Power per turbine (MW)', x: tt.map((v) => round(v, 4)), y: pShow.map((v) => round(v, 4)) }],
          explain: 'Fifteen days of the M2 + S2 tide and the power one turbine delivers. Output drops to zero four times a day at slack water and is capped at the rating around springs.' },
        { type: 'line', title: 'Current-speed exceedance', x: { label: 'Time exceeded', unit: '%' }, y: { label: 'Surface current speed', unit: 'm/s' },
          series: [{ name: 'Speed', x: ex.map((v) => round(v, 3)), y: sp.map((v) => round(v, 4)) }, { name: 'Rated speed (surface equivalent)', x: [0, 100], y: [ratedSpeed / avg.powerVelocity, ratedSpeed / avg.powerVelocity], dash: true }, { name: 'Cut-in', x: [0, 100], y: [inp.cutIn / avg.powerVelocity, inp.cutIn / avg.powerVelocity], dash: true }],
          explain: 'Fraction of the year each current speed is exceeded. Where the curve crosses the rated-speed line is the time spent at full power; below the cut-in line the turbine is idle.' },
        { type: 'line', title: 'Power coefficient of a disc in the channel', x: { label: 'Wake velocity ratio α₄', unit: '-' }, y: { label: 'Cp', unit: '-' },
          series: [{ name: `Blocked channel (B = ${fmt(B)})`, x: a4.filter((_, i) => cpB[i] != null).map((v) => round(v, 4)), y: cpB.filter((v) => v != null).map((v) => round(v, 5)) }, { name: 'Open water', x: a4.map((v) => round(v, 4)), y: cp0.map((v) => round(v, 5)) }, { name: 'Rotor (channel)', x: [round(lim.alpha4, 3)], y: [round(cpH * blk.gainCp, 4)] }],
          explain: 'Actuator-disc Cp against wake velocity ratio with and without blockage. The peak of the blocked curve is the Garrett–Cummins limit; the marker is this rotor at its channel operating point.' },
        { type: 'line', title: 'Current profile and rotor disc at spring peak', x: { label: 'Current speed', unit: 'm/s' }, y: { label: 'Height above seabed', unit: 'm' },
          series: [{ name: 'Spring peak profile', x: zz.map((z) => round(inp.springPeak * hubProfile(z), 4)), y: zz.map((z) => round(z, 3)) }, { name: 'Rotor disc', x: [inp.springPeak * hubProfile(inp.hubZ), inp.springPeak * hubProfile(inp.hubZ)], y: [inp.hubZ - R, inp.hubZ + R] }],
          explain: 'Power-law current profile over the depth with the vertical extent of the rotor. A hub high in the column sees more flow but less submergence for cavitation and more wave action.' },
        ...(sol ? [{ type: 'line', title: 'Spanwise blade loading at 2 m/s', x: { label: 'Radius', unit: 'm' }, y: { label: 'Force per unit length', unit: 'kN/m' },
          series: [{ name: 'Normal', x: sol.spanwise.map((q) => round(q.r, 4)), y: sol.spanwise.map((q) => round(q.Np / 1e3, 4)) }, { name: 'Tangential', x: sol.spanwise.map((q) => round(q.r, 4)), y: sol.spanwise.map((q) => round(q.Tp / 1e3, 4)) }],
          explain: 'Hydrodynamic load per metre of blade at a 2 m/s current (loads scale with speed squared). The steep outboard rise is why tidal blades are short, thick and heavily built.' }] : []),
      ],
      tables: [
        { title: 'Power curve', columns: [{ key: 'u', label: 'Current speed', unit: 'm/s' }, { key: 'p', label: 'Power', unit: 'kW' }, { key: 'cp', label: 'Cp (electrical)', unit: '-' }, { key: 'ct', label: 'Ct', unit: '-' }, { key: 'T', label: 'Thrust', unit: 'kN' }],
          rows: curve.u.filter((_, i) => i % 2 === 0).map((u, j) => { const i = 2 * j; return { u: round(u, 3), p: round(curve.p[i] / 1e3, 5), cp: round(curve.cp[i], 4), ct: round(curve.ct[i], 4), T: round(curve.thrust[i] / 1e3, 5) }; }),
          explain: 'Tabulated electrical power, coefficients and thrust per current-speed bin, as required for IEC TS 62600-200 reporting. Speeds are rotor-averaged.' },
        ...(rotor ? [{ title: 'Blade planform (optimum rotor)', columns: [{ key: 'r', label: 'Radius', unit: 'm' }, { key: 'c', label: 'Chord', unit: 'm' }, { key: 'tw', label: 'Twist', unit: 'deg' }, { key: 'a', label: 'a', unit: '-' }, { key: 'al', label: 'α', unit: 'deg' }, { key: 'cl', label: 'Cl', unit: '-' }],
          rows: thin([rotor.stations.map((_, i) => i)], 25)[0].map((i) => ({ r: round(rotor.stations[i].r, 4), c: round(rotor.stations[i].chord, 4), tw: round(rotor.stations[i].twistDeg, 3), a: round(sol.spanwise[i].a, 3), al: round(sol.spanwise[i].alpha / DEG, 3), cl: round(sol.spanwise[i].cl, 3) })),
          explain: 'Chord and twist of the Glauert–Schmitz optimum planform with the BEM induction and incidence at each station. A starting geometry for blade design, not a manufactured blade.' }] : []),
      ],
      notes: [
        `Rotor ${fmt(inp.D)} m, ${inp.rotorModel === 'bem' ? `optimum planform at λ = ${inp.tsr}, BEM Cp ${fmt(cpH)} / Ct ${fmt(ct)}` : `actuator disc a = ${fmt(inp.induction)}`}; blockage ${fmt(B * 100)} %, Froude number ${fmt(Fr)}, power gain ${fmt((blk.gainCp - 1) * 100)} %.`,
        `Tide: M2 ${fmt(sn.constituents[0].amplitude)} m/s + S2 ${fmt(sn.constituents[1]?.amplitude ?? 0)} m/s at the surface; rotor-averaged speed is ${fmt(avg.powerVelocity)} of the surface value; ${fmt(en.ratedFraction * 100)} % of the time at rated power.`,
        `Governing check: ${gov.name} (utilisation ${fmt(gov.uc)}). Cavitation: σ = ${fmt(cav.sigma)} against −Cp,min = ${fmt(-cpMin)} (screening estimate of the suction peak; vapour pressure ${fmt(pv, 4)} Pa at ${inp.Tsea} °C).`,
        `Extreme: design velocity ${fmt(ext.designVelocity)} m/s (wave orbital ${fmt(ext.waveVelocity)} m/s at hub), base shear ${fmt(Vbase / 1e3)} kN, base moment ${fmt(Mbase / 1e6)} MN·m unfactored.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `The turbine reaches ${fmt(inp.ratedMW)} MW at ${fmt(ratedSpeed)} m/s with a peak thrust of ${fmt(thrustPeak / 1e3)} kN; cavitation utilisation ${fmt(cav.utilisation)} and column stress utilisation ${fmt(checks[1].uc)}.`,
          gov.uc > 1 ? (gov === checks[0] ? `Cavitation margin is not met: reduce the tip-speed ratio to ${fmt(Math.max(2, cav.maxTipSpeed / Math.max(Ucav, 0.1)), 2)} or lower the hub by ${fmt(Math.max(1, (cav.sigmaRequired - cav.sigma) * 0.5 * rho * cav.relativeSpeed ** 2 / (rho * G)), 2)} m.` : `Increase the column wall to about ${fmt(inp.colT * gov.uc * 1.05, 3)} mm or the diameter by ${fmt((Math.sqrt(gov.uc) - 1) * 100, 2)} %.`) : 'Margins are met: confirm the suction peak with section data (the estimate ignores leading-edge peaks), and the thrust with basin or supplier measurements before sizing the foundation.',
          'IEC TS 62600-2:2019 (design requirements for marine energy systems) and DNV-ST-0164 (2015) for tidal turbines; power performance per IEC TS 62600-200:2013.',
          `Support columns: ${fmt(steelT)} t of steel for ${inp.nUnits} turbine(s). ${steelCo2Text(steelT)} The seabed footprint is limited to the foundations, and energy is predictable — valuable for system integration.`),
        adv('economic', 'info', `The array yields ${fmt(aep / 1e9, 4)} GWh/yr at a capacity factor of ${fmt(en.capacityFactor * 100)} %; each point of availability is worth ${fmt(dAvail)} MWh/yr.`,
          `Test the rating against the resource: the rated speed is ${fmt(ratedSpeed / (inp.springPeak * avg.powerVelocity) * 100)} % of the rotor-averaged spring peak — a lower rating raises the capacity factor and lowers thrust (${fmt(thrustPeak / 1e3)} kN now) for a small energy loss; retrieval for maintenance needs a vessel in slack-water windows, so availability dominates the yield.`,
          'LCOE conventions (IEA-OES / IRENA) with whole-life costing per ISO 15686-5: rating, availability and retrieval strategy move the cost of tidal energy more than rotor efficiency.',
          energyCo2Text(aep / 1e6)),
      ],
      publish: { 'tidal.ratedPower': P * inp.nUnits, 'tidal.aep': aep, 'tidal.thrustMax': thrustPeak, 'tidal.baseMoment': Mbase, 'tidal.capacityFactor': en.capacityFactor },
      quantities: { steel_t: round(steelT, 4) },
    };
  },
  verify: [
    { name: 'Ideal actuator disc at a = 1/3: Cp = 16/27, Ct = 8/9', source: 'Lanchester–Betz limit; Burton et al., Wind Energy Handbook §3.2',
      inputs: { rotorModel: 'disc', induction: 1 / 3 }, expect: [{ kpi: 'cpRotor', value: 0.592593, tol: 1e-4 }, { kpi: 'ctRotor', value: 0.888889, tol: 1e-4 }] },
    { name: 'Rigid-lid blocked channel: Cp,max = (16/27)(1 − B)⁻² at B = 0.1 (low Froude number)', source: 'Garrett & Cummins (2007) J. Fluid Mech. 588:243–251, Eq. 3.5',
      // B = π·10² / (78.54 × 40) = 0.1000; 0.592593 / 0.81 = 0.731596. Spring peak 0.5 m/s gives Fr ≈ 0.02, i.e. the rigid-lid limit within 0.3 %.
      inputs: { rotorModel: 'disc', D: 20, width: 78.54, depth: 40, springPeak: 0.5, surge: 0.6, cutIn: 0.2, ratedMW: 0.05 }, expect: [{ kpi: 'cpMaxBlocked', value: 0.731596, tol: 0.004 }, { kpi: 'blockage', value: 0.1, tol: 0.001 }] },
    { name: 'M2-only tide: mean kinetic power density = ½ρ·a³·4/(3π) at hub height', source: 'Hand calculation: time-mean of |cos|³ = 4/(3π)',
      // hub speed amplitude = 3 × (18/40)^(1/7) = 3 × 0.892190 = 2.676570 m/s; ½ × 1025 × 2.676570³ × 0.424413 = 4170.9 W/m² = 4.171 kW/m²
      inputs: { springPeak: 3, neapRatio: 1, depth: 40, hubZ: 18, profileN: 7, surge: 3.5 }, expect: [{ kpi: 'powerDensity', value: 4.171, tol: 0.005 }] },
  ],
};
