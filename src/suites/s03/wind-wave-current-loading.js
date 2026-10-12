// Coupled wind, wave and current loading on a monopile-supported turbine: steady rotor thrust from blade element
// momentum theory, tower drag in a sheared wind, Morison wave loading with diffraction correction and current,
// combined into shear and overturning moment at the mudline.
import { bemSolve, inductionFromCt, actuatorDisc } from '../../engine/aero/bem.js';
import { pileForceClosedForm } from '../../engine/hydro/morison.js';
import { DEG, need, run, fmt, round, adv, check, steelCo2Text, tube, turbineInputs, supportInputs, turbine, operatingPoint, LIMIT_AERO, LIMIT_HIFI, HELP_SUPPLIER } from './_common.js';
import { waveCurrentLoad } from './_hydro.js';

export default {
  id: 's03-wind-wave-current-loading', suite: 3, title: 'Coupled wind, wave and current loading',
  summary: 'Rotor thrust, tower wind drag, wave and current forces on the monopile and their combination into shear force and overturning moment at the mudline for a chosen operating condition.',
  covers: [3], fidelity: 'analytical',
  standards: ['IEC 61400-3-1:2019 §6.4 (marine conditions), §7.3–7.4', 'IEC 61400-1:2019 §6.3.2.2 (normal wind profile)', 'DNV-RP-C205 (2019) §3 (waves), §4 (current), §5 (wind), §6 (Morison, MacCamy–Fuchs)', 'DNV-ST-0437 (2016) §4', 'API RP 2A-WSD (21st ed.) §2.3.1 (wave–current kinematics)'],
  models: ['Blade element momentum theory with Prandtl tip-loss and hub-loss and Glauert/Buhl high-thrust correction', 'Actuator-disc momentum theory (induction from thrust coefficient)', 'Aerodynamic thrust and torque from the steady operating schedule (variable speed, pitch to feather)', 'Power-law wind shear and drag on the tower', 'Airy wave theory and Morison equation integrated over depth', 'MacCamy–Fuchs diffraction for large diameter', 'Power-law current profile with nonlinear wave–current drag interaction', 'Force and moment equilibrium about the mudline'],
  limits: `Steady rotor loads (no turbulence, gusts or dynamics — see the coupled-simulation and extreme-loading modules); linear waves integrated to the still-water level; wind, wave and current forces are combined as collinear vectors rotated by the stated misalignment. ${LIMIT_AERO} ${LIMIT_HIFI}`,
  publishes: ['owt.thrustRated', 'owt.rna'],
  searchTerms: 'offshore wind turbine thrust wave current load monopile Morison MacCamy-Fuchs',
  inputs: [
    ...turbineInputs(),
    ...supportInputs({ embed: false }).filter((i) => i.key !== 'towerMassFactor' && i.key !== 't' && i.key !== 'towerTb' && i.key !== 'towerTt'),
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'level', label: 'Water level above mean (tide + surge)', unit: 'm', type: 'number', default: 1.5, min: -5, max: 10, group: 'Site', role: 'bc', help: 'Raises the wave loading zone and the lever arm of the wave force.' },
    { key: 'Vhub', label: 'Wind speed at hub height', unit: 'm/s', type: 'number', default: 11.5, min: 2, max: 70, group: 'Wind', role: 'bc', uncertain: { dist: 'normal', cov: 0.05 }, help: 'Thrust peaks near rated wind speed. Above cut-out the rotor is parked and feathered.' },
    { key: 'shear', label: 'Wind shear exponent', unit: '-', type: 'number', default: 0.14, min: 0, max: 0.4, group: 'Wind', role: 'bc', help: 'IEC 61400-3-1 normal wind profile offshore: 0.14.' },
    { key: 'yaw', label: 'Yaw misalignment', unit: 'deg', type: 'number', default: 0, min: -45, max: 45, group: 'Wind', role: 'bc' },
    { key: 'rhoAir', label: 'Air density', unit: 'kg/m³', type: 'number', default: 1.225, min: 1.0, max: 1.4, group: 'Wind', role: 'environment' },
    { key: 'cdTower', label: 'Tower drag coefficient', unit: '-', type: 'number', default: 0.6, min: 0.3, max: 1.2, group: 'Wind', role: 'load', help: 'Smooth cylinder at high Reynolds number 0.5–0.7 (DNV-RP-C205 §5).' },
    { key: 'Hs', label: 'Significant wave height', unit: 'm', type: 'number', default: 4, min: 0, max: 16, group: 'Waves and current', role: 'bc', bind: 'site.Hs1', uncertain: { dist: 'lognormal', cov: 0.1 }, help: 'Sea state acting with the wind condition; the design wave is taken as H = 1.86 Hs.' },
    { key: 'Tp', label: 'Wave period', unit: 's', type: 'number', default: 8.5, min: 3, max: 22, group: 'Waves and current', role: 'bc', bind: 'site.Tp1' },
    { key: 'Uc', label: 'Surface current speed', unit: 'm/s', type: 'number', default: 0.8, min: 0, max: 3, group: 'Waves and current', role: 'bc', bind: 'site.Uc1', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'alphaC', label: 'Current profile exponent', unit: '-', type: 'number', default: 0.1429, min: 0, max: 0.5, group: 'Waves and current', role: 'bc', help: '1/7 for tidal current (DNV-RP-C205 §4.1.4); 0 gives a uniform profile.' },
    { key: 'misDeg', label: 'Wind–wave misalignment', unit: 'deg', type: 'number', default: 0, min: 0, max: 90, group: 'Waves and current', role: 'bc' },
    { key: 'mg', label: 'Marine growth thickness', unit: 'mm', type: 'number', default: 50, min: 0, max: 200, group: 'Hydrodynamic coefficients', role: 'environment', help: 'Added to the radius below the splash zone; raises drag diameter and inertia volume.' },
    { key: 'Cd', label: 'Drag coefficient', unit: '-', type: 'number', default: 0.9, min: 0.4, max: 1.6, group: 'Hydrodynamic coefficients', role: 'environment', uncertain: { dist: 'normal', cov: 0.15 } },
    { key: 'Cm', label: 'Inertia coefficient', unit: '-', type: 'number', default: 2, min: 1, max: 2.5, group: 'Hydrodynamic coefficients', role: 'environment' },
    { key: 'mcf', label: 'MacCamy–Fuchs diffraction correction', type: 'bool', default: true, group: 'Hydrodynamic coefficients', role: 'environment', help: 'Replaces Cm by the diffraction value when the pile is large against the wavelength.' },
    { key: 'nZ', label: 'Integration stations over depth', type: 'int', default: 40, min: 8, max: 400, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nZ', label: 'Depth stations', levels: [10, 20, 40, 80], qoi: 'waveMoment' },
  compute(inp) {
    const tb = turbine(inp, { rhoAir: inp.rhoAir });
    const d = inp.depth + inp.level;
    need(d > 3, 'The water level correction leaves less than 3 m of water: check the level and the depth.');
    const op = operatingPoint(tb, inp.Vhub);
    const hubAbove = inp.hubHeight + inp.depth; // lever arms are taken about the mudline
    // Rotor: BEM at the steady operating point (or parked and feathered).
    const bem = run(() => bemSolve({ rotor: tb.rotor, U: inp.Vhub, omega: op.parked ? 0 : op.omega, pitchDeg: op.parked ? 90 : op.pitchDeg, yawDeg: inp.yaw, rho: inp.rhoAir, shearExp: inp.shear, hubHeight: inp.hubHeight }), 'Rotor aerodynamics');
    const thrust = Math.abs(bem.thrust);
    const ct = thrust / (0.5 * inp.rhoAir * tb.area * inp.Vhub ** 2);
    const aInd = op.parked ? 0 : inductionFromCt(Math.min(ct, 1.9));
    const ideal = actuatorDisc(Math.min(aInd, 0.5)); // ideal-disc Cp for the same induction (upper bound of the BEM value)
    // Tower and dry monopile: drag in the sheared wind, integrated in slices.
    const zTop = tb.zTowerTop, zLow = Math.max(inp.level, 0.1);
    need(zTop > inp.zInt + 5, 'The tower top must be above the interface level.');
    let Ft = 0, Mt = 0;
    const nS = 60;
    for (let i = 0; i < nS; i++) {
      const z = zLow + ((i + 0.5) * (zTop - zLow)) / nS, dz = (zTop - zLow) / nS;
      const D = z < inp.zInt ? inp.D : inp.towerDb + ((inp.towerDt - inp.towerDb) * (z - inp.zInt)) / (zTop - inp.zInt);
      const U = inp.Vhub * (z / inp.hubHeight) ** inp.shear;
      const f = 0.5 * inp.rhoAir * inp.cdTower * D * U * U * dz;
      Ft += f; Mt += f * (z + inp.depth);
    }
    // Waves and current on the wet pile.
    const Dh = inp.D + (2 * inp.mg) / 1e3;
    const H = 1.86 * inp.Hs;
    const hyd = run(() => waveCurrentLoad({ H, T: inp.Tp, d, D: Dh, Cd: inp.Cd, Cm: inp.Cm, mcf: inp.mcf, rho: inp.rhoW, uc0: inp.Uc, alpha: inp.alphaC, nZ: inp.nZ }), 'Wave and current loading');
    const waveOnly = H > 0 ? waveCurrentLoad({ H, T: inp.Tp, d, D: Dh, Cd: inp.Cd, Cm: inp.Cm, mcf: inp.mcf, rho: inp.rhoW, uc0: 0, nZ: inp.nZ }) : { shearMax: 0, momentMax: 0 };
    const curOnly = inp.Uc > 0 ? waveCurrentLoad({ H: 0, T: inp.Tp, d, D: Dh, Cd: inp.Cd, Cm: inp.Cm, mcf: false, rho: inp.rhoW, uc0: inp.Uc, alpha: inp.alphaC, nZ: inp.nZ }) : { shearMax: 0, momentMax: 0 };
    const closed = H > 0 ? pileForceClosedForm({ H, T: inp.Tp, depth: d, D: Dh, Cd: inp.Cd, Cm: hyd.cmEff, rho: inp.rhoW }) : { Fmax: 0, Mmax: 0, FD: 0, FI: 0 };
    // Combination: aerodynamic loads along the wind, hydrodynamic loads rotated by the misalignment. The hydro lever is
    // referred to the mean mudline (water level only shifts the loaded zone).
    const c = Math.cos(inp.misDeg * DEG), s = Math.sin(inp.misDeg * DEG);
    const Fa = thrust + Ft, Ma = thrust * hubAbove + Mt;
    const Fh = hyd.shearMax, Mh = hyd.momentMax;
    const shear = Math.hypot(Fa + Fh * c, Fh * s), moment = Math.hypot(Ma + Mh * c, Mh * s);
    const interaction = hyd.momentMax / Math.max(waveOnly.momentMax + curOnly.momentMax, 1e-9);
    const vs = tb.schedule.U, sched = tb.schedule;
    const pile = { D: inp.D };
    const regime = hyd.ka > 0.2 * Math.PI ? 'diffraction (D/L > 0.2)' : 'slender-body Morison (D/L < 0.2)';
    const KC = H > 0 ? (Math.PI * H / inp.Tp / Math.tanh(hyd.k * d)) * inp.Tp / Dh : 0;
    const aeroShare = Ma / Math.max(Ma + Mh, 1e-9);
    const checks = [
      check('Wave height against the breaking limit', H, hyd.Hb, 'm', 'DNV-RP-C205 (2019) §3.4.6 (Miche / depth-limited breaking)', 'Linear-wave Morison loading is only valid for non-breaking waves. Above the limit the wave breaks at the site and slamming loads must be added (see the extreme-loading module).'),
      check('Rotor thrust coefficient against the momentum-theory limit', ct, 1, '-', 'Glauert (1935); Buhl (2005) NREL/TP-500-36834', 'Above a thrust coefficient of about 0.96 the rotor is in the turbulent-wake state where simple momentum theory fails and the empirical correction governs; thrust predictions carry larger uncertainty.'),
    ];
    const span = bem.spanwise || [];
    const dThrustPerMs = (operatingPoint(tb, Math.min(inp.Vhub + 0.5, 60)).thrust - operatingPoint(tb, Math.max(inp.Vhub - 0.5, 0.5)).thrust);
    const Wneed = moment * 1.35 / (355e6 / 1.1), tNeed = Wneed / (Math.PI * inp.D * inp.D / 4) * 1e3;
    return {
      kpis: [
        { key: 'thrust', label: 'Rotor thrust at this wind speed', unit: 'kN', value: thrust / 1e3, explain: `Steady aerodynamic thrust on the rotor from blade element momentum theory (${op.parked ? 'parked, blades feathered' : `rotor at ${fmt(op.omega * 30 / Math.PI)} rpm, pitch ${fmt(op.pitchDeg)}°`}). ${HELP_SUPPLIER}` },
        { key: 'thrustRated', label: 'Maximum steady thrust (near rated wind)', unit: 'kN', value: tb.thrustRated / 1e3, explain: `Peak of the steady thrust curve, at ${fmt(tb.windAtMaxThrust)} m/s. It is the reference rotor load for foundation sizing; gusts add 20–40 % dynamically.` },
        { key: 'ct', label: 'Thrust coefficient', unit: '-', value: ct, explain: 'Thrust divided by ½ρAV². About 0.75–0.9 below rated wind, falling quickly above rated as the blades pitch.' },
        { key: 'cp', label: 'Power coefficient', unit: '-', value: op.parked ? 0 : bem.Cp, explain: `Aerodynamic power over the power in the wind through the rotor disc. The Betz limit is 0.593 and an ideal disc at this induction would give ${fmt(ideal.cp ?? ideal.Cp ?? 0)}; modern rotors reach 0.45–0.50.` },
        { key: 'aInd', label: 'Axial induction factor (momentum theory)', unit: '-', value: aInd, explain: 'Fractional slow-down of the wind at the rotor plane implied by the thrust coefficient. One third is the Betz optimum; above 0.4 the turbulent-wake correction applies.' },
        { key: 'towerDrag', label: 'Wind drag on tower and dry pile', unit: 'kN', value: Ft / 1e3, explain: 'Drag of the sheared wind on the tower and the monopile above the water. Small in production, dominant for the parked turbine in a storm.' },
        { key: 'waveShear', label: 'Maximum wave force (waves only)', unit: 'kN', value: waveOnly.shearMax / 1e3, explain: 'Largest horizontal Morison force over the wave cycle on the wet pile without current. On a large monopile it is mainly inertia force, a quarter period before the crest.' },
        { key: 'waveMoment', label: 'Wave + current moment at the mudline', unit: 'MN·m', value: hyd.momentMax / 1e6, explain: 'Largest overturning moment from waves and current together over the wave cycle, about the mudline.' },
        { key: 'currentForce', label: 'Current drag force (current only)', unit: 'kN', value: curOnly.shearMax / 1e3, explain: 'Steady drag of the current alone on the wet pile. Small by itself but it amplifies the wave drag because velocities add before squaring.' },
        { key: 'cmEff', label: 'Effective inertia coefficient', unit: '-', value: hyd.cmEff, explain: 'Inertia coefficient used: 2.0 for a slender pile, reduced by diffraction (MacCamy–Fuchs) when the diameter is large against the wavelength.' },
        { key: 'dOverL', label: 'Diameter / wavelength', unit: '-', value: Dh / hyd.L, explain: 'Above about 0.2 the pile disturbs the wave field and diffraction theory replaces the Morison inertia term.' },
        { key: 'shear', label: 'Combined shear at the mudline', unit: 'MN', value: shear / 1e6, explain: 'Rotor thrust, tower drag and the maximum wave and current force added as vectors with the stated misalignment (characteristic, unfactored, steady wind).' },
        { key: 'moment', label: 'Combined overturning moment at the mudline', unit: 'MN·m', value: moment / 1e6, explain: 'Sum of the aerodynamic and hydrodynamic overturning moments about the mudline. The foundation is sized on this quantity and its fatigue-driving fluctuation.' },
        { key: 'aeroShare', label: 'Aerodynamic share of the moment', unit: '%', value: aeroShare * 100, explain: 'Share of the mudline moment caused by wind loads. High in production; in a storm with the rotor parked the wave share dominates.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Steady thrust and power curves', x: { label: 'Hub wind speed', unit: 'm/s' }, y: { label: 'Thrust (kN) / power (10 kW)', unit: '' },
          series: [{ name: 'Thrust, kN', x: vs, y: sched.thrust.map((v) => v / 1e3) }, { name: 'Electrical power, ×10 kW', x: vs, y: sched.power.map((v) => v / 1e4) }, { name: 'This condition', x: [inp.Vhub], y: [thrust / 1e3] }],
          explain: 'Thrust rises with the square of wind speed up to rated, then falls as the blades pitch to hold rated power. The peak near rated wind sets the operating design load.' },
        { type: 'line', title: 'Blade loading along the span', x: { label: 'Radius', unit: 'm' }, y: { label: 'Force per length', unit: 'kN/m' },
          series: [{ name: 'Normal (thrust direction)', x: span.map((q) => q.r), y: span.map((q) => q.Np / 1e3) }, { name: 'Tangential (torque direction)', x: span.map((q) => q.r), y: span.map((q) => q.Tp / 1e3) }, { name: 'Tip/hub loss factor ×1', x: span.map((q) => q.r), y: span.map((q) => q.F), dash: true }],
          explain: 'Aerodynamic force per metre of one blade. Loading grows towards the tip and drops to zero at the tip through the Prandtl loss factor (dashed).' },
        { type: 'line', title: 'Wave and current force over depth', x: { label: 'Force per length', unit: 'kN/m' }, y: { label: 'Elevation', unit: 'm' },
          series: [{ name: `Phase of maximum moment (${fmt(hyd.phaseDeg)}°)`, x: hyd.f.map((v) => v / 1e3), y: hyd.z }],
          explain: 'Distribution of the hydrodynamic force on the pile at the wave phase giving the largest mudline moment. Inertia loading decays with depth; the concentration near the surface raises the lever arm.' },
        { type: 'bar', title: 'Contributions to the mudline moment', categories: ['Rotor thrust', 'Tower drag', 'Waves', 'Current', 'Wave–current interaction'], series: [{ name: 'Moment, MN·m', y: [thrust * hubAbove / 1e6, Mt / 1e6, waveOnly.momentMax / 1e6, curOnly.momentMax / 1e6, Math.max(hyd.momentMax - waveOnly.momentMax - curOnly.momentMax, 0) / 1e6] }], x: { label: 'Source' }, y: { label: 'Moment', unit: 'MN·m' },
          explain: 'Breakdown of the overturning moment by source. The interaction bar is the extra drag caused by adding wave and current velocities before squaring.' },
      ],
      tables: [
        { title: 'Load summary at the mudline', columns: [{ key: 'k', label: 'Source' }, { key: 'F', label: 'Force', unit: 'kN' }, { key: 'M', label: 'Moment', unit: 'MN·m' }, { key: 'h', label: 'Lever arm', unit: 'm' }],
          rows: [['Rotor thrust', thrust, thrust * hubAbove], ['Tower and dry-pile drag', Ft, Mt], ['Waves and current (maximum over the cycle)', Fh, Mh], ['Waves only', waveOnly.shearMax, waveOnly.momentMax], ['Current only', curOnly.shearMax, curOnly.momentMax], ['Morison closed form, waves only (Airy)', closed.Fmax, closed.Mmax], ['Combined', shear, moment]]
            .map(([k, F, M]) => ({ k, F: round(F / 1e3), M: round(M / 1e6), h: F > 0 ? round(M / F) : 0 })),
          explain: 'Force, moment and lever arm of each load source about the mudline. The closed-form row is an independent check of the numerical depth integration.' },
        { title: 'Hydrodynamic regime', columns: [{ key: 'k', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['Design wave height H = 1.86 Hs', H, 'm'], ['Wavelength', hyd.L, 'm'], ['ka = πD/L', hyd.ka, '-'], ['Keulegan–Carpenter number at the surface', KC, '-'], ['Breaking wave height at the site', hyd.Hb, 'm'], ['Recommended wave theory', hyd.recommended, ''], ['Loading regime', regime, ''], ['Drag part of the wave force', closed.FD / 1e3, 'kN'], ['Inertia part of the wave force', closed.FI / 1e3, 'kN'], ['Wave–current interaction factor on moment', interaction, '-']].map(([k, v, u]) => ({ k, v: typeof v === 'number' ? round(v) : v, u })),
          explain: 'Parameters that decide which wave-force model applies: D/L for diffraction, KC for drag against inertia, and the wave theory recommended for this height, period and depth.' },
      ],
      notes: [`${tb.name}; ${op.parked ? 'wind outside the operating range — rotor parked and feathered' : `steady operating point at ${fmt(inp.Vhub)} m/s`}.`, `Loading regime: ${regime}; the recommended wave theory is ${hyd.recommended}${hyd.recommended !== 'airy' ? ' — use the extreme-loading module for nonlinear kinematics' : ''}.`,
        `Aerodynamic loads give ${fmt(aeroShare * 100)} % of the mudline moment.`, 'Rotor loads come from generic polars and a geometrically scaled planform: replace them with the supplier\'s thrust curve for design.'],
      advice: [
        adv('engineering', checks.some((q) => q.uc > 1) ? 'warn' : 'ok', `Combined characteristic loads at the mudline: shear ${fmt(shear / 1e6)} MN, moment ${fmt(moment / 1e6)} MN·m (${fmt(aeroShare * 100)} % aerodynamic). Maximum steady thrust is ${fmt(tb.thrustRated / 1e3)} kN at ${fmt(tb.windAtMaxThrust)} m/s.`,
          `Treat the rotor loads as preliminary (generic polars, scaled planform) and replace them with the turbine supplier's load document; apply dynamic amplification and gust effects (coupled simulation) and the partial load factors before sizing; near this wind speed the thrust changes by ${fmt(dThrustPerMs / 1e3)} kN per m/s. ${checks[0].uc > 1 ? 'The design wave exceeds the breaking limit: add slamming loads.' : ''}`,
          'IEC 61400-3-1:2019 §7.4 and §7.6; DNV-RP-C205 (2019) §6',
          `With a load factor of 1.35 this moment needs a monopile section modulus of ${fmt(Wneed)} m³, about ${fmt(tNeed)} mm of wall at D = ${fmt(pile.D)} m for yield alone (${fmt(tube(inp.D, Math.max(tNeed, 10) / 1e3).A * 7.85)} t/m). Per metre of pile: ${steelCo2Text(tube(inp.D, Math.max(tNeed, 10) / 1e3).A * 7.85)}`),
        adv('economic', 'info', `The rotor contributes ${fmt(thrust * hubAbove / 1e6)} MN·m and waves with current ${fmt(Mh / 1e6)} MN·m; every metre of extra hub height adds ${fmt(thrust / 1e6)} MN·m to the foundation moment.`,
          'Trade hub height and rotor size against foundation steel in the layout optimisation: yield rises with hub height through wind shear, but foundation moment and steel rise almost linearly with it.',
          'LCOE-based design optimisation (energy yield against CAPEX); quantity-based estimating per AACE RP 18R-97 Class 4'),
      ],
      publish: { 'owt.thrustRated': tb.thrustRated, 'owt.rna': tb.rnaMass },
    };
  },
  verify: [
    { name: 'Current drag on the pile with a 1/7 power-law profile: F = ½ρ·Cd·D·U0²·d·7/9', source: 'Morison drag integrated over depth, ∫(z′/d)^(2/7) dz′ = 7d/9 (hand calculation): 0.5 × 1025 × 0.9 × 6 × 1² × 25 × 7/9 = 53.81 kN',
      inputs: { Hs: 0, Uc: 1, alphaC: 1 / 7, level: 0, mg: 0, Cd: 0.9, D: 6, depth: 25, nZ: 80 }, expect: [{ kpi: 'currentForce', value: 53.81, tol: 0.002 }] },
    { name: 'Inertia-dominated Airy wave force: F = ρ·Cm·(πD²/4)·g·(H/2)·tanh(kd)', source: 'Morison et al. (1950) closed form (Dean & Dalrymple 1991, §8): H = 2 m, T = 10 s, d = 25 m → k = 0.04820 1/m, tanh(kd) = 0.8352; F = 1025 × 2 × 28.274 × 9.80665 × 1 × 0.8352 = 474.7 kN',
      inputs: { Hs: 2 / 1.86, Tp: 10, Uc: 0, level: 0, mg: 0, Cd: 0, Cm: 2, mcf: false, D: 6, depth: 25, nZ: 80 }, expect: [{ kpi: 'waveShear', value: 474.7, tol: 0.002 }] },
  ],
};
