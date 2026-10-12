// Dynamic response and resonance of the turbine–tower–foundation system: frequency-domain response of the first
// fore-aft mode to turbulent thrust and wave loading, damping budget, 1P / 3P harmonic amplification and nacelle
// acceleration.
import { ntmSigma, peakFactor } from '../../engine/aero/turbulence.js';
import { unsteadyBem, harmonicAmplitudes } from '../../engine/aero/unsteady.js';
import { dynamicAmplification } from '../../engine/structural/framedynamics.js';
import { rayleighDamping } from '../../engine/structural/framemodal.js';
import { run, fmt, round, adv, check, worst, sev, steelCo2Text, thin, turbineInputs, supportInputs, foundationInputs, steelInputs, steelOf, turbine, steelMass, LIMIT_AERO, LIMIT_HIFI, HELP_SUPPLIER } from './_common.js';
import { dynamicsInputs, supportStructure, modalSpectra } from './_dynamics.js';

const TWO_PI = 2 * Math.PI;
const daf = (r, z) => dynamicAmplification({ ratio: r, zeta: z }).daf;

export default {
  id: 's03-dynamic-response', suite: 3, title: 'Dynamic response and resonance',
  summary: 'Displacement, acceleration and mudline-moment response of the first support-structure mode to turbulence and waves, the damping budget including aerodynamic damping, and the amplification at rotor (1P) and blade-passing (3P) frequencies.',
  covers: [4], fidelity: 'reduced-order',
  standards: ['DNV-ST-0126 (2021) §4.5 (dynamic analysis, resonance)', 'IEC 61400-3-1:2019 §7.5', 'IEC 61400-1:2019 §6.3.2.3 (normal turbulence model), Annex C (Kaimal spectrum)', 'DNV-RP-C205 (2019) §3.5 (JONSWAP), §6 (Morison, MacCamy–Fuchs)', 'Salzmann & van der Tempel (2005) aerodynamic damping'],
  models: ['Modal analysis of the stick model (first fore-aft mode)', 'Modal damping model: structural + soil + hydrodynamic + aerodynamic', 'Rayleigh damping coefficients for time-domain use', 'Kaimal turbulence spectrum with rotor-averaging filter', 'JONSWAP wave spectrum, Airy wave theory, Morison inertia with MacCamy–Fuchs diffraction', 'Aerodynamic thrust slope from the BEM coefficient surface', 'Harmonic resonance screening (dynamic amplification factor)', 'Davenport peak factor', 'Unsteady blade element momentum with Øye dynamic inflow, Beddoes–Leishman dynamic stall and tower shadow (optional 3P load)'],
  limits: `Single-mode, linear, frequency-domain model: wave drag is not linearised, the controller is represented by a first-order blend between the steady-schedule and frozen-control thrust slopes at the pitch-loop frequency, and side-side response is not computed. ${LIMIT_AERO} ${LIMIT_HIFI}`,
  searchTerms: 'offshore wind turbine monopile dynamic response aerodynamic damping resonance 1P 3P',
  inputs: [
    ...turbineInputs(),
    ...supportInputs({ embed: false }),
    ...steelInputs(),
    ...foundationInputs(),
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'Vhub', label: 'Mean wind speed at hub height', unit: 'm/s', type: 'number', default: 12, min: 2, max: 60, group: 'Environment', role: 'bc', uncertain: { dist: 'normal', cov: 0.05 }, help: 'Ten-minute mean. Above cut-out the rotor is parked and aerodynamic damping vanishes.' },
    { key: 'Iref', label: 'Reference turbulence intensity', unit: '-', type: 'number', default: 0.14, min: 0.05, max: 0.25, group: 'Environment', role: 'bc', help: 'IEC class A 0.16, B 0.14, C 0.12; the normal turbulence model gives σ1 = Iref (0.75 V + 5.6).' },
    { key: 'Hs', label: 'Significant wave height', unit: 'm', type: 'number', default: 2.5, min: 0, max: 15, group: 'Environment', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.1 } },
    { key: 'Tp', label: 'Peak wave period', unit: 's', type: 'number', default: 7, min: 2, max: 22, group: 'Environment', role: 'bc' },
    { key: 'gamma', label: 'JONSWAP peak enhancement', unit: '-', type: 'number', default: 3.3, min: 1, max: 7, group: 'Environment', role: 'bc' },
    ...dynamicsInputs(),
    { key: 'imb', label: 'Rotor mass imbalance moment', unit: 'kg·m', type: 'number', default: 1800, min: 0, max: 1e5, group: 'Rotor excitation', role: 'load', help: 'Residual imbalance (mass × radius) after balancing; about 0.5 % of blade mass at its centre of gravity. ' + HELP_SUPPLIER },
    { key: 'p3Mode', label: 'Blade-passing (3P) thrust fluctuation', type: 'select', default: 'input', group: 'Rotor excitation', role: 'load', options: [{ value: 'input', label: 'Entered as a share of mean thrust' }, { value: 'bem', label: 'Computed: unsteady BEM with tower shadow, shear and dynamic stall (slow)' }],
      help: 'The computed option runs a blade-resolved unsteady blade-element-momentum simulation over three rotor revolutions at rated wind (Øye dynamic inflow, Beddoes–Leishman dynamic stall, potential-flow tower shadow).' },
    { key: 'shear', label: 'Wind shear exponent', unit: '-', type: 'number', default: 0.14, min: 0, max: 0.4, group: 'Rotor excitation', role: 'bc', showIf: { key: 'p3Mode', in: ['bem'] } },
    { key: 'p3Frac', label: '3P thrust fluctuation / mean thrust', unit: '%', type: 'number', default: 3, min: 0, max: 20, group: 'Rotor excitation', role: 'load', showIf: { key: 'p3Mode', in: ['input'] }, help: 'Amplitude of the blade-passing thrust variation from tower shadow and wind shear as a share of the mean thrust.' },
    { key: 'accLim', label: 'Allowable peak nacelle acceleration', unit: 'm/s²', type: 'number', default: 2.0, min: 0.2, max: 10, group: 'Acceptance criteria', role: 'criteria', help: 'Operational limit of the turbine supplier (typically 0.15–0.3 g fore-aft).' },
    { key: 'dafLim', label: 'Allowable amplification at 1P / 3P', unit: '-', type: 'number', default: 2.5, min: 1.1, max: 20, group: 'Acceptance criteria', role: 'criteria', help: 'Screening limit on the dynamic amplification of rotor harmonics at any speed in the operating range.' },
    { key: 'nFreq', label: 'Frequency points', type: 'int', default: 400, min: 60, max: 3200, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nFreq', label: 'Frequency points', levels: [200, 400, 800, 1600], qoi: 'sigmaM' },
  compute(inp, ctx) {
    const steel = steelOf(inp, ctx);
    const tb = turbine(inp);
    const str = supportStructure(inp, tb, steel, { nElem: 30, nModes: 2 });
    const md = str.modes[0], f1 = md.frequency;
    const sigma1 = ntmSigma({ Vhub: inp.Vhub, Iref: inp.Iref });
    const sp = run(() => modalSpectra({ str, tb, inp, V: inp.Vhub, sigma1, Hs: inp.Hs, Tp: inp.Tp, gamma: inp.gamma, nFreq: inp.nFreq }), 'Spectral response');
    const zeta = sp.zetaTotal;
    const T10 = 600;
    const nuM = Math.max(sp.nu0, 2 / T10);
    const g = peakFactor(nuM, T10);
    const peakM = sp.meanMoment + g * sp.sigmaM, peakX = sp.meanThrust / md.stiffness + g * sp.sigmaX, peakA = g * sp.sigmaA;
    // Rotor harmonics: amplification at rated speed and the worst case within the speed range.
    const p1lo = tb.rpmMin / 60, p1hi = tb.rpmRated / 60;
    const worstIn = (lo, hi) => (f1 >= lo && f1 <= hi ? daf(1, zeta) : Math.max(daf(lo / f1, zeta), daf(hi / f1, zeta)));
    const daf1 = daf(p1hi / f1, zeta), daf3 = daf((3 * p1hi) / f1, zeta);
    const daf1w = worstIn(p1lo, p1hi), daf3w = worstIn(3 * p1lo, 3 * p1hi);
    const omR = tb.rotor.ratedRotorSpeed;
    const F1 = inp.imb * omR * omR, x1 = (F1 / md.stiffness) * daf1;
    // 3P thrust harmonic: entered, or from a blade-resolved unsteady BEM run at the wind speed of maximum thrust.
    let F3 = (inp.p3Frac / 100) * tb.thrustRated, p3Note = `3P thrust amplitude entered as ${inp.p3Frac} % of the maximum steady thrust.`;
    if (inp.p3Mode === 'bem') {
      const Trev = TWO_PI / omR, dtB = Trev / 96, nRev = 3;
      const ub = run(() => unsteadyBem({ rotor: tb.rotor, dt: dtB, nSteps: nRev * 96, wind: tb.windAtMaxThrust, omega: omR, pitchDeg: tb.controller.design.pitchFineDeg, rho: tb.rhoAir, shearExp: inp.shear, hubHeight: inp.hubHeight,
        towerShadow: { topDiameter: inp.towerDt, baseDiameter: inp.towerDb, height: tb.zTowerTop - inp.zInt, overhang: tb.rotor.overhang }, dynamicInflow: 'oye', dynamicStall: 'bl' }), 'Unsteady rotor aerodynamics');
      const h = harmonicAmplitudes(Array.from(ub.thrust.slice((nRev - 1) * 96, nRev * 96)), 6);
      const h3 = h.harmonics.find((q) => q.n === 3);
      F3 = h3 ? h3.amp : 0;
      p3Note = `3P thrust amplitude from unsteady BEM: ${fmt(F3 / 1e3)} kN (${fmt(F3 / Math.max(h.mean, 1) * 100)} % of the mean ${fmt(h.mean / 1e3)} kN); 6P ${fmt((h.harmonics.find((q) => q.n === 6)?.amp ?? 0) / 1e3)} kN.`;
    }
    const x3 = (F3 / md.stiffness) * daf3;
    const a1 = x1 * (TWO_PI * p1hi) ** 2, a3 = x3 * (TWO_PI * 3 * p1hi) ** 2;
    const accTot = peakA + a1 + a3;
    const rd = rayleighDamping({ f1, f2: str.modes[1].frequency, zeta1: md.zeta, zeta2: md.zeta });
    const checks = [
      check('Peak nacelle acceleration', accTot, inp.accLim, 'm/s²', 'Turbine supplier operational limit; DNV-ST-0126 (2021) §4.5', 'Guards the nacelle equipment and the drive train against vibration. Above the limit the turbine alarms and shuts down, losing production.'),
      check('Amplification of 1P over the speed range', daf1w, inp.dafLim, '-', 'DNV-ST-0126 (2021) §4.5.4', 'Guards against resonance with rotor imbalance at any operating speed. A large factor means small imbalances produce large tower sway and fatigue.'),
      check('Amplification of 3P over the speed range', daf3w, inp.dafLim, '-', 'DNV-ST-0126 (2021) §4.5.4', 'Guards against resonance with blade-passing loads during run-up and low-wind operation. Exceedance calls for a speed-exclusion window in the controller.'),
    ];
    const gov = worst(checks), ok = gov.uc <= 1;
    const [ff, sm, sw, sv] = thin([sp.f, sp.SM, sp.SMwind, sp.SMwave], 500);
    const fr = Array.from({ length: 160 }, (_, i) => 0.02 + (i * (Math.max(1, 3.4 * p1hi) - 0.02)) / 159);
    const dMax = Math.min(daf(1, zeta), 60);
    // Response against wind speed (coarser grid to stay fast).
    const vs = [5, 8, 11, 14, 18, 22, 25].filter((v) => v >= tb.rotor.cutIn && v <= tb.rotor.cutOut);
    const sweep = vs.map((v) => modalSpectra({ str, tb, inp, V: v, sigma1: ntmSigma({ Vhub: v, Iref: inp.Iref }), Hs: inp.Hs, Tp: inp.Tp, gamma: inp.gamma, nFreq: Math.max(100, Math.round(inp.nFreq / 4)) }));
    const sm0 = steelMass({ ...inp, Lp: 0 }, tb, steel.rho);
    // Wave-resonance share: energy near f1 against total wave response.
    const zetaNeed = accTot > inp.accLim ? zeta * (accTot / inp.accLim) ** 2 : zeta;
    return {
      kpis: [
        { key: 'f1', label: 'First natural frequency', unit: 'Hz', value: f1, explain: 'Fore-aft bending frequency of the stick model on the foundation springs; all amplification factors are referred to it.' },
        { key: 'zetaAero', label: 'Aerodynamic damping ratio', unit: '%', value: sp.zetaAero * 100, explain: 'Damping of the fore-aft mode produced by the operating rotor (thrust slope dT/dV at frozen speed and pitch). Typically 3–8 % in production and zero when parked, which is why idling storms govern wave-induced fatigue.' },
        { key: 'zetaTot', label: 'Total damping ratio', unit: '%', value: zeta * 100, explain: 'Structural, soil and hydrodynamic damping plus aerodynamic damping. Response at resonance is inversely proportional to it.' },
        { key: 'sigmaM', label: 'Std. deviation of mudline moment', unit: 'MN·m', value: sp.sigmaM / 1e6, explain: 'Root-mean-square fluctuation of the overturning moment at the mudline from turbulence and waves together. It drives fatigue; compare the wind and wave parts in the table.' },
        { key: 'meanM', label: 'Mean mudline moment', unit: 'MN·m', value: sp.meanMoment / 1e6, explain: 'Steady overturning moment from the mean rotor thrust times the hub height above the mudline.' },
        { key: 'peakM', label: 'Expected 10-minute peak mudline moment', unit: 'MN·m', value: peakM / 1e6, explain: 'Mean plus peak factor times the standard deviation (Davenport, Gaussian process). A characteristic operating extreme, before load factors.' },
        { key: 'sigmaX', label: 'Std. deviation of hub displacement', unit: 'mm', value: sp.sigmaX * 1e3, explain: 'Fluctuating fore-aft displacement at hub height.' },
        { key: 'peakX', label: 'Expected peak hub displacement', unit: 'm', value: peakX, explain: 'Static deflection under mean thrust plus the dynamic peak. Check against blade–tower clearance and cable hang-off assumptions.' },
        { key: 'peakA', label: 'Peak nacelle acceleration', unit: 'm/s²', value: accTot, limit: inp.accLim, sense: 'max', explain: 'Stochastic peak plus the 1P and 3P harmonic amplitudes at rated speed. Turbines alarm at roughly 0.15–0.3 g.' },
        { key: 'daf1P', label: 'Amplification of 1P at rated speed', unit: '-', value: daf1, explain: 'Dynamic amplification factor of a once-per-revolution load at rated rotor speed. Close to 1 is quasi-static; large values mean resonance.' },
        { key: 'daf3P', label: 'Amplification of 3P at rated speed', unit: '-', value: daf3, explain: 'Dynamic amplification of the blade-passing load at rated speed. Below 1 the structure isolates the excitation (soft side of 3P).' },
        { key: 'dafWave', label: 'Amplification at the wave peak frequency', unit: '-', value: daf(1 / inp.Tp / f1, zeta), explain: 'Amplification of loads at the spectral peak of the sea state. It grows as the wave period approaches the natural period, typical for large monopiles in short steep seas.' },
        { key: 'x1P', label: 'Hub displacement amplitude at 1P', unit: 'mm', value: x1 * 1e3, explain: 'Sway caused by the residual rotor imbalance at rated speed.' },
        { key: 'x3P', label: 'Hub displacement amplitude at 3P', unit: 'mm', value: x3 * 1e3, explain: 'Sway caused by the blade-passing thrust variation at rated speed.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Mudline moment response spectrum', x: { label: 'Frequency', unit: 'Hz' }, y: { label: 'Spectral density', unit: '(MN·m)²/Hz' },
          series: [{ name: 'Total', x: ff, y: sm.map((v) => v / 1e12) }, { name: 'Turbulence', x: ff, y: sw.map((v) => v / 1e12), dash: true }, { name: 'Waves', x: ff, y: sv.map((v) => v / 1e12), dash: true }],
          explain: 'Where the fluctuating moment comes from: low-frequency turbulence, the wave band around 1/Tp and the resonance peak at f1. A dominant resonance peak means damping controls fatigue.' },
        { type: 'line', title: 'Dynamic amplification and excitation ranges', x: { label: 'Excitation frequency', unit: 'Hz' }, y: { label: 'Amplification factor', unit: '-' },
          series: [{ name: `DAF (ζ = ${fmt(zeta * 100)} %)`, x: fr, y: fr.map((f) => Math.min(daf(f / f1, zeta), 60)) }, { name: '1P range', x: [p1lo, p1lo, p1hi, p1hi], y: [0, dMax, dMax, 0], dash: true }, { name: '3P range', x: [3 * p1lo, 3 * p1lo, 3 * p1hi, 3 * p1hi], y: [0, dMax, dMax, 0], dash: true }, { name: 'Wave peak', x: [1 / inp.Tp], y: [daf(1 / inp.Tp / f1, zeta)] }],
          explain: 'Amplification of a harmonic load against its frequency, with the rotor-speed (1P) and blade-passing (3P) ranges boxed. The resonance peak must fall in the gap between the boxes.' },
        { type: 'line', title: 'Response against wind speed', x: { label: 'Hub wind speed', unit: 'm/s' }, y: { label: 'Mudline moment', unit: 'MN·m' },
          series: [{ name: 'Mean', x: vs, y: sweep.map((s) => s.meanMoment / 1e6) }, { name: 'Standard deviation', x: vs, y: sweep.map((s) => s.sigmaM / 1e6) }, { name: 'Aerodynamic damping ×10, %', x: vs, y: sweep.map((s) => s.zetaAero * 1000), dash: true }],
          explain: 'Mean and fluctuating moment across the operating range in the same sea state. The mean peaks at rated wind; the fluctuation depends on turbulence and on how much aerodynamic damping is available.' },
      ],
      tables: [
        { title: 'Response and damping budget', columns: [{ key: 'k', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['Turbulence standard deviation σ1', sigma1, 'm/s'], ['Thrust slope, frozen control', sp.slopes.frozen / 1e3, 'kN·s/m'], ['Thrust slope along the schedule', sp.slopes.schedule / 1e3, 'kN·s/m'],
            ['Structural + soil + hydrodynamic damping', md.zeta * 100, '%'], ['Aerodynamic damping', sp.zetaAero * 100, '%'], ['Modal mass (unit hub displacement)', md.mass / 1e3, 't'], ['Modal stiffness', md.stiffness / 1e6, 'MN/m'],
            ['σ of mudline moment, turbulence', sp.sigmaMwind / 1e6, 'MN·m'], ['σ of mudline moment, waves', sp.sigmaMwave / 1e6, 'MN·m'], ['σ of external shear at mudline', sp.sigmaV / 1e6, 'MN'], ['Mean-crossing frequency of the moment', sp.nu0, 'Hz'], ['Peak factor (10 min)', g, '-'],
            ['Stochastic peak acceleration', peakA, 'm/s²'], ['1P acceleration amplitude', a1, 'm/s²'], ['3P thrust amplitude', F3 / 1e3, 'kN'], ['3P acceleration amplitude', a3, 'm/s²'], ['Worst 1P amplification in speed range', daf1w, '-'], ['Worst 3P amplification in speed range', daf3w, '-'],
            ['Rayleigh mass coefficient α', rd.alpha, '1/s'], ['Rayleigh stiffness coefficient β', rd.beta, 's']].map(([k, v, u]) => ({ k, v: round(v), u })),
          explain: 'Intermediate quantities of the response calculation, including the Rayleigh coefficients that reproduce the non-aerodynamic damping at the first two modes in a time-domain model.' },
      ],
      notes: [`${tb.name}. ${sp.slopes.op.parked ? 'Rotor parked: no aerodynamic damping.' : `Operating at ${fmt(inp.Vhub)} m/s with ${fmt(sp.zetaAero * 100)} % aerodynamic damping.`}`,
        `${sp.sigmaMwind > sp.sigmaMwave ? 'Turbulence dominates' : 'Waves dominate'} the fluctuating mudline moment (${fmt(sp.sigmaMwind / 1e6)} against ${fmt(sp.sigmaMwave / 1e6)} MN·m).`, `${gov.name} governs (utilisation ${fmt(gov.uc)}).`, p3Note],
      advice: [
        ok ? adv('engineering', sev(gov.uc), `Dynamic response is acceptable: peak nacelle acceleration ${fmt(accTot)} m/s² against ${inp.accLim}, worst 1P / 3P amplification ${fmt(daf1w)} / ${fmt(daf3w)} against ${inp.dafLim}.`,
          `Aerodynamic damping here comes from a generic rotor model: confirm it with the supplier. Re-run with the rotor parked above cut-out and in the idling 1-year storm: without the ${fmt(sp.zetaAero * 100)} % aerodynamic damping the resonant wave response grows by about ${fmt(Math.sqrt(zeta / md.zeta))} times.`,
          'DNV-ST-0126 (2021) §4.5; IEC 61400-3-1:2019 §7.5',
          `Fatigue-driving moment fluctuation is ${fmt(sp.sigmaM / 1e6)} MN·m rms; a tower damper adding 1 % damping would cut the resonant part by about ${fmt((1 - Math.sqrt(zeta / (zeta + 0.01))) * 100)} % and avoid adding steel (1 % more wall is ${fmt(sm0.total / 1e5)} t ≈ ${fmt(sm0.total / 1e5 * 1.9)} t CO₂e at the 1.9 t/t screening factor).`)
          : adv('engineering', 'fail', `${gov.name}: ${fmt(gov.demand)} against ${fmt(gov.capacity)} ${gov.unit}.`,
            gov.name.includes('acceleration') ? `Raise total damping from ${fmt(zeta * 100)} % to about ${fmt(zetaNeed * 100)} % (tuned mass or sloshing damper, active tower damping in the pitch controller) or move the natural frequency away from the wave band.` : `Move f1 = ${fmt(f1)} Hz out of the excitation range (stiffen or soften the support structure) or agree a rotor-speed exclusion window around ${fmt(f1 * 60 / (gov.name.includes('3P') ? 3 : 1))} rpm.`,
            'DNV-ST-0126 (2021) §4.5.4; IEC 61400-3-1:2019 §7.5',
            `A damper or controller change adds almost no material; solving the same problem with stiffness costs steel — ${steelCo2Text(sm0.total / 1e4)} for a 10 % wall increase.`),
        adv('economic', ok ? 'info' : 'warn', `Fluctuating mudline moment is ${fmt(sp.sigmaM / sp.meanMoment * 100)} % of the mean; damping totals ${fmt(zeta * 100)} % of which ${fmt(sp.zetaAero * 100)} % depends on the turbine running.`,
          'Value availability in storms: every hour the turbine idles with waves near the natural period consumes fatigue life several times faster than production, so high availability and a tower-damping controller mode protect the asset as well as the yield.',
          'Whole-life costing per ISO 15686-5; fatigue life consumption as a condition-based operating cost (DNV-ST-0262 lifetime extension practice)'),
      ],
    };
  },
  verify: [
    { name: 'Amplification of 1P and 3P for a uniform cantilever tower with the rotor parked: DAF = 1/√((1−r²)² + (2ζr)²)', source: 'Hand calculation: f1 = 0.4794 Hz (Blevins Table 8-1 for D = 6 m, t = 60 mm, L = 112.6 m, no top mass); rated speed 12.1 rpm → 1P = 0.2017 Hz, r = 0.4207 → DAF = 1.2149; 3P = 0.605 Hz, r = 1.2619, ζ = 1.2 % → DAF = 1.6858 (Chopra, Dynamics of Structures, Eq. 3.2.11)',
      inputs: { fnd: 'fixed', D: 6, t: 60, towerDb: 6, towerDt: 6, towerTb: 60, towerTt: 60, towerMassFactor: 1, Ca: 0, flooded: false, rnaT: 0.001, depth: 25, hubHeight: 90, ratedMW: 5, Vhub: 30, zetaS: 1.2 },
      expect: [{ kpi: 'daf1P', value: 1.2149, tol: 0.005 }, { kpi: 'daf3P', value: 1.6858, tol: 0.01 }, { kpi: 'zetaAero', value: 0, abs: 1e-9 }] },
  ],
};
