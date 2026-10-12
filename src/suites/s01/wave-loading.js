// Wave, current and wind loading on the jacket: regular design wave by Airy / Stokes 5th / stream function,
// Morison drag + inertia on every wetted member stepped through the wave phase and around the compass.
import { fmt, round, adv, check, linspace, MORISON_CAVEAT, CONFIRM, BASIS_COST, steelCo2 } from './_common.js';
import { geomInputs, envInputs, meshInputs, headingInput, headingList, buildJacket, waveField, envLoader, windLoad, maxLegD, steelMass, testFrame } from './_model.js';
import { breakingLimits } from '../../engine/metocean/waves.js';
import { jonswap, piersonMoskowitz, jonswapGamma, jonswapPeriodRatios, maxWaveHeight } from '../../engine/metocean/spectra.js';
import { powerLawProfile, logProfile } from '../../engine/metocean/wind.js';
import { currentPowerLaw } from '../../engine/metocean/current.js';
import { keuleganCarpenter, dragCoefficientSteady, inertiaCoefficient } from '../../engine/hydro/morison.js';

const THEORIES = [['airy', 'Airy to still water'], ['wheeler', 'Airy + Wheeler'], ['stokes5', 'Stokes 5th'], ['stream', 'Stream function']];

export default {
  id: 's01-wave-loading', suite: 1, title: 'Wave, current and wind loading',
  summary: 'Design-wave kinematics, Morison base shear and overturning moment by heading and wave theory, drag–inertia split, load distribution over the water column, wind force and air gap.',
  covers: [12], fidelity: 'numerical',
  standards: ['ISO 19901-1:2015 §8, §A.8–A.9 (metocean actions)', 'ISO 19902:2020 §9.5 (hydrodynamic actions)', 'API RP 2A-WSD 22nd ed. §5.3.1–5.3.2', 'DNV-RP-C205 (2019) §3, §6'],
  models: ['Airy linear wave theory', 'Wheeler stretching', 'Stokes 5th-order wave theory (Fenton)', 'Stream-function wave theory', 'Morison equation (drag + inertia)', 'Keulegan–Carpenter and roughness-dependent coefficients (DNV-RP-C205)', 'Power-law current profile with blockage', 'Power-law and logarithmic wind profiles', 'JONSWAP and Pierson–Moskowitz spectra', 'Rayleigh maximum wave height'],
  limits: `One regular design wave per heading on a rigid structure; irregular-sea, wave-in-deck, slam and vortex-shedding loads are not computed. ${MORISON_CAVEAT} Hydrodynamic coefficients, kinematics and blockage factors are criteria to be confirmed against the governing edition; conductors and appurtenances enter only through the coefficients.`,
  searchTerms: 'Morison wave loading jacket stream function Stokes fifth base shear',
  inputs: [
    ...geomInputs(), ...envInputs({ heading: false }).slice(0, 5), headingInput(),
    { key: 'Hs', label: 'Significant wave height of the design sea state', unit: 'm', type: 'number', default: 9.7, min: 0.5, max: 22, group: 'Environment', role: 'environment', bind: 'site.Hs100', help: 'Used for the spectrum plot and the Rayleigh check of the design wave height (Hmax ≈ 1.86 Hs in a 3-hour storm).' },
    { key: 'stormH', label: 'Storm duration', unit: 'h', type: 'number', default: 3, min: 0.5, max: 12, group: 'Environment', role: 'environment' },
    ...envInputs({ heading: false }).slice(5),
    { key: 'rough', label: 'Surface roughness k/D', unit: '-', type: 'number', default: 0.02, min: 0, max: 0.1, group: 'Hydrodynamic coefficients', role: 'environment', help: 'Used only for the DNV-RP-C205 coefficient comparison: below 1e-4 smooth, above 1e-2 rough (marine growth).' },
    { key: 'windA', label: 'Topsides wind area', unit: 'm²', type: 'number', default: 600, min: 0, max: 2e4, group: 'Wind', role: 'bc', help: 'Projected area normal to the wind.' },
    { key: 'windCs', label: 'Wind shape coefficient', unit: '-', type: 'number', default: 1.0, min: 0.3, max: 2, group: 'Wind', role: 'criteria', help: `1.0 for the overall projected area of a platform (API RP 2A). ${CONFIRM}` },
    { key: 'windZ', label: 'Wind area centroid above still water', unit: 'm', type: 'number', default: 30, min: 5, max: 150, group: 'Wind', role: 'bc' },
    { key: 'deckEl', label: 'Underside of the lowest deck above still water', unit: 'm', type: 'number', default: 21, min: 3, max: 60, group: 'Wind', role: 'geometry', help: 'For the air-gap check against the crest of the design wave.' },
    { key: 'airGap', label: 'Required air gap', unit: 'm', type: 'number', default: 1.5, min: 0, max: 5, group: 'Wind', role: 'criteria', help: `1.5 m above the design crest with surge and tide in ISO 19902 / API RP 2A. ${CONFIRM}` },
    ...meshInputs({ nSub: 2 }),
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'FhWave' },
  meshTime: { key: 'nPhase', label: 'Wave phase steps', levels: [8, 16, 32], qoi: 'FhWave' },
  compute(inp, ctx) {
    // SI conversion: member sizes mm → m inside the shared model; forces reported in MN.
    const J = buildJacket(inp, ctx, { found: 'fixed' }), heads = headingList(inp.nDir), wf = waveField(inp), d = wf.d;
    const Fw = windLoad(inp), armW = inp.windZ - J.zMud;
    const dirs = heads.map((b) => { const ld = envLoader(J, inp, wf, b), sc = ld.scan(inp.nPhase); return { b, ld, sc, F: Math.abs(sc.FhMax), M: Math.abs(sc.MhMax) }; });
    const gov = dirs.reduce((a, r) => (r.F > a.F ? r : a)), tG = gov.sc.t[gov.sc.iF];
    const FhWave = gov.F, baseShear = FhWave + Fw, otm = Math.max(...dirs.map((r) => r.M)) + Fw * armW;
    // drag / inertia split and current-only load at the governing phase
    const fD = envLoader(J, inp, wf, gov.b, { Cm: 0 }).at(tG).Fh, fI = envLoader(J, inp, wf, gov.b, { Cd: 0 }).at(tG).Fh;
    const still = waveField({ ...inp, H: 0 }), fCur = inp.Uc > 0 ? envLoader(J, inp, still, gov.b).at(0).Fh : 0;
    const nCmp = Math.max(8, Math.round(inp.nPhase / 2));   // coarser stepping for the comparison runs
    const noCur = wf.H > 0 ? Math.abs(envLoader(J, inp, wf, gov.b, { Uc: 0 }).scan(nCmp).FhMax) : 0;
    // other wave theories at the governing heading
    const th = [];
    if (wf.H > 0) for (const [k, label] of THEORIES) {
      try { const w2 = k === wf.theory ? wf : waveField(inp, { theory: k, fallback: false }); th.push({ label, F: k === wf.theory ? FhWave : Math.abs(envLoader(J, inp, w2, gov.b).scan(nCmp).FhMax), crest: w2.crest, used: k === wf.theory }); } catch { /* theory outside its range for this wave */ }
    }
    // load per unit height at the governing phase
    const st = gov.ld.at(tG, true), nb = 14, z0 = J.zMud, z1 = Math.min(J.zTop, inp.surge + wf.crest), dz = (z1 - z0) / nb, bins = new Array(nb).fill(0);
    J.elems.forEach((el) => { const zm = 0.5 * (el.a[2] + el.b[2]), i = Math.min(nb - 1, Math.max(0, Math.floor((zm - z0) / dz))); bins[i] += st.elemF[el.e]; });
    // kinematics under the crest
    const zs = linspace(-d, wf.crest, 30), uW = zs.map((z) => (wf.H > 0 ? wf.kin(0, z, 0).u : 0)), uC = zs.map((z) => (inp.Uc > 0 ? currentPowerLaw(Math.min(0, z), { v0: inp.Uc, d }) : 0));
    const uSurf = (wf.H > 0 ? wf.kin(0, Math.min(0, wf.crest), 0).u : 0) * inp.kinFac, D = maxLegD(J) + (inp.mg > 0 ? 2 * inp.mg / 1e3 : 0);
    const KC = wf.H > 0 ? keuleganCarpenter(uSurf, wf.T, D) : 0, Re = (Math.abs(uSurf) + inp.Uc) * D / 1.19e-6, DL = wf.H > 0 ? D / wf.L : 0;
    const cds = dragCoefficientSteady(inp.rough), cmDnv = inertiaCoefficient(KC, cds).Cm;
    const br = wf.H > 0 ? breakingLimits({ T: wf.T, d }) : { governing: Infinity };
    // sea state: spectrum and most probable maximum wave
    const gam = jonswapGamma(inp.Hs, inp.T), Tz = inp.T * jonswapPeriodRatios(Math.min(gam, 7)).TzOverTp, mw = maxWaveHeight({ Hs: inp.Hs, Tz, duration: inp.stormH * 3600 });
    const om = linspace(0.5 * 2 * Math.PI / inp.T, 3.5 * 2 * Math.PI / inp.T, 80);
    const crestEl = inp.surge + wf.crest, checks = [
      check('Morison slenderness D/λ', DL, 0.2, '-', 'DNV-RP-C205 §6.1.2; ISO 19902 §9.5.2', 'Morison loading assumes the member does not disturb the wave. Above D/λ = 0.2 diffraction reduces the inertia load and the result here becomes conservative and inaccurate.'),
      check('Wave height against breaking limit', wf.H, br.governing, 'm', 'Miche (1944) / McCowan (1894) limits; ISO 19901-1 §A.8.4', 'A wave higher than the breaking limit cannot exist at this depth and period; regular-wave theories fail and a breaking-wave (slam) load model is needed instead.'),
      check('Air gap under the lowest deck', crestEl + inp.airGap, inp.deckEl, 'm', 'ISO 19902:2020 §6.3.3.2; API RP 2A-WSD §5.3.4.3 (1.5 m air gap)', 'Crest elevation with surge and tide plus the required clearance, against the deck underside. A wave reaching the deck multiplies the base shear and is the classic cause of jacket collapse in hurricanes.'),
    ];
    const steel = steelMass(J, inp), dragShare = Math.abs(fD) / (Math.abs(fD) + Math.abs(fI) || 1), ok = checks.every((c) => c.uc <= 1);
    return {
      kpis: [
        { key: 'baseShear', label: 'Total base shear (wave + current + wind)', unit: 'MN', value: baseShear / 1e6, explain: 'Largest horizontal environmental force over the headings and wave phases, unfactored. It is the reference load for member sizing, pile lateral design and the reserve strength ratio.' },
        { key: 'otm', label: 'Overturning moment at the mudline', unit: 'MN·m', value: otm / 1e6, explain: 'Largest moment of the environmental loads about the mudline centre. Divided by the leg spacing at the base it gives the push–pull on the piles.' },
        { key: 'FhWave', label: 'Hydrodynamic base shear (wave + current)', unit: 'MN', value: FhWave / 1e6, explain: 'Morison force on all wetted members at the governing wave phase and heading, without wind.' },
        { key: 'Fwind', label: 'Wind force on the topsides', unit: 'MN', value: Fw / 1e6, explain: 'Sustained wind force on the projected deck area with the height coefficient at the centroid elevation.' },
        { key: 'dragShare', label: 'Drag share of the peak load', unit: '-', value: dragShare, explain: 'Fraction of the peak hydrodynamic force that is drag. Near 1 the load grows with H² and current matters; near 0 it is inertia-dominated and linear in H.' },
        { key: 'currentEffect', label: 'Increase of base shear due to current', unit: '-', value: noCur > 0 ? FhWave / noCur - 1 : 0, explain: 'Relative increase of the peak wave force when the current is added. In drag-dominated loading a modest current adds a large share because of the (u + Uc)² term.' },
        { key: 'crestEl', label: 'Crest elevation above still water', unit: 'm', value: crestEl, explain: 'Wave crest of the selected theory plus surge and tide; compare with the deck underside for air gap.' },
        { key: 'waveLen', label: 'Wavelength', unit: 'm', value: wf.H > 0 ? wf.L : 0, explain: 'Wavelength of the design wave at this depth. When it is close to twice the leg spacing, loads on opposite legs cancel partly.' },
        { key: 'KC', label: 'Keulegan–Carpenter number at the surface (leg)', unit: '-', value: KC, explain: 'Orbital excursion over member diameter. Above about 30 the flow is drag dominated and the steady-flow drag coefficient applies; below 6 inertia dominates.' },
        { key: 'HmaxRayleigh', label: 'Most probable maximum wave of the sea state', unit: 'm', value: mw.mostProbable, explain: 'Rayleigh estimate Hs·√(½ ln N) for the storm duration. The design wave height should not be smaller than this.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Hydrodynamic base shear through the wave cycle', x: { label: 'Time', unit: 's' }, y: { label: 'Base shear', unit: 'MN' }, series: dirs.map((r) => ({ name: `${r.b}°`, x: r.sc.t.map((v) => round(v)), y: r.sc.Fh.map((v) => round(v / 1e6)) })),
          explain: 'Force along each heading as the wave passes; the crest is over the centre of the jacket at t = 0. The peak slightly ahead of the crest shows the inertia contribution; a smooth, well-resolved peak confirms the phase step.' },
        { type: 'polar', title: 'Base shear and overturning moment by heading', series: [{ name: 'Base shear (MN)', thetaDeg: dirs.map((r) => r.b), r: dirs.map((r) => round((r.F + Fw) / 1e6)) }, { name: 'OTM (10 MN·m)', thetaDeg: dirs.map((r) => r.b), r: dirs.map((r) => round((r.M + Fw * armW) / 1e7)) }],
          explain: 'Total load against the direction of wave travel. Diagonal approach often gives the largest load on a square jacket because more braces are exposed.' },
        { type: 'bar', title: 'Peak hydrodynamic base shear by wave theory', categories: th.map((t) => t.label + (t.used ? ' (used)' : '')), series: [{ name: 'Base shear (MN)', y: th.map((t) => round(t.F / 1e6)) }], y: { label: 'Base shear', unit: 'MN' },
          explain: 'Same wave and structure with each theory valid for this case. Airy loaded to still water misses the crest loading; the spread between the nonlinear theories indicates the theory uncertainty.' },
        { type: 'line', title: 'Load intensity over the water column at the governing phase', x: { label: 'Horizontal load per metre height', unit: 'kN/m' }, y: { label: 'Elevation', unit: 'm' }, series: [{ name: 'Wave + current', x: bins.map((v) => round(v / dz / 1e3)), y: bins.map((v, i) => round(z0 + (i + 0.5) * dz)) }],
          explain: 'Horizontal Morison load collected in elevation bands. Most of the load sits in the top third of the water column, which is why marine growth and appurtenances near the surface are so costly.' },
        { type: 'line', title: 'Horizontal particle velocity under the crest', x: { label: 'Velocity', unit: 'm/s' }, y: { label: 'Elevation above still water', unit: 'm' }, series: [{ name: `Wave (${wf.label})`, x: uW.map((v) => round(v)), y: zs.map((v) => round(v)) }, { name: 'Current', x: uC.map((v) => round(v)), y: zs.map((v) => round(v)) }],
          explain: 'Velocity profile of the selected wave theory and of the current. Drag follows the square of their sum, so the upper part of the profile controls the load.' },
        { type: 'line', title: 'Wave spectrum of the design sea state', x: { label: 'Angular frequency', unit: 'rad/s' }, y: { label: 'S(ω)', unit: 'm²·s/rad' }, series: [{ name: `JONSWAP γ = ${fmt(gam)}`, x: om.map((v) => round(v)), y: om.map((w) => round(jonswap(w, { Hs: inp.Hs, Tp: inp.T }))) }, { name: 'Pierson–Moskowitz', x: om.map((v) => round(v)), y: om.map((w) => round(piersonMoskowitz(w, inp.Hs, inp.T)), ), dash: true }],
          explain: 'Energy of the sea state over frequency with the peak at the wave period. Compare the peak with the natural frequency of the platform in the dynamic module.' },
        { type: 'line', title: 'Wind speed profile', x: { label: 'Mean wind speed', unit: 'm/s' }, y: { label: 'Elevation', unit: 'm' }, series: [{ name: 'Power law (α = 0.125)', x: linspace(5, 100, 20).map((z) => round(powerLawProfile(z, inp.Uw, 10, 0.125))), y: linspace(5, 100, 20) }, { name: 'Logarithmic (z0 = 0.2 mm)', x: linspace(5, 100, 20).map((z) => round(logProfile(z, inp.Uw))), y: linspace(5, 100, 20), dash: true }],
          explain: 'One-hour mean wind speed with height by the two usual profiles. The force uses the power-law height coefficient at the centroid of the wind area.' },
      ],
      tables: [
        { title: 'Loads by heading', columns: [{ key: 'b', label: 'Heading', unit: 'deg' }, { key: 'F', label: 'Hydrodynamic base shear', unit: 'MN' }, { key: 'M', label: 'Hydrodynamic OTM', unit: 'MN·m' }, { key: 'tF', label: 'Phase of max shear', unit: 's' }, { key: 'Fz', label: 'Vertical load at that phase', unit: 'MN' }, { key: 'Ft', label: 'Total base shear with wind', unit: 'MN' }],
          rows: dirs.map((r) => ({ b: r.b, F: round(r.F / 1e6), M: round(r.M / 1e6), tF: round(r.sc.t[r.sc.iF]), Fz: round(r.sc.Fz[r.sc.iF] / 1e6), Ft: round((r.F + Fw) / 1e6) })), explain: 'One row per wave heading with the peak wave-and-current load, when in the cycle it occurs and the total including wind. Use the largest row for global design and all rows for pile loads.' },
        { title: 'Load build-up at the governing phase', columns: [{ key: 'c', label: 'Component' }, { key: 'F', label: 'Force', unit: 'MN' }],
          rows: [{ c: 'Drag (wave + current)', F: round(fD / 1e6) }, { c: 'Inertia', F: round(fI / 1e6) }, { c: 'Current alone, still water', F: round(fCur / 1e6) }, { c: 'Wave without current (peak)', F: round(noCur / 1e6) }, { c: 'Wind on topsides', F: round(Fw / 1e6) }],
          explain: 'Components of the base shear. Drag and inertia are evaluated at the phase of maximum total force, so they add up to the hydrodynamic base shear.' },
        { title: 'Hydrodynamic coefficients', columns: [{ key: 'c', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'n', label: 'Note' }],
          rows: [{ c: 'Cd used', v: inp.Cd, n: 'input (API RP 2A: 0.65 smooth, 1.05 rough)' }, { c: 'Cm used', v: inp.Cm, n: 'input (API RP 2A: 1.6 smooth, 1.2 rough)' }, { c: 'Steady-flow Cd for the roughness', v: round(cds), n: 'DNV-RP-C205 §6.7.1.5' }, { c: 'Cm for the KC number', v: round(cmDnv), n: 'DNV-RP-C205 §6.9.1' }, { c: 'Reynolds number (leg, surface)', v: round(Re, 3), n: 'post-critical above about 1e6' }, { c: 'D/λ (leg)', v: round(DL), n: 'Morison valid below 0.2' }],
          explain: 'Coefficients applied and the DNV-RP-C205 values implied by roughness and Keulegan–Carpenter number, for comparison. A large difference is a reason to review the inputs.' },
      ],
      notes: [`Governing heading ${gov.b}°: ${fmt(FhWave / 1e6)} MN hydrodynamic + ${fmt(Fw / 1e6)} MN wind. Wave theory: ${wf.label}; H/Hb = ${fmt(wf.H / br.governing)}.`,
        `Drag provides ${fmt(dragShare * 100)} % of the peak; the current raises the peak wave force by ${fmt((noCur > 0 ? FhWave / noCur - 1 : 0) * 100)} %.`,
        inp.H < 0.95 * mw.mostProbable ? `The design wave height (${fmt(inp.H)} m) is below the most probable maximum of the sea state (${fmt(mw.mostProbable)} m): check the metocean basis.` : `The design wave height is consistent with the sea state (most probable maximum ${fmt(mw.mostProbable)} m in ${fmt(inp.stormH)} h).`],
      advice: [
        adv('engineering', ok ? 'ok' : 'fail', ok ? `Loading basis is consistent: base shear ${fmt(baseShear / 1e6)} MN, overturning moment ${fmt(otm / 1e6)} MN·m, air gap ${fmt(inp.deckEl - crestEl)} m against ${fmt(inp.airGap)} m required.` : `Loading basis check failed: ${checks.filter((c) => c.uc > 1).map((c) => c.name).join('; ')}.`,
          ok ? 'Carry the base shear and overturning moment into the in-place and pile analyses. Confirm Cd, Cm, the kinematics factor and the blockage factor against the governing edition, and add conductors, risers and boat landings through equivalent diameters before final design.' : (crestEl + inp.airGap > inp.deckEl ? `Raise the deck underside to at least ${fmt(crestEl + inp.airGap)} m above still water or reassess with a wave-in-deck load.` : 'Reduce the wave to the breaking limit or use a diffraction / breaking-wave load model.'),
          'ISO 19901-1:2015 §8; ISO 19902:2020 §9.5; API RP 2A-WSD 22nd ed. §5.3', `Wave load scales with projected diameter: every 100 mm of marine growth avoided on the splash-zone members cuts the drag load there by about ${fmt(200 / (maxLegD(J) * 1e3) * 100)} % on the legs, and steel with it. The jacket holds ${steelCo2(steel)}`),
        adv('economic', 'info', `The structure is ${dragShare > 0.6 ? 'drag' : 'inertia'} dominated; ${fmt(bins.slice(Math.floor(nb * 2 / 3)).reduce((s, v) => s + v, 0) / (st.Fh || 1) * 100)} % of the wave load acts on the top third of the water column.`,
          'Spend the metocean budget where it pays: a site-specific hindcast that trims the 100-year wave height by 5 % cuts drag loading by about 10 % and the steel that follows. Keep appurtenances out of the wave zone and plan marine-growth cleaning instead of designing for the thickest growth.',
          `${BASIS_COST}; value of metocean data per ISO 19901-1 §5`),
      ],
      quantities: { steel_t: round(steel) },
    };
  },
  verify: [
    // Inertia force on a vertical surface-piercing cylinder in an Airy wave integrated to still water:
    // F = ρ Cm (πD²/4)·(H/2)·g·tanh(kd). d = 80 m, T = 13.5 s → k = 0.023196 rad/m, tanh(kd) = 0.95227;
    // D = 1.5 m, Cm = 2, H = 2 m: F = 1025·2·1.76715·1·9.80665·0.95227 = 33 831 N = 0.03383 MN.
    { name: 'Vertical cylinder, Airy wave, inertia only: F = ρ·Cm·(πD²/4)·(H/2)·g·tanh(kd)', source: 'Morison, O\'Brien, Johnson & Schaaf (1950); Dean & Dalrymple (1991) §8.3 closed-form pile force',
      inputs: { geom: testFrame('column'), H: 2, T: 13.5, Uc: 0, Uw: 0, surge: 0, theory: 'airy', Cd: 0, Cm: 2, mg: 0, kinFac: 1, appF: 1, nDir: '1', nSub: 8, nPhase: 48, nInt: 3 },
      expect: [{ kpi: 'FhWave', value: 0.03383, tol: 0.005 }, { kpi: 'dragShare', value: 0, abs: 1e-9 }] },
  ],
};
