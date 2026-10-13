// Extreme storm (survival) condition: most probable maximum and fractile motions in the 100-year sea, a sweep of
// peak period along the contour, air gap or freeboard exceedance with wave-in-deck load, static wind heel and
// extreme accelerations against survival limits.
import { pointMotion, relativeMotion, airGap } from '../../engine/hydro/floater.js';
import { waveInDeckForce } from '../../engine/hydro/impact.js';
import { hullInputs, hullFrom, motionInputs, motionFrom, seaInputs, seaFrom, spectrum, omegaGrid, stats, pointInputs, dofScale,
  G, DEG, LIMIT_HYDRO, LIMIT_CFD, need, run, round, rnd, fmt, adv, check, sev, worst, linspace, steelCo2Text } from './_common.js';

export default {
  id: 's11-survival',
  suite: 11,
  title: 'Extreme storm and survival condition',
  summary: 'Extreme motions, accelerations, inclination and air gap of the unit in the 100-year storm, scanned over the peak-period range, with the wave-in-deck load if the air gap is lost.',
  covers: [11],
  fidelity: 'reduced-order',
  standards: ['DNV-OS-C103 (2020) Ch.2 Sec.2 (air gap of column-stabilised units)', 'DNV-RP-C205 (2019) §3.7 (contours), §7, §8.2 (air gap)', 'DNV-OS-C101 (2021) Ch.2 Sec.1 (ULS environmental loads)', 'ISO 19904-1 (2019) §6, §9', 'API RP 2A-WSD (22nd ed.) §C17.6.2 (wave-in-deck)', 'IMO MODU Code 2009 §3.2.3 (wind)'],
  models: ['Linear spectral response with Rayleigh maxima: most probable maximum and fractile extremes', 'Peak-period scan along the extreme sea-state contour', 'Relative wave elevation and air gap with an asymmetry factor', 'Wave-in-deck silhouette drag force', 'Static wind heel from GM', 'Short-term probability of air-gap exceedance'],
  limits: `${LIMIT_HYDRO} Extremes assume a linear Gaussian response; wave run-up along columns, diffraction amplification and nonlinear crest kinematics are represented only by the user asymmetry factor. Static wind heel is the small-angle GM value. Wave-in-deck load is a silhouette drag estimate. ${LIMIT_CFD}`,
  searchTerms: 'air gap extreme response survival condition semi-submersible 100-year storm',
  inputs: [
    ...hullInputs({ type: 'semi', withMesh: false, draft: 20, kg: 21.5 }),
    ...seaInputs({ Hs: 14, Tp: 16, ret: '100', heading: 0 }),
    { key: 'TpLow', label: 'Lowest peak period on the contour', unit: 's', type: 'number', default: 13, min: 4, max: 30, group: 'Sea state', role: 'environment', help: 'Short-period end of the 100-year contour at this Hs (steeper waves).' },
    { key: 'TpHigh', label: 'Highest peak period on the contour', unit: 's', type: 'number', default: 19, min: 5, max: 35, group: 'Sea state', role: 'environment', help: 'Long-period end of the contour (swell-like storms that excite heave).' },
    { key: 'Uw', label: 'Storm wind speed (1-hour mean at 10 m)', unit: 'm/s', type: 'number', default: 40, min: 0, max: 80, group: 'Wind', role: 'environment', bind: 'site.Uw100', help: 'Wind for the static heel.' },
    { key: 'AwSide', label: 'Windage area', unit: 'm²', type: 'number', default: 3200, min: 0, max: 40000, group: 'Wind', role: 'load', help: 'Projected area above water in the wind direction.' },
    { key: 'Cwind', label: 'Wind force coefficient', unit: '-', type: 'number', default: 1.1, min: 0.4, max: 2, group: 'Wind', role: 'load', help: 'Combined shape and height coefficient: a user input from the MODU Code tables or a wind-tunnel test.' },
    { key: 'zWind', label: 'Windage centroid above waterline', unit: 'm', type: 'number', default: 32, min: 1, max: 200, group: 'Wind', role: 'load', help: 'Height of the centre of wind pressure.' },
    ...motionInputs(),
    ...pointInputs({ x: 35, y: 35, z: 40, group: 'Deck point', what: 'lowest deck point checked for air gap (deck-box underside at a column, or the bow deck edge of a ship)' }),
    { key: 'waveFactor', label: 'Wave asymmetry / amplification factor', unit: '-', type: 'number', default: 1.2, min: 1, max: 2, group: 'Deck point', role: 'bc', help: 'Scales the incident wave elevation for crest asymmetry and local amplification near columns (1.0 = undisturbed linear wave; 1.2–1.3 is customary without model tests).' },
    { key: 'deckWidth', label: 'Deck width exposed to the wave', unit: 'm', type: 'number', default: 70, min: 1, max: 300, group: 'Deck point', role: 'geometry', help: 'Width of deck normal to the waves for the wave-in-deck estimate.' },
    { key: 'CdDeck', label: 'Wave-in-deck drag coefficient', unit: '-', type: 'number', default: 2.0, min: 1, max: 3, group: 'Deck point', role: 'load', help: 'API RP 2A silhouette coefficient: about 1.2–1.6 for bare decks, 2.0–2.5 for heavily equipped decks. User input.' },
    { key: 'airGapMin', label: 'Required minimum air gap / freeboard', unit: 'm', type: 'number', default: 0, min: 0, max: 10, group: 'Acceptance criteria', role: 'criteria', help: 'Zero means the wave must not touch the deck; some projects require a positive margin (commonly 1.5 m).' },
    { key: 'inclMax', label: 'Allowable maximum inclination', unit: 'deg', type: 'number', default: 15, min: 2, max: 40, group: 'Acceptance criteria', role: 'criteria', help: 'Survival limit on static wind heel plus dynamic roll or pitch.' },
    { key: 'accMax', label: 'Allowable maximum acceleration at the point', unit: 'm/s²', type: 'number', default: 3, min: 0.2, max: 20, group: 'Acceptance criteria', role: 'criteria', help: 'Structural or equipment limit on the extreme vertical acceleration.' },
    { key: 'heaveMax', label: 'Allowable maximum heave', unit: 'm', type: 'number', default: 12, min: 0.5, max: 40, group: 'Acceptance criteria', role: 'criteria', help: 'Riser tensioner or tendon limit in the survival condition.' },
    { key: 'fractile', label: 'Design fractile of the 3-hour maximum', unit: '-', type: 'number', default: 0.9, min: 0.37, max: 0.99, group: 'Acceptance criteria', role: 'criteria', help: '0.37 is the most probable maximum; DNV practice takes the 90 % fractile for ULS response.' },
  ],
  mesh: { key: 'nFreq', label: 'Frequency points', levels: [50, 100, 200], qoi: 'heaveMpm' },
  compute(inp) {
    need(inp.TpHigh > inp.TpLow, 'The highest peak period of the contour must exceed the lowest.');
    const hull = hullFrom(inp); const mot = motionFrom(inp, hull);
    const sea = seaFrom(inp); const om = omegaGrid(inp.nFreq); const dur = inp.durationH * 3600;
    const r = mot.rao(om, inp.headingDeg, sea);
    const zDeck = inp.pz - hull.T;
    need(zDeck > 0, `The deck point (${fmt(inp.pz)} m above keel) is below the waterline at ${fmt(hull.T)} m draft: raise the point or reduce the draft.`);
    const pt = [inp.px, inp.py, zDeck];
    const pm = pointMotion({ rao: r, point: pt });
    const rel = relativeMotion({ rao: r, point: pt, waveFactor: inp.waveFactor });
    const frac = (s) => (s.sigma > 0 ? s.sigma * Math.sqrt(2 * Math.log(s.nCycles / -Math.log(inp.fractile))) : 0);
    const evalSea = (Tp) => {
      const S = spectrum(om, { ...sea, Tp });
      const h = stats(om, r.amp[2], S, dur); const ro = stats(om, r.amp[3], S, dur); const pi = stats(om, r.amp[4], S, dur);
      const ac = stats(om, pm.accAmp[2], S, dur); const rl = stats(om, rel.amp, S, dur);
      return { Tp, S, h, ro, pi, ac, rl, gap: zDeck - frac(rl) };
    };
    const tps = [...new Set([...linspace(inp.TpLow, inp.TpHigh, 9), inp.Tp].map((v) => round(v, 6)))].sort((a, b) => a - b);
    const sweep = tps.map(evalSea);
    const d = evalSea(inp.Tp);
    const ag = d.rl.m0 > 0 ? run(() => airGap({ rao: r, point: pt, staticAirGap: zDeck, S: d.S, duration: dur, waveFactor: inp.waveFactor }), 'Air gap') : { minAirGap: zDeck, exceedancePerCycle: 0, mpmRelative: 0 };
    const wGap = sweep.reduce((a, s) => (s.gap < a.gap ? s : a), sweep[0]);
    const wHeave = sweep.reduce((a, s) => (frac(s.h) > frac(a.h) ? s : a), sweep[0]);
    const wRot = sweep.reduce((a, s) => (Math.max(frac(s.ro), frac(s.pi)) > Math.max(frac(a.ro), frac(a.pi)) ? s : a), sweep[0]);
    const wAcc = sweep.reduce((a, s) => (frac(s.ac) > frac(a.ac) ? s : a), sweep[0]);
    // static wind heel, small-angle: M = F·(zWind + T/2), heel = M/(Δ g GM)
    const Fw = 0.5 * 1.222 * inp.Cwind * inp.AwSide * inp.Uw ** 2;
    const GM = hull.type === 'tlp' ? null : Math.min(hull.hs.GMt, hull.hs.GMl);
    const Krot = hull.type === 'tlp' ? Math.min(mot.Kext[3][3], mot.Kext[4][4]) : hull.mass * G * GM;
    const heelWind = Krot > 0 ? Fw * (inp.zWind + hull.T / 2) / Krot / DEG : 0;
    const rotDyn = Math.max(frac(wRot.ro), frac(wRot.pi)) * dofScale(3);
    const incl = heelWind + rotDyn;
    const pEx = wGap.rl.sigma > 0 ? 1 - Math.exp(-wGap.rl.nCycles * Math.exp(-(zDeck * zDeck) / (2 * wGap.rl.m0))) : 0;
    let wid = { force: 0, inundation: 0 };
    if (wGap.gap < 0) {
      // crest particle velocity of the equivalent regular wave (deep water): u = π·H/T with H = 2·relative extreme, T = Tz
      const u = Math.PI * 2 * frac(wGap.rl) / Math.max(wGap.rl.Tz, 1);
      wid = run(() => waveInDeckForce({ crestElevation: frac(wGap.rl), deckBottom: zDeck, width: inp.deckWidth, u, Cd: inp.CdDeck, rho: hull.rho }), 'Wave-in-deck force');
    }
    const checks = [
      check('Air gap / freeboard at the deck point (extreme relative elevation + required margin)', frac(wGap.rl) + inp.airGapMin, zDeck, 'm', 'DNV-OS-C103 (2020) Ch.2 Sec.2 [3] (positive air gap or design for impact); DNV-RP-C205 (2019) §8.2',
        'Guards the deck box against wave impact: a lost air gap brings slamming and horizontal wave loads onto structure and equipment not designed for them.'),
      check('Maximum heave', frac(wHeave.h), inp.heaveMax, 'm', 'Riser tensioner / tendon stroke in survival (project criterion)', 'Guards the riser and tendon interfaces, which stroke out or go slack beyond their heave range.'),
      check('Maximum inclination (wind heel + dynamic)', incl, inp.inclMax, 'deg', 'DNV-OS-C101 (2021) Ch.2 Sec.1; equipment survival inclination (project criterion)', 'Guards topsides equipment, drilling and process systems and the stability range, all designed for a limiting inclination in the survival condition.'),
      check('Maximum vertical acceleration at the point', frac(wAcc.ac), inp.accMax, 'm/s²', 'DNV-OS-C101 (2021) Ch.2 Sec.1 (inertia loads on topsides)', 'Guards module supports and sea-fastening against inertia loads beyond their design acceleration.'),
    ];
    const gov = worst(checks);
    const drise = Math.max(0, inp.airGapMin - wGap.gap);
    return {
      kpis: [
        { key: 'heaveMpm', label: 'Most probable maximum heave (design Tp)', value: round(d.h.mpm, 4), unit: 'm', explain: 'Heave amplitude most likely to be the largest in the storm duration at the design peak period.' },
        { key: 'heaveExt', label: 'Design extreme heave (worst Tp)', value: round(frac(wHeave.h), 4), unit: 'm', limit: inp.heaveMax, sense: 'max', explain: `Heave amplitude at the chosen fractile of the 3-hour maximum, at the worst peak period of the contour (${fmt(wHeave.Tp)} s).` },
        { key: 'rollExt', label: 'Design extreme roll', value: Math.abs(frac(wRot.ro) * dofScale(3)) < 1e-6 ? 0 : round(frac(wRot.ro) * dofScale(3), 4), unit: 'deg', explain: 'Roll amplitude at the chosen fractile, worst peak period. In head or following seas it vanishes by symmetry: a value below 1e-6 deg is arithmetic round-off and is reported as 0.' },
        { key: 'pitchExt', label: 'Design extreme pitch', value: round(frac(wRot.pi) * dofScale(4), 4), unit: 'deg', explain: 'Pitch amplitude at the chosen fractile, worst peak period.' },
        { key: 'heelWind', label: 'Static wind heel', value: round(heelWind, 4), unit: 'deg', explain: 'Steady inclination from the storm wind on the windage area, resisted by the hydrostatic (or tendon) stiffness.' },
        { key: 'inclination', label: 'Maximum inclination', value: round(incl, 4), unit: 'deg', limit: inp.inclMax, sense: 'max', explain: 'Static wind heel plus extreme dynamic roll or pitch: the angle topsides equipment must survive.' },
        { key: 'accExt', label: 'Design extreme vertical acceleration at the point', value: round(frac(wAcc.ac), 4), unit: 'm/s²', limit: inp.accMax, sense: 'max', explain: 'Extreme vertical acceleration of the deck point; multiply by supported mass for the survival inertia load.' },
        { key: 'staticGap', label: 'Still-water air gap / freeboard', value: round(zDeck, 4), unit: 'm', explain: 'Height of the deck point above the still waterline at the survival draft.' },
        { key: 'relExt', label: 'Design extreme relative wave elevation', value: round(frac(wGap.rl), 4), unit: 'm', explain: 'Extreme upward wave elevation relative to the moving deck point, including the asymmetry factor.' },
        { key: 'minAirGap', label: 'Minimum air gap (worst Tp)', value: round(wGap.gap, 4), unit: 'm', limit: inp.airGapMin, sense: 'min', explain: `Still-water air gap minus the extreme relative elevation at the worst peak period (${fmt(wGap.Tp)} s). Negative means the wave reaches the deck.` },
        { key: 'minAirGapMpm', label: 'Minimum air gap, most probable (design Tp)', value: round(ag.minAirGap, 4), unit: 'm', explain: 'Air gap using the most probable maximum relative motion at the design peak period, for comparison with the fractile value.' },
        { key: 'pExceed', label: 'Probability of deck contact in the storm', value: round(pEx, 5), unit: '-', explain: 'Probability that at least one wave reaches the deck point during the storm, from Rayleigh statistics of the relative elevation.' },
        { key: 'deckForce', label: 'Wave-in-deck force', value: round(wid.force / 1e6, 5), unit: 'MN', explain: 'Horizontal force on the inundated deck silhouette when the air gap is negative (zero otherwise).' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Extreme responses along the contour', x: { label: 'Peak period', unit: 's' }, y: { label: 'Extreme response', unit: 'm, deg, m/s²' },
          series: [{ name: 'Heave (m)', x: tps, y: rnd(sweep.map((s) => frac(s.h))) }, { name: 'Pitch (deg)', x: tps, y: rnd(sweep.map((s) => frac(s.pi) * dofScale(4))) }, { name: 'Roll (deg)', x: tps, y: rnd(sweep.map((s) => frac(s.ro) * dofScale(3))) }, { name: 'Vertical acceleration (m/s²)', x: tps, y: rnd(sweep.map((s) => frac(s.ac))) }],
          explain: 'Design extremes for the same Hs at each peak period of the contour. The worst period differs by response: heave grows towards its natural period, accelerations towards short periods.' },
        { type: 'line', title: 'Air gap along the contour', x: { label: 'Peak period', unit: 's' }, y: { label: 'Height', unit: 'm' },
          series: [{ name: 'Minimum air gap', x: tps, y: rnd(sweep.map((s) => s.gap)) }, { name: 'Still-water air gap', x: [tps[0], tps[tps.length - 1]], y: [round(zDeck), round(zDeck)], dash: true }, { name: 'Required', x: [tps[0], tps[tps.length - 1]], y: [inp.airGapMin, inp.airGapMin], dash: true }],
          explain: 'Remaining clearance between the extreme wave crest and the deck point. Where the curve drops below the required line the deck is hit; steep short-period storms are usually worst.' },
        { type: 'line', title: 'Relative-elevation and heave spectra (design sea)', x: { label: 'Frequency', unit: 'rad/s' }, y: { label: 'Spectral density', unit: 'm²·s/rad' },
          series: [{ name: 'Wave', x: rnd(om), y: rnd(d.S) }, { name: 'Relative elevation at the deck point', x: rnd(om), y: rnd(rel.amp.map((a, i) => a * a * d.S[i])) }, { name: 'Heave', x: rnd(om), y: rnd(r.amp[2].map((a, i) => a * a * d.S[i])) }],
          explain: 'Where the unit follows the waves the relative-elevation spectrum falls below the wave spectrum; where it does not (short waves) or moves out of phase (near resonance) it rises above it.' },
      ],
      tables: [
        { title: 'Contour scan', columns: [{ key: 'Tp', label: 'Tp', unit: 's' }, { key: 'h', label: 'Heave', unit: 'm' }, { key: 'r', label: 'Roll', unit: 'deg' }, { key: 'p', label: 'Pitch', unit: 'deg' }, { key: 'a', label: 'Vert. acc.', unit: 'm/s²' }, { key: 'rel', label: 'Relative elevation', unit: 'm' }, { key: 'gap', label: 'Air gap', unit: 'm' }],
          rows: sweep.map((s) => ({ Tp: s.Tp, h: round(frac(s.h), 4), r: round(frac(s.ro) * dofScale(3), 4), p: round(frac(s.pi) * dofScale(4), 4), a: round(frac(s.ac), 4), rel: round(frac(s.rl), 4), gap: round(s.gap, 4) })),
          explain: `One row per peak period on the contour, all at Hs ${fmt(inp.Hs)} m and the ${fmt(100 * inp.fractile)} % fractile of the storm maximum.` },
      ],
      notes: [
        `Survival draft ${fmt(hull.T)} m: still-water air gap ${fmt(zDeck)} m, minimum ${fmt(wGap.gap)} m at Tp ${fmt(wGap.Tp)} s (asymmetry factor ${fmt(inp.waveFactor)}).`,
        `Extreme heave ${fmt(frac(wHeave.h))} m at Tp ${fmt(wHeave.Tp)} s; inclination ${fmt(incl)}° = wind heel ${fmt(heelWind)}° + dynamic ${fmt(rotDyn)}°.`,
        `Drag damping linearised at the design sea is used for the whole contour. Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at utilisation ${fmt(gov.uc)}; minimum air gap ${fmt(wGap.gap)} m, probability of deck contact ${fmt(pEx)}.`,
          wGap.gap < inp.airGapMin ? `Recover ${fmt(drise)} m of air gap: deballast to a shallower survival draft by that amount (check stability at the new draft), raise the deck, or design the exposed deck for the ${fmt(wid.force / 1e6)} MN wave-in-deck load and local slamming.`
            : `Clearance margin is ${fmt(wGap.gap - inp.airGapMin)} m. Confirm the asymmetry factor and run-up at columns by model test or diffraction analysis: they can add 20–30 % to the relative elevation.`,
          'DNV-OS-C103 (2020) Ch.2 Sec.2 [3]; DNV-RP-C205 (2019) §8.2 and §3.7 (contour method, fractile of the short-term extreme).'),
        adv('economic', 'info', `Air gap is bought with column height: one metre more on four columns is about ${fmt(hull.type === 'semi' ? 4 * Math.PI * inp.Dcol * 0.03 * 7.85 * 1.6 : 60)} t of steel and raises the centre of gravity of the whole deck.`,
          'Use the survival draft as the first lever (it is free in steel), then the asymmetry factor documented by model tests, and only then column height.',
          'Value-of-information for a model-test campaign versus added steel; AACE RP 18R-97 class 4 estimating.',
          `${steelCo2Text(hull.type === 'semi' ? 4 * Math.PI * inp.Dcol * 0.03 * 7.85 * 1.6 : 60)} per metre of extra column height; a deck designed for wave impact needs heavier plating over the whole exposed underside instead.`),
      ],
    };
  },
  verify: [
    { name: 'Flat calm: minimum air gap equals the still-water air gap', source: 'Limiting case Hs = 0 (hydrostatic clearance = point height − draft = 40 − 20 m)', inputs: { Hs: 0 },
      expect: [{ kpi: 'minAirGap', value: 20, tol: 1e-9 }, { kpi: 'heaveExt', value: 0, abs: 1e-9 }, { kpi: 'pExceed', value: 0, abs: 1e-12 }] },
    { // Small barge 20 × 8 × 0.8 m in Hs = 4 m, Tp = 20 s (γ = 1): heave follows the wave, σ = Hs/4 = 1.0 m, Tz ≈ Tp·√(6/12)
      // = 14.14 s, N = 10 800/14.14 = 764, most probable maximum = σ·√(2 ln N) = 1.0·3.644 = 3.64 m.
      name: 'Rayleigh most probable maximum of a wave-following body: σ √(2 ln N)',
      source: 'Longuet-Higgins (1952) / DNV-RP-C205 (2019) §3.5.11 — hand calculation',
      inputs: { hullType: 'barge', L: 20, B: 8, Dhull: 4, draft: 0.8, KG: 1.5, Hs: 4, Tp: 20, gamma: 1, headingDeg: 180, TpLow: 19, TpHigh: 21, px: 0, py: 0, pz: 4, nFreq: 200 },
      expect: [{ kpi: 'heaveMpm', value: 3.64, tol: 0.03 }] },
    { // Wind heel by hand: F = ½·1.222·1.0·1000·40² = 977.6 kN; arm = 20 + 5/2 = 22.5 m; Δ g GM = 10 250 000·9.80665·5.1667
      // = 519.35 MN·m/rad → heel = 21.996/519.35 = 0.04235 rad = 2.427°.
      name: 'Static wind heel M/(Δ g GM) of a box barge',
      source: 'Small-angle hydrostatics, hand calculation',
      inputs: { hullType: 'barge', L: 100, B: 20, Dhull: 12, draft: 5, KG: 4, Hs: 0, Uw: 40, AwSide: 1000, Cwind: 1, zWind: 20, pz: 12 },
      expect: [{ kpi: 'heelWind', value: 2.427, tol: 2e-3 }] },
  ],
};
