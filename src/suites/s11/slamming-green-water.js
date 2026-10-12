// Slamming, green water and sloshing screening: relative motion at a hull section from the motion model, Ochi
// slamming and deck-wetness probabilities, water-entry pressures (von Kármán, Wagner and the DNV space-average
// coefficient), dam-break green-water flow against a deck structure, and tank sloshing periods against the motion
// periods.
import { relativeMotion } from '../../engine/hydro/floater.js';
import { wedgeWaterEntry, slammingPressureCoefficient, slammingProbability, deckWetnessProbability, damBreakDeckFlow, sloshingPeriods, sloshingResonanceCheck } from '../../engine/hydro/impact.js';
import { hullInputs, hullFrom, motionInputs, motionFrom, seaInputs, seaFrom, spectrum, omegaGrid, trapz,
  G, DEG, LIMIT_HYDRO, LIMIT_CFD, need, run, round, rnd, fmt, adv, check, checkMin, sev, worst, linspace, steelCo2Text } from './_common.js';

// DNV-RP-C205 §8 gives the space-average coefficient for deadrise < 4° (2π) and ≥ 15° (2.5/tan^1.1 β) only; the gap is
// bridged linearly between the two end values (screening assumption, stated in the notes).
function averageCoefficient(betaDeg) {
  if (betaDeg < 4) return slammingPressureCoefficient({ shape: 'flat', deadriseDeg: betaDeg });
  if (betaDeg >= 15) return slammingPressureCoefficient({ shape: 'wedge', deadriseDeg: betaDeg });
  const lo = 2 * Math.PI; const hi = slammingPressureCoefficient({ shape: 'wedge', deadriseDeg: 15 });
  return lo + (hi - lo) * (betaDeg - 4) / 11;
}

export default {
  id: 's11-slamming-green-water',
  suite: 11,
  title: 'Slamming, green water and sloshing',
  summary: 'How often the bow section slams and ships green water in a sea state, the water-entry pressures, the green-water load on a deck structure behind the bow, and whether a partly filled tank sloshes in resonance with the motions.',
  covers: [13],
  fidelity: 'screening',
  standards: ['DNV-RP-C205 (2019) §8 (slamming, wave impact)', 'DNV-OS-C102 (2021) Ch.2 Sec.4 (green water, bow impact, sloshing on ship-shaped units)', 'DNV-CG-0158 (2021) (sloshing analysis of tanks)', 'IACS UR S11 / CSR-H (2022) Pt 1 Ch 4 Sec 5 (bottom and bow impact pressure concepts)'],
  models: ['Relative vertical motion and velocity spectra at a hull section', 'Ochi slamming criterion (emergence and threshold velocity)', 'Deck-wetness probability (Rayleigh exceedance of freeboard)', 'von Kármán and Wagner wedge water entry', 'DNV space-average slamming pressure coefficient', 'Dam-break (Ritter) shallow-water deck flow', 'Linear sloshing natural periods of rectangular and cylindrical tanks'],
  limits: `Slamming and green water are analytical screening: rigid two-dimensional wedge entry at constant velocity without air cushion, hydroelasticity or three-dimensional relief; green water is the frictionless dam-break analogy on a dry horizontal deck; sloshing is the linear natural period only, with no impact pressure. Bow swell-up and wave asymmetry enter through one user factor. ${LIMIT_HYDRO} ${LIMIT_CFD}`,
  searchTerms: 'green water FPSO bow slamming Wagner water entry sloshing',
  inputs: [
    ...hullInputs({ type: 'ship', only: ['ship', 'barge'], withMesh: false, draft: 19, kg: 14 }),
    ...seaInputs({ Hs: 12, Tp: 13.5, ret: '100', heading: 180 }),
    ...motionInputs(),
    { key: 'xSec', label: 'Section forward of midship', unit: 'm', type: 'number', default: 122, min: -250, max: 250, group: 'Exposed section', role: 'geometry', help: 'Longitudinal position of the section checked for slamming and deck wetness (the forward perpendicular region governs for head seas).' },
    { key: 'swellUp', label: 'Swell-up / wave asymmetry factor', unit: '-', type: 'number', default: 1.2, min: 1, max: 2, group: 'Exposed section', role: 'bc', help: 'Scales the incident wave elevation at the section for the bow wave, diffraction and crest asymmetry (1.0 = undisturbed linear wave; 1.1–1.3 without model tests).' },
    { key: 'slamDepth', label: 'Depth of the slamming surface below still water (0 = keel at the draft)', unit: 'm', type: 'number', default: 0, min: -30, max: 60, group: 'Exposed section', role: 'geometry', help: 'Bottom slamming: the local draft. Bow-flare slamming: enter the height of the flare knuckle above the waterline as a negative number.' },
    { key: 'deadriseDeg', label: 'Deadrise / flare angle to the water surface', unit: 'deg', type: 'number', default: 12, min: 1, max: 60, group: 'Exposed section', role: 'geometry', design: true, help: 'Angle between the hull surface and the water at impact: a few degrees for a flat bottom, 20–40° for a flared bow.' },
    { key: 'halfBeam', label: 'Half-breadth of the impacting surface', unit: 'm', type: 'number', default: 6, min: 0.3, max: 40, group: 'Exposed section', role: 'geometry', help: 'Width from the first contact point to the chine or knuckle; the impact force stops growing when it is wetted.' },
    { key: 'vThreshold', label: 'Threshold relative velocity (0 = Ochi 0.093·√(gL))', unit: 'm/s', type: 'number', default: 0, min: 0, max: 15, group: 'Exposed section', role: 'criteria', help: 'Re-entry velocity below which the impact is not counted as a slam.' },
    { key: 'bulwark', label: 'Bulwark height above the deck edge', unit: 'm', type: 'number', default: 0, min: 0, max: 6, group: 'Green water', role: 'geometry', design: true, help: 'A solid bulwark adds to the freeboard seen by the wave.' },
    { key: 'h0Manual', label: 'Water height above deck at the edge (0 = from the relative-motion extreme)', unit: 'm', type: 'number', default: 0, min: 0, max: 20, group: 'Green water', role: 'bc', help: 'Leave 0 to use the most probable maximum freeboard exceedance of the sea state; enter a model-test value to override.' },
    { key: 'xStruct', label: 'Distance from the deck edge to the structure', unit: 'm', type: 'number', default: 15, min: 0.5, max: 150, group: 'Green water', role: 'geometry', help: 'Run of open deck between the bow (or side) and the first exposed structure or equipment.' },
    { key: 'wStruct', label: 'Exposed width of the structure', unit: 'm', type: 'number', default: 12, min: 0.2, max: 80, group: 'Green water', role: 'geometry', help: 'Width normal to the flow.' },
    { key: 'pDeckAllow', label: 'Allowable pressure on the deck structure', unit: 'kPa', type: 'number', default: 150, min: 5, max: 2000, group: 'Acceptance criteria', role: 'criteria', help: 'Design lateral pressure of the exposed front wall or protection (green-water breakwater).' },
    { key: 'pSlamAllow', label: 'Allowable slamming pressure on the plating', unit: 'kPa', type: 'number', default: 1200, min: 20, max: 6000, group: 'Acceptance criteria', role: 'criteria', help: 'Space-average impact pressure the local panel is designed for (from the local strength module or class rules).' },
    { key: 'pSlamMax', label: 'Allowable slamming probability per wave cycle', unit: '-', type: 'number', default: 0.03, min: 0.001, max: 0.5, group: 'Acceptance criteria', role: 'criteria', help: 'Ochi & Motter operability criterion: 3 slams per 100 cycles.' },
    { key: 'pWetMax', label: 'Allowable deck-wetness probability per wave cycle', unit: '-', type: 'number', default: 0.05, min: 0.001, max: 0.5, group: 'Acceptance criteria', role: 'criteria', help: 'Customary seakeeping criterion: 5 wetnesses per 100 cycles (lower for manned areas).' },
    { key: 'tankShape', label: 'Tank shape', type: 'select', options: [{ value: 'rectangular', label: 'Rectangular (prismatic)' }, { value: 'cylindrical', label: 'Upright cylinder' }], default: 'rectangular', group: 'Sloshing tank', role: 'geometry', help: 'Smooth tank without internal structure.' },
    { key: 'tankL', label: 'Free-surface length in the direction of motion (diameter for a cylinder)', unit: 'm', type: 'number', default: 38, min: 1, max: 120, group: 'Sloshing tank', role: 'geometry', design: true, help: 'Tank breadth for roll-induced sloshing, tank length for pitch-induced sloshing.' },
    { key: 'tankFill', label: 'Liquid depth', unit: 'm', type: 'number', default: 9, min: 0.1, max: 40, group: 'Sloshing tank', role: 'ic', help: 'Fill depth of the partly filled tank in this loading condition.' },
    { key: 'tankDir', label: 'Exciting motion', type: 'select', options: [{ value: 'roll', label: 'Roll (transverse sloshing)' }, { value: 'pitch', label: 'Pitch (longitudinal sloshing)' }], default: 'roll', group: 'Sloshing tank', role: 'bc', help: 'The motion whose natural period is compared with the first sloshing period.' },
    { key: 'sloshBand', label: 'Resonance band', unit: '-', type: 'number', default: 0.2, min: 0.05, max: 0.5, group: 'Acceptance criteria', role: 'criteria', help: 'Relative separation required between the sloshing period and the motion period (screening value, ±20 %).' },
  ],
  mesh: { key: 'nFreq', label: 'Frequency points', levels: [50, 100, 200], qoi: 'relSig' },
  compute(inp) {
    const hull = hullFrom(inp); const mot = motionFrom(inp, hull);
    const sea = seaFrom(inp); const om = omegaGrid(inp.nFreq); const dur = inp.durationH * 3600; const rho = hull.rho;
    need(Math.abs(inp.xSec) <= hull.Lc / 2 + 1e-9, `The section at ${fmt(inp.xSec)} m lies beyond the hull ends (±${fmt(hull.Lc / 2)} m).`);
    const r = mot.rao(om, inp.headingDeg, sea);
    const rel = relativeMotion({ rao: r, point: [inp.xSec, 0, 0], waveFactor: inp.swellUp });
    const S = spectrum(om, sea);
    const m0r = trapz(om, Array.from(rel.amp, (a, i) => a * a * S[i]));
    const m0v = trapz(om, Array.from(rel.amp, (a, i) => (a * om[i]) ** 2 * S[i]));
    const calm = !(m0r > 1e-12);
    const Tz = calm ? 0 : 2 * Math.PI * Math.sqrt(m0r / m0v); const nCyc = calm ? 0 : dur / Tz;
    const relSig = Math.sqrt(Math.max(m0r, 0)); const relMpm = calm ? 0 : relSig * Math.sqrt(2 * Math.log(Math.max(nCyc, 2)));
    // --- slamming -------------------------------------------------------------------------------------------
    const d = inp.slamDepth === 0 ? hull.T : Math.abs(inp.slamDepth);
    const vThr = inp.vThreshold > 0 ? inp.vThreshold : 0.093 * Math.sqrt(G * hull.Lc);
    const sl = calm ? { probability: 0, ratePerHour: 0, emergenceProbability: 0 } : run(() => slammingProbability({ m0r, m0v, d, vThreshold: vThr }), 'Slamming probability');
    const nSlam = sl.ratePerHour * inp.durationH;
    // Ochi & Motter (1973): slam velocities follow a truncated Rayleigh law, so the most probable extreme of N slams is
    // v = √(v_cr² + 2·m0v·ln N); with fewer than one expected slam the threshold velocity itself is the design value.
    const vSlam = calm ? 0 : Math.sqrt(vThr * vThr + 2 * m0v * Math.log(Math.max(nSlam, 1)));
    const Cpa = averageCoefficient(inp.deadriseDeg);
    const tb = Math.tan(inp.deadriseDeg * DEG);
    const tChine = vSlam > 0 ? inp.halfBeam * tb / ((Math.PI / 2) * vSlam) : 0;
    const we = wedgeWaterEntry({ deadriseDeg: inp.deadriseDeg, V: vSlam, t: tChine, halfBeam: inp.halfBeam, rho });
    const pAvg = 0.5 * rho * Cpa * vSlam * vSlam; const pPeak = we.wagner.peakPressure;
    const tt = linspace(0, Math.max(tChine, 1e-6), 40);
    const hist = tt.map((t) => wedgeWaterEntry({ deadriseDeg: inp.deadriseDeg, V: vSlam, t, rho }));
    // --- green water ----------------------------------------------------------------------------------------
    const freeboard = hull.depth - hull.T + inp.bulwark;
    const wet = calm ? { probability: 0, ratePerHour: 0 } : run(() => deckWetnessProbability({ m0r, m2r: m0v, freeboard }), 'Deck wetness');
    const h0 = inp.h0Manual > 0 ? inp.h0Manual : Math.max(relMpm - freeboard, 0);
    let gw = { front: 0, fluxMax: 0, pWall: 0, tArr: 0, hs: [], us: [], ts: [], hMax: 0 };
    if (h0 > 1e-6) {
      const c0 = Math.sqrt(G * h0); const tArr = inp.xStruct / (2 * c0);
      const ts = linspace(tArr, 12 * tArr, 60); const fl = ts.map((t) => damBreakDeckFlow({ h0, x: inp.xStruct, t, rho }));
      // the momentum flux ρ·h·u² at a fixed station peaks when x/t = √(g·h0)/2, where h = h0/4 and u = √(g·h0)
      const pk = damBreakDeckFlow({ h0, x: inp.xStruct, t: 2 * inp.xStruct / c0, rho });
      gw = { front: 2 * c0, fluxMax: pk.momentumFlux, pWall: rho * pk.u * pk.u, tArr, ts, hs: fl.map((f) => f.h), us: fl.map((f) => f.u), hMax: Math.max(...fl.map((f) => f.h)) };
    }
    const Fgw = gw.fluxMax * inp.wStruct;
    // --- sloshing -------------------------------------------------------------------------------------------
    need(inp.tankFill > 0, 'Enter a positive liquid depth for the sloshing tank.');
    const slosh = run(() => sloshingPeriods(inp.tankShape === 'rectangular' ? { shape: 'rectangular', L: inp.tankL, h: inp.tankFill } : { shape: 'cylindrical', R: inp.tankL / 2, h: inp.tankFill }), 'Sloshing periods');
    const Tmot = inp.tankDir === 'roll' ? mot.Tn[3] : mot.Tn[4];
    need(Tmot > 0, 'The natural period of the exciting motion could not be evaluated for this hull.');
    const res = sloshingResonanceCheck({ sloshingPeriod: slosh.periods[0], motionPeriods: [Tmot, inp.Tp], band: inp.sloshBand });
    const sep = Math.abs(res.ratios[0] - 1);
    const checks = [
      check('Slamming probability per cycle', sl.probability, inp.pSlamMax, '-', 'Ochi (1964); Ochi & Motter (1973) operability criterion (3 per 100 cycles)', 'Limits how often the section leaves the water and re-enters hard. Frequent slamming causes whipping stresses, local plate damage and forces a heading or draft change.'),
      check('Slamming pressure on the plating (space average)', pAvg / 1e3, inp.pSlamAllow, 'kPa', 'DNV-RP-C205 (2019) §8.6 (space-average slamming pressure)', 'Guards the bottom or flare panels against permanent set: the space-average impact pressure acts as a short lateral load on plate and stiffeners.'),
      check('Deck-wetness probability per cycle', wet.probability, inp.pWetMax, '-', 'Ochi (1964); seakeeping criterion (5 per 100 cycles); DNV-OS-C102 (2021) Ch.2 Sec.4', 'Limits how often solid water comes over the deck edge, which endangers personnel, piping and equipment on the exposed deck.'),
      check('Green-water pressure on the deck structure', gw.pWall / 1e3, inp.pDeckAllow, 'kPa', 'DNV-OS-C102 (2021) Ch.2 Sec.4 (green-water design pressures); Buchner (2002) dam-break analogy', 'Guards the first structure behind the deck edge against the horizontal jet of green water, which can flatten unprotected walls and equipment.'),
      checkMin(`Sloshing period separation from ${inp.tankDir}`, sep, inp.sloshBand, '-', 'DNV-CG-0158 (2021) (sloshing resonance screening); linear sloshing theory', 'Keeps the first sloshing period away from the motion period. In resonance the liquid surge grows until it hits the tank roof and bulkheads with impact pressures the linear model does not predict.'),
    ];
    const gov = worst(checks);
    const fbNeed = calm ? 0 : Math.sqrt(-2 * m0r * Math.log(inp.pWetMax));
    return {
      kpis: [
        { key: 'relSig', label: 'Relative motion at the section (standard deviation)', value: round(relSig, 5), unit: 'm', explain: 'Standard deviation of the wave elevation relative to the moving hull at the section, including the swell-up factor. It drives both emergence (slamming) and freeboard exceedance (green water).' },
        { key: 'relMpm', label: 'Most probable maximum relative motion', value: round(relMpm, 4), unit: 'm', explain: 'Largest relative elevation expected in the sea-state duration; compare with the local draft and the freeboard.' },
        { key: 'relVelSig', label: 'Relative velocity (standard deviation)', value: round(Math.sqrt(Math.max(m0v, 0)), 4), unit: 'm/s', explain: 'Standard deviation of the relative vertical velocity between hull and water; slams occur when re-entry exceeds the threshold velocity.' },
        { key: 'pSlam', label: 'Slamming probability per cycle', value: round(sl.probability, 5), unit: '-', limit: inp.pSlamMax, sense: 'max', explain: 'Probability that a wave cycle produces a slam (section emerges and re-enters above the threshold velocity). 0.03 is the usual operability limit.' },
        { key: 'slamsPerHour', label: 'Slams per hour', value: round(sl.ratePerHour, 4), unit: '1/h', explain: 'Expected number of slams per hour in this sea state and heading.' },
        { key: 'vSlam', label: 'Design slamming velocity', value: round(vSlam, 4), unit: 'm/s', explain: 'Most probable extreme re-entry velocity among the slams of the sea state (the threshold velocity when fewer than one slam is expected).' },
        { key: 'CpWagner', label: 'Wagner peak pressure coefficient', value: round(we.wagner.Cp, 5), unit: '-', explain: 'Peak pressure at the spray root divided by ½ρV²: π²/(4·tan²β). It rises steeply as the deadrise falls and is very localised.' },
        { key: 'pPeak', label: 'Wagner peak slamming pressure', value: round(pPeak / 1e3, 4), unit: 'kPa', explain: 'Local peak pressure travelling with the spray root; it acts on a narrow strip for milliseconds and is not a panel design pressure by itself.' },
        { key: 'pSlamAvg', label: 'Space-average slamming pressure', value: round(pAvg / 1e3, 4), unit: 'kPa', limit: inp.pSlamAllow, sense: 'max', explain: 'Pressure ½ρ·C_Pa·v² averaged over a panel, the value to carry into plate and stiffener checks.' },
        { key: 'FslamMax', label: 'Impact force per metre at chine wetting (Wagner)', value: round(we.wagner.force / 1e3, 4), unit: 'kN/m', explain: 'Vertical force per unit length of hull when the wetted width reaches the chine; multiply by the affected length for the global slam load.' },
        { key: 'freeboard', label: 'Effective freeboard', value: round(freeboard, 4), unit: 'm', explain: 'Deck edge plus bulwark above the still waterline at the section.' },
        { key: 'pWet', label: 'Deck-wetness probability per cycle', value: round(wet.probability, 6), unit: '-', limit: inp.pWetMax, sense: 'max', explain: 'Probability that the relative elevation exceeds the freeboard in a wave cycle.' },
        { key: 'wetPerHour', label: 'Green-water events per hour', value: round(wet.ratePerHour ?? 0, 4), unit: '1/h', explain: 'Expected number of freeboard exceedances per hour.' },
        { key: 'h0', label: 'Water height above the deck edge', value: round(h0, 4), unit: 'm', explain: 'Freeboard exceedance of the most probable maximum (or the entered value): the reservoir height of the dam-break flow.' },
        { key: 'frontSpeed', label: 'Green-water front speed', value: round(gw.front, 4), unit: 'm/s', explain: 'Speed of the water front running along the deck, 2·√(g·h0).' },
        { key: 'fluxMax', label: 'Peak momentum flux at the structure', value: round(gw.fluxMax / 1e3, 4), unit: 'kN/m', explain: 'Largest horizontal force per metre width on a wall that stops the flow, ρ·h·u².' },
        { key: 'pGreen', label: 'Green-water pressure on the structure', value: round(gw.pWall / 1e3, 4), unit: 'kPa', limit: inp.pDeckAllow, sense: 'max', explain: 'Momentum pressure ρ·u² over the wetted height at the time of peak momentum flux. Impact peaks measured in model tests can be higher.' },
        { key: 'Fgreen', label: 'Green-water force on the structure', value: round(Fgw / 1e3, 4), unit: 'kN', explain: 'Peak horizontal force on the exposed width.' },
        { key: 'sloshT1', label: 'First sloshing period', value: round(slosh.periods[0], 5), unit: 's', explain: 'Natural period of the first sloshing mode of the tank at this fill depth.' },
        { key: 'sloshRatio', label: 'Motion period / sloshing period', value: round(res.ratios[0], 4), unit: '-', explain: 'Ratio of the natural period of the exciting motion to the sloshing period; values near 1.0 mean resonance.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Wave and relative-motion spectra at the section', x: { label: 'Frequency', unit: 'rad/s' }, y: { label: 'Spectral density', unit: 'm²·s/rad' },
          series: [{ name: 'Wave', x: rnd(om), y: rnd(S) }, { name: 'Relative motion', x: rnd(om), y: rnd(Array.from(rel.amp, (a, i) => a * a * S[i])) }],
          explain: 'Where the relative-motion spectrum rises above the wave spectrum the hull moves against the wave (pitch out of phase), which is what exposes the bow to slamming and green water.' },
        { type: 'line', title: 'Water-entry force during the design slam', x: { label: 'Time after contact', unit: 'ms' }, y: { label: 'Force per unit length', unit: 'kN/m' },
          series: [{ name: 'Wagner (with pile-up)', x: rnd(tt.map((t) => t * 1e3)), y: rnd(hist.map((h) => h.wagner.force / 1e3)) }, { name: 'von Kármán', x: rnd(tt.map((t) => t * 1e3)), y: rnd(hist.map((h) => h.vonKarman.force / 1e3)) }],
          explain: 'Force grows linearly until the chine is wetted. Wagner includes the pile-up of water at the spray root and is (π/2)² times the von Kármán momentum value; the rise time shown is the load duration for a dynamic panel check.' },
        { type: 'line', title: 'Probability versus height of the exposed surface', x: { label: 'Draft or freeboard', unit: 'm' }, y: { label: 'Probability per cycle', unit: '-' },
          series: (() => { const hs = linspace(0, Math.max(3 * relMpm / 2, freeboard * 1.2, d * 1.2, 1), 50); return [{ name: 'Exceedance of height (wetness / emergence)', x: rnd(hs), y: rnd(hs.map((h) => (calm ? 0 : Math.exp(-h * h / (2 * m0r))))) }, { name: 'Freeboard', x: [round(freeboard), round(freeboard)], y: [0, 1], dash: true }, { name: 'Slamming depth', x: [round(d), round(d)], y: [0, 1], dash: true }]; })(),
          explain: 'Rayleigh probability that the relative motion exceeds a given height. Read the deck-wetness probability at the freeboard and the emergence probability at the slamming depth; a metre of extra freeboard moves you down this curve.' },
        { type: 'line', title: 'Green water at the structure', x: { label: 'Time after the water comes aboard', unit: 's' }, y: { label: 'Depth (m) / velocity (m/s)', unit: '-' },
          series: [{ name: 'Water depth', x: rnd(gw.ts || []), y: rnd(gw.hs) }, { name: 'Flow velocity', x: rnd(gw.ts || []), y: rnd(gw.us) }],
          explain: 'Dam-break flow at the structure position: the thin fast front arrives first, then the depth builds while the velocity falls. The load peaks in between, when ρ·h·u² is largest. Empty when no water comes aboard.' },
      ],
      tables: [
        { title: 'Sloshing modes', columns: [{ key: 'n', label: 'Mode' }, { key: 'T', label: 'Period', unit: 's' }, { key: 'rMot', label: `${inp.tankDir} period / T` }, { key: 'rTp', label: 'Tp / T' }],
          rows: slosh.periods.map((T, i) => ({ n: i + 1, T: round(T, 5), rMot: round(Tmot / T, 4), rTp: round(inp.Tp / T, 4) })),
          explain: 'First three sloshing periods with their ratio to the motion natural period and to the wave peak period. A ratio within the resonance band of 1.0 for mode 1 calls for a fill restriction or a detailed sloshing assessment.' },
      ],
      notes: [
        `Relative motion at x = ${fmt(inp.xSec)} m: σ ${fmt(relSig)} m, most probable maximum ${fmt(relMpm)} m against a slamming depth of ${fmt(d)} m and a freeboard of ${fmt(freeboard)} m.`,
        `Threshold velocity ${fmt(vThr)} m/s; ${fmt(nSlam)} slams expected in ${fmt(inp.durationH)} h. Space-average coefficient C_Pa = ${fmt(Cpa)}${inp.deadriseDeg >= 4 && inp.deadriseDeg < 15 ? ' (interpolated between the 4° and 15° expressions of DNV-RP-C205, which gives none in between)' : ''}.`,
        `Sloshing: T₁ = ${fmt(slosh.periods[0])} s against a ${inp.tankDir} natural period of ${fmt(Tmot)} s (ratio ${fmt(res.ratios[0])})${res.resonant[1] ? '; the wave peak period is also inside the resonance band' : ''}. Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at utilisation ${fmt(gov.uc)}: slamming probability ${fmt(sl.probability)}, deck-wetness probability ${fmt(wet.probability)}, green-water pressure ${fmt(gw.pWall / 1e3)} kPa.`,
          wet.probability > inp.pWetMax || gw.pWall / 1e3 > inp.pDeckAllow ? `Raise the effective freeboard to about ${fmt(fbNeed)} m (${fmt(Math.max(fbNeed - freeboard, 0))} m more by bulwark, forecastle or a lighter draft) to bring wetness to ${fmt(inp.pWetMax)}, or protect the structure with a breakwater designed for ${fmt(gw.pWall / 1e3)} kPa.`
            : sl.probability > inp.pSlamMax ? `Increase the forward draft to about ${fmt(Math.sqrt(Math.max(-2 * m0r * Math.log(inp.pSlamMax), 0)))} m in this sea state, or change heading; verify the panel for ${fmt(pAvg / 1e3)} kPa.`
              : sep < inp.sloshBand ? `Restrict the fill range of this tank: move the liquid depth so that T₁ leaves ${fmt(Tmot * (1 - inp.sloshBand))}–${fmt(Tmot * (1 + inp.sloshBand))} s, or fit a swash bulkhead that halves the free-surface length.`
                : 'All screening criteria are met. Confirm the swell-up factor and the green-water height by model test before fixing the breakwater and bow plating, and carry the space-average pressure into the local strength module.',
          'DNV-RP-C205 (2019) §8; DNV-OS-C102 (2021) Ch.2 Sec.4; Ochi & Motter (1973); Buchner (2002).'),
        adv('economic', 'info', `Each green-water or slamming exceedance in an operating sea state means a heading or draft restriction; a breakwater ${fmt(inp.wStruct)} m wide for ${fmt(gw.pWall / 1e3)} kPa is about ${fmt(inp.wStruct * 4 * 0.02 * 7.85 * 1.8)} t of stiffened plate (4 m high, 20 mm, 80 % stiffening allowance).`,
          'Compare the steel of a bulwark or breakwater with the lost payload of a reduced draft; decide before the topsides layout is frozen because moving equipment back from the bow is free at that stage.',
          'Whole-life costing per ISO 15686-5; model-test value-of-information against rule-based design pressures.',
          `${steelCo2Text(inp.wStruct * 4 * 0.02 * 7.85 * 1.8)} Avoided green-water damage also avoids offshore repair campaigns and the vessel fuel they burn.`),
      ],
    };
  },
  verify: [
    { // Wagner: Cp = π²/(4·tan²β); β = 10° → tan = 0.176327, Cp = 9.8696/(4·0.031091) = 79.36.
      // Dam break with h0 = 4 m: front speed 2·√(g·h0) = 2·√39.2266 = 12.526 m/s; peak momentum flux ρ·g·h0²/4
      // = 1025·9.80665·16/4 = 40 207 N/m = 40.21 kN/m; wall pressure ρ·g·h0 = 40.21 kPa.
      name: 'Wagner peak pressure coefficient and Ritter dam-break flow',
      source: 'Wagner (1932) ZAMM 12; Ritter (1892) — hand calculation',
      inputs: { deadriseDeg: 10, h0Manual: 4 },
      expect: [{ kpi: 'CpWagner', value: 79.36, tol: 1e-3 }, { kpi: 'frontSpeed', value: 12.526, tol: 1e-3 }, { kpi: 'fluxMax', value: 40.21, tol: 1e-3 }, { kpi: 'pGreen', value: 40.21, tol: 1e-3 }] },
    { // Rectangular tank L = 40 m, h = 10 m: k = π/40 = 0.07854, tanh(k·h) = tanh(0.7854) = 0.65579,
      // ω² = 9.80665·0.07854·0.65579 = 0.50510 → ω = 0.71070 rad/s, T = 8.841 s.
      name: 'First sloshing period of a rectangular tank',
      source: 'Lamb (1932) Hydrodynamics §257; Faltinsen & Timokha (2009) Eq. 4.20 — hand calculation',
      inputs: { tankL: 40, tankFill: 10 },
      expect: [{ kpi: 'sloshT1', value: 8.841, tol: 1e-3 }] },
    { name: 'Flat calm: no slamming, no deck wetness',
      source: 'Limiting case Hs = 0',
      inputs: { Hs: 0 },
      expect: [{ kpi: 'pSlam', value: 0, abs: 1e-12 }, { kpi: 'pWet', value: 0, abs: 1e-12 }, { kpi: 'h0', value: 0, abs: 1e-12 }] },
  ],
};
