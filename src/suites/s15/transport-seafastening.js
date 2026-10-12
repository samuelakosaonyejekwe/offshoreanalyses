// Sea transportation: motion criteria → inertia forces at the cargo centre of gravity → support reactions,
// seafastening brace forces, grillage beam bending (frame finite elements) and voyage fatigue of the welds.
import { analyzeFrame } from '../../engine/structural/frame3d.js';
import { iSection } from '../../engine/structural/sections.js';
import { snCurve, snConstants, SN_CURVE_NAMES } from '../../engine/fatigue/sn.js';
import { minerWeibull } from '../../engine/fatigue/damage.js';
import { G, DEG, need, run, round, fmt, kpi, check, adv, bar, line, tbl, worst, chiEC3, tube, steelCo2Text, CONFIRM, N001 } from './_common.js';

export default {
  id: 's15-transport-seafastening', suite: 15, title: 'Transportation and sea-fastening',
  summary: 'Turns the voyage motion criteria into inertia forces at the cargo centre of gravity and checks support reactions, uplift, seafastening braces and welds, the grillage beam and voyage fatigue.',
  covers: [5], fidelity: 'numerical',
  standards: [`${N001} Sec. 11 (sea voyages: default motion criteria, load combinations, seafastening and grillage design)`, 'EN 1993-1-1 (2005) §6.2–6.3 (member resistance), EN 1993-1-8 §4.5 (fillet welds)', 'DNV-RP-C203 (2021) §2.4, §5.1 (S–N curves, closed-form Weibull damage)', 'IMO CSS Code Annex 13 (acceleration approach for non-standardised cargo)'],
  models: ['Transportation inertia forces (roll/pitch about the waterline plus heave, d\'Alembert)', 'Rigid-body equilibrium of cargo on discrete supports', 'Coulomb friction (optional credit)', 'Nonlinear seafastening restraint simplified to axial braces', 'Structural grillage analysis by 3-D beam finite elements', 'Column buckling (EN 1993-1-1 curves)', 'Fatigue under repeated vessel motions (Miner sum, Rayleigh-distributed ranges)'],
  limits: 'Quasi-static design-wave approach: roll and pitch are single harmonic motions about an axis in the waterplane, combined with heave but not with each other, and the cargo is rigid. Default motion criteria are recalled DNV-ST-N001 values exposed as inputs — replace them with the voyage motion analysis. Barge deck and frame strength are represented by one allowable reaction; local plate, stiffener and nonlinear contact (lift-off of cribbing) are outside the module.',
  publishes: ['install.seafasteningForce', 'install.transportAccT', 'install.transportAccL'],
  searchTerms: 'sea transportation seafastening grillage barge motion criteria inertia forces DNV-ST-N001',
  inputs: [
    { key: 'W', label: 'Cargo weight', unit: 't', type: 'number', default: 6000, min: 5, max: 80000, group: 'Cargo', role: 'load', bind: 'install.liftWeight', uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'fWeight', label: 'Weight contingency factor', unit: '-', type: 'number', default: 1.05, min: 1, max: 1.5, group: 'Cargo', role: 'criteria', help: CONFIRM },
    { key: 'zc', label: 'Cargo CoG height above the deck', unit: 'm', type: 'number', default: 14, min: 0.1, max: 150, group: 'Cargo', role: 'geometry', uncertain: { dist: 'normal', cov: 0.05 } },
    { key: 'zDeck', label: 'Deck height above the motion centre (waterline)', unit: 'm', type: 'number', default: 3, min: -10, max: 40, group: 'Cargo', role: 'geometry', help: 'Barge depth minus draft. The motion centre is taken at the waterline, as in the simplified criteria.' },
    { key: 'yc', label: 'Cargo CoG offset from the barge centreline', unit: 'm', type: 'number', default: 0, min: -60, max: 60, group: 'Cargo', role: 'ic' },
    { key: 'xc', label: 'Cargo CoG distance from the pitch centre (midship)', unit: 'm', type: 'number', default: 10, min: -200, max: 200, group: 'Cargo', role: 'ic', help: 'Initial condition of the stowage: where the cargo stands on the deck.' },
    { key: 'b', label: 'Transverse spacing of the cargo supports', unit: 'm', type: 'number', default: 20, min: 1, max: 100, group: 'Cargo', role: 'geometry', design: true },
    { key: 'a', label: 'Longitudinal spacing of the cargo supports', unit: 'm', type: 'number', default: 30, min: 1, max: 300, group: 'Cargo', role: 'geometry' },
    { key: 'rollDeg', label: 'Roll single amplitude', unit: 'deg', type: 'number', default: 20, min: 0, max: 45, group: 'Motion criteria', role: 'bc', bind: 'install.transportRollDeg', uncertain: { dist: 'normal', cov: 0.15 }, help: `Default for a large barge on an unrestricted voyage: 20° in 10 s. ${CONFIRM}` },
    { key: 'rollT', label: 'Roll period', unit: 's', type: 'number', default: 10, min: 3, max: 40, group: 'Motion criteria', role: 'bc', bind: 'install.transportRollT' },
    { key: 'pitchDeg', label: 'Pitch single amplitude', unit: 'deg', type: 'number', default: 12.5, min: 0, max: 30, group: 'Motion criteria', role: 'bc', bind: 'install.transportPitchDeg', help: `Default for a large barge: 12.5° in 10 s. ${CONFIRM}` },
    { key: 'pitchT', label: 'Pitch period', unit: 's', type: 'number', default: 10, min: 3, max: 40, group: 'Motion criteria', role: 'bc' },
    { key: 'heaveG', label: 'Heave acceleration', unit: 'g', type: 'number', default: 0.2, min: 0, max: 1, group: 'Motion criteria', role: 'bc', bind: 'install.transportHeaveG', help: `Default 0.2 g, combined with roll or pitch. ${CONFIRM}` },
    { key: 'Awind', label: 'Cargo windage area (transverse)', unit: 'm²', type: 'number', default: 900, min: 0, max: 30000, group: 'Wind', role: 'geometry' },
    { key: 'Vwind', label: 'Design wind speed (1-min mean at cargo height)', unit: 'm/s', type: 'number', default: 35, min: 0, max: 80, group: 'Wind', role: 'environment', bind: 'site.Uw1' },
    { key: 'Cs', label: 'Wind shape coefficient', unit: '-', type: 'number', default: 1.1, min: 0.3, max: 2.5, group: 'Wind', role: 'load' },
    { key: 'mu', label: 'Friction credited between cargo and grillage', unit: '-', type: 'number', default: 0, min: 0, max: 0.5, group: 'Seafastening', role: 'bc', help: 'Friction is normally not relied on for steel-to-steel seafastening; credit timber cribbing only where the standard allows.' },
    { key: 'mat', label: 'Steel of seafastening and grillage', type: 'select', source: 'materials:steel', default: 's355', group: 'Seafastening', role: 'material' },
    { key: 'nT', label: 'Transverse braces (total, both sides effective)', unit: '-', type: 'int', default: 16, min: 1, max: 80, group: 'Seafastening', role: 'bc', design: true, help: 'Braces restraining the cargo in the transverse direction. Count only those that can take the load in the sense considered if they are tension-only or compression-only.' },
    { key: 'nL', label: 'Longitudinal braces (total)', unit: '-', type: 'int', default: 10, min: 1, max: 80, group: 'Seafastening', role: 'bc', design: true },
    { key: 'braceDeg', label: 'Brace angle to the deck', unit: 'deg', type: 'number', default: 45, min: 0, max: 75, group: 'Seafastening', role: 'geometry' },
    { key: 'Dbr', label: 'Brace tube diameter', unit: 'mm', type: 'number', default: 610, min: 60, max: 2000, group: 'Seafastening', role: 'geometry', design: true },
    { key: 'tbr', label: 'Brace wall thickness', unit: 'mm', type: 'number', default: 25, min: 4, max: 80, group: 'Seafastening', role: 'geometry', design: true },
    { key: 'Lbr', label: 'Brace length', unit: 'm', type: 'number', default: 6, min: 0.5, max: 40, group: 'Seafastening', role: 'geometry' },
    { key: 'aw', label: 'Fillet weld throat at the brace ends', unit: 'mm', type: 'number', default: 20, min: 3, max: 40, group: 'Seafastening', role: 'geometry', design: true },
    { key: 'hG', label: 'Grillage beam depth', unit: 'mm', type: 'number', default: 2200, min: 200, max: 5000, group: 'Grillage', role: 'geometry', design: true },
    { key: 'bG', label: 'Grillage flange width', unit: 'mm', type: 'number', default: 800, min: 100, max: 2000, group: 'Grillage', role: 'geometry' },
    { key: 'tfG', label: 'Grillage flange thickness', unit: 'mm', type: 'number', default: 50, min: 8, max: 150, group: 'Grillage', role: 'geometry', design: true },
    { key: 'twG', label: 'Grillage web thickness', unit: 'mm', type: 'number', default: 50, min: 6, max: 100, group: 'Grillage', role: 'geometry' },
    { key: 'sFrame', label: 'Spacing of the barge frames supporting the grillage', unit: 'm', type: 'number', default: 2.4, min: 0.5, max: 12, group: 'Grillage', role: 'bc' },
    { key: 'nSpan', label: 'Number of frame spacings covered by the grillage beam', unit: '-', type: 'int', default: 4, min: 1, max: 12, group: 'Grillage', role: 'bc', design: true, help: 'The cargo leg stands at the middle of the span nearest mid-length; the beam spreads its reaction to nSpan + 1 barge frames (rigid supports).' },
    { key: 'frameCap', label: 'Allowable reaction on one barge frame', unit: 'kN', type: 'number', default: 24000, min: 10, max: 2e5, group: 'Grillage', role: 'bc', help: 'From the barge deck-strength documentation at a transverse or longitudinal bulkhead/frame.' },
    { key: 'tieCap', label: 'Tie-down capacity per cargo support', unit: 'kN', type: 'number', default: 8000, min: 0, max: 2e5, group: 'Seafastening', role: 'bc', design: true, help: 'Design resistance of the uplift restraint (welded plates, clamps or lashings) at one support; zero when none is fitted.' },
    { key: 'fLoad', label: 'Load factor on environmental inertia loads', unit: '-', type: 'number', default: 1.3, min: 1, max: 2, group: 'Criteria', role: 'criteria', help: CONFIRM },
    { key: 'fG', label: 'Load factor on self-weight in the same combination', unit: '-', type: 'number', default: 1.0, min: 0.8, max: 1.5, group: 'Criteria', role: 'criteria', help: CONFIRM },
    { key: 'gM', label: 'Material factor', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.5, group: 'Criteria', role: 'criteria', help: CONFIRM },
    { key: 'days', label: 'Voyage duration', unit: 'days', type: 'number', default: 14, min: 0.2, max: 120, group: 'Voyage', role: 'operation' },
    { key: 'snName', label: 'S–N curve of the seafastening weld', type: 'select', default: 'F', group: 'Voyage', role: 'material', options: SN_CURVE_NAMES.map((n) => ({ value: n, label: `DNV-RP-C203 curve ${n} (air)` })), help: 'Load-carrying fillet welds are typically F to W3 depending on detail and throat.' },
    { key: 'dff', label: 'Design fatigue factor for the voyage', unit: '-', type: 'number', default: 3, min: 1, max: 10, group: 'Criteria', role: 'criteria', help: CONFIRM },
  ],
  compute(inp, ctx) {
    // t → kg, mm → m, deg → rad, g → m/s².
    const steel = ctx.material(inp.mat); need(steel, `Steel '${inp.mat}' is not in the materials library.`);
    const fy = steel.fy, fu = steel.fu || 1.3 * fy, E = steel.E || 2.1e11, fd = fy / inp.gM;
    const m = inp.W * inp.fWeight * 1e3, z = inp.zc + inp.zDeck, ah = inp.heaveG * G;
    const thR = inp.rollDeg * DEG, thP = inp.pitchDeg * DEG, wR = 2 * Math.PI / inp.rollT, wP = 2 * Math.PI / inp.pitchT;
    // Accelerations in the deck frame at the cargo CoG: tangential rotary term + gravity/heave component.
    const aT = thR * wR * wR * z + (G + ah) * Math.sin(thR), aL = thP * wP * wP * z + (G + ah) * Math.sin(thP);
    const aVr = [(G + ah) * Math.cos(thR) + thR * wR * wR * Math.abs(inp.yc), (G - ah) * Math.cos(thR) - thR * wR * wR * Math.abs(inp.yc)];
    const aVp = [(G + ah) * Math.cos(thP) + thP * wP * wP * Math.abs(inp.xc), (G - ah) * Math.cos(thP) - thP * wP * wP * Math.abs(inp.xc)];
    const Fwind = 0.5 * 1.225 * inp.Cs * inp.Awind * inp.Vwind ** 2;
    const FT = m * aT + Fwind, FL = m * aL;
    // Factored design actions (inertia part by the environmental factor, static gravity part by the weight factor).
    const fac = (aTot, aStat) => inp.fG * aStat + inp.fLoad * (aTot - aStat);
    const HT = m * fac(aT, G * Math.sin(thR)) + inp.fLoad * Fwind, HL = m * fac(aL, G * Math.sin(thP));
    const MT = HT * inp.zc, ML = HL * inp.zc; // wind taken to act at the cargo CoG height
    const Vmax = (aV, st) => m * fac(aV[0], st), Vmin = (aV, st) => m * (inp.fG * st - inp.fLoad * (st - aV[1]));
    const sR = G * Math.cos(thR), sP = G * Math.cos(thP);
    // Four supports: two per side; overturning couple taken by the two supports on each side.
    const Rroll = [Vmax(aVr, sR) / 4 + MT / (2 * inp.b), Vmin(aVr, sR) / 4 - MT / (2 * inp.b)];
    const Rpitch = [Vmax(aVp, sP) / 4 + ML / (2 * inp.a), Vmin(aVp, sP) / 4 - ML / (2 * inp.a)];
    const Rmax = Math.max(Rroll[0], Rpitch[0]), Rmin = Math.min(Rroll[1], Rpitch[1]);
    const uplift = Math.max(-Rmin, 0);
    // Sliding force left for the braces after any credited friction on the minimum vertical load.
    const ST = Math.max(HT - inp.mu * Math.max(Vmin(aVr, sR), 0), 0), SL = Math.max(HL - inp.mu * Math.max(Vmin(aVp, sP), 0), 0);
    const cb = Math.cos(inp.braceDeg * DEG);
    const NT = ST / (inp.nT * cb), NL = SL / (inp.nL * cb), Nbr = Math.max(NT, NL);
    const sec = tube(inp.Dbr / 1e3, inp.tbr / 1e3);
    need(inp.tbr < inp.Dbr / 2, 'Brace wall thickness must be smaller than the radius.');
    const lam = (inp.Lbr / sec.r) / (Math.PI * Math.sqrt(E / fy)), chi = chiEC3(lam, 0.49), NbRd = chi * sec.A * fd;
    // Weld: brace force on the perimeter fillet, directional method simplified to shear on the throat.
    const Lw = Math.PI * inp.Dbr / 1e3 / cb, tauW = Nbr / (inp.aw / 1e3 * Lw), fvw = fu / (Math.sqrt(3) * 0.9 * 1.25);
    // Grillage beam: continuous over nSpan frame spacings on rigid frames, the leg reaction at mid-length.
    const gs = run(() => iSection({ h: inp.hG / 1e3, b: inp.bG / 1e3, tw: inp.twG / 1e3, tf: inp.tfG / 1e3 }), 'Grillage section');
    // The leg stands at the middle of the span nearest mid-length (worst position for bending).
    const nSp = inp.nSpan, Lg = nSp * inp.sFrame, xL = (Math.floor((nSp - 1) / 2) + 0.5) * inp.sFrame, xs = [];
    for (let i = 0; i <= nSp; i++) xs.push(i * inp.sFrame);
    xs.push(xL);
    const order = xs.map((x, i) => i).sort((p, q) => xs[p] - xs[q]);
    const nodes = order.map((i) => [xs[i], 0, 0]), isSup = order.map((i) => i <= nSp), kLoad = order.indexOf(nSp + 1);
    const members = []; for (let i = 0; i < nodes.length - 1; i++) members.push({ a: i, b: i + 1, section: 'g', material: 'st' });
    const supports = []; let first = true;
    isSup.forEach((s, i) => { if (s) { supports.push({ node: i, fix: first ? [1, 1, 1, 1, 0, 0] : [0, 1, 1, 0, 0, 0] }); first = false; } });
    const model = { nodes, members, supports, sections: { g: { A: gs.A, Iy: gs.Iy, Iz: gs.Iz, J: gs.J } }, materials: { st: { E, nu: 0.3, rho: 7850 } }, loadCases: [{ name: 'leg', nodal: [{ node: kLoad, F: [0, 0, -Rmax, 0, 0, 0] }] }] };
    const fe = run(() => analyzeFrame(model, { shearDeformation: false, cases: ['leg'] }), 'Grillage frame analysis').results.leg;
    let Mg = 0, Vg = 0; for (const mf of fe.memberForces) { Mg = Math.max(Mg, Math.abs(mf.endA.My), Math.abs(mf.endB.My)); Vg = Math.max(Vg, Math.abs(mf.endA.Vz), Math.abs(mf.endB.Vz)); }
    const Rfr = fe.reactions.map((r) => r.R[2]), RfrMax = Math.max(...Rfr.map(Math.abs));
    const defl = Math.abs(fe.displacements[6 * kLoad + 2]);
    const sigG = Mg / gs.Wely, tauG = Vg / (gs.hw * inp.twG / 1e3);
    // Voyage fatigue of the brace weld: ranges Rayleigh distributed, extreme range = 2 × the unfactored design amplitude.
    const NbrChar = Math.max(m * aT / (inp.nT * cb), m * aL / (inp.nL * cb));
    const dS = 2 * NbrChar / (inp.aw / 1e3 * Lw), nCyc = inp.days * 86400 / inp.rollT;
    const dmg = dS <= 0 ? 0 : run(() => minerWeibull({ sn: snConstants(snCurve(inp.snName, 'air')), shape: 2, cycles: nCyc, maxStressRange: dS }).damage, 'Fatigue damage');
    const checks = [
      check('Seafastening brace: buckling / yield', Nbr / 1e3, NbRd / 1e3, 'kN', `${N001} Sec. 11 (seafastening design loads); EN 1993-1-1 §6.3.1 (curve c)`, 'Axial force in the most loaded brace against its buckling resistance. A buckled brace sheds its load to the others and the cargo can shift.'),
      check('Seafastening weld: shear on the throat', tauW / 1e6, fvw / 1e6, 'MPa', 'EN 1993-1-8 §4.5.3.3 (simplified method, fvw,d = fu/(√3·βw·γM2))', 'Brace end weld to deck or cargo. Weld failure is brittle and gives no warning; it is the usual weak point of seafastening.'),
      check('Grillage beam: bending', sigG / 1e6, fd / 1e6, 'MPa', `${N001} Sec. 11 (grillage); EN 1993-1-1 §6.2.5`, 'Bending of the grillage beam spreading the cargo leg reaction to the barge frames. Yielding lets the cargo settle and overloads the nearest frame.'),
      check('Grillage beam: web shear', tauG / 1e6, fd / Math.sqrt(3) / 1e6, 'MPa', 'EN 1993-1-1 §6.2.6', 'Shear in the grillage web next to the leg. Short stiff grillages are usually governed by shear rather than bending.'),
      check('Barge frame reaction', RfrMax / 1e3, inp.frameCap, 'kN', `${N001} Sec. 11 (vessel strength: local deck loads); barge deck-load documentation`, 'Largest reaction delivered to one barge frame. Exceedance damages the deck structure under the grillage.'),
      check('Voyage fatigue of the brace weld', dmg * inp.dff, 1, 'Miner sum × DFF', 'DNV-RP-C203 (2021) §5.1 (closed-form Weibull damage, shape 2)', 'Every roll cycle of the voyage stresses the weld. A long tow with a poor weld detail can crack the seafastening before arrival.'),
    ];
    if (uplift > 0) checks.push(check('Uplift at a cargo support against tie-down capacity', uplift / 1e3, Math.max(inp.tieCap, 1e-6), 'kN', `${N001} Sec. 11 (uplift)`, 'The overturning moment exceeds the stabilising weight on one side: without a tie-down designed for this force the cargo rocks on its supports.'));
    const gov = worst(checks), ok = gov.uc <= 1;
    const braceMass = (inp.nT + inp.nL) * sec.mass * inp.Lbr / 1e3, grillMass = 4 * gs.A * 7850 * Lg / 1e3, steel_t = braceMass + grillMass;
    const tReq = inp.tbr * Nbr / NbRd;
    const rolls = [5, 10, 15, 20, 25, 30];
    return {
      kpis: [
        kpi('aT', 'Transverse acceleration at the cargo CoG', aT / G, 'g', 'Parallel to the deck in the roll case: rotary acceleration θ·ω²·z plus the gravity and heave component sin θ·(g + a_heave).'),
        kpi('aL', 'Longitudinal acceleration at the cargo CoG', aL / G, 'g', 'Parallel to the deck in the pitch case.'),
        kpi('aV', 'Maximum acceleration normal to the deck', Math.max(aVr[0], aVp[0]) / G, 'g', 'Largest deck-normal acceleration including gravity; the minimum value drives uplift.'),
        kpi('FT', 'Transverse force incl. wind (unfactored)', FT / 1e3, 'kN', 'Characteristic transverse force on the cargo: inertia plus wind.'),
        kpi('FL', 'Longitudinal force (unfactored)', FL / 1e3, 'kN', 'Characteristic longitudinal inertia force on the cargo in the pitch case.'),
        kpi('Rmax', 'Maximum design support reaction', Rmax / 1e3, 'kN', 'Largest factored vertical reaction on one cargo support; the design load of the grillage under that leg.'),
        kpi('Rmin', 'Minimum design support reaction', Rmin / 1e3, 'kN', 'Smallest factored reaction; a negative value is uplift and needs a tie-down of that capacity.'),
        kpi('Nbr', 'Design force in the most loaded brace', Nbr / 1e3, 'kN', 'Factored sliding force shared by the braces of the governing direction, resolved along the brace.', NbRd / 1e3, 'max'),
        kpi('Mg', 'Grillage beam design bending moment', Mg / 1e3, 'kN·m', 'Largest moment in the grillage beam from the beam finite-element model on rigid frames.'),
        kpi('sigG', 'Grillage bending stress', sigG / 1e6, 'MPa', 'Extreme-fibre stress under the design leg reaction.', fd / 1e6, 'max'),
        kpi('defl', 'Grillage deflection under the leg', defl * 1e3, 'mm', 'Elastic deflection at the leg; large values tilt the cargo and change the load sharing between legs.'),
        kpi('Rfr', 'Largest barge frame reaction', RfrMax / 1e3, 'kN', 'Load delivered to the most loaded barge frame under the grillage.', inp.frameCap, 'max'),
        kpi('damage', 'Voyage fatigue damage of the brace weld', dmg, '-', 'Miner sum for the voyage with Rayleigh-distributed stress ranges; multiplied by the design fatigue factor it must stay below 1.', 1 / inp.dff, 'max'),
        kpi('steel', 'Seafastening and grillage steel', steel_t, 't', 'Braces plus four grillage beams of the size entered: temporary steel that is cut off after arrival.'),
      ],
      checks,
      plots: [
        bar('Design support reactions', ['Roll: max', 'Roll: min', 'Pitch: max', 'Pitch: min'], 'Reaction per support', 'kN', [{ name: 'Factored reaction', y: [Rroll[0], Rroll[1], Rpitch[0], Rpitch[1]].map((v) => round(v / 1e3)) }],
          'Largest and smallest reaction on one of four supports in the roll and pitch cases. A negative bar is uplift.'),
        line('Transverse force against roll amplitude', 'Roll single amplitude', 'deg', 'Transverse force on the cargo', 'kN', [{ name: `Period ${inp.rollT} s`, x: rolls, y: rolls.map((r) => round((m * (r * DEG * wR * wR * z + (G + ah) * Math.sin(r * DEG)) + Fwind) / 1e3)) }, { name: 'Period × 1.3', x: rolls, y: rolls.map((r) => round((m * (r * DEG * (wR / 1.3) ** 2 * z + (G + ah) * Math.sin(r * DEG)) + Fwind) / 1e3)) }],
          'How the seafastening load grows with the roll criterion; the second curve shows the relief from a longer roll period. Use it to judge what a voyage-specific motion analysis could save.'),
        line('Grillage beam: frame reactions', 'Position along the beam', 'm', 'Frame reaction', 'kN', [{ name: 'Reaction on barge frame', x: fe.reactions.map((r) => round(nodes[r.node][0])), y: Rfr.map((v) => round(v / 1e3)) }],
          'How the leg load spreads to the barge frames. A long flexible beam loads only the frames next to the leg; negative values are hold-down forces at the ends.', 'scatter'),
      ],
      tables: [tbl('Load cases at the cargo centre of gravity', [['c', 'Case'], ['ah', 'Along-deck acceleration', 'g'], ['av1', 'Normal acceleration max', 'g'], ['av2', 'Normal acceleration min', 'g'], ['H', 'Design horizontal force', 'kN'], ['M', 'Design overturning moment', 'kN·m'], ['N', 'Force per brace', 'kN']], [
        { c: 'Roll ± heave + wind', ah: round(aT / G), av1: round(aVr[0] / G), av2: round(aVr[1] / G), H: round(HT / 1e3), M: round(MT / 1e3), N: round(NT / 1e3) },
        { c: 'Pitch ± heave', ah: round(aL / G), av1: round(aVp[0] / G), av2: round(aVp[1] / G), H: round(HL / 1e3), M: round(ML / 1e3), N: round(NL / 1e3) },
      ], 'The two quasi-static voyage load cases with the accelerations they produce and the factored actions on seafastening and supports.')],
      notes: [`Roll ${inp.rollDeg}° / ${inp.rollT} s and pitch ${inp.pitchDeg}° / ${inp.pitchT} s with ±${inp.heaveG} g heave; lever to the motion centre ${fmt(z)} m.`, uplift > 0 ? `Uplift of ${fmt(uplift / 1e3)} kN occurs at a support: tie-downs are required.` : 'No uplift: all supports stay in compression.', `"${gov.name}" governs at ${fmt(gov.uc)}.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : uplift > 0 ? 'warn' : 'ok', `Seafastening and grillage pass; "${gov.name}" governs at ${fmt(gov.uc)}. Transverse acceleration ${fmt(aT / G)} g, brace force ${fmt(Nbr / 1e3)} kN.${uplift > 0 ? ` Uplift of ${fmt(uplift / 1e3)} kN must be tied down.` : ''}`,
          `Replace the default motion criteria with the voyage motion analysis for the actual barge and route, and design the brace-end welds and the deck under-structure for ${fmt(Nbr / 1e3)} kN per brace.`, `${N001} Sec. 11`)
          : adv('engineering', 'fail', `"${gov.name}" is exceeded (${fmt(gov.uc)}).`,
            gov.name.startsWith('Seafastening brace') ? `Increase the brace wall to about ${fmt(Math.ceil(tReq))} mm or add braces: ${fmt(Math.ceil(inp.nT * gov.uc))} transverse braces of this size would pass.` : gov.name.startsWith('Seafastening weld') ? `Increase the weld throat to ${fmt(Math.ceil(inp.aw * gov.uc))} mm or use a full-penetration weld.` : gov.name.startsWith('Grillage') ? `Deepen or thicken the grillage beam (section modulus × ${fmt(gov.uc)}), or shorten its span onto more frames.` : gov.name.startsWith('Uplift') ? `Provide tie-downs for ${fmt(uplift / 1e3)} kN per support, or widen the support spacing to about ${fmt(inp.b * Math.sqrt(gov.uc))} m.` : gov.name.startsWith('Barge') ? 'Lengthen the grillage to engage more frames, or reposition the cargo legs over bulkheads.' : 'Improve the weld detail category, increase the throat, or shorten the exposure with a weather-routed voyage.',
            `${N001} Sec. 11; EN 1993-1-1; DNV-RP-C203`),
        adv('economic', 'info', `Temporary steel for this arrangement is about ${fmt(steel_t)} t; the roll criterion alone accounts for ${fmt(100 * (m * (thR * wR * wR * z + (G + ah) * Math.sin(thR))) / FT)} % of the transverse force.`,
          'A voyage-specific motion analysis or weather-routed tow typically cuts the default roll criterion by a third, and the seafastening steel and offshore cutting time with it; weigh that study against the fabrication, installation and removal of the steel saved.',
          'Seafastening is priced per tonne fabricated, fitted and later removed offshore; AACE RP 18R-97 estimate classes',
          `Seafastening and grillage are scrapped or recycled after one voyage: ${steelCo2Text(steel_t)} Reusable grillages and bolted or clamped seafastening avoid most of it.`),
      ],
      publish: { 'install.seafasteningForce': Math.max(ST, SL), 'install.transportAccT': aT, 'install.transportAccL': aL },
      quantities: { steel_t: round(steel_t) },
    };
  },
  verify: [
    { name: 'Transverse acceleration: a = θ·(2π/T)²·z + (g + a_heave)·sin θ', source: 'DNV-ST-N001 simplified motion criteria, hand calculation: θ = 20° = 0.34907 rad, T = 10 s, z = 17 m, heave 0.2 g → 0.34907·0.39478·17 + 1.2·9.80665·0.34202 = 2.3427 + 4.0249 = 6.3676 m/s² = 0.6493 g',
      inputs: {}, expect: [{ kpi: 'aT', value: (20 * Math.PI / 180 * (2 * Math.PI / 10) ** 2 * 17 + 1.2 * 9.80665 * Math.sin(20 * Math.PI / 180)) / 9.80665, tol: 1e-9 }] },
    { name: 'Grillage beam on two frames: M = P·L/4, δ = P·L³/(48·E·I)', source: 'Simply supported beam with a central point load (Euler–Bernoulli); one frame spacing of 2.4 m; checked through the ratio M/Rmax = L/4 = 0.6 m with zero motion (Rmax = W·g/4)',
      inputs: { nSpan: 1, rollDeg: 0, pitchDeg: 0, heaveG: 0, Vwind: 0, fG: 1, fLoad: 1 }, expect: [{ kpi: 'Rmax', value: 6.3e6 * 9.80665 / 4 / 1e3, tol: 1e-9 }, { kpi: 'Mg', value: 6.3e6 * 9.80665 / 4 * 2.4 / 4 / 1e3, tol: 1e-6 }] },
    { name: 'Two-span continuous grillage, load at the middle of one span: M = 13·P·L/64, centre reaction 11·P/16', source: 'Three-moment equation for a two-span continuous beam with a point load at the middle of one span (Roark, Formulas for Stress and Strain, continuous beams): support moment 3PL/32, moment under the load 13PL/64, reactions 13P/32, 11P/16, −3P/32',
      inputs: { nSpan: 2, rollDeg: 0, pitchDeg: 0, heaveG: 0, Vwind: 0, fG: 1, fLoad: 1 }, expect: [{ kpi: 'Mg', value: 13 * (6.3e6 * 9.80665 / 4) * 2.4 / 64 / 1e3, tol: 1e-6 }, { kpi: 'Rfr', value: 11 / 16 * 6.3e6 * 9.80665 / 4 / 1e3, tol: 1e-6 }] },
  ],
};
