// Floating platform motion and stability: hydrostatics and intact stability of a three-column semi-submersible
// wind platform under rotor thrust, natural periods, six-degree-of-freedom response amplitude operators and the
// motion statistics (platform motions, nacelle acceleration) in the operating and extreme sea states.
import { multiColumnHydrostatics } from '../../engine/hydro/hydrostatics.js';
import { buildFloater, naturalPeriods, waveSpectrum, floaterRAO, pointMotion, responseStatistics } from '../../engine/hydro/floater.js';
import { site, adv, check, worst, sev, fmt, round, run, need, linspace, G, DEG, LIMIT_CODE, CONFIRM, steelCo2Text } from './_common.js';

const DOF = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'];

export default {
  id: 's16-floater-stability', suite: 16, title: 'Floating platform motions and stability',
  summary: 'Hydrostatics, metacentric height, heel under rotor thrust, intact-stability area ratio, natural periods, RAOs and sea-state motion statistics of a semi-submersible floating wind platform.',
  covers: [7], fidelity: 'reduced-order',
  standards: ['DNV-ST-0119 (2021) Floating wind turbine structures §7 (floating stability) and §4 (motions)', 'IEC TS 61400-3-2:2019 §11 (floating stability)', 'DNV-RP-C205 (2021) §3.5, §6 and §7 (spectra, Morison, floater response)', 'IMO MODU Code 2009 §3.3 (intact stability of column-stabilised units, area ratio 1.3)', 'DNV-OS-C301 (stability and watertight integrity)'],
  models: ['Hydrostatic restoring of a multi-column body (parallel-axis waterplane inertia)', 'Offshore structural stability: metacentric height, righting lever, wind-heeling area ratio', 'Six-degree-of-freedom rigid floating-body dynamics', 'Strip (Morison) added mass and long-wave Froude–Krylov excitation', 'Stochastically linearised viscous drag (Borgman)', 'JONSWAP spectrum and Rayleigh most-probable-maximum statistics', 'Uncoupled and coupled natural periods'],
  limits: `Slender-member long-wave hydrodynamics (valid for column diameter below a fifth of the wavelength): no diffraction solution, frequency-independent strip added mass, radiation damping from the Haskind relation with the same long-wave excitation plus a user fraction of critical; a sea state whose peak period is below 1.4 times the period at which the largest member reaches diameter / wavelength = 0.2 carries a warning (on a single column the surge force was 6 % high at that limit and 47 % high at 0.32 against a panel solution, and surge and pitch differ by up to 25 % within 0.9–1.4 times the pitch natural period); no second-order (slow-drift) response and no aerodynamic or controller coupling (see the coupled dynamics module). The righting lever uses the wall-sided formula, valid until a column emerges or the deck immerses; damaged stability is not covered. ${LIMIT_CODE}`,
  publishes: ['vessel.disp', 'vessel.GM', 'vessel.heaveSig', 'vessel.rollSigDeg', 'vessel.Tn3', 'vessel.Tn4', 'quantities.steel_t', 'vessel.heelDeg'],
  searchTerms: 'semi-submersible floating wind platform stability RAO natural period DNV-ST-0119',
  inputs: [
    { key: 'colR', label: 'Radius of the offset columns from the centre', unit: 'm', type: 'number', default: 28.87, min: 10, max: 80, group: 'Hull geometry', role: 'geometry', design: true, help: 'Three columns at 120°. Column spacing is the main lever on waterplane inertia and stability.' },
    { key: 'colD', label: 'Offset-column diameter', unit: 'm', type: 'number', default: 12, min: 3, max: 30, group: 'Hull geometry', role: 'geometry', design: true },
    { key: 'baseD', label: 'Base-column (heave-plate) diameter', unit: 'm', type: 'number', default: 24, min: 3, max: 50, group: 'Hull geometry', role: 'geometry', design: true },
    { key: 'baseH', label: 'Base-column height', unit: 'm', type: 'number', default: 6, min: 0.5, max: 20, group: 'Hull geometry', role: 'geometry' },
    { key: 'cenD', label: 'Centre-column diameter', unit: 'm', type: 'number', default: 6.5, min: 0, max: 20, group: 'Hull geometry', role: 'geometry', help: '0 = no centre column (turbine on one of the offset columns).' },
    { key: 'draft', label: 'Operating draft', unit: 'm', type: 'number', default: 20, min: 5, max: 60, group: 'Hull geometry', role: 'ic', design: true, help: 'Initial floating position; the platform mass equals the displacement at this draft (ballast adjusts it).' },
    { key: 'freeboard', label: 'Freeboard to the deck', unit: 'm', type: 'number', default: 12, min: 3, max: 30, group: 'Hull geometry', role: 'geometry' },
    { key: 'KG', label: 'Centre of gravity above the keel (platform, tower, turbine, ballast)', unit: 'm', type: 'number', default: 10.11, min: 1, max: 60, group: 'Mass distribution', role: 'ic', uncertain: { dist: 'normal', cov: 0.05 }, design: true },
    { key: 'rxx', label: 'Roll and pitch radius of gyration', unit: 'm', type: 'number', default: 28, min: 5, max: 80, group: 'Mass distribution', role: 'geometry', help: 'About the centre of gravity, turbine and tower included.' },
    { key: 'rzz', label: 'Yaw radius of gyration', unit: 'm', type: 'number', default: 30, min: 5, max: 80, group: 'Mass distribution', role: 'geometry' },
    { key: 'steelFrac', label: 'Hull steel mass / displacement', unit: '-', type: 'number', default: 0.28, min: 0.1, max: 0.6, group: 'Mass distribution', role: 'geometry', help: 'Primary and secondary hull steel as a share of displacement (the rest is ballast, turbine and mooring load): 0.2–0.35 for steel semi-submersibles.' },
    { key: 'thrust', label: 'Rotor thrust', unit: 'kN', type: 'number', default: 750, min: 0, max: 10000, group: 'Turbine loads', role: 'bc', bind: 'owt.thrustRated', uncertain: { dist: 'normal', cov: 0.1 } },
    { key: 'hubHeight', label: 'Hub height above still water', unit: 'm', type: 'number', default: 90, min: 30, max: 220, group: 'Turbine loads', role: 'geometry', bind: 'owt.hubHeight' },
    { key: 'fairleadDepth', label: 'Mooring fairlead depth (centre of lateral resistance)', unit: 'm', type: 'number', default: 14, min: 0, max: 60, group: 'Turbine loads', role: 'bc', help: 'The thrust is reacted at the fairleads: the heeling lever is hub height plus this depth.' },
    { key: 'Kxx', label: 'Mooring surge stiffness', unit: 'kN/m', type: 'number', default: 120, min: 0, max: 50000, group: 'Turbine loads', role: 'bc', bind: 'mooring.Kxx', help: 'Horizontal restraint from the mooring module (yaw stiffness is taken as Kxx × fairlead radius²).' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 200, min: 30, max: 2000, group: 'Sea states', role: 'environment', bind: 'site.depth' },
    { key: 'Hs', label: 'Operating significant wave height', unit: 'm', type: 'number', default: 4, min: 0.5, max: 12, group: 'Sea states', role: 'bc', bind: 'site.Hs1', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Sea state in which the turbine must keep producing.' },
    { key: 'Tp', label: 'Operating peak period', unit: 's', type: 'number', default: 9, min: 3, max: 20, group: 'Sea states', role: 'bc', bind: 'site.Tp1' },
    site('Hs100', { group: 'Sea states' }), site('Tp100', { group: 'Sea states' }),
    { key: 'heading', label: 'Wave heading', unit: 'deg', type: 'number', default: 0, min: 0, max: 180, group: 'Sea states', role: 'bc' },
    { key: 'zetaRad', label: 'Added damping ratio per mode', unit: '-', type: 'number', default: 0.03, min: 0.005, max: 0.3, group: 'Sea states', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.4 }, help: 'Fraction of critical damping added in every restored mode on top of the computed wave-radiation damping (Haskind relation) and the linearised drag: it stands for mooring, cable and appendage damping. Calibrate against decay tests or a diffraction analysis.' },
    { key: 'heelLimit', label: 'Allowable static heel in operation', unit: 'deg', type: 'number', default: 6, min: 1, max: 15, group: 'Acceptance', role: 'criteria', help: `Turbine suppliers typically allow 5–7° mean and 10° maximum tilt. ${CONFIRM}` },
    { key: 'gmMin', label: 'Minimum metacentric height', unit: 'm', type: 'number', default: 1, min: 0.1, max: 10, group: 'Acceptance', role: 'criteria', help: `DNV-ST-0119 / DNV-OS-C301: GM ≥ 1.0 m for column-stabilised units in operation. ${CONFIRM}` },
    { key: 'areaRatio', label: 'Required righting / heeling area ratio', unit: '-', type: 'number', default: 1.3, min: 1, max: 2, group: 'Acceptance', role: 'criteria', help: '1.3 for column-stabilised units (IMO MODU Code / DNV-ST-0119).' },
    { key: 'downflood', label: 'Downflooding or limiting heel angle', unit: 'deg', type: 'number', default: 25, min: 8, max: 50, group: 'Acceptance', role: 'criteria' },
    { key: 'periodMin', label: 'Minimum heave and pitch natural period', unit: 's', type: 'number', default: 16, min: 5, max: 40, group: 'Acceptance', role: 'criteria', help: 'Above the peak periods of the energetic sea states (about 1.25 × the extreme peak period); semi-submersibles reach 17–25 s in heave.' },
    { key: 'accLimit', label: 'Allowable nacelle acceleration (3-h maximum)', unit: 'm/s²', type: 'number', default: 3, min: 0.5, max: 8, group: 'Acceptance', role: 'criteria', help: 'About 0.2–0.4 g by turbine supplier.' },
    { key: 'nOmega', label: 'Frequency points', type: 'int', default: 80, min: 20, max: 640, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nOmega', label: 'Frequency points', levels: [40, 80, 160, 320], qoi: 'heaveSig' },
  compute(inp) {
    const T = inp.draft, R = inp.colR;
    need(inp.baseH < T - 1, 'The base columns must be shorter than the draft.');
    need(inp.baseD >= inp.colD, 'The base-column diameter cannot be smaller than the upper column.');
    need(R > inp.baseD / 2 + inp.cenD / 2, 'The base columns overlap the centre column: increase the column radius.');
    const pos = [0, 120, 240].map((a) => [R * Math.cos(a * DEG), R * Math.sin(a * DEG)]);
    const vBase = Math.PI * inp.baseD ** 2 / 4 * inp.baseH;
    const hs = run(() => multiColumnHydrostatics({ columns: [...pos.map(([x, y]) => ({ x, y, D: inp.colD, zBottom: inp.baseH })), ...(inp.cenD > 0 ? [{ x: 0, y: 0, D: inp.cenD, zBottom: 0 }] : [])],
      submerged: pos.map(([x, y]) => ({ volume: vBase, x, y, z: inp.baseH / 2 })), T, KG: inp.KG }), 'Hydrostatics');
    const GM = Math.min(hs.GMt, hs.GMl), disp = hs.displacement, zG = inp.KG - T;
    need(GM > 0, `The platform is unstable upright: GM = ${fmt(GM)} m. Lower the centre of gravity (more ballast low down) or widen the column spacing.`);
    // static heel under thrust: moment about the mooring fairlead level against the hydrostatic pitch stiffness
    const Mheel = inp.thrust * 1e3 * (inp.hubHeight + inp.fairleadDepth), C55 = disp * G * GM;
    const heel = Mheel / C55;
    // wall-sided righting lever and area ratio to the limiting angle (heeling lever decays with cos²)
    const BM = Math.min(hs.BMt, hs.BMl), lever0 = Mheel / (disp * G);
    const phis = linspace(0, inp.downflood, 51), gz = phis.map((p) => Math.sin(p * DEG) * (GM + 0.5 * BM * Math.tan(p * DEG) ** 2)), hl = phis.map((p) => lever0 * Math.cos(p * DEG) ** 2);
    let aR = 0, aH = 0, second = inp.downflood;
    for (let i = 1; i < phis.length; i++) {
      const dphi = (phis[i] - phis[i - 1]) * DEG;
      aR += 0.5 * (gz[i] + gz[i - 1]) * dphi; aH += 0.5 * (hl[i] + hl[i - 1]) * dphi;
    }
    const ratio = aH > 0 ? aR / aH : 99;
    // slender-member dynamic model (origin on the still-water plane)
    const members = [
      ...pos.flatMap(([x, y]) => [{ type: 'cylinder', p1: [x, y, -T + inp.baseH], p2: [x, y, inp.freeboard], D: inp.colD, Ca: 1, Cd: 0.7 }, { type: 'cylinder', p1: [x, y, -T], p2: [x, y, -T + inp.baseH], D: inp.baseD, Ca: 1, Cd: 0.8, end1: 'half-disc', end2: 'half-disc', CdAxial: 4 }]),
      ...(inp.cenD > 0 ? [{ type: 'cylinder', p1: [0, 0, -T], p2: [0, 0, inp.freeboard], D: inp.cenD, Ca: 1, Cd: 0.7, end1: 'half-disc' }] : []),
    ];
    const Kext = Array.from({ length: 6 }, () => new Array(6).fill(0));
    Kext[0][0] = Kext[1][1] = inp.Kxx * 1e3; Kext[5][5] = inp.Kxx * 1e3 * (R + inp.colD / 2) ** 2;
    const model = run(() => buildFloater({ members, mass: disp, cog: [0, 0, zG], gyradii: [inp.rxx, inp.rxx, inp.rzz], depth: inp.depth, Kext, nStrips: 10 }), 'Floater model');
    const Tn = naturalPeriods(model).uncoupled;
    const omega = linspace(0.12, 2.2, inp.nOmega);
    const resp = (Hs, Tp, dur) => {
      const sea = { Hs, Tp, gamma: 3.3 };
      const rao = run(() => floaterRAO({ model, omega, headingDeg: inp.heading, seaState: sea, extraDampingRatio: inp.zetaRad }), 'Response amplitude operators');
      const S = waveSpectrum(omega, sea);
      const st = rao.amp.map((amp) => responseStatistics({ omega, amp, S, duration: dur }));
      const nac = pointMotion({ rao, point: [0, 0, inp.hubHeight] });
      const acc = responseStatistics({ omega, amp: nac.accAmp[0], S, duration: dur });
      return { rao, S, st, acc };
    };
    const op = resp(inp.Hs, inp.Tp, 10800), ex = resp(inp.Hs100, inp.Tp100, 10800);
    const pitchSig = op.st[4].significant / DEG, pitchMax = ex.st[4].mpm / DEG, heaveSig = op.st[2].significant;
    const dynTilt = heel / DEG;   // operating mean tilt; storm pitch is reported separately (turbine parked)
    const checks = [
      check('Metacentric height', inp.gmMin, GM, 'm', 'DNV-ST-0119 (2021) §7; DNV-OS-C301 (GM ≥ 1.0 m, column-stabilised units)', 'Demand is the required GM, capacity the calculated one. Low GM means small restoring moment: large heel under thrust and slow, large-amplitude pitch.'),
      check('Static heel under rotor thrust', heel / DEG, inp.heelLimit, 'deg', 'DNV-ST-0119 (2021) §4 / turbine supplier tilt limit', 'Mean tilt of the platform at rated thrust against the allowable operating tilt. Larger tilt loses energy (cos³), raises tower-base bending and wears bearings; active ballast is the usual remedy.'),
      check('Intact stability: righting / heeling area ratio', inp.areaRatio, ratio, '-', 'IMO MODU Code 2009 §3.3.2 (ratio ≥ 1.3); DNV-ST-0119 §7.2', 'Required against actual ratio of the area under the righting-lever curve to that under the heeling-lever curve up to the limiting angle. Below the requirement a gust could heel the platform to downflooding.'),
      check('Nacelle fore-aft acceleration, operating sea state (3-h maximum)', op.acc.mpm, inp.accLimit, 'm/s²', 'Turbine supplier limit; DNV-ST-0119 §4 (motion criteria)', 'Most probable maximum acceleration at hub height from wave-induced platform motion against the allowable for nacelle components.'),
      check('Heave and pitch natural periods above the wave range', inp.periodMin, Math.min(Tn[2] ?? 99, Tn[4] ?? 99), 's', 'DNV-RP-C205 (2021) §7.1 (natural periods of column-stabilised units)', 'Required minimum against the shorter of the heave and pitch natural periods. Periods inside the range of energetic waves mean first-order resonance in common sea states.'),
    ];
    const gov = worst(checks);
    const steel = inp.steelFrac * disp / 1e3;
    // KG that would bring the heel to the limit: GM_req = M / (Δ g φ_lim)
    const gmReq = Mheel / (disp * G * inp.heelLimit * DEG), kgReq = hs.KB + BM - gmReq;
    const idx = linspace(0, omega.length - 1, Math.min(omega.length, 120)).map(Math.round);
    const per = idx.map((i) => round(2 * Math.PI / omega[i], 3)).reverse();
    const ser = (arr, f = 1) => idx.map((i) => round(arr[i] * f, 5)).reverse();
    return {
      kpis: [
        { key: 'disp', label: 'Displacement', value: round(disp / 1e3, 5), unit: 't', explain: 'Mass of water displaced at the operating draft; it equals the total mass of hull, ballast, tower and turbine plus the vertical mooring load.' },
        { key: 'GM', label: 'Metacentric height GM', value: round(GM, 4), unit: 'm', limit: inp.gmMin, sense: 'min', explain: 'Height of the metacentre above the centre of gravity: KB + BM − KG. The initial stability of a semi-submersible comes almost entirely from the column spacing (waterplane inertia).' },
        { key: 'KB', label: 'Centre of buoyancy above keel', value: round(hs.KB, 4), unit: 'm', explain: 'Vertical centroid of the displaced volume.' },
        { key: 'BM', label: 'Metacentric radius BM', value: round(BM, 4), unit: 'm', explain: 'Waterplane second moment over displaced volume: the stabilising contribution of the columns.' },
        { key: 'heel', label: 'Static heel under rotor thrust', value: round(heel / DEG, 4), unit: 'deg', limit: inp.heelLimit, sense: 'max', explain: 'Mean inclination at rated thrust from the hydrostatic pitch stiffness ρg∇·GM. It adds to wave-induced pitch.' },
        { key: 'ratio', label: 'Righting / heeling area ratio', value: round(Math.min(ratio, 99), 4), unit: '-', limit: inp.areaRatio, sense: 'min', explain: `Energy margin of the righting moment over the thrust heeling moment up to ${inp.downflood}°.` },
        { key: 'TnHeave', label: 'Heave natural period', value: round(Tn[2] ?? 0, 4), unit: 's', explain: 'Vertical oscillation period 2π√((M + A₃₃)/ρgA_w). Semi-submersibles aim for 17–25 s, above the bulk of wave energy.' },
        { key: 'TnPitch', label: 'Pitch natural period', value: round(Tn[4] ?? 0, 4), unit: 's', explain: 'Pitch period from inertia (with added inertia) over hydrostatic stiffness; 20–35 s is typical. The pitch controller of the turbine must be slower than this.' },
        { key: 'TnSurge', label: 'Surge natural period', value: round(Tn[0] ?? 0, 4), unit: 's', explain: 'Horizontal period on the mooring stiffness (0 if no mooring stiffness is entered): far above wave periods, excited by slow-drift forces.' },
        { key: 'heaveSig', label: 'Significant heave amplitude, operating sea', value: round(heaveSig, 4), unit: 'm', explain: 'Twice the standard deviation of heave in the operating sea state.' },
        { key: 'pitchSig', label: 'Significant pitch amplitude, operating sea', value: round(pitchSig, 4), unit: 'deg', explain: 'Twice the standard deviation of wave-induced pitch in the operating sea state (added to the mean thrust heel).' },
        { key: 'surgeSig', label: 'Significant surge amplitude, operating sea', value: round(op.st[0].significant, 4), unit: 'm', explain: 'Wave-frequency surge motion; slow-drift motion on the mooring comes on top and is treated in the mooring module.' },
        { key: 'accNac', label: 'Nacelle acceleration, 3-h maximum, operating sea', value: round(op.acc.mpm, 4), unit: 'm/s²', limit: inp.accLimit, sense: 'max', explain: 'Most probable maximum horizontal acceleration at hub height from surge and pitch.' },
        { key: 'pitchMaxStorm', label: 'Pitch, 3-h maximum, extreme sea', value: round(pitchMax, 4), unit: 'deg', explain: 'Most probable maximum wave-induced pitch in the extreme sea state (turbine parked): a survival quantity for tower and nacelle.' },
        { key: 'heaveMaxStorm', label: 'Heave, 3-h maximum, extreme sea', value: round(ex.st[2].mpm, 4), unit: 'm', explain: 'Most probable maximum heave in the extreme sea state: compare with the air gap under the deck and the cable hang-off design.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Response amplitude operators', x: { label: 'Wave period', unit: 's' }, y: { label: 'm/m · deg/m', unit: '-' }, series: [{ name: 'Surge (m/m)', x: per, y: ser(op.rao.amp[0]) }, { name: 'Heave (m/m)', x: per, y: ser(op.rao.amp[2]) }, { name: 'Pitch (deg/m)', x: per, y: ser(op.rao.amp[4], 1 / DEG) }],
          explain: 'Motion per metre of wave amplitude against wave period for the operating sea state. Peaks mark the natural periods; they should lie to the right of the periods where the sea has energy.' },
        { type: 'line', title: 'Wave spectra and natural periods', x: { label: 'Wave period', unit: 's' }, y: { label: 'Spectral density', unit: 'm²·s/rad' }, series: [{ name: 'Operating sea', x: per, y: ser(op.S) }, { name: 'Extreme sea', x: per, y: ser(ex.S) }, { name: 'Heave natural period', x: [Tn[2] ?? 0, Tn[2] ?? 0], y: [0, Math.max(...ex.S)], dash: true }, { name: 'Pitch natural period', x: [Tn[4] ?? 0, Tn[4] ?? 0], y: [0, Math.max(...ex.S)], dash: true }],
          explain: 'Energy of the two sea states against period with the heave and pitch natural periods. Overlap of a dashed line with the extreme spectrum means resonant storm response.' },
        { type: 'line', title: 'Righting and heeling levers', x: { label: 'Heel angle', unit: 'deg' }, y: { label: 'Lever', unit: 'm' }, series: [{ name: 'Righting lever GZ (wall-sided)', x: phis.map((v) => round(v, 3)), y: gz.map((v) => round(v, 4)) }, { name: 'Thrust heeling lever', x: phis.map((v) => round(v, 3)), y: hl.map((v) => round(v, 4)) }],
          explain: 'Righting lever of the hull and heeling lever of the rotor thrust against heel angle. Their first crossing is the static heel; the area between the curves beyond it is the reserve against gusts.' },
        { type: 'model3d', title: 'Hull columns', nodes: [...pos.flatMap(([x, y]) => [[round(x, 2), round(y, 2), -T], [round(x, 2), round(y, 2), -T + inp.baseH], [round(x, 2), round(y, 2), inp.freeboard]]), [0, 0, -T], [0, 0, inp.freeboard], [0, 0, inp.hubHeight]], lines: [[0, 1], [1, 2], [3, 4], [4, 5], [6, 7], [7, 8], [9, 10], [10, 11], [1, 9], [4, 9], [7, 9], [2, 10], [5, 10], [8, 10]], tris: [],
          nodeValues: [inp.baseD, inp.baseD, inp.colD, inp.baseD, inp.baseD, inp.colD, inp.baseD, inp.baseD, inp.colD, inp.cenD, inp.cenD, 0], valueLabel: 'Member diameter', valueUnit: 'm',
          explain: 'Centre lines of the three offset columns with their base columns, the centre column and the tower to the hub, with the bracing shown schematically. Values are member diameters.' },
      ],
      tables: [
        { title: 'Natural periods and motion statistics', columns: [{ key: 'dof', label: 'Motion' }, { key: 'Tn', label: 'Natural period', unit: 's' }, { key: 'sigOp', label: 'Significant amplitude, operating' }, { key: 'maxOp', label: '3-h maximum, operating' }, { key: 'maxEx', label: '3-h maximum, extreme' }, { key: 'u', label: 'Unit' }],
          rows: DOF.map((d, i) => { const f = i >= 3 ? 1 / DEG : 1; return { dof: d, Tn: Tn[i] ? round(Tn[i], 4) : '–', sigOp: round(op.st[i].significant * f, 4), maxOp: round(op.st[i].mpm * f, 4), maxEx: round(ex.st[i].mpm * f, 4), u: i >= 3 ? 'deg' : 'm' }; }),
          explain: 'All six rigid-body motions: natural period, significant amplitude in the operating sea and most probable 3-hour maxima in both sea states. Zero rows are motions not excited at this wave heading.' },
        { title: 'Hydrostatics', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [{ q: 'Displaced volume', v: round(hs.volume, 5), u: 'm³' }, { q: 'Waterplane area', v: round(hs.waterplaneArea, 5), u: 'm²' }, { q: 'Waterplane second moment', v: round(Math.min(hs.Ixx, hs.Iyy), 6), u: 'm⁴' }, { q: 'KB / BM / KG', v: `${fmt(hs.KB)} / ${fmt(BM)} / ${fmt(inp.KG)}`, u: 'm' }, { q: 'Heave stiffness', v: round(hs.C33 / 1e3, 5), u: 'kN/m' }, { q: 'Pitch stiffness ρg∇·GM', v: round(C55 / 1e6, 5), u: 'MN·m/rad' }, { q: 'Thrust heeling moment', v: round(Mheel / 1e6, 5), u: 'MN·m' }, { q: 'Valid frequency limit of the slender-member model', v: round(model.validOmegaMax, 4), u: 'rad/s' }],
          explain: 'Upright hydrostatic particulars behind the stability results and the stiffness terms used in the motion model.' },
      ],
      notes: [
        `Semi-submersible: three ${fmt(inp.colD)} m columns at ${fmt(R)} m radius on ${fmt(inp.baseD)} m base columns, draft ${fmt(T)} m, displacement ${fmt(disp / 1e3, 5)} t, GM ${fmt(GM)} m.`,
        `Thrust ${fmt(inp.thrust)} kN at ${inp.hubHeight} m (lever ${fmt(inp.hubHeight + inp.fairleadDepth)} m) heels the platform ${fmt(heel / DEG)}°; with the significant wave pitch the operating tilt reaches about ${fmt(heel / DEG + pitchSig)}°.`,
        `Natural periods: heave ${fmt(Tn[2] ?? 0)} s, pitch ${fmt(Tn[4] ?? 0)} s, surge ${fmt(Tn[0] ?? 0)} s; wave peak periods ${inp.Tp} s (operating) and ${inp.Tp100} s (extreme).`,
        `Governing check: ${gov.name} (utilisation ${fmt(gov.uc)}). The slender-member model is valid up to ${fmt(model.validOmegaMax)} rad/s (periods above ${fmt(2 * Math.PI / model.validOmegaMax)} s).`,
      ],
      advice: [
        ...(Math.min(inp.Tp, inp.Tp100) * model.validOmegaMax < 1.4 * 2 * Math.PI ? [adv('engineering', 'warn',
          `Peak period ${fmt(Math.min(inp.Tp, inp.Tp100))} s is too close to the limit of the slender-member model: its largest member reaches diameter / wavelength = 0.2 at ${fmt(2 * Math.PI / model.validOmegaMax)} s, and a sea with this peak has much of its energy below that period (the limit asks for Tp ≥ ${fmt(1.4 * 2 * Math.PI / model.validOmegaMax)} s). Diffraction, which the model omits, lowers the true loads there: on a single column the surge force was 6 % high at the limit and 47 % high at diameter / wavelength 0.32.`,
          `Take wave-frequency loads and motions for this sea state from a radiation–diffraction analysis, or treat the present values as upper-bound estimates for periods below ${fmt(2 * Math.PI / model.validOmegaMax)} s.`,
          'DNV-RP-C205 (2019) §6.1.2 (slender-member limit); comparison of the slender-member model with Capytaine 3.0.0 on a single column (cross-solver verification of this application).')] : []),
        adv('engineering', sev(gov.uc), `GM = ${fmt(GM)} m, static heel ${fmt(dynTilt)}° at rated thrust (limit ${inp.heelLimit}°), area ratio ${fmt(Math.min(ratio, 99))}; heave and pitch periods ${fmt(Tn[2] ?? 0)} s and ${fmt(Tn[4] ?? 0)} s; nacelle acceleration ${fmt(op.acc.mpm)} m/s².`,
          gov.uc > 1 ? (gov === checks[1] || gov === checks[0] || gov === checks[2] ? `Raise the stability: GM of ${fmt(gmReq)} m is needed for the heel limit — lower KG to about ${fmt(kgReq)} m, widen the column radius by about ${fmt((Math.sqrt(Math.max(gmReq - hs.KB + inp.KG, 0.1) / BM) - 1) * 100)} %, or use active ballast to cancel the mean thrust moment.` : gov === checks[4] ? 'Lengthen the natural periods: larger heave plates (added mass) or a smaller waterplane for heave; more pitch inertia or lower GM for pitch.' : 'Reduce nacelle acceleration: lengthen the pitch period or add heave-plate damping; check with the coupled model including aerodynamic damping.') : 'Stability and motions meet the criteria: confirm the hydrodynamic coefficients with a diffraction–radiation analysis and model test, and check damaged stability (one compartment flooded).',
          'DNV-ST-0119 (2021) §7 (floating stability) and IEC TS 61400-3-2:2019; IMO MODU Code area-ratio criterion for column-stabilised units.',
          `Hull steel at ${fmt(inp.steelFrac * 100)} % of displacement: ${fmt(steel)} t. ${steelCo2Text(steel)} Stability bought with column spacing rather than ballast mass keeps both steel and displacement down.`),
        adv('economic', 'info', `The hull carries about ${fmt(steel)} t of steel for a displacement of ${fmt(disp / 1e3, 5)} t; static heel costs energy as roughly 1 − cos³(tilt) = ${fmt((1 - Math.cos(heel) ** 3) * 100)} % at rated thrust.`,
          `Active ballast that trims out the mean heel recovers that energy and reduces tower fatigue for the cost of pumps and controls; a 10 % smaller column spacing would cut bracing steel but lowers BM by about 19 % (BM ∝ spacing²).`,
          'Whole-life costing per ISO 15686-5 and LCOE convention: hull steel tonnage is the dominant floating-substructure CAPEX quantity, energy lost to tilt enters the denominator.'),
      ],
      publish: { 'vessel.disp': disp, 'vessel.GM': GM, 'vessel.heaveSig': heaveSig, 'vessel.rollSigDeg': op.st[3].significant / DEG, 'vessel.Tn3': Tn[2] ?? 0, 'vessel.Tn4': Tn[4] ?? 0, 'quantities.steel_t': steel, 'vessel.heelDeg': heel / DEG },
      quantities: { steel_t: round(steel, 4) },
    };
  },
  verify: [
    { name: 'Three-column semi-submersible: GM = KB + I_wp/∇ − KG and static heel = M/(ρg∇·GM)', source: 'Hydrostatics by the parallel-axis theorem (Biran & López-Pulido 2014, Ship Hydrostatics and Stability, ch. 2): hand calculation',
      // A_col = π·6² = 113.097 m², A_c = π·3.25² = 33.183 m²; I = 3·1017.88 + 113.097·1.5·28.87² + 87.62 = 144 537 m⁴
      // ∇ = 3·(113.097·14 + π·12²·6) + 33.183·20 = 13 556.8 m³; KB = [3·(1583.36·13 + 2714.34·3) + 663.66·10]/13 556.8 = 6.8465 m
      // BM = 10.6616 m; GM = 6.8465 + 10.6616 − 10.11 = 7.398 m; displacement = 1.025 × 13 556.8 = 13 895.7 t
      // heel = 800 kN × 104 m / (1025·9.80665·13 556.8·7.398) = 0.08253 rad = 4.728°
      inputs: { colR: 28.87, colD: 12, baseD: 24, baseH: 6, cenD: 6.5, draft: 20, KG: 10.11, thrust: 800, hubHeight: 90, fairleadDepth: 14 },
      expect: [{ kpi: 'GM', value: 7.398, tol: 0.002 }, { kpi: 'disp', value: 13895.7, tol: 0.001 }, { kpi: 'heel', value: 4.728, tol: 0.003 }] },
  ],
};
