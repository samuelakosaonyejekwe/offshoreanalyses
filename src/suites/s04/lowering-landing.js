// Lifting, lowering through the splash zone and water column, and landing on the seabed.
import { staticTension, axialResonance, heaveResponse, landingVelocity, waveKinematicsN103, splashZoneN103, currentOffset } from '../../engine/lines/lowering.js';
import { landingImpact, landingVelocityLimit } from '../../engine/geotech/mudmat.js';
import { G, DEG, need, run, fmt, round, adv, check, vesselCo2Text, CONFIRM } from './_common.js';

/** Touchdown as a single mass on an elastic–perfectly-plastic soil spring, velocity-Verlet in time. */
function touchdown({ M, W, k, Qult, v0, nSteps }) {
  const om = Math.sqrt(k / M), dt = (Math.PI / om) / nSteps; // half a natural period resolved by nSteps
  let x = 0, v = v0, Fmax = 0, xmax = 0;
  const Fs = (y) => (y > 0 ? Math.min(k * y, Qult) : 0);
  const ts = [0], Fh = [0], xh = [0];
  let a = (W - Fs(x)) / M;
  const every = Math.max(1, Math.floor(nSteps / 100));
  for (let i = 1; i <= 4 * nSteps; i++) {
    x += v * dt + 0.5 * a * dt * dt;
    const a2 = (W - Fs(x)) / M;
    v += 0.5 * (a + a2) * dt; a = a2;
    Fmax = Math.max(Fmax, Fs(x)); xmax = Math.max(xmax, x);
    if (i % every === 0) { ts.push(i * dt); Fh.push(Fs(x) / 1e3); xh.push(x * 1e3); }
    if (v <= 0) break;
  }
  return { Fmax, xmax, ts, Fh, xh, plastic: Fmax >= Qult * (1 - 1e-9) };
}

export default {
  id: 's04-lowering-landing', suite: 4, title: 'Lifting, lowering and landing',
  summary: 'Sling and hook loads in air, splash-zone forces and slack-sling check, line tension and axial resonance through the water column, current offset, and touchdown load and allowable set-down speed on the seabed.',
  covers: [13], fidelity: 'reduced-order',
  standards: ['DNV-RP-N103 (2017) §4 (splash zone), §5 (deep-water lowering), §6 (landing)', 'DNV-ST-N001 (2021) §16 (lifting operations)', 'DNV-RP-C205 (2019) App. D (added mass)', 'ISO 13628-1 / API RP 17A (installation)'],
  models: ['Archimedes submerged weight', 'Lowering-line axial dynamics (rod with end mass)', 'Single-degree-of-freedom heave response with linearised drag', 'DNV-RP-N103 Simplified Method for the splash zone', 'Added-mass model', 'Seabed touchdown contact (elastic–plastic soil spring, time domain)', 'Energy-based landing impact', 'Cable offset in current'],
  limits: 'Vertical single-degree-of-freedom response of the load on an elastic wire with regular crane-tip motion; the DNV-RP-N103 Simplified Method expressions are recalled and their coefficients are inputs to be confirmed against the governing edition. Coupled vessel–load multibody dynamics, wave-frequency response in irregular seas, pendulum motion in air and water-exit of entrapped air are not modelled; the soil at touchdown is a single elastic–plastic spring.',
  publishes: ['install.liftWeight', 'install.hookLoad', 'install.vesselDays', 'subsea.landingLoad'],
  searchTerms: 'subsea lift lowering splash zone DNV-RP-N103 landing velocity mudmat',
  inputs: [
    { key: 'mass', label: 'Mass of the lifted structure', unit: 't', type: 'number', default: 120, min: 1, max: 20000, group: 'Lifted object', role: 'ic', bind: 'subsea.mass', uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'vol', label: 'Displaced volume', unit: 'm³', type: 'number', default: 18, min: 0, max: 20000, group: 'Lifted object', role: 'ic', bind: 'subsea.volume' },
    { key: 'A33', label: 'Heave added mass (incl. entrapped water)', unit: 't', type: 'number', default: 350, min: 0, max: 1e5, group: 'Lifted object', role: 'load', uncertain: { dist: 'lognormal', cov: 0.25 },
      help: 'For a solid rectangular plate a × b: ρ·CA·(π/4)·a²·b with CA ≈ 0.58–1.0 (DNV-RP-C205 App. D); perforation of 15–25 % typically halves it.' },
    { key: 'Aplan', label: 'Plan area (drag, slamming, waterplane)', unit: 'm²', type: 'number', default: 90, min: 0.5, max: 5000, group: 'Lifted object', role: 'geometry', help: 'Solid projected area on a horizontal plane.' },
    { key: 'CdZ', label: 'Vertical drag coefficient', unit: '-', type: 'number', default: 2.5, min: 0.5, max: 8, group: 'Lifted object', role: 'criteria', help: `DNV-RP-N103 §4.3.8: not less than 2.5 for typical subsea structures in oscillatory flow unless documented. ${CONFIRM}` },
    { key: 'Cs', label: 'Slamming coefficient', unit: '-', type: 'number', default: 5.0, min: 1, max: 10, group: 'Lifted object', role: 'criteria', help: `DNV-RP-N103 §4.5: not less than 5.0 for structures of complex shape. ${CONFIRM}` },
    { key: 'zetaFac', label: 'Characteristic wave amplitude / Hs', unit: '-', type: 'number', default: 0.9, min: 0.5, max: 1.2, group: 'Lifted object', role: 'criteria', help: `DNV-RP-N103 §4.3.2: ζa = 0.9·Hs for operations shorter than 30 minutes (1.0·Hs otherwise). ${CONFIRM}` },
    { key: 'slackFac', label: 'Slack-sling margin on the static load', unit: '-', type: 'number', default: 0.9, min: 0.5, max: 1, group: 'Lifted object', role: 'criteria', help: `DNV-RP-N103 §4.4.3: Fhyd ≤ 0.9·Fstatic,min. ${CONFIRM}` },
    { key: 'slamFrac', label: 'Slamming area / plan area', unit: '-', type: 'number', default: 0.3, min: 0, max: 1, group: 'Lifted object', role: 'geometry', help: 'Share of the plan area that meets the surface at once (flat bottom plating).' },
    { key: 'nSling', label: 'Number of slings', type: 'int', default: 4, min: 1, max: 8, group: 'Rigging', role: 'operation' },
    { key: 'slingAng', label: 'Sling angle from vertical', unit: 'deg', type: 'number', default: 30, min: 0, max: 60, group: 'Rigging', role: 'operation' },
    { key: 'skl', label: 'Skew load factor', unit: '-', type: 'number', default: 1.25, min: 1, max: 2, group: 'Rigging', role: 'criteria', help: `DNV-ST-N001: 1.25 for a statically indeterminate four-sling lift. ${CONFIRM}` },
    { key: 'slingMBL', label: 'Sling minimum breaking load', unit: 'kN', type: 'number', default: 2500, min: 10, max: 1e5, group: 'Rigging', role: 'material', design: true },
    { key: 'slingSF', label: 'Sling nominal safety factor', unit: '-', type: 'number', default: 3.0, min: 1.5, max: 6, group: 'Rigging', role: 'criteria', help: `DNV-ST-N001 nominal safety factor for steel-wire slings is of this order. ${CONFIRM}` },
    { key: 'dafAir', label: 'Dynamic amplification factor in air', unit: '-', type: 'number', default: 1.15, min: 1, max: 2, group: 'Rigging', role: 'criteria', help: `DNV-ST-N001 Table for offshore lifts 100–300 t. ${CONFIRM}` },
    { key: 'craneCap', label: 'Crane capacity at the lift radius', unit: 'kN', type: 'number', default: 2500, min: 10, max: 1e6, group: 'Crane and wire', role: 'criteria' },
    { key: 'EA', label: 'Wire axial stiffness EA', unit: 'MN', type: 'number', default: 280, min: 5, max: 5000, group: 'Crane and wire', role: 'material', uncertain: { dist: 'normal', cov: 0.1 } },
    { key: 'wireMass', label: 'Wire mass per metre', unit: 'kg/m', type: 'number', default: 26, min: 0.5, max: 300, group: 'Crane and wire', role: 'material' },
    { key: 'wireD', label: 'Wire diameter', unit: 'mm', type: 'number', default: 76, min: 10, max: 250, group: 'Crane and wire', role: 'geometry' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 150, min: 5, max: 3500, group: 'Site', role: 'environment', bind: 'site.depth' },
    { key: 'Hs', label: 'Installation significant wave height', unit: 'm', type: 'number', default: 1.5, min: 0, max: 6, group: 'Operation limits', role: 'bc', design: true, help: 'Operational limiting sea state (design criterion OPLIM).' },
    { key: 'Tz', label: 'Zero-up-crossing period', unit: 's', type: 'number', default: 6.5, min: 2, max: 18, group: 'Operation limits', role: 'bc' },
    { key: 'ctAmp', label: 'Crane-tip heave amplitude', unit: 'm', type: 'number', default: 0.7, min: 0, max: 6, group: 'Operation limits', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Prescribed single-amplitude vertical motion of the crane tip in the limiting sea state, from the vessel RAOs.' },
    { key: 'ctT', label: 'Crane-tip motion period', unit: 's', type: 'number', default: 8, min: 2, max: 25, group: 'Operation limits', role: 'bc' },
    { key: 'vc', label: 'Lowering speed', unit: 'm/s', type: 'number', default: 0.5, min: 0, max: 2, group: 'Operation limits', role: 'ic' },
    { key: 'vLand', label: 'Lowering speed at set-down', unit: 'm/s', type: 'number', default: 0.1, min: 0, max: 2, group: 'Operation limits', role: 'ic', design: true, help: 'Winch pay-out speed in the last metres; heave compensation reduces the superimposed motion.' },
    { key: 'ahc', label: 'Residual motion with heave compensation', unit: '-', type: 'number', default: 0.3, min: 0, max: 1, group: 'Operation limits', role: 'operation', help: '1.0 without compensation; 0.1–0.3 with active heave compensation at landing.' },
    { key: 'Uc', label: 'Current during lowering (uniform)', unit: 'm/s', type: 'number', default: 0.5, min: 0, max: 3, group: 'Operation limits', role: 'bc', bind: 'site.Uc1' },
    { key: 'AdragH', label: 'Horizontal drag area Cd·A of the object', unit: 'm²', type: 'number', default: 40, min: 0, max: 5000, group: 'Operation limits', role: 'load' },
    { key: 'ksoil', label: 'Vertical soil stiffness at touchdown', unit: 'MN/m', type: 'number', default: 50, min: 0.1, max: 1e5, group: 'Seabed', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.4 }, help: 'Subgrade modulus × base area, or secant stiffness from the bearing-capacity curve.' },
    { key: 'Qult', label: 'Ultimate bearing capacity of the seabed', unit: 'kN', type: 'number', default: 6000, min: 1, max: 1e7, group: 'Seabed', role: 'soil', bind: 'subsea.Vult' },
    { key: 'fsLand', label: 'Factor on bearing at landing', unit: '-', type: 'number', default: 1.5, min: 1, max: 4, group: 'Seabed', role: 'criteria', help: CONFIRM },
    { key: 'tolLand', label: 'Allowable landing position tolerance', unit: 'm', type: 'number', default: 2.5, min: 0.05, max: 50, group: 'Seabed', role: 'criteria', help: 'Target-box half-width; the current offset must be compensated by vessel positioning within it.' },
    { key: 'nSteps', label: 'Time steps per half period of touchdown', type: 'int', default: 400, min: 20, max: 20000, group: 'Solver', role: 'control' },
  ],
  meshTime: { key: 'nSteps', label: 'Time steps', levels: [50, 100, 200, 400, 800], qoi: 'FlandTD' },
  compute(inp) {
    // SI: t → kg, kN → N, MN → N, mm → m, deg → rad.
    const rho = 1025, m = inp.mass * 1e3, A33 = inp.A33 * 1e3, EA = inp.EA * 1e6, d = inp.depth;
    const Wair = m * G, Wsub = (m - rho * inp.vol) * G;
    need(Wsub > 0, 'The structure is buoyant: it cannot be lowered on a wire without ballast.');
    const wSub = inp.wireMass * G * (1 - rho / 7850);
    // --- in air: hook and sling loads
    const hookAir = inp.dafAir * Wair;
    const sling = hookAir * inp.skl / (inp.nSling * Math.cos(inp.slingAng * DEG));
    // --- splash zone (object just submerged)
    const wk = run(() => waveKinematicsN103({ Hs: Math.max(inp.Hs, 1e-6), Tz: inp.Tz, d: 1.0, amplitudeFactor: inp.zetaFac }), 'Wave kinematics');
    const omC = 2 * Math.PI / inp.ctT, vCt = inp.ctAmp * omC, aCt = inp.ctAmp * omC * omC;
    const sp = run(() => splashZoneN103({ mass: m, volume: inp.vol, addedMass: A33, slamArea: inp.slamFrac * inp.Aplan, Cs: inp.Cs, dragArea: inp.Aplan, CD: inp.CdZ, waterplaneArea: 0.1 * inp.Aplan,
      vc: inp.vc, etaCt: inp.ctAmp, vCt, aCt, zetaA: inp.Hs > 0 ? wk.zetaA : 0, vW: inp.Hs > 0 ? wk.velocity : 0, aW: inp.Hs > 0 ? wk.acceleration : 0, staticMin: Wsub, staticMax: Wair, rho }), 'Splash zone');
    const slackUc = sp.staticMin > 0 ? sp.Fhyd / (inp.slackFac * sp.staticMin) : 99;
    // --- deep lowering: static tension, resonance and dynamic response along the depth
    const depths = Array.from({ length: 30 }, (_, i) => Math.max(d * (i + 1) / 30, 1));
    const st = run(() => staticTension({ mass: m, volume: inp.vol, cableWeight: wSub, depths, rho }), 'Static tension');
    const prof = depths.map((Lw, i) => {
      const res = axialResonance({ mass: m, addedMass: A33, EA, length: Lw, cableMass: inp.wireMass });
      let hr = null;
      if (inp.ctAmp > 0) hr = heaveResponse({ mass: m, addedMass: A33, EA, length: Lw, cableMass: inp.wireMass, heaveAmp: inp.ctAmp, period: inp.ctT, Cd: inp.CdZ, area: inp.Aplan, rho, staticTension: Wsub });
      return { L: Lw, Ttop: st.rows[i].tension, T0: res.periodLumped, dyn: hr ? hr.dynamicForce : 0, minT: hr ? hr.minTension : Wsub, amp: hr ? hr.motionAmp : 0, slack: hr ? hr.slack : false };
    });
    const worstDyn = prof.reduce((a, p) => (p.dyn > a.dyn ? p : a));
    const atBed = prof[prof.length - 1];
    const TtopMax = Math.max(...prof.map((p) => p.Ttop + p.dyn));
    const hookMax = Math.max(hookAir, sp.Ftotal, TtopMax);
    const minT = Math.min(...prof.map((p) => p.minT));
    // --- current offset at full depth
    const off = inp.Uc > 0 ? run(() => currentOffset({ length: d, objectWeight: Wsub, cableWeight: wSub, cableDiameter: inp.wireD / 1e3, objectDragArea: inp.AdragH, current: inp.Uc, rho }), 'Current offset') : { offset: 0, topAngleDeg: 0 };
    // --- landing
    const lv = landingVelocity({ loweringSpeed: inp.vLand, objectMotionAmp: atBed.amp * inp.ahc, period: inp.ctT });
    const k = inp.ksoil * 1e6, Qult = inp.Qult * 1e3, Qallow = Qult / inp.fsLand;
    const li = landingImpact({ mass: m, addedMass: A33, W: Wsub, v: lv.velocity, k });
    const vLim = landingVelocityLimit({ mass: m, addedMass: A33, W: Wsub, k, Qallow });
    const td = touchdown({ M: m + A33, W: Wsub, k, Qult, v0: lv.velocity, nSteps: inp.nSteps });
    const hours = 4 + d / Math.max(inp.vc, 0.05) / 3600 + 3; // rigging and overboarding, lowering, landing and release
    const days = hours / 24;
    const checks = [
      check('Crane capacity', hookMax / 1e3, inp.craneCap, 'kN', 'DNV-ST-N001 (2021) §16 (dynamic hook load); DNV-RP-N103 §4.4.4', 'Guards against overloading the crane at any stage: in air, in the splash zone or at the depth of largest dynamic tension.'),
      check('Sling load', sling / 1e3, inp.slingMBL / inp.slingSF, 'kN', 'DNV-ST-N001 (2021) §16 (sling design load with skew load factor)', 'Guards against sling failure. The sling load includes dynamic amplification, the skew load factor and the sling angle.'),
      { name: 'Slack sling in the splash zone', demand: round(sp.Fhyd / 1e3), capacity: round(inp.slackFac * sp.staticMin / 1e3), unit: 'kN', uc: Math.min(slackUc, 99), ref: `DNV-RP-N103 (2017) §4.4.3: Fhyd ≤ ${inp.slackFac}·Fstatic,min`, explain: 'Guards against the slings going slack as waves lift the structure, which is followed by a snap load when they re-tension.' },
      { name: 'Slack wire during deep lowering', demand: round(worstDyn.dyn / 1e3), capacity: round(inp.slackFac * Wsub / 1e3), unit: 'kN', uc: Math.min(worstDyn.dyn / (inp.slackFac * Wsub), 99), ref: 'DNV-RP-N103 (2017) §5.3 (dynamic force against static submerged weight)', explain: 'Guards against the dynamic tension amplitude exceeding the submerged weight near axial resonance, where the wire goes slack and snaps.' },
      check('Landing reaction against seabed bearing', td.Fmax / 1e3, Qallow / 1e3, 'kN', 'DNV-RP-N103 (2017) §6.2 (landing on seabed); ISO 19901-4 §7', 'Guards against bearing failure and excessive penetration at touchdown, when kinetic energy of structure and added mass is absorbed by the soil.'),
      check('Landing position offset by current', off.offset, inp.tolLand, 'm', 'DNV-RP-N103 (2017) §5.2.2 (horizontal offset); project installation tolerance', 'The current pushes the load away from the crane tip. Beyond the tolerance the vessel must be re-positioned or the operation delayed.'),
    ];
    const ok = checks.every((c) => c.uc <= 1);
    const w = checks.reduce((a, c) => (c.uc > a.uc ? c : a));
    const resDepth = prof.find((p) => p.T0 >= inp.ctT);
    const advice = [
      ok ? adv('engineering', 'ok', `Maximum hook load ${fmt(hookMax / 1e3)} kN of ${fmt(inp.craneCap)} kN; splash-zone hydrodynamic force ${fmt(sp.Fhyd / 1e3)} kN (slack criterion ${fmt(slackUc)}); landing at ${fmt(lv.velocity)} m/s gives ${fmt(td.Fmax / 1e3)} kN against ${fmt(Qallow / 1e3)} kN allowable.`,
        `Operate at Hs ≤ ${fmt(inp.Hs)} m with set-down speed below ${fmt(vLim)} m/s. ${resDepth ? `Axial resonance with the crane-tip period occurs near ${fmt(resDepth.L)} m of wire: pass that depth without stopping.` : 'The axial natural period stays below the crane-tip period over the whole depth.'} Confirm the splash-zone coefficients against the governing edition of DNV-RP-N103.`,
        'DNV-RP-N103 (2017) §4–6; DNV-ST-N001 (2021) §16',
        `A workable limiting sea state avoids waiting on weather: each day of standby is ${vesselCo2Text(1)}`)
        : adv('engineering', 'fail', `"${w.name}" is exceeded (${fmt(w.uc)}): demand ${fmt(w.demand)} against ${fmt(w.capacity)} ${w.unit}.`,
          w.name.startsWith('Slack sling') ? `Reduce the limiting sea state (hydrodynamic force scales roughly with Hs to the power 1–2), increase perforation to cut added mass and slamming area, or lower faster through the surface (now ${fmt(inp.vc)} m/s).`
            : w.name.startsWith('Landing') ? `Reduce the set-down speed to below ${fmt(vLim)} m/s (heave compensation), or land at reduced weight by keeping hook load on until the skirts are home.`
              : w.name.startsWith('Slack wire') ? 'Avoid dwelling near the resonance depth, use a softer passive compensator to shift the natural period, or reduce added mass with perforation.'
                : w.name.startsWith('Sling') ? `Use slings of at least ${fmt(sling * inp.slingSF / 1e3)} kN breaking load or reduce the sling angle from vertical.` : `Use a crane of at least ${fmt(hookMax / 1e3)} kN at this radius or reduce the lift weight by installing modules separately.`,
          'DNV-RP-N103 (2017) §4–6; DNV-ST-N001 (2021) §16',
          w.name.startsWith('Crane') || w.name.startsWith('Sling') ? `A larger crane vessel burns more fuel per day — ${vesselCo2Text(1)} Splitting the lift adds a second operation of similar duration.` : `Each day of weather standby caused by a lower sea-state limit: ${vesselCo2Text(1)} Perforating the mat cuts added mass and slamming area and removes a little steel as well.`),
      adv('economic', ok ? 'info' : 'warn', `The operation takes about ${fmt(hours)} h of vessel time (${fmt(days)} days) at ${fmt(inp.vc)} m/s lowering speed, excluding transit and waiting on weather; the limiting sea state is Hs ${fmt(inp.Hs)} m.`,
        'Operability drives cost more than duration: raising the installable Hs by 0.5 m typically adds tens of percentage points of workable time in winter. Spend engineering on perforation, heave compensation and landing speed before accepting a low limit.',
        'DNV-ST-N001 (2021) §2 operational limits and weather windows (alpha factors); AACE RP 18R-97'),
    ];
    return {
      kpis: [
        { key: 'hookMax', label: 'Maximum dynamic hook load', unit: 'kN', value: hookMax / 1e3, limit: inp.craneCap, sense: 'max', explain: 'Largest hook load over the three phases: in air with dynamic factor, in the splash zone, and at the depth of largest dynamic line tension.' },
        { key: 'Wair', label: 'Lift weight in air', unit: 'kN', value: Wair / 1e3, explain: 'Static weight of the structure in air, the basis for crane and rigging selection.' },
        { key: 'Wsub', label: 'Submerged weight', unit: 'kN', value: Wsub / 1e3, explain: 'Weight in water by Archimedes; it is the static tension at the hook once submerged and the load landed on the seabed.' },
        { key: 'sling', label: 'Sling design load', unit: 'kN', value: sling / 1e3, limit: inp.slingMBL / inp.slingSF, sense: 'max', explain: 'Load per sling with dynamic amplification, skew load factor and sling angle.' },
        { key: 'Fhyd', label: 'Splash-zone hydrodynamic force', unit: 'kN', value: sp.Fhyd / 1e3, explain: 'Characteristic hydrodynamic force combining slamming, drag, inertia and varying buoyancy by the DNV-RP-N103 Simplified Method.' },
        { key: 'slackUc', label: 'Slack-sling criterion', unit: '-', value: Math.min(slackUc, 99), limit: 1, sense: 'max', explain: 'Hydrodynamic force over the margin factor (0.9 by default) × minimum static load. Above 1 the slings can go slack and snap loads must be assessed.' },
        { key: 'TtopSeabed', label: 'Static top tension at the seabed', unit: 'kN', value: atBed.Ttop / 1e3, explain: 'Submerged weight of the structure plus the submerged wire paid out to the seabed.' },
        { key: 'T0bed', label: 'Axial natural period at the seabed', unit: 's', value: atBed.T0, explain: 'Natural period of the load bouncing on the wire at full depth. Close to the crane-tip period means resonance and large dynamic tension.' },
        { key: 'dynMax', label: 'Largest dynamic line force', unit: 'kN', value: worstDyn.dyn / 1e3, explain: 'Largest dynamic tension amplitude over the lowering depth; it must stay below the submerged weight to keep the wire taut.' },
        { key: 'minT', label: 'Minimum line tension', unit: 'kN', value: minT / 1e3, limit: 0, sense: 'min', explain: 'Lowest tension at the object over the depth; zero or negative means a slack wire.' },
        { key: 'offset', label: 'Horizontal offset in current', unit: 'm', value: off.offset, explain: 'Drift of the load relative to the crane tip at full depth; the vessel is positioned up-current by this amount.' },
        { key: 'vTouch', label: 'Touchdown velocity', unit: 'm/s', value: lv.velocity, limit: vLim || 1e-9, sense: 'max', explain: 'Set-down speed plus the residual heave velocity of the load at the seabed.' },
        { key: 'vLim', label: 'Allowable touchdown velocity', unit: 'm/s', value: vLim, explain: 'Largest set-down velocity for which the landing reaction stays within the allowable bearing load (energy balance).' },
        { key: 'Fland', label: 'Landing reaction (energy method)', unit: 'kN', value: li.Fmax / 1e3, explain: 'Peak soil reaction from the closed-form energy balance on a linear soil spring.' },
        { key: 'FlandTD', label: 'Landing reaction (time domain)', unit: 'kN', value: td.Fmax / 1e3, limit: Qallow / 1e3, sense: 'max', explain: 'Peak soil reaction from the time-domain touchdown simulation, capped at the ultimate bearing capacity when the soil yields.' },
        { key: 'pen', label: 'Touchdown penetration', unit: 'mm', value: td.xmax * 1e3, explain: 'Largest vertical movement into the soil at touchdown; large values mean plastic bearing failure during landing.' },
        { key: 'days', label: 'Vessel time for the operation', unit: 'days', value: days, explain: 'Rigging, overboarding, lowering and landing time at the stated lowering speed, without transit or weather standby.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Line tension against wire length', x: { label: 'Tension at the crane tip', unit: 'kN' }, y: { label: 'Depth of the load', unit: 'm' },
          series: [{ name: 'Static', x: prof.map((p) => p.Ttop / 1e3), y: prof.map((p) => -p.L) }, { name: 'Static + dynamic', x: prof.map((p) => (p.Ttop + p.dyn) / 1e3), y: prof.map((p) => -p.L) }, { name: 'Static − dynamic', x: prof.map((p) => (p.Ttop - p.dyn) / 1e3), y: prof.map((p) => -p.L) }],
          explain: 'Envelope of top tension as the load descends. A widening envelope marks the approach to axial resonance; the lower branch reaching the wire weight line means slack at the object.' },
        { type: 'line', title: 'Axial natural period against depth', x: { label: 'Natural period', unit: 's' }, y: { label: 'Depth of the load', unit: 'm' },
          series: [{ name: 'Natural period', x: prof.map((p) => p.T0), y: prof.map((p) => -p.L) }, { name: 'Crane-tip period', x: [inp.ctT, inp.ctT], y: [0, -d], dash: true }],
          explain: 'The natural period lengthens with paid-out wire. Where it crosses the crane-tip period the load resonates; do not hold the load at that depth.' },
        { type: 'line', title: 'Touchdown: soil reaction in time', x: { label: 'Time after contact', unit: 's' }, y: { label: 'Soil reaction', unit: 'kN' },
          series: [{ name: 'Reaction', x: td.ts, y: td.Fh }, { name: 'Allowable', x: [0, td.ts[td.ts.length - 1]], y: [Qallow / 1e3, Qallow / 1e3], dash: true }, { name: 'Submerged weight', x: [0, td.ts[td.ts.length - 1]], y: [Wsub / 1e3, Wsub / 1e3], dash: true }],
          explain: 'Build-up of soil reaction as the structure decelerates. The peak above the submerged weight is the landing amplification; a flat top means the soil has reached its ultimate capacity.' },
        { type: 'bar', title: 'Splash-zone force components', categories: ['Slamming', 'Varying buoyancy', 'Mass (inertia)', 'Drag', 'Combined'], series: [{ name: 'Force (kN)', y: [sp.Fslam, sp.Frho, sp.FM, sp.FD, sp.Fhyd].map((v) => round(v / 1e3)) }], y: { label: 'Force', unit: 'kN' },
          explain: 'Components of the Simplified Method. The largest bar shows what to change: slamming area, perforation (added mass, drag) or sea state.' },
      ],
      tables: [{ title: 'Lowering profile', columns: [{ key: 'L', label: 'Wire in water', unit: 'm' }, { key: 'T', label: 'Static top tension', unit: 'kN' }, { key: 'T0', label: 'Natural period', unit: 's' }, { key: 'dyn', label: 'Dynamic force', unit: 'kN' }, { key: 'amp', label: 'Load motion amplitude', unit: 'm' }, { key: 'slack', label: 'Slack' }],
        rows: prof.filter((_, i) => i % 3 === 2).map((p) => ({ L: round(p.L), T: round(p.Ttop / 1e3), T0: round(p.T0), dyn: round(p.dyn / 1e3), amp: round(p.amp), slack: p.slack ? 'yes' : 'no' })),
        explain: 'Static tension, axial natural period and dynamic response at ten depths during lowering. Use it to mark hold points and no-stop zones in the installation procedure.' }],
      notes: [`"${w.name}" is the governing check (${fmt(w.uc)}).`, td.plastic ? 'The soil reaches its ultimate capacity at touchdown: penetration is plastic.' : 'The soil stays elastic at touchdown.', `Crane-tip velocity ${fmt(vCt)} m/s and acceleration ${fmt(aCt)} m/s²; wave amplitude ζa = ${fmt(inp.Hs > 0 ? wk.zetaA : 0)} m.`],
      advice,
      publish: { 'install.liftWeight': Wair, 'install.hookLoad': hookMax, 'install.vesselDays': days, 'subsea.landingLoad': td.Fmax },
      quantities: { vessel_days: { installation: round(days) } },
    };
  },
  verify: [
    { name: 'Static top tension: T = (m − ρV)g + w′·L', source: 'Hand calculation: (100 000 − 1025·10)·9.80665 = 880 147 N; wire 20 kg/m → w′ = 20·9.80665·(1 − 1025/7850) = 170.52 N/m; at 1000 m T = 1050.67 kN',
      inputs: { mass: 100, vol: 10, wireMass: 20, depth: 1000, ctAmp: 0, Uc: 0 }, expect: [{ kpi: 'TtopSeabed', value: 1050.671, tol: 1e-5 }, { kpi: 'Wsub', value: 880.147, tol: 1e-5 }] },
    { name: 'Landing on a linear soil spring: F = W′ + √(W′² + k(m + ma)v²)', source: 'Energy balance ½(m+ma)v² + W′δ = ½kδ² (cf. DNV-RP-N103 §6.2): W′ = 880 147 N, k = 50 MN/m, m + ma = 300 t, v = 0.5 m/s → F = 880 147 + √(7.7466e11 + 3.75e12) = 3007.3 kN',
      inputs: { mass: 100, vol: 10, A33: 200, ksoil: 50, vLand: 0.5, ctAmp: 0, Qult: 1e6, nSteps: 800, Uc: 0 },
      expect: [{ kpi: 'Fland', value: 3007.28, tol: 1e-4 }, { kpi: 'FlandTD', value: 3007.28, tol: 0.002 }] },
    { name: 'Lumped axial natural period: T0 = 2π√((M + A33 + m_wire·L/3)·L/EA)', source: 'DNV-RP-N103 §5.3.3 with θ = 1/3: M + A33 = 470 t, wire 26 kg/m × 150 m, EA = 280 MN',
      inputs: { mass: 120, A33: 350, EA: 280, wireMass: 26, depth: 150 }, expect: [{ kpi: 'T0bed', value: 2 * Math.PI * Math.sqrt((470e3 + 26 * 150 / 3) * 150 / 280e6), tol: 1e-3 }] },
  ],
};
