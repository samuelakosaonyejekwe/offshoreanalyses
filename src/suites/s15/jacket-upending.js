// Jacket upending by controlled flooding and hook assistance: planar rigid-body time simulation (heave + pitch).
import { rk4 } from '../../engine/math/ode.js';
import { G, DEG, need, round, fmt, kpi, check, adv, line, tbl, thin, worst, vesselCo2Text, LIMIT_ROM } from './_common.js';
import { jacketInputs, buildJacket, hydro, addForce, addPointMass, solve3 } from './_jacket.js';

export default {
  id: 's15-jacket-upending', suite: 15, title: 'Jacket upending',
  summary: 'Simulates a floating jacket rotating to the vertical as its lower legs are flooded and the crane holds the top, and reports the rotation history, hook load, bottom clearance, set-down or floating draft and upright stability.',
  covers: [3], fidelity: 'reduced-order',
  standards: ['DNV-ST-N001 (2021) Sec. 13 (upending: stability, clearance, buoyancy and flooding contingencies)', 'ISO 19901-6 (2009) §15–16 (upending and set-down)', 'API RP 2A-WSD (2014) §2.4.5 (upending forces)'],
  models: ['Newton–Euler planar rigid-body dynamics (heave and pitch)', 'Controlled flooding as a time-varying point mass', 'Time-varying buoyancy by strip integration', 'Hydrostatic equilibrium and metacentric stability of the upright jacket', 'Added mass and quadratic hydrodynamic damping (Morison strips)', 'Prescribed hook load', 'Seabed contact (penalty spring)'],
  limits: `${LIMIT_ROM} Flood water is a point mass at one station with a prescribed rate, the hook load is prescribed rather than solved with the crane wire, surge and out-of-plane motion are ignored, and the jacket buoyancy is a linear distribution along the axis. Compartment-by-compartment flooding with free surfaces and the three-dimensional rotation of an asymmetric jacket need a dedicated upending analysis.`,
  publishes: ['install.upendHookLoad', 'install.upendFloodMass'],
  searchTerms: 'jacket upending controlled flooding crane assisted simulation stability',
  inputs: [
    ...jacketInputs({ M: 14000, L: 165, xiG: 92, xiB: 76, reserve: 30 }),
    { key: 'th0', label: 'Initial inclination of the axis', unit: 'deg', type: 'number', default: 2, min: -30, max: 90, group: 'Initial conditions', role: 'ic', help: 'Floating attitude after launch (lower end down positive). The initial draft is found from vertical equilibrium.' },
    { key: 'mFlood', label: 'Flood water admitted', unit: 't', type: 'number', default: 1200, min: 0, max: 60000, group: 'Flooding', role: 'bc', design: true, uncertain: { dist: 'normal', cov: 0.05 }, help: 'Total water let into the lower leg compartments by the end of the flooding sequence.' },
    { key: 'xiF', label: 'Station of the flooded compartments from the upper end', unit: 'm', type: 'number', default: 158, min: 0, max: 500, group: 'Flooding', role: 'bc', help: 'Centroid of the flood water along the axis; close to the lower end for bottom leg compartments.' },
    { key: 'qFlood', label: 'Flooding rate', unit: 't/min', type: 'number', default: 200, min: 1, max: 20000, group: 'Flooding', role: 'bc', help: 'Controlled by the flood valves; real sequences take 20–60 min, a faster rate here shortens the simulation without changing the end state.' },
    { key: 'hook', label: 'Hook load held at the upper end', unit: 't', type: 'number', default: 0, min: 0, max: 20000, group: 'Crane', role: 'bc', help: 'Prescribed vertical force at the hook attachment (constant-tension mode). Zero for a self-upending jacket.' },
    { key: 'xiH', label: 'Hook attachment from the upper end', unit: 'm', type: 'number', default: 5, min: 0, max: 500, group: 'Crane', role: 'bc' },
    { key: 'hookCap', label: 'Crane capacity at the upending radius', unit: 't', type: 'number', default: 3000, min: 0, max: 20000, group: 'Crane', role: 'criteria' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 150, min: 10, max: 3000, group: 'Site', role: 'environment', bind: 'site.depth' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'clearReq', label: 'Required bottom clearance while rotating', unit: 'm', type: 'number', default: 5, min: 0, max: 100, group: 'Criteria', role: 'criteria', help: 'Applies until the jacket is within the upright tolerance; afterwards touching down is the intended set-down. Confirm the value against DNV-ST-N001 and the project procedure.' },
    { key: 'tolDeg', label: 'Upright tolerance', unit: 'deg', type: 'number', default: 3, min: 0.2, max: 20, group: 'Criteria', role: 'criteria' },
    { key: 'gmReq', label: 'Required upright GM', unit: 'm', type: 'number', default: 1, min: 0, max: 20, group: 'Criteria', role: 'criteria', help: 'Minimum metacentric height of the free-floating upright jacket before set-down; DNV-ST-N001 asks for positive stability with margin in every stage.' },
    { key: 'tEnd', label: 'Simulated time', unit: 's', type: 'number', default: 600, min: 60, max: 7200, group: 'Numerics', role: 'control' },
    { key: 'rate', label: 'Time steps per second', unit: '1/s', type: 'int', default: 8, min: 2, max: 200, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'rate', label: 'Time steps per second', levels: [4, 8, 16, 32], qoi: 'thEnd' },
  compute(inp) {
    const j = buildJacket(inp, inp.rhoW);
    const mF = inp.mFlood * 1e3, q = inp.qFlood * 1e3 / 60, xiF = inp.xiF - j.La, xiH = inp.xiH - j.La, H = inp.hook * 1e3 * G;
    need(inp.xiF <= j.L && inp.xiH <= j.L, 'The flooding and hook stations must lie within the jacket length.');
    const zBed = -inp.depth, kBed = j.M * G / 0.1 / 4, cBed = 0.8 * Math.sqrt(kBed * j.M);
    const flood = (t) => Math.min(q * t, mF);
    const rhs = (t, y, out) => {
      const z = y[0], th = y[1], vz = y[2], om = y[3];
      const Q = [0, -j.M * G, 0], Mm = [j.M, 0, 0, 0, j.M, 0, 0, 0, j.I];
      const Vsub = hydro(j, z, th, 0, vz, om, Q, Mm);
      addPointMass(Q, Mm, th, om, xiF, flood(t));
      if (H > 0) addForce(Q, th, xiH, 0, H);
      // Seabed under either end.
      let Rbed = 0;
      const s = Math.sin(th), c = Math.cos(th);
      for (const xi of [-j.La, j.Lf]) {
        const pz = z - xi * s, pen = zBed - pz;
        if (pen > 0) { const N = Math.max(kBed * pen + cBed * -(vz - xi * om * c), 0); addForce(Q, th, xi, 0, N); Rbed += N; }
      }
      if (out) { out.Rbed = Rbed; out.Vsub = Vsub; }
      // Surge is not a degree of freedom: solve heave and pitch only.
      const a = Mm[4], b = Mm[5], d = Mm[8], det = a * d - b * b;
      return [(Q[1] * d - b * Q[2]) / det, (a * Q[2] - b * Q[1]) / det];
    };
    // Initial draft: vertical equilibrium at the initial inclination (bisection on z).
    const th0 = inp.th0 * DEG;
    const net = (z) => { const Q = [0, -j.M * G, 0], Mm = new Array(9).fill(0); hydro(j, z, th0, 0, 0, 0, Q, Mm); return Q[1] + H; };
    let lo = -j.L - 50, hi = j.L + 50;
    need(net(lo) > 0, 'The jacket does not float: buoyancy plus hook load is less than its weight. Increase the reserve buoyancy or the hook load.');
    for (let i = 0; i < 80; i++) { const m = 0.5 * (lo + hi); if (net(m) > 0) lo = m; else hi = m; }
    const nSteps = Math.round(inp.tEnd * inp.rate);
    const sol = rk4((t, y) => { const a = rhs(t, y); return [y[2], y[3], a[0], a[1]]; }, 0, [0.5 * (lo + hi), th0, 0, 0], inp.tEnd, nSteps);
    const n = sol.t.length, tt = sol.t, o = {};
    const thD = new Float64Array(n), zLow = new Float64Array(n), zTop = new Float64Array(n), mf = new Float64Array(n), rb = new Float64Array(n);
    let clearMin = Infinity, omMax = 0, tUp = NaN, RbedMax = 0, mUp = NaN;
    const tol = inp.tolDeg * DEG;
    for (let k = 0; k < n; k++) {
      const y = sol.y[k], th = y[1];
      need(Number.isFinite(y[0]) && Number.isFinite(th), 'The upending simulation diverged. Increase the number of time steps per second.');
      rhs(tt[k], y, o);
      const s = Math.sin(th);
      thD[k] = th / DEG; zLow[k] = Math.min(y[0] - j.Lf * s, y[0] + j.La * s); zTop[k] = y[0] + j.La * s; mf[k] = flood(tt[k]) / 1e3; rb[k] = o.Rbed;
      const upright = Math.abs(th - Math.PI / 2) <= tol;
      if (upright && Number.isNaN(tUp)) { tUp = tt[k]; mUp = mf[k]; }
      if (!upright) clearMin = Math.min(clearMin, zLow[k] - zBed);
      omMax = Math.max(omMax, Math.abs(y[3]) / DEG); RbedMax = Math.max(RbedMax, o.Rbed);
    }
    if (!Number.isFinite(clearMin)) clearMin = zLow[0] - zBed;
    // End state: mean over the last fifth of the record to remove residual oscillation.
    const k0 = Math.floor(0.8 * (n - 1)); let thE = 0, topE = 0, lowE = 0, rbE = 0;
    for (let k = k0; k < n; k++) { thE += thD[k]; topE += zTop[k]; lowE += zLow[k]; rbE += rb[k]; }
    const cnt = n - k0; thE /= cnt; topE /= cnt; lowE /= cnt; rbE /= cnt;
    // Upright free-floating statics with the final flood mass (Archimedes along the axis, lower end down).
    const mTot = j.M + mF - H / G;
    let Ls = NaN, xiBsub = NaN, cum = 0, mom = 0;
    for (let i = j.n - 1; i >= 0; i--) {
      const need1 = mTot / j.rho - cum;
      if (j.vol[i] >= need1) { const f = need1 / j.vol[i]; Ls = (j.n - 1 - i + f) * j.dxi; mom += need1 * (j.xi[i] + j.dxi / 2 - f * j.dxi / 2); cum += need1; xiBsub = mom / cum; break; }
      cum += j.vol[i]; mom += j.vol[i] * j.xi[i];
    }
    const floats = Number.isFinite(Ls);
    const xiGtot = (mF * xiF - (H / G) * xiH) / mTot; // CoG of jacket + flood (hook as negative mass at its station)
    // Upright GM = BG (B above G when ξG > ξB, ξ positive downwards) + waterplane inertia of the legs, neglected (conservative).
    const GM = floats ? xiGtot - xiBsub : NaN;
    const freeboard = floats ? j.L - Ls : 0;
    const onBottom = (j.M + mF) * G - H - j.rho * G * j.V;
    const needsBed = !floats;
    const checks = [
      check('Bottom clearance while rotating', inp.clearReq, Math.max(clearMin, 1e-6), 'm (required / achieved)', 'DNV-ST-N001 (2021) Sec. 13 (upending: bottom clearance)', 'The lower end swings down as the jacket rotates. Touching the seabed before the jacket is upright loads the legs sideways and can damage mudmats.'),
      check('Final inclination from the vertical', Math.abs(90 - thE), inp.tolDeg, 'deg', 'DNV-ST-N001 (2021) Sec. 13 (upending: final attitude); project set-down tolerance', 'The jacket must end close to vertical before set-down or pile stabbing; a residual list means the flooding or hook load does not balance the weight distribution.'),
    ];
    if (floats) checks.push(check('Upright floating stability (GM)', inp.gmReq, Math.max(GM, 1e-6), 'm (required / achieved)', 'DNV-ST-N001 (2021) Sec. 13 (stability in all upending stages)', 'With the centre of gravity above the centre of buoyancy the upright jacket is unstable and will fall back over; flood lower or add hook load.'));
    if (H > 0) checks.push(check('Hook load against crane capacity', inp.hook, Math.max(inp.hookCap, 1e-6), 't', 'DNV-ST-N001 (2021) Sec. 16 (crane capacity)', 'The prescribed hook load must be within the crane chart at the upending radius.'));
    const gov = worst(checks), ok = gov.uc <= 1;
    const [tp, a1, a2, a3, a4] = thin([tt, thD, zLow, zTop, mf], 400);
    const mFloatMax = j.rho * j.V + H / G - j.M; // flood mass at which the jacket just stops floating
    return {
      kpis: [
        kpi('thEnd', 'Final inclination of the axis from the horizontal', thE, 'deg', '90° is upright. Averaged over the last fifth of the simulation to remove residual oscillation.'),
        kpi('tUp', 'Time to reach the upright tolerance', Number.isFinite(tUp) ? tUp : inp.tEnd, 's', 'First time the axis is within the tolerance of vertical. Equal to the simulated time when the jacket never gets there.'),
        kpi('mUp', 'Flood water on board when upright is first reached', Number.isFinite(mUp) ? mUp : inp.mFlood, 't', 'The buoyancy requirement of the upending: flood water needed at the given station to bring the jacket within tolerance of vertical.'),
        kpi('clearMin', 'Minimum bottom clearance while rotating', clearMin, 'm', 'Smallest gap between the lowest point and the seabed before the jacket is upright.', inp.clearReq, 'min'),
        kpi('omMax', 'Maximum rotation rate', omMax, 'deg/s', 'Peak pitch velocity. Fast rotation means the jacket passes through an unstable attitude: slow the flooding or hold back with the hook.'),
        kpi('topEnd', 'Final elevation of the upper end', topE, 'm', 'Height of the top of the jacket above still water at the end; the freeboard available for rigging access and pile stabbing.'),
        kpi('freeboard', 'Upright free-floating freeboard (statics)', freeboard, 'm', 'Length of jacket above water when floating upright with the final flood mass, from Archimedes. Zero when the jacket cannot float and rests on the seabed.'),
        kpi('GM', 'Upright metacentric height (BG, waterplane neglected)', floats ? GM : 0, 'm', 'Centre of buoyancy above centre of gravity of the upright floating jacket. Positive is stable; reported as zero when the jacket is bottom-founded.', floats ? inp.gmReq : null, 'min'),
        kpi('onBottom', 'On-bottom weight after flooding', Math.max(onBottom, 0) / 1e3, 'kN', 'Weight less buoyancy and hook load with the jacket fully submerged: the load on the mudmats at set-down. Zero while the jacket still floats.'),
        kpi('Rbed', 'Seabed reaction at the end of the simulation', rbE / 1e3, 'kN', 'Contact force on the lower end averaged over the last fifth of the simulation.'),
        kpi('mFloatMax', 'Flood mass at which the jacket stops floating', mFloatMax / 1e3, 't', 'Reserve buoyancy plus hook load expressed as water: admitting more than this sets the jacket down.'),
        kpi('Mup', 'Upending moment of the flood water about the CoG', mF * G * Math.abs(xiF) / 1e6, 'MN·m', 'Flood weight times its lever to the centre of gravity with the jacket horizontal: the moment that starts the rotation.'),
      ],
      checks,
      plots: [
        line('Rotation and flooding history', 'Time', 's', 'Inclination (deg) · flood water (100 t)', '', [{ name: 'Axis inclination (deg)', x: tp, y: a1 }, { name: 'Flood water (×100 t)', x: tp, y: a4.map((v) => v / 100) }, { name: 'Upright', x: [0, inp.tEnd], y: [90, 90], dash: true }],
          'The jacket rotates as water is admitted. A sudden jump in angle shows the attitude at which it becomes unstable and self-rotates; a plateau short of 90° shows insufficient flooding.'),
        line('Elevation of the jacket ends', 'Time', 's', 'Elevation', 'm', [{ name: 'Lowest point', x: tp, y: a2 }, { name: 'Upper end', x: tp, y: a3 }, { name: 'Seabed', x: [0, inp.tEnd], y: [zBed, zBed], dash: true }, { name: 'Still water', x: [0, inp.tEnd], y: [0, 0], dash: true }],
          'How deep the lower end swings and how much of the top stays above water. The lowest-point curve must stay above the seabed until the jacket is upright.'),
      ],
      tables: [tbl('Upending end state', [['k', 'Quantity'], ['v', 'Value'], ['u', 'Unit']], [
        ['Total buoyant volume', j.V, 'm³'], ['Flood water admitted', inp.mFlood, 't'], ['Hook load', inp.hook, 't'], ['Final inclination', thE, 'deg'], ['Final elevation of the lower end', lowE, 'm'], ['Final elevation of the upper end', topE, 'm'],
        ['Upright submerged length (statics)', floats ? Ls : j.L, 'm'], ['Upright BG', floats ? GM : 0, 'm'], ['On-bottom weight', Math.max(onBottom, 0) / 1e3, 'kN'],
      ].map(([k, v, u]) => ({ k, v: round(v), u })), 'The configuration reached at the end of the flooding sequence, with the upright statics for cross-checking the simulated end state.')],
      notes: [needsBed ? `With ${inp.mFlood} t of flood water the jacket no longer floats: it sets down with ${fmt(onBottom / 1e3)} kN on bottom.` : `The jacket floats upright with ${fmt(freeboard)} m freeboard and BG = ${fmt(GM)} m.`, Number.isFinite(tUp) ? `Upright within ${inp.tolDeg}° after ${fmt(tUp)} s with ${fmt(mUp)} t on board.` : 'The jacket does not reach the upright tolerance within the simulated time.', `"${gov.name}" governs at ${fmt(gov.uc)}.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `The jacket upends to ${fmt(thE)}° with ${fmt(clearMin)} m bottom clearance while rotating and a peak rotation rate of ${fmt(omMax)} deg/s.`,
          `Check the same sequence with one flooded buoyancy compartment and with the flood valves stuck open (rate × 2); keep at least ${fmt(inp.clearReq)} m clearance in both. Stage the flooding so that the unstable attitude is passed under hook control.`, 'DNV-ST-N001 (2021) Sec. 13; ISO 19901-6 §15')
          : adv('engineering', 'fail', `"${gov.name}" is not met (${fmt(gov.uc)}): final inclination ${fmt(thE)}°, bottom clearance ${fmt(clearMin)} m.`,
            gov.name.startsWith('Final') ? `Admit more water lower down (the jacket stops floating above ${fmt(mFloatMax / 1e3)} t), move the flooded compartments towards the lower end, or apply hook load at the top.` : gov.name.startsWith('Bottom') ? `Upend in deeper water, hold more hook load to keep the jacket high during rotation, or flood in stages; the lower end needs ${fmt(inp.clearReq - clearMin)} m more.` : gov.name.startsWith('Upright') ? 'Lower the centre of gravity by flooding the bottom compartments first, or keep the hook attached until set-down.' : 'Reduce the prescribed hook load or use a crane with more capacity at this radius.',
            'DNV-ST-N001 (2021) Sec. 13; ISO 19901-6 §15'),
        adv('economic', 'info', `Self-upending needs ${fmt(inp.mFlood)} t of controlled flooding; the crane contributes ${fmt(inp.hook)} t.`,
          inp.hook > 0 ? 'Crane-assisted upending ties up the heavy-lift vessel for the whole sequence: compare the hook time against extra flooding compartments and valves that would let the jacket upend by itself.' : 'A self-upending design removes the heavy-lift vessel from this step; the cost is the flooding system (valves, umbilicals, closures) and its testing. Keep the option of hook assistance as a contingency.',
          'Installation method comparison, ISO 19901-6 §6; AACE RP 18R-97 for the estimate class',
          `Every hour of crane-vessel time avoided saves about ${vesselCo2Text(1 / 24, 60)}`),
      ],
      publish: { 'install.upendHookLoad': H, 'install.upendFloodMass': mF },
      quantities: { vessel_days: { upending: round(inp.tEnd / 86400 + 0.25) } },
    };
  },
  verify: [
    { name: 'Upright free-floating freeboard of a uniform buoyant column (Archimedes)', source: 'Uniform column L = 100 m, mass 6000 t, reserve 70.8333 % → V = 10 000 m³ (100 m²); with 3000 t flood water at the base the submerged length is 9.0e6/(1025·100) = 87.80 m, freeboard 12.20 m; G is 16.67 m below mid-length and B 6.10 m below it → BG = 10.57 m',
      inputs: { M: 6000, L: 100, xiG: 50, xiB: 50, kGyr: 29, reserve: 70.8333333, mFlood: 3000, xiF: 100, qFlood: 3000, th0: 1, depth: 500, Dn: 10, Aax: 100 },
      expect: [{ kpi: 'freeboard', value: 100 - 9.0e6 / (1025 * 100), tol: 0.005 }, { kpi: 'GM', value: (3000 * 50 / 9000) - (50 - 9.0e6 / (1025 * 100) / 2), tol: 0.01 }, { kpi: 'topEnd', value: 100 - 9.0e6 / (1025 * 100), tol: 0.03 }, { kpi: 'thEnd', value: 90, tol: 0.01 }] },
  ],
};
