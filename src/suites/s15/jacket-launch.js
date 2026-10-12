// Jacket launch from a trimmed barge over a rocker pin: planar rigid-body time simulation with penalty contact.
import { rk4 } from '../../engine/math/ode.js';
import { G, DEG, need, round, fmt, kpi, check, adv, line, tbl, thin, worst, vesselCo2Text, LIMIT_ROM } from './_common.js';
import { jacketInputs, buildJacket, hydro, addForce, solve3 } from './_jacket.js';

export default {
  id: 's15-jacket-launch', suite: 15, title: 'Jacket launch from a barge',
  summary: 'Simulates the jacket sliding down the skid beams of a trimmed launch barge, tipping over the rocker, separating and diving, and reports trajectory, rocker reaction, accelerations and seabed clearance.',
  covers: [3], fidelity: 'reduced-order',
  standards: ['DNV-ST-N001 (2021) Sec. 13 (launching: trim, friction, bottom clearance, reserve buoyancy)', 'API RP 2A-WSD (2014) §2.4.4 (launching forces)', 'ISO 19901-6 (2009) §15 (launching)'],
  models: ['Newton–Euler planar rigid-body dynamics (3 degrees of freedom)', 'Launch cradle Coulomb friction', 'Rocker-arm load transfer (penalty contact at the rocker pin)', 'Time-varying buoyancy by strip integration', 'Added mass and quadratic hydrodynamic damping (Morison strips)', 'Explicit Runge–Kutta time integration'],
  limits: `${LIMIT_ROM} The barge is held at its prescribed trim and draft throughout (its change of trim under the moving jacket is only estimated afterwards from the rocker reaction), the rocker is a pin without beam length, the jacket buoyancy is a linear distribution along the axis, and water-entry slamming and member-level loads are not resolved. Use it to size trim, friction and clearance before a coupled launch analysis.`,
  publishes: ['install.launchRockerReaction', 'install.launchDiveDepth'],
  searchTerms: 'jacket launch barge rocker arm trajectory simulation dive depth',
  inputs: [
    ...jacketInputs(),
    { key: 's0', label: 'Initial distance of the CoG ahead of the rocker pin', unit: 'm', type: 'number', default: 30, min: 0.5, max: 300, group: 'Initial conditions', role: 'ic', help: 'Position on the skid beams when the jacket starts to move (after cutting the sea-fastening and jacking to break out).' },
    { key: 'v0', label: 'Initial sliding velocity', unit: 'm/s', type: 'number', default: 0.2, min: 0, max: 5, group: 'Initial conditions', role: 'ic', help: 'Velocity given by the jacks or winches at break-out.' },
    { key: 'trimDeg', label: 'Barge trim (skid-beam slope)', unit: 'deg', type: 'number', default: 4, min: 0.5, max: 12, group: 'Launch barge', role: 'bc', design: true, help: 'Prescribed boundary condition: stern-down trim set by ballast, typically 2–5°.' },
    { key: 'hPin', label: 'Rocker pin height above still water', unit: 'm', type: 'number', default: 3, min: -5, max: 40, group: 'Launch barge', role: 'bc', help: 'At the stern with the launch trim and draft applied.' },
    { key: 'mu', label: 'Kinetic friction on skid beams and rocker', unit: '-', type: 'number', default: 0.04, min: 0, max: 0.4, group: 'Launch barge', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'About 0.03–0.07 for greased timber or PTFE on steel. Check both ends of the range.' },
    { key: 'muS', label: 'Static (break-out) friction', unit: '-', type: 'number', default: 0.12, min: 0, max: 0.6, group: 'Launch barge', role: 'bc', help: 'Typically 0.10–0.20 after transportation; sets the jacking force to start the launch.' },
    { key: 'jackCap', label: 'Jack / winch capacity for break-out', unit: 'kN', type: 'number', default: 12000, min: 0, max: 2e5, group: 'Launch barge', role: 'bc' },
    { key: 'rockerCap', label: 'Rocker capacity (both arms)', unit: 'MN', type: 'number', default: 140, min: 1, max: 2000, group: 'Launch barge', role: 'bc', help: 'Rated vertical load of the rocker arms from the barge documentation.' },
    { key: 'Lb', label: 'Barge length', unit: 'm', type: 'number', default: 180, min: 30, max: 400, group: 'Launch barge', role: 'geometry' },
    { key: 'Bb', label: 'Barge beam', unit: 'm', type: 'number', default: 42, min: 8, max: 120, group: 'Launch barge', role: 'geometry' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 150, min: 10, max: 3000, group: 'Site', role: 'environment', bind: 'site.depth' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'clearReq', label: 'Required seabed clearance during the dive', unit: 'm', type: 'number', default: 10, min: 0, max: 200, group: 'Criteria', role: 'criteria', help: 'DNV-ST-N001 asks for a documented bottom clearance, commonly the larger of 10 m or 10 % of the water depth. Confirm against the governing edition.' },
    { key: 'tEnd', label: 'Simulated time', unit: 's', type: 'number', default: 70, min: 10, max: 400, group: 'Numerics', role: 'control' },
    { key: 'rate', label: 'Time steps per second', unit: '1/s', type: 'int', default: 100, min: 25, max: 2000, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'rate', label: 'Time steps per second', levels: [50, 100, 200, 400], qoi: 'tTip' },
  compute(inp) {
    const j = buildJacket(inp, inp.rhoW);
    const al = inp.trimDeg * DEG, sa = Math.sin(al), ca = Math.cos(al), mu = inp.mu, hp = inp.hPin;
    // Penalty contact: static penetration of about 20 mm under the full weight, near-critical damping.
    const kc = j.M * G / 0.02, cc = 1.4 * Math.sqrt(kc * j.M);
    const kStrip = 2 * kc / j.n, cStrip = 2 * cc / j.n;
    const reg = (v) => Math.tanh(v / 0.02);
    // Returns accelerations and, through `out`, the contact reactions.
    const rhs = (y, out) => {
      const x = y[0], z = y[1], th = y[2], vx = y[3], vz = y[4], om = y[5];
      const s = Math.sin(th), c = Math.cos(th), ux = c, uz = -s, nx = s, nz = c;
      const Q = [0, -j.M * G, 0], Mm = [j.M, 0, 0, 0, j.M, 0, 0, 0, j.I];
      hydro(j, z, th, vx, vz, om, Q, Mm);
      let Rdeck = 0, Rpin = 0;
      // Skid beams: every strip still forward of the pin bears on the inclined deck line through the pin.
      for (let i = 0; i < j.n; i++) {
        const xi = j.xi[i], px = x + xi * ux, pz = z + xi * uz;
        if (px >= 0) continue;
        const pen = (hp - px * sa / ca - pz) * ca;
        if (pen <= 0) continue;
        const pvx = vx - xi * om * nx, pvz = vz - xi * om * nz;
        const N = Math.max(kStrip * pen - cStrip * (pvx * sa + pvz * ca), 0);
        const Ft = -mu * N * reg(pvx * ca - pvz * sa);
        addForce(Q, th, xi, N * sa + Ft * ca, N * ca - Ft * sa); Rdeck += N;
      }
      // Rocker pin against the jacket runner line.
      const xp = (0 - x) * ux + (hp - z) * uz;
      if (xp >= -j.La && xp <= j.Lf) {
        const pen = (0 - x) * nx + (hp - z) * nz;
        if (pen > 0) {
          const vn = vx * nx + vz * nz - xp * om;
          const N = Math.max(kc * pen - cc * vn, 0);
          const Ft = -mu * N * reg(vx * ux + vz * uz);
          addForce(Q, th, xp, N * nx + Ft * ux, N * nz + Ft * uz); Rpin = N;
        }
      }
      if (out) { out.Rdeck = Rdeck; out.Rpin = Rpin; }
      return solve3(Mm, Q);
    };
    // Start bedded into the contact springs (static penetration under the weight component normal to the deck).
    let nOn = 0; for (let i = 0; i < j.n; i++) if (-inp.s0 * ca + j.xi[i] * ca < 0) nOn++;
    const pen0 = j.M * G * ca / (kStrip * Math.max(nOn, 1) + kc);
    const y0 = [-inp.s0 * ca - pen0 * sa, hp + inp.s0 * sa - pen0 * ca, al, inp.v0 * ca, -inp.v0 * sa, 0];
    const nSteps = Math.round(inp.tEnd * inp.rate);
    const sol = rk4((t, y) => { const a = rhs(y); return [y[3], y[4], y[5], a[0], a[1], a[2]]; }, 0, y0, inp.tEnd, nSteps);
    const n = sol.t.length, tt = sol.t;
    let tRmax = 0, tTip = NaN, tSep = 0, vSep = 0, Rmax = 0, Rdmax = 0, low = Infinity, tLow = 0, aMax = 0, thMax = 0, vMax = 0;
    const lowZ = new Float64Array(n), Rp = new Float64Array(n), sp = new Float64Array(n), thD = new Float64Array(n), cx = new Float64Array(n), cz = new Float64Array(n);
    const o = {};
    for (let k = 0; k < n; k++) {
      const y = sol.y[k], th = y[2];
      need(Number.isFinite(y[0]) && Number.isFinite(th), 'The launch simulation diverged. Increase the number of time steps per second.');
      const a = rhs(y, o);
      const zLead = y[1] - j.Lf * Math.sin(th) - j.rad[j.n - 1], zTrail = y[1] + j.La * Math.sin(th) - j.rad[0];
      lowZ[k] = Math.min(zLead, zTrail); Rp[k] = o.Rpin + o.Rdeck; thD[k] = th / DEG; cx[k] = y[0]; cz[k] = y[1];
      sp[k] = Math.hypot(y[3], y[4]);
      if (lowZ[k] < low) { low = lowZ[k]; tLow = tt[k]; }
      if (o.Rpin > Rmax) { Rmax = o.Rpin; tRmax = tt[k]; }
      if (o.Rdeck > Rdmax) Rdmax = o.Rdeck;
      if (o.Rpin + o.Rdeck > 1e-3 * j.M * G) { tSep = tt[k]; vSep = sp[k]; }
      if (k > 0 && Number.isNaN(tTip) && y[0] >= 0) { const y1 = sol.y[k - 1]; tTip = tt[k - 1] + (tt[k] - tt[k - 1]) * (0 - y1[0]) / (y[0] - y1[0]); }
      if (tt[k] > 1) aMax = Math.max(aMax, Math.hypot(a[0], a[1])); thMax = Math.max(thMax, Math.abs(th) / DEG); vMax = Math.max(vMax, sp[k]);
    }
    need(Number.isFinite(tTip), `The centre of gravity does not reach the rocker pin within ${inp.tEnd} s: with ${inp.trimDeg}° trim (tan = ${fmt(Math.tan(al))}) and kinetic friction ${mu} the jacket ${Math.tan(al) > mu ? 'needs a longer simulated time' : 'does not keep sliding — increase the trim or reduce friction'}.`);
    const separated = tSep < inp.tEnd - 2 / inp.rate;
    const dive = Math.max(-low, 0), clear = inp.depth - dive;
    const Fbreak = Math.max(j.M * G * (inp.muS * ca - sa), 0);
    // Barge trim change if the peak rocker reaction acted statically at the stern: Δθ = R·(Lb/2)/(ρ·g·B·Lb³/12).
    const trimEst = Rmax * (inp.Lb / 2) / (inp.rhoW * G * inp.Bb * inp.Lb ** 3 / 12) / DEG;
    const checks = [
      check('Seabed clearance during the dive', inp.clearReq, Math.max(clear, 1e-6), 'm (required / achieved)', 'DNV-ST-N001 (2021) Sec. 13 (launch: bottom clearance)', 'The leading end plunges below its floating draft after separation. If the clearance is used up the jacket strikes the seabed and can damage legs, mudmats or pile sleeves.'),
      check('Rocker reaction', Rmax / 1e6, inp.rockerCap, 'MN', 'DNV-ST-N001 (2021) Sec. 13 (rocker arm loads); barge operating manual', 'When the centre of gravity passes the pin the whole jacket weight is carried by the rocker arms. Overload damages the arms, the barge stern or the jacket launch trusses.'),
      check('Break-out force', Fbreak / 1e3, Math.max(inp.jackCap, 1e-6), 'kN', 'DNV-ST-N001 (2021) Sec. 13 (launch: friction range, push/pull capacity)', 'Force needed to start sliding against static friction less the weight component down the slope. Without it the jacket stays on the barge.'),
    ];
    const gov = worst(checks), ok = gov.uc <= 1;
    const [tp, zp, rp, spd, thp] = thin([tt, lowZ, Rp, sp, thD], 400);
    // Stick figures of the jacket axis at eight instants.
    const snaps = [];
    for (let q = 0; q < 8; q++) { const k = Math.round(q * (n - 1) / 7), y = sol.y[k], s = Math.sin(y[2]), c = Math.cos(y[2]); snaps.push({ name: `t = ${fmt(tt[k], 3)} s`, x: [round(y[0] - j.La * c), round(y[0] + j.Lf * c)], y: [round(y[1] + j.La * s), round(y[1] - j.Lf * s)] }); }
    const aSlide = G * (sa - mu * ca);
    return {
      kpis: [
        kpi('tTip', 'Time for the CoG to reach the rocker pin', tTip, 's', 'End of the sliding phase: from here the jacket rotates on the rocker. Closed form for steady friction: s0 = v0·t + ½·g·(sin α − μ·cos α)·t².'),
        kpi('aSlide', 'Sliding acceleration on the skid beams', aSlide, 'm/s²', 'g·(sin α − μ·cos α). Negative means the jacket decelerates and may stop on the barge: more trim or less friction is needed.', 0, 'min'),
        kpi('tSep', 'Time of separation from the barge', tSep, 's', 'Last instant with contact on the skid beams or rocker. If it equals the simulated time the jacket has not left the barge.'),
        kpi('vSep', 'Velocity at separation', vSep, 'm/s', 'Speed of the centre of gravity when the trailing end leaves the rocker; with the dive it sets the stand-off distance to the barge.'),
        kpi('Rmax', 'Maximum rocker reaction', Rmax / 1e6, 'MN', 'Peak normal force on the rocker pin. It approaches the jacket weight less the buoyancy already acting on the immersed leading end.', inp.rockerCap, 'max'),
        kpi('dive', 'Maximum dive depth of the lowest point', dive, 'm', 'Deepest position below still water reached by either end of the jacket during the launch.'),
        kpi('clear', 'Minimum seabed clearance', clear, 'm', 'Water depth less the dive depth. Compare with the required clearance.', inp.clearReq, 'min'),
        kpi('aMax', 'Maximum acceleration of the centre of gravity', aMax, 'm/s²', 'Peak rigid-body acceleration during tipping and water entry; multiply member masses by it for inertia loads on appurtenances.'),
        kpi('thMax', 'Maximum pitch angle', thMax, 'deg', 'Largest inclination of the jacket axis from the horizontal during tipping and the dive.'),
        kpi('Fbreak', 'Break-out force', Fbreak / 1e3, 'kN', 'Push or pull needed to overcome static friction at the launch trim. Zero when the slope alone starts the jacket.', inp.jackCap, 'max'),
        kpi('trimEst', 'Barge trim change under the peak rocker load (static estimate)', trimEst, 'deg', 'Indicative only: the barge is held at fixed trim in the simulation. A value comparable with the launch trim means barge response must be modelled.'),
      ],
      checks,
      plots: [
        line('Launch trajectory', 'Horizontal position astern of the rocker pin', 'm', 'Elevation', 'm', [...snaps, { name: 'CoG path', x: Array.from(thin([cx], 300)[0]), y: Array.from(thin([cz], 300)[0]), dash: true }, { name: 'Still water', x: [-j.L, 2 * j.L], y: [0, 0], dash: true }],
          'Stick figures of the jacket axis at eight instants and the path of the centre of gravity. Look for the deepest excursion of the leading end and how far the jacket travels astern.'),
        line('Elevation of the lowest point', 'Time', 's', 'Elevation', 'm', [{ name: 'Lowest point of the jacket', x: tp, y: zp }, { name: 'Seabed + required clearance', x: [0, inp.tEnd], y: [-inp.depth + inp.clearReq, -inp.depth + inp.clearReq], dash: true }],
          'The dive of the jacket after separation. The curve must stay above the dashed clearance line.'),
        line('Skid and rocker reaction', 'Time', 's', 'Normal reaction', 'MN', [{ name: 'Total contact reaction', x: tp, y: rp.map((v) => v / 1e6) }, { name: 'Rocker capacity', x: [0, inp.tEnd], y: [inp.rockerCap, inp.rockerCap], dash: true }],
          'Load transfer from the skid beams to the rocker and its release as buoyancy picks up the leading end. The peak is the rocker design load.'),
        line('Speed and pitch', 'Time', 's', 'Speed (m/s) · pitch (deg)', '', [{ name: 'CoG speed (m/s)', x: tp, y: spd }, { name: 'Pitch angle (deg)', x: tp, y: thp }],
          'Sliding speed builds up until tipping; the pitch angle peaks during the dive and then decays towards the floating attitude.'),
      ],
      tables: [tbl('Launch stages', [['k', 'Stage'], ['t', 'Time', 's'], ['v', 'Value'], ['u', 'Unit']], [
        { k: 'CoG at the rocker pin (tipping starts)', t: round(tTip), v: round(aSlide * tTip + inp.v0), u: 'm/s sliding speed (closed form)' },
        { k: 'Peak rocker reaction', t: round(tRmax), v: round(Rmax / 1e6), u: 'MN' },
        { k: 'Separation', t: round(tSep), v: round(vSep), u: 'm/s' },
        { k: 'Deepest point', t: round(tLow), v: round(dive), u: 'm below still water' },
      ], 'The sequence of the launch with the time and governing value of each stage; use it to brief the marine crew and to set the camera and survey plan.')],
      notes: [separated ? `The jacket separates after ${fmt(tSep)} s at ${fmt(vSep)} m/s.` : 'The jacket has not separated within the simulated time — extend it.', `Reserve buoyancy ${inp.reserve} %; total buoyant volume ${fmt(j.V)} m³.`, `"${gov.name}" governs at ${fmt(gov.uc)}.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `The launch completes with ${fmt(clear)} m seabed clearance (dive ${fmt(dive)} m) and a peak rocker reaction of ${fmt(Rmax / 1e6)} MN (${fmt(100 * Rmax / (j.M * G))} % of the jacket weight).`,
          `Repeat with friction ${fmt(mu * 0.6)} and ${fmt(mu * 1.6)}, with one buoyancy compartment flooded, and with the trim ±0.5°; then confirm with a coupled barge–jacket launch analysis before fixing the ballast plan.`, 'DNV-ST-N001 (2021) Sec. 13; ISO 19901-6 §15')
          : adv('engineering', 'fail', `"${gov.name}" is not met (${fmt(gov.uc)}): clearance ${fmt(clear)} m against ${inp.clearReq} m, rocker ${fmt(Rmax / 1e6)} MN against ${inp.rockerCap} MN, break-out ${fmt(Fbreak / 1e3)} kN.`,
            gov.name.startsWith('Seabed') ? `Add buoyancy (reserve above ${inp.reserve} %), reduce the launch trim to lower the entry speed, or launch in deeper water: about ${fmt(inp.clearReq + dive)} m is needed for this configuration.` : gov.name.startsWith('Rocker') ? 'Lower the rocker pin towards the water or increase the trim so that buoyancy relieves the rocker earlier, or use a barge with stronger rocker arms.' : `Provide at least ${fmt(Fbreak / 1e3)} kN of jacking capacity or increase the trim to ${fmt(Math.atan(inp.muS) / DEG)}° to self-start.`,
            'DNV-ST-N001 (2021) Sec. 13; ISO 19901-6 §15'),
        adv('economic', 'info', `Launch avoids a heavy-lift vessel for a ${fmt(inp.M)} t jacket but needs ${fmt(j.V - j.M / inp.rhoW)} m³ of reserve buoyancy and a launch barge with ${fmt(Rmax / 1e6)} MN rocker capacity.`,
          'Compare the steel of launch trusses and auxiliary buoyancy tanks (often 5–8 % of jacket weight, removed after upending) against the day rate and availability of a crane vessel able to lift the jacket.',
          'Concept selection by installation method, ISO 19901-6 §6; AACE RP 18R-97 for the estimate class of the comparison',
          `Launch steel at 6 % of the jacket weight is about ${fmt(0.06 * inp.M)} t ≈ ${fmt(0.06 * inp.M * 1.9)} t CO₂e (screening factor 1.9 t CO₂e/t); a two-day heavy-lift alternative burns about ${vesselCo2Text(2, 60)}`),
      ],
      publish: { 'install.launchRockerReaction': Rmax, 'install.launchDiveDepth': dive },
      quantities: { vessel_days: { launchBarge: 1 } },
    };
  },
  verify: [
    { name: 'Sliding phase on a dry incline: s0 = v0·t + ½·g·(sin α − μ·cos α)·t²', source: 'Particle on a rough incline (Newton): α = 4°, μ = 0.04, s0 = 30 m, v0 = 0.2 m/s, pin 30 m above water so that no buoyancy acts → a = 9.80665·(0.069756 − 0.04·0.997564) = 0.29276 m/s², t = (−0.2 + √(0.04 + 2·0.29276·30))/0.29276 = 13.65 s',
      inputs: { hPin: 30, tEnd: 20 }, expect: [{ kpi: 'tTip', value: (-0.2 + Math.sqrt(0.04 + 2 * 9.80665 * (Math.sin(4 * Math.PI / 180) - 0.04 * Math.cos(4 * Math.PI / 180)) * 30)) / (9.80665 * (Math.sin(4 * Math.PI / 180) - 0.04 * Math.cos(4 * Math.PI / 180))), tol: 0.01 }] },
  ],
};
