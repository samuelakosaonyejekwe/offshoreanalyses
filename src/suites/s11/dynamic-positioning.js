// Dynamic positioning and station keeping: mean wind, current and wave-drift force on the unit, static thrust
// allocation intact and after the worst single thruster failure, the capability polar, PID or LQR controller gains
// and a closed-loop three-degree-of-freedom station-keeping simulation with a wave filter and rate-limited thrusters.
import { allocateThrust, dpPidGains, dpLqrGains, simulateDp } from '../../engine/hydro/dp.js';
import { meanDriftForce } from '../../engine/hydro/drift.js';
import { hullInputs, hullFrom, seaInputs, seaFrom, spectrum, omegaGrid, windCurrentInputs, windCurrentForce, driftCoefficient,
  DEG, CO2_MGO, LIMIT_CFD, need, run, round, rnd, fmt, adv, check, sev, worst, thin } from './_common.js';

const thr = (name, x, y, type = 'azimuth', TmaxkN = 550) => ({ name, x, y, type, TmaxkN });

export default {
  id: 's11-dynamic-positioning',
  suite: 11,
  title: 'Dynamic positioning and station keeping',
  summary: 'Whether the thrusters can hold the unit on station in a given wind, current and sea — intact and after the worst single failure — with the capability polar, the controller gains and the simulated position footprint.',
  covers: [16],
  fidelity: 'reduced-order',
  standards: ['IMO MSC.1/Circ.1580 (2017) (DP equipment classes)', 'DNV-ST-0111 (2021) (assessment of station-keeping capability of DP vessels)', 'IMCA M 140 Rev. 1 (2000) (DP capability plots)', 'DNV-RU-SHIP Pt.6 Ch.3 Sec.1 (2022) (DYNPOS notations)', 'API RP 2SK (2005) §7 (thruster-assisted and DP station keeping)'],
  models: ['Mean wind and current force on projected areas with user coefficients', 'Mean wave-drift force from a screening drift coefficient and the wave spectrum', 'Thrust allocation by bounded (redistributed) pseudo-inverse', 'Capability polar by bisection on wind speed with the bounded allocation, intact and worst single failure', 'PID gains by pole placement and LQR gains from the algebraic Riccati equation', 'Closed-loop low-speed 3-DOF simulation with steady-state Kalman wave filter, bias estimation and thruster rate limits'],
  limits: `DP is a three-degree-of-freedom (surge, sway, yaw) linear low-speed model: no Coriolis or nonlinear damping, no thrust losses, thruster–thruster or thruster–hull interaction, no power-plant or blackout logic. Wind, current and wave-drift coefficients are user inputs; the wave drift is a short-wave screening value, and in the capability polar wind and wave forces scale with the square of wind speed while current stays constant. Model-predictive and extended-Kalman schemes are not implemented. ${LIMIT_CFD}`,
  searchTerms: 'dynamic positioning capability plot thrust allocation station keeping',
  inputs: [
    ...hullInputs({ type: 'semi', only: ['semi', 'ship', 'barge', 'spar'], withMesh: false }),
    ...seaInputs({ Hs: 5, Tp: 10, ret: '1', heading: 90, withDuration: false }),
    ...windCurrentInputs({ Uw: 20, Uc: 0.8, ret: '1' }),
    { key: 'thrusters', label: 'Thrusters', type: 'table', group: 'Thrusters', role: 'geometry',
      columns: [{ key: 'name', label: 'Thruster', type: 'text' }, { key: 'x', label: 'x (fwd of midship)', unit: 'm', type: 'number' }, { key: 'y', label: 'y (to port)', unit: 'm', type: 'number' }, { key: 'type', label: 'Type (azimuth / tunnel)', type: 'text' }, { key: 'TmaxkN', label: 'Maximum thrust', unit: 'kN', type: 'number' }],
      default: [thr('T1 fwd port', 50, 35), thr('T2 fwd port', 38, 35), thr('T3 fwd stbd', 50, -35), thr('T4 fwd stbd', 38, -35), thr('T5 aft port', -38, 35), thr('T6 aft port', -50, 35), thr('T7 aft stbd', -38, -35), thr('T8 aft stbd', -50, -35)],
      help: 'Position and bollard thrust of each thruster. Azimuth thrusters turn freely; tunnel thrusters act athwartships in both directions. Enter thrust net of losses.' },
    { key: 'kNperKW', label: 'Thrust per unit power at full load', unit: 'kN/kW', type: 'number', default: 0.17, min: 0.05, max: 0.35, group: 'Thrusters', role: 'operation', help: 'Bollard-pull efficiency of the thrusters (0.15–0.20 for ducted azimuth units); power is taken to vary with thrust^1.5.' },
    { key: 'thrustRate', label: 'Thrust slew rate', unit: '1/s', type: 'number', default: 0.1, min: 0.01, max: 1, group: 'Thrusters', role: 'bc', help: 'Fraction of maximum thrust the thruster can gain or shed per second.' },
    { key: 'aziRate', label: 'Azimuth slew rate', unit: 'deg/s', type: 'number', default: 12, min: 1, max: 60, group: 'Thrusters', role: 'bc', help: 'Turning rate of the azimuth thrusters (typically 12°/s, two revolutions per minute).' },
    { key: 'kYaw', label: 'Yaw gyradius / overall length', unit: '-', type: 'number', default: 0.35, min: 0.15, max: 0.8, group: 'Vessel dynamics', role: 'ic', help: 'Yaw radius of gyration as a fraction of the overall length.' },
    { key: 'amSurge', label: 'Surge added mass / mass', unit: '-', type: 'number', default: 0.15, min: 0, max: 2, group: 'Vessel dynamics', role: 'bc', help: 'Low-frequency added mass in surge (0.05–0.15 ships, 0.2–0.6 column-stabilised units).' },
    { key: 'amSway', label: 'Sway added mass / mass', unit: '-', type: 'number', default: 0.6, min: 0, max: 2, group: 'Vessel dynamics', role: 'bc', help: 'Low-frequency added mass in sway (0.5–1.0 for ships beam-on).' },
    { key: 'amYaw', label: 'Yaw added inertia / inertia', unit: '-', type: 'number', default: 0.5, min: 0, max: 2, group: 'Vessel dynamics', role: 'bc', help: 'Low-frequency added inertia in yaw.' },
    { key: 'Tdamp', label: 'Low-frequency damping time constant', unit: 's', type: 'number', default: 150, min: 20, max: 2000, group: 'Vessel dynamics', role: 'bc', help: 'Linear damping = (mass + added mass) / time constant in each degree of freedom (100–250 s is typical).' },
    { key: 'controller', label: 'Controller', type: 'select', options: [{ value: 'pid', label: 'PID by pole placement' }, { value: 'lqr', label: 'LQR (Riccati) with integral action' }], default: 'pid', group: 'Control system', role: 'control', help: 'Both act on the wave-filtered low-frequency position and velocity estimates.' },
    { key: 'Tctrl', label: 'Closed-loop natural period (PID)', unit: 's', type: 'number', default: 100, min: 20, max: 400, group: 'Control system', role: 'control', design: true, showIf: { key: 'controller', in: ['pid'] }, help: 'Bandwidth of the position loop: 60–150 s keeps the controller away from wave frequencies.' },
    { key: 'zetaCtrl', label: 'Closed-loop damping ratio (PID)', unit: '-', type: 'number', default: 0.7, min: 0.3, max: 1.5, group: 'Control system', role: 'control', showIf: { key: 'controller', in: ['pid'] }, help: '0.7–1.0 gives a well-damped response.' },
    { key: 'lqrPos', label: 'LQR position error scale', unit: 'm', type: 'number', default: 2, min: 0.1, max: 50, group: 'Control system', role: 'control', showIf: { key: 'controller', in: ['lqr'] }, help: 'Bryson weighting: the position error considered as costly as half the total thrust. Smaller = stiffer control.' },
    { key: 'lqrPsi', label: 'LQR heading error scale', unit: 'deg', type: 'number', default: 2, min: 0.1, max: 30, group: 'Control system', role: 'control', showIf: { key: 'controller', in: ['lqr'] }, help: 'Bryson weighting for the heading error.' },
    { key: 'x0', label: 'Initial offset north (x)', unit: 'm', type: 'number', default: 0, min: -100, max: 100, group: 'Simulation', role: 'ic', help: 'Initial position error along x; thrusters start at zero thrust.' },
    { key: 'y0', label: 'Initial offset east (y)', unit: 'm', type: 'number', default: 0, min: -100, max: 100, group: 'Simulation', role: 'ic', help: 'Initial position error along y.' },
    { key: 'psi0', label: 'Initial heading error', unit: 'deg', type: 'number', default: 0, min: -45, max: 45, group: 'Simulation', role: 'ic', help: 'Initial deviation from the heading setpoint.' },
    { key: 'slowFrac', label: 'Slowly varying force (fraction of mean)', unit: '-', type: 'number', default: 0.3, min: 0, max: 1, group: 'Simulation', role: 'bc', help: 'Amplitude of a sinusoidal gust / slow-drift variation added to the mean environmental force.' },
    { key: 'slowT', label: 'Period of the slowly varying force', unit: 's', type: 'number', default: 80, min: 10, max: 600, group: 'Simulation', role: 'bc', help: 'Typical wind-gust and wave-group periods are 30–200 s.' },
    { key: 'wfSig', label: 'Wave-frequency motion in the position measurement (σ)', unit: 'm', type: 'number', default: 0.5, min: 0, max: 5, group: 'Simulation', role: 'bc', help: 'First-order surge/sway motion the wave filter must remove so that the thrusters do not chase the waves.' },
    { key: 'noise', label: 'Position sensor noise (σ)', unit: 'm', type: 'number', default: 0.2, min: 0, max: 5, group: 'Simulation', role: 'bc', help: 'Standard deviation of the position reference noise (DGNSS ≈ 0.1–0.5 m).' },
    { key: 'duration', label: 'Simulated duration', unit: 's', type: 'number', default: 600, min: 60, max: 7200, group: 'Simulation', role: 'control', help: 'The environment is applied as a step at t = 0; statistics are taken over the last 75 %.' },
    { key: 'dt', label: 'Time step', unit: 's', type: 'number', default: 1, min: 0.05, max: 2, group: 'Simulation', role: 'control', help: 'Controller and plant step (DP controllers typically run at 1 Hz).' },
    { key: 'utilIntact', label: 'Allowable thrust utilisation, intact', unit: '-', type: 'number', default: 0.8, min: 0.3, max: 1, group: 'Acceptance criteria', role: 'criteria', help: '80 % leaves a 20 % allowance for dynamics (IMCA M 140 practice).' },
    { key: 'utilFail', label: 'Allowable thrust utilisation after single failure', unit: '-', type: 'number', default: 1.0, min: 0.3, max: 1, group: 'Acceptance criteria', role: 'criteria', help: 'After the worst single failure the remaining thrusters may be used fully (DP class 2 and 3).' },
    { key: 'footMax', label: 'Allowable position excursion', unit: 'm', type: 'number', default: 5, min: 0.2, max: 200, group: 'Acceptance criteria', role: 'criteria', help: 'Watch-circle radius set by the riser, gangway or drilling operation.' },
    { key: 'psiMax', label: 'Allowable heading deviation', unit: 'deg', type: 'number', default: 3, min: 0.2, max: 45, group: 'Acceptance criteria', role: 'criteria', help: 'Largest acceptable heading error.' },
  ],
  meshTime: { key: 'dt', label: 'Time step', levels: [2, 1, 0.5], qoi: 'maxRadius', isSize: true },
  compute(inp) {
    const hull = hullFrom(inp); const rho = hull.rho;
    need(hull.type !== 'tlp', 'A tension-leg platform is not dynamically positioned: choose another hull type.');
    const rows = (inp.thrusters || []).filter((t) => t.TmaxkN > 0);
    need(rows.length >= 2, 'At least two thrusters with a positive maximum thrust are needed to control surge, sway and yaw.');
    const thrusters = rows.map((t) => {
      const type = String(t.type || 'azimuth').trim().toLowerCase();
      need(type === 'azimuth' || type === 'tunnel', `Thruster "${t.name}": type must be "azimuth" or "tunnel".`);
      return { name: t.name, x: t.x, y: t.y, type, Tmax: t.TmaxkN * 1e3 };
    });
    const span = Math.max(...thrusters.map((t) => Math.hypot(t.x, t.y)));
    need(span > 0 && new Set(thrusters.map((t) => `${t.x},${t.y}`)).size >= 2, 'All thrusters are at one point: yaw cannot be controlled. Spread the thrusters.');
    const Ttot = thrusters.reduce((a, t) => a + t.Tmax, 0);
    // --- mean environmental force (collinear wind, waves and current travelling towards headingDeg) --------
    const dir = inp.headingDeg; const om = omegaGrid(120); const sea = seaFrom(inp); const S = spectrum(om, sea);
    const Dw = driftCoefficient(inp, hull, om, dir);
    const drift = inp.Hs > 0 ? meanDriftForce({ omega: om, S, D: Dw }) : 0;
    const wc0 = windCurrentForce(inp, hull, dir, inp.Uw, inp.Uc);
    const e0 = { Fx: wc0.wind.Fx + wc0.cur.Fx + drift * Math.cos(dir * DEG), Fy: wc0.wind.Fy + wc0.cur.Fy + drift * Math.sin(dir * DEG), Mz: wc0.wind.Mz + wc0.cur.Mz, wind: wc0.wind, cur: wc0.cur, drift };
    const Fenv = Math.hypot(e0.Fx, e0.Fy);
    const demand = { Fx: -e0.Fx, Fy: -e0.Fy, Mz: -e0.Mz };
    // --- static allocation, intact and each single failure ---------------------------------------------------
    const alloc = run(() => allocateThrust({ thrusters, demand }), 'Thrust allocation');
    const fails = thrusters.map((_, i) => ({ i, a: allocateThrust({ thrusters, demand, disabled: [i] }) }));
    const wf = fails.reduce((a, f) => (f.a.demandFraction < a.a.demandFraction - 1e-9 || (Math.abs(f.a.demandFraction - a.a.demandFraction) <= 1e-9 && f.a.maxUtilisation > a.a.maxUtilisation) ? f : a), fails[0]);
    const uInt = alloc.feasible ? alloc.maxUtilisation : 1 / Math.max(alloc.demandFraction, 1e-3);
    const uFail = wf.a.feasible ? wf.a.maxUtilisation : 1 / Math.max(wf.a.demandFraction, 1e-3);
    // --- capability polar: limiting wind speed with wind and wave drift ∝ V², current fixed --------------------
    const dirs = Array.from({ length: 24 }, (_, i) => 15 * i);
    const driftDir = new Map(); // drift force per direction at the design sea state, cached for the bisection
    const envPolar = (V, d) => {
      if (!driftDir.has(d)) driftDir.set(d, inp.Hs > 0 ? meanDriftForce({ omega: om, S, D: driftCoefficient(inp, hull, om, d) }) : 0);
      const wc = windCurrentForce(inp, hull, d, V, inp.Uc); const dr = driftDir.get(d) * (inp.Uw > 0 ? (V / inp.Uw) ** 2 : 0);
      return { Fx: wc.wind.Fx + wc.cur.Fx + dr * Math.cos(d * DEG), Fy: wc.wind.Fy + wc.cur.Fy + dr * Math.sin(d * DEG), Mz: wc.wind.Mz + wc.cur.Mz };
    };
    const vTop = 60;
    // Limiting wind speed (IMCA M 140 static balance). The force resultant can never exceed the summed available thrust,
    // which gives an upper bound without calling the allocator. When the bounded allocation is infeasible there, it
    // reports the deliverable fraction of the demand; the wind speed whose resultant equals that fraction is the limit
    // (the demand direction hardly changes with wind speed), confirmed by a feasible allocation 0.1 % below it.
    const limitV = (d, disabled, util) => {
      const th = thrusters.map((t) => ({ ...t, Tmax: t.Tmax * util }));
      const avail = th.reduce((a, t, i) => a + (disabled.includes(i) ? 0 : t.Tmax), 0);
      const dem = (V) => { const e = envPolar(V, d); return { Fx: -e.Fx, Fy: -e.Fy, Mz: -e.Mz }; };
      const mag = (V) => { const q = dem(V); return Math.hypot(q.Fx, q.Fy); };
      const vFor = (F, top) => { if (mag(top) <= F) return top; let a = 0; let b = top; while (b - a > 0.005) { const mid = 0.5 * (a + b); if (mag(mid) <= F) a = mid; else b = mid; } return a; };
      const alloc = (V) => allocateThrust({ thrusters: th, demand: dem(V), disabled });
      if (!(avail > 0) || mag(0) > avail) return 0;
      let hi = vFor(avail, vTop);
      for (let k = 0; k < 6; k++) {
        const r = alloc(hi);
        if (r.feasible) return hi;
        const target = r.demandFraction * mag(hi);
        if (!(target > mag(0))) break;
        hi = vFor(target, hi) * 0.999;
      }
      // fallback: plain bisection (current alone close to the thrust limit, or a strongly turning demand)
      if (!alloc(0).feasible) return 0;
      let lo = 0;
      while (hi - lo > 0.1) { const mid = 0.5 * (lo + hi); if (alloc(mid).feasible) lo = mid; else hi = mid; }
      return lo;
    };
    const capI = { intact: dirs.map((d) => limitV(d, [], inp.utilIntact)) };
    // the single-failure polar needs one solve per thruster and direction: it is drawn on every second direction
    const dirsF = dirs.filter((_, i) => i % 2 === 0);
    const capF = { failure: dirsF.map((d) => Math.min(...thrusters.map((_, i) => limitV(d, [i], inp.utilFail)))) };
    const vInt = Math.min(...capI.intact); const vFail = Math.min(...capF.failure);
    const vHead = capI.intact[0];
    // --- controller ------------------------------------------------------------------------------------------
    const m = hull.mass; const Izz = m * (inp.kYaw * hull.Lc) ** 2;
    const Md = [m * (1 + inp.amSurge), m * (1 + inp.amSway), Izz * (1 + inp.amYaw)];
    const M = [[Md[0], 0, 0], [0, Md[1], 0], [0, 0, Md[2]]]; const D = M.map((r) => r.map((v) => v / inp.Tdamp));
    let gains;
    if (inp.controller === 'pid') {
      const wn = 2 * Math.PI / inp.Tctrl;
      gains = dpPidGains({ M, D, omegaN: [wn, wn, wn], zeta: [inp.zetaCtrl, inp.zetaCtrl, inp.zetaCtrl] });
      need(gains.Kd.every((v) => v > 0), 'The requested closed-loop damping is below the natural damping of the vessel: shorten the closed-loop period or raise the damping ratio.');
    } else {
      // Bryson's rule: state weights 1/(allowable error)², control weights 1/(half the available thrust or moment)²
      const fS = 0.5 * Ttot; const mS = fS * span;
      const Q = Array.from({ length: 6 }, (_, i) => Array.from({ length: 6 }, (_, j) => (i === j && i < 3 ? 1 / (i < 2 ? inp.lqrPos : inp.lqrPsi * DEG) ** 2 : 0)));
      const R = [[1 / fS ** 2, 0, 0], [0, 1 / fS ** 2, 0], [0, 0, 1 / mS ** 2]];
      const g = run(() => dpLqrGains({ M, D, Q, R }), 'LQR gains');
      const Kp = [0, 1, 2].map((i) => g.Kp[i][i]); const Kd = [0, 1, 2].map((i) => g.Kd[i][i]);
      gains = { Kp, Kd, Ki: Kp.map((k, i) => 0.1 * Math.sqrt(k / Md[i]) * k) };
    }
    const Tcl = gains.Kp.map((k, i) => 2 * Math.PI * Math.sqrt(Md[i] / k));
    // --- closed-loop simulation ----------------------------------------------------------------------------
    need(inp.duration / inp.dt <= 40000, `The simulation needs ${fmt(inp.duration / inp.dt)} steps: raise the time step or shorten the duration (limit 40 000).`);
    const sim = run(() => simulateDp({ M, D, thrusters, gains, initial: { x: inp.x0, y: inp.y0, psiDeg: inp.psi0 },
      envForce: (t) => { const f = 1 + inp.slowFrac * Math.sin(2 * Math.PI * t / inp.slowT); return { Fx: e0.Fx * f, Fy: e0.Fy * f, Mz: e0.Mz * f }; },
      waveMotion: inp.wfSig > 0 ? { sigma: [inp.wfSig, inp.wfSig, inp.wfSig / Math.max(hull.Lc, 1)], period: inp.Tp } : undefined,
      noiseStd: [inp.noise, inp.noise, inp.noise / Math.max(hull.Lc, 1)], thrustRate: inp.thrustRate, azimuthRateDeg: inp.aziRate, dt: inp.dt, duration: inp.duration, seed: 11 }), 'DP simulation');
    const fp = sim.footprint;
    let peak = 0; for (let i = 0; i < sim.t.length; i++) peak = Math.max(peak, Math.hypot(sim.x[i], sim.y[i]));
    const uSim = Math.max(...sim.utilisation.map((u) => u.max));
    // mean power with P ∝ T^1.5 and the static allocation as the mean operating point
    const pFull = (t) => t.Tmax / 1e3 / inp.kNperKW;
    const power = alloc.thrusters.reduce((a, t, i) => a + pFull(thrusters[i]) * Math.min(t.utilisation, 1) ** 1.5, 0);
    const fuelDay = power * 24 * 0.21 / 1e3;
    const checks = [
      check('Thrust utilisation, intact (static balance)', uInt, inp.utilIntact, '-', 'DNV-ST-0111 (2021) Sec.4 (DP capability level 1); IMCA M 140 Rev. 1 (80 % of thrust for mean loads)', 'Guards the margin for dynamics: if the mean environment already uses most of the thrust, gusts and wave groups push the unit off station.'),
      check('Thrust utilisation after the worst single failure', uFail, inp.utilFail, '-', 'IMO MSC.1/Circ.1580 (2017) §2.2 (class 2/3: position kept after a single failure); DNV-ST-0111 (2021)', 'Guards the redundancy concept: the unit must hold position with the worst thruster (or its switchboard group) lost.'),
      check('Position excursion in the simulation', Math.max(fp.maxRadius, 0), inp.footMax, 'm', 'Watch-circle limit of the operation (project criterion); API RP 2SK (2005) §7', 'Guards the riser, gangway or well connection, which are released or damaged when the unit leaves the watch circle.'),
      check('Heading deviation in the simulation', fp.maxPsiDevDeg, inp.psiMax, 'deg', 'Operational heading limit (project criterion)', 'Guards heading-sensitive operations (offloading hose, gangway, weathervaning) and keeps the environmental force from growing as the unit turns beam-on.'),
    ];
    const gov = worst(checks);
    const [ts, xs, ys, ps] = thin([sim.t, sim.x, sim.y, sim.psiDeg], 600);
    const need1 = Fenv / 1e3 / inp.utilIntact;
    return {
      kpis: [
        { key: 'Fenv', label: 'Mean environmental force', value: round(Fenv / 1e3, 5), unit: 'kN', explain: 'Resultant of the mean wind, current and wave-drift forces the thrusters must balance.' },
        { key: 'Fwind', label: 'Wind force', value: round(Math.hypot(e0.wind.Fx, e0.wind.Fy) / 1e3, 5), unit: 'kN', explain: 'Mean wind force on the projected windage with the user coefficient.' },
        { key: 'Fcur', label: 'Current force', value: round(Math.hypot(e0.cur.Fx, e0.cur.Fy) / 1e3, 5), unit: 'kN', explain: 'Current drag on the projected underwater area with the user coefficient.' },
        { key: 'Fdrift', label: 'Mean wave-drift force', value: round(e0.drift / 1e3, 5), unit: 'kN', explain: 'Second-order mean wave force from the screening drift coefficient: the least certain of the three components.' },
        { key: 'Menv', label: 'Environmental yaw moment', value: round(e0.Mz / 1e6, 5), unit: 'MN·m', explain: 'Yaw moment of the environment about midship (zero for symmetric areas without lever arms).' },
        { key: 'thrustTotal', label: 'Installed thrust', value: round(Ttot / 1e3, 5), unit: 'kN', explain: 'Sum of the maximum thrust of all thrusters.' },
        { key: 'utilMax', label: 'Highest thruster utilisation, intact', value: round(uInt, 5), unit: '-', limit: inp.utilIntact, sense: 'max', explain: 'Largest thrust / maximum thrust in the static allocation. Above 1 the demand cannot be met and the value is the factor by which thrust falls short.' },
        { key: 'utilFail', label: 'Highest utilisation after the worst single failure', value: round(uFail, 5), unit: '-', limit: inp.utilFail, sense: 'max', explain: `Same after losing the most critical thruster (${thrusters[wf.i].name}).` },
        { key: 'vLimIntact', label: 'Limiting wind speed, intact (worst direction)', value: round(vInt, 4), unit: 'm/s', explain: `Lowest value of the capability polar over all directions with ${fmt(100 * inp.utilIntact)} % of thrust; ${vTop} m/s means the limit was not reached.` },
        { key: 'vLimFail', label: 'Limiting wind speed after single failure (worst direction)', value: round(vFail, 4), unit: 'm/s', explain: 'Lowest value of the worst-single-failure polar: the weather in which the unit keeps position with one thruster lost.' },
        { key: 'vLimHead', label: 'Limiting wind speed, environment along +x', value: round(vHead, 4), unit: 'm/s', explain: 'Intact capability with wind, waves and current travelling towards the bow direction.' },
        { key: 'KpSurge', label: 'Surge position gain Kp', value: round(gains.Kp[0] / 1e3, 5), unit: 'kN/m', explain: 'Restoring thrust per metre of surge error: the DP system acts like a spring of this stiffness.' },
        { key: 'KdSurge', label: 'Surge velocity gain Kd', value: round(gains.Kd[0] / 1e3, 5), unit: 'kN·s/m', explain: 'Damping thrust per m/s of low-frequency surge velocity.' },
        { key: 'Tclosed', label: 'Closed-loop surge period', value: round(Tcl[0], 5), unit: 's', explain: 'Natural period of the controlled surge motion; it should stay well above the wave periods (> 5·Tp is comfortable).' },
        { key: 'maxRadius', label: 'Largest excursion after settling', value: round(fp.maxRadius, 5), unit: 'm', limit: inp.footMax, sense: 'max', explain: 'Largest distance from the setpoint in the statistics window (last 75 % of the run).' },
        { key: 'peakTransient', label: 'Peak excursion including the start-up transient', value: round(peak, 5), unit: 'm', explain: 'Largest distance over the whole run, including the response to the environment applied as a step with thrusters at zero.' },
        { key: 'r95', label: '95 % footprint radius', value: round(fp.r95, 5), unit: 'm', explain: 'Radius containing 95 % of the simulated positions after settling.' },
        { key: 'headingDev', label: 'Largest heading deviation', value: round(fp.maxPsiDevDeg, 5), unit: 'deg', limit: inp.psiMax, sense: 'max', explain: 'Largest heading error in the statistics window.' },
        { key: 'saturated', label: 'Share of time with saturated allocation', value: round(sim.saturatedFraction, 5), unit: '-', explain: 'Fraction of steps in which the commanded force could not be delivered. Anything above a few percent means position is being lost.' },
        { key: 'power', label: 'Mean thruster power', value: round(power, 5), unit: 'kW', explain: 'Shaft power at the mean operating point, with power proportional to thrust^1.5.' },
      ],
      checks,
      plots: [
        { type: 'polar', title: 'DP capability polar (limiting wind speed)',
          series: [{ name: 'Intact', thetaDeg: dirs, r: rnd(capI.intact) }, { name: 'Worst single failure', thetaDeg: dirsF, r: rnd(capF.failure) }, { name: 'Design wind', thetaDeg: dirs, r: dirs.map(() => inp.Uw) }],
          explain: 'Limiting wind speed (m/s) against the direction the environment travels towards. The design-wind circle must lie inside both envelopes; the dents show the beam-on directions where windage and drift are largest.' },
        { type: 'scatter', title: 'Position footprint', x: { label: 'x (north)', unit: 'm' }, y: { label: 'y (east)', unit: 'm' },
          series: [{ name: 'Track', x: rnd(xs), y: rnd(ys) }, { name: 'Watch circle', x: rnd(Array.from({ length: 61 }, (_, i) => inp.footMax * Math.cos(i * 6 * DEG))), y: rnd(Array.from({ length: 61 }, (_, i) => inp.footMax * Math.sin(i * 6 * DEG))) }],
          explain: 'Low-frequency track of the unit around the setpoint with the allowable watch circle. The initial excursion is the start-up transient; the cluster shows steady station keeping.' },
        { type: 'line', title: 'Position and heading history', x: { label: 'Time', unit: 's' }, y: { label: 'Offset (m) / heading (deg)', unit: '-' },
          series: [{ name: 'x', x: rnd(ts), y: rnd(xs) }, { name: 'y', x: rnd(ts), y: rnd(ys) }, { name: 'Heading', x: rnd(ts), y: rnd(ps) }],
          explain: 'Response to the environment applied at t = 0. A good tuning returns to the setpoint within two or three closed-loop periods without sustained oscillation; the remaining ripple follows the slowly varying force.' },
        { type: 'bar', title: 'Thruster utilisation', categories: thrusters.map((t) => t.name), y: { label: 'Thrust / maximum', unit: '-' },
          series: [{ name: 'Intact (static)', y: rnd(alloc.thrusters.map((t) => t.utilisation)) }, { name: `After loss of ${thrusters[wf.i].name}`, y: rnd(wf.a.thrusters.map((t) => t.utilisation)) }, { name: 'Simulation maximum', y: rnd(sim.utilisation.map((u) => u.max)) }],
          explain: 'Load on each thruster for the mean environment, after the worst single failure, and the peak during the simulation. Bars at 1.0 are saturated.' },
      ],
      tables: [
        { title: 'Thrust allocation', columns: [{ key: 'name', label: 'Thruster' }, { key: 'T', label: 'Thrust', unit: 'kN' }, { key: 'ang', label: 'Direction', unit: 'deg' }, { key: 'u', label: 'Utilisation' }, { key: 'uf', label: 'After worst failure' }, { key: 'P', label: 'Power', unit: 'kW' }],
          rows: alloc.thrusters.map((t, i) => ({ name: t.name, T: round(t.T / 1e3, 5), ang: round(t.angleDeg ?? 0, 4), u: round(t.utilisation, 4), uf: round(wf.a.thrusters[i].utilisation, 4), P: round(pFull(thrusters[i]) * Math.min(t.utilisation, 1) ** 1.5, 5) })),
          explain: 'Static allocation for the mean environmental force: thrust and direction of the force each thruster applies to the unit, its utilisation intact and after the worst single failure, and the shaft power.' },
        { title: 'Controller gains', columns: [{ key: 'dof', label: 'Degree of freedom' }, { key: 'Kp', label: 'Kp', unit: 'N/m, N·m/rad' }, { key: 'Kd', label: 'Kd', unit: 'N·s/m, N·m·s/rad' }, { key: 'Ki', label: 'Ki', unit: 'N/(m·s), N·m/(rad·s)' }, { key: 'T', label: 'Closed-loop period', unit: 's' }],
          rows: ['Surge', 'Sway', 'Yaw'].map((dof, i) => ({ dof, Kp: round(gains.Kp[i], 5), Kd: round(gains.Kd[i], 5), Ki: round(gains.Ki[i], 5), T: round(Tcl[i], 4) })),
          explain: `${inp.controller === 'pid' ? 'PID gains from pole placement' : 'Diagonal of the LQR gain matrix with an integral term at one tenth of the bandwidth'} on the wave-filtered estimates. Use them as starting values for sea-trial tuning.` },
      ],
      notes: [
        `Environment travelling towards ${fmt(dir)}°: wind ${fmt(Math.hypot(e0.wind.Fx, e0.wind.Fy) / 1e3)} kN, current ${fmt(Math.hypot(e0.cur.Fx, e0.cur.Fy) / 1e3)} kN, wave drift ${fmt(e0.drift / 1e3)} kN; resultant ${fmt(Fenv / 1e3)} kN against ${fmt(Ttot / 1e3)} kN installed.`,
        `Static utilisation ${fmt(uInt)} intact and ${fmt(uFail)} after losing ${thrusters[wf.i].name}. Capability: ${fmt(vInt)} m/s intact, ${fmt(vFail)} m/s after single failure (worst direction).`,
        `Simulation (${inp.controller.toUpperCase()}, closed-loop period ${fmt(Tcl[0])} s): excursion ${fmt(fp.maxRadius)} m after settling, ${fmt(peak)} m including the start-up transient; allocation saturated ${fmt(100 * sim.saturatedFraction)} % of the time. Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at utilisation ${fmt(gov.uc)}: thrust ${fmt(uInt)} intact, ${fmt(uFail)} after single failure, excursion ${fmt(fp.maxRadius)} m against ${fmt(inp.footMax)} m.`,
          uInt > inp.utilIntact || uFail > inp.utilFail ? `Install about ${fmt(Math.max(need1 - Ttot / 1e3, Ttot / 1e3 * (uFail / inp.utilFail - 1), 0))} kN more thrust (total ≈ ${fmt(Math.max(need1, Ttot / 1e3 * uFail / inp.utilFail))} kN), turn the unit head-on to cut windage and drift, or lower the operating limit to ${fmt(Math.min(vInt, vFail))} m/s of wind.`
            : fp.maxRadius > inp.footMax || fp.maxPsiDevDeg > inp.psiMax ? `Tighten the control: shorten the closed-loop period towards ${fmt(Math.max(Tcl[0] * inp.footMax / Math.max(fp.maxRadius, 1e-6), 5 * inp.Tp))} s (not below about five wave periods), add wind feed-forward, and check the thrust slew rate, which limits how fast the slowly varying force can be followed.`
              : `The unit holds position with ${fmt(100 * (1 - uInt / inp.utilIntact))} % thrust margin intact. Confirm the wind, current and drift coefficients by wind-tunnel and model-basin data, include thrust losses, and repeat for the worst-case failure group from the FMEA rather than a single thruster.`,
          'DNV-ST-0111 (2021); IMCA M 140 Rev. 1; IMO MSC.1/Circ.1580 (2017); Fossen (2011) ch. 12.'),
        adv('economic', 'info', `Holding station in this environment takes about ${fmt(power)} kW of thruster power, i.e. ${fmt(fuelDay)} t of fuel per day; the DP weather limit of ${fmt(Math.min(vInt, vFail))} m/s wind sets the weather downtime of the operation.`,
          'Compare the downtime between the intact and the single-failure limit with the cost of an extra thruster: one more unit raises both limits and lowers the fuel burn because thrusters run at lower, more efficient load.',
          'Operability and weather-downtime assessment per DNV-ST-N001 (2021); whole-life costing per ISO 15686-5 (fuel as the dominant DP operating cost).',
          `${fmt(fuelDay)} t of marine gas oil per day ≈ ${fmt(fuelDay * CO2_MGO)} t CO₂ per day at the screening factor of ${CO2_MGO} t CO₂ per tonne of fuel (0.21 kg fuel per kWh assumed). A moored or thruster-assisted alternative removes most of it for long-duration station keeping.`),
      ],
    };
  },
  verify: [
    { // Wind only, travelling towards +x: F = ½·1.226·1.0·2400·20² = 588.48 kN. Eight equal azimuth thrusters placed symmetrically
      // share it equally in the minimum-norm allocation: 73.56 kN each = 0.1471 of 500 kN.
      // Capability along +x with 80 % of 8 × 500 kN: V = √(0.8·4.0e6/(½·1.226·2400)) = √2175.1 = 46.64 m/s.
      name: 'Wind-only balance: equal sharing among symmetric thrusters and the limiting wind speed',
      source: 'Static equilibrium, hand calculation; minimum-norm allocation per Fossen (2011) §12.3',
      inputs: { Hs: 0, Uc: 0, Uw: 20, headingDeg: 0, AwFront: 2400, Cwind: 1, slowFrac: 0, wfSig: 0, noise: 0,
        thrusters: [thr('1', 50, 35, 'azimuth', 500), thr('2', 38, 35, 'azimuth', 500), thr('3', 50, -35, 'azimuth', 500), thr('4', 38, -35, 'azimuth', 500), thr('5', -38, 35, 'azimuth', 500), thr('6', -50, 35, 'azimuth', 500), thr('7', -38, -35, 'azimuth', 500), thr('8', -50, -35, 'azimuth', 500)] },
      expect: [{ kpi: 'Fenv', value: 588.48, tol: 1e-4 }, { kpi: 'utilMax', value: 0.1471, tol: 2e-3 }, { kpi: 'vLimHead', value: 46.64, tol: 2e-3 }] },
    { // Barge 100 × 20 × 5 m: m = 10 250 t; surge added mass 10 % → 11 275 t; ωn = 2π/100 = 0.062832 rad/s;
      // Kp = 11.275e6·0.0039478 = 44.51 kN/m; Kd = 2·0.7·0.062832·11.275e6 − 11.275e6/150 = 991.8 − 75.2 = 916.6 kN·s/m.
      name: 'PID gains by pole placement: Kp = m·ωn², Kd = 2ζωn·m − d',
      source: 'Fossen (2011) Handbook of Marine Craft Hydrodynamics and Motion Control §12.2 — hand calculation',
      inputs: { hullType: 'barge', L: 100, B: 20, Dhull: 12, draft: 5, KG: 4, amSurge: 0.1, Tctrl: 100, zetaCtrl: 0.7, Tdamp: 150,
        thrusters: [thr('1', 40, 8, 'azimuth', 300), thr('2', 40, -8, 'azimuth', 300), thr('3', -40, 8, 'azimuth', 300), thr('4', -40, -8, 'azimuth', 300)] },
      expect: [{ kpi: 'KpSurge', value: 44.51, tol: 1e-3 }, { kpi: 'KdSurge', value: 916.6, tol: 1e-3 }, { kpi: 'Tclosed', value: 100, tol: 1e-6 }] },
    { name: 'No environment, no noise: the unit stays on the setpoint', source: 'Limiting case — zero disturbance gives zero excursion and zero thrust',
      inputs: { Hs: 0, Uc: 0, Uw: 0, slowFrac: 0, wfSig: 0, noise: 0 },
      expect: [{ kpi: 'maxRadius', value: 0, abs: 1e-9 }, { kpi: 'utilMax', value: 0, abs: 1e-9 }, { kpi: 'Fenv', value: 0, abs: 1e-9 }] },
  ],
};
