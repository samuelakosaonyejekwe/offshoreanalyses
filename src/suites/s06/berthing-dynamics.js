// Time-domain berthing impact: effective vessel mass striking a nonlinear hysteretic fender mounted on a
// flexible dolphin. Gives the contact history, peak reaction, rebound and an energy / momentum audit.
import { GENERIC_FENDERS, fenderCurve, interp } from './_berthing.js';
import { need, fmt, round, adv, check, CONFIRM, NOT_HIFI, ECON_BASIS } from './_common.js';

export default {
  id: 's06-berthing-dynamics', suite: 6, title: 'Berthing impact dynamics (vessel–fender–dolphin)',
  summary: 'Transient response of vessel, fender and dolphin during a berthing: reaction history, peak fender compression, dolphin displacement, contact duration, rebound velocity and the energy balance.',
  covers: [1, 10], fidelity: 'reduced-order',
  standards: ['PIANC WG33 Guidelines for the Design of Fender Systems (2002) §4', 'PIANC WG211 (2024)', 'BS 6349-4:2014 §5–6', 'DNV-RP-C204 (2019) §3 (impact energy balance)'],
  models: ['Conservation of linear momentum and energy in impact', 'Effective berthing mass (added mass and eccentricity)', 'Nonlinear fender compression with hysteresis loop', 'Rate (velocity) correction factor', 'Two-degree-of-freedom coupled vessel motion–fender–structure response', 'Viscous structural damping', 'Explicit velocity-Verlet time integration'],
  limits: 'Reduced-order two-mass model along the line of impact: the vessel is the effective mass at the contact point (M·Cm·Ce), the dolphin a single linear spring–mass–damper. Hysteresis is rate-independent with a constant ratio between the unloading and loading branches; rate dependence is a constant factor. Vessel yaw after contact, hull flexibility, multi-fender contact and soil nonlinearity during the event are not resolved. ' + NOT_HIFI,
  publishes: ['dolphin.fenderReactionDynamic', 'dolphin.reboundVelocity'],
  searchTerms: 'berthing impact simulation fender hysteresis vessel dolphin dynamic response',
  inputs: [
    { key: 'disp', label: 'Vessel displacement', unit: 't', type: 'number', default: 75000, min: 100, max: 700000, group: 'Vessel', role: 'load', bind: 'vessel.disp' },
    { key: 'Cm', label: 'Added-mass factor Cm', unit: '-', type: 'number', default: 1.78, min: 1, max: 3, group: 'Vessel', role: 'load', help: 'From the berthing-energy module.' },
    { key: 'Ce', label: 'Eccentricity factor Ce', unit: '-', type: 'number', default: 0.58, min: 0.1, max: 1, group: 'Vessel', role: 'load', help: 'From the berthing-energy module; the effective mass at the contact point is M·Cm·Ce.' },
    { key: 'v0', label: 'Approach velocity normal to the berth', unit: 'm/s', type: 'number', default: 0.115, min: 0.005, max: 1.5, group: 'Initial conditions', role: 'ic', bind: 'dolphin.berthVelocity', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Vessel velocity at first contact; fender compression and dolphin displacement start from zero.' },
    { key: 'type', label: 'Fender type', type: 'select', default: 'cone', group: 'Fender', role: 'geometry', options: Object.entries(GENERIC_FENDERS).map(([value, g]) => ({ value, label: `${g.label} — generic curve` })) },
    { key: 'H', label: 'Fender height (or diameter)', unit: 'm', type: 'number', default: 1.3, min: 0.1, max: 4.5, group: 'Fender', role: 'geometry', design: true },
    { key: 'RR', label: 'Rated reaction per unit', unit: 'kN', type: 'number', default: 1300, min: 5, max: 20000, group: 'Fender', role: 'geometry', design: true },
    { key: 'nUnits', label: 'Units sharing the impact', unit: '-', type: 'int', default: 1, min: 1, max: 12, group: 'Fender', role: 'geometry' },
    { key: 'curve', label: 'Certified manufacturer curve (replaces the generic shape)', type: 'table', group: 'Fender', role: 'geometry', columns: [{ key: 'defl', label: 'Deflection', unit: '% of height', type: 'number' }, { key: 'reaction', label: 'Reaction', unit: 'kN', type: 'number' }], default: [], help: 'Generic curves are not manufacturer data; enter the certified curve for design.' },
    { key: 'VF', label: 'Velocity (rate) factor on reaction', unit: '-', type: 'number', default: 1.0, min: 0.8, max: 1.5, group: 'Fender', role: 'criteria', help: `Rubber reacts more stiffly at real compression rates than in slow tests; manufacturer data. ${CONFIRM}` },
    { key: 'hyst', label: 'Hysteresis: unloading reaction lower by the fraction', unit: '-', type: 'number', default: 0.5, min: 0, max: 0.95, group: 'Fender', role: 'material', help: 'Fraction of the absorbed energy dissipated in the rubber; the unloading reaction is (1 − fraction) of the loading reaction.' },
    { key: 'flex', label: 'Dolphin flexibility included', type: 'bool', default: true, group: 'Dolphin', role: 'bc', help: 'Off: the fender is backed by a rigid structure.' },
    { key: 'Kd', label: 'Dolphin lateral stiffness', unit: 'MN/m', type: 'number', default: 8.5, min: 0.05, max: 5000, group: 'Dolphin', role: 'bc', bind: 'dolphin.Kh', showIf: { key: 'flex', in: [true] }, uncertain: { dist: 'lognormal', cov: 0.25 } },
    { key: 'md', label: 'Dolphin mass moving with the head', unit: 't', type: 'number', default: 250, min: 1, max: 20000, group: 'Dolphin', role: 'bc', showIf: { key: 'flex', in: [true] }, help: 'Cap and equipment plus about a quarter of the pile and contained-water mass above fixity.' },
    { key: 'zeta', label: 'Dolphin damping ratio', unit: '%', type: 'number', default: 5, min: 0, max: 30, group: 'Dolphin', role: 'material', showIf: { key: 'flex', in: [true] } },
    { key: 'nSteps', label: 'Time steps over the estimated contact duration', unit: '-', type: 'int', default: 2000, min: 100, max: 64000, group: 'Numerical', role: 'control' },
  ],
  meshTime: { key: 'nSteps', label: 'Time steps', levels: [500, 1000, 2000, 4000], qoi: 'deflPeak' },
  compute(inp) {
    // SI: t → kg, kN → N, MN/m → N/m
    const m1 = inp.disp * 1e3 * inp.Cm * inp.Ce, v0 = inp.v0, n = inp.nUnits;
    const c = fenderCurve({ type: inp.type, H: inp.H, RR: inp.RR * 1e3, table: inp.curve });
    const kEnd = Math.max((c.R[c.R.length - 1] - c.R[c.R.length - 2]) / (c.d[c.d.length - 1] - c.d[c.d.length - 2]), c.RRated / c.dRated);
    const Rload = (d) => (d <= 0 ? 0 : n * inp.VF * (d <= c.dMax ? interp(c.d, c.R, d) : c.R[c.R.length - 1] + kEnd * (d - c.dMax)));
    const flex = !!inp.flex, k2 = inp.Kd * 1e6, m2 = inp.md * 1e3, c2 = flex ? 2 * (inp.zeta / 100) * Math.sqrt(k2 * m2) : 0;
    const E0 = 0.5 * m1 * v0 * v0;
    // contact duration estimate from the impulse 2·m·v against half the rated reaction, plus the dolphin half period
    const Rref = Math.max(0.5 * n * inp.VF * c.RRated, 1);
    const Test = 2 * m1 * v0 / Rref + (flex ? Math.PI * Math.sqrt(m1 / k2) : 0);
    // The form control sets the output step; it is subdivided when the dolphin period needs a finer step for stability.
    const dt0 = Test / inp.nSteps;
    const sub = flex ? Math.max(1, Math.ceil(dt0 * Math.sqrt(k2 / m2) / 0.4)) : 1;
    need(sub * inp.nSteps <= 4e5, `The dolphin natural period (${fmt(2 * Math.PI * Math.sqrt(m2 / k2))} s) is too short against the contact duration (about ${fmt(Test)} s) to integrate in reasonable time: increase the dolphin mass, or switch the dolphin flexibility off (it responds quasi-statically).`);
    const dt = dt0 / sub;
    // Rate-independent hysteresis: virgin loading follows R(δ), unloading follows (1 − h)·R(δ); the fender moves between
    // the two branches along a stiff elastic line of slope kT, so the force is continuous in δ and reloading retraces it.
    const kT = Math.min(100 * n * inp.VF * c.RRated / c.dRated, 0.5 * (flex ? m2 : m1) / (dt * dt));
    // Velocity-Verlet (leapfrog) integration: second order, symplectic for the conservative part.
    let x1 = 0, v1 = v0, x2 = 0, v2 = 0, F = 0, t = 0, contact = false, separated = false;
    const acc = () => [-F / m1, flex ? (F - k2 * x2 - c2 * v2) / m2 : 0];
    let Rpeak = 0, dPeak = 0, x2Peak = 0, Wf = 0, Dc = 0, imp = 0, eDolPeak = 0;
    const tt = [0], FF = [0], dd = [0], xx = [0], eK = [round(E0 / 1e3)], eF = [0], eD = [0];
    const maxSteps = 6 * inp.nSteps * sub, stride = Math.max(1, Math.ceil(maxSteps / 600));
    let [a1, a2] = acc();
    for (let k = 1; k <= maxSteps; k++) {
      const d0 = x1 - x2, F0 = F, v2a = v2;
      v1 += 0.5 * dt * a1; v2 += 0.5 * dt * a2;
      x1 += dt * v1; if (flex) x2 += dt * v2;
      const d = x1 - x2;
      if (d > 0) { const up = Rload(d); F = Math.min(Math.max(F0 + kT * (d - d0), (1 - inp.hyst) * up), up); contact = true; } else F = 0;
      [a1, a2] = acc();
      v1 += 0.5 * dt * a1; v2 += 0.5 * dt * a2;
      t += dt;
      Wf += 0.5 * (F + F0) * (d - d0); Dc += 0.5 * c2 * (v2 * v2 + v2a * v2a) * dt; imp += 0.5 * (F + F0) * dt;
      Rpeak = Math.max(Rpeak, F); dPeak = Math.max(dPeak, d); x2Peak = Math.max(x2Peak, Math.abs(x2)); eDolPeak = Math.max(eDolPeak, 0.5 * k2 * x2 * x2);
      const done = contact && d <= 0;
      if (k % stride === 0 || done) { tt.push(round(t)); FF.push(round(F / 1e3)); dd.push(round(Math.max(d, 0) * 1e3)); xx.push(round(x2 * 1e3)); eK.push(round(0.5 * m1 * v1 * v1 / 1e3)); eF.push(round(Wf / 1e3)); eD.push(round((flex ? 0.5 * k2 * x2 * x2 + 0.5 * m2 * v2 * v2 : 0) / 1e3)); }
      if (done) { separated = true; break; }
    }
    const yE = [x1, v1, x2, v2], Tc = t;
    const Eend = 0.5 * m1 * yE[1] * yE[1] + (flex ? 0.5 * m2 * yE[3] * yE[3] + 0.5 * k2 * yE[2] * yE[2] : 0) + Wf + Dc;
    // Balance errors are reported to three figures, and below 1e-9 of the initial energy or momentum (arithmetic
    // round-off, not a property of the simulation) as exactly zero, so that the figures do not depend on the platform.
    const audit = (e) => (e < 1e-9 ? 0 : Number(e.toPrecision(3)));
    const eErr = audit(Math.abs(E0 - Eend) / E0);
    const pErr = audit(Math.abs(m1 * (v0 - yE[1]) - imp) / (m1 * v0));
    const vReb = separated ? Math.max(-v1, 0) : 0;
    // peak compression energy from the loading curve, for the efficiency split
    const share = flex ? eDolPeak / E0 : 0;
    const dRated = c.dRated;
    const checks = [
      check('Peak fender compression against rated deflection', dPeak * 1e3, dRated * 1e3, 'mm', 'PIANC WG33 (2002) §4.3; BS 6349-4:2014 §6', 'Largest compression reached in the simulated berthing against the rated deflection of the unit. Beyond it the rubber is overstrained and the reaction climbs steeply.'),
      check('Energy balance error of the simulation', eErr * 100, 1, '%', 'Energy-balance diagnostic for nonlinear impact simulation (conservation of energy)', 'Initial kinetic energy against the sum of remaining kinetic energy, strain energy, fender work and damping work. An error above 1 % means the time step is too coarse for the result to be trusted.'),
      check('Momentum balance error of the simulation', pErr * 100, 1, '%', 'Impulse–momentum theorem: ∫F dt = m·Δv', 'Change of vessel momentum against the time integral of the contact force. It verifies the reaction history independently of the energy audit.'),
    ];
    const w = checks.reduce((a, q) => (q.uc > a.uc ? q : a)), ok = w.uc <= 1;
    const advice = [
      adv('engineering', ok ? 'ok' : 'fail', `Peak reaction ${fmt(Rpeak / 1e3)} kN at ${fmt(dPeak * 1e3)} mm compression (${fmt(dPeak / c.H * 100)} % of height); dolphin head moves ${fmt(x2Peak * 1e3)} mm; contact lasts ${fmt(Tc)} s${separated ? ` and the vessel rebounds at ${fmt(vReb)} m/s` : ' (vessel still in contact at the end of the simulated time)'}.`,
        ok ? `The dolphin absorbs ${fmt(share * 100)} % of the kinetic energy at its peak displacement, so the fender works ${share > 0.1 ? 'noticeably' : 'only slightly'} less than the rigid-backing calculation assumes. Use the reaction history as the load function for the structural dynamic check if the contact duration is close to the dolphin period (${flex ? fmt(2 * Math.PI * Math.sqrt(m2 / k2)) : '–'} s).`
          : (w.name.startsWith('Peak') ? `Select a fender with at least ${fmt(E0 / 1e3 / n)} kJ per unit at this velocity, or reduce the approach velocity to about ${fmt(v0 * Math.sqrt(Math.min(n * c.ERated * inp.VF / E0, 1)))} m/s.` : `Increase the number of time steps (try ${inp.nSteps * 2}).`),
        'PIANC WG33 (2002) §4; PIANC WG211 (2024); BS 6349-4:2014; DNV-RP-C204 (2019) §3',
        `Of ${fmt(E0 / 1e3)} kJ arriving, ${fmt(Wf / 1e3)} kJ is turned to heat in the rubber each berthing: hysteresis that protects the lines from rebound but ages the compound — fender life, and with it polymer replacement, follows the energy dissipated per year.`),
      adv('economic', 'info', `Rebound velocity ${fmt(vReb)} m/s returns ${fmt(0.5 * m1 * vReb * vReb / 1e3)} kJ to the vessel, which tugs or lines must control.`,
        'A fender with higher hysteresis shortens the berthing manoeuvre and tug time; compare that operational saving with the unit price and service life of high-damping compounds.', `${ECON_BASIS}; berth-cycle time assessment`),
    ];
    return {
      kpis: [
        { key: 'Rpeak', label: 'Peak reaction', unit: 'kN', value: Rpeak / 1e3, explain: 'Largest contact force during the berthing, transmitted through the fender to the dolphin and to the hull.' },
        { key: 'deflPeak', label: 'Peak fender compression', unit: 'mm', value: dPeak * 1e3, limit: dRated * 1e3, sense: 'max', explain: 'Maximum compression of the fender in the event. It should stay below the rated deflection.' },
        { key: 'xDolphin', label: 'Peak dolphin head displacement', unit: 'mm', value: x2Peak * 1e3, explain: 'Largest movement of the dolphin head during the impact (zero for rigid backing).' },
        { key: 'Tc', label: 'Contact duration', unit: 's', value: Tc, explain: 'Time from first contact to separation (or to the end of the simulation). Compare with the natural period of the dolphin: similar values mean dynamic amplification.' },
        { key: 'vReb', label: 'Rebound velocity', unit: 'm/s', value: vReb, explain: 'Velocity at which the vessel leaves the fender. High rebound loads the mooring lines and prolongs the manoeuvre.' },
        { key: 'E0', label: 'Kinetic energy at contact', unit: 'kJ', value: E0 / 1e3, explain: 'Effective kinetic energy ½·(M·Cm·Ce)·v² arriving at the fender: the normal berthing energy.' },
        { key: 'Ediss', label: 'Energy dissipated in the fender', unit: 'kJ', value: Wf / 1e3, explain: 'Net work done on the fender over the event: the area of the hysteresis loop.' },
        { key: 'shareDol', label: 'Peak energy share of the dolphin', unit: '%', value: share * 100, explain: 'Strain energy stored in the dolphin at its peak displacement as a fraction of the arriving energy. The more flexible the dolphin, the less the fender must absorb.' },
        { key: 'eErr', label: 'Energy balance error', unit: '%', value: eErr * 100, limit: 1, sense: 'max', explain: 'Relative difference between the initial kinetic energy and the audited energies at the end. A numerical-quality indicator: it should be far below 1 %. An error below 1e-9 of the initial energy is round-off and is reported as 0.' },
        { key: 'pErr', label: 'Momentum balance error', unit: '%', value: pErr * 100, limit: 1, sense: 'max', explain: 'Relative difference between the vessel momentum change and the impulse of the contact force. It should be far below 1 %; an error below 1e-9 of the initial momentum is round-off and is reported as 0.' },
        { key: 'impulse', label: 'Impulse of the contact force', unit: 'kN·s', value: imp / 1e3, explain: 'Time integral of the reaction, equal to the change of vessel momentum (arrest plus rebound).' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Reaction, fender compression and dolphin displacement in time', x: { label: 'Time from first contact', unit: 's' }, y: { label: 'Reaction (kN) / displacement (mm)', unit: '' },
          series: [{ name: 'Reaction (kN)', x: tt, y: FF }, { name: 'Fender compression (mm)', x: tt, y: dd }, { name: 'Dolphin displacement (mm)', x: tt, y: xx }],
          explain: 'Time histories through the impact. A buckling fender shows a flat-topped reaction; the drop at maximum compression is the hysteresis; the ripple on the curves is the dolphin vibrating on its own period.' },
        { type: 'line', title: 'Force–compression loop of the fender', x: { label: 'Fender compression', unit: 'mm' }, y: { label: 'Reaction', unit: 'kN' }, series: [{ name: 'Simulated path', x: dd, y: FF }],
          explain: 'Loading and unloading path followed in the event. The enclosed area is the dissipated energy; the upper branch is the catalogue curve scaled by the rate factor.' },
        { type: 'line', title: 'Energy partition in time', x: { label: 'Time', unit: 's' }, y: { label: 'Energy', unit: 'kJ' }, series: [{ name: 'Vessel kinetic', x: tt, y: eK }, { name: 'Work on fender', x: tt, y: eF }, { name: 'Dolphin (strain + kinetic)', x: tt, y: eD }],
          explain: 'Where the kinetic energy goes: into the fender, into the dolphin, and back to the vessel on rebound. The three curves plus damping losses always sum to the initial energy.' },
      ],
      tables: [{ title: 'Energy and momentum audit', columns: [{ key: 'item', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
        rows: [['Initial kinetic energy', E0 / 1e3, 'kJ'], ['Vessel kinetic energy at the end', 0.5 * m1 * yE[1] * yE[1] / 1e3, 'kJ'], ['Net work on the fender (dissipated)', Wf / 1e3, 'kJ'], ['Dolphin strain + kinetic energy at the end', (flex ? 0.5 * m2 * yE[3] * yE[3] + 0.5 * k2 * yE[2] * yE[2] : 0) / 1e3, 'kJ'], ['Dolphin damping work', Dc / 1e3, 'kJ'], ['Initial momentum', m1 * v0 / 1e3, 'kN·s'], ['Impulse of the reaction', imp / 1e3, 'kN·s'], ['Effective mass M·Cm·Ce', m1 / 1e3, 't'], ['Time step', dt, 's'], ['Branch-transition stiffness', kT / 1e6, 'MN/m']].map(([item, v, u]) => ({ item, v: round(v), u })),
        explain: 'Terms of the energy balance and of the impulse–momentum balance, to audit the simulation against conservation laws.' }],
      notes: [c.generic ? 'Generic fender curve used — replace with the certified manufacturer curve for design.' : 'Certified manufacturer curve used.', separated ? `Vessel separated after ${fmt(Tc)} s.` : 'The vessel had not separated at the end of the simulated time (it stays against the fender).', `${w.name} governs (utilisation ${fmt(w.uc)}).`],
      advice,
      publish: { 'dolphin.fenderReactionDynamic': Rpeak, 'dolphin.reboundVelocity': vReb },
      quantities: {},
    };
  },
  verify: [
    // Mass on a linear spring: m = 1e7 kg, k = 1000 kN / 0.5 m = 2e6 N/m, v = 0.2 m/s.
    // F_peak = v·√(k·m) = 0.2·√(2e13) = 894.4 kN; δ_max = v·√(m/k) = 447.2 mm; contact time π·√(m/k) = 7.025 s; elastic rebound at 0.2 m/s.
    { name: 'Linear fender on rigid backing: F = v√(km), contact time π√(m/k)', source: 'Closed-form half-cycle of a mass–spring oscillator (energy and momentum conservation)',
      inputs: { disp: 10000, Cm: 1, Ce: 1, v0: 0.2, curve: [{ defl: 0, reaction: 0 }, { defl: 50, reaction: 1000 }], H: 1, nUnits: 1, VF: 1, hyst: 0, flex: false, nSteps: 4000 },
      expect: [{ kpi: 'Rpeak', value: 894.43, tol: 0.002 }, { kpi: 'deflPeak', value: 447.21, tol: 0.002 }, { kpi: 'Tc', value: 7.0248, tol: 0.005 }, { kpi: 'vReb', value: 0.2, tol: 0.002 }, { kpi: 'eErr', value: 0, abs: 0.1 }] },
    // With the unloading reaction halved, the energy returned is half the energy absorbed: v_rebound = v·√(1 − 0.5) = 0.14142 m/s
    // (the stiff transition between the branches returns about 0.5 % more energy than the ideal loop).
    { name: 'Hysteretic fender: rebound velocity v√(1 − h)', source: 'Energy balance of a loading–unloading loop with constant force ratio',
      inputs: { disp: 10000, Cm: 1, Ce: 1, v0: 0.2, curve: [{ defl: 0, reaction: 0 }, { defl: 50, reaction: 1000 }], H: 1, nUnits: 1, VF: 1, hyst: 0.5, flex: false, nSteps: 4000 },
      expect: [{ kpi: 'vReb', value: 0.141421, tol: 0.005 }, { kpi: 'Ediss', value: 100, tol: 0.01 }] },
  ],
};
