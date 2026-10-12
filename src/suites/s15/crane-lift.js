// Heavy lift: hook-load build-up per DNV-ST-N001 logic and crane capacity against the capacity–radius curve.
import { G, need, round, fmt, kpi, check, adv, line, tbl, interp, dafN001, liftWeightInputs, vesselCo2Text, CONFIRM, N001, LIMIT_N001, worst } from './_common.js';

const DEFAULT_CURVE = [
  { radius: 30, capacity: 5000, hook: 95 }, { radius: 40, capacity: 4200, hook: 92 }, { radius: 50, capacity: 3300, hook: 88 },
  { radius: 60, capacity: 2500, hook: 82 }, { radius: 70, capacity: 1900, hook: 75 }, { radius: 80, capacity: 1400, hook: 66 }, { radius: 90, capacity: 1000, hook: 55 },
];

export default {
  id: 's15-crane-lift', suite: 15, title: 'Heavy lift and crane capacity',
  summary: 'Builds the static and dynamic hook load from weight, contingency, rigging and dynamic amplification, and checks it against the crane capacity–radius curve and the available hook height.',
  covers: [1], fidelity: 'analytical',
  standards: [`${N001} Sec. 16 (lifting operations)`, 'DNV-RP-N103 (2017) §9.2 (lift-off dynamics in air)', 'EN 13852-1 (offshore cranes, capacity ratings)'],
  models: ['Rigid-body equilibrium', 'Crane load-radius model (tabulated curve)', 'Dynamic amplification by weight class', 'Hook-load dynamics: DAF = 1 + vr·√(K/M)/g', 'Pendulum motion of the hook load'],
  limits: `${LIMIT_N001} The crane is represented by its tabulated capacity–radius curve only: boom and pedestal structure, slew-bearing loads, vessel heel under the lift and multi-crane load sharing are not modelled. The lift-off DAF is a single-degree-of-freedom estimate, not a coupled vessel–crane–load simulation.`,
  publishes: ['install.liftWeight', 'install.hookLoad', 'install.daf', 'install.craneUtilisation'],
  searchTerms: 'offshore heavy lift crane capacity radius dynamic amplification factor DNV-ST-N001',
  inputs: [
    ...liftWeightInputs(),
    { key: 'curve', label: 'Crane capacity–radius curve', type: 'table', group: 'Crane', role: 'bc',
      columns: [{ key: 'radius', label: 'Radius', unit: 'm', type: 'number' }, { key: 'capacity', label: 'Rated capacity', unit: 't', type: 'number' }, { key: 'hook', label: 'Hook height above deck', unit: 'm', type: 'number' }],
      default: DEFAULT_CURVE, help: 'From the crane load chart for the boom configuration and sea state used. The default is a generic 5000 t-class revolving crane, indicative only.' },
    { key: 'radius', label: 'Operating radius', unit: 'm', type: 'number', default: 55, min: 5, max: 200, group: 'Crane', role: 'bc', design: true, help: 'Largest radius reached during the lift (pick-up, slew or set-down), measured from the slew axis to the hook.' },
    { key: 'dafCurve', label: 'DAF already included in the crane curve', unit: '-', type: 'number', default: 1.1, min: 1, max: 2, group: 'Crane', role: 'criteria', help: 'Offshore crane charts are rated with a stated dynamic factor. A larger operational DAF reduces the usable static capacity in proportion.' },
    { key: 'hObj', label: 'Object height (base to lift points)', unit: 'm', type: 'number', default: 22, min: 0, max: 150, group: 'Geometry', role: 'geometry' },
    { key: 'hRig', label: 'Rigging height (lift points to hook)', unit: 'm', type: 'number', default: 36, min: 0, max: 150, group: 'Geometry', role: 'geometry', design: true, help: 'Sets the sling angle; follows from the rigging design.' },
    { key: 'hLand', label: 'Height the object base must clear (above crane deck)', unit: 'm', type: 'number', default: 18, min: -50, max: 150, group: 'Geometry', role: 'bc', help: 'Highest obstruction passed or the landing elevation, relative to the crane vessel deck.' },
    { key: 'clear', label: 'Required vertical clearance', unit: 'm', type: 'number', default: 3, min: 0, max: 20, group: 'Criteria', role: 'criteria', help: `DNV-ST-N001 asks for about 3 m between the lifted object and obstructions during an offshore lift. ${CONFIRM}` },
    { key: 'dynOn', label: 'Check lift-off DAF from hoist dynamics', type: 'bool', default: true, group: 'Lift-off dynamics', role: 'control' },
    { key: 'Khoist', label: 'Hoisting system stiffness (wire, boom, slings in series)', unit: 'kN/m', type: 'number', default: 60000, min: 100, max: 5e6, group: 'Lift-off dynamics', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.2 }, showIf: { key: 'dynOn', in: [true] } },
    { key: 'vHoist', label: 'Hoisting speed at lift-off', unit: 'm/s', type: 'number', default: 0.1, min: 0, max: 2, group: 'Lift-off dynamics', role: 'ic', showIf: { key: 'dynOn', in: [true] } },
    { key: 'vTip', label: 'Relative vertical velocity crane tip – object', unit: 'm/s', type: 'number', default: 0.25, min: 0, max: 5, group: 'Lift-off dynamics', role: 'ic', bind: 'install.craneTipVelocity', showIf: { key: 'dynOn', in: [true] },
      help: 'Single amplitude of the relative velocity between the crane tip and the barge deck the object is lifted from, in the limiting sea state (from the vessel-motion module).' },
    { key: 'hoistLen', label: 'Hoist-wire length, crane tip to hook', unit: 'm', type: 'number', default: 40, min: 2, max: 300, group: 'Lift-off dynamics', role: 'geometry', help: 'Sets the pendulum period of the suspended load.' },
    { key: 'liftHours', label: 'Duration of the lift operation incl. rigging', unit: 'h', type: 'number', default: 18, min: 1, max: 240, group: 'Operation', role: 'operation' },
  ],
  compute(inp) {
    // t → kg, kN/m → N/m. Hook loads are reported in tonnes (mass-equivalent), forces in kN.
    const Wc = inp.W * inp.fWeight, SHL = Wc + inp.Wrig;
    const dafTab = inp.dafUser > 0 ? inp.dafUser : dafN001(SHL, inp.where);
    // Lift-off from a moving deck: DAF = 1 + vr·√(K/M)/g with vr = hoist speed + relative crane-tip velocity (DNV-RP-N103 §9.2 form).
    const vr = inp.vHoist + inp.vTip;
    const dafDyn = inp.dynOn ? 1 + vr * Math.sqrt(inp.Khoist * 1e3 / (SHL * 1e3)) / G : 1;
    const daf = Math.max(dafTab, dafDyn), DHL = daf * SHL;
    const rows = (inp.curve || []).map((r) => ({ r: Number(r.radius), c: Number(r.capacity), h: Number(r.hook) })).filter((r) => Number.isFinite(r.r) && Number.isFinite(r.c)).sort((a, b) => a.r - b.r);
    need(rows.length >= 2, 'The crane curve needs at least two rows (radius, capacity).');
    need(rows.every((r) => r.c > 0 && r.r > 0), 'Crane curve radii and capacities must be positive.');
    need(inp.radius >= rows[0].r && inp.radius <= rows[rows.length - 1].r, `The operating radius ${inp.radius} m is outside the crane curve (${rows[0].r}–${rows[rows.length - 1].r} m): the lift cannot be rated there.`);
    const rs = rows.map((r) => r.r), cs = rows.map((r) => r.c), hs = rows.map((r) => (Number.isFinite(r.h) ? r.h : 0));
    const capR = interp(rs, cs, inp.radius), hookR = interp(rs, hs, inp.radius);
    // Demand on the chart: static hook load scaled when the operational DAF exceeds the DAF the chart was rated for.
    const demand = SHL * Math.max(daf / inp.dafCurve, 1);
    const ucCrane = demand / capR;
    // Largest radius at which the lift is still within the chart.
    let rMax = rows[0].r;
    if (cs[0] >= demand) { rMax = rs[rs.length - 1]; for (let i = 1; i < rs.length; i++) if (cs[i] < demand) { rMax = rs[i - 1] + (rs[i] - rs[i - 1]) * (cs[i - 1] - demand) / (cs[i - 1] - cs[i]); break; } }
    const hookReq = inp.hLand + inp.clear + inp.hObj + inp.hRig;
    const Tpend = 2 * Math.PI * Math.sqrt((inp.hoistLen + inp.hRig + inp.hObj / 2) / G);
    const checks = [
      check('Crane capacity at the operating radius', demand, capR, 't', `${N001} Sec. 16 (crane capacity); crane load chart`, 'Guards against overloading the crane boom, hoist and slew bearing. Above 1.0 the lift is outside the certified chart and must not proceed at this radius.'),
    ];
    if (hookR > 0) checks.push(check('Hook height', hookReq, hookR, 'm', `${N001} Sec. 16 (clearances)`, 'The hook must be high enough to carry object, rigging and the clearance over the highest obstruction; otherwise the object strikes the structure or the block two-blocks.'));
    const w = worst(checks), ok = w.uc <= 1;
    const days = inp.liftHours / 24;
    const spare = capR - demand;
    return {
      kpis: [
        kpi('SHL', 'Static hook load', SHL, 't', 'Weight with contingency plus rigging: the load the crane holds in still conditions.'),
        kpi('daf', 'Dynamic amplification factor used', daf, '-', 'Larger of the weight-class table value (or override) and the lift-off estimate. Multiplies the static hook load to cover vessel motion and hoisting dynamics.'),
        kpi('dafTable', 'DAF by weight class', dafTab, '-', 'Recalled DNV-ST-N001 value for this weight class and location, or the override entered.'),
        kpi('dafDyn', 'DAF from lift-off dynamics', dafDyn, '-', 'One-degree-of-freedom snatch estimate 1 + vr·√(K/M)/g for the relative velocity at lift-off. A value above the table value means the sea state or hoist speed, not the weight class, governs.'),
        kpi('DHL', 'Dynamic hook load', DHL, 't', 'Maximum hook load during the lift: static hook load × DAF. Rigging and lift points are designed from it.'),
        kpi('capR', 'Crane capacity at the operating radius', capR, 't', 'Rated capacity interpolated from the crane curve at the operating radius.'),
        kpi('ucCrane', 'Crane capacity utilisation', ucCrane, '-', 'Chart demand over rated capacity. Keep below 1.0; values above about 0.9 leave little margin for weight growth.', 1, 'max'),
        kpi('rMax', 'Largest feasible radius', rMax, 'm', 'Radius at which the chart capacity equals the demand: the reach available for positioning the crane vessel.'),
        kpi('hookReq', 'Required hook height above deck', hookReq, 'm', 'Landing or obstruction height plus clearance, object height and rigging height.', hookR > 0 ? hookR : null, 'max'),
        kpi('Tpend', 'Pendulum period of the suspended load', Tpend, 's', 'Natural period of horizontal swing about the crane tip. If it lies in the wave-period band the load will swing and tugger lines or a lower limiting sea state are needed.'),
      ],
      checks,
      plots: [
        line('Crane capacity against radius', 'Radius', 'm', 'Hook load', 't', [{ name: 'Rated capacity', x: rs, y: cs }, { name: 'Chart demand', x: [rs[0], rs[rs.length - 1]], y: [demand, demand], dash: true }, { name: 'Operating point', x: [inp.radius], y: [demand] }],
          'The operating point must lie below the capacity curve. Where the demand line crosses the curve is the largest radius the lift can reach.'),
        line('Utilisation against radius', 'Radius', 'm', 'Utilisation', '-', [{ name: 'Demand / capacity', x: rs, y: cs.map((c) => demand / c) }, { name: 'Limit', x: [rs[0], rs[rs.length - 1]], y: [1, 1], dash: true }],
          'How fast the margin is consumed as the boom is lowered. Use it to choose the vessel stand-off position.'),
      ],
      tables: [tbl('Hook-load build-up', [['k', 'Item'], ['v', 'Value'], ['u', 'Unit']], [
        ['Net weight', inp.W, 't'], ['Weight with contingency', Wc, 't'], ['Rigging weight', inp.Wrig, 't'], ['Static hook load', SHL, 't'], ['DAF', daf, '-'], ['Dynamic hook load', DHL, 't'], ['Dynamic hook load', DHL * G, 'kN'], ['Spare chart capacity', spare, 't'],
      ].map(([k, v, u]) => ({ k, v: round(v), u })), 'Each row adds one factor of the DNV-ST-N001 load build-up, so the origin of the design hook load can be traced.')],
      notes: [`${dafDyn > dafTab ? 'The lift-off dynamic estimate governs the DAF' : 'The weight-class table (or override) governs the DAF'}: ${fmt(daf)}.`, `"${w.name}" governs at ${fmt(w.uc)}.`],
      advice: [
        ok ? adv('engineering', w.uc > 0.9 ? 'warn' : 'ok', `The crane is used to ${fmt(ucCrane * 100)} % at ${inp.radius} m radius (demand ${fmt(demand)} t, capacity ${fmt(capR)} t); the lift stays within the chart out to ${fmt(rMax)} m.`,
          `Hold the weight report and rigging weight to the values used; ${fmt(spare)} t of chart capacity remains for weight growth. Confirm the DAF (${fmt(daf)}) against the governing edition and the crane vessel's motion analysis.`, `${N001} Sec. 16; crane load chart`)
          : adv('engineering', 'fail', `"${w.name}" is exceeded (${fmt(w.uc)}): demand ${fmt(demand)} t against ${fmt(capR)} t at ${inp.radius} m, required hook height ${fmt(hookReq)} m.`,
            ucCrane > 1 ? `Reduce the radius to ${fmt(rMax)} m or less, lower the limiting sea state to bring the DAF down, remove ${fmt(demand - capR)} t of weight, or use a larger crane.` : `Shorten the rigging by ${fmt(hookReq - hookR)} m (check the sling angle), reduce the radius to gain hook height, or lower the landing clearance requirement.`, `${N001} Sec. 16; crane load chart`),
        adv('economic', 'info', `The lift occupies the crane vessel for about ${fmt(days)} day(s); the chart reserve is ${fmt(spare)} t at the chosen radius.`,
          ucCrane < 0.6 ? 'The crane is lightly used: a smaller crane vessel class may be able to do the lift — check the chart of the next class down before committing to this spread.' : 'Fix the weight budget early: every tonne of late weight growth near the chart limit risks a re-rig or a change of vessel, which costs vessel days rather than steel.',
          'AACE RP 18R-97 estimate classes for marine spread day rates; schedule risk per DNV-ST-N001 Sec. 2 (operation reference period)',
          `Heavy-lift vessel time dominates the footprint of the operation: ${vesselCo2Text(days, 60)}`),
      ],
      publish: { 'install.liftWeight': Wc * 1e3 * G, 'install.hookLoad': DHL * 1e3 * G, 'install.daf': daf, 'install.craneUtilisation': ucCrane },
      quantities: { vessel_days: { heavyLift: round(days) } },
    };
  },
  verify: [
    { name: 'Hook-load build-up: SHL = W·f + Wrig, DHL = DAF·SHL',
      source: 'Hand calculation: 1200 t × 1.05 + 60 t = 1320 t; weight class 1000–2500 t offshore → DAF 1.15 (recalled DNV-ST-N001 table); DHL = 1.15 × 1320 = 1518 t. Capacity at 55 m by linear interpolation between 3300 t (50 m) and 2500 t (60 m) = 2900 t.',
      inputs: { dynOn: false }, expect: [{ kpi: 'SHL', value: 1320, tol: 1e-9 }, { kpi: 'DHL', value: 1518, tol: 1e-9 }, { kpi: 'capR', value: 2900, tol: 1e-9 }, { kpi: 'ucCrane', value: 1320 * (1.15 / 1.1) / 2900, tol: 1e-9 }] },
    { name: 'Light lift offshore: DAF = 1 + 0.25·√(100/W)', source: 'Recalled DNV-ST-N001 expression for 3–100 t: W = 76 t + 4 t rigging = 80 t → 1 + 0.25·√1.25 = 1.2795',
      inputs: { W: 76, fWeight: 1, Wrig: 4, dynOn: false }, expect: [{ kpi: 'daf', value: 1 + 0.25 * Math.sqrt(1.25), tol: 1e-9 }] },
    { name: 'Lift-off snatch: DAF = 1 + vr·√(K/M)/g', source: 'SDOF energy balance: K = 60 MN/m, M = 1320 t, vr = 0.35 m/s → 1 + 0.35·√(6e7/1.32e6)/9.80665 = 1.2406',
      inputs: {}, expect: [{ kpi: 'dafDyn', value: 1 + 0.35 * Math.sqrt(6e7 / 1.32e6) / 9.80665, tol: 1e-9 }] },
  ],
};
