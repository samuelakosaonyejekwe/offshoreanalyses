// Subsea equipment support and load distribution: weight take-off, centre of gravity and support reactions.
import { contactPressure } from '../../engine/geotech/mudmat.js';
import { G, need, run, fmt, round, adv, check, steelCo2Text, matInputs } from './_common.js';

const DEFAULT_EQUIP = [
  { name: 'Support structure and mudmat', mass: 55, vol: 7.0, x: 0, y: 0, z: 1.2 },
  { name: 'Manifold headers and piping', mass: 32, vol: 4.6, x: 0.4, y: 0.5, z: 1.8 },
  { name: 'Valves and actuators (6 off)', mass: 14, vol: 1.9, x: -0.6, y: 1.5, z: 2.1 },
  { name: 'Connectors and hubs (4 off)', mass: 8, vol: 1.0, x: 3.8, y: 0, z: 1.5 },
  { name: 'Subsea control module', mass: 3, vol: 0.9, x: -2.5, y: -3.0, z: 2.6 },
  { name: 'Protection roof and hatches', mass: 8, vol: 2.6, x: 0, y: 0, z: 3.8 },
];
const DEFAULT_SUPPORTS = [
  { name: 'S1', x: -3.5, y: -5.0, k: 1 }, { name: 'S2', x: 3.5, y: -5.0, k: 1 }, { name: 'S3', x: 3.5, y: 5.0, k: 1 }, { name: 'S4', x: -3.5, y: 5.0, k: 1 },
  { name: 'S5', x: 0, y: -5.0, k: 1 }, { name: 'S6', x: 0, y: 5.0, k: 1 },
];

export default {
  id: 's04-equipment-load-distribution', suite: 4, title: 'Equipment support and load distribution',
  summary: 'Weight and buoyancy take-off of structure and equipment, the overall centre of gravity, the reactions at each equipment support and the resulting mudmat corner pressures.',
  covers: [10], fidelity: 'analytical',
  standards: ['ISO 13628-1:2005 / API RP 17A (subsea production systems, structures)', 'ISO 13628-15:2011 / API 17P (subsea structures and manifolds)', 'DNV-ST-N001 (2021) §5 (weight and centre-of-gravity control)', 'ISO 19901-5:2021 (weight management)'],
  models: ['Force and moment equilibrium', 'Rigid-body load distribution to elastic supports', 'Archimedes buoyancy', 'Bearing-pressure distribution under a rigid mat'],
  limits: 'The supported equipment is a rigid body on linear vertical supports of given relative stiffness, so reactions follow a plane; flexibility of the equipment skid and of the support beams is not modelled (use the in-place frame module). Horizontal loads and dynamic amplification are not included.',
  publishes: ['subsea.mass', 'subsea.volume', 'subsea.submergedWeight', 'subsea.cogX', 'subsea.cogY', 'subsea.cogZ'],
  searchTerms: 'subsea manifold weight control centre of gravity support reactions',
  inputs: [
    { key: 'items', label: 'Structure and equipment', type: 'table', group: 'Equipment', role: 'ic',
      columns: [{ key: 'name', label: 'Item', type: 'text' }, { key: 'mass', label: 'Dry mass', unit: 't', type: 'number' }, { key: 'vol', label: 'Displaced volume', unit: 'm³', type: 'number' }, { key: 'x', label: 'x', unit: 'm', type: 'number' }, { key: 'y', label: 'y', unit: 'm', type: 'number' }, { key: 'z', label: 'z above base', unit: 'm', type: 'number' }],
      default: DEFAULT_EQUIP, help: 'Centre-of-gravity coordinates relative to the mat centre (x along B, y along L). Buoyancy modules are entered with their small mass and large volume.' },
    { key: 'supported', label: 'Rows carried by the equipment supports', type: 'text', default: '2,3,4,5', group: 'Equipment', role: 'bc', help: 'Comma-separated row numbers of the table (1-based) that sit on the support points; all rows load the mudmat.' },
    { key: 'supports', label: 'Equipment support points', type: 'table', group: 'Supports', role: 'bc',
      columns: [{ key: 'name', label: 'Support', type: 'text' }, { key: 'x', label: 'x', unit: 'm', type: 'number' }, { key: 'y', label: 'y', unit: 'm', type: 'number' }, { key: 'k', label: 'Relative stiffness', type: 'number' }],
      default: DEFAULT_SUPPORTS, help: 'At least three non-collinear supports. Relative vertical stiffness 1 for equal supports.' },
    { key: 'cont', label: 'Weight contingency', unit: '%', type: 'number', default: 5, min: 0, max: 30, group: 'Equipment', role: 'criteria', help: 'DNV-ST-N001 / ISO 19901-5 weight contingency on the calculated mass.' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'inWater', label: 'Condition', type: 'select', default: 'water', group: 'Site', role: 'ic', options: [{ value: 'water', label: 'Submerged (in place)' }, { value: 'air', label: 'In air (load-out, lift)' }] },
    ...matInputs(),
    { key: 'Rallow', label: 'Allowable load per support', unit: 'kN', type: 'number', default: 250, min: 1, max: 1e5, group: 'Acceptance criteria', role: 'criteria', design: true, help: 'Design capacity of each landing pad, guide or support beam seat.' },
    { key: 'cogTol', label: 'Allowable centre-of-gravity offset', unit: 'm', type: 'number', default: 0.5, min: 0.01, max: 10, group: 'Acceptance criteria', role: 'criteria', help: 'Offset from the mat centre (and lift point) accepted for landing tilt and sling loads.' },
  ],
  compute(inp) {
    // SI: t → kg, kN → N.
    const items = (inp.items || []).filter((r) => r.mass > 0 || r.vol > 0);
    need(items.length > 0, 'Give at least one item with a mass.');
    const sups = (inp.supports || []).filter((s) => (s.k ?? 1) > 0);
    need(sups.length >= 3, 'At least three support points are needed: fewer leave the equipment free to rock (under-constrained).');
    const f = 1 + inp.cont / 100, rho = inp.inWater === 'water' ? inp.rhoW : 0;
    const w = items.map((r) => ({ ...r, W: (r.mass * 1e3 * f - rho * (r.vol || 0)) * G }));
    const sum = (arr, fn) => arr.reduce((a, r) => a + fn(r), 0);
    const mass = sum(items, (r) => r.mass * 1e3) * f, vol = sum(items, (r) => r.vol || 0);
    const W = sum(w, (r) => r.W);
    need(W > 0, 'The assembly is buoyant in this condition: check masses and displaced volumes.');
    const cog = { x: sum(w, (r) => r.W * r.x) / W, y: sum(w, (r) => r.W * r.y) / W, z: sum(items, (r) => r.mass * r.z) / sum(items, (r) => r.mass) };
    const cogAir = { x: sum(items, (r) => r.mass * r.x) / sum(items, (r) => r.mass), y: sum(items, (r) => r.mass * r.y) / sum(items, (r) => r.mass) };
    // Items on the equipment supports.
    const idx = String(inp.supported || '').split(/[,;\s]+/).map(Number).filter((n) => Number.isInteger(n) && n >= 1 && n <= w.length);
    const sel = idx.length ? idx.map((n) => w[n - 1]) : w;
    const Ws = sum(sel, (r) => r.W);
    need(Ws > 0, 'The supported items have no net weight in this condition.');
    const xs = sum(sel, (r) => r.W * r.x) / Ws, ys = sum(sel, (r) => r.W * r.y) / Ws;
    // Rigid body on vertical springs: R_i = k_i (a + b·x_i + c·y_i); a, b, c from ΣR = W, ΣR·x = W·x, ΣR·y = W·y.
    const K = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (const s of sups) { const v = [1, s.x, s.y], k = s.k ?? 1; for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) K[i][j] += k * v[i] * v[j]; }
    const det = K[0][0] * (K[1][1] * K[2][2] - K[1][2] * K[2][1]) - K[0][1] * (K[1][0] * K[2][2] - K[1][2] * K[2][0]) + K[0][2] * (K[1][0] * K[2][1] - K[1][1] * K[2][0]);
    const scale = K[0][0] * K[1][1] * K[2][2] || 1;
    need(Math.abs(det) > 1e-9 * Math.abs(scale) && det !== 0, 'The support points are collinear (or coincide): the equipment can rotate about that line. Add a support off the line.');
    // 3 × 3 system by Cramer's rule.
    const rhs = [Ws, Ws * xs, Ws * ys];
    const det3 = (m) => m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    const abc = [0, 1, 2].map((c) => det3(K.map((row, i) => row.map((v, j) => (j === c ? rhs[i] : v)))) / det);
    const reacts = sups.map((s) => ({ name: s.name, x: s.x, y: s.y, R: (s.k ?? 1) * (abc[0] + abc[1] * s.x + abc[2] * s.y) }));
    const Rmax = Math.max(...reacts.map((r) => r.R)), Rmin = Math.min(...reacts.map((r) => r.R));
    const even = Ws / sups.length;
    // Mudmat pressure from the whole assembly at its centre of gravity.
    const off = Math.hypot(cog.x, cog.y);
    need(Math.abs(cog.x) < inp.B / 2 && Math.abs(cog.y) < inp.L / 2, 'The centre of gravity lies outside the mudmat plan: the structure would topple. Check the coordinates.');
    const cp = run(() => contactPressure({ B: inp.B, L: inp.L, V: W, ex: cog.x, ey: cog.y, n: 80 }), 'Mat contact pressure');
    const checks = [
      check('Highest support reaction', Rmax / 1e3, inp.Rallow, 'kN', 'ISO 13628-15 / API 17P §5 (structural interfaces); project support design load', 'Guards against overloading one landing pad or support seat when the equipment weight is not centred on the support pattern.'),
      check('Centre-of-gravity offset', off, inp.cogTol, 'm', 'DNV-ST-N001 (2021) §5.6 weight and CoG control; ISO 19901-5', 'Keeps the centre of gravity within the envelope assumed for lifting, landing tilt and foundation eccentricity.'),
      { name: 'Support uplift', demand: round(Math.max(-Rmin, 0) / 1e3), capacity: round(even / 1e3), unit: 'kN', uc: Rmin < 0 ? Math.min(-Rmin / even + 1, 99) : 0, ref: 'Statics: supports are compression-only unless bolted or locked', explain: 'A negative reaction means the equipment lifts off that support; the rigid-plane distribution is then no longer valid and the remaining supports carry more.' },
    ];
    const ok = checks.every((c) => c.uc <= 1);
    const heavy = w.reduce((a, r) => (r.W > a.W ? r : a));
    const worstS = reacts.reduce((a, r) => (r.R > a.R ? r : a));
    const advice = [
      ok ? adv('engineering', 'ok', `Submerged weight ${fmt(W / 1e3)} kN (${fmt(mass / 1e3)} t dry incl. ${inp.cont} % contingency); centre of gravity at (${fmt(cog.x)}, ${fmt(cog.y)}, ${fmt(cog.z)}) m; support reactions ${fmt(Rmin / 1e3)}–${fmt(Rmax / 1e3)} kN against ${fmt(even / 1e3)} kN if evenly shared.`,
        'Freeze the centre-of-gravity envelope in the weight report and confirm by weighing; re-run when vendor masses change. Check the in-air condition for load-out and lift, where buoyant items shift the centre.',
        'DNV-ST-N001 (2021) §5.6; ISO 19901-5:2021; ISO 13628-15 / API 17P',
        `Keeping the contingency at ${inp.cont} % rather than a blanket 10 % avoids designing the mat and rigging for ${fmt(mass / 1e3 / f * (0.10 - inp.cont / 100))} t that may never exist; mat and rigging steel scale with it.`)
        : adv('engineering', 'fail', `${checks.reduce((a, c) => (c.uc > a.uc ? c : a)).name} is outside its limit: support ${worstS.name} carries ${fmt(Rmax / 1e3)} kN against ${fmt(inp.Rallow)} kN allowable; centre-of-gravity offset ${fmt(off)} m against ${fmt(inp.cogTol)} m.`,
          Rmax / 1e3 > inp.Rallow ? `Move the supported equipment ${fmt(Math.hypot(xs, ys))} m towards the centre of the support pattern, add a support near (${fmt(xs)}, ${fmt(ys)}) m, or design the seats for ${fmt(Rmax / 1e3 * 1.1)} kN.` : `Shift "${heavy.name}" or add ${fmt(W * off / (Math.max(inp.B, inp.L) / 2) / 1e3)} kN of trim ballast at the opposite edge to bring the centre of gravity back.`,
          'DNV-ST-N001 (2021) §5.6; ISO 19901-5:2021; ISO 13628-15 / API 17P',
          `Re-arranging equipment costs no material; trim ballast of ${fmt(W * off / (Math.max(inp.B, inp.L) / 2) / 1e3 / G)} t would add ${steelCo2Text(W * off / (Math.max(inp.B, inp.L) / 2) / 1e3 / G)}`),
      adv('economic', 'info', `Lift weight in air is ${fmt(mass * G / 1e3)} kN (${fmt(mass / 1e3)} t). "${heavy.name}" is the heaviest contributor in this condition.`,
        'Crane-vessel selection steps at capacity thresholds: keep the dry mass with contingency below the next threshold of the intended vessel and track it in the weight report at each design revision. A late overweight finding typically means a vessel change or offshore outfitting of modules.',
        'ISO 19901-5:2021 weight management; AACE RP 18R-97 estimate classes'),
    ];
    return {
      kpis: [
        { key: 'Wsub', label: inp.inWater === 'water' ? 'Submerged weight' : 'Weight in air', unit: 'kN', value: W / 1e3, explain: 'Total weight in the selected condition with contingency: (mass − ρ·displaced volume)·g summed over all items.' },
        { key: 'mass', label: 'Dry mass with contingency', unit: 't', value: mass / 1e3, explain: 'Sum of item masses times the weight contingency: the lift weight in air.' },
        { key: 'vol', label: 'Displaced volume', unit: 'm³', value: vol, explain: 'Total volume displacing water; multiplied by the water density it is the buoyancy.' },
        { key: 'cogX', label: 'Centre of gravity x', unit: 'm', value: cog.x, explain: 'Offset of the weight resultant from the mat centre along B in the selected condition. It is the load eccentricity seen by the foundation.' },
        { key: 'cogY', label: 'Centre of gravity y', unit: 'm', value: cog.y, explain: 'Offset of the weight resultant from the mat centre along L.' },
        { key: 'cogZ', label: 'Centre of gravity height', unit: 'm', value: cog.z, explain: 'Height of the mass centre above the base; it sets the overturning arm on a slope and the stability during lowering.' },
        { key: 'Rmax', label: 'Highest support reaction', unit: 'kN', value: Rmax / 1e3, limit: inp.Rallow, sense: 'max', explain: 'Largest reaction among the equipment supports from the rigid-body distribution.' },
        { key: 'Rmin', label: 'Lowest support reaction', unit: 'kN', value: Rmin / 1e3, limit: 0, sense: 'min', explain: 'Smallest support reaction; negative means that support would lift off.' },
        { key: 'pMax', label: 'Peak mat pressure under weight', unit: 'kPa', value: cp.pMax / 1e3, explain: 'Largest contact pressure under the rigid mat from the weight resultant alone.' },
        { key: 'pMin', label: 'Lowest corner pressure', unit: 'kPa', value: cp.pMin / 1e3, explain: 'Pressure at the least loaded corner; a value near zero means the weight is close to the edge of the kern.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Support reactions', categories: reacts.map((r) => r.name), series: [{ name: 'Reaction (kN)', y: reacts.map((r) => round(r.R / 1e3)) }, { name: 'Even share (kN)', y: reacts.map(() => round(even / 1e3)) }], y: { label: 'Reaction', unit: 'kN' },
          explain: 'Reaction at each equipment support against the even share. Large differences show an eccentric equipment weight or a poorly placed support.' },
        { type: 'scatter', title: 'Plan of weights, supports and centre of gravity', x: { label: 'x along B', unit: 'm' }, y: { label: 'y along L', unit: 'm' },
          series: [{ name: 'Items', x: w.map((r) => r.x), y: w.map((r) => r.y) }, { name: 'Supports', x: sups.map((s) => s.x), y: sups.map((s) => s.y) }, { name: 'Centre of gravity', x: [round(cog.x)], y: [round(cog.y)] },
            { name: 'Mat outline', x: [-inp.B / 2, inp.B / 2, inp.B / 2, -inp.B / 2, -inp.B / 2], y: [-inp.L / 2, -inp.L / 2, inp.L / 2, inp.L / 2, -inp.L / 2] }],
          explain: 'Where each item and support sits in plan and where the resultant falls. The centre of gravity should lie near the centre of the support pattern and of the mat.' },
        { type: 'bar', title: 'Weight by item', categories: w.map((r) => r.name), series: [{ name: 'Weight in this condition (kN)', y: w.map((r) => round(r.W / 1e3)) }], y: { label: 'Weight', unit: 'kN' },
          explain: 'Contribution of each item to the total weight. Negative bars are net buoyant items.' },
      ],
      tables: [
        { title: 'Weight take-off', columns: [{ key: 'name', label: 'Item' }, { key: 'mass', label: 'Dry mass', unit: 't' }, { key: 'vol', label: 'Volume', unit: 'm³' }, { key: 'W', label: 'Weight in condition', unit: 'kN' }, { key: 'x', label: 'x', unit: 'm' }, { key: 'y', label: 'y', unit: 'm' }, { key: 'z', label: 'z', unit: 'm' }, { key: 'on', label: 'On supports' }],
          rows: w.map((r, i) => ({ name: r.name, mass: round(r.mass * f), vol: r.vol || 0, W: round(r.W / 1e3), x: r.x, y: r.y, z: r.z, on: !idx.length || idx.includes(i + 1) ? 'yes' : 'no' })),
          explain: 'One row per item with contingency applied: the basis of the weight report and of the centre-of-gravity calculation.' },
        { title: 'Support reactions', columns: [{ key: 'name', label: 'Support' }, { key: 'x', label: 'x', unit: 'm' }, { key: 'y', label: 'y', unit: 'm' }, { key: 'R', label: 'Reaction', unit: 'kN' }, { key: 'uc', label: 'Utilisation' }],
          rows: reacts.map((r) => ({ name: r.name, x: r.x, y: r.y, R: round(r.R / 1e3), uc: round(r.R / 1e3 / inp.Rallow) })),
          explain: 'Vertical reaction at each support and its utilisation against the allowable seat load.' },
      ],
      notes: [`Supported weight ${fmt(Ws / 1e3)} kN acts at (${fmt(xs)}, ${fmt(ys)}) m.`, `Centre of gravity in air: (${fmt(cogAir.x)}, ${fmt(cogAir.y)}) m; in the selected condition: (${fmt(cog.x)}, ${fmt(cog.y)}) m.`, cp.fullContact ? 'The mat stays in full contact under weight alone.' : 'The weight resultant is outside the kern: part of the mat lifts off under weight alone.'],
      advice,
      publish: { 'subsea.mass': mass, 'subsea.volume': vol, 'subsea.submergedWeight': W, 'subsea.cogX': cog.x, 'subsea.cogY': cog.y, 'subsea.cogZ': cog.z },
      quantities: {},
    };
  },
  verify: [
    { name: 'Four corner supports, eccentric load: R = W/4 ± W·ex·x/Σx²', source: 'Rigid-body statics: W = 100 t·g = 980.665 kN at x = 1 m on supports at x = ±2 m, y = ±3 m → Rmax = W(1/4 + 1·2/16) = 367.75 kN, Rmin = 122.58 kN',
      inputs: { items: [{ name: 'load', mass: 100, vol: 0, x: 1, y: 0, z: 1 }], supported: '1', cont: 0, inWater: 'air',
        supports: [{ name: 'a', x: -2, y: -3, k: 1 }, { name: 'b', x: 2, y: -3, k: 1 }, { name: 'c', x: 2, y: 3, k: 1 }, { name: 'd', x: -2, y: 3, k: 1 }] },
      expect: [{ kpi: 'Rmax', value: 367.749, tol: 1e-5 }, { kpi: 'Rmin', value: 122.583, tol: 1e-5 }, { kpi: 'cogX', value: 1, tol: 1e-9 }] },
    { name: 'Archimedes: W′ = (m − ρV)·g', source: 'Hand calculation: (120 000 − 1025·18)·9.80665 = 995 866 N',
      inputs: { items: [{ name: 'all', mass: 120, vol: 18, x: 0, y: 0, z: 1 }], supported: '1', cont: 0, inWater: 'water', Rallow: 1000 },
      expect: [{ kpi: 'Wsub', value: 995.866, tol: 1e-5 }] },
  ],
};
