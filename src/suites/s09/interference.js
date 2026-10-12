// Riser-to-riser interference: current deflection of two neighbouring risers with wake shielding of the
// downstream one, clearance along the length, the current at which they touch and a collision energy estimate.
import { clearanceUnderCurrent, impactVelocityEstimate } from '../../engine/lines/interference.js';
import {
  need, run, fmt, round, thin, currentInputs, currentFn, adv, check, kpi, plot, finish, steelCo2Text, vesselCo2Text, CONFIRM,
} from './_common.js';

const tt = { key: 'mode', in: ['top-tensioned'] }, fh = { key: 'mode', in: ['free-hanging'] };
const riser = (p, label, d, w, T, group) => [
  { key: `D${p}`, label: `${label}: drag diameter`, unit: 'mm', type: 'number', default: d, min: 50, max: 3000, group, role: 'geometry', help: 'Outer diameter including coating, buoyancy or strakes.' },
  { key: `w${p}`, label: `${label}: apparent weight`, unit: 'N/m', type: 'number', default: w, min: 0, max: 20000, group, role: 'load', uncertain: { dist: 'normal', cov: 0.05 } },
  { key: `T${p}`, label: `${label}: top tension`, unit: 'kN', type: 'number', default: T, min: 10, max: 50000, group, role: 'bc', showIf: tt, design: true, bind: p === 'Up' ? 'riser.ttrTopTension' : undefined },
  { key: `We${p}`, label: `${label}: end weight (submerged)`, unit: 'kN', type: 'number', default: 500, min: 1, max: 20000, group, role: 'bc', showIf: fh },
  { key: `Cd${p}`, label: `${label}: drag coefficient`, unit: '-', type: 'number', default: 1.2, min: 0.3, max: 3, group, role: 'load', uncertain: { dist: 'normal', cov: 0.15 }, help: 'Multiply by the VIV drag amplification where the riser vibrates.' },
  { key: `m${p}`, label: `${label}: mass incl. added mass`, unit: 'kg/m', type: 'number', default: 250, min: 5, max: 5000, group, role: 'load' },
];

export default {
  id: 's09-interference', suite: 9, title: 'Riser interference and clashing',
  summary: 'Clearance between two neighbouring risers in steady current with wake shielding, the location of closest approach, the current that causes contact and the impact energy.',
  covers: [9], fidelity: 'analytical',
  standards: ['DNV-RP-F203 (2009) riser interference', 'DNV-RP-C205 (2019) §6.2, §6.10 wake effects', 'API RP 2RD (1998, R2013) §4.6 interference'],
  models: ['Tensioned-string current deflection (Sparks)', 'Huse wake velocity deficit', 'Blevins wake drag model', 'Morison drag', 'Work–energy collision estimate', 'Fixed-point coupling of deflection and shielding'],
  limits: 'Quasi-static, in the plane of the current, small-angle tensioned-string deflection without bending stiffness; no wake-induced lift, wake-induced oscillation or VIV (apply drag amplification through Cd). Both risers have the same length and end elevations. The impact estimate is an undamped upper bound, not a contact–impact simulation; arrays of more than two risers and 3-D current directions need a dedicated model.',
  publishes: ['riser.minClearance'],
  searchTerms: 'riser interference clashing wake shielding Huse clearance',
  inputs: [
    { key: 'mode', label: 'Riser type', type: 'select', default: 'top-tensioned', group: 'Configuration', role: 'geometry', options: [{ value: 'top-tensioned', label: 'Top-tensioned (pinned both ends)' }, { value: 'free-hanging', label: 'Free-hanging (pinned at the top, weight at the bottom)' }] },
    { key: 'L', label: 'Riser length', unit: 'm', type: 'number', default: 980, min: 20, max: 3600, group: 'Configuration', role: 'geometry', bind: 'site.depth' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 1000, min: 20, max: 3600, group: 'Configuration', role: 'environment', bind: 'site.depth' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Configuration', role: 'environment', bind: 'site.rhoW' },
    ...riser('Up', 'Upstream riser', 340, 920, 1260, 'Upstream riser'),
    ...riser('Dn', 'Downstream riser', 340, 920, 1260, 'Downstream riser'),
    { key: 'spaceTop', label: 'Centre spacing at the top', unit: 'm', type: 'number', default: 4.5, min: 0.3, max: 100, group: 'Support arrangement', role: 'bc', design: true },
    { key: 'spaceBot', label: 'Centre spacing at the bottom', unit: 'm', type: 'number', default: 4.5, min: 0.3, max: 200, group: 'Support arrangement', role: 'bc', design: true },
    { key: 'transverse', label: 'Transverse offset from the wake centreline', unit: 'm', type: 'number', default: 0, min: 0, max: 50, group: 'Support arrangement', role: 'bc', help: '0 = downstream riser directly in line with the current (most shielding).' },
    ...currentInputs({ Us: 0.8, Ub: 0.1 }),
    { key: 'wakeModel', label: 'Wake model', type: 'select', default: 'huse', group: 'Wake', role: 'criteria', options: [{ value: 'huse', label: 'Huse (velocity deficit)' }, { value: 'blevins', label: 'Blevins (wake drag)' }, { value: 'none', label: 'No shielding' }] },
    { key: 'k1', label: 'Huse wake constant k1', unit: '-', type: 'number', default: 0.25, min: 0.1, max: 0.5, group: 'Wake', role: 'criteria', showIf: { key: 'wakeModel', in: ['huse'] }, help: `Wake half-width constant (0.25 for a smooth cylinder). ${CONFIRM}` },
    { key: 'k2', label: 'Huse wake constant k2', unit: '-', type: 'number', default: 1.0, min: 0.5, max: 1.5, group: 'Wake', role: 'criteria', showIf: { key: 'wakeModel', in: ['huse'] }, help: `Centreline deficit constant. ${CONFIRM}` },
    { key: 'clearMin', label: 'Required minimum clearance', unit: 'm', type: 'number', default: 0.3, min: 0, max: 20, group: 'Acceptance criteria', role: 'criteria', help: 'Surface-to-surface gap to keep in the design current (0 = no-contact criterion). Allow for VIV amplitude and installation tolerance.' },
    { key: 'n', label: 'Stations along the riser', type: 'int', default: 200, min: 40, max: 1600, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'n', label: 'Stations', levels: [100, 200, 400, 800], qoi: 'minClearance' },

  compute(inp) {
    // engineering units → SI: mm → m, kN → N
    const L = inp.L; need(L <= inp.depth + 40, `Riser length ${L} m exceeds the water depth plus deck elevation.`);
    const cf = currentFn(inp, inp.depth), n = inp.n, sT = [], spd = [];
    for (let i = 0; i <= 40; i++) { sT.push(L * i / 40); spd.push(cf(L * i / 40)); }
    const def = (p) => ({ L, w: inp[`w${p}`], d: inp[`D${p}`] / 1e3, Cd: inp[`Cd${p}`], mode: inp.mode, ...(inp.mode === 'top-tensioned' ? { topTension: inp[`T${p}`] * 1e3 } : { endWeight: inp[`We${p}`] * 1e3 }) });
    if (inp.mode === 'top-tensioned') for (const p of ['Up', 'Dn']) need(inp[`T${p}`] * 1e3 > 1.02 * inp[`w${p}`] * L, `${p === 'Up' ? 'Upstream' : 'Downstream'} top tension ${inp[`T${p}`]} kN does not exceed the apparent weight ${fmt(inp[`w${p}`] * L / 1e3)} kN: the base would be in compression.`);
    const spacing = Array.from({ length: n + 1 }, (_, i) => inp.spaceBot + (inp.spaceTop - inp.spaceBot) * i / n);
    const rMin = 0.5 * (inp.DUp + inp.DDn) / 1e3;
    need(Math.min(inp.spaceTop, inp.spaceBot) > rMin, `Centre spacing must exceed the sum of the radii (${fmt(rMin)} m): the risers overlap as installed.`);
    const wake = inp.wakeModel === 'huse' ? { k1: inp.k1, k2: inp.k2 } : undefined;
    const solve = (scale) => run(() => clearanceUnderCurrent({ upstream: def('Up'), downstream: def('Dn'), current: { s: sT, speed: spd.map((v) => v * scale) }, spacing, transverseOffset: inp.transverse, wakeModel: inp.wakeModel, wake, n, rho: inp.rhoW, maxIter: 400 }), 'Clearance analysis failed');
    const r = solve(1);
    need(r.converged, 'The deflection–shielding iteration did not converge; the risers are probably in hard contact. Increase spacing or tension.');
    need(r.upstream.maxAngle < 0.35 && r.downstream.maxAngle < 0.35, `Deflection angles reach ${fmt(Math.max(r.upstream.maxAngle, r.downstream.maxAngle) * 57.3)}°, beyond the small-angle range of the tensioned-string model. Increase tension or reduce current.`);
    // current sweep for the contact speed
    const scales = [0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6], sweep = scales.map((f) => { try { const q = f === 1 ? r : solve(f); return { f, c: q.minClearance, ok: q.converged }; } catch { return { f, c: NaN, ok: false }; } }).filter((q) => q.ok && Number.isFinite(q.c));
    let uCrit = NaN;
    for (let i = 1; i < sweep.length; i++) if ((sweep[i - 1].c - inp.clearMin) * (sweep[i].c - inp.clearMin) <= 0 && sweep[i].c !== sweep[i - 1].c) { uCrit = inp.Us * (sweep[i - 1].f + (inp.clearMin - sweep[i - 1].c) * (sweep[i].f - sweep[i - 1].f) / (sweep[i].c - sweep[i - 1].c)); break; }
    const uCap = inp.Us * scales[scales.length - 1];
    const k = Math.round(r.sAtMin / L * n), gap0 = spacing[k] - rMin;
    // collision screening: lateral stiffness per metre k ≈ 8T/L² (taut string at mid-span) and the uniform closing
    // load that gives the computed static relative deflection, applied suddenly (current onset / wake switching)
    const Tloc = Math.max(1, Math.min(r.upstream.Te[k], r.downstream.Te[k])), kLat = 8 * Tloc / (L * L);
    const fNet = kLat * Math.max(0, gap0 - r.minClearance);
    const imp = impactVelocityEstimate({ force: fNet, gap: Math.max(gap0, 0), mass: Math.min(inp.mUp, inp.mDn), stiffness: kLat });
    const minDF = Math.min(...r.dragFactor);

    const kpis = [
      kpi('minClearance', 'Minimum clearance (surface to surface)', r.minClearance, 'm', 'Smallest gap between the two riser surfaces along the length in the design current. Zero or negative means the risers touch.', inp.clearMin, 'min'),
      kpi('sAtMin', 'Height of closest approach above the lower end', r.sAtMin, 'm', 'Elevation where the risers come closest. Bumpers, clamps or local coating protection belong here.'),
      kpi('yMaxUp', 'Upstream riser maximum deflection', r.upstream.yMax, 'm', 'Largest downstream deflection of the upstream riser, which sees the undisturbed current.'),
      kpi('yMaxDn', 'Downstream riser maximum deflection', r.downstream.yMax, 'm', 'Largest deflection of the downstream riser. It is smaller than the upstream one because it sits in the wake, and that difference closes the gap.'),
      kpi('shielding', 'Lowest drag shielding factor', minDF, '-', 'Smallest ratio of downstream drag to free-stream drag along the length. 1.0 is no shielding; around 0.5–0.8 is typical a few diameters behind a cylinder.'),
      kpi('angleUp', 'Upstream riser maximum angle', r.upstream.maxAngle * 57.2958, 'deg', 'Largest inclination from vertical of the upstream riser; the string model is valid for small angles (below about 12°).'),
      kpi('uContact', 'Surface current at which the clearance limit is reached', Number.isFinite(uCrit) ? uCrit : uCap, 'm/s', Number.isFinite(uCrit) ? 'Surface speed of the same profile shape at which the clearance falls to the required minimum. Compare with the extreme current to see the margin.' : `The clearance limit is not crossed anywhere in the swept range; the value shown is the top of the sweep (${fmt(uCap)} m/s).`),
      kpi('impactVelocity', 'Impact velocity estimate', imp.velocity, 'm/s', 'Undamped work–energy estimate of the closing velocity if the load that causes the static relative deflection is applied suddenly (dynamic overshoot up to twice the static value). Zero means the risers do not reach each other even with overshoot.'),
      kpi('impactEnergy', 'Impact energy per metre', imp.energy, 'J/m', 'Kinetic energy per metre of riser at contact from the same estimate, for comparison with the denting or coating-damage resistance of the pipe.'),
    ];
    const checks = [check('Clearance in the design current', inp.clearMin + Math.max(0, gap0 - r.minClearance), Math.max(gap0, 1e-9), 'm', 'DNV-RP-F203 (2009) §2 design principles (no collision / acceptable collision)', 'Compares the gap consumed by the relative deflection plus the required margin with the installed gap at the point of closest approach. Above 1.0 the risers touch or come closer than allowed.')];

    const [sx, cl, df] = thin([r.s, r.clearance, r.dragFactor], 300);
    const yu = thin([r.s, r.upstream.y], 300)[1], yd = thin([r.s, r.downstream.y], 300)[1], sp = thin([r.s, spacing], 300)[1];
    const plots = [
      plot('Deflected risers', 'Horizontal position (current direction)', 'm', 'Height above lower end', 'm', [{ name: 'Upstream riser', x: yu, y: sx }, { name: 'Downstream riser', x: yd.map((v, i) => v + sp[i]), y: sx }], 'Both riser axes in the design current, with the downstream riser plotted at its installed spacing. The gap closes where the upstream curve bows towards the downstream one.'),
      plot('Clearance along the length', 'Height above lower end', 'm', 'Surface-to-surface clearance', 'm', [{ name: 'Clearance', x: sx, y: cl }, { name: 'Required minimum', x: [0, L], y: [inp.clearMin, inp.clearMin], dash: true }], 'Gap between the riser surfaces against elevation. The minimum of the curve relative to the dashed line is the interference margin.'),
      plot('Wake shielding of the downstream riser', 'Height above lower end', 'm', 'Drag factor', '-', [{ name: 'Downstream drag / free-stream drag', x: sx, y: df }], 'Reduction of the drag on the downstream riser caused by the wake of the upstream one. Lower values mean more shielding and therefore more relative deflection.'),
      plot('Clearance against current speed', 'Surface current speed', 'm/s', 'Minimum clearance', 'm', [{ name: 'Minimum clearance', x: sweep.map((q) => q.f * inp.Us), y: sweep.map((q) => q.c) }, { name: 'Required minimum', x: [sweep[0].f * inp.Us, sweep[sweep.length - 1].f * inp.Us], y: [inp.clearMin, inp.clearMin], dash: true }], 'Minimum clearance as the whole current profile is scaled. Where the curve crosses the dashed line is the current speed that exhausts the margin.'),
    ];
    const tables = [{ title: 'Clearance profile', columns: [{ key: 's', label: 'Height', unit: 'm' }, { key: 'U', label: 'Current', unit: 'm/s' }, { key: 'yu', label: 'Upstream deflection', unit: 'm' }, { key: 'yd', label: 'Downstream deflection', unit: 'm' }, { key: 'df', label: 'Drag factor', unit: '-' }, { key: 'c', label: 'Clearance', unit: 'm' }],
      rows: Array.from({ length: 21 }, (_, i) => { const j = Math.round(i * n / 20); return { s: round(r.s[j], 5), U: round(cf(r.s[j]), 3), yu: round(r.upstream.y[j], 4), yd: round(r.downstream.y[j], 4), df: round(r.dragFactor[j], 3), c: round(r.clearance[j], 4) }; }),
      explain: 'Deflections, shielding factor and clearance at 21 elevations, for the interference report and for locating bumpers or spacers.' }];
    const uc = checks[0].uc;
    // tension increase on the upstream riser that would recover the shortfall: relative deflection ∝ 1/T
    const rel = gap0 - r.minClearance, tReq = rel > 0 && gap0 - inp.clearMin > 0 ? inp.TUp * rel / (gap0 - inp.clearMin) : NaN;
    const advice = [
      adv('engineering', uc > 1 ? 'fail' : uc > 0.8 ? 'warn' : 'ok',
        `Minimum clearance ${fmt(r.minClearance)} m at ${fmt(r.sAtMin)} m above the lower end (installed gap ${fmt(gap0)} m, required ${inp.clearMin} m); shielding factor down to ${fmt(minDF)}.`,
        uc > 1 ? (inp.mode === 'top-tensioned' && Number.isFinite(tReq) ? `Raise the upstream top tension to about ${fmt(tReq)} kN (relative deflection scales with 1/tension), increase spacing by ${fmt(inp.clearMin - r.minClearance)} m, or stagger the risers across the current.` : `Increase spacing by at least ${fmt(inp.clearMin - r.minClearance)} m, add weight to the downstream line, or accept contact and design bumpers for ${fmt(imp.energy)} J/m.`)
          : `Margin exists up to about ${fmt(Number.isFinite(uCrit) ? uCrit : uCap)} m/s surface current. Re-run with the VIV-amplified drag on the upstream riser and with the lightest-contents downstream case, which is the usual governing combination.`,
        `DNV-RP-F203 (2009) riser interference; wake constants per Huse (1993) / Blevins (2005). ${CONFIRM}`),
      adv('economic', 'info',
        `Spacing is ${inp.spaceTop} m at the top. Each extra metre of well-bay or porch spacing is deck and hull steel for the life of the facility, whereas ${Number.isFinite(tReq) ? `${fmt(Math.max(0, tReq - inp.TUp))} kN of extra tension` : 'extra tension or ballast'} is tensioner capacity or payload.`,
        'Settle interference by tension and layout during concept selection; retrofitting spacers or bumpers offshore needs a vessel campaign and riser downtime.',
        'Whole-life costing (ISO 15686-5); downtime valuation per the field production forecast.',
        `A retrofit campaign of about 3 vessel days: ${vesselCo2Text(3)} A collision that dents a riser means replacing joints: ${steelCo2Text(5)} for a 5 t joint string.`),
    ];
    return finish({
      kpis, checks, plots, tables,
      notes: [`${inp.wakeModel === 'none' ? 'No wake shielding' : `${inp.wakeModel === 'huse' ? 'Huse' : 'Blevins'} wake model`}; coupling converged in ${r.iterations} iterations.`, r.collision ? 'Contact is predicted in the design current.' : `No contact in the design current; closest approach at ${fmt(r.sAtMin)} m.`],
      advice,
      publish: { 'riser.minClearance': r.minClearance },
    });
  },

  verify: [
    {
      // Weightless tensioned string in uniform current: parabola with mid-span deflection y = f·L²/(8T),
      // f = ½ρ·Cd·D·U² = 0.5·1025·1.2·0.34·1² = 209.1 N/m, L = 980 m, T = 1260 kN → y = 19.92 m.
      name: 'Tensioned string under uniform drag: y_mid = f·L²/(8T)', source: 'Sparks, Fundamentals of Marine Riser Mechanics (2007) Ch. 6; classical taut-string solution',
      inputs: { wUp: 0, wDn: 0, Us: 1, Ub: 1, wakeModel: 'none', spaceTop: 30, spaceBot: 30, n: 400 },
      expect: [{ kpi: 'yMaxUp', value: 0.5 * 1025 * 1.2 * 0.34 * 980 ** 2 / (8 * 1260e3), tol: 0.002 },
        // identical unshielded risers deflect together: clearance = spacing − diameter
        { kpi: 'minClearance', value: 30 - 0.34, tol: 1e-4 }],
    },
  ],
};
