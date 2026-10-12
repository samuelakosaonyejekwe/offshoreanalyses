// Pipeline crossing and support analysis: the crossing pipe as a beam on discrete supports and a unilateral elastic
// seabed, giving support reactions, bending moment and shear diagrams, the spans created, the vertical separation
// from the crossed line, and the mattress and rock quantities.
import { beamOnSupports, crossingSeparation, crossingQuantities, rockSizeShields } from '../../engine/pipeline/crossing.js';
import { pipeInputs, pipeFrom, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, sev, worst, thin, arr,
  vesselCo2, CONFIRM, LIMIT_NO_FE, ECON_BASIS, CO2_CONCRETE, MPA } from './_common.js';

const G_M = 'Model';
const G_X = 'Crossed line and separation';
const G_Q = 'Protection quantities';
const DEFAULT_SUPPORTS = [
  { x: 82, height: 0.55, stiffness: 0 }, { x: 100, height: 0.95, stiffness: 0 }, { x: 118, height: 0.55, stiffness: 0 },
];

export default {
  id: 's08-crossing',
  suite: 8,
  title: 'Pipeline crossing and supports',
  summary: 'Lays the pipe over crossing supports and returns the support loads, bending stress, spans created and the clearance to the crossed line.',
  covers: [8],
  fidelity: 'numerical',
  standards: ['DNV-ST-F101 (2021) §5 (crossings, minimum separation 0.3 m)', 'DNV-RP-F105 (2017) (spans created at crossings)', 'CIRIA C683 Rock Manual (2007) §5.2 (rock sizing)'],
  models: ['Euler–Bernoulli beam finite elements (Hermite cubic)', 'Unilateral support and Winkler seabed contact with lift-off (active set)', 'Shear force and bending moment diagrams', 'Shields rock stability'],
  limits: `Small-deflection elastic beam without axial force: tension stiffening during laying and effective compression in operation are not included, so spans under high compression need the free-span module. Supports act at the nearest node. Settlement is modelled as a reduced support height. ${LIMIT_NO_FE}`,
  publishes: ['pipeline.crossingMoment', 'quantities.concrete_m3', 'quantities.rock_t'],
  searchTerms: 'pipeline crossing design sleeper mattress support reaction separation',
  inputs: [
    ...pipeInputs(),
    num('Lmodel', 'Modelled length', 'm', 200, 20, 2000, G_M, 'geometry', 'Length of pipe modelled, centred on the crossing; long enough that both ends rest on the seabed.'),
    sel('ends', 'End condition of the model', 'free', [['free', 'Free (resting on the seabed)'], ['pinned', 'Pinned at seabed level'], ['fixed', 'Fixed at seabed level']], G_M, 'bc', 'Free ends on an elastic seabed represent the continuing pipeline.'),
    num('kSeabed', 'Seabed vertical stiffness (0 = no seabed)', 'kN/m²', 500, 0, 100000, G_M, 'soil', 'Winkler stiffness per metre of pipe; contact in compression only.', { uncertain: { dist: 'lognormal', cov: 0.4 } }),
    { key: 'supports', label: 'Supports', type: 'table', group: G_M, role: 'bc',
      columns: [{ key: 'x', label: 'Position along the model', unit: 'm', type: 'number' }, { key: 'height', label: 'Top elevation above seabed', unit: 'm', type: 'number' }, { key: 'stiffness', label: 'Vertical stiffness (0 = rigid)', unit: 'kN/m', type: 'number' }],
      default: DEFAULT_SUPPORTS, help: 'Sleepers, mattress stacks or berm crests carrying the pipe; lift-off is allowed. Height is where the pipe bottom rests.' },
    { key: 'seabedProfile', label: 'Imported seabed / support profile (optional)', type: 'geometry', kinds: ['polyline', 'points'], default: null, group: G_M, role: 'geometry',
      help: 'Survey profile along the pipe (x along the route, z elevation). Each point becomes a unilateral support with the seabed stiffness over its tributary length and replaces the flat seabed and the table.' },
    num('settle', 'Support settlement', 'm', 0.05, 0, 1, G_M, 'bc', 'Expected settlement of the supports, deducted from their height.'),
    num('xCross', 'Position of the crossed line', 'm', 100, 0, 2000, G_X, 'geometry', 'Where the existing pipeline or cable passes under the model.'),
    num('topExisting', 'Top of the crossed line above seabed', 'm', 0.45, -2, 5, G_X, 'geometry', 'Elevation of the top of the existing line including its coating or protection (negative if buried).'),
    num('minSep', 'Minimum vertical separation', 'm', 0.3, 0.05, 2, G_X, 'criteria', `DNV-ST-F101: 0.3 m unless otherwise agreed with the owner of the crossed line. ${CONFIRM}`),
    num('spanAllow', 'Allowable span length', 'm', 30, 2, 200, 'Acceptance criteria', 'criteria', 'From the free-span module for this pipe and location.', { bind: 'pipeline.spanMax' }),
    num('etaBend', 'Allowable bending stress / SMYS', '-', 0.5, 0.1, 1, 'Acceptance criteria', 'criteria', 'Simplified limit on the static bending stress; use the combined-loading module with the moment below for the code check.'),
    num('nMat', 'Number of concrete mattresses', '-', 8, 0, 500, G_Q, 'geometry', 'Mattresses for separation and protection.'),
    num('matL', 'Mattress length', 'm', 6, 1, 12, G_Q, 'geometry', 'Plan length of one mattress.'),
    num('matW', 'Mattress width', 'm', 3, 1, 6, G_Q, 'geometry', 'Plan width of one mattress.'),
    num('matT', 'Mattress thickness', 'm', 0.3, 0.1, 1, G_Q, 'geometry', 'Thickness of one mattress.'),
    num('rockLen', 'Rock berm length over the crossing (0 = none)', 'm', 40, 0, 2000, G_Q, 'geometry', 'Length of post-lay rock cover.'),
    num('rockCover', 'Rock cover above the pipe crown', 'm', 0.5, 0.1, 3, G_Q, 'geometry', 'Cover height of the post-lay berm.'),
    { key: 'rock', label: 'Rock type', type: 'select', source: 'materials:rock', default: 'rock-granite', group: G_Q, role: 'material', help: 'Particle density of the berm rock.' },
    num('Ubed', 'Near-bed design velocity for rock sizing', 'm/s', 1.0, 0, 5, G_Q, 'environment', 'Combined wave and current velocity at the berm for hydraulic stability.', { bind: 'site.UcBed100' }),
    { key: 'nElem', label: 'Beam elements', type: 'int', default: 200, min: 20, max: 4000, group: 'Numerics', role: 'control', help: 'Number of beam elements; supports snap to the nearest node.' },
  ],
  mesh: { key: 'nElem', label: 'Beam elements', levels: [100, 200, 400, 800], qoi: 'maxMoment' },
  compute(inp, ctx) {
    const p = pipeFrom(inp, ctx);
    need(p.ws > 0, 'The pipe is buoyant in operation and would not rest on the supports.');
    const L = inp.Lmodel, D = p.Do;
    let supports = (inp.supports || []).filter((s) => Number.isFinite(Number(s.x)) && Number.isFinite(Number(s.height)))
      .map((s) => ({ x: Number(s.x), height: Number(s.height) - inp.settle, ...(Number(s.stiffness) > 0 ? { stiffness: Number(s.stiffness) * 1e3 } : {}), unilateral: true }));
    let seabed = inp.kSeabed > 0 ? { stiffness: inp.kSeabed * 1e3, elevation: 0 } : undefined;
    const g = inp.seabedProfile;
    const P = g && (g.points || g.positions);
    if (P && P.length >= 9) {
      // survey points: x along the route (relative to the first point), z elevation relative to the lowest point
      const n = P.length / 3; const pts = [];
      for (let i = 0; i < n; i++) pts.push({ x: Math.hypot(P[3 * i] - P[0], P[3 * i + 1] - P[1]), z: P[3 * i + 2] });
      pts.sort((a, b) => a.x - b.x);
      const zMin = Math.min(...pts.map((q) => q.z));
      need(pts[n - 1].x > 1, 'The imported profile has no length along the route.');
      const k = Math.max(inp.kSeabed, 50) * 1e3;
      supports = pts.filter((q) => q.x <= L).map((q, i, a) => ({ x: q.x, height: q.z - zMin, stiffness: k * 0.5 * ((a[i + 1]?.x ?? q.x) - (a[i - 1]?.x ?? q.x)) || k, unilateral: true }));
      seabed = undefined;
    }
    for (const s of supports) need(s.x >= 0 && s.x <= L, `Support at ${fmt(s.x)} m lies outside the modelled length of ${fmt(L)} m.`);
    need(supports.length > 0 || seabed || inp.ends !== 'free', 'The model is under-constrained: with free ends it needs a seabed or at least one support to carry the pipe weight.');
    need(!(supports.length === 1 && !seabed && inp.ends === 'free'), 'The model is under-constrained: a single support with free ends and no seabed cannot hold the pipe in rotation.');
    const b = run(() => beamOnSupports({ L, EI: p.EI, w: p.ws, supports, seabed, ends: inp.ends, nElem: inp.nElem, D: p.D, E: p.E }), 'Beam solution');
    // pipe axis elevation: deflection is the pipe position relative to resting on the seabed (bottom at 0)
    const iX = Math.round(inp.xCross / L * inp.nElem);
    const bottomAtX = b.deflection[Math.max(0, Math.min(inp.nElem, iX))];
    const sep = crossingSeparation({ axisElevation: bottomAtX + D / 2, Dcross: D, topOfExisting: inp.topExisting, minSeparation: inp.minSep, settlement: 0 });
    const spans = b.spans.filter((s) => s.length > 0.5);
    const maxSpan = spans.reduce((a, s) => Math.max(a, s.length), 0);
    const active = b.reactions.filter((r) => r.active);
    const Rmax = active.reduce((a, r) => Math.max(a, r.R), 0);
    let Vmax = 0; for (const v of b.shear) Vmax = Math.max(Vmax, Math.abs(v));
    const sigAllow = inp.etaBend * p.smys;
    const rockM = (ctx?.material && ctx.material(inp.rock)) || { rho: 2650 };
    const q = crossingQuantities({ mattresses: inp.nMat > 0 ? [{ count: inp.nMat, length: inp.matL, width: inp.matW, thickness: inp.matT }] : [],
      rockBerms: inp.rockLen > 0 ? [{ D, H: inp.rockCover, crestWidth: 1, sideSlope: 2.5, length: inp.rockLen, rhoSolid: rockM.rho, rhoWater: p.rhoW }] : [], rhoWater: p.rhoW });
    const d50 = inp.Ubed > 0 ? rockSizeShields({ thetaCr: 0.035, U: inp.Ubed, f: 0.02, rhoRock: rockM.rho, rhoW: p.rhoW }).D50 : 0;
    const checks = [
      chk('Vertical separation from the crossed line', inp.minSep, Math.max(sep.separation, 1e-6), 'm', 'DNV-ST-F101 (2021) §5 (crossings)',
        'Guards against contact between the two lines after settlement: contact damages coatings, couples the cathodic protection systems and loads the crossed line.'),
      chk('Static bending stress', b.maxStress, sigAllow, 'Pa', 'Beam statics; simplified allowable (see combined-loading module)',
        'Guards against overstressing the pipe where it bends over the supports and sags into the spans either side.'),
      chk('Longest span created', maxSpan, inp.spanAllow, 'm', 'DNV-RP-F105 (2017) allowable span from the free-span assessment',
        'Guards against the crossing creating free spans long enough to vibrate in the current; if exceeded add intermediate supports.'),
    ];
    const gov = worst(checks);
    const [xs, defl, mom, shr] = thin([b.x, b.deflection, b.moment, b.shear], 400);
    const days = Math.max(0.5, inp.nMat / 16) + (q.rockMass > 0 ? q.rockMass / 1e3 / 8000 : 0);
    const heightNeed = inp.topExisting + inp.minSep; // pipe bottom elevation needed at the crossing
    return {
      kpis: [
        kpi('maxMoment', 'Maximum bending moment', b.maxMoment / 1e3, 'kN·m', 'Largest bending moment magnitude in the crossing pipe, normally hogging over the highest support. Use it as the functional moment in the combined-loading check.'),
        kpi('maxStress', 'Maximum bending stress', b.maxStress / MPA, 'MPa', 'Bending stress at the outer steel fibre at the location of the largest moment.', sigAllow / MPA, 'max'),
        kpi('maxShear', 'Maximum shear force', Vmax / 1e3, 'kN', 'Largest shear force in the pipe, next to the most heavily loaded support.'),
        kpi('Rmax', 'Largest support reaction', Rmax / 1e3, 'kN', 'Vertical load on the most heavily loaded support; design the sleeper, mattress stack and its foundation for it.'),
        kpi('nActive', 'Supports in contact', active.length, '-', 'Number of supports actually carrying the pipe; a support that is not in contact is too low to work.'),
        kpi('separation', 'Vertical separation at the crossing', sep.separation, 'm', 'Clear gap between the underside of the crossing pipe and the top of the crossed line, after settlement.', inp.minSep, 'min'),
        kpi('maxSpan', 'Longest span created', maxSpan, 'm', 'Longest length without seabed or support contact; it must be within the allowable span for vibration and stress.', inp.spanAllow, 'max'),
        kpi('liftHeight', 'Pipe bottom elevation at the crossing', bottomAtX, 'm', 'Height of the pipe underside above the seabed at the crossing point.'),
        kpi('totalReaction', 'Total vertical reaction', b.totalReaction / 1e3, 'kN', 'Sum of support and seabed reactions; it equals the weight of the modelled pipe when the solution is in equilibrium.'),
        kpi('matVolume', 'Concrete in mattresses', q.mattressVolume, 'm³', 'Total concrete volume of the mattresses.'),
        kpi('rockMass', 'Rock in the berm', q.rockMass / 1e3, 't', 'Placed mass of the post-lay rock berm.'),
        kpi('rockD50', 'Stable rock size D50', d50 * 1e3, 'mm', 'Median rock size that resists the design near-bed velocity by the Shields criterion (θ_cr = 0.035, friction factor 0.02).'),
      ],
      checks,
      plots: [
        plot('line', 'Pipe profile over the crossing', 'Distance along the model', 'm', 'Elevation of the pipe underside', 'm', [
          { name: 'Pipe underside', x: arr(xs), y: arr(defl, 1, 4) }, { name: 'Supports', x: supports.map((s) => round(s.x)), y: supports.map((s) => round(s.height)) },
          { name: 'Top of crossed line', x: [round(inp.xCross)], y: [inp.topExisting] },
        ], 'Deflected shape of the pipe over the supports with the crossed line. Look for supports the pipe does not touch and for the touchdown points that end the spans.'),
        plot('line', 'Bending moment diagram', 'Distance along the model', 'm', 'Bending moment (sagging +)', 'kN·m', [{ name: 'Moment', x: arr(xs), y: arr(mom, 1e-3) }],
          'Hogging peaks over the supports and sagging in the spans. The largest magnitude drives the bending stress and the local buckling check.'),
        plot('line', 'Shear force diagram', 'Distance along the model', 'm', 'Shear force', 'kN', [{ name: 'Shear', x: arr(xs), y: arr(shr, 1e-3) }],
          'Shear jumps by the reaction at each support; the jump size is the load that support carries.'),
      ],
      tables: [
        tbl('Support reactions', [col('x', 'Position', 'm'), col('h', 'Top elevation', 'm'), col('R', 'Reaction', 'kN'), col('contact', 'In contact')],
          b.reactions.slice(0, 60).map((r, i) => ({ x: round(r.x), h: round(supports[i]?.height ?? 0), R: round(r.R / 1e3), contact: r.active ? 'yes' : 'no' })),
          'One row per support: its position, height after settlement, the vertical load it carries and whether the pipe touches it.'),
        tbl('Spans created', [col('from', 'From', 'm'), col('to', 'To', 'm'), col('len', 'Length', 'm')], spans.map((s) => ({ from: round(s.from), to: round(s.to), len: round(s.length) })),
          'Each row is a length of pipe without contact between two touchdown or support points; compare with the allowable span.'),
      ],
      notes: [`${active.length} of ${supports.length} supports carry load; equilibrium: reactions ${fmt(b.totalReaction / 1e3)} kN against weight ${fmt(b.totalLoad / 1e3)} kN.`,
        `${gov.name} governs at ${fmt(gov.uc)}. Contact iteration converged in ${b.iterations} passes.`,
        'No axial force is included: check the spans for the operating effective compression in the free-span module.'],
      advice: [
        adv('engineering', sev(gov.uc), `Separation ${fmt(sep.separation)} m (minimum ${fmt(inp.minSep)} m), bending stress ${fmt(b.maxStress / MPA)} MPa, longest span ${fmt(maxSpan)} m; governing utilisation ${fmt(gov.uc)} (${gov.name}).`,
          gov.uc > 1 ? (gov.name.startsWith('Vertical') ? `Raise the central support so the pipe underside is at least ${fmt(heightNeed)} m above the seabed at the crossing (now ${fmt(bottomAtX)} m).` : gov.name.startsWith('Longest') ? 'Add intermediate supports or grade the approach with mattresses to shorten the spans.' : 'Lower the central support or add graded approach supports to reduce the curvature over the crossing.')
            : 'The arrangement works statically. Confirm the support settlement with the geotechnical data and check on-bottom stability and trawl loads on the raised section.',
          'DNV-ST-F101 (2021) §5 crossings; DNV-RP-F105 (2017) for the spans'),
        adv('economic', 'info', `The crossing uses ${fmt(inp.nMat)} mattresses (${fmt(q.mattressVolume)} m³ of concrete) and ${fmt(q.rockMass / 1e3)} t of rock, about ${fmt(days)} vessel days of construction support.`,
          active.length < supports.length ? `${supports.length - active.length} support(s) carry no load and can be removed or lowered, saving their installation.` : 'All supports work; the saving left is in the rock volume, which scales with the berm length and cover height.',
          ECON_BASIS, `${fmt(q.mattressVolume)} m³ of mattress concrete ≈ ${fmt(q.mattressVolume * CO2_CONCRETE)} t CO₂e (screening factor ${CO2_CONCRETE} t/m³). ${vesselCo2(days)} Mattresses and rock cover ${fmt(q.mattressArea)} m² of seabed.`),
      ],
      publish: { 'pipeline.crossingMoment': b.maxMoment, 'quantities.concrete_m3': q.mattressVolume, 'quantities.rock_t': q.rockMass / 1e3 },
      quantities: { concrete_m3: round(q.mattressVolume), rock_t: round(q.rockMass / 1e3), vessel_days: { crossing: round(days, 3) } },
    };
  },
  verify: [{
    // Simply supported beam under uniform load: M_max = w L²/8 = 1000·40²/8 = 200 kN·m at mid-span; end reactions w L/2 = 20 kN each
    // (carried by the pinned ends, so no support in the table carries load).
    name: 'Simply supported beam, M = wL²/8',
    source: 'Beam statics (Roark\'s Formulas for Stress and Strain, Table 8.1 case 2e)',
    inputs: { D: 323.9, t: 14.3, wsGiven: 1000, Lmodel: 40, ends: 'pinned', kSeabed: 0, supports: [], xCross: 20, topExisting: -5, nElem: 200, settle: 0 },
    expect: [{ kpi: 'maxMoment', value: 200, tol: 0.002 }],
  }, {
    // Vertical equilibrium on the default crossing: total reaction = submerged weight × modelled length = 1000·200 = 200 kN.
    name: 'Vertical equilibrium of the pipe on supports and seabed',
    source: 'Statics: Σ reactions = w·L',
    inputs: { wsGiven: 1000 },
    expect: [{ kpi: 'totalReaction', value: 200, tol: 0.001 }],
  }],
};
