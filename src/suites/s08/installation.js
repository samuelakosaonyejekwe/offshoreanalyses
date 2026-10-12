// Pipeline installation: static S-lay over a circular stinger or J-lay from a tower, with overbend and sagbend strain,
// lay tension, touchdown geometry, a quasi-dynamic tension bound for vessel motion and the DNV-ST-F101 simplified
// laying criteria.
import { sLay, jLay, stiffenedCatenary, layDynamicTension, layCriteriaF101 } from '../../engine/lines/pipelay.js';
import { pipeInputs, pipeFrom, lengthInput, depthInput, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, linspace, sev, worst, thin,
  vesselCo2, CONFIRM, ECON_BASIS, G, MPA, KM_PER_DAY } from './_common.js';

const G_V = 'Lay vessel';
const G_M = 'Vessel motion';
const G_C = 'Acceptance criteria';
const GRADES = [[360e6, 'X52'], [415e6, 'X60'], [450e6, 'X65'], [485e6, 'X70']];

export default {
  id: 's08-installation',
  suite: 8,
  title: 'Pipeline installation and laying',
  summary: 'Solves the pipe shape between the lay vessel and the seabed and checks overbend strain, sagbend stress and lay tension for S-lay or J-lay.',
  covers: [7],
  fidelity: 'numerical',
  standards: ['DNV-ST-F101 (2021) Sec. 10 (installation) and simplified laying criteria', 'Bai & Bai, Subsea Pipelines and Risers (2005) ch. 34', 'Dixon & Rutledge (1968) stiffened catenary'],
  models: ['Nonlinear beam-column (rod) model with stinger and seabed contact', 'Natural and stiffened catenary', 'Curved beam overbend strain D/(2R) + T/EA', 'Quasi-dynamic tension bounds for harmonic top motion', 'von Mises equivalent stress in the sagbend'],
  limits: 'Static in-plane analysis with a continuous frictionless circular stinger (discrete rollers, stinger flexibility and current are not resolved) and a constant apparent pipe weight above and below the waterline. Vessel motion enters only through quasi-dynamic tension bounds, not a time-domain lay simulation. The overbend strain limits by grade are recalled values to be confirmed; reel-lay plasticity and explicit non-linear finite-element lay analysis are outside this module.',
  publishes: ['pipeline.layTension', 'install.vesselDays'],
  searchTerms: 'pipeline S-lay J-lay installation analysis overbend sagbend strain lay tension',
  inputs: [
    ...pipeInputs({ omit: ['rhoCont', 'wsGiven'] }),
    lengthInput(10, 'Pipe'),
    depthInput(120, 'Site'),
    sel('method', 'Lay method', 'slay', [['slay', 'S-lay over a stinger'], ['jlay', 'J-lay from a tower']], G_V, 'operation', 'S-lay bends the pipe over a stinger (overbend) then through the sagbend; J-lay leaves near-vertically and has a sagbend only.'),
    num('hinge', 'Height of the tensioner pipe axis above sea level', 'm', 5, 0, 30, G_V, 'geometry', 'Elevation of the stinger hinge (S-lay) or tower exit (J-lay) above the waterline.'),
    num('tension', 'Lay tension at the tensioner', 'kN', 1000, 20, 20000, G_V, 'bc', 'Tensioner force held during laying (S-lay). Must exceed the weight of the suspended pipe.', { design: true, showIf: { key: 'method', in: ['slay'] } }),
    num('Rst', 'Stinger radius', 'm', 150, 30, 600, G_V, 'geometry', 'Radius of the roller bed; sets the overbend strain D/(2R).', { design: true, showIf: { key: 'method', in: ['slay'] } }),
    num('Lst', 'Stinger length (arc)', 'm', 90, 10, 200, G_V, 'geometry', 'Arc length of the supported overbend; the pipe should lift off before the last roller.', { showIf: { key: 'method', in: ['slay'] } }),
    num('tower', 'Tower angle from vertical', 'deg', 8, 0.5, 45, G_V, 'bc', 'J-lay departure angle; the catenary is matched to it.', { showIf: { key: 'method', in: ['jlay'] } }),
    num('Tcap', 'Tensioner capacity', 'kN', 2000, 50, 30000, G_V, 'operation', 'Rated holding capacity of the vessel tensioners.'),
    num('kSeabed', 'Seabed normal stiffness', 'kPa', 1000, 20, 100000, 'Site', 'soil', 'Contact stiffness per unit length and penetration at the touchdown zone.'),
    num('heave', 'Heave amplitude at the lay point', 'm', 0.8, 0, 6, G_M, 'bc', 'Vertical motion amplitude of the stinger tip or tower top in the limiting sea state (heave plus pitch arm).'),
    num('surge', 'Surge amplitude at the lay point', 'm', 0.5, 0, 6, G_M, 'bc', 'Horizontal motion amplitude in the lay direction.'),
    num('period', 'Motion period', 's', 8, 3, 25, G_M, 'bc', 'Period of the vessel motion.'),
    sel('criterion', 'Overbend strain criterion', 1, [[1, 'Criterion I (static loading)'], [2, 'Criterion II (static + dynamic)']], G_C, 'criteria', 'DNV simplified laying criteria: limit on total overbend strain by steel grade.'),
    num('epsLimit', 'Overbend strain limit (0 = by grade)', '%', 0, 0, 2, G_C, 'criteria', `Recalled defaults: X70 0.270/0.325 %, X65 0.250/0.305 %, X60 0.230/0.290 %, X52 0.205/0.260 % for criterion I/II. ${CONFIRM}`),
    num('sagFactor', 'Sagbend allowable equivalent stress / yield', '-', 0.87, 0.5, 1, G_C, 'criteria', `DNV simplified criterion for the sagbend and stinger tip. ${CONFIRM}`),
    num('layRate', 'Lay rate', 'km/day', KM_PER_DAY, 0.2, 10, 'Campaign', 'operation', 'Average progress including welding and downtime, for vessel days.'),
    { key: 'nElem', label: 'Rod elements', type: 'int', default: 100, min: 30, max: 600, group: 'Numerics', role: 'control', help: 'Elements along the suspended and grounded pipe in the rod model.' },
  ],
  mesh: { key: 'nElem', label: 'Rod elements', levels: [50, 100, 200], qoi: 'sagStrain' },
  compute(inp, ctx) {
    const p = pipeFrom({ ...inp, rhoCont: 0 }, ctx);
    const w = p.wsEmpty;
    need(w > 0, 'The empty pipe is buoyant: it cannot be laid without added weight (concrete coating or flooding).');
    need(inp.depth > 5, 'The water depth must be more than 5 m for a lay analysis.');
    const h = inp.depth + inp.hinge;
    const pipe = { EA: p.EA, EI: p.EI, w };
    const seabed = { kn: inp.kSeabed * 1e3 };
    const isS = inp.method === 'slay';
    let r, topAngle;
    if (isS) {
      const T = inp.tension * 1e3;
      need(T > 1.05 * w * h, `The lay tension (${fmt(inp.tension)} kN) does not carry the suspended pipe: it must exceed about ${fmt(1.05 * w * h / 1e3)} kN (weight × height) for this depth.`);
      r = run(() => sLay({ pipe, D: p.D, hingeHeight: h, stingerRadius: inp.Rst, stingerLength: inp.Lst, tension: T, seabed, nElem: inp.nElem }), 'S-lay solution (try a higher tension, a longer stinger or more elements)');
      topAngle = r.departureAngleDeg;
    } else {
      r = run(() => jLay({ pipe, D: p.D, depth: h, towerAngleDeg: inp.tower, seabed, nElem: inp.nElem }), 'J-lay solution');
      topAngle = 90 - inp.tower;
    }
    const rod = r.rod;
    const overbend = isS ? r.overbendStrain : 0;
    const grade = GRADES.reduce((a, g) => (Math.abs(g[0] - p.smys) < Math.abs(a[0] - p.smys) ? g : a), GRADES[0])[1];
    const pe = p.rhoW * G * inp.depth;
    const crit = layCriteriaF101({ grade, criterion: Number(inp.criterion), strainLimit: inp.epsLimit > 0 ? inp.epsLimit / 100 : undefined, overbendStrain: Math.max(overbend, 1e-9), D: p.D, t: p.t, fy: p.smys,
      sagbendTe: r.bottomTension, sagbendMoment: r.sagbendMoment, pe, stressFactor: inp.sagFactor });
    const mAx = p.W.massSteel + p.W.massCoatings;
    const dyn = layDynamicTension({ staticTension: r.topTension, w, EA: p.EA, m: mAx, suspendedLength: r.suspendedLength, topAngleDeg: topAngle, heaveAmp: inp.heave, surgeAmp: inp.surge, period: inp.period });
    const Tdyn = r.topTension + dyn.lower; // quasi-static bound; the drag-locked bound is reported, not used for acceptance
    const checks = [
      chk('Sagbend equivalent stress', crit.sagbendStress, crit.sagbendAllowable, 'Pa', 'DNV-ST-F101 (2021) Sec. 10 simplified laying criteria',
        'Guards against yielding or local buckling in the sagbend, where bending combines with external pressure on the empty pipe. Exceeding it risks a buckle at touchdown.'),
      chk('Tensioner capacity (static + quasi-static dynamic tension)', Tdyn, inp.Tcap * 1e3, 'N', 'DNV-ST-F101 (2021) Sec. 10; vessel operating limits',
        'Guards against the lay tension, including the motion-induced part, exceeding what the tensioners can hold; beyond it the pipe slips or must be abandoned.'),
      chk('No slack at the lay point in the motion cycle', dyn.lower, r.topTension, 'N', 'Quasi-static dynamic tension (first principles)',
        'Guards against the dynamic tension amplitude exceeding the static tension, which would slacken the pipe and cause compression and buckling in the sagbend.'),
    ];
    if (isS) checks.unshift(chk('Overbend total strain', overbend, crit.overbendLimit, '-', `DNV simplified laying criterion ${Number(inp.criterion) === 1 ? 'I' : 'II'} (${grade})`,
      'Guards against excessive bending strain on the stinger, which causes residual curvature, coating damage and local buckling over the rollers.'));
    const gov = worst(checks);
    // sagbend strain against bottom tension from the stiffened catenary (independent of the rod model)
    const Hs = linspace(0.3 * r.bottomTension, 2.5 * r.bottomTension, 20);
    const sc = Hs.map((H) => stiffenedCatenary({ EI: p.EI, H, w, h, D: p.D }));
    const scNow = stiffenedCatenary({ EI: p.EI, H: r.bottomTension, w, h, D: p.D });
    const [xs, zs, ss, mom, te] = thin([rod.x, rod.z, rod.s, rod.moment, rod.Te], 300);
    const strainAlong = mom.map((m) => m * p.D / (2 * p.EI) * 100);
    const days = inp.Lkm / inp.layRate;
    const L = inp.Lkm * 1e3;
    return {
      kpis: [
        kpi('topTension', 'Top tension', r.topTension / 1e3, 'kN', 'Effective tension at the tensioner or tower: what the vessel must hold statically.', inp.Tcap, 'max'),
        kpi('bottomTension', 'Bottom tension', r.bottomTension / 1e3, 'kN', 'Horizontal tension at the touchdown point; it controls the sagbend curvature and the residual lay tension left in the pipe.'),
        kpi('overbendStrain', 'Overbend strain', overbend * 100, '%', 'Maximum total strain over the stinger (bending plus axial). Zero for J-lay.', isS ? crit.overbendLimit * 100 : undefined, 'max'),
        kpi('sagStrain', 'Sagbend bending strain', r.sagbendStrain * 100, '%', 'Maximum bending strain in the sagbend from the rod model.'),
        kpi('sagStress', 'Sagbend equivalent stress', crit.sagbendStress / MPA, 'MPa', 'von Mises stress in the sagbend from bending, effective tension and external pressure on the empty pipe.', crit.sagbendAllowable / MPA, 'max'),
        kpi('sagMoment', 'Sagbend bending moment', r.sagbendMoment / 1e3, 'kN·m', 'Largest bending moment in the sagbend; input to the local buckling check for installation.'),
        kpi('sagStrainCat', 'Sagbend strain, stiffened catenary', scNow.sagbendStrain * 100, '%', 'Independent boundary-layer catenary estimate of the sagbend strain for cross-checking the rod model.'),
        kpi('departure', 'Departure angle below horizontal', topAngle, 'deg', 'Angle at which the pipe leaves the stinger or tower.'),
        kpi('layBack', 'Lay-back to touchdown', r.layBack, 'm', 'Horizontal distance from the lay point to the touchdown point; needed for route curves and start-up.'),
        kpi('suspended', 'Suspended length', r.suspendedLength, 'm', 'Length of pipe hanging between the vessel and the seabed.'),
        kpi('daf', 'Dynamic amplification of top tension (quasi-static)', dyn.dafLower, '-', 'Peak over static top tension when the suspended pipe follows the vessel motion quasi-statically. A value of 2 means the pipe goes slack in the trough of the cycle.', 2, 'max'),
        kpi('dafUpper', 'Dynamic amplification, drag-locked upper bound', dyn.dafUpper, '-', 'Bound for a line that cannot change shape within the motion cycle, so the motion is taken by axial stretch. It is far above the quasi-static value in shallow water; a time-domain analysis lies between the two.'),
        kpi('Tdyn', 'Peak top tension with vessel motion', Tdyn / 1e3, 'kN', 'Static tension plus the quasi-static dynamic amplitude.'),
        kpi('vesselDays', 'Lay duration', days, 'days', 'Route length over the average lay rate.'),
      ],
      checks,
      plots: [
        plot('line', 'Pipe configuration during laying', 'Horizontal distance from the lay point', 'm', 'Height above seabed', 'm', [
          { name: 'Pipe axis', x: xs.map((v) => round(v)), y: zs.map((v) => round(v)) }, { name: 'Sea surface', x: [round(xs[0]), round(xs[xs.length - 1])], y: [inp.depth, inp.depth], dash: true },
        ], 'Static shape of the pipe from the vessel to the seabed: the overbend on the stinger (S-lay), the inflection, the sagbend and the touchdown point.'),
        plot('line', 'Bending strain and tension along the pipe', 'Arc length from the lay point', 'm', 'Bending strain', '%', [
          { name: 'Bending strain (M·D/2EI)', x: ss.map((v) => round(v)), y: strainAlong.map((v) => round(v, 3)) }, { name: 'Effective tension (MN ÷ 10)', x: ss.map((v) => round(v)), y: te.map((v) => round(v / 1e7, 3)) },
        ], 'Signed bending strain shows the overbend (one sign) and the sagbend (the other); tension falls from the vessel to the touchdown point by the weight of the suspended pipe.'),
        plot('line', 'Sagbend strain vs bottom tension', 'Bottom tension', 'kN', 'Sagbend bending strain', '%', [
          { name: 'Stiffened catenary', x: Hs.map((v) => round(v / 1e3)), y: sc.map((s) => round(s.sagbendStrain * 100)) }, { name: 'Rod model at the lay tension', x: [round(r.bottomTension / 1e3)], y: [round(r.sagbendStrain * 100)] },
        ], 'More tension flattens the sagbend. The curve shows how far the tension can be reduced before the sagbend strain becomes critical.'),
      ],
      tables: [tbl('Lay analysis summary', [col('q', 'Quantity'), col('v', 'Value'), col('u', 'Unit')], [
        { q: 'Lay method', v: isS ? 'S-lay' : 'J-lay', u: '' }, { q: 'Empty submerged weight', v: round(w), u: 'N/m' },
        { q: 'Overbend strain, analytic D/(2R) + T/EA', v: isS ? round(r.overbendStrainAnalytic * 100) : 0, u: '%' }, { q: 'Pipe lifts off before the stinger tip', v: isS ? (r.bearsOnTip ? 'no — bears on the tip' : 'yes') : 'n/a', u: '' },
        { q: 'Maximum roller reaction', v: isS ? round(r.maxRollerReaction / 1e3) : 0, u: 'kN/m' }, { q: 'Dynamic tension amplitude, lower / upper bound', v: `${fmt(dyn.lower / 1e3)} / ${fmt(dyn.upper / 1e3)}`, u: 'kN' },
        { q: 'Minimum tension in the cycle (upper bound)', v: round(dyn.minTensionUpper / 1e3), u: 'kN' }, { q: 'Strain limit applied', v: round(crit.overbendLimit * 100), u: '%' },
      ], 'Key figures of the lay analysis in one place, including the stinger contact condition and the bounds on dynamic tension.')],
      notes: [isS ? (r.bearsOnTip ? 'The pipe bears on the stinger tip: the curvature there is concentrated — raise the tension or lengthen the stinger.' : 'The pipe lifts off the stinger before the tip.') : 'J-lay: no overbend; the tower angle is matched to the catenary.',
        `${gov.name} governs at ${fmt(gov.uc)}. Rod model converged: ${rod.converged ? 'yes' : 'no'}.`,
        'Dynamic tension is bracketed by a quasi-static and a drag-locked bound for one harmonic motion; acceptance uses the quasi-static value, so confirm with a time-domain lay analysis.'],
      advice: [
        adv('engineering', sev(gov.uc), `${isS ? 'S-lay' : 'J-lay'} in ${fmt(inp.depth)} m: top tension ${fmt(r.topTension / 1e3)} kN, ${isS ? `overbend strain ${fmt(overbend * 100)} % (limit ${fmt(crit.overbendLimit * 100)} %), ` : ''}sagbend stress ${fmt(crit.sagbendStress / MPA)} MPa; governing utilisation ${fmt(gov.uc)}.`,
          gov.uc > 1 ? (gov.name.startsWith('Overbend') ? `Increase the stinger radius to at least ${fmt(p.D / (2 * Math.max(crit.overbendLimit - r.topTension / p.EA, 1e-5)))} m or use a strain-based engineering criticality assessment.` : 'Raise the lay tension to flatten the sagbend, or reduce the limiting sea state for laying.')
            : `The configuration is acceptable statically. Confirm with a dynamic lay analysis for the limiting sea state; the tensioners have ${fmt(inp.Tcap - Tdyn / 1e3)} kN in hand.`,
          'DNV-ST-F101 (2021) Sec. 10 simplified laying criteria (recalled limits — confirm)'),
        adv('economic', 'info', `Laying ${fmt(inp.Lkm)} km at ${fmt(inp.layRate)} km/day is ${fmt(days)} vessel days with a peak tension of ${fmt(Tdyn / 1e3)} kN.`,
          Tdyn > 0.8 * inp.Tcap * 1e3 ? 'The tension is close to the vessel capacity, which narrows the weather window; a lighter coating or a higher-capacity vessel should be compared on total campaign days.' : `The tension needs only ${fmt(Tdyn / (inp.Tcap * 1e3) * 100)} % of the tensioner capacity, so a smaller lay vessel could be considered.`,
          ECON_BASIS, `${vesselCo2(days)} Lay rate is the main lever: each 10 % gain in progress saves about ${fmt(0.1 * days)} vessel days of fuel.`),
      ],
      publish: { 'pipeline.layTension': r.topTension, 'install.vesselDays': days },
      quantities: { pipe_km: inp.Lkm, steel_t: round(p.steelKgM * L / 1e3), concrete_m3: round(p.concM3M * L), vessel_days: { pipelay: round(days, 3) } },
    };
  },
  verify: [{
    // J-lay in deep water tends to the natural catenary: T_top = H + w h and H = T_top cos α, so T_top = w h/(1 − cos α),
    // α = 90° − tower angle from the horizontal. Bare 323.9 × 14.3 mm pipe, empty: w = (0.013909·7850 − 0.082397·1025)·9.80665
    // = 242.5 N/m. h = 1000 m, tower 8°: T_top = 242.5·1000/(1 − cos 82°) = 242.5e3/0.86083 = 281.7 kN.
    name: 'J-lay top tension against the natural catenary, T = w·h/(1 − cos α)',
    source: 'Catenary statics (Irvine, Cable Structures, §1.2); hand calculation',
    inputs: { D: 323.9, t: 14.3, tCoat: 0, tConc: 0, method: 'jlay', depth: 1000, hinge: 0, tower: 8, Tcap: 2000, nElem: 100 },
    expect: [{ kpi: 'topTension', value: 281.7, tol: 0.01 }],
  }],
};
