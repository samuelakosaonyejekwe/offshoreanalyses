// Load-out by skidding or self-propelled trailers: friction and pull, support reactions, quay loading and the
// ballast sequence that keeps the barge level with the quay while the load transfers.
import { boxHydrostatics } from '../../engine/hydro/hydrostatics.js';
import { ballastPlan } from '../../engine/hydro/stability.js';
import { G, DEG, need, run, round, fmt, kpi, check, adv, line, bar, tbl, worst, vesselCo2Text, CONFIRM, N001 } from './_common.js';

const DEFAULT_SUPPORTS = [
  { name: 'A1', x: -15, y: -10 }, { name: 'A2', x: -5, y: -10 }, { name: 'A3', x: 5, y: -10 }, { name: 'A4', x: 15, y: -10 },
  { name: 'B1', x: -15, y: 10 }, { name: 'B2', x: -5, y: 10 }, { name: 'B3', x: 5, y: 10 }, { name: 'B4', x: 15, y: 10 },
];

export default {
  id: 's15-loadout', suite: 15, title: 'Load-out and skidway',
  summary: 'Gives the break-out and running pull, the support reactions on skid shoes or trailer axles, the quay loading and the barge ballast and pumping sequence needed to keep the deck level during load transfer.',
  covers: [4], fidelity: 'analytical',
  standards: [`${N001} Sec. 10 (load-out: friction, pull capacity, ballast and pump capacity, trailer stability)`, 'ISO 19901-6 (2009) §11 (load-out)', 'Biran & López-Pulido (2014) Ship Hydrostatics and Stability ch. 7 (weight shifts and trim)'],
  models: ['Rigid-body load transfer', 'Skid friction (Coulomb, static and kinetic)', 'Support reaction redistribution (rigid body on discrete supports)', 'Buoyancy equilibrium and small-angle hydrostatics', 'Bounded least-squares ballast plan', 'Trailer hydraulic three-point stability'],
  limits: 'The structure and the barge are rigid and the supports equally stiff, so reactions follow rigid-body statics: load redistribution from barge deflection, quay settlement or a step between quay and barge is not computed. The barge is a box with small-angle hydrostatics; friction coefficients, pull and pump margins are recalled DNV-ST-N001 practice exposed as inputs to be confirmed. A grounded load-out and link-span bridges are outside the module.',
  publishes: ['install.loadoutPull', 'install.loadoutBallast'],
  searchTerms: 'load-out skidding trailer SPMT ballast sequence barge quay friction',
  inputs: [
    { key: 'W', label: 'Structure weight', unit: 't', type: 'number', default: 6000, min: 10, max: 80000, group: 'Structure', role: 'load', bind: 'install.liftWeight', uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'fWeight', label: 'Weight contingency factor', unit: '-', type: 'number', default: 1.05, min: 1, max: 1.5, group: 'Structure', role: 'criteria', help: CONFIRM },
    { key: 'Ls', label: 'Structure length along the load-out direction', unit: 'm', type: 'number', default: 40, min: 2, max: 400, group: 'Structure', role: 'geometry' },
    { key: 'ex', label: 'CoG offset from the support centroid, x', unit: 'm', type: 'number', default: 1.0, min: -50, max: 50, group: 'Structure', role: 'ic', help: 'Initial condition: where the weight sits on its supports before the move starts.' },
    { key: 'ey', label: 'CoG offset from the support centroid, y', unit: 'm', type: 'number', default: 0.5, min: -50, max: 50, group: 'Structure', role: 'ic' },
    { key: 'hCog', label: 'CoG height above the support plane', unit: 'm', type: 'number', default: 14, min: 0, max: 150, group: 'Structure', role: 'geometry' },
    { key: 'supports', label: 'Support points (skid shoes or trailer groups)', type: 'table', group: 'Supports', role: 'bc',
      columns: [{ key: 'name', label: 'Name', type: 'text' }, { key: 'x', label: 'x', unit: 'm', type: 'number' }, { key: 'y', label: 'y', unit: 'm', type: 'number' }],
      default: DEFAULT_SUPPORTS, help: 'Positions of the skid shoes or trailer support groups in the structure frame (x along the load-out direction). At least three, not collinear.' },
    { key: 'method', label: 'Load-out method', type: 'select', default: 'skid', group: 'Operation', role: 'operation', options: [{ value: 'skid', label: 'Skidded on skid beams (strand jacks or winches)' }, { value: 'trailer', label: 'Self-propelled modular trailers' }] },
    { key: 'muS', label: 'Static (break-out) friction', unit: '-', type: 'number', default: 0.15, min: 0.005, max: 0.6, group: 'Skidway', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.25 }, showIf: { key: 'method', in: ['skid'] }, help: `Steel on PTFE about 0.12–0.20 at break-out, steel on greased timber 0.20–0.25, steel on steel 0.30. ${CONFIRM}` },
    { key: 'muK', label: 'Moving friction', unit: '-', type: 'number', default: 0.08, min: 0.005, max: 0.5, group: 'Skidway', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.25 }, showIf: { key: 'method', in: ['skid'] }, help: `Steel on PTFE about 0.05–0.10 once moving. ${CONFIRM}` },
    { key: 'rollRes', label: 'Trailer rolling resistance', unit: '-', type: 'number', default: 0.03, min: 0.005, max: 0.15, group: 'Trailers', role: 'bc', showIf: { key: 'method', in: ['trailer'] }, help: 'About 2–3 % on a firm level surface, more on link plates.' },
    { key: 'grade', label: 'Gradient against the direction of travel', unit: '%', type: 'number', default: 0.5, min: -5, max: 10, group: 'Skidway', role: 'bc', help: 'Includes the residual slope from barge trim and quay level difference.' },
    { key: 'pullCap', label: 'Installed pull capacity', unit: 'kN', type: 'number', default: 14000, min: 1, max: 5e5, group: 'Skidway', role: 'bc', design: true, help: 'Strand jacks, winches or trailer tractive effort.' },
    { key: 'fPull', label: 'Required margin on the break-out force', unit: '-', type: 'number', default: 1.3, min: 1, max: 3, group: 'Criteria', role: 'criteria', help: `Installed capacity over the computed break-out force. ${CONFIRM}` },
    { key: 'shoeCap', label: 'Capacity of one support (shoe or trailer group)', unit: 't', type: 'number', default: 1100, min: 1, max: 20000, group: 'Supports', role: 'bc', design: true },
    { key: 'shoeLen', label: 'Bearing length of one support on the skid beam / quay', unit: 'm', type: 'number', default: 6, min: 0.2, max: 60, group: 'Supports', role: 'geometry' },
    { key: 'quayCap', label: 'Quay / skid-beam line-load capacity', unit: 'kN/m', type: 'number', default: 2500, min: 10, max: 50000, group: 'Quay', role: 'bc', help: 'Per skid line or trailer train, from the quay structural assessment.' },
    { key: 'nAxle', label: 'Number of trailer axle lines', unit: '-', type: 'int', default: 160, min: 2, max: 2000, group: 'Trailers', role: 'geometry', showIf: { key: 'method', in: ['trailer'] } },
    { key: 'axleCap', label: 'Capacity per axle line', unit: 't', type: 'number', default: 44, min: 5, max: 80, group: 'Trailers', role: 'bc', showIf: { key: 'method', in: ['trailer'] }, help: 'Payload plus self-weight per axle line at the load-out speed (typically 36–48 t).' },
    { key: 'axleTare', label: 'Trailer self-weight per axle line', unit: 't', type: 'number', default: 4, min: 1, max: 10, group: 'Trailers', role: 'load', showIf: { key: 'method', in: ['trailer'] } },
    { key: 'hTrailer', label: 'Trailer deck height', unit: 'm', type: 'number', default: 1.5, min: 0.5, max: 3, group: 'Trailers', role: 'geometry', showIf: { key: 'method', in: ['trailer'] } },
    { key: 'tipReq', label: 'Required trailer stability angle', unit: 'deg', type: 'number', default: 7, min: 1, max: 20, group: 'Criteria', role: 'criteria', showIf: { key: 'method', in: ['trailer'] }, help: `Angle at which the CoG passes the tipping line of the hydraulic suspension triangle. ${CONFIRM}` },
    { key: 'Lb', label: 'Barge length', unit: 'm', type: 'number', default: 120, min: 30, max: 400, group: 'Barge', role: 'geometry' },
    { key: 'Bb', label: 'Barge beam', unit: 'm', type: 'number', default: 36, min: 8, max: 120, group: 'Barge', role: 'geometry' },
    { key: 'Db', label: 'Barge depth', unit: 'm', type: 'number', default: 7.6, min: 2, max: 30, group: 'Barge', role: 'geometry' },
    { key: 'T0', label: 'Barge draft at the start of load-out', unit: 'm', type: 'number', default: 4.6, min: 0.5, max: 25, group: 'Barge', role: 'ic', help: 'Ballasted so that the deck is level with the quay at the tide when the move starts.' },
    { key: 'fill0', label: 'Initial filling of the ballast tanks', unit: '%', type: 'number', default: 80, min: 0, max: 100, group: 'Barge', role: 'ic', help: 'Eight wing tanks of equal size are assumed, together 80 % of the hull volume, operated as a forward and an aft group with port–starboard transfer for heel.' },
    { key: 'xFinal', label: 'Final position of the structure CoG from barge midship', unit: 'm', type: 'number', default: 0, min: -150, max: 150, group: 'Barge', role: 'bc', help: 'Positive towards the bow; the load comes over the stern.' },
    { key: 'speed', label: 'Load-out speed', unit: 'm/h', type: 'number', default: 15, min: 1, max: 300, group: 'Operation', role: 'operation', design: true },
    { key: 'tideRate', label: 'Tide rate during the move (rising positive)', unit: 'm/h', type: 'number', default: 0.1, min: -2, max: 2, group: 'Operation', role: 'bc', bind: 'install.tideRate', help: 'The barge must be ballasted against the tide as well as against the incoming load.' },
    { key: 'pumpCap', label: 'Installed ballast pump capacity (all pumps)', unit: 'm³/h', type: 'number', default: 8000, min: 10, max: 2e5, group: 'Barge', role: 'bc', design: true },
    { key: 'fPump', label: 'Required pump capacity margin', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Criteria', role: 'criteria', help: `Installed over required pumping rate, covering a pump failure. ${CONFIRM}` },
    { key: 'nStep', label: 'Number of load-transfer steps', unit: '-', type: 'int', default: 10, min: 4, max: 40, group: 'Operation', role: 'control' },
  ],
  compute(inp) {
    // t → kg, % → fraction, m/h and m³/h kept per hour.
    const W = inp.W * inp.fWeight * 1e3, Wg = W * G, rho = 1025;
    const sup = (inp.supports || []).map((s, i) => ({ name: s.name || `S${i + 1}`, x: Number(s.x), y: Number(s.y) })).filter((s) => Number.isFinite(s.x) && Number.isFinite(s.y));
    need(sup.length >= 3, 'At least three support points are needed to carry the structure.');
    // Rigid body on equally stiff supports: R_i = W/N + My·x_i/Σx² + Mx·y_i/Σy² about the support centroid.
    const xc = sup.reduce((s, p) => s + p.x, 0) / sup.length, yc = sup.reduce((s, p) => s + p.y, 0) / sup.length;
    const Sxx = sup.reduce((s, p) => s + (p.x - xc) ** 2, 0), Syy = sup.reduce((s, p) => s + (p.y - yc) ** 2, 0), Sxy = sup.reduce((s, p) => s + (p.x - xc) * (p.y - yc), 0);
    const det = Sxx * Syy - Sxy * Sxy;
    need(det > 1e-9 * Math.max(Sxx, Syy, 1) ** 2, 'The support points are collinear: the structure is not restrained against overturning about that line. Add supports off the line.');
    const a = Wg * (inp.ex * Syy - inp.ey * Sxy) / det, b = Wg * (inp.ey * Sxx - inp.ex * Sxy) / det;
    const R = sup.map((p) => Wg / sup.length + a * (p.x - xc) + b * (p.y - yc));
    const Rmax = Math.max(...R), Rmin = Math.min(...R);
    need(Rmin > 0, 'A support lifts off (negative reaction): the centre of gravity is outside the kern of the support pattern. Rearrange the supports or correct the CoG offsets.');
    const grade = inp.grade / 100, skid = inp.method === 'skid';
    const Fbreak = Wg * ((skid ? inp.muS : inp.rollRes) + grade), Frun = Wg * ((skid ? inp.muK : inp.rollRes) + grade);
    const checks = [
      check('Pull capacity against break-out force', Fbreak * inp.fPull / 1e3, inp.pullCap, 'kN', `${N001} Sec. 10 (load-out: push/pull capacity and friction)`, 'The jacks or tractors must start the structure against static friction and gradient with margin; otherwise the move stalls half-way across the quay edge.'),
      check('Maximum support reaction', Rmax / G / 1e3, inp.shoeCap, 't', `${N001} Sec. 10 (load-out: support reactions)`, 'The most heavily loaded skid shoe or trailer group. Overload crushes the shoe, the skid beam or the structure node above it.'),
      check('Quay / skid-beam line load', Rmax / inp.shoeLen / 1e3, inp.quayCap, 'kN/m', `${N001} Sec. 10 (quay and skidway capacity); quay structural assessment`, 'Support reaction spread over its bearing length. Exceedance risks local failure of the quay deck or skid-beam foundation during the move.'),
    ];
    let axle = 0, tipDeg = 0;
    if (!skid) {
      axle = (W / 1e3 + inp.nAxle * inp.axleTare) / inp.nAxle;
      // Three-point hydraulic suspension: tipping lines form a triangle; the nearest line is about one third of the support width from the CoG.
      const half = Math.max(...sup.map((p) => Math.abs(p.y - yc)));
      tipDeg = Math.atan((2 * half / 3 - Math.abs(inp.ey)) / (inp.hCog + inp.hTrailer)) / DEG;
      checks.push(check('Trailer axle-line load', axle, inp.axleCap, 't', `${N001} Sec. 10 (trailers); trailer manufacturer load chart`, 'Average load per axle line including trailer self-weight. Above the chart value tyres and hydraulic cylinders are overloaded.'),
        check('Trailer stability angle', inp.tipReq, Math.max(tipDeg, 1e-6), 'deg (required / achieved)', `${N001} Sec. 10 (trailer stability)`, 'Tilt at which the weight passes over the tipping line of the hydraulic three-point support. A small angle means camber, wind or a hydraulic fault can topple the load.'));
    }
    // Barge hydrostatics at the starting draft (box barge), then the ballast needed at each transfer step.
    const hs = run(() => boxHydrostatics({ L: inp.Lb, B: inp.Bb, T: inp.T0, KG: 0.6 * inp.Db, rho }), 'Barge hydrostatics');
    need(inp.T0 < inp.Db, 'The starting draft exceeds the barge depth.');
    // Tanks are grouped fore and aft (four wing tanks each): two unknowns for the two targets draft and trim, which keeps the
    // bounded least-squares plan full rank. Heel from the transverse CoG offset is corrected by a port–starboard transfer.
    const capGroup = 0.8 * inp.Lb * inp.Bb * inp.Db / 2;
    const tanks = [{ name: 'Aft tanks', x: -inp.Lb / 4, y: 0, capacity: capGroup, volume: capGroup * inp.fill0 / 100, rho }, { name: 'Forward tanks', x: inp.Lb / 4, y: 0, capacity: capGroup, volume: capGroup * inp.fill0 / 100, rho }];
    const hydro = { displacement: hs.displacement, waterplaneArea: hs.waterplaneArea, lcf: 0, GMt: Math.max(hs.GMt, 0.5), GMl: Math.max(hs.GMl, 1), rho };
    const dist = inp.Ls + (inp.xFinal + inp.Lb / 2 - inp.Ls / 2); // travel from first contact with the stern to the final position
    need(inp.xFinal + inp.Lb / 2 >= inp.Ls / 2 - 1e-9 && inp.xFinal + inp.Ls / 2 <= inp.Lb / 2 + 1e-9, 'The final position leaves part of the structure off the barge deck.');
    const steps = [];
    const nT = inp.nStep, nTot = nT + Math.max(2, Math.round(nT * (dist - inp.Ls) / inp.Ls));
    let prev = null, qMax = 0, moved = 0, feasible = true, trimMax = 0;
    for (let k = 0; k <= nTot; k++) {
      const s = dist * k / nTot; // travel of the leading edge beyond the stern
      const f = Math.min(s / inp.Ls, 1), P = f * W;
      // Centroid of the part on the barge (uniform weight along Ls), measured from midship, positive forward.
      const xP = -inp.Lb / 2 + (f < 1 ? s / 2 : s - inp.Ls / 2);
      const t = s / inp.speed, tide = inp.tideRate * t;
      const current = { draft: inp.T0 + P / (rho * hs.waterplaneArea), heelDeg: 0, trimDeg: P * xP / (hs.displacement * hydro.GMl) / DEG };
      const plan = run(() => ballastPlan({ tanks, hydro, current, target: { draft: inp.T0 + tide, trimDeg: 0 }, mode: 'fill-discharge' }), 'Ballast plan');
      // Heel moment P·ey balanced by moving water between wing tanks whose centroids are Bb/2 apart.
      const vols = [...plan.tanks.map((q) => q.newVolume), P * Math.abs(inp.ey) / (rho * inp.Bb / 2)];
      let dV = 0; if (prev) for (let i = 0; i < vols.length; i++) dV += Math.abs(vols[i] - prev.vols[i]);
      const dt = prev ? t - prev.t : 0, rate = dt > 0 ? dV / dt : 0;
      qMax = Math.max(qMax, rate); moved += dV; feasible = feasible && plan.feasible; trimMax = Math.max(trimMax, Math.abs(current.trimDeg));
      steps.push({ k, s, f, P, xP, t, tide, aft: vols[0], fwd: vols[1], cross: vols[2], ballast: vols[0] + vols[1], rate, feasible: plan.feasible, errT: plan.error.draft ?? 0, trimFree: current.trimDeg });
      prev = { vols, t };
    }
    const dBallast = steps[steps.length - 1].ballast - steps[0].ballast;
    checks.push(check('Ballast pump capacity', qMax * inp.fPump, inp.pumpCap, 'm³/h', `${N001} Sec. 10 (load-out: ballast system capacity and redundancy)`, 'Pumping rate needed to follow the incoming load and the tide at the chosen speed, with margin for a pump failure. If the pumps cannot keep up the barge deck drops below the quay and the structure bridges the gap.'));
    checks.push(check('Ballast available in the tanks', feasible ? 0.5 : 1.5, 1, '- (plan feasible ≤ 1)', 'Bounded least-squares ballast plan (tank limits 0 … capacity)', 'Whether the assumed tanks can hold the deck level and the barge even at every step. If not, tanks run empty or full before the transfer is complete.'));
    const gov = worst(checks), ok = gov.uc <= 1;
    const hours = dist / inp.speed, days = (hours + 12) / 24;
    const speedMax = inp.speed * inp.pumpCap / Math.max(qMax * inp.fPump, 1e-9);
    return {
      kpis: [
        kpi('Fbreak', 'Break-out force', Fbreak / 1e3, 'kN', 'Force to start the structure moving: weight × (static friction or rolling resistance + gradient).', inp.pullCap / inp.fPull, 'max'),
        kpi('Frun', 'Running force', Frun / 1e3, 'kN', 'Force to keep the structure moving once friction has dropped to its kinetic value.'),
        kpi('Rmax', 'Maximum support reaction', Rmax / G / 1e3, 't', 'Largest reaction on a skid shoe or trailer group with the centre of gravity offset.', inp.shoeCap, 'max'),
        kpi('Rmin', 'Minimum support reaction', Rmin / G / 1e3, 't', 'Smallest reaction; it must stay positive or the support lifts off and the others take its share.'),
        kpi('qLine', 'Line load on the quay / skid beam', Rmax / inp.shoeLen / 1e3, 'kN/m', 'Maximum support reaction over its bearing length.', inp.quayCap, 'max'),
        kpi('axle', 'Trailer axle-line load', axle, 't', 'Average load per axle line; zero for a skidded load-out.'),
        kpi('tipDeg', 'Trailer stability angle', tipDeg, 'deg', 'Tilt to the tipping line of the hydraulic support triangle; zero for a skidded load-out.'),
        kpi('dBallast', 'Net ballast change over the load-out', dBallast, 'm³', 'Negative means water discharged. It equals the incoming weight less the tide effect on the required draft.'),
        kpi('moved', 'Total ballast pumped (all tanks)', moved, 'm³', 'Sum of the absolute volume changes of every tank over all steps: discharge plus the transfers that hold trim and heel.'),
        kpi('qMax', 'Peak pumping rate required', qMax, 'm³/h', 'Largest volume pumped per hour between two consecutive steps at the chosen load-out speed.', inp.pumpCap / inp.fPump, 'max'),
        kpi('speedMax', 'Largest load-out speed the pumps can follow', speedMax, 'm/h', 'Speed at which the required rate with margin equals the installed pump capacity.'),
        kpi('trimFree', 'Largest barge trim if no ballast were moved', trimMax, 'deg', 'Shows how much the ballast sequence has to correct; a large value means small pumping errors give a visible step at the quay edge.'),
        kpi('hours', 'Duration of the move', hours, 'h', 'Travel distance over load-out speed, excluding preparation and seafastening.'),
      ],
      checks,
      plots: [
        line('Ballast sequence', 'Travel of the leading edge past the stern', 'm', 'Total ballast on board', 'm³', [{ name: 'Ballast on board', x: steps.map((q) => round(q.s)), y: steps.map((q) => round(q.ballast)) }, { name: 'Load on barge (t)', x: steps.map((q) => round(q.s)), y: steps.map((q) => round(q.P / 1e3)) }],
          'Ballast on board falls as the structure comes onto the barge; a rising tide adds to the discharge. The slope is the pumping rate needed per metre of travel.'),
        line('Pumping rate along the move', 'Travel of the leading edge past the stern', 'm', 'Pumping rate', 'm³/h', [{ name: 'Required', x: steps.slice(1).map((q) => round(q.s)), y: steps.slice(1).map((q) => round(q.rate)) }, { name: 'Installed / margin', x: [0, round(dist)], y: [inp.pumpCap / inp.fPump, inp.pumpCap / inp.fPump], dash: true }],
          'Required pumping rate at each step against the usable capacity. Peaks occur while the weight transfers and the trim correction is largest.'),
        bar('Support reactions', sup.map((p) => p.name), 'Reaction', 't', [{ name: 'Reaction', y: R.map((v) => round(v / G / 1e3)) }, { name: 'Capacity', y: R.map(() => inp.shoeCap) }],
          'Share of the weight on each support with the CoG offset. Uneven bars show the effect of the offset; none may exceed the capacity.'),
      ],
      tables: [tbl('Load-transfer steps', [['k', 'Step'], ['s', 'Travel', 'm'], ['f', 'Share on barge', '%'], ['xP', 'Load centroid from midship', 'm'], ['t', 'Time', 'h'], ['aft', 'Aft tanks', 'm³'], ['fwd', 'Forward tanks', 'm³'], ['cross', 'Port–starboard transfer', 'm³'], ['rate', 'Pumping rate', 'm³/h'], ['trim', 'Trim without ballast', 'deg']],
        steps.map((q) => ({ k: q.k, s: round(q.s), f: round(q.f * 100), xP: round(q.xP), t: round(q.t, 3), aft: round(q.aft), fwd: round(q.fwd), cross: round(q.cross), rate: round(q.rate), trim: round(q.trimFree, 3) })),
        'One row per step of the move: how much of the structure is on the barge, the ballast that keeps the deck level with the quay and the barge even, and the pumping rate to get there from the previous step.')],
      notes: [`Break-out ${fmt(Fbreak / 1e3)} kN, running ${fmt(Frun / 1e3)} kN; ${skid ? `friction ${inp.muS} / ${inp.muK}` : `rolling resistance ${inp.rollRes}`} with ${inp.grade} % gradient.`, `Peak pumping ${fmt(qMax)} m³/h at ${inp.speed} m/h; tide ${inp.tideRate} m/h.`, `"${gov.name}" governs at ${fmt(gov.uc)}.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `The load-out is feasible: pull ${fmt(Fbreak / 1e3)} kN against ${fmt(inp.pullCap)} kN installed, peak support ${fmt(Rmax / G / 1e3)} t, peak pumping ${fmt(qMax)} m³/h against ${fmt(inp.pumpCap)} m³/h.`,
          `Hold the speed below ${fmt(speedMax)} m/h so the pumps keep up with margin, and re-run with the upper-bound friction (${fmt((skid ? inp.muS : inp.rollRes) * 1.3)}) and the largest tide rate of the planned window.`, `${N001} Sec. 10`)
          : adv('engineering', 'fail', `"${gov.name}" is not met (${fmt(gov.uc)}).`,
            gov.name.startsWith('Pull') ? `Install at least ${fmt(Fbreak * inp.fPull / 1e3)} kN of pull or reduce friction (PTFE pads, fresh lubrication).` : gov.name.startsWith('Ballast pump') ? `Reduce the load-out speed to ${fmt(speedMax)} m/h, add ${fmt(qMax * inp.fPump - inp.pumpCap)} m³/h of portable pumps, or plan the move around slack tide.` : gov.name.startsWith('Ballast avail') ? 'Start with more ballast on board (higher initial filling or deeper starting draft) so that the tanks do not run empty before the transfer is complete.' : `Add supports or spread the load: the governing reaction must come down by ${fmt((gov.uc - 1) / gov.uc * 100)} %.`,
            `${N001} Sec. 10`),
        adv('economic', 'info', `The move takes about ${fmt(hours)} h and pumps ${fmt(moved)} m³ of ballast; the barge and tugs are committed for roughly ${fmt(days)} day(s) including preparation.`,
          skid ? 'Skidding needs skid beams, shoes and strand jacks but little hire; trailers avoid the skid beams and cut preparation time at a hire cost per axle line per day — compare for this weight.' : `The spread is ${inp.nAxle} axle lines: each line removed by better weight control saves its hire for the whole mobilisation.`,
          'Load-out spread costs are time-related hire (barge, tugs, pumps, jacks or trailers): AACE RP 18R-97 estimate classes; ISO 19901-6 §11 planning',
          `Marine spread on station for the load-out: ${vesselCo2Text(days, 8)}`),
      ],
      publish: { 'install.loadoutPull': Fbreak, 'install.loadoutBallast': Math.abs(dBallast) },
      quantities: { vessel_days: { loadoutBarge: round(days) } },
    };
  },
  verify: [
    { name: 'Break-out and running pull: F = W·g·(μ + gradient)', source: 'Coulomb friction on a slope: 6000 t × 1.05 = 6300 t; break-out 6.3e6·9.80665·(0.15 + 0.005) = 9576.2 kN; running 6.3e6·9.80665·(0.08 + 0.005) = 5251.5 kN',
      inputs: {}, expect: [{ kpi: 'Fbreak', value: 6.3e6 * 9.80665 * 0.155 / 1e3, tol: 1e-9 }, { kpi: 'Frun', value: 6.3e6 * 9.80665 * 0.085 / 1e3, tol: 1e-9 }] },
    { name: 'Rigid-body support reactions with a CoG offset', source: 'Statics: 8 supports at x = ±5, ±15 m and y = ±10 m; W = 6300 t, ex = 1 m, ey = 0.5 m → R = W/8 + W·ex·x/Σx² + W·ey·y/Σy² = 787.5 + 6300·15/1000 + 6300·0.5·10/800 = 921.375 t at (15, 10)',
      inputs: {}, expect: [{ kpi: 'Rmax', value: 787.5 + 94.5 + 39.375, tol: 1e-9 }, { kpi: 'Rmin', value: 787.5 - 94.5 - 39.375, tol: 1e-9 }] },
    { name: 'Net ballast discharge equals the weight taken on at constant draft', source: 'Archimedes: with no tide the barge must discharge the mass it receives: 6.3e6 kg / 1025 kg/m³ = 6146.3 m³',
      inputs: { tideRate: 0 }, expect: [{ kpi: 'dBallast', value: -6.3e6 / 1025, tol: 0.002 }] },
  ],
};
