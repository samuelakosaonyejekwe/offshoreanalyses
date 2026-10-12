// Mooring line tensions and bollard loads of a ship alongside a continuous quay under wind and current: static
// three-degree-of-freedom equilibrium with tension-only elastic lines and compression-only fenders, the compass
// sweep and the limiting wind speed.
import { berthEquilibrium, genericCoefficientCurves } from '../../engine/lines/berthmooring.js';
import { need, run, r, f, check, eng, eco, maxWhere, CONFIRM, ECON_BASIS, CO2_STEEL } from './_common.js';

const LINES_DEFAULT = [
  { name: 'Head lines', bol: 'B1', fx: 112, fy: -8, fz: 19, bx: 150, by: -3, bz: 4, n: 3, MBL: 1300, material: 'synthetic', EA: 40000, pre: 80 },
  { name: 'Fwd breast', bol: 'B2', fx: 100, fy: -14, fz: 19, bx: 110, by: -3, bz: 4, n: 3, MBL: 1300, material: 'synthetic', EA: 40000, pre: 80 },
  { name: 'Fwd spring', bol: 'B3', fx: 80, fy: -18.5, fz: 19, bx: 40, by: -3, bz: 4, n: 3, MBL: 1300, material: 'synthetic', EA: 40000, pre: 80 },
  { name: 'Aft spring', bol: 'B4', fx: -80, fy: -18.5, fz: 19, bx: -40, by: -3, bz: 4, n: 3, MBL: 1300, material: 'synthetic', EA: 40000, pre: 80 },
  { name: 'Aft breast', bol: 'B5', fx: -100, fy: -14, fz: 19, bx: -110, by: -3, bz: 4, n: 3, MBL: 1300, material: 'synthetic', EA: 40000, pre: 80 },
  { name: 'Stern lines', bol: 'B6', fx: -112, fy: -8, fz: 19, bx: -150, by: -3, bz: 4, n: 3, MBL: 1300, material: 'synthetic', EA: 40000, pre: 80 },
];
const FENDERS_DEFAULT = [{ name: 'F1', x: -75, RR: 1700, defl: 0.9, h: 1.3 }, { name: 'F2', x: -25, RR: 1700, defl: 0.9, h: 1.3 }, { name: 'F3', x: 25, RR: 1700, defl: 0.9, h: 1.3 }, { name: 'F4', x: 75, RR: 1700, defl: 0.9, h: 1.3 }];

function build(inp) {
  const vessel = { LBP: inp.LBP, beam: inp.beam, draft: inp.draft, AT: inp.AT, AL: inp.AL, type: inp.vtype, side: 'starboard', parallelBody: [-0.35 * inp.LBP, 0.35 * inp.LBP] };
  const rows = inp.lines || [];
  need(rows.length >= 1, 'Give at least one mooring line.');
  const lines = rows.map((q, i) => {
    const n = Math.max(1, Math.round(q.n || 1)), mat = String(q.material || 'synthetic').trim().toLowerCase();
    need(['steel', 'synthetic', 'polyamide'].includes(mat), `Line ${i + 1}: material must be steel, synthetic or polyamide.`);
    need(q.MBL > 0 && q.EA > 0, `Line ${i + 1}: breaking strength and axial stiffness must be positive.`);
    need((q.pre || 0) < q.MBL, `Line ${i + 1}: the pretension must be below the breaking strength.`);
    // n parallel lines between the same fairlead and bollard act as one line of n-fold strength and stiffness (kN → N)
    return { name: String(q.name || `Line ${i + 1}`), fairlead: [Number(q.fx), Number(q.fy), Number(q.fz)], bollard: [Number(q.bx), Number(q.by), Number(q.bz)], MBL: n * q.MBL * 1e3, material: mat, pretension: n * (q.pre || 0) * 1e3, EA: n * q.EA * 1e3, bollardId: String(q.bol || `B${i + 1}`), count: n };
  });
  const names = new Set(lines.map((l) => l.name));
  need(names.size === lines.length, 'Mooring line names must be unique.');
  const fenders = (inp.fenders || []).map((q, i) => {
    need(q.RR > 0 && q.defl > 0 && q.h > q.defl, `Fender ${i + 1}: rated reaction and deflection must be positive and the height must exceed the rated deflection.`);
    return { name: String(q.name || `Fender ${i + 1}`), position: [Number(q.x), -q.h], normal: [0, 1], height: Number(q.h), curve: { deflection: [0, q.defl], reaction: [0, q.RR * 1e3] }, ratedReaction: q.RR * 1e3 };
  });
  const ids = [...new Set(lines.map((l) => l.bollardId))];
  const windCurves = genericCoefficientCurves({ type: inp.vtype, kind: 'wind', ...(inp.CyW > 0 ? { Cy0: inp.CyW } : {}), ...(inp.CxW > 0 ? { Cx0: inp.CxW } : {}), ...(inp.CxyW >= 0 ? { Cxy0: inp.CxyW } : {}) });
  const curCurves = genericCoefficientCurves({ type: inp.vtype, kind: 'current', shallowFactor: inp.shallow });
  const model = { vessel, pose0: { x: 0, y: inp.beam / 2, headingDeg: 0 }, lines, fenders, reference: { waterLevel: 0, draft: inp.draft }, condition: { waterLevel: inp.tide, draft: inp.draftCond > 0 ? inp.draftCond : inp.draft },
    limits: { steel: inp.limSteel, synthetic: inp.limSyn, polyamide: inp.limPa }, bollards: ids.map((id) => ({ id, SWL: inp.swl * 1e3 })) };
  const env = (Uw, from) => ({ wind: { speed: Uw, fromDeg: from, curves: windCurves, rho: inp.rhoAir }, current: { speed: inp.Uc, fromDeg: inp.curFrom, curves: curCurves, rho: inp.rhoW }, extra: { surge: inp.extraSurge * 1e3, sway: inp.extraSway * 1e3, yaw: 0 } });
  return { model, env };
}

export default {
  id: 's05-mooring-bollard', suite: 5, title: 'Mooring forces and bollard loads',
  summary: 'Line tensions, bollard and hook loads, fender reactions and vessel excursion for a ship moored alongside under wind, current and passing-vessel loads, with the limiting wind speed for the installed equipment.',
  covers: [3], fidelity: 'reduced-order',
  standards: ['OCIMF Mooring Equipment Guidelines, 4th ed. (2018)', 'BS 6349-4:2014 §7 (mooring)', 'PIANC WG 153 (2016) marine oil and petrochemical terminals', 'BS 6349-1-2:2016 (actions)'],
  models: ['Static equilibrium in surge, sway and yaw', 'Line tension equilibrium with nonlinear (tension-only) elastic lines', 'Compression-only fender springs', 'Wind and current force coefficients on the vessel', 'Vector resolution of line tension at bollards'],
  limits: 'Static three-degree-of-freedom model: no wave-induced motion, no dynamic amplification of gusts, long waves or passing ships (their peak forces can be entered as additional loads). The built-in wind and current coefficient curves are generic shapes, not OCIMF data: enter vessel-specific coefficients for design. Dynamic mooring analysis is outside the module.',
  publishes: ['wharf.bollardLoad', 'wharf.bollardLoadBent', 'wharf.mooringLineMax'],
  searchTerms: 'mooring analysis bollard load OCIMF berth line tension',
  inputs: [
    { key: 'vtype', label: 'Vessel type (generic coefficient shape)', type: 'select', default: 'container', group: 'Vessel', role: 'operation', options: [{ value: 'tanker', label: 'Tanker' }, { value: 'lng', label: 'LNG carrier' }, { value: 'container', label: 'Container ship' }, { value: 'bulk', label: 'Bulk carrier' }] },
    { key: 'LBP', label: 'Length between perpendiculars', unit: 'm', type: 'number', default: 240, min: 30, max: 420, group: 'Vessel', role: 'geometry' },
    { key: 'beam', label: 'Beam', unit: 'm', type: 'number', default: 37, min: 5, max: 70, group: 'Vessel', role: 'geometry' },
    { key: 'draft', label: 'Draft when the lines were tended', unit: 'm', type: 'number', default: 12, min: 1, max: 25, group: 'Vessel', role: 'ic' },
    { key: 'draftCond', label: 'Draft in the condition analysed (0 = same)', unit: 'm', type: 'number', default: 0, min: 0, max: 25, group: 'Vessel', role: 'ic', help: 'Cargo operations change the draft and therefore the fairlead height and line tensions.' },
    { key: 'AL', label: 'Lateral windage area', unit: 'm²', type: 'number', default: 6500, min: 50, max: 25000, group: 'Vessel', role: 'geometry', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Broadside area above water including deck cargo.' },
    { key: 'AT', label: 'Transverse (head-on) windage area', unit: 'm²', type: 'number', default: 1100, min: 20, max: 5000, group: 'Vessel', role: 'geometry' },
    { key: 'lines', label: 'Mooring lines', type: 'table', group: 'Mooring arrangement', role: 'bc',
      columns: [{ key: 'name', label: 'Line group', type: 'text' }, { key: 'bol', label: 'Bollard', type: 'text' }, { key: 'fx', label: 'Fairlead x (fwd +)', unit: 'm', type: 'number' }, { key: 'fy', label: 'Fairlead y (berth side −)', unit: 'm', type: 'number' }, { key: 'fz', label: 'Fairlead above keel', unit: 'm', type: 'number' },
        { key: 'bx', label: 'Bollard x', unit: 'm', type: 'number' }, { key: 'by', label: 'Bollard y (behind fender face −)', unit: 'm', type: 'number' }, { key: 'bz', label: 'Bollard level', unit: 'm', type: 'number' }, { key: 'n', label: 'Lines', type: 'number' }, { key: 'MBL', label: 'MBL per line', unit: 'kN', type: 'number' }, { key: 'material', label: 'steel / synthetic / polyamide', type: 'text' }, { key: 'EA', label: 'EA per line', unit: 'kN', type: 'number' }, { key: 'pre', label: 'Pretension per line', unit: 'kN', type: 'number' }],
      default: LINES_DEFAULT, help: 'Ship axes: origin amidships on the centreline, x forward, y to port; the starboard side lies against the fenders. Berth axes: x along the quay, y = 0 at the fender face, bollards at negative y; levels from the water-level datum.' },
    { key: 'fenders', label: 'Fenders in contact', type: 'table', group: 'Mooring arrangement', role: 'bc', columns: [{ key: 'name', label: 'Fender', type: 'text' }, { key: 'x', label: 'x along the quay', unit: 'm', type: 'number' }, { key: 'RR', label: 'Rated reaction', unit: 'kN', type: 'number' }, { key: 'defl', label: 'Rated deflection', unit: 'm', type: 'number' }, { key: 'h', label: 'Height', unit: 'm', type: 'number' }],
      default: FENDERS_DEFAULT, help: 'Each fender is a linear spring to its rated point (secant stiffness); the berthing module holds the full curve.' },
    { key: 'swl', label: 'Bollard safe working load', unit: 'kN', type: 'number', default: 2000, min: 50, max: 5000, group: 'Mooring arrangement', role: 'criteria', design: true, help: 'Rated load of the bollard and its anchorage; several lines may share one bollard.' },
    { key: 'nShare', label: 'Bents sharing one bollard load', unit: '-', type: 'number', default: 3, min: 1, max: 12, group: 'Mooring arrangement', role: 'bc' },
    { key: 'Uw', label: 'Wind speed (30 s gust at 10 m)', unit: 'm/s', type: 'number', default: 22, min: 0, max: 60, group: 'Environment', role: 'environment', uncertain: { dist: 'gumbel', cov: 0.12 }, help: 'Operating limit wind for remaining alongside; OCIMF uses a 30-second gust.' },
    { key: 'windFrom', label: 'Wind from (anticlockwise from the quay x axis)', unit: 'deg', type: 'number', default: 270, min: 0, max: 360, group: 'Environment', role: 'environment', help: '270° blows the vessel off the quay and loads the breast lines; 90° presses it on the fenders; 0° and 180° load the springs.' },
    { key: 'Uc', label: 'Current speed', unit: 'm/s', type: 'number', default: 0.5, min: 0, max: 4, group: 'Environment', role: 'environment', bind: 'site.Uc1' },
    { key: 'curFrom', label: 'Current from', unit: 'deg', type: 'number', default: 180, min: 0, max: 360, group: 'Environment', role: 'environment' },
    { key: 'shallow', label: 'Shallow-water factor on lateral current force', unit: '-', type: 'number', default: 2, min: 1, max: 6, group: 'Environment', role: 'environment', help: 'Lateral current load grows several-fold as depth/draft falls towards 1.1; take from OCIMF charts for the actual ratio.' },
    { key: 'tide', label: 'Water level in the condition analysed', unit: 'm', type: 'number', default: 0, min: -6, max: 8, group: 'Environment', role: 'bc', help: 'Relative to the level at which the lines were tended; a rising tide slackens and a falling tide tightens untended lines.' },
    { key: 'extraSurge', label: 'Additional surge force (passing vessel)', unit: 'kN', type: 'number', default: 0, min: -5000, max: 5000, group: 'Environment', role: 'bc', help: 'Peak value from a passing-ship calculation, in ship axes.' },
    { key: 'extraSway', label: 'Additional sway force (off the berth −)', unit: 'kN', type: 'number', default: 0, min: -5000, max: 5000, group: 'Environment', role: 'bc' },
    { key: 'rhoAir', label: 'Air density', unit: 'kg/m³', type: 'number', default: 1.28, min: 1.1, max: 1.4, group: 'Environment', role: 'environment', help: '1.28 kg/m³ is the value used with the OCIMF wind coefficients.' },
    { key: 'rhoW', label: 'Water density', unit: 'kg/m³', type: 'number', default: 1025, min: 995, max: 1035, group: 'Environment', role: 'environment', bind: 'site.rhoW' },
    { key: 'CyW', label: 'Lateral wind coefficient (0 = generic)', unit: '-', type: 'number', default: 0, min: 0, max: 2, group: 'Coefficients', role: 'criteria', help: 'Peak of the lateral coefficient curve on the lateral windage area. Enter the vessel-specific value (OCIMF, wind-tunnel); the generic shapes are not OCIMF data.' },
    { key: 'CxW', label: 'Longitudinal wind coefficient (0 = generic)', unit: '-', type: 'number', default: 0, min: 0, max: 2, group: 'Coefficients', role: 'criteria' },
    { key: 'CxyW', label: 'Wind yaw coefficient (−1 = generic)', unit: '-', type: 'number', default: -1, min: -1, max: 0.5, group: 'Coefficients', role: 'criteria' },
    { key: 'limSteel', label: 'Line limit, steel wire (fraction of MBL)', unit: '-', type: 'number', default: 0.55, min: 0.2, max: 1, group: 'Acceptance', role: 'criteria', help: `OCIMF MEG. ${CONFIRM}` },
    { key: 'limSyn', label: 'Line limit, synthetic (fraction of MBL)', unit: '-', type: 'number', default: 0.5, min: 0.2, max: 1, group: 'Acceptance', role: 'criteria', help: `OCIMF MEG. ${CONFIRM}` },
    { key: 'limPa', label: 'Line limit, polyamide (fraction of MBL)', unit: '-', type: 'number', default: 0.45, min: 0.2, max: 1, group: 'Acceptance', role: 'criteria', help: `OCIMF MEG. ${CONFIRM}` },
    { key: 'surgeAllow', label: 'Allowable surge alongside', unit: 'm', type: 'number', default: 1.0, min: 0.1, max: 10, group: 'Acceptance', role: 'criteria', help: 'PIANC WG24 (1995) motion criteria for cargo handling: about 1 m for container cranes at 100 % efficiency, 2–3 m for loading arms and conveyors.' },
    { key: 'sweep', label: 'Sweep the wind around the compass', type: 'bool', default: true, group: 'Acceptance', role: 'control', help: 'Runs the wind from every 30° to find the governing direction.' },
  ],
  compute(inp) {
    const { model, env } = build(inp);
    const solve = (Uw, from) => berthEquilibrium({ ...model, ...env(Uw, from) });
    const res = run(() => solve(inp.Uw, inp.windFrom), 'Mooring equilibrium');
    const count = Object.fromEntries(model.lines.map((l) => [l.name, l.count]));
    const mbl = Object.fromEntries(model.lines.map((l) => [l.name, l.MBL]));
    const lines = res.lines.map((l) => ({ ...l, MBL: mbl[l.name], per: l.tension / count[l.name] }));
    const gov = lines.reduce((a, l) => (l.uc > a.uc ? l : a), lines[0]);
    const bGov = res.bollards.reduce((a, b) => (b.load > a.load ? b : a), res.bollards[0]);
    const fGov = res.fenders.length ? res.fenders.reduce((a, q) => (q.reaction > a.reaction ? q : a), res.fenders[0]) : null;
    const envF = res.environment.reduce((a, e) => ({ Fx: a.Fx + e.Fx, Fy: a.Fy + e.Fy }), { Fx: 0, Fy: 0 });
    const ucOf = (q) => Math.max(q.checks.maxLineUc, q.checks.maxBollardUc ?? 0);
    const dirs = [], ucDir = [];
    let govDir = inp.windFrom, govUc = ucOf(res);
    if (inp.sweep) {
      for (let d = 0; d < 360; d += 30) {
        let u = 99; try { u = ucOf(solve(inp.Uw, d)); } catch { u = 99; }
        dirs.push(d); ucDir.push(Math.min(u, 3));
        if (u > govUc) { govUc = u; govDir = d; }
      }
    }
    // limiting wind speed from the governing direction: utilisation of lines and bollards reaches 1
    const okU = (U) => { try { return ucOf(solve(U, govDir)) <= 1; } catch { return false; } };
    const Ulim = maxWhere(okU, 0, 60, 10);
    const surge = Math.abs(res.excursion.surge), sway = Math.abs(res.excursion.sway);
    const ok = govUc <= 1;
    const swlNeed = bGov.load / 1e3;
    return {
      kpis: [
        { key: 'Fx', label: 'Longitudinal environmental force', value: r(envF.Fx / 1e3), unit: 'kN', explain: 'Wind plus current force along the vessel (bow positive), resisted by the spring lines.' },
        { key: 'Fy', label: 'Lateral environmental force', value: r(Math.abs(envF.Fy) / 1e3), unit: 'kN', explain: 'Wind plus current force across the vessel, resisted by the breast lines when off the berth and by the fenders when onto it.' },
        { key: 'Tmax', label: 'Largest line tension (per line)', value: r(gov.per / 1e3), unit: 'kN', explain: `Tension in each line of the most loaded group (${gov.name}). Compare with the working limit, a fraction of the minimum breaking load.` },
        { key: 'lineUc', label: 'Line utilisation', value: r(gov.uc), unit: '-', limit: 1, sense: 'max', explain: 'Tension over the permitted fraction of the breaking load. Above 1 the line is overloaded: a parted line whips and sheds its load onto the others.' },
        { key: 'bollard', label: 'Largest bollard load', value: r(bGov.load / 1e3), unit: 'kN', limit: inp.swl, sense: 'max', explain: `Vector sum of the lines on bollard ${bGov.id}, against its safe working load.` },
        { key: 'bollardH', label: 'Horizontal component at that bollard', value: r(bGov.horizontal / 1e3), unit: 'kN', explain: 'Horizontal pull passed into the deck and piles; the vertical component lifts the bollard foundation.' },
        { key: 'bollardV', label: 'Uplift at that bollard', value: r(bGov.vertical / 1e3), unit: 'kN', explain: 'Vertical component of the line pull (up positive). Steep lines to a high-freeboard vessel give large uplift on the holding-down bolts.' },
        { key: 'fender', label: 'Largest fender reaction', value: r((fGov?.reaction ?? 0) / 1e3), unit: 'kN', explain: 'Reaction of the most compressed fender when the vessel is pressed onto the berth; zero when all fenders are clear.' },
        { key: 'surge', label: 'Surge excursion', value: r(surge), unit: 'm', limit: inp.surgeAllow, sense: 'max', explain: 'Movement of the vessel along the quay from the tended position. Cargo-handling equipment tolerates only limited surge.' },
        { key: 'sway', label: 'Sway excursion', value: r(sway), unit: 'm', explain: 'Movement off (or onto) the quay from the tended position.' },
        { key: 'govUc', label: 'Governing utilisation over all wind directions', value: r(govUc), unit: '-', limit: 1, sense: 'max', explain: `Largest line or bollard utilisation found in the compass sweep (wind from ${f(govDir)}°).` },
        { key: 'Ulim', label: 'Maximum permissible wind speed alongside', value: r(Ulim), unit: 'm/s', explain: 'Wind speed from the governing direction at which a line or bollard reaches its limit, with the stated current. The vessel should leave the berth, or lines be doubled, before this is forecast. Reported at 60 m/s when not reached.' },
      ],
      checks: [
        check(`Line tension (${gov.name})`, gov.per / 1e3, gov.limitFraction * gov.MBL / count[gov.name] / 1e3, 'kN', 'OCIMF MEG4 (2018) §1 and §3 (line working load limits)', 'Guards against a line parting. A parted line is a snap-back hazard and overloads its neighbours in sequence.'),
        check(`Bollard load (${bGov.id})`, bGov.load / 1e3, inp.swl, 'kN', 'BS 6349-4:2014 §7; PIANC WG 153 (2016)', 'Guards against pulling the bollard or its anchorage out of the deck.'),
        check('Surge excursion', surge, inp.surgeAllow, 'm', 'PIANC WG24 (1995) Criteria for Movements of Moored Ships in Harbours, Table 4.1', 'Guards against interrupting cargo handling: beyond the limit cranes or loading arms must stop or disconnect.'),
      ],
      plots: [
        { type: 'bar', title: 'Line group tensions', categories: lines.map((l) => l.name), series: [{ name: 'Tension per line (kN)', y: lines.map((l) => r(l.per / 1e3)) }, { name: 'Limit per line (kN)', y: lines.map((l) => r(l.limitFraction * l.MBL / count[l.name] / 1e3)) }],
          x: { label: 'Line group', unit: '' }, y: { label: 'Tension', unit: 'kN' }, explain: 'Tension in each line against its working limit. A good arrangement shares the load evenly between lines of the same service; one bar far above the rest marks a short, stiff line.' },
        ...(inp.sweep ? [{ type: 'polar', title: 'Utilisation against wind direction', series: [{ name: 'Line / bollard utilisation', thetaDeg: [...dirs, 360], r: [...ucDir, ucDir[0]] }, { name: 'Limit', thetaDeg: [...dirs, 360], r: [...dirs.map(() => 1), 1] }],
          explain: 'Largest line or bollard utilisation for the wind from each direction (capped at 3). Lobes outside the unit circle are directions in which the vessel cannot stay at this wind speed.' }] : []),
        { type: 'scatter', title: 'Mooring plan', x: { label: 'Along the quay', unit: 'm' }, y: { label: 'Off the quay', unit: 'm' },
          series: [{ name: 'Fairleads', x: model.lines.map((l) => r(l.fairlead[0])), y: model.lines.map((l) => r(inp.beam / 2 + l.fairlead[1])) }, { name: 'Bollards', x: model.lines.map((l) => l.bollard[0]), y: model.lines.map((l) => l.bollard[1]) }, { name: 'Fenders', x: model.fenders.map((q) => q.position[0]), y: model.fenders.map(() => 0) }],
          explain: 'Plan of fairleads, bollards and fenders to check the entered geometry: breast lines should be near-perpendicular to the hull and springs near-parallel.' },
      ],
      tables: [
        { title: 'Mooring lines', columns: [{ key: 'name', label: 'Group' }, { key: 'n', label: 'Lines' }, { key: 'T', label: 'Tension per line', unit: 'kN' }, { key: 'pct', label: '% MBL' }, { key: 'uc', label: 'Utilisation' }, { key: 'vert', label: 'Vertical angle', unit: 'deg' }],
          rows: lines.map((l) => ({ name: l.name, n: count[l.name], T: r(l.per / 1e3), pct: r(l.utilisation * 100, 3), uc: r(l.uc, 3), vert: r(l.verticalAngle * 180 / Math.PI, 3) })), explain: 'One row per line group. Lines steeper than about 25° lose efficiency and load the bollard in uplift.' },
        { title: 'Bollard loads', columns: [{ key: 'id', label: 'Bollard' }, { key: 'load', label: 'Resultant', unit: 'kN' }, { key: 'H', label: 'Horizontal', unit: 'kN' }, { key: 'V', label: 'Vertical', unit: 'kN' }, { key: 'uc', label: 'Utilisation' }],
          rows: res.bollards.map((b) => ({ id: b.id, load: r(b.load / 1e3), H: r(b.horizontal / 1e3), V: r(b.vertical / 1e3), uc: r(b.load / (inp.swl * 1e3), 3) })), explain: 'Vector sum of the line pulls on each bollard for the wind direction entered. Use the horizontal value as the mooring action on the deck and piles.' },
      ],
      notes: [`${gov.name} governs at ${f(gov.utilisation * 100)} % of MBL; bollard ${bGov.id} carries ${f(bGov.load / 1e3)} kN. Equilibrium converged in ${res.iterations} iteration(s).`,
        `Generic ${inp.vtype} coefficient shapes${inp.CyW > 0 ? ' with the lateral wind coefficient overridden' : ''}; the static model has no gust or wave dynamics.`],
      advice: [
        ok ? eng('ok', `Lines and bollards are within limits for ${f(inp.Uw)} m/s wind: governing utilisation ${f(govUc)} with the wind from ${f(govDir)}°; the arrangement reaches its limit at ${f(Ulim)} m/s.`,
          `Set the berth operating limit below ${f(Ulim)} m/s with a margin for gust dynamics (the model is static), and confirm the wind and current coefficients for the actual vessel.`, 'OCIMF MEG4 (2018) §3; BS 6349-4:2014 §7',
          'Keeping within the installed bollard rating avoids replacing bollards and strengthening their anchorage, which would mean breaking out and recasting deck concrete.')
          : eng('fail', `Utilisation ${f(govUc)} with the wind from ${f(govDir)}° at ${f(inp.Uw)} m/s: ${gov.uc > 1 ? `${gov.name} at ${f(gov.utilisation * 100)} % of MBL` : `bollard ${bGov.id} at ${f(bGov.load / 1e3)} kN against ${f(inp.swl)} kN`}.`,
            `Limit the wind alongside to ${f(Ulim)} m/s, or add lines in the governing service (breast lines for off-berth wind, springs for along-berth), or fit bollards of at least ${f(Math.ceil(swlNeed / 50) * 50)} kN.`, 'OCIMF MEG4 (2018) §3; BS 6349-4:2014 §7',
            `Upgrading a bollard anchorage adds roughly ${f(0.004 * Math.max(swlNeed - inp.swl, 0))} t of anchor steel per bollard (${f(0.004 * Math.max(swlNeed - inp.swl, 0) * CO2_STEEL)} t CO₂e, screening factor); an operating limit achieves the same safety with no material.`),
        eco('info', `The berth can hold this vessel up to ${f(Ulim)} m/s; above that the vessel must sail or take extra lines. Surge is ${f(surge)} m against ${f(inp.surgeAllow)} m for cargo handling.`,
          'Convert the limiting wind speed into expected weather downtime with the site wind statistics, and compare the downtime hours with the cost of stronger bollards or shore-tension units.', `${ECON_BASIS}; PIANC WG24 (1995) operability criteria`),
      ],
      publish: { 'wharf.bollardLoad': bGov.load, 'wharf.bollardLoadBent': bGov.horizontal / inp.nShare, 'wharf.mooringLineMax': gov.per },
    };
  },
  verify: [
    // Beam wind off the berth: Fy = ½·1.28·1.0·20²·3000 = 768 kN; two perpendicular horizontal breast lines share it: T = 384 kN each.
    { name: 'Two horizontal breast lines under beam wind: T = F/2', source: 'Statics; wind force ½ ρ Cy V² AL with ρ = 1.28 kg/m³',
      inputs: { vtype: 'tanker', LBP: 200, beam: 30, draft: 10, AL: 3000, AT: 500, Uw: 20, windFrom: 270, Uc: 0, CyW: 1, CxW: 0, CxyW: 0, sweep: false, fenders: [], swl: 5000, surgeAllow: 10,
        lines: [{ name: 'B1', bol: 'A', fx: 50, fy: -15, fz: 15, bx: 50, by: -30, bz: 5, n: 1, MBL: 5000, material: 'steel', EA: 1e6, pre: 0 }, { name: 'B2', bol: 'B', fx: -50, fy: -15, fz: 15, bx: -50, by: -30, bz: 5, n: 1, MBL: 5000, material: 'steel', EA: 1e6, pre: 0 },
          { name: 'S1', bol: 'C', fx: 60, fy: -15, fz: 15, bx: 0, by: -1, bz: 5, n: 1, MBL: 5000, material: 'steel', EA: 1e6, pre: 0 }, { name: 'S2', bol: 'D', fx: -60, fy: -15, fz: 15, bx: 0, by: -1, bz: 5, n: 1, MBL: 5000, material: 'steel', EA: 1e6, pre: 0 }] },
      expect: [{ kpi: 'Fy', value: 768, tol: 1e-6 }, { kpi: 'Tmax', value: 384, tol: 0.01 }] },
    // Bollards 10 m below the fairleads over 30 m horizontal: cos φ = 30/√(30² + 10²) = 0.94868; T = 384/0.94868 = 404.8 kN;
    // horizontal component on each bollard 384 kN, uplift 384·tan φ = 128 kN.
    { name: 'Inclined breast lines: T = F/(2 cos φ), bollard uplift F/2·tan φ', source: 'Resolution of line tension for a known fairlead and bollard geometry (statics)',
      inputs: { vtype: 'tanker', LBP: 200, beam: 30, draft: 10, AL: 3000, AT: 500, Uw: 20, windFrom: 270, Uc: 0, CyW: 1, CxW: 0, CxyW: 0, sweep: false, fenders: [], swl: 5000, surgeAllow: 10,
        lines: [{ name: 'B1', bol: 'A', fx: 50, fy: -15, fz: 15, bx: 50, by: -30, bz: -5, n: 1, MBL: 5000, material: 'steel', EA: 1e6, pre: 0 }, { name: 'B2', bol: 'B', fx: -50, fy: -15, fz: 15, bx: -50, by: -30, bz: -5, n: 1, MBL: 5000, material: 'steel', EA: 1e6, pre: 0 },
          { name: 'S1', bol: 'C', fx: 60, fy: -15, fz: 15, bx: 0, by: -1, bz: 5, n: 1, MBL: 5000, material: 'steel', EA: 1e6, pre: 0 }, { name: 'S2', bol: 'D', fx: -60, fy: -15, fz: 15, bx: 0, by: -1, bz: 5, n: 1, MBL: 5000, material: 'steel', EA: 1e6, pre: 0 }] },
      expect: [{ kpi: 'Tmax', value: 404.77, tol: 0.01 }, { kpi: 'bollardH', value: 384, tol: 0.01 }, { kpi: 'bollardV', value: 128, tol: 0.02 }] },
  ],
};
