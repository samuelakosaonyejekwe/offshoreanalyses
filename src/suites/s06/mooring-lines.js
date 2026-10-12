// Mooring line tensions, bollard / hook reactions and fender loads of a ship moored at dolphins under wind and
// current: static three-degree-of-freedom equilibrium with tension-only lines and compression-only fenders.
import { berthEquilibrium, runBerthCases, genericCoefficientCurves } from '../../engine/lines/berthmooring.js';
import { inextensibleCatenaryShape } from '../../engine/lines/catenary.js';
import { need, run, fmt, round, adv, check, DEG, CONFIRM, NOT_HIFI, ECON_BASIS, vesselCo2Text } from './_common.js';

const LINES_DEFAULT = [
  { name: 'Head', dol: 'MD1', fx: 108, fy: -6, fz: 20, bx: 150, by: -35, bz: 6, n: 2, MBL: 900, material: 'steel', EA: 70000, pre: 80, tailL: 11, tailEA: 6000 },
  { name: 'Fwd breast', dol: 'MD2', fx: 92, fy: -13, fz: 20, bx: 95, by: -45, bz: 6, n: 2, MBL: 900, material: 'steel', EA: 70000, pre: 80, tailL: 11, tailEA: 6000 },
  { name: 'Fwd spring', dol: 'BD1', fx: 75, fy: -16.1, fz: 20, bx: 40, by: -3, bz: 6, n: 2, MBL: 900, material: 'steel', EA: 70000, pre: 80, tailL: 11, tailEA: 6000 },
  { name: 'Aft spring', dol: 'BD2', fx: -75, fy: -16.1, fz: 20, bx: -40, by: -3, bz: 6, n: 2, MBL: 900, material: 'steel', EA: 70000, pre: 80, tailL: 11, tailEA: 6000 },
  { name: 'Aft breast', dol: 'MD3', fx: -92, fy: -13, fz: 20, bx: -95, by: -45, bz: 6, n: 2, MBL: 900, material: 'steel', EA: 70000, pre: 80, tailL: 11, tailEA: 6000 },
  { name: 'Stern', dol: 'MD4', fx: -108, fy: -6, fz: 20, bx: -150, by: -35, bz: 6, n: 2, MBL: 900, material: 'steel', EA: 70000, pre: 80, tailL: 11, tailEA: 6000 },
];
const FENDERS_DEFAULT = [{ name: 'BD1', x: 40, RR: 2500, defl: 0.9, h: 1.3 }, { name: 'BD2', x: -40, RR: 2500, defl: 0.9, h: 1.3 }];

function buildModel(inp) {
  const vessel = { LBP: inp.LBP, beam: inp.beam, draft: inp.draft, AT: inp.AT, AL: inp.AL, type: inp.vtype, side: 'starboard', parallelBody: [-inp.pbFrac * inp.LBP, inp.pbFrac * inp.LBP] };
  const rows = inp.lines || [];
  need(rows.length >= 1, 'Give at least one mooring line.');
  const lines = rows.map((r, i) => {
    const n = Math.max(1, Math.round(r.n || 1));
    const mat = String(r.material || 'synthetic').trim().toLowerCase();
    need(['steel', 'synthetic', 'polyamide'].includes(mat), `Line ${i + 1}: material must be steel, synthetic or polyamide.`);
    need(r.MBL > 0 && r.EA > 0, `Line ${i + 1}: breaking strength and axial stiffness must be positive.`);
    need((r.pre || 0) < r.MBL, `Line ${i + 1}: the pretension must be below the breaking strength.`);
    // n parallel lines on the same fairlead and hook act as one line of n-fold strength and stiffness (SI: kN → N)
    const ln = { name: r.name || `Line ${i + 1}`, fairlead: [r.fx, r.fy, r.fz], bollard: [r.bx, r.by, r.bz], MBL: n * r.MBL * 1e3, material: mat, pretension: n * (r.pre || 0) * 1e3, EA: n * r.EA * 1e3, bollardId: String(r.dol || `B${i + 1}`), count: n };
    if (r.tailL > 0 && r.tailEA > 0) ln.tail = { L: r.tailL, EA: n * r.tailEA * 1e3 };
    return ln;
  });
  const fenders = (inp.fenders || []).map((f, i) => {
    need(f.RR > 0 && f.defl > 0 && f.h > f.defl, `Fender ${i + 1}: rated reaction and deflection must be positive and the height must exceed the rated deflection.`);
    return { name: f.name || `Fender ${i + 1}`, position: [f.x, -f.h], normal: [0, 1], height: f.h, curve: { deflection: [0, f.defl], reaction: [0, f.RR * 1e3] }, ratedReaction: f.RR * 1e3 };
  });
  const ids = [...new Set(lines.map((l) => l.bollardId))];
  const windCurves = genericCoefficientCurves({ type: inp.vtype, kind: 'wind', Cx0: inp.CxW, Cy0: inp.CyW, Cxy0: inp.CxyW });
  const curCurves = genericCoefficientCurves({ type: inp.vtype, kind: 'current', Cx0: inp.CxC, Cy0: inp.CyC, Cxy0: inp.CxyC, shallowFactor: inp.shallow });
  const model = { vessel, pose0: { x: 0, y: inp.beam / 2, headingDeg: 0 }, lines, fenders, reference: { waterLevel: 0, draft: inp.draft }, condition: { waterLevel: inp.tide, draft: inp.draftCond > 0 ? inp.draftCond : inp.draft },
    limits: { steel: inp.limSteel, synthetic: inp.limSyn, polyamide: inp.limPa }, bollards: ids.map((id) => ({ id, SWL: inp.swl * 1e3 })), armEnvelope: { surge: inp.armSurge, sway: inp.armSway } };
  const env = (windFrom) => ({
    wind: { speed: inp.Uw, fromDeg: windFrom, curves: windCurves, rho: inp.rhoAir },
    current: { speed: inp.Uc, fromDeg: inp.curFrom, curves: curCurves, rho: inp.rhoW },
    extra: { surge: inp.extraSurge * 1e3, sway: inp.extraSway * 1e3, yaw: inp.extraYaw * 1e3 },
  });
  return { model, env };
}

export default {
  id: 's06-mooring-lines', suite: 6, title: 'Mooring line forces and bollard loads',
  summary: 'Line tensions, bollard and quick-release-hook reactions, fender loads and vessel excursion for a ship moored at dolphins under wind, current and additional vessel loads.',
  covers: [2], fidelity: 'reduced-order',
  standards: ['OCIMF Mooring Equipment Guidelines, 4th ed. (MEG4, 2018) §3–5', 'BS 6349-4:2014 §7–8', 'PIANC WG153 (2016) Recommendations for the design of marine oil and petrochemical terminals', 'ROM 2.0-11'],
  models: ['Static equilibrium in surge, sway and yaw', 'Elastic tension-only mooring lines with tails', 'Catenary line equation (sag of the governing line)', 'Compression-only fenders', 'OCIMF-form wind and current force coefficients', 'Vector resolution of line tension at bollards', 'Quasi-static mooring analysis'],
  limits: 'Static three-degree-of-freedom model: rigid ship at level trim, no wave-frequency or low-frequency dynamic amplification, no passing-ship regression (enter such loads as additional forces). The built-in wind and current coefficient curves are generic smooth shapes of realistic magnitude, NOT OCIMF or SIGTTO data — supply project coefficients for design. Line-load limits (OCIMF MEG) are overridable defaults. Catenary sag of the lines is neglected (taut lines); dynamic line equations and coupled vessel motion are outside this module. ' + NOT_HIFI,
  publishes: ['dolphin.lineTmax', 'dolphin.mooringLoad', 'dolphin.mooringUplift', 'dolphin.mooredFenderReaction'],
  searchTerms: 'mooring analysis berth OCIMF line tension bollard load dolphin',
  inputs: [
    { key: 'vtype', label: 'Vessel type (generic coefficient shape)', type: 'select', default: 'tanker', group: 'Vessel', role: 'load', options: [{ value: 'tanker', label: 'Tanker' }, { value: 'lng', label: 'Gas carrier' }, { value: 'container', label: 'Container ship' }, { value: 'bulk', label: 'Bulk carrier' }] },
    { key: 'LBP', label: 'Length between perpendiculars', unit: 'm', type: 'number', default: 220, min: 20, max: 420, group: 'Vessel', role: 'geometry' },
    { key: 'beam', label: 'Beam', unit: 'm', type: 'number', default: 32.2, min: 4, max: 70, group: 'Vessel', role: 'geometry' },
    { key: 'draft', label: 'Draft when the lines were tended', unit: 'm', type: 'number', default: 8, min: 1, max: 26, group: 'Vessel', role: 'ic', help: 'Reference condition at which the pretension is set. Ballast draft gives the largest windage.' },
    { key: 'draftCond', label: 'Draft analysed (0 = same)', unit: 'm', type: 'number', default: 0, min: 0, max: 26, group: 'Vessel', role: 'ic', help: 'A different draft (loading or discharging without re-tending) changes fairlead height and line tensions.' },
    { key: 'AT', label: 'Head-on windage area', unit: 'm²', type: 'number', default: 900, min: 10, max: 5000, group: 'Vessel', role: 'load' },
    { key: 'AL', label: 'Broadside windage area', unit: 'm²', type: 'number', default: 4200, min: 20, max: 20000, group: 'Vessel', role: 'load', uncertain: { dist: 'normal', cov: 0.05 } },
    { key: 'pbFrac', label: 'Parallel mid-body half-length / LBP', unit: '-', type: 'number', default: 0.35, min: 0.1, max: 0.5, group: 'Vessel', role: 'geometry', help: 'Fenders outside the parallel body are not in contact.' },
    { key: 'lines', label: 'Mooring lines', type: 'table', group: 'Mooring arrangement', role: 'bc',
      columns: [{ key: 'name', label: 'Line', type: 'text' }, { key: 'dol', label: 'Dolphin / hook', type: 'text' },
        { key: 'fx', label: 'Fairlead x (ship, fwd +)', unit: 'm', type: 'number' }, { key: 'fy', label: 'Fairlead y (ship, port +)', unit: 'm', type: 'number' }, { key: 'fz', label: 'Fairlead above keel', unit: 'm', type: 'number' },
        { key: 'bx', label: 'Hook x', unit: 'm', type: 'number' }, { key: 'by', label: 'Hook y', unit: 'm', type: 'number' }, { key: 'bz', label: 'Hook elevation', unit: 'm', type: 'number' },
        { key: 'n', label: 'Lines', unit: '-', type: 'number' }, { key: 'MBL', label: 'Breaking strength each', unit: 'kN', type: 'number' }, { key: 'material', label: 'Material (steel / synthetic / polyamide)', type: 'text' },
        { key: 'EA', label: 'Axial stiffness EA each', unit: 'kN', type: 'number' }, { key: 'pre', label: 'Pretension each', unit: 'kN', type: 'number' }, { key: 'tailL', label: 'Tail length', unit: 'm', type: 'number' }, { key: 'tailEA', label: 'Tail EA each', unit: 'kN', type: 'number' }],
      default: LINES_DEFAULT,
      help: 'Berth axes: x along the fender line (towards the bow), y from the fender line towards the vessel, elevations above the water level at tending. The vessel lies starboard side to, midships at x = 0.' },
    { key: 'fenders', label: 'Breasting fenders', type: 'table', group: 'Mooring arrangement', role: 'bc',
      columns: [{ key: 'name', label: 'Dolphin', type: 'text' }, { key: 'x', label: 'Position x', unit: 'm', type: 'number' }, { key: 'RR', label: 'Rated reaction', unit: 'kN', type: 'number' }, { key: 'defl', label: 'Rated deflection', unit: 'm', type: 'number' }, { key: 'h', label: 'Fender height', unit: 'm', type: 'number' }],
      default: FENDERS_DEFAULT, help: 'Represented by a secant-linear spring up to the rated point; the fender module gives the full curve.' },
    { key: 'wLine', label: 'Line weight in air', unit: 'kg/m', type: 'number', default: 7, min: 0.1, max: 60, group: 'Mooring arrangement', role: 'material', help: 'About 7 kg/m for 42 mm steel wire, 1–2 kg/m for HMPE or polyester of the same strength. Used for the catenary sag of the governing line only.' },
    { key: 'swl', label: 'Bollard / hook safe working load', unit: 'kN', type: 'number', default: 1500, min: 10, max: 20000, group: 'Mooring arrangement', role: 'criteria', design: true },
    { key: 'Uw', label: 'Wind speed (30 s gust at 10 m)', unit: 'm/s', type: 'number', default: 30.9, min: 0, max: 70, group: 'Environment', role: 'environment', bind: 'site.Uw100', uncertain: { dist: 'gumbel', cov: 0.12 }, help: 'OCIMF standard criterion: 60 knots (30.9 m/s) from any direction.' },
    { key: 'windFrom', label: 'Wind from (0 = from the bow side along +x, 270 = off the berth)', unit: 'deg', type: 'number', default: 270, min: 0, max: 360, group: 'Environment', role: 'environment', help: 'Direction the wind comes from, anticlockwise from the berth x axis. 270° blows the vessel off the dolphins and loads the breast lines; 90° presses it on the fenders.' },
    { key: 'Uc', label: 'Current speed', unit: 'm/s', type: 'number', default: 1.0, min: 0, max: 4, group: 'Environment', role: 'environment', bind: 'site.Uc100', uncertain: { dist: 'normal', cov: 0.15 } },
    { key: 'curFrom', label: 'Current from', unit: 'deg', type: 'number', default: 0, min: 0, max: 360, group: 'Environment', role: 'environment', help: 'Same convention as the wind. OCIMF: 3 kn at 0°/180°, 2 kn at 10°/170°, 0.75 kn broadside.' },
    { key: 'tide', label: 'Water level relative to tending', unit: 'm', type: 'number', default: 0, min: -8, max: 8, group: 'Environment', role: 'bc', help: 'Rise (+) or fall (−) since the lines were tended.' },
    { key: 'rhoAir', label: 'Air density', unit: 'kg/m³', type: 'number', default: 1.28, min: 1.1, max: 1.4, group: 'Environment', role: 'environment', help: '1.28 kg/m³ is the value paired with the OCIMF wind coefficients.' },
    { key: 'rhoW', label: 'Water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Environment', role: 'environment', bind: 'site.rhoW' },
    { key: 'CxW', label: 'Wind coefficient, longitudinal', unit: '-', type: 'number', default: 0.9, min: 0, max: 2, group: 'Force coefficients (user data)', role: 'load', help: 'Peak magnitude of the generic curve. Generic shape, not OCIMF data: enter the project coefficients for the vessel and loading condition.' },
    { key: 'CyW', label: 'Wind coefficient, lateral', unit: '-', type: 'number', default: 0.9, min: 0, max: 2, group: 'Force coefficients (user data)', role: 'load', uncertain: { dist: 'normal', cov: 0.1 } },
    { key: 'CxyW', label: 'Wind coefficient, yaw moment', unit: '-', type: 'number', default: 0.12, min: 0, max: 0.5, group: 'Force coefficients (user data)', role: 'load' },
    { key: 'CxC', label: 'Current coefficient, longitudinal', unit: '-', type: 'number', default: 0.04, min: 0, max: 0.5, group: 'Force coefficients (user data)', role: 'load' },
    { key: 'CyC', label: 'Current coefficient, lateral (deep water)', unit: '-', type: 'number', default: 0.6, min: 0, max: 5, group: 'Force coefficients (user data)', role: 'load' },
    { key: 'CxyC', label: 'Current coefficient, yaw moment', unit: '-', type: 'number', default: 0.15, min: 0, max: 1, group: 'Force coefficients (user data)', role: 'load' },
    { key: 'shallow', label: 'Shallow-water multiplier on lateral current load', unit: '-', type: 'number', default: 2.0, min: 1, max: 8, group: 'Force coefficients (user data)', role: 'load', help: 'Lateral current load grows several-fold as depth/draft falls towards 1.1; take from the coefficient source.' },
    { key: 'extraSurge', label: 'Additional surge force (ship axes)', unit: 'kN', type: 'number', default: 0, min: -1e5, max: 1e5, group: 'Additional vessel loads', role: 'bc', help: 'Peak passing-ship, wave-drift or seiche force from a dedicated calculation.' },
    { key: 'extraSway', label: 'Additional sway force (to port +)', unit: 'kN', type: 'number', default: 0, min: -1e5, max: 1e5, group: 'Additional vessel loads', role: 'bc', bind: 'vessel.Fmean' },
    { key: 'extraYaw', label: 'Additional yaw moment', unit: 'kN·m', type: 'number', default: 0, min: -1e7, max: 1e7, group: 'Additional vessel loads', role: 'bc' },
    { key: 'limSteel', label: 'Line limit, steel wire (fraction of MBL)', unit: '-', type: 'number', default: 0.55, min: 0.2, max: 1, group: 'Acceptance', role: 'criteria', help: `OCIMF MEG. ${CONFIRM}` },
    { key: 'limSyn', label: 'Line limit, synthetic (fraction of MBL)', unit: '-', type: 'number', default: 0.5, min: 0.2, max: 1, group: 'Acceptance', role: 'criteria', help: `OCIMF MEG. ${CONFIRM}` },
    { key: 'limPa', label: 'Line limit, polyamide (fraction of MBL)', unit: '-', type: 'number', default: 0.45, min: 0.2, max: 1, group: 'Acceptance', role: 'criteria', help: `OCIMF MEG. ${CONFIRM}` },
    { key: 'armSurge', label: 'Allowable surge at the manifold', unit: 'm', type: 'number', default: 3, min: 0.1, max: 20, group: 'Acceptance', role: 'criteria' },
    { key: 'armSway', label: 'Allowable sway at the manifold', unit: 'm', type: 'number', default: 3, min: 0.1, max: 20, group: 'Acceptance', role: 'criteria' },
    { key: 'sweep', label: 'Sweep the wind around the compass', type: 'bool', default: true, group: 'Acceptance', role: 'control', help: 'Runs the wind from every 30° to find the governing direction (adds twelve equilibrium solutions).' },
  ],
  compute(inp) {
    const { model, env } = buildModel(inp);
    const r = run(() => berthEquilibrium({ ...model, ...env(inp.windFrom) }), 'Mooring equilibrium');
    const count = Object.fromEntries(model.lines.map((l) => [l.name, l.count]));
    const lines = r.lines.map((l) => ({ ...l, per: l.tension / count[l.name] }));
    const gov = lines.reduce((a, l) => (l.uc > a.uc ? l : a), lines[0]);
    const bGov = r.bollards.reduce((a, b) => (b.load > a.load ? b : a), r.bollards[0]);
    const fGov = r.fenders.length ? r.fenders.reduce((a, f) => (f.reaction > a.reaction ? f : a), r.fenders[0]) : null;
    const envF = r.environment.reduce((a, e) => ({ Fx: a.Fx + e.Fx, Fy: a.Fy + e.Fy, Mz: a.Mz + e.Mz }), { Fx: 0, Fy: 0, Mz: 0 });
    let dirs = [], ucDir = [], govDir = inp.windFrom, govUc = gov.uc, fendDir = [];
    if (inp.sweep) {
      const cases = []; for (let d = 0; d < 360; d += 30) cases.push({ name: String(d), ...env(d) });
      const sw = runBerthCases(model, cases);
      for (const c of sw.cases) { dirs.push(Number(c.name)); ucDir.push(c.converged ? round(c.maxLineUc) : 0); fendDir.push(c.converged ? round(c.maxFenderUc) : 0); if (c.converged && c.maxLineUc > govUc) { govUc = c.maxLineUc; govDir = Number(c.name); } }
    }
    const maxUplift = Math.max(...r.bollards.map((b) => b.vertical), 0);
    // Catenary sag of the governing line between level supports: z(X/2) = (H/w)·[cosh(wX/2H) − 1] (taut-line model check)
    const Hline = Math.max(gov.horizontalTension / count[gov.name], 1), span = gov.length * Math.cos(gov.verticalAngle);
    const sag = inextensibleCatenaryShape({ H: Hline, w: inp.wLine * 9.80665, x: span / 2 }).z;
    const checks = [
      check(`Line load, ${gov.name}`, gov.per / 1e3, gov.limitFraction * model.lines.find((l) => l.name === gov.name).MBL / count[gov.name] / 1e3, 'kN', 'OCIMF MEG4 (2018) §3 line load limits (55 % steel, 50 % synthetic, 45 % polyamide of MBL)', 'Keeps the most loaded line below its working limit. A parted line under load recoils with lethal energy and transfers its load to the neighbours.'),
      check(`Bollard / hook load, ${bGov.id}`, bGov.load / 1e3, inp.swl, 'kN', 'BS 6349-4:2014 §8; OCIMF MEG4 §4', 'Vector sum of the lines on one hook assembly against its safe working load. Overload fails the hook, its foundation bolts or the dolphin deck.'),
      check('Manifold surge', Math.abs(r.excursion.manifoldX), inp.armSurge, 'm', 'OCIMF MEG4 §3; loading-arm envelope', 'Movement along the berth at the manifold against the working envelope of loading arms or hoses. Beyond it the arms disconnect or are damaged.'),
      check('Manifold sway', Math.abs(r.excursion.manifoldY), inp.armSway, 'm', 'OCIMF MEG4 §3; loading-arm envelope', 'Movement off the berth at the manifold against the loading-arm envelope.'),
    ];
    if (fGov) checks.push(check(`Fender reaction, ${fGov.name}`, fGov.reaction / 1e3, fGov.ratedReaction / 1e3, 'kN', 'BS 6349-4:2014 §6', 'Reaction of the breasting fender under on-berth wind and current against its rated reaction. Beyond it the fender is over-compressed.'));
    if (inp.sweep) checks.push({ name: `Line load, worst wind direction (${govDir}°)`, demand: round(govUc, 4), capacity: 1, unit: 'fraction of limit', uc: Math.min(govUc, 99), ref: 'OCIMF MEG4 §3 (wind from any direction)', explain: 'The largest line utilisation found with the wind swept around the compass. The mooring pattern must hold for every direction, not only the one entered.' });
    const w = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = w.uc <= 1;
    const mblNeed = gov.per / gov.limitFraction / 1e3;
    const advice = [
      adv('engineering', ok ? 'ok' : 'fail', `Most loaded line ${gov.name}: ${fmt(gov.per / 1e3)} kN per line (${fmt(gov.utilisation * 100)} % of MBL, limit ${fmt(gov.limitFraction * 100)} %); largest hook load ${fmt(bGov.load / 1e3)} kN on ${bGov.id}.`,
        ok ? 'Re-run with the project wind and current coefficients for ballast and laden draft and with the tide at both extremes; check the governing direction from the polar plot.' : (w.name.startsWith('Line') ? `Raise the breaking strength to at least ${fmt(mblNeed)} kN per line, add a line to the ${gov.name} group, or improve the lead (flatter, more perpendicular breast lines).` : w.name.startsWith('Bollard') ? `Provide hooks of at least ${fmt(bGov.load / 1e3)} kN or split the lines over two dolphins.` : 'Stiffen the spring or breast lines (shorter tails, higher pretension) to bring the excursion inside the loading-arm envelope.'),
        'OCIMF MEG4 (2018) §3–5; BS 6349-4:2014 §7–8',
        'Lines sized to the working limit rather than oversized are lighter to handle and replaced on a normal 5–8 year cycle: every unnecessary line is roughly a quarter tonne of wire or fibre per replacement.'),
      adv('economic', ok ? 'info' : 'warn', `Mooring dolphin design load from this case: ${fmt(bGov.horizontal / 1e3)} kN horizontal and ${fmt(maxUplift / 1e3)} kN uplift per dolphin; vessel excursion ${fmt(r.excursion.surge)} m surge, ${fmt(r.excursion.sway)} m sway.`,
        'Wind limits for staying alongside are the lever on downtime: derive the wind speed at which the governing line reaches its limit and compare the lost berth hours with the cost of stronger hooks and dolphins.', `${ECON_BASIS}; berth-availability (downtime) assessment per PIANC WG153`,
        `Each avoided emergency departure saves tug and vessel fuel: ${vesselCo2Text(0.25)}`),
    ];
    const pub = { 'dolphin.lineTmax': gov.per, 'dolphin.mooringLoad': bGov.horizontal, 'dolphin.mooringUplift': maxUplift };
    if (fGov) pub['dolphin.mooredFenderReaction'] = fGov.reaction;
    return {
      kpis: [
        { key: 'Tmax', label: 'Largest line tension (per line)', unit: 'kN', value: gov.per / 1e3, explain: 'Tension in the most utilised single line. Compare with the working limit, not the breaking strength.' },
        { key: 'ucLine', label: 'Largest line utilisation', unit: '-', value: gov.uc, limit: 1, sense: 'max', explain: 'Line tension over its allowable load (fraction of MBL by material). Above 1 the line exceeds the OCIMF working limit.' },
        { key: 'bollard', label: 'Largest bollard / hook load', unit: 'kN', value: bGov.load / 1e3, limit: inp.swl, sense: 'max', explain: 'Resultant of all lines on the most loaded hook assembly: the mooring design load of that dolphin.' },
        { key: 'bollardH', label: 'Horizontal component on that dolphin', unit: 'kN', value: bGov.horizontal / 1e3, explain: 'Horizontal part of the hook load, carried by the piles in bending and push–pull.' },
        { key: 'uplift', label: 'Largest uplift on a dolphin', unit: 'kN', value: maxUplift / 1e3, explain: 'Upward vertical component of the line pull (high fairleads at low water or light draft). It reduces pile compression and can govern tension piles.' },
        { key: 'sag', label: 'Catenary sag of the governing line', unit: 'm', value: sag, explain: 'Mid-span sag of the most loaded line under its own weight at the computed tension. A sag that is small against the line length confirms the taut elastic-line model; a large sag (slack heavy wire) means the real line is softer than modelled.' },
        { key: 'surge', label: 'Vessel surge', unit: 'm', value: r.excursion.surge, explain: 'Movement along the berth from the tended position; positive towards the bow.' },
        { key: 'sway', label: 'Vessel sway', unit: 'm', value: r.excursion.sway, explain: 'Movement normal to the berth; positive away from the dolphins.' },
        { key: 'yaw', label: 'Vessel yaw', unit: 'deg', value: r.excursion.yawDeg, explain: 'Rotation in plan from the tended heading. Large yaw indicates unbalanced bow and stern restraint.' },
        { key: 'fender', label: 'Largest fender reaction', unit: 'kN', value: fGov ? fGov.reaction / 1e3 : 0, explain: 'Reaction on the most loaded breasting dolphin with the vessel pressed on by wind and current: the sustained (not impact) lateral load.' },
        { key: 'Fx', label: 'Environmental surge force', unit: 'kN', value: envF.Fx / 1e3, explain: 'Sum of wind, current and additional force along the ship, forward positive.' },
        { key: 'Fy', label: 'Environmental sway force', unit: 'kN', value: envF.Fy / 1e3, explain: 'Sum of wind, current and additional force across the ship, positive to port (off the berth for a starboard-side-to vessel).' },
        { key: 'Mz', label: 'Environmental yaw moment', unit: 'kN·m', value: envF.Mz / 1e3, explain: 'Moment about midships, bow to port positive.' },
        { key: 'ucSweep', label: 'Line utilisation, worst wind direction', unit: '-', value: govUc, limit: 1, sense: 'max', explain: 'Largest line utilisation over all wind directions (equal to the entered direction when the sweep is off).' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Line tension as a fraction of breaking strength', categories: lines.map((l) => l.name), series: [{ name: 'Tension / MBL', y: lines.map((l) => round(l.utilisation)) }, { name: 'Limit', y: lines.map((l) => l.limitFraction) }], y: { label: 'Tension / MBL', unit: '-' },
          explain: 'Load sharing between the lines. Bars of very different height show an unbalanced pattern: short stiff lines attract load, long soft ones do little.' },
        ...(inp.sweep ? [{ type: 'polar', title: 'Largest line utilisation against wind direction', series: [{ name: 'Line load / limit', thetaDeg: [...dirs, dirs[0]], r: [...ucDir, ucDir[0]] }, { name: 'Limit', thetaDeg: [...dirs, dirs[0]], r: dirs.map(() => 1).concat(1) }],
          explain: 'Utilisation of the most loaded line for the design wind from each direction. Lobes outside the unit circle are directions for which the pattern is inadequate.' }] : []),
        { type: 'scatter', title: 'Mooring plan', x: { label: 'Along berth', unit: 'm' }, y: { label: 'From fender line', unit: 'm' },
          series: [{ name: 'Fairleads', x: model.lines.map((l) => round(l.fairlead[0])), y: model.lines.map((l) => round(inp.beam / 2 + l.fairlead[1])) }, { name: 'Hooks', x: model.lines.map((l) => l.bollard[0]), y: model.lines.map((l) => l.bollard[1]) }, { name: 'Fenders', x: model.fenders.map((f) => f.position[0]), y: model.fenders.map(() => 0) }],
          explain: 'Plan of fairleads, hooks and fenders in berth axes at the tended position, to check the entered geometry: breast lines should be near-perpendicular, springs near-parallel to the hull.' },
      ],
      tables: [
        { title: 'Line results', columns: [{ key: 'name', label: 'Line' }, { key: 'T', label: 'Tension per line', unit: 'kN' }, { key: 'pct', label: 'of MBL', unit: '%' }, { key: 'uc', label: 'Utilisation' }, { key: 'va', label: 'Vertical angle', unit: 'deg' }, { key: 'ha', label: 'Plan angle', unit: 'deg' }, { key: 'len', label: 'Outboard length', unit: 'm' }],
          rows: lines.map((l) => ({ name: l.name, T: round(l.per / 1e3), pct: round(l.utilisation * 100), uc: round(l.uc), va: round(l.verticalAngle / DEG), ha: round(l.horizontalAngle / DEG), len: round(l.length) })),
          explain: 'One row per line group: tension in each line, its share of breaking strength and its lead angles (vertical angle negative when the hook is below the fairlead).' },
        { title: 'Dolphin (hook) reactions', columns: [{ key: 'id', label: 'Dolphin' }, { key: 'R', label: 'Resultant', unit: 'kN' }, { key: 'Fx', label: 'Along berth', unit: 'kN' }, { key: 'Fy', label: 'Normal to berth', unit: 'kN' }, { key: 'Fz', label: 'Vertical (up +)', unit: 'kN' }, { key: 'uc', label: 'Load / SWL' }],
          rows: r.bollards.map((b) => ({ id: b.id, R: round(b.load / 1e3), Fx: round(b.force[0] / 1e3), Fy: round(b.force[1] / 1e3), Fz: round(b.force[2] / 1e3), uc: round(b.uc ?? 0) })),
          explain: 'Force components each dolphin receives from its lines, in berth axes: the load case to apply at hook level in the structural modules.' },
      ],
      notes: [`${w.name} governs (utilisation ${fmt(w.uc)}).`, 'Wind and current coefficients are generic curve shapes scaled by the entered magnitudes — not OCIMF data.', ...(inp.sweep ? [`Worst wind direction for the lines: from ${govDir}°.`] : [])],
      advice, publish: pub, quantities: {},
    };
  },
  verify: [
    // Beam wind off the berth: Fy = ½·1.28·1.0·20²·3000 = 768 kN; two perpendicular horizontal breast lines share it: T = 384 kN each.
    { name: 'Two horizontal breast lines under beam wind: T = F/2', source: 'Statics; wind force ½·ρ·Cy·V²·A_L with ρ = 1.28 kg/m³ (OCIMF force definition)',
      inputs: { beam: 30, draft: 10, AL: 3000, AT: 500, Uw: 20, windFrom: 270, Uc: 0, CyW: 1, CxW: 0, CxyW: 0, sweep: false, fenders: [], swl: 5000, armSurge: 20, armSway: 20,
        lines: [{ name: 'B1', dol: 'A', fx: 50, fy: -15, fz: 15, bx: 50, by: -30, bz: 5, n: 1, MBL: 5000, material: 'steel', EA: 1e6, pre: 0, tailL: 0, tailEA: 0 },
          { name: 'B2', dol: 'B', fx: -50, fy: -15, fz: 15, bx: -50, by: -30, bz: 5, n: 1, MBL: 5000, material: 'steel', EA: 1e6, pre: 0, tailL: 0, tailEA: 0 },
          { name: 'S1', dol: 'C', fx: 60, fy: -15, fz: 15, bx: 0, by: -1, bz: 5, n: 1, MBL: 5000, material: 'steel', EA: 1e6, pre: 0, tailL: 0, tailEA: 0 },
          { name: 'S2', dol: 'D', fx: -60, fy: -15, fz: 15, bx: 0, by: -1, bz: 5, n: 1, MBL: 5000, material: 'steel', EA: 1e6, pre: 0, tailL: 0, tailEA: 0 }] },
      expect: [{ kpi: 'Fy', value: 768, tol: 1e-6 }, { kpi: 'Tmax', value: 384, tol: 0.01 }] },
    // Same with the hooks 10 m below the fairleads over 30 m horizontal: cos φ = 30/√(30² + 10²) = 0.94868; T = 384/0.94868 = 404.8 kN;
    // the horizontal component on each dolphin stays 384 kN.
    { name: 'Inclined breast lines: T = F/(2·cos φ)', source: 'Resolution of line tension for a known fairlead and hook geometry (statics)',
      inputs: { beam: 30, draft: 10, AL: 3000, AT: 500, Uw: 20, windFrom: 270, Uc: 0, CyW: 1, CxW: 0, CxyW: 0, sweep: false, fenders: [], swl: 5000, armSurge: 20, armSway: 20,
        lines: [{ name: 'B1', dol: 'A', fx: 50, fy: -15, fz: 15, bx: 50, by: -30, bz: -5, n: 1, MBL: 5000, material: 'steel', EA: 1e6, pre: 0, tailL: 0, tailEA: 0 },
          { name: 'B2', dol: 'B', fx: -50, fy: -15, fz: 15, bx: -50, by: -30, bz: -5, n: 1, MBL: 5000, material: 'steel', EA: 1e6, pre: 0, tailL: 0, tailEA: 0 },
          { name: 'S1', dol: 'C', fx: 60, fy: -15, fz: 15, bx: 0, by: -1, bz: 5, n: 1, MBL: 5000, material: 'steel', EA: 1e6, pre: 0, tailL: 0, tailEA: 0 },
          { name: 'S2', dol: 'D', fx: -60, fy: -15, fz: 15, bx: 0, by: -1, bz: 5, n: 1, MBL: 5000, material: 'steel', EA: 1e6, pre: 0, tailL: 0, tailEA: 0 }] },
      expect: [{ kpi: 'Tmax', value: 404.77, tol: 0.01 }, { kpi: 'bollardH', value: 384, tol: 0.01 }] },
  ],
};
