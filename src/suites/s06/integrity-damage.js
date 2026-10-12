// Structural integrity and damage assessment of an existing dolphin: residual strength with measured wall loss,
// dents and defective connections (plastic pushover of the pile group), predicted damage from an accidental
// impact, and the operating envelope the damaged structure still supports.
import { pushoverAnalysis } from '../../engine/structural/pushover.js';
import { need, run, fmt, round, cap, adv, check, linspace, pileInputs, soilInputs, layoutInputs, capInputs, frameInputs, fixityLength, pileProps, profileFrom, pilesFrom, frameModel, steelCo2Text, vesselCo2Text, CONFIRM, NOT_HIFI, ECON_BASIS } from './_common.js';

const COND_DEFAULT = [{ tMeas: 23, dent: 0, knock: 1 }, { tMeas: 25, dent: 0, knock: 1 }, { tMeas: 20, dent: 70, knock: 1 }, { tMeas: 24, dent: 0, knock: 0.8 }];

/**
 * Plastic capacities of a dented tube, from first principles (thin ring, radius r, wall t).
 * A dent of depth δ flattens the arc of half-angle α = acos(1 − 2δ/D). Taking that arc as ineffective and the
 * dent on the compression side, the remaining arc spans θ = α … 2π − α (θ from the dent centre). The plastic
 * neutral axis halves the remaining area: β − α = π − β → β = (π + α)/2. With lever arm r·cos θ,
 *   Mp,d = 2·fy·t·r²·[∫α..β cos θ dθ − ∫β..π cos θ dθ] = 2·fy·t·r²·(2 sin β − sin α) = 4·fy·t·r²·[cos(α/2) − ½ sin α],
 * so Mp,d/Mp = cos(α/2) − ½ sin α and Np,d/Np = 1 − α/π. (Lower bound: the flattened wall still carries some stress.)
 */
export function dentFactors(dent, D) {
  if (!(dent > 0)) return { m: 1, n: 1, alpha: 0 };
  const alpha = Math.acos(Math.max(1 - 2 * dent / D, -1));
  return { m: Math.max(Math.cos(alpha / 2) - 0.5 * Math.sin(alpha), 0.02), n: Math.max(1 - alpha / Math.PI, 0.02), alpha };
}

function pushover(inp, P, piles, Le, fixed, states) {
  // states: per pile { t (m), fm, fn } — wall thickness and capacity factors on Mp and Np
  const fm = frameModel({ piles, P, Le, nSub: inp.nSub, fixedHead: fixed });
  const model = fm.model;
  piles.forEach((p, k) => {
    const s = states[k], Di = P.D - 2 * s.t;
    const A = Math.PI / 4 * (P.D ** 2 - Di ** 2), Zp = (P.D ** 3 - Di ** 3) / 6;
    model.sections[`pile${k}`] = { type: 'tube', D: P.D, t: s.t };
    for (let e = 0; e < inp.nSub; e++) {
      const m = model.members[k * inp.nSub + e];
      m.section = `pile${k}`; m.K = fixed ? 1 : 2; m.bucklingLength = Le;
      if (s.fm < 1 || s.fn < 1) m.capacity = { Np: A * P.fy * s.fn, Mpy: Zp * P.fy * s.fm, Mpz: Zp * P.fy * s.fm, Pt: A * P.fy * s.fn };
    }
  });
  // the cap is not allowed to hinge: give its tie members a capacity far above the piles
  const big = 1e3 * P.A * P.fy, bigM = 1e3 * P.Mp;
  for (const m of model.members) if (m.group === 'cap') m.capacity = { Np: big, Mpy: bigM, Mpz: bigM, Pc: big, Pt: big };
  const W = (inp.capB * inp.capL * inp.capT * inp.gammaConc + inp.Wequip) * 1e3, F = inp.Fd * 1e3;
  model.loadCases = [{ name: 'grav', nodal: [{ node: fm.centre, F: [0, 0, -W, 0, 0, 0] }] }, { name: 'lat', nodal: [{ node: fm.centre, F: [0, F, 0, -F * inp.zLoad, 0, 0] }] }];
  const r = run(() => pushoverAnalysis(model, { pattern: 'lat', base: 'grav', control: { node: fm.centre, dof: 1 }, pDelta: !!inp.pDelta, shearDeformation: !!inp.shear, maxEvents: 60 }), 'Pushover analysis');
  // energy absorbed up to each point of the curve: ∫ H dδ
  const H = r.curve.loadFactor.map((l) => l * F), d = r.curve.displacement.map(Math.abs), E = [0];
  for (let i = 1; i < H.length; i++) E.push(E[i - 1] + 0.5 * (H[i] + H[i - 1]) * (d[i] - d[i - 1]));
  if (d[0] !== 0 || H[0] !== 0) { H.unshift(0); d.unshift(0); E.unshift(0); for (let i = 1; i < E.length; i++) E[i] = E[i - 1] + 0.5 * (H[i] + H[i - 1]) * (d[i] - d[i - 1]); }
  const iFirst = H.findIndex((h) => h >= r.firstEventLoadFactor * F * (1 - 1e-9));
  return { r, H, d, E, Hu: r.ultimateLoadFactor * F, H1: r.firstEventLoadFactor * F, iFirst: iFirst < 0 ? H.length - 1 : iFirst, model };
}

export default {
  id: 's06-integrity-damage', suite: 6, title: 'Structural integrity and damage assessment',
  summary: 'Residual lateral strength of a corroded or damaged dolphin, predicted damage and permanent set from an accidental impact, repair requirement and the vessel size, berthing velocity and mooring load it can still accept.',
  covers: [12], fidelity: 'numerical',
  standards: ['ISO 19902:2020 §12 (assessment of existing structures), §13.7 (dented members)', 'API RP 2SIM (2014) structural integrity management', 'BS 6349-2:2019 §9 (assessment and repair)', 'PIANC WG103 Life Cycle Management of Port Structures (2008)', 'DNV-RP-C204 (2019) accidental loads (energy absorption)'],
  models: ['Nonlinear frame analysis: event-to-event plastic pushover with P–Δ', 'Plastic hinge interaction of tubular sections', 'Dented-tube plastic capacity (first-principles thin-ring derivation)', 'Residual strength ratio and damaged strength ratio', 'Energy balance of impact: fender, elastic and plastic work', 'Kinetic-energy berthing limit for the damaged structure'],
  limits: 'Piles are columns fixed at an equivalent depth; plastic hinges form at member ends with elastic–perfectly-plastic behaviour and a rotation limit set by the user. The dent model treats the flattened arc as ineffective (lower bound) and ignores dent growth; weld and connection defects enter as a strength factor from inspection judgement or fracture assessment. Local denting energy, rate effects and hull crushing are not modelled. ' + NOT_HIFI,
  publishes: ['dolphin.lateralCapacity', 'dolphin.residualStrengthRatio', 'dolphin.energyCapacityDamaged'],
  searchTerms: 'residual strength corroded pile dolphin damage assessment pushover dent',
  inputs: [
    ...layoutInputs(),
    ...pileInputs().filter((i) => i.key !== 'corr'),
    { key: 'condition', label: 'Inspection record per pile (same order as the layout)', type: 'table', group: 'Structural condition', role: 'ic',
      columns: [{ key: 'tMeas', label: 'Measured wall at the critical section', unit: 'mm', type: 'number' }, { key: 'dent', label: 'Dent depth', unit: 'mm', type: 'number' }, { key: 'knock', label: 'Strength factor for weld / connection defects (0 = ineffective)', unit: '-', type: 'number' }],
      default: COND_DEFAULT, help: 'Ultrasonic wall readings near the section of maximum moment, measured dent depth, and a factor from weld inspection (1 sound, 0 severed or cracked through). Missing rows take the as-built wall.' },
    { key: 'corrRate', label: 'Corrosion rate from here on', unit: 'mm/yr', type: 'number', default: 0.15, min: 0, max: 2, group: 'Structural condition', role: 'ic', uncertain: { dist: 'lognormal', cov: 0.4 }, help: 'About 0.1 mm/yr immersed with working cathodic protection lost, 0.2–0.4 mm/yr in the splash zone unprotected, more with accelerated low-water corrosion.' },
    { key: 'fenderFactor', label: 'Fender energy remaining (deterioration factor)', unit: '-', type: 'number', default: 0.8, min: 0, max: 1, group: 'Structural condition', role: 'ic', help: 'Observed fender deterioration: 1 as new, 0.7–0.9 aged rubber, 0 missing.' },
    { key: 'Efender', label: 'Fender energy capacity as new', unit: 'kJ', type: 'number', default: 755, min: 0, max: 1e5, group: 'Structural condition', role: 'ic', bind: 'dolphin.energyCapacity' },
    { key: 'y0', label: 'Observed residual head deflection', unit: 'mm', type: 'number', default: 0, min: 0, max: 3000, group: 'Structural condition', role: 'ic', help: 'Permanent lean recorded by survey; it uses up part of the allowable plastic displacement.' },
    ...capInputs(),
    ...soilInputs(),
    ...frameInputs().map((i) => (i.key === 'nSub' ? { ...i, default: 2, help: 'Hinges form at member ends; two elements per pile capture the head and base hinges, more refine the P–Δ shape.' } : i)),
    { key: 'pDelta', label: 'Include P–Δ in the pushover', type: 'bool', default: true, group: 'Frame model', role: 'control' },
    { key: 'Fd', label: 'Design lateral load', unit: 'kN', type: 'number', default: 1575, min: 1, max: 1e5, group: 'Loading', role: 'bc', bind: 'dolphin.fenderReaction', help: 'Reference load pattern scaled to collapse: the berthing reaction or line pull normal to the berth.' },
    { key: 'zLoad', label: 'Load height above the pile heads', unit: 'm', type: 'number', default: 1.5, min: -5, max: 20, group: 'Loading', role: 'bc' },
    { key: 'Eacc', label: 'Accidental impact energy', unit: 'kJ', type: 'number', default: 1200, min: 0, max: 1e6, group: 'Loading', role: 'bc', help: 'Collision or grossly abnormal berthing to be survived without collapse (e.g. twice the abnormal berthing energy, or a drifting-vessel scenario).' },
    { key: 'disp', label: 'Design vessel displacement', unit: 't', type: 'number', default: 75000, min: 100, max: 700000, group: 'Operating envelope', role: 'operation', bind: 'vessel.disp' },
    { key: 'CeCm', label: 'Product Ce·Cm·Cs·Cc of the berthing calculation', unit: '-', type: 'number', default: 1.03, min: 0.2, max: 2.5, group: 'Operating envelope', role: 'operation', help: 'From the berthing-energy module (eccentricity × added mass × softness × configuration).' },
    { key: 'Fab', label: 'Abnormal berthing factor', unit: '-', type: 'number', default: 1.25, min: 1, max: 3, group: 'Operating envelope', role: 'criteria', help: CONFIRM },
    { key: 'rsrReq', label: 'Required residual strength ratio', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Acceptance', role: 'criteria', help: `Collapse load over the design lateral load for continued operation. ISO 19902 / API RP 2SIM use 1.6–2.0 for manned or high-consequence platforms; 1.3–1.5 is common for berth structures with operational controls. ${CONFIRM}` },
    { key: 'thetaPl', label: 'Allowable plastic drift (permanent set / length to fixity)', unit: '-', type: 'number', default: 0.03, min: 0.002, max: 0.1, group: 'Acceptance', role: 'criteria', help: 'Rotation capacity of the hinges before local buckling or fracture: about 0.02–0.04 for compact tubes (D/t below 40–50), much less for slender or dented walls.' },
  ],
  mesh: { key: 'nSub', label: 'Elements per pile', levels: [1, 2, 4], qoi: 'Hu' },
  compute(inp, ctx) {
    const P = pileProps({ ...inp, corr: 0 }, ctx);
    const all = pilesFrom(inp);
    const fixedAll = inp.headFix === 'fixed' && all.length > 1;
    const profile = inp.LeUser > 0 ? null : profileFrom(inp, P.D);
    const Le = fixityLength(inp, P, profile, fixedAll);
    // condition per pile (mm → m)
    const cond = all.map((p, k) => {
      const c = (inp.condition || [])[k] || {};
      const t = Math.min(c.tMeas > 0 ? c.tMeas / 1e3 : P.t, P.t), knock = c.knock == null || c.knock === '' ? 1 : Math.min(Math.max(+c.knock, 0), 1);
      const df = dentFactors((c.dent || 0) / 1e3, P.D);
      need(c.dent == null || c.dent < 0.6 * P.D * 1e3, `Pile ${k + 1}: a dent of ${fmt(c.dent)} mm is more than 60 % of the diameter — treat the pile as ineffective (strength factor 0).`);
      return { id: k + 1, t, dent: c.dent || 0, knock, fm: knock * df.m, fn: knock * df.n, loss: 1 - t / P.t };
    });
    const live = cond.filter((c) => c.knock > 0 && c.t >= 0.004);
    need(live.length >= 1, 'No effective pile remains: every pile is flagged ineffective or has less than 4 mm of wall.');
    const pilesD = all.filter((p, k) => cond[k].knock > 0 && cond[k].t >= 0.004);
    const fixedD = inp.headFix === 'fixed' && pilesD.length > 1;
    const intact = pushover(inp, P, all, Le, fixedAll, all.map(() => ({ t: P.t, fm: 1, fn: 1 })));
    const dam = pushover(inp, P, pilesD, Le, fixedD, live);
    need(intact.Hu > 0 && dam.Hu > 0, 'The pile group cannot carry its own cap weight over this length to fixity (it is unstable under gravity with P–Δ before any lateral load is applied). Check the pile diameter, the length to fixity and the cap weight.');
    const F = inp.Fd * 1e3;
    const rsrI = intact.Hu / F, rsrD = dam.Hu / F, dsr = dam.Hu / intact.Hu;
    // Energy balance of the accidental impact (J)
    const Ef = inp.Efender * 1e3 * inp.fenderFactor;
    const Eel = dam.E[dam.iFirst];                         // to first yield
    const Eult = dam.E[dam.E.length - 1], dUlt = dam.d[dam.d.length - 1];
    const dAllow = Math.max(inp.thetaPl * Le - inp.y0 / 1e3, 0);  // further plastic displacement available
    const Ecap = Ef + Eult + dam.Hu * dAllow;               // collapse-prevention capacity
    const Es = Math.max(inp.Eacc * 1e3 - Ef, 0);
    let set = 0, level = 0;
    if (Es > Eel) { level = 1; if (Es > Eult) { set = (Es - Eult) / dam.Hu; level = 2; } else { // between first yield and mechanism: unload elastically from the curve point
      let i = 1; while (i < dam.E.length - 1 && dam.E[i] < Es) i++;
      const f = (Es - dam.E[i - 1]) / Math.max(dam.E[i] - dam.E[i - 1], 1e-9), dNow = dam.d[i - 1] + f * (dam.d[i] - dam.d[i - 1]), hNow = dam.H[i - 1] + f * (dam.H[i] - dam.H[i - 1]);
      const k0 = dam.H[1] / Math.max(dam.d[1], 1e-12); set = Math.max(dNow - hNow / k0, 0); } }
    if (level === 2 && set > dAllow) level = 3;
    const LEVELS = ['No structural damage: energy taken by the fender and elastic bending', 'Local yielding at the most loaded section: inspect, no repair expected', 'Plastic hinges with permanent set: straightening or pile replacement required', 'Beyond the hinge rotation capacity: collapse or loss of the dolphin'];
    // Operating envelope of the damaged dolphin
    const Hallow = dam.Hu / inp.rsrReq;
    const Eop = Ef + Eel;                                   // energy without permanent damage
    const vAllow = Math.sqrt(2 * Eop / (inp.Fab * inp.disp * 1e3 * inp.CeCm));
    const ms = linspace(0.25, 2, 15).map((f) => f * inp.disp);
    // Years until the residual strength ratio falls to the requirement (bending-governed: capacity ∝ wall)
    const tAvg = live.reduce((s, c) => s + c.t, 0) / live.length;
    const years = inp.corrRate > 0 && rsrD > inp.rsrReq ? (tAvg * 1e3) * (1 - inp.rsrReq / rsrD) / inp.corrRate : (rsrD > inp.rsrReq ? 200 : 0);
    const checks = [
      { name: 'Residual strength ratio of the damaged dolphin', demand: round(inp.rsrReq, 4), capacity: round(cap(rsrD), 4), unit: 'ratio (required / achieved)', uc: Math.min(inp.rsrReq / rsrD, 99), ref: 'ISO 19902:2020 §12; API RP 2SIM (2014) §6 (overridable acceptance value)', explain: 'Lateral collapse load of the structure as inspected over the design load. Below the required ratio the dolphin should be load-restricted, strengthened or repaired.' },
      check('Accidental impact energy against collapse-prevention capacity', inp.Eacc, Ecap / 1e3, 'kJ', 'DNV-RP-C204 (2019) §3 (energy absorption); ISO 19902:2020 §10', 'Energy of the accidental event against deteriorated fender, elastic work and plastic hinge work up to the rotation limit. If exceeded the dolphin collapses and the vessel reaches the structure behind.'),
      check('Largest wall loss against 50 % of as-built wall', Math.max(...cond.map((c) => c.loss)) * 100, 50, '%', 'API RP 2SIM (2014) §5 (anomaly criteria, project-specific); BS 6349-2:2019 §9', 'Local thinning beyond about half the wall invalidates the member model (holing, local buckling). Such sections need sleeving or replacement regardless of global strength.'),
    ];
    const w = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = w.uc <= 1;
    const worstPile = cond.reduce((a, c) => (c.fm * c.t < a.fm * a.t ? c : a), cond[0]);
    const sleeveT = Math.PI * P.D * 0.02 * 6 * 7.85; // 20 mm grouted sleeve, 6 m long
    const repair = level >= 2 ? 'Replace or straighten the hinged piles and re-survey the cap' : (rsrD < inp.rsrReq ? `Strengthen: sleeve pile ${worstPile.id} (and any pile below ${fmt(P.t * 0.7 * 1e3)} mm) or restrict the berth` : (dsr < 0.9 ? `Schedule sleeve or coating repair of pile ${worstPile.id}; restore cathodic protection` : 'No repair required: continue inspection'));
    const advice = [
      adv('engineering', ok ? 'ok' : 'fail', `Collapse load ${fmt(dam.Hu / 1e3)} kN as inspected against ${fmt(intact.Hu / 1e3)} kN as built (damaged strength ratio ${fmt(dsr)}); residual strength ratio ${fmt(rsrD)} on the design load of ${fmt(inp.Fd)} kN. Accidental impact of ${fmt(inp.Eacc)} kJ: ${LEVELS[level].toLowerCase()}${set > 0 ? `, permanent set about ${fmt(set * 1e3)} mm` : ''}.`,
        `${repair}. Until then limit the lateral load to ${fmt(Hallow / 1e3)} kN and berthing of the ${fmt(inp.disp)} t vessel to ${fmt(vAllow)} m/s.`,
        'ISO 19902:2020 §12 and §13.7; API RP 2SIM (2014); BS 6349-2:2019 §9; DNV-RP-C204 (2019)',
        `A grouted sleeve repair of one pile is of the order of ${steelCo2Text(sleeveT)} Replacing the pile instead would use about ${fmt(Math.PI * P.D * P.t * (inp.Lfree + inp.Lemb) * 7.85)} t of steel plus piling-spread days.`),
      adv('economic', ok ? 'info' : 'warn', `At ${fmt(inp.corrRate)} mm/yr the residual strength ratio reaches the required ${fmt(inp.rsrReq)} in about ${fmt(Math.min(years, 200))} years; the dolphin currently retains ${fmt(dsr * 100)} % of its as-built strength.`,
        years < 10 ? 'Plan the repair within the next maintenance window and restore corrosion protection now: the cost of anodes or coating is small against a berth outage for pile replacement.' : 'Set the next detailed inspection at about one third of the remaining period and compare retrofit cathodic protection against accepting the wall loss to end of life.',
        `Risk-based inspection per API RP 2SIM / DNV-RP-G101; ${ECON_BASIS}`, `One diver or ROV inspection campaign is about ${vesselCo2Text(1.5)}`),
    ];
    const dMax = Math.max(dam.d[dam.d.length - 1], intact.d[intact.d.length - 1]) + dAllow;
    return {
      kpis: [
        { key: 'Hu', label: 'Collapse load, as inspected', unit: 'kN', value: dam.Hu / 1e3, explain: 'Lateral load at which the damaged pile group forms a plastic mechanism (pushover with the measured walls, dents and defects).' },
        { key: 'HuIntact', label: 'Collapse load, as built', unit: 'kN', value: intact.Hu / 1e3, explain: 'Same analysis with the as-built wall and no damage: the reference for the loss of strength.' },
        { key: 'rsr', label: 'Residual strength ratio', unit: '-', value: cap(rsrD), limit: inp.rsrReq, sense: 'min', explain: 'Collapse load as inspected over the design lateral load. The reserve the structure still has against the event it was designed for.' },
        { key: 'dsr', label: 'Damaged strength ratio', unit: '-', value: dsr, explain: 'Collapse load as inspected over collapse load as built. 1.0 means no loss; below about 0.9 the damage is structurally significant.' },
        { key: 'H1', label: 'Load at first hinge, as inspected', unit: 'kN', value: dam.H1 / 1e3, explain: 'Lateral load at which the first plastic hinge forms. Above it the dolphin takes a permanent set.' },
        { key: 'redund', label: 'Redundancy (collapse / first hinge)', unit: '-', value: dam.Hu / Math.max(dam.H1, 1), explain: 'Reserve between first yielding and collapse. Close to 1.0 the group has no redistribution and fails as soon as one pile hinges.' },
        { key: 'Ecap', label: 'Energy capacity before collapse', unit: 'kJ', value: Ecap / 1e3, limit: inp.Eacc, sense: 'min', explain: 'Deteriorated fender energy plus elastic and plastic work of the piles up to the allowable hinge rotation.' },
        { key: 'Eop', label: 'Energy capacity without permanent damage', unit: 'kJ', value: Eop / 1e3, explain: 'Deteriorated fender energy plus elastic work to first hinge: the berthing energy the dolphin can take repeatedly.' },
        { key: 'set', label: 'Predicted permanent set from the accidental impact', unit: 'mm', value: set * 1e3 + inp.y0, explain: 'Residual head displacement after the accidental event (plus any existing lean). Zero when the event is absorbed elastically.' },
        { key: 'level', label: 'Predicted damage level (0 none … 3 collapse)', unit: '-', value: level, limit: 2, sense: 'max', explain: '0 no damage, 1 local yielding, 2 plastic hinges with permanent set (repairable), 3 beyond the rotation capacity (loss of the dolphin).' },
        { key: 'Hallow', label: 'Allowable lateral / mooring load', unit: 'kN', value: Hallow / 1e3, explain: 'Collapse load as inspected divided by the required residual strength ratio: the working limit for fender reaction or line pull on the damaged dolphin.' },
        { key: 'vAllow', label: 'Allowable berthing velocity (design vessel)', unit: 'm/s', value: vAllow, explain: 'Velocity at which the abnormal berthing energy of the design vessel equals the no-damage energy capacity of the dolphin as inspected.' },
        { key: 'dispAllow', label: 'Largest vessel at 0.10 m/s', unit: 't', value: 2 * Eop / (inp.Fab * 0.01 * inp.CeCm) / 1e3, explain: 'Displacement whose abnormal berthing energy at 0.10 m/s equals the no-damage energy capacity: a size limit for the berth in its present state.' },
        { key: 'years', label: 'Years to the minimum residual strength', unit: 'years', value: Math.min(years, 200), explain: 'Time for further corrosion at the stated rate to bring the residual strength ratio down to the required value (capacity taken proportional to wall).' },
        { key: 'nEff', label: 'Effective piles', unit: '-', value: live.length, explain: 'Piles still contributing to the lateral resistance; piles flagged ineffective are removed from the model.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Pushover curves: as built and as inspected', x: { label: 'Cap displacement', unit: 'mm' }, y: { label: 'Lateral load', unit: 'kN' },
          series: [{ name: 'As built', x: intact.d.map((v) => round(v * 1e3)), y: intact.H.map((v) => round(v / 1e3)) }, { name: 'As inspected', x: dam.d.map((v) => round(v * 1e3)), y: dam.H.map((v) => round(v / 1e3)) }, { name: 'Design load', x: [0, round(dMax * 1e3)], y: [inp.Fd, inp.Fd], dash: true }],
          explain: 'Lateral load against cap displacement up to the plastic mechanism. Each kink is a hinge forming; the plateau height is the collapse load, and its ratio to the dashed design load is the residual strength ratio.' },
        { type: 'bar', title: 'Condition of each pile', categories: cond.map((c) => `Pile ${c.id}`), series: [{ name: 'Wall remaining (%)', y: cond.map((c) => round((1 - c.loss) * 100)) }, { name: 'Bending capacity remaining (%)', y: cond.map((c) => round(c.fm * (1 - c.loss) * 100)) }], y: { label: 'Share of as-built', unit: '%' },
          explain: 'Remaining wall and remaining plastic bending capacity of each pile (wall loss, dent and defect factors combined). The shortest bar is the pile to repair first.' },
        { type: 'line', title: 'Operating envelope of the dolphin as inspected', x: { label: 'Vessel displacement', unit: 't' }, y: { label: 'Allowable berthing velocity', unit: 'm/s' },
          series: [{ name: 'No permanent damage', x: ms.map((m) => round(m)), y: ms.map((m) => round(Math.sqrt(2 * Eop / (inp.Fab * m * 1e3 * inp.CeCm)))) }, { name: 'Collapse prevention (single event)', x: ms.map((m) => round(m)), y: ms.map((m) => round(Math.sqrt(2 * Ecap / (m * 1e3 * inp.CeCm)))), dash: true }],
          explain: 'Vessel size against approach velocity. Below the solid line berthings leave no damage (abnormal factor included); between the lines a single mishap is survived with permanent set; above the dashed line the dolphin is lost.' },
      ],
      tables: [
        { title: 'Pile condition and capacity factors', columns: [{ key: 'id', label: 'Pile' }, { key: 't', label: 'Wall', unit: 'mm' }, { key: 'loss', label: 'Wall loss', unit: '%' }, { key: 'dent', label: 'Dent', unit: 'mm' }, { key: 'knock', label: 'Defect factor' }, { key: 'fm', label: 'Bending factor' }, { key: 'fn', label: 'Axial factor' }, { key: 'state', label: 'State' }],
          rows: cond.map((c) => ({ id: c.id, t: round(c.t * 1e3), loss: round(c.loss * 100), dent: c.dent, knock: c.knock, fm: round(c.fm), fn: round(c.fn), state: c.knock > 0 && c.t >= 0.004 ? 'effective' : 'removed' })),
          explain: 'Inspection data as used: measured wall, dent depth and defect factor, with the resulting factors on plastic bending and axial capacity (applied on top of the reduced wall).' },
        { title: 'Hinge sequence, as inspected', columns: [{ key: 'n', label: 'Event' }, { key: 'lf', label: 'Lateral load', unit: 'kN' }, { key: 'd', label: 'Displacement', unit: 'mm' }, { key: 'pile', label: 'Pile' }, { key: 'type', label: 'Type' }],
          rows: dam.r.events.slice(0, 40).map((e, i) => ({ n: i + 1, lf: round(e.loadFactor * F / 1e3), d: round(Math.abs(e.displacement) * 1e3), pile: e.member < pilesD.length * inp.nSub ? live[Math.floor(e.member / inp.nSub)].id : 'cap', type: e.type })),
          explain: 'Order in which hinges or member failures form as the load grows: the first rows identify the weakest piles and where to strengthen.' },
      ],
      notes: [`${w.name} governs.`, `Damage prediction for ${fmt(inp.Eacc)} kJ: ${LEVELS[level]}.`, `Repair requirement: ${repair}.`, `Pushover ended with status "${dam.r.status}"; length to fixity ${fmt(Le)} m.`],
      advice,
      publish: { 'dolphin.lateralCapacity': dam.Hu, 'dolphin.residualStrengthRatio': rsrD, 'dolphin.energyCapacityDamaged': Eop },
      quantities: level >= 2 || rsrD < inp.rsrReq ? { steel_t: round(sleeveT), vessel_days: { repair: 4 } } : { vessel_days: { inspection: 1.5 } },
    };
  },
  verify: [
    // Cantilever pile, L = 20 m, undamaged: plastic collapse when the base moment reaches Mp = fy(D³ − Di³)/6
    // = 355e6·(1 − 0.96³)/6 = 6819.8 kN·m → Hu = Mp/L = 341.0 kN (negligible axial load, no P–Δ).
    { name: 'Plastic collapse of a cantilever pile: Hu = Mp/L', source: 'Plastic hinge theory (Neal, The Plastic Methods of Structural Analysis)',
      inputs: { piles: [{ x: 0, y: 0, batter: 0, azimuth: 0 }], condition: [{ tMeas: 20, dent: 0, knock: 1 }], D: 1, t: 20, steel: 's355', LeUser: 20, shear: false, pDelta: false, nSub: 4, capB: 0.5, capL: 0.5, capT: 0.2, gammaConc: 10, Wequip: 0, Fd: 100, zLoad: 0 },
      expect: [{ kpi: 'Hu', value: 341.0, tol: 0.005 }, { kpi: 'dsr', value: 1, tol: 1e-6 }] },
    // Four fixed-head piles, rigid cap, widely spaced: sway mechanism with hinges at head and base of each pile,
    // Hu = n·2Mp/L = 4·2·6819.8/20 = 2727.9 kN (push–pull axial force N/Np ≈ 0.03 lowers Mp by about 0.1 %).
    { name: 'Sway mechanism of a fixed-head group: Hu = 2·n·Mp/L', source: 'Plastic hinge theory, portal sway mechanism (Neal)',
      inputs: { piles: [{ x: -10, y: -10, batter: 0, azimuth: 0 }, { x: 10, y: -10, batter: 0, azimuth: 0 }, { x: 10, y: 10, batter: 0, azimuth: 0 }, { x: -10, y: 10, batter: 0, azimuth: 0 }], condition: [], headFix: 'fixed', D: 1, t: 20, steel: 's355', LeUser: 20, shear: false, pDelta: false, nSub: 4, capB: 0.5, capL: 0.5, capT: 0.2, gammaConc: 10, Wequip: 0, Fd: 1000, zLoad: 0 },
      expect: [{ kpi: 'Hu', value: 2727.9, tol: 0.01 }] },
    // Dented tube (thin-ring derivation in the module): δ/D = 0.1 → α = acos(0.8) = 36.87°; Mp,d/Mp = cos(18.43°) − ½ sin(36.87°) = 0.6487.
    { name: 'Dented cantilever: Hu = [cos(α/2) − ½ sin α]·Mp/L', source: 'First-principles plastic modulus of a thin ring with an ineffective arc (derived in the module source)',
      inputs: { piles: [{ x: 0, y: 0, batter: 0, azimuth: 0 }], condition: [{ tMeas: 20, dent: 100, knock: 1 }], D: 1, t: 20, steel: 's355', LeUser: 20, shear: false, pDelta: false, nSub: 4, capB: 0.5, capL: 0.5, capT: 0.2, gammaConc: 10, Wequip: 0, Fd: 100, zLoad: 0 },
      expect: [{ kpi: 'dsr', value: 0.6487, tol: 0.003 }] },
  ],
};
