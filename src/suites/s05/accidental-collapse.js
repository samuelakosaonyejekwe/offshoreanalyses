// Accidental damage and progressive collapse of a piled wharf bent: plastic pushover to collapse and the reserve
// over the design lateral load, notional removal of each pile in turn (alternate load path), residual lateral
// strength with a pile lost, and the energy a pile can absorb when struck by a vessel.
import { pushoverAnalysis, removeMembers } from '../../engine/structural/pushover.js';
import { analyzeFrame } from '../../engine/structural/frame3d.js';
import { collisionEnergy, pileImpact, tubeDentDepth } from '../../engine/structural/impact.js';
import { momentCapacity } from '../../engine/structural/concrete.js';
import { need, run, r, f, cap, check, sfCheck, eng, eco, bentInputs, bentLoadInputs, pileInputs, rcInputs, rcMember, bentModel, framePlot, steelCo2, CONFIRM, NOT_HIFI, RC_NOTE, ECON_BASIS, MGO_T_PER_DAY, CO2_MGO } from './_common.js';

/** Bent model prepared for plastic analysis: capacities on deck and piles, one gravity and one lateral load case. */
function plasticBent(inp, ctx) {
  const psi = inp.psiA;
  const B = bentModel({ ...inp, qLive: inp.qLive * psi, crane: (inp.crane || []).map((c) => ({ ...c, V: Number(c.V) * psi, H: 0 })) }, ctx);
  const rc = rcMember(inp, ctx, inp.capB, inp.capH);
  const mat = { ...rc.material, gammaC: inp.gammaCacc, gammaS: inp.gammaSacc };
  const MpS = momentCapacity({ section: rc.sag, material: mat }).MRd, MpH = momentCapacity({ section: rc.hog, material: mat }).MRd;
  const big = 1e3 * inp.capB * inp.capH * rc.fck;
  // the hinge model has one moment capacity per member: the smaller of sagging and hogging is used for the deck beam
  for (const i of B.meta.deckMembers) B.model.members[i].capacity = { Np: big, Mpy: Math.min(MpS, MpH), Mpz: big, Pc: big, Pt: big };
  for (const i of B.meta.pileMembers) B.model.members[i].K = inp.headFix === 'pinned' ? 2.1 : 1.2;
  const lc = Object.fromEntries(B.model.loadCases.map((c) => [c.name, c]));
  const Href = inp.Hd * 1e3;
  B.model.loadCases = [
    { name: 'grav', selfWeight: lc.G.selfWeight, member: [...(lc.G.member || []), ...(lc.Q.member || []), ...(lc.C.member || [])], nodal: [...(lc.C.nodal || [])] },
    { name: 'lat', nodal: [{ node: B.headSea, F: [Href, 0, 0, 0, 0, 0] }] },
  ];
  return { B, rc, MpDeck: Math.min(MpS, MpH), Href };
}

export default {
  id: 's05-accidental-collapse', suite: 5, title: 'Accidental damage and progressive collapse',
  summary: 'Collapse load of the bent and its reserve over the accidental design load, whether the deck survives the loss of any one pile, the lateral strength that remains, and whether a pile struck by a vessel absorbs the impact energy.',
  covers: [14], fidelity: 'numerical',
  standards: ['EN 1991-1-7:2006 §3.3, Annex A (accidental actions, notional removal)', 'ISO 19902:2020 §10, §12 (accidental situations, reserve strength)', 'DNV-RP-C204 (2019) Structural design against accidental loads', 'UFC 4-023-03 (2016) alternate path method', 'BS 6349-2:2019 §6 (accidental berthing)', 'ASCE 61-14 §3 (pushover of piers and wharves)'],
  models: ['Nonlinear frame analysis: event-to-event plastic hinge pushover', 'Member removal benchmarks (alternate load path)', 'Column buckling with post-buckling residual strength', 'P–Δ geometric stiffness', 'Impact energy balance: elastic cantilever plus plastic hinge rotation', 'Local denting of a tube (plastic mechanism)', 'Momentum transfer (added-mass kinetic energy)'],
  limits: `Static plastic-hinge analysis of one bent in its plane: dynamic effects of sudden member loss enter through an amplification factor, and the deck spanning along the wharf to neighbouring bents after a pile is lost enters only through the load share entered. Soil is represented by the fixity depth. Explicit nonlinear contact finite-element simulation of the collision is outside the module. ${RC_NOTE} ${NOT_HIFI}`,
  publishes: ['wharf.RSR', 'wharf.collapseLoad', 'wharf.damagedStrengthRatio'],
  searchTerms: 'wharf progressive collapse pushover pile loss vessel collision accidental design',
  inputs: [
    ...bentInputs(), ...pileInputs(), ...bentLoadInputs().filter((q) => ['gSdl', 'qLive', 'crane'].includes(q.key)), ...rcInputs().filter((q) => !['gammaC', 'gammaS'].includes(q.key)),
    { key: 'gammaC', label: 'γc concrete (persistent situation, for bar fit only)', unit: '-', type: 'number', default: 1.5, min: 1, max: 2, group: 'Accidental criteria', role: 'criteria', help: CONFIRM },
    { key: 'gammaS', label: 'γs reinforcement (persistent situation)', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.5, group: 'Accidental criteria', role: 'criteria', help: CONFIRM },
    { key: 'gammaCacc', label: 'γc concrete, accidental situation', unit: '-', type: 'number', default: 1.2, min: 1, max: 1.5, group: 'Accidental criteria', role: 'criteria', help: `EN 1992-1-1 Table 2.1N. ${CONFIRM}` },
    { key: 'gammaSacc', label: 'γs reinforcement, accidental situation', unit: '-', type: 'number', default: 1.0, min: 1, max: 1.15, group: 'Accidental criteria', role: 'criteria', help: CONFIRM },
    { key: 'psiA', label: 'Share of cargo and crane load present in the accident', unit: '-', type: 'number', default: 0.5, min: 0, max: 1, group: 'Accidental criteria', role: 'criteria', help: `Frequent or quasi-permanent combination factor of EN 1990 §6.4.3.3. ${CONFIRM}` },
    { key: 'Hd', label: 'Accidental design lateral load on the bent', unit: 'kN', type: 'number', default: 900, min: 10, max: 50000, group: 'Accidental criteria', role: 'bc', bind: 'wharf.fenderReactionBent', help: 'For example the abnormal berthing reaction shared to one bent, or the impact force of an out-of-control vessel.' },
    { key: 'rsrMin', label: 'Required reserve strength ratio', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Accidental criteria', role: 'criteria', help: `Collapse load over the accidental design load; ISO 19902 quotes 1.85 for manned jackets under the 100-year environment, port structures typically 1.25–1.6 on the accidental load. ${CONFIRM}` },
    { key: 'daf', label: 'Dynamic amplification on gravity loads after sudden pile loss', unit: '-', type: 'number', default: 1.5, min: 1, max: 2, group: 'Accidental criteria', role: 'criteria', help: `2.0 is the upper bound for instantaneous removal of an elastic structure; about 1.3–1.5 where plastic redistribution is counted, as here; 1.0 for a pile lost gradually (corrosion) or removed for repair. ${CONFIRM}` },
    { key: 'kLong', label: 'Share of a lost pile reaction carried along the wharf to neighbouring bents', unit: '-', type: 'number', default: 0.8, min: 0, max: 0.95, group: 'Accidental criteria', role: 'bc', help: 'The rail beam and slab span two bent spacings over a lost pile. 0 treats the bent alone (no longitudinal continuity, for example at an expansion joint); 0.6–0.85 for a continuous rail beam and slab detailed as a tie. Confirm with a longitudinal analysis of the rail beam over the double span.' },
    { key: 'pDelta', label: 'Include P–Δ effects', type: 'bool', default: true, group: 'Accidental criteria', role: 'control' },
    { key: 'postBuck', label: 'Residual compression after buckling (fraction)', unit: '-', type: 'number', default: 0.4, min: 0, max: 1, group: 'Accidental criteria', role: 'criteria' },
    { key: 'mVessel', label: 'Mass of the striking vessel', unit: 't', type: 'number', default: 1500, min: 10, max: 300000, group: 'Vessel impact on a pile', role: 'ic', help: 'A tug or workboat striking an exposed pile is the usual accidental case; large ships are kept off the piles by the fender line.' },
    { key: 'vVessel', label: 'Impact speed', unit: 'm/s', type: 'number', default: 0.8, min: 0, max: 6, group: 'Vessel impact on a pile', role: 'ic', uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'CmA', label: 'Added-mass factor', unit: '-', type: 'number', default: 1.1, min: 1, max: 1.8, group: 'Vessel impact on a pile', role: 'criteria', help: '1.1 bow or stern impact, 1.4 broadside (DNV-RP-C204).' },
    { key: 'shareE', label: 'Share of the energy absorbed by the pile', unit: '-', type: 'number', default: 0.5, min: 0.05, max: 1, group: 'Vessel impact on a pile', role: 'criteria', help: 'The rest is absorbed by the vessel hull and by other piles through the deck; 1.0 is the strength-design upper bound.' },
    { key: 'zImpact', label: 'Impact level above datum', unit: 'm', type: 'number', default: 1.0, min: -10, max: 10, group: 'Vessel impact on a pile', role: 'bc' },
    { key: 'thetaCap', label: 'Plastic rotation capacity of the pile hinge', unit: 'rad', type: 'number', default: 0.03, min: 0, max: 0.15, group: 'Vessel impact on a pile', role: 'criteria', help: 'Depends on D/t and steel grade: about 0.02–0.05 rad for compact tubes (D/t below 9000/fy), near zero for slender ones.' },
  ],
  compute(inp, ctx) {
    const { B, MpDeck, Href } = plasticBent(inp, ctx), model = B.model;
    const opts = { control: { node: B.headSea, dof: 0 }, pDelta: !!inp.pDelta, postBucklingRatio: inp.postBuck, maxEvents: 80 };
    const push = run(() => pushoverAnalysis(model, { pattern: 'lat', base: 'grav', ...opts }), 'Pushover analysis');
    need(push.status !== 'collapse-under-base-load', 'The bent collapses under its gravity loads alone in the accidental combination: check the deck beam and pile sizes before assessing accidental actions.');
    const Hu = push.ultimateLoadFactor * Href, H1 = push.firstEventLoadFactor * Href, RSR = push.ultimateLoadFactor;
    // notional removal of each pile member in turn under the amplified gravity loads. The rail beam and slab span along
    // the wharf to the neighbouring bents: the share entered of the lost pile's gravity reaction is taken there, applied
    // as an upward force at the head of the removed pile that scales with the gravity load.
    const piles = B.meta.pileMembers;
    const grav0 = run(() => analyzeFrame(model), 'Frame analysis').results.grav;
    const damagedModel = (i) => {
      const m = removeMembers(model, [i]).model, head = model.members[i].b, Nc = Math.max(-grav0.memberForces[i].endB.N, 0);
      const dz = Math.abs(model.nodes[model.members[i].b][2] - model.nodes[model.members[i].a][2]), len = Math.hypot(...[0, 1, 2].map((q) => model.nodes[model.members[i].b][q] - model.nodes[model.members[i].a][q]));
      m.loadCases = m.loadCases.map((c) => (c.name === 'grav' ? { ...c, nodal: [...(c.nodal || []), { node: head, F: [0, 0, inp.kLong * Nc * dz / len, 0, 0, 0] }] } : c));
      return m;
    };
    const cases = piles.length > 1 ? piles.map((i) => {
      const p = run(() => pushoverAnalysis(damagedModel(i), { pattern: 'grav', pDelta: !!inp.pDelta, postBucklingRatio: inp.postBuck, maxEvents: 80, maxLoadFactor: 20 }), 'Member removal analysis');
      return { member: i, ultimateLoadFactor: p.ultimateLoadFactor, uc: inp.daf / Math.max(p.ultimateLoadFactor, 1e-6), survives: p.ultimateLoadFactor >= inp.daf, status: p.status };
    }) : [];
    const pc = { cases };
    const crit = pc.cases.length ? pc.cases.reduce((m, c) => (c.uc > m.uc ? c : m), pc.cases[0]) : null;
    // lateral strength with the critical pile removed (if the bent still stands under gravity)
    let ds = null;
    if (crit && crit.ultimateLoadFactor > 1) { try { const d = pushoverAnalysis(damagedModel(crit.member), { pattern: 'lat', base: 'grav', ...opts }); if (d.status !== 'collapse-under-base-load') ds = { damaged: d, damagedStrengthRatio: d.ultimateLoadFactor / push.ultimateLoadFactor }; } catch { ds = null; } }
    const DSR = ds ? ds.damagedStrengthRatio : 0, HuDam = ds ? ds.damaged.ultimateLoadFactor * Href : 0;
    // energy absorbed by the bent up to collapse: ∫ H dδ along the pushover curve
    const Hc = push.curve.loadFactor.map((l) => l * Href), dc = push.curve.displacement.map(Math.abs); let Ebent = 0; for (let i = 1; i < Hc.length; i++) Ebent += 0.5 * (Hc[i] + Hc[i - 1]) * (dc[i] - dc[i - 1]);
    // vessel impact on the seaward pile
    const P = B.P, Ek = collisionEnergy({ mass: inp.mVessel * 1e3, velocity: inp.vVessel, addedMassFactor: inp.CmA }).energy, Epile = inp.shareE * Ek;
    const hImp = inp.zImpact - (B.meta.seabedZ ? B.meta.seabedZ[0] : -inp.depth) + (B.meta.fix ?? 5 * P.D);
    need(hImp > 1, 'The impact level must be at least 1 m above the point of fixity of the seaward pile.');
    const imp = pileImpact({ energy: Epile, Mp: P.Mp, EI: P.EI, height: hImp, rotationCapacity: inp.thetaCap });
    const dent = tubeDentDepth({ D: P.D, t: P.t, fy: P.fy, force: imp.force });
    const vPerm = Math.sqrt(2 * imp.capacity / (inp.shareE * inp.mVessel * 1e3 * inp.CmA));
    const names = (i) => { const k = piles.indexOf(i); return k >= 0 ? `pile ${k + 1}${model.members[i].group === 'raker' ? ' (raker)' : ''}` : `deck beam member ${i}`; };
    const checks = [
      sfCheck('Reserve strength under the accidental lateral load', RSR, inp.rsrMin, 'ISO 19902:2020 §12.6 (RSR); BS 6349-2:2019 §6', 'Collapse load from the pushover over the accidental design load. Below the required ratio a plastic mechanism forms in the design accident and the deck is lost.'),
      ...(crit ? [check(`Survival after loss of ${names(crit.member)}`, inp.daf, cap(crit.ultimateLoadFactor), '× gravity load', 'EN 1991-1-7:2006 Annex A; UFC 4-023-03 (2016) §3', 'Gravity load factor the bent can carry with one pile removed, against the amplified load at the moment of loss. Below it the deck beam hinges over the gap and collapse spreads to the next bay.')] : []),
      check('Vessel impact energy on the seaward pile', Epile / 1e3, imp.capacity / 1e3, 'kJ', 'DNV-RP-C204 (2019) §3; EN 1991-1-7:2006 §4.6', 'Energy to be absorbed by the pile against its elastic energy plus the plastic work of the hinge. If exceeded the hinge runs out of rotation capacity: the pile tears or buckles locally and stops carrying the deck.'),
      check('Local dent of the pile wall (depth / diameter)', dent.dentRatio, 0.1, '-', 'DNV-RP-C204 (2019) §3.6; API RP 2A-WSD §C17', 'Dent depth under the contact force over the diameter. Beyond about 10 % the bending resistance of the tube falls quickly and the hinge capacity assumed above is not available.'),
    ];
    const w = checks.reduce((m, c) => (c.uc > m.uc ? c : m)), ok = w.uc <= 1;
    const fracs = pc.cases.map((c) => r(Math.min(c.ultimateLoadFactor, 20), 3));
    const tNeed = inp.corr + (inp.t - inp.corr) * Math.max(inp.rsrMin / Math.max(RSR, 0.05), 0.2), steelBent = B.pileLen * P.massPerM / 1e3;
    const evt = push.events.slice(0, 20).map((e, i) => ({ n: i + 1, H: r(e.loadFactor * Href / 1e3), member: names(e.member), type: e.type, d: r(Math.abs(e.displacement) * 1e3) }));
    const hingeVal = model.nodes.map(() => 0); push.events.forEach((e, i) => { const m = model.members[e.member]; if (m) { hingeVal[m.a] = Math.max(hingeVal[m.a], i + 1); hingeVal[m.b] = Math.max(hingeVal[m.b], i + 1); } });
    let dmax = 0; for (let i = 0; i < model.nodes.length; i++) dmax = Math.max(dmax, Math.hypot(push.displacements[6 * i], push.displacements[6 * i + 2]));
    const xs = model.nodes.map((p) => p[0]), zs = model.nodes.map((p) => p[2]), size = Math.max(Math.max(...xs) - Math.min(...xs), Math.max(...zs) - Math.min(...zs)), sc = dmax > 0 ? r(0.08 * size / dmax, 2) : 1;
    return {
      kpis: [
        { key: 'Hu', label: 'Collapse lateral load of the bent', value: r(Hu / 1e3), unit: 'kN', explain: `Lateral load at which a plastic mechanism forms (${push.status}), with gravity loads held. It is the largest berthing or impact force the bent can resist once.` },
        { key: 'H1', label: 'Lateral load at first yield', value: r(H1 / 1e3), unit: 'kN', explain: 'Load at the first plastic hinge or member buckling. Up to here the bent is undamaged; beyond it permanent deformation remains.' },
        { key: 'RSR', label: 'Reserve strength ratio', value: r(RSR), unit: '-', limit: inp.rsrMin, sense: 'min', explain: 'Collapse load over the accidental design lateral load. It measures how far the design accident is from collapse.' },
        { key: 'redundancy', label: 'Redundancy factor', value: r(Number.isFinite(push.redundancyFactor) ? push.redundancyFactor : 1), unit: '-', explain: 'Collapse load over first-yield load. Close to 1 means failure follows first yield almost at once; higher means the bent redistributes load between piles.' },
        { key: 'dU', label: 'Deck displacement at collapse', value: r(dc[dc.length - 1] * 1e3), unit: 'mm', explain: 'Sway at the last equilibrium state of the pushover: the deformation the deck, crane rails and services would see in a near-collapse event.' },
        { key: 'Ebent', label: 'Energy absorbed by the bent up to collapse', value: r(Ebent / 1e3), unit: 'kJ', explain: 'Area under the load–displacement curve: the impact energy one bent can dissipate through yielding of its piles.' },
        { key: 'lamRem', label: 'Gravity load factor with the critical pile removed', value: r(crit ? Math.min(crit.ultimateLoadFactor, 20) : 20), unit: '-', limit: inp.daf, sense: 'min', explain: `Multiple of the accidental gravity combination the bent carries after losing ${crit ? names(crit.member) : 'a pile'} (capped at 20). It must exceed the dynamic amplification factor for a sudden loss.` },
        { key: 'DSR', label: 'Damaged strength ratio (lateral)', value: r(DSR), unit: '-', explain: 'Lateral collapse load with the critical pile removed over the intact value. It tells how much berthing load the damaged bent can still take; zero means it cannot stand at all.' },
        { key: 'HuDam', label: 'Maximum permissible lateral load with one pile lost', value: r(HuDam / inp.rsrMin / 1e3), unit: 'kN', explain: 'Damaged collapse load over the required reserve ratio: the berthing or mooring load limit to impose on this bent until it is repaired.' },
        { key: 'Eimpact', label: 'Kinetic energy of the striking vessel', value: r(Ek / 1e3), unit: 'kJ', explain: '½·m·Cm·v² of the vessel with its added mass. The pile is checked for the share entered.' },
        { key: 'Eel', label: 'Elastic energy capacity of the pile', value: r(imp.elasticEnergy / 1e3), unit: 'kJ', explain: 'Energy the pile absorbs as a cantilever before a hinge forms at fixity, Mp²h/(6EI). Impacts below it leave no permanent damage.' },
        { key: 'Ecap', label: 'Total energy capacity of the pile', value: r(imp.capacity / 1e3), unit: 'kJ', limit: r(Epile / 1e3), sense: 'min', explain: 'Elastic energy plus plastic work Mp·θ of the hinge up to its rotation capacity. It must not be less than the energy the pile has to absorb.' },
        { key: 'Fimpact', label: 'Contact force on the pile', value: r(imp.force / 1e3), unit: 'kN', explain: `Force between vessel and pile during the impact; it produces a dent of ${f(dent.dent * 1e3)} mm and is passed into the deck through the pile head.` },
        { key: 'vPerm', label: 'Maximum permissible impact speed', value: r(vPerm), unit: 'm/s', explain: 'Speed of this vessel at which the pile just absorbs its share of the energy: a limit for workboat operations around unprotected piles.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Pushover curve of the bent', x: { label: 'Deck displacement', unit: 'mm' }, y: { label: 'Lateral load', unit: 'kN' },
          series: [{ name: 'Intact', x: dc.map((v) => r(v * 1e3)), y: Hc.map((v) => r(v / 1e3)) }, ...(ds ? [{ name: `Without ${names(crit.member)}`, x: ds.damaged.curve.displacement.map((v) => r(Math.abs(v) * 1e3)), y: ds.damaged.curve.loadFactor.map((l) => r(l * Href / 1e3)), dash: true }] : []), { name: 'Accidental design load', x: [0, r(dc[dc.length - 1] * 1e3)], y: [inp.Hd, inp.Hd], dash: true }],
          explain: 'Lateral load against deck sway as hinges form one after another. Each kink is a member yielding or buckling; the plateau is the collapse load. The dashed curve is the bent with one pile lost.' },
        framePlot(model, { disp: push.displacements, scale: sc, title: `Deformed shape at collapse (× ${sc})`, values: hingeVal, valueLabel: 'Order of yielding at adjacent members', valueUnit: '-',
          explain: 'The bent at the end of the pushover. Coloured nodes mark the members that yielded, in sequence: the pattern shows the collapse mechanism and where ductile detailing matters.' }),
        ...(pc.cases.length ? [{ type: 'bar', title: 'Gravity load factor after removal of each pile', categories: pc.cases.map((c) => names(c.member)), series: [{ name: 'Ultimate gravity load factor', y: fracs }, { name: 'Required (amplification factor)', y: pc.cases.map(() => inp.daf) }], x: { label: 'Removed member', unit: '' }, y: { label: 'Load factor', unit: '-' },
          explain: 'How many times the accidental gravity combination the bent carries with each pile removed in turn. Bars below the required factor identify key elements whose loss leads to collapse.' }] : []),
      ],
      tables: [
        { title: 'Sequence of yielding in the pushover', columns: [{ key: 'n', label: 'Event' }, { key: 'H', label: 'Lateral load', unit: 'kN' }, { key: 'member', label: 'Member' }, { key: 'type', label: 'Type' }, { key: 'd', label: 'Deck displacement', unit: 'mm' }],
          rows: evt, explain: 'Each plastic event in order with the load and sway at which it occurs. The first rows locate the damage to look for after an overload; the last ones complete the mechanism.' },
        ...(pc.cases.length ? [{ title: 'Notional member removal', columns: [{ key: 'member', label: 'Removed' }, { key: 'lam', label: 'Gravity load factor' }, { key: 'uc', label: 'Utilisation' }, { key: 'res', label: 'Result' }],
          rows: pc.cases.map((c) => ({ member: names(c.member), lam: r(Math.min(c.ultimateLoadFactor, 20), 3), uc: r(cap(c.uc), 3), res: c.survives ? 'alternate load path holds' : 'collapse spreads' })), explain: 'Alternate-load-path screening: one row per pile removed. A pile whose removal fails the check is a key element and should be protected or designed for the accidental action itself.' }] : []),
      ],
      notes: [`Pushover ended with status '${push.status}' after ${push.events.length} plastic event(s); first event in ${push.events.length ? names(push.events[0].member) : '–'}. Deck beam hinge capacity ${f(MpDeck / 1e3)} kN·m (accidental factors).`,
        `Impact: ${f(Epile / 1e3)} kJ on the pile at ${f(hImp)} m above fixity; the pile ${imp.elastic ? 'stays elastic' : `forms a hinge and rotates ${f(imp.plasticRotation)} rad`}.`],
      advice: [
        ok ? eng('ok', `Reserve strength ratio ${f(RSR)} (collapse at ${f(Hu / 1e3)} kN, first yield at ${f(H1 / 1e3)} kN); ${crit ? `the bent carries ${f(Math.min(crit.ultimateLoadFactor, 20))} times gravity with ${names(crit.member)} removed; ` : ''}pile impact uses ${f(checks[checks.length - 2].uc * 100)} % of the energy capacity.`,
          `Record ${f(H1 / 1e3)} kN as the load above which the bent must be inspected after a hard berthing, and ${f(HuDam / inp.rsrMin / 1e3)} kN as the operating limit if a pile is lost. Detail pile heads for the hinge rotations the mechanism needs.`, 'EN 1991-1-7:2006 §3.3; ISO 19902:2020 §12; DNV-RP-C204 (2019)',
          `Robustness here comes from redundancy rather than added material: the ${f(steelBent)} t of pile steel per bent already provides an alternate load path, avoiding sacrificial protection piles (${steelCo2(P.massPerM * 30 / 1e3)} for each 30 m protection pile avoided).`)
          : eng('fail', `${w.name}: utilisation ${f(w.uc)}. Reserve strength ratio ${f(RSR)} against ${f(inp.rsrMin)}; ${crit ? `gravity load factor ${f(Math.min(crit.ultimateLoadFactor, 20))} with ${names(crit.member)} removed against ${f(inp.daf)}.` : ''}`,
            w.name.startsWith('Reserve') ? `Increase the pile wall to about ${f(Math.ceil(tNeed))} mm or add raker piles; alternatively limit the accidental lateral load to ${f(Hu / inp.rsrMin / 1e3)} kN by fender selection.` : w.name.startsWith('Survival') ? `Strengthen the deck beam to span the gap (hinge capacity of about ${f(MpDeck * w.uc / 1e3)} kN·m), tie the deck along the wharf so that neighbouring bents share the load, or protect that pile as a key element.` : `Keep vessels off the piles (fender or protection pile line), limit the speed of workboats to ${f(vPerm)} m/s, or use a thicker, more compact pile section in the impact zone.`,
            'EN 1991-1-7:2006 §3.3, Annex A; DNV-RP-C204 (2019) §3; UFC 4-023-03 (2016)',
            `Thicker piles to ${f(Math.ceil(tNeed))} mm: ${steelCo2(steelBent * (Math.ceil(tNeed) - inp.t) / inp.t)} per bent. An operating speed limit or a fender line change adds no structural material.`),
        eco(ok ? 'info' : 'warn', `A collapse of one bent would close the berth for months; the damaged bent can still take ${f(HuDam / 1e3)} kN laterally (${f(DSR * 100)} % of intact), which decides whether the berth stays open during repair.`,
          `Price the accident, not only the structure: compare the cost of robustness (redundant piles, deck ties, protection) with the expected loss of berth revenue and the repair spread, about ${f(MGO_T_PER_DAY * CO2_MGO)} t CO₂ per marine-spread day (screening). Hold an emergency response plan with the inspection trigger of ${f(H1 / 1e3)} kN.`, `Risk assessment per EN 1991-1-7 Annex B and ISO 2394; ${ECON_BASIS}`),
      ],
      publish: { 'wharf.RSR': RSR, 'wharf.collapseLoad': Hu, 'wharf.damagedStrengthRatio': DSR },
    };
  },
  verify: [
    // Two vertical piles with pinned heads under a stiff deck beam, level bed: each pile is a cantilever built in at fixity, so the
    // mechanism needs one hinge at the base of each pile: Hu = 2·Mp/L. Tube 1.0 m × 20 mm, S355: Mp = 355e6·(1 − 0.96³)/6 = 6819.8 kN·m;
    // L = 4 + 12 + 5 = 21 m → Hu = 649.5 kN. First yield coincides (both piles reach Mp together): redundancy 1.0. P–Δ switched off.
    { name: 'Sway mechanism of two pinned-head piles: Hu = 2 Mp / L', source: 'Plastic collapse by the mechanism method (Neal, Plastic Methods of Structural Analysis)',
      inputs: { nPiles: 2, D: 1, t: 20, corr: 0, bedRise: 0, headFix: 'pinned', fixDepth: 5, deckLevel: 4, depth: 12, qLive: 0, gSdl: 0, crane: [], pDelta: false, Hd: 500, kLong: 0.9 },
      expect: [{ kpi: 'Hu', value: 649.5, tol: 0.01 }, { kpi: 'RSR', value: 1.299, tol: 0.01 }] },
    // Vessel of 2000 t at 1.5 m/s broadside (Cm = 1.4): E = ½·2e6·1.4·1.5² = 3150 kJ.
    // Pile as above struck 21 m above fixity: elastic energy Mp²·h/(6EI) = (6.8198e6)²·21/(6·1.5530e9) = 104.8 kJ; with θ = 0.03 rad the capacity is 104.8 + 6819.8·0.03 = 309.4 kJ.
    { name: 'Collision energy ½ m Cm v² and pile energy capacity Mp²h/(6EI) + Mp θ', source: 'DNV-RP-C204 (2019) §3.3; elastic cantilever and plastic hinge work (hand calculation)',
      inputs: { nPiles: 2, D: 1, t: 20, corr: 0, bedRise: 0, fixDepth: 5, depth: 12, zImpact: 4, mVessel: 2000, vVessel: 1.5, CmA: 1.4, thetaCap: 0.03, crane: [] },
      expect: [{ kpi: 'Eimpact', value: 3150, tol: 1e-6 }, { kpi: 'Eel', value: 104.8, tol: 0.003 }, { kpi: 'Ecap', value: 309.4, tol: 0.003 }] },
  ],
};
