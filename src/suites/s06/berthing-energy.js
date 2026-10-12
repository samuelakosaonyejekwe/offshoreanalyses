// Berthing impact energy by the kinetic-energy method (PIANC WG33 2002 / WG211 2024 practice, BS 6349-4).
import { blockCoefficient, radiusOfGyration, eccentricityFactor, contactGeometry, addedMassFactor, berthingEnergy, velocityFromTable, exceedanceLognormal, normCdf, VELOCITY_TABLE_DEFAULT } from './_berthing.js';
import { DEG, need, fmt, round, adv, check, linspace, CONFIRM, NOT_HIFI, ECON_BASIS } from './_common.js';

function evaluate(inp, over = {}) {
  const p = { ...inp, ...over };
  // SI: t → kg, deg → rad
  const M = p.disp * 1e3, L = p.LBP, B = p.beam, D = p.draft;
  need(p.depth > D, `The water depth (${fmt(p.depth)} m) must exceed the vessel draft (${fmt(D)} m): the vessel is aground.`);
  const Cb = blockCoefficient({ M, L, B, D, rho: p.rhoW });
  need(Cb > 0.3 && Cb <= 1.0, `Displacement, length, beam and draft give a block coefficient of ${fmt(Cb)}; a ship hull lies between about 0.4 and 0.9. Check the vessel data.`);
  const v = p.velMode === 'table' ? velocityFromTable(p.velTable, p.dwt, p.cond) : p.vUser;
  need(v > 0, 'The berthing velocity must be positive.');
  const K = radiusOfGyration({ Cb, L });
  const { R, gamma } = contactGeometry({ L, B, x: p.xFrac * L, alpha: p.alpha * DEG });
  // near-midships contact at an angle gives a small negative γ from the formula; Ce depends on cos²γ only
  const Ce = eccentricityFactor({ K, R, gamma });
  const Kc = p.depth - D;
  const Cm = addedMassFactor({ method: p.cmMethod, D, B, Cb, Kc, user: p.CmUser });
  const e = berthingEnergy({ M, v, Ce, Cm, Cs: p.Cs, Cc: p.Cc, Fab: p.Fab });
  return { M, v, Cb, K, R, gamma, Ce, Cm, Kc, ...e, E0: 0.5 * M * v * v };
}

export default {
  id: 's06-berthing-energy', suite: 6, title: 'Berthing impact energy',
  summary: 'Normal and abnormal berthing energy a dolphin must absorb, from vessel mass, approach velocity, eccentricity, added mass, softness and berth configuration.',
  covers: [1], fidelity: 'analytical',
  standards: ['PIANC WG33 Guidelines for the Design of Fender Systems (2002) §4.2', 'PIANC WG211 (2024)', 'BS 6349-4:2014 §5', 'ROM 2.0-11'],
  models: ['Vessel kinetic energy', 'Conservation of angular momentum about the contact point (eccentricity factor)', 'Added mass: Vasco Costa, Ueda, PIANC under-keel clearance', 'Softness and berth-configuration factors', 'Abnormal berthing factor', 'Lognormal berthing-velocity exceedance'],
  limits: 'Rigid vessel, single point of contact arrested without sliding, velocity normal to the berthing line. The velocity guidance table and the abnormal-berthing factor are indicative defaults to be confirmed with the PIANC edition in force and with port records. Coupled vessel–fender time-domain dynamics and tug or wave-induced approach are not modelled. ' + NOT_HIFI,
  publishes: ['dolphin.berthEnergyNormal', 'dolphin.berthEnergyAbnormal', 'dolphin.berthVelocity'],
  searchTerms: 'berthing energy fender design eccentricity added mass PIANC',
  inputs: [
    { key: 'disp', label: 'Vessel displacement at berthing', unit: 't', type: 'number', default: 75000, min: 100, max: 700000, group: 'Vessel', role: 'load', bind: 'vessel.disp', uncertain: { dist: 'normal', cov: 0.05 },
      help: 'Displacement (not deadweight) in the berthing condition: laden for import berths, ballast for export.' },
    { key: 'dwt', label: 'Deadweight (for the velocity table)', unit: 't', type: 'number', default: 60000, min: 100, max: 600000, group: 'Vessel', role: 'load' },
    { key: 'LBP', label: 'Length between perpendiculars', unit: 'm', type: 'number', default: 220, min: 10, max: 420, group: 'Vessel', role: 'geometry' },
    { key: 'beam', label: 'Beam', unit: 'm', type: 'number', default: 32.2, min: 3, max: 70, group: 'Vessel', role: 'geometry' },
    { key: 'draft', label: 'Draft at berthing', unit: 'm', type: 'number', default: 12.5, min: 1, max: 26, group: 'Vessel', role: 'geometry' },
    { key: 'velMode', label: 'Berthing velocity from', type: 'select', default: 'table', group: 'Approach', role: 'ic',
      options: [{ value: 'table', label: 'Guidance table (by deadweight and condition)' }, { value: 'user', label: 'Specified value' }] },
    { key: 'cond', label: 'Navigation condition', type: 'select', default: 'b', group: 'Approach', role: 'ic', showIf: { key: 'velMode', in: ['table'] },
      options: [{ value: 'a', label: 'a — easy berthing, sheltered' }, { value: 'b', label: 'b — difficult berthing, sheltered' }, { value: 'c', label: 'c — easy berthing, exposed' }, { value: 'd', label: 'd — good berthing, exposed' }, { value: 'e', label: 'e — difficult berthing, exposed' }] },
    { key: 'velTable', label: 'Berthing-velocity guidance (indicative)', type: 'table', group: 'Approach', role: 'ic', showIf: { key: 'velMode', in: ['table'] },
      columns: [{ key: 'dwt', label: 'Deadweight', unit: 't', type: 'number' }, { key: 'a', label: 'a', unit: 'm/s', type: 'number' }, { key: 'b', label: 'b', unit: 'm/s', type: 'number' }, { key: 'c', label: 'c', unit: 'm/s', type: 'number' }, { key: 'd', label: 'd', unit: 'm/s', type: 'number' }, { key: 'e', label: 'e', unit: 'm/s', type: 'number' }],
      default: VELOCITY_TABLE_DEFAULT, help: 'Indicative Brolsma-type values as reproduced in fender design manuals after PIANC WG33 (2002). Replace with the curves of the PIANC edition in force or with measured berthing records of the port.' },
    { key: 'vUser', label: 'Berthing velocity normal to the berth', unit: 'm/s', type: 'number', default: 0.15, min: 0.01, max: 1.5, group: 'Approach', role: 'ic', showIf: { key: 'velMode', in: ['user'] }, uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'alpha', label: 'Berthing angle', unit: 'deg', type: 'number', default: 6, min: 0, max: 25, group: 'Approach', role: 'ic', uncertain: { dist: 'normal', cov: 0.3 }, help: 'Angle between hull and berthing line at first contact; 5–6° for large tankers on dolphins, up to 10–15° for small vessels.' },
    { key: 'xFrac', label: 'Contact point from the bow / LBP', unit: '-', type: 'number', default: 0.25, min: 0.05, max: 0.5, group: 'Approach', role: 'ic', help: '0.25 quarter-point berthing (typical for breasting dolphins), 0.33 third-point, 0.5 midships (Ce → 1).' },
    { key: 'depth', label: 'Water depth at the berth', unit: 'm', type: 'number', default: 15, min: 1.5, max: 100, group: 'Site', role: 'environment', bind: 'site.depth' },
    { key: 'rhoW', label: 'Water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'cmMethod', label: 'Added-mass factor', type: 'select', default: 'vasco-costa', group: 'Factors', role: 'criteria',
      options: [{ value: 'vasco-costa', label: 'Vasco Costa: 1 + 2D/B' }, { value: 'ueda', label: 'Ueda: 1 + πD/(2·Cb·B)' }, { value: 'pianc', label: 'PIANC 2002 by under-keel clearance' }, { value: 'user', label: 'Specified value' }] },
    { key: 'CmUser', label: 'Added-mass factor Cm', unit: '-', type: 'number', default: 1.8, min: 1, max: 3, group: 'Factors', role: 'criteria', showIf: { key: 'cmMethod', in: ['user'] } },
    { key: 'Cs', label: 'Softness factor Cs', unit: '-', type: 'number', default: 1.0, min: 0.8, max: 1, group: 'Factors', role: 'criteria', help: '1.0 for soft fenders (deflection above about 150 mm), 0.9 for hard fenders where the hull itself deforms.' },
    { key: 'Cc', label: 'Berth configuration factor Cc', unit: '-', type: 'number', default: 1.0, min: 0.7, max: 1, group: 'Factors', role: 'criteria', help: '1.0 for open piled structures such as dolphins; 0.8–0.9 for closed quay walls with parallel approach (water cushion).' },
    { key: 'Fab', label: 'Abnormal berthing factor', unit: '-', type: 'number', default: 1.25, min: 1, max: 3, group: 'Factors', role: 'criteria',
      help: `PIANC WG33 Table 4.2.5: 1.25 largest tankers and bulk carriers, 1.75 smallest; containers 1.5–2.0; general cargo 1.75; ferries and tugs 2.0 or more. WG211 (2024) moves to partial factors on the characteristic energy. ${CONFIRM}` },
    { key: 'Ecap', label: 'Energy capacity of fender and dolphin', unit: 'kJ', type: 'number', default: 800, min: 1, max: 1e5, group: 'Acceptance', role: 'criteria', bind: 'dolphin.energyCapacity', design: true, help: 'Minimum-tolerance energy the fender system (and flexible piles, if any) can absorb; published by the fender and flexible-dolphin modules.' },
    { key: 'vCov', label: 'Coefficient of variation of berthing velocity', unit: '-', type: 'number', default: 0.35, min: 0.05, max: 1, group: 'Acceptance', role: 'criteria', help: 'Scatter of measured berthing velocities for this vessel class; 0.3–0.5 in published port records.' },
    { key: 'vFractile', label: 'Design velocity is the fractile', unit: '-', type: 'number', default: 0.98, min: 0.5, max: 0.9999, group: 'Acceptance', role: 'criteria', help: 'Non-exceedance probability per berthing attributed to the design velocity (the guidance curves are commonly read as about a 1-in-50 event).' },
  ],
  compute(inp) {
    const r = evaluate(inp);
    const EN = r.EN / 1e3, EA = r.EA / 1e3;
    // Normal score of the design fractile by bisection on the distribution function.
    let lo = -6, hi = 6; for (let i = 0; i < 60; i++) { const m = 0.5 * (lo + hi); if (normCdf(m) < inp.vFractile) lo = m; else hi = m; }
    const zk = 0.5 * (lo + hi);
    const vCap = r.v * Math.sqrt(inp.Ecap / EN);
    const pEx = exceedanceLognormal({ vk: r.v, zk, cov: inp.vCov, vCap });
    const vs = linspace(0.25 * r.v, 2.5 * r.v, 30);
    const eN = vs.map((v) => round(EN * (v / r.v) ** 2)), eA = eN.map((e) => round(e * inp.Fab));
    const xs = linspace(0.05, 0.5, 19), ces = xs.map((x) => { try { return round(evaluate(inp, { xFrac: x }).Ce); } catch { return null; } });
    const xOk = xs.filter((_, i) => ces[i] != null), ceOk = ces.filter((c) => c != null);
    // Envelope: allowable velocity against displacement for geometrically similar vessels (Ce, Cm unchanged).
    const ms = linspace(0.25, 2, 15).map((f) => f * inp.disp);
    const vAllow = ms.map((m) => round(Math.sqrt(2 * inp.Ecap * 1e3 / (inp.Fab * m * 1e3 * r.Ce * r.Cm * inp.Cs * inp.Cc))));
    const dispMax = inp.disp * inp.Ecap / EA;
    const ok = EA <= inp.Ecap;
    const dE10 = EA * (1 - 0.81);
    const checks = [check('Abnormal berthing energy against capacity', EA, inp.Ecap, 'kJ', 'PIANC WG33 (2002) §4.2.8; BS 6349-4:2014 §5.4',
      'Guards against the fender bottoming out and the surplus energy going into the hull and the piles. If exceeded the reaction rises steeply and the structure or hull plating is damaged.')];
    const advice = [
      adv('engineering', ok ? 'ok' : 'fail', `Normal berthing energy ${fmt(EN)} kJ, abnormal ${fmt(EA)} kJ (v = ${fmt(r.v)} m/s, Ce = ${fmt(r.Ce)}, Cm = ${fmt(r.Cm)}); capacity ${fmt(inp.Ecap)} kJ.`,
        ok ? `Reserve of ${fmt(inp.Ecap - EA)} kJ. Confirm the velocity with port berthing records and repeat for the ballast and the smallest design vessel, which can govern hull pressure.` : `Provide at least ${fmt(EA)} kJ at minimum tolerance (after velocity, temperature and angle corrections), or restrict the vessel to ${fmt(dispMax)} t displacement, or the approach to ${fmt(vCap / Math.sqrt(inp.Fab))} m/s.`,
        'PIANC WG33 (2002) §4.2; PIANC WG211 (2024); BS 6349-4:2014',
        `Energy scales with velocity squared: a 10 % lower approach velocity removes ${fmt(dE10)} kJ of demand, which carries directly into fender rubber volume and pile steel.`),
      adv('economic', 'info', `The demand of ${fmt(EA)} kJ sets the fender size and the lateral design load of the dolphin; velocity is the dominant driver (±10 % in velocity is ±${fmt(EA * 0.21)} kJ).`,
        'Before enlarging fenders or piles, compare the cost of a berthing-aid (laser docking) system and enforced approach limits, which lower both the design velocity and its scatter.', `${ECON_BASIS}; value-of-information on measured berthing velocities`),
    ];
    return {
      kpis: [
        { key: 'EN', label: 'Normal berthing energy', unit: 'kJ', value: EN, explain: 'Energy the fender system absorbs in a normal berthing: ½·M·v² reduced by eccentricity and increased by added mass. The fender must absorb it at minimum manufacturing tolerance without exceeding its rated deflection.' },
        { key: 'EA', label: 'Abnormal berthing energy', unit: 'kJ', value: EA, limit: inp.Ecap, sense: 'max', explain: 'Normal energy times the abnormal-berthing factor, covering mishandling, engine failure or squalls. It must not exceed the energy capacity of fender plus structure.' },
        { key: 'v', label: 'Berthing velocity', unit: 'm/s', value: r.v, explain: 'Velocity of the vessel normal to the berthing line at first contact. Energy grows with its square, so it is the most influential and most uncertain input.' },
        { key: 'Ce', label: 'Eccentricity factor Ce', unit: '-', value: r.Ce, explain: 'Fraction of the kinetic energy not converted into yaw rotation about the contact point. 1.0 for midships contact; about 0.5 for quarter-point berthing.' },
        { key: 'Cm', label: 'Added-mass factor Cm', unit: '-', value: r.Cm, explain: 'Ratio of vessel plus entrained water mass to vessel mass in sway. Typically 1.5–1.8; it rises as the under-keel clearance shrinks.' },
        { key: 'Cb', label: 'Block coefficient', unit: '-', value: r.Cb, explain: 'Displaced volume over length × beam × draft. A check on the vessel data: 0.8–0.85 for tankers and bulk carriers, 0.6–0.7 for container ships.' },
        { key: 'K', label: 'Radius of gyration', unit: 'm', value: r.K, explain: 'Yaw radius of gyration, (0.19·Cb + 0.11)·L. Together with the contact distance it fixes how much energy goes into rotation.' },
        { key: 'R', label: 'Contact point to centre of mass', unit: 'm', value: r.R, explain: 'Distance from the point of contact to the centre of mass. The larger it is, the more the vessel pivots and the less the fender must absorb.' },
        { key: 'gammaDeg', label: 'Velocity-vector angle γ', unit: 'deg', value: r.gamma / DEG, explain: 'Angle between the velocity vector and the line from the contact point to the centre of mass. 90° gives the smallest eccentricity factor.' },
        { key: 'momentum', label: 'Vessel momentum normal to the berth', unit: 'kN·s', value: r.M * r.v / 1e3, explain: 'Linear momentum the dolphin arrests or redirects. Divided by the peak fender reaction it gives the order of the contact duration.' },
        { key: 'Erot', label: 'Energy remaining in rotation', unit: 'kJ', value: r.E0 * r.Cm * (1 - r.Ce) / 1e3, explain: 'Kinetic energy left in the vessel as yaw after first contact (energy balance: ½·M·v²·Cm = absorbed + rotation). It is delivered to the second dolphin when the stern swings in.' },
        { key: 'pExceed', label: 'Probability a berthing exceeds the capacity', unit: '-', value: pEx, explain: 'Chance per berthing that the normal-berthing energy exceeds the energy capacity, for a lognormal velocity. Multiply by berthings per year for the annual rate; values above about 1e-4 per berthing warrant attention.' },
        { key: 'vAllow', label: 'Allowable berthing velocity', unit: 'm/s', value: vCap / Math.sqrt(inp.Fab), explain: 'Largest normal velocity for which the abnormal energy stays within the capacity for this vessel: the operating limit to give to pilots.' },
        { key: 'dispMax', label: 'Largest vessel at the design velocity', unit: 't', value: dispMax, explain: 'Displacement at which the abnormal energy equals the capacity, for a vessel of similar proportions approaching at the design velocity.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Berthing energy against approach velocity', x: { label: 'Velocity normal to berth', unit: 'm/s' }, y: { label: 'Energy', unit: 'kJ' },
          series: [{ name: 'Normal', x: vs.map((v) => round(v)), y: eN }, { name: 'Abnormal', x: vs.map((v) => round(v)), y: eA }, { name: 'Capacity', x: [round(vs[0]), round(vs[29])], y: [inp.Ecap, inp.Ecap], dash: true }],
          explain: 'Quadratic growth of energy with velocity. The crossing of the abnormal curve with the dashed capacity line is the velocity limit of the berth.' },
        { type: 'line', title: 'Eccentricity factor against contact position', x: { label: 'Contact point from bow / LBP', unit: '-' }, y: { label: 'Ce', unit: '-' },
          series: [{ name: 'Ce', x: xOk.map((x) => round(x)), y: ceOk }], explain: 'How the absorbed fraction changes as the contact moves from the bow towards midships. Dolphins placed nearer midships see more energy.' },
        { type: 'line', title: 'Allowable berthing velocity against vessel size', x: { label: 'Displacement', unit: 't' }, y: { label: 'Allowable velocity', unit: 'm/s' },
          series: [{ name: 'Capacity envelope', x: ms.map((m) => round(m)), y: vAllow }, { name: 'Design vessel', x: [inp.disp], y: [round(r.v)] }],
          explain: 'Operating envelope of the berth: combinations of vessel size and velocity below the curve keep the abnormal energy within the capacity. The design point should lie under it.' },
      ],
      tables: [{ title: 'Energy build-up', columns: [{ key: 'item', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
        rows: [['Kinetic energy ½·M·v²', r.E0 / 1e3, 'kJ'], ['× added mass Cm', r.E0 * r.Cm / 1e3, 'kJ'], ['× eccentricity Ce', r.E0 * r.Cm * r.Ce / 1e3, 'kJ'], ['× softness Cs and configuration Cc (normal energy)', EN, 'kJ'], ['× abnormal factor', EA, 'kJ'], ['Under-keel clearance', r.Kc, 'm'], ['Velocity at capacity (normal berthing)', vCap, 'm/s']].map(([item, v, u]) => ({ item, v: round(v), u })),
        explain: 'Each row applies one factor of the kinetic-energy method in turn, for checking against a hand calculation.' }],
      notes: [`Velocity ${inp.velMode === 'table' ? `read from the guidance table (condition ${inp.cond}, ${fmt(inp.dwt)} t deadweight)` : 'as specified'}: ${fmt(r.v)} m/s. The guidance table is indicative.`,
        ok ? 'The abnormal berthing energy is within the stated capacity.' : 'The abnormal berthing energy exceeds the stated capacity.'],
      advice,
      publish: { 'dolphin.berthEnergyNormal': r.EN, 'dolphin.berthEnergyAbnormal': r.EA, 'dolphin.berthVelocity': r.v },
      quantities: {},
    };
  },
  verify: [
    { name: 'Midships parallel contact: Ce = 1, E = ½·M·v²', source: 'Hand calculation: ½ · 1.0e7 kg · (0.2 m/s)² = 200 kJ; R = B/2 and γ = 0 give Ce = 1',
      inputs: { disp: 10000, LBP: 120, beam: 20, draft: 7, depth: 12, velMode: 'user', vUser: 0.2, alpha: 0, xFrac: 0.5, cmMethod: 'user', CmUser: 1, Cs: 1, Cc: 1, Fab: 1 },
      expect: [{ kpi: 'EN', value: 200, tol: 1e-6 }, { kpi: 'Ce', value: 1, tol: 1e-9 }] },
    // Cb = 40e6/(1025·200·30·10) = 0.65041; K = (0.19·0.65041 + 0.11)·200 = 46.715 m; R = √(50² + 15²) = 52.202 m;
    // γ = 90 − 10 − asin(15/52.202) = 63.301°; Ce = (2182.3 + 2725·cos²γ)/4907.3 = 0.5568; Cm = 1 + 20/30 = 1.6667;
    // E = 450 kJ · 0.5568 · 1.6667 = 417.6 kJ.
    { name: 'Quarter-point berthing of a 40 000 t vessel at 10°', source: 'Hand calculation with the PIANC WG33 (2002) §4.2 expressions (written out in the module source)',
      inputs: { disp: 40000, LBP: 200, beam: 30, draft: 10, depth: 20, velMode: 'user', vUser: 0.15, alpha: 10, xFrac: 0.25, cmMethod: 'vasco-costa', Cs: 1, Cc: 1, Fab: 1.5 },
      expect: [{ kpi: 'Ce', value: 0.5568, tol: 0.002 }, { kpi: 'Cm', value: 1.66667, tol: 1e-5 }, { kpi: 'EN', value: 417.6, tol: 0.003 }, { kpi: 'EA', value: 626.4, tol: 0.003 }] },
  ],
};
