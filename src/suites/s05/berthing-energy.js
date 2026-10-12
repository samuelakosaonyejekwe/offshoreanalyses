// Berthing energy by the kinetic-energy method, fender selection from the force–deflection curve, the impact
// force passed into the wharf, the one-degree-of-freedom impact history and the permissible berthing envelope.
import { blockCoefficient, radiusOfGyration, contactGeometry, eccentricityFactor, addedMassFactor, berthingEnergy, VELOCITY_TABLE_DEFAULT, velocityFromTable,
  exceedanceLognormal, GENERIC_FENDERS, fenderCurve, fenderAtEnergy, interp } from './_berthing.js';
import { erfinv } from '../../engine/math/special.js';
import { DEG, need, r, f, check, eng, eco, linspace, maxWhere, CONFIRM, NOT_HIFI, ECON_BASIS, CO2_STEEL } from './_common.js';

/** Effective energy coefficients and energies for a vessel of mass M [kg] at velocity v [m/s]. */
function energyOf(inp, M, v) {
  const Cb = blockCoefficient({ M, L: inp.LBP, B: inp.beam, D: inp.draft, rho: inp.rhoW });
  const K = radiusOfGyration({ Cb, L: inp.LBP });
  const g = contactGeometry({ L: inp.LBP, B: inp.beam, x: inp.contact * inp.LBP, alpha: inp.alpha * DEG });
  const Ce = inp.CeUser > 0 ? inp.CeUser : eccentricityFactor({ K, R: g.R, gamma: g.gamma });
  const Cm = addedMassFactor({ method: inp.cmMethod, D: inp.draft, B: inp.beam, Cb, Kc: inp.ukc, user: inp.CmUser });
  const e = berthingEnergy({ M, v, Ce, Cm, Cs: inp.Cs, Cc: inp.Cc, Fab: inp.Fab });
  return { Cb, K, R: g.R, gamma: g.gamma, Ce, Cm, ...e, mEff: M * Ce * Cm * inp.Cs * inp.Cc };
}

export default {
  id: 's05-berthing-energy', suite: 5, title: 'Berthing energy and fender selection',
  summary: 'Normal and abnormal berthing energy of the design vessel, fender deflection and reaction from the performance curve, the force passed into the wharf and the largest permissible approach velocity and vessel.',
  covers: [2], fidelity: 'analytical',
  standards: ['PIANC WG33 Guidelines for the Design of Fender Systems (2002) §4', 'PIANC WG211 Guidelines for the Design, Manufacturing and Testing of Fender Systems (2024)', 'BS 6349-4:2014 §5', 'ROM 2.0-11'],
  models: ['Vessel kinetic energy', 'Effective berthing mass (added mass: Vasco Costa, Ueda, PIANC)', 'Eccentricity factor from rigid-body impact', 'Nonlinear fender force–deflection curve', 'Energy absorption by integration of the curve', 'One-degree-of-freedom impact dynamics', 'Lognormal berthing-velocity exceedance'],
  limits: `Kinetic-energy method for a single side-on contact: the vessel is a rigid body, the fender arrests the contact point and the wharf is rigid compared with the fender. Berthing velocities and the generic fender shapes are indicative table values, not PIANC or manufacturer data: enter the certified curve and the velocities of the edition in force. ${NOT_HIFI}`,
  publishes: ['wharf.berthingEnergy', 'wharf.fenderReaction', 'wharf.fenderReactionBent'],
  searchTerms: 'berthing energy fender design PIANC WG211 eccentricity added mass',
  inputs: [
    { key: 'dwt', label: 'Deadweight of the design vessel', unit: 't', type: 'number', default: 60000, min: 500, max: 550000, group: 'Vessel', role: 'operation', help: 'Used only to read the berthing-velocity table.' },
    { key: 'disp', label: 'Displacement at berthing', unit: 't', type: 'number', default: 80000, min: 500, max: 700000, group: 'Vessel', role: 'ic', bind: 'vessel.disp', uncertain: { dist: 'normal', cov: 0.05 }, help: 'Mass of the vessel in the berthing condition (laden for an import berth).' },
    { key: 'LBP', label: 'Length between perpendiculars', unit: 'm', type: 'number', default: 240, min: 20, max: 420, group: 'Vessel', role: 'geometry' },
    { key: 'beam', label: 'Beam', unit: 'm', type: 'number', default: 37, min: 4, max: 70, group: 'Vessel', role: 'geometry' },
    { key: 'draft', label: 'Draft at berthing', unit: 'm', type: 'number', default: 12.5, min: 1, max: 25, group: 'Vessel', role: 'ic' },
    { key: 'ukc', label: 'Under-keel clearance', unit: 'm', type: 'number', default: 1.5, min: 0.1, max: 30, group: 'Vessel', role: 'bc', help: 'Used by the PIANC added-mass option.' },
    { key: 'rhoW', label: 'Water density', unit: 'kg/m³', type: 'number', default: 1025, min: 995, max: 1035, group: 'Vessel', role: 'environment', bind: 'site.rhoW' },
    { key: 'vMode', label: 'Approach velocity from', type: 'select', default: 'table', group: 'Approach', role: 'ic', options: [{ value: 'table', label: 'Velocity table by navigation condition' }, { value: 'user', label: 'Value entered below' }] },
    { key: 'navCond', label: 'Navigation condition', type: 'select', default: 'b', group: 'Approach', role: 'ic', showIf: { key: 'vMode', in: ['table'] },
      options: [{ value: 'a', label: 'a  Easy berthing, sheltered' }, { value: 'b', label: 'b  Difficult berthing, sheltered' }, { value: 'c', label: 'c  Easy berthing, exposed' }, { value: 'd', label: 'd  Good berthing, exposed' }, { value: 'e', label: 'e  Difficult berthing, exposed' }] },
    { key: 'vUser', label: 'Approach velocity normal to the berth', unit: 'm/s', type: 'number', default: 0.12, min: 0.01, max: 1.5, group: 'Approach', role: 'ic', showIf: { key: 'vMode', in: ['user'] }, uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Characteristic value at first contact, from measured berthings or pilot data where available.' },
    { key: 'vTable', label: 'Berthing-velocity table (indicative)', type: 'table', group: 'Approach', role: 'ic',
      columns: [{ key: 'dwt', label: 'Deadweight', unit: 't', type: 'number' }, { key: 'a', label: 'a', unit: 'm/s', type: 'number' }, { key: 'b', label: 'b', unit: 'm/s', type: 'number' }, { key: 'c', label: 'c', unit: 'm/s', type: 'number' }, { key: 'd', label: 'd', unit: 'm/s', type: 'number' }, { key: 'e', label: 'e', unit: 'm/s', type: 'number' }],
      default: VELOCITY_TABLE_DEFAULT, help: 'Indicative Brolsma-type values as tabulated in fender manuals after PIANC WG33; WG211 (2024) revises the basis using measured berthings. Replace with the values of the edition and port in force.' },
    { key: 'alpha', label: 'Berthing angle', unit: 'deg', type: 'number', default: 6, min: 0, max: 20, group: 'Approach', role: 'ic', uncertain: { dist: 'normal', cov: 0.25 } },
    { key: 'contact', label: 'Contact point from the bow (fraction of length)', unit: '-', type: 'number', default: 0.25, min: 0.05, max: 0.5, group: 'Approach', role: 'ic', help: '0.25 quarter-point (continuous quay), 0.33 third-point, 0.5 midships.' },
    { key: 'cmMethod', label: 'Added-mass factor Cm', type: 'select', default: 'vasco-costa', group: 'Coefficients', role: 'criteria',
      options: [{ value: 'vasco-costa', label: 'Vasco Costa: 1 + 2D/B' }, { value: 'ueda', label: 'Ueda: 1 + πD/(2 Cb B)' }, { value: 'pianc', label: 'PIANC WG33 by under-keel clearance' }, { value: 'user', label: 'User value' }] },
    { key: 'CmUser', label: 'Cm (user)', unit: '-', type: 'number', default: 1.7, min: 1, max: 3, group: 'Coefficients', role: 'criteria', showIf: { key: 'cmMethod', in: ['user'] } },
    { key: 'CeUser', label: 'Eccentricity factor Ce (0 = from geometry)', unit: '-', type: 'number', default: 0, min: 0, max: 1, group: 'Coefficients', role: 'criteria', help: 'Computed as (K² + R² cos²γ)/(K² + R²) when zero.' },
    { key: 'Cs', label: 'Softness factor Cs', unit: '-', type: 'number', default: 1, min: 0.8, max: 1, group: 'Coefficients', role: 'criteria', help: `1.0 for soft fenders (deflection > 150 mm), 0.9 for hard fenders. ${CONFIRM}` },
    { key: 'Cc', label: 'Berth configuration factor Cc', unit: '-', type: 'number', default: 1, min: 0.7, max: 1, group: 'Coefficients', role: 'criteria', help: `1.0 for an open piled wharf, 0.8–0.9 for a solid quay wall with small berthing angle. ${CONFIRM}` },
    { key: 'Fab', label: 'Abnormal berthing factor', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Coefficients', role: 'criteria', help: `PIANC WG33 Table 4.2.5: 1.25–1.75 tankers and bulk carriers, 1.5–2.0 container ships, 1.75 general cargo, 2.0 ferries and tugs. ${CONFIRM}` },
    { key: 'fType', label: 'Fender type (generic curve shape)', type: 'select', default: 'cone', group: 'Fender', role: 'bc', options: Object.entries(GENERIC_FENDERS).map(([value, g]) => ({ value, label: g.label })) },
    { key: 'fH', label: 'Fender height / diameter', unit: 'm', type: 'number', default: 1.3, min: 0.2, max: 4.5, group: 'Fender', role: 'geometry', design: true },
    { key: 'fRR', label: 'Rated reaction of one unit', unit: 'kN', type: 'number', default: 1700, min: 20, max: 12000, group: 'Fender', role: 'bc', design: true },
    { key: 'fCurve', label: 'Certified fender curve (replaces the generic shape)', type: 'table', group: 'Fender', role: 'bc', columns: [{ key: 'defl', label: 'Deflection', unit: '% of height', type: 'number' }, { key: 'reaction', label: 'Reaction', unit: 'kN', type: 'number' }], default: [],
      help: 'Rated performance data of one unit from the certified test (PIANC protocol). Leave empty to use the generic shape scaled to the rated reaction; the last row is the rated point.' },
    { key: 'nF', label: 'Fender units sharing the impact', unit: '-', type: 'int', default: 1, min: 1, max: 6, group: 'Fender', role: 'bc', help: 'Units compressed together by the hull at the berthing angle and fender spacing.' },
    { key: 'fE', label: 'Correction factor on energy (velocity × temperature × angle × tolerance)', unit: '-', type: 'number', default: 0.9, min: 0.4, max: 1.3, group: 'Fender', role: 'criteria', help: 'Product of the manufacturer correction factors applied to the catalogue energy, low-end tolerance.' },
    { key: 'fR', label: 'Correction factor on reaction (high-end tolerance)', unit: '-', type: 'number', default: 1.1, min: 0.6, max: 1.6, group: 'Fender', role: 'criteria' },
    { key: 'panelA', label: 'Fender panel contact area', unit: 'm²', type: 'number', default: 12, min: 0.2, max: 80, group: 'Fender', role: 'geometry', design: true },
    { key: 'pHull', label: 'Permissible hull pressure', unit: 'kPa', type: 'number', default: 200, min: 50, max: 1000, group: 'Fender', role: 'criteria', help: `Indicative: 150–250 kPa tankers and bulk carriers, up to 400 kPa container ships (PIANC WG33 Table 4.4.1). ${CONFIRM}` },
    { key: 'mu', label: 'Friction coefficient at the panel face', unit: '-', type: 'number', default: 0.2, min: 0, max: 0.6, group: 'Fender', role: 'bc', help: '0.2 for UHMW-PE facing on steel.' },
    { key: 'Rallow', label: 'Design horizontal resistance of the wharf at one fender', unit: 'kN', type: 'number', default: 3000, min: 50, max: 50000, group: 'Structure', role: 'criteria', help: 'Factored resistance of the fender support and deck diaphragm to a single fender reaction, from the global strength module.' },
    { key: 'nShare', label: 'Bents sharing one fender reaction', unit: '-', type: 'number', default: 3, min: 1, max: 12, group: 'Structure', role: 'bc', help: 'Load spread through the deck diaphragm; used for the reaction published to the bent analysis.' },
    { key: 'covV', label: 'Coefficient of variation of the berthing velocity', unit: '-', type: 'number', default: 0.4, min: 0.1, max: 1, group: 'Traffic', role: 'operation', help: 'Scatter of measured berthings about the mean; 0.3–0.5 in published records.' },
    { key: 'pChar', label: 'Non-exceedance probability of the characteristic velocity', unit: '-', type: 'number', default: 0.98, min: 0.5, max: 0.999, group: 'Traffic', role: 'operation' },
    { key: 'calls', label: 'Berthings per year', unit: '1/yr', type: 'number', default: 150, min: 1, max: 5000, group: 'Traffic', role: 'operation' },
    { key: 'nStep', label: 'Time steps of the impact history', unit: '-', type: 'int', default: 400, min: 50, max: 6400, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nStep', label: 'Time steps', levels: [100, 200, 400, 800], qoi: 'dDyn' },
  compute(inp) {
    const M = inp.disp * 1e3;
    need(inp.disp >= 0.3 * inp.rhoW / 1e3 * inp.LBP * inp.beam * inp.draft, `A displacement of ${f(inp.disp)} t is too small for a hull of ${f(inp.LBP)} × ${f(inp.beam)} × ${f(inp.draft)} m (block coefficient below 0.3). Check the vessel data.`);
    const v = inp.vMode === 'table' ? velocityFromTable((inp.vTable || []).map((q) => ({ dwt: Number(q.dwt), a: Number(q.a), b: Number(q.b), c: Number(q.c), d: Number(q.d), e: Number(q.e) })), inp.dwt, inp.navCond) : inp.vUser;
    const e = energyOf(inp, M, v);
    need(e.Cb <= 1.0, `The displacement of ${f(inp.disp)} t exceeds the box volume of the hull (block coefficient ${f(e.Cb)}). Check the vessel data.`);
    const rows = (inp.fCurve || []).map((q) => ({ defl: Number(q.defl), reaction: Number(q.reaction) })).filter((q) => Number.isFinite(q.defl) && Number.isFinite(q.reaction));
    const curve = fenderCurve({ type: inp.fType, H: inp.fH, RR: inp.fRR * 1e3, table: rows });
    const nF = Math.round(inp.nF), fac = { fE: inp.fE, fR: inp.fR };
    const sN = fenderAtEnergy(curve, e.EN / nF, fac), sA = fenderAtEnergy(curve, e.EA / nF, fac);
    const Ecap = nF * curve.ERated * inp.fE, EcapMax = nF * sA.energyCapacity;
    const Rn = nF * sN.peakReaction, Ra = nF * sA.peakReaction;
    const pHull = Rn / nF / inp.panelA;
    // permissible conditions: the abnormal energy must stay within the rated energy, and the peak reaction within the wharf resistance
    const okAt = (Mx, vx) => { const q = energyOf(inp, Mx, vx); const s = fenderAtEnergy(curve, q.EA / nF, fac); return q.EA <= Ecap && nF * s.peakReaction <= inp.Rallow * 1e3 && !s.exceeded; };
    const vPerm = maxWhere((x) => okAt(M, x), 0.005, 1.5);
    const dispPerm = maxWhere((x) => okAt(x * 1e3, v), 0.2 * inp.disp, 5 * inp.disp);
    const dispGrid = linspace(0.3 * inp.disp, 1.6 * inp.disp, 14), vEnv = dispGrid.map((d) => maxWhere((x) => okAt(d * 1e3, x), 0.005, 1.5));
    // probability that a berthing exceeds the fender's rated energy (energy ∝ v²)
    const zk = Math.sqrt(2) * erfinv(2 * inp.pChar - 1);
    const vCap = v * Math.sqrt(Ecap / e.EN);
    const pEx = exceedanceLognormal({ vk: v, zk, cov: inp.covV, vCap });
    const pAnnual = 1 - (1 - pEx) ** inp.calls;
    // impact history: m_eff ẍ = −nF·R(x), velocity-Verlet until the vessel stops (rate effects are in the correction factors)
    const nS = Math.round(inp.nStep), force = (x) => nF * inp.fE * interp(curve.d, curve.R, x); // energy-consistent corrected curve
    const dt = (2.5 * sN.deflection / v) / nS;
    let x = 0, vel = v, a = -force(0) / e.mEff, tt = 0, xPrev = 0, vPrev = v, th = [0], xh = [0], fh = [0];
    for (let i = 0; i < 4 * nS && vel > 0; i++) {
      xPrev = x; vPrev = vel;
      x += vel * dt + 0.5 * a * dt * dt;
      const an = -force(Math.min(x, curve.dMax)) / e.mEff;
      vel += 0.5 * (a + an) * dt; a = an; tt += dt;
      th.push(tt); xh.push(x); fh.push(force(Math.min(x, curve.dMax)));
    }
    // deflection where the velocity passes zero (linear interpolation within the last step)
    const dDyn = vel < 0 ? xPrev + (x - xPrev) * vPrev / (vPrev - vel) : x;
    const tStop = vel < 0 ? tt - dt + dt * vPrev / (vPrev - vel) : tt;
    const okE = e.EA <= Ecap, okR = Ra <= inp.Rallow * 1e3, okP = pHull <= inp.pHull * 1e3;
    const needRR = inp.fRR * (e.EA / Ecap);
    const dPct = curve.d.map((q) => r(q / inp.fH * 100, 4));
    const stride = Math.max(1, Math.floor(th.length / 250));
    const pick = (aa) => aa.filter((_, i) => i % stride === 0 || i === aa.length - 1);
    return {
      kpis: [
        { key: 'v', label: 'Approach velocity normal to the berth', value: r(v), unit: 'm/s', explain: 'Characteristic velocity of the vessel at first contact. Energy grows with its square, so it is the most influential and least certain input.' },
        { key: 'Cb', label: 'Block coefficient', value: r(e.Cb), unit: '-', explain: 'Displaced volume over the box L × B × D. 0.6–0.7 for container ships, 0.8–0.85 for tankers and bulk carriers; a value outside 0.5–0.9 suggests inconsistent vessel data.' },
        { key: 'Ce', label: 'Eccentricity factor', value: r(e.Ce), unit: '-', explain: 'Share of the kinetic energy that the fender must absorb; the rest stays in the vessel as rotation about the contact point. 0.4–0.7 for quarter-point contact, 1.0 for midships.' },
        { key: 'Cm', label: 'Added-mass factor', value: r(e.Cm), unit: '-', explain: 'Vessel plus entrained water over vessel mass. Typically 1.5–1.9; it rises as the under-keel clearance shrinks.' },
        { key: 'EN', label: 'Normal berthing energy', value: r(e.EN / 1e3), unit: 'kJ', explain: 'Energy of a normal berthing, ½ M v² Ce Cm Cs Cc. The fender absorbs it within its normal working range.' },
        { key: 'EA', label: 'Abnormal berthing energy', value: r(e.EA / 1e3), unit: 'kJ', limit: r(Ecap / 1e3), sense: 'max', explain: 'Normal energy times the abnormal-berthing factor: the energy the fender system must absorb without exceeding its rated deflection. Above the limit the fender bottoms out and the reaction rises steeply.' },
        { key: 'Efender', label: 'Rated energy capacity of the fender system', value: r(Ecap / 1e3), unit: 'kJ', explain: 'Area under the corrected force–deflection curve up to the rated deflection, for all units sharing the impact.' },
        { key: 'deflN', label: 'Fender deflection, normal berthing', value: r(sN.deflection / inp.fH * 100), unit: '%', explain: 'Compression as a percentage of fender height in a normal berthing.' },
        { key: 'deflA', label: 'Fender deflection, abnormal berthing', value: r(sA.deflection / inp.fH * 100), unit: '%', limit: r(curve.dRated / inp.fH * 100), sense: 'max', explain: 'Compression in the abnormal berthing against the rated deflection of the unit. Beyond the rated deflection the rubber is over-compressed and may be damaged.' },
        { key: 'RN', label: 'Fender reaction, normal berthing', value: r(Rn / 1e3), unit: 'kN', explain: 'Largest reaction passed on the way to the normal-berthing deflection, with the high-end tolerance factor. This is the characteristic berthing force on the wharf.' },
        { key: 'RA', label: 'Fender reaction, abnormal berthing', value: r(Ra / 1e3), unit: 'kN', limit: r(inp.Rallow), sense: 'max', explain: 'Largest reaction in the abnormal berthing: the impact force transmitted into the supporting structure, compared with its design resistance.' },
        { key: 'pHull', label: 'Hull pressure', value: r(pHull / 1e3), unit: 'kPa', limit: inp.pHull, sense: 'max', explain: 'Normal-berthing reaction over the panel contact area. Above the permissible value the hull plating between frames may dent.' },
        { key: 'shear', label: 'Friction force along the berth', value: r(inp.mu * Rn / 1e3), unit: 'kN', explain: 'Friction between hull and panel, μ × reaction, carried by the fender chains and the deck in the longitudinal direction.' },
        { key: 'vPerm', label: 'Maximum permissible approach velocity', value: r(vPerm), unit: 'm/s', explain: 'Largest characteristic velocity of this vessel for which the abnormal energy stays within the rated fender energy and the reaction within the wharf resistance. A berthing-aid alarm can be set from it.' },
        { key: 'dispPerm', label: 'Maximum permissible displacement at the design velocity', value: r(dispPerm), unit: 't', explain: 'Largest vessel (same dimensions) that can berth at the design velocity within the fender and wharf limits. Above five times the design vessel it is reported at that cap.' },
        { key: 'pExceed', label: 'Annual probability of exceeding the rated fender energy', value: r(pAnnual, 3), unit: '-', explain: 'From a lognormal berthing-velocity distribution through the characteristic value, over the yearly number of berthings. It measures how often the abnormal reserve is called on.' },
        { key: 'dDyn', label: 'Peak deflection from the impact history', value: r(dDyn * 1e3, 7), unit: 'mm', explain: 'Peak fender compression from time integration of the vessel decelerating against the fender. It should match the energy-method deflection: the difference is the time-step error.' },
        { key: 'tStop', label: 'Time to arrest the vessel', value: r(tStop), unit: 's', explain: 'Duration of the first compression stroke. Several seconds means the wharf responds quasi-statically; it is far from the structural natural period.' },
      ],
      checks: [
        check('Fender energy capacity (abnormal berthing)', e.EA / 1e3, Ecap / 1e3, 'kJ', 'PIANC WG33 (2002) §4.2.8; BS 6349-4:2014 §5.6', 'Guards against bottoming the fender. If the energy exceeds the rated capacity the reaction climbs steeply and the hull and wharf take the excess.'),
        check('Impact force on the wharf', Ra / 1e3, inp.Rallow, 'kN', 'BS 6349-2:2019 §6 (berthing actions); PIANC WG33 §4.3', 'Guards against overloading the fender support, deck edge and piles. Exceeding it means local cracking at the cope or overstressed front piles.'),
        check('Hull pressure', pHull / 1e3, inp.pHull, 'kPa', 'PIANC WG33 (2002) §4.4, Table 4.4.1', 'Guards against denting the ship side shell. A larger panel spreads the reaction.'),
      ],
      plots: [
        { type: 'line', title: 'Fender reaction and energy against deflection', x: { label: 'Deflection', unit: '% of height' }, y: { label: 'Reaction (kN) / energy (kJ)', unit: '-' },
          series: [{ name: 'Reaction, all units (kN)', x: dPct, y: curve.R.map((q) => r(nF * q * inp.fR / 1e3)) }, { name: 'Energy, all units (kJ)', x: dPct, y: curve.E.map((q) => r(nF * q * inp.fE / 1e3)) },
            { name: 'Normal berthing', x: [r(sN.deflection / inp.fH * 100), r(sN.deflection / inp.fH * 100)], y: [0, r(e.EN / 1e3)], dash: true }, { name: 'Abnormal berthing', x: [r(sA.deflection / inp.fH * 100), r(sA.deflection / inp.fH * 100)], y: [0, r(e.EA / 1e3)], dash: true }],
          explain: 'Corrected performance curve of the fender system with the two working points. A buckling fender reaches its plateau reaction early; look for the abnormal point staying left of the rated deflection.' },
        { type: 'line', title: 'Permissible approach velocity against displacement', x: { label: 'Displacement', unit: 't' }, y: { label: 'Velocity', unit: 'm/s' },
          series: [{ name: 'Permissible', x: dispGrid.map((d) => r(d)), y: vEnv.map((q) => r(q)) }, { name: 'Design point', x: [inp.disp], y: [r(v)] }],
          explain: 'Berthing envelope of the installed fender and wharf: any vessel and velocity below the curve is acceptable. The design point should sit below it with margin.' },
        { type: 'line', title: 'Impact history (normal berthing)', x: { label: 'Time', unit: 's' }, y: { label: 'Deflection (mm) / reaction (kN)', unit: '-' },
          series: [{ name: 'Deflection (mm)', x: pick(th).map((q) => r(q)), y: pick(xh).map((q) => r(q * 1e3)) }, { name: 'Reaction (kN)', x: pick(th).map((q) => r(q)), y: pick(fh).map((q) => r(q / 1e3)) }],
          explain: 'Compression stroke of the fender as the effective vessel mass decelerates. The reaction history is the load-time function to use if a dynamic check of the wharf is wanted.' },
      ],
      tables: [{ title: 'Berthing energy build-up', columns: [{ key: 'item', label: 'Quantity' }, { key: 'value', label: 'Value' }, { key: 'unit', label: 'Unit' }],
        rows: [{ item: 'Displacement M', value: r(inp.disp), unit: 't' }, { item: 'Velocity v', value: r(v), unit: 'm/s' }, { item: '½ M v²', value: r(0.5 * M * v * v / 1e3), unit: 'kJ' }, { item: 'Radius of gyration K', value: r(e.K), unit: 'm' }, { item: 'Contact distance R', value: r(e.R), unit: 'm' }, { item: 'Velocity angle γ', value: r(e.gamma / DEG), unit: 'deg' }, { item: 'Ce', value: r(e.Ce), unit: '-' }, { item: 'Cm', value: r(e.Cm), unit: '-' }, { item: 'Cs × Cc', value: r(inp.Cs * inp.Cc), unit: '-' }, { item: 'Normal energy', value: r(e.EN / 1e3), unit: 'kJ' }, { item: 'Abnormal factor', value: inp.Fab, unit: '-' }, { item: 'Abnormal energy', value: r(e.EA / 1e3), unit: 'kJ' }],
        explain: 'Each factor of the kinetic-energy method in turn, so the calculation can be checked by hand and compared with the fender supplier\'s sheet.' }],
      notes: [`${curve.generic ? 'Generic curve shape scaled to the rated reaction (not manufacturer data)' : 'Certified curve from the table'}; ${nF} unit(s) share the impact. Energy method deflection ${f(sN.deflection * 1e3)} mm against ${f(dDyn * 1e3)} mm from the impact history.`,
        `${okE ? 'The fender absorbs the abnormal energy within its rated deflection' : 'The abnormal energy exceeds the rated fender energy'}; ${Ra > Rn * 1.001 ? 'the abnormal reaction governs the wharf' : 'the reaction is on the plateau, so the abnormal case adds no force'}.`],
      advice: [
        okE && okR && okP
          ? eng('ok', `Abnormal energy ${f(e.EA / 1e3)} kJ against ${f(Ecap / 1e3)} kJ rated (${f(e.EA / Ecap * 100)} %); reaction ${f(Ra / 1e3)} kN against ${f(inp.Rallow)} kN; hull pressure ${f(pHull / 1e3)} kPa.`,
            `The fender system is adequate. Confirm the velocity (${f(v)} m/s) with pilot or berthing-aid records and replace the generic curve with the certified one before ordering; the permissible velocity is ${f(vPerm)} m/s.`, 'PIANC WG33 (2002) §4; PIANC WG211 (2024); BS 6349-4:2014',
            `Pile moments and fender-support forces scale directly with the fender reaction: every 10 % of reaction avoided by choosing a taller, softer unit for the same energy is 10 % less bending demand on the front piles, with no extra material.`)
          : eng('fail', `${!okE ? `Abnormal energy ${f(e.EA / 1e3)} kJ exceeds the rated ${f(Ecap / 1e3)} kJ. ` : ''}${!okR ? `Reaction ${f(Ra / 1e3)} kN exceeds the wharf resistance of ${f(inp.Rallow)} kN. ` : ''}${!okP ? `Hull pressure ${f(pHull / 1e3)} kPa exceeds ${f(inp.pHull)} kPa.` : ''}`,
            `${!okE ? `Select a unit of about ${f(needRR)} kN rated reaction at the same height, or a taller unit with the same reaction (energy scales with height). ` : ''}${!okR ? `Use a taller, softer fender or limit the approach velocity to ${f(vPerm)} m/s. ` : ''}${!okP ? `Increase the panel to at least ${f(Rn / nF / (inp.pHull * 1e3))} m².` : ''}`, 'PIANC WG33 (2002) §4; BS 6349-4:2014 §5',
            `A taller fender of the same reaction absorbs more energy without adding load to the wharf; raising the rated reaction to ${f(needRR)} kN instead raises the pile and fender-support demand by ${f(Math.max(needRR / inp.fRR - 1, 0) * 100)} %, which is paid for in steel at about ${CO2_STEEL} t CO₂e per tonne (screening factor).`),
        eco(pAnnual > 0.05 ? 'warn' : 'info', `With ${f(inp.calls)} berthings a year the chance of at least one berthing above the rated fender energy is ${f(pAnnual * 100)} % per year; the permissible envelope is ${f(vPerm)} m/s for this vessel and ${f(dispPerm)} t at the design velocity.`,
          'Weigh a larger fender against operational control: a berthing-aid system with a velocity alarm at the permissible value costs far less than strengthening the wharf, and fender replacement after an overload means berth downtime of days to weeks.', `${ECON_BASIS}; PIANC WG211 (2024) reliability-based partial factors`),
      ],
      publish: { 'wharf.berthingEnergy': e.EA, 'wharf.fenderReaction': Rn, 'wharf.fenderReactionBent': Rn / inp.nShare },
    };
  },
  verify: [
    // Cb = 1e8/(1025·250·40·14) = 0.6969; K = (0.19·0.6969 + 0.11)·250 = 60.60 m; R = √(62.5² + 20²) = 65.62 m;
    // γ = 90 − 6 − asin(20/65.62) = 66.26°; Ce = (60.60² + 65.62²·cos²66.26°)/(60.60² + 65.62²) = (3672.5 + 698.3)/7978.7 = 0.5478;
    // Cm = 1 + 2·14/40 = 1.7; E = ½·1e8·0.15²·0.5478·1.7 = 1.0477e6 J.
    { name: 'Kinetic-energy method, quarter-point berthing of a 100 000 t vessel', source: 'Hand calculation with the PIANC WG33 (2002) §4.2 expressions (Vasco Costa added mass)',
      inputs: { disp: 100000, LBP: 250, beam: 40, draft: 14, vMode: 'user', vUser: 0.15, alpha: 6, contact: 0.25, cmMethod: 'vasco-costa', Cs: 1, Cc: 1, Fab: 1, fRR: 3000, fH: 2 },
      expect: [{ kpi: 'Ce', value: 0.5478, tol: 0.002 }, { kpi: 'Cm', value: 1.7, tol: 1e-6 }, { kpi: 'EN', value: 1047.7, tol: 0.002 }] },
    // Linear spring fender, 1000 kN at 50 % of a 1 m unit: E_rated = ½·1000·0.5 = 250 kJ. Midships contact, Ce = 1, Cm = 1:
    // E = ½·2e7·0.1² = 100 kJ → δ = 0.5·√(100/250) = 0.3162 m = 31.62 %, R = 632.5 kN.
    { name: 'Linear fender: energy ½Rδ and deflection δ = δr √(E/Er)', source: 'Closed form for a linear spring',
      inputs: { disp: 20000, LBP: 150, beam: 25, draft: 8, vMode: 'user', vUser: 0.1, contact: 0.5, alpha: 0, CeUser: 1, cmMethod: 'user', CmUser: 1, Cs: 1, Cc: 1, Fab: 1, fE: 1, fR: 1, fH: 1, fCurve: [{ defl: 0, reaction: 0 }, { defl: 50, reaction: 1000 }], nStep: 800 },
      expect: [{ kpi: 'Efender', value: 250, tol: 1e-3 }, { kpi: 'deflN', value: 31.62, tol: 0.002 }, { kpi: 'RN', value: 632.5, tol: 0.002 }, { kpi: 'dDyn', value: 316.2, tol: 0.005 }] },
  ],
};

