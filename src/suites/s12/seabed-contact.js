// Umbilical–seabed contact: touchdown-zone reactions and curvature, touchdown movement with vessel offset, axial
// friction anchoring towards the termination, embedment and lateral on-bottom stability of the grounded length.
import { riserOffsetCases, boundaryLayerTdp } from '../../engine/lines/riser.js';
import { embedment } from '../../engine/pipeline/pipesoil.js';
import { cableOnBottomStability } from '../../engine/xsection/cablemech.js';
import { lineInputs, configInputs, solveConfig, adv, check, worst, sev, fmt, round, run, need, thin, linspace, G, LIMIT_2D, STD_UMB, vesselCo2Text } from './_common.js';

export default {
  id: 's12-seabed-contact', suite: 12, title: 'Umbilical–seabed contact',
  summary: 'Contact reactions, touchdown movement, friction anchoring, embedment and lateral stability of the umbilical on the seabed.',
  covers: [7], fidelity: 'numerical',
  standards: ['DNV-RP-F114 (2017) §4 (pipe–soil interaction)', 'DNV-RP-F109 (2021) §3 (on-bottom stability, absolute method in simplified form)', `${STD_UMB} §7.3`, 'Plunkett (1967) J. Eng. Ind. 89 (touchdown boundary layer)'],
  models: ['Touchdown-zone contact on linear lift-off springs', 'Nonlinear seabed friction (Coulomb, axial)', 'Stiffened-catenary boundary layer', 'Bearing-capacity embedment', "Morison drag and lift on a seabed cylinder"],
  limits: `${LIMIT_2D} Trenching, soil suction and cyclic degradation of the touchdown soil are not modelled; lateral stability is a static force balance with a combined near-bed velocity, not a dynamic stability analysis.`,
  publishes: ['umbilical.tdpTravel', 'umbilical.groundedLength', 'umbilical.embedment'],
  searchTerms: 'umbilical touchdown seabed interaction friction embedment stability',
  inputs: [
    ...lineInputs(), ...configInputs(),
    { key: 'offsetPct', label: 'Vessel offset (near and far)', unit: '% depth', type: 'number', default: 8, min: 0, max: 20, group: 'Vessel offset', role: 'bc', uncertain: { dist: 'normal', cov: 0.15 } },
    { key: 'soil', label: 'Seabed soil', type: 'select', source: 'materials:soil', default: 'clay-soft', group: 'Seabed', role: 'soil', help: 'Strength and unit weight for the embedment estimate (typical values — site investigation data should replace them).' },
    { key: 'kDyn', label: 'Dynamic embedment factor', unit: '-', type: 'number', default: 1.5, min: 1, max: 5, group: 'Seabed', role: 'soil', help: 'Extra penetration from touchdown motion, typically 1–3 on clay (DNV-RP-F114).' },
    { key: 'muLat', label: 'Lateral seabed friction coefficient', unit: '-', type: 'number', default: 0.5, min: 0.05, max: 1.5, group: 'Seabed', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'DNV-RP-F109 uses 0.2 on clay and 0.6 on sand for pipes; embedment adds passive resistance not counted here.' },
    { key: 'Ubed', label: 'Near-bed steady current (normal to the line)', unit: 'm/s', type: 'number', default: 0.4, min: 0, max: 3, group: 'Near-bed flow', role: 'bc', bind: 'site.UcBed100', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'Uwave', label: 'Wave-induced near-bed velocity amplitude', unit: 'm/s', type: 'number', default: 0.1, min: 0, max: 3, group: 'Near-bed flow', role: 'bc', help: 'Significant orbital velocity at the seabed in the design storm (negligible in deep water).' },
    { key: 'CdBed', label: 'Drag coefficient on the seabed', unit: '-', type: 'number', default: 1.2, min: 0.3, max: 3, group: 'Near-bed flow', role: 'load' },
    { key: 'ClBed', label: 'Lift coefficient on the seabed', unit: '-', type: 'number', default: 0.9, min: 0, max: 3, group: 'Near-bed flow', role: 'load' },
    { key: 'sfStab', label: 'Stability safety factor', unit: '-', type: 'number', default: 1.1, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'termRating', label: 'Termination axial pull rating', unit: 'kN', type: 'number', default: 60, min: 1, max: 5000, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'minGrounded', label: 'Minimum grounded length in the far position', unit: 'm', type: 'number', default: 30, min: 0, max: 2000, group: 'Acceptance criteria', role: 'criteria', help: 'Keeps the touchdown point away from the termination so that it is not lifted or bent.' },
    { key: 'nElem', label: 'Rod elements', type: 'int', default: 120, min: 40, max: 1200, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Rod elements', levels: [60, 120, 240, 480], qoi: 'seabedForce' },
  compute(inp, ctx) {
    const { res, line } = solveConfig(inp, { nElem: inp.nElem });
    need(res.tdp, 'The umbilical does not touch the seabed before the termination: there is no touchdown zone to assess. Lengthen the line or move the termination.');
    const n = res.s.length - 1, L = inp.length;
    // touchdown located by interpolating the elevation to the seabed between the last free and first grounded node
    const td = (r) => {
      if (!r || !r.tdp) return null;
      const k = r.tdp.node, z0 = r.z[k - 1], z1 = r.z[k], f = z0 > 0 && z1 < z0 ? Math.min(1, z0 / (z0 - z1)) : 0.5;
      return { s: r.s[k - 1] + f * (r.s[k] - r.s[k - 1]), x: r.x[k - 1] + f * (r.x[k] - r.x[k - 1]), node: k };
    };
    const off = inp.offsetPct / 100 * inp.depth;
    const cases = off > 0 ? run(() => riserOffsetCases(res, [-off, off]), 'Offset cases did not solve') : [];
    const near = cases[0]?.result, far = cases[1]?.result;
    // seabed force = ∫ reaction ds (trapezoid); peak reaction and its ratio to the submerged weight
    let Fbed = 0, rMax = 0, sPeak = 0;
    for (let k = 1; k <= n; k++) { const ds = res.s[k] - res.s[k - 1]; Fbed += 0.5 * (res.reaction[k] + res.reaction[k - 1]) * ds; }
    for (let k = 0; k <= n; k++) if (res.reaction[k] > rMax) { rMax = res.reaction[k]; sPeak = res.s[k]; }
    const tdM = td(res), grounded = L - tdM.s, kLay = Math.max(1, rMax / line.w);
    const Htdp = Math.max(res.Te[res.tdp.node], 1);
    const bl = boundaryLayerTdp({ EI: line.EI, H: Htdp, w: line.w });
    const tdN = td(near), tdF = td(far);
    const groundedFar = tdF ? L - tdF.s : (far ? 0 : grounded);
    const tdpS = [tdN?.s, tdM.s, tdF?.s].filter((v) => v != null), tdpX = [tdN?.x, tdM.x, tdF?.x].filter((v) => v != null);
    const travelS = Math.max(...tdpS) - Math.min(...tdpS), travelX = Math.max(...tdpX) - Math.min(...tdpX);
    // axial friction anchoring: length needed to shed the touchdown tension, and what reaches the termination
    const HtdpFar = far?.tdp ? Math.max(far.Te[far.tdp.node], 0) : Htdp;
    const anchorLen = inp.muAx > 0 ? HtdpFar / (inp.muAx * line.w) : Infinity;
    const fricAvail = inp.muAx * line.w * groundedFar, termResidual = Math.max(0, HtdpFar - fricAvail);
    const termLoad = Math.max(Math.abs((far || res).endB.Fx), Math.abs(res.endB.Fx));
    // embedment from bearing capacity under the concentrated touchdown reaction
    const soil = ctx.material(inp.soil); need(soil, `Soil '${inp.soil}' is not in the materials library.`);
    const clay = soil.type === 'clay';
    let emb = null, embNote = '';
    try {
      emb = embedment({ D: line.D, W: line.w, model: clay ? 'undrained' : 'drained', kLay, kDyn: inp.kDyn, su0: soil.su, suGrad: soil.suGrad ?? 0, gammaSub: soil.gammaSub, phiDeg: soil.phiDeg });
    } catch (e) { embNote = `Embedment could not be evaluated (${e.message}); the soil may be too weak for the touchdown load within two diameters.`; }
    const zEmb = emb ? emb.z : 2 * line.D;
    const Uc = inp.Ubed + inp.Uwave;
    const stab = Uc > 0 ? run(() => cableOnBottomStability({ submergedWeight: line.w, outerDiameter: line.D, velocity: Uc, mu: inp.muLat, Cd: inp.CdBed, Cl: inp.ClBed, rho: line.rhoW, safetyFactor: inp.sfStab, mass: line.m }), 'On-bottom stability') : null;
    const checks = [
      check('Lateral on-bottom stability', stab ? stab.requiredWeight : 0, line.w, 'N/m', 'DNV-RP-F109 (2021) §3.6 (absolute stability, simplified static balance)', 'Demand is the submerged weight needed to resist drag and lift with the safety factor, capacity the actual weight. Above 1.0 the grounded umbilical slides sideways in the design near-bed flow and can be dragged into free spans or other assets.'),
      check('Axial load on the termination', termLoad / 1e3, inp.termRating, 'kN', 'Project termination rating; DNV-RP-F114 (2017) §4.3 (axial resistance)', 'Horizontal tension left at the termination after seabed friction has shed part of the touchdown tension. Exceeding the rating drags the termination unit or overloads its connector.'),
      check('Grounded length in the far position', inp.minGrounded, Math.max(groundedFar, 1e-6), 'm', 'Project criterion (touchdown clear of the termination)', 'Demand is the required grounded length, capacity the length still on the seabed in the far position. Too little and the touchdown point reaches the termination, bending it and loading it vertically.'),
    ];
    const gov = worst(checks);
    // plots
    const i0 = Math.max(0, res.tdp.node - Math.round(n * 0.12)), i1 = n;
    const sl = (a, f = 1) => Array.from(a.slice(i0, i1 + 1), (v) => v * f);
    const [sR, rR, fR, kR] = thin([sl(res.s), sl(res.reaction), sl(res.friction).map(Math.abs), sl(res.curvature).map(Math.abs)], 300);
    const xi = linspace(0, Math.min(12 * bl.lambda, res.tdp.s), 60), kb = boundaryLayerTdp({ EI: line.EI, H: Htdp, w: line.w, xi }).curvature;
    const shapes = [['Near', near], ['Mean', res], ['Far', far]].filter((c) => c[1]).map(([name, r]) => {
      const j0 = Math.max(0, (r.tdp?.node ?? n) - Math.round(n * 0.25)); const [xx, zz] = thin([Array.from(r.x.slice(j0)), Array.from(r.z.slice(j0))], 150);
      return { name, x: xx, y: zz, dash: name !== 'Mean' };
    });
    const vs = linspace(0, Math.max(1.5, Uc * 2), 25);
    const ucV = vs.map((v) => (v > 0 ? cableOnBottomStability({ submergedWeight: line.w, outerDiameter: line.D, velocity: v, mu: inp.muLat, Cd: inp.CdBed, Cl: inp.ClBed, rho: line.rhoW, safetyFactor: inp.sfStab }).requiredWeight / line.w : 0));
    const extraW = stab ? Math.max(0, stab.requiredWeight - line.w) : 0;
    return {
      kpis: [
        { key: 'seabedForce', label: 'Total seabed reaction', value: round(Fbed / 1e3, 6), unit: 'kN', explain: 'Integral of the normal contact reaction: the weight of the grounded length carried by the seabed in the mean position.' },
        { key: 'grounded', label: 'Grounded length (mean)', value: round(grounded), unit: 'm', explain: 'Length resting on the seabed between the touchdown point and the termination (resolved to one element).' },
        { key: 'groundedFar', label: 'Grounded length (far)', value: round(groundedFar), unit: 'm', limit: inp.minGrounded, sense: 'min', explain: 'Length still on the seabed with the vessel at its far offset; it must not shrink to the point where the termination is lifted.' },
        { key: 'peakReaction', label: 'Peak contact reaction', value: round(rMax), unit: 'N/m', explain: 'Largest seabed reaction per metre, found just behind the touchdown point where the bending boundary layer concentrates the load.' },
        { key: 'kLay', label: 'Touchdown load concentration', value: round(kLay, 3), unit: '-', explain: 'Peak reaction over submerged weight. Values above 1 mean the touchdown zone presses harder than the static weight and embeds deeper.' },
        { key: 'lambda', label: 'Bending boundary-layer length', value: round(bl.lambda, 3), unit: 'm', explain: 'Decay length √(EI/H) over which curvature rises from zero on the seabed to the catenary value; the element size should be of this order to resolve touchdown curvature.' },
        { key: 'kTdp', label: 'Catenary curvature at touchdown', value: round(bl.catenaryCurvature, 3), unit: '1/m', explain: 'Curvature w/H just above the touchdown point (radius = H/w). Compare its inverse with the dynamic minimum bend radius.' },
        { key: 'tdpTravel', label: 'Touchdown travel along the line (near → far)', value: round(travelS), unit: 'm', explain: 'Length of umbilical that alternately lifts off and lands as the vessel moves between near and far: the sheath-abrasion and trenching zone.' },
        { key: 'tdpTravelX', label: 'Touchdown movement on the seabed', value: round(travelX), unit: 'm', explain: 'Horizontal shift of the touchdown point between the near and far positions.' },
        { key: 'frictionForce', label: 'Available axial friction on the grounded length (far)', value: round(fricAvail / 1e3), unit: 'kN', explain: 'Coulomb resistance μ·w·L of the grounded length in the far position: the tension it can shed before the termination once the line has moved enough to mobilise it.' },
        { key: 'termResidual', label: 'Termination load with friction fully mobilised', value: round(termResidual / 1e3), unit: 'kN', explain: 'Touchdown tension minus the available friction. Zero means friction alone can anchor the line; the check uses the pinned-end reaction, which assumes no friction relief.' },
        { key: 'anchorLen', label: 'Friction anchoring length (far)', value: round(Math.min(anchorLen, 1e6)), unit: 'm', explain: 'Grounded length needed for Coulomb friction μ·w to absorb the whole touchdown tension. If the grounded length is shorter, the termination carries the remainder.' },
        { key: 'termLoad', label: 'Horizontal load on the termination', value: round(termLoad / 1e3), unit: 'kN', limit: inp.termRating, sense: 'max', explain: 'Support reaction at the pinned seabed end from the rod model (largest of mean and far positions). Conservative: axial friction is mobilised only by movement, which a pinned end prevents.' },
        { key: 'embedment', label: 'Embedment at touchdown', value: round(zEmb * 1e3), unit: 'mm', explain: 'Penetration at which soil bearing capacity balances the concentrated touchdown load, including the dynamic factor. Deeper embedment raises lateral and axial resistance but hinders recovery.' },
        { key: 'stabUc', label: 'Lateral stability utilisation', value: round(checks[0].uc), unit: '-', limit: 1, sense: 'max', explain: 'Required over actual submerged weight in the design near-bed flow. Above 1.0 the grounded section is predicted to slide.' },
        { key: 'vCrit', label: 'Critical near-bed velocity', value: round(stab ? stab.criticalVelocity : 0, 3), unit: 'm/s', explain: 'Flow speed normal to the line at which sliding starts (without the safety factor).' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Seabed reaction and friction through the touchdown zone', x: { label: 'Arc length from hang-off', unit: 'm' }, y: { label: 'Force per unit length', unit: 'N/m' },
          series: [{ name: 'Normal reaction', x: sR, y: rR }, { name: 'Axial friction', x: sR, y: fR }, { name: 'Submerged weight', x: [sR[0], sR[sR.length - 1]], y: [line.w, line.w], dash: true }],
          explain: 'Contact reaction rises from zero at lift-off to a peak just behind the touchdown point and settles to the submerged weight. Friction is mobilised where the grounded line is being pulled; where it falls to zero the line is anchored.' },
        { type: 'line', title: 'Curvature at touchdown: rod model and boundary-layer solution', x: { label: 'Distance above the touchdown point', unit: 'm' }, y: { label: 'Curvature', unit: '1/m' },
          series: [{ name: 'Boundary layer κ = (w/H)(1 − e^(−ξ/λ))', x: xi, y: kb }, { name: 'Rod model', x: sR.map((v) => tdM.s - v).filter((v) => v >= 0 && v <= xi[xi.length - 1]), y: kR.filter((_, i) => tdM.s - sR[i] >= 0 && tdM.s - sR[i] <= xi[xi.length - 1]) }],
          explain: 'Curvature builds up from zero at the seabed to the catenary value over a few boundary-layer lengths. Agreement between the two curves shows that the mesh resolves the touchdown bend.' },
        { type: 'line', title: 'Touchdown region in the near, mean and far positions', x: { label: 'Horizontal position', unit: 'm' }, y: { label: 'Elevation above seabed', unit: 'm' }, series: shapes,
          explain: 'Lower part of the line for each vessel position. The spread of the points where the curves meet the seabed is the touchdown travel.' },
        { type: 'line', title: 'Lateral stability vs near-bed velocity', x: { label: 'Velocity normal to the line', unit: 'm/s' }, y: { label: 'Required / actual submerged weight', unit: '-' },
          series: [{ name: 'Utilisation', x: vs, y: ucV }, { name: 'Limit', x: [0, vs[vs.length - 1]], y: [1, 1], dash: true }, { name: 'Design velocity', x: [Uc, Uc], y: [0, Math.max(1, checks[0].uc)], dash: true }],
          explain: 'Utilisation grows with the square of the flow speed. The crossing with the limit line is the largest near-bed velocity the as-laid umbilical tolerates without burial, mattresses or rock.' },
      ],
      tables: [{ title: 'Touchdown by vessel position', explain: 'One row per position: arc length and horizontal position of the touchdown point, grounded length, tension at touchdown and horizontal load on the termination.',
        columns: [{ key: 'pos', label: 'Position' }, { key: 's', label: 'Touchdown arc length', unit: 'm' }, { key: 'x', label: 'Touchdown x', unit: 'm' }, { key: 'g', label: 'Grounded length', unit: 'm' }, { key: 'T', label: 'Tension at touchdown', unit: 'kN' }, { key: 'H', label: 'Termination H', unit: 'kN' }],
        rows: [['Near', near], ['Mean', res], ['Far', far]].filter((c) => c[1]).map(([pos, r]) => ({ pos, s: r.tdp ? round(td(r).s) : null, x: r.tdp ? round(td(r).x) : null, g: r.tdp ? round(L - td(r).s) : 0, T: r.tdp ? round(r.Te[r.tdp.node] / 1e3) : null, H: round(Math.abs(r.endB.Fx) / 1e3) })) }],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        `Peak reaction ${fmt(rMax)} N/m at ${fmt(sPeak, 4)} m from the hang-off is ${fmt(kLay)} × the submerged weight; the point reaction of the rigid-seabed boundary layer is ${fmt(bl.pointReaction)} N.`,
        embNote || `Embedment ${fmt(zEmb * 1e3)} mm (${fmt(zEmb / line.D * 100)} % of the diameter) on ${soil.name.toLowerCase()} with a dynamic factor of ${inp.kDyn}.`,
        Number.isFinite(anchorLen) && anchorLen < groundedFar ? `Once mobilised, friction can anchor the line ${fmt(anchorLen)} m behind the touchdown point (grounded ${fmt(groundedFar)} m); the pinned-end reaction of ${fmt(termLoad / 1e3)} kN used in the check is therefore an upper bound.` : `The grounded length (${fmt(groundedFar)} m) is shorter than the friction anchoring length (${fmt(Math.min(anchorLen, 1e6))} m): the termination carries at least ${fmt(termResidual / 1e3)} kN.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}. Touchdown travels ${fmt(travelS)} m along the line and the far position leaves ${fmt(groundedFar)} m grounded.`,
          gov.uc > 1 ? (gov.name.startsWith('Lateral') ? `Add ${fmt(extraW)} N/m of submerged weight over the exposed length (ballast, heavier armour) or protect it by trenching, mattresses or rock; alternatively orient the route along the dominant near-bed flow.` : 'Lengthen the grounded section (move the termination away or add line), add a hold-back anchor or clump weight behind the touchdown zone, or raise the termination rating.')
            : 'Specify abrasion protection (sacrificial sheath or wear sleeve) over the touchdown travel plus margin, and confirm soil parameters from the route survey before relying on the embedment and friction values.',
          'DNV-RP-F114 (2017) §4; DNV-RP-F109 (2021) §3; site data per ISO 19901-8'),
        adv('economic', checks[0].uc > 1 ? 'warn' : 'info', checks[0].uc > 1 ? `The grounded section is ${fmt((checks[0].uc - 1) * 100)} % short of the weight needed for stability: stabilisation is required over the exposed route.` : `The as-laid umbilical is stable up to ${fmt(stab ? stab.criticalVelocity : 0)} m/s near-bed flow, so no stabilisation is needed for the specified ${fmt(Uc)} m/s.`,
          checks[0].uc > 1 ? 'Compare post-lay trenching (typically several hundred metres per hour for a jet trencher in soft clay) with spot mattresses on vessel-day basis; trenching also removes the fishing-interaction risk.' : 'Avoid specifying trenching or rock for stability alone on this section; reserve the intervention budget for crossings and third-party protection.',
          'AACE RP 18R-97 (quantity and vessel-day based estimate); DNV-RP-F109 stabilisation options',
          `${vesselCo2Text(2)} Each avoided intervention spread also avoids seabed disturbance over the trench or berm footprint; trenching in soft clay re-suspends sediment along the whole route.`),
      ],
      publish: { 'umbilical.tdpTravel': travelS, 'umbilical.groundedLength': grounded, 'umbilical.embedment': zEmb },
      quantities: { cable_km: grounded / 1e3 },
    };
  },
  verify: [
    // Frictionless catenary on a flat seabed (Irvine §1.2): h = 340 m, L = 600 m, X = 330 m → a = 32.00 m,
    // suspended S = sqrt(h(h + 2a)) = sqrt(340 × 404) = 370.62 m, grounded = 229.38 m, seabed force = w × 229.38 =
    // 235.11 × 229.38 = 53.93 kN.
    { name: 'Catenary with touchdown: seabed reaction = w·(L − S), S = √(h(h + 2a))',
      source: 'Irvine, Cable Structures (1981) §1.2; hand solution a = 32.00 m',
      inputs: { config: 'catenary', length: 600, span: 330, offsetPct: 0, EI: 200, muAx: 0, nElem: 320, minGrounded: 1 },
      expect: [{ kpi: 'seabedForce', value: 53.93, tol: 0.015 }, { kpi: 'grounded', value: 229.38, tol: 0.02 }] },
  ],
};
