// Pile foundation and soil–structure interaction: axial capacity (API), lateral response on p–y springs, axial
// load–settlement on t–z / Q–z springs, the pile-head stiffness matrix handed to the jacket model, and the load
// sharing of the pile group under base shear and overturning moment.
import { need, run, fmt, round, adv, check, steelCo2, CONFIRM, BASIS_COST } from './_common.js';
import { createProfile, applyScour, validateLayers } from '../../engine/geotech/soilprofile.js';
import { axialCapacity, capacityVsPenetration, requiredPenetration } from '../../engine/geotech/axialcapacity.js';
import { solveLateralPile, pileHeadStiffness, axialLoadSettlement, steelUtilisation } from '../../engine/geotech/pile.js';
import { analysePileGroup, pMultiplier } from '../../engine/geotech/pilegroup.js';

export default {
  id: 's01-pile-foundation', suite: 1, title: 'Pile foundation and soil–structure interaction',
  summary: 'Axial capacity and required penetration, lateral deflection, bending moment and soil reaction along the pile, load–settlement response, pile-head stiffness for the jacket model and pile-group load sharing.',
  covers: [13], fidelity: 'numerical',
  standards: ['API RP 2GEO (2011, add. 2014) §8.1–8.5', 'API RP 2A-WSD 22nd ed. §9 (21st ed. §6.4–6.8)', 'ISO 19901-4:2016 §8', 'ISO 19902:2020 §17'],
  models: ['Beam on nonlinear Winkler foundation (p–y curves: Matlock soft clay, O\'Neill–Murchison sand)', 'Axial shaft friction (α method for clay, β method for sand) and end bearing', 't–z and Q–z load transfer curves', 'Pile-head stiffness by static condensation (secant at working load)', 'Rigid-cap pile group stiffness method', 'p-multipliers for group shadowing', 'Local and general scour'],
  limits: 'Single vertical pipe pile with soil springs uncoupled between depths (Winkler); soil layering corrections, cyclic degradation beyond the cyclic p–y curves, set-up, pile driving and installation stresses are not modelled, and continuum soil models (Mohr–Coulomb, Cam-Clay finite elements) are outside the module. Safety factors and soil parameters are criteria to be confirmed against the governing edition and the site investigation.',
  publishes: ['foundation.Kh', 'foundation.Kr', 'foundation.Khr', 'foundation.Kv'],
  searchTerms: 'offshore pile p-y curves axial capacity API RP 2GEO pile head stiffness',
  inputs: [
    { key: 'pileD', label: 'Pile outer diameter', unit: 'mm', type: 'number', default: 1829, min: 300, max: 5000, group: 'Pile', role: 'geometry', design: true },
    { key: 'pileT', label: 'Pile wall thickness', unit: 'mm', type: 'number', default: 60, min: 10, max: 150, group: 'Pile', role: 'geometry', design: true },
    { key: 'pileL', label: 'Penetration below mudline', unit: 'm', type: 'number', default: 70, min: 5, max: 200, group: 'Pile', role: 'geometry', design: true },
    { key: 'mat', label: 'Pile steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Pile', role: 'material' },
    { key: 'layers', label: 'Soil layers (from the mudline down)', type: 'table', group: 'Soil', role: 'soil',
      columns: [{ key: 'type', label: 'Type (clay / sand)', type: 'text' }, { key: 'thickness', label: 'Thickness', unit: 'm', type: 'number' }, { key: 'gammaSub', label: 'Submerged unit weight', unit: 'kN/m³', type: 'number' }, { key: 'su', label: 'su at top (clay)', unit: 'kPa', type: 'number' }, { key: 'suGrad', label: 'su gradient (clay)', unit: 'kPa/m', type: 'number' }, { key: 'phi', label: 'φ′ (sand)', unit: 'deg', type: 'number' }, { key: 'eps50', label: 'ε50 (clay)', unit: '-', type: 'number' }],
      default: [{ type: 'clay', thickness: 12, gammaSub: 8, su: 50, suGrad: 4, phi: 0, eps50: 0.008 }, { type: 'sand', thickness: 20, gammaSub: 9.5, su: 0, suGrad: 0, phi: 35, eps50: 0 }, { type: 'clay', thickness: 68, gammaSub: 9, su: 120, suGrad: 2.5, phi: 0, eps50: 0.006 }],
      help: 'Stratigraphy from the site investigation. Clay rows use su and ε50, sand rows φ′; the last layer must reach below the pile tip.' },
    { key: 'scour', label: 'General scour', unit: 'm', type: 'number', default: 1.0, min: 0, max: 15, group: 'Soil', role: 'ic', help: 'Soil removed over the whole site; local scour of 1.5 pile diameters is added around the pile (API RP 2A commentary).' },
    { key: 'localScour', label: 'Include local scour of 1.5 D', type: 'bool', default: true, group: 'Soil', role: 'ic' },
    { key: 'loading', label: 'Lateral soil curves', type: 'select', default: 'cyclic', group: 'Soil', role: 'criteria', options: [{ value: 'cyclic', label: 'Cyclic (storm loading)' }, { value: 'static', label: 'Static' }] },
    { key: 'Nc', label: 'Pile-head compression', unit: 'MN', type: 'number', default: 25, min: 0, max: 500, group: 'Pile-head loads (boundary condition)', role: 'bc', bind: 'jacket.pileAxialMax', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Largest unfactored pile-head compression from the in-place analysis.' },
    { key: 'Nt', label: 'Pile-head tension', unit: 'MN', type: 'number', default: 2, min: 0, max: 500, group: 'Pile-head loads (boundary condition)', role: 'bc' },
    { key: 'Hh', label: 'Pile-head shear', unit: 'MN', type: 'number', default: 3, min: 0.001, max: 100, group: 'Pile-head loads (boundary condition)', role: 'bc', bind: 'jacket.pileShearMax', uncertain: { dist: 'normal', cov: 0.1 } },
    { key: 'Mh', label: 'Pile-head moment (free head only)', unit: 'MN·m', type: 'number', default: 0, min: -500, max: 500, group: 'Pile-head loads (boundary condition)', role: 'bc', showIf: { key: 'head', in: ['free'] } },
    { key: 'head', label: 'Pile-head rotation', type: 'select', default: 'fixed', group: 'Pile-head loads (boundary condition)', role: 'bc', options: [{ value: 'fixed', label: 'Restrained by the jacket leg (fixed rotation)' }, { value: 'free', label: 'Free head (shear and moment given)' }], help: 'A pile grouted or welded into a leg is close to rotation-fixed at the mudline.' },
    { key: 'FSc', label: 'Safety factor on axial capacity', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Design criteria', role: 'criteria', help: `API RP 2A-WSD: 1.5 for extreme, 2.0 for operating conditions. ${CONFIRM}` },
    { key: 'yLim', label: 'Allowable mudline deflection as a share of diameter', unit: '%', type: 'number', default: 10, min: 1, max: 30, group: 'Design criteria', role: 'criteria', help: 'Serviceability screen; large deflection means the p–y curves are near their ultimate over the upper pile.' },
    { key: 'usage', label: 'Allowable steel stress as a share of yield', unit: '-', type: 'number', default: 0.8, min: 0.4, max: 1, group: 'Design criteria', role: 'criteria', help: '0.6·1.33 = 0.8 for extreme conditions in working-stress design.' },
    { key: 'nPiles', label: 'Piles in the group', type: 'int', default: 4, min: 1, max: 16, group: 'Pile group', role: 'geometry', help: 'Placed on a square or rectangular grid at the leg spacing (one pile per leg for 4).' },
    { key: 'spacing', label: 'Pile spacing (leg spacing at the mudline)', unit: 'm', type: 'number', default: 37, min: 2, max: 150, group: 'Pile group', role: 'geometry' },
    { key: 'Vbase', label: 'Base shear on the group', unit: 'MN', type: 'number', default: 11, min: 0, max: 500, group: 'Pile group', role: 'bc', bind: 'jacket.baseShear' },
    { key: 'Mbase', label: 'Overturning moment on the group', unit: 'MN·m', type: 'number', default: 700, min: 0, max: 5e4, group: 'Pile group', role: 'bc', bind: 'jacket.otm' },
    { key: 'Wtot', label: 'Vertical load on the group', unit: 'MN', type: 'number', default: 50, min: 0, max: 5000, group: 'Pile group', role: 'bc', help: 'Topsides plus buoyant jacket weight.' },
    { key: 'nElem', label: 'Pile elements', type: 'int', default: 60, min: 10, max: 480, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Pile elements', levels: [30, 60, 120], qoi: 'headDefl' },
  compute(inp, ctx) {
    // SI conversion: mm → m, kPa → Pa, kN/m³ → N/m³, MN → N.
    const D = inp.pileD / 1e3, t = inp.pileT / 1e3, L = inp.pileL, mat = ctx.material(inp.mat);
    need(mat, `Steel '${inp.mat}' is not in the materials library.`);
    need(t < D / 2, 'Pile wall thickness must be smaller than the pile radius.');
    let top = 0;
    const layers = (inp.layers || []).map((r, i) => {
      const type = String(r.type || '').trim().toLowerCase(), th = Number(r.thickness), g = Number(r.gammaSub) * 1e3;
      need(type === 'clay' || type === 'sand', `Soil layer ${i + 1}: type must be "clay" or "sand".`);
      need(th > 0 && g > 0, `Soil layer ${i + 1}: thickness and submerged unit weight must be positive.`);
      const lay = type === 'clay' ? { type, top, thickness: th, gammaSub: g, suTop: Number(r.su) * 1e3, suGrad: (Number(r.suGrad) || 0) * 1e3, eps50: Number(r.eps50) > 0 ? Number(r.eps50) : 0.01 } : { type, top, thickness: th, gammaSub: g, phiDeg: Number(r.phi), Dr: Math.min(0.95, Math.max(0.36, 0.04 * (Number(r.phi) - 20))) };
      need(type !== 'clay' || lay.suTop > 0, `Soil layer ${i + 1}: clay needs a positive undrained shear strength.`);
      need(type !== 'sand' || (lay.phiDeg >= 29 && lay.phiDeg <= 45), `Soil layer ${i + 1}: sand friction angle must be between 29° and 45°.`);
      top += th; return lay;
    });
    need(layers.length, 'Give at least one soil layer.');
    need(top >= L + 2, `The soil layers reach ${fmt(top)} m but the pile tip is at ${fmt(L)} m: extend the last layer at least 2 m below the tip.`);
    const v = validateLayers(layers); need(v.valid, `Soil layers: ${v.errors.join('; ')}`);
    let profile = createProfile(layers);
    const loc = inp.localScour ? 1.5 * D : 0;
    if (inp.scour > 0 || loc > 0) { need(inp.scour + loc < 0.5 * L, 'Scour removes more than half the penetration: increase the penetration.'); profile = run(() => applyScour(profile, { general: inp.scour, local: loc, D }), 'Scour'); }
    // axial capacity and penetration
    const ax = run(() => axialCapacity({ profile, L, D, t }), 'Axial capacity'), Nc = inp.Nc * 1e6, Nt = inp.Nt * 1e6, H = inp.Hh * 1e6;
    const cv = run(() => capacityVsPenetration({ profile, D, t, Lmin: Math.max(5, Math.ceil(inp.scour + loc + 2)), Lmax: Math.min(top - 2, 1.5 * L), step: Math.max(1, Math.round(L / 40)) }), 'Capacity curve');
    const rq = Nc > 0 ? run(() => requiredPenetration({ profile, D, t, load: Nc, mode: 'compression', safetyFactor: inp.FSc, Lmin: 5, Lmax: top - 2 }), 'Required penetration') : { L: 0, found: true };
    // lateral response at the working load
    const o = { L, D, t, E: mat.E, nElem: inp.nElem, profile, loading: inp.loading, axialLoad: Nc, head: inp.head === 'fixed' ? { H, fixedRotation: true } : { H, M: inp.Mh * 1e6 } };
    const lat = run(() => solveLateralPile(o), 'Lateral pile analysis (the head load may exceed the lateral soil capacity)');
    need(lat.convergence.converged, 'The lateral pile analysis did not converge: the head load is close to the lateral capacity. Increase diameter, wall or penetration.');
    const hs = run(() => pileHeadStiffness(o), 'Pile-head stiffness'), K = { Kh: hs.secant.KHH, Kr: hs.secant.KMM, Khr: Math.abs(hs.secant.KHM) };
    const ls = run(() => axialLoadSettlement({ profile, L, D, t, E: mat.E, nElem: Math.max(20, Math.round(inp.nElem / 2)) }, { steps: 40 }), 'Axial load–settlement'), Kv = ls.headStiffness;
    need([K.Kh, K.Kr, K.Khr, Kv].every((x) => Number.isFinite(x) && x > 0), 'The pile-head stiffness could not be evaluated: check soil strength and pile dimensions.');
    const y0 = Math.abs(lat.head.y), Mmax = Math.abs(lat.maxMoment), su = steelUtilisation({ N: Nc, M: Mmax, D, t, fy: mat.fy, usageFactor: inp.usage });
    // settlement at the working compression from the load–settlement curve
    let wN = 0; for (let i = 1; i < ls.Q.length; i++) if (ls.Q[i] >= Nc) { wN = ls.w[i - 1] + (ls.w[i] - ls.w[i - 1]) * (Nc - ls.Q[i - 1]) / (ls.Q[i] - ls.Q[i - 1]); break; } if (!wN && Nc > 0) wN = ls.wAtCapacity;
    // pile group on a rigid cap: near-square grid at the given spacing, load along x
    const nx = Math.max(1, Math.round(Math.sqrt(inp.nPiles))), ny = Math.ceil(inp.nPiles / nx), pos = [];
    for (let i = 0; i < nx; i++) for (let j = 0; j < ny && pos.length < inp.nPiles; j++) pos.push({ x: (i - (nx - 1) / 2) * inp.spacing, y: (j - (ny - 1) / 2) * inp.spacing });
    const sD = inp.spacing / D, pm = inp.nPiles > 1 ? [1, 2].map((row) => pMultiplier({ spacingOverD: sD, row })) : [1, 1];
    const grp = run(() => analysePileGroup({ piles: pos.map((p) => ({ ...p, stiffness: { KHH: K.Kh, KHM: -K.Khr, KMM: K.Kr, Kaxial: Kv, Ktorsion: 0.25 * K.Kr, head: 'fixed' } })), load: { Fx: inp.Vbase * 1e6, Fz: -inp.Wtot * 1e6, My: inp.Mbase * 1e6 } }), 'Pile group analysis');
    const gAx = Math.max(...grp.piles.map((p) => p.axial)), gTen = Math.max(0, ...grp.piles.map((p) => -p.axial)), gSh = Math.max(...grp.piles.map((p) => p.shear));
    const steel = Math.PI * (D - t) * t * L * (mat.rho ?? 7850) / 1e3 * inp.nPiles, perM = steel / L;
    const checks = [
      check('Axial compression capacity', Nc * inp.FSc / 1e6, ax.Qc / 1e6, 'MN', `API RP 2GEO §8.1; API RP 2A-WSD §6.3.4 (FS = ${inp.FSc})`, 'Pile-head compression times the safety factor against shaft friction plus end bearing. If exceeded the pile plunges: settlement grows without limit under the storm load.'),
      check('Axial tension (pull-out) capacity', Nt * inp.FSc / 1e6, ax.Qt / 1e6, 'MN', `API RP 2GEO §8.1; shaft friction only (FS = ${inp.FSc})`, 'Pile-head uplift times the safety factor against shaft friction. If exceeded the pile is pulled out and the jacket overturns.'),
      check('Pile steel stress at maximum moment', su.sigmaMax / 1e6, inp.usage * mat.fy / 1e6, 'MPa', 'API RP 2A-WSD §6.9–6.10 (axial + bending, first yield with usage factor)', 'Axial plus bending stress in the pile wall at the depth of maximum moment. If exceeded the pile yields below the mudline, where it cannot be inspected or repaired.'),
      check('Mudline lateral deflection', y0 * 1e3, inp.yLim / 100 * D * 1e3, 'mm', 'Serviceability criterion (project); p–y curves per API RP 2GEO §8.5', 'Lateral movement of the pile head. Large values mean the upper soil is at its ultimate resistance and the foundation stiffness assumed in the jacket analysis is degrading.'),
      check('Lateral soil mobilisation (peak p/pu)', lat.maxMobilisation, 1, '-', 'API RP 2GEO §8.5 (ultimate lateral resistance)', 'Largest ratio of soil reaction to ultimate lateral resistance along the pile. A value of 1.0 over a growing depth is the approach to lateral failure.'),
    ];
    const ok = checks.every((c) => c.uc <= 1), gov = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), stp = Math.max(1, Math.ceil(lat.z.length / 150)), ix = Array.from(lat.z, (z, i) => i).filter((i) => i % stp === 0 || i === lat.z.length - 1);
    return {
      kpis: [
        { key: 'Qc', label: 'Ultimate compression capacity', unit: 'MN', value: ax.Qc / 1e6, explain: 'Outer shaft friction plus the smaller of plugged and unplugged end bearing at the given penetration, before any safety factor.' },
        { key: 'Qshaft', label: 'Outer shaft friction', unit: 'MN', value: ax.Qfo / 1e6, explain: 'Friction on the outside of the pile integrated over the penetration: the whole of the tension capacity and usually most of the compression capacity.' },
        { key: 'Qt', label: 'Ultimate tension capacity', unit: 'MN', value: ax.Qt / 1e6, explain: 'Pull-out resistance from shaft friction (pile and plug weight not included).' },
        { key: 'fsAxial', label: 'Factor of safety in compression', unit: '-', value: Nc > 0 ? Math.min(ax.Qc / Nc, 999) : 999, limit: inp.FSc, sense: 'min', explain: 'Ultimate capacity over the pile-head compression. API RP 2A asks for 1.5 in extreme and 2.0 in operating conditions.' },
        { key: 'Lreq', label: 'Penetration required for compression', unit: 'm', value: rq.L, explain: `Shortest penetration at which capacity divided by ${fmt(inp.FSc)} equals the compression load${rq.found ? '' : ' (not reached within the soil profile)'}.` },
        { key: 'headDefl', label: 'Mudline lateral deflection', unit: 'mm', value: y0 * 1e3, limit: inp.yLim / 100 * D * 1e3, sense: 'max', explain: 'Lateral displacement of the pile at the original mudline under the head shear (and moment), on nonlinear p–y springs.' },
        { key: 'Mmax', label: 'Maximum pile bending moment', unit: 'MN·m', value: Mmax / 1e6, explain: `Largest bending moment in the pile, at ${fmt(lat.depthMaxMoment)} m below the mudline; it sizes the wall of the upper pile section.` },
        { key: 'settlement', label: 'Pile-head settlement at working load', unit: 'mm', value: wN * 1e3, explain: 'Axial displacement of the pile head under the compression load from the t–z / Q–z analysis.' },
        { key: 'Kh', label: 'Lateral pile-head stiffness', unit: 'MN/m', value: K.Kh / 1e6, explain: 'Secant lateral stiffness at the working load: the Kh of the pile-head spring in the jacket model.' },
        { key: 'Kr', label: 'Rotational pile-head stiffness', unit: 'GN·m/rad', value: K.Kr / 1e9, explain: 'Secant rocking stiffness of the pile head.' },
        { key: 'Khr', label: 'Lateral–rotational coupling', unit: 'MN', value: K.Khr / 1e6, explain: 'Magnitude of the off-diagonal term of the pile-head stiffness matrix (negative in the pile sign convention).' },
        { key: 'Kv', label: 'Axial pile-head stiffness', unit: 'MN/m', value: Kv / 1e6, explain: 'Initial secant stiffness of the load–settlement curve: the vertical spring under each leg.' },
        { key: 'groupAxial', label: 'Largest pile compression in the group', unit: 'MN', value: gAx / 1e6, explain: 'Compression of the most loaded pile of a rigid-cap group under the vertical load, base shear and overturning moment. Compare with the pile-head compression used above.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Lateral deflection along the pile', x: { label: 'Deflection', unit: 'mm' }, y: { label: 'Depth below mudline', unit: 'm' }, series: [{ name: 'y', x: ix.map((i) => round(lat.y[i] * 1e3)), y: ix.map((i) => round(-lat.z[i])) }], explain: 'Deflected shape of the pile. The depth of the first zero crossing shows how much of the pile works laterally; anything below about two such lengths only adds axial capacity.' },
        { type: 'line', title: 'Bending moment and shear along the pile', x: { label: 'Moment (MN·m) / shear (MN)', unit: '-' }, y: { label: 'Depth below mudline', unit: 'm' }, series: [{ name: 'Moment', x: ix.map((i) => round(lat.moment[i] / 1e6)), y: ix.map((i) => round(-lat.z[i])) }, { name: 'Shear', x: ix.map((i) => round(lat.shear[i] / 1e6)), y: ix.map((i) => round(-lat.z[i])) }], explain: 'Moment and shear in the pile. A fixed head has its peak moment at the mudline; a free head a few diameters below. The moment diagram defines where the heavy-wall section can stop.' },
        { type: 'line', title: 'Soil reaction and ultimate resistance', x: { label: 'Lateral soil reaction', unit: 'kN/m' }, y: { label: 'Depth below mudline', unit: 'm' }, series: [{ name: 'Mobilised p', x: ix.map((i) => round(lat.p[i] / 1e3)), y: ix.map((i) => round(-lat.z[i])) }, { name: 'Ultimate pu', x: ix.map((i) => round(lat.pu[i] / 1e3)), y: ix.map((i) => round(-lat.z[i])), dash: true }], explain: 'Soil reaction along the pile against the ultimate lateral resistance. Where the two curves touch the soil has yielded; the scoured zone at the top gives no support.' },
        { type: 'line', title: 'Axial capacity against penetration', x: { label: 'Capacity', unit: 'MN' }, y: { label: 'Penetration', unit: 'm' }, series: [{ name: 'Compression', x: cv.Qc.map((q) => round(q / 1e6)), y: cv.L.map((l) => round(-l)) }, { name: 'Tension', x: cv.Qt.map((q) => round(q / 1e6)), y: cv.L.map((l) => round(-l)) }, { name: 'Required (load × FS)', x: [round(Nc * inp.FSc / 1e6), round(Nc * inp.FSc / 1e6)], y: [round(-cv.L[0]), round(-cv.L[cv.L.length - 1])], dash: true }], explain: 'Ultimate compression and tension capacity as the pile is driven deeper. The penetration is adequate where the compression curve passes the dashed line; steps mark stronger layers.' },
        { type: 'line', title: 'Pile-head load–settlement', x: { label: 'Settlement', unit: 'mm' }, y: { label: 'Head load', unit: 'MN' }, series: [{ name: 't–z / Q–z', x: ls.w.map((w) => round(w * 1e3)), y: ls.Q.map((q) => round(q / 1e6)) }, { name: 'Working load', x: [0, round(ls.w[ls.w.length - 1] * 1e3)], y: [round(Nc / 1e6), round(Nc / 1e6)], dash: true }], explain: 'Axial response of the pile head. The initial slope is the axial stiffness used in the jacket model; the plateau is the capacity.' },
      ],
      tables: [
        { title: 'Pile-head stiffness matrix for the jacket model', columns: [{ key: 'c', label: 'Term' }, { key: 'sec', label: 'Secant at working load' }, { key: 'tan', label: 'Tangent' }, { key: 'u', label: 'Unit' }],
          rows: [{ c: 'Kh (lateral)', sec: round(K.Kh / 1e6), tan: round(hs.tangent.KHH / 1e6), u: 'MN/m' }, { c: 'Kr (rotational)', sec: round(K.Kr / 1e9), tan: round(hs.tangent.KMM / 1e9), u: 'GN·m/rad' }, { c: 'Khr (coupling, magnitude)', sec: round(K.Khr / 1e6), tan: round(Math.abs(hs.tangent.KHM) / 1e6), u: 'MN' }, { c: 'Kv (axial)', sec: round(Kv / 1e6), tan: round(Kv / 1e6), u: 'MN/m' }],
          explain: 'Linearised pile-head springs. Use the secant values with the load level they were derived at; iterate with the in-place analysis if the pile-head loads change by more than about 10 %.' },
        { title: 'Axial capacity build-up', columns: [{ key: 'c', label: 'Component' }, { key: 'v', label: 'Value', unit: 'MN' }], rows: [{ c: 'Outer shaft friction', v: round(ax.Qfo / 1e6) }, { c: 'Inner shaft friction (plug)', v: round(ax.Qfi / 1e6) }, { c: 'End bearing, plugged', v: round(ax.QpPlugged / 1e6) }, { c: 'End bearing, unplugged (annulus + inner friction)', v: round(ax.QpUnplugged / 1e6) }, { c: `End bearing used (${ax.plugged ? 'plugged' : 'unplugged'})`, v: round(ax.Qp / 1e6) }, { c: 'Compression capacity', v: round(ax.Qc / 1e6) }, { c: 'Tension capacity', v: round(ax.Qt / 1e6) }],
          explain: 'Components of the ultimate axial capacity. The pile behaves as plugged when the inner friction exceeds the end bearing on the soil plug.' },
        { title: 'Pile group under base shear and overturning', columns: [{ key: 'id', label: 'Pile' }, { key: 'x', label: 'x', unit: 'm' }, { key: 'y', label: 'y', unit: 'm' }, { key: 'ax', label: 'Axial (compression +)', unit: 'MN' }, { key: 'sh', label: 'Shear', unit: 'MN' }, { key: 'm', label: 'Head moment', unit: 'MN·m' }],
          rows: grp.piles.map((p, i) => ({ id: i + 1, x: round(pos[i].x), y: round(pos[i].y), ax: round(p.axial / 1e6), sh: round(p.shear / 1e6), m: round(p.moment / 1e6) })), explain: 'Load on each pile of a rigid-cap group with the stiffness above, load along x. Negative axial values are uplift. It is an independent check of the reactions from the jacket model.' },
      ],
      notes: [`Governing: ${gov.name} at ${fmt(gov.uc)}. Capacity ${fmt(ax.Qc / 1e6)} MN (${ax.plugged ? 'plugged' : 'unplugged'}) against ${fmt(Nc / 1e6)} MN × ${fmt(inp.FSc)}; required penetration ${rq.found ? `${fmt(rq.L)} m` : 'not reached in this profile'}.`,
        `Lateral: ${fmt(y0 * 1e3)} mm at the mudline, maximum moment ${fmt(Mmax / 1e6)} MN·m at ${fmt(lat.depthMaxMoment)} m, converged in ${lat.convergence.iterations} iterations; scour ${fmt(inp.scour)} m general + ${fmt(loc)} m local.`,
        `Group of ${inp.nPiles}: spacing ${fmt(sD)} diameters, p-multipliers ${fmt(pm[0])} (leading) and ${fmt(pm[1])} (trailing) — ${sD > 8 ? 'no group reduction is needed at this spacing' : 'apply these to the p–y curves of the trailing piles'}; largest group compression ${fmt(gAx / 1e6)} MN, tension ${fmt(gTen / 1e6)} MN, shear ${fmt(gSh / 1e6)} MN.`],
      advice: [
        ok ? adv('engineering', 'ok', `The pile is adequate: axial factor of safety ${fmt(Nc > 0 ? ax.Qc / Nc : 999)} (≥ ${fmt(inp.FSc)}), steel stress ${fmt(su.sigmaMax / 1e6)} MPa, mudline deflection ${fmt(y0 * 1e3)} mm.`, rq.found && rq.L < 0.85 * L ? `Compression needs only ${fmt(rq.L)} m of penetration: up to ${fmt(L - Math.max(rq.L, 0))} m could be saved if tension, lateral fixity and drivability allow. Feed the pile-head stiffness back into the in-place analysis and confirm the loads.` : 'Penetration is well matched to the load. Feed the pile-head stiffness back into the in-place analysis, repeat with the updated loads, and confirm drivability to this depth.', 'API RP 2GEO (2011) §8; API RP 2A-WSD §6.3.4; ISO 19901-4:2016 §8', `Each metre of penetration on ${inp.nPiles} piles is about ${steelCo2(perM)} Driving energy and vessel time scale with it as well.`)
          : adv('engineering', 'fail', `The pile does not pass: ${gov.name} at ${fmt(gov.uc)}.`, gov.name.startsWith('Axial') ? `Increase penetration to about ${rq.found ? fmt(Math.max(rq.L, L * 1.05)) : fmt(L * gov.uc)} m or the diameter; shaft friction grows with both.` : gov.name.startsWith('Pile steel') ? `Increase the wall in the upper ${fmt(2 * lat.depthMaxMoment + 5 * D)} m to about ${fmt(inp.pileT * gov.uc / 0.9)} mm.` : 'Increase the pile diameter (lateral stiffness ∝ D⁴·t) or protect against scour, which removes the stiffest support.', 'API RP 2GEO (2011) §8; API RP 2A-WSD §6.3.4', `Extending ${inp.nPiles} piles by ${fmt(Math.max((rq.found ? rq.L : L * 1.2) - L, 1))} m is about ${steelCo2(perM * Math.max((rq.found ? rq.L : L * 1.2) - L, 1))}`),
        adv('economic', 'info', `${inp.nPiles} piles of ${fmt(D * 1e3)} × ${fmt(t * 1e3)} mm and ${fmt(L)} m penetration weigh ${fmt(steel)} t below the mudline.`, 'Pile steel and driving time are a large share of installation cost: a site-specific soil investigation with cone tests at each leg typically allows shorter piles than conservative strength profiles, and avoids refusal or re-drive days of the installation vessel.', `${BASIS_COST}; ISO 19901-8:2014 (marine soil investigations)`),
      ],
      publish: { 'foundation.Kh': K.Kh, 'foundation.Kr': K.Kr, 'foundation.Khr': K.Khr, 'foundation.Kv': Kv },
      quantities: { steel_t: round(steel) },
    };
  },
  verify: [
    // Uniform clay, su = 100 kPa, γ′ = 8 kN/m³, D = 1.524 m, L = 60 m, no scour. API α method with ψ = su/σ′v:
    // 0–12.5 m (ψ > 1): f = 0.5·su·ψ^-0.25 → ∫ f dz = f(12.5)·12.5/1.25 = 50 kPa·10 m = 500 kN/m;
    // 12.5–50 m (0.25 ≤ ψ ≤ 1): f = 0.5·√(su·σ′v) = 14.142·√z kPa → ∫ = 14.142·(2/3)·(50^1.5 − 12.5^1.5) = 2916.7 kN/m;
    // 50–60 m (α capped at 1): f = su → 1000 kN/m. Qfo = π·1.524·(500 + 2916.7 + 1000) kN = 21.15 MN.
    { name: 'Shaft friction in uniform clay by the API α method (closed-form integral)', source: 'API RP 2GEO (2011) §8.1.3, Eq. 17–18 (hand integration)',
      inputs: { layers: [{ type: 'clay', thickness: 80, gammaSub: 8, su: 100, suGrad: 0, phi: 0, eps50: 0.01 }], scour: 0, localScour: false, pileL: 60, pileD: 1524, pileT: 50 },
      expect: [{ kpi: 'Qshaft', value: 21.146, tol: 0.005 }] },
  ],
};
