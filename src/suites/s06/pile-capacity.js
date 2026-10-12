// Axial (compression, tension) and lateral capacity of a dolphin pile.
import { axialCapacity, capacityVsPenetration } from '../../engine/geotech/axialcapacity.js';
import { solveLateralPile, bromsClay, bromsSand } from '../../engine/geotech/pile.js';
import { propertiesAt } from '../../engine/geotech/soilprofile.js';
import { need, run, fmt, round, cap, adv, check, pileInputs, soilInputs, pileProps, profileFrom, latOpts, pileSteelT, steelCo2Text, CONFIRM, NOT_HIFI, ECON_BASIS } from './_common.js';

/**
 * Lateral load at which the pile forms a plastic hinge (maximum moment = Mp) or the soil fails (no equilibrium):
 * load stepping on the nonlinear pile, with the last step interpolated on moment.
 */
function lateralLimit(base, P, fixedHead, e) {
  const arm0 = (base.stickUp || 0) + e + 3 * P.D;
  const dH = (fixedHead ? 2 : 1) * P.Mp / arm0 / 12;
  const H = [0], y = [0], M = [0];
  let mode = 'plastic hinge in the pile', Hult = null;
  for (let k = 1; k <= 60; k++) {
    const h = k * dH;
    let s;
    try { s = solveLateralPile({ ...base, head: fixedHead ? { H: h, fixedRotation: true } : { H: h, M: h * e } }); } catch { mode = 'soil failure (no lateral equilibrium)'; Hult = H[H.length - 1]; break; }
    const m = Math.abs(s.maxMoment);
    if (!(s.head.y > y[y.length - 1]) || !Number.isFinite(m)) { mode = 'soil failure (no lateral equilibrium)'; Hult = H[H.length - 1]; break; }
    if (m >= P.Mp) { const f = (P.Mp - M[M.length - 1]) / (m - M[M.length - 1]); Hult = H[H.length - 1] + f * (h - H[H.length - 1]); H.push(Hult); y.push(y[y.length - 1] + f * (s.head.y - y[y.length - 1])); M.push(P.Mp); break; }
    H.push(h); y.push(s.head.y); M.push(m);
  }
  if (Hult === null) { Hult = H[H.length - 1]; mode = 'load stepping limit reached without failure'; }
  return { Hult, mode, H, y, M };
}

export default {
  id: 's06-pile-capacity', suite: 6, title: 'Pile lateral and axial capacity',
  summary: 'Ultimate compression, tension and lateral capacity of a dolphin pile with factors of safety, and the penetration needed for the design loads.',
  covers: [4], fidelity: 'numerical',
  standards: ['API RP 2GEO (2011, add. 2014) §8.1–8.5', 'API RP 2A-WSD 22nd ed. §6.3–6.8', 'ISO 19901-4:2016 §8', 'Broms (1964) lateral resistance of piles'],
  models: ['API α method (clay) and β method (sand) shaft friction', 'End bearing with plug check', 'Lateral pile capacity by load stepping to a plastic hinge on p–y springs', 'Broms short / long pile limit for reference'],
  limits: 'Static capacity of a driven open-ended steel pipe pile in siliceous soils; pile and plug weight, cyclic axial degradation, set-up and group block failure are not included (apply them through the capacity factors). Lateral capacity assumes an elastic pile up to the first plastic hinge; post-hinge redistribution is in the integrity module. ' + NOT_HIFI,
  publishes: ['dolphin.pileQc', 'dolphin.pileQt', 'dolphin.pileHult', 'dolphin.LembRequired'],
  searchTerms: 'pile axial capacity API alpha beta method lateral capacity Broms',
  inputs: [
    ...pileInputs(),
    { key: 'headFix', label: 'Pile head', type: 'select', default: 'fixed', group: 'Pile', role: 'bc', options: [{ value: 'fixed', label: 'Restrained against rotation by the cap' }, { value: 'free', label: 'Free (single pile or pinned)' }] },
    { key: 'eLoad', label: 'Load height above the pile head (free head)', unit: 'm', type: 'number', default: 1.5, min: 0, max: 20, group: 'Pile', role: 'bc', showIf: { key: 'headFix', in: ['free'] } },
    ...soilInputs(),
    { key: 'plugMode', label: 'End bearing', type: 'select', default: 'auto', group: 'Soil', role: 'soil', options: [{ value: 'auto', label: 'Smaller of plugged and coring' }, { value: 'plugged', label: 'Plugged' }, { value: 'unplugged', label: 'Coring (annulus + internal friction)' }] },
    { key: 'tensionFactor', label: 'Tension / compression shaft friction ratio', unit: '-', type: 'number', default: 1.0, min: 0.3, max: 1, group: 'Soil', role: 'soil', help: 'API uses the same friction in tension as in compression; 0.7–0.8 is often applied in sand.' },
    { key: 'cyclicFactor', label: 'Capacity factor for cyclic loading and set-up', unit: '-', type: 'number', default: 1.0, min: 0.3, max: 1.5, group: 'Soil', role: 'soil', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Multiplier on shaft friction: below 1 for degradation under repeated berthing loads, above 1 for demonstrated set-up.' },
    { key: 'Nc', label: 'Design pile compression', unit: 'kN', type: 'number', default: 3200, min: 0, max: 2e5, group: 'Design loads', role: 'bc', bind: 'dolphin.pileAxialMax', uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'Nt', label: 'Design pile tension', unit: 'kN', type: 'number', default: 2400, min: 0, max: 2e5, group: 'Design loads', role: 'bc', bind: 'dolphin.pileTensionMax', uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'Hd', label: 'Design pile shear', unit: 'kN', type: 'number', default: 480, min: 0, max: 1e5, group: 'Design loads', role: 'bc', bind: 'dolphin.pileShearMax' },
    { key: 'fsC', label: 'Factor of safety, compression', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Acceptance', role: 'criteria', help: `API RP 2A-WSD §6.3.4: 1.5 extreme, 2.0 operating. ${CONFIRM}` },
    { key: 'fsT', label: 'Factor of safety, tension', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Acceptance', role: 'criteria', help: CONFIRM },
    { key: 'fsH', label: 'Factor of safety, lateral', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Acceptance', role: 'criteria', help: `On the load at first plastic hinge or soil failure. ${CONFIRM}` },
    { key: 'nElem', label: 'Elements over the embedded length', unit: '-', type: 'int', default: 40, min: 10, max: 320, group: 'Numerical', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Pile elements (embedded length)', levels: [20, 40, 80], qoi: 'Hult' },
  compute(inp, ctx) {
    const P = pileProps(inp, ctx);
    const profile = profileFrom(inp, P.D);
    const axO = { profile, D: P.D, t: P.t, plugMode: inp.plugMode, tensionFactor: inp.tensionFactor, setupFactor: inp.cyclicFactor };
    const ax = run(() => axialCapacity({ ...axO, L: inp.Lemb }), 'Axial capacity');
    const Lmax = Math.min(profile.depth - 0.5, Math.max(inp.Lemb * 1.6, inp.Lemb + 10));
    const cv = run(() => capacityVsPenetration({ ...axO, Lmin: Math.max(3, Math.min(5, inp.Lemb / 2)), Lmax, step: Math.max(0.5, (Lmax - 3) / 40) }), 'Capacity against penetration');
    const fixed = inp.headFix === 'fixed';
    const lat = lateralLimit(latOpts(inp, P, profile), P, fixed, fixed ? 0 : inp.eLoad);
    const fsC = inp.Nc > 0 ? ax.Qc / (inp.Nc * 1e3) : Infinity, fsT = inp.Nt > 0 ? ax.Qt / (inp.Nt * 1e3) : Infinity, fsH = inp.Hd > 0 ? lat.Hult / (inp.Hd * 1e3) : Infinity;
    let Lreq = NaN;
    for (let i = 0; i < cv.L.length; i++) if (cv.Qc[i] >= inp.fsC * inp.Nc * 1e3 && cv.Qt[i] >= inp.fsT * inp.Nt * 1e3) { Lreq = cv.L[i]; break; }
    // Broms reference with the soil over the top eight diameters below the scoured surface
    const zs = (inp.scour ?? 0) + 4 * P.D, top = propertiesAt(profile, Math.min(zs, profile.depth - 0.1));
    let broms = null;
    try {
      broms = top.type === 'clay' ? bromsClay({ su: top.su, D: P.D, L: inp.Lemb - (inp.scour ?? 0), e: inp.Lfree + (inp.scour ?? 0) + (fixed ? 0 : inp.eLoad), My: P.Mp })
        : bromsSand({ gammaSub: top.gammaSub, phiDeg: top.phiDeg, D: P.D, L: inp.Lemb - (inp.scour ?? 0), e: inp.Lfree + (inp.scour ?? 0) + (fixed ? 0 : inp.eLoad), My: P.Mp });
    } catch { broms = null; }
    const checks = [
      check('Axial compression', inp.Nc, ax.Qc / inp.fsC / 1e3, 'kN', 'API RP 2GEO (2011) §8.1; API RP 2A-WSD §6.3.4', 'Pile compression against shaft friction plus end bearing over the factor of safety. If exceeded the pile plunges and the cap tilts.'),
      check('Axial tension (pull-out)', inp.Nt, ax.Qt / inp.fsT / 1e3, 'kN', 'API RP 2GEO (2011) §8.1; API RP 2A-WSD §6.3.4', 'Pile tension against shaft friction alone over the factor of safety. If exceeded the rear piles of the group pull out progressively under repeated berthing.'),
      check('Lateral load', inp.Hd, lat.Hult / inp.fsH / 1e3, 'kN', 'API RP 2GEO (2011) §8.5; Broms (1964)', 'Pile shear against the lateral load that forms a plastic hinge or fails the soil. If exceeded the pile is permanently bent.'),
    ];
    const w = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = w.uc <= 1;
    const dL = Number.isFinite(Lreq) ? Lreq - inp.Lemb : NaN;
    const tPerM = pileSteelT({ ...inp, Lfree: 0, Lemb: 1 }, 1);
    const advice = [
      adv('engineering', ok ? 'ok' : 'fail', `Capacity: compression ${fmt(ax.Qc / 1e3)} kN (FS ${fmt(cap(fsC))}), tension ${fmt(ax.Qt / 1e3)} kN (FS ${fmt(cap(fsT))}), lateral ${fmt(lat.Hult / 1e3)} kN by ${lat.mode} (FS ${fmt(cap(fsH))}).`,
        ok ? (Number.isFinite(dL) && dL < -1 ? `The axial loads would be carried at ${fmt(Lreq)} m penetration; before shortening confirm lateral fixity (critical length), drivability and the cyclic capacity factor.` : 'Confirm the shaft friction with CPT-based methods and a drivability study; verify the tension piles for cyclic degradation.')
          : (w.name.startsWith('Axial') ? (Number.isFinite(Lreq) ? `Increase the penetration from ${fmt(inp.Lemb)} m to about ${fmt(Lreq)} m.` : `The axial loads cannot be carried within ${fmt(Lmax)} m: increase the diameter, add piles or widen the pile spacing to reduce push–pull.`) : `Increase the wall or diameter: the plastic moment needed is about ${fmt(P.Mp * w.uc / 1e3)} kN·m against ${fmt(P.Mp / 1e3)} kN·m.`),
        'API RP 2GEO (2011, add. 2014) §8; API RP 2A-WSD 22nd ed. §6.3–6.8; ISO 19901-4:2016',
        `Each metre of pile is ${fmt(tPerM)} t of steel ≈ ${fmt(tPerM * 1.9)} t CO₂e (screening factor 1.9 t CO₂e/t)${Number.isFinite(dL) ? `; the required penetration differs from the present one by ${fmt(dL)} m per pile.` : '.'}`),
      adv('economic', ok ? 'info' : 'warn', `Shaft friction supplies ${fmt(ax.Qfo / ax.Qc * 100)} % of the compression capacity; penetration costs ${fmt(tPerM)} t of steel per metre per pile plus driving time.`,
        'A site-specific CPT campaign and one instrumented pile (dynamic testing with signal matching) typically justify a higher capacity than the API defaults and shorten every pile: compare its cost with the steel and driving time saved.', `${ECON_BASIS}; value-of-information for site investigation`),
    ];
    return {
      kpis: [
        { key: 'Qc', label: 'Ultimate compression capacity', unit: 'kN', value: ax.Qc / 1e3, explain: 'Outer shaft friction plus end bearing (the smaller of plugged and coring unless forced). Divide by the factor of safety for the allowable load.' },
        { key: 'Qt', label: 'Ultimate tension capacity', unit: 'kN', value: ax.Qt / 1e3, explain: 'Outer shaft friction times the tension ratio; pile weight is not counted. Governs the rear piles of a dolphin under overturning.' },
        { key: 'Qs', label: 'Shaft friction (outer)', unit: 'kN', value: ax.Qfo / 1e3, explain: 'Integrated unit skin friction over the embedded outer surface: α·su in clay, β·σ′v (capped) in sand.' },
        { key: 'Qp', label: 'End bearing', unit: 'kN', value: ax.Qp / 1e3, explain: 'Tip resistance: 9·su or Nq·σ′v (capped) on the gross area if plugged, on the steel annulus plus internal friction if coring.' },
        { key: 'fsC', label: 'Factor of safety, compression', unit: '-', value: cap(fsC), limit: inp.fsC, sense: 'min', explain: 'Ultimate compression capacity over the design compression. Below the required value the pile is too short or too small.' },
        { key: 'fsT', label: 'Factor of safety, tension', unit: '-', value: cap(fsT), limit: inp.fsT, sense: 'min', explain: 'Ultimate tension capacity over the design tension.' },
        { key: 'Hult', label: 'Ultimate lateral capacity', unit: 'kN', value: lat.Hult / 1e3, explain: 'Head load at which the pile first reaches its plastic moment, or at which the soil can no longer balance the load. A long flexible dolphin pile is governed by the hinge.' },
        { key: 'fsH', label: 'Factor of safety, lateral', unit: '-', value: cap(fsH), limit: inp.fsH, sense: 'min', explain: 'Ultimate lateral capacity over the design pile shear.' },
        { key: 'Lreq', label: 'Penetration required for the axial loads', unit: 'm', value: Number.isFinite(Lreq) ? Lreq : Lmax, explain: 'Shallowest penetration at which both compression and tension meet their factors of safety (the largest depth examined if never met).' },
        { key: 'Mp', label: 'Plastic moment of the pile', unit: 'kN·m', value: P.Mp / 1e3, explain: 'Full plastic bending resistance of the corroded tube, fy·(D³ − Di³)/6. It caps the lateral capacity of a long pile.' },
        { key: 'plugged', label: 'Plugged at the tip (1 = yes)', unit: '-', value: ax.plugged ? 1 : 0, explain: 'Whether the soil plug carries the tip load (1) or the pile cores (0); it decides which end-bearing term applies.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Axial capacity against penetration', x: { label: 'Capacity', unit: 'kN' }, y: { label: 'Penetration', unit: 'm' },
          series: [{ name: 'Compression', x: cv.Qc.map((q) => round(q / 1e3)), y: cv.L.map((l) => -round(l)) }, { name: 'Tension', x: cv.Qt.map((q) => round(q / 1e3)), y: cv.L.map((l) => -round(l)) },
            { name: 'Required compression (FS × load)', x: [inp.fsC * inp.Nc, inp.fsC * inp.Nc], y: [-cv.L[0], -cv.L[cv.L.length - 1]], dash: true }, { name: 'Required tension (FS × load)', x: [inp.fsT * inp.Nt, inp.fsT * inp.Nt], y: [-cv.L[0], -cv.L[cv.L.length - 1]], dash: true }],
          explain: 'Ultimate capacity as the pile is driven deeper. Where each curve crosses its dashed line is the penetration required; steps mark layer boundaries and changes between plugged and coring behaviour.' },
        { type: 'line', title: 'Unit shaft friction with depth', x: { label: 'Unit friction', unit: 'kPa' }, y: { label: 'Depth', unit: 'm' }, series: [{ name: 'f', x: Array.from(ax.f, (f) => round(f / 1e3)), y: Array.from(ax.z, (z) => -round(z)) }],
          explain: 'Skin friction available at each depth. Weak upper layers contribute little; the plot shows which stratum the capacity depends on.' },
        { type: 'line', title: 'Lateral load against head deflection up to the limit', x: { label: 'Head deflection', unit: 'mm' }, y: { label: 'Head load', unit: 'kN' },
          series: [{ name: 'Pile', x: lat.y.map((v) => round(v * 1e3)), y: lat.H.map((v) => round(v / 1e3)) }, { name: 'Design shear', x: [0, round(lat.y[lat.y.length - 1] * 1e3)], y: [inp.Hd, inp.Hd], dash: true }],
          explain: 'Nonlinear load–deflection curve of the pile head to the point of first plastic hinge or soil failure. The curvature shows soil yielding; the end point is the lateral capacity.' },
      ],
      tables: [{ title: 'Capacity components', columns: [{ key: 'item', label: 'Component' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
        rows: [['Outer shaft friction', ax.Qfo / 1e3, 'kN'], ['Inner shaft friction', ax.Qfi / 1e3, 'kN'], ['End bearing, plugged', ax.QpPlugged / 1e3, 'kN'], ['End bearing, coring', ax.QpUnplugged / 1e3, 'kN'], ['Unit end bearing', ax.qTip / 1e3, 'kPa'],
          ['Lateral limit mode', lat.mode, ''], ...(broms ? [['Broms reference capacity (free head, uniform soil at 4 D)', broms.Hu / 1e3, 'kN'], ['Broms failure mode', broms.mode, '']] : [])].map(([item, v, u]) => ({ item, v: typeof v === 'number' ? round(v) : v, u })),
        explain: 'Breakdown of the axial capacity and the lateral limit, with the Broms closed-form value as an order-of-magnitude cross-check of the numerical lateral capacity.' }],
      notes: [`${w.name} governs (utilisation ${fmt(w.uc)}).`, `Lateral limit reached by ${lat.mode}.`, ax.plugged ? 'End bearing governed by the plugged condition.' : 'End bearing governed by coring (annulus plus internal friction).'],
      advice,
      publish: { 'dolphin.pileQc': ax.Qc, 'dolphin.pileQt': ax.Qt, 'dolphin.pileHult': lat.Hult, 'dolphin.LembRequired': Number.isFinite(Lreq) ? Lreq : Lmax },
      quantities: { steel_t: round(pileSteelT(inp, 1)) },
    };
  },
  verify: [
    // Uniform clay su = 20 kPa, γ′ = 10 kN/m³, D = 1 m, L = 20 m, α = 0.5ψ^−0.5 (ψ ≤ 1), 0.5ψ^−0.25 (ψ > 1), α ≤ 1, ψ = su/σ′v:
    //  0–2 m (ψ > 1): f = 10·(z/2)^0.25 kPa → ∫ = 10·2/1.25 = 16.00 kPa·m
    //  2–8 m: f = 0.5·√(20·10z) = 7.0711·√z → ∫ = 7.0711·(2/3)(8^1.5 − 2^1.5) = 93.33
    //  8–20 m (α = 1): 20·12 = 240     Σ = 349.33 kPa·m × πD = 1097.5 kN;  plugged tip: 9·20·π/4 = 141.37 kN.
    { name: 'API α method in uniform clay: shaft friction and plugged end bearing', source: 'Hand integration of the API RP 2GEO §8.1 α expressions (written out in the module source)',
      inputs: { D: 1, t: 25, corr: 0, Lemb: 20, Lfree: 10, scour: 0, plugMode: 'plugged', tensionFactor: 1, cyclicFactor: 1, layers: [{ type: 'clay', thick: 40, gamma: 10, su: 20, suGrad: 0, phi: 0, Dr: 0, eps50: 0.01 }] },
      expect: [{ kpi: 'Qs', value: 1097.5, tol: 0.01 }, { kpi: 'Qp', value: 141.37, tol: 1e-3 }, { kpi: 'Qt', value: 1097.5, tol: 0.01 }] },
    // Free-head cantilever built in 5 m below the mudline, 15 m above it: Hult = Mp/L with Mp = fy(D³ − Di³)/6
    // = 355e6·(1 − 0.96³)/6 = 6.820 MN·m → Hult = 6820/20 = 341.0 kN.
    { name: 'Plastic hinge of a cantilever pile: Hult = Mp/L', source: 'Plastic section modulus of a tube and statics of a cantilever',
      inputs: { D: 1, t: 20, corr: 0, steel: 's355', Lemb: 20, Lfree: 15, scour: 0, soilModel: 'fixed', fixDepth: 5, headFix: 'free', eLoad: 0, nElem: 80 },
      expect: [{ kpi: 'Hult', value: 341.0, tol: 0.005 }, { kpi: 'Mp', value: 6819.8, tol: 1e-3 }] },
  ],
};
