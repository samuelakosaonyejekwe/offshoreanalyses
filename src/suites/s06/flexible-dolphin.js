// Deflection and energy absorption of a flexible dolphin: the piles themselves absorb the berthing energy
// through bending, E = ∫ F dy on the nonlinear pile-head load–deflection curve.
import { solveLateralPile } from '../../engine/geotech/pile.js';
import { interp } from './_berthing.js';
import { DEG, need, fmt, round, adv, check, pileInputs, soilInputs, pileProps, profileFrom, latOpts, pileSteelT, steelCo2Text, CONFIRM, NOT_HIFI, ECON_BASIS } from './_common.js';

/** Head load–deflection curve of one pile up to the stress limit, with cumulative energy (trapezoidal ∫F dy). */
function headCurve(base, P, fixed, e, sigLim, Naxial) {
  const arm0 = (base.stickUp || 0) + e + 3 * P.D;
  const Mlim = Math.max(sigLim - Naxial / P.A, 0.05 * sigLim) * P.Z;
  const dH = (fixed ? 2 : 1) * Mlim / arm0 / 10;
  const c = { H: [0], y: [0], th: [0], M: [0], E: [0], yMud: [0] };
  let limited = 'stress limit';
  for (let k = 1; k <= 40; k++) {
    const h = k * dH;
    let s;
    try { s = solveLateralPile({ ...base, head: fixed ? { H: h, fixedRotation: true } : { H: h, M: h * e } }); } catch { limited = 'soil failure'; break; }
    const i = c.H.length - 1, m = Math.abs(s.maxMoment);
    // a head that does not move further in the load direction means the pile has lost lateral equilibrium
    if (!(s.head.y > c.y[i]) || !Number.isFinite(m)) { limited = 'soil failure'; break; }
    const im = s.z.findIndex((z) => z >= -1e-9);
    let f = 1, stop = false;
    if (m >= Mlim) { f = (Mlim - c.M[i]) / (m - c.M[i]); stop = true; }
    const Hn = c.H[i] + f * (h - c.H[i]), yn = c.y[i] + f * (s.head.y - c.y[i]);
    c.H.push(Hn); c.y.push(yn); c.th.push(c.th[i] + f * (s.head.theta - c.th[i])); c.M.push(c.M[i] + f * (m - c.M[i])); c.yMud.push(c.yMud[i] + f * (s.y[Math.max(im, 0)] - c.yMud[i]));
    c.E.push(c.E[i] + 0.5 * (Hn + c.H[i]) * (yn - c.y[i]));
    if (stop) break;
    if (k === 40) limited = 'load stepping limit';
  }
  need(c.H.length >= 3, 'The pile cannot sustain even a small lateral load with this soil and embedment: increase the embedded length.');
  return { ...c, limited, Mlim };
}

export default {
  id: 's06-flexible-dolphin', suite: 6, title: 'Deflection and energy absorption of a flexible dolphin',
  summary: 'Head deflection, rotation, stiffness and energy absorbed in bending by the piles of a flexible dolphin for a given berthing energy, with the stress and deflection limits.',
  covers: [6], fidelity: 'numerical',
  standards: ['BS 6349-2:2019 §7.9 (flexible dolphins)', 'EAU 2012 (Recommendations of the Committee for Waterfront Structures) R 69, R 218', 'PIANC WG33 (2002) §4', 'API RP 2GEO (2011) §8.5'],
  models: ['Elastic beam bending on nonlinear p–y springs', 'Strain energy as the area under the head load–deflection curve', 'Timoshenko shear deformation (optional)', 'Fender in series as an additional energy absorber', 'Energy-based impact assessment'],
  limits: 'Piles act in parallel with equal head deflection (rows tied by a cap or walings) and remain elastic: the curve stops at the allowable stress. Dynamic amplification, vessel hull flexibility and torsion from eccentric contact are not included. High-strength steel raises the energy capacity in proportion to fy² only if local buckling (D/t) permits. ' + NOT_HIFI,
  publishes: ['dolphin.energyCapacity', 'dolphin.fenderReaction', 'dolphin.capDisp', 'dolphin.Kh', 'dolphin.stressPerLoad'],
  searchTerms: 'flexible dolphin energy absorption tubular pile berthing deflection',
  inputs: [
    { key: 'nPiles', label: 'Piles acting together', unit: '-', type: 'int', default: 1, min: 1, max: 24, group: 'Dolphin', role: 'geometry', design: true },
    { key: 'headFix', label: 'Pile heads', type: 'select', default: 'free', group: 'Dolphin', role: 'bc', options: [{ value: 'free', label: 'Free (mono-pile, or piles linked by pinned walings)' }, { value: 'fixed', label: 'Rotation restrained by a rigid cap' }] },
    { key: 'eLoad', label: 'Contact height above the pile head', unit: 'm', type: 'number', default: 0, min: 0, max: 15, group: 'Dolphin', role: 'bc', showIf: { key: 'headFix', in: ['free'] }, help: 'Height of the fender contact above the top of the pile; zero when the vessel bears at pile-head level.' },
    ...pileInputs({ D: 2.4, t: 45, Lfree: 20, Lemb: 30 }),
    { key: 'timo', label: 'Include shear deformation (Timoshenko)', type: 'bool', default: false, group: 'Pile', role: 'control' },
    ...soilInputs(),
    { key: 'E', label: 'Berthing energy to absorb', unit: 'kJ', type: 'number', default: 515, min: 0.1, max: 1e5, group: 'Demand', role: 'bc', bind: 'dolphin.berthEnergyAbnormal', uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'Efender', label: 'Energy absorbed by a fender in series', unit: 'kJ', type: 'number', default: 0, min: 0, max: 1e5, group: 'Demand', role: 'bc', help: 'Zero for a bare flexible dolphin with rubbing strips. With a fender, enter the energy it absorbs at the reaction the piles reach (iterate with the fender module).' },
    { key: 'N', label: 'Axial load per pile', unit: 'kN', type: 'number', default: 300, min: 0, max: 1e5, group: 'Demand', role: 'ic', help: 'Self-weight of pile head, platform and equipment: the sustained state before impact.' },
    { key: 'y0', label: 'Residual head deflection before impact', unit: 'mm', type: 'number', default: 0, min: 0, max: 2000, group: 'Demand', role: 'ic', help: 'Permanent set observed on an existing dolphin; it is added to the computed deflection.' },
    { key: 'usage', label: 'Allowable stress / yield', unit: '-', type: 'number', default: 1.0, min: 0.5, max: 1, group: 'Acceptance', role: 'criteria', help: `1.0 of yield for abnormal berthing, about 0.67–0.8 for normal berthing in working-stress practice (EAU, BS 6349-2). ${CONFIRM}` },
    { key: 'yLimit', label: 'Allowable head deflection', unit: 'mm', type: 'number', default: 1500, min: 10, max: 5000, group: 'Acceptance', role: 'criteria', help: 'Set by vessel flare and clearance to the structure behind, loading-arm reach and walkway articulation.' },
    { key: 'nElem', label: 'Elements over the embedded length', unit: '-', type: 'int', default: 40, min: 10, max: 320, group: 'Numerical', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Pile elements (embedded length)', levels: [20, 40, 80], qoi: 'yHead' },
  compute(inp, ctx) {
    const P = pileProps(inp, ctx);
    const profile = profileFrom(inp, P.D);
    const fixed = inp.headFix === 'fixed', n = inp.nPiles;
    const base = latOpts(inp, P, profile, { axialLoad: inp.N * 1e3, timoshenko: !!inp.timo });
    const c = headCurve(base, P, fixed, fixed ? 0 : inp.eLoad, inp.usage * P.fy, inp.N * 1e3);
    // SI: kJ → J; energy per pile
    const Edem = Math.max(inp.E - inp.Efender, 0) * 1e3 / n;
    const Ecap1 = c.E[c.E.length - 1];
    const over = Edem > Ecap1;
    // deflection, load and moment vary nearly as √E: interpolate on √E so that coarse load steps stay accurate
    const rE = c.E.map(Math.sqrt);
    const at = (arr2) => (over ? arr2[arr2.length - 1] : interp(rE, arr2, Math.sqrt(Edem)));
    const y = at(c.y), H = at(c.H), M = at(c.M), th = at(c.th), yM = at(c.yMud);
    const sig = inp.N * 1e3 / P.A + M / P.Z;
    const k0 = c.H[1] / c.y[1], ksec = H / Math.max(y, 1e-9);
    const yTot = y * 1e3 + inp.y0;
    const checks = [
      check('Berthing energy against elastic energy capacity', Edem * n / 1e3, Ecap1 * n / 1e3, 'kJ', 'BS 6349-2:2019 §7.9; EAU 2012 R 218', 'Energy to be stored in pile bending against the energy at the allowable stress. If exceeded the pile yields: a permanent set remains and the capacity for the next berthing is reduced.'),
      check('Head deflection', yTot, inp.yLimit, 'mm', 'BS 6349-2:2019 §7.9 (project limit)', 'Movement of the dolphin head at the design energy against the clearance available. Too much deflection lets the hull flare or the loading arms reach their limits.'),
      check('Bending plus axial stress', sig / 1e6, inp.usage * P.fy / 1e6, 'MPa', 'EAU 2012 R 218; ISO 19902:2020 §13.3', 'Extreme-fibre stress at the section of maximum moment against the allowable fraction of yield.'),
    ];
    const w = checks.reduce((a, q) => (q.uc > a.uc ? q : a)), ok = w.uc <= 1;
    // elastic energy ∝ fy²·(section) : factor on fy to close an energy gap, or diameter by the same token
    const fyReq = P.fy * Math.sqrt(Math.max(Edem / Ecap1, 1e-9));
    const steelT = pileSteelT(inp, n);
    const eff = Ecap1 * n / 1e3 / steelT;
    const advice = [
      adv('engineering', ok ? 'ok' : 'fail', `At ${fmt(inp.E - inp.Efender)} kJ the head deflects ${fmt(yTot)} mm under ${fmt(H * n / 1e3)} kN, stress ${fmt(sig / 1e6)} MPa; elastic energy capacity ${fmt(Ecap1 * n / 1e3)} kJ at ${fmt(inp.usage * 100)} % of yield.`,
        ok ? `Energy reserve ${fmt((1 - Edem / Ecap1) * 100)} %. Check local buckling of the tube (D/t = ${fmt(P.D / P.t)}) at the section of maximum moment and repeat at lowest water level, where the contact point is lowest and the dolphin stiffest.`
          : (over ? `Raise the yield strength to about ${fmt(fyReq / 1e6)} MPa (energy capacity grows with fy² at unchanged stiffness), increase the wall, or add a fender in series absorbing at least ${fmt((Edem - Ecap1) * n / 1e3)} kJ.` : `Stiffen the dolphin (larger diameter, shorter free length) to bring the deflection within ${fmt(inp.yLimit)} mm; the present secant stiffness is ${fmt(ksec * n / 1e3)} kN/m.`),
        'BS 6349-2:2019 §7.9; EAU 2012 R 69 and R 218; PIANC WG33 (2002)',
        `This dolphin stores ${fmt(eff)} kJ per tonne of pile steel. High-strength steel raises that figure without extra mass: the piles hold ${steelCo2Text(steelT)}`),
      adv('economic', ok ? 'info' : 'warn', `${n} pile(s) of ${fmt(steelT / n)} t each absorb ${fmt(Ecap1 * n / 1e3)} kJ; the reaction passed to the hull at the design energy is ${fmt(H * n / 1e3)} kN.`,
        'Compare the flexible dolphin (steel only, no rubber unit to replace) with a rigid dolphin plus fender on steel tonnage, piling set-ups and fender replacement over the berth life; flexible dolphins win where soil is firm and water deep.', ECON_BASIS),
    ];
    return {
      kpis: [
        { key: 'yHead', label: 'Head deflection at the design energy', unit: 'mm', value: yTot, limit: inp.yLimit, sense: 'max', explain: 'Lateral movement of the pile head when the piles have absorbed the berthing energy (plus any residual set). The larger it is, the gentler the reaction on the hull.' },
        { key: 'H', label: 'Reaction on the vessel', unit: 'kN', value: H * n / 1e3, explain: 'Total force between hull and dolphin at maximum deflection. Check it against the hull pressure and the rubbing-strip area.' },
        { key: 'Ecap', label: 'Elastic energy capacity', unit: 'kJ', value: Ecap1 * n / 1e3, limit: Math.max(inp.E - inp.Efender, 0), sense: 'min', explain: 'Area under the head load–deflection curve up to the allowable stress, for all piles. It must cover the berthing energy not taken by a fender.' },
        { key: 'sigma', label: 'Stress at the design energy', unit: 'MPa', value: sig / 1e6, limit: inp.usage * P.fy / 1e6, sense: 'max', explain: 'Axial plus bending stress at the section of maximum moment. Reaching the allowable value marks the end of the usable elastic range.' },
        { key: 'k0', label: 'Initial lateral stiffness', unit: 'kN/m', value: k0 * n / 1e3, explain: 'Head load per unit deflection at small load. A soft dolphin absorbs energy at low reaction but moves more.' },
        { key: 'ksec', label: 'Secant stiffness at the design energy', unit: 'kN/m', value: ksec * n / 1e3, explain: 'Reaction over deflection at the design point; lower than the initial stiffness because the soil softens.' },
        { key: 'theta', label: 'Head rotation', unit: 'deg', value: th / DEG, explain: 'Tilt of the pile head at maximum deflection. It sets the contact angle on the hull and the travel needed in walkway bearings.' },
        { key: 'yMud', label: 'Deflection at the mudline', unit: 'mm', value: yM * 1e3, explain: 'Movement of the pile at the seabed; large values mean the soil near the surface is yielding and will degrade under repeated berthings.' },
        { key: 'M', label: 'Maximum bending moment per pile', unit: 'kN·m', value: M / 1e3, explain: 'Peak moment in the pile at the design energy, a few diameters below the mudline for a free head.' },
        { key: 'Hyield', label: 'Reaction at the allowable stress', unit: 'kN', value: c.H[c.H.length - 1] * n / 1e3, explain: 'Total head load at which the allowable stress is reached: the largest static lateral load (e.g. mooring or wind-on-vessel) the dolphin takes elastically.' },
        { key: 'eff', label: 'Energy stored per tonne of steel', unit: 'kJ/t', value: eff, explain: 'Elastic energy capacity over pile steel mass: the figure of merit for comparing diameters, walls and grades.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Head load against deflection', x: { label: 'Head deflection', unit: 'mm' }, y: { label: 'Total head load', unit: 'kN' },
          series: [{ name: 'Dolphin', x: c.y.map((v) => round(v * 1e3)), y: c.H.map((v) => round(v * n / 1e3)) }, { name: 'Design point', x: [round(y * 1e3)], y: [round(H * n / 1e3)] }],
          explain: 'Nonlinear load–deflection curve of the dolphin head up to the allowable stress. The area under it to the design point is the absorbed energy; the end of the curve is the elastic limit.' },
        { type: 'line', title: 'Absorbed energy against deflection', x: { label: 'Head deflection', unit: 'mm' }, y: { label: 'Energy', unit: 'kJ' },
          series: [{ name: 'Strain energy ∫F dy', x: c.y.map((v) => round(v * 1e3)), y: c.E.map((v) => round(v * n / 1e3)) }, { name: 'Energy to absorb', x: [0, round(c.y[c.y.length - 1] * 1e3)], y: [inp.E - inp.Efender, inp.E - inp.Efender], dash: true }],
          explain: 'Energy stored in the piles as the head deflects. The crossing with the dashed demand line gives the design deflection; no crossing means the piles yield first.' },
        { type: 'line', title: 'Maximum pile stress against head load', x: { label: 'Total head load', unit: 'kN' }, y: { label: 'Stress', unit: 'MPa' },
          series: [{ name: 'Axial + bending', x: c.H.map((v) => round(v * n / 1e3)), y: c.M.map((m) => round((inp.N * 1e3 / P.A + m / P.Z) / 1e6)) }, { name: 'Allowable', x: [0, round(c.H[c.H.length - 1] * n / 1e3)], y: [inp.usage * P.fy / 1e6, inp.usage * P.fy / 1e6], dash: true }],
          explain: 'Stress at the critical section as the load grows. The slight upward curvature comes from the moment peak moving deeper as the soil yields.' },
      ],
      tables: [{ title: 'Load–deflection–energy curve', columns: [{ key: 'H', label: 'Head load', unit: 'kN' }, { key: 'y', label: 'Head deflection', unit: 'mm' }, { key: 'E', label: 'Energy', unit: 'kJ' }, { key: 'M', label: 'Max. moment per pile', unit: 'kN·m' }, { key: 's', label: 'Stress', unit: 'MPa' }],
        rows: c.H.map((h, i) => ({ H: round(h * n / 1e3), y: round(c.y[i] * 1e3), E: round(c.E[i] * n / 1e3), M: round(c.M[i] / 1e3), s: round((inp.N * 1e3 / P.A + c.M[i] / P.Z) / 1e6) })),
        explain: 'Each row is one load step of the nonlinear pile analysis: use it as the force–deflection characteristic of the dolphin in a berthing or mooring model.' }],
      notes: [`${w.name} governs (utilisation ${fmt(w.uc)}).`, `Curve ended by ${c.limited}.`, ...(over ? ['The demand exceeds the elastic capacity: values are reported at the allowable stress.'] : [])],
      advice,
      publish: { 'dolphin.energyCapacity': Ecap1 * n + inp.Efender * 1e3, 'dolphin.fenderReaction': H * n, 'dolphin.capDisp': yTot / 1e3, 'dolphin.Kh': ksec * n, 'dolphin.stressPerLoad': H > 0 ? (M / P.Z) / (H * n) : 0 },
      quantities: { steel_t: round(steelT), vessel_days: { piling: round(n * 0.75 + 1.5, 3) } },
    };
  },
  verify: [
    // Cantilever built in 5 m below the mudline, 15 m free: Lc = 20 m, EI = 210e9·π/64·(1 − 0.96⁴) = 1.5530e9 N·m².
    // Strain energy E = H²Lc³/(6EI) → H = √(6·EI·E/Lc³) = √(6·1.5530e9·1e5/8000) = 341.3 kN; y = H·Lc³/(3EI) = 586.0 mm;
    // stress M/Z = H·Lc/(2I/D) = 6.826e6/0.014790 = 461.5 MPa (fy 690 MPa so the pile stays elastic).
    { name: 'Cantilever pile: E = H²L³/(6EI), y = HL³/(3EI)', source: 'Elementary beam theory (strain energy of a tip-loaded cantilever)',
      inputs: { nPiles: 1, headFix: 'free', eLoad: 0, D: 1, t: 20, corr: 0, steel: 's690', Lfree: 15, Lemb: 20, scour: 0, soilModel: 'fixed', fixDepth: 5, E: 100, Efender: 0, N: 0, y0: 0, usage: 1, yLimit: 5000, nElem: 160, timo: false },
      expect: [{ kpi: 'yHead', value: 586.0, tol: 0.01 }, { kpi: 'H', value: 341.3, tol: 0.01 }, { kpi: 'sigma', value: 461.5, tol: 0.01 }] },
  ],
};
