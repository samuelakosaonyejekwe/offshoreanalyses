// Sliding and overturning stability of a mudmat on a sloping seabed, with skirts and torsion.
import { slidingResistance, overturningCheck } from '../../engine/geotech/bearing.js';
import { torsionEquivalentH, torsionalCapacityRect } from '../../engine/geotech/mudmat.js';
import { rankine } from '../../engine/geotech/earthpressure.js';
import { DEG, need, run, fmt, round, cap, adv, fsCheck, worst, steelCo2Text, matInputs, soilInputs, soilFrom, CONFIRM } from './_common.js';

function evaluate(inp, over = {}) {
  const p = { ...inp, ...over };
  const s = soilFrom(p);
  const bx = p.slopeX * DEG, by = p.slopeY * DEG;
  const W = p.W * 1e3;
  // Weight resolved on the inclined base: normal component carries friction, the down-slope components add to H.
  const N = W * Math.cos(bx) * Math.cos(by) + p.Vint * 1e3;
  need(N > 0, 'The net normal load on the base is not compressive: the structure lifts off. Reduce the uplift reaction or add weight.');
  const Hx = p.Hx * 1e3 + W * Math.sin(bx), Hy = p.Hy * 1e3 + W * Math.sin(by);
  const Hres = Math.hypot(Hx, Hy);
  const A = p.B * p.L, T = p.T * 1e3;
  const Heq = T > 0 ? torsionEquivalentH({ H: Hres, T, Leff: Math.max(p.B, p.L) }) : Hres;
  // Skirt face normal to the resultant: projected width of the rectangle seen along the load direction.
  const dir = Hres > 0 ? Math.atan2(Math.abs(Hy), Math.abs(Hx)) : 0;
  const wSkirt = p.L * Math.cos(dir) + p.B * Math.sin(dir);
  const sr = run(() => (s.clay
    ? slidingResistance({ mode: 'undrained', A, suBase: p.alphaInt * s.suAt(p.D), skirtDepth: p.D, skirtWidth: wSkirt, suSkirt: s.suAt(p.D / 2), gamma: s.gamma, skirtMobilisation: p.mob })
    : slidingResistance({ mode: 'drained', A, V: N, delta: p.deltaDeg * DEG, c: s.c, skirtDepth: p.D, skirtWidth: wSkirt, phi: s.phi, gamma: s.gamma, skirtMobilisation: p.mob })), 'Sliding resistance');
  const fsS = Heq > 0 ? sr.H / Heq : Infinity;
  // Overturning about each base edge: horizontal loads at their lever arm, weight component at the centre of gravity.
  const Mx = p.Hx * 1e3 * p.zH + W * Math.sin(bx) * p.zG + p.Mx * 1e3;  // overturning in the x–z plane (eccentricity along B)
  const My = p.Hy * 1e3 * p.zH + W * Math.sin(by) * p.zG + p.My * 1e3;  // eccentricity along L
  const ox = overturningCheck({ V: N, M: Math.abs(Mx), B: p.B });
  const oy = overturningCheck({ V: N, M: Math.abs(My), B: p.L });
  const Tult = s.clay ? torsionalCapacityRect({ B: Math.min(p.B, p.L), L: Math.max(p.B, p.L), su: p.alphaInt * s.suAt(p.D) }) : NaN;
  return { s, N, W, Hx, Hy, Hres, Heq, sr, fsS, Mx, My, ox, oy, A, T, Tult, wSkirt };
}

/** Best-fit plane z = a + b·x + c·y through the grid soundings inside the mat footprint centred on the grid (whole grid if too few). */
function planeFit(g, B, L) {
  need(g.kind === 'grid' && g.nx >= 2 && g.ny >= 2 && g.z && g.z.length >= g.nx * g.ny, 'The seabed surface must be a bathymetry grid with at least 2 × 2 soundings.');
  const cx = g.x0 + g.dx * (g.nx - 1) / 2, cy = g.y0 + g.dy * (g.ny - 1) / 2;
  const collect = (all) => { const pts = []; for (let j = 0; j < g.ny; j++) for (let i = 0; i < g.nx; i++) { const x = g.x0 + i * g.dx - cx, y = g.y0 + j * g.dy - cy, z = g.z[j * g.nx + i]; if (Number.isFinite(z) && (all || (Math.abs(x) <= B / 2 && Math.abs(y) <= L / 2))) pts.push([x, y, z]); } return pts; };
  let pts = collect(false); if (pts.length < 4) pts = collect(true);
  need(pts.length >= 3, 'The bathymetry grid holds fewer than three valid soundings.');
  let sx = 0, sy = 0, sz = 0; for (const p of pts) { sx += p[0]; sy += p[1]; sz += p[2]; }
  const n = pts.length, mx = sx / n, my = sy / n, mz = sz / n;
  let sxx = 0, sxy = 0, syy = 0, sxz = 0, syz = 0;
  for (const p of pts) { const x = p[0] - mx, y = p[1] - my, z = p[2] - mz; sxx += x * x; sxy += x * y; syy += y * y; sxz += x * z; syz += y * z; }
  const det = sxx * syy - sxy * sxy;
  need(det > 1e-12, 'The soundings inside the footprint are collinear: a seabed plane cannot be fitted.');
  const b = (sxz * syy - syz * sxy) / det, c = (syz * sxx - sxz * sxy) / det;
  let ss = 0; for (const p of pts) { const e = p[2] - mz - b * (p[0] - mx) - c * (p[1] - my); ss += e * e; }
  return { slopeX: Math.abs(Math.atan(b)) / DEG, slopeY: Math.abs(Math.atan(c)) / DEG, rms: Math.sqrt(ss / n), n };
}

export default {
  id: 's04-sliding-overturning', suite: 4, title: 'Sliding and overturning stability',
  summary: 'Factors of safety against sliding and overturning of a skirted or flat mudmat on a sloping seabed under horizontal load, moment and torsion.',
  covers: [4], fidelity: 'analytical',
  standards: ['ISO 19901-4:2016 §7.4', 'API RP 2GEO (2011) §7.3.3', 'DNV-RP-C212 (2019) §5.4', 'API RP 2A-WSD (21st ed.) §6.13–6.14'],
  models: ['Rigid-body sliding equilibrium', 'Overturning moment equilibrium', 'Coulomb friction', 'Interface adhesion', 'Undrained shear sliding resistance', 'Rankine passive and active skirt resistance', 'Torsion-equivalent horizontal load'],
  limits: 'Rigid mat, planar seabed, uniform interface strength; passive skirt resistance from Rankine limit states scaled by a mobilisation factor. Torsion is combined through the DNV equivalent-load expression (recalled, to be confirmed). Deep-seated failure below the skirts, cyclic degradation and scour are not included here.',
  publishes: ['subsea.fsSliding', 'subsea.fsOverturning'],
  searchTerms: 'mudmat sliding overturning skirt passive resistance seabed slope',
  inputs: [
    ...matInputs(),
    { key: 'D', label: 'Skirt depth', unit: 'm', type: 'number', default: 0.5, min: 0, max: 10, group: 'Mudmat geometry', role: 'geometry', design: true },
    { key: 'zG', label: 'Centre of gravity above the base', unit: 'm', type: 'number', default: 1.6, min: 0, max: 30, group: 'Mudmat geometry', role: 'geometry', bind: 'subsea.cogZ' },
    { key: 'bathy', label: 'Surveyed seabed surface (optional)', type: 'geometry', kinds: ['grid'], default: null, group: 'Seabed', role: 'ic', help: 'Bathymetry grid (GeoTIFF, ESRI ASCII, XYZ). When given, the slopes along B and L are the best-fit plane over the mat footprint centred on the grid, and replace the two values below.' },
    { key: 'slopeX', label: 'Seabed slope along B', unit: 'deg', type: 'number', default: 2, min: 0, max: 25, group: 'Seabed', role: 'ic', uncertain: { dist: 'normal', cov: 0.3 }, help: 'Longitudinal slope; the down-slope weight component is added to the load along B.' },
    { key: 'slopeY', label: 'Seabed slope along L', unit: 'deg', type: 'number', default: 0, min: 0, max: 25, group: 'Seabed', role: 'ic' },
    ...soilInputs(),
    { key: 'alphaInt', label: 'Base interface factor', unit: '-', type: 'number', default: 1.0, min: 0.05, max: 1, group: 'Soil', role: 'soil', uncertain: { dist: 'normal', cov: 0.15 }, showIf: { key: 'soil', in: ['clay'] },
      help: '1.0 when skirts force shearing in intact soil at tip level; about 1/sensitivity (0.2–0.5) for a flat steel plate on clay.' },
    { key: 'deltaDeg', label: 'Interface friction angle δ', unit: 'deg', type: 'number', default: 27, min: 5, max: 45, group: 'Soil', role: 'soil', uncertain: { dist: 'normal', cov: 0.1 }, showIf: { key: 'soil', in: ['sand'] }, help: 'φ′ − 5° for steel on sand; φ′ for a skirted base.' },
    { key: 'mob', label: 'Skirt resistance mobilised', unit: '-', type: 'number', default: 0.5, min: 0, max: 1, group: 'Soil', role: 'soil', help: 'Passive resistance needs more displacement than base shear; 0.5 is a common allowance, 0 ignores the skirts.' },
    { key: 'W', label: 'Submerged weight', unit: 'kN', type: 'number', default: 1000, min: 1, max: 2e5, group: 'Loads', role: 'ic', bind: 'subsea.submergedWeight', uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'Vint', label: 'Vertical interface reaction (down +)', unit: 'kN', type: 'number', default: 0, min: -1e5, max: 1e5, group: 'Loads', role: 'bc', help: 'Net vertical pipeline / jumper reaction and hydrodynamic lift (negative).' },
    { key: 'Hx', label: 'Horizontal load along B', unit: 'kN', type: 'number', default: 120, min: 0, max: 1e5, group: 'Loads', role: 'bc', bind: 'subsea.designH', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'Hy', label: 'Horizontal load along L', unit: 'kN', type: 'number', default: 30, min: 0, max: 1e5, group: 'Loads', role: 'bc' },
    { key: 'zH', label: 'Lever arm of the horizontal loads', unit: 'm', type: 'number', default: 2.0, min: 0, max: 40, group: 'Loads', role: 'bc', bind: 'subsea.hydroArm' },
    { key: 'Mx', label: 'Additional moment (eccentricity along B)', unit: 'kN·m', type: 'number', default: 150, min: 0, max: 1e6, group: 'Loads', role: 'bc', help: 'Moments not already contained in load × lever arm, e.g. connector moments.' },
    { key: 'My', label: 'Additional moment (eccentricity along L)', unit: 'kN·m', type: 'number', default: 0, min: 0, max: 1e6, group: 'Loads', role: 'bc' },
    { key: 'T', label: 'Torsion about the vertical axis', unit: 'kN·m', type: 'number', default: 200, min: 0, max: 1e6, group: 'Loads', role: 'bc', bind: 'subsea.designT' },
    { key: 'fsSlide', label: 'Required sliding factor of safety', unit: '-', type: 'number', default: 1.5, min: 1, max: 4, group: 'Acceptance criteria', role: 'criteria', help: `API RP 2A-WSD 6.13: 1.5 for sliding. ${CONFIRM}` },
    { key: 'fsOT', label: 'Required overturning factor of safety', unit: '-', type: 'number', default: 1.5, min: 1, max: 4, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
  ],
  compute(inp0) {
    // SI: kN → N, kN·m → N·m, kPa → Pa, deg → rad (inside evaluate).
    const fit = inp0.bathy ? planeFit(inp0.bathy, inp0.B, inp0.L) : null;
    const inp = fit ? { ...inp0, slopeX: fit.slopeX, slopeY: fit.slopeY } : inp0;
    need(inp.slopeX <= 25 && inp.slopeY <= 25, `The surveyed seabed slope (${fmt(inp.slopeX)}° / ${fmt(inp.slopeY)}°) exceeds 25°: outside the range of this rigid-block check.`);
    const r = evaluate(inp);
    const fsOx = r.ox.fsOverturning, fsOy = r.oy.fsOverturning;
    const checks = [
      fsCheck('Sliding (with torsion)', r.fsS, inp.fsSlide, 'ISO 19901-4:2016 §7.4 / DNV-RP-C212 §5.4', 'Guards against translation of the mat along the seabed. If exceeded the structure moves until skirts or connected pipework arrest it, overloading connectors.'),
      fsCheck('Overturning about the edge along L', fsOx, inp.fsOT, 'Statics: N·B/2 over the overturning moment', 'Guards against rotation about the long edge. If exceeded the heel lifts and bearing pressure concentrates at the toe.'),
      fsCheck('Overturning about the edge along B', fsOy, inp.fsOT, 'Statics: N·L/2 over the overturning moment', 'Guards against rotation about the short edge under transverse load.'),
    ];
    if (r.T > 0 && Number.isFinite(r.Tult)) checks.push({ name: 'Pure torsion on the base', demand: round(r.T / 1e3), capacity: round(r.Tult / 1e3), unit: 'kN·m', uc: Math.min(r.T / r.Tult, 99), ref: 'Plastic torsional capacity of a rectangular base in full contact', explain: 'Torque from eccentric pulls against the plastic twisting resistance of the base. If exceeded the mat rotates in plan.' });
    const w = worst(checks), ok = w.uc <= 1;
    // Sweeps: slope angle and horizontal load.
    const sl = [], fsSl = [], fsOl = [];
    for (let i = 0; i <= 15; i++) { try { const q = evaluate(inp, { slopeX: i }); sl.push(i); fsSl.push(cap(q.fsS)); fsOl.push(cap(q.ox.fsOverturning)); } catch { /* lift-off */ } }
    const hx = [], fsH = [];
    for (let i = 1; i <= 16; i++) { const h = inp.Hx * (0.25 + i * 0.125) || i * 10; try { hx.push(round(h)); fsH.push(cap(evaluate(inp, { Hx: h }).fsS)); } catch { hx.pop(); } }
    // Skirt depth that restores the sliding factor (bisection on D).
    let Dreq = inp.D;
    if (r.fsS < inp.fsSlide && inp.mob > 0) {
      let lo = inp.D, hi = 10;
      if (evaluate(inp, { D: hi }).fsS >= inp.fsSlide) { for (let i = 0; i < 40; i++) { const mid = 0.5 * (lo + hi); if (evaluate(inp, { D: mid }).fsS >= inp.fsSlide) hi = mid; else lo = mid; } Dreq = hi; } else Dreq = 99;
    }
    const perim = 2 * (inp.B + inp.L);
    const skirtSteel = (d) => perim * d * 0.02 * 7.85; // 20 mm skirt plate, t
    const Ka = !r.s.clay ? rankine(r.s.phi) : null;
    const advice = [];
    if (ok) advice.push(adv('engineering', 'ok', `Sliding factor ${fmt(cap(r.fsS))} and overturning factors ${fmt(cap(fsOx))} / ${fmt(cap(fsOy))} meet the criteria; ${fmt(r.sr.Hbase / (r.sr.H || 1) * 100)} % of the sliding resistance comes from base shear.`,
      'Re-run with the lower-bound interface strength and with the as-surveyed seabed slope; confirm that skirt resistance is not double-counted with scour protection assumptions.',
      'ISO 19901-4:2016 §7.4; API RP 2GEO (2011) §7.3.3; DNV-RP-C212 (2019) §5.4',
      `The reserve allows the skirts to stay at ${fmt(inp.D)} m: each additional 0.1 m of 20 mm skirt on this perimeter is ${steelCo2Text(skirtSteel(0.1))}`));
    else advice.push(adv('engineering', 'fail', `${w.name}: achieved factor ${fmt(w.capacity)} against ${fmt(w.demand)} required.`,
      w.name.startsWith('Sliding') ? (Dreq < 99 && Dreq > inp.D ? `Deepen the skirts from ${fmt(inp.D)} m to about ${fmt(Dreq)} m, or reduce the horizontal resultant ${fmt(r.Heq / 1e3)} kN (torsion included).` : `Skirts alone do not close the gap: enlarge the base area (resistance ∝ area on clay) or add weight on sand; required resistance is ${fmt(inp.fsSlide * r.Heq / 1e3)} kN against ${fmt(r.sr.H / 1e3)} kN.`)
        : `Widen the base in the direction of overturning or lower the load line: the required restoring moment is ${fmt(inp.fsOT * Math.max(Math.abs(r.Mx), Math.abs(r.My)) / 1e3)} kN·m.`,
      'ISO 19901-4:2016 §7.4; API RP 2GEO (2011) §7.3.3; DNV-RP-C212 (2019) §5.4',
      w.name.startsWith('Sliding') && Dreq < 99 && Dreq > inp.D ? `Deeper skirts add ${steelCo2Text(skirtSteel(Dreq - inp.D))}` : 'A larger base increases both steel mass (about 0.16 t per m² of mat) and the disturbed seabed footprint; compare against seabed levelling to remove the slope component.'));
    advice.push(adv('economic', ok ? 'info' : 'warn', `Skirt plate on this ${fmt(perim)} m perimeter weighs about ${fmt(skirtSteel(1))} t per metre of depth (20 mm plate); the seabed slope alone contributes ${fmt(r.W * Math.sin(inp.slopeX * DEG) / 1e3)} kN of driving force.`,
      ok ? 'Keep the skirt depth; a pre-installation survey confirming the slope avoids carrying a conservative slope allowance into fabrication.' : 'Compare deeper skirts (fabrication tonnage, higher penetration resistance) with seabed preparation (about one vessel day of survey and levelling) before committing.',
      'AACE RP 18R-97 estimate classes; whole-life costing per ISO 15686-5'));
    return {
      kpis: [
        { key: 'fsSlide', label: 'Sliding factor of safety', unit: '-', value: cap(r.fsS), limit: inp.fsSlide, sense: 'min', explain: 'Base shear plus mobilised skirt resistance over the horizontal resultant (slope component and torsion-equivalent load included). Below the required value the mat may translate.' },
        { key: 'fsOTx', label: 'Overturning factor (about the L edge)', unit: '-', value: cap(fsOx), limit: inp.fsOT, sense: 'min', explain: 'Restoring moment N·B/2 over the overturning moment in the B direction. Below 1 the mat tips; below about 3 the heel lifts off an unskirted mat.' },
        { key: 'fsOTy', label: 'Overturning factor (about the B edge)', unit: '-', value: cap(fsOy), limit: inp.fsOT, sense: 'min', explain: 'Restoring moment N·L/2 over the overturning moment in the L direction.' },
        { key: 'Hres', label: 'Sliding resistance', unit: 'kN', value: r.sr.H / 1e3, explain: 'Total horizontal resistance: base shear (adhesion or friction) plus the mobilised net passive resistance of the skirts.' },
        { key: 'Hdrive', label: 'Driving horizontal force', unit: 'kN', value: r.Heq / 1e3, explain: 'Resultant of applied horizontal loads, the down-slope weight component and the torsion-equivalent load.' },
        { key: 'Mrestore', label: 'Restoring moment (B direction)', unit: 'kN·m', value: r.N * inp.B / 2 / 1e3, explain: 'Normal load times half the base width: the moment available to resist overturning about the long edge.' },
        { key: 'Mover', label: 'Overturning moment (B direction)', unit: 'kN·m', value: Math.abs(r.Mx) / 1e3, explain: 'Horizontal load times its lever arm, the slope component of weight at the centre of gravity and applied moments.' },
        { key: 'ecc', label: 'Eccentricity along B', unit: 'm', value: r.ox.e, limit: inp.B / 6, sense: 'max', explain: 'Offset of the base reaction from the centre. Inside B/6 the whole base stays in compression.' },
        { key: 'margin', label: 'Overturning stability margin', unit: 'kN·m', value: (r.N * inp.B / 2 - Math.abs(r.Mx)) / 1e3, explain: 'Restoring minus overturning moment in the B direction; negative means the mat rotates about its edge.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Sliding resistance and driving force', categories: ['Base shear', 'Skirt passive (mobilised)', 'Skirt active (mobilised)', 'Driving force'],
          series: [{ name: 'Force (kN)', y: [r.sr.Hbase / 1e3, inp.mob * r.sr.Hpassive / 1e3, -inp.mob * r.sr.Hactive / 1e3, r.Heq / 1e3].map((v) => round(v)) }], y: { label: 'Force', unit: 'kN' },
          explain: 'Where the horizontal resistance comes from, next to the load it has to carry. A small base-shear bar means the result depends on the assumed skirt mobilisation.' },
        { type: 'line', title: 'Factors of safety against seabed slope', x: { label: 'Slope along B', unit: 'deg' }, y: { label: 'Factor of safety', unit: '-' },
          series: [{ name: 'Sliding', x: sl, y: fsSl }, { name: 'Overturning', x: sl, y: fsOl }, { name: 'Required sliding', x: [0, 15], y: [inp.fsSlide, inp.fsSlide], dash: true }],
          explain: 'How both factors fall as the seabed steepens. The crossing with the dashed line is the installation tolerance on slope.' },
        { type: 'line', title: 'Sliding factor against horizontal load', x: { label: 'Horizontal load along B', unit: 'kN' }, y: { label: 'Sliding factor of safety', unit: '-' },
          series: [{ name: 'Sliding', x: hx, y: fsH }, { name: 'Required', x: [hx[0], hx[hx.length - 1]], y: [inp.fsSlide, inp.fsSlide], dash: true }],
          explain: 'Sliding reserve for a range of applied load. Read off the largest horizontal load the foundation tolerates at the required factor.' },
      ],
      tables: [{ title: 'Force and moment balance', columns: [{ key: 'item', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
        rows: [['Normal load on the base', r.N / 1e3, 'kN'], ['Down-slope weight component (B)', r.W * Math.sin(inp.slopeX * DEG) / 1e3, 'kN'], ['Horizontal resultant', r.Hres / 1e3, 'kN'], ['Torsion-equivalent horizontal load', r.Heq / 1e3, 'kN'],
          ['Base shear resistance', r.sr.Hbase / 1e3, 'kN'], ['Skirt passive resistance (full)', r.sr.Hpassive / 1e3, 'kN'], ['Skirt active thrust (full)', r.sr.Hactive / 1e3, 'kN'], ['Skirt width normal to load', r.wSkirt, 'm'],
          ['Kern limit B/6', r.ox.kernLimit, 'm'], ['Contact width B − 2e', r.ox.contactWidth, 'm'], ...(Ka ? [['Rankine Ka', Ka.Ka, '-'], ['Rankine Kp', Ka.Kp, '-']] : []), ...(Number.isFinite(r.Tult) ? [['Torsional capacity of the base', r.Tult / 1e3, 'kN·m']] : [])]
          .map(([item, v, u]) => ({ item, v: round(v), u })),
        explain: 'Each row is one term of the horizontal force balance or the moment balance about the base edge, for checking against a hand calculation.' }],
      notes: [...(fit ? [`Slopes from the surveyed grid over the footprint: ${fmt(fit.slopeX)}° along B, ${fmt(fit.slopeY)}° along L; seabed unevenness (rms deviation from the plane) ${fmt(fit.rms * 1e3)} mm from ${fit.n} soundings.`] : []), `${w.name} governs (required/achieved = ${fmt(w.uc)}).`, r.ox.fullContact ? 'The base reaction stays inside the kern: full contact along B.' : 'The base reaction is outside the kern along B: part of an unskirted base lifts off.'],
      advice,
      publish: { 'subsea.fsSliding': cap(r.fsS), 'subsea.fsOverturning': cap(Math.min(fsOx, fsOy)) },
      quantities: {},
    };
  },
  verify: [
    { name: 'Drained sliding on a level bed: FS = V·tan δ / H', source: 'Hand calculation: 1000·tan 30° = 577.35 kN; FS = 577.35/200 = 2.887',
      inputs: { soil: 'sand', deltaDeg: 30, W: 1000, Vint: 0, Hx: 200, Hy: 0, T: 0, D: 0, slopeX: 0, slopeY: 0, Mx: 0, zH: 5, B: 10 },
      expect: [{ kpi: 'fsSlide', value: 2.88675, tol: 1e-4 }, { kpi: 'fsOTx', value: 5, tol: 1e-6 }] },
    { name: 'Undrained sliding on a slope: FS = su·A / (W sin β)', source: 'Hand calculation: su·A = 5·140 = 700 kN; W sin 10° = 173.65 kN; FS = 4.031',
      inputs: { soil: 'clay', su0: 5, suGrad: 0, alphaInt: 1, W: 1000, Hx: 0, Hy: 0, T: 0, D: 0, slopeX: 10, Mx: 0, B: 10, L: 14 },
      expect: [{ kpi: 'fsSlide', value: 4.0311, tol: 1e-4 }] },
  ],
};
