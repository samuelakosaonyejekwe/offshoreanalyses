// Wake interaction and array performance: farm energy over the wind rose with an engineering wake model, per-turbine
// and per-sector array efficiency, wake-added turbulence for fatigue, layout checks and the inter-array cable length.
import { farmAep } from '../../engine/aero/aep.js';
import { farmFlow, jensenDeficit, gaussianDeficit, wakeGrowthFromTI, effectiveTurbulence, deepArrayFrandsen, superpose, rotorPoints } from '../../engine/aero/wake.js';
import { cableLayout, capacityDensity, checkLayout } from '../../engine/aero/layout.js';
import { curveInputs, getCurve } from './_curve.js';
import { costInput, money, site, adv, check, worst, sev, fmt, round, run, need, linspace, interp, DEG, LIMIT_WAKE, LIMIT_POLAR, LIMIT_CODE, CONFIRM, energyCo2Text } from './_common.js';

// Indicative 12-sector rose with south-westerly prevailing winds (frequency %, sector Weibull scale / omnidirectional scale).
const ROSE = [[0, 5, 0.85], [30, 5, 0.85], [60, 6, 0.88], [90, 6, 0.9], [120, 6, 0.92], [150, 7, 0.95], [180, 9, 1.02], [210, 13, 1.1], [240, 15, 1.14], [270, 12, 1.1], [300, 9, 1.02], [330, 7, 0.92]];

function layoutPoints(inp, D) {
  const g = inp.layout;
  if (g && g.positions && g.positions.length >= 3) {
    const n = Math.floor(g.positions.length / 3), pts = [];
    for (let i = 0; i < n; i++) { const x = g.positions[3 * i], y = g.positions[3 * i + 1]; need(Number.isFinite(x) && Number.isFinite(y), 'The imported layout contains non-numeric coordinates.'); pts.push([x, y]); }
    need(n <= 100, `The imported layout has ${n} positions; the farm solver in this module is limited to 100 turbines per run (run time grows with the square of the number: split the site or decimate).`);
    return { pts, imported: true };
  }
  const a = inp.rowAngle * DEG, ex = [Math.cos(a), Math.sin(a)], ey = [-Math.sin(a), Math.cos(a)], pts = [];
  need(inp.nRows * inp.nCols <= 100, `The parametric grid is limited to 100 turbines per run (${inp.nRows} × ${inp.nCols} entered); run time grows with the square of the number.`);
  for (let r = 0; r < inp.nRows; r++) for (let c = 0; c < inp.nCols; c++) {
    const s = (c + (inp.stagger && r % 2 ? 0.5 : 0)) * inp.sx * D, t = r * inp.sy * D;
    pts.push([s * ex[0] + t * ey[0], s * ex[1] + t * ey[1]]);
  }
  return { pts, imported: false };
}

export default {
  id: 's16-wake-array', suite: 16, title: 'Wake interaction and array performance',
  summary: 'Wake losses, array efficiency and array-level energy of a wind farm over the wind rose, with per-turbine results, wake-added turbulence and layout checks.',
  covers: [14], fidelity: 'reduced-order',
  standards: ['IEC 61400-1:2019 Annex E (wake and wind-farm turbulence, Frandsen effective turbulence)', 'IEC 61400-15-1 (site energy yield assessment)', 'Jensen (1983) Risø-M-2411 / Katić, Højstrup & Jensen (1986)', 'Bastankhah & Porté-Agel (2014) Renewable Energy 70', 'Niayifar & Porté-Agel (2016) Energies 9'],
  models: ['Jensen / Park top-hat wake model', 'Bastankhah–Porté-Agel Gaussian wake deficit', 'Wake growth from local turbulence intensity', 'Frandsen and Crespo–Hernández wake-added turbulence', 'Wake superposition (root-sum-square, linear, energy, maximum)', 'Summed-induction turbine blockage', 'Frandsen deep-array (infinite farm) roughness model', 'Weibull wind rose integration', 'Esau–Williams capacitated cable strings'],
  limits: `${LIMIT_WAKE} Blockage is a summed-induction estimate, not a flow solution; turbines share one hub height and one curve. ${LIMIT_POLAR} Layout optimisation is left to the generic design studies (spacing and orientation are design variables). ${LIMIT_CODE}`,
  publishes: ['energy.wakeLossFrac', 'energy.aepArray', 'array.cableLength', 'array.count', 'array.effectiveTI'],
  searchTerms: 'wind farm wake model Jensen Gaussian Bastankhah array efficiency layout',
  inputs: [
    ...curveInputs(),
    { key: 'hubHeight', label: 'Hub height', unit: 'm', type: 'number', default: 90, min: 20, max: 220, group: 'Generating unit', role: 'geometry', bind: 'owt.hubHeight' },
    { key: 'layout', label: 'Imported turbine coordinates (optional)', type: 'geometry', kinds: ['points'], default: null, group: 'Array arrangement', role: 'geometry',
      help: 'CSV / XYZ / GeoJSON point set in metres (x east, y north), one point per turbine. Leave empty to use the parametric grid below.' },
    { key: 'nRows', label: 'Rows', type: 'int', default: 4, min: 1, max: 30, group: 'Array arrangement', role: 'geometry', design: true },
    { key: 'nCols', label: 'Turbines per row', type: 'int', default: 5, min: 1, max: 30, group: 'Array arrangement', role: 'geometry', design: true },
    { key: 'sx', label: 'In-row spacing', unit: 'D', type: 'number', default: 7, min: 2, max: 20, group: 'Array arrangement', role: 'geometry', design: true, help: 'Rotor diameters between turbines along a row.' },
    { key: 'sy', label: 'Row spacing', unit: 'D', type: 'number', default: 5, min: 2, max: 20, group: 'Array arrangement', role: 'geometry', design: true },
    { key: 'rowAngle', label: 'Row direction (from east, counter-clockwise)', unit: 'deg', type: 'number', default: 30, min: -90, max: 90, group: 'Array arrangement', role: 'geometry', design: true, help: 'Rows along the prevailing wind maximise wake losses; turn them away from the prevailing sectors.' },
    { key: 'stagger', label: 'Stagger alternate rows', type: 'bool', default: false, group: 'Array arrangement', role: 'geometry' },
    site('A', { role: 'bc' }), site('k', { role: 'bc' }),
    { key: 'rose', label: 'Wind rose', type: 'table', group: 'Wind resource', role: 'bc',
      columns: [{ key: 'dir', label: 'Direction (from)', unit: 'deg', type: 'number' }, { key: 'freq', label: 'Frequency', unit: '%', type: 'number' }, { key: 'factor', label: 'Sector A / omnidirectional A', unit: '-', type: 'number' }],
      default: ROSE.map(([dir, freq, factor]) => ({ dir, freq, factor })), help: 'Inflow boundary condition by direction sector (clockwise from north, direction the wind comes from). Frequencies are normalised to 100 %.' },
    { key: 'TI', label: 'Ambient turbulence intensity', unit: '-', type: 'number', default: 0.06, min: 0.02, max: 0.25, group: 'Wind resource', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Hub-height ambient TI; offshore 0.05–0.08. Low TI means slow wake recovery and higher losses.' },
    { key: 'model', label: 'Wake model', type: 'select', default: 'gaussian', group: 'Wake model', role: 'control', options: [{ value: 'gaussian', label: 'Gaussian (Bastankhah & Porté-Agel 2014)' }, { value: 'jensen', label: 'Jensen / Park top-hat' }], help: 'Engineering wake models: no CFD or large-eddy simulation and no dynamic wake meandering. Calibrate against farm operating data where available.' },
    { key: 'kJensen', label: 'Jensen wake decay constant k', unit: '-', type: 'number', default: 0.04, min: 0.02, max: 0.1, group: 'Wake model', role: 'control', help: `0.04–0.05 offshore, 0.075 onshore. ${CONFIRM}`, uncertain: { dist: 'normal', cov: 0.15 } },
    { key: 'superposition', label: 'Wake superposition', type: 'select', default: 'rss', group: 'Wake model', role: 'control', options: [{ value: 'rss', label: 'Root-sum-square of deficits (Katić)' }, { value: 'linear', label: 'Linear sum' }, { value: 'energy', label: 'Energy balance' }, { value: 'max', label: 'Largest deficit only' }] },
    { key: 'blockage', label: 'Include array blockage (summed upstream induction)', type: 'bool', default: false, group: 'Wake model', role: 'control' },
    { key: 'nBins', label: 'Wind-speed bins between cut-in and cut-out', type: 'int', default: 32, min: 6, max: 128, group: 'Solver', role: 'control', help: 'Discretisation of the Weibull distribution in each direction sector (32 bins ≈ 0.7 m/s for a 3–25 m/s machine). Run time grows with bins × directions × turbines².' },
    { key: 'subSectors', label: 'Directions evaluated per sector', type: 'int', default: 3, min: 1, max: 9, group: 'Solver', role: 'control', help: 'One direction per 30° sector over-states row-aligned losses; 3 or more smooths them.' },
    { key: 'maxPerString', label: 'Turbines per array-cable string', type: 'int', default: 6, min: 1, max: 15, group: 'Electrical collection', role: 'geometry', help: 'About 6 × 5 MW or 5 × 15 MW per 66 kV string.' },
    { key: 'wohler', label: 'Wöhler exponent for effective turbulence', unit: '-', type: 'number', default: 4, min: 3, max: 14, group: 'Acceptance', role: 'criteria', help: '3–5 for welded steel (tower, monopile), 10–12 for blade composites.' },
    { key: 'tiClass', label: 'Design turbulence intensity of the turbine class', unit: '-', type: 'number', default: 0.14, min: 0.08, max: 0.2, group: 'Acceptance', role: 'criteria', help: `Representative TI of the type certificate near rated wind speed: IEC class B Iref = 0.14 (A 0.16, C 0.12). ${CONFIRM}` },
    { key: 'minSpacing', label: 'Minimum centre spacing', unit: 'D', type: 'number', default: 4, min: 2, max: 10, group: 'Acceptance', role: 'criteria', help: 'Supplier site-suitability limit; below about 3–4 D wake-induced loads usually exceed the type class.' },
    { key: 'lossTarget', label: 'Allowable wake loss', unit: '%', type: 'number', default: 12, min: 1, max: 40, group: 'Acceptance', role: 'criteria' },
    costInput('energyPrice', 'Value of delivered energy', 'USD/MWh', 80, 1000, 'Used to express the wake loss as revenue.'),
    costInput('cableRate', 'Array cable supplied and installed', 'USD/m', 500, 5000, 'Supply, lay and burial per metre of route.'),
  ],
  mesh: { key: 'nBins', label: 'Wind-speed bins', levels: [16, 32, 64], qoi: 'wakeLoss' },
  compute(inp) {
    const c = getCurve(inp), D = c.D;
    need(inp.hubHeight > D / 2 + 5, 'The hub height leaves no clearance under the rotor.');
    const { pts, imported } = layoutPoints(inp, D), n = pts.length;
    const rows = (inp.rose || []).filter((r) => r.freq > 0);
    need(rows.length >= 1, 'The wind rose needs at least one sector with a positive frequency.');
    need(rows.every((r) => r.factor > 0.2 && r.factor < 3 && r.dir >= 0 && r.dir <= 360), 'Wind rose: directions must be 0–360° and sector factors between 0.2 and 3.');
    const fSum = rows.reduce((a, r) => a + r.freq, 0);
    const rose = rows.map((r) => ({ dirDeg: r.dir, frequency: r.freq / fSum, A: inp.weibA * r.factor, k: inp.weibK }));
    const curve = { u: c.u, p: c.p, ct: c.ct };
    const turbines = pts.map(([x, y]) => ({ x, y, D, hubHeight: inp.hubHeight, curve }));
    // 12-point rotor quadrature: within 0.1 % of the 24-point default on farm efficiency at an eighth of the cost
    const wake = { model: inp.model, k: inp.kJensen, TI: inp.TI, superposition: inp.superposition, blockage: inp.blockage, turbulenceModel: 'crespo', points: rotorPoints(2, 6) };
    const du = (c.cutOut - c.cutIn) / inp.nBins;
    const lay = n > 1 ? checkLayout({ points: pts, minSpacing: inp.minSpacing * D }) : { minDistance: Infinity, spacingViolations: 0 };
    need(n === 1 || lay.minDistance > 1.5 * D, `Two turbines are only ${fmt(lay.minDistance / D)} D apart: rotors would nearly touch. Check the layout coordinates or spacing.`);
    const far = run(() => farmAep({ turbines, rose, du, uMin: c.cutIn + du / 2, uMax: c.cutOut - du / 4, subSectors: n > 1 ? inp.subSectors : 1, wake }), 'Farm energy calculation');
    need(far.grossWh > 0, 'The gross farm energy is zero: the power curve does not overlap the wind climate.');
    const eff = far.perTurbine.map((t) => t.efficiency), effMin = Math.min(...eff), iMin = eff.indexOf(effMin);
    // prevailing direction by energy content (frequency × A³) for the flow-field illustration at a below-rated speed
    const prev = rose.reduce((a, r) => (r.frequency * r.A ** 3 > a.frequency * a.A ** 3 ? r : a), rose[0]);
    const Ushow = Math.min(9, 0.5 * (c.cutIn + c.cutOut));
    const ff = run(() => farmFlow({ turbines: turbines.map((t) => ({ ...t, ct: { u: c.u, v: c.ct }, power: { u: c.u, v: c.p } })), U: Ushow, wdDeg: prev.dirDeg, ...wake }), 'Farm flow field');
    // Frandsen effective turbulence at the most-waked turbine from its eight nearest neighbours (IEC 61400-1 Annex E)
    const ctRated = interp(c.u, c.ct, Math.max(c.cutIn, 0.85 * (c.u[c.p.indexOf(Math.max(...c.p))] || 11)));
    let ieff = inp.TI, nearest = Infinity;
    if (n > 1) {
      const d = pts.map((p, j) => (j === iMin ? Infinity : Math.hypot(p[0] - pts[iMin][0], p[1] - pts[iMin][1]))).sort((a, b) => a - b);
      nearest = lay.minDistance / D;
      const neigh = d.slice(0, Math.min(8, n - 1)).filter((x) => x / D < 10).map((x) => ({ distanceD: x / D, Ct: ctRated }));
      if (neigh.length) ieff = effectiveTurbulence({ I0: inp.TI, neighbours: neigh, m: inp.wohler }).Ieff;
    }
    const cx = pts.reduce((a, p) => a + p[0], 0) / n, cy = pts.reduce((a, p) => a + p[1], 0) / n;
    const cab = n > 1 ? cableLayout({ turbines: pts, substation: [cx, cy], maxPerString: inp.maxPerString, routeFactor: 1.1 }) : { length: 0, strings: [], mstLength: 0 };
    let dens = n > 2 ? capacityDensity({ points: pts, ratedPower: c.ratedPower }) : { hullDensity: 0, hullArea: 0 };
    if (!Number.isFinite(dens.hullDensity)) dens = { hullDensity: 0, hullArea: 0 };   // a single row has no enclosed area
    const deep = imported ? null : deepArrayFrandsen({ Ct: ctRated, sx: inp.sx, sy: inp.sy, hubHeight: inp.hubHeight, z0: 2e-4, H: Math.max(5 * inp.hubHeight, 500) });
    // single-wake centre-line deficit versus distance, both models, at the illustration speed
    const ct8 = interp(c.u, c.ct, Ushow), kStar = wakeGrowthFromTI(inp.TI);
    const xD = linspace(2, 20, 37);
    const jen = xD.map((s) => jensenDeficit({ Ct: ct8, x: s * D, D, k: inp.kJensen }).deficit);
    const gau = xD.map((s) => gaussianDeficit({ Ct: ct8, x: s * D, D, kStar }).centreDeficit);
    const jensenAtSx = jensenDeficit({ Ct: ct8, x: inp.sx * D, D, k: inp.kJensen }).deficit;
    // hub-height speed field for the prevailing direction
    const w = [-Math.sin(prev.dirDeg * DEG), -Math.cos(prev.dirDeg * DEG)];
    const xsAll = pts.map((p) => p[0]), ysAll = pts.map((p) => p[1]);
    const pad = 4 * D, gx = linspace(Math.min(...xsAll) - pad, Math.max(...xsAll) + pad, 56), gy = linspace(Math.min(...ysAll) - pad, Math.max(...ysAll) + pad, 44);
    const field = gy.map((y) => gx.map((x) => {
      const ds = [];
      for (let j = 0; j < n; j++) {
        const dx = x - pts[j][0], dy = y - pts[j][1], xd = dx * w[0] + dy * w[1];
        if (xd <= 0.5 * D) continue;
        const r = Math.abs(dx * w[1] - dy * w[0]);
        if (inp.model === 'jensen') { const q = jensenDeficit({ Ct: ff.Ct[j], x: xd, D, k: inp.kJensen }); if (r <= q.wakeDiameter / 2) ds.push(q.deficit); } else {
          const q = gaussianDeficit({ Ct: ff.Ct[j], x: xd, D, kStar: wakeGrowthFromTI(ff.TI[j]), r }); if (q.deficit > 1e-4) ds.push(q.deficit);
        }
      }
      return round(Ushow * (1 - Math.min(ds.length ? superpose(ds, inp.superposition) : 0, 0.95)), 4);
    }));
    const wakeLoss = far.wakeLoss;
    const checks = [
      check('Wake loss against the project allowance', wakeLoss * 100, inp.lossTarget, '%', 'Project criterion; IEC 61400-15-1 loss category "wake effect"', 'Annual energy lost to wakes inside the array against the allowance in the business case. Above it the layout gives away more energy than planned.'),
      ...(n > 1 ? [check('Minimum turbine spacing', inp.minSpacing, nearest, 'D', 'Supplier site-suitability practice; IEC 61400-1:2019 Annex E', 'Demand is the required minimum spacing, capacity the smallest centre distance in the layout. Closer spacing raises wake-induced fatigue loads beyond the turbine class.'),
        check('Effective turbulence at the most-waked turbine', ieff, inp.tiClass, '-', 'IEC 61400-1:2019 Annex E (Frandsen effective turbulence)', 'Fatigue-equivalent turbulence intensity including neighbour wakes against the design value of the turbine class. Exceeding it shortens the fatigue life of blades, drivetrain and tower unless a site-specific load assessment shows margin.')] : []),
    ];
    const gov = worst(checks);
    const gwh = (v) => round(v / 1e9, 5);
    const lostMwh = (far.grossWh - far.aepWh) / 1e6;
    // sensitivity of the wake loss to a wider in-row spacing is left to the design study; quote the energy value of one point of loss
    const onePoint = far.grossWh * 0.01 / 1e6;
    return {
      kpis: [
        { key: 'wakeLoss', label: 'Wake (array) loss', value: round(wakeLoss * 100, 4), unit: '%', limit: inp.lossTarget, sense: 'max', explain: 'Share of the gross farm energy lost because downstream rotors see slower, more turbulent wind. Large offshore arrays typically lose 5–15 %.' },
        { key: 'arrayEff', label: 'Array efficiency', value: round(far.arrayEfficiency * 100, 4), unit: '%', explain: 'Farm energy with wakes over the energy of the same turbines standing alone. It multiplies the gross yield in the loss cascade.' },
        { key: 'aepGross', label: 'Gross farm energy (no interaction)', value: gwh(far.grossWh), unit: 'GWh/yr', explain: 'Annual energy of all turbines in the undisturbed wind rose, before wakes and other losses.' },
        { key: 'aepArray', label: 'Farm energy with wakes', value: gwh(far.aepWh), unit: 'GWh/yr', explain: 'Array-level annual energy after wake (and, if selected, blockage) effects, before electrical and availability losses.' },
        { key: 'effMin', label: 'Efficiency of the most-waked turbine', value: round(effMin * 100, 4), unit: '%', explain: `Lowest turbine energy relative to its free-stream energy (turbine ${iMin + 1}). Interior positions in line with the prevailing wind lose most.` },
        { key: 'ieff', label: 'Effective turbulence intensity (most-waked turbine)', value: round(ieff, 4), unit: '-', limit: inp.tiClass, sense: 'max', explain: 'Fatigue-equivalent turbulence including wake-added turbulence from the nearest neighbours (Frandsen, IEC 61400-1 Annex E). It must stay below the design turbulence of the turbine class.' },
        { key: 'spacingMin', label: 'Smallest centre spacing', value: n > 1 ? round(nearest, 4) : 0, unit: 'D', explain: 'Shortest distance between two turbines in rotor diameters. 0 is shown for a single turbine.' },
        { key: 'density', label: 'Capacity density', value: round(dens.hullDensity, 4), unit: 'MW/km²', explain: 'Installed capacity over the convex hull of the positions. Offshore farms are built at about 4–8 MW/km²; higher density raises wake loss (0 for a single row, which encloses no area).' },
        { key: 'cableLength', label: 'Inter-array cable length (estimate)', value: round(cab.length / 1e3, 4), unit: 'km', explain: 'Total length of radial strings to a substation at the layout centroid, with 10 % route allowance. It is passed to the array-cable module.' },
        { key: 'jensenAtSx', label: 'Single-wake centre-line deficit at the in-row spacing (Jensen)', value: round(jensenAtSx, 5), unit: '-', explain: 'Fractional speed reduction directly behind one turbine at the in-row spacing by the Jensen model at the illustration wind speed. A quick sanity value: 0.15–0.25 at 7 D offshore.' },
        { key: 'count', label: 'Turbines in the layout', value: n, unit: '-', explain: 'Number of positions used in the farm calculation (imported points or the parametric grid).' },
      ],
      checks,
      plots: [
        { type: 'model3d', title: 'Layout with turbine efficiency', nodes: pts.map((p) => [round(p[0], 2), round(p[1], 2), 0]), lines: cab.edges ? cab.edges.filter((e) => e.from >= 0 && e.to >= 0).map((e) => [e.from, e.to]) : [], tris: [], nodeValues: eff.map((e) => round(e * 100, 3)), valueLabel: 'Energy relative to free stream', valueUnit: '%',
          explain: 'Plan view of the turbine positions coloured by annual energy relative to an isolated turbine, with the array-cable strings between turbines. Look for the low-efficiency positions in the interior and along rows aligned with the prevailing wind.' },
        { type: 'heatmap', title: `Hub-height wind speed, wind from ${fmt(prev.dirDeg)}° at ${fmt(Ushow)} m/s`, xs: gx.map((v) => round(v, 1)), ys: gy.map((v) => round(v, 1)), z: field, zLabel: 'Wind speed', zUnit: 'm/s', x: { label: 'Easting', unit: 'm' }, y: { label: 'Northing', unit: 'm' },
          explain: 'Wake field for the energy-prevailing direction at a below-rated wind speed (highest thrust coefficient). Dark streaks are wakes: turbines standing in them lose power and see added turbulence.' },
        { type: 'polar', title: 'Array efficiency by wind direction', series: [{ name: 'Array efficiency (%)', thetaDeg: far.perSector.map((s) => s.dirDeg), r: far.perSector.map((s) => round(s.efficiency * 100, 3)) }, { name: 'Frequency × 5 (%)', thetaDeg: far.perSector.map((s) => s.dirDeg), r: far.perSector.map((s) => round(s.frequency * 500, 3)) }],
          explain: 'Farm efficiency for each direction sector with the sector frequency. A dip that coincides with a frequent direction is what the row orientation should avoid.' },
        { type: 'bar', title: 'Energy per turbine relative to free stream', categories: eff.map((_, i) => `T${i + 1}`), series: [{ name: 'Efficiency (%)', y: eff.map((e) => round(e * 100, 3)) }], y: { label: 'Efficiency', unit: '%' },
          explain: 'One bar per turbine. The spread between perimeter and interior turbines shows how unevenly the wake loss (and wake-induced fatigue) is shared.' },
        { type: 'line', title: 'Single-wake centre-line deficit versus distance', x: { label: 'Downstream distance', unit: 'D' }, y: { label: 'Velocity deficit', unit: '-' },
          series: [{ name: `Jensen (k = ${inp.kJensen})`, x: xD, y: jen.map((v) => round(v, 4)) }, { name: `Gaussian (k* = ${fmt(kStar)})`, x: xD, y: gau.map((v) => round(v, 4)) }, { name: 'In-row spacing', x: [inp.sx, inp.sx], y: [0, Math.max(...gau, ...jen)], dash: true }],
          explain: 'Speed deficit on the wake axis behind one turbine for both models. The Gaussian centre-line value is higher than the top-hat average close in; both recover with distance — the knee shows where extra spacing stops paying.' },
      ],
      tables: [
        { title: 'Turbines', columns: [{ key: 'i', label: 'Turbine' }, { key: 'x', label: 'Easting', unit: 'm' }, { key: 'y', label: 'Northing', unit: 'm' }, { key: 'e', label: 'Energy with wakes', unit: 'GWh/yr' }, { key: 'eff', label: 'Efficiency', unit: '%' }, { key: 'u', label: `Speed at ${fmt(Ushow)} m/s, prevailing`, unit: 'm/s' }, { key: 'ti', label: 'Local TI (prevailing)', unit: '-' }],
          rows: pts.map((p, i) => ({ i: i + 1, x: round(p[0], 2), y: round(p[1], 2), e: gwh(far.perTurbine[i].aepWh), eff: round(eff[i] * 100, 3), u: round(ff.U[i], 4), ti: round(ff.TI[i], 4) })),
          explain: 'One row per turbine: position, annual energy with wakes, efficiency against free stream, and rotor-averaged speed and turbulence for the prevailing direction. Use it for per-position site-suitability checks.' },
        { title: 'Direction sectors', columns: [{ key: 'dir', label: 'Direction', unit: 'deg' }, { key: 'f', label: 'Frequency', unit: '%' }, { key: 'A', label: 'Weibull A', unit: 'm/s' }, { key: 'g', label: 'Gross', unit: 'GWh/yr' }, { key: 'e', label: 'With wakes', unit: 'GWh/yr' }, { key: 'eff', label: 'Efficiency', unit: '%' }],
          rows: far.perSector.map((s, i) => ({ dir: s.dirDeg, f: round(s.frequency * 100, 3), A: round(rose[i]?.A ?? 0, 4), g: gwh(s.grossWh), e: gwh(s.aepWh), eff: round(s.efficiency * 100, 3) })),
          explain: 'Energy and efficiency per direction sector. The sectors with large gross energy and low efficiency drive the annual wake loss.' },
      ],
      notes: [
        `${n} × ${fmt(c.ratedPower / 1e6)} MW (${c.label}), rotor ${fmt(D, 4)} m, ${imported ? 'imported coordinates' : `grid ${inp.nRows} × ${inp.nCols} at ${inp.sx} D × ${inp.sy} D, rows at ${inp.rowAngle}°`}; ${inp.model === 'jensen' ? 'Jensen' : 'Gaussian'} wake model, ${inp.superposition} superposition${inp.blockage ? ', with blockage' : ''}.`,
        `Wake loss ${fmt(wakeLoss * 100)} %; governing check: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        ...(deep ? [`Frandsen infinite-array asymptote for this spacing: hub-height speed ratio ${fmt(deep.speedRatio)} (power ratio ${fmt(deep.powerRatio)}) — the level a very large farm tends to in its interior.`] : []),
        `Array cables: ${cab.strings.length} strings, ${fmt(cab.length / 1e3, 4)} km including 10 % route allowance (unconstrained minimum spanning tree ${fmt(cab.mstLength / 1e3, 4)} km).`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `Wake loss is ${fmt(wakeLoss * 100)} % (array efficiency ${fmt(far.arrayEfficiency * 100, 4)} %); the most-waked turbine produces ${fmt(effMin * 100)} % of its free-stream energy and sees an effective turbulence of ${fmt(ieff)}.`,
          gov.uc > 1 ? `'${gov.name}' is exceeded (utilisation ${fmt(gov.uc)}): widen the spacing in the prevailing direction, rotate the rows away from ${fmt(prev.dirDeg)}° or remove interior positions, then re-run.` : `Keep the layout; test row direction and spacing with the design study (the prevailing energy direction is ${fmt(prev.dirDeg)}°) and confirm the turbine class against the effective turbulence with the supplier.`,
          'IEC 61400-1:2019 Annex E (wind-farm turbulence) and IEC 61400-15-1 energy-yield practice; engineering wake models are calibrated against farm SCADA and lidar data before financial close.',
          `Wakes cost ${fmt(lostMwh, 4)} MWh a year. ${energyCo2Text(lostMwh)} Wider spacing recovers energy but lengthens array cables (now ${fmt(cab.length / 1e3, 4)} km) and enlarges the seabed footprint.`),
        adv('economic', 'info', `One percentage point of wake loss is ${fmt(onePoint, 4)} MWh a year for this array; the cable estimate is ${fmt(cab.length / 1e3, 4)} km. At the indicative rates entered the wake loss costs ${money(lostMwh * inp.energyPrice)} a year (${money(onePoint * inp.energyPrice)} per point) and the array cable ${money(cab.length * inp.cableRate)}.`,
          'Trade spacing against cable length and lease area with the layout as design variables: each extra rotor diameter of in-row spacing adds roughly one diameter of cable per turbine but recovers energy every year of operation.' + ` Here one point of wake loss a year equals the installed cost of about ${fmt(onePoint * inp.energyPrice / Math.max(inp.cableRate, 1e-9))} m of cable at the indicative rates.`,
          'LCOE-based layout optimisation (IEA Wind Task 37 practice): maximise net energy per unit of whole-life cost, not energy alone; ISO 15686-5 whole-life costing for the cable and loss trade.'),
      ],
      publish: { 'energy.wakeLossFrac': wakeLoss, 'energy.aepArray': far.aepWh, 'array.cableLength': cab.length, 'array.count': n, 'array.effectiveTI': ieff },
      quantities: { cable_km: round(cab.length / 1e3, 5) },
    };
  },
  verify: [
    { name: 'Two turbines in line, power ∝ U, Ct = 0.8, Jensen k = 0.05 at 7 D: deficit 0.1913 → array efficiency (1 + 0.8087)/2',
      source: 'Jensen (1983) / Katić et al. (1986) top-hat wake, hand calculation',
      // deficit = (1 − √(1 − 0.8)) / (1 + 2·0.05·7)² = 0.552786 / 2.89 = 0.191275; wake diameter 1.7 D covers the downstream rotor fully.
      // With P = c·U the downstream turbine yields (1 − 0.191275) of the upstream one at every speed → wake loss = 0.191275 / 2 = 9.564 %
      inputs: { curveSource: 'table', curve: [{ u: 0, p: 0, ct: 0.8 }, { u: 20, p: 2000, ct: 0.8 }, { u: 45, p: 4500, ct: 0.8 }], D: 100, layout: null, nRows: 1, nCols: 2, sx: 7, sy: 5, rowAngle: 0, stagger: false,
        rose: [{ dir: 270, freq: 100, factor: 1 }], subSectors: 1, model: 'jensen', kJensen: 0.05, superposition: 'rss', blockage: false, weibA: 8, weibK: 2 },
      expect: [{ kpi: 'wakeLoss', value: 9.564, tol: 0.005 }, { kpi: 'jensenAtSx', value: 0.191275, tol: 0.002 }] },
    { name: 'Single turbine: no wake interaction (limiting case)', source: 'Limiting case — an isolated rotor has array efficiency 1',
      inputs: { nRows: 1, nCols: 1 }, expect: [{ kpi: 'wakeLoss', value: 0, abs: 1e-9 }, { kpi: 'arrayEff', value: 100, tol: 1e-9 }] },
  ],
};
