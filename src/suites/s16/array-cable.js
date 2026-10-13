// Subsea array-cable thermal and dynamic analysis: IEC 60287 current rating and operating temperature of the string
// cable in its installation, string losses at rated power and over the year, and the static and wave-induced
// tension and curvature of the hanging section (J-tube entry of a fixed foundation or the catenary under a floater).
import { genericArrayCable, ARRAY_CABLE_SIZES } from '../../engine/electrical/collection.js';
import { cableGeometry, rateCable, temperaturesAtCurrent } from '../../engine/electrical/ampacity.js';
import { solveCatenary } from '../../engine/lines/catenary.js';
import { weibull } from '../../engine/aero/aep.js';
import { curveInputs, getCurve } from './_curve.js';
import { costInput, money, site, adv, check, worst, sev, fmt, round, run, need, linspace, interp, G, HOURS_YEAR, LIMIT_CODE, CONFIRM, CO2_COPPER, CO2_ALUMINIUM, CO2_STEEL, CO2_POLYMER, energyCo2Text, vesselCo2Text } from './_common.js';

const cat = { key: 'config', in: ['catenary'] };
const RHO = { Cu: 8890, Al: 2703 };

export default {
  id: 's16-array-cable', suite: 16, title: 'Array cable: thermal rating, losses and hanging section',
  summary: 'Current rating and conductor temperature of the inter-array cable, electrical losses of a turbine string at rated power and over the year, and the tension and bending of the cable where it hangs from the structure or floater.',
  covers: [15], fidelity: 'analytical',
  standards: ['IEC 60287-1-1:2023 and IEC 60287-2-1:2023 (current rating and thermal resistances)', 'IEC 63026:2019 (submarine power cables 6–60 kV) / CIGRE TB 490', 'DNV-ST-0359 (2021) Subsea power cables for wind power plants', 'DNV-RP-0360 (2016) Subsea power cables in shallow water', 'CIGRE TB 862 (2022) Recommendations for mechanical testing of submarine cables for dynamic applications', 'IEC 61400-3-1:2019 / IEC TS 61400-3-2:2019 (electrical collection interface)'],
  models: ['IEC 60287 thermal resistance network with conductor, screen and armour losses', 'Skin and proximity effect, dielectric loss', 'Cable resistance and thermal loss models (I²R with temperature-dependent resistance)', 'Radial string power flow at unity voltage (current proportional to the number of turbines fed)', 'Loss-load factor from the Weibull wind distribution and the power curve', 'Elastic catenary with seabed contact', 'Quasi-static wave-induced tension from heave acceleration and vessel offset'],
  limits: `Generic three-core array-cable construction (indicative: replace with the manufacturer data sheet and type-test results). IEC 60287 over-estimates the armour loss of three-core submarine cables (CIGRE TB 880), so the rating is conservative. The hanging section is analysed quasi-statically in one plane: no hydrodynamic drag, vortex-induced vibration, bend stiffener or fatigue of the dynamic cable (see the cable suite for dynamic and fatigue analysis). ${LIMIT_CODE}`,
  publishes: ['cable.length', 'cable.ampacity', 'cable.lossesFrac', 'cable.topTension', 'quantities.cable_km'],
  searchTerms: 'inter-array cable 66 kV ampacity IEC 60287 string losses dynamic cable floating wind',
  inputs: [
    ...curveInputs({ group: 'Turbine and string' }),
    { key: 'nString', label: 'Turbines on the string', type: 'int', default: 5, min: 1, max: 15, group: 'Turbine and string', role: 'geometry', design: true, help: 'The cable section next to the substation carries all of them.' },
    { key: 'spacing', label: 'Cable length between turbines', unit: 'm', type: 'number', default: 950, min: 100, max: 5000, group: 'Turbine and string', role: 'geometry', help: 'Route length per section including J-tube or hang-off lengths and burial slack.' },
    { key: 'totalKm', label: 'Total array-cable length of the farm', unit: 'km', type: 'number', default: 60, min: 0.1, max: 2000, group: 'Turbine and string', role: 'geometry', bind: 'array.cableLength', help: 'From the wake and array module (string layout) when it has been run.' },
    { key: 'route', label: 'Imported cable routes (optional)', type: 'geometry', kinds: ['polyline'], default: null, group: 'Turbine and string', role: 'geometry', help: 'CSV / DXF / GeoJSON polyline(s) of the array-cable routes in metres. The summed length replaces the total length above; leave empty to use the entered value.' },
    { key: 'U', label: 'Array voltage (line-to-line)', unit: 'kV', type: 'select', default: '33', group: 'Electrical system', role: 'bc', options: [{ value: '33', label: '33 kV' }, { value: '66', label: '66 kV' }], help: 'Receiving-network boundary condition at the offshore substation busbar.' },
    { key: 'pf', label: 'Power factor', unit: '-', type: 'number', default: 0.95, min: 0.8, max: 1, group: 'Electrical system', role: 'bc' },
    { key: 'area', label: 'Conductor cross-section', unit: 'mm²', type: 'select', default: '240', group: 'Cable', role: 'geometry', design: true, options: ARRAY_CABLE_SIZES.map((a) => ({ value: String(Math.round(a * 1e6)), label: `${Math.round(a * 1e6)} mm²` })) },
    { key: 'metal', label: 'Conductor metal', type: 'select', default: 'Cu', group: 'Cable', role: 'material', options: [{ value: 'Cu', label: 'Copper' }, { value: 'Al', label: 'Aluminium' }] },
    { key: 'thetaMax', label: 'Maximum conductor temperature', unit: '°C', type: 'number', default: 90, min: 60, max: 105, group: 'Cable', role: 'criteria', help: `90 °C for XLPE (IEC 63026). ${CONFIRM}` },
    { key: 'burial', label: 'Burial depth to top of cable', unit: 'm', type: 'number', default: 1.5, min: 0.2, max: 6, group: 'Installation', role: 'ic', design: true },
    { key: 'rhoSoil', label: 'Soil thermal resistivity', unit: 'K·m/W', type: 'number', default: 0.7, min: 0.3, max: 3, group: 'Installation', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Saturated sand 0.4–0.7, clay 0.7–1.0.' },
    { key: 'Tsea', label: 'Seabed / sea-water temperature', unit: '°C', type: 'number', default: 12, min: -2, max: 32, group: 'Installation', role: 'bc', bind: 'site.Tsea', help: 'Thermal boundary condition (ambient).' },
    site('A', { group: 'Wind resource', role: 'bc' }), site('k', { group: 'Wind resource', role: 'bc' }),
    { key: 'config', label: 'Hanging section', type: 'select', default: 'jtube', group: 'Hanging section', role: 'bc', options: [{ value: 'jtube', label: 'Fixed foundation: vertical hang in the J-tube' }, { value: 'catenary', label: 'Floater: free-hanging catenary to the seabed' }] },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 30, min: 5, max: 1500, group: 'Hanging section', role: 'environment', bind: 'site.depth' },
    { key: 'hang', label: 'Hang-off height above still water', unit: 'm', type: 'number', default: 15, min: -30, max: 40, group: 'Hanging section', role: 'bc', help: 'Level of the hang-off flange (negative = below the water line, as on a floater).' },
    { key: 'span', label: 'Horizontal distance hang-off to seabed anchor point', unit: 'm', type: 'number', default: 120, min: 5, max: 3000, group: 'Hanging section', role: 'bc', showIf: cat },
    { key: 'hangLength', label: 'Cable length hang-off to anchor point', unit: 'm', type: 'number', default: 0, min: 0, max: 5000, group: 'Hanging section', role: 'geometry', showIf: cat, help: '0 = 15 % longer than the straight distance.' },
    { key: 'offset', label: 'Floater offset (± about the mean position)', unit: 'm', type: 'number', default: 15, min: 0, max: 200, group: 'Hanging section', role: 'bc', showIf: cat, bind: 'vessel.offsetMax' },
    { key: 'heave', label: 'Hang-off heave amplitude', unit: 'm', type: 'number', default: 2.5, min: 0, max: 15, group: 'Hanging section', role: 'bc', showIf: cat, bind: 'vessel.heaveSig' },
    { key: 'Tw', label: 'Wave period for the dynamic tension', unit: 's', type: 'number', default: 10, min: 3, max: 25, group: 'Hanging section', role: 'bc', showIf: cat, bind: 'site.Tp1' },
    { key: 'massOver', label: 'Cable mass in air (0 = from the construction)', unit: 'kg/m', type: 'number', default: 0, min: 0, max: 150, group: 'Hanging section', role: 'geometry', help: 'Use the data-sheet value when available.' },
    { key: 'odOver', label: 'Cable outer diameter (0 = from the construction)', unit: 'mm', type: 'number', default: 0, min: 0, max: 300, group: 'Hanging section', role: 'geometry' },
    { key: 'Tallow', label: 'Allowable tension of the cable', unit: 'kN', type: 'number', default: 120, min: 5, max: 2000, group: 'Acceptance', role: 'criteria', help: 'Maximum permissible tension from the data sheet (installation and operation), typically 100–250 kN for armoured array cables.' },
    { key: 'mbrFactor', label: 'Minimum bending radius / outer diameter', unit: '-', type: 'number', default: 15, min: 8, max: 30, group: 'Acceptance', role: 'criteria', help: 'About 15 × OD under tension for three-core armoured cables (data-sheet value governs).' },
    costInput('energyPrice', 'Value of delivered energy', 'USD/MWh', 80, 1000, 'Used to capitalise the cable losses.'),
    costInput('metalPrice', 'Conductor metal price', 'USD/t', 9500, 1e5, 'Per tonne of conductor metal as cabled: order of 9 500 for copper, 2 600 for aluminium.'),
    costInput('layDayRate', 'Cable-lay vessel spread day rate', 'USD/day', 120000, 2e6),
  ],
  compute(inp) {
    const c = getCurve(inp), U = Number(inp.U) * 1e3, area = Number(inp.area) * 1e-6;
    const desc = genericArrayCable({ voltage: U, area, material: inp.metal });
    const g = run(() => cableGeometry(desc), 'Cable construction');
    need(inp.thetaMax > inp.Tsea + 10, 'The maximum conductor temperature must exceed the ambient temperature by at least 10 K.');
    const installation = { type: 'buried', depth: inp.burial + g.De / 2, rhoSoil: inp.rhoSoil, ambientC: inp.Tsea, formation: 'single' }, system = { f: 50, U, thetaMaxC: inp.thetaMax };
    const rate = run(() => rateCable({ cable: g, installation, system }), 'Current rating');
    const I1 = c.ratedPower / (Math.sqrt(3) * U * inp.pf), Imax = I1 * inp.nString;
    let load;
    try { load = temperaturesAtCurrent({ cable: g, installation, system, I: Imax }); } catch { throw new Error(`No thermal equilibrium at ${fmt(Imax)} A: the string current is far above the cable rating of ${fmt(rate.I)} A. Choose a larger conductor or fewer turbines per string.`); }
    // string losses: section k (counted from the far end) carries k turbines; resistance and loss factors at the temperature of each section
    const secs = [];
    let lossRated = 0;
    for (let kk = 1; kk <= inp.nString; kk++) {
      const tk = temperaturesAtCurrent({ cable: g, installation, system, I: I1 * kk });
      const P = tk.losses.total * inp.spacing; lossRated += P;
      secs.push({ k: kk, I: I1 * kk, T: tk.conductorC, wpm: tk.losses.total, P });
    }
    const Pstring = c.ratedPower * inp.nString, lossFracRated = lossRated / Pstring;
    // annual loss: current-dependent losses scale with (P/Prated)²; dielectric loss is constant
    const wb = weibull(inp.weibA, inp.weibK), us = linspace(0, 35, 351);
    let e1 = 0, e2 = 0;
    for (let i = 1; i < us.length; i++) { const u = 0.5 * (us[i] + us[i - 1]), pr = (u < c.cutIn || u > c.cutOut ? 0 : interp(c.u, c.p, u) / c.ratedPower), w = wb.cdf(us[i]) - wb.cdf(us[i - 1]); e1 += pr * w; e2 += pr * pr * w; }
    const Wd3 = (load.losses.dielectric ?? 0) * inp.spacing * inp.nString;
    const lossAnnualWh = ((lossRated - Wd3) * e2 + Wd3) * HOURS_YEAR, genWh = Pstring * e1 * HOURS_YEAR;
    const lossFracAnnual = genWh > 0 ? lossAnnualWh / genWh : 0;
    // mass per metre from the construction
    const A = Object.fromEntries(g.layers.map((L) => [L.name, L.area]));
    const mCond = 3 * (A.conductor || 0) * RHO[inp.metal], mScreen = 3 * (A.sheath || 0) * RHO.Cu, mArm = (A.armour || 0) * 7850;
    const mPoly = 3 * ((A.insulation || 0) * 1000 + (A.coreJacket || 0) * 950) + (A.bedding || 0) * 600 + (A.serving || 0) * 950;
    const mass = inp.massOver > 0 ? inp.massOver : mCond + mScreen + mArm + mPoly, OD = inp.odOver > 0 ? inp.odOver / 1e3 : g.De;
    const wSub = (mass - 1025 * Math.PI * OD * OD / 4) * G;
    need(wSub > 0, 'The cable is buoyant with these mass and diameter values: check the overrides.');
    // hanging section
    let Ttop, Tdyn = 0, mbr = Infinity, shape = null, near = null, far = null, Lh = 0;
    if (inp.config === 'jtube') {
      need(inp.hang >= 0, 'A J-tube hang-off is above the water line: enter a positive hang-off height or select the catenary configuration.');
      Ttop = wSub * inp.depth + mass * G * inp.hang;
    } else {
      const Z = inp.depth + inp.hang;
      need(Z > 5, 'The hang-off must be at least 5 m above the seabed.');
      Lh = inp.hangLength > 0 ? inp.hangLength : 1.15 * Math.hypot(inp.span, Z);
      const solve = (X) => { need(Lh > Math.hypot(X, Z) * 1.002, `The cable (${fmt(Lh, 4)} m) is too short for a hang-off ${fmt(Z, 4)} m above the seabed and ${fmt(X, 4)} m away from the anchor point (straight distance ${fmt(Math.hypot(X, Z), 4)} m): lengthen it or reduce the offset.`); need(Lh < X + Z, `The cable (${fmt(Lh, 4)} m) is so long that it piles up on the seabed (limit ${fmt(X + Z, 4)} m): shorten it.`); return run(() => solveCatenary({ X, Z, L: Lh, w: wSub, EA: 4e8, CB: 0.3 }), 'Catenary solution'); };
      const mean = solve(inp.span); far = solve(inp.span + inp.offset); near = inp.span - inp.offset > 5 ? solve(inp.span - inp.offset) : mean;
      // wave-induced tension: suspended weight times heave acceleration / g at the far position
      const acc = inp.heave * (2 * Math.PI / inp.Tw) ** 2;
      Tdyn = far.T * acc / G; Ttop = far.T + Tdyn;
      // smallest radius of curvature of a catenary is at the touchdown: R = H / w (near position governs)
      mbr = Math.min(mean.H, near.H, far.H) / wSub;
      const a = mean.H / wSub, xs = linspace(0, a * Math.acosh(1 + Z / a), 40);
      shape = { x: xs.map((x) => round(x + mean.groundedLength, 3)), z: xs.map((x) => round(a * (Math.cosh(x / a) - 1) - inp.depth, 3)), grounded: mean.groundedLength, Tmean: mean.T };
    }
    const mbrAllow = inp.mbrFactor * OD;
    const checks = [
      check('String current against the continuous rating', Imax, rate.I, 'A', 'IEC 60287-1-1:2023 §1.4 (permissible current rating); DNV-ST-0359 §4', 'Current of the most loaded section (all turbines at rated power) against the current that brings the conductor to its maximum temperature in this burial. Above it the insulation ages faster than its design life assumes.'),
      check('Conductor temperature at rated power', load.conductorC - inp.Tsea, inp.thetaMax - inp.Tsea, 'K', 'IEC 60287-1-1 / IEC 63026 (maximum conductor temperature of XLPE)', 'Temperature rise of the conductor above the seabed at the full string current against the permissible rise.'),
      check('Top tension of the hanging section', Ttop / 1e3, inp.Tallow, 'kN', 'DNV-ST-0359 (2021) §4.6 (mechanical design); CIGRE TB 862', 'Largest tension at the hang-off (static weight plus wave-induced part) against the allowable tension of the cable. Over-tension damages the armour termination and stretches the conductors.'),
      ...(inp.config === 'catenary' ? [check('Bending radius at the touchdown against the minimum bending radius', mbrAllow, mbr, 'm', 'DNV-ST-0359 (2021) §4.6; cable data sheet (MBR under tension)', 'Demand is the minimum permissible bending radius, capacity the smallest radius of the catenary at the seabed in the near position. Tighter bending buckles the armour and cracks the lead-free screens over time.')] : []),
    ];
    const gov = worst(checks);
    // imported route: 3-D length of the polyline replaces the entered total
    let km = inp.totalKm;
    const rp = inp.route?.points;
    if (rp && rp.length >= 6) { let Lr = 0; for (let i = 3; i + 2 < rp.length; i += 3) Lr += Math.hypot(rp[i] - rp[i - 3], rp[i + 1] - rp[i - 2], rp[i + 2] - rp[i - 1]); need(Number.isFinite(Lr) && Lr > 10, 'The imported route has no usable length: check the coordinates and units (metres).'); km = Lr / 1e3; }
    const condT = mCond * km, steelT = mArm * km, polyT = mPoly * km, cuScreenT = mScreen * km;
    const co2 = condT * (inp.metal === 'Cu' ? CO2_COPPER : CO2_ALUMINIUM) + cuScreenT * CO2_COPPER + steelT * CO2_STEEL + polyT * CO2_POLYMER;
    // next larger section for the loss trade-off
    const sizes = ARRAY_CABLE_SIZES.map((a) => Math.round(a * 1e6)), iS = sizes.indexOf(Number(inp.area)), nextA = sizes[Math.min(iS + 1, sizes.length - 1)];
    const lossSave = lossAnnualWh * (1 - Number(inp.area) / nextA) / 1e6, extraCu = 3 * (nextA - Number(inp.area)) * 1e-6 * RHO[inp.metal] * inp.spacing * inp.nString / 1e3;
    const Is = linspace(0, Math.max(rate.I * 1.1, Imax * 1.05), 23), Ts = Is.map((I) => { try { return temperaturesAtCurrent({ cable: g, installation, system, I }).conductorC; } catch { return null; } });
    const layDays = km / 4;
    return {
      kpis: [
        { key: 'ampacity', label: 'Continuous current rating', value: round(rate.I, 5), unit: 'A', explain: 'Current that brings the conductor to its maximum temperature in this burial at 100 % load factor (IEC 60287). Wind-farm load is cyclic, so the dynamic rating is somewhat higher.' },
        { key: 'stringCurrent', label: 'Current of the most loaded section at rated power', value: round(Imax, 5), unit: 'A', limit: round(rate.I, 4), sense: 'max', explain: 'All turbines of the string at rated power: n·P/(√3·U·cos φ).' },
        { key: 'conductorT', label: 'Conductor temperature at rated power', value: round(load.conductorC, 4), unit: '°C', limit: inp.thetaMax, sense: 'max', explain: 'Steady conductor temperature of the most loaded section. Every 8–10 K above the limit roughly halves insulation life.' },
        { key: 'lossRated', label: 'String loss at rated power', value: round(lossRated / 1e3, 5), unit: 'kW', explain: 'Conductor, screen, armour and dielectric losses of all sections with every turbine at rated power.' },
        { key: 'lossFracRated', label: 'Loss at rated power / string power', value: round(lossFracRated * 100, 4), unit: '%', explain: 'Share of the string output dissipated in the cable at full load; array systems are usually designed for 1–2 %.' },
        { key: 'lossFracAnnual', label: 'Annual energy loss of the string', value: round(lossFracAnnual * 100, 4), unit: '%', explain: 'Energy lost over the year relative to the energy generated, using the wind distribution and power curve (losses scale with the square of power). This is the electrical array loss for the yield assessment.' },
        { key: 'lossAnnual', label: 'Annual energy loss of the string', value: round(lossAnnualWh / 1e6, 5), unit: 'MWh/yr', explain: 'Energy dissipated in the string cable per year.' },
        { key: 'wpm', label: 'Heat dissipated by the most loaded section', value: round(load.losses.total, 4), unit: 'W/m', explain: 'Total loss per metre of cable at rated power: the heat the seabed must carry away.' },
        { key: 'massPerM', label: 'Cable mass in air', value: round(mass, 4), unit: 'kg/m', explain: 'Mass per metre from the construction (or the entered data-sheet value); with the diameter it gives the submerged weight.' },
        { key: 'wSub', label: 'Submerged weight', value: round(wSub, 4), unit: 'N/m', explain: 'Weight per metre in sea water: it sets the tension of the hanging section and the on-bottom stability.' },
        { key: 'topTension', label: 'Maximum top tension of the hanging section', value: round(Ttop / 1e3, 5), unit: 'kN', limit: inp.Tallow, sense: 'max', explain: inp.config === 'jtube' ? 'Weight of the cable hanging in the J-tube from the hang-off to the seabed.' : 'Hang-off tension at the far offset plus the wave-induced part from heave acceleration.' },
        { key: 'mbr', label: 'Smallest bending radius in the hanging section', value: round(Number.isFinite(mbr) ? mbr : 0, 4), unit: 'm', explain: 'Radius of curvature at the touchdown of the catenary in the near position (0 for the straight J-tube hang; the bellmouth and bend restrictor govern there).' },
        { key: 'cableKm', label: 'Array-cable length of the farm', value: round(km, 5), unit: 'km', explain: 'Total route length carried to the quantity and cost take-off.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Current and conductor temperature along the string', categories: secs.map((s) => `Section ${s.k}`), series: [{ name: 'Current (A)', y: secs.map((s) => round(s.I, 4)) }, { name: 'Conductor temperature × 5 (°C)', y: secs.map((s) => round(s.T * 5, 3)) }], y: { label: 'A · °C×5', unit: '-' },
          explain: 'Each section from the far end to the substation carries one more turbine. Only the last sections run warm: tapering the conductor size along the string saves metal.' },
        { type: 'line', title: 'Conductor temperature versus current', x: { label: 'Current', unit: 'A' }, y: { label: 'Conductor temperature', unit: '°C' }, series: [{ name: 'Conductor', x: Is.filter((_, i) => Ts[i] != null).map((v) => round(v, 3)), y: Ts.filter((v) => v != null).map((v) => round(v, 3)) }, { name: 'Limit', x: [0, Is[Is.length - 1]], y: [inp.thetaMax, inp.thetaMax], dash: true }, { name: 'String current', x: [Imax, Imax], y: [inp.Tsea, inp.thetaMax], dash: true }],
          explain: 'Steady conductor temperature against load current for this burial. The curve steepens because resistance rises with temperature; its crossing of the limit is the rating.' },
        { type: 'bar', title: 'Loss breakdown of the most loaded section', categories: ['Conductor', 'Screens', 'Armour', 'Dielectric'], series: [{ name: 'W/m', y: [load.losses.conductor, load.losses.sheath, load.losses.armour, load.losses.dielectric].map((v) => round(v ?? 0, 4)) }], y: { label: 'Loss', unit: 'W/m' },
          explain: 'Where the heat is generated in the cable at rated power. Armour loss is over-estimated by IEC 60287 for three-core cables, which makes the rating conservative.' },
        ...(shape ? [{ type: 'line', title: 'Hanging cable (mean position)', x: { label: 'Distance from the seabed anchor point', unit: 'm' }, y: { label: 'Elevation', unit: 'm' }, series: [{ name: 'Suspended cable', x: shape.x, y: shape.z }, { name: 'On the seabed', x: [0, round(shape.grounded, 3)], y: [-inp.depth, -inp.depth] }],
          explain: 'Catenary of the cable from the touchdown to the hang-off. The sharpest curvature is at the touchdown; a lazy-wave shape with buoyancy modules is used when tension or curvature limits cannot be met.' }] : []),
      ],
      tables: [
        { title: 'String sections', columns: [{ key: 'k', label: 'Section (from the far end)' }, { key: 'I', label: 'Current', unit: 'A' }, { key: 'T', label: 'Conductor temperature', unit: '°C' }, { key: 'w', label: 'Loss', unit: 'W/m' }, { key: 'P', label: 'Section loss', unit: 'kW' }],
          rows: secs.map((s) => ({ k: s.k, I: round(s.I, 4), T: round(s.T, 3), w: round(s.wpm, 4), P: round(s.P / 1e3, 4) })),
          explain: 'Current, conductor temperature and loss of each cable section at rated power. The last row is the section that sizes the cable.' },
        { title: 'Cable construction and quantities', columns: [{ key: 'q', label: 'Item' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [{ q: 'Construction', v: `3 × ${inp.area} mm² ${inp.metal}, ${inp.U} kV XLPE, steel-wire armoured`, u: '' }, { q: 'Outer diameter', v: round(OD * 1e3, 4), u: 'mm' }, { q: 'AC resistance at maximum temperature', v: round(rate.R * 1e6, 5), u: 'mΩ/km' }, { q: 'Screen / armour loss factors λ1 / λ2', v: `${fmt(rate.lambda1)} / ${fmt(rate.lambda2)}`, u: '-' },
            { q: `Conductor ${inp.metal === 'Cu' ? 'copper' : 'aluminium'} (farm)`, v: round(condT, 4), u: 't' }, { q: 'Screen copper (farm)', v: round(cuScreenT, 4), u: 't' }, { q: 'Armour steel (farm)', v: round(steelT, 4), u: 't' }, { q: 'Polymers (farm)', v: round(polyT, 4), u: 't' }],
          explain: 'Generic construction used for the rating and the material take-off for the whole array-cable length.' },
      ],
      notes: [
        `${inp.nString} × ${fmt(c.ratedPower / 1e6)} MW at ${inp.U} kV: ${fmt(Imax)} A on a 3 × ${inp.area} mm² ${inp.metal} cable rated ${fmt(rate.I)} A buried ${inp.burial} m in soil of ${inp.rhoSoil} K·m/W at ${inp.Tsea} °C; conductor ${fmt(load.conductorC)} °C.`,
        `Losses: ${fmt(lossRated / 1e3)} kW at rated power (${fmt(lossFracRated * 100)} %), ${fmt(lossFracAnnual * 100)} % of the annual energy of the string.`,
        inp.config === 'jtube' ? `J-tube hang: ${fmt(Ttop / 1e3)} kN at the hang-off (${fmt(inp.depth)} m of water, ${inp.hang} m above the water line).` : `Catenary of ${fmt(Lh, 4)} m: mean top tension ${fmt(shape.Tmean / 1e3)} kN, far position ${fmt(far.T / 1e3)} kN plus ${fmt(Tdyn / 1e3)} kN wave-induced; touchdown radius ${fmt(mbr)} m against ${fmt(mbrAllow)} m allowed.`,
        `Governing check: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `The most loaded section carries ${fmt(Imax)} A against a rating of ${fmt(rate.I)} A (conductor ${fmt(load.conductorC)} °C); top tension ${fmt(Ttop / 1e3)} kN against ${inp.Tallow} kN allowed.`,
          checks[0].uc > 1 || checks[1].uc > 1 ? `Use the next conductor size (${nextA} mm²) for the sections next to the substation, reduce the string to ${Math.max(1, Math.floor(rate.I / I1))} turbines, or check the hot spots (J-tube, landfall) with a dynamic rating.` : gov.uc > 1 ? (inp.config === 'catenary' ? 'Change the hanging configuration: a longer cable or a lazy wave with buoyancy modules lowers tension and opens the touchdown radius; add a bend stiffener at the hang-off.' : 'Reduce the hanging weight (support in the J-tube) or select a cable with a higher allowable tension.') : `The cable is thermally and mechanically adequate with ${fmt((1 - checks[0].uc) * 100)} % current margin: taper the conductor along the string (the first ${Math.max(1, inp.nString - 2)} sections need less) and verify the J-tube and hang-off hot spots.`,
          'IEC 60287 (2023) rating with DNV-ST-0359 (2021) mechanical criteria; CIGRE TB 880 on armour-loss conservatism; dynamic cables per CIGRE TB 862.',
          `Array cable for the farm: ${fmt(condT)} t of conductor metal, ${fmt(steelT)} t of armour steel and ${fmt(polyT)} t of polymers ≈ ${fmt(co2)} t CO₂e (screening factors: copper ${CO2_COPPER}, aluminium ${CO2_ALUMINIUM}, steel ${CO2_STEEL}, polymer ${CO2_POLYMER} t CO₂e/t). Aluminium conductors cut conductor mass by about half at equal rating.`),
        adv('economic', 'info', `The string loses ${fmt(lossAnnualWh / 1e6, 4)} MWh a year (${fmt(lossFracAnnual * 100)} % of its energy); the farm needs ${fmt(km, 4)} km of cable, about ${fmt(layDays)} cable-lay days at 4 km/day. At the indicative rates entered the string losses are worth ${money(lossAnnualWh / 1e6 * inp.energyPrice)} a year and the lay campaign ${money(layDays * inp.layDayRate)}.`,
          `Economic conductor sizing: stepping to ${nextA} mm² would save about ${fmt(lossSave)} MWh a year on this string for ${fmt(extraCu)} t more conductor metal — worth it when the capitalised value of the losses over the design life exceeds the metal cost. At the indicative rates that is ${money(lossSave * inp.energyPrice)} a year saved for ${money(extraCu * inp.metalPrice)} of metal${lossSave * inp.energyPrice > 0 ? `: a simple payback of ${fmt(extraCu * inp.metalPrice / (lossSave * inp.energyPrice))} years` : ''}.`,
          'IEC 60287-3-2 (economic optimisation of cable size) and whole-life costing per ISO 15686-5: capitalised loss value against conductor cost.',
          `${energyCo2Text(lossAnnualWh / 1e6)} ${vesselCo2Text(layDays, 25)}`),
      ],
      publish: { 'cable.length': km * 1e3, 'cable.ampacity': rate.I, 'cable.lossesFrac': lossFracAnnual, 'cable.topTension': Ttop, 'quantities.cable_km': km },
      quantities: { cable_km: round(km, 5), ...(inp.metal === 'Cu' ? { copper_t: round(condT + cuScreenT, 4) } : { aluminium_t: round(condT, 4), copper_t: round(cuScreenT, 4) }), steel_t: round(steelT, 4), polymer_t: round(polyT, 4), vessel_days: { cable_lay: round(layDays, 3) } },
    };
  },
  verify: [
    { name: 'String current n·P/(√3·U·cos φ) and J-tube hang-off tension w_sub·d + m·g·h', source: 'Three-phase power relation and statics of a vertically hanging cable (hand calculation)',
      // 6 × 5 MW at 66 kV, cos φ = 0.95: 30e6/(1.73205 × 66e3 × 0.95) = 276.24 A
      // m = 30 kg/m, OD = 120 mm: w_sub = (30 − 1025·π·0.12²/4)·9.80665 = 180.51 N/m; T = 180.51 × 30 + 30 × 9.80665 × 15 = 9.828 kN
      inputs: { curveSource: 'cubic', ratedCubicMW: 5, D: 126, nString: 6, U: '66', pf: 0.95, config: 'jtube', depth: 30, hang: 15, massOver: 30, odOver: 120 }, expect: [{ kpi: 'stringCurrent', value: 276.24, tol: 0.001 }, { kpi: 'topTension', value: 9.828, tol: 0.001 }] },
    { name: 'Catenary with grounded length: top tension = H + w·h and touchdown radius = H / w', source: 'Inextensible catenary (Faltinsen 1990, Eq. 8.20); elastic stretch below 0.1 % here',
      // m = 30 kg/m, OD = 120 mm → w = 180.51 N/m; h = 85 m, X = 120 m, L = 160 m. Catenary parameter a = H/w from L − s + a·asinh(s/a) = X with s = √(h² + 2ah):
      // a = 60.9 m → s = √(7225 + 10353) = 132.58 m, a·asinh(s/a) = 60.9 × 1.5201 = 92.57 m, grounded 27.42 m, sum 119.99 m ✓
      // touchdown radius R = a = 60.9 m; top tension T = w·(a + h) = 180.51 × 145.9 = 26.34 kN (no heave, no offset)
      inputs: { config: 'catenary', depth: 100, hang: -15, span: 120, hangLength: 160, offset: 0, heave: 0, massOver: 30, odOver: 120, Tallow: 200 }, expect: [{ kpi: 'topTension', value: 26.34, tol: 0.01 }, { kpi: 'mbr', value: 60.9, tol: 0.01 }] },
  ],
};
