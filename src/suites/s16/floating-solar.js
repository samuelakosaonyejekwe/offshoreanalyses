// Floating photovoltaic plant: hourly energy yield from the site irradiance (transposition, cell temperature with
// floating-specific coefficients, single-diode module model, inverter and losses), and the float platform —
// buoyancy and freeboard, island stability, wind, wave and current loads, mooring load sharing and connectors.
import { syntheticYear, simulatePv, extractSingleDiode, singleDiodeCurve, singleDiodeAt, cellTempFaiman } from '../../engine/aero/pv.js';
import { floatBuoyancy, islandStability, panelWindLoad, floatWaveLoad, meanDriftForce, currentDrag, mooringLoadSharing, anchorCount, hingedRaftResponse, coolingYieldGain, exposureClass } from '../../engine/aero/floatingpv.js';
import { dispersion } from '../../engine/aero/wec.js';
import { adv, check, worst, sev, fmt, round, run, need, linspace, DEG, G, LIMIT_PV, LIMIT_CODE, CONFIRM, CO2_POLYMER, CO2_ALUMINIUM, CO2_GRID, energyCo2Text } from './_common.js';

const MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], DAYS = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], MID = [17, 47, 75, 105, 135, 162, 198, 228, 258, 288, 318, 344];

/** Monthly mean daily extraterrestrial insolation on a horizontal plane [Wh/m²/day] (Duffie & Beckman, Eq. 1.10.3). */
function h0Monthly(latDeg) {
  const phi = latDeg * DEG;
  return MID.map((nDay) => {
    const dec = 23.45 * DEG * Math.sin(2 * Math.PI * (284 + nDay) / 365), x = -Math.tan(phi) * Math.tan(dec), ws = x >= 1 ? 0 : x <= -1 ? Math.PI : Math.acos(x);
    return 24 / Math.PI * 1361 * (1 + 0.033 * Math.cos(2 * Math.PI * nDay / 365)) * (Math.cos(phi) * Math.cos(dec) * Math.sin(ws) + ws * Math.sin(phi) * Math.sin(dec));
  });
}

export default {
  id: 's16-floating-solar', suite: 16, title: 'Floating solar: energy yield and float platform',
  summary: 'Annual, monthly and daily energy of a floating photovoltaic plant with module and thermal models, and the buoyancy, stability, environmental loads and mooring of the float island.',
  covers: [13, 7], fidelity: 'numerical',
  standards: ['DNV-RP-0584 (2021) Design, development and operation of floating solar photovoltaic systems', 'IEC 61853-1/-2 (PV module performance testing and energy rating)', 'IEC 61724-1:2021 (PV system performance: performance ratio)', 'IEC TS 62738:2018 (design guidelines for PV power plants)', 'EN 1991-1-4:2005 §7.3 (wind on monopitch canopies)', 'DNV-RP-C205 (2021) §6 (wave and current loads)'],
  models: ['Solar position and extraterrestrial irradiance', 'Erbs diffuse decomposition and Hay–Davies / isotropic transposition', 'Faiman cell-temperature model with floating heat-loss coefficients', 'Photovoltaic single-diode equivalent circuit (De Soto five-parameter model)', 'Photovoltaic thermal derating (PVWatts DC model)', 'Inverter clipping and efficiency, solar array DC/AC loss chain', 'Buoyancy stability: Archimedes draft and freeboard, metacentric height and heel of a float island', 'Wind load on tilted panel rows with shelter factor', 'Long-wave Morison force and Maruo mean drift on floats', 'Elastic mooring load sharing with slack-line removal', 'Hinged-raft wave response (connector rotation and shear)'],
  limits: `${LIMIT_PV} The synthetic year uses one clearness index per month (no day-to-day variability or measured diffuse fraction); the double-diode model, module-level mismatch from wave-induced tilt and soiling by birds or salt are represented only through loss inputs. The island is rigid for stability and mooring load sharing; the hinged-raft response is a quasi-static screening; wind-pressure and shelter coefficients are inputs from code or wind-tunnel data. ${LIMIT_CODE}`,
  publishes: ['fpv.ratedPower', 'fpv.aep', 'fpv.capacityFactor', 'fpv.specificYield', 'fpv.mooringLoad'],
  searchTerms: 'floating photovoltaic energy yield single-diode Faiman cooling float mooring DNV-RP-0584',
  inputs: [
    { key: 'pdcMW', label: 'Installed DC capacity', unit: 'MWp', type: 'number', default: 20, min: 0.05, max: 1000, group: 'PV plant', role: 'geometry', design: true },
    { key: 'Isc', label: 'Module short-circuit current', unit: 'A', type: 'number', default: 13.9, min: 1, max: 25, group: 'PV module (datasheet, STC)', role: 'material' },
    { key: 'Voc', label: 'Module open-circuit voltage', unit: 'V', type: 'number', default: 49.6, min: 10, max: 120, group: 'PV module (datasheet, STC)', role: 'material' },
    { key: 'Imp', label: 'Current at maximum power', unit: 'A', type: 'number', default: 13.2, min: 1, max: 25, group: 'PV module (datasheet, STC)', role: 'material' },
    { key: 'Vmp', label: 'Voltage at maximum power', unit: 'V', type: 'number', default: 41.7, min: 8, max: 100, group: 'PV module (datasheet, STC)', role: 'material', help: 'Default: a 550 Wp monocrystalline module (144 half-cells, 72 in series).' },
    { key: 'Ns', label: 'Cells in series', type: 'int', default: 72, min: 30, max: 160, group: 'PV module (datasheet, STC)', role: 'material' },
    { key: 'gammaP', label: 'Power temperature coefficient', unit: '%/K', type: 'number', default: -0.35, min: -0.6, max: -0.15, group: 'PV module (datasheet, STC)', role: 'material', help: 'Photovoltaic thermal derating: −0.26 (heterojunction) to −0.40 (PERC).' },
    { key: 'betaVoc', label: 'Voc temperature coefficient', unit: '%/K', type: 'number', default: -0.27, min: -0.45, max: -0.15, group: 'PV module (datasheet, STC)', role: 'material' },
    { key: 'modArea', label: 'Module area', unit: 'm²', type: 'number', default: 2.58, min: 0.5, max: 4, group: 'PV module (datasheet, STC)', role: 'geometry' },
    { key: 'tilt', label: 'Module tilt', unit: 'deg', type: 'number', default: 12, min: 0, max: 35, group: 'PV plant', role: 'geometry', design: true, help: 'Floating systems use 5–15° to limit wind load.' },
    { key: 'orient', label: 'Module orientation', type: 'select', default: 'equator', group: 'PV plant', role: 'geometry', options: [{ value: 'equator', label: 'Facing the equator' }, { value: 'custom', label: 'Entered azimuth' }] },
    { key: 'azimuth', label: 'Module azimuth (clockwise from north)', unit: 'deg', type: 'number', default: 180, min: 0, max: 360, group: 'PV plant', role: 'geometry', showIf: { key: 'orient', in: ['custom'] } },
    { key: 'dcAc', label: 'DC/AC ratio', unit: '-', type: 'number', default: 1.2, min: 0.9, max: 1.8, group: 'PV plant', role: 'geometry', design: true, help: 'Electrical operating constraint: the inverter clips above P_dc/ratio.' },
    { key: 'etaInv', label: 'Inverter efficiency', unit: '-', type: 'number', default: 0.975, min: 0.9, max: 0.995, group: 'PV plant', role: 'operation' },
    { key: 'lossDc', label: 'DC losses (mismatch, wiring, wave-induced tilt)', unit: '%', type: 'number', default: 3, min: 0, max: 15, group: 'PV plant', role: 'operation' },
    { key: 'lossAc', label: 'AC losses (transformer, floating and export cable)', unit: '%', type: 'number', default: 2, min: 0, max: 15, group: 'PV plant', role: 'operation' },
    { key: 'soiling', label: 'Soiling loss (salt, birds)', unit: '%', type: 'number', default: 3, min: 0, max: 20, group: 'PV plant', role: 'operation' },
    { key: 'avail', label: 'Availability', unit: '%', type: 'number', default: 98, min: 70, max: 100, group: 'PV plant', role: 'ic', help: 'Operational availability state of the plant.' },
    { key: 'ghi', label: 'Annual global horizontal irradiation', unit: 'kWh/m²', type: 'number', default: 1500, min: 600, max: 2600, group: 'Solar resource', role: 'bc', bind: 'site.ghiAnnual', uncertain: { dist: 'normal', cov: 0.04 }, help: 'Prescribed irradiance boundary condition (long-term mean). North Sea about 1000–1100, Mediterranean 1600–1900, tropics 1700–2100.' },
    { key: 'lat', label: 'Latitude', unit: 'deg', type: 'number', default: 40, min: -60, max: 65, group: 'Solar resource', role: 'environment', bind: 'site.lat' },
    { key: 'Tair', label: 'Annual mean air temperature', unit: '°C', type: 'number', default: 16, min: -5, max: 35, group: 'Solar resource', role: 'bc', bind: 'site.Tair' },
    { key: 'Tamp', label: 'Seasonal amplitude of the monthly temperature', unit: 'K', type: 'number', default: 7, min: 0, max: 20, group: 'Solar resource', role: 'bc' },
    { key: 'Tsea', label: 'Water temperature (annual mean)', unit: '°C', type: 'number', default: 15, min: 0, max: 32, group: 'Solar resource', role: 'bc', bind: 'site.Tsea', help: 'Thermal exchange with the water: a water body cooler than the air lowers the module temperature through the offset below.' },
    { key: 'wind', label: 'Mean wind speed at module height', unit: 'm/s', type: 'number', default: 4.5, min: 0.5, max: 12, group: 'Solar resource', role: 'bc' },
    { key: 'u0', label: 'Faiman constant heat-loss coefficient U0 (floating)', unit: 'W/m²K', type: 'number', default: 31, min: 15, max: 60, group: 'Thermal model', role: 'bc', uncertain: { dist: 'normal', cov: 0.15 }, help: 'Floating-specific input: published fits range from about 25 (closed membrane, small footprint) to 45 W/m²K (open, well-ventilated floats); 25 is the free-standing land value.' },
    { key: 'u1', label: 'Faiman wind heat-loss coefficient U1', unit: 'W·s/m³K', type: 'number', default: 6.84, min: 0, max: 15, group: 'Thermal model', role: 'bc' },
    { key: 'coolFrac', label: 'Share of the air–water temperature difference felt by the modules', unit: '-', type: 'number', default: 0.3, min: 0, max: 1, group: 'Thermal model', role: 'bc', help: 'Offset = share × (water − air temperature), applied when the water is cooler. 0 for air-cooled floats, towards 1 for modules in contact with water.' },
    { key: 'floatVol', label: 'Float volume per module', unit: 'm³', type: 'number', default: 0.22, min: 0.02, max: 5, group: 'Float platform', role: 'geometry', design: true },
    { key: 'floatH', label: 'Float height', unit: 'm', type: 'number', default: 0.35, min: 0.1, max: 3, group: 'Float platform', role: 'geometry' },
    { key: 'deadMass', label: 'Permanent mass per module (float, module, frame, cables)', unit: 'kg', type: 'number', default: 58, min: 10, max: 2000, group: 'Float platform', role: 'load' },
    { key: 'liveMass', label: 'Maintenance live load per module position', unit: 'kg', type: 'number', default: 40, min: 0, max: 500, group: 'Float platform', role: 'load', help: 'Personnel and tools shared over neighbouring floats.' },
    { key: 'fbMin', label: 'Required freeboard', unit: 'm', type: 'number', default: 0.1, min: 0, max: 1, group: 'Float platform', role: 'criteria', help: `Project criterion (DNV-RP-0584 gives guidance). ${CONFIRM}` },
    { key: 'islandMW', label: 'Capacity of one float island', unit: 'MWp', type: 'number', default: 5, min: 0.05, max: 100, group: 'Float platform', role: 'geometry' },
    { key: 'packing', label: 'Module area / island deck area', unit: '-', type: 'number', default: 0.6, min: 0.2, max: 0.95, group: 'Float platform', role: 'geometry', help: 'Ground-cover ratio including walkways.' },
    { key: 'kg', label: 'Centre of gravity above the float bottom', unit: 'm', type: 'number', default: 0.45, min: 0.05, max: 3, group: 'Float platform', role: 'geometry' },
    { key: 'Vext', label: 'Extreme wind speed at panel height (gust)', unit: 'm/s', type: 'number', default: 35, min: 10, max: 80, group: 'Environmental loads', role: 'bc', uncertain: { dist: 'normal', cov: 0.08 } },
    { key: 'cN', label: 'Net normal-force coefficient of the panels', unit: '-', type: 'number', default: 1.0, min: 0.3, max: 2.5, group: 'Environmental loads', role: 'bc', help: `EN 1991-1-4 §7.3 for monopitch canopies at this tilt, or wind-tunnel data for the float system. ${CONFIRM}` },
    { key: 'shelter', label: 'Shelter factor of the inner rows', unit: '-', type: 'number', default: 0.25, min: 0.05, max: 1, group: 'Environmental loads', role: 'bc', help: 'Load on inner rows relative to the windward row (wind-tunnel value; 0.15–0.4 for dense low-tilt arrays).' },
    { key: 'Hs', label: 'Extreme significant wave height at the site', unit: 'm', type: 'number', default: 1.0, min: 0, max: 6, group: 'Environmental loads', role: 'bc', bind: 'site.Hs100', help: 'Sheltered inland or nearshore water: < 1 m; offshore floating PV needs dedicated structures.' },
    { key: 'Tp', label: 'Associated peak period', unit: 's', type: 'number', default: 4, min: 1.5, max: 14, group: 'Environmental loads', role: 'bc' },
    { key: 'Uc', label: 'Current speed', unit: 'm/s', type: 'number', default: 0.5, min: 0, max: 3, group: 'Environmental loads', role: 'bc', bind: 'site.Uc100' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 15, min: 1.5, max: 300, group: 'Environmental loads', role: 'environment', bind: 'site.depth' },
    { key: 'refl', label: 'Wave reflection coefficient of the island', unit: '-', type: 'number', default: 0.4, min: 0, max: 1, group: 'Environmental loads', role: 'bc' },
    { key: 'nLines', label: 'Mooring lines per island', type: 'int', default: 48, min: 4, max: 200, group: 'Mooring', role: 'bc', design: true, help: 'Equally spaced around the perimeter; the island is restrained elastically.' },
    { key: 'kLine', label: 'Horizontal stiffness per line', unit: 'kN/m', type: 'number', default: 60, min: 1, max: 5000, group: 'Mooring', role: 'bc' },
    { key: 'pre', label: 'Pretension per line', unit: 'kN', type: 'number', default: 8, min: 0, max: 500, group: 'Mooring', role: 'ic', help: 'Initial mooring equilibrium of the island.' },
    { key: 'lineCap', label: 'Design capacity of one line and anchor', unit: 'kN', type: 'number', default: 180, min: 5, max: 5000, group: 'Mooring', role: 'criteria', design: true },
    { key: 'gammaF', label: 'Load factor on environmental loads', unit: '-', type: 'number', default: 1.35, min: 1, max: 2, group: 'Mooring', role: 'criteria', help: CONFIRM },
    { key: 'hingeLimit', label: 'Allowable connector rotation', unit: 'deg', type: 'number', default: 12, min: 2, max: 45, group: 'Float platform', role: 'criteria', help: 'From the float supplier: fatigue of connector pins and ears is the usual failure mode.' },
  ],
  compute(inp) {
    const pdc0 = inp.pdcMW * 1e6;
    need(inp.Vmp < inp.Voc && inp.Imp < inp.Isc, 'Datasheet values are inconsistent: Vmp must be below Voc and Imp below Isc.');
    // ---- module: five-parameter single-diode model from the datasheet ----
    const sd = run(() => extractSingleDiode({ Isc: inp.Isc, Voc: inp.Voc, Imp: inp.Imp, Vmp: inp.Vmp, Ns: inp.Ns, alphaIsc: 0.0005 * inp.Isc, betaVoc: inp.betaVoc / 100 * inp.Voc }), 'Single-diode parameter extraction (the datasheet values admit no physical solution — check Isc, Voc, Imp, Vmp and cells in series)');
    const stc = singleDiodeCurve(sd, 60), pModule = stc.pmp, nMod = Math.ceil(pdc0 / pModule), etaStc = pModule / (1000 * inp.modArea);
    // ---- resource: monthly insolation from the annual total with the seasonality of extraterrestrial radiation ----
    const h0 = h0Monthly(inp.lat), wsum = h0.reduce((a, h, i) => a + h * DAYS[i], 0);
    need(wsum > 0, 'No sun at this latitude.');
    const monthlyGhi = h0.map((h) => inp.ghi * 1e3 * h / wsum);
    const north = inp.lat >= 0, monthlyTemp = MONTHS.map((_, i) => inp.Tair + inp.Tamp * Math.cos(2 * Math.PI * (i - (north ? 6.5 : 0.5)) / 12));
    const yr = run(() => syntheticYear({ monthlyGhi, monthlyTempC: monthlyTemp, diurnalRangeC: 6, lat: inp.lat, lon: 0 }), 'Synthetic irradiance year');
    const offset = inp.coolFrac * Math.min(0, inp.Tsea - inp.Tair);
    const base = { times: yr.times, ghi: yr.ghi, tempAirC: yr.tempAirC, windSpeed: inp.wind, lat: inp.lat, lon: 0, tiltDeg: inp.tilt, surfaceAzimuthDeg: inp.orient === 'custom' ? inp.azimuth : (inp.lat >= 0 ? 180 : 0), pdc0, albedo: 0.08, soiling: inp.soiling / 100,
      gammaPdc: inp.gammaP / 100, dcLoss: inp.lossDc / 100, dcAcRatio: inp.dcAc, inverterEfficiency: inp.etaInv, acLoss: inp.lossAc / 100, availability: inp.avail / 100, degradationPerYear: 0 };
    const sim = run(() => simulatePv({ ...base, thermal: 'faiman', thermalParams: { u0: inp.u0, u1: inp.u1, waterCoolingOffsetC: offset } }), 'PV energy simulation');
    const land = simulatePv({ ...base, thermal: 'faiman', thermalParams: { u0: 25, u1: 6.84 } });
    need(sim.energyWh > 0, 'The plant produces no energy with these inputs.');
    // irradiance-weighted cell temperatures
    let sw = 0, st = 0, sl = 0, tMax = -99, pMax = 0, iMax = 0;
    for (let i = 0; i < sim.poa.length; i++) { const w = sim.poa[i]; sw += w; st += w * sim.tempCellC[i]; sl += w * land.tempCellC[i]; if (sim.tempCellC[i] > tMax) tMax = sim.tempCellC[i]; if (sim.pac[i] > pMax) { pMax = sim.pac[i]; iMax = i; } }
    const tCell = st / sw, tLand = sl / sw, gainCool = coolingYieldGain({ deltaTempC: tCell - tLand, gammaPdc: inp.gammaP / 100 });
    // ---- float platform ----
    const buoy = floatBuoyancy({ floatVolume: inp.floatVol, floatHeight: inp.floatH, deadMass: inp.deadMass, liveMass: inp.liveMass, rho: 1025, minFreeboard: inp.fbMin });
    const cDead = buoy.cases[0], cLive = buoy.cases[1];
    const nIsl = Math.max(1, Math.round(inp.pdcMW / inp.islandMW)), modIsl = Math.ceil(nMod / nIsl);
    const deck = modIsl * inp.modArea / inp.packing, side = Math.sqrt(deck), massIsl = modIsl * inp.deadMass;
    const fill = Math.min(0.95, modIsl * buoy.waterplaneArea / deck);
    const rows = Math.max(1, Math.round(side / (inp.modArea ** 0.5 * 1.3 / Math.sqrt(inp.packing)))), rowArea = modIsl * inp.modArea / rows;
    const windL = panelWindLoad({ windSpeed: inp.Vext, panelArea: rowArea, tiltDeg: Math.max(inp.tilt, 1), cN: inp.cN, rows, shelterFactor: inp.shelter });
    // stability of an edge strip of the island (four float rows wide) under the windward-row load and one person at the edge
    const stripW = Math.min(side, 4 * Math.sqrt(inp.modArea / inp.packing)), stripMass = massIsl * stripW / side;
    const stab = islandStability({ mass: stripMass, kg: inp.kg, draft: cDead.draft, length: side, width: stripW, fillRatio: fill, windForce: windL.dragFirstRow, windLever: inp.floatH + 0.5, personMass: 100, personOffset: stripW / 2, rho: 1025 });
    const wl = dispersion(2 * Math.PI / inp.Tp, inp.depth).wavelength;
    const Hmax = Math.min(1.86 * inp.Hs, 0.78 * inp.depth);
    const subVol = massIsl / 1025;
    const wave = inp.Hs > 0 ? floatWaveLoad({ H: Hmax, T: inp.Tp, depth: inp.depth, submergedVolume: subVol * Math.min(1, wl / (2 * side)), projectedArea: side * cDead.draft, rho: 1025 }).max : 0;
    const drift = inp.Hs > 0 ? meanDriftForce({ Hs: inp.Hs, width: side, reflection: inp.refl, rho: 1025 }) : 0;
    const cur = currentDrag({ speed: inp.Uc, area: side * cDead.draft + 0.004 * deck, Cd: 1.0, rho: 1025 });
    const Fenv = windL.drag + wave + drift + cur, Fd = Fenv * inp.gammaF;
    const lines = Array.from({ length: inp.nLines }, (_, i) => ({ azimuthDeg: 360 * i / inp.nLines, stiffness: inp.kLine * 1e3, pretension: inp.pre * 1e3 }));
    const moor = mooringLoadSharing({ Fx: Fd, Fy: 0, lines });
    const anc = anchorCount({ load: Fd, anchorCapacity: inp.lineCap * 1e3, sharing: Math.max(0.05, Math.min(1, moor.sharing * moor.loadedLines / inp.nLines)) });
    const hinge = inp.Hs > 0 ? hingedRaftResponse({ H: Hmax, wavelength: wl, moduleLength: Math.sqrt(inp.modArea / inp.packing) * 1.3, moduleWidth: Math.sqrt(inp.modArea / inp.packing), fillRatio: fill, rho: 1025 }) : { relativeRotationDeg: 0, connectorShear: 0, regime: 'contouring', lengthRatio: Infinity };
    const expo = exposureClass({ Hs: inp.Hs });
    const checks = [
      check('Freeboard with maintenance load', inp.fbMin, Math.max(cLive.freeboard, 1e-6), 'm', 'DNV-RP-0584 (2021) §5 (buoyancy and freeboard); Archimedes', 'Demand is the required freeboard, capacity the freeboard with dead and live load. Too little freeboard lets waves wash over the floats and wets connectors and cables.'),
      check('Buoyancy utilisation (weight / total buoyancy)', cLive.utilisation, 0.7, '-', 'DNV-RP-0584 (2021) §5 (reserve buoyancy); project criterion 70 %', 'Share of the float buoyancy used by dead plus live load. Above about 0.7 there is little reserve for snow, marine growth, water ingress or a punctured float.'),
      check('Heel of an edge strip under wind and a person at the edge', Math.abs(stab.heelDeg), 10, 'deg', 'DNV-RP-0584 (2021) §5 (stability); small-angle hydrostatics', 'Heel of a four-row edge strip against 10°. Larger heel submerges the float edge, makes walkways unsafe and changes the wind load on the panels.'),
      check('Maximum mooring-line tension (factored loads)', moor.maxTension / 1e3, inp.lineCap, 'kN', 'DNV-RP-0584 (2021) §7 (station keeping); DNV-OS-E301 principles', 'Tension of the most loaded line under factored wind, wave, drift and current load against the design capacity of line and anchor. An overloaded line fails first and sheds its load to the neighbours (unzipping).'),
      check('Connector rotation in the extreme wave', hinge.relativeRotationDeg, inp.hingeLimit, 'deg', 'DNV-RP-0584 (2021) §6 (float connections); quasi-static hinged-raft kinematics', 'Relative rotation between neighbouring floats when the island follows the extreme wave against the allowable of the connector. Repeated large rotations fatigue the connector ears and pins.'),
    ];
    const gov = worst(checks);
    const polymerT = nMod * inp.floatVol * 60 / 1e3, alT = nMod * 4 / 1e3;   // ~60 kg HDPE per m³ of float volume, 4 kg aluminium frame per module (indicative)
    const co2 = polymerT * CO2_POLYMER + alT * CO2_ALUMINIUM;
    const mwh = sim.energyWh / 1e6;
    // I–V curves at four irradiances and the operating cell temperature
    const ivs = [1000, 800, 600, 200].map((Gi) => { const p = singleDiodeAt(sd, { G: Gi, tempCellC: cellTempFaiman({ poa: Gi, tempAirC: inp.Tair, windSpeed: inp.wind, u0: inp.u0, u1: inp.u1 }) }); return { G: Gi, c: singleDiodeCurve(p, 40) }; });
    const d0 = Math.floor(iMax / 24) * 24, hrs = linspace(0, 23, 24);
    return {
      kpis: [
        { key: 'aep', label: 'Annual energy (AC, year 1)', value: round(sim.energyWh / 1e9, 5), unit: 'GWh/yr', explain: 'Energy delivered at the grid connection in the first year, after all losses and availability, from the long-term mean irradiation (a P50 value: exceeded with 50 % probability in a given year).' },
        { key: 'yield', label: 'Specific yield', value: round(sim.specificYield, 5), unit: 'kWh/kWp', explain: 'Annual AC energy per installed kWp. It allows comparison between sites and layouts independent of plant size.' },
        { key: 'pr', label: 'Performance ratio', value: round(sim.performanceRatio * 100, 4), unit: '%', explain: 'Actual yield over the yield of a loss-free plant at standard conditions under the same in-plane irradiation (IEC 61724-1). 78–86 % is normal; floating plants gain from cooling and lose on soiling and mismatch.' },
        { key: 'cf', label: 'Capacity factor (AC energy / DC rating)', value: round(sim.capacityFactor * 100, 4), unit: '%', explain: 'Mean output over the DC rating: 11–20 % depending on latitude.' },
        { key: 'ghiYear', label: 'Annual global horizontal irradiation of the simulated year', value: round(sim.ghiWhPerM2 / 1e3, 5), unit: 'kWh/m²', explain: 'Sum of the hourly irradiance used in the simulation: it must reproduce the entered annual value (resource conservation check).' },
        { key: 'poaYear', label: 'Annual in-plane irradiation', value: round(sim.poaWhPerM2 / 1e3, 5), unit: 'kWh/m²', explain: 'Irradiation on the tilted module plane after transposition. The gain over the horizontal value is small at the low tilts used on floats.' },
        { key: 'pmpStc', label: 'Module power at STC from the single-diode model', value: round(pModule, 5), unit: 'W', explain: 'Maximum power of the five-parameter model at 1000 W/m² and 25 °C: it must equal Vmp × Imp of the datasheet.' },
        { key: 'etaStc', label: 'Module efficiency at STC', value: round(etaStc * 100, 4), unit: '%', explain: 'Module power over irradiance times module area: the photovoltaic conversion efficiency.' },
        { key: 'tCell', label: 'Irradiance-weighted cell temperature', value: round(tCell, 4), unit: '°C', explain: `Mean cell temperature weighted by irradiance (maximum ${fmt(tMax)} °C). The same plant on land with the reference coefficients would run at ${fmt(tLand)} °C.` },
        { key: 'gainCool', label: 'Yield gain from cooling over water', value: round(gainCool * 100, 4), unit: '%', explain: 'Energy gain relative to a land-based plant from the lower cell temperature: temperature coefficient × temperature difference. Measured gains are typically 1–5 %.' },
        { key: 'clipped', label: 'Energy clipped by the inverters', value: round(sim.clippedWh / Math.max(sim.energyWh + sim.clippedWh, 1) * 100, 4), unit: '%', explain: 'Energy lost when the DC power exceeds the inverter rating. A few percent is the economic optimum of the DC/AC ratio.' },
        { key: 'draft', label: 'Float draft under permanent load', value: round(cDead.draft, 5), unit: 'm', explain: 'Immersion of the float under its permanent load, m/(ρ·A_wp).' },
        { key: 'freeboard', label: 'Freeboard with maintenance load', value: round(cLive.freeboard, 5), unit: 'm', limit: inp.fbMin, sense: 'min', explain: 'Float height minus draft with a person on the float.' },
        { key: 'GM', label: 'Metacentric height of the edge strip', value: round(Math.min(stab.GM, 999), 4), unit: 'm', explain: 'Initial stability of a four-row edge strip of the island: large for wide rafts, the critical case is the narrow strip at the perimeter.' },
        { key: 'Fenv', label: 'Extreme horizontal load on one island', value: round(Fenv / 1e3, 5), unit: 'kN', explain: 'Characteristic wind drag on the panels plus wave, drift and current forces on the floats, acting together.' },
        { key: 'Tmoor', label: 'Maximum mooring-line tension (factored)', value: round(moor.maxTension / 1e3, 5), unit: 'kN', limit: inp.lineCap, sense: 'max', explain: 'Tension in the most loaded line; lines on the leeward side go slack and carry nothing.' },
        { key: 'offset', label: 'Island offset under the factored load', value: round(moor.offset, 4), unit: 'm', explain: 'Horizontal excursion of the island; the floating power cable must have this much slack.' },
        { key: 'hingeRot', label: 'Connector rotation in the extreme wave', value: round(hinge.relativeRotationDeg, 4), unit: 'deg', limit: inp.hingeLimit, sense: 'max', explain: `Relative rotation between adjacent floats (wavelength ${fmt(wl)} m, regime: ${hinge.regime}).` },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Monthly energy', categories: MONTHS, series: [{ name: 'AC energy (MWh)', y: sim.monthly.map((m) => round(m.energyWh / 1e6, 4)) }], y: { label: 'Energy', unit: 'MWh' },
          explain: 'AC energy per calendar month. Solar output peaks in summer and complements winter-peaking wind on a shared grid connection.' },
        { type: 'line', title: 'Power on the day of maximum output', x: { label: 'Hour (UTC, longitude 0°)', unit: 'h' }, y: { label: 'Power (MW) · irradiance (kW/m² × 10)', unit: '-' }, series: [{ name: 'AC power (MW)', x: hrs, y: hrs.map((h) => round((sim.pac[d0 + h] ?? 0) / 1e6, 4)) }, { name: 'In-plane irradiance × 10 (kW/m²)', x: hrs, y: hrs.map((h) => round((sim.poa[d0 + h] ?? 0) / 100, 4)) }, { name: 'Inverter limit (MW)', x: [0, 23], y: [pdc0 / inp.dcAc / 1e6, pdc0 / inp.dcAc / 1e6], dash: true }],
          explain: 'Instantaneous output over the best day with the in-plane irradiance. A flat top at the dashed line is inverter clipping.' },
        { type: 'line', title: 'Module I–V curves at operating temperature', x: { label: 'Voltage', unit: 'V' }, y: { label: 'Current', unit: 'A' }, series: ivs.map((q) => ({ name: `${q.G} W/m² (Pmp ${fmt(q.c.pmp)} W)`, x: q.c.v.map((v) => round(v, 3)), y: q.c.i.map((v) => round(v, 4)) })),
          explain: 'Current–voltage characteristics from the single-diode model at four irradiance levels and the corresponding cell temperatures. The knee of each curve is the maximum-power point the inverter tracks.' },
        { type: 'bar', title: 'Environmental loads on one island (characteristic)', categories: ['Wind on panels', 'Wave (Morison)', 'Mean wave drift', 'Current'], series: [{ name: 'Force (kN)', y: [windL.drag, wave, drift, cur].map((v) => round(v / 1e3, 4)) }], y: { label: 'Force', unit: 'kN' },
          explain: 'Contribution of each environmental action to the horizontal load on the mooring. Wind on the panel rows normally dominates on sheltered water.' },
        { type: 'polar', title: 'Mooring-line tensions around the island', series: [{ name: 'Tension (kN)', thetaDeg: lines.map((l) => round(l.azimuthDeg, 2)), r: moor.tensions.map((v) => round(v / 1e3, 4)) }],
          explain: 'Tension of each line by azimuth under the factored load (acting towards 90°). Windward lines carry the load; zero means a slack line.' },
      ],
      tables: [
        { title: 'Monthly irradiation and energy', columns: [{ key: 'm', label: 'Month' }, { key: 'ghi', label: 'GHI', unit: 'kWh/m²' }, { key: 'poa', label: 'In-plane', unit: 'kWh/m²' }, { key: 'T', label: 'Air temperature', unit: '°C' }, { key: 'e', label: 'AC energy', unit: 'MWh' }, { key: 'd', label: 'Daily mean', unit: 'MWh/day' }, { key: 'pr', label: 'PR', unit: '%' }],
          rows: MONTHS.map((m, i) => ({ m, ghi: round(monthlyGhi[i] * DAYS[i] / 1e3, 4), poa: round(sim.monthly[i].poaWhPerM2 / 1e3, 4), T: round(monthlyTemp[i], 3), e: round(sim.monthly[i].energyWh / 1e6, 4), d: round(sim.monthly[i].energyWh / DAYS[i] / 1e6, 4), pr: round(sim.monthly[i].poaWhPerM2 > 0 ? sim.monthly[i].energyWh / (pdc0 * sim.monthly[i].poaWhPerM2 / 1000) * 100 : 0, 3) })),
          explain: 'Horizontal and in-plane irradiation, temperature, energy, mean daily energy and performance ratio for each month. The summer drop in PR is the thermal derating.' },
        { title: 'Float island', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [{ q: 'Modules (plant / island)', v: `${nMod} / ${modIsl}`, u: '-' }, { q: 'Islands', v: nIsl, u: '-' }, { q: 'Island deck area and side', v: `${fmt(deck, 4)} / ${fmt(side)}`, u: 'm² / m' }, { q: 'Draft dead / with live load', v: `${fmt(cDead.draft)} / ${fmt(cLive.draft)}`, u: 'm' }, { q: 'Reserve buoyancy (dead load)', v: round(cDead.reserveBuoyancy, 4), u: '-' }, { q: 'Wind drag, windward row / all rows', v: `${fmt(windL.dragFirstRow / 1e3)} / ${fmt(windL.drag / 1e3)}`, u: 'kN' },
            { q: 'Wind uplift on the island', v: round(windL.lift / 1e3, 4), u: 'kN' }, { q: 'Loaded / slack lines', v: `${moor.loadedLines} / ${moor.slackLines}`, u: '-' }, { q: 'Anchors needed for the design load', v: anc.count, u: '-' }, { q: 'Connector shear in the extreme wave', v: round(hinge.connectorShear / 1e3, 4), u: 'kN' }, { q: 'Exposure class', v: expo.class, u: '' }],
          explain: 'Geometry, buoyancy, wind load and mooring particulars of one float island.' },
      ],
      notes: [
        `${fmt(inp.pdcMW)} MWp with ${nMod} modules of ${fmt(pModule)} W (single-diode: n = ${fmt(sd.n)}, Rs = ${fmt(sd.Rs)} Ω, Rsh = ${fmt(sd.Rsh)} Ω); tilt ${inp.tilt}°, DC/AC ${inp.dcAc}.`,
        `Yield ${fmt(sim.specificYield, 4)} kWh/kWp at ${inp.ghi} kWh/m² and latitude ${inp.lat}°: PR ${fmt(sim.performanceRatio * 100)} %, cooling gain ${fmt(gainCool * 100)} % against a land plant (U0 = ${inp.u0} against 25 W/m²K, water offset ${fmt(offset)} K). Annual energy is a P50 value: exceeded with 50 % probability.`,
        `Island of ${fmt(inp.islandMW)} MWp: ${fmt(side)} m square, ${fmt(massIsl / 1e3)} t; extreme load ${fmt(Fenv / 1e3)} kN (wind ${fmt(windL.drag / 1e3)} kN), ${inp.nLines} lines, maximum factored tension ${fmt(moor.maxTension / 1e3)} kN.`,
        `Governing check: ${gov.name} (utilisation ${fmt(gov.uc)}). Site exposure class: ${expo.class}.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `The plant yields ${fmt(sim.energyWh / 1e9, 4)} GWh/yr (${fmt(sim.specificYield, 4)} kWh/kWp, PR ${fmt(sim.performanceRatio * 100)} %); freeboard ${fmt(cLive.freeboard)} m, maximum mooring tension ${fmt(moor.maxTension / 1e3)} kN against ${inp.lineCap} kN, connector rotation ${fmt(hinge.relativeRotationDeg)}°.`,
          gov.uc > 1 ? (gov === checks[3] ? `Add mooring lines (about ${Math.ceil(inp.nLines * gov.uc * 1.1)}) or raise the line and anchor capacity to ${fmt(moor.maxTension / 1e3 * 1.05)} kN; softer lines share load more evenly.` : gov === checks[0] || gov === checks[1] ? `Increase the float volume to about ${fmt(inp.floatVol * Math.max(gov.uc, 1.05) * 1.05)} m³ per module or reduce the supported mass.` : gov === checks[4] ? 'The site is too exposed for this float type: choose larger modules (longer than half the wavelength), a membrane or a semi-submersible concept, or a more sheltered location.' : 'Widen the float arrangement at the perimeter or lower the centre of gravity to reduce heel.') : 'The float platform meets the screening criteria. Confirm wind coefficients by wind-tunnel tests of the float system and the cell-temperature coefficients from an operating plant or the supplier.',
          'DNV-RP-0584 (2021) for the floating system; IEC 61853 and IEC 61724-1 for module rating and performance ratio; yield modelling per IEC TS 62738.',
          `Floats and frames: about ${fmt(polymerT)} t of HDPE and ${fmt(alT)} t of aluminium ≈ ${fmt(co2)} t CO₂e (screening factors ${CO2_POLYMER} and ${CO2_ALUMINIUM} t CO₂e/t); the plant covers ${fmt(nIsl * deck / 1e4)} ha of water, which reduces evaporation and light for aquatic life — keep coverage of a water body moderate. Yield displaces about ${fmt(mwh * CO2_GRID)} t CO₂ a year.`),
        adv('economic', 'info', `Specific yield ${fmt(sim.specificYield, 4)} kWh/kWp; the cooling gain is worth ${fmt(mwh * gainCool / (1 + gainCool))} MWh a year and inverter clipping costs ${fmt(sim.clippedWh / 1e6)} MWh.`,
          `Optimise the DC/AC ratio (now ${inp.dcAc}) against clipping and tilt against wind load: a lower tilt cuts mooring and float cost for a small energy loss. One percent of soiling is ${fmt(mwh * 0.01)} MWh a year — the cleaning interval follows from that.`,
          'LCOE conventions (IEA PVPS / IRENA): specific yield and performance ratio in the denominator, float and mooring CAPEX per kWp in the numerator; yield reporting at P50 with P90 for financing.',
          energyCo2Text(mwh)),
      ],
      publish: { 'fpv.ratedPower': pdc0, 'fpv.aep': sim.energyWh, 'fpv.capacityFactor': sim.capacityFactor, 'fpv.specificYield': sim.specificYield, 'fpv.mooringLoad': Fenv },
      quantities: { polymer_t: round(polymerT, 4), aluminium_t: round(alT, 4) },
    };
  },
  verify: [
    { name: 'Single-diode model reproduces the datasheet maximum power: Pmp = Vmp × Imp', source: 'De Soto, Klein & Beckman (2006) Solar Energy 80: five-parameter model fitted to the datasheet points',
      // 41.7 V × 13.2 A = 550.44 W
      inputs: { Isc: 13.9, Voc: 49.6, Imp: 13.2, Vmp: 41.7, Ns: 72 }, expect: [{ kpi: 'pmpStc', value: 550.44, tol: 0.002 }] },
    { name: 'Float draft = m / (ρ·V/h) and conservation of the annual irradiation', source: 'Archimedes principle (hand calculation); synthetic year built from monthly insolation',
      // V = 0.12 m³, h = 0.30 m → A_wp = 0.40 m²; m = 60 kg, ρ = 1025 kg/m³ → draft = 60/(1025 × 0.40) = 0.14634 m; annual GHI must return 1500 kWh/m²
      inputs: { floatVol: 0.12, floatH: 0.3, deadMass: 60, ghi: 1500 }, expect: [{ kpi: 'draft', value: 0.14634, tol: 0.001 }, { kpi: 'ghiYear', value: 1500, tol: 0.01 }] },
  ],
};
