// Thermal expansion and temperature response: steady core and sheath temperatures from conductor losses, cooling
// of a warm service fluid along a tube, differential thermal strain and stress between the warm tube and the rest
// of the bundle, axial expansion against seabed friction and the warm-up / cool-down time constants.
import { umbilicalCoreTemperature } from '../../engine/xsection/umbilical.js';
import { internalFilmCoefficient, externalFilmCoefficient, overallU, cooldownLumped } from '../../engine/flow/thermalpipe.js';
import { endExpansionUniform } from '../../engine/pipeline/expansion.js';
import { adv, check, worst, sev, fmt, round, run, need, linspace, G, STD_UMB, CO2_POLYMER } from './_common.js';

export default {
  id: 's12-thermal-response', suite: 12, title: 'Thermal expansion and temperature response',
  summary: 'Core, sheath and fluid temperatures, thermal strain and stress, end expansion or restraint force on the seabed, and warm-up and cool-down times.',
  covers: [12], fidelity: 'analytical',
  standards: ['IEC 60287-1-1 / -2-1 (thermal-resistance method, simplest series form)', `${STD_UMB} §7 (thermal analysis, temperature limits of materials)`, 'Incropera et al., Fundamentals of Heat and Mass Transfer §3.3 (composite cylinder)', 'Palmer & King, Subsea Pipeline Engineering ch. 12 (end expansion against friction)'],
  models: ['Fourier heat conduction through multilayer cylindrical resistances', "Newton's law of cooling at the outer surface", 'Internal convective heat transfer (Gnielinski / laminar)', 'Transient thermal diffusion (lumped capacity)', 'Thermal expansion and coupled thermal stress between elements sharing axial strain', 'Coulomb friction anchoring of axial expansion'],
  limits: 'Series thermal network: all core losses cross the same bundle resistance, which must be estimated or measured; no solar heating in the I-tube, no burial (use the soil resistance input of the power-cable suite for buried sections). The warm tube and the rest of the bundle share one axial strain (no slip); temperature-dependent material properties and a transient radial finite-difference solution are outside this module.',
  publishes: ['umbilical.thermalForce', 'umbilical.endExpansion', 'umbilical.coreTemp'],
  searchTerms: 'umbilical thermal analysis conductor temperature expansion',
  inputs: [
    { key: 'OD', label: 'Overall diameter', unit: 'mm', type: 'number', default: 145, min: 30, max: 400, group: 'Umbilical', role: 'geometry', bind: 'umbilical.OD' },
    { key: 'sheathT', label: 'Outer sheath thickness', unit: 'mm', type: 'number', default: 6, min: 0.5, max: 30, group: 'Umbilical', role: 'geometry' },
    { key: 'sheathMat', label: 'Outer sheath material', type: 'select', source: 'materials:polymer', default: 'hdpe', group: 'Umbilical', role: 'material' },
    { key: 'EA', label: 'Axial stiffness EA', unit: 'MN', type: 'number', default: 720, min: 1, max: 5000, group: 'Umbilical', role: 'material', bind: 'umbilical.EA', uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'wSub', label: 'Submerged weight', unit: 'N/m', type: 'number', default: 235, min: 1, max: 3000, group: 'Umbilical', role: 'geometry', bind: 'umbilical.wSub' },
    { key: 'mass', label: 'Mass in air', unit: 'kg/m', type: 'number', default: 40.9, min: 1, max: 300, group: 'Umbilical', role: 'geometry', bind: 'umbilical.massPerM' },
    { key: 'cpEff', label: 'Effective specific heat of the bundle', unit: 'J/(kg·K)', type: 'number', default: 750, min: 300, max: 2500, group: 'Umbilical', role: 'material', help: 'Mass-weighted: steel ≈ 480, copper ≈ 385, polyethylene ≈ 1900, water ≈ 4000.' },
    { key: 'alpha', label: 'Thermal expansion coefficient of the load-bearing steel', unit: '10⁻⁶/K', type: 'number', default: 13, min: 5, max: 30, group: 'Umbilical', role: 'material', help: 'Super duplex 13, carbon-steel armour 11.7, copper 17.' },
    { key: 'length', label: 'Length on the seabed', unit: 'm', type: 'number', default: 5000, min: 10, max: 3e5, group: 'Umbilical', role: 'geometry', bind: 'cable.length' },
    { key: 'cores', label: 'Electrical cores (heat sources)', type: 'table', group: 'Heat sources', role: 'bc',
      columns: [{ key: 'name', label: 'Core', type: 'text' }, { key: 'loss', label: 'Losses', unit: 'W/m', type: 'number' }, { key: 'dCond', label: 'Conductor diameter', unit: 'mm', type: 'number' }, { key: 'dIns', label: 'Insulation diameter', unit: 'mm', type: 'number' }, { key: 'kIns', label: 'Insulation conductivity', unit: 'W/(m·K)', type: 'number' }],
      default: [{ name: 'Power core L1', loss: 6, dCond: 8, dIns: 16, kIns: 0.286 }, { name: 'Power core L2', loss: 6, dCond: 8, dIns: 16, kIns: 0.286 }, { name: 'Power core L3', loss: 6, dCond: 8, dIns: 16, kIns: 0.286 }, { name: 'Signal quad', loss: 0.3, dCond: 3, dIns: 7, kIns: 0.3 }],
      help: 'I²R (and dielectric) losses per metre of each core at the design current, from the electrical analysis.' },
    { key: 'Rbundle', label: 'Bundle (filler) thermal resistance', unit: 'K·m/W', type: 'number', default: 0.35, min: 0, max: 5, group: 'Heat sources', role: 'material', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Between the cores and the inside of the sheath; depends on lay-up and whether interstices are flooded. Estimate or measure.' },
    { key: 'tubeID', label: 'Warm-fluid tube bore', unit: 'mm', type: 'number', default: 26.6, min: 3, max: 100, group: 'Warm service fluid', role: 'geometry' },
    { key: 'tubeWall', label: 'Tube wall', unit: 'mm', type: 'number', default: 2.6, min: 0.3, max: 20, group: 'Warm service fluid', role: 'geometry' },
    { key: 'tubeMat', label: 'Tube material', type: 'select', source: 'materials:steel', default: 'sdss-25cr', group: 'Warm service fluid', role: 'material' },
    { key: 'fluid', label: 'Service fluid', type: 'select', source: 'materials:fluid', default: 'meg', group: 'Warm service fluid', role: 'material' },
    { key: 'Tin', label: 'Fluid temperature entering the umbilical', unit: '°C', type: 'number', default: 45, min: -20, max: 150, group: 'Warm service fluid', role: 'bc', help: 'Topside supply temperature (hydraulic boundary condition).' },
    { key: 'flow', label: 'Flow rate', unit: 'L/min', type: 'number', default: 15, min: 0, max: 2000, group: 'Warm service fluid', role: 'bc', help: '0 = stagnant line at ambient temperature.' },
    { key: 'kFill', label: 'Conductivity between the tube and the sheath', unit: 'W/(m·K)', type: 'number', default: 0.35, min: 0.05, max: 5, group: 'Warm service fluid', role: 'material', help: 'Smeared fillers and flooded interstices.' },
    { key: 'Tamb', label: 'Ambient sea-water temperature', unit: '°C', type: 'number', default: 4, min: -2, max: 35, group: 'Environment', role: 'bc', bind: 'site.Tbed' },
    { key: 'Tinstall', label: 'Temperature at installation (strain-free state)', unit: '°C', type: 'number', default: 15, min: -10, max: 45, group: 'Environment', role: 'ic', help: 'Initial condition: temperature of the umbilical when laid and locked in by seabed friction.' },
    { key: 'Ubed', label: 'Near-bed current', unit: 'm/s', type: 'number', default: 0.1, min: 0.005, max: 3, group: 'Environment', role: 'bc', bind: 'site.UcBed100', help: 'Sets the external film coefficient; a low value is conservative for temperature.' },
    { key: 'muAx', label: 'Axial seabed friction coefficient', unit: '-', type: 'number', default: 0.4, min: 0.01, max: 1.5, group: 'Restraint', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'ends', label: 'End condition', type: 'select', default: 'free', group: 'Restraint', role: 'bc', options: [{ value: 'free', label: 'Free to slide (termination on a sliding mudmat or slack loop)' }, { value: 'fixed', label: 'Fixed terminations (fully restrained)' }] },
    { key: 'TcoreMax', label: 'Allowable conductor temperature', unit: '°C', type: 'number', default: 90, min: 40, max: 150, group: 'Acceptance criteria', role: 'criteria', help: 'XLPE / EPR 90 °C continuous (IEC 60502-2).' },
    { key: 'TsheathMax', label: 'Allowable sheath / filler temperature', unit: '°C', type: 'number', default: 60, min: 30, max: 130, group: 'Acceptance criteria', role: 'criteria', help: 'Polyethylene sheath and fillers: typically 60–70 °C continuous; hose liners may be lower.' },
    { key: 'eta', label: 'Allowable thermal stress / yield in the tube', unit: '-', type: 'number', default: 0.67, min: 0.3, max: 1, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'termRating', label: 'Termination axial rating', unit: 'kN', type: 'number', default: 150, min: 1, max: 5000, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'Tthreshold', label: 'Minimum acceptable fluid temperature after shut-in', unit: '°C', type: 'number', default: 10, min: -20, max: 100, group: 'Acceptance criteria', role: 'criteria', help: 'E.g. viscosity or wax limit of the service fluid; used for the cool-down time.' },
  ],
  compute(inp, ctx) {
    const D = inp.OD / 1e3, ro = D / 2, ri = ro - inp.sheathT / 1e3;
    need(ri > 0.005, 'Sheath thickness leaves no room for the bundle.');
    const sm = ctx.material(inp.sheathMat), tm = ctx.material(inp.tubeMat), fl = ctx.material(inp.fluid), sw = ctx.material('seawater');
    need(sm?.k, `Sheath material '${inp.sheathMat}' has no thermal conductivity in the library.`);
    need(tm?.E && fl?.cp && fl?.k && fl?.mu, 'Tube material or service fluid lacks the properties needed (E, cp, k, viscosity).');
    const cores = (inp.cores || []).filter((c) => c.loss >= 0 && c.dIns > 0);
    need(cores.length > 0, 'Give at least one electrical core (set its losses to 0 if the umbilical carries no power).');
    for (const c of cores) need(c.dIns > c.dCond && c.dCond > 0 && c.kIns > 0, `Core '${c.name}': insulation diameter must exceed conductor diameter and conductivity must be positive.`);
    // external film in cross flow (Newton's law of cooling at the surface)
    const hExt = run(() => externalFilmCoefficient({ v: inp.Ubed, D, rho: sw.rho, mu: sw.mu, k: sw.k, cp: sw.cp }), 'External film coefficient').h;
    const ct = run(() => umbilicalCoreTemperature({ cores: cores.map((c) => ({ loss: c.loss, conductorDiameter: c.dCond / 1e3, insulationDiameter: c.dIns / 1e3, kInsulation: c.kIns })), ambientC: inp.Tamb, bundleResistance: inp.Rbundle, sheath: { ri, ro, k: sm.k }, hExternal: hExt, outerDiameter: D }), 'Core temperature');
    const TcoreMax = Math.max(...ct.coreTempsC), bundleC = ct.sheathInnerTempC + 0.5 * ct.totalLoss * inp.Rbundle;   // mean bundle temperature from core losses
    // warm fluid tube: overall U from bore film + tube wall + smeared filler + sheath + external film
    const di = inp.tubeID / 1e3, tw = inp.tubeWall / 1e3, At = Math.PI * ((di + 2 * tw) ** 2 - di ** 2) / 4;
    const Q = inp.flow / 60000, v = Q / (Math.PI * di * di / 4), mDot = fl.rho * Q;
    const hi = v > 0 ? run(() => internalFilmCoefficient({ v, D: di, rho: fl.rho, mu: fl.mu, k: fl.k, cp: fl.cp, L: inp.length }), 'Internal film coefficient') : null;
    const tFill = ri - (di / 2 + tw);
    need(tFill > 0, 'The warm tube does not fit inside the sheath.');
    const ou = run(() => overallU({ Di: di, layers: [{ t: tw, k: tm.k ?? 15 }, { t: tFill, k: inp.kFill }, { t: inp.sheathT / 1e3, k: sm.k }], hi: hi ? hi.h : Infinity, ho: hExt * 1 }), 'Overall heat-transfer coefficient');
    // T(x) = Tamb' + (Tin − Tamb')·exp(−UA·x/(ṁ·cp)), Tamb' = bundle temperature from core losses
    const Lth = mDot > 0 ? mDot * fl.cp / ou.UA : 0;
    const Tf = (x) => (mDot > 0 ? bundleC + (inp.Tin - bundleC) * Math.exp(-x / Lth) : bundleC);
    const L = inp.length, Tout = Tf(L), TfMean = mDot > 0 ? bundleC + (inp.Tin - bundleC) * Lth / L * (1 - Math.exp(-L / Lth)) : bundleC;
    // shared axial strain of the warm tube (EA_h) and the rest of the bundle (EA_r): ε = Σ EA·α·ΔT / Σ EA
    const al = inp.alpha * 1e-6, EA = inp.EA * 1e6, EAh = Math.min(tm.E * At, 0.9 * EA), EAr = EA - EAh;
    const dTr = bundleC - inp.Tinstall, dThIn = inp.Tin - inp.Tinstall, dThMean = TfMean - inp.Tinstall;
    const epsAt = (dTh) => (EAh * al * dTh + EAr * al * dTr) / EA;
    const epsIn = epsAt(mDot > 0 ? dThIn : dTr), epsMean = epsAt(mDot > 0 ? dThMean : dTr);
    const sigTube = tm.E * (epsIn - al * (mDot > 0 ? dThIn : dTr));       // + tension, − compression (free bundle, warm inlet)
    const sigTubeRestr = -tm.E * al * (mDot > 0 ? dThIn : dTr);            // tube stress when the whole bundle is held at zero strain
    // axial response against seabed friction
    const P = EA * Math.abs(epsMean), f = inp.muAx * inp.wSub;
    const ex = P > 0 ? endExpansionUniform({ P, f, EA, L }) : { anchorLength: 0, expansion: 0, regime: 'long' };
    const fixed = inp.ends === 'fixed';
    const endMove = fixed ? 0 : ex.expansion * Math.sign(epsMean), termForce = fixed ? P : 0, Nmid = ex.regime === 'long' || fixed ? P : f * L / 2;
    // lumped warm-up / cool-down
    const C = inp.mass * inp.cpEff, Rsurf = Math.log(ro / ri) / (2 * Math.PI * sm.k) + 1 / (Math.PI * D * hExt) + 0.5 * inp.Rbundle, tau = C * Rsurf;
    const cd = inp.Tin > inp.Tthreshold && inp.Tthreshold > inp.Tamb ? cooldownLumped({ T0: inp.Tin + 273.15, Tamb: inp.Tamb + 273.15, U: ou.U, D: ou.Dref, heatCapacity: C, Tthreshold: inp.Tthreshold + 273.15 }) : null;
    const tCool = cd && Number.isFinite(cd.timeTo) ? cd.timeTo / 3600 : 0;
    const sigUse = fixed ? sigTubeRestr : sigTube, fy = tm.fy ?? 550e6;
    const checks = [
      check('Conductor temperature', TcoreMax, inp.TcoreMax, '°C', 'IEC 60287-1-1 (permissible conductor temperature); IEC 60502-2', 'Hottest conductor against the continuous rating of its insulation. Exceeding it ages the insulation exponentially (roughly halving life per 8–10 K).'),
      check('Sheath / filler temperature', Math.max(ct.sheathInnerTempC, mDot > 0 ? inp.Tin : -273), inp.TsheathMax, '°C', `${STD_UMB} §7 (temperature limits of polymer components)`, 'Highest temperature seen by polymer sheath and fillers — from core losses or from the warm tube at the inlet. Above the limit the polymers creep and soften, loosening the lay-up.'),
      check('Thermal stress in the warm tube', Math.abs(sigUse) / 1e6, inp.eta * fy / 1e6, 'MPa', `${STD_UMB} §7.2 (stress utilisation); thermal stress σ = E(ε − αΔT)`, 'Axial stress locked into the warm tube because the cooler bundle (or the end restraint) prevents its free expansion. It adds to pressure and tension stress and cycles with every start-up.'),
      check('Thermal force on the termination', termForce / 1e3, inp.termRating, 'kN', 'Project termination rating; restrained force EA·αΔT', 'Axial force on fixed terminations from restrained expansion or contraction (zero when the ends are free to slide). Exceeding the rating overloads the connector or foundation.'),
    ];
    const gov = worst(checks);
    const xs = linspace(0, L, 60), ts = linspace(0, 5 * tau, 50);
    const rS = [0, 0.4 * ri, ri, ro, ro * 1.15], Trad = [TcoreMax, bundleC, ct.sheathInnerTempC, ct.surfaceTempC, inp.Tamb];
    const strainX = xs.map((x) => epsAt(mDot > 0 ? Tf(x) - inp.Tinstall : dTr) * 100);
    const lossTotal = ct.totalLoss * L / 1e3, energyYr = lossTotal * 8766 / 1e3;   // kW, MWh/yr
    return {
      kpis: [
        { key: 'coreTemp', label: 'Hottest conductor temperature', value: round(TcoreMax), unit: '°C', limit: inp.TcoreMax, sense: 'max', explain: 'Steady conductor temperature of the most heavily loaded core: ambient plus the temperature rises across the external film, sheath, bundle and its own insulation.' },
        { key: 'sheathTemp', label: 'Sheath inner temperature', value: round(ct.sheathInnerTempC), unit: '°C', explain: 'Temperature at the inside of the outer sheath from core losses; fillers and the sheath must tolerate it continuously.' },
        { key: 'surfaceTemp', label: 'Surface temperature', value: round(ct.surfaceTempC), unit: '°C', explain: 'Outer-surface temperature above ambient; a small difference shows that external cooling is effective and the internal resistances govern.' },
        { key: 'riseInternal', label: 'Conductor-to-surface temperature rise', value: round(TcoreMax - ct.surfaceTempC, 5), unit: 'K', explain: 'Temperature drop across insulation, bundle and sheath for the hottest core. It depends only on the losses and the internal resistances, not on the sea conditions.' },
        { key: 'loss', label: 'Total heat loss', value: round(ct.totalLoss), unit: 'W/m', explain: 'Sum of core losses per metre, all of which leaves through the sheath to the sea.' },
        { key: 'hExt', label: 'External film coefficient', value: round(hExt), unit: 'W/(m²·K)', explain: 'Convective coefficient on the outer surface in the near-bed current (Newton cooling). Low current gives a lower value and a warmer umbilical.' },
        { key: 'UA', label: 'Tube-to-sea conductance', value: round(ou.UA), unit: 'W/(m·K)', explain: 'Heat lost per metre of warm tube per kelvin of fluid-to-sea temperature difference, through bore film, wall, fillers, sheath and external film.' },
        { key: 'Tout', label: 'Fluid temperature at the far end', value: round(Tout), unit: '°C', explain: 'Service-fluid temperature after cooling along the tube. In a long umbilical it reaches the bundle temperature within the first few hundred metres.' },
        { key: 'coolLength', label: 'Thermal decay length of the flow', value: round(Lth), unit: 'm', explain: 'Distance over which the fluid loses 63 % of its excess temperature (ṁ·cp/UA). Thermal effects of the warm fluid are confined to a few of these lengths from the inlet.' },
        { key: 'strain', label: 'Free thermal strain (mean over length)', value: round(epsMean * 100, 3), unit: '%', explain: 'Stiffness-weighted αΔT relative to the installation temperature, averaged along the line. Negative means the umbilical wants to shorten in cold water.' },
        { key: 'strainInlet', label: 'Free thermal strain at the warm inlet', value: round(epsIn * 100, 3), unit: '%', explain: 'Shared axial strain of tube and bundle where the fluid is warmest.' },
        { key: 'restrainedForce', label: 'Fully restrained thermal force', value: round(P / 1e3 * Math.sign(epsMean || 1) * -1), unit: 'kN', explain: 'Axial force if the strain is fully prevented, EA·ε: negative = compression (warming), positive = tension (cooling below the installation temperature).' },
        { key: 'endExpansion', label: 'End movement (per end)', value: round(endMove * 1e3), unit: 'mm', explain: 'Axial movement of each free end after seabed friction has anchored the middle; positive = outwards. The termination or its slack loop must accommodate it.' },
        { key: 'anchorLength', label: 'Friction anchor length', value: round(Math.min(ex.anchorLength, L / 2)), unit: 'm', explain: 'Length from each end over which the umbilical slides; beyond it friction holds the line fully restrained.' },
        { key: 'tubeStress', label: 'Thermal stress in the warm tube', value: round(sigUse / 1e6), unit: 'MPa', explain: 'Axial stress in the warm tube at the inlet caused by differential expansion against the cooler bundle (or by end restraint). Negative = compression.' },
        { key: 'tau', label: 'Thermal time constant', value: round(tau / 3600, 3), unit: 'h', explain: 'Lumped heat capacity times thermal resistance to the sea: the bundle reaches 63 % of a temperature step in this time and about 95 % in three time constants.' },
        { key: 'coolTime', label: 'Cool-down time to the fluid limit', value: round(tCool, 3), unit: 'h', explain: 'Time after shut-in for the stagnant fluid at the inlet to fall to the minimum acceptable temperature; 0 when the limit is outside the range between inlet and ambient.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Radial temperature profile (core losses)', x: { label: 'Radius', unit: 'mm' }, y: { label: 'Temperature', unit: '°C' }, series: [{ name: 'Steady state', x: rS.map((r) => r * 1e3), y: Trad }],
          explain: 'Temperature from the hottest conductor through insulation, bundle and sheath to the sea. The steepest segment is the resistance that controls the conductor temperature.' },
        { type: 'line', title: 'Fluid temperature along the umbilical', x: { label: 'Distance from the topside end', unit: 'm' }, y: { label: 'Temperature', unit: '°C' }, series: [{ name: 'Service fluid', x: xs, y: xs.map(Tf) }, { name: 'Bundle (core losses)', x: [0, L], y: [bundleC, bundleC], dash: true }, { name: 'Installation temperature', x: [0, L], y: [inp.Tinstall, inp.Tinstall], dash: true }],
          explain: 'Exponential cooling of the warm fluid towards the bundle temperature. Where the curve has flattened the tube is in thermal equilibrium with the bundle and carries no differential thermal stress.' },
        { type: 'line', title: 'Free thermal strain along the umbilical', x: { label: 'Distance from the topside end', unit: 'm' }, y: { label: 'Strain', unit: '%' }, series: [{ name: 'Stiffness-weighted αΔT', x: xs, y: strainX }],
          explain: 'Strain the umbilical would take if unrestrained, relative to the installation state. Its integral is the free elongation; seabed friction converts part of it into axial force.' },
        { type: 'line', title: 'Warm-up of the bundle after energising', x: { label: 'Time', unit: 'h' }, y: { label: 'Mean bundle temperature', unit: '°C' }, series: [{ name: 'Lumped response', x: ts.map((t) => t / 3600), y: ts.map((t) => bundleC - (bundleC - inp.Tamb) * Math.exp(-t / tau)) }, { name: 'Steady state', x: [0, 5 * tau / 3600], y: [bundleC, bundleC], dash: true }],
          explain: 'First-order rise of the bundle temperature from ambient when the cores are energised. Short overloads shorter than the time constant heat the bundle far less than the steady value.' },
      ],
      tables: [{ title: 'Temperature rises (series resistances)', explain: 'Each row is one thermal resistance on the path from the hottest conductor to the sea, with the temperature at its inner side. The largest drop identifies where a design change is most effective.',
        columns: [{ key: 'node', label: 'Location' }, { key: 'T', label: 'Temperature', unit: '°C' }, { key: 'dT', label: 'Rise above the next node', unit: 'K' }],
        rows: [
          { node: 'Hottest conductor', T: round(TcoreMax), dT: round(TcoreMax - (ct.sheathInnerTempC + ct.totalLoss * inp.Rbundle)) },
          { node: 'Bundle at the cores', T: round(ct.sheathInnerTempC + ct.totalLoss * inp.Rbundle), dT: round(ct.totalLoss * inp.Rbundle) },
          { node: 'Sheath inner surface', T: round(ct.sheathInnerTempC), dT: round(ct.sheathInnerTempC - ct.surfaceTempC) },
          { node: 'Outer surface', T: round(ct.surfaceTempC), dT: round(ct.surfaceTempC - inp.Tamb) },
          { node: 'Sea water', T: inp.Tamb, dT: 0 },
        ] }],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        `Mean bundle temperature ${fmt(bundleC)} °C against ${inp.Tinstall} °C at installation: ${epsMean < 0 ? 'the umbilical contracts' : 'the umbilical expands'} by ${fmt(Math.abs(epsMean) * 1e6)} microstrain; ${fixed ? `fixed ends carry ${fmt(P / 1e3)} kN` : `each free end moves ${fmt(Math.abs(endMove) * 1e3)} mm with ${fmt(Math.min(ex.anchorLength, L / 2))} m sliding`}.`,
        mDot > 0 ? `The warm fluid cools from ${inp.Tin} °C to within 1 K of the bundle temperature in about ${fmt(Lth * Math.log(Math.max(Math.abs(inp.Tin - bundleC), 1.0001)))} m.` : 'No flow: the tube is at bundle temperature.',
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}. Conductor ${fmt(TcoreMax)} °C, sheath ${fmt(ct.sheathInnerTempC)} °C, tube thermal stress ${fmt(sigUse / 1e6)} MPa, restrained force ${fmt(P / 1e3)} kN.`,
          gov.uc > 1 ? (gov.name.startsWith('Conductor') ? `Cut the losses by ${fmt((1 - 1 / gov.uc) * 100)} % (larger conductor) or reduce the bundle resistance (flooded interstices, cores nearer the sheath).` : gov.name.includes('termination') ? 'Release the ends (sliding mudmat or slack loop sized for the end movement) or raise the termination rating.' : 'Lower the fluid supply temperature or pre-cool it topside; check the combined stress in the tension-and-bending module with the thermal stress added.')
            : 'Add the thermal tube stress to the pressure and tension stresses for the start-up case, and size the termination slack loop for the end movement with margin for friction uncertainty.',
          `${STD_UMB} §7; IEC 60287-1-1 / -2-1; Palmer & King (2008) ch. 12`),
        adv('economic', 'info', `Core losses of ${fmt(ct.totalLoss)} W/m dissipate ${fmt(lossTotal)} kW over ${fmt(L / 1e3)} km — about ${fmt(energyYr)} MWh a year if loaded continuously.`,
          `Conductor losses fall in proportion to the copper area: one size up on the power cores typically cuts them by 20–35 % and lowers conductor temperature (${fmt(TcoreMax)} °C against ${inp.TcoreMax} °C), which extends insulation life. Compare the added copper with the energy saved over the field life.`,
          'Whole-life costing ISO 15686-5 (capitalised losses, as in IEC 60287-3-2 economic optimisation of cable size)',
          `${fmt(energyYr)} MWh/yr of losses is supplied by platform power generation; at a gas-turbine screening intensity of about 0.5 t CO₂/MWh that is ≈ ${fmt(energyYr * 0.5)} t CO₂ a year, far more than the embodied carbon of the extra copper (screening factor, state the site value).`),
      ],
      publish: { 'umbilical.thermalForce': P, 'umbilical.endExpansion': Math.abs(endMove), 'umbilical.coreTemp': TcoreMax },
      quantities: {},
    };
  },
  verify: [
    // No flow, no losses: uniform ΔT = Tamb − Tinstall = 44 − 4 = +40 K. Free strain αΔT = 13e-6 × 40 = 0.052 %;
    // restrained force EA·αΔT = 720e6 × 5.2e-4 = 374.4 kN (compression → −374.4); friction f = 0.4 × 235 = 94 N/m;
    // anchor length P/f = 3983 m > L/2 → short line: δ = (P·L/2 − f·L²/8)/EA = (374 400 × 2500 − 94 × 25e6/8)/720e6 = 0.8920 m.
    { name: 'Uniform temperature change: ε = αΔT, N = EA·αΔT, end expansion of a short line against friction',
      source: 'Thermal expansion relationship; Palmer & King, Subsea Pipeline Engineering ch. 12 (short pipeline, anchor at mid-length)',
      inputs: { flow: 0, Tamb: 44, Tinstall: 4, cores: [{ name: 'none', loss: 0, dCond: 3, dIns: 7, kIns: 0.3 }], Tthreshold: 50 },
      expect: [{ kpi: 'strain', value: 0.052, tol: 0.002 }, { kpi: 'restrainedForce', value: -374.4, tol: 0.002 }, { kpi: 'endExpansion', value: 892.0, tol: 0.003 }] },
    // Series thermal resistances with a prescribed film (very fast current → surface ≈ ambient is avoided by checking the
    // insulation + bundle + sheath rises): one core W = 10 W/m, d = 8/16 mm, k = 0.3 → ΔT_ins = 10·ln2/(2π·0.3) = 3.677 K;
    // bundle 10 × 0.35 = 3.5 K; sheath 10·ln(72.5/66.5)/(2π·0.41) = 0.3353 K → conductor − surface = 7.512 K.
    { name: 'Composite cylinder: conductor-to-surface rise = W·(R_ins + R_bundle + R_sheath)',
      source: 'Incropera et al. §3.3 (radial conduction, R = ln(r2/r1)/(2πk)); IEC 60287-2-1 series form',
      inputs: { cores: [{ name: 'one', loss: 10, dCond: 8, dIns: 16, kIns: 0.3 }], flow: 0 },
      expect: [{ kpi: 'riseInternal', value: 7.512, tol: 0.003 }] },
  ],
};
