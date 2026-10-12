// Electrical loading and power loss of an AC cable link: resistance, inductance, capacitance and sequence impedances,
// exact distributed-parameter load flow with shunt compensation, voltage and current profile, losses and efficiency,
// charging current and critical length, annual energy loss from a load-duration curve, insulation electric stress
// and short-circuit thermal withstand.
import { temperaturesAtCurrent } from '../../engine/electrical/ampacity.js';
import { inductance, sequenceImpedance, chargingCurrent, criticalLength, surgeImpedance, compensatedLink, abcdDistributed, solveTwoPort, lineProfile, lossLoadFactor, threePhasePower } from '../../engine/electrical/circuit.js';
import { gridImpedance, threePhaseFault, faultWithstand } from '../../engine/electrical/shortcircuit.js';
import { stressCheck, acStress, impulseLevels } from '../../engine/electrical/insulation.js';
import { cabs } from '../../engine/electrical/complex.js';
import { cableInputs, buildElec, massesPerMetre, cableQuantities, embodiedCo2, polylineLength, adv, check, worst, sev, fmt, round, run, need, linspace, LIMIT_GENERIC, LIMIT_FE, TB880, CONFIRM, HOURS_YEAR, CO2_COPPER, lossCo2Text } from './_common.js';

const METAL_KEY = { cu: 'Cu', al: 'Al', 'lead-sheath': 'lead', 'steel-armour': 'steel' };
const DEFAULT_LDC = [{ p: 100, h: 2200 }, { p: 75, h: 1500 }, { p: 50, h: 1700 }, { p: 25, h: 2000 }, { p: 5, h: 1366 }];

export default {
  id: 's13-electrical-losses', suite: 13, title: 'Electrical loading and power loss',
  summary: 'Cable electrical constants, voltage drop, current profile, transmission losses and efficiency of the link, annual energy loss, insulation stress and short-circuit withstand.',
  covers: [9], fidelity: 'analytical',
  standards: ['IEC 60287-1-1:2023 §2 (AC resistance, dielectric, sheath and armour losses)', 'IEC 60909-0:2016 (short-circuit currents)', 'IEC 60949:1988+A1 (adiabatic short-circuit withstand)', 'IEC 60840 / IEC 63026 (rated voltages, test levels)', 'IEC 60071-1:2019 (insulation co-ordination)', 'Grainger & Stevenson, Power System Analysis ch. 6 (distributed-parameter line)', 'CIGRE TB 880 (2022) (armour losses)'],
  models: ["Ohm's law and Kirchhoff's laws on the per-phase equivalent", 'Distributed-parameter (hyperbolic ABCD) AC transmission-line model', 'Conductor resistance with temperature, skin and proximity effects', 'Inductance (GMD/GMR) and coaxial capacitance (Gauss)', 'Sequence impedances with Carson earth return', 'Sheath and armour loss factors', 'Dielectric loss', 'Shunt reactive compensation', 'Loss-load factor from a load-duration curve', 'Coaxial insulation electric stress', 'IEC 60909 three-phase fault and IEC 60949 adiabatic withstand'],
  limits: `${LIMIT_GENERIC} Balanced positive-sequence steady state at power frequency: harmonics, resonance, switching transients and converter control are not modelled. ${TB880} The armour-loss reduction factor is an input. ${LIMIT_FE} IEC 60909 / 60949 constants and the stress limits are overridable defaults to confirm against the governing edition.`,
  publishes: ['cable.length', 'cable.lossesFrac', 'cable.voltageDropFrac', 'cable.chargingPower', 'cable.loadCurrent'],
  searchTerms: 'submarine cable transmission losses charging current reactive compensation voltage drop',
  inputs: [
    ...cableInputs(),
    { key: 'route', label: 'Imported route (optional)', type: 'geometry', kinds: ['polyline'], default: null, group: 'Link', role: 'geometry', help: 'Route polyline; its length replaces the route length. Leave empty to use the form value.' },
    { key: 'L', label: 'Route length', unit: 'km', type: 'number', default: 25, min: 0.05, max: 400, group: 'Link', role: 'geometry', design: true },
    { key: 'U', label: 'Sending-end voltage (line-to-line)', unit: 'kV', type: 'number', default: 66, min: 1, max: 550, group: 'Link', role: 'bc', help: 'Electrical boundary condition at the offshore (sending) end.' },
    { key: 'Um', label: 'Highest voltage for equipment Um', unit: 'kV', type: 'number', default: 72.5, min: 1.2, max: 550, group: 'Link', role: 'criteria' },
    { key: 'f', label: 'Frequency', unit: 'Hz', type: 'number', default: 50, min: 16, max: 60, group: 'Link', role: 'bc' },
    { key: 'P', label: 'Transmitted active power (sending end)', unit: 'MW', type: 'number', default: 60, min: 0.1, max: 2500, group: 'Link', role: 'bc', bind: 'energy.ratedPower', design: true },
    { key: 'pf', label: 'Power factor at the sending end', unit: '-', type: 'number', default: 0.98, min: 0.7, max: 1, group: 'Link', role: 'bc' },
    { key: 'leading', label: 'Leading power factor', type: 'bool', default: false, group: 'Link', role: 'bc' },
    { key: 'compS', label: 'Shunt reactor at the sending end', unit: '% of charging power', type: 'number', default: 0, min: 0, max: 100, group: 'Link', role: 'bc' },
    { key: 'compR', label: 'Shunt reactor at the receiving end', unit: '% of charging power', type: 'number', default: 50, min: 0, max: 100, group: 'Link', role: 'bc', help: 'Sending + receiving must not exceed 100 % (over-compensation is rejected).' },
    { key: 'ambient', label: 'Ambient temperature', unit: '°C', type: 'number', default: 12, min: -2, max: 35, group: 'Thermal state', role: 'bc', bind: 'site.Tbed' },
    { key: 'T4', label: 'External thermal resistance', unit: 'K·m/W', type: 'number', default: 0.42, min: 0.005, max: 3, group: 'Thermal state', role: 'bc', help: 'From the thermal-rating module (T4); it sets the conductor temperature and therefore the resistance at the load current.' },
    { key: 'armourLossFactor', label: 'Armour-loss reduction factor', unit: '-', type: 'number', default: 1.0, min: 0.05, max: 1, group: 'Thermal state', role: 'criteria', help: `Multiplier on the IEC 60287 λ2. ${TB880} 1.0 = IEC value.` },
    { key: 'rating', label: 'Continuous current rating', unit: 'A', type: 'number', default: 630, min: 10, max: 5000, group: 'Thermal state', role: 'criteria', bind: 'cable.ampacity' },
    { key: 'rhoEarth', label: 'Resistivity of the return medium', unit: 'Ω·m', type: 'number', default: 1, min: 0.1, max: 1000, group: 'Thermal state', role: 'environment', help: 'Sea water ≈ 0.25, saturated seabed ≈ 1, soil 10–100.' },
    { key: 'ldc', label: 'Load-duration curve', type: 'table', group: 'Load profile', role: 'load', columns: [{ key: 'p', label: 'Power', unit: '% of rated', type: 'number' }, { key: 'h', label: 'Hours per year', type: 'number' }],
      default: DEFAULT_LDC.map((r) => ({ ...r })), help: 'Generic wind-farm duration curve; replace with the project production profile.' },
    { key: 'energyValue', label: 'Indicative value of lost energy', unit: 'USD/MWh', type: 'number', default: 70, min: 0, max: 1000, group: 'Load profile', role: 'cost', help: 'Indicative unit rate used only to express the loss economically.' },
    { key: 'Sk', label: 'Grid short-circuit power at the receiving end', unit: 'MVA', type: 'number', default: 2500, min: 10, max: 1e5, group: 'Short circuit', role: 'bc' },
    { key: 'tk', label: 'Fault clearing time', unit: 's', type: 'number', default: 1, min: 0.05, max: 5, group: 'Short circuit', role: 'bc' },
    { key: 'Ief', label: 'Earth-fault current in one metallic sheath', unit: 'kA', type: 'number', default: 6, min: 0, max: 100, group: 'Short circuit', role: 'bc', help: 'Share of the single-phase fault current returning in the sheath of the faulted core.' },
    { key: 'cFactor', label: 'Voltage factor c', unit: '-', type: 'number', default: 1.1, min: 0.95, max: 1.1, group: 'Short circuit', role: 'criteria', help: `IEC 60909-0 Table 1: 1.1 for maximum short-circuit currents above 1 kV, 1.0 for minimum currents. ${CONFIRM}` },
    { key: 'rOverX', label: 'R/X ratio of the grid feeder', unit: '-', type: 'number', default: 0.1, min: 0.01, max: 1, group: 'Short circuit', role: 'criteria', help: `IEC 60909-0 §6.2: 0.1 for feeders above 35 kV when not known. ${CONFIRM}` },
    { key: 'Kcond', label: 'Adiabatic constant K of the conductor metal', unit: 'A·s½/mm²', type: 'number', default: 0, min: 0, max: 400, group: 'Short circuit', role: 'criteria', help: `0 = IEC 60949 Table I value of the selected metal (copper 226, aluminium 148). ${CONFIRM}` },
    { key: 'betaCond', label: 'Constant β of the conductor metal', unit: 'K', type: 'number', default: 0, min: 0, max: 400, group: 'Short circuit', role: 'criteria', help: '0 = IEC 60949 Table I value (copper 234.5, aluminium 228).' },
    { key: 'thetaSc', label: 'Permissible conductor temperature after a short circuit', unit: '°C', type: 'number', default: 250, min: 130, max: 350, group: 'Acceptance criteria', role: 'criteria', help: `250 °C for XLPE and EPR. ${CONFIRM}` },
    { key: 'thetaScSheath', label: 'Permissible sheath temperature after a short circuit', unit: '°C', type: 'number', default: 170, min: 120, max: 350, group: 'Acceptance criteria', role: 'criteria', help: `170–200 °C for lead alloy depending on the covering. ${CONFIRM}` },
    { key: 'dUmax', label: 'Allowable voltage difference between ends', unit: '%', type: 'number', default: 5, min: 0.5, max: 15, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'Einner', label: 'Allowable stress at the conductor screen', unit: 'kV/mm', type: 'number', default: 7, min: 1, max: 20, group: 'Acceptance criteria', role: 'criteria', help: `Indicative for wet-design 66 kV array cables; dry-design HV cables use higher values. ${CONFIRM}` },
    { key: 'Eouter', label: 'Allowable stress at the insulation screen', unit: 'kV/mm', type: 'number', default: 5, min: 0.5, max: 12, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
  ],
  compute(inp, ctx) {
    const b = buildElec(inp, ctx), g = b.geom, Lm = polylineLength(inp.route) ?? inp.L * 1e3;
    const U = inp.U * 1e3, P = inp.P * 1e6, w = 2 * Math.PI * inp.f;
    need(inp.compS + inp.compR <= 100, 'Over-compensated: the sending and receiving reactors together exceed 100 % of the cable charging power.');
    const sp = threePhasePower({ U, P, powerFactor: inp.pf, leading: inp.leading });
    const system = { f: inp.f, U, thetaMaxC: 90, armourLossReduction: inp.armourLossFactor };
    const installation = { type: 'custom', T4: inp.T4, ambientC: inp.ambient };
    const lineAt = (I) => {
      let t;
      try { t = temperaturesAtCurrent({ cable: g, installation, system, I }); } catch (e) { throw new Error(`No thermal equilibrium at ${fmt(I)} A: the current is far above the cable rating for this thermal resistance.`); }
      const gmd = b.cores === 3 ? g.s : 2 * g.De;      // three-core: core spacing; single-core: indicative spacing of two diameters
      const z = g.Rs20 != null && g.sheathMeanDiameter ? sequenceImpedance({ Rc: t.R, Rs: g.Rs20 * (1 + g.alphaS * (t.sheathC - 20)), dc: g.dc, ds: g.sheathMeanDiameter, gmd, f: inp.f, rhoEarth: inp.rhoEarth }) : null;
      const l = z ? z.Z1.im / w : inductance({ dc: g.dc, gmd });
      // series resistance per phase carrying conductor, sheath and armour losses; shunt conductance carrying the dielectric loss
      const r = t.R * (1 + t.lambda1 + t.lambda2), c = g.C, gSh = w * g.C * g.tanDelta;
      return { t, r, l, c, g: gSh, z, gmd };
    };
    // the cable current depends on the losses and charging current: two passes
    let ln = lineAt(sp.I), link;
    const comp = { sending: inp.compS / 100, receiving: inp.compR / 100 };
    for (let it = 0; it < 2; it++) {
      link = run(() => compensatedLink({ r: ln.r, l: ln.l, c: ln.c, g: ln.g, f: inp.f, length: Lm, U, P, Q: sp.Q, compensation: comp, nProfile: 21 }), 'Load flow');
      ln = lineAt(link.ImaxCable);
    }
    const abcd = abcdDistributed({ r: ln.r, l: ln.l, c: ln.c, g: ln.g, f: inp.f, length: Lm });
    // cable alone (no reactors) from the sending end: profile of voltage and current
    const tp = run(() => solveTwoPort({ abcd, known: 'sending', U, P, Q: sp.Q - comp.sending * link.Qcharging }), 'Two-port solution');
    const prof = lineProfile({ r: ln.r, l: ln.l, c: ln.c, g: ln.g, f: inp.f, length: Lm, Vr: tp.Vr, Ir: tp.Ir, n: 41 });
    const ch = chargingCurrent({ C: g.C, U, f: inp.f, length: Lm });
    const Lcrit = criticalLength({ C: g.C, U, f: inp.f, I: inp.rating });
    const sz = surgeImpedance({ L: ln.l, C: g.C, U, f: inp.f });
    const lossFrac = link.losses / P, dU = (U - link.Ur) / U;
    // annual energy loss: series losses scale with I² (load-duration curve), dielectric loss is constant while energised
    const rows = (inp.ldc || []).filter((r) => r.h > 0);
    need(rows.length > 0 && rows.some((r) => r.p > 0), 'The load-duration curve needs at least one row with positive power and hours.');
    const hours = rows.reduce((a, r) => a + r.h, 0);
    need(hours <= 8784 + 1e-6, `The load-duration curve covers ${fmt(hours)} h — more than a year.`);
    const Wd3 = 3 * ln.t.Wd * Lm, seriesLossPeak = Math.max(link.losses - Wd3, 0);
    const llf = lossLoadFactor({ load: rows.map((r) => r.p), weights: rows.map((r) => r.h), lossAtPeak: seriesLossPeak, hours });
    const eLoss = (llf.energyLoss + Wd3 * HOURS_YEAR) / 1e6, eSent = P * rows.reduce((a, r) => a + r.p / 100 * r.h, 0) / 1e6;   // MWh/yr
    // insulation stress at Um/√3
    const ri = g.dScreen / 2, ro = g.Di / 2, U0m = inp.Um * 1e3 / Math.sqrt(3);
    const es = stressCheck({ U0: U0m, ri, ro, limitInner: inp.Einner * 1e6, limitOuter: inp.Eouter * 1e6 });
    const bil = impulseLevels({ Um: inp.Um * 1e3 });
    // short circuit at the offshore end fed from the grid through the cable (IEC 60909 far-from-generator)
    const zq = gridImpedance({ Un: U, Sk: inp.Sk * 1e6, c: inp.cFactor, rOverX: inp.rOverX });
    const zk = { re: zq.re + ln.t.R * Lm, im: zq.im + w * ln.l * Lm };
    const scFar = threePhaseFault({ Un: U, Zk: zk, c: inp.cFactor, f: inp.f, Tk: inp.tk }), scNear = threePhaseFault({ Un: U, Zk: zq, c: inp.cFactor, f: inp.f, Tk: inp.tk });
    const fwC = faultWithstand({ Ifault: scNear.Ith, material: METAL_KEY[b.m.cond.id], area: inp.area * 1e-6, thetaInitialC: 90, thetaFinalC: inp.thetaSc, duration: inp.tk, ...(inp.Kcond > 0 ? { K: inp.Kcond } : {}), ...(inp.betaCond > 0 ? { beta: inp.betaCond } : {}) });
    const sheathLayer = g.layers.find((q) => q.name === 'sheath');
    const fwS = sheathLayer && inp.Ief > 0 ? faultWithstand({ Ifault: inp.Ief * 1e3, material: METAL_KEY[b.m.sh.id], area: sheathLayer.area, thetaInitialC: Math.min(ln.t.sheathC + 20, 80), thetaFinalC: inp.thetaScSheath, duration: inp.tk }) : null;
    const checks = [
      check('Cable current against the continuous rating', link.ImaxCable, inp.rating, 'A', 'IEC 60287-1-1 §1.4 (rating from the thermal module)', 'Largest current along the cable — load current plus charging current at the more heavily loaded end — against the thermal rating.'),
      check('Voltage difference between the ends', Math.abs(dU) * 100, inp.dUmax, '%', 'Grid-code voltage band; project criterion', 'Difference between sending and receiving voltage at rated power. Outside the band the transformers run out of tap range; a cable can also show a voltage rise at light load (Ferranti effect).'),
      check('Route length against the critical length', Lm / 1e3, Lcrit / 1e3, 'km', 'Charging-current limit l = I/(ω·C·U0)', 'Length at which charging current alone uses the whole rating of a cable fed from one end. Close to it the link needs compensation at both ends or a mid-point reactor — or HVDC.'),
      check('Electric stress at the conductor screen', es.Emax / 1e6, inp.Einner, 'kV/mm', 'Coaxial field E = U0/(r·ln(R/r)); IEC 60840 / IEC 63026 design practice', 'Maximum AC stress in the insulation at Um/√3. Higher stress shortens electrical life (inverse power law) and, in wet designs, accelerates water treeing.'),
      check('Electric stress at the insulation screen', es.Emin / 1e6, inp.Eouter, 'kV/mm', 'Coaxial field at the outer radius; joint and termination design', 'Stress at the insulation screen, which governs accessories (joints and terminations).'),
      check('Conductor short-circuit withstand', scNear.Ith / 1e3, fwC.Ipermissible / 1e3, 'kA', 'IEC 60909-0 §8–10 (thermal equivalent current); IEC 60949 §3 (adiabatic)', 'Thermal equivalent fault current at the grid end against the adiabatic capacity of the conductor for the clearing time. Beyond it the conductor exceeds its permissible temperature and damages the insulation.'),
      ...(fwS ? [check('Metallic-sheath earth-fault withstand', inp.Ief, fwS.Ipermissible / 1e3, 'kA', 'IEC 60949 §3 (adiabatic, sheath)', 'Earth-fault current returning in one sheath against its adiabatic capacity. An undersized sheath or screen melts at the fault and the cable cannot be re-energised.')] : []),
    ];
    const gov = worst(checks);
    const Ps = linspace(0.1, 1.2, 12).map((k) => k * P), lossP = Ps.map((p) => { try { return compensatedLink({ r: ln.r, l: ln.l, c: ln.c, g: ln.g, f: inp.f, length: Lm, U, P: p, Q: sp.Q * p / P, compensation: comp, nProfile: 5 }).losses / p * 100; } catch { return null; } });
    const rr = linspace(ri, ro, 30), Er = rr.map((r) => acStress({ U0: U0m, ri, ro, r }).E / 1e6);
    const mass = massesPerMetre(b), co2 = embodiedCo2(mass, Lm);
    const z1 = ln.z ? cabs(ln.z.Z1) : Math.hypot(ln.t.R, w * ln.l), z0 = ln.z ? cabs(ln.z.Z0) : null;
    const xProf = Array.from(prof.x, (x) => (Lm - x) / 1e3);
    return {
      kpis: [
        { key: 'Rdc20', label: 'Conductor DC resistance at 20 °C', value: round(g.R20 * 1e6, 5), unit: 'mΩ/km', explain: 'IEC 60228 maximum resistance of the conductor; check against the routine test certificate.' },
        { key: 'Rac', label: 'Conductor AC resistance at operating temperature', value: round(ln.t.R * 1e6, 5), unit: 'mΩ/km', explain: 'DC resistance at the conductor temperature increased by skin and proximity effects.' },
        { key: 'Reff', label: 'Effective series resistance incl. sheath and armour losses', value: round(ln.r * 1e6, 5), unit: 'mΩ/km', explain: 'Resistance that reproduces all current-dependent losses of one phase: R·(1 + λ1 + λ2).' },
        { key: 'Lind', label: 'Inductance', value: round(ln.l * 1e6, 5), unit: 'mH/km', explain: 'Positive-sequence inductance per phase (sheaths bonded at both ends).' },
        { key: 'X', label: 'Reactance', value: round(w * ln.l * 1e3, 5), unit: 'Ω/km', explain: 'Positive-sequence series reactance ωL per phase.' },
        { key: 'C', label: 'Capacitance', value: round(g.C * 1e9, 5), unit: 'µF/km', explain: 'Capacitance of one core to its screen; it sets the charging current and reactive power of the cable.' },
        { key: 'Z1', label: 'Positive-sequence impedance', value: round(z1 * 1e3, 5), unit: 'Ω/km', explain: 'Magnitude of the series impedance per phase used in load-flow and fault studies.' },
        ...(z0 != null ? [{ key: 'Z0', label: 'Zero-sequence impedance', value: round(z0 * 1e3, 5), unit: 'Ω/km', explain: 'Series impedance for earth-fault currents returning in sheaths and sea / seabed in parallel.' }] : []),
        { key: 'Zc', label: 'Surge impedance', value: round(sz.Zc), unit: 'Ω', explain: 'Characteristic impedance √(L/C): cables are far below overhead lines (≈ 300–400 Ω), which is why they generate reactive power.' },
        { key: 'conductorTemp', label: 'Conductor temperature at the cable current', value: round(ln.t.conductorC), unit: '°C', explain: 'Electro-thermal operating point used for the resistance (losses and temperature solved together).' },
        { key: 'IcKm', label: 'Charging current per km', value: round(ch.IcPerMetre * 1e3, 5), unit: 'A/km', explain: 'Capacitive current drawn per phase per kilometre at rated voltage, ω·C·U0.' },
        { key: 'Qc', label: 'Total charging power', value: round(ch.Qc / 1e6, 5), unit: 'MVAr', explain: 'Reactive power generated by the cable at rated voltage; the shunt reactors are sized as a fraction of it.' },
        { key: 'Lcrit', label: 'Critical length', value: round(Lcrit / 1e3), unit: 'km', explain: 'Length at which the charging current of an uncompensated cable equals its rating.' },
        { key: 'Iload', label: 'Load current (sending end, before compensation)', value: round(sp.I), unit: 'A', explain: 'Current corresponding to the transmitted power and power factor at the sending end.' },
        { key: 'Imax', label: 'Largest current along the cable', value: round(link.ImaxCable), unit: 'A', limit: inp.rating, sense: 'max', explain: 'Load and charging current combined; in a cable the maximum is at one of the ends and depends on where the reactors are.' },
        { key: 'Ur', label: 'Receiving-end voltage', value: round(link.Ur / 1e3), unit: 'kV', explain: 'Voltage at the onshore (receiving) end for the given sending-end voltage and power.' },
        { key: 'dU', label: 'Voltage drop sending → receiving', value: round(dU * 100, 3), unit: '%', explain: 'Positive = drop, negative = rise towards the receiving end.' },
        { key: 'loss', label: 'Transmission losses at rated power', value: round(link.losses / 1e6), unit: 'MW', explain: 'Conductor, sheath, armour and dielectric losses of the three phases over the route.' },
        { key: 'lossFrac', label: 'Losses as a fraction of transmitted power', value: round(lossFrac * 100), unit: '%', explain: 'Loss fraction at rated power; published for the energy-yield calculation.' },
        { key: 'efficiency', label: 'Power-delivery efficiency', value: round(link.efficiency * 100, 5), unit: '%', explain: 'Received over sent active power at rated load.' },
        { key: 'eLoss', label: 'Annual energy loss', value: round(eLoss), unit: 'MWh/yr', explain: 'Series losses weighted by the square of the load over the duration curve plus dielectric losses for the whole year.' },
        { key: 'eLossFrac', label: 'Annual loss as a fraction of energy sent', value: round(eSent > 0 ? eLoss / eSent * 100 : 0), unit: '%', explain: 'Energy-weighted loss fraction; lower than the rated-power value because losses fall with the square of the load.' },
        { key: 'Emax', label: 'Stress at the conductor screen (at Um)', value: round(es.Emax / 1e6), unit: 'kV/mm', limit: inp.Einner, sense: 'max', explain: 'Maximum electric stress in the insulation at the highest system voltage.' },
        { key: 'Ik', label: 'Three-phase fault current at the grid end', value: round(scNear.Ik / 1e3), unit: 'kA', explain: 'Initial symmetrical short-circuit current I″k fed by the grid into a fault at the receiving end of the cable.' },
        { key: 'IkFar', label: 'Fault current at the far (offshore) end', value: round(scFar.Ik / 1e3), unit: 'kA', explain: 'The cable impedance limits the grid contribution to a fault at the offshore end.' },
        { key: 'Isc', label: 'Conductor adiabatic withstand', value: round(fwC.Ipermissible / 1e3), unit: 'kA', explain: `Current the conductor can carry for ${inp.tk} s from 90 °C without exceeding ${inp.thetaSc} °C (IEC 60949).` },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Voltage and current along the cable', x: { label: 'Distance from the sending end', unit: 'km' }, y: { label: 'Per cent of sending-end value', unit: '%' },
          series: [{ name: 'Voltage', x: xProf, y: Array.from(prof.U, (v) => v / U * 100) }, { name: 'Current', x: xProf, y: Array.from(prof.I, (v) => v / Math.max(tp.Is ? cabs(tp.Is) : 1, 1e-9) * 100) }],
          explain: 'Distributed-parameter profile of the cable alone (reactors at the ends). The current changes along the route as charging current adds vectorially to the load current; the larger end governs the thermal rating.' },
        { type: 'line', title: 'Loss fraction against transmitted power', x: { label: 'Power', unit: 'MW' }, y: { label: 'Losses', unit: '% of power' }, series: [{ name: 'Loss fraction', x: Ps.map((p) => p / 1e6), y: lossP }],
          explain: 'At low power the constant charging-current and dielectric losses dominate (rising percentage towards zero load); at high power the I²R term grows linearly in percentage. The minimum marks the most efficient loading.' },
        { type: 'line', title: 'Electric stress across the insulation at Um', x: { label: 'Radius', unit: 'mm' }, y: { label: 'Stress', unit: 'kV/mm' },
          series: [{ name: 'E(r)', x: rr.map((r) => r * 1e3), y: Er }, { name: 'Limit at the conductor screen', x: [ri * 1e3, ro * 1e3], y: [inp.Einner, inp.Einner], dash: true }],
          explain: 'The stress falls as 1/r from the conductor screen to the insulation screen. The value at the left end is the design stress; a thicker insulation or larger conductor lowers it.' },
        { type: 'bar', title: 'Loss breakdown at the cable current', categories: ['Conductor', 'Sheath', 'Armour', 'Dielectric'], series: [{ name: 'W/m (three phases)', y: [ln.t.losses.conductor, ln.t.losses.sheath, ln.t.losses.armour, ln.t.losses.dielectric].map((v) => round(v)) }],
          explain: 'Per-metre losses of the complete cable. The armour bar carries the IEC 60287 over-estimate for three-core cables unless a reduction factor is applied.' },
      ],
      tables: [
        { title: 'Electrical constants per kilometre', explain: 'Per-phase constants of the cable for load-flow, fault and harmonic studies by others.',
          columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [{ q: 'R (AC, operating)', v: round(ln.t.R * 1e3, 5), u: 'Ω/km' }, { q: 'R effective (with sheath and armour losses)', v: round(ln.r * 1e3, 5), u: 'Ω/km' }, { q: 'X', v: round(w * ln.l * 1e3, 5), u: 'Ω/km' }, { q: 'C', v: round(g.C * 1e9, 5), u: 'µF/km' }, { q: 'B = ωC', v: round(w * g.C * 1e9, 5), u: 'µS/km' },
            ...(ln.z ? [{ q: 'Z1', v: `${fmt(ln.z.Z1.re * 1e3, 4)} + j${fmt(ln.z.Z1.im * 1e3, 4)}`, u: 'Ω/km' }, { q: 'Z0', v: `${fmt(ln.z.Z0.re * 1e3, 4)} + j${fmt(ln.z.Z0.im * 1e3, 4)}`, u: 'Ω/km' }] : []),
            { q: 'Skin / proximity factors', v: `${fmt(ln.t.ys)} / ${fmt(ln.t.yp)}`, u: '-' }, { q: 'λ1 / λ2', v: `${fmt(ln.t.lambda1)} / ${fmt(ln.t.lambda2)}`, u: '-' }, { q: 'Lightning impulse level for Um', v: round(bil.lightning / 1e3), u: 'kV peak' }] },
        { title: 'Load flow at rated power', explain: 'Sending- and receiving-end quantities of the compensated link, with reactor ratings.',
          columns: [{ key: 'q', label: 'Quantity' }, { key: 's', label: 'Sending' }, { key: 'r', label: 'Receiving' }, { key: 'u', label: 'Unit' }],
          rows: [{ q: 'Voltage', s: round(U / 1e3), r: round(link.Ur / 1e3), u: 'kV' }, { q: 'Active power', s: round(P / 1e6), r: round(link.Pr / 1e6), u: 'MW' }, { q: 'Reactive power', s: round(sp.Q / 1e6), r: round(link.Qr / 1e6), u: 'MVAr' }, { q: 'Cable current', s: round(link.Is), r: round(link.Ir), u: 'A' }, { q: 'Shunt reactor', s: round(link.reactors.sending / 1e6), r: round(link.reactors.receiving / 1e6), u: 'MVAr' }] },
      ],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        `Charging power ${fmt(ch.Qc / 1e6)} MVAr over ${fmt(Lm / 1e3)} km (${fmt(Lm / Lcrit * 100)} % of the critical length); compensation ${inp.compS} % / ${inp.compR} % at sending / receiving end.`,
        `Loss-load factor ${fmt(llf.lossLoadFactor)} (load factor ${fmt(llf.loadFactor)}): annual loss ${fmt(eLoss)} MWh of ${fmt(eSent)} MWh sent.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name}: ${fmt(gov.uc)}. Cable current ${fmt(link.ImaxCable)} A of ${inp.rating} A, losses ${fmt(lossFrac * 100)} % at rated power, voltage change ${fmt(dU * 100)} %, stress ${fmt(es.Emax / 1e6)} kV/mm.`,
          gov.uc > 1 ? (gov.name.startsWith('Cable current') ? `Move compensation towards the more loaded end or increase the conductor to about ${fmt(inp.area * gov.uc ** 2)} mm² (losses and heating scale with I²/A).` : gov.name.startsWith('Route length') ? 'Add mid-point compensation, raise the voltage level or assess HVDC for this distance.' : gov.name.includes('stress') ? `Increase the insulation thickness to at least ${fmt(es.thicknessRequired * 1e3)} mm on this conductor.` : gov.name.includes('short-circuit') || gov.name.includes('fault') ? 'Shorten the clearing time or increase the cross-section of the conductor / sheath.' : 'Adjust transformer taps or the reactive compensation split to bring the voltage back into band.')
            : 'Issue the per-kilometre constants to the grid-study team, and confirm resonance and energisation transients (not covered here) for long cables with shunt reactors.',
          'IEC 60287-1-1:2023 §2; IEC 60909-0:2016; IEC 60949; Grainger & Stevenson (1994) ch. 6'),
        adv('economic', 'info', `Annual energy loss ≈ ${fmt(eLoss)} MWh (${fmt(eSent > 0 ? eLoss / eSent * 100 : 0)} % of the energy sent) — indicatively ${fmt(eLoss * inp.energyValue / 1e3)} k USD a year at the entered ${inp.energyValue} USD/MWh.`,
          `One conductor size up cuts the I²R part by roughly ${fmt((1 - inp.area / (inp.area * 1.27)) * 100)} % (≈ ${fmt(llf.energyLoss / 1e6 * (1 - 1 / 1.27))} MWh/yr) for about ${fmt(mass.copper * 0.27 * Lm / 1e3)} t more conductor metal; compare the capitalised losses over the project life with that metal (Kelvin's law, IEC 60287-3-2).`,
          'IEC 60287-3-2 (economic optimisation of power cable size); LCOE conventions (losses reduce net energy yield)',
          `${lossCo2Text(eLoss)} Against that, ${fmt(mass.copper * 0.27 * Lm / 1e3)} t of extra copper embodies ≈ ${fmt(mass.copper * 0.27 * Lm / 1e3 * CO2_COPPER)} t CO₂e once (screening ${CO2_COPPER} t CO₂e/t); the cable as designed embodies ≈ ${fmt(co2)} t CO₂e.`),
      ],
      publish: { 'cable.length': Lm, 'cable.lossesFrac': lossFrac, 'cable.voltageDropFrac': dU, 'cable.chargingPower': ch.Qc, 'cable.loadCurrent': link.ImaxCable },
      quantities: cableQuantities(mass, Lm),
    };
  },
  verify: [
    // Default cable core: dc = 29.854 mm, screen 1 mm → ri = 15.927 mm; insulation 9 mm → ro = 24.927 mm; ln(ro/ri) = 0.44795.
    // C = 2.5e-9/(18 × 0.44795) = 3.1006e-10 F/m = 0.31006 µF/km (IEC 60287-1-1 §2.2).
    // Charging current ω·C·U0 = 314.159 × 3.1006e-10 × 38105 = 3.712e-3 A/m = 3.712 A/km; Qc = ω·C·U²·l = 314.159 × 3.1006e-10 × 66000² × 25000 = 10.61 MVAr.
    // Critical length = I/(ω·C·U0) = 630/3.712 = 169.7 km. Stress at Um: (72.5/√3)/(15.927 × 0.44795) = 5.867 kV/mm.
    // Adiabatic: 226 × 630 × sqrt(ln((250 + 234.5)/(90 + 234.5))/1 s) = 142380 × 0.6331 = 90.14 kA (IEC 60949).
    { name: 'Coaxial capacitance, charging current and power, critical length, insulation stress and adiabatic short-circuit current',
      source: 'IEC 60287-1-1 §2.2; Gauss law for a coaxial field; IEC 60949 §3 (hand calculation)',
      inputs: {},
      expect: [{ kpi: 'C', value: 0.31006, tol: 0.002 }, { kpi: 'IcKm', value: 3.712, tol: 0.002 }, { kpi: 'Qc', value: 10.61, tol: 0.002 }, { kpi: 'Lcrit', value: 169.7, tol: 0.003 }, { kpi: 'Emax', value: 5.867, tol: 0.002 }, { kpi: 'Isc', value: 90.14, tol: 0.003 }] },
    // Short link limit: at 100 m the cable is electrically short, so losses → 3·I²·R_eff·l and the efficiency → 1 − loss/P;
    // the load current is P/(√3·U·cosφ) = 60e6/(1.7321 × 66e3 × 0.98) = 535.6 A.
    { name: 'Balanced three-phase load current I = P/(√3·U·cosφ)',
      source: 'Three-phase power relations (Grainger & Stevenson §1)',
      inputs: { L: 0.1, compR: 0 },
      expect: [{ kpi: 'Iload', value: 535.6, tol: 0.002 }, { kpi: 'Imax', value: 535.6, tol: 0.01 }] },
  ],
};
