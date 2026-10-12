// Cable failure and reliability: failure frequency of the circuit from internal, external, joint and termination
// rates, repair time and availability, fault tree of loss of transmission with importance measures, energy not
// delivered, Weibull analysis of service records, electro-thermal ageing of the insulation with the size effect on
// breakdown probability, an FMEA ranking and the installation / operating / service-life compliance summary.
import { evaluateFaultTree, fmeca } from '../../engine/ram/faulttree.js';
import { fitWeibullMle, weibullFunctions, bLife, sparesForServiceLevel } from '../../engine/ram/lifedata.js';
import { combinedLife, thermalLife, weibullBreakdownProbability, weibullEnlargement } from '../../engine/electrical/insulation.js';
import { adv, check, worst, sev, fmt, round, run, need, linspace, CONFIRM, HOURS_YEAR, CO2_COPPER, vesselCo2Text } from './_common.js';

const DEFAULT_RECORDS = [
  { age: 0.4, status: 'F' }, { age: 1.1, status: 'F' }, { age: 2.6, status: 'F' }, { age: 5.5, status: 'F' }, { age: 9, status: 'F' }, { age: 14, status: 'F' },
  { age: 4, status: 'S' }, { age: 6, status: 'S' }, { age: 8, status: 'S' }, { age: 10, status: 'S' }, { age: 12, status: 'S' }, { age: 12, status: 'S' },
  { age: 15, status: 'S' }, { age: 16, status: 'S' }, { age: 18, status: 'S' }, { age: 20, status: 'S' },
];
const DEFAULT_FMEA = [
  { item: 'Cable (route)', mode: 'Anchor strike or drag', cause: 'external', severity: 9, detection: 3 },
  { item: 'Cable (route)', mode: 'Fishing-gear hooking or abrasion', cause: 'external', severity: 7, detection: 5 },
  { item: 'Cable (route)', mode: 'Insulation breakdown at a defect', cause: 'internal', severity: 9, detection: 7 },
  { item: 'Cable (route)', mode: 'Thermal overload at a hotspot', cause: 'internal', severity: 6, detection: 4 },
  { item: 'Joint', mode: 'Interface breakdown / water ingress', cause: 'joint', severity: 9, detection: 7 },
  { item: 'Termination', mode: 'Tracking, moisture, hang-off damage', cause: 'termination', severity: 6, detection: 3 },
  { item: 'Dynamic section', mode: 'Armour or sheath fatigue', cause: 'internal', severity: 8, detection: 6 },
];
const CAUSES = ['internal', 'external', 'joint', 'termination'];

export default {
  id: 's13-failure-reliability', suite: 13, title: 'Cable failure and reliability',
  summary: 'Failure frequency, repair time, availability and energy not delivered of the cable circuit, with the fault tree of loss of transmission, the Weibull trend of service records, insulation ageing and the compliance summary.',
  covers: [14], fidelity: 'analytical',
  standards: ['CIGRE TB 379 (2009) Update of service experience of HV underground and submarine cable systems', 'CIGRE TB 815 (2020) Update of service experience of HV underground and submarine cable systems', 'IEC 61025:2006 (fault tree analysis)', 'IEC 60812:2018 (FMEA / FMECA)', 'IEC 61649:2008 (Weibull analysis)', 'IEC 60216-1 / IEC 60505 (thermal and electrical endurance of insulation)', 'IEC 62539 (statistical analysis of breakdown data)', 'DNV-RP-0360 (2016) §8 (operation, repair preparedness)'],
  models: ['Homogeneous Poisson failure process with constant rates', 'Two-state repairable components (failure rate, repair rate)', 'Fault-tree analysis by binary decision diagram with minimal cut sets and importance measures', 'β-factor common-cause model for parallel circuits', 'Weibull reliability distribution fitted by maximum likelihood to censored life data', 'Arrhenius thermal ageing and inverse-power-law electrical ageing (insulation ageing models)', 'Weibull weakest-link enlargement law for breakdown probability', 'Failure mode and effects analysis with risk priority number', 'Probabilistic external damage rate from the burial risk assessment', 'Poisson spares provisioning'],
  limits: 'Failure rates are fleet statistics with wide scatter and are overridable defaults: they do not know this cable, its installation quality or its route — use project data where it exists. Rates are constant in time for the availability figures (the Weibull fit shows whether that is fair), repairs are independent with one mean repair time, and weather windows enter only through the waiting time. The ageing model multiplies thermal and electrical acceleration without synergy and ignores water treeing; condition data (partial discharge, DTS, sheath tests) must update it in service.',
  publishes: ['cable.availability', 'cable.failureRate', 'cable.energyNotDelivered', 'cable.insulationLife'],
  searchTerms: 'submarine cable failure rate statistics repair time availability CIGRE TB 815',
  inputs: [
    { key: 'L', label: 'Route length per circuit', unit: 'km', type: 'number', default: 25, min: 0.05, max: 800, group: 'Circuit', role: 'geometry', bind: 'cable.length' },
    { key: 'layout', label: 'Circuit arrangement', type: 'select', default: 'single', group: 'Circuit', role: 'bc', options: [{ value: 'single', label: 'Single circuit (no redundancy)' }, { value: 'dual', label: 'Two parallel circuits' }], help: 'Boundary of the system whose loss is assessed: the cable connection between the two terminations.' },
    { key: 'share', label: 'Capacity of one circuit when the other is out', unit: '% of total', type: 'number', default: 60, min: 0, max: 100, group: 'Circuit', role: 'bc', showIf: { key: 'layout', in: ['dual'] }, help: 'Share of the full output the surviving circuit can carry.' },
    { key: 'betaCC', label: 'Common-cause fraction of external failures', unit: '-', type: 'number', default: 0.1, min: 0, max: 1, group: 'Circuit', role: 'operation', showIf: { key: 'layout', in: ['dual'] }, help: 'Share of external damage events (anchor drag) that take out both circuits; it falls with circuit spacing.' },
    { key: 'nJoints', label: 'Joints per circuit', type: 'int', default: 4, min: 0, max: 400, group: 'Circuit', role: 'geometry' },
    { key: 'nTerm', label: 'Terminations per circuit', type: 'int', default: 2, min: 0, max: 20, group: 'Circuit', role: 'geometry' },
    { key: 'life', label: 'Design life', unit: 'years', type: 'number', default: 30, min: 1, max: 60, group: 'Circuit', role: 'criteria', bind: 'life.designYears' },
    { key: 'age0', label: 'Age of the cable today', unit: 'years', type: 'number', default: 0, min: 0, max: 60, group: 'Circuit', role: 'ic', help: 'Initial condition for the degradation assessment: years already in service (0 = new).' },
    { key: 'lamInt', label: 'Internal failure rate of the cable', unit: 'per 100 km·yr', type: 'number', default: 0.03, min: 0, max: 5, group: 'Failure rates', role: 'operation', uncertain: { dist: 'lognormal', cov: 0.6 }, help: `CIGRE TB 379 / TB 815 report internal failures of extruded submarine cables of the order of 0.02–0.1 per 100 circuit-km per year. ${CONFIRM}` },
    { key: 'lamExt', label: 'External damage rate over the route', unit: 'per year', type: 'number', default: 0.02, min: 0, max: 5, group: 'Failure rates', role: 'operation', bind: 'cable.strikeRate', uncertain: { dist: 'lognormal', cov: 0.8 }, help: 'Anchor and fishing damage per year for the whole route — from the burial risk assessment. Service statistics give 0.05–0.3 per 100 km·yr for unburied or shallow cables.' },
    { key: 'lamJoint', label: 'Failure rate per joint', unit: 'per year', type: 'number', default: 0.002, min: 0, max: 1, group: 'Failure rates', role: 'operation', uncertain: { dist: 'lognormal', cov: 0.6 }, help: `Indicative; field (repair) joints are several times worse than factory joints. ${CONFIRM}` },
    { key: 'lamTerm', label: 'Failure rate per termination', unit: 'per year', type: 'number', default: 0.002, min: 0, max: 1, group: 'Failure rates', role: 'operation', uncertain: { dist: 'lognormal', cov: 0.6 } },
    { key: 'records', label: 'Service records for the Weibull trend', type: 'table', group: 'Failure rates', role: 'operation', columns: [{ key: 'age', label: 'Age at failure or today', unit: 'years', type: 'number' }, { key: 'status', label: 'F = failed, S = still in service', type: 'text' }],
      default: DEFAULT_RECORDS.map((r) => ({ ...r })), help: 'Historical failures and surviving sections of comparable cables (operator fleet data). The default set is illustrative.' },
    { key: 'mob', label: 'Fault location and vessel mobilisation', unit: 'days', type: 'number', default: 30, min: 0, max: 365, group: 'Repair', role: 'operation', uncertain: { dist: 'lognormal', cov: 0.4 } },
    { key: 'wait', label: 'Waiting on weather', unit: 'days', type: 'number', default: 15, min: 0, max: 200, group: 'Repair', role: 'operation', uncertain: { dist: 'lognormal', cov: 0.6 } },
    { key: 'rep', label: 'Offshore repair, test and re-burial', unit: 'days', type: 'number', default: 20, min: 1, max: 200, group: 'Repair', role: 'operation' },
    { key: 'spare', label: 'Spare cable and joint kits in store', type: 'bool', default: true, group: 'Repair', role: 'operation' },
    { key: 'lead', label: 'Lead time without a spare', unit: 'days', type: 'number', default: 180, min: 0, max: 900, group: 'Repair', role: 'operation', showIf: { key: 'spare', in: [false] } },
    { key: 'spareL', label: 'Spare cable length held', unit: 'm', type: 'number', default: 1000, min: 0, max: 20000, group: 'Repair', role: 'operation', showIf: { key: 'spare', in: [true] }, help: 'A repair needs about 2–3 water depths of new cable plus two joints.' },
    { key: 'P', label: 'Rated transmitted power', unit: 'MW', type: 'number', default: 60, min: 0.1, max: 5000, group: 'Energy', role: 'bc', bind: 'energy.ratedPower' },
    { key: 'cf', label: 'Capacity factor', unit: '-', type: 'number', default: 0.45, min: 0.05, max: 1, group: 'Energy', role: 'bc', bind: 'energy.capacityFactor' },
    { key: 'energyValue', label: 'Indicative value of energy', unit: 'USD/MWh', type: 'number', default: 70, min: 0, max: 1000, group: 'Energy', role: 'cost', help: 'Indicative unit rate used only to express lost production economically.' },
    { key: 'dayRate', label: 'Indicative repair spread day rate', unit: 'USD/day', type: 'number', default: 150000, min: 0, max: 2e6, group: 'Energy', role: 'cost', help: 'Indicative; repair vessel with jointing spread.' },
    { key: 'thetaOp', label: 'Conductor temperature in service (energy-weighted)', unit: '°C', type: 'number', default: 75, min: 0, max: 130, group: 'Insulation ageing', role: 'ic', bind: 'cable.conductorTemp', help: 'Thermal state the insulation lives at; the peak value is conservative.' },
    { key: 'thetaRef', label: 'Reference temperature of the design life', unit: '°C', type: 'number', default: 90, min: 40, max: 110, group: 'Insulation ageing', role: 'criteria' },
    { key: 'lifeRef', label: 'Insulation life at the reference conditions', unit: 'years', type: 'number', default: 30, min: 5, max: 80, group: 'Insulation ageing', role: 'criteria', help: `Qualification basis: life at the rated temperature and design stress (IEC 60840 / IEC 63026 prequalification intent). ${CONFIRM}` },
    { key: 'Ea', label: 'Activation energy of thermal ageing', unit: 'eV', type: 'number', default: 1.1, min: 0.5, max: 2, group: 'Insulation ageing', role: 'material', uncertain: { dist: 'normal', cov: 0.1 }, help: '1.0–1.3 eV for XLPE (material specific).' },
    { key: 'E', label: 'Operating electric stress', unit: 'kV/mm', type: 'number', default: 5.5, min: 0.5, max: 30, group: 'Insulation ageing', role: 'bc', help: 'Stress at the conductor screen at the operating voltage (electrical module).' },
    { key: 'Eref', label: 'Design stress of the qualification', unit: 'kV/mm', type: 'number', default: 5.9, min: 0.5, max: 30, group: 'Insulation ageing', role: 'criteria' },
    { key: 'nLife', label: 'Electrical life exponent', unit: '-', type: 'number', default: 12, min: 5, max: 25, group: 'Insulation ageing', role: 'material', help: '9–15 for extruded AC insulation (IEC 60505; CIGRE TB 496 for DC).' },
    { key: 'E0', label: '63 % breakdown stress of the test sample', unit: 'kV/mm', type: 'number', default: 40, min: 5, max: 200, group: 'Insulation ageing', role: 'material', help: 'Weibull scale of AC breakdown tests on model or full-size samples.' },
    { key: 'bW', label: 'Weibull shape of breakdown stress', unit: '-', type: 'number', default: 10, min: 3, max: 30, group: 'Insulation ageing', role: 'material' },
    { key: 'L0', label: 'Test sample length', unit: 'm', type: 'number', default: 10, min: 0.1, max: 1000, group: 'Insulation ageing', role: 'material' },
    { key: 'I', label: 'Load current', unit: 'A', type: 'number', default: 580, min: 0, max: 6000, group: 'Results of other modules', role: 'bc', bind: 'cable.loadCurrent' },
    { key: 'ampacity', label: 'Continuous rating', unit: 'A', type: 'number', default: 630, min: 1, max: 6000, group: 'Results of other modules', role: 'criteria', bind: 'cable.ampacity' },
    { key: 'Tmax', label: 'Largest tension (installation or service)', unit: 'kN', type: 'number', default: 90, min: 0, max: 5000, group: 'Results of other modules', role: 'bc', bind: 'cable.topTension' },
    { key: 'Tallow', label: 'Allowable tension', unit: 'kN', type: 'number', default: 250, min: 1, max: 5000, group: 'Results of other modules', role: 'criteria', bind: 'cable.allowTension' },
    { key: 'mbrMin', label: 'Tightest bend radius found', unit: 'm', type: 'number', default: 6, min: 0.1, max: 1e4, group: 'Results of other modules', role: 'bc', bind: 'cable.minBendRadius' },
    { key: 'mbrAllow', label: 'Allowable bend radius', unit: 'm', type: 'number', default: 3.2, min: 0.1, max: 60, group: 'Results of other modules', role: 'criteria', bind: 'cable.mbrDynamic' },
    { key: 'fatigueYears', label: 'Calculated fatigue life', unit: 'years', type: 'number', default: 400, min: 0.1, max: 1e6, group: 'Results of other modules', role: 'ic', bind: 'life.fatigueYears' },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 10, min: 1, max: 20, group: 'Acceptance criteria', role: 'criteria', help: `10 is usual for dynamic cables that cannot be inspected (DNV-ST-0359 / DNV-RP-0360). ${CONFIRM}` },
    { key: 'availTarget', label: 'Required availability of the connection', unit: '%', type: 'number', default: 99, min: 80, max: 99.99, group: 'Acceptance criteria', role: 'criteria', help: 'Energy-weighted availability from the business case or the grid connection agreement.' },
    { key: 'serviceLevel', label: 'Confidence level for spares provisioning', unit: '%', type: 'number', default: 95, min: 50, max: 99.9, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'fmea', label: 'Failure modes', type: 'table', group: 'FMEA', role: 'operation', columns: [{ key: 'item', label: 'Item', type: 'text' }, { key: 'mode', label: 'Failure mode', type: 'text' }, { key: 'cause', label: 'Rate class (internal / external / joint / termination)', type: 'text' }, { key: 'severity', label: 'Severity 1–10', type: 'number' }, { key: 'detection', label: 'Detection 1–10', type: 'number' }],
      default: DEFAULT_FMEA.map((r) => ({ ...r })), help: 'Severity and detection are judgement ranks (10 = worst / undetectable); occurrence is derived from the failure rate of the class.' },
  ],
  compute(inp) {
    const T = inp.life, mttrDays = inp.mob + inp.wait + inp.rep + (inp.spare ? 0 : inp.lead);
    need(inp.age0 < T, 'The cable is already older than the design life: enter the extended life as the design life to assess a life extension.');
    const mu = 365.25 / mttrDays;                                   // repairs per year
    const rate = { internal: inp.lamInt * inp.L / 100, external: inp.lamExt, joint: inp.lamJoint * inp.nJoints, termination: inp.lamTerm * inp.nTerm };
    const lam = rate.internal + rate.external + rate.joint + rate.termination;
    need(lam > 0, 'All failure rates are zero: nothing to assess. Enter at least one non-zero rate.');
    const dual = inp.layout === 'dual', beta = dual ? inp.betaCC : 0;
    // basic events of one circuit: repairable two-state units, q = λ/(λ + μ)
    const ev = (id, l) => ({ type: 'basic', id, lambda: l, mu });
    const circuit = (c) => ({ type: 'or', children: [ev(`${c}: internal fault`, rate.internal), ev(`${c}: external damage`, rate.external * (1 - beta)), ev(`${c}: joint failure`, rate.joint), ev(`${c}: termination failure`, rate.termination)] });
    const tree = dual ? { type: 'or', children: [{ type: 'and', children: [circuit('A'), circuit('B')] }, ev('Common-cause external damage', rate.external * beta)] } : circuit('A');
    const ft = run(() => evaluateFaultTree(tree), 'Fault tree');
    const uOne = run(() => evaluateFaultTree(dual ? { type: 'or', children: [circuit('A'), ev('Common-cause external damage', rate.external * beta)] } : circuit('A'), { cutSets: false }), 'Fault tree').probability;
    const uTop = ft.probability;                                      // probability that no power can be transmitted
    // energy not delivered: total loss + (dual) partial loss while exactly one circuit is out
    const Eyear = inp.P * inp.cf * HOURS_YEAR;                        // MWh/yr
    const pOne = dual ? Math.max(0, 2 * (uOne - uTop)) : 0;
    const lostFrac = uTop + pOne * (1 - inp.share / 100);
    const end = lostFrac * Eyear, avail = 1 - lostFrac;
    const lamSys = dual ? 2 * lam : lam;                              // repair events per year, all circuits
    const Nlife = lamSys * (T - inp.age0), pNone = Math.exp(-Nlife);
    const sp = sparesForServiceLevel({ failureRate: lamSys, period: T - inp.age0, serviceLevel: inp.serviceLevel / 100 });   // no replenishment assumed
    // Weibull trend of the service records
    const recs = (inp.records || []).filter((r) => r.age > 0), nF = recs.filter((r) => String(r.status).trim().toUpperCase().startsWith('F')).length;
    let wb = null;
    if (nF >= 3) { try { wb = fitWeibullMle({ times: recs.map((r) => r.age), censored: recs.map((r) => !String(r.status).trim().toUpperCase().startsWith('F')) }); } catch { wb = null; } }
    if (wb && !(Number.isFinite(wb.shape) && Number.isFinite(wb.scale) && wb.shape > 0)) wb = null;
    const trend = !wb ? 'not fitted' : wb.shape < 0.9 ? 'early failures (decreasing hazard)' : wb.shape > 1.2 ? 'wear-out (increasing hazard)' : 'random failures (constant hazard)';
    // insulation ageing
    const lifeIns = combinedLife({ E: inp.E, Eref: inp.Eref, n: inp.nLife, thetaC: inp.thetaOp, thetaRefC: inp.thetaRef, lifeRef: inp.lifeRef, activationEnergyEv: inp.Ea });
    const afT = inp.lifeRef / thermalLife({ thetaC: inp.thetaOp, thetaRefC: inp.thetaRef, lifeRef: inp.lifeRef, activationEnergyEv: inp.Ea });
    const sizeRatio = inp.L * 1e3 * (dual ? 2 : 1) * 3 / inp.L0;      // three cores per circuit
    const pBd = weibullBreakdownProbability({ E: inp.E * 1e6, E0: inp.E0 * 1e6, b: inp.bW, sizeRatio });
    const enl = weibullEnlargement({ E0: inp.E0 * 1e6, b: inp.bW, sizeRatio, probability: 0.01 });
    // FMEA: occurrence rank from the probability of at least one event of the class in the design life
    const occ = (cause) => { const p = 1 - Math.exp(-(rate[cause] ?? 0) * T * (dual ? 2 : 1)); return Math.min(10, Math.max(1, Math.round(1 + 9 * Math.sqrt(p)))); };
    const fmRows = (inp.fmea || []).filter((r) => r.mode).map((r, i) => {
      const cause = String(r.cause || '').trim().toLowerCase();
      need(CAUSES.includes(cause), `FMEA row ${i + 1}: the rate class must be one of ${CAUSES.join(', ')}.`);
      need(r.severity >= 1 && r.severity <= 10 && r.detection >= 1 && r.detection <= 10, `FMEA row ${i + 1}: severity and detection must be between 1 and 10.`);
      return { item: String(r.item || ''), mode: String(r.mode), cause, severity: Number(r.severity), detection: Number(r.detection), occurrence: occ(cause) };
    });
    const fm = fmRows.length ? fmeca(fmRows) : { byRpn: [] };
    // compliance summary
    const comp = [
      { stage: 'Installation', item: 'Tension', demand: inp.Tmax, capacity: inp.Tallow, unit: 'kN', uc: inp.Tmax / inp.Tallow },
      { stage: 'Installation', item: 'Minimum bend radius', demand: inp.mbrAllow, capacity: inp.mbrMin, unit: 'm', uc: inp.mbrAllow / inp.mbrMin },
      { stage: 'Operation', item: 'Thermal rating', demand: inp.I, capacity: inp.ampacity, unit: 'A', uc: inp.I / inp.ampacity },
      { stage: 'Operation', item: 'Availability', demand: (1 - avail) * 100, capacity: 100 - inp.availTarget, unit: '% unavailable', uc: (1 - avail) / Math.max(1 - inp.availTarget / 100, 1e-9) },
      { stage: 'Service life', item: 'Fatigue life', demand: T * inp.dff, capacity: inp.fatigueYears, unit: 'years', uc: T * inp.dff / inp.fatigueYears },
      { stage: 'Service life', item: 'Insulation life', demand: T, capacity: lifeIns, unit: 'years', uc: T / lifeIns },
    ];
    const checks = [
      check('Availability of the connection', (1 - avail) * 100, 100 - inp.availTarget, '% unavailable', 'Frequency–duration method (IEC 61078 / IEC 61025); CIGRE TB 379 / TB 815 statistics', 'Energy-weighted unavailability against the allowance of the business case. Every failure costs the whole repair time, so availability is driven as much by repair preparedness as by failure rate.'),
      check('Insulation electro-thermal life', T, lifeIns, 'years', 'IEC 60216-1 (Arrhenius); IEC 60505 (inverse power law)', 'Design life against the life of the insulation at the service temperature and stress. A shorter life means the cable ages out before the asset it serves.'),
      check('Fatigue life with the design fatigue factor', T * inp.dff, inp.fatigueYears, 'years', 'DNV-ST-0359 (2021) / DNV-RP-0360 (design fatigue factor); fatigue module of this suite', 'Design life × design fatigue factor against the calculated fatigue life of the governing component (armour, sheath or conductor).'),
      check('Thermal rating utilisation', inp.I, inp.ampacity, 'A', 'IEC 60287-1-1 §1.4 (rating from the thermal module)', 'Load current against the continuous rating: sustained operation above it is what turns the ageing model against the cable.'),
      check('Tension utilisation', inp.Tmax, inp.Tallow, 'kN', 'CIGRE TB 623 (tensile capacity); mechanical modules of this suite', 'Largest tension seen in installation or service against the allowable: overstrained cables fail early and are over-represented in the infant-mortality statistics.'),
      check('Bend radius compliance', inp.mbrAllow, inp.mbrMin, 'm', 'DNV-RP-0360 (MBR); CIGRE TB 623', 'Allowable bend radius against the tightest radius found; over-bending damages screens and sheath and shows up later as an internal fault.'),
    ];
    const gov = worst(checks);
    const yrs = linspace(0, T, 31);
    const repairDaysLife = Nlife * inp.rep, vesselDaysYear = lamSys * (inp.rep + 0.2 * inp.mob);     // offshore spread days, 20 % of mobilisation at sea
    const cut = (ft.cutSets || []).slice(0, 8);
    const imp = [...ft.importance].sort((a, c) => (c.fussellVesely ?? c.criticality) - (a.fussellVesely ?? a.criticality));
    const sharePct = (v) => round(v / lam * 100);
    const spareCu = inp.spare ? 17 * inp.spareL / 1e3 : 0;           // ≈ 17 kg/m of copper in a 3 × 630 mm² cable
    return {
      kpis: [
        { key: 'lambda', label: 'Failure frequency per circuit', value: round(lam, 5), unit: '1/yr', explain: 'Expected failures per year of one circuit: internal faults, external damage, joints and terminations added as a series system.' },
        { key: 'mtbf', label: 'Mean time between failures per circuit', value: round(1 / lam), unit: 'years', explain: 'Reciprocal of the failure frequency. Compare it with the design life: a value of the same order means at least one repair should be planned for.' },
        { key: 'mttr', label: 'Mean time to restore', unit: 'days', value: round(mttrDays), explain: 'Outage per failure: fault location and mobilisation, waiting on weather, offshore repair — plus the manufacturing lead time when no spare is held.' },
        { key: 'failuresLife', label: 'Expected failures over the remaining life', value: round(Nlife, 5), unit: '-', explain: 'Failure frequency of all circuits × remaining years: the number of repair campaigns to budget for.' },
        { key: 'pNoFailure', label: 'Probability of no failure in the remaining life', value: round(pNone, 5), unit: '-', explain: 'exp(−expected failures) for a Poisson process. Even well-protected export cables seldom exceed 0.5 over 25–30 years.' },
        { key: 'availability', label: 'Energy-weighted availability', value: round(avail * 100, 6), unit: '%', limit: inp.availTarget, sense: 'min', explain: 'Share of the annual energy the cable connection can deliver. With two circuits the partial capacity of the surviving one is credited.' },
        { key: 'unavailTotal', label: 'Probability of total loss of transmission', value: round(uTop * 100, 5), unit: '%', explain: 'Top event of the fault tree: no circuit available. For two circuits it is dominated by common-cause damage.' },
        { key: 'downtime', label: 'Equivalent full outage', value: round(lostFrac * HOURS_YEAR), unit: 'h/yr', explain: 'Unavailability expressed as hours per year of complete outage.' },
        { key: 'end', label: 'Energy not delivered', value: round(end), unit: 'MWh/yr', explain: 'Expected annual production lost to cable outages at the entered capacity factor.' },
        { key: 'spares', label: 'Repair sets to cover the remaining life', value: sp.spares, unit: '-', explain: `Number of repair sets (spare length and two joints) that covers the failures of the remaining life with ${inp.serviceLevel} % confidence if none is re-ordered; hold the first in store and re-order after each use.` },
        { key: 'weibullShape', label: 'Weibull shape of the service records', value: round(wb ? wb.shape : 1), unit: '-', explain: 'Below 1: failures concentrate early in life (installation damage, defects); about 1: random (external); above 1: wear-out. 1 is shown when too few failures are recorded to fit.' },
        { key: 'weibullScale', label: 'Weibull characteristic life of the records', value: round(wb ? wb.scale : 1 / lam), unit: 'years', explain: 'Age by which 63 % of comparable sections have failed according to the fitted distribution.' },
        { key: 'b10', label: 'B10 life of the records', value: round(wb ? bLife(wb, 0.1) : -Math.log(0.9) / lam), unit: 'years', explain: 'Age by which 10 % of comparable sections have failed.' },
        { key: 'insLife', label: 'Insulation life at service conditions', value: round(Math.min(lifeIns, 1e6)), unit: 'years', limit: T, sense: 'min', explain: 'Reference life scaled by the Arrhenius factor for the service temperature and the inverse power law for the service stress.' },
        { key: 'ageingRate', label: 'Thermal ageing rate relative to the reference', value: round(afT, 5), unit: '-', explain: 'Rate at which thermal life is consumed compared with continuous operation at the reference temperature; 0.5 means the life is used half as fast.' },
        { key: 'pBreakdown', label: 'Breakdown probability at the operating stress (size effect)', value: round(pBd, 4), unit: '-', explain: 'Weibull weakest-link probability for the whole insulated length at the operating stress, scaled from the test-sample statistics. It shows why long cables need a large ratio of test to service stress.' },
        { key: 'E1pct', label: 'Stress for 1 % breakdown probability of the route', value: round(enl.Eprobability / 1e6), unit: 'kV/mm', explain: 'Stress at which the full length has a 1 % chance of breakdown according to the enlargement law.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Contribution to the failure frequency', categories: ['Internal', 'External', 'Joints', 'Terminations'], series: [{ name: 'Failures per year', y: [rate.internal, rate.external, rate.joint, rate.termination].map((v) => round(v, 4)) }],
          explain: 'Which cause produces the failures. External damage is addressed by burial and routing, joints by reducing their number and by factory testing, internal faults by qualification and after-installation tests.' },
        { type: 'line', title: 'Probability of surviving without a failure', x: { label: 'Years from today', unit: 'yr' }, y: { label: 'Probability', unit: '-' },
          series: [{ name: 'No failure on any circuit', x: yrs, y: yrs.map((t) => Math.exp(-lamSys * t)) }, ...(wb ? [{ name: 'Weibull fit of the records (one section)', x: yrs, y: yrs.map((t) => weibullFunctions(Math.max(t, 1e-9), wb).reliability), dash: true }] : [])],
          explain: 'Reliability function of the cable system for constant rates, with the fitted Weibull curve of the service records for comparison. Read off the probability at the design life.' },
        ...(wb ? [{ type: 'line', title: 'Hazard rate of the service records', x: { label: 'Age', unit: 'yr' }, y: { label: 'Hazard', unit: '1/yr' },
          series: [{ name: 'Weibull hazard', x: yrs.slice(1), y: yrs.slice(1).map((t) => weibullFunctions(t, wb).hazard) }, { name: 'Constant rate used', x: [yrs[1], T], y: [lam, lam], dash: true }],
          explain: 'A falling curve means early-life failures dominate (test and warranty matter most); a rising curve means ageing (plan condition monitoring and replacement). The dashed line is the constant rate of the availability calculation.' }] : []),
        { type: 'bar', title: 'Importance of the basic events', categories: imp.slice(0, 8).map((e) => e.id), series: [{ name: 'Fussell–Vesely importance', y: imp.slice(0, 8).map((e) => round(e.fussellVesely ?? e.criticality, 4)) }],
          explain: 'Fraction of the loss-of-transmission probability that involves each basic event. The tallest bar is where a reduction in rate or repair time pays most.' },
        { type: 'line', title: 'Insulation life against conductor temperature', x: { label: 'Conductor temperature', unit: '°C' }, y: { label: 'log₁₀ life', unit: 'years' },
          series: (() => { const th = linspace(40, 110, 29); return [{ name: 'Electro-thermal life', x: th, y: th.map((t) => Math.log10(combinedLife({ E: inp.E, Eref: inp.Eref, n: inp.nLife, thetaC: t, thetaRefC: inp.thetaRef, lifeRef: inp.lifeRef, activationEnergyEv: inp.Ea }))) }, { name: 'Design life', x: [40, 110], y: [Math.log10(T), Math.log10(T)], dash: true }]; })(),
          explain: 'Arrhenius line at the service stress. Where it crosses the design life is the highest sustained conductor temperature the insulation can live at; each 8–10 K above it halves the life.' },
      ],
      tables: [
        { title: 'Installation, operating and service-life compliance', explain: 'One row per governing criterion carried from the other modules of the suite, with its utilisation. Rows above 1.0 are non-compliant and must be closed before the cable is accepted for the stage shown.',
          columns: [{ key: 'stage', label: 'Stage' }, { key: 'item', label: 'Criterion' }, { key: 'demand', label: 'Demand' }, { key: 'capacity', label: 'Capacity' }, { key: 'unit', label: 'Unit' }, { key: 'uc', label: 'Utilisation' }, { key: 'status', label: 'Status' }],
          rows: comp.map((r) => ({ stage: r.stage, item: r.item, demand: round(r.demand), capacity: round(r.capacity), unit: r.unit, uc: round(r.uc), status: r.uc <= 1 ? 'complies' : 'does not comply' })) },
        { title: 'Failure frequency by cause', explain: 'Rate, share and mean time between events for each cause on one circuit.',
          columns: [{ key: 'cause', label: 'Cause' }, { key: 'rate', label: 'Rate', unit: '1/yr' }, { key: 'share', label: 'Share', unit: '%' }, { key: 'mtbf', label: 'Mean time between', unit: 'years' }],
          rows: CAUSES.map((c) => ({ cause: c, rate: round(rate[c], 4), share: sharePct(rate[c]), mtbf: rate[c] > 0 ? round(1 / rate[c]) : '—' })) },
        { title: 'Minimal cut sets of loss of transmission', explain: 'Smallest combinations of basic events that interrupt transmission, most probable first. Single-event cut sets are single points of failure.',
          columns: [{ key: 'events', label: 'Events' }, { key: 'order', label: 'Order' }, { key: 'p', label: 'Probability' }],
          rows: cut.map((c) => ({ events: c.events.join(' + '), order: c.order, p: Number(c.probability.toExponential(3)) })) },
        ...(fm.byRpn.length ? [{ title: 'Failure modes ranked by risk priority number', explain: 'FMEA rows sorted by RPN = severity × occurrence × detection. Occurrence comes from the failure rate of the class over the design life; act on the top rows first.',
          columns: [{ key: 'rank', label: 'Rank' }, { key: 'item', label: 'Item' }, { key: 'mode', label: 'Failure mode' }, { key: 'S', label: 'Severity' }, { key: 'O', label: 'Occurrence' }, { key: 'D', label: 'Detection' }, { key: 'rpn', label: 'RPN' }],
          rows: fm.byRpn.map((r) => ({ rank: r.rpnRank, item: r.item, mode: r.mode, S: r.severity, O: r.occurrence, D: r.detection, rpn: r.rpn })) }] : []),
      ],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        `Failure frequency ${fmt(lam)} per circuit-year (${sharePct(rate.external)} % external, ${sharePct(rate.internal)} % internal, ${sharePct(rate.joint)} % joints, ${sharePct(rate.termination)} % terminations); mean time to restore ${fmt(mttrDays)} days.`,
        wb ? `Service records: Weibull shape ${fmt(wb.shape)}, characteristic life ${fmt(wb.scale)} years from ${nF} failures and ${recs.length - nF} survivors — ${trend}.` : 'Service records: fewer than three failures (or no convergent fit) — the Weibull trend is not evaluated and a constant rate is shown.',
        `Insulation: ${fmt(lifeIns)} years at ${inp.thetaOp} °C and ${inp.E} kV/mm (thermal ageing rate ${fmt(afT)} × reference).`,
        dual ? `Two circuits: total loss ${fmt(uTop * 100)} % of the time, one circuit out ${fmt(pOne * 100)} % of the time at ${inp.share} % capacity.` : 'Single circuit: every failure interrupts transmission for the full repair time.',
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}. Availability ${fmt(avail * 100, 5)} % against ${inp.availTarget} %; ${fmt(Nlife)} failure(s) expected in ${fmt(T - inp.age0)} years (probability of none ${fmt(pNone)}); insulation life ${fmt(lifeIns)} years.`,
          gov.uc > 1 ? (gov.name.startsWith('Availability') ? `Cut the outage per failure first: ${inp.spare ? 'a repair framework agreement and pre-engineered repair procedure shorten mobilisation' : `holding a spare length removes ${inp.lead} days per failure`}; then reduce the dominant cause (${CAUSES.reduce((a, c) => (rate[c] > rate[a] ? c : a), 'internal')}) — deeper burial or rock for external damage, fewer field joints, or a second circuit.` : gov.name.startsWith('Insulation') ? `Lower the sustained conductor temperature: the life meets ${T} years below about ${fmt(inp.thetaOp - 10 * Math.log2(T / lifeIns))} °C (rule of thumb, 10 K per halving), or reduce the design stress.` : gov.name.startsWith('Fatigue') ? 'Return to the configuration and fatigue modules: reduce curvature ranges at the hot spot (bend stiffener, buoyancy) or plan a mid-life replacement of the dynamic section.' : 'Close the non-compliant row in the compliance table in the module that produced it before accepting the cable.')
            : `Hold the repair preparedness that the availability assumes (${inp.spare ? `${inp.spareL} m of spare cable with jointing kits in store; ${sp.spares} repair set(s) cover the remaining life at the ${inp.serviceLevel} % level` : 'a spare length, which is not held now'}), keep fibre-based monitoring (DTS for hotspots, DAS for anchor activity) to move detection earlier, and update the rates with the as-built burial survey.`,
          'CIGRE TB 815 (2020) / TB 379 (2009); IEC 61025; IEC 60812; DNV-RP-0360 (2016) §8'),
        adv('economic', 'info', `Expected energy not delivered ${fmt(end)} MWh/yr — indicatively ${fmt(end * inp.energyValue / 1e3)} k USD a year at ${inp.energyValue} USD/MWh — plus ${fmt(vesselDaysYear)} repair-spread day(s) a year on average (≈ ${fmt(vesselDaysYear * inp.dayRate / 1e3)} k USD at the entered day rate).`,
          `One avoided failure saves about ${fmt(mttrDays * 24 * inp.P * inp.cf)} MWh and ${fmt(inp.rep)} offshore days. ${inp.spare ? 'The spare already held is worth' : 'A spare length would be worth'} ${fmt(lamSys * inp.lead * 24 * inp.P * inp.cf * (dual ? 1 - inp.share / 100 : 1))} MWh a year in avoided waiting: compare that with its storage and capital cost. Value burial and monitoring options by the change in failure rate × this cost per failure.`,
          'Reliability-centred maintenance (IEC 60300-3-11); risk-based inspection DNV-RP-G101 principles; LCOE conventions (availability losses); whole-life costing ISO 15686-5',
          `Each repair consumes a spare length and joints${spareCu > 0 ? ` (the ${inp.spareL} m spare holds about ${fmt(spareCu)} t of copper ≈ ${fmt(spareCu * CO2_COPPER)} t CO₂e, screening ${CO2_COPPER} t CO₂e/t, for a 3 × 630 mm² cable)` : ''} and a vessel campaign: ${vesselCo2Text(inp.rep)} Lost wind energy of ${fmt(end)} MWh/yr is replaced by other generation — prevention by burial depth is the lower-carbon choice.`),
      ],
      publish: { 'cable.availability': avail, 'cable.failureRate': lamSys, 'cable.energyNotDelivered': end * 3.6e9, 'cable.insulationLife': Math.min(lifeIns, 1e6) },
      quantities: { vessel_days: { cable_repair: round(repairDaysLife) }, ...(inp.spare ? { cable_km: round(inp.spareL / 1e3, 5) } : {}) },
    };
  },
  verify: [
    // One repairable cause only: λ = 0.1 /yr, restore time 36.525 days → μ = 10 /yr. Steady-state unavailability q = λ/(λ + μ)
    // = 0.1/10.1 = 0.0099010 → availability 99.00990 %; expected failures in 30 years 3.0; P(no failure) = e⁻³ = 0.049787;
    // MTBF = 10 years; energy not delivered = 0.009901 × 60 MW × 0.45 × 8766 h = 2343.4 MWh/yr.
    { name: 'Two-state repairable unit: U = λ/(λ + μ), Poisson count and energy not delivered',
      source: 'Rausand & Høyland, System Reliability Theory (2nd ed.) §8.3 and §9.4; hand calculation',
      inputs: { lamInt: 0, lamJoint: 0, lamTerm: 0, lamExt: 0.1, mob: 16.525, wait: 0, rep: 20 },
      expect: [{ kpi: 'availability', value: 99.0099, tol: 1e-5 }, { kpi: 'failuresLife', value: 3, tol: 1e-6 }, { kpi: 'pNoFailure', value: 0.049787, tol: 1e-4 }, { kpi: 'mtbf', value: 10, tol: 1e-6 }, { kpi: 'end', value: 2343.4, tol: 0.001 }] },
    // Arrhenius: L = L0·exp[(Ea/k)(1/T − 1/T0)], Ea = 1.1 eV, k = 8.617333e-5 eV/K → Ea/k = 12764.9 K;
    // T = 353.15 K (80 °C), T0 = 363.15 K: 1/T − 1/T0 = 10/(353.15 × 363.15) = 7.7975e-5 → exponent 0.99534 → factor 2.7056.
    // At the reference stress (E = Eref) the life is 30 × 2.7056 = 81.17 years; the ageing rate is 1/2.7056 = 0.3696.
    // Size effect: P = 1 − exp[−(E/E0)^b·L/L0] = 1 − exp[−(5.9/40)^10 × 75000/10] = 1 − exp(−4.8698e-9 × 7500) = 3.652e-5.
    { name: 'Arrhenius insulation life and Weibull enlargement law',
      source: 'IEC 60216-1 (Arrhenius relation); IEC 62539 (Weibull size effect); hand calculation',
      inputs: { thetaOp: 80, E: 5.9, Eref: 5.9 },
      expect: [{ kpi: 'insLife', value: 81.17, tol: 0.003 }, { kpi: 'ageingRate', value: 0.3696, tol: 0.003 }, { kpi: 'pBreakdown', value: 3.652e-5, tol: 0.01 }] },
    // Two identical circuits without common cause: P(both out) = q², q = 0.0099010 → 9.8030e-5 = 0.0098030 %.
    { name: 'Two independent parallel circuits: P(both out) = q²',
      source: 'Parallel structure of independent repairable components (Rausand & Høyland §9.4); hand calculation',
      inputs: { layout: 'dual', betaCC: 0, lamInt: 0, lamJoint: 0, lamTerm: 0, lamExt: 0.1, mob: 16.525, wait: 0, rep: 20 },
      expect: [{ kpi: 'unavailTotal', value: 0.0098030, tol: 1e-4 }] },
  ],
};
