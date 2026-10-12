// Electromagnetic fields and induced currents: magnetic flux density around the cable (Biot–Savart line conductors),
// distance to exposure thresholds, induced electric field in sea water and fauna, sheath and armour induced currents
// and losses, sheath standing voltage and voltage-limiter duty, induction on a parallel pipeline with AC-corrosion
// screening, the DC field with compass deviation, and the electric field in the insulation.
import { fluxDensity, arrangement, twistAttenuation, fluxDensityProfile, distanceToThreshold, dcFieldWithGeomagnetic, motionalElectricField, faradayElectricField, inducedVoltage, sheathStandingVoltage, svlSelection, acCorrosionCurrentDensity } from '../../engine/electrical/emf.js';
import { temperaturesAtCurrent } from '../../engine/electrical/ampacity.js';
import { acStress, dcField } from '../../engine/electrical/insulation.js';
import { cableInputs, buildElec, adv, check, worst, sev, fmt, round, run, need, linspace, MU0, LIMIT_GENERIC, LIMIT_FE, TB880, CONFIRM, HOURS_YEAR, lossCo2Text } from './_common.js';

export default {
  id: 's13-emf-induced-currents', suite: 13, title: 'Electromagnetic fields and induced currents',
  summary: 'Magnetic field at and above the seabed, induced sheath and armour currents and losses, sheath standing voltage, induction on a neighbouring pipeline and the electric stress in the insulation.',
  covers: [10], fidelity: 'analytical',
  standards: ['ICNIRP (2010) Guidelines for limiting exposure to time-varying electric and magnetic fields (1 Hz–100 kHz)', 'IEC 60287-1-1:2023 §2.3–2.4 (sheath and armour losses)', 'IEEE Std 575-2014 (sheath bonding and induced voltages)', 'EN 50443:2011 / ISO 18086:2019 (AC interference and AC corrosion on pipelines)', 'CIGRE TB 880 (2022) (armour losses)', 'Carson (1926) (earth-return mutual impedance)'],
  models: ["Maxwell's equations in the quasi-static limit: Ampère's law (Biot–Savart superposition of line conductors)", "Faraday's law (induced sheath voltage, induced electric field in a conducting loop)", "Gauss's law (coaxial electric field in the insulation)", 'Motional induction E = v × B', 'Carson earth-return mutual impedance', 'Sheath circulating-current and armour loss factors', 'Helical-lay field attenuation (exponential term)', 'DC resistive field with temperature- and field-dependent conductivity', 'Spreading resistance of a coating defect (AC corrosion)'],
  limits: `${LIMIT_GENERIC} Fields come from straight line-conductor (Biot–Savart) models with a user shielding factor: the screening by lead sheath and steel armour, eddy currents in the armour and the exact twisted-core field need 2-D / 3-D electromagnetic finite elements, which are outside the suite. The lay-twist attenuation keeps only the exponential term and is a screening factor. ${TB880} Exposure and interference limits are overridable defaults to confirm against the governing edition.`,
  publishes: ['cable.Bseabed', 'cable.sheathLossFrac'],
  searchTerms: 'submarine cable magnetic field EMF induced sheath voltage ICNIRP pipeline interference',
  inputs: [
    ...cableInputs(),
    { key: 'arr', label: 'Circuit arrangement', type: 'select', default: 'three-core', group: 'Circuit', role: 'geometry', options: [{ value: 'three-core', label: 'One three-core AC cable' }, { value: 'trefoil', label: 'Three single-core AC cables, trefoil' }, { value: 'flat', label: 'Three single-core AC cables, flat' }, { value: 'bipole', label: 'DC bipole (two single-core cables)' }, { value: 'monopole', label: 'DC monopole (single conductor, remote return)' }] },
    { key: 'spacing', label: 'Axial spacing of single-core cables', unit: 'm', type: 'number', default: 0.5, min: 0.03, max: 100, group: 'Circuit', role: 'geometry', design: true, showIf: { key: 'arr', in: ['trefoil', 'flat', 'bipole'] } },
    { key: 'I', label: 'Current (rms phase or DC pole)', unit: 'A', type: 'number', default: 580, min: 0, max: 6000, group: 'Circuit', role: 'bc', bind: 'cable.loadCurrent', uncertain: { dist: 'normal', cov: 0.05 }, help: 'Energised state: conductor current prescribed as the electrical boundary condition.' },
    { key: 'U', label: 'System voltage (line-to-line, or DC pole voltage)', unit: 'kV', type: 'number', default: 66, min: 1, max: 640, group: 'Circuit', role: 'bc' },
    { key: 'f', label: 'Frequency (AC arrangements)', unit: 'Hz', type: 'number', default: 50, min: 16, max: 60, group: 'Circuit', role: 'bc' },
    { key: 'depth', label: 'Burial depth to top of cable', unit: 'm', type: 'number', default: 1.5, min: 0, max: 30, group: 'Circuit', role: 'ic', bind: 'cable.burialDepth', design: true },
    { key: 'corePitch', label: 'Lay length of the cores (three-core)', unit: 'm', type: 'number', default: 2.2, min: 0.5, max: 10, group: 'Circuit', role: 'geometry', showIf: { key: 'arr', in: ['three-core'] } },
    { key: 'useTwist', label: 'Credit the lay-twist attenuation', type: 'bool', default: false, group: 'Circuit', role: 'criteria', showIf: { key: 'arr', in: ['three-core'] }, help: 'Screening-level exponential factor; leave off for a conservative (untwisted trefoil) field.' },
    { key: 'shield', label: 'Shielding factor of sheath and armour', unit: '-', type: 'number', default: 1.0, min: 0.05, max: 1, group: 'Circuit', role: 'criteria', help: 'Ratio of shielded to unshielded field; 0.8–1.0 at 50 Hz for steel-wire armour from measurement or FE. 1.0 = no credit.' },
    { key: 'hProf', label: 'Height of the field profile above the seabed', unit: 'm', type: 'number', default: 0, min: 0, max: 100, group: 'Circuit', role: 'geometry' },
    { key: 'bond', label: 'Sheath bonding', type: 'select', default: 'both', group: 'Sheath bonding', role: 'bc', options: [{ value: 'both', label: 'Solidly bonded and earthed at both ends' }, { value: 'single', label: 'Single-point bonded' }], help: 'Electrical boundary condition of the metallic sheaths. Submarine sheaths are normally bonded at both ends (and often along the route through semiconductive jackets).' },
    { key: 'Lsec', label: 'Section length for the standing voltage', unit: 'km', type: 'number', default: 1, min: 0.01, max: 400, group: 'Sheath bonding', role: 'geometry', showIf: { key: 'bond', in: ['single'] } },
    { key: 'Ifault', label: 'Through-fault current', unit: 'kA', type: 'number', default: 20, min: 0.1, max: 100, group: 'Sheath bonding', role: 'bc' },
    { key: 'ambient', label: 'Ambient temperature', unit: '°C', type: 'number', default: 12, min: -2, max: 35, group: 'Sheath bonding', role: 'bc', bind: 'site.Tbed' },
    { key: 'T4', label: 'External thermal resistance', unit: 'K·m/W', type: 'number', default: 0.42, min: 0.005, max: 3, group: 'Sheath bonding', role: 'bc', help: 'From the thermal-rating module; sets the temperatures at which the loss factors are evaluated.' },
    { key: 'armourLossFactor', label: 'Armour-loss reduction factor', unit: '-', type: 'number', default: 1.0, min: 0.05, max: 1, group: 'Sheath bonding', role: 'criteria', help: `Multiplier on the IEC 60287 λ2. ${TB880}` },
    { key: 'L', label: 'Route length', unit: 'km', type: 'number', default: 25, min: 0.05, max: 800, group: 'Sheath bonding', role: 'geometry', bind: 'cable.length' },
    { key: 'pipeD', label: 'Distance to a parallel pipeline or cable', unit: 'm', type: 'number', default: 30, min: 0.5, max: 5000, group: 'Neighbouring pipeline', role: 'geometry' },
    { key: 'pipeL', label: 'Length of parallelism', unit: 'km', type: 'number', default: 5, min: 0, max: 400, group: 'Neighbouring pipeline', role: 'geometry', help: '0 = no parallel pipeline.' },
    { key: 'rhoEarth', label: 'Resistivity of the return medium', unit: 'Ω·m', type: 'number', default: 1, min: 0.1, max: 1000, group: 'Neighbouring pipeline', role: 'environment' },
    { key: 'gamma', label: 'Pipeline propagation constant', unit: '1/km', type: 'number', default: 0.5, min: 0.01, max: 10, group: 'Neighbouring pipeline', role: 'material', help: 'Magnitude for the long-line limit: about 0.1–1 per km for well-coated pipelines.' },
    { key: 'rhoHoliday', label: 'Soil resistivity at a coating defect', unit: 'Ω·m', type: 'number', default: 1, min: 0.1, max: 1000, group: 'Neighbouring pipeline', role: 'soil' },
    { key: 'fauna', label: 'Swimming / current speed through the field', unit: 'm/s', type: 'number', default: 1, min: 0, max: 10, group: 'Marine environment', role: 'environment' },
    { key: 'bodyR', label: 'Half-width of an animal (induction loop radius)', unit: 'm', type: 'number', default: 0.15, min: 0.01, max: 3, group: 'Marine environment', role: 'environment' },
    { key: 'heading', label: 'Cable heading from magnetic north (DC)', unit: 'deg', type: 'number', default: 45, min: 0, max: 360, group: 'Marine environment', role: 'geometry', showIf: { key: 'arr', in: ['bipole', 'monopole'] } },
    { key: 'geoB', label: 'Geomagnetic field intensity', unit: 'µT', type: 'number', default: 50, min: 20, max: 70, group: 'Marine environment', role: 'environment' },
    { key: 'dTins', label: 'Temperature drop across the insulation (DC)', unit: 'K', type: 'number', default: 12, min: 0, max: 40, group: 'Marine environment', role: 'bc', showIf: { key: 'arr', in: ['bipole', 'monopole'] } },
    { key: 'Blim', label: 'Magnetic-field reference level', unit: 'µT', type: 'number', default: 200, min: 1, max: 5e5, group: 'Acceptance criteria', role: 'criteria', help: `ICNIRP 2010 general public at 50 Hz: 200 µT (100 µT in the 1998 guidelines); 40 000 µT for static fields (ICNIRP 2009: 400 mT). ${CONFIRM}` },
    { key: 'Vsheath', label: 'Allowable sheath standing voltage', unit: 'V', type: 'number', default: 150, min: 10, max: 1000, group: 'Acceptance criteria', role: 'criteria', help: `Utility practice ranges from 65 V to several hundred volts under full load. ${CONFIRM}` },
    { key: 'Vpipe', label: 'Allowable induced voltage on the pipeline', unit: 'V', type: 'number', default: 60, min: 5, max: 2000, group: 'Acceptance criteria', role: 'criteria', help: `EN 50443 long-term interference limit 60 V (rms). ${CONFIRM}` },
    { key: 'Jlim', label: 'AC current-density limit at a coating defect', unit: 'A/m²', type: 'number', default: 30, min: 1, max: 300, group: 'Acceptance criteria', role: 'criteria', help: `ISO 18086 / EN 15280: corrosion unlikely below about 30 A/m² on a 1 cm² defect. ${CONFIRM}` },
    { key: 'compassLim', label: 'Allowable compass deviation at the profile height (DC)', unit: 'deg', type: 'number', default: 3, min: 0.1, max: 45, group: 'Acceptance criteria', role: 'criteria', showIf: { key: 'arr', in: ['bipole', 'monopole'] } },
  ],
  compute(inp, ctx) {
    const b = buildElec(inp, ctx), g = b.geom, De = g.De;
    const dc = inp.arr === 'bipole' || inp.arr === 'monopole';
    need(!(inp.arr === 'three-core' && b.cores !== 3), 'Incompatible combination: the arrangement is one three-core cable but the cable construction is single-core. Change the number of cores or the arrangement.');
    need(!(inp.arr !== 'three-core' && b.cores === 3), 'Incompatible combination: this arrangement uses single-core cables but the cable construction is three-core. Set the cores to single-core.');
    const s = inp.arr === 'three-core' ? g.s : inp.spacing;
    need(inp.arr === 'three-core' || inp.arr === 'monopole' || s >= De, `The axial spacing (${s} m) is smaller than the cable diameter (${fmt(De)} m): the cables would overlap.`);
    const zc = -(inp.depth + De / 2), kind = inp.arr === 'three-core' ? 'trefoil' : inp.arr;
    const cond = arrangement({ arrangement: kind, I: inp.I, s, z: zc });
    // shielding and (optionally) the lay-twist attenuation evaluated at the seabed distance
    const rSea = Math.abs(zc), tw = inp.arr === 'three-core' && inp.useTwist ? twistAttenuation({ r: rSea + inp.hProf, pitch: inp.corePitch, coreRadius: g.c }) : 1;
    const sf = Math.max(inp.shield * tw, 1e-6);
    const at = (x, z) => fluxDensity({ conductors: cond, point: { x, z }, shieldingFactor: sf }).rms;
    const Bsea = Math.max(at(0, 0), at(s / 2, 0), at(-s / 2, 0));
    const span = Math.max(10, 6 * rSea);
    const prof = fluxDensityProfile({ conductors: cond, height: inp.hProf, xMin: -span, xMax: span, n: 161, shieldingFactor: sf });
    const prof1 = fluxDensityProfile({ conductors: cond, height: inp.hProf + 1, xMin: -span, xMax: span, n: 161, shieldingFactor: sf });
    const dLim = inp.I > 0 ? distanceToThreshold({ conductors: cond, threshold: inp.Blim * 1e-6, shieldingFactor: sf }) : 0;
    const dGeo = inp.I > 0 ? distanceToThreshold({ conductors: cond, threshold: 5e-6, shieldingFactor: sf }) : 0;   // 10 % of a 50 µT background
    const Eind = dc ? motionalElectricField({ velocity: inp.fauna, B: prof.Bmax }) : faradayElectricField({ B: prof.Bmax, f: inp.f, radius: inp.bodyR });
    const Emot = motionalElectricField({ velocity: inp.fauna, B: prof.Bmax });
    // induced currents and losses in sheath and armour at the operating point
    const system = { ...(dc ? { dc: true } : { f: inp.f }), U: inp.U * 1e3, thetaMaxC: 90, armourLossReduction: inp.armourLossFactor };
    const installation = { type: 'custom', T4: inp.T4, ambientC: inp.ambient, ...(inp.arr === 'trefoil' || inp.arr === 'flat' ? { formation: inp.arr, spacing: s, bonding: inp.bond === 'both' ? 'solid' : 'singlePoint' } : {}) };
    let t;
    try { t = temperaturesAtCurrent({ cable: g, installation, system, I: inp.I }); } catch (e) { throw new Error(`The loss calculation did not reach thermal equilibrium at ${inp.I} A: ${e.message}`); }
    const Rs = g.Rs20 != null ? g.Rs20 * (1 + g.alphaS * (t.sheathC - 20)) : null, RA = g.RA20 != null ? g.RA20 * (1 + g.alphaA * (t.armourC - 20)) : null;
    // loss-equivalent currents: Is²·Rs = λ1·I²·R per core; IA²·RA = n·λ2·I²·R for the armour
    const Is = Rs ? inp.I * Math.sqrt(Math.max(t.lambda1, 0) * t.R / Rs) : 0, IA = RA ? inp.I * Math.sqrt(b.cores * Math.max(t.lambda2, 0) * t.R / RA) : 0;
    // sheath standing voltage (open sheath circuit) under load and under a through fault on one phase
    const ds = g.sheathMeanDiameter;
    const ssv = !dc && ds ? sheathStandingVoltage({ conductors: cond, sheathMeanDiameter: ds, f: inp.f, length: inp.Lsec * 1e3 }) : { max: 0, voltage: 0 };
    const faultCond = cond.map((c, i) => ({ x: c.x, z: c.z, I: i === 0 ? inp.Ifault * 1e3 : 0 }));
    const ssvF = !dc && ds ? sheathStandingVoltage({ conductors: faultCond, sheathMeanDiameter: ds, f: inp.f, length: inp.Lsec * 1e3 }) : { max: 0 };
    const svl = !dc && ds && inp.bond === 'single' ? svlSelection({ faultGradient: ssvF.max, length: inp.Lsec * 1e3, ratings: [1e3, 3e3, 6e3, 9e3, 12e3, 18e3, 24e3] }) : null;
    // parallel pipeline: balanced load (only the differences in distance remain) and the single-phase fault
    const victims = cond.map((c) => ({ d: Math.max(Math.hypot(inp.pipeD - c.x, c.z - zc), 0.05), I: c.I }));
    // Carson's simplified mutual impedance holds inside the equivalent earth-return depth De = 658.87·√(ρ/f); beyond it the coupling is negligible
    const DeCarson = 658.87 * Math.sqrt(inp.rhoEarth / (dc ? 50 : inp.f)), farPipe = inp.pipeL > 0 && !dc && (inp.pipeD >= DeCarson || victims.some((v) => v.d >= DeCarson));
    const hasPipe = inp.pipeL > 0 && !dc && !farPipe;
    const ind = hasPipe ? inducedVoltage({ conductors: victims, f: inp.f, rhoEarth: inp.rhoEarth, length: inp.pipeL * 1e3, gamma: inp.gamma / 1e3 }) : { E: 0, voltage: 0, voltageInsulated: 0 };
    const indF = hasPipe ? inducedVoltage({ conductors: [{ d: inp.pipeD, I: inp.Ifault * 1e3 }], f: inp.f, rhoEarth: inp.rhoEarth, length: inp.pipeL * 1e3, gamma: inp.gamma / 1e3 }) : { voltage: 0 };
    const acc = hasPipe ? acCorrosionCurrentDensity({ Vac: ind.voltage, rhoSoil: inp.rhoHoliday, limit: inp.Jlim }) : { J: 0, utilisation: 0, VacAllowable: Infinity };
    // DC: compass deviation and resistive field in the insulation
    const comp = dc ? dcFieldWithGeomagnetic({ conductors: cond.map((c) => ({ x: c.x, z: c.z, I: typeof c.I === 'number' ? c.I : c.I.re })), point: { x: 0, z: inp.hProf }, headingDeg: inp.heading, geomagneticTotal: inp.geoB * 1e-6 }) : null;
    const ri = g.dScreen / 2, ro = g.Di / 2, U0 = dc ? inp.U * 1e3 : inp.U * 1e3 / Math.sqrt(3);
    const rr = linspace(ri, ro, 40);
    const dcf = dc ? run(() => dcField({ U: U0, ri, ro, alpha: 0.08, beta: 0.05e-6, deltaT: inp.dTins, n: 81 }), 'DC field in the insulation') : null;
    const Eac = rr.map((r) => acStress({ U0, ri, ro, r }).E / 1e6);
    const checks = [
      check('Magnetic flux density at the seabed surface', Bsea * 1e6, inp.Blim, 'µT', dc ? 'ICNIRP (2009) static magnetic fields' : 'ICNIRP (2010) Table 4 (general-public reference level)', 'Largest rms flux density on the seabed directly above the cable against the exposure reference level. It falls with the square of distance for a balanced AC circuit.'),
      ...(!dc && inp.bond === 'single' && ds ? [check('Sheath standing voltage under load', ssv.voltage, inp.Vsheath, 'V', 'IEEE Std 575-2014 (induced sheath voltage); IEC 60287-1-1 §2.3.1 (X·I)', 'Voltage at the open end of a single-point-bonded sheath section at the load current. Above the limit the jacket and link boxes are overstressed and touch safety is compromised.')] : []),
      ...(hasPipe ? [
        check('Induced voltage on the parallel pipeline (steady state)', ind.voltage, inp.Vpipe, 'V', 'EN 50443:2011 (long-term interference); Carson mutual impedance', 'Longitudinal voltage induced on the neighbouring pipeline by the load currents. Above the limit it is a touch hazard at above-water appurtenances and drives AC corrosion.'),
        check('AC corrosion current density at a coating defect', acc.J, inp.Jlim, 'A/m²', 'ISO 18086:2019 / EN 15280 (AC corrosion likelihood)', 'AC current density through a 1 cm² coating holiday for the induced voltage. Above about 30 A/m² AC corrosion becomes likely despite cathodic protection.'),
      ] : []),
      ...(dc ? [check('Compass deviation at the profile height', Math.abs(comp.compassDeviationDeg), inp.compassLim, 'deg', 'Navigation authority criterion; Biot–Savart field superposed on the geomagnetic field', 'Rotation of the horizontal magnetic field direction caused by the DC cable. Large deviations affect magnetic compasses and migratory species near the cable.'),
        check('DC stress at the insulation screen (stress inversion)', dcf.Eouter / 1e6, Math.max(dcf.Elaplace / 1e6, 1e-6), 'kV/mm', 'Eoll (1975); Jeroense & Morshuis (1998); CIGRE TB 496', 'Loaded DC stress at the insulation screen against the no-load stress at the conductor screen. Above 1.0 the field has inverted and the outer insulation and accessories see more than the nominal design stress.')] : []),
    ];
    const gov = worst(checks);
    const sheathLossKw = (t.losses.sheath + t.losses.armour) * inp.L, sheathMwh = sheathLossKw * HOURS_YEAR / 1e3 * 0.4;   // indicative loss-load factor 0.4
    const heights = linspace(0, Math.max(20, 3 * rSea), 40);
    const Bh = heights.map((h) => at(0, h) * 1e6);
    return {
      kpis: [
        { key: 'Bseabed', label: 'Flux density at the seabed above the cable', value: round(Bsea * 1e6, 5), unit: 'µT', limit: inp.Blim, sense: 'max', explain: 'Largest rms magnetic flux density on the seabed surface over the cable. Compare with the exposure reference level and with the geomagnetic background of 25–65 µT.' },
        { key: 'Bprofile', label: 'Peak flux density at the profile height', value: round(prof.Bmax * 1e6, 5), unit: 'µT', explain: 'Maximum of the lateral profile at the chosen height above the seabed.' },
        { key: 'dLimit', label: 'Distance from the cable to the reference level', value: round(dLim), unit: 'm', explain: 'Distance from the cable centre (vertically) beyond which the field is below the reference level; 0 means it is below everywhere outside the cable group.' },
        { key: 'd5uT', label: 'Distance to 5 µT', value: round(dGeo), unit: 'm', explain: 'Distance at which the cable field falls to about a tenth of the geomagnetic background — an indicative range of detectability by magneto-sensitive species.' },
        { key: 'Einduced', label: 'Induced electric field in an animal at the profile height', value: round(Eind * 1e6, 4), unit: 'µV/m', explain: dc ? 'Motional field v × B for an animal or water moving through the static field.' : 'Faraday induction π·f·a·B in a body of the given half-width; electro-sensitive fish respond to fields of the order of µV/m.' },
        { key: 'Emotional', label: 'Motional electric field in moving sea water', value: round(Emot * 1e6, 4), unit: 'µV/m', explain: 'E = v·B for water or an animal moving at the entered speed across the field.' },
        { key: 'lambda1', label: 'Sheath loss factor λ1', value: round(t.lambda1), unit: '-', explain: 'Sheath loss over conductor loss from circulating and eddy currents induced by the conductor currents (Faraday).' },
        { key: 'lambda2', label: 'Armour loss factor λ2 (as used)', value: round(t.lambda2), unit: '-', explain: 'Armour loss over conductor loss after the reduction factor; the IEC value is known to be high for three-core cables.' },
        { key: 'Isheath', label: 'Loss-equivalent sheath current', value: round(Is), unit: 'A', explain: 'Sheath current that would dissipate the calculated sheath loss in the sheath resistance: the induced circulating current for solid bonding.' },
        { key: 'Iarmour', label: 'Loss-equivalent armour current', value: round(IA), unit: 'A', explain: 'Equivalent current for the armour loss (the true loss is hysteresis and eddy currents in the wires, not a net current, in a balanced three-core cable).' },
        { key: 'Wsheath', label: 'Sheath and armour losses', value: round(t.losses.sheath + t.losses.armour), unit: 'W/m', explain: 'Heat generated in the metallic sheaths and armour of the cable at this current: induced-current losses that add to conductor heating.' },
        { key: 'gradLoad', label: 'Sheath voltage gradient under load', value: round(ssv.max * 1e3, 5), unit: 'V/km', explain: 'Open-circuit voltage induced per kilometre of sheath at the load current (ω·μ0/2π·I·ln(2s/d) for a trefoil).' },
        { key: 'gradFault', label: 'Sheath voltage gradient under through fault', value: round(ssvF.max * 1e3), unit: 'V/km', explain: 'Gradient with the fault current in one phase: the duty of the sheath voltage limiters.' },
        { key: 'Vstanding', label: 'Sheath standing voltage of the section', value: round(inp.bond === 'single' ? ssv.voltage : 0), unit: 'V', explain: 'Standing voltage at the open end of a single-point-bonded section (0 when both ends are bonded: the voltage drives a circulating current instead).' },
        { key: 'Vpipeline', label: 'Induced voltage on the pipeline (load)', value: round(ind.voltage, 4), unit: 'V', explain: 'Longitudinal voltage on the parallel pipeline in normal operation; small for a balanced circuit because the three phase contributions nearly cancel.' },
        { key: 'VpipeFault', label: 'Induced voltage on the pipeline (earth fault)', value: round(indF.voltage), unit: 'V', explain: 'Voltage during a single-phase fault with earth return, before any credit for sheath screening: the short-term touch and coating-stress case.' },
        { key: 'Jac', label: 'AC current density at a coating defect', value: round(acc.J, 4), unit: 'A/m²', explain: 'Current density through a 1 cm² holiday at the induced steady-state voltage.' },
        { key: 'Emax', label: 'Electric stress at the conductor screen', value: round((dc ? dcf.Einner : acStress({ U0, ri, ro }).Emax) / 1e6), unit: 'kV/mm', explain: 'Electric field in the insulation next to the conductor screen. The field is confined inside the earthed screen: there is no external electric field from a screened cable.' },
        ...(dc ? [{ key: 'compass', label: 'Compass deviation at the profile height', value: round(comp.compassDeviationDeg), unit: 'deg', limit: inp.compassLim, sense: 'max', explain: 'Change of magnetic heading above the cable; it depends strongly on cable heading and conductor spacing.' },
          { key: 'Eouter', label: 'DC stress at the insulation screen (loaded)', value: round(dcf.Eouter / 1e6), unit: 'kV/mm', explain: 'Resistive field at the outer screen under load; when it exceeds the inner value the stress has inverted.' }] : []),
      ],
      checks,
      plots: [
        { type: 'line', title: 'Magnetic flux density across the route', x: { label: 'Horizontal distance from the cable', unit: 'm' }, y: { label: 'Flux density (rms)', unit: 'µT' },
          series: [{ name: `At ${inp.hProf} m above the seabed`, x: Array.from(prof.x), y: Array.from(prof.B, (v) => v * 1e6) }, { name: `At ${inp.hProf + 1} m`, x: Array.from(prof1.x), y: Array.from(prof1.B, (v) => v * 1e6) }],
          explain: 'Lateral profile of the field. The peak sits over the cable and the width at half height is about the burial depth: deeper burial lowers and spreads the peak.' },
        { type: 'line', title: 'Flux density above the cable', x: { label: 'Height above the seabed', unit: 'm' }, y: { label: 'log₁₀ flux density', unit: 'µT' },
          series: [{ name: 'On the cable axis', x: heights, y: Bh.map((v) => Math.log10(Math.max(v, 1e-9))) }, { name: 'Reference level', x: [0, heights[39]], y: [Math.log10(inp.Blim), Math.log10(inp.Blim)], dash: true }, { name: 'Geomagnetic background (50 µT)', x: [0, heights[39]], y: [Math.log10(50), Math.log10(50)], dash: true }],
          explain: 'Decay with height: 1/r² for a balanced three-phase group or a bipole, 1/r for a single DC conductor. Where the curve crosses the background line the cable is no longer the dominant field.' },
        { type: 'line', title: 'Electric field in the insulation', x: { label: 'Radius', unit: 'mm' }, y: { label: 'Field', unit: 'kV/mm' },
          series: dc ? [{ name: 'DC, loaded (resistive)', x: Array.from(dcf.r, (r) => r * 1e3), y: Array.from(dcf.E, (e) => e / 1e6) }, { name: 'No load (capacitive)', x: rr.map((r) => r * 1e3), y: Eac }] : [{ name: 'AC (capacitive)', x: rr.map((r) => r * 1e3), y: Eac }],
          explain: dc ? 'Under load the warm inner insulation conducts better and the DC field shifts outwards; a rising curve means stress inversion. The no-load curve is the Laplace field.' : 'Coaxial field falling as 1/r from the conductor screen to the insulation screen; outside the earthed screen the electric field is zero.' },
        { type: 'bar', title: 'Loss split in the metallic layers', categories: ['Conductor', 'Sheath', 'Armour'], series: [{ name: 'W/m', y: [t.losses.conductor, t.losses.sheath, t.losses.armour].map((v) => round(v)) }],
          explain: 'Induced-current losses in sheath and armour beside the conductor loss. They are the price of bonding both ends (no standing voltage) and of magnetic armour.' },
      ],
      tables: [{ title: 'Conductor positions and currents', explain: 'Line conductors used in the Biot–Savart superposition (x horizontal, z relative to the seabed, negative when buried) with their phasor currents.',
        columns: [{ key: 'i', label: 'Conductor' }, { key: 'x', label: 'x', unit: 'm' }, { key: 'z', label: 'z', unit: 'm' }, { key: 'mag', label: 'Current', unit: 'A' }, { key: 'ang', label: 'Phase', unit: 'deg' }],
        rows: cond.map((c, i) => { const I = typeof c.I === 'number' ? { re: c.I, im: 0 } : c.I; return { i: i + 1, x: round(c.x), z: round(c.z), mag: round(Math.hypot(I.re, I.im)), ang: round(Math.atan2(I.im, I.re) * 180 / Math.PI) }; }) }],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        `Field model: ${cond.length} line conductor(s), spacing ${fmt(s * 1e3)} mm, axis ${fmt(-zc)} m below the seabed; shielding ${inp.shield}${tw < 1 ? `, twist attenuation ${fmt(tw)}` : ''}.`,
        ...(farPipe ? [`The parallel pipeline at ${inp.pipeD} m lies beyond the equivalent earth-return depth of ${fmt(DeCarson)} m for a return medium of ${inp.rhoEarth} Ω·m: inductive coupling is negligible and the pipeline checks are not evaluated.`] : []),
        inp.bond === 'both' ? `Sheaths bonded at both ends: circulating current ≈ ${fmt(Is)} A per sheath, no standing voltage.` : `Single-point bonding: standing voltage ${fmt(ssv.voltage)} V over ${inp.Lsec} km${svl ? `; sheath voltage limiter ≥ ${fmt(svl.required / 1e3)} kV for the through fault${svl.selected ? ` (select ${svl.selected / 1e3} kV)` : ''}` : ''}.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name}: ${fmt(gov.uc)}. Field ${fmt(Bsea * 1e6)} µT at the seabed, below ${inp.Blim} µT beyond ${fmt(dLim)} m; sheath and armour losses ${fmt(t.losses.sheath + t.losses.armour)} W/m (λ1 ${fmt(t.lambda1)}, λ2 ${fmt(t.lambda2)}).`,
          gov.uc > 1 ? (gov.name.startsWith('Magnetic') ? 'Bury deeper, reduce the phase spacing (bundle the cables) or use a three-core cable: the field of a balanced group scales with spacing / distance².' : gov.name.startsWith('Sheath') ? `Shorten the single-point-bonded section to ${fmt(inp.Lsec / gov.uc)} km or bond both ends and accept the circulating-current loss.` : gov.name.includes('pipeline') || gov.name.includes('corrosion') ? `Increase the separation, shorten the parallelism or install gradient-control earthing / decoupling devices on the pipeline; the induced voltage must fall below ${fmt(Math.min(inp.Vpipe, acc.VacAllowable))} V.` : 'Reduce the conductor spacing (bundle the poles) or change the cable heading; review the DC insulation design for the loaded stress.')
            : 'State the field at the seabed and at 1 m above it in the environmental statement with the burial depth assumed, and confirm the sheath and armour loss basis by measurement or 2-D/3-D FE where the rating depends on it.',
          'ICNIRP (2010); IEEE Std 575-2014; IEC 60287-1-1:2023 §2.3–2.4; EN 50443:2011; ISO 18086:2019'),
        adv('economic', 'info', `Induced losses in sheath and armour are ${fmt(t.losses.sheath + t.losses.armour)} W/m — ${fmt(sheathLossKw)} kW over ${inp.L} km, of the order of ${fmt(sheathMwh)} MWh a year (indicative loss-load factor 0.4).`,
          `These losses are ${fmt((t.lambda1 + t.lambda2) / (1 + t.lambda1 + t.lambda2) * 100)} % of the current-dependent loss. A justified armour-loss reduction factor (measurement or FE) usually removes most of the armour share on paper and can release rating without more copper; non-magnetic or reduced-coverage armour removes it physically at higher cable cost.`,
          'IEC 60287-3-2 (capitalised losses); whole-life costing ISO 15686-5',
          `${lossCo2Text(sheathMwh)} Deeper burial to reduce the seabed field costs trenching energy and rating; a three-core design achieves a low field with no extra seabed footprint.`),
      ],
      publish: { 'cable.Bseabed': Bsea, 'cable.sheathLossFrac': (t.lambda1 + t.lambda2) / (1 + t.lambda1 + t.lambda2) },
      quantities: {},
    };
  },
  verify: [
    // Single DC conductor (Ampère / Biot–Savart): B = μ0·I/(2πr). Single-core cable diameter 81.454 mm, top 1.5 m below the seabed
    // → axis at 1.5407 m; I = 1000 A → B = 2e-7 × 1000 / 1.5407 = 129.81 µT at the seabed.
    { name: 'Long straight conductor: B = μ0·I/(2π·r)',
      source: "Ampère's law / Biot–Savart (any electromagnetics text); hand calculation",
      inputs: { cores: '1', arr: 'monopole', I: 1000, depth: 1.5, Blim: 40000 },
      expect: [{ kpi: 'Bseabed', value: 129.81, tol: 0.002 }] },
    // Balanced three-phase flat formation, far field (r ≫ s): B ≈ √3·μ0·I·s/(2π·r²). s = 0.5 m, I = 1000 A, point 10 m above the
    // centre conductor (profile height 8.4593 m + 1.5407 m): B = 1.7321 × 2e-7 × 1000 × 0.5 / 100 = 1.732 µT (within 1 % at r/s = 20).
    // Trefoil sheath gradient check via the flat arrangement is avoided: ω·μ0/(2π)·I·ln(2s/d) applies to a trefoil.
    { name: 'Flat three-phase group, far-field dipole limit B = √3·μ0·I·s/(2π·r²)',
      source: 'Biot–Savart superposition, far-field expansion for a balanced flat formation (e.g. CIGRE TB 373 Annex)',
      inputs: { cores: '1', arr: 'flat', spacing: 0.5, I: 1000, depth: 1.5, hProf: 8.4593, pipeL: 0 },
      expect: [{ kpi: 'Bprofile', value: 1.732, tol: 0.01 }] },
    // Trefoil, single-point bonded: sheath gradient E = ω·μ0/(2π)·I·ln(2s/d). s = 0.5 m, sheath mean diameter d = 51.854 + 1.8 = 53.654 mm,
    // I = 1000 A, 50 Hz: 314.159 × 2e-7 × 1000 × ln(18.638) = 0.062832 × 2.9252 = 0.18380 V/m = 183.8 V/km.
    { name: 'Trefoil sheath standing voltage gradient E = ω·μ0/(2π)·I·ln(2s/d)',
      source: 'IEC 60287-1-1 §2.3.1 (X = 2ω·10⁻⁷·ln(2s/d)); IEEE Std 575 Annex',
      inputs: { cores: '1', arr: 'trefoil', spacing: 0.5, I: 1000, bond: 'single', pipeL: 0 },
      expect: [{ kpi: 'gradLoad', value: 183.8, tol: 0.003 }] },
  ],
};
