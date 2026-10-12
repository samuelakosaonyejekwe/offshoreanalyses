// Service-life and failure assessment: combined damage budget, reliability-based fatigue life (lognormal format),
// armour-wire fretting wear and corrosion, end-of-life tension capacity, functional redundancy of the lines, and the
// integrity / design-compliance summary of the utilisations produced by the other modules.
import { combineDamage, fatigueReliabilityLognormal, medianLifeFromCharacteristic } from '../../engine/fatigue/damage.js';
import { archardWear, wireSlipPerCycle } from '../../engine/xsection/flexfatigue.js';
import { corrosionWastage } from '../../engine/fatigue/degradation.js';
import { adv, check, worst, sev, fmt, round, run, need, linspace, DEG, CO2_STEEL, STD_UMB, vesselCo2Text } from './_common.js';

// probability that at least k of n identical independent lines are still working when each has failed with probability q
function kOfN(n, k, q) {
  let p = 0;
  const C = (a, b) => { let c = 1; for (let i = 1; i <= b; i++) c = c * (a - b + i) / i; return c; };
  for (let j = k; j <= n; j++) p += C(n, j) * (1 - q) ** j * q ** (n - j);
  return p;
}

export default {
  id: 's12-service-life', suite: 12, title: 'Service-life and failure assessment',
  summary: 'Remaining life, failure probability over time, wear and corrosion margins, functional redundancy and a compliance summary of all integrity checks.',
  covers: [15], fidelity: 'reduced-order',
  standards: [`${STD_UMB} §7.4 and Annex (service life, design fatigue factor)`, 'DNV-RP-C203 (2021) §2.2 and Appendix F (Miner sum, fatigue reliability)', 'Wirsching (1984) J. Struct. Eng. 110(10) (lognormal fatigue reliability)', 'DNV-RP-G101 (2021) (risk-based inspection)', 'ISO 14224 (2016) (reliability data, failure rates)'],
  models: ['Palmgren–Miner damage accumulation over several sources with design fatigue factors', 'Reliability-based service life (lognormal fatigue format)', 'Archard wear model for fretting between armour wires', 'Linear corrosion wastage after coating breakdown', 'k-out-of-n functional redundancy with constant failure rates', 'Failure-risk indicators and compliance ranking'],
  limits: 'Damage sources are added linearly without sequence effects; the reliability model is a single-slope lognormal format with generic uncertainty values, so the failure probabilities are notional (for ranking and inspection planning, not actuarial). Wear uses one contact pressure and slip range; fretting-fatigue interaction and hydrogen effects are not modelled. Failure rates of functional lines must come from the operator database.',
  publishes: ['life.fatigueYears', 'life.designYears', 'umbilical.serviceLifeYears', 'umbilical.complianceMax'],
  searchTerms: 'umbilical service life assessment reliability life extension integrity',
  inputs: [
    { key: 'life', label: 'Design service life', unit: 'years', type: 'number', default: 25, min: 1, max: 60, group: 'Service conditions', role: 'criteria', bind: 'life.designYears' },
    { key: 'age', label: 'Present age', unit: 'years', type: 'number', default: 0, min: 0, max: 60, group: 'Service conditions', role: 'ic', help: 'Years already in service (0 for a new design). Sets the starting point of the remaining-life estimate.' },
    { key: 'sources', label: 'Damage sources', type: 'table', group: 'Fatigue damage budget', role: 'load',
      columns: [{ key: 'name', label: 'Source', type: 'text' }, { key: 'rate', label: 'Damage per year', unit: '1/yr', type: 'number' }, { key: 'once', label: 'One-off damage', unit: '-', type: 'number' }, { key: 'dff', label: 'DFF', unit: '-', type: 'number' }],
      default: [{ name: 'Wave-frequency tension and bending', rate: 0.0025, once: 0, dff: 10 }, { name: 'Vortex-induced vibration', rate: 0.0003, once: 0, dff: 10 }, { name: 'Pressure cycles', rate: 0.0002, once: 0, dff: 10 }, { name: 'Reeling, load-out and installation', rate: 0, once: 0.004, dff: 10 }],
      help: 'Characteristic Miner damage of the governing element from each source (wave fatigue, VIV, pressure, installation), without the design fatigue factor; one-off damage is consumed before service starts.' },
    { key: 'snM', label: 'S–N inverse slope of the governing detail', unit: '-', type: 'number', default: 3, min: 2.5, max: 6, group: 'Reliability model', role: 'material' },
    { key: 'sLogN', label: 'Standard deviation of log10 N', unit: '-', type: 'number', default: 0.2, min: 0.05, max: 0.5, group: 'Reliability model', role: 'material', help: '0.20 for the DNV-RP-C203 design curves (mean minus two standard deviations).' },
    { key: 'deltaCov', label: 'CoV of the Miner sum at failure', unit: '-', type: 'number', default: 0.3, min: 0, max: 1, group: 'Reliability model', role: 'material' },
    { key: 'bCov', label: 'CoV of the stress calculation', unit: '-', type: 'number', default: 0.2, min: 0, max: 0.6, group: 'Reliability model', role: 'load', help: 'Uncertainty of the load spectrum and stress model (0.15–0.30 typical).' },
    { key: 'pfTarget', label: 'Target annual failure probability', unit: '-', type: 'number', default: 1e-4, min: 1e-6, max: 1e-2, group: 'Acceptance criteria', role: 'criteria', help: '10⁻⁴ per year for components whose failure has serious but not catastrophic consequences (DNV safety-class convention); a project decision.' },
    { key: 'kWear', label: 'Wear coefficient', unit: '10⁻¹⁶ m²/N', type: 'number', default: 2, min: 0, max: 1000, group: 'Armour wear and corrosion', role: 'material', uncertain: { dist: 'lognormal', cov: 0.6 }, help: 'Archard coefficient of galvanised wire on wire, lubricated by sea water or filling compound: order 10⁻¹⁶–10⁻¹⁵ m²/N; calibrate against fretting tests.' },
    { key: 'pContact', label: 'Contact pressure on the wire footprint', unit: 'MPa', type: 'number', default: 15, min: 0, max: 500, group: 'Armour wear and corrosion', role: 'load' },
    { key: 'R', label: 'Armour pitch-circle radius', unit: 'mm', type: 'number', default: 59, min: 5, max: 300, group: 'Armour wear and corrosion', role: 'geometry' },
    { key: 'alpha', label: 'Armour lay angle', unit: 'deg', type: 'number', default: 15, min: 3, max: 55, group: 'Armour wear and corrosion', role: 'geometry' },
    { key: 'dWire', label: 'Armour wire diameter', unit: 'mm', type: 'number', default: 5, min: 1, max: 12, group: 'Armour wear and corrosion', role: 'geometry' },
    { key: 'dKappa', label: 'Mean curvature range per cycle', unit: '1/m', type: 'number', default: 0.002, min: 0, max: 0.2, group: 'Armour wear and corrosion', role: 'load', bind: 'umbilical.curvatureRange', help: 'Long-term average curvature range at the hang-off (not the storm value).' },
    { key: 'cyclesYr', label: 'Bending cycles per year', unit: '-', type: 'number', default: 5e6, min: 0, max: 5e7, group: 'Armour wear and corrosion', role: 'load' },
    { key: 'corrRate', label: 'Corrosion rate after coating breakdown', unit: 'mm/yr', type: 'number', default: 0.01, min: 0, max: 1, group: 'Armour wear and corrosion', role: 'environment', uncertain: { dist: 'lognormal', cov: 0.5 }, help: 'Radial wastage of the wire in a flooded, oxygen-depleted annulus is low (≈ 0.01–0.1 mm/yr); a breached sheath in the splash zone is much higher.' },
    { key: 'coatLife', label: 'Galvanising / sheath protection life', unit: 'years', type: 'number', default: 15, min: 0, max: 60, group: 'Armour wear and corrosion', role: 'ic' },
    { key: 'lossAllow', label: 'Allowable loss of wire section', unit: '%', type: 'number', default: 20, min: 1, max: 60, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'Tmax', label: 'Maximum tension in service', unit: 'kN', type: 'number', default: 120, min: 0, max: 5000, group: 'End-of-life strength', role: 'bc', bind: 'umbilical.topTension' },
    { key: 'Tcap', label: 'Allowable tension when new', unit: 'kN', type: 'number', default: 320, min: 1, max: 20000, group: 'End-of-life strength', role: 'criteria', bind: 'umbilical.allowTension' },
    { key: 'armShare', label: 'Share of tension carried by the armour', unit: '-', type: 'number', default: 0.7, min: 0, max: 1, group: 'End-of-life strength', role: 'material', help: 'From the cross-section module; only this share of the capacity degrades with wire loss.' },
    { key: 'nLines', label: 'Functional lines of one service (incl. spares)', type: 'int', default: 5, min: 1, max: 40, group: 'Functional redundancy', role: 'geometry', help: 'E.g. hydraulic tubes of the same service laid up in the umbilical.' },
    { key: 'nReq', label: 'Lines required for operation', type: 'int', default: 4, min: 1, max: 40, group: 'Functional redundancy', role: 'criteria' },
    { key: 'lambda', label: 'Failure rate per line', unit: '1/yr', type: 'number', default: 0.003, min: 0, max: 1, group: 'Functional redundancy', role: 'load', uncertain: { dist: 'lognormal', cov: 0.5 }, help: 'Leak / blockage / loss of insulation rate per line from the operator reliability database (ISO 14224 format).' },
    { key: 'relTarget', label: 'Required probability of retaining the function', unit: '-', type: 'number', default: 0.95, min: 0.5, max: 0.9999, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'ucStress', label: 'Tension–bending stress utilisation', unit: '-', type: 'number', default: 0.88, min: 0, max: 20, group: 'Design-compliance summary', role: 'load', bind: 'umbilical.stressUtilisation' },
    { key: 'ucTube', label: 'Tube pressure-containment utilisation', unit: '-', type: 'number', default: 0.94, min: 0, max: 20, group: 'Design-compliance summary', role: 'load', bind: 'umbilical.tubeUtilisation' },
    { key: 'ucCollapse', label: 'Collapse utilisation', unit: '-', type: 'number', default: 0.06, min: 0, max: 20, group: 'Design-compliance summary', role: 'load', bind: 'umbilical.collapseUtilisation' },
    { key: 'mbrAllow', label: 'Dynamic minimum bend radius (allowable)', unit: 'm', type: 'number', default: 10.3, min: 0.1, max: 100, group: 'Design-compliance summary', role: 'criteria', bind: 'umbilical.mbrDynamic' },
    { key: 'mbrActual', label: 'Tightest operating bend radius', unit: 'm', type: 'number', default: 16.6, min: 0.1, max: 1e5, group: 'Design-compliance summary', role: 'load', bind: 'umbilical.minBendRadius' },
    { key: 'dayRate', label: 'Indicative replacement-campaign day rate', unit: 'USD', type: 'number', default: 180000, min: 0, max: 2e6, group: 'Economics', role: 'cost', help: 'Indicative unit rate, used only to express the value of deferring a replacement.' },
    { key: 'campaignDays', label: 'Replacement campaign duration', unit: 'days', type: 'number', default: 12, min: 1, max: 120, group: 'Economics', role: 'operation' },
    { key: 'steelT', label: 'Steel in the dynamic section', unit: 't', type: 'number', default: 23, min: 0, max: 5000, group: 'Economics', role: 'operation', bind: 'quantities.steel_t' },
  ],
  compute(inp) {
    const src = (inp.sources || []).filter((s) => (s.rate || 0) > 0 || (s.once || 0) > 0);
    need(src.length > 0, 'Enter at least one damage source with a yearly or one-off damage.');
    for (const s of src) need(s.rate >= 0 && s.once >= 0 && s.dff >= 1, `Damage source '${s.name}': damage cannot be negative and the DFF must be at least 1.`);
    need(inp.nReq <= inp.nLines, `Over-constrained redundancy: ${inp.nReq} lines are required but only ${inp.nLines} exist.`);
    const L = inp.life;
    const comb = run(() => combineDamage(src.map((s) => ({ name: s.name, damage: s.once + s.rate * L, dff: s.dff }))), 'Damage combination');
    const rate = src.reduce((a, s) => a + s.rate, 0), once = src.reduce((a, s) => a + s.once, 0);
    const rateF = src.reduce((a, s) => a + s.rate * s.dff, 0), onceF = src.reduce((a, s) => a + s.once * s.dff, 0);
    // fatigue design life: factored damage reaches 1
    const lifeDff = rateF > 0 ? Math.max(0, (1 - onceF) / rateF) : 1e6;
    const remaining = Math.max(0, lifeDff - inp.age);
    // reliability: median life from the characteristic (design-curve) life, lognormal format
    const charLife = rate > 0 ? Math.max(1e-3, (1 - once) / rate) : 1e9;
    const med = medianLifeFromCharacteristic(charLife, { sLogN: inp.sLogN });
    const tMax = Math.max(L * 2, 10), times = linspace(1, tMax, Math.min(120, Math.max(20, Math.round(tMax))));
    const rel = run(() => fatigueReliabilityLognormal({ medianLife: med.medianLife, m: inp.snM, times, deltaCov: inp.deltaCov, aCov: med.aCov, bCov: inp.bCov }), 'Fatigue reliability');
    const atL = run(() => fatigueReliabilityLognormal({ medianLife: med.medianLife, m: inp.snM, times: [L], deltaCov: inp.deltaCov, aCov: med.aCov, bCov: inp.bCov }), 'Fatigue reliability');
    let relLife = tMax; for (let i = 0; i < times.length; i++) if (rel.annualPf[i] > inp.pfTarget) { relLife = i ? times[i - 1] : 0; break; }
    // armour wear (Archard) and corrosion wastage of the wire
    const slip = wireSlipPerCycle({ R: inp.R / 1e3, alpha: inp.alpha * DEG, kappaRange: inp.dKappa });
    const d0 = inp.dWire / 1e3;
    const lossAt = (t) => {
      const w = archardWear({ wearCoefficient: inp.kWear * 1e-16, contactPressure: inp.pContact * 1e6, slidingPerCycle: slip.slidingPerCycle, cycles: inp.cyclesYr * t, thickness: d0 });
      const c = corrosionWastage({ time: t, coatingLife: inp.coatLife, rate: inp.corrRate / 1e3 });
      // wear flattens one side of the wire, corrosion removes a uniform radial layer: remaining section ≈ (d − 2c)·(d − 2c − h)/d² of the original
      const dc = Math.max(0, d0 - 2 * c.depth), rem = Math.max(0, dc * Math.max(0, dc - w.wearDepth) / (d0 * d0));
      return { wear: w.wearDepth, corr: c.depth, areaLoss: 1 - rem };
    };
    const eol = lossAt(L);
    const capEol = inp.Tcap * (1 - inp.armShare * eol.areaLoss);
    // functional redundancy
    const qL = 1 - Math.exp(-inp.lambda * L), funcRel = kOfN(inp.nLines, inp.nReq, qL);
    const ucMbr = inp.mbrAllow / inp.mbrActual;
    const checks = [
      check('Combined fatigue damage with DFF over the design life', comb.factoredDamage, 1, '-', `${STD_UMB} §7.4; DNV-RP-C203 (2021) §2.2`, 'Sum of the Miner damage of all sources, each times its design fatigue factor. Above 1.0 the fatigue budget is exceeded before the end of the design life.'),
      check('Annual failure probability at the design life', atL.annualPf[0], inp.pfTarget, '1/yr', 'DNV-RP-C203 (2021) App. F; Wirsching (1984) lognormal format; target per project safety class', 'Notional probability that fatigue failure of the governing element occurs in the last year of the design life against the target. Above it, inspection, monitoring or replacement is needed before that year.'),
      check('Armour wire section loss (wear + corrosion)', eol.areaLoss * 100, inp.lossAllow, '%', 'Archard (1953); API RP 17B (2014) §8 (wear and corrosion of armour) — allowable loss is a project criterion', 'Loss of wire cross-section at the end of the design life from fretting wear and corrosion. Beyond the allowable the remaining wires are overstressed and break progressively.'),
      check('End-of-life tension capacity', inp.Tmax, capEol, 'kN', `${STD_UMB} §7.3 with degraded armour`, 'Maximum service tension against the allowable tension reduced for the wire section lost by the end of life.'),
      check('Function retained with the available spare lines', 1 - funcRel, 1 - inp.relTarget, '-', 'k-out-of-n reliability with exponential line life (ISO 14224 data format)', 'Probability of having fewer working lines than required at the end of the design life, against the accepted probability. Above it more spares are needed or an intervention must be planned.'),
      check('Design-compliance summary (largest utilisation)', Math.max(inp.ucStress, inp.ucTube, inp.ucCollapse, ucMbr), 1, '-', `${STD_UMB} §7 (compliance of each design check)`, 'Largest utilisation among the stress, pressure, collapse and bend-radius checks reported by the other modules. Above 1.0 at least one design criterion is not met.'),
    ];
    const gov = worst(checks);
    const serviceLife = Math.min(lifeDff, relLife);
    const tt = linspace(0, tMax, 41);
    const risk = (uc) => (uc > 1 ? 'High' : uc > 0.7 ? 'Medium' : 'Low');
    const comp = [
      { item: 'Tension and bending stress', uc: inp.ucStress, mode: 'Yield / overstress of tubes and armour' },
      { item: 'Tube pressure containment', uc: inp.ucTube, mode: 'Leak or burst of a tube' },
      { item: 'Hydrostatic collapse', uc: inp.ucCollapse, mode: 'Collapse of a depressurised tube' },
      { item: 'Minimum bend radius', uc: ucMbr, mode: 'Overbending / kinking' },
      { item: 'Fatigue (all sources, with DFF)', uc: checks[0].uc, mode: 'Fatigue cracking of a tube or wire' },
      { item: 'Armour wear and corrosion', uc: checks[2].uc, mode: 'Armour-wire damage, loss of tensile capacity' },
      { item: 'End-of-life tension capacity', uc: checks[3].uc, mode: 'Tensile failure with degraded armour' },
      { item: 'Functional redundancy', uc: checks[4].uc, mode: 'Loss of a hydraulic / electrical function' },
    ];
    const deferYears = Math.max(0, serviceLife - L);
    return {
      kpis: [
        { key: 'damageFactored', label: 'Factored fatigue damage at the design life', value: round(comb.factoredDamage), unit: '-', limit: 1, sense: 'max', explain: 'Miner sum of all sources including their design fatigue factors at the end of the design life. It must not exceed 1.0.' },
        { key: 'lifeDff', label: 'Fatigue design life (with DFF)', value: round(Math.min(lifeDff, 1e6)), unit: 'years', limit: L, sense: 'min', explain: 'Time at which the factored Miner sum reaches 1.0, after deducting the one-off installation damage.' },
        { key: 'remaining', label: 'Remaining fatigue life from the present age', value: round(Math.min(remaining, 1e6)), unit: 'years', explain: 'Fatigue design life minus the years already served: the basis of a life-extension case.' },
        { key: 'pfLife', label: 'Cumulative failure probability at the design life', value: round(atL.pf[0], 3), unit: '-', explain: 'Notional probability that the governing element has failed in fatigue by the end of the design life (lognormal format; for ranking, not an actuarial value).' },
        { key: 'pfAnnual', label: 'Annual failure probability in the last year', value: round(atL.annualPf[0], 3), unit: '1/yr', limit: inp.pfTarget, sense: 'max', explain: 'Increase of the failure probability during the final year of the design life; compared with the target annual value for the safety class.' },
        { key: 'beta', label: 'Reliability index at the design life', value: round(atL.beta[0], 3), unit: '-', explain: 'Cumulative reliability index β = −Φ⁻¹(Pf). About 3.7 corresponds to 10⁻⁴; lower values call for inspection or monitoring.' },
        { key: 'relLife', label: 'Reliability-based service life', value: round(relLife), unit: 'years', explain: 'Time at which the annual failure probability first exceeds the target: the latest date for inspection, re-assessment or replacement on a risk basis.' },
        { key: 'serviceLife', label: 'Predicted service life (governing)', value: round(Math.min(serviceLife, 1e6)), unit: 'years', limit: L, sense: 'min', explain: 'Smaller of the DFF-based fatigue life and the reliability-based life.' },
        { key: 'wearDepth', label: 'Armour wear depth at the design life', value: round(eol.wear * 1e3, 5), unit: 'mm', explain: 'Material worn from the wire at the contact with its neighbour by fretting slip (Archard: coefficient × pressure × sliding distance).' },
        { key: 'corrDepth', label: 'Corrosion wastage at the design life', value: round(eol.corr * 1e3, 3), unit: 'mm', explain: 'Radial metal loss of the wire after the galvanising / sheath protection is exhausted.' },
        { key: 'areaLoss', label: 'Armour wire section lost at the design life', value: round(eol.areaLoss * 100), unit: '%', limit: inp.lossAllow, sense: 'max', explain: 'Combined loss of wire cross-section from wear and corrosion; it raises wire stress in proportion.' },
        { key: 'capEol', label: 'Allowable tension at end of life', value: round(capEol), unit: 'kN', limit: inp.Tmax, sense: 'min', explain: 'Allowable tension reduced for the armour section lost; it must still exceed the maximum service tension.' },
        { key: 'funcRel', label: 'Probability of retaining the function at the design life', value: round(funcRel), unit: '-', limit: inp.relTarget, sense: 'min', explain: `Probability that at least ${inp.nReq} of the ${inp.nLines} lines still work at the end of the design life. Spare lines are what keeps this high.` },
        { key: 'complianceMax', label: 'Largest design utilisation', value: round(Math.max(...comp.map((c) => c.uc))), unit: '-', limit: 1, sense: 'max', explain: 'Highest utilisation across all integrity checks in the compliance summary — the single number for the design-compliance report.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Fatigue damage budget by source (factored, at the design life)', categories: comb.shares.map((s) => s.name), series: [{ name: 'Factored damage [-]', y: comb.shares.map((s) => round(s.factoredDamage, 3)) }],
          explain: 'Contribution of each damage source to the Miner sum including its design fatigue factor. The largest bar is where analysis refinement or mitigation pays most; the bars must total no more than 1.0.' },
        { type: 'line', title: 'Fatigue failure probability over time (log scale)', x: { label: 'Service time', unit: 'years' }, y: { label: 'log10 probability', unit: '-' },
          series: [{ name: 'Cumulative Pf', x: Array.from(times), y: Array.from(rel.pf, (p) => Math.log10(Math.max(p, 1e-12))) }, { name: 'Annual Pf', x: Array.from(times), y: Array.from(rel.annualPf, (p) => Math.log10(Math.max(p, 1e-12))) }, { name: 'Target annual Pf', x: [times[0], tMax], y: [Math.log10(inp.pfTarget), Math.log10(inp.pfTarget)], dash: true }, { name: 'Design life', x: [L, L], y: [-12, 0], dash: true }],
          explain: 'Notional probability of fatigue failure against time. The year in which the annual curve crosses the target line is the reliability-based service life; if that is before the design-life marker, inspection or mitigation is required.' },
        { type: 'line', title: 'Armour degradation over time', x: { label: 'Service time', unit: 'years' }, y: { label: 'Loss', unit: 'mm / %' },
          series: [{ name: 'Wear depth [mm]', x: tt, y: tt.map((t) => lossAt(t).wear * 1e3) }, { name: 'Corrosion depth [mm]', x: tt, y: tt.map((t) => lossAt(t).corr * 1e3) }, { name: 'Section loss [%]', x: tt, y: tt.map((t) => lossAt(t).areaLoss * 100) }, { name: 'Allowable section loss [%]', x: [0, tMax], y: [inp.lossAllow, inp.lossAllow], dash: true }],
          explain: 'Wear grows linearly from the start; corrosion starts when the protection is exhausted. The section-loss curve crossing the allowable line marks the end of life of the armour.' },
        { type: 'line', title: 'Functional reliability over time', x: { label: 'Service time', unit: 'years' }, y: { label: 'Probability of retaining the function', unit: '-' },
          series: [{ name: `${inp.nReq} of ${inp.nLines} lines`, x: tt, y: tt.map((t) => kOfN(inp.nLines, inp.nReq, 1 - Math.exp(-inp.lambda * t))) }, { name: 'No spare lines', x: tt, y: tt.map((t) => Math.exp(-inp.lambda * t * inp.nReq)), dash: true }, { name: 'Required', x: [0, tMax], y: [inp.relTarget, inp.relTarget], dash: true }],
          explain: 'Probability that enough lines of the service are still working. The gap between the solid and dashed curves is what the spare lines buy; where the solid curve crosses the requirement an intervention should be planned.' },
        { type: 'bar', title: 'Integrity and design-compliance summary', categories: comp.map((c) => c.item), series: [{ name: 'Utilisation [-]', y: comp.map((c) => round(c.uc)) }],
          explain: 'One bar per integrity check, including those passed in from the other modules. Bars above 1.0 are non-compliances; the tallest compliant bar shows where the least margin remains.' },
      ],
      tables: [
        { title: 'Integrity and design-compliance report', explain: 'One row per failure mode: the check that guards against it, its utilisation, compliance status and a qualitative risk indicator (Low ≤ 0.7 < Medium ≤ 1.0 < High). Use it as the summary sheet of the integrity file.',
          columns: [{ key: 'item', label: 'Check' }, { key: 'mode', label: 'Failure mode' }, { key: 'uc', label: 'Utilisation', unit: '-' }, { key: 'status', label: 'Status' }, { key: 'risk', label: 'Risk indicator' }],
          rows: comp.map((c) => ({ item: c.item, mode: c.mode, uc: round(c.uc), status: c.uc <= 1 ? 'Complies' : 'Does not comply', risk: risk(c.uc) })) },
        { title: 'Damage sources', explain: 'Each source of fatigue damage with its unfactored damage at the design life, design fatigue factor, factored damage and share of the budget.',
          columns: [{ key: 'name', label: 'Source' }, { key: 'd', label: 'Damage at design life', unit: '-' }, { key: 'dff', label: 'DFF', unit: '-' }, { key: 'fd', label: 'Factored damage', unit: '-' }, { key: 'share', label: 'Share', unit: '%' }],
          rows: comb.shares.map((s, i) => ({ name: s.name, d: round(s.damage, 3), dff: src[i].dff, fd: round(s.factoredDamage, 3), share: round(s.share * 100, 3) })) },
      ],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        `Fatigue design life ${fmt(Math.min(lifeDff, 1e6))} years; reliability-based life ${fmt(relLife)} years at a target annual failure probability of ${inp.pfTarget}; design life ${L} years.`,
        `At ${L} years the armour wire has lost ${fmt(eol.areaLoss * 100)} % of its section (wear ${fmt(eol.wear * 1e3)} mm, corrosion ${fmt(eol.corr * 1e3)} mm) and the allowable tension is ${fmt(capEol)} kN.`,
        `Function retained with probability ${fmt(funcRel)} (${inp.nReq} of ${inp.nLines} lines, ${inp.lambda} failures per line-year).`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}. Predicted service life ${fmt(Math.min(serviceLife, 1e6))} years against ${L} years; annual failure probability in the last year ${fmt(atL.annualPf[0])}.`,
          gov.uc > 1 ? (gov.name.includes('fatigue') || gov.name.includes('probability') ? `Reduce the dominant damage source (${comb.shares.reduce((a, s) => (s.factoredDamage > a.factoredDamage ? s : a), comb.shares[0]).name}) or plan monitoring and a re-assessment before year ${fmt(Math.max(1, Math.min(lifeDff, relLife)))}; without inspection access the DFF cannot be relaxed.` : gov.name.includes('Function') ? `Add spare lines: ${inp.nLines + 1} lines give a retention probability of ${fmt(kOfN(inp.nLines + 1, inp.nReq, qL))}.` : 'Address the non-complying check in its own module and update the summary.')
            : `Record the load history (tension, hang-off angle, pressure cycles) from first use so that the remaining life can be re-assessed with measured data; first re-assessment no later than year ${fmt(Math.min(relLife, L))}.`,
          `${STD_UMB} §7.4; DNV-RP-C203 (2021) App. F; DNV-RP-G101 (2021); NORSOK U-009 (2011) life extension of subsea systems`),
        adv('economic', deferYears > 0 ? 'info' : 'warn', deferYears > 0 ? `The predicted life exceeds the design life by ${fmt(deferYears)} years.` : `The predicted life is ${fmt(L - serviceLife)} years short of the design life: one replacement of the dynamic section falls inside the field life.`,
          deferYears > 0 ? `That margin is the basis of a life-extension case: deferring a replacement campaign of ${inp.campaignDays} days (about ${fmt(inp.campaignDays * inp.dayRate / 1e6)} million USD of vessel time at the indicative rate, plus the umbilical itself) is worth far more than the monitoring needed to justify it.` : `Budget a replacement campaign of ${inp.campaignDays} days (about ${fmt(inp.campaignDays * inp.dayRate / 1e6)} million USD of vessel time at the indicative rate, plus the new section and deferred production), or spend on the mitigation above now.`,
          'Whole-life costing ISO 15686-5; value of information and risk-based inspection DNV-RP-G101 / API RP 580',
          `Each replacement avoided saves ≈ ${fmt(inp.steelT * CO2_STEEL)} t CO₂e of steel re-manufacture (screening ${CO2_STEEL} t CO₂e/t for ${fmt(inp.steelT)} t) and the campaign fuel: ${vesselCo2Text(inp.campaignDays)} A recovered umbilical can be separated into steel and thermoplastics for recycling.`),
      ],
      publish: { 'life.fatigueYears': Math.min(lifeDff, 1e6), 'life.designYears': L, 'umbilical.serviceLifeYears': Math.min(serviceLife, 1e6), 'umbilical.complianceMax': Math.max(...comp.map((c) => c.uc)) },
      quantities: { vessel_days: { umbilicalReplacement: deferYears > 0 ? 0 : inp.campaignDays } },
    };
  },
  verify: [
    // One source, rate 0.004 /yr, DFF 10, no one-off damage: factored damage at 25 yr = 10 × 0.004 × 25 = 1.000;
    // fatigue design life = 1/(10 × 0.004) = 25 yr.
    { name: 'Miner budget: DFF·ḋ·T = 1 at the design life',
      source: 'DNV-RP-C203 (2021) §2.2 (linear damage accumulation with design fatigue factor); hand calculation',
      inputs: { sources: [{ name: 'Wave', rate: 0.004, once: 0, dff: 10 }] },
      expect: [{ kpi: 'damageFactored', value: 1.0, tol: 1e-6 }, { kpi: 'lifeDff', value: 25, tol: 1e-6 }] },
    // 2-out-of-3 lines, λ = 0.01 /yr, T = 25 yr: r = e^(−0.25) = 0.77880, q = 0.22120;
    // P = r³ + 3r²q = 0.47237 + 3 × 0.60653 × 0.22120 = 0.87486.
    { name: 'k-out-of-n redundancy: P = Σ C(n,j)·r^j·q^(n−j)',
      source: 'Binomial reliability of identical independent components (Rausand & Høyland, System Reliability Theory §4.3); hand calculation',
      inputs: { nLines: 3, nReq: 2, lambda: 0.01 },
      expect: [{ kpi: 'funcRel', value: 0.87486, tol: 0.0005 }] },
    // Archard: slip range = R²cos²α·Δκ/sinα = 0.059² × 0.93301 × 0.002/0.25882 = 2.5097e-5 m; sliding per cycle = 2 × = 5.0194e-5 m;
    // h = k·p·s·N = 2e-16 × 15e6 × 5.0194e-5 × (5e6 × 25) = 0.018823 mm.
    { name: 'Archard wear depth h = k·p·s·N with helical slip per bending cycle',
      source: 'Archard (1953) J. Appl. Phys. 24; Féret & Bournazel (1987) wire slip; hand calculation',
      inputs: {},
      expect: [{ kpi: 'wearDepth', value: 0.018823, tol: 0.003 }] },
  ],
};
