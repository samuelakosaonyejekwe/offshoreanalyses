// Structural integrity and service life of a critical hull detail: corrosion wastage against the renewal
// thickness, fatigue damage accumulating faster as the plate thins, fracture-mechanics growth of a surface crack
// to the critical size on the failure assessment diagram, the inspection interval that follows, and the fatigue
// reliability index over the life.
import { corrosionWastage } from '../../engine/fatigue/degradation.js';
import { parisLawBs7910, surfaceCrackGrowth, assessSurfaceCrack, criticalCrackSize, fadOption1 } from '../../engine/fatigue/fracture.js';
import { form } from '../../engine/reliability/form.js';
import { need, run, round, rnd, fmt, adv, check, checkMin, sev, worst, linspace, thin, steelCo2Text } from './_common.js';

const CAP = 1000; // years reported when a limit is never reached
const ENV = [{ value: 'air', label: 'Air / dry internal space' }, { value: 'marine-cp850', label: 'Sea water with cathodic protection' }, { value: 'marine-free', label: 'Sea water, free corrosion' }];

export default {
  id: 's11-integrity-life',
  suite: 11,
  title: 'Structural integrity and service life',
  summary: 'How long a critical hull detail lasts: when corrosion reaches the renewal thickness, when fatigue damage is used up, how fast a crack grows to its critical size, what inspection interval that allows and how the reliability falls with age.',
  covers: [18],
  fidelity: 'analytical',
  standards: ['BS 7910:2019 clauses 7, 8 and Annexes M, P (fracture and fatigue crack growth assessment)', 'DNV-RP-C203 (2021) §2, App. F (fatigue design, fracture mechanics)', 'DNV-RP-C210 (2019) (probabilistic methods for planning of inspection for fatigue cracks)', 'DNV-RP-C101 (2021) / IACS UR Z10 (corrosion additions and renewal thickness)', 'DNV-RP-G101 (2021) (risk-based inspection)', 'ISO 19904-1 (2019) §18 (in-service inspection, monitoring and maintenance)', 'ISO 2394 (2015) (reliability of structures)'],
  models: ['Corrosion wastage after coating breakdown (linear or Guedes Soares–Garbatov)', 'Miner fatigue accumulation with stress growing as the plate thins (damage rate ∝ (t₀/t)^m)', 'Paris–Erdogan crack growth of a semi-elliptical surface crack with Newman–Raju stress intensities (BS 7910 laws)', 'BS 7910 Option 1 failure assessment diagram and critical crack size', 'Inspection interval from the detectable-to-critical growth time', 'First-order reliability method on the Miner limit state'],
  limits: 'One detail in a flat plate under membrane and bending stress: the equivalent constant-amplitude stress range, the weld-toe magnification (not applied) and the residual stress are user inputs; corrosion is uniform and thins the plate without changing the crack geometry solution. Crack growth runs on the thickness at mid-life with a zero growth threshold, because the equivalent range stands for a variable-amplitude spectrum. The reliability index covers the Miner sum only, with lognormal uncertainties supplied by the user, and is not updated by inspection results. Detailed finite-element hot-spot analysis and hull monitoring data are outside the module.',
  publishes: ['vessel.serviceLife'],
  searchTerms: 'hull structural integrity life extension corrosion fatigue crack growth inspection interval',
  inputs: [
    { key: 'steel', label: 'Steel', type: 'select', source: 'materials:steel', default: 'ah36', group: 'Detail', role: 'material', help: 'Grade of the plate at the detail; sets yield and tensile strength.' },
    { key: 't0', label: 'As-built plate thickness', unit: 'mm', type: 'number', default: 22, min: 6, max: 120, group: 'Detail', role: 'geometry', design: true, help: 'Gross thickness at the detail when new.' },
    { key: 'W', label: 'Plate width at the detail', unit: 'm', type: 'number', default: 0.85, min: 0.1, max: 10, group: 'Detail', role: 'geometry', help: 'Width available for the crack to grow along the weld (stiffener spacing or bracket length).' },
    { key: 'sigMax', label: 'Extreme tensile stress at the detail', unit: 'MPa', type: 'number', default: 230, min: 0, max: 700, group: 'Loading (boundary conditions)', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Largest primary membrane stress in the design storm (still-water plus wave), from the global or local strength modules.' },
    { key: 'Qres', label: 'Residual stress / yield', unit: '-', type: 'number', default: 0.5, min: 0, max: 1, group: 'Loading (boundary conditions)', role: 'bc', help: 'Secondary (welding residual) membrane stress as a fraction of yield: 1.0 as-welded, about 0.2–0.3 after stress relief; it relaxes under proof loading.' },
    { key: 'dsEq', label: 'Equivalent constant-amplitude stress range', unit: 'MPa', type: 'number', default: 22, min: 0, max: 300, group: 'Loading (boundary conditions)', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Stress range that gives the same crack growth as the long-term distribution at the yearly cycle count (m-th root of the mean m-th power, m = 3).' },
    { key: 'bendFrac', label: 'Bending share of the stress range', unit: '-', type: 'number', default: 0.2, min: 0, max: 1, group: 'Loading (boundary conditions)', role: 'bc', help: 'Fraction of the hot-spot stress range that is through-thickness bending; bending slows growth through the wall.' },
    { key: 'nYear', label: 'Stress cycles per year', unit: '1/yr', type: 'number', default: 5e6, min: 1e4, max: 5e7, group: 'Loading (boundary conditions)', role: 'bc', help: 'About 5 million wave cycles per year for a permanently moored unit.' },
    { key: 'L0', label: 'Calculated fatigue life, as built (design S-N curve)', unit: 'years', type: 'number', default: 150, min: 1, max: 5000, group: 'Fatigue', role: 'ic', bind: 'life.fatigueYears', help: 'S-N fatigue life of the detail at as-built thickness without the design fatigue factor, from the fatigue module.' },
    { key: 'mSN', label: 'S-N slope m', unit: '-', type: 'number', default: 3, min: 2.5, max: 5, group: 'Fatigue', role: 'material', help: 'Negative inverse slope of the governing part of the S-N curve; stress grows as the plate thins and damage with stress^m.' },
    { key: 'DFF', label: 'Design fatigue factor', unit: '-', type: 'number', default: 2, min: 1, max: 10, group: 'Acceptance criteria', role: 'criteria', help: '1 to 10 depending on consequence and access for inspection and repair (DNV-OS-C101).' },
    { key: 'designLife', label: 'Required service life', unit: 'years', type: 'number', default: 25, min: 1, max: 100, group: 'Acceptance criteria', role: 'criteria', bind: 'life.designYears', help: 'Design or extended service life of the unit on site.' },
    { key: 'age', label: 'Present age', unit: 'years', type: 'number', default: 0, min: 0, max: 80, group: 'Initial condition', role: 'ic', help: 'Years already in service; damage and wastage accumulated so far are computed with the same models.' },
    { key: 'a0', label: 'Initial crack depth', unit: 'mm', type: 'number', default: 0.5, min: 0.05, max: 30, group: 'Initial condition', role: 'ic', uncertain: { dist: 'lognormal', cov: 0.5 }, help: 'Depth of the pre-existing weld-toe flaw (0.1–1 mm for sound welds) or of a crack found by inspection.' },
    { key: 'aspect', label: 'Initial crack aspect ratio a/c', unit: '-', type: 'number', default: 0.2, min: 0.05, max: 1, group: 'Initial condition', role: 'ic', help: 'Depth over half-length: 0.2 for long shallow weld-toe cracks, 1.0 for a semicircular flaw.' },
    { key: 'corrModel', label: 'Corrosion model', type: 'select', options: [{ value: 'linear', label: 'Constant rate after coating life' }, { value: 'guedes-soares-garbatov', label: 'Guedes Soares–Garbatov (levelling off)' }], default: 'linear', group: 'Corrosion', role: 'environment', help: 'Linear is the customary conservative choice; the nonlinear model levels off at a long-term wastage.' },
    { key: 'coatLife', label: 'Coating life', unit: 'years', type: 'number', default: 10, min: 0, max: 60, group: 'Corrosion', role: 'environment', help: 'Years before the coating breaks down and wastage starts (10–15 for ballast-tank epoxy to PSPC).' },
    { key: 'corrRate', label: 'Corrosion rate per exposed side', unit: 'mm/yr', type: 'number', default: 0.1, min: 0, max: 1, group: 'Corrosion', role: 'environment', uncertain: { dist: 'lognormal', cov: 0.4 }, showIf: { key: 'corrModel', in: ['linear'] }, help: '0.05–0.1 in cargo tanks, 0.1–0.2 in ballast tanks, 0.3–0.5 in the splash zone. Project values govern.' },
    { key: 'dInf', label: 'Long-term wastage per exposed side', unit: 'mm', type: 'number', default: 2, min: 0.05, max: 15, group: 'Corrosion', role: 'environment', showIf: { key: 'corrModel', in: ['guedes-soares-garbatov'] }, help: 'Asymptotic wastage of the nonlinear model.' },
    { key: 'tauT', label: 'Transition time', unit: 'years', type: 'number', default: 12, min: 1, max: 60, group: 'Corrosion', role: 'environment', showIf: { key: 'corrModel', in: ['guedes-soares-garbatov'] }, help: 'Time constant of the nonlinear model after coating breakdown.' },
    { key: 'nSides', label: 'Corroding sides', unit: '-', type: 'int', default: 1, min: 0, max: 2, group: 'Corrosion', role: 'environment', help: '2 for a bulkhead between two ballast tanks, 1 for shell plating with cathodic protection outside, 0 for dry coated spaces.' },
    { key: 'renewFrac', label: 'Allowable thickness diminution', unit: '-', type: 'number', default: 0.2, min: 0.05, max: 0.4, group: 'Acceptance criteria', role: 'criteria', help: 'Fraction of the as-built thickness that may be lost before renewal (class rules: typically 20 % for plating, less for strength-deck and bottom).' },
    { key: 'cgEnv', label: 'Crack-growth environment', type: 'select', options: ENV, default: 'marine-cp850', group: 'Fracture mechanics', role: 'environment', help: 'Selects the BS 7910 growth law (mean + 2 standard deviations, R ≥ 0.5 as for as-welded joints).' },
    { key: 'Kmat', label: 'Fracture toughness K_mat', unit: 'MPa·√m', type: 'number', default: 120, min: 20, max: 400, group: 'Fracture mechanics', role: 'material', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Toughness at the minimum design temperature, from CTOD tests or a Charpy correlation (BS 7910 Annex J). 60–100 for older hull steel, 150+ for modern TMCP plate.' },
    { key: 'aDet', label: 'Detectable crack depth', unit: 'mm', type: 'number', default: 3, min: 0.2, max: 40, group: 'Inspection', role: 'operation', help: 'Depth found with high probability by the chosen method: about 1–2 mm for eddy current or MPI on a cleaned weld, 5 mm and more for close visual inspection.' },
    { key: 'inspFactor', label: 'Inspections required within the growth period', unit: '-', type: 'number', default: 2, min: 1, max: 5, group: 'Inspection', role: 'criteria', help: 'The interval is the detectable-to-critical growth time divided by this factor, so that a crack is seen at least this many times before it becomes critical.' },
    { key: 'inspInterval', label: 'Planned inspection interval', unit: 'years', type: 'number', default: 5, min: 0.5, max: 30, group: 'Inspection', role: 'operation', help: 'Class renewal surveys come every 5 years; in-water or continuous survey schemes differ.' },
    { key: 'covDelta', label: 'Uncertainty of the Miner sum at failure (CoV)', unit: '-', type: 'number', default: 0.3, min: 0.05, max: 1, group: 'Reliability', role: 'criteria', help: 'Lognormal with mean 1.0 (DNV-RP-C210 uses CoV 0.3).' },
    { key: 'covB', label: 'Uncertainty of the stress calculation (CoV)', unit: '-', type: 'number', default: 0.2, min: 0.02, max: 0.6, group: 'Reliability', role: 'criteria', help: 'Lognormal with mean 1.0 on the stress range: 0.1–0.15 for refined analysis with monitoring, 0.2–0.3 for simplified load and stress models.' },
    { key: 'betaTarget', label: 'Target reliability index (accumulated)', unit: '-', type: 'number', default: 2.0, min: 0.5, max: 5, group: 'Reliability', role: 'criteria', help: 'Project target for the accumulated probability of a fatigue crack at the end of life (β = 2.3 ≈ 10⁻², 3.1 ≈ 10⁻³); higher for details that cannot be inspected or whose failure is critical.' },
  ],
  compute(inp, ctx) {
    const mat = ctx.material(inp.steel); need(mat && mat.fy > 0, 'Select a steel from the materials library.');
    const fy = mat.fy; const fu = mat.fu || 1.3 * mat.fy; const E = mat.E || 206e9;
    const t0 = inp.t0 / 1e3; const m = inp.mSN; const Tdes = inp.designLife;
    need(inp.a0 < inp.t0, `The initial crack depth (${fmt(inp.a0)} mm) must be smaller than the plate thickness (${fmt(inp.t0)} mm).`);
    // --- corrosion -------------------------------------------------------------------------------------------
    const waste = (yr) => inp.nSides * run(() => corrosionWastage(inp.corrModel === 'linear' ? { time: yr, model: 'linear', coatingLife: inp.coatLife, rate: inp.corrRate / 1e3 } : { time: yr, model: 'guedes-soares-garbatov', coatingLife: inp.coatLife, dInf: inp.dInf / 1e3, transitionTime: inp.tauT }), 'Corrosion wastage').depth;
    const tMin = t0 * (1 - inp.renewFrac);
    const horizon = Math.max(2 * Tdes, Tdes + 20, 60); const nT = 400; const yrs = linspace(0, horizon, nT + 1);
    const thick = yrs.map((y) => Math.max(t0 - waste(y), 1e-4));
    const firstYear = (arr, pred) => { for (let i = 1; i < arr.length; i++) if (pred(arr[i])) { const a = arr[i - 1]; const b = arr[i]; const f = pred(a) ? 0 : (arr === thick ? (a - tMin) / (a - b) : (1 - a) / (b - a)); return yrs[i - 1] + f * (yrs[i] - yrs[i - 1]); } return null; };
    const yRenew = firstYear(thick, (v) => v <= tMin);
    // --- fatigue with thinning plate: dD/dt = (1/L0)·(t0/t)^m ---------------------------------------------------
    const dmg = [0]; for (let i = 1; i <= nT; i++) dmg.push(dmg[i - 1] + 0.5 * ((t0 / thick[i]) ** m + (t0 / thick[i - 1]) ** m) * (yrs[i] - yrs[i - 1]) / inp.L0);
    const at = (arr, y) => { const x = Math.min(Math.max(y, 0), horizon) / horizon * nT; const i = Math.min(Math.floor(x), nT - 1); return arr[i] + (arr[i + 1] - arr[i]) * (x - i); };
    const dDes = at(dmg, Tdes); const dNow = at(dmg, inp.age);
    let yFat = firstYear(dmg.map((d) => d * inp.DFF), (v) => v >= 1);
    if (yFat === null) { const rate = (dmg[nT] - dmg[nT - 1]) / (yrs[nT] - yrs[nT - 1]); yFat = Math.min(horizon + (1 / inp.DFF - dmg[nT]) / Math.max(rate, 1e-12), CAP); }
    // --- fracture mechanics ------------------------------------------------------------------------------------
    const tCg = Math.max(t0 - waste(Tdes / 2), 0.3 * t0);
    const a0 = Math.min(inp.a0 / 1e3, 0.95 * tCg); const c0 = a0 / inp.aspect;
    const eca = (a, c) => run(() => assessSurfaceCrack({ a, c, t: tCg, W: inp.W, Pm: inp.sigMax * 1e6, Qm: inp.Qres * fy, yieldStrength: fy, tensileStrength: fu, E, toughness: inp.Kmat * 1e6 }), 'Fracture assessment');
    const crit = run(() => criticalCrackSize({ aspectRatio: inp.aspect, t: tCg, W: inp.W, Pm: inp.sigMax * 1e6, Qm: inp.Qres * fy, yieldStrength: fy, tensileStrength: fu, E, toughness: inp.Kmat * 1e6 }), 'Critical crack size');
    const e0 = eca(a0, Math.min(c0, 0.45 * inp.W));
    let yCrack = CAP; let yDet = CAP; let status = 'no cyclic stress'; let hist = { y: [0], a: [a0 * 1e3], c: [c0 * 1e3] }; let path = [{ Lr: e0.Lr, Kr: e0.Kr }];
    if (inp.dsEq > 0 && e0.acceptable) {
      // zero threshold: the equivalent range stands for a variable-amplitude spectrum whose larger cycles exceed ΔK0
      const law = parisLawBs7910({ environment: inp.cgEnv, curve: 'mean+2sd', rRatio: 'high', threshold: 0 });
      const ds = inp.dsEq * 1e6;
      // growth is followed to wall break-through; the life ends earlier, when the growing crack leaves the failure assessment diagram
      const g = run(() => surfaceCrackGrowth({ a0, c0: Math.min(c0, 0.39 * inp.W), t: tCg, W: inp.W, membraneRange: ds * (1 - inp.bendFrac), bendingRange: ds * inp.bendFrac, law, maxCycles: CAP * inp.nYear }), 'Crack growth');
      const hy = g.history.cycles.map((n) => n / inp.nYear);
      status = g.status;
      if (g.status === 'arrested') { status = 'below the growth threshold'; } else {
        let iEnd = g.history.a.length - 1; path = [];
        const stepI = Math.max(1, Math.floor(g.history.a.length / 60));
        for (let i = 0; i < g.history.a.length; i += stepI) {
          const r = eca(Math.min(g.history.a[i], 0.999 * tCg), g.history.c[i]); path.push({ Lr: r.Lr, Kr: r.Kr });
          if (!r.acceptable) { iEnd = i; status = 'fracture (failure assessment diagram)'; break; }
        }
        yCrack = Math.min(hy[iEnd], CAP);
        const iDet = g.history.a.findIndex((a) => a >= inp.aDet / 1e3);
        yDet = iDet < 0 || iDet > iEnd ? yCrack : hy[iDet];
        const [ty, ta, tc] = thin([hy.slice(0, iEnd + 1), g.history.a.slice(0, iEnd + 1), g.history.c.slice(0, iEnd + 1)], 300);
        hist = { y: ty, a: ta.map((v) => v * 1e3), c: tc.map((v) => v * 1e3) };
      }
    } else if (!e0.acceptable) { yCrack = 0; yDet = 0; status = 'initial crack already unacceptable'; }
    const window = Math.max(yCrack - yDet, 0);
    const interval = yCrack >= CAP ? CAP : window / inp.inspFactor;
    // --- reliability of the Miner sum ----------------------------------------------------------------------------
    const vars = { Delta: { dist: 'lognormal', mean: 1, cov: inp.covDelta }, Bs: { dist: 'lognormal', mean: 1, cov: inp.covB } };
    const beta = (D) => (D > 1e-12 ? run(() => form((x) => Math.log(x.Delta) - m * Math.log(x.Bs) - Math.log(D), vars), 'FORM').beta : 9);
    const by = linspace(Math.max(1, horizon / 40), horizon, 40); const bs = by.map((y) => Math.min(beta(at(dmg, y)), 9));
    const bDes = Math.min(beta(dDes), 9);
    let yBeta = CAP; for (let i = 0; i < by.length; i++) if (bs[i] < inp.betaTarget) { yBeta = i === 0 ? by[0] : by[i - 1] + (bs[i - 1] - inp.betaTarget) / (bs[i - 1] - bs[i]) * (by[i] - by[i - 1]); break; }
    const pf = 0.5 * erfc(bDes / Math.SQRT2);
    const yR = yRenew ?? CAP;
    const life = Math.min(yR, yFat, yBeta);
    const tDes = at(thick, Tdes);
    const checks = [
      checkMin('Thickness at the end of the service life against the renewal thickness', tDes * 1e3, tMin * 1e3, 'mm', 'IACS UR Z10 series / DNV-RU-OU-0300 (2021) (allowable diminution); DNV-RP-C101 (2021)', 'Guards the scantling: below the renewal thickness the plate no longer meets the yield and buckling checks it was designed to, and class requires it to be cropped and renewed.'),
      check('Fatigue damage at the end of the service life, with the design fatigue factor', dDes * inp.DFF, 1, '-', 'DNV-RP-C203 (2021) §2.2 (Miner sum); DNV-OS-C101 (2021) Ch.2 Sec.5 (DFF)', 'Guards against fatigue cracking within the service life. Corrosion raises the stress in the thinning plate, so damage accumulates faster late in life.'),
      checkMin('Crack-growth window against the planned inspection interval', interval, inp.inspInterval, 'years', 'DNV-RP-C210 (2019) §5; BS 7910:2019 clause 8; ISO 19904-1 (2019) §18', 'Guards the inspection strategy: a crack must be detectable for long enough before it becomes critical that the planned surveys see it at least the required number of times.'),
      checkMin('Fatigue reliability index at the end of the service life', bDes, inp.betaTarget, '-', 'DNV-RP-C210 (2019) §4 (target reliability); ISO 2394 (2015)', 'Guards the probability of a fatigue crack accumulated over the life, given the uncertainty of loads, stress analysis and the Miner rule.'),
    ];
    if (yCrack >= CAP) checks[2].uc = 0;
    const gov = worst(checks);
    const fadL = linspace(0, (fy + fu) / (2 * fy), 60);
    const tAdd = Math.max(tMin - tDes, 0) * 1e3;
    const steelPerM2 = (mm) => mm * 7.85 / 1e3 * 1000; // t per 1000 m² of plating
    return {
      kpis: [
        { key: 'tEnd', label: 'Thickness at the end of the service life', value: round(tDes * 1e3, 5), unit: 'mm', limit: round(tMin * 1e3, 5), sense: 'min', explain: `Remaining plate thickness after ${fmt(Tdes)} years; the renewal thickness is ${fmt(100 * (1 - inp.renewFrac))} % of as-built.` },
        { key: 'yRenew', label: 'Year the renewal thickness is reached', value: round(yR, 5), unit: 'years', explain: yRenew === null ? `The renewal thickness is not reached within ${fmt(horizon)} years (reported as ${CAP}).` : 'Age at which wastage equals the allowable diminution: steel renewal is due about then.' },
        { key: 'wasteEnd', label: 'Corrosion wastage at the end of the service life', value: round((t0 - tDes) * 1e3, 5), unit: 'mm', explain: 'Total thickness lost from all corroding sides after the coating life.' },
        { key: 'dmgNow', label: 'Fatigue damage accumulated to date', value: round(dNow, 5), unit: '-', explain: 'Miner sum at the present age (without the design fatigue factor).' },
        { key: 'dmgEnd', label: 'Fatigue damage at the end of the service life', value: round(dDes, 5), unit: '-', limit: round(1 / inp.DFF, 5), sense: 'max', explain: 'Miner sum after the required life including the effect of thinning; it must stay below 1/DFF.' },
        { key: 'fatigueLife', label: 'Fatigue life with corrosion and DFF', value: round(yFat, 5), unit: 'years', limit: Tdes, sense: 'min', explain: 'Age at which damage × DFF reaches 1.0. Without corrosion it equals the as-built life divided by DFF.' },
        { key: 'aCrit', label: 'Critical crack depth', value: round(crit.a * 1e3, 5), unit: 'mm', explain: crit.limitedBy === 'wall' ? 'The crack stays acceptable up to full thickness at the extreme stress: failure is by leakage, not fracture (leak-before-break in the plate).' : 'Depth at which the extreme stress drives the assessment point onto the failure assessment line: unstable fracture or plastic collapse of the ligament.' },
        { key: 'crackYears', label: 'Crack growth life from the initial flaw', value: round(yCrack, 5), unit: 'years', explain: `Years for the initial flaw to grow to the critical size or through the wall (${status}); ${CAP} means no growth.` },
        { key: 'detYears', label: 'Years until the crack is detectable', value: round(yDet, 5), unit: 'years', explain: 'Time for the flaw to reach the detectable depth of the inspection method.' },
        { key: 'inspMax', label: 'Longest admissible inspection interval', value: round(interval, 5), unit: 'years', limit: inp.inspInterval, sense: 'min', explain: 'Detectable-to-critical growth time divided by the required number of inspections; the planned interval must not be longer.' },
        { key: 'KrNow', label: 'Fracture ratio Kr of the initial flaw', value: round(e0.Kr, 5), unit: '-', explain: 'Applied stress intensity (primary plus residual) over toughness at the extreme stress, for the initial flaw.' },
        { key: 'LrNow', label: 'Load ratio Lr of the initial flaw', value: round(e0.Lr, 5), unit: '-', explain: 'Reference stress over yield: proximity to plastic collapse of the cracked section.' },
        { key: 'betaEnd', label: 'Reliability index at the end of the service life', value: round(bDes, 5), unit: '-', limit: inp.betaTarget, sense: 'min', explain: 'First-order reliability index of the Miner limit state accumulated over the life; higher is safer.' },
        { key: 'pfEnd', label: 'Accumulated probability of fatigue failure', value: round(pf, 6), unit: '-', explain: 'Φ(−β): probability that a fatigue crack has developed by the end of the service life, given the stated uncertainties.' },
        { key: 'serviceLife', label: 'Service life of the detail', value: round(life, 5), unit: 'years', limit: Tdes, sense: 'min', explain: 'Earliest of: renewal thickness reached, fatigue damage with DFF used up, reliability below target. Crack growth sets the inspection interval rather than the life.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Thickness and fatigue damage over the life', x: { label: 'Age', unit: 'years' }, y: { label: 'Thickness (mm) / damage × DFF × 10', unit: '-' },
          series: (() => { const [y, th, d] = thin([yrs, thick, dmg], 200); return [{ name: 'Plate thickness (mm)', x: rnd(y), y: rnd(th.map((v) => v * 1e3)) }, { name: 'Renewal thickness (mm)', x: [0, round(horizon)], y: [round(tMin * 1e3), round(tMin * 1e3)], dash: true }, { name: 'Damage × DFF × 10', x: rnd(y), y: rnd(d.map((v) => v * inp.DFF * 10)) }, { name: 'Damage limit (× 10)', x: [0, round(horizon)], y: [10, 10], dash: true }]; })(),
          explain: 'Thickness falls after the coating life while the factored fatigue damage (scaled by 10 to share the axis) curves upwards as the stress grows. The life ends where either curve meets its dashed limit.' },
        { type: 'line', title: 'Crack growth', x: { label: 'Time from the initial flaw', unit: 'years' }, y: { label: 'Crack size', unit: 'mm' },
          series: [{ name: 'Depth a', x: rnd(hist.y), y: rnd(hist.a) }, { name: 'Half-length c', x: rnd(hist.y), y: rnd(hist.c) }, { name: 'Detectable depth', x: [0, round(hist.y[hist.y.length - 1] || 1)], y: [inp.aDet, inp.aDet], dash: true }, { name: 'Critical depth', x: [0, round(hist.y[hist.y.length - 1] || 1)], y: [round(crit.a * 1e3), round(crit.a * 1e3)], dash: true }],
          explain: 'Depth and surface half-length of the crack against time. Growth is slow for most of the life and accelerates sharply: the gap between the detectable and critical lines is the window in which inspection can find it.' },
        { type: 'line', title: 'Failure assessment diagram', x: { label: 'Load ratio Lr', unit: '-' }, y: { label: 'Fracture ratio Kr', unit: '-' },
          series: [{ name: 'BS 7910 Option 1 line', x: rnd(fadL), y: rnd(fadL.map((l) => fadOption1(l, { yieldStrength: fy, tensileStrength: fu, E }).Kr)) }, { name: 'Crack growing at extreme stress', x: rnd(path.map((p) => p.Lr)), y: rnd(path.map((p) => p.Kr)) }],
          explain: 'Assessment points of the growing crack at the extreme stress. Inside the line the crack is tolerable; the point where the path crosses it is the critical size. A path leaving through the top is fracture-controlled, through the right side collapse-controlled.' },
        { type: 'line', title: 'Fatigue reliability over the life', x: { label: 'Age', unit: 'years' }, y: { label: 'Reliability index β', unit: '-' },
          series: [{ name: 'β (accumulated)', x: rnd(by), y: rnd(bs) }, { name: 'Target', x: [round(by[0]), round(horizon)], y: [inp.betaTarget, inp.betaTarget], dash: true }],
          explain: 'Reliability index of the Miner limit state against age. Where it crosses the target, an inspection (or a life-extension reassessment with measured loads) is needed to restore confidence.' },
      ],
      tables: [
        { title: 'Life-limiting mechanisms', columns: [{ key: 'mech', label: 'Mechanism' }, { key: 'year', label: 'Limit reached at', unit: 'years' }, { key: 'margin', label: 'Margin on the required life', unit: 'years' }, { key: 'note', label: 'Basis' }],
          rows: [
            { mech: 'Corrosion to renewal thickness', year: round(yR, 4), margin: round(yR - Tdes, 4), note: `${fmt(100 * inp.renewFrac)} % diminution of ${fmt(inp.t0)} mm` },
            { mech: 'Fatigue (Miner × DFF = 1)', year: round(yFat, 4), margin: round(yFat - Tdes, 4), note: `as-built life ${fmt(inp.L0)} yr, DFF ${fmt(inp.DFF)}, m = ${fmt(m)}` },
            { mech: 'Reliability below target', year: round(yBeta, 4), margin: round(yBeta - Tdes, 4), note: `β target ${fmt(inp.betaTarget)}` },
            { mech: 'Crack growth from the initial flaw', year: round(yCrack, 4), margin: round(yCrack - Tdes, 4), note: status },
          ],
          explain: `One row per degradation mechanism with the age at which its limit is reached (${CAP} = not reached). The smallest year of the first three is the service life; the crack-growth row drives the inspection interval.` },
      ],
      notes: [
        `Corrosion: ${fmt((t0 - tDes) * 1e3)} mm lost after ${fmt(Tdes)} years (coating life ${fmt(inp.coatLife)} yr, ${inp.nSides} side${inp.nSides === 1 ? '' : 's'}); renewal thickness ${fmt(tMin * 1e3)} mm ${yRenew === null ? 'not reached' : `reached at ${fmt(yRenew)} years`}.`,
        `Fatigue: damage ${fmt(dDes)} at ${fmt(Tdes)} years (${fmt(100 * (dDes / (Tdes / inp.L0) - 1))} % more than without thinning); with DFF ${fmt(inp.DFF)} the life is ${fmt(yFat)} years. Reliability index ${fmt(bDes)} (pf ${fmt(pf)}).`,
        `Crack growth on ${fmt(tCg * 1e3)} mm (mid-life thickness): ${status}; ${yCrack >= CAP ? 'no growth' : `${fmt(yCrack)} years from ${fmt(inp.a0)} mm to the end state, detectable after ${fmt(yDet)} years`}; critical depth ${fmt(crit.a * 1e3)} mm. Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at utilisation ${fmt(gov.uc)}; service life of the detail ${fmt(life)} years against ${fmt(Tdes)} required, admissible inspection interval ${fmt(interval)} years.`,
          checks[0].uc > 1 ? `Add ${fmt(tAdd)} mm of corrosion addition (as-built ${fmt(inp.t0 + tAdd / (1 - inp.renewFrac))} mm), or maintain the coating so that its life reaches ${fmt(Math.min(Tdes, inp.coatLife + (yR < CAP ? Tdes - yR : 0)))} years; plan steel renewal at about year ${fmt(yR)} otherwise.`
            : checks[1].uc > 1 ? `Reduce the stress range by ${fmt(100 * (1 - (1 / (dDes * inp.DFF)) ** (1 / m)))} % (thicker insert, softer bracket toe, weld-toe grinding, which is worth a factor of about 2 on life when documented), or accept a lower DFF with a matching inspection programme.`
              : checks[2].uc > 1 ? `Shorten the inspection interval to ${fmt(interval)} years or use a method that detects ${fmt(inp.aDet / 2)} mm cracks; a tougher plate or lower residual stress enlarges the critical crack and the window.`
                : checks[3].uc > 1 ? `Schedule a detailed inspection of this detail by year ${fmt(yBeta)}; a clean result restores the reliability index, as does reducing the stress uncertainty with strain monitoring.`
                  : `The detail meets the required life with ${fmt(life - Tdes)} years of margin. Confirm the coating life at each survey: the fatigue and renewal margins both assume wastage starts at year ${fmt(inp.coatLife)}.`,
          'DNV-RP-C203 (2021); BS 7910:2019; DNV-RP-C210 (2019); ISO 19904-1 (2019) §18; IACS UR Z10 series.'),
        adv('economic', 'info', `${yCrack >= CAP ? 'The initial flaw does not grow, so crack growth sets no inspection interval' : `The crack-growth window supports an inspection interval of up to ${fmt(interval)} years`}; each additional millimetre of corrosion addition is about ${fmt(steelPerM2(1))} t per 1000 m² of plating.`,
          'Put the inspection effort where the admissible interval is shortest and extend it elsewhere (risk-based inspection); compare added steel now with steel renewal on station, which costs many times more per tonne and interrupts production.',
          'Risk-based inspection per DNV-RP-G101 (2021) and API RP 580 (2016); whole-life costing per ISO 15686-5; life-extension assessment per NORSOK N-006 (2015).',
          `${steelCo2Text(steelPerM2(1))} per 1000 m² and millimetre. Keeping the coating intact avoids both the added steel and later renewal; extending the life of an existing hull by inspection avoids the embodied carbon of a replacement unit.`),
      ],
      publish: { 'vessel.serviceLife': life },
    };
  },
  verify: [
    { // Linear corrosion: 0.1 mm/yr on one side after 10 years; 20 mm plate with 20 % allowable diminution → 4 mm → 10 + 4/0.1 = 50 yr.
      // After 25 years: 1.5 mm lost → 18.5 mm. Fatigue with thinning, m = 3, L0 = 200 yr:
      // D(25) = (1/200)·[10 + ∫₀¹⁵ (20/(20 − 0.1τ))³ dτ] = (1/200)·[10 + (8000/0.2)·(1/18.5² − 1/20²)] = (10 + 16.874)/200 = 0.13437.
      name: 'Linear wastage to the renewal thickness and Miner damage of a thinning plate',
      source: 'Hand calculation: d = r·(t − Tc); dD/dt = (1/L0)·(t0/t)^m integrated in closed form',
      inputs: { t0: 20, corrModel: 'linear', coatLife: 10, corrRate: 0.1, nSides: 1, renewFrac: 0.2, L0: 200, mSN: 3, designLife: 25 },
      expect: [{ kpi: 'yRenew', value: 50, tol: 2e-3 }, { kpi: 'tEnd', value: 18.5, tol: 1e-4 }, { kpi: 'dmgEnd', value: 0.13437, tol: 2e-3 }] },
    { // No corrosion: life = L0/DFF = 200/2 = 100 yr; damage at 25 yr = 0.125. Lognormal limit state Δ = B^m·D:
      // s_Δ² = ln(1 + 0.3²) = 0.086178, μ_Δ = −0.043089; s_B² = ln(1 + 0.2²) = 0.039221, μ_B = −0.019610;
      // β = (μ_Δ − 3·μ_B − ln 0.125)/√(s_Δ² + 9·s_B²) = (−0.043089 + 0.058831 + 2.079442)/0.662697 = 3.1616.
      name: 'Fatigue life without corrosion and the closed-form lognormal reliability index',
      source: 'DNV-RP-C203 (2021) §2.2; lognormal R–S format (Madsen, Krenk & Lind, Methods of Structural Safety, 1986, ch. 4) — hand calculation',
      inputs: { nSides: 0, L0: 200, DFF: 2, mSN: 3, designLife: 25, covDelta: 0.3, covB: 0.2 },
      expect: [{ kpi: 'fatigueLife', value: 100, tol: 2e-3 }, { kpi: 'dmgEnd', value: 0.125, tol: 1e-6 }, { kpi: 'betaEnd', value: 3.1616, tol: 2e-3 }] },
    { name: 'No cyclic stress: the crack does not grow', source: 'Limiting case Δσ = 0 (Paris law)', inputs: { dsEq: 0 },
      expect: [{ kpi: 'crackYears', value: 1000, tol: 1e-9 }] },
  ],
};

// complementary error function (Abramowitz & Stegun 7.1.26, |error| < 1.5e-7), for Φ(−β)
function erfc(x) {
  const z = Math.abs(x); const t = 1 / (1 + 0.3275911 * z);
  const y = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))) * Math.exp(-z * z);
  return x >= 0 ? y : 2 - y;
}
