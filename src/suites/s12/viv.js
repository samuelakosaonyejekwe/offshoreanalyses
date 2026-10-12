// Vortex-induced vibration and hydrodynamic loading: Morison drag on the static shape in a sheared current, modal
// cross-flow VIV screening of the suspended length, VIV drag amplification and VIV fatigue of the governing tube.
import { tensionedBeamModes, vivModalScreening, vivFatigue, dragAmplification, sheddingFrequency } from '../../engine/lines/viv.js';
import { snCurve, snConstants, cyclesFromConstants, SN_CURVE_NAMES } from '../../engine/fatigue/sn.js';
import { lineInputs, configInputs, currentInputs, solveConfig, rodCurrent, currentFn, adv, check, worst, sev, fmt, round, run, need, thin, LIMIT_2D, STD_UMB, CO2_POLYMER } from './_common.js';

export default {
  id: 's12-viv-hydrodynamic', suite: 12, title: 'Vortex-induced vibration and hydrodynamic loading',
  summary: 'Current drag on the umbilical, the modes that vortex shedding can excite, response amplitude, drag amplification and VIV fatigue life.',
  covers: [10], fidelity: 'screening',
  standards: ['DNV-RP-C205 (2021) §6.2 and §9 (Morison load, VIV)', 'DNV-RP-F204 (2021) §4 (riser VIV fatigue)', `${STD_UMB} §7.3–7.4`, 'DNV-RP-C203 (2021) §2 (S–N curves)'],
  models: ["Morison's equation: normal and tangential drag in a sheared current", 'Strouhal relation and reduced-velocity lock-in window', 'Modal power-in VIV screening (Vandiver 1993)', 'Drag amplification under VIV (Vandiver 1983)', 'S–N fatigue of the tube from modal curvature'],
  limits: `${LIMIT_2D} VIV modes are those of an equivalent pinned tensioned beam over the suspended length at its mean effective tension — a screening substitute for the cross-flow modes of the real shape. Single-mode lock-in per mode, cross-flow only, no in-line VIV, no wake interference; a dedicated VIV code is needed if the screening life is short.`,
  publishes: ['umbilical.vivLifeYears', 'umbilical.cdAmplified', 'umbilical.currentDrag'],
  searchTerms: 'umbilical vortex induced vibration fatigue drag amplification',
  inputs: [
    ...lineInputs(), ...configInputs(), ...currentInputs({ Us: 0.8, Ub: 0.15 }),
    { key: 'St', label: 'Strouhal number', unit: '-', type: 'number', default: 0.2, min: 0.1, max: 0.3, group: 'VIV parameters', role: 'load' },
    { key: 'zeta', label: 'Structural damping ratio', unit: '-', type: 'number', default: 0.02, min: 0.001, max: 0.2, group: 'VIV parameters', role: 'material', uncertain: { dist: 'lognormal', cov: 0.4 },
      help: 'Umbilicals damp more than steel pipe because of interlayer friction: 1–5 % is typical; calibrate from a bending hysteresis test.' },
    { key: 'Cl', label: 'Lift coefficient (energy balance)', unit: '-', type: 'number', default: 0.7, min: 0.1, max: 1.2, group: 'VIV parameters', role: 'load' },
    { key: 'strakeFrac', label: 'Suppression coverage (fraction of suspended length from the top)', unit: '-', type: 'number', default: 0, min: 0, max: 1, group: 'VIV parameters', role: 'geometry', design: true },
    { key: 'strakeEff', label: 'Suppression efficiency', unit: '-', type: 'number', default: 0.8, min: 0, max: 1, group: 'VIV parameters', role: 'geometry' },
    { key: 'pCurrent', label: 'Fraction of time this current occurs', unit: '-', type: 'number', default: 0.05, min: 0.001, max: 1, group: 'VIV parameters', role: 'load', help: 'Probability of the current case in the long-term distribution; the damage rate is weighted by it.' },
    { key: 'tubeOD', label: 'Outside diameter of the governing steel tube', unit: 'mm', type: 'number', default: 31.8, min: 3, max: 120, group: 'Fatigue resistance', role: 'geometry', help: 'Bending stress in a tube bent about its own axis is E·(d/2)·κ: the largest tube governs.' },
    { key: 'Etube', label: 'Tube Young’s modulus', unit: 'GPa', type: 'number', default: 200, min: 50, max: 220, group: 'Fatigue resistance', role: 'material' },
    { key: 'sn', label: 'S–N curve (tube)', type: 'select', default: 'C1', options: SN_CURVE_NAMES.map((n) => ({ value: n, label: `DNV-RP-C203 ${n}` })), group: 'Fatigue resistance', role: 'material' },
    { key: 'snEnv', label: 'Environment', type: 'select', default: 'seawater-cp', options: [{ value: 'air', label: 'Air' }, { value: 'seawater-cp', label: 'Sea water with cathodic protection' }, { value: 'free-corrosion', label: 'Free corrosion' }], group: 'Fatigue resistance', role: 'material' },
    { key: 'scf', label: 'Stress concentration factor', unit: '-', type: 'number', default: 1.15, min: 1, max: 4, group: 'Fatigue resistance', role: 'material' },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 10, min: 1, max: 20, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'life', label: 'Design service life', unit: 'years', type: 'number', default: 25, min: 1, max: 60, group: 'Acceptance criteria', role: 'criteria', bind: 'life.designYears' },
    { key: 'Tallow', label: 'Allowable tension', unit: 'kN', type: 'number', default: 350, min: 5, max: 5000, group: 'Acceptance criteria', role: 'criteria', bind: 'umbilical.allowTension' },
    { key: 'nElem', label: 'Rod elements', type: 'int', default: 120, min: 30, max: 1200, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Rod elements', levels: [60, 120, 240], qoi: 'dragForce' },
  compute(inp) {
    need(inp.Us !== 0 || inp.Ub !== 0, 'Both current speeds are zero: there is no hydrodynamic loading or vortex shedding to assess.');
    const Uof = currentFn(inp, inp.depth);
    // still-water shape first; the current case and the amplified-drag case are warm-started from it
    const still = solveConfig(inp, { nElem: inp.nElem }).res;
    const { res, line, buoySec, s0, s1 } = solveConfig(inp, { nElem: inp.nElem, current: rodCurrent(inp, inp.rhoW), initial: still });
    const n = res.s.length - 1, rho = line.rhoW;
    const inB = (s) => buoySec && s >= s0 && s <= s1;
    // Morison drag per unit length on the deflected shape (DNV-RP-C205 §6.2): normal ½ρ·Cdn·D·|un|·un, tangential on πD
    const fx = new Float64Array(n + 1), Un = new Float64Array(n + 1), Dn = new Float64Array(n + 1), mn = new Float64Array(n + 1);
    for (let i = 0; i <= n; i++) {
      const sec = inB(res.s[i]) ? buoySec : line.section, U = res.z[i] >= -1e-6 ? Uof(res.z[i]) : 0, th = res.angle[i];
      const sn = Math.sin(th), cs = Math.cos(th);
      Un[i] = Math.abs(U * sn); Dn[i] = sec.D; mn[i] = sec.m + sec.ma;
      fx[i] = Math.sign(U) * (0.5 * rho * sec.Cdn * sec.D * U * U * Math.abs(sn) * sn * sn + 0.5 * rho * sec.Cdt * Math.PI * sec.D * U * U * Math.abs(cs) * cs * cs);
    }
    let Fdrag = 0;
    for (let i = 1; i <= n; i++) Fdrag += 0.5 * (fx[i] + fx[i - 1]) * (res.s[i] - res.s[i - 1]);
    // equivalent tensioned beam over the suspended length for the VIV modes
    const kEnd = res.tdp ? res.tdp.node : n, Ls = res.s[kEnd];
    let Tm = 0, mm = 0; for (let i = 1; i <= kEnd; i++) { const ds = res.s[i] - res.s[i - 1]; Tm += 0.5 * (res.Te[i] + res.Te[i - 1]) * ds; mm += 0.5 * (mn[i] + mn[i - 1]) * ds; }
    Tm /= Ls; mm /= Ls;
    need(Tm > 0, 'Mean effective tension of the suspended length is not positive: the VIV screening needs a tensioned line.');
    let UnMax = 0, fsMax = 0;
    for (let i = 0; i <= kEnd; i++) { UnMax = Math.max(UnMax, Un[i]); fsMax = Math.max(fsMax, sheddingFrequency({ U: Un[i], D: Dn[i], St: inp.St })); }
    // mode number whose string frequency reaches the highest shedding frequency: n = 2·L·f·√(m/T)
    const nNeed = Math.ceil(1.25 * 2 * Ls * fsMax * Math.sqrt(mm / Tm)) + 2, nModes = Math.max(6, Math.min(100, nNeed));
    const nPts = Math.min(1601, Math.max(401, 8 * nModes + 1));
    const tb = run(() => tensionedBeamModes({ L: Ls, EI: line.EI, m: mm, T: Tm, nModes, nPoints: nPts }), 'Tensioned-beam modes');
    const sv = tb.s, at = (arr) => { const out = new Float64Array(sv.length); let k = 1; for (let j = 0; j < sv.length; j++) { while (k < kEnd && res.s[k] < sv[j]) k++; const f = (sv[j] - res.s[k - 1]) / Math.max(res.s[k] - res.s[k - 1], 1e-12); out[j] = arr[k - 1] + Math.min(1, Math.max(0, f)) * (arr[k] - arr[k - 1]); } return out; };
    const Uv = at(Un), Dv = at(Dn), mv = at(mn);
    const suppression = inp.strakeFrac > 0 ? { zones: [[0, inp.strakeFrac * Ls]], efficiency: inp.strakeEff } : undefined;
    const scr = run(() => vivModalScreening({ s: sv, modes: tb.modes, U: Uv, D: Dv, mass: mv, zeta: inp.zeta, rho, St: inp.St, Cl: inp.Cl, suppression }), 'VIV screening');
    const sn = snConstants(run(() => snCurve(inp.sn, inp.snEnv), 'S–N curve'), { scf: 1 });
    const fat = run(() => vivFatigue({ s: sv, modes: tb.modes, screening: scr, E: inp.Etube * 1e9, stressDiameter: inp.tubeOD / 1e3, scf: inp.scf, snCyclesToFailure: (S) => cyclesFromConstants(sn, S) }), 'VIV fatigue');
    const dom = scr.dominant >= 0 ? scr.modes[scr.dominant] : null;
    const aD = dom ? dom.aOverD : 0, amp = dragAmplification({ Cd0: inp.Cdn, aRmsOverD: aD / Math.SQRT2 });
    const rate = fat.maxDamagePerYear * inp.pCurrent, life = rate > 0 ? 1 / rate : 1e6;
    const sigMax = Math.max(...fat.stressAmplitude);
    // static response with the VIV-amplified drag coefficient
    const resA = dom && amp.factor > 1.001 ? solveConfig({ ...inp, Cdn: amp.Cd }, { nElem: inp.nElem, current: rodCurrent(inp, inp.rhoW), initial: res }).res : res;
    const xShift = (r) => { let d = 0; for (let i = 0; i <= n; i++) d = Math.max(d, Math.abs(r.x[i] - still.x[i])); return d; };
    const checks = [
      check('VIV fatigue damage over the design life', inp.dff * rate * inp.life, 1, '-', 'DNV-RP-F204 (2021) §4; DNV-RP-C203 (2021) §2.2 (Miner sum × DFF ≤ 1)', 'Miner damage from vortex-induced bending cycles of the governing tube, weighted by the occurrence of the current and multiplied by the design fatigue factor. Above 1.0 VIV alone consumes the fatigue budget and suppression is required.'),
      check('Top tension with VIV-amplified drag', resA.topTension / 1e3, inp.Tallow, 'kN', `${STD_UMB} §7.3`, 'Hang-off tension in the design current with the drag coefficient increased for VIV, against the allowable tension.'),
    ];
    const gov = worst(checks);
    const exc = scr.excited.slice(0, 12).map((i) => scr.modes[i]);
    const [sp, fxp, unp] = thin([res.s, fx, Un], 300);
    const [svp, sgp, dmp] = thin([sv, Array.from(fat.stressAmplitude, (v) => v / 1e6), Array.from(fat.damagePerYear, (v) => v * inp.pCurrent)], 400);
    const strakeLen = inp.strakeFrac * Ls, strakeT = strakeLen * 6 / 1e3;   // ≈ 6 kg/m of polyurethane strake shells (indicative)
    return {
      kpis: [
        { key: 'dragForce', label: 'Total current drag on the line', value: round(Fdrag / 1e3, 6), unit: 'kN', explain: 'Horizontal Morison drag integrated over the deflected shape with the basic drag coefficient. It is shared between the hang-off and the termination as horizontal reaction.' },
        { key: 'topTension', label: 'Top tension in the current', value: round(res.topTension / 1e3), unit: 'kN', explain: 'Hang-off tension with current drag; compare with the still-water value to see how much the current loads the line.' },
        { key: 'dTcurrent', label: 'Top-tension change due to current', value: round((res.topTension - still.topTension) / 1e3, 3), unit: 'kN', explain: 'Difference between top tension with and without current. Negative values mean the current pushes the line towards the termination and slackens it.' },
        { key: 'deflection', label: 'Largest in-plane deflection by current', value: round(xShift(res), 3), unit: 'm', explain: 'Largest horizontal movement of any point of the line from its still-water position; it consumes clearance to neighbouring lines.' },
        { key: 'fsMax', label: 'Highest vortex-shedding frequency', value: round(fsMax), unit: 'Hz', explain: 'St·U/D at the fastest normal flow. Modes with natural frequencies up to about this value can lock in.' },
        { key: 'nExcited', label: 'Modes with a lock-in zone', value: scr.excited.length, unit: '-', explain: 'Number of modes whose reduced velocity falls in the lock-in window somewhere along the line. Many excited modes in a sheared current mean each has only a short power-in length.' },
        { key: 'domMode', label: 'Dominant mode number', value: dom ? dom.index + 1 : 0, unit: '-', explain: 'Mode receiving the most power from the flow. High mode numbers mean short wavelengths and high curvature per unit amplitude.' },
        { key: 'domFreq', label: 'Dominant mode frequency', value: round(dom ? dom.fHz : 0), unit: 'Hz', explain: 'Vibration frequency of the dominant mode — the cycle rate of the VIV stress.' },
        { key: 'aOverD', label: 'Response amplitude A/D (dominant mode)', value: round(aD, 3), unit: '-', explain: 'Cross-flow amplitude at the antinode in diameters. Above roughly 0.1–0.2 VIV matters for fatigue and drag; the upper bound for a flexible cylinder is about 1.' },
        { key: 'cdAmp', label: 'Drag amplification factor', value: round(amp.factor, 3), unit: '-', explain: 'Increase of the mean drag coefficient caused by cross-flow vibration (Vandiver 1983). Apply it to the drag coefficient of the global analyses.' },
        { key: 'sigmaViv', label: 'Maximum VIV stress amplitude (incl. SCF)', value: round(sigMax / 1e6), unit: 'MPa', explain: 'Bending stress amplitude of the governing tube at the worst point for the dominant mode. Twice this value is the stress range cycling at the dominant frequency.' },
        { key: 'vivLife', label: 'VIV fatigue life (no DFF)', value: round(Math.min(life, 1e6)), unit: 'years', limit: inp.dff * inp.life, sense: 'min', explain: 'Years to a Miner sum of 1.0 from VIV alone, weighted by the occurrence of this current. It must exceed the service life times the design fatigue factor.' },
        { key: 'topTensionAmp', label: 'Top tension with amplified drag', value: round(resA.topTension / 1e3), unit: 'kN', limit: inp.Tallow, sense: 'max', explain: 'Hang-off tension re-computed with the VIV-amplified drag coefficient.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Normal current velocity and drag load along the line', x: { label: 'Arc length from hang-off', unit: 'm' }, y: { label: 'Velocity [m/s] / drag [N/m ÷ 100]', unit: '-' },
          series: [{ name: 'Normal velocity [m/s]', x: sp, y: unp }, { name: 'Horizontal drag [N/m ÷ 100]', x: sp, y: fxp.map((v) => v / 100) }],
          explain: 'Flow component normal to the axis and the resulting Morison drag per metre. The buoyant section stands out because of its larger diameter; steeply inclined parts see less normal flow.' },
        { type: 'bar', title: 'Excited modes: response amplitude', categories: exc.map((m) => `mode ${m.index + 1} (${fmt(m.fHz)} Hz)`), series: [{ name: 'A/D [-]', y: exc.map((m) => round(m.aOverD, 3)) }, { name: 'Power-in fraction of length [-]', y: exc.map((m) => round(m.powerInFraction, 3)) }],
          explain: 'The twelve modes receiving most power, with their antinode amplitude and the share of the length over which the flow feeds them. The first bar is the dominant mode used for fatigue.' },
        { type: 'line', title: 'VIV stress amplitude along the suspended length', x: { label: 'Arc length from hang-off', unit: 'm' }, y: { label: 'Stress amplitude', unit: 'MPa' }, series: [{ name: 'Governing tube, dominant mode', x: svp, y: sgp }],
          explain: 'Bending stress amplitude of the governing tube from the dominant mode. The pinned-end model puts zero stress at the ends; real hang-off and touchdown boundary conditions concentrate curvature there.' },
        { type: 'line', title: 'VIV damage rate along the suspended length', x: { label: 'Arc length from hang-off', unit: 'm' }, y: { label: 'Damage per year', unit: '1/yr' }, series: [{ name: 'Weighted by current occurrence', x: svp, y: dmp }, { name: `Allowable 1/(DFF·life) = ${fmt(1 / (inp.dff * inp.life))}`, x: [0, Ls], y: [1 / (inp.dff * inp.life), 1 / (inp.dff * inp.life)], dash: true }],
          explain: 'Yearly Miner damage from VIV. Any part above the dashed line fails the fatigue criterion on VIV alone and locates where strakes or fairings are needed.' },
      ],
      tables: [{ title: 'Excited modes', explain: 'One row per excited mode, ranked by input power: frequency, reduced-velocity range along the line, power-in length, mass-damping parameter, amplitude and drag factor.',
        columns: [{ key: 'mode', label: 'Mode' }, { key: 'f', label: 'Frequency', unit: 'Hz' }, { key: 'vr', label: 'Vr range', unit: '-' }, { key: 'Lin', label: 'Power-in length', unit: 'm' }, { key: 'Ks', label: 'Ks', unit: '-' }, { key: 'aD', label: 'A/D', unit: '-' }, { key: 'cd', label: 'Drag factor', unit: '-' }],
        rows: exc.map((m) => ({ mode: m.index + 1, f: round(m.fHz), vr: `${fmt(m.vrMin)}–${fmt(m.vrMax)}`, Lin: round(m.powerInLength), Ks: round(m.Ks, 3), aD: round(m.aOverD, 3), cd: round(m.dragFactor, 3) })) }],
      notes: [
        `Current drag ${fmt(Fdrag / 1e3)} kN moves the line up to ${fmt(xShift(res))} m and changes the top tension by ${fmt((res.topTension - still.topTension) / 1e3)} kN.`,
        dom ? `Dominant VIV mode ${dom.index + 1} at ${fmt(dom.fHz)} Hz, A/D = ${fmt(aD)}; ${scr.excited.length} of ${nModes} modes have a lock-in zone.` : 'No mode has a lock-in zone in this current: VIV is not predicted.',
        `Modes of an equivalent pinned beam of ${fmt(Ls, 4)} m at mean tension ${fmt(Tm / 1e3)} kN${nNeed > 100 ? '; the shedding band extends beyond the 100 modes retained, so the highest-frequency response is not covered' : ''}.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), dom ? `VIV screening gives A/D = ${fmt(aD)} in mode ${dom.index + 1} (${fmt(dom.fHz)} Hz), a tube stress range of ${fmt(2 * sigMax / 1e6)} MPa and a VIV life of ${fmt(Math.min(life, 1e6))} years; required ${fmt(inp.dff * inp.life)} years.` : 'No lock-in is predicted for this current profile.',
          checks[0].uc > 1 ? `Fit strakes or fairings over the high-current upper part of the line (start with ${fmt(Math.min(1, Math.max(0.3, inp.strakeFrac + 0.3)) * 100)} % of the suspended length) and re-run, and use the amplified drag coefficient ${fmt(amp.Cd)} in the global analyses.` : `Use Cd = ${fmt(amp.Cd)} (factor ${fmt(amp.factor)}) in the global and interference analyses for this current, and add the VIV damage to the wave fatigue damage before judging the total life.`,
          'DNV-RP-F204 (2021) §4; DNV-RP-C205 (2021) §9'),
        adv('economic', 'info', inp.strakeFrac > 0 ? `Suppression over ${fmt(strakeLen)} m (${fmt(inp.strakeFrac * 100)} % of the suspended length) is included.` : 'No VIV suppression is fitted in this run.',
          checks[0].uc > 1 ? 'Strakes are fitted on the lay vessel as the umbilical is paid out and slow the lay of the dynamic section (typically by a few hours per 100 m); weigh that against the replacement of a section that fails early in fatigue.' : 'With the VIV criterion met without (more) suppression, leave strakes out: they raise drag by 30–60 %, add offshore fitting time and collect marine growth.',
          'Whole-life costing ISO 15686-5; installation time on vessel-day basis (AACE RP 18R-97)',
          `Strake shells are polyurethane: about ${fmt(strakeT)} t for the covered length ≈ ${fmt(strakeT * CO2_POLYMER)} t CO₂e (screening ${CO2_POLYMER} t/t). Unneeded strakes also raise drag and therefore top tension and hang-off steel.`),
      ],
      publish: { 'umbilical.vivLifeYears': Math.min(life, 1e6), 'umbilical.cdAmplified': amp.Cd, 'umbilical.currentDrag': Fdrag },
      quantities: { polymer_t: strakeT },
    };
  },
  verify: [
    // Near-vertical taut line (L = 341 m between points 340 m apart vertically and 10 m horizontally) in a uniform 1 m/s
    // current: F = ½·ρ·Cd·D·U²·L = 0.5 × 1025 × 1.2 × 0.145 × 1 × 341 = 30.41 kN (reduced by ~1–2 % by the bow of the line).
    // Strouhal: fs = St·U/D = 0.2 × 1/0.145 = 1.379 Hz.
    { name: 'Uniform current on a near-vertical line: F = ½ρ·Cd·D·U²·L and fs = St·U/D',
      source: 'Morison equation (DNV-RP-C205 §6.2) and Strouhal relation; hand calculation',
      inputs: { config: 'catenary', length: 341, span: 10, Us: 1, Ub: 1, Cdt: 0, kn: 200 },
      expect: [{ kpi: 'dragForce', value: 30.41, tol: 0.03 }, { kpi: 'fsMax', value: 1.379, tol: 0.01 }] },
  ],
};
