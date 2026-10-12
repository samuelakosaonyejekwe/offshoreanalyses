// Vortex-induced vibration of a riser in sheared current: natural modes, modal power-in screening, response
// amplitude, drag amplification, strake coverage, VIV fatigue over a current scatter, and a wake-oscillator check.
import { tensionedBeamModes, vivModalScreening, vivFatigue, wakeOscillator, suppressionCoverage } from '../../engine/lines/viv.js';
import { rodModes } from '../../engine/lines/rod.js';
import { cyclesFromConstants } from '../../engine/fatigue/sn.js';
import {
  need, run, fmt, round, thin, pipeInputs, waterInputs, hydroInputs, currentInputs, currentFn, pipeFrom, snInputs, snFrom, adv, check, kpi, plot, meshInput, finish,
  handWeight, steelTonnes, steelCo2Text, CONFIRM,
} from './_common.js';
import { buildScr } from './_scr.js';

const vert = { key: 'config', in: ['vertical'] }, scr = { key: 'config', in: ['scr'] };
// hand values for the verification case (273.1 × 20.6 mm, 3 mm coating, contents 700 kg/m³, Ca = 1)
const Q = Math.PI / 4, M_HAND = Q * (0.2731 ** 2 - 0.2319 ** 2) * 7850 + Q * (0.2791 ** 2 - 0.2731 ** 2) * 900 + Q * 0.2319 ** 2 * 700 + Q * 0.2791 ** 2 * 1025;
const EI_HAND = 207e9 * Math.PI / 64 * (0.2731 ** 4 - 0.2319 ** 4), T_HAND = (1.4 - 0.5) * handWeight() * 980;

export default {
  id: 's09-viv', suite: 9, title: 'Vortex-induced vibration',
  summary: 'Natural frequencies, excited modes, cross-flow VIV amplitude, drag amplification and VIV fatigue life of a riser in sheared current, with strake coverage.',
  covers: [5], fidelity: 'screening',
  standards: ['DNV-RP-F204 (2019) §4 VIV fatigue', 'DNV-RP-C205 (2019) §9 vortex-induced oscillations', 'DNV-RP-C203 (2021) §2.4', 'API RP 2RD (1998, R2013) §6.5'],
  models: ['Strouhal vortex-shedding relation', 'Reduced velocity and lock-in window', 'Tensioned-beam modes (closed form / Rayleigh–Ritz) and rod finite-element modes', 'Modal power-in screening in sheared current (Vandiver)', 'Griffin / Sarpkaya response amplitude vs mass-damping', 'Hydrodynamic damping outside the power-in zone', 'Drag amplification (Vandiver)', 'van der Pol wake oscillator (Facchinetti et al. 2004)', 'S–N fatigue with Miner summation over a current scatter'],
  limits: 'Screening-level modal method: cross-flow response only, one mode locked in at a time, constant Strouhal number and lift coefficient, no in-line VIV, higher harmonics, travelling waves or mode competition in time. For the catenary option the in-plane modes of the 2-D rod model stand in for the out-of-plane modes and the full current speed is taken as normal to the pipe. Response and damping coefficients are overridable defaults — confirm them against the governing project edition or a dedicated VIV program; CFD-derived coefficients can be entered but CFD itself is outside the module.',
  publishes: ['riser.vivDamagePerYear', 'riser.vivDragFactor', 'riser.vivAOverD'],
  searchTerms: 'riser vortex induced vibration fatigue strakes sheared current',
  inputs: [
    { key: 'config', label: 'Riser model', type: 'select', default: 'vertical', group: 'Configuration', role: 'geometry',
      options: [{ value: 'vertical', label: 'Vertical tensioned riser (TTR / drilling / tower)' }, { value: 'scr', label: 'Steel catenary riser (finite-element modes, slower)' }] },
    ...pipeInputs(),
    ...waterInputs(),
    ...hydroInputs().filter((i) => i.key === 'Ca' || i.key === 'Cdn'),
    { key: 'length', label: 'Riser length', unit: 'm', type: 'number', default: 980, min: 30, max: 3600, group: 'Configuration', role: 'geometry', showIf: vert },
    { key: 'tensionFactor', label: 'Top tension factor', unit: '-', type: 'number', default: 1.4, min: 1.02, max: 3, group: 'End conditions', role: 'bc', showIf: vert, design: true, help: 'Top tension / apparent weight. Both ends are pinned in the modal model.' },
    { key: 'tensionProfile', label: 'Tension along the riser', type: 'select', default: 'weight', group: 'End conditions', role: 'bc', showIf: vert,
      options: [{ value: 'weight', label: 'Linear with apparent weight' }, { value: 'uniform', label: 'Uniform at the mid-length value (closed-form modes)' }] },
    { key: 'hangoffDepth', label: 'Hang-off depth below still water', unit: 'm', type: 'number', default: 15, min: 0, max: 200, group: 'End conditions', role: 'bc', showIf: scr },
    { key: 'hangoffAngle', label: 'Hang-off angle from vertical', unit: 'deg', type: 'number', default: 12, min: 3, max: 30, group: 'End conditions', role: 'bc', showIf: scr },
    ...currentInputs({ Us: 0.7, Ub: 0.1 }),
    { key: 'curScatter', label: 'Current scatter (scale on the design profile)', type: 'table', group: 'Current', role: 'bc',
      columns: [{ key: 'scale', label: 'Profile scale', unit: '-', type: 'number' }, { key: 'prob', label: 'Occurrence', unit: '%', type: 'number' }],
      default: [{ scale: 0.3, prob: 50 }, { scale: 0.55, prob: 35 }, { scale: 0.8, prob: 12 }, { scale: 1, prob: 3 }],
      help: 'Long-term current statistics as multiples of the design profile. Occurrences should sum to 100 %; any shortfall is treated as no VIV.' },
    { key: 'strakeFrac', label: 'Strake coverage from the top', unit: '-', type: 'number', default: 0.5, min: 0, max: 1, group: 'VIV suppression', role: 'geometry', design: true, help: 'Fraction of the (suspended) length fitted with strakes or fairings, measured down from the top where the current is strongest.' },
    { key: 'strakeEff', label: 'Suppression efficiency', unit: '-', type: 'number', default: 0.9, min: 0, max: 0.99, group: 'VIV suppression', role: 'criteria', help: `Fraction of excitation removed in covered zones (helical strakes 0.8–0.95 when clean; marine growth reduces it). ${CONFIRM}` },
    { key: 'solveStrakes', label: 'Solve the coverage needed for a target amplitude (slower)', type: 'bool', default: false, group: 'VIV suppression', role: 'control' },
    { key: 'targetAD', label: 'Target amplitude A/D', unit: '-', type: 'number', default: 0.1, min: 0.01, max: 1, group: 'VIV suppression', role: 'criteria', showIf: { key: 'solveStrakes', in: [true] } },
    { key: 'St', label: 'Strouhal number', unit: '-', type: 'number', default: 0.2, min: 0.12, max: 0.3, group: 'VIV coefficients', role: 'criteria', uncertain: { dist: 'normal', cov: 0.08 }, help: `0.2 for sub-critical flow on a bare cylinder; 0.17–0.24 depending on Reynolds number and roughness. ${CONFIRM}` },
    { key: 'vrLo', label: 'Power-in reduced velocity, lower bound', unit: '-', type: 'number', default: 4, min: 2, max: 6, group: 'VIV coefficients', role: 'criteria', help: `Lock-in band used to find where each mode is excited. ${CONFIRM}` },
    { key: 'vrHi', label: 'Power-in reduced velocity, upper bound', unit: '-', type: 'number', default: 8, min: 6, max: 16, group: 'VIV coefficients', role: 'criteria' },
    { key: 'Cl', label: 'Lift coefficient (energy balance)', unit: '-', type: 'number', default: 0.7, min: 0.1, max: 1.2, group: 'VIV coefficients', role: 'criteria', uncertain: { dist: 'normal', cov: 0.2 } },
    { key: 'respModel', label: 'Amplitude model', type: 'select', default: 'griffin', group: 'VIV coefficients', role: 'criteria', options: [{ value: 'griffin', label: 'Griffin & Ramberg' }, { value: 'sarpkaya', label: 'Sarpkaya' }] },
    { key: 'zeta', label: 'Structural damping ratio', unit: '-', type: 'number', default: 0.003, min: 0.0005, max: 0.05, group: 'VIV coefficients', role: 'material', uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'Csw', label: 'Still-water damping coefficient', unit: '-', type: 'number', default: 0.25, min: 0, max: 1, group: 'VIV coefficients', role: 'criteria', help: `Venugopal damping model coefficients (still water, low and high reduced velocity). ${CONFIRM}` },
    { key: 'Cvl', label: 'Low-reduced-velocity damping coefficient', unit: '-', type: 'number', default: 0.18, min: 0, max: 1, group: 'VIV coefficients', role: 'criteria' },
    { key: 'Cvh', label: 'High-reduced-velocity damping coefficient', unit: '-', type: 'number', default: 0.2, min: 0, max: 1, group: 'VIV coefficients', role: 'criteria' },
    { key: 'nModes', label: 'Modes computed', type: 'int', default: 20, min: 3, max: 60, group: 'Solver', role: 'control', help: 'Must reach beyond the highest shedding frequency; the result warns when the top mode is excited.' },
    ...snInputs({ curve: 'E' }),
    { key: 'dff', label: 'Design fatigue factor for VIV', unit: '-', type: 'number', default: 10, min: 1, max: 30, group: 'Acceptance criteria', role: 'criteria', help: `VIV fatigue predictions are uncertain: factors of 10–20 are customary. ${CONFIRM}` },
    { key: 'designLife', label: 'Design life', unit: 'yr', type: 'number', default: 25, min: 1, max: 60, group: 'Acceptance criteria', role: 'criteria', bind: 'life.designYears' },
    meshInput(200, '100 / 200 / 400'),
  ],
  mesh: { key: 'nElem', label: 'Stations along the riser', levels: [100, 200, 400], qoi: 'aOverD' },

  compute(inp, ctx) {
    // engineering units → SI: mm → m in pipeFrom; table occurrences % → fraction
    const pp = pipeFrom(inp, ctx), sec = pp.section;
    need(inp.vrHi > inp.vrLo, 'The upper bound of the power-in window must exceed the lower bound.');
    const mass = sec.m + sec.ma, cf = currentFn(inp, inp.depth);
    const isV = inp.config === 'vertical';
    let s, modes, Uof, Lact, zone, from;
    if (isV) {
      const L = inp.length; need(L <= inp.depth + 40, `Riser length ${L} m exceeds the water depth plus deck elevation.`);
      need(sec.w > 0, 'The riser is buoyant; the vertical tensioned model needs positive apparent weight.');
      const T1 = inp.tensionFactor * sec.w * L, T0 = T1 - sec.w * L;
      const args = inp.tensionProfile === 'uniform' ? { T: 0.5 * (T0 + T1) } : { T0, T1 };
      const md = run(() => tensionedBeamModes({ L, EI: sec.EI, m: mass, ...args, nModes: inp.nModes, nPoints: inp.nElem + 1 }), 'Modal analysis failed');
      s = md.s; modes = md.modes; Uof = (x) => cf(x); Lact = L; zone = [L * (1 - inp.strakeFrac), L]; from = 'end';   // s runs upwards from the base
    } else {
      const { base } = buildScr({ ...inp, topType: 'pinned', kn: 1000, mu: 0, pInt: 0, groundedLength: 300 }, ctx, { nElem: inp.nElem });
      const md = run(() => rodModes(base, { nModes: Math.min(inp.nModes, 30) }), 'Modal analysis of the catenary failed');
      need(md.converged, 'The eigenvalue solution did not converge; reduce the number of modes.');
      s = md.s; modes = md.modes; const zs = base.z; Lact = base.tdp ? base.tdp.s : s[s.length - 1];
      const zAt = (x) => { let lo = 0, hi = s.length - 1; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (s[m] <= x) lo = m; else hi = m; } return zs[lo] + (zs[hi] - zs[lo]) * (x - s[lo]) / (s[hi] - s[lo] || 1); };
      Uof = (x) => (x > Lact ? 0 : cf(zAt(x))); zone = [0, Lact * inp.strakeFrac]; from = 'start';
    }
    const { sn, label } = snFrom(inp, pp.t);
    const supp = inp.strakeFrac > 0 && inp.strakeEff > 0 ? { zones: [zone], efficiency: inp.strakeEff } : undefined;
    const scArgs = (scale, sup) => ({ s, modes, U: (x) => scale * Uof(x), D: sec.D, mass, zeta: inp.zeta, rho: pp.rhoW, St: inp.St, vrWindow: [inp.vrLo, inp.vrHi], Cl: inp.Cl, responseModel: inp.respModel, suppression: sup, damping: { Csw: inp.Csw, Cvl: inp.Cvl, Cvh: inp.Cvh } });
    const Nf = (S) => cyclesFromConstants(sn, S);
    const evalCase = (scale, sup) => {
      const sc = run(() => vivModalScreening(scArgs(scale, sup)), 'VIV screening failed (a mode has a vanishing power-in zone for this window; widen the reduced-velocity window or change the number of stations)');
      if (sc.dominant < 0) return { sc, fat: null, dom: null };
      const fat = vivFatigue({ s, modes, screening: sc, E: pp.mat.E, stressDiameter: pp.D, scf: 1, snCyclesToFailure: Nf });
      return { sc, fat, dom: sc.modes.find((m) => m.index === sc.dominant) ?? sc.modes[sc.dominant] };
    };
    const design = evalCase(1, supp), bare = supp ? evalCase(1, undefined) : design;
    const cases = (inp.curScatter || []).filter((c) => c.scale > 0 && c.prob > 0).map((c) => ({ scale: c.scale, p: c.prob / 100 }));
    need(cases.length > 0, 'The current scatter table needs at least one row with positive scale and occurrence.');
    need(cases.reduce((a, c) => a + c.p, 0) <= 1.0001, 'Current occurrences sum to more than 100 %.');
    const dmg = new Float64Array(s.length);
    for (const c of cases) {
      c.res = Math.abs(c.scale - 1) < 1e-9 ? design : evalCase(c.scale, supp);
      if (c.res.fat) for (let k = 0; k < s.length; k++) dmg[k] += c.p * c.res.fat.damagePerYear[k];
    }
    let dMax = 0, kMax = 0; for (let k = 0; k < s.length; k++) if (dmg[k] > dMax) { dMax = dmg[k]; kMax = k; }
    const CAP = 1e6;   // years; reported when no damage is predicted
    const life = dMax > 0 ? Math.min(CAP, 1 / dMax) : CAP, need_ = inp.dff * inp.designLife;
    const dom = design.dom, sigA = design.fat ? Math.max(...design.fat.stressAmplitude) * inp.scf : 0;
    const Umax = Math.max(...Array.from(s, (x) => Uof(x)));
    const Re = Umax * sec.D / 1.19e-6, fs = inp.St * Umax / sec.D, fTop = modes[modes.length - 1].fHz;
    const mStar = sec.m / (pp.rhoW * Math.PI / 4 * sec.D ** 2);
    const urs = []; for (let u = 3; u <= 11.001; u += 0.4) urs.push(u);
    const wo = run(() => wakeOscillator({ massRatio: mStar, zeta: inp.zeta, St: inp.St, Ur: urs, cycles: 120, stepsPerCycle: 30 }), 'Wake-oscillator model failed');
    let cov = null;
    if (inp.solveStrakes) cov = run(() => suppressionCoverage({ ...scArgs(1, undefined), efficiency: inp.strakeEff, targetAOverD: inp.targetAD, from }), 'Strake coverage search failed');
    const aD = dom ? dom.aOverD : 0, dragF = dom ? dom.dragFactor : 1;

    const kpis = [
      kpi('f1', 'Fundamental natural frequency', modes[0].fHz, 'Hz', 'Lowest transverse natural frequency of the riser in water (with added mass). Long risers have very low fundamentals, so currents excite high mode numbers.'),
      kpi('fShed', 'Peak vortex-shedding frequency', fs, 'Hz', 'Strouhal frequency St·U/D at the fastest current along the riser. Modes with natural frequencies up to about this value can lock in.'),
      kpi('modeDom', 'Dominant excited mode number', dom ? dom.index + 1 : 0, '-', 'Mode with the largest input power in the design current (0 = none excited). Higher modes mean shorter wavelengths, more curvature and more cycles.'),
      kpi('fDom', 'Dominant response frequency', dom ? dom.fHz : 0, 'Hz', 'Natural frequency of the dominant mode: the frequency at which the riser vibrates and accumulates fatigue cycles.'),
      kpi('nExcited', 'Number of potentially excited modes', design.sc.excited.length, '-', 'Modes whose lock-in band overlaps the current profile somewhere along the riser. Many candidates indicate multi-mode response in sheared flow.'),
      kpi('aOverD', 'Cross-flow amplitude A/D (design current)', aD, '-', 'Single amplitude of the dominant mode at its antinode in diameters, with the fitted suppression. Bare risers reach about 0.5–1.0; good strake coverage brings it towards 0.1.'),
      kpi('aOverDBare', 'Amplitude A/D without suppression', bare.dom ? bare.dom.aOverD : 0, '-', 'Amplitude the same riser would reach with no strakes or fairings, for comparison with the suppressed value.'),
      kpi('powerIn', 'Power-in length of the dominant mode', dom ? dom.powerInLength : 0, 'm', 'Length over which the local reduced velocity lies in the lock-in band for the dominant mode. The rest of the riser provides hydrodynamic damping.'),
      kpi('zetaEff', 'Effective damping ratio (structural + hydrodynamic)', dom ? dom.zetaEff : inp.zeta, '-', 'Modal damping including hydrodynamic damping outside the power-in zone; it limits the amplitude in sheared current.'),
      kpi('dragFactor', 'Drag amplification factor', dragF, '-', 'Multiplier on the drag coefficient caused by the vibration. Use the amplified drag in the static and interference analyses of the same current case.'),
      kpi('stressAmp', 'Maximum VIV stress amplitude (with SCF)', sigA / 1e6, 'MPa', 'Largest hot-spot bending stress amplitude along the riser for the dominant mode in the design current; the stress range per cycle is twice this.'),
      kpi('damageYr', 'Annual VIV damage (current scatter)', dMax, '1/yr', 'Miner damage per year at the worst section, weighted over the current scatter. Pass it to the fatigue-life module to combine with wave fatigue.'),
      kpi('life', 'VIV fatigue life (unfactored)', life, 'yr', `Years to a Miner sum of 1.0 from VIV alone at the worst section (capped at ${CAP} when no damage is predicted). It must exceed the design life times the VIV fatigue factor.`, need_, 'min'),
      kpi('sCrit', isV ? 'Critical section height above base' : 'Critical section arc length from hang-off', s[kMax], 'm', 'Location of the largest VIV damage along the riser; target it for strakes, weld quality and monitoring.'),
      kpi('wakeAD', 'Wake-oscillator peak amplitude A/D', wo.peakAOverD, '-', 'Peak amplitude of a spring-mounted rigid cylinder with the same mass ratio and damping from the van der Pol wake-oscillator model: an independent order-of-magnitude check on the empirical amplitude.'),
      kpi('reynolds', 'Reynolds number at peak current', Re, '-', 'U·D/ν at the fastest current. Near and above 2×10⁵ (critical regime) the Strouhal number and lift become sensitive to roughness.'),
    ];
    if (cov) kpis.push(kpi('strakeLength', 'Strake coverage for the target amplitude', cov.coveredLength, 'm', `Covered length, grown from the top, that brings the dominant-mode amplitude to the target A/D of ${inp.targetAD}${cov.achievable ? '' : ' (not reached even with full coverage at this efficiency)'}.`));
    const checks = [check('VIV fatigue life with design fatigue factor', need_, life, 'yr', 'DNV-RP-F204 (2019) §4, §6; DNV-ST-F201 (2018) Sec. 5 F', 'Requires the VIV fatigue life to exceed the design life times the fatigue factor. If it does not, vibration at the dominant mode can crack a girth weld within the service life.')];

    const act = design.sc.modes.filter((m) => m.excited).slice(0, 40);
    const [sx, dy] = thin([s, dmg], 300);
    const shape = dom ? thin([s, modes[dom.index]?.phi ?? modes[0].phi], 300) : thin([s, modes[0].phi], 300);
    const up = thin([s, Array.from(s, (x) => Uof(x))], 200);
    const plots = [
      plot('Current profile and suppression zone', isV ? 'Height above base' : 'Arc length from hang-off', 'm', 'Current speed', 'm/s', [{ name: 'Design current', x: up[0], y: up[1] }, ...(supp ? [{ name: 'Strake zone', x: zone, y: [0, 0], dash: true }] : [])], 'Current speed normal to the riser along its length; the dashed bar marks the straked length. Excitation comes from the fast part of the profile, damping from the slow part.'),
      plot('Modal response in the design current', 'Mode natural frequency', 'Hz', 'A/D', '-', [{ name: 'Empirical amplitude', x: act.map((m) => m.fHz), y: act.map((m) => m.aOverD) }, { name: 'Energy-balance amplitude', x: act.map((m) => m.fHz), y: act.map((m) => m.aOverDEnergy), dash: true }], 'Amplitude each candidate mode would reach if it alone locked in. The dominant mode is the one with the highest input power, usually near the shedding frequency of the fastest current.', 'scatter'),
      plot('Dominant mode shape', isV ? 'Height above base' : 'Arc length from hang-off', 'm', 'Mode shape', '-', [{ name: dom ? `Mode ${dom.index + 1}` : 'Mode 1', x: shape[0], y: shape[1] }], 'Shape of the dominant mode (unit amplitude). Bending stress peaks at the points of greatest curvature of this shape.'),
      plot('Annual VIV damage along the riser', isV ? 'Height above base' : 'Arc length from hang-off', 'm', 'Damage per year', '1/yr', [{ name: 'Scatter-weighted damage', x: sx, y: dy }, { name: 'Allowable 1/(DFF·life)', x: [s[0], s[s.length - 1]], y: [1 / need_, 1 / need_], dash: true }], 'Miner damage per year from VIV, summed over the current scatter. Sections above the dashed line fail the factored-life requirement.'),
      plot('Wake-oscillator response', 'Reduced velocity', '-', 'A/D', '-', [{ name: 'van der Pol wake oscillator', x: wo.Ur, y: wo.aOverD }, { name: 'Empirical amplitude (bare)', x: [wo.Ur[0], wo.Ur[wo.Ur.length - 1]], y: [bare.dom ? bare.dom.aOverD : 0, bare.dom ? bare.dom.aOverD : 0], dash: true }], 'Amplitude of a rigid cylinder on springs against reduced velocity from the wake-oscillator model. The lock-in hump should sit inside the power-in window and its height should be of the same order as the empirical amplitude.'),
    ];
    const tables = [
      { title: 'Candidate modes (design current)', columns: [{ key: 'n', label: 'Mode' }, { key: 'f', label: 'Frequency', unit: 'Hz' }, { key: 'Lin', label: 'Power-in length', unit: 'm' }, { key: 'zeta', label: 'Damping ratio', unit: '-' }, { key: 'SG', label: 'Skop–Griffin SG', unit: '-' }, { key: 'aD', label: 'A/D', unit: '-' }, { key: 'P', label: 'Input power', unit: 'W' }, { key: 'Cd', label: 'Drag factor', unit: '-' }],
        rows: act.map((m) => ({ n: m.index + 1, f: round(m.fHz, 4), Lin: round(m.powerInLength, 4), zeta: round(m.zetaEff, 3), SG: round(m.SG, 3), aD: round(m.aOverD, 3), P: round(m.power, 4), Cd: round(m.dragFactor, 3) })),
        explain: 'Every mode that can lock in somewhere along the riser in the design current, with its power-in length, damping, amplitude and input power. The dominant mode has the largest power.' },
      { title: 'Current scatter cases', columns: [{ key: 'scale', label: 'Profile scale', unit: '-' }, { key: 'p', label: 'Occurrence', unit: '%' }, { key: 'mode', label: 'Dominant mode' }, { key: 'aD', label: 'A/D', unit: '-' }, { key: 'd', label: 'Max damage rate', unit: '1/yr' }, { key: 'share', label: 'Weighted share at critical section', unit: '%' }],
        rows: cases.map((c) => ({ scale: c.scale, p: round(c.p * 100, 4), mode: c.res.dom ? c.res.dom.index + 1 : 0, aD: round(c.res.dom ? c.res.dom.aOverD : 0, 3), d: Number((c.res.fat ? c.res.fat.maxDamagePerYear : 0).toExponential(3)), share: round(dMax > 0 && c.res.fat ? 100 * c.p * c.res.fat.damagePerYear[kMax] / dMax : 0, 3) })),
        explain: 'One row per current condition: which mode it excites, the amplitude, the damage rate while it lasts and its share of the long-term damage at the critical section.' },
    ];
    const uc = checks[0].uc, topExcited = dom && dom.index >= modes.length - 2;
    const strakeLen = inp.strakeFrac * Lact, perM = Math.PI / 4 * ((sec.D * 1.5) ** 2 - sec.D ** 2) * 0.35 * 1150 / 1e3;   // t of polymer per metre of strake shell (screening)
    const advice = [
      adv('engineering', uc > 1 ? 'fail' : uc > 0.5 || topExcited ? 'warn' : 'ok',
        `Design current excites mode ${dom ? dom.index + 1 : '–'} at ${fmt(dom ? dom.fHz : 0)} Hz with A/D = ${fmt(aD)} (bare ${fmt(bare.dom ? bare.dom.aOverD : 0)}); VIV life ${fmt(life)} yr against ${fmt(need_)} yr required; drag amplified by ${fmt(dragF)}.`,
        topExcited ? `The highest computed mode (${modes.length}, ${fmt(fTop)} Hz) is excited while shedding reaches ${fmt(fs)} Hz: increase "Modes computed" until the dominant mode sits below the top of the set.`
          : uc > 1 ? `Extend strake coverage beyond ${fmt(inp.strakeFrac * 100)} % of the length (tick the coverage solver for the length that reaches a target amplitude), raise top tension to push modes up, or upgrade the weld class at ${fmt(s[kMax])} m.`
            : `Carry the amplified drag coefficient ${fmt(inp.Cdn * dragF)} into the static, interference and wave analyses for this current, and confirm the result with a dedicated VIV program before deleting any strakes.`,
        `DNV-RP-F204 (2019) §4 (VIV fatigue assessment) and DNV-RP-C205 (2019) §9. ${CONFIRM}`),
      adv('economic', 'info',
        `${fmt(strakeLen)} m of strakes are assumed (${fmt(inp.strakeFrac * 100)} % of the length), about ${fmt(strakeLen * perM)} t of polymer shells that must be fitted joint by joint on the lay vessel and kept clean of marine growth. They also raise drag (Cd ≈ 1.4–2.0), which increases vessel offset loads.`,
        'Optimise coverage against VIV life rather than covering by rule of thumb: strakes are paid for per metre, slow the lay rate, and need cleaning or replacement when fouled. Compare with fairings where the drag penalty matters.',
        'Whole-life costing (ISO 15686-5) including installation rate and in-service cleaning; inspection planning per DNV-RP-G101.',
        `Each 100 m of strake avoided saves about ${fmt(100 * perM)} t of polyurethane/polyethylene (no established recycling route offshore). A riser that needs a thicker wall instead: ${steelCo2Text(steelTonnes(pp.D, 0.001, Lact, pp.mat.rho))} per millimetre.`),
    ];
    return finish({
      kpis, checks, plots, tables,
      notes: [`${modes.length} modes from ${isV ? (inp.tensionProfile === 'uniform' ? 'the closed-form tensioned beam' : 'Rayleigh–Ritz with linearly varying tension') : 'the catenary rod model (in-plane modes)'}; ${design.sc.excited.length} candidates in the design current; ${label}.`,
        dMax > 0 ? `Critical section at ${fmt(s[kMax])} m; strakes over ${fmt(inp.strakeFrac * 100)} % at ${fmt(inp.strakeEff * 100)} % efficiency.` : 'No VIV fatigue damage is predicted for the given current scatter.'],
      advice,
      publish: { 'riser.vivDamagePerYear': dMax, 'riser.vivDragFactor': dragF, 'riser.vivAOverD': aD },
      quantities: { polymer_t: strakeLen * perM },
    });
  },

  verify: [
    {
      // Pinned–pinned tensioned beam, uniform tension T = (TF − ½)·w·L = 0.9·w·980 N, mass per length incl. added mass
      // m = steel + coating + contents + Ca·ρ·π/4·Do² (hand formula above):
      //   f1 = (1/2L)·sqrt(T/m)·sqrt(1 + π²·EI/(T·L²)).
      name: 'Tensioned-beam fundamental f1 = (1/2L)·√(T/m)·√(1 + π²EI/(TL²))', source: 'Blevins, Formulas for Natural Frequency and Mode Shape (1979), tensioned beam; DNV-RP-F204 App. A',
      inputs: { config: 'vertical', tensionProfile: 'uniform' },
      expect: [{ kpi: 'f1', value: 1 / (2 * 980) * Math.sqrt(T_HAND / M_HAND) * Math.sqrt(1 + Math.PI ** 2 * EI_HAND / (T_HAND * 980 ** 2)), tol: 0.002 }],
    },
    {
      // Strouhal relation: fs = St·U/D with U = 0.7 m/s at the surface (riser top within the water column),
      // D = 0.2791 m hydrodynamic diameter → fs = 0.2·0.7/0.2791 = 0.5016 Hz.
      name: 'Shedding frequency fs = St·U/D', source: 'DNV-RP-C205 (2019) §9.1.1',
      inputs: { config: 'vertical' },
      expect: [{ kpi: 'fShed', value: 0.2 * (0.1 + 0.6 * (980 / 1000) ** 1.5) / 0.2791, tol: 0.002 }],
    },
  ],
};
