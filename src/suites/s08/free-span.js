// Free-span assessment in the DNV-RP-F105 format: natural frequencies, static response, onset screening, in-line and
// cross-flow vortex-induced vibration amplitudes and fatigue, a simplified direct-wave-loading fatigue estimate, the
// ultimate bending stress and the allowable span length.
import { concreteStiffnessFactor, addedMassCoefficient, spanResponse, stabilityParameter, inlineResponse, crossflowResponse,
  vivScreening, vivFatigue, spanUls, maxSpanLength, boundaryCoefficients } from '../../engine/pipeline/span.js';
import { soilStiffness } from '../../engine/pipeline/pipesoil.js';
import { snCurve, cyclesToFailure } from '../../engine/fatigue/sn.js';
import { wakeOscillator } from '../../engine/lines/viv.js';
import { pipeInputs, pipeFrom, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, linspace, sev, worst, clamp,
  vesselCo2, CONFIRM, LIMIT_NO_FE, RBI_BASIS, MPA, YEAR } from './_common.js';

const G_G = 'Span geometry and support';
const G_F = 'Flow at pipe level';
const G_D = 'Damping and fatigue';
const G_C = 'Acceptance criteria';
const CLASSES = ['sand-loose', 'sand-medium', 'sand-dense', 'clay-very-soft', 'clay-soft', 'clay-firm', 'clay-stiff', 'clay-very-stiff', 'clay-hard'];
const ETA = { low: 1.0, normal: 0.5, high: 0.25 }, GK = { low: 1.0, normal: 1.15, high: 1.3 };

export default {
  id: 's08-free-span',
  suite: 8,
  title: 'Free-span and vortex-induced vibration',
  summary: 'Finds the natural frequencies of a free span, whether it will vibrate in the current and waves, the fatigue life if it does, and the longest span that can be left uncorrected.',
  covers: [4],
  fidelity: 'analytical',
  standards: ['DNV-RP-F105 (2006/2017) §2.3–2.5, §4.3–4.4, §6.7, Table 6-1', 'DNV-RP-C203 (2021) S–N curves', 'DNV-RP-C205 (2019) §6 (Morison loads)'],
  models: ['Free-span beam vibration on elastic shoulders (effective span length)', 'Strouhal / reduced-velocity onset of vortex shedding', 'In-line VIV response model', 'Cross-flow VIV response model', 'Direct wave loading (quasi-static Morison with dynamic amplification)', 'Wake-oscillator model (Facchinetti et al.) as a cross-check of the cross-flow amplitude', 'S–N fatigue with Palmgren–Miner summation over a Weibull current distribution'],
  limits: `Single span, first mode in each direction, with the DNV-RP-F105 (2006) response-model points and soil stiffness tables as recalled — expose and confirm them against the edition in force. Cross-flow-induced in-line motion, multi-mode and multi-span interaction are not included; direct wave loading is a simplified quasi-static Morison estimate, not the RP force model; the three fatigue contributions are summed. The wake oscillator is a rigid-cylinder single-degree-of-freedom cross-check, not used for acceptance; CFD is outside this module. ${LIMIT_NO_FE}`,
  publishes: ['pipeline.spanMax', 'life.fatigueYears'],
  searchTerms: 'pipeline free span vortex induced vibration DNV-RP-F105 fatigue',
  inputs: [
    ...pipeInputs(),
    num('L', 'Span length', 'm', 24, 3, 300, G_G, 'geometry', 'Free length between the touchdown points from the survey.', { design: true, uncertain: { dist: 'normal', cov: 0.1 } }),
    num('gap', 'Gap under the pipe', 'm', 0.5, 0.01, 20, G_G, 'geometry', 'Average clearance between the pipe and the seabed over the central third of the span.'),
    sel('boundary', 'Span boundary condition', 'seabed', [['seabed', 'Single span on elastic seabed shoulders'], ['pinned', 'Pinned–pinned'], ['fixed', 'Fixed–fixed']], G_G, 'bc', 'Seabed shoulders use the effective span length concept of DNV-RP-F105 §6.7.'),
    sel('soilClass', 'Shoulder soil class', 'sand-medium', CLASSES, G_G, 'soil', `Selects the dynamic and static soil stiffness (F105 Tables 7-5/7-6). ${CONFIRM}`),
    num('Seff', 'Effective axial force (tension positive)', 'kN', -150, -20000, 20000, G_G, 'bc', 'Compression lowers the frequencies and can buckle the span; take it from the expansion analysis at the span location.', { bind: 'pipeline.effectiveForce', uncertain: { dist: 'normal', cov: 0.15 } }),
    num('kc', 'Concrete stiffness factor constant k_c', '-', 0.25, 0, 0.5, G_G, 'material', '0.33 for asphalt, 0.25 for PP/PE anti-corrosion coating under the concrete; zero ignores the concrete stiffness.'),
    num('Uc100', '100-year current at pipe level', 'm/s', 0.6, 0, 3, G_F, 'bc', 'Extreme current normal to the pipe at the span elevation.', { bind: 'site.UcBed100', uncertain: { dist: 'gumbel', cov: 0.15 } }),
    num('ucScale', 'Weibull scale of the long-term current', 'm/s', 0.18, 0.01, 2, G_F, 'bc', 'Scale parameter of the current speed distribution at pipe level, normal to the pipe.'),
    num('ucShape', 'Weibull shape of the long-term current', '-', 1.8, 0.5, 5, G_F, 'bc', 'Shape parameter of the current distribution.'),
    num('Uw1', '1-year significant wave-induced velocity at pipe level', 'm/s', 0.25, 0, 3, G_F, 'bc', 'Significant near-bed wave velocity normal to the pipe in the 1-year sea state.'),
    num('Uwf', 'Fatigue sea-state wave-induced velocity', 'm/s', 0.08, 0, 2, G_F, 'bc', 'Representative significant wave-induced velocity for fatigue (long-term weighted).'),
    num('Tu', 'Wave-induced velocity period', 's', 10, 3, 25, G_F, 'bc', 'Zero-crossing period of the near-bed wave velocity.'),
    num('Uw100', '100-year wave-induced velocity at pipe level', 'm/s', 0.5, 0, 4, G_F, 'bc', 'For the ultimate limit state.'),
    num('zetaStr', 'Structural damping ratio', '-', 0.005, 0.001, 0.05, G_D, 'material', '0.005 without data; 0.01–0.02 with concrete coating.'),
    num('zetaSoil', 'Soil damping ratio', '-', 0.01, 0, 0.05, G_D, 'soil', 'Typically 0.005–0.02 depending on soil and span length.'),
    sel('sn', 'S–N curve of the girth weld', 'F1', ['C1', 'D', 'E', 'F', 'F1', 'F3'], G_D, 'material', 'DNV-RP-C203 curve for the weld root or cap; F1 is common for single-sided girth welds.'),
    sel('snEnv', 'S–N environment', 'seawater-cp', [['seawater-cp', 'Sea water with cathodic protection'], ['air', 'Air (internal, dry)'], ['free-corrosion', 'Free corrosion']], G_D, 'environment', 'Environment at the crack initiation site.'),
    num('scf', 'Stress concentration factor', '-', 1.2, 1, 3, G_D, 'material', 'Girth-weld misalignment and geometry.'),
    num('life', 'Exposure time of the span', 'years', 25, 0.1, 100, G_D, 'operation', 'Design life for a permanent span, or the duration of a temporary phase.', { bind: 'life.designYears' }),
    sel('safetyClass', 'Safety class', 'normal', [['low', 'Low'], ['normal', 'Normal'], ['high', 'High']], G_C, 'criteria', 'Sets the allowable damage ratio η (1.0 / 0.5 / 0.25) and γ_k (1.0 / 1.15 / 1.3).'),
    num('gammaS', 'Stress range factor γ_s', '-', 1.3, 1, 1.6, G_C, 'criteria', `F105 Table 2-2. ${CONFIRM}`),
    num('gammaOnIL', 'Onset factor in-line γ_on,IL', '-', 1.1, 1, 1.5, G_C, 'criteria', `F105 Table 2-2. ${CONFIRM}`),
    num('gammaOnCF', 'Onset factor cross-flow γ_on,CF', '-', 1.2, 1, 1.5, G_C, 'criteria', `F105 Table 2-2. ${CONFIRM}`),
    num('gammaScr', 'Screening factor on frequency', '-', 1.4, 1, 2, G_C, 'criteria', `F105 Table 2-1 (γ_IL = γ_CF = 1.4). ${CONFIRM}`),
    num('etaUls', 'Allowable static + environmental bending stress / SMYS', '-', 0.72, 0.3, 1, G_C, 'criteria', 'Simplified ULS acceptance; use the combined-loading module for the code check.'),
  ],
  compute(inp, ctx) {
    const p = pipeFrom(inp, ctx);
    need(p.ws > 0, 'The pipe is buoyant in operation: a free span cannot sag under its own weight.');
    const D = p.Do, L = inp.L, eD = inp.gap / D;
    const Adisp = Math.PI * D * D / 4;
    const Ca = addedMassCoefficient(eD);
    const me = p.massOp + Ca * p.rhoW * Adisp;
    const massRatio = p.massOp / (p.rhoW * Adisp);
    const conc = (ctx?.material && ctx.material(inp.concrete)) || { E: 30e9 };
    let CSF = 0;
    if (inp.tConc > 0 && inp.kc > 0) {
      const ri = D / 2 - inp.tConc / 1e3;
      CSF = concreteStiffnessFactor({ EIconc: (conc.E || 30e9) * Math.PI / 4 * ((D / 2) ** 4 - ri ** 4), EIsteel: p.EI, kc: inp.kc });
    }
    const ks = run(() => soilStiffness({ D, massRatio: clamp(massRatio, 1.2, 2.0), soilType: inp.soilClass }), 'Soil stiffness');
    const spanArgs = { D, Ds: p.D, t: p.t, EI: p.EI, me, E: p.E, Seff: inp.Seff * 1e3, q: p.ws, CSF, boundary: inp.boundary, KV: ks.KV, KL: ks.KL, KVS: ks.KVS ?? ks.KV };
    const rs = run(() => spanResponse({ ...spanArgs, L }), 'Span response (the span buckles as a column under this compression — shorten it or reduce the axial force)');
    const zetaT = inp.zetaStr + inp.zetaSoil;
    const stab = stabilityParameter({ me, zetaT, D, rhoW: p.rhoW, gammaK: GK[inp.safetyClass] });
    const VRil = (U) => U / (rs.f1IL * D), VRcf = (U) => U / (rs.f1CF * D);
    const alpha100 = inp.Uc100 / Math.max(inp.Uc100 + inp.Uw1, 1e-9);
    const KC = inp.Uw1 > 0 ? inp.Uw1 * inp.Tu / D : Infinity;
    const cfOpt = { Ksd: stab.Ksd, alpha: alpha100, KC, gapRatio: eD, massRatio, gammaOn: inp.gammaOnCF };
    const il0 = inlineResponse({ VR: VRil(inp.Uc100 + inp.Uw1), Ksd: stab.Ksd, gammaOn: inp.gammaOnIL, alpha: alpha100 });
    const cf0 = crossflowResponse({ VR: VRcf(inp.Uc100 + inp.Uw1), ...cfOpt });
    const scr = vivScreening({ f1IL: rs.f1IL, f1CF: rs.f1CF, Uc100: inp.Uc100, Uw1: inp.Uw1, D, L, VRonsetIL: il0.VRonset, VRonsetCF: cf0.VRonset, gammaIL: inp.gammaScr, gammaCF: inp.gammaScr });
    const curve = snCurve(inp.sn, inp.snEnv);
    const N = (S) => cyclesToFailure({ curve, stressRange: Math.max(S, 1), thickness: p.t, scf: inp.scf });
    const fatArgs = { D, Ksd: stab.Ksd, current: { scale: inp.ucScale, shape: inp.ucShape }, sn: N, gammaS: inp.gammaS, gammaOnIL: inp.gammaOnIL, gammaOnCF: inp.gammaOnCF, gapRatio: eD, massRatio };
    if (inp.Uwf > 0) fatArgs.seaStates = [{ prob: 1, Uw: inp.Uwf, Tu: inp.Tu }];
    const fat = run(() => vivFatigue({ ...fatArgs, f1IL: rs.f1IL, f1CF: rs.f1CF, AIL: rs.AIL, ACF: rs.ACF }), 'VIV fatigue');
    // direct wave loading: Morison force range on the span as a quasi-static beam load with single-mode amplification
    const waveStress = (r, Uw, Uc) => {
      if (!(Uw > 0)) return 0;
      const CD = 1.0, CM = 1 + Ca, om = 2 * Math.PI / inp.Tu;
      const drag = (u) => 0.5 * p.rhoW * D * CD * u * Math.abs(u);
      const dF = Math.hypot(drag(Uc + Uw) - drag(Uc - Uw), 2 * p.rhoW * Adisp * CM * om * Uw);
      const bc = boundaryCoefficients(inp.boundary, r.LeffIL / L);
      const daf = 1 / Math.max(Math.abs(1 - (1 / inp.Tu / r.f1IL) ** 2), 0.05);
      return bc.C5 * dF * r.LeffIL ** 2 * p.D / (2 * p.I) * daf;
    };
    const ucMean = inp.ucScale * 0.89;
    const Swave = waveStress(rs, inp.Uwf, ucMean) * inp.gammaS;
    const dWave = Swave > 0 ? (YEAR / inp.Tu) / N(Swave) : 0;
    const dTot = fat.damagePerYearIL + fat.damagePerYearCF + dWave;
    const lifeTot = dTot > 0 ? 1 / dTot : 1e9;
    const eta = ETA[inp.safetyClass];
    const ucFat = dTot * inp.life / eta;
    const sigAllow = inp.etaUls * p.smys;
    const uls = run(() => spanUls({ response: rs, D, Uextreme: inp.Uc100 + inp.Uw100, alpha: inp.Uc100 / Math.max(inp.Uc100 + inp.Uw100, 1e-9), Ksd: stab.Ksd, KC: inp.Uw100 > 0 ? inp.Uw100 * inp.Tu / D : Infinity,
      sigmaDrag: 0.5 * waveStress(rs, inp.Uw100, inp.Uc100), gammaS: inp.gammaS, sigmaAllow: sigAllow }), 'Ultimate limit state');
    const maxL = (criterion, extra) => { try { return maxSpanLength({ criterion, span: spanArgs, Lmin: 2, Lmax: 300, ...extra }).L; } catch { return NaN; } };
    const LmaxStatic = maxL('static', { sigmaAllow: sigAllow });
    const LmaxIL = maxL('onset-il', { U: inp.Uc100 + inp.Uw1, VRonset: il0.VRonset, gammaF: inp.gammaScr });
    const LmaxCF = maxL('onset-cf', { U: inp.Uc100 + inp.Uw1, VRonset: cf0.VRonset, gammaF: inp.gammaScr });
    // span length sweep: frequencies and fatigue life
    const Ls = linspace(Math.max(0.4 * L, 5), 1.8 * L, 15);
    const sweep = Ls.map((l) => {
      try {
        const r = spanResponse({ ...spanArgs, L: l });
        const f = vivFatigue({ ...fatArgs, f1IL: r.f1IL, f1CF: r.f1CF, AIL: r.AIL, ACF: r.ACF });
        const sw = waveStress(r, inp.Uwf, ucMean) * inp.gammaS;
        const d = f.damagePerYearIL + f.damagePerYearCF + (sw > 0 ? (YEAR / inp.Tu) / N(sw) : 0);
        return { l, fIL: r.f1IL, fCF: r.f1CF, life: Math.min(d > 0 ? 1 / d : 1e6, 1e6) };
      } catch { return null; }
    }).filter(Boolean);
    let LmaxFat = sweep.length ? sweep[sweep.length - 1].l : L;
    const target = inp.life / eta;
    for (let i = 0; i < sweep.length; i++) if (sweep[i].life < target) { LmaxFat = i === 0 ? sweep[0].l : sweep[i - 1].l + (sweep[i].l - sweep[i - 1].l) * (Math.log(sweep[i - 1].life) - Math.log(target)) / (Math.log(sweep[i - 1].life) - Math.log(sweep[i].life)); break; }
    const Lallow = Math.min(...[LmaxStatic, LmaxFat].filter(Number.isFinite));
    const vr = linspace(0, 17, 86);
    const ilC = vr.map((v) => inlineResponse({ VR: v, Ksd: stab.Ksd, gammaOn: inp.gammaOnIL }).AyOverD);
    const cfC = vr.map((v) => { const c = crossflowResponse({ VR: v, ...cfOpt }); return c.AzOverD * c.Rk; });
    // wake-oscillator cross-check of the cross-flow lock-in amplitude (reduced velocity on the still-water frequency, St = 0.2)
    let wo = null;
    try { wo = wakeOscillator({ massRatio: Math.max(massRatio, 0.3), zeta: zetaT, Ur: linspace(2, 12, 21), cycles: 120, stepsPerCycle: 30 }); } catch { wo = null; }
    const checks = [
      chk('Fatigue damage over the exposure time', dTot * inp.life, eta, '-', 'DNV-RP-F105 §2.4 (η by safety class); DNV-RP-C203',
        'Guards against a fatigue crack at the girth weld from in-line VIV, cross-flow VIV and direct wave loading accumulated over the exposure time.'),
      chk('Ultimate bending stress (static + environmental)', uls.sigmaTotal, sigAllow, 'Pa', 'DNV-RP-F105 §2.5 (simplified stress format)',
        'Guards against overstressing the span shoulders or mid-span under sagging plus the extreme environmental bending; exceeding it calls for the local buckling check with the combined moment.'),
      chk('Span slenderness within validity (L/D ≤ 140)', L / p.D, 140, '-', 'DNV-RP-F105 §1.6 / §6.7',
        'The beam response expressions are valid for spans up to about 140 diameters; longer spans behave as cables and need a finite-element model.'),
    ];
    const gov = worst(checks);
    const extra = Math.max(0, L - Lallow);
    return {
      kpis: [
        kpi('f1IL', 'In-line natural frequency', rs.f1IL, 'Hz', 'First horizontal bending frequency of the span including the effective axial force and shoulder flexibility.'),
        kpi('f1CF', 'Cross-flow natural frequency', rs.f1CF, 'Hz', 'First vertical bending frequency; slightly higher than in-line because the sagged shape stiffens the span.'),
        kpi('lifeTotal', 'Fatigue life (all mechanisms, unfactored)', Math.min(lifeTot, 1e9), 'years', 'One over the summed annual damage from in-line VIV, cross-flow VIV and direct wave loading. It must exceed the exposure time divided by η.', target, 'min'),
        kpi('lifeIL', 'In-line VIV fatigue life', Math.min(fat.lifeILYears, 1e9), 'years', 'Life from in-line vortex-induced vibration alone.'),
        kpi('lifeCF', 'Cross-flow VIV fatigue life', Math.min(fat.lifeCFYears, 1e9), 'years', 'Life from cross-flow vortex-induced vibration alone.'),
        kpi('lifeWave', 'Direct wave loading fatigue life', dWave > 0 ? Math.min(1 / dWave, 1e9) : 1e9, 'years', 'Life from the simplified quasi-static Morison estimate of direct wave loading alone.'),
        kpi('ucFatigue', 'Fatigue utilisation', ucFat, '-', 'Damage over the exposure time divided by the allowable damage ratio η.', 1, 'max'),
        kpi('scrIL', 'In-line screening ratio', scr.ucIL, '-', 'Required over available in-line frequency in the DNV-RP-F105 screening criterion (Eq. 2.1). Below 1.0 in-line VIV fatigue need not be calculated; above 1.0 the fatigue result decides.'),
        kpi('scrCF', 'Cross-flow screening ratio', scr.ucCF, '-', 'Required over available cross-flow frequency in the screening criterion (Eq. 2.2). Above 1.0 cross-flow VIV is possible at the extreme flow and the fatigue result decides.'),
        kpi('AyD', 'In-line amplitude at the 100-year flow', il0.AyOverD * il0.psiAlpha, 'D', 'In-line vibration amplitude in pipe diameters at the 100-year current plus 1-year wave velocity.'),
        kpi('AzD', 'Cross-flow amplitude at the 100-year flow', cf0.AzOverD * cf0.Rk, 'D', 'Cross-flow vibration amplitude in pipe diameters at the same flow, after the damping reduction.'),
        kpi('AzWake', 'Peak cross-flow amplitude, wake oscillator', wo ? wo.peakAOverD : 0, 'D', 'Largest lock-in amplitude of a rigid cylinder with the same mass ratio and damping from the wake-oscillator model: an independent check on the response-model plateau, not used for acceptance.'),
        kpi('VRil', 'In-line reduced velocity at the 100-year flow', VRil(inp.Uc100 + inp.Uw1), '-', 'U/(f·D). In-line VIV starts near 1.0 and cross-flow near 3.'),
        kpi('Ksd', 'Design stability parameter', stab.Ksd, '-', 'Mass-damping parameter 4π·m_e·ζ/(ρD²) divided by γ_k; higher values suppress in-line VIV.'),
        kpi('sigStatic', 'Static bending stress', rs.stressStatic / MPA, 'MPa', 'Bending stress at the outer steel fibre from the sag of the span under submerged weight.'),
        kpi('ucUls', 'Ultimate stress utilisation', uls.uc ?? 0, '-', 'Static plus extreme environmental bending stress over the simplified allowable.', 1, 'max'),
        kpi('deflection', 'Static mid-span deflection', rs.deflection, 'm', 'Sag under the operating submerged weight; if it exceeds the gap the span touches down and shortens.'),
        kpi('Leff', 'Effective span length (cross-flow)', rs.LeffCF, 'm', 'Length of the fixed-fixed beam with the same frequency; longer than the visible span because the shoulders yield.'),
        kpi('Lallow', 'Allowable span length', Lallow, 'm', 'Longest span that meets both the fatigue life with the safety factor and the static stress limit; spans beyond it need support or correction.'),
      ],
      checks,
      plots: [
        plot('line', 'VIV response models', 'Reduced velocity', '-', 'Amplitude / D', '-', [
          { name: 'In-line A_Y/D', x: vr.map((v) => round(v)), y: ilC.map((v) => round(v)) }, { name: 'Cross-flow A_Z/D', x: vr.map((v) => round(v)), y: cfC.map((v) => round(v)) },
          ...(wo ? [{ name: 'Wake oscillator (cross-flow, rigid cylinder)', x: wo.Ur.map((v) => round(v)), y: wo.aOverD.map((v) => round(v)), dash: true }] : []),
          { name: 'In-line at 100-year flow', x: [round(VRil(inp.Uc100 + inp.Uw1))], y: [round(il0.AyOverD)] }, { name: 'Cross-flow at 100-year flow', x: [round(VRcf(inp.Uc100 + inp.Uw1))], y: [round(cf0.AzOverD * cf0.Rk)] },
        ], 'Amplitude response of the span against reduced velocity. Points on the zero part of a curve mean no vibration of that type at the extreme flow; points on the plateau mean full lock-in.'),
        plot('line', 'Natural frequency vs span length', 'Span length', 'm', 'Frequency', 'Hz', [
          { name: 'In-line', x: sweep.map((s) => round(s.l)), y: sweep.map((s) => round(s.fIL)) }, { name: 'Cross-flow', x: sweep.map((s) => round(s.l)), y: sweep.map((s) => round(s.fCF)) },
        ], 'Frequencies fall roughly with the square of the span length; where they drop below the shedding frequency of the design current the span starts to vibrate.'),
        plot('line', 'Fatigue life vs span length', 'Span length', 'm', 'log10 of fatigue life in years', '-', [
          { name: 'All mechanisms', x: sweep.map((s) => round(s.l)), y: sweep.map((s) => round(Math.log10(s.life))) }, { name: 'Required (exposure / η)', x: [round(Ls[0]), round(Ls[Ls.length - 1])], y: [round(Math.log10(target)), round(Math.log10(target))], dash: true },
        ], 'Fatigue life collapses over a narrow range of span length. The crossing with the required life is the allowable span for fatigue; values of 6 mean no damage.'),
      ],
      tables: [tbl('Allowable span length by criterion', [col('crit', 'Criterion'), col('L', 'Allowable span', 'm'), col('basis', 'Basis')], [
        { crit: 'Static bending stress', L: round(LmaxStatic), basis: `stress ≤ ${fmt(inp.etaUls)} SMYS` }, { crit: 'Fatigue life', L: round(LmaxFat), basis: `life ≥ ${fmt(target)} years` },
        { crit: 'No in-line VIV onset', L: round(LmaxIL), basis: 'screening frequency criterion' }, { crit: 'No cross-flow VIV onset', L: round(LmaxCF), basis: 'screening frequency criterion' },
      ].filter((r) => Number.isFinite(r.L)), 'Longest span meeting each criterion. The onset rows are conservative screening limits; the fatigue and static rows are the acceptance limits.')],
      notes: [`Effective mass ${fmt(me)} kg/m (added-mass coefficient ${fmt(Ca)}), concrete stiffness factor ${fmt(CSF)}, Seff/Pcr = ${fmt(inp.Seff * 1e3 / rs.PcrCF)}.`,
        `Screening: in-line ${scr.ilOk ? 'passes' : 'fails'}, cross-flow ${scr.cfOk ? 'passes' : 'fails'}. Governing acceptance check: ${gov.name} (${fmt(gov.uc)}).`,
        ...rs.warnings.map((w) => `Validity: ${w}.`), 'Direct wave loading is a simplified quasi-static estimate; the three damage contributions are added.'],
      advice: [
        adv('engineering', sev(gov.uc), `The ${fmt(L)} m span has a fatigue life of ${fmt(lifeTot)} years against ${fmt(target)} required and an ultimate stress utilisation of ${fmt(uls.uc ?? 0)}; allowable span ${fmt(Lallow)} m.`,
          gov.uc > 1 ? `Shorten the span to ${fmt(Lallow)} m or less with a mid-span support (grout bag or mechanical support), or trench the shoulders; ${fmt(extra)} m must be removed.` : `The span can stay. Re-check after each survey if it grows beyond ${fmt(Lallow)} m, and confirm the current distribution at pipe level.`,
          'DNV-RP-F105 (2017) §2.3–2.5 with safety factors of Table 2-2 (recalled — confirm)'),
        adv('economic', 'info', gov.uc > 1 ? 'The span needs one intervention to bring it within the allowable length.' : `The span is ${fmt(Lallow - L)} m shorter than the allowable length, so no correction is needed now.`,
          gov.uc > 1 ? 'One grout-bag support is typically half a vessel day; batch all span corrections along the route into a single campaign after the as-laid survey.' : `Set the inspection interval so that scour growth of ${fmt(Lallow - L)} m is detected: the cost is survey time, not intervention.`,
          RBI_BASIS, `${vesselCo2(0.5)} Avoiding unnecessary supports also keeps grout and rock off the seabed.`),
      ],
      publish: { 'pipeline.spanMax': Lallow, 'life.fatigueYears': Math.min(lifeTot, 1e9) },
      quantities: { vessel_days: { spanCorrection: gov.uc > 1 ? 0.5 : 0 } },
    };
  },
  verify: [{
    // Pinned–pinned span, no axial force: f1 = C1·√(EI/(m_e L⁴)), C1 = 1.57 (≈ π/2). Bare 400 × 20 mm pipe, empty:
    // EI = 207e9·π/64·(0.4⁴ − 0.36⁴) = 8.9457e7 N·m²; steel 187.4 kg/m + added mass 1025·π/4·0.4² = 128.8 kg/m (C_a = 1 for e/D ≥ 0.8)
    // → m_e = 316.2 kg/m. L = 40 m: f1 = 1.57·√(8.9457e7/(316.2·40⁴)) = 1.57·0.3324 = 0.522 Hz.
    name: 'Pinned–pinned span, in-line frequency f1 = C1·√(EI/(m_e L⁴))',
    source: 'DNV-RP-F105 Eq. 6.16 with Table 6-1; Blevins, Formulas for Natural Frequency and Mode Shape, Table 8-1',
    inputs: { D: 400, t: 20, tCoat: 0, tConc: 0, rhoCont: 0, L: 40, gap: 1, boundary: 'pinned', Seff: 0, kc: 0 },
    expect: [{ kpi: 'f1IL', value: 0.522, tol: 0.005 }],
  }],
};
