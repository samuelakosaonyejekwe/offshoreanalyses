// Cable free span: static sag and curvature of a tensioned beam, natural frequencies, vortex-induced vibration
// onset and amplitude (response models in the DNV-RP-F105 form), allowable span length and strain-based fatigue
// of the metallic sheath from the vibration.
import { tensionedBeamModes } from '../../engine/lines/viv.js';
import { stabilityParameter, inlineResponse, crossflowResponse, vivScreening, addedMassCoefficient } from '../../engine/pipeline/span.js';
import { nearBedSpectralVelocity } from '../../engine/pipeline/stability.js';
import { lineInputs, lineFrom, adv, check, worst, sev, fmt, round, run, need, linspace, CONFIRM, STD_CABLE, CO2_ROCK, vesselCo2Text } from './_common.js';

/** First natural frequency of a pinned tensioned beam: f1 = (1/2L)·sqrt(T/m + π²EI/(m·L²)). */
const f1Of = (L, T, EI, me) => Math.sqrt(T / me + Math.PI ** 2 * EI / (me * L * L)) / (2 * L);

export default {
  id: 's13-free-span-viv', suite: 13, title: 'Cable free span and vibration',
  summary: 'Sag, curvature, natural frequencies and vortex-induced vibration susceptibility of an unsupported cable length, with the allowable span and the fatigue life of the metallic sheath.',
  covers: [12], fidelity: 'analytical',
  standards: ['DNV-RP-F105 (2017) Free spanning pipelines (response models, screening, safety factors)', 'DNV-RP-0360 (2016) §4.7 (free spans of subsea power cables)', 'DNV-RP-C205 (2021) §9 (vortex-induced oscillations)', 'Blevins, Formulas for Natural Frequency and Mode Shape (1979) Table 8'],
  models: ['Euler–Bernoulli tensioned beam (static sag, curvature and natural frequencies)', 'Vortex shedding: reduced velocity and stability parameter', 'In-line and cross-flow VIV response models', 'Wave-induced near-bed velocity (JONSWAP transfer)', 'Strain-based fatigue of the metallic sheath with Miner summation'],
  limits: `Single span with pinned ends (lower-bound frequency) and uniform tension; the DNV-RP-F105 response models were calibrated for steel pipelines and are applied to the cable as a screening method with the cable's own mass, damping and bending stiffness. Modes above the fifth, multi-span interaction, span shoulders settling, and stick–slip bending hysteresis (which adds damping) are not modelled. The ε–N parameters are indicative.`,
  publishes: ['cable.spanAllowable', 'cable.spanF1'],
  searchTerms: 'submarine cable free span vortex induced vibration fatigue',
  inputs: [
    ...lineInputs({ hydro: false }),
    { key: 'Lspan', label: 'Span length', unit: 'm', type: 'number', default: 3, min: 0.5, max: 200, group: 'Span', role: 'geometry', design: true, help: 'Unsupported length from survey (scour, sandwave trough, crossing shoulder, rock outcrop).' },
    { key: 'gap', label: 'Gap between cable and seabed', unit: 'm', type: 'number', default: 0.3, min: 0.01, max: 20, group: 'Span', role: 'geometry' },
    { key: 'T', label: 'Residual effective tension in the span', unit: 'kN', type: 'number', default: 8, min: 0.01, max: 2000, group: 'Span', role: 'ic', bind: 'cable.tdpTension', design: true, uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Initial pretension left by laying; the span shoulders are taken as pinned supports.' },
    { key: 'Ca', label: 'Added-mass coefficient', unit: '-', type: 'number', default: 1.0, min: 0.5, max: 3, group: 'Span', role: 'load', help: 'About 1.0 far from the seabed, rising as the gap closes (DNV-RP-F105 Eq. 6.21 value is reported).' },
    { key: 'zeta', label: 'Total modal damping ratio', unit: '-', type: 'number', default: 0.03, min: 0.002, max: 0.3, group: 'Span', role: 'material', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Structural (armour stick–slip gives a few per cent in cables), soil and hydrodynamic damping together.' },
    { key: 'Uc', label: 'Current velocity at the span', unit: 'm/s', type: 'number', default: 0.35, min: 0, max: 4, group: 'Environment', role: 'bc', bind: 'site.UcBed100', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Velocity normal to the cable at its level (100-year for screening).' },
    { key: 'Hs', label: 'Significant wave height', unit: 'm', type: 'number', default: 1.5, min: 0, max: 20, group: 'Environment', role: 'bc', bind: 'site.Hs1' },
    { key: 'Tp', label: 'Peak period', unit: 's', type: 'number', default: 7, min: 3, max: 25, group: 'Environment', role: 'bc', bind: 'site.Tp1' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 40, min: 3, max: 3000, group: 'Environment', role: 'environment', bind: 'site.depth' },
    { key: 'occ', label: 'Fraction of time with vibration-inducing flow', unit: '-', type: 'number', default: 0.2, min: 0, max: 1, group: 'Environment', role: 'load', help: 'Share of the year the flow at the span exceeds the onset velocity (from the current statistics).' },
    { key: 'rSheath', label: 'Radius of the metallic sheath of one core', unit: 'mm', type: 'number', default: 27.7, min: 3, max: 120, group: 'Fatigue resistance', role: 'geometry', help: 'Each core bends about its own axis: strain = curvature × this radius.' },
    { key: 'epsRef', label: 'Sheath strain range at 10⁶ cycles', unit: '%', type: 'number', default: 0.1, min: 0.005, max: 2, group: 'Fatigue resistance', role: 'material', help: `ε–N power law N = 10⁶·(ε_ref/Δε)^m; indicative for lead alloy. ${CONFIRM}` },
    { key: 'epsM', label: 'ε–N slope m', unit: '-', type: 'number', default: 3.5, min: 1.5, max: 10, group: 'Fatigue resistance', role: 'material' },
    { key: 'gammaF', label: 'Safety factor on natural frequency', unit: '-', type: 'number', default: 1.4, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: `Screening factor of DNV-RP-F105 §2.3 (1.4 in-line and cross-flow). ${CONFIRM}` },
    { key: 'mbrStatic', label: 'Allowable static bend radius', unit: 'm', type: 'number', default: 1.9, min: 0.2, max: 60, group: 'Acceptance criteria', role: 'criteria', bind: 'cable.mbrStatic' },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 10, min: 1, max: 20, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'life', label: 'Design life', unit: 'years', type: 'number', default: 25, min: 1, max: 60, group: 'Acceptance criteria', role: 'criteria', bind: 'life.designYears' },
  ],
  compute(inp) {
    const line = lineFrom({ ...inp, Ca: inp.Ca }), D = line.D, w = line.w, EI = line.EI, T = inp.T * 1e3, L = inp.Lspan;
    need(w > 0, `The cable is buoyant (apparent weight ${fmt(w)} N/m).`);
    const me = line.m + line.ma;
    const modes = run(() => tensionedBeamModes({ L, EI, m: me, T, nModes: 5, nPoints: 41 }), 'Natural frequencies');
    const f1 = modes.modes[0].fHz;
    // static tensioned beam, pinned ends, uniform load w: u = (L/2)·sqrt(T/EI)
    const u = L / 2 * Math.sqrt(T / EI), sech = 1 / Math.cosh(Math.min(u, 300));
    const kStat = w / T * (1 - sech), sag = w * L * L / (8 * T) - w * EI / (T * T) * (1 - sech);
    const sagX = (x) => { const k = Math.sqrt(T / EI), c = Math.min(u, 300); return w / (2 * T) * x * (L - x) - w * EI / (T * T) * (1 - Math.cosh(k * (x - L / 2)) / Math.cosh(c)); };
    // environment
    const nb = inp.Hs > 0 ? nearBedSpectralVelocity({ Hs: inp.Hs, Tp: inp.Tp, depth: inp.depth }) : { Us: 0, Tu: inp.Tp };
    const Uw = nb.Us, U = inp.Uc + Uw, alpha = U > 0 ? inp.Uc / U : 1;
    const KC = Uw > 0 ? Uw * nb.Tu / D : Infinity;
    const st = stabilityParameter({ me, zetaT: inp.zeta, D, rhoW: line.rhoW });
    const massRatio = line.m / (line.rhoW * line.Ae), gapRatio = inp.gap / D;
    const ilAt = (VR) => inlineResponse({ VR, Ksd: st.Ksd, alpha });
    const cfAt = (VR) => crossflowResponse({ VR, Ksd: st.Ksd, alpha, KC, gapRatio, massRatio });
    const VR = U / (f1 * D), il = ilAt(VR), cf = cfAt(VR);
    const scr = vivScreening({ f1IL: f1, f1CF: f1, Uc100: inp.Uc, Uw1: Uw, D, L, VRonsetIL: il.VRonset, VRonsetCF: cf.VRonset, gammaIL: inp.gammaF, gammaCF: inp.gammaF });
    // every mode n can lock in: curvature amplitude of a sine mode κ = A·(nπ/L)², strain range = 2·κ·r, cycles at f_n
    const rS = inp.rSheath / 1e3, cyc = 31557600 * inp.occ;
    const nOf = (eps) => (eps > 0 ? 1e6 * (inp.epsRef / 100 / eps) ** inp.epsM : Infinity);
    const perMode = modes.modes.map((m) => {
      const vr = U / (m.fHz * D), a = ilAt(vr).AyOverD, c = cfAt(vr).AzOverD, kk = (m.n * Math.PI / L) ** 2;
      const eIL = 2 * a * D * kk * rS, eCF = 2 * c * D * kk * rS;
      return { n: m.n, f: m.fHz, vr, a, c, eps: Math.max(eIL, eCF), d: cyc * m.fHz * (1 / nOf(eIL) + 1 / nOf(eCF)) };
    });
    const dYr = perMode.reduce((s2, q) => s2 + q.d, 0), epsMax = Math.max(...perMode.map((q) => q.eps));
    const ampIL = Math.max(...perMode.map((q) => q.a)), ampCF = Math.max(...perMode.map((q) => q.c));
    const lifeYr = dYr > 0 ? Math.min(1 / dYr, 1e6) : 1e6;
    // allowable span lengths: frequency criterion f1 ≥ γ·U/(VRonset·D), and static curvature
    const solveL = (fReq) => { if (!(fReq > 0)) return 500; let lo = 0.2, hi = 500; if (f1Of(hi, T, EI, me) > fReq) return hi; for (let i = 0; i < 60; i++) { const m = 0.5 * (lo + hi); if (f1Of(m, T, EI, me) > fReq) lo = m; else hi = m; } return lo; };
    const LallowIL = solveL(inp.gammaF * U / (il.VRonset * D)), LallowCF = solveL(inp.gammaF * U / (cf.VRonset * D));
    const checks = [
      check('In-line VIV onset (screening)', scr.ucIL, 1, '-', 'DNV-RP-F105 (2017) §2.3 (in-line screening criterion)', 'Required over available natural frequency for no in-line vortex-induced vibration. Above 1.0 the span vibrates in line with the flow at small amplitude but high cycle count.'),
      check('Cross-flow VIV onset (screening)', scr.ucCF, 1, '-', 'DNV-RP-F105 (2017) §2.3 (cross-flow screening criterion)', 'Required over available natural frequency for no cross-flow vibration. Cross-flow amplitudes approach one diameter and fatigue the sheath quickly.'),
      check('Static curvature in the span', kStat, 1 / inp.mbrStatic, '1/m', `${STD_CABLE} (static MBR); tensioned-beam statics`, 'Mid-span curvature under the submerged weight against the allowable static curvature. A long span at low tension bends the cable beyond its installed bend radius.'),
      check('Sag against the seabed gap', sag, inp.gap, 'm', 'Tensioned-beam deflection', 'Mid-span sag against the available gap. When the sag exceeds the gap the cable touches down and the real span is shorter — re-run with the reduced length.'),
      check('Sheath fatigue from vibration (Miner × DFF over life)', dYr * inp.dff * inp.life, 1, '-', 'DNV-RP-F105 §2.4 (fatigue criterion) with ε–N of the sheath; DNV-RP-C203 §2.2', 'Accumulated strain-fatigue damage of the metallic sheath from the predicted vibration amplitudes. Above 1.0 the water barrier is expected to crack within the design life.'),
    ];
    const gov = worst(checks);
    const Ls = linspace(Math.max(1, 0.25 * L), 3 * L, 40), VRs = linspace(0.5, 12, 60);
    const xs = linspace(0, L, 41);
    return {
      kpis: [
        { key: 'f1', label: 'First natural frequency', value: round(f1, 5), unit: 'Hz', explain: 'Lowest natural frequency of the span as a pinned tensioned beam in water (with added mass). Lower frequency means earlier vortex lock-in.' },
        { key: 'f2', label: 'Second natural frequency', value: round(modes.modes[1].fHz), unit: 'Hz', explain: 'Second mode; it can be excited at higher velocities.' },
        { key: 'VR', label: 'Reduced velocity', value: round(VR), unit: '-', explain: 'Flow velocity over (frequency × diameter). In-line vibration starts near 1, cross-flow near 2–3, with the largest cross-flow amplitude around 5–7.' },
        { key: 'Ksd', label: 'Design stability parameter', value: round(st.Ksd), unit: '-', explain: 'Mass-damping parameter 4π·mₑ·ζ/(ρD²) with its safety factor. High values (heavy, well-damped cables) suppress in-line vibration.' },
        { key: 'onsetIL', label: 'In-line onset reduced velocity', value: round(il.VRonset), unit: '-', explain: 'Reduced velocity at which in-line vibration starts for this stability parameter.' },
        { key: 'onsetCF', label: 'Cross-flow onset reduced velocity', value: round(cf.VRonset), unit: '-', explain: 'Reduced velocity at which cross-flow vibration starts, including seabed proximity and mass-ratio effects.' },
        { key: 'ampIL', label: 'In-line amplitude A/D', value: round(ampIL), unit: '-', explain: 'Largest predicted in-line vibration amplitude over diameter among the first five modes at the design velocity (0 = no vibration).' },
        { key: 'ampCF', label: 'Cross-flow amplitude A/D', value: round(ampCF), unit: '-', explain: 'Largest predicted cross-flow amplitude over diameter among the first five modes at the design velocity.' },
        { key: 'sag', label: 'Static sag at mid-span', value: round(sag, 5), unit: 'm', limit: inp.gap, sense: 'max', explain: 'Deflection of the span under its submerged weight with the residual tension; it must stay below the gap for the span to be free.' },
        { key: 'curvStatic', label: 'Static curvature at mid-span', value: round(kStat, 5), unit: '1/m', limit: round(1 / inp.mbrStatic), sense: 'max', explain: 'Bending of the cable in the span, (w/T)(1 − sech u); it tends to w/T for long spans.' },
        { key: 'strainRange', label: 'Sheath strain range from vibration', value: round(epsMax * 100, 4), unit: '%', explain: 'Largest cyclic strain of the metallic sheath from the modal curvature of any excited mode.' },
        { key: 'fatigueLife', label: 'Sheath fatigue life (no factor)', value: round(lifeYr), unit: 'years', explain: 'Characteristic life under the predicted vibration for the stated occurrence; 10⁶ means no vibration is predicted.' },
        { key: 'LallowIL', label: 'Allowable span, in-line criterion', value: round(LallowIL), unit: 'm', explain: 'Longest span for which the natural frequency stays above the in-line onset frequency with the safety factor (500 = not limiting).' },
        { key: 'LallowCF', label: 'Allowable span, cross-flow criterion', value: round(LallowCF), unit: 'm', explain: 'Longest span that avoids cross-flow vibration; spans between the two limits vibrate in line only and need the fatigue check.' },
        { key: 'CaF105', label: 'Added-mass coefficient near the seabed (reference)', value: round(addedMassCoefficient(gapRatio)), unit: '-', explain: 'DNV-RP-F105 value for this gap ratio, for comparison with the entered coefficient.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Natural frequency against span length', x: { label: 'Span length', unit: 'm' }, y: { label: 'Frequency', unit: 'Hz' },
          series: [{ name: 'First natural frequency', x: Ls, y: Ls.map((l) => f1Of(l, T, EI, me)) }, { name: 'Required, in-line', x: [Ls[0], Ls[39]], y: [inp.gammaF * U / (il.VRonset * D), inp.gammaF * U / (il.VRonset * D)], dash: true }, { name: 'Required, cross-flow', x: [Ls[0], Ls[39]], y: [inp.gammaF * U / (cf.VRonset * D), inp.gammaF * U / (cf.VRonset * D)], dash: true }, { name: 'This span', x: [L], y: [f1] }],
          explain: 'The frequency falls with span length. The crossings with the dashed lines are the allowable span lengths; the marker shows where this span sits.' },
        { type: 'line', title: 'VIV response model', x: { label: 'Reduced velocity', unit: '-' }, y: { label: 'Amplitude / diameter', unit: '-' },
          series: [{ name: 'In-line', x: VRs, y: VRs.map((v) => ilAt(v).AyOverD) }, { name: 'Cross-flow', x: VRs, y: VRs.map((v) => cfAt(v).AzOverD) }, { name: 'Design point', x: [VR, VR], y: [0, Math.max(ampCF, ampIL, 0.1)], dash: true }],
          explain: 'Amplitude envelopes of the response models for this stability parameter. The dashed line is the design reduced velocity: left of both curves means no vibration.' },
        { type: 'line', title: 'Static deflected shape', x: { label: 'Distance along the span', unit: 'm' }, y: { label: 'Elevation', unit: 'm' },
          series: [{ name: 'Cable', x: xs, y: xs.map((x) => -sagX(x)) }, { name: 'Seabed', x: [0, L], y: [-inp.gap, -inp.gap], dash: true }],
          explain: 'Sag of the span under its weight between pinned shoulders. If the curve reaches the seabed line the cable touches down and the span divides.' },
      ],
      tables: [{ title: 'Span modes', explain: 'Natural frequencies of the first five modes with the reduced velocity each sees at the design flow, the predicted amplitudes and the sheath damage rate each contributes.',
        columns: [{ key: 'n', label: 'Mode' }, { key: 'f', label: 'Frequency', unit: 'Hz' }, { key: 'VR', label: 'Reduced velocity', unit: '-' }, { key: 'a', label: 'In-line A/D', unit: '-' }, { key: 'c', label: 'Cross-flow A/D', unit: '-' }, { key: 'd', label: 'Sheath damage', unit: '1/yr' }], rows: perMode.map((q) => ({ n: q.n, f: round(q.f), VR: round(q.vr), a: round(q.a), c: round(q.c), d: round(q.d, 3) })) }],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        `Flow at the span: current ${inp.Uc} m/s + significant wave velocity ${fmt(Uw)} m/s (current ratio ${fmt(alpha)}); reduced velocity ${fmt(VR)}.`,
        `Tension parameter u = (L/2)·√(T/EI) = ${fmt(u)}: the span behaves ${u > 3 ? 'as a tensioned string' : 'as a beam'}.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name}: ${fmt(gov.uc)}. f₁ = ${fmt(f1)} Hz, reduced velocity ${fmt(VR)}, amplitudes ${fmt(ampIL)} D in-line and ${fmt(ampCF)} D cross-flow; allowable span ${fmt(Math.min(LallowIL, LallowCF))} m (cross-flow ${fmt(LallowCF)} m).`,
          gov.uc > 1 ? `Shorten the span to below ${fmt(gov.name.includes('Cross') ? LallowCF : gov.name.includes('In-line') ? LallowIL : L / Math.sqrt(gov.uc))} m with grout bags or rock supports at mid-span, or re-bury the section; raising the residual tension helps only marginally and loads the adjacent joints.`
            : 'Record the span in the survey register, re-survey after the first storms (spans grow by scour), and re-run with the measured length and gap; cables are more tolerant than pipes because of their damping, but lead sheaths are not.',
          'DNV-RP-F105 (2017) §2.3–2.4; DNV-RP-0360 (2016) §4.7'),
        adv('economic', 'info', `Rectifying one ${fmt(L)} m span needs of the order of ${fmt(Math.max(1, Math.ceil(L / Math.max(Math.min(LallowIL, LallowCF), 1)) - 1))} support(s) or ${fmt(2.5 * L)} t of rock and about one vessel day; leaving it requires the fatigue life of ${fmt(lifeYr)} years to exceed ${inp.dff} × ${inp.life} years.`,
          gov.uc <= 1 ? 'No intervention is needed now: put the span on a risk-based survey interval instead of rectifying it, and intervene only if it grows beyond the allowable length.' : 'Combine span rectification with other route interventions in one campaign: mobilisation dominates the cost of a single rock dump or grout-bag operation.',
          'Risk-based inspection planning (DNV-RP-G101 principles); whole-life costing ISO 15686-5',
          `Rock of ${fmt(2.5 * L)} t ≈ ${fmt(2.5 * L * CO2_ROCK)} t CO₂e at the quarry (screening ${CO2_ROCK} t CO₂e/t); the vessel dominates: ${vesselCo2Text(1)}`),
      ],
      publish: { 'cable.spanAllowable': Math.min(LallowIL, LallowCF), 'cable.spanF1': f1 },
      quantities: {},
    };
  },
  verify: [
    // Pinned tensioned beam (Blevins Table 8-1 with axial tension): f1 = (1/2L)·sqrt(T/mₑ + π²EI/(mₑL²)).
    // mₑ = 63 + 1.0 × 1025 × π × 0.16²/4 = 83.609 kg/m; L = 20 m, T = 5 kN, EI = 12 000 N·m²:
    // T/mₑ = 59.80, π²EI/(mₑL²) = 3.541 → f1 = sqrt(63.34)/40 = 0.19897 Hz.
    // Static: w = 415.73 N/m, u = 10·sqrt(5000/12000) = 6.455, sech u = 0.003146 → κ = (w/T)(1 − sech u) = 0.08288 1/m;
    // sag = wL²/(8T) − (w·EI/T²)(1 − sech u) = 4.1573 − 0.1989 = 3.9584 m.
    { name: 'Tensioned pinned beam: first natural frequency, mid-span curvature and sag',
      source: 'Blevins (1979) Table 8-1 (tensioned beam); Timoshenko, Strength of Materials II (tie rod with lateral load); hand calculation',
      inputs: { Lspan: 20, T: 5, gap: 10, Hs: 0 },
      expect: [{ kpi: 'f1', value: 0.19897, tol: 0.003 }, { kpi: 'curvStatic', value: 0.08288, tol: 0.003 }, { kpi: 'sag', value: 3.9584, tol: 0.003 }] },
  ],
};
