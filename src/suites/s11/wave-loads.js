// First-order wave excitation on the floater: force and moment transfer functions, design-wave amplitudes,
// irregular-sea force statistics, mean wave drift, and a diffraction check on the governing column or beam.
import { maccamyFuchs } from '../../engine/hydro/morison.js';
import { meanDriftForce } from '../../engine/hydro/drift.js';
import { responseStatistics } from '../../engine/hydro/floater.js';
import { hullInputs, hullFrom, motionInputs, motionFrom, motionRangeAdvice, seaInputs, seaFrom, spectrum, omegaGrid, driftCoefficient, windCurrentInputs,
  G, DEG, LIMIT_HYDRO, LIMIT_CFD, need, round, rnd, fmt, adv, check, sev, steelCo2Text } from './_common.js';

const sinc = (u) => (Math.abs(u) < 1e-8 ? 1 : Math.sin(u) / u);

/** Froude–Krylov force amplitudes per unit wave amplitude on a box L × B × T in deep water (exact pressure integral). */
export function boxFroudeKrylov({ L, B, T, omega, headingDeg, rho }) {
  const k = omega * omega / G; const kx = k * Math.cos(headingDeg * DEG); const ky = k * Math.sin(headingDeg * DEG);
  const side = (1 - Math.exp(-k * T)) / k; const bot = Math.exp(-k * T);
  const ux = kx * L / 2; const uy = ky * B / 2;
  const mom = (kk, half) => (Math.abs(kk) < 1e-9 ? 0 : 2 * (Math.sin(kk * half) - kk * half * Math.cos(kk * half)) / (kk * kk));
  return [
    rho * G * side * B * sinc(uy) * 2 * Math.abs(Math.sin(ux)),
    rho * G * side * L * sinc(ux) * 2 * Math.abs(Math.sin(uy)),
    rho * G * bot * L * B * Math.abs(sinc(ux) * sinc(uy)),
    rho * G * bot * L * Math.abs(sinc(ux)) * Math.abs(mom(ky, B / 2)),
    rho * G * bot * B * Math.abs(sinc(uy)) * Math.abs(mom(kx, L / 2)),
    0,
  ];
}

const NAMES = ['Surge force', 'Sway force', 'Heave force', 'Roll moment', 'Pitch moment', 'Yaw moment'];

export default {
  id: 's11-wave-loads',
  suite: 11,
  title: 'Hydrodynamic wave loads on the hull',
  summary: 'Wave excitation forces and moments on the floater per unit wave amplitude, their design-wave and irregular-sea values, the mean drift force, and whether the slender-body assumption holds for the governing member.',
  covers: [2],
  fidelity: 'reduced-order',
  standards: ['DNV-RP-C205 (2019) §6.2, §7.2–7.4', 'DNV-RP-C103 (2015) §4 (column-stabilised units, design-wave loads)', 'API RP 2FPS (2nd ed., 2011) §7', 'ISO 19904-1 (2019) §8'],
  models: ['Froude–Krylov pressure integration (long-wave / exact box)', 'Morison inertia with strip added mass', 'Stochastically linearised Morison drag (Borgman)', 'MacCamy–Fuchs linear diffraction on a vertical cylinder (validity check)', 'Wave spectral response (JONSWAP)', 'Mean second-order drift: Maruo and short-wave cylinder limit'],
  limits: `${LIMIT_HYDRO} Box-shaped hulls are loaded by the undisturbed (Froude–Krylov) pressure only, moments from the bottom pressure. Mean drift is a screening estimate with a user reflection coefficient; quadratic transfer functions are not computed. ${LIMIT_CFD}`,
  searchTerms: 'wave excitation force floater Froude-Krylov Morison diffraction',
  inputs: [
    ...hullInputs({ type: 'semi', withMesh: false }),
    ...seaInputs({ Hs: 12, Tp: 15, ret: '100', heading: 0 }),
    { key: 'Hdes', label: 'Design regular wave height', unit: 'm', type: 'number', default: 22, min: 0.1, max: 40, group: 'Design wave', role: 'environment', bind: 'site.Hmax100', help: 'Height of the regular design wave (about 1.86·Hs for the 3-hour maximum).' },
    { key: 'Tdes', label: 'Design regular wave period', unit: 's', type: 'number', default: 14, min: 3, max: 30, group: 'Design wave', role: 'environment', help: 'Period of the design wave; scan it across the range allowed by the steepness limit to find the governing load.' },
    ...motionInputs().filter((i) => ['kRoll', 'kPitch', 'kYaw', 'zetaRad', 'zetaRoll', 'CdCol', 'CdPon', 'Kmoor', 'tendonEA', 'nFreq'].includes(i.key)),
    windCurrentInputs().find((i) => i.key === 'driftR'),
    { key: 'dOverLambdaMax', label: 'Slender-body limit D/λ', unit: '-', type: 'number', default: 0.2, min: 0.05, max: 0.5, group: 'Acceptance criteria', role: 'criteria', help: 'Above this ratio diffraction changes the wave load and a panel (radiation–diffraction) analysis is needed. DNV-RP-C205 uses 0.2.' },
  ],
  mesh: { key: 'nFreq', label: 'Frequency points', levels: [50, 100, 200], qoi: 'FzSig' },
  compute(inp) {
    const hull = hullFrom(inp);
    const mot = motionFrom(inp, hull);
    const sea = seaFrom(inp); const om = omegaGrid(inp.nFreq); const S = spectrum(om, sea);
    const rho = hull.rho; const beta = inp.headingDeg;
    const forceAmp = (omega, s) => {
      if (mot.kind === 'slender') { const r = mot.rao(omega, beta, s); return [0, 1, 2, 3, 4, 5].map((d) => Array.from(omega, (_, i) => Math.hypot(r.force.re[d][i], r.force.im[d][i]))); }
      const per = Array.from(omega, (w) => boxFroudeKrylov({ L: hull.Lc, B: hull.Bc, T: hull.T, omega: w, headingDeg: beta, rho }));
      return [0, 1, 2, 3, 4, 5].map((d) => per.map((p) => p[d]));
    };
    const F = forceAmp(om, sea);
    const wd = 2 * Math.PI / inp.Tdes; const a = inp.Hdes / 2;
    const lamDes = 2 * Math.PI * G / (wd * wd);
    need(inp.Hdes / lamDes < 0.143, `The design wave is steeper than the breaking limit H/λ = 1/7 (H/λ = ${fmt(inp.Hdes / lamDes)}): increase the period or reduce the height.`);
    const Fd = forceAmp([wd], null).map((v) => v[0] * a);
    const dur = inp.durationH * 3600;
    const st = F.map((f) => { const r = responseStatistics({ omega: om, amp: f, S, duration: dur }); return r.m0 > 0 ? r : { significant: 0, mpm: 0, Tz: 0 }; });
    const D = hull.type === 'semi' ? inp.Dcol : hull.type === 'spar' ? inp.Dspar : hull.type === 'tlp' ? inp.DcolT : hull.Bc;
    const lamP = 2 * Math.PI * G / (2 * Math.PI / inp.Tp) ** 2;
    const dl = D / lamP; const dlDes = D / lamDes;
    const mf = maccamyFuchs({ H: inp.Hdes, T: inp.Tdes, depth: inp.depth, D, rho });
    const KC = Math.PI * inp.Hdes / D;
    const dragInertia = (inp.CdCol ?? 1) * KC / (Math.PI * Math.PI * 2);
    const Dw = driftCoefficient(inp, hull, om, beta);
    const drift = meanDriftForce({ omega: om, S, D: Dw });
    const chk = [
      check('Slender-body validity at the peak period', dl, inp.dOverLambdaMax, '-', 'DNV-RP-C205 (2019) §6.1.2 and §7.2 (D/λ > 0.2: diffraction regime)',
        'Not a strength check: above the limit the member scatters the wave and the Froude–Krylov plus Morison load used here over-predicts the excitation; a radiation–diffraction analysis is then required.'),
      check('Slender-body validity at the design-wave period', dlDes, inp.dOverLambdaMax, '-', 'DNV-RP-C205 (2019) §6.1.2',
        'Same criterion for the regular design wave. Short design waves on wide columns or beams need diffraction theory.'),
    ];
    const gov = chk[0].uc >= chk[1].uc ? chk[0] : chk[1];
    const T = om.map((w) => round(2 * Math.PI / w, 4));
    const sc = (d) => (d < 3 ? 1e-6 : 1e-6);
    const cmErr = 100 * (2 - mf.CmEff) / mf.CmEff;
    return {
      kpis: [
        { key: 'FxDes', label: 'Surge force amplitude, design wave', value: round(Fd[0] / 1e6, 5), unit: 'MN', explain: 'Horizontal wave force amplitude along x in the regular design wave (inertia and Froude–Krylov; drag excluded). It loads the station-keeping system and the column–pontoon connections.' },
        { key: 'FyDes', label: 'Sway force amplitude, design wave', value: round(Fd[1] / 1e6, 5), unit: 'MN', explain: 'Horizontal wave force amplitude along y in the design wave; zero for head or following seas.' },
        { key: 'FzDes', label: 'Heave force amplitude, design wave', value: round(Fd[2] / 1e6, 5), unit: 'MN', explain: 'Vertical wave force amplitude. For a semi-submersible it passes through zero at the cancellation period, where column and pontoon forces balance.' },
        { key: 'MyDes', label: 'Pitch moment amplitude, design wave', value: round(Fd[4] / 1e6, 5), unit: 'MN·m', explain: 'Wave overturning moment about the transverse axis at the waterline origin.' },
        { key: 'MxDes', label: 'Roll moment amplitude, design wave', value: Math.abs(Fd[3]) < 1e-9 * Math.max(Math.abs(Fd[4]), Math.abs(Fd[5]), 1) ? 0 : round(Fd[3] / 1e6, 5), unit: 'MN·m', explain: 'Wave overturning moment about the longitudinal axis; needs a beam-sea component. In head or following seas it vanishes by symmetry: a value below 1e-9 of the pitch or yaw moment is arithmetic round-off and is reported as 0.' },
        { key: 'FxSig', label: 'Significant surge force amplitude', value: round(st[0].significant / 1e6, 5), unit: 'MN', explain: 'Twice the standard deviation of the surge wave force in the irregular sea state (with linearised drag).' },
        { key: 'FzSig', label: 'Significant heave force amplitude', value: round(st[2].significant / 1e6, 5), unit: 'MN', explain: 'Twice the standard deviation of the vertical wave force in the sea state.' },
        { key: 'FxMpm', label: 'Most probable maximum surge force', value: round(st[0].mpm / 1e6, 5), unit: 'MN', explain: 'Expected largest surge force amplitude in the storm duration, σ·√(2 ln N).' },
        { key: 'MyMpm', label: 'Most probable maximum pitch moment', value: round(st[4].mpm / 1e6, 5), unit: 'MN·m', explain: 'Expected largest wave pitching moment in the storm duration.' },
        { key: 'drift', label: 'Mean wave-drift force', value: round(drift / 1e3, 5), unit: 'kN', explain: 'Steady second-order force in the wave direction, 2∫S·D dω with a screening drift coefficient. It adds to wind and current in the mean offset.' },
        { key: 'dOverLambda', label: 'Member size / wavelength at Tp', value: round(dl, 4), unit: '-', limit: inp.dOverLambdaMax, sense: 'max', explain: 'Governing diameter (or beam) over the peak wavelength. Below 0.2 the member is hydrodynamically transparent and the slender-body loads are valid.' },
        { key: 'CmEff', label: 'MacCamy–Fuchs effective inertia coefficient', value: round(mf.CmEff, 4), unit: '-', explain: 'Inertia coefficient that reproduces the linear diffraction force on a vertical cylinder of the governing diameter at the design period. Compare with the slender value of 2.0: the further below 2, the more diffraction matters.' },
        { key: 'KC', label: 'Keulegan–Carpenter number', value: round(KC, 4), unit: '-', explain: 'π·H/D at the surface in the design wave. Below about 2 the load is inertia-dominated and drag is a small correction.' },
        { key: 'dragShare', label: 'Drag / inertia force ratio', value: round(dragInertia, 4), unit: '-', explain: 'Peak Morison drag over peak inertia force per unit length at the surface, Cd·KC/(π²·Cm). Small values justify linear frequency-domain analysis.' },
      ],
      checks: chk,
      plots: [
        { type: 'line', title: 'Force transfer functions', x: { label: 'Wave period', unit: 's' }, y: { label: 'Force per metre amplitude', unit: 'MN/m' },
          series: [0, 1, 2].map((d) => ({ name: NAMES[d], x: T, y: rnd(F[d].map((v) => v * sc(d))) })).filter((s) => s.y.some((v) => v > 1e-9)),
          explain: 'Wave force amplitude per metre of wave amplitude against period. Look for the heave cancellation (a zero near 20 s on a semi) and for the short-period range where the slender-body model stops being valid.' },
        { type: 'line', title: 'Moment transfer functions', x: { label: 'Wave period', unit: 's' }, y: { label: 'Moment per metre amplitude', unit: 'MN·m/m' },
          series: [3, 4, 5].map((d) => ({ name: NAMES[d], x: T, y: rnd(F[d].map((v) => v * 1e-6)) })).filter((s) => s.y.some((v) => v > 1e-9)),
          explain: 'Wave moment amplitude per metre of wave amplitude. Peaks occur where half a wavelength matches the column or hull spacing.' },
        { type: 'line', title: 'Wave spectrum and surge-force spectrum', x: { label: 'Frequency', unit: 'rad/s' }, y: { label: 'Normalised spectral density', unit: '-' },
          series: [{ name: 'Wave spectrum', x: rnd(om), y: rnd(Array.from(S, (v) => v / Math.max(...S, 1e-12))) }, (() => { const sf = F[0].map((f, i) => f * f * S[i]); const m = Math.max(...sf, 1e-12); return { name: 'Surge force spectrum', x: rnd(om), y: rnd(sf.map((v) => v / m)), dash: true }; })()],
          explain: 'Both spectra scaled to a peak of one: where the force spectrum sits relative to the wave spectrum shows which wave periods load the hull most.' },
        { type: 'line', title: 'Mean drift coefficient', x: { label: 'Frequency', unit: 'rad/s' }, y: { label: 'F / a²', unit: 'kN/m²' }, series: [{ name: 'Drift coefficient', x: rnd(om), y: rnd(Dw.map((v) => v / 1e3)) }],
          explain: 'Screening mean drift force per squared wave amplitude. It rises with frequency as shorter waves are reflected; the mean drift in the sea state is the spectrum-weighted integral.' },
      ],
      tables: [
        { title: 'Wave loads summary', columns: [{ key: 'c', label: 'Component' }, { key: 'des', label: 'Design-wave amplitude' }, { key: 'sig', label: 'Significant amplitude' }, { key: 'mpm', label: 'Most probable maximum' }, { key: 'tz', label: 'Tz', unit: 's' }, { key: 'u', label: 'Unit' }],
          rows: NAMES.map((n, d) => ({ c: n, des: round(Fd[d] / 1e6, 5), sig: round(st[d].significant / 1e6, 5), mpm: round(st[d].mpm / 1e6, 5), tz: round(st[d].Tz, 4), u: d < 3 ? 'MN' : 'MN·m' })),
          explain: 'One row per load component about the waterline origin: regular design-wave amplitude and the irregular-sea statistics for the same heading.' },
      ],
      notes: [
        `${mot.kind === 'slender' ? 'Slender-member model' : 'Box Froude–Krylov model'}; governing dimension ${fmt(D)} m, D/λ = ${fmt(dl)} at Tp and ${fmt(dlDes)} at the design period.`,
        `MacCamy–Fuchs effective Cm = ${fmt(mf.CmEff)} (ka = ${fmt(mf.ka)}): the slender value of 2.0 ${Math.abs(cmErr) < 5 ? 'is adequate' : `overstates the column inertia force by about ${fmt(cmErr)} %`}.`,
        mot.kind === 'slender' ? `Loads are valid for periods above ${fmt(2 * Math.PI / mot.validOmegaMax)} s.` : 'Diffraction and radiation forces on the box are not included: values are the undisturbed-wave pressure integrals.',
      ],
      advice: [
        ...motionRangeAdvice(mot, Math.min(inp.Tp, inp.Tdes), inp.headingDeg),
        adv('engineering', gov.uc > 1 ? 'warn' : 'ok', `Design-wave loads: surge ${fmt(Fd[0] / 1e6)} MN, heave ${fmt(Fd[2] / 1e6)} MN, pitch moment ${fmt(Fd[4] / 1e6)} MN·m; D/λ = ${fmt(dl)} against the slender-body limit ${inp.dOverLambdaMax}.`,
          gov.uc > 1 ? 'The hull is in the diffraction regime for this sea: use these loads for sizing only and obtain excitation from a panel analysis, imported through the RAO module.'
            : `Scan the design-wave period between ${fmt(Math.sqrt(inp.Hdes * 7 * 2 * Math.PI / G))} s (steepness limit) and ${fmt(inp.Tp * 1.3)} s and the heading in 15° steps: split forces and racking moments on column-stabilised units peak at specific period–heading pairs.`,
          'DNV-RP-C103 (2015) §4 (characteristic global responses, design-wave approach); DNV-RP-C205 (2019) §6–7.'),
        adv('economic', 'info', `The mean drift of ${fmt(drift / 1e3)} kN and the surge force maximum of ${fmt(st[0].mpm / 1e6)} MN set the station-keeping capacity; wave loads scale roughly with displaced volume near the surface.`,
          'Before adding steel to the connections, reduce the load at source: deeper pontoons and smaller waterplane cut the wave-frequency force, which reduces mooring size, fatigue demand and inspection scope together.',
          'Design optimisation practice of DNV-RP-C103 §2; whole-life costing per ISO 15686-5 (operating cost of station keeping and inspection).',
          `Load-driven steel is carbon: every 100 t of connection reinforcement avoided saves ${fmt(190)} t CO₂e at the screening factor. ${steelCo2Text(100)}`),
      ],
    };
  },
  verify: [
    { // Spar D = 30 m, T = 90 m, ω = 0.4 rad/s (T = 15.708 s), unit amplitude (H = 2 m): k = ω²/g = 0.016315, e^(−kT) = 0.23030.
      // Heave: (ρ g Aw − (2/3) π ρ r³ ω²)·e^(−kT) = (7.1052e6 − 7.2453e6·0.16)·0.23030 = (7.1052e6 − 1.1592e6)·0.23030 = 1.3693 MN.
      // Surge: (1 + Ca)·ρ·A·g·(1 − e^(−kT)) = 2·1025·706.858·9.80665·0.76970 = 10.938 MN, less the added mass lost around the
      // free keel end at e^(−kT): the column and its image are 180 m long, l/D = 6, C_A = 0.90 + 0.06·(6 − 5)/(9 − 5) = 0.915
      // (DNV-RP-C205 Table D-2), loss = 0.085·ρ·A·T = 0.085·65.208e6 = 5.5427e6 kg, ω²·loss·e^(−kT) = 0.16·5.5427e6·0.23030
      // = 0.2042 MN → 10.734 MN.
      name: 'Spar: heave and surge excitation of a deep vertical cylinder in long waves',
      source: 'Faltinsen (1990) Sea Loads on Ships and Offshore Structures ch. 3 — long-wave force on a vertical cylinder, hand calculation',
      inputs: { hullType: 'spar', Dspar: 30, Hspar: 105, draft: 90, KG: 38, Hdes: 2, Tdes: 15.70796, headingDeg: 0, depth: 4000 },
      expect: [{ kpi: 'FzDes', value: 1.3693, tol: 0.005 }, { kpi: 'FxDes', value: 10.734, tol: 0.01 }] },
    { // Box 100 × 20 × 5 m draft, beam... head sea, ω = 0.5: k = 0.025493, kL/2 = 1.27466, sinc = 0.75034, e^(−kT) = 0.88033.
      // Fz = ρ g L B e^(−kT) sinc(kL/2) = 1025·9.80665·2000·0.88033·0.75034 = 13.279 MN per metre amplitude.
      name: 'Box barge: Froude–Krylov heave force in a head sea',
      source: 'Exact integral of the undisturbed deep-water wave pressure over the bottom of a rectangular box (hand calculation)',
      inputs: { hullType: 'barge', L: 100, B: 20, Dhull: 12, draft: 5, KG: 4, Hdes: 2, Tdes: 12.56637, headingDeg: 180 },
      expect: [{ kpi: 'FzDes', value: 13.279, tol: 0.003 }] },
  ],
};
