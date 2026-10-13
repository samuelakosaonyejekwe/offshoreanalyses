// Vessel-motion-induced riser response: hang-off motion from the host RAOs drives the solved catenary riser in
// the frequency domain; dynamic tension, bending stress and amplification relative to the quasi-static response.
import { riserOffsetCases, quasiStaticHeave } from '../../engine/lines/riser.js';
import {
  DEG, need, run, fmt, round, thin, scaled, adv, check, kpi, plot, meshInput, finish, handWeight, steelTonnes, steelCo2Text, LIMIT_2D, CONFIRM,
} from './_common.js';
import { scrInputs, criteriaInputs, buildScr, stressOf, tdpNode } from './_scr.js';
import { vesselInputs, hangoffRao, vesselRangeAdvice, LIMIT_VESSEL } from './_vessel.js';
import { dynInputs, harmonic, stressAmp, sweep, seaGrid, responseStats, mpm, mvOf, tdpTravel } from './_dyn.js';

export default {
  id: 's09-vessel-motion', suite: 9, title: 'Vessel-motion-induced response',
  summary: 'Dynamic top tension, touchdown bending stress and dynamic amplification of a catenary riser driven by host-vessel heave, pitch and surge in a design wave and in an irregular sea.',
  covers: [4], fidelity: 'reduced-order',
  standards: ['API RP 2RD (1998, R2013) §5.5 dynamic analysis', 'DNV-ST-F201 (2018) Sec. 4, App. A', 'DNV-RP-C205 (2019) §3.5, §6.2, §7'],
  models: ['Response amplitude operators (closed-form box vessel or user table)', 'Frequency-domain rod dynamics with consistent mass and added mass', 'Harmonically linearised Morison drag', 'JONSWAP spectrum', 'Rayleigh extreme-value statistics', 'Quasi-static offset method'],
  limits: `${LIMIT_2D} Dynamics are linearised about the static shape: drag damping is equivalent-linear, seabed contact status is frozen, and wave loads act only through vessel motion (no direct wave kinematics on the riser). The vessel transfer functions are closed-form estimates, not a diffraction analysis; six-degree-of-freedom motions reduce to in-plane surge and vertical motion at the porch. A nonlinear time-domain run is needed where compression or snap loading is predicted.` + LIMIT_VESSEL,
  publishes: ['riser.dynTensionAmp', 'riser.daf', 'riser.hangoffHeaveAmp'],
  searchTerms: 'steel catenary riser vessel motion heave dynamic response touchdown',
  inputs: [
    ...scrInputs(),
    ...vesselInputs(),
    { key: 'Hmax', label: 'Design regular wave height', unit: 'm', type: 'number', default: 18, min: 0, max: 40, group: 'Waves', role: 'bc', bind: 'site.Hmax100', uncertain: { dist: 'normal', cov: 0.1 },
      help: 'Height of the regular design wave (e.g. 100-year Hmax). Set 0 to obtain the static result.' },
    { key: 'Twave', label: 'Design wave period', unit: 's', type: 'number', default: 13, min: 3, max: 30, group: 'Waves', role: 'bc', bind: 'site.Tp100' },
    { key: 'Hs', label: 'Significant wave height (irregular sea)', unit: 'm', type: 'number', default: 10, min: 0.1, max: 22, group: 'Waves', role: 'bc', bind: 'site.Hs100', uncertain: { dist: 'normal', cov: 0.1 } },
    { key: 'Tp', label: 'Spectral peak period', unit: 's', type: 'number', default: 14, min: 3, max: 30, group: 'Waves', role: 'bc', bind: 'site.Tp100' },
    { key: 'gamma', label: 'JONSWAP peak factor', unit: '-', type: 'number', default: 2.5, min: 1, max: 7, group: 'Waves', role: 'bc' },
    ...dynInputs(),
    { key: 'scf', label: 'Stress concentration factor', unit: '-', type: 'number', default: 1.2, min: 1, max: 5, group: 'Dynamic model', role: 'material', help: 'Applied to the dynamic stress amplitude reported for the welds.' },
    ...criteriaInputs().map((c) => (c.key === 'Cf' ? { ...c, default: 1.2 } : c)),
    meshInput(160),
  ],
  mesh: { key: 'nElem', label: 'Beam elements', levels: [80, 160, 320, 640], qoi: 'TtopMax' },

  compute(inp, ctx) {
    // engineering units → SI in buildScr; wave heights m, periods s
    const { pp, base, pressure } = buildScr(inp, ctx, { nElem: inp.nElem });
    const n = base.s.length - 1, kT = tdpNode(base), Ts = base.topTension;
    const omega = 2 * Math.PI / inp.Twave, a = inp.Hmax / 2;
    const rao = hangoffRao(inp, [omega]);
    const xAmp = rao.surge[0] * a, zAmp = rao.heave[0] * a;
    need(Math.hypot(xAmp, zAmp) < 0.1 * (inp.depth - inp.hangoffDepth), `Hang-off motion amplitude ${fmt(Math.hypot(xAmp, zAmp))} m exceeds 10 % of the hang-off height: outside the range of a linearised dynamic analysis.`);
    const reg = harmonic(base, { omega, xAmp, zAmp, rho: pp.rhoW, zeta: inp.zeta });
    need(reg.converged, 'The drag linearisation did not converge for the design wave; reduce the wave height.');
    const dT = reg.tensionAmp[0];
    // quasi-static reference: heave and surge solved statically and combined as 90° out of phase
    let qs = 0;
    if (a > 0) {
      const qh = zAmp > 1e-6 ? run(() => quasiStaticHeave(base, { heaveAmp: zAmp }), 'Quasi-static heave case').tensionAmplitude : 0;
      let qx = 0;
      if (xAmp > 1e-6) { const oc = run(() => riserOffsetCases(base, [-xAmp, xAmp]), 'Quasi-static surge case'); qx = 0.5 * Math.abs(oc[1].topTension - oc[0].topTension); }
      qs = Math.hypot(qh, qx);
    }
    const sAmp = stressAmp(reg, pp, inp.scf, mvOf(inp, base, xAmp, zAmp));
    const travel = tdpTravel(base, pp, xAmp, zAmp);
    const st = stressOf(base, pp, pressure);
    let vmDyn = 0, kV = 0, minTe = Infinity, kMin = 0;
    for (let k = 0; k <= n; k++) {
      const v = st.vonMises[k] + sAmp[k] / inp.scf; if (v > vmDyn) { vmDyn = v; kV = k; }
      const t = base.Te[k] - reg.tensionAmp[k]; if (t < minTe) { minTe = t; kMin = k; }
    }
    // touchdown zone = suspended part within a quarter of the suspended length of the touchdown point
    const sT = base.s[kT]; let sTdz = 0, kTdz = kT, sTop = 0;
    for (let k = 0; k <= n; k++) { if (Math.abs(base.s[k] - sT) < 0.25 * sT && sAmp[k] > sTdz) { sTdz = sAmp[k]; kTdz = k; } if (base.s[k] < 0.1 * sT) sTop = Math.max(sTop, sAmp[k]); }

    // irregular sea: transfer functions linearised at the significant wave amplitude
    const sw = sweep(base, pp, inp, inp.Hs / 2);
    const grid = seaGrid(sw, inp);
    const hT = sw.rows.map((r) => r.resp.tensionAmp[0] / r.a), hS = sw.rows.map((r) => stressAmp(r.resp, pp, inp.scf, r.mv)[kTdz] / r.a);
    const hZ = sw.rows.map((r) => r.heave), hX = sw.rows.map((r) => r.surge);
    const stT = responseStats(sw, grid, hT), stS = responseStats(sw, grid, hS), stZ = responseStats(sw, grid, hZ);
    const stV = responseStats(sw, grid, sw.rows.map((r) => r.heave * r.omega)), stA = responseStats(sw, grid, sw.rows.map((r) => r.heave * r.omega ** 2));
    const mpmT = mpm(stT.sigma, stT.Tz);
    const TtopMax = Ts + Math.max(dT, 0), daf = TtopMax / Ts, dafDyn = qs > 1e-6 ? dT / qs : 1;
    const fy = pp.mat.fy, allow = inp.Cf * inp.CaAllow * fy;
    const TeTdp = base.Te[kT];

    const kpis = [
      kpi('TtopStatic', 'Static top tension', Ts / 1e3, 'kN', 'Top effective tension of the riser at rest in the mean position; the reference for the dynamic amplification.'),
      kpi('TdynAmp', 'Dynamic top-tension amplitude (design wave)', dT / 1e3, 'kN', 'Amplitude of the cyclic top tension when the hang-off moves with the vessel in the regular design wave. It is driven mostly by vertical (heave + pitch) motion at the porch.'),
      kpi('TtopMax', 'Maximum top tension (design wave)', TtopMax / 1e3, 'kN', 'Static tension plus the dynamic amplitude in the design wave: the extreme axial load on the hang-off.'),
      kpi('TtopMin', 'Minimum top tension (design wave)', (Ts - dT) / 1e3, 'kN', 'Static tension minus the dynamic amplitude. If it approaches zero the riser goes slack at the top during downward heave.'),
      kpi('daf', 'Total amplification Tmax / Tstatic', daf, '-', 'Ratio of the extreme to the static top tension. Catenary risers on ship-shaped hosts typically show 1.1–1.5 in the 100-year wave; higher values point to heave-induced compression waves.'),
      kpi('dafDyn', 'Dynamic / quasi-static tension amplitude', dafDyn, '-', 'How much inertia and drag change the tension amplitude compared with moving the hang-off slowly. For a deep catenary riser the axial inertia of the suspended pipe dominates, so values of 10 or more are normal; the quasi-static method alone would badly under-predict tension variation.'),
      kpi('minTe', 'Minimum effective tension along the riser', minTe / 1e3, 'kN', 'Lowest value of static minus dynamic effective tension anywhere along the riser (usually near the touchdown). Negative means dynamic compression and possible buckling.', 0, 'min'),
      kpi('stressTdz', 'Dynamic stress amplitude at the touchdown zone', sTdz / 1e6, 'MPa', 'Amplitude of the cyclic hot-spot stress (tension plus bending, with SCF) at the worst point of the touchdown zone in the design wave. Twice this is the stress range of one wave cycle.'),
      kpi('stressTop', 'Dynamic stress amplitude near the hang-off', sTop / 1e6, 'MPa', 'Amplitude of the cyclic stress in the top 10 % of the riser; it is governed by axial tension variation and, with a stiff joint, by end bending.'),
      kpi('vmMax', 'Maximum von Mises stress incl. dynamics', vmDyn / 1e6, 'MPa', 'Static von Mises stress plus the dynamic stress amplitude at the same section (added in phase, without SCF), against the allowable for the extreme case.', allow / 1e6, 'max'),
      kpi('sigT', 'Top-tension standard deviation (irregular sea)', stT.sigma / 1e3, 'kN', 'Standard deviation of the dynamic top tension in the irregular sea state, from the tension transfer function and the JONSWAP spectrum.'),
      kpi('mpmT', 'Most probable maximum dynamic tension in 3 h', mpmT / 1e3, 'kN', 'Most probable largest dynamic top tension in a three-hour storm, σ·√(2·ln N). Compare with the design-wave amplitude: the larger governs.'),
      kpi('sigStressTdz', 'Touchdown stress standard deviation (irregular sea)', stS.sigma / 1e6, 'MPa', 'Standard deviation of the hot-spot stress at the critical touchdown node in the sea state; input to spectral fatigue.'),
      kpi('tdpTravel', 'Touchdown travel amplitude (design wave)', travel, 'm', 'Distance the touchdown point moves along the pipe either side of its mean position in the design wave (catenary geometry). The sag-bend curvature profile travels with it, which is what cycles the bending stress there.'),
      kpi('heaveAmp', 'Vertical hang-off motion amplitude (design wave)', zAmp, 'm', 'Vertical motion of the porch in the design wave from heave plus lever arm × pitch. It is the main driver of catenary riser dynamics.'),
      kpi('surgeAmp', 'Horizontal hang-off motion amplitude (design wave)', xAmp, 'm', 'In-plane horizontal first-order motion of the porch in the design wave.'),
      kpi('velHang', 'Vertical hang-off velocity amplitude', zAmp * omega, 'm/s', 'Peak vertical velocity of the porch. When the downward velocity approaches the terminal sinking velocity of the pipe, tension drops sharply and compression can develop at touchdown.'),
      kpi('accHang', 'Vertical hang-off acceleration amplitude', zAmp * omega * omega, 'm/s²', 'Peak vertical acceleration of the porch in the design wave; it scales the inertial part of the dynamic tension.'),
      kpi('sigAcc', 'Vertical hang-off acceleration, standard deviation', stA.sigma, 'm/s²', 'Standard deviation of porch vertical acceleration in the irregular sea state, for comparison with motion monitoring data.'),
    ];
    const checks = [
      check('Combined stress incl. dynamics (von Mises)', vmDyn / 1e6, allow / 1e6, 'MPa', 'API RP 2RD (1998) §5.2.3, Table 2 (extreme case Cf)', 'Guards against yielding under the static stress plus the dynamic increment in the design wave.'),
      check('Dynamic tension at touchdown vs static tension', reg.tensionAmp[kT] / 1e3, Math.max(TeTdp, 1) / 1e3, 'kN', 'DNV-ST-F201 (2018) Sec. 5 — effective compression to be avoided or assessed', 'Compares the dynamic tension amplitude at the touchdown point with the static tension there. Above 1.0 the pipe goes into effective compression each cycle and can buckle in the sag bend.'),
    ];

    const Ts_ = sw.rows.map((r) => r.T);
    const [sx, tAmp, mAmp, nAmp, sA] = thin([base.s, reg.tensionAmp, reg.momentAmp, reg.normalAmp, sAmp], 300);
    const plots = [
      plot('Top-tension transfer function', 'Wave period', 's', 'Dynamic top tension per metre of wave amplitude', 'kN/m', [{ name: 'Top tension RAO', x: Ts_, y: hT.map((v) => v / 1e3) }], 'Dynamic top tension per unit wave amplitude against wave period (drag linearised at the significant amplitude). Peaks show where vessel motion and riser dynamics combine most unfavourably.'),
      plot('Hang-off motion RAOs', 'Wave period', 's', 'Motion per metre of wave amplitude', 'm/m', [{ name: 'Vertical (heave + arm·pitch)', x: Ts_, y: hZ }, { name: 'Surge', x: Ts_, y: hX }], `Motion of the porch per metre of wave amplitude (${sw.source}). The vertical curve multiplied by the wave spectrum gives the motion that drives the riser.`),
      plot('Dynamic tension amplitude along the riser (design wave)', 'Arc length from hang-off', 'm', 'Tension amplitude', 'kN', [{ name: 'Dynamic amplitude', x: sx, y: scaled(tAmp, 1e-3) }, { name: 'Static effective tension', x: thin([base.s], 300)[0], y: scaled(thin([base.Te], 300)[0], 1e-3), dash: true }], 'Dynamic effective-tension amplitude compared with the static tension. Where the solid curve reaches the dashed one the pipe loses all tension during the cycle.'),
      plot('Dynamic bending moment amplitude (design wave)', 'Arc length from hang-off', 'm', 'Moment amplitude', 'kN·m', [{ name: 'Moment amplitude', x: sx, y: scaled(mAmp, 1e-3) }], 'Amplitude of the cyclic bending moment from the linearised (frozen-contact) solution. Its sharp spike at the static touchdown node is an artefact of freezing contact; the stress plot uses the selected touchdown bending model.'),
      plot('Dynamic stress amplitude (design wave)', 'Arc length from hang-off', 'm', 'Stress amplitude', 'MPa', [{ name: 'Hot-spot stress amplitude', x: sx, y: scaled(sA, 1e-6) }], 'Cyclic hot-spot stress amplitude including the SCF. The two fatigue-critical regions — hang-off and touchdown — stand out as peaks.'),
      plot('Motion envelope along the riser (design wave)', 'Arc length from hang-off', 'm', 'Normal displacement amplitude', 'm', [{ name: 'Normal displacement', x: sx, y: nAmp }], 'Amplitude of motion normal to the riser axis. It shows how the imposed hang-off motion travels down the riser and is damped by drag.'),
      plot('Wave and top-tension spectra', 'Frequency', 'rad/s', 'Spectral density (normalised to peak)', '-', (() => { const sT2 = grid.omega.map((w, j) => { let lo = 0; while (lo < sw.omega.length - 2 && sw.omega[lo + 1] < w) lo++; const h = hT[lo] + (hT[lo + 1] - hT[lo]) * Math.min(1, Math.max(0, (w - sw.omega[lo]) / (sw.omega[lo + 1] - sw.omega[lo]))); return grid.S[j] * h * h; }); const m1 = Math.max(...grid.S, 1e-30), m2 = Math.max(...sT2, 1e-30); return [{ name: 'Wave spectrum', x: grid.omega, y: grid.S.map((v) => v / m1) }, { name: 'Top-tension response spectrum', x: grid.omega, y: sT2.map((v) => v / m2) }]; })(), 'JONSWAP wave spectrum and the resulting top-tension spectrum, each scaled to its own peak. A shift between the peaks shows which wave periods actually load the riser.'),
    ];
    const tables = [{ title: 'Response sweep', columns: [{ key: 'T', label: 'Wave period', unit: 's' }, { key: 'heave', label: 'Vertical motion RAO', unit: 'm/m' }, { key: 'surge', label: 'Surge RAO', unit: 'm/m' }, { key: 'T1', label: 'Top tension RAO', unit: 'kN/m' }, { key: 'S1', label: 'Touchdown stress RAO', unit: 'MPa/m' }, { key: 'a', label: 'Linearisation amplitude', unit: 'm' }],
      rows: sw.rows.map((r, i) => ({ T: round(r.T, 4), heave: round(r.heave, 3), surge: round(r.surge, 3), T1: round(hT[i] / 1e3, 4), S1: round(hS[i] / 1e6, 4), a: round(r.a, 3) })),
      explain: 'Transfer functions at each solved wave period: porch motion, dynamic top tension and touchdown hot-spot stress per metre of wave amplitude, with the amplitude used to linearise drag.' }];
    const ucMax = Math.max(...checks.map((c) => c.uc));
    const dAng = 2;
    const advice = [
      ...vesselRangeAdvice(inp, Math.min(inp.Twave, inp.Tp)),
      adv('engineering', ucMax > 1 ? 'fail' : ucMax > 0.85 ? 'warn' : 'ok',
        `Design wave: dynamic tension ±${fmt(dT / 1e3)} kN on ${fmt(Ts / 1e3)} kN static (amplification ${fmt(daf)}); dynamic/quasi-static ratio ${fmt(dafDyn)}; minimum effective tension ${fmt(minTe / 1e3)} kN at s = ${fmt(base.s[kMin])} m; touchdown stress amplitude ${fmt(sTdz / 1e6)} MPa.`,
        checks[1].uc > 1 ? `Dynamic compression is predicted at touchdown. Move the porch closer to the pitch centre (lever arm now ${inp.hangArm} m), increase the hang-off angle by about ${dAng}° to raise bottom tension, or consider a lazy-wave configuration that decouples vessel heave; confirm with a nonlinear time-domain analysis.`
          : checks[0].uc > 1 ? `Reduce the dynamic stress: the governing section is at s = ${fmt(base.s[kV])} m. Reduce the vertical porch motion or upgrade the wall locally to about ${fmt(pp.t * 1e3 * vmDyn / (0.9 * allow))} mm.`
            : 'Use the tension and stress transfer functions in the wave-fatigue analysis, and repeat with the vessel RAOs from the diffraction analysis once available — the closed-form vessel model is only a starting point.',
        `API RP 2RD (1998, R2013) §5.5; DNV-ST-F201 (2018) App. A on frequency-domain methods. ${CONFIRM}`),
      adv('economic', 'info',
        `Porch vertical motion of ±${fmt(zAmp)} m produces ${fmt(dT / Math.max(zAmp, 1e-6) / 1e3)} kN of tension per metre. Hang-off location is a layout decision made early and at almost no cost; correcting dynamics later means a different riser concept or a larger hull.`,
        `Before buying pipe, check a porch position ${fmt(Math.max(0, inp.hangArm - 30))} m from midship instead of ${inp.hangArm} m: vertical motion from pitch falls in proportion to the lever arm and with it fatigue damage (roughly with the cube of stress range).`,
        'Value-of-information and front-end loading practice (AACE RP 18R-97 class 4 → 3); whole-life costing ISO 15686-5.',
        `Lower dynamic stress allows a thinner touchdown section and fewer upgraded welds: each millimetre of wall over the suspended length is ${fmt(steelTonnes(pp.D, 0.001, base.s[kT], pp.mat.rho))} t of steel. ${steelCo2Text(steelTonnes(pp.D, 0.001, base.s[kT], pp.mat.rho))}`),
    ];
    return finish({
      kpis, checks, plots, tables,
      notes: [`Vessel motions from the ${rao.source}; regular wave H = ${inp.Hmax} m, T = ${inp.Twave} s; irregular sea Hs = ${inp.Hs} m, Tp = ${inp.Tp} s.`,
        `Drag linearisation converged in ${reg.iterations} iterations. ${checks[0].uc >= checks[1].uc ? 'Combined stress' : 'Touchdown compression'} is the governing check (${fmt(ucMax)}).`],
      advice,
      publish: { 'riser.dynTensionAmp': dT, 'riser.daf': daf, 'riser.hangoffHeaveAmp': zAmp },
      quantities: { steel_t: steelTonnes(pp.D, pp.t, base.s[n], pp.mat.rho) },
    });
  },

  verify: [
    {
      // Zero wave height: no hang-off motion, so the maximum top tension is the static catenary value
      // T = w·h/(1 − sinθ) (h = 985 m, θ = 12°) and the dynamic amplitude vanishes.
      name: 'Zero-wave limit: dynamic amplitude → 0, top tension → static catenary w·h/(1 − sin θ)', source: 'Limiting case; Irvine, Cable Structures (1981) §1.2',
      inputs: { Hmax: 0, topType: 'pinned', pInt: 0, nElem: 320 },
      expect: [{ kpi: 'TtopMax', value: handWeight() * 985 / (1 - Math.sin(12 * DEG)) / 1e3, tol: 0.005 }, { kpi: 'TdynAmp', value: 0, abs: 1e-6 }, { kpi: 'daf', value: 1, tol: 1e-9 }],
    },
  ],
};
