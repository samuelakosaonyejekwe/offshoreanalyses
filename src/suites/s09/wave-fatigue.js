// Wave- and current-induced (first-order vessel motion) fatigue of a catenary riser by the spectral method:
// stress transfer functions along the riser, a sea-state scatter table, S–N curves and Miner summation.
import { narrowBandDamageRate, dirlikDamageRate } from '../../engine/fatigue/spectral.js';
import { cyclesFromConstants } from '../../engine/fatigue/sn.js';
import { jonswap } from '../../engine/metocean/spectra.js';
import {
  YEAR_S, need, fmt, round, thin, adv, check, kpi, plot, meshInput, finish, snInputs, snFrom, currentInputs, steelTonnes, steelCo2Text, vesselCo2Text, LIMIT_2D, CONFIRM,
} from './_common.js';
import { scrInputs, buildScr, currentOf, tdpNode } from './_scr.js';
import { vesselInputs, vesselRangeAdvice, LIMIT_VESSEL } from './_vessel.js';
import { dynInputs, stressAmp, sweep, interpAt } from './_dyn.js';

export default {
  id: 's09-wave-fatigue', suite: 9, title: 'Wave- and current-induced fatigue',
  summary: 'Long-term fatigue damage and life along a catenary riser from vessel motion in a sea-state scatter table, by spectral analysis with Dirlik and narrow-band cycle counting.',
  covers: [6], fidelity: 'reduced-order',
  standards: ['DNV-RP-C203 (2021) §2.4 S–N curves, §5 simplified and spectral fatigue', 'DNV-ST-F201 (2018) Sec. 5 fatigue limit state', 'DNV-RP-F204 (2019) riser fatigue', 'API RP 2RD (1998, R2013) §5.6'],
  models: ['Spectral fatigue (stress transfer function × JONSWAP spectrum)', "Dirlik's rainflow-range distribution", 'Narrow-band (Rayleigh) damage', 'S–N curves with thickness correction and SCF', 'Palmgren–Miner summation', 'Response amplitude operators', 'Linearised Morison drag'],
  limits: `${LIMIT_2D} Stress transfer functions are linearised once, at a damage-weighted wave amplitude, and applied to every sea state; tension and bending stress are added in phase (conservative). Mean current enters only through the static shape; wave directionality, second-order (slow-drift) motion, VIV and installation damage are not included here — combine them in the fatigue-life and integrity module. Rainflow counting of nonlinear time histories is outside this module.` + LIMIT_VESSEL,
  publishes: ['riser.waveDamagePerYear', 'riser.waveFatigueLife', 'riser.fatigueStressRange', 'riser.fatigueCyclesPerYear'],
  searchTerms: 'steel catenary riser wave fatigue spectral touchdown Dirlik',
  inputs: [
    ...scrInputs(),
    ...currentInputs({ Us: 0.3, Ub: 0.05, bind: 'site.Uc1' }),
    ...vesselInputs(),
    { key: 'scatter', label: 'Sea-state scatter table', type: 'table', group: 'Wave climate', role: 'bc',
      columns: [{ key: 'Hs', label: 'Hs', unit: 'm', type: 'number' }, { key: 'Tp', label: 'Tp', unit: 's', type: 'number' }, { key: 'prob', label: 'Occurrence', unit: '%', type: 'number' }],
      default: [{ Hs: 0.5, Tp: 7, prob: 20 }, { Hs: 1, Tp: 8, prob: 34 }, { Hs: 1.5, Tp: 9, prob: 24 }, { Hs: 2, Tp: 10, prob: 12 }, { Hs: 2.5, Tp: 11, prob: 6 }, { Hs: 3.5, Tp: 12, prob: 3 }, { Hs: 5, Tp: 13, prob: 0.8 }, { Hs: 7, Tp: 14, prob: 0.2 }],
      help: 'Long-term distribution of sea states (lumped blocks). Occurrences should sum to 100 %; any shortfall is treated as calm.' },
    { key: 'gamma', label: 'JONSWAP peak factor', unit: '-', type: 'number', default: 2, min: 1, max: 7, group: 'Wave climate', role: 'bc' },
    ...dynInputs(),
    ...snInputs({ curve: 'E' }),
    { key: 'cycleMethod', label: 'Cycle-counting method', type: 'select', default: 'dirlik', group: 'Fatigue resistance', role: 'control', options: [{ value: 'dirlik', label: 'Dirlik (wide-band)' }, { value: 'narrow', label: 'Narrow-band Rayleigh (conservative)' }] },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 10, min: 1, max: 20, group: 'Acceptance criteria', role: 'criteria',
      help: `DNV-ST-F201 Table 5-11: 3 / 6 / 10 for low / normal / high safety class; API RP 2RD uses 10 for non-inspectable risers. ${CONFIRM}` },
    { key: 'designLife', label: 'Design life', unit: 'yr', type: 'number', default: 25, min: 1, max: 60, group: 'Acceptance criteria', role: 'criteria', bind: 'life.designYears' },
    meshInput(160),
  ],
  mesh: { key: 'nElem', label: 'Beam elements', levels: [80, 160, 320, 640], qoi: 'sigRef' },

  compute(inp, ctx) {
    // engineering units → SI in buildScr; occurrences in % → fraction
    const states = (inp.scatter || []).filter((r) => r.Hs > 0 && r.Tp > 0 && r.prob > 0).map((r) => ({ Hs: r.Hs, Tp: r.Tp, p: r.prob / 100 }));
    need(states.length > 0, 'The scatter table needs at least one sea state with positive Hs, Tp and occurrence.');
    const pSum = states.reduce((a, r) => a + r.p, 0);
    need(pSum <= 1.0001, `Sea-state occurrences sum to ${fmt(pSum * 100)} %; they cannot exceed 100 %.`);
    const { pp, base } = buildScr(inp, ctx, { nElem: inp.nElem, current: currentOf(inp, 1) });
    const n = base.s.length - 1, kT = tdpNode(base);
    const { sn, label } = snFrom(inp, pp.t);
    // damage-weighted linearisation height: (Σ p·Hs^m / Σ p)^(1/m) with the S–N slope m
    const HsLin = (states.reduce((a, r) => a + r.p * r.Hs ** sn.m1, 0) / pSum) ** (1 / sn.m1);
    const sw = sweep(base, pp, inp, HsLin / 2);
    const H = sw.rows.map((r) => { const a = stressAmp(r.resp, pp, 1, r.mv); for (let k = 0; k <= n; k++) a[k] /= r.a; return a; });   // Pa per m, nominal
    // wave spectra on a common grid
    const m = 80, w0 = sw.omega[0], w1 = sw.omega[sw.omega.length - 1], dw = (w1 - w0) / m, wg = [];
    for (let j = 0; j < m; j++) wg.push(w0 + (j + 0.5) * dw);
    const Sg = states.map((s) => wg.map((w) => jonswap(w, { Hs: s.Hs, Tp: s.Tp, gamma: inp.gamma })));
    const mom = (k, i) => { let m0 = 0, m1 = 0, m2 = 0, m4 = 0; const col = H.map((h) => h[k]); for (let j = 0; j < m; j++) { const h = interpAt(sw.omega, col, wg[j]), e = Sg[i][j] * h * h * dw, f = wg[j] / (2 * Math.PI); m0 += e; m1 += e * f; m2 += e * f * f; m4 += e * f ** 4; } return { m0, m1, m2, m4 }; };
    const rateOf = (mo) => (mo.m0 > 0 ? (inp.cycleMethod === 'dirlik' ? dirlikDamageRate({ moments: mo, sn }) : narrowBandDamageRate({ sigma: Math.sqrt(mo.m0), rate: Math.sqrt(mo.m2 / mo.m0), sn })) : 0);
    // narrow-band damage at every node (cheap closed form) to locate the hot spots
    const dNb = new Float64Array(n + 1);
    for (let k = 0; k <= n; k++) for (let i = 0; i < states.length; i++) { const mo = mom(k, i); if (mo.m0 > 0) dNb[k] += states[i].p * narrowBandDamageRate({ sigma: Math.sqrt(mo.m0), rate: Math.sqrt(mo.m2 / mo.m0), sn }) * YEAR_S; }
    const sT = base.s[kT]; let kTdz = kT, kTop = 0;
    for (let k = 0; k <= n; k++) { if (base.s[k] > 0.5 * sT && dNb[k] > dNb[kTdz]) kTdz = k; if (base.s[k] <= 0.5 * sT && dNb[k] > dNb[kTop]) kTop = k; }
    const at = (k) => { const per = states.map((s, i) => { const mo = mom(k, i); return { ...s, sigma: Math.sqrt(mo.m0), nu0: mo.m0 > 0 ? Math.sqrt(mo.m2 / mo.m0) : 0, d: s.p * rateOf(mo) * YEAR_S }; }); return { per, d: per.reduce((a, r) => a + r.d, 0) }; };
    const tdz = at(kTdz), top = at(kTop);
    const dMax = Math.max(tdz.d, top.d), crit = tdz.d >= top.d ? tdz : top, kCrit = tdz.d >= top.d ? kTdz : kTop;
    need(dMax > 0, 'No fatigue damage is predicted: every stress range is below the fatigue limit of the selected S–N curve. Check the scatter table and the vessel motions.');
    const life = 1 / dMax, need_ = inp.dff * inp.designLife;
    const iRef = states.reduce((b, s, i) => (crit.per[i].d > crit.per[b].d ? i : b), 0);
    // equivalent constant-amplitude stress range and cycles per year at the critical node (m = m1), for crack growth
    const cyc = crit.per.reduce((a, r) => a + r.p * r.nu0 * YEAR_S, 0);
    const sEq = cyc > 0 ? (dMax * sn.a1 / cyc) ** (1 / sn.m1) : 0;
    const ratio = dNb[kCrit] > 0 ? dMax / dNb[kCrit] : 1;

    const kpis = [
      kpi('lifeMin', 'Fatigue life (unfactored)', life, 'yr', 'Years to a Miner sum of 1.0 at the most damaged location with the characteristic (mean minus two standard deviations) S–N curve, before the design fatigue factor.', need_, 'min'),
      kpi('lifeTdz', 'Fatigue life at the touchdown zone', tdz.d > 0 ? 1 / tdz.d : 1e9, 'yr', 'Unfactored life at the worst node in the lower half of the suspended length, where the moving touchdown point bends the pipe every wave cycle.'),
      kpi('lifeTop', 'Fatigue life near the hang-off', top.d > 0 ? 1 / top.d : 1e9, 'yr', 'Unfactored life at the worst node in the upper half of the riser, driven by tension variation and end bending below the joint.'),
      kpi('damageYr', 'Annual damage at the critical location', dMax, '1/yr', 'Miner damage accumulated per year at the critical location. Multiply by the service years and the design fatigue factor to compare with 1.0.'),
      kpi('factoredLife', 'Life divided by the design fatigue factor', life / inp.dff, 'yr', 'Fatigue life after dividing by the design fatigue factor; it must be at least the design life.', inp.designLife, 'min'),
      kpi('sCrit', 'Critical location (arc length from hang-off)', base.s[kCrit], 'm', 'Arc length of the most damaged section. Use it to place the high-quality welds and to target inspection.'),
      kpi('sigRef', 'Stress standard deviation in the most damaging sea state', crit.per[iRef].sigma * inp.scf / 1e6, 'MPa', `Hot-spot stress standard deviation (with SCF) at the critical location in the sea state that contributes most damage (Hs = ${states[iRef].Hs} m).`),
      kpi('sEq', 'Equivalent constant stress range', sEq * inp.scf / 1e6, 'MPa', 'Constant-amplitude hot-spot stress range that would give the same annual damage with the same number of cycles on the first S–N slope; used for crack-growth assessment.'),
      kpi('cyclesYr', 'Stress cycles per year', cyc, '1/yr', 'Number of zero-upcrossing stress cycles per year at the critical location, weighted over the scatter table.'),
      kpi('nRef', 'S–N cycles to failure at 100 MPa nominal range', cyclesFromConstants(sn, 100e6), '-', 'Cycles to failure of the selected S–N curve (with SCF and thickness correction) at a 100 MPa nominal stress range: a quick check that the intended curve is in use.'),
      kpi('HsLin', 'Linearisation wave height', HsLin, 'm', 'Damage-weighted significant wave height at which the drag damping was linearised. Sea states far from it are slightly over- or under-damped.'),
    ];
    const checks = [check('Fatigue life with design fatigue factor', need_, life, 'yr', 'DNV-ST-F201 (2018) Sec. 5 F; DNV-RP-C203 (2021) §2.2; API RP 2RD §5.6', 'Requires the calculated life to exceed the design life times the design fatigue factor. If it does not, a through-wall fatigue crack is credible within the service life at the critical weld.')];

    const [sx, dy] = thin([base.s, dNb], 300);
    const plots = [
      plot('Annual fatigue damage along the riser', 'Arc length from hang-off', 'm', 'Damage per year (narrow-band)', '1/yr', [{ name: 'Annual damage', x: sx, y: dy.map((v) => v * ratio) }, { name: 'Allowable 1/(DFF·life)', x: [0, base.s[n]], y: [1 / need_, 1 / need_], dash: true }], 'Miner damage per year along the riser. The peaks mark the fatigue-critical zones (touchdown and hang-off); anything above the dashed line fails the factored-life requirement.'),
      plot('Damage contribution by sea state', 'Sea state', '', 'Share of damage', '%', [], 'Share of the annual damage at the critical location from each scatter block. Moderate, frequent sea states usually dominate, not the rare storms.', 'bar'),
      plot('Hot-spot stress transfer function at the critical location', 'Wave period', 's', 'Stress amplitude per metre of wave amplitude', 'MPa/m', [{ name: 'Critical location', x: sw.rows.map((r) => r.T), y: H.map((h) => h[kCrit] * inp.scf / 1e6) }], 'Stress amplitude at the critical node per metre of wave amplitude against wave period. Sea states whose energy sits under its peak do the damage.'),
      plot('S–N curve and operating stress ranges', 'Cycles to failure', '-', 'Stress range', 'MPa', (() => { const xs = [], ys = []; for (let e = 4; e <= 9; e += 0.25) { const N = 10 ** e; const S = N <= (sn.sKnee > 0 ? sn.a1 / sn.sKnee ** sn.m1 : Infinity) ? (sn.a1 / N) ** (1 / sn.m1) : (sn.a2 / N) ** (1 / sn.m2); xs.push(e); ys.push(S / 1e6); } return [{ name: `${label} incl. SCF (log10 N)`, x: xs, y: ys }, { name: 'Equivalent range at annual cycles × DFF × life', x: [Math.log10(Math.max(cyc * need_, 10))], y: [sEq / 1e6] }]; })(), 'Selected S–N curve (nominal stress range against log10 of cycles) with the equivalent operating point. The point must lie below the curve for the factored life to be met.', 'scatter'),
    ];
    plots[1] = { type: 'bar', title: 'Damage contribution by sea state', categories: states.map((s) => `${s.Hs} m / ${s.Tp} s`), x: { label: 'Sea state (Hs / Tp)', unit: '' }, y: { label: 'Share of damage', unit: '%' }, series: [{ name: 'Critical location', y: crit.per.map((r) => 100 * r.d / dMax) }], explain: plots[1].explain };
    const tables = [{ title: 'Sea-state damage at the critical location', columns: [{ key: 'Hs', label: 'Hs', unit: 'm' }, { key: 'Tp', label: 'Tp', unit: 's' }, { key: 'p', label: 'Occurrence', unit: '%' }, { key: 'sigma', label: 'Hot-spot stress σ', unit: 'MPa' }, { key: 'Tz', label: 'Stress Tz', unit: 's' }, { key: 'd', label: 'Annual damage', unit: '1/yr' }, { key: 'share', label: 'Share', unit: '%' }],
      rows: crit.per.map((r) => ({ Hs: r.Hs, Tp: r.Tp, p: round(r.p * 100, 4), sigma: round(r.sigma * inp.scf / 1e6, 4), Tz: round(r.nu0 > 0 ? 1 / r.nu0 : 0, 3), d: Number(r.d.toExponential(3)), share: round(100 * r.d / dMax, 3) })),
      explain: 'One row per scatter block: hot-spot stress standard deviation, cycle period and the damage it contributes per year at the critical location.' }];
    const uc = checks[0].uc;
    // stress reduction needed: damage ∝ S^m, so S must fall by uc^(1/m); for bending-dominated stress σ ∝ D/(2I) ~ 1/t
    const sFac = uc ** (1 / sn.m1), tNew = pp.t * sFac, dSteel = steelTonnes(pp.D, tNew - pp.t, 0.3 * sT, pp.mat.rho);
    const advice = [
      ...vesselRangeAdvice(inp, crit.per.reduce((a, r) => (r.d > a.d ? r : a), crit.per[0]).Tp),
      adv('engineering', uc > 1 ? 'fail' : uc > 0.7 ? 'warn' : 'ok',
        `Unfactored life ${fmt(life)} yr at s = ${fmt(base.s[kCrit])} m (${kCrit === kTdz ? 'touchdown zone' : 'hang-off region'}); required ${fmt(need_)} yr (DFF ${inp.dff} × ${inp.designLife} yr). ${label}, SCF ${inp.scf}.`,
        uc > 1 ? `Stress ranges must fall by ${fmt((1 - 1 / sFac) * 100)} % (damage ∝ Δσ^${sn.m1}). Options: a better weld class or flush-ground welds in the critical zone, an SCF reduced from ${inp.scf} by tighter fit-up, a wall of about ${fmt(tNew * 1e3)} mm locally, a larger hang-off angle, or a porch nearer the vessel pitch centre.`
          : `Margin on life is ${fmt(1 / uc)}. Before relying on it, add VIV and installation damage in the fatigue-life module and confirm the S–N class with the welding procedure qualification.`,
        `DNV-RP-C203 (2021) §2.4 and §5; DNV-ST-F201 (2018) Sec. 5 F. ${CONFIRM}`),
      adv('economic', uc > 1 ? 'warn' : 'info',
        `Damage is concentrated near s = ${fmt(base.s[kCrit])} m. Upgrading the wall to ${fmt(Math.max(tNew, pp.t) * 1e3)} mm over the critical 30 % of the suspended length changes steel by ${fmt(dSteel)} t; a higher weld class costs inspection and grinding time per weld instead.`,
        'Compare a local weld-quality upgrade (per-weld cost, no added mass) with a local wall upgrade (added steel and welding time) and with a lower design fatigue factor justified by an inspection programme.',
        'Whole-life costing (ISO 15686-5) and risk-based inspection planning (DNV-RP-G101 / API RP 580).',
        `Local wall upgrade: ${steelCo2Text(dSteel)} One additional in-service inspection campaign: ${vesselCo2Text(2)}`),
    ];
    return finish({
      kpis, checks, plots, tables,
      notes: [`${states.length} sea states (${fmt(pSum * 100)} % occurrence); transfer functions at ${inp.nFreq} periods from the ${sw.source}, linearised at Hs = ${fmt(HsLin)} m.`,
        `${inp.cycleMethod === 'dirlik' ? 'Dirlik' : 'Narrow-band'} counting at the critical nodes; the along-riser curve is the narrow-band profile scaled to the critical value (factor ${fmt(ratio)}).`],
      advice,
      publish: { 'riser.waveDamagePerYear': dMax, 'riser.waveFatigueLife': life, 'riser.fatigueStressRange': sEq * inp.scf, 'riser.fatigueCyclesPerYear': cyc },
      quantities: { steel_t: Math.max(0, dSteel) },
    });
  },

  verify: [
    {
      // DNV-RP-C203 Table 2-1, curve D in air: log a = 12.164, m = 3 (N ≤ 1e7), no thickness correction for
      // t = 20.6 mm < t_ref = 25 mm. At Δσ = 100 MPa: N = 10^12.164 / 100³ = 1.4588e6 cycles.
      name: 'S–N curve D in air: N(100 MPa) = 10^12.164 / 100³', source: 'DNV-RP-C203 (2021) Table 2-1',
      inputs: { snCurve: 'D', snEnv: 'air', scf: 1 },
      expect: [{ kpi: 'nRef', value: 10 ** 12.164 / 1e6, tol: 0.001 }],
    },
    {
      // Same curve with SCF = 1.5: the hot-spot range is 150 MPa, N = 10^12.164 / 150³ = 4.3224e5 cycles.
      name: 'SCF folded into the S–N constants: N = 10^12.164 / (1.5·100)³', source: 'DNV-RP-C203 (2021) §2.4, Table 2-1',
      inputs: { snCurve: 'D', snEnv: 'air', scf: 1.5 },
      expect: [{ kpi: 'nRef', value: 10 ** 12.164 / 150 ** 3, tol: 0.001 }],
    },
  ],
};
