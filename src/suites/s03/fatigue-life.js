// Fatigue damage and remaining useful life of the support structure: long-term damage over the joint wind–wave
// climate by frequency-domain response and Dirlik rainflow statistics (default) or by rainflow counting of coupled
// time-domain simulations (opt-in), S–N curves with thickness and stress-concentration effects, Palmgren–Miner sum.
import { snCurve, snConstants, cyclesToFailure, SN_CURVE_NAMES } from '../../engine/fatigue/sn.js';
import { spectralDamage } from '../../engine/fatigue/spectral.js';
import { minerSum, fatigueLife } from '../../engine/fatigue/damage.js';
import { ntmSigma } from '../../engine/aero/turbulence.js';
import { windWaveBins, runDesignLoadCases } from '../../engine/aero/dlc.js';
import { rainflowCount } from '../../engine/aero/coupled.js';
import { YEAR_S, need, run, fmt, round, adv, check, sev, steelCo2Text, thin, tube, turbineInputs, supportInputs, foundationInputs, steelInputs, steelOf, turbine, LIMIT_AERO, LIMIT_COUPLED, LIMIT_CODE, LIMIT_HIFI, CONFIRM } from './_common.js';
import { dynamicsInputs, supportStructure, modalSpectra } from './_dynamics.js';

/** Sea state conditional on hub wind speed from the user table (linear interpolation, held beyond the ends). */
function seaTable(rows) {
  const tab = (rows || []).map((r) => ({ U: Number(r.U), Hs: Number(r.Hs), Tp: Number(r.Tp) })).filter((r) => Number.isFinite(r.U) && r.Hs >= 0 && r.Tp > 0).sort((a, b) => a.U - b.U);
  need(tab.length >= 2, 'The wind–wave table needs at least two rows (hub wind speed, Hs, Tp).');
  for (let i = 1; i < tab.length; i++) need(tab[i].U > tab[i - 1].U, 'The wind–wave table has repeated wind speeds.');
  return tab;
}
const seaAt = (tab, U) => {
  if (U <= tab[0].U) return tab[0];
  if (U >= tab[tab.length - 1].U) return tab[tab.length - 1];
  let i = 0;
  while (tab[i + 1].U < U) i++;
  const w = (U - tab[i].U) / (tab[i + 1].U - tab[i].U);
  return { Hs: tab[i].Hs + w * (tab[i + 1].Hs - tab[i].Hs), Tp: tab[i].Tp + w * (tab[i + 1].Tp - tab[i].Tp) };
};

export default {
  id: 's03-fatigue-life', suite: 3, title: 'Fatigue and remaining useful life',
  summary: 'Lifetime fatigue damage at the mudline or tower-base weld from the joint wind–wave climate, the fatigue life with the design fatigue factor, and the remaining useful life given the years already in service.',
  covers: [5], fidelity: 'reduced-order',
  standards: ['DNV-RP-C203 (2021) §2 (S–N curves, thickness effect), §5', 'DNV-ST-0126 (2021) §4.11 (design fatigue factors)', 'IEC 61400-3-1:2019 §7.4 DLC 1.2, 6.4, 7.2', 'IEC 61400-1:2019 §7.6.3', 'DNV-ST-0262 (2016) lifetime extension of wind turbines', 'Dirlik (1985)'],
  models: ['Wind–wave joint probability: Weibull wind-speed bins with conditional sea states', 'Spectral fatigue analysis: Dirlik rainflow-range distribution on the response spectrum', 'Rainflow cycle counting of time-domain simulations (opt-in)', 'Palmgren–Miner damage accumulation', 'DNV-RP-C203 S–N relationships with thickness correction and stress concentration factor', 'Kaimal turbulence, JONSWAP waves, MacCamy–Fuchs inertia loading, aerodynamic damping', 'Idling / unavailable time without aerodynamic damping'],
  limits: `Single fore-aft mode with wind and waves aligned (conservative for monopiles), no 1P/3P harmonics, no start/stop transients and no driving-induced damage (see the pile-driving module). Wave drag is not linearised. The optional time-domain check uses a simplified DLC 1.2 with one short seed per bin. ${LIMIT_COUPLED} ${LIMIT_AERO} ${LIMIT_CODE} ${LIMIT_HIFI}`,
  publishes: ['life.fatigueYears'],
  searchTerms: 'monopile fatigue life wind wave scatter Dirlik S-N curve lifetime extension',
  inputs: [
    ...turbineInputs(),
    ...supportInputs({ embed: false }),
    ...steelInputs(),
    ...foundationInputs(),
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'weibullA', label: 'Weibull scale of hub wind speed', unit: 'm/s', type: 'number', default: 10.5, min: 4, max: 16, group: 'Wind and wave climate', role: 'bc', bind: 'site.windWeibullA', uncertain: { dist: 'normal', cov: 0.05 } },
    { key: 'weibullK', label: 'Weibull shape', unit: '-', type: 'number', default: 2.1, min: 1.2, max: 3.5, group: 'Wind and wave climate', role: 'bc', bind: 'site.windWeibullK' },
    { key: 'Iref', label: 'Reference turbulence intensity', unit: '-', type: 'number', default: 0.14, min: 0.05, max: 0.25, group: 'Wind and wave climate', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 } },
    { key: 'seaStates', label: 'Sea state against hub wind speed', type: 'table', group: 'Wind and wave climate', role: 'bc',
      columns: [{ key: 'U', label: 'Hub wind speed', unit: 'm/s', type: 'number' }, { key: 'Hs', label: 'Hs', unit: 'm', type: 'number' }, { key: 'Tp', label: 'Tp', unit: 's', type: 'number' }],
      default: [{ U: 2, Hs: 0.6, Tp: 5.0 }, { U: 6, Hs: 0.9, Tp: 5.3 }, { U: 10, Hs: 1.4, Tp: 5.8 }, { U: 14, Hs: 2.1, Tp: 6.6 }, { U: 18, Hs: 3.0, Tp: 7.4 }, { U: 22, Hs: 3.9, Tp: 8.2 }, { U: 26, Hs: 4.8, Tp: 9.0 }, { U: 34, Hs: 6.5, Tp: 10.3 }],
      help: 'Expected Hs and Tp for each wind speed from the site wind–wave scatter diagram (lumped sea states). Replace the generic values with site data.' },
    { key: 'gamma', label: 'JONSWAP peak enhancement', unit: '-', type: 'number', default: 3.3, min: 1, max: 7, group: 'Wind and wave climate', role: 'bc' },
    { key: 'avail', label: 'Turbine availability', unit: '%', type: 'number', default: 95, min: 50, max: 100, group: 'Wind and wave climate', role: 'operation', help: 'Share of the time within the operating wind range that the rotor is producing. The rest is idling without aerodynamic damping.' },
    ...dynamicsInputs(),
    { key: 'section', label: 'Section assessed', type: 'select', default: 'mudline', group: 'Detail', role: 'geometry', options: [{ value: 'mudline', label: 'Monopile girth weld at the mudline' }, { value: 'towerBase', label: 'Tower-base girth weld / flange' }] },
    { key: 'curve', label: 'S–N curve (DNV-RP-C203)', type: 'select', default: 'C1', group: 'Detail', role: 'material', options: SN_CURVE_NAMES.map((n) => ({ value: n, label: n })), help: 'C1 for girth welds ground flush on both sides (usual for monopile cans), D for as-welded double-sided girth welds, E or F for attachments.' },
    { key: 'env', label: 'Environment', type: 'select', default: 'seawater-cp', group: 'Detail', role: 'material', options: [{ value: 'air', label: 'In air / coated' }, { value: 'seawater-cp', label: 'Sea water with cathodic protection' }, { value: 'free-corrosion', label: 'Free corrosion' }] },
    { key: 'scf', label: 'Stress concentration factor', unit: '-', type: 'number', default: 1.15, min: 1, max: 5, group: 'Detail', role: 'geometry', uncertain: { dist: 'lognormal', cov: 0.1 }, help: 'Girth-weld misalignment SCF = 1 + 3δ/t·exp(−√(t/D)) (DNV-RP-C203 §3.3.7); typically 1.05–1.3 for monopile cans.' },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 3, min: 1, max: 10, group: 'Acceptance criteria', role: 'criteria', help: `DNV-ST-0126: 3 for non-inspectable submerged structure, 2 in the splash zone when inspectable, 1 accessible in air. ${CONFIRM}` },
    { key: 'life', label: 'Required service life', unit: 'years', type: 'number', default: 25, min: 5, max: 60, group: 'Acceptance criteria', role: 'criteria', bind: 'life.designYears' },
    { key: 'yearsUsed', label: 'Years already in service', unit: 'years', type: 'number', default: 0, min: 0, max: 60, group: 'Existing condition', role: 'ic', help: 'Initial damage state: damage is accumulated at the computed rate for this period.' },
    { key: 'd0', label: 'Additional damage already consumed', unit: '-', type: 'number', default: 0.02, min: 0, max: 1, group: 'Existing condition', role: 'ic', help: 'Miner damage from transport, pile driving and installation (characteristic, before the fatigue factor).' },
    { key: 'sRef', label: 'Reference stress range for the S–N readout', unit: 'MPa', type: 'number', default: 100, min: 5, max: 500, group: 'Detail', role: 'control' },
    { key: 'method', label: 'Method', type: 'select', default: 'spectral', group: 'Numerics', role: 'control', options: [{ value: 'spectral', label: 'Frequency domain (Dirlik)' }, { value: 'time', label: 'Also run time-domain rainflow on the coupled simulator (slow)' }],
      help: 'The time-domain option runs one coupled simulation per wind bin (design load case 1.2) and reports its damage beside the spectral value.' },
    { key: 'tdDur', label: 'Time-domain record per bin', unit: 's', type: 'number', default: 150, min: 60, max: 3600, group: 'Numerics', role: 'control', showIf: { key: 'method', in: ['time'] } },
    { key: 'binW', label: 'Wind-speed bin width', unit: 'm/s', type: 'number', default: 2, min: 0.5, max: 4, group: 'Numerics', role: 'control' },
    { key: 'nFreq', label: 'Frequency points', type: 'int', default: 240, min: 60, max: 2000, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nFreq', label: 'Frequency points', levels: [120, 240, 480, 960], qoi: 'damageYr' },
  compute(inp, ctx) {
    const steel = steelOf(inp, ctx);
    const tb = turbine(inp);
    const str = supportStructure(inp, tb, steel, { nElem: 30, nModes: 2 });
    const md = str.modes[0];
    const tab = seaTable(inp.seaStates);
    const atMud = inp.section === 'mudline';
    const sec = atMud ? tube(inp.D, inp.t / 1e3) : tube(inp.towerDb, inp.towerTb / 1e3);
    // Section moment relative to the mudline moment: quasi-static lever ratio for the tower base (thrust-dominated).
    const mRatio = atMud ? 1 : Math.abs(md.momentBase / md.momentMudline);
    const curve = run(() => snCurve(inp.curve, inp.env), 'S–N curve');
    const sn = snConstants(curve, { thickness: sec.t, scf: inp.scf });
    const cutIn = tb.rotor.cutIn, cutOut = tb.rotor.cutOut;
    const prodBins = run(() => windWaveBins({ weibullA: inp.weibullA, weibullK: inp.weibullK, vIn: cutIn, vOut: cutOut, binWidth: inp.binW, hubHeight: inp.hubHeight, seaState: tab }), 'Wind–wave bins');
    const cdf = (v) => 1 - Math.exp(-((v / inp.weibullA) ** inp.weibullK));
    const av = inp.avail / 100;
    // Idle conditions: calm below cut-in and storm above cut-out (parked rotor, extreme-wind turbulence 0.11 V).
    const idle = [{ lo: 0, hi: cutIn }, { lo: cutOut, hi: cutOut + 5 }, { lo: cutOut + 5, hi: cutOut + 12 }, { lo: cutOut + 12, hi: 80 }]
      .map((b) => ({ ...b, U: b.hi === 80 ? b.lo + 4 : 0.5 * (b.lo + b.hi), prob: cdf(b.hi) - cdf(b.lo) }));
    const states = [];
    const rate = (V, sigma1, ss, parked) => {
      const sp = modalSpectra({ str, tb, inp, V: Math.max(V, 1), sigma1, Hs: ss.Hs, Tp: ss.Tp, gamma: inp.gamma, nFreq: inp.nFreq, parked });
      const S = Float64Array.from(sp.SM, (v) => v * (mRatio / sec.W) ** 2);
      const d = spectralDamage({ f: sp.f, S, sn, duration: 1, method: 'dirlik' });
      return { rate: d.damageRate, sigma: sp.sigmaM * mRatio / sec.W, nu0: sp.nu0, zeta: sp.zetaTotal, sp };
    };
    let dYr = 0, sum4 = 0;
    const add = (label, U, prob, r, ss, mode) => {
      const dy = prob * r.rate * YEAR_S;
      dYr += dy; sum4 += prob * YEAR_S * r.nu0 * (2 * Math.SQRT2 * r.sigma) ** 4 * 2; // narrow-band E[S⁴] = (2√2σ)⁴ Γ(3)
      states.push({ label, U, prob, Hs: ss.Hs, Tp: ss.Tp, mode, sigma: r.sigma, zeta: r.zeta, dy });
    };
    let rated = null;
    for (const b of prodBins) {
      const ss = { Hs: b.Hs, Tp: b.Tp };
      const rp = rate(b.U, ntmSigma({ Vhub: b.U, Iref: inp.Iref }), ss, false);
      if (!rated || Math.abs(b.U - tb.windAtMaxThrust) < Math.abs(rated.U - tb.windAtMaxThrust)) rated = { U: b.U, sp: rp.sp };
      add(`${fmt(b.lo)}–${fmt(b.hi)} m/s`, b.U, b.prob * av, rp, ss, 'production');
      if (av < 1) add(`${fmt(b.lo)}–${fmt(b.hi)} m/s`, b.U, b.prob * (1 - av), rate(b.U, ntmSigma({ Vhub: b.U, Iref: inp.Iref }), ss, true), ss, 'idling (unavailable)');
    }
    for (const b of idle) {
      if (b.prob < 1e-9) continue;
      const ss = seaAt(tab, b.U);
      add(`${fmt(b.lo)}–${b.hi === 80 ? '' : fmt(b.hi)} m/s`, b.U, b.prob, rate(b.U, b.U > cutOut ? 0.11 * b.U : ntmSigma({ Vhub: b.U, Iref: inp.Iref }), ss, true), ss, b.U > cutOut ? 'parked (storm)' : 'idling (calm)');
    }
    need(Number.isFinite(dYr) && dYr >= 0, 'Fatigue damage could not be evaluated: check the climate table and the S–N inputs.');
    const dUsed = inp.d0 + dYr * inp.yearsUsed;
    const fl = fatigueLife({ damageRate: Math.max(dYr, 1e-30), designLife: Math.max(inp.life - inp.yearsUsed, 0), dff: inp.dff, accumulatedDamage: dUsed });
    const lifeDesign = Math.min(fl.designLife, 1e4), remaining = Math.max(Math.min(fl.remainingLife, 1e4), 0);
    const totalLife = Math.min((1 / inp.dff - inp.d0) / Math.max(dYr, 1e-30), 1e4);
    const ucLife = inp.dff * (inp.d0 + dYr * inp.life);
    const delS = (sum4 * inp.life / 1e7) ** 0.25; // stress range, nominal
    const prodFrac = cdf(cutOut) - cdf(cutIn);
    // Optional time-domain cross-check through the design-load-case driver (DLC 1.2, production only).
    let td = null;
    if (inp.method === 'time') {
      const res = run(() => runDesignLoadCases({
        model: { surface: tb.surface, turbine: { inertia: tb.rotor.drivetrainInertia, hubHeight: inp.hubHeight }, controller: tb.controller, structure: str, rho: tb.rhoAir, rhoWater: inp.rhoW, dt: 0.05 },
        site: { hubHeight: inp.hubHeight, rotorDiameter: 2 * tb.R, waterDepth: inp.depth, Vref: 50, Iref: inp.Iref, ratedWind: tb.schedule.ratedWind || 11.4, cutIn, cutOut, weibullA: inp.weibullA, weibullK: inp.weibullK, seaState: tab, gamma: inp.gamma, macCamyFuchs: inp.mcf },
        cases: ['1.2'], nSeeds: 1, duration: inp.tdDur, tTransient: 20, binWidth: Math.max(inp.binW, 2), lifetimeYears: 1, availability: av, wohlerM: 4, rainflow: rainflowCount,
        snDamage: ({ ranges, counts }) => minerSum({ ranges: ranges.map((x) => x / sec.W), counts, sn }).damage,
      }), 'Time-domain fatigue (design load case 1.2)');
      const ch = atMud ? 'mudlineMoment' : 'towerBaseMoment';
      // The driver returns the damage of the channel moment divided by the section modulus of the assessed section.
      td = { damage: res.fatigue.damage[ch], del: res.fatigue.del[ch] / sec.W, notes: res.notes };
    }
    const prodOnly = states.filter((s) => s.mode === 'production').reduce((a, s) => a + s.dy, 0);
    const checks = [
      check('Fatigue damage over the service life with the design fatigue factor', ucLife, 1, '-', 'DNV-RP-C203 (2021) §2.2 (Miner sum); DNV-ST-0126 (2021) §4.11', 'Miner sum over the required life times the design fatigue factor. Above 1 a fatigue crack is predicted at the weld before the end of life with the required safety margin.'),
    ];
    const ok = ucLife <= 1;
    // Wall thickness for unit utilisation: damage ∝ stress^m with m between 3 and 5; nominal stress ∝ 1/t (thin wall) — use m = 4.
    const tNow = sec.t * 1e3, tReq = tNow * ucLife ** 0.25;
    const lenSec = atMud ? 12 : 10; // length of the thick can affected, m
    const dMass = (tube(sec.D, Math.min(tReq, 0.4 * sec.D * 1e3) / 1e3).A - sec.A) * steel.rho * lenSec / 1e3;
    const cats = states.map((s) => `${s.label} ${s.mode === 'production' ? '' : '(' + s.mode.split(' ')[0] + ')'}`.trim());
    const [rf, rS] = thin([rated.sp.f, Float64Array.from(rated.sp.SM, (v) => v * (mRatio / sec.W) ** 2 / 1e12)], 400);
    const sr = Array.from({ length: 40 }, (_, i) => 10 * 10 ** ((i * 1.6) / 39));
    const top = states.reduce((a, s) => (s.dy > a.dy ? s : a), states[0]);
    return {
      kpis: [
        { key: 'damageYr', label: 'Fatigue damage per year (characteristic)', unit: '1/yr', value: dYr, explain: 'Miner damage accumulated in an average year at the assessed weld, without the design fatigue factor. Its inverse is the characteristic fatigue life.' },
        { key: 'lifeYears', label: 'Fatigue life with the design fatigue factor', unit: 'years', value: Math.max(totalLife, 0), limit: inp.life, sense: 'min', explain: 'Years from new until the factored Miner sum reaches 1, allowing for the damage consumed before operation. It must exceed the required service life.' },
        { key: 'remaining', label: 'Remaining useful life', unit: 'years', value: remaining, explain: 'Years left from today at the present damage rate, including the design fatigue factor and the damage already consumed. The basis for a life-extension decision.' },
        { key: 'ucLife', label: 'Fatigue utilisation at the required life', unit: '-', value: ucLife, limit: 1, sense: 'max', explain: 'Design fatigue factor times the Miner sum at the end of the required service life.' },
        { key: 'dUsed', label: 'Damage consumed to date', unit: '-', value: dUsed, explain: 'Characteristic Miner sum already used: installation damage plus the years in service at the computed rate.' },
        { key: 'delStress', label: 'Damage-equivalent stress range (1e7 cycles, m = 4)', unit: 'MPa', value: delS / 1e6, explain: 'Constant-amplitude nominal stress range that, repeated ten million times, gives the lifetime damage of a slope-4 curve (narrow-band estimate). Useful for comparing designs; not used in the Miner sum.' },
        { key: 'prodShare', label: 'Share of damage from power production', unit: '%', value: dYr > 0 ? prodOnly / dYr * 100 : 0, explain: 'How much of the damage arises while the rotor is producing. A low share means idling and storm conditions without aerodynamic damping dominate — typical of large monopiles in exposed sites.' },
        { key: 'prodFrac', label: 'Time within the operating wind range', unit: '%', value: prodFrac * 100, explain: 'Probability that the hub wind speed lies between cut-in and cut-out for the given Weibull climate.' },
        { key: 'nRef', label: 'Cycles to failure at the reference stress range', unit: '-', value: cyclesToFailure({ curve, stressRange: inp.sRef * 1e6, thickness: 0.025, scf: 1 }), explain: 'Characteristic S–N life at the reference range for the basic curve (reference thickness, SCF 1): a readout to check the curve and environment selected.' },
        { key: 'f1', label: 'First natural frequency', unit: 'Hz', value: md.frequency, explain: 'Frequency of the mode whose resonance amplifies wave loading; the closer it is to the wave band, the larger the damage.' },
        ...(td ? [{ key: 'damageYrTime', label: 'Damage per year, time-domain rainflow (production)', unit: '1/yr', value: td.damage, explain: 'Miner damage per year from rainflow counting of coupled simulations, one short record per wind bin, production only. Agreement with the spectral production value within a factor of about 2 is normal for single short seeds.' }] : []),
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Annual damage by condition', categories: cats, series: [{ name: 'Damage per year', y: states.map((s) => s.dy) }], x: { label: 'Wind bin' }, y: { label: 'Miner damage per year', unit: '1/yr' },
          explain: 'Contribution of each wind bin and turbine state to the yearly damage. Look for the dominant bins: around rated wind (thrust fluctuation) or the idling and storm bins (wave resonance without aerodynamic damping).' },
        { type: 'line', title: 'Stress spectrum near rated wind', x: { label: 'Frequency', unit: 'Hz' }, y: { label: 'Nominal stress PSD', unit: 'MPa²/Hz' }, series: [{ name: `V = ${fmt(rated.U)} m/s`, x: rf, y: rS }],
          explain: 'Frequency content of the nominal bending stress in the production bin nearest rated wind. The area is the stress variance; the resonance peak shows the role of damping.' },
        { type: 'line', title: 'S–N curve and equivalent stress range', x: { label: 'log10 cycles to failure', unit: '-' }, y: { label: 'Stress range', unit: 'MPa' },
          series: [{ name: `Curve ${inp.curve}, ${inp.env}, t = ${fmt(tNow)} mm, SCF ${inp.scf}`, x: sr.map((s) => Math.log10(cyclesToFailure({ curve, stressRange: s * 1e6, thickness: sec.t, scf: inp.scf }))), y: sr }, { name: 'Equivalent range at 1e7 cycles', x: [7], y: [delS / 1e6] }],
          explain: 'Design S–N curve including thickness correction and stress concentration, with the damage-equivalent nominal stress range plotted at ten million cycles. The point must lie below the curve with margin for the fatigue factor.' },
      ],
      tables: [
        { title: 'Damage by wind bin and turbine state', columns: [{ key: 'bin', label: 'Wind bin' }, { key: 'mode', label: 'State' }, { key: 'p', label: 'Probability', unit: '%' }, { key: 'Hs', label: 'Hs', unit: 'm' }, { key: 'Tp', label: 'Tp', unit: 's' }, { key: 'sig', label: 'σ of nominal stress', unit: 'MPa' }, { key: 'z', label: 'Damping', unit: '%' }, { key: 'd', label: 'Damage per year' }, { key: 'sh', label: 'Share', unit: '%' }],
          rows: states.map((s) => ({ bin: s.label, mode: s.mode, p: round(s.prob * 100), Hs: round(s.Hs), Tp: round(s.Tp), sig: round(s.sigma / 1e6), z: round(s.zeta * 100), d: round(s.dy), sh: round(dYr > 0 ? s.dy / dYr * 100 : 0) })),
          explain: 'One row per lumped environmental state: probability, sea state, stress standard deviation, total damping and yearly damage. Use it to see which conditions to simulate in detail.' },
        ...(td ? [{ title: 'Frequency-domain against time-domain (production bins)', columns: [{ key: 'k', label: 'Quantity' }, { key: 's', label: 'Spectral (Dirlik)' }, { key: 't', label: 'Time-domain rainflow' }],
          rows: [{ k: 'Damage per year at full availability share', s: round(prodOnly), t: round(td.damage) }, { k: 'Ratio time / spectral', s: 1, t: round(td.damage / Math.max(prodOnly, 1e-30)) }],
          explain: 'Cross-check of the two methods on the production bins. The time-domain value comes from one short seed per bin and includes controller dynamics; the spectral value includes the full low-frequency turbulence.' }] : []),
      ],
      notes: [`Assessed: ${atMud ? 'monopile at the mudline' : 'tower base'}, D = ${fmt(sec.D)} m, t = ${fmt(tNow)} mm, curve ${inp.curve} (${inp.env}), SCF ${inp.scf}, DFF ${inp.dff}.`,
        `Largest contribution: ${top.label} ${top.mode} (${fmt(top.dy / Math.max(dYr, 1e-30) * 100)} % of the damage).`, `Characteristic life ${fmt(Math.min(1 / Math.max(dYr, 1e-30), 1e4))} years; with DFF ${inp.dff} and initial damage ${inp.d0}: ${fmt(Math.max(totalLife, 0))} years.`,
        ...(td ? td.notes : []), 'Wind and waves are taken as aligned and every cycle acts on the same weld point: directional spreading typically reduces monopile damage by 20–40 %.'],
      advice: [
        ok ? adv('engineering', sev(ucLife), `Fatigue utilisation is ${fmt(ucLife)} at ${inp.life} years (life ${fmt(Math.max(totalLife, 0))} years with DFF ${inp.dff}); ${fmt(remaining)} years remain after ${inp.yearsUsed} years in service.`,
          `Confirm with the supplier's load simulations (this estimate uses a generic rotor and a single-mode model) and the site scatter diagram, and check the governing ${top.mode} condition at ${top.label}; a wall of about ${fmt(tReq)} mm at this section would use the life fully.`,
          'DNV-RP-C203 (2021) §2; DNV-ST-0126 (2021) §4.11; DNV-ST-0262 for life extension',
          `Fatigue life scales with about the fourth power of wall thickness: going from ${fmt(tNow)} to ${fmt(tReq)} mm over a ${lenSec} m can changes steel by ${fmt(dMass)} t — ${steelCo2Text(dMass)} Weld grinding (curve C1) gains more life than thickness with no extra steel.`)
          : adv('engineering', 'fail', `Fatigue utilisation ${fmt(ucLife)} exceeds 1: life is ${fmt(Math.max(totalLife, 0))} years against ${inp.life} required; ${top.label} ${top.mode} contributes ${fmt(top.dy / dYr * 100)} %.`,
            `Increase the wall at this section to about ${fmt(tReq * 1.03)} mm, improve the weld class (grinding to C1), reduce the SCF by tighter fit-up tolerance, or add damping if idling states dominate.`,
            'DNV-RP-C203 (2021) §2.4, §3.3.7 and Appendix D (weld improvement); DNV-ST-0126 (2021) §4.11',
            `Thickening the ${lenSec} m can adds ${fmt(dMass)} t — ${steelCo2Text(dMass)} Weld improvement or a damper achieves the same life with negligible material.`),
        adv('economic', ok ? 'info' : 'warn', `At ${fmt(dYr)} damage per year the structure has ${fmt(Math.max(totalLife - inp.life, 0))} years of fatigue life beyond the ${inp.life}-year requirement.`,
          totalLife > inp.life + 5 ? `The reserve supports a life-extension case of roughly ${fmt(Math.min(totalLife - inp.life, 15))} years of additional production, subject to inspection of the governing weld; alternatively trim the wall towards ${fmt(tReq)} mm at design stage.` : 'Budget for strain-gauge monitoring of the governing section from year one: measured damage rates are usually below design and are the evidence that later unlocks a life extension.',
          'Life-extension economics per DNV-ST-0262 / DNV-SE-0263; risk-based inspection planning per DNV-RP-C210; LCOE reduction through extended operating life'),
      ],
      publish: { 'life.fatigueYears': Math.max(totalLife, 0) },
    };
  },
  verify: [
    { name: 'S–N curve D in air: N = 10^12.164 / Δσ³ at Δσ = 100 MPa', source: 'DNV-RP-C203 Table 2-1, curve D (log a1 = 12.164, m1 = 3, t ≤ 25 mm): N = 1.4588e12 / 1e6 = 1.459e6 cycles', inputs: { curve: 'D', env: 'air', sRef: 100 }, expect: [{ kpi: 'nRef', value: 1.4588e6, tol: 0.002 }] },
    { name: 'Time within the operating range for a Weibull climate: exp(−(Vin/A)^k) − exp(−(Vout/A)^k)', source: 'Weibull distribution, hand calculation: A = 10.5 m/s, k = 2.1, 3–25 m/s → exp(−0.07202) − exp(−6.1827) = 0.93051 − 0.00207 = 0.9284', inputs: { weibullA: 10.5, weibullK: 2.1 }, expect: [{ kpi: 'prodFrac', value: 92.84, tol: 0.001 }] },
  ],
};
