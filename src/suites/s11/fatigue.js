// Spectral fatigue of a hull detail: stress transfer function from the wave-induced hull-girder moment, from
// motion-induced inertia, or a direct stress RAO; short-term damage per sea state and long-term life over the
// scatter diagram with DNV S–N curves and a design fatigue factor.
import { pointMotion } from '../../engine/hydro/floater.js';
import { snCurve, snConstants, SN_CURVE_NAMES } from '../../engine/fatigue/sn.js';
import { longTermSpectralDamage } from '../../engine/fatigue/spectral.js';
import { fatigueLife } from '../../engine/fatigue/damage.js';
import { hullInputs, hullFrom, motionInputs, motionFrom, seaFrom, spectrum, omegaGrid, pointInputs,
  G, DEG, YEAR_S, LIMIT_HYDRO, need, run, round, rnd, fmt, adv, check, sev, steelCo2Text } from './_common.js';
import { scatterInput, scatterFrom } from './frequency-domain.js';

/** Midship vertical wave bending moment per unit wave amplitude of a heave-balanced box (N·m/m), with Smith factor e^(−kT). */
export function boxMomentRao(L, B, T, omega, headingDeg, rho) {
  const k = omega * omega / G; const ke = Math.max(Math.abs(k * Math.cos(headingDeg * DEG)), 1e-9); const u = ke * L / 2;
  const sinc = Math.sin(u) / u;
  return rho * G * B * Math.exp(-k * T) * Math.abs((Math.cos(u) - 1) / (ke * ke) + (L / 2) * Math.sin(u) / ke - sinc * L * L / 8);
}

export default {
  id: 's11-fatigue',
  suite: 11,
  title: 'Wave-induced and motion-induced fatigue',
  summary: 'Long-term spectral fatigue damage and life of a hull detail from the scatter diagram, with the stress coming from hull-girder wave bending, from motion-induced inertia loads or from a given stress transfer coefficient.',
  covers: [9],
  fidelity: 'reduced-order',
  standards: ['DNV-RP-C203 (2021) §2.4, §5 and App. F', 'DNV-RP-C206 (2012) (fatigue methodology of offshore ships)', 'DNV-OS-C101 (2021) Ch.2 Sec.5 (design fatigue factors)', 'IACS Rec. 56 (fatigue assessment of ship structures)'],
  models: ['Stress transfer function × wave spectrum → stress spectrum', 'Closed-form box hull-girder wave moment with Smith correction', 'Motion-induced inertia stress from point acceleration RAO', 'Dirlik, narrow-band, Wirsching–Light and Tovo–Benasciutti spectral damage', 'Bilinear S–N curves with thickness correction (DNV-RP-C203)', 'Palmgren–Miner accumulation over the scatter diagram'],
  limits: `Nominal or hot-spot stress is a linear function of wave amplitude from one heading (no spreading, no heading distribution). The hull-girder transfer function is a quasi-static box-girder moment amidships with the Smith factor; local pressure, torsion, springing and whipping stresses are not included. ${LIMIT_HYDRO} Low-cycle fatigue from loading and unloading must be added separately.`,
  publishes: ['life.fatigueYears', 'life.designYears'],
  searchTerms: 'spectral fatigue FPSO hull detail S-N curve scatter diagram',
  inputs: [
    { key: 'mode', label: 'Source of cyclic stress', type: 'select', default: 'hull-girder', group: 'Stress transfer function', role: 'bc',
      options: [{ value: 'hull-girder', label: 'Hull-girder wave bending (ship, barge)' }, { value: 'motion', label: 'Motion-induced inertia at a point' }, { value: 'direct', label: 'Constant stress per metre wave amplitude' }],
      help: 'Hull-girder: deck or bottom detail of a ship-shaped unit. Motion-induced: supports of topsides, flare tower, crane pedestal. Direct: a coefficient from a finite-element load transfer.' },
    { key: 'Zsec', label: 'Hull-girder section modulus at the detail', unit: 'm³', type: 'number', default: 58, min: 0.1, max: 1000, group: 'Stress transfer function', role: 'geometry', showIf: { key: 'mode', in: ['hull-girder'] }, design: true, help: 'Section modulus at the deck or bottom from the global strength module.' },
    { key: 'kAcc', label: 'Stress per unit vertical acceleration', unit: 'MPa/(m/s²)', type: 'number', default: 25, min: 0, max: 1000, group: 'Stress transfer function', role: 'load', showIf: { key: 'mode', in: ['motion'] }, help: 'Nominal stress at the detail caused by 1 m/s² of vertical acceleration of the supported mass (from a static load case).' },
    { key: 'kDirect', label: 'Stress per metre of wave amplitude', unit: 'MPa/m', type: 'number', default: 12, min: 0, max: 500, group: 'Stress transfer function', role: 'load', showIf: { key: 'mode', in: ['direct'] }, help: 'Frequency-independent stress amplitude per metre of wave amplitude.' },
    ...hullInputs({ type: 'ship', withMesh: false }),
    { key: 'headingDeg', label: 'Wave heading', unit: 'deg', type: 'number', default: 180, min: 0, max: 360, group: 'Sea states', role: 'environment', help: '180° = head sea. A weather-vaning unit sees mostly head seas; a spread-moored unit needs each heading weighted separately.' },
    { key: 'gamma', label: 'JONSWAP peak factor γ', unit: '-', type: 'number', default: 2.0, min: 1, max: 7, group: 'Sea states', role: 'environment', help: 'Peakedness used for every sea state of the scatter diagram.' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 350, min: 20, max: 4000, group: 'Sea states', role: 'environment', bind: 'site.depth', help: 'Still-water depth.' },
    scatterInput('Sea states', [{ Hs: 1, Tp: 7, prob: 0.1 }, { Hs: 1.5, Tp: 8, prob: 0.2 }, { Hs: 2, Tp: 9, prob: 0.25 }, { Hs: 2.5, Tp: 10, prob: 0.18 }, { Hs: 3, Tp: 11, prob: 0.12 }, { Hs: 3.5, Tp: 12, prob: 0.08 }, { Hs: 4, Tp: 12, prob: 0.04 }, { Hs: 5, Tp: 13, prob: 0.02 }, { Hs: 6, Tp: 14, prob: 0.01 }]),
    ...motionInputs(),
    ...pointInputs({ x: 60, y: 0, z: 45 }),
    { key: 'curve', label: 'S–N curve', type: 'select', default: 'D', options: SN_CURVE_NAMES.map((n) => ({ value: n, label: `DNV-RP-C203 curve ${n}` })), group: 'Detail and S–N data', role: 'material', help: 'Joint classification of the welded detail (D for butt welds and hot-spot stress, F–F3 for attachments and fillet welds, B1–C2 for base material).' },
    { key: 'env', label: 'Environment', type: 'select', default: 'air', options: [{ value: 'air', label: 'In air / coated' }, { value: 'seawater-cp', label: 'Sea water with cathodic protection' }, { value: 'free-corrosion', label: 'Free corrosion' }], group: 'Detail and S–N data', role: 'material', help: 'Selects the S–N table of DNV-RP-C203.' },
    { key: 'thick', label: 'Plate thickness at the detail', unit: 'mm', type: 'number', default: 22, min: 5, max: 150, group: 'Detail and S–N data', role: 'geometry', help: 'For the thickness correction above the 25 mm reference.' },
    { key: 'scf', label: 'Stress concentration factor', unit: '-', type: 'number', default: 1.0, min: 1, max: 10, group: 'Detail and S–N data', role: 'geometry', design: true, uncertain: { dist: 'lognormal', cov: 0.1 }, help: 'Hot-spot over nominal stress for the detail (bracket toe, cut-out, misalignment).' },
    { key: 'method', label: 'Spectral damage method', type: 'select', default: 'dirlik', options: [{ value: 'dirlik', label: 'Dirlik' }, { value: 'narrow-band', label: 'Narrow-band (Rayleigh)' }, { value: 'wirsching-light', label: 'Wirsching–Light' }, { value: 'tovo-benasciutti', label: 'Tovo–Benasciutti' }], group: 'Detail and S–N data', role: 'control', help: 'Narrow-band is conservative for broad spectra; Dirlik is the usual choice.' },
    { key: 'D0', label: 'Damage already accumulated', unit: '-', type: 'number', default: 0, min: 0, max: 1, group: 'Design life', role: 'ic', help: 'Miner sum consumed before this assessment (previous service, conversion, tow).' },
    { key: 'designLife', label: 'Required service life', unit: 'yr', type: 'number', default: 25, min: 1, max: 100, group: 'Design life', role: 'criteria', bind: 'life.designYears', help: 'Remaining service life to be demonstrated.' },
    { key: 'dff', label: 'Design fatigue factor', unit: '-', type: 'number', default: 2, min: 1, max: 10, group: 'Design life', role: 'criteria', help: 'DNV-OS-C101: 1 for inspectable non-critical details, 2–3 for details below water or not accessible, up to 10 for non-inspectable critical ones.' },
  ],
  mesh: { key: 'nFreq', label: 'Frequency points', levels: [50, 100, 200], qoi: 'damageYear' },
  compute(inp) {
    const om = omegaGrid(inp.nFreq); const f = om.map((w) => w / (2 * Math.PI));
    let H; let what;
    if (inp.mode === 'hull-girder') {
      need(['ship', 'barge'].includes(inp.hullType), 'Hull-girder fatigue applies to ship-shaped hulls and barges: choose one of those hull types, or use the motion-induced or direct stress option.');
      const hull = hullFrom(inp);
      // effective breadth = waterplane area / L, so the tapered ends are not counted as full-breadth buoyancy
      const Beff = hull.hs.waterplaneArea / hull.Lc;
      H = om.map((w) => boxMomentRao(hull.Lc, Beff, hull.T, w, inp.headingDeg, hull.rho) / inp.Zsec);
      what = `hull-girder bending, Z = ${fmt(inp.Zsec)} m³`;
    } else if (inp.mode === 'motion') {
      const hull = hullFrom(inp); const mot = motionFrom(inp, hull);
      // drag linearised for a moderate fatigue sea (the cells that dominate damage)
      const r = mot.rao(om, inp.headingDeg, { Hs: 3, Tp: 10, gamma: inp.gamma });
      const pm = pointMotion({ rao: r, point: [inp.px, inp.py, inp.pz - hull.T] });
      H = Array.from(pm.accAmp[2], (a) => a * inp.kAcc * 1e6);
      what = `motion-induced inertia, ${fmt(inp.kAcc)} MPa per m/s²`;
    } else { H = om.map(() => inp.kDirect * 1e6); what = `direct coefficient ${fmt(inp.kDirect)} MPa/m`; }
    const curve = run(() => snCurve(inp.curve, inp.env), 'S–N curve');
    const sn = snConstants(curve, { thickness: inp.thick / 1e3, scf: inp.scf });
    const sc = scatterFrom(inp);
    const states = sc.map((s) => { const S = spectrum(om, { Hs: s.Hs, Tp: s.Tp, gamma: inp.gamma }); return { probability: s.prob, f, S: Array.from(S, (v, i) => H[i] * H[i] * v * 2 * Math.PI), label: `Hs ${s.Hs} / Tp ${s.Tp}` }; });
    const lt = run(() => longTermSpectralDamage({ states, sn, duration: YEAR_S, method: inp.method }), 'Spectral fatigue');
    const dYear = lt.damage;
    need(Number.isFinite(dYear) && dYear >= 0, 'The fatigue damage could not be evaluated: check the stress transfer coefficient and the S–N data.');
    const fl = fatigueLife({ damageRate: Math.max(dYear, 1e-30), designLife: inp.designLife, dff: inp.dff, accumulatedDamage: inp.D0 });
    const cap = (v) => Math.min(v, 1e5);
    const contrib = lt.contributions.map((c, i) => ({ ...c, Hs: sc[i].Hs, Tp: sc[i].Tp, share: dYear > 0 ? c.damage / dYear : 0 }));
    const top = contrib.reduce((a, c) => (c.damage > a.damage ? c : a), contrib[0]);
    const sigMax = Math.max(...contrib.map((c) => c.sigma));
    const util = fl.utilisation ?? inp.dff * (inp.D0 + dYear * inp.designLife);
    const chk = [check('Fatigue damage with design fatigue factor', util, 1, '-', 'DNV-RP-C203 (2021) §2.2; DNV-OS-C101 (2021) Ch.2 Sec.5 (DFF)',
      'Guards against fatigue cracking within the service life: the Miner sum times the design fatigue factor must not exceed 1.0, otherwise a through-thickness crack is expected before the end of life at the target reliability.')];
    // stress reduction needed: damage ∝ S^m with the governing slope (m between 3 and 5)
    const mEff = sn.m2 && sigMax * 2.83 < sn.sKnee ? sn.m2 : sn.m1;
    const scfAllow = inp.scf * Math.min((1 / Math.max(util, 1e-12)) ** (1 / mEff), 50);
    const scfReq = Math.min(scfAllow, inp.scf);
    const hsVals = [...new Set(sc.map((s) => s.Hs))].sort((a, b) => a - b); const tpVals = [...new Set(sc.map((s) => s.Tp))].sort((a, b) => a - b);
    const z = hsVals.map((Hh) => tpVals.map((Tp) => { const c = contrib.find((s) => s.Hs === Hh && s.Tp === Tp); return c ? round(100 * c.share, 3) : 0; }));
    const sr = [20, 30, 50, 80, 120, 200, 300].map((s) => s * 1e6);
    const nOf = (s) => (s >= sn.sKnee || !sn.m2 ? sn.a1 * s ** -sn.m1 : sn.a2 * s ** -sn.m2);
    return {
      kpis: [
        { key: 'damageYear', label: 'Fatigue damage per year', value: round(dYear, 5), unit: '1/yr', explain: 'Characteristic Miner sum accumulated in one year over the scatter diagram, without the design fatigue factor.' },
        { key: 'life', label: 'Characteristic fatigue life', value: round(cap(fl.characteristicLife), 5), unit: 'yr', explain: 'Years to a Miner sum of 1.0 on the design S–N curve (mean minus two standard deviations), from new.' },
        { key: 'lifeDesign', label: 'Design fatigue life (with DFF)', value: round(cap(fl.remainingLife), 5), unit: 'yr', limit: inp.designLife, sense: 'min', explain: 'Remaining life after dividing by the design fatigue factor and deducting the damage already consumed. It must cover the required service life.' },
        { key: 'utilisation', label: 'Fatigue utilisation', value: round(util, 5), unit: '-', limit: 1, sense: 'max', explain: 'DFF × (consumed damage + annual damage × service life). Above 1.0 the detail does not reach its life at the required reliability.' },
        { key: 'sigmaMax', label: 'Largest sea-state stress standard deviation', value: round(sigMax / 1e6, 4), unit: 'MPa', explain: 'Standard deviation of nominal stress in the most severe sea state; the significant stress range is four times this (before SCF).' },
        { key: 'cycles', label: 'Stress cycles per year', value: round(lt.cycles, 4), unit: '-', explain: 'Number of stress cycles in a year over all sea states; about 4–5 million for wave loading.' },
        { key: 'topShare', label: 'Share of damage from the worst sea state', value: round(100 * top.share, 4), unit: '%', explain: 'How concentrated the damage is. Fatigue is normally driven by frequent moderate seas, not by the storm.' },
        { key: 'scfAllow', label: 'Stress concentration factor that meets the life', value: round(scfAllow, 4), unit: '-', explain: 'Largest SCF (or equivalent stress scaling) at which utilisation is 1.0, estimated with the governing S–N slope: the target for detail improvement or the reserve available.' },
      ],
      checks: chk,
      plots: [
        { type: 'line', title: 'Stress transfer function', x: { label: 'Wave period', unit: 's' }, y: { label: 'Stress amplitude per metre wave amplitude', unit: 'MPa/m' }, series: [{ name: 'Nominal stress RAO', x: rnd(om.map((w) => 2 * Math.PI / w)), y: rnd(H.map((v) => v / 1e6)) }],
          explain: 'Nominal stress at the detail per metre of wave amplitude. For hull-girder bending it peaks where the wavelength is close to the hull length; sea states with peak periods there do the most damage.' },
        { type: 'heatmap', title: 'Damage share over the scatter diagram', xs: tpVals, ys: hsVals, z, zLabel: 'Share of annual damage', zUnit: '%', x: { label: 'Tp', unit: 's' }, y: { label: 'Hs', unit: 'm' },
          explain: 'Percentage of the annual fatigue damage from each sea state. The hot cells are the sea states to get right in the metocean data and the transfer function.' },
        { type: 'line', title: 'S–N curve with sea-state stress ranges', x: { label: 'log10 cycles to failure', unit: '-' }, y: { label: 'Stress range', unit: 'MPa' },
          series: [{ name: `Curve ${inp.curve} with SCF and thickness`, x: rnd(sr.map((s) => Math.log10(nOf(s)))), y: sr.map((s) => s / 1e6) },
            { name: 'Significant stress range per sea state (at annual cycles)', x: rnd(contrib.map((c) => Math.log10(Math.max(c.nu0 * c.probability * YEAR_S, 1)))), y: rnd(contrib.map((c) => 4 * c.sigma / 1e6)) }],
          explain: 'The design S–N curve (nominal stress, SCF folded in) with one point per sea state at its significant stress range and yearly cycle count. Points approaching the curve are the damaging sea states.' },
      ],
      tables: [
        { title: 'Damage by sea state', columns: [{ key: 'Hs', label: 'Hs', unit: 'm' }, { key: 'Tp', label: 'Tp', unit: 's' }, { key: 'p', label: 'Probability', unit: '%' }, { key: 'sig', label: 'Stress std dev', unit: 'MPa' }, { key: 'tz', label: 'Stress Tz', unit: 's' }, { key: 'd', label: 'Damage per year' }, { key: 'sh', label: 'Share', unit: '%' }],
          rows: contrib.map((c) => ({ Hs: c.Hs, Tp: c.Tp, p: round(100 * c.probability, 4), sig: round(c.sigma / 1e6, 4), tz: round(c.nu0 > 0 ? 1 / c.nu0 : 0, 4), d: round(c.damage, 4), sh: round(100 * c.share, 4) })),
          explain: 'One row per sea state: nominal stress standard deviation, stress zero-crossing period and its contribution to the annual damage.' },
      ],
      notes: [
        `Stress from ${what}; curve ${inp.curve} (${inp.env}), SCF ${fmt(inp.scf)}, thickness ${fmt(inp.thick)} mm, method ${inp.method}.`,
        `Annual damage ${fmt(dYear)}: characteristic life ${fmt(cap(fl.characteristicLife))} yr, ${fmt(cap(fl.designLife))} yr with DFF ${fmt(inp.dff)}. Worst sea state ${top.label} gives ${fmt(100 * top.share)} % of the damage.`,
      ],
      advice: [
        adv('engineering', sev(util), `Fatigue utilisation is ${fmt(util)} for ${fmt(inp.designLife)} years with DFF ${fmt(inp.dff)}; design life ${fmt(cap(fl.remainingLife))} yr.`,
          util > 1 ? `Reduce the hot-spot stress by ${fmt(100 * (1 - (1 / util) ** (1 / mEff)))} % (SCF from ${fmt(inp.scf)} to about ${fmt(scfReq)}): soft-toe brackets, ground weld toes (one S–N class better when documented), thicker insert plate, or a larger section modulus. Halving stress multiplies life by ${fmt(2 ** mEff)}.`
            : `Margin corresponds to a ${fmt(100 * ((1 / Math.max(util, 1e-9)) ** (1 / mEff) - 1))} % higher stress. Confirm the SCF by hot-spot finite-element analysis and include every heading with its probability for a spread-moored unit.`,
          'DNV-RP-C203 (2021) §2.4 (S–N curves), §5 (simplified and spectral fatigue), §7 (weld improvement); DNV-RP-C206 (2012).'),
        adv('economic', util > 1 ? 'warn' : 'info', `A detail with ${fmt(cap(fl.remainingLife))} years of design life against ${fmt(inp.designLife)} required ${util > 1 ? 'will need in-service inspection and probably repair offshore' : 'needs no fatigue-driven intervention in the design life'}.`,
          util > 1 ? 'Fix it in the yard: a bracket re-shape or weld-toe grinding at fabrication costs hours; an offshore crack repair costs tank entry, hot work, and often production deferral.' : 'Use the margin to set a longer inspection interval for this detail in the risk-based inspection plan rather than removing steel.',
          'Risk-based inspection planning per DNV-RP-C210 (2019) and DNV-RP-G101; whole-life costing per ISO 15686-5.',
          `Life extension by detail improvement uses almost no material, whereas replacing a cracked 2 t insert offshore repeats its footprint and adds vessel and hot-work emissions: ${steelCo2Text(2)}`),
      ],
      publish: { 'life.fatigueYears': cap(fl.designLife), 'life.designYears': inp.designLife },
    };
  },
  verify: [
    { // One sea state Hs = 4 m, Tp = 10 s (γ = 1), stress 10 MPa per metre amplitude → σ = 10·Hs/4 = 10 MPa;
      // Tz ≈ Tp·√((5 + γ)/(11 + γ)) = 7.071 s (DNV-RP-C205 §3.5.5.6) → ν0 = 0.14142 Hz. Curve D, free corrosion: m = 3, log a = 11.687.
      // Narrow band: D = ν0·T·(2√2·σ)^m·Γ(1 + m/2)/a = 0.14142·31 557 600·(28.284)³·1.3293/4.864e11 = 0.2760 per year → life 3.62 yr.
      name: 'Narrow-band closed form, single-slope S–N curve, one sea state',
      source: 'DNV-RP-C203 (2021) App. F / Rayleigh damage formula; curve D free corrosion (Table 2-4) — hand calculation',
      inputs: { mode: 'direct', kDirect: 10, scatter: [{ Hs: 4, Tp: 10, prob: 1 }], gamma: 1, curve: 'D', env: 'free-corrosion', thick: 20, scf: 1, method: 'narrow-band', nFreq: 200 },
      expect: [{ kpi: 'damageYear', value: 0.276, tol: 0.05 }, { kpi: 'sigmaMax', value: 10, tol: 0.01 }] },
    { // Flat transfer function of 74.39 MPa per metre amplitude in Hs = 2 m: σ = 74.39·2/4 = 37.195 MPa (m0 of the wave = Hs²/16).
      name: 'Stress standard deviation = transfer coefficient × Hs/4',
      source: 'Definition of significant wave height, Hs = 4·√m0 (DNV-RP-C205 §3.5)',
      inputs: { mode: 'direct', kDirect: 74.39, scatter: [{ Hs: 2, Tp: 9, prob: 1 }], nFreq: 200 },
      expect: [{ kpi: 'sigmaMax', value: 37.195, tol: 0.01 }] },
  ],
};
