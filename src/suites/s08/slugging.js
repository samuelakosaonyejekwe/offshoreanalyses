// Slugging: hydrodynamic slug characteristics from the unit-cell model (frequency, length, velocity, volume, arrival
// history, surge volume), stratified-flow stability (Kelvin–Helmholtz), riser-induced severe slugging criteria and
// cycle, and slug loads on a bend.
import { flowlineRiserProfile, trackSlugs } from '../../engine/flow/transient.js';
import { marchPipeline } from '../../engine/flow/multiphase.js';
import { flowRegime, taitelDukler } from '../../engine/flow/flowregime.js';
import { slugUnitCell, slugFrequencyGregoryScott, slugFrequencyZabaras, meanSlugLength, slugLengthDistribution, slugTrain, kelvinHelmholtzStability,
  boeCriterion, potsNumber, taitelStability, severeSluggingCycle, severeSluggingMitigation, surgeVolume, slugBodyHoldupGregory, slugTranslationalVelocity } from '../../engine/flow/slug.js';
import { slugBendLoad, dynamicLoadFactor, allowableSlugVelocity, slugFatigueCycles } from '../../engine/flow/slugloads.js';
import { pipeInputs, pipeFrom, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, sev, worst, linspace, thin, clamp,
  steelCo2, LIMIT_NO_CFD, ECON_BASIS, CO2_STEEL, BAR, KELVIN, DEG, G } from './_common.js';
import { fluidInputs, fluidFrom } from './_flow.js';

const G_S = 'Flowline and riser';
const G_H = 'Slug model';
const G_L = 'Slug loads and receiving facility';

export default {
  id: 's08-slugging',
  suite: 8,
  title: 'Slugging and slug loads',
  summary: 'Predicts whether the line slugs, how often, how long and how fast the slugs are, the surge the separator must absorb, and whether the riser will severe-slug.',
  covers: [15],
  fidelity: 'reduced-order',
  standards: ['Taitel & Barnea (1990) slug unit-cell model', 'Gregory & Scott (1969); Zabaras (2000) slug frequency', 'Scott, Shoham & Brill (1989); Brill et al. (1981) slug length', 'Bøe (1981); Taitel (1986); Pots et al. (1987) severe slugging', 'DNV-RP-D101 (2017) §3.12 (slug loads)'],
  models: ['Mechanistic slug-flow unit cell', 'Slug frequency correlations', 'Lognormal slug-length distribution and slug train (slug arrival profile)', 'Lagrangian slug tracking over a uniform film (growth, decay and merging)', 'Kelvin–Helmholtz instability (inviscid and viscous) — slug initiation', 'Terrain-induced severe slugging criteria and cycle', 'Momentum force of a slug on a bend with a single-degree-of-freedom dynamic load factor'],
  limits: `Fully developed hydrodynamic slug flow from a steady unit cell with empirical frequency and length correlations calibrated on small-diameter or specific field data (scatter of a factor two or more is normal); slug lengths are sampled from a lognormal distribution with a fixed seed. Severe slugging is assessed by quasi-static criteria and a cycle estimate, not a dynamic simulation (see the transient module). The slug tracker uses a prescribed uniform film, so liquid is not conserved between slugs. Slug tracking with film dynamics and gas compressibility, and two-fluid six/seven-equation models are outside this module. ${LIMIT_NO_CFD}`,
  publishes: ['pipeline.slugFrequency', 'pipeline.slugVolume', 'pipeline.slugVelocity'],
  searchTerms: 'slug flow frequency slug length severe slugging riser slug catcher surge volume',
  inputs: [
    ...pipeInputs({ omit: ['tCoat', 'rhoCoat', 'tConc', 'concrete', 'rhoCont', 'rhoW', 'wsGiven'] }),
    num('rough', 'Wall roughness', 'mm', 0.045, 0.001, 3, 'Pipe', 'geometry', 'Absolute roughness of the bore.'),
    num('Lflow', 'Flowline length', 'm', 1500, 50, 200000, G_S, 'geometry', 'Length of the flowline upstream of the riser base.'),
    num('angle', 'Flowline inclination (negative = downhill to the riser)', 'deg', -1, -10, 10, G_S, 'geometry', 'A flowline falling towards the riser base lets liquid block the gas and is the classic severe-slugging geometry.'),
    num('hRiser', 'Riser height', 'm', 120, 0, 3000, G_S, 'geometry', 'Vertical height from the riser base to the separator; zero for a line arriving at seabed level.', { bind: 'site.depth' }),
    num('pSep', 'Separator (outlet) pressure', 'bar', 20, 1.5, 300, G_S, 'bc', 'Absolute pressure at the riser top. With the rate it fixes the hydraulic state: the inlet pressure is solved.', { design: true }),
    num('Tflow', 'Flowing temperature', '°C', 40, 0, 150, G_S, 'bc', 'Mean fluid temperature in the flowline (isothermal assessment).'),
    ...fluidInputs({ qo: 3000, GOR: 160, WC: 20 }),
    num('lenFactor', 'Mean slug length in diameters (0 = Scott correlation)', '-', 32, 0, 600, G_H, 'operation', 'About 32 D for fully developed horizontal slugs; zero uses the Scott–Shoham–Brill large-diameter correlation.'),
    num('sigmaLn', 'Log-standard deviation of slug length', '-', 0.5, 0.1, 1.2, G_H, 'operation', 'Brill et al. found about 0.5 at Prudhoe Bay; the 0.1 % exceedance length is the design maximum.'),
    num('duration', 'Length of the slug arrival history', 's', 600, 60, 7200, G_H, 'control', 'Duration of the synthetic slug train at the outlet.'),
    { key: 'seed', label: 'Random seed of the slug train', type: 'int', default: 1, min: 1, max: 99999, group: G_H, role: 'control', help: 'Fixed seed so the slug train is reproducible.' },
    num('drain', 'Liquid drain rate / mean liquid rate', '-', 1.2, 0.5, 5, G_L, 'operation', 'Liquid handling capacity of the separator relative to the mean liquid arrival rate.'),
    num('vSurge', 'Available surge volume', 'm³', 25, 0.1, 5000, G_L, 'criteria', 'Liquid volume between normal and high level in the separator or slug catcher.'),
    num('bend', 'Bend angle', 'deg', 90, 10, 180, G_L, 'geometry', 'Angle the flow turns through at the bend checked (riser base, topsides).'),
    num('fn', 'Natural frequency of the bend support', 'Hz', 4, 0.2, 50, G_L, 'bc', 'First natural frequency of the piping span or support carrying the bend.'),
    num('Fallow', 'Allowable dynamic bend force', 'kN', 20, 0.5, 5000, G_L, 'criteria', 'Force range the support and piping can take, from the piping stress analysis.'),
  ],
  compute(inp, ctx) {
    const p = pipeFrom({ ...inp, rhoCont: 0 }, ctx);
    const D = p.Di, A = Math.PI * D * D / 4, rough = inp.rough / 1e3;
    const fluid = fluidFrom(inp);
    const T = inp.Tflow + KELVIN, theta = inp.angle * DEG;
    const route = inp.hRiser > 1 ? flowlineRiserProfile({ flowlineLength: inp.Lflow, flowlineAngleDeg: inp.angle, riserHeight: inp.hRiser }) : [{ x: 0, z: 0 }, { x: inp.Lflow * Math.cos(theta), z: inp.Lflow * Math.sin(theta) }];
    const Ltot = inp.Lflow + inp.hRiser;
    const march = run(() => marchPipeline({ profile: route, D, roughness: rough, fluidAt: fluid.at, temperature: T, pOut: inp.pSep * BAR, massRate: fluid.massRate, method: 'beggs-brill', maxStep: Ltot / 30 }), 'Steady hydraulic state');
    const stF = march.stations.filter((s) => s.s <= inp.Lflow + 1e-6);
    const base = stF[stF.length - 1], mid = stF[Math.floor(stF.length / 2)];
    const s = fluid.state(mid.P, T, D), f = fluid.at(mid.P, T);
    need(s.vsl > 1e-4 && s.vsg > 1e-4, `Slug flow needs both gas and liquid: at flowline conditions (${fmt(mid.P / BAR)} bar) the stream is single-phase (v_sl ${fmt(s.vsl)} m/s, v_sg ${fmt(s.vsg)} m/s). Lower the pressure or raise the gas–oil ratio.`);
    const vm = s.vsl + s.vsg;
    const args = { vsl: s.vsl, vsg: s.vsg, D, theta, rhoL: f.rhoL, rhoG: f.rhoG, muL: f.muL, muG: f.muG };
    const reg = run(() => flowRegime({ ...args, sigma: f.sigma }), 'Flow regime');
    const Lmean = inp.lenFactor > 0 ? meanSlugLength({ D, method: 'diameters', factor: inp.lenFactor }) : meanSlugLength({ D, method: 'scott' });
    const fGS = slugFrequencyGregoryScott({ vsl: s.vsl, vm, D }), fZ = slugFrequencyZabaras({ vsl: s.vsl, vm, D, theta: Math.max(theta, 0) });
    let uc = run(() => slugUnitCell({ ...args, roughness: rough, slugLength: Lmean }), 'Slug unit cell');
    const cellOk = uc.feasible;
    if (!uc.feasible) {
      // the flowline is stratified at this angle: characterise the slugs it would carry if they form (turndown, start-up, upward
      // sections) with the horizontal unit cell, or failing that with the stratified film and an empirical frequency
      const ucH = slugUnitCell({ ...args, theta: 0, roughness: rough, slugLength: Lmean });
      if (ucH.feasible) uc = ucH;
      else {
        const HLS0 = slugBodyHoldupGregory(vm), vt0 = slugTranslationalVelocity({ vm, D, theta: 0 }).vt;
        const HLF = clamp(Math.min(mid.holdup, 0.6 * HLS0), 0.02, HLS0 - 0.05), vLF = s.vsl / HLF * 0.5;
        const beta = clamp((s.vsl - vLF * HLF) / Math.max(vm * HLS0 - vLF * HLF, 1e-9), 0.02, 0.9);
        uc = { ...uc, slugHoldup: HLS0, vt: vt0, filmHoldup: HLF, vLF, vGF: (vm - vLF * HLF) / (1 - HLF), slugFraction: beta, frequency: fZ };
      }
    }
    const dist = slugLengthDistribution({ meanLength: Lmean, sigmaLn: inp.sigmaLn });
    const HLS = uc.slugHoldup, vt = uc.vt;
    const volMean = A * Lmean * HLS, volMax = A * dist.maxLength * HLS;
    const slugFlow = cellOk && /intermittent|slug/.test(reg.regime);
    const train = run(() => slugTrain({ D, vm, vt, slugHoldup: HLS, filmHoldup: clamp(uc.filmHoldup, 0.01, HLS - 0.01), vLF: uc.vLF, vGF: uc.vGF, slugFraction: clamp(uc.slugFraction, 0.02, 1), meanSlugLength: Lmean, sigmaLn: inp.sigmaLn, duration: inp.duration, nSamples: 600, seed: inp.seed }), 'Slug train');
    // Lagrangian tracking of slugs initiated at the inlet at the unit-cell frequency: do they grow, decay or merge on the way?
    let trk = null;
    try {
      const fIn = Math.max(uc.frequency, 1e-4), HLFt = clamp(uc.filmHoldup, 0.01, HLS - 0.02);
      trk = trackSlugs({ length: inp.Lflow, D, vm, vt, slugHoldup: HLS, filmHoldup: HLFt, vLF: uc.vLF, inletPeriod: 1 / fIn, inletSlugLength: 0.5 * Lmean, wakeA: 0, tEnd: Math.min(Math.max(3 * inp.Lflow / vt, 60), 7200), dt: Math.max(0.02 * D / vt, Math.min(3 * inp.Lflow / vt, 7200) / 4000) });
    } catch { trk = null; }
    const qLmean = s.qLiquid;
    const surge = surgeVolume({ slugVolume: volMax, inflowRate: A * vm * HLS, drainRate: inp.drain * qLmean });
    // stratified-flow stability at the equilibrium level (slug initiation)
    let kh = null;
    try { const td = taitelDukler(args); kh = kelvinHelmholtzStability({ h: td.hL, ...args, roughness: rough }); } catch { kh = null; }
    // severe slugging at the riser base
    const hasRiser = inp.hRiser > 1;
    const sb = fluid.state(base.P, T, D), fb = fluid.at(base.P, T);
    const alphaF = clamp(1 - mid.holdup, 0.05, 0.98);
    const boe = hasRiser ? boeCriterion({ vsl: sb.vsl, vsg: Math.max(sb.vsg, 1e-6), P: base.P, rhoL: fb.rhoL, alphaPipe: alphaF, L: inp.Lflow }) : { piSS: 99, severeSlugging: false, vslCritical: 0 };
    const gasRT = base.P / fb.rhoG;
    const pots = hasRiser && sb.mGas > 0 ? potsNumber({ massRateG: sb.mGas, massRateL: sb.mTotal - sb.mGas, T, M: 8.314462618 * T / gasRT, Z: 1, alphaPipe: alphaF, L: inp.Lflow }) : { piSS: 99, severeSlugging: false };
    const tai = hasRiser ? taitelStability({ Psep: inp.pSep * BAR, rhoL: fb.rhoL, alphaPipe: alphaF, L: inp.Lflow, riserHeight: inp.hRiser }) : { stable: true, PsepMin: 0 };
    const downhill = inp.angle < 0;
    const severe = hasRiser && downhill && boe.severeSlugging && !tai.stable;
    const cyc = hasRiser && sb.mGas > 0 ? run(() => severeSluggingCycle({ vsl: sb.vsl, massRateG: sb.mGas, gasRT, D, L: inp.Lflow, riserHeight: inp.hRiser, Psep: inp.pSep * BAR, rhoL: fb.rhoL, alphaPipe: alphaF }), 'Severe-slugging cycle') : null;
    const mit = hasRiser ? severeSluggingMitigation({ vsl: sb.vsl, vsg: Math.max(sb.vsg, 1e-6), Psep: inp.pSep * BAR, rhoL: fb.rhoL, alphaPipe: alphaF, L: inp.Lflow, riserHeight: inp.hRiser, D }) : { chokeDpForStability: 0 };
    // slug load on a bend
    const rhoSlug = HLS * f.rhoL + (1 - HLS) * f.rhoG, rhoFilm = uc.filmHoldup * f.rhoL + (1 - uc.filmHoldup) * f.rhoG;
    const load = slugBendLoad({ rhoSlug, rhoFilm, v: vm, D, bendAngle: inp.bend * DEG });
    const dlf = dynamicLoadFactor({ naturalFrequency: inp.fn, duration: Lmean / vt, riseTime: 1.5 * D * (inp.bend * DEG) / vt, damping: 0.02, stepsPerPeriod: 60 }).dlf;
    const Fdyn = dlf * load.forceRange;
    const vAllow = allowableSlugVelocity({ allowableForce: inp.Fallow * 1e3, rhoSlug, rhoFilm, D, bendAngle: inp.bend * DEG, dlf });
    const cycles = slugFatigueCycles({ forceRange: load.forceRange, stressPerForce: 1, slugFrequency: uc.frequency, dlf });
    const governingVol = severe && cyc ? Math.max(volMax, cyc.slugVolume) : volMax;
    const surgeGov = severe && cyc ? Math.max(surge.surgeVolume, surgeVolume({ slugVolume: cyc.slugVolume, inflowRate: Math.max(cyc.peakLiquidRate, inp.drain * qLmean * 1.01), drainRate: inp.drain * qLmean }).surgeVolume) : surge.surgeVolume;
    const checks = [
      chk('Slug surge volume against the separator capacity', surgeGov, inp.vSurge, 'm³', 'Liquid volume balance (GPSA Engineering Data Book §17 practice)',
        'Guards against the largest slug overfilling the separator between normal and high-high level, which trips production.'),
      chk('Dynamic slug force on the bend', Fdyn, inp.Fallow * 1e3, 'N', 'Momentum theorem with SDOF dynamic load factor; DNV-RP-D101 (2017) §3.12',
        'Guards against the force step as a slug front passes a bend overloading supports and fatiguing the piping.'),
      chk('Severe slugging (Bøe criterion with Taitel stability)', severe ? 1 / Math.max(boe.piSS, 1e-6) : Math.min(1 / Math.max(boe.piSS, 1e-6), tai.stable || !downhill ? 0.99 : 1 / Math.max(boe.piSS, 1e-6)), 1, '-', 'Bøe (1981); Taitel (1986) Eq. 6',
        'Flags riser-induced severe slugging: liquid blocks the riser base and is expelled in long slugs with large pressure and rate swings. Above 1.0 with an unstable Taitel criterion, mitigation is needed.'),
    ];
    const gov = worst(checks);
    const [tt, ql] = thin([train.t, train.qL], 600);
    // frequency against rate (turndown)
    const frs = linspace(0.3, 1.8, 10);
    const fSweep = frs.map((k) => { const fl = fluid.scaled(k).state(mid.P, T, D); return { k, vm: fl.vsl + fl.vsg, f: fl.vsl > 0 ? slugFrequencyZabaras({ vsl: fl.vsl, vm: fl.vsl + fl.vsg, D, theta: Math.max(theta, 0) }) : 0 }; });
    const vg = linspace(0.05, Math.max(3 * sb.vsg, 1), 20);
    const lens = train.slugs.map((q) => q.length);
    const steelSC = 2 * Math.max(surgeGov, 1) * 0.35; // rough vessel steel: 0.35 t per m³ of gross volume, twice the surge volume
    return {
      kpis: [
        kpi('slugFreq', 'Slug frequency (unit cell)', uc.frequency, '1/s', 'Slugs per second passing a point, from the translational velocity over the slug-unit length of the unit-cell model.'),
        kpi('slugFreqGS', 'Slug frequency (Gregory–Scott)', fGS, '1/s', 'Empirical frequency correlation for horizontal flow; compare with the unit-cell value to see the model spread.'),
        kpi('slugFreqZ', 'Slug frequency (Zabaras)', fZ, '1/s', 'Gregory–Scott with an inclination correction for upward flow.'),
        kpi('slugLen', 'Mean slug length', Lmean, 'm', 'Mean length of the liquid slug body.'),
        kpi('slugLenMax', 'Design maximum slug length (0.1 % exceedance)', dist.maxLength, 'm', 'Length exceeded by one slug in a thousand in the lognormal distribution; used for sizing.'),
        kpi('slugVel', 'Slug translational velocity', vt, 'm/s', 'Speed of the slug front and tail (Bendiksen): faster than the mixture, it sets arrival times and the rise time of loads.'),
        kpi('slugHoldup', 'Liquid holdup in the slug body', HLS, '-', 'Liquid fraction of the aerated slug body (Gregory et al.).'),
        kpi('slugVol', 'Liquid volume of the design slug', governingVol, 'm³', 'Liquid delivered by the largest slug: the hydrodynamic design slug, or the severe-slugging blow-out volume when that governs.'),
        kpi('surge', 'Surge volume required', surgeGov, 'm³', 'Volume the separator must absorb while the design slug arrives faster than it is drained.', inp.vSurge, 'max'),
        kpi('trackGrowth', 'Slug front minus tail velocity (tracking)', trk ? trk.frontVelocity - vt : 0, 'm/s', 'From the slug tracker: positive means slugs pick up more film liquid at the front than they shed at the tail and grow along the line; negative means they decay before the outlet.'),
        kpi('trackArrivals', 'Slugs reaching the outlet in the tracked period', trk ? trk.arrivals.length : 0, '-', 'Number of tracked slugs that survive to the outlet; fewer than were launched means decay or merging on the way.'),
        kpi('peakQL', 'Peak liquid arrival rate', train.peakQL * 3600, 'm³/h', 'Largest liquid rate at the outlet while a slug body passes, against the mean rate below.'),
        kpi('meanQL', 'Mean liquid rate', qLmean * 3600, 'm³/h', 'Time-averaged in-situ liquid rate at flowline conditions.'),
        kpi('vm', 'Mixture velocity in the flowline', vm, 'm/s', 'Sum of superficial gas and liquid velocities at mid-flowline conditions.'),
        kpi('piSS', 'Severe-slugging number (Bøe)', Math.min(boe.piSS, 99), '-', 'Ratio of flowline gas pressure build-up to riser head build-up. Below 1.0 severe slugging is possible if the flowline is stratified and the Taitel criterion is unstable.'),
        kpi('piPots', 'Severe-slugging number (Pots)', Math.min(pots.piSS, 99), '-', 'Equivalent mass-rate form of the same ratio; severe slugging for values below 1.'),
        kpi('chokeDp', 'Choke pressure drop for stability', mit.chokeDpForStability / BAR, 'bar', 'Topside choking that just satisfies the Taitel stability criterion; zero when already stable.'),
        kpi('ssPeriod', 'Severe-slugging cycle period', cyc && cyc.severeSlugging ? cyc.period : 0, 's', 'Estimated time for one build-up, production and blow-out cycle; zero when no full blockage cycle forms.'),
        kpi('Fbend', 'Slug force range on the bend', load.forceRange / 1e3, 'kN', 'Step in the resultant bend force between film and slug body passing, before dynamic amplification.'),
        kpi('Fdyn', 'Dynamic slug force', Fdyn / 1e3, 'kN', 'Force range times the dynamic load factor of the support.', inp.Fallow, 'max'),
        kpi('dlf', 'Dynamic load factor', dlf, '-', 'Amplification of the slug force pulse by the support dynamics; 2 is the undamped upper bound for a sudden load.'),
      ],
      checks,
      plots: [
        plot('line', 'Liquid arrival at the outlet (slug arrival profile)', 'Time', 's', 'Liquid rate', 'm³/h', [
          { name: 'Outlet liquid rate', x: tt.map((v) => round(v)), y: ql.map((v) => round(v * 3600)) }, { name: 'Mean rate', x: [0, round(inp.duration)], y: [round(train.meanQL * 3600), round(train.meanQL * 3600)], dash: true },
        ], 'Synthetic slug train at the outlet: plateaus are slug bodies, troughs the film zones between them. Peak height and duration size the separator inlet and level control.'),
        { type: 'hist', title: 'Slug lengths in the train', values: lens.map((v) => round(v)), bins: 15, x: { label: 'Slug length', unit: 'm' }, y: { label: 'Count', unit: '-' },
          explain: 'Distribution of the slug lengths sampled for the arrival history. The long right-hand tail is why the design slug is several times the mean.' },
        plot('line', 'Slug frequency vs production rate', 'Mixture velocity', 'm/s', 'Slug frequency', '1/s', [
          { name: 'Zabaras correlation', x: fSweep.map((q) => round(q.vm)), y: fSweep.map((q) => round(q.f)) }, { name: 'Operating point', x: [round(vm)], y: [round(fZ)] },
        ], 'How slug frequency moves with rate at the same pressure. The minimum of the curve marks the transition from few long slugs at low rate to many short ones at high rate.'),
        ...(trk ? [plot('line', 'Slugs in the flowline (tracking)', 'Time', 's', 'Number of slugs', '-', [(() => { const k = Math.max(1, Math.floor(trk.t.length / 300)); const ix = Array.from({ length: Math.floor(trk.t.length / k) }, (_, i) => i * k); return { name: 'Slugs present', x: ix.map((i) => round(trk.t[i])), y: ix.map((i) => trk.count[i]) }; })()],
          'Number of slugs travelling in the flowline as tracked from the inlet. A steady plateau means slugs survive to the outlet; a count that stays near zero means they decay in the stratified film.')] : []),
        ...(hasRiser ? [plot('line', 'Severe-slugging boundary (Bøe)', 'Superficial gas velocity at the riser base', 'm/s', 'Superficial liquid velocity', 'm/s', [
          { name: 'Bøe boundary (severe slugging above)', x: vg.map((v) => round(v)), y: vg.map((v) => round(boe.vslCritical / Math.max(sb.vsg, 1e-6) * v)) }, { name: 'Operating point', x: [round(sb.vsg)], y: [round(sb.vsl)] },
        ], 'Operating point against the Bøe line. Points above the line can severe-slug: the riser head builds faster than the flowline gas pressure.')] : []),
      ],
      tables: [tbl('Slug unit cell and stability', [col('q', 'Quantity'), col('v', 'Value'), col('u', 'Unit')], [
        { q: 'Predicted flow regime in the flowline', v: reg.regime, u: '' }, { q: 'Unit-cell solution at the flowline angle (slug flow)', v: cellOk ? 'yes' : 'no — horizontal cell or stratified film used', u: '' },
        { q: 'Superficial liquid / gas velocity', v: `${fmt(s.vsl)} / ${fmt(s.vsg)}`, u: 'm/s' }, { q: 'Film holdup / film velocity', v: `${fmt(uc.filmHoldup)} / ${fmt(uc.vLF)}`, u: '- / m/s' },
        { q: 'Slug fraction of the unit', v: round(uc.slugFraction), u: '-' }, { q: 'Unit-cell pressure gradient', v: round(uc.gradTotal), u: 'Pa/m' },
        { q: 'Stratified flow: viscous Kelvin–Helmholtz stable', v: kh ? (kh.vkhStable ? 'yes' : 'no — waves grow into slugs') : 'n/a', u: '' }, { q: 'Stratified flow: inviscid Kelvin–Helmholtz stable', v: kh ? (kh.ikhStable ? 'yes' : 'no') : 'n/a', u: '' },
        { q: 'Riser-base pressure', v: round(base.P / BAR), u: 'bar' }, { q: 'Taitel stability criterion', v: hasRiser ? (tai.stable ? 'stable' : 'unstable') : 'n/a', u: '' },
        { q: 'Slug force cycles per year at the bend', v: round(cycles.cyclesPerYear, 3), u: '1/yr' }, { q: 'Allowable slug velocity for the bend', v: round(vAllow), u: 'm/s' },
      ], 'Intermediate results: the flow regime and unit-cell state, the stratified-flow stability that governs slug initiation, the riser stability criteria and the load cycle count for the fatigue module.')],
      notes: [slugFlow ? `Hydrodynamic slug flow is predicted in the flowline (${reg.regime}).` : `The regime map gives ${reg.regime}: the slug characteristics below are what the line would show if it slugs (for example at turndown).`,
        hasRiser ? (severe ? `Severe slugging is predicted: Bøe number ${fmt(boe.piSS)}, Taitel unstable; cycle about ${fmt(cyc?.period ?? 0)} s with ${fmt(cyc?.slugVolume ?? 0)} m³ per blow-out.` : `No severe slugging at this rate (Bøe number ${fmt(boe.piSS)}, Taitel ${tai.stable ? 'stable' : 'unstable'}${downhill ? '' : ', flowline not downhill'}).`) : 'No riser: severe slugging is not assessed.',
        'Frequency and length correlations scatter widely; use the spread between the three frequency values as the uncertainty.'],
      advice: [
        adv('engineering', sev(gov.uc), `Slugs of ${fmt(Lmean)} m mean and ${fmt(dist.maxLength)} m design length arrive every ${fmt(1 / Math.max(uc.frequency, 1e-9))} s at ${fmt(vt)} m/s; required surge volume ${fmt(surgeGov)} m³, dynamic bend force ${fmt(Fdyn / 1e3)} kN.`,
          gov.uc > 1 ? (gov.name.startsWith('Slug surge') ? `Provide at least ${fmt(surgeGov)} m³ of surge volume, raise the drain capacity, or control the slugs with topside choking.` : gov.name.startsWith('Dynamic') ? `Keep the slug velocity below ${fmt(vAllow)} m/s or stiffen and strengthen the support for ${fmt(Fdyn / 1e3)} kN.` : mit.chokeDpForStability / BAR < 2 * inp.pSep ? `Choke the riser top by about ${fmt(mit.chokeDpForStability / BAR)} bar, add riser-base gas lift, or keep the rate above the Bøe boundary.`
              : `The quasi-static Taitel criterion asks for ${fmt(mit.chokeDpForStability / BAR)} bar of choking, which is impractical for this flowline length: confirm the behaviour with the transient simulation and plan riser-base gas lift or active choke control.`)
            : 'The facility absorbs the predicted slugs. Repeat at turndown and start-up rates, where slugs are longer, and confirm severe slugging with the transient module.',
          'Taitel & Barnea (1990); Bøe (1981) / Taitel (1986); DNV-RP-D101 (2017) §3.12'),
        adv('economic', 'info', `The design slug needs ${fmt(surgeGov)} m³ of surge capacity against ${fmt(inp.vSurge)} m³ available.`,
          surgeGov > inp.vSurge ? `A larger slug catcher is roughly ${fmt(steelSC)} t of vessel steel and topsides weight; slug control by choking costs a few bar of back-pressure (deferred production) instead — compare both over the field life.` : 'No extra surge volume is needed; a smaller separator could be considered if late-life rates do not lengthen the slugs.',
          ECON_BASIS, `Slug-catcher steel of about ${fmt(steelSC)} t ≈ ${fmt(steelSC * CO2_STEEL)} t CO₂e (screening factor ${CO2_STEEL} t/t). Stable flow also avoids flaring from separator trips, usually the larger emission.`),
      ],
      publish: { 'pipeline.slugFrequency': uc.frequency, 'pipeline.slugVolume': governingVol, 'pipeline.slugVelocity': vt },
    };
  },
  verify: [{
    // v_sl = 1 m/s and v_sg = 3 m/s in a horizontal 0.3 m bore (constant properties: 60.08 kg/s liquid at 850 kg/m³, 12.72 kg/s gas at 60 kg/m³).
    // Gregory et al.: H_LS = 1/(1 + (v_m/8.66)^1.39) = 1/(1 + (4/8.66)^1.39) = 0.7453.
    // Bendiksen (Fr = 4/√(9.81·0.3) = 2.33 < 3.5): v_t = 1.05 v_m + 0.54 √(gD) = 4.2 + 0.926 = 5.126 m/s.
    name: 'Slug-body holdup (Gregory) and translational velocity (Bendiksen)',
    source: 'Gregory, Nicholson & Aziz (1978); Bendiksen (1984); hand calculation',
    inputs: { D: 320, t: 10, angle: 0, hRiser: 0, fluidModel: 'constant', mRate: 72.806, xGas: 0.17476, rhoLc: 850, rhoGc: 60, muLc: 2, Lflow: 2000 },
    expect: [{ kpi: 'slugHoldup', value: 0.7453, tol: 0.003 }, { kpi: 'slugVel', value: 5.126, tol: 0.005 }, { kpi: 'vm', value: 4, tol: 0.003 }],
  }],
};
