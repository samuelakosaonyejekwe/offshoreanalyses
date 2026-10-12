// Rotor, blade and support-structure fatigue: lifetime Miner damage at the monopile mudline weld and the tower base
// from rainflow-counted coupled simulations over the operating wind-speed bins (joint wind–wave climate), blade-root
// flapwise fatigue from turbulence and edgewise fatigue from gravity on a Goodman diagram, with blade deflection,
// root strain and tip-to-tower clearance at rated thrust.
import { simulateCoupled, buildSupportStructure, rainflowCount } from '../../engine/aero/coupled.js';
import { windWaveBins } from '../../engine/aero/dlc.js';
import { bemSolve } from '../../engine/aero/bem.js';
import { ntmSigma, rotorEffectiveTurbulence } from '../../engine/aero/turbulence.js';
import { genericBladeDistributions, bladeStaticDeflection, tipClearance, bendingStrain, goodmanCycles, gravityEdgewiseCycles, genericLayup, laminateABD, firstPlyFailure } from '../../engine/aero/blade.js';
import { snCurve, cyclesToFailure, SN_CURVE_NAMES } from '../../engine/fatigue/sn.js';
import { fatigueLife } from '../../engine/fatigue/damage.js';
import { turbineInputs, getTurbine, steadyAt, supportInputs, supportSegments, sectionAt, foundationInputs, foundationOf, segmentsMass, site, adv, check, worst, sev, fmt, round, run, need, G, RPM,
  LIMIT_POLAR, LIMIT_COUPLED, LIMIT_CODE, CONFIRM, SUPPLIER, steelCo2Text } from './_common.js';

const YEAR_S = 365.25 * 86400;
const curves = ['B1', 'B2', 'C', 'C1', 'C2', 'D', 'E', 'F', 'F1', 'F3'].filter((c) => SN_CURVE_NAMES.includes(c));

export default {
  id: 's16-fatigue', suite: 16, title: 'Fatigue of support structure and blades',
  summary: 'Lifetime fatigue damage and life of the monopile mudline weld, the tower base and the blade root from the wind–wave climate, with damage-equivalent loads, blade deflection, strain and tip clearance.',
  covers: [6], fidelity: 'reduced-order',
  standards: ['DNV-RP-C203 (2021) Fatigue design of offshore steel structures (S–N curves, thickness effect, Miner sum)', 'DNV-ST-0126 (2021) §4.11 (design fatigue factors for support structures)', 'IEC 61400-1:2019 §7.6.3 (fatigue) and DLC 1.2', 'IEC 61400-3-1:2019 (joint wind–wave fatigue climate)', 'DNV-ST-0376 (2015) Rotor blades (composite fatigue, constant-life diagram)', 'ASTM E1049-85 (rainflow counting)'],
  models: ['Reduced-order coupled aero-hydro-servo-elastic simulation per wind-speed bin', 'Normal turbulence model with rotor averaging', 'Joint wind–wave bins (Weibull wind, conditional sea state)', 'Rainflow cycle counting', 'Palmgren–Miner linear damage with DNV-RP-C203 two-slope S–N curves and thickness correction', 'Damage-equivalent loads', 'Euler–Bernoulli blade beam with centrifugal stiffening', 'Blade element momentum loads on the blade', 'Composite laminate theory (ABD matrices) and Tsai–Wu first-ply failure', 'Composite fatigue on a linear Goodman constant-life diagram', 'Gravity edgewise reversal per revolution'],
  limits: `${LIMIT_COUPLED} Fatigue is extrapolated from one short simulation per wind bin (IEC requires six 10-minute seeds per bin); only the first fore-aft mode is retained; idling, start/stop transients, side-side response, wind–wave misalignment and wave loading outside the operating wind range are not counted, so the support-structure damage is a lower-bound screening value. Blade properties are generic distributions scaled to mass and frequency, the root section is an equivalent tube, and the flapwise load spectrum is a quasi-static narrow-band estimate. ${LIMIT_POLAR} Fracture-mechanics crack growth is outside the module. ${LIMIT_CODE}`,
  publishes: ['life.fatigueYears', 'life.designYears', 'owt.delMudlineLife', 'owt.fatigueDamage'],
  searchTerms: 'offshore wind monopile fatigue damage equivalent load S-N curve blade root fatigue Goodman',
  inputs: [
    ...turbineInputs(),
    ...supportInputs({ steel: false }),
    ...foundationInputs(),
    site('A', { role: 'bc' }), site('k', { role: 'bc' }),
    { key: 'Iref', label: 'Reference turbulence intensity', unit: '-', type: 'number', default: 0.12, min: 0.05, max: 0.2, group: 'Wind resource', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 }, help: 'IEC category C 0.12 (B 0.14, A 0.16); wake-added turbulence in an array raises the effective value (see the wake module).', bind: 'array.effectiveTI' },
    { key: 'seaTable', label: 'Sea state by wind speed (empty = generic fetch-limited wind sea)', type: 'table', group: 'Sea state', role: 'bc', columns: [{ key: 'U', label: 'Hub wind speed', unit: 'm/s', type: 'number' }, { key: 'Hs', label: 'Hs', unit: 'm', type: 'number' }, { key: 'Tp', label: 'Tp', unit: 's', type: 'number' }], default: [],
      help: 'Conditional normal sea state from the site scatter data. When empty, a generic fetch-limited wind-sea correlation (Shore Protection Manual, 100 km fetch, no swell) is used and flagged.' },
    { key: 'avail', label: 'Time in power production within the operating wind range', unit: '%', type: 'number', default: 95, min: 50, max: 100, group: 'Operation', role: 'ic', help: 'Operational availability state: share of the time between cut-in and cut-out that the turbine actually runs.' },
    { key: 'D0', label: 'Fatigue damage already accumulated (Miner sum, characteristic)', unit: '-', type: 'number', default: 0, min: 0, max: 1, group: 'Operation', role: 'ic', help: 'Initial condition for an existing asset (0 for a new structure).' },
    { key: 'designLife', label: 'Design life', unit: 'yr', type: 'number', default: 25, min: 5, max: 50, group: 'Acceptance', role: 'criteria', bind: 'life.designYears' },
    { key: 'snMono', label: 'S–N curve of the monopile girth weld (seawater with cathodic protection)', type: 'select', default: 'D', group: 'Fatigue resistance', role: 'material', options: curves.map((c) => ({ value: c, label: `DNV-RP-C203 curve ${c}` })), help: 'D for a girth weld made from both sides; C1 if ground flush. Free corrosion is not assumed.' },
    { key: 'snTower', label: 'S–N curve of the tower-base weld (in air)', type: 'select', default: 'D', group: 'Fatigue resistance', role: 'material', options: curves.map((c) => ({ value: c, label: `DNV-RP-C203 curve ${c}` })) },
    { key: 'scf', label: 'Stress concentration factor at the welds', unit: '-', type: 'number', default: 1.15, min: 1, max: 3, group: 'Fatigue resistance', role: 'material', uncertain: { dist: 'lognormal', cov: 0.1 }, help: 'Thickness transition and fabrication misalignment at the girth weld (DNV-RP-C203 §3.3.7); attachments need their own detail category.' },
    { key: 'dffMono', label: 'Design fatigue factor, monopile', unit: '-', type: 'number', default: 3, min: 1, max: 10, group: 'Acceptance', role: 'criteria', help: `DNV-ST-0126 Table 4-13: 3 where inspection and repair are not planned (submerged, below mudline), 2 in the splash zone with access, 1 in the atmospheric zone. ${CONFIRM}` },
    { key: 'dffTower', label: 'Design fatigue factor, tower', unit: '-', type: 'number', default: 2, min: 1, max: 10, group: 'Acceptance', role: 'criteria', help: CONFIRM },
    { key: 'flapHz', label: 'Blade first flapwise frequency at standstill', unit: 'Hz', type: 'number', default: 0.70, min: 0.2, max: 3, group: 'Blade', role: 'material' },
    { key: 'edgeHz', label: 'Blade first edgewise frequency at standstill', unit: 'Hz', type: 'number', default: 1.08, min: 0.3, max: 5, group: 'Blade', role: 'material' },
    { key: 'rootD', label: 'Blade-root diameter / blade length', unit: '-', type: 'number', default: 0.057, min: 0.03, max: 0.1, group: 'Blade', role: 'geometry', help: 'NREL 5-MW: 3.5 m root on a 61.5 m blade. The extreme fibre is half this diameter from the neutral axis.' },
    { key: 'cgFrac', label: 'Blade centre of mass / blade length', unit: '-', type: 'number', default: 0.333, min: 0.2, max: 0.5, group: 'Blade', role: 'geometry', help: 'NREL 5-MW blade: 20.5 m from the root on 61.5 m.' },
    { key: 'laminate', label: 'Root laminate', type: 'select', source: 'materials:polymer', default: 'gfrp-ud', group: 'Blade', role: 'material', help: 'Generic unidirectional laminate from the library (E, tensile and compressive strength) — replace with tested characteristic values.' },
    { key: 'mBlade', label: 'Wöhler exponent of the laminate', unit: '-', type: 'number', default: 10, min: 6, max: 16, group: 'Blade', role: 'material', help: 'About 10 for glass/epoxy, 14 for carbon/epoxy; confirm by coupon tests.' },
    { key: 'gammaBlade', label: 'Partial factor on blade fatigue stress', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Acceptance', role: 'criteria', help: `Combined material and consequence factor for composite fatigue (DNV-ST-0376 / IEC 61400-5); knock-downs for the multiaxial root laminate are on top. ${CONFIRM}` },
    { key: 'clearMin', label: 'Minimum residual tip-to-tower clearance / unloaded clearance', unit: '-', type: 'number', default: 0.3, min: 0, max: 0.8, group: 'Acceptance', role: 'criteria', help: 'IEC 61400-1 asks for no contact under factored loads; 30 % residual is a common design target.' },
    { key: 'binWidth', label: 'Wind-speed bin width', unit: 'm/s', type: 'number', default: 4, min: 1, max: 6, group: 'Solver', role: 'control', help: 'IEC uses 2 m/s bins; wider bins are faster and coarser.' },
    { key: 'simT', label: 'Simulated time per wind bin', unit: 's', type: 'number', default: 120, min: 60, max: 660, group: 'Solver', role: 'control', help: 'After a 30 s transient. 600 s with several seeds is the design standard; short runs under-sample the largest cycles.' },
    { key: 'dt', label: 'Time step', unit: 's', type: 'number', default: 0.2, min: 0.0125, max: 0.2, group: 'Solver', role: 'control', help: 'The damage-equivalent load changes by about 1 % between 0.2 s and 0.05 s for a first mode near 0.25 Hz; reduce the step for stiffer structures.' },
    { key: 'seed', label: 'Random seed', type: 'int', default: 11, min: 1, max: 99999, group: 'Solver', role: 'ic' },
  ],
  meshTime: { key: 'dt', label: 'Time step (s)', levels: [0.2, 0.1, 0.05], qoi: 'delMud', isSize: true },
  compute(inp, ctx) {
    const tb = getTurbine(inp), rotor = tb.rotor, R = rotor.tipRadius;
    need(inp.hubHeight > R + 5, 'The hub height leaves no clearance under the blade tip.');
    const segs = supportSegments(inp, { zHub: inp.hubHeight, hubOffset: 2.4 * tb.scale });
    const structure = run(() => buildSupportStructure({ segments: segs, rnaMass: rotor.rnaMass, zHub: inp.hubHeight, zTowerBase: inp.transZ, foundation: foundationOf(inp), nElem: 30, nModes: 1, zeta: 0.01 }), 'Support-structure model');   // first fore-aft mode: it carries the fatigue response
    const table = (inp.seaTable || []).filter((r) => r.U > 0 && r.Hs >= 0 && r.Tp > 0).sort((a, b) => a.U - b.U);
    const bins = run(() => windWaveBins({ weibullA: inp.weibA, weibullK: inp.weibK, vIn: rotor.cutIn, vOut: rotor.cutOut, binWidth: inp.binWidth, hubHeight: inp.hubHeight, ...(table.length >= 2 ? { seaState: table } : {}) }), 'Wind–wave bins');
    // section properties and S–N curves
    const secM = { D: inp.monoD, t: inp.monoT / 1e3 }, secT = sectionAt(segs, inp.transZ + 0.01);
    const W = (s) => Math.PI / 32 * (s.D ** 4 - (s.D - 2 * s.t) ** 4) / s.D;
    const cvM = snCurve(inp.snMono, 'seawater-cp'), cvT = snCurve(inp.snTower, 'air');
    const dmgOf = (rf, Wsec, curve, t) => { let d = 0; for (let i = 0; i < rf.ranges.length; i++) { const sr = rf.ranges[i] / Wsec; if (sr > 1e4) d += rf.counts[i] / cyclesToFailure({ curve, stressRange: sr, thickness: t, scf: inp.scf }); } return d; };
    const tTrans = 30, av = inp.avail / 100;
    let dMud = 0, dTow = 0, sumM = 0, sumT = 0, spec = null, worstBin = -1;
    const rows = bins.map((b, ib) => {
      const sigma1 = ntmSigma({ Vhub: b.U, Iref: inp.Iref });
      // turbulence on a fixed 0.1 s grid so that the realisation does not depend on the solver step
      const ws = rotorEffectiveTurbulence({ Vhub: b.U, sigma1, hubHeight: inp.hubHeight, R, duration: 0.1 * 2 ** Math.ceil(Math.log2((inp.simT + tTrans + 2) / 0.1)), dt: 0.1, seed: inp.seed + ib }).u;
      const wind = (t) => { const x = Math.max(0, t) / 0.1, i = Math.min(ws.length - 2, Math.floor(x)); return ws[i] + (ws[i + 1] - ws[i]) * (x - i); };
      const st0 = steadyAt(tb, b.U);
      const sim = run(() => simulateCoupled({ surface: tb.surface, turbine: { inertia: rotor.drivetrainInertia, hubHeight: inp.hubHeight }, controller: tb.controller, structure, wind,
        ...(b.Hs > 0.05 ? { waves: { Hs: Math.min(b.Hs, 0.5 * inp.depth), Tp: b.Tp, gamma: 3.3, seed: inp.seed + 100 + ib, depth: inp.depth, macCamyFuchs: true } } : {}),
        dt: inp.dt, duration: inp.simT + tTrans, tTransient: tTrans, shearExp: 0.14, towerWindDrag: true, initial: { omega: st0.omega, pitchDeg: st0.pitchDeg, mode: 'production' }, psdNfft: 256 }), `Simulation of the ${fmt(b.U)} m/s bin`);
      const i0 = Math.round(tTrans / inp.dt);
      const rfM = rainflowCount(sim.mudlineMoment.subarray(i0)), rfT = rainflowCount(sim.towerBaseMoment.subarray(i0));
      const perS = (v) => v / inp.simT;
      const dm = perS(dmgOf(rfM, W(secM), cvM, secM.t)) * YEAR_S * b.prob * av, dt = perS(dmgOf(rfT, W(secT), cvT, secT.t)) * YEAR_S * b.prob * av;
      const s4 = (rf) => { let s = 0; for (let i = 0; i < rf.ranges.length; i++) s += rf.counts[i] * rf.ranges[i] ** 4; return s; };
      sumM += perS(s4(rfM)) * YEAR_S * b.prob * av; sumT += perS(s4(rfT)) * YEAR_S * b.prob * av;
      if (worstBin < 0 || dm > spec.dm) { worstBin = ib; spec = { dm, rf: rfM, U: b.U }; }
      dMud += dm; dTow += dt;
      return { U: b.U, prob: b.prob, Hs: b.Hs, Tp: b.Tp, generic: b.generic, dm, dt, std: sim.stats.mudlineMoment.std, mean: sim.stats.mudlineMoment.mean, del: sim.del.mudlineMoment };
    });
    need(Number.isFinite(dMud) && Number.isFinite(dTow), 'The fatigue simulations diverged: reduce the time step.');
    // lifetime damage-equivalent moments at 1e7 cycles, Wöhler exponent 4
    const nEq = 1e7, delMud = (sumM * inp.designLife / nEq) ** 0.25, delTow = (sumT * inp.designLife / nEq) ** 0.25;
    const lifeM = fatigueLife({ damageRate: Math.max(dMud, 1e-12), designLife: inp.designLife, dff: inp.dffMono, accumulatedDamage: inp.D0 });
    const lifeT = fatigueLife({ damageRate: Math.max(dTow, 1e-12), designLife: inp.designLife, dff: inp.dffTower, accumulatedDamage: inp.D0 });
    // ---- blade ----
    const lam = ctx.material(inp.laminate);
    need(lam && lam.E > 0 && lam.Xt > 0 && lam.Xc > 0, `Laminate '${inp.laminate}' lacks modulus or strength data in the materials library: choose a laminate entry (glass or carbon UD).`);
    const gb = run(() => genericBladeDistributions({ length: rotor.bladeLength, mass: rotor.bladeMass, flapFrequency: inp.flapHz, edgeFrequency: inp.edgeHz }), 'Blade distributions');
    const c = inp.rootD * rotor.bladeLength / 2;
    const epsT = lam.Xt / lam.E, epsC = lam.Xc / lam.E;
    // rated operating point: distributed normal load → deflection, root moment and strain
    const Vr = tb.ratedWind, opR = steadyAt(tb, Vr);
    const bemR = run(() => bemSolve({ rotor, U: Vr, omega: opR.omega, pitchDeg: opR.pitchDeg }), 'BEM at rated wind');
    const xs = gb.x, rs = bemR.spanwise.map((q) => q.r - rotor.hubRadius), np = bemR.spanwise.map((q) => q.Np);
    const loadAt = (x) => { if (x <= rs[0]) return np[0] * x / Math.max(rs[0], 1e-6); if (x >= rs[rs.length - 1]) return np[np.length - 1] * Math.max(0, (rotor.bladeLength - x) / Math.max(rotor.bladeLength - rs[rs.length - 1], 1e-6)); let i = 1; while (rs[i] < x) i++; return np[i - 1] + (np[i] - np[i - 1]) * (x - rs[i - 1]) / (rs[i] - rs[i - 1]); };
    const defl = run(() => bladeStaticDeflection({ x: xs, EI: gb.EIflap, massPerLength: gb.massPerLength, load: xs.map(loadAt), omega: opR.omega, hubRadius: rotor.hubRadius }), 'Blade deflection');
    const tipZ = inp.hubHeight - R;
    const clr = tipClearance({ tipDeflection: defl.tipDeflection * 1.35, overhang: rotor.overhang, tiltDeg: rotor.tiltDeg, preconeDeg: rotor.preconeDeg, tipRadius: R, towerRadius: sectionAt(segs, Math.max(tipZ, inp.transZ)).D / 2 });
    const epsRated = bendingStrain(defl.rootMoment, gb.EIflap[0], c);
    // root laminate: triaxial lay-up (0/±45) of the selected ply by classical laminate theory; first-ply failure under the factored rated strain
    const ply = { E1: lam.E, E2: lam.E2 ?? lam.E / 4, G12: lam.G12 ?? lam.E / 10, nu12: lam.nu ?? 0.28, Xt: lam.Xt, Xc: lam.Xc, Yt: lam.Yt ?? lam.Xt / 25, Yc: lam.Yc ?? lam.Xc / 4, S: lam.S ?? lam.Xt / 15, density: lam.rho ?? 1900, thickness: 0.9e-3 };
    const layup = genericLayup('triax', ply, 8), abd = run(() => laminateABD(layup), 'Laminate stiffness');
    const fpf = run(() => firstPlyFailure(layup, { N: [-1.35 * inp.gammaBlade * epsRated * abd.Ex * abd.thickness, 0, 0] }, { criterion: 'tsai-wu' }), 'First-ply failure');
    // flapwise fatigue: quasi-static narrow-band response to turbulence, cycles at the rotor frequency (rotational sampling)
    let dFlap = 0; const bladeRows = [];
    for (const b of bins) {
      const op = steadyAt(tb, b.U); if (!(op.omega > 0)) continue;
      const h = 0.5, m1 = bemSolve({ rotor, U: b.U + h, omega: op.omega, pitchDeg: op.pitchDeg }).flapMoment, m0 = bemSolve({ rotor, U: b.U - h, omega: op.omega, pitchDeg: op.pitchDeg }).flapMoment;
      const mMean = 0.5 * (m0 + m1), sigM = Math.abs(m1 - m0) / (2 * h) * ntmSigma({ Vhub: b.U, Iref: inp.Iref });
      const epsMean = bendingStrain(mMean, gb.EIflap[0], c) * inp.gammaBlade, sigEps = bendingStrain(sigM, gb.EIflap[0], c) * inp.gammaBlade;
      // Rayleigh-distributed amplitudes: Miner sum by quadrature over the amplitude density a/σ²·exp(−a²/2σ²)
      let dPerCycle = 0; const na = 40;
      for (let i = 0; i < na; i++) { const a = (i + 0.5) / na * 5 * sigEps, pdf = a / sigEps ** 2 * Math.exp(-a * a / (2 * sigEps ** 2)), N = goodmanCycles({ amplitude: Math.max(a, 1e-9), mean: epsMean, Rt: epsT, Rc: epsC, m: inp.mBlade }).N; dPerCycle += (N > 0 ? 1 / N : 1) * pdf * 5 * sigEps / na; }
      const cyc = op.omega / (2 * Math.PI) * YEAR_S * b.prob * av, d = dPerCycle * cyc;
      dFlap += d; bladeRows.push({ U: b.U, mMean, sigM, epsMean, sigEps, d });
    }
    // edgewise: gravity reversal every revolution about the mean torque moment
    const meanOmega = bins.reduce((a, b) => a + steadyAt(tb, b.U).omega * b.prob, 0) / Math.max(bins.reduce((a, b) => a + b.prob, 0), 1e-9);
    const grav = gravityEdgewiseCycles({ bladeMass: rotor.bladeMass, cgRadius: inp.cgFrac * rotor.bladeLength, meanRotorSpeed: meanOmega, years: 1, availability: av * bins.reduce((a, b) => a + b.prob, 0) });
    const epsEdgeAmp = bendingStrain(grav.momentRange / 2, gb.EIedge[0], c) * inp.gammaBlade, epsEdgeMean = bendingStrain(bemR.edgeMoment, gb.EIedge[0], c) * inp.gammaBlade;
    const nEdge = goodmanCycles({ amplitude: epsEdgeAmp, mean: epsEdgeMean, Rt: epsT, Rc: epsC, m: inp.mBlade }).N;
    const dEdge = nEdge > 0 ? grav.cycles / nEdge : 99;
    const lifeFlap = 1 / Math.max(dFlap, 1e-12), lifeEdge = 1 / Math.max(dEdge, 1e-12);
    const cap = (v) => Math.min(v, 1e5);
    const lifeMin = Math.min(lifeM.designLife, lifeT.designLife, lifeFlap, lifeEdge);
    const checks = [
      check('Fatigue damage at the monopile mudline weld over the design life', inp.dffMono * (inp.D0 + dMud * inp.designLife), 1, '-', `DNV-RP-C203 (2021) §2.2 (Miner sum) curve ${inp.snMono} in seawater with cathodic protection; DNV-ST-0126 DFF = ${inp.dffMono}`, 'Miner sum of the rainflow-counted stress ranges times the design fatigue factor. Above 1 a fatigue crack is expected to initiate at the girth weld within the design life; the weld below water cannot be inspected or repaired economically.'),
      check('Fatigue damage at the tower base over the design life', inp.dffTower * (inp.D0 + dTow * inp.designLife), 1, '-', `DNV-RP-C203 (2021) curve ${inp.snTower} in air; DFF = ${inp.dffTower}`, 'Miner sum at the tower-base weld times the design fatigue factor. The tower base carries nearly the full rotor-induced moment range at a thinner wall than the pile.'),
      check('Blade-root flapwise fatigue damage over the design life', dFlap * inp.designLife, 1, '-', 'DNV-ST-0376 (2015) §2 / IEC 61400-5 (linear Goodman constant-life diagram, Miner sum)', 'Damage of the root laminate from turbulence-driven flapwise strain cycles with the partial factor applied. Above 1 the laminate is predicted to crack (matrix cracking, delamination) before the end of life.'),
      check('Blade-root edgewise fatigue damage over the design life (gravity)', dEdge * inp.designLife, 1, '-', 'DNV-ST-0376 (2015); gravity reversal once per revolution', 'Damage from the weight of the blade reversing the edgewise root moment every revolution — about 10⁸ fully reversed cycles in a lifetime, which governs the trailing-edge and root design of large blades.'),
      check('Blade-root laminate, first-ply failure at factored rated load (Tsai–Wu)', 1 / fpf.strengthRatio, 1, '-', 'Classical laminate theory with the Tsai–Wu criterion (Jones 1999); DNV-ST-0376 (2015) §2 (fibre failure, inter-fibre failure)', 'Failure index of the most stressed ply of a triaxial root laminate on the compression side under the rated flapwise strain times the load and material factors. Above 1 a ply cracks (first-ply failure), the start of progressive laminate damage.'),
      check('Blade-tip deflection against the tower clearance (load factor 1.35)', 1 - clr.clearance / clr.undeflected, 1 - inp.clearMin, '-', 'IEC 61400-1:2019 §7.6.5 (critical deflection analysis)', 'Share of the unloaded tip-to-tower distance used by the factored flapwise deflection at rated thrust against the share allowed. A blade strike on the tower is a total loss of rotor and often of the tower.'),
    ];
    const gov = worst(checks);
    // exceedance spectrum of mudline stress ranges in the most damaging bin and the S–N curve
    const srs = Array.from(spec.rf.ranges, (r, i) => [r / W(secM) * inp.scf / 1e6, spec.rf.counts[i]]).filter((q) => q[0] > 0.5).sort((a, b) => b[0] - a[0]);
    let cum = 0; const perYear = YEAR_S / inp.simT * (bins[worstBin].prob * av);
    const ex = srs.map(([s, n]) => { cum += n * perYear * inp.designLife; return [s, cum]; });
    const stepE = Math.max(1, Math.floor(ex.length / 120)), exT = ex.filter((_, i) => i % stepE === 0);
    const snS = [300, 200, 150, 100, 70, 50, 35, 25, 15, 10], snN = snS.map((s) => cyclesToFailure({ curve: cvM, stressRange: s * 1e6, thickness: secM.t }));
    const steel = segmentsMass(segs) / 1e3;
    // thickness change that brings the mudline damage to 1/DFF: stress ∝ 1/t, damage ∝ stress^m (m ≈ 4 between the slopes)
    const tScale = (inp.dffMono * (inp.D0 + dMud * inp.designLife)) ** (1 / 4);
    const nRef = cyclesToFailure({ curve: cvT, stressRange: 100e6 });
    return {
      kpis: [
        { key: 'fatigueYears', label: 'Governing design fatigue life', value: round(cap(lifeMin), 4), unit: 'yr', limit: inp.designLife, sense: 'min', explain: 'Shortest of the lives of the monopile weld, the tower base (both divided by their design fatigue factor) and the blade root. It must exceed the design life.' },
        { key: 'lifeMud', label: 'Design fatigue life at the mudline weld', value: round(cap(lifeM.designLife), 4), unit: 'yr', explain: `Years to a Miner sum of 1/DFF at the monopile girth weld (characteristic life ${fmt(cap(lifeM.characteristicLife))} years without the factor).` },
        { key: 'lifeTower', label: 'Design fatigue life at the tower base', value: round(cap(lifeT.designLife), 4), unit: 'yr', explain: 'Years to a Miner sum of 1/DFF at the tower-base weld.' },
        { key: 'lifeFlap', label: 'Blade-root flapwise fatigue life', value: round(cap(lifeFlap), 4), unit: 'yr', explain: 'Life of the root laminate under turbulence-driven flapwise cycles with the partial factor on strain.' },
        { key: 'lifeEdge', label: 'Blade-root edgewise fatigue life (gravity)', value: round(cap(lifeEdge), 4), unit: 'yr', explain: 'Life under the gravity reversal of the edgewise moment each revolution. It shortens rapidly as blades grow, since the moment scales with length to the fourth power.' },
        { key: 'delMud', label: 'Lifetime damage-equivalent mudline moment (m = 4, 10⁷ cycles)', value: round(delMud / 1e6, 5), unit: 'MN·m', explain: 'Constant moment range that would cause the lifetime damage in ten million cycles: the single number used to compare designs and sites.' },
        { key: 'delTower', label: 'Lifetime damage-equivalent tower-base moment', value: round(delTow / 1e6, 5), unit: 'MN·m', explain: 'The same measure at the tower base.' },
        { key: 'dmgMud', label: 'Annual damage at the mudline weld (characteristic)', value: round(dMud, 5), unit: '1/yr', explain: 'Miner damage accumulated per year of operation before the design fatigue factor.' },
        { key: 'tipDefl', label: 'Blade-tip deflection at rated thrust', value: round(defl.tipDeflection, 4), unit: 'm', explain: 'Out-of-plane deflection of the blade tip under the rated aerodynamic load with centrifugal stiffening (characteristic).' },
        { key: 'clearance', label: 'Residual tip-to-tower clearance (factored deflection)', value: round(clr.clearance, 4), unit: 'm', explain: `Distance left between blade tip and tower wall when the deflection is multiplied by 1.35 (unloaded clearance ${fmt(clr.undeflected)} m).` },
        { key: 'epsRoot', label: 'Blade-root flapwise strain at rated', value: round(epsRated * 1e6, 5), unit: 'µε', explain: `Extreme-fibre strain at the blade root under the rated load; the generic laminate fails at about ${fmt(epsT * 1e6)} µε in tension.` },
        { key: 'sigRoot', label: 'Blade-root flapwise stress at rated', value: round(epsRated * lam.E / 1e6, 5), unit: 'MPa', explain: 'Root strain times the laminate modulus: the stress in the load-carrying plies.' },
        { key: 'exLam', label: 'Axial modulus of the triaxial root laminate', value: round(abd.Ex / 1e9, 4), unit: 'GPa', explain: 'Effective in-plane modulus of a 0/±45 lay-up of the selected ply from classical laminate theory: the off-axis plies lower it against the unidirectional value.' },
        { key: 'fpf', label: 'First-ply failure strength ratio at factored rated load', value: round(Math.min(fpf.strengthRatio, 999), 4), unit: '-', limit: 1, sense: 'min', explain: `Load multiplier to first-ply failure (Tsai–Wu) of the root laminate; the critical ply is at ${fmt(fpf.criticalAngleDeg)}°.` },
        { key: 'edgeRange', label: 'Gravity edgewise moment range at the blade root', value: round(grav.momentRange / 1e3, 5), unit: 'kN·m', explain: 'Twice blade weight times the distance of its centre of mass from the root: the moment swing every revolution.' },
        { key: 'nRef', label: 'Tower S–N curve: cycles to failure at a 100 MPa range', value: round(nRef, 5), unit: 'cycles', explain: 'Reference point of the selected S–N curve (no stress concentration, reference thickness) for checking against DNV-RP-C203 Table 2-1.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Annual fatigue damage by wind-speed bin', categories: rows.map((r) => `${fmt(r.U)} m/s`), series: [{ name: 'Mudline weld', y: rows.map((r) => round(r.dm, 6)) }, { name: 'Tower base', y: rows.map((r) => round(r.dt, 6)) }], y: { label: 'Damage per year', unit: '1/yr' },
          explain: 'Contribution of each wind-speed bin (with its sea state and probability) to the yearly Miner sum. Bins around rated wind usually dominate: high thrust variation and frequent occurrence.' },
        { type: 'line', title: `Stress-range spectrum at the mudline, ${fmt(spec.U)} m/s bin, against the S–N curve`, x: { label: 'log₁₀ cycles in the design life', unit: '-' }, y: { label: 'Stress range', unit: 'MPa' }, series: [{ name: 'Cumulative cycles exceeding the range', x: exT.map((q) => round(Math.log10(Math.max(q[1], 1)), 4)), y: exT.map((q) => round(q[0], 4)) }, { name: `S–N curve ${inp.snMono} (cycles to failure)`, x: snN.map((n) => round(Math.log10(n), 4)), y: snS, dash: true }],
          explain: 'Rainflow-counted hot-spot stress ranges of the most damaging bin as an exceedance curve, with the S–N curve. Damage accumulates where the spectrum comes closest to the S–N line — typically in the many small cycles of the lower slope.' },
        { type: 'line', title: 'Blade deflection and bending moment at rated thrust', x: { label: 'Span from the root', unit: 'm' }, y: { label: 'Deflection (m) · moment (MN·m)', unit: '-' }, series: [{ name: 'Flapwise deflection (m)', x: Array.from(defl.x, (v) => round(v, 3)), y: Array.from(defl.deflection, (v) => round(v, 4)) }, { name: 'Flapwise bending moment (MN·m)', x: Array.from(defl.x, (v) => round(v, 3)), y: Array.from(defl.moment, (v) => round(v / 1e6, 4)) }],
          explain: 'Out-of-plane deflection and bending moment along the blade under the rated aerodynamic load. The moment is largest at the root; deflection grows towards the flexible outboard part.' },
        { type: 'line', title: 'Goodman constant-life diagram of the root laminate', x: { label: 'Mean strain', unit: 'µε' }, y: { label: 'Strain amplitude', unit: 'µε' },
          series: [{ name: 'Static envelope', x: [-epsC * 1e6, (epsT - epsC) / 2 * 1e6, epsT * 1e6], y: [0, (epsT + epsC) / 2 * 1e6, 0] }, { name: 'Edgewise gravity cycle (factored)', x: [round(epsEdgeMean * 1e6, 3)], y: [round(epsEdgeAmp * 1e6, 3)] }, { name: 'Flapwise 2σ amplitude per bin (factored)', x: bladeRows.map((r) => round(r.epsMean * 1e6, 3)), y: bladeRows.map((r) => round(2 * r.sigEps * 1e6, 3)) }],
          explain: 'Mean strain against strain amplitude with the static-strength triangle. Constant-life lines are scaled copies of the triangle: points close to the envelope have short lives, points near the origin practically infinite life.' },
      ],
      tables: [
        { title: 'Wind–wave bins and support-structure damage', columns: [{ key: 'U', label: 'Wind speed', unit: 'm/s' }, { key: 'p', label: 'Probability', unit: '%' }, { key: 'Hs', label: 'Hs', unit: 'm' }, { key: 'Tp', label: 'Tp', unit: 's' }, { key: 'mean', label: 'Mean mudline moment', unit: 'MN·m' }, { key: 'std', label: 'Std. dev.', unit: 'MN·m' }, { key: 'del', label: '1-Hz DEL', unit: 'MN·m' }, { key: 'dm', label: 'Damage mudline', unit: '1/yr' }, { key: 'dt', label: 'Damage tower base', unit: '1/yr' }],
          rows: rows.map((r) => ({ U: round(r.U, 3), p: round(r.prob * 100, 3), Hs: round(r.Hs, 3), Tp: round(r.Tp, 3), mean: round(r.mean / 1e6, 4), std: round(r.std / 1e6, 4), del: round(r.del / 1e6, 4), dm: round(r.dm, 6), dt: round(r.dt, 6) })),
          explain: 'One row per wind-speed bin: probability, associated sea state, moment statistics of the simulation and the yearly damage it contributes at both welds.' },
        { title: 'Blade-root flapwise fatigue by wind bin', columns: [{ key: 'U', label: 'Wind speed', unit: 'm/s' }, { key: 'M', label: 'Mean flap moment', unit: 'MN·m' }, { key: 's', label: 'Std. dev. of moment', unit: 'MN·m' }, { key: 'em', label: 'Mean strain (factored)', unit: 'µε' }, { key: 'es', label: 'Std. dev. of strain (factored)', unit: 'µε' }, { key: 'd', label: 'Damage', unit: '1/yr' }],
          rows: bladeRows.map((r) => ({ U: round(r.U, 3), M: round(r.mMean / 1e6, 4), s: round(r.sigM / 1e6, 4), em: round(r.epsMean * 1e6, 4), es: round(r.sigEps * 1e6, 4), d: Number(r.d.toExponential(3)) })),
          explain: 'Mean and fluctuating flapwise root moment and strain per wind bin with the resulting yearly damage of the root laminate.' },
      ],
      notes: [
        `${rows.length} wind bins of ${inp.binWidth} m/s, ${inp.simT} s each (seed ${inp.seed}); sea states ${rows.some((r) => r.generic) ? 'from the generic fetch-limited wind-sea correlation — replace with the site table' : 'from the entered table'}; first natural frequency ${fmt(structure.modes[0].frequency)} Hz.`,
        `Support structure: annual damage ${fmt(dMud)} at the mudline (curve ${inp.snMono}, SCF ${inp.scf}, t = ${inp.monoT} mm) and ${fmt(dTow)} at the tower base; the ${fmt(spec.U)} m/s bin contributes most.`,
        `Blade: tip deflection ${fmt(defl.tipDeflection)} m at rated, root strain ${fmt(epsRated * 1e6)} µε; flapwise life ${fmt(cap(lifeFlap))} yr, edgewise (gravity) life ${fmt(cap(lifeEdge))} yr with γ = ${inp.gammaBlade}.`,
        `Governing check: ${gov.name} (utilisation ${fmt(gov.uc)}). Short single-seed runs: treat lives as indicative within a factor of about two.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `Governing design fatigue life is ${fmt(cap(lifeMin))} years against ${inp.designLife} required (mudline ${fmt(cap(lifeM.designLife))}, tower base ${fmt(cap(lifeT.designLife))}, blade flap ${fmt(cap(lifeFlap))}, blade edge ${fmt(cap(lifeEdge))} years).`,
          checks[0].uc > 1 ? `Increase the monopile wall at the mudline to about ${fmt(inp.monoT * tScale, 3)} mm, or grind the weld flush (curve C1), or reduce the load ranges with a tower damper; ${fmt((1 - 1 / tScale) * 100)} % lower stress ranges are needed.` : gov.uc > 1 ? `'${gov.name}' is exceeded (utilisation ${fmt(gov.uc)}): ${gov === checks[1] ? `increase the tower-base wall by about ${fmt((gov.uc ** 0.25 - 1) * 100)} %` : gov === checks[5] ? 'stiffen the blade (carbon spar caps), add prebend or increase overhang, tilt or precone' : 'thicken the root laminate or reduce the partial factor by testing'}.` : `Fatigue is satisfied with a reserve factor of ${fmt(lifeMin / inp.designLife)} on life. Run six 10-minute seeds per 2 m/s bin and add idling and transient cases before releasing steel; ${SUPPLIER}`,
          'DNV-RP-C203 (2021) and DNV-ST-0126 (2021) §4.11 for steel; DNV-ST-0376 (2015) / IEC 61400-5:2020 for blades; IEC 61400-1:2019 DLC 1.2 for the load basis.',
          `Support structure ${fmt(steel)} t. ${steelCo2Text(steel)} Fatigue life beyond the design life is the basis for lifetime extension — the lowest-carbon energy a wind farm can deliver.`),
        adv('economic', 'info', `The fatigue reserve is ${fmt(lifeMin / inp.designLife)} × the design life; the mudline weld uses ${fmt(inp.dffMono * dMud * inp.designLife * 100)} % of its allowable damage in ${inp.designLife} years.`,
          lifeMin > 1.3 * inp.designLife ? `The reserve supports a life-extension case of about ${fmt(Math.min(lifeMin - inp.designLife, 15))} years, or a wall reduction of roughly ${fmt((1 - (inp.designLife / Math.min(lifeM.designLife, 1e4)) ** 0.25) * 100)} % at the mudline (${fmt(steel * 0.4 * (1 - (inp.designLife / Math.min(lifeM.designLife, 1e4)) ** 0.25))} t) if no other limit state governs.` : `Fatigue governs the steel: each 10 % of extra wall thickness buys about ${fmt((1.1 ** 4 - 1) * 100)} % more life (damage ∝ stress⁴); a weld-improvement step (grinding) is usually cheaper than steel.`,
          'Whole-life costing per ISO 15686-5 and risk-based inspection per DNV-RP-G101 / DNV-RP-C210: fatigue reserve is valued through inspection intervals and life extension (DNV-ST-0262).'),
      ],
      publish: { 'life.fatigueYears': cap(lifeMin), 'life.designYears': inp.designLife, 'owt.delMudlineLife': delMud, 'owt.fatigueDamage': inp.D0 + dMud * inp.designLife },
      quantities: { steel_t: round(steel, 4) },
    };
  },
  verify: [
    { name: 'Gravity edgewise moment range at the blade root: 2·m·g·r_cg', source: 'Statics of the rotating blade; NREL 5-MW blade (17 740 kg, centre of mass 20.475 m from the root; Jonkman et al. 2009, Table 2-2)',
      // 2 × 17 740 kg × 9.80665 m/s² × (0.33293 × 61.5 m = 20.475 m) = 7 124.1 kN·m
      inputs: { cgFrac: 0.33293 }, expect: [{ kpi: 'edgeRange', value: 7124.1, tol: 0.001 }] },
    { name: 'DNV-RP-C203 curve D in air at Δσ = 100 MPa: N = 10^12.164 / 100³', source: 'DNV-RP-C203 (2021) Table 2-1: curve D, m₁ = 3, log ā₁ = 12.164',
      // 10^12.164 = 1.4588e12 → N = 1.4588e6 cycles
      inputs: { snTower: 'D' }, expect: [{ kpi: 'nRef', value: 1.4588e6, tol: 0.002 }] },
  ],
};
