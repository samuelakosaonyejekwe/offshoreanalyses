// Pipe–soil interaction: embedment and contact pressure, axial and lateral breakout and residual resistance with
// low / best / high estimates, force–displacement mobilisation curves and the soil stiffnesses used by span analysis.
import { contactGeometry, verticalBearingUndrained, verticalBearingDrained, embedment, axialResistance, lateralResistance, mobilisationCurve,
  frictionRange, soilStiffness } from '../../engine/pipeline/pipesoil.js';
import { pipeInputs, pipeFrom, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, linspace, sev,
  ECON_BASIS, CONFIRM, LIMIT_NO_FE, MGO_T_PER_DAY, CO2_MGO } from './_common.js';

const G_S = 'Seabed soil';
const G_E = 'Embedment';
const G_A = 'Axial response';
const G_T = 'Lateral response';
const CLASSES = ['sand-loose', 'sand-medium', 'sand-dense', 'clay-very-soft', 'clay-soft', 'clay-firm', 'clay-stiff', 'clay-very-stiff', 'clay-hard'];

export default {
  id: 's08-pipe-soil',
  suite: 8,
  title: 'Pipe–soil interaction',
  summary: 'Derives embedment, axial and lateral soil resistance with their uncertainty range, and the mobilisation curves other analyses need.',
  covers: [6],
  fidelity: 'analytical',
  standards: ['DNV-RP-F114 (2017) §4.2–4.4 (pipe–soil interaction for submarine pipelines)', 'DNV-RP-F109 (2021) §3.4.6 (passive resistance)', 'DNV-RP-F105 (2006) §7.4 (soil stiffness)'],
  models: ['Vertical bearing resistance of a partially embedded pipe (undrained and drained)', 'Coulomb friction with wedging factor', 'Interface adhesion (total stress) axial model', 'Rankine passive resistance / Verley lateral models', 'Tri-linear force–displacement mobilisation', 'Lognormal low / best / high estimates'],
  limits: `Monotonic breakout from the as-laid embedment with tri-linear mobilisation curves. Cyclic berm growth, consolidation hardening of axial friction and rate effects enter only through the factors the user types. Adhesion, interface friction, lay factors and the F105 stiffness coefficients are recalled or site-specific values to be confirmed against pipe–soil tests. ${LIMIT_NO_FE}`,
  publishes: ['pipeline.embedment', 'pipeline.axialFriction', 'pipeline.lateralFriction', 'pipeline.soilKV', 'pipeline.soilKL'],
  searchTerms: 'pipe-soil interaction embedment axial friction lateral breakout DNV-RP-F114',
  inputs: [
    ...pipeInputs({ omit: ['steel'] }),
    sel('soil', 'Soil type', 'clay', [['clay', 'Clay (undrained)'], ['sand', 'Sand (drained)']], G_S, 'soil', 'Selects the bearing, axial and lateral models.'),
    num('su0', 'Undrained strength at mudline', 'kPa', 3, 0.2, 500, G_S, 'soil', 'Intact strength at the seabed surface.', { bind: 'site.su0', showIf: { key: 'soil', in: ['clay'] }, uncertain: { dist: 'lognormal', cov: 0.25 } }),
    num('suGrad', 'Strength gradient', 'kPa/m', 1.5, 0, 50, G_S, 'soil', 'Increase of strength with depth.', { bind: 'site.suGrad', showIf: { key: 'soil', in: ['clay'] } }),
    num('phi', 'Friction angle', 'deg', 32, 25, 45, G_S, 'soil', 'Drained friction angle of the sand.', { bind: 'site.phiDeg', showIf: { key: 'soil', in: ['sand'] } }),
    num('gammaSub', 'Submerged unit weight', 'N/m³', 6000, 2000, 13000, G_S, 'soil', 'Soil unit weight in water.', { bind: 'site.gammaSub' }),
    sel('soilClass', 'Soil class for dynamic stiffness', 'clay-soft', CLASSES, G_S, 'soil', `Selects the DNV-RP-F105 stiffness coefficients C_V, C_L. ${CONFIRM}`),
    sel('zMode', 'Embedment', 'computed', [['computed', 'From bearing resistance with lay effects'], ['given', 'Specified (survey)']], G_E, 'ic', 'Measured embedment from the as-laid survey overrides the prediction.'),
    num('zGiven', 'Specified embedment', 'mm', 100, 0, 1500, G_E, 'ic', 'Invert depth below the seabed.', { showIf: { key: 'zMode', in: ['given'] } }),
    num('kLay', 'Lay factor on contact force', '-', 2, 1, 5, G_E, 'ic', 'Touchdown contact force over submerged weight (1–3, higher for stiff pipe in shallow water).', { showIf: { key: 'zMode', in: ['computed'] } }),
    num('kDyn', 'Dynamic embedment factor', '-', 1.5, 1, 4, G_E, 'ic', 'Multiplier on depth for lay-vessel motion (1–3 on clay).', { showIf: { key: 'zMode', in: ['computed'] } }),
    num('alphaAd', 'Interface adhesion factor α', '-', 0.5, 0.05, 1, G_A, 'soil', 'Clay: axial shear over undrained strength (0.2–1.0 depending on sensitivity and coating roughness).', { showIf: { key: 'soil', in: ['clay'] }, uncertain: { dist: 'lognormal', cov: 0.3 } }),
    num('muAx', 'Interface friction coefficient tan δ', '-', 0.5, 0.1, 1.2, G_A, 'soil', 'Sand: axial interface friction (0.2–0.7).', { showIf: { key: 'soil', in: ['sand'] }, uncertain: { dist: 'lognormal', cov: 0.2 } }),
    num('resA', 'Axial residual / breakout ratio', '-', 0.8, 0.2, 1, G_A, 'soil', 'Residual resistance after large displacement over the peak.'),
    num('xA', 'Axial mobilisation displacement', 'mm', 3, 0.5, 100, G_A, 'soil', 'Displacement at which axial breakout is reached (a few millimetres).'),
    sel('latModel', 'Lateral model', 'coulomb-passive', [['coulomb-passive', 'Coulomb friction + Rankine passive resistance'], ['verley', 'Coulomb friction + Verley passive resistance (DNV-RP-F109)']], G_T, 'control', 'Both add a passive term from the embedment to Coulomb friction.'),
    num('muLat', 'Lateral Coulomb friction coefficient', '-', 0.2, 0.05, 1.5, G_T, 'soil', 'DNV-RP-F109: 0.2 on clay, 0.6 on sand.', { uncertain: { dist: 'lognormal', cov: 0.2 } }),
    num('muLatRes', 'Residual lateral friction coefficient', '-', 0.2, 0.05, 1.5, G_T, 'soil', 'Friction after the pipe has ridden out of its embedment.'),
    num('berm', 'Berm factor on passive resistance at breakout', '-', 1, 0.5, 4, G_T, 'soil', 'Above 1 for berms built by earlier cycles.'),
    num('xL', 'Lateral breakout displacement / D', '-', 0.2, 0.02, 1, G_T, 'soil', 'Typically 0.1–0.5 D on soft clay.'),
    num('xLres', 'Lateral residual displacement / D', '-', 2, 0.3, 6, G_T, 'soil', 'Travel after which the residual resistance applies (0.5–3 D).'),
    num('kRange', 'High / best estimate ratio', '-', 1.5, 1.05, 4, 'Uncertainty', 'criteria', 'Low estimate = best / ratio, high estimate = best × ratio, taken as the 5 % and 95 % fractiles of a lognormal distribution.'),
  ],
  compute(inp, ctx) {
    const p = pipeFrom({ ...inp, steel: 'api5l-x65' }, ctx);
    need(p.ws > 0 && p.wsEmpty > 0, 'The pipe must be heavier than water, empty and in operation, to rest on the seabed.');
    const D = p.Do, clay = inp.soil === 'clay';
    const soil = clay ? { model: 'undrained', su0: inp.su0 * 1e3, suGrad: inp.suGrad * 1e3, gammaSub: inp.gammaSub } : { model: 'drained', phiDeg: inp.phi, gammaSub: inp.gammaSub };
    const Qv = (z) => (clay ? verticalBearingUndrained({ D, z, su0: soil.su0, suGrad: soil.suGrad, gammaSub: soil.gammaSub }).Qv : verticalBearingDrained({ D, z, phiDeg: inp.phi, gammaSub: inp.gammaSub }).Qv);
    let z, zLay = 0, zStat = 0;
    if (inp.zMode === 'given') z = inp.zGiven / 1e3;
    else {
      zLay = run(() => embedment({ D, W: p.wsEmpty, kLay: inp.kLay, kDyn: inp.kDyn, ...soil }), 'As-laid embedment (the seabed cannot carry the pipe within two diameters)').z;
      zStat = run(() => embedment({ D, W: Math.max(p.ws, p.wsFlood), ...soil }), 'Static embedment').z;
      z = Math.max(zLay, zStat);
    }
    need(z < 2 * D, 'Embedment exceeds two diameters: the pipe is buried, not resting on the seabed.');
    const W = p.ws, geo = contactGeometry({ D, z: Math.max(z, 1e-6) });
    const suZ = clay ? soil.su0 + soil.suGrad * z : 0;
    const ax = run(() => axialResistance(clay ? { model: 'undrained', D, z: Math.max(z, 1e-6), W, alpha: inp.alphaAd, su: suZ, residualRatio: inp.resA } : { model: 'drained', D, z: Math.max(z, 1e-6), W, mu: inp.muAx, residualRatio: inp.resA }), 'Axial resistance');
    const latArgs = { model: inp.latModel, soil: inp.soil, D, W, mu: inp.muLat, muResidual: inp.muLatRes, su: suZ || undefined, gammaSub: inp.gammaSub, gammaDry: 18000, phiDeg: clay ? undefined : inp.phi, bermFactor: inp.berm };
    const lat = run(() => lateralResistance({ ...latArgs, z }), 'Lateral resistance');
    const muA = ax.breakout / W, muL = lat.breakout / W;
    const rng = (be) => { const r = frictionRange({ LE: be / inp.kRange, HE: be * inp.kRange }); return { LE: be / inp.kRange, BE: r.BE, HE: be * inp.kRange, cov: r.cov }; };
    const rA = rng(muA), rL = rng(muL), rAr = rng(ax.residual / W), rLr = rng(lat.residual / W);
    const mA = mobilisationCurve({ breakout: ax.breakout, residual: ax.residual, xBreakout: inp.xA / 1e3, xResidual: Math.max(10 * inp.xA / 1e3, inp.xA / 1e3 * 1.01) });
    const mL = mobilisationCurve({ breakout: lat.breakout, residual: lat.residual, xBreakout: inp.xL * D, xResidual: Math.max(inp.xLres, inp.xL * 1.01) * D });
    const massRatio = Math.min(Math.max(p.massOp / (p.rhoW * Math.PI * D * D / 4), 1.2), 2.0);
    const ks = run(() => soilStiffness({ D, massRatio, soilType: inp.soilClass }), 'Soil stiffness');
    const QvHalf = Qv(0.5 * D), Vmax = Math.max(p.ws, p.wsFlood, inp.zMode === 'given' ? 0 : inp.kLay * p.wsEmpty);
    const zs = linspace(0.005 * D, 0.6 * D, 25);
    const curveA = zs.map((zz) => { const g2 = clay ? soil.su0 + soil.suGrad * zz : 0; return axialResistance(clay ? { model: 'undrained', D, z: zz, W, alpha: inp.alphaAd, su: g2 } : { model: 'drained', D, z: zz, W, mu: inp.muAx }).breakout; });
    const curveL = zs.map((zz) => { try { return lateralResistance({ ...latArgs, z: zz, su: clay ? soil.su0 + soil.suGrad * zz : undefined }).breakout; } catch { return 0; } });
    const xa = [0, ...mA.x.slice(1), mA.x[mA.x.length - 1] * 2], xl = [0, ...mL.x.slice(1), mL.x[mL.x.length - 1] * 1.5];
    const checks = [
      chk('Embedment within model validity', z / D, 0.5, '-', 'DNV-RP-F114 (2017) §4.2',
        'Flags embedment beyond half a diameter, where the bearing, wedging and passive-resistance expressions are extrapolated and the pipe behaves as partly buried.'),
      chk('Vertical bearing at half-diameter embedment', Vmax, QvHalf, 'N/m', 'DNV-RP-F114 (2017) §4.2',
        'Guards against the seabed being too weak to carry the heaviest contact force before the pipe is half buried; if exceeded the pipe sinks until soil strength or buoyancy stops it.'),
    ];
    return {
      kpis: [
        kpi('zEmb', 'Embedment', z * 1e3, 'mm', 'Depth of the pipe invert below the seabed used for all resistances.'),
        kpi('zOverD', 'Embedment / diameter', z / D, '-', 'Relative embedment; the models are calibrated up to about 0.5.', 0.5, 'max'),
        kpi('B', 'Contact width', geo.B * 1e3, 'mm', 'Chord of the pipe in contact with the soil at this embedment.'),
        kpi('qContact', 'Mean contact pressure', W / geo.B / 1e3, 'kPa', 'Operating submerged weight over the contact width.'),
        kpi('zeta', 'Wedging factor', ax.zeta, '-', 'Ratio of the integrated normal contact force to the vertical load: 1.0 at zero embedment, 1.27 at half a diameter. It raises drained axial friction.'),
        kpi('Qv', 'Vertical bearing resistance at this embedment', Qv(Math.max(z, 1e-6)), 'N/m', 'Soil bearing resistance at the embedment used; equal to the governing vertical load when the embedment is computed.'),
        kpi('FA', 'Axial breakout resistance', ax.breakout, 'N/m', 'Peak axial resistance per metre of pipe.'),
        kpi('FAres', 'Axial residual resistance', ax.residual, 'N/m', 'Axial resistance after large displacement, used for expansion and walking.'),
        kpi('muA', 'Equivalent axial friction factor', muA, '-', 'Axial breakout resistance over submerged weight: the friction factor for expansion, walking and buckle feed-in.'),
        kpi('FL', 'Lateral breakout resistance', lat.breakout, 'N/m', 'Peak lateral resistance: Coulomb friction plus passive resistance of the soil ahead of the embedded pipe.'),
        kpi('FLres', 'Lateral residual resistance', lat.residual, 'N/m', 'Lateral resistance once the pipe has broken out and slides on the surface.'),
        kpi('muL', 'Equivalent lateral friction factor', muL, '-', 'Lateral breakout resistance over submerged weight: the factor for lateral buckling and stability.'),
        kpi('KV', 'Dynamic vertical soil stiffness', ks.KV / 1e6, 'MN/m²', 'Vertical stiffness per metre of pipe for free-span dynamics.'),
        kpi('KL', 'Dynamic lateral soil stiffness', ks.KL / 1e6, 'MN/m²', 'Lateral stiffness per metre of pipe for free-span dynamics.'),
      ],
      checks,
      plots: [
        plot('line', 'Mobilisation curves', 'Displacement', 'mm', 'Resistance', 'N/m', [
          { name: 'Axial', x: xa.map((x) => round(x * 1e3)), y: xa.map((x) => round(mA.force(x))) }, { name: 'Lateral', x: xl.map((x) => round(x * 1e3)), y: xl.map((x) => round(mL.force(x))) },
        ], 'Force–displacement curves for axial and lateral movement: rise to breakout, fall to residual. These are the springs to use in expansion, walking and buckling models.'),
        plot('line', 'Resistance vs embedment', 'Embedment / D', '-', 'Breakout resistance', 'N/m', [
          { name: 'Axial breakout', x: zs.map((v) => round(v / D)), y: curveA.map((v) => round(v)) }, { name: 'Lateral breakout', x: zs.map((v) => round(v / D)), y: curveL.map((v) => round(v)) },
          { name: 'Bearing resistance', x: zs.map((v) => round(v / D)), y: zs.map((v) => round(Qv(v))) }, { name: 'Embedment used', x: [round(z / D), round(z / D)], y: [0, round(Math.max(lat.breakout, ax.breakout))], dash: true },
        ], 'How resistances grow with embedment. The steepness of the lateral curve shows how sensitive lateral buckling and stability are to the embedment estimate.'),
      ],
      tables: [tbl('Friction factors: low, best and high estimates', [col('q', 'Quantity'), col('LE', 'Low estimate', '-'), col('BE', 'Best estimate', '-'), col('HE', 'High estimate', '-'), col('cov', 'Coefficient of variation', '-')], [
        { q: 'Axial breakout', LE: round(rA.LE, 3), BE: round(rA.BE, 3), HE: round(rA.HE, 3), cov: round(rA.cov, 3) }, { q: 'Axial residual', LE: round(rAr.LE, 3), BE: round(rAr.BE, 3), HE: round(rAr.HE, 3), cov: round(rAr.cov, 3) },
        { q: 'Lateral breakout', LE: round(rL.LE, 3), BE: round(rL.BE, 3), HE: round(rL.HE, 3), cov: round(rL.cov, 3) }, { q: 'Lateral residual', LE: round(rLr.LE, 3), BE: round(rLr.BE, 3), HE: round(rLr.HE, 3), cov: round(rLr.cov, 3) },
      ], 'Equivalent friction factors as low, best and high estimates. Use the high axial with the low lateral value for buckling, the low axial for end expansion, and the high lateral for buckle strain.')],
      notes: [inp.zMode === 'given' ? `Embedment specified at ${fmt(z * 1e3)} mm.` : `As-laid embedment ${fmt(zLay * 1e3)} mm (lay factor ${fmt(inp.kLay)}, dynamic factor ${fmt(inp.kDyn)}); static ${fmt(zStat * 1e3)} mm.`,
        `Passive resistance contributes ${fmt(lat.passive)} N/m of the lateral breakout; ${clay ? `strength at the invert ${fmt(suZ / 1e3)} kPa` : `K_p from φ = ${fmt(inp.phi)}°`}.`,
        `Static vertical stiffness ${ks.KVS ? fmt(ks.KVS / 1e3) + ' kN/m²' : 'not tabulated for this class'}; coefficients are the DNV-RP-F105 (2006) table values.`],
      advice: [
        adv('engineering', sev(Math.max(...checks.map((c) => c.uc))), `Embedment ${fmt(z * 1e3)} mm (${fmt(z / D)} D) gives equivalent friction factors of ${fmt(muA)} axial and ${fmt(muL)} lateral, with a range of ×/÷ ${fmt(inp.kRange)}.`,
          z / D > 0.5 ? 'The pipe is more than half buried: treat it as partly buried and derive resistances from a site-specific model or tests.' : 'Carry the low, best and high estimates into expansion, walking and buckling; calibrate adhesion and lay factors against the as-laid survey and pipe–soil tests.',
          'DNV-RP-F114 (2017) §3.3 (LE/BE/HE) and §4.2–4.4; DNV-RP-F110 (2018) §5'),
        adv('economic', 'info', `The friction range is ×/÷ ${fmt(inp.kRange)}; design quantities such as end expansion scale with 1/friction and buckle-trigger spacing with friction.`,
          `A site-specific pipe–soil test programme that narrows the range from ${fmt(inp.kRange)} to about 1.25 typically removes sleepers, anchors or spool length worth far more than the survey; value it before fixing the buckle-management design.`,
          'Value-of-information reasoning for site investigation; whole-life costing per ISO 15686-5',
          `A geotechnical vessel burns about ${MGO_T_PER_DAY} t of fuel per day (≈ ${fmt(MGO_T_PER_DAY * CO2_MGO)} t CO₂, screening); narrower friction bounds avoid over-built restraint, rock and steel for the life of the line.`),
      ],
      publish: { 'pipeline.embedment': z, 'pipeline.axialFriction': muA, 'pipeline.lateralFriction': muL, 'pipeline.soilKV': ks.KV, 'pipeline.soilKL': ks.KL },
    };
  },
  verify: [{
    // Drained axial resistance at half-diameter embedment: wedging factor ζ = 2 sin β/(β + sin β cos β) with β = π/2 → 4/π = 1.2732;
    // F_A = μ W ζ = 0.5·1000·1.2732 = 636.6 N/m.
    name: 'Wedging factor and drained axial resistance of a half-embedded pipe',
    source: 'White & Randolph (2007) wedging factor as used in DNV-RP-F114 §4.3; hand calculation',
    inputs: { D: 400, t: 20, tCoat: 0, tConc: 0, wsGiven: 1000, soil: 'sand', zMode: 'given', zGiven: 200, muAx: 0.5, soilClass: 'sand-medium', muLat: 0.6, muLatRes: 0.6 },
    expect: [{ kpi: 'zeta', value: 1.2732, tol: 0.001 }, { kpi: 'FA', value: 636.6, tol: 0.002 }],
  }, {
    // Undrained bearing at z/D = 0.1: N = min[6·0.1^0.25; 3.4·(10·0.1)^0.5] = min(3.374; 3.4) = 3.374 → 3.374·D·s_u = 3.374·0.4·5000 = 6748 N/m.
    // Soil buoyancy: cos β = 0.8, A_bm = (D²/4)(β − sin β cos β) = 0.04·(0.6435 − 0.48) = 0.00654 m²; 1.5·γ'·A_bm = 1.5·2000·0.00654 = 19.6 N/m.
    // Q_v = 6748 + 19.6 = 6767.6 N/m.
    name: 'Undrained vertical bearing resistance at z/D = 0.1',
    source: 'DNV-RP-F114 (2017) §4.2 closed form; hand calculation',
    inputs: { D: 400, t: 20, tCoat: 0, tConc: 0, wsGiven: 1000, soil: 'clay', su0: 5, suGrad: 0, gammaSub: 2000, zMode: 'given', zGiven: 40 },
    expect: [{ kpi: 'Qv', value: 6767.6, tol: 0.002 }],
  }],
};
