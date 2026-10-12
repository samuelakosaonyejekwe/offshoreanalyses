// Upheaval buckling of buried or trenched pipe and seabed embedment: required download over an imperfection
// (Palmer et al.), uplift resistance of the cover, cover height and rock quantities, and as-laid embedment.
import { restrainedAxialForce } from '../../engine/pipeline/expansion.js';
import { upheavalDownload, hobbsMinimum, palmerPhiW } from '../../engine/pipeline/buckling.js';
import { upliftResistanceSand, upliftResistanceClay, coverQuantities, embedment, contactGeometry } from '../../engine/pipeline/pipesoil.js';
import { pipeInputs, pipeFrom, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, linspace, sev,
  vesselCo2, CONFIRM, LIMIT_NO_FE, ECON_BASIS, BAR } from './_common.js';

const G_L = 'Operating state';
const G_I = 'Imperfection';
const G_C = 'Cover';
const G_S = 'Seabed and embedment';
const ROCK_T_PER_DAY = 8000; // fall-pipe vessel placement rate, t/day (screening)

export default {
  id: 's08-upheaval-embedment',
  suite: 8,
  title: 'Upheaval buckling and seabed embedment',
  summary: 'Sizes the cover that holds a buried hot pipe down over a seabed imperfection and estimates how deep the pipe embeds when laid.',
  covers: [17],
  fidelity: 'analytical',
  standards: ['DNV-RP-F110 (2018) §8 and App. B (upheaval buckling of buried pipelines)', 'Palmer et al. (1990), OTC 6335', 'Schaminée et al. (1990), OTC 6486', 'DNV-RP-F114 (2017) §4.2 (embedment)'],
  models: ['Upheaval buckling universal design curve (Palmer)', 'Prop imperfection beam model', 'Uplift resistance of cohesionless and cohesive cover', 'Hobbs vertical buckling mode', 'Vertical bearing capacity of a partially embedded pipe', 'Nonlinear seabed contact (bearing resistance vs depth)'],
  limits: `Palmer's design curve for a single symmetric imperfection in a small-slope elastic beam-column; the safety factor on download is project specific (DNV-RP-F110 uses partial factors on load and resistance, not built in). Uplift coefficients and bearing factors are recalled and site dependent. Ratcheting under cyclic loading, liquefaction of the backfill and finite-element upheaval analysis with the surveyed profile are outside the module. ${LIMIT_NO_FE}`,
  publishes: ['pipeline.embedment', 'quantities.rock_t'],
  searchTerms: 'upheaval buckling buried pipeline cover height uplift resistance embedment',
  inputs: [
    ...pipeInputs(),
    sel('forceMode', 'Effective compressive force', 'restrained', [['restrained', 'Fully restrained, from pressure and temperature'], ['given', 'Specified']], G_L, 'bc', 'A buried pipe is normally fully restrained.'),
    num('pi', 'Operating internal pressure', 'bar', 150, 0, 1500, G_L, 'bc', 'Internal pressure rise since trenching.', { showIf: { key: 'forceMode', in: ['restrained'] } }),
    num('Tinst', 'Temperature when buried', '°C', 5, -10, 40, G_L, 'ic', 'Pipe temperature when the cover was placed.', { bind: 'site.Tbed', showIf: { key: 'forceMode', in: ['restrained'] } }),
    num('Top', 'Operating temperature', '°C', 80, -20, 200, G_L, 'bc', 'Pipe wall temperature in the buried section.', { showIf: { key: 'forceMode', in: ['restrained'] }, uncertain: { dist: 'normal', cov: 0.05 } }),
    num('Hlay', 'Residual lay tension', 'kN', 0, 0, 5000, G_L, 'ic', 'Tension locked in at burial.', { showIf: { key: 'forceMode', in: ['restrained'] } }),
    num('Pgiven', 'Effective compressive force', 'kN', 2000, 0, 50000, G_L, 'bc', 'Compression positive.', { showIf: { key: 'forceMode', in: ['given'] } }),
    num('delta', 'Imperfection height', 'm', 0.3, 0.02, 2, G_I, 'ic', 'Height of the out-of-straightness (overbend) of the as-trenched profile; the result is very sensitive to it, so use the survey accuracy as a minimum.', { uncertain: { dist: 'lognormal', cov: 0.3 } }),
    sel('lengthMode', 'Imperfection length', 'prop', [['prop', 'Prop shape from the installation weight'], ['given', 'Specified']], G_I, 'ic', 'A pipe draped over a point high spot takes the prop length; a surveyed length can be typed instead.'),
    num('Limp', 'Specified imperfection length', 'm', 40, 2, 500, G_I, 'ic', 'Full length of the imperfection.', { showIf: { key: 'lengthMode', in: ['given'] } }),
    sel('cover', 'Cover material', 'rock', [['rock', 'Rock berm'], ['sand', 'Sand backfill'], ['clay', 'Clay backfill']], G_C, 'soil', 'Rock and sand use the Schaminée uplift model; clay uses the vertical-slip / flow-around model.'),
    num('Hcover', 'Cover above the pipe crown', 'm', 0.8, 0, 5, G_C, 'geometry', 'Height of cover provided over the top of the pipe.', { design: true }),
    num('gammaCover', 'Submerged unit weight of sand or clay cover', 'N/m³', 8500, 2000, 13000, G_C, 'soil', 'Rock uses the particle density and porosity below.', { showIf: { key: 'cover', in: ['sand', 'clay'] } }),
    num('fUp', 'Uplift coefficient f', '-', 0.5, 0.05, 1, G_C, 'soil', 'About 0.1 loose sand, 0.5 dense sand, 0.6–0.7 rock. Site specific.', { uncertain: { dist: 'lognormal', cov: 0.2 } }),
    num('suCover', 'Remoulded strength of clay backfill', 'kPa', 3, 0.2, 100, G_C, 'soil', 'Undrained strength of the backfill, usually well below the intact value.', { showIf: { key: 'cover', in: ['clay'] } }),
    { key: 'rock', label: 'Rock type', type: 'select', source: 'materials:rock', default: 'rock-granite', group: G_C, role: 'material', help: 'Particle density of the berm rock.', showIf: { key: 'cover', in: ['rock'] } },
    num('porosity', 'Berm porosity', '-', 0.4, 0.25, 0.5, G_C, 'soil', 'Void fraction of the placed rock.'),
    num('crest', 'Berm crest width', 'm', 1, 0.3, 10, G_C, 'geometry', 'Width of the flat top of the berm.'),
    num('slope', 'Berm side slope (horizontal per vertical)', '-', 2.5, 1, 6, G_C, 'geometry', 'Rock berms are typically 1:2.5 to 1:3.'),
    num('Lcover', 'Length to be covered', 'km', 2, 0.01, 500, G_C, 'geometry', 'Length of pipeline receiving this cover.'),
    num('sf', 'Safety factor on required download', '-', 1.3, 1, 3, 'Acceptance criteria', 'criteria', `Project specific. ${CONFIRM}`),
    sel('soilModel', 'Seabed response at laying', 'undrained', [['undrained', 'Clay (undrained)'], ['drained', 'Sand (drained)']], G_S, 'soil', 'Bearing model for as-laid embedment.'),
    num('su0', 'Seabed undrained strength at mudline', 'kPa', 3, 0.2, 500, G_S, 'soil', 'Clay strength at the seabed surface.', { bind: 'site.su0', showIf: { key: 'soilModel', in: ['undrained'] } }),
    num('suGrad', 'Strength gradient', 'kPa/m', 1.5, 0, 50, G_S, 'soil', 'Increase of undrained strength with depth.', { bind: 'site.suGrad', showIf: { key: 'soilModel', in: ['undrained'] } }),
    num('phi', 'Friction angle of the seabed sand', 'deg', 32, 25, 45, G_S, 'soil', 'Drained friction angle.', { bind: 'site.phiDeg', showIf: { key: 'soilModel', in: ['drained'] } }),
    num('gammaBed', 'Submerged unit weight of the seabed soil', 'N/m³', 6000, 2000, 13000, G_S, 'soil', 'Soil buoyancy term of the bearing resistance.', { bind: 'site.gammaSub' }),
    num('kLay', 'Lay factor on contact force', '-', 2, 1, 5, G_S, 'ic', 'Touchdown-zone contact force over submerged weight during laying (1–3).'),
    num('kDyn', 'Dynamic embedment factor', '-', 1.5, 1, 4, G_S, 'ic', 'Extra embedment from vessel-motion-induced pipe movement at touchdown (1–3 on clay).'),
  ],
  compute(inp, ctx) {
    const p = pipeFrom(inp, ctx);
    need(p.ws > 0 && p.wsEmpty > 0, 'The pipe must be negatively buoyant both empty and in operation for this assessment.');
    let P;
    if (inp.forceMode === 'given') P = inp.Pgiven * 1e3;
    else P = -restrainedAxialForce({ H: inp.Hlay * 1e3, deltaPi: inp.pi * BAR, Ai: p.Ai, As: p.As, deltaT: inp.Top - inp.Tinst, E: p.E, alpha: p.alpha, nu: p.nu }).S0;
    need(P > 0, 'The effective force is tensile: upheaval buckling cannot occur. Raise the temperature or pressure, or lower the lay tension.');
    const D = p.Do;
    const rockM = (ctx?.material && ctx.material(inp.rock)) || { rho: 2650 };
    const cq = (H, len) => coverQuantities({ D, H: Math.max(H, 1e-3), crestWidth: inp.crest, sideSlope: inp.slope, length: len, rhoSolid: inp.cover === 'rock' ? rockM.rho : 2650, porosity: inp.porosity, rhoWater: p.rhoW });
    const gCov = inp.cover === 'rock' ? cq(1, 1).gammaSub : inp.gammaCover;
    const model = inp.cover === 'clay' ? { soil: 'clay', D, gammaSub: gCov, su: inp.suCover * 1e3 } : { soil: 'sand', D, gammaSub: gCov, f: inp.fUp };
    const upArgs = { P, EI: p.EI, delta: inp.delta, wInstall: p.wsEmpty, wPipe: p.ws, safetyFactor: inp.sf, cover: model };
    if (inp.lengthMode === 'given') upArgs.L = inp.Limp;
    const up = run(() => upheavalDownload(upArgs), 'Upheaval download');
    const Rof = (H) => (H <= 0 ? 0 : inp.cover === 'clay' ? upliftResistanceClay({ D, H, su: inp.suCover * 1e3, gammaSub: gCov }).R : upliftResistanceSand({ D, H, gammaSub: gCov, f: inp.fUp }));
    const R = Rof(inp.Hcover);
    const hReq = up.coverHeight ?? (up.stableWithoutCover ? 0 : NaN);
    need(Number.isFinite(hReq), `Clay backfill of ${fmt(inp.suCover)} kPa cannot supply the required uplift resistance at any depth (flow-around limit): use rock or a stronger backfill.`);
    const uc = up.wRequired / (p.ws + R);
    const vert = run(() => hobbsMinimum({ mode: 'vertical', EI: p.EI, EA: p.EA, w: p.ws + R, muAxial: 0.5, D: p.D }), 'Hobbs vertical mode');
    // as-laid embedment: empty pipe with lay and dynamic factors, then static settlement under the operating weight
    const soil = inp.soilModel === 'undrained' ? { model: 'undrained', su0: inp.su0 * 1e3, suGrad: inp.suGrad * 1e3, gammaSub: inp.gammaBed } : { model: 'drained', phiDeg: inp.phi, gammaSub: inp.gammaBed };
    const eLay = run(() => embedment({ D, W: p.wsEmpty, kLay: inp.kLay, kDyn: inp.kDyn, ...soil }), 'As-laid embedment (the seabed cannot carry the pipe within two diameters — check the soil strength)');
    const eOp = run(() => embedment({ D, W: Math.max(p.ws, p.wsFlood), ...soil }), 'Embedment under the heaviest static weight');
    const z = Math.max(eLay.z, eOp.z);
    const geo = contactGeometry({ D, z });
    const qContact = Math.max(p.ws, p.wsFlood) / geo.B;
    const L = inp.Lcover * 1e3;
    const qNow = cq(inp.Hcover, L), qReq = cq(hReq, L);
    const isRock = inp.cover === 'rock';
    const tNow = inp.Hcover > 0 ? qNow.mass / 1e3 : 0, tReq = hReq > 0 ? qReq.mass / 1e3 : 0;
    const ds = linspace(0.05, Math.max(2 * inp.delta, 0.6), 16);
    const hOfDelta = ds.map((d) => { try { const u = upheavalDownload({ ...upArgs, delta: d }); return u.coverHeight ?? 0; } catch { return null; } });
    const hs = linspace(0, Math.max(2 * inp.Hcover, 2 * hReq, 1), 21);
    const days = (isRock ? Math.max(tNow, tReq) : 0) / ROCK_T_PER_DAY;
    return {
      kpis: [
        kpi('ucUpheaval', 'Upheaval utilisation', uc, '-', 'Required download (with the safety factor) over the pipe weight plus the uplift resistance of the cover. Above 1.0 the pipe can lift through the cover at the imperfection.', 1, 'max'),
        kpi('wReq', 'Required download', up.wRequired, 'N/m', 'Total downward force per metre needed to hold the pipe on the imperfection at the operating compressive force, including the safety factor.'),
        kpi('Rcover', 'Uplift resistance of the cover provided', R, 'N/m', 'Resistance of the cover to the pipe moving upwards, excluding the pipe weight.'),
        kpi('coverReq', 'Cover height required', hReq, 'm', 'Height of cover above the pipe crown that supplies the missing download. Zero means the pipe weight alone is enough.'),
        kpi('P', 'Effective compressive force', P / 1e3, 'kN', 'Compressive effective force in the buried section in operation.'),
        kpi('Limp', 'Imperfection length', up.L, 'm', 'Full length of the overbend; the prop length is where a pipe draped over the high spot touches down.'),
        kpi('phiL', 'Dimensionless imperfection length', up.phiL, '-', 'L·√(P/EI). Short imperfections (below 4.49) need the maximum dimensionless download of 0.0646.'),
        kpi('phiW', 'Dimensionless download', up.phiW, '-', 'w·EI/(δ·P²) from Palmer\'s universal curve.'),
        kpi('P0vert', 'Vertical buckling force with cover (straight pipe)', vert.P0min / 1e3, 'kN', 'Hobbs vertical-mode safe force with the cover resistance added to the weight; the imperfection lowers the real capacity, which is why the download criterion governs.'),
        kpi('zEmb', 'As-laid embedment', z * 1e3, 'mm', 'Depth of the pipe invert below the seabed after laying and settlement under the heaviest static weight.'),
        kpi('zOverD', 'Embedment / diameter', z / D, '-', 'Relative embedment. Above about 0.5 the pipe is half buried and lateral models lose validity.'),
        kpi('qContact', 'Mean contact pressure', qContact / 1e3, 'kPa', 'Heaviest static weight over the contact width at this embedment.'),
        kpi('rock', isRock ? 'Rock in the berm provided' : 'Cover volume provided', isRock ? tNow : (inp.Hcover > 0 ? qNow.volume : 0), isRock ? 't' : 'm³', 'Quantity of cover material for the specified height and length.'),
        kpi('rockReq', isRock ? 'Rock for the required cover' : 'Cover volume required', isRock ? tReq : (hReq > 0 ? qReq.volume : 0), isRock ? 't' : 'm³', 'Quantity of cover material if the cover is set to the required height.'),
      ],
      checks: [
        chk('Upheaval buckling over the imperfection', up.wRequired, p.ws + R, 'N/m', 'Palmer et al. (1990) OTC 6335; DNV-RP-F110 (2018) §8',
          'Guards against the hot pipe lifting out of its trench at an overbend. Once it starts the cover fails progressively and the pipe arches out of the seabed.'),
        chk('Embedment within model validity', z / D, 0.5, '-', 'DNV-RP-F114 (2017) §4.2; DNV-RP-F109 §3.4.6',
          'Flags embedment beyond half a diameter, where the bearing and lateral-resistance models are extrapolated and the pipe may sink further.'),
      ],
      plots: [
        plot('line', 'Required cover vs imperfection height', 'Imperfection height', 'm', 'Cover above crown', 'm', [
          { name: 'Required cover', x: ds.filter((_, i) => hOfDelta[i] !== null).map((d) => round(d)), y: hOfDelta.filter((h) => h !== null).map((h) => round(h)) },
          { name: 'Cover provided', x: [round(ds[0]), round(ds[ds.length - 1])], y: [inp.Hcover, inp.Hcover], dash: true }, { name: 'Design imperfection', x: [inp.delta], y: [round(hReq)] },
        ], 'Cover height needed as the as-trenched imperfection grows. The crossing with the provided cover is the largest imperfection that can be tolerated, which sets the survey accuracy required.'),
        plot('line', 'Download vs cover height', 'Cover above crown', 'm', 'Download', 'N/m', [
          { name: 'Pipe weight + uplift resistance', x: hs.map((h) => round(h)), y: hs.map((h) => round(p.ws + Rof(h))) },
          { name: 'Required download', x: [0, round(hs[hs.length - 1])], y: [round(up.wRequired), round(up.wRequired)], dash: true },
        ], 'Available hold-down force as the cover is raised against the required download. The crossing is the required cover height.'),
      ],
      tables: [tbl('Cover quantities', [col('case', 'Case'), col('H', 'Cover height', 'm'), col('area', 'Section area', 'm²'), col('vol', 'Volume', 'm³'), col('mass', 'Placed mass', 't'), col('base', 'Base width', 'm')], [
        { case: 'As specified', H: inp.Hcover, area: round(inp.Hcover > 0 ? qNow.area : 0), vol: round(inp.Hcover > 0 ? qNow.volume : 0), mass: round(tNow), base: round(qNow.baseWidth) },
        { case: 'Required', H: round(hReq), area: round(hReq > 0 ? qReq.area : 0), vol: round(hReq > 0 ? qReq.volume : 0), mass: round(tReq), base: round(qReq.baseWidth) },
      ], 'Berm or backfill cross-section, volume and placed mass over the covered length for the specified and for the required cover height.')],
      notes: [`Cover submerged unit weight ${fmt(gCov)} N/m³; uplift model: ${inp.cover === 'clay' ? 'vertical slip / flow-around (clay)' : 'Schaminée (cohesionless)'}.`,
        up.stableWithoutCover ? 'The pipe weight alone provides the required download.' : `The cover must supply ${fmt(up.coverResistance)} N/m of uplift resistance.`,
        `As-laid embedment ${fmt(eLay.z * 1e3)} mm (lay factor ${fmt(inp.kLay)}, dynamic factor ${fmt(inp.kDyn)}); static embedment under the heaviest weight ${fmt(eOp.z * 1e3)} mm.`],
      advice: [
        adv('engineering', sev(uc), `Upheaval utilisation ${fmt(uc)}: required download ${fmt(up.wRequired)} N/m against ${fmt(p.ws + R)} N/m available with ${fmt(inp.Hcover)} m of cover.`,
          uc > 1 ? `Raise the cover to at least ${fmt(hReq)} m over the imperfection, or reduce the imperfection by re-trenching; the requirement scales with the imperfection height.` : `${fmt(inp.Hcover - hReq)} m of cover is in reserve. Confirm the imperfection height against the out-of-straightness survey accuracy and repeat with the as-trenched profile.`,
          'Palmer et al. (1990) OTC 6335; DNV-RP-F110 (2018) §8 (partial-factor format to be applied by the project)'),
        adv('economic', 'info', isRock ? `The specified berm is ${fmt(tNow)} t of rock over ${fmt(inp.Lcover)} km; the required cover needs ${fmt(tReq)} t (${fmt(tReq - tNow)} t difference).` : `The specified backfill is ${fmt(inp.Hcover > 0 ? qNow.volume : 0)} m³ over ${fmt(inp.Lcover)} km.`,
          isRock ? `At about ${ROCK_T_PER_DAY} t per day for a fall-pipe vessel this is ${fmt(days)} vessel days. Spot-dumping only at surveyed overbends instead of continuous cover typically cuts the rock by more than half.` : 'Natural or mechanical backfill avoids quarried rock; verify its strength after consolidation before taking credit for it.',
          ECON_BASIS, `${vesselCo2(Math.max(days, 0.1))} Rock berms occupy a seabed strip ${fmt(qNow.baseWidth)} m wide; minimising cover by a better out-of-straightness survey reduces quarrying, shipping and footprint together.`),
      ],
      publish: { 'pipeline.embedment': z, 'quantities.rock_t': isRock ? tNow : 0 },
      quantities: { rock_t: round(isRock ? tNow : 0), vessel_days: { rockPlacement: round(days, 3) } },
    };
  },
  verify: [{
    // Schaminée uplift: R = γ' H D (1 + f H/D) = 9000·1·0.4·(1 + 0.5·1/0.4) = 8100 N/m (bare 400 mm pipe, sand cover).
    // Palmer download for a short imperfection (Φ_L = L√(P/EI) = 20·√(2e6/8.946e7) = 2.99 < 4.49): Φ_w = 0.0646,
    // w = Φ_w δ P²/EI = 0.0646·0.3·(2e6)²/8.946e7 = 866.6 N/m (EI = 207e9·π/64·(0.4⁴ − 0.36⁴) = 8.946e7 N·m²).
    name: 'Schaminée uplift resistance and Palmer required download for a short imperfection',
    source: 'Schaminée et al. (1990) OTC 6486; Palmer et al. (1990) OTC 6335; hand calculation',
    inputs: { D: 400, t: 20, tCoat: 0, tConc: 0, forceMode: 'given', Pgiven: 2000, delta: 0.3, lengthMode: 'given', Limp: 20, cover: 'sand', Hcover: 1, gammaCover: 9000, fUp: 0.5, sf: 1 },
    expect: [{ kpi: 'Rcover', value: 8100, tol: 0.002 }, { kpi: 'wReq', value: 866.6, tol: 0.003 }],
  }],
};
