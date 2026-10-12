// Global buckling: restrained effective force, Hobbs lateral and vertical thermal-buckling modes, infinite-mode
// screening, post-buckling amplitude and strain, buckle interaction and trigger spacing, and propagation buckling.
import { restrainedAxialForce } from '../../engine/pipeline/expansion.js';
import { hobbsLateral, hobbsVertical, hobbsMinimum, hobbsCurve, infiniteModeForce, buckleSpacing } from '../../engine/pipeline/buckling.js';
import { propagationPressure, localBucklingDisplacementControlled, dnvResistanceFactors } from '../../engine/pipeline/strength.js';
import { characteristicStrength, fabricationFactor } from '../../engine/pipeline/pipeprops.js';
import { pipeInputs, pipeFrom, depthInput, safetyClassInput, num, sel, kpi, chk, chkUc, adv, plot, bar, tbl, col, run, need, fmt, round, sev, arr,
  vesselCo2, CONFIRM, LIMIT_F101, LIMIT_NO_FE, STD_F101, ECON_BASIS, G, BAR } from './_common.js';

const G_L = 'Operating state';
const G_S = 'Pipe–soil friction';
const MODES = [1, 2, 3, 4, 'inf'];

export default {
  id: 's08-global-buckling',
  suite: 8,
  title: 'Global buckling — lateral, upheaval and propagation',
  summary: 'Tells whether the restrained hot pipeline can buckle sideways or upwards, how severe the buckle would be, and whether a collapse would propagate.',
  covers: [3],
  fidelity: 'analytical',
  standards: ['DNV-RP-F110 (2018) §5–7 (global buckling of exposed pipelines)', 'Hobbs (1984), J. Transp. Eng. ASCE 110(2)', `${STD_F101} §5.4.5 (propagation buckling), §5.4.6 (strain criterion)`],
  models: ['Hobbs lateral buckling modes 1–4 and infinite mode', 'Hobbs vertical (upheaval) mode on a rigid seabed', 'Fully restrained effective axial force', 'Euler beam-column buckling with Coulomb friction', 'Axial friction feed-in and virtual anchor spacing', 'Buckle propagation pressure'],
  limits: `Hobbs analysis of a perfectly straight elastic pipe on a rigid frictional seabed: it gives the lowest force at which a buckled shape can exist and the post-buckling equilibrium, not the snap load of a real imperfection. Imperfection-sensitive finite-element lateral buckling, trawl-triggered buckles, sleepers and buoyancy sections are outside this module. ${LIMIT_F101} ${LIMIT_NO_FE}`,
  publishes: ['pipeline.effectiveForce', 'pipeline.buckleSafeForce'],
  searchTerms: 'pipeline lateral buckling Hobbs upheaval thermal buckling DNV-RP-F110',
  inputs: [
    ...pipeInputs(),
    depthInput(120, G_L),
    num('pi', 'Operating internal pressure', 'bar', 150, 0, 1500, G_L, 'bc', 'Internal pressure in operation; the pressure at installation is taken as zero.', { uncertain: { dist: 'normal', cov: 0.03 } }),
    num('Tinst', 'Installation temperature', '°C', 5, -10, 40, G_L, 'ic', 'Pipe temperature when laid and restrained.', { bind: 'site.Tbed' }),
    num('Top', 'Operating temperature', '°C', 80, -20, 200, G_L, 'bc', 'Pipe wall temperature in the restrained section.', { design: true, uncertain: { dist: 'normal', cov: 0.05 } }),
    num('Hlay', 'Residual lay tension', 'kN', 100, 0, 5000, G_L, 'ic', 'Effective tension left in the pipe after laying; it delays buckling.'),
    num('muA', 'Axial friction coefficient', '-', 0.5, 0.05, 1.5, G_S, 'soil', 'Use the pipe–soil module for a site-specific low / best / high estimate.', { bind: 'pipeline.axialFriction', uncertain: { dist: 'lognormal', cov: 0.25 } }),
    num('muL', 'Lateral friction coefficient', '-', 0.6, 0.05, 2.5, G_S, 'soil', 'Breakout lateral resistance over submerged weight. A low value makes buckling easier but gentler.', { bind: 'pipeline.lateralFriction', uncertain: { dist: 'lognormal', cov: 0.25 } }),
    num('spacing', 'Spacing of planned buckle triggers (0 = none)', 'm', 2000, 0, 50000, 'Buckle management', 'geometry', 'Distance between sleepers, snake-lay curves or buoyancy sections; sets the feed-in to each buckle.'),
    sel('process', 'Pipe manufacturing process', 'seamless', ['seamless', 'UO', 'TRB', 'ERW', 'HFW', 'UOE'], 'Pipe', 'material', 'Fabrication factor for propagation and collapse.'),
    safetyClassInput('Acceptance criteria'),
    num('alphaGw', 'Girth-weld factor α_gw', '-', 1, 0.5, 1, 'Acceptance criteria', 'criteria', `For the strain criterion in the buckle crown. ${CONFIRM}`),
    num('scfStrain', 'Strain concentration factor', '-', 1.2, 1, 3, 'Acceptance criteria', 'criteria', 'On the buckle bending strain for field joints and weight-coating discontinuities.'),
  ],
  compute(inp, ctx) {
    const p = pipeFrom(inp, ctx);
    const pe = p.rhoW * G * inp.depth, pi = inp.pi * BAR, dT = inp.Top - inp.Tinst;
    need(p.ws > 0, 'The pipe is buoyant in operation: lateral and upheaval resistance cannot be evaluated.');
    const r = restrainedAxialForce({ H: inp.Hlay * 1e3, deltaPi: pi, Ai: p.Ai, As: p.As, deltaT: dT, E: p.E, alpha: p.alpha, nu: p.nu });
    const P0 = -r.S0; // compressive positive
    const com = { EI: p.EI, EA: p.EA, w: p.ws, muLateral: inp.muL, muAxial: inp.muA, D: p.D, alpha: p.alpha, deltaPi: pi, Ai: p.Ai, nu: p.nu, H: inp.Hlay * 1e3 };
    const mins = MODES.map((mode) => ({ mode, ...run(() => hobbsMinimum({ ...com, mode }), `Hobbs mode ${mode}`) }));
    const vert = run(() => hobbsMinimum({ mode: 'vertical', EI: p.EI, EA: p.EA, w: p.ws, muAxial: inp.muA, D: p.D, alpha: p.alpha, deltaPi: pi, Ai: p.Ai, nu: p.nu, H: inp.Hlay * 1e3 }), 'Hobbs vertical mode');
    const gov = mins.reduce((a, m) => (m.P0min < a.P0min ? m : a), mins[0]);
    const inf = infiniteModeForce({ EA: p.EA, EI: p.EI, fL: inp.muL * p.ws, S0: Math.max(P0, 1) });
    const ucLat = Math.max(P0, 0) / gov.P0min, ucVert = Math.max(P0, 0) / vert.P0min;
    // post-buckling equilibrium of the governing lateral mode at the operating force (stable branch, L > L at the minimum)
    let post = null;
    if (ucLat > 1) {
      const c = hobbsCurve({ ...com, mode: gov.mode, Lmin: gov.L, Lmax: 12 * gov.L, n: 400 });
      for (let i = 1; i < c.L.length; i++) {
        if (c.P0[i] >= P0) {
          const L = c.L[i - 1] + (c.L[i] - c.L[i - 1]) * (P0 - c.P0[i - 1]) / (c.P0[i] - c.P0[i - 1]);
          post = { L, ...hobbsLateral({ ...com, mode: gov.mode, L }) }; break;
        }
      }
    }
    const str = characteristicStrength({ smys: p.smys, smts: p.smts, tempC: Math.max(inp.Top, 20) });
    const sc = inp.safetyClass, tab = dnvResistanceFactors(sc);
    const epsSd = post ? post.strain * inp.scfStrain : 0;
    const dcc = post ? run(() => localBucklingDisplacementControlled({ D: p.D, t2: p.t, fy: str.fy, fu: str.fu, epsSd, pMin: pi, pe, alphaGw: inp.alphaGw, safetyClass: sc, E: p.E, nu: p.nu }), 'Strain criterion in the buckle') : null;
    const f = inp.muA * p.ws;
    const sp = post ? buckleSpacing({ P0, Ppost: post.P, f, EA: p.EA, spacing: inp.spacing > 0 ? inp.spacing : undefined, PcritUnplanned: gov.P0min }) : null;
    const prop = run(() => propagationPressure({ D: p.D, t2: p.t, fy: str.fy, alphaFab: fabricationFactor(inp.process), pe, pMin: 0, safetyClass: sc }), 'Propagation buckling');
    const Tcrit = inp.Tinst + gov.deltaT, TcritV = inp.Tinst + vert.deltaT;
    const checks = [
      chk('Propagating buckling', prop.demand, prop.capacity, 'Pa', `${STD_F101} Eq. 5.16`,
        'Guards against a local collapse (for example at a buckle crown) running along the line under external pressure; above 1.0 arrestors are needed.'),
    ];
    if (dcc) checks.push(chkUc('Local buckling in the buckle crown (strain criterion)', dcc.uc, `${STD_F101} Eq. 5.30`,
      'Guards against the wall wrinkling at the crown of the lateral buckle, where the bending strain is displacement controlled. Above 1.0 the uncontrolled buckle is not acceptable.'));
    const Ls = Array.from({ length: 60 }, (_, i) => gov.L * (0.45 + i * 0.06));
    const curves = [1, 2, 3, 4].map((mode) => { const c = hobbsCurve({ ...com, mode, Lmin: Ls[0], Lmax: Ls[Ls.length - 1], n: 60 }); return { name: `Mode ${mode}`, x: arr(c.L), y: Array.from(c.P0, (v) => (Number.isFinite(v) ? round(v / 1e3) : null)) }; })
      .map((s) => { const k = s.y.map((v, i) => (v === null ? -1 : i)).filter((i) => i >= 0); return { name: s.name, x: k.map((i) => s.x[i]), y: k.map((i) => s.y[i]) }; });
    const worstUc = Math.max(dcc ? dcc.uc : 0, prop.uc);
    const status = ucLat <= 1 ? 'ok' : dcc && dcc.uc > 1 ? 'fail' : 'warn';
    const nTrig = inp.spacing > 0 ? 1e4 / inp.spacing : 0;
    return {
      kpis: [
        kpi('S0', 'Fully restrained effective force', r.S0 / 1e3, 'kN', 'Effective axial force where the pipe cannot expand: lay tension minus the pressure and thermal terms. Negative is compression and drives buckling.'),
        kpi('P0crit', 'Lowest lateral buckling force', gov.P0min / 1e3, 'kN', `Minimum of the Hobbs curves (mode ${gov.mode}). Below this compressive force no lateral buckle can exist on an even seabed.`),
        kpi('ucLat', 'Lateral buckling susceptibility', ucLat, '-', 'Restrained compressive force over the lowest lateral buckling force. Above 1.0 the pipeline is susceptible to lateral buckling and the buckle must be shown acceptable (strain check) or be managed; it is a screening ratio, not a failure criterion.'),
        kpi('Tcrit', 'Operating temperature at which buckling becomes possible', Tcrit, '°C', 'Wall temperature at which the restrained force reaches the lowest lateral buckling force at the stated pressure.'),
        kpi('P0vert', 'Lowest vertical buckling force', vert.P0min / 1e3, 'kN', 'Minimum of the Hobbs vertical-mode curve for a pipe resting on the seabed without cover.'),
        kpi('ucVert', 'Upheaval susceptibility (exposed pipe)', ucVert, '-', 'Restrained force over the lowest vertical buckling force. An exposed pipe that exceeds both limits buckles laterally, the lower of the two.'),
        kpi('Sinf', 'Infinite-mode screening force', inf.Sinf / 1e3, 'kN', 'Closed-form minimum of the Hobbs infinite mode, the even-seabed screening force of DNV-RP-F110.'),
        kpi('Lbuckle', 'Buckle length at the minimum', gov.L, 'm', 'Length of the central lobe of the governing mode at its lowest force.'),
        kpi('amplitude', 'Post-buckling lateral amplitude', post ? post.amplitude : 0, 'm', 'Sideways deflection of the buckle crown at the operating force. Zero when the line is not susceptible.'),
        kpi('strain', 'Post-buckling bending strain', post ? post.strain * 100 : 0, '%', 'Bending strain at the buckle crown before the strain concentration factor.'),
        kpi('ucStrain', 'Strain utilisation in the buckle', dcc ? dcc.uc : 0, '-', 'Design strain in the crown over the factored critical strain. Above 1.0 an uncontrolled buckle would wrinkle the wall.', 1, 'max'),
        kpi('feedIn', 'Axial feed-in to the buckle', sp ? sp.feedIn : 0, 'm', 'Length of pipe that slides into each buckle from both sides; the quantity that buckle management limits.'),
        kpi('slip', 'Slip length each side of an isolated buckle', sp ? sp.slipLength : 0, 'm', 'Distance over which friction rebuilds the force after a buckle; buckles closer than twice this share the feed-in.'),
        kpi('ucProp', 'Propagation utilisation', prop.uc, '-', 'External pressure over the factored propagation pressure. Above 1.0 a collapse would run and arrestors are required.', 1, 'max'),
      ],
      checks,
      plots: [
        plot('line', 'Hobbs buckling curves', 'Buckle length', 'm', 'Pre-buckle compressive force P0', 'kN', [
          ...curves, { name: 'Operating restrained force', x: [round(Ls[0]), round(Ls[Ls.length - 1])], y: [round(P0 / 1e3), round(P0 / 1e3)], dash: true },
        ], 'Far-field compressive force needed to hold each buckle mode at a given length. If the dashed operating force lies above the minimum of a curve, that buckle can exist; the right-hand crossing is the stable buckled state.'),
        bar('Critical temperature by mode', [...MODES.map((m) => `Lateral ${m}`), 'Vertical'], 'Operating temperature', '°C',
          [{ name: 'Temperature at the mode minimum', y: [...mins.map((m) => round(inp.Tinst + m.deltaT)), round(TcritV)] }, { name: 'Operating temperature', y: [...MODES.map(() => inp.Top), inp.Top] }],
          'Operating temperature at which each mode first becomes possible, against the actual operating temperature. Modes whose bar is lower than the operating temperature can occur.'),
      ],
      tables: [tbl('Hobbs modes at their minimum', [col('mode', 'Mode'), col('L', 'Buckle length', 'm'), col('P0', 'Safe force', 'kN'), col('T', 'Critical temperature', '°C'), col('amp', 'Amplitude', 'm'), col('eps', 'Bending strain', '%')],
        [...mins.map((m) => ({ mode: `Lateral ${m.mode}`, L: round(m.L), P0: round(m.P0min / 1e3), T: round(inp.Tinst + m.deltaT), amp: round(m.amplitude), eps: round((m.strain ?? 0) * 100) })),
          { mode: 'Vertical', L: round(vert.L), P0: round(vert.P0min / 1e3), T: round(TcritV), amp: round(vert.amplitude), eps: round((vert.strain ?? 0) * 100) }],
        'One row per buckle mode: the buckle length, far-field force, operating temperature, amplitude and bending strain at the lowest point of its curve. The lowest safe force governs.')],
      notes: [
        `Restrained force ${fmt(P0 / 1e3)} kN compression (thermal ${fmt(-r.thermalTerm / 1e3)} kN, pressure ${fmt(-r.pressureTerm / 1e3)} kN, lay tension ${fmt(inp.Hlay)} kN).`,
        ucLat > 1 ? `Lateral mode ${gov.mode} is possible above ${fmt(Tcrit)} °C; at operating conditions the buckle is ${fmt(post ? post.L : 0)} m long with ${fmt(post ? post.amplitude : 0)} m amplitude.` : `The line is below the lowest lateral buckling force; first buckling at about ${fmt(Tcrit)} °C.`,
        'Hobbs curves assume a straight pipe: real lines buckle at imperfections at forces between the mode minimum and the perfect-pipe value.',
      ],
      advice: [
        adv('engineering', status, ucLat > 1 ? `The pipeline is susceptible to lateral buckling (force ratio ${fmt(ucLat)}); the uncontrolled buckle has ${fmt((post?.strain ?? 0) * 100)} % bending strain and a strain utilisation of ${fmt(dcc ? dcc.uc : 0)}.` : `The restrained force is ${fmt(ucLat * 100)} % of the lowest lateral buckling force: no global buckling is expected on an even seabed.`,
          ucLat > 1 ? (dcc && dcc.uc > 1 ? `Control the buckling: plan triggers at no more than ${fmt(sp?.maxSpacing ?? inp.spacing)} m spacing so each buckle takes limited feed-in, or reduce the operating temperature below ${fmt(Tcrit)} °C.` : 'The buckle is acceptable on strain; confirm with an imperfection-based finite-element analysis and check fatigue of the buckle under shutdown cycles.') : `Buckling starts near ${fmt(Tcrit)} °C; keep the margin when the design temperature or pressure is revised.`,
          'DNV-RP-F110 (2018) §6–7; Hobbs (1984)'),
        adv('economic', 'info', inp.spacing > 0 ? `Planned triggers at ${fmt(inp.spacing)} m mean about ${fmt(nTrig)} sleepers or buoyancy sections per 10 km of line.` : 'No buckle triggers are planned.',
          ucLat > 1 ? `Each sleeper adds roughly half a vessel day of pre-lay work; ${fmt(nTrig)} per 10 km is about ${fmt(0.5 * nTrig)} vessel days. A wider spacing is the saving if strain permits.` : 'No buckle-management hardware is needed, which removes pre-lay sleepers and their vessel time from the installation campaign.',
          ECON_BASIS, `${vesselCo2(0.5 * Math.max(nTrig, 1))} Avoiding rock-dump restraint in favour of controlled buckling also keeps the seabed footprint to the sleepers alone.`),
      ],
      publish: { 'pipeline.effectiveForce': r.S0, 'pipeline.buckleSafeForce': gov.P0min },
      quantities: { vessel_days: { sleepers: round(ucLat > 1 ? 0.5 * nTrig : 0, 3) } },
    };
  },
  verify: [{
    // Restrained force: S0 = H − Δp A_i (1 − 2ν) − A_s E α ΔT with A_i = 0.068488 m², A_s = 0.013909 m² (323.9 × 14.3 mm),
    // Δp = 15 MPa, ΔT = 75 K, α = 1.17e-5, E = 207 GPa, H = 0: −410.9 kN − 2526.5 kN = −2937.4 kN.
    // Infinite mode: S∞ = 2.2889 (EA·EI)^¼ f_L^½, EA = 2.8792e9 N, EI = 3.4571e7 N·m², f_L = 0.6·1000 N/m → 995.9 kN.
    name: 'Fully restrained effective force and Hobbs infinite-mode force',
    source: 'DNV-RP-F110 (2018) Eq. 6; Hobbs (1984) mode ∞ minimised analytically; hand calculation',
    inputs: { D: 323.9, t: 14.3, pi: 150, Tinst: 5, Top: 80, Hlay: 0, wsGiven: 1000, muL: 0.6 },
    expect: [{ kpi: 'S0', value: -2937.4, tol: 0.003 }, { kpi: 'Sinf', value: 995.9, tol: 0.01 }],
  }],
};
