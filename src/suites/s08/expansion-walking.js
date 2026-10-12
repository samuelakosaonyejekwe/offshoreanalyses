// Pipeline walking: axial ratcheting of a short pipeline over start-up / shutdown cycles from a seabed slope, an end
// tension (riser) or the travelling thermal front, by a numerical 1-D friction model and the closed-form solutions.
import { walkingSlope, walkingEndTension, walkingArrestForce, thermalTransientSteps, uniformCycleSteps } from '../../engine/pipeline/expansion.js';
import { walkingModel } from './_walking.js';
import { pipeInputs, pipeFrom, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, sev, arr,
  vesselCo2, LIMIT_NO_FE, ECON_BASIS, BAR, DEG } from './_common.js';

const G_O = 'Operating cycle';
const G_R = 'Route and restraint';

export default {
  id: 's08-expansion-walking',
  suite: 8,
  title: 'Pipeline walking and cyclic expansion',
  summary: 'Predicts how far a short hot pipeline creeps along the seabed with every shutdown and restart, and the anchor force that would stop it.',
  covers: [5],
  fidelity: 'numerical',
  standards: ['DNV-RP-F110 (2018) §9 (pipeline walking)', 'Carr, Sinclair & Bruton (2006), OTC 17945 (SAFEBUCK walking solutions)'],
  models: ['Nonlinear pipeline walking model (1-D elastic bar on Coulomb friction, cycled)', 'Elastic–perfectly-plastic axial pipe–soil friction', 'Travelling thermal front transient', 'Closed-form slope and end-tension walking (Carr et al.)', 'Axial restraint / anchor spring'],
  limits: `Walking is computed with a numerical one-dimensional friction model of a straight pipe: no lateral buckles (which interact with walking), no route curve pull-out, rate-independent friction, and an idealised start-up transient with a sharp thermal front and uniform cool-down. The closed forms apply only to fully mobilised (short) pipelines. ${LIMIT_NO_FE}`,
  publishes: ['pipeline.walkPerCycle', 'pipeline.anchorForce'],
  searchTerms: 'pipeline walking axial ratcheting thermal transient SAFEBUCK',
  inputs: [
    ...pipeInputs(),
    num('L', 'Pipeline length', 'm', 2000, 100, 30000, G_R, 'geometry', 'Walking affects short lines that are fully mobilised axially (typically below a few kilometres).'),
    num('slope', 'Seabed slope (falling from hot to cold end)', 'deg', 0.5, -10, 10, G_R, 'geometry', 'Positive when the seabed falls from the inlet towards the outlet; the pipe walks downhill.'),
    num('Tcold', 'Tension applied at the cold end', 'kN', 0, 0, 5000, G_R, 'bc', 'Constant tension from a steel catenary riser or lazy-wave jumper at the outlet; it pulls the pipe towards that end.', { bind: 'riser.tdpTension' }),
    num('Thot', 'Tension applied at the hot end', 'kN', 0, 0, 5000, G_R, 'bc', 'Constant tension at the inlet end.'),
    sel('anchor', 'Hold-back anchor', 'hot', [['none', 'None'], ['hot', 'At the hot end'], ['mid', 'At mid-length'], ['cold', 'At the cold end']], G_R, 'bc', 'A hold-back anchor (suction pile and chain) modelled as a linear spring to ground.'),
    num('kAnchor', 'Anchor stiffness', 'kN/m', 20000, 100, 1e7, G_R, 'bc', 'Stiffness of the anchor and its chain.', { showIf: { key: 'anchor', in: ['hot', 'mid', 'cold'] } }),
    num('Tin', 'Inlet temperature in operation', '°C', 90, 0, 200, G_O, 'bc', 'Hot-end temperature at steady state.', { uncertain: { dist: 'normal', cov: 0.05 } }),
    num('Tamb', 'Ambient (shutdown) temperature', '°C', 5, -2, 35, G_O, 'ic', 'Temperature the whole line returns to on cool-down.', { bind: 'site.Tbed' }),
    num('lambda', 'Temperature decay length', 'km', 15, 0.2, 1000, G_O, 'operation', 'Decay length of the steady temperature profile, ṁc_p/(UπD); a short value gives a steep gradient and more walking.'),
    num('pi', 'Operating pressure', 'bar', 150, 0, 1500, G_O, 'bc', 'Internal pressure applied with the first heating step and removed on shutdown.'),
    sel('transient', 'Heat-up transient', 'front', [['front', 'Thermal front travelling from the hot end'], ['uniform', 'Uniform heating (no gradient)']], G_O, 'ic', 'A travelling front makes the pipe walk towards the cold end even on a flat seabed.'),
    num('muA', 'Axial friction coefficient', '-', 0.5, 0.05, 1.5, 'Pipe–soil friction', 'soil', 'Residual axial friction; walking is larger for lower friction on a slope.', { bind: 'pipeline.axialFriction', uncertain: { dist: 'lognormal', cov: 0.25 } }),
    num('uMob', 'Friction mobilisation displacement', 'mm', 5, 0.1, 100, 'Pipe–soil friction', 'soil', 'Axial displacement at which full friction is reached.'),
    num('nCycles', 'Full shutdown / restart cycles in the design life', '-', 100, 1, 5000, 'Design life', 'operation', 'Number of full thermal cycles; walking accumulates with each.'),
    num('allow', 'Allowable accumulated end displacement', 'm', 1.5, 0.05, 20, 'Acceptance criteria', 'criteria', 'Displacement capacity of the end connections (spool, jumper or riser) over the life.'),
    { key: 'simCycles', label: 'Cycles simulated', type: 'int', default: 3, min: 2, max: 12, group: 'Numerics', role: 'control', help: 'Walking per cycle settles after the first one or two cycles.' },
    { key: 'nElem', label: 'Elements along the pipe', type: 'int', default: 60, min: 20, max: 640, group: 'Numerics', role: 'control', help: 'Number of bar elements in the walking model.' },
  ],
  mesh: { key: 'nElem', label: 'Elements along the pipe', levels: [30, 60, 120], qoi: 'walkRate' },
  compute(inp, ctx) {
    const p = pipeFrom(inp, ctx);
    need(p.ws > 0, 'The pipe is buoyant in operation: it has no axial friction and the walking model does not apply.');
    const L = inp.L, n = inp.nElem, W = p.ws, mu = inp.muA, phi = inp.slope * DEG;
    need(Math.abs(Math.tan(phi)) < mu, `The seabed slope (${fmt(inp.slope)}°) is steeper than friction can hold (tan φ ≥ μ): the pipe slides without any thermal cycling.`);
    need(inp.Tin > inp.Tamb, 'The inlet temperature must be above ambient for a heat-up cycle.');
    const aDT = p.alpha * (inp.Tin - inp.Tamb);
    const ePress = inp.pi * BAR * p.Ai * (1 - 2 * p.nu) / p.EA;
    const lam = inp.lambda * 1e3;
    const steps = inp.transient === 'front' ? thermalTransientSteps({ L, nElem: n, alphaDeltaT: aDT, lambda: lam, pressureStrain: ePress, nHeat: 12, nCool: 4 })
      : uniformCycleSteps({ nElem: n, freeStrain: aDT + ePress, nUp: 4, nDown: 4 });
    const args = { L, EA: p.EA, W, mu, steps, cycles: inp.simCycles, slopeRad: phi, tensionAtEnd: inp.Tcold * 1e3, tensionAtStart: inp.Thot * 1e3, uMob: inp.uMob / 1e3 };
    if (inp.anchor !== 'none') args.anchor = { x: inp.anchor === 'hot' ? 0 : inp.anchor === 'mid' ? L / 2 : L, stiffness: inp.kAnchor * 1e3 };
    need(Math.abs(inp.Tcold - inp.Thot) * 1e3 + W * L * Math.abs(Math.sin(phi)) < mu * W * L * Math.cos(phi), 'The net end tension plus the downslope weight exceeds the total axial friction μ·W·L: the pipe is pulled along the seabed without thermal cycling. Add an anchor or reduce the tension.');
    const m = walkingModel(args);
    const dP = p.EA * (aDT + ePress);
    const cfSlope = walkingSlope({ deltaP: dP, W, L, mu, slopeRad: Math.abs(phi), EA: p.EA });
    const cfTens = inp.Tcold > 0 ? walkingEndTension({ deltaP: dP, SR: inp.Tcold * 1e3, W, L, mu, EA: p.EA }) : { walkPerCycle: 0, fullyMobilised: cfSlope.fullyMobilised };
    const arrest = walkingArrestForce({ W, L, mu, slopeRad: Math.abs(phi), SR: inp.Tcold * 1e3 });
    const rate = m.walkRate;
    const total = Math.abs(rate) * inp.nCycles;
    const anchorF = Math.max(Math.abs(m.anchorForceMin), Math.abs(m.anchorForceMax));
    const xe = Array.from({ length: n }, (_, i) => (i + 0.5) * L / n);
    const cyc = m.walkPerCycle.map((_, i) => i + 1);
    const lifeCycles = Math.abs(rate) > 1e-9 ? inp.allow / Math.abs(rate) : 1e9;
    const uc = total / inp.allow;
    const anchored = inp.anchor !== 'none';
    const hotRange = Math.max(...Array.from(m.hotForce, (v) => -v));
    return {
      kpis: [
        kpi('walkRate', 'Walk per cycle', rate, 'm', 'Net axial displacement of the whole pipeline in each full shutdown and restart cycle (positive towards the cold end), from the numerical model after the start-up cycles.'),
        kpi('walkLife', 'Accumulated walk over the design life', total, 'm', 'Walk per cycle times the number of full cycles. The end connections must tolerate it or an anchor is needed.', inp.allow, 'max'),
        kpi('cyclesAllow', 'Cycles to reach the allowable displacement', Math.min(lifeCycles, 1e9), '-', 'Number of full thermal cycles before the allowable end displacement is used up.'),
        kpi('walkSlopeCF', 'Closed-form walk per cycle, slope', cfSlope.walkPerCycle, 'm', 'Carr et al. solution for a fully mobilised pipeline on a slope with no thermal gradient, for comparison.'),
        kpi('walkTensionCF', 'Closed-form walk per cycle, end tension', cfTens.walkPerCycle, 'm', 'Carr et al. solution for a constant tension at the cold end alone.'),
        kpi('holdBack', 'Hold-back tension to arrest walking', arrest.holdBackTension / 1e3, 'kN', 'Constant tension at the upper end that makes the heat-up and cool-down anchors coincide, stopping slope and tension walking.'),
        kpi('anchorForce', 'Peak anchor force (modelled anchor)', anchorF / 1e3, 'kN', 'Largest force in the modelled hold-back anchor over the last cycle. Zero when no anchor is modelled.'),
        kpi('rigidAnchor', 'Peak tension in a rigid end anchor', arrest.rigidAnchorTension / 1e3, 'kN', 'Upper bound of the cool-down tension a rigid anchor at the end would see for a fully mobilised line.'),
        kpi('endHot', 'Hot-end displacement after the simulated cycles', m.endDisplacementStart[m.endDisplacementStart.length - 1], 'm', 'Cold-state position of the inlet end relative to as-laid.'),
        kpi('endCold', 'Cold-end displacement after the simulated cycles', m.endDisplacementEnd[m.endDisplacementEnd.length - 1], 'm', 'Cold-state position of the outlet end relative to as-laid.'),
        kpi('Pmax', 'Peak effective compression in the cycle', hotRange / 1e3, 'kN', 'Largest compressive effective force reached during heat-up; compare with the fully restrained force to judge mobilisation.'),
        kpi('dP', 'Fully constrained force change', dP / 1e3, 'kN', 'Change in fully restrained force between shutdown and operation at the inlet.'),
      ],
      checks: [
        chk('Accumulated walking vs end-connection capacity', total, inp.allow, 'm', 'DNV-RP-F110 (2018) §9',
          'Guards against the pipeline creeping so far over its life that spools, jumpers or the riser are overstressed or pulled off their supports.'),
      ],
      plots: [
        plot('line', 'Effective force along the pipe', 'Distance from the hot end', 'm', 'Effective force (tension +)', 'kN', [
          { name: 'Most compressive state (hot)', x: arr(xe), y: arr(m.hotForce, 1e-3) }, { name: 'After cool-down', x: arr(xe), y: arr(m.coldForce, 1e-3) },
        ], 'Force profiles at the hot and cold ends of the cycle. The peaks are the virtual anchors; when the hot and cold anchors sit at different places the pipe walks by the strain between them.'),
        plot('line', 'Walk per cycle', 'Cycle', '-', 'Displacement', 'm', [
          { name: 'Walk in the cycle', x: cyc, y: m.walkPerCycle.map((v) => round(v, 5)) }, { name: 'Mean pipe displacement', x: cyc, y: m.meanDisplacement.map((v) => round(v, 5)) },
        ], 'Net movement in each simulated cycle and the running total. A constant walk per cycle after the first cycle shows steady ratcheting.'),
        plot('line', 'Axial displacement along the pipe (cold state)', 'Distance from the hot end', 'm', 'Displacement', 'm', [
          { name: 'After the simulated cycles', x: arr(m.x), y: arr(m.u, 1, 5) },
        ], 'Residual axial displacement after the last cool-down. A uniform offset is pure walking; the slope shows locked-in strain.'),
      ],
      tables: [tbl('Walking mechanisms', [col('mech', 'Mechanism'), col('walk', 'Walk per cycle', 'm'), col('note', 'Basis')], [
        { mech: 'Numerical model (all effects)', walk: round(rate, 5), note: `${n} elements, ${inp.simCycles} cycles, ${inp.transient === 'front' ? 'travelling thermal front' : 'uniform heating'}` },
        { mech: 'Seabed slope alone', walk: round(cfSlope.walkPerCycle, 5), note: cfSlope.fullyMobilised ? 'closed form, fully mobilised' : 'closed form not applicable (not fully mobilised)' },
        { mech: 'Cold-end tension alone', walk: round(cfTens.walkPerCycle, 5), note: inp.Tcold > 0 ? 'closed form' : 'no end tension applied' },
      ], 'Walk per cycle from the numerical model against the closed-form contributions, to show which mechanism drives the walking.')],
      notes: [`${cfSlope.fullyMobilised ? 'The line is fully mobilised (short pipeline): walking is possible.' : 'The line is not fully mobilised in the closed-form sense; any walking comes from the thermal transient.'}`,
        `Numerical model converged in ${m.iterations} iterations; friction ${fmt(mu * W)} N/m, mobilisation ${fmt(inp.uMob)} mm.`,
        'Walking is computed with a numerical 1-D friction model with an idealised transient; calibrate the transient against the flow-assurance warm-up profile.'],
      advice: [
        adv('engineering', sev(uc), `The pipe walks ${fmt(Math.abs(rate) * 1e3)} mm per cycle towards the ${rate >= 0 ? 'cold' : 'hot'} end, ${fmt(total)} m over ${fmt(inp.nCycles)} cycles against ${fmt(inp.allow)} m allowed.`,
          uc > 1 ? `${anchored ? 'Stiffen or relocate the anchor' : `Install a hold-back anchor for about ${fmt(Math.max(arrest.holdBackTension, arrest.rigidAnchorTension) / 1e3)} kN at the uphill end`}, or design the end connections for ${fmt(total)} m; the allowable is reached after ${fmt(lifeCycles)} cycles.`
            : anchored ? `The modelled anchor holds the line with a peak force of ${fmt(anchorF / 1e3)} kN; design the pile, chain and clamp for it with the high-estimate friction, and monitor end positions in service.` : 'No anchor is needed for the stated cycles. Re-check with low-estimate friction and with the measured start-up temperature transient, and monitor end positions in service.',
          'DNV-RP-F110 (2018) §9; Carr, Sinclair & Bruton (2006) OTC 17945'),
        adv('economic', 'info', anchored ? `The hold-back anchor must carry about ${fmt(Math.max(anchorF, arrest.holdBackTension) / 1e3)} kN (suction pile, chain and clamp); with it walking uses ${fmt(uc * 100)} % of the end-connection capacity.`
          : uc > 1 ? `Stopping the walking needs a hold-back anchor of roughly ${fmt(arrest.rigidAnchorTension / 1e3)} kN capacity (suction pile, chain and clamp).` : `Walking uses ${fmt(uc * 100)} % of the end-connection capacity over the life, so no anchor pile is required.`,
          anchored || uc > 1 ? 'An anchor pile is about one to two vessel days plus the pile; limiting full shutdowns or accepting periodic spool re-setting are the operational alternatives to compare over the field life. Run the case without the anchor to see whether it can be omitted.' : 'Avoiding the anchor saves the pile, its chain and one to two installation vessel days; keep an inspection allowance to confirm end positions at each ROV survey.',
          ECON_BASIS, `${vesselCo2(1.5)} An anchor pile of about ${fmt(Math.max(arrest.rigidAnchorTension / 1e3 / 50, 5))} t of steel would add its own embodied carbon; operating with fewer full thermal cycles avoids both.`),
      ],
      publish: { 'pipeline.walkPerCycle': rate, 'pipeline.anchorForce': Math.max(arrest.holdBackTension, arrest.rigidAnchorTension) },
    };
  },
  verify: [{
    // Fully mobilised line on a slope, uniform heating, no pressure: Δ = [ΔP + W L sin φ − μ W L cos φ]·L tan φ/(E A μ).
    // E A = 207e9·0.013909 = 2.8792e9 N; ΔP = E A α ΔT = 2.8792e9·1.17e-5·80 = 2.6949 MN; W = 1000 N/m, L = 2000 m, μ = 0.5,
    // tan φ = 0.05 (φ = 2.8624°): bracket = 2.6949e6 + 99 876 − 998 752 = 1.7961e6 N; Δ = 1.7961e6·2000·0.05/(2.8792e9·0.5) = 0.1248 m.
    // The virtual anchors (0.45 L and 0.55 L) fall on element nodes for 60 elements; the small mobilisation displacement
    // approaches the rigid-plastic friction of the closed form.
    name: 'Slope-driven walking of a fully mobilised pipeline',
    source: 'Carr, Sinclair & Bruton (2006), OTC 17945, slope walking solution; hand calculation',
    inputs: { D: 323.9, t: 14.3, wsGiven: 1000, L: 2000, slope: 2.8624, Tcold: 0, Thot: 0, anchor: 'none', Tin: 85, Tamb: 5, pi: 0, transient: 'uniform', muA: 0.5, uMob: 0.2, nElem: 60, simCycles: 3 },
    expect: [{ kpi: 'walkRate', value: 0.1248, tol: 0.01 }, { kpi: 'walkSlopeCF', value: 0.1248, tol: 0.003 }],
  }, {
    // No slope, no end tension and uniform heating: loading is symmetric about mid-length, so the net displacement per cycle is zero.
    name: 'Symmetric loading does not walk',
    source: 'Symmetry (limiting case)',
    inputs: { slope: 0, Tcold: 0, Thot: 0, anchor: 'none', transient: 'uniform' },
    expect: [{ kpi: 'walkRate', value: 0, abs: 1e-6 }],
  }],
};
