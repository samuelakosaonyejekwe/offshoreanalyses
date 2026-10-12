// Blast response of a deck beam, blast-wall stiffener or secondary member: equivalent single-degree-of-freedom
// system with elastic–perfectly-plastic resistance under a given overpressure–time history, with strain-rate
// enhancement, ductility, support rotation, plastic strain, dynamic reactions and the pressure–impulse diagram.
import { iSection } from '../../engine/structural/sections.js';
import { beamBlastResponse, equivalentSdof, pressureImpulseCurve, meanStrainRate, dynamicIncreaseFactor, blastDragLoad, loadHistory } from '../../engine/structural/blast.js';
import { johnsonCook } from '../../engine/structural/materials-nl.js';
import { steelOf, need, run, finite, fmt, round, adv, check, thin, steelCo2Text, CONFIRM, NOT_FE } from './_common.js';

const BAR = 1e5;

export default {
  id: 's02-blast-response', suite: 2, title: 'Blast response (equivalent SDOF and pressure–impulse)',
  summary: 'Peak deflection, ductility, support rotation, plastic strain and dynamic reactions of a beam or stiffener under an explosion overpressure–time history, with its pressure–impulse diagram and the drag load on pipework.',
  covers: [10], fidelity: 'reduced-order',
  standards: ['ISO 19901-3:2014 §8.6 (accidental actions, explosion)', 'NORSOK N-004 (2022) Annex A.6 (explosion response)', 'DNV-RP-C204 (2019) §6 (design against explosions)', 'FABIG Technical Note 7 / API RP 2FB (2006) §C6 (SDOF methods)', 'Biggs (1964) Introduction to Structural Dynamics, ch. 5'],
  models: ['Single-degree-of-freedom blast response (Biggs transformation factors)', 'Elastic–perfectly-plastic resistance (plastic-hinge formulation)', 'Pressure–impulse diagram', 'Cowper–Symonds strain-rate enhancement', 'Johnson–Cook rate term', 'Transient dynamic analysis (explicit time integration)', 'Morison-type blast drag'],
  limits: `Equivalent SDOF of one prismatic member in bending with the overpressure–time history as an input: the pressure is not derived from gas-cloud size or scaled distance, and no CFD explosion load is computed. Membrane action, local buckling of the compression flange, connection failure modes other than end shear, and interaction with adjacent members are not modelled. ${NOT_FE}`,
  searchTerms: 'explosion blast response SDOF Biggs pressure impulse diagram offshore blast wall ductility',
  inputs: [
    { key: 'support', label: 'End conditions', type: 'select', default: 'simple', group: 'Member', role: 'bc', options: [{ value: 'simple', label: 'Simply supported' }, { value: 'fixed', label: 'Fixed at both ends' }, { value: 'cantilever', label: 'Cantilever' }], help: 'Fixed ends need connections able to develop the plastic moment of the member.' },
    { key: 'L', label: 'Span', unit: 'm', type: 'number', default: 6, min: 0.5, max: 40, group: 'Member', role: 'geometry', design: true },
    { key: 's', label: 'Loaded width (member spacing)', unit: 'm', type: 'number', default: 2.0, min: 0.1, max: 20, group: 'Member', role: 'geometry', design: true, help: 'Width of blast wall or deck plating that delivers its pressure to this member.' },
    { key: 'h', label: 'Section depth', unit: 'mm', type: 'number', default: 400, min: 80, max: 3000, group: 'Member', role: 'geometry', design: true },
    { key: 'b', label: 'Flange width', unit: 'mm', type: 'number', default: 180, min: 40, max: 1200, group: 'Member', role: 'geometry' },
    { key: 'tw', label: 'Web thickness', unit: 'mm', type: 'number', default: 10, min: 3, max: 80, group: 'Member', role: 'geometry' },
    { key: 'tf', label: 'Flange thickness', unit: 'mm', type: 'number', default: 16, min: 4, max: 120, group: 'Member', role: 'geometry', design: true },
    { key: 'mArea', label: 'Mass of plating and attachments', unit: 'kg/m²', type: 'number', default: 100, min: 0, max: 2000, group: 'Member', role: 'ic', help: 'Cladding, plate, insulation and supported items moving with the member (they add mass but no strength).' },
    { key: 'mat', label: 'Steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Member', role: 'material' },
    { key: 'sig0', label: 'Bending stress before the explosion', unit: 'MPa', type: 'number', default: 40, min: 0, max: 400, group: 'Member', role: 'ic', help: 'Stress from permanent loads at the hinge location: it uses up part of the plastic moment available for the blast.' },
    { key: 'pt', label: 'Overpressure–time history', type: 'series', xLabel: 'Time', xUnit: 'ms', yLabel: 'Overpressure', yUnit: 'bar', default: { x: [0, 40, 120], y: [0, 0.7, 0] }, group: 'Explosion load', role: 'bc',
      help: 'Design overpressure on the loaded surface from the explosion study (CFD or the project design accidental load). Gas explosions rise and decay over 50–200 ms; the history is used as given.' },
    { key: 'pScale', label: 'Scale factor on the overpressure', unit: '-', type: 'number', default: 1, min: 0.05, max: 10, group: 'Explosion load', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.25 }, help: 'Multiplies every pressure ordinate: use it for sensitivity studies and exceedance levels.' },
    { key: 'qDyn', label: 'Dynamic (drag) pressure of the blast wind', unit: 'kPa', type: 'number', default: 15, min: 0, max: 500, group: 'Explosion load', role: 'bc', help: '½ρu² behind the front, from the explosion study; it loads pipework and small equipment in drag.' },
    { key: 'dPipe', label: 'Diameter of drag-loaded pipe', unit: 'mm', type: 'number', default: 324, min: 10, max: 3000, group: 'Explosion load', role: 'geometry' },
    { key: 'Cd', label: 'Drag coefficient of the pipe', unit: '-', type: 'number', default: 1.2, min: 0.3, max: 3, group: 'Explosion load', role: 'load' },
    { key: 'rate', label: 'Strain-rate enhancement of yield', type: 'select', default: 'cs', group: 'Material rate effects', role: 'material', options: [{ value: 'none', label: 'None (static yield)' }, { value: 'cs', label: 'Cowper–Symonds' }, { value: 'jc', label: 'Johnson–Cook rate term' }] },
    { key: 'csD', label: 'Cowper–Symonds D', unit: '1/s', type: 'number', default: 40.4, min: 1, max: 1e5, group: 'Material rate effects', role: 'material', showIf: { key: 'rate', in: ['cs'] }, help: 'Commonly quoted for mild steel yield: D = 40.4 1/s, q = 5; high-tensile steel 3200 1/s. Material data — confirm for the grade.' },
    { key: 'csQ', label: 'Cowper–Symonds q', unit: '-', type: 'number', default: 5, min: 1, max: 20, group: 'Material rate effects', role: 'material', showIf: { key: 'rate', in: ['cs'] } },
    { key: 'jcC', label: 'Johnson–Cook C', unit: '-', type: 'number', default: 0.022, min: 0, max: 0.2, group: 'Material rate effects', role: 'material', showIf: { key: 'rate', in: ['jc'] }, help: 'Rate constant for a reference strain rate of 1e-3 1/s. Material test data — values of 0.01–0.03 are reported for structural steels.' },
    { key: 'muLim', label: 'Ductility limit', unit: '-', type: 'number', default: 5, min: 1, max: 40, group: 'Acceptance criteria', role: 'criteria', help: `Allowed ratio of peak to yield deflection. Depends on section class, restraint and the function to be kept (1 for elastic response; 5–10 for compact restrained beams). ${CONFIRM}` },
    { key: 'rotLim', label: 'Support rotation limit', unit: 'deg', type: 'number', default: 2, min: 0.2, max: 12, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'epsLim', label: 'Plastic strain limit', unit: '%', type: 'number', default: 5, min: 0.2, max: 20, group: 'Acceptance criteria', role: 'criteria', help: `Critical strain at the hinge for the grade and weld quality (NORSOK N-004 Table A.3-4 / DNV-RP-C208). ${CONFIRM}` },
    { key: 'nStep', label: 'Time steps per natural period', type: 'int', default: 400, min: 50, max: 6400, group: 'Solver', role: 'control', help: 'Explicit integration step = period / value (and never coarser than 1/400 of the load duration).' },
  ],
  mesh: { key: 'nStep', label: 'Time steps per natural period', levels: [100, 200, 400, 800], qoi: 'uMax' },
  compute(inp, ctx) {
    // SI: mm → m, ms → s, bar → Pa, kPa → Pa, MPa → Pa.
    const st = steelOf(ctx, inp.mat);
    need(2 * inp.tf < inp.h && inp.tw < inp.b, 'The section is not a valid I-shape: check the flange and web thicknesses against the depth and width.');
    const sec = run(() => iSection({ h: inp.h / 1e3, b: inp.b / 1e3, tw: inp.tw / 1e3, tf: inp.tf / 1e3, welded: true }), 'Section');
    const tx = (inp.pt?.x || []).map(Number), py = (inp.pt?.y || []).map(Number);
    need(tx.length >= 2 && tx.length === py.length, 'The overpressure history needs at least two (time, pressure) points.');
    for (let i = 1; i < tx.length; i++) need(tx[i] > tx[i - 1], 'The times of the overpressure history must increase.');
    need(tx.every(Number.isFinite) && py.every(Number.isFinite), 'The overpressure history contains a non-numeric value.');
    const L = inp.L, t = tx.map((v) => (v - tx[0]) / 1e3), area = inp.s * L, F = py.map((p) => p * BAR * inp.pScale * area);
    const pPeak = Math.max(...py.map(Math.abs)) * BAR * inp.pScale;
    need(pPeak > 0, 'The overpressure history is zero everywhere.');
    const hist = loadHistory({ type: 'table', t, F }), td = hist.duration, impulse = hist.impulse / area;
    const mass = (sec.A * st.rho + inp.mArea * inp.s) * L, EI = st.E * sec.Iy;
    const avail = 1 - (inp.sig0 * 1e6) / st.fy;
    need(avail > 0.05, 'The stress before the explosion uses up the yield strength: nothing is left to resist the blast.');
    // Dynamic yield: iterate the mean strain rate fy/(E tE) with the time to first yield of the previous pass.
    let dif = 1, resp = null, rateEps = 0;
    const mpOf = (d, mpFac) => sec.Wply * st.fy * d * avail * mpFac;
    // The period does not depend on the strength, so one evaluation fixes the time step: period / nStep, refined to
    // resolve short pulses (duration / 400) but never finer than period / (10 nStep) to bound the run time.
    const T0 = run(() => equivalentSdof({ support: inp.support, L, EI, Mp: mpOf(1, 1), mass }), 'Equivalent system').period;
    const dt = Math.max(Math.min(T0 / inp.nStep, td / 400), T0 / (10 * inp.nStep));
    const solveDt = (d, mpFac = 1) => run(() => beamBlastResponse({ support: inp.support, L, EI, Mp: mpOf(d, mpFac), mass, history: { type: 'table', t, F }, dt }), 'Blast response');
    for (let it = 0; it < 4; it++) {
      resp = solveDt(dif);
      if (inp.rate === 'none') break;
      const R = resp.response, ue = R.ue;
      let tE = R.tMax;
      for (let i = 0; i < R.u.length; i++) if (Math.abs(R.u[i]) >= ue) { tE = R.t[i]; break; }
      rateEps = meanStrainRate({ fy: st.fy * dif, E: st.E, timeToYield: Math.max(tE, 1e-5) });
      const next = inp.rate === 'cs' ? dynamicIncreaseFactor({ strainRate: rateEps, D: inp.csD, q: inp.csQ })
        : johnsonCook({ epsP: 0, strainRate: Math.max(rateEps, 1e-3), A: st.fy, B: 0, n: 1, C: inp.jcC, refRate: 1e-3 }).rateFactor;
      if (Math.abs(next - dif) < 1e-3) { dif = next; break; }
      dif = next;
    }
    resp = solveDt(dif);
    const R = resp.response, sd = resp.sdof, T = sd.period;
    const uMax = finite(resp.uMax, 'Peak deflection'), mu = resp.ductility, rot = resp.rotationDeg;
    const up = Math.max(uMax - R.ue, 0);
    // Hinge strain estimate: plastic hinge rotation over a hinge length equal to the section depth, extreme fibre at h/2.
    const thetaHinge = inp.support === 'cantilever' ? up / L : 2 * up / (L / 2), epsP = thetaHinge / 2;
    const Vpl = (inp.h / 1e3) * (inp.tw / 1e3) * st.fy * dif / Math.sqrt(3);
    const permanent = Math.abs(R.plasticDisplacement || 0);
    const drag = blastDragLoad({ dynamicPressure: inp.qDyn * 1e3 * inp.pScale, Cd: inp.Cd, D: inp.dPipe / 1e3, length: 1 });
    // Section-strength factor that would bring the ductility to its limit (bisection on the plastic moment).
    let fReq = 1;
    if (mu > inp.muLim) { let lo = 1, hi = 2; while (solveDt(dif, hi).ductility > inp.muLim && hi < 64) hi *= 2; for (let i = 0; i < 9; i++) { const mid = 0.5 * (lo + hi); if (solveDt(dif, mid).ductility > inp.muLim) lo = mid; else hi = mid; } fReq = hi; }
    const klm = sd.KLM;
    const pi = run(() => pressureImpulseCurve({ mass, k: sd.k, Rm: sd.Rm, klm, ductility: inp.muLim, shape: 'triangular', nPoints: 16, durationRatio: [0.05, 20] }), 'Pressure–impulse curve');
    const piEl = run(() => pressureImpulseCurve({ mass, k: sd.k, Rm: sd.Rm, klm, ductility: 1, shape: 'triangular', nPoints: 16, durationRatio: [0.05, 20] }), 'Pressure–impulse curve');
    const checks = [
      check('Ductility demand', mu, inp.muLim, 'μ = u_max / u_yield', 'NORSOK N-004 Annex A.6.9; DNV-RP-C204 §6 (ductility limits by section class)', 'Peak deflection over the deflection at first full yield. Above the limit the hinge loses moment capacity through local buckling or tearing and the member no longer carries the blast or the permanent loads.'),
      check('Support rotation', rot, inp.rotLim, 'deg', 'FABIG TN 7 / UFC 3-340-02 (response limits)', 'Angle of the member at its supports at peak deflection. Large rotations strain the end connections and the attached plating beyond what a bending model captures.'),
      check('Plastic strain at the hinge', epsP * 100, inp.epsLim, '%', 'NORSOK N-004 Table A.3-4; DNV-RP-C208 §5 (critical strain)', 'Extreme-fibre plastic strain estimated from the hinge rotation over a hinge length of one section depth. Above the critical strain the flange or its weld tears.'),
      check('End shear from the dynamic reaction', resp.reactionMax / 1e3, Vpl / 1e3, 'kN', 'Biggs §5.5 dynamic reaction V = aR + bF against the plastic shear resistance of the web (EN 1993-1-1 §6.2.6, dynamic yield)', 'The reaction during the blast must be carried by the web and the end connection; shear failure is brittle and prevents the ductile bending response this analysis assumes.'),
    ];
    const gov = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = gov.uc <= 1;
    const [ht, hu, hR, hF] = thin([R.t, R.u, R.R, R.F], 400);
    const mBeam = sec.A * st.rho * L / 1e3;
    return {
      kpis: [
        { key: 'uMax', label: 'Peak deflection', unit: 'mm', value: uMax * 1e3, explain: 'Largest mid-span (or tip) deflection during the blast. Compare with clearances to equipment, piping and escape routes behind the member.' },
        { key: 'mu', label: 'Ductility demand', unit: '-', value: mu, limit: inp.muLim, sense: 'max', explain: 'Peak deflection over the yield deflection. Below 1 the member stays elastic; above 1 it deforms permanently.' },
        { key: 'rot', label: 'Support rotation', unit: 'deg', value: rot, limit: inp.rotLim, sense: 'max', explain: 'Chord rotation at the supports at peak deflection: the usual damage measure with ductility.' },
        { key: 'epsP', label: 'Plastic strain at the hinge', unit: '%', value: epsP * 100, limit: inp.epsLim, sense: 'max', explain: 'Estimated extreme-fibre plastic strain at the hinge; zero while the response is elastic.' },
        { key: 'permanent', label: 'Permanent deflection', unit: 'mm', value: permanent * 1e3, explain: 'Residual set after the blast: the damage left to repair, and the pre-deformation for any subsequent fire.' },
        { key: 'dlf', label: 'Dynamic load factor', unit: '-', value: resp.dlf, explain: 'Peak deflection over the static deflection under the peak load. Up to 2 for long pulses on an elastic member; below 1 for impulsive loads.' },
        { key: 'T', label: 'Natural period', unit: 'ms', value: T * 1e3, explain: 'Elastic period of the equivalent system. The ratio of load duration to this period decides whether the response is impulsive, dynamic or quasi-static.' },
        { key: 'tdT', label: 'Load duration over period', unit: '-', value: td / T, explain: 'Below about 0.3 only the impulse matters; above about 3 the peak pressure governs and the member feels up to twice the static load.' },
        { key: 'pPeak', label: 'Peak overpressure', unit: 'bar', value: pPeak / BAR, explain: 'Largest ordinate of the applied history including the scale factor.' },
        { key: 'impulse', label: 'Specific impulse', unit: 'bar·ms', value: impulse / BAR * 1e3, explain: 'Area under the overpressure–time history.' },
        { key: 'dif', label: 'Dynamic increase factor on yield', unit: '-', value: dif, explain: `Yield strength enhancement at the mean strain rate of ${fmt(rateEps)} 1/s. Typically 1.1–1.3 for gas explosions.` },
        { key: 'Vmax', label: 'Peak dynamic reaction', unit: 'kN', value: resp.reactionMax / 1e3, explain: 'Largest support reaction during the response: the design force for the end connection and the supporting structure.' },
        { key: 'dragLoad', label: 'Blast drag on the pipe', unit: 'kN/m', value: drag.forcePerLength / 1e3, explain: 'Drag force per metre of pipe from the blast wind, to be applied to pipe supports with its own dynamic factor.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Deflection history', x: { label: 'Time', unit: 'ms' }, y: { label: 'Deflection', unit: 'mm' }, series: [{ name: 'Deflection', x: ht.map((v) => round(v * 1e3)), y: hu.map((v) => round(v * 1e3)) }, { name: 'Yield deflection', x: [0, round(ht[ht.length - 1] * 1e3)], y: [round(R.ue * 1e3), round(R.ue * 1e3)], dash: true }],
          explain: 'Deflection against time. Crossing the dashed yield line means plastic deformation; the level about which the member vibrates afterwards is the permanent set.' },
        { type: 'line', title: 'Load and resistance histories', x: { label: 'Time', unit: 'ms' }, y: { label: 'Force on the member', unit: 'kN' }, series: [{ name: 'Blast load', x: ht.map((v) => round(v * 1e3)), y: hF.map((v) => round(v / 1e3)) }, { name: 'Resistance', x: ht.map((v) => round(v * 1e3)), y: hR.map((v) => round(v / 1e3)) }],
          explain: 'Total blast load on the member and the resistance it mobilises. A flat top on the resistance curve is the plastic plateau; the difference between the curves accelerates the mass.' },
        { type: 'line', title: 'Pressure–impulse diagram', x: { label: 'Specific impulse', unit: 'bar·ms' }, y: { label: 'Peak overpressure', unit: 'bar' },
          series: [{ name: `Ductility ${inp.muLim} (limit)`, x: pi.points.map((p) => round(p.I / area / BAR * 1e3)), y: pi.points.map((p) => round(p.F / area / BAR)) }, { name: 'First yield (μ = 1)', x: piEl.points.map((p) => round(p.I / area / BAR * 1e3)), y: piEl.points.map((p) => round(p.F / area / BAR)), dash: true }, { name: 'Design explosion', x: [round(impulse / BAR * 1e3)], y: [round(pPeak / BAR)], mode: 'markers' }],
          explain: 'Iso-damage curves for a zero-rise triangular pulse: combinations of peak pressure and impulse giving first yield and the ductility limit. A design point below and to the left of a curve is safe for that damage level; the vertical branch is the impulsive asymptote, the horizontal one the quasi-static asymptote.' },
      ],
      tables: [{ title: 'Equivalent system', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
        rows: [{ q: 'Stiffness (bilinear)', v: round(sd.k / 1e6), u: 'kN/mm' }, { q: 'Ultimate resistance', v: round(sd.Rm / 1e3), u: 'kN' }, { q: 'Ultimate resistance as pressure', v: round(sd.Rm / area / BAR), u: 'bar' }, { q: 'Yield deflection', v: round(sd.ue * 1e3), u: 'mm' }, { q: 'Mass', v: round(mass), u: 'kg' },
          { q: 'Load–mass factor (elastic / plastic)', v: `${fmt(klm.elastic)} / ${fmt(klm.plastic)}`, u: '-' }, { q: 'Impulsive asymptote at the ductility limit', v: round(pi.impulsiveAsymptote / area / BAR * 1e3), u: 'bar·ms' }, { q: 'Quasi-static asymptote at the ductility limit', v: round(pi.quasiStaticAsymptote / area / BAR), u: 'bar' }, { q: 'Strain rate', v: round(rateEps), u: '1/s' }, { q: 'Time step', v: round(dt * 1e6), u: 'µs' }],
        explain: 'Properties of the equivalent single-degree-of-freedom system and the asymptotes of its pressure–impulse diagram. Use the ultimate resistance as a pressure for a first comparison with the design overpressure.' }],
      notes: [`${gov.name} governs at ${fmt(gov.uc)}. Response is ${td / T < 0.3 ? 'impulsive' : td / T > 3 ? 'quasi-static' : 'dynamic'} (td/T = ${fmt(td / T)}).`, `Dynamic yield ${fmt(st.fy * dif / 1e6)} MPa (factor ${fmt(dif)}); ${fmt(avail * 100)} % of the plastic moment is available after the permanent-load stress.`,
        'The overpressure history is an input: confirm it against the explosion risk analysis for the module and its exceedance frequency.'],
      advice: [
        ok ? adv('engineering', gov.uc > 0.85 ? 'warn' : 'ok', `The member survives the design explosion: ductility ${fmt(mu)} (limit ${inp.muLim}), support rotation ${fmt(rot)}°, peak deflection ${fmt(uMax * 1e3)} mm, permanent set ${fmt(permanent * 1e3)} mm.`,
          `Design the end connections for the dynamic reaction of ${fmt(resp.reactionMax / 1e3)} kN and keep ${fmt(uMax * 1e3 * 1.2, 2)} mm clear behind the member. Confirm the section stays compact (class 1) at the ductility used and check the pipe supports for ${fmt(drag.forcePerLength / 1e3)} kN/m drag.`,
          'NORSOK N-004 (2022) Annex A.6; DNV-RP-C204 (2019) §6; API RP 2FB (2006)', mu > 1 ? `Allowing ductile response instead of elastic design is what keeps the member light: an elastic design for the same pulse would need roughly ${fmt(mu, 2)} times the section strength, in the order of ${steelCo2Text(mBeam * (Math.sqrt(mu) - 1))}` : `The member stays elastic under this explosion and holds ${steelCo2Text(mBeam)} A lighter section working at a ductility of 2–3 would do the same duty with less steel where permanent set is acceptable.`)
          : adv('engineering', 'fail', `${gov.name} is exceeded (utilisation ${fmt(gov.uc)}): ductility ${fmt(mu)} against ${inp.muLim}, rotation ${fmt(rot)}° against ${inp.rotLim}°, plastic strain ${fmt(epsP * 100)} %.`,
            gov.name.startsWith('End shear') ? `Thicken the web at the ends to about ${fmt(inp.tw * gov.uc, 2)} mm or add end stiffeners and a connection rated for ${fmt(resp.reactionMax / 1e3)} kN.` : `Increase the plastic section modulus by a factor of about ${fmt(fReq, 3)} (deeper section or thicker flanges), reduce the spacing to about ${fmt(inp.s / Math.max(fReq, 1), 2)} m, fix the ends, or reduce the load by venting or relocating the member.`,
            'NORSOK N-004 (2022) Annex A.6; DNV-RP-C204 (2019) §6',
            `Raising the section strength by ${fmt(fReq, 3)} through depth costs about ${steelCo2Text(mBeam * Math.max(Math.sqrt(fReq) - 1, 0.05))} per member; venting area or a relief panel lowers the load for every member at once.`),
        adv('economic', 'info', `This member weighs ${fmt(mBeam)} t and carries ${fmt(area)} m² of blast-loaded surface (${fmt(mBeam * 1e3 / area)} kg of primary steel per m²); permanent set after the event is ${fmt(permanent * 1e3)} mm.`,
          permanent > 0 ? 'Decide with the operator whether the member is sacrificial: a permanent set means replacement after the event, which is usually acceptable for a rare accidental load and much cheaper than elastic design.' : 'The member stays elastic and is reusable after the event. If it is not on an escape route or supporting safety-critical items, accepting ductile response would save steel.',
          'ISO 19901-3:2014 §8.6 and ISO 13702:2015 (risk-based accidental load levels); AACE RP 18R-97'),
      ],
      quantities: { steel_t: mBeam },
    };
  },
  verify: [
    { // Long rectangular pulse on an elastic system: the peak deflection is twice the static value (Biggs Fig. 2.7).
      name: 'Suddenly applied constant load on an elastic member: DLF = 2', source: 'Biggs (1964) Introduction to Structural Dynamics, §2.3 (Fig. 2.7)',
      inputs: { pt: { x: [0, 2000], y: [0.02, 0.02] }, rate: 'none', sig0: 0, mArea: 0 }, expect: [{ kpi: 'dlf', value: 2, tol: 0.01 }, { kpi: 'mu', value: 0, abs: 0.999 }] },
    { // Bare beam 400 × 180 × 10 × 16, L = 6 m, simply supported: I = (180·400³ − 170·368³)/12 = 2.5399e8 mm⁴,
      // A = 9440 mm² → m = 74.10 kg/m. Exact first mode ω = π²√(EI/(m L⁴)) = 9.8696·√(210e9·2.5399e-4/(74.10·1296)) = 232.6 rad/s
      // → T = 27.01 ms. The equivalent SDOF (K_LM = 0.78, k = 384EI/5L³) is within 0.5 % of this.
      name: 'Natural period of a simply supported beam: T = (2/π)·L²·√(m/EI)', source: 'Blevins, Formulas for Natural Frequency and Mode Shape, Table 8-1',
      inputs: { pt: { x: [0, 2000], y: [0.02, 0.02] }, rate: 'none', sig0: 0, mArea: 0 }, expect: [{ kpi: 'T', value: 27.01, tol: 0.01 }] },
  ],
};
