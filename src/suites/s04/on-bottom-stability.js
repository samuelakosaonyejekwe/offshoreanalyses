// On-bottom stability of a subsea structure under operating and extreme near-bed wave and current loading.
import { onBottomStability } from '../../engine/geotech/mudmat.js';
import { nearBedWaveCurrent, currentLogNearBed } from '../../engine/metocean/current.js';
import { G, DEG, need, run, fmt, round, cap, adv, fsCheck, worst, concreteCo2Text, waterInputs, CONFIRM, CO2_STEEL } from './_common.js';

const NPH = 72;

function caseLoads(c, p) {
  // Regular design wave H = hFac·Hs at T = Tp; orbital velocity at the load centroid by linear theory.
  const H = p.hFac * c.Hs;
  const Uw = H > 0 ? run(() => nearBedWaveCurrent({ H, T: c.Tp, d: p.depth, zBed: p.zF }), `Wave kinematics (${c.name})`).Uw : 0;
  const Uc = c.Uc > 0 ? currentLogNearBed(p.zF, { vRef: c.Uc, zRef: p.zRef, z0: p.z0 }) : 0;
  const om = 2 * Math.PI / c.Tp;
  const ph = [], FHs = [], FLs = [], FDs = [], FIs = [], fsS = [], fsO = [], fsU = [];
  let min = { s: Infinity, o: Infinity, u: Infinity }, FHmax = 0, FLmax = 0, ballast = 0, areaReq = 0;
  for (let i = 0; i < NPH; i++) {
    const th = 2 * Math.PI * i / NPH;
    const ux = Uc * Math.cos(p.ang) + Uw * Math.cos(th), uy = Uc * Math.sin(p.ang);
    const sp = Math.hypot(ux, uy);
    const kD = 0.5 * p.rho * p.Cd * p.Ap * sp;                 // drag acts along the instantaneous velocity
    const FI = -p.rho * p.Cm * p.volH * om * Uw * Math.sin(th); // inertia along the wave direction
    let FH = Math.hypot(kD * ux + FI, kD * uy);
    let FL = 0.5 * p.rho * p.Cl * p.Aplan * sp * sp;
    if (c.ext) { if (p.FHext > 0) FH = p.FHext; if (p.FLext > 0) FL = p.FLext; }
    const ob = onBottomStability({ W: p.W, FH, FL, arm: p.zF, B: p.B, mu: p.mu, su: p.suInt, A: p.A, Hpassive: p.Hp });
    ph.push(round(th / DEG, 4)); FHs.push(FH / 1e3); FLs.push(FL / 1e3); FDs.push(kD * sp / 1e3); FIs.push(FI / 1e3);
    fsS.push(cap(ob.fsSliding)); fsO.push(cap(ob.fsOverturning)); fsU.push(cap(ob.fsUplift));
    min = { s: Math.min(min.s, ob.fsSliding), o: Math.min(min.o, ob.fsOverturning), u: Math.min(min.u, ob.fsUplift) };
    FHmax = Math.max(FHmax, FH); FLmax = Math.max(FLmax, FL);
    // Extra submerged weight that would restore each required factor at this phase.
    const dS = p.mu > 0 ? (p.fsS * FH - p.suInt * p.A - p.Hp) / p.mu + FL - p.W : 0;
    const dO = 2 * p.fsO * FH * p.zF / p.B + FL - p.W;
    const dU = p.fsU * FL - p.W;
    ballast = Math.max(ballast, dS, dO, dU);
    if (p.mu === 0 && p.suInt > 0) areaReq = Math.max(areaReq, (p.fsS * FH - p.Hp) / p.suInt);
  }
  return { ...c, H, Uw, Uc, KC: Uw * c.Tp / Math.sqrt(p.Ap), ph, FHs, FLs, FDs, FIs, fsS, fsO, fsU, min, FHmax, FLmax, ballast, areaReq };
}

export default {
  id: 's04-on-bottom-stability', suite: 4, title: 'On-bottom stability',
  summary: 'Sliding, overturning and uplift factors of safety of a subsea structure in the operating and extreme near-bed wave and current climate.',
  covers: [1], fidelity: 'analytical',
  standards: ['DNV-RP-C205 (2019) §6.2, §6.7', 'DNV-RP-F109 (2021) §3.4 (near-bed velocities)', 'ISO 19901-4:2016 §7', 'DNV-RP-N103 (2017) §6'],
  models: ['Morison wave-current loading', 'Near-bed logarithmic current profile', 'Linear wave theory', 'Coulomb friction', 'Undrained interface adhesion', 'Rigid-body sliding and overturning equilibrium', 'Archimedes buoyancy'],
  limits: 'Rigid body on a flat seabed with lumped drag, inertia and lift coefficients on the whole structure; regular design wave by linear theory, no wave-current boundary-layer interaction, no diffraction, shielding or CFD-derived loads. Passive skirt resistance is a user value. Cyclic soil degradation and scour are treated in their own modules.',
  publishes: ['subsea.obsFsSliding', 'subsea.ballastRequired'],
  searchTerms: 'subsea structure mudmat on-bottom stability hydrodynamic loading',
  inputs: [
    { key: 'mass', label: 'Dry mass of structure and equipment', unit: 't', type: 'number', default: 120, min: 1, max: 20000, group: 'Weight and buoyancy', role: 'ic', bind: 'subsea.mass', uncertain: { dist: 'normal', cov: 0.03 },
      help: 'As-weighed mass in air including equipment installed at the time considered.' },
    { key: 'vol', label: 'Displaced volume', unit: 'm³', type: 'number', default: 18, min: 0, max: 20000, group: 'Weight and buoyancy', role: 'ic', bind: 'subsea.volume',
      help: 'Steel volume plus closed buoyant volumes. Submerged weight = (mass − ρ·volume)·g.' },
    { key: 'ballast', label: 'Added ballast (submerged weight)', unit: 'kN', type: 'number', default: 0, min: 0, max: 1e5, group: 'Weight and buoyancy', role: 'ic', design: true },
    { key: 'B', label: 'Base width in the load direction', unit: 'm', type: 'number', default: 10, min: 0.5, max: 60, group: 'Geometry', role: 'geometry', design: true, bind: 'subsea.matB' },
    { key: 'L', label: 'Base length across the load', unit: 'm', type: 'number', default: 14, min: 0.5, max: 80, group: 'Geometry', role: 'geometry', design: true, bind: 'subsea.matL' },
    { key: 'Ap', label: 'Projected solid area normal to the flow', unit: 'm²', type: 'number', default: 30, min: 0.1, max: 2000, group: 'Geometry', role: 'geometry',
      help: 'Sum of member, equipment and panel areas seen by the flow (not the envelope area unless the structure is closed).' },
    { key: 'zF', label: 'Height of the load centroid above the base', unit: 'm', type: 'number', default: 1.8, min: 0.05, max: 30, group: 'Geometry', role: 'geometry', bind: 'subsea.hydroArm' },
    { key: 'Aplan', label: 'Plan area exposed to lift', unit: 'm²', type: 'number', default: 60, min: 0, max: 5000, group: 'Geometry', role: 'geometry', help: 'Solid plan area of roof panels, mat plating and equipment.' },
    { key: 'volH', label: 'Hydrodynamic volume for inertia', unit: 'm³', type: 'number', default: 25, min: 0, max: 20000, group: 'Geometry', role: 'geometry', help: 'Displaced volume plus entrapped water moving with the structure.' },
    { key: 'Cd', label: 'Drag coefficient', unit: '-', type: 'number', default: 1.3, min: 0.3, max: 3, group: 'Hydrodynamic coefficients', role: 'load', uncertain: { dist: 'normal', cov: 0.15 },
      help: 'On the projected solid area. Calibrate against flume tests for the actual shape, seabed proximity and shielding.' },
    { key: 'Cm', label: 'Inertia coefficient', unit: '-', type: 'number', default: 2.0, min: 1, max: 4, group: 'Hydrodynamic coefficients', role: 'load', uncertain: { dist: 'normal', cov: 0.1 } },
    { key: 'Cl', label: 'Lift coefficient', unit: '-', type: 'number', default: 0.3, min: 0, max: 2, group: 'Hydrodynamic coefficients', role: 'load', uncertain: { dist: 'normal', cov: 0.25 },
      help: 'On the plan area. Open frames 0–0.3; solid roofs close to the bed can be much higher.' },
    ...waterInputs(),
    { key: 'z0mm', label: 'Seabed roughness z0', unit: 'mm', type: 'number', default: 0.04, min: 0.001, max: 50, group: 'Site', role: 'environment', help: 'DNV-RP-F109 Table 3-1: silt 0.005, fine sand 0.01, medium sand 0.04, coarse sand 0.1, gravel 0.3 mm.' },
    { key: 'zRef', label: 'Reference height of the current', unit: 'm', type: 'number', default: 1, min: 0.1, max: 20, group: 'Site', role: 'environment' },
    { key: 'Hs1', label: 'Operating significant wave height', unit: 'm', type: 'number', default: 6, min: 0, max: 20, group: 'Operating condition', role: 'bc', bind: 'site.Hs1' },
    { key: 'Tp1', label: 'Operating peak period', unit: 's', type: 'number', default: 11, min: 3, max: 25, group: 'Operating condition', role: 'bc', bind: 'site.Tp1' },
    { key: 'Uc1', label: 'Operating near-bed current', unit: 'm/s', type: 'number', default: 0.4, min: 0, max: 3, group: 'Operating condition', role: 'bc' },
    { key: 'Hs100', label: 'Extreme significant wave height', unit: 'm', type: 'number', default: 13, min: 0, max: 25, group: 'Extreme condition', role: 'bc', bind: 'site.Hs100', uncertain: { dist: 'normal', cov: 0.08 } },
    { key: 'Tp100', label: 'Extreme peak period', unit: 's', type: 'number', default: 15.5, min: 3, max: 25, group: 'Extreme condition', role: 'bc', bind: 'site.Tp100' },
    { key: 'Uc100', label: 'Extreme near-bed current', unit: 'm/s', type: 'number', default: 0.7, min: 0, max: 3, group: 'Extreme condition', role: 'bc', bind: 'site.UcBed100', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'hFac', label: 'Design wave height / Hs', unit: '-', type: 'number', default: 1.86, min: 1, max: 2, group: 'Extreme condition', role: 'bc', help: 'Most probable maximum wave in a 3-hour sea state ≈ 1.86 Hs.' },
    { key: 'angDeg', label: 'Angle between current and waves', unit: 'deg', type: 'number', default: 0, min: 0, max: 90, group: 'Extreme condition', role: 'bc' },
    { key: 'FHext', label: 'Extreme horizontal load from the member-level analysis', unit: 'kN', type: 'number', default: 0, min: 0, max: 1e5, group: 'Extreme condition', role: 'bc', bind: 'subsea.hydroFH',
      help: 'If greater than zero it replaces the lumped-coefficient load of the extreme case (all phases).' },
    { key: 'FLext', label: 'Extreme lift from the member-level analysis', unit: 'kN', type: 'number', default: 0, min: 0, max: 1e5, group: 'Extreme condition', role: 'bc', bind: 'subsea.hydroFL' },
    { key: 'mode', label: 'Base resistance', type: 'select', default: 'undrained', group: 'Seabed interface', role: 'soil',
      options: [{ value: 'undrained', label: 'Undrained adhesion (clay)' }, { value: 'drained', label: 'Drained friction (sand)' }] },
    { key: 'su0', label: 'Undrained strength at base level', unit: 'kPa', type: 'number', default: 5, min: 0.2, max: 500, group: 'Seabed interface', role: 'soil', bind: 'site.su0', uncertain: { dist: 'lognormal', cov: 0.25 }, showIf: { key: 'mode', in: ['undrained'] } },
    { key: 'alphaInt', label: 'Interface adhesion factor', unit: '-', type: 'number', default: 0.6, min: 0.05, max: 1, group: 'Seabed interface', role: 'soil', uncertain: { dist: 'normal', cov: 0.2 }, showIf: { key: 'mode', in: ['undrained'] },
      help: 'Ratio of interface to intact strength: about 1/sensitivity for smooth steel, 1.0 for a skirted base shearing in soil.' },
    { key: 'deltaDeg', label: 'Interface friction angle δ', unit: 'deg', type: 'number', default: 25, min: 5, max: 45, group: 'Seabed interface', role: 'soil', uncertain: { dist: 'normal', cov: 0.1 }, showIf: { key: 'mode', in: ['drained'] } },
    { key: 'Hp', label: 'Passive resistance of skirts / embedment', unit: 'kN', type: 'number', default: 0, min: 0, max: 1e5, group: 'Seabed interface', role: 'soil', design: true },
    { key: 'fsSlide', label: 'Required sliding factor of safety', unit: '-', type: 'number', default: 1.5, min: 1, max: 4, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'fsOT', label: 'Required overturning factor of safety', unit: '-', type: 'number', default: 1.5, min: 1, max: 4, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'fsUp', label: 'Required uplift factor of safety', unit: '-', type: 'number', default: 1.2, min: 1, max: 4, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
  ],
  compute(inp) {
    // SI conversion: t → kg, kN → N, kPa → Pa, mm → m, deg → rad.
    const rho = inp.rhoW;
    const W = (inp.mass * 1e3 - rho * inp.vol) * G + inp.ballast * 1e3;
    need(W > 0, `The structure is buoyant: mass ${fmt(inp.mass)} t against ${fmt(rho * inp.vol / 1e3)} t of displaced water. Reduce the displaced volume or add ballast.`);
    need(inp.zF < inp.depth, `The load centroid height (${fmt(inp.zF)} m) must be below the water surface (depth ${fmt(inp.depth)} m).`);
    const und = inp.mode === 'undrained';
    const p = {
      rho, W, depth: inp.depth, zF: inp.zF, zRef: inp.zRef, z0: inp.z0mm / 1e3, ang: inp.angDeg * DEG, hFac: inp.hFac,
      Cd: inp.Cd, Cm: inp.Cm, Cl: inp.Cl, Ap: inp.Ap, Aplan: inp.Aplan, volH: inp.volH, B: inp.B, A: inp.B * inp.L,
      mu: und ? 0 : Math.tan(inp.deltaDeg * DEG), suInt: und ? inp.alphaInt * inp.su0 * 1e3 : 0, Hp: inp.Hp * 1e3,
      fsS: inp.fsSlide, fsO: inp.fsOT, fsU: inp.fsUp, FHext: inp.FHext * 1e3, FLext: inp.FLext * 1e3,
    };
    const cases = [
      caseLoads({ name: 'Operating', Hs: inp.Hs1, Tp: inp.Tp1, Uc: inp.Uc1, ext: false }, p),
      caseLoads({ name: 'Extreme', Hs: inp.Hs100, Tp: inp.Tp100, Uc: inp.Uc100, ext: true }, p),
    ];
    const gov = (k) => cases.reduce((a, c) => (c.min[k] < a.min[k] ? c : a));
    const fsS = Math.min(...cases.map((c) => c.min.s)), fsO = Math.min(...cases.map((c) => c.min.o)), fsU = Math.min(...cases.map((c) => c.min.u));
    const ex = cases[1];
    const ballast = Math.max(...cases.map((c) => c.ballast), 0);
    const areaReq = Math.max(...cases.map((c) => c.areaReq), 0);
    const checks = [
      fsCheck(`Sliding (${gov('s').name.toLowerCase()})`, fsS, inp.fsSlide, 'ISO 19901-4 §7 / DNV-RP-N103 §6 (rigid-body sliding)', 'Guards against the structure being pushed along the seabed. If exceeded the structure walks, overstressing jumpers and connectors.'),
      fsCheck(`Overturning (${gov('o').name.toLowerCase()})`, fsO, inp.fsOT, 'Statics about the downstream edge', 'Guards against rotation about the downstream edge. If exceeded the upstream edge lifts off and the base pressure concentrates under the toe.'),
      fsCheck(`Uplift (${gov('u').name.toLowerCase()})`, fsU, inp.fsUp, 'Statics, vertical equilibrium', 'Guards against hydrodynamic lift exceeding the submerged weight, at which point the structure loses all friction.'),
    ];
    const w = worst(checks);
    const ok = w.uc <= 1;
    // Sensitivity of the sliding factor to the extreme sea state.
    const hs = [], fss = [];
    for (let i = 0; i <= 12; i++) {
      const f = 0.4 + i * 0.1;
      const c = caseLoads({ name: 'sweep', Hs: inp.Hs100 * f, Tp: inp.Tp100, Uc: inp.Uc100, ext: false }, p);
      hs.push(round(inp.Hs100 * f, 4)); fss.push(cap(c.min.s));
    }
    const ballastM3 = ballast / ((2400 - rho) * G);
    const advice = [];
    if (ok) {
      advice.push(adv('engineering', 'ok', `All on-bottom stability factors are met; the governing check is ${w.name.toLowerCase()} at ${fmt(1 / w.uc * w.demand)} against ${fmt(w.demand)} required.`,
        'Confirm the lumped drag, inertia and lift coefficients against the member-level load module or flume data, and re-run with the scoured and cyclically degraded soil state.',
        'DNV-RP-C205 (2019) §6; ISO 19901-4:2016 §7; DNV-RP-N103 (2017) §6',
        `Reserve of ${fmt((1 / w.uc - 1) * 100)} % on the governing check: do not add ballast — each tonne of avoided ballast steel saves about ${fmt(CO2_STEEL)} t CO₂e (screening factor).`));
    } else {
      const fix = p.mu > 0 || w.name.startsWith('Over') || w.name.startsWith('Upl')
        ? `Add about ${fmt(ballast / 1e3)} kN of submerged ballast (≈ ${fmt(ballastM3)} m³ of concrete), or reduce the exposed area.`
        : `Ballast does not raise undrained sliding resistance: increase the bearing area to about ${fmt(areaReq)} m² (now ${fmt(p.A)} m²) or add skirts providing ${fmt((inp.fsSlide * ex.FHmax - p.suInt * p.A - p.Hp) / 1e3)} kN of passive resistance.`;
      advice.push(adv('engineering', 'fail', `${w.name} is not satisfied: factor of safety ${fmt(w.capacity)} against ${fmt(w.demand)} required.`, fix,
        'ISO 19901-4:2016 §7; DNV-RP-N103 (2017) §6; DNV-RP-C205 (2019) §6',
        ballast > 0 ? `Ballast option: ${concreteCo2Text(ballastM3)} Skirts reach the same resistance with far less material.` : 'Skirts mobilise soil that is already there: a few tonnes of plate instead of tens of tonnes of ballast.'));
    }
    advice.push(adv('economic', ok ? 'info' : 'warn',
      ok ? 'No stabilisation measures are needed, so no ballast, rock or additional vessel time enters the estimate.' : (ballast > 0 ? `Stabilisation adds material and offshore time: ${fmt(ballast / 1e3 / G)} t of submerged ballast (or the equivalent in skirts or mat area) has to be fabricated, lifted and landed.` : `Stabilisation adds material: the base area has to grow from ${fmt(p.A)} m² to about ${fmt(areaReq)} m², or skirts have to be added; ballast would not help on an undrained base.`),
      ok ? 'Hold the current weight budget; track weight growth against the stability reserve during detail design.' : 'Compare ballast (heavier lift, possibly a larger crane vessel) with skirts or a larger mat (fabrication only) before fixing the concept; a post-installation fix by rock or mattresses costs roughly one extra vessel day.',
      'AACE RP 18R-97 estimate classes; whole-life costing per ISO 15686-5'));
    return {
      kpis: [
        { key: 'fsSlide', label: 'Sliding factor of safety (minimum)', value: cap(fsS), unit: '-', limit: inp.fsSlide, sense: 'min', explain: 'Horizontal resistance of the base over the peak horizontal hydrodynamic load, minimum over the wave cycle and both sea states. Below the required value the structure may slide. Values are capped at 99.' },
        { key: 'fsOT', label: 'Overturning factor of safety (minimum)', value: cap(fsO), unit: '-', limit: inp.fsOT, sense: 'min', explain: 'Restoring moment of the net submerged weight about the downstream edge over the hydrodynamic overturning moment. Below 1 the structure tips.' },
        { key: 'fsUp', label: 'Uplift factor of safety (minimum)', value: cap(fsU), unit: '-', limit: inp.fsUp, sense: 'min', explain: 'Submerged weight over peak hydrodynamic lift. Close to 1 the base carries no friction and every other check collapses.' },
        { key: 'Wsub', label: 'Submerged weight', unit: 'kN', value: W / 1e3, explain: 'Weight in water by Archimedes: (mass − ρ·displaced volume)·g plus ballast. This is the only stabilising action on a drained base.' },
        { key: 'FHmax', label: 'Peak horizontal load (extreme)', unit: 'kN', value: ex.FHmax / 1e3, explain: 'Largest combined drag and inertia load over the extreme wave cycle. Drag dominates when it peaks under the crest; inertia peaks a quarter cycle earlier.' },
        { key: 'FLmax', label: 'Peak lift (extreme)', unit: 'kN', value: ex.FLmax / 1e3, explain: 'Largest vertical hydrodynamic lift in the extreme case. It unloads the base in phase with the drag load.' },
        { key: 'UbExt', label: 'Peak near-bed velocity (extreme)', unit: 'm/s', value: ex.Uc * Math.cos(p.ang) + ex.Uw, explain: 'Current plus wave orbital velocity at the load centroid. Loads scale with its square, so this is the most influential input.' },
        { key: 'ballastReq', label: 'Additional submerged ballast required', unit: 'kN', value: ballast / 1e3, explain: 'Extra submerged weight that brings every factor up to its required value where weight helps (friction, overturning, uplift). Zero when all checks pass.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Factors of safety through the extreme wave cycle', x: { label: 'Wave phase', unit: 'deg' }, y: { label: 'Factor of safety', unit: '-' },
          series: [{ name: 'Sliding', x: ex.ph, y: ex.fsS }, { name: 'Overturning', x: ex.ph, y: ex.fsO }, { name: 'Uplift', x: ex.ph, y: ex.fsU }, { name: 'Required sliding', x: [0, 360], y: [inp.fsSlide, inp.fsSlide], dash: true }],
          explain: 'Each factor of safety at every phase of the design wave (capped at 99). Look for the trough of each curve and whether it dips below the required line.' },
        { type: 'line', title: 'Hydrodynamic load components through the extreme wave cycle', x: { label: 'Wave phase', unit: 'deg' }, y: { label: 'Force', unit: 'kN' },
          series: [{ name: 'Drag', x: ex.ph, y: ex.FDs }, { name: 'Inertia', x: ex.ph, y: ex.FIs }, { name: 'Resultant horizontal', x: ex.ph, y: ex.FHs }, { name: 'Lift', x: ex.ph, y: ex.FLs }],
          explain: 'Drag, inertia, their resultant and lift over one wave period. Drag and lift peak together under the crest, which is why sliding usually governs there.' },
        { type: 'line', title: 'Sliding factor of safety against extreme wave height', x: { label: 'Hs', unit: 'm' }, y: { label: 'Sliding factor of safety', unit: '-' },
          series: [{ name: 'Sliding', x: hs, y: fss }, { name: 'Required', x: [hs[0], hs[hs.length - 1]], y: [inp.fsSlide, inp.fsSlide], dash: true }],
          explain: 'How fast the sliding reserve erodes as the sea state grows at the extreme period and current. The crossing with the dashed line is the limiting sea state.' },
      ],
      tables: [{ title: 'Sea-state summary', columns: [{ key: 'name', label: 'Case' }, { key: 'H', label: 'Design wave', unit: 'm' }, { key: 'Uw', label: 'Orbital velocity', unit: 'm/s' }, { key: 'Uc', label: 'Current at centroid', unit: 'm/s' }, { key: 'KC', label: 'KC', unit: '-' }, { key: 'FH', label: 'Peak horizontal', unit: 'kN' }, { key: 'FL', label: 'Peak lift', unit: 'kN' }, { key: 'fsS', label: 'FS sliding' }, { key: 'fsO', label: 'FS overturning' }, { key: 'fsU', label: 'FS uplift' }],
        rows: cases.map((c) => ({ name: c.name, H: round(c.H), Uw: round(c.Uw), Uc: round(c.Uc), KC: round(c.KC), FH: round(c.FHmax / 1e3), FL: round(c.FLmax / 1e3), fsS: round(cap(c.min.s)), fsO: round(cap(c.min.o)), fsU: round(cap(c.min.u)) })),
        explain: 'One row per sea state: kinematics at the load centroid, peak loads and the minimum factors over the cycle. KC is based on the square root of the projected area.' }],
      notes: [`${w.name} governs with required/achieved = ${fmt(w.uc)}.`, und ? `Undrained base: resistance = α·su·A = ${fmt(p.suInt * p.A / 1e3)} kN, independent of weight.` : `Drained base: resistance = tan δ × (W′ − lift) with tan δ = ${fmt(p.mu)}.`,
        p.FHext > 0 ? 'The extreme-case horizontal load was taken from the member-level hydrodynamic analysis.' : 'Loads use lumped coefficients on the whole structure.'],
      advice,
      publish: { 'subsea.obsFsSliding': cap(fsS), 'subsea.ballastRequired': ballast },
      quantities: ballast > 0 ? { concrete_m3: round(ballastM3) } : {},
    };
  },
  verify: [
    { name: 'Current only, drained base: FS = tan δ · W′ / (½ρ·Cd·A·U²)', source: 'Hand calculation: W′ = (100 000 − 1025·10)·9.80665 = 880 147 N; FH = ½·1025·1.2·20·1² = 12 300 N; FS = 0.5·880 147/12 300 = 35.78',
      inputs: { mass: 100, vol: 10, Hs1: 0, Hs100: 0, Uc1: 1, Uc100: 1, zRef: 2, zF: 2, Cd: 1.2, Ap: 20, Cl: 0, mode: 'drained', deltaDeg: 26.565051, B: 10, angDeg: 0 },
      expect: [{ kpi: 'fsSlide', value: 35.778, tol: 0.002 }, { kpi: 'Wsub', value: 880.147, tol: 1e-4 }, { kpi: 'FHmax', value: 12.3, tol: 1e-4 }] },
    { name: 'Uplift: FS = W′ / (½ρ·Cl·A·U²)', source: 'Hand calculation: FL = ½·1025·1.0·50·1² = 25 625 N; FS = 880 147/25 625 = 34.35; overturning FS = (880 147 − 25 625)·5/(12 300·2) = 173.7 → capped at 99',
      inputs: { mass: 100, vol: 10, Hs1: 0, Hs100: 0, Uc1: 1, Uc100: 1, zRef: 2, zF: 2, Cd: 1.2, Ap: 20, Cl: 1, Aplan: 50, mode: 'drained', deltaDeg: 26.565051 },
      expect: [{ kpi: 'fsUp', value: 34.347, tol: 0.002 }, { kpi: 'fsSlide', value: 34.737, tol: 0.002 }] },
  ],
};

