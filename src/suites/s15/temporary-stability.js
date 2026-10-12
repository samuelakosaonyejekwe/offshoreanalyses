// Temporary conditions during installation: on-bottom stability of the unpiled jacket and elastic stability of a
// free-standing pile with the hammer on top (frame finite-element buckling with self-weight).
import { pileForceClosedForm } from '../../engine/hydro/morison.js';
import { bucklingAnalysis } from '../../engine/structural/framebuckling.js';
import { G, DEG, need, run, round, fmt, kpi, check, fsCheck, adv, line, bar, tbl, worst, chiEC3, tube, linspace, steelCo2Text, vesselCo2Text, seaInputs, CONFIRM } from './_common.js';

export default {
  id: 's15-temporary-stability', suite: 15, title: 'Temporary stability during installation',
  summary: 'Checks the jacket standing unpiled on its mudmats against sliding, overturning and bearing failure in the installation sea state, and the free-standing pile with the hammer on top against buckling and bending.',
  covers: [17], fidelity: 'numerical',
  standards: ['API RP 2A-WSD (2014) §6.4, §6.13 (mudmat stability of unpiled jackets) and §12.5 (pile stick-up)', 'ISO 19902 (2020) §22 (installation: on-bottom stability, free-standing piles)', 'DNV-ST-N001 (2021) Sec. 13 and 18 (jacket set-down, pile installation)', 'EN 1993-1-1 (2005) §6.3 (buckling resistance, interaction)', 'Timoshenko & Gere, Theory of Elastic Stability §2.12 (column under its own weight)'],
  models: ['Rigid-body equilibrium on four mudmats', 'Morison wave loading (closed-form vertical cylinder, drag and inertia)', 'Undrained bearing capacity and base sliding', 'Linear eigenvalue buckling by beam finite elements with geometric stiffness (self-weight and tip load)', 'Column buckling curve and beam–column interaction', 'Second-order amplification of bending'],
  limits: 'The jacket wave load is that of equivalent vertical cylinders in a linear wave without current blockage or wave spreading, mudmat capacity uses a simple undrained or drained expression without embedment, inclination or eccentricity interaction, and all four mudmats are equally stiff. The pile is a prismatic cantilever fixed at its support level; guide clearances, sleeve flexibility and vortex-induced vibration are not modelled. Safety factors are recalled values exposed as inputs.',
  publishes: ['install.onBottomFS', 'install.stickupLambda'],
  searchTerms: 'unpiled jacket on-bottom stability mudmat free-standing pile stick-up buckling hammer',
  inputs: [
    { key: 'Wob', label: 'On-bottom weight of the jacket', unit: 'kN', type: 'number', default: 30000, min: 10, max: 1e6, group: 'Unpiled jacket', role: 'load', uncertain: { dist: 'normal', cov: 0.05 }, help: 'Weight less buoyancy after set-down and flooding, before piles are stabbed.' },
    { key: 'Bj', label: 'Base width between mudmat centres (in the load direction)', unit: 'm', type: 'number', default: 45, min: 5, max: 150, group: 'Unpiled jacket', role: 'geometry' },
    { key: 'Amat', label: 'Area of one mudmat', unit: 'm²', type: 'number', default: 120, min: 2, max: 2000, group: 'Unpiled jacket', role: 'geometry', design: true },
    { key: 'nLeg', label: 'Equivalent vertical cylinders for wave load', unit: '-', type: 'int', default: 4, min: 1, max: 40, group: 'Unpiled jacket', role: 'geometry' },
    { key: 'Deq', label: 'Equivalent cylinder diameter (legs with braces and appurtenances)', unit: 'm', type: 'number', default: 4.5, min: 0.3, max: 20, group: 'Unpiled jacket', role: 'geometry', help: 'Chosen so that nLeg cylinders of this diameter give the drag and inertia of the jacket; calibrate against the in-place wave-load analysis.' },
    { key: 'Cd', label: 'Drag coefficient', unit: '-', type: 'number', default: 1.05, min: 0.5, max: 2, group: 'Unpiled jacket', role: 'load', uncertain: { dist: 'normal', cov: 0.15 } },
    { key: 'Cm', label: 'Inertia coefficient', unit: '-', type: 'number', default: 1.6, min: 1, max: 2.5, group: 'Unpiled jacket', role: 'load' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 90, min: 5, max: 500, group: 'Site', role: 'environment', bind: 'site.depth' },
    ...seaInputs({ Hs: 2.5, Tp: 8, group: 'Installation sea state' }),
    { key: 'fHmax', label: 'Design wave height / Hs', unit: '-', type: 'number', default: 1.86, min: 1, max: 2.2, group: 'Installation sea state', role: 'criteria', help: `Most probable maximum wave in about 1000 waves; the unpiled condition is usually designed for a seasonal 1-year storm. ${CONFIRM}` },
    { key: 'Fwind', label: 'Wind and current force on the jacket', unit: 'kN', type: 'number', default: 250, min: 0, max: 1e5, group: 'Installation sea state', role: 'bc', help: 'Additional horizontal force acting with the wave.' },
    { key: 'zWind', label: 'Lever of that force above the mudline', unit: 'm', type: 'number', default: 80, min: 0, max: 400, group: 'Installation sea state', role: 'bc' },
    { key: 'soil', label: 'Seabed soil', type: 'select', default: 'clay', group: 'Soil', role: 'soil', options: [{ value: 'clay', label: 'Clay (undrained)' }, { value: 'sand', label: 'Sand (drained)' }] },
    { key: 'su', label: 'Undrained shear strength at the mudmats', unit: 'kPa', type: 'number', default: 25, min: 0.5, max: 500, group: 'Soil', role: 'soil', bind: 'site.su0', uncertain: { dist: 'lognormal', cov: 0.25 }, showIf: { key: 'soil', in: ['clay'] } },
    { key: 'phiDeg', label: 'Friction angle', unit: 'deg', type: 'number', default: 32, min: 20, max: 45, group: 'Soil', role: 'soil', bind: 'site.phiDeg', showIf: { key: 'soil', in: ['sand'] } },
    { key: 'gammaSub', label: 'Submerged unit weight', unit: 'kN/m³', type: 'number', default: 9, min: 3, max: 13, group: 'Soil', role: 'soil', bind: 'site.gammaSub', showIf: { key: 'soil', in: ['sand'] } },
    { key: 'fsBear', label: 'Required factor of safety, bearing', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Criteria', role: 'criteria', help: `API RP 2A: 2.0 for dead plus live loads, 1.5 with the installation environment. ${CONFIRM}` },
    { key: 'fsSlide', label: 'Required factor of safety, sliding and overturning', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Criteria', role: 'criteria', help: CONFIRM },
    { key: 'Dp', label: 'Pile outer diameter', unit: 'mm', type: 'number', default: 2134, min: 300, max: 6000, group: 'Free-standing pile', role: 'geometry' },
    { key: 'tp', label: 'Pile wall thickness', unit: 'mm', type: 'number', default: 60, min: 10, max: 150, group: 'Free-standing pile', role: 'geometry', design: true },
    { key: 'Ls', label: 'Stick-up above the support (sleeve or mudline fixity)', unit: 'm', type: 'number', default: 55, min: 2, max: 150, group: 'Free-standing pile', role: 'bc', design: true, help: 'Boundary condition: the pile is taken as fixed at this level and free at the top where the hammer sits.' },
    { key: 'mHammer', label: 'Hammer, sleeve and cap mass on the pile top', unit: 't', type: 'number', default: 250, min: 0, max: 2000, group: 'Free-standing pile', role: 'load', uncertain: { dist: 'normal', cov: 0.05 } },
    { key: 'batterDeg', label: 'Pile inclination from the vertical', unit: 'deg', type: 'number', default: 5.7, min: 0, max: 20, group: 'Free-standing pile', role: 'ic', help: 'Batter of the jacket leg or allowable out-of-verticality; it makes the hammer weight bend the pile. 1:10 batter = 5.7°.' },
    { key: 'qLat', label: 'Lateral wave and current load on the stick-up', unit: 'kN/m', type: 'number', default: 2, min: 0, max: 200, group: 'Free-standing pile', role: 'bc' },
    { key: 'mat', label: 'Pile steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Free-standing pile', role: 'material' },
    { key: 'gM', label: 'Material factor', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.5, group: 'Criteria', role: 'criteria' },
    { key: 'fLoad', label: 'Load factor on hammer and pile weight', unit: '-', type: 'number', default: 1.3, min: 1, max: 2, group: 'Criteria', role: 'criteria', help: CONFIRM },
    { key: 'nElem', label: 'Beam elements along the stick-up', unit: '-', type: 'int', default: 16, min: 2, max: 200, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Beam elements', levels: [4, 8, 16, 32], qoi: 'lambdaCr' },
  compute(inp, ctx) {
    // kN → N, kPa → Pa, mm → m, t → kg.
    const steel = ctx.material(inp.mat); need(steel, `Steel '${inp.mat}' is not in the materials library.`);
    const rho = 1025, W = inp.Wob * 1e3;
    // --- Unpiled jacket: wave load on equivalent cylinders with the design wave of the installation sea state.
    const Hd = inp.fHmax * inp.Hs, Tw = 0.9 * inp.Tp;
    need(Hd < 0.78 * inp.depth, 'The design wave exceeds the breaking limit for this water depth: reduce the sea state.');
    const wv = inp.Hs > 0 ? run(() => pileForceClosedForm({ H: Hd, T: Tw, depth: inp.depth, D: inp.Deq, Cd: inp.Cd, Cm: inp.Cm, rho }), 'Wave load') : { Fmax: 0, Mmax: 0, FD: 0, FI: 0 };
    const Fh = inp.nLeg * wv.Fmax + inp.Fwind * 1e3, Mo = inp.nLeg * wv.Mmax + inp.Fwind * 1e3 * inp.zWind;
    const clay = inp.soil === 'clay', phi = inp.phiDeg * DEG, A = inp.Amat, Bm = Math.sqrt(A);
    // Two mudmats on each side share the overturning couple.
    const Vmax = W / 4 + Mo / (2 * inp.Bj), Vmin = W / 4 - Mo / (2 * inp.Bj);
    let qult, Hres;
    if (clay) { qult = 5.14 * 1.2 * inp.su * 1e3; Hres = 4 * A * inp.su * 1e3; } // square footing shape factor 1.2 on Nc = 5.14
    else { const Nq = Math.exp(Math.PI * Math.tan(phi)) * Math.tan(Math.PI / 4 + phi / 2) ** 2, Ng = 1.5 * (Nq - 1) * Math.tan(phi); qult = 0.5 * inp.gammaSub * 1e3 * Bm * Ng * 0.6; Hres = W * Math.tan(phi); } // surface footing: ½γ′BNγ·sγ, sγ = 0.6 (square)
    const fsBear = qult * A / Math.max(Vmax, 1e-9), fsSl = Hres / Math.max(Fh, 1e-9), fsOt = (W * inp.Bj / 2) / Math.max(Mo, 1e-9);
    // --- Free-standing pile: fixed-base cantilever with hammer at the top, self-weight along the length.
    const D = inp.Dp / 1e3, t = inp.tp / 1e3;
    need(t < D / 2, 'Pile wall thickness must be smaller than the radius.');
    const sec = tube(D, t), E = steel.E || 2.1e11, fy = steel.fy, L = inp.Ls, n = inp.nElem;
    const P = inp.mHammer * 1e3 * G * inp.fLoad, q = sec.mass * G * inp.fLoad; // factored tip load [N] and self-weight [N/m]
    need(P + q * L > 0, 'Enter a hammer mass or a stick-up length above zero.');
    const nodes = []; for (let i = 0; i <= n; i++) nodes.push([0, 0, L * i / n]);
    const members = []; for (let i = 0; i < n; i++) members.push({ a: i, b: i + 1, section: 'p', material: 'st' });
    const model = { nodes, members, supports: [{ node: 0, fix: 'fixed' }], sections: { p: { A: sec.A, Iy: sec.I, Iz: sec.I, J: 2 * sec.I } }, materials: { st: { E, nu: 0.3, rho: 7850 } },
      loadCases: [{ name: 'axial', nodal: [{ node: n, F: [0, 0, -P, 0, 0, 0] }], member: members.map((m, i) => ({ member: i, type: 'uniform', q: [0, 0, -q], axes: 'global' })) }] };
    const bk = run(() => bucklingAnalysis(model, { loadCase: 'axial', nModes: 1, shearDeformation: false }), 'Buckling analysis');
    const lam = bk.factors[0];
    need(Number.isFinite(lam) && lam > 0, 'The buckling analysis did not return a positive load factor: check the pile data.');
    // Equivalent slenderness from the elastic critical load of the whole column: λ̄ = √(N_pl / N_cr) with N at the base.
    const Nbase = P + q * L, Ncr = lam * Nbase, Npl = sec.A * fy, lamBar = Math.sqrt(Npl / Ncr), chi = chiEC3(lamBar, 0.49), NbRd = chi * Npl / inp.gM;
    const beta = inp.batterDeg * DEG;
    // First-order base moment: hammer and self-weight components normal to the axis plus lateral load, amplified by 1/(1 − 1/λ).
    const M1 = P * Math.sin(beta) * L + q * Math.sin(beta) * L * L / 2 + inp.qLat * 1e3 * inp.fLoad * L * L / 2;
    const ampl = lam > 1 ? 1 / (1 - 1 / lam) : 99, M2 = M1 * ampl, MRd = fy * sec.Wp / inp.gM;
    const ucInt = Nbase * Math.cos(beta) / NbRd + M2 / MRd;
    const tipDefl = (P * Math.sin(beta) * L ** 3 / 3 + (q * Math.sin(beta) + inp.qLat * 1e3 * inp.fLoad) * L ** 4 / 8) / (E * sec.I) * ampl;
    const dOverT = D / t;
    const checks = [
      fsCheck('Unpiled jacket: mudmat bearing', fsBear, inp.fsBear, 'API RP 2A-WSD (2014) §6.13 (shallow foundations: bearing) with §6.4 safety factors', 'The leeward mudmats carry the weight plus the overturning couple. Bearing failure tilts the jacket before the piles can be stabbed.'),
      fsCheck('Unpiled jacket: sliding', fsSl, inp.fsSlide, 'API RP 2A-WSD (2014) §6.13.3 (sliding stability)', 'Wave, wind and current push the unpiled jacket sideways on its mudmats; sliding moves it off location and off level.'),
      fsCheck('Unpiled jacket: overturning', fsOt, inp.fsSlide, 'ISO 19902 (2020) §22 (on-bottom stability); rigid-body moment equilibrium about the leeward mudmat line', 'Stabilising moment of the on-bottom weight over the environmental overturning moment.'),
      fsCheck('Free-standing pile: elastic buckling factor', lam, 2.0, 'Timoshenko & Gere §2.12; ISO 19902 (2020) §22 (free-standing pile); factor of 2 on the elastic critical load as a screen', 'Load factor at which the pile with hammer buckles as a flagpole under factored weight. Close to 1 the deflections grow without limit.'),
      check('Free-standing pile: axial + bending interaction at the support', ucInt, 1, '-', 'EN 1993-1-1 (2005) §6.3.3 in simplified form: N/Nb,Rd + M·(1/(1 − 1/λcr))/M,Rd; API RP 2A-WSD §12.5.3', 'Compression from hammer and self-weight with second-order bending from batter and lateral load. Exceedance yields the pile at the sleeve or guide.'),
    ];
    if (Vmin < 0) checks.push(check('Unpiled jacket: uplift on the windward mudmats', -Vmin / 1e3, 1e-6, 'kN (tension not available)', 'Rigid-body equilibrium; API RP 2A-WSD §6.13', 'Mudmats cannot take tension: the windward side lifts and the jacket rocks.'));
    const gov = worst(checks), ok = gov.uc <= 1;
    const WobReq = Math.max(inp.fsSlide * Mo * 2 / inp.Bj, clay ? 0 : inp.fsSlide * Fh / Math.tan(phi));
    const LsMax = L * Math.sqrt(lam / 2);
    const hsS = linspace(0.5, Math.max(6, 1.5 * inp.Hs), 12);
    const fsVs = hsS.map((hs) => { const hd = Math.min(inp.fHmax * hs, 0.77 * inp.depth); const w2 = pileForceClosedForm({ H: hd, T: Tw, depth: inp.depth, D: inp.Deq, Cd: inp.Cd, Cm: inp.Cm, rho }); const f2 = inp.nLeg * w2.Fmax + inp.Fwind * 1e3, m2 = inp.nLeg * w2.Mmax + inp.Fwind * 1e3 * inp.zWind; return [Hres / f2, W * inp.Bj / 2 / m2, qult * A / (W / 4 + m2 / (2 * inp.Bj))]; });
    return {
      kpis: [
        kpi('Fh', 'Environmental base shear on the unpiled jacket', Fh / 1e3, 'kN', 'Maximum wave force on the equivalent cylinders plus the wind and current force.'),
        kpi('Mo', 'Overturning moment about the mudline', Mo / 1e3, 'kN·m', 'Wave moment plus wind and current force times its lever.'),
        kpi('fsBear', 'Mudmat bearing factor of safety', Math.min(fsBear, 99), '-', 'Ultimate bearing load of one mudmat over its largest reaction.', inp.fsBear, 'min'),
        kpi('fsSlide', 'Sliding factor of safety', Math.min(fsSl, 99), '-', clay ? 'Undrained base shear resistance su·A of the four mudmats over the base shear.' : 'Frictional resistance W′·tan φ over the base shear.', inp.fsSlide, 'min'),
        kpi('fsOt', 'Overturning factor of safety', Math.min(fsOt, 99), '-', 'Stabilising moment W′·B/2 over the overturning moment.', inp.fsSlide, 'min'),
        kpi('Vmax', 'Largest mudmat reaction', Vmax / 1e3, 'kN', 'Leeward mudmat load: W′/4 + M/(2B).'),
        kpi('Vmin', 'Smallest mudmat reaction', Vmin / 1e3, 'kN', 'Windward mudmat load; negative means the mudmat lifts off.', 0, 'min'),
        kpi('WobReq', 'On-bottom weight needed for sliding and overturning', WobReq / 1e3, 'kN', 'Weight at which the governing of the two just meets the required factor (in clay only overturning depends on weight).'),
        kpi('lambdaCr', 'Elastic buckling load factor of the free-standing pile', lam, '-', 'Multiplier on the factored hammer and self-weight at which the cantilever pile buckles (finite-element eigenvalue).', 2, 'min'),
        kpi('Ncr', 'Elastic critical load at the support', Ncr / 1e3, 'kN', 'Buckling factor times the axial force at the base.'),
        kpi('ucInt', 'Pile axial + bending utilisation', ucInt, '-', 'Beam–column interaction at the support with second-order amplified moment.', 1, 'max'),
        kpi('M2', 'Second-order bending moment at the support', M2 / 1e3, 'kN·m', 'First-order moment from batter and lateral load times the amplification 1/(1 − 1/λcr).'),
        kpi('tipDefl', 'Lateral deflection of the pile top', tipDefl, 'm', 'Elastic deflection at hammer level including second-order amplification; it must fit the hammer sleeve and guide clearances.'),
        kpi('LsMax', 'Longest stick-up with a buckling factor of 2', LsMax, 'm', 'Scaled from the present stick-up with the critical load varying roughly as 1/L².'),
      ],
      checks,
      plots: [
        line('On-bottom factors of safety against sea state', 'Significant wave height', 'm', 'Factor of safety', '-', [{ name: 'Sliding', x: hsS.map((v) => round(v)), y: fsVs.map((v) => round(Math.min(v[0], 20))) }, { name: 'Overturning', x: hsS.map((v) => round(v)), y: fsVs.map((v) => round(Math.min(v[1], 20))) }, { name: 'Bearing', x: hsS.map((v) => round(v)), y: fsVs.map((v) => round(Math.min(v[2], 20))) }, { name: 'Required', x: [hsS[0], hsS[hsS.length - 1]], y: [inp.fsSlide, inp.fsSlide], dash: true }],
          'How fast the stability of the unpiled jacket erodes with wave height. Where a curve meets the required factor is the limiting sea state for leaving the jacket unpiled.'),
        bar('Mudmat reactions', ['Leeward (max)', 'Windward (min)', 'Static W′/4'], 'Reaction', 'kN', [{ name: 'Reaction', y: [Vmax, Vmin, W / 4].map((v) => round(v / 1e3)) }, { name: 'Ultimate bearing', y: [qult * A, qult * A, qult * A].map((v) => round(v / 1e3)) }],
          'Reaction on the most and least loaded mudmats against the ultimate bearing load. A negative windward bar means lift-off.'),
      ],
      tables: [tbl('Temporary-condition summary', [['k', 'Quantity'], ['v', 'Value'], ['u', 'Unit']], [
        ['Design wave height', Hd, 'm'], ['Wave force per cylinder (drag amplitude)', wv.FD / 1e3, 'kN'], ['Wave force per cylinder (inertia amplitude)', wv.FI / 1e3, 'kN'], ['Ultimate bearing pressure', qult / 1e3, 'kPa'], ['Sliding resistance', Hres / 1e3, 'kN'],
        ['Pile D/t', dOverT, '-'], ['Factored axial force at the support', Nbase / 1e3, 'kN'], ['Relative slenderness', lamBar, '-'], ['Buckling reduction χ', chi, '-'], ['First-order moment at the support', M1 / 1e3, 'kN·m'], ['Moment amplification', ampl, '-'],
      ].map(([k, v, u]) => ({ k, v: round(v), u })), 'Intermediate values of both temporary conditions for checking against the installation engineering report.')],
      notes: [`Unpiled jacket in Hs ${inp.Hs} m (design wave ${fmt(Hd)} m): FS bearing ${fmt(fsBear)}, sliding ${fmt(fsSl)}, overturning ${fmt(fsOt)}.`, `Pile stick-up ${L} m with ${inp.mHammer} t on top: buckling factor ${fmt(lam)}, interaction ${fmt(ucInt)}, D/t = ${fmt(dOverT)}.`, `"${gov.name}" governs at ${fmt(gov.uc)}.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `Both temporary conditions are stable: the unpiled jacket has factors of ${fmt(Math.min(fsBear, fsSl, fsOt))} or more in Hs ${inp.Hs} m, and the free-standing pile has a buckling factor of ${fmt(lam)} with ${fmt(ucInt)} interaction.`,
          `State the limiting sea state for the unpiled condition in the procedure (see the factor-of-safety curves) and the weather window needed to stab and drive enough piles; limit the stick-up to ${fmt(LsMax)} m with this hammer.`, 'API RP 2A-WSD (2014) §6.13, §12.5; ISO 19902 (2020) §22')
          : adv('engineering', 'fail', `"${gov.name}" is not met (${fmt(gov.uc)}): FS bearing ${fmt(fsBear)}, sliding ${fmt(fsSl)}, overturning ${fmt(fsOt)}; pile buckling factor ${fmt(lam)}, interaction ${fmt(ucInt)}.`,
            gov.name.startsWith('Unpiled') ? (gov.name.includes('bearing') ? `Increase the mudmat area to about ${fmt(A * inp.fsBear / fsBear)} m² each, or restrict the unpiled condition to a lower sea state.` : `Add on-bottom weight (flood more legs) to at least ${fmt(WobReq / 1e3)} kN, widen the base, or stab piles before the weather exceeds the limit shown by the curves.`) : `Shorten the stick-up to ${fmt(LsMax)} m (add-on sections), use a lighter hammer, add a pile guide or gripper near the top, or increase the wall thickness.`,
            'API RP 2A-WSD (2014) §6.13, §12.5; ISO 19902 (2020) §22'),
        adv('economic', 'info', `The unpiled jacket depends on ${fmt(4 * A)} m² of mudmat and ${fmt(inp.Wob)} kN on-bottom weight; the pile is driven with ${L} m of stick-up.`,
          'Mudmats are installation-only steel: size them for the seasonal sea state of the actual installation window rather than an all-year storm, and trade pile add-on welds offshore (vessel hours each) against a longer stick-up with a guide frame.',
          'Installation-phase steel and offshore welding hours versus weather-window risk (DNV-ST-N001 Sec. 2); AACE RP 18R-97 for the estimate class',
          `Mudmat steel at about 0.25 t/m² is roughly ${steelCo2Text(4 * A * 0.25)} Each offshore pile add-on avoided saves several vessel hours: one day is about ${vesselCo2Text(1, 40)}`),
      ],
      publish: { 'install.onBottomFS': Math.min(fsBear, fsSl, fsOt, 99), 'install.stickupLambda': lam },
      quantities: { steel_t: round(4 * A * 0.25) },
    };
  },
  verify: [
    { name: 'Euler cantilever with tip load: Pcr = π²·E·I/(4·L²)', source: 'Euler (Timoshenko & Gere §2.1): self-weight removed by a negligible wall mass is not possible, so the tip load is made dominant (hammer 20 000 t on a 20 m stick-up; self-weight is 0.3 % of it and lowers λ by about 0.1 %): Pcr = π²·210e9·0.20973/(4·400) = 271.7 MN → λ = Pcr/(1.3·20 000 t·g)',
      inputs: { mHammer: 20000, Ls: 20, batterDeg: 0, qLat: 0 }, expect: [{ kpi: 'lambdaCr', value: Math.PI ** 2 * 210e9 * (Math.PI / 64 * (2.134 ** 4 - 2.014 ** 4)) / (4 * 400) / (1.3 * 2e7 * 9.80665), tol: 0.003 }] },
    { name: 'Column buckling under its own weight: (q·L)cr = 7.837·E·I/L²', source: 'Greenhill (1881); Timoshenko & Gere, Theory of Elastic Stability §2.12: no hammer, L = 150 m → λ = 7.837·E·I/(q·L³) with q = 1.3·ρ·g·A',
      inputs: { mHammer: 0, Ls: 150, batterDeg: 0, qLat: 0, nElem: 32 }, expect: [{ kpi: 'lambdaCr', value: 7.837 * 210e9 * (Math.PI / 64 * (2.134 ** 4 - 2.014 ** 4)) / (1.3 * 7850 * 9.80665 * Math.PI / 4 * (2.134 ** 2 - 2.014 ** 2) * 150 ** 3), tol: 0.005 }] },
    { name: 'Rigid-body on-bottom factors with no waves: FS sliding = 4·su·A/F, overturning = W′·B/(2·M)', source: 'Statics: 250 kN at 80 m lever only → sliding 4·25·120/250 = 48, overturning 12 000·45/2/(250·80) = 13.5', inputs: { Hs: 0, Wob: 12000 }, expect: [{ kpi: 'fsSlide', value: 48, tol: 1e-9 }, { kpi: 'fsOt', value: 13.5, tol: 1e-9 }] },
  ],
};
