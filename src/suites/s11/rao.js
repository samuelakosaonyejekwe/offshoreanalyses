// Response amplitude operators: six-degree-of-freedom motion RAOs with phases, and displacement, velocity and
// acceleration RAOs at a point, from the built-in models or from an imported hydrodynamic database.
import { pointMotion } from '../../engine/hydro/floater.js';
import { parseHydroJson, parseHydroCsv, raoSolve, checkHydroDb } from '../../engine/hydro/hydrodb.js';
import { hullInputs, hullFrom, motionInputs, motionFrom, motionRangeAdvice, seaInputs, seaFrom, omegaGrid, pointInputs, massMatrix, DOF, DOF_UNIT, dofScale,
  G, DEG, LIMIT_HYDRO, LIMIT_CFD, need, run, round, rnd, fmt, adv, check, steelCo2Text } from './_common.js';

const at = (xs, ys, x) => { if (x <= xs[0]) return ys[0]; for (let i = 1; i < xs.length; i++) if (x <= xs[i]) return ys[i - 1] + (ys[i] - ys[i - 1]) * (x - xs[i - 1]) / (xs[i] - xs[i - 1]); return ys[ys.length - 1]; };
const peak = (om, a) => { let j = 0; for (let i = 1; i < a.length; i++) if (a[i] > a[j]) j = i; return { value: a[j], period: 2 * Math.PI / om[j], edge: j === 0 }; };

function parseDb(text, rho, depth) {
  const t = String(text).trim();
  const db = run(() => (t.startsWith('{') ? parseHydroJson(JSON.parse(t)) : parseHydroCsv(t, { rho, g: G, depth })), 'Imported hydrodynamic database');
  need(db.omega.length >= 2, 'The imported hydrodynamic database needs at least two frequencies.');
  need(db.headingsDeg.length >= 1, 'The imported hydrodynamic database contains no wave excitation (F rows): RAOs cannot be solved.');
  return db;
}

export default {
  id: 's11-rao',
  suite: 11,
  title: 'Response amplitude operators',
  summary: 'Motion RAOs and phases for all six degrees of freedom and displacement, velocity and acceleration RAOs at a point, from the slender-member or box model or from imported radiation–diffraction coefficients.',
  covers: [4],
  fidelity: 'reduced-order',
  standards: ['DNV-RP-C205 (2019) §7.1–7.3', 'ISO 19904-1 (2019) §9', 'ITTC Recommended Procedure 7.5-02-07-03.1 (floating offshore platform experiments)', 'WAMIT v7 User Manual (numeric output conventions)'],
  models: ['Linear frequency-domain equation of motion [−ω²(M + A) + iωB + C + K]ξ = F', 'Imported added mass, radiation damping and excitation (WAMIT-style JSON or CSV)', 'Long-wave Froude–Krylov and Morison excitation on slender members', 'Box-barge heave, pitch and roll from a three-dimensional panel solution with the free-surface Green function (frequency-dependent added mass, wave damping and Haskind excitation)', 'Haskind / symmetry consistency checks of imported coefficients', 'Rigid-body point kinematics (displacement, velocity, acceleration)'],
  limits: `${LIMIT_HYDRO} Imported coefficients are used as given, interpolated linearly in frequency and heading, about the database origin, which must be the waterline point above the keel origin; viscous damping is then the user fraction of critical added on the diagonal. RAOs are per unit wave amplitude without drag unless a sea state is supplied for the linearisation. ${LIMIT_CFD}`,
  searchTerms: 'response amplitude operator RAO floating platform WAMIT import',
  inputs: [
    ...hullInputs({ type: 'semi', withMesh: false }),
    { key: 'source', label: 'Hydrodynamic model', type: 'select', default: 'model', group: 'Hydrodynamic model', role: 'bc',
      options: [{ value: 'model', label: 'Built-in (slender members / box barge)' }, { value: 'imported', label: 'Imported coefficients' }], help: 'Use imported coefficients whenever a radiation–diffraction analysis exists; the built-in models are for concept work.' },
    { key: 'hydroDb', label: 'Hydrodynamic database (JSON or CSV)', type: 'text', default: '', group: 'Hydrodynamic model', role: 'bc', showIf: { key: 'source', in: ['imported'] },
      help: 'JSON: { rho, g, depth, omega[], headingsDeg[], A[nw][6][6], B[nw][6][6], Fre[nb][nw][6], Fim[nb][nw][6] } in SI. CSV with header kind,omega,headingDeg,i,j,re,im (kind A | B | F). Convention exp(+iωt), unit wave amplitude, heading = propagation direction.' },
    ...seaInputs({ Hs: 6, Tp: 11, ret: '1', heading: 0, withDuration: false }),
    { key: 'dragSea', label: 'Linearise drag for this sea state', type: 'bool', default: true, group: 'Sea state', role: 'bc', help: 'Viscous drag depends on amplitude: with this on, the RAOs include drag damping linearised for the sea state above; off gives undamped-by-drag RAOs.' },
    ...motionInputs(),
    ...pointInputs(),
  ],
  mesh: { key: 'nFreq', label: 'Frequency points', levels: [50, 100, 200], qoi: 'heavePeak' },
  compute(inp) {
    const hull = hullFrom(inp); const mot = motionFrom(inp, hull);
    const om = omegaGrid(inp.nFreq); const sea = seaFrom(inp);
    let r; let dbInfo = null; let src;
    if (inp.source === 'imported') {
      need(String(inp.hydroDb || '').trim().length > 0, 'Paste the hydrodynamic database (JSON or CSV) or switch the hydrodynamic model to "Built-in".');
      const db = parseDb(inp.hydroDb, hull.rho, inp.depth);
      const M = massMatrix(hull.mass, [hull.cog[0], hull.cog[1], hull.KG - hull.T], mot.gyr);
      const C = hull.hs.stiffness.map((row) => row.slice());
      // user damping fraction on the restored degrees of freedom, in addition to the imported radiation damping
      const Bext = Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]);
      for (let d = 0; d < 6; d++) { const k = C[d][d] + mot.Kext[d][d]; if (k > 0) Bext[d][d] = 2 * inp.zetaRad * Math.sqrt(k * (M[d][d] + db.A[0][d][d])); }
      const w = om.filter((x) => x >= db.omega[0] && x <= db.omega[db.omega.length - 1]);
      need(w.length >= 2, `The database covers ${fmt(db.omega[0])}–${fmt(db.omega[db.omega.length - 1])} rad/s, which does not overlap the analysis grid.`);
      r = run(() => raoSolve({ db, M, C, Kext: mot.Kext, Bext, omega: w, headingDeg: inp.headingDeg }), 'RAO solution from imported coefficients');
      r.depth = inp.depth; r.g = G;
      const ck = checkHydroDb(db, {});
      dbInfo = { n: db.omega.length, nb: db.headingsDeg.length, ok: ck.ok, neg: ck.negativeDamping, asym: Math.max(ck.asymmetryA ?? 0, ck.asymmetryB ?? 0) };
      src = `imported coefficients (${db.omega.length} frequencies, ${db.headingsDeg.length} heading(s))`;
    } else {
      r = mot.rao(om, inp.headingDeg, inp.dragSea ? sea : null);
      src = mot.kind === 'slender' ? 'slender-member model' : 'box-barge model';
    }
    const w = Array.from(r.omega); const T = w.map((x) => round(2 * Math.PI / x, 4));
    for (const a of r.amp) need(a.every(Number.isFinite), 'The RAO solution is singular at one or more frequencies: a degree of freedom has neither restoring nor inertia. Check the mass data and the station-keeping stiffness.');
    const pm = pointMotion({ rao: r, point: [inp.px, inp.py, inp.pz - hull.T] });
    const pk = r.amp.map((a) => peak(w, a));
    const wp = 2 * Math.PI / inp.Tp;
    const dOL = hull.type === 'semi' ? inp.Dcol : hull.type === 'spar' ? inp.Dspar : hull.type === 'tlp' ? inp.DcolT : hull.Bc;
    const Tvalid = mot.kind === 'slender' && inp.source !== 'imported' ? 2 * Math.PI / mot.validOmegaMax : 0;
    const chk = [check('Wave period range covered by the hydrodynamic model', Tvalid, Math.max(inp.Tp / 1.4, 1e-3), 's', 'DNV-RP-C205 (2019) §6.1.2 (slender-body limit D/λ < 0.2)',
      'Model-validity check: the shortest period at which the slender-member loads are valid should lie below the bulk of the wave energy (taken as Tp/1.4 ≈ Tz); otherwise the short-wave part of the RAO is over-predicted and diffraction coefficients should be imported.')];
    if (dbInfo) chk.push(check('Imported damping non-negative and matrices symmetric', dbInfo.ok ? 0 : 1.5, 1, '-', 'Newman (1977) Marine Hydrodynamics ch. 6 (properties of radiation coefficients)', 'Consistency check of the imported database: negative radiation damping or asymmetric matrices indicate irregular frequencies, unit errors or a wrong reference point.'));
    const thinIdx = w.map((_, i) => i).filter((i) => i % Math.ceil(w.length / 50) === 0);
    return {
      kpis: [
        { key: 'heavePeak', label: 'Peak heave RAO', value: round(pk[2].value, 4), unit: 'm/m', explain: 'Largest heave amplitude per metre of wave amplitude. Values above 1 mean resonance amplifies the wave; the height of the peak is set by the damping.' },
        { key: 'heavePeakT', label: 'Period of the heave peak', value: round(pk[2].period, 4), unit: 's', explain: 'Wave period at which heave is largest: normally the heave natural period. Compare it with the storm peak periods at the site.' },
        { key: 'rollPeak', label: 'Peak roll RAO', value: Math.abs(pk[3].value * dofScale(3)) < 1e-6 ? 0 : round(pk[3].value * dofScale(3), 4), unit: 'deg/m', explain: 'Largest roll angle per metre of wave amplitude for this heading (zero in head or following seas). A response below 1e-6 deg/m is round-off of a symmetric case and is reported as 0.' },
        { key: 'pitchPeak', label: 'Peak pitch RAO', value: round(pk[4].value * dofScale(4), 4), unit: 'deg/m', explain: 'Largest pitch angle per metre of wave amplitude for this heading.' },
        { key: 'pitchPeakT', label: 'Period of the pitch peak', value: round(pk[4].period, 4), unit: 's', explain: 'Wave period at which pitch is largest.' },
        { key: 'surgeAtTp', label: 'Surge RAO at Tp', value: round(at(w, r.amp[0], wp), 4), unit: 'm/m', explain: 'Surge amplitude per metre of wave amplitude at the sea-state peak period; in long waves it tends to the horizontal particle motion.' },
        { key: 'raoHeaveAtTp', label: 'Heave RAO at Tp', value: round(at(w, r.amp[2], wp), 4), unit: 'm/m', explain: 'Heave amplitude per metre of wave amplitude at the peak period: the number that scales most of the heave response in that sea.' },
        { key: 'raoHeaveLong', label: 'Heave RAO at the longest period', value: round(r.amp[2][0], 4), unit: 'm/m', explain: `Heave RAO at ${fmt(T[0])} s. It should approach 1.0 for a free-floating body in very long waves (the unit rides the surface) — a quick sanity check of any RAO set.` },
        { key: 'ptDispPeak', label: 'Peak vertical displacement RAO at the point', value: round(Math.max(...pm.amp[2]), 4), unit: 'm/m', explain: 'Largest vertical motion of the point per metre of wave amplitude, combining heave with the roll and pitch levers in phase.' },
        { key: 'ptVelPeak', label: 'Peak vertical velocity RAO at the point', value: round(Math.max(...pm.velAmp[2]), 4), unit: '(m/s)/m', explain: 'Largest vertical velocity of the point per metre of wave amplitude.' },
        { key: 'ptAccPeak', label: 'Peak vertical acceleration RAO at the point', value: round(Math.max(...pm.accAmp[2]), 4), unit: '(m/s²)/m', explain: 'Largest vertical acceleration of the point per metre of wave amplitude; short periods dominate because acceleration scales with ω².' },
        { key: 'Tvalid', label: 'Shortest valid wave period', value: round(Tvalid, 4), unit: 's', explain: 'Below this period the governing member exceeds D/λ = 0.2 and the slender-member excitation is no longer valid (zero for the box model and imported coefficients, which have no such limit).' },
      ],
      checks: chk,
      plots: [
        { type: 'line', title: 'Translation RAOs', x: { label: 'Wave period', unit: 's' }, y: { label: 'Amplitude', unit: 'm/m' }, series: [0, 1, 2].map((d) => ({ name: DOF[d], x: T, y: rnd(r.amp[d]) })),
          explain: 'Surge, sway and heave amplitude per metre of wave amplitude. Look for the heave resonance peak, the cancellation trough of a semi-submersible and the long-period asymptote of 1.0 in heave.' },
        { type: 'line', title: 'Rotation RAOs', x: { label: 'Wave period', unit: 's' }, y: { label: 'Amplitude', unit: 'deg/m' }, series: [3, 4, 5].map((d) => ({ name: DOF[d], x: T, y: rnd(r.amp[d].map((v) => v / DEG)) })),
          explain: 'Roll, pitch and yaw amplitude per metre of wave amplitude. The peaks mark the roll and pitch natural periods; in long waves pitch tends to the wave slope.' },
        { type: 'line', title: 'RAO phases (heave, roll, pitch)', x: { label: 'Wave period', unit: 's' }, y: { label: 'Phase lead', unit: 'deg' }, series: [2, 3, 4].map((d) => ({ name: DOF[d], x: T, y: rnd(r.phase[d].map((v) => v / DEG)) })),
          explain: 'Phase lead relative to the wave crest at the origin. A 180° swing through a narrow band marks a resonance; phases are needed to combine motions at a point (box-model heave and pitch carry no phase information and show zero).' },
        { type: 'line', title: 'Point RAOs: vertical displacement, velocity, acceleration', x: { label: 'Wave period', unit: 's' }, y: { label: 'Per metre wave amplitude', unit: 'm, m/s, m/s²' },
          series: [{ name: 'Displacement (m/m)', x: T, y: rnd(pm.amp[2]) }, { name: 'Velocity ((m/s)/m)', x: T, y: rnd(pm.velAmp[2]) }, { name: 'Acceleration ((m/s²)/m)', x: T, y: rnd(pm.accAmp[2]) }],
          explain: 'Vertical motion of the chosen point. Velocity and acceleration weight the short periods (×ω and ×ω²), so their peaks sit at shorter periods than the displacement peak.' },
      ],
      tables: [
        { title: 'RAO table', columns: [{ key: 'T', label: 'Period', unit: 's' }, ...DOF.map((n, d) => ({ key: `a${d}`, label: n, unit: d < 3 ? 'm/m' : 'deg/m' })), { key: 'pz', label: 'Point vertical', unit: 'm/m' }, { key: 'pa', label: 'Point vert. acc.', unit: '(m/s²)/m' }],
          rows: thinIdx.map((i) => ({ T: T[i], ...Object.fromEntries(DOF.map((_, d) => [`a${d}`, round(r.amp[d][i] * dofScale(d), 4)])), pz: round(pm.amp[2][i], 4), pa: round(pm.accAmp[2][i], 4) })),
          explain: 'Amplitude operators by wave period for this heading, in the layout of a motion report; each row is one frequency.' },
        { title: 'RAO peaks', columns: [{ key: 'dof', label: 'DOF' }, { key: 'pk', label: 'Peak RAO' }, { key: 'T', label: 'At period', unit: 's' }, { key: 'u', label: 'Unit' }],
          rows: DOF.map((n, d) => ({ dof: n, pk: round(pk[d].value * dofScale(d), 4), T: round(pk[d].period, 4), u: d < 3 ? 'm/m' : 'deg/m' })),
          explain: 'Largest RAO of each degree of freedom and where it occurs. A peak at the longest period of the grid means the response is still rising towards a resonance outside the range.' },
      ],
      notes: [
        `RAOs from the ${src} for heading ${fmt(inp.headingDeg)}°${inp.source !== 'imported' && inp.dragSea && mot.kind === 'slender' ? `, drag linearised for Hs ${fmt(inp.Hs)} m / Tp ${fmt(inp.Tp)} s (${r.converged ? 'converged' : 'not converged'})` : ''}.`,
        `Heave peaks at ${fmt(pk[2].period)} s with ${fmt(pk[2].value)} m/m; pitch peaks at ${fmt(pk[4].period)} s.`,
        dbInfo ? `Database check: ${dbInfo.ok ? 'damping non-negative and matrices symmetric' : `problems found (negative damping entries: ${dbInfo.neg ? 'yes' : 'no'}, asymmetry ${fmt(dbInfo.asym)})`}.` : `Governing member dimension ${fmt(dOL)} m; the added damping of ${fmt(100 * inp.zetaRad)} % of critical is an assumption to calibrate.`,
      ],
      advice: [
        ...(inp.source === 'imported' ? [] : motionRangeAdvice(mot, inp.Tp, inp.headingDeg)),
        adv('engineering', chk.some((c) => c.uc > 1) ? 'warn' : 'ok', `Heave RAO peaks at ${fmt(pk[2].value)} m/m at ${fmt(pk[2].period)} s; at the sea-state peak period it is ${fmt(at(w, r.amp[2], wp))} m/m. The long-period heave asymptote is ${fmt(r.amp[2][0])}.`,
          inp.source === 'imported' ? 'Keep the reference point, units and phase convention of the database consistent with this suite (SI, exp(+iωt), origin at the waterline above the keel origin) and add viscous damping from a decay test.'
            : `Treat these as concept-level RAOs: resonant peaks scale with 1/damping (${fmt(100 * inp.zetaRad)} % of critical added to the computed wave radiation and drag). Run a radiation–diffraction analysis, import its coefficients here and compare before design decisions.`,
          'DNV-RP-C205 (2019) §7.1–7.3; ITTC 7.5-02-07-03.1 for model-test calibration of damping.'),
        adv('economic', 'info', `Moving the heave natural period (now ${fmt(pk[2].period)} s) by a few seconds changes the response in every sea state of the scatter diagram, and with it operability and fatigue.`,
          'Settle natural periods at concept stage, when it costs only draft and pontoon proportions; correcting motions later means added steel or heave plates on a built hull.',
          'Front-end loading practice (AACE RP 18R-97 class 5–4 estimates); value-of-information for commissioning a diffraction analysis or model test.',
          `A heave-plate or pontoon modification of 500 t to shift a period after hull fabrication: ${steelCo2Text(500)}`),
      ],
    };
  },
  verify: [
    { // A 30 × 10 × 1 m barge at the longest period of the grid (2π/0.08 = 78.5 s, λ ≈ 9.6 km): heave RAO → 1.
      name: 'Long-wave limit: heave RAO of a free-floating body tends to 1',
      source: 'Newman (1977) Marine Hydrodynamics ch. 6 — long-wavelength limit of heave response',
      inputs: { hullType: 'barge', L: 30, B: 10, Dhull: 4, draft: 1, KG: 1.5, headingDeg: 180 },
      expect: [{ kpi: 'raoHeaveLong', value: 1.0, tol: 0.01 }] },
    { // Imported database with A = B = 0 and heave excitation equal to the hydrostatic stiffness ρ g Aw (a body that does not
      // disturb the wave): heave RAO = 1/(1 − ω²·T/g) for a box of draft T. T = 6 m, Tp = 15.708 s (ω = 0.4):
      // 1/(1 − 0.16·6/9.80665) = 1.1085.
      name: 'Imported coefficients: undamped heave oscillator RAO = 1/(1 − ω²T/g)',
      source: 'Single-degree-of-freedom frequency response, hand calculation',
      inputs: { hullType: 'barge', L: 100, B: 20, Dhull: 12, draft: 6, KG: 5, source: 'imported', Tp: 15.70796, headingDeg: 0, zetaRad: 0.005, Kmoor: 250, nFreq: 200,
        hydroDb: JSON.stringify((() => { const omega = [0.05, 0.2, 0.4, 0.8, 1.6, 2.6]; const z = () => Array.from({ length: 6 }, () => [0, 0, 0, 0, 0, 0]); const c33 = 1025 * 9.80665 * 2000;
          return { rho: 1025, g: 9.80665, depth: null, omega, headingsDeg: [0], A: omega.map(z), B: omega.map(z), Fre: [omega.map(() => [0, 0, c33, 0, 0, 0])], Fim: [omega.map(() => [0, 0, 0, 0, 0, 0])] }; })()) },
      expect: [{ kpi: 'raoHeaveAtTp', value: 1.1085, tol: 0.01 }] },
  ],
};
