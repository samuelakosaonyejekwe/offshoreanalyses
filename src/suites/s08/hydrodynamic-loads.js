// Hydrodynamic loading under a regular design wave and current: near-bed kinematics by Airy or Stokes second-order
// theory, Morison load history on the pipe, Keulegan–Carpenter number, seabed sediment mobility (Shields) and a
// sliding-block estimate of the lateral displacement through a wave train.
import { airyKinematics, stokes2Kinematics, waveProperties, selectWaveTheory } from '../../engine/metocean/waves.js';
import { grainRoughness, bedShearCurrent, bedShearWave, combinedShearSoulsby, mobilityCheck } from '../../engine/coastal/sediment.js';
import { addedMassCoefficient } from '../../engine/pipeline/span.js';
import { pipeInputs, pipeFrom, depthInput, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, linspace, sev, thin,
  ECON_BASIS, CO2_CONCRETE, MGO_T_PER_DAY, CO2_MGO } from './_common.js';

const G_ENV = 'Design wave and current';
const G_HYD = 'Hydrodynamic coefficients';
const G_SOIL = 'Seabed';

export default {
  id: 's08-hydrodynamic-loads',
  suite: 8,
  title: 'Hydrodynamic loads and stability under waves and current',
  summary: 'Builds the wave-and-current load history on the pipe, tests seabed sediment mobility and estimates how far the pipe slides in a wave train.',
  covers: [9],
  fidelity: 'reduced-order',
  standards: ['DNV-RP-C205 (2019) §3 (wave kinematics), §6 (Morison loads)', 'DNV-RP-F109 (2021) §3.5 (dynamic lateral stability, principle)', 'Soulsby, Dynamics of Marine Sands (1997) ch. 4–5'],
  models: ['Airy wave theory', 'Stokes second-order wave theory', 'Morison equation (drag, inertia, lift) with relative velocity', 'Keulegan–Carpenter number', 'Shields sediment mobility', 'Coulomb friction sliding block'],
  limits: 'Regular-wave, single-degree-of-freedom sliding block with constant Morison coefficients and Coulomb friction plus a constant passive resistance: it shows the order of displacement, not the calibrated generalised method of DNV-RP-F109 (whose peak-load coefficient tables are not built in). No wake model, no penetration build-up during sliding, no pipe bending or axial restraint. 3-D CFD is outside the module.',
  publishes: ['pipeline.hydroForceY'],
  searchTerms: 'pipeline hydrodynamic forces Morison Keulegan-Carpenter seabed lateral displacement',
  inputs: [
    ...pipeInputs({ omit: ['steel'] }),
    depthInput(50, G_ENV),
    num('H', 'Design wave height', 'm', 10, 0.1, 35, G_ENV, 'bc', 'Height of the regular design wave (e.g. the 100-year maximum wave).', { bind: 'site.Hmax100', uncertain: { dist: 'gumbel', cov: 0.1 } }),
    num('T', 'Wave period', 's', 12, 3, 25, G_ENV, 'bc', 'Period associated with the design wave.'),
    sel('theory', 'Wave theory', 'airy', [['airy', 'Airy (linear)'], ['stokes2', 'Stokes second order']], G_ENV, 'control', 'The recommended theory for the H, T, d combination is reported in the notes.'),
    num('Uc', 'Steady current at pipe level', 'm/s', 0.5, 0, 3, G_ENV, 'bc', 'Current normal to the pipe at pipe height, in the wave direction.', { bind: 'site.UcBed100', uncertain: { dist: 'gumbel', cov: 0.15 } }),
    num('heading', 'Angle between wave direction and pipe axis', 'deg', 90, 10, 90, G_ENV, 'bc', '90° = waves crossing the pipe at right angles; only the normal velocity component loads the pipe.'),
    num('gap', 'Gap between pipe and seabed', 'mm', 0, 0, 3000, G_ENV, 'geometry', 'Zero for a pipe resting on the bed. With a gap the pipe is treated as a span section and friction is not mobilised.'),
    num('CD', 'Drag coefficient', '-', 0.7, 0.3, 2.5, G_HYD, 'load', 'Effective value for a pipe on the bed 0.7; free-stream rough cylinder 1.05–1.2.', { uncertain: { dist: 'lognormal', cov: 0.15 } }),
    num('CM', 'Inertia coefficient', '-', 3.29, 1, 5, G_HYD, 'load', 'On the bed 3.29; free stream 2.0.'),
    num('CL', 'Lift coefficient', '-', 0.9, 0, 2, G_HYD, 'load', 'On the bed 0.9; falls quickly with a gap.'),
    num('mu', 'Pipe–soil friction coefficient', '-', 0.6, 0.05, 1.5, G_SOIL, 'soil', '0.6 sand, 0.2 clay for concrete-coated pipe (DNV-RP-F109).', { uncertain: { dist: 'lognormal', cov: 0.15 } }),
    num('FR', 'Passive soil resistance', 'N/m', 0, 0, 20000, G_SOIL, 'soil', 'Constant passive resistance from embedment (take it from the on-bottom stability or pipe–soil module).'),
    num('d50', 'Median grain size of the seabed', 'mm', 0.25, 0.01, 100, G_SOIL, 'soil', 'For the Shields mobility check of the seabed around the pipe.', { bind: 'site.d50' }),
    num('nWaves', 'Number of waves simulated', '-', 10, 2, 200, 'Simulation', 'control', 'Waves in the regular train; displacement accumulates wave by wave once sliding starts.'),
    { key: 'nSteps', label: 'Time steps per wave period', type: 'int', default: 180, min: 30, max: 2880, group: 'Simulation', role: 'control', help: 'Time resolution of the load history and of the sliding integration.' },
    { key: 'x0', label: 'Initial lateral offset', unit: 'm', type: 'number', default: 0, min: -50, max: 50, group: 'Simulation', role: 'ic', help: 'Starting lateral position of the pipe; it starts at rest.' },
  ],
  mesh: { key: 'nSteps', label: 'Time steps per wave period', levels: [90, 180, 360, 720], qoi: 'FYmax' },
  compute(inp, ctx) {
    const p = pipeFrom({ ...inp, steel: 'api5l-x65' }, ctx);
    const D = p.Do, d = inp.depth, e = inp.gap / 1e3;
    need(d > 1, 'Water depth must be more than 1 m.');
    need(e + D < d, 'The pipe (with its gap) does not fit in the water depth.');
    const wp = run(() => waveProperties({ T: inp.T, d, H: inp.H }), 'Wave properties');
    need(inp.H < 0.78 * d, `A ${fmt(inp.H)} m wave breaks in ${fmt(d)} m of water (H/d > 0.78): reduce the wave height or use a breaking-wave assessment.`);
    const rec = selectWaveTheory({ H: inp.H, T: inp.T, d });
    const zc = -d + e + D / 2;
    const sinH = Math.sin(inp.heading * Math.PI / 180);
    const kin = (t) => {
      const k = inp.theory === 'stokes2' ? stokes2Kinematics({ H: inp.H, T: inp.T, d, z: zc, t, rho: p.rhoW }) : airyKinematics({ H: inp.H, T: inp.T, d, z: zc, t, rho: p.rhoW, k: wp.k });
      return { u: k.u * sinH + inp.Uc, a: k.ax * sinH, eta: k.eta };
    };
    const A = Math.PI * D * D / 4, rho = p.rhoW;
    const Ca = addedMassCoefficient(e / D);
    const mTot = p.massOp + Ca * rho * A;
    const onBed = e <= 1e-6;
    const ws = p.ws;
    // explicit time stepping of the sliding block; friction acts only when the pipe touches the bed
    const n = inp.nSteps, dt = inp.T / n, N = Math.round(inp.nWaves) * n;
    const ts = new Float64Array(N + 1), FY = new Float64Array(N + 1), FZ = new Float64Array(N + 1), X = new Float64Array(N + 1), ETA = new Float64Array(N + 1), U = new Float64Array(N + 1);
    let x = inp.x0, v = 0, FYmax = 0, FZmax = 0, Umax = 0, ucMin = Infinity, liftOff = false, xMax = 0;
    for (let i = 0; i <= N; i++) {
      const t = i * dt, k = kin(t);
      const ur = k.u - v;
      const fy = 0.5 * rho * D * inp.CD * ur * Math.abs(ur) + rho * A * inp.CM * k.a;
      const fz = 0.5 * rho * D * inp.CL * ur * ur;
      const R = onBed ? inp.mu * Math.max(0, ws - fz) + inp.FR : 0;
      if (fz > ws) liftOff = true;
      ts[i] = t; FY[i] = fy; FZ[i] = fz; X[i] = x; ETA[i] = k.eta; U[i] = k.u;
      FYmax = Math.max(FYmax, Math.abs(fy)); FZmax = Math.max(FZmax, fz); Umax = Math.max(Umax, Math.abs(k.u));
      if (onBed && Math.abs(fy) > 0) ucMin = Math.min(ucMin, R / Math.abs(fy));
      if (onBed) {
        if (v === 0 && Math.abs(fy) <= R) { /* stuck */ } else {
          const dir = v !== 0 ? Math.sign(v) : Math.sign(fy);
          const acc = (fy - dir * R) / mTot;
          const vNew = v + acc * dt;
          if (v !== 0 && vNew * v < 0) { x += 0.5 * v * dt * Math.abs(v / (v - vNew)); v = 0; } // stops within the step
          else { x += 0.5 * (v + vNew) * dt; v = vNew; }
        }
      }
      xMax = Math.max(xMax, Math.abs(x - inp.x0));
    }
    const slide = Math.abs(x - inp.x0);
    const Uw = Umax - inp.Uc;
    const KC = Math.abs(Uw) * inp.T / D;
    const sf = Number.isFinite(ucMin) ? ucMin : 99;
    // sediment mobility of the surrounding seabed: wave + current skin friction (Soulsby)
    const d50 = inp.d50 / 1e3, z0 = grainRoughness(d50).z0;
    const bedU = Math.abs(airyKinematics({ H: inp.H, T: inp.T, d, z: -d, t: 0, k: wp.k }).u);
    const tw = bedU > 1e-6 ? bedShearWave({ Uw: bedU, T: inp.T, z0, rho }).tau : 0;
    const tc = inp.Uc > 1e-6 ? bedShearCurrent({ U: inp.Uc, h: d, z0, rho }).tau : 0;
    const tMax = tw > 0 && tc > 0 ? combinedShearSoulsby({ tauC: tc, tauW: tw, phiDeg: 0 }).tauMax : Math.max(tw, tc);
    const mob = run(() => mobilityCheck({ tau: Math.max(tMax, 1e-9), d: d50, rho }), 'Sediment mobility');
    const wsReq = onBed ? (FYmax + inp.mu * FZmax - inp.FR) / inp.mu : 0; // conservative: peaks taken as simultaneous
    const [tt, fyp, fzp, xp, up] = thin([ts, FY, FZ, X, U], 600);
    const res = fzp.map((fz) => round(inp.mu * Math.max(0, ws - fz) + inp.FR));
    const dWs = Math.max(0, wsReq - ws);
    const dConc = dWs / ((p.concRho - rho) * 9.80665); // m³ of concrete per metre that adds dWs of submerged weight
    return {
      kpis: [
        kpi('Uw', 'Wave-induced velocity amplitude at the pipe', Uw, 'm/s', 'Largest wave-induced horizontal velocity normal to the pipe at pipe level over the wave cycle (current excluded).'),
        kpi('KC', 'Keulegan–Carpenter number', KC, '-', 'U_w·T/D. Below about 5 the load is inertia dominated; above about 20 drag and lift dominate and wake effects grow.'),
        kpi('FYmax', 'Peak horizontal load', FYmax, 'N/m', 'Largest Morison drag-plus-inertia load on the pipe during the wave train.'),
        kpi('FZmax', 'Peak lift', FZmax, 'N/m', 'Largest vertical hydrodynamic load; it reduces the contact force and with it the friction.', ws, 'max'),
        kpi('sfSlide', 'Minimum resistance / load ratio', Math.min(sf, 99), '-', 'Lowest ratio of sliding resistance to horizontal load during the wave cycle. Below 1.0 the pipe slides at that instant (99 = not on the bed or no load).', 1, 'min'),
        kpi('slide', 'Net lateral displacement after the wave train', slide, 'm', 'Distance the pipe has moved from its starting position after the simulated waves. Zero means the pipe is absolutely stable in this wave.'),
        kpi('slideMax', 'Largest lateral excursion', xMax, 'm', 'Maximum distance from the starting position reached at any time; compare with the allowable displacement of the corridor (often 10 pipe diameters).', 10 * D, 'max'),
        kpi('wsReq', 'Submerged weight for no sliding', Math.max(wsReq, 0), 'N/m', 'Weight at which peak drag plus inertia is resisted with peak lift acting, taking the two peaks as simultaneous (conservative).'),
        kpi('shields', 'Shields parameter of the seabed', mob.theta, '-', 'Dimensionless bed shear stress under the design wave and current. Above the critical value (see the mobility ratio) the seabed is mobile, so scour, self-burial or loss of support can develop.'),
        kpi('mobility', 'Seabed mobility ratio θ / θcr', mob.theta / mob.thetaCr, '-', 'Shields parameter over its critical value (Soulsby 1997 Eq. 77). Above 1 the bed sediment moves in the design condition. That is usual for sand in a design storm and is a screening result for scour and embedment, not a pass/fail criterion.'),
        kpi('wavelength', 'Wavelength', wp.L, 'm', 'Length of the design wave at this depth from the linear dispersion relation.'),
        kpi('ursell', 'Ursell number', wp.ursell, '-', 'H·L²/d³. Above about 26 linear theory loses accuracy and a higher-order or stream-function theory is appropriate.'),
      ],
      checks: [
        chk('Sliding in the design wave', onBed ? 1 : 0, onBed ? Math.max(sf, 1e-6) : 1, '-', 'DNV-RP-F109 (2021) §3.6 (static criterion, unfactored)',
          'Guards against the horizontal wave-and-current load exceeding friction plus passive resistance at any instant of the wave; if exceeded the pipe starts to slide.'),
        chk('Lift against submerged weight', FZmax, ws, 'N/m', 'DNV-RP-F109 (2021) Eq. 3.39 (unfactored)',
          'Guards against the lift force unloading the seabed contact completely, after which the pipe has no lateral resistance.'),
      ],
      plots: [
        plot('line', 'Load history on the pipe', 'Time', 's', 'Load per unit length', 'N/m', [
          { name: 'Horizontal load', x: tt.map((v) => round(v)), y: fyp.map((v) => round(v)) }, { name: 'Lift', x: tt.map((v) => round(v)), y: fzp.map((v) => round(v)) },
          { name: 'Sliding resistance', x: tt.map((v) => round(v)), y: res, dash: true },
        ], 'Morison horizontal load and lift through the wave train with the instantaneous sliding resistance. Sliding occurs where the horizontal load rises above the dashed line.'),
        plot('line', 'Lateral displacement', 'Time', 's', 'Displacement', 'm', [{ name: 'Pipe position', x: tt.map((v) => round(v)), y: xp.map((v) => round(v, 5)) }],
          'Position of the pipe during the wave train. A staircase shows displacement accumulating each wave; a flat line shows a stable pipe.'),
        plot('line', 'Flow velocity at pipe level', 'Time', 's', 'Velocity', 'm/s', [{ name: 'Wave + current', x: tt.map((v) => round(v)), y: up.map((v) => round(v)) }],
          'Horizontal velocity normal to the pipe. Stokes theory sharpens the crest velocity and flattens the trough; the current shifts the whole curve.'),
      ],
      tables: [tbl('Flow and load summary', [col('item', 'Quantity'), col('value', 'Value'), col('unit', 'Unit')], [
        { item: 'Recommended wave theory for H, T, d', value: rec.theory, unit: '' },
        { item: 'Near-bed orbital velocity (linear)', value: round(bedU), unit: 'm/s' },
        { item: 'Wave skin-friction shear stress', value: round(tw), unit: 'Pa' }, { item: 'Current skin-friction shear stress', value: round(tc), unit: 'Pa' },
        { item: 'Maximum combined bed shear stress', value: round(tMax), unit: 'Pa' }, { item: 'Critical Shields parameter', value: round(mob.thetaCr), unit: '-' },
        { item: 'Added mass coefficient', value: round(Ca), unit: '-' }, { item: 'Mass incl. added mass', value: round(mTot), unit: 'kg/m' },
      ], 'Intermediate quantities behind the KPIs: the recommended wave theory, the bed shear stresses used in the Shields check and the inertia used in the sliding integration.')],
      notes: [
        `Wave theory used: ${inp.theory === 'stokes2' ? 'Stokes second order' : 'Airy'}; the regime chart suggests ${rec.theory} for this wave.`,
        liftOff ? 'Lift exceeds the submerged weight during part of the cycle: the pipe leaves the seabed.' : 'The pipe stays in contact with the seabed throughout the cycle.',
        mob.mobile ? 'The seabed is mobile in the design condition: expect scour and embedment changes.' : 'The seabed sediment is below the threshold of motion in the design condition.',
      ],
      advice: [
        adv('engineering', sf < 1 ? (xMax > 10 * D ? 'fail' : 'warn') : 'ok',
          sf < 1 ? `The pipe slides: the resistance falls to ${fmt(sf)} × the load and the pipe moves ${fmt(slide)} m in ${fmt(inp.nWaves)} waves (largest excursion ${fmt(xMax)} m).` : `The pipe is stable in the design wave: resistance is at least ${fmt(sf)} × the horizontal load.`,
          sf < 1 ? `Raise the submerged weight to about ${fmt(wsReq)} N/m, trench the section, or accept the displacement after a generalised-method assessment of the full storm.` : 'Confirm with the code-format check in the on-bottom stability module using the project γ_SC and with the irregular sea state.',
          'DNV-RP-F109 (2021) §3.5–3.6; DNV-RP-C205 (2019) §6 Morison loads'),
        ...(mob.mobile ? [adv('engineering', 'warn', `The seabed is mobile in the design condition: the Shields parameter is ${fmt(mob.theta / mob.thetaCr)} times its critical value.`,
          'Do not rely on the as-laid embedment: assess scour, free-span development and self-burial along the route, and use the embedment expected after sediment transport in the stability check.',
          'DNV-RP-F109 (2021) §3.4.4 (soil resistance and embedment); Soulsby (1997) Eq. 77; DNV-RP-F105 (2021) for the resulting free spans')] : []),
        adv('economic', 'info', dWs > 0 ? `Absolute stability needs ${fmt(dWs)} N/m more submerged weight, about ${fmt(dConc * 1e3)} m³ of extra weight coating per kilometre.` : `The pipe carries ${fmt(ws - Math.max(wsReq, 0))} N/m more submerged weight than this wave requires.`,
          dWs > 0 ? 'Compare the extra coating with one trenching pass; allowing a limited displacement usually removes most of the extra weight.' : 'A lighter coating may be possible for this section; check the installation and empty cases before reducing it.',
          ECON_BASIS,
          `Each cubic metre of weight coating is about ${CO2_CONCRETE} t CO₂e (screening factor); a trenching spread burns about ${MGO_T_PER_DAY} t of fuel per day ≈ ${fmt(MGO_T_PER_DAY * CO2_MGO)} t CO₂ per day. A mobile seabed that buries the pipe naturally avoids both.`),
      ],
      publish: { 'pipeline.hydroForceY': FYmax },
    };
  },
  verify: [{
    // Airy near-bed velocity amplitude: U = πH/(T sinh kd). d = 30 m, T = 10 s: ω² = g k tanh(kd) gives
    // k = 0.04578 rad/m (kd = 1.3734, sinh kd = 1.8478). H = 5 m: U = π·5/(10·1.8478) = 0.850 m/s at the bed.
    // The pipe axis is 0.2 m above the bed: cosh(k·0.2) = 1.00004, so the same value to four figures.
    name: 'Airy near-bed orbital velocity, U = πH/(T sinh kd)',
    source: 'Dean & Dalrymple, Water Wave Mechanics for Engineers and Scientists (1991) §4.4; hand calculation',
    inputs: { D: 400, t: 20, tCoat: 0, tConc: 0, depth: 30, H: 5, T: 10, Uc: 0, theory: 'airy', heading: 90 },
    expect: [{ kpi: 'Uw', value: 0.850, tol: 0.005 }, { kpi: 'KC', value: 21.25, tol: 0.006 }],
  }],
};
