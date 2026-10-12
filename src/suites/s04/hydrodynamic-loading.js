// Member-by-member hydrodynamic loading of a subsea structure in the near-bed wave and current field.
import { nearBedWaveCurrent, currentLogNearBed } from '../../engine/metocean/current.js';
import { dragCoefficientSteady, wakeAmplification } from '../../engine/hydro/morison.js';
import { G, DEG, NU_SEA, need, run, fmt, round, adv, steelCo2Text, waterInputs } from './_common.js';

const DEFAULT_MEMBERS = [
  { name: 'Base beams across the flow', D: 0.5, len: 14, z: 0.6, n: 2, Cd: 1.2, Cm: 2.0 },
  { name: 'Corner posts', D: 0.4, len: 3.4, z: 2.0, n: 4, Cd: 1.05, Cm: 2.0 },
  { name: 'Top frame across the flow', D: 0.4, len: 14, z: 3.8, n: 2, Cd: 1.05, Cm: 2.0 },
  { name: 'Manifold headers and valves (equivalent)', D: 0.8, len: 8, z: 2.2, n: 1, Cd: 1.2, Cm: 2.0 },
  { name: 'Roof panel edges', D: 0.15, len: 14, z: 4.0, n: 1, Cd: 2.0, Cm: 1.5 },
];

export default {
  id: 's04-hydrodynamic-loading', suite: 4, title: 'Hydrodynamic loading and seabed current',
  summary: 'Drag, inertia and lift on every member and equipment item from the near-bed current profile and wave orbital motion, with the resultant load and its lever arm.',
  covers: [8], fidelity: 'analytical',
  standards: ['DNV-RP-C205 (2019) §6.2, §6.7, App. D', 'DNV-RP-F109 (2021) §3.4 Eq. 3.2', 'ISO 13628-1 / API RP 17A (subsea structures, environmental loads)'],
  models: ['Morison wave-current loading', 'Near-bed logarithmic velocity profile', 'Linear wave theory (near-bed orbital velocity)', 'Keulegan–Carpenter number', 'Reynolds-number / roughness dependent drag', 'Added-mass inertia'],
  limits: 'Each item is a slender Morison body normal to the flow at one elevation; shielding and blockage are a single user factor, the wave is a regular linear wave and wave-current boundary-layer interaction is not modelled. Diffraction, proximity effects between members and CFD-derived loads are outside the module.',
  publishes: ['subsea.hydroFH', 'subsea.hydroFL', 'subsea.hydroArm'],
  searchTerms: 'Morison near-bed wave current load subsea template drag coefficient',
  inputs: [
    { key: 'members', label: 'Members and equipment exposed to the flow', type: 'table', group: 'Structure', role: 'geometry',
      columns: [{ key: 'name', label: 'Item', type: 'text' }, { key: 'D', label: 'Width normal to flow', unit: 'm', type: 'number' }, { key: 'len', label: 'Length', unit: 'm', type: 'number' }, { key: 'z', label: 'Centroid height above bed', unit: 'm', type: 'number' }, { key: 'n', label: 'Number', type: 'number' }, { key: 'Cd', label: 'Cd', type: 'number' }, { key: 'Cm', label: 'Cm', type: 'number' }],
      default: DEFAULT_MEMBERS, help: 'One row per group of identical items. Width × length is the projected area; the inertia volume is taken as π·D²/4 × length.' },
    { key: 'cdMode', label: 'Drag coefficient', type: 'select', default: 'table', group: 'Hydrodynamic coefficients', role: 'load',
      options: [{ value: 'table', label: 'As tabulated per item' }, { value: 'c205', label: 'DNV-RP-C205 roughness and KC (circular members, KC < 12)' }],
      help: 'The C205 option uses C_DS(k/D)·ψ(KC) where the closed expression exists (KC < 12) and the tabulated value elsewhere.' },
    { key: 'kRough', label: 'Surface roughness k', unit: 'mm', type: 'number', default: 5, min: 0.001, max: 100, group: 'Hydrodynamic coefficients', role: 'load', showIf: { key: 'cdMode', in: ['c205'] }, help: 'New painted steel ≈ 0.005 mm; marine growth 5–50 mm.' },
    { key: 'shield', label: 'Shielding / blockage factor', unit: '-', type: 'number', default: 1.0, min: 0.4, max: 1.2, group: 'Hydrodynamic coefficients', role: 'load', uncertain: { dist: 'normal', cov: 0.1 },
      help: 'Multiplier on all member loads. 1.0 ignores shielding; calibrate against flume tests of the actual arrangement.' },
    { key: 'Aplan', label: 'Solid plan area exposed to lift', unit: 'm²', type: 'number', default: 60, min: 0, max: 5000, group: 'Hydrodynamic coefficients', role: 'geometry' },
    { key: 'Cl', label: 'Lift coefficient on the plan area', unit: '-', type: 'number', default: 0.3, min: 0, max: 2, group: 'Hydrodynamic coefficients', role: 'load', uncertain: { dist: 'normal', cov: 0.25 } },
    { key: 'zLift', label: 'Height used for the lift velocity', unit: 'm', type: 'number', default: 3.8, min: 0.05, max: 40, group: 'Hydrodynamic coefficients', role: 'geometry' },
    ...waterInputs(),
    { key: 'z0mm', label: 'Seabed roughness z0', unit: 'mm', type: 'number', default: 0.04, min: 0.001, max: 50, group: 'Site', role: 'environment', help: 'DNV-RP-F109 Table 3-1: silt 0.005, fine sand 0.01, medium sand 0.04, coarse sand 0.1, gravel 0.3 mm.' },
    { key: 'Hs', label: 'Significant wave height', unit: 'm', type: 'number', default: 13, min: 0, max: 25, group: 'Wave and current', role: 'bc', bind: 'site.Hs100', uncertain: { dist: 'normal', cov: 0.08 } },
    { key: 'Tp', label: 'Peak period', unit: 's', type: 'number', default: 15.5, min: 3, max: 25, group: 'Wave and current', role: 'bc', bind: 'site.Tp100' },
    { key: 'hFac', label: 'Design wave height / Hs', unit: '-', type: 'number', default: 1.86, min: 1, max: 2, group: 'Wave and current', role: 'bc' },
    { key: 'Uref', label: 'Current speed at the reference height', unit: 'm/s', type: 'number', default: 0.7, min: 0, max: 3, group: 'Wave and current', role: 'bc', bind: 'site.UcBed100', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'zRef', label: 'Reference height of the current', unit: 'm', type: 'number', default: 1, min: 0.1, max: 20, group: 'Wave and current', role: 'bc' },
    { key: 'angDeg', label: 'Angle between current and waves', unit: 'deg', type: 'number', default: 0, min: 0, max: 90, group: 'Wave and current', role: 'bc', help: 'Loads are resolved along the wave direction and across it.' },
    { key: 'nPhase', label: 'Phase steps per wave cycle', type: 'int', default: 72, min: 8, max: 720, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nPhase', label: 'Phase steps', levels: [18, 36, 72, 144], qoi: 'FH' },
  compute(inp) {
    // SI: mm → m, deg → rad; forces reported in kN.
    const rows = (inp.members || []).filter((m) => m.D > 0 && m.len > 0 && m.n > 0);
    need(rows.length > 0, 'Give at least one member with positive width, length and number.');
    for (const m of rows) need(m.z > 0 && m.Cd >= 0 && m.Cm >= 0, `Member "${m.name}": height above bed must be positive and coefficients non-negative.`);
    need(inp.zLift < inp.depth && inp.zRef < inp.depth, 'The lift and current reference heights must be below the water surface.');
    const rho = inp.rhoW, z0 = inp.z0mm / 1e3, H = inp.hFac * inp.Hs, T = inp.Tp, om = 2 * Math.PI / T, ang = inp.angDeg * DEG;
    const uwAt = (z) => (H > 0 ? run(() => nearBedWaveCurrent({ H, T, d: inp.depth, zBed: z }), 'Wave kinematics').Uw : 0);
    const ucAt = (z) => (inp.Uref > 0 ? currentLogNearBed(z, { vRef: inp.Uref, zRef: inp.zRef, z0 }) : 0);
    for (const m of rows) need(m.z < inp.depth, `Member "${m.name}" lies above the water surface (height ${fmt(m.z)} m in ${fmt(inp.depth)} m of water).`);
    const zTop = Math.min(Math.max(...rows.map((m) => m.z)) * 1.3, inp.depth * 0.999);
    const items = rows.map((m) => {
      const Uw = uwAt(m.z), Uc = ucAt(m.z);
      const KC = Uw * T / m.D, Re = (Uc + Uw) * m.D / NU_SEA;
      let Cd = m.Cd, cdSrc = 'table';
      if (inp.cdMode === 'c205' && KC > 0 && KC < 12) { const cds = dragCoefficientSteady(inp.kRough / 1e3 / m.D); Cd = cds * wakeAmplification(KC, cds); cdSrc = 'C205 ψ(KC)'; }
      else if (inp.cdMode === 'c205' && KC === 0) { Cd = dragCoefficientSteady(inp.kRough / 1e3 / m.D); cdSrc = 'C205 steady'; }
      return { ...m, Uw, Uc, KC, Re, Cd, cdSrc, kD: 0.5 * rho * Cd * m.D * m.len * m.n * inp.shield, kI: rho * m.Cm * Math.PI * m.D * m.D / 4 * m.len * m.n * inp.shield * om * Uw, FDmax: 0, Fmax: 0 };
    });
    const n = inp.nPhase, ph = [], Fx = [], Fy = [], Fd = [], Fi = [], Fres = [], Mom = [];
    let iPk = 0;
    for (let i = 0; i < n; i++) {
      const th = 2 * Math.PI * i / n;
      let fx = 0, fy = 0, fd = 0, fi = 0, mo = 0;
      for (const m of items) {
        const ux = m.Uc * Math.cos(ang) + m.Uw * Math.cos(th), uy = m.Uc * Math.sin(ang), sp = Math.hypot(ux, uy);
        const dx = m.kD * sp * ux, dy = m.kD * sp * uy, I = -m.kI * Math.sin(th);
        const f = Math.hypot(dx + I, dy);
        m.Fmax = Math.max(m.Fmax, f); m.FDmax = Math.max(m.FDmax, m.kD * sp * sp);
        fx += dx + I; fy += dy; fd += dx; fi += I; mo += (dx + I) * m.z;
      }
      ph.push(round(th / DEG, 5)); Fx.push(fx / 1e3); Fy.push(fy / 1e3); Fd.push(fd / 1e3); Fi.push(fi / 1e3); Fres.push(Math.hypot(fx, fy) / 1e3); Mom.push(mo / 1e3);
      if (Fres[i] > Fres[iPk]) iPk = i;
    }
    const FH = Fres[iPk] * 1e3;
    const arm = Math.abs(Fx[iPk]) > 1e-9 ? Math.abs(Mom[iPk] / Fx[iPk]) : items.reduce((a, m) => a + m.z, 0) / items.length;
    const uL = ucAt(inp.zLift) + uwAt(inp.zLift);
    const FL = 0.5 * rho * inp.Cl * inp.Aplan * uL * uL;
    const FDpk = Math.max(...Fd.map(Math.abs)), FIpk = Math.max(...Fi.map(Math.abs));
    const zs = [], ucs = [], uws = [];
    for (let i = 1; i <= 40; i++) { const z = zTop * i / 40; zs.push(round(z)); ucs.push(ucAt(z)); uws.push(uwAt(z)); }
    const big = items.reduce((a, m) => (m.Fmax > a.Fmax ? m : a));
    const dragDominated = FDpk >= FIpk;
    const advice = [
      adv('engineering', 'info', `Peak horizontal load ${fmt(FH / 1e3)} kN acting ${fmt(arm)} m above the bed; ${dragDominated ? 'drag' : 'inertia'} dominates (drag ${fmt(FDpk)} kN, inertia ${fmt(FIpk)} kN). "${big.name}" contributes most (${fmt(big.Fmax / 1e3)} kN).`,
        dragDominated ? `Reduce projected area or streamline "${big.name}" first; a 10 % cut in its drag area removes about ${fmt(0.1 * big.FDmax / 1e3)} kN. Calibrate Cd and the shielding factor against flume or basin tests before final design.` : `Inertia scales with displaced plus entrained volume: open up closed boxes around "${big.name}" and confirm Cm (added mass) against DNV-RP-C205 App. D for the actual shapes.`,
        'DNV-RP-C205 (2019) §6.2, §6.7 and App. D; DNV-RP-F109 (2021) §3.4',
        'Lower hydrodynamic load means a smaller mat and less ballast: every kilonewton of peak load avoided on a friction base with tan δ = 0.5 and FS 1.5 saves about 0.3 t of submerged ballast.'),
      adv('economic', 'info', `The resultant ${fmt(FH / 1e3)} kN horizontal and ${fmt(FL / 1e3)} kN lift size the foundation; on a drained base at a sliding factor of 1.5 and tan δ = 0.5 they call for about ${fmt((1.5 * FH / 0.5 + FL) / 1e3 / G)} t of submerged weight.`,
        'Spend on a flume test or CFD study when the coefficients are uncertain: a 15 % reduction in Cd translates directly into foundation steel and lift weight. Carry the coefficient uncertainty as a range in the estimate.',
        'Value-of-information reasoning; AACE RP 18R-97 estimate classes',
        steelCo2Text(0.15 * (1.5 * FH / 0.5) / 1e3 / G)),
    ];
    return {
      kpis: [
        { key: 'FH', label: 'Peak resultant horizontal load', unit: 'kN', value: FH / 1e3, explain: 'Largest vector sum of drag and inertia on all items over the wave cycle. This is the load the foundation must resist in sliding.' },
        { key: 'arm', label: 'Lever arm of the peak load', unit: 'm', value: arm, explain: 'Height above the seabed at which the peak horizontal load acts; load × arm is the hydrodynamic overturning moment.' },
        { key: 'Mot', label: 'Overturning moment at peak load', unit: 'kN·m', value: Math.abs(Mom[iPk]), explain: 'Moment of the member loads about the seabed at the phase of peak horizontal load.' },
        { key: 'FL', label: 'Peak lift', unit: 'kN', value: FL / 1e3, explain: 'Vertical lift on the solid plan area at the peak velocity. It reduces the weight available for friction.' },
        { key: 'FDpk', label: 'Peak drag component', unit: 'kN', value: FDpk, explain: 'Largest total drag along the wave direction. Compare with the inertia component to see which coefficient matters.' },
        { key: 'FIpk', label: 'Peak inertia component', unit: 'kN', value: FIpk, explain: 'Largest total inertia (added-mass plus Froude–Krylov) load; it leads the drag peak by a quarter cycle.' },
        { key: 'Ubed', label: 'Peak velocity at the reference height', unit: 'm/s', value: ucAt(inp.zRef) + uwAt(inp.zRef), explain: 'Current plus orbital velocity amplitude at the current reference height: the basic environmental intensity near the bed.' },
        { key: 'KCmax', label: 'Largest Keulegan–Carpenter number', unit: '-', value: Math.max(...items.map((m) => m.KC)), explain: 'KC = Uw·T/D. Below about 5 the load is inertia dominated; above about 30 it is quasi-steady drag.' },
      ],
      checks: [],
      plots: [
        { type: 'line', title: 'Total load through the wave cycle', x: { label: 'Wave phase', unit: 'deg' }, y: { label: 'Force', unit: 'kN' },
          series: [{ name: 'Drag (along waves)', x: ph, y: Fd }, { name: 'Inertia', x: ph, y: Fi }, { name: 'Along waves', x: ph, y: Fx }, { name: 'Resultant', x: ph, y: Fres }],
          explain: 'Sum over all members at each wave phase. The resultant peak and the phase shift between drag and inertia show which component governs.' },
        { type: 'line', title: 'Near-bed velocity profile', x: { label: 'Velocity', unit: 'm/s' }, y: { label: 'Height above bed', unit: 'm' },
          series: [{ name: 'Current (log law)', x: ucs, y: zs }, { name: 'Wave orbital amplitude', x: uws, y: zs }, { name: 'Combined peak', x: ucs.map((u, i) => u + uws[i]), y: zs }],
          explain: 'Logarithmic boundary-layer current and linear-theory orbital amplitude against height. Members close to the bed see markedly less current.' },
        { type: 'bar', title: 'Peak load by item', categories: items.map((m) => m.name), series: [{ name: 'Peak load (kN)', y: items.map((m) => round(m.Fmax / 1e3)) }], y: { label: 'Force', unit: 'kN' },
          explain: 'Peak Morison load of each group of items. The tallest bar is where a change in geometry or coefficient pays most.' },
      ],
      tables: [{ title: 'Member loads', columns: [{ key: 'name', label: 'Item' }, { key: 'Uc', label: 'Current', unit: 'm/s' }, { key: 'Uw', label: 'Orbital', unit: 'm/s' }, { key: 'KC', label: 'KC' }, { key: 'Re', label: 'Re' }, { key: 'Cd', label: 'Cd used' }, { key: 'src', label: 'Cd source' }, { key: 'FD', label: 'Peak drag', unit: 'kN' }, { key: 'FI', label: 'Peak inertia', unit: 'kN' }, { key: 'F', label: 'Peak total', unit: 'kN' }],
        rows: items.map((m) => ({ name: m.name, Uc: round(m.Uc), Uw: round(m.Uw), KC: round(m.KC), Re: round(m.Re, 3), Cd: round(m.Cd), src: m.cdSrc, FD: round(m.FDmax / 1e3), FI: round(m.kI / 1e3), F: round(m.Fmax / 1e3) })),
        explain: 'One row per group of items: local velocities, flow regime numbers, the drag coefficient actually used and the peak load components (all members of the group together).' }],
      notes: [`Design wave H = ${fmt(H)} m, T = ${fmt(T)} s in ${fmt(inp.depth)} m of water.`, `${dragDominated ? 'Drag' : 'Inertia'} governs the resultant; the peak occurs at a phase of ${fmt(ph[iPk])}°.`],
      advice,
      publish: { 'subsea.hydroFH': FH, 'subsea.hydroFL': FL, 'subsea.hydroArm': arm },
      quantities: {},
    };
  },
  verify: [
    { name: 'Steady current on one cylinder: F = ½ρ·Cd·D·L·U²', source: 'Hand calculation: ½·1025·1.0·0.5·10·1² = 2562.5 N',
      inputs: { members: [{ name: 'c', D: 0.5, len: 10, z: 1, n: 1, Cd: 1.0, Cm: 2 }], Hs: 0, Uref: 1, zRef: 1, shield: 1, Cl: 0 },
      expect: [{ kpi: 'FH', value: 2.5625, tol: 1e-6 }, { kpi: 'arm', value: 1, tol: 1e-6 }] },
    { name: 'Waves only, inertia: F = ρ·Cm·(πD²/4)·L·ω·Uw with Uw = (πH/T)·cosh(kz)/sinh(kd)', source: 'Linear wave theory; k from ω² = g·k·tanh(kd) solved independently by fixed-point iteration',
      inputs: { members: [{ name: 'c', D: 1, len: 10, z: 2, n: 1, Cd: 0, Cm: 2 }], Hs: 2, hFac: 1, Tp: 8, Uref: 0, depth: 30, shield: 1, Cl: 0, nPhase: 72 },
      expect: [{ kpi: 'FH', value: (() => { const om = 2 * Math.PI / 8; let k = om * om / G; for (let i = 0; i < 200; i++) k = om * om / (G * Math.tanh(k * 30)); const Uw = Math.PI * 2 / 8 * Math.cosh(k * 2) / Math.sinh(k * 30); return 1025 * 2 * Math.PI / 4 * 10 * om * Uw / 1e3; })(), tol: 0.002 }] },
  ],
};
