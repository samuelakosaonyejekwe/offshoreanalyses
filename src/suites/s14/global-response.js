// Global static and dynamic response of a flexible riser: free-hanging catenary or lazy-wave configuration on
// the planar rod solver, near/far vessel offsets, current, and the response to harmonic vessel motion at the
// hang-off.
import { solveRiser, solveWaveRiser, catenaryRiser, riserOffsetCases, quasiStaticHeave } from '../../engine/lines/riser.js';
import { rodHarmonicResponse } from '../../engine/lines/rod.js';
import { equivalentBeamProperties } from '../../engine/xsection/flexible.js';
import { stackInputs, stackFrom, massSplit, muInput, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, thin, linspace, sev, worst, LIMIT_2D, UNBONDED, GENERIC, STD_17J, STD_17B, CONFIRM, CO2_STEEL, BAR, G } from './_common.js';

export default {
  id: 's14-global-response',
  suite: 14,
  title: 'Global static and dynamic response',
  summary: 'Solves the installed riser configuration and gives effective tension and curvature along it for static, offset and wave-frequency vessel-motion cases.',
  covers: [1],
  fidelity: 'numerical',
  standards: ['API RP 17B (5th ed.) §7 (riser configurations), §8', 'API Spec 17J (4th ed.) §5.3', 'DNV-RP-C205 §6 (Morison loading)', 'API RP 2RD §5 (load cases)'],
  models: ['Geometrically nonlinear beam / rod finite elements (corotational)', 'Nonlinear catenary start solution', 'Effective tension formulation', 'Morison drag with harmonic linearisation', 'Added mass', 'Seabed contact springs at touchdown', 'Vessel-motion response in the frequency domain'],
  limits: LIMIT_2D + ' Dynamics are a single-harmonic frequency-domain response about the static shape with linearised drag (or a quasi-static heave shift): no irregular-wave time domain, no wave kinematics on the riser, and bending hysteresis enters only through the constant slip stiffness. ' + UNBONDED + ' ' + GENERIC,
  publishes: ['flexible.topTension', 'riser.topTension', 'riser.hangoffAngleDeg', 'riser.tdpTension', 'riser.reactionH', 'flexible.maxCurvature', 'flexible.minBendRadius', 'flexible.tensionRange', 'flexible.curvatureRange', 'quantities.pipe_km'],
  searchTerms: 'flexible riser lazy wave global analysis effective tension curvature',
  inputs: [
    { key: 'propSource', label: 'Riser section properties from', type: 'select', default: 'stack', options: [{ value: 'stack', label: 'Layer stack (computed)' }, { value: 'datasheet', label: 'Manufacturer data sheet values' }], group: 'Riser section', role: 'geometry',
      help: 'Global analysis only needs axial and bending stiffness, weight and diameter; use the data sheet when it is available.' },
    { key: 'dsEA', label: 'Axial stiffness', unit: 'MN', type: 'number', default: 400, min: 1, max: 5000, group: 'Riser section', role: 'geometry', showIf: { key: 'propSource', in: ['datasheet'] }, bind: 'flexible.EA' },
    { key: 'dsEI', label: 'Bending stiffness', unit: 'kN·m²', type: 'number', default: 50, min: 0.5, max: 50000, group: 'Riser section', role: 'geometry', showIf: { key: 'propSource', in: ['datasheet'] }, bind: 'flexible.EIslip' },
    { key: 'dsW', label: 'Submerged weight, operating', unit: 'N/m', type: 'number', default: 650, min: 10, max: 20000, group: 'Riser section', role: 'load', showIf: { key: 'propSource', in: ['datasheet'] }, bind: 'flexible.submergedWeight' },
    { key: 'dsOD', label: 'Outer diameter', unit: 'mm', type: 'number', default: 236, min: 30, max: 1000, group: 'Riser section', role: 'geometry', showIf: { key: 'propSource', in: ['datasheet'] } },
    { key: 'dsMass', label: 'Mass in air incl. contents', unit: 'kg/m', type: 'number', default: 113, min: 1, max: 2000, group: 'Riser section', role: 'geometry', showIf: { key: 'propSource', in: ['datasheet'] } },
    ...stackInputs(),
    muInput('Riser section'),
    { key: 'rhoContent', label: 'Bore fluid density', unit: 'kg/m³', type: 'number', default: 700, min: 0, max: 2000, group: 'Internal fluid', role: 'load', uncertain: { dist: 'normal', cov: 0.1 } },
    { key: 'pInt', label: 'Internal pressure', unit: 'bar', type: 'number', default: 200, min: 0, max: 1500, group: 'Internal fluid', role: 'load', help: 'Used for the contact state of the layer stack (stiffness); effective tension is solved directly.' },
    { key: 'config', label: 'Configuration', type: 'select', default: 'lazy-wave', options: [{ value: 'lazy-wave', label: 'Lazy wave (distributed buoyancy)' }, { value: 'free-hanging', label: 'Free-hanging catenary' }], group: 'Installed configuration', role: 'ic' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 1000, min: 30, max: 3500, group: 'Installed configuration', role: 'environment', bind: 'site.depth' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Installed configuration', role: 'environment', bind: 'site.rhoW' },
    { key: 'zTop', label: 'Hang-off depth below the surface', unit: 'm', type: 'number', default: 15, min: 0, max: 200, group: 'Installed configuration', role: 'bc', help: 'The top end is pinned at this point (the bend stiffener is analysed in the end-fitting module).' },
    { key: 'hangAngle', label: 'Hang-off angle from vertical', unit: 'deg', type: 'number', default: 8, min: 1, max: 40, group: 'Installed configuration', role: 'ic', showIf: { key: 'config', in: ['free-hanging'] }, design: true, help: 'Sets the suspended length and lay-back of the free-hanging catenary.' },
    { key: 'grounded', label: 'Length on the seabed to the anchored end', unit: 'm', type: 'number', default: 300, min: 20, max: 5000, group: 'Installed configuration', role: 'bc', showIf: { key: 'config', in: ['free-hanging'] } },
    { key: 'upperLength', label: 'Upper catenary length', unit: 'm', type: 'number', default: 900, min: 20, max: 6000, group: 'Installed configuration', role: 'ic', showIf: { key: 'config', in: ['lazy-wave'] }, design: true },
    { key: 'buoyLength', label: 'Buoyancy section length', unit: 'm', type: 'number', default: 300, min: 10, max: 3000, group: 'Installed configuration', role: 'ic', showIf: { key: 'config', in: ['lazy-wave'] }, design: true },
    { key: 'lowerLength', label: 'Lower length incl. seabed part', unit: 'm', type: 'number', default: 600, min: 20, max: 6000, group: 'Installed configuration', role: 'ic', showIf: { key: 'config', in: ['lazy-wave'] } },
    { key: 'span', label: 'Horizontal distance hang-off to seabed end', unit: 'm', type: 'number', default: 1100, min: 20, max: 10000, group: 'Installed configuration', role: 'bc', showIf: { key: 'config', in: ['lazy-wave'] }, help: 'The seabed end is pinned at this distance.' },
    { key: 'buoyNet', label: 'Net uplift of the buoyancy section', unit: 'N/m', type: 'number', default: 500, min: 10, max: 10000, group: 'Installed configuration', role: 'load', showIf: { key: 'config', in: ['lazy-wave'] }, design: true, uncertain: { dist: 'normal', cov: 0.05 },
      help: 'Upward force per metre of the buoyant section after deducting the riser weight (allow for water absorption over life).' },
    { key: 'buoyOD', label: 'Buoyancy module outer diameter', unit: 'm', type: 'number', default: 0.6, min: 0.1, max: 3, group: 'Installed configuration', role: 'geometry', showIf: { key: 'config', in: ['lazy-wave'] } },
    { key: 'kSeabed', label: 'Seabed vertical stiffness', unit: 'kN/m/m', type: 'number', default: 100, min: 1, max: 1e5, group: 'Installed configuration', role: 'soil', help: 'Linear lift-off spring per metre of pipe.' },
    { key: 'Us', label: 'Surface current speed (in plane)', unit: 'm/s', type: 'number', default: 0.8, min: 0, max: 3, group: 'Hydrodynamic loading', role: 'bc', bind: 'site.Uc100', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'Ub', label: 'Near-bed current speed', unit: 'm/s', type: 'number', default: 0.1, min: 0, max: 2, group: 'Hydrodynamic loading', role: 'bc', bind: 'site.UcBed100' },
    { key: 'curExp', label: 'Current profile exponent', unit: '-', type: 'number', default: 1.5, min: 0.1, max: 8, group: 'Hydrodynamic loading', role: 'bc', help: 'U = Ub + (Us − Ub)·(height/depth)^n.' },
    { key: 'Cdn', label: 'Normal drag coefficient', unit: '-', type: 'number', default: 1.2, min: 0.3, max: 3, group: 'Hydrodynamic loading', role: 'load', uncertain: { dist: 'normal', cov: 0.15 } },
    { key: 'Ca', label: 'Added-mass coefficient', unit: '-', type: 'number', default: 1, min: 0, max: 2.5, group: 'Hydrodynamic loading', role: 'load' },
    { key: 'offset', label: 'Vessel offset, near and far', unit: 'm', type: 'number', default: 60, min: 0, max: 500, group: 'Vessel motion', role: 'bc', bind: 'vessel.offsetMax', help: 'In-plane hang-off excursion applied towards (near) and away from (far) the seabed end.' },
    { key: 'dynMethod', label: 'Dynamic method', type: 'select', default: 'harmonic', options: [{ value: 'harmonic', label: 'Harmonic response (frequency domain, drag linearised)' }, { value: 'quasi', label: 'Quasi-static heave' }], group: 'Vessel motion', role: 'control' },
    { key: 'heaveAmp', label: 'Hang-off vertical motion amplitude', unit: 'm', type: 'number', default: 2, min: 0, max: 15, group: 'Vessel motion', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'surgeAmp', label: 'Hang-off horizontal motion amplitude', unit: 'm', type: 'number', default: 1, min: 0, max: 15, group: 'Vessel motion', role: 'bc' },
    { key: 'period', label: 'Motion period', unit: 's', type: 'number', default: 10, min: 3, max: 30, group: 'Vessel motion', role: 'bc', bind: 'site.Tp100' },
    { key: 'allowTension', label: 'Allowable effective tension', unit: 'kN', type: 'number', default: 2400, min: 10, max: 50000, group: 'Acceptance criteria', role: 'criteria', bind: 'flexible.allowableTension', help: 'From the tension-and-bending module or the pipe data sheet.' },
    { key: 'mbr', label: 'Dynamic minimum bend radius', unit: 'm', type: 'number', default: 2.34, min: 0.3, max: 50, group: 'Acceptance criteria', role: 'criteria', bind: 'flexible.mbrDynamic', help: 'From the MBR module or the data sheet. ' + CONFIRM },
    { key: 'minClear', label: 'Minimum sag-bend clearance above seabed', unit: 'm', type: 'number', default: 10, min: 0, max: 500, group: 'Acceptance criteria', role: 'criteria', showIf: { key: 'config', in: ['lazy-wave'] } },
    { key: 'centreline', label: 'Imported centreline for comparison (optional)', type: 'geometry', kinds: ['polyline'], default: null, group: 'Installed configuration', role: 'geometry', help: 'A surveyed or third-party configuration (x horizontal, z up) overlaid on the computed shape.' },
    { key: 'nElem', label: 'Rod elements', type: 'int', default: 120, min: 30, max: 600, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Rod elements', levels: [60, 120, 240], qoi: 'topTension' },
  compute(inp, ctx) {
    const rho = inp.rhoW, depth = inp.depth;
    need(inp.zTop < depth - 10, 'The hang-off must be at least 10 m above the seabed.');
    let sec, od, steelT = 0, polyT = 0;
    if (inp.propSource === 'datasheet') {
      od = inp.dsOD / 1e3;
      sec = { EA: inp.dsEA * 1e6, EI: inp.dsEI * 1e3, w: inp.dsW, m: inp.dsMass };
    } else {
      const S = stackFrom(inp, ctx);
      const pMid = rho * G * depth / 2;
      const beam = run(() => equivalentBeamProperties(S.stack, { loads: { tension: 3e5, pInt: inp.pInt * BAR, pExt: pMid }, mu: inp.mu, weight: { rhoContent: inp.rhoContent, rhoSeawater: rho } }), 'Equivalent beam properties');
      od = S.od;
      sec = { EA: beam.EA, EI: beam.EIslip, w: beam.mass.submergedOperating, m: beam.mass.massOperating };
      const ms = massSplit(S.stack); steelT = ms.steel; polyT = ms.polymer;
    }
    need(sec.w > 0, `The riser section has a submerged weight of ${fmt(sec.w)} N/m — it floats; a catenary configuration is impossible. Check the contents density and layer stack.`);
    const pipe = { ...sec, D: od, Cdn: inp.Cdn, Cdt: 0.01, ma: inp.Ca * rho * Math.PI / 4 * od * od };
    const zs = linspace(-depth, 0, 11);
    const current = inp.Us > 0 || inp.Ub > 0 ? { rho, z: zs, U: zs.map((z) => inp.Ub + (inp.Us - inp.Ub) * ((z + depth) / depth) ** inp.curExp) } : undefined;
    const seabed = { z: -depth, kn: inp.kSeabed * 1e3 };
    const top = { x: 0, z: -inp.zTop };
    const lazy = inp.config === 'lazy-wave';
    let base, total;
    if (lazy) {
      total = inp.upperLength + inp.buoyLength + inp.lowerLength;
      const chord = Math.hypot(inp.span, depth - inp.zTop);
      need(total > chord * 1.01, `The riser (${fmt(total)} m) is too short to reach from the hang-off to the seabed end (${fmt(chord)} m straight line). Lengthen it or reduce the horizontal distance.`);
      need(total < inp.span + (depth - inp.zTop) + 0.9 * inp.buoyLength, `The riser (${fmt(total)} m) is too long for a ${fmt(inp.span)} m horizontal distance: the slack would pile up. Shorten it or move the seabed end out.`);
      const buoyant = { ...pipe, w: -inp.buoyNet, D: inp.buoyOD, m: sec.m + 0.45 * rho * Math.PI / 4 * (inp.buoyOD ** 2 - od ** 2), ma: inp.Ca * rho * Math.PI / 4 * inp.buoyOD ** 2 };
      base = run(() => solveWaveRiser({ pipe, buoyant, upperLength: inp.upperLength, buoyancyLength: inp.buoyLength, lowerLength: inp.lowerLength, top, bottom: { x: inp.span, z: -depth }, seabed, current, nElem: inp.nElem }), 'Lazy-wave static solution (if it does not converge, reduce the current or adjust the length split and buoyancy)');
    } else {
      const cat = run(() => catenaryRiser({ h: depth - inp.zTop, w: sec.w, hangoffAngleDeg: inp.hangAngle }), 'Catenary layout');
      total = cat.suspendedLength + inp.grounded;
      base = run(() => solveRiser({ sections: [{ ...pipe, length: total }], top, bottom: { x: cat.horizontalSpan + inp.grounded, z: -depth }, seabed, current, nElem: inp.nElem }), 'Free-hanging static solution (if it does not converge, reduce the current or increase the hang-off angle)');
    }
    need(base.converged, 'The static riser solution did not converge. Check lengths, buoyancy and the horizontal distance — the configuration may be geometrically impossible.');
    const cases = [{ name: 'Mean position', offset: 0, r: base }];
    if (inp.offset > 0) {
      const off = run(() => riserOffsetCases(base, [-inp.offset, inp.offset]), 'Vessel-offset cases');
      cases.push({ name: 'Near offset', offset: -inp.offset, r: off[0].result, tdpX: off[0].tdpX }, { name: 'Far offset', offset: inp.offset, r: off[1].result, tdpX: off[1].tdpX });
    }
    for (const c of cases) need(c.r.converged !== false, `${c.name}: the static solution did not converge — reduce the offset.`);
    // dynamic response about the mean position
    const n = base.s.length;
    let tAmp = new Array(n).fill(0), kAmp = new Array(n).fill(0), dynNote;
    if (inp.heaveAmp > 0 || inp.surgeAmp > 0) {
      if (inp.dynMethod === 'harmonic') {
        const hr = run(() => rodHarmonicResponse(base, { omega: 2 * Math.PI / inp.period, end: 'A', xAmp: inp.surgeAmp, zAmp: inp.heaveAmp, rho }), 'Harmonic response');
        tAmp = Array.from(hr.tensionAmp); kAmp = Array.from(hr.curvatureAmp);
        // the end elements of a pinned rod carry a numerical curvature spike: exclude 2 % of the length at each end
        const skip = Math.max(2, Math.round(0.02 * n));
        for (let i = 0; i < n; i++) if (i < skip || i >= n - skip) kAmp[i] = 0;
        dynNote = `Harmonic response at T = ${fmt(inp.period)} s ${hr.converged ? 'converged' : 'did not fully converge'} in ${hr.iterations} drag iterations.`;
      } else {
        const q = run(() => quasiStaticHeave(base, { heaveAmp: Math.max(inp.heaveAmp, 1e-3) }), 'Quasi-static heave');
        tAmp = tAmp.map(() => q.tensionAmplitude);
        kAmp = kAmp.map(() => Math.max(0, q.maxCurvature - base.maxCurvature));
        dynNote = 'Quasi-static heave: inertia and drag ignored (lower bound for fast motions).';
      }
    } else dynNote = 'No vessel motion applied: static results only.';
    for (const v of tAmp) need(Number.isFinite(v), 'The dynamic response is not finite — try the quasi-static method or a different motion period.');
    const Te0 = Array.from(base.Te), k0 = Array.from(base.curvature, Math.abs);
    const tAmpTop = tAmp[0];
    const maxT = Math.max(...cases.map((c) => c.r.maxTension)) + Math.max(...tAmp);
    let minT = Infinity, ampAtMin = 0, staticMin = Infinity;
    for (let i = 0; i < n; i++) { const v = Te0[i] - tAmp[i]; if (v < minT) { minT = v; ampAtMin = tAmp[i]; staticMin = Te0[i]; } }
    minT = Math.min(minT, ...cases.map((c) => c.r.minTension));
    const kDyn = Math.max(...kAmp);
    const maxK = Math.max(...cases.map((c) => c.r.maxCurvature)) + kDyn;
    need(maxK > 0, 'The riser is perfectly straight — check the configuration inputs.');
    const checks = [
      chk('Maximum effective tension', maxT / 1e3, inp.allowTension, 'kN', `${STD_17J} Table 6 via the allowable tension of the cross-section`, 'Guards against overstressing the tensile armour at the most loaded point (normally the hang-off) in the worst offset with vessel motion.'),
      chk('Minimum bend radius (static + dynamic curvature)', maxK, 1 / inp.mbr, '1/m', `${STD_17J} §5.3.1 (dynamic MBR)`, 'Guards against over-bending at the sag bend, hog bend or touchdown. Curvature above 1/MBR overstrains the sheaths and can lock the interlocked layers.'),
      chk('No effective compression (dynamic amplitude vs static tension)', ampAtMin, Math.max(staticMin, 1e-6), 'N', `${STD_17B} §8 (bird-caging; compression to be avoided or assessed)`, 'Guards against the riser going into effective compression during the motion cycle, which can buckle the line and bird-cage the armour.'),
    ];
    if (lazy) checks.push(chk('Sag-bend clearance above the seabed', inp.minClear, Math.max(Math.min(...cases.map((c) => c.r.sagClearance ?? base.sagClearance)), 1e-6), 'm', `${STD_17B} §7 (configuration design: clearance)`,
      'Guards against the sag bend touching the seabed in the near position, which abrades the outer sheath and changes the configuration.'));
    const gov = worst(checks);
    const tdpX = cases.map((c) => c.r.tdp?.x).filter(Number.isFinite);
    const tdpExc = tdpX.length > 1 ? Math.max(...tdpX) - Math.min(...tdpX) : 0;
    const shapeSeries = cases.map((c) => { const [x, z] = thin([c.r.x, c.r.z], 200); return { name: c.name, x: x.map((v) => round(v)), y: z.map((v) => round(v)) }; });
    const cl = inp.centreline;
    if (cl && cl.points && cl.points.length >= 6) {
      const np = cl.points.length / 3, px = [], pz = [];
      for (let i = 0; i < np; i++) { px.push(cl.points[3 * i]); pz.push(cl.points[3 * i + 2]); }
      const [x, z] = thin([px, pz], 300);
      shapeSeries.push({ name: cl.name || 'Imported centreline', x: x.map((v) => round(v)), y: z.map((v) => round(v)), dash: true });
    }
    const [sT, T1, Tup, Tlo] = thin([base.s, Te0, Te0.map((v, i) => v + tAmp[i]), Te0.map((v, i) => v - tAmp[i])], 250);
    const [sK, K1, Kup] = thin([base.s, k0, k0.map((v, i) => v + kAmp[i])], 250);
    const km = total / 1e3;
    const ucT = checks[0].uc, ucK = checks[1].uc;
    return {
      kpis: [
        kpi('topTension', 'Static top effective tension', base.topTension / 1e3, 'kN', 'Effective tension at the hang-off in the mean position. It is the suspended weight less buoyancy, plus the horizontal pull.'),
        kpi('maxTension', 'Maximum effective tension', maxT / 1e3, 'kN', 'Largest tension over all offset cases with the dynamic amplitude added.', inp.allowTension, 'max'),
        kpi('minTension', 'Minimum effective tension', minT / 1e3, 'kN', 'Smallest tension along the riser including the dynamic amplitude. Negative means effective compression.', 0, 'min'),
        kpi('hangoffAngle', 'Hang-off angle from vertical', base.hangoffAngleDeg, 'deg', 'Static top angle in the mean position; with the offset and motion ranges it sizes the bend stiffener.'),
        kpi('maxCurvature', 'Maximum curvature', maxK, '1/m', 'Largest static curvature of any case plus the dynamic curvature amplitude.', 1 / inp.mbr, 'max'),
        kpi('minBendRadius', 'Minimum bend radius reached', 1 / maxK, 'm', 'Tightest radius along the riser; must stay above the dynamic MBR.', inp.mbr, 'min'),
        kpi('tdpTension', 'Tension at the seabed end', base.bottomTension / 1e3, 'kN', 'Effective tension delivered to the seabed end or flowline — the horizontal pull on the riser base.'),
        kpi('tensionRange', 'Top tension range from vessel motion', 2 * tAmpTop / 1e3, 'kN', 'Peak-to-peak dynamic tension at the hang-off for the applied harmonic motion; input to armour fatigue.'),
        kpi('curvatureRange', 'Largest curvature range from vessel motion', 2 * kDyn, '1/m', 'Peak-to-peak dynamic curvature away from the end fittings; input to armour fatigue.'),
        kpi('tdpExcursion', 'Touchdown excursion near to far', tdpExc, 'm', 'Travel of the touchdown point between the near and far vessel positions; it sets the length of abrasion protection and the trench extent.'),
        kpi('sagClearance', 'Sag-bend clearance (mean position)', lazy ? base.sagClearance : 0, 'm', 'Height of the sag bend above the seabed for a lazy wave (zero for a free-hanging catenary).'),
        kpi('hogElevation', 'Hog-bend height above seabed', lazy ? base.hogElevation : 0, 'm', 'Height of the top of the wave above the seabed (zero for a free-hanging catenary).'),
        kpi('length', 'Riser length', total, 'm', 'Total length from hang-off to the seabed end.'),
        kpi('maxUc', 'Governing global utilisation', gov.uc, '-', `Highest utilisation of the global checks (${gov.name}).`, 1, 'max'),
      ],
      checks,
      plots: [
        plot('line', 'Riser configuration', 'Horizontal distance', 'm', 'Elevation', 'm', shapeSeries, 'Shape of the riser in the mean, near and far vessel positions. For a lazy wave look at the sag-bend clearance in the near case and the hog height in the far case.'),
        plot('line', 'Effective tension along the riser', 'Arc length from hang-off', 'm', 'Effective tension', 'kN', [
          { name: 'Static (mean)', x: sT.map((v) => round(v)), y: T1.map((v) => round(v / 1e3)) }, { name: 'Static + dynamic', x: sT.map((v) => round(v)), y: Tup.map((v) => round(v / 1e3)) }, { name: 'Static − dynamic', x: sT.map((v) => round(v)), y: Tlo.map((v) => round(v / 1e3)) }],
        'Tension from the hang-off down to the seabed end with the dynamic envelope. The minimum must stay positive; the maximum is compared with the allowable tension.'),
        plot('line', 'Curvature along the riser', 'Arc length from hang-off', 'm', 'Curvature', '1/m', [
          { name: 'Static (mean)', x: sK.map((v) => round(v)), y: K1.map((v) => round(v)) }, { name: 'Static + dynamic amplitude', x: sK.map((v) => round(v)), y: Kup.map((v) => round(v)) }, { name: '1 / MBR', x: [0, round(total)], y: [round(1 / inp.mbr), round(1 / inp.mbr)], dash: true }],
        'Curvature along the riser. Peaks mark the sag bend, hog bend and touchdown; none may reach the dashed 1/MBR line.'),
      ],
      tables: [tbl('Static cases', [col('case', 'Case'), col('offset', 'Offset', 'm'), col('top', 'Top tension', 'kN'), col('min', 'Minimum tension', 'kN'), col('angle', 'Hang-off angle', 'deg'), col('kappa', 'Max curvature', '1/m'), col('mbr', 'Min bend radius', 'm'), col('tdp', 'Touchdown x', 'm')],
        cases.map((c) => ({ case: c.name, offset: c.offset, top: round(c.r.topTension / 1e3), min: round(c.r.minTension / 1e3), angle: round(c.r.hangoffAngleDeg), kappa: round(c.r.maxCurvature), mbr: round(1 / Math.max(c.r.maxCurvature, 1e-9)), tdp: c.r.tdp ? round(c.r.tdp.x) : null })),
        'Static results for the mean position and the near and far vessel offsets (in-plane only).')],
      notes: [`${gov.name} governs (${fmt(gov.uc)}).`, dynNote, `Section used: EA ${fmt(sec.EA / 1e6)} MN, EI ${fmt(sec.EI / 1e3)} kN·m², submerged weight ${fmt(sec.w)} N/m (${inp.propSource === 'stack' ? 'from the layer stack' : 'data sheet'}).`, LIMIT_2D],
      advice: [
        adv('engineering', sev(gov.uc), `Maximum tension ${fmt(maxT / 1e3)} kN (${fmt(ucT * 100)} % of allowable); tightest radius ${fmt(1 / maxK)} m against an MBR of ${fmt(inp.mbr)} m; minimum tension ${fmt(minT / 1e3)} kN.`,
          gov.uc > 1 ? (minT < 0 ? 'Remove the compression: move buoyancy towards the touchdown, lengthen the lower catenary or limit the vessel motion.' : (ucK > ucT ? `Open the tight bend: ${lazy ? 'lengthen the buoyancy section or reduce its uplift by about ' + fmt((1 - 1 / ucK) * 100) + ' %' : 'increase the hang-off angle or lay-back'}.` : 'Reduce top tension with more distributed buoyancy or a lighter (gas-filled check) contents case.')) : 'Run the near, far and (in a 3-D tool) cross cases for the full load-case matrix, and pass the tension and curvature ranges to the armour-fatigue module.',
          `${STD_17B} §7–8; ${STD_17J} §5.3; DNV-RP-C205 §6`),
        adv('economic', 'info', `The configuration uses ${fmt(total)} m of riser${lazy ? ' and ' + fmt(inp.buoyLength) + ' m of buoyancy modules providing about ' + fmt((inp.buoyNet + sec.w) * inp.buoyLength / 1e3) + ' kN of gross lift' : ''}.`,
          lazy ? `Top tension is ${fmt(base.topTension / 1e3)} kN; each 10 % less buoyancy saves about ${fmt(0.1 * (inp.buoyNet + sec.w) * inp.buoyLength / G / 1e3 / 0.55)} t of syntactic foam but raises curvature and seabed contact — optimise buoyancy against the clearance check rather than adding margin everywhere.` : `A free-hanging catenary avoids buoyancy modules altogether; check it stays compression-free in the near position before paying for a wave configuration.`,
          'Concept selection and quantity-based estimating per AACE RP 18R-97 (Class 4); whole-life costing per ISO 15686-5',
          steelT > 0 ? `The riser contains about ${fmt(steelT * km)} t of steel ≈ ${fmt(steelT * km * CO2_STEEL)} t CO₂e (screening factor ${CO2_STEEL} t CO₂e per tonne) and ${fmt(polyT * km)} t of polymer; every 100 m of length saved removes ${fmt(steelT * 0.1 * CO2_STEEL)} t CO₂e.` : `Riser length drives embodied carbon: every 100 m less is roughly ${fmt(sec.m * 0.1 * 0.7 * CO2_STEEL)} t CO₂e at an assumed 70 % steel content and ${CO2_STEEL} t CO₂e per tonne (screening).`),
      ],
      publish: { 'flexible.topTension': maxT, 'riser.topTension': base.topTension, 'riser.hangoffAngleDeg': base.hangoffAngleDeg, 'riser.tdpTension': base.bottomTension, 'riser.reactionH': base.bottomTension,
        'flexible.maxCurvature': maxK, 'flexible.minBendRadius': 1 / maxK, 'flexible.tensionRange': 2 * tAmpTop, 'flexible.curvatureRange': 2 * kDyn, 'quantities.pipe_km': km },
      quantities: steelT > 0 ? { pipe_km: km, steel_t: steelT * km, polymer_t: polyT * km } : { pipe_km: km },
    };
  },
  verify: [{
    // Free-hanging inextensible catenary: vertical equilibrium T·cosθ = w·S and T = H + w·h give T_top = w·h/(1 − sinθ)
    // and H = T·sinθ (θ from vertical). w = 600 N/m, h = 1000 m, θ = 10°: T = 600 000/(1 − 0.173648) = 726.08 kN,
    // H = 126.08 kN. Bending stiffness 10 kN·m² and no current, so the rod solution should match the cable.
    name: 'Free-hanging catenary: top tension and horizontal tension',
    source: 'Irvine, Cable Structures, §1.2; Sparks (2007) catenary equations; hand calculation',
    inputs: { propSource: 'datasheet', dsEA: 400, dsEI: 10, dsW: 600, dsOD: 236, dsMass: 113, config: 'free-hanging', depth: 1000, zTop: 0, hangAngle: 10, grounded: 300, Us: 0, Ub: 0, offset: 0, heaveAmp: 0, surgeAmp: 0, nElem: 240 },
    expect: [{ kpi: 'topTension', value: 726.08, tol: 0.01 }, { kpi: 'tdpTension', value: 126.08, tol: 0.02 }],
  }],
};
