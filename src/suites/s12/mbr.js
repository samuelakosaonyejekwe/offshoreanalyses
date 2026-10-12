// Minimum bend radius and overbending: storage, installation and operating MBR from the stress capacity of the
// cross-section, bend-stiffener response at the hang-off, bend-restrictor locking radius at the static termination.
import { umbilicalMbr } from '../../engine/xsection/umbilical.js';
import { bendStiffener, bendRestrictor } from '../../engine/xsection/ancillaries.js';
import { sectionInputs, criteriaInputs, buildSection, adv, check, worst, sev, fmt, round, run, need, linspace, thin, G, DEG, CO2_POLYMER, LIMIT_RIGID, STD_UMB, CONFIRM } from './_common.js';

export default {
  id: 's12-minimum-bend-radius', suite: 12, title: 'Minimum bend radius and overbending',
  summary: 'Allowable bend radii for storage, installation and operation, and the curvature actually reached at the reel, chute, bend stiffener and bend restrictor.',
  covers: [4], fidelity: 'reduced-order',
  standards: [`${STD_UMB} §7.3 (MBR for storage, installation and operation)`, 'API RP 17B (2014) §4 (bend stiffeners and restrictors)', 'Boef & Out (1990) OTC 6436 (bend-stiffener beam model)'],
  models: ['Minimum bend-radius relationship from the von Mises capacity of tubes and wires', 'Elastic rod theory: planar elastica with varying EI (bend stiffener)', 'Bend-restrictor locking geometry', 'Nonlinear bending with tension stiffening'],
  limits: `${LIMIT_RIGID} MBR is an elastic stress criterion: acceptance of plastic reeling strain in steel tubes is a manufacturer qualification and is not modelled. The stiffener is a linear-elastic polyurethane cone without a gap; its modulus is temperature and rate dependent.`,
  publishes: ['umbilical.mbrDynamic', 'umbilical.mbrStorage', 'umbilical.mbrInstall', 'umbilical.hangoffMoment', 'umbilical.hangoffShear'],
  searchTerms: 'umbilical minimum bend radius bend stiffener design',
  inputs: [
    ...sectionInputs(),
    { key: 'Top', label: 'Operating tension at the hang-off', unit: 'kN', type: 'number', default: 120, min: 1, max: 5000, group: 'Loading', role: 'bc', bind: 'umbilical.topTension', uncertain: { dist: 'lognormal', cov: 0.1 } },
    { key: 'Tinst', label: 'Installation tension over the chute', unit: 'kN', type: 'number', default: 100, min: 0, max: 5000, group: 'Loading', role: 'bc', bind: 'umbilical.layTension' },
    { key: 'angle', label: 'Hang-off angle deviation (static + dynamic)', unit: 'deg', type: 'number', default: 12, min: 0.5, max: 60, group: 'Loading', role: 'bc', uncertain: { dist: 'normal', cov: 0.15 },
      help: 'Angle between the line tension and the stiffener root axis in the extreme case (vessel pitch/roll plus line angle change).' },
    { key: 'zSection', label: 'Water depth at the hang-off', unit: 'm', type: 'number', default: 10, min: 0, max: 3500, group: 'Loading', role: 'bc' },
    { key: 'globalMbr', label: 'Tightest radius from the global analysis', unit: 'm', type: 'number', default: 25, min: 0.5, max: 1e5, group: 'Loading', role: 'bc', bind: 'umbilical.minBendRadius', help: 'Sag, hog or touchdown bend from the static / dynamic configuration.' },
    { key: 'reelR', label: 'Storage reel / carousel hub radius', unit: 'm', type: 'number', default: 6.5, min: 0.5, max: 30, group: 'Handling radii', role: 'bc', design: true },
    { key: 'chuteR', label: 'Overboarding chute radius', unit: 'm', type: 'number', default: 9, min: 0.5, max: 30, group: 'Handling radii', role: 'bc', design: true },
    { key: 'stiffL', label: 'Bend-stiffener length', unit: 'm', type: 'number', default: 4.5, min: 0.3, max: 15, group: 'Bend stiffener', role: 'geometry', design: true },
    { key: 'rootD', label: 'Stiffener root diameter', unit: 'm', type: 'number', default: 0.7, min: 0.03, max: 3, group: 'Bend stiffener', role: 'geometry', design: true },
    { key: 'tipD', label: 'Stiffener tip diameter', unit: 'm', type: 'number', default: 0.19, min: 0.03, max: 3, group: 'Bend stiffener', role: 'geometry' },
    { key: 'puMat', label: 'Stiffener material', type: 'select', source: 'materials:polymer', default: 'pu-stiffener', group: 'Bend stiffener', role: 'material' },
    { key: 'puStrain', label: 'Allowable polyurethane strain', unit: '%', type: 'number', default: 8, min: 1, max: 20, group: 'Bend stiffener', role: 'criteria', help: 'Typical design values are a few per cent; take it from the stiffener supplier.' },
    { key: 'eiBasis', label: 'Umbilical bending stiffness inside the stiffener', type: 'select', default: 'slip', group: 'Bend stiffener', role: 'material', options: [{ value: 'slip', label: 'Slip stiffness (large curvature, conservative for curvature)' }, { value: 'stick', label: 'Stick stiffness (small cycles)' }] },
    { key: 'resPitch', label: 'Bend-restrictor element pitch', unit: 'm', type: 'number', default: 0.35, min: 0.05, max: 2, group: 'Bend restrictor (static end)', role: 'geometry' },
    { key: 'resLock', label: 'Locking angle per element', unit: 'deg', type: 'number', default: 2.0, min: 0.2, max: 15, group: 'Bend restrictor (static end)', role: 'geometry', design: true },
    { key: 'resN', label: 'Number of restrictor elements', type: 'int', default: 30, min: 2, max: 200, group: 'Bend restrictor (static end)', role: 'geometry' },
    ...criteriaInputs(),
    { key: 'etaStorage', label: 'Utilisation factor — storage (unpressurised, no tension)', unit: '-', type: 'number', default: 1.0, min: 0.3, max: 1, group: 'Acceptance criteria', role: 'criteria', help: `Storage on the reel is usually allowed to approach yield in bending. ${CONFIRM}` },
    { key: 'nSteps', label: 'Integration steps along the stiffener', type: 'int', default: 200, min: 20, max: 4000, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nSteps', label: 'Integration steps', levels: [50, 100, 200, 400], qoi: 'stiffKmax' },
  compute(inp, ctx) {
    const b = buildSection(inp, ctx), s = b.sec;
    need(s.elements.some((e) => e.component.yieldStress), 'No component has a yield stress: the bend radius cannot be derived. Enter the yield of the tubes or wires.');
    const b0 = buildSection(inp, ctx, { pressureScale: 0 });        // unpressurised section for storage on the reel
    const pExt = inp.rhoW * G * inp.zSection, mu = inp.muInt;
    const mbrOf = (sec, T, eta, pe) => run(() => umbilicalMbr(sec, { tension: T, eta, mu, pExt: pe }), 'Minimum bend radius').mbr;
    const mbrSt = mbrOf(b0.sec, 0, inp.etaStorage, 0), mbrIn = mbrOf(s, inp.Tinst * 1e3, inp.etaAbnormal, 0), mbrOp = mbrOf(s, inp.Top * 1e3, inp.etaNormal, pExt);
    for (const [n, v] of [['storage', mbrSt], ['installation', mbrIn], ['operation', mbrOp]]) need(Number.isFinite(v), `No bend radius is admissible in ${n}: the section is already over-utilised by tension and pressure alone. Reduce the tension or raise the utilisation factor.`);
    // bend stiffener at the hang-off
    const pu = ctx.material(inp.puMat); need(pu && pu.E, `Stiffener material '${inp.puMat}' has no modulus in the materials library.`);
    need(inp.tipD >= b.OD, `Stiffener tip diameter (${inp.tipD} m) is smaller than the umbilical (${fmt(b.OD, 4)} m).`);
    need(inp.rootD >= inp.tipD, 'Stiffener root diameter must not be smaller than its tip diameter.');
    const EIp = inp.eiBasis === 'slip' ? s.EIslip : s.EIstick;
    const bsArgs = { tension: inp.Top * 1e3, length: inp.stiffL, rootDiameter: inp.rootD, tipDiameter: inp.tipD, innerDiameter: b.OD, Epu: pu.E, EIpipe: EIp, nSteps: inp.nSteps, mbr: mbrOp, strainLimit: inp.puStrain / 100 };
    const bs = run(() => bendStiffener({ ...bsArgs, angle: inp.angle * DEG }), 'Bend stiffener');
    // unprotected exit for comparison: root curvature of a tensioned beam, κ = 2·sqrt(T/EI)·sin(θ/2) (elastica first integral)
    const kBare = 2 * Math.sqrt(inp.Top * 1e3 / EIp) * Math.sin(inp.angle * DEG / 2);
    const res = run(() => bendRestrictor({ pitch: inp.resPitch, lockAngle: inp.resLock * DEG, nElements: inp.resN }), 'Bend restrictor');
    const mbrStatic = mbrOf(s, 0, inp.etaNormal, inp.rhoW * G * 0);
    const checks = [
      check('Storage: reel hub radius', mbrSt, inp.reelR, 'm', `${STD_UMB} §7.3 (storage MBR)`, 'Demand is the storage MBR, capacity the hub radius. Winding tighter than the MBR yields the tubes in bending and can ovalise or kink them.'),
      check('Installation: chute radius', mbrIn, inp.chuteR, 'm', `${STD_UMB} §7.3 (installation MBR at lay tension)`, 'Demand is the MBR at the lay tension, capacity the chute radius. Overbending under tension at the chute is the classic installation damage mode.'),
      check('Operation: curvature in the bend stiffener', bs.maxCurvature, 1 / mbrOp, '1/m', `${STD_UMB} §7.3 (dynamic MBR); API RP 17B §4`, 'Largest curvature of the umbilical inside or just beyond the stiffener against the allowable dynamic curvature. Exceeding it overstresses tubes at the most fatigue-sensitive location.'),
      check('Bend-stiffener polyurethane strain', bs.maxPuStrain * 100, inp.puStrain, '%', 'API RP 17B (2014) §4 / supplier design strain', 'Outer-fibre strain of the polyurethane cone. Overstrain cracks the stiffener at its root, after which it no longer protects the umbilical.'),
      check('Operation: global configuration radius', mbrOp, inp.globalMbr, 'm', `${STD_UMB} §7.3`, 'Demand is the dynamic MBR, capacity the tightest radius of the sag, hog or touchdown bend from the global analysis.'),
      check('Static end: bend-restrictor locking radius', mbrStatic, res.lockingRadius, 'm', 'API RP 17B (2014) §4 (bend restrictors)', 'Demand is the static MBR, capacity the radius at which the restrictor elements lock. A restrictor that locks tighter than the MBR does not protect the umbilical at the termination.'),
    ];
    const gov = worst(checks);
    const angles = linspace(Math.max(1, inp.angle * 0.2), Math.min(60, inp.angle * 1.6), 9);
    const kOfAngle = angles.map((a) => { try { return bendStiffener({ ...bsArgs, nSteps: Math.min(inp.nSteps, 120), angle: a * DEG }).maxCurvature; } catch { return null; } });
    const tens = linspace(0, Math.max(inp.Top, inp.Tinst) * 2, 8);
    const mbrCurve = (eta, pe) => tens.map((t) => { const v = umbilicalMbr(s, { tension: t * 1e3, eta, mu, pExt: pe }).mbr; return Number.isFinite(v) ? round(v) : null; });
    const [ss, kk] = thin([bs.s, bs.curvature], 200);
    // polyurethane mass of the cone (frustum minus bore)
    const puVol = Math.PI * inp.stiffL / 12 * (inp.rootD ** 2 + inp.rootD * inp.tipD + inp.tipD ** 2) - Math.PI * b.OD ** 2 / 4 * inp.stiffL;
    const puT = puVol * (pu.rho ?? 1150) / 1e3;
    return {
      kpis: [
        { key: 'mbrDynamic', label: 'MBR — operation (dynamic)', value: round(mbrOp), unit: 'm', explain: `Tightest radius at which the governing tube or wire reaches ${inp.etaNormal} × yield at the operating tension and pressure. Every operating bend must stay above it.` },
        { key: 'mbrInstall', label: 'MBR — installation', value: round(mbrIn), unit: 'm', explain: `Tightest radius at the lay tension with the abnormal-condition factor ${inp.etaAbnormal}; it sizes the chute, lay wheel and tensioner exit.` },
        { key: 'mbrStorage', label: 'MBR — storage', value: round(mbrSt), unit: 'm', explain: `Tightest radius with no tension or pressure at factor ${inp.etaStorage}; it sizes the reel or carousel hub.` },
        { key: 'mbrOverD', label: 'Dynamic MBR / diameter', value: round(mbrOp / b.OD), unit: '-', explain: 'MBR as a multiple of overall diameter — a quick comparison with data-sheet values (steel-tube umbilicals are typically several tens of diameters).' },
        { key: 'stiffKmax', label: 'Maximum curvature at the hang-off', value: round(bs.maxCurvature, 6), unit: '1/m', limit: round(1 / mbrOp, 5), sense: 'max', explain: 'Largest curvature of the umbilical along the bend stiffener (or just beyond its tip) at the extreme angle. Its inverse is the operating bend radius at the hang-off.' },
        { key: 'stiffRmin', label: 'Minimum radius at the hang-off', value: round(1 / bs.maxCurvature), unit: 'm', limit: round(mbrOp), sense: 'min', explain: 'Inverse of the maximum curvature in the stiffener region; it must exceed the dynamic MBR.' },
        { key: 'bareRmin', label: 'Radius without a stiffener', value: round(1 / kBare), unit: 'm', explain: 'Radius the umbilical would take at a rigid clamped exit with no stiffener (κ = 2√(T/EI)·sin(θ/2)): shows how much the stiffener is needed.' },
        { key: 'rootMoment', label: 'Stiffener root moment', value: round(bs.rootMoment / 1e3), unit: 'kN·m', explain: 'Bending moment the stiffener flange and the I-tube or hang-off structure must react.' },
        { key: 'rootShear', label: 'Stiffener root shear', value: round(bs.rootShear / 1e3), unit: 'kN', explain: 'Lateral support reaction at the stiffener root, T·sin(θ).' },
        { key: 'puStrainMax', label: 'Maximum polyurethane strain', value: round(bs.maxPuStrain * 100), unit: '%', limit: inp.puStrain, sense: 'max', explain: 'Outer-fibre strain of the stiffener cone; governs stiffener fatigue and cracking.' },
        { key: 'lockR', label: 'Bend-restrictor locking radius', value: round(res.lockingRadius), unit: 'm', limit: round(mbrStatic), sense: 'min', explain: 'Radius at which the restrictor vertebrae lock (p / 2·sin(δ/2)); it must not be tighter than the static MBR.' },
        { key: 'utilMax', label: 'Governing bend-radius utilisation', value: round(gov.uc), unit: '-', limit: 1, sense: 'max', explain: 'Largest ratio of allowable to actual radius (or actual to allowable curvature) over storage, installation and operation.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Curvature along the bend stiffener', x: { label: 'Distance from the root', unit: 'm' }, y: { label: 'Curvature', unit: '1/m' },
          series: [{ name: 'Umbilical curvature', x: ss, y: kk }, { name: 'Allowable (1/MBR dynamic)', x: [0, inp.stiffL], y: [1 / mbrOp, 1 / mbrOp], dash: true }],
          explain: 'A well-proportioned stiffener spreads curvature evenly along its length. A peak at the tip means the stiffener is too short or too stiff; a peak at the root means it is too soft. The curve must stay below the dashed limit.' },
        { type: 'line', title: 'Maximum hang-off curvature vs angle', x: { label: 'Angle deviation', unit: 'deg' }, y: { label: 'Maximum curvature', unit: '1/m' },
          series: [{ name: 'With stiffener', x: angles, y: kOfAngle }, { name: 'Allowable', x: [angles[0], angles[angles.length - 1]], y: [1 / mbrOp, 1 / mbrOp], dash: true }],
          explain: 'How the curvature grows with the angle between the line and the hang-off axis at the operating tension. The crossing with the dashed line is the largest angle the stiffener can accept.' },
        { type: 'line', title: 'Minimum bend radius vs tension', x: { label: 'Effective tension', unit: 'kN' }, y: { label: 'MBR', unit: 'm' },
          series: [{ name: `Normal (η = ${inp.etaNormal})`, x: tens, y: mbrCurve(inp.etaNormal, pExt) }, { name: `Abnormal (η = ${inp.etaAbnormal})`, x: tens, y: mbrCurve(inp.etaAbnormal, 0) }, { name: `Accidental (η = ${inp.etaAccidental})`, x: tens, y: mbrCurve(inp.etaAccidental, 0) }],
          explain: 'The allowable radius grows with tension because the tension share uses up stress. Where a curve ends, no bending at all is admissible at that tension. Read the MBR for any handling step from the curve of its load condition.' },
        { type: 'bar', title: 'Bend-radius utilisation by stage', categories: checks.map((c) => c.name), series: [{ name: 'Utilisation [-]', y: checks.map((c) => round(c.uc)) }],
          explain: 'Allowable over actual radius for each stage of life. The tallest bar is where overbending is most likely; anything above 1.0 needs a larger radius or a lower tension.' },
      ],
      tables: [{ title: 'Minimum bend radii and actual radii', explain: 'One row per life stage: the tension and utilisation factor used, the resulting MBR and the radius actually imposed by the equipment or configuration. Margin = actual / allowable.',
        columns: [{ key: 'stage', label: 'Stage' }, { key: 'T', label: 'Tension', unit: 'kN' }, { key: 'eta', label: 'η', unit: '-' }, { key: 'mbr', label: 'MBR', unit: 'm' }, { key: 'actual', label: 'Actual radius', unit: 'm' }, { key: 'margin', label: 'Margin', unit: '-' }],
        rows: [
          { stage: 'Storage (reel hub)', T: 0, eta: inp.etaStorage, mbr: round(mbrSt), actual: inp.reelR, margin: round(inp.reelR / mbrSt) },
          { stage: 'Installation (chute)', T: inp.Tinst, eta: inp.etaAbnormal, mbr: round(mbrIn), actual: inp.chuteR, margin: round(inp.chuteR / mbrIn) },
          { stage: 'Operation (bend stiffener)', T: inp.Top, eta: inp.etaNormal, mbr: round(mbrOp), actual: round(1 / bs.maxCurvature), margin: round(1 / bs.maxCurvature / mbrOp) },
          { stage: 'Operation (global bends)', T: inp.Top, eta: inp.etaNormal, mbr: round(mbrOp), actual: inp.globalMbr, margin: round(inp.globalMbr / mbrOp) },
          { stage: 'Static end (restrictor)', T: 0, eta: inp.etaNormal, mbr: round(mbrStatic), actual: round(res.lockingRadius), margin: round(res.lockingRadius / mbrStatic) },
        ] }],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        `The stiffener raises the hang-off radius from ${fmt(1 / kBare)} m (bare exit) to ${fmt(1 / bs.maxCurvature)} m; the maximum is ${fmt(bs.maxCurvatureAt)} m from the root.`,
        `Storage MBR uses the unpressurised section at η = ${inp.etaStorage}. ${CONFIRM}`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}. MBR: storage ${fmt(mbrSt)} m, installation ${fmt(mbrIn)} m, operation ${fmt(mbrOp)} m (${fmt(mbrOp / b.OD)} × diameter).`,
          gov.uc > 1 ? (gov.name.includes('stiffener') ? `Lengthen the stiffener (try ${fmt(inp.stiffL * Math.min(2, Math.max(1.15, gov.uc)))} m) and re-check the curvature distribution, or reduce the hang-off angle range by re-orienting the I-tube.` : `Increase the governing radius to at least ${fmt(gov.demand * 1.05)} m or lower the tension in that stage.`)
            : 'Hold these radii as limits in the installation procedure and the stiffener specification; confirm the storage radius against the manufacturer-qualified reeling strain.',
          `${STD_UMB} §7.3; API RP 17B (2014) §4`),
        adv('economic', 'info', `The stiffener cone contains about ${fmt(puT)} t of polyurethane; the reel hub radius ${inp.reelR} m and chute radius ${inp.chuteR} m set the class of lay spread that can be used.`,
          `A storage MBR of ${fmt(mbrSt)} m allows a hub of that radius: if it fits standard reels the umbilical can be shipped and laid from a reel-drive spread instead of a carousel vessel, which usually shortens mobilisation by days.`,
          'Installation-method selection on vessel-day basis (AACE RP 18R-97 Class 4); DNV-ST-N001 marine operations planning',
          `Polyurethane embodied carbon ≈ ${fmt(puT * CO2_POLYMER)} t CO₂e (screening ${CO2_POLYMER} t CO₂e/t). An overbend during lay scraps the damaged length — avoiding it avoids re-manufacture of tubes with a far larger footprint.`),
      ],
      publish: { 'umbilical.mbrDynamic': mbrOp, 'umbilical.mbrStorage': mbrSt, 'umbilical.mbrInstall': mbrIn, 'umbilical.hangoffMoment': bs.rootMoment, 'umbilical.hangoffShear': bs.rootShear },
      quantities: { polymer_t: puT },
    };
  },
  verify: [
    { name: 'Known curvature limit of a single tube: MBR = E·(D/2)/(η·fy)',
      source: 'Euler–Bernoulli bending, σ = E·r/R. 200 GPa × 15.9 mm / (0.67 × 550 MPa) = 8.630 m; storage at η = 1: 5.782 m',
      inputs: { sheathT: 0, layers: [], Top: 1, Tinst: 0, zSection: 0, tipD: 0.032, rootD: 0.032, components: [{ name: 'Tube', kind: 'tube', layer: 0, count: 1, od: 31.8, wall: 2.6, E: 200, rho: 7800, fy: 550, p: 0, rhoF: 0 }] },
      // operating tension 1 kN adds σ = 1000/238.5e-6 = 4.19 MPa → MBR = 200e9 × 0.0159/(368.5e6 − 4.19e6) = 8.729 m
      expect: [{ kpi: 'mbrStorage', value: 5.782, tol: 0.003 }, { kpi: 'mbrInstall', value: 7.227, tol: 0.003 }, { kpi: 'mbrDynamic', value: 8.729, tol: 0.003 }] },
    { name: 'Tensioned beam leaving a clamp at a small angle: κ = θ·√(T/EI)',
      source: 'Sparks, Fundamentals of Marine Riser Mechanics (2007) Ch. 13. EI = 200 GPa × 25 622 mm⁴ = 5124 N·m²; T = 20 kN; θ = 1° → κ = 0.017453 × √(20000/5124) = 0.03448 1/m',
      inputs: { sheathT: 0, layers: [], Top: 20, Tinst: 0, zSection: 0, angle: 1, tipD: 0.032, rootD: 0.032, stiffL: 1, nSteps: 400, components: [{ name: 'Tube', kind: 'tube', layer: 0, count: 1, od: 31.8, wall: 2.6, E: 200, rho: 7800, fy: 550, p: 0, rhoF: 0 }] },
      expect: [{ kpi: 'stiffKmax', value: 0.03448, tol: 0.01 }] },
  ],
};
