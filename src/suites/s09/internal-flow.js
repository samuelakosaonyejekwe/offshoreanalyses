// Internal pressure, temperature and flow-induced loading of a riser: pressure containment, wall stresses with
// depth, effective versus true wall tension, thermal expansion, internal-flow momentum effects and slug loading.
import { wallStresses, catenaryRiser } from '../../engine/lines/riser.js';
import { burstPressureDnv, burstCheckDnv, burstPressureApi1111, hoopStress } from '../../engine/lines/risercheck.js';
import { bendForce } from '../../engine/flow/slugloads.js';
import { cyclesFromConstants } from '../../engine/fatigue/sn.js';
import {
  G, DEG, need, fmt, round, pipeInputs, waterInputs, pipeFrom, snInputs, snFrom, steelTonnes, adv, check, kpi, plot, finish, steelCo2Text, CONFIRM, CO2_STEEL,
} from './_common.js';
import { criteriaInputs } from './_scr.js';

export default {
  id: 's09-internal-flow', suite: 9, title: 'Internal pressure, temperature and flow loading',
  summary: 'Pressure containment, wall stress along depth, pressure and thermal axial effects, and loads from internal flow and slugging in a production riser.',
  covers: [12], fidelity: 'analytical',
  standards: ['DNV-ST-F201 (2018) Sec. 5 C (bursting)', 'API RP 1111 (2015) §4.3.1 burst design', 'API RP 2RD (1998, R2013) §5.2', 'DNV-ST-F101 (2021) Sec. 5 (temperature derating)', 'DNV-RP-C203 (2021) §2.4'],
  models: ['Lamé thick-wall stresses', 'Effective tension Te = Tw − pi·Ai + pe·Ae', 'Burst resistance (DNV and API forms)', 'Thermal expansion and restrained thermal force', 'Internal-flow momentum (ρAv²) reduction of effective tension', 'Centrifugal load of flow in a curved pipe', 'Slug weight variation and bend momentum force', 'S–N fatigue of slug cycles'],
  limits: 'Steady, single-density contents with a hydrostatic pressure profile; friction pressure drop, multiphase flow regime prediction, hydrate / wax and transient thermo-hydraulics belong to the flow-assurance suite. Slugs are idealised as liquid plugs of given length, speed and frequency; the riser is a planar catenary for the geometry terms. Temperature derating and design factors are overridable defaults.',
  publishes: ['riser.slugDamagePerYear', 'riser.burstUc', 'riser.thermalExpansion'],
  searchTerms: 'riser internal pressure burst slug loading thermal expansion effective tension',
  inputs: [
    ...pipeInputs(),
    ...waterInputs(),
    { key: 'tCorr', label: 'Corrosion allowance', unit: 'mm', type: 'number', default: 3, min: 0, max: 10, group: 'Riser pipe', role: 'geometry' },
    { key: 'tFab', label: 'Fabrication under-tolerance', unit: '%', type: 'number', default: 8, min: 0, max: 15, group: 'Riser pipe', role: 'geometry' },
    { key: 'pd', label: 'Design pressure at the top', unit: 'bar', type: 'number', default: 250, min: 0, max: 1500, group: 'Internal fluid', role: 'bc', uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'gInc', label: 'Incidental / design pressure ratio', unit: '-', type: 'number', default: 1.1, min: 1, max: 1.3, group: 'Internal fluid', role: 'criteria', help: `Typically 1.10 (DNV-ST-F201). ${CONFIRM}` },
    { key: 'Tin', label: 'Contents temperature', unit: '°C', type: 'number', default: 80, min: -20, max: 200, group: 'Internal fluid', role: 'bc', bind: 'pipeline.Tout' },
    { key: 'Tinst', label: 'Installation temperature', unit: '°C', type: 'number', default: 8, min: -5, max: 40, group: 'Internal fluid', role: 'ic', bind: 'site.Tsea' },
    { key: 'restraint', label: 'Axial restraint at the seabed end', unit: '-', type: 'number', default: 0.5, min: 0, max: 1, group: 'End conditions', role: 'bc', help: '0 = flowline free to expand into the riser, 1 = fully anchored. Scales the restrained thermal and pressure force.' },
    { key: 'Qm', label: 'Mass flow rate', unit: 'kg/s', type: 'number', default: 80, min: 0, max: 2000, group: 'Flow', role: 'bc' },
    { key: 'hangoffAngle', label: 'Hang-off angle from vertical', unit: 'deg', type: 'number', default: 12, min: 3, max: 30, group: 'End conditions', role: 'bc', bind: 'riser.hangoffAngleDeg' },
    { key: 'hangoffDepth', label: 'Hang-off depth below still water', unit: 'm', type: 'number', default: 15, min: 0, max: 200, group: 'End conditions', role: 'bc' },
    { key: 'slugRhoL', label: 'Slug liquid density', unit: 'kg/m³', type: 'number', default: 800, min: 300, max: 1300, group: 'Slugging', role: 'load' },
    { key: 'slugRhoG', label: 'Gas bubble density', unit: 'kg/m³', type: 'number', default: 150, min: 1, max: 500, group: 'Slugging', role: 'load' },
    { key: 'slugLength', label: 'Slug length', unit: 'm', type: 'number', default: 60, min: 0, max: 2000, group: 'Slugging', role: 'load', uncertain: { dist: 'lognormal', cov: 0.4 } },
    { key: 'slugV', label: 'Slug velocity', unit: 'm/s', type: 'number', default: 6, min: 0, max: 30, group: 'Slugging', role: 'load', uncertain: { dist: 'normal', cov: 0.2 } },
    { key: 'slugFreq', label: 'Slug frequency', unit: '1/h', type: 'number', default: 30, min: 0, max: 600, group: 'Slugging', role: 'load' },
    ...snInputs({ curve: 'E' }),
    ...criteriaInputs(),
    { key: 'gammaSC', label: 'Safety-class factor γSC (burst)', unit: '-', type: 'number', default: 1.14, min: 1, max: 1.4, group: 'Acceptance criteria', role: 'criteria', help: `1.04 / 1.14 / 1.26 low / normal / high. ${CONFIRM}` },
    { key: 'gammaM', label: 'Material factor γm', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.3, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'fdApi', label: 'API RP 1111 riser design factor fd', unit: '-', type: 'number', default: 0.75, min: 0.5, max: 0.9, group: 'Acceptance criteria', role: 'criteria', help: `0.75 for risers, 0.90 for pipelines; design pressure ≤ 0.80·fd·Pb. ${CONFIRM}` },
    { key: 'fyDerate', label: 'Yield derating at temperature', unit: 'MPa', type: 'number', default: 18, min: 0, max: 150, group: 'Acceptance criteria', role: 'criteria', help: `C-Mn steel: none below 50 °C, about 0.6 MPa/°C above (18 MPa at 80 °C, 30 MPa at 100 °C; DNV-ST-F101 Fig. 5-2). ${CONFIRM}` },
  ],

  compute(inp, ctx) {
    // engineering units → SI: mm → m, bar → Pa, MPa → Pa
    const pp = pipeFrom(inp, ctx), { D, t, mat, props } = pp;
    const t1 = t * (1 - inp.tFab / 100) - inp.tCorr / 1e3, t2 = t - inp.tCorr / 1e3;
    need(t < D / 2.2, `Wall thickness ${inp.t} mm is too large for a ${inp.D} mm pipe. Reduce the wall or increase the diameter.`);
    need(t1 > 0.002, `Wall after tolerance and corrosion allowance is ${fmt(t1 * 1e3)} mm: too thin. Increase the wall.`);
    const fy = Math.max(mat.fy - inp.fyDerate * 1e6, 0.5 * mat.fy), fu = Math.max(mat.fu - inp.fyDerate * 1e6, 0.5 * mat.fu);
    const h = inp.depth - inp.hangoffDepth; need(h > 20, 'The hang-off must be at least 20 m above the seabed.');
    need(pp.section.w > 0, 'The riser is buoyant with these contents; the catenary geometry terms need positive apparent weight.');
    const cat = catenaryRiser({ h, w: pp.section.w, hangoffAngleDeg: inp.hangoffAngle });
    const pTop = inp.pd * 1e5, rhoI = inp.rhoInt, rhoW = pp.rhoW;
    const piAt = (z) => pTop + rhoI * G * z, peAt = (z) => rhoW * G * (inp.hangoffDepth + z);   // z = depth below the hang-off
    // profile with depth: effective tension of the catenary Te = H + w·(height above touchdown)
    const zs = [], piA = [], peA = [], hoop = [], vm = [], ax = [], ucB = [];
    let vmMax = 0, zVm = 0, bMax = null, zB = 0;
    for (let i = 0; i <= 40; i++) {
      const z = h * i / 40, pi = piAt(z), pe = peAt(z), Te = cat.H + pp.section.w * (h - z);
      const ws = wallStresses({ Te, M: 0, pi, pe, D, t: t2 });
      const bc = burstCheckDnv({ pli: inp.gInc * pTop + rhoI * G * z, pe, D, t: t1, fy, fu, gammaSC: inp.gammaSC, gammaM: inp.gammaM });
      zs.push(z + inp.hangoffDepth); piA.push(pi / 1e6); peA.push(pe / 1e6); hoop.push(ws.hoopInner / 1e6); vm.push(ws.vonMises / 1e6); ax.push(ws.axial / 1e6); ucB.push(bc.uc);
      if (ws.vonMises > vmMax) { vmMax = ws.vonMises; zVm = z; }
      if (!bMax || bc.uc > bMax.uc) { bMax = bc; zB = z; }
    }
    const top = wallStresses({ Te: cat.topTension, M: 0, pi: pTop, pe: peAt(0), D, t: t2 }), bot = wallStresses({ Te: cat.H, M: 0, pi: piAt(h), pe: peAt(h), D, t: t2 });
    const lameTop = wallStresses({ Te: 0, M: 0, pi: pTop, pe: 0, D, t }).hoopInner;
    const pb = burstPressureDnv({ D, t: t1, fy, fu }), pbApi = burstPressureApi1111({ D, t: t2, S: fy, U: fu });
    const netTop = pTop - peAt(0), netBot = piAt(h) - peAt(h), netMax = Math.max(netTop, netBot);
    const barlow = hoopStress({ pi: piAt(zB), pe: peAt(zB), D, t: t1 });
    const allow = inp.Cf * inp.CaAllow * fy;
    // thermal and pressure axial effects
    const dT = inp.Tin - inp.Tinst, S = cat.suspendedLength;
    const expTh = mat.alpha * dT * S;
    const pMid = 0.5 * (pTop + piAt(h)), peMid = 0.5 * (peAt(0) + peAt(h));
    const expP = (1 - 2 * mat.nu) * (pMid * props.Ai - peMid * props.Ae) / props.EA * S;
    const Frestr = inp.restraint * (mat.E * props.As * mat.alpha * dT + (1 - 2 * mat.nu) * (piAt(h) * props.Ai - peAt(h) * props.Ae));
    // flow
    const v = rhoI > 0 ? inp.Qm / (rhoI * props.Ai) : 0, mom = rhoI * props.Ai * v * v;
    const vCrit = rhoI > 0 ? Math.sqrt(cat.H / (rhoI * props.Ai)) : 1e6;
    const qCent = mom * cat.tdpCurvature;
    // slugging
    const Di = D - 2 * t, dW = (inp.slugRhoL - inp.slugRhoG) * props.Ai * G;
    const slugT = dW * Math.min(inp.slugLength, S) * Math.cos(0.5 * inp.hangoffAngle * DEG);    // weight carried mainly on the steep upper part
    const bend = bendForce({ rho: inp.slugRhoL, v: inp.slugV, D: Di, bendAngle: (90 - inp.hangoffAngle) * DEG });
    const momSlug = inp.slugRhoL * props.Ai * inp.slugV ** 2;
    const dSig = (slugT + momSlug) / props.As;
    const { sn, label } = snFrom(inp, t);
    const nYr = inp.slugFreq * 8766, dSlug = dSig > 0 && nYr > 0 ? nYr / cyclesFromConstants(sn, dSig) : 0;
    const dKap = inp.slugLength > 0 ? dW / cat.H : 0;   // change of sag-bend curvature w/H while a slug sits in the bend
    const dSigB = mat.E * dKap * D / 2;

    const kpis = [
      kpi('piBottom', 'Internal pressure at the seabed', piAt(h) / 1e6, 'MPa', 'Design pressure at the top plus the head of the contents. Dense contents make the bottom of the riser the highest-pressure point.'),
      kpi('netMax', 'Largest net internal overpressure', netMax / 1e6, 'MPa', `Internal minus external pressure where it is greatest (${netTop >= netBot ? 'at the top: the contents are lighter than sea water, so the surface end governs burst' : 'at the seabed: the contents are heavier than sea water'}).`),
      kpi('hoopTop', 'Lamé hoop stress at the bore, top (nominal wall, no external pressure)', lameTop / 1e6, 'MPa', 'Thick-wall hoop stress at the inner surface from design pressure alone on the nominal wall: the reference value for the verification case.'),
      kpi('hoopMax', 'Maximum hoop stress (corroded wall)', Math.max(...hoop), 'MPa', 'Largest hoop stress at the bore along the riser on the corroded wall, including external pressure.'),
      kpi('pb', 'Burst resistance pb (DNV)', pb / 1e6, 'MPa', 'Pressure-containment resistance of the minimum wall, (2/√3)·(2t/(D − t))·min(fy, fu/1.15), with the temperature-derated strength.'),
      kpi('burstUc', 'Burst utilisation (DNV)', Math.max(0, bMax.uc), '-', `Local incidental net pressure over the factored burst resistance, at its worst (${fmt(zB + inp.hangoffDepth)} m water depth). Above 1.0 the wall is too thin for pressure containment.`, 1, 'max'),
      kpi('burstApiUc', 'Burst utilisation (API RP 1111)', netMax / (0.8 * inp.fdApi * pbApi), '-', 'Design net pressure over 0.80·fd·Pb with the API burst pressure 0.45(S + U)·ln(D/Di). An independent check on the wall for pressure.', 1, 'max'),
      kpi('vmMax', 'Maximum von Mises stress (pressure + tension)', vmMax / 1e6, 'MPa', `Combined stress from pressure and the catenary effective tension without bending, greatest at ${fmt(zVm + inp.hangoffDepth)} m water depth, against the API RP 2RD allowable with derated yield.`, allow / 1e6, 'max'),
      kpi('TwTop', 'True wall tension at the top', top.Tw / 1e3, 'kN', 'Axial force in the steel at the hang-off: effective tension plus pi·Ai − pe·Ae. Internal pressure raises it well above the effective tension.'),
      kpi('TwBottom', 'True wall tension at the touchdown', bot.Tw / 1e3, 'kN', 'Axial force in the steel at the seabed. It can differ in sign from the effective tension; only effective tension governs buckling.'),
      kpi('expansion', 'Free thermal expansion of the suspended riser', expTh, 'm', 'Lengthening of the suspended pipe from installation to operating temperature if unrestrained. It lowers the sag bend and moves the touchdown point.'),
      kpi('expansionP', 'Pressure-induced elongation', expP, 'm', 'Change of length from the pressure end-cap and Poisson effects, (1 − 2ν)(pi·Ai − pe·Ae)/EA over the suspended length.'),
      kpi('restrainedForce', 'Restrained expansion force at the seabed end', Frestr / 1e3, 'kN', 'Compressive effective force generated where thermal and pressure expansion are restrained (scaled by the restraint factor). It pushes the flowline into the touchdown zone and reduces bottom tension.'),
      kpi('velocity', 'Mean flow velocity', v, 'm/s', 'Mass flow over density and bore area. It drives erosion, momentum loads at bends and the flow-induced reduction of effective tension.'),
      kpi('momentum', 'Internal-flow momentum flux ρAv²', mom / 1e3, 'kN', 'Momentum flux of the contents. It acts like an extra compression in the effective-tension equation for a curved pipe (Te,flow = Te − ρAv²).'),
      kpi('vCrit', 'Flow velocity that cancels the bottom tension', vCrit, 'm/s', 'Velocity at which ρAv² equals the effective tension at touchdown; flow-induced buckling of the sag bend would start near it. Normal production velocities are far below.'),
      kpi('centrifugal', 'Centrifugal load in the sag bend', qCent, 'N/m', 'Distributed outward load ρAv²·κ on the pipe where the flow follows the sag-bend curvature, to compare with the submerged weight.'),
      kpi('slugTension', 'Top-tension variation per slug', slugT / 1e3, 'kN', 'Change of top tension when a liquid slug replaces gas over the slug length: the cyclic axial load from slugging.'),
      kpi('slugBendForce', 'Slug momentum force at the sag bend', bend.resultant / 1e3, 'kN', 'Resultant force on the sag bend when a slug front turns through it, ρAv²·√(2(1 − cos θ)). It excites the riser at the slug frequency.'),
      kpi('slugStress', 'Stress range per slug (axial)', dSig / 1e6, 'MPa', 'Nominal axial stress range in the wall from slug weight and momentum variation; the bending range at the sag bend is reported in the table.'),
      kpi('slugDamage', 'Annual fatigue damage from slugging', dSlug, '1/yr', 'Miner damage per year from slug cycles on the selected S–N curve. Add it to wave and VIV damage in the fatigue-life module.'),
    ];
    const checks = [
      check('Pressure containment (burst), DNV', Math.max(0, bMax.demand) / 1e6, bMax.capacity / 1e6, 'MPa', 'DNV-ST-F201 (2018) Sec. 5 C200 (Eq. 5.14–5.15)', 'Guards against rupture of the wall under the local incidental pressure. Failure means loss of containment.'),
      check('Burst design pressure, API RP 1111', netMax / 1e6, 0.8 * inp.fdApi * pbApi / 1e6, 'MPa', 'API RP 1111 (2015) §4.3.1: Pd ≤ 0.80·fd·fe·ft·Pb', 'Independent pressure-design criterion with the riser design factor; it limits design pressure to a fraction of the specified burst pressure.'),
      check('Combined stress, pressure + tension (von Mises)', vmMax / 1e6, allow / 1e6, 'MPa', 'API RP 2RD (1998) §5.2.3 with temperature-derated yield', 'Guards against yielding under hoop, radial and axial stress together, before any bending is added by the global analysis.'),
      check('Flow momentum vs bottom tension', mom / 1e3, cat.H / 1e3, 'kN', 'Internal-flow effective tension (Païdoussis; Sparks 2007 Ch. 3)', 'Compares ρAv² with the effective tension at touchdown. If the momentum flux approached the tension, the sag bend would lose its restoring stiffness and could buckle.'),
    ];

    const tt = [], Tt = []; const per = inp.slugFreq > 0 ? 3600 / inp.slugFreq : 120, pass = inp.slugV > 0 ? Math.min(per, S / inp.slugV) : per;
    for (let c = 0; c < 4; c++) { tt.push(c * per, c * per + 0.15 * pass, c * per + pass, c * per + 1.15 * pass); Tt.push(cat.topTension / 1e3, (cat.topTension + slugT) / 1e3, (cat.topTension + slugT) / 1e3, cat.topTension / 1e3); }
    const qs = [], tf = []; for (let i = 0; i <= 30; i++) { const q = Math.max(inp.Qm, 50) * 6 * i / 30; qs.push(q); tf.push((cat.H - (rhoI > 0 ? q * q / (rhoI * props.Ai) : 0)) / 1e3); }
    const plots = [
      plot('Pressure with depth', 'Pressure', 'MPa', 'Water depth', 'm', [{ name: 'Internal', x: piA, y: zs }, { name: 'External', x: peA, y: zs }, { name: 'Net internal', x: piA.map((p, i) => p - peA[i]), y: zs }], 'Internal, external and net pressure down the riser. The slope difference shows whether burst governs at the top (light contents) or at the bottom (heavy contents).'),
      plot('Wall stress with depth', 'Stress', 'MPa', 'Water depth', 'm', [{ name: 'Hoop at bore', x: hoop, y: zs }, { name: 'Axial (true wall)', x: ax, y: zs }, { name: 'von Mises', x: vm, y: zs }, { name: 'Allowable', x: [allow / 1e6, allow / 1e6], y: [zs[0], zs[40]], dash: true }], 'Hoop, axial and von Mises stress in the wall against depth, without bending. The gap to the dashed allowable is what remains for bending in the global analysis.'),
      plot('Burst utilisation with depth', 'Utilisation', '-', 'Water depth', 'm', [{ name: 'DNV burst', x: ucB, y: zs }, { name: 'Limit', x: [1, 1], y: [zs[0], zs[40]], dash: true }], 'Pressure-containment utilisation along the riser; the governing depth is where the curve is furthest right.'),
      plot('Top tension during slugging (operating transient)', 'Time', 's', 'Top tension', 'kN', [{ name: 'Top tension', x: tt, y: Tt }], 'Idealised top-tension history for four slugs: tension rises as each liquid slug fills the riser and falls as it leaves. The step height is the cyclic load; the spacing is the slug period.'),
      plot('Bottom tension against flow rate', 'Mass flow rate', 'kg/s', 'Effective tension at touchdown incl. flow', 'kN', [{ name: 'Te − ρAv²', x: qs, y: tf }, { name: 'Operating point', x: [inp.Qm], y: [(cat.H - mom) / 1e3] }], 'Effective tension at the touchdown after subtracting the flow momentum flux. The curve reaching zero marks the flow rate at which flow-induced buckling would begin.'),
    ];
    const tables = [{ title: 'Pressure, temperature and flow effects', columns: [{ key: 'item', label: 'Quantity' }, { key: 'value', label: 'Value' }, { key: 'unit', label: 'Unit' }],
      rows: [['Design pressure at top', pTop / 1e6, 'MPa'], ['Internal pressure at seabed', piAt(h) / 1e6, 'MPa'], ['External pressure at seabed', peAt(h) / 1e6, 'MPa'], ['Hoop stress, mean-diameter form at governing depth', barlow.meanDiameter / 1e6, 'MPa'], ['Burst pressure DNV', pb / 1e6, 'MPa'], ['Burst pressure API', pbApi / 1e6, 'MPa'], ['Effective tension at top (catenary)', cat.topTension / 1e3, 'kN'], ['Effective tension at touchdown', cat.H / 1e3, 'kN'], ['Temperature rise', dT, '°C'], ['Derated yield', fy / 1e6, 'MPa'], ['Free thermal strain', mat.alpha * dT * 100, '%'], ['Slug weight change', dW, 'N/m'], ['Slug cycles per year', nYr, '1/yr'], ['Sag-bend curvature change per slug', dKap * 1e3, '1/km'], ['Sag-bend bending stress range per slug', dSigB / 1e6, 'MPa']].map(([item, value, unit]) => ({ item, value: round(value, 5), unit })),
      explain: 'Intermediate values behind the pressure, thermal and slug results, for the calculation note and for hand checking.' }];
    const ucMax = Math.max(...checks.map((c) => c.uc)), gov = checks.find((c) => c.uc === ucMax);
    // nominal wall for burst uc = 1: pb ∝ 2t1/(D − t1)
    const pReq = bMax.demand * inp.gammaSC * inp.gammaM, k = pReq / ((2 / Math.sqrt(3)) * Math.min(fy, fu / 1.15)), t1Req = k * D / (2 + k), tNom = (t1Req + inp.tCorr / 1e3) / (1 - inp.tFab / 100);
    const perMm = steelTonnes(D, 0.001, 1000, mat.rho), dWall = (tNom - t) * 1e3;
    const advice = [
      adv('engineering', ucMax > 1 ? 'fail' : ucMax > 0.85 ? 'warn' : 'ok',
        `${gov.name} governs at ${fmt(ucMax)}. Net pressure ${fmt(netMax / 1e6)} MPa; burst needs a nominal wall of ${fmt(tNom * 1e3)} mm (selected ${inp.t} mm). Slugging cycles the top tension by ${fmt(slugT / 1e3)} kN ${fmt(nYr)} times a year (damage ${fmt(dSlug)} per year).`,
        ucMax > 1 ? (gov.name.startsWith('Flow') ? `Reduce the flow velocity below ${fmt(0.5 * vCrit)} m/s or raise the bottom tension (larger hang-off angle).` : `Increase the nominal wall to at least ${fmt(tNom * 1e3)} mm or use a higher grade; check the corrosion allowance of ${inp.tCorr} mm against the corrosion assessment, because it is paid for in wall everywhere.`)
          : `Pass the true wall tension and the slug stress range to the local strength and fatigue modules; ${dSigB > dSig ? `the sag-bend bending range per slug (${fmt(dSigB / 1e6)} MPa) exceeds the axial range and should be included in touchdown fatigue.` : 'the axial slug range governs over the sag-bend bending range.'}`,
        `DNV-ST-F201 (2018) Sec. 5 C; API RP 1111 (2015) §4.3.1; temperature derating per DNV-ST-F101 (2021) Sec. 5. ${CONFIRM}`),
      adv('economic', dWall > 0 ? 'warn' : 'info',
        dWall > 0 ? `Pressure containment needs ${fmt(dWall)} mm more wall: ${fmt(perMm * dWall)} t of steel per kilometre.` : `The selected wall exceeds the burst minimum by ${fmt(-dWall)} mm (${fmt(perMm * -dWall)} t of steel per kilometre). Of the nominal wall, ${fmt(inp.tCorr + inp.t * inp.tFab / 100)} mm is corrosion allowance and tolerance that carry no design pressure.`,
        'Weigh the corrosion allowance against inhibition or a corrosion-resistant liner on whole-life terms, and check whether slug mitigation (topsides choking, gas lift) is cheaper than designing the riser welds for the slug cycles.',
        'Whole-life costing (ISO 15686-5): pipe tonnage and lay time against chemical inhibition and production deferment.',
        `Each millimetre of wall is ${fmt(perMm)} t of steel per kilometre ≈ ${fmt(perMm * CO2_STEEL)} t CO₂e (screening factor). ${dWall > 0 ? steelCo2Text(perMm * dWall) : 'A thinner, higher-grade or lined pipe reduces both steel mass and hull payload.'}`),
    ];
    return finish({
      kpis, checks, plots, tables,
      notes: [`Contents ${rhoI} kg/m³ at ${inp.pd} bar and ${inp.Tin} °C; yield derated by ${inp.fyDerate} MPa. Burst is assessed on t1 = ${fmt(t1 * 1e3)} mm, stresses on t2 = ${fmt(t2 * 1e3)} mm.`, `${gov.name} governs (${fmt(ucMax)}). Slug fatigue uses ${label}.`],
      advice,
      publish: { 'riser.slugDamagePerYear': dSlug, 'riser.burstUc': Math.max(0, bMax.uc), 'riser.thermalExpansion': expTh },
      quantities: { steel_t: Math.max(0, perMm * dWall) },
    });
  },

  verify: [
    {
      // Lamé thick cylinder under internal pressure only: σθ(ri) = p·(ro² + ri²)/(ro² − ri²).
      // p = 25 MPa, ro = 136.55 mm, ri = 115.95 mm → σθ = 25·(18646 + 13444)/(18646 − 13444) = 154.2 MPa.
      name: 'Lamé hoop stress at the bore σθ = p·(ro² + ri²)/(ro² − ri²)', source: 'Timoshenko & Goodier, Theory of Elasticity §28',
      inputs: {},
      expect: [{ kpi: 'hoopTop', value: 25 * (0.13655 ** 2 + 0.11595 ** 2) / (0.13655 ** 2 - 0.11595 ** 2), tol: 1e-4 }],
    },
    {
      // DNV burst resistance with no tolerance, corrosion or derating: pb = (2/√3)·(2t/(D − t))·min(fy, fu/1.15)
      // = 1.1547·(2·20.6/252.5)·min(450, 535/1.15 = 465.2) = 84.78 MPa.
      name: 'Burst resistance pb = (2/√3)·(2t/(D − t))·min(fy, fu/1.15)', source: 'DNV-ST-F201 (2018) Eq. 5.15; DNV-ST-F101 (2021) Eq. 5.8',
      inputs: { tCorr: 0, tFab: 0, fyDerate: 0 },
      expect: [{ kpi: 'pb', value: 2 / Math.sqrt(3) * (2 * 20.6 / (273.1 - 20.6)) * Math.min(450, 535 / 1.15), tol: 1e-4 }],
    },
  ],
};
