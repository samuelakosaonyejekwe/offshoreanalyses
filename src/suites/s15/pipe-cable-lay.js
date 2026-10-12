// Pipeline, cable and umbilical laying: catenary tensions and lay-back, sagbend and overbend strain, bend
// radius and sidewall pressure on the chute, reeling strain, and quasi-dynamic tension bounds.
import { stiffenedCatenary, cableLay, reelLay, layDynamicTension, layCriteriaF101 } from '../../engine/lines/pipelay.js';
import { G, DEG, need, run, round, fmt, kpi, check, adv, line, tbl, worst, linspace, steelCo2Text, vesselCo2Text, CONFIRM } from './_common.js';

const METHODS = [{ value: 'slay', label: 'S-lay (rigid pipe over a stinger)' }, { value: 'jlay', label: 'J-lay (rigid pipe from a tower)' }, { value: 'reel', label: 'Reel-lay (rigid pipe, reeled)' }, { value: 'cable', label: 'Cable / umbilical / flexible over a chute' }];
const rigid = ['slay', 'jlay', 'reel'];

export default {
  id: 's15-pipe-cable-lay', suite: 15, title: 'Pipeline, cable and umbilical laying',
  summary: 'Gives the lay tension, lay-back and suspended length, the sagbend and overbend strain of a rigid pipe or the bend radius and sidewall pressure of a cable, the reeling strain, and the dynamic tension range for the lay vessel motion.',
  covers: [8], fidelity: 'analytical',
  standards: ['DNV-ST-F101 (2021) Sec. 10 (installation: simplified laying criteria) and Sec. 5 (ovalisation)', 'DNV-RP-0360 (2016) §4 (subsea power cables: installation, minimum bend radius, sidewall pressure)', 'API RP 17B (2014) §9 / ISO 13628-5 (flexible pipe and umbilical installation limits)', 'Dixon & Rutledge (1968) Stiffened catenary calculations in pipeline laying, J. Eng. Ind. 90'],
  models: ['Nonlinear catenary equations', 'Stiffened catenary (boundary-layer bending correction)', 'Overbend and sagbend analysis', 'Tensioner stiffness / quasi-dynamic tension bounds', 'Seabed touchdown (catenary touchdown radius)', 'Minimum bend radius and sidewall pressure', 'Elastic–perfectly-plastic reel bending and spring-back'],
  limits: 'Static in-plane catenary with a first-order bending correction: no current, no seabed slope, no roller-by-roller stinger contact and no pipe rotation or residual curvature tracking. Dynamic tension is bounded between a quasi-static and a drag-locked estimate, not simulated; a rod finite-element lay analysis with vessel RAOs is needed for the installation procedure. Strain limits of DNV-ST-F101 are recalled defaults to be confirmed.',
  publishes: ['install.layTension', 'install.layBack'],
  searchTerms: 'pipelay S-lay J-lay reel-lay cable lay catenary sagbend overbend minimum bend radius',
  inputs: [
    { key: 'method', label: 'Lay method', type: 'select', default: 'slay', options: METHODS, group: 'Method', role: 'operation' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 120, min: 5, max: 3500, group: 'Site', role: 'environment', bind: 'site.depth' },
    { key: 'D', label: 'Steel outer diameter', unit: 'mm', type: 'number', default: 508, min: 50, max: 1500, group: 'Pipe', role: 'geometry', bind: 'pipeline.D', showIf: { key: 'method', in: rigid } },
    { key: 't', label: 'Wall thickness', unit: 'mm', type: 'number', default: 15.9, min: 3, max: 60, group: 'Pipe', role: 'geometry', bind: 'pipeline.t', design: true, showIf: { key: 'method', in: rigid } },
    { key: 'mat', label: 'Line-pipe grade', type: 'select', source: 'materials:steel', default: 'api5l-x65', group: 'Pipe', role: 'material', showIf: { key: 'method', in: rigid } },
    { key: 'tCoat', label: 'Concrete / coating thickness', unit: 'mm', type: 'number', default: 50, min: 0, max: 200, group: 'Pipe', role: 'geometry', showIf: { key: 'method', in: rigid } },
    { key: 'rhoCoat', label: 'Coating density', unit: 'kg/m³', type: 'number', default: 3040, min: 500, max: 3600, group: 'Pipe', role: 'material', showIf: { key: 'method', in: rigid } },
    { key: 'wCable', label: 'Cable submerged weight', unit: 'N/m', type: 'number', default: 450, min: 5, max: 5000, group: 'Cable', role: 'load', showIf: { key: 'method', in: ['cable'] } },
    { key: 'mCable', label: 'Cable mass in air', unit: 'kg/m', type: 'number', default: 75, min: 1, max: 600, group: 'Cable', role: 'load', showIf: { key: 'method', in: ['cable'] } },
    { key: 'eaCable', label: 'Cable axial stiffness', unit: 'MN', type: 'number', default: 600, min: 1, max: 5000, group: 'Cable', role: 'material', showIf: { key: 'method', in: ['cable'] } },
    { key: 'mbr', label: 'Minimum bend radius of the product', unit: 'm', type: 'number', default: 3.0, min: 0.2, max: 30, group: 'Cable', role: 'criteria', showIf: { key: 'method', in: ['cable'] }, help: 'Installation (dynamic) MBR from the manufacturer, under tension.' },
    { key: 'swpLimit', label: 'Allowable sidewall pressure', unit: 'kN/m', type: 'number', default: 30, min: 1, max: 500, group: 'Cable', role: 'criteria', showIf: { key: 'method', in: ['cable'] } },
    { key: 'tAllow', label: 'Allowable installation tension of the product', unit: 'kN', type: 'number', default: 250, min: 1, max: 20000, group: 'Cable', role: 'criteria', showIf: { key: 'method', in: ['cable'] } },
    { key: 'chuteR', label: 'Chute / lay-wheel radius', unit: 'm', type: 'number', default: 4.5, min: 0.3, max: 30, group: 'Lay equipment', role: 'geometry', design: true, showIf: { key: 'method', in: ['cable'] } },
    { key: 'Ttop', label: 'Lay tension at the tensioner', unit: 'kN', type: 'number', default: 900, min: 1, max: 30000, group: 'Lay equipment', role: 'bc', design: true, help: 'Boundary condition of the catenary: the tension held by the tensioners (or the cable engine).' },
    { key: 'Tcap', label: 'Tensioner capacity', unit: 'kN', type: 'number', default: 1500, min: 1, max: 40000, group: 'Lay equipment', role: 'bc' },
    { key: 'hTop', label: 'Height of the departure point above still water', unit: 'm', type: 'number', default: -18, min: -200, max: 80, group: 'Lay equipment', role: 'bc', help: 'S-lay: elevation of the stinger tip (negative below water). J-lay, reel and cable: elevation of the tower exit or chute.' },
    { key: 'Rst', label: 'Stinger radius', unit: 'm', type: 'number', default: 180, min: 30, max: 1000, group: 'Lay equipment', role: 'geometry', design: true, showIf: { key: 'method', in: ['slay'] } },
    { key: 'reelR', label: 'Reel hub radius', unit: 'm', type: 'number', default: 8.2, min: 2, max: 30, group: 'Lay equipment', role: 'geometry', showIf: { key: 'method', in: ['reel'] } },
    { key: 'reelStrainMax', label: 'Allowable reeling strain', unit: '%', type: 'number', default: 2.5, min: 0.3, max: 5, group: 'Criteria', role: 'criteria', showIf: { key: 'method', in: ['reel'] }, help: `Qualified by the engineering critical assessment of the girth welds; commonly 2–3 %. ${CONFIRM}` },
    { key: 'grade', label: 'Overbend strain criterion grade', type: 'select', default: 'X65', options: ['X52', 'X60', 'X65', 'X70'].map((g) => ({ value: g, label: g })), group: 'Criteria', role: 'criteria', showIf: { key: 'method', in: ['slay'] }, help: `Simplified laying criteria, static (Criterion I): X65 0.25 %, X70 0.27 %. ${CONFIRM}` },
    { key: 'heave', label: 'Vertical motion amplitude of the departure point', unit: 'm', type: 'number', default: 0.8, min: 0, max: 8, group: 'Vessel motion', role: 'bc', bind: 'install.craneTipHeave', help: 'Single amplitude in the limiting lay sea state.' },
    { key: 'Tmot', label: 'Motion period', unit: 's', type: 'number', default: 8, min: 3, max: 25, group: 'Vessel motion', role: 'bc' },
    { key: 'layRate', label: 'Lay rate', unit: 'km/day', type: 'number', default: 3, min: 0.05, max: 30, group: 'Operation', role: 'operation' },
    { key: 'length', label: 'Length to lay', unit: 'km', type: 'number', default: 40, min: 0.1, max: 2000, group: 'Operation', role: 'operation', bind: 'pipeline.length' },
  ],
  compute(inp, ctx) {
    // mm → m, kN → N, MN → N, km kept.
    const rho = 1025, h = inp.depth + inp.hTop;
    need(h > 1, 'The departure point is at or below the seabed.');
    const T = inp.Ttop * 1e3, isCable = inp.method === 'cable';
    const kpis = [], checks = [], rows = [];
    let w, EA, mLin, H, layBack, susp, Tbot, angDeg, steelPerKm = 0, sigma = [], curv = [], tReqNote = '';
    if (isCable) {
      w = inp.wCable; EA = inp.eaCable * 1e6; mLin = inp.mCable;
      need(T > w * h, `The lay tension of ${inp.Ttop} kN is less than the weight of the suspended cable (${fmt(w * h / 1e3)} kN): the cable cannot be held. Increase the tension.`);
      const c = run(() => cableLay({ w, h, topTension: T, chuteRadius: inp.chuteR, mbr: inp.mbr, swpLimit: inp.swpLimit * 1e3 }), 'Cable lay');
      H = c.H; layBack = c.layBack; susp = c.suspendedLength; Tbot = H; angDeg = c.departureAngleDeg;
      checks.push(
        check('Chute radius against MBR', inp.mbr, inp.chuteR, 'm (required / provided)', 'DNV-RP-0360 (2016) §4 (minimum bend radius during installation)', 'Bending the cable tighter than its MBR under tension damages armour, insulation or fibre elements.'),
        check('Touchdown bend radius against MBR', inp.mbr, c.tdpBendRadius, 'm (required / achieved)', 'DNV-RP-0360 (2016) §4; catenary touchdown radius a = H/w', 'Too little bottom tension lets the cable kink or loop at the touchdown point.'),
        check('Sidewall pressure on the chute', c.sidewallPressure / 1e3, inp.swpLimit, 'kN/m', 'DNV-RP-0360 (2016) §4 (SWP = T/R)', 'Radial crushing load where the tensioned cable bears on the chute; exceedance deforms the cable cross-section.'),
        check('Top tension against the allowable installation tension', T / 1e3, inp.tAllow, 'kN', 'Manufacturer installation limits; DNV-RP-0360 §4', 'Tension the product can take while being bent over the chute.'));
      kpis.push(kpi('tdpR', 'Touchdown bend radius', c.tdpBendRadius, 'm', 'Radius of curvature at the seabed touchdown, a = H/w. It must stay above the product MBR.', inp.mbr, 'min'),
        kpi('swp', 'Sidewall pressure on the chute', c.sidewallPressure / 1e3, 'kN/m', 'Top tension divided by chute radius.', inp.swpLimit, 'max'));
      for (const s of linspace(0, susp, 40)) { sigma.push(s); curv.push(1 / ((H / w) * (1 + (s * w / H) ** 2))); }
    } else {
      const steel = ctx.material(inp.mat); need(steel, `Steel '${inp.mat}' is not in the materials library.`);
      const D = inp.D / 1e3, t = inp.t / 1e3, E = steel.E || 2.07e11, fy = steel.fy;
      need(t < D / 2, 'Wall thickness must be smaller than the radius.');
      const As = Math.PI / 4 * (D * D - (D - 2 * t) ** 2), I = Math.PI / 64 * (D ** 4 - (D - 2 * t) ** 4);
      const Dc = D + 2 * inp.tCoat / 1e3, Ac = Math.PI / 4 * (Dc * Dc - D * D);
      mLin = As * 7850 + Ac * inp.rhoCoat; EA = E * As;
      w = (mLin - rho * Math.PI / 4 * Dc * Dc) * G; steelPerKm = As * 7850;
      need(w > 0, `The empty pipe floats (submerged weight ${fmt(w)} N/m): add weight coating or wall thickness before it can be laid.`);
      // S-lay: tensioner tension less the weight of pipe between the tensioner level (taken at the surface) and the stinger tip.
      const Tdep = inp.method === 'slay' ? T - w * Math.max(-inp.hTop, 0) : T;
      H = Tdep - w * h;
      need(H > 0, `The lay tension of ${inp.Ttop} kN does not carry the suspended pipe (weight ${fmt(w * h / 1e3)} kN to the departure point): no catenary exists. Increase the tension.`);
      const sc = run(() => stiffenedCatenary({ EI: E * I, H, w, h, D, EA, topFree: inp.method !== 'jlay' }), 'Stiffened catenary');
      susp = sc.suspendedLength; Tbot = H; angDeg = Math.acos(Math.min(H / Tdep, 1)) / DEG; layBack = (H / w) * Math.asinh(susp * w / H);
      sigma = Array.from(sc.sigma); curv = Array.from(sc.curvature);
      const overStrain = inp.method === 'slay' ? D / (2 * inp.Rst) + T / EA : 0;
      const pe = rho * G * inp.depth;
      const crit = run(() => layCriteriaF101({ grade: inp.grade, criterion: 1, overbendStrain: Math.max(overStrain, 1e-9), D, t, fy, sagbendTe: H, sagbendMoment: sc.sagbendMoment, pe }), 'Laying criteria');
      checks.push(check('Sagbend equivalent stress', crit.sagbendStress / 1e6, crit.sagbendAllowable / 1e6, 'MPa', 'DNV-ST-F101 (2021) Sec. 10 (simplified laying criteria: 0.87·fy in the sagbend)', 'Bending at the sagbend combined with tension and external pressure. Too little lay tension tightens the sagbend until the pipe yields or buckles.'));
      if (inp.method === 'slay') checks.push(check('Overbend strain on the stinger', overStrain * 100, crit.overbendLimit * 100, '%', 'DNV-ST-F101 (2021) Sec. 10 (simplified laying criteria, Criterion I by grade)', 'Strain from bending over the stinger plus tension. Exceedance damages the concrete coating and accumulates plastic strain in the pipe.'));
      kpis.push(kpi('sagStrain', 'Sagbend bending strain', sc.sagbendStrain * 100, '%', 'Bending plus axial strain at the point of largest sagbend curvature from the stiffened catenary.'),
        kpi('sagStress', 'Sagbend equivalent stress', crit.sagbendStress / 1e6, 'MPa', 'Von Mises stress from tension, bending and external pressure at the sagbend.', crit.sagbendAllowable / 1e6, 'max'),
        kpi('sagMoment', 'Sagbend bending moment', sc.sagbendMoment / 1e3, 'kN·m', 'EI × maximum sagbend curvature.'),
        kpi('overStrain', 'Overbend strain', overStrain * 100, '%', 'D/(2R) on the stinger plus the axial strain from the lay tension; zero for methods without an overbend.', inp.method === 'slay' ? crit.overbendLimit * 100 : null, 'max'),
        kpi('lambda', 'Bending boundary-layer length at touchdown', sc.lambda0, 'm', '√(EI/H): the length over which bending stiffness modifies the catenary near the seabed. The model is valid while it is much shorter than the suspended length.'));
      if (inp.method === 'reel') {
        const rl = run(() => reelLay({ D, t, reelRadius: inp.reelR, fy, E }), 'Reel bending');
        checks.push(check('Reeling strain', rl.reelStrain * 100, inp.reelStrainMax, '%', 'DNV-ST-F101 (2021) Sec. 10 (installation by reeling; strain qualified by ECA)', 'Plastic bending strain on the reel hub. Above the qualified strain, girth-weld defects can grow and the pipe may buckle locally on the reel.'),
          check('Ovality after reeling', rl.ovalityF101 * 100, 3, '%', 'DNV-ST-F101 (2021) Sec. 5 (ovalisation limit 3 %)', 'Flattening of the cross-section from bending; it reduces collapse resistance in deep water.'));
        kpis.push(kpi('reelStrain', 'Reeling strain on the hub', rl.reelStrain * 100, '%', '(D/2)/(R + D/2): nominal bending strain of the pipe wrapped on the reel.', inp.reelStrainMax, 'max'),
          kpi('ovality', 'Ovality after reeling', rl.ovalityF101 * 100, '%', 'DNV-ST-F101 bending ovalisation expression evaluated at the reeling strain.', 3, 'max'));
        rows.push(['Residual radius after spring-back', rl.residualRadius, 'm'], ['Reel bending moment', rl.reelMoment / 1e3, 'kN·m']);
      }
      // Tension needed to bring the sagbend stress to its limit (bisection on H; stress falls as H rises).
      let lo = 0.02 * w * h, hi = 50 * w * h;
      const ucAt = (HH) => { const s2 = stiffenedCatenary({ EI: E * I, H: HH, w, h, D, EA, topFree: inp.method !== 'jlay' }); const c2 = layCriteriaF101({ grade: inp.grade, criterion: 1, overbendStrain: 1e-9, D, t, fy, sagbendTe: HH, sagbendMoment: s2.sagbendMoment, pe }); return c2.ucSagbend; };
      try { if (ucAt(lo) > 1 && ucAt(hi) < 1) { for (let i = 0; i < 40; i++) { const mid = Math.sqrt(lo * hi); if (ucAt(mid) > 1) lo = mid; else hi = mid; } const Tmin = hi + w * h + (inp.method === 'slay' ? w * Math.max(-inp.hTop, 0) : 0); kpis.push(kpi('TminSag', 'Minimum lay tension for the sagbend criterion', Tmin / 1e3, 'kN', 'Tensioner tension at which the sagbend equivalent stress just reaches its limit. Laying below it over-bends the sagbend.')); tReqNote = `Minimum tension for the sagbend: ${fmt(Tmin / 1e3)} kN.`; } } catch { /* optional figure */ }
    }
    checks.push(check('Tensioner capacity', T / 1e3, inp.Tcap, 'kN', 'Lay equipment rating; DNV-ST-F101 Sec. 10 (installation equipment)', 'The tensioners must hold the static lay tension with margin for dynamics; if they slip the pipe or cable runs away.'));
    const dyn = run(() => layDynamicTension({ staticTension: T, w, EA, m: mLin, suspendedLength: susp, topAngleDeg: angDeg, heaveAmp: inp.heave, period: inp.Tmot }), 'Dynamic tension');
    checks.push(check('Dynamic tension (quasi-static bound) against tensioner capacity', (T + dyn.lower) / 1e3, inp.Tcap, 'kN', 'Quasi-dynamic lower bound: ΔT = w·a_v + m·Ls·ω²·a_t (first principles; DNV-ST-F101 Sec. 10 requires dynamics to be included)', 'Least peak tension to expect when the suspended span follows the vessel motion. The drag-locked upper bound is reported as a KPI: if that one exceeds the rating, only a dynamic lay analysis can clear the sea state.'));
    if (dyn.lower >= T) checks.push(check('Loss of tension at the touchdown (quasi-static bound)', dyn.lower / 1e3, T / 1e3, 'kN', 'Loss of tension criterion (DNV-ST-F101 Sec. 10; DNV-RP-0360 §4 for cables)', 'The dynamic variation exceeds the static tension: the line goes slack and can buckle or loop at the touchdown.'));
    const straddle = T + dyn.upper > inp.Tcap * 1e3 || dyn.minTensionUpper < 0;
    const gov = worst(checks), ok = gov.uc <= 1;
    const days = inp.length / inp.layRate, steel_t = steelPerKm * inp.length;
    const quantities = isCable ? { cable_km: round(inp.length), vessel_days: { cableLay: round(days) } } : { pipe_km: round(inp.length), steel_t: round(steel_t), vessel_days: { pipelay: round(days) } };
    return {
      kpis: [
        kpi('Ttop', 'Top (tensioner) tension', T / 1e3, 'kN', 'Static tension held by the lay equipment.', inp.Tcap, 'max'),
        kpi('Hbottom', 'Bottom (horizontal) tension', Tbot / 1e3, 'kN', 'Tension at the touchdown point: top tension less the submerged weight of the vertical height carried.', 0, 'min'),
        kpi('layBack', 'Lay-back', layBack, 'm', 'Horizontal distance from the departure point to the touchdown point.'),
        kpi('susp', 'Suspended length', susp, 'm', 'Length of product hanging between the departure point and the seabed.'),
        kpi('angle', 'Departure angle from the horizontal', angDeg, 'deg', 'Angle of the catenary at the departure point; it must match the stinger tip, tower or chute geometry.'),
        kpi('w', 'Submerged weight', w, 'N/m', 'Weight in water per metre, empty pipe with coating or cable.'),
        ...kpis,
        kpi('dTlow', 'Dynamic tension amplitude, lower bound', dyn.lower / 1e3, 'kN', 'The span follows the motion quasi-statically: weight of the extra length lifted plus axial inertia.'),
        kpi('dTup', 'Dynamic tension amplitude, upper bound', dyn.upper / 1e3, 'kN', 'The span is locked by drag and the motion is taken by elastic stretch. The real amplitude lies between the bounds.'),
        kpi('days', 'Lay duration', days, 'days', 'Length over lay rate, without weather downtime.'),
      ],
      checks,
      plots: [
        line('Curvature along the suspended span', 'Arc length from touchdown', 'm', 'Curvature', '1/km', [{ name: 'Curvature', x: sigma.map((v) => round(v)), y: curv.map((v) => round(v * 1e3, 5)) }],
          'Curvature rises from zero at the touchdown (bending boundary layer) to its sagbend peak and falls towards the top. The peak sets the sagbend strain or the bend radius.'),
        line('Tension components against lay tension', 'Tensioner tension', 'kN', 'Bottom tension (kN) · touchdown radius (m)', '', (() => { const ts = linspace(1.05, 3, 14).map((f) => f * w * h + (T - H - w * h)); return [{ name: 'Bottom tension (kN)', x: ts.map((v) => round(v / 1e3)), y: ts.map((v) => round((v - (T - H)) / 1e3)) }, { name: 'Touchdown radius H/w (m)', x: ts.map((v) => round(v / 1e3)), y: ts.map((v) => round((v - (T - H)) / w)) }]; })(),
          'More lay tension raises the bottom tension one-for-one and opens the touchdown radius, at the price of longer lay-back and higher residual tension on the seabed.'),
      ],
      tables: [tbl('Lay configuration', [['k', 'Quantity'], ['v', 'Value'], ['u', 'Unit']], [
        ['Height of the departure point above the seabed', h, 'm'], ['Catenary parameter a = H/w', H / w, 'm'], ['Axial stiffness EA', EA / 1e6, 'MN'], ['Mass in air', mLin, 'kg/m'], ['Minimum tension with upper-bound dynamics', dyn.minTensionUpper / 1e3, 'kN'], ...rows,
      ].map(([k, v, u]) => ({ k, v: round(v), u })), 'Geometry and stiffness values behind the lay analysis, for comparison with the lay vessel configuration sheet.')],
      notes: [`${METHODS.find((m) => m.value === inp.method).label}: top ${fmt(T / 1e3)} kN, bottom ${fmt(Tbot / 1e3)} kN, lay-back ${fmt(layBack)} m in ${inp.depth} m of water.`, ...(tReqNote ? [tReqNote] : []), `"${gov.name}" governs at ${fmt(gov.uc)}.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 || straddle ? 'warn' : 'ok', `The lay configuration is within its limits; "${gov.name}" governs at ${fmt(gov.uc)}. Dynamic tension lies between ±${fmt(dyn.lower / 1e3)} and ±${fmt(dyn.upper / 1e3)} kN for ${inp.heave} m of motion.`,
          `${straddle ? 'The drag-locked upper bound exceeds the tensioner rating or unloads the touchdown, so the two bounds straddle the limit: run a dynamic lay analysis' : 'Confirm with a dynamic lay analysis'} for the vessel RAOs before fixing the limiting sea state, and set the tensioner dead-band so that the bottom tension stays between ${fmt(0.8 * Tbot / 1e3)} and ${fmt(1.2 * Tbot / 1e3)} kN.`, isCable ? 'DNV-RP-0360 (2016) §4' : 'DNV-ST-F101 (2021) Sec. 10')
          : adv('engineering', 'fail', `"${gov.name}" is exceeded (${fmt(gov.uc)}).`,
            gov.name.startsWith('Sagbend') ? `Increase the lay tension${tReqNote ? ` (${tReqNote.toLowerCase()})` : ''} or the wall thickness.` : gov.name.startsWith('Overbend') ? `Use a stinger radius of at least ${fmt(inp.Rst * gov.uc)} m or reduce the tension.` : gov.name.startsWith('Chute') ? `Use a chute of at least ${inp.mbr} m radius.` : gov.name.startsWith('Touchdown') ? `Raise the lay tension so that the bottom tension is at least ${fmt(inp.mbr * w / 1e3)} kN.` : gov.name.startsWith('Sidewall') ? `Enlarge the chute radius to ${fmt(inp.chuteR * gov.uc)} m or lower the tension.` : 'Reduce the tension demand (shallower departure angle, lighter product) or use equipment with more capacity; lower the limiting sea state to cut the dynamic part.',
            isCable ? 'DNV-RP-0360 (2016) §4' : 'DNV-ST-F101 (2021) Sec. 10'),
        adv('economic', 'info', `Laying ${inp.length} km at ${inp.layRate} km/day takes about ${fmt(days)} vessel days before weather downtime.`,
          'Lay rate, not tension, drives the cost: check that the limiting sea state implied by the dynamic tension does not cut workability below the planned rate, and use the weather-window module for the downtime allowance.',
          'Lay spread day rate × duration (AACE RP 18R-97 estimate classes); weather downtime per DNV-ST-N001 Sec. 2',
          `${isCable ? 'Lay campaign' : `Line pipe steel ≈ ${steelCo2Text(steel_t)} Lay campaign`}: ${vesselCo2Text(days, isCable ? 25 : 45)}`),
      ],
      publish: { 'install.layTension': T, 'install.layBack': layBack },
      quantities,
    };
  },
  verify: [
    { name: 'Cable catenary: H = T − w·h, lay-back = a·acosh(1 + h/a)', source: 'Inextensible catenary (Irvine 1981 §1.2): w = 450 N/m, chute 8 m above water in 120 m depth → h = 128 m, T = 120 kN → H = 62.4 kN, a = 138.67 m, lay-back = 138.67·acosh(1 + 128/138.67) = 176.3 m, S = √(h·(h + 2a)) = 227.8 m',
      inputs: { method: 'cable', Ttop: 120, hTop: 8, Tcap: 400 }, expect: [{ kpi: 'Hbottom', value: 62.4, tol: 1e-9 }, { kpi: 'layBack', value: (62400 / 450) * Math.acosh(1 + 128 / (62400 / 450)), tol: 1e-6 }, { kpi: 'susp', value: Math.sqrt(128 * (128 + 2 * 62400 / 450)), tol: 1e-6 }, { kpi: 'tdpR', value: 62400 / 450, tol: 1e-9 }] },
    { name: 'Reeling strain: ε = (D/2)/(R + D/2)', source: 'Geometry: 323.9 mm pipe on an 8.2 m hub → 0.16195/8.36195 = 1.937 %',
      inputs: { method: 'reel', D: 323.9, t: 17.5, tCoat: 3, rhoCoat: 900, hTop: 10, Ttop: 600 }, expect: [{ kpi: 'reelStrain', value: 100 * 0.16195 / 8.36195, tol: 1e-6 }] },
    { name: 'S-lay overbend strain: ε = D/(2R) + T/EA', source: 'DNV-ST-F101 Sec. 10 simplified overbend strain: D = 508 mm, R = 180 m → 0.1411 %; T/EA = 900e3/(207e9·0.024581) = 0.0177 % → 0.1588 %',
      inputs: {}, expect: [{ kpi: 'overStrain', value: 100 * (0.508 / 360 + 900e3 / (207e9 * Math.PI / 4 * (0.508 ** 2 - 0.4762 ** 2))), tol: 1e-6 }] },
  ],
};
