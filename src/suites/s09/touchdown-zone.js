// Touchdown-zone and seabed interaction of a steel catenary riser: touchdown-point movement, sag-bend curvature,
// seabed reaction, penetration, trench effect, cyclic stiffness degradation and pipe–soil resistance.
import { riserOffsetCases, boundaryLayerTdp, trenchProfile } from '../../engine/lines/riser.js';
import { embedment, axialResistance, lateralResistance } from '../../engine/pipeline/pipesoil.js';
import {
  DEG, need, run, fmt, round, scaled, adv, check, kpi, plot, meshInput, finish, handWeight, steelTonnes, steelCo2Text, vesselCo2Text, LIMIT_2D, CONFIRM,
} from './_common.js';
import { scrInputs, buildScr, tdpNode } from './_scr.js';

const clay = { key: 'soilModel', in: ['undrained'] }, sand = { key: 'soilModel', in: ['drained'] };
const A12 = 985 * Math.sin(12 * DEG) / (1 - Math.sin(12 * DEG)); // catenary parameter of the default riser [m]

/** Chainage / elevation arrays from an imported polyline (x, y, z points) or the first row of a bathymetry grid. */
function surveyProfile(model) {
  if (!model) return null;
  const xs = [], zs = [];
  if (model.kind === 'grid' && model.z?.length >= 2) for (let i = 0; i < model.nx; i++) { xs.push(i * model.dx); zs.push(model.z[i]); }
  else if (model.points?.length >= 6) { const p = model.points; let c = 0; for (let i = 0; i < p.length / 3; i++) { if (i > 0) c += Math.hypot(p[3 * i] - p[3 * i - 3], p[3 * i + 1] - p[3 * i - 2]); xs.push(c); zs.push(p[3 * i + 2]); } }
  if (xs.length < 2 || !zs.every(Number.isFinite)) throw new Error('The imported seabed profile needs at least two points with finite elevations.');
  for (let i = 1; i < xs.length; i++) if (!(xs[i] > xs[i - 1])) throw new Error('The imported seabed profile must advance in chainage (no repeated or vertical points).');
  return { xs, zs };
}

export default {
  id: 's09-touchdown-zone', suite: 9, title: 'Touchdown zone and seabed interaction',
  summary: 'Touchdown-point movement, touchdown-zone curvature, seabed contact force, penetration and friction for a catenary riser on a flat or trenched seabed.',
  covers: [7], fidelity: 'numerical',
  standards: ['DNV-RP-F114 (2017) §4 pipe–soil interaction', 'DNV-ST-F201 (2018) Sec. 4', 'API RP 2RD (1998, R2013) §5'],
  models: ['Winkler seabed springs with lift-off (contact–separation)', 'Coulomb friction with elastic shear stiffness', 'Stiffened-catenary boundary layer (Plunkett 1967)', 'Pipe penetration resistance (undrained / drained bearing)', 'Raised-cosine trench profile', 'Secant stiffness degradation'],
  limits: `${LIMIT_2D} Soil response is linear-elastic with lift-off: hysteretic load–unload loops, suction during uplift and trench formation are not simulated — degradation enters as a factor on the secant stiffness and suction as a screening force. The trench shape is an assumption unless survey data sets it.`,
  publishes: ['riser.tdzCurvature', 'riser.tdpRange', 'riser.embedment'],
  searchTerms: 'steel catenary riser touchdown zone seabed interaction trench soil stiffness',
  inputs: [
    ...scrInputs(),
    { key: 'ktAx', label: 'Axial seabed shear stiffness', unit: 'kPa', type: 'number', default: 100, min: 0, max: 1e4, group: 'Bottom end (seabed)', role: 'soil',
      help: 'Axial force per metre of pipe per metre of slip before sliding. 0 switches axial friction off (the far-end anchor then carries the tension).' },
    { key: 'soilModel', label: 'Soil response', type: 'select', default: 'undrained', group: 'Seabed soil', role: 'soil', options: [{ value: 'undrained', label: 'Clay (undrained)' }, { value: 'drained', label: 'Sand (drained)' }] },
    { key: 'su0', label: 'Undrained shear strength at mudline', unit: 'kPa', type: 'number', default: 2, min: 0.2, max: 200, group: 'Seabed soil', role: 'soil', bind: 'site.su0', showIf: clay, uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'suGrad', label: 'Strength gradient', unit: 'kPa/m', type: 'number', default: 1.5, min: 0, max: 20, group: 'Seabed soil', role: 'soil', showIf: clay },
    { key: 'alphaAdh', label: 'Axial adhesion factor α', unit: '-', type: 'number', default: 0.5, min: 0.1, max: 1, group: 'Seabed soil', role: 'soil', showIf: clay, help: 'Interface shear / undrained strength; 0.2–1.0 depending on coating roughness and sensitivity.' },
    { key: 'phiDeg', label: 'Friction angle', unit: 'deg', type: 'number', default: 32, min: 20, max: 45, group: 'Seabed soil', role: 'soil', bind: 'site.phiDeg', showIf: sand },
    { key: 'gammaSub', label: 'Submerged unit weight', unit: 'kN/m³', type: 'number', default: 6, min: 2, max: 12, group: 'Seabed soil', role: 'soil', bind: 'site.gammaSub' },
    { key: 'muLat', label: 'Lateral friction coefficient', unit: '-', type: 'number', default: 0.5, min: 0.1, max: 1.5, group: 'Seabed soil', role: 'soil' },
    { key: 'kDyn', label: 'Dynamic embedment factor', unit: '-', type: 'number', default: 1.5, min: 1, max: 4, group: 'Seabed soil', role: 'soil', help: 'Multiplies static embedment for riser motion at the touchdown (1–3 on clay, site-specific).' },
    { key: 'knDegr', label: 'Cyclic stiffness degradation factor', unit: '-', type: 'number', default: 0.3, min: 0.02, max: 1, group: 'Seabed soil', role: 'soil',
      help: 'Degraded secant stiffness / initial stiffness after cyclic loading. The riser is re-solved with the degraded value to bracket the response.' },
    { key: 'fSuction', label: 'Suction factor on uplift', unit: '-', type: 'number', default: 0.2, min: 0, max: 1, group: 'Seabed soil', role: 'soil', help: 'Uplift resistance as a fraction of the penetration resistance during fast lift-off (screening).' },
    { key: 'trenchDepth', label: 'Trench depth', unit: 'm', type: 'number', default: 0, min: 0, max: 5, group: 'Seabed profile', role: 'bc', help: '0 = flat seabed. Observed SCR trenches are typically 2–5 pipe diameters deep.' },
    { key: 'trenchLength', label: 'Trench length', unit: 'm', type: 'number', default: 150, min: 20, max: 600, group: 'Seabed profile', role: 'bc', help: 'Centred on the nominal touchdown point.' },
    { key: 'seabedProfile', label: 'Surveyed seabed profile along the riser plane (optional)', type: 'geometry', kinds: ['polyline', 'grid'], default: null, group: 'Seabed profile', role: 'bc',
      help: 'A polyline (chainage from the hang-off, elevation) or a bathymetry grid sampled along its first row replaces the parametric trench. Elevations are taken relative to the level at the nominal touchdown point. Leave empty for the flat or trenched seabed above.' },
    { key: 'offset', label: 'Vessel offset (near / far)', unit: 'm', type: 'number', default: 40, min: 1, max: 300, group: 'Top end (hang-off)', role: 'bc', bind: 'vessel.offsetMax', uncertain: { dist: 'normal', cov: 0.2 } },
    { key: 'flowlineLength', label: 'Flowline length available for anchoring', unit: 'm', type: 'number', default: 2000, min: 50, max: 2e5, group: 'Bottom end (seabed)', role: 'geometry', bind: 'pipeline.length',
      help: 'Length of grounded pipe beyond the touchdown point whose axial friction can react the bottom tension (to the PLET or first anchor).' },
    { key: 'strainAllow', label: 'Allowable bending strain', unit: '%', type: 'number', default: 0.15, min: 0.05, max: 0.5, group: 'Acceptance criteria', role: 'criteria', help: `Elastic design limit for the sag bend. ${CONFIRM}` },
    meshInput(160),
  ],
  mesh: { key: 'nElem', label: 'Beam elements', levels: [80, 160, 320, 640], qoi: 'kappaMax' },

  compute(inp, ctx) {
    // engineering units → SI inside buildScr (mm, bar, kPa) and below (kPa → Pa, kN/m³ → N/m³)
    const hh = inp.depth - inp.hangoffDepth;
    need(inp.offset <= 0.25 * hh, `A vessel offset of ${inp.offset} m exceeds 25 % of the hang-off height (${fmt(hh)} m). Reduce the offset.`);
    let seabed = {};
    const survey = surveyProfile(inp.seabedProfile);
    if (survey || inp.trenchDepth > 0) {
      const x0 = buildScr(inp, ctx, { nElem: 80 }).base.catenary.horizontalSpan;
      if (survey) {
        // imported profile: linear interpolation in chainage, datum shifted to the level at the nominal touchdown
        const at = (x) => { const { xs, zs } = survey, n = xs.length; if (x <= xs[0]) return zs[0]; if (x >= xs[n - 1]) return zs[n - 1]; let k = 0; while (xs[k + 1] < x) k++; return zs[k] + (zs[k + 1] - zs[k]) * (x - xs[k]) / (xs[k + 1] - xs[k]); };
        const z0 = at(x0);
        need(survey.xs[survey.xs.length - 1] > x0, `The imported seabed profile ends at chainage ${fmt(survey.xs[survey.xs.length - 1])} m, short of the touchdown at about ${fmt(x0)} m from the hang-off.`);
        // the far end of the model is pinned at the datum level: fade the imported relief out over the last 30 % of
        // the grounded length so the boundary does not sit above or below the seabed
        const xEnd = x0 + inp.groundedLength, xFade = xEnd - 0.3 * inp.groundedLength;
        const fade = (x) => (x <= xFade ? 1 : x >= xEnd ? 0 : 0.5 * (1 + Math.cos(Math.PI * (x - xFade) / (xEnd - xFade))));
        seabed = { profile: (x) => (at(x) - z0) * fade(x) };
      } else seabed = { profile: trenchProfile({ zSeabed: 0, xStart: x0 - inp.trenchLength / 2, length: inp.trenchLength, depth: inp.trenchDepth }) };
    }
    const { pp, base, h } = buildScr(inp, ctx, { nElem: inp.nElem, seabed });
    const w = pp.section.w, EI = pp.section.EI, kn = inp.kn * 1e3, n = base.s.length - 1;
    const offs = [-inp.offset, -inp.offset / 2, inp.offset / 2, inp.offset];
    const oc = run(() => riserOffsetCases(base, offs), 'Vessel-offset cases did not solve');
    oc.forEach((o) => need(o.result.converged && o.result.tdp, `The offset case ${o.offset} m did not converge or lost seabed contact; reduce the offset or lengthen the grounded pipe.`));
    need(base.tdp, 'The riser does not touch the seabed: lengthen the pipe on the seabed.');
    const degr = buildScr(inp, ctx, { nElem: inp.nElem, seabed, knFactor: inp.knDegr }).base;
    const kT = tdpNode(base), H = base.Te[kT];
    need(H > 0, 'The touchdown point is in effective compression in the mean position; increase the hang-off angle.');
    const all = [...oc.slice(0, 2).map((o) => o.result), base, ...oc.slice(2).map((o) => o.result)];
    const kappaMax = Math.max(...all.map((r) => r.maxCurvature));
    const maxOf = (a) => { let m = 0; for (const v of a) if (v > m) m = v; return m; };
    const rMax = Math.max(...all.map((r) => maxOf(r.reaction)));
    const bl = boundaryLayerTdp({ EI, H, w });
    const Do = pp.D + 2 * inp.coatT / 1e3;
    const und = inp.soilModel === 'undrained';
    const kLay = Math.max(1, rMax / w);
    const embArgs = { D: Do, W: w, model: inp.soilModel, kDyn: inp.kDyn, su0: inp.su0 * 1e3, suGrad: inp.suGrad * 1e3, gammaSub: inp.gammaSub * 1e3, phiDeg: inp.phiDeg };
    let emb, embNote = '';
    try { emb = embedment({ ...embArgs, kLay }); } catch {
      // the concentrated reaction exceeds the bearing resistance within two diameters: report the static-weight
      // embedment and flag that the pipe will trench in until the reaction spreads
      emb = run(() => embedment({ ...embArgs, kLay: 1 }), 'Pipe embedment could not be found even under the static pipe weight; check the soil strength');
      embNote = `The peak reaction (${fmt(kLay)} × weight) exceeds the soil bearing resistance within two diameters: the embedment shown is for the static weight only and the pipe will trench in further.`;
    }
    const suZ = inp.su0 * 1e3 + inp.suGrad * 1e3 * emb.z;
    const ax = axialResistance({ model: inp.soilModel, D: Do, z: emb.z, W: w, mu: inp.mu, alpha: inp.alphaAdh, su: suZ });
    const lat = lateralResistance({ soil: und ? 'clay' : 'sand', D: Do, z: emb.z, W: w, mu: inp.muLat, su: suZ, gammaSub: inp.gammaSub * 1e3, phiDeg: inp.phiDeg });
    const uplift = inp.fSuction * emb.V;
    const Lg = base.s[n] - base.tdp.s;
    const strain = kappaMax * pp.D / 2;
    const tdpRange = oc[3].tdpX - oc[0].tdpX;
    const pen = rMax / kn;
    const kappaDeg = degr.maxCurvature;
    const sigB = pp.mat.E * strain;

    const kpis = [
      kpi('kappaMax', 'Maximum touchdown-zone curvature', kappaMax * 1e3, '1/km', 'Largest curvature in the sag bend over the near, mean and far positions. It sets bending strain and the fatigue stress range at the touchdown; the near case governs.'),
      kpi('kappaCat', 'Catenary curvature w/H at touchdown', bl.catenaryCurvature * 1e3, '1/km', 'Curvature a perfectly flexible catenary would have at touchdown. Bending stiffness and seabed compliance keep the real maximum slightly below it.'),
      kpi('strain', 'Maximum bending strain', strain * 100, '%', 'Outer-fibre strain κ·D/2 at the tightest point of the sag bend, compared with the elastic design limit.', inp.strainAllow, 'max'),
      kpi('bendStress', 'Bending stress at the sag bend', sigB / 1e6, 'MPa', 'Elastic bending stress E·κ·D/2 at the tightest bend. Its variation between vessel positions is what drives touchdown fatigue.'),
      kpi('tdpTension', 'Effective tension at touchdown (mean)', H / 1e3, 'kN', 'Horizontal tension where the pipe meets the seabed. The sag-bend radius is about H/w, so low tension means a tight bend.'),
      kpi('tdpX', 'Touchdown distance from hang-off (mean)', base.tdp.x, 'm', 'Horizontal distance from the hang-off to the first point of seabed contact in the mean vessel position.'),
      kpi('tdpRange', 'Touchdown movement near → far', tdpRange, 'm', 'Length of seabed over which the touchdown point travels between the near and far offsets. Fatigue damage and trenching are spread over this zone.'),
      kpi('lambda', 'Bending boundary-layer length', bl.lambda, 'm', 'Decay length √(EI/H) over which bending stiffness matters at the touchdown. The mesh must resolve it and the seabed reaction peak sits within it.'),
      kpi('reactionMax', 'Peak seabed reaction', rMax / 1e3, 'kN/m', 'Largest contact force per metre over all positions. It exceeds the pipe weight just behind the touchdown point because the seabed has to turn the shear force.'),
      kpi('kLay', 'Reaction concentration factor', kLay, '-', 'Peak seabed reaction divided by the submerged weight. It is the load factor used for embedment (k_lay in DNV-RP-F114).'),
      kpi('penetration', 'Elastic penetration at the peak', pen * 1e3, 'mm', 'Peak reaction divided by the seabed stiffness: the spring compression in the model, to be compared with the bearing-capacity embedment.'),
      kpi('embedment', 'Pipe embedment (bearing capacity)', emb.z * 1e3, 'mm', 'Depth at which the soil bearing resistance equals the concentrated touchdown load, including the dynamic factor. Deeper embedment raises lateral and uplift resistance.'),
      kpi('axialRes', 'Axial soil resistance', ax.breakout / 1e3, 'kN/m', 'Axial friction available per metre of grounded pipe at the calculated embedment; it anchors the horizontal tension along the flowline.'),
      kpi('lateralRes', 'Lateral soil resistance (breakout)', lat.breakout / 1e3, 'kN/m', 'Lateral resistance per metre including passive soil ahead of the embedded pipe. It restrains out-of-plane touchdown motion, which this 2-D model does not solve.'),
      kpi('upliftRes', 'Uplift (suction) resistance', uplift / 1e3, 'kN/m', 'Screening estimate of the suction that resists fast lift-off, as a fraction of the penetration resistance. Not applied in the solver; it adds bending during uplift.'),
      kpi('kappaDegraded', 'Curvature with degraded seabed stiffness', kappaDeg * 1e3, '1/km', 'Maximum mean-position curvature when the seabed stiffness is reduced by the cyclic degradation factor. A softer seabed spreads the reaction and usually lowers the curvature peak.'),
    ];
    const anchor = ax.breakout * Math.max(Lg, inp.flowlineLength);
    const checks = [
      check('Sag-bend bending strain', strain * 100, inp.strainAllow, '%', 'DNV-ST-F201 (2018) Sec. 5 (elastic design); project strain limit', 'Guards against yielding and local buckling of the pipe in the sag bend. Above the limit the bend is no longer elastic and ratcheting or wrinkling can develop.'),
      check('Axial anchoring of the bottom tension', H / 1e3, anchor / 1e3, 'kN', 'DNV-RP-F114 (2017) §4.3 axial resistance', 'Compares the horizontal tension at touchdown with the friction available along the grounded flowline. If friction is insufficient the flowline walks towards the riser unless an anchor or PLET holds it.'),
    ];

    // touchdown-zone window for the plots
    const sT = base.tdp.s, win = Math.max(12 * bl.lambda, 0.6 * Math.abs(tdpRange) + 4 * bl.lambda);
    const pick = (r, f, ref = sT) => { const xs = [], ys = []; for (let k = 0; k <= n; k++) if (Math.abs(r.s[k] - ref) <= win) { xs.push(r.s[k] - ref); ys.push(f(r, k)); } return { xs, ys }; };
    const names = [`Near −${inp.offset} m`, `Near −${inp.offset / 2} m`, 'Mean', `Far +${inp.offset / 2} m`, `Far +${inp.offset} m`];
    const curv = [0, 2, 4].map((i) => { const q = pick(all[i], (r, k) => Math.abs(r.curvature[k]) * 1e3); return { name: names[i], x: q.xs, y: q.ys }; });
    const xiA = []; for (let i = 0; i <= 40; i++) xiA.push(win * i / 40);
    const blc = boundaryLayerTdp({ EI, H, w, xi: xiA }).curvature;
    curv.push({ name: 'Boundary-layer solution (rigid seabed, mean)', x: xiA.map((v) => -v), y: scaled(blc, 1e3), dash: true });
    const reac = [0, 2, 4].map((i) => { const q = pick(all[i], (r, k) => r.reaction[k] / 1e3); return { name: names[i], x: q.xs, y: q.ys }; });
    reac.push({ name: 'Submerged weight', x: [-win, win], y: [w / 1e3, w / 1e3], dash: true });
    const shp = [0, 2, 4].map((i) => { const xs = [], ys = []; const r = all[i]; for (let k = 0; k <= n; k++) if (r.z[k] < 12 && Math.abs(r.x[k] - base.tdp.x) < 3 * win) { xs.push(r.x[k]); ys.push(r.z[k]); } return { name: names[i], x: xs, y: ys }; });
    const bx = []; const by = []; for (let i = 0; i <= 60; i++) { const x = base.tdp.x - 3 * win + 6 * win * i / 60; bx.push(x); by.push(seabed.profile ? seabed.profile(x) : 0); }
    shp.push({ name: 'Seabed', x: bx, y: by, dash: true });
    const plots = [
      plot('Touchdown-zone shape', 'Horizontal distance from mean hang-off', 'm', 'Height above seabed', 'm', shp, 'Pipe axis close to the seabed for the near, mean and far positions (vertical scale exaggerated). Look at where each curve meets the seabed and how the trench, if any, changes the contact.'),
      plot('Curvature through the touchdown zone', 'Arc length from the mean touchdown point (suspended side negative)', 'm', 'Curvature', '1/km', curv, 'Curvature magnitude either side of the touchdown point. It rises from zero on the seabed to the catenary value over a few boundary-layer lengths; the dashed curve is the analytical stiffened-catenary solution.'),
      plot('Seabed reaction', 'Arc length from the mean touchdown point', 'm', 'Reaction', 'kN/m', reac, 'Contact force per metre. The spike just behind the touchdown point is the concentrated reaction (≈ w·λ on a rigid seabed); further along it settles to the submerged weight.'),
      plot('Touchdown position against vessel offset', 'Vessel offset (far positive)', 'm', 'Touchdown x relative to mean', 'm', [{ name: 'Touchdown excursion', x: [offs[0], offs[1], 0, offs[2], offs[3]], y: [oc[0].tdpExcursion, oc[1].tdpExcursion, 0, oc[2].tdpExcursion, oc[3].tdpExcursion] }], 'How far the touchdown point moves along the seabed per metre of vessel offset. The slope (typically above 1) shows how vessel motion is amplified at the seabed.'),
      plot('Maximum curvature against vessel offset', 'Vessel offset (far positive)', 'm', 'Maximum curvature', '1/km', [{ name: 'Maximum curvature', x: [offs[0], offs[1], 0, offs[2], offs[3]], y: [oc[0], oc[1], { maxCurvature: base.maxCurvature }, oc[2], oc[3]].map((o) => o.maxCurvature * 1e3) }], 'Sag-bend curvature for each vessel position. The steep rise towards the near side is why near offsets and low tension govern touchdown strength and fatigue.'),
    ];
    const tables = [
      { title: 'Vessel positions', columns: [{ key: 'case', label: 'Position' }, { key: 'off', label: 'Offset', unit: 'm' }, { key: 'tdp', label: 'Touchdown excursion', unit: 'm' }, { key: 'H', label: 'Bottom tension', unit: 'kN' }, { key: 'kappa', label: 'Max curvature', unit: '1/km' }, { key: 'R', label: 'Peak reaction', unit: 'kN/m' }],
        rows: all.map((r, i) => ({ case: names[i], off: [offs[0], offs[1], 0, offs[2], offs[3]][i], tdp: round(i === 2 ? 0 : oc[i < 2 ? i : i - 1].tdpExcursion, 4), H: round(r.Te[tdpNode(r)] / 1e3, 4), kappa: round(r.maxCurvature * 1e3, 4), R: round(maxOf(r.reaction) / 1e3, 4) })),
        explain: 'One row per vessel position: where the touchdown point sits, the tension left at the seabed and the resulting curvature and contact force.' },
      { title: 'Pipe–soil resistance at the touchdown', columns: [{ key: 'item', label: 'Quantity' }, { key: 'value', label: 'Value' }, { key: 'unit', label: 'Unit' }],
        rows: [{ item: 'Embedment', value: round(emb.z * 1e3, 4), unit: 'mm' }, { item: 'Embedment / diameter', value: round(emb.zOverD, 3), unit: '-' }, { item: 'Design vertical load', value: round(emb.V / 1e3, 4), unit: 'kN/m' },
          { item: 'Axial resistance (breakout)', value: round(ax.breakout / 1e3, 4), unit: 'kN/m' }, { item: 'Lateral resistance (breakout)', value: round(lat.breakout / 1e3, 4), unit: 'kN/m' }, { item: 'Lateral resistance (residual)', value: round(lat.residual / 1e3, 4), unit: 'kN/m' },
          { item: 'Uplift (suction) resistance', value: round(uplift / 1e3, 4), unit: 'kN/m' }, { item: 'Frictional resistance over the flowline length', value: round(anchor / 1e3, 4), unit: 'kN' }],
        explain: 'Soil resistances at the calculated embedment, for use as spring limits in a detailed touchdown model and for the flowline anchoring check.' },
    ];
    const ucMax = Math.max(...checks.map((c) => c.uc));
    const steel = steelTonnes(pp.D, pp.t, Math.abs(tdpRange) + 2 * win, pp.mat.rho);
    const advice = [
      adv('engineering', ucMax > 1 ? 'fail' : ucMax > 0.85 ? 'warn' : 'ok',
        `Sag-bend strain ${fmt(strain * 100)} % against ${inp.strainAllow} %; touchdown point travels ${fmt(tdpRange)} m; peak reaction ${fmt(kLay)} × weight; degraded seabed changes the mean curvature from ${fmt(base.maxCurvature * 1e3)} to ${fmt(kappaDeg * 1e3)} 1/km.`,
        checks[0].uc > 1 ? `Raise the hang-off angle (now ${inp.hangoffAngle}°): bottom tension must rise to about ${fmt(w * pp.D / 2 / (inp.strainAllow / 100) / 1e3)} kN to keep κ = w/H within the strain limit, or limit the near offset.`
          : checks[1].uc > 1 ? `Provide a hold-back anchor or PLET rated for at least ${fmt((H - anchor) / 1e3)} kN, or confirm at least ${fmt(H / ax.breakout)} m of grounded flowline to mobilise enough friction.`
            : 'Carry the range of seabed stiffness (initial to degraded) and the trench profile into the fatigue analysis: touchdown fatigue is more sensitive to them than static strength is.',
        'DNV-RP-F114 (2017) §4 for pipe–soil parameters; DNV-ST-F201 (2018) Sec. 5. Soil factors are site-specific and must be confirmed by the geotechnical report.'),
      adv('economic', 'info',
        `Touchdown-zone damage concentrates in about ${fmt(Math.abs(tdpRange) + 2 * win)} m of pipe (${fmt(steel)} t). Upgrading only this length (thicker wall, ground flush welds or clad pipe) is far cheaper than upgrading the whole riser.`,
        'Specify the fatigue-critical weld quality and any wall upset for the touchdown window only, and budget one ROV survey of the trench after the first year to confirm the seabed assumptions before fixing the inspection interval.',
        'Risk-based inspection (DNV-RP-G101 / API RP 580) and whole-life costing (ISO 15686-5).',
        `Local upgrade affects ${steelCo2Text(steel)} A trench survey by ROV support vessel: ${vesselCo2Text(1)}`),
    ];
    if (embNote) advice.push(adv('engineering', 'warn', embNote, 'Use a softer (degraded) seabed stiffness consistent with trench formation, or obtain cyclic pipe–soil test data for the touchdown zone.', 'DNV-RP-F114 (2017) §4.2 (embedment under lay and dynamic loads).'));
    return finish({
      kpis, checks, plots, tables,
      notes: [...(embNote ? [embNote] : []), `Seabed stiffness ${inp.kn} kPa, ${survey ? 'imported seabed profile' : inp.trenchDepth > 0 ? `trench ${inp.trenchDepth} m deep over ${inp.trenchLength} m` : 'flat seabed'}; ${inp.ktAx > 0 ? `axial friction μ = ${inp.mu} on a ${inp.ktAx} kPa shear spring` : 'no axial friction'}.`,
        `Analytical boundary layer: λ = ${fmt(bl.lambda)} m, concentrated reaction w·λ = ${fmt(bl.pointReaction / 1e3)} kN. ${checks[0].uc >= checks[1].uc ? 'Bending strain' : 'Axial anchoring'} governs.`],
      advice,
      publish: { 'riser.tdzCurvature': kappaMax, 'riser.tdpRange': Math.abs(tdpRange), 'riser.embedment': emb.z },
      quantities: { steel_t: steelTonnes(pp.D, pp.t, base.s[n], pp.mat.rho), vessel_days: { survey: 1 } },
    });
  },

  verify: [
    {
      // Catenary on a flat seabed: a = h·sinθ/(1 − sinθ) = 985·sin12°/(1 − sin12°) = 258.5 m,
      // bottom tension H = w·a, touchdown distance X = a·asinh(S/a) with S = sqrt(h(h + 2a)).
      // The stiff pipe touches down slightly further out than the cable, hence the 2 % tolerance on X.
      name: 'Flat-seabed catenary: H = w·a and touchdown distance a·asinh(S/a)', source: 'Irvine, Cable Structures (1981) §1.2',
      inputs: { topType: 'pinned', pInt: 0, ktAx: 0, trenchDepth: 0, nElem: 320 },
      expect: [{ kpi: 'tdpTension', value: handWeight() * A12 / 1e3, tol: 0.01 },
        { kpi: 'tdpX', value: A12 * Math.asinh(Math.sqrt(985 * (985 + 2 * A12)) / A12), tol: 0.02 },
        { kpi: 'kappaCat', value: 1e3 / A12, tol: 0.01 }],
    },
  ],
};
