// Foundation settlement: immediate, primary consolidation (layered, time-dependent), secondary compression and tilt.
import { elasticSettlement, footingRotation, consolidationSettlement, consolidationFD, consolidationU, consolidationTv, secondaryCompression, biotColumn } from '../../engine/geotech/settlement.js';
import { janbuModulus } from '../../engine/geotech/soilmodels.js';
import { DEG, need, run, fmt, round, adv, check, vesselCo2Text, matInputs, NOT_FE } from './_common.js';

const YEAR = 31557600;

const DEFAULT_LAYERS = [
  { name: 'Soft clay', thickness: 4, gamma: 6, Cc: 0.45, Cr: 0.07, e0: 1.6, ocr: 1.5, cv: 1.5 },
  { name: 'Firm clay', thickness: 8, gamma: 7.5, Cc: 0.30, Cr: 0.05, e0: 1.1, ocr: 2.5, cv: 3 },
  { name: 'Stiff clay', thickness: 12, gamma: 9, Cc: 0.20, Cr: 0.03, e0: 0.8, ocr: 5, cv: 6 },
];

export default {
  id: 's04-settlement', suite: 4, title: 'Foundation settlement and differential settlement',
  summary: 'Immediate, consolidation and secondary settlement of the mudmat through layered soil, its development in time, and the rotation and differential settlement caused by moment.',
  covers: [6], fidelity: 'numerical',
  standards: ['ISO 19901-4:2016 §7.5', 'API RP 2GEO (2011) §7.4', 'DNV-RP-C212 (2019) §5.6'],
  models: ['Elastic settlement theory (rigid rectangle on half-space)', 'Boussinesq stress distribution', 'One-dimensional Terzaghi consolidation (layered finite-volume)', 'e–log σ′ compression with pre-consolidation', 'Secondary compression', 'Biot poroelastic column (cross-check)', 'Janbu stress-dependent stiffness', 'Rigid footing rotation on elastic half-space'],
  limits: `One-dimensional vertical drainage and strain below the mat centre; stress spread by Boussinesq or the 2:1 rule; rotation from an elastic half-space with an equivalent secant modulus. Three-dimensional consolidation, creep coupling and lateral soil flow are not modelled. ${NOT_FE}`,
  publishes: ['subsea.settlementTotal', 'subsea.tiltDeg'],
  searchTerms: 'mudmat consolidation settlement clay subsea foundation tilt',
  inputs: [
    ...matInputs(),
    { key: 'D', label: 'Base level below mudline (skirt depth)', unit: 'm', type: 'number', default: 0.5, min: 0, max: 10, group: 'Mudmat geometry', role: 'geometry' },
    { key: 'V', label: 'Sustained vertical load', unit: 'kN', type: 'number', default: 1000, min: 1, max: 2e5, group: 'Sustained loads', role: 'bc', bind: 'subsea.designV', uncertain: { dist: 'normal', cov: 0.05 }, help: 'Submerged weight plus sustained interface reactions (not storm peaks).' },
    { key: 'M', label: 'Sustained moment (eccentricity along B)', unit: 'kN·m', type: 'number', default: 400, min: 0, max: 1e6, group: 'Sustained loads', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Long-term moment from eccentric equipment weight and pipeline reactions; it tilts the mat.' },
    { key: 'layers', label: 'Soil layers below the base', type: 'table', group: 'Soil', role: 'soil',
      columns: [{ key: 'name', label: 'Layer', type: 'text' }, { key: 'thickness', label: 'Thickness', unit: 'm', type: 'number' }, { key: 'gamma', label: 'γ′', unit: 'kN/m³', type: 'number' }, { key: 'Cc', label: 'Cc', type: 'number' }, { key: 'Cr', label: 'Cr', type: 'number' }, { key: 'e0', label: 'e0', type: 'number' }, { key: 'ocr', label: 'OCR', type: 'number' }, { key: 'cv', label: 'cv', unit: 'm²/yr', type: 'number' }],
      default: DEFAULT_LAYERS, help: 'From the base downward. Compression indices and cv from oedometer tests; OCR from the pre-consolidation pressure.' },
    { key: 'Eu', label: 'Undrained Young’s modulus', unit: 'MPa', type: 'number', default: 3, min: 0.1, max: 2000, group: 'Soil', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.4 }, help: 'For immediate (distortion) settlement with ν = 0.5; typically 200–500 × su.' },
    { key: 'Calpha', label: 'Secondary compression index Cα', unit: '-', type: 'number', default: 0.008, min: 0, max: 0.1, group: 'Soil', role: 'soil', help: 'Δe per log cycle of time after end of primary; about 0.03–0.05 × Cc for clays.' },
    { key: 'stress', label: 'Stress distribution', type: 'select', default: 'boussinesq', group: 'Soil', role: 'control', options: [{ value: 'boussinesq', label: 'Boussinesq (below centre)' }, { value: '2:1', label: '2:1 load spread' }] },
    { key: 'top', label: 'Drainage at the base of the mat', type: 'select', default: 'drained', group: 'Drainage boundaries', role: 'bc', options: [{ value: 'drained', label: 'Drained (vented / permeable interface)' }, { value: 'undrained', label: 'Impervious (sealed plate)' }] },
    { key: 'bottom', label: 'Drainage at the bottom of the profile', type: 'select', default: 'undrained', group: 'Drainage boundaries', role: 'bc', options: [{ value: 'drained', label: 'Drained (sand layer below)' }, { value: 'undrained', label: 'Impervious' }] },
    { key: 'U0', label: 'Consolidation already achieved at the start', unit: '%', type: 'number', default: 0, min: 0, max: 95, group: 'Initial state', role: 'ic', help: 'Start state of the simulation: 0 for a new installation; use the degree reached when assessing an existing structure. Settlement reported is the remaining part.' },
    { key: 'life', label: 'Design life', unit: 'yr', type: 'number', default: 25, min: 0.5, max: 100, group: 'Time', role: 'operation', bind: 'life.designYears' },
    { key: 'tEval', label: 'Time for the reported degree of consolidation', unit: 'yr', type: 'number', default: 1, min: 0.01, max: 100, group: 'Time', role: 'operation' },
    { key: 'sAllow', label: 'Allowable total settlement', unit: 'mm', type: 'number', default: 250, min: 5, max: 2000, group: 'Acceptance criteria', role: 'criteria', help: 'Set by jumper and connector stroke; project-specific.' },
    { key: 'tiltAllow', label: 'Allowable tilt', unit: 'deg', type: 'number', default: 1.0, min: 0.05, max: 10, group: 'Acceptance criteria', role: 'criteria', help: 'Typical connector and tree landing tolerance is 0.5–1.5°; project-specific.' },
    { key: 'nCells', label: 'Finite-volume cells per layer', type: 'int', default: 20, min: 4, max: 320, group: 'Solver', role: 'control' },
    { key: 'nSteps', label: 'Time steps', type: 'int', default: 400, min: 50, max: 6400, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nCells', label: 'Cells per layer', levels: [5, 10, 20, 40], qoi: 'Ueval' },
  meshTime: { key: 'nSteps', label: 'Time steps', levels: [100, 200, 400, 800], qoi: 'Ueval' },
  compute(inp) {
    // SI: kN → N, kPa/MPa → Pa, kN/m³ → N/m³, m²/yr → m²/s, yr → s.
    const layers = (inp.layers || []).filter((l) => l.thickness > 0);
    need(layers.length > 0, 'Give at least one soil layer with a positive thickness.');
    for (const l of layers) need(l.gamma > 0 && l.Cc >= 0 && l.Cr >= 0 && l.e0 > 0 && l.cv > 0 && (l.ocr ?? 1) >= 1, `Layer "${l.name}": γ′, e0 and cv must be positive, Cc and Cr non-negative and OCR ≥ 1.`);
    need(!(inp.top === 'undrained' && inp.bottom === 'undrained'), 'Both drainage boundaries are impervious: excess pore pressure cannot dissipate and no consolidation occurs. Make at least one boundary drained.');
    need(inp.tEval <= inp.life, `The evaluation time (${inp.tEval} yr) must not exceed the design life (${inp.life} yr).`);
    const b = Math.min(inp.B, inp.L), l = Math.max(inp.B, inp.L), A = b * l;
    const q = inp.V * 1e3 / A;
    const Eu = inp.Eu * 1e6;
    const imm = run(() => elasticSettlement({ q, B: b, L: l, E: Eu, nu: 0.5, shape: 'rect', point: 'rigid' }), 'Immediate settlement');
    const cs = run(() => consolidationSettlement({ q, B: b, L: l, stressMethod: inp.stress, sigmaV0Base: layers[0].gamma * 1e3 * inp.D,
      layers: layers.map((y) => ({ thickness: y.thickness, gamma: y.gamma * 1e3, Cc: y.Cc, Cr: y.Cr, e0: y.e0, ocr: y.ocr ?? 1 })), nSub: 12 }), 'Consolidation settlement');
    const sPrim = cs.s;
    // Equivalent compressibility per layer so that the column model returns the layer settlements: mv = s / (h·q).
    const Htot = layers.reduce((a, y) => a + y.thickness, 0);
    const fvLayers = layers.map((y, i) => { const mv = Math.max(cs.layers[i].s / (y.thickness * q), 1e-12); return { thickness: y.thickness, mv, k: (y.cv / YEAR) * mv * 9810 }; });
    const tEnd = inp.life * YEAR;
    const fd = run(() => consolidationFD({ layers: fvLayers, load: () => q, tEnd, nSteps: inp.nSteps, nCells: inp.nCells, top: inp.top, bottom: inp.bottom }), 'Consolidation in time');
    const U0 = inp.U0 / 100;
    const tY = Array.from(fd.t, (t) => t / YEAR), Uraw = Array.from(fd.U);
    // Start state: time axis shifted so that the simulation begins at the achieved degree of consolidation.
    let iStart = 0; while (iStart < Uraw.length - 1 && Uraw[iStart] < U0) iStart++;
    const interpU = (ty) => { const tt = ty + tY[iStart]; if (tt >= tY[tY.length - 1]) return Uraw[Uraw.length - 1]; let i = 1; while (tY[i] < tt) i++; const f = (tt - tY[i - 1]) / (tY[i] - tY[i - 1]); return Uraw[i - 1] + f * (Uraw[i] - Uraw[i - 1]); };
    const Ueval = interpU(inp.tEval), Uend = interpU(inp.life);
    const timeTo = (U) => { for (let i = 1; i < Uraw.length; i++) if (Uraw[i] >= U) { const f = (U - Uraw[i - 1]) / (Uraw[i] - Uraw[i - 1] || 1); return Math.max(tY[i - 1] + f * (tY[i] - tY[i - 1]) - tY[iStart], 0); } return NaN; };
    const t50 = timeTo(0.5), t90 = timeTo(0.9);
    // Secondary compression after end of primary (t at U = 95 % of the equivalent single layer).
    const Hdr = inp.top === 'drained' && inp.bottom === 'drained' ? Htot / 2 : Htot;
    const cvEq = Htot * Htot / layers.reduce((a, y) => a + y.thickness / Math.sqrt(y.cv / YEAR), 0) ** 2; // equivalent-thickness rule
    const tp = consolidationTv(0.95) * Hdr * Hdr / cvEq;
    const sSec = (ty) => (inp.Calpha > 0 && ty * YEAR > tp ? secondaryCompression({ Calpha: inp.Calpha, e0: layers[0].e0, H: Htot, t1: tp, t2: ty * YEAR }) : 0);
    const sAt = (ty) => (U0 > 0 ? 0 : imm.s) + (interpU(ty) - U0) * sPrim + sSec(ty);
    const sLife = sAt(inp.life);
    const sFinal = (U0 > 0 ? 0 : imm.s) + (1 - U0) * sPrim + sSec(inp.life);
    // Rotation: immediate with Eu, long-term with the secant modulus that reproduces immediate + primary settlement.
    const M = inp.M * 1e3;
    const Esec = Eu * imm.s / (imm.s + sPrim) * (1 - 0.3 * 0.3) / (1 - 0.25);
    const rot0 = footingRotation({ M, B: inp.B, L: inp.L, E: Eu, nu: 0.5 });
    const rotL = footingRotation({ M, B: inp.B, L: inp.L, E: Esec, nu: 0.3 });
    const tilt = rotL.theta / DEG;
    // Cross-check at tEval: equivalent single layer by Terzaghi and by the Biot column.
    const mvAvg = sPrim / (Htot * q), Kv = 1 / Math.max(mvAvg, 1e-12), nuD = 0.3;
    const UT = consolidationU(cvEq * inp.tEval * YEAR / (Hdr * Hdr));
    const biot = run(() => biotColumn({ q, h: Hdr, K: Kv * (1 + nuD) / (3 * (1 - nuD)), G: Kv * (1 - 2 * nuD) / (2 * (1 - nuD)), k: cvEq * mvAvg * 9810, gammaW: 9810, t: inp.tEval * YEAR }), 'Biot column');
    const janbu = janbuModulus({ m: 12, j: 0, sigma: Math.max(layers[0].gamma * 1e3 * (inp.D + layers[0].thickness / 2), 1e3) });
    const ts = [], ss = [], us = [];
    for (let i = 0; i <= 60; i++) { const ty = inp.life * (i / 60) ** 2; ts.push(round(ty, 5)); ss.push(sAt(ty) * 1e3); us.push(interpU(ty) * 100); }
    const checks = [
      check('Total settlement at end of design life', sLife * 1e3, inp.sAllow, 'mm', 'ISO 19901-4:2016 §7.5 (serviceability); project connector stroke', 'Guards against settlement exceeding what jumpers and connectors can follow. Beyond it rigid spools are overstressed at their hubs.'),
      check('Tilt from sustained moment', Math.abs(tilt), inp.tiltAllow, 'deg', 'ISO 13628-1 / API RP 17A installation and operating tolerances (project-specific)', 'Guards against rotation of the structure beyond connector alignment and tooling tolerance.'),
    ];
    const ok = checks.every((c) => c.uc <= 1);
    const worstLayer = cs.layers.reduce((a, c, i) => (c.s > a.s ? { s: c.s, i } : a), { s: -1, i: 0 });
    const advice = [
      adv('engineering', ok ? 'ok' : 'fail', `Settlement at ${fmt(inp.life)} years is ${fmt(sLife * 1e3)} mm (immediate ${fmt(imm.s * 1e3)}, primary ${fmt(sPrim * 1e3)} of which ${fmt(Uend * 100)} % has developed, secondary ${fmt(sSec(inp.life) * 1e3)}); tilt ${fmt(tilt)}°. Layer "${layers[worstLayer.i].name}" contributes most (${fmt(worstLayer.s * 1e3)} mm).`,
        ok ? `Give jumpers a vertical stroke of at least ${fmt(sLife * 1e3 * 1.3)} mm (30 % allowance) and schedule metrology after ${fmt(Number.isFinite(t50) ? t50 : inp.life)} years, when half the consolidation is complete.` : `Reduce the bearing pressure (${fmt(q / 1e3)} kPa): a mat ${fmt(Math.sqrt(sLife * 1e3 / inp.sAllow))} times larger in area brings settlement near the limit; or re-centre the load to cut the sustained moment (${fmt(inp.M)} kN·m) that drives the tilt.`,
        'ISO 19901-4:2016 §7.5; API RP 2GEO (2011) §7.4; Terzaghi one-dimensional consolidation',
        'Designing flexibility into the jumpers avoids an oversized mat; an intervention to re-level or re-make a connection costs vessel time: ' + vesselCo2Text(2)),
      adv('economic', 'info', `Half of the consolidation settlement takes ${fmt(Number.isFinite(t50) ? t50 : inp.life)} years and 90 % ${Number.isFinite(t90) ? fmt(t90) + ' years' : 'longer than the design life'}: connections made at installation see most of the movement.`,
        'Where settlement uses up connector stroke, compare a larger mat (steel at about 0.16 t/m²) with a planned re-metrology and spool adjustment campaign (about 2 vessel days) on a whole-life basis; oedometer tests that firm up cv and Cc narrow the choice.',
        'Whole-life costing per ISO 15686-5; risk-based inspection planning per DNV-RP-G101'),
    ];
    return {
      kpis: [
        { key: 'sLife', label: 'Settlement at end of design life', unit: 'mm', value: sLife * 1e3, limit: inp.sAllow, sense: 'max', explain: 'Immediate plus developed consolidation plus secondary settlement below the mat centre at the end of the design life (remaining settlement when a start state is given).' },
        { key: 'sImm', label: 'Immediate settlement', unit: 'mm', value: imm.s * 1e3, explain: 'Undrained elastic distortion settlement of the rigid rectangle, occurring as the load is applied.' },
        { key: 'sPrim', label: 'Final primary consolidation settlement', unit: 'mm', value: sPrim * 1e3, explain: 'Settlement from dissipation of excess pore pressure, summed over the layers from the e–log σ′ curves.' },
        { key: 'sFinal', label: 'Ultimate settlement', unit: 'mm', value: sFinal * 1e3, explain: 'Settlement once consolidation is complete, with secondary compression up to the design life.' },
        { key: 'Ueval', label: 'Degree of consolidation at the evaluation time', unit: '-', value: Ueval, explain: 'Share of the primary settlement developed at the chosen time, from the layered finite-volume solution.' },
        { key: 'Uend', label: 'Degree of consolidation at end of life', unit: '-', value: Uend, explain: 'Close to 1 means settlement is essentially complete within the design life; lower values mean movement continues throughout operation.' },
        { key: 't50', label: 'Time to 50 % consolidation', unit: 'yr', value: Number.isFinite(t50) ? t50 : inp.life, explain: 'Time for half the consolidation settlement. Reported as the design life when it is not reached within it.' },
        { key: 'tilt', label: 'Long-term tilt', unit: 'deg', value: tilt, limit: inp.tiltAllow, sense: 'max', explain: 'Rotation of the rigid mat under the sustained moment using the secant modulus that reproduces the long-term settlement.' },
        { key: 'sDiff', label: 'Differential settlement edge to edge', unit: 'mm', value: rotL.differential * 1e3, explain: 'Difference in settlement between opposite edges of the mat from the long-term rotation.' },
        { key: 'tilt0', label: 'Immediate rotation', unit: 'deg', value: rot0.theta / DEG, explain: 'Rotation on application of the moment with the undrained modulus; the part visible at landing.' },
        { key: 'q', label: 'Net bearing pressure', unit: 'kPa', value: q / 1e3, explain: 'Sustained vertical load over the gross mat area: the stress increment that drives settlement.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Settlement against time', x: { label: 'Time', unit: 'yr' }, y: { label: 'Settlement', unit: 'mm' },
          series: [{ name: 'Settlement', x: ts, y: ss }, { name: 'Allowable', x: [0, inp.life], y: [inp.sAllow, inp.sAllow], dash: true }],
          explain: 'Development of settlement over the design life. Look at how much occurs in the first year (before tie-in) and whether the curve crosses the allowable line.' },
        { type: 'line', title: 'Degree of consolidation against time', x: { label: 'Time', unit: 'yr' }, y: { label: 'Degree of consolidation', unit: '%' }, series: [{ name: 'Layered finite-volume', x: ts, y: us }],
          explain: 'Share of the primary consolidation completed. A curve still rising at the end of life means long-term movement must be accommodated.' },
        { type: 'line', title: 'Stress and strain below the mat centre', x: { label: 'Stress', unit: 'kPa' }, y: { label: 'Depth below base', unit: 'm' },
          series: [{ name: 'Initial effective stress', x: cs.profile.map((p) => p.sigma0 / 1e3), y: cs.profile.map((p) => -p.z) }, { name: 'Stress increase', x: cs.profile.map((p) => p.dSigma / 1e3), y: cs.profile.map((p) => -p.z) }],
          explain: 'In-situ vertical effective stress and the increment from the mat against depth. Settlement comes from the zone where the increment is a significant fraction of the in-situ stress.' },
        { type: 'line', title: 'Excess pore pressure at end of life', x: { label: 'Excess pore pressure', unit: 'kPa' }, y: { label: 'Depth below base', unit: 'm' }, series: [{ name: 'u', x: Array.from(fd.u, (u) => u / 1e3), y: Array.from(fd.z, (z) => -z) }],
          explain: 'Remaining excess pore pressure through the column at the end of the simulation. A bulge marks the layer that is still consolidating.' },
      ],
      tables: [
        { title: 'Settlement by layer', columns: [{ key: 'name', label: 'Layer' }, { key: 'h', label: 'Thickness', unit: 'm' }, { key: 's', label: 'Primary settlement', unit: 'mm' }, { key: 'mv', label: 'Equivalent mv', unit: '1/MPa' }, { key: 'k', label: 'Permeability', unit: 'm/s' }],
          rows: layers.map((y, i) => ({ name: y.name, h: y.thickness, s: round(cs.layers[i].s * 1e3), mv: round(fvLayers[i].mv * 1e6), k: round(fvLayers[i].k, 3) })),
          explain: 'Contribution of each layer to the primary settlement with the equivalent compressibility and permeability used in the time-domain column model.' },
        { title: 'Cross-check of the degree of consolidation at the evaluation time', columns: [{ key: 'model', label: 'Model' }, { key: 'U', label: 'U', unit: '-' }],
          rows: [{ model: 'Layered finite-volume column', U: round(Ueval) }, { model: 'Terzaghi single equivalent layer', U: round(UT) }, { model: 'Biot poroelastic column (incompressible constituents)', U: round(biot.U) }],
          explain: 'Three independent estimates of the same quantity. Agreement of the last two is exact by theory; their difference from the first shows the effect of layering.' },
      ],
      notes: [`Drainage: top ${inp.top}, bottom ${inp.bottom}; drainage path ${fmt(Hdr)} m; equivalent cv ${fmt(cvEq * YEAR)} m²/yr.`, `Long-term secant modulus ${fmt(Esec / 1e6)} MPa; Janbu modulus of a normally consolidated clay (m = 12) at mid-depth of the first layer would be ${fmt(janbu / 1e6)} MPa for comparison.`,
        checks[0].uc >= checks[1].uc ? 'Total settlement governs over tilt.' : 'Tilt governs over total settlement.'],
      advice,
      publish: { 'subsea.settlementTotal': sLife, 'subsea.tiltDeg': tilt },
      quantities: { vessel_days: { metrology: 1 } },
    };
  },
  verify: [
    { name: 'Terzaghi single layer, one-way drainage: U = 50 % at Tv = 0.197', source: 'Terzaghi (1943) series solution: H = 10 m, cv = 1 m²/yr → t = 0.197·100/1 = 19.7 yr gives U = 0.5007',
      inputs: { layers: [{ name: 'clay', thickness: 10, gamma: 6, Cc: 0.4, Cr: 0.05, e0: 1.5, ocr: 1, cv: 1 }], top: 'drained', bottom: 'undrained', life: 50, tEval: 19.7, nSteps: 1600, nCells: 80, U0: 0 },
      expect: [{ kpi: 'Ueval', value: 0.5007, tol: 0.01 }] },
    { name: 'Immediate settlement of a rigid square: s = q·B·(1 − ν²)·I/E with I = 0.93 × 0.946', source: 'Schleicher (1926) average factor 0.946 for a flexible square; rigid ≈ 0.93 × average (Bowles 1996 §5.6): 10e3·10·0.75·0.880/5e6 = 13.20 mm',
      inputs: { B: 10, L: 10, V: 1000, Eu: 5 }, expect: [{ kpi: 'sImm', value: 13.20, tol: 0.01 }] },
  ],
};
