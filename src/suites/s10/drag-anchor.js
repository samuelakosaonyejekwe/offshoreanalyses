// Drag-embedment anchor in clay: embedment trajectory, ultimate depth and capacity, drag distance and the
// embedded-chain attenuation (Neubecker & Randolph), cross-checked against the design-chart power law.
import { dragAnchorTrajectory, dragAnchorUhc, DRAG_ANCHOR_INDICATIVE } from '../../engine/geotech/anchors.js';
import { DEG, G, need, run, round, fmt, adv, check, worst, sev, thin, interp1, steelCo2Text, vesselCo2Text } from './_common.js';
import { anchorLoadInputs, clayInputs, capacityFactorInputs, requiredSf, anchorLoad, LIMIT_ANCHOR } from './_anchor.js';

const SOILS = Object.keys(DRAG_ANCHOR_INDICATIVE);

export default {
  id: 's10-drag-anchor',
  suite: 10,
  title: 'Drag anchor embedment and capacity',
  summary: 'Predicts how deep and how far a drag-embedment anchor travels in clay as it is tensioned, the capacity it reaches, and whether the design load is met without uplift.',
  covers: [10],
  fidelity: 'reduced-order',
  standards: ['DNV-RP-E301 (2017) Design and installation of fluke anchors', 'API RP 2SK (3rd ed., 2005) §6.5 and anchor annex', 'ISO 19901-7 (2013) §10.4'],
  models: ['Drag-anchor kinematic trajectory model (Neubecker & Randolph 1996)', 'Embedded anchor-line inverse catenary and frictional capacity (Neubecker & Randolph 1995)', 'Undrained bearing on the fluke, Ta = Ne·Af·su', 'Design-chart power law UHC = a·W^b', 'Consolidation set-up and cyclic degradation factors'],
  limits: `${LIMIT_ANCHOR} The trajectory model is the simplified kinematic form: the anchor travels parallel to its fluke, the resultant acts at a fixed angle to the fluke, and the anchor line enters the seabed horizontally. Capacity factor Ne and angle θw are anchor-specific and must come from the manufacturer or calibrated installation records; sand, layered soil and anchor structural strength are not covered.`,
  searchTerms: 'drag embedment anchor trajectory clay Neubecker Randolph',
  inputs: [
    ...anchorLoadInputs(),
    { key: 'Af', label: 'Fluke plan area', unit: 'm²', type: 'number', default: 15, min: 1, max: 150, group: 'Anchor', role: 'geometry', design: true, help: 'Projected area of the fluke; roughly 1.2–1.5 m² per tonne^(2/3) for high-holding-power anchors.' },
    { key: 'mass', label: 'Anchor mass', unit: 't', type: 'number', default: 22, min: 0.5, max: 200, group: 'Anchor', role: 'geometry', help: 'Used for the chart cross-check and the steel quantity.' },
    { key: 'Ne', label: 'Anchor capacity factor Ne', unit: '-', type: 'number', default: 13.5, min: 3, max: 30, group: 'Anchor', role: 'soil', uncertain: { dist: 'normal', cov: 0.15 },
      help: 'Ta = Ne·Af·su at the fluke. Form factor 1.4–1.6 times Nc = 9 is reported for fluke anchors; anchor-specific.' },
    { key: 'thetaWDeg', label: 'Resultant-to-fluke angle θw', unit: 'deg', type: 'number', default: 25, min: 10, max: 60, group: 'Anchor', role: 'geometry', help: 'Depends on the fluke–shank angle (about 25–45°); larger for the mud setting (50° shank) than the sand setting (32°).' },
    { key: 'dLine', label: 'Anchor-line diameter (forerunner)', unit: 'mm', type: 'number', default: 95, min: 20, max: 250, group: 'Anchor line', role: 'geometry', help: 'Chain bar diameter or wire diameter of the embedded forerunner. Wire cuts deeper than chain.' },
    { key: 'lineKind', label: 'Forerunner type', type: 'select', default: 'chain', options: [{ value: 'chain', label: 'Chain (En = 2.5, μ = 0.4)' }, { value: 'wire', label: 'Wire rope (En = 1.0, μ = 0.2)' }], group: 'Anchor line', role: 'geometry' },
    ...clayInputs(),
    { key: 'dragSoil', label: 'Soil class for chart cross-check', type: 'select', default: SOILS[0], options: SOILS.map((v) => ({ value: v, label: v })), group: 'Soil', role: 'soil' },
    { key: 'installT', label: 'Installation (proof) tension at mudline', unit: 'kN', type: 'number', default: 3200, min: 10, max: 60000, group: 'Installation', role: 'ic',
      help: 'Tension applied when the anchor is set. Initial condition of the anchor: it fixes the installed depth and drag before service.' },
    ...capacityFactorInputs({ sf: 1.5 }),
    { key: 'nSteps', label: 'Trajectory integration steps', type: 'int', default: 400, min: 50, max: 3200, group: 'Numerical', role: 'control', help: 'Steps along the embedment path. 400 gives the drag distance to well under 1 %.' },
  ],
  mesh: { key: 'nSteps', label: 'Trajectory steps', levels: [100, 200, 400, 800], qoi: 'dragDistance' },
  compute(inp) {
    const ld = anchorLoad(inp);
    const sf = requiredSf(inp, 'dragAnchor');
    const su0 = inp.su0 * 1e3, k = inp.suGrad * 1e3;
    need(su0 > 0 || k > 0, 'The clay has no strength: enter a mudline strength or a gradient.');
    const En = inp.lineKind === 'chain' ? 2.5 : 1.0, mu = inp.lineKind === 'chain' ? 0.4 : 0.2;
    const tr = run(() => dragAnchorTrajectory({ Af: inp.Af, Ne: inp.Ne, thetaW: inp.thetaWDeg * DEG, d: inp.dLine / 1e3, En, Nc: 7.6, mu, su0, k, fraction: 0.95, n: inp.nSteps }), 'Drag-anchor trajectory');
    const chart = dragAnchorUhc({ mass: inp.mass * 1e3, ...DRAG_ANCHOR_INDICATIVE[inp.dragSoil] });
    const fac = inp.setUp * inp.cyclic;
    // installed state: depth and drag at which the mudline tension equals the installation tension
    const Ti = inp.installT * 1e3;
    const reached = Ti >= tr.Tm[tr.Tm.length - 1];
    const zInst = reached ? tr.z[tr.z.length - 1] : interp1(tr.Tm, tr.z, Ti);
    const xInst = reached ? tr.x[tr.x.length - 1] : interp1(tr.Tm, tr.x, Ti);
    // capacity available without further drag is the installation tension (times set-up and cyclic factors);
    // the ultimate capacity needs additional drag to mobilise
    const capInstalled = Math.min(Ti, tr.TmUlt) * fac, capUlt = tr.TmUlt * fac;
    const demand = ld.T * sf;
    const zReq = demand / fac >= tr.Tm[tr.Tm.length - 1] ? null : interp1(tr.Tm, tr.z, demand / fac);
    const xReq = zReq == null ? null : interp1(tr.Tm, tr.x, demand / fac);
    const checks = [
      check('Ultimate holding capacity', demand / 1e3, capUlt / 1e3, 'kN', `DNV-RP-E301 (2017) §4; API RP 2SK (3rd ed.) factor ${sf}`,
        'Guards against the anchor dragging without limit: the factored mudline load must not exceed the ultimate capacity at full embedment (after set-up and cyclic factors).'),
      check('Capacity proven by installation tension', ld.T / 1e3, capInstalled / 1e3, 'kN', 'DNV-RP-E301 (2017) §5 (installation tension vs design load)',
        'Guards against additional drag in service: if the unfactored design load exceeds what was proven at installation the anchor will move further before it holds, slackening the line and shifting the mooring pattern.'),
      { name: 'Uplift at the mudline', demand: round(ld.angleDeg, 4), capacity: 0, unit: 'deg', uc: ld.V > 1e-6 * ld.T ? 99 : 0, ref: 'API RP 2SK (3rd ed.) §6.5',
        explain: 'The model assumes the line enters the seabed horizontally; uplift at the mudline reduces embedment and can pull a fluke anchor out.' },
    ];
    const g = worst(checks);
    const [x, z, Tm, Ta] = thin([tr.x, tr.z, tr.Tm, tr.Ta], 200);
    const dAH = Math.max(0.5, xInst / 600 + 0.5);                 // screening: drag at ~25 m/h plus half a day handling per anchor
    return {
      kpis: [
        { key: 'zUlt', label: 'Ultimate fluke embedment', value: round(tr.zUlt, 5), unit: 'm', explain: 'Depth at which the anchor stops diving (the padeye line angle equals the resultant angle) and travels horizontally at constant capacity.' },
        { key: 'TaUlt', label: 'Ultimate capacity at the padeye', value: round(tr.TaUlt / 1e3, 5), unit: 'kN', explain: 'Fluke bearing resistance Ne·Af·su at the ultimate depth; the load the anchor body itself can take.' },
        { key: 'TmUlt', label: 'Ultimate holding capacity at the mudline', value: round(tr.TmUlt / 1e3, 5), unit: 'kN', explain: 'Padeye capacity plus the friction carried by the embedded anchor line; this is the UHC to compare with the mudline design load.' },
        { key: 'dragDistance', label: 'Drag to 95 % of ultimate depth', value: round(tr.dragDistance, 5), unit: 'm', explain: 'Horizontal travel needed to approach the ultimate depth. Long drag (tens of fluke lengths) means the UHC is rarely reached in practice.' },
        { key: 'zInstalled', label: 'Embedment at installation tension', value: round(zInst, 5), unit: 'm', explain: 'Fluke depth reached when the installation tension is applied; the starting state of the anchor in service.' },
        { key: 'xInstalled', label: 'Drag at installation tension', value: round(xInst, 5), unit: 'm', explain: 'Distance the anchor drags during setting; allow for it when placing the anchor relative to its target position.' },
        { key: 'zRequired', label: 'Embedment required for the factored load', value: round(zReq ?? 999, 5), unit: 'm', limit: round(0.95 * tr.zUlt, 5), sense: 'max', explain: 'Depth at which the mudline capacity equals the design load times the safety factor. 999 means the load exceeds what this anchor can reach.' },
        { key: 'fs', label: 'Factor of safety on ultimate capacity', value: round(capUlt / ld.T, 4), unit: '-', limit: sf, sense: 'min', explain: 'Ultimate mudline capacity (with set-up and cyclic factors) over the design load.' },
        { key: 'chartUhc', label: 'Design-chart UHC (power law)', value: round(chart.UHC / 1e3, 5), unit: 'kN', explain: 'Indicative chart value for the entered anchor mass and soil class, as an independent order-of-magnitude check on the trajectory model.' },
        { key: 'lineShare', label: 'Share of capacity carried by the embedded line', value: round(100 * (1 - tr.TaUlt / tr.TmUlt), 4), unit: '%', explain: 'Portion of the mudline capacity that comes from soil friction on the buried forerunner rather than the fluke.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Embedment trajectory', x: { label: 'Drag distance', unit: 'm' }, y: { label: 'Fluke depth below mudline', unit: 'm' },
          series: [{ name: 'Trajectory', x: x.map((v) => round(v)), y: z.map((v) => round(-v)) }, { name: 'Ultimate depth', x: [0, x[x.length - 1]], y: [-tr.zUlt, -tr.zUlt], dash: true }, { name: 'Installed', x: [round(xInst)], y: [round(-zInst)] }],
          explain: 'Path of the fluke as the anchor is dragged. It dives steeply at first and flattens towards the ultimate depth; the marker shows where installation leaves it.' },
        { type: 'line', title: 'Capacity mobilised with drag', x: { label: 'Drag distance', unit: 'm' }, y: { label: 'Tension', unit: 'kN' },
          series: [{ name: 'Mudline tension', x: x.map((v) => round(v)), y: Tm.map((v) => round(v / 1e3)) }, { name: 'Padeye tension', x: x.map((v) => round(v)), y: Ta.map((v) => round(v / 1e3)) },
            { name: 'Design load × SF / factors', x: [0, x[x.length - 1]], y: [demand / fac / 1e3, demand / fac / 1e3], dash: true }, { name: 'Installation tension', x: [0, x[x.length - 1]], y: [Ti / 1e3, Ti / 1e3], dash: true }],
          explain: 'Tension needed to keep the anchor moving, at the mudline and at the padeye. Where the curve crosses a dashed line the anchor has dragged far enough to hold that load.' },
      ],
      tables: [
        { title: 'Trajectory stations', columns: [{ key: 'x', label: 'Drag', unit: 'm' }, { key: 'z', label: 'Depth', unit: 'm' }, { key: 'Tm', label: 'Mudline tension', unit: 'kN' }, { key: 'Ta', label: 'Padeye tension', unit: 'kN' }, { key: 'th', label: 'Line angle at padeye', unit: 'deg' }],
          rows: thin([tr.x, tr.z, tr.Tm, tr.Ta, tr.thetaA], 16)[0].map((_, i, a) => { const c = thin([tr.x, tr.z, tr.Tm, tr.Ta, tr.thetaA], 16); return { x: round(c[0][i]), z: round(c[1][i]), Tm: round(c[2][i] / 1e3), Ta: round(c[3][i] / 1e3), th: round(c[4][i] / DEG, 3) }; }),
          explain: 'Depth, tensions and padeye line angle at intervals along the drag path; use it to read the installation tension needed for a target embedment.' },
      ],
      notes: [`Ultimate depth ${fmt(tr.zUlt)} m (${fmt(tr.zUlt / Math.sqrt(inp.Af))} fluke lengths) with mudline capacity ${fmt(tr.TmUlt / 1e3)} kN; chart value for ${inp.mass} t is ${fmt(chart.UHC / 1e3)} kN.`,
        `Installation at ${inp.installT} kN leaves the fluke at ${fmt(zInst)} m after ${fmt(xInst)} m of drag.`,
        'Ne, θw and the chart coefficients are indicative and anchor-specific; calibrate against installation records of the actual anchor type.'],
      advice: [
        adv('engineering', sev(g.uc), `Ultimate mudline capacity ${fmt(capUlt / 1e3)} kN gives a factor of safety of ${fmt(capUlt / ld.T)} against ${sf} required; installation tension proves ${fmt(capInstalled / 1e3)} kN against a design load of ${fmt(ld.T / 1e3)} kN (governing: ${g.name}, ${fmt(g.uc)}).`,
          g.uc > 1 ? (g.name.startsWith('Ultimate') ? `Increase the fluke area to about ${fmt(inp.Af * g.uc)} m² (capacity ∝ Af and deeper embedment), or use a wire forerunner to dive deeper.` : g.name.startsWith('Capacity') ? `Raise the installation tension to at least ${fmt(ld.T / fac / 1e3)} kN at the mudline, or accept and quantify the additional drag (${xReq != null ? fmt(Math.max(0, xReq - xInst)) : 'not reachable'} m to reach the factored load).` : 'Add ground chain so the line stays horizontal at the anchor, or change anchor type.')
            : `Record installation tension, drag and (where possible) depth for every anchor and compare with ${fmt(zInst)} m / ${fmt(xInst)} m predicted; a shortfall signals harder soil or a non-diving anchor.`,
          'DNV-RP-E301 (2017) §4–5; API RP 2SK (3rd ed., 2005) §6.5.'),
        adv('economic', 'info', `Setting each anchor needs about ${fmt(xInst)} m of drag at ${inp.installT} kN — roughly ${fmt(dAH)} anchor-handling vessel day(s) per anchor including handling (screening).`,
          `A proof tension of ${inp.installT} kN corresponds to about ${fmt(inp.installT / 9.81)} t bollard pull at the mudline; where that exceeds the vessel, cross-tension opposite lines from the platform instead of chartering a larger tug.`,
          'Marine-operations planning per DNV-ST-N001 (2021); AACE RP 18R-97 Class 4 quantities.',
          `${vesselCo2Text(dAH)} ${steelCo2Text(inp.mass)} The drag path disturbs a seabed strip about ${fmt(xInst * Math.sqrt(inp.Af))} m² per anchor.`),
      ],
      quantities: { steel_t: round(inp.mass), vessel_days: { anchor_setting: round(dAH, 3) } },
    };
  },
  verify: [
    { // Uniform clay su = 10 kPa. At the ultimate depth θa = θw, and the embedded-chain equation Ta·θa²/2 = z·En·d·Nc·su with
      // Ta = Ne·Af·su gives z_ult = Ne·Af·θw²/(2·En·d·Nc) = 13.5·20·0.61087²/(2·2.5·0.095·7.6) = 100.75/3.61 = 27.91 m.
      // Ta,ult = 13.5·20·10 = 2700 kN; Tm,ult = Ta·exp(μ·θw) = 2700·exp(0.4·0.61087) = 3447 kN.
      name: 'Ultimate embedment and capacity in uniform clay (closed form)',
      source: 'Neubecker & Randolph (1995) J. Geotech. Eng. 121(11) Eqs. for the embedded chain; (1996) Can. Geotech. J. 33 for θa = θw at ultimate depth',
      inputs: { Af: 20, Ne: 13.5, thetaWDeg: 35, dLine: 95, lineKind: 'chain', su0: 10, suGrad: 0 },
      expect: [{ kpi: 'zUlt', value: 27.91, tol: 0.005 }, { kpi: 'TaUlt', value: 2700, tol: 0.003 }, { kpi: 'TmUlt', value: 3447, tol: 0.005 }] },
  ],
};
