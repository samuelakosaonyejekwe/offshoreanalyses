// Mooring interference and clearance: line-to-line and line-to-asset distances at the calm-water and loaded
// vessel positions, anchor stand-off from seabed assets, crossings in plan, and the loss of clearance between
// neighbouring lines in a cross-current with wake shielding.
import { clearances, polylineMinDistance, lineForces } from '../../engine/lines/mooringsystem.js';
import { clearanceUnderCurrent } from '../../engine/lines/interference.js';
import { DEG, need, run, round, fmt, adv, check, worst, sev, systemInputs, loadInputs, meanState, vesselCo2Text, LIMIT_QS, LIMIT_PROPS } from './_common.js';

const ASSETS = [
  { name: 'Export pipeline', x1: 260, y1: -300, z1: -299.7, x2: 260, y2: 900, z2: -299.7, radius: 0.3 },
  { name: 'Production riser', x1: 0, y1: -15, z1: -22, x2: -160, y2: -230, z2: -300, radius: 0.25 },
  { name: 'Subsea template', x1: -420, y1: -620, z1: -294, x2: -400, y2: -600, z2: -294, radius: 12 },
];
const sysInputs = systemInputs().map((i) => (i.key === 'layout' ? { ...i, default: 'clustered' } : i));

/** Do two plan segments (p1–p2, p3–p4) intersect? */
function crosses(p1, p2, p3, p4) {
  const d = (a, b, c) => (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
  const d1 = d(p3, p4, p1), d2 = d(p3, p4, p2), d3 = d(p1, p2, p3), d4 = d(p1, p2, p4);
  return d1 * d2 < 0 && d3 * d4 < 0;
}

export default {
  id: 's10-interference-clearance',
  suite: 10,
  title: 'Mooring interference and clearance',
  summary: 'Measures how close the mooring lines come to each other and to pipelines, risers and structures as the vessel moves, and flags crossings and insufficient anchor stand-off.',
  covers: [14],
  fidelity: 'numerical',
  standards: ['API RP 2SK (3rd ed., 2005) §6.6 (clearances)', 'DNV-OS-E301 (2018) Ch.2 Sec.2 [5] (mooring line clearance)', 'DNV-RP-F203 (2009) riser interference (wake models)', 'ISO 19901-7 (2013) §10.6'],
  models: ['Three-dimensional catenary line profiles from the quasi-static system', 'Minimum distance between polylines (segment–segment closest approach)', 'Tensioned-string deflection of a line in steady current (Sparks 2007)', 'Wake shielding of the downstream line (Huse 1993; Blevins 2005)', 'Plan-view crossing detection'],
  limits: `${LIMIT_QS} Lines and assets are polylines with an exclusion radius: clearances are static at the two evaluated vessel positions, without wave-frequency line motion, vortex-induced vibration or contact mechanics. The current check treats the suspended part of the closest pair as tensioned strings in a uniform cross-flow. Required clearances are project values to confirm with the pipeline and field operators. ${LIMIT_PROPS}`,
  searchTerms: 'mooring line clearance pipeline crossing interference wake shielding',
  inputs: [
    ...sysInputs,
    ...loadInputs(),
    { key: 'assets', label: 'Subsea assets (straight segments)', type: 'table', group: 'Assets and obstacles', role: 'geometry',
      columns: [{ key: 'name', label: 'Asset', type: 'text' }, { key: 'x1', label: 'x start', unit: 'm', type: 'number' }, { key: 'y1', label: 'y start', unit: 'm', type: 'number' }, { key: 'z1', label: 'z start', unit: 'm', type: 'number' },
        { key: 'x2', label: 'x end', unit: 'm', type: 'number' }, { key: 'y2', label: 'y end', unit: 'm', type: 'number' }, { key: 'z2', label: 'z end', unit: 'm', type: 'number' }, { key: 'radius', label: 'Exclusion radius', unit: 'm', type: 'number' }],
      default: ASSETS,
      help: 'Pipelines, risers, umbilicals and structures as straight segments in global coordinates (z up from the still-water level, seabed at −depth). The exclusion radius is the asset half-width plus any protection zone.' },
    { key: 'assetGeom', label: 'Asset route (imported)', type: 'geometry', kinds: ['polyline'], default: null, group: 'Assets and obstacles', role: 'geometry',
      help: 'Optional surveyed pipeline or cable route as a 3-D polyline (x, y, z with z up). It is checked in addition to the table.' },
    { key: 'assetGeomRadius', label: 'Exclusion radius of the imported route', unit: 'm', type: 'number', default: 0.5, min: 0, max: 100, group: 'Assets and obstacles', role: 'geometry' },
    { key: 'Ucross', label: 'Cross-current on the lines', unit: 'm/s', type: 'number', default: 1.0, min: 0, max: 3.5, group: 'Current on lines', role: 'bc', bind: 'site.Uc100',
      help: 'Uniform current normal to the plane of the closest pair of lines (boundary condition for the interference check).' },
    { key: 'CdLine', label: 'Line drag coefficient', unit: '-', type: 'number', default: 2.4, min: 0.5, max: 4, group: 'Current on lines', role: 'load', help: 'On the nominal diameter: studless chain 2.4, wire rope 1.2–1.8, fibre rope 1.6.' },
    { key: 'wake', label: 'Wake model', type: 'select', default: 'huse', options: [{ value: 'huse', label: 'Huse (1993)' }, { value: 'blevins', label: 'Blevins (2005)' }, { value: 'none', label: 'No shielding' }], group: 'Current on lines', role: 'control' },
    { key: 'reqLineLine', label: 'Required line-to-line clearance', unit: 'm', type: 'number', default: 2, min: 0, max: 100, group: 'Acceptance criteria', role: 'criteria', help: 'Minimum surface distance between two mooring lines of the same unit away from the fairleads (contact must be avoided; project value).' },
    { key: 'reqLineAsset', label: 'Required line-to-asset clearance', unit: 'm', type: 'number', default: 10, min: 0, max: 200, group: 'Acceptance criteria', role: 'criteria', help: 'API RP 2SK guidance: about 10 m vertical clearance of a line over a pipeline in the intact condition; no contact in the damaged condition. Confirm with the asset owner.' },
    { key: 'reqAnchor', label: 'Required anchor stand-off from seabed assets', unit: 'm', type: 'number', default: 100, min: 0, max: 1000, group: 'Acceptance criteria', role: 'criteria', help: 'Horizontal distance of drag anchors from pipelines and structures (API RP 2SK guidance: 100 m, more where the line pulls towards the asset).' },
    { key: 'nProfile', label: 'Profile points per line segment', type: 'int', default: 20, min: 5, max: 160, group: 'Numerical', role: 'control', help: 'Resolution of each line polyline used for the distance search.' },
  ],
  mesh: { key: 'nProfile', label: 'Profile points per segment', levels: [10, 20, 40], qoi: 'minLineLine' },
  compute(inp) {
    const { b, sys, eq, calm } = meanState(inp);
    const assets = (inp.assets || []).filter((a) => a && a.name).map((a, i) => {
      for (const k of ['x1', 'y1', 'z1', 'x2', 'y2', 'z2']) need(Number.isFinite(a[k]), `Asset '${a.name}' (row ${i + 1}): coordinate ${k} is missing.`);
      return { name: a.name, points: [[a.x1, a.y1, a.z1], [a.x2, a.y2, a.z2]], radius: Math.max(a.radius || 0, 0) };
    });
    if (inp.assetGeom) {
      const p = inp.assetGeom.points || inp.assetGeom.positions;
      need(p && p.length >= 6, 'The imported asset route needs at least two points.');
      const pts = []; for (let i = 0; i + 2 < p.length; i += 3) pts.push([p[i], p[i + 1], p[i + 2]]);
      assets.push({ name: inp.assetGeom.name || 'Imported route', points: pts, radius: inp.assetGeomRadius });
    }
    const poseC = [calm.x, calm.y, 0, 0, 0, 0], poseL = [eq.x, eq.y, 0, 0, 0, 0];
    const cC = run(() => clearances(sys, { pose: poseC, nProfile: inp.nProfile, assets }), 'Clearance search (calm water)');
    const cL = run(() => clearances(sys, { pose: poseL, nProfile: inp.nProfile, assets }), 'Clearance search (loaded position)');
    const minOf = (arr, key) => arr.reduce((a, r) => (!a || r[key] < a[key] ? r : a), null);
    const llC = minOf(cC.lineLine, 'distance'), llL = minOf(cL.lineLine, 'distance');
    const laC = minOf(cC.lineAsset, 'clearance'), laL = minOf(cL.lineAsset, 'clearance');
    const llMin = Math.min(llC.distance, llL.distance), llGov = llL.distance < llC.distance ? llL : llC;
    const laMin = assets.length ? Math.min(laC.clearance, laL.clearance) : 9999, laGov = assets.length ? (laL.clearance < laC.clearance ? laL : laC) : null;
    // nominal line diameter for surface clearance
    const dLine = Math.max(...b.segs.map((s) => s.d));
    // anchors to seabed assets (assets with any point within 15 m of the seabed)
    let anchorMin = 9999, anchorAt = null;
    for (const a of assets) {
      if (!a.points.some((p) => p[2] < -b.depth + 15)) continue;
      sys.lines.forEach((ln, i) => {
        const dist = polylineMinDistance([[ln.anchor[0], ln.anchor[1], 0], [ln.anchor[0], ln.anchor[1], 0]], a.points.map((p) => [p[0], p[1], 0])).distance - a.radius;
        if (dist < anchorMin) { anchorMin = dist; anchorAt = { line: ln.name, asset: a.name }; }
      });
    }
    // crossings in plan between lines (anchor–fairlead chords at the loaded position)
    const ends = sys.lines.map((ln) => [[ln.anchor[0], ln.anchor[1]], [eq.x + ln.fairlead[0], eq.y + ln.fairlead[1]]]);
    let nCross = 0; const crossPairs = [];
    for (let i = 0; i < ends.length; i++) for (let j = i + 1; j < ends.length; j++) if (crosses(ends[i][0], ends[i][1], ends[j][0], ends[j][1])) { nCross++; crossPairs.push(`${sys.lines[i].name} × ${sys.lines[j].name}`); }
    // --- cross-current on the closest pair: suspended parts as tensioned strings, downstream line shielded by the wake
    const iA = llGov.a, iB = llGov.b;
    const lA = eq.lines[iA], lB = eq.lines[iB];
    const sus = Math.max(Math.min(b.Ltot - lA.groundedLength, b.Ltot - lB.groundedLength), 10);
    const az = (ln) => Math.atan2(ln.anchor[1] - ln.fairlead[1], ln.anchor[0] - ln.fairlead[0]);
    const dth = Math.abs(((az(sys.lines[iA]) - az(sys.lines[iB]) + 3 * Math.PI) % (2 * Math.PI)) - Math.PI);
    const rf = inp.fairleadRadius, Xs = Math.max(Math.sqrt(Math.max(sus * sus - b.hF * b.hF, 0)), 1);
    const nS = 60;
    const sep = (r) => Math.max(2 * r * Math.sin(dth / 2), 1.2 * dLine);
    const spacing = Array.from({ length: nS + 1 }, (_, i) => sep(rf + Xs * (1 - i / nS)));      // station 0 = touchdown end, n = fairlead
    const top = b.segs[b.segs.length - 1];
    let cur = null;
    if (inp.Ucross > 0) {
      const def = (T) => ({ L: sus, w: Math.max((T - Math.max(T - top.w * b.hF, 0.05 * T)) / sus, 0), d: dLine, Cd: inp.CdLine, mode: 'top-tensioned', topTension: T });
      try { cur = clearanceUnderCurrent({ upstream: def(lA.T), downstream: def(lB.T), current: inp.Ucross, spacing, wakeModel: inp.wake, n: nS, rho: inp.rhoW }); } catch (e) { cur = null; }
    }
    const curMin = cur ? cur.minClearance : Math.min(...spacing) - dLine;
    const curLoss = cur ? Math.max(0, ...spacing.map((v, i) => v - dLine - cur.clearance[i])) : 0;
    const checks = [
      check('Line-to-line clearance', inp.reqLineLine, Math.max(llMin - dLine, 1e-9), 'm', 'DNV-OS-E301 Ch.2 Sec.2 [5]; API RP 2SK (3rd ed.) §6.6',
        'Required over available surface distance between the closest pair of lines. Contact between lines chafes wire and rope and can foul chain; the smallest gap is normally just below the fairleads.'),
      check('Line-to-asset clearance', inp.reqLineAsset, Math.max(laMin, 1e-9), 'm', 'API RP 2SK (3rd ed.) §6.6 (lines over pipelines and near structures)',
        'Required over available clearance to the nearest pipeline, riser or structure (after its exclusion radius). A line touching a pipeline can abrade its coating or drag it sideways; a line touching a riser can cut it.'),
      check('Anchor stand-off from seabed assets', inp.reqAnchor, Math.max(anchorMin, 1e-9), 'm', 'API RP 2SK (3rd ed.) §6.6 (anchor to pipeline distance)',
        'Required over available horizontal distance from an anchor to a seabed asset. An anchor that drags during setting or in a storm must not reach the pipeline.'),
      { name: 'Lines crossing in plan', demand: nCross, capacity: 0, unit: 'pairs', uc: nCross > 0 ? 99 : 0, ref: 'DNV-OS-E301 Ch.2 Sec.2 [5] (crossing lines to be avoided or shown clear)',
        explain: 'Two lines whose plan projections cross can touch when the vessel moves; crossings need an explicit vertical-separation study in all conditions.' },
      check('Line-to-line clearance in cross-current', inp.reqLineLine, Math.max(curMin, 1e-9), 'm', 'DNV-RP-F203 (2009) §4 (wake-induced interference)',
        'Required over available gap between the closest pair when a cross-current pushes the upstream line into its shielded neighbour, which deflects less.'),
    ];
    const g = worst(checks);
        const lf = run(() => lineForces(sys, poseL, { nProfile: 10 }), 'Line profiles').lines.filter((l) => !l.broken).map((l) => ({ name: l.name, points: l.points }));
    return {
      kpis: [
        { key: 'minLineLine', label: 'Minimum line-to-line distance (centres)', value: round(llMin, 5), unit: 'm', limit: round(inp.reqLineLine + dLine, 5), sense: 'min', explain: 'Smallest centre-to-centre distance between any two lines over both vessel positions. Lines of one cluster are closest at the fairleads.' },
        { key: 'minLineLineCalm', label: 'Line-to-line distance, calm water', value: round(llC.distance, 5), unit: 'm', explain: 'Closest approach of two lines with the vessel at its calm-water position.' },
        { key: 'minLineLineLoaded', label: 'Line-to-line distance, loaded position', value: round(llL.distance, 5), unit: 'm', explain: 'Closest approach of two lines with the vessel at its mean offset under the environmental load.' },
        { key: 'minLineAsset', label: 'Minimum line-to-asset clearance', value: round(laMin, 5), unit: 'm', limit: inp.reqLineAsset, sense: 'min', explain: 'Smallest distance from any line to any asset, less the asset exclusion radius, over both vessel positions (9999 when no assets are defined).' },
        { key: 'minLineAssetCalm', label: 'Line-to-asset clearance, calm water', value: round(assets.length ? laC.clearance : 9999, 5), unit: 'm', explain: 'Closest line-to-asset clearance at the calm-water position.' },
        { key: 'minLineAssetLoaded', label: 'Line-to-asset clearance, loaded position', value: round(assets.length ? laL.clearance : 9999, 5), unit: 'm', explain: 'Closest line-to-asset clearance at the loaded position: leeward lines sag towards the seabed, windward lines lift.' },
        { key: 'anchorStandoff', label: 'Smallest anchor stand-off from seabed assets', value: round(anchorMin, 5), unit: 'm', limit: inp.reqAnchor, sense: 'min', explain: 'Horizontal distance from the nearest anchor to a pipeline or seabed structure (9999 when none is defined).' },
        { key: 'crossings', label: 'Line pairs crossing in plan', value: nCross, unit: '-', limit: 0, sense: 'max', explain: 'Number of line pairs whose anchor-to-fairlead chords intersect when seen from above.' },
        { key: 'currentDeflection', label: 'Line deflection in the cross-current', value: round(cur ? cur.upstream.yMax : 0, 4), unit: 'm', explain: 'Sideways bow of the upstream line of the closest pair under the cross-current, as a tensioned string over its suspended length.' },
        { key: 'currentClearance', label: 'Closest-pair clearance in the cross-current', value: round(curMin, 4), unit: 'm', limit: inp.reqLineLine, sense: 'min', explain: 'Surface gap between the closest pair after the upstream line has deflected towards the wake-shielded downstream line.' },
        { key: 'wakeLoss', label: 'Largest clearance lost to wake shielding', value: round(curLoss, 4), unit: 'm', explain: 'Largest local reduction of the gap along the suspended length, caused by the downstream line sitting in the wake and deflecting less than its neighbour.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Plan view: lines and assets at the loaded position', x: { label: 'x', unit: 'm' }, y: { label: 'y', unit: 'm' },
          series: [...lf.map((l) => ({ name: l.name, x: l.points.map((p) => round(p[0], 5)), y: l.points.map((p) => round(p[1], 5)) })), ...assets.map((a) => ({ name: a.name, x: a.points.map((p) => round(p[0], 5)), y: a.points.map((p) => round(p[1], 5)), dash: true }))],
          explain: 'Lines (solid) and assets (dashed) seen from above with the vessel at its mean offset. Look for lines passing over assets and for anchors close to them.' },
        { type: 'model3d', title: 'Three-dimensional view of lines and assets', nodes: [...lf.flatMap((l) => l.points.map((p) => p.map((v) => round(v, 5)))), ...assets.flatMap((a) => a.points.map((p) => p.map((v) => round(v, 5))))],
          lines: (() => { const out = []; let k = 0; for (const l of lf) { for (let i = 0; i + 1 < l.points.length; i++) out.push([k + i, k + i + 1]); k += l.points.length; } for (const a of assets) { for (let i = 0; i + 1 < a.points.length; i++) out.push([k + i, k + i + 1]); k += a.points.length; } return out; })(),
          explain: 'Spatial arrangement of the catenaries and the subsea assets (z up, still-water level at zero). Rotate to judge vertical separation where a line passes over a pipeline or beside a riser.' },
        { type: 'bar', title: 'Clearance by asset', categories: assets.map((a) => a.name),
          series: [{ name: 'Calm water (m)', y: assets.map((a) => round(Math.min(...cC.lineAsset.filter((r) => r.asset === a.name).map((r) => r.clearance)))) }, { name: 'Loaded position (m)', y: assets.map((a) => round(Math.min(...cL.lineAsset.filter((r) => r.asset === a.name).map((r) => r.clearance)))) }, { name: 'Required (m)', y: assets.map(() => inp.reqLineAsset) }],
          explain: 'Smallest line clearance to each asset at the two vessel positions against the requirement. A bar that shrinks from calm to loaded shows which assets are approached as the vessel offsets.' },
        ...(cur ? [{ type: 'line', title: 'Closest pair in cross-current', x: { label: 'Arc length from touchdown', unit: 'm' }, y: { label: 'Distance', unit: 'm' },
          series: [{ name: 'Clearance with current', x: Array.from(cur.s, (v) => round(v)), y: Array.from(cur.clearance, (v) => round(v)) }, { name: 'Clearance without current', x: Array.from(cur.s, (v) => round(v)), y: spacing.map((v) => round(v - dLine)) },
            { name: 'Upstream deflection', x: Array.from(cur.s, (v) => round(v)), y: Array.from(cur.upstream.y, (v) => round(v)) }, { name: 'Downstream deflection', x: Array.from(cur.s, (v) => round(v)), y: Array.from(cur.downstream.y, (v) => round(v)) }],
          explain: 'Gap between the closest two lines along their suspended length, with and without the cross-current, and the deflection of each. The difference between the two deflection curves is the clearance lost to wake shielding.' }] : []),
      ],
      tables: [
        { title: 'Line-to-asset clearances', columns: [{ key: 'line', label: 'Line' }, { key: 'asset', label: 'Asset' }, { key: 'calm', label: 'Calm water', unit: 'm' }, { key: 'loaded', label: 'Loaded position', unit: 'm' }, { key: 'ok', label: 'Meets requirement' }],
          rows: cC.lineAsset.map((r, i) => ({ line: sys.lines[r.line].name, asset: r.asset, calm: round(r.clearance), loaded: round(cL.lineAsset[i].clearance), ok: Math.min(r.clearance, cL.lineAsset[i].clearance) >= inp.reqLineAsset ? 'yes' : 'NO' })).sort((p, q) => Math.min(p.calm, p.loaded) - Math.min(q.calm, q.loaded)).slice(0, 60),
          explain: 'Every line against every asset, sorted by smallest clearance (exclusion radius deducted). The top rows are the interference cases to report.' },
        { title: 'Closest line pairs', columns: [{ key: 'a', label: 'Line' }, { key: 'b', label: 'Line' }, { key: 'calm', label: 'Calm water', unit: 'm' }, { key: 'loaded', label: 'Loaded position', unit: 'm' }],
          rows: cC.lineLine.map((r, i) => ({ a: sys.lines[r.a].name, b: sys.lines[r.b].name, calm: round(r.distance), loaded: round(cL.lineLine[i].distance) })).sort((p, q) => Math.min(p.calm, p.loaded) - Math.min(q.calm, q.loaded)).slice(0, 20),
          explain: 'Centre-to-centre distance of the closest line pairs at both positions. Pairs within a cluster head the list; their gap is set by fairlead spacing.' },
      ],
      notes: [`Closest lines: ${sys.lines[llGov.a].name} and ${sys.lines[llGov.b].name} at ${fmt(llMin)} m centre to centre (line diameter ${fmt(dLine * 1e3)} mm).`,
        laGov ? `Closest asset: ${laGov.asset}, ${fmt(laMin)} m clear of ${sys.lines[laGov.line].name}.` : 'No assets defined.',
        nCross ? `Crossing pairs: ${crossPairs.slice(0, 6).join('; ')}.` : 'No lines cross in plan.',
        anchorAt ? `Nearest anchor to a seabed asset: ${anchorAt.line} anchor, ${fmt(anchorMin)} m from ${anchorAt.asset}.` : 'No seabed assets near the anchors.'],
      advice: [
        adv('engineering', sev(g.uc), `Smallest line-to-line distance ${fmt(llMin)} m, line-to-asset clearance ${fmt(laMin)} m${laGov ? ` (${laGov.asset})` : ''}, anchor stand-off ${fmt(anchorMin)} m, ${nCross} crossing(s). Governing: ${g.name} (${fmt(g.uc)}).`,
          g.uc > 1 ? (g.name.startsWith('Line-to-asset') ? `Raise the line over ${laGov.asset} with a buoy or a lighter (wire/fibre) section there, shorten the grounded length at that crossing, or re-head the line by a few degrees; ${fmt(inp.reqLineAsset - laMin)} m of clearance is missing.`
            : g.name.startsWith('Anchor') ? `Move the ${anchorAt.line} anchor at least ${fmt(inp.reqAnchor - anchorMin)} m further from ${anchorAt.asset}, or use a non-dragging anchor (suction/pile) there.`
              : g.name.startsWith('Lines crossing') ? 'Re-arrange headings so that no two lines cross in plan, or prove vertical separation in every intact and damaged condition.'
                : `Widen the fairlead spacing or the angle between lines in a group (now ${fmt(dth / DEG)}° for the closest pair) to open the gap by ${fmt(inp.reqLineLine - Math.min(llMin - dLine, curMin))} m.`)
            : 'Clearances are met at both positions; repeat the check at the damaged-condition offset and with the as-laid anchor and pipeline survey positions.',
          'API RP 2SK (3rd ed., 2005) §6.6; DNV-OS-E301 (2018) Ch.2 Sec.2 [5]; DNV-RP-F203 (2009).'),
        adv('economic', 'info', g.uc > 1 ? 'A clearance violation found at design stage costs a re-heading or a buoyancy element; found offshore it costs a re-lay.' : 'No re-routing or protection is needed for the assets listed.',
          g.uc > 1 ? 'Budget for either one buoy and its connectors per affected line, or pipeline protection (mattresses) at the crossing plus about one vessel day to install it, and compare with moving the anchor.' : 'Hold the anchor-placement tolerance (typically ±10–15 m) in the installation contract so the as-built clearances stay within these results; an as-laid survey avoids a later inspection campaign.',
          'Marine-operations planning per DNV-ST-N001 (2021); crossing agreements with third-party asset owners; AACE RP 18R-97 Class 4.',
          `${vesselCo2Text(1)} Avoiding line contact with pipelines and structures also avoids coating damage, spill risk and unplanned intervention on the seabed.`),
      ],
      quantities: { steel_t: round(b.steel_t), vessel_days: { clearance_survey: 1 } },
    };
  },
  verify: [
    { // Spread of 8 lines with line 1 along +x (anchor at x = 1130 m, seabed z = −300 m). A horizontal asset segment at x = 1080 m,
      // 5 m above the seabed, crosses over the grounded part of line 1: distance 5.0 m, less the 0.5 m exclusion radius = 4.5 m.
      name: 'Asset crossing 5 m above a grounded line',
      source: 'Geometry (hand calculation)',
      inputs: { layout: 'spread', headingDeg: 0, Fmean: 0, Ucross: 0, assets: [{ name: 'test pipe', x1: 1080, y1: -100, z1: -295, x2: 1080, y2: 100, z2: -295, radius: 0.5 }] },
      expect: [{ kpi: 'minLineAsset', value: 4.5, tol: 0.002 }] },
    { // Same spread: adjacent fairleads sit on a 40 m radius 45° apart, so the closest two lines are 2·40·sin(22.5°) = 30.61 m apart.
      name: 'Closest approach of adjacent lines at the fairleads',
      source: 'Geometry (chord of the fairlead circle)',
      inputs: { layout: 'spread', headingDeg: 0, Fmean: 0, Ucross: 0, assets: [] },
      expect: [{ kpi: 'minLineLine', value: 30.615, tol: 0.002 }] },
  ],
};

