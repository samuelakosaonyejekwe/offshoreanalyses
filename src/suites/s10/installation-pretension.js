// Mooring installation and pretensioning: line length or anchor position for the target pretension, tolerance
// band, sensitivity of tension to paid-out length, stepwise tensioning programme and proof-load check.
import { setPretension, lineForces } from '../../engine/lines/mooringsystem.js';
import { DEG, need, run, round, fmt, adv, check, worst, sev, systemInputs, loadInputs, buildSystem, fastSystem, elevationPlot, vesselCo2Text, steelCo2Text, LIMIT_QS, LIMIT_PROPS } from './_common.js';

export default {
  id: 's10-installation-pretension',
  suite: 10,
  title: 'Mooring installation and pretension',
  summary: 'Finds the line length or anchor radius that gives the target pretension, how tightly it must be controlled, and the tensioning and proof-loading steps to reach it.',
  covers: [13],
  fidelity: 'analytical',
  standards: ['API RP 2SK (3rd ed., 2005) §6.5 (anchor test load), §8 (installation)', 'DNV-OS-E301 (2018) Ch.2 Sec.4 and Ch.3 (installation, proof tensioning)', 'DNV-ST-N001 (2021) Marine operations and marine warranty', 'ISO 19901-7 (2013) §12–13'],
  models: ['Multi-segment elastic catenary inverse solve for target fairlead tension (Brent)', 'Tension–length sensitivity dT/dL', 'Linearised offset from pretension imbalance, ΔF / K', 'Quasi-static equilibrium at the tolerance limits', 'Stepwise cross-tensioning programme'],
  limits: `${LIMIT_QS} Installation is treated as a sequence of static states: winch and chain-jack dynamics, vessel motion during tensioning, anchor drag during proof loading and seabed friction hysteresis (the line remembers its lay tension) are not modelled. ${LIMIT_PROPS}`,
  publishes: ['mooring.pretension', 'mooring.lineLength'],
  searchTerms: 'mooring pretension installation tensioning proof load tolerance',
  inputs: [
    ...systemInputs(),
    { key: 'targetT', label: 'Target pretension at the fairlead', unit: 'kN', type: 'number', default: 1300, min: 10, max: 30000, group: 'Mooring pretension', role: 'ic', design: true,
      help: 'Initial condition of the installed system: calm-water fairlead tension of every line at the as-designed vessel position (typically 10–20 % of MBL for catenary systems, 15–25 % for taut polyester).' },
    { key: 'tolPct', label: 'Installation tension tolerance', unit: '± %', type: 'number', default: 5, min: 0.5, max: 30, group: 'Mooring pretension', role: 'criteria', help: 'Allowed deviation of each line from the target (±5 % is common; ±10 % for mobile units).' },
    { key: 'mode', label: 'Pretension set by', type: 'select', default: 'length', group: 'Mooring pretension', role: 'bc',
      options: [{ value: 'length', label: 'Line length (winch / chain jack at the fairlead)' }, { value: 'anchorRadius', label: 'Anchor position (fixed line length)' }],
      help: 'Boundary condition used to establish pretension: prescribed line length paid out at the fairlead, or anchors placed at the radius that gives the tension with the as-built length.' },
    { key: 'adjustSeg', label: 'Segment adjusted (0 = top segment)', type: 'int', default: 0, min: 0, max: 12, group: 'Mooring pretension', role: 'control', showIf: { key: 'mode', in: ['length'] }, help: 'Number of the segment in the make-up table whose length is changed; normally the top chain pulled through the fairlead.' },
    { key: 'layPct', label: 'As-laid tension before tensioning', unit: '% of target', type: 'number', default: 50, min: 2, max: 100, group: 'Mooring pretension', role: 'ic', help: 'Tension left in the line after hook-up, before final tensioning; sets how much line must be hauled in.' },
    { key: 'nStages', label: 'Tensioning stages', type: 'int', default: 4, min: 1, max: 10, group: 'Tensioning procedure', role: 'operation', help: 'Number of equal tension increments applied alternately to opposite lines so the vessel stays centred.' },
    { key: 'proofT', label: 'Anchor proof tension at the fairlead', unit: 'kN', type: 'number', default: 3100, min: 0, max: 60000, group: 'Tensioning procedure', role: 'operation', bind: 'mooring.Tmax',
      help: 'Test load held (typically 15–30 minutes) to prove the anchor: API RP 2SK asks for the maximum intact design tension for permanent systems unless the anchor type justifies less.' },
    { key: 'winchPull', label: 'Tensioning capacity at the fairlead', unit: 'kN', type: 'number', default: 4000, min: 10, max: 60000, group: 'Tensioning procedure', role: 'operation', help: 'Stall pull of the winch or chain jack, or the tension obtainable by cross-tensioning with an installation vessel.' },
    { key: 'hoursPerStage', label: 'Time per line per stage', unit: 'h', type: 'number', default: 1.5, min: 0.1, max: 24, group: 'Tensioning procedure', role: 'operation' },
    ...loadInputs(),
  ],
  compute(inp) {
    const b = buildSystem(inp);
    const sys = b.system, n = b.nLines;
    const target = inp.targetT * 1e3;
    need(target < 0.6 * b.MBL, `The target pretension ${inp.targetT} kN is above 60 % of the weakest segment's break load (${fmt(b.MBL / 1e3)} kN): choose a lower pretension.`);
    const iSeg = inp.adjustSeg > 0 ? inp.adjustSeg - 1 : b.segs.length - 1;
    need(iSeg < b.segs.length, `Segment ${inp.adjustSeg} does not exist: the make-up has ${b.segs.length} segments.`);
    // identical lines (generated layout): solve one line and reuse; imported anchors: solve all
    const identical = !inp.layoutGeom;
    const sub = identical ? { ...sys, lines: [sys.lines[0]] } : sys;
    const setAt = (T) => run(() => setPretension(sub, { target: T, mode: inp.mode, segment: inp.mode === 'length' ? iSeg : undefined }),
      `Pretension solve (${inp.mode === 'length' ? 'the segment length' : 'the anchor radius'} cannot produce ${fmt(T / 1e3)} kN: check the make-up, depth and target)`);
    const set = setAt(target);
    const tol = inp.tolPct / 100;
    const hi = setAt(target * (1 + tol)), lo = setAt(target * (1 - tol));
    // the as-laid tension cannot be below the weight of line hanging at the fairlead: step up until a state exists
    let lay = null, layFrac = inp.layPct / 100;
    for (const f of [inp.layPct / 100, 0.5, 0.65, 0.8, 0.9, 0.97].filter((v) => v >= inp.layPct / 100)) {
      try { lay = setPretension(sub, { target: target * f, mode: inp.mode, segment: inp.mode === 'length' ? iSeg : undefined }); layFrac = f; break; } catch (e) { lay = null; }
    }
    if (!lay) { lay = set; layFrac = 1; }
    const adj = set.adjustments[0], adjHi = hi.adjustments[0], adjLo = lo.adjustments[0], adjLay = lay.adjustments[0];
    // d(tension)/d(length or span) from the tolerance band
    const dTdA = (2 * tol * target) / Math.abs(adjHi - adjLo);
    const band = Math.abs(adjHi - adjLo) / 2;
    const Lseg0 = b.segs[iSeg].L, span0 = inp.anchorRadius - inp.fairleadRadius;
    const Ltot = inp.mode === 'length' ? b.Ltot - Lseg0 + adj : b.Ltot;
    const radius = inp.mode === 'length' ? inp.anchorRadius : adj + inp.fairleadRadius;
    // full adjusted system for equilibrium checks
    const full = (s) => (identical ? { ...sys, lines: sys.lines.map((ln) => {
      const src = s.system.lines[0];
      if (inp.mode === 'length') return { ...ln, segments: src.segments.map((g) => ({ ...g })) };
      const az = Math.atan2(ln.anchor[1] - ln.fairlead[1], ln.anchor[0] - ln.fairlead[0]), sp = s.adjustments[0];
      return { ...ln, anchor: [ln.fairlead[0] + sp * Math.cos(az), ln.fairlead[1] + sp * Math.sin(az), ln.anchor[2]] };
    }) } : s.system);
    const F = inp.Fmean * 1e3, a = inp.loadDirDeg * DEG, Fext = [F * Math.cos(a), F * Math.sin(a), 0, 0, 0, 0];
    const evalSys = (s) => {
      const fs = fastSystem(full(s), { maxOffset: Math.max(0.15 * b.depth, 30) });
      const calm = fs.solve([0, 0, 0, 0, 0, 0]);
      const eq = fs.solve(Fext, [calm.x, calm.y]);
      const f1 = fs.forces(calm.x + 0.05, calm.y), f0 = fs.forces(calm.x - 0.05, calm.y);
      return { fs, calm, eq, K: -(f1.F[0] - f0.F[0]) / 0.1 };
    };
    const mid = evalSys(set), up = evalSys(hi), dn = evalSys(lo);
    // worst pretension imbalance: lines on one side at +tol, the other side at −tol (linearised)
    const dF = tol * mid.calm.lines.reduce((s, l, i) => s + l.H * Math.abs(Math.cos(Math.atan2(sys.lines[i].anchor[1], sys.lines[i].anchor[0]))), 0);
    const imbalance = dF / mid.K;
    // proof-load state of one line (vessel held on station by the opposite line)
    const proof = inp.proofT * 1e3;
    let proofSol = null;
    if (proof > 0) { try { proofSol = setPretension({ ...sub, lines: [full(set).lines[0]] }, { target: proof, mode: 'length', segment: iSeg }); } catch (e) { proofSol = null; } }
    let proofLine = null;
    if (proofSol) proofLine = lineForces(proofSol.system, [0, 0, 0, 0, 0, 0]).lines[0];
    const pitch = b.segs[iSeg].family === 'chain' ? 4 * b.segs[iSeg].d : 0.05;     // one chain link advances 4d; winches are continuous
    const checks = [
      check('Tensioning capacity for the proof load', Math.max(proof, target * (1 + tol)) / 1e3, inp.winchPull, 'kN', 'API RP 2SK (3rd ed.) §6.5 (anchor test load); DNV-OS-E301 Ch.2 Sec.4',
        'Guards against an installation spread that cannot reach the proof tension: the anchors would then enter service unproven and the pretension could not be set at the top of its tolerance band.'),
      check('Length resolution against the tolerance band', pitch / 2, Math.max(band, 1e-9), 'm', 'Installation practice: adjustment increment (one chain link = 4 d at a chain stopper) vs allowable length error',
        'Guards against a pretension tolerance tighter than the hardware can deliver: locking the chain one link early or late must keep the tension inside the band.'),
      { name: 'Anchor uplift at proof tension', demand: round(proofLine ? proofLine.upliftAngle / DEG : 0, 4), capacity: 0, unit: 'deg', uc: proofLine && proofLine.upliftAngle > 1e-6 ? 99 : 0, ref: 'API RP 2SK (3rd ed.) §6.5',
        explain: 'During proof loading the line must still pull horizontally on a drag anchor; uplift at the test load would pull the anchor out instead of setting it.' },
      check('Maximum tension with pretension at the upper tolerance', up.eq.maxTension / 1e3, b.MBL / 2.0 / 1e3, 'kN', 'API RP 2SK (3rd ed.) §6.3 quasi-static factor 2.0 on mean tension (screen)',
        'Guards against the as-installed scatter eating into strength: with every line at the top of the band the mean storm tension is highest.'),
    ];
    const g = worst(checks);
    const stages = [];
    for (let s = 0; s <= inp.nStages; s++) {
      const T = target * (layFrac + (1 - layFrac) * s / inp.nStages);
      const st = s === 0 ? lay : s === inp.nStages ? set : setAt(T);
      stages.push({ stage: s === 0 ? 'As laid' : `Stage ${s}`, T, adj: st.adjustments[0] });
    }
    if (proof > 0 && proofSol) stages.push({ stage: 'Proof load (hold)', T: proof, adj: proofSol.adjustments[0] }, { stage: 'Lock-off at target', T: target, adj });
    const haul = Math.abs(adjLay - adj);
    const days = n * (inp.nStages + (proof > 0 ? 2 : 0)) * inp.hoursPerStage / 24 + 1;
    const unit = inp.mode === 'length' ? 'segment length' : 'anchor–fairlead span';
    return {
      kpis: [
        { key: 'pretension', label: 'Pretension achieved', value: round(mid.calm.lines[0].T / 1e3, 5), unit: 'kN', explain: 'Calm-water fairlead tension of the adjusted system; it should equal the target.' },
        { key: 'lineLength', label: 'Line length for the target pretension', value: round(Ltot, 6), unit: 'm', explain: 'Total unstretched length between anchor and fairlead that gives the target pretension at the design vessel position.' },
        { key: 'segLength', label: `Length of adjusted segment (${b.segs[iSeg].name})`, value: round(inp.mode === 'length' ? adj : Lseg0, 6), unit: 'm', explain: 'Length of the segment used for adjustment after tensioning; the difference from the supplied length is what stays in the chain locker or is cut.' },
        { key: 'anchorRadiusSet', label: 'Anchor radius for the target pretension', value: round(radius, 6), unit: 'm', explain: 'Horizontal distance from the mooring centre to the anchors that gives the target pretension with the entered line length.' },
        { key: 'dTdL', label: 'Tension sensitivity to adjustment', value: round(dTdA / 1e3, 4), unit: 'kN/m', explain: `Change of pretension per metre of ${unit}. A high value means a stiff, taut system where small length errors cause large tension errors.` },
        { key: 'lengthTol', label: 'Adjustment tolerance for the tension band', value: round(band, 4), unit: 'm', explain: `Allowed error in ${unit} that keeps the pretension within ± the tolerance. Compare with one chain link (4 d) and with anchor placement accuracy.` },
        { key: 'haul', label: 'Line to haul in from the as-laid state', value: round(haul, 4), unit: 'm', explain: 'Length of line to recover (or distance to move the anchor) between the as-laid tension and the target.' },
        { key: 'Kxx', label: 'System stiffness at the target pretension', value: round(mid.K / 1e3, 5), unit: 'kN/m', explain: 'Surge stiffness of the installed system; pretension is the main control on it.' },
        { key: 'imbalance', label: 'Offset from worst pretension imbalance', value: round(imbalance, 4), unit: 'm', explain: 'Calm-water offset if all lines on one side sit at the upper tolerance and those opposite at the lower (linearised ΔF/K). It consumes part of the riser offset budget.' },
        { key: 'offsetHi', label: 'Mean offset, pretension at +tolerance', value: round(up.eq.offset, 4), unit: 'm', explain: 'Offset under the mean load when every line is at the top of the band (stiffest system).' },
        { key: 'offsetLo', label: 'Mean offset, pretension at −tolerance', value: round(dn.eq.offset, 4), unit: 'm', explain: 'Offset under the mean load when every line is at the bottom of the band (softest system); the governing case for excursion.' },
        { key: 'TmaxHi', label: 'Largest mean tension, pretension at +tolerance', value: round(up.eq.maxTension / 1e3, 5), unit: 'kN', explain: 'Mean storm tension with every line at the top of the band; the governing case for strength.' },
        { key: 'proofGrounded', label: 'Grounded length at proof tension', value: round(proofLine ? proofLine.groundedLength : 0, 4), unit: 'm', explain: 'Line still on the seabed while the proof load is held. It must stay positive for a drag anchor.' },
        { key: 'days', label: 'Tensioning campaign duration', value: round(days, 4), unit: 'days', explain: 'Vessel time for staged tensioning and proof loading of all lines at the entered rate, plus one day of set-up.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Pretension vs adjustment', x: { label: inp.mode === 'length' ? 'Adjusted segment length' : 'Anchor–fairlead span', unit: 'm' }, y: { label: 'Fairlead tension', unit: 'kN' },
          series: [{ name: 'Tension', x: stages.slice(0, inp.nStages + 1).map((s) => round(s.adj, 6)), y: stages.slice(0, inp.nStages + 1).map((s) => round(s.T / 1e3)) },
            { name: 'Upper tolerance', x: [Math.min(adjLay, adjHi), Math.max(adjLay, adjHi)], y: [target * (1 + tol) / 1e3, target * (1 + tol) / 1e3], dash: true }, { name: 'Lower tolerance', x: [Math.min(adjLay, adjHi), Math.max(adjLay, adjHi)], y: [target * (1 - tol) / 1e3, target * (1 - tol) / 1e3], dash: true }],
          explain: 'Tension build-up as the line is hauled in (or the anchor moved out) through the tensioning stages. The curve steepens near the target: the last metres need the finest control.' },
        elevationPlot(full(set), [0, 0, 0, 0, 0, 0], [0], 'Line profile at the target pretension'),
        { type: 'bar', title: 'Effect of the tolerance band under the mean load', categories: ['−tolerance', 'target', '+tolerance'],
          series: [{ name: 'Mean offset (m × 100)', y: [dn, mid, up].map((e) => round(e.eq.offset * 100)) }, { name: 'Largest mean tension (kN)', y: [dn, mid, up].map((e) => round(e.eq.maxTension / 1e3)) }],
          explain: 'How installation scatter moves the design results: low pretension gives more offset, high pretension more tension. Offsets are scaled by 100 to share the axis.' },
      ],
      tables: [
        { title: 'Tensioning programme (per line; alternate opposite lines at each stage)', columns: [{ key: 'stage', label: 'Step' }, { key: 'T', label: 'Fairlead tension', unit: 'kN' }, { key: 'adj', label: inp.mode === 'length' ? 'Adjusted segment length' : 'Anchor–fairlead span', unit: 'm' }, { key: 'd', label: 'Change from previous step', unit: 'm' }],
          rows: stages.map((s, i) => ({ stage: s.stage, T: round(s.T / 1e3), adj: round(s.adj, 6), d: i ? round(s.adj - stages[i - 1].adj, 4) : 0 })),
          explain: 'Target tension and corresponding length (or span) at each step. Tension opposite lines in pairs through each stage to keep the vessel centred; hold the proof load before locking off at the target.' },
        { title: 'Installed line results at the target pretension', columns: [{ key: 'line', label: 'Line' }, { key: 'T', label: 'Pretension', unit: 'kN' }, { key: 'lo', label: 'Lower limit', unit: 'kN' }, { key: 'hi', label: 'Upper limit', unit: 'kN' }, { key: 'aH', label: 'Anchor load', unit: 'kN' }, { key: 'gr', label: 'Grounded length', unit: 'm' }],
          rows: mid.calm.lines.map((l) => ({ line: l.name, T: round(l.T / 1e3), lo: round(target * (1 - tol) / 1e3), hi: round(target * (1 + tol) / 1e3), aH: round(l.anchorH / 1e3), gr: round(l.groundedLength) })),
          explain: 'Individual line pretensions with their acceptance band, and the equilibrium anchor load and grounded length in calm water: the as-installed record to sign off against.' },
      ],
      notes: [`Target ${inp.targetT} kN reached with ${unit} ${fmt(adj, 6)} m; ±${inp.tolPct} % corresponds to ±${fmt(band)} m (${fmt(dTdA / 1e3)} kN per metre).`,
        `From the as-laid state at ${fmt(100 * layFrac)} % of target${layFrac > inp.layPct / 100 + 1e-9 ? ' (raised from the entered value, which is below the weight of line hanging at the fairlead)' : ''}, ${fmt(haul)} m must be ${inp.mode === 'length' ? 'hauled in' : 'gained by moving the anchor'}.`,
        proofLine ? `At the ${inp.proofT} kN proof load ${fmt(proofLine.groundedLength)} m stays grounded and the anchor sees ${fmt(proofLine.anchorH / 1e3)} kN.` : 'No proof-load state (zero proof tension or not reachable with this make-up).'],
      advice: [
        adv('engineering', sev(g.uc), `Target pretension needs ${unit} ${fmt(adj, 6)} m with a tolerance of ±${fmt(band)} m; worst imbalance shifts the vessel ${fmt(imbalance)} m. Governing: ${g.name} (${fmt(g.uc)}).`,
          g.uc > 1 ? (g.name.startsWith('Tensioning') ? `Provide ${fmt(Math.max(proof, target * (1 + tol)) / 1e3)} kN at the fairlead: cross-tension opposite lines with the platform winches, use a subsea tensioner, or agree a reduced proof load with the certifier for the anchor type.`
            : g.name.startsWith('Length') ? `Relax the tolerance to about ±${fmt(100 * (pitch / 2) * dTdA / target)} % or use a continuous tensioner; one chain link (${fmt(pitch)} m) changes tension by ${fmt(pitch * dTdA / 1e3)} kN.`
              : g.name.startsWith('Anchor') ? 'Add ground chain or reduce the proof load so the line stays tangent to the seabed at the anchor during the test.' : 'Lower the target pretension or strengthen the line: the upper tolerance case over-uses the break load.')
            : `Record achieved tension and paid-out length per line; re-check pretension after the first storm, since seabed friction and anchor settlement typically relax it by a few per cent.`,
          'API RP 2SK (3rd ed., 2005) §6.5 and §8; DNV-OS-E301 (2018) Ch.2 Sec.4; DNV-ST-N001 (2021) §2.'),
        adv('economic', 'info', `Staged tensioning and proof loading of ${n} lines takes about ${fmt(days)} vessel days at ${inp.hoursPerStage} h per line per stage.`,
          inp.nStages > 3 ? `Dropping from ${inp.nStages} to 3 stages would save about ${fmt(n * (inp.nStages - 3) * inp.hoursPerStage / 24)} days; acceptable where the vessel can tolerate the temporary offset of ${fmt(target / inp.nStages / mid.K * 1.5)} m per unbalanced stage.` : 'The programme is already short; savings would come from tensioning opposite lines simultaneously.',
          'Marine-operations planning with weather-restricted durations per DNV-ST-N001 (2021) §2.6; AACE RP 18R-97 Class 4 quantities.',
          `${vesselCo2Text(days)} A tighter pretension tolerance also avoids over-sizing lines for installation scatter: ${steelCo2Text(b.steel_t * tol)}`),
      ],
      publish: { 'mooring.pretension': mid.calm.lines[0].T, 'mooring.lineLength': Ltot },
      quantities: { steel_t: round(b.steel_t * Ltot / b.Ltot), vessel_days: { mooring_tensioning: round(days, 3) } },
    };
  },
  verify: [
    { // Inextensible chain w = 1500 N/m, fairlead 280 m above the seabed, span 910.986 m, target T = 1220 kN → H = T − w·h = 800 kN,
      // a = 533.33 m, suspended s = sqrt(h² + 2ha) = 614.06 m, x = a·acosh(1 + h/a) = 525.04 m → L = s + (X − x) = 1000.0 m.
      name: 'Line length for a target pretension (inextensible catenary)',
      source: 'Faltinsen (1990) Sea Loads, Eqs. 8.14–8.20 inverted for the line length',
      inputs: { mode: 'length', targetT: 1220, depth: 300, fairleadDraft: 20, fairleadRadius: 40, anchorRadius: 950.986, CB: 0, proofT: 0, Fmean: 0,
        segments: [{ name: 'chain', type: 'custom', d: 0, grade: '', L: 1050, w: 1500, EA: 1e7, MBL: 9000 }] },
      expect: [{ kpi: 'lineLength', value: 1000, tol: 0.0005 }, { kpi: 'pretension', value: 1220, tol: 0.002 }] },
    { // Same line with L = 1000 m fixed: the span for T = 1220 kN is 910.986 m, i.e. an anchor radius of 950.986 m.
      name: 'Anchor radius for a target pretension (inextensible catenary)',
      source: 'Faltinsen (1990) Sea Loads, Eqs. 8.14–8.20',
      inputs: { mode: 'anchorRadius', targetT: 1220, depth: 300, fairleadDraft: 20, fairleadRadius: 40, anchorRadius: 930, CB: 0, proofT: 0, Fmean: 0,
        segments: [{ name: 'chain', type: 'custom', d: 0, grade: '', L: 1000, w: 1500, EA: 1e7, MBL: 9000 }] },
      expect: [{ kpi: 'anchorRadiusSet', value: 950.986, tol: 0.0005 }] },
  ],
};
