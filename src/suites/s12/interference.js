// Mechanical interference and clearance: deflection of the umbilical and a neighbouring riser or mooring line in a
// steady current with wake shielding of the downstream line, clearance along the length and contact screening.
import { clearanceUnderCurrent, impactVelocityEstimate } from '../../engine/lines/interference.js';
import { adv, check, worst, sev, fmt, round, run, need, linspace, thin, G, CO2_POLYMER, STD_UMB, vesselCo2Text } from './_common.js';

export default {
  id: 's12-interference-clearance', suite: 12, title: 'Mechanical interference and clearance',
  summary: 'Clearance between the umbilical and a neighbouring line under current, where they come closest and whether they can touch.',
  covers: [14], fidelity: 'reduced-order',
  standards: ['DNV-RP-F203 (2009) §3–5 (riser interference: wake models, clearance and collision assessment)', 'Huse (1993) OTC 7237 (wake velocity deficit)', 'Blevins (2005) JOMAE 127 (wake drag)', `${STD_UMB} §7.3 (interference with adjacent structures)`],
  models: ["Morison drag on a tensioned string (Sparks 2007 ch. 6–7)", 'Huse wake velocity deficit and Blevins wake drag (shielding)', 'Fixed-point coupling of two deflected lines', 'Work–energy impact velocity bound'],
  limits: 'Both lines are near-vertical tensioned strings of equal length sharing end elevations (valid for angles below about 0.2 rad): use it for the upper catenary next to risers or moorings, not for the sag bend or touchdown. In-line arrangement with a constant transverse offset; no wake-induced lift or oscillation, no wave-frequency relative motion and no 3-D catenary geometry — a time-domain interference analysis is needed where the static clearance is small.',
  publishes: ['umbilical.minClearance'],
  searchTerms: 'riser umbilical interference clearance wake shielding current',
  inputs: [
    { key: 'L', label: 'Assessed length (common vertical extent)', unit: 'm', type: 'number', default: 250, min: 10, max: 3000, group: 'Geometry', role: 'geometry', help: 'Near-vertical part over which the two lines run side by side, measured from its lower end.' },
    { key: 'gapTop', label: 'Centre spacing at the top', unit: 'm', type: 'number', default: 5, min: 0.2, max: 200, group: 'Geometry', role: 'bc', design: true, help: 'Hang-off spacing in the current direction (upper termination positions).' },
    { key: 'gapBot', label: 'Centre spacing at the lower end', unit: 'm', type: 'number', default: 12, min: 0.2, max: 500, group: 'Geometry', role: 'bc', design: true },
    { key: 'yOff', label: 'Transverse offset from the wake centreline', unit: 'm', type: 'number', default: 0, min: 0, max: 100, group: 'Geometry', role: 'bc' },
    { key: 'order', label: 'Current direction', type: 'select', default: 'umbilical-up', group: 'Geometry', role: 'bc', options: [{ value: 'umbilical-down', label: 'Neighbour upstream, umbilical in its wake' }, { value: 'umbilical-up', label: 'Umbilical upstream, neighbour in its wake' }] },
    { key: 'OD', label: 'Umbilical diameter', unit: 'mm', type: 'number', default: 145, min: 30, max: 400, group: 'Umbilical', role: 'geometry', bind: 'umbilical.OD' },
    { key: 'w', label: 'Umbilical submerged weight', unit: 'N/m', type: 'number', default: 235, min: 0, max: 3000, group: 'Umbilical', role: 'geometry', bind: 'umbilical.wSub' },
    { key: 'T', label: 'Umbilical tension at the top of the assessed length', unit: 'kN', type: 'number', default: 80, min: 1, max: 5000, group: 'Umbilical', role: 'bc', bind: 'umbilical.topTension', uncertain: { dist: 'lognormal', cov: 0.1 } },
    { key: 'Cd', label: 'Umbilical drag coefficient', unit: '-', type: 'number', default: 1.2, min: 0.3, max: 3.5, group: 'Umbilical', role: 'load', bind: 'umbilical.cdAmplified', uncertain: { dist: 'normal', cov: 0.15 }, help: 'Include the VIV drag amplification from the VIV module where lock-in is predicted.' },
    { key: 'mass', label: 'Umbilical mass in air', unit: 'kg/m', type: 'number', default: 40.9, min: 1, max: 300, group: 'Umbilical', role: 'geometry', bind: 'umbilical.massPerM' },
    { key: 'nD', label: 'Neighbour diameter', unit: 'mm', type: 'number', default: 320, min: 20, max: 2000, group: 'Neighbouring line', role: 'geometry', help: 'Flexible riser, steel riser or mooring line next to the umbilical.' },
    { key: 'nW', label: 'Neighbour submerged weight', unit: 'N/m', type: 'number', default: 900, min: 0, max: 20000, group: 'Neighbouring line', role: 'geometry' },
    { key: 'nT', label: 'Neighbour tension at the top of the assessed length', unit: 'kN', type: 'number', default: 500, min: 1, max: 50000, group: 'Neighbouring line', role: 'bc', bind: 'riser.topTension' },
    { key: 'nCd', label: 'Neighbour drag coefficient', unit: '-', type: 'number', default: 1.1, min: 0.3, max: 3.5, group: 'Neighbouring line', role: 'load' },
    { key: 'Us', label: 'Current speed at the top', unit: 'm/s', type: 'number', default: 1.0, min: 0, max: 3, group: 'Current', role: 'bc', bind: 'site.Uc100', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'Ul', label: 'Current speed at the lower end', unit: 'm/s', type: 'number', default: 0.5, min: 0, max: 3, group: 'Current', role: 'bc' },
    { key: 'curExp', label: 'Profile exponent', unit: '-', type: 'number', default: 1, min: 0.1, max: 8, group: 'Current', role: 'bc' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Current', role: 'environment', bind: 'site.rhoW' },
    { key: 'wake', label: 'Wake model', type: 'select', default: 'huse', group: 'Current', role: 'control', options: [{ value: 'huse', label: 'Huse (1993) velocity deficit' }, { value: 'blevins', label: 'Blevins (2005) wake drag' }, { value: 'none', label: 'No shielding' }] },
    { key: 'minClear', label: 'Required surface-to-surface clearance', unit: 'm', type: 'number', default: 1.0, min: 0, max: 50, group: 'Acceptance criteria', role: 'criteria', help: 'Project criterion; DNV-RP-F203 accepts either no collision (with margin for dynamics and VIV amplitudes) or a collision-damage assessment.' },
    { key: 'n', label: 'Intervals along the length', type: 'int', default: 100, min: 20, max: 2000, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'n', label: 'Intervals', levels: [50, 100, 200, 400], qoi: 'minClearance' },
  compute(inp) {
    const n = inp.n, L = inp.L, umbUp = inp.order === 'umbilical-up';
    const umb = { L, w: inp.w, d: inp.OD / 1e3, Cd: inp.Cd, mode: 'top-tensioned', topTension: inp.T * 1e3 };
    const nb = { L, w: inp.nW, d: inp.nD / 1e3, Cd: inp.nCd, mode: 'top-tensioned', topTension: inp.nT * 1e3 };
    for (const [name, l] of [['umbilical', umb], ['neighbouring line', nb]]) need(l.topTension - l.w * L > 0, `The ${name} would be in compression at the lower end of the assessed length (top tension ${fmt(l.topTension / 1e3)} kN is less than its weight ${fmt(l.w * L / 1e3)} kN): shorten the length or raise the tension.`);
    const sG = linspace(0, L, 41), current = { s: sG, speed: sG.map((s) => inp.Ul + (inp.Us - inp.Ul) * (s / L) ** inp.curExp) };
    const spacing = Float64Array.from({ length: n + 1 }, (_, i) => inp.gapBot + (inp.gapTop - inp.gapBot) * i / n);
    const up = umbUp ? umb : nb, down = umbUp ? nb : umb;
    const solve = (wakeModel) => run(() => clearanceUnderCurrent({ upstream: up, downstream: down, current, spacing, transverseOffset: inp.yOff, wakeModel, n, rho: inp.rhoW }), 'Clearance under current');
    const r = solve(inp.wake), rNo = inp.wake === 'none' ? r : solve('none');
    need(r.converged !== false, 'The wake-shielding iteration did not converge. The lines probably cross: increase the spacing or the downstream tension.');
    const rSum = (umb.d + nb.d) / 2;
    let static0 = Infinity; for (let i = 0; i <= n; i++) static0 = Math.min(static0, Math.hypot(spacing[i], inp.yOff) - rSum);
    const minC = r.minClearance, shieldMin = Math.min(...r.dragFactor), angMax = Math.max(r.upstream.maxAngle, r.downstream.maxAngle);
    // if the gap can close: velocity bound from the net lateral load working across the still-water gap
    const iMin = Math.max(0, Math.min(n, Math.round(r.sAtMin / L * n)));
    const dF = Math.abs(r.upstream.f[iMin] - r.downstream.f[iMin]);
    // the umbilical is the lighter body: its mass plus added mass (Ca = 1) is accelerated across the gap
    const imp = impactVelocityEstimate({ force: dF, gap: Math.max(static0, 0.01), mass: inp.mass + inp.rhoW * Math.PI * umb.d ** 2 / 4 });
    const checks = [
      check('Minimum clearance under current', inp.minClear, Math.max(minC, 1e-6), 'm', 'DNV-RP-F203 (2009) §4 (no-collision criterion with project margin)', 'Demand is the required surface-to-surface clearance, capacity the smallest clearance along the length in the design current. Below it the lines can touch once wave-frequency motion and VIV amplitudes are added; contact abrades the sheath and can snag buoyancy or clamps.'),
      check('Small-angle validity of the string model', angMax, 0.2, 'rad', 'Sparks (2007) ch. 6 (tensioned-string deflection)', 'Largest inclination of either line from vertical against the validity limit of the model. Beyond it the deflections are under-predicted and a catenary / rod interference analysis is required.'),
    ];
    const gov = worst(checks);
    const [sp, cl, cl0, yU, yD, sh] = thin([r.s, r.clearance, rNo.clearance, r.upstream.y, r.downstream.y, r.dragFactor], 200);
    const nodes = [], lines = [], vals = [];
    for (const [yArr, base] of [[r.upstream.y, 0], [r.downstream.y, 1]]) {
      const [ss, yy, cc, gg] = thin([r.s, yArr, r.clearance, spacing], 80); const b0 = nodes.length;
      for (let i = 0; i < ss.length; i++) { nodes.push([yy[i] + base * gg[i], base * inp.yOff, ss[i]]); vals.push(cc[i]); if (i) lines.push([b0 + i - 1, b0 + i]); }
    }
    const collide = r.collision || minC <= 0;
    const gapNeed = Math.max(0, inp.minClear - minC);
    return {
      kpis: [
        { key: 'minClearance', label: 'Minimum clearance under current', value: round(minC, 6), unit: 'm', limit: inp.minClear, sense: 'min', explain: 'Smallest surface-to-surface distance between the two lines in the design current with wake shielding. Negative means the static shapes already cross.' },
        { key: 'sAtMin', label: 'Height of the closest point above the lower end', value: round(r.sAtMin), unit: 'm', explain: 'Potential contact location: where spacers, wear protection or a clash-guard would be needed.' },
        { key: 'clearStill', label: 'Minimum clearance in still water', value: round(static0), unit: 'm', explain: 'Smallest clearance of the undeflected lines; the difference to the value under current is what the current consumes.' },
        { key: 'clearLost', label: 'Clearance consumed by the current', value: round(static0 - minC), unit: 'm', explain: 'Reduction of the minimum clearance caused by unequal deflection of the two lines. The lighter, lower-tensioned line normally deflects more.' },
        { key: 'clearNoWake', label: 'Minimum clearance without wake shielding', value: round(rNo.minClearance), unit: 'm', explain: 'Same case with both lines in the free stream. Shielding lowers the drag on the downstream line and lets the upstream line catch it up.' },
        { key: 'deflUp', label: 'Largest deflection, upstream line', value: round(r.upstream.yMax), unit: 'm', explain: `Downstream deflection of the upstream line (${umbUp ? 'umbilical' : 'neighbour'}) in the free stream; for uniform load and tension it is f·L²/(8T).` },
        { key: 'deflDown', label: 'Largest deflection, downstream line', value: round(r.downstream.yMax), unit: 'm', explain: `Deflection of the downstream line (${umbUp ? 'neighbour' : 'umbilical'}) with its drag reduced by the wake.` },
        { key: 'shield', label: 'Lowest shielding factor on downstream drag', value: round(shieldMin, 3), unit: '-', explain: 'Ratio of drag in the wake to free-stream drag at the most shielded point; 1 means no shielding.' },
        { key: 'angleMax', label: 'Largest inclination from vertical', value: round(angMax, 3), unit: 'rad', limit: 0.2, sense: 'max', explain: 'Validity indicator of the string model; above about 0.2 rad use a catenary model.' },
        { key: 'impactV', label: 'Contact velocity bound', value: round(collide ? imp.velocity : 0, 3), unit: 'm/s', explain: 'Upper bound of the closing velocity if the net current load accelerated the umbilical across the still-water gap (work–energy). Zero when no contact is predicted.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Clearance along the length', x: { label: 'Height above the lower end', unit: 'm' }, y: { label: 'Surface-to-surface clearance', unit: 'm' },
          series: [{ name: `With wake (${inp.wake})`, x: sp, y: cl }, { name: 'No shielding', x: sp, y: cl0, dash: true }, { name: 'Required', x: [0, L], y: [inp.minClear, inp.minClear], dash: true }],
          explain: 'Clearance at every height. The minimum of the solid curve and its height locate the potential contact; a curve dipping below the required line fails the criterion and below zero means the lines cross.' },
        { type: 'line', title: 'Deflected shapes', x: { label: 'Deflection in the current direction', unit: 'm' }, y: { label: 'Height above the lower end', unit: 'm' },
          series: [{ name: `Upstream (${umbUp ? 'umbilical' : 'neighbour'})`, x: yU, y: sp }, { name: `Downstream (${umbUp ? 'neighbour' : 'umbilical'}) incl. spacing`, x: yD.map((v, i) => v + (inp.gapBot + (inp.gapTop - inp.gapBot) * sp[i] / L)), y: sp }],
          explain: 'Both lines in the plane of the current, the downstream one drawn at its spacing. Where the two curves approach each other is the clearance pinch point.' },
        { type: 'line', title: 'Wake shielding of the downstream line', x: { label: 'Height above the lower end', unit: 'm' }, y: { label: 'Drag factor', unit: '-' }, series: [{ name: 'Drag in wake / free-stream drag', x: sp, y: sh }],
          explain: 'How much of its free-stream drag the downstream line keeps at each height. Low values where the lines are close show strong shielding, which drives the lines together.' },
        { type: 'model3d', title: 'Lines in the current (deflected)', nodes, lines, tris: [], nodeValues: vals, valueLabel: 'Clearance', valueUnit: 'm',
          explain: 'Three-dimensional view of the two deflected lines over the assessed length (z up from the lower end), coloured by local clearance. The darkest region is the potential contact location.' },
      ],
      tables: [{ title: 'Clearance at stations', explain: 'Clearance, shielding factor and both deflections at ten heights. Use it to place spacers or wear protection around the minimum.',
        columns: [{ key: 's', label: 'Height', unit: 'm' }, { key: 'gap', label: 'Still-water spacing', unit: 'm' }, { key: 'yU', label: 'Upstream deflection', unit: 'm' }, { key: 'yD', label: 'Downstream deflection', unit: 'm' }, { key: 'sh', label: 'Shielding', unit: '-' }, { key: 'cl', label: 'Clearance', unit: 'm' }],
        rows: Array.from({ length: 11 }, (_, k) => { const i = Math.round(k * n / 10); return { s: round(r.s[i]), gap: round(spacing[i]), yU: round(r.upstream.y[i], 3), yD: round(r.downstream.y[i], 3), sh: round(r.dragFactor[i], 3), cl: round(r.clearance[i], 3) }; }) }],
      notes: [
        `Minimum clearance ${fmt(minC)} m at ${fmt(r.sAtMin, 4)} m above the lower end (still water ${fmt(static0)} m; without shielding ${fmt(rNo.minClearance)} m).`,
        collide ? 'The static shapes cross: contact is predicted in the design current.' : `No static contact; the margin over the required clearance is ${fmt(minC - inp.minClear)} m before wave-frequency motion and VIV amplitudes.`,
        `Upstream line: ${umbUp ? 'umbilical' : 'neighbour'}; largest inclination ${fmt(angMax)} rad.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `Minimum clearance ${fmt(minC)} m against ${inp.minClear} m required; the current consumes ${fmt(static0 - minC)} m and wake shielding accounts for ${fmt(rNo.minClearance - minC)} m of that.`,
          gov.uc > 1 ? (gov.name.startsWith('Small') ? 'The lines incline too much for the string model: run a catenary / rod interference analysis for this pair.' : `Open the spacing by at least ${fmt(gapNeed)} m at the pinch point (hang-off azimuth or separation), equalise the deflections (weight-to-diameter ratio of the umbilical ${fmt(inp.w / (inp.OD / 1e3) / 1e3)} kN/m² against ${fmt(inp.nW / (inp.nD / 1e3) / 1e3)} kN/m² for the neighbour), or fit spacers / a clash-guard at ${fmt(r.sAtMin, 4)} m.`)
            : 'Add the wave-frequency relative motion and the VIV amplitudes of both lines to the static clearance, and repeat for the reversed current direction and for near / far offsets before accepting the layout.',
          'DNV-RP-F203 (2009) §4–5'),
        adv('economic', 'info', collide || gov.uc > 1 ? 'The layout does not keep the required clearance, so either the hang-off arrangement or added hardware has to change.' : 'The layout keeps the required clearance without spacers or clash protection.',
          collide || gov.uc > 1 ? 'Changing hang-off slots at design stage is essentially free; retro-fitted spacers or clash-guards need diver or ROV campaigns and periodic inspection for the rest of the field life — settle it in layout.' : 'Record the clearance margin as a constraint for later tie-ins: adding a line in the gap later will need this analysis repeated, which is far cheaper than discovering sheath wear at inspection.',
          'Risk-based inspection planning DNV-RP-G101; whole-life costing ISO 15686-5',
          `A sheath worn through by line contact floods the armour and typically ends in early replacement of the dynamic section; avoiding one retrofit campaign saves ${vesselCo2Text(3)} Clash-guard shells are polyurethane (≈ ${CO2_POLYMER} t CO₂e/t screening).`),
      ],
      publish: { 'umbilical.minClearance': minC },
    };
  },
  verify: [
    // Two identical lines in a uniform current with no shielding deflect equally: clearance = spacing − diameter = 6 − 0.145 = 5.855 m.
    // Deflection of a weightless tensioned string under uniform load (parabola): y = f·L²/(8T), f = ½ρ·Cd·D·U² =
    // 0.5 × 1025 × 1.2 × 0.145 × 1² = 89.175 N/m → y = 89.175 × 250²/(8 × 80 000) = 8.709 m.
    { name: 'Identical lines without shielding keep their spacing; string deflection y = f·L²/(8T)',
      source: 'Sparks, Fundamentals of Marine Riser Mechanics (2007) ch. 6 (uniform load on a tensioned string); hand calculation',
      inputs: { wake: 'none', gapTop: 6, gapBot: 6, w: 0, nW: 0, nD: 145, nT: 80, nCd: 1.2, Us: 1, Ul: 1, n: 200, minClear: 0 },
      expect: [{ kpi: 'minClearance', value: 5.855, tol: 0.001 }, { kpi: 'deflUp', value: 8.709, tol: 0.003 }] },
  ],
};
