// Line–seabed contact and friction of a single catenary line: grounded length, touchdown, frictional resistance,
// anchor tension and uplift, minimum length for zero uplift, and the embedded-chain attenuation to a buried padeye.
import { solveCatenary, catenaryProfile, minLengthNoUplift } from '../../engine/lines/catenary.js';
import { checkAnchorUplift, checkGroundedLength } from '../../engine/lines/mooringsystem.js';
import { embeddedChainFromMudline } from '../../engine/geotech/anchors.js';
import { DEG, need, run, round, fmt, adv, check, worst, sev, linspace, lineInputs, lineFrom, checkSpan, steelCo2Text, LIMIT_PROPS } from './_common.js';

/** Least-squares mean slope [rad] of an imported seabed profile (polyline along the line, or the central row of a grid). */
function slopeFrom(model) {
  const sArr = [], zArr = [];
  if (model.kind === 'grid' || (model.nx && model.z)) {
    const j = Math.floor((model.ny || 1) / 2);
    for (let i = 0; i < model.nx; i++) { const z = model.z[j * model.nx + i]; if (Number.isFinite(z)) { sArr.push(i * (model.dx || 1)); zArr.push(z); } }
  } else {
    const p = model.points || model.positions;
    need(p && p.length >= 6, 'The imported seabed profile needs at least two points.');
    let s = 0;
    for (let i = 0; i + 2 < p.length; i += 3) { if (i > 0) s += Math.hypot(p[i] - p[i - 3], p[i + 1] - p[i - 2]); sArr.push(s); zArr.push(p[i + 2]); }
  }
  need(sArr.length >= 2, 'The imported seabed profile has fewer than two valid points.');
  const n = sArr.length, ms = sArr.reduce((a, v) => a + v, 0) / n, mz = zArr.reduce((a, v) => a + v, 0) / n;
  let num = 0, den = 0;
  for (let i = 0; i < n; i++) { num += (sArr[i] - ms) * (zArr[i] - mz); den += (sArr[i] - ms) ** 2; }
  need(den > 0, 'The imported seabed profile has no horizontal extent.');
  return Math.atan(num / den);
}

export default {
  id: 's10-seabed-contact',
  suite: 10,
  title: 'Line–seabed contact and friction',
  summary: 'Works out how much line rests on the seabed, the friction it mobilises, the tension and uplift angle reaching the anchor, and the length needed to keep the anchor load horizontal.',
  covers: [8],
  fidelity: 'analytical',
  standards: ['API RP 2SK (3rd ed., 2005) §5.3, §6.5', 'DNV-OS-E301 (2018) Ch.2 Sec.2 [2.5]', 'DNV-RP-E301 (2017) App. A (embedded line)', 'ISO 19901-7 (2013) §8.4'],
  models: ['Elastic catenary with frictional seabed (Jonkman 2007)', 'Coulomb friction on the grounded length', 'Sloping seabed', 'Minimum length for zero uplift (Faltinsen 1990 Eq. 8.20)', 'Embedded chain inverse catenary (Neubecker & Randolph 1995)'],
  limits: `Static, two-dimensional, uniform line on a plane seabed of constant slope; friction is fully mobilised towards the anchor (the upper-bound relief of anchor tension) and lateral friction, trenching and local bathymetry are not modelled. The embedded-chain solution is the small-angle closed form for clay. ${LIMIT_PROPS}`,
  searchTerms: 'mooring chain seabed friction grounded length touchdown anchor uplift',
  inputs: [
    ...lineInputs({ L: 1150, X: 1055 }),
    { key: 'CB', label: 'Seabed friction coefficient', unit: '-', type: 'number', default: 0.6, min: 0, max: 2, group: 'Seabed', role: 'bc', uncertain: { dist: 'uniform', min: 0.4, max: 1.0 },
      help: 'Axial friction between line and seabed: chain on clay or sand 0.6–1.0 (starting), 0.4–0.7 (sliding); wire rope 0.3–0.6. Calibrate against site tests where available.' },
    { key: 'slopeDeg', label: 'Seabed slope (rising towards the fairlead)', unit: 'deg', type: 'number', default: 0, min: -10, max: 10, group: 'Seabed', role: 'bc',
      help: 'Positive when the seabed rises from the anchor towards the vessel. From bathymetry along the line corridor.' },
    { key: 'bathy', label: 'Seabed profile along the line (imported)', type: 'geometry', kinds: ['polyline', 'grid'], default: null, group: 'Seabed', role: 'bc',
      help: 'Optional surveyed seabed: a polyline of points from the anchor towards the fairlead (x, y, z with z up) or a bathymetry grid, whose central row along x is used. The best-fit plane slope replaces the slope entered above; local undulations are not modelled.' },
    { key: 'Tdesign', label: 'Design fairlead tension', unit: 'kN', type: 'number', default: 3100, min: 10, max: 60000, group: 'Design condition', role: 'load', bind: 'mooring.Tmax',
      help: 'Maximum tension from the tension-and-offset module; used for the minimum line length that keeps the anchor free of uplift.' },
    { key: 'su0', label: 'Undrained shear strength at mudline', unit: 'kPa', type: 'number', default: 5, min: 0, max: 500, group: 'Embedded line to padeye', role: 'soil', bind: 'site.su0' },
    { key: 'suGrad', label: 'Strength gradient', unit: 'kPa/m', type: 'number', default: 1.5, min: 0, max: 20, group: 'Embedded line to padeye', role: 'soil', bind: 'site.suGrad' },
    { key: 'zPadeye', label: 'Padeye depth below mudline', unit: 'm', type: 'number', default: 10, min: 0, max: 80, group: 'Embedded line to padeye', role: 'ic',
      help: 'Installed embedment of the anchor attachment point (0 for an anchor on the surface). Initial condition of the anchor: the line cuts through the soil in an inverse catenary down to this depth.' },
    { key: 'muSoil', label: 'Embedded-line friction ratio', unit: '-', type: 'number', default: 0.4, min: 0, max: 1, group: 'Embedded line to padeye', role: 'soil', help: 'Tangential over normal soil resistance: chain 0.4–0.6, wire 0.1–0.3.' },
    { key: 'minGrounded', label: 'Required grounded length', unit: 'm', type: 'number', default: 50, min: 0, max: 2000, group: 'Acceptance criteria', role: 'criteria',
      help: 'Reserve of line on the seabed required in the governing condition (project value; zero if the anchor accepts uplift).' },
    { key: 'upliftAllowDeg', label: 'Allowable uplift angle at mudline', unit: 'deg', type: 'number', default: 0, min: 0, max: 60, group: 'Acceptance criteria', role: 'criteria' },
  ],
  compute(inp, ctx) {
    const p = lineFrom(inp, ctx);
    const { L, X, h } = inp;
    checkSpan(L, X, h, p.EA, p.MBL);
    const slope = inp.bathy ? slopeFrom(inp.bathy) : inp.slopeDeg * DEG;
    need(Math.abs(slope) <= 15 * DEG, `The imported seabed has a mean slope of ${fmt(slope / DEG)}°: outside the ±15° range of the plane-seabed catenary.`);
    const sol = run(() => solveCatenary({ X, Z: h, L, w: p.w, EA: p.EA, CB: inp.CB, slope, nProfile: 0 }), 'Catenary');
    need(!sol.slack, 'The line is fully slack (fairlead too close to the anchor): there is no tension to analyse. Increase the horizontal distance.');
    const pr = catenaryProfile(sol, 120);
    const Lg = sol.groundedLength;
    const wn = p.w * Math.cos(slope);                     // normal contact force per unit length [N/m]
    const friction = Math.max(0, sol.H / Math.cos(slope) - sol.anchorT > 0 ? Math.min(inp.CB * wn * Lg, sol.H / Math.cos(slope)) : 0);
    const frictionCap = inp.CB * wn * Lg;
    const Td = inp.Tdesign * 1e3;
    need(Td > p.w * h, `The design tension ${fmt(Td / 1e3)} kN is below the weight of line hanging straight down (${fmt(p.w * h / 1e3)} kN): it cannot be a fairlead tension for this depth.`);
    const noUp = run(() => minLengthNoUplift({ Z: h, w: p.w, T: Td, EA: p.EA }), 'Minimum length');
    const up = checkAnchorUplift({ upliftAngle: sol.upliftAngle, allowable: inp.upliftAllowDeg * DEG, anchorV: sol.anchorV });
    const gl = checkGroundedLength({ groundedLength: Lg, required: inp.minGrounded });
    // embedded line from mudline to padeye (clay), evaluated in the design condition: mudline tension = horizontal design
    // tension less the friction on whatever length is still grounded at that tension
    const TmDesign = Math.max(noUp.H - inp.CB * wn * Math.max(L - noUp.L, 0), 0);
    let emb = null;
    if (inp.zPadeye > 0 && TmDesign > 0 && (inp.su0 > 0 || inp.suGrad > 0)) {
      try {
        emb = embeddedChainFromMudline({ Tm: TmDesign, za: inp.zPadeye, d: p.d, En: p.family === 'chain' ? 2.5 : 1, Nc: 7.6, su0: inp.su0 * 1e3, k: inp.suGrad * 1e3, mu: inp.muSoil, thetaM: 0 });
        if (!(emb.thetaA < 60 * DEG)) emb = null;       // outside the small-angle range of the closed form
      } catch (e) { emb = null; }                        // tension too small to reach the padeye: held by the soil
    }
    // sweep of fairlead position: grounded length and anchor tension against fairlead tension
    const sw = { T: [], gr: [], Ta: [], x: [] };
    for (const x of linspace(Math.max(1, X - 0.06 * L), X + 0.06 * L, 41)) {
      let s; try { s = solveCatenary({ X: x, Z: h, L, w: p.w, EA: p.EA, CB: inp.CB, slope }); } catch (e) { break; }
      if (s.slack || s.T > p.MBL) continue;
      sw.x.push(round(x, 6)); sw.T.push(round(s.T / 1e3)); sw.gr.push(round(s.groundedLength)); sw.Ta.push(round(s.anchorT / 1e3));
    }
    const checks = [
      { name: 'Uplift angle at the mudline', demand: round(sol.upliftAngle / DEG, 4), capacity: inp.upliftAllowDeg, unit: 'deg', uc: inp.upliftAllowDeg > 0 ? sol.upliftAngle / DEG / inp.upliftAllowDeg : (up.ok ? 0 : 99),
        ref: 'API RP 2SK (3rd ed.) §6.5', explain: 'A drag anchor needs a horizontal pull. Any uplift angle with zero allowable means the whole line is lifted and the anchor can be pulled out of the seabed.' },
      check('Grounded length reserve', inp.minGrounded, Math.max(Lg, 1e-9), 'm', 'Design practice; API RP 2SK (3rd ed.) §6.5',
        'Required over available grounded length. The grounded line is the reserve against uplift and adds frictional holding; too little means small extra offsets start lifting the anchor.'),
      check('Line length for zero uplift at design tension', noUp.L, L, 'm', 'Faltinsen (1990) Eq. 8.20; API RP 2SK (3rd ed.) §6.5',
        'Minimum length for the line to remain tangent to the seabed at the anchor when the design tension is applied; a shorter line lifts the anchor at that tension.'),
    ];
    const g = worst(checks);
    const mass = p.mass * L / 1e3;
    const extra = Math.max(0, noUp.L + inp.minGrounded - L);
    const spare = Math.max(0, L - noUp.L - inp.minGrounded);
    return {
      kpis: [
        { key: 'T', label: 'Fairlead tension', value: round(sol.T / 1e3, 5), unit: 'kN', explain: 'Tension at the vessel end for the entered geometry; friction on the seabed does not change it, only the anchor tension.' },
        { key: 'H', label: 'Horizontal tension', value: round(sol.H / 1e3, 5), unit: 'kN', explain: 'Horizontal component at the fairlead, constant along the suspended part; the force the grounded line and anchor together must resist.' },
        { key: 'grounded', label: 'Grounded (contact) length', value: round(Lg, 5), unit: 'm', limit: inp.minGrounded, sense: 'min', explain: 'Unstretched line length resting on the seabed between anchor and touchdown.' },
        { key: 'touchdownX', label: 'Touchdown distance from anchor', value: round(sol.touchdownX, 5), unit: 'm', explain: 'Horizontal position where the line lifts off the seabed; the zone of repeated lifting, wear and trenching.' },
        { key: 'contactForce', label: 'Seabed contact force', value: round(wn * Lg / 1e3, 5), unit: 'kN', explain: 'Total normal force between grounded line and seabed (submerged weight of the grounded length); the distributed contact force is the submerged weight per metre.' },
        { key: 'contactPerLength', label: 'Contact force per unit length', value: round(wn, 5), unit: 'N/m', explain: 'Normal seabed reaction per metre of grounded line; with the line diameter it gives the bearing pressure that drives self-burial.' },
        { key: 'frictionCap', label: 'Available friction on grounded length', value: round(frictionCap / 1e3, 5), unit: 'kN', explain: 'Friction coefficient × contact force: the most tension the seabed can take out of the line before it reaches the anchor.' },
        { key: 'anchorT', label: 'Tension at the anchor (mudline)', value: round(sol.anchorT / 1e3, 5), unit: 'kN', explain: 'Line tension reaching the anchor shackle after friction; zero means the grounded line alone holds the vessel and the anchor is unloaded.' },
        { key: 'frictionUsed', label: 'Friction mobilised', value: round(friction / 1e3, 5), unit: 'kN', explain: 'Tension actually removed by seabed friction between touchdown and anchor (cannot exceed the horizontal tension).' },
        { key: 'uplift', label: 'Uplift angle at the anchor', value: round(sol.upliftAngle / DEG, 4), unit: 'deg', explain: 'Angle of the line above the seabed at the anchor; zero while any line remains grounded.' },
        { key: 'Lmin', label: 'Minimum length for zero uplift at design tension', value: round(noUp.L, 5), unit: 'm', limit: round(L, 5), sense: 'max', explain: 'Shortest line that stays tangent to the seabed at the anchor under the design tension; the actual length should exceed it by the required grounded reserve.' },
        { key: 'mudlineTd', label: 'Mudline tension at design tension', value: round(TmDesign / 1e3, 5), unit: 'kN', explain: 'Tension arriving at the anchor point on the seabed when the design fairlead tension is applied, after friction on the length still grounded.' },
        { key: 'padeyeT', label: 'Padeye tension at design tension', value: round((emb ? emb.Ta : TmDesign) / 1e3, 5), unit: 'kN', explain: 'Tension left at the buried padeye after the embedded line has shed load into the soil; the load the anchor body itself must hold.' },
        { key: 'padeyeAngle', label: 'Line angle at the padeye', value: round((emb ? emb.thetaA : 0) / DEG, 4), unit: 'deg', explain: 'Inclination of the embedded line at the padeye (inverse catenary). It sets the split into horizontal and vertical anchor load even when the mudline angle is zero.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Line profile and seabed', x: { label: 'Horizontal distance from anchor', unit: 'm' }, y: { label: 'Elevation above anchor', unit: 'm' },
          series: [{ name: 'Line', x: Array.from(pr.x, (v) => round(v)), y: Array.from(pr.z, (v) => round(v)) }, { name: 'Seabed', x: [0, X], y: [0, round(X * Math.tan(slope))], dash: true }],
          explain: 'Catenary shape over the seabed. The part lying on the dashed seabed line is the grounded length; where the curve leaves it is the touchdown point.' },
        { type: 'line', title: 'Tension along the line', x: { label: 'Arc length from anchor', unit: 'm' }, y: { label: 'Tension', unit: 'kN' },
          series: [{ name: 'Tension', x: Array.from(pr.s, (v) => round(v)), y: Array.from(pr.T, (v) => round(v / 1e3)) }],
          explain: 'Tension from anchor to fairlead. The rising ramp on the grounded part is seabed friction (slope = friction × weight per metre); a flat zero stretch means the anchor is not loaded at all.' },
        { type: 'line', title: 'Grounded length and anchor tension vs fairlead tension', x: { label: 'Fairlead tension', unit: 'kN' }, y: { label: 'Grounded length (m) / anchor tension (kN)', unit: '' },
          series: [{ name: 'Grounded length (m)', x: sw.T, y: sw.gr }, { name: 'Anchor tension (kN)', x: sw.T, y: sw.Ta }],
          explain: 'How the contact length shrinks and the anchor tension grows as the vessel moves away. The tension at which the grounded length reaches zero is the lift-off limit of this line.' },
      ],
      tables: [
        { title: 'Contact and friction summary', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['Suspended length', sol.suspendedLength, 'm'], ['Grounded length', Lg, 'm'], ['Stretched length', sol.stretchedLength, 'm'], ['Contact force per metre', wn, 'N/m'], ['Friction capacity', frictionCap / 1e3, 'kN'],
            ['Tension at touchdown', sol.H / Math.cos(slope) / 1e3, 'kN'], ['Tension at anchor (mudline)', sol.anchorT / 1e3, 'kN'], ['Anchor horizontal load', sol.anchorH / 1e3, 'kN'], ['Anchor vertical load', sol.anchorV / 1e3, 'kN'],
            ['Mudline tension at design tension', TmDesign / 1e3, 'kN'], ['Padeye tension at design tension', (emb ? emb.Ta : TmDesign) / 1e3, 'kN'], ['Horizontal stiffness dH/dX', sol.stiffness[0][0] / 1e3, 'kN/m']].map(([q, v, u]) => ({ q, v: round(v, 5), u })),
          explain: 'Geometry and force balance of the line on the seabed. Friction capacity minus tension at touchdown tells whether the anchor is loaded at all.' },
      ],
      notes: [`${fmt(Lg)} m of ${fmt(L)} m rests on the seabed; friction can remove up to ${fmt(frictionCap / 1e3)} kN of the ${fmt(sol.H / 1e3)} kN horizontal tension.`,
        emb ? `Embedded line to a padeye at ${inp.zPadeye} m in the design condition: tension falls from ${fmt(TmDesign / 1e3)} to ${fmt(emb.Ta / 1e3)} kN and the line arrives at ${fmt(emb.thetaA / DEG)}° to the horizontal.` : 'No embedded-line result: padeye at mudline, no soil strength, or the mudline tension is too small for the small-angle inverse-catenary solution (the soil holds the line).',
        'Friction is taken as fully mobilised towards the anchor; for anchor design also check the zero-friction case.'],
      advice: [
        adv('engineering', sev(g.uc), `Grounded length ${fmt(Lg)} m, uplift angle ${fmt(sol.upliftAngle / DEG)}°, anchor tension ${fmt(sol.anchorT / 1e3)} kN; at the design tension of ${inp.Tdesign} kN the line needs ${fmt(noUp.L)} m to avoid uplift (governing: ${g.name}, ${fmt(g.uc)}).`,
          g.uc > 1 ? `Add about ${fmt(Math.max(extra, noUp.L - L))} m of ground chain (or move the anchor out by a similar distance with the same pretension), or select an anchor rated for a ${fmt(Math.max(sol.upliftAngle / DEG, 5))}° uplift angle.`
            : `The line has ${fmt(spare)} m more than needed for zero uplift plus the ${inp.minGrounded} m reserve. Design the anchor for the no-friction tension ${fmt(sol.H / 1e3)} kN as an upper bound and ${fmt(sol.anchorT / 1e3)} kN as best estimate.`,
          'API RP 2SK (3rd ed., 2005) §6.5; DNV-OS-E301 (2018) Ch.2 Sec.2 [2.5] (friction to be bounded: high for line tension, low for anchor load).'),
        adv('economic', 'info', `Ground chain costs mass: each 100 m of this line is ${fmt(p.mass * 100 / 1e3)} t.`,
          spare > 20 ? `About ${fmt(spare)} m (${fmt(p.mass * spare / 1e3)} t) could be removed while keeping zero uplift and the required reserve at the design tension — verify against the damaged-condition offset first.` : extra > 0 ? `The shortfall of ${fmt(extra)} m is ${fmt(p.mass * extra / 1e3)} t of extra chain; an uplift-rated anchor (suction, pile, VLA) may be the cheaper answer in deep water.` : 'Length is matched to the design tension; no saving available here.',
          'Whole-life costing per ISO 15686-5; quantity-based Class 4 estimate per AACE RP 18R-97.',
          `${steelCo2Text(spare > 20 ? p.mass * spare / 1e3 : extra > 0 ? p.mass * extra / 1e3 : mass)} Grounded chain also sweeps and trenches the seabed over the touchdown zone; shorter ground chain reduces that footprint.`),
      ],
      quantities: p.family === 'fibre' ? { polymer_t: round(mass) } : { steel_t: round(mass) },
    };
  },
  verify: [
    { // Inextensible chain, w = 1500 N/m, h = 280 m, H = 800 kN: a = 533.33 m, suspended s = sqrt(h² + 2ha) = 614.06 m,
      // x = a·acosh(1 + h/a) = 525.04 m. With L = 1000 m: grounded 385.94 m, X = 910.99 m, T = H + w·h = 1220 kN.
      // Friction CB = 0.5: anchor tension = H − CB·w·Lg = 800 − 0.5·1.5·385.94 = 510.5 kN. Minimum length for no uplift at
      // T = 1220 kN equals the suspended length 614.06 m.
      name: 'Inextensible catenary with seabed friction',
      source: 'Faltinsen (1990) Sea Loads, Eqs. 8.14–8.20; Coulomb friction on the grounded length',
      inputs: { lineType: 'custom', wCustom: 1500, EAcustom: 1e7, MBLcustom: 9000, massCustom: 176, L: 1000, X: 910.986, h: 280, CB: 0.5, slopeDeg: 0, Tdesign: 1220, zPadeye: 0 },
      expect: [{ kpi: 'T', value: 1220, tol: 0.002 }, { kpi: 'grounded', value: 385.94, tol: 0.003 }, { kpi: 'anchorT', value: 510.5, tol: 0.005 }, { kpi: 'Lmin', value: 614.06, tol: 0.002 }, { kpi: 'contactForce', value: 578.9, tol: 0.003 }] },
  ],
};
