// Accidental impact and dropped objects: impact energy at the seabed, local denting, plastic deformation and residual capacity.
import { tubeSection } from '../../engine/structural/sections.js';
import { G, DEG, need, fmt, round, adv, check, steelCo2Text, vesselCo2Text, CONFIRM } from './_common.js';

/** Abramowitz & Stegun 7.1.26 error function (|error| < 1.5e-7). */
function erf(x) {
  const t = 1 / (1 + 0.3275911 * Math.abs(x));
  const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
  return x >= 0 ? y : -y;
}

/** Plastic section modulus of a thin ring with the dented arc (half-angle alpha at the crown) ineffective, relative to the intact ring. */
function residualPlasticRatio(alpha) {
  if (alpha <= 0) return 1;
  const n = 1440, ys = [];
  for (let i = 0; i < n; i++) { const th = (i + 0.5) * 2 * Math.PI / n; const a = Math.min(Math.abs(th), 2 * Math.PI - Math.abs(th)); if (a > alpha) ys.push(Math.cos(th)); } // y = cos θ, crown at θ = 0
  ys.sort((p, q) => p - q);
  const yPna = ys[Math.floor(ys.length / 2)];
  const zp = ys.reduce((s, y) => s + Math.abs(y - yPna), 0) * (2 * Math.PI / n);
  return zp / 4; // intact thin ring: ∮|cos θ| dθ = 4
}

export default {
  id: 's04-dropped-object', suite: 4, title: 'Accidental impact and dropped objects',
  summary: 'Velocity and kinetic energy of an object dropped from the vessel when it reaches the structure, the dent and plastic deformation it causes in a protection member, roof plate or pipe, and the residual capacity.',
  covers: [15], fidelity: 'screening',
  standards: ['DNV-RP-F107 (2019) §4–5 (dropped objects, impact capacity)', 'DNV-RP-C204 (2019) §3–4 (accidental loads, energy absorption)', 'NORSOK N-004 (2013) Annex A.3–A.4', 'NORSOK U-001 (2015) §5.3 (dropped-object protection)'],
  models: ['Energy-based dropped-object method', 'Terminal velocity in water with added mass', 'Local denting of a tube (Wierzbicki & Suh, as in DNV-RP-F107)', 'Rigid-plastic beam mechanism', 'Rigid-plastic plate collapse under a patch load', 'Dented-section residual capacity (ineffective dented arc)', 'Normal distribution of lateral excursion'],
  limits: 'Energy balance with rigid-plastic resistance: the denting expression, drag and added-mass coefficients, excursion angles and the dent damage bands of DNV-RP-F107 are recalled and exposed as inputs to be confirmed against the governing edition. Membrane action in plates and beams is neglected (conservative), the object is rigid, and strain-rate effects, fracture and contact stresses are not modelled; explicit nonlinear finite-element impact analysis is outside the module.',
  publishes: ['subsea.impactEnergy', 'subsea.dentRatio'],
  searchTerms: 'dropped object subsea protection structure impact energy DNV-RP-F107 dent',
  inputs: [
    { key: 'mObj', label: 'Object mass', unit: 'kg', type: 'number', default: 1000, min: 1, max: 5e5, group: 'Dropped object', role: 'load', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Typical crane-handled items: 1–8 t (containers, baskets, casing joints, BOP parts).' },
    { key: 'volObj', label: 'Object volume', unit: 'm³', type: 'number', default: 0.13, min: 0, max: 500, group: 'Dropped object', role: 'geometry' },
    { key: 'Aobj', label: 'Projected area in the fall direction', unit: 'm²', type: 'number', default: 1.0, min: 0.001, max: 200, group: 'Dropped object', role: 'geometry', help: 'Flat or broadside for boxes and tumbling tubulars; end-on area gives the highest velocity.' },
    { key: 'CdObj', label: 'Drag coefficient', unit: '-', type: 'number', default: 1.0, min: 0.2, max: 3, group: 'Dropped object', role: 'load', uncertain: { dist: 'normal', cov: 0.2 }, help: `DNV-RP-F107 Table: slender 0.7–1.5, box 1.2–1.3, complex 0.6–2.0. ${CONFIRM}` },
    { key: 'CaObj', label: 'Added-mass coefficient', unit: '-', type: 'number', default: 1.0, min: 0, max: 2.5, group: 'Dropped object', role: 'load', help: `DNV-RP-F107 Table: slender 0.1–1.0, box 0.6–1.5. ${CONFIRM}` },
    { key: 'hAir', label: 'Drop height in air', unit: 'm', type: 'number', default: 20, min: 0, max: 150, group: 'Initial conditions', role: 'ic', help: 'Height of the load above the sea surface when released; sets the water-entry velocity.' },
    { key: 'depth', label: 'Water depth above the target', unit: 'm', type: 'number', default: 150, min: 1, max: 3500, group: 'Initial conditions', role: 'ic', bind: 'site.depth' },
    { key: 'devDeg', label: 'Excursion angle (one standard deviation)', unit: 'deg', type: 'number', default: 5, min: 0.5, max: 30, group: 'Initial conditions', role: 'criteria', help: `DNV-RP-F107 Table: flat / long objects up to 15°, box-shaped about 3–5°, heavy compact about 2°. ${CONFIRM}` },
    { key: 'rTarget', label: 'Equivalent radius of the structure in plan', unit: 'm', type: 'number', default: 6.7, min: 0.2, max: 100, group: 'Initial conditions', role: 'geometry', help: '√(B·L/π) for a rectangular template.' },
    { key: 'target', label: 'Impacted component', type: 'select', default: 'tube', group: 'Target', role: 'bc',
      options: [{ value: 'tube', label: 'Tubular protection member (dent + beam mechanism)' }, { value: 'plate', label: 'Roof plate panel' }, { value: 'pipe', label: 'Pressure-containing pipe or spool (dent only)' }], help: 'Sets the contact location and the mechanism by which the impact energy is absorbed.' },
    { key: 'mat', label: 'Steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Target', role: 'material' },
    { key: 'Dt', label: 'Tube / pipe outer diameter', unit: 'mm', type: 'number', default: 406.4, min: 50, max: 3000, group: 'Target', role: 'geometry', design: true, showIf: { key: 'target', in: ['tube', 'pipe'] } },
    { key: 'tt', label: 'Tube / pipe wall thickness', unit: 'mm', type: 'number', default: 20, min: 3, max: 100, group: 'Target', role: 'geometry', design: true, showIf: { key: 'target', in: ['tube', 'pipe'] } },
    { key: 'Lt', label: 'Member span', unit: 'm', type: 'number', default: 5, min: 0.5, max: 60, group: 'Target', role: 'geometry', showIf: { key: 'target', in: ['tube'] } },
    { key: 'endFix', label: 'Member end restraint', type: 'select', default: 'fixed', group: 'Target', role: 'bc', showIf: { key: 'target', in: ['tube'] }, options: [{ value: 'fixed', label: 'Fixed ends (collapse load 8·Mp/L)' }, { value: 'pinned', label: 'Pinned ends (collapse load 4·Mp/L)' }] },
    { key: 'tp', label: 'Roof plate thickness', unit: 'mm', type: 'number', default: 12, min: 3, max: 80, group: 'Target', role: 'geometry', design: true, showIf: { key: 'target', in: ['plate'] } },
    { key: 'Rp', label: 'Panel radius (half the shorter span)', unit: 'm', type: 'number', default: 0.6, min: 0.1, max: 10, group: 'Target', role: 'geometry', showIf: { key: 'target', in: ['plate'] } },
    { key: 'aPatch', label: 'Contact patch radius', unit: 'm', type: 'number', default: 0.1, min: 0.005, max: 5, group: 'Target', role: 'bc', showIf: { key: 'target', in: ['plate'] }, help: 'Radius of the contact area of the object on the plate.' },
    { key: 'eta', label: 'Share of energy absorbed by the target', unit: '-', type: 'number', default: 1.0, min: 0.1, max: 1, group: 'Acceptance criteria', role: 'criteria', help: 'Less than 1 when the object itself deforms or coating / soil absorbs part; 1.0 is conservative.' },
    { key: 'dentAllow', label: 'Allowable dent depth / diameter', unit: '%', type: 'number', default: 5, min: 1, max: 50, group: 'Acceptance criteria', role: 'criteria', help: `DNV-RP-F107 (2017) Table 4-1 for steel pipes: below 5 % is minor damage; 5–10 % is already major damage with leakage anticipated; above 10 % leakage and rupture are anticipated. The default is the minor-damage limit, which is the one to use for pressure-containing pipes. For protection structures and sacrificial frame members a larger dent may be accepted — raise it deliberately. ${CONFIRM}` },
    { key: 'dMinor', label: 'Dent limit of the minor-damage class', unit: '%', type: 'number', default: 5, min: 1, max: 20, group: 'Acceptance criteria', role: 'criteria', help: `DNV-RP-F107 damage classification: dents up to this depth need no immediate repair. ${CONFIRM}` },
    { key: 'dMajor', label: 'Dent depth above which rupture is anticipated', unit: '%', type: 'number', default: 10, min: 2, max: 40, group: 'Acceptance criteria', role: 'criteria', help: `DNV-RP-F107 (2017) Table 4-1: between the minor limit and this depth the damage is major with leakage anticipated; above it both leakage and rupture are anticipated. ${CONFIRM}` },
    { key: 'defAllow', label: 'Allowable permanent deflection / span', unit: '%', type: 'number', default: 5, min: 0.5, max: 30, group: 'Acceptance criteria', role: 'criteria', help: 'Limits beam or plate deformation so the protected equipment is not struck and the strain stays ductile.' },
    { key: 'Edesign', label: 'Specified design impact energy', unit: 'kJ', type: 'number', default: 0, min: 0, max: 1e5, group: 'Acceptance criteria', role: 'criteria', help: `If greater than zero it replaces the computed energy (e.g. NORSOK U-001: 50 kJ over a 700 mm diameter for multi-well structures, 20 kJ / 5 kJ classes). ${CONFIRM}` },
  ],
  compute(inp, ctx) {
    // SI: mm → m, kJ → J, % → fraction.
    const steel = ctx.material(inp.mat);
    need(steel, `Steel '${inp.mat}' is not in the materials library.`);
    const rho = 1025, fy = steel.fy;
    const m = inp.mObj, V = inp.volObj, ma = rho * inp.CaObj * V;
    const Wsub = (m - rho * V) * G;
    need(Wsub > 0, 'The object floats (mass below displaced water): it does not reach the seabed.');
    // Free fall in air, then (m + ma)·v·dv/dz = W′ − ½ρ·Cd·A·v² in water.
    const v0 = Math.sqrt(2 * G * inp.hAir);
    const vT = Math.sqrt(2 * Wsub / (rho * inp.CdObj * inp.Aobj));
    const beta = rho * inp.CdObj * inp.Aobj / (m + ma);
    const vAt = (z) => Math.sqrt(vT * vT + (v0 * v0 - vT * vT) * Math.exp(-beta * z));
    const v = vAt(inp.depth);
    const Ecalc = 0.5 * (m + ma) * v * v;
    const E = (inp.Edesign > 0 ? inp.Edesign * 1e3 : Ecalc) * inp.eta;
    const sig = inp.depth * Math.tan(inp.devDeg * DEG);
    const pHit = erf(inp.rTarget / (sig * Math.SQRT2));
    const kpis = [], checks = [], rows = [];
    let dentRatio = 0, defl = 0, cap = 0, Fpeak = 0, resN = 1, resM = 1, tReq = 0, steelDelta = 0, what = '', sigImp = 0;
    const xs = [], Es = [];
    if (inp.target === 'plate') {
      const t = inp.tp / 1e3, R = inp.Rp, a = Math.min(inp.aPatch, 0.9 * R);
      const mp = fy * t * t / 4;
      // Simply supported circular plate under a central patch load (Hopkins & Prager): Pc = 2π·mp / (1 − 2a/3R).
      const Pc = 2 * Math.PI * mp / (1 - 2 * a / (3 * R));
      const wAllow = inp.defAllow / 100 * 2 * R;
      cap = Pc * wAllow; defl = E / Pc; Fpeak = Pc;
      sigImp = 6 * mp / (t * t); // the collapse mechanism means the full plastic moment is reached: extreme fibre at 1.5·fy nominal elastic
      tReq = t * Math.sqrt(E / cap) * 1e3; steelDelta = Math.max(tReq - inp.tp, 0) * 7.85e-3 * 140;
      what = 'roof plate';
      checks.push(check('Roof plate: energy against plastic capacity', E / 1e3, cap / 1e3, 'kJ', 'DNV-RP-C204 (2019) §3 (energy dissipation by plastic deformation); rigid-plastic plate collapse', 'Guards against the object punching through or deflecting the roof onto the equipment. Membrane action is neglected, so the capacity is a lower bound.'));
      checks.push(check('Roof plate: permanent deflection', defl * 1e3, wAllow * 1e3, 'mm', 'Project clearance and ductility limit on plate deflection', 'Permanent dishing of the panel after impact; it must stay clear of the protected equipment.'));
      rows.push(['Plate plastic moment mp', mp / 1e3, 'kN·m/m'], ['Plate collapse load', Pc / 1e3, 'kN'], ['Allowable deflection', wAllow * 1e3, 'mm']);
      for (let i = 0; i <= 20; i++) { const w = wAllow * 2 * i / 20; xs.push(round(w * 1e3)); Es.push(Pc * w / 1e3); }
      kpis.push({ key: 'defl', label: 'Permanent deflection of the plate', unit: 'mm', value: defl * 1e3, limit: wAllow * 1e3, sense: 'max', explain: 'Rigid-plastic deflection needed to absorb the impact energy at the plate collapse load, without membrane stiffening.' });
    } else {
      const D = inp.Dt / 1e3, t = inp.tt / 1e3;
      need(t < D / 2, 'Wall thickness must be smaller than the radius.');
      const mp = fy * t * t / 4;
      // DNV-RP-F107 denting: E = 16·(2π/9)^½·mp·(D/t)^½·D·(δ/D)^(3/2).
      const cD = 16 * Math.sqrt(2 * Math.PI / 9) * mp * Math.sqrt(D / t) * D;
      const Edent = (r) => cD * r ** 1.5, Fdent = (r) => 1.5 * cD * Math.sqrt(r) / D;
      const sec = tubeSection({ D, t });
      const Mp = fy * sec.Wply;
      let Ebeam = 0, Pc = Infinity, rYield = 1;
      if (inp.target === 'tube') {
        Pc = (inp.endFix === 'fixed' ? 8 : 4) * Mp / inp.Lt;
        rYield = Math.min((Pc * D / (1.5 * cD)) ** 2, 1); // dent ratio at which the denting force reaches the beam collapse load
      }
      const dAllow = inp.dentAllow / 100;
      if (E <= Edent(rYield)) dentRatio = (E / cD) ** (2 / 3);
      else { dentRatio = rYield; Ebeam = E - Edent(rYield); defl = Ebeam / Pc; }
      Fpeak = Fdent(dentRatio);
      // Nominal elastic bending stress of the member under the peak contact force (span moment FL/8 fixed, FL/4 pinned).
      sigImp = inp.target === 'tube' ? Math.min(Fpeak, Pc) * inp.Lt / (inp.endFix === 'fixed' ? 8 : 4) / sec.Wely : Fpeak / (Math.PI * D * t);
      const wAllow = inp.target === 'tube' ? inp.defAllow / 100 * inp.Lt : 0;
      cap = Edent(Math.min(dAllow, rYield)) + (inp.target === 'tube' && dAllow >= rYield ? Pc * wAllow : 0);
      // Dented section: flattened arc treated as ineffective.
      const alpha = Math.acos(Math.max(1 - 2 * Math.min(dentRatio, 1), -1));
      resN = 1 - alpha / Math.PI; resM = residualPlasticRatio(alpha);
      tReq = t * (E / Math.max(Edent(dAllow), 1e-9)) ** 0.4 * 1e3; // dent energy ∝ t^2.5 at a given dent ratio
      steelDelta = inp.target === 'tube' ? Math.max(tReq - inp.tt, 0) / 1e3 * Math.PI * D * inp.Lt * 7.85 : 0;
      what = inp.target === 'tube' ? 'protection member' : 'pipe';
      checks.push(check('Dent depth', dentRatio * 100, inp.dentAllow, '% of diameter', 'DNV-RP-F107 (2019) §4 (dent–energy relation and damage classes)', 'Guards against a dent deep enough to rupture a pipe or cripple a member. For pressure-containing pipes dents above about 5 % need assessment.'));
      if (inp.target === 'tube') checks.push(check('Member permanent deflection', defl * 1e3, wAllow * 1e3, 'mm', 'DNV-RP-C204 (2019) §3; NORSOK N-004 A.3 (rigid-plastic beam mechanism)', 'Once the denting force reaches the collapse load the member forms plastic hinges; the deflection must stay clear of the protected equipment.'));
      checks.push(check('Impact energy against absorption capacity', E / 1e3, cap / 1e3, 'kJ', 'DNV-RP-C204 (2019) §3 energy method; DNV-RP-F107 (2019) §4', 'Energy that can be absorbed within the allowable dent and deflection. Above it the protection does not stop the object.'));
      rows.push(['Wall plastic moment mp', mp / 1e3, 'kN·m/m'], ['Energy to a 5 % dent', Edent(0.05) / 1e3, 'kJ'], ['Energy to the allowable dent', Edent(dAllow) / 1e3, 'kJ'], ['Section plastic moment Mp', Mp / 1e3, 'kN·m']);
      if (inp.target === 'tube') rows.push(['Beam collapse load', Pc / 1e3, 'kN'], ['Dent ratio at onset of the beam mechanism', rYield * 100, '%'], ['Energy absorbed in beam bending', Ebeam / 1e3, 'kJ']);
      for (let i = 0; i <= 30; i++) { const r = 0.3 * i / 30; xs.push(round(r * 100)); Es.push(Edent(r) / 1e3); }
      kpis.push(
        { key: 'dent', label: 'Dent depth / diameter', unit: '%', value: dentRatio * 100, limit: inp.dentAllow, sense: 'max', explain: 'Permanent local dent from the absorbed energy by the DNV-RP-F107 relation. Compare with the damage classes for pipes: below 5 % minor; 5–10 % major, leakage anticipated; above 10 % rupture anticipated.' },
        { key: 'Edent5', label: 'Energy to produce a 5 % dent', unit: 'kJ', value: Edent(0.05) / 1e3, explain: 'Impact energy that produces a dent of 5 % of the diameter: the usual threshold between minor and moderate damage for pipes.' },
        { key: 'defl', label: 'Permanent deflection of the member', unit: 'mm', value: defl * 1e3, explain: 'Rigid-plastic deflection of the member after the denting force has reached its collapse load. Zero when the energy is absorbed by denting alone.' },
        { key: 'resN', label: 'Residual axial capacity', unit: '%', value: resN * 100, explain: 'Squash-load fraction remaining when the dented arc is treated as ineffective. Use it to re-check the damaged member for the post-accident condition.' },
        { key: 'resM', label: 'Residual bending capacity', unit: '%', value: resM * 100, explain: 'Plastic moment fraction remaining with the dented arc ineffective, bending with the dent in compression.' });
    }
    const ok = checks.every((c) => c.uc <= 1);
    const w = checks.reduce((a, c) => (c.uc > a.uc ? c : a));
    const zs = [], vs = [], ek = [];
    for (let i = 0; i <= 40; i++) { const z = inp.depth * i / 40; zs.push(round(-z)); vs.push(vAt(z)); ek.push(0.5 * (m + ma) * vAt(z) ** 2 / 1e3); }
    const advice = [
      ok ? adv('engineering', 'ok', `Impact energy ${fmt(E / 1e3)} kJ (velocity ${fmt(v)} m/s, terminal ${fmt(vT)} m/s) is absorbed by the ${what} within the limits: capacity ${fmt(cap / 1e3)} kJ, governing check "${w.name}" at ${fmt(w.uc)}.`,
        `Re-check the damaged member with ${fmt(resN * 100)} % axial and ${fmt(resM * 100)} % bending capacity in the post-accident condition, and confirm the object list, drop frequencies and the F107 coefficients against the governing edition.`,
        'DNV-RP-F107 (2019) §4–5; DNV-RP-C204 (2019) §3; NORSOK U-001 (2015) §5.3',
        'The protection is adequate without extra steel; over-sizing every roof member by one wall-thickness step would add several tonnes for no reduction in risk.')
        : adv('engineering', 'fail', `"${w.name}" is exceeded (${fmt(w.uc)}): the ${what} cannot absorb ${fmt(E / 1e3)} kJ within the allowable deformation (capacity ${fmt(cap / 1e3)} kJ).`,
          inp.target === 'plate' ? `Increase the plate to about ${fmt(Math.ceil(tReq))} mm (collapse load ∝ t²), reduce the panel span with an extra stiffener, or restrict lifting over the structure.` : `Increase the wall to about ${fmt(Math.ceil(tReq))} mm (dent energy ∝ t^2.5), add a sacrificial cover or GRP protection, or restrict the lifting zone so the object cannot land here.`,
          'DNV-RP-F107 (2019) §4–5; DNV-RP-C204 (2019) §3; NORSOK U-001 (2015) §5.3',
          steelDelta > 0 ? `Strengthening adds about ${steelCo2Text(steelDelta)} Operational restriction of lifts over the structure achieves the risk reduction with no material.` : 'A lifting-zone restriction or a lighter composite cover reduces risk with far less material than thicker steel.'),
      adv('economic', ok ? 'info' : 'warn', `${fmt(pHit * 100)} % of drops from directly above land within the structure radius (one-sigma excursion ${fmt(sig)} m at this depth).`,
        'Dropped-object protection is a risk trade: compare the steel and installation weight of protection against the frequency × consequence of a hit (repair spread of several vessel days plus deferred production). Moving the lifting zone or the vessel heading often lowers the hit probability at no capital cost.',
        'Risk assessment per DNV-RP-F107 (2019) §2 and ISO 17776; whole-life costing per ISO 15686-5',
        `A repair campaign after a hit would cost at least ${vesselCo2Text(5)}`),
    ];
    return {
      kpis: [
        { key: 'Eimp', label: 'Impact energy at the target', unit: 'kJ', value: Ecalc / 1e3, explain: 'Kinetic energy of the object including its hydrodynamic added mass when it reaches the structure.' },
        { key: 'Eabs', label: 'Energy to be absorbed', unit: 'kJ', value: E / 1e3, limit: cap / 1e3, sense: 'max', explain: 'Design energy (computed, or the specified value if given) times the share taken by the target, against the capacity within allowable deformation.' },
        { key: 'vImp', label: 'Impact velocity', unit: 'm/s', value: v, explain: 'Velocity at the depth of the structure; it approaches the terminal velocity within a few tens of metres for most objects.' },
        { key: 'vT', label: 'Terminal velocity in water', unit: 'm/s', value: vT, explain: 'Velocity at which submerged weight balances drag. It depends strongly on the orientation (projected area) of the falling object.' },
        { key: 'vEntry', label: 'Water-entry velocity', unit: 'm/s', value: v0, explain: 'Velocity after the free fall in air, before any loss at water entry (conservative).' },
        { key: 'Fpeak', label: 'Peak contact force', unit: 'kN', value: Fpeak / 1e3, explain: 'Largest contact force during the impact from the force–deformation relation of the target.' },
        { key: 'sigImp', label: 'Nominal impact stress', unit: 'MPa', value: sigImp / 1e6, explain: 'Elastic stress under the peak contact force: member bending stress (protection member), plate extreme-fibre stress at collapse (roof plate) or contact force over the pipe wall section (pipe). A value above yield confirms that the impact is absorbed plastically.' },
        { key: 'capacity', label: 'Energy absorption capacity', unit: 'kJ', value: cap / 1e3, explain: 'Energy the component can absorb within the allowable dent and deflection.' },
        { key: 'sigExc', label: 'Lateral excursion (one sigma)', unit: 'm', value: sig, explain: 'Standard deviation of the landing point from the vertical below the drop point: depth × tan(excursion angle).' },
        { key: 'pHit', label: 'Share of drops landing within the structure', unit: '%', value: pHit * 100, explain: 'Probability that a drop from directly above lands within the equivalent radius of the structure, for a normally distributed excursion.' },
        ...kpis,
      ],
      checks,
      plots: [
        { type: 'line', title: 'Velocity and kinetic energy against depth', x: { label: 'Velocity (m/s) · energy (kJ)', unit: '' }, y: { label: 'Depth below surface', unit: 'm' },
          series: [{ name: 'Velocity (m/s)', x: vs, y: zs }, { name: 'Kinetic energy incl. added mass (kJ)', x: ek, y: zs }],
          explain: 'How the object slows from its entry velocity (or accelerates) to terminal velocity. Below the depth where the curves turn vertical, depth no longer matters.' },
        { type: 'line', title: inp.target === 'plate' ? 'Absorbed energy against plate deflection' : 'Dent energy against dent depth', x: { label: inp.target === 'plate' ? 'Deflection' : 'Dent depth / diameter', unit: inp.target === 'plate' ? 'mm' : '%' }, y: { label: 'Absorbed energy', unit: 'kJ' },
          series: [{ name: 'Capacity curve', x: xs, y: Es }, { name: 'Demand', x: [xs[0], xs[xs.length - 1]], y: [E / 1e3, E / 1e3], dash: true }],
          explain: 'Energy absorbed as deformation grows. The deformation at which the curve reaches the demand line is the predicted damage.' },
      ],
      tables: [{ title: 'Impact and resistance parameters', columns: [{ key: 'k', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
        rows: [['Object submerged weight', Wsub / 1e3, 'kN'], ['Object added mass', ma, 'kg'], ['Velocity decay length 1/β', 1 / beta, 'm'], ['Energy from the fall in air only', 0.5 * m * v0 * v0 / 1e3, 'kJ'], ...rows].map(([k, vv, u]) => ({ k, v: round(vv), u })),
        explain: 'The object properties and the resistance parameters of the target that enter the energy balance.' }],
      notes: [...(inp.target !== 'plate' ? [`Damage class by dent depth: ${dentRatio * 100 <= inp.dMinor ? 'minor' : dentRatio * 100 <= inp.dMajor ? 'major, leakage anticipated' : 'major, leakage and rupture anticipated'} (limits ${inp.dMinor} % and ${inp.dMajor} %).`] : []), inp.Edesign > 0 ? `The specified design energy of ${inp.Edesign} kJ was used instead of the computed ${fmt(Ecalc / 1e3)} kJ.` : `Computed impact energy ${fmt(Ecalc / 1e3)} kJ; ${v0 > vT ? 'the object slows down in water' : 'the object accelerates in water'} towards ${fmt(vT)} m/s.`, `"${w.name}" governs (${fmt(w.uc)}).`],
      advice,
      publish: { 'subsea.impactEnergy': Ecalc, 'subsea.dentRatio': dentRatio },
      quantities: steelDelta > 0 && !ok ? { steel_t: round(steelDelta) } : {},
    };
  },
  verify: [
    { name: 'Terminal velocity and energy: vT = √(2W′/(ρ·Cd·A)), E = ½(m + ma)vT²', source: 'Hand calculation: W′ = (1000 − 1025·0.13)·9.80665 = 8500.0 N; vT = √(2·8500.0/1025) = 4.0725 m/s; E = ½·(1000 + 133.25)·16.585 = 9.398 kJ (deep water, terminal)',
      inputs: { mObj: 1000, volObj: 0.13, Aobj: 1, CdObj: 1, CaObj: 1, hAir: 0, depth: 1000 },
      expect: [{ kpi: 'vT', value: Math.sqrt(2 * (1000 - 133.25) * 9.80665 / 1025), tol: 1e-6 }, { kpi: 'Eimp', value: 0.5 * 1133.25 * (2 * 866.75 * 9.80665 / 1025) / 1e3, tol: 1e-4 }] },
    { name: 'Free fall in air: v = √(2gh)', source: 'Galileo: h = 20 m → 19.806 m/s', inputs: { hAir: 20 }, expect: [{ kpi: 'vEntry', value: 19.8057, tol: 1e-4 }] },
    { name: 'Dent energy: E = 16·(2π/9)^½·(fy·t²/4)·(D/t)^½·D·(δ/D)^1.5', source: 'DNV-RP-F107 dent–energy relation (Wierzbicki & Suh 1988): D = 406.4 mm, t = 16 mm, fy = 355 MPa, δ/D = 5 %',
      inputs: { target: 'pipe', Dt: 406.4, tt: 16, mat: 's355' },
      expect: [{ kpi: 'Edent5', value: 16 * Math.sqrt(2 * Math.PI / 9) * (355e6 * 0.016 ** 2 / 4) * Math.sqrt(406.4 / 16) * 0.4064 * 0.05 ** 1.5 / 1e3, tol: 1e-6 }] },
  ],
};
