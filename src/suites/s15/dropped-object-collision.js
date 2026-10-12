// Accidental events during installation: a load dropped from the crane onto a deck or subsea asset, and vessel
// collision with a leg or brace. Energy demand, local denting, plastic beam response and residual capacity.
import { erf } from '../../engine/math/special.js';
import { G, DEG, need, round, fmt, kpi, check, adv, line, tbl, worst, linspace, tube, steelCo2Text, vesselCo2Text, CONFIRM } from './_common.js';

export default {
  id: 's15-dropped-object-collision', suite: 15, title: 'Dropped object and collision',
  summary: 'Computes the impact energy of a dropped lift or a colliding vessel, the dent and plastic deflection of the struck tubular member, its residual capacity, and how likely a dropped object is to hit a subsea asset.',
  covers: [16], fidelity: 'screening',
  standards: ['DNV-RP-C204 (2019) §3–4 (design against accidental loads: ship collision, dropped objects, energy dissipation)', 'NORSOK N-003 (2017) §8 / N-004 (2013) Annex A (collision energies, tubular denting and beam mechanisms)', 'DNV-RP-F107 (2019) §4–5 (dropped objects over subsea assets: excursion, dent–energy relation)', 'DNV-ST-N001 (2021) Sec. 16 (lifting over live assets, dropped-object protection)'],
  models: ['Conservation of energy and momentum (impact energy with added mass)', 'Free fall in air and terminal velocity in water', 'Local denting of a tube (Wierzbicki & Suh / Ellinas & Walker form used in DNV-RP-F107)', 'Rigid-plastic beam mechanism', 'Dented-section residual capacity (dented arc ineffective)', 'Normal distribution of lateral excursion through the water column'],
  limits: 'Energy method with rigid-plastic resistance: the striking object is rigid unless an energy share is entered, membrane action, strain-rate effects, fracture and joint failure are not modelled, and the dent expression, added-mass factors, excursion angles and design collision energies are recalled values exposed as inputs to be confirmed. Explicit nonlinear finite-element impact analysis and a quantitative risk assessment of drop frequencies are outside the module.',
  publishes: ['install.impactEnergy', 'install.dentRatio'],
  searchTerms: 'dropped object vessel collision impact energy tubular dent residual strength DNV-RP-C204',
  inputs: [
    { key: 'scenario', label: 'Accidental scenario', type: 'select', default: 'drop', group: 'Scenario', role: 'operation', options: [{ value: 'drop', label: 'Load dropped from the crane' }, { value: 'collision', label: 'Vessel collision with the structure' }] },
    { key: 'mObj', label: 'Dropped object mass', unit: 't', type: 'number', default: 5, min: 0.05, max: 5000, group: 'Dropped object', role: 'load', uncertain: { dist: 'lognormal', cov: 0.2 }, showIf: { key: 'scenario', in: ['drop'] } },
    { key: 'hAir', label: 'Drop height in air', unit: 'm', type: 'number', default: 8, min: 0, max: 150, group: 'Dropped object', role: 'ic', showIf: { key: 'scenario', in: ['drop'] }, help: 'Initial condition: height of the load above the deck it strikes, or above the sea surface for a subsea target.' },
    { key: 'subsea', label: 'Target is below water', type: 'bool', default: false, group: 'Dropped object', role: 'bc', showIf: { key: 'scenario', in: ['drop'] } },
    { key: 'depth', label: 'Water depth above the target', unit: 'm', type: 'number', default: 100, min: 1, max: 3500, group: 'Dropped object', role: 'environment', bind: 'site.depth', showIf: { key: 'scenario', in: ['drop'] } },
    { key: 'volObj', label: 'Object volume', unit: 'm³', type: 'number', default: 1.1, min: 0, max: 5000, group: 'Dropped object', role: 'geometry', showIf: { key: 'scenario', in: ['drop'] } },
    { key: 'Aobj', label: 'Projected area in the fall direction', unit: 'm²', type: 'number', default: 4, min: 0.01, max: 2000, group: 'Dropped object', role: 'geometry', showIf: { key: 'scenario', in: ['drop'] } },
    { key: 'CdObj', label: 'Drag coefficient in water', unit: '-', type: 'number', default: 1.0, min: 0.2, max: 3, group: 'Dropped object', role: 'load', showIf: { key: 'scenario', in: ['drop'] }, help: `Slender 0.7–1.5, box-shaped 1.2–1.3, complex 0.6–2.0. ${CONFIRM}` },
    { key: 'CaObj', label: 'Added-mass coefficient in water', unit: '-', type: 'number', default: 1.0, min: 0, max: 2.5, group: 'Dropped object', role: 'load', showIf: { key: 'scenario', in: ['drop'] } },
    { key: 'devDeg', label: 'Excursion angle (one standard deviation)', unit: 'deg', type: 'number', default: 5, min: 0.5, max: 30, group: 'Dropped object', role: 'criteria', showIf: { key: 'scenario', in: ['drop'] }, help: `Flat or long objects up to 15°, box-shaped 3–5°, heavy compact about 2°. ${CONFIRM}` },
    { key: 'rTarget', label: 'Radius of the subsea asset in plan', unit: 'm', type: 'number', default: 8, min: 0.1, max: 200, group: 'Dropped object', role: 'geometry', showIf: { key: 'scenario', in: ['drop'] } },
    { key: 'mVessel', label: 'Vessel displacement', unit: 't', type: 'number', default: 5000, min: 50, max: 2e5, group: 'Collision', role: 'load', showIf: { key: 'scenario', in: ['collision'] }, help: `Supply-vessel class for the design collision; 5000 t is the customary minimum, larger for modern field vessels. ${CONFIRM}` },
    { key: 'vVessel', label: 'Impact velocity', unit: 'm/s', type: 'number', default: 2, min: 0.05, max: 10, group: 'Collision', role: 'ic', uncertain: { dist: 'lognormal', cov: 0.3 }, showIf: { key: 'scenario', in: ['collision'] }, help: `2 m/s for a drifting supply vessel in the accidental limit state; 0.5 m/s for an operational bump. ${CONFIRM}` },
    { key: 'aVessel', label: 'Hydrodynamic added mass fraction', unit: '-', type: 'number', default: 0.4, min: 0, max: 1, group: 'Collision', role: 'load', showIf: { key: 'scenario', in: ['collision'] }, help: `0.4 for sideways (broadside) impact, 0.1 for bow or stern impact. ${CONFIRM}` },
    { key: 'eta', label: 'Share of the energy absorbed by the struck member', unit: '-', type: 'number', default: 1.0, min: 0.05, max: 1, group: 'Scenario', role: 'criteria', help: 'Less than 1 when the striking object or vessel deforms as well (shared-energy design); 1.0 is conservative for the member.' },
    { key: 'mat', label: 'Steel of the struck member', type: 'select', source: 'materials:steel', default: 's355', group: 'Struck member', role: 'material' },
    { key: 'D', label: 'Member outer diameter', unit: 'mm', type: 'number', default: 1200, min: 100, max: 6000, group: 'Struck member', role: 'geometry', design: true },
    { key: 't', label: 'Wall thickness', unit: 'mm', type: 'number', default: 40, min: 5, max: 150, group: 'Struck member', role: 'geometry', design: true },
    { key: 'Lm', label: 'Member span between joints', unit: 'm', type: 'number', default: 18, min: 1, max: 100, group: 'Struck member', role: 'geometry' },
    { key: 'fix', label: 'End restraint', type: 'select', default: 'fixed', group: 'Struck member', role: 'bc', options: [{ value: 'fixed', label: 'Fixed ends (collapse load 8·Mp/L)' }, { value: 'pinned', label: 'Pinned ends (collapse load 4·Mp/L)' }], help: 'Boundary condition of the struck member: joint flexibility lies between the two.' },
    { key: 'dentAllow', label: 'Allowable dent depth / diameter', unit: '%', type: 'number', default: 15, min: 1, max: 50, group: 'Criteria', role: 'criteria', help: `Structural members tolerate larger dents than pressure-containing pipe (5–10 %) if the residual capacity is checked. ${CONFIRM}` },
    { key: 'duct', label: 'Allowable plastic deflection / span', unit: '%', type: 'number', default: 4, min: 0.2, max: 20, group: 'Criteria', role: 'criteria', help: 'Limits the hinge rotation so that the strain stays ductile; depends on section class and material toughness.' },
    { key: 'resReq', label: 'Required residual axial capacity after impact', unit: '%', type: 'number', default: 60, min: 0, max: 100, group: 'Criteria', role: 'criteria', help: 'What the damaged member must still carry in the post-accident load case (from the structural analysis).' },
  ],
  compute(inp, ctx) {
    // t → kg, mm → m, % → fraction.
    const steel = ctx.material(inp.mat); need(steel, `Steel '${inp.mat}' is not in the materials library.`);
    const fy = steel.fy, rho = 1025, D = inp.D / 1e3, t = inp.t / 1e3;
    need(t < D / 2, 'Wall thickness must be smaller than the radius.');
    const drop = inp.scenario === 'drop';
    let E0, v, vEntry = 0, vT = 0, pHit = 1, sig = 0, mEff;
    const zs = [], vs = [];
    if (drop) {
      const m = inp.mObj * 1e3; vEntry = Math.sqrt(2 * G * inp.hAir);
      if (inp.subsea) {
        const Wsub = (m - rho * inp.volObj) * G, ma = rho * inp.CaObj * inp.volObj;
        need(Wsub > 0, 'The object floats (mass below displaced water): it does not reach the subsea target.');
        vT = Math.sqrt(2 * Wsub / (rho * inp.CdObj * inp.Aobj));
        // (m + ma)·v·dv/dz = W′ − ½ρ·Cd·A·v² → v² relaxes exponentially from the entry value to the terminal value.
        const beta = rho * inp.CdObj * inp.Aobj / (m + ma), vAt = (z) => Math.sqrt(vT * vT + (vEntry * vEntry - vT * vT) * Math.exp(-beta * z));
        v = vAt(inp.depth); mEff = m + ma;
        sig = inp.depth * Math.tan(inp.devDeg * DEG); pHit = erf(inp.rTarget / (sig * Math.SQRT2));
        for (const z of linspace(0, inp.depth, 30)) { zs.push(-round(z)); vs.push(round(vAt(z))); }
      } else { v = vEntry; mEff = m; }
      E0 = 0.5 * mEff * v * v;
    } else {
      mEff = inp.mVessel * 1e3 * (1 + inp.aVessel); v = inp.vVessel; E0 = 0.5 * mEff * v * v;
    }
    const E = E0 * inp.eta;
    // Tube: dent energy E = 16·(2π/9)^½·mp·(D/t)^½·D·(δ/D)^1.5, contact force F = dE/dδ.
    const mp = fy * t * t / 4, cD = 16 * Math.sqrt(2 * Math.PI / 9) * mp * Math.sqrt(D / t) * D;
    const Ed = (r) => cD * r ** 1.5, Fd = (r) => 1.5 * cD * Math.sqrt(r) / D;
    const sec = tube(D, t), Mp = fy * sec.Wp, Pc = (inp.fix === 'fixed' ? 8 : 4) * Mp / inp.Lm;
    // Denting proceeds until its contact force reaches the beam collapse load; the rest goes into the beam mechanism.
    const rY = Math.min((Pc * D / (1.5 * cD)) ** 2, 1);
    let dent, defl = 0, Ebeam = 0;
    if (E <= Ed(rY)) dent = (E / cD) ** (2 / 3); else { dent = rY; Ebeam = E - Ed(rY); defl = Ebeam / Pc; }
    const dAllow = inp.dentAllow / 100, wAllow = inp.duct / 100 * inp.Lm;
    const cap = Ed(Math.min(dAllow, rY)) + (dAllow >= rY ? Pc * wAllow : 0);
    const Fpeak = Math.min(Fd(dent), Pc);
    // Residual capacity with the flattened arc ineffective: half-angle α from cos α = 1 − 2δ/D.
    const alpha = Math.acos(Math.max(1 - 2 * Math.min(dent, 1), -1)), resN = 1 - alpha / Math.PI;
    const checks = [
      check('Impact energy against absorption capacity', E / 1e3, cap / 1e3, 'kJ', 'DNV-RP-C204 (2019) §3 (energy dissipation); NORSOK N-004 (2013) A.3', 'Energy the member can absorb by denting and plastic bending within the allowable deformation. Above it the member ruptures or is pushed through and the load path is lost.'),
      check('Dent depth', dent * 100, inp.dentAllow, '% of diameter', 'DNV-RP-F107 (2019) §4 dent–energy relation; NORSOK N-004 A.3.10 (dented tubulars)', 'Local flattening of the wall under the contact point. Deep dents destroy the axial and bending capacity of the section.'),
      check('Plastic deflection of the member', defl, wAllow, 'm', 'DNV-RP-C204 (2019) §3.10 (ductility limits for plastic hinges)', 'Global bending after the hinge mechanism forms. The rotation must stay within the ductility of the section or the wall tears.'),
      check('Residual axial capacity after impact', inp.resReq, Math.max(resN * 100, 1e-6), '% (required / remaining)', 'NORSOK N-004 (2013) A.3.10.3 (resistance of dented tubulars); DNV-RP-C204 §3 (post-damage check)', 'What is left of the member for the post-accident load case. If it is below the demand, the frame must redistribute or the member be repaired before the next storm.'),
    ];
    const gov = worst(checks), ok = gov.uc <= 1;
    const tReq = t * Math.max((E / Math.max(cap, 1e-9)), 1e-9) ** 0.5 * 1e3; // capacity grows roughly with t² (dent ∝ t^2.5, beam ∝ t)
    const steelAdd = Math.max(tReq / 1e3 - t, 0) * Math.PI * D * inp.Lm * 7.85;
    const rs = linspace(0, Math.max(0.3, 1.3 * dent), 31);
    const cls = dent * 100 <= 5 ? 'minor' : dent * 100 <= 10 ? 'moderate' : 'major';
    return {
      kpis: [
        kpi('E0', 'Impact energy', E0 / 1e3, 'kJ', drop ? 'Kinetic energy of the object when it strikes: m·g·h in air, or ½·(m + added mass)·v² at depth for a subsea target.' : 'Kinetic energy of the vessel including hydrodynamic added mass: ½·(1 + a)·m·v².'),
        kpi('Eabs', 'Energy to be absorbed by the member', E / 1e3, 'kJ', 'Impact energy times the share assigned to the struck member.', cap / 1e3, 'max'),
        kpi('v', 'Impact velocity', v, 'm/s', drop ? 'Velocity at the target: free fall in air, relaxed towards terminal velocity in water.' : 'Design drift or approach velocity of the vessel.'),
        kpi('vT', 'Terminal velocity in water', vT, 'm/s', 'Velocity at which submerged weight balances drag; zero when the target is above water or for a collision.'),
        kpi('pHit', 'Share of drops landing on the subsea asset', pHit * 100, '%', 'Probability that a drop from directly above lands within the asset radius, for a normally distributed excursion; 100 % for a target in air.'),
        kpi('dent', 'Dent depth / diameter', dent * 100, '%', 'Permanent local dent from the absorbed energy.', inp.dentAllow, 'max'),
        kpi('defl', 'Plastic deflection of the member', defl, 'm', 'Global permanent deflection once the denting force reaches the beam collapse load; zero if denting alone absorbs the energy.', wAllow, 'max'),
        kpi('Fpeak', 'Peak contact force', Fpeak / 1e6, 'MN', 'Largest force between object and member: the denting force, capped by the beam collapse load.'),
        kpi('Pc', 'Beam collapse load of the member', Pc / 1e6, 'MN', 'Lateral point load at mid-span that forms a plastic mechanism (8·Mp/L fixed, 4·Mp/L pinned).'),
        kpi('cap', 'Energy absorption capacity', cap / 1e3, 'kJ', 'Dent energy to the allowable dent plus beam plastic work to the allowable deflection.'),
        kpi('resN', 'Residual axial capacity', resN * 100, '%', 'Squash-load fraction remaining with the dented arc treated as ineffective.', inp.resReq, 'min'),
        kpi('tReq', 'Wall thickness to absorb the energy within limits', tReq, 'mm', 'Indicative thickness at which capacity equals demand (capacity taken to grow with t²).'),
      ],
      checks,
      plots: [
        line('Absorbed energy against dent depth', 'Dent depth / diameter', '%', 'Energy', 'kJ', [{ name: 'Dent energy', x: rs.map((r) => round(r * 100)), y: rs.map((r) => round(Ed(Math.min(r, rY)) / 1e3)) }, { name: 'Demand', x: [0, round(rs[rs.length - 1] * 100)], y: [E / 1e3, E / 1e3], dash: true }, { name: 'Capacity incl. beam mechanism', x: [0, round(rs[rs.length - 1] * 100)], y: [cap / 1e3, cap / 1e3], dash: true }],
          'Energy absorbed by local denting as the dent deepens; it flattens once the member starts to bend as a beam. The demand line shows where the impact ends up on this curve.'),
        ...(drop && inp.subsea ? [line('Velocity of the falling object against depth', 'Velocity', 'm/s', 'Depth', 'm', [{ name: 'Velocity', x: vs, y: zs }], 'The object slows (or accelerates) from its entry velocity to terminal velocity within the first tens of metres; below that depth does not change the impact energy.')] : []),
      ],
      tables: [tbl('Impact and resistance parameters', [['k', 'Quantity'], ['v', 'Value'], ['u', 'Unit']], [
        ['Effective mass incl. added mass', mEff / 1e3, 't'], ['Wall plastic moment mp', mp / 1e3, 'kN·m/m'], ['Section plastic moment Mp', Mp / 1e6, 'MN·m'], ['Dent ratio at onset of the beam mechanism', rY * 100, '%'], ['Energy absorbed by denting', Math.min(E, Ed(rY)) / 1e3, 'kJ'], ['Energy absorbed by beam bending', Ebeam / 1e3, 'kJ'], ['Water-entry velocity', vEntry, 'm/s'], ['Lateral excursion (one sigma)', sig, 'm'],
      ].map(([k, vv, u]) => ({ k, v: round(vv), u })), 'How the impact energy is built up and where it is absorbed in the member.')],
      notes: [`${drop ? `Dropped ${inp.mObj} t from ${inp.hAir} m${inp.subsea ? ` onto a target at ${inp.depth} m depth` : ''}` : `${inp.mVessel} t vessel at ${inp.vVessel} m/s with ${inp.aVessel} added mass`}: ${fmt(E0 / 1e3)} kJ.`, `Damage class by dent depth: ${cls}; residual axial capacity ${fmt(resN * 100)} %.`, `"${gov.name}" governs at ${fmt(gov.uc)}.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `The ${inp.D} × ${inp.t} mm member absorbs ${fmt(E / 1e3)} kJ with a ${fmt(dent * 100)} % dent and ${fmt(defl * 1e3)} mm deflection; ${fmt(resN * 100)} % of its axial capacity remains.`,
          `Run the post-accident structural case with this member at ${fmt(resN * 100)} % capacity, check the joints for the peak force of ${fmt(Fpeak / 1e6)} MN, and confirm the design energy against the governing edition and the field-specific vessel sizes.`, 'DNV-RP-C204 (2019) §3–4; NORSOK N-004 (2013) Annex A')
          : adv('engineering', 'fail', `"${gov.name}" is not met (${fmt(gov.uc)}): ${fmt(E / 1e3)} kJ against ${fmt(cap / 1e3)} kJ capacity, dent ${fmt(dent * 100)} %, residual capacity ${fmt(resN * 100)} %.`,
            drop ? `Restrict the lifting zone or lift height (the energy scales with height), add a dropped-object protection frame, or thicken the member to about ${fmt(Math.ceil(tReq))} mm.` : `Thicken the member to about ${fmt(Math.ceil(tReq))} mm, add fendering or a sacrificial boat-impact member, or demonstrate shared energy absorption with the vessel (strength design).`,
            'DNV-RP-C204 (2019) §3–4; NORSOK N-004 (2013) Annex A'),
        adv('economic', 'info', `${drop ? `${fmt(pHit * 100)} % of drops from above reach the target` : 'The design collision is an accidental-limit-state event'}; the damage class is ${cls} with ${fmt(resN * 100)} % residual capacity.`,
          'Weigh protection steel against operational controls: lifting-zone restrictions, vessel approach procedures and heading control cut the frequency at no capital cost, while a repair after a major dent means a clamp or member replacement campaign and possibly restricted operation until it is done.',
          'Risk-based decision per ISO 17776 / DNV-RP-C204 §2 (frequency × consequence); repair campaign priced in vessel days (AACE RP 18R-97)',
          steelAdd > 0 && !ok ? `Thickening this member adds about ${steelCo2Text(steelAdd)}` : `A subsea or splash-zone repair campaign of about ten vessel days would cost roughly ${vesselCo2Text(10, 30)} Avoiding the event is the lowest-carbon option.`),
      ],
      publish: { 'install.impactEnergy': E0, 'install.dentRatio': dent },
      quantities: !ok && steelAdd > 0 ? { steel_t: round(steelAdd) } : {},
    };
  },
  verify: [
    { name: 'Drop in air: E = m·g·h, v = √(2·g·h)', source: 'Energy conservation: 8 t from 15 m → 8000·9.80665·15 = 1176.8 kJ, v = 17.153 m/s', inputs: { mObj: 8, hAir: 15 }, expect: [{ kpi: 'E0', value: 8000 * 9.80665 * 15 / 1e3, tol: 1e-9 }, { kpi: 'v', value: Math.sqrt(2 * 9.80665 * 15), tol: 1e-9 }] },
    { name: 'Supply-vessel collision energies: 14 MJ sideways, 11 MJ bow/stern', source: 'NORSOK N-003 / DNV-RP-C204 design values: 5000 t at 2 m/s with added mass 0.4 → ½·1.4·5e6·4 = 14.0 MJ; with 0.1 → 11.0 MJ',
      inputs: { scenario: 'collision' }, expect: [{ kpi: 'E0', value: 14000, tol: 1e-9 }] },
    { name: 'Dent energy relation: E = 16·(2π/9)^½·(fy·t²/4)·(D/t)^½·D·(δ/D)^1.5 inverted for the dent', source: 'DNV-RP-F107 (2019) dent–energy relation (Wierzbicki & Suh 1988): D = 1200 mm, t = 40 mm, fy = 355 MPa, E = 200 kJ, short fixed span so that no beam mechanism forms → δ/D = (E/c)^(2/3) with c = 16·0.83555·142 000·5.4772·1.2 = 12 478 kJ → 6.36 %',
      inputs: { mObj: 200 / (9.80665 * 15), hAir: 15, Lm: 4, fix: 'fixed' }, expect: [{ kpi: 'dent', value: 100 * (200e3 / (16 * Math.sqrt(2 * Math.PI / 9) * (355e6 * 0.04 ** 2 / 4) * Math.sqrt(30) * 1.2)) ** (2 / 3), tol: 1e-6 }] },
  ],
};
