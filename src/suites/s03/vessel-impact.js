// Accidental impact and structural damage: service-vessel collision with a monopile or a jacket member —
// kinetic energy, sharing of strain energy between vessel and structure, local denting, global bending response
// and the residual strength of the damaged section.
import { collisionEnergy, tubeDentResistance, tubeDentDepth, beamMechanismResistance, pileImpact, sharedEnergy, classifyCollisionDesign } from '../../engine/structural/impact.js';
import { need, run, fmt, round, adv, check, worst, sev, steelCo2Text, vesselCo2Text, tube, steelInputs, steelOf, LIMIT_CODE, CONFIRM } from './_common.js';

/** Remaining plastic moment of a dented tube relative to the intact value (DNV-RP-C204 §3.10.2; Ellinas 1984). */
function residualMoment(dentRatio) {
  const x = Math.min(Math.max(dentRatio, 0), 0.5);
  const th = Math.acos(1 - 2 * x);
  return Math.max(Math.cos(th / 2) - 0.5 * Math.sin(th), 0);
}

export default {
  id: 's03-vessel-impact', suite: 3, title: 'Vessel impact and structural damage',
  summary: 'Operational and accidental collision of a service vessel with the foundation: impact energy, contact force, energy shared between vessel and structure, dent depth, global response and the bending capacity left in the damaged section.',
  covers: [15], fidelity: 'analytical',
  standards: ['DNV-ST-0126 (2021) §4.3.6 (ship impact: operational and accidental)', 'DNV-RP-C204 (2019) §3 (ship collisions: energy, denting, residual strength)', 'IEC 61400-3-1:2019 §7.4 DLC 8.x and accidental situations', 'NORSOK N-004 (2013) Annex A.3', 'ISO 19902:2020 §10 (accidental situations)'],
  models: ['Kinetic energy with hydrodynamic added mass', 'Strain-energy sharing from force–deformation curves in series', 'Local denting of tubular members (DNV-RP-C204 resistance curve)', 'Elastic cantilever and plastic-hinge energy absorption of a pile', 'Three-hinge beam mechanism with membrane action for braced members', 'Residual plastic moment of a dented tube', 'Force and moment equilibrium'],
  limits: `Quasi-static energy method with a rigid–plastic or elastic structure and a bilinear vessel curve; strain-rate effects, fracture, fender and boat-landing energy absorption, the dynamic response of the turbine to the impulse and soil nonlinearity at the fixity point are not modelled. ${LIMIT_CODE} Explicit nonlinear finite-element collision analysis is outside the module.`,
  searchTerms: 'service vessel collision monopile dent energy absorption DNV-RP-C204 offshore wind',
  inputs: [
    { key: 'member', label: 'Struck member', type: 'select', default: 'monopile', group: 'Structure', role: 'geometry', options: [{ value: 'monopile', label: 'Monopile (cantilever from the soil)' }, { value: 'span', label: 'Jacket leg or brace (span between joints)' }] },
    { key: 'D', label: 'Outer diameter at the impact level', unit: 'm', type: 'number', default: 6, min: 0.3, max: 14, group: 'Structure', role: 'geometry', design: true },
    { key: 't', label: 'Wall thickness at the impact level', unit: 'mm', type: 'number', default: 70, min: 8, max: 160, group: 'Structure', role: 'geometry', design: true },
    ...steelInputs(),
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 25, min: 5, max: 80, group: 'Structure', role: 'environment', bind: 'site.depth', showIf: { key: 'member', in: ['monopile'] } },
    { key: 'zImp', label: 'Impact elevation above still water', unit: 'm', type: 'number', default: 3, min: -5, max: 20, group: 'Structure', role: 'bc', showIf: { key: 'member', in: ['monopile'] }, help: 'Height of the vessel fender or bow contact, taking tide and wave motion into account.' },
    { key: 'fixD', label: 'Apparent fixity depth / diameter', unit: '-', type: 'number', default: 3.5, min: 0, max: 8, group: 'Structure', role: 'bc', showIf: { key: 'member', in: ['monopile'] }, help: 'Depth below the mudline at which the pile is taken as clamped: 3–4 D in dense sand or stiff clay, 6–8 D in soft clay.' },
    { key: 'span', label: 'Member span between joints', unit: 'm', type: 'number', default: 18, min: 2, max: 80, group: 'Structure', role: 'geometry', showIf: { key: 'member', in: ['span'] } },
    { key: 'ends', label: 'End fixity of the member', type: 'select', default: 'clamped', group: 'Structure', role: 'bc', showIf: { key: 'member', in: ['span'] }, options: [{ value: 'clamped', label: 'Clamped (stiff joints)' }, { value: 'pinned', label: 'Pinned' }] },
    { key: 'dent0', label: 'Existing dent depth', unit: 'mm', type: 'number', default: 0, min: 0, max: 3000, group: 'Structure', role: 'ic', help: 'Damage already present from inspection; the residual strength is evaluated for the larger of this and the computed dent.' },
    { key: 'mass', label: 'Vessel displacement', unit: 't', type: 'number', default: 3500, min: 20, max: 30000, group: 'Vessel', role: 'bc', help: 'Largest service vessel authorised to approach: about 150–400 t for crew-transfer vessels, 3000–8000 t for service operation vessels.' },
    { key: 'vOp', label: 'Operational impact speed (bow or stern)', unit: 'm/s', type: 'number', default: 0.5, min: 0.05, max: 3, group: 'Vessel', role: 'bc', help: 'DNV-ST-0126: not less than 0.5 m/s for boat landing in normal operation.' },
    { key: 'vAcc', label: 'Accidental impact speed (drifting sideways)', unit: 'm/s', type: 'number', default: 2, min: 0.2, max: 6, group: 'Vessel', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'DNV-ST-0126: not less than 2 m/s for an unintended drifting collision.' },
    { key: 'cmBow', label: 'Added-mass factor, surge', unit: '-', type: 'number', default: 1.1, min: 1, max: 1.5, group: 'Vessel', role: 'bc' },
    { key: 'cmSide', label: 'Added-mass factor, sway', unit: '-', type: 'number', default: 1.4, min: 1, max: 2, group: 'Vessel', role: 'bc' },
    { key: 'kShip', label: 'Vessel crushing stiffness', unit: 'MN/m', type: 'number', default: 12, min: 0.5, max: 500, group: 'Vessel', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.4 }, help: 'Initial slope of the vessel force–indentation curve (hull side or fendered bow). Stiffer vessels put more energy into the structure.' },
    { key: 'fShip', label: 'Vessel crushing plateau force', unit: 'MN', type: 'number', default: 15, min: 0.5, max: 200, group: 'Vessel', role: 'bc', help: 'Force at which the vessel structure crushes at roughly constant load (DNV-RP-C204 Fig. 3-4 gives curves for supply vessels).' },
    { key: 'dShip', label: 'Vessel deformation limit', unit: 'm', type: 'number', default: 2.5, min: 0.2, max: 10, group: 'Vessel', role: 'bc' },
    { key: 'dentOp', label: 'Allowable dent / diameter, operational', unit: '-', type: 'number', default: 0.01, min: 0.001, max: 0.1, group: 'Acceptance criteria', role: 'criteria', help: 'Operational impact must leave no damage needing repair.' },
    { key: 'dentAcc', label: 'Allowable dent / diameter, accidental', unit: '-', type: 'number', default: 0.1, min: 0.01, max: 0.5, group: 'Acceptance criteria', role: 'criteria', help: `Limit of the denting formula and of the residual-strength expression. ${CONFIRM}` },
    { key: 'rotCap', label: 'Plastic hinge rotation capacity', unit: 'rad', type: 'number', default: 0.01, min: 0, max: 0.1, group: 'Acceptance criteria', role: 'criteria', help: 'Small for slender (class 3–4) shells such as monopiles; up to 0.05 for compact jacket members.' },
    { key: 'Mpost', label: 'Bending moment to survive after the accident', unit: 'MN·m', type: 'number', default: 120, min: 0, max: 5000, group: 'Acceptance criteria', role: 'load', help: 'Post-damage environmental moment at the impact section (for example the 1-year condition), with load factor 1.0.' },
  ],
  compute(inp, ctx) {
    const steel = steelOf(inp, ctx);
    const D = inp.D, t = inp.t / 1e3, sec = tube(D, t);
    const mono = inp.member === 'monopile';
    const Mp = steel.fy * sec.Zp, EI = steel.E * sec.I;
    const m = inp.mass * 1e3;
    const eOp = collisionEnergy({ mass: m, velocity: inp.vOp, addedMassFactor: inp.cmBow }).energy;
    const eAcc = collisionEnergy({ mass: m, velocity: inp.vAcc, addedMassFactor: inp.cmSide }).energy;
    // Vessel curve: linear to the plateau, then constant to the deformation limit.
    const dPl = inp.fShip / inp.kShip;
    need(inp.dShip > dPl, 'The vessel deformation limit must exceed the indentation at which the plateau force is reached.');
    const shipCurve = { d: [0, dPl, inp.dShip], F: [0, inp.fShip * 1e6, inp.fShip * 1e6 * 1.0001] };
    // Structure curve: local dent and global deflection add at equal force (springs in series).
    const h = mono ? inp.zImp + inp.depth + inp.fixD * D : inp.span;
    need(h > 2 * D || !mono, 'The impact point is too close to the fixity point for beam theory: check the impact elevation and fixity depth.');
    const Fglob = mono ? Mp / h : beamMechanismResistance({ section: 'tube', D, t, fy: steel.fy, L: inp.span, ends: inp.ends, w: 0 }).R0;
    const kGlob = mono ? (3 * EI) / h ** 3 : (inp.ends === 'clamped' ? 192 : 48) * EI / inp.span ** 3;
    const dentAt = (F) => tubeDentDepth({ D, t, fy: steel.fy, force: F }).dent;
    const Fdmax = tubeDentResistance({ D, t, fy: steel.fy, dent: 0.3 * D }).R;
    const Fmax = Math.min(Fglob, Fdmax);
    const sd = [0], sF = [0];
    for (let i = 1; i <= 24; i++) { const F = (Fmax * i) / 24; sd.push(dentAt(F) + F / kGlob); sF.push(F); }
    // beyond the global collapse load: plastic hinge rotation (monopile) or mechanism deflection (span) at constant force
    const dPlast = mono ? inp.rotCap * h : 0.5 * D;
    if (Fglob <= Fdmax && dPlast > 0) { sd.push(sd[sd.length - 1] + dPlast); sF.push(Fmax * 1.0001); }
    const structureCurve = { d: sd, F: sF };
    const share = (E) => run(() => sharedEnergy({ shipCurve, structureCurve, energy: E }), 'Energy sharing');
    const op = share(eOp), acc = share(eAcc);
    const dentOf = (r) => dentAt(Math.min(r.force, Fdmax));
    const dOp = dentOf(op), dAcc = dentOf(acc);
    const globE = (r, dnt) => Math.max(r.structureEnergy - tubeDentResistance({ D, t, fy: steel.fy, dent: dnt }).energy, 0);
    const pi = mono ? pileImpact({ energy: globE(acc, dAcc), Mp, EI, height: h, rotationCapacity: inp.rotCap }) : null;
    const cat = classifyCollisionDesign({ structureEnergy: acc.structureEnergy, shipEnergy: acc.shipEnergy });
    const dentRes = Math.max(dAcc, inp.dent0 / 1e3);
    const res = residualMoment(dentRes / D);
    const Mres = res * Mp / 1.1;
    const lever = mono ? inp.zImp + inp.depth : inp.span / 4;
    const Mimp = acc.force * lever, My = steel.fy * sec.W;
    const checks = [
      check('Operational impact: dent depth', dOp / D, inp.dentOp, 'dent / D', 'DNV-ST-0126 (2021) §4.3.6 (boat impact in normal operation: no permanent damage)', 'Routine boat landings must not leave damage. A visible dent after an operational impact means the fendering or the wall thickness is inadequate.'),
      check('Operational impact: contact force against first-yield resistance of the member', op.force / 1e6, (mono ? My / h : Fglob * 0.75) / 1e6, 'MN', 'DNV-ST-0126 (2021) §4.3.6; elastic response required', 'Keeps the global response elastic in normal operation: the force at which the extreme fibre yields at the fixity point (or the member starts forming hinges).'),
      check('Accidental impact: dent depth', dAcc / D, inp.dentAcc, 'dent / D', 'DNV-RP-C204 (2019) §3.10 (denting, validity of the residual-strength expression)', 'Limits local damage in the accidental event to a dent the section can survive. Deeper dents tear the wall or flatten the tube.'),
      check('Accidental impact: energy to be absorbed against the capacity of vessel and structure', eAcc / 1e6, acc.capacity / 1e6, 'MJ', 'DNV-RP-C204 (2019) §3.3 (strain-energy dissipation)', 'The kinetic energy must be dissipated before either the vessel or the structure reaches its deformation limit. If not, the member collapses or the hull is penetrated.'),
      check('Post-accident bending resistance of the dented section', inp.Mpost, Mres / 1e6, 'MN·m', 'DNV-RP-C204 (2019) §3.10.2; DNV-ST-0126 (2021) §4.3.6 (damaged condition, load factor 1.0)', 'After the accident the structure must still resist the environmental loads until repair. A dent removes part of the compression flange of the tube.'),
    ];
    const gov = worst(checks), ok = gov.uc <= 1;
    const dents = Array.from({ length: 31 }, (_, i) => (0.2 * D * i) / 30);
    // Wall thickness for the accidental dent limit: denting force ∝ t² approximately (Rc ∝ fy t² √(D/t)) → scan.
    let tReq = t;
    if (dAcc / D > inp.dentAcc) { for (let k = 0; k < 60 && tReq < 0.2; k++) { tReq *= 1.03; if (tubeDentDepth({ D, t: tReq, fy: steel.fy, force: Math.min(acc.force, tubeDentResistance({ D, t: tReq, fy: steel.fy, dent: 0.3 * D }).R) }).dent / D <= inp.dentAcc) break; } }
    const dM = (tube(D, tReq).A - sec.A) * steel.rho * 8 / 1e3; // thicker can over an 8 m splash-zone height
    return {
      kpis: [
        { key: 'eOp', label: 'Operational impact energy', unit: 'MJ', value: eOp / 1e6, explain: 'Kinetic energy of the vessel in a normal boat landing, including added mass in surge.' },
        { key: 'eAcc', label: 'Accidental impact energy', unit: 'MJ', value: eAcc / 1e6, explain: 'Kinetic energy of the drifting vessel striking sideways, including added mass in sway. This energy must be dissipated by deformation of vessel and structure.' },
        { key: 'forceAcc', label: 'Contact force, accidental', unit: 'MN', value: acc.force / 1e6, explain: 'Force between vessel and structure at the end of the accidental impact, from the two force–deformation curves in series.' },
        { key: 'dentAcc', label: 'Dent depth, accidental', unit: 'mm', value: dAcc * 1e3, explain: `Local indentation of the tube wall at the contact force (${fmt(dAcc / D * 100)} % of the diameter).` },
        { key: 'dentOp', label: 'Dent depth, operational', unit: 'mm', value: dOp * 1e3, explain: 'Indentation in a routine landing; it should be negligible.' },
        { key: 'shareStruct', label: 'Energy absorbed by the structure, accidental', unit: '%', value: acc.structureFraction * 100, explain: `Share of the strain energy taken by the foundation (${cat.category} design in DNV-RP-C204 terms). A stiff structure pushes the energy into the vessel.` },
        { key: 'deflAcc', label: 'Deflection of the member at the impact point', unit: 'mm', value: Math.max(acc.structureDeformation - dAcc, 0) * 1e3, explain: 'Global bending deflection of the pile or member at the contact point, excluding the local dent.' },
        { key: 'Mimp', label: mono ? 'Mudline moment from the contact force' : 'Member bending moment from the contact force', unit: 'MN·m', value: Mimp / 1e6, explain: mono ? 'Contact force times its height above the mudline: compare with the environmental design moment.' : 'Midspan bending moment of the struck member under the contact force (fixed-ended estimate).' },
        { key: 'resM', label: 'Residual bending capacity of the dented section', unit: '%', value: res * 100, explain: 'Plastic moment of the dented tube relative to the intact tube, for the larger of the computed and the existing dent.' },
        { key: 'ucEnergy', label: 'Energy utilisation, accidental', unit: '-', value: eAcc / Math.max(acc.capacity, 1), limit: 1, sense: 'max', explain: 'Impact energy over the energy the two bodies can absorb before one reaches its deformation limit.' },
        { key: 'Fcollapse', label: 'Global collapse force of the member', unit: 'MN', value: Fglob / 1e6, explain: mono ? 'Horizontal force at the impact level that forms a plastic hinge at the fixity point (Mp/h).' : 'Force that turns the member into a three-hinge mechanism (bending only, before membrane action).' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Force–deformation curves of vessel and structure', x: { label: 'Deformation', unit: 'm' }, y: { label: 'Force', unit: 'MN' },
          series: [{ name: 'Vessel', x: shipCurve.d, y: shipCurve.F.map((v) => v / 1e6) }, { name: 'Structure (dent + deflection)', x: sd, y: sF.map((v) => v / 1e6) }, { name: 'Accidental impact', x: [acc.shipDeformation, acc.structureDeformation], y: [acc.force / 1e6, acc.force / 1e6] }, { name: 'Operational impact', x: [op.shipDeformation, op.structureDeformation], y: [op.force / 1e6, op.force / 1e6] }],
          explain: 'The two bodies carry the same force; the areas under the curves up to that force are the energies each absorbs. The softer body takes most of the energy.' },
        { type: 'line', title: 'Denting resistance and residual capacity', x: { label: 'Dent depth / diameter', unit: '-' }, y: { label: 'Force (MN) / residual moment (%)', unit: '' },
          series: [{ name: 'Denting force, MN', x: dents.map((d) => d / D), y: dents.map((d) => tubeDentResistance({ D, t, fy: steel.fy, dent: d }).R / 1e6) }, { name: 'Residual plastic moment, %', x: dents.map((d) => d / D), y: dents.map((d) => residualMoment(d / D) * 100), dash: true }, { name: 'Accidental dent', x: [dAcc / D], y: [Math.min(acc.force, Fdmax) / 1e6] }],
          explain: 'Force needed to dent the wall against dent depth, and the bending capacity that remains. The marked point is the accidental impact; its position on the dashed curve is the strength left for the storm that follows.' },
      ],
      tables: [
        { title: 'Impact cases', columns: [{ key: 'k', label: 'Quantity' }, { key: 'o', label: 'Operational' }, { key: 'a', label: 'Accidental' }, { key: 'u', label: 'Unit' }],
          rows: [['Speed', inp.vOp, inp.vAcc, 'm/s'], ['Kinetic energy', eOp / 1e6, eAcc / 1e6, 'MJ'], ['Contact force', op.force / 1e6, acc.force / 1e6, 'MN'], ['Vessel indentation', op.shipDeformation, acc.shipDeformation, 'm'], ['Energy in the vessel', op.shipEnergy / 1e6, acc.shipEnergy / 1e6, 'MJ'], ['Energy in the structure', op.structureEnergy / 1e6, acc.structureEnergy / 1e6, 'MJ'], ['Dent depth', dOp * 1e3, dAcc * 1e3, 'mm'], ['Dent / diameter', dOp / D, dAcc / D, '-']]
            .map(([k, o, a, u]) => ({ k, o: round(o), a: round(a), u })),
          explain: 'Energy balance of both design impacts: how the kinetic energy divides between vessel and structure and the damage that results.' },
      ],
      notes: [`${mono ? `Monopile clamped ${fmt(inp.fixD * D)} m below the mudline, impact ${fmt(h)} m above the fixity point` : `Member span ${inp.span} m, ${inp.ends} ends`}; D/t = ${fmt(D / t)}.`,
        `Accidental impact: ${cat.category} design, the structure absorbs ${fmt(acc.structureFraction * 100)} % of ${fmt(eAcc / 1e6)} MJ.${pi ? ` Global bending stays ${pi.elastic ? 'elastic' : 'plastic (hinge rotation ' + fmt(pi.plasticRotation) + ' rad)'}.` : ''}`,
        `${gov.name} governs (utilisation ${fmt(gov.uc)}).`, acc.exceeded ? 'The energy cannot be dissipated within the deformation limits entered: one of the bodies fails.' : 'The energy is dissipated within the deformation limits of both bodies.'],
      advice: [
        ok ? adv('engineering', sev(gov.uc), `Both impacts are survivable: accidental dent ${fmt(dAcc * 1e3)} mm (${fmt(dAcc / D * 100)} % of D), contact force ${fmt(acc.force / 1e6)} MN, residual bending capacity ${fmt(res * 100)} %.`,
          `Document ${fmt(inp.mass)} t and ${inp.vAcc} m/s as the limits of the marine operations manual, add the impact moment ${fmt(Mimp / 1e6)} MN·m to the accidental load combinations, and inspect for dents after any reported contact.`,
          'DNV-ST-0126 (2021) §4.3.6; DNV-RP-C204 (2019) §3',
          `No strengthening is needed. Fendering on the boat landing absorbs operational energy without steel in the primary structure; a dent that needs a clamp repair costs several vessel days — ${vesselCo2Text(3)}`)
          : adv('engineering', 'fail', `${gov.name}: ${fmt(gov.demand)} against ${fmt(gov.capacity)} ${gov.unit}.`,
            gov.name.includes('dent') ? `Increase the wall in the collision zone to about ${fmt(tReq * 1e3)} mm, add a sacrificial fender or boat-landing structure, or restrict the vessel size and approach speed.` : gov.name.includes('Post') ? `The dented section cannot carry ${inp.Mpost} MN·m: thicken the collision zone to about ${fmt(tReq * 1e3)} mm or plan immediate repair and restricted operation after an impact.` : 'Restrict vessel displacement or speed in the operations manual, or add energy-absorbing fendering so that the structure stays within its deformation limit.',
            'DNV-RP-C204 (2019) §3.3 and §3.10; DNV-ST-0126 (2021) §4.3.6',
            `Thickening an 8 m collision zone to ${fmt(tReq * 1e3)} mm adds ${fmt(dM)} t — ${steelCo2Text(dM)} Operational limits on vessel speed achieve the same safety with no material.`),
        adv('economic', 'info', `The accidental case puts ${fmt(acc.structureEnergy / 1e6)} MJ into the foundation and leaves ${fmt(res * 100)} % of its bending capacity; the operational case leaves a ${fmt(dOp * 1e3)} mm dent.`,
          'Weigh a thicker collision zone or fender (one-off cost per foundation) against the consequence of a damaged monopile: weeks of lost production, a repair campaign and possible derating. Limiting the authorised vessel size is free but constrains the O&M strategy.',
          'Risk-based design for accidental loads (ISO 19900 exposure levels; DNV-RP-C204 §2); O&M logistics costing'),
      ],
    };
  },
  verify: [
    { name: 'Collision energy with added mass: E = ½·m·Cm·v²', source: 'DNV-RP-C204 (2019) §3.2: 0.5 × 3.5e6 kg × 1.4 × (2 m/s)² = 9.8 MJ (accidental) and 0.5 × 3.5e6 × 1.1 × 0.5² = 0.481 MJ (operational)', inputs: { mass: 3500, vAcc: 2, cmSide: 1.4, vOp: 0.5, cmBow: 1.1 }, expect: [{ kpi: 'eAcc', value: 9.8, tol: 1e-6 }, { kpi: 'eOp', value: 0.48125, tol: 1e-6 }] },
    { name: 'Residual plastic moment of a dented tube: M/Mp = cos(θ/2) − ½ sin θ, θ = arccos(1 − 2·wd/D)', source: 'DNV-RP-C204 (2019) §3.10.2 hand calculation for wd/D = 0.1: θ = arccos(0.8) = 0.6435 rad → cos(0.3218) − 0.5 × 0.6 = 0.9487 − 0.3 = 0.6487', inputs: { mass: 20, vAcc: 0.2, vOp: 0.05, dent0: 600, D: 6 }, expect: [{ kpi: 'resM', value: 64.87, tol: 0.001 }] },
    { name: 'Collapse force of a cantilever pile: F = Mp / h with Mp = fy·(D³ − Di³)/6', source: 'Plastic hinge at the fixity point: D = 6 m, t = 70 mm, fy = 355 MPa → Zp = (216 − 201.23)/6 = 2.4617 m³, Mp = 873.9 MN·m; h = 3 + 25 + 3.5 × 6 = 49 m → F = 17.83 MN', inputs: { D: 6, t: 70, depth: 25, zImp: 3, fixD: 3.5 }, expect: [{ kpi: 'Fcollapse', value: 17.83, tol: 0.003 }] },
  ],
};
