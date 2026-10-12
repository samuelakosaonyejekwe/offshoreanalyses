// Structural collapse and accidental damage: ship-collision energy with a compliant floating unit, sharing of the
// strain energy between the striking bow and the struck side structure from their force–deformation curves,
// indentation against rupture and against the inner boundary, and the residual ultimate hull-girder (or pontoon
// box-girder) strength with the damaged shell removed.
import { collisionEnergy, collisionStrainEnergy, plateImpactCapacity, sharedEnergy, classifyCollisionDesign } from '../../engine/structural/impact.js';
import { hullGirderSection, hullGirderUltimate } from '../../engine/structural/plates.js';
import { LIMIT_CFD, need, run, round, rnd, fmt, adv, check, sev, worst, linspace, steelCo2Text } from './_common.js';

const trap = (d, F) => { let e = 0; for (let i = 1; i < d.length; i++) e += 0.5 * (F[i] + F[i - 1]) * (d[i] - d[i - 1]); return e; };

export default {
  id: 's11-accidental-damage',
  suite: 11,
  title: 'Collision damage and residual strength',
  summary: 'How much energy a ship collision puts into the floating unit, how deep the side structure is dented and whether it ruptures, and how much hull-girder bending strength is left with the damaged shell removed.',
  covers: [17],
  fidelity: 'analytical',
  standards: ['DNV-RP-C204 (2019) ch. 3 (ship collision: energy, force–deformation, strain-energy sharing)', 'NORSOK N-004 Rev. 3 (2013) Annex A.3', 'NORSOK N-003 (2017) §9.3 (collision energies)', 'DNV-OS-A101 (2021) Ch.2 Sec.1 (accidental loads)', 'DNV-OS-C101 (2021) Ch.2 Sec.6 (accidental limit state)', 'IACS CSR-H (2022) Pt 1 Ch 5 Sec 3 (residual hull-girder strength)'],
  models: ['Collision-impact energy with added mass and conservation of momentum (compliant installation)', 'Strain-energy sharing from force–deformation curves in series', 'Rigid-plastic plate with membrane action under a blunt punch (yield-line bending then membrane)', 'Rupture limit from a critical membrane strain', 'Smith-type progressive collapse of the hull girder (elastoplastic elements, buckling-limited compression)', 'Residual strength index of the damaged section'],
  limits: `Central, perfectly plastic impact; the struck panel is an equivalent circular plate with stiffeners smeared into the thickness and a strain-based rupture limit, and the bow curve is a user table. Residual strength removes the damaged shell strips and keeps a horizontal neutral axis (no rotation for asymmetric damage, no shear or torsion interaction). Nonlinear explicit finite-element collision analysis, fracture propagation and dynamic collapse are outside the module. ${LIMIT_CFD}`,
  searchTerms: 'ship collision offshore strain energy dissipation residual hull girder strength',
  inputs: [
    { key: 'shipMassT', label: 'Displacement of the striking vessel', unit: 't', type: 'number', default: 7500, min: 100, max: 300000, group: 'Collision event', role: 'ic', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Supply vessels: 5 000 t was the historic design value; NORSOK N-003 (2017) raises it to about 10 000 t. Shuttle tankers and passing ships are far larger.' },
    { key: 'vShip', label: 'Impact speed', unit: 'm/s', type: 'number', default: 2.0, min: 0, max: 12, group: 'Collision event', role: 'ic', uncertain: { dist: 'lognormal', cov: 0.25 }, help: 'Speed of the vessel at contact: 0.5 m/s for a drifting vessel in a small sea state, 2 m/s for the traditional accidental case, 3 m/s and more for powered head-on impact.' },
    { key: 'CmShip', label: 'Added-mass factor of the vessel', unit: '-', type: 'number', default: 1.1, min: 1, max: 2, group: 'Collision event', role: 'bc', help: '1.1 for bow or stern impact (surge), 1.4 for broadside impact (sway).' },
    { key: 'unitMassT', label: 'Displacement of the floating unit', unit: 't', type: 'number', default: 45000, min: 500, max: 2e6, group: 'Collision event', role: 'ic', bind: 'vessel.disp', help: 'Mass of the struck unit. A compliant unit recoils and dissipates less strain energy than a fixed one.' },
    { key: 'CmUnit', label: 'Added-mass factor of the unit', unit: '-', type: 'number', default: 1.4, min: 1, max: 3, group: 'Collision event', role: 'bc', help: 'Added mass of the unit in the direction of impact (1.2–2 for column-stabilised units and ships struck broadside).' },
    { key: 'fixed', label: 'Treat the unit as fixed', type: 'bool', default: false, group: 'Collision event', role: 'bc', help: 'Conservative: all the kinetic energy of the vessel is dissipated as strain energy.' },
    { key: 'shipCurve', label: 'Force–indentation curve of the striking bow', type: 'table', group: 'Striking vessel', role: 'bc',
      columns: [{ key: 'd', label: 'Indentation of the vessel', unit: 'm', type: 'number' }, { key: 'F', label: 'Force', unit: 'MN', type: 'number' }],
      default: [{ d: 0, F: 0 }, { d: 0.25, F: 8 }, { d: 0.5, F: 14 }, { d: 1, F: 20 }, { d: 1.5, F: 24 }, { d: 2, F: 27 }, { d: 3, F: 30 }],
      help: 'Crushing force of the vessel against a rigid wall, starting at zero and non-decreasing. The default is an indicative supply-vessel bow shape: replace it with the DNV-RP-C204 recommended curve for the actual vessel, impact location and contact width.' },
    { key: 'steel', label: 'Hull steel', type: 'select', source: 'materials:steel', default: 'ah36', group: 'Struck structure', role: 'material', help: 'Grade of the struck shell; sets the yield stress.' },
    { key: 'tPlate', label: 'Shell plate thickness', unit: 'mm', type: 'number', default: 20, min: 6, max: 80, group: 'Struck structure', role: 'geometry', design: true, help: 'Outer shell plating at the impact location.' },
    { key: 'sStiff', label: 'Stiffener spacing', unit: 'mm', type: 'number', default: 700, min: 300, max: 1500, group: 'Struck structure', role: 'geometry', help: 'Spacing of the shell stiffeners, used to smear their area into the plate.' },
    { key: 'aStiff', label: 'Stiffener cross-section area', unit: 'cm²', type: 'number', default: 60, min: 0, max: 400, group: 'Struck structure', role: 'geometry', design: true, help: 'Area of one stiffener (0 = unstiffened plate). Smeared stiffeners add membrane capacity in their own direction only: the equivalent thickness is an upper estimate.' },
    { key: 'lSpan', label: 'Panel span between supports', unit: 'm', type: 'number', default: 3, min: 0.6, max: 12, group: 'Struck structure', role: 'geometry', design: true, help: 'Shorter distance between the frames, decks or stringers that restrain the panel; the equivalent circular plate has half this as its radius.' },
    { key: 'aPunch', label: 'Contact radius of the striking object', unit: 'm', type: 'number', default: 0.5, min: 0.05, max: 5, group: 'Struck structure', role: 'bc', help: 'Radius of the contact patch (bulb tip, stem, stern roller). Must be smaller than half the panel span.' },
    { key: 'edge', label: 'Panel edge restraint', type: 'select', options: [{ value: 'clamped', label: 'Clamped (continuous over frames)' }, { value: 'simple', label: 'Simply supported' }], default: 'clamped', group: 'Struck structure', role: 'bc', help: 'Rotational restraint at the panel boundary; in-plane restraint (membrane action) is assumed in both cases.' },
    { key: 'epsCr', label: 'Critical membrane strain at rupture', unit: '-', type: 'number', default: 0.12, min: 0.02, max: 0.3, group: 'Struck structure', role: 'criteria', uncertain: { dist: 'lognormal', cov: 0.25 }, help: 'Average strain over the deformed span at which the shell tears. 0.10–0.15 is customary for mild and high-tensile hull steel in simplified methods; calibrate against DNV-RP-C208 for fine meshes.' },
    { key: 'dInner', label: 'Distance to the inner boundary', unit: 'm', type: 'number', default: 1.5, min: 0.1, max: 10, group: 'Struck structure', role: 'criteria', help: 'Depth of the double side, or clearance to the first tank, riser or equipment behind the shell that must not be reached.' },
    { key: 'B', label: 'Girder breadth', unit: 'm', type: 'number', default: 46, min: 4, max: 100, group: 'Residual strength section', role: 'geometry', help: 'Breadth of the hull (or of a pontoon box girder).' },
    { key: 'Dg', label: 'Girder depth', unit: 'm', type: 'number', default: 27, min: 2, max: 60, group: 'Residual strength section', role: 'geometry', help: 'Depth from bottom to strength deck.' },
    { key: 'tDeck', label: 'Deck thickness (equivalent, with stiffeners)', unit: 'mm', type: 'number', default: 30, min: 5, max: 120, group: 'Residual strength section', role: 'geometry', design: true, help: 'Plate plus smeared longitudinals.' },
    { key: 'tBottom', label: 'Bottom thickness (equivalent)', unit: 'mm', type: 'number', default: 32, min: 5, max: 120, group: 'Residual strength section', role: 'geometry', help: 'Plate plus smeared longitudinals (include the inner bottom).' },
    { key: 'tSide', label: 'Side thickness (equivalent, each side)', unit: 'mm', type: 'number', default: 26, min: 5, max: 120, group: 'Residual strength section', role: 'geometry', help: 'Side shell plus smeared longitudinals (include an inner skin).' },
    { key: 'sigUfrac', label: 'Compressive strength of stiffened panels / yield', unit: '-', type: 'number', default: 0.8, min: 0.3, max: 1, group: 'Residual strength section', role: 'material', help: 'Ultimate compressive stress of the stiffened deck, bottom and side panels as a fraction of yield (from the local strength module).' },
    { key: 'zDamLow', label: 'Damage extent, lower edge above base', unit: 'm', type: 'number', default: 12, min: 0, max: 60, group: 'Residual strength section', role: 'ic', help: 'Lower edge of the side shell assumed ineffective after the collision.' },
    { key: 'zDamHigh', label: 'Damage extent, upper edge above base', unit: 'm', type: 'number', default: 27, min: 0, max: 60, group: 'Residual strength section', role: 'ic', help: 'Upper edge of the ineffective side shell; reaching the deck also removes a deck strip.' },
    { key: 'bDeckDam', label: 'Deck width lost with the damage', unit: 'm', type: 'number', default: 3, min: 0, max: 50, group: 'Residual strength section', role: 'ic', help: 'Width of deck (stringer plate and adjacent strakes) taken as ineffective when the damage reaches the deck.' },
    { key: 'Mdam', label: 'Bending moment in the damaged condition', unit: 'MN·m', type: 'number', default: 6500, min: 0, max: 100000, group: 'Residual strength section', role: 'bc', help: 'Still-water moment of the damaged condition plus the reduced wave moment (CSR-H uses 67 % of the rule wave moment), from the global strength module.' },
    { key: 'gammaR', label: 'Partial factor on residual capacity', unit: '-', type: 'number', default: 1.0, min: 1, max: 1.5, group: 'Acceptance criteria', role: 'criteria', help: 'Accidental limit state: 1.0 (CSR-H residual strength uses 1.0).' },
  ],
  compute(inp, ctx) {
    const mat = ctx.material(inp.steel); need(mat && mat.fy > 0, 'Select a hull steel from the materials library.');
    const fy = mat.fy; const E = mat.E || 206e9;
    // --- energy ---------------------------------------------------------------------------------------------
    const ke = collisionEnergy({ mass: inp.shipMassT * 1e3, velocity: inp.vShip, addedMassFactor: inp.CmShip });
    const mi = inp.fixed ? Infinity : inp.unitMassT * 1e3 * inp.CmUnit;
    const se = collisionStrainEnergy({ shipMass: ke.effectiveMass, shipVelocity: inp.vShip, installationMass: mi });
    const Es = se.strainEnergy;
    // --- curves ----------------------------------------------------------------------------------------------
    const pts = (inp.shipCurve || []).filter((r) => Number.isFinite(r.d) && Number.isFinite(r.F)).sort((a, b) => a.d - b.d);
    need(pts.length >= 2 && pts[0].d === 0, 'The bow curve needs at least two points and must start at zero indentation.');
    for (let i = 1; i < pts.length; i++) need(pts[i].d > pts[i - 1].d && pts[i].F >= pts[i - 1].F, 'The bow curve must have increasing indentation and non-decreasing force: enter the envelope of a softening curve.');
    const shipCurve = { d: pts.map((p) => p.d), F: pts.map((p) => p.F * 1e6) };
    const Rp = inp.lSpan / 2;
    need(inp.aPunch < Rp, `The contact radius (${fmt(inp.aPunch)} m) must be smaller than half the panel span (${fmt(Rp)} m): the object would bear directly on the supports.`);
    const tEq = inp.tPlate / 1e3 + (inp.aStiff * 1e-4) / (inp.sStiff / 1e3);
    // average radial strain of the cone between punch and support: ε ≈ ½·(δ/(Rp − a))² → δ at rupture
    const dRupt = (Rp - inp.aPunch) * Math.sqrt(2 * inp.epsCr);
    const plate = (d) => run(() => plateImpactCapacity({ t: tEq, fy, radius: Rp, punchRadius: inp.aPunch, edge: inp.edge, membrane: true, deflection: d }), 'Plate impact resistance');
    const ds = linspace(0, dRupt, 41); const st = ds.map(plate);
    const structureCurve = { d: ds, F: st.map((s) => s.resistance) };
    const sh = run(() => sharedEnergy({ shipCurve, structureCurve, energy: Es }), 'Strain-energy sharing');
    const cat = sh.structureEnergy + sh.shipEnergy > 0 ? classifyCollisionDesign({ structureEnergy: sh.structureEnergy, shipEnergy: sh.shipEnergy }) : { category: 'strength (no impact energy)' };
    const capStruct = trap(structureCurve.d, structureCurve.F);
    const ruptured = sh.exceeded || sh.structureDeformation >= dRupt * (1 - 1e-9);
    // --- residual strength ------------------------------------------------------------------------------------
    need(inp.zDamHigh >= inp.zDamLow, 'The upper edge of the damage must not be below its lower edge.');
    const D = inp.Dg; const B = inp.B; const su = inp.sigUfrac * fy; const nv = 24; const dz = D / nv;
    need(inp.bDeckDam < B, 'The lost deck width must be smaller than the girder breadth.');
    const section = (damaged) => {
      const zL = Math.min(inp.zDamLow, D); const zH = Math.min(inp.zDamHigh, D);
      const toDeck = damaged && zH >= D - 1e-9 && zH > zL;
      const el = [{ area: (B - (toDeck ? inp.bDeckDam : 0)) * inp.tDeck / 1e3, z: D, fy, sigmaU: su }, { area: B * inp.tBottom / 1e3, z: 0, fy, sigmaU: su }];
      for (let i = 0; i < nv; i++) {
        const z0 = i * dz; const z1 = z0 + dz; const z = 0.5 * (z0 + z1);
        el.push({ area: dz * inp.tSide / 1e3, z, i0: inp.tSide / 1e3 * dz ** 3 / 12, fy, sigmaU: su });
        // struck side: only the part of each strip outside the damage extent remains effective
        const lost = damaged ? Math.max(0, Math.min(z1, zH) - Math.max(z0, zL)) : 0;
        if (dz - lost > 1e-9) el.push({ area: (dz - lost) * inp.tSide / 1e3, z, i0: inp.tSide / 1e3 * (dz - lost) ** 3 / 12, fy, sigmaU: su });
      }
      return el;
    };
    const el0 = section(false); const el1 = section(true);
    const s0 = hullGirderSection({ elements: el0 }); const s1 = hullGirderSection({ elements: el1 });
    const u0 = run(() => hullGirderUltimate({ elements: el0, E, nSteps: 120 }), 'Intact ultimate strength');
    const u1 = run(() => hullGirderUltimate({ elements: el1, E, nSteps: 120 }), 'Damaged ultimate strength');
    const Mu0 = Math.min(u0.sagging.Mu, u0.hogging.Mu); const Mu1 = Math.min(u1.sagging.Mu, u1.hogging.Mu);
    const rsi = Mu1 / Mu0; const areaLost = 1 - s1.A / s0.A;
    const checks = [
      check('Strain energy against the absorption capacity up to shell rupture', Es / 1e6, sh.capacity / 1e6, 'MJ', 'DNV-RP-C204 (2019) §3.4 (dissipation of strain energy); NORSOK N-004 (2013) A.3.5', 'Guards against tearing of the outer shell: if the two bodies cannot absorb the collision energy before the panel reaches its rupture strain, the shell opens and the compartment floods.'),
      check('Shell indentation against the rupture deflection', sh.structureDeformation, dRupt, 'm', 'DNV-RP-C204 (2019) §3.10 (ductility limits, critical strain); DNV-RP-C208 (2019)', 'Guards the watertight integrity of the struck panel: beyond this deflection the membrane strain exceeds what the plating and its welds can carry.'),
      check('Shell indentation against the inner boundary', sh.structureDeformation, inp.dInner, 'm', 'DNV-OS-A101 (2021) Ch.2 Sec.1 (accidental loads not to impair the main safety functions); NORSOK N-004 (2013) A.3', 'Guards what lies behind the shell: an indentation deeper than the double side reaches cargo tanks, risers or the inner watertight boundary.'),
      check('Bending moment against the residual hull-girder capacity', inp.Mdam, Mu1 / inp.gammaR / 1e6, 'MN·m', 'IACS CSR-H (2022) Pt 1 Ch 5 Sec 3 (residual strength); DNV-OS-C101 (2021) Ch.2 Sec.6 (ALS, post-damage condition)', 'Guards against progressive collapse after the accident: the damaged girder must still carry the still-water and reduced wave bending until the unit is repaired or made safe.'),
    ];
    const gov = worst(checks);
    // plate thickness that just avoids rupture: membrane energy scales linearly with thickness
    const eNeedStruct = Math.max(Es - sh.shipEnergy, 0); const tNeed = sh.exceeded ? inp.tPlate * 1.25 : tEq * 1e3 * Math.min(eNeedStruct / Math.max(capStruct, 1), 5);
    const panelSteel = (dt) => Math.PI * Rp * Rp * 4 * (dt / 1e3) * 7.85; // four panels around the impact zone
    return {
      kpis: [
        { key: 'Ekin', label: 'Kinetic energy of the vessel', value: round(ke.energy / 1e6, 5), unit: 'MJ', explain: '½·(m + added mass)·v² of the striking vessel just before contact.' },
        { key: 'Estrain', label: 'Strain energy to dissipate', value: round(Es / 1e6, 5), unit: 'MJ', explain: 'Part of the kinetic energy that must be absorbed by deformation of the two bodies; the rest remains as motion of the vessel and the recoiling unit.' },
        { key: 'vCommon', label: 'Common velocity after impact', value: round(se.commonVelocity, 5), unit: 'm/s', explain: 'Velocity of vessel and unit moving together at the end of a perfectly plastic impact (zero for a fixed unit).' },
        { key: 'PcPlate', label: 'Plastic collapse load of the panel', value: round(st[0].Pc / 1e3, 5), unit: 'kN', explain: 'Force at which the plate forms a yield-line mechanism under the contact patch; above it the panel deflects plastically and membrane action takes over.' },
        { key: 'dRupture', label: 'Deflection at shell rupture', value: round(dRupt, 5), unit: 'm', explain: 'Indentation at which the average membrane strain of the panel reaches the critical strain.' },
        { key: 'capStruct', label: 'Energy capacity of the panel up to rupture', value: round(capStruct / 1e6, 5), unit: 'MJ', explain: 'Area under the panel force–deflection curve up to rupture: what the struck structure can absorb by itself.' },
        { key: 'Fcontact', label: 'Contact force at the end of impact', value: round(sh.force / 1e6, 5), unit: 'MN', explain: 'Common force in bow and panel when the strain energy has been absorbed; the load for the supporting frames and the global check.' },
        { key: 'Estruct', label: 'Energy absorbed by the unit', value: round(sh.structureEnergy / 1e6, 5), unit: 'MJ', explain: 'Share of the strain energy dissipated in the struck panel.' },
        { key: 'Eship', label: 'Energy absorbed by the vessel', value: round(sh.shipEnergy / 1e6, 5), unit: 'MJ', explain: 'Share dissipated by crushing of the striking bow.' },
        { key: 'structShare', label: 'Share absorbed by the unit', value: round(sh.structureFraction, 5), unit: '-', explain: 'Below about 0.1 the design is "strength" (the ship does the work), above 0.9 "ductile" (the unit does), otherwise "shared".' },
        { key: 'indent', label: 'Indentation of the shell', value: round(sh.structureDeformation, 5), unit: 'm', limit: round(Math.min(dRupt, inp.dInner), 5), sense: 'max', explain: 'Permanent deflection of the struck panel at the end of the impact, to compare with the rupture deflection and the distance to the inner boundary.' },
        { key: 'ruptured', label: 'Shell ruptured (1 = yes)', value: ruptured ? 1 : 0, unit: '-', limit: 0, sense: 'max', explain: 'Whether the panel reaches its rupture deflection before the collision energy is absorbed. A ruptured shell means the compartment behind it floods.' },
        { key: 'indentShip', label: 'Crushing of the vessel bow', value: round(sh.shipDeformation, 5), unit: 'm', explain: 'Indentation of the striking vessel at the end of the impact.' },
        { key: 'MuIntact', label: 'Ultimate bending capacity, intact', value: round(Mu0 / 1e6, 5), unit: 'MN·m', explain: 'Lower of the sagging and hogging ultimate moments of the intact section (progressive collapse).' },
        { key: 'MuDamaged', label: 'Ultimate bending capacity, damaged', value: round(Mu1 / 1e6, 5), unit: 'MN·m', explain: 'Same with the damaged side shell and deck strip removed.' },
        { key: 'rsi', label: 'Residual strength index', value: round(rsi, 5), unit: '-', explain: 'Damaged over intact ultimate moment. Typical collision damage leaves 0.85–0.95; much lower values mean the damage has reached the deck or bottom flange.' },
        { key: 'areaLost', label: 'Section area lost', value: round(areaLost, 5), unit: '-', explain: 'Fraction of the longitudinally effective steel area removed by the assumed damage extent.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Force–deformation curves and the end-of-impact point', x: { label: 'Deformation', unit: 'm' }, y: { label: 'Force', unit: 'MN' },
          series: [{ name: 'Striking bow', x: rnd(shipCurve.d), y: rnd(shipCurve.F.map((f) => f / 1e6)) }, { name: 'Struck panel (to rupture)', x: rnd(ds), y: rnd(structureCurve.F.map((f) => f / 1e6)) }, { name: 'Contact force', x: [0, round(Math.max(shipCurve.d[shipCurve.d.length - 1], dRupt))], y: [round(sh.force / 1e6), round(sh.force / 1e6)], dash: true }],
          explain: 'Both bodies carry the same contact force (dashed). Where that line cuts each curve is the deformation of that body, and the areas under the curves up to there add up to the strain energy. A panel curve that ends below the line has ruptured.' },
        { type: 'bar', title: 'Energy balance', categories: ['Kinetic energy', 'Strain energy', 'Absorbed by vessel', 'Absorbed by unit', 'Panel capacity to rupture'], y: { label: 'Energy', unit: 'MJ' },
          series: [{ name: 'Energy', y: rnd([ke.energy / 1e6, Es / 1e6, sh.shipEnergy / 1e6, sh.structureEnergy / 1e6, capStruct / 1e6]) }],
          explain: 'Where the collision energy goes. The difference between kinetic and strain energy stays as motion; the unit share must stay below the panel capacity for the shell to remain intact.' },
        { type: 'line', title: 'Moment–curvature of the girder, intact and damaged (sagging)', x: { label: 'Curvature', unit: '1/km' }, y: { label: 'Bending moment', unit: 'MN·m' },
          series: [{ name: 'Intact', x: rnd(u0.sagging.curvature.map((c) => Math.abs(c) * 1e3)), y: rnd(u0.sagging.moment.map((m) => Math.abs(m) / 1e6)) }, { name: 'Damaged', x: rnd(u1.sagging.curvature.map((c) => Math.abs(c) * 1e3)), y: rnd(u1.sagging.moment.map((m) => Math.abs(m) / 1e6)) }, { name: 'Moment in the damaged condition', x: [0, round(Math.abs(u0.sagging.curvature[u0.sagging.curvature.length - 1]) * 1e3)], y: [inp.Mdam, inp.Mdam], dash: true }],
          explain: 'Progressive-collapse curves: the peak of each is the ultimate moment. The damaged curve must still peak above the dashed demand line; the drop between the two peaks is the strength lost to the collision.' },
      ],
      tables: [
        { title: 'Collision summary', columns: [{ key: 'item', label: 'Quantity' }, { key: 'value', label: 'Value' }, { key: 'unit', label: 'Unit' }],
          rows: [
            { item: 'Vessel mass incl. added mass', value: round(ke.effectiveMass / 1e3, 5), unit: 't' }, { item: 'Unit mass incl. added mass', value: inp.fixed ? 'fixed' : round(mi / 1e3, 5), unit: 't' },
            { item: 'Equivalent shell thickness', value: round(tEq * 1e3, 4), unit: 'mm' }, { item: 'Equivalent plate radius', value: round(Rp, 4), unit: 'm' },
            { item: 'Design category', value: cat.category, unit: '' }, { item: 'Shell ruptured', value: ruptured ? 'yes' : 'no', unit: '' },
            { item: 'Neutral axis, intact / damaged', value: `${fmt(s0.zNA, 4)} / ${fmt(s1.zNA, 4)}`, unit: 'm' }, { item: 'Sagging Mu, intact / damaged', value: `${fmt(u0.sagging.Mu / 1e6, 4)} / ${fmt(u1.sagging.Mu / 1e6, 4)}`, unit: 'MN·m' },
            { item: 'Hogging Mu, intact / damaged', value: `${fmt(u0.hogging.Mu / 1e6, 4)} / ${fmt(u1.hogging.Mu / 1e6, 4)}`, unit: 'MN·m' },
          ],
          explain: 'Intermediate values of the collision and residual-strength calculation, for the accidental limit state report.' },
      ],
      notes: [
        `Kinetic energy ${fmt(ke.energy / 1e6)} MJ; ${inp.fixed ? 'the unit is taken as fixed, so all of it' : `${fmt(Es / 1e6)} MJ (${fmt(100 * Es / Math.max(ke.energy, 1e-9))} %)`} is strain energy. Design category: ${cat.category} (unit absorbs ${fmt(100 * sh.structureFraction)} %).`,
        ruptured ? `The shell ruptures: the panel reaches its rupture deflection of ${fmt(dRupt)} m with ${fmt(sh.capacity / 1e6)} MJ absorbed against ${fmt(Es / 1e6)} MJ required.` : `The shell is dented ${fmt(sh.structureDeformation)} m (rupture at ${fmt(dRupt)} m, inner boundary at ${fmt(inp.dInner)} m) at a contact force of ${fmt(sh.force / 1e6)} MN.`,
        `Damage from ${fmt(inp.zDamLow)} to ${fmt(inp.zDamHigh)} m removes ${fmt(100 * areaLost)} % of the section area and leaves ${fmt(100 * rsi)} % of the ultimate moment. Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at utilisation ${fmt(gov.uc)}: indentation ${fmt(sh.structureDeformation)} m, residual strength index ${fmt(rsi)}.`,
          ruptured ? `Raise the panel energy capacity from ${fmt(capStruct / 1e6)} to about ${fmt(eNeedStruct / 1e6)} MJ: an equivalent shell thickness near ${fmt(tNeed)} mm, a longer span between hard supports, or fendering that spreads the contact; otherwise design the compartment behind as flooded (damage stability and progressive flooding modules).`
            : sh.structureDeformation > inp.dInner ? `Stiffen the shell so that the indentation stays below ${fmt(inp.dInner)} m, or widen the double side; the contact force of ${fmt(sh.force / 1e6)} MN should then be carried by the frames.`
              : checks[3].uc > 1 ? `Restore residual capacity: ${fmt(inp.Mdam * inp.gammaR - Mu1 / 1e6)} MN·m is missing. Thicken the deck stringer and sheer strake (they carry most of the lost moment), or restrict the damaged-condition loading to lower the still-water moment.`
                : `The collision is absorbed with ${fmt(100 * (1 - sh.structureDeformation / dRupt))} % deflection margin to rupture and ${fmt(100 * (1 - checks[3].uc))} % reserve on residual strength. Check the supporting frames for the ${fmt(sh.force / 1e6)} MN contact force, and confirm the bow curve for the actual vessel class.`,
          'DNV-RP-C204 (2019) ch. 3; NORSOK N-004 (2013) Annex A.3; IACS CSR-H (2022) Pt 1 Ch 5 Sec 3; DNV-OS-C101 (2021) Ch.2 Sec.6.'),
        adv('economic', 'info', `Raising the shell by 5 mm over the four panels around the collision zone is about ${fmt(panelSteel(5))} t of steel per zone; a ruptured shell instead means flooding, production shut-in and an offshore or dry-dock repair.`,
          'Decide the collision design energy with the operator early: moving from the 2 m/s supply-vessel case to the 10 000 t / 3 m/s case roughly triples the energy and is far cheaper to accommodate in the frame spacing and fendering before the hull is ordered.',
          'Accidental limit state design to a 10⁻⁴ annual event (NORSOK N-003 (2017) §9; DNV-OS-A101 (2021)); risk-based inspection after impact per DNV-RP-G101 (2021).',
          `${steelCo2Text(panelSteel(5))} per collision zone; repairing a torn shell offshore typically takes weeks of support-vessel time, which outweighs the embodied carbon of the extra plate many times over.`),
      ],
    };
  },
  verify: [
    { // Vessel 7500 t × 1.1 = 8250 t at 2 m/s: E = ½·8.25e6·2² = 16.50 MJ. Unit 45 000 t × 1.4 = 63 000 t:
      // Es = 16.50/(1 + 8250/63 000) = 16.50/1.13095 = 14.589 MJ; common velocity = 8250·2/71 250 = 0.2316 m/s.
      // Plate 20 mm (no stiffeners), fy = 355 MPa: mp = 355e6·0.02²/4 = 35.5 kN·m/m; clamped circular plate R = 1.5 m, punch 0.5 m:
      // Pc = 4π·mp·R/(R − a) = 4π·35 500·1.5 = 669.2 kN. Rupture deflection = (1.5 − 0.5)·√(2·0.12) = 0.4899 m.
      name: 'Collision strain energy of a compliant unit and plastic collapse load of the struck plate',
      source: 'DNV-RP-C204 (2019) Eq. 3.1–3.2; yield-line solution for a clamped circular plate (Jones, Structural Impact, 2012, ch. 2) — hand calculation',
      inputs: { aStiff: 0 },
      expect: [{ kpi: 'Ekin', value: 16.5, tol: 1e-6 }, { kpi: 'Estrain', value: 14.589, tol: 1e-4 }, { kpi: 'vCommon', value: 0.2316, tol: 1e-3 }, { kpi: 'PcPlate', value: 669.2, tol: 1e-3 }, { kpi: 'dRupture', value: 0.4899, tol: 1e-3 }] },
    { name: 'Fixed unit dissipates the whole kinetic energy; no damage extent leaves the full strength', source: 'Limiting cases: installation mass → ∞ (DNV-RP-C204 Eq. 3.3) and zero damage',
      inputs: { fixed: true, zDamLow: 10, zDamHigh: 10 },
      expect: [{ kpi: 'Estrain', value: 16.5, tol: 1e-6 }, { kpi: 'rsi', value: 1, tol: 1e-9 }, { kpi: 'areaLost', value: 0, abs: 1e-12 }] },
    { // Fully plastic check of the progressive-collapse routine through the module: with σU = fy and symmetric flanges the
      // ultimate moment of a thin box B × D is Mp = fy·(B·t_f·D + 2·t_s·D²/4). B = 20, D = 10 m, t_f = t_s = 20 mm, fy = 355 MPa:
      // Mp = 355e6·(20·0.02·10 + 2·0.02·100/4) = 355e6·(4 + 1) = 1775 MN·m.
      name: 'Plastic moment of a thin-walled box girder',
      source: 'Plastic section theory, hand calculation',
      inputs: { B: 20, Dg: 10, tDeck: 20, tBottom: 20, tSide: 20, sigUfrac: 1, zDamLow: 5, zDamHigh: 5, Mdam: 500 },
      expect: [{ kpi: 'MuIntact', value: 1775, tol: 0.02 }] },
  ],
};
