// Dropped-object impact on a deck: impact energy of an object falling in air, plastic energy absorption of the
// deck plate panel (bending then membrane) up to a strain-based deflection limit, the plastic mechanism of the
// supporting beam, and the allowable drop height against mass for lifting restrictions.
import { dropInAir, plateImpactCapacity } from '../../engine/structural/impact.js';
import { iSection } from '../../engine/structural/sections.js';
import { steelOf, need, run, fmt, round, linspace, adv, check, steelCo2Text, CONFIRM, NOT_FE, G } from './_common.js';

export default {
  id: 's02-dropped-object', suite: 2, title: 'Dropped-object impact on the deck',
  summary: 'Impact energy of a dropped load against the plastic energy absorption of the deck plate and its supporting beam, with plate deflection, penetration margin and the allowable drop height for each lifted mass.',
  covers: [10], fidelity: 'analytical',
  standards: ['DNV-RP-C204 (2019) §5 (dropped objects)', 'NORSOK N-004 (2022) Annex A.4 (dropped objects)', 'ISO 19901-3:2014 §8.5 (accidental actions, impact)', 'DNV-ST-N001 (2021) §16 (lifting over live plant)'],
  models: ['Newton’s second law (free fall)', 'Rigid-plastic yield-line mechanism of a plate', 'Plastic membrane action', 'Three-hinge plastic mechanism of a beam (plastic-hinge formulation)', 'Energy balance (plastic collapse)'],
  limits: `Rigid–perfectly-plastic energy balance for a blunt object striking the centre of a plate panel idealised as a circular plate: the object is rigid, strain-rate and inertia effects of the deck are neglected, and rupture is represented only by a limit on the average membrane strain. Sharp or rotating objects, impact on an edge or stiffener, and progressive tearing need explicit analysis. ${NOT_FE}`,
  searchTerms: 'dropped object impact deck plate energy absorption membrane plastic offshore',
  inputs: [
    { key: 'mass', label: 'Mass of the dropped object', unit: 't', type: 'number', default: 2, min: 0.01, max: 500, group: 'Dropped object', role: 'load', uncertain: { dist: 'normal', cov: 0.05 }, help: 'Heaviest item lifted over this deck area (container, drill collar, pipe bundle, tool basket).' },
    { key: 'height', label: 'Drop height above the deck', unit: 'm', type: 'number', default: 3, min: 0.05, max: 100, group: 'Dropped object', role: 'ic', design: true, help: 'Largest height of the load above the deck along its lifting path.' },
    { key: 'dContact', label: 'Contact diameter', unit: 'mm', type: 'number', default: 300, min: 10, max: 3000, group: 'Dropped object', role: 'geometry', help: 'Diameter of the part of the object that strikes the plate (corner fitting, pipe end, shackle). Smaller is more severe.' },
    { key: 'bPanel', label: 'Plate panel width (stiffener spacing)', unit: 'mm', type: 'number', default: 700, min: 100, max: 4000, group: 'Deck', role: 'geometry', design: true },
    { key: 'tPlate', label: 'Deck plate thickness', unit: 'mm', type: 'number', default: 12, min: 3, max: 80, group: 'Deck', role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.03 } },
    { key: 'edge', label: 'Plate edge rotation', type: 'select', default: 'clamped', group: 'Deck', role: 'bc', options: [{ value: 'clamped', label: 'Clamped (continuous plating)' }, { value: 'simple', label: 'Simply supported' }] },
    { key: 'membrane', label: 'Plate edges restrained in plane (membrane action)', type: 'bool', default: true, group: 'Deck', role: 'bc', help: 'True for continuous welded deck plating surrounded by other panels; false for an isolated hatch cover or panel at a free edge.' },
    { key: 'mat', label: 'Steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Deck', role: 'material' },
    { key: 'Lb', label: 'Span of the supporting beam', unit: 'm', type: 'number', default: 5, min: 0.5, max: 30, group: 'Supporting beam', role: 'geometry' },
    { key: 'hb', label: 'Beam depth', unit: 'mm', type: 'number', default: 300, min: 60, max: 2500, group: 'Supporting beam', role: 'geometry', design: true },
    { key: 'bb', label: 'Beam flange width', unit: 'mm', type: 'number', default: 150, min: 30, max: 1000, group: 'Supporting beam', role: 'geometry' },
    { key: 'twb', label: 'Beam web thickness', unit: 'mm', type: 'number', default: 8, min: 3, max: 60, group: 'Supporting beam', role: 'geometry' },
    { key: 'tfb', label: 'Beam flange thickness', unit: 'mm', type: 'number', default: 12, min: 4, max: 100, group: 'Supporting beam', role: 'geometry' },
    { key: 'beamEnds', label: 'Beam end fixity', type: 'select', default: 'clamped', group: 'Supporting beam', role: 'bc', options: [{ value: 'clamped', label: 'Continuous / fixed' }, { value: 'pinned', label: 'Pinned' }] },
    { key: 'epsCr', label: 'Critical average membrane strain of the plate', unit: '%', type: 'number', default: 5, min: 0.5, max: 20, group: 'Acceptance criteria', role: 'criteria', help: `Strain at which the plate or its welds are taken to tear. NORSOK N-004 gives 15–20 % for base material over a long gauge length; welded plating is commonly limited to about 5 %. ${CONFIRM}` },
    { key: 'rotLim', label: 'Plastic rotation limit of the beam', unit: 'deg', type: 'number', default: 4, min: 0.5, max: 15, group: 'Acceptance criteria', role: 'criteria', help: `Chord rotation the beam can reach before local buckling or connection failure. ${CONFIRM}` },
    { key: 'beamShare', label: 'Count the supporting beam as an energy absorber', type: 'bool', default: true, group: 'Acceptance criteria', role: 'criteria', help: 'Off when the plate alone must stop the object (items directly below the plate).' },
  ],
  compute(inp, ctx) {
    // SI: t → kg, mm → m, % → fraction.
    const st = steelOf(ctx, inp.mat), m = inp.mass * 1e3, t = inp.tPlate / 1e3, Rp = inp.bPanel / 2e3, a = inp.dContact / 2e3;
    need(a < 0.9 * Rp, 'The contact diameter must be smaller than 90 % of the panel width: a wider object loads the stiffeners directly.');
    const drop = dropInAir({ mass: m, height: inp.height });
    // Deflection at which the average radial membrane strain of the cone between punch and edge reaches the limit:
    // (1 + ε)² = 1 + (δ/(Rp − a))².
    const eps = inp.epsCr / 100, dLim = (Rp - a) * Math.sqrt((1 + eps) ** 2 - 1);
    const plate = (tt, d) => run(() => plateImpactCapacity({ t: tt, fy: st.fy, radius: Rp, punchRadius: a, edge: inp.edge, membrane: inp.membrane, deflection: d }), 'Plate capacity');
    const pl = plate(t, dLim);
    need(2 * inp.tfb < inp.hb && inp.twb < inp.bb, 'The supporting beam is not a valid I-shape: check its thicknesses.');
    const sec = run(() => iSection({ h: inp.hb / 1e3, b: inp.bb / 1e3, tw: inp.twb / 1e3, tf: inp.tfb / 1e3, welded: true }), 'Beam section');
    // Three-hinge mechanism under a midspan point load, bending only: R0 = 4 c Mp / L (c = 2 fixed, 1 pinned); E = R0 · w.
    const Mp = sec.Wply * st.fy, R0 = 4 * (inp.beamEnds === 'clamped' ? 2 : 1) * Mp / inp.Lb, wLim = (inp.Lb / 2) * Math.tan(inp.rotLim * Math.PI / 180), Ebeam = inp.beamShare ? R0 * wLim : 0;
    const E = drop.energy, Ecap = pl.energy + Ebeam;
    // Plate deflection at which the impact energy is absorbed (bisection on the plate energy), capped at the limit.
    let dAct = dLim;
    if (pl.energy > E) { let lo = 0, hi = dLim; for (let i = 0; i < 40; i++) { const mid = 0.5 * (lo + hi); if (plate(t, mid).energy > E) hi = mid; else lo = mid; } dAct = hi; }
    const Ebm = Math.max(Math.min(E - pl.energy, Ebeam), 0), wBeam = R0 > 0 && inp.beamShare ? Ebm / R0 : 0;
    // Plate thickness that alone absorbs the energy within the strain limit.
    let tReq = t;
    if (pl.energy < E) { let lo = t, hi = 2 * t; while (plate(hi, dLim).energy < E && hi < 0.3) hi *= 2; for (let i = 0; i < 30; i++) { const mid = 0.5 * (lo + hi); if (plate(mid, dLim).energy < E) lo = mid; else hi = mid; } tReq = hi; }
    const hAllowPlate = pl.energy / (m * G), hAllow = Ecap / (m * G);
    const checks = [
      check('Plate energy absorption (no tearing)', E / 1e3, pl.energy / 1e3, 'kJ', 'DNV-RP-C204 §5.3; NORSOK N-004 A.4.4 (plastic energy of plating up to the critical strain)', 'Impact energy against the plastic work the plate panel can absorb before its average membrane strain reaches the limit. Above 1.0 the plate is predicted to tear and the object reaches whatever lies below.'),
      check(inp.beamShare ? 'Plate and supporting beam together' : 'Plate alone (beam not counted)', E / 1e3, Ecap / 1e3, 'kJ', 'DNV-RP-C204 §5.2 (energy dissipation by the struck structure)', 'Impact energy against the combined plastic work of the plate and the three-hinge mechanism of the beam beneath it. Above 1.0 the deck cannot stop the object.'),
    ];
    const ok = checks[1].uc <= 1, plateOk = checks[0].uc <= 1;
    const ds = linspace(0, dLim, 41), masses = linspace(0.2, 2.5, 24).map((f) => f * inp.mass);
    const panelT = (inp.bPanel / 1e3) ** 2 * t * st.rho / 1e3;
    return {
      kpis: [
        { key: 'energy', label: 'Impact energy', unit: 'kJ', value: E / 1e3, explain: 'Kinetic energy at impact, m g h, for a fall in air without drag.' },
        { key: 'velocity', label: 'Impact velocity', unit: 'm/s', value: drop.velocity, explain: 'Velocity at impact, √(2 g h).' },
        { key: 'Pc', label: 'Plate bending collapse load', unit: 'kN', value: pl.Pc / 1e3, explain: 'Yield-line collapse load of the plate under the contact patch before membrane action develops.' },
        { key: 'Eplate', label: 'Plate energy capacity', unit: 'kJ', value: pl.energy / 1e3, explain: `Plastic work absorbed by the plate up to the limiting deflection of ${fmt(dLim * 1e3)} mm (average membrane strain ${inp.epsCr} %).` },
        { key: 'Ebeam', label: 'Beam energy capacity', unit: 'kJ', value: Ebeam / 1e3, explain: 'Plastic work of the supporting beam mechanism up to its rotation limit (zero when the beam is not counted).' },
        { key: 'ucEnergy', label: 'Impact energy over total capacity', unit: '-', value: round(E / Ecap, 5), limit: 1, sense: 'max', explain: 'Fraction of the available plastic energy absorption used by the impact. Above 1.0 the object is not stopped.' },
        { key: 'dPlate', label: 'Plate deflection under the object', unit: 'mm', value: dAct * 1e3, limit: dLim * 1e3, sense: 'max', explain: 'Permanent dish of the plate needed to absorb the energy (the limit value when the plate is predicted to tear).' },
        { key: 'wBeam', label: 'Beam plastic deflection', unit: 'mm', value: wBeam * 1e3, explain: 'Mid-span set of the supporting beam needed to absorb the energy the plate could not.' },
        { key: 'hAllow', label: 'Allowable drop height for this mass', unit: 'm', value: hAllow, limit: inp.height, sense: 'min', explain: 'Height from which this mass can fall and still be stopped by the deck: the lifting height restriction over this area.' },
        { key: 'tReq', label: 'Plate thickness to absorb the energy alone', unit: 'mm', value: tReq * 1e3, explain: 'Thickness at which the plate panel alone absorbs the impact within the strain limit (equal to the actual thickness when it already does).' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Plate resistance and absorbed energy', x: { label: 'Plate deflection', unit: 'mm' }, y: { label: 'Resistance (kN) / energy (kJ)', unit: 'kN, kJ' },
          series: [{ name: 'Resistance', x: ds.map((d) => round(d * 1e3)), y: ds.map((d) => round(plate(t, Math.max(d, 1e-9)).resistance / 1e3)) }, { name: 'Absorbed energy', x: ds.map((d) => round(d * 1e3)), y: ds.map((d) => round(plate(t, Math.max(d, 1e-9)).energy / 1e3)) }, { name: 'Impact energy', x: [0, round(dLim * 1e3)], y: [round(E / 1e3), round(E / 1e3)], dash: true }],
          explain: 'Punch resistance of the plate (flat bending plateau, then the rising membrane branch) and the energy absorbed, up to the strain-limited deflection. The plate stops the object where the energy curve crosses the dashed impact energy.' },
        { type: 'line', title: 'Allowable drop height against mass', x: { label: 'Mass of the object', unit: 't' }, y: { label: 'Drop height', unit: 'm' },
          series: [{ name: 'Plate and beam', x: masses.map((v) => round(v)), y: masses.map((v) => round(Ecap / (v * 1e3 * G))) }, { name: 'Plate alone (no tearing)', x: masses.map((v) => round(v)), y: masses.map((v) => round(pl.energy / (v * 1e3 * G))) }, { name: 'Design case', x: [inp.mass], y: [inp.height], mode: 'markers' }],
          explain: 'Heights from which a given mass can fall and still be stopped: use it as the lifting restriction chart for this deck area. The design case must lie below the curve relied upon.' },
        { type: 'bar', title: 'Energy balance', categories: ['Impact energy', 'Plate capacity', 'Beam capacity'], series: [{ name: 'Energy', y: [round(E / 1e3), round(pl.energy / 1e3), round(Ebeam / 1e3)] }], y: { label: 'Energy', unit: 'kJ' }, explain: 'Impact energy beside the energy each structural element can dissipate plastically; the capacities add when the beam is counted.' },
      ],
      tables: [{ title: 'Impact summary', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
        rows: [{ q: 'Fall time', v: round(drop.time), u: 's' }, { q: 'Momentum', v: round(m * drop.velocity / 1e3), u: 'kN·s' }, { q: 'Limiting plate deflection', v: round(dLim * 1e3), u: 'mm' }, { q: 'Deflection at which membrane action takes over', v: round(pl.transitionDeflection * 1e3), u: 'mm' }, { q: 'Plate resistance at the limit', v: round(pl.resistance / 1e3), u: 'kN' },
          { q: 'Beam mechanism load', v: round(R0 / 1e3), u: 'kN' }, { q: 'Beam deflection at the rotation limit', v: round(wLim * 1e3), u: 'mm' }, { q: 'Allowable drop height, plate alone', v: round(hAllowPlate), u: 'm' }],
        explain: 'Supporting numbers of the energy balance: the resistance levels can be compared with the capacity of equipment or piping directly under the plate.' }],
      notes: [`Impact energy ${fmt(E / 1e3)} kJ against ${fmt(pl.energy / 1e3)} kJ in the plate and ${fmt(Ebeam / 1e3)} kJ in the beam.`, plateOk ? `The plate stops the object with a ${fmt(dAct * 1e3)} mm dish.` : `The plate is predicted to tear; ${ok ? 'the beam absorbs the rest' : 'the deck does not stop the object'}.`],
      advice: [
        ok ? adv('engineering', plateOk ? (checks[0].uc > 0.85 ? 'warn' : 'ok') : 'warn', plateOk ? `The ${inp.tPlate} mm plate stops the ${inp.mass} t object from ${inp.height} m with a ${fmt(dAct * 1e3)} mm dish (${fmt(checks[0].uc * 100)} % of its energy capacity).` : `The plate tears (energy ${fmt(checks[0].uc)} times its capacity) but the supporting beam stops the object with ${fmt(wBeam * 1e3)} mm of plastic deflection.`,
          plateOk ? `Restrict lifts over this area to ${fmt(hAllowPlate, 2)} m for this mass, and plan a plate repair after any drop.` : `Do not route this lift over items directly below the plate; either limit the height to ${fmt(hAllowPlate, 2)} m or thicken the plate to ${fmt(Math.ceil(tReq * 1e3), 2)} mm in the drop zone.`,
          'DNV-RP-C204 (2019) §5; NORSOK N-004 (2022) Annex A.4', `Local protection is light: a ${inp.bPanel} mm square panel at ${inp.tPlate} mm weighs ${fmt(panelT * 1e3)} kg, so designated drop zones with thicker plate cost far less steel than a thicker deck overall (each extra mm over 100 m² is ${steelCo2Text(100 * 7.85 / 1e3)}).`)
          : adv('engineering', 'fail', `The deck does not stop the ${inp.mass} t object from ${inp.height} m: impact energy ${fmt(E / 1e3)} kJ against a capacity of ${fmt(Ecap / 1e3)} kJ.`,
            `Limit the lifting height to ${fmt(hAllow, 2)} m for this mass, or thicken the plate in the drop zone to about ${fmt(Math.ceil(tReq * 1e3), 2)} mm, or fit a sacrificial protection deck or energy-absorbing grating over the area.`,
            'DNV-RP-C204 (2019) §5; NORSOK N-004 (2022) Annex A.4', `Thickening a 100 m² drop zone from ${inp.tPlate} to ${fmt(Math.ceil(tReq * 1e3), 2)} mm adds about ${steelCo2Text(100 * Math.max(tReq - t, 0.001) * st.rho / 1e3)} A lifting restriction costs no material.`),
        adv('economic', 'info', `The allowable drop height for ${inp.mass} t is ${fmt(hAllow, 3)} m (plate alone ${fmt(hAllowPlate, 3)} m): every lift above that needs protection below or a shutdown of the equipment underneath.`,
          'Use the height–mass chart to set crane operating restrictions and laydown-area protection early; a protected laydown zone is paid once, while shutting in production for unprotected heavy lifts is paid at every lift.', 'DNV-ST-N001 (2021) lifting over live plant; risk-based decision practice per ISO 17776:2016'),
      ],
    };
  },
  verify: [
    { // E = m g h = 5000·9.80665·10 = 490.33 kJ; v = √(2·9.80665·10) = 14.005 m/s.
      name: 'Free fall in air: E = m g h, v = √(2 g h)', source: 'Newtonian mechanics (hand calculation)', inputs: { mass: 5, height: 10 }, expect: [{ kpi: 'energy', value: 490.33, tol: 1e-4 }, { kpi: 'velocity', value: 14.005, tol: 1e-4 }] },
    { // Clamped circular plate, rigid punch: Pc = 4π mp Rp/(Rp − a), mp = fy t²/4 = 355e6·0.01²/4 = 8875 N·m/m;
      // Rp = 0.35 m, a = 0.15 m → Pc = 4π·8875·0.35/0.20 = 195.2 kN.
      name: 'Yield-line collapse load of a clamped circular plate under a rigid punch', source: 'Upper-bound theorem (Johansen yield-line theory); Jones, Structural Impact, §2.6',
      inputs: { tPlate: 10, bPanel: 700, dContact: 300, edge: 'clamped' }, expect: [{ kpi: 'Pc', value: 195.2, tol: 0.002 }] },
  ],
};
