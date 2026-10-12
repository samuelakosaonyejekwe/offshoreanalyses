// Accidental impact on a marine structure: vessel collision with a bridge pier, tower or dolphin (AASHTO and
// Eurocode equivalent static forces, energy balance, dynamic amplification, protection sizing) and floating
// debris impact.
import { G, CO2_STEEL, CO2_CONC, need, finite, r, f, linspace, check, eng, eco, slrInputs, slrRows, SLR_BASIS } from './_common.js';

/** Peak dynamic load factor of an undamped single-degree-of-freedom system under a rectangular pulse of duration td. */
const dafRect = (ratio) => (ratio >= 0.5 ? 2 : 2 * Math.sin(Math.PI * ratio));

export default {
  id: 's07-vessel-impact', suite: 7, title: 'Vessel and debris impact',
  summary: 'Collision force and energy of a vessel striking a pier or tower, the dynamic amplification, the demand on the structure and the protection needed to absorb it; plus floating-debris impact.',
  covers: [11], fidelity: 'screening',
  standards: ['AASHTO LRFD Bridge Design Specifications (9th ed., 2020) §3.14', 'AASHTO Guide Specification for Vessel Collision Design of Highway Bridges (2nd ed., 2009)', 'EN 1991-1-7:2006 §4.6 and Annex C', 'ASCE 7-22 §5.4.5 and C5.4.5 (debris impact)', 'PIANC WG 33 Guidelines for the Design of Fender Systems (2002)', 'DNV-RP-C204 (2019) (accidental loads, energy dissipation)'],
  models: ['Conservation of momentum and kinetic energy', 'AASHTO ship collision force (Woisin)', 'Eurocode hard-impact force F = v √(k m)', 'Impulse–momentum debris impact', 'Single-degree-of-freedom dynamic load factor', 'Plastic energy absorption', 'Beam statics for the pier demand'],
  limits: 'Equivalent static forces from head-on impact of a single vessel on a rigid, fixed-base pier; oblique and glancing blows, vessel–pier interaction, local crushing of the pier section, pile-group redistribution and collision probability (annual frequency of collapse) are not computed. Explicit nonlinear finite-element collision analysis is outside the module.',
  publishes: ['marine.impactForce', 'marine.impactEnergy'],
  searchTerms: 'vessel collision bridge pier impact force AASHTO Eurocode fender energy debris impact',
  inputs: [
    { key: 'dwt', label: 'Vessel deadweight tonnage', unit: 't', type: 'number', default: 8000, min: 50, max: 400000, group: 'Vessel', role: 'bc', help: 'Design vessel from the fleet using the waterway.' },
    { key: 'disp', label: 'Vessel displacement', unit: 't', type: 'number', default: 11000, min: 50, max: 500000, group: 'Vessel', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 } },
    { key: 'V', label: 'Impact speed', unit: 'm/s', type: 'number', default: 3.0, min: 0.1, max: 12, group: 'Vessel', role: 'ic', uncertain: { dist: 'lognormal', cov: 0.25 }, help: 'Initial condition of the collision: transit speed in the channel, reducing to the drift speed away from it (AASHTO Fig. 3.14.6-1).' },
    { key: 'CH', label: 'Hydrodynamic mass coefficient', unit: '-', type: 'number', default: 1.05, min: 1, max: 1.5, group: 'Vessel', role: 'criteria', help: 'AASHTO: 1.05 with large under-keel clearance, 1.25 when it is less than 10 % of draught.' },
    { key: 'kShip', label: 'Equivalent vessel stiffness (Eurocode)', unit: 'MN/m', type: 'number', default: 15, min: 1, max: 100, group: 'Vessel', role: 'criteria', help: 'EN 1991-1-7 Annex C.4.3: 5 MN/m inland vessels, 15 MN/m sea-going vessels.' },
    { key: 'zImp', label: 'Height of impact above the pier base', unit: 'm', type: 'number', default: 12, min: 0.5, max: 80, group: 'Structure', role: 'geometry', help: 'Bow contact level above the foundation at the present design water level.' },
    ...slrInputs({ group: 'Structure' }),
    { key: 'Rlat', label: 'Lateral resistance of the pier at impact level', unit: 'MN', type: 'number', default: 45, min: 0.1, max: 2000, group: 'Structure', role: 'criteria', design: true, help: 'Plastic or ultimate lateral capacity from the structural and foundation design (extreme event limit state).' },
    { key: 'Mres', label: 'Moment resistance at the base', unit: 'MN·m', type: 'number', default: 650, min: 1, max: 1e5, group: 'Structure', role: 'criteria', design: true },
    { key: 'Tn', label: 'Natural period of the pier', unit: 's', type: 'number', default: 0.6, min: 0.02, max: 10, group: 'Structure', role: 'geometry' },
    { key: 'support', label: 'Pier head condition', type: 'select', default: 'free', options: [{ value: 'free', label: 'Free-standing (cantilever from the base)' }, { value: 'propped', label: 'Restrained by the superstructure at the top' }], group: 'Structure', role: 'bc' },
    { key: 'Hp', label: 'Pier height (base to top restraint)', unit: 'm', type: 'number', default: 30, min: 1, max: 200, group: 'Structure', role: 'geometry' },
    { key: 'protect', label: 'Protection system', type: 'select', default: 'none', options: [{ value: 'none', label: 'None (pier takes the impact)' }, { value: 'fender', label: 'Fender / sacrificial dolphin' }], group: 'Protection', role: 'geometry' },
    { key: 'Ffend', label: 'Protection crushing force', unit: 'MN', type: 'number', default: 20, min: 0.1, max: 500, group: 'Protection', role: 'geometry', showIf: { key: 'protect', in: ['fender'] }, design: true },
    { key: 'stroke', label: 'Protection crush stroke available', unit: 'm', type: 'number', default: 3, min: 0.05, max: 30, group: 'Protection', role: 'geometry', showIf: { key: 'protect', in: ['fender'] }, design: true },
    { key: 'mDeb', label: 'Debris mass', unit: 'kg', type: 'number', default: 450, min: 10, max: 100000, group: 'Debris', role: 'bc', help: 'ASCE 7: 450 kg log or pole as a minimum; shipping containers 2–30 t.' },
    { key: 'vDeb', label: 'Debris speed (flow velocity)', unit: 'm/s', type: 'number', default: 3, min: 0.1, max: 15, group: 'Debris', role: 'ic' },
    { key: 'tDeb', label: 'Debris impact duration', unit: 's', type: 'number', default: 0.03, min: 0.005, max: 1, group: 'Debris', role: 'criteria', help: 'ASCE 7-22 C5.4.5 recommends 0.03 s.' },
    { key: 'Rlocal', label: 'Local resistance to debris impact', unit: 'kN', type: 'number', default: 400, min: 1, max: 1e6, group: 'Debris', role: 'criteria' },
  ],
  compute(inp) {
    const m = inp.disp * 1000, V = inp.V;
    const KE = 0.5 * inp.CH * m * V * V;                                  // J  (AASHTO Eq. 3.14.7-1: 500 CH W V²)
    const Ps = 1.2e5 * V * Math.sqrt(inp.dwt);                            // N  (AASHTO Eq. 3.14.8-1, SI form)
    const aB = 1.54 * KE / Ps;                                            // m  (AASHTO Eq. 3.14.9-1, SI form)
    const Fec = V * Math.sqrt(inp.kShip * 1e6 * m);                       // N  (EN 1991-1-7 Eq. C.1)
    const Fdes = Math.max(Ps, Fec);
    const td = inp.CH * m * V / Fdes;                                     // s  momentum over force: equivalent rectangular pulse
    const daf = dafRect(td / inp.Tn);
    finite(Fdes, 'Collision force');
    // protection: absorbs Ffend × stroke; the pier then sees the fender force while it crushes, and the residual energy if it bottoms out
    const fender = inp.protect === 'fender';
    const Eabs = fender ? inp.Ffend * 1e6 * inp.stroke : 0, Eres = Math.max(0, KE - Eabs);
    const strokeReq = fender ? KE / (inp.Ffend * 1e6) : 0;
    const Fpier = fender ? (Eres > 0 ? Fdes * Math.sqrt(Eres / KE) : inp.Ffend * 1e6) : Fdes; // residual force scales with √energy (F ∝ V)
    const demand = (rise) => {
      const a = Math.min(inp.zImp + rise, inp.Hp), L = inp.Hp;
      // fixed base; free top: M = F a.  Propped top (point load at a from the fixed end): M_base = F a b (L + b)/(2 L²)
      const b = L - a;
      return { a, M: inp.support === 'propped' ? Fpier * a * b * (L + b) / (2 * L * L) : Fpier * a, Vb: inp.support === 'propped' ? Fpier * b * (3 * L * L - b * b) / (2 * L ** 3) : Fpier };
    };
    need(inp.zImp < inp.Hp, 'The impact level must lie below the top of the pier.');
    const d = demand(inp.slr);
    const ucF = Fpier / (inp.Rlat * 1e6), ucM = d.M / (inp.Mres * 1e6);
    // debris: half-sine impulse, F = π m v / (2 Δt)
    const Fdeb = Math.PI * inp.mDeb * inp.vDeb / (2 * inp.tDeb), ucDeb = Fdeb / (inp.Rlocal * 1e3);
    const ok = ucF <= 1 && ucM <= 1, okDeb = ucDeb <= 1;
    const scen = slrRows(inp).map((s) => { const q = demand(s.rise); return { name: s.name, rise: s.rise, a: r(q.a, 3), M: r(q.M / 1e6), uc: r(q.M / (inp.Mres * 1e6), 3) }; });
    const vs = linspace(0.25, Math.max(2 * V, 2), 30);
    // steel that must yield to absorb the whole energy: σy ε / ρ ≈ 355 MPa × 10 % / 7850 kg/m³ ≈ 5 kJ/kg (theoretical minimum, screening)
    const steelT = KE / 5e3 / 1000;
    const vAllow = V / Math.max(Math.max(ucF, ucM), 1e-9);

    const advice = [
      ok
        ? eng('ok', `Design collision force is ${f(Fpier / 1e6)} MN (AASHTO ${f(Ps / 1e6)}, Eurocode ${f(Fec / 1e6)} MN; energy ${f(KE / 1e6)} MJ): utilisation ${f(ucF, 2)} in shear and ${f(ucM, 2)} in base moment.`,
          `The pier resists the design vessel at ${f(V)} m/s; the limiting speed is about ${f(vAllow, 2)} m/s. Confirm the design vessel and speed from traffic data and check the superstructure connections for the prop reaction.`,
          'AASHTO LRFD (2020) §3.14.8–3.14.14; EN 1991-1-7:2006 §4.6 and Annex C', 'A pier that survives the design collision by its own strength needs no separate protection structures in the waterway, avoiding their steel, piling and seabed footprint.')
        : eng('fail', `The pier does not resist the design collision: force ${f(Fpier / 1e6)} MN against ${f(inp.Rlat)} MN (utilisation ${f(ucF, 2)}), base moment ${f(d.M / 1e6)} MN·m against ${f(inp.Mres)} MN·m (${f(ucM, 2)}).`,
          `Either strengthen to ${f(Fpier / 1e6)} MN and ${f(d.M / 1e6)} MN·m, or protect: an independent dolphin or fender absorbing ${f(KE / 1e6)} MJ (for example ${f(Math.min(inp.Rlat, Fdes / 1e6) * 0.8)} MN crushing force over ${f(KE / (Math.min(inp.Rlat * 1e6, Fdes) * 0.8), 2)} m), or an artificial island; a speed limit of ${f(vAllow, 2)} m/s would also bring the pier within capacity.`,
          'AASHTO LRFD (2020) §3.14.14–3.14.15 (protection); EN 1991-1-7:2006 Annex C; PIANC WG 33 (2002)',
          `Absorbing ${f(KE / 1e6)} MJ plastically takes at least ${f(steelT, 3)} t of yielding steel in a sacrificial dolphin (≈ ${f(steelT * CO2_STEEL, 3)} t CO₂e; theoretical minimum at about 5 kJ/kg of plastic work — real dolphins need several times more; screening factors); a rock island does it with low-carbon fill but a large seabed footprint.`),
      eng('info', `Impact lasts about ${f(td, 2)} s against a pier period of ${f(inp.Tn)} s: dynamic load factor ${f(daf, 2)} for an undamped rectangular pulse; bow crush depth ${f(aB, 2)} m.`,
        daf > 1.5 ? 'The pulse is long compared with the period, so the pier responds almost statically to the crushing force — the equivalent static force already represents this; do not add the dynamic factor on top of the AASHTO force, but do apply it in a time-history check.' : 'The pulse is short compared with the pier period: the response is impulsive and a dynamic analysis will give a lower demand than the equivalent static force.',
        'Biggs (1964) Introduction to Structural Dynamics ch. 2; EN 1991-1-7 Annex C.2'),
      eng(okDeb ? 'info' : 'fail', `Debris impact: ${f(inp.mDeb)} kg at ${f(inp.vDeb)} m/s gives ${f(Fdeb / 1e3)} kN over ${f(inp.tDeb * 1000)} ms (utilisation ${f(ucDeb, 2)} of the local resistance). Sea-level rise lifts the impact level by ${f(inp.slr)} m, raising the base moment to ${f(d.M / 1e6)} MN·m.`,
        okDeb ? 'Local elements withstand the debris strike; keep parapets and services out of the debris zone at the raised water level.' : `Strengthen the exposed element to ${f(Fdeb / 1e3)} kN, add a debris deflector or sacrificial fender, or raise the element above the flood level.`, `ASCE 7-22 §5.4.5; ${SLR_BASIS}`),
      eco(ok ? 'info' : 'warn', `The collision energy of ${f(KE / 1e6)} MJ must go somewhere: into the pier (strengthening), a protection structure (at least ${f(steelT, 3)} t of fully yielding steel, in practice several times more) or be prevented (speed limits, traffic management, pilotage).`,
        'Rank the options by annualised cost per unit of risk removed: operational controls are cheapest, protection structures next, strengthening of an existing pier usually the most expensive. Use the annual frequency of collapse to justify the level of spend.',
        'AASHTO LRFD (2020) §3.14.5 (annual frequency of collapse); ALARP cost–benefit practice; whole-life costing per ISO 15686-5',
        `Each MN·m of extra base moment capacity costs foundation concrete and steel (${f(CO2_CONC)} t CO₂e/m³ and ${f(CO2_STEEL)} t CO₂e/t, screening factors); preventing the collision has no embodied carbon.`),
    ];
    return {
      kpis: [
        { key: 'KE', label: 'Collision kinetic energy', value: r(KE / 1e6), unit: 'MJ', explain: 'Energy of the vessel with its hydrodynamic added mass at the impact speed. It must be dissipated by crushing of the bow, the protection and the structure.' },
        { key: 'Ps', label: 'AASHTO equivalent static force', value: r(Ps / 1e6), unit: 'MN', explain: 'Head-on ship collision force from Woisin\'s tests, 1.2×10⁵ V √DWT: the mean crushing force of the bow, used as an equivalent static load.' },
        { key: 'Fec', label: 'Eurocode hard-impact force', value: r(Fec / 1e6), unit: 'MN', explain: 'F = v √(k m) with the code\'s equivalent vessel stiffness, assuming a rigid structure that absorbs no energy.' },
        { key: 'Fpier', label: 'Design force on the pier', value: r(Fpier / 1e6), unit: 'MN', limit: inp.Rlat, sense: 'max', explain: 'Larger of the two code forces, reduced where a protection system absorbs energy first. Compared with the lateral resistance at impact level.' },
        { key: 'Mbase', label: 'Base moment from the impact', value: r(d.M / 1e6), unit: 'MN·m', limit: inp.Mres, sense: 'max', explain: 'Moment at the pier base from the design force at the impact level, for the stated head restraint. Rises with the water level.' },
        { key: 'aB', label: 'Bow crush depth', value: r(aB, 3), unit: 'm', explain: 'Length of bow damaged in the collision (AASHTO). The superstructure and pier cap must be set back by at least this from the vessel\'s path.' },
        { key: 'td', label: 'Impact duration', value: r(td, 3), unit: 's', explain: 'Momentum divided by the design force — the length of the equivalent rectangular force pulse.' },
        { key: 'daf', label: 'Dynamic load factor', value: r(daf, 3), unit: '-', explain: 'Peak dynamic over static response of the pier to the pulse. 2 for long pulses; less than 1 when the impact is short against the natural period.' },
        { key: 'strokeReq', label: 'Protection stroke to absorb the energy', value: r(strokeReq, 3), unit: 'm', explain: 'Crush distance the fender or dolphin needs at its crushing force to absorb all the collision energy. Zero without protection.' },
        { key: 'Eres', label: 'Energy reaching the pier', value: r(Eres / 1e6), unit: 'MJ', explain: 'Collision energy left after the protection has crushed through its available stroke.' },
        { key: 'vAllow', label: 'Impact speed the pier can resist', value: r(vAllow, 3), unit: 'm/s', explain: 'Speed at which the governing utilisation would be 1.0 (force is proportional to speed). A basis for speed limits in the waterway.' },
        { key: 'Fdeb', label: 'Debris impact force', value: r(Fdeb / 1e3), unit: 'kN', limit: inp.Rlocal, sense: 'max', explain: 'Impulse–momentum force of a floating object stopped in the stated duration (half-sine pulse). Governs parapets, screens and slender piles.' },
      ],
      checks: [
        check('Collision force against lateral resistance', Fpier / 1e6, inp.Rlat, 'MN', 'AASHTO LRFD §3.14.14; EN 1991-1-7 §4.6', 'Guards against the pier shearing or the foundation failing laterally at impact. Exceeded, the pier is displaced and the span it carries can drop.'),
        check('Base moment against moment resistance', d.M / 1e6, inp.Mres, 'MN·m', 'AASHTO LRFD §3.14.14 (extreme event II)', 'Guards against a plastic hinge or foundation overturning at the base under the impact.'),
        check('Debris impact against local resistance', Fdeb / 1e3, inp.Rlocal, 'kN', 'ASCE 7-22 §5.4.5', 'Guards against local failure of slender or thin elements struck by floating debris during floods and storm surges.'),
        ...(fender ? [check('Protection stroke', strokeReq, inp.stroke, 'm', 'PIANC WG 33 (2002); AASHTO Guide Specification §7', 'Guards against the fender or dolphin bottoming out before the vessel is stopped, after which the pier takes the remaining energy.')] : []),
      ],
      plots: [
        { type: 'line', title: 'Collision force versus impact speed', x: { label: 'Impact speed', unit: 'm/s' }, y: { label: 'Force', unit: 'MN' },
          series: [{ name: 'AASHTO', x: vs.map((v) => r(v, 3)), y: vs.map((v) => r(1.2e5 * v * Math.sqrt(inp.dwt) / 1e6, 3)) }, { name: 'Eurocode', x: vs.map((v) => r(v, 3)), y: vs.map((v) => r(v * Math.sqrt(inp.kShip * 1e6 * m) / 1e6, 3)) }, { name: 'Pier resistance', x: [vs[0], vs[29]], y: [inp.Rlat, inp.Rlat], dash: true }],
          explain: 'Both code forces grow linearly with speed. The crossing with the resistance line is the highest speed the unprotected pier survives.' },
        { type: 'line', title: 'Collision energy versus impact speed', x: { label: 'Impact speed', unit: 'm/s' }, y: { label: 'Energy', unit: 'MJ' },
          series: [{ name: 'Kinetic energy', x: vs.map((v) => r(v, 3)), y: vs.map((v) => r(0.5 * inp.CH * m * v * v / 1e6, 3)) }, ...(fender ? [{ name: 'Protection capacity', x: [vs[0], vs[29]], y: [r(Eabs / 1e6, 3), r(Eabs / 1e6, 3)], dash: true }] : [])],
          explain: 'Energy to be absorbed rises with the square of speed, so halving the speed removes three quarters of the problem. The dashed line is what the protection can take.' },
      ],
      tables: [
        { title: 'Sea-level-rise scenarios: impact level and base moment', columns: [{ key: 'name', label: 'Scenario' }, { key: 'rise', label: 'Rise', unit: 'm' }, { key: 'a', label: 'Impact height', unit: 'm' }, { key: 'M', label: 'Base moment', unit: 'MN·m' }, { key: 'uc', label: 'Utilisation', unit: '-' }],
          rows: scen, explain: 'A higher sea lifts the point of impact up the pier, lengthening the lever arm of a free-standing pier. Check the row for the end of the design life.' },
        { title: 'Impact summary', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['Vessel mass with added mass', inp.CH * m / 1000, 't'], ['Momentum', inp.CH * m * V / 1e6, 'MN·s'], ['Base shear', d.Vb / 1e6, 'MN'], ['Energy absorbed by protection', Math.min(Eabs, KE) / 1e6, 'MJ'], ['Pulse duration / pier period', td / inp.Tn, '-'], ['Debris momentum', inp.mDeb * inp.vDeb, 'N·s']].map(([q, v, u]) => ({ q, v: r(v), u })),
          explain: 'Momentum and energy bookkeeping of the collision, for comparison with a hand check or a detailed collision analysis.' },
      ],
      notes: [`${Ps >= Fec ? 'AASHTO' : 'Eurocode'} force governs; ${fender ? (Eres > 0 ? 'the protection bottoms out and the pier takes the residual energy.' : 'the protection absorbs the whole collision.') : 'no protection: the pier takes the full force.'}`, `${ucM >= ucF ? 'Base moment' : 'Lateral force'} is the governing demand on the pier.`],
      advice,
      publish: { 'marine.impactForce': Fpier, 'marine.impactEnergy': KE },
      quantities: { steel_t: r(ok ? 0 : 3 * steelT) },
    };
  },
  verify: [
    { name: 'AASHTO ship collision: PS = 1.2×10⁵ V √DWT, KE = 500 CH W V², aB = 1.54 KE/PS', source: 'AASHTO LRFD (2020) Eq. 3.14.8-1, 3.14.7-1, 3.14.9-1 (SI). DWT = 10 000 t, W = 13 000 t, V = 4 m/s, CH = 1.05: PS = 1.2e5 × 4 × 100 = 48.0 MN; KE = 500 × 1.05 × 13 000 × 16 = 109.2 MJ; aB = 1.54 × 109.2/48.0 = 3.503 m',
      inputs: { dwt: 10000, disp: 13000, V: 4, CH: 1.05 }, expect: [{ kpi: 'Ps', value: 48.0, tol: 0.001 }, { kpi: 'KE', value: 109.2, tol: 0.001 }, { kpi: 'aB', value: 3.503, tol: 0.002 }] },
    { name: 'Eurocode hard impact F = v √(k m) and debris impulse F = π m v/(2 Δt)', source: 'EN 1991-1-7 Eq. C.1: v = 4 m/s, k = 15 MN/m, m = 13 000 t → F = 4 × √(15e6 × 13e6) = 55.86 MN. Half-sine impulse: 450 kg at 3 m/s in 0.03 s → π × 1350/0.06 = 70.69 kN',
      inputs: { dwt: 10000, disp: 13000, V: 4, kShip: 15, mDeb: 450, vDeb: 3, tDeb: 0.03 }, expect: [{ kpi: 'Fec', value: 55.86, tol: 0.001 }, { kpi: 'Fdeb', value: 70.69, tol: 0.001 }] },
  ],
};
