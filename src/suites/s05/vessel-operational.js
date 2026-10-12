// Vessel-induced operational loading at a berth: the drawdown and return current of a passing ship and the surge
// force they put on the moored vessel, ship wash waves on the piles, propeller and thruster jet loads on piles and
// on the quay face, the number of load events per year, and the passing speed the moorings can tolerate.
import { pileForceClosedForm } from '../../engine/hydro/morison.js';
import { wavelength } from '../../engine/metocean/waves.js';
import { propellerJet, returnCurrent, washWave } from './_vessel.js';
import { G, need, r, f, check, eng, eco, linspace, maxWhere, CONFIRM, ECON_BASIS, CO2_MGO } from './_common.js';

export default {
  id: 's05-vessel-operational', suite: 5, title: 'Vessel-induced operational loading',
  summary: 'Loads that ship traffic and manoeuvring put on the berth in service: passing-ship drawdown and surge force on the moored vessel, wash-wave and propeller-jet forces on the piles, thruster jet force on the quay face, and the permissible passing speed.',
  covers: [13], fidelity: 'screening',
  standards: ['PIANC WG180 (2015) scour and loads from ship propulsion', 'PIANC MarCom WG 121 Harbour Approach Channels Design Guidelines (2014)', 'BAW (2010) Principles for the Design of Bank and Bottom Protection for Inland Waterways (GBB)', 'OCIMF Mooring Equipment Guidelines 4th ed. (2018) §3 (passing ships)', 'BS 6349-1-2:2016 §30 (vessel-induced actions)'],
  models: ['One-dimensional return current and drawdown of a passing ship (Schijf energy method)', 'Surge force from the water-surface slope (hydrostatic slope method)', 'Ship wash wave height and period (Kelvin wake with empirical coefficient)', "Morison's equation on piles", 'Propeller jet: efflux, Gaussian decay, drag on a pile and momentum force on a wall'],
  limits: 'Screening level. The passing-ship force uses a one-dimensional drawdown and a water-surface slope over a ramp length entered by the user; it does not resolve the lateral separation, the sway force and yaw moment, or the dynamic response of the moored ship, for which a passing-ship potential-flow model with a dynamic mooring analysis is needed. Wash and jet coefficients are indicative and editable. CFD of the jet and of ship waves is outside the module.',
  publishes: ['wharf.passingSurgeForce', 'wharf.washForcePile', 'wharf.jetForcePile'],
  searchTerms: 'passing ship moored vessel forces ship wash propeller jet load quay',
  inputs: [
    { key: 'h', label: 'Water depth in the channel and berth', unit: 'm', type: 'number', default: 14, min: 3, max: 40, group: 'Waterway', role: 'bc', bind: 'site.depth' },
    { key: 'bw', label: 'Waterway width at the water surface', unit: 'm', type: 'number', default: 350, min: 30, max: 5000, group: 'Waterway', role: 'geometry', help: 'Width of the basin or channel in which the passing ship displaces water (quay to opposite bank).' },
    { key: 'rhoW', label: 'Water density', unit: 'kg/m³', type: 'number', default: 1025, min: 995, max: 1035, group: 'Waterway', role: 'environment', bind: 'site.rhoW' },
    { key: 'Vp', label: 'Speed of the passing ship through the water', unit: 'm/s', type: 'number', default: 3.6, min: 0.3, max: 12, group: 'Passing ship', role: 'bc', design: true, uncertain: { dist: 'normal', cov: 0.15 }, help: '3.6 m/s is 7 knots. The force grows with roughly the square of the speed.' },
    { key: 'Bp', label: 'Beam of the passing ship', unit: 'm', type: 'number', default: 43, min: 5, max: 65, group: 'Passing ship', role: 'geometry' },
    { key: 'Tp', label: 'Draft of the passing ship', unit: 'm', type: 'number', default: 12.5, min: 1, max: 24, group: 'Passing ship', role: 'geometry' },
    { key: 'Lp', label: 'Length of the passing ship', unit: 'm', type: 'number', default: 300, min: 30, max: 420, group: 'Passing ship', role: 'geometry' },
    { key: 'cm', label: 'Midship section coefficient', unit: '-', type: 'number', default: 0.97, min: 0.7, max: 1, group: 'Passing ship', role: 'geometry' },
    { key: 'sPass', label: 'Passing distance, track to berth line', unit: 'm', type: 'number', default: 120, min: 15, max: 2000, group: 'Passing ship', role: 'bc' },
    { key: 'kBerth', label: 'Share of the drawdown reaching the berth', unit: '-', type: 'number', default: 0.7, min: 0.1, max: 1, group: 'Passing ship', role: 'criteria', help: 'The one-dimensional method gives the section-average depression; it is larger beside the hull and smaller at the far bank. 1.0 is the conservative choice for a narrow waterway.' },
    { key: 'ramp', label: 'Length over which the water level falls (fraction of the passing ship length)', unit: '-', type: 'number', default: 0.35, min: 0.1, max: 1, group: 'Passing ship', role: 'criteria', help: `The level drops from the bow wave to the trough over part of the fore body; the slope sets the force on the moored ship. 0.25–0.5 is a screening range; calibrate against a passing-ship model or measurements. ${CONFIRM}` },
    { key: 'dispM', label: 'Displacement of the moored ship', unit: 't', type: 'number', default: 80000, min: 500, max: 600000, group: 'Moored ship', role: 'ic', bind: 'vessel.disp' },
    { key: 'Lm', label: 'Length of the moored ship', unit: 'm', type: 'number', default: 240, min: 30, max: 420, group: 'Moored ship', role: 'geometry' },
    { key: 'surgeCap', label: 'Surge restraint capacity of the mooring', unit: 'kN', type: 'number', default: 1500, min: 50, max: 20000, group: 'Moored ship', role: 'criteria', help: 'Sum of the spring-line working limits resolved along the berth (from the mooring module), for example three springs at 50 % of 1300 kN at 30° to the berth line.' },
    { key: 'alphaW', label: 'Wash coefficient α', unit: '-', type: 'number', default: 1, min: 0.2, max: 1.5, group: 'Wash waves', role: 'criteria', help: `H = α h (s/h)^(−1/3) Fr⁴. About 1.0 tugs, patrol boats and loaded conventional vessels, 0.5 empty push units, 0.35 empty conventional vessels. ${CONFIRM}` },
    { key: 'Vw', label: 'Speed of the wash-generating vessel', unit: 'm/s', type: 'number', default: 5.5, min: 0.5, max: 15, group: 'Wash waves', role: 'bc', help: 'Fast craft (tugs, pilot boats, ferries) usually govern wash; they may differ from the large passing ship.' },
    { key: 'sWash', label: 'Distance from that vessel to the berth', unit: 'm', type: 'number', default: 80, min: 5, max: 2000, group: 'Wash waves', role: 'bc' },
    { key: 'D', label: 'Pile diameter', unit: 'm', type: 'number', default: 1.067, min: 0.2, max: 3, group: 'Piles', role: 'geometry', bind: 'wharf.pileD' },
    { key: 'Cd', label: 'Drag coefficient', unit: '-', type: 'number', default: 1.05, min: 0.5, max: 2, group: 'Piles', role: 'criteria' },
    { key: 'Cm', label: 'Inertia coefficient', unit: '-', type: 'number', default: 1.6, min: 1, max: 2.2, group: 'Piles', role: 'criteria' },
    { key: 'pileCap', label: 'Lateral design resistance of one pile for local loads', unit: 'kN', type: 'number', default: 300, min: 10, max: 10000, group: 'Piles', role: 'criteria', bind: 'wharf.pileHult', help: 'From the pile capacity module, with its resistance factor.' },
    { key: 'P', label: 'Propulsion or thruster power applied', unit: 'kW', type: 'number', default: 2500, min: 50, max: 90000, group: 'Propeller and thruster jets', role: 'operation', design: true, help: 'Power actually used at the berth: a bow thruster at full power, or the main engine share used when unberthing.' },
    { key: 'Dp', label: 'Propeller or thruster diameter', unit: 'm', type: 'number', default: 2.8, min: 0.5, max: 11, group: 'Propeller and thruster jets', role: 'geometry' },
    { key: 'ducted', label: 'Ducted (thruster tunnel or nozzle)', type: 'bool', default: true, group: 'Propeller and thruster jets', role: 'geometry' },
    { key: 'cV0', label: 'Efflux velocity coefficient', unit: '-', type: 'number', default: 1.17, min: 0.8, max: 1.8, group: 'Propeller and thruster jets', role: 'criteria', help: `V0 = c (P/(ρ Dp²))^(1/3): about 1.15 open, 1.17 ducted. ${CONFIRM}` },
    { key: 'xJet', label: 'Distance from the jet outlet to the pile or wall', unit: 'm', type: 'number', default: 12, min: 1, max: 200, group: 'Propeller and thruster jets', role: 'bc' },
    { key: 'calls', label: 'Ship calls per year at this berth', unit: '1/yr', type: 'number', default: 150, min: 1, max: 5000, group: 'Traffic', role: 'operation' },
    { key: 'passes', label: 'Ships passing per year while the berth is occupied', unit: '1/yr', type: 'number', default: 2500, min: 0, max: 200000, group: 'Traffic', role: 'operation' },
  ],
  compute(inp) {
    const rho = inp.rhoW, h = inp.h, Ac = inp.bw * h, Am = inp.cm * inp.Bp * inp.Tp;
    need(inp.Tp < h, 'The draft of the passing ship must be less than the water depth.');
    need(Am < 0.5 * Ac, `The passing ship blocks ${f(Am / Ac * 100)} % of the waterway cross-section: above 50 % the one-dimensional method does not apply.`);
    const sch = returnCurrent({ V: inp.Vp, Ac, Am, bw: inp.bw });
    const zB = inp.kBerth * sch.z, Lramp = inp.ramp * inp.Lp;
    // surge force on the moored ship: weight times the water-surface slope, limited by the share of the hull on the ramp
    const M = inp.dispM * 1e3, share = Math.min(1, Lramp / inp.Lm);
    const surgeOf = (V) => { const q = returnCurrent({ V, Ac, Am, bw: inp.bw }); return M * G * (inp.kBerth * q.z / Lramp) * share; };
    const Fsurge = M * G * (zB / Lramp) * share;
    const Vperm = maxWhere((V) => surgeOf(V) <= inp.surgeCap * 1e3, 0.2, sch.Vlim * 0.99, 30);
    // wash waves at the berth and their force on a pile
    const wash = washWave({ V: inp.Vw, h, s: inp.sWash, alpha: inp.alphaW });
    need(wash.Fr < 1, `The wash-generating vessel is at a depth Froude number of ${f(wash.Fr)}: above 1 the wash relation does not apply (planing or supercritical regime).`);
    const Hw = Math.min(wash.H, 0.78 * h), Tw = Math.max(wash.T, 1);
    const pw = Hw > 1e-4 ? pileForceClosedForm({ H: Hw, T: Tw, depth: h, D: inp.D, Cd: inp.Cd, Cm: inp.Cm, rho }) : { Fmax: 0, Mmax: 0 };
    const steep = Hw > 1e-4 ? Hw / wavelength(Tw, h) : 0;
    // return current drag on a pile over the depth (uniform return flow)
    const Fret = 0.5 * rho * inp.Cd * inp.D * sch.Ur * sch.Ur * h;
    // jet: velocity on the axis at the pile, drag on the pile integrated across the Gaussian jet, and wall force
    const jet = propellerJet({ P: inp.P * 1e3, Dp: inp.Dp, hp: 1, rho, ducted: !!inp.ducted, cV0: inp.cV0 });
    const bRad = 15.4, Vax = jet.axis(inp.xJet);
    // drag ∫ ½ρ Cd D V² dr with V = V_axis exp(−b r²/x²): ∫ exp(−2b r²/x²) dr = x √(π/(2b))
    const jetOf = (P) => { const j = propellerJet({ P, Dp: inp.Dp, hp: 1, rho, ducted: !!inp.ducted, cV0: inp.cV0 }); const v = j.axis(inp.xJet); return 0.5 * rho * inp.Cd * inp.D * v * v * Math.min(inp.xJet * Math.sqrt(Math.PI / (2 * bRad)), h); };
    const Fjet = jetOf(inp.P * 1e3);
    const Pperm = maxWhere((P) => jetOf(P) <= inp.pileCap * 1e3, 1e3, 9e7, 40);
    const Fwall = jet.thrust, pStag = 0.5 * rho * Vax * Vax;
    const FpileMax = Math.max(pw.Fmax, Fjet, Fret);
    const checks = [
      check('Passing-ship surge force against the mooring restraint', Fsurge / 1e3, inp.surgeCap, 'kN', 'OCIMF MEG4 (2018) §3; PIANC WG 121 (2014) §3', 'Force along the berth from the water-level slope of the passing ship against the working capacity of the spring lines. If exceeded the moored ship surges, lines part and loading arms or gangways are at risk.'),
      check('Largest vessel-induced force on a pile', FpileMax / 1e3, inp.pileCap, 'kN', 'BS 6349-1-2:2016 §30; PIANC WG180 (2015) §8', 'Largest of the wash-wave, return-current and jet forces on one pile against its lateral design resistance. These loads recur at every ship movement and add to fatigue at the pile head.'),
      check('Blockage of the waterway by the passing ship', Am / Ac, 0.25, '-', 'PIANC WG 121 (2014) §2 (blockage ratio)', 'Midship section over waterway section. Above about 0.25 the return flow and drawdown rise steeply and passing speeds must be restricted.'),
      check('Passing speed against the limiting speed', inp.Vp, sch.Vlim, 'm/s', 'Schijf (1949) limiting speed', 'Above the limiting speed the water cannot return past the hull: the ship squats heavily and pushes a surge ahead of it that loads every berth it passes.'),
    ];
    const w = checks.slice(0, 2).reduce((a, c) => (c.uc > a.uc ? c : a)), ok = w.uc <= 1;
    const vs = linspace(0.5, Math.min(sch.Vlim * 0.98, 8), 30), xs = linspace(Math.max(2 * inp.Dp, 2), Math.max(4 * inp.xJet, 20), 40);
    return {
      kpis: [
        { key: 'z', label: 'Water-level depression beside the passing ship', value: r(sch.z), unit: 'm', explain: 'Section-average drawdown of the primary wave from continuity and Bernoulli. It travels with the ship and is felt at the berth as a slow fall and rise of the water level.' },
        { key: 'Ur', label: 'Return current', value: r(sch.Ur), unit: 'm/s', explain: 'Mean velocity of the water flowing back past the hull, opposite to the ship. It adds to the tidal current on piles and on the bed.' },
        { key: 'Vlim', label: 'Limiting speed of the passing ship', value: r(sch.Vlim), unit: 'm/s', explain: 'Highest speed at which the displaced water can pass the hull in this cross-section. Ships normally sail below 80–90 % of it.' },
        { key: 'Fsurge', label: 'Surge force on the moored ship', value: r(Fsurge / 1e3), unit: 'kN', limit: inp.surgeCap, sense: 'max', explain: 'Weight of the moored ship times the water-surface slope reaching the berth. It reverses as the ship passes, so both fore and aft springs are loaded in turn.' },
        { key: 'Vperm', label: 'Maximum permissible passing speed', value: r(Vperm), unit: 'm/s', explain: `Passing speed at which the surge force equals the mooring restraint (${f(Vperm / 0.5144)} knots): the basis for a harbour speed limit past the occupied berth.` },
        { key: 'Hwash', label: 'Wash wave height at the berth', value: r(Hw), unit: 'm', explain: `Height of the largest secondary waves from a vessel at ${f(inp.Vw)} m/s (depth Froude number ${f(wash.Fr)}), with a period of about ${f(Tw)} s.` },
        { key: 'Fwash', label: 'Wash-wave force on one pile', value: r(pw.Fmax / 1e3), unit: 'kN', explain: 'Morison force of the wash wave on a pile, from linear theory. Small per event but repeated at every vessel movement.' },
        { key: 'Fret', label: 'Return-current drag on one pile', value: r(Fret / 1e3), unit: 'kN', explain: 'Drag of the return current over the water depth on one pile.' },
        { key: 'V0', label: 'Jet efflux velocity', value: r(jet.V0), unit: 'm/s', explain: 'Velocity of the propeller or thruster jet at its outlet for the power applied.' },
        { key: 'Vjet', label: 'Jet axis velocity at the pile or wall', value: r(Vax), unit: 'm/s', explain: 'Centreline velocity after the jet has spread over the distance entered.' },
        { key: 'Fjet', label: 'Jet drag on one pile', value: r(Fjet / 1e3), unit: 'kN', explain: 'Drag of the jet on a pile crossing its axis, integrated over the Gaussian velocity profile. It acts near the propeller level, often well below the deck.' },
        { key: 'Fwall', label: 'Jet force on a quay face or fender panel', value: r(Fwall / 1e3), unit: 'kN', explain: `Momentum flux of the jet, the upper bound of the force when it strikes a wall square on; the stagnation pressure on the axis is ${f(pStag / 1e3)} kPa.` },
        { key: 'Pperm', label: 'Maximum permissible jet power for the pile', value: r(Pperm / 1e3), unit: 'kW', explain: 'Power at which the jet drag on a pile at this distance reaches its lateral design resistance (reported at 90 MW when not reached).' },
        { key: 'events', label: 'Vessel-induced load events per year', value: r(inp.passes + 2 * inp.calls), unit: '1/yr', explain: 'Passing ships plus an arrival and a departure per call: the cycle count for fatigue of pile heads, fender fixings and mooring lines.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Surge force on the moored ship against passing speed', x: { label: 'Passing speed', unit: 'm/s' }, y: { label: 'Surge force', unit: 'kN' },
          series: [{ name: 'Surge force', x: vs.map((v) => r(v)), y: vs.map((v) => r(surgeOf(v) / 1e3)) }, { name: 'Mooring restraint', x: [r(vs[0]), r(vs[vs.length - 1])], y: [inp.surgeCap, inp.surgeCap], dash: true }],
          explain: 'The force rises faster than the square of the speed as the limiting speed is approached. The crossing with the dashed line is the permissible passing speed.' },
        { type: 'line', title: 'Jet velocity against distance', x: { label: 'Distance from the outlet', unit: 'm' }, y: { label: 'Axis velocity', unit: 'm/s' },
          series: [{ name: 'Jet axis', x: xs.map((x) => r(x)), y: xs.map((x) => r(jet.axis(x))) }, { name: 'Pile or wall', x: [inp.xJet, inp.xJet], y: [0, r(Vax)], dash: true }],
          explain: 'Decay of the jet with distance (inversely proportional beyond the establishment zone). Doubling the stand-off halves the velocity and quarters the load.' },
        { type: 'bar', title: 'Vessel-induced forces on one pile', categories: ['Wash wave', 'Return current', 'Propeller / thruster jet'], series: [{ name: 'Force (kN)', y: [r(pw.Fmax / 1e3), r(Fret / 1e3), r(Fjet / 1e3)] }], x: { label: 'Source', unit: '' }, y: { label: 'Force', unit: 'kN' },
          explain: 'The three ship-related loads on a pile side by side. The jet usually dominates close to thrusters; wash and return current matter by repetition.' },
      ],
      tables: [{ title: 'Vessel-induced loads on the berth', columns: [{ key: 'load', label: 'Load' }, { key: 'value', label: 'Value', unit: 'kN' }, { key: 'on', label: 'Acts on' }, { key: 'n', label: 'Events per year' }],
        rows: [{ load: 'Passing-ship surge', value: r(Fsurge / 1e3), on: 'Spring lines and their bollards', n: r(inp.passes) }, { load: 'Wash wave', value: r(pw.Fmax / 1e3), on: 'Each pile (near the surface)', n: r(inp.passes + 2 * inp.calls) }, { load: 'Return current', value: r(Fret / 1e3), on: 'Each pile (full depth)', n: r(inp.passes) }, { load: 'Jet drag', value: r(Fjet / 1e3), on: 'Piles in the jet path', n: r(2 * inp.calls) }, { load: 'Jet on a wall', value: r(Fwall / 1e3), on: 'Quay face or fender panel', n: r(2 * inp.calls) }],
        explain: 'Each ship-related load with where it acts and how often: enter the forces in the strength and mooring modules and the counts in the fatigue module.' }],
      notes: [`Blockage ${f(Am / Ac * 100)} %; passing at ${f(inp.Vp / sch.Vlim * 100)} % of the limiting speed${sch.limited ? ' (capped at the limit: the speed entered cannot be reached in this section)' : ''}. Drawdown at the berth ${f(zB)} m over a ramp of ${f(Lramp)} m.`,
        `Wash steepness ${f(steep)}; jet outlet ${f(jet.D0)} m in diameter with ${f(jet.thrust / 1e3)} kN of thrust.`],
      advice: [
        ok ? eng('ok', `Surge force ${f(Fsurge / 1e3)} kN against ${f(inp.surgeCap)} kN of spring capacity; the largest ship-related pile load is ${f(FpileMax / 1e3)} kN against ${f(inp.pileCap)} kN.`,
          `Ships may pass at up to ${f(Vperm)} m/s (${f(Vperm / 0.5144)} knots) with this mooring. Confirm with a passing-ship model if large vessels pass closer than about two beams, and include the ${f(inp.passes + 2 * inp.calls)} events a year in the fatigue check.`, 'OCIMF MEG4 (2018) §3; PIANC WG 121 (2014); PIANC WG180 (2015)',
          `A speed limit is the cheapest mitigation: slowing by 20 % cuts the surge force by about ${f((1 - surgeOf(0.8 * inp.Vp) / Math.max(Fsurge, 1)) * 100)} %, and because propulsion power goes with roughly the cube of speed the passing ship burns about ${f((1 - 0.8 ** 3) * 100)} % less fuel per unit time along the berth frontage (${CO2_MGO} t CO₂ per tonne of marine gas oil saved).`)
          : eng('fail', `${w.name}: ${f(w.demand)} kN against ${f(w.capacity)} kN.`,
            w.name.startsWith('Passing') ? `Limit passing speed to ${f(Vperm)} m/s (${f(Vperm / 0.5144)} knots), add spring lines or shore-tension units, or increase the passing distance.` : `Limit the power used at this stand-off to ${f(Pperm / 1e3)} kW, move the thruster position away from the piles (load falls with the square of the distance), or protect the piles with deflector panels.`,
            'OCIMF MEG4 (2018) §3; PIANC WG180 (2015) §8–9',
            `An operating limit needs no material. Extra mooring hardware or pile strengthening would add steel at about 1.9 t CO₂e per tonne (screening factor), and a parted line is a safety hazard as well as lost berth time.`),
        eco('info', `About ${f(inp.passes)} passing events and ${f(2 * inp.calls)} manoeuvres a year act on this berth; the permissible passing speed is ${f(Vperm / 0.5144)} knots and the permissible jet power ${f(Pperm / 1e3)} kW.`,
          'Weigh a harbour speed limit (longer transits for every passing ship) against stronger moorings at this berth (a one-off cost and slower line handling per call). Line-tension monitoring settles the question with data at low cost.', `${ECON_BASIS}; port operations planning (PIANC WG 121)`),
      ],
      publish: { 'wharf.passingSurgeForce': Fsurge, 'wharf.washForcePile': pw.Fmax, 'wharf.jetForcePile': Fjet },
    };
  },
  verify: [
    // Schijf: Ac = 100·10 = 1000 m², Am = 1·10·10 = 100 m², V = 3 m/s. Fixed-point iteration by hand on
    //   v = V·Ac/(Ac − Am − bw·z), z = (v² − V²)/(2g):  z = 0.1 → 0.12045 → 0.12310 → 0.12346 → 0.12351 m; v = 3.3797 → Ur = 0.3797 m/s.
    { name: 'Return current and drawdown of a ship in a confined waterway', source: 'Schijf (1949) continuity and Bernoulli, solved by hand iteration (written out in the module source)',
      inputs: { h: 10, bw: 100, Vp: 3, Bp: 10, Tp: 10 * 0.9999, cm: 1 / 0.9999, kBerth: 1 }, expect: [{ kpi: 'z', value: 0.1235, tol: 0.003 }, { kpi: 'Ur', value: 0.3797, tol: 0.003 }] },
    // Wash: V = 6 m/s, h = 10 m, s = 100 m, α = 1: Fr = 6/√98.07 = 0.6059, Fr⁴ = 0.1348, (s/h)^(−1/3) = 0.4642 → H = 10·0.4642·0.1348 = 0.6255 m.
    { name: 'Wash wave height H = α h (s/h)^(−1/3) Fr⁴', source: 'Delft wash relation as given in PIANC and the Rock Manual (hand calculation)',
      inputs: { h: 10, Tp: 8, Vw: 6, sWash: 100, alphaW: 1 }, expect: [{ kpi: 'Hwash', value: 0.6255, tol: 0.003 }] },
    // Jet: P = 5 MW, Dp = 5 m open: V0 = 6.670 m/s, D0 = 3.536 m; at x = 20 m V = 2.8·6.670·3.536/20 = 3.302 m/s;
    // drag on a 1 m pile, Cd = 1: ½·1025·3.302²·20·√(π/30.8) = 5588·6.387 = 35.69 kN; wall force ρ(πD0²/4)V0² = 1025·9.817·44.49 = 447.7 kN.
    { name: 'Jet axis velocity, drag on a pile and momentum force on a wall', source: 'Blaauw & van de Kaa (1978) jet decay; momentum flux of a round jet (hand calculation)',
      inputs: { P: 5000, Dp: 5, ducted: false, cV0: 1.15, xJet: 20, D: 1, Cd: 1, h: 20, Tp: 12 }, expect: [{ kpi: 'Vjet', value: 3.302, tol: 0.003 }, { kpi: 'Fjet', value: 35.69, tol: 0.005 }, { kpi: 'Fwall', value: 447.7, tol: 0.005 }] },
  ],
};
