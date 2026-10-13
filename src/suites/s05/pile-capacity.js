// Pile and foundation capacity of a wharf: axial compression and tension capacity of the driven tubular piles,
// group efficiency, lateral capacity with the free length above the seabed, the penetration needed for the design
// loads, the permissible pile loads, and the bearing resistance of a shallow footing (landside rail beam or pad).
import { axialCapacity, capacityVsPenetration, groupBlockCapacityClay } from '../../engine/geotech/axialcapacity.js';
import { converseLabarre } from '../../engine/geotech/pilegroup.js';
import { solveLateralPile } from '../../engine/geotech/pile.js';
import { propertiesAt } from '../../engine/geotech/soilprofile.js';
import { bearingCapacity, terzaghiBearing } from '../../engine/geotech/bearing.js';
import { elasticSettlement } from '../../engine/geotech/settlement.js';
import { DEG, need, run, r, f, check, eng, eco, pileInputs, pileProps, soilInputs, profileFrom, pileSoilModelInputs, latOpts, lateralLimit, steelCo2, CONFIRM, ECON_BASIS, NOT_HIFI } from './_common.js';

export default {
  id: 's05-pile-capacity', suite: 5, title: 'Pile and foundation capacity',
  summary: 'Axial compression, tension and lateral capacity of the wharf piles with group effects, the penetration required for the design loads, the permissible pile loads, and the bearing resistance of a shallow landside footing.',
  covers: [5], fidelity: 'numerical',
  standards: ['API RP 2GEO (2011, add. 2014) §8.1–8.5', 'ISO 19901-4:2016 §8', 'EN 1997-1:2004 §7 (pile foundations), Annex D (bearing resistance)', 'BS 6349-2:2019 §7', 'Terzaghi (1943); Meyerhof (1963) bearing capacity'],
  models: ['Axial pile shaft resistance (API α and β methods)', 'Pile toe resistance, plugged and unplugged', 'Pile-group interaction (Converse–Labarre efficiency, block failure)', 'Lateral pile resistance curves (p–y) on a Winkler foundation', 'Terzaghi bearing capacity', 'Meyerhof bearing capacity with shape, depth and inclination factors', 'Elastic settlement of a rigid footing'],
  limits: `Static capacity of driven open-ended steel pipe piles in siliceous clay and sand from the API main-text methods; pile weight, drivability, cyclic degradation beyond the factor entered and negative skin friction are not modelled. The lateral limit is first hinge of an elastic pile on p–y springs. The footing check assumes homogeneous soil. ${NOT_HIFI}`,
  publishes: ['wharf.pileQc', 'wharf.pileQt', 'wharf.pileHult', 'wharf.pileLemb'],
  searchTerms: 'driven pile axial capacity API wharf pile group lateral capacity bearing capacity',
  inputs: [
    ...pileInputs(), ...soilInputs(), ...pileSoilModelInputs(),
    { key: 'plugMode', label: 'End bearing', type: 'select', default: 'auto', group: 'Axial capacity', role: 'bc', options: [{ value: 'auto', label: 'Smaller of plugged and unplugged' }, { value: 'plugged', label: 'Plugged' }, { value: 'unplugged', label: 'Unplugged (coring)' }] },
    { key: 'tensionFactor', label: 'Tension / compression shaft friction ratio', unit: '-', type: 'number', default: 0.8, min: 0.4, max: 1, group: 'Axial capacity', role: 'criteria', help: `Reduction of shaft friction in tension; 0.7–0.9 in sand, 1.0 in clay by the API main text. ${CONFIRM}` },
    { key: 'cyclicFactor', label: 'Set-up / cyclic factor on shaft friction', unit: '-', type: 'number', default: 1, min: 0.5, max: 1.5, group: 'Axial capacity', role: 'criteria', help: 'Below 1 for cyclic degradation under repeated berthing; above 1 only with re-strike evidence of set-up.' },
    { key: 'Nc', label: 'Design compression on the pile', unit: 'kN', type: 'number', default: 4500, min: 0, max: 60000, group: 'Design loads', role: 'bc', bind: 'wharf.pileAxialMax', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Largest factored pile compression from the strength module.' },
    { key: 'Nt', label: 'Design tension on the pile', unit: 'kN', type: 'number', default: 400, min: 0, max: 60000, group: 'Design loads', role: 'bc', bind: 'wharf.pileTensionMax' },
    { key: 'Hd', label: 'Design shear at the pile head', unit: 'kN', type: 'number', default: 400, min: 0, max: 20000, group: 'Design loads', role: 'bc', bind: 'wharf.pileShearMax' },
    { key: 'gRc', label: 'Resistance factor, compression', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Design loads', role: 'criteria', help: `Divides the ultimate capacity; used with factored loads. 1.25–1.5 in limit-state codes with load testing, 2.0 for working loads (API RP 2A-WSD operating). ${CONFIRM}` },
    { key: 'gRt', label: 'Resistance factor, tension', unit: '-', type: 'number', default: 1.75, min: 1, max: 3, group: 'Design loads', role: 'criteria', help: CONFIRM },
    { key: 'gRh', label: 'Resistance factor, lateral', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Design loads', role: 'criteria', help: CONFIRM },
    { key: 'rows', label: 'Pile rows across the wharf', unit: '-', type: 'int', default: 4, min: 1, max: 12, group: 'Pile group', role: 'geometry' },
    { key: 'cols', label: 'Bents acting together', unit: '-', type: 'int', default: 5, min: 1, max: 40, group: 'Pile group', role: 'geometry', help: 'Bents within one deck segment between expansion joints (for the block-failure check).' },
    { key: 'sGroup', label: 'Pile spacing (smaller of the two directions)', unit: 'm', type: 'number', default: 6, min: 1, max: 15, group: 'Pile group', role: 'geometry' },
    { key: 'footB', label: 'Footing width', unit: 'm', type: 'number', default: 2.5, min: 0.3, max: 20, group: 'Shallow footing', role: 'geometry', design: true, help: 'Landside crane rail beam on grade, a relieving-platform base or an equipment pad.' },
    { key: 'footL', label: 'Footing length (0 = strip)', unit: 'm', type: 'number', default: 0, min: 0, max: 100, group: 'Shallow footing', role: 'geometry' },
    { key: 'footD', label: 'Founding depth', unit: 'm', type: 'number', default: 1.2, min: 0, max: 10, group: 'Shallow footing', role: 'geometry' },
    { key: 'footPhi', label: 'Friction angle of the founding soil', unit: 'deg', type: 'number', default: 34, min: 0, max: 45, group: 'Shallow footing', role: 'soil', bind: 'site.phiDeg', uncertain: { dist: 'normal', cov: 0.07 } },
    { key: 'footC', label: 'Cohesion / undrained strength', unit: 'kPa', type: 'number', default: 0, min: 0, max: 500, group: 'Shallow footing', role: 'soil' },
    { key: 'footGamma', label: 'Effective unit weight of the founding soil', unit: 'kN/m³', type: 'number', default: 18, min: 5, max: 23, group: 'Shallow footing', role: 'soil', help: 'Bulk weight above the water table, submerged weight below.' },
    { key: 'footE', label: 'Drained modulus of the founding soil', unit: 'MPa', type: 'number', default: 40, min: 2, max: 500, group: 'Shallow footing', role: 'soil', help: 'For the elastic settlement of the footing: about 10–30 MPa loose sand or firm clay, 40–80 MPa dense sand or compacted fill.' },
    { key: 'footV', label: 'Design vertical load on the footing', unit: 'kN/m', type: 'number', default: 700, min: 0, max: 20000, group: 'Shallow footing', role: 'bc', help: 'Per metre run for a strip (for example landside rail load spread by the beam); total load in kN for a pad when a length is given.' },
    { key: 'footH', label: 'Design horizontal load on the footing', unit: 'kN/m', type: 'number', default: 40, min: 0, max: 5000, group: 'Shallow footing', role: 'bc' },
    { key: 'gRv', label: 'Resistance factor, bearing', unit: '-', type: 'number', default: 1.4, min: 1, max: 4, group: 'Shallow footing', role: 'criteria', help: `EN 1997-1 set R2 uses 1.4 with factored loads; 2.5–3 with working loads. ${CONFIRM}` },
  ],
  // Levels put nodes on the layer interfaces (6 m, 14 m) and the scour depth (1 m) of the default 28 m penetration.
  mesh: { key: 'nElem', label: 'Pile elements (embedded length)', levels: [14, 28, 56, 112], qoi: 'Hult' },
  compute(inp, ctx) {
    const P = pileProps(inp, ctx), profile = profileFrom(inp, P.D);
    const axO = { profile, D: P.D, t: P.t, plugMode: inp.plugMode, tensionFactor: inp.tensionFactor, setupFactor: inp.cyclicFactor };
    const ax = run(() => axialCapacity({ ...axO, L: inp.Lemb }), 'Axial capacity');
    const Lmax = Math.min(profile.depth - 0.5, Math.max(inp.Lemb * 1.6, inp.Lemb + 10));
    const cv = run(() => capacityVsPenetration({ ...axO, Lmin: Math.max(3, Math.min(5, inp.Lemb / 2)), Lmax, step: Math.max(0.5, (Lmax - 3) / 40) }), 'Capacity against penetration');
    // group effects: efficiency formula for friction piles in clay, and block failure of one deck segment
    const rows = Math.round(inp.rows), cols = Math.round(inp.cols);
    const toe = propertiesAt(profile, Math.min(inp.Lemb, profile.depth - 0.1)), mid = propertiesAt(profile, Math.min(inp.Lemb / 2, profile.depth - 0.1));
    const clayShaft = mid.type === 'clay';
    const Eg = clayShaft && rows * cols > 1 ? converseLabarre({ m: rows, n: cols, D: P.D, s: inp.sGroup }) : 1;
    let block = null;
    if (clayShaft && toe.type === 'clay' && rows * cols > 1) block = groupBlockCapacityClay({ B: (rows - 1) * inp.sGroup + P.D, W: (cols - 1) * inp.sGroup + P.D, L: inp.Lemb, suShaft: mid.su, suBase: toe.su, nPiles: rows * cols, Qsingle: ax.Qc });
    const eff = Math.min(Eg, block ? block.efficiency : 1);
    const QcG = ax.Qc * eff, QtG = ax.Qt * eff;
    const fixed = inp.headFix === 'fixed';
    const lat = lateralLimit(solveLateralPile, latOpts(inp, P, profile), P, fixed);
    need(lat.Hult > 0, 'The pile has no lateral capacity with this soil and free length: the first load step already fails. Increase the penetration or the diameter.');
    const RdC = QcG / inp.gRc, RdT = QtG / inp.gRt, RdH = lat.Hult / inp.gRh;
    let Lreq = NaN;
    for (let i = 0; i < cv.L.length; i++) if (cv.Qc[i] * eff / inp.gRc >= inp.Nc * 1e3 && cv.Qt[i] * eff / inp.gRt >= inp.Nt * 1e3) { Lreq = cv.L[i]; break; }
    // shallow footing: Meyerhof with inclination, Terzaghi for reference
    const strip = !(inp.footL > 0), phi = inp.footPhi * DEG, c = inp.footC * 1e3, gam = inp.footGamma * 1e3;
    need(phi > 0 || c > 0, 'The founding soil of the footing needs a friction angle or a cohesion.');
    const V = inp.footV * 1e3, Hh = inp.footH * 1e3;
    const bc = run(() => bearingCapacity({ method: 'meyerhof', B: inp.footB, L: strip ? Infinity : inp.footL, D: inp.footD, c, phi, gamma: gam, V: Math.max(V, 1), HB: Hh }), 'Bearing capacity');
    const tz = terzaghiBearing({ shape: strip ? 'strip' : 'square', B: inp.footB, c, phi, gamma: gam, q: gam * inp.footD });
    const Qfoot = bc.qult * inp.footB * (strip ? 1 : inp.footL), RdV = Qfoot / inp.gRv;
    // elastic settlement of the footing under the load entered: rigid base on a layer three widths thick (the depth that matters);
    // a strip is taken as ten widths long
    const Bs = strip ? inp.footB : Math.min(inp.footB, inp.footL), Ls = strip ? 10 * inp.footB : Math.max(inp.footB, inp.footL);
    const sFoot = V > 0 ? run(() => elasticSettlement({ q: V / (inp.footB * (strip ? 1 : inp.footL)), B: Bs, L: Ls, E: inp.footE * 1e6, nu: 0.3, point: 'rigid', H: 3 * Bs }), 'Footing settlement').s : 0;
    const checks = [
      check('Pile axial compression', inp.Nc, RdC / 1e3, 'kN', 'API RP 2GEO (2011) §8.1; EN 1997-1:2004 §7.6.2', 'Factored pile compression against shaft friction plus end bearing over the resistance factor, with the group efficiency. If exceeded the pile plunges and the deck settles differentially.'),
      check('Pile axial tension (pull-out)', inp.Nt, RdT / 1e3, 'kN', 'API RP 2GEO (2011) §8.1; EN 1997-1:2004 §7.6.3', 'Factored pile tension against shaft friction alone. If exceeded the pile ratchets upward under repeated berthing or mooring loads.'),
      check('Pile lateral load', inp.Hd, RdH / 1e3, 'kN', 'API RP 2GEO (2011) §8.5; EN 1997-1:2004 §7.7', 'Pile-head shear against the load that forms the first plastic hinge or exhausts the soil. If exceeded the pile is permanently bent and the wharf keeps a residual sway.'),
      check('Footing bearing resistance', inp.footV, RdV / 1e3, strip ? 'kN/m' : 'kN', 'EN 1997-1:2004 Annex D; Meyerhof (1963)', 'Vertical load on the shallow footing against the bearing resistance of the soil wedge below it, reduced for load inclination. If exceeded the footing punches or rotates and the rail loses level.'),
    ];
    const w = checks.reduce((a, q) => (q.uc > a.uc ? q : a)), ok = w.uc <= 1;
    const dL = Number.isFinite(Lreq) ? Lreq - inp.Lemb : NaN, tPerM = P.massPerM / 1e3;
    const zz = ax.z.map((z) => r(z)), stride = Math.max(1, Math.floor(zz.length / 150)), pk = (a) => a.filter((_, i) => i % stride === 0);
    return {
      kpis: [
        { key: 'Qs', label: 'Shaft friction (compression)', value: r(ax.Qfo / 1e3), unit: 'kN', explain: 'Outer shaft friction integrated over the penetration with the α method in clay and the β method in sand, for a single pile.' },
        { key: 'Qp', label: 'End bearing', value: r(ax.Qp / 1e3), unit: 'kN', explain: `Toe resistance; the pile behaves ${ax.plugged ? 'plugged (full base area)' : 'unplugged (steel annulus plus inner friction)'} because that mode gives the smaller value.` },
        { key: 'Qc', label: 'Ultimate compression capacity, single pile', value: r(ax.Qc / 1e3), unit: 'kN', explain: 'Shaft friction plus end bearing. The design resistance divides this by the resistance factor and multiplies by the group efficiency.' },
        { key: 'Qt', label: 'Ultimate tension capacity, single pile', value: r(ax.Qt / 1e3), unit: 'kN', explain: 'Shaft friction times the tension ratio; no end bearing acts in uplift.' },
        { key: 'eff', label: 'Group efficiency', value: r(eff), unit: '-', explain: 'Capacity of a pile in the group over that of an isolated pile: the smaller of the Converse–Labarre factor and block failure of the deck segment. 1.0 in sand and at wide spacing.' },
        { key: 'RdC', label: 'Maximum permissible factored compression per pile', value: r(RdC / 1e3), unit: 'kN', limit: inp.Nc, sense: 'min', explain: 'Design compression resistance of a pile in the group. It must not be less than the factored pile load; the difference is the reserve for heavier cargo or cranes.' },
        { key: 'RdT', label: 'Maximum permissible factored tension per pile', value: r(RdT / 1e3), unit: 'kN', limit: inp.Nt, sense: 'min', explain: 'Design uplift resistance. It limits the berthing and mooring force the bent can take once pile tension develops.' },
        { key: 'Hult', label: 'Ultimate lateral capacity, single pile', value: r(lat.Hult / 1e3, 6), unit: 'kN', explain: `Pile-head shear at ${lat.mode}, with the free length above the seabed and the scour entered.` },
        { key: 'Mp', label: 'Plastic moment of the pile', value: r(P.Mp / 1e3), unit: 'kN·m', explain: 'Full plastic bending resistance of the corroded tube, fy·(D³ − Di³)/6. The lateral limit is reached when the largest moment equals it.' },
        { key: 'Lreq', label: 'Penetration required for the axial loads', value: r(Number.isFinite(Lreq) ? Lreq : Lmax), unit: 'm', explain: `Shallowest penetration at which both the compression and tension design resistances cover the loads${Number.isFinite(Lreq) ? '' : ' (not reached within the depth examined: the value shown is the search limit)'}.` },
        { key: 'qult', label: 'Footing ultimate bearing pressure (Meyerhof)', value: r(bc.qult / 1e3), unit: 'kPa', explain: 'Gross ultimate pressure under the footing with shape, depth and load-inclination factors.' },
        { key: 'qultT', label: 'Footing ultimate bearing pressure (Terzaghi)', value: r(tz.qult / 1e3), unit: 'kPa', explain: 'Classical value for a vertical central load without depth or inclination factors: a cross-check on the Meyerhof result.' },
        { key: 'sFoot', label: 'Elastic settlement of the footing', value: r(sFoot * 1e3), unit: 'mm', explain: 'Immediate settlement of the footing under the load entered, for a rigid base on an elastic layer three footing widths thick. It sets the level difference between a rail on grade and a rail on the piled deck.' },
        { key: 'RdV', label: 'Maximum permissible footing load', value: r(RdV / 1e3), unit: strip ? 'kN/m' : 'kN', explain: 'Design bearing resistance of the footing; the landside rail load per metre (or pad load) must stay below it.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Capacity against penetration', x: { label: 'Capacity', unit: 'kN' }, y: { label: 'Penetration below the seabed', unit: 'm' },
          series: [{ name: 'Compression', x: cv.Qc.map((v) => r(v * eff / 1e3)), y: cv.L.map((v) => r(-v)) }, { name: 'Tension', x: cv.Qt.map((v) => r(v * eff / 1e3)), y: cv.L.map((v) => r(-v)) },
            { name: 'Required compression (γR·N)', x: [inp.Nc * inp.gRc, inp.Nc * inp.gRc], y: [r(-cv.L[0]), r(-Lmax)], dash: true }, { name: 'Present toe level', x: [0, r(Math.max(...cv.Qc) * eff / 1e3)], y: [-inp.Lemb, -inp.Lemb], dash: true }],
          explain: 'Ultimate capacity of a pile in the group as the toe goes deeper. The pile is long enough where the compression curve passes the dashed load line; steps mark stronger layers worth founding in.' },
        { type: 'line', title: 'Unit shaft friction with depth', x: { label: 'Unit friction', unit: 'kPa' }, y: { label: 'Depth below the seabed', unit: 'm' }, series: [{ name: 'f', x: pk(ax.f.map((v) => r(v / 1e3))), y: pk(zz.map((z) => -z)) }],
          explain: 'Friction available on the pile wall at each depth. Zero in the scoured zone; jumps show the layer boundaries of the borehole log.' },
        { type: 'line', title: 'Lateral load against head deflection', x: { label: 'Head deflection', unit: 'mm' }, y: { label: 'Pile-head shear', unit: 'kN' },
          series: [{ name: fixed ? 'Fixed head' : 'Free head', x: lat.y.map((v) => r(v * 1e3)), y: lat.H.map((v) => r(v / 1e3)) }, { name: 'Design shear', x: [0, r(lat.y[lat.y.length - 1] * 1e3)], y: [inp.Hd, inp.Hd], dash: true }],
          explain: 'Nonlinear response of the pile in soil up to first hinge. The curve softens as the upper soil yields; the design shear should sit on the stiff initial part.' },
      ],
      tables: [{ title: 'Capacity summary', columns: [{ key: 'item', label: 'Resistance' }, { key: 'ult', label: 'Ultimate', unit: 'kN' }, { key: 'factor', label: 'γR' }, { key: 'design', label: 'Design', unit: 'kN' }, { key: 'load', label: 'Load', unit: 'kN' }, { key: 'uc', label: 'Utilisation' }],
        rows: [{ item: 'Compression (group pile)', ult: r(QcG / 1e3), factor: inp.gRc, design: r(RdC / 1e3), load: inp.Nc, uc: r(checks[0].uc, 3) }, { item: 'Tension (group pile)', ult: r(QtG / 1e3), factor: inp.gRt, design: r(RdT / 1e3), load: inp.Nt, uc: r(checks[1].uc, 3) },
          { item: 'Lateral (single pile)', ult: r(lat.Hult / 1e3), factor: inp.gRh, design: r(RdH / 1e3), load: inp.Hd, uc: r(checks[2].uc, 3) }, { item: `Footing bearing${strip ? ' (per metre)' : ''}`, ult: r(Qfoot / 1e3), factor: inp.gRv, design: r(RdV / 1e3), load: inp.footV, uc: r(checks[3].uc, 3) }],
        explain: 'Each foundation resistance with its factor and the load it is compared with. The design column is the largest factored load permitted.' }],
      notes: [`Shaft friction supplies ${f(ax.Qfo / ax.Qc * 100)} % of the compression capacity; the pile is ${ax.plugged ? 'plugged' : 'unplugged'} (unit end bearing ${f(ax.qTip / 1e6)} MPa).`,
        `Group efficiency ${f(eff)}${block ? ` (Converse–Labarre ${f(Eg)}, block failure ${f(block.efficiency)})` : clayShaft ? ` (Converse–Labarre)` : ' (sand at mid-depth: no reduction)'}; lateral limit by ${lat.mode}.`],
      advice: [
        eng(ok ? 'ok' : 'fail', `${w.name} governs at ${f(w.uc)}. Design resistances per pile: compression ${f(RdC / 1e3)} kN, tension ${f(RdT / 1e3)} kN, lateral ${f(RdH / 1e3)} kN; footing ${f(RdV / 1e3)} ${strip ? 'kN/m' : 'kN'}.`,
          ok ? (Number.isFinite(dL) && dL < -1 ? `The axial loads would be carried at ${f(Lreq)} m penetration. Before shortening, confirm lateral fixity, drivability and the cyclic factor.` : 'Confirm the shaft friction with CPT-based methods and dynamic testing with signal matching on the first piles driven.')
            : w.name.startsWith('Pile axial') ? (Number.isFinite(Lreq) ? `Drive to about ${f(Lreq)} m penetration (now ${f(inp.Lemb)} m).` : `The axial loads cannot be carried within ${f(Lmax)} m: use a larger diameter or more piles per bent.`)
              : w.name.startsWith('Pile lateral') ? `Increase the wall or diameter (plastic moment needed about ${f(P.Mp * w.uc / 1e3)} kN·m against ${f(P.Mp / 1e3)} kN·m) or add raker piles.` : `Widen the footing to about ${f(inp.footB * Math.sqrt(w.uc))} m or found it deeper.`,
          'API RP 2GEO (2011, add. 2014) §8; ISO 19901-4:2016; EN 1997-1:2004 §7',
          `Each metre of pile is ${f(tPerM)} t of steel${Number.isFinite(dL) ? `; the required penetration differs from the present one by ${f(dL)} m per pile, ${steelCo2(dL * tPerM)}` : ` (${steelCo2(tPerM)})`}.`),
        eco(ok ? 'info' : 'warn', `The permissible factored pile loads are ${f(RdC / 1e3)} kN compression and ${f(RdT / 1e3)} kN tension against ${f(inp.Nc)} and ${f(inp.Nt)} kN applied; penetration costs ${f(tPerM)} t of steel per metre per pile plus driving time.`,
          'A static or dynamic pile load test on the first bent usually justifies a lower resistance factor than the default and shortens every pile that follows: compare the test cost with the steel and hammer time saved over the whole wharf.', `${ECON_BASIS}; value-of-information for pile testing`),
      ],
      publish: { 'wharf.pileQc': QcG, 'wharf.pileQt': QtG, 'wharf.pileHult': lat.Hult, 'wharf.pileLemb': inp.Lemb },
      quantities: { steel_t: r(tPerM * (inp.Lemb + inp.Lfree)) },
    };
  },
  verify: [
    // Uniform clay su = 20 kPa, γ′ = 10 kN/m³, D = 1 m, L = 20 m, α = 0.5ψ^−0.5 (ψ ≤ 1), 0.5ψ^−0.25 (ψ > 1), α ≤ 1, ψ = su/σ′v:
    //  0–2 m (ψ > 1): f = 10·(z/2)^0.25 kPa → ∫ = 10·2/1.25 = 16.00 kPa·m
    //  2–8 m: f = 0.5·√(20·10z) = 7.0711·√z → ∫ = 7.0711·(2/3)(8^1.5 − 2^1.5) = 93.33
    //  8–20 m (α = 1): 20·12 = 240     Σ = 349.33 kPa·m × πD = 1097.5 kN;  plugged tip: 9·20·π/4 = 141.37 kN.
    { name: 'API α method in uniform clay: shaft friction and plugged end bearing', source: 'Hand integration of the API RP 2GEO §8.1 α expressions',
      inputs: { D: 1, t: 25, corr: 0, Lemb: 20, Lfree: 10, scour: 0, plugMode: 'plugged', tensionFactor: 1, cyclicFactor: 1, rows: 1, cols: 1, layers: [{ type: 'clay', thick: 40, gamma: 10, su: 20, suGrad: 0, phi: 0, Dr: 0, eps50: 0.01 }] },
      expect: [{ kpi: 'Qs', value: 1097.5, tol: 0.01 }, { kpi: 'Qp', value: 141.37, tol: 1e-3 }, { kpi: 'Qt', value: 1097.5, tol: 0.01 }] },
    // Terzaghi strip footing on the surface of a purely cohesive soil: qult = c·Nc = 100·5.7 = 570 kPa (Nc = 1.5π + 1 = 5.712).
    { name: 'Terzaghi strip footing on clay: qult = 5.71·su', source: 'Terzaghi (1943) Theoretical Soil Mechanics, Nc = 1.5π + 1 for φ = 0',
      inputs: { footPhi: 0, footC: 100, footD: 0, footB: 2, footL: 0, footH: 0 }, expect: [{ kpi: 'qultT', value: 571.2, tol: 0.003 }] },
    // Meyerhof / Prandtl strip on clay at the surface, vertical load: qult = (2 + π)·su = 514.2 kPa.
    { name: 'Prandtl strip footing on clay: qult = (2 + π)·su', source: 'Prandtl (1921) exact plasticity solution, as used by Meyerhof (1963)',
      inputs: { footPhi: 0, footC: 100, footD: 0, footB: 2, footL: 0, footH: 0 }, expect: [{ kpi: 'qult', value: 514.2, tol: 0.003 }] },
  ],
};
