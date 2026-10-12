// Effective tension and bending moment of vertical tensioned risers: top-tensioned (production / completion)
// risers on a tensioner and free-standing hybrid riser towers held up by a buoyancy can.
import { solveTtr, solveFreeStandingRiser, stressAlong } from '../../engine/lines/riser.js';
import {
  G, need, run, fmt, round, thin, scaled, pipeInputs, waterInputs, hydroInputs, currentInputs, currentFn, pipeFrom, steelTonnes, adv, check, kpi, plot,
  meshInput, finish, handWeight, stations, steelCo2Text, maxAbs, DEG, CONFIRM,
} from './_common.js';
import { criteriaInputs } from './_scr.js';

const ttr = { key: 'riserType', in: ['ttr'] }, hyb = { key: 'riserType', in: ['hybrid'] };
const W_TTR = handWeight({ D: 339.7, t: 15.9, coatT: 0, rhoInt: 800 });

export default {
  id: 's09-tension-moment', suite: 9, title: 'Effective tension and bending moment',
  summary: 'Effective and true wall tension, bending moment, shear and joint angles along a top-tensioned riser or a free-standing hybrid riser tower under offset and current.',
  covers: [3], fidelity: 'numerical',
  standards: ['API RP 2RD (1998, R2013) §5', 'API RP 16Q (2017) §3 riser analysis', 'DNV-ST-F201 (2018) Sec. 4–5'],
  models: ['Tensioned Euler–Bernoulli beam with geometric stiffness', 'Effective tension formulation (Sparks)', 'Morison current drag', 'Tensioner stiffness', 'Catenary jumper (hybrid riser)', 'Lamé thick-wall stresses'],
  limits: 'In-plane (2-D) quasi-static analysis: out-of-plane current, torsion and tensioner dynamics are not modelled. Concentric casing/tubing strings are represented as one equivalent pipe; the hybrid buoyancy can is a point uplift and the jumper an inextensible catenary without bending stiffness.',
  publishes: ['riser.ttrTopTension', 'riser.ttrStroke', 'riser.ttrLowerAngleDeg'],
  searchTerms: 'top tensioned riser effective tension tensioner hybrid riser tower',
  inputs: [
    { key: 'riserType', label: 'Riser type', type: 'select', default: 'ttr', group: 'Configuration', role: 'geometry',
      options: [{ value: 'ttr', label: 'Top-tensioned riser (tensioner)' }, { value: 'hybrid', label: 'Free-standing hybrid riser (buoyancy can + jumper)' }] },
    ...pipeInputs({ D: 339.7, t: 15.9, rhoInt: 800, coat: 0 }),
    ...waterInputs(),
    ...hydroInputs(),
    { key: 'length', label: 'Riser length, wellhead to tension ring', unit: 'm', type: 'number', default: 980, min: 30, max: 3600, group: 'Configuration', role: 'geometry', showIf: ttr },
    { key: 'tensionFactor', label: 'Top tension factor (top tension / apparent weight)', unit: '-', type: 'number', default: 1.4, min: 1.02, max: 3, group: 'Top end', role: 'bc', showIf: ttr, design: true,
      help: 'Production TTRs typically run 1.3–1.6; below about 1.1 the base approaches compression.' },
    { key: 'tensionerK', label: 'Tensioner stiffness', unit: 'kN/m', type: 'number', default: 0, min: 0, max: 5000, group: 'Top end', role: 'bc', showIf: ttr,
      help: 'Tension gained per metre of stroke-out. 0 = ideal constant-tension system; hydro-pneumatic tensioners are typically 20–200 kN/m.' },
    { key: 'kUpper', label: 'Upper joint rotational stiffness', unit: 'kN·m/deg', type: 'number', default: 0, min: 0, max: 5000, group: 'Top end', role: 'bc', showIf: ttr, help: '0 = free rotation at the tension ring / keel joint.' },
    { key: 'towerLength', label: 'Tower length, seabed to buoyancy can', unit: 'm', type: 'number', default: 850, min: 30, max: 3400, group: 'Configuration', role: 'geometry', showIf: hyb },
    { key: 'canUplift', label: 'Net uplift of the buoyancy can', unit: 'kN', type: 'number', default: 2500, min: 50, max: 30000, group: 'Top end', role: 'bc', showIf: hyb, design: true,
      help: 'Can buoyancy minus can weight; it is the top tension of the tower.' },
    { key: 'canDragArea', label: 'Can drag area Cd·A', unit: 'm²', type: 'number', default: 60, min: 0, max: 600, group: 'Top end', role: 'load', showIf: hyb },
    { key: 'jumperLength', label: 'Flexible jumper length', unit: 'm', type: 'number', default: 450, min: 20, max: 3000, group: 'Jumper', role: 'geometry', showIf: hyb },
    { key: 'jumperW', label: 'Jumper apparent weight', unit: 'N/m', type: 'number', default: 800, min: 20, max: 5000, group: 'Jumper', role: 'load', showIf: hyb },
    { key: 'vesselX', label: 'Vessel hang-off horizontal distance from tower base', unit: 'm', type: 'number', default: 300, min: 5, max: 2500, group: 'Jumper', role: 'bc', showIf: hyb },
    { key: 'jumperHangDepth', label: 'Jumper hang-off depth below still water', unit: 'm', type: 'number', default: 20, min: 0, max: 200, group: 'Jumper', role: 'bc', showIf: hyb },
    { key: 'kLower', label: 'Lower joint rotational stiffness', unit: 'kN·m/deg', type: 'number', default: 50, min: 0, max: 50000, group: 'Bottom end', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.25 },
      help: '0 = free ball joint. Flex joints 20–100 kN·m/deg; a tapered stress joint behaves like several thousand kN·m/deg.' },
    { key: 'offset', label: 'Vessel offset', unit: 'm', type: 'number', default: 30, min: 0, max: 400, group: 'Top end', role: 'bc', showIf: ttr, bind: 'vessel.offsetMax', uncertain: { dist: 'normal', cov: 0.2 } },
    ...currentInputs({ Us: 1.0 }),
    { key: 'pInt', label: 'Internal pressure at the top', unit: 'bar', type: 'number', default: 200, min: 0, max: 1500, group: 'Internal fluid', role: 'bc' },
    ...criteriaInputs(),
    { key: 'minTF', label: 'Minimum tension factor at the base', unit: '-', type: 'number', default: 1.1, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria',
      help: `Top tension must exceed this multiple of the apparent weight so the base stays in effective tension with margin. ${CONFIRM}` },
    { key: 'angleAllow', label: 'Allowable lower joint angle', unit: 'deg', type: 'number', default: 5, min: 0.5, max: 15, group: 'Acceptance criteria', role: 'criteria', help: 'Equipment rating of the lower flex / stress joint for the load case considered.' },
    meshInput(100, '50 / 100 / 200'),
  ],
  mesh: { key: 'nElem', label: 'Beam elements', levels: [50, 100, 200, 400], qoi: 'maxMoment' },

  compute(inp, ctx) {
    // engineering units → SI: mm → m (pipeFrom), kN·m/deg → N·m/rad, kN → N, bar → Pa
    const pp = pipeFrom(inp, ctx);
    const sec = pp.section, isT = inp.riserType === 'ttr';
    const L = isT ? inp.length : inp.towerLength;
    need(sec.w > 0, `The riser string is buoyant (apparent weight ${Math.round(sec.w)} N/m); a tensioned riser needs positive apparent weight for this model.`);
    need(isT ? L <= inp.depth + 40 : L <= inp.depth - 30, isT ? `Riser length ${L} m exceeds the water depth plus a 40 m deck elevation.` : `Tower length ${L} m must end at least 30 m below the surface (depth ${inp.depth} m).`);
    const cf = currentFn(inp, inp.depth);
    const current = inp.Us > 0 || inp.Ub > 0 ? { rho: pp.rhoW, profile: (z) => cf(z) } : undefined;
    const kRad = (k) => (k > 0 ? k * 1e3 / DEG : undefined);
    const pressure = { rhoInt: inp.rhoInt, pIntRef: inp.pInt * 1e5, zIntRef: L, rhoExt: pp.rhoW, zSurface: inp.depth, nu: pp.mat.nu };
    const W = sec.w * L;
    let r, extra = {};
    if (isT) {
      need(inp.offset <= 0.2 * L, `An offset of ${inp.offset} m is more than 20 % of the riser length; a top-tensioned riser cannot follow it.`);
      r = run(() => solveTtr({ length: L, pipe: sec, tensionFactor: inp.tensionFactor, offset: inp.offset, current, kFlexLower: kRad(inp.kLower), kFlexUpper: kRad(inp.kUpper), tensionerStiffness: inp.tensionerK * 1e3, nElem: inp.nElem, pressure }), 'The top-tensioned riser did not solve');
    } else {
      need(inp.canUplift * 1e3 > W, `Can uplift ${inp.canUplift} kN is below the tower apparent weight ${fmt(W / 1e3)} kN: the tower would not stand. Increase the uplift.`);
      const vz = inp.depth - inp.jumperHangDepth;
      need(inp.jumperLength > Math.hypot(inp.vesselX, vz - L) * 1.03, `The jumper (${inp.jumperLength} m) is too short to span from the can to the vessel (${fmt(Math.hypot(inp.vesselX, vz - L))} m apart).`);
      const Utop = cf(L);
      const h = run(() => solveFreeStandingRiser({ length: L, pipe: sec, canNetUplift: inp.canUplift * 1e3, jumper: { length: inp.jumperLength, w: inp.jumperW, vesselX: inp.vesselX, vesselZ: vz }, kFlexBase: kRad(inp.kLower), current, canDragForce: 0.5 * pp.rhoW * inp.canDragArea * Utop * Utop, nElem: inp.nElem }), 'The hybrid riser did not solve');
      const nn = h.rod.s.length - 1;
      r = { topTension: h.topTension, bottomTension: h.baseTension, tensionFactor: h.topTension / W, lowerAngleDeg: h.baseAngleDeg, upperAngleDeg: (Math.PI / 2 - h.rod.angle[nn]) / DEG, stroke: h.topSetDown, maxMoment: h.maxMoment, maxCurvature: maxAbs(h.rod.curvature), maxDeflection: h.topOffset, lowerJointMoment: Math.abs(h.rod.moment[0]), rod: h.rod };
      extra = h;
    }
    const rod = r.rod;
    need(rod.converged !== false, 'The static solution did not converge; raise the tension or reduce offset and current.');
    const n = rod.s.length - 1;
    const st = stressAlong(rod, { D: pp.D, t: pp.t, pIntRef: pressure.pIntRef, zIntRef: L, rhoInt: inp.rhoInt, rhoExt: pp.rhoW, zSurface: inp.depth });
    // true wall tension from effective tension and the local pressures (Sparks 2007 Eq. 2.5)
    const Tw = Array.from(rod.Te, (te, k) => { const pi = Math.max(0, pressure.pIntRef + inp.rhoInt * G * (L - rod.z[k])), pe = Math.max(0, pp.rhoW * G * (inp.depth - rod.z[k])); return te + pi * pp.props.Ai - pe * pp.props.Ae; });
    const fy = pp.mat.fy, allow = inp.Cf * inp.CaAllow * fy;
    const minTe = Math.min(...rod.Te);
    const bend = r.maxMoment * (pp.D / 2) / pp.props.I;

    const kpis = [
      kpi('topTension', 'Top effective tension', r.topTension / 1e3, 'kN', 'Tension the tensioner or buoyancy can applies at the top. It is the payload the host or can must supply and the main control on riser bending.'),
      kpi('bottomTension', 'Effective tension at the base', r.bottomTension / 1e3, 'kN', 'Effective tension at the wellhead or tower base: top tension less the apparent weight below. It governs lateral stiffness and must stay positive.', 0, 'min'),
      kpi('TwTop', 'True wall tension at the top', Tw[n] / 1e3, 'kN', 'Axial force actually carried by the steel at the top (effective tension plus the internal and minus the external pressure end-cap terms). Use it for axial wall stress and connector loads.'),
      kpi('TwBottom', 'True wall tension at the base', Tw[0] / 1e3, 'kN', 'Axial force in the steel at the base. It differs from the effective tension by the pressure terms pi·Ai − pe·Ae and can be compressive even when the effective tension is positive.'),
      kpi('tensionFactor', 'Tension factor', r.tensionFactor, '-', 'Top tension over total apparent weight of the string. Higher values straighten the riser and reduce joint angles at the cost of tensioner capacity and hull payload.', inp.minTF, 'min'),
      kpi('lowerAngle', 'Lower joint angle', Math.abs(r.lowerAngleDeg), 'deg', 'Inclination of the riser from vertical at the wellhead / tower base. It drives wear, the bending moment into the wellhead and the stress-joint design.', inp.angleAllow, 'max'),
      kpi('upperAngle', 'Upper joint angle', Math.abs(r.upperAngleDeg), 'deg', 'Inclination of the riser at the top end relative to vertical; it loads the keel joint or the can connection.'),
      kpi('stroke', isT ? 'Tensioner stroke (set-down)' : 'Tower-top set-down', r.stroke, 'm', 'Vertical movement of the top end caused by offset and current. For a TTR it must stay inside the tensioner stroke range with margin for heave and tide.'),
      kpi('maxDeflection', isT ? 'Largest deflection from the straight chord' : 'Tower-top horizontal offset', r.maxDeflection, 'm', 'Lateral bow of the riser under current (TTR) or the excursion of the can (hybrid). It feeds the clearance check against neighbouring risers.'),
      kpi('maxMoment', 'Maximum bending moment', r.maxMoment / 1e3, 'kN·m', 'Largest bending moment along the string, normally at the lower joint when it has rotational stiffness, otherwise in the high-current zone.'),
      kpi('bendStress', 'Maximum bending stress', bend / 1e6, 'MPa', 'Outer-fibre bending stress M·(D/2)/I at the point of largest moment; it is the cyclic part of the stress when the vessel moves.'),
      kpi('vmMax', 'Maximum von Mises stress', st.max / 1e6, 'MPa', 'Largest combined stress from true wall tension, bending and pressure along the riser, compared with the API RP 2RD allowable Cf·Ca·σy.', allow / 1e6, 'max'),
    ];
    if (!isT) kpis.push(kpi('jumperTension', 'Jumper tension at the vessel', extra.jumperVesselTension / 1e3, 'kN', 'Tension of the flexible jumper at the vessel hang-off; it is the only riser load the host sees in a hybrid system.'),
      kpi('jumperH', 'Jumper horizontal pull on the tower', extra.jumperH / 1e3, 'kN', 'Horizontal force the jumper applies to the tower top; with current drag on the can it sets the tower inclination.'));

    const checks = [
      check('Combined stress (von Mises)', st.max / 1e6, allow / 1e6, 'MPa', 'API RP 2RD (1998) §5.2.3, Table 2', 'Guards against yielding of the pipe wall under true wall tension, bending and pressure. Beyond the allowable the code margin against plastic deformation is lost.'),
      check('Minimum top tension', inp.minTF * W / 1e3, r.topTension / 1e3, 'kN', 'API RP 16Q (2017) §3 minimum top tension; API RP 2RD §5', 'Keeps the whole string in positive effective tension with margin. If the base goes into effective compression the riser buckles sideways like a column.'),
      check('Lower joint angle', Math.abs(r.lowerAngleDeg), inp.angleAllow, 'deg', 'Equipment rating / API RP 16Q operating limits (project value)', 'Limits rotation of the lower flex or stress joint. Exceeding it overloads the joint and transfers excessive bending into the wellhead and conductor.'),
    ];

    const zs = Array.from(rod.z);
    const tt = (a, f) => { const [z, y] = thin([zs, a], 300); return { z, y: scaled(y, f) }; };
    const te = tt(rod.Te, 1e-3), tw = tt(Tw, 1e-3), mo = tt(rod.moment, 1e-3), sh = tt(rod.shear, 1e-3), vm = tt(st.vonMises, 1e-6), xx = tt(rod.x, 1);
    const [x3, z3, v3] = thin([rod.x, rod.z, st.vonMises], 150);
    const plots = [
      plot('Effective and true wall tension', 'Tension', 'kN', 'Height above base', 'm', [{ name: 'Effective tension', x: te.y, y: te.z }, { name: 'True wall tension', x: tw.y, y: tw.z }], 'Both tensions against elevation. Effective tension controls shape and buckling; true wall tension gives the axial stress. The gap between them is the pressure end-cap term, largest at depth.'),
      plot('Deflected shape', 'Horizontal position', 'm', 'Height above base', 'm', [{ name: 'Riser', x: xx.y, y: xx.z }, { name: 'Straight chord', x: [0, rod.x[n]], y: [0, rod.z[n]], dash: true }], 'Riser axis under offset and current. The bow away from the chord shows where current drag acts; the slope at the bottom is the lower joint angle.'),
      plot('Bending moment', 'Bending moment', 'kN·m', 'Height above base', 'm', [{ name: 'Moment', x: mo.y, y: mo.z }], 'Bending moment against elevation. A peak at the base means the lower joint stiffness governs; a mid-depth peak means current drag governs.'),
      plot('Shear force', 'Shear force', 'kN', 'Height above base', 'm', [{ name: 'Shear', x: sh.y, y: sh.z }], 'Transverse shear along the riser; the base value is the lateral load delivered to the wellhead or tower foundation.'),
      plot('von Mises stress', 'von Mises stress', 'MPa', 'Height above base', 'm', [{ name: 'von Mises', x: vm.y, y: vm.z }, { name: 'Allowable', x: [allow / 1e6, allow / 1e6], y: [0, rod.z[n]], dash: true }], 'Combined stress against the allowable (dashed). The governing elevation is where the curve is closest to the limit.'),
      { type: 'model3d', title: 'Riser model (coloured by von Mises stress)', nodes: z3.map((z) => [0, 0, z]), deformed: x3.map((x, i) => [x, 0, z3[i]]), lines: x3.slice(1).map((_, i) => [i, i + 1]), nodeValues: scaled(v3, 1e-6), valueLabel: 'von Mises stress', valueUnit: 'MPa',
        explain: 'Undeformed vertical riser and its deflected shape under offset and current, coloured by combined wall stress.' },
    ];
    const tables = [{ title: 'Section forces along the riser', columns: [{ key: 'z', label: 'Height', unit: 'm' }, { key: 'x', label: 'Deflection', unit: 'm' }, { key: 'Te', label: 'Effective tension', unit: 'kN' }, { key: 'Tw', label: 'True wall tension', unit: 'kN' }, { key: 'M', label: 'Moment', unit: 'kN·m' }, { key: 'V', label: 'Shear', unit: 'kN' }, { key: 'vm', label: 'von Mises', unit: 'MPa' }],
      rows: stations({ z: rod.z, x: rod.x, Te: scaled(rod.Te, 1e-3), Tw: scaled(Tw, 1e-3), M: scaled(rod.moment, 1e-3), V: scaled(rod.shear, 1e-3), vm: scaled(st.vonMises, 1e-6) }, 30),
      explain: 'Axial force (effective and true wall), bending moment and shear at stations from the base upwards; use the rows at joints and connectors for component checks.' }];

    const steel = steelTonnes(pp.D, pp.t, L, pp.mat.rho);
    const ucs = checks.map((c) => c.uc), ucMax = Math.max(...ucs), gov = checks[ucs.indexOf(ucMax)];
    const dT = 0.1 * W;
    const advice = [
      adv('engineering', ucMax > 1 ? 'fail' : ucMax > 0.9 ? 'warn' : 'ok',
        `Governing check: ${gov.name}, utilisation ${fmt(ucMax)}. Tension factor ${fmt(r.tensionFactor)}, lower joint angle ${fmt(Math.abs(r.lowerAngleDeg))}°, base effective tension ${fmt(r.bottomTension / 1e3)} kN.`,
        ucMax > 1 ? (gov.name.startsWith('Lower') ? `Raise the ${isT ? 'tension factor' : 'can uplift'} or reduce the allowed offset; the angle falls roughly in proportion to 1/tension.` : gov.name.startsWith('Minimum') ? `Increase top tension to at least ${fmt(inp.minTF * W / 1e3)} kN.` : `Reduce bending by softening the lower joint (currently ${inp.kLower} kN·m/deg) or thicken the wall to about ${fmt(pp.t * 1e3 * st.max / (0.9 * allow))} mm.`)
          : 'Confirm the tension setting against the dynamic stroke range and the fatigue of the lower joint before fixing it; static margins here do not cover vessel-motion cycling.',
        `API RP 2RD (1998, R2013) §5 and API RP 16Q (2017) §3. ${CONFIRM}`),
      adv('economic', 'info',
        `Top tension of ${fmt(r.topTension / 1e3)} kN (${fmt(r.topTension / 1e3 / G)} t) is permanent ${isT ? 'tensioner load on the hull' : 'buoyancy that the can must provide'}. Each 0.1 on the tension factor costs ${fmt(dT / 1e3)} kN (${fmt(dT / 1e3 / G)} t) of ${isT ? 'tensioner capacity and deck payload' : 'net can buoyancy'}.`,
        `Set the tension factor at the lowest value that meets the joint-angle and minimum-tension checks with the heaviest contents; ${isT ? 'tensioner cylinders and the hull displacement that carries them scale with it' : 'can volume and its installation scale with it'}.`,
        'Whole-life costing (ISO 15686-5): payload-driven hull and equipment quantities; AACE RP 18R-97 class 4 until detailed design.',
        `Riser string: ${steelCo2Text(steel)} Lower tension also reduces wall stress, which allows a thinner or lower-grade string.`),
    ];
    return finish({
      kpis, checks, plots, tables,
      notes: [`Apparent weight ${fmt(sec.w)} N/m over ${L} m = ${fmt(W / 1e3)} kN. ${gov.name} governs (utilisation ${fmt(ucMax)}).`,
        isT ? `Tensioner modelled as constant tension${inp.tensionerK > 0 ? ` plus ${inp.tensionerK} kN/m of stroke stiffness` : ''}.` : `Jumper catenary iterated with the tower in ${extra.iterations} passes.`],
      advice,
      publish: { 'riser.ttrTopTension': r.topTension, 'riser.ttrStroke': r.stroke, 'riser.ttrLowerAngleDeg': Math.abs(r.lowerAngleDeg) },
      quantities: { steel_t: steel, pipe_km: L / 1e3 },
    });
  },

  verify: [
    {
      // Vertical riser, no offset, no current: effective tension at the base = top tension − apparent weight
      //   Te,base = (TF − 1)·w·L with w from the hand weight formula (339.7 × 15.9 mm, contents 800 kg/m³).
      // True wall tension at the base (pi = ρi·g·L at zero top pressure, pe = ρw·g·d):
      //   Tw = Te + pi·Ai − pe·Ae (Sparks 2007 Eq. 2.5).
      name: 'Vertical tensioned string: Te,base = (TF − 1)·w·L and Tw = Te + pi·Ai − pe·Ae', source: 'Sparks, Fundamentals of Marine Riser Mechanics (2007) Ch. 2, Eq. 2.5',
      inputs: { riserType: 'ttr', offset: 0, Us: 0, Ub: 0, pInt: 0, tensionerK: 0, kLower: 0 },
      expect: [{ kpi: 'bottomTension', value: 0.4 * W_TTR * 980 / 1e3, tol: 0.002 },
        { kpi: 'TwBottom', value: (0.4 * W_TTR * 980 + 800 * G * 980 * Math.PI / 4 * 0.3079 ** 2 - 1025 * G * 1000 * Math.PI / 4 * 0.3397 ** 2) / 1e3, tol: 0.005 }],
    },
  ],
};
