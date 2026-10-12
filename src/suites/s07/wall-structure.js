// Structural strength of a reinforced-concrete wall stem (seawall, caisson front wall, intake chamber wall) under
// hydrostatic, wave and earth pressure: beam action over the height with the chosen supports, section
// resistance to EN 1992-1-1, deflection, and plate action of a wall panel between cross walls.
import { navierPlate } from '../../engine/structural/plates.js';
import { CO2_CONC, CO2_STEEL, need, run, finite, r, f, linspace, thin, check, eng, eco, slrInputs, slrRows, SLR_BASIS } from './_common.js';

/** Pressure on the stem at height z above the base [Pa] (all components act in the same direction). */
function pressureAt(z, o) {
  let p = 0;
  if (z < o.hw) p += o.gw * (o.hw - z) * o.headShare;                      // net hydrostatic head
  if (o.p1 > 0) {                                                            // wave: base → still water → η*
    if (z <= o.hw) p += o.hw > 0 ? o.p3 + (o.p1 - o.p3) * z / o.hw : o.p1;
    else if (z < o.hw + o.eta) p += o.p1 * (1 - (z - o.hw) / o.eta);
  }
  if (z < o.he) p += o.K * (o.gs * (o.he - z) + o.q);                        // earth and surcharge
  return p;
}

/** Bending moment, shear and deflection of the stem by numerical integration (trapezoidal rule on n elements). */
function solveStem({ H, n, p, EI, base, top }) {
  const dz = H / n, z = linspace(0, H, n + 1), pz = z.map(p);
  // cantilever from the free top: V(z) = ∫_z^H p ds, M(z) = ∫_z^H V ds
  const V0 = new Array(n + 1).fill(0), M0 = new Array(n + 1).fill(0);
  for (let i = n - 1; i >= 0; i--) { V0[i] = V0[i + 1] + 0.5 * (pz[i] + pz[i + 1]) * dz; M0[i] = M0[i + 1] + 0.5 * (V0[i] + V0[i + 1]) * dz; }
  const defl = (M, theta0) => { const th = [theta0], w = [0]; for (let i = 0; i < n; i++) { const k0 = M[i] / EI(z[i]), k1 = M[i + 1] / EI(z[i + 1]); th.push(th[i] + 0.5 * (k0 + k1) * dz); w.push(w[i] + 0.5 * (th[i] + th[i + 1]) * dz); } return { th, w }; };
  let R = 0, Mb0 = 0;
  if (base === 'fixed' && top === 'propped') {
    // force method: top prop reaction R removes the tip deflection; unit prop gives M = −(H − z)
    const d0 = defl(M0, 0).w[n], d1 = defl(z.map((zz) => H - zz), 0).w[n];
    R = d0 / d1;
  } else if (base === 'pinned' && top === 'propped') {
    R = M0[0] / H; // moment equilibrium about the pinned base
  }
  const M = M0.map((m, i) => m - R * (H - z[i])), V = V0.map((v) => v - R);
  let w;
  if (base === 'pinned') { const a = defl(M, 0); const th0 = -a.w[n] / H; w = a.w.map((v, i) => v + th0 * z[i]); Mb0 = 0; } else { w = defl(M, 0).w; Mb0 = M[0]; }
  return { z, pz, V, M, w, R, Mbase: Mb0 };
}

export default {
  id: 's07-wall-structure', suite: 7, title: 'Wall stem strength and deflection',
  summary: 'Bending, shear and deflection of a concrete wall stem under water, wave and earth pressure, checked against reinforced-concrete resistance, with plate action of a panel between cross walls.',
  covers: [1], fidelity: 'numerical',
  standards: ['EN 1992-1-1:2004 §6.1, §6.2.2, §7.2–7.4', 'EN 1990:2002 Annex A1 (load factors)', 'BS 6349-2:2019 / BS 6349-1-4:2021 (maritime concrete)', 'Timoshenko & Woinowsky-Krieger, Theory of Plates and Shells (1959) §28–30'],
  models: ['Euler–Bernoulli beam theory (numerical integration, force method)', 'Hydrostatic pressure', 'Goda-type linear wave pressure diagram', 'Rankine earth pressure', 'Rectangular concrete stress block', 'EN 1992 shear resistance without links', 'Kirchhoff plate theory (Navier series)'],
  limits: 'One-way vertical strip of unit width with linear-elastic stiffness (a cracked-stiffness factor is applied uniformly), plus a simply supported plate check of one panel; two-way action near corners, dynamic amplification under impulsive wave loads, crack-width calculation, fatigue and detailing are not covered. Nonlinear and 3-D finite-element analysis are outside the module.',
  publishes: ['marine.stemMoment', 'marine.stemShear'],
  searchTerms: 'reinforced concrete seawall stem bending shear wave pressure cantilever wall design',
  inputs: [
    { key: 'H', label: 'Stem height', unit: 'm', type: 'number', default: 7, min: 0.5, max: 40, group: 'Geometry', role: 'geometry' },
    { key: 'tb', label: 'Thickness at the base', unit: 'mm', type: 'number', default: 1000, min: 150, max: 5000, group: 'Geometry', role: 'geometry', design: true },
    { key: 'tt', label: 'Thickness at the top', unit: 'mm', type: 'number', default: 500, min: 150, max: 5000, group: 'Geometry', role: 'geometry', design: true },
    { key: 'base', label: 'Support at the base', type: 'select', default: 'fixed', options: [{ value: 'fixed', label: 'Fixed (monolithic with the base slab)' }, { value: 'pinned', label: 'Pinned (hinged or keyed)' }], group: 'Supports', role: 'bc' },
    { key: 'top', label: 'Support at the top', type: 'select', default: 'free', options: [{ value: 'free', label: 'Free (cantilever)' }, { value: 'propped', label: 'Propped (deck slab, tie or strut)' }], group: 'Supports', role: 'bc' },
    { key: 'conc', label: 'Concrete', type: 'select', source: 'materials:concrete', default: 'c40', group: 'Materials', role: 'material' },
    { key: 'fyk', label: 'Reinforcement yield strength', unit: 'MPa', type: 'number', default: 500, min: 250, max: 700, group: 'Materials', role: 'material', uncertain: { dist: 'normal', cov: 0.05 } },
    { key: 'dia', label: 'Main bar diameter', unit: 'mm', type: 'number', default: 32, min: 8, max: 50, group: 'Materials', role: 'geometry', design: true },
    { key: 'spacing', label: 'Main bar spacing', unit: 'mm', type: 'number', default: 140, min: 75, max: 400, group: 'Materials', role: 'geometry', design: true },
    { key: 'cover', label: 'Cover to main bars', unit: 'mm', type: 'number', default: 65, min: 25, max: 120, group: 'Materials', role: 'geometry', help: 'XS3 splash-zone exposure typically needs 55–75 mm.' },
    { key: 'stiff', label: 'Cracked / gross stiffness', unit: '-', type: 'number', default: 0.5, min: 0.2, max: 1, group: 'Materials', role: 'material', help: '1.0 uncracked; about 0.5 for a cracked wall at service load.' },
    { key: 'hw', label: 'Water level above the base on the loaded side', unit: 'm', type: 'number', default: 4.5, min: 0, max: 40, group: 'Loads', role: 'bc', help: 'Design still-water level before sea-level rise.' },
    { key: 'headShare', label: 'Share of the head that is unbalanced', unit: '-', type: 'number', default: 0.4, min: 0, max: 1, group: 'Loads', role: 'bc', help: '1 when the other side is dry (dewatered chamber, wave trough); 0 when water levels are equal.' },
    { key: 'p1', label: 'Wave pressure at still-water level', unit: 'kPa', type: 'number', default: 40, min: 0, max: 1000, group: 'Loads', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'p1 from the wave-pressure module (Goda).' },
    { key: 'p3', label: 'Wave pressure at the base of the stem', unit: 'kPa', type: 'number', default: 30, min: 0, max: 1000, group: 'Loads', role: 'bc' },
    { key: 'eta', label: 'Height above still water where wave pressure vanishes', unit: 'm', type: 'number', default: 5, min: 0.1, max: 40, group: 'Loads', role: 'bc', help: 'η* of the Goda diagram.' },
    { key: 'he', label: 'Retained soil height on the loaded side', unit: 'm', type: 'number', default: 0, min: 0, max: 40, group: 'Loads', role: 'bc' },
    { key: 'K', label: 'Earth-pressure coefficient', unit: '-', type: 'number', default: 0.33, min: 0.1, max: 1, group: 'Loads', role: 'soil' },
    { key: 'gs', label: 'Soil unit weight (effective)', unit: 'kN/m³', type: 'number', default: 18, min: 5, max: 23, group: 'Loads', role: 'soil' },
    { key: 'q', label: 'Surcharge on the retained soil', unit: 'kPa', type: 'number', default: 10, min: 0, max: 200, group: 'Loads', role: 'bc', help: 'Traffic or storage load.' },
    ...slrInputs({ group: 'Loads' }),
    { key: 'gF', label: 'Load factor (ULS)', unit: '-', type: 'number', default: 1.35, min: 1, max: 2, group: 'Criteria', role: 'criteria' },
    { key: 'gC', label: 'Concrete partial factor γc', unit: '-', type: 'number', default: 1.5, min: 1, max: 2, group: 'Criteria', role: 'criteria' },
    { key: 'gS', label: 'Steel partial factor γs', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.5, group: 'Criteria', role: 'criteria' },
    { key: 'acc', label: 'Long-term factor αcc', unit: '-', type: 'number', default: 0.85, min: 0.8, max: 1, group: 'Criteria', role: 'criteria', help: 'EN 1992 recommends 1.0; several national annexes use 0.85.' },
    { key: 'sLim', label: 'Service steel-stress limit (crack control)', unit: 'MPa', type: 'number', default: 200, min: 100, max: 400, group: 'Criteria', role: 'criteria', help: 'EN 1992-1-1 Table 7.2N / 7.3N: about 200 MPa for 0.2–0.3 mm cracks with 25 mm bars.' },
    { key: 'deflLim', label: 'Deflection limit (height / …)', unit: '-', type: 'number', default: 250, min: 50, max: 1000, group: 'Criteria', role: 'criteria' },
    { key: 'panelW', label: 'Panel width between cross walls', unit: 'm', type: 'number', default: 4.5, min: 0.5, max: 30, group: 'Geometry', role: 'geometry', help: 'Caisson cell or counterfort spacing for the plate check.' },
    { key: 'length', label: 'Length of wall', unit: 'm', type: 'number', default: 300, min: 1, max: 20000, group: 'Geometry', role: 'geometry' },
    { key: 'n', label: 'Elements over the height', type: 'int', default: 40, min: 4, max: 1280, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'n', label: 'Elements', levels: [10, 20, 40, 80], qoi: 'defl' },
  compute(inp, ctx) {
    need(!(inp.base === 'pinned' && inp.top === 'free'), 'A stem pinned at the base and free at the top is a mechanism (under-constrained): fix the base or prop the top.');
    const c = ctx.material(inp.conc); need(c && c.fck > 0, 'The selected concrete is not in the materials library.');
    const H = inp.H, n = Math.round(inp.n), tb = inp.tb / 1000, tt = inp.tt / 1000;
    const tAt = (z) => tb + (tt - tb) * z / H;
    const EI = (z) => inp.stiff * c.E * tAt(z) ** 3 / 12;
    const loadOf = (rise) => ({ hw: Math.min(inp.hw + rise, H), gw: 10055, headShare: inp.headShare, p1: inp.p1 * 1e3, p3: inp.p3 * 1e3, eta: inp.eta, he: Math.min(inp.he, H), K: inp.K, gs: inp.gs * 1e3, q: inp.he > 0 ? inp.q * 1e3 : 0 });
    const solve = (rise, nn = n) => solveStem({ H, n: nn, p: (z) => pressureAt(z, loadOf(rise)), EI, base: inp.base, top: inp.top });
    const s = run(() => solve(inp.slr), 'Stem analysis');
    need(Math.max(...s.pz) > 0, 'No pressure acts on the stem: give a water level, wave pressure or retained soil.');
    // section resistance per metre width (EN 1992-1-1)
    const fcd = inp.acc * c.fck / inp.gC, fyd = inp.fyk * 1e6 / inp.gS, As = Math.PI * (inp.dia / 1000) ** 2 / 4 / (inp.spacing / 1000);
    const sect = (z) => {
      const d = tAt(z) - inp.cover / 1000 - inp.dia / 2000;
      need(d > 0.05, 'The section is too thin for the cover and bar diameter.');
      const x = As * fyd / (0.8 * fcd), MRd = x < 0.45 * d ? As * fyd * (d - 0.4 * x) : 0.8 * 0.45 * d * fcd * (d - 0.18 * d); // over-reinforced sections limited to x = 0.45 d
      const k = Math.min(2, 1 + Math.sqrt(0.2 / d)), rho = Math.min(0.02, As / d), fck = c.fck / 1e6;
      const vRd = Math.max(0.18 / inp.gC * k * Math.cbrt(100 * rho * fck), 0.035 * k ** 1.5 * Math.sqrt(fck)) * 1e6; // EN 1992-1-1 Eq. 6.2a, 6.3N
      return { d, MRd, VRd: vRd * d, x };
    };
    // required steel for a moment: 0.8 x fcd (d − 0.4 x) = M
    const asReq = (M, d) => { const a = 0.32 * fcd, b = -0.8 * fcd * d, disc = b * b - 4 * a * M; if (disc < 0) return NaN; const x = (-b - Math.sqrt(disc)) / (2 * a); return 0.8 * x * fcd / fyd; };
    let ucM = 0, ucV = 0, zM = 0, zV = 0;
    s.z.forEach((z, i) => { const q = sect(z), um = inp.gF * Math.abs(s.M[i]) / q.MRd, uv = inp.gF * Math.abs(s.V[i]) / q.VRd; if (um > ucM) { ucM = um; zM = z; } if (uv > ucV) { ucV = uv; zV = z; } });
    const iM = s.M.reduce((b, m, i) => (Math.abs(m) > Math.abs(s.M[b]) ? i : b), 0), Mmax = Math.abs(s.M[iM]), Vmax = Math.max(...s.V.map(Math.abs));
    const qM = sect(s.z[iM]);
    const AsReq = asReq(inp.gF * Mmax, qM.d);
    const sigS = Mmax / (As * 0.9 * qM.d); // service steel stress with lever arm ≈ 0.9 d
    const defl = Math.max(...s.w.map(Math.abs));
    finite(defl, 'Deflection');
    const dLim = H / inp.deflLim;
    // plate action: panel between cross walls, simply supported on four edges, under the peak pressure
    const pPeak = Math.max(...s.pz), tPanel = tAt(Math.min(inp.hw + inp.slr, H));
    const pl = run(() => navierPlate({ a: inp.panelW, b: H, t: tPanel, E: c.E, nu: c.nu ?? 0.2, load: { type: 'uniform', q: pPeak }, nTerms: 39 }), 'Plate check');
    const mPlate = Math.max(Math.abs(pl.mx), Math.abs(pl.my)), ucPlate = inp.gF * mPlate / sect(Math.min(inp.hw + inp.slr, H)).MRd;
    // thickness at the base that would bring the governing utilisation to 1 (resistance ∝ d for shear, ∝ d for lightly reinforced bending)
    const ucGov = Math.max(ucM, ucV), tReq = ucGov > 1 ? tb * ucGov : tb;
    const vol = 0.5 * (tb + tt) * H, concAll = vol * inp.length, steelT = As * H * inp.length * 7850 * 1.6 / 1000;
    const scen = slrRows(inp).map((sc) => { const q = solve(sc.rise, Math.min(n, 40)); const Mx = Math.max(...q.M.map(Math.abs)); return { name: sc.name, rise: sc.rise, M: r(Mx / 1e3), V: r(Math.max(...q.V.map(Math.abs)) / 1e3), uc: r(inp.gF * Mx / qM.MRd, 3) }; });
    const ok = ucM <= 1 && ucV <= 1 && sigS / 1e6 <= inp.sLim && defl <= dLim;
    const [pz, pp, pV, pMm, pw] = thin([s.z, s.pz, s.V, s.M, s.w], 200);

    const advice = [
      ok
        ? eng('ok', `The stem carries M = ${f(Mmax / 1e3)} kN·m/m and V = ${f(Vmax / 1e3)} kN/m: bending utilisation ${f(ucM, 2)}, shear ${f(ucV, 2)}, service steel stress ${f(sigS / 1e6)} MPa, deflection ${f(defl * 1000)} mm (H/${f(H / defl)}).`,
          `Steel required is ${f(AsReq * 1e6)} mm²/m against ${f(As * 1e6)} mm²/m provided. ${ucGov < 0.7 ? `The base could be thinned towards ${f(tb * Math.max(ucGov, 0.5) * 1000)} mm if crack control and durability cover allow.` : 'Little reserve remains: keep the section.'}`,
          'EN 1992-1-1:2004 §6.1, §6.2.2, §7.3; EN 1990 Annex A1', `The stem holds ${f(vol, 3)} m³ of concrete per metre run; each 100 mm off the base thickness saves ${f(0.05 * H * inp.length, 3)} m³ over the wall (≈ ${f(0.05 * H * inp.length * CO2_CONC, 3)} t CO₂e, screening factor).`)
        : eng('fail', `The stem is overstressed: bending ${f(ucM, 2)} (at ${f(zM)} m), shear ${f(ucV, 2)} (at ${f(zV)} m), service steel stress ${f(sigS / 1e6)} MPa (limit ${f(inp.sLim)}), deflection ${f(defl * 1000)} mm (limit ${f(dLim * 1000)}).`,
          `${ucM > 1 ? `Provide at least ${Number.isFinite(AsReq) ? `${f(AsReq * 1e6)} mm²/m` : 'a deeper section — the compression zone is exhausted'} of main steel (now ${f(As * 1e6)}). ` : ''}${ucV > 1 ? `Thicken the base to about ${f(tReq * 1000)} mm or add shear links. ` : ''}${sigS / 1e6 > inp.sLim ? `For crack control use about ${f(As * sigS / (inp.sLim * 1e6) * 1e6)} mm²/m. ` : ''}${defl > dLim ? 'Prop the top or thicken the stem to limit deflection.' : ''}`,
          'EN 1992-1-1:2004 §6.1, §6.2.2, §7.3; EN 1990 Annex A1', `Thickening the base to ${f(tReq * 1000)} mm adds ${f(0.5 * (tReq - tb) * H * inp.length, 3)} m³ of concrete (≈ ${f(0.5 * (tReq - tb) * H * inp.length * CO2_CONC, 3)} t CO₂e); propping the top with a deck or tie cuts the base moment several-fold with far less material (screening factors).`),
      eng(ucPlate > 1 ? 'warn' : 'info', `Panel between cross walls (${f(inp.panelW)} × ${f(H)} m, ${f(tPanel * 1000)} mm): plate moment ${f(mPlate / 1e3)} kN·m/m under the ${f(pPeak / 1e3)} kPa peak pressure, utilisation ${f(ucPlate, 2)} with the vertical steel.`,
        ucPlate > 1 ? 'Provide horizontal reinforcement for the plate moment or reduce the cross-wall spacing; check the supports of the panel for the reactions.' : 'Horizontal steel of similar area to the vertical steel covers two-way action; confirm with a slab analysis where the panel is continuous over cross walls.', 'Timoshenko & Woinowsky-Krieger (1959) §30; EN 1992-1-1 §9.3'),
      eng(scen[scen.length - 1].uc > 1 ? 'warn' : 'info', `Sea-level rise: the base moment goes from ${f(scen[0].M)} to ${f(scen[scen.length - 1].M)} kN·m/m across the scenarios (bending utilisation ${f(scen[scen.length - 1].uc, 2)} at +${f(scen[scen.length - 1].rise)} m).`,
        'Hydrostatic moment grows with the cube of the head: design the stem for the end-of-life water level, and detail starter bars so that a future wave wall can be added.', SLR_BASIS),
      eco('info', `The stem takes ${f(concAll, 3)} m³ of concrete and roughly ${f(steelT, 3)} t of reinforcement over ${f(inp.length)} m (screening take-off with laps and distribution steel).`,
        'Optimise thickness against reinforcement: a thinner, more heavily reinforced stem saves concrete but raises steel stress and crack width; in the splash zone durability cover usually fixes the minimum thickness.', 'AACE RP 18R-97 quantity-based estimating; whole-life costing per ISO 15686-5',
        `Embodied carbon ≈ ${f(concAll * CO2_CONC + steelT * CO2_STEEL, 3)} t CO₂e (0.3 t/m³ concrete, 1.9 t/t steel, screening factors); cement replacement lowers the concrete share and improves chloride resistance.`),
    ];

    return {
      kpis: [
        { key: 'Mmax', label: 'Maximum bending moment (characteristic)', value: r(Mmax / 1e3), unit: 'kN·m/m', explain: 'Largest moment in the stem per metre width, before the load factor. At the base for a cantilever; shared between base and span for a propped stem.' },
        { key: 'Mbase', label: 'Moment at the base', value: r(Math.abs(s.Mbase) / 1e3), unit: 'kN·m/m', explain: 'Fixed-end moment transferred into the base slab; zero for a pinned base.' },
        { key: 'Vmax', label: 'Maximum shear force (characteristic)', value: r(Vmax / 1e3), unit: 'kN/m', explain: 'Largest shear in the stem per metre width, normally at the base.' },
        { key: 'Rtop', label: 'Prop reaction at the top', value: r(s.R / 1e3), unit: 'kN/m', explain: 'Force the deck, tie or strut must provide at the top of the stem; zero for a free top.' },
        { key: 'defl', label: 'Maximum deflection', value: r(defl * 1000), unit: 'mm', limit: r(dLim * 1000), sense: 'max', explain: 'Largest lateral displacement of the stem with the cracked stiffness factor. It governs joint sealing and the fit of gates or parapets.' },
        { key: 'ucM', label: 'Bending utilisation MEd / MRd', value: r(ucM, 3), unit: '-', limit: 1, sense: 'max', explain: 'Factored moment over the bending resistance of the reinforced section, worst over the height. Above 1 the tension steel yields and the stem hinges.' },
        { key: 'ucV', label: 'Shear utilisation VEd / VRd,c', value: r(ucV, 3), unit: '-', limit: 1, sense: 'max', explain: 'Factored shear over the resistance without links. Above 1 a brittle diagonal shear failure is possible unless links are provided.' },
        { key: 'sigS', label: 'Service steel stress', value: r(sigS / 1e6), unit: 'MPa', limit: inp.sLim, sense: 'max', explain: 'Tension in the main bars under unfactored load. It controls crack width, which controls chloride ingress and the life of the wall.' },
        { key: 'AsReq', label: 'Main steel required (ULS)', value: r(Number.isFinite(AsReq) ? AsReq * 1e6 : 99999), unit: 'mm²/m', explain: 'Reinforcement area that just resists the factored maximum moment; 99999 means the concrete compression zone is exhausted and the section must be deeper.' },
        { key: 'AsProv', label: 'Main steel provided', value: r(As * 1e6), unit: 'mm²/m', explain: 'Area of the specified bar diameter and spacing per metre width.' },
        { key: 'ucPlate', label: 'Panel plate-bending utilisation', value: r(ucPlate, 3), unit: '-', explain: 'Navier plate moment of a panel between cross walls over the section resistance. Shows whether two-way action needs its own horizontal steel.' },
        { key: 'pPeak', label: 'Peak design pressure on the stem', value: r(pPeak / 1e3), unit: 'kPa', explain: 'Largest combined hydrostatic, wave and earth pressure ordinate, unfactored.' },
      ],
      checks: [
        check('Bending resistance', ucM, 1, '-', 'EN 1992-1-1 §6.1', 'Guards against yielding of the tension reinforcement and crushing of the compression zone. Exceeded, the stem forms a hinge and rotates.'),
        check('Shear resistance without links', ucV, 1, '-', 'EN 1992-1-1 §6.2.2 Eq. 6.2a', 'Guards against sudden diagonal shear failure of a wall without shear reinforcement — the most dangerous mode for thick stems.'),
        check('Service steel stress (crack control)', sigS / 1e6, inp.sLim, 'MPa', 'EN 1992-1-1 §7.3.3 Tables 7.2N, 7.3N', 'Guards against wide cracks that let chlorides reach the reinforcement. Exceeded, durability, not strength, is lost first.'),
        check('Deflection', defl * 1000, dLim * 1000, 'mm', 'EN 1992-1-1 §7.4; project serviceability criterion', 'Guards against visible lean and opening of joints. Exceeded, seals fail and attached equipment loses alignment.'),
      ],
      plots: [
        { type: 'line', title: 'Pressure, shear and moment over the height', x: { label: 'Pressure (kPa) / shear (kN/m) / moment (kN·m/m)', unit: '' }, y: { label: 'Height above base', unit: 'm' },
          series: [{ name: 'Pressure', x: pp.map((v) => r(v / 1e3, 3)), y: pz.map((v) => r(v, 3)) }, { name: 'Shear', x: pV.map((v) => r(v / 1e3, 3)), y: pz.map((v) => r(v, 3)) }, { name: 'Moment', x: pMm.map((v) => r(v / 1e3, 3)), y: pz.map((v) => r(v, 3)) }],
          explain: 'Loading and internal forces from base to top. For a cantilever both shear and moment peak at the base; a propped stem shows a sign change in the moment.' },
        { type: 'line', title: 'Moment against resistance', x: { label: 'Height above base', unit: 'm' }, y: { label: 'Moment', unit: 'kN·m/m' },
          series: [{ name: 'Factored |MEd|', x: pz.map((v) => r(v, 3)), y: pMm.map((v) => r(inp.gF * Math.abs(v) / 1e3, 3)) }, { name: 'MRd', x: pz.map((v) => r(v, 3)), y: pz.map((v) => r(sect(v).MRd / 1e3, 3)), dash: true }],
          explain: 'Design moment envelope against the resistance of the tapering section. Where the demand curve approaches the resistance line is where bars can or cannot be curtailed.' },
        { type: 'line', title: 'Deflected shape', x: { label: 'Deflection', unit: 'mm' }, y: { label: 'Height above base', unit: 'm' },
          series: [{ name: 'Deflection', x: pw.map((v) => r(v * 1000, 3)), y: pz.map((v) => r(v, 3)) }],
          explain: 'Lateral displacement of the stem under unfactored load with the cracked stiffness. The tip value is checked against the limit.' },
      ],
      tables: [
        { title: 'Sea-level-rise scenarios', columns: [{ key: 'name', label: 'Scenario' }, { key: 'rise', label: 'Rise', unit: 'm' }, { key: 'M', label: 'Max moment', unit: 'kN·m/m' }, { key: 'V', label: 'Max shear', unit: 'kN/m' }, { key: 'uc', label: 'Bending utilisation', unit: '-' }],
          rows: scen, explain: 'Internal forces with the water level raised by each scenario and the bending utilisation at the critical section. A value above 1 marks the scenario at which strengthening is needed.' },
        { title: 'Critical section', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['Height of maximum moment', s.z[iM], 'm'], ['Effective depth', qM.d * 1000, 'mm'], ['Neutral-axis depth', qM.x * 1000, 'mm'], ['MRd', qM.MRd / 1e3, 'kN·m/m'], ['VRd,c', qM.VRd / 1e3, 'kN/m'], ['fcd', fcd / 1e6, 'MPa'], ['fyd', fyd / 1e6, 'MPa'], ['Reinforcement ratio', As / qM.d * 100, '%']].map(([q, v, u]) => ({ q, v: r(v), u })),
          explain: 'Section properties and resistances at the governing section, for comparison with a hand check.' },
      ],
      notes: [`${inp.base === 'fixed' ? 'Fixed' : 'Pinned'} base, ${inp.top === 'free' ? 'free' : 'propped'} top; ${ucM >= ucV ? 'bending' : 'shear'} governs the strength.`, `Water ${f(Math.min(inp.hw + inp.slr, H))} m above the base including ${f(inp.slr)} m of sea-level rise; ${f(inp.headShare * 100)} % of the head unbalanced.`],
      advice,
      publish: { 'marine.stemMoment': Mmax, 'marine.stemShear': Vmax },
      quantities: { concrete_m3: r(concAll), steel_t: r(steelT) },
    };
  },
  verify: [
    { name: 'Cantilever under triangular hydrostatic load: M = γw h³/6, V = γw h²/2, δ = q0 L⁴/(30 EI)', source: 'Roark\'s Formulas for Stress and Strain, Table 8.1 case 2d. h = H = 6 m, γw = 10.055 kN/m³: q0 = 60.33 kPa, M = 361.98 kN·m/m, V = 180.99 kN/m; C40 (E = 35 GPa), t = 0.6 m, EI = 6.30e8 N·m² → δ = 60330 × 1296/(30 × 6.30e8) = 4.137 mm',
      inputs: { H: 6, tb: 600, tt: 600, base: 'fixed', top: 'free', conc: 'c40', stiff: 1, hw: 6, headShare: 1, p1: 0, p3: 0, he: 0, slr: 0, n: 80 },
      expect: [{ kpi: 'Mmax', value: 361.98, tol: 0.002 }, { kpi: 'Vmax', value: 180.99, tol: 0.002 }, { kpi: 'defl', value: 4.137, tol: 0.005 }] },
    { name: 'Propped cantilever under triangular load: prop reaction q0 L/10, fixed-end moment q0 L²/15', source: 'Roark\'s Formulas for Stress and Strain, Table 8.1 case 2c. q0 = 60.33 kPa, L = 6 m → R = 36.20 kN/m, Mbase = 144.79 kN·m/m',
      inputs: { H: 6, tb: 600, tt: 600, base: 'fixed', top: 'propped', conc: 'c40', stiff: 1, hw: 6, headShare: 1, p1: 0, p3: 0, he: 0, slr: 0, n: 80 },
      expect: [{ kpi: 'Rtop', value: 36.20, tol: 0.005 }, { kpi: 'Mbase', value: 144.79, tol: 0.005 }] },
  ],
};
