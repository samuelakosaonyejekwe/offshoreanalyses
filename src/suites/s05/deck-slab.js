// Deck, beam and slab loading of a wharf: a slab panel under uniform cargo, container corner loads and truck or
// reach-stacker wheel patches (Navier plate solution), punching under the heaviest patch, cracking of the slab,
// the crane rail beam under a moving wheel group (influence lines by the stiffness method), and the permissible
// cargo, wheel and crane wheel loads.
import { navierPlate } from '../../engine/structural/plates.js';
import { analyzeFrame, memberDiagram } from '../../engine/structural/frame3d.js';
import { momentCapacity, shearResistanceNoLinks, shearResistanceLinks, crackedSection, crackWidth } from '../../engine/structural/concrete.js';
import { need, run, r, f, check, eng, eco, clamp, rcInputs, rcMember, concCo2, CONFIRM, NOT_HIFI, RC_NOTE, ECON_BASIS } from './_common.js';

const PATCH_DEFAULT = [
  { name: 'Truck axle, wheel 1', P: 75, u: 0.4, v: 0.4, x: 2.1, y: 3.0, dyn: 1.25 },
  { name: 'Truck axle, wheel 2', P: 75, u: 0.4, v: 0.4, x: 3.9, y: 3.0, dyn: 1.25 },
  { name: 'Container corner (4 stacks meet)', P: 300, u: 0.36, v: 0.36, x: 3.0, y: 1.2, dyn: 1.0 },
];

export default {
  id: 's05-deck-slab', suite: 5, title: 'Deck slab and crane rail beam loading',
  summary: 'Slab moments, stresses and reinforcement demand under cargo, container and wheel loads, punching and crack width, the crane rail beam under a moving wheel group, and the permissible cargo, wheel and crane wheel loads.',
  covers: [7], fidelity: 'analytical',
  standards: ['EN 1992-1-1:2004 §6.1, §6.2, §6.4, §7.3', 'EN 1991-2:2003 (traffic loads)', 'BS 6349-1-2:2016 §8 (cargo and equipment actions)', 'ISO 1496-1 (container corner fitting loads)', 'FEM 1.001 / EN 13001 (crane wheel loads from the supplier)'],
  models: ['Kirchhoff–Love plate theory (Navier double series)', 'Reinforced-concrete sectional analysis (parabola–rectangle)', 'Cracked-section analysis and crack width', 'Punching shear on a control perimeter', 'Influence lines of a continuous beam (matrix stiffness method)', 'Timoshenko beam elements'],
  limits: `The slab panel is analysed as a simply supported thin elastic plate with factors for continuity over the beams; thick-plate (Mindlin–Reissner) shear deformation, membrane action and yield-line reserve are not included. The rail beam is a continuous beam on unyielding supports. ${RC_NOTE} ${NOT_HIFI}`,
  publishes: ['wharf.cargoPermissible', 'wharf.craneWheelPermissible', 'wharf.wheelPermissible'],
  searchTerms: 'wharf deck slab wheel load container corner load crane rail beam punching',
  inputs: [
    { key: 'a', label: 'Slab panel span across the wharf', unit: 'm', type: 'number', default: 6, min: 1, max: 15, group: 'Slab panel', role: 'geometry', design: true },
    { key: 'b', label: 'Slab panel span along the wharf', unit: 'm', type: 'number', default: 6, min: 1, max: 15, group: 'Slab panel', role: 'geometry' },
    { key: 'slabT', label: 'Slab thickness', unit: 'm', type: 'number', default: 0.45, min: 0.15, max: 1.5, group: 'Slab panel', role: 'geometry', design: true },
    { key: 'conc', label: 'Concrete class', type: 'select', source: 'materials:concrete', default: 'c40', group: 'Slab panel', role: 'material' },
    { key: 'gConc', label: 'Unit weight of reinforced concrete', unit: 'kN/m³', type: 'number', default: 25, min: 0, max: 28, group: 'Slab panel', role: 'material' },
    { key: 'nu', label: "Poisson's ratio of the slab", unit: '-', type: 'number', default: 0.2, min: 0, max: 0.35, group: 'Slab panel', role: 'material', help: '0.2 for uncracked concrete, 0 for cracked (EN 1992-1-1 §3.1.3).' },
    { key: 'kSpan', label: 'Continuity factor on the span moment', unit: '-', type: 'number', default: 0.8, min: 0.5, max: 1, group: 'Slab panel', role: 'bc', help: 'Span moment of the continuous slab over that of the simply supported panel; 1.0 for an edge panel or precast planks without continuity.' },
    { key: 'kSup', label: 'Support (hogging) moment as a fraction of the simple span moment', unit: '-', type: 'number', default: 0.65, min: 0, max: 1, group: 'Slab panel', role: 'bc', help: 'Moment over the beams of a continuous slab; 0 for simply supported planks.' },
    { key: 'sBar', label: 'Bottom bar spacing', unit: 'mm', type: 'number', default: 150, min: 75, max: 400, group: 'Slab reinforcement', role: 'geometry', design: true },
    { key: 'dBar', label: 'Bottom bar diameter', unit: 'mm', type: 'number', default: 25, min: 10, max: 40, group: 'Slab reinforcement', role: 'geometry', design: true },
    { key: 'sBarTop', label: 'Top bar spacing', unit: 'mm', type: 'number', default: 150, min: 75, max: 400, group: 'Slab reinforcement', role: 'geometry' },
    { key: 'dBarTop', label: 'Top bar diameter', unit: 'mm', type: 'number', default: 25, min: 10, max: 40, group: 'Slab reinforcement', role: 'geometry' },
    { key: 'gSdl', label: 'Surfacing and services', unit: 'kPa', type: 'number', default: 2.5, min: 0, max: 30, group: 'Operational loads', role: 'ic' },
    { key: 'qCargo', label: 'Uniform cargo / container stacking load', unit: 'kPa', type: 'number', default: 40, min: 0, max: 300, group: 'Operational loads', role: 'bc', design: true, uncertain: { dist: 'gumbel', cov: 0.2 }, help: 'For example 20 kPa general cargo and personnel areas, 40–60 kPa containers stacked 3–5 high (averaged), 100+ kPa bulk and heavy lift.' },
    { key: 'patches', label: 'Wheel, axle and container corner loads on the panel', type: 'table', group: 'Operational loads', role: 'bc',
      columns: [{ key: 'name', label: 'Load', type: 'text' }, { key: 'P', label: 'Load on the patch', unit: 'kN', type: 'number' }, { key: 'u', label: 'Patch size across', unit: 'm', type: 'number' }, { key: 'v', label: 'Patch size along', unit: 'm', type: 'number' }, { key: 'x', label: 'x in the panel', unit: 'm', type: 'number' }, { key: 'y', label: 'y in the panel', unit: 'm', type: 'number' }, { key: 'dyn', label: 'Dynamic factor', unit: '-', type: 'number' }],
      default: PATCH_DEFAULT, help: 'One row per wheel or corner casting, positioned for the worst effect (near mid-panel for bending, near a beam for shear). Typical: truck axle 150 kN on two wheels 1.8 m apart; reach stacker front axle 4 × 250 kN; a container corner where four 30 t stacks meet about 300 kN per tier group.' },
    { key: 'patchWithCargo', label: 'Uniform cargo acts together with the patch loads', type: 'bool', default: false, group: 'Operational loads', role: 'bc', help: 'Normally the stacking area and the traffic lane are alternatives; tick only where both can load the same panel.' },
    { key: 'railB', label: 'Rail beam width', unit: 'm', type: 'number', default: 1.2, min: 0.3, max: 4, group: 'Crane rail beam', role: 'geometry' },
    { key: 'railH', label: 'Rail beam depth', unit: 'm', type: 'number', default: 1.6, min: 0.4, max: 4, group: 'Crane rail beam', role: 'geometry', design: true },
    { key: 'railL', label: 'Rail beam span (bent spacing)', unit: 'm', type: 'number', default: 6, min: 2, max: 20, group: 'Crane rail beam', role: 'geometry' },
    { key: 'nSpans', label: 'Continuous spans', unit: '-', type: 'int', default: 3, min: 1, max: 5, group: 'Crane rail beam', role: 'bc', help: 'Spans between expansion joints modelled; supports are the pile bents, taken as unyielding.' },
    { key: 'railTrib', label: 'Slab width carried by the rail beam', unit: 'm', type: 'number', default: 3, min: 0, max: 15, group: 'Crane rail beam', role: 'geometry' },
    { key: 'wheelP', label: 'Crane wheel load', unit: 'kN', type: 'number', default: 450, min: 0, max: 2500, group: 'Crane rail beam', role: 'bc', design: true, uncertain: { dist: 'normal', cov: 0.1 }, help: 'Largest characteristic wheel load from the crane supplier (operating with load and wind, or stowed storm case).' },
    { key: 'nWheels', label: 'Wheels in the corner group', unit: '-', type: 'int', default: 8, min: 1, max: 16, group: 'Crane rail beam', role: 'bc' },
    { key: 'wheelS', label: 'Wheel spacing', unit: 'm', type: 'number', default: 1.2, min: 0.5, max: 3, group: 'Crane rail beam', role: 'bc' },
    { key: 'dynCrane', label: 'Dynamic factor on wheel loads', unit: '-', type: 'number', default: 1.1, min: 1, max: 1.6, group: 'Crane rail beam', role: 'criteria', help: `Travelling and hoisting dynamics (EN 13001 / FEM). ${CONFIRM}` },
    ...rcInputs({ group: 'Rail beam reinforcement', nBot: 12, nTop: 12, links: 16, linkLegs: 6, linkS: 150 }),
    { key: 'gG', label: 'γG permanent actions', unit: '-', type: 'number', default: 1.35, min: 0.9, max: 1.6, group: 'Partial factors', role: 'criteria', help: CONFIRM },
    { key: 'gQ', label: 'γQ cargo, wheel and crane actions', unit: '-', type: 'number', default: 1.5, min: 1, max: 2, group: 'Partial factors', role: 'criteria', help: CONFIRM },
    { key: 'psi2', label: 'ψ2 quasi-permanent share of the cargo load', unit: '-', type: 'number', default: 0.8, min: 0, max: 1, group: 'Partial factors', role: 'criteria', help: `Storage areas (EN 1990 category E): 0.8. Wheel loads are taken as not quasi-permanent. ${CONFIRM}` },
    { key: 'wLim', label: 'Crack width limit', unit: 'mm', type: 'number', default: 0.3, min: 0.05, max: 0.5, group: 'Partial factors', role: 'criteria', help: `EN 1992-1-1 Table 7.1N for XS exposure under the quasi-permanent combination; some owners ask 0.2 mm. ${CONFIRM}` },
    { key: 'nTerms', label: 'Harmonics of the plate series', unit: '-', type: 'int', default: 48, min: 8, max: 400, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nTerms', label: 'Series harmonics per direction', levels: [12, 24, 48, 96], qoi: 'mEd' },
  compute(inp, ctx) {
    const a = inp.a, b = inp.b, t = inp.slabT, nT = Math.round(inp.nTerms);
    // slab strip of 1 m width with the bar spacing converted to bars per metre
    const slabInp = { ...inp, nBot: 1000 / inp.sBar, dBot: inp.dBar, nTop: 1000 / inp.sBarTop, dTop: inp.dBarTop, dLink: 0, linkLegs: 0 };
    const rc = rcMember(slabInp, ctx, 1, t), Ecm = rc.props.Ecm;
    const patches = (inp.patches || []).map((q, i) => {
      const p = { name: String(q.name || `Patch ${i + 1}`), P: Number(q.P) * 1e3, u: Number(q.u), v: Number(q.v), x: Number(q.x), y: Number(q.y), dyn: Number(q.dyn) > 0 ? Number(q.dyn) : 1 };
      need(p.P >= 0 && p.u > 0 && p.v > 0, `Patch load ${i + 1}: the load must not be negative and the patch sizes must be positive.`);
      need(p.x - p.u / 2 >= -1e-9 && p.x + p.u / 2 <= a + 1e-9 && p.y - p.v / 2 >= -1e-9 && p.y + p.v / 2 <= b + 1e-9, `Patch load '${p.name}' does not lie within the ${f(a)} × ${f(b)} m panel.`);
      return p;
    });
    // unit solutions: uniform 1 Pa and each patch at 1 N, evaluated at the panel centre and under every patch
    const pts = [{ x: a / 2, y: b / 2, name: 'Panel centre' }, ...patches.map((p) => ({ x: p.x, y: p.y, name: p.name }))];
    const plate = (load, x, y) => run(() => navierPlate({ a, b, t, E: Ecm, nu: inp.nu, load, x, y, nTerms: nT }), 'Plate analysis');
    const uni = pts.map((p) => plate({ type: 'uniform', q: 1 }, p.x, p.y));
    const pat = patches.map((p) => pts.map((q) => plate({ type: 'patch', P: 1, x0: p.x, y0: p.y, u: p.u, v: p.v }, q.x, q.y)));
    const gk = (t * inp.gConc + inp.gSdl) * 1e3, qk = inp.qCargo * 1e3;
    // effect at point j for factors on permanent, uniform cargo and patch loads
    const eff = (j, fg, fq, fp) => { let mx = (fg * gk + fq * qk) * uni[j].mx, my = (fg * gk + fq * qk) * uni[j].my, w = (fg * gk + fq * qk) * uni[j].w; patches.forEach((p, i) => { mx += fp * p.P * p.dyn * pat[i][j].mx; my += fp * p.P * p.dyn * pat[i][j].my; w += fp * p.P * p.dyn * pat[i][j].w; }); return { mx, my, w }; };
    const worstOf = (fg, fq, fp) => pts.map((_, j) => eff(j, fg, fq, fp)).reduce((m, e, j) => (Math.max(e.mx, e.my) > m.m ? { m: Math.max(e.mx, e.my), j, e } : m), { m: -Infinity, j: 0, e: null });
    const both = !!inp.patchWithCargo;
    // ULS: stacking case and traffic case (alternatives unless combined)
    const cases = both ? [{ name: 'Cargo with wheel and corner loads', fq: inp.gQ, fp: inp.gQ }] : [{ name: 'Uniform cargo', fq: inp.gQ, fp: 0 }, { name: 'Wheel and corner loads', fq: 0, fp: inp.gQ }];
    const uls = cases.map((c) => ({ ...c, ...worstOf(inp.gG, c.fq, c.fp) }));
    const gov = uls.reduce((m, c) => (c.m > m.m ? c : m), uls[0]);
    const mSS = gov.m, mEd = inp.kSpan * mSS, mSup = inp.kSup * mSS;
    const MRd = momentCapacity({ section: rc.sag, material: rc.material }).MRd, MRdH = momentCapacity({ section: rc.hog, material: rc.material }).MRd;
    const fyd = rc.fyk / rc.material.gammaS, AsReq = mEd / (0.9 * rc.d * fyd), AsReqTop = mSup / (0.9 * rc.dHog * fyd);
    const sigma = 6 * (inp.kSpan * worstOf(1, both ? 1 : 0, 1).m) / (t * t); // uncracked extreme-fibre stress, characteristic wheel case
    // punching under the heaviest patch (EN 1992-1-1 §6.4): basic control perimeter at 2d, mean of the two bar directions
    const dP = rc.d - rc.db / 2, heavy = patches.length ? patches.reduce((m, p) => (p.P * p.dyn > m.P * m.dyn ? p : m), patches[0]) : null;
    const vRdc = shearResistanceNoLinks({ bw: 1, d: dP, fck: rc.fck, Asl: rc.AsB, gammaC: inp.gammaC, alphaCC: inp.alphaCC }).VRdc / dP; // Pa on bw·d
    let vEd = 0, u1 = 0;
    if (heavy) { u1 = 2 * (heavy.u + heavy.v) + 4 * Math.PI * dP; const inside = (heavy.u + 4 * dP) * (heavy.v + 4 * dP) - (4 - Math.PI) * 4 * dP * dP; vEd = Math.max(inp.gQ * heavy.P * heavy.dyn - (inp.gG * gk) * inside, 0) / (u1 * dP); }
    // SLS quasi-permanent: cracking
    const mQp = inp.kSpan * worstOf(1, inp.psi2, 0).m;
    const cs = crackedSection({ section: rc.sag, Ec: Ecm, Es: rc.Es, M: Math.max(mQp, 1) });
    const cw = crackWidth({ sigmaS: cs.sigmaS, As: rc.AsB, b: 1, h: t, d: rc.d, x: cs.x, cover: rc.cover, barDiameter: rc.db, fctEff: rc.props.fctm, Ecm, Es: rc.Es, barSpacing: inp.sBar / 1e3 });
    const wDefl = inp.kSpan * Math.max(...pts.map((_, j) => eff(j, 1, both ? 1 : 0, 1).w), ...pts.map((_, j) => eff(j, 1, 1, 0).w));
    // permissible loads on the slab: bending and (for the patches) punching
    const mG = inp.kSpan * inp.gG * gk * uni[0].mx, mGy = inp.kSpan * inp.gG * gk * uni[0].my;
    const mQ1 = inp.kSpan * inp.gQ * 1e3 * Math.max(uni[0].mx, uni[0].my);
    const supLim = inp.kSup > 0 ? MRdH * inp.kSpan / inp.kSup : Infinity, Mlim = Math.min(MRd, supLim);
    const qPerm = Math.max((Mlim - Math.max(mG, mGy)) / mQ1, 0);
    let kPatch = 20;
    if (patches.length) {
      const okK = (k) => { const m = inp.kSpan * worstOf(inp.gG, both ? inp.gQ : 0, inp.gQ * k).m; const v = heavy ? Math.max(inp.gQ * k * heavy.P * heavy.dyn, 0) / (u1 * dP) : 0; return m <= Mlim && v <= vRdc; };
      let lo = 0, hi = 20; if (!okK(0)) kPatch = 0; else if (okK(hi)) kPatch = hi; else { for (let i = 0; i < 30; i++) { const mid = 0.5 * (lo + hi); if (okK(mid)) lo = mid; else hi = mid; } kPatch = lo; }
    }
    // crane rail beam: continuous beam, moving wheel group
    const rb = rcMember(inp, ctx, inp.railB, inp.railH), L = inp.railL, nSp = Math.round(inp.nSpans), nW = Math.round(inp.nWheels);
    const sec = { A: inp.railB * inp.railH, Iy: inp.railB * inp.railH ** 3 / 12, Iz: inp.railH * inp.railB ** 3 / 12, J: 0.2 * Math.min(inp.railB, inp.railH) ** 3 * Math.max(inp.railB, inp.railH), Asy: 5 / 6 * inp.railB * inp.railH, Asz: 5 / 6 * inp.railB * inp.railH };
    const mat = { E: rb.props.Ecm, nu: 0.2, rho: 0, fy: rb.fck };
    const nodes = Array.from({ length: nSp + 1 }, (_, i) => [i * L, 0, 0]);
    const members = Array.from({ length: nSp }, (_, i) => ({ a: i, b: i + 1, section: sec, material: mat }));
    const supports = nodes.map((_, i) => ({ node: i, fix: i === 0 ? [1, 1, 1, 1, 0, 1] : [0, 1, 1, 1, 0, 1] }));
    const nPos = 20 * nSp, groupLen = (nW - 1) * inp.wheelS, total = nSp * L;
    const loadCases = [{ name: 'g', member: members.map((_, i) => ({ member: i, type: 'uniform', q: [0, 0, -1], axes: 'global' })) }];
    for (let k = 0; k <= nPos; k++) {
      const x0 = -groupLen / 2 + total * k / nPos, pl = [];
      for (let wI = 0; wI < nW; wI++) { const x = x0 + wI * inp.wheelS; if (x < 0 || x > total) continue; const m = Math.min(Math.floor(x / L), nSp - 1); pl.push({ member: m, type: 'point', s: clamp(x - m * L, 0, L), P: [0, 0, -1], axes: 'global' }); }
      loadCases.push({ name: `p${k}`, member: pl });
    }
    const fr = run(() => analyzeFrame({ nodes, members, supports, loadCases }), 'Rail beam analysis');
    const nSt = 21, env = (name) => members.map((_, i) => { const d = memberDiagram(fr.results[name], i, nSt); return { My: Array.from(d.My), Vz: Array.from(d.Vz) }; });
    const gU = env('g');
    const wRail = (inp.railB * inp.railH * inp.gConc + inp.railTrib * (t * inp.gConc + inp.gSdl)) * 1e3; // permanent line load, N/m
    // unit-wheel envelopes over all group positions (frame convention for z up: My > 0 hogging)
    const sagU = members.map(() => new Float64Array(nSt)), hogU = members.map(() => new Float64Array(nSt)), shU = members.map(() => new Float64Array(nSt));
    for (let k = 0; k <= nPos; k++) { const e = env(`p${k}`); for (let i = 0; i < nSp; i++) for (let j = 0; j < nSt; j++) { sagU[i][j] = Math.max(sagU[i][j], -e[i].My[j]); hogU[i][j] = Math.max(hogU[i][j], e[i].My[j]); shU[i][j] = Math.max(shU[i][j], Math.abs(e[i].Vz[j])); } }
    const railAt = (Pw) => { let Ms = 0, Mh = 0, V = 0; const xs = [], sag = [], hog = []; for (let i = 0; i < nSp; i++) for (let j = 0; j < nSt; j++) { const s = inp.gG * wRail * -gU[i].My[j] + inp.gQ * Pw * sagU[i][j], h = inp.gG * wRail * gU[i].My[j] + inp.gQ * Pw * hogU[i][j], v = inp.gG * wRail * Math.abs(gU[i].Vz[j]) + inp.gQ * Pw * shU[i][j]; Ms = Math.max(Ms, s); Mh = Math.max(Mh, h); V = Math.max(V, v); xs.push(r(i * L + j * L / (nSt - 1))); sag.push(r(Math.max(s, 0) / 1e3)); hog.push(r(-Math.max(h, 0) / 1e3)); } return { Ms, Mh, V, xs, sag, hog }; };
    const Pw = inp.wheelP * 1e3 * inp.dynCrane, rail = railAt(Pw);
    const MRdR = momentCapacity({ section: rb.sag, material: rb.material }).MRd, MRdRH = momentCapacity({ section: rb.hog, material: rb.material }).MRd;
    const VRdR = rb.AswS > 0 ? shearResistanceLinks({ bw: rb.b, d: rb.d, fck: rb.fck, fywk: rb.fyk, AswOverS: rb.AswS, VEd: rail.V, gammaC: inp.gammaC, gammaS: inp.gammaS, alphaCC: inp.alphaCC }).VRd : shearResistanceNoLinks({ bw: rb.b, d: rb.d, fck: rb.fck, Asl: rb.AsB, gammaC: inp.gammaC, alphaCC: inp.alphaCC }).VRdc;
    const railOk = (P) => { const q = railAt(P); return q.Ms <= MRdR && q.Mh <= MRdRH && q.V <= VRdR; };
    let wheelPerm = 0; if (railOk(0)) { let lo = 0, hi = 2e7; if (railOk(hi)) wheelPerm = hi; else { for (let i = 0; i < 40; i++) { const mid = 0.5 * (lo + hi); if (railOk(mid)) lo = mid; else hi = mid; } wheelPerm = lo; } }
    wheelPerm /= inp.dynCrane;
    const checks = [
      check('Slab span moment', mEd / 1e3, MRd / 1e3, 'kN·m/m', 'EN 1992-1-1:2004 §6.1', 'Sagging moment in the panel against the resistance of the bottom bars. If exceeded the bottom steel yields under the wheel or stack and the slab sags permanently.'),
      ...(inp.kSup > 0 ? [check('Slab support moment', mSup / 1e3, MRdH / 1e3, 'kN·m/m', 'EN 1992-1-1:2004 §6.1', 'Hogging moment over the beams against the resistance of the top bars. If exceeded wide cracks open on the deck surface along the beam lines.')] : []),
      ...(heavy ? [check(`Punching under '${heavy.name}'`, vEd / 1e6, vRdc / 1e6, 'MPa', 'EN 1992-1-1:2004 §6.4.3, §6.4.4', 'Shear stress on the control perimeter two effective depths from the patch against the resistance without shear reinforcement. Punching is brittle: the wheel or corner casting pushes a cone out of the slab.')] : []),
      check('Slab crack width (quasi-permanent)', cw.wk * 1e3, inp.wLim, 'mm', 'EN 1992-1-1:2004 §7.3.4, Table 7.1N', 'Calculated crack width under the long-term share of the cargo load. Wider cracks let chlorides reach the bars faster and shorten the service life.'),
      check('Rail beam sagging', rail.Ms / 1e3, MRdR / 1e3, 'kN·m', 'EN 1992-1-1:2004 §6.1', 'Largest span moment under the moving wheel group against the bottom reinforcement of the rail beam.'),
      check('Rail beam hogging', rail.Mh / 1e3, MRdRH / 1e3, 'kN·m', 'EN 1992-1-1:2004 §6.1', 'Largest moment over the bents against the top reinforcement. It peaks with the wheel group straddling two spans.'),
      check('Rail beam shear', rail.V / 1e3, VRdR / 1e3, 'kN', 'EN 1992-1-1:2004 §6.2.3', 'Largest shear beside a bent with the wheel group next to the support, against the link resistance. Shear failure is brittle.'),
    ];
    const w = checks.reduce((m, c) => (c.uc > m.uc ? c : m)), ok = w.uc <= 1;
    const tNeed = t * Math.sqrt(Math.max(mEd / MRd, 0.01)), dConc = (tNeed - t) * a * b;
    // moment along the panel centre line for the governing case
    const xs = Array.from({ length: 17 }, (_, i) => a * (i + 0.5) / 17), fgov = gov;
    const plateP = (load, x, y) => navierPlate({ a, b, t, E: Ecm, nu: inp.nu, load, x, y, nTerms: Math.min(nT, 32) }); // the plotted line needs fewer harmonics than the design value
    const mLine = xs.map((x) => { const yy = pts[gov.j].y; let m = (inp.gG * gk + fgov.fq * qk) * plateP({ type: 'uniform', q: 1 }, x, yy).mx; for (const p of patches) m += fgov.fp * p.P * p.dyn * plateP({ type: 'patch', P: 1, x0: p.x, y0: p.y, u: p.u, v: p.v }, x, yy).mx; return r(inp.kSpan * m / 1e3); });
    return {
      kpis: [
        { key: 'mEd', label: 'Slab design span moment', value: r(mEd / 1e3, 6), unit: 'kN·m/m', limit: r(MRd / 1e3), sense: 'max', explain: `Largest factored sagging moment per metre width (${gov.name}, at ${pts[gov.j].name}), with the continuity factor, against the resistance of the bottom bars.` },
        { key: 'mField', label: 'Simply supported panel moment (before continuity)', value: r(mSS / 1e3), unit: 'kN·m/m', explain: 'Factored plate moment from the Navier solution with all edges simply supported: the reference from which span and support moments are taken.' },
        { key: 'mSup', label: 'Slab design support moment', value: r(mSup / 1e3), unit: 'kN·m/m', limit: r(MRdH / 1e3), sense: 'max', explain: 'Hogging moment over the beams taken as a fraction of the simple span moment, against the resistance of the top bars.' },
        { key: 'AsReq', label: 'Bottom reinforcement demand', value: r(AsReq * 1e6), unit: 'mm²/m', limit: r(rc.AsB * 1e6), sense: 'max', explain: 'Steel area needed for the span moment with a lever arm of 0.9 d, against the area provided. The margin is what a heavier cargo rating could use.' },
        { key: 'AsReqTop', label: 'Top reinforcement demand', value: r(AsReqTop * 1e6), unit: 'mm²/m', limit: r(rc.AsT * 1e6), sense: 'max', explain: 'Steel area needed over the beams for the support moment, against the area provided.' },
        { key: 'sigma', label: 'Slab bending stress (uncracked, characteristic)', value: r(sigma / 1e6), unit: 'MPa', explain: `Extreme-fibre stress 6m/t² under unfactored loads. Above the tensile strength of about ${f(rc.props.fctm / 1e6)} MPa the slab cracks and the reinforcement takes the tension.` },
        { key: 'vEd', label: 'Punching shear stress', value: r(vEd / 1e6), unit: 'MPa', limit: r(vRdc / 1e6), sense: 'max', explain: 'Factored shear stress on the perimeter two effective depths from the heaviest patch, against the concrete resistance without shear reinforcement.' },
        { key: 'wk', label: 'Crack width, quasi-permanent', value: r(cw.wk * 1e3), unit: 'mm', limit: inp.wLim, sense: 'max', explain: `Characteristic crack width with the bottom steel at ${f(cs.sigmaS / 1e6)} MPa. It is a durability measure for the splash zone, not a strength measure.` },
        { key: 'defl', label: 'Elastic panel deflection (characteristic)', value: r(wDefl * 1e3), unit: 'mm', explain: 'Uncracked plate deflection under unfactored loads; multiply by about 3–4 for cracking and creep when ponding or rail level matters.' },
        { key: 'railMsag', label: 'Rail beam span moment', value: r(rail.Ms / 1e3), unit: 'kN·m', limit: r(MRdR / 1e3), sense: 'max', explain: 'Envelope maximum of the factored sagging moment as the wheel group travels along the beam.' },
        { key: 'railMhog', label: 'Rail beam support moment', value: r(rail.Mh / 1e3), unit: 'kN·m', limit: r(MRdRH / 1e3), sense: 'max', explain: 'Envelope maximum of the factored hogging moment over the bents.' },
        { key: 'railV', label: 'Rail beam shear', value: r(rail.V / 1e3), unit: 'kN', limit: r(VRdR / 1e3), sense: 'max', explain: 'Envelope maximum of the factored shear force; also the largest reaction share passed to one bent divided between its two sides.' },
        { key: 'qPerm', label: 'Maximum permissible uniform cargo load', value: r(qPerm), unit: 'kPa', explain: 'Characteristic stacking load at which the slab span or support moment reaches its resistance. Post it as the deck rating for storage areas.' },
        { key: 'patchPerm', label: 'Maximum permissible wheel / corner load (heaviest patch)', value: r(heavy ? kPatch * heavy.P / 1e3 : 0), unit: 'kN', explain: 'Characteristic load on the heaviest patch, all patches scaled together, at which bending or punching reaches its limit. Reported at 20 times the entered load when not reached.' },
        { key: 'wheelPerm', label: 'Maximum permissible crane wheel load', value: r(wheelPerm / 1e3), unit: 'kN', explain: 'Characteristic wheel load of this wheel group at which the rail beam reaches its bending or shear resistance. Compare with the wheel loads of a replacement or upgraded crane.' },
      ],
      checks,
      plots: [
        { type: 'line', title: `Slab moment across the panel (${gov.name})`, x: { label: 'x across the panel', unit: 'm' }, y: { label: 'Moment mx', unit: 'kN·m/m' },
          series: [{ name: 'Design moment', x: xs.map((x) => r(x)), y: mLine }, { name: 'Resistance', x: [0, a], y: [r(MRd / 1e3), r(MRd / 1e3)], dash: true }],
          explain: 'Factored span moment along a line through the governing point. Peaks sit under the wheel patches; the curve must stay below the dashed resistance.' },
        { type: 'line', title: 'Rail beam moment envelope', x: { label: 'Distance along the beam', unit: 'm' }, y: { label: 'Moment (sagging +)', unit: 'kN·m' },
          series: [{ name: 'Sagging envelope', x: rail.xs, y: rail.sag }, { name: 'Hogging envelope', x: rail.xs, y: rail.hog }, { name: 'Sagging resistance', x: [0, r(total)], y: [r(MRdR / 1e3), r(MRdR / 1e3)], dash: true }, { name: 'Hogging resistance', x: [0, r(total)], y: [r(-MRdRH / 1e3), r(-MRdRH / 1e3)], dash: true }],
          explain: 'Largest sagging and hogging moment at each section for any position of the wheel group, with the permanent load. Both envelopes must stay inside the dashed resistances.' },
        { type: 'bar', title: 'Utilisation by check', categories: checks.map((c) => c.name), series: [{ name: 'Utilisation', y: checks.map((c) => r(c.uc, 3)) }], x: { label: 'Check', unit: '' }, y: { label: 'Utilisation', unit: '-' },
          explain: 'Demand over capacity for every slab and rail beam check. The tallest bar is the element that sets the deck rating.' },
      ],
      tables: [
        { title: 'Slab load cases', columns: [{ key: 'name', label: 'Case' }, { key: 'at', label: 'Governing point' }, { key: 'mx', label: 'mx', unit: 'kN·m/m' }, { key: 'my', label: 'my', unit: 'kN·m/m' }],
          rows: uls.map((c) => ({ name: c.name, at: pts[c.j].name, mx: r(inp.kSpan * c.e.mx / 1e3), my: r(inp.kSpan * c.e.my / 1e3) })), explain: 'Factored span moments in the two directions for each load arrangement, at the point where they are largest. Bars across the panel resist mx, bars along it my.' },
        { title: 'Permissible operational loads', columns: [{ key: 'item', label: 'Load' }, { key: 'applied', label: 'Applied' }, { key: 'perm', label: 'Permissible' }, { key: 'unit', label: 'Unit' }],
          rows: [{ item: 'Uniform cargo / container stacking', applied: inp.qCargo, perm: r(qPerm), unit: 'kPa' }, ...(heavy ? [{ item: `Heaviest patch (${heavy.name})`, applied: r(heavy.P / 1e3), perm: r(kPatch * heavy.P / 1e3), unit: 'kN' }] : []), { item: `Crane wheel (group of ${nW} at ${f(inp.wheelS)} m)`, applied: inp.wheelP, perm: r(wheelPerm / 1e3), unit: 'kN' }],
          explain: 'Characteristic loads the deck can carry at the ultimate limit state with the factors entered: the basis of the deck load plan posted for operations.' },
      ],
      notes: [`${gov.name} governs the slab at ${pts[gov.j].name}; resistances ${f(MRd / 1e3)} (span) and ${f(MRdH / 1e3)} kN·m/m (support) with ${f(inp.dBar)} mm bars at ${f(inp.sBar)} mm.`,
        `Rail beam: permanent load ${f(wRail / 1e3)} kN/m; wheel group of ${nW} × ${f(inp.wheelP)} kN (× ${f(inp.dynCrane)}) over ${f(groupLen)} m on ${nSp} span(s) of ${f(L)} m.`],
      advice: [
        ok ? eng('ok', `${w.name} governs at ${f(w.uc)}. The deck can carry ${f(qPerm)} kPa of cargo${heavy ? `, ${f(kPatch * heavy.P / 1e3)} kN on the heaviest patch` : ''} and ${f(wheelPerm / 1e3)} kN per crane wheel.`,
          `Issue the deck load plan from the permissible values and mark traffic lanes and stacking areas on the deck. A slab of about ${f(Math.max(tNeed, 0.2))} m would bring the span moment to unity; confirm punching and cracking before thinning.`, 'EN 1992-1-1:2004 §6–7; BS 6349-1-2:2016 §8',
          `Thinning the slab to ${f(Math.max(tNeed, 0.2))} m changes the panel by ${f((Math.max(tNeed, 0.2) - t) * a * b)} m³: ${concCo2((Math.max(tNeed, 0.2) - t) * a * b)}. The crack width of ${f(cw.wk * 1e3)} mm governs durability in the splash zone.`)
          : eng('fail', `${w.name}: ${f(w.demand)} against ${f(w.capacity)} ${w.unit} (${f(w.uc)}).`,
            w.name.startsWith('Punching') ? `Enlarge the bearing area (spreader plates or a thicker wearing course), add shear reinforcement, or limit the patch load to ${f(heavy ? kPatch * heavy.P / 1e3 : 0)} kN.` : w.name.startsWith('Rail') ? `Limit the crane wheel load to ${f(wheelPerm / 1e3)} kN, or strengthen the rail beam (depth, bars or links as the failed check shows).` : w.name.includes('crack') ? 'Use smaller bars at closer spacing for the same area, or more steel: crack width falls with steel stress and bar spacing.' : `Thicken the slab to about ${f(tNeed)} m or reduce the bar spacing; otherwise limit the cargo load to ${f(qPerm)} kPa.`,
            'EN 1992-1-1:2004 §6.1, §6.2, §6.4, §7.3',
            `A slab of ${f(tNeed)} m adds ${f(dConc)} m³ per panel: ${concCo2(dConc)}. An operating limit on the load costs no material.`),
        eco('info', `Ratings: ${f(qPerm)} kPa cargo against ${f(inp.qCargo)} kPa applied; ${f(wheelPerm / 1e3)} kN against ${f(inp.wheelP)} kN per crane wheel.`,
          'Unused rating is commercial capacity: a higher stack or a larger crane without civil works. Where the rating is short, compare an operating restriction (lost stacking height, slower handling) with local strengthening before committing to either.', `${ECON_BASIS}; terminal capacity planning (PIANC WG 135)`),
      ],
      publish: { 'wharf.cargoPermissible': qPerm * 1e3, 'wharf.craneWheelPermissible': wheelPerm, 'wharf.wheelPermissible': heavy ? kPatch * heavy.P : 0 },
      quantities: { concrete_m3: r(a * b * t + inp.railB * inp.railH * total) },
    };
  },
  verify: [
    // Simply supported square plate, uniform load, ν = 0.3: M = 0.0479·q·a² (Timoshenko & Woinowsky-Krieger Table 8).
    // q = 10 kPa, a = 4 m → 7.664 kN·m/m. Self-weight and surfacing removed, γ = 1.
    { name: 'Simply supported square plate under uniform load: M = 0.0479 q a²', source: 'Timoshenko & Woinowsky-Krieger, Theory of Plates and Shells, Table 8 (ν = 0.3)',
      inputs: { a: 4, b: 4, nu: 0.3, gConc: 0, gSdl: 0, qCargo: 10, patches: [], gG: 1, gQ: 1, kSpan: 1, nTerms: 96 }, expect: [{ kpi: 'mField', value: 7.664, tol: 0.003 }, { kpi: 'mEd', value: 7.664, tol: 0.003 }] },
    // One wheel at mid-span of a single simply supported span: M = P·L/4 = 1000·6/4 = 1500 kN·m; V = P = 1000 kN at the support.
    { name: 'Single wheel on a simple span: M = PL/4, V = P', source: 'Statics of a simply supported beam under a moving point load',
      inputs: { gConc: 0, gSdl: 0, nSpans: 1, railL: 6, nWheels: 1, wheelP: 1000, dynCrane: 1, gG: 1, gQ: 1 }, expect: [{ kpi: 'railMsag', value: 1500, tol: 0.002 }, { kpi: 'railV', value: 1000, tol: 0.002 }] },
    // Three equal spans under uniform load: support moment 0.100·w·L². Beam 1 × 1 m at 25 kN/m³: w = 25 kN/m, L = 6 m → 90 kN·m.
    { name: 'Three-span continuous beam under uniform load: M = 0.100 w L²', source: 'Three-moment equation (Clapeyron); standard beam tables',
      inputs: { gConc: 25, gSdl: 0, railB: 1, railH: 1, railTrib: 0, nSpans: 3, railL: 6, wheelP: 0, gG: 1, gQ: 1, nBot: 8, nTop: 8 }, expect: [{ kpi: 'railMhog', value: 90, tol: 0.01 }] },
  ],
};
