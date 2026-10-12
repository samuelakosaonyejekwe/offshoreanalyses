// Global hull-girder strength of a ship-shaped unit or barge: still-water and quasi-static wave shear force and
// bending moment from the weight and buoyancy distributions, rule wave moment, section modulus and stresses,
// oblique-wave torsion and ultimate bending capacity by progressive collapse.
import { hullGirderSection, hullGirderUltimate, stiffenerColumnBuckling } from '../../engine/structural/plates.js';
import { G, DEG, need, run, round, rnd, fmt, adv, check, sev, worst, linspace, steelCo2Text } from './_common.js';

const halfBreadth = (x, L, B, taper) => { if (!(taper > 0)) return B / 2; const a = taper * L; const d = Math.abs(x) - (L / 2 - a); return d <= 0 ? B / 2 : Math.max(0, (B / 2) * (1 - (d / a) ** 2)); };
// IACS UR S4 higher-tensile steel factor k by yield strength (235 → 1.0, 315 → 0.78, 355 → 0.72, 390 → 0.68)
const kMat = (fyMPa) => { const t = [[235, 1], [315, 0.78], [355, 0.72], [390, 0.68]]; if (fyMPa <= 235) return 1; for (let i = 1; i < t.length; i++) if (fyMPa <= t[i][0]) return t[i - 1][1] + (t[i][1] - t[i - 1][1]) * (fyMPa - t[i - 1][0]) / (t[i][0] - t[i - 1][0]); return 0.68; };
// IACS UR S11 wave coefficient C(L)
const waveC = (L) => { const l = Math.min(Math.max(L, 90), 500); return l <= 300 ? 10.75 - ((300 - l) / 100) ** 1.5 : l <= 350 ? 10.75 : 10.75 - ((l - 350) / 150) ** 1.5; };

/** Shear force and bending moment (hogging positive) of a wall-sided hull balanced on a static wave. */
function girder({ xs, b2, w, D, rho, W, xG, wave }) {
  const n = xs.length; const dx = xs[1] - xs[0];
  const integ = (f) => { let s = 0; for (let i = 1; i < n; i++) s += 0.5 * (f[i] + f[i - 1]) * dx; return s; };
  const buoy = (T0, th) => xs.map((x, i) => rho * G * b2[i] * Math.min(Math.max(T0 + th * x + (wave ? wave[i] : 0), 0), D));
  let T0 = W / (rho * G * integ(b2)); let th = 0;
  for (let it = 0; it < 40; it++) {
    const res = (a, t) => { const q = buoy(a, t); return [integ(q) - W, integ(q.map((v, i) => v * xs[i])) - W * xG]; };
    const r0 = res(T0, th); const ea = 1e-4; const et = 1e-7;
    const ra = res(T0 + ea, th); const rt = res(T0, th + et);
    const J = [[(ra[0] - r0[0]) / ea, (rt[0] - r0[0]) / et], [(ra[1] - r0[1]) / ea, (rt[1] - r0[1]) / et]];
    const det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    need(Math.abs(det) > 0, 'The weight cannot be balanced by buoyancy: check the weight table against the hull size.');
    const dT = (r0[0] * J[1][1] - r0[1] * J[0][1]) / det; const dth = (J[0][0] * r0[1] - J[1][0] * r0[0]) / det;
    T0 -= dT; th -= dth;
    if (Math.abs(dT) < 1e-9 && Math.abs(dth) < 1e-12) break;
  }
  const q = buoy(T0, th); const p = q.map((v, i) => v - w[i]);
  const V = new Array(n).fill(0); const M = new Array(n).fill(0);
  for (let i = 1; i < n; i++) V[i] = V[i - 1] + 0.5 * (p[i] + p[i - 1]) * dx;
  // closure: residual end shear from discretisation is distributed linearly
  const vEnd = V[n - 1]; for (let i = 0; i < n; i++) V[i] -= vEnd * i / (n - 1);
  for (let i = 1; i < n; i++) M[i] = M[i - 1] + 0.5 * (V[i] + V[i - 1]) * dx;
  const mEnd = M[n - 1]; for (let i = 0; i < n; i++) M[i] -= mEnd * i / (n - 1);
  return { T0, trim: th, q, V, hog: M.map((v) => -v) };
}

export default {
  id: 's11-global-strength',
  suite: 11,
  title: 'Global hull-girder strength',
  summary: 'Still-water and wave-induced shear force, bending moment and torsion along a ship-shaped hull or barge, the midship section modulus and stresses against rule allowables, and the ultimate hull-girder bending capacity.',
  covers: [7],
  fidelity: 'numerical',
  standards: ['IACS UR S11 Rev.10 (longitudinal strength standard)', 'IACS UR S4 (material factor k)', 'IACS CSR-H (2022) Pt 1 Ch 5 Sec 2 (hull-girder ultimate strength)', 'DNV-OS-C102 (2021) (ship-shaped offshore units)', 'DNV-RP-C201 (2010) (buckling of plated structures)'],
  models: ['Classical beam theory: shear force and bending moment equilibrium of weight and buoyancy', 'Hull balanced on a static sinusoidal wave (hogging and sagging)', 'Rule wave bending moment (IACS UR S11)', 'Hull-girder section properties from lumped longitudinal elements', 'Shear flow τ = F·S/(I·Σt) and Bredt closed-cell torsion', 'Stiffened-panel ultimate compressive strength', 'Smith-type progressive collapse for ultimate moment'],
  limits: 'Wall-sided hull with a parabolic plan taper; quasi-static wave balance without Smith correction or dynamic effects, so the direct wave moment is an estimate to compare with the rule value, which should be used for design. Section elements are plate strips with smeared longitudinals; no shear lag, openings, or warping torsion of open decks. The progressive-collapse analysis holds the ultimate stress after buckling (upper estimate). Hydroelastic response (springing, whipping) and finite-element global models are outside the module.',
  searchTerms: 'hull girder longitudinal strength still water wave bending moment FPSO ultimate',
  inputs: [
    { key: 'L', label: 'Length', unit: 'm', type: 'number', default: 260, min: 40, max: 500, group: 'Hull geometry', role: 'geometry', help: 'Rule length of the hull.' },
    { key: 'B', label: 'Beam', unit: 'm', type: 'number', default: 46, min: 8, max: 100, group: 'Hull geometry', role: 'geometry', help: 'Moulded breadth.' },
    { key: 'Dhull', label: 'Depth', unit: 'm', type: 'number', default: 27, min: 3, max: 60, group: 'Hull geometry', role: 'geometry', help: 'Moulded depth to the strength deck.' },
    { key: 'endTaper', label: 'Bow / stern taper length (fraction of L; 0 = box barge)', unit: '-', type: 'number', default: 0.12, min: 0, max: 0.45, group: 'Hull geometry', role: 'geometry', help: 'Each end tapers parabolically in plan over this share of the length.' },
    { key: 'weights', label: 'Weight distribution', type: 'table', group: 'Loading condition', role: 'ic',
      columns: [{ key: 'name', label: 'Item', type: 'text' }, { key: 'x1', label: 'From (aft +ve fwd of midship)', unit: 'm', type: 'number' }, { key: 'x2', label: 'To', unit: 'm', type: 'number' }, { key: 'mass', label: 'Mass', unit: 't', type: 'number' }],
      default: [{ name: 'Lightship hull', x1: -130, x2: 130, mass: 30000 }, { name: 'Machinery and accommodation', x1: -130, x2: -95, mass: 7000 }, { name: 'Topsides modules', x1: -85, x2: 105, mass: 24000 }, { name: 'Cargo oil', x1: -90, x2: 95, mass: 105000 }, { name: 'Turret and mooring', x1: 95, x2: 125, mass: 5000 }],
      help: 'Each item is spread uniformly between its two stations (metres from midship, positive forward). Enter one loading condition: full, ballast or an intermediate offloading state.' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1040, group: 'Loading condition', role: 'environment', bind: 'site.rhoW', help: 'Density of the water.' },
    { key: 'Hw', label: 'Design wave height for the direct wave moment', unit: 'm', type: 'number', default: 9.8, min: 0, max: 30, group: 'Wave loading', role: 'bc', help: 'Height of the static wave of length L on which the hull is balanced. The classic estimate 0.607·√L gives 9.8 m for 260 m; the rule moment is evaluated independently.' },
    { key: 'chiDeg', label: 'Oblique wave heading for torsion', unit: 'deg', type: 'number', default: 60, min: 10, max: 80, group: 'Wave loading', role: 'bc', help: 'Angle between the wave direction and the hull axis for the quasi-static torsional moment.' },
    { key: 'steel', label: 'Hull steel', type: 'select', source: 'materials:steel', default: 'ah36', group: 'Midship section', role: 'material', help: 'Grade of the deck and bottom plating; sets the yield stress and the rule material factor k.' },
    { key: 'tDeck', label: 'Deck plate thickness (net)', unit: 'mm', type: 'number', default: 22, min: 5, max: 80, group: 'Midship section', role: 'geometry', design: true, help: 'Strength-deck plating, net of corrosion addition.' },
    { key: 'tBottom', label: 'Bottom plate thickness (net)', unit: 'mm', type: 'number', default: 22, min: 5, max: 80, group: 'Midship section', role: 'geometry', design: true, help: 'Bottom shell plating.' },
    { key: 'tSide', label: 'Side shell thickness (net)', unit: 'mm', type: 'number', default: 18, min: 5, max: 80, group: 'Midship section', role: 'geometry', design: true, help: 'Side shell plating (each side).' },
    { key: 'tInner', label: 'Inner-bottom thickness (0 = single bottom)', unit: 'mm', type: 'number', default: 16, min: 0, max: 60, group: 'Midship section', role: 'geometry', help: 'Inner-bottom plating of a double bottom.' },
    { key: 'hDB', label: 'Double-bottom height', unit: 'm', type: 'number', default: 2.8, min: 0.5, max: 6, group: 'Midship section', role: 'geometry', help: 'Height of the inner bottom above the keel.' },
    { key: 'nLbhd', label: 'Longitudinal bulkheads', unit: '-', type: 'int', default: 2, min: 0, max: 6, group: 'Midship section', role: 'geometry', help: 'Number of continuous longitudinal bulkheads (inner hull counts as two).' },
    { key: 'tLbhd', label: 'Longitudinal bulkhead thickness', unit: 'mm', type: 'number', default: 15, min: 5, max: 60, group: 'Midship section', role: 'geometry', help: 'Plating of each longitudinal bulkhead.' },
    { key: 'sStiff', label: 'Longitudinal spacing', unit: 'mm', type: 'number', default: 850, min: 300, max: 1500, group: 'Midship section', role: 'geometry', help: 'Spacing of the longitudinal stiffeners.' },
    { key: 'lFrame', label: 'Web-frame spacing', unit: 'm', type: 'number', default: 4.5, min: 1, max: 10, group: 'Midship section', role: 'geometry', help: 'Span of the longitudinals between transverse web frames.' },
    { key: 'hw', label: 'Longitudinal web height', unit: 'mm', type: 'number', default: 420, min: 80, max: 1000, group: 'Midship section', role: 'geometry', help: 'Web of the tee or angle longitudinals.' },
    { key: 'tw', label: 'Longitudinal web thickness', unit: 'mm', type: 'number', default: 12, min: 5, max: 40, group: 'Midship section', role: 'geometry', help: 'Web thickness.' },
    { key: 'bf', label: 'Longitudinal flange width', unit: 'mm', type: 'number', default: 130, min: 0, max: 400, group: 'Midship section', role: 'geometry', help: 'Flange width (0 for flat bars).' },
    { key: 'tf', label: 'Longitudinal flange thickness', unit: 'mm', type: 'number', default: 18, min: 0, max: 50, group: 'Midship section', role: 'geometry', help: 'Flange thickness.' },
    { key: 'gammaS', label: 'Partial factor, still-water moment', unit: '-', type: 'number', default: 1.0, min: 1, max: 1.5, group: 'Acceptance criteria', role: 'criteria', help: 'γS for the ultimate check (CSR-H: 1.0).' },
    { key: 'gammaW', label: 'Partial factor, wave moment', unit: '-', type: 'number', default: 1.2, min: 1, max: 1.7, group: 'Acceptance criteria', role: 'criteria', help: 'γW for the ultimate check (CSR-H: 1.2). Site-specific offshore units use the 100-year moment with the factors of the governing offshore standard.' },
    { key: 'gammaR', label: 'Partial factor, ultimate capacity', unit: '-', type: 'number', default: 1.1, min: 1, max: 1.5, group: 'Acceptance criteria', role: 'criteria', help: 'γR = γM·γDB (CSR-H: 1.1 for sagging; 1.21 with the double-bottom factor in hogging).' },
    { key: 'nSec', label: 'Stations along the hull', unit: '-', type: 'int', default: 80, min: 20, max: 800, group: 'Numerical controls', role: 'control', help: 'Number of intervals for the load, shear and moment integration.' },
  ],
  mesh: { key: 'nSec', label: 'Stations', levels: [40, 80, 160, 320], qoi: 'MswMax' },
  compute(inp, ctx) {
    const { L, B, Dhull: D } = inp; const rho = inp.rhoW; const mat = ctx.material(inp.steel);
    need(mat && mat.fy > 0, 'Select a hull steel from the materials library.');
    const fy = mat.fy; const E = mat.E; const k = kMat(fy / 1e6);
    const n = Math.max(20, Math.round(inp.nSec)); const xs = linspace(-L / 2, L / 2, n + 1); const dx = L / n;
    const b2 = xs.map((x) => 2 * halfBreadth(x, L, B, inp.endTaper));
    const items = (inp.weights || []).filter((r) => r.mass > 0);
    need(items.length > 0, 'The weight table is empty: add the lightship and deadweight items.');
    const w = new Array(n + 1).fill(0); let W = 0; let Wx = 0;
    for (const it of items) {
      const a = Math.max(Math.min(it.x1, it.x2), -L / 2); const b = Math.min(Math.max(it.x1, it.x2), L / 2);
      need(b > a, `Weight item "${it.name}" lies outside the hull or has zero length: its stations must be within ±${fmt(L / 2)} m and differ.`);
      const q = it.mass * 1e3 * G / (b - a);
      // nodal load by exact overlap of the item with each node's tributary strip, so totals and centroids are preserved
      for (let i = 0; i <= n; i++) { const lo = Math.max(a, xs[i] - dx / 2, -L / 2); const hi = Math.min(b, xs[i] + dx / 2, L / 2); if (hi > lo) w[i] += q * (hi - lo) / (i === 0 || i === n ? dx / 2 : dx); }
      W += it.mass * 1e3 * G; Wx += it.mass * 1e3 * G * 0.5 * (a + b);
    }
    const xG = Wx / W;
    const base = { xs, b2, w, D, rho, W, xG };
    const still = girder(base);
    need(still.T0 < D && still.T0 > 0, `The loading condition gives a draft of ${fmt(still.T0)} m against a depth of ${fmt(D)} m: the hull is overloaded.`);
    const a = inp.Hw / 2; const hogW = girder({ ...base, wave: xs.map((x) => a * Math.cos(2 * Math.PI * x / L)) }); const sagW = girder({ ...base, wave: xs.map((x) => -a * Math.cos(2 * Math.PI * x / L)) });
    const absMax = (v) => v.reduce((m, x) => (Math.abs(x) > Math.abs(m) ? x : m), 0);
    const Msw = absMax(still.hog); const Fsw = absMax(still.V);
    const MwHog = Math.max(...hogW.hog.map((v, i) => v - still.hog[i])); const MwSag = Math.min(...sagW.hog.map((v, i) => v - still.hog[i]));
    const FwMax = Math.max(...hogW.V.map((v, i) => Math.abs(v - still.V[i])), ...sagW.V.map((v, i) => Math.abs(v - still.V[i])));
    const vol = W / (rho * G); const Cb = Math.max(vol / (L * B * still.T0), 0.6); const C = waveC(L);
    const MrHog = 0.19 * C * L * L * B * Cb * 1e3; const MrSag = -0.11 * C * L * L * B * (Cb + 0.7) * 1e3; // N·m
    // midship section: plate strips with longitudinals smeared into the plate area
    const As = (inp.hw * inp.tw + inp.bf * inp.tf) * 1e-6; const s = inp.sStiff / 1e3;
    const smear = (t) => t / 1e3 + As / s;
    const suOf = (t) => run(() => stiffenerColumnBuckling({ L: inp.lFrame, s, tp: t / 1e3, hw: inp.hw / 1e3, tw: inp.tw / 1e3, bf: inp.bf / 1e3, tf: inp.tf / 1e3, E, fy }), 'Stiffened-panel strength').sigmaU;
    const el = [];
    el.push({ area: B * smear(inp.tDeck), z: D, fy, sigmaU: suOf(inp.tDeck), part: 'deck' });
    el.push({ area: B * smear(inp.tBottom), z: 0, fy, sigmaU: suOf(inp.tBottom), part: 'bottom' });
    if (inp.tInner > 0) { need(inp.hDB < D, 'The double-bottom height must be less than the hull depth.'); el.push({ area: B * smear(inp.tInner), z: inp.hDB, fy, sigmaU: suOf(inp.tInner), part: 'inner bottom' }); }
    const nv = 12; const dz = D / nv;
    for (let i = 0; i < nv; i++) {
      const z = (i + 0.5) * dz;
      el.push({ area: 2 * dz * smear(inp.tSide), z, i0: 2 * smear(inp.tSide) * dz ** 3 / 12, fy, sigmaU: suOf(inp.tSide), part: 'side' });
      if (inp.nLbhd > 0) el.push({ area: inp.nLbhd * dz * smear(inp.tLbhd), z, i0: inp.nLbhd * smear(inp.tLbhd) * dz ** 3 / 12, fy, sigmaU: suOf(inp.tLbhd), part: 'bulkhead' });
    }
    const sec = hullGirderSection({ elements: el });
    const ult = run(() => hullGirderUltimate({ elements: el, E, nSteps: 120 }), 'Ultimate hull-girder strength');
    const Mhog = Math.max(...still.hog) + MrHog; const Msag = Math.min(...still.hog) + MrSag;
    const Mtot = Math.max(Math.abs(Mhog), Math.abs(Msag));
    const sigDeck = Mtot / sec.Ztop; const sigKeel = Mtot / sec.Zbottom;
    const Sq = el.filter((e) => e.z > sec.zNA).reduce((acc, e) => acc + e.area * (e.z - sec.zNA), 0);
    const tSum = (2 * inp.tSide + inp.nLbhd * inp.tLbhd) / 1e3;
    // rule wave shear is not evaluated (its distribution factors are tabulated per station); the direct quasi-static value is used
    const Ftot = Math.abs(Fsw) + FwMax; const tau = Ftot * Sq / (sec.I * tSum);
    // quasi-static torsion in an oblique wave with λ·/cos χ = L (see header)
    const chi = inp.chiDeg * DEG; const kx = 2 * Math.PI / L; const kw = kx / Math.cos(chi); const ky = kw * Math.sin(chi); const hb = B / 2;
    const Tors = rho * G * a * 2 * (Math.sin(ky * hb) - ky * hb * Math.cos(ky * hb)) / (ky * ky) * 2 / kx;
    const tauT = Math.abs(Tors) / (2 * B * D * Math.min(inp.tDeck, inp.tSide, inp.tBottom) / 1e3);
    const sagDemand = inp.gammaS * Math.abs(Math.min(...still.hog, 0)) + inp.gammaW * Math.abs(MrSag);
    const hogDemand = inp.gammaS * Math.max(...still.hog, 0) + inp.gammaW * MrHog;
    const sigAllow = 175e6 / k; const tauAllow = 110e6 / k;
    const checks = [
      check('Deck bending stress (still water + rule wave)', sigDeck / 1e6, sigAllow / 1e6, 'MPa', 'IACS UR S11 §S11.3.1 (permissible 175/k N/mm²); k per UR S4',
        'Guards against yielding and low-cycle cracking of the strength deck under the combined hull-girder moment; exceeding it means the deck area is insufficient.'),
      check('Bottom bending stress (still water + rule wave)', sigKeel / 1e6, sigAllow / 1e6, 'MPa', 'IACS UR S11 §S11.3.1', 'Same check at the keel, which also carries lateral sea pressure; exceeding it calls for more bottom area or a lower still-water moment.'),
      check('Hull-girder shear stress (still water + direct wave)', tau / 1e6, tauAllow / 1e6, 'MPa', 'IACS UR S11 §S11.4 (permissible 110/k N/mm²)', 'Guards the side shell and longitudinal bulkheads against shear yielding and shear buckling near the quarter lengths, where the shear force peaks.'),
      check('Ultimate strength, sagging', sagDemand / 1e9, ult.sagging.Mu / inp.gammaR / 1e9, 'GN·m', 'IACS CSR-H (2022) Pt 1 Ch 5 Sec 2 [2] (γS·Msw + γW·Mw ≤ MU/γR)', 'Guards against collapse of the hull girder by deck buckling in sagging: the failure that breaks a ship\'s back.'),
      check('Ultimate strength, hogging', hogDemand / 1e9, ult.hogging.Mu / inp.gammaR / 1e9, 'GN·m', 'IACS CSR-H (2022) Pt 1 Ch 5 Sec 2 [2]', 'Guards against collapse by bottom buckling in hogging.'),
      check('Torsional shear stress (oblique wave)', tauT / 1e6, tauAllow / 1e6, 'MPa', 'Bredt closed-cell torsion; DNV-OS-C102 (2021) global strength', 'Guards against shear overstress from twisting in quartering seas; significant only for hulls with large deck openings, where warping stresses (not modelled) add.'),
    ];
    const gov = worst(checks);
    const Zreq = Mtot / sigAllow; const steelPerM = sec.A * 7.85; // t per metre of hull, longitudinal material
    const dT = Math.max(0, (Zreq / sec.Ztop - 1)) * inp.tDeck;
    const xsR = rnd(xs);
    return {
      kpis: [
        { key: 'draft', label: 'Mean draft of the loading condition', value: round(still.T0, 5), unit: 'm', explain: 'Draft amidships at which buoyancy balances the weight table.' },
        { key: 'trim', label: 'Trim (bow down positive)', value: round(still.trim * L, 4), unit: 'm', explain: 'Difference of draft between bow and stern from the longitudinal position of the weights.' },
        { key: 'MswMax', label: 'Still-water bending moment (hogging +)', value: round(Msw / 1e6, 5), unit: 'MN·m', explain: 'Largest still-water moment along the hull. Hogging (positive) puts the deck in tension; a full cargo amidships gives sagging.' },
        { key: 'FswMax', label: 'Still-water shear force', value: round(Math.abs(Fsw) / 1e6, 5), unit: 'MN', explain: 'Largest still-water shear force, normally near the quarter lengths or at a tank boundary.' },
        { key: 'MwHog', label: 'Direct wave moment, hogging', value: round(MwHog / 1e6, 5), unit: 'MN·m', explain: 'Extra moment when the hull is balanced on a wave of its own length with the crest amidships (quasi-static, no Smith correction).' },
        { key: 'MwSag', label: 'Direct wave moment, sagging', value: round(MwSag / 1e6, 5), unit: 'MN·m', explain: 'Extra moment with the trough amidships (negative = sagging).' },
        { key: 'MwRuleHog', label: 'Rule wave moment, hogging', value: round(MrHog / 1e6, 5), unit: 'MN·m', explain: 'IACS UR S11 vertical wave bending moment amidships for unrestricted (North Atlantic) service, 0.19·C·L²·B·Cb.' },
        { key: 'MwRuleSag', label: 'Rule wave moment, sagging', value: round(MrSag / 1e6, 5), unit: 'MN·m', explain: 'IACS UR S11 sagging wave moment, −0.11·C·L²·B·(Cb + 0.7).' },
        { key: 'Zdeck', label: 'Section modulus at deck', value: round(sec.Ztop, 5), unit: 'm³', explain: 'Hull-girder section modulus at the strength deck from the plating and smeared longitudinals.' },
        { key: 'Zkeel', label: 'Section modulus at keel', value: round(sec.Zbottom, 5), unit: 'm³', explain: 'Hull-girder section modulus at the bottom.' },
        { key: 'sigDeck', label: 'Deck stress', value: round(sigDeck / 1e6, 4), unit: 'MPa', limit: round(sigAllow / 1e6, 4), sense: 'max', explain: 'Nominal hull-girder bending stress at the deck under still-water plus rule wave moment.' },
        { key: 'tau', label: 'Hull-girder shear stress', value: round(tau / 1e6, 4), unit: 'MPa', limit: round(tauAllow / 1e6, 4), sense: 'max', explain: 'Shear stress at the neutral axis from the largest total shear force, shared by side shell and longitudinal bulkheads.' },
        { key: 'MuSag', label: 'Ultimate moment, sagging', value: round(ult.sagging.Mu / 1e9, 5), unit: 'GN·m', explain: 'Peak of the moment–curvature curve with the deck in compression (stiffened-panel buckling limits the deck stress).' },
        { key: 'MuHog', label: 'Ultimate moment, hogging', value: round(ult.hogging.Mu / 1e9, 5), unit: 'GN·m', explain: 'Peak of the moment–curvature curve with the bottom in compression.' },
        { key: 'torsion', label: 'Wave torsional moment (oblique wave)', value: round(Math.abs(Tors) / 1e6, 5), unit: 'MN·m', explain: 'Quasi-static twisting moment amidships in a quartering wave whose projected length equals the hull length.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Weight and buoyancy distribution (still water)', x: { label: 'Station from midship', unit: 'm' }, y: { label: 'Load', unit: 'MN/m' },
          series: [{ name: 'Weight', x: xsR, y: rnd(w.map((v) => v / 1e6)) }, { name: 'Buoyancy', x: xsR, y: rnd(still.q.map((v) => v / 1e6)) }],
          explain: 'The two distributions whose difference loads the hull girder. Where weight exceeds buoyancy the hull is pushed down; large local differences produce shear and moment.' },
        { type: 'line', title: 'Shear force', x: { label: 'Station from midship', unit: 'm' }, y: { label: 'Shear force', unit: 'MN' },
          series: [{ name: 'Still water', x: xsR, y: rnd(still.V.map((v) => v / 1e6)) }, { name: 'Wave crest amidships', x: xsR, y: rnd(hogW.V.map((v) => v / 1e6)), dash: true }, { name: 'Wave trough amidships', x: xsR, y: rnd(sagW.V.map((v) => v / 1e6)), dash: true }],
          explain: 'Shear force along the hull; it must close to zero at both ends. Peaks near the quarter lengths are checked against the shear capacity of the side shell and bulkheads.' },
        { type: 'line', title: 'Bending moment (hogging positive)', x: { label: 'Station from midship', unit: 'm' }, y: { label: 'Bending moment', unit: 'MN·m' },
          series: [{ name: 'Still water', x: xsR, y: rnd(still.hog.map((v) => v / 1e6)) }, { name: 'Total, crest amidships', x: xsR, y: rnd(hogW.hog.map((v) => v / 1e6)), dash: true }, { name: 'Total, trough amidships', x: xsR, y: rnd(sagW.hog.map((v) => v / 1e6)), dash: true }],
          explain: 'Still-water moment and the envelope on the hogging and sagging wave. The largest ordinate, normally within 0.4 L amidships, drives the deck and bottom scantlings.' },
        { type: 'line', title: 'Moment–curvature (progressive collapse)', x: { label: 'Curvature', unit: '1/km' }, y: { label: 'Moment', unit: 'GN·m' },
          series: [{ name: 'Sagging', x: rnd(ult.sagging.curvature.map((c) => Math.abs(c) * 1e3)), y: rnd(ult.sagging.moment.map((m) => Math.abs(m) / 1e9)) }, { name: 'Hogging', x: rnd(ult.hogging.curvature.map((c) => Math.abs(c) * 1e3)), y: rnd(ult.hogging.moment.map((m) => Math.abs(m) / 1e9)) },
            { name: 'Design sagging demand', x: [0, Math.abs(ult.sagging.curvature[ult.sagging.curvature.length - 1]) * 1e3], y: [sagDemand / 1e9, sagDemand / 1e9], dash: true }],
          explain: 'Bending capacity as curvature grows. The plateau is the ultimate moment; the dashed line is the factored sagging demand, which must stay below the capacity divided by γR.' },
      ],
      tables: [
        { title: 'Moments and shear forces', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['Total weight', W / G / 1e3, 't'], ['LCG from midship', xG, 'm'], ['Block coefficient (rule minimum 0.6)', Cb, '-'], ['Wave coefficient C', C, '-'], ['Still-water moment (max hog)', Math.max(...still.hog) / 1e6, 'MN·m'], ['Still-water moment (max sag)', Math.min(...still.hog) / 1e6, 'MN·m'],
            ['Rule wave moment hog', MrHog / 1e6, 'MN·m'], ['Rule wave moment sag', MrSag / 1e6, 'MN·m'], ['Direct wave moment hog', MwHog / 1e6, 'MN·m'], ['Direct wave moment sag', MwSag / 1e6, 'MN·m'], ['Still-water shear', Math.abs(Fsw) / 1e6, 'MN'], ['Direct wave shear', FwMax / 1e6, 'MN'], ['Torsional moment', Math.abs(Tors) / 1e6, 'MN·m'],
            ['Neutral axis above keel', sec.zNA, 'm'], ['Moment of inertia', sec.I, 'm⁴'], ['Required deck section modulus', Zreq, 'm³'], ['Fully plastic moment', ult.Mp / 1e9, 'GN·m']].map(([q, v, u]) => ({ q, v: round(v, 5), u })),
          explain: 'Summary of the global loads and section properties in the layout of a longitudinal-strength report.' },
      ],
      notes: [
        `Loading condition: ${fmt(W / G / 1e3, 5)} t at ${fmt(still.T0)} m draft, trim ${fmt(still.trim * L)} m; still-water moment ${fmt(Msw / 1e6)} MN·m (${Msw >= 0 ? 'hogging' : 'sagging'}).`,
        `Rule wave moments ${fmt(MrHog / 1e6)} / ${fmt(MrSag / 1e6)} MN·m against the direct quasi-static values ${fmt(MwHog / 1e6)} / ${fmt(MwSag / 1e6)} MN·m for H = ${fmt(inp.Hw)} m.`,
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}); material factor k = ${fmt(k)} for fy = ${fmt(fy / 1e6)} MPa.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at utilisation ${fmt(gov.uc)}; deck stress ${fmt(sigDeck / 1e6)} MPa against ${fmt(sigAllow / 1e6)} MPa, ultimate sagging margin ${fmt(ult.sagging.Mu / inp.gammaR / Math.max(sagDemand, 1))}.`,
          gov.uc > 1 ? `Increase the deck section modulus to at least ${fmt(Zreq)} m³ (about ${fmt(dT)} mm more deck plate or heavier deck longitudinals), or cut the still-water moment by redistributing cargo and ballast: it is ${fmt(100 * Math.abs(Msw) / Math.max(Mtot, 1))} % of the total.`
            : `Reserve is ${fmt(100 * (1 - gov.uc))} %. Check every loading and offloading sequence (the intermediate states often give the largest still-water moment) and replace the rule wave moment by the site-specific 100-year value for a permanently moored unit.`,
          'IACS UR S11 Rev.10; IACS CSR-H (2022) Pt 1 Ch 5; DNV-OS-C102 (2021) for site-specific wave loads on ship-shaped units.'),
        adv('economic', 'info', `Longitudinal material is ${fmt(steelPerM)} t per metre of hull (about ${fmt(steelPerM * 0.6 * L)} t over the 0.6 L mid-body). A loading manual that limits the still-water moment is free; steel is not.`,
          `Set permissible still-water moments in the loading computer at the value used here (${fmt(Math.abs(Msw) / 1e6)} MN·m) plus a margin, and optimise the offloading sequence before adding scantlings.`,
          'Whole-life costing per ISO 15686-5; AACE RP 18R-97 (hull steel weight drives class 4 estimates).',
          `Each millimetre of deck plate over the mid-body adds ${fmt(B * 0.6 * L * 0.001 * 7.85)} t: ${steelCo2Text(B * 0.6 * L * 0.001 * 7.85)} Higher-tensile steel reaches the same utilisation with less mass but lowers fatigue and buckling margins.`),
      ],
    };
  },
  verify: [
    { // Box barge 100 × 20 m, weight 10 000 t spread over the middle half (−25…25 m): buoyancy W/L uniform, weight 2W/L.
      // Midship moment = (W/2)·(L/4) − (W/2)·(L/8) = W·L/16 = 98.0665e6·100/16 = 612.9 MN·m, sagging (negative).
      // Shear at x = −25 m: buoyancy of the aft quarter = W/4 = 24.52 MN (the load step is smeared over one station, hence 1.5 %).
      name: 'Box barge with the weight over the middle half: M = W·L/16 sagging, F = W/4',
      source: 'Beam statics, hand calculation (weight and buoyancy as uniform loads)',
      inputs: { L: 100, B: 20, Dhull: 12, endTaper: 0, Hw: 0, weights: [{ name: 'cargo', x1: -25, x2: 25, mass: 10000 }], nSec: 160 },
      expect: [{ kpi: 'MswMax', value: -612.9, tol: 0.003 }, { kpi: 'FswMax', value: 24.517, tol: 0.015 }, { kpi: 'draft', value: 4.878, tol: 0.001 }] },
    { // Uniformly loaded box barge on a wave of its own length, H = 4 m, crest amidships: extra buoyancy ρ g B a cos(2πx/L);
      // midship hogging moment = ρ g B a L²/(2π²) = 1025·9.80665·20·2·10⁴/19.739 = 203.7 MN·m; sagging −203.7 MN·m.
      name: 'Box barge on a sinusoidal wave λ = L: Mw = ρ g B a L²/(2π²)',
      source: 'Static balance on a wave, closed form for a rectangular waterplane (hand calculation)',
      inputs: { L: 100, B: 20, Dhull: 12, endTaper: 0, Hw: 4, weights: [{ name: 'uniform', x1: -50, x2: 50, mass: 10000 }], nSec: 160 },
      expect: [{ kpi: 'MwHog', value: 203.7, tol: 0.003 }, { kpi: 'MwSag', value: -203.7, tol: 0.003 }, { kpi: 'MswMax', value: 0, abs: 0.5 }] },
    { // L = 100 m: C = 10.75 − ((300 − 100)/100)^1.5 = 10.75 − 2.8284 = 7.9216; M = 0.19·7.9216·100²·20·1.0 = 301 020 kN·m.
      name: 'IACS wave coefficient and rule hogging moment for a 100 m barge (Cb = 1)',
      source: 'IACS UR S11 Rev.10 §S11.2.2.1 — formula evaluated by hand',
      inputs: { L: 100, B: 20, Dhull: 12, endTaper: 0, weights: [{ name: 'uniform', x1: -50, x2: 50, mass: 10000 }] },
      expect: [{ kpi: 'MwRuleHog', value: 301.02, tol: 0.001 }] },
  ],
};
