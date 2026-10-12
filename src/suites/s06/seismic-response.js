// Seismic and extreme loading of the dolphin: modal analysis of the equivalent-fixity frame with cap mass and
// hydrodynamic added mass, response-spectrum analysis, pile checks, drift and P–Δ sensitivity.
import { modalAnalysis } from '../../engine/structural/framemodal.js';
import { ec8Spectrum, responseSpectrumAnalysis } from '../../engine/structural/framedynamics.js';
import { isoTubularCheck } from '../../engine/structural/tubular.js';
import { G, DEG, need, run, fmt, round, cap, adv, check, ucCheck, linspace, pileInputs, soilInputs, layoutInputs, capInputs, frameInputs, fixityLength, pileProps, profileFrom, pilesFrom, frameModel, framePlot, govUc, pileSteelT, steelCo2Text, CONFIRM, NOT_HIFI, ECON_BASIS } from './_common.js';

export default {
  id: 's06-seismic-response', suite: 6, title: 'Seismic and extreme loading response',
  summary: 'Natural periods, response-spectrum base shear, cap displacement, pile forces and member utilisation of the dolphin under earthquake ground motion.',
  covers: [11], fidelity: 'numerical',
  standards: ['EN 1998-1:2004 §3.2.2 (response spectra), §4.3.3.3 (modal response spectrum analysis), §4.4.2.2 (P–Δ)', 'EN 1998-5:2004 §5 (foundations)', 'ISO 19901-2:2017 (seismic design procedures)', 'PIANC WG34 Seismic Design Guidelines for Port Structures (2001)', 'ISO 19902:2020 §13'],
  models: ['Modal vibration analysis (consistent mass, subspace iteration)', 'Seismic response spectrum analysis with CQC modal combination', 'Hydrodynamic added mass and contained water', 'Equivalent depth of fixity from p–y pile-head stiffness', 'Rayleigh quotient check of the fundamental period', 'P–Δ stability coefficient'],
  limits: 'Linear response-spectrum analysis on equivalent-fixity piles: soil nonlinearity enters only through the secant stiffness and the behaviour factor. Kinematic pile bending from soil displacement, liquefaction and lateral spreading are not modelled — screen them with the geotechnical suite. Nonlinear time-history analysis and vessel-moored seismic cases are outside this module. ' + NOT_HIFI,
  publishes: ['dolphin.f1', 'dolphin.seismicShear', 'dolphin.seismicDisp'],
  searchTerms: 'seismic design pile supported dolphin response spectrum wharf',
  inputs: [
    ...layoutInputs(),
    ...pileInputs(),
    ...capInputs(),
    ...soilInputs(),
    ...frameInputs(),
    { key: 'Ca', label: 'Added-mass coefficient of the submerged piles', unit: '-', type: 'number', default: 1.0, min: 0, max: 2, group: 'Frame model', role: 'load', help: '1.0 for a circular cylinder oscillating in still water; 0 to ignore the water.' },
    { key: 'flooded', label: 'Piles filled with water or soil plug', type: 'bool', default: true, group: 'Frame model', role: 'ic' },
    { key: 'nModes', label: 'Modes included', unit: '-', type: 'int', default: 6, min: 1, max: 20, group: 'Frame model', role: 'control' },
    { key: 'ag', label: 'Design ground acceleration on rock', unit: 'm/s²', type: 'number', default: 2.0, min: 0.05, max: 10, group: 'Seismic action', role: 'bc', bind: 'site.pga475', uncertain: { dist: 'lognormal', cov: 0.3 }, help: '475-year reference peak ground acceleration times the importance factor.' },
    { key: 'ground', label: 'Ground type', type: 'select', default: 'C', group: 'Seismic action', role: 'soil', options: ['A', 'B', 'C', 'D', 'E'].map((v) => ({ value: v, label: `EN 1998-1 ground type ${v}` })) },
    { key: 'specType', label: 'Spectrum type', type: 'select', default: '1', group: 'Seismic action', role: 'bc', options: [{ value: '1', label: 'Type 1 (surface-wave magnitude above 5.5)' }, { value: '2', label: 'Type 2 (magnitude 5.5 or less)' }] },
    { key: 'zeta', label: 'Damping ratio', unit: '%', type: 'number', default: 5, min: 0.5, max: 20, group: 'Seismic action', role: 'material', help: 'Structural plus soil radiation damping; 2–5 % for steel pile dolphins.' },
    { key: 'q', label: 'Behaviour factor q', unit: '-', type: 'number', default: 1.5, min: 1, max: 4, group: 'Seismic action', role: 'criteria', help: `1.0 for elastic response; 1.5 limited ductility; up to 2–3 for vertical steel piles designed for ductility (EN 1998, PIANC WG34). Displacements are multiplied back by q. ${CONFIRM}` },
    { key: 'dirRule', label: 'Directional combination', type: 'select', default: '100-30-30', group: 'Seismic action', role: 'control', options: [{ value: '100-30-30', label: '100 % + 30 % rule' }, { value: 'SRSS', label: 'SRSS of the two horizontal directions' }] },
    { key: 'uLimit', label: 'Allowable cap displacement', unit: 'mm', type: 'number', default: 400, min: 5, max: 3000, group: 'Acceptance', role: 'criteria', help: 'Set by walkway bearings, loading arms and pipework crossing to the dolphin.' },
    { key: 'Hcap', label: 'Lateral capacity of the dolphin', unit: 'kN', type: 'number', default: 4000, min: 1, max: 1e6, group: 'Acceptance', role: 'criteria', bind: 'dolphin.lateralCapacity' },
  ],
  mesh: { key: 'nSub', label: 'Elements per pile', levels: [1, 2, 4, 8], qoi: 'T1' },
  compute(inp, ctx) {
    const P = pileProps(inp, ctx);
    const piles = pilesFrom(inp), n = piles.length;
    const fixed = inp.headFix === 'fixed' && n > 1;
    const profile = inp.LeUser > 0 ? null : profileFrom(inp, P.D);
    const Le = fixityLength(inp, P, profile, fixed);
    // SI: kN → N → kg
    const W = (inp.capB * inp.capL * inp.capT * inp.gammaConc + inp.Wequip) * 1e3, mCap = W / G;
    need(mCap > 0, 'The cap has no mass: give its dimensions or the equipment weight.');
    const Di = P.D - 2 * P.t;
    // water moving with the pile (added mass) and inside it, per metre; applied over the whole equivalent column
    const extraMass = 1025 * Math.PI / 4 * (inp.Ca * P.D * P.D + (inp.flooded ? Di * Di : 0));
    const fm = frameModel({ piles, P, Le, nSub: inp.nSub, capMass: mCap, fixedHead: fixed, extraMass });
    const opt = { shearDeformation: !!inp.shear };
    const nFree = fm.model.nodes.length * 6 - n * 6;
    const modal = run(() => modalAnalysis(fm.model, { nModes: Math.max(1, Math.min(inp.nModes, Math.floor(nFree / 2))), ...opt }), 'Modal analysis');
    const zeta = inp.zeta / 100;
    const sp = run(() => ec8Spectrum({ ag: inp.ag, groundType: inp.ground, type: Number(inp.specType), zeta, ...(inp.q > 1 ? { q: inp.q } : {}) }), 'Response spectrum');
    const rsa = run(() => responseSpectrumAnalysis(modal, { spectrum: sp.Sa, directions: ['x', 'y'], zeta, combination: 'CQC', directional: inp.dirRule }), 'Response-spectrum analysis');
    const c = fm.centre, qd = Math.max(inp.q, 1);
    const ux = rsa.displacements[6 * c] * qd, uy = rsa.displacements[6 * c + 1] * qd;
    const u = Math.max(Math.abs(ux), Math.abs(uy));
    const rot = Math.hypot(rsa.displacements[6 * c + 3], rsa.displacements[6 * c + 4]) * qd;
    const T1 = modal.periods[0], Sa1 = sp.Sa(T1);
    const Vb = Math.max(rsa.baseShear.x, rsa.baseShear.y);
    const mTot = modal.totalMass.x;
    const rows = piles.map((p, k) => {
      let N = 0, M = 0, V = 0;
      for (let e = 0; e < inp.nSub; e++) { const f = rsa.memberForces[k * inp.nSub + e]; N = Math.max(N, f[0], f[6]); V = Math.max(V, Math.hypot(f[1], f[2]), Math.hypot(f[7], f[8])); M = Math.max(M, Math.hypot(f[4], f[5]), Math.hypot(f[10], f[11])); }
      const Ng = W / n;                                   // gravity share, compression
      const chk = isoTubularCheck({ D: P.D, t: P.t, fy: P.fy, E: P.E, L: Le, Ky: fixed ? 1 : 2, Kz: fixed ? 1 : 2, N: -(Ng + N), My: M, V });
      const g = govUc(chk);
      return { id: k + 1, Ns: N, Nc: Ng + N, Nt: Math.max(N - Ng, 0), M, V, uc: g.uc, gov: g.name };
    });
    const gov = rows.reduce((a, r) => (r.uc > a.uc ? r : a), rows[0]);
    // P–Δ sensitivity: θ = W·d / (V·h) with d the inelastic displacement (EN 1998-1 §4.4.2.2)
    const theta = W * u / (Math.max(Vb, 1) * Le);
    const part = Math.max(rsa.massParticipation.x, rsa.massParticipation.y);
    const checks = [
      ucCheck(`Tubular member check under gravity + earthquake, pile ${gov.id} (${gov.gov || 'combined'})`, gov.uc, 'ISO 19902:2020 §13.3; EN 1998-1 §4.4.2', 'Axial force and bending from the seismic combination against the member strength and stability expressions. Above 1 the pile yields; with q above 1 it must then be detailed for ductility.'),
      check('Seismic base shear against lateral capacity', Vb / 1e3, inp.Hcap, 'kN', 'EN 1998-5:2004 §5.4; PIANC WG34 (2001)', 'Design base shear against the ultimate horizontal resistance of piles and soil. If exceeded a plastic mechanism forms in the design earthquake.'),
      check('Cap displacement', u * 1e3, inp.uLimit, 'mm', 'EN 1998-1:2004 §4.3.4 (ds = q·de); project limit', 'Peak displacement of the dolphin head (elastic value times q) against what walkways, pipework and loading arms tolerate. Exceedance pulls walkways off their seats.'),
      check('P–Δ sensitivity coefficient θ', theta, 0.3, '-', 'EN 1998-1:2004 §4.4.2.2 (θ ≤ 0.1 negligible, ≤ 0.2 amplify by 1/(1 − θ), never above 0.3)', 'Overturning from the cap weight acting on the displaced dolphin, relative to the seismic shear. Above 0.3 the structure is too flexible to be stable in the earthquake.'),
    ];
    const w = checks.reduce((a, q) => (q.uc > a.uc ? q : a)), ok = w.uc <= 1;
    const Ts = linspace(0.02, Math.max(4, 1.5 * T1), 60);
    const mode = modal.modes[0];
    let mx = 0; for (let i = 0; i < fm.model.nodes.length; i++) mx = Math.max(mx, Math.hypot(mode[6 * i], mode[6 * i + 1], mode[6 * i + 2]));
    const steelT = pileSteelT(inp, n);
    const plateau = T1 >= sp.TB && T1 <= sp.TC;
    const advice = [
      adv('engineering', ok ? 'ok' : 'fail', `Fundamental period ${fmt(T1)} s (spectral acceleration ${fmt(Sa1)} m/s²); base shear ${fmt(Vb / 1e3)} kN = ${fmt(Vb / (mTot * G) * 100)} % of the seismic weight; cap displacement ${fmt(u * 1e3)} mm; member utilisation ${fmt(gov.uc)}.`,
        ok ? `${plateau ? 'The dolphin sits on the spectral plateau: any stiffening will not reduce the force.' : T1 > sp.TC ? 'The dolphin is on the descending branch: flexibility is protecting it, displacement governs.' : 'The dolphin is stiffer than the plateau.'} Confirm liquefaction and lateral-spreading potential of the upper soils and add the kinematic pile moments before finalising; check ${fmt(part * 100)} % mass participation is sufficient (≥ 90 %).`
          : (w.name.startsWith('Cap') ? `Reduce the displacement by stiffening (diameter, rake) or accept it by articulating walkways and pipework for ±${fmt(u * 1e3)} mm.` : w.name.startsWith('P–Δ') ? 'Reduce cap mass or stiffen the piles: the dolphin is too flexible for its top weight.' : `Reduce the cap mass (${fmt(mCap / 1e3)} t), which drives the inertia force, or detail the pile heads for ductility and justify a higher behaviour factor.`),
        'EN 1998-1:2004 §3.2.2, §4.3.3.3, §4.4.2; EN 1998-5:2004; PIANC WG34 (2001); ISO 19901-2:2017',
        `Seismic force scales with the cap mass: a steel deck instead of the ${fmt(inp.capB * inp.capL * inp.capT)} m³ concrete cap would cut both the inertia load and about ${fmt(inp.capB * inp.capL * inp.capT * 0.3)} t CO₂e of concrete (screening factor 0.3 t/m³), at the price of a steel deck of roughly 0.25 t/m²: ${steelCo2Text(inp.capB * inp.capL * 0.25)}`),
      adv('economic', ok ? 'info' : 'warn', `Seismic design adds ${fmt(Vb / 1e3)} kN of lateral demand; the berthing load is usually larger for breasting dolphins, so earthquake governs mainly mooring dolphins with heavy caps.`,
        'Decide the performance level per earthquake return period (serviceable after the 75-year event, repairable after the 475-year event) with the terminal operator: it sets q and the repair cost and downtime carried as risk rather than as steel.',
        `Performance-based seismic design per PIANC WG34; ${ECON_BASIS}`),
    ];
    return {
      kpis: [
        { key: 'T1', label: 'Fundamental period', unit: 's', value: T1, explain: 'Period of the first sway mode with cap mass, pile mass, contained and added water. It places the dolphin on the response spectrum; long periods mean low force but large displacement.' },
        { key: 'f1', label: 'Fundamental frequency', unit: 'Hz', value: modal.frequencies[0], explain: 'Inverse of the fundamental period. Keep it away from wave periods and from vessel surge / sway periods on the mooring.' },
        { key: 'Sa1', label: 'Design spectral acceleration at T1', unit: 'm/s²', value: Sa1, explain: 'Ordinate of the design spectrum at the fundamental period (already divided by the behaviour factor when q > 1).' },
        { key: 'Vb', label: 'Seismic base shear', unit: 'kN', value: Vb / 1e3, limit: inp.Hcap, sense: 'max', explain: 'Combined (CQC) horizontal shear at the level of fixity in the more loaded direction: the earthquake design load of the foundation.' },
        { key: 'uCap', label: 'Cap displacement (q × elastic)', unit: 'mm', value: u * 1e3, limit: inp.uLimit, sense: 'max', explain: 'Peak horizontal movement of the pile cap in the design earthquake, including the inelastic part through q.' },
        { key: 'rotCap', label: 'Cap rotation', unit: 'deg', value: rot / DEG, explain: 'Peak tilt of the cap (pile-head rotation for a mono-pile), with q applied.' },
        { key: 'Mmax', label: 'Largest seismic pile moment', unit: 'kN·m', value: Math.max(...rows.map((r) => r.M)) / 1e3, explain: 'Peak bending moment in any pile from the response-spectrum analysis, normally at the cap for fixed heads.' },
        { key: 'NsMax', label: 'Largest seismic axial force', unit: 'kN', value: Math.max(...rows.map((r) => r.Ns)) / 1e3, explain: 'Push–pull axial force from overturning, to be added to and subtracted from the gravity load.' },
        { key: 'NtMax', label: 'Largest net pile tension', unit: 'kN', value: Math.max(...rows.map((r) => r.Nt)) / 1e3, explain: 'Seismic axial force less the gravity share: positive means uplift on the pile, to be checked against pull-out capacity.' },
        { key: 'uc', label: 'Governing member utilisation', unit: '-', value: gov.uc, limit: 1, sense: 'max', explain: 'ISO 19902 utilisation of the most loaded pile under gravity plus earthquake.' },
        { key: 'theta', label: 'P–Δ sensitivity coefficient', unit: '-', value: theta, limit: 0.3, sense: 'max', explain: 'Second-order overturning over first-order overturning. Below 0.1 P–Δ can be ignored; between 0.1 and 0.2 amplify seismic effects by 1/(1 − θ).' },
        { key: 'part', label: 'Mass participation of the modes used', unit: '%', value: part * 100, explain: 'Share of the total mass captured by the modes included, in the better-captured horizontal direction. Add modes if it is below 90 %.' },
        { key: 'mass', label: 'Seismic mass', unit: 't', value: mTot / 1e3, explain: 'Total mass moving horizontally: cap and equipment, piles to the level of fixity, contained and added water.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Design response spectrum and modal periods', x: { label: 'Period', unit: 's' }, y: { label: 'Spectral acceleration', unit: 'm/s²' },
          series: [{ name: `Spectrum (ground ${inp.ground}, q = ${fmt(qd)})`, x: Ts.map((t) => round(t)), y: Ts.map((t) => round(sp.Sa(t))) }, { name: 'Modes', x: modal.periods.map((t) => round(t)), y: modal.periods.map((t) => round(sp.Sa(t))) }],
          explain: 'The design spectrum with a marker at each natural period. The position of the first mode relative to the plateau (between TB and TC) tells whether stiffening raises or lowers the seismic force.' },
        framePlot(fm.model, { disp: mode, scale: 0.15 * Le / (mx || 1), title: 'First mode shape', explain: 'Equivalent-fixity frame with the fundamental mode overlaid (arbitrary amplitude): a sway of the cap. A twisting first mode would indicate an eccentric mass or unbalanced pile layout.' }),
        { type: 'bar', title: 'Effective modal mass', categories: modal.periods.map((t, i) => `Mode ${i + 1} (${fmt(t)} s)`), series: [{ name: 'x', y: modal.massFraction.map((m) => round(m.x * 100)) }, { name: 'y', y: modal.massFraction.map((m) => round(m.y * 100)) }], y: { label: 'Share of total mass', unit: '%' },
          explain: 'Fraction of the mass responding in each mode and direction. One dominant bar per direction means the dolphin behaves as a single-degree-of-freedom oscillator.' },
      ],
      tables: [
        { title: 'Pile forces in the seismic combination', columns: [{ key: 'id', label: 'Pile' }, { key: 'Ns', label: 'Seismic axial ±', unit: 'kN' }, { key: 'Nc', label: 'Max compression', unit: 'kN' }, { key: 'Nt', label: 'Net tension', unit: 'kN' }, { key: 'V', label: 'Shear', unit: 'kN' }, { key: 'M', label: 'Moment', unit: 'kN·m' }, { key: 'uc', label: 'Utilisation' }],
          rows: rows.map((r) => ({ id: r.id, Ns: round(r.Ns / 1e3), Nc: round(r.Nc / 1e3), Nt: round(r.Nt / 1e3), V: round(r.V / 1e3), M: round(r.M / 1e3), uc: round(r.uc) })),
          explain: 'One row per pile: CQC-combined seismic forces, their combination with the gravity share and the member utilisation.' },
        { title: 'Modes', columns: [{ key: 'i', label: 'Mode' }, { key: 'T', label: 'Period', unit: 's' }, { key: 'f', label: 'Frequency', unit: 'Hz' }, { key: 'Sa', label: 'Sa', unit: 'm/s²' }, { key: 'mx', label: 'Mass x', unit: '%' }, { key: 'my', label: 'Mass y', unit: '%' }],
          rows: modal.periods.map((t, i) => ({ i: i + 1, T: round(t), f: round(modal.frequencies[i]), Sa: round(sp.Sa(t)), mx: round(modal.massFraction[i].x * 100), my: round(modal.massFraction[i].y * 100) })),
          explain: 'Natural periods with their spectral ordinate and effective mass, for checking the modal truncation.' },
      ],
      notes: [`${w.name} governs (utilisation ${fmt(w.uc)}).`, `Length to fixity ${fmt(Le)} m; spectrum corner periods TB = ${fmt(sp.TB)} s, TC = ${fmt(sp.TC)} s, TD = ${fmt(sp.TD)} s, soil factor S = ${fmt(sp.S)}.`, 'Kinematic soil–pile interaction and liquefaction are not included.'],
      advice,
      publish: { 'dolphin.f1': modal.frequencies[0], 'dolphin.seismicShear': Vb, 'dolphin.seismicDisp': u },
      quantities: { steel_t: round(steelT) },
    };
  },
  verify: [
    // Cantilever pile L = 20 m with tip mass M = 100 t (cap 4×4×2 m at 30.6458 kN/m³ = 980.665 kN) and pile mass
    // m·L = 7850·0.061575·20 = 9667 kg. Rayleigh with the static deflected shape: M_eff = M + (33/140)·m·L = 102 279 kg,
    // k = 3EI/L³ = 582.4 kN/m → T = 2π√(M_eff/k) = 2.633 s. EN 1998-1 Type 1, ground A (S = 1, TC = 0.4, TD = 2 s), 5 %:
    // T > TD → Se = ag·S·2.5·TC·TD/T² = 2·2.5·0.4·2/2.633² = 0.5769 m/s². Tip displacement = Γ·Se/ω² with
    // Γ = (M + 0.375 m·L)/(M + 0.2357 m·L) = 1.0132: 1.0132·0.5769/(2π/2.633)² = 102.6 mm.
    { name: 'Cantilever with tip mass: T = 2π√((M + 33/140·mL)/(3EI/L³)) and EC8 spectrum ordinate', source: 'Blevins, Formulas for Natural Frequency and Mode Shape (Rayleigh, cantilever with end mass); EN 1998-1:2004 Eq. 3.5',
      inputs: { piles: [{ x: 0, y: 0, batter: 0, azimuth: 0 }], D: 1, t: 20, corr: 0, steel: 's355', LeUser: 20, shear: false, nSub: 8, capB: 4, capL: 4, capT: 2, gammaConc: 30.6458, Wequip: 0, Ca: 0, flooded: false, ag: 2, ground: 'A', specType: '1', zeta: 5, q: 1, dirRule: 'SRSS', nModes: 6 },
      expect: [{ kpi: 'T1', value: 2.633, tol: 0.003 }, { kpi: 'Sa1', value: 0.5769, tol: 0.006 }, { kpi: 'uCap', value: 102.6, tol: 0.01 }] },
  ],
};
