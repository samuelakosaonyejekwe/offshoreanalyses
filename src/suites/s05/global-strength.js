// Global and local strength of a piled wharf bent: elastic frame analysis for permanent, cargo, crane, berthing,
// mooring and environmental actions, tubular pile checks, deck-beam bending and shear, and the permissible cargo
// load, crane load and fender reaction.
import { momentCapacity, shearResistanceLinks, shearResistanceNoLinks } from '../../engine/structural/concrete.js';
import { r, f, check, ucCheck, eng, eco, maxWhere, bentInputs, bentLoadInputs, pileInputs, factorInputs, rcInputs, rcMember, bentModel, solveBent, bentCombinations, framePlot, deformScale, bentQuantities,
  steelCo2, CONFIRM, NOT_HIFI, RC_NOTE, ECON_BASIS } from './_common.js';

export default {
  id: 's05-global-strength', suite: 5, title: 'Global and local strength of the wharf bent',
  summary: 'Frame analysis of one transverse bent under cargo, crane, berthing, mooring and environmental actions: pile forces and utilisation, deck-beam moments and shears, deck sway, and the permissible cargo, crane and fender loads.',
  covers: [1], fidelity: 'numerical',
  standards: ['BS 6349-2:2019 (quay walls, jetties and dolphins)', 'EN 1990:2002+A1 Annex A1 (combinations)', 'ISO 19902:2020 §13 (tubular members)', 'EN 1992-1-1:2004 §6.1, §6.2', 'ASCE 61-14 (piers and wharves)'],
  models: ['Matrix stiffness method (3-D frame, direct stiffness)', 'Timoshenko beam elements', 'Static force and moment equilibrium', 'Equivalent-fixity pile model', 'Column buckling (ISO 19902 column curve)', 'Reinforced-concrete sectional analysis', 'Ultimate limit state combinations'],
  limits: `One bent analysed in its own plane with linear elastic members, piles built-in at an apparent fixity depth and a deck beam of reduced (cracked) stiffness; load spread between bents by the deck diaphragm is taken by the shares entered. ${RC_NOTE} ${NOT_HIFI}`,
  publishes: ['wharf.pileAxialMax', 'wharf.pileTensionMax', 'wharf.pileMomentMax', 'wharf.pileShearMax', 'wharf.deckSway', 'wharf.deckMsag', 'wharf.deckMhog', 'wharf.deckShear', 'wharf.pileD', 'wharf.pileT', 'quantities.steel_t', 'quantities.concrete_m3'],
  searchTerms: 'piled wharf structural analysis pile bent crane load berthing',
  inputs: [
    ...bentInputs(), ...pileInputs(),
    { key: 'Kpile', label: 'Pile effective-length factor (0 = automatic)', unit: '-', type: 'number', default: 0, min: 0, max: 3, group: 'Piles', role: 'criteria', help: 'Applied to the length from deck to fixity. Automatic: 1.2 for a moment-connected sway frame, 2.1 for pinned heads.' },
    ...bentLoadInputs(), ...factorInputs(), ...rcInputs(),
    { key: 'swayLim', label: 'Deck sway limit (height / value)', unit: '-', type: 'number', default: 300, min: 50, max: 1000, group: 'Serviceability', role: 'criteria', help: `Sway under the serviceability berthing combination against the height from fixity to deck divided by this value; crane rail tolerances may govern. ${CONFIRM}` },
  ],
  compute(inp, ctx) {
    const B = bentModel(inp, ctx), S = solveBent(B, inp), combos = bentCombinations(inp);
    const rc = rcMember(inp, ctx, inp.capB, inp.capH);
    const MRdSag = momentCapacity({ section: rc.sag, material: rc.material }).MRd, MRdHog = momentCapacity({ section: rc.hog, material: rc.material }).MRd;
    const results = combos.map((c) => ({ ...c, ...S.combo(c.factors) }));
    const uls = results.filter((q) => q.ls === 'ULS'), sls = results.filter((q) => q.ls === 'SLS');
    const gov = uls.reduce((a, q) => (q.ucMax > a.ucMax ? q : a), uls[0]);
    const govPile = gov.piles.reduce((a, p) => (p.uc > a.uc ? p : a), gov.piles[0]);
    const allP = uls.flatMap((q) => q.piles);
    const Ncomp = Math.max(0, ...allP.map((p) => -p.N)), Ntens = Math.max(0, ...allP.map((p) => p.N)), Mmax = Math.max(...allP.map((p) => p.M)), Vmax = Math.max(...allP.map((p) => p.V));
    const Msag = Math.max(...uls.map((q) => q.Msag)), Mhog = Math.max(...uls.map((q) => q.Mhog)), Vd = Math.max(...uls.map((q) => q.Vd));
    const shear = (VEd) => (rc.AswS > 0 ? shearResistanceLinks({ bw: rc.b, d: rc.d, fck: rc.fck, fywk: rc.fyk, AswOverS: rc.AswS, VEd, gammaC: inp.gammaC, gammaS: inp.gammaS, alphaCC: inp.alphaCC }).VRd
      : shearResistanceNoLinks({ bw: rc.b, d: rc.d, fck: rc.fck, Asl: rc.AsB, gammaC: inp.gammaC, alphaCC: inp.alphaCC }).VRdc);
    const VRd = shear(Vd);
    const deckUc = (q) => Math.max(q.Msag / MRdSag, q.Mhog / MRdHog, q.Vd / VRd);
    const swaySls = Math.max(...sls.map((q) => q.sway)), swayUls = Math.max(...uls.map((q) => q.sway));
    const hFix = Math.max(...gov.piles.map((p) => p.L)), swayAllow = hFix / inp.swayLim;
    const ucAll = (q) => Math.max(q.ucMax, deckUc(q));
    const ucMax = Math.max(...uls.map(ucAll));
    // permissible actions: scale one action in its governing ULS combination until a pile or the deck beam reaches unity
    const scaled = (c, key, k) => { const fac = { ...c.factors }; fac[key] = (fac[key] ?? 0) * k; return ucAll(S.combo(fac)); };
    const perm = (ci, key) => maxWhere((k) => scaled(combos[ci], key, k) <= 1, 0, 20, 30);
    const kQ = perm(0, 'Q'), kC = perm(0, 'C'), kB = perm(1, 'B');
    const Vcrane = (inp.crane || []).reduce((s, c) => Math.max(s, Number(c.V) || 0), 0);
    const q = bentQuantities(B, inp);
    const sc = deformScale(B.model, gov.disp);
    const deckX = [], deckM = [], deckV = [];
    const order = gov.diag.map((d) => ({ d, x0: Math.min(B.model.nodes[B.model.members[d.member].a][0], B.model.nodes[B.model.members[d.member].b][0]), flip: B.model.nodes[B.model.members[d.member].a][0] > B.model.nodes[B.model.members[d.member].b][0] })).sort((a, b) => a.x0 - b.x0);
    for (const o of order) { const L = o.d.s[o.d.s.length - 1]; for (let i = 0; i < o.d.s.length; i++) { const k = o.flip ? o.d.s.length - 1 - i : i; deckX.push(r(o.x0 + (o.flip ? L - o.d.s[k] : o.d.s[k]) - B.x0)); deckM.push(r(-o.d.My[k] / 1e3)); deckV.push(r(o.d.Vz[k] / 1e3)); } }
    const ok = ucMax <= 1, okSway = swaySls <= swayAllow;
    // wall thickness that would bring the governing pile to unity, assuming stress scales inversely with thickness
    const tNeed = inp.corr + (inp.t - inp.corr) * govPile.uc, dSteel = q.steel_t * (tNeed - inp.t) / inp.t;
    return {
      kpis: [
        { key: 'ucMax', label: 'Governing utilisation (piles and deck beam)', value: r(ucMax), unit: '-', limit: 1, sense: 'max', explain: 'Largest demand over resistance of any pile or of the deck beam in the ultimate limit state combinations. Above 1 the member does not meet the code resistance.' },
        { key: 'pileUc', label: 'Governing pile utilisation', value: r(gov.ucMax), unit: '-', limit: 1, sense: 'max', explain: `ISO 19902 interaction of axial force and bending in the most loaded pile (pile ${govPile.id}, ${gov.name}).` },
        { key: 'Ncomp', label: 'Largest pile compression', value: r(Ncomp / 1e3), unit: 'kN', explain: 'Factored axial compression at the pile head or fixity; the demand for the geotechnical pile-capacity check.' },
        { key: 'Ntens', label: 'Largest pile tension', value: r(Ntens / 1e3), unit: 'kN', explain: 'Factored uplift in any pile. Tension appears in the seaward row under berthing and in rakers; it needs shaft friction and a tension connection.' },
        { key: 'Mpile', label: 'Largest pile bending moment', value: r(Mmax / 1e3), unit: 'kN·m', explain: 'Factored moment at the pile head or at the fixity point. Short landward piles attract the most because they are stiffest.' },
        { key: 'Vpile', label: 'Largest pile shear', value: r(Vmax / 1e3), unit: 'kN', explain: 'Factored horizontal force taken by one pile; the lateral demand on the soil.' },
        { key: 'Msag', label: 'Deck beam sagging moment', value: r(Msag / 1e3), unit: 'kN·m', limit: r(MRdSag / 1e3), sense: 'max', explain: 'Largest factored span moment (tension at the bottom) against the moment resistance of the bottom reinforcement.' },
        { key: 'Mhog', label: 'Deck beam hogging moment', value: r(Mhog / 1e3), unit: 'kN·m', limit: r(MRdHog / 1e3), sense: 'max', explain: 'Largest factored support moment (tension at the top) against the moment resistance of the top reinforcement.' },
        { key: 'Vdeck', label: 'Deck beam shear', value: r(Vd / 1e3), unit: 'kN', limit: r(VRd / 1e3), sense: 'max', explain: 'Largest factored shear force in the deck beam against the resistance of the links (or of the concrete alone if there are none).' },
        { key: 'sway', label: 'Deck sway, serviceability', value: r(swaySls * 1e3), unit: 'mm', limit: r(swayAllow * 1e3), sense: 'max', explain: 'Horizontal movement of the deck under unfactored actions. It governs crane rail alignment, expansion joints and the feel of the berth.' },
        { key: 'swayUls', label: 'Deck sway, ultimate', value: r(swayUls * 1e3), unit: 'mm', explain: 'Sway under factored actions: check that P–Δ effects stay small (sway × vertical load against the lateral load × height).' },
        { key: 'Rz', label: 'Total vertical reaction (ULS operating)', value: r(Math.abs(results[0].Rz) / 1e3), unit: 'kN', explain: 'Sum of the factored vertical loads on the bent: an equilibrium check against the applied loads and the total demand on the piles.' },
        { key: 'qPerm', label: 'Maximum permissible uniform cargo load', value: r(kQ * inp.qLive), unit: 'kPa', explain: 'Characteristic cargo load at which a pile or the deck beam reaches unity in the operating combination, the crane loads unchanged. Reported at 20 times the entered load when not reached.' },
        { key: 'cranePerm', label: 'Maximum permissible crane rail load', value: r(kC * Vcrane), unit: 'kN', explain: 'Largest characteristic rail load (all crane loads scaled together) the bent can carry with the entered cargo load.' },
        { key: 'fenderPerm', label: 'Maximum permissible fender reaction on the bent', value: r(kB * inp.Ffender), unit: 'kN', explain: 'Characteristic fender reaction at which the bent reaches unity in the berthing combination. Compare with the abnormal reaction from the berthing module.' },
      ],
      checks: [
        ucCheck(`Pile ${govPile.id}: ${govPile.gov || 'axial + bending'}`, govPile.uc, govPile.ref, 'Guards against yielding or buckling of the pile between deck and fixity. Exceedance means permanent bending of the pile and loss of lateral stiffness of the wharf.'),
        check('Deck beam sagging', Msag / 1e3, MRdSag / 1e3, 'kN·m', 'EN 1992-1-1:2004 §6.1', 'Guards against flexural failure of the deck beam in the span: the bottom bars yield and the concrete crushes at the top.'),
        check('Deck beam hogging', Mhog / 1e3, MRdHog / 1e3, 'kN·m', 'EN 1992-1-1:2004 §6.1', 'Guards against flexural failure over the piles, where the pile-head moment enters the beam.'),
        check('Deck beam shear', Vd / 1e3, VRd / 1e3, 'kN', 'EN 1992-1-1:2004 §6.2.2, §6.2.3', 'Guards against diagonal shear failure next to the piles and under the crane rails, which is brittle.'),
        check('Deck sway (serviceability)', swaySls * 1e3, swayAllow * 1e3, 'mm', 'BS 6349-2:2019 §5 (serviceability); crane supplier rail tolerances', 'Guards against movements that misalign crane rails and open joints. It is a stiffness criterion: more or stiffer piles, or rakers, reduce it.'),
      ],
      plots: [
        framePlot(B.model, { disp: gov.disp, scale: sc, title: `Bent and deformed shape (${gov.name}, × ${sc})`, values: B.model.nodes.map((p, i) => Math.hypot(gov.disp[6 * i], gov.disp[6 * i + 2]) * 1e3), valueLabel: 'Displacement', valueUnit: 'mm',
          explain: 'The bent with its displacement in the governing combination, magnified. Look for the sway mode and for the short landward piles doing most of the lateral work.' }),
        { type: 'line', title: `Deck beam bending moment and shear (${gov.name})`, x: { label: 'Distance from the seaward pile', unit: 'm' }, y: { label: 'Moment (kN·m, sagging +) / shear (kN)', unit: '-' },
          series: [{ name: 'Bending moment (kN·m)', x: deckX, y: deckM }, { name: 'Shear force (kN)', x: deckX, y: deckV, dash: true }], explain: 'Moment (sagging positive) and shear along the deck beam. Peaks under the crane rails and jumps at the pile heads show where reinforcement is needed.' },
        { type: 'bar', title: 'Pile utilisation by combination', categories: gov.piles.map((p) => `Pile ${p.id}${p.group === 'raker' ? ' (raker)' : ''}`), series: uls.map((qq) => ({ name: qq.name, y: qq.piles.map((p) => r(p.uc, 3)) })), x: { label: 'Pile (seaward first)', unit: '' }, y: { label: 'Utilisation', unit: '-' },
          explain: 'ISO 19902 utilisation of each pile in each ultimate combination. Bars near 1 are fully used; a large spread suggests the pile sizes or lengths could be rebalanced.' },
      ],
      tables: [
        { title: `Pile forces (${gov.name})`, columns: [{ key: 'id', label: 'Pile' }, { key: 'x', label: 'x', unit: 'm' }, { key: 'N', label: 'Axial (tension +)', unit: 'kN' }, { key: 'M', label: 'Moment', unit: 'kN·m' }, { key: 'V', label: 'Shear', unit: 'kN' }, { key: 'L', label: 'Length to fixity', unit: 'm' }, { key: 'uc', label: 'Utilisation' }, { key: 'gov', label: 'Governing check' }],
          rows: gov.piles.map((p) => ({ id: p.group === 'raker' ? `${p.id} raker` : p.id, x: r(p.x - B.x0, 3), N: r(p.N / 1e3), M: r(p.M / 1e3), V: r(p.V / 1e3), L: r(p.L, 3), uc: r(p.uc, 3), gov: p.gov })), explain: 'Factored member forces of each pile in the governing combination; use the axial forces for the geotechnical check and the moments for the pile-head connection.' },
        { title: 'Combination summary', columns: [{ key: 'name', label: 'Combination' }, { key: 'uc', label: 'Pile utilisation' }, { key: 'deck', label: 'Deck beam utilisation' }, { key: 'Msag', label: 'Sagging', unit: 'kN·m' }, { key: 'Mhog', label: 'Hogging', unit: 'kN·m' }, { key: 'V', label: 'Shear', unit: 'kN' }, { key: 'sway', label: 'Sway', unit: 'mm' }],
          rows: results.map((qq) => ({ name: qq.name, uc: r(qq.ucMax, 3), deck: r(deckUc(qq), 3), Msag: r(qq.Msag / 1e3), Mhog: r(qq.Mhog / 1e3), V: r(qq.Vd / 1e3), sway: r(qq.sway * 1e3, 3) })), explain: 'One row per load combination with the factors entered. The row with the largest utilisation governs the member sizes; the serviceability rows govern stiffness.' },
      ],
      notes: [`${gov.name} governs the piles (pile ${govPile.id}, ${govPile.gov || 'combined check'}); effective-length factor ${f(S.K)} on the length to fixity.`,
        `Deck beam resistances: sagging ${f(MRdSag / 1e3)} kN·m, hogging ${f(MRdHog / 1e3)} kN·m, shear ${f(VRd / 1e3)} kN with γc = ${f(inp.gammaC)}, γs = ${f(inp.gammaS)}, αcc = ${f(inp.alphaCC)}.`],
      advice: [
        ok ? eng(okSway ? 'ok' : 'warn', `Governing utilisation ${f(ucMax)} (${gov.name}); deck sway ${f(swaySls * 1e3)} mm against ${f(swayAllow * 1e3)} mm. The bent can carry up to ${f(kQ * inp.qLive)} kPa of cargo, ${f(kC * Vcrane)} kN on a crane rail or ${f(kB * inp.Ffender)} kN of fender reaction.`,
          okSway ? `Strength reserve of ${f((1 - ucMax) * 100)} %: a wall of about ${f(Math.max(tNeed, 10))} mm would bring the governing pile to unity. Confirm the fixity depth with the soil–structure module before trimming.` : 'Strength is adequate but the deck is too flexible: add a pair of raker piles, increase the pile diameter or reduce the free length.', 'BS 6349-2:2019; ISO 19902:2020 §13; EN 1992-1-1:2004',
          `Reducing the wall to ${f(Math.max(tNeed, 10))} mm would change the bent steel by ${f(q.steel_t * (Math.max(tNeed, 10) - inp.t) / inp.t)} t: ${steelCo2(q.steel_t * (Math.max(tNeed, 10) - inp.t) / inp.t)}.`)
          : eng('fail', `Utilisation ${f(ucMax)} in ${gov.name}: ${gov.ucMax >= deckUc(gov) ? `pile ${govPile.id} at ${f(govPile.uc)} (${govPile.gov})` : `deck beam at ${f(deckUc(gov))}`}.`,
            gov.ucMax >= deckUc(gov) ? `Increase the pile wall to about ${f(Math.ceil(tNeed))} mm (or a higher steel grade), add rakers to take the horizontal load in axial force, or limit the cargo load to ${f(kQ * inp.qLive)} kPa.` : `Add reinforcement or depth to the deck beam: sagging ${f(Msag / MRdSag)}, hogging ${f(Mhog / MRdHog)}, shear ${f(Vd / VRd)} of resistance.`, 'ISO 19902:2020 §13; EN 1992-1-1:2004 §6',
            `Thickening the piles to ${f(Math.ceil(tNeed))} mm adds about ${f(dSteel)} t per bent: ${steelCo2(dSteel)}. Rakers carry the same horizontal load with far less steel than bending of vertical piles.`),
        eco('info', `One bent with its deck uses about ${f(q.steel_t)} t of pile steel and ${f(q.concrete_m3)} m³ of concrete; the permissible loads are ${f(kQ * inp.qLive)} kPa cargo and ${f(kC * Vcrane)} kN per crane rail.`,
          'Fix the cargo and crane ratings from the permissible values rather than the design values: unused capacity is revenue (heavier cranes, higher stacks) at no construction cost, while a shortfall is far cheaper to resolve by an operating limit than by strengthening.', ECON_BASIS),
      ],
      publish: { 'wharf.pileAxialMax': Ncomp, 'wharf.pileTensionMax': Ntens, 'wharf.pileMomentMax': Mmax, 'wharf.pileShearMax': Vmax, 'wharf.deckSway': swaySls, 'wharf.deckMsag': Msag, 'wharf.deckMhog': Mhog, 'wharf.deckShear': Vd, 'wharf.pileD': inp.D, 'wharf.pileT': inp.t / 1e3, 'quantities.steel_t': q.steel_t, 'quantities.concrete_m3': q.concrete_m3 },
      quantities: q,
    };
  },
  verify: [
    // Single moment-connected pile, free to sway: cantilever from fixity. L = 4 + 12 + 5 = 21 m; tube 914 × 22: I = π/64·(0.914⁴ − 0.870⁴) = 6.1357e-3 m⁴,
    // E = 210 GPa → EI = 1.2885e9 N·m²; δ = H L³/(3EI) = 5e5·9261/(3·1.2885e9) = 1.1979 m bending, plus shear H L/(G As) ≈ 0.4 % → about 1.203 m.
    { name: 'Single-pile bent as a cantilever: sway = H L³/(3EI)', source: 'Euler–Bernoulli cantilever (Roark Table 8.1 case 1a) with the Timoshenko shear term',
      inputs: { D: 0.914, t: 22, nPiles: 1, corr: 0, fixDepth: 5, gSdl: 0, qLive: 0, crane: [], Ffender: 500, Fmoor: 0, Fenv: 0, psi: 0, slabT: 0.2 }, expect: [{ kpi: 'sway', value: 1203, tol: 0.006 }] },
    // Two equal piles, level bed, rigid-enough beam, pinned heads: each pile is a cantilever carrying H/2: M = (H/2)·L = 250e3·21 = 5250 kN·m (γ = 1.5 → 7875).
    { name: 'Two pinned-head piles share the fender load: M = γ (H/2) L', source: 'Statics of a pinned portal with equal columns',
      inputs: { D: 0.914, t: 22, nPiles: 2, bedRise: 0, headFix: 'pinned', corr: 0, fixDepth: 5, gSdl: 0, qLive: 0, crane: [], Ffender: 500, Fmoor: 0, Fenv: 0, psi: 0, gG: 1, gB: 1.5 }, expect: [{ kpi: 'Mpile', value: 7875, tol: 0.005 }, { kpi: 'Vpile', value: 375, tol: 0.005 }] },
  ],
};
