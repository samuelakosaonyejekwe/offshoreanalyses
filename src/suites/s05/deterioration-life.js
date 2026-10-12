// Structural deterioration and remaining life of a wharf: chloride ingress into the concrete deck and the time to
// corrosion initiation (deterministic and probabilistic), loss of bar section and of moment resistance with time,
// wall loss of the steel piles by exposure zone, growth of a recorded weld crack, and the remaining service life
// with the repair or protection that extends it.
import { rcSection, chlorideProfile, chlorideInitiationTime, corrosionRate, barSectionLoss, residualCapacityHistory, pileCorrosionLoss } from '../../engine/structural/concrete.js';
import { parisLawBs7910, parisLifeClosedForm } from '../../engine/fatigue/fracture.js';
import { form } from '../../engine/reliability/form.js';
import { erfc } from '../../engine/math/special.js';
import { YEAR_S, need, run, r, f, cap, check, sfCheck, eng, eco, linspace, CONFIRM, ECON_BASIS, CO2_CONC, CO2_STEEL } from './_common.js';

const ZONES_DEFAULT = [
  { name: 'Atmospheric (under deck)', rate: 0.05, protect: 15, faces: 1 }, { name: 'Splash zone', rate: 0.3, protect: 20, faces: 1 }, { name: 'Tidal zone', rate: 0.1, protect: 20, faces: 1 },
  { name: 'Low-water zone', rate: 0.2, protect: 25, faces: 1 }, { name: 'Permanent immersion', rate: 0.1, protect: 25, faces: 1 }, { name: 'Below the seabed', rate: 0.02, protect: 0, faces: 1 },
];
const normCdf = (x) => 0.5 * erfc(-x / Math.SQRT2);

export default {
  id: 's05-deterioration-life', suite: 5, title: 'Deterioration and remaining service life',
  summary: 'When chlorides reach the deck reinforcement, how fast bars and pile walls lose section, how the moment resistance and reliability fall with age, how long a recorded crack takes to grow through, and how many years of service remain before repair or strengthening.',
  covers: [15], fidelity: 'analytical',
  standards: ['fib Bulletin 34 Model Code for Service Life Design (2006) §B2 (chloride ingress)', 'fib Model Code 2010 §7.8', 'EN 1992-1-1:2004 §4 (durability and cover)', 'EN 1993-5:2007 Tables 4-1, 4-2 (corrosion of steel piling)', 'BS 6349-1-4:2021 §4 (durability of maritime structures)', 'BS 7910:2019 clause 8 (fatigue crack growth)', 'ISO 13822:2010; ISO 2394:2015 (assessment of existing structures, reliability)'],
  models: ["Fick's second law of diffusion with an ageing diffusion coefficient", 'Chloride-induced corrosion initiation (critical chloride content)', "Corrosion propagation by Faraday's law", 'Bar section loss with pitting factor', 'Residual moment resistance of the corroding section', 'Zone-wise wall loss of steel piles after coating or cathodic protection life', 'Paris crack-growth law (closed form)', 'Time-dependent structural reliability (first-order reliability method)'],
  limits: 'Chloride ingress is one-dimensional Fickian diffusion with an ageing factor and user-supplied surface content, diffusion coefficient and critical content: convection in the splash zone, cracks, carbonation and temperature cycles enter only through those parameters. Bond loss, spalling and loss of bar ductility are not modelled, so the strength curve is optimistic once cracking starts. Pile corrosion rates and the crack-growth law are editable defaults. Calibrate against chloride profiles, half-cell and wall-thickness surveys.',
  publishes: ['life.remainingYears', 'wharf.initiationYears', 'wharf.pileWallRemaining'],
  searchTerms: 'chloride ingress service life marine concrete corrosion steel pile wall loss remaining life',
  inputs: [
    { key: 'age', label: 'Present age of the structure', unit: 'years', type: 'number', default: 20, min: 0, max: 150, group: 'Service history', role: 'ic', help: 'Zero for a new design. For reassessment, the time since construction or since the last full repair of the element considered.' },
    { key: 'life', label: 'Required service life from construction', unit: 'years', type: 'number', default: 50, min: 5, max: 200, group: 'Service history', role: 'criteria', bind: 'life.designYears' },
    { key: 'b', label: 'Deck beam width', unit: 'm', type: 'number', default: 1.2, min: 0.2, max: 5, group: 'Concrete member', role: 'geometry' },
    { key: 'h', label: 'Deck beam depth', unit: 'm', type: 'number', default: 1.4, min: 0.2, max: 5, group: 'Concrete member', role: 'geometry' },
    { key: 'nBot', label: 'Bottom bars', unit: '-', type: 'int', default: 12, min: 2, max: 80, group: 'Concrete member', role: 'geometry' },
    { key: 'dBot', label: 'Bottom bar diameter', unit: 'mm', type: 'number', default: 32, min: 8, max: 50, group: 'Concrete member', role: 'geometry' },
    { key: 'cover', label: 'Cover to the main bars (mean, as built)', unit: 'mm', type: 'number', default: 70, min: 15, max: 150, group: 'Concrete member', role: 'geometry', design: true, uncertain: { dist: 'normal', cov: 0.12 }, help: 'From the cover survey for an existing structure; nominal cover for a new one.' },
    { key: 'coverStd', label: 'Standard deviation of the cover', unit: 'mm', type: 'number', default: 8, min: 1, max: 30, group: 'Concrete member', role: 'ic', help: '6–10 mm for ordinary site control, less with a cover survey of the actual element.' },
    { key: 'conc', label: 'Concrete class', type: 'select', source: 'materials:concrete', default: 'c40', group: 'Concrete member', role: 'material' },
    { key: 'rebar', label: 'Reinforcing steel', type: 'select', source: 'materials:steel', default: 'reinf-b500', group: 'Concrete member', role: 'material' },
    { key: 'MEd', label: 'Design moment the member must keep resisting', unit: 'kN·m', type: 'number', default: 3300, min: 0, max: 1e5, group: 'Concrete member', role: 'bc', help: 'Governing ULS moment from the strength module.' },
    { key: 'gammaC', label: 'γc concrete', unit: '-', type: 'number', default: 1.5, min: 1, max: 2, group: 'Concrete member', role: 'criteria', help: CONFIRM },
    { key: 'gammaS', label: 'γs reinforcement', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.5, group: 'Concrete member', role: 'criteria', help: CONFIRM },
    { key: 'alphaCC', label: 'αcc', unit: '-', type: 'number', default: 0.85, min: 0.7, max: 1, group: 'Concrete member', role: 'criteria', help: CONFIRM },
    { key: 'Cs', label: 'Surface chloride content', unit: '% binder', type: 'number', default: 3.0, min: 0.2, max: 8, group: 'Chloride exposure', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Chloride content at the surface (or at the end of the convection zone) by mass of binder: about 2–4 % in the splash and tidal zones, 1–2 % in the atmospheric zone under a deck. Use the measured profile where there is one.' },
    { key: 'Ccrit', label: 'Critical chloride content', unit: '% binder', type: 'number', default: 0.6, min: 0.1, max: 2.5, group: 'Chloride exposure', role: 'criteria', uncertain: { dist: 'lognormal', cov: 0.25 }, help: `fib Bulletin 34 uses a mean of 0.6 % by mass of binder; 0.4 % is a common deterministic threshold. ${CONFIRM}` },
    { key: 'Ci', label: 'Initial chloride content of the mix', unit: '% binder', type: 'number', default: 0.05, min: 0, max: 0.4, group: 'Chloride exposure', role: 'ic' },
    { key: 'D0', label: 'Apparent diffusion coefficient at the reference age', unit: '×10⁻¹² m²/s', type: 'number', default: 5, min: 0.1, max: 60, group: 'Chloride exposure', role: 'material', design: true, uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'From rapid chloride migration tests at the reference age (times the environmental and test factors): about 2–5 for slag or fly-ash blends, 8–15 for plain Portland cement at w/c 0.45.' },
    { key: 't0', label: 'Reference age of that coefficient', unit: 'days', type: 'number', default: 28, min: 1, max: 3650, group: 'Chloride exposure', role: 'material' },
    { key: 'mAge', label: 'Ageing exponent', unit: '-', type: 'number', default: 0.4, min: 0, max: 0.8, group: 'Chloride exposure', role: 'material', help: 'Da(t) = D0·(t0/t)^m: about 0.3 for Portland cement, 0.45 slag, 0.6 fly ash (fib Bulletin 34). Zero keeps the coefficient constant.' },
    { key: 'dxConv', label: 'Depth of the convection zone', unit: 'mm', type: 'number', default: 0, min: 0, max: 20, group: 'Chloride exposure', role: 'bc', help: 'About 8–10 mm in the splash zone, 0 where the surface stays wet.' },
    { key: 'profile', label: 'Measured chloride profile (optional)', type: 'table', group: 'Chloride exposure', role: 'ic', columns: [{ key: 'x', label: 'Depth', unit: 'mm', type: 'number' }, { key: 'C', label: 'Chloride content', unit: '% binder', type: 'number' }], default: [],
      help: 'Dust samples from the present structure. When given, the apparent diffusion coefficient that fits them at the present age is reported next to the one entered.' },
    { key: 'icorr', label: 'Corrosion current density after initiation', unit: 'µA/cm²', type: 'number', default: 1.0, min: 0.05, max: 10, group: 'Corrosion propagation', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.5 }, help: '0.1–0.5 low, 0.5–1 moderate, above 1 high (tidal and splash zones); 1 µA/cm² removes 11.6 µm of steel a year.' },
    { key: 'pit', label: 'Pitting factor', unit: '-', type: 'number', default: 5, min: 1, max: 12, group: 'Corrosion propagation', role: 'criteria', help: 'Deepest pit over mean penetration: 4–8 for chloride attack. Reported as pit depth; the strength curve uses the mean loss.' },
    { key: 'pCrack', label: 'Mean penetration at cover cracking', unit: 'µm', type: 'number', default: 100, min: 10, max: 500, group: 'Corrosion propagation', role: 'criteria', help: 'Rust expansion cracks the cover at about 50–100 µm of bar loss for ordinary cover-to-diameter ratios; the visible end of the maintenance-free life.' },
    { key: 'betaT', label: 'Target reliability index against depassivation', unit: '-', type: 'number', default: 1.3, min: 0.5, max: 3.8, group: 'Corrosion propagation', role: 'criteria', help: `fib Bulletin 34 uses about 1.3 (a 10 % chance of initiation) as a serviceability target. ${CONFIRM}` },
    { key: 'D', label: 'Pile outer diameter', unit: 'm', type: 'number', default: 1.067, min: 0.2, max: 3, group: 'Steel piles', role: 'geometry', bind: 'wharf.pileD' },
    { key: 't', label: 'Pile wall thickness as built', unit: 'mm', type: 'number', default: 25, min: 6, max: 80, group: 'Steel piles', role: 'geometry', bind: 'wharf.pileT' },
    { key: 'tMin', label: 'Minimum wall needed for strength', unit: 'mm', type: 'number', default: 15, min: 3, max: 80, group: 'Steel piles', role: 'criteria', design: true, help: 'Wall at which the governing pile reaches unity: as-built wall times its utilisation from the strength module (corrosion allowance removed).' },
    { key: 'zones', label: 'Corrosion zones of the piles', type: 'table', group: 'Steel piles', role: 'bc', columns: [{ key: 'name', label: 'Zone', type: 'text' }, { key: 'rate', label: 'Rate per face', unit: 'mm/yr', type: 'number' }, { key: 'protect', label: 'Life of coating or cathodic protection', unit: 'years', type: 'number' }, { key: 'faces', label: 'Corroding faces', type: 'number' }],
      default: ZONES_DEFAULT, help: 'Indicative rates for unprotected carbon steel in temperate sea water; EN 1993-5 Tables 4-1 and 4-2 and BS 6349-1-4 give design values, and accelerated low-water corrosion can exceed 0.5 mm/yr. Replace with site measurements.' },
    { key: 'a0', label: 'Depth of a recorded weld crack (0 = none)', unit: 'mm', type: 'number', default: 0, min: 0, max: 40, group: 'Crack record', role: 'ic', help: 'Surface-breaking crack found at a pile-head or splice weld by inspection.' },
    { key: 'dSig', label: 'Equivalent stress range at the crack', unit: 'MPa', type: 'number', default: 40, min: 1, max: 400, group: 'Crack record', role: 'bc', help: 'Constant-amplitude equivalent of the service spectrum at the detail (from the fatigue module).' },
    { key: 'nCyc', label: 'Stress cycles per year', unit: '1/yr', type: 'number', default: 500000, min: 1, max: 1e8, group: 'Crack record', role: 'bc' },
    { key: 'Ygeo', label: 'Crack geometry factor Y', unit: '-', type: 'number', default: 1.12, min: 0.5, max: 3, group: 'Crack record', role: 'criteria', help: '1.12 for a shallow edge crack; weld-toe magnification raises it for small cracks.' },
  ],
  compute(inp, ctx) {
    const cm = ctx?.material ? ctx.material(inp.conc) : null, sm = ctx?.material ? ctx.material(inp.rebar) : null;
    need(cm && cm.fck > 0, `Concrete '${inp.conc}' is not in the materials library.`); need(sm && sm.fy > 0, `Reinforcing steel '${inp.rebar}' is not in the materials library.`);
    need(inp.Cs > inp.Ci, 'The surface chloride content must exceed the initial content of the mix.');
    need(inp.tMin < inp.t, 'The minimum wall for strength must be less than the as-built wall.');
    // SI: mm → m, 1e-12 m²/s, days and years → s
    const c = inp.cover / 1e3, D0 = inp.D0 * 1e-12, t0 = inp.t0 * 86400, dx = inp.dxConv / 1e3, life = inp.life, age = inp.age;
    need(c > dx, 'The cover must exceed the depth of the convection zone.');
    const chl = { Cs: inp.Cs, Ci: inp.Ci, D0, t0, m: inp.mAge, dx };
    const tiS = chlorideInitiationTime({ cover: c, Ccrit: inp.Ccrit, ...chl }), ti = Number.isFinite(tiS) ? tiS / YEAR_S : 1e4;
    const Cat = (x, yr) => (yr > 0 ? chlorideProfile({ x: Math.max(x, dx), t: yr * YEAR_S, ...chl }).C : inp.Ci);
    const Cnow = Cat(c, Math.max(age, 1e-6)), Clife = Cat(c, life);
    // propagation
    const rate = corrosionRate({ icorr: inp.icorr * 0.01 }), rateYr = rate * YEAR_S, db = inp.dBot / 1e3;
    const tCrack = ti + inp.pCrack * 1e-6 / rateYr;
    const lossAt = (yr) => barSectionLoss({ diameter: db, rate, t: yr * YEAR_S, initiationTime: Math.min(ti, 1e4) * YEAR_S, pittingFactor: inp.pit });
    const bl = lossAt(life);
    // residual moment resistance of the beam (bottom bars corroding)
    const AsB = inp.nBot * Math.PI * db * db / 4, dEff = inp.h - c - db / 2;
    need(dEff > 0.1, 'The section is too shallow for the cover and bar diameter given.');
    const section = rcSection({ shape: 'rect', b: inp.b, h: inp.h, bars: [{ depth: dEff, area: AsB, diameter: db }] });
    const material = { fck: cm.fck, fyk: sm.fy, gammaC: inp.gammaC, gammaS: inp.gammaS, alphaCC: inp.alphaCC };
    const horizon = Math.max(1.6 * life, age + 30, 60), yrs = linspace(0, horizon, 49);
    const hist = run(() => residualCapacityHistory({ section, material, times: yrs.map((y) => y * YEAR_S), initiationTime: Math.min(ti, 1e4) * YEAR_S, rate }), 'Residual capacity');
    const MRd0 = hist[0].MRd, MEd = inp.MEd * 1e3;
    const MRdAt = (yr) => { const k = Math.min(yrs.length - 2, Math.max(0, Math.floor(yr / horizon * (yrs.length - 1)))); const q = (yr - yrs[k]) / (yrs[k + 1] - yrs[k]); return hist[k].MRd + q * (hist[k + 1].MRd - hist[k].MRd); };
    let tStr = 1e4; if (MRd0 < MEd) tStr = 0; else for (let i = 1; i < hist.length; i++) if (hist[i].MRd < MEd) { tStr = yrs[i - 1] + (yrs[i] - yrs[i - 1]) * (hist[i - 1].MRd - MEd) / (hist[i - 1].MRd - hist[i].MRd); break; }
    if (tStr === 1e4 && MEd > 0 && Number.isFinite(tiS)) { // beyond the plotted horizon: bar area needed scales with the moment (lever arm nearly constant)
      const need2 = Math.sqrt(Math.min(MEd / MRd0, 1)); tStr = ti + (db * (1 - need2) / 2) / rateYr;
    }
    // probability of depassivation with time (FORM on cover, diffusion coefficient, surface and critical content)
    const vars = { cover: { dist: 'normal', mean: c, std: inp.coverStd / 1e3 }, D0: { dist: 'lognormal', mean: D0 * 1e12, cov: 0.3 }, Cs: { dist: 'lognormal', mean: inp.Cs, cov: 0.3 }, Cc: { dist: 'lognormal', mean: inp.Ccrit, cov: 0.25 } };
    const betaAt = (yr) => { if (yr <= 0) return 8; const g = (x) => x.Cc - chlorideProfile({ x: Math.max(x.cover, dx + 1e-4), t: yr * YEAR_S, Cs: Math.max(x.Cs, inp.Ci + 1e-6), Ci: inp.Ci, D0: Math.max(x.D0, 1e-4) * 1e-12, t0, m: inp.mAge, dx }).C; try { return Math.max(-8, Math.min(8, form(g, vars, { maxIter: 100 }).beta)); } catch { return NaN; } };
    const ty = linspace(horizon / 12, horizon, 12), betas = ty.map(betaAt);
    need(betas.every(Number.isFinite), 'The reliability analysis of chloride ingress did not converge: check that the chloride parameters are physically consistent.');
    const bNow = age > 0 ? betaAt(age) : 8, bLife = betaAt(life);
    let tRel = horizon; for (let i = 0; i < ty.length; i++) if (betas[i] < inp.betaT) { tRel = i === 0 ? ty[0] * Math.max(0, (8 - inp.betaT) / (8 - betas[0])) : ty[i - 1] + (ty[i] - ty[i - 1]) * (betas[i - 1] - inp.betaT) / (betas[i - 1] - betas[i]); break; }
    // steel piles by zone
    const zrows = (inp.zones || []).map((z, i) => { need(z.rate >= 0 && z.protect >= 0, `Corrosion zone ${i + 1}: rate and protection life must not be negative.`); return { name: String(z.name || `Zone ${i + 1}`), rate: Number(z.rate) / 1e3 / YEAR_S, protectionLife: Number(z.protect) * YEAR_S, faces: Math.max(1, Math.round(Number(z.faces) || 1)), rateMm: Number(z.rate), protect: Number(z.protect) }; });
    need(zrows.length >= 1, 'Give at least one corrosion zone for the piles.');
    const zl = pileCorrosionLoss({ zones: zrows, t: life * YEAR_S, thickness: inp.t / 1e3, D: inp.D }), zn = pileCorrosionLoss({ zones: zrows, t: Math.max(age, 1e-6) * YEAR_S, thickness: inp.t / 1e3, D: inp.D });
    const zLife = zrows.map((z, i) => ({ ...z, lossLife: zl[i].loss, remLife: zl[i].remaining, remNow: zn[i].remaining, modulus: zl[i].modulusRatio ?? 1, end: z.rateMm > 0 ? z.protect + (inp.t - inp.tMin) / (z.rateMm * z.faces) : 1e4 }));
    const govZ = zLife.reduce((m, z) => (z.end < m.end ? z : m), zLife[0]), tPile = govZ.end;
    // recorded crack: Paris growth from a0 to the wall thickness
    let tCrk = 1e4;
    if (inp.a0 > 0) { need(inp.a0 < inp.t, 'The recorded crack depth must be less than the wall thickness.'); const law = parisLawBs7910({ environment: 'marine', curve: 'simplified' }).stages[0]; const N = parisLifeClosedForm({ a0: inp.a0 / 1e3, af: inp.t / 1e3, stressRange: inp.dSig * 1e6, Y: inp.Ygeo, A: law.A, m: law.m }); tCrk = age + N / inp.nCyc; }
    // service life: first of the limit states, measured from construction
    const lim = [{ name: 'Cover cracking of the deck (corrosion visible)', t: tCrack }, { name: 'Deck beam moment resistance below the design moment', t: tStr }, { name: `Pile wall at its minimum in the ${govZ.name.toLowerCase()}`, t: tPile }, ...(inp.a0 > 0 ? [{ name: 'Recorded crack through the pile wall', t: tCrk }] : [])];
    const gov = lim.reduce((m, q) => (q.t < m.t ? q : m), lim[0]), endOfLife = gov.t, remaining = Math.max(endOfLife - age, 0);
    // measured profile: apparent diffusion coefficient that fits it at the present age (least squares on the erfc argument)
    let Dfit = 0; const prof = (inp.profile || []).filter((q) => Number(q.x) > inp.dxConv && Number(q.C) > inp.Ci && Number(q.C) < inp.Cs);
    if (prof.length >= 2 && age > 0) { const inv = (y) => { let a = 0, b = 6; for (let i = 0; i < 60; i++) { const m = 0.5 * (a + b); if (erfc(m) > y) a = m; else b = m; } return 0.5 * (a + b); }; let sxy = 0, sxx = 0; for (const q of prof) { const z = inv((q.C - inp.Ci) / (inp.Cs - inp.Ci)), x = (q.x / 1e3 - dx); sxy += x * z; sxx += z * z; } const s = sxy / sxx; /* x = 2√(Da t)·z */ Dfit = (s / 2) ** 2 / (age * YEAR_S); }
    const DaNow = age > 0 ? D0 * (t0 / (age * YEAR_S)) ** inp.mAge : D0;
    // cover that would delay initiation to the required life (bisection)
    let cNeed = c; if (ti < life) { let a = c, b = 0.2; for (let i = 0; i < 40; i++) { const m = 0.5 * (a + b); if (chlorideInitiationTime({ cover: m, Ccrit: inp.Ccrit, ...chl }) / YEAR_S < life) a = m; else b = m; } cNeed = b; }
    const checks = [
      sfCheck('Service life against the required life', endOfLife, life, 'ISO 13822:2010 §6; fib Model Code 2010 §7.8', 'Age at which the first deterioration limit state is reached over the required service life. Below 1 the structure needs repair, protection or strengthening before the end of its intended life.'),
      check('Chloride content at the reinforcement at the end of life', Clife, inp.Ccrit, '% binder', 'fib Bulletin 34 (2006) §B2; EN 206 / BS 8500 durability provisions', 'Chloride content at the bar depth against the content that breaks down the passive film. Above it corrosion starts and the maintenance-free period is over.'),
      sfCheck('Reliability against depassivation at the end of life', Math.max(bLife, 0.01), inp.betaT, 'fib Bulletin 34 (2006) §3.2; ISO 2394:2015', 'Reliability index of the event that chlorides reach the critical content at the bars, with scatter in cover, diffusion and chloride contents. Below the target more than the accepted share of the surface is expected to be corroding.'),
      check('Deck beam moment resistance at the end of life', MEd / 1e3, MRdAt(Math.min(life, horizon)) / 1e3, 'kN·m', 'EN 1992-1-1:2004 §6.1 with reduced bar area; ISO 13822:2010', 'Design moment against the resistance of the section with corroded bars. Once it is exceeded the deck needs a load restriction or strengthening.'),
      check(`Pile wall loss at the end of life (${govZ.name})`, inp.tMin, Math.max(...zLife.map(() => 0), Math.min(...zLife.map((z) => z.remLife)) * 1e3), 'mm', 'EN 1993-5:2007 §4.4; BS 6349-1-4:2021', 'Wall needed for strength against the wall left in the worst zone. Below it the pile is overstressed at that level: local buckling in the splash or low-water zone is the usual failure.'),
    ];
    const w = checks.reduce((m, q) => (q.uc > m.uc ? q : m)), ok = w.uc <= 1;
    const xs = linspace(0, Math.max(1.6 * c, 0.1), 33);
    const repairArea = 2 * (inp.b + inp.h), concRepair = repairArea * 0.1; // m² and m³ of soffit and side repair per metre of beam, 100 mm deep
    const tPerM = Math.PI * inp.D * (inp.t / 1e3) * 7.85;
    return {
      kpis: [
        { key: 'ti', label: 'Time to corrosion initiation', value: r(cap(ti, 9999)), unit: 'years', explain: 'Age at which the chloride content at the bars reaches the critical value (mean parameters). Until then the reinforcement is protected by the alkaline concrete.' },
        { key: 'Cnow', label: 'Chloride content at the bars today', value: r(Cnow), unit: '% binder', limit: inp.Ccrit, sense: 'max', explain: `Predicted content at ${f(inp.cover)} mm after ${f(age)} years, against the critical content. Compare with dust samples from the structure.` },
        { key: 'Clife', label: 'Chloride content at the bars at the end of life', value: r(Clife), unit: '% binder', explain: 'Predicted content at the bar depth at the required service life.' },
        { key: 'corrRate', label: 'Bar corrosion rate after initiation', value: r(rateYr * 1e6), unit: 'µm/yr', explain: "Mean loss of radius from the corrosion current by Faraday's law. The deepest pits advance several times faster." },
        { key: 'tCrack', label: 'Time to cracking of the cover', value: r(cap(tCrack, 9999)), unit: 'years', explain: 'Initiation plus the time for rust to crack the cover: when rust stains and delamination appear and repair is needed to stop the decay accelerating.' },
        { key: 'lossAtLife', label: 'Bar section lost at the end of life', value: r(bl.lossFraction * 100), unit: '%', explain: `Mean loss of bar area at the required life; the deepest pit would be about ${f(bl.pitDepth * 1e3)} mm on a ${f(inp.dBot)} mm bar.` },
        { key: 'MRdLife', label: 'Moment resistance at the end of life', value: r(MRdAt(Math.min(life, horizon)) / 1e3), unit: 'kN·m', limit: r(MEd / 1e3), sense: 'min', explain: `Resistance with the corroded bar area (${f(MRd0 / 1e3)} kN·m as built). It must stay above the design moment.` },
        { key: 'betaNow', label: 'Reliability index against depassivation today', value: r(bNow), unit: '-', explain: `Probability that corrosion has already started somewhere on the surface: ${f(normCdf(-bNow) * 100)} %. It indicates how much of the soffit a condition survey should expect to find active.` },
        { key: 'betaLife', label: 'Reliability index against depassivation at the end of life', value: r(bLife), unit: '-', limit: inp.betaT, sense: 'min', explain: `Reliability index at the required life (probability of initiation ${f(normCdf(-bLife) * 100)} %), against the serviceability target.` },
        { key: 'tRel', label: 'Age at which the reliability target is reached', value: r(tRel), unit: 'years', explain: 'Probabilistic end of the initiation period: the age at which the chance of depassivation exceeds the accepted level.' },
        { key: 'tRemain', label: 'Least pile wall at the end of life', value: r(Math.min(...zLife.map((z) => z.remLife)) * 1e3), unit: 'mm', limit: inp.tMin, sense: 'min', explain: `Remaining wall in the worst zone (${zLife.reduce((m, z) => (z.remLife < m.remLife ? z : m), zLife[0]).name}) after coating or cathodic protection has expired and the steel has corroded freely.` },
        { key: 'tPile', label: 'Age at which the pile wall reaches its minimum', value: r(cap(tPile, 9999)), unit: 'years', explain: `End of life of the piles by wall loss, governed by the ${govZ.name.toLowerCase()}.` },
        { key: 'tCrk', label: 'Age at which the recorded crack is through the wall', value: r(cap(tCrk, 9999)), unit: 'years', explain: inp.a0 > 0 ? 'Present age plus the Paris-law growth time from the recorded depth to the wall thickness. Repair before then.' : 'No crack recorded: reported at the cap value.' },
        { key: 'life', label: 'Estimated service life', value: r(cap(endOfLife, 9999)), unit: 'years', limit: life, sense: 'min', explain: `Age at which the first limit state is reached: ${gov.name}.` },
        { key: 'remaining', label: 'Estimated remaining service life', value: r(cap(remaining, 9999)), unit: 'years', explain: 'Years from today until that limit state, without intervention. It is the window for planning repair, protection or strengthening.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Chloride profile in the concrete', x: { label: 'Depth from the surface', unit: 'mm' }, y: { label: 'Chloride content', unit: '% binder' },
          series: [...(age > 0 ? [{ name: `Today (${f(age)} years)`, x: xs.map((x) => r(x * 1e3)), y: xs.map((x) => r(Cat(x, age))) }] : []), { name: `End of life (${f(life)} years)`, x: xs.map((x) => r(x * 1e3)), y: xs.map((x) => r(Cat(x, life))) }, { name: 'Critical content', x: [0, r(xs[xs.length - 1] * 1e3)], y: [inp.Ccrit, inp.Ccrit], dash: true }, { name: 'Bar depth', x: [inp.cover, inp.cover], y: [0, inp.Cs], dash: true }, ...(prof.length ? [{ name: 'Measured', x: prof.map((q) => Number(q.x)), y: prof.map((q) => Number(q.C)) }] : [])],
          explain: 'Chloride content against depth now and at the end of life. Corrosion starts when a profile crosses the critical content at the bar depth (the intersection of the two dashed lines).' },
        { type: 'line', title: 'Moment resistance and reliability with age', x: { label: 'Age', unit: 'years' }, y: { label: 'Resistance ratio / reliability index', unit: '-' },
          series: [{ name: 'Moment resistance / as built', x: yrs.map((y) => r(y)), y: hist.map((q) => r(q.capacityRatio)) }, { name: 'Design moment / as built', x: [0, r(horizon)], y: [r(MEd / MRd0), r(MEd / MRd0)], dash: true }, { name: 'Reliability index against depassivation', x: ty.map((y) => r(y)), y: betas.map((b2) => r(b2)) }, { name: 'Target index', x: [0, r(horizon)], y: [inp.betaT, inp.betaT], dash: true }],
          explain: 'Two views of ageing: the reliability index falls from the start as chlorides advance; the strength stays flat until initiation and then declines with bar loss. Each end of life is where a curve meets its dashed line.' },
        { type: 'bar', title: 'Pile wall by corrosion zone', categories: zLife.map((z) => z.name), series: [{ name: 'Today (mm)', y: zLife.map((z) => r(z.remNow * 1e3)) }, { name: 'End of life (mm)', y: zLife.map((z) => r(z.remLife * 1e3)) }, { name: 'Minimum for strength (mm)', y: zLife.map(() => inp.tMin) }], x: { label: 'Zone', unit: '' }, y: { label: 'Remaining wall', unit: 'mm' },
          explain: 'Remaining wall in each exposure zone now and at the end of life against the wall needed for strength. The splash and low-water zones usually govern once the coating has gone.' },
      ],
      tables: [
        { title: 'Service-life limit states', columns: [{ key: 'name', label: 'Limit state' }, { key: 't', label: 'Reached at age', unit: 'years' }, { key: 'rem', label: 'Years from today' }, { key: 'ok', label: 'Meets the required life' }],
          rows: lim.map((q) => ({ name: q.name, t: r(cap(q.t, 9999)), rem: r(cap(Math.max(q.t - age, 0), 9999)), ok: q.t >= life ? 'yes' : 'no' })), explain: 'Each deterioration mechanism with the age at which it ends the service life. The earliest row governs and names what to repair or protect first.' },
        { title: 'Pile corrosion by zone', columns: [{ key: 'name', label: 'Zone' }, { key: 'rate', label: 'Rate', unit: 'mm/yr' }, { key: 'protect', label: 'Protection life', unit: 'years' }, { key: 'loss', label: 'Loss at end of life', unit: 'mm' }, { key: 'rem', label: 'Remaining wall', unit: 'mm' }, { key: 'W', label: 'Section modulus ratio' }, { key: 'end', label: 'Minimum wall reached at', unit: 'years' }],
          rows: zLife.map((z) => ({ name: z.name, rate: z.rateMm, protect: z.protect, loss: r(z.lossLife * 1e3), rem: r(z.remLife * 1e3), W: r(z.modulus, 3), end: r(cap(z.end, 9999)) })), explain: 'Wall loss of the piles zone by zone after the protection has expired. Use the remaining wall as the corrosion deduction in the strength module for a reassessment.' },
      ],
      notes: [`${gov.name} governs the service life at ${f(cap(endOfLife, 9999))} years. Apparent diffusion coefficient today ${f(DaNow * 1e12)} ×10⁻¹² m²/s${Dfit > 0 ? `; the measured profile fits ${f(Dfit * 1e12)} ×10⁻¹² m²/s` : ''}.`,
        `Initiation at ${f(cap(ti, 9999))} years with the mean parameters; the reliability target of ${f(inp.betaT)} is reached at ${f(tRel)} years. A cover of ${f(cNeed * 1e3)} mm would delay initiation to the required life.`],
      advice: [
        ok ? eng('ok', `Estimated service life ${f(cap(endOfLife, 9999))} years against ${f(life)} required (${f(cap(remaining, 9999))} years remain); corrosion initiation at ${f(cap(ti, 9999))} years, least pile wall ${f(Math.min(...zLife.map((z) => z.remLife)) * 1e3)} mm at the end of life.`,
          `Verify the prediction with chloride profiles and a cover survey at the next inspection, and renew the pile coating or anodes before year ${f(Math.min(...zrows.filter((z) => z.protect > 0).map((z) => z.protect), 9999))} when the first protection expires.`, 'fib Bulletin 34 (2006); ISO 13822:2010; EN 1993-5:2007 §4',
          `Durability is the low-carbon option: reaching the full life avoids rebuilding the deck soffit (${f(concRepair)} m³ of repair concrete per metre of beam ≈ ${f(concRepair * CO2_CONC * 1e3)} kg CO₂e) and replacing piles (${f(tPerM)} t of steel per metre ≈ ${f(tPerM * CO2_STEEL)} t CO₂e, screening factors).`)
          : eng('fail', `Estimated service life ${f(cap(endOfLife, 9999))} years against ${f(life)} required: ${gov.name.toLowerCase()}. ${w.name} is the governing check at ${f(w.uc)}.`,
            gov.name.startsWith('Cover') || gov.name.startsWith('Deck') ? `Before initiation (year ${f(cap(ti, 9999))}): apply a surface treatment or coating to slow ingress. After it: install cathodic protection of the reinforcement, or break out and repair with ${f(Math.max(cNeed * 1e3, inp.cover))} mm of low-permeability cover; strengthen in flexure if the resistance has fallen below ${f(MEd / 1e3)} kN·m.` : gov.name.startsWith('Pile') ? `Renew the coating or fit sacrificial anodes before year ${f(govZ.protect)}; where the wall is already below ${f(inp.tMin)} mm, sleeve the piles through the ${govZ.name.toLowerCase()} with welded plates or a grouted jacket.` : 'Grind out and re-weld the cracked detail, then re-inspect; reduce the stress range if the crack recurs.',
            'fib Bulletin 34 (2006); EN 1504 (repair of concrete); EN ISO 12696 (cathodic protection of steel in concrete); BS 6349-1-4:2021',
            `Cathodic protection or timely coating renewal keeps the existing material in service: replacing instead costs about ${f(tPerM)} t of steel per metre of pile (${f(tPerM * CO2_STEEL)} t CO₂e) and ${f(concRepair)} m³ of concrete per metre of beam (${f(concRepair * CO2_CONC * 1e3)} kg CO₂e, screening factors).`),
        eco(ok ? 'info' : 'warn', `Without intervention ${f(cap(remaining, 9999))} years of service remain; the maintenance-free period of the deck ends at about year ${f(cap(Math.min(tRel, tCrack), 9999))} and pile protection expires from year ${f(Math.min(...zrows.filter((z) => z.protect > 0).map((z) => z.protect), 9999))}.`,
          'Plan interventions on whole-life cost: preventive measures before corrosion starts (coating, anodes) cost a fraction of patch repair after cracking, and repair costs roughly quintuple with each stage of delay. Set inspection intervals from the predicted dates and update the model with each survey.', `${ECON_BASIS}; de Sitter's law of fives; risk-based inspection per DNV-RP-G101`),
      ],
      publish: { 'life.remainingYears': Math.min(remaining, 9999), 'wharf.initiationYears': Math.min(ti, 9999), 'wharf.pileWallRemaining': Math.min(...zLife.map((z) => z.remLife)) },
    };
  },
  verify: [
    // Constant diffusion coefficient (m = 0): c/(2√(D t)) = erfc⁻¹[(Ccrit − Ci)/(Cs − Ci)]. Ccrit/Cs = 0.6/3.0 = 0.2 → erfc⁻¹(0.2) = 0.9062;
    // c = 60 mm, D = 1e-12 m²/s: t = (0.06/(2·0.9062))²/1e-12 = 1.0960e9 s = 34.73 years.
    // At 54.73 years (1.7271e9 s): C = 3·erfc[0.06/(2·√(1.7271e-3))] = 3·erfc(0.7219) = 3·0.3073 = 0.922 %.
    { name: "Fick's second law: time to corrosion initiation", source: 'Crank (1975) error-function solution; Collepardi et al. (1972) — hand calculation',
      inputs: { cover: 60, Cs: 3, Ccrit: 0.6, Ci: 0, D0: 1, mAge: 0, dxConv: 0, age: 10, life: 54.73, icorr: 1, pit: 1 }, expect: [{ kpi: 'ti', value: 34.73, tol: 0.003 }, { kpi: 'Clife', value: 0.922, tol: 0.005 }] },
    // Faraday: 1 µA/cm² = 0.01 A/m² → 0.01·0.055845/(2·96485·7850) m/s = 3.687e-13 m/s = 11.63 µm/yr.
    // Twenty years after initiation (age 54.73): penetration 0.2327 mm, d = 32 − 0.465 = 31.535 mm, area lost = 1 − (31.535/32)² = 2.89 %.
    { name: "Faraday's law: 1 µA/cm² = 11.6 µm/yr, and the bar section lost after 20 years", source: "Faraday's law for iron (M = 55.845 g/mol, z = 2, ρ = 7850 kg/m³) — hand calculation",
      inputs: { cover: 60, Cs: 3, Ccrit: 0.6, Ci: 0, D0: 1, mAge: 0, dxConv: 0, age: 10, life: 54.73, icorr: 1, pit: 1, dBot: 32 }, expect: [{ kpi: 'corrRate', value: 11.63, tol: 0.003 }, { kpi: 'lossAtLife', value: 2.89, tol: 0.02 }] },
    // Splash zone at 0.3 mm/yr after 10 years of coating life, 50-year life: loss = 0.3·40 = 12 mm; 25 mm wall → 13 mm remains;
    // the 15 mm minimum is reached after 10 + (25 − 15)/0.3 = 43.3 years.
    { name: 'Pile wall loss after the coating life: loss = rate × (life − protection life)', source: 'EN 1993-5 §4.4 method (hand calculation)',
      inputs: { life: 50, t: 25, tMin: 15, zones: [{ name: 'Splash zone', rate: 0.3, protect: 10, faces: 1 }] }, expect: [{ kpi: 'tRemain', value: 13, tol: 1e-6 }, { kpi: 'tPile', value: 43.33, tol: 1e-3 }] },
  ],
};
