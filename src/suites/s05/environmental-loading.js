// Wave, current, wind and tidal loading on one bent of a piled wharf or jetty: Morison loading on the pile row
// with the phase lag between rows, the crest against the deck soffit (air gap), wave-in-deck and uplift screening,
// wind on the deck equipment, and the wave height at which the deck is first reached.
import { airyKinematics, stokes5Wave, stokes5Kinematics, selectWaveTheory, wavelength } from '../../engine/metocean/waves.js';
import { streamFunctionWave, streamFunctionKinematics } from '../../engine/metocean/streamfunction.js';
import { currentPowerLaw } from '../../engine/metocean/current.js';
import { windForce, heightCoefficient } from '../../engine/metocean/wind.js';
import { pileForceClosedForm, keuleganCarpenter } from '../../engine/hydro/morison.js';
import { waveInDeckForce } from '../../engine/hydro/impact.js';
import { G, need, run, r, f, check, eng, eco, linspace, maxWhere, CONFIRM, NOT_HIFI, ECON_BASIS, CO2_CONC } from './_common.js';

/** Wave solution of the selected theory with one kinematics call signature: (x, z, t) → { eta, u, ax }. */
function waveModel(theory, { H, T, d, rho }) {
  if (theory === 'stokes5') {
    const sol = run(() => stokes5Wave({ H, d, T }), 'Stokes fifth-order wave');
    return { L: sol.L, valid: sol.valid, at: (x, z, t) => { const e = stokes5Kinematics(sol, { x, z: 0, t }).eta; const k = stokes5Kinematics(sol, { x, z: Math.min(z, e - 1e-9), t, rho }); return { eta: e, u: k.u, ax: k.ax }; } };
  }
  if (theory === 'stream') {
    const sol = run(() => streamFunctionWave({ H, d, T, order: 10 }), 'Stream-function wave');
    need(sol.report.converged, 'The stream-function wave did not converge: the wave is too close to breaking for this depth and period.');
    return { L: sol.L, valid: true, at: (x, z, t) => { const k0 = streamFunctionKinematics(sol, { x, z: -d, t, rho }); const k = streamFunctionKinematics(sol, { x, z: Math.min(z, k0.eta - 1e-9), t, rho }); return { eta: k0.eta, u: k.u, ax: k.ax }; } };
  }
  // linear theory with Wheeler stretching so that the crest kinematics stay bounded
  return { L: wavelength(T, d), valid: true, at: (x, z, t) => { const k = airyKinematics({ H, T, d, x, z, t, stretching: 'wheeler', clipAboveSurface: false, rho }); return { eta: k.eta, u: k.u, ax: k.ax }; } };
}

export default {
  id: 's05-environmental-loading', suite: 5, title: 'Wave, current, wind and tidal loading',
  summary: 'Wave and current force on the pile rows of one bent with their phase lag, crest level against the deck soffit, wave-in-deck and uplift screening, wind on deck equipment, and the limiting wave height for the air gap.',
  covers: [4], fidelity: 'analytical',
  standards: ['BS 6349-1-2:2016 (assessment of actions)', 'DNV-RP-C205 (2019) §3, §6', 'ISO 21650:2007 (actions from waves and currents on coastal structures)', 'API RP 2A-WSD 22nd ed. §5.3', 'ISO 19902:2020 §A.9.4 (air gap)'],
  models: ['Airy wave theory with Wheeler stretching', 'Stokes fifth-order wave theory (Fenton)', 'Stream-function wave theory', "Morison's equation", 'Power-law tidal current profile', 'Wind pressure with height coefficient', 'Wave-in-deck silhouette drag', 'Hydrostatic uplift on a submerged soffit'],
  limits: `Slender-pile Morison loading (diameter below one fifth of the wavelength) with regular-wave kinematics evaluated for the depth at the berthing face; no diffraction, no pile-group interference beyond the shielding factor entered. Wave-in-deck force and uplift are quasi-static screening values without slam dynamics. ${NOT_HIFI}`,
  publishes: ['wharf.envForceBent', 'wharf.windForceBent', 'wharf.waveCrestLevel'],
  searchTerms: 'wave load piled jetty Morison air gap wave-in-deck uplift',
  inputs: [
    { key: 'depth', label: 'Water depth at the berthing face (below datum)', unit: 'm', type: 'number', default: 12, min: 2, max: 40, group: 'Water levels', role: 'bc', bind: 'site.depth' },
    { key: 'tide', label: 'High-water level above datum', unit: 'm', type: 'number', default: 1.2, min: -3, max: 8, group: 'Water levels', role: 'bc', help: 'Design high water (for example HAT or MHWS) relative to the datum of the levels below. Use the low-water level to check the longest wetted pile.' },
    { key: 'surge', label: 'Storm surge', unit: 'm', type: 'number', default: 0.6, min: 0, max: 5, group: 'Water levels', role: 'bc', bind: 'site.surge100' },
    { key: 'soffit', label: 'Deck soffit level above datum', unit: 'm', type: 'number', default: 3.3, min: 0.5, max: 15, group: 'Water levels', role: 'geometry', design: true, help: 'Underside of the lowest deck beam.' },
    { key: 'deckTop', label: 'Deck surface level above datum', unit: 'm', type: 'number', default: 4.7, min: 1, max: 18, group: 'Water levels', role: 'geometry' },
    { key: 'Hmax', label: 'Design maximum wave height at the berth', unit: 'm', type: 'number', default: 2.2, min: 0, max: 15, group: 'Waves and current', role: 'environment', uncertain: { dist: 'gumbel', cov: 0.15 }, help: 'Individual maximum wave inside the harbour after diffraction and refraction (about 1.8 to 2 Hs). Zero for current alone.' },
    { key: 'T', label: 'Associated wave period', unit: 's', type: 'number', default: 7, min: 2, max: 25, group: 'Waves and current', role: 'environment', bind: 'site.Tp100' },
    { key: 'Uc', label: 'Tidal current speed at the surface', unit: 'm/s', type: 'number', default: 1.0, min: 0, max: 4, group: 'Waves and current', role: 'environment', bind: 'site.Uc100', uncertain: { dist: 'normal', cov: 0.15 }, help: 'In line with the wave (conservative); a one-seventh power profile is applied over the depth.' },
    { key: 'theory', label: 'Wave theory', type: 'select', default: 'auto', group: 'Waves and current', role: 'control', options: [{ value: 'auto', label: 'Automatic from the validity chart' }, { value: 'airy', label: 'Airy with Wheeler stretching' }, { value: 'stokes5', label: 'Stokes fifth order' }, { value: 'stream', label: 'Stream function (order 10)' }] },
    { key: 'rhoW', label: 'Water density', unit: 'kg/m³', type: 'number', default: 1025, min: 995, max: 1035, group: 'Waves and current', role: 'environment', bind: 'site.rhoW' },
    { key: 'nPiles', label: 'Pile rows in the bent', unit: '-', type: 'int', default: 4, min: 1, max: 12, group: 'Piles', role: 'geometry' },
    { key: 'sx', label: 'Row spacing across the wharf', unit: 'm', type: 'number', default: 6, min: 1.5, max: 15, group: 'Piles', role: 'geometry' },
    { key: 'sy', label: 'Bent spacing along the wharf', unit: 'm', type: 'number', default: 6, min: 1.5, max: 15, group: 'Piles', role: 'geometry', help: 'Tributary deck length of one bent for wave-in-deck, uplift and wind.' },
    { key: 'bedRise', label: 'Rise of the seabed to the landward row', unit: 'm', type: 'number', default: 6, min: 0, max: 30, group: 'Piles', role: 'bc' },
    { key: 'D', label: 'Pile outer diameter', unit: 'm', type: 'number', default: 1.067, min: 0.2, max: 3, group: 'Piles', role: 'geometry', design: true, bind: 'wharf.pileD' },
    { key: 'growth', label: 'Marine growth thickness', unit: 'mm', type: 'number', default: 50, min: 0, max: 200, group: 'Piles', role: 'ic', help: 'Added to the radius for the hydrodynamic diameter.' },
    { key: 'Cd', label: 'Drag coefficient', unit: '-', type: 'number', default: 1.05, min: 0.5, max: 2, group: 'Piles', role: 'criteria', help: `1.05 rough, 0.65 smooth (DNV-RP-C205 §6.7; API RP 2A). ${CONFIRM}` },
    { key: 'Cm', label: 'Inertia coefficient', unit: '-', type: 'number', default: 1.6, min: 1, max: 2.2, group: 'Piles', role: 'criteria', help: `1.2 rough, 1.6 smooth by API RP 2A; 2.0 potential-flow value. ${CONFIRM}` },
    { key: 'shield', label: 'Shielding / interference factor on the row total', unit: '-', type: 'number', default: 1, min: 0.5, max: 1.5, group: 'Piles', role: 'criteria', help: 'Below 1 only with model-test evidence; above 1 for closely spaced piles (spacing under about 4 D).' },
    { key: 'CdDeck', label: 'Deck silhouette drag coefficient', unit: '-', type: 'number', default: 2.5, min: 1, max: 4, group: 'Deck', role: 'criteria', help: `API RP 2A-WSD §C17.6.2: 2.0 to 2.5 for a heavily equipped or solid deck edge. ${CONFIRM}` },
    { key: 'deckWidth', label: 'Deck width across the wharf', unit: 'm', type: 'number', default: 21, min: 2, max: 80, group: 'Deck', role: 'geometry' },
    { key: 'Uw', label: 'Wind speed (3 s gust at 10 m)', unit: 'm/s', type: 'number', default: 45, min: 0, max: 90, group: 'Wind', role: 'environment', uncertain: { dist: 'gumbel', cov: 0.1 }, help: 'Out-of-service gust for the stowed crane and stacked cargo.' },
    { key: 'alphaW', label: 'Wind profile exponent', unit: '-', type: 'number', default: 0.12, min: 0.07, max: 0.3, group: 'Wind', role: 'environment' },
    { key: 'wind', label: 'Wind-exposed items carried by one bent', type: 'table', group: 'Wind', role: 'bc',
      columns: [{ key: 'name', label: 'Item', type: 'text' }, { key: 'A', label: 'Projected area per bent', unit: 'm²', type: 'number' }, { key: 'Cs', label: 'Shape coefficient', unit: '-', type: 'number' }, { key: 'z', label: 'Centroid above datum', unit: 'm', type: 'number' }],
      default: [{ name: 'Crane (stowed), share per bent', A: 60, Cs: 1.6, z: 30 }, { name: 'Container stack, 4 high', A: 62, Cs: 1.3, z: 10 }, { name: 'Deck edge and services', A: 9, Cs: 1.5, z: 4 }],
      help: 'Projected area facing the wind that one bent carries. Crane wind reactions from the supplier can instead be entered as rail loads in the strength module.' },
    { key: 'nSeg', label: 'Integration segments over the water column', unit: '-', type: 'int', default: 40, min: 5, max: 640, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nSeg', label: 'Segments over the water column', levels: [10, 20, 40, 80], qoi: 'Fwave' },
  compute(inp) {
    const d = inp.depth + inp.tide + inp.surge, swl = inp.tide + inp.surge, rho = inp.rhoW, De = inp.D + 2 * inp.growth / 1e3;
    need(d > 1.5, 'The water depth at the berthing face including tide and surge must exceed 1.5 m.');
    need(inp.deckTop > inp.soffit, 'The deck surface must be above the soffit.');
    const n = Math.round(inp.nPiles), nSeg = Math.round(inp.nSeg), H = Math.max(inp.Hmax, 1e-6), hasWave = inp.Hmax > 0;
    const rec = hasWave ? selectWaveTheory({ H, T: inp.T, d }) : { theory: 'airy', HoverHb: 0, ursell: 0, dOverL: 0, steepness: 0 };
    need(rec.theory !== 'breaking', `A ${f(inp.Hmax)} m, ${f(inp.T)} s wave breaks in ${f(d)} m of water (H/Hb = ${f(rec.HoverHb)}): use the depth-limited height of about ${f(inp.Hmax / rec.HoverHb)} m, or a breaking-wave (slam) assessment.`);
    let theory = inp.theory === 'auto' ? ({ airy: 'airy', stokes5: 'stokes5' }[rec.theory] ?? 'stream') : inp.theory, note = '';
    let wm;
    try { wm = waveModel(theory, { H, T: inp.T, d, rho }); if (!wm.valid) throw new Error('outside validity'); } catch (e) {
      if (inp.theory !== 'auto' && theory !== 'stokes5') throw e;
      note = `${theory === 'stokes5' ? 'Stokes fifth-order' : 'Stream-function'} theory is not usable for this wave (${e.message}); ${theory === 'stokes5' ? 'the stream-function solution' : 'Airy theory with Wheeler stretching'} is used instead.`;
      theory = theory === 'stokes5' ? 'stream' : 'airy';
      try { wm = waveModel(theory, { H, T: inp.T, d, rho }); } catch { theory = 'airy'; wm = waveModel('airy', { H, T: inp.T, d, rho }); note = 'Nonlinear wave solutions failed for this wave; Airy theory with Wheeler stretching is used.'; }
    }
    // pile rows: x landward from the seaward row, local bed level rising linearly (kinematics of the front depth)
    const piles = Array.from({ length: n }, (_, i) => ({ id: i + 1, x: i * inp.sx, dLoc: Math.max(d - (n > 1 ? inp.bedRise * i / (n - 1) : 0), 0) }));
    const nPh = hasWave ? 36 : 1, phases = [], Ftot = [], Mtot = [], etaFront = [];
    const perPile = piles.map(() => ({ Fmax: 0, Mmax: 0 }));
    let uMax = 0, crest = 0;
    for (let ip = 0; ip < nPh; ip++) {
      const t = inp.T * ip / nPh; let F = 0, M = 0;
      piles.forEach((p, j) => {
        if (p.dLoc <= 0.05) return;
        // waves travel landward (+x); phase lag kx between rows
        const eta = wm.at(p.x, 0, t).eta, top = eta, bot = -p.dLoc, dz = (top - bot) / nSeg;
        let Fp = 0, Mp = 0;
        for (let i = 0; i < nSeg; i++) {
          const z = bot + (i + 0.5) * dz, k = wm.at(p.x, z, t);
          const zc = Math.min(0, (z - eta) * d / (d + eta));  // current profile stretched with the surface
          const u = k.u + (inp.Uc > 0 ? currentPowerLaw(Math.max(zc, -d), { v0: inp.Uc, d }) : 0);
          const q = 0.5 * rho * inp.Cd * De * u * Math.abs(u) + rho * inp.Cm * Math.PI * De * De / 4 * k.ax;   // Morison, N/m
          Fp += q * dz; Mp += q * dz * (z - bot);
          if (Math.abs(u) > uMax) uMax = Math.abs(u);
        }
        F += Fp; M += Mp;
        if (Math.abs(Fp) > perPile[j].Fmax) perPile[j].Fmax = Math.abs(Fp);
        if (Math.abs(Mp) > perPile[j].Mmax) perPile[j].Mmax = Math.abs(Mp);
        if (j === 0) { etaFront.push(eta); crest = Math.max(crest, eta); }
      });
      phases.push(360 * ip / nPh); Ftot.push(inp.shield * F); Mtot.push(inp.shield * M);
    }
    let iMax = 0; Ftot.forEach((v, i) => { if (Math.abs(v) > Math.abs(Ftot[iMax])) iMax = i; });
    const Fwave = Math.abs(Ftot[iMax]), Mwave = Math.abs(Mtot[iMax]), sumPeaks = inp.shield * perPile.reduce((s, p) => s + p.Fmax, 0);
    if (!hasWave) crest = 0;
    const cf = hasWave ? pileForceClosedForm({ H, T: inp.T, depth: d, D: De, Cd: inp.Cd, Cm: inp.Cm, rho }) : { Fmax: 0, FD: 0, FI: 0 };
    const Um = Math.max(uMax - inp.Uc, 0), KC = hasWave && Um > 0 ? keuleganCarpenter(Um, inp.T, De) : 0;
    // air gap and deck loads
    const crestLevel = swl + crest, gap = inp.soffit - crestLevel;
    const uCrest = hasWave ? wm.at(0, crest - 1e-6, 0).u : 0;
    const wid = waveInDeckForce({ crestElevation: crestLevel, deckBottom: inp.soffit, deckTop: inp.deckTop, width: inp.sy, u: uCrest, U: inp.Uc, Cd: inp.CdDeck, rho });
    // quasi-static uplift: head of the crest above the soffit integrated over the wetted length of soffit
    let upl = 0, wetLen = 0;
    if (gap < 0 && hasWave) { const xs = linspace(-wm.L / 2, wm.L / 2, 121), dx = wm.L / 120; for (const x of xs) { const h = swl + wm.at(x, 0, 0).eta - inp.soffit; if (h > 0 && wetLen < inp.deckWidth) { upl += rho * G * h * dx * inp.sy; wetLen += dx; } } }
    const deckSelf = 25e3 * 0.6 * inp.sy * Math.max(wetLen, 1e-9); // indicative 0.6 m equivalent concrete thickness resisting uplift
    // limiting wave height for a dry soffit (crest scales with the same theory)
    const crestOf = (h) => { if (h <= 0) return 0; try { const s = selectWaveTheory({ H: h, T: inp.T, d }); if (s.theory === 'breaking') return 1e9; const m = waveModel(theory, { H: h, T: inp.T, d, rho }); return m.at(0, 0, 0).eta; } catch { return 1e9; } };
    const Hlim = inp.soffit > swl ? maxWhere((h) => swl + crestOf(h) <= inp.soffit, 0.01, Math.min(0.78 * d, 15), 12) : 0;
    // wind on deck items
    const items = (inp.wind || []).map((q, i) => { need(q.A >= 0 && q.Cs >= 0 && q.z > 0, `Wind item ${i + 1}: area and shape coefficient must not be negative and the centroid must be above datum.`); const w = windForce({ U: inp.Uw, A: Number(q.A), Cs: Number(q.Cs), Ch: heightCoefficient(Math.max(Number(q.z), 1), 10, inp.alphaW) }); return { name: String(q.name || `Item ${i + 1}`), F: w.F, z: Number(q.z), A: Number(q.A) }; });
    const Fwind = items.reduce((s, q) => s + q.F, 0), MwindDeck = items.reduce((s, q) => s + q.F * (q.z - inp.deckTop), 0);
    const Fenv = Fwave + wid.force, lever = Fwave > 0 ? Mwave / Fwave : 0;
    const DoverL = hasWave ? De / wm.L : 0;
    const checks = [
      check('Wave crest against the deck soffit (air gap)', Math.max(crestLevel, 0), inp.soffit, 'm', 'ISO 19902:2020 §A.9.4; BS 6349-1-2:2016 §27', 'Crest level including tide and surge against the underside of the deck. Once the crest reaches the soffit the deck takes slam, drag and uplift that a pile-only calculation misses.'),
      check('Slender-member validity (D / wavelength)', DoverL, 0.2, '-', 'DNV-RP-C205 (2019) §6.1.2', 'Morison loading assumes the pile does not disturb the wave. Above 0.2 diffraction governs and the inertia force is over-estimated.'),
      check('Wave height against the breaking limit', rec.HoverHb, 1, '-', 'DNV-RP-C205 (2019) §3.4.6 (Miche and depth limits)', 'A wave near breaking has steep, fast crest kinematics and may slam on the piles and deck; regular-wave theory then under-predicts the peak load.'),
      ...(gap < 0 ? [check('Uplift on the wetted soffit against deck weight', upl / 1e3, deckSelf / 1e3, 'kN', 'BS 6349-1-2:2016 §28 (wave uplift, quasi-static screening)', 'Upward wave pressure on the underside against an indicative deck weight over the same length. If exceeded, slab panels can lift off their bearings unless they are tied down.')] : []),
    ];
    const w = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = w.uc <= 1;
    const dRaise = Math.max(0.5 - gap, 0), concRaise = dRaise * 0.9 * n; // pile-head/plinth concrete per bent if raised, m³ (indicative 0.9 m² plinth)
    const zs = linspace(-d, 0, 21);
    return {
      kpis: [
        { key: 'Fwave', label: 'Wave and current force on the pile rows of one bent', value: r(Fwave / 1e3, 6), unit: 'kN', explain: 'Largest simultaneous horizontal force on all pile rows during the wave cycle, with the phase lag between rows and the shielding factor. This is the environmental load on the bent.' },
        { key: 'Mwave', label: 'Moment of that force about the local seabed', value: r(Mwave / 1e3), unit: 'kN·m', explain: 'Sum over the rows of each pile force times its height above the local bed, at the phase of maximum force. Divide by the force for the lever arm.' },
        { key: 'lever', label: 'Lever arm above the seabed', value: r(lever), unit: 'm', explain: 'Height of the resultant above the bed. A high resultant (near the surface) is typical of drag-dominated loading with current.' },
        { key: 'FwaveSum', label: 'Sum of the individual row peaks', value: r(sumPeaks / 1e3), unit: 'kN', explain: 'Upper bound that ignores the phase lag between rows. The difference from the simultaneous force is the benefit of the rows not peaking together.' },
        { key: 'FwaveLin', label: 'Single pile, linear closed form (no current)', value: r(cf.Fmax / 1e3), unit: 'kN', explain: 'Morison force on one pile at the berthing face from linear theory integrated to still water: an independent check on the numerical integration.' },
        { key: 'KC', label: 'Keulegan–Carpenter number', value: r(KC), unit: '-', explain: 'Wave orbital excursion relative to the pile diameter. Below about 6 inertia dominates; above about 20 drag dominates and the force follows velocity squared.' },
        { key: 'crest', label: 'Crest level above datum', value: r(crestLevel), unit: 'm', limit: inp.soffit, sense: 'max', explain: 'Tide plus surge plus crest height from the wave theory used. It must stay below the soffit unless the deck is designed for wave loading.' },
        { key: 'airGap', label: 'Air gap below the soffit', value: r(gap), unit: 'm', explain: 'Clearance between crest and deck underside. Negative means the deck is inundated; a margin of at least 0.5 to 1.5 m is normal practice.' },
        { key: 'Hlim', label: 'Maximum wave height for a dry soffit', value: r(Hlim), unit: 'm', explain: 'Largest wave of this period at this water level whose crest just reaches the soffit: an operating and design limit for the deck level.' },
        { key: 'Fdeck', label: 'Wave-in-deck horizontal force on one bent', value: r(wid.force / 1e3), unit: 'kN', explain: 'Drag on the inundated deck silhouette over one bent spacing. Zero while the crest is below the soffit.' },
        { key: 'uplift', label: 'Quasi-static uplift on the soffit of one bent', value: r(upl / 1e3), unit: 'kN', explain: 'Hydrostatic head of the crest above the soffit integrated over the wetted length. Slam pressures are several times higher and need model tests or specialist analysis.' },
        { key: 'Fwind', label: 'Wind force on the items carried by one bent', value: r(Fwind / 1e3), unit: 'kN', explain: 'Sum of ½ρU²·Cs·Ch·A over the items entered, with the gust speed referred to each centroid height.' },
        { key: 'Fenv', label: 'Wave, current and deck force passed to the bent', value: r(Fenv / 1e3), unit: 'kN', explain: 'Pile-row force plus wave-in-deck force: the value published for the strength analysis of the bent.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Force on the bent through the wave cycle', x: { label: 'Phase', unit: 'deg' }, y: { label: 'Force', unit: 'kN' },
          series: [{ name: 'All rows', x: phases, y: Ftot.map((v) => r(v / 1e3)) }, ...(hasWave ? [{ name: 'Front-row surface elevation × 100 (m)', x: phases, y: etaFront.map((v) => r(v * 100)), dash: true }] : [])],
          explain: 'Total horizontal force against wave phase. Inertia-dominated loading peaks a quarter cycle before the crest; drag-dominated loading peaks at the crest.' },
        { type: 'line', title: 'Horizontal velocity under the crest', x: { label: 'Velocity', unit: 'm/s' }, y: { label: 'Elevation above still water', unit: 'm' },
          series: [{ name: 'Wave + current', x: zs.map((z) => r(wm.at(0, z * (1 - 1e-9), 0).u * (hasWave ? 1 : 0) + (inp.Uc > 0 ? currentPowerLaw(z, { v0: inp.Uc, d }) : 0))), y: zs.map((z) => r(z)) }, { name: 'Current', x: zs.map((z) => r(inp.Uc > 0 ? currentPowerLaw(z, { v0: inp.Uc, d }) : 0)), y: zs.map((z) => r(z)), dash: true }],
          explain: 'Velocity profile below still water at the crest phase. The current raises the whole profile and, because drag goes with velocity squared, adds more force than its own drag alone.' },
        { type: 'bar', title: 'Peak force per pile row', categories: piles.map((p) => `Row ${p.id}`), series: [{ name: 'Peak force (kN)', y: perPile.map((p) => r(inp.shield * p.Fmax / 1e3)) }], x: { label: 'Row (seaward first)', unit: '' }, y: { label: 'Force', unit: 'kN' },
          explain: 'Largest force on each row during the cycle. Landward rows stand in shallower water on the under-wharf slope and take less.' },
      ],
      tables: [
        { title: 'Pile rows', columns: [{ key: 'id', label: 'Row' }, { key: 'x', label: 'x', unit: 'm' }, { key: 'd', label: 'Local depth', unit: 'm' }, { key: 'F', label: 'Peak force', unit: 'kN' }, { key: 'M', label: 'Peak moment about the bed', unit: 'kN·m' }],
          rows: piles.map((p, j) => ({ id: p.id, x: r(p.x), d: r(p.dLoc), F: r(inp.shield * perPile[j].Fmax / 1e3), M: r(inp.shield * perPile[j].Mmax / 1e3) })), explain: 'Each row with its water depth at the design level and its own peak force and moment; apply these as line loads over the wetted length in a detailed frame model.' },
        { title: 'Wind items', columns: [{ key: 'name', label: 'Item' }, { key: 'A', label: 'Area', unit: 'm²' }, { key: 'z', label: 'Centroid', unit: 'm' }, { key: 'F', label: 'Force', unit: 'kN' }],
          rows: items.map((q) => ({ name: q.name, A: r(q.A), z: r(q.z), F: r(q.F / 1e3) })), explain: 'Wind force on each item at its centroid height. The moment about deck level is the sum of force times height above the deck.' },
      ],
      notes: [`${theory === 'airy' ? 'Airy theory with Wheeler stretching' : theory === 'stokes5' ? 'Stokes fifth-order theory' : 'Stream-function theory (order 10)'} at ${f(d)} m depth (H/Hb = ${f(rec.HoverHb)}, Ursell number ${f(rec.ursell)}); hydrodynamic diameter ${f(De)} m.${note ? ' ' + note : ''}`,
        `Wind moment about deck level ${f(MwindDeck / 1e3)} kN·m; the simultaneous pile force is ${f(sumPeaks > 0 ? Fwave / sumPeaks * 100 : 100)} % of the sum of row peaks.`],
      advice: [
        ok ? eng('ok', `Environmental force on the bent ${f(Fenv / 1e3)} kN with a lever arm of ${f(lever)} m; crest at ${f(crestLevel)} m leaves an air gap of ${f(gap)} m; wind adds ${f(Fwind / 1e3)} kN at deck level and above.`,
          `Carry ${f(Fenv / 1e3)} kN into the bent strength check. ${gap < 0.5 ? `The air gap is below 0.5 m: allow for sea-level rise and wave setup or raise the soffit by ${f(dRaise)} m.` : `The soffit stays dry up to a ${f(Hlim)} m wave at this water level.`} Repeat at low water for the longest wetted piles.`, 'BS 6349-1-2:2016; DNV-RP-C205 (2019) §6; ISO 19902:2020 §A.9.4',
          `Marine growth of ${f(inp.growth)} mm accounts for ${f((1 - inp.D / De) * 100)} % of the drag area: cleaning in the splash zone lowers the load without new material.`)
          : eng('fail', `${w.name}: ${f(w.demand)} against ${f(w.capacity)} ${w.unit}. Crest at ${f(crestLevel)} m, soffit at ${f(inp.soffit)} m; wave-in-deck force ${f(wid.force / 1e3)} kN and uplift ${f(upl / 1e3)} kN per bent.`,
            gap < 0 ? `Raise the soffit by at least ${f(dRaise)} m, or design the deck for wave loading (tie-downs for ${f(upl / 1e3)} kN uplift per bent, venting, and physical model tests for slam).` : 'Use a diffraction solution for the piles or a breaking-wave assessment as the failed check indicates.', 'BS 6349-1-2:2016 §27–28; ISO 19902:2020 §A.9.4',
            `Raising the deck by ${f(dRaise)} m needs about ${f(concRaise)} m³ more pile-head concrete per bent (${f(concRaise * CO2_CONC)} t CO₂e, screening factor) but avoids recurring storm damage to soffit, services and fender fixings.`),
        eco('info', `The deck stays clear of waves up to ${f(Hlim)} m at the design water level; the design wave is ${f(inp.Hmax)} m.`,
          'Set the deck level from the whole-life balance: each 0.1 m of freeboard costs longer piles and ramps once, while an inundated soffit costs repairs and berth closure after every storm that exceeds the limit. Check the return period of the limiting wave with sea-level rise over the design life.', `${ECON_BASIS}; PIANC WG 178 climate change adaptation planning`),
      ],
      publish: { 'wharf.envForceBent': Fenv, 'wharf.windForceBent': Fwind, 'wharf.waveCrestLevel': crestLevel },
    };
  },
  verify: [
    // Inertia-dominated pile in deep water: k = ω²/g = 0.2516 1/m, tanh(kd) = 1;
    // F_I = ρ·Cm·(πD²/4)·g·(H/2) = 1025·2·3.1416·9.80665·0.5 = 31.58 kN; drag ≈ 1.3 kN is out of phase, so Fmax = F_I.
    { name: 'Inertia force on a vertical pile, deep water: F = ρ·Cm·(πD²/4)·g·H/2', source: 'Morison et al. (1950) with Airy kinematics; Dean & Dalrymple (1991) §8',
      inputs: { nPiles: 1, D: 2, growth: 0, shield: 1, depth: 38, tide: 0, surge: 0, bedRise: 0, Hmax: 1, T: 4, Uc: 0, Cd: 1, Cm: 2, theory: 'airy', nSeg: 160, soffit: 10, deckTop: 12 },
      expect: [{ kpi: 'Fwave', value: 31.58, tol: 0.01 }, { kpi: 'FwaveLin', value: 31.58, tol: 0.002 }] },
    // Steady current only: F = ½ρCdD∫u²dz with u = U(1 + z/d)^(1/7): ∫ = U²·d·7/9 → 0.5·1025·1·1·1²·20·7/9 = 7.972 kN;
    // lever arm ∫u² s ds / ∫u² ds = d·(7/16)/(7/9) = 0.5625 d = 11.25 m.
    { name: 'Current drag with a 1/7 power profile: F = ½ρCdD·U²·d·7/9', source: 'Closed-form integral of the power-law profile (DNV-RP-C205 §4.1.4)',
      inputs: { nPiles: 1, D: 1, growth: 0, shield: 1, depth: 20, tide: 0, surge: 0, bedRise: 0, Hmax: 0, T: 8, Uc: 1, Cd: 1, Cm: 2, nSeg: 160, soffit: 10, deckTop: 12 },
      expect: [{ kpi: 'Fwave', value: 7.972, tol: 0.005 }, { kpi: 'lever', value: 11.25, tol: 0.005 }] },
    // Wind on one item at 10 m: F = ½·1.225·40²·1.5·100 = 147.0 kN.
    { name: 'Wind force F = ½ρU²·Cs·A at the reference height', source: 'API RP 2A-WSD §5.3.2.2 Eq. 5.6 (hand calculation)',
      inputs: { Uw: 40, wind: [{ name: 'Panel', A: 100, Cs: 1.5, z: 10 }] }, expect: [{ kpi: 'Fwind', value: 147.0, tol: 0.002 }] },
  ],
};
