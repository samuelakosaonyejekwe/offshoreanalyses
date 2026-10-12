// Hydrostatic pressure, uplift, seepage and flotation of a submerged or buried marine chamber (intake pump
// chamber, outfall shaft, dry structure behind a seawall): Archimedes buoyancy, confined seepage under the base
// with a cut-off wall, exit gradient and piping / heave safety.
import { flotationCheck } from '../../engine/coastal/caisson.js';
import { confinedSeepage2D, criticalGradient } from '../../engine/geotech/seepage.js';
import { G, CO2_CONC, need, run, finite, r, f, linspace, check, sfCheck, eng, eco, slrInputs, slrRows, SLR_BASIS } from './_common.js';

export default {
  id: 's07-uplift-flotation', suite: 7, title: 'Hydrostatic uplift, seepage and flotation',
  summary: 'Buoyancy safety of an emptied chamber, water pressure on its walls and base, and the seepage, uplift and piping risk under the structure with a cut-off wall.',
  covers: [10], fidelity: 'numerical',
  standards: ['EN 1997-1:2004 §10 (hydraulic failure: UPL, HYD)', 'USACE EM 1110-2-2100 (2005) §3-10 flotation', 'USACE EM 1110-2-1901 Seepage Analysis and Control', 'BS 6349-2:2019'],
  models: ['Hydrostatic pressure', 'Archimedes buoyancy principle', 'Darcy seepage (2-D finite volumes)', 'Continuity (Laplace equation for head)', 'Terzaghi heave prism', 'Critical exit gradient'],
  limits: 'Steady confined seepage in a homogeneous (optionally anisotropic) layer over an impervious base in a vertical 2-D section; transient tidal seepage, 3-D flow round the ends, soil layering and wall friction along the chamber sides (unless entered as resisting force) are not modelled. The exit gradient next to a cut-off converges slowly: refine the grid before relying on it.',
  publishes: ['marine.upliftForce', 'marine.flotationSF'],
  searchTerms: 'flotation uplift buoyancy seepage cut-off piping heave intake chamber',
  inputs: [
    { key: 'Lc', label: 'Chamber length (outside)', unit: 'm', type: 'number', default: 14, min: 1, max: 300, group: 'Chamber', role: 'geometry' },
    { key: 'Bc', label: 'Chamber width (outside)', unit: 'm', type: 'number', default: 8, min: 1, max: 100, group: 'Chamber', role: 'geometry' },
    { key: 'Hc', label: 'Chamber height (outside)', unit: 'm', type: 'number', default: 7, min: 1, max: 60, group: 'Chamber', role: 'geometry' },
    { key: 'tw', label: 'Wall thickness', unit: 'm', type: 'number', default: 0.6, min: 0.15, max: 5, group: 'Chamber', role: 'geometry', design: true },
    { key: 'tbs', label: 'Base slab thickness', unit: 'm', type: 'number', default: 1.4, min: 0.2, max: 8, group: 'Chamber', role: 'geometry', design: true },
    { key: 'roof', label: 'Roof slab thickness (0 = open top)', unit: 'm', type: 'number', default: 0.4, min: 0, max: 3, group: 'Chamber', role: 'geometry' },
    { key: 'conc', label: 'Concrete', type: 'select', source: 'materials:concrete', default: 'c35', group: 'Chamber', role: 'material' },
    { key: 'extraW', label: 'Permanent equipment and ballast weight', unit: 'kN', type: 'number', default: 300, min: 0, max: 1e6, group: 'Chamber', role: 'ic', help: 'Only weight that is certain to be present in the governing case (emptied for maintenance).' },
    { key: 'resist', label: 'Additional resisting force (anchors, side friction, toe overburden)', unit: 'kN', type: 'number', default: 0, min: 0, max: 1e6, group: 'Chamber', role: 'bc' },
    { key: 'inside', label: 'Water depth inside the chamber', unit: 'm', type: 'number', default: 0, min: 0, max: 60, group: 'Water levels', role: 'ic', help: 'Initial flooded or dry state: 0 for the dewatered maintenance case.' },
    { key: 'hw', label: 'External water level above the underside of the base', unit: 'm', type: 'number', default: 5.2, min: 0, max: 60, group: 'Water levels', role: 'bc', help: 'Design high water or groundwater before the sea-level-rise allowance.' },
    ...slrInputs({ group: 'Water levels' }),
    { key: 'rhoW', label: 'Water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Water levels', role: 'environment', bind: 'site.rhoW' },
    { key: 'dH', label: 'Head difference driving seepage under the base', unit: 'm', type: 'number', default: 3.5, min: 0, max: 30, group: 'Seepage', role: 'bc', help: 'Sea or river side above the protected / dewatered side.' },
    { key: 'T', label: 'Thickness of the permeable layer below the base', unit: 'm', type: 'number', default: 10, min: 1, max: 80, group: 'Seepage', role: 'soil' },
    { key: 'k', label: 'Horizontal permeability', unit: 'm/s', type: 'number', default: 1e-4, min: 1e-9, max: 1e-1, group: 'Seepage', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.8 } },
    { key: 'aniso', label: 'Vertical / horizontal permeability', unit: '-', type: 'number', default: 0.5, min: 0.01, max: 1, group: 'Seepage', role: 'soil' },
    { key: 'cut', label: 'Cut-off wall depth at the downstream edge', unit: 'm', type: 'number', default: 3, min: 0, max: 60, group: 'Seepage', role: 'geometry', design: true, help: '0 for none. Must stop short of the impervious base.' },
    { key: 'gamSub', label: 'Submerged unit weight of the soil', unit: 'kN/m³', type: 'number', default: 9.5, min: 4, max: 13, group: 'Seepage', role: 'soil', bind: 'site.gammaSub' },
    { key: 'sfF', label: 'Required flotation factor', unit: '-', type: 'number', default: 1.1, min: 1, max: 2, group: 'Criteria', role: 'criteria', help: 'USACE: 1.1 extreme to 1.3–1.5 usual conditions; EN 1997 UPL uses 0.9 on stabilising weight.' },
    { key: 'sfP', label: 'Required factor against piping (exit gradient)', unit: '-', type: 'number', default: 3, min: 1.5, max: 8, group: 'Criteria', role: 'criteria' },
    { key: 'sfH', label: 'Required factor against heave', unit: '-', type: 'number', default: 1.5, min: 1.1, max: 4, group: 'Criteria', role: 'criteria' },
    { key: 'nx', label: 'Seepage grid columns', type: 'int', default: 120, min: 30, max: 960, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nx', label: 'Seepage grid columns', levels: [60, 120, 240], qoi: 'upliftSeep' },
  compute(inp, ctx) {
    const c = ctx.material(inp.conc); need(c && c.rho > 0, 'The selected concrete is not in the materials library.');
    need(2 * inp.tw < Math.min(inp.Lc, inp.Bc), 'The walls are thicker than the chamber is wide.');
    need(inp.tbs + inp.roof < inp.Hc, 'Base and roof slabs are thicker than the chamber is high.');
    need(inp.cut < inp.T, 'The cut-off reaches the impervious base: seepage is then blocked and this analysis does not apply. Shorten the cut-off or treat the base as sealed.');
    const rhoC = c.rhoReinforced || c.rho, A = inp.Lc * inp.Bc;
    const Vc = A * inp.Hc - (inp.Lc - 2 * inp.tw) * (inp.Bc - 2 * inp.tw) * (inp.Hc - inp.tbs - inp.roof);
    const Wc = Vc * rhoC * G, Win = inp.rhoW * G * (inp.Lc - 2 * inp.tw) * (inp.Bc - 2 * inp.tw) * Math.min(inp.inside, inp.Hc - inp.tbs - inp.roof);
    const W = Wc + inp.extraW * 1e3 + Win;
    const flo = (rise) => flotationCheck({ weight: W, displacedVolume: Math.max(1e-9, A * Math.min(inp.hw + rise, inp.Hc)), resisting: inp.resist * 1e3, rho: inp.rhoW });
    const fl = flo(inp.slr), hDes = Math.min(inp.hw + inp.slr, inp.Hc);
    const pBase = inp.rhoW * G * hDes, wallThrust = 0.5 * inp.rhoW * G * Math.max(hDes - inp.inside, 0) ** 2;
    const ballast = Math.max(0, fl.ballastForSf(inp.sfF)), slabExtra = ballast / ((rhoC - inp.rhoW) * G * A);
    // water level at which the chamber just floats (FS = 1)
    const hFloat = (W + inp.resist * 1e3) / (inp.rhoW * G * A);
    // seepage under the base: domain five base widths wide, structure centred, cut-off at the downstream edge
    let sp = null;
    if (inp.dH > 0) {
      const width = 5 * inp.Bc, x1 = 2 * inp.Bc, x2 = 3 * inp.Bc, nx = Math.round(inp.nx), nz = Math.max(10, Math.round(nx * inp.T / width));
      sp = run(() => confinedSeepage2D({ width, thickness: inp.T, nx, nz: Math.min(nz, 160), kx: inp.k, kz: inp.k * inp.aniso, headUp: inp.dH, headDown: 0, structure: { x1, x2 }, cutoffs: inp.cut > 0 ? [{ x: x2, depth: inp.cut }] : [], gammaSub: inp.gamSub * 1e3 }), 'Seepage');
    }
    const upliftSeep = sp ? sp.uplift.force : 0, iExit = sp ? sp.exitGradient : 0, iCr = criticalGradient({ gammaSub: inp.gamSub * 1e3 });
    const fsPipe = iExit > 1e-9 ? iCr / iExit : 99, fsHeave = sp && Number.isFinite(sp.fsHeave) ? sp.fsHeave : 99, qSeep = sp ? sp.q : 0;
    finite(fl.sf, 'Flotation factor');
    const scen = slrRows(inp).map((s) => { const q = flo(s.rise); return { name: s.name, rise: s.rise, h: r(Math.min(inp.hw + s.rise, inp.Hc), 3), U: r(q.buoyancy / 1e3), sf: r(q.sf, 3), ballast: r(Math.max(0, q.ballastForSf(inp.sfF)) / 1e3) }; });
    const ok = fl.sf >= inp.sfF, okP = fsPipe >= inp.sfP && fsHeave >= inp.sfH;
    const hs = linspace(0, inp.Hc, 30);

    const advice = [
      ok
        ? eng('ok', `Flotation factor is ${f(fl.sf)} (need ${f(inp.sfF)}) with the water ${f(hDes)} m above the base including ${f(inp.slr)} m of sea-level rise; net downward force ${f(fl.netDownward / 1e3)} kN.`,
          `The chamber would float when the water reaches ${f(hFloat)} m above the base — ${f(hFloat - hDes, 2)} m above the design level. Post that level as the limit for dewatering.`, 'USACE EM 1110-2-2100 (2005) §3-10; EN 1997-1:2004 §10.2 (UPL)',
          'A chamber that is safe by its own weight needs no tension anchors or pumped pressure relief, avoiding steel tendons and lifetime pumping energy.')
        : eng('fail', `The emptied chamber is not safe against flotation: factor ${f(fl.sf)} against ${f(inp.sfF)} required (buoyancy ${f(fl.buoyancy / 1e3)} kN, weight ${f(W / 1e3)} kN).`,
          `Add ${f(ballast / 1e3)} kN of permanent weight — about ${f(slabExtra, 2)} m more base slab — or extend the base as a toe to mobilise overburden, install tension piles or anchors, or prohibit dewatering above ${f(hFloat)} m water level.`, 'USACE EM 1110-2-2100 (2005) §3-10; EN 1997-1:2004 §10.2 (UPL)',
          `Thickening the slab by ${f(slabExtra, 2)} m takes ${f(slabExtra * A, 3)} m³ of concrete (≈ ${f(slabExtra * A * CO2_CONC, 3)} t CO₂e, screening factor); a base toe that engages soil weight or a dewatering limit achieves the same with little or no concrete.`),
      eng(okP ? 'info' : 'fail', `Seepage under the base: ${f(qSeep * 86400, 3)} m³/day per metre, uplift ${f(upliftSeep / 1e3)} kN/m, exit gradient ${f(iExit, 3)} (critical ${f(iCr, 3)}): piping factor ${f(fsPipe)} (need ${f(inp.sfP)}), heave factor ${f(fsHeave)} (need ${f(inp.sfH)}).`,
        okP ? 'Seepage is controlled. Confirm permeability and layering by pumping tests; an inverted filter at the exit adds robustness.' : 'Deepen the cut-off, add an upstream blanket, or place a weighted filter berm at the exit; piping progresses backwards unseen and fails without warning.', 'EN 1997-1:2004 §10.3 (HYD) and §10.4; USACE EM 1110-2-1901 ch. 4'),
      eng(scen[scen.length - 1].sf < inp.sfF ? 'warn' : 'info', `Sea-level rise: the flotation factor falls from ${f(scen[0].sf)} to ${f(scen[scen.length - 1].sf)} at +${f(scen[scen.length - 1].rise)} m (${scen[scen.length - 1].name}).`,
        'Size the permanent weight for the end-of-life groundwater and sea level; where that is uneconomic, provide pressure-relief or a controlled-flooding procedure with a monitored trigger level.', SLR_BASIS),
      eco('info', `The chamber contains ${f(Vc, 3)} m³ of concrete (${f(Wc / 1e3)} kN); flotation safety needs ${f(ballast / 1e3)} kN more and seepage costs ${f(qSeep * 86400 * inp.Lc, 3)} m³/day of pumping when dewatered.`,
        'Compare mass concrete ballast with tension anchors and with an operational dewatering limit on whole-life cost, including the pumping energy for seepage during maintenance.', 'Whole-life costing per ISO 15686-5; AACE RP 18R-97 quantity-based estimating',
        `The chamber concrete embodies about ${f(Vc * CO2_CONC, 3)} t CO₂e (screening factor); a deeper cut-off reduces both uplift and the pumped seepage volume.`),
    ];
    return {
      kpis: [
        { key: 'sfFlot', label: 'Flotation safety factor', value: r(Math.min(fl.sf, 999)), unit: '-', limit: inp.sfF, sense: 'min', explain: 'Weight plus resisting forces over buoyancy for the stated internal water state. Below 1 the chamber lifts; 1.1 is a usual minimum for the emptied case.' },
        { key: 'buoy', label: 'Buoyancy', value: r(fl.buoyancy / 1e3), unit: 'kN', explain: 'Weight of water displaced by the chamber below the external water level (Archimedes).' },
        { key: 'weight', label: 'Total stabilising weight', value: r(W / 1e3), unit: 'kN', explain: 'Concrete, permanent equipment and any water inside. Only weight certain to be present in the governing case should be counted.' },
        { key: 'ballast', label: 'Extra weight for the required factor', value: r(ballast / 1e3), unit: 'kN', explain: 'Additional permanent weight (or anchor force) needed to reach the required flotation factor; zero when already satisfied.' },
        { key: 'hFloat', label: 'Water level at which the chamber floats', value: r(hFloat), unit: 'm', explain: 'External level above the underside of the base at which buoyancy equals the stabilising forces. The operating limit for dewatering.' },
        { key: 'pBase', label: 'Hydrostatic pressure at the base', value: r(pBase / 1e3), unit: 'kPa', explain: 'Water pressure ρ g h on the underside of the base slab; the uplift pressure the slab spans against when the chamber is empty.' },
        { key: 'wallThrust', label: 'Net water thrust on a wall', value: r(wallThrust / 1e3), unit: 'kN/m', explain: 'Resultant ½ ρ g h² of the unbalanced head on each metre of wall, acting at one third of the head above its foot.' },
        { key: 'upliftSeep', label: 'Seepage uplift under the base', value: r(upliftSeep / 1e3), unit: 'kN/m', explain: 'Resultant of the excess pore pressure under the base from seepage, per metre run, additional to the downstream hydrostatic level.' },
        { key: 'qSeep', label: 'Seepage discharge', value: r(qSeep * 86400, 4), unit: 'm³/day/m', explain: 'Flow under the structure per metre run; the pumping demand of a dewatered excavation or chamber.' },
        { key: 'iExit', label: 'Exit hydraulic gradient', value: r(iExit, 4), unit: '-', explain: 'Upward gradient where seepage emerges. Near the critical gradient (about 1) the soil boils and piping starts.' },
        { key: 'fsPipe', label: 'Factor against piping', value: r(Math.min(fsPipe, 99)), unit: '-', limit: inp.sfP, sense: 'min', explain: 'Critical over exit gradient. Practice asks for 3 or more because local variations of permeability concentrate the flow.' },
        { key: 'fsHeave', label: 'Factor against heave', value: r(Math.min(fsHeave, 99)), unit: '-', limit: inp.sfH, sense: 'min', explain: 'Submerged weight of the soil prism beside the cut-off over the mean seepage uplift on it (Terzaghi). 99 when there is no cut-off or no seepage.' },
      ],
      checks: [
        sfCheck('Flotation of the chamber', fl.sf, inp.sfF, 'USACE EM 1110-2-2100 §3-10; EN 1997-1 §10.2', 'Guards against the chamber lifting when emptied at high water — sudden, and it breaks every connected pipe and joint.'),
        sfCheck('Piping at the seepage exit', fsPipe, inp.sfP, 'EN 1997-1 §10.3; USACE EM 1110-2-1901', 'Guards against soil being carried out by emerging seepage; a pipe then works back under the structure and undermines it.'),
        sfCheck('Heave beside the cut-off', fsHeave, inp.sfH, 'Terzaghi (1943); EN 1997-1 §10.3', 'Guards against the soil block on the downstream side of the cut-off being lifted by seepage pressure.'),
      ],
      plots: [
        { type: 'line', title: 'Flotation factor versus water level', x: { label: 'External water level above base', unit: 'm' }, y: { label: 'Flotation factor', unit: '-' },
          series: [{ name: 'Factor', x: hs.slice(1).map((h) => r(h, 3)), y: hs.slice(1).map((h) => r(Math.min((W + inp.resist * 1e3) / (inp.rhoW * G * A * h), 5), 3)) }, { name: 'Required', x: [0, inp.Hc], y: [inp.sfF, inp.sfF], dash: true }, { name: 'Design level', x: [hDes, hDes], y: [0, Math.min(fl.sf, 5)], dash: true }],
          explain: 'How the margin against flotation falls as the water rises. The crossing with the required line is the highest level at which the chamber may be emptied.' },
        ...(sp ? [{ type: 'line', title: 'Uplift pressure under the base from seepage', x: { label: 'Distance across the base', unit: 'm' }, y: { label: 'Excess pressure', unit: 'kPa' },
          series: [{ name: 'Uplift', x: sp.uplift.x.map((x) => r(x - 2 * inp.Bc, 3)), y: sp.uplift.pressure.map((p) => r(p / 1e3, 3)) }],
          explain: 'Seepage pressure along the underside of the base from the high-head side to the cut-off. A downstream cut-off raises the pressure under the base while lowering the exit gradient.' }] : []),
      ],
      tables: [
        { title: 'Sea-level-rise scenarios: flotation', columns: [{ key: 'name', label: 'Scenario' }, { key: 'rise', label: 'Rise', unit: 'm' }, { key: 'h', label: 'Water above base', unit: 'm' }, { key: 'U', label: 'Buoyancy', unit: 'kN' }, { key: 'sf', label: 'Factor', unit: '-' }, { key: 'ballast', label: 'Extra weight needed', unit: 'kN' }],
          rows: scen, explain: 'Buoyancy and flotation factor with the external level raised by each scenario (capped at the top of the chamber), and the weight each would need.' },
        { title: 'Weights and volumes', columns: [{ key: 'q', label: 'Item' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['Concrete volume', Vc, 'm³'], ['Concrete weight', Wc / 1e3, 'kN'], ['Equipment and ballast', inp.extraW, 'kN'], ['Water inside', Win / 1e3, 'kN'], ['Resisting force', inp.resist, 'kN'], ['Displaced volume', A * hDes, 'm³'], ['Critical gradient', iCr, '-']].map(([q, v, u]) => ({ q, v: r(v), u })),
          explain: 'The terms of the flotation balance, for an independent hand check.' },
      ],
      notes: [`Reinforced-concrete density ${f(rhoC)} kg/m³; ${inp.inside > 0 ? `${f(inp.inside)} m of water inside` : 'chamber empty'}.`, sp ? `Seepage grid ${sp.nx} × ${sp.nz} cells${sp.converged ? '' : ' (solver not fully converged)'}; cut-off ${f(inp.cut)} m.` : 'No head difference: seepage not analysed.'],
      advice,
      publish: { 'marine.upliftForce': fl.buoyancy, 'marine.flotationSF': Math.min(fl.sf, 999) },
      quantities: { concrete_m3: r(Vc + (ok ? 0 : slabExtra * A)) },
    };
  },
  verify: [
    { name: 'Archimedes: flotation factor = weight / (ρ g V)', source: 'Hand calculation: open chamber 10 × 6 × 5 m, walls 0.5 m, base 0.8 m → concrete 300 − 9 × 5 × 4.2 = 111 m³ at 2500 kg/m³ = 2721.3 kN; water 4 m above base → buoyancy 1025 × 9.80665 × 240 = 2412.4 kN; factor 1.128',
      inputs: { Lc: 10, Bc: 6, Hc: 5, tw: 0.5, tbs: 0.8, roof: 0, extraW: 0, resist: 0, inside: 0, hw: 4, slr: 0, rhoW: 1025, conc: 'c35' }, expect: [{ kpi: 'sfFlot', value: 1.128, tol: 0.002 }, { kpi: 'buoy', value: 2412.4, tol: 0.001 }] },
    { name: 'Symmetric seepage without a cut-off: mean uplift head under the base = ΔH/2', source: 'Antisymmetry of the Laplace solution about the centre of the base (Harr 1962, Groundwater and Seepage §3). ΔH = 4 m, B = 10 m, γw = 9.81 kN/m³ → uplift = 9.81 × 2 × 10 = 196.2 kN/m',
      inputs: { Bc: 10, dH: 4, cut: 0, T: 10, nx: 240 }, expect: [{ kpi: 'upliftSeep', value: 196.2, tol: 0.01 }] },
  ],
};
