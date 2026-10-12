// Seismic and liquefaction assessment of a monopile foundation: SPT-based liquefaction triggering with depth,
// potential index and post-liquefaction settlement, modal response-spectrum loads of the turbine–tower–pile system,
// and the pile response with degraded soil reactions in the liquefied layers.
import { sptProfile, liquefactionPotentialIndex, postLiquefactionSettlement } from '../../engine/geotech/liquefaction.js';
import { ec8Spectrum, combineModal } from '../../engine/structural/framedynamics.js';
import { solveLateralPile } from '../../engine/geotech/pile.js';
import { monopileStiffness } from '../../engine/geotech/monopile.js';
import { G, DEG, need, run, fmt, round, cap, adv, check, worst, sev, steelCo2Text, turbineInputs, supportInputs, steelInputs, steelOf, soilInputs, soilProfile, turbine, LIMIT_PY, LIMIT_CODE, CONFIRM } from './_common.js';
import { pileOptions } from './_pile.js';
import { supportStructure } from './_dynamics.js';

export default {
  id: 's03-seismic-liquefaction', suite: 3, title: 'Seismic and liquefaction assessment',
  summary: 'Factor of safety against liquefaction with depth, liquefaction potential index and settlement, earthquake base shear and overturning moment from a modal response-spectrum analysis, and the pile response with liquefied layers.',
  covers: [11], fidelity: 'reduced-order',
  standards: ['IEC 61400-1:2019 §11.6 and Annex D (earthquake)', 'DNV-RP-0585 (2021) seismic design of wind power plants', 'EN 1998-1:2004 §3.2.2 (elastic response spectrum), EN 1998-5:2004 §4.1.4 and Annex B (liquefaction)', 'ISO 19901-2:2017 (seismic design procedures)', 'Youd et al. (2001) NCEER liquefaction procedure', 'Iwasaki et al. (1982); Idriss & Boulanger (2008)'],
  models: ['Simplified cyclic stress ratio (Seed & Idriss) with stress reduction and magnitude scaling', 'SPT clean-sand cyclic resistance curve (NCEER)', 'Liquefaction potential index and post-liquefaction volumetric strain', 'Modal response-spectrum analysis with SRSS/CQC combination', 'EN 1998-1 elastic response spectrum', 'p-multiplier degradation of soil reactions in liquefied layers', 'Three-spring foundation stiffness and beam eigenfrequencies'],
  limits: `Horizontal ground motion in one direction on a planar stick model; vertical excitation, soil–structure kinematic interaction, lateral spreading loads and site-specific ground response are not modelled (the spectrum's ground type stands in for site amplification). Liquefaction triggering is total-stress and SPT-based. ${LIMIT_PY} ${LIMIT_CODE} Effective-stress site-response and 3-D soil–structure analysis are outside the module.`,
  searchTerms: 'offshore wind turbine monopile earthquake liquefaction response spectrum',
  inputs: [
    ...turbineInputs(),
    ...supportInputs(),
    ...steelInputs(),
    ...soilInputs(),
    { key: 'pga', label: 'Peak ground acceleration on rock (475-year)', unit: 'g', type: 'number', default: 0.12, min: 0.01, max: 1, group: 'Earthquake', role: 'bc', bind: 'site.pga475', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Reference acceleration on type A ground for the design return period.' },
    { key: 'Mw', label: 'Moment magnitude', unit: '-', type: 'number', default: 6.5, min: 5, max: 9, group: 'Earthquake', role: 'bc', help: 'Controls the number of loading cycles through the magnitude scaling factor.' },
    { key: 'ground', label: 'Ground type (EN 1998-1)', type: 'select', default: 'C', group: 'Earthquake', role: 'bc', options: ['A', 'B', 'C', 'D', 'E'].map((v) => ({ value: v, label: `Type ${v}` })), help: 'A rock, B very dense sand / stiff clay, C dense to medium sand, D loose sand / soft clay, E shallow soft layer on rock.' },
    { key: 'specType', label: 'Spectrum type', type: 'select', default: 1, group: 'Earthquake', role: 'bc', options: [{ value: 1, label: 'Type 1 (Ms > 5.5)' }, { value: 2, label: 'Type 2 (Ms ≤ 5.5)' }] },
    { key: 'zetaE', label: 'Damping in the earthquake', unit: '%', type: 'number', default: 2, min: 0.5, max: 10, group: 'Earthquake', role: 'ic', help: 'Structural plus soil damping; 1–2 % for a parked turbine, about 5 % if aerodynamic damping in production is counted.' },
    { key: 'T', label: 'Rotor thrust acting with the earthquake', unit: 'kN', type: 'number', default: 715, min: 0, max: 20000, group: 'Earthquake', role: 'bc', bind: 'owt.thrustRated', help: 'IEC 61400-1 combines the earthquake with normal operation (or with an emergency stop triggered by it).' },
    { key: 'spt', label: 'SPT profile', type: 'table', group: 'Liquefaction', role: 'soil',
      columns: [{ key: 'z', label: 'Depth', unit: 'm', type: 'number' }, { key: 'N160', label: '(N1)60', unit: '-', type: 'number' }, { key: 'FC', label: 'Fines content', unit: '%', type: 'number' }, { key: 'gamma', label: 'Saturated unit weight', unit: 'kN/m³', type: 'number' }],
      default: [{ z: 1.5, N160: 16, FC: 8, gamma: 18.5 }, { z: 3, N160: 18, FC: 8, gamma: 18.5 }, { z: 4.5, N160: 20, FC: 6, gamma: 19 }, { z: 6, N160: 22, FC: 5, gamma: 19 }, { z: 8, N160: 25, FC: 5, gamma: 19.5 }, { z: 10, N160: 27, FC: 5, gamma: 19.5 }, { z: 13, N160: 32, FC: 5, gamma: 20 }, { z: 16, N160: 36, FC: 5, gamma: 20 }, { z: 20, N160: 40, FC: 5, gamma: 20 }],
      help: 'Overburden- and energy-corrected blow counts from the boreholes (or converted from CPT), with fines content and total unit weight of the soil above each depth.' },
    { key: 'fsLim', label: 'Required factor of safety against liquefaction', unit: '-', type: 'number', default: 1.25, min: 1, max: 2, group: 'Liquefaction', role: 'criteria', help: `EN 1998-5 §4.1.4 uses 1.25 on the cyclic resistance. ${CONFIRM}` },
    { key: 'pLiq', label: 'p-multiplier in liquefied soil', unit: '-', type: 'number', default: 0.1, min: 0, max: 0.5, group: 'Liquefaction', role: 'soil', help: 'Residual lateral resistance of liquefied sand as a fraction of the non-liquefied p–y curve (0.05–0.3 depending on density).' },
    { key: 'gammaF', label: 'Load factor on seismic action', unit: '-', type: 'number', default: 1.0, min: 1, max: 1.5, group: 'Acceptance criteria', role: 'criteria', help: 'IEC 61400-1 §11.6: 1.0 on the earthquake combined with operational loads.' },
    { key: 'gammaM', label: 'Material factor on yield', unit: '-', type: 'number', default: 1.1, min: 1, max: 1.5, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'nElem', label: 'Pile elements', type: 'int', default: 40, min: 10, max: 160, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nElem', label: 'Pile elements', levels: [20, 40, 80], qoi: 'yIntact' },
  compute(inp, ctx) {
    const steel = steelOf(inp, ctx);
    const tb = turbine(inp);
    const amax = inp.pga * G;
    const spec = run(() => ec8Spectrum({ ag: amax, groundType: inp.ground, type: Number(inp.specType), zeta: inp.zetaE / 100 }), 'Response spectrum');
    // Liquefaction triggering on the site-amplified surface acceleration ag·S.
    const rows = (inp.spt || []).map((r) => ({ z: Number(r.z), N160: Number(r.N160), FC: Number(r.FC) || 0, gamma: Number(r.gamma) * 1e3 })).filter((r) => r.z > 0 && r.N160 > 0 && r.gamma > 0).sort((a, b) => a.z - b.z);
    need(rows.length >= 1, 'The SPT table is empty: enter at least one depth with a blow count and unit weight.');
    for (let i = 1; i < rows.length; i++) need(rows[i].z > rows[i - 1].z, 'The SPT table has repeated depths.');
    const aSurf = amax * spec.S;
    const prof = run(() => sptProfile({ points: rows, waterTableDepth: 0, amax: aSurf, Mw: inp.Mw }), 'Liquefaction profile');
    const fs = prof.map((p) => cap(p.FS, 5));
    const fsMin = Math.min(...fs), iMin = fs.indexOf(fsMin);
    const lpi = liquefactionPotentialIndex(prof.map((p) => p.z), fs);
    const thick = prof.map((p, i) => (i === 0 ? p.z : p.z - prof[i - 1].z));
    const liq = prof.map((p, i) => fs[i] < 1);
    const zLiq = prof.reduce((a, p, i) => (liq[i] ? p.z : a), 0);
    const tLiq = thick.reduce((a, t, i) => a + (liq[i] ? t : 0), 0);
    const settle = postLiquefactionSettlement(prof.map((p, i) => ({ thickness: thick[i], N160cs: p.N160cs, FS: Math.max(fs[i], 0.05) }))).settlement;
    // Foundation stiffness: intact and with liquefied layers (p-multiplier applied down to the deepest liquefied depth).
    const profile = soilProfile(inp);
    const optBase = pileOptions({ ...inp, pyMode: 'api', loading: 'cyclic' }, profile, steel);
    const isLiq = (z) => { for (let i = 0; i < prof.length; i++) if (z <= prof[i].z) return liq[i]; return false; };
    const optLiq = { ...optBase, pMult: (z) => (isLiq(z) ? Math.max(inp.pLiq, 1e-3) : 1) };
    const strOf = (K) => supportStructure({ ...inp, fnd: 'springs', Kh: K.KLL / 1e9, Kr: K.KRR / 1e9, Khr: K.KLR / 1e6, zetaS: inp.zetaE }, tb, steel, { nElem: 30, nModes: 3 });
    const lever = inp.hubHeight + inp.depth;
    const Tn = inp.T * 1e3;
    const seismic = (opt) => {
      // secant stiffness at a representative seismic load level (iterated once on the response-spectrum result)
      let head = { H: 1.5e6, M: 1.0e8 };
      let K = run(() => monopileStiffness({ ...opt, head }), 'Foundation stiffness').secant, str = strOf(K), res = null;
      for (let it = 0; it < 2; it++) {
        const om = str.modes.map((m) => 2 * Math.PI * m.frequency);
        const q = str.modes.map((m, i) => (m.participation * spec.Sa(1 / m.frequency)) / (om[i] * om[i]));
        const Mi = str.modes.map((m, i) => m.momentMudline * q[i]), Vi = str.modes.map((m, i) => m.shearMudline * q[i]), Xi = str.modes.map((m, i) => m.phiHub * q[i]);
        const cq = (p) => combineModal({ peaks: p, omegas: om, zeta: inp.zetaE / 100, method: 'CQC' });
        res = { M: cq(Mi), V: cq(Vi), x: cq(Xi), acc: cq(Xi.map((x, i) => x * om[i] * om[i])), f: str.modes.map((m) => m.frequency), meff: str.modes.map((m) => m.effectiveMass), str, K };
        head = { H: Math.max(res.V + Tn, 1e5), M: Math.max(res.M + Tn * lever, 1e6) };
        if (it === 0) { K = monopileStiffness({ ...opt, head }).secant; str = strOf(K); }
      }
      return res;
    };
    const s0 = seismic(optBase);
    const anyLiq = tLiq > 0;
    const sL = anyLiq ? seismic(optLiq) : s0;
    const Md = inp.gammaF * sL.M + Tn * lever, Hd = inp.gammaF * sL.V + Tn;
    const r0 = run(() => solveLateralPile({ ...optBase, head: { H: inp.gammaF * s0.V + Tn, M: inp.gammaF * s0.M + Tn * lever } }), 'Pile response, intact soil');
    const rL = run(() => solveLateralPile({ ...optLiq, head: { H: Hd, M: Md }, throwOnFail: false }), 'Pile response, liquefied soil');
    need(rL.convergence.converged, 'The pile is unstable with the liquefied layers at the seismic load: the remaining embedment below the liquefied zone is insufficient. Lengthen the pile or improve the ground.');
    const sec = optBase.sec, fyd = steel.fy / inp.gammaM;
    const sig = rL.maxMoment / sec.W;
    const checks = [
      check(`Liquefaction: required over available factor of safety (z = ${fmt(prof[iMin].z)} m)`, inp.fsLim, Math.max(fsMin, 1e-9), '-', 'EN 1998-5:2004 §4.1.4 (11); Youd et al. (2001)', 'Guards against loss of strength of saturated sand during the earthquake. A liquefied layer gives almost no lateral support to the pile and settles afterwards.'),
      check('Pile bending stress with liquefied layers (seismic + operating)', sig / 1e6, fyd / 1e6, 'MPa', 'IEC 61400-1:2019 §11.6; DNV-RP-0585 (2021) §5', 'Bending stress in the pile when the upper soil has liquefied and the fixity point has moved down. Yielding here means permanent tilt after the event.'),
      check('Mudline rotation in the earthquake with liquefied layers', Math.abs(rL.head.theta) / DEG, 0.5, 'deg', 'DNV-ST-0126 (2021) §7.6.2 (rotation tolerance, applied as a post-event serviceability criterion)', 'Peak rotation of the pile at the mudline during shaking. Large values indicate residual tilt and loss of operability after the earthquake.'),
    ];
    checks[0].uc = Math.min(inp.fsLim / Math.max(fsMin, 1e-9), 99);
    const gov = worst(checks), ok = gov.uc <= 1;
    const Ts = Array.from({ length: 60 }, (_, i) => 0.05 * 1.085 ** i);
    const WpilePerM = sec.A * steel.rho / 1e3;
    const extra = anyLiq ? zLiq : 0;
    return {
      kpis: [
        { key: 'fsMin', label: 'Minimum factor of safety against liquefaction', unit: '-', value: fsMin, limit: inp.fsLim, sense: 'min', explain: `Cyclic resistance over cyclic stress at the most critical depth (${fmt(prof[iMin].z)} m), capped at 5. Below 1 the layer is expected to liquefy in the design earthquake.` },
        { key: 'csrRef', label: 'Cyclic stress ratio at the first test depth', unit: '-', value: prof[0].CSR, explain: 'Earthquake-induced shear stress relative to effective overburden at the shallowest SPT depth: 0.65·(amax/g)·(σv/σ′v)·rd.' },
        { key: 'lpi', label: 'Liquefaction potential index', unit: '-', value: lpi, explain: 'Depth-weighted severity over the top 20 m (Iwasaki): below 5 low, 5–15 high, above 15 severe surface effects.' },
        { key: 'tLiq', label: 'Thickness of liquefied soil', unit: 'm', value: tLiq, explain: 'Total thickness of layers with a factor of safety below 1; these layers lose their lateral support to the pile.' },
        { key: 'settle', label: 'Post-liquefaction seabed settlement', unit: 'mm', value: settle * 1e3, explain: 'Free-field reconsolidation settlement of the liquefied layers. It loads cables and J-tubes at the seabed and may expose scour protection.' },
        { key: 'saPlateau', label: 'Spectral acceleration at the plateau', unit: 'm/s²', value: spec.Sa(0.5 * (spec.TB + spec.TC)), explain: 'Maximum of the elastic response spectrum, 2.5·ag·S·η. The first tower mode sits far to the right of the plateau; higher modes can reach it.' },
        { key: 'f1', label: 'First natural frequency, intact soil', unit: 'Hz', value: s0.f[0], explain: 'Fundamental frequency with the cyclic foundation stiffness at the seismic load level.' },
        { key: 'f1Liq', label: 'First natural frequency with liquefied layers', unit: 'Hz', value: sL.f[0], explain: 'Fundamental frequency once the liquefied layers have lost their stiffness. A large drop may bring the structure into the rotor-speed range after the event.' },
        { key: 'Vseis', label: 'Seismic shear at the mudline', unit: 'MN', value: sL.V / 1e6, explain: 'Combined modal (CQC) base shear from the horizontal earthquake. Higher modes contribute most of it because the first mode is long-period.' },
        { key: 'Mseis', label: 'Seismic moment at the mudline', unit: 'MN·m', value: sL.M / 1e6, explain: 'Combined modal overturning moment at the mudline from the earthquake alone, before adding the operating thrust.' },
        { key: 'Mcomb', label: 'Combined design moment (earthquake + thrust)', unit: 'MN·m', value: Md / 1e6, explain: 'Seismic moment times its load factor plus the operating rotor thrust moment: compare with the storm design moment to see whether the earthquake governs.' },
        { key: 'accHub', label: 'Peak nacelle acceleration', unit: 'm/s²', value: sL.acc, explain: 'Combined modal acceleration at hub height. Nacelle equipment is typically qualified to a few m/s²; higher values need supplier confirmation.' },
        { key: 'yIntact', label: 'Mudline displacement, intact soil', unit: 'mm', value: r0.head.y * 1e3, explain: 'Pile-head displacement under the combined seismic and operating load without liquefaction.' },
        { key: 'yLiq', label: 'Mudline displacement with liquefied layers', unit: 'mm', value: rL.head.y * 1e3, explain: 'Pile-head displacement once the liquefied layers carry only the residual soil reaction.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Factor of safety against liquefaction with depth', x: { label: 'Factor of safety', unit: '-' }, y: { label: 'Depth below seabed', unit: 'm' },
          series: [{ name: 'FS (capped at 5)', x: fs, y: prof.map((p) => -p.z) }, { name: 'FS = 1', x: [1, 1], y: [0, -prof[prof.length - 1].z], dash: true }, { name: `Required ${inp.fsLim}`, x: [inp.fsLim, inp.fsLim], y: [0, -prof[prof.length - 1].z], dash: true }],
          explain: 'Liquefaction safety at each test depth. Points left of the FS = 1 line are layers that liquefy; the deepest of them sets how far down the pile loses support.' },
        { type: 'line', title: 'Cyclic stress and resistance ratios', x: { label: 'Ratio', unit: '-' }, y: { label: 'Depth below seabed', unit: 'm' },
          series: [{ name: 'Cyclic stress ratio (demand)', x: prof.map((p) => p.CSR), y: prof.map((p) => -p.z) }, { name: 'Cyclic resistance ratio, M 7.5', x: prof.map((p) => Math.min(p.CRR75, 0.6)), y: prof.map((p) => -p.z) }],
          explain: 'Earthquake demand against soil resistance (resistance capped at 0.6 for display, where sand is too dense to liquefy). Demand is nearly constant with depth; resistance follows the blow count.' },
        { type: 'line', title: 'Elastic response spectrum and structural periods', x: { label: 'Period', unit: 's' }, y: { label: 'Spectral acceleration', unit: 'm/s²' },
          series: [{ name: `EN 1998-1 type ${inp.specType}, ground ${inp.ground}, ζ = ${inp.zetaE} %`, x: Ts, y: Ts.map((t) => spec.Sa(t)) }, { name: 'Modes, intact soil', x: s0.f.map((f) => 1 / f), y: s0.f.map((f) => spec.Sa(1 / f)) }, { name: 'Modes, liquefied layers', x: sL.f.map((f) => 1 / f), y: sL.f.map((f) => spec.Sa(1 / f)) }],
          explain: 'Design spectrum with the periods of the first three modes marked. The long first period attracts little acceleration; the second and third modes near the plateau drive shear and tower-top acceleration.' },
        { type: 'line', title: 'Pile deflection in the earthquake', x: { label: 'Deflection', unit: 'mm' }, y: { label: 'Depth below mudline', unit: 'm' },
          series: [{ name: 'Intact soil', x: Array.from(r0.y, (v) => v * 1e3), y: Array.from(r0.z, (v) => -v) }, { name: 'Liquefied layers', x: Array.from(rL.y, (v) => v * 1e3), y: Array.from(rL.z, (v) => -v), dash: true }],
          explain: 'Deflected pile under the combined seismic and operating load. With liquefied layers the pile bends over a longer unsupported length and the fixity point moves below the liquefied zone.' },
      ],
      tables: [
        { title: 'Liquefaction assessment by depth', columns: [{ key: 'z', label: 'Depth', unit: 'm' }, { key: 'sv', label: "σ′v", unit: 'kPa' }, { key: 'rd', label: 'rd' }, { key: 'csr', label: 'CSR' }, { key: 'n', label: '(N1)60cs' }, { key: 'crr', label: 'CRR7.5' }, { key: 'fs', label: 'FS' }, { key: 's', label: 'State' }],
          rows: prof.map((p, i) => ({ z: p.z, sv: round(p.sigmaVeff / 1e3), rd: round(p.rd), csr: round(p.CSR), n: round(p.N160cs), crr: Number.isFinite(p.CRR75) ? round(p.CRR75) : 'non-liquefiable', fs: round(fs[i]), s: liq[i] ? 'liquefies' : fs[i] < inp.fsLim ? 'marginal' : 'safe' })),
          explain: 'One row per SPT depth with the stress, demand, resistance and factor of safety, for comparison with the geotechnical report.' },
        { title: 'Modal response-spectrum results', columns: [{ key: 'm', label: 'Mode' }, { key: 'f', label: 'Frequency, intact', unit: 'Hz' }, { key: 'fl', label: 'Frequency, liquefied', unit: 'Hz' }, { key: 'me', label: 'Effective mass', unit: 't' }, { key: 'sa', label: 'Sa', unit: 'm/s²' }],
          rows: s0.f.map((f, i) => ({ m: i + 1, f: round(f), fl: round(sL.f[i]), me: round(sL.meff[i] / 1e3), sa: round(spec.Sa(1 / sL.f[i])) })),
          explain: 'Frequencies, effective masses and spectral accelerations of the modes used. The effective masses should add up to most of the total mass; if not, more modes matter for base shear.' },
      ],
      notes: [anyLiq ? `Layers down to ${fmt(zLiq)} m liquefy (${fmt(tLiq)} m in total): lateral soil reactions there are reduced to ${fmt(inp.pLiq * 100)} %.` : 'No layer liquefies in the design earthquake.',
        `Surface acceleration ${fmt(aSurf / G)} g (soil factor ${spec.S}); the three modes carry ${fmt(sL.meff.reduce((a, b) => a + b, 0) / sL.str.totalMass * 100)} % of the mass.`, `${gov.name} governs (utilisation ${fmt(gov.uc)}).`,
        `Combined design moment ${fmt(Md / 1e6)} MN·m: compare with the storm moment from the extreme-loading module.`],
      advice: [
        ok ? adv('engineering', sev(gov.uc), `Minimum liquefaction factor of safety is ${fmt(fsMin)} (required ${inp.fsLim}); the earthquake adds ${fmt(sL.M / 1e6)} MN·m at the mudline and the pile stress is ${fmt(sig / fyd * 100)} % of design yield.`,
          `Check whether the combined moment ${fmt(Md / 1e6)} MN·m exceeds the storm design moment; confirm the SPT/CPT correlation and magnitude with a site-specific seismic hazard study, and add the emergency-stop case triggered by the earthquake.`,
          'IEC 61400-1:2019 §11.6 and Annex D; DNV-RP-0585 (2021); EN 1998-5:2004 §4.1.4',
          `No ground improvement or extra embedment is needed: each metre of additional pile that is avoided saves ${fmt(WpilePerM)} t of steel ≈ ${fmt(WpilePerM * 1.9)} t CO₂e (screening factor 1.9 t CO₂e/t).`)
          : adv('engineering', 'fail', `${gov.name}: ${fmt(gov.demand)} against ${fmt(gov.capacity)} ${gov.unit}. ${anyLiq ? `Soil liquefies down to ${fmt(zLiq)} m.` : ''}`,
            gov.name.startsWith('Liquefaction') ? `Design the pile to work without lateral support above ${fmt(zLiq || prof[iMin].z)} m (the checks below show stress and rotation in that state), extend the embedment by about ${fmt(extra)} m, or densify the upper sand (vibro-compaction, stone columns).` : `Increase the wall thickness or diameter over the liquefied depth and extend the pile about ${fmt(extra)} m into non-liquefiable soil.`,
            'EN 1998-5:2004 §4.1.4 and §5.4.2; DNV-RP-0585 (2021) §6; IEC 61400-1:2019 §11.6',
            `Extending the pile by ${fmt(extra)} m adds ${fmt(extra * WpilePerM)} t — ${steelCo2Text(extra * WpilePerM)} Ground improvement avoids the steel but needs a separate marine campaign and disturbs ${fmt(Math.PI * (2 * inp.D) ** 2)} m² of seabed per position.`),
        adv('economic', anyLiq ? 'warn' : 'info', anyLiq ? `Liquefaction to ${fmt(zLiq)} m would need about ${fmt(extra * WpilePerM)} t of extra pile steel per position, and ${fmt(settle * 1e3)} mm of seabed settlement affects cable entries.` : `The seismic case adds ${fmt(sL.M / 1e6)} MN·m at the mudline, ${fmt(sL.M / Math.max(Tn * lever, 1) * 100)} % of the operating thrust moment.`,
          'Spend on site-specific probabilistic seismic hazard and CPT-based liquefaction assessment before adding steel: generic code spectra and SPT correlations are conservative and the investigation cost is small against steel across the array.',
          'Value-of-information for geotechnical and seismic studies (ISO 19901-2 §6, ISO 19901-8); whole-life costing per ISO 15686-5'),
      ],
    };
  },
  verify: [
    { name: 'Cyclic stress ratio: CSR = 0.65·(amax/g)·(σv/σ′v)·rd', source: 'Seed & Idriss (1971) with rd = 1 − 0.00765 z (Liao & Whitman 1986): z = 5 m, γ = 19 kN/m³ → σv = 95 kPa, σ′v = 95 − 49.05 = 45.95 kPa, rd = 0.96175; amax = 0.2 g (ground type A, S = 1) → CSR = 0.65 × 0.2 × 2.0675 × 0.96175 = 0.2585',
      inputs: { pga: 0.2, ground: 'A', spt: [{ z: 5, N160: 15, FC: 5, gamma: 19 }, { z: 10, N160: 30, FC: 5, gamma: 19 }] }, expect: [{ kpi: 'csrRef', value: 0.2585, tol: 0.001 }] },
    { name: 'Plateau of the EN 1998-1 elastic spectrum: Se = 2.5·ag·S·η', source: 'EN 1998-1:2004 Eq. (3.3), type 1, ground type C (S = 1.15), ζ = 5 % (η = 1): 2.5 × 0.2 × 9.80665 × 1.15 = 5.639 m/s²', inputs: { pga: 0.2, ground: 'C', specType: 1, zetaE: 5 }, expect: [{ kpi: 'saPlateau', value: 5.639, tol: 0.001 }] },
  ],
};
