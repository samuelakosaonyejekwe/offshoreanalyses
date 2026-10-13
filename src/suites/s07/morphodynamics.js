// Morphodynamic bed evolution near a marine structure: Exner sediment continuity driven by a current over a
// trench or mound (outfall trench infill, intake channel siltation, disposal mound migration), and the
// equilibrium-profile response of the adjacent shoreline to sea-level rise (Dean profile, Bruun rule).
import { exnerEvolve, deanParameter, deanProfile, bruunRule, depthOfClosure } from '../../engine/coastal/morpho.js';
import { engelundHansen, bedShearCurrent, grainRoughness, settlingVelocitySoulsby, criticalShearStress } from '../../engine/coastal/sediment.js';
import { CO2_MGO, MGO_T_PER_DAY, need, run, finite, r, f, linspace, thin, check, eng, eco, slrInputs, slrRows, SLR_BASIS } from './_common.js';

export default {
  id: 's07-morphodynamics', suite: 7, title: 'Seabed evolution and shoreline response',
  summary: 'How fast a dredged trench or channel silts up (or a mound migrates) under the current, and how far the shoreline retreats as the sea rises.',
  covers: [6], fidelity: 'numerical',
  standards: ['Soulsby, Dynamics of Marine Sands (1997) ch. 10–11', 'USACE Coastal Engineering Manual Part III ch. 3 and Part V ch. 4', 'CIRIA C683 Rock Manual (2007) §4.6', 'PIANC Report 102 (2008) minimising harbour siltation'],
  models: ['Exner sediment continuity', 'Engelund–Hansen total-load transport', 'Bed-slope transport correction', 'Quasi-steady flow continuity (q = U h)', 'Flux-limited Lax–Wendroff finite volumes', 'Dean equilibrium beach profile', 'Bruun rule', 'Hallermeier depth of closure'],
  limits: 'One horizontal dimension with quasi-steady flow that adapts instantly to the bed (no flow separation, no lag of suspended sediment, no waves stirring the bed, no tidal reversal); the infill of deep, narrow trenches is therefore over-predicted in rate but conservative in volume. Shoreline retreat is the equilibrium Bruun estimate, not a process-based forecast; coupled wave–sediment–morphology (2-D/3-D) models are outside the module.',
  publishes: ['marine.shorelineRetreat', 'marine.infillRate'],
  searchTerms: 'trench infill siltation Exner morphodynamics Bruun rule shoreline retreat',
  inputs: [
    { key: 'feature', label: 'Initial bed feature', type: 'select', default: 'trench', options: [{ value: 'trench', label: 'Dredged trench / channel (below bed)' }, { value: 'mound', label: 'Mound / bar / disposal site (above bed)' }], group: 'Initial seabed', role: 'ic' },
    { key: 'amp', label: 'Feature depth or height', unit: 'm', type: 'number', default: 2, min: 0.01, max: 15, group: 'Initial seabed', role: 'ic', design: true },
    { key: 'width', label: 'Feature width (Gaussian, ±2σ)', unit: 'm', type: 'number', default: 80, min: 5, max: 2000, group: 'Initial seabed', role: 'ic', design: true },
    { key: 'h0', label: 'Undisturbed water depth', unit: 'm', type: 'number', default: 8, min: 1, max: 60, group: 'Flow boundary', role: 'bc', bind: 'site.depth' },
    { key: 'U0', label: 'Representative current over the undisturbed bed', unit: 'm/s', type: 'number', default: 0.9, min: 0.1, max: 4, group: 'Flow boundary', role: 'bc', bind: 'site.Uc1', uncertain: { dist: 'lognormal', cov: 0.15 }, help: 'Morphologically representative (transport-weighted) tidal or river current, not the extreme.' },
    { key: 'supply', label: 'Upstream sediment supply / transport capacity', unit: '-', type: 'number', default: 1, min: 0, max: 3, group: 'Flow boundary', role: 'bc', help: '1 = equilibrium inflow; < 1 starved (downdrift of a breakwater or dam); > 1 over-supplied.' },
    { key: 'days', label: 'Simulated duration', unit: 'days', type: 'number', default: 60, min: 0.1, max: 3650, group: 'Flow boundary', role: 'bc' },
    { key: 'd50', label: 'Median grain size d50', unit: 'mm', type: 'number', default: 0.25, min: 0.063, max: 5, group: 'Sediment', role: 'soil', bind: 'site.d50', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'rhoS', label: 'Grain density', unit: 'kg/m³', type: 'number', default: 2650, min: 1800, max: 3200, group: 'Sediment', role: 'soil' },
    { key: 'por', label: 'Bed porosity', unit: '-', type: 'number', default: 0.4, min: 0.25, max: 0.6, group: 'Sediment', role: 'soil' },
    { key: 'beta', label: 'Bed-slope transport coefficient β', unit: '-', type: 'number', default: 1.6, min: 0, max: 5, group: 'Sediment', role: 'criteria', help: 'qs is multiplied by (1 − β ∂zb/∂x): Soulsby–Van Rijn use 1.6. It smooths the bed and is the calibration parameter for trench infill against surveys; 0 gives pure migration.' },
    { key: 'allowFill', label: 'Tolerable infill before maintenance dredging', unit: '%', type: 'number', default: 50, min: 5, max: 100, group: 'Criteria', role: 'criteria' },
    { key: 'He', label: 'Wave height exceeded 12 h per year', unit: 'm', type: 'number', default: 3.5, min: 0.3, max: 12, group: 'Shoreline', role: 'bc', help: 'Effective wave height for the depth of closure (Hallermeier).' },
    { key: 'Te', label: 'Associated wave period', unit: 's', type: 'number', default: 9, min: 3, max: 20, group: 'Shoreline', role: 'bc' },
    { key: 'berm', label: 'Beach berm height above mean sea level', unit: 'm', type: 'number', default: 2, min: 0.2, max: 10, group: 'Shoreline', role: 'geometry' },
    { key: 'beachW', label: 'Beach width available in front of the structure', unit: 'm', type: 'number', default: 60, min: 0, max: 1000, group: 'Shoreline', role: 'criteria', help: 'Retreat beyond this exposes the seawall toe.' },
    ...slrInputs({ group: 'Shoreline', slr: 0.5 }),
    { key: 'nx', label: 'Grid cells', type: 'int', default: 400, min: 40, max: 3200, group: 'Solver', role: 'control' },
    { key: 'dredgeRate', label: 'Dredger production', unit: 'm³/day', type: 'number', default: 8000, min: 200, max: 200000, group: 'Maintenance', role: 'operation' },
    { key: 'reachW', label: 'Width of trench / channel across the flow', unit: 'm', type: 'number', default: 300, min: 1, max: 20000, group: 'Maintenance', role: 'geometry' },
  ],
  mesh: { key: 'nx', label: 'Grid cells', levels: [50, 100, 200, 400], qoi: 'infill' },
  compute(inp) {
    const d = inp.d50 / 1000, n = Math.round(inp.nx), sgn = inp.feature === 'trench' ? -1 : 1;
    need(inp.feature === 'trench' || inp.amp < 0.7 * inp.h0, 'A mound higher than about 70 % of the water depth chokes the flow: the quasi-steady transport model no longer applies.');
    const sig = inp.width / 4, Ldom = Math.max(12 * inp.width, 40 * inp.h0), dx = Ldom / n, x0 = 0.3 * Ldom;
    const xs = Array.from({ length: n }, (_, i) => (i + 0.5) * dx);
    const zb0 = Float64Array.from(xs, (x) => -inp.h0 + sgn * inp.amp * Math.exp(-0.5 * ((x - x0) / sig) ** 2));
    const q = inp.U0 * inp.h0;
    // constant Chézy roughness from the grain log law at the undisturbed depth; τ = ρ g U²/C²
    const C = bedShearCurrent({ U: inp.U0, h: inp.h0, z0: grainRoughness(d).z0 }).chezy;
    const qsOf = (h) => { const U = q / h; return engelundHansen({ U, tau: 1025 * 9.80665 * U * U / (C * C), d50: d, rhoS: inp.rhoS }).qt; };
    const qs0 = finite(qsOf(inp.h0), 'Transport rate');
    need(qs0 > 0, 'The transport rate is zero for this current and grain size: nothing will move.');
    const tEnd = inp.days * 86400;
    const qsMax = Math.max(qs0, qsOf(Math.max(inp.h0 - (sgn > 0 ? inp.amp : 0), 0.05 * inp.h0)));
    // work estimate before marching: advective and diffusive step limits
    const dtAdv = 0.8 * dx * inp.h0 * (1 - inp.por) / (5 * qsMax), dtDif = inp.beta > 0 ? 0.2 * dx * dx * (1 - inp.por) / (inp.beta * qsMax) : Infinity;
    const estSteps = tEnd / Math.min(dtAdv, dtDif, tEnd / 20);
    need(estSteps * n < 2e6, `The run needs about ${Math.round(estSteps)} time steps on ${n} cells: the bed reworks itself in roughly ${f(Math.min(dtAdv, dtDif) * 3000 / 86400, 2)} days at this transport rate. Shorten the duration or coarsen the grid.`);
    const res = run(() => exnerEvolve({ dx, zb0, tEnd, porosity: inp.por, scheme: 'laxWendroff', limiter: 'vanleer',
      // bed-slope effect: transport reduced going up-slope, increased down-slope, qs (1 − β ∂zb/∂x) — the physical
      // diffusion that fills a trench; the explicit step is limited to the diffusion stability bound
      transport: (zb) => Float64Array.from(zb, (z, i) => qsOf(Math.max(-z, 0.05 * inp.h0)) * (1 - inp.beta * (zb[Math.min(i + 1, n - 1)] - zb[Math.max(i - 1, 0)]) / (2 * dx))),
      dtMax: Math.max(1, Math.min(tEnd / 20, inp.beta > 0 ? 0.2 * dx * dx * (1 - inp.por) / (inp.beta * qsMax) : tEnd)),
      left: Math.abs(inp.supply - 1) < 1e-9 ? 'zeroGradient' : { q: inp.supply * qs0 }, right: 'zeroGradient', maxSteps: 400000 }), 'Bed evolution');
    const dev0 = Array.from(zb0, (z) => z + inp.h0), dev = Array.from(res.zb, (z) => z + inp.h0);
    const ext0 = sgn < 0 ? -Math.min(...dev0) : Math.max(...dev0), ext = sgn < 0 ? -Math.min(...dev) : Math.max(...dev);
    const infill = (1 - ext / ext0) * 100; // % of the initial trench depth (or mound height) removed at the deepest / highest point
    // migration of the feature: centroid of the bed perturbation (same sign as the feature)
    const cen = (a) => { let m = 0, mx = 0; a.forEach((v, i) => { const w = Math.max(sgn * v, 0); m += w; mx += w * xs[i]; }); return m > 0 ? mx / m : x0; };
    const celNum = (cen(dev) - cen(dev0)) / tEnd;
    // linear kinematic-wave celerity of a small bed perturbation: c = (dqs/dzb)/(1 − p) = 5 qs/(h (1 − p)) for qs ∝ U⁵
    const celTh = 5 * qs0 / (inp.h0 * (1 - inp.por));
    const volIn = (ext0 - ext) / ext0 * inp.amp * sig * Math.sqrt(2 * Math.PI); // m³ per m width (Gaussian area change)
    // time for the tolerable infill by rerunning is costly: scale from the simulated rate (first-order estimate)
    const tFill = infill > 0.5 ? inp.days * inp.allowFill / infill : 9999;
    const dredgeVol = Math.max(0, volIn) * inp.reachW * (365 / inp.days), dredgeDays = dredgeVol / inp.dredgeRate;

    // shoreline: Dean profile to the closure depth and Bruun retreat for each scenario
    const ws = settlingVelocitySoulsby({ d, rhoS: inp.rhoS });
    const A = deanParameter(ws), hc = depthOfClosure({ He: inp.He, Te: inp.Te });
    need(hc > 0, 'The depth of closure is not positive for this wave height and period: check the effective wave climate.');
    const Lp = deanProfile({ y: [0, 1], A }).widthToDepth(hc);
    const retreatOf = (S) => (S > 0 ? bruunRule({ S, hClosure: hc, B: inp.berm, L: Lp }).retreat : 0);
    const retreat = retreatOf(inp.slr);
    const scen = slrRows(inp).map((s) => ({ name: s.name, rise: s.rise, retreat: r(retreatOf(s.rise), 3), left: r(inp.beachW - retreatOf(s.rise), 3) }));
    const ys = linspace(0, Lp, 40), prof = deanProfile({ y: ys, A });
    const [px, p0, p1] = thin([xs, Array.from(zb0), Array.from(res.zb)], 400);
    const thetaCr = criticalShearStress({ d, rhoS: inp.rhoS });
    const okFill = infill <= inp.allowFill, okBeach = retreat <= inp.beachW || inp.beachW === 0;

    const advice = [
      okFill
        ? eng('ok', `After ${f(inp.days)} days the ${inp.feature} has lost ${f(infill, 2)} % of its ${f(inp.amp)} m ${sgn < 0 ? 'depth' : 'height'} (tolerable ${f(inp.allowFill)} %); it migrates at ${f(celNum * 86400, 2)} m/day.`,
          `Maintenance is due after roughly ${f(tFill, 2)} days at this rate. Confirm with repeated bathymetric surveys and calibrate the transport rate (${f(qs0 * 86400, 2)} m³/m/day) against the measured volume change.`,
          'Soulsby (1997) ch. 11; PIANC Report 102 (2008)', 'A surveyed, calibrated infill rate lets dredging be scheduled only when needed, avoiding unnecessary campaigns and their seabed disturbance.')
        : eng('fail', `After ${f(inp.days)} days the ${inp.feature} has lost ${f(infill, 2)} % of its ${f(inp.amp)} m ${sgn < 0 ? 'depth' : 'height'}, beyond the tolerable ${f(inp.allowFill)} %; the limit is reached after about ${f(tFill, 2)} days.`,
          `Over-dredge by the expected infill, widen the trench so that the flow slows less abruptly, backfill promptly after pipe-lay, or add a sediment trap upstream; plan maintenance dredging every ${f(tFill, 2)} days otherwise.`,
          'Soulsby (1997) ch. 11; PIANC Report 102 (2008); CIRIA C683 §4.6', `Sustained maintenance would move about ${f(dredgeVol, 3)} m³ a year; each dredging day burns of the order of ${f(MGO_T_PER_DAY, 2)} t of fuel (≈ ${f(MGO_T_PER_DAY * CO2_MGO, 2)} t CO₂) and re-disturbs the benthic habitat.`),
      eng(okBeach ? 'info' : 'warn', `Shoreline: with ${f(inp.slr)} m of sea-level rise the Bruun rule gives ${f(retreat, 3)} m of retreat (closure depth ${f(hc, 2)} m, active profile ${f(Lp, 3)} m wide) against ${f(inp.beachW)} m of beach in front of the structure.`,
        okBeach ? 'The beach remains in front of the structure for the design allowance; monitor profiles yearly and compare with the scenario table to time any nourishment.' : `The beach is lost before the design allowance is reached: allow for toe scour and higher waves at the wall, or nourish about ${f((retreat - inp.beachW) * (hc + inp.berm), 3)} m³ per metre of frontage to hold the line.`,
        `Bruun (1962); Hallermeier (1978); ${SLR_BASIS}`, 'Beach nourishment keeps the natural wave-dissipating buffer and habitat in front of the wall; losing the beach means heavier armour and a higher crest instead.'),
      eco('info', `At the simulated rate the ${inp.feature} takes in about ${f(Math.max(0, volIn), 3)} m³ per metre width in ${f(inp.days)} days — about ${f(dredgeVol, 3)} m³/year over the ${f(inp.reachW)} m reach, or ${f(dredgeDays, 2)} dredger days a year.`,
        'Weigh capital over-dredging against recurring maintenance on a whole-life basis, and schedule campaigns from survey triggers rather than fixed intervals.', 'Whole-life costing per ISO 15686-5; PIANC Report 102 (2008)',
        `Annual maintenance dredging would burn roughly ${f(dredgeDays * MGO_T_PER_DAY, 2)} t of marine gas oil (≈ ${f(dredgeDays * MGO_T_PER_DAY * CO2_MGO, 2)} t CO₂, screening factors).`),
    ];

    return {
      kpis: [
        { key: 'infill', label: `${sgn < 0 ? 'Trench infill' : 'Mound lowering'} at the end of the run`, value: r(infill), unit: '%', limit: inp.allowFill, sense: 'max', explain: 'Share of the initial depth (or height) of the feature lost at its deepest (highest) point. 100 % means the bed has returned to its undisturbed level.' },
        { key: 'ext', label: 'Remaining depth / height of the feature', value: r(ext), unit: 'm', explain: 'Largest remaining deviation from the undisturbed bed. For a pipe trench it is the cover depth still available.' },
        { key: 'tFill', label: 'Time to the tolerable infill', value: r(tFill, 3), unit: 'days', explain: 'First-order estimate of when maintenance is needed, scaled from the simulated infill rate. 9999 means no measurable infill.' },
        { key: 'qs0', label: 'Transport rate over the undisturbed bed', value: r(qs0 * 86400), unit: 'm³/m/day', explain: 'Solid volume of sand carried per metre width per day (Engelund–Hansen). This is the sediment supply that fills the trench; calibrate it against surveys.' },
        { key: 'cel', label: 'Migration speed of the feature', value: r(celNum * 86400), unit: 'm/day', explain: 'Speed at which the centroid of the trench or mound moves downstream. A pipe laid in a migrating trench can be left exposed on the upstream side.' },
        { key: 'celRatio', label: 'Numerical / theoretical bed-wave celerity', value: r(celNum / celTh), unit: '-', explain: 'Check of the solver against the kinematic-wave celerity 5 qs/[h(1 − p)]. Near 1 for small features; larger features deform and travel at a different speed.' },
        { key: 'massErr', label: 'Sediment volume balance error', value: Math.abs(res.mass.error) < 1e-9 ? 0 : r(Math.abs(res.mass.error), 3), unit: 'm³/m', explain: 'Change of bed volume minus the net boundary flux. Zero confirms the scheme conserves sediment; an error below 1e-9 m³ per metre of shoreline (a cubic millimetre) is arithmetic round-off and is reported as 0.' },
        { key: 'hc', label: 'Depth of closure', value: r(hc), unit: 'm', explain: 'Seaward limit of significant profile change (Hallermeier). Structures and nourishment designs use it as the base of the active beach.' },
        { key: 'retreat', label: 'Shoreline retreat for the design sea-level rise', value: r(retreat), unit: 'm', limit: inp.beachW, sense: 'max', explain: 'Equilibrium landward shift of the shoreline as the profile rises with the sea (Bruun rule). Beyond the available beach width the structure toe is exposed.' },
        { key: 'deanA', label: 'Dean profile parameter A', value: r(A), unit: 'm^⅓', explain: 'Scale of the equilibrium profile h = A y^⅔; coarser sand gives a larger A and a steeper beach.' },
      ],
      checks: [
        check('Infill against the maintenance trigger', Math.max(infill, 0), inp.allowFill, '%', 'PIANC Report 102 (2008); project dredging tolerance', 'Guards against a trench or channel silting beyond its tolerance before the pipe is laid or before the next survey; exceeded, the pipe cannot reach grade or the intake draws sediment.'),
        check('Shoreline retreat against the available beach', retreat, Math.max(inp.beachW, 1e-6), 'm', 'Bruun (1962); USACE CEM Part III-3', 'Guards against loss of the beach in front of the structure, which lets larger waves reach the wall and starts toe scour.'),
      ],
      plots: [
        { type: 'line', title: 'Bed level before and after', x: { label: 'Distance along the flow', unit: 'm' }, y: { label: 'Bed level', unit: 'm' },
          series: [{ name: 'Initial', x: px.map((v) => r(v)), y: p0.map((v) => r(v)), dash: true }, { name: `After ${f(inp.days)} days`, x: px.map((v) => r(v)), y: p1.map((v) => r(v)) }],
          explain: 'The bed along the flow at the start and end of the run. A trench fills from the upstream slope and migrates downstream; look at the remaining depth and at where the upstream edge has moved.' },
        { type: 'line', title: 'Equilibrium beach profile and closure depth', x: { label: 'Distance offshore', unit: 'm' }, y: { label: 'Depth', unit: 'm' },
          series: [{ name: 'Dean profile', x: ys.map((v) => r(v)), y: Array.from(prof.h, (v) => -r(v)) }, { name: 'Depth of closure', x: [0, Lp], y: [-hc, -hc], dash: true }],
          explain: 'The h = A y^⅔ profile out to the depth of closure. Its width and the berm height set how many metres the shoreline retreats per metre of sea-level rise.' },
        { type: 'line', title: 'Shoreline retreat versus sea-level rise', x: { label: 'Sea-level rise', unit: 'm' }, y: { label: 'Retreat', unit: 'm' },
          series: [{ name: 'Bruun rule', x: scen.map((s) => s.rise), y: scen.map((s) => s.retreat) }, { name: 'Beach width available', x: [scen[0].rise, scen[scen.length - 1].rise], y: [inp.beachW, inp.beachW], dash: true }],
          explain: 'Retreat for each scenario against the beach width in front of the structure. The crossing is the sea level at which the wall becomes the shoreline.' },
      ],
      tables: [
        { title: 'Sea-level-rise scenarios: shoreline', columns: [{ key: 'name', label: 'Scenario' }, { key: 'rise', label: 'Rise', unit: 'm' }, { key: 'retreat', label: 'Retreat', unit: 'm' }, { key: 'left', label: 'Beach remaining', unit: 'm' }],
          rows: scen, explain: 'Equilibrium retreat and the beach left in front of the structure for each scenario; negative remaining width means the toe is exposed.' },
        { title: 'Run diagnostics', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }],
          rows: [{ q: 'Time steps', v: res.steps }, { q: 'Maximum Courant number', v: f(res.courantMax) }, { q: 'Cell size (m)', v: f(dx) }, { q: 'Chézy coefficient (m^½/s)', v: f(C) }, { q: 'Critical Shields number', v: f(thetaCr.thetaCr) }, { q: 'Boundary net inflow (m³/m)', v: f(res.mass.boundaryNet) }, { q: 'Theoretical celerity (m/day)', v: f(celTh * 86400) }],
          explain: 'Numerical and physical diagnostics of the bed-evolution run, for the verification record (conservation, stability, resolution).' },
      ],
      notes: ['Transport ∝ U⁵ (Engelund–Hansen) with flow continuity and a bed-slope correction: deeper water over the trench slows the flow and drops sediment, and the side slopes relax.', `Sediment balance error ${f(Math.abs(res.mass.error), 2)} m³/m over ${res.steps} steps.`],
      advice,
      publish: { 'marine.shorelineRetreat': retreat, 'marine.infillRate': Math.max(0, volIn) / tEnd },
      quantities: { vessel_days: { 'maintenance dredging per year': r(dredgeDays, 3) } },
    };
  },
  verify: [
    { name: 'Small bed wave migrates at the kinematic celerity c = (dqs/dzb)/(1 − p) = 5 qs/[h (1 − p)] and conserves sediment',
      source: 'Linearised Exner equation with qs ∝ U⁵ and U = q/h (Soulsby 1997 §11; LeVeque 2002 §6): a small perturbation translates unchanged at c. Independent of the numerical scheme',
      inputs: { feature: 'mound', amp: 0.02, width: 200, h0: 8, U0: 0.9, days: 200, supply: 1, nx: 400 },
      expect: [{ kpi: 'celRatio', value: 1, tol: 0.03 }, { kpi: 'massErr', value: 0, abs: 1e-6 }] },
    { name: 'Bruun rule R = S L/(B + h*)', source: 'Bruun (1962). With the module\'s Dean profile width L and closure depth h*: Hallermeier h* = 2.28 He − 68.5 He²/(g Te²) = 2.28 × 3.5 − 68.5 × 12.25/(9.80665 × 81) = 6.924 m',
      inputs: { He: 3.5, Te: 9 }, expect: [{ kpi: 'hc', value: 6.924, tol: 0.002 }] },
  ],
};
