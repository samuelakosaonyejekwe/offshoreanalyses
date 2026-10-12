// Flooding behind a coastal defence: weir overflow and wave overtopping over the crest routed across the
// hinterland with the 1-D nonlinear shallow-water equations (wetting and drying), with volume balance.
import { solveSwe1d } from '../../engine/coastal/swe.js';
import { G, need, run, finite, r, f, linspace, check, eng, eco, slrInputs, slrRows, SLR_BASIS, CO2_CONC } from './_common.js';

export default {
  id: 's07-flood-inundation', suite: 7, title: 'Overflow flooding of the hinterland',
  summary: 'Routes the water that overflows and overtops a defence across the land behind it: flood depth, velocity, arrival time and extent at a receptor, with a closed volume balance.',
  covers: [10], fidelity: 'numerical',
  standards: ['EurOtop Manual (2018) ch. 3 (overtopping consequences)', 'EN 1991-1-6 / ISO 21650:2007', 'Defra/EA FD2320 Flood Risks to People (2006) hazard rating', 'USACE EM 1110-2-1100 Part VI-5-2'],
  models: ['Nonlinear shallow-water equations (continuity and momentum)', 'HLL finite-volume scheme with wetting and drying', 'Broad-crested weir flow (Bernoulli, critical depth)', 'Manning friction', 'Level-pool storage (continuity)'],
  limits: 'One-dimensional (a flow strip normal to the defence): no lateral spreading, no street-scale 2-D flow paths, hydrostatic and depth-averaged, so breaking-wave jets at the crest are represented only by their volume. Use it to bound depths and arrival times; a 2-D flood model or 3-D CFD is outside the module and is needed for mapping.',
  publishes: ['marine.floodDepth', 'marine.floodVolume'],
  searchTerms: 'coastal flooding overflow overtopping inundation shallow water equations hazard',
  inputs: [
    { key: 'crest', label: 'Defence crest level above datum', unit: 'm', type: 'number', default: 4.2, min: 0, max: 30, group: 'Defence', role: 'geometry', design: true },
    { key: 'swl', label: 'Still-water level at the start of the event', unit: 'm', type: 'number', default: 2.6, min: -5, max: 20, group: 'Sea boundary', role: 'bc', help: 'Tide level above datum on which the surge rides.' },
    { key: 'surge', label: 'Peak storm surge', unit: 'm', type: 'number', default: 1.5, min: 0, max: 10, group: 'Sea boundary', role: 'bc', bind: 'site.surge100', uncertain: { dist: 'lognormal', cov: 0.2 } },
    ...slrInputs({ group: 'Sea boundary' }),
    { key: 'durMin', label: 'Duration of the surge peak', unit: 'min', type: 'number', default: 40, min: 2, max: 360, group: 'Sea boundary', role: 'bc', help: 'The sea level follows swl + surge·sin²(πt/T) over this duration.' },
    { key: 'q0', label: 'Peak wave-overtopping discharge', unit: 'l/s/m', type: 'number', default: 20, min: 0, max: 2000, group: 'Sea boundary', role: 'bc', help: 'From the overtopping module at the peak water level; follows the same time shape as the surge.' },
    { key: 'Cw', label: 'Weir coefficient q/(√g h^1.5)', unit: '-', type: 'number', default: 0.544, min: 0.3, max: 0.7, group: 'Defence', role: 'criteria', help: 'Broad-crested weir with critical flow: (2/3)^1.5 = 0.544. A calibration parameter against measured overflow.' },
    { key: 'L', label: 'Length of the flow path inland', unit: 'm', type: 'number', default: 500, min: 20, max: 20000, group: 'Hinterland', role: 'geometry' },
    { key: 'z0', label: 'Ground level just behind the defence', unit: 'm', type: 'number', default: 1.5, min: -10, max: 30, group: 'Hinterland', role: 'ic' },
    { key: 'slope', label: 'Ground slope inland (+ rising, − falling)', unit: 'm/km', type: 'number', default: -1, min: -50, max: 50, group: 'Hinterland', role: 'geometry', help: 'Falling ground behind a defence (a polder) ponds water at the landward end.' },
    { key: 'h0', label: 'Initial standing water depth', unit: 'm', type: 'number', default: 0, min: 0, max: 5, group: 'Hinterland', role: 'ic', help: 'Initial condition: 0 for dry ground, or the depth already ponded by rainfall.' },
    { key: 'manning', label: 'Manning roughness n', unit: 's/m^⅓', type: 'number', default: 0.035, min: 0.01, max: 0.2, group: 'Hinterland', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.25 }, help: 'Bed boundary: 0.02 paved, 0.035 grass, 0.06–0.1 built-up or vegetated.' },
    { key: 'landward', label: 'Landward boundary', type: 'select', default: 'wall', options: [{ value: 'wall', label: 'Closed (rising ground, second defence): water accumulates' }, { value: 'transmissive', label: 'Open: water drains away inland' }], group: 'Hinterland', role: 'bc' },
    { key: 'xR', label: 'Receptor distance from the defence', unit: 'm', type: 'number', default: 300, min: 1, max: 20000, group: 'Receptor', role: 'geometry' },
    { key: 'dTol', label: 'Tolerable flood depth at the receptor', unit: 'm', type: 'number', default: 0.3, min: 0.01, max: 5, group: 'Receptor', role: 'criteria', help: 'Typically the threshold level of buildings or the wading limit of vehicles.' },
    { key: 'hrTol', label: 'Tolerable hazard rating d (v + 0.5)', unit: 'm²/s', type: 'number', default: 0.75, min: 0.1, max: 5, group: 'Receptor', role: 'criteria', help: 'FD2320: below 0.75 low hazard, 0.75–1.25 danger for some, above 1.25 danger for most (debris factor not included).' },
    { key: 'front', label: 'Length of defence overflowing', unit: 'm', type: 'number', default: 300, min: 1, max: 50000, group: 'Defence', role: 'geometry' },
    { key: 'nx', label: 'Grid cells along the flow path', type: 'int', default: 50, min: 10, max: 800, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nx', label: 'Grid cells', levels: [25, 50, 100], qoi: 'depthR' },
  compute(inp) {
    const T = inp.durMin * 60, n = Math.round(inp.nx), dx = inp.L / n;
    need(inp.xR < inp.L, 'The receptor must lie within the modelled flow path: increase the path length or move the receptor.');
    need(inp.z0 < inp.crest, 'The ground behind the defence is above its crest: there is nothing to flood.');
    const shape = (t) => (t <= 0 || t >= T ? 0 : Math.sin(Math.PI * t / T) ** 2);
    const qOf = (t, rise) => { const head = inp.swl + rise + inp.surge * shape(t) - inp.crest; return (head > 0 ? inp.Cw * Math.sqrt(G) * head ** 1.5 : 0) + inp.q0 / 1000 * shape(t); };
    const volOf = (rise) => { let v = 0; const m = 400; for (let i = 0; i < m; i++) v += qOf((i + 0.5) * T / m, rise) * T / m; return v; };
    const peakOf = (rise) => qOf(T / 2, rise);
    const Vin = volOf(inp.slr);
    need(inp.swl + inp.slr < inp.crest, `The still-water level including sea-level rise (${f(inp.swl + inp.slr)} m) is above the crest before the surge arrives: the defence is permanently overflowed. Raise the crest.`);
    const zOf = (x) => inp.z0 + inp.slope / 1000 * x;
    const tEnd = 1.5 * T; // run on after the peak so that the flood wave reaches the far end

    let res = null, x = [], scale = 1, tR = [], hR = [], iR = 0;
    if (Vin > 1e-6) {
      const solve = (cells, qScale, outputs) => {
        const d = inp.L / cells;
        const zb = Float64Array.from({ length: cells }, (_, i) => zOf((i + 0.5) * d));
        return solveSwe1d({ dx: d, zb, h0: new Float64Array(cells).fill(inp.h0), tEnd, manning: inp.manning, left: { type: 'discharge', q: (t) => qScale * qOf(t, inp.slr) }, right: { type: inp.landward }, outputTimes: outputs, hShore: 0.005, maxSteps: 60000 });
      };
      // the characteristic discharge boundary delivers a few per cent more than prescribed on an initially dry bed:
      // a pre-run measures the ratio and the inflow is rescaled so that the routed volume matches the hydrograph
      if (inp.landward === 'wall') { const pre = run(() => solve(n, 1, linspace(tEnd / 1200, tEnd, 1200)), 'Flood routing'); if (pre.volume.boundaryNet > 0) scale = Vin / pre.volume.boundaryNet; }
      // output times double as a step cap: on a dry bed the solver has no wave speed to limit its time step, so
      // without them the first step would swallow the whole hydrograph
      const outT = linspace(tEnd / 1200, tEnd, 1200);
      res = run(() => solve(n, scale, outT), 'Flood routing (shorten the event, coarsen the grid or raise the roughness if the step limit is hit)');
      x = Array.from(res.x);
      iR = Math.min(n - 1, Math.floor(inp.xR / dx));
      for (let i = 19; i < res.snapshots.length; i += 20) { tR.push(res.snapshots[i].t); hR.push(res.snapshots[i].h[iR]); }
    }
    const depthR = res ? Math.max(res.envelope.hMax[iR] - inp.h0, 0) + inp.h0 : inp.h0;
    const vR = res ? res.envelope.uMax[iR] : 0;
    // hazard rating from the snapshot with the largest d(v + 0.5) at the receptor
    let hr = depthR > 0.005 ? depthR * 0.5 : 0, arrive = -1;
    if (res) res.snapshots.forEach((s) => { const d = s.h[iR], v = Math.abs(s.u[iR]); if (d > 0.005) hr = Math.max(hr, d * (v + 0.5)); if (arrive < 0 && d - inp.h0 > 0.01) arrive = s.t; });
    const extent = res ? Math.min(inp.L, Math.max(0, res.envelope.inundationLimit)) : 0;
    const stored = res ? res.volume.final - res.volume.initial : 0;
    const delivered = res ? res.volume.boundaryNet : 0;
    const volErr = res ? Math.abs(res.volume.error) : 0;
    const depthMaxAll = res ? Math.max(...res.envelope.hMax) : inp.h0;
    finite(depthR, 'Flood depth');
    // crest that stops overflow at the surge peak (overtopping remains)
    const crestNoOverflow = inp.swl + inp.slr + inp.surge;
    const raise = Math.max(0, crestNoOverflow - inp.crest);
    const scen = slrRows(inp).map((s) => { const sub = inp.swl + s.rise >= inp.crest; const v = sub ? 0 : volOf(s.rise); return { name: s.name, rise: s.rise, head: r(Math.max(0, inp.swl + s.rise + inp.surge - inp.crest), 3), qPeak: r(peakOf(s.rise) * 1000, 3), vol: r(v), pool: sub ? 'permanently overflowed' : f(v / inp.L, 3) }; });
    const ok = depthR <= inp.dTol && hr <= inp.hrTol;
    const concRaise = raise * 0.5 * inp.front;

    const advice = [];
    advice.push(ok
      ? eng('ok', `At the receptor ${f(inp.xR)} m inland the flood reaches ${f(depthR, 2)} m (tolerable ${f(inp.dTol)} m) with hazard rating ${f(hr, 2)} m²/s; ${f(Vin, 3)} m³ per metre of defence enters over the event.`,
        'The receptor stays within tolerable limits for this event. Check the higher sea-level scenarios in the table and confirm flow paths with a 2-D model where streets or channels concentrate the flow.', 'Defra/EA FD2320 (2006) hazard rating; EurOtop (2018) ch. 3',
        'Keeping floodwater shallow at the receptor avoids contaminated-water damage and the waste from stripping out flooded buildings.')
      : eng('fail', `At the receptor ${f(inp.xR)} m inland the flood reaches ${f(depthR, 2)} m (tolerable ${f(inp.dTol)} m), velocity ${f(vR, 2)} m/s, hazard rating ${f(hr, 2)} m²/s (tolerable ${f(inp.hrTol)}); it arrives ${arrive >= 0 ? `${f(arrive / 60, 2)} min` : 'within the event'} after the start.`,
        `Raise the crest by ${f(raise, 2)} m to ${f(crestNoOverflow)} m to stop overflow at the surge peak, or provide ${f(Vin * inp.front, 3)} m³ of flood storage or equivalent pumping, or set a secondary defence line before the receptor; use the arrival time for the warning and evacuation plan.`,
        'Defra/EA FD2320 (2006); EurOtop (2018) ch. 3; EN 1997-1 / CIRIA C731 International Levee Handbook (2013) §8',
        `A ${f(raise, 2)} m crest wall over ${f(inp.front)} m takes about ${f(concRaise, 3)} m³ of concrete (≈ ${f(concRaise * CO2_CONC, 3)} t CO₂e, screening factor); set against repeated flood damage and waste, raising is usually the lower-carbon path.`));
    const backLevel = res ? zOf(0.5 * dx) + res.envelope.hMax[0] : inp.z0;
    if (backLevel > Math.max(inp.crest, inp.swl + inp.slr + inp.surge) + 0.05) advice.push(eng('warn', `The water behind the defence rises to ${f(backLevel)} m, above the crest and the peak sea level: in reality the inflow would be drowned and water would return seaward.`, 'Treat the depths as an upper bound; lengthen the flow path, open the landward boundary or reduce the inflow to represent the drowned weir.', 'Weir submergence (Chow 1959 §14); level-pool continuity'));
    const worst = scen[scen.length - 1];
    advice.push(eng(worst.head > 0 ? 'warn' : 'info', `Sea-level rise: at +${f(worst.rise)} m (${worst.name}) the peak overflow head is ${f(worst.head, 2)} m and the peak inflow ${f(worst.qPeak)} l/s/m; overflow volume grows with head to the power 1.5.`,
      'Fix the crest for the adopted end-of-life scenario and keep the option of a demountable or added wall for the higher ones; review at each update of the projections.', SLR_BASIS));
    advice.push(eco(ok ? 'info' : 'warn', `The event puts ${f(Vin * inp.front, 3)} m³ over ${f(inp.front)} m of defence; ${f(stored * inp.front, 3)} m³ remains on the strip at the end of the run.`,
      'Value the options by expected annual damage avoided: crest raising removes the overflow term entirely, storage and pumping only manage it, and warning time reduces loss of life but not property damage.', 'Expected-annual-damage appraisal (national flood-risk appraisal guidance); whole-life costing per ISO 15686-5',
      'Pumped drainage of the stored volume costs energy on every event; a passive crest raise has no operating energy.'));

    const snaps = res ? [299, 599, 899, 1199].map((i) => res.snapshots[Math.min(i, res.snapshots.length - 1)]) : [];
    return {
      kpis: [
        { key: 'depthR', label: 'Maximum flood depth at the receptor', value: r(depthR), unit: 'm', limit: inp.dTol, sense: 'max', explain: 'Largest water depth reached at the receptor during the event. Compare with building thresholds and wading limits.' },
        { key: 'vR', label: 'Maximum velocity at the receptor', value: r(vR, 3), unit: 'm/s', explain: 'Largest depth-averaged flow speed at the receptor. Above about 1–2 m/s people and vehicles are swept away even in shallow water.' },
        { key: 'hazard', label: 'Hazard rating d (v + 0.5) at the receptor', value: r(hr, 3), unit: 'm²/s', limit: inp.hrTol, sense: 'max', explain: 'Flood-hazard-to-people index combining depth and velocity (debris factor excluded). Below 0.75 low hazard; above 1.25 danger for most.' },
        { key: 'arrival', label: 'Flood arrival time at the receptor', value: r(arrive >= 0 ? arrive / 60 : 0, 3), unit: 'min', explain: 'Time from the start of the event until the water at the receptor has risen by 10 mm: the available warning time. Zero means the flood does not arrive.' },
        { key: 'extent', label: 'Inundation extent inland', value: r(extent), unit: 'm', explain: 'Furthest point reached by the wet front along the flow path. Equal to the path length when the whole strip is flooded.' },
        { key: 'volIn', label: 'Volume entering per metre of defence', value: r(Vin), unit: 'm³/m', explain: 'Time integral of weir overflow plus wave overtopping over the event — the inflow boundary condition of the flood.' },
        { key: 'qPeak', label: 'Peak inflow', value: r(peakOf(inp.slr) * 1000), unit: 'l/s/m', explain: 'Largest discharge over the crest at the surge peak. Once the sea is above the crest, weir flow dwarfs wave overtopping.' },
        { key: 'head', label: 'Peak overflow head above the crest', value: r(Math.max(0, inp.swl + inp.slr + inp.surge - inp.crest), 3), unit: 'm', explain: 'How far the peak still-water level stands above the crest. Zero means no overflow, only wave overtopping.' },
        { key: 'depthMax', label: 'Maximum depth anywhere on the path', value: r(depthMaxAll), unit: 'm', explain: 'Deepest water along the strip, usually at the low point or against the landward boundary.' },
        { key: 'volRouted', label: 'Volume delivered by the solver', value: r(delivered), unit: 'm³/m', explain: 'Net volume that crossed the boundaries in the simulation. With a closed landward boundary it should equal the prescribed inflow.' },
        { key: 'volErr', label: 'Water-volume balance error', value: r(volErr, 3), unit: 'm³/m', explain: 'Stored volume minus net boundary inflow. Machine-zero confirms that the routing neither creates nor loses water.' },
        { key: 'crestReq', label: 'Crest level that stops overflow', value: r(crestNoOverflow), unit: 'm', explain: 'Peak still-water level including surge and the sea-level-rise allowance. A crest at this level leaves only wave overtopping.' },
      ],
      checks: [
        check('Flood depth at the receptor', depthR, inp.dTol, 'm', 'Project flood-protection criterion; EurOtop (2018) ch. 3', 'Guards against water entering buildings or stopping vehicles at the receptor. Exceeded, property floods and access is cut.'),
        check('Hazard to people at the receptor', hr, inp.hrTol, 'm²/s', 'Defra/EA FD2320 (2006) Table 13.1', 'Guards against people being knocked over or swept away by the combination of depth and velocity.'),
      ],
      plots: res ? [
        { type: 'line', title: 'Water surface along the flow path', x: { label: 'Distance inland', unit: 'm' }, y: { label: 'Level above datum', unit: 'm' },
          series: [{ name: 'Ground', x: x.map((v) => r(v)), y: x.map((v) => r(zOf(v), 3)), dash: true }, ...snaps.map((s) => ({ name: `t = ${f(s.t / 60, 2)} min`, x: x.map((v) => r(v)), y: x.map((v, i) => r(zOf(v) + s.h[i], 3)) }))],
          explain: 'Flood surface at four times. Look for the advancing front, the bore reflecting from the landward boundary and where the water finally ponds.' },
        { type: 'area', title: 'Maximum depth and velocity envelope', x: { label: 'Distance inland', unit: 'm' }, y: { label: 'Depth (m) / velocity (m/s)', unit: '' },
          series: [{ name: 'Maximum depth', x: x.map((v) => r(v)), y: Array.from(res.envelope.hMax, (v) => r(v, 3)) }, { name: 'Maximum velocity', x: x.map((v) => r(v)), y: Array.from(res.envelope.uMax, (v) => r(v, 3)) }, { name: 'Tolerable depth', x: [0, inp.L], y: [inp.dTol, inp.dTol], dash: true }],
          explain: 'Largest depth and speed reached at each point over the whole event. Velocity is highest just behind the defence; depth is greatest where the water ponds.' },
        { type: 'line', title: 'Inflow hydrograph and receptor depth', x: { label: 'Time', unit: 'min' }, y: { label: 'Inflow (l/s/m ÷ 100) / depth (m)', unit: '' },
          series: [{ name: 'Inflow ÷ 100', x: linspace(0, tEnd, 60).map((t) => r(t / 60, 3)), y: linspace(0, tEnd, 60).map((t) => r(qOf(t, inp.slr) * 10, 3)) }, { name: 'Depth at receptor', x: tR.map((t) => r(t / 60, 3)), y: hR.map((v) => r(v, 3)) }],
          explain: 'The discharge over the crest and the response at the receptor. The lag between the inflow peak and the depth peak is the travel time of the flood wave.' },
      ] : [{ type: 'line', title: 'Inflow hydrograph', x: { label: 'Time', unit: 'min' }, y: { label: 'Inflow', unit: 'l/s/m' }, series: [{ name: 'Inflow', x: [0, inp.durMin], y: [0, 0] }], explain: 'No water crosses the crest in this event: the hydrograph is zero and the hinterland stays dry.' }],
      tables: [
        { title: 'Sea-level-rise scenarios: inflow', columns: [{ key: 'name', label: 'Scenario' }, { key: 'rise', label: 'Rise', unit: 'm' }, { key: 'head', label: 'Peak head', unit: 'm' }, { key: 'qPeak', label: 'Peak inflow', unit: 'l/s/m' }, { key: 'vol', label: 'Volume', unit: 'm³/m' }, { key: 'pool', label: 'Level-pool depth (m)' }],
          rows: scen, explain: 'Overflow head, peak discharge and event volume for each scenario, with the mean depth if that volume ponded over the flow path. The steep growth marks the scenario at which the defence stops working.' },
        { title: 'Volume balance', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'm³ per m' }],
          rows: [{ q: 'Prescribed inflow (hydrograph integral)', v: r(Vin) }, { q: 'Net boundary inflow in the solver', v: r(delivered) }, { q: 'Change of stored volume', v: r(stored) }, { q: 'Balance error', v: r(volErr, 3) }, { q: 'Inflow scale factor applied', v: r(scale, 4) }, { q: 'Time steps', v: res ? res.steps : 0 }],
          explain: 'Water-volume conservation of the run: inflow, storage and the residual. The scale factor corrects the small over-delivery of the characteristic discharge boundary on a dry bed.' },
      ],
      notes: [Vin > 0 ? `Peak head over the crest ${f(Math.max(0, inp.swl + inp.slr + inp.surge - inp.crest), 2)} m; weir flow provides ${f((1 - Math.min(1, inp.q0 / 1000 * T / 2 / Vin)) * 100, 2)} % of the volume.` : 'No overflow or overtopping in this event.',
        'One-dimensional strip: depths are an upper bound where the flow can spread sideways and a lower bound where it is funnelled.'],
      advice,
      publish: { 'marine.floodDepth': depthR, 'marine.floodVolume': Vin * inp.front },
      quantities: { concrete_m3: r(ok ? 0 : concRaise) },
    };
  },
  verify: [
    { name: 'Volume conservation: closed basin stores the hydrograph integral ∫ q0 sin²(πt/T) dt = q0 T/2',
      source: 'Continuity (analytic integral). Overtopping only, q0 = 100 l/s/m, T = 30 min: V = 0.1 × 1800/2 = 90 m³/m',
      inputs: { crest: 8, swl: 2, surge: 1, slr: 0, durMin: 30, q0: 100, landward: 'wall', slope: 0, h0: 0 },
      expect: [{ kpi: 'volIn', value: 90, tol: 0.001 }, { kpi: 'volRouted', value: 90, tol: 0.01 }, { kpi: 'volErr', value: 0, abs: 1e-6 }] },
    { name: 'Broad-crested weir: q = (2/3)^1.5 √g h^1.5 at the surge peak', source: 'Bernoulli with critical depth on the crest (Chow 1959 §3-7). Head 0.5 m: q = 0.5443 × 3.1316 × 0.35355 = 0.6026 m²/s = 602.6 l/s/m',
      inputs: { crest: 4, swl: 3, surge: 1.5, slr: 0, q0: 0, Cw: 0.5443, durMin: 10 }, expect: [{ kpi: 'qPeak', value: 602.6, tol: 0.002 }] },
  ],
};
