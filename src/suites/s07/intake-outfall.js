// Hydraulics of a sea outfall or intake conduit: head balance against the sea level (Bernoulli with friction
// and minor losses), gravity capacity, conduit and port velocities, screen approach velocity, pumping power
// and pump-trip surge, with the loss of gravity head under sea-level rise.
import { liquidLinePressureDrop, pumpPower, waveSpeed, joukowsky } from '../../engine/flow/singlephase.js';
import { G, CO2_CONC, need, run, finite, r, f, linspace, check, eng, eco, slrInputs, slrRows, SLR_BASIS } from './_common.js';

export default {
  id: 's07-intake-outfall', suite: 7, title: 'Intake and outfall hydraulics',
  summary: 'Head needed to drive the design flow through an outfall (or the drawdown of an intake) against the sea level, the velocities in the conduit, ports and screens, and how sea-level rise erodes gravity capacity.',
  covers: [10], fidelity: 'analytical',
  standards: ['EN 752:2017 / EN 16933-2 (hydraulic design of drains and outfalls)', 'WRc Design Guide for Marine Treatment Schemes (1990)', 'US EPA 40 CFR 125 Subpart J (intake through-screen velocity 0.15 m/s)', 'ANSI/HI 9.8 (pump intake design)', 'Wylie & Streeter, Fluid Transients in Systems (1993)'],
  models: ['Continuity equation', 'Bernoulli equation with losses', 'Darcy–Weisbach friction (Colebrook–White)', 'Minor-loss coefficients', 'Hydrostatic density head', 'Densimetric Froude number', 'Joukowsky water hammer'],
  limits: 'Steady full-pipe flow in a single conduit with lumped minor losses and an evenly distributing diffuser; port-by-port manifold hydraulics, air entrainment, sediment deposition, near-field and far-field dilution, wave-induced oscillation in the conduit and full transient (method-of-characteristics) surge analysis are not modelled. 3-D CFD of chambers and plumes is outside the module.',
  publishes: ['marine.outfallHead', 'marine.outfallCapacity'],
  searchTerms: 'sea outfall intake hydraulics head loss diffuser Froude number screen velocity sea level rise',
  inputs: [
    { key: 'mode', label: 'Structure', type: 'select', default: 'outfall', options: [{ value: 'outfall', label: 'Outfall (discharge to sea)' }, { value: 'intake', label: 'Intake (abstraction from sea)' }], group: 'Conduit', role: 'geometry' },
    { key: 'Q', label: 'Design flow', unit: 'm³/s', type: 'number', default: 2.5, min: 0.01, max: 200, group: 'Flow boundary', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Prescribed discharge or abstraction rate (peak wet-weather flow for an outfall).' },
    { key: 'route', label: 'Conduit centre-line (optional import)', type: 'geometry', kinds: ['polyline'], default: null, group: 'Conduit', role: 'geometry', help: 'A 3-D polyline of the alignment gives the length; leave empty to use the length below.' },
    { key: 'L', label: 'Conduit length', unit: 'm', type: 'number', default: 1200, min: 5, max: 30000, group: 'Conduit', role: 'geometry' },
    { key: 'D', label: 'Internal diameter', unit: 'm', type: 'number', default: 1.4, min: 0.1, max: 8, group: 'Conduit', role: 'geometry', design: true },
    { key: 'rough', label: 'Wall roughness', unit: 'mm', type: 'number', default: 1.5, min: 0.01, max: 30, group: 'Conduit', role: 'material', uncertain: { dist: 'lognormal', cov: 0.4 }, help: 'New concrete or lined steel 0.1–0.6; slimed or marine-fouled conduits 1.5–6. Calibrate against measured head loss.' },
    { key: 'K', label: 'Sum of minor-loss coefficients', unit: '-', type: 'number', default: 2.5, min: 0, max: 50, group: 'Conduit', role: 'criteria', help: 'Entry 0.5, bends 0.2–0.4 each, valves, chamber and riser losses; port exit loss is added separately.' },
    { key: 'nPorts', label: 'Number of diffuser ports / intake openings', type: 'int', default: 16, min: 1, max: 500, group: 'Terminal', role: 'geometry', design: true },
    { key: 'dPort', label: 'Port / opening diameter', unit: 'm', type: 'number', default: 0.3, min: 0.02, max: 5, group: 'Terminal', role: 'geometry', design: true },
    { key: 'Cdp', label: 'Port discharge coefficient', unit: '-', type: 'number', default: 0.65, min: 0.4, max: 1, group: 'Terminal', role: 'criteria' },
    { key: 'hOut', label: 'Water depth at the terminal', unit: 'm', type: 'number', default: 15, min: 0.5, max: 200, group: 'Receiving water boundary', role: 'bc', bind: 'site.depth' },
    { key: 'sea', label: 'Design sea level above datum', unit: 'm', type: 'number', default: 2.4, min: -10, max: 15, group: 'Receiving water boundary', role: 'bc', help: 'High water plus surge for an outfall (worst for gravity discharge); lowest water for an intake.' },
    ...slrInputs({ group: 'Receiving water boundary' }),
    { key: 'rhoSea', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 995, max: 1050, group: 'Receiving water boundary', role: 'environment', bind: 'site.rhoW' },
    { key: 'rhoE', label: 'Density of the conveyed water', unit: 'kg/m³', type: 'number', default: 999, min: 950, max: 1100, group: 'Flow boundary', role: 'bc', help: 'Fresh effluent ≈ 999; brine up to 1050; sea water for an intake.' },
    { key: 'zAvail', label: 'Available water level at the landward chamber', unit: 'm', type: 'number', default: 6.5, min: -20, max: 60, group: 'Flow boundary', role: 'bc', help: 'Outfall: highest level the headworks can stand at without flooding upstream. Intake: minimum level the pumps need in the sump.' },
    { key: 'screenA', label: 'Net open area of the intake screens', unit: 'm²', type: 'number', default: 20, min: 0.1, max: 2000, group: 'Terminal', role: 'geometry', showIf: { key: 'mode', in: ['intake'] } },
    { key: 'vScreen', label: 'Tolerable through-screen velocity', unit: 'm/s', type: 'number', default: 0.15, min: 0.05, max: 1, group: 'Criteria', role: 'criteria', showIf: { key: 'mode', in: ['intake'] } },
    { key: 'vMin', label: 'Minimum self-cleansing velocity', unit: 'm/s', type: 'number', default: 0.8, min: 0.2, max: 2, group: 'Criteria', role: 'criteria', help: 'About 0.6–1.0 m/s at least daily to keep grit and marine growth moving.' },
    { key: 'vMax', label: 'Maximum conduit velocity', unit: 'm/s', type: 'number', default: 3, min: 1, max: 8, group: 'Criteria', role: 'criteria' },
    { key: 'eff', label: 'Pump efficiency (wire to water)', unit: '-', type: 'number', default: 0.75, min: 0.3, max: 0.92, group: 'Pumping', role: 'operation' },
    { key: 'hours', label: 'Pumping hours per year', unit: 'h', type: 'number', default: 1500, min: 0, max: 8760, group: 'Pumping', role: 'operation' },
    { key: 'tClose', label: 'Pump run-down / valve closure time', unit: 's', type: 'number', default: 15, min: 0.1, max: 600, group: 'Pumping', role: 'operation' },
    { key: 'wall', label: 'Conduit wall thickness', unit: 'mm', type: 'number', default: 120, min: 3, max: 800, group: 'Conduit', role: 'geometry' },
    { key: 'Ewall', label: 'Wall Young\'s modulus', unit: 'GPa', type: 'number', default: 35, min: 0.5, max: 210, group: 'Conduit', role: 'material', help: 'Concrete 30–40, steel 207, HDPE about 1.' },
  ],
  compute(inp) {
    const outfall = inp.mode === 'outfall';
    let L = inp.L;
    if (inp.route && inp.route.kind === 'polyline' && inp.route.points && inp.route.points.length >= 6) { const p = inp.route.points; L = 0; for (let i = 3; i < p.length; i += 3) L += Math.hypot(p[i] - p[i - 3], p[i + 1] - p[i - 2], p[i + 2] - p[i - 1]); }
    need(L > 1, 'The imported centre-line has no length.');
    const Aport = inp.nPorts * Math.PI * inp.dPort ** 2 / 4, Apipe = Math.PI * inp.D ** 2 / 4;
    const mu = 1.1e-3;
    /** Total head loss [m of conveyed water] at flow Q. */
    const losses = (Q) => {
      const dp = liquidLinePressureDrop({ Q, D: inp.D, rho: inp.rhoE, mu, rough: inp.rough / 1000, L, K: inp.K });
      const vp = Q / (inp.Cdp * Aport);
      return { hf: dp.dPfriction / (inp.rhoE * G), hm: dp.dPminor / (inp.rhoE * G), hp: vp * vp / (2 * G), v: dp.v, fr: dp.f, Re: dp.Re, vPort: Q / Aport };
    };
    need(inp.Q / Apipe < 15, `The conduit velocity would be ${f(inp.Q / Apipe)} m/s: the diameter is far too small for this flow.`);
    const lo = run(() => losses(inp.Q), 'Conduit hydraulics');
    const hLoss = lo.hf + lo.hm + lo.hp;
    // density head: a column of lighter effluent must stand higher than the sea to balance it at the port depth
    const hDens = outfall ? (inp.rhoSea / inp.rhoE - 1) * inp.hOut : 0;
    const need_ = (rise, Q = inp.Q) => { const l = losses(Q); return inp.sea + rise + (outfall ? 1 : -1) * (l.hf + l.hm + l.hp) + hDens; };
    const zReq = finite(need_(inp.slr), 'Required head');
    // outfall: required chamber level ≤ available; intake: sump level (sea − losses) ≥ minimum
    const marginH = outfall ? inp.zAvail - zReq : zReq - inp.zAvail;
    const capacity = (rise) => {
      const m = (Q) => (outfall ? inp.zAvail - need_(rise, Q) : need_(rise, Q) - inp.zAvail);
      if (m(1e-4) <= 0) return 0;
      let a = 1e-4, b = Math.max(inp.Q, 0.1); let guard = 0;
      while (m(b) > 0 && guard++ < 40) b *= 1.6;
      for (let i = 0; i < 26; i++) { const mid = 0.5 * (a + b); if (m(mid) > 0) a = mid; else b = mid; }
      return a;
    };
    const Qcap = capacity(inp.slr);
    const pumpHead = Math.max(0, -marginH), pw = pumpHead > 0 ? pumpPower({ Q: inp.Q, head: pumpHead, rho: inp.rhoE, efficiency: inp.eff }).shaft : 0;
    const energy = pw * inp.hours / 1e6; // MWh per year
    const gPrime = G * Math.abs(inp.rhoSea - inp.rhoE) / inp.rhoSea;
    const Frd = gPrime > 1e-6 ? lo.vPort / Math.sqrt(gPrime * inp.dPort) : 99;
    const vScr = inp.Q / inp.screenA;
    const a = waveSpeed({ K: 2.2e9, rho: inp.rhoE, D: inp.D, e: inp.wall / 1000, E: inp.Ewall * 1e9 });
    const wh = joukowsky({ rho: inp.rhoE, a, dv: lo.v, L, tClose: inp.tClose });
    const surge = (wh.rapid ? wh.dP : wh.dPslow) / (inp.rhoE * G);
    // diameter that passes the design flow by gravity (outfall) or within the sump limit (intake): bisection on D
    let Dreq = inp.D;
    if (marginH < 0) {
      const mD = (D) => { const dp = liquidLinePressureDrop({ Q: inp.Q, D, rho: inp.rhoE, mu, rough: inp.rough / 1000, L, K: inp.K }); const h = (dp.dPfriction + dp.dPminor) / (inp.rhoE * G) + lo.hp; return outfall ? inp.zAvail - (inp.sea + inp.slr + h + hDens) : inp.sea + inp.slr - h - inp.zAvail; };
      let x0 = inp.D, x1 = 12; if (mD(x1) > 0) { for (let i = 0; i < 26; i++) { const mid = 0.5 * (x0 + x1); if (mD(mid) > 0) x1 = mid; else x0 = mid; } Dreq = x1; } else Dreq = 12;
    }
    const scen = slrRows(inp).map((s) => { const z = need_(s.rise), cap = capacity(s.rise); return { name: s.name, rise: s.rise, z: r(z, 3), margin: r(outfall ? inp.zAvail - z : z - inp.zAvail, 3), cap: r(cap, 3), capPct: r(cap / inp.Q * 100, 3) }; });
    const qs = linspace(0.1 * inp.Q, 1.6 * Math.max(inp.Q, Qcap), 30);
    const okH = marginH >= 0, okV = lo.v >= inp.vMin && lo.v <= inp.vMax, okScr = outfall || vScr <= inp.vScreen, okFr = !outfall || Frd >= 1;
    const dConc = Math.PI * ((Dreq + 2 * inp.wall / 1000) ** 2 - Dreq ** 2 - (inp.D + 2 * inp.wall / 1000) ** 2 + inp.D ** 2) / 4 * L;
    const trigger = scen.find((s) => s.margin < 0);

    const advice = [
      okH
        ? eng('ok', outfall ? `The design flow of ${f(inp.Q)} m³/s needs the chamber at ${f(zReq)} m (sea ${f(inp.sea + inp.slr)} m + losses ${f(hLoss)} m + density head ${f(hDens, 2)} m) against ${f(inp.zAvail)} m available: ${f(marginH, 2)} m of gravity head in hand.` : `At ${f(inp.Q)} m³/s the sump draws down to ${f(zReq)} m (sea ${f(inp.sea + inp.slr)} m − losses ${f(hLoss)} m), ${f(marginH, 2)} m above the minimum pump level.`,
          `Gravity capacity is ${f(Qcap)} m³/s (${f(Qcap / inp.Q * 100, 3)} % of design). ${trigger ? `The margin is used up at about +${f(trigger.rise)} m of sea-level rise (${trigger.name}).` : 'The margin survives every scenario listed.'} Calibrate the roughness against a measured head-loss test before relying on the margin.`,
          'EN 752:2017 §8; WRc Design Guide for Marine Treatment Schemes (1990) ch. 6', 'Gravity discharge has no operating energy; keeping the margin through the design life avoids a pumping station and its lifetime electricity.')
        : eng('fail', outfall ? `The outfall cannot pass ${f(inp.Q)} m³/s by gravity: the chamber would need ${f(zReq)} m against ${f(inp.zAvail)} m available (short by ${f(-marginH, 2)} m). Gravity capacity is ${f(Qcap)} m³/s.` : `The intake draws the sump down to ${f(zReq)} m, ${f(-marginH, 2)} m below the minimum pump level; capacity within the limit is ${f(Qcap)} m³/s.`,
          `Enlarge the conduit to about ${f(Dreq, 3)} m, or cut losses (smoother lining, fewer bends, more or larger ports), or ${outfall ? `pump the peak with ${f(pw / 1e3)} kW at ${f(pumpHead, 2)} m head, or attenuate the peak flow in storage` : 'lower the pump setting'}.`,
          'EN 752:2017 §8; WRc Design Guide for Marine Treatment Schemes (1990) ch. 6; ANSI/HI 9.8',
          `Pumping ${f(inp.hours)} h a year uses about ${f(energy, 3)} MWh; a larger conduit avoids it for roughly ${f(Math.max(0, dConc), 3)} m³ of extra concrete (≈ ${f(Math.max(0, dConc) * CO2_CONC, 3)} t CO₂e, screening factor) — compare the two on whole-life carbon.`),
      eng(okV && okFr && okScr ? 'info' : 'warn', `Conduit velocity ${f(lo.v, 2)} m/s (target ${f(inp.vMin)}–${f(inp.vMax)}), port velocity ${f(lo.vPort, 2)} m/s${outfall ? `, port densimetric Froude number ${f(Frd, 2)}` : `, through-screen velocity ${f(vScr, 2)} m/s (tolerable ${f(inp.vScreen)})`}; pump-trip surge ±${f(surge, 2)} m.`,
        `${lo.v < inp.vMin ? 'The conduit is too slow to keep itself clean: reduce the diameter, twin the pipe for staged flow, or schedule flushing. ' : ''}${lo.v > inp.vMax ? 'The conduit velocity is high: check abrasion and surge. ' : ''}${outfall && Frd < 1 ? 'Ports run with a densimetric Froude number below 1: sea water intrudes into the diffuser — close ports, fit duckbill valves or reduce port size. ' : ''}${!outfall && vScr > inp.vScreen ? `Increase the net screen area to ${f(inp.Q / inp.vScreen)} m² to protect fish and limit debris impingement. ` : ''}${okV && okFr && okScr ? 'Velocities are within the usual ranges; confirm the surge with a transient analysis if the closure is faster than assumed.' : ''}`,
        'WRc Design Guide (1990) ch. 6–7; US EPA 40 CFR 125.84 (0.15 m/s); Wylie & Streeter (1993) ch. 1'),
      eng(trigger ? 'warn' : 'info', `Sea-level rise takes head directly from a gravity ${outfall ? 'outfall' : 'intake margin'}: capacity falls from ${f(scen[0].cap)} to ${f(scen[scen.length - 1].cap)} m³/s between the first and last scenario.`,
        outfall ? 'Fix the scenario at which pumping or storage becomes necessary and reserve space and power supply for it now; fit non-return valves so that high sea levels cannot back-flood the network.' : 'A rising sea helps intake submergence but raises flood risk to the pump house: set floor and electrical levels for the end-of-life scenario.', SLR_BASIS),
      eco(okH ? 'info' : 'warn', `Head loss is ${f(hLoss)} m at design flow (friction ${f(lo.hf)} m, fittings ${f(lo.hm)} m, ports ${f(lo.hp)} m)${pw > 0 ? `; pumping the shortfall takes ${f(pw / 1e3)} kW, about ${f(energy, 3)} MWh a year` : ''}.`,
        'Choose the diameter on whole-life cost: capital rises roughly with D^1.5 while friction head falls with D^5, so a modest increase often removes decades of pumping energy and defers the sea-level-rise trigger.', 'Whole-life costing per ISO 15686-5; EN 16933-2 economic pipe sizing',
        `Each metre of avoidable pumping head at this flow and duty costs about ${f(inp.rhoE * G * inp.Q / inp.eff * inp.hours / 1e6, 3)} MWh a year of electricity.`),
    ];

    return {
      kpis: [
        { key: 'zReq', label: outfall ? 'Chamber level needed for the design flow' : 'Sump level at the design flow', value: r(zReq), unit: 'm', explain: outfall ? 'Water level the landward chamber must reach to push the design flow out against the sea: sea level plus losses plus density head.' : 'Level in the pump sump when abstracting the design flow: sea level less the losses through screens, openings and conduit.' },
        { key: 'margin', label: 'Head margin', value: r(marginH, 3), unit: 'm', limit: 0, sense: 'min', explain: 'Available level minus required level (outfall), or sump level minus the minimum pump level (intake). Negative means the flow cannot be passed without pumping or a larger conduit.' },
        { key: 'hf', label: 'Friction head loss', value: r(lo.hf), unit: 'm', explain: 'Darcy–Weisbach loss along the conduit. Grows with the square of flow and falls with the fifth power of diameter; sensitive to fouling.' },
        { key: 'hm', label: 'Minor head losses', value: r(lo.hm), unit: 'm', explain: 'Losses at entry, bends, valves and chambers, K v²/2g.' },
        { key: 'hp', label: 'Port / opening head loss', value: r(lo.hp), unit: 'm', explain: 'Head to drive the flow through the diffuser ports or intake openings with their discharge coefficient.' },
        { key: 'hDens', label: 'Density head', value: r(hDens, 3), unit: 'm', explain: 'Extra level a column of lighter effluent needs to balance sea water at the port depth. Zero for an intake.' },
        { key: 'v', label: 'Conduit velocity', value: r(lo.v, 3), unit: 'm/s', explain: 'Mean velocity in the conduit at the design flow. Too slow lets sediment and growth settle; too fast erodes and makes surge severe.' },
        { key: 'vPort', label: 'Port / opening velocity', value: r(lo.vPort, 3), unit: 'm/s', explain: 'Mean velocity through the diffuser ports or intake openings — the local flow field divers, fish and seabed sediment see.' },
        { key: 'Frd', label: 'Port densimetric Froude number', value: r(Math.min(Frd, 99), 3), unit: '-', explain: 'Port velocity over √(g′ d). Above 1 the ports run full and keep sea water out; higher values improve initial mixing.' },
        { key: 'Qcap', label: 'Capacity within the available head', value: r(Qcap), unit: 'm³/s', explain: 'Largest flow that can pass with the available chamber level (or sump limit) at the design sea level.' },
        { key: 'power', label: 'Pump shaft power for the shortfall', value: r(pw / 1e3), unit: 'kW', explain: 'Power to lift the design flow through the missing head; zero when gravity suffices.' },
        { key: 'surge', label: 'Pump-trip surge head', value: r(surge), unit: 'm', explain: 'Joukowsky pressure rise or fall when the flow is stopped in the stated time. Compare with the pressure rating and with column separation on the down-surge.' },
        { key: 'vScreen', label: 'Through-screen velocity', value: r(vScr, 3), unit: 'm/s', explain: 'Design flow over the net open screen area (intake). Above about 0.15 m/s fish are impinged and debris loading rises; informative only for an outfall.' },
      ],
      checks: [
        check(outfall ? 'Gravity head for the design flow' : 'Sump drawdown at the design flow', inp.Q, Math.max(Qcap, 1e-6), 'm³/s', 'EN 752:2017 §8; Bernoulli with Darcy–Weisbach losses', outfall ? 'Guards against the outfall backing up and flooding the upstream network when high sea level coincides with peak flow.' : 'Guards against the sump falling below the level the pumps need, which causes vortexing, air entrainment and cavitation.'),
        check('Self-cleansing velocity', inp.vMin, Math.max(lo.v, 1e-6), 'm/s', 'WRc Design Guide (1990) §6.3', 'Guards against grit, sand and marine growth settling in the conduit and progressively throttling it.'),
        check('Maximum conduit velocity', lo.v, inp.vMax, 'm/s', 'EN 752:2017 / project criterion', 'Guards against lining abrasion, cavitation at fittings and violent surge on pump trip.'),
        outfall ? check('Port purging (densimetric Froude ≥ 1)', 1, Math.max(Frd, 1e-6), '-', 'WRc Design Guide (1990) §7; Wilkinson (1988)', 'Guards against sea water entering the diffuser through low-flow ports, which traps sediment and marine growth and blocks ports.')
          : check('Through-screen velocity', vScr, inp.vScreen, 'm/s', 'US EPA 40 CFR 125.84(b)(2)', 'Guards against impingement of fish and excessive debris loading on the screens.'),
      ],
      plots: [
        { type: 'line', title: 'System curve against the available level', x: { label: 'Flow', unit: 'm³/s' }, y: { label: outfall ? 'Chamber level needed' : 'Sump level', unit: 'm' },
          series: [{ name: `Sea +${f(inp.slr)} m (design)`, x: qs.map((q) => r(q, 3)), y: qs.map((q) => r(need_(inp.slr, q), 3)) }, { name: `Sea +${f(scen[scen.length - 1].rise)} m`, x: qs.map((q) => r(q, 3)), y: qs.map((q) => r(need_(scen[scen.length - 1].rise, q), 3)), dash: true }, { name: outfall ? 'Available level' : 'Minimum pump level', x: [qs[0], qs[29]], y: [inp.zAvail, inp.zAvail], dash: true }, { name: 'Design flow', x: [inp.Q], y: [r(zReq, 3)] }],
          explain: 'Level required as a function of flow. Where the curve meets the available level is the gravity capacity; the dashed curve shows how the highest sea-level scenario shifts it.' },
        { type: 'bar', title: 'Head budget at design flow', categories: ['Friction', 'Fittings', 'Ports', 'Density head'], series: [{ name: 'm', y: [lo.hf, lo.hm, lo.hp, hDens].map((v) => r(v, 3)) }],
          explain: 'Where the head goes. A dominant friction bar points to diameter or roughness; a dominant port bar points to the diffuser design.' },
        { type: 'line', title: 'Capacity versus sea-level rise', x: { label: 'Sea-level rise', unit: 'm' }, y: { label: 'Capacity', unit: 'm³/s' },
          series: [{ name: 'Capacity', x: scen.map((s) => s.rise), y: scen.map((s) => s.cap) }, { name: 'Design flow', x: [scen[0].rise, scen[scen.length - 1].rise], y: [inp.Q, inp.Q], dash: true }],
          explain: 'Gravity capacity as the sea rises. The crossing with the design flow is the trigger for pumping, storage or a larger conduit.' },
      ],
      tables: [
        { title: 'Sea-level-rise scenarios', columns: [{ key: 'name', label: 'Scenario' }, { key: 'rise', label: 'Rise', unit: 'm' }, { key: 'z', label: outfall ? 'Level needed' : 'Sump level', unit: 'm' }, { key: 'margin', label: 'Margin', unit: 'm' }, { key: 'cap', label: 'Capacity', unit: 'm³/s' }, { key: 'capPct', label: 'Share of design flow', unit: '%' }],
          rows: scen, explain: 'Head balance and capacity for each scenario. The first row with a negative margin is the adaptation trigger.' },
        { title: 'Hydraulic parameters', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['Conduit length', L, 'm'], ['Reynolds number', lo.Re, '-'], ['Darcy friction factor', lo.fr, '-'], ['Total port area', Aport, 'm²'], ['Port area / conduit area', Aport / Apipe, '-'], ['Pressure-wave speed', a, 'm/s'], ['Pipe period 2L/a', wh.period, 's'], ['Diameter for gravity flow', Dreq, 'm']].map(([q, v, u]) => ({ q, v: r(v), u })),
          explain: 'Intermediate hydraulic values for checking against a hand calculation or a measured head-loss test. Port area should normally be 0.5–0.7 of the conduit area for even distribution.' },
      ],
      notes: [`${outfall ? 'Outfall' : 'Intake'}: total loss ${f(hLoss)} m at ${f(inp.Q)} m³/s; ${wh.rapid ? 'closure is rapid (shorter than 2L/a): full Joukowsky surge' : 'closure is slow: reduced surge'}.`, `Sea level in the design case ${f(inp.sea + inp.slr)} m including ${f(inp.slr)} m of rise.`],
      advice,
      publish: { 'marine.outfallHead': zReq, 'marine.outfallCapacity': Qcap },
      quantities: { concrete_m3: r(Math.PI * ((inp.D + 2 * inp.wall / 1000) ** 2 - inp.D ** 2) / 4 * L), pipe_km: r(L / 1000) },
    };
  },
  verify: [
    { name: 'Darcy–Weisbach with Colebrook–White friction', source: 'Hand calculation (Moody 1944; Colebrook 1939): Q = 2 m³/s, D = 1.2 m → v = 1.7684 m/s, Re = 1.929e6 (ρ = 1000, μ = 1.1e-3), ε/D = 0.00125 → f = 0.02088 (two Colebrook iterations), hf = f (L/D) v²/2g = 0.02088 × 833.3 × 0.15945 = 2.774 m; fittings K = 2: hm = 2 × 0.15945 = 0.319 m',
      inputs: { Q: 2, D: 1.2, L: 1000, rough: 1.5, K: 2, rhoE: 1000, route: null }, expect: [{ kpi: 'hf', value: 2.774, tol: 0.005 }, { kpi: 'hm', value: 0.319, tol: 0.005 }, { kpi: 'v', value: 1.7684, tol: 0.001 }] },
    { name: 'Zero-flow limit: hydrostatic balance of fresh effluent against sea water', source: 'Hydrostatics: the chamber stands (ρsea/ρe − 1) h above the sea. ρsea = 1025, ρe = 1000, h = 20 m → 0.5 m; with sea level 2 m and negligible flow the chamber level is 2.5 m',
      inputs: { Q: 0.01, hOut: 20, sea: 2, slr: 0, rhoSea: 1025, rhoE: 1000, K: 0, nPorts: 50, dPort: 0.5 }, expect: [{ kpi: 'hDens', value: 0.5, tol: 0.001 }, { kpi: 'zReq', value: 2.5, abs: 0.002 }] },
  ],
};
