// Marine towing: environmental resistance of the tow, required bollard pull, towline strength, catenary and
// dynamic tension, and the tow speed and voyage time the tug can deliver.
import { solveCatenary } from '../../engine/lines/catenary.js';
import { maruoDrift, meanDriftForce } from '../../engine/hydro/drift.js';
import { waveSpectrum } from '../../engine/hydro/floater.js';
import { G, need, run, round, fmt, kpi, check, adv, line, bar, tbl, worst, linspace, CO2_MGO, CONFIRM, N001 } from './_common.js';

/** Towline minimum breaking load over continuous bollard pull, as recalled from DNV-ST-N001 Sec. 11 (BP in tonnes). */
const mblFactor = (bp) => (bp <= 40 ? 3 : bp <= 90 ? 3.8 - bp / 50 : 2);

export default {
  id: 's15-towing', suite: 15, title: 'Towing and towline',
  summary: 'Computes the wind, current and wave-drift resistance of a towed barge, the bollard pull needed to hold it in the design weather, the towline strength, catenary sag and dynamic tension, and the achievable tow speed and voyage duration.',
  covers: [7], fidelity: 'analytical',
  standards: [`${N001} Sec. 11 (towing: holding criteria, tug efficiency, towline minimum breaking load)`, 'DNV-RP-C205 (2019) §5 (wind force), §6 (current drag), §7.4 (mean drift)', 'Irvine (1981) Cable Structures ch. 1 (elastic catenary)', 'IMO MSC/Circ.884 (guidelines for safe ocean towing)'],
  models: ['Wind force from projected area', 'Hydrodynamic drag of the hull', 'Mean wave drift force (Maruo far-field formula, reflection coefficient as input)', 'JONSWAP spectrum', 'Elastic catenary towline (cable tension relationships)', 'Towline elasticity and geometric stiffness in series', 'Rigid-body equilibrium of the tow'],
  limits: 'Steady resistance and a quasi-static towline: course stability (yaw and sway of the tow), tug propulsion characteristics, shallow-water effects on hull drag and wave-frequency towline dynamics beyond a stiffness-times-motion estimate are not modelled. The holding sea state, tug efficiency and towline breaking-load factors are recalled DNV-ST-N001 values exposed as inputs to be confirmed; the drift force uses a constant reflection coefficient rather than a diffraction analysis.',
  publishes: ['install.towBollardPull', 'install.towlineTension', 'install.towDays'],
  searchTerms: 'ocean towing bollard pull towline catenary tension barge tow resistance',
  inputs: [
    { key: 'Lb', label: 'Towed barge length', unit: 'm', type: 'number', default: 120, min: 20, max: 400, group: 'Tow', role: 'geometry' },
    { key: 'Bb', label: 'Barge beam', unit: 'm', type: 'number', default: 36, min: 6, max: 120, group: 'Tow', role: 'geometry' },
    { key: 'Tb', label: 'Barge draft', unit: 'm', type: 'number', default: 4.5, min: 0.5, max: 25, group: 'Tow', role: 'ic', bind: 'install.bargeDraft', help: 'Initial condition from the loading condition (stability module).' },
    { key: 'CdHull', label: 'Hull drag coefficient on B × T', unit: '-', type: 'number', default: 1.0, min: 0.2, max: 2.5, group: 'Tow', role: 'load', uncertain: { dist: 'normal', cov: 0.2 }, help: 'About 0.8–1.2 for a raked box barge including skegs; higher for a square bow.' },
    { key: 'Awind', label: 'Head-on windage area (hull above water + cargo)', unit: 'm²', type: 'number', default: 1100, min: 0, max: 50000, group: 'Tow', role: 'geometry' },
    { key: 'Cs', label: 'Wind shape coefficient', unit: '-', type: 'number', default: 1.1, min: 0.3, max: 2, group: 'Tow', role: 'load' },
    { key: 'Rrefl', label: 'Wave reflection coefficient of the hull', unit: '-', type: 'number', default: 0.6, min: 0, max: 1, group: 'Tow', role: 'load', uncertain: { dist: 'uniform', min: 0.4, max: 0.9 }, help: 'Mean drift per unit width is ½·ρ·g·a²·R² in deep water. 1.0 is the upper bound for a wall-sided hull in short waves; shallow raked bows reflect less.' },
    { key: 'Vw', label: 'Holding condition: wind speed', unit: 'm/s', type: 'number', default: 20, min: 0, max: 60, group: 'Holding condition', role: 'bc', help: `Zero forward speed against this wind, sea and current acting together. Unrestricted tows: 20 m/s, Hs 5 m, 0.5 m/s. ${CONFIRM}` },
    { key: 'Hs', label: 'Holding condition: significant wave height', unit: 'm', type: 'number', default: 5, min: 0, max: 15, group: 'Holding condition', role: 'bc', help: CONFIRM },
    { key: 'Tp', label: 'Holding condition: peak period', unit: 's', type: 'number', default: 10, min: 3, max: 22, group: 'Holding condition', role: 'bc' },
    { key: 'Uc', label: 'Holding condition: current', unit: 'm/s', type: 'number', default: 0.5, min: 0, max: 4, group: 'Holding condition', role: 'bc', help: CONFIRM },
    { key: 'depth', label: 'Least water depth on the route', unit: 'm', type: 'number', default: 60, min: 5, max: 6000, group: 'Route', role: 'environment', help: 'Governs towline clearance to the seabed.' },
    { key: 'BP', label: 'Tug continuous bollard pull', unit: 't', type: 'number', default: 150, min: 5, max: 600, group: 'Tug', role: 'bc', design: true },
    { key: 'Te', label: 'Tug efficiency in the holding sea state', unit: '-', type: 'number', default: 0.75, min: 0.2, max: 1, group: 'Tug', role: 'criteria', help: `Fraction of bollard pull available as towline pull in waves; about 0.75 for a large tug in Hs 5 m, much less for small tugs. ${CONFIRM}` },
    { key: 'vFree', label: 'Tug free-running speed', unit: 'kn', type: 'number', default: 14, min: 6, max: 22, group: 'Tug', role: 'operation', help: 'Available pull is taken to fall linearly from the bollard pull at rest to zero at this speed.' },
    { key: 'mat', label: 'Towline construction', type: 'select', source: 'materials:wire', default: 'wire-sixstrand', group: 'Towline', role: 'material' },
    { key: 'd', label: 'Towline diameter', unit: 'mm', type: 'number', default: 76, min: 20, max: 160, group: 'Towline', role: 'geometry', design: true },
    { key: 'Lline', label: 'Towline length deployed', unit: 'm', type: 'number', default: 900, min: 100, max: 2500, group: 'Towline', role: 'bc', design: true, help: 'Longer lines give a softer catenary and lower dynamic tension but more sag.' },
    { key: 'fMbl', label: 'Required MBL / bollard pull (0 = recalled rule)', unit: '-', type: 'number', default: 0, min: 0, max: 6, group: 'Criteria', role: 'criteria', help: `Recalled rule: 3.0 up to 40 t, 3.8 − BP/50 for 40–90 t, 2.0 above 90 t. ${CONFIRM}` },
    { key: 'xDyn', label: 'Relative surge amplitude tug – tow at wave frequency', unit: 'm', type: 'number', default: 1.0, min: 0, max: 15, group: 'Towline', role: 'bc', help: 'Single amplitude of the relative motion along the towline in the design sea state; multiplied by the line stiffness for the dynamic tension.' },
    { key: 'clearReq', label: 'Required towline clearance above the seabed', unit: 'm', type: 'number', default: 5, min: 0, max: 100, group: 'Criteria', role: 'criteria' },
    { key: 'dist', label: 'Tow distance', unit: 'nm', type: 'number', default: 900, min: 1, max: 15000, group: 'Route', role: 'operation' },
    { key: 'fuel', label: 'Tug fuel consumption at tow power', unit: 't/day', type: 'number', default: 28, min: 1, max: 150, group: 'Tug', role: 'operation' },
    { key: 'vTransit', label: 'Wind and sea during normal transit (fraction of holding condition)', unit: '-', type: 'number', default: 0.4, min: 0, max: 1, group: 'Route', role: 'environment', help: 'Scales wind speed and wave height of the holding condition for the average tow-speed estimate.' },
  ],
  compute(inp, ctx) {
    // t → N (bollard pull as force), mm → m, nm → m.
    const wire = ctx.material(inp.mat); need(wire, `Towline '${inp.mat}' is not in the materials library.`);
    const rho = 1025, d = inp.d / 1e3;
    const MBL = wire.mblCoeff * d * d, mLin = wire.massCoeff * d * d, EA = wire.eaCoeff * d * d, w = mLin * G * 0.87; // steel wire: submerged ≈ 0.87 × dry weight
    const Fwind = (V) => 0.5 * 1.225 * inp.Cs * inp.Awind * V * V;
    const Fcur = (U) => 0.5 * rho * inp.CdHull * inp.Bb * inp.Tb * U * U;
    const om = linspace(0.2, 2.5, 120);
    const Fdrift = (Hs) => {
      if (Hs <= 0) return 0;
      const S = waveSpectrum(om, { Hs, Tp: inp.Tp, gamma: 3.3 });
      const D = om.map((o) => maruoDrift({ a: 1, R: inp.Rrefl, T: 2 * Math.PI / o, depth: Math.max(inp.depth, 5), rho }) * inp.Bb);
      return meanDriftForce({ omega: om, S, D });
    };
    const Fw = Fwind(inp.Vw), Fc = Fcur(inp.Uc), Fd = run(() => Fdrift(inp.Hs), 'Wave drift force');
    const TPR = Fw + Fc + Fd; // towline pull required
    const BPreq = TPR / inp.Te / (G * 1e3);
    const pull = inp.BP * 1e3 * G * inp.Te; // towline pull available in the holding sea state
    const fM = inp.fMbl > 0 ? inp.fMbl : mblFactor(inp.BP), MBLreq = fM * inp.BP * 1e3 * G;
    // Static towline at full bollard pull: symmetric elastic catenary between two points at the same level.
    const H = inp.BP * 1e3 * G, a = H / w, half = inp.Lline / 2;
    const sag = Math.sqrt(a * a + half * half) - a, span = 2 * a * Math.asinh(half / a) + H * inp.Lline / EA;
    const Tend = Math.hypot(H, w * half);
    let kLine = 1 / (inp.Lline / EA + (w * w * inp.Lline ** 3) / (12 * H ** 3)), Tcat = Tend;
    try {
      const sol = solveCatenary({ X: span, Z: 0, L: inp.Lline, w, EA, seabed: false });
      if (Number.isFinite(sol.stiffness?.[0]?.[0]) && sol.stiffness[0][0] > 0) kLine = sol.stiffness[0][0];
      if (Number.isFinite(sol.T)) Tcat = sol.T;
    } catch { /* closed-form values stand */ }
    const dT = kLine * inp.xDyn, Tmax = Tcat + dT, Tmin = Tcat - dT;
    // Sag at reduced pull (e.g. slowing down): the line drops as tension falls.
    const pulls = linspace(0.1, 1, 19).map((f) => f * H), sags = pulls.map((h) => { const aa = h / w; return Math.sqrt(aa * aa + half * half) - aa; });
    // Tow speed: the pull a tug can deliver falls about linearly from bollard pull at rest to zero at its free-running
    // speed; the tow speed is where that pull, less wind and wave drift, equals the hull drag (bisection).
    const vFree = inp.vFree * 0.514444;
    const speedFor = (frac) => {
      const env = Fwind(inp.Vw * frac) + run(() => Fdrift(inp.Hs * frac), 'Wave drift force');
      const eff = frac > 0 ? Math.min(1, inp.Te + (1 - inp.Te) * (1 - frac)) : 1;
      const net = (v) => inp.BP * 1e3 * G * eff * (1 - v / vFree) - env - 0.5 * rho * inp.CdHull * inp.Bb * inp.Tb * v * v;
      if (net(0) <= 0) return 0;
      let lo = 0, hi = vFree;
      for (let i = 0; i < 60; i++) { const mid = 0.5 * (lo + hi); if (net(mid) > 0) lo = mid; else hi = mid; }
      return 0.5 * (lo + hi);
    };
    const vCalm = speedFor(0), vTow = speedFor(inp.vTransit);
    need(vTow > 0.05, 'The tug cannot make headway in the transit weather entered: increase the bollard pull or lower the transit weather fraction.');
    const days = inp.dist * 1852 / vTow / 86400, fuelT = days * inp.fuel;
    const checks = [
      check('Bollard pull to hold the tow', BPreq, inp.BP, 't', `${N001} Sec. 11 (towing: bollard pull requirement, tug efficiency)`, 'The tug must hold the tow at zero speed against the design wind, sea and current together; otherwise the tow drifts astern in a storm.'),
      check('Towline minimum breaking load', MBLreq / 1e3, MBL / 1e3, 'kN', `${N001} Sec. 11 (towline MBL as a multiple of bollard pull)`, 'The towline must be stronger than the tug by the stated factor so that it survives peak loads in a seaway.'),
      check('Peak towline tension (static + dynamic)', Tmax / 1e3, MBL / 2 / 1e3, 'kN', 'Working-load practice: peak tension not above 50 % of MBL in the design sea state (confirm with the towing manual)', 'Mean tension at full pull plus the wave-frequency variation from the relative motion of tug and tow. High peaks fatigue the line at the fairlead and can part it.'),
      check('Towline clearance above the seabed', inp.clearReq + sag, inp.depth, 'm (sag + clearance against depth)', `${N001} Sec. 11 (towline catenary in shallow water)`, 'A long heavy towline at full pull sags below the tug stern. In shallow water it drags on the seabed, wears and can snag.'),
    ];
    if (Tmin <= 0) checks.push(check('Towline goes slack in the design sea', dT / 1e3, Tcat / 1e3, 'kN', 'Slack–snatch avoidance (DNV-RP-N103 §4.7 principle applied to a towline)', 'Dynamic variation larger than the mean tension unloads the line; the following snatch load can be several times the mean.'));
    const gov = worst(checks), ok = gov.uc <= 1;
    const dReq = Math.sqrt(MBLreq / wire.mblCoeff) * 1e3;
    return {
      kpis: [
        kpi('Fwind', 'Wind force on the tow', Fw / 1e3, 'kN', '½·ρair·Cs·A·V² on the head-on windage area in the holding wind.'),
        kpi('Fcur', 'Current drag on the hull', Fc / 1e3, 'kN', '½·ρ·Cd·B·T·U² on the immersed frontal area in the holding current.'),
        kpi('Fdrift', 'Mean wave drift force', Fd / 1e3, 'kN', 'Second-order mean force from wave reflection in the holding sea state; it usually dominates for a beamy barge.'),
        kpi('TPR', 'Towline pull required', TPR / 1e3, 'kN', 'Sum of wind, current and wave drift: the pull the towline must deliver to hold the tow.'),
        kpi('BPreq', 'Required continuous bollard pull', BPreq, 't', 'Towline pull required divided by tug efficiency in the sea state.', inp.BP, 'max'),
        kpi('MBL', 'Towline minimum breaking load', MBL / 1e3, 'kN', 'From the construction coefficient and diameter; replace with the certificate value.', MBLreq / 1e3, 'min'),
        kpi('dReq', 'Required towline diameter', dReq, 'mm', 'Diameter that just meets the MBL requirement for this construction.'),
        kpi('Tstat', 'Static towline tension at the ends (full pull)', Tcat / 1e3, 'kN', 'Tension at tug and tow with the full bollard pull as horizontal component: √(H² + (w·L/2)²).'),
        kpi('sag', 'Towline sag at full bollard pull', sag, 'm', 'Depth of the lowest point of the catenary below the towing points: √(a² + (L/2)²) − a with a = H/w.', inp.depth - inp.clearReq, 'max'),
        kpi('kLine', 'Towline stiffness (elastic + catenary)', kLine / 1e3, 'kN/m', 'Change of horizontal tension per metre of relative surge at full pull. Stiff short lines transmit wave motion as tension.'),
        kpi('Tmax', 'Peak towline tension in the design sea', Tmax / 1e3, 'kN', 'Static end tension plus stiffness × relative surge amplitude.', MBL / 2 / 1e3, 'max'),
        kpi('vCalm', 'Tow speed in calm water', vCalm / 0.514444, 'kn', 'Speed at which hull drag equals the pull the tug can still deliver at that speed, with no wind or waves.'),
        kpi('vTow', 'Average tow speed in transit weather', vTow / 0.514444, 'kn', 'Speed with wind and wave drift of the transit weather deducted from the available pull.'),
        kpi('days', 'Tow duration', days, 'days', 'Distance over average tow speed; add weather routing, bunkering and port time.'),
        kpi('fuelT', 'Tug fuel for the tow', fuelT, 't', 'Duration times the daily consumption entered.'),
      ],
      checks,
      plots: [
        bar('Resistance of the tow in the holding condition', ['Wind', 'Current', 'Wave drift', 'Total', 'Available pull'], 'Force', 'kN', [{ name: 'Force', y: [Fw, Fc, Fd, TPR, pull].map((v) => round(v / 1e3)) }],
          'The three environmental components against the towline pull the tug can deliver in the sea state. The total must not exceed the available pull.'),
        line('Towline sag against pull', 'Horizontal towline pull', 'kN', 'Sag below the towing points', 'm', [{ name: 'Sag', x: pulls.map((v) => round(v / 1e3)), y: sags.map((v) => round(v)) }, { name: 'Depth − clearance', x: [round(pulls[0] / 1e3), round(H / 1e3)], y: [inp.depth - inp.clearReq, inp.depth - inp.clearReq], dash: true }],
          'The line hangs deeper as the tug eases off. Where the curve crosses the dashed line is the least pull at which this length can be used in the shallowest part of the route.'),
      ],
      tables: [tbl('Towline properties and loads', [['k', 'Quantity'], ['v', 'Value'], ['u', 'Unit']], [
        ['Mass per metre', mLin, 'kg/m'], ['Submerged weight', w, 'N/m'], ['Axial stiffness EA', EA / 1e6, 'MN'], ['MBL', MBL / 1e3, 'kN'], ['Required MBL factor', fM, '× BP'], ['Horizontal span at full pull', span, 'm'], ['Dynamic tension amplitude', dT / 1e3, 'kN'], ['Minimum tension', Tmin / 1e3, 'kN'],
      ].map(([k, v, u]) => ({ k, v: round(v), u })), 'Derived towline data and the tensions at full bollard pull; compare with the towline certificate and the towing winch settings.')],
      notes: [`Holding condition ${inp.Vw} m/s, Hs ${inp.Hs} m, ${inp.Uc} m/s: wind ${fmt(Fw / 1e3)} kN, current ${fmt(Fc / 1e3)} kN, drift ${fmt(Fd / 1e3)} kN.`, `Towline ${inp.d} mm × ${inp.Lline} m: MBL ${fmt(MBL / 1e3)} kN against ${fmt(MBLreq / 1e3)} kN required (factor ${fmt(fM)}).`, `"${gov.name}" governs at ${fmt(gov.uc)}.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `A ${inp.BP} t tug holds the tow (required ${fmt(BPreq)} t) and the ${inp.d} mm towline is adequate; peak tension ${fmt(Tmax / 1e3)} kN, sag ${fmt(sag)} m in ${inp.depth} m of water.`,
          `Shorten the towline in water shallower than ${fmt(sag + inp.clearReq)} m, keep the catenary by not exceeding ${fmt(100 * Tmax / MBL)} % of MBL in peaks, and confirm the tug efficiency and holding sea state against the governing edition and the warranty surveyor's requirements.`, `${N001} Sec. 11; IMO MSC/Circ.884`)
          : adv('engineering', 'fail', `"${gov.name}" is not met (${fmt(gov.uc)}): required bollard pull ${fmt(BPreq)} t, required MBL ${fmt(MBLreq / 1e3)} kN, peak tension ${fmt(Tmax / 1e3)} kN, sag ${fmt(sag)} m.`,
            gov.name.startsWith('Bollard') ? `Use a tug of at least ${fmt(Math.ceil(BPreq / 5) * 5)} t continuous bollard pull, add a second tug, or restrict the tow to a weather window with a lower holding sea state.` : gov.name.startsWith('Towline min') ? `Increase the towline to at least ${fmt(Math.ceil(dReq / 2) * 2)} mm.` : gov.name.startsWith('Towline clear') ? `Shorten the line to about ${fmt(inp.Lline * Math.sqrt(Math.max(inp.depth - inp.clearReq, 1) / sag))} m in the shallow section or keep more tension on it.` : 'Lengthen the towline or add a stretcher / chain pendant to soften it, and reduce power in heavy weather.',
            `${N001} Sec. 11; IMO MSC/Circ.884`),
        adv('economic', 'info', `The tow takes about ${fmt(days)} days at ${fmt(vTow / 0.514444)} kn and burns about ${fmt(fuelT)} t of fuel; the tug is used to ${fmt(100 * BPreq / inp.BP)} % of its bollard pull in the holding condition.`,
          BPreq < 0.6 * inp.BP ? `The tug is larger than the holding criterion needs: a tug of about ${fmt(Math.ceil(BPreq * 1.15 / 5) * 5)} t would hold the tow at a lower day rate, at the price of a slower passage — compare day rate × days for both.` : 'Wave drift sets the bollard pull: a seasonal or weather-routed tow with a lower holding sea state allows a smaller tug class, traded against waiting on weather.',
          'Tug day rate × tow duration plus fuel; weather-routed versus unrestricted tow per DNV-ST-N001 Sec. 2; AACE RP 18R-97 for the estimate class',
          `Fuel for the tow ≈ ${fmt(fuelT)} t marine gas oil ≈ ${fmt(fuelT * CO2_MGO)} t CO₂ (${CO2_MGO} t CO₂ per t fuel, screening factor); a 10 % lower speed cuts the hull-drag fuel per mile by about 19 %.`),
      ],
      publish: { 'install.towBollardPull': BPreq * 1e3 * G, 'install.towlineTension': Tmax, 'install.towDays': days },
      quantities: { vessel_days: { tug: round(days) } },
    };
  },
  verify: [
    { name: 'Wind and current resistance: F = ½·ρ·C·A·V²', source: 'DNV-RP-C205 §5.3 and §6.2, hand calculation: wind ½·1.225·1.1·1100·20² = 296.45 kN; current ½·1025·1.0·36·4.5·0.5² = 20.76 kN',
      inputs: {}, expect: [{ kpi: 'Fwind', value: 0.5 * 1.225 * 1.1 * 1100 * 400 / 1e3, tol: 1e-9 }, { kpi: 'Fcur', value: 0.5 * 1025 * 36 * 4.5 * 0.25 / 1e3, tol: 1e-9 }] },
    { name: 'Deep-water mean drift on a reflecting wall: F = ρ·g·Hs²·R²·B/16', source: 'Maruo (1960) / Longuet-Higgins (1977): ½·ρ·g·a²·R² per unit width with the mean square amplitude 2·m0 = Hs²/8; R = 1, B = 36 m, Hs = 5 m, deep water → 1025·9.80665·25·36/16 = 565.4 kN (spectrum truncated at 0.2–2.5 rad/s, hence 2 % tolerance)',
      inputs: { Rrefl: 1, depth: 3000, Tp: 8 }, expect: [{ kpi: 'Fdrift', value: 1025 * 9.80665 * 25 * 36 / 16 / 1e3, tol: 0.02 }] },
    { name: 'Catenary sag: √(a² + (L/2)²) − a, a = H/w', source: 'Inextensible catenary (Irvine 1981 §1.2): 76 mm six-strand wire, m = 4000·0.076² = 23.104 kg/m, w = 0.87·23.104·9.80665 = 197.12 N/m; H = 150 t·g = 1471.0 kN → a = 7462.5 m; L/2 = 450 m → sag = 13.556 m',
      inputs: {}, expect: [{ kpi: 'sag', value: Math.sqrt((150e3 * 9.80665 / (0.87 * 4000 * 0.076 ** 2 * 9.80665)) ** 2 + 450 ** 2) - 150e3 * 9.80665 / (0.87 * 4000 * 0.076 ** 2 * 9.80665), tol: 1e-6 }] },
  ],
};
