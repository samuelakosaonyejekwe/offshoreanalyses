// Skirt penetration, suction installation and uplift / breakout resistance of mudmat skirts and suction buckets.
import { skirtPenetrationClay, skirtPenetrationCpt, selfWeightPenetration, breakoutForce } from '../../engine/geotech/mudmat.js';
import { suctionInstallClay, suctionPenetrationCurve, caissonVerticalCapacityClay, ringBucklingPressure } from '../../engine/geotech/suction.js';
import { G, need, run, fmt, round, adv, check, steelCo2Text, concreteCo2Text, vesselCo2Text, CONFIRM, NOT_FE } from './_common.js';

export default {
  id: 's04-skirt-suction-uplift', suite: 4, title: 'Suction, penetration and uplift resistance',
  summary: 'Skirt penetration resistance and self-weight penetration, the suction or ballast needed to reach full penetration with its plug-heave limit, and uplift and breakout resistance for retrieval.',
  covers: [9], fidelity: 'analytical',
  standards: ['DNV-RP-C212 (2019) §6.2 (skirt penetration)', 'DNV-RP-E303 (2021) suction anchors in clay', 'DNV-RP-N103 (2017) §6.2–6.3 (landing, retrieval)', 'ISO 19901-4:2016 §7.6', 'Houlsby & Byrne (2005) Proc. ICE Geotech. Eng. 158'],
  models: ['Skirt penetration resistance (wall friction and tip bearing)', 'CPT-based penetration resistance (kp, kf)', 'Suction-assisted penetration and plug-heave limit (Houlsby & Byrne)', 'Uplift resistance (reverse end bearing, wall friction)', 'Suction breakout factor model', 'Hydrostatic pressure and cavitation limit', 'Ring buckling under suction (Bresse–Bryan)'],
  limits: `Clay resistance from su = su0 + k·z with an adhesion factor; sand through CPT coefficients only. Breakout suction is a simple factor model (a user fraction of reverse end bearing), not a rate-dependent consolidation analysis. Seepage-assisted penetration in sand, set-up with time and soil plug stability beyond the reverse-bearing check are not modelled. ${NOT_FE}`,
  publishes: ['subsea.skirtPenetrationForce', 'subsea.breakoutForce', 'subsea.upliftResistance'],
  searchTerms: 'mudmat skirt penetration resistance suction bucket installation clay breakout retrieval',
  inputs: [
    { key: 'type', label: 'Foundation type', type: 'select', default: 'mat', group: 'Foundation', role: 'geometry', options: [{ value: 'mat', label: 'Skirted mudmat (rectangular, compartments)' }, { value: 'bucket', label: 'Suction bucket (circular)' }] },
    { key: 'B', label: 'Mudmat width', unit: 'm', type: 'number', default: 10, min: 0.5, max: 60, group: 'Foundation', role: 'geometry', bind: 'subsea.matB', showIf: { key: 'type', in: ['mat'] } },
    { key: 'L', label: 'Mudmat length', unit: 'm', type: 'number', default: 14, min: 0.5, max: 80, group: 'Foundation', role: 'geometry', bind: 'subsea.matL', showIf: { key: 'type', in: ['mat'] } },
    { key: 'ribs', label: 'Internal skirt (rib) length', unit: 'm', type: 'number', default: 28, min: 0, max: 2000, group: 'Foundation', role: 'geometry', showIf: { key: 'type', in: ['mat'] }, help: 'Total plan length of internal skirts, added to the perimeter.' },
    { key: 'Dc', label: 'Bucket diameter', unit: 'm', type: 'number', default: 5, min: 0.5, max: 40, group: 'Foundation', role: 'geometry', design: true, showIf: { key: 'type', in: ['bucket'] } },
    { key: 'nB', label: 'Number of buckets', type: 'int', default: 4, min: 1, max: 12, group: 'Foundation', role: 'geometry', showIf: { key: 'type', in: ['bucket'] } },
    { key: 'hs', label: 'Skirt height', unit: 'm', type: 'number', default: 0.6, min: 0.05, max: 30, group: 'Foundation', role: 'geometry', design: true },
    { key: 'ts', label: 'Skirt thickness', unit: 'mm', type: 'number', default: 15, min: 4, max: 80, group: 'Foundation', role: 'geometry', design: true },
    { key: 'mat', label: 'Skirt steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Foundation', role: 'material' },
    { key: 'W', label: 'Submerged weight available at landing', unit: 'kN', type: 'number', default: 1000, min: 1, max: 2e5, group: 'Loads', role: 'ic', bind: 'subsea.submergedWeight', uncertain: { dist: 'normal', cov: 0.03 }, help: 'Weight on the skirts once the hook load is released.' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 150, min: 5, max: 3500, group: 'Site', role: 'bc', bind: 'site.depth', help: 'Sets the external hydrostatic pressure and thus the cavitation limit on suction.' },
    { key: 'soil', label: 'Soil model', type: 'select', default: 'clay', group: 'Soil', role: 'soil', options: [{ value: 'clay', label: 'Clay — su profile' }, { value: 'cpt', label: 'CPT-based (sand or clay)' }] },
    { key: 'su0', label: 'Undrained shear strength at mudline', unit: 'kPa', type: 'number', default: 5, min: 0.2, max: 500, group: 'Soil', role: 'soil', bind: 'site.su0', uncertain: { dist: 'lognormal', cov: 0.25 }, showIf: { key: 'soil', in: ['clay'] } },
    { key: 'suGrad', label: 'Strength gradient', unit: 'kPa/m', type: 'number', default: 1.5, min: 0, max: 50, group: 'Soil', role: 'soil', showIf: { key: 'soil', in: ['clay'] } },
    { key: 'St', label: 'Sensitivity (intact / remoulded)', unit: '-', type: 'number', default: 3, min: 1, max: 20, group: 'Soil', role: 'soil', uncertain: { dist: 'lognormal', cov: 0.2 }, showIf: { key: 'soil', in: ['clay'] }, help: 'Penetration adhesion factor α = 1/sensitivity: the wall shears remoulded soil.' },
    { key: 'NcTip', label: 'Skirt tip bearing factor', unit: '-', type: 'number', default: 7.5, min: 5, max: 12, group: 'Soil', role: 'soil', showIf: { key: 'soil', in: ['clay'] } },
    { key: 'qc0', label: 'Cone resistance at mudline', unit: 'MPa', type: 'number', default: 0.5, min: 0.01, max: 60, group: 'Soil', role: 'soil', showIf: { key: 'soil', in: ['cpt'] } },
    { key: 'qcGrad', label: 'Cone resistance gradient', unit: 'MPa/m', type: 'number', default: 2, min: 0, max: 30, group: 'Soil', role: 'soil', showIf: { key: 'soil', in: ['cpt'] } },
    { key: 'kp', label: 'CPT tip coefficient kp', unit: '-', type: 'number', default: 0.6, min: 0.05, max: 1.5, group: 'Soil', role: 'soil', showIf: { key: 'soil', in: ['cpt'] }, help: `DNV: clay 0.4 / 0.6, sand 0.3 / 0.6 (most probable / highest expected). ${CONFIRM}` },
    { key: 'kf', label: 'CPT friction coefficient kf', unit: '-', type: 'number', default: 0.003, min: 0.0005, max: 0.1, group: 'Soil', role: 'soil', showIf: { key: 'soil', in: ['cpt'] }, help: `DNV: clay 0.03 / 0.05, sand 0.001 / 0.003. ${CONFIRM}` },
    { key: 'gammaSub', label: 'Submerged unit weight γ′', unit: 'kN/m³', type: 'number', default: 6, min: 2, max: 13, group: 'Soil', role: 'soil', bind: 'site.gammaSub' },
    { key: 'setup', label: 'Long-term adhesion factor (set-up)', unit: '-', type: 'number', default: 0.65, min: 0.1, max: 1, group: 'Uplift and retrieval', role: 'soil', help: 'Wall adhesion after consolidation, used for in-place uplift and retrieval.' },
    { key: 'beta', label: 'Breakout suction factor β', unit: '-', type: 'number', default: 0.5, min: 0, max: 1, group: 'Uplift and retrieval', role: 'criteria',
      help: 'Fraction of reverse end bearing mobilised as base suction during a lift: 0 vented or very slow, up to 1 for a rapid lift of a sealed base. Simple factor model — DNV-RP-N103 §6.3 allows suction up to the bearing capacity. Confirm against the governing edition.' },
    { key: 'Uup', label: 'Design uplift load in place', unit: 'kN', type: 'number', default: 200, min: 0, max: 2e5, group: 'Uplift and retrieval', role: 'bc', help: 'Net upward load on the foundation from overturning, lift or accidental snag.' },
    { key: 'fsUp', label: 'Required uplift factor of safety', unit: '-', type: 'number', default: 1.5, min: 1, max: 4, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'fsPlug', label: 'Required factor on plug heave', unit: '-', type: 'number', default: 1.5, min: 1, max: 4, group: 'Acceptance criteria', role: 'criteria', help: `Allowable over required suction. ${CONFIRM}` },
    { key: 'craneCap', label: 'Crane capacity available for retrieval', unit: 'kN', type: 'number', default: 4000, min: 10, max: 1e6, group: 'Acceptance criteria', role: 'criteria' },
  ],
  compute(inp, ctx) {
    // SI: kN → N, kPa/MPa → Pa, mm → m.
    const steel = ctx.material(inp.mat);
    need(steel, `Steel '${inp.mat}' is not in the materials library.`);
    const mat = inp.type === 'mat', clay = inp.soil === 'clay';
    need(mat || clay, 'The suction-bucket branch uses the clay strength profile: choose the clay soil model, or analyse the skirts as a mudmat with CPT coefficients.');
    const t = inp.ts / 1e3, hs = inp.hs, W = inp.W * 1e3, gamma = inp.gammaSub * 1e3;
    const su0 = inp.su0 * 1e3, k = inp.suGrad * 1e3, alpha = 1 / inp.St;
    const A = mat ? inp.B * inp.L : inp.nB * Math.PI * inp.Dc * inp.Dc / 4;
    const wall = mat ? 2 * (inp.B + inp.L) + inp.ribs : inp.nB * Math.PI * (inp.Dc - t);
    need(mat || t < inp.Dc / 2, 'Skirt thickness must be smaller than the bucket radius.');
    const qc = (z) => (inp.qc0 + inp.qcGrad * z) * 1e6;
    const Rz = (z) => (clay ? skirtPenetrationClay({ z, wallLength: wall, t, su0, k, alpha, Nc: inp.NcTip, gamma }) : skirtPenetrationCpt({ z, qc, wallLength: wall, t, kp: inp.kp, kf: inp.kf }));
    const full = run(() => Rz(hs), 'Skirt penetration resistance');
    // Self-weight penetration (closed solver for clay, bisection for CPT).
    let zSelf;
    if (clay) zSelf = run(() => selfWeightPenetration({ W, skirtHeight: hs, wallLength: wall, t, su0, k, alpha, Nc: inp.NcTip, gamma }), 'Self-weight penetration').z;
    else if (full.R <= W) zSelf = hs; else if (Rz(0).R >= W) zSelf = 0; else { let lo = 0, hi = hs; for (let i = 0; i < 50; i++) { const m = 0.5 * (lo + hi); if (Rz(m).R < W) lo = m; else hi = m; } zSelf = 0.5 * (lo + hi); }
    const deficit = Math.max(full.R - W, 0);
    // Suction (under-pressure in the skirt compartments) that supplies the deficit, and its limits.
    const sReq = deficit / A;
    const suTip = clay ? su0 + k * hs : qc(hs) / 15; // CPT: su ≈ qc/Nkt with Nkt ≈ 15 for the plug-heave estimate only
    const suAvg = clay ? su0 + 0.5 * k * hs : suTip;
    let sAllow;
    if (!mat) {
      const r = run(() => suctionInstallClay({ Do: inp.Dc, t, h: hs, su: { su0, k }, alphaO: alpha, alphaI: alpha, gamma, Vprime: W / inp.nB }), 'Suction installation');
      sAllow = r.sAllow;
    } else {
      // Reverse end bearing of the plug in one compartment plus inside wall friction (Houlsby & Byrne 2005a form).
      sAllow = 6.2 * suTip + alpha * suAvg * hs * wall / A;
    }
    const pHyd = 1025 * G * inp.depth + 101325;
    const sCav = pHyd - 2300;
    const sLim = Math.min(sAllow, sCav);
    const pCr = !mat ? ringBucklingPressure({ E: steel.E, nu: steel.nu ?? 0.3, t, D: inp.Dc - t }) : NaN;
    // Uplift in place and breakout on retrieval, with set-up adhesion.
    const fr = inp.setup * suAvg * hs * wall * 2; // both faces of every skirt
    const cap = !mat ? run(() => caissonVerticalCapacityClay({ Do: inp.Dc, t, L: hs, su: { su0, k }, alphaO: inp.setup, alphaI: inp.setup, W: W / inp.nB, gamma }), 'Caisson pull-out capacity') : null;
    const Nrev = mat ? 1.2 * (2 + Math.PI) : 9;
    const upSealed = cap ? cap.sealed * inp.nB : W + 0.5 * fr + Nrev * suTip * A;
    const upVented = cap ? cap.vented * inp.nB : W + fr;
    const bo = breakoutForce({ W, A, su: suTip, breakoutFactor: inp.beta, Nc: Nrev, skirtFriction: mat ? fr : 0.5 * fr });
    const boVented = breakoutForce({ W, A, su: suTip, breakoutFactor: 0, Nc: Nrev, skirtFriction: mat ? fr : 0.5 * fr });
    const zs = [], Rs = [], Fr = [], Tp = [], sR = [], sA = [];
    for (let i = 0; i <= 40; i++) { const z = hs * i / 40; const r = Rz(z); zs.push(round(-z)); Rs.push(r.R / 1e3); Fr.push(r.friction / 1e3); Tp.push(r.tip / 1e3); sR.push(Math.max(r.R - W, 0) / A / 1e3); sA.push((mat ? 6.2 * (clay ? su0 + k * z : qc(z) / 15) + alpha * (clay ? su0 + 0.5 * k * z : qc(z) / 15) * z * wall / A : sAllow * (z / hs || 0)) / 1e3); }
    let curve = null;
    if (!mat) try { curve = suctionPenetrationCurve({ soil: 'clay', L: hs, n: 40, params: { Do: inp.Dc, t, su: { su0, k }, alphaO: alpha, alphaI: alpha, gamma, Vprime: W / inp.nB } }); } catch { curve = null; }
    const checks = [
      check('Penetration: resistance against weight plus allowable suction', full.R / 1e3, (W + sLim / inp.fsPlug * A) / 1e3, 'kN', 'DNV-RP-C212 (2019) §6.2; Houlsby & Byrne (2005a)', 'Guards against refusal before the skirts are fully home. If exceeded the mat stands proud, the base is not in contact and capacity is not achieved.'),
      check('Suction against plug heave and cavitation', sReq / 1e3, sLim / inp.fsPlug / 1e3, 'kPa', 'Houlsby & Byrne (2005a) reverse-bearing limit; DNV-RP-E303 §5', 'Guards against the soil plug being sucked up inside the skirt (or the water cavitating) instead of the skirt going down.'),
      check('Uplift in place (vented base)', inp.Uup * inp.fsUp, upVented / 1e3, 'kN', 'ISO 19901-4:2016 §7.6; DNV-RP-E303 §6 (vented, wall friction only)', 'Guards against the foundation being lifted by sustained upward load, when no base suction can be relied on.'),
      check('Retrieval: breakout force against crane capacity', bo.F / 1e3, inp.craneCap, 'kN', 'DNV-RP-N103 (2017) §6.3 (retrieval, suction up to bearing capacity)', 'Guards against the lift force at breakout exceeding the crane or rigging capacity; the load releases suddenly once the base vents.'),
    ];
    if (!mat) checks.push(check('Skirt ring buckling under suction', sReq / 1e3, pCr / 1.5 / 1e3, 'kPa', 'Bresse–Bryan long-cylinder lower bound; DNV-RP-C202 for the code check', 'Guards against collapse of the thin skirt wall under the installation under-pressure.'));
    const ok = checks.every((c) => c.uc <= 1);
    const ballastM3 = deficit / ((2400 - 1025) * G);
    const skirtSteel = wall * hs * t * 7.85;
    const w = checks.reduce((a, c) => (c.uc > a.uc ? c : a));
    const advice = [
      ok ? adv('engineering', 'ok', `Skirts ${zSelf >= hs ? 'penetrate fully under self-weight' : `penetrate ${fmt(zSelf)} m of ${fmt(hs)} m under self-weight and need ${fmt(sReq / 1e3)} kPa of suction (or ${fmt(deficit / 1e3)} kN of ballast) for the rest`}; uplift resistance ${fmt(upVented / 1e3)} kN vented, ${fmt(upSealed / 1e3)} kN sealed.`,
        zSelf >= hs ? 'Confirm the upper-bound strength profile does not prevent penetration, and provide vent area so the water can escape at the landing speed.' : `Specify a suction system for at least ${fmt(sReq * inp.fsPlug / 1e3)} kPa under-pressure and monitor penetration against pressure; plug-heave limit is ${fmt(sAllow / 1e3)} kPa.`,
        'DNV-RP-C212 (2019) §6.2; DNV-RP-E303 (2021); DNV-RP-N103 (2017) §6',
        `Suction replaces ballast: reaching full penetration by weight alone would take ${concreteCo2Text(ballastM3)}`)
        : adv('engineering', 'fail', `"${w.name}" is not satisfied (utilisation ${fmt(w.uc)}); penetration resistance at full depth is ${fmt(full.R / 1e3)} kN against ${fmt(W / 1e3)} kN of weight.`,
          w.name.startsWith('Retrieval') ? `Vent the base and lift slowly: the vented breakout force is ${fmt(boVented.F / 1e3)} kN against ${fmt(bo.F / 1e3)} kN with β = ${inp.beta}; or inject water under the base.`
            : w.name.startsWith('Uplift') ? `Deepen the skirts or seal the base to mobilise reverse end bearing (${fmt(upSealed / 1e3)} kN sealed); alternatively add ${fmt((inp.Uup * inp.fsUp * 1e3 - upVented) / 1e3)} kN of ballast.`
              : w.name.startsWith('Skirt ring') ? `Thicken the skirt wall to about ${fmt(inp.ts * (w.uc * 1.05) ** (1 / 3))} mm (buckling pressure ∝ t³), add ring stiffeners, or reduce the required suction by adding ballast.`
                : `Reduce skirt wall length or thickness (tip resistance is ${fmt(full.tip / full.R * 100)} % of the total), shorten the skirts, or add ${fmt(Math.max(full.R - W - sLim / inp.fsPlug * A, 0) / 1e3)} kN of ballast.`,
          'DNV-RP-C212 (2019) §6.2; DNV-RP-E303 (2021); DNV-RP-N103 (2017) §6',
          w.name.startsWith('Retrieval') ? `Venting or water injection needs no material; mobilising a larger crane vessel instead would add vessel time — ${vesselCo2Text(1)}`
            : w.name.startsWith('Uplift') ? `Ballast is the heavy fix — ${concreteCo2Text(Math.max(inp.Uup * inp.fsUp * 1e3 - upVented, 0) / ((2400 - 1025) * G))} Deeper skirts or a sealed base reach the same resistance with little material.`
              : `Ballast is the heavy fix — ${concreteCo2Text(Math.max(ballastM3, 1))} Thinner or shorter skirts remove resistance and steel at the same time.`),
      adv('economic', ok ? 'info' : 'warn', `Skirt steel is ${fmt(skirtSteel)} t. ${zSelf >= hs ? 'No suction spread is required.' : 'Suction-assisted penetration needs a pump skid and ROV support, typically a few hours per foundation.'} Retrieval at end of life needs ${fmt(bo.F / 1e3)} kN of lift (β = ${inp.beta}).`,
        'Size the retrieval lift now: if the breakout force exceeds the crane class used for installation, decommissioning needs a larger vessel or a jetting / over-pressure spread. Record the as-installed penetration and pressure log for the decommissioning plan.',
        'Whole-life costing per ISO 15686-5 (decommissioning included); AACE RP 18R-97',
        `Designing for retrievability keeps the skirt steel recoverable for recycling, where an abandoned foundation is lost material — ${steelCo2Text(skirtSteel)} Each extra vessel day for a difficult recovery: ${vesselCo2Text(1)}`),
    ];
    return {
      kpis: [
        { key: 'Rfull', label: 'Penetration resistance at full skirt depth', unit: 'kN', value: full.R / 1e3, explain: 'Wall friction on both faces plus tip bearing with the skirts fully penetrated: the installation force that must be supplied by weight, ballast or suction.' },
        { key: 'zSelf', label: 'Self-weight penetration', unit: 'm', value: zSelf, explain: 'Depth reached when the hook load is released and only the submerged weight acts. Equal to the skirt height means no assistance is needed.' },
        { key: 'ballast', label: 'Additional force for full penetration', unit: 'kN', value: deficit / 1e3, explain: 'Shortfall between resistance and weight at full depth: the ballast or suction force required.' },
        { key: 'sReq', label: 'Required suction', unit: 'kPa', value: sReq / 1e3, limit: sLim / inp.fsPlug / 1e3, sense: 'max', explain: 'Under-pressure inside the skirt compartments that supplies the shortfall over the base area.' },
        { key: 'sAllow', label: 'Allowable suction (plug heave)', unit: 'kPa', value: sAllow / 1e3, explain: 'Under-pressure at which the soil plug fails upward (reverse end bearing plus inside friction). Suction above this lifts the soil instead of lowering the skirt.' },
        { key: 'sCav', label: 'Cavitation limit', unit: 'kPa', value: sCav / 1e3, explain: 'Largest under-pressure the water column allows (hydrostatic plus atmospheric minus vapour pressure); it only matters in shallow water.' },
        { key: 'upVented', label: 'Uplift resistance, vented', unit: 'kN', value: upVented / 1e3, limit: inp.Uup * inp.fsUp, sense: 'min', explain: 'Submerged weight plus wall friction with set-up: the resistance to sustained uplift when no base suction can be relied on.' },
        { key: 'upSealed', label: 'Uplift resistance, sealed', unit: 'kN', value: upSealed / 1e3, explain: 'Short-term resistance with a sealed base mobilising reverse end bearing; available only for transient loads.' },
        { key: 'breakout', label: 'Breakout force on retrieval', unit: 'kN', value: bo.F / 1e3, limit: inp.craneCap, sense: 'max', explain: 'Hook load needed to lift the foundation: weight, skirt friction and the assumed fraction of base suction.' },
        { key: 'suction', label: 'Breakout suction pressure', unit: 'kPa', value: bo.suctionPressure / 1e3, explain: 'Base under-pressure implied by the breakout factor; it vanishes if the base is vented or water is injected.' },
        { key: 'friction', label: 'Skirt wall friction at full depth', unit: 'kN', value: full.friction / 1e3, explain: 'Friction part of the penetration resistance; the remainder is tip bearing, which thinner plate reduces.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Penetration resistance against depth', x: { label: 'Force', unit: 'kN' }, y: { label: 'Tip depth below mudline', unit: 'm' },
          series: [{ name: 'Total resistance', x: Rs, y: zs }, { name: 'Wall friction', x: Fr, y: zs }, { name: 'Tip bearing', x: Tp, y: zs }, { name: 'Submerged weight', x: [W / 1e3, W / 1e3], y: [0, -hs], dash: true }],
          explain: 'Resistance to pushing the skirts in against tip depth. Where the total crosses the weight line, self-weight penetration stops and suction or ballast must take over.' },
        { type: 'line', title: 'Required and allowable suction against depth', x: { label: 'Under-pressure', unit: 'kPa' }, y: { label: 'Tip depth below mudline', unit: 'm' },
          series: curve ? [{ name: 'Required', x: curve.s.map((v) => v / 1e3), y: curve.h.map((v) => -v) }, { name: 'Plug-heave limit', x: curve.sLimit.map((v) => v / 1e3), y: curve.h.map((v) => -v) }] : [{ name: 'Required', x: sR, y: zs }, { name: 'Plug-heave limit', x: sA, y: zs }],
          explain: 'Under-pressure needed to keep the skirt moving and the pressure at which the plug would heave. The two curves must not touch before full depth.' },
        { type: 'bar', title: 'Uplift and breakout forces', categories: ['Submerged weight', 'Uplift, vented', 'Uplift, sealed', 'Breakout (vented)', `Breakout (β = ${inp.beta})`, 'Crane capacity'],
          series: [{ name: 'Force (kN)', y: [W, upVented, upSealed, boVented.F, bo.F, inp.craneCap * 1e3].map((v) => round(v / 1e3)) }], y: { label: 'Force', unit: 'kN' },
          explain: 'Resistance to uplift in service next to the force needed to recover the foundation. A breakout bar above the crane bar means retrieval needs venting or a larger vessel.' },
      ],
      tables: [{ title: 'Geometry and resistance components', columns: [{ key: 'k', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
        rows: [['Base area', A, 'm²'], ['Total skirt length in plan', wall, 'm'], ['Tip bearing at full depth', full.tip / 1e3, 'kN'], ['Wall friction at full depth', full.friction / 1e3, 'kN'], ['Penetration adhesion factor', clay ? alpha : inp.kf, '-'], ['Strength at tip level', suTip / 1e3, 'kPa'],
          ['Set-up wall friction', fr / 1e3, 'kN'], ['Base suction force at breakout', bo.suction / 1e3, 'kN'], ['External hydrostatic pressure (absolute)', pHyd / 1e3, 'kPa'], ...(Number.isFinite(pCr) ? [['Ring buckling pressure of the skirt', pCr / 1e3, 'kPa']] : [])].map(([kk, v, u]) => ({ k: kk, v: round(v), u })),
        explain: 'The geometric quantities and each resistance component behind the headline numbers, for comparison with installation records (penetration against pressure).' }],
      notes: [zSelf >= hs ? 'Full self-weight penetration is predicted.' : `Self-weight penetration stops at ${fmt(zSelf)} m; ${fmt(deficit / 1e3)} kN more is needed.`, `Tip bearing is ${fmt(full.tip / full.R * 100)} % of the penetration resistance.`, `"${w.name}" is the governing check (${fmt(w.uc)}).`],
      advice,
      publish: { 'subsea.skirtPenetrationForce': full.R, 'subsea.breakoutForce': bo.F, 'subsea.upliftResistance': upVented },
      quantities: deficit > 0 && !ok ? { concrete_m3: round(ballastM3) } : {},
    };
  },
  verify: [
    { name: 'Skirt resistance in uniform clay: R = α·su·2·z·l + Nc·su·l·t', source: 'Hand calculation: friction 0.5·5·0.5·48·2 = 120 kN; tip 7.5·5·48·0.02 = 36 kN; R = 156 kN',
      inputs: { type: 'mat', B: 10, L: 14, ribs: 0, hs: 0.5, ts: 20, soil: 'clay', su0: 5, suGrad: 0, St: 2, NcTip: 7.5, gammaSub: 6, W: 1000 },
      expect: [{ kpi: 'friction', value: 120, tol: 1e-6 }, { kpi: 'Rfull', value: 156 + 6 * 0.5 * 48 * 0.02, tol: 1e-6 }, { kpi: 'zSelf', value: 0.5, tol: 1e-6 }] },
    { name: 'Breakout: F = W′ + friction + β·Nc·su·A', source: 'Hand calculation: friction 0.65·5·0.5·48·2 = 156 kN; suction 0.5·6.17·5·140 = 2159.6 kN; F = 1000 + 156 + 2159.6 = 3315.6 kN',
      inputs: { type: 'mat', B: 10, L: 14, ribs: 0, hs: 0.5, ts: 20, soil: 'clay', su0: 5, suGrad: 0, St: 2, W: 1000, setup: 0.65, beta: 0.5 },
      expect: [{ kpi: 'breakout', value: 1000 + 156 + 0.5 * 1.2 * (2 + Math.PI) * 5 * 140, tol: 1e-5 }] },
  ],
};
