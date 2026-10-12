// Fire response of the topsides frame: steel temperature histories of the members in the fire zone by the step
// equations of EN 1993-1-2 (unprotected or with passive fire protection), strength and stiffness reduction,
// frame analysis with the degraded stiffness, residual capacity, critical temperature and fire resistance time.
import { fireCurve, steelTemperatureUnprotected, steelTemperatureProtected, steelReductionFactors, criticalTemperatureC, tubeSectionFactor, timeToTemperature, pfpThickness, restrainedThermalStress, steelThermalStrain } from '../../engine/structural/thermal.js';
import { deckInputs, solverInput, buildDeck, loadCase, solveDeck, memberForces, memberChecks, groupSummary, framePlot, equipmentOf, need, run, fmt, round, clamp, adv, check, fsCheck, thin, steelCo2Text, CONFIRM, DECK_LIMITS, NOT_FE } from './_common.js';

const interp = (xs, ys, x) => { if (x <= xs[0]) return ys[0]; for (let i = 1; i < xs.length; i++) if (x <= xs[i]) return ys[i - 1] + ((ys[i] - ys[i - 1]) * (x - xs[i - 1])) / (xs[i] - xs[i - 1]); return ys[ys.length - 1]; };

export default {
  id: 's02-fire-response', suite: 2, title: 'Fire response and passive fire protection',
  summary: 'Steel temperatures in the fire zone, strength and stiffness degradation, deformation and residual capacity of the heated frame, critical temperature, fire resistance time and required protection thickness.',
  covers: [10], fidelity: 'numerical',
  standards: ['EN 1993-1-2:2005 §3.2 (Table 3.1), §4.2.4, §4.2.5 (step equations)', 'EN 1991-1-2:2002 §3.1–3.2 (nominal fire curves, heat flux)', 'ISO 19901-3:2014 §8.7 (fire)', 'ISO 13702:2015 (control and mitigation of fires and explosions)', 'NORSOK N-004 (2022) Annex A.7'],
  models: ['Transient heat transfer (lumped-capacitance step equation)', 'Temperature-dependent material degradation (EN 1993-1-2 Table 3.1)', 'Thermal expansion stress model', 'Matrix stiffness method with reduced modulus', 'Critical temperature method', 'Residual capacity at temperature'],
  limits: `${DECK_LIMITS} Steel temperatures use the EN 1993-1-2 step equations with uniform temperature per member and the Table 3.1 reduction factors; the nominal ISO 834, hydrocarbon and external curves are built in and any other exposure (jet fire, CFD fire) is entered as a gas temperature history — no jet-fire curve is built in. The hot frame analysis is linear with reduced modulus: plastic redistribution, creep, catenary action and restrained-expansion buckling are not modelled. ${NOT_FE}`,
  searchTerms: 'offshore topsides fire response passive fire protection critical temperature hydrocarbon fire steel',
  inputs: [
    ...deckInputs(),
    { key: 'fire', label: 'Fire exposure', type: 'select', default: 'hydrocarbon', group: 'Fire', role: 'bc', options: [{ value: 'hydrocarbon', label: 'Hydrocarbon pool fire (EN 1991-1-2)' }, { value: 'iso834', label: 'Standard cellulosic (ISO 834)' }, { value: 'external', label: 'External fire curve' }, { value: 'custom', label: 'Gas temperature history (table)' }],
      help: 'Jet fires are not built in: enter the gas temperature history (typically 1100–1200 °C within a minute) and a convection coefficient for the project jet-fire definition.' },
    { key: 'gasT', label: 'Gas temperature history', type: 'series', xLabel: 'Time', xUnit: 'min', yLabel: 'Gas temperature', yUnit: '°C', default: { x: [0, 1, 60, 120], y: [20, 1150, 1150, 1150] }, group: 'Fire', role: 'bc', showIf: { key: 'fire', in: ['custom'] } },
    { key: 'alphaC', label: 'Convection coefficient (custom history)', unit: 'W/m²K', type: 'number', default: 50, min: 5, max: 500, group: 'Fire', role: 'bc', showIf: { key: 'fire', in: ['custom'] }, help: '25 for cellulosic, 50 for hydrocarbon pool fires; jet fires are higher (project value).' },
    { key: 'tReq', label: 'Required fire resistance time', unit: 'min', type: 'number', default: 60, min: 5, max: 240, group: 'Fire', role: 'criteria', help: `Time the structure must keep its function (escape, shutdown, blowdown). ${CONFIRM}` },
    { key: 'zFire0', label: 'Fire zone from elevation', unit: 'm', type: 'number', default: 4.5, min: -1, max: 100, group: 'Fire', role: 'bc', help: 'Members whose mid-point lies between the two elevations (above the jacket interface) are heated. The default band is the columns and braces standing on the lowest deck and the deck above it.' },
    { key: 'zFire1', label: 'Fire zone to elevation', unit: 'm', type: 'number', default: 12, min: 0, max: 120, group: 'Fire', role: 'bc' },
    { key: 'theta0', label: 'Steel temperature before the fire', unit: '°C', type: 'number', default: 20, min: -20, max: 80, group: 'Fire', role: 'ic', bind: 'site.Tair' },
    { key: 'protect', label: 'Passive fire protection', type: 'select', default: 'pfp', group: 'Protection', role: 'geometry', options: [{ value: 'none', label: 'Unprotected steel' }, { value: 'pfp', label: 'Insulating coating or cladding' }] },
    { key: 'dp', label: 'Protection thickness', unit: 'mm', type: 'number', default: 15, min: 1, max: 200, group: 'Protection', role: 'geometry', design: true, showIf: { key: 'protect', in: ['pfp'] } },
    { key: 'lamP', label: 'Protection thermal conductivity', unit: 'W/mK', type: 'number', default: 0.15, min: 0.01, max: 2, group: 'Protection', role: 'material', showIf: { key: 'protect', in: ['pfp'] }, uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Effective value at fire temperature from the product type approval (lightweight cementitious 0.1–0.2; char of an epoxy intumescent is product-specific).' },
    { key: 'rhoP', label: 'Protection density', unit: 'kg/m³', type: 'number', default: 550, min: 0, max: 3000, group: 'Protection', role: 'material', showIf: { key: 'protect', in: ['pfp'] } },
    { key: 'cP', label: 'Protection specific heat', unit: 'J/kgK', type: 'number', default: 1100, min: 0, max: 3000, group: 'Protection', role: 'material', showIf: { key: 'protect', in: ['pfp'] } },
    { key: 'psi', label: 'Fraction of live load present in the fire', unit: '-', type: 'number', default: 0.5, min: 0, max: 1, group: 'Acceptance criteria', role: 'criteria', help: `Accidental combination G + ψ Q with all partial factors 1.0. ${CONFIRM}` },
    { key: 'expand', label: 'Include restrained thermal expansion in the frame analysis', type: 'bool', default: false, group: 'Acceptance criteria', role: 'control', help: 'Elastic restraint forces are an upper bound that real members shed by yielding and bowing; leave off for the strength check and read the restrained-stress indicator instead.' },
    { key: 'dtFire', label: 'Time step of the temperature integration', unit: 's', type: 'number', default: 5, min: 0.25, max: 5, group: 'Solver', role: 'control', help: 'EN 1993-1-2 §4.2.5.1(4) limits the step to 5 s for unprotected steel.' },
    solverInput(1),
  ],
  mesh: { key: 'dtFire', label: 'Time step', levels: [5, 2.5, 1.25], qoi: 'thetaMax', isSize: true },
  compute(inp, ctx) {
    // SI: min → s, mm → m. Temperatures stay in °C as the engine arguments require.
    const deck = buildDeck(inp, ctx), eqs = equipmentOf(deck, inp), st = deck.steel;
    need(inp.zFire1 > inp.zFire0, 'The top of the fire zone must be above its bottom.');
    const heated = deck.orig.map((o) => o.mid[2] >= inp.zFire0 && o.mid[2] <= inp.zFire1);
    need(heated.some(Boolean), `No member mid-point lies between ${inp.zFire0} and ${inp.zFire1} m: move the fire zone onto the structure (deck levels at ${deck.levelsZ.map((z) => fmt(z)).join(', ')} m).`);
    const dur = inp.tReq * 60;
    let fire = inp.fire, alphaC;
    if (inp.fire === 'custom') {
      const xs = (inp.gasT?.x || []).map((v) => Number(v) * 60), ys = (inp.gasT?.y || []).map(Number);
      need(xs.length >= 2 && xs.length === ys.length && xs.every(Number.isFinite) && ys.every(Number.isFinite), 'The gas temperature history needs at least two numeric (time, temperature) points.');
      for (let i = 1; i < xs.length; i++) need(xs[i] > xs[i - 1], 'The times of the gas temperature history must increase.');
      fire = (t) => interp(xs, ys, t); alphaC = inp.alphaC;
    }
    // Section factor Am/V: I-sections heated on three sides (the deck plate shields the top flange), tubes on the outside.
    const sf = (o) => { const s = o.sec; if (o.kind === 'tube') return tubeSectionFactor({ D: s.D, t: s.t }); if (o.kind === 'I') return (2 * s.h + (o.horizontal ? 3 : 4) * s.b - 2 * s.tw) / s.A; return 4 / Math.sqrt(s.A); };
    const classes = new Map();
    deck.orig.forEach((o, i) => { if (!heated[i]) return; const f = clamp(Math.round(sf(o) / 5) * 5, 10, 1000), key = `${f}`; if (!classes.has(key)) classes.set(key, { f, members: [], group: o.group }); classes.get(key).members.push(i); });
    const pfp = inp.protect === 'pfp';
    const history = (f, duration) => run(() => (pfp
      ? steelTemperatureProtected({ sectionFactor: f, dp: inp.dp / 1e3, lambdaP: inp.lamP, rhoP: inp.rhoP, cP: inp.cP, fire, duration, dt: Math.min(inp.dtFire * 6, 30), rho: st.rho, theta0C: inp.theta0 })
      : steelTemperatureUnprotected({ sectionFactor: f, fire, duration, dt: inp.dtFire, alphaC, rho: st.rho, theta0C: inp.theta0 })), 'Steel temperature');
    // Histories run to twice the required time so that the resistance time can be read beyond the requirement.
    for (const c of classes.values()) { c.h = history(c.f, 2 * dur); c.theta = interp(Array.from(c.h.t), Array.from(c.h.steelC), dur); c.red = steelReductionFactors(c.theta); }
    const cls = [...classes.values()], hot = cls.reduce((a, c) => (c.theta > a.theta ? c : a));
    const memberClass = new Array(deck.orig.length).fill(null); cls.forEach((c) => c.members.forEach((i) => { memberClass[i] = c; }));
    // Accidental combination, cold and hot (reduced modulus in the heated members).
    const lc = (name, thermal) => { const l = loadCase(deck, inp, { name, gLive: inp.psi, eqs }); if (thermal) deck.orig.forEach((o, i) => { if (memberClass[i]) l.member.push({ member: i, type: 'thermal', dT: memberClass[i].theta - inp.theta0 }); }); return l; };
    const cold = solveDeck(deck, [lc('cold', false)], inp.nSub);
    const materials = { ...deck.model.materials };
    cls.forEach((c, k) => { c.mat = `hot${k}`; materials[c.mat] = { ...deck.model.materials.steel, E: st.E * Math.max(c.red.kE, 1e-3), G: undefined }; });
    const hotDeck = { ...deck, model: { ...deck.model, materials, members: deck.model.members.map((m, i) => (memberClass[i] ? { ...m, material: memberClass[i].mat } : m)) } };
    const hotSol = solveDeck(hotDeck, [lc('hot', inp.expand)], inp.nSub);
    const rc = cold.sol.results.cold, rh = hotSol.sol.results.hot;
    const fc = memberForces(deck, cold.ref, rc), fh = memberForces(deck, hotSol.ref, rh);
    const coldRows = memberChecks(deck, fc), rows = [];
    fh.forEach((f) => { const c = memberClass[f.i]; rows.push(memberChecks(deck, [f], c ? { fyFac: Math.max(c.red.ky, 1e-3), EFac: Math.max(c.red.kE, 1e-3) } : {})[0]); });
    const govHot = rows.reduce((a, c) => (c.uc > a.uc ? c : a)), groups = groupSummary(rows);
    // Critical temperature of each heated member from its cold utilisation, and the time at which it is reached.
    let crit = null;
    coldRows.forEach((r) => { const c = memberClass[r.i]; if (!c) return; const mu0 = clamp(r.uc, 0.013, 1), thCr = r.uc >= 1 ? inp.theta0 : criticalTemperatureC(mu0, 'table'), tF = r.uc >= 1 ? 0 : timeToTemperature(c.h, thCr); if (!crit || tF < crit.tF || (tF === crit.tF && thCr < crit.thCr)) crit = { i: r.i, mu0: r.uc, thCr, tF, c, group: r.group }; });
    const tFail = Number.isFinite(crit.tF) ? crit.tF / 60 : 2 * inp.tReq;
    const beyond = !Number.isFinite(crit.tF);
    const req = crit.mu0 < 1 ? run(() => pfpThickness({ sectionFactor: crit.c.f, lambdaP: inp.lamP, rhoP: inp.rhoP, cP: inp.cP, fire, duration: dur, limitC: crit.thCr, dpMax: 0.5 }), 'Protection thickness') : { dp: 0.5, required: true };
    const nN = hotSol.ref.model.nodes.length, uH = rh.displacements, uC = rc.displacements;
    let dHot = 0, dCold = 0; for (let i = 0; i < nN; i++) { dHot = Math.max(dHot, Math.abs(uH[6 * i + 2])); dCold = Math.max(dCold, Math.abs(uC[6 * i + 2])); }
    const nodeT = new Array(nN).fill(inp.theta0); hotSol.ref.nodeMap.forEach((ns, i) => { if (memberClass[i]) for (const n of ns) nodeT[n] = Math.max(nodeT[n], memberClass[i].theta); });
    const eps = steelThermalStrain(Math.max(hot.theta, 20)) - steelThermalStrain(Math.max(inp.theta0, 20));
    const rs = restrainedThermalStress({ E: st.E * Math.max(hot.red.kE, 1e-3), thermalStrain: eps, restraint: 1, fy: st.fy * Math.max(hot.red.ky, 1e-3) });
    const [ht, hg, hs] = thin([hot.h.t, hot.h.gasC, hot.h.steelC], 300);
    const gasEnd = typeof fire === 'function' ? fire(dur) : fireCurve(fire, dur).gasC;
    // Protected surface of the heated members (for quantities and advice).
    const surf = deck.orig.reduce((s, o, i) => (heated[i] ? s + sf(o) * o.sec.A * o.L : s), 0);
    const checks = [
      check('Member strength at temperature', govHot.uc, 1, 'interaction', `${govHot.ref} with fy,θ = ky,θ fy and Ea,θ = kE,θ E (EN 1993-1-2 §4.2.3, γM,fi = 1.0)`, `Strength and buckling utilisation of the most loaded member (${govHot.group.toLowerCase()} ${govHot.i + 1}) at the end of the required fire duration. Above 1.0 the member yields or buckles and the load must find another path or the deck collapses locally.`),
      fsCheck('Fire resistance time', Math.max(tFail, 1e-3), inp.tReq, 'EN 1993-1-2 §4.2.4 (critical temperature) with the §4.2.5 temperature history', `Time at which the first heated member (${crit.group.toLowerCase()} ${crit.i + 1}, cold utilisation ${fmt(crit.mu0)}) reaches its critical temperature of ${fmt(crit.thCr)} °C. It must exceed the time needed for escape, shutdown and blowdown.`),
    ];
    if (pfp) checks.push(check('Protection thickness', req.dp * 1e3, inp.dp, 'mm', 'EN 1993-1-2 §4.2.5.2 Eq. (4.27), thickness for θa ≤ θa,cr at the required time', 'Thickness of insulation that keeps the governing member at or below its critical temperature for the required duration, against the thickness provided.'));
    const gov = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = gov.uc <= 1;
    return {
      kpis: [
        { key: 'thetaMax', label: 'Highest steel temperature at the required time', unit: '°C', value: hot.theta, explain: `Temperature of the fastest-heating member class (section factor ${hot.f} 1/m) after ${inp.tReq} minutes. Carbon steel keeps full yield strength to 400 °C and about half at 590 °C.` },
        { key: 'gasEnd', label: 'Gas temperature at the required time', unit: '°C', value: gasEnd, explain: 'Temperature of the fire gases from the selected curve at the end of the required duration.' },
        { key: 'ky', label: 'Yield strength retained (ky,θ)', unit: '-', value: hot.red.ky, explain: 'EN 1993-1-2 Table 3.1 reduction factor on yield strength for the hottest member class.' },
        { key: 'kE', label: 'Stiffness retained (kE,θ)', unit: '-', value: hot.red.kE, explain: 'Reduction factor on Young’s modulus: stiffness falls faster than strength, so deflection and buckling degrade first.' },
        { key: 'ucHot', label: 'Highest utilisation at temperature', unit: '-', value: round(govHot.uc, 5), limit: 1, sense: 'max', explain: 'Member utilisation under the accidental combination with reduced yield strength and modulus. Its inverse is the residual capacity margin.' },
        { key: 'residual', label: 'Residual capacity factor', unit: '-', value: govHot.uc > 0 ? Math.min(1 / govHot.uc, 99) : 99, limit: 1, sense: 'min', explain: 'Factor by which the accidental loads could be scaled before the hottest member reaches its reduced resistance.' },
        { key: 'thetaCr', label: 'Critical temperature of the governing member', unit: '°C', value: crit.thCr, explain: 'Temperature at which the retained yield strength equals the utilisation of the member under the accidental loads.' },
        { key: 'tFail', label: 'Fire resistance time', unit: 'min', value: tFail, limit: inp.tReq, sense: 'min', explain: `Time to reach the critical temperature in the first member${beyond ? ` (not reached within ${2 * inp.tReq} minutes — reported at the end of the computed history)` : ''}.` },
        { key: 'dpReq', label: 'Required protection thickness', unit: 'mm', value: req.required ? req.dp * 1e3 : 0, explain: 'Insulation thickness (with the entered conductivity) that keeps the governing member at its critical temperature for the required time; zero means unprotected steel is sufficient.' },
        { key: 'dHot', label: 'Largest vertical deflection in the fire', unit: 'mm', value: dHot * 1e3, explain: `Vertical displacement with the reduced modulus (cold value ${fmt(dCold * 1e3)} mm). The ratio shows the stiffness lost.` },
        { key: 'sigTh', label: 'Fully restrained thermal stress', unit: 'MPa', value: Math.abs(rs.stress) / 1e6, explain: `Compressive stress in a heated member whose expansion is fully prevented (elastic value ${fmt(Math.abs(rs.elasticStress) / 1e6)} MPa, capped at the reduced yield). It shows why restrained members bow or yield early in a fire.` },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Gas and steel temperature', x: { label: 'Time', unit: 'min' }, y: { label: 'Temperature', unit: '°C' },
          series: [{ name: 'Fire gases', x: ht.map((v) => round(v / 60)), y: hg.map((v) => round(v)) }, ...cls.slice(0, 5).map((c) => { const [a, b] = thin([c.h.t, c.h.steelC], 300); return { name: `Steel, ${c.group.toLowerCase()} (Am/V ${c.f} 1/m)`, x: a.map((v) => round(v / 60)), y: b.map((v) => round(v)) }; }),
            { name: 'Critical temperature', x: [0, 2 * inp.tReq], y: [round(crit.thCr), round(crit.thCr)], dash: true }, { name: 'Required time', x: [inp.tReq, inp.tReq], y: [0, round(Math.max(...hg))], dash: true }],
          explain: 'Temperature of the fire gases and of each heated member class. The fire resistance time is where the governing steel curve crosses the dashed critical temperature; it should lie to the right of the required-time line.' },
        { type: 'line', title: 'Strength and stiffness retained by the hottest member', x: { label: 'Time', unit: 'min' }, y: { label: 'Reduction factor', unit: '-' },
          series: [{ name: 'ky,θ (yield)', x: ht.map((v) => round(v / 60)), y: hs.map((v) => round(steelReductionFactors(v).ky)) }, { name: 'kE,θ (modulus)', x: ht.map((v) => round(v / 60)), y: hs.map((v) => round(steelReductionFactors(v).kE)) }, { name: 'Utilisation of the governing member', x: [0, 2 * inp.tReq], y: [round(Math.min(crit.mu0, 1)), round(Math.min(crit.mu0, 1))], dash: true }],
          explain: 'Yield strength and modulus retained over time. The member fails when the yield curve drops to its utilisation (dashed line); the modulus curve shows that deflections grow well before that.' },
        framePlot({ title: 'Steel temperature and deformed shape at the required time', ref: hotSol.ref, u: uH, values: nodeT, valueLabel: 'Steel temperature', valueUnit: '°C', explain: 'Frame coloured by steel temperature with the deflected shape under the accidental loads and reduced stiffness. Check that the heated members are the ones the fire scenario reaches.' }),
        { type: 'bar', title: 'Utilisation at temperature by member group', categories: groups.map((g) => g.group), series: [{ name: 'In the fire', y: groups.map((g) => round(g.uc)) }, { name: 'Before the fire', y: groups.map((g) => round(Math.max(...coldRows.filter((r) => r.group === g.group).map((r) => r.uc)))) }], y: { label: 'Utilisation', unit: '-' }, explain: 'Governing utilisation of each member group before the fire and at the required time. The growth of a bar is the strength lost to temperature (and force redistribution toward cooler members).' },
      ],
      tables: [{ title: 'Heated member classes', columns: [{ key: 'group', label: 'Typical member' }, { key: 'n', label: 'Members' }, { key: 'f', label: 'Section factor', unit: '1/m' }, { key: 'theta', label: 'Steel temperature', unit: '°C' }, { key: 'ky', label: 'ky,θ' }, { key: 'kE', label: 'kE,θ' }, { key: 'uc', label: 'Utilisation at temperature' }],
        rows: cls.map((c) => ({ group: c.group, n: c.members.length, f: c.f, theta: round(c.theta), ky: round(c.red.ky), kE: round(c.red.kE), uc: round(Math.max(...rows.filter((r) => memberClass[r.i] === c).map((r) => r.uc))) })),
        explain: 'One row per section-factor class in the fire zone at the required time: thin, high section-factor members heat fastest and normally govern the protection thickness.' }],
      notes: [`${heated.filter(Boolean).length} member(s) between ${inp.zFire0} and ${inp.zFire1} m are heated; ${govHot.group.toLowerCase()} ${govHot.i + 1} governs at ${fmt(govHot.uc)} (${govHot.ref}).`,
        `Governing critical temperature ${fmt(crit.thCr)} °C (${crit.group.toLowerCase()} ${crit.i + 1}, utilisation before the fire ${fmt(crit.mu0)}), ${beyond ? `not reached within ${2 * inp.tReq} min` : `reached after ${fmt(tFail)} min`}.`,
        inp.expand ? 'Restrained thermal expansion is included elastically: member forces are an upper bound.' : 'Thermal expansion forces are not included in the member checks (see the restrained thermal stress indicator).'],
      advice: [
        ok ? adv('engineering', gov.uc > 0.85 ? 'warn' : 'ok', `The frame keeps its function for ${fmt(tFail)} min${beyond ? ' or more' : ''} against ${inp.tReq} min required; hottest steel ${fmt(hot.theta)} °C, utilisation at temperature ${fmt(govHot.uc)}.`,
          pfp ? `The governing member needs ${fmt(req.dp * 1e3)} mm of protection against ${inp.dp} mm provided: ${req.dp * 1e3 < 0.7 * inp.dp ? `the thickness could be reduced to about ${fmt(Math.ceil(req.dp * 1e3 * 1.1), 2)} mm on this member class, subject to the product type approval.` : 'keep the thickness and confirm it with the product type approval for the member section factor.'}` : 'Unprotected steel is sufficient for this duration; confirm the fire scenario and duration with the fire risk analysis.',
          'EN 1993-1-2:2005 §4.2.4–4.2.5; ISO 13702:2015; ISO 19901-3:2014 §8.7', `Protected surface in the fire zone about ${fmt(surf)} m²; at ${inp.dp} mm that is ${fmt(surf * inp.dp / 1e3 * inp.rhoP / 1e3)} t of coating carried for the life of the platform — rating members by their own critical temperature avoids over-coating lightly loaded steel.`)
          : adv('engineering', 'fail', `${gov.name} is not met (utilisation ${fmt(gov.uc)}): fire resistance ${fmt(tFail)} min against ${inp.tReq} min, hottest steel ${fmt(hot.theta)} °C, utilisation at temperature ${fmt(govHot.uc)}.`,
            req.dp < 0.5 ? `Apply at least ${fmt(Math.ceil(req.dp * 1e3), 2)} mm of protection (λ = ${inp.lamP} W/mK) to the ${crit.group.toLowerCase()} class, or lower the utilisation of ${crit.group.toLowerCase()} ${crit.i + 1} below ${fmt(Math.max(steelReductionFactors(hot.theta).ky, 0.02))} with a heavier section; deluge or a fire wall that shortens the exposure is the alternative.` : 'No practical insulation thickness holds the critical temperature: the governing member is overstressed before the fire or the exposure is too long — strengthen it, provide an alternative load path, or shield it with a fire wall.',
            'EN 1993-1-2:2005 §4.2.4–4.2.5; ISO 13702:2015', `Protection of ${fmt(surf)} m² at ${fmt(Math.ceil(req.dp * 1e3), 2)} mm is about ${fmt(surf * req.dp * inp.rhoP / 1e3)} t of coating; a section with a lower section factor (thicker, more compact) heats more slowly and can need less coating for the same rating — heavier steel at ${steelCo2Text(0.1 * deck.frameSteelT)} per 10 % of frame weight is the trade-off.`),
        adv('economic', 'info', `Passive fire protection is priced by area and thickness: ${fmt(surf)} m² of steel surface lies in this fire zone and each millimetre of thickness is ${fmt(surf * 1e-3 * inp.rhoP / 1e3)} t of material plus application and inspection access for life.`,
          'Optimise by member: protect primary members to their own critical temperature, leave redundant secondary steel unprotected where the frame survives without it, and weigh PFP against active protection and layout changes in the fire risk analysis.', 'ISO 13702:2015 (performance-based fire protection); whole-life costing per ISO 15686-5 (coating inspection and corrosion under PFP)'),
      ],
      quantities: { polymer_t: pfp ? surf * inp.dp / 1e3 * inp.rhoP / 1e3 : 0 },
    };
  },
  verify: [
    { // Constant gas temperature of 600 °C held for 4 h on unprotected steel: the steel reaches the gas temperature,
      // and EN 1993-1-2 Table 3.1 gives ky,θ = 0.47 and kE,θ = 0.31 at 600 °C.
      name: 'Steel soaked at 600 °C: ky,θ = 0.47, kE,θ = 0.31', source: 'EN 1993-1-2:2005 Table 3.1',
      inputs: { fire: 'custom', gasT: { x: [0, 0.1, 500], y: [600, 600, 600] }, alphaC: 50, tReq: 240, protect: 'none' }, expect: [{ kpi: 'thetaMax', value: 600, tol: 0.005 }, { kpi: 'ky', value: 0.47, tol: 0.01 }, { kpi: 'kE', value: 0.31, tol: 0.01 }] },
    { // ISO 834: θg = 20 + 345 log10(8t + 1), t in minutes → at 60 min 20 + 345·log10(481) = 945.3 °C.
      name: 'ISO 834 gas temperature at 60 minutes: 945 °C', source: 'EN 1991-1-2:2002 Eq. (3.4)',
      inputs: { fire: 'iso834', tReq: 60 }, expect: [{ kpi: 'gasEnd', value: 945.3, tol: 0.001 }] },
  ],
};
