// Equipment support and machinery foundation: a machine on rigid mounts or isolators carried by deck support
// beams or by a block footing on an elastic half-space — natural frequencies, tuning against the running speed,
// vibration amplitudes, transmitted force, support reactions, holding-down bolts and foundation stresses.
import { unbalanceForce, transmissibility, sdofForcedResponse, isolatorSelection, twoDofFoundation, equivalentFootingRadius, footingImpedance, vibrationSeverity, campbellData } from '../../engine/structural/machinery.js';
import { iSection } from '../../engine/structural/sections.js';
import { boltResistance, BOLT_STRESS_AREA, BOLT_GRADES } from '../../engine/structural/steelcheck.js';
import { steelOf, need, run, finite, fmt, round, linspace, adv, check, fsCheck, steelCo2Text, CONFIRM, G } from './_common.js';

export default {
  id: 's02-machinery-foundation', suite: 2, title: 'Equipment support and machinery foundation',
  summary: 'Tuning, vibration amplitudes and transmitted forces of a rotating machine on its mounts and supporting structure, with support reactions, holding-down bolt loads and stresses in the support beams or bearing pressure under a block foundation.',
  covers: [7], fidelity: 'reduced-order',
  standards: ['ISO 20816-1:2016 / ISO 20816-3:2022 (vibration severity zones)', 'ISO 21940-11:2016 (balance quality)', 'API RP 686 (2009) ch. 4 (machinery foundations)', 'ACI 351.3R-18 (foundations for dynamic equipment)', 'EN 1993-1-8:2005 Table 3.4 (bolts)', 'ISO 19901-3:2014 §7.6', 'DIN 4024-1 (machine foundations)'],
  models: ['Newton’s second law (two-degree-of-freedom machine–foundation model)', 'Rotor unbalance forcing model', 'Vibration isolation (transmissibility)', 'Forced harmonic response', 'Machinery foundation impedance (elastic half-space lumped parameters)', 'Campbell diagram', 'Euler–Bernoulli beam stiffness', 'Bolted joint resistance'],
  limits: 'Vertical translation only, with the machine as one rigid mass and the support as one mass on one spring: rocking, coupled horizontal modes, skid flexibility and the higher modes of the deck are not modelled (use the deck vibration module for the framing). Footing impedances are the static half-space values of Richart, Hall & Woods with geometric damping — no frequency dependence, embedment, layering or soil material damping. Vendor-certified dynamic forces should replace the balance-grade estimate when available.',
  searchTerms: 'machinery foundation vibration isolation transmissibility rotating equipment support offshore skid',
  inputs: [
    { key: 'mMach', label: 'Machine and skid mass', unit: 't', type: 'number', default: 40, min: 0.1, max: 2000, group: 'Machine', role: 'load' },
    { key: 'mRotor', label: 'Rotor mass', unit: 't', type: 'number', default: 4, min: 0.01, max: 500, group: 'Machine', role: 'load' },
    { key: 'speed', label: 'Operating speed', unit: 'rpm', type: 'number', default: 1500, min: 30, max: 30000, group: 'Machine', role: 'bc', design: true },
    { key: 'speedMin', label: 'Lowest continuous speed', unit: 'rpm', type: 'number', default: 1200, min: 0, max: 30000, group: 'Machine', role: 'bc' },
    { key: 'grade', label: 'Balance quality grade G', unit: 'mm/s', type: 'number', default: 6.3, min: 0.4, max: 630, group: 'Machine', role: 'load', help: 'ISO 21940-11: G 2.5 turbines and compressors, G 6.3 pumps and general machinery, G 16 and above reciprocating.' },
    { key: 'service', label: 'Service factor on the unbalance', unit: '-', type: 'number', default: 2.5, min: 1, max: 10, group: 'Machine', role: 'load', help: 'Allowance for the balance deteriorating in service (API RP 686 and ACI 351.3R use factors of 2 and more).' },
    { key: 'Fextra', label: 'Additional harmonic force at running speed', unit: 'kN', type: 'number', default: 0, min: 0, max: 5000, group: 'Machine', role: 'bc', help: 'Certified vendor force (reciprocating unbalance, gas forces) acting vertically.' },
    { key: 'torque', label: 'Fault or start-up torque', unit: 'kN·m', type: 'number', default: 150, min: 0, max: 1e5, group: 'Machine', role: 'bc', help: 'Short-circuit or locked-rotor torque reacted by the supports as a couple across the footprint width.' },
    { key: 'hCog', label: 'Centre of gravity above the supports', unit: 'm', type: 'number', default: 1.4, min: 0, max: 15, group: 'Machine', role: 'geometry' },
    { key: 'aH', label: 'Horizontal design acceleration', unit: 'g', type: 'number', default: 0.1, min: 0, max: 2, group: 'Machine', role: 'bc', help: 'Platform motion, seismic or transport acceleration acting on the machine with the operating loads.' },
    { key: 'bFoot', label: 'Footprint width between support lines', unit: 'm', type: 'number', default: 2.4, min: 0.3, max: 20, group: 'Machine', role: 'geometry' },
    { key: 'nSup', label: 'Number of support points', type: 'int', default: 8, min: 3, max: 40, group: 'Machine', role: 'geometry', help: 'Holding-down positions, half on each support line.' },
    { key: 'mount', label: 'Mounting', type: 'select', default: 'iso', group: 'Mounts', role: 'bc', options: [{ value: 'rigid', label: 'Rigid (bolted or grouted)' }, { value: 'iso', label: 'Anti-vibration mounts' }] },
    { key: 'kIso', label: 'Total vertical stiffness of the mounts', unit: 'kN/mm', type: 'number', default: 60, min: 0.1, max: 1e5, group: 'Mounts', role: 'bc', design: true, showIf: { key: 'mount', in: ['iso'] }, uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'zIso', label: 'Damping ratio of the mounts', unit: '%', type: 'number', default: 5, min: 0, max: 30, group: 'Mounts', role: 'ic', showIf: { key: 'mount', in: ['iso'] } },
    { key: 'support', label: 'Supporting structure', type: 'select', default: 'deck', group: 'Supporting structure', role: 'bc', options: [{ value: 'deck', label: 'Deck support beams' }, { value: 'block', label: 'Block foundation on soil (half-space)' }] },
    { key: 'nBeams', label: 'Support beams', type: 'int', default: 2, min: 1, max: 12, group: 'Supporting structure', role: 'geometry', showIf: { key: 'support', in: ['deck'] } },
    { key: 'Lb', label: 'Support beam span', unit: 'm', type: 'number', default: 6, min: 1, max: 30, group: 'Supporting structure', role: 'geometry', showIf: { key: 'support', in: ['deck'] }, design: true },
    { key: 'ends', label: 'Support beam ends', type: 'select', default: 'fixed', group: 'Supporting structure', role: 'bc', showIf: { key: 'support', in: ['deck'] }, options: [{ value: 'simple', label: 'Simply supported' }, { value: 'fixed', label: 'Continuous / fixed' }] },
    { key: 'h', label: 'Beam depth', unit: 'mm', type: 'number', default: 800, min: 100, max: 3000, group: 'Supporting structure', role: 'geometry', showIf: { key: 'support', in: ['deck'] }, design: true },
    { key: 'b', label: 'Beam flange width', unit: 'mm', type: 'number', default: 300, min: 50, max: 1200, group: 'Supporting structure', role: 'geometry', showIf: { key: 'support', in: ['deck'] } },
    { key: 'tw', label: 'Beam web thickness', unit: 'mm', type: 'number', default: 15, min: 4, max: 80, group: 'Supporting structure', role: 'geometry', showIf: { key: 'support', in: ['deck'] } },
    { key: 'tf', label: 'Beam flange thickness', unit: 'mm', type: 'number', default: 30, min: 6, max: 120, group: 'Supporting structure', role: 'geometry', showIf: { key: 'support', in: ['deck'] } },
    { key: 'mSup', label: 'Mass moving with the support', unit: 't', type: 'number', default: 12, min: 0.01, max: 5000, group: 'Supporting structure', role: 'ic', help: 'Deck: about half the mass of the support beams plus plating and items on them. Block: the concrete block.' },
    { key: 'zSup', label: 'Structural damping ratio of the deck support', unit: '%', type: 'number', default: 2, min: 0.2, max: 10, group: 'Supporting structure', role: 'ic', showIf: { key: 'support', in: ['deck'] } },
    { key: 'kUser', label: 'Support stiffness override (0 = computed)', unit: 'kN/mm', type: 'number', default: 0, min: 0, max: 1e12, group: 'Supporting structure', role: 'bc', help: 'Measured or finite-element dynamic stiffness at the machine position; replaces the beam or footing formula when above zero.' },
    { key: 'Bf', label: 'Block width', unit: 'm', type: 'number', default: 4, min: 0.5, max: 40, group: 'Supporting structure', role: 'geometry', showIf: { key: 'support', in: ['block'] } },
    { key: 'Lf', label: 'Block length', unit: 'm', type: 'number', default: 8, min: 0.5, max: 60, group: 'Supporting structure', role: 'geometry', showIf: { key: 'support', in: ['block'] } },
    { key: 'Gs', label: 'Soil shear modulus', unit: 'MPa', type: 'number', default: 50, min: 2, max: 2000, group: 'Supporting structure', role: 'soil', showIf: { key: 'support', in: ['block'] }, uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'Small-strain value from shear-wave velocity or correlation.' },
    { key: 'nuS', label: 'Soil Poisson’s ratio', unit: '-', type: 'number', default: 0.35, min: 0.1, max: 0.49, group: 'Supporting structure', role: 'soil', showIf: { key: 'support', in: ['block'] } },
    { key: 'rhoS', label: 'Soil density', unit: 'kg/m³', type: 'number', default: 1900, min: 1200, max: 2500, group: 'Supporting structure', role: 'soil', showIf: { key: 'support', in: ['block'] } },
    { key: 'qAllow', label: 'Allowable bearing pressure', unit: 'kPa', type: 'number', default: 150, min: 20, max: 2000, group: 'Supporting structure', role: 'criteria', showIf: { key: 'support', in: ['block'] }, help: 'For machine foundations commonly limited to half the static allowable value.' },
    { key: 'mat', label: 'Steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Supporting structure', role: 'material' },
    { key: 'bolt', label: 'Holding-down bolt', type: 'select', default: 'M30', group: 'Holding-down bolts', role: 'geometry', options: Object.keys(BOLT_STRESS_AREA).map((v) => ({ value: v, label: v })) },
    { key: 'boltGrade', label: 'Bolt grade', type: 'select', default: '8.8', group: 'Holding-down bolts', role: 'material', options: Object.keys(BOLT_GRADES).map((v) => ({ value: String(v), label: String(v) })) },
    { key: 'nBolt', label: 'Bolts per support point', type: 'int', default: 2, min: 1, max: 12, group: 'Holding-down bolts', role: 'geometry' },
    { key: 'marginReq', label: 'Required separation margin', unit: '%', type: 'number', default: 20, min: 5, max: 50, group: 'Acceptance criteria', role: 'criteria', help: `Between each system natural frequency and the running-speed range. ${CONFIRM}` },
    { key: 'vAB', label: 'Vibration zone A/B boundary (rms)', unit: 'mm/s', type: 'number', default: 2.3, min: 0.1, max: 50, group: 'Acceptance criteria', role: 'criteria', help: `ISO 20816-3 values by machine group and support flexibility. ${CONFIRM}` },
    { key: 'vBC', label: 'Vibration zone B/C boundary (rms)', unit: 'mm/s', type: 'number', default: 4.5, min: 0.2, max: 60, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'vCD', label: 'Vibration zone C/D boundary (rms)', unit: 'mm/s', type: 'number', default: 7.1, min: 0.3, max: 80, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'gM', label: 'Material factor on the support steel', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.5, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'effTarget', label: 'Target isolation efficiency', unit: '%', type: 'number', default: 90, min: 50, max: 99, group: 'Acceptance criteria', role: 'criteria' },
  ],
  compute(inp, ctx) {
    // SI: t → kg, rpm → Hz, kN/mm → N/m, mm → m, MPa → Pa, % → fraction.
    const st = steelOf(ctx, inp.mat), m1 = inp.mMach * 1e3, m2 = inp.mSup * 1e3, deckSup = inp.support === 'deck', iso = inp.mount === 'iso';
    need(inp.speedMin <= inp.speed, 'The lowest continuous speed must not exceed the operating speed.');
    need(inp.vAB < inp.vBC && inp.vBC < inp.vCD, 'The vibration zone boundaries must increase from A/B to C/D.');
    need(inp.mRotor <= inp.mMach, 'The rotor mass cannot exceed the machine mass.');
    const fe = inp.speed / 60, w = 2 * Math.PI * fe;
    const ub = unbalanceForce({ mass: inp.mRotor * 1e3, gradeMmS: inp.grade, speedRpm: inp.speed, serviceFactor: inp.service }), F0 = ub.force + inp.Fextra * 1e3;
    let sec = null, kCalc, c2, foot = null;
    if (deckSup) {
      need(2 * inp.tf < inp.h && inp.tw < inp.b, 'The support beam is not a valid I-shape: check its thicknesses.');
      sec = run(() => iSection({ h: inp.h / 1e3, b: inp.b / 1e3, tw: inp.tw / 1e3, tf: inp.tf / 1e3, welded: true }), 'Support beam section');
      kCalc = inp.nBeams * (inp.ends === 'fixed' ? 192 : 48) * st.E * sec.Iy / inp.Lb ** 3;   // midspan point load
    } else {
      const r = equivalentFootingRadius({ B: inp.Bf, L: inp.Lf });
      foot = run(() => footingImpedance({ G: inp.Gs * 1e6, nu: inp.nuS, rho: inp.rhoS, radius: r.translation, radiusRocking: r.rocking, radiusTorsion: r.torsion, mass: m1 + m2 }), 'Footing impedance');
      kCalc = foot.kz;
    }
    const k2 = inp.kUser > 0 ? inp.kUser * 1e6 : kCalc, zeta2 = deckSup ? inp.zSup / 100 : Math.min(foot.Dz ?? 0.3, 0.7);
    c2 = 2 * zeta2 * Math.sqrt(k2 * (m1 + m2));
    const k1 = inp.kIso * 1e6, c1 = iso ? 2 * (inp.zIso / 100) * Math.sqrt(k1 * m1) : 0;
    // Response at a forcing frequency f with the unbalance force growing as f² (plus the constant vendor force).
    const forceAt = (f) => ub.force * (f / fe) ** 2 + inp.Fextra * 1e3;
    const resp = (f) => {
      const F = Math.max(forceAt(f), 1e-9);
      if (iso) { const r = twoDofFoundation({ m1, k1, c1, m2, k2, c2, F0: F, frequencyHz: f }); return { x1: r.x1, x2: r.x2, Ft: r.forceToGround, Fi: r.forceInIsolators, TR: r.transmissibility, fn: r.naturalHz }; }
      const r = sdofForcedResponse({ mass: m1 + m2, k: k2, c: c2, frequencyHz: f, F0: F }); return { x1: r.amplitude, x2: r.amplitude, Ft: r.transmittedForce, Fi: F, TR: r.transmissibility, fn: [r.naturalHz] };
    };
    const op = resp(fe), fns = Array.from(op.fn).filter((v) => Number.isFinite(v) && v > 0);
    const fnIso = Math.sqrt(k1 / m1) / (2 * Math.PI), fnSup = Math.sqrt(k2 / (m1 + m2)) / (2 * Math.PI);
    const v1 = finite(op.x1 * w / Math.SQRT2, 'Machine vibration velocity'), v2 = op.x2 * w / Math.SQRT2;
    const sev = vibrationSeverity({ velocityRms: v1, limits: { AB: inp.vAB / 1e3, BC: inp.vBC / 1e3, CD: inp.vCD / 1e3 } });
    const lo = inp.speedMin / 60, margins = fns.map((f) => (f >= lo && f <= fe ? 0 : Math.min(Math.abs(lo - f), Math.abs(fe - f)) / f)), margin = Math.min(...margins), fCrit = fns[margins.indexOf(margin)];
    const sel = isolatorSelection({ efficiency: inp.effTarget / 100, frequencyHz: fe, mass: m1 });
    // Support reactions: weight, fault torque couple across the footprint, horizontal inertia at the centre of gravity
    // and the dynamic force through the mounts, shared by the support points (half on each line).
    const n = Math.max(3, Math.round(inp.nSup)), half = n / 2, W = m1 * G, H = inp.aH * W;
    const couple = (inp.torque * 1e3 + H * inp.hCog) / inp.bFoot / half, dyn = op.Fi / n;
    const Rmax = W / n + couple + dyn, Rmin = W / n - couple - dyn, shearPt = (H + F0) / n;
    const A = BOLT_STRESS_AREA[inp.bolt], dB = Number(inp.bolt.slice(1)) / 1e3;
    const br = run(() => boltResistance({ d: dB, As: A, grade: inp.boltGrade, FvEd: shearPt / inp.nBolt, FtEd: Math.max(-Rmin, 0) / inp.nBolt }), 'Bolt resistance');
    const checks = [
      fsCheck('Separation of natural frequencies from the running speed', Math.max(margin, 1e-4), inp.marginReq / 100, 'API RP 686 ch. 4 / ACI 351.3R §6 (frequency separation); project margin', `Distance between the nearest system natural frequency (${fmt(fCrit)} Hz) and the operating speed range ${fmt(lo)}–${fmt(fe)} Hz. Inside the margin the response is governed by damping alone and small changes in mass or stiffness cause large changes in vibration.`),
      check(`Machine vibration severity (zone ${sev.zone})`, v1 * 1e3, inp.vBC, 'mm/s rms', 'ISO 20816-1 / -3 zone B/C boundary', 'Predicted rms vibration velocity of the machine against the limit for unrestricted long-term operation. Above it, bearings, seals and couplings wear quickly and the protection system will alarm or trip.'),
      check('Holding-down bolts (shear and tension)', br.ucCombined, 1, 'interaction', `EN 1993-1-8 Table 3.4, ${inp.nBolt} × ${inp.bolt} grade ${inp.boltGrade} per support point`, 'Shear from horizontal inertia and the rotating force, with tension where the fault torque and overturning lift a support point. Bolt failure lets the machine walk on its supports and breaks the alignment.'),
    ];
    let sigma = 0, defl = 0, qBear = 0;
    if (deckSup) {
      const P = (W + m2 * G * 0) / inp.nBeams, Pd = op.Ft / inp.nBeams, cM = inp.ends === 'fixed' ? 8 : 4, M = (P + Pd) * inp.Lb / cM;
      sigma = M / sec.Wely; defl = W / k2;
      checks.push(check('Support beam bending stress', sigma / 1e6, st.fy / inp.gM / 1e6, 'MPa', `Elastic bending, fy/γM with γM = ${inp.gM} (EN 1993-1-1 §6.2.5; NORSOK N-004 §6)`, 'Bending stress in a support beam under the machine weight plus the transmitted dynamic force at mid-span. Machinery supports are normally sized by stiffness, so this should be low; a high value signals an undersized beam.'));
      checks.push(check('Static deflection of the support under the machine', defl * 1e3, inp.Lb * 1e3 / 1000, 'mm', 'Span / 1000 (alignment practice for rotating machinery supports)', 'Deflection of the support beams under the machine weight. Large deflection upsets shaft alignment between driver and driven machine when they sit on different beams.'));
    } else {
      qBear = ((m1 + m2) * G + op.Ft) / (inp.Bf * inp.Lf);
      checks.push(check('Bearing pressure under the block', qBear / 1e3, inp.qAllow, 'kPa', 'ACI 351.3R §4.3 (static plus dynamic bearing pressure against a reduced allowable)', 'Soil pressure under the block from the total weight plus the transmitted dynamic force. Excess pressure causes settlement and tilting that the machine alignment cannot tolerate.'));
      checks.push(check('Block mass ratio', 3 * m1 / 1e3, m2 / 1e3, 't', 'API RP 686 ch. 4 §4.3 (foundation mass at least 3 × machine mass for rotating machinery)', 'A heavy block keeps the amplitude low and the centre of gravity near the soil. A light block vibrates more and is sensitive to the soil stiffness assumed.'));
    }
    const gov = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = gov.uc <= 1;
    const fs = linspace(0.02 * fe, 1.6 * fe, 240).concat(fns.filter((f) => f < 1.6 * fe)).sort((a, b) => a - b), sw = fs.map(resp);
    const rs = linspace(0.1, 4, 80), zi = iso ? inp.zIso / 100 : zeta2;
    const camp = run(() => campbellData({ modes: fns, orders: [1, 2], speedMinRpm: 0, speedMaxRpm: 1.3 * inp.speed, operatingRangeRpm: [inp.speedMin, inp.speed], required: inp.marginReq / 100, nPoints: 41 }), 'Campbell diagram');
    const beamT = deckSup ? inp.nBeams * sec.A * inp.Lb * st.rho / 1e3 : 0, kNeed = (m1 + m2) * (2 * Math.PI * fe * (1 + inp.marginReq / 100)) ** 2;
    return {
      kpis: [
        { key: 'F0', label: 'Dynamic force at running speed', unit: 'kN', value: F0 / 1e3, explain: 'Rotating unbalance m·G·Ω with the service factor plus the vendor force: the amplitude that excites the foundation once per revolution.' },
        { key: 'fnIso', label: 'Natural frequency of the machine on its mounts', unit: 'Hz', value: iso ? fnIso : fnSup, explain: iso ? 'Vertical frequency of the machine on the anti-vibration mounts with a rigid support. Isolation needs the running frequency to be at least √2, preferably 3, times this value.' : 'Vertical frequency of the machine and support mass on the support stiffness (rigid mounting).' },
        { key: 'fn1', label: 'Lowest system natural frequency', unit: 'Hz', value: fns[0], explain: 'First vertical natural frequency of the coupled machine–mount–support system.' },
        { key: 'fn2', label: 'Second system natural frequency', unit: 'Hz', value: fns[1] ?? fns[0], explain: 'Second vertical natural frequency (the support mass on its stiffness, shifted by the machine). Equal to the first for rigid mounting.' },
        { key: 'margin', label: 'Separation margin', unit: '%', value: margin * 100, limit: inp.marginReq, sense: 'min', explain: 'Smallest relative distance between a natural frequency and the operating speed range. Zero means a resonance lies inside the range.' },
        { key: 'TR', label: 'Force transmissibility', unit: '-', value: op.TR, explain: 'Force passed into the supporting structure over the force generated by the machine. Below 1 the mounts isolate; above 1 they amplify.' },
        { key: 'Ft', label: 'Dynamic force transmitted to the structure', unit: 'kN', value: op.Ft / 1e3, explain: 'Amplitude of the alternating force entering the deck or soil at the running speed: the input for the deck vibration and fatigue checks.' },
        { key: 'vMach', label: 'Machine vibration velocity (rms)', unit: 'mm/s', value: v1 * 1e3, limit: inp.vBC, sense: 'max', explain: `Predicted rms velocity of the machine body: zone ${sev.zone} of the entered boundaries (A new, B acceptable, C restricted, D damaging).` },
        { key: 'vSup', label: 'Support vibration velocity (rms)', unit: 'mm/s', value: v2 * 1e3, explain: 'Velocity of the supporting structure under the machine; compare with work-area and adjacent-equipment limits.' },
        { key: 'xMach', label: 'Machine displacement amplitude', unit: 'µm', value: op.x1 * 1e6, explain: 'Zero-to-peak displacement of the machine at the running speed.' },
        { key: 'kSup', label: 'Support stiffness', unit: 'kN/mm', value: k2 / 1e6, explain: deckSup ? 'Vertical stiffness of the support beams at mid-span (or the override value).' : 'Static vertical spring of the block on the elastic half-space, 4Gr/(1 − ν).' },
        { key: 'Rmax', label: 'Largest support-point reaction', unit: 'kN', value: Rmax / 1e3, explain: 'Weight share plus the couple from fault torque and horizontal inertia plus the dynamic force, at the most loaded support point.' },
        { key: 'Rmin', label: 'Smallest support-point reaction', unit: 'kN', value: Rmin / 1e3, explain: 'Least loaded support point; a negative value is uplift carried by the holding-down bolts.' },
        { key: 'stress', label: deckSup ? 'Support beam bending stress' : 'Bearing pressure', unit: deckSup ? 'MPa' : 'kPa', value: deckSup ? sigma / 1e6 : qBear / 1e3, explain: deckSup ? 'Bending stress in a support beam from the machine weight and transmitted dynamic force.' : 'Soil pressure under the block from weight plus transmitted dynamic force.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Vibration velocity against speed', x: { label: 'Speed', unit: 'rpm' }, y: { label: 'Velocity (rms)', unit: 'mm/s' }, series: [{ name: 'Machine', x: fs.map((f) => round(f * 60)), y: sw.map((r, i) => round(r.x1 * 2 * Math.PI * fs[i] / Math.SQRT2 * 1e3)) }, { name: 'Support', x: fs.map((f) => round(f * 60)), y: sw.map((r, i) => round(r.x2 * 2 * Math.PI * fs[i] / Math.SQRT2 * 1e3)) }, { name: 'Zone B/C', x: [0, round(1.6 * inp.speed)], y: [inp.vBC, inp.vBC], dash: true }, { name: 'Operating speed', x: [inp.speed, inp.speed], y: [0, inp.vBC], dash: true }],
          explain: 'Run-up curve: rms velocity of the machine and of its support as the speed rises, with the unbalance force growing with speed squared. Peaks are the critical speeds passed during start-up; the operating speed should sit on a flat part below the zone boundary.' },
        { type: 'line', title: 'Transmissibility of the mounts', x: { label: 'Frequency ratio f / fn', unit: '-' }, y: { label: 'Transmissibility', unit: '-' }, series: [{ name: `Damping ${fmt(zi * 100)} %`, x: rs.map((r) => round(r)), y: rs.map((r) => round(transmissibility({ r, zeta: zi }).transmissibility)) }, { name: 'Undamped', x: rs.filter((r) => Math.abs(r - 1) > 0.08).map((r) => round(r)), y: rs.filter((r) => Math.abs(r - 1) > 0.08).map((r) => round(Math.min(transmissibility({ r, zeta: 0 }).transmissibility, 20))), dash: true }, { name: 'Operating point', x: [round(fe / (iso ? fnIso : fnSup))], y: [round(transmissibility({ r: fe / (iso ? fnIso : fnSup), zeta: zi }).transmissibility)], mode: 'markers' }],
          explain: 'Single-degree-of-freedom transmissibility of the mounting on a rigid base. Isolation starts beyond a ratio of √2; the operating point should be well to the right (ratio 3 gives about 88 % isolation undamped).' },
        { type: 'line', title: 'Campbell diagram', x: { label: 'Speed', unit: 'rpm' }, y: { label: 'Frequency', unit: 'Hz' }, series: [...camp.modeHz.map((mm, k) => ({ name: `Natural frequency ${k + 1}`, x: Array.from(camp.speedRpm, (v) => round(v)), y: Array.from(mm, (v) => round(v)) })), ...camp.orderHz.map((o, k) => ({ name: `${k + 1}× order`, x: Array.from(camp.speedRpm, (v) => round(v)), y: Array.from(o, (v) => round(v)), dash: true }))],
          explain: 'System natural frequencies against the 1× and 2× excitation lines. Crossings are critical speeds; none should lie between the lowest continuous speed and the operating speed.' },
        { type: 'bar', title: 'Support-point reactions', categories: ['Weight share', 'Torque and inertia couple', 'Dynamic', 'Largest', 'Smallest'], series: [{ name: 'Reaction', y: [round(W / n / 1e3), round(couple / 1e3), round(dyn / 1e3), round(Rmax / 1e3), round(Rmin / 1e3)] }], y: { label: 'Force per support point', unit: 'kN' }, explain: 'Build-up of the reaction at one support point: static share of the weight, the couple from fault torque and horizontal inertia, and the alternating force. The smallest value shows whether a point lifts.' },
      ],
      tables: [{ title: 'Foundation dynamic properties', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
        rows: [{ q: 'Running frequency', v: round(fe), u: 'Hz' }, { q: 'Rotor eccentricity (permissible)', v: round(ub.eccentricity * 1e6), u: 'µm' }, { q: 'Mount stiffness', v: iso ? inp.kIso : 'rigid', u: 'kN/mm' }, { q: 'Support stiffness', v: round(k2 / 1e6), u: 'kN/mm' }, { q: 'Support damping ratio', v: round(zeta2 * 100), u: '%' }, { q: 'Static deflection of the mounts', v: iso ? round(W / k1 * 1e3) : 0, u: 'mm' },
          { q: `Mount natural frequency for ${inp.effTarget} % isolation`, v: round(sel.naturalHz), u: 'Hz' }, { q: `Mount stiffness for ${inp.effTarget} % isolation`, v: round(sel.stiffness / 1e6), u: 'kN/mm' }, { q: 'Mount static deflection for that isolation', v: round(sel.staticDeflection * 1e3), u: 'mm' }, { q: 'Support stiffness for a high-tuned rigid mounting', v: round(kNeed / 1e6), u: 'kN/mm' }, { q: 'Bolt shear / tension per bolt', v: `${fmt(shearPt / inp.nBolt / 1e3)} / ${fmt(Math.max(-Rmin, 0) / inp.nBolt / 1e3)}`, u: 'kN' }],
        explain: 'Dynamic properties of the machine–mount–support system with the mount selection that would give the target isolation and the support stiffness a rigid, high-tuned mounting would need.' }],
      notes: [`${gov.name} governs at ${fmt(gov.uc)}. Natural frequencies ${fns.map((f) => fmt(f)).join(' and ')} Hz against a running frequency of ${fmt(fe)} Hz (${fe > fns[fns.length - 1] ? 'low-tuned: the machine runs above all system frequencies and passes them at start-up' : fe < fns[0] ? 'high-tuned: all system frequencies lie above the running speed' : 'the running speed lies between the two system frequencies'}).`,
        deckSup ? `Support: ${inp.nBeams} beam(s) of ${inp.Lb} m, stiffness ${fmt(k2 / 1e6)} kN/mm.` : `Block ${inp.Bf} × ${inp.Lf} m on soil with G = ${inp.Gs} MPa: kz = ${fmt(k2 / 1e6)} kN/mm, geometric damping ${fmt(zeta2 * 100)} % (static half-space values).`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.85 ? 'warn' : 'ok', `The foundation is ${fe > fns[fns.length - 1] ? 'low-tuned' : 'tuned'} with ${fmt(margin * 100)} % separation; machine vibration ${fmt(v1 * 1e3)} mm/s rms (zone ${sev.zone}); ${fmt(op.TR * 100)} % of the dynamic force reaches the structure.`,
          `Keep the arrangement. Give the vendor the support stiffness of ${fmt(k2 / 1e6)} kN/mm, confirm the certified unbalanced forces, and verify the natural frequencies by a bump test after installation.`, 'API RP 686 (2009) ch. 4; ISO 20816-3:2022; ACI 351.3R-18 §6',
          deckSup ? `The support beams weigh ${fmt(beamT)} t (${steelCo2Text(beamT)}); isolating the machine lets the deck steel be sized for strength rather than for a high-tuned stiffness of ${fmt(kNeed / 1e6)} kN/mm.` : `The block is ${fmt(m2 / 1e3)} t of concrete ≈ ${fmt(m2 / 2400)} m³ ≈ ${fmt(m2 / 2400 * 0.3)} t CO₂e at 0.3 t CO₂e/m³ (screening factor); soil stiffness measurements can justify a smaller block.`)
          : adv('engineering', 'fail', `${gov.name} is not met (utilisation ${fmt(gov.uc)}): separation ${fmt(margin * 100)} %, machine vibration ${fmt(v1 * 1e3)} mm/s rms (zone ${sev.zone}), smallest support reaction ${fmt(Rmin / 1e3)} kN.`,
            gov.name.startsWith('Separation') ? (iso ? `Move the natural frequency of ${fmt(fCrit)} Hz out of the speed range: mounts of about ${fmt(sel.stiffness / 1e6)} kN/mm give ${inp.effTarget} % isolation, or stiffen the support above ${fmt(kNeed / 1e6)} kN/mm and mount rigidly.` : `Stiffen the support to at least ${fmt(kNeed / 1e6)} kN/mm (deeper beams, shorter span or an extra column), or fit mounts of about ${fmt(sel.stiffness / 1e6)} kN/mm.`)
              : gov.name.startsWith('Machine vibration') ? `Balance to a finer grade than G ${inp.grade}, increase the damping, or change the tuning as for the separation margin.` : gov.name.startsWith('Holding-down') ? `Use ${Math.ceil(inp.nBolt * br.ucCombined)} bolts per point or a larger size, and add shear keys so the bolts carry tension only.` : deckSup ? `Deepen the support beams to about ${fmt(inp.h * Math.sqrt(gov.uc), 3)} mm or shorten the span.` : 'Enlarge the block footprint or mass, or improve the soil below it.',
            'API RP 686 (2009) ch. 4; ISO 20816-3:2022; EN 1993-1-8:2005 Table 3.4', deckSup ? `Stiffness by depth costs little steel: ${fmt((Math.sqrt(Math.max(gov.uc, 1.05)) - 1) * 100, 2)} % more depth at constant area, against ${steelCo2Text(beamT * Math.max(gov.uc - 1, 0.1))} if the same stiffness were bought with area.` : `Each extra cubic metre of block is about 2.4 t of concrete ≈ 0.3 t CO₂e (screening factor); soil improvement or isolators may reach the same result with less material.`),
        adv('economic', 'info', `Vibration zone ${sev.zone} at ${fmt(v1 * 1e3)} mm/s: machines in zone C or D run with shortened bearing and seal life and unplanned trips, while zone A/B supports the planned maintenance interval.`,
          'Agree the support stiffness and tuning with the machine vendor before the deck steel is ordered and include a site vibration acceptance test in the purchase order; retuning a foundation after start-up is paid for in lost production.', 'ISO 20816-1:2016 acceptance practice; ISO 14224:2016 reliability data and reliability-centred maintenance; whole-life costing per ISO 15686-5'),
      ],
      quantities: deckSup ? { steel_t: beamT } : { concrete_m3: m2 / 2400 },
    };
  },
  verify: [{
    // 10 t machine on undamped mounts of 10 kN/mm on a rigid support: fn = √(k/m)/(2π) = √(1e7/1e4)/(2π) = 5.0329 Hz.
    // Running at three times that (905.9 rpm) the undamped transmissibility is 1/(r² − 1) = 1/8 = 0.125.
    name: 'Undamped isolator on a rigid base: fn = √(k/m)/2π and TR = 1/(r² − 1) at r = 3', source: 'Den Hartog, Mechanical Vibrations, §2.12 (vibration isolation)',
    inputs: { mMach: 10, mRotor: 1, mount: 'iso', kIso: 10, zIso: 0, speed: 905.92, speedMin: 800, kUser: 1e9, mSup: 1, zSup: 0.2 },
    expect: [{ kpi: 'fnIso', value: 5.0329, tol: 1e-3 }, { kpi: 'TR', value: 0.125, tol: 0.005 }],
  }],
};
