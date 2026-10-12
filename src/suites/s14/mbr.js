// Minimum bend radius and over-bending: locking radius of the interlocked layers, polymer sheath strain, storage /
// static / dynamic MBR and compliance of the radii actually applied on the reel and in service.
import { axisymmetricResponse, minimumBendRadius, bendingProperties } from '../../engine/xsection/flexible.js';
import { overbendCheck } from '../../engine/xsection/flexcollapse.js';
import { sheathStrainCheck } from '../../engine/xsection/flexfatigue.js';
import { stackInputs, stackFrom, loadInputs, loadsFrom, muInput, kpi, chk, adv, plot, bar, tbl, col, run, fmt, round, linspace, sev, worst, UNBONDED, GENERIC, STD_17J, STD_17B, CONFIRM, CO2_MGO, MGO_T_PER_DAY } from './_common.js';

export default {
  id: 's14-mbr',
  suite: 14,
  title: 'Minimum bend radius and over-bending',
  summary: 'Derives the storage, static and dynamic minimum bend radius from layer locking and sheath strain, and checks the radii applied on the reel and in service.',
  covers: [10],
  fidelity: 'analytical',
  standards: ['API Spec 17J (4th ed.) §5.3.1 (MBR factors, sheath strain limits)', 'API RP 17B (5th ed.) §8 (over-bending)'],
  models: ['Locking radius from gap closure of interlocked layers', 'Polymer sheath bending strain r/ρ', 'Minimum bend-radius factors (storage, static, dynamic)', 'Bilinear stick–slip moment–curvature relation'],
  limits: 'Geometric model: locking radius = outer radius of the interlocked layer / locking strain; ovalisation and layer interaction at the bend are ignored. MBR factors (1.1 / 1.0 / 1.5) and sheath strain limits are recalled defaults to be confirmed against the governing edition, and the locking strain must come from the manufacturer\'s profile geometry. ' + UNBONDED + ' ' + GENERIC,
  publishes: ['flexible.mbrStorage', 'flexible.mbrDynamic'],
  searchTerms: 'flexible pipe minimum bend radius locking radius sheath strain',
  inputs: [
    ...stackInputs(),
    { key: 'service', label: 'Service', type: 'select', default: 'dynamic', options: [{ value: 'dynamic', label: 'Dynamic (riser, jumper)' }, { value: 'static', label: 'Static (flowline)' }], group: 'Applied bending', role: 'operation' },
    { key: 'reelRadius', label: 'Reel hub / carousel inner radius', unit: 'm', type: 'number', default: 2.0, min: 0.3, max: 20, group: 'Applied bending', role: 'bc', design: true, help: 'Smallest radius the pipe is wound to for storage and transport (prescribed contact geometry).' },
    { key: 'opRadius', label: 'Smallest bend radius in service', unit: 'm', type: 'number', default: 6, min: 0.3, max: 5000, group: 'Applied bending', role: 'bc', bind: 'flexible.minBendRadius', design: true, help: 'From global analysis: touchdown, sag bend or inside the bend stiffener.' },
    { key: 'cycRadius', label: 'Bend radius of the cyclic curvature range', unit: 'm', type: 'number', default: 20, min: 0.3, max: 5000, group: 'Applied bending', role: 'bc', help: '1 / (peak-to-peak curvature range) in service, for the cyclic sheath strain.' },
    ...loadInputs({ Te: 300 }),
    muInput(),
    { key: 'storageFactor', label: 'Storage MBR / locking radius', unit: '-', type: 'number', default: 1.1, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'staticFactor', label: 'Static MBR / storage MBR', unit: '-', type: 'number', default: 1.0, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'dynamicFactor', label: 'Dynamic MBR / storage MBR', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'strainStatic', label: 'Allowable sheath strain, static and storage', unit: '%', type: 'number', default: 7.7, min: 1, max: 15, group: 'Acceptance criteria', role: 'criteria', help: 'PE and PA 7.7 %, PVDF 7.0 %. ' + CONFIRM },
    { key: 'strainDynamic', label: 'Allowable sheath strain, dynamic', unit: '%', type: 'number', default: 7.7, min: 0.5, max: 15, group: 'Acceptance criteria', role: 'criteria', help: 'PE and PA 7.7 %, PVDF 3.5 %. ' + CONFIRM },
    { key: 'strainRange', label: 'Allowable cyclic sheath strain range', unit: '%', type: 'number', default: 3.5, min: 0.1, max: 15, group: 'Acceptance criteria', role: 'criteria', help: 'Project value for cyclic strain of the sheath material; from the material qualification.' },
  ],
  compute(inp, ctx) {
    const S = stackFrom(inp, ctx, { strain: { static: inp.strainStatic / 100, dynamic: inp.strainDynamic / 100 } });
    const m = run(() => minimumBendRadius(S.stack, { storageFactor: inp.storageFactor, staticFactor: inp.staticFactor, dynamicFactor: inp.dynamicFactor }), 'Minimum bend radius');
    const dyn = inp.service === 'dynamic';
    const mbrService = dyn ? m.dynamic : m.static;
    const O = S.stack.layers[S.idx.outer];
    const sOp = sheathStrainCheck({ rOut: O.rOut, kappaMax: 1 / inp.opRadius, strainLimit: (dyn ? inp.strainDynamic : inp.strainStatic) / 100, kappaRange: 1 / inp.cycRadius, strainRangeLimit: inp.strainRange / 100 });
    const sReel = sheathStrainCheck({ rOut: O.rOut, kappaMax: 1 / inp.reelRadius, strainLimit: inp.strainStatic / 100 });
    const checks = [
      chk('Storage: MBR vs reel radius', m.storage, inp.reelRadius, 'm', `${STD_17J} §5.3.1 (storage MBR ≥ ${inp.storageFactor} × locking radius)`,
        'Guards against winding the pipe tighter than its storage bend radius, which closes the interlock gaps and can unlock the carcass or pressure armour.'),
      chk(`Service (${inp.service}): MBR vs smallest radius`, mbrService, inp.opRadius, 'm', `${STD_17J} §5.3.1 (${dyn ? 'dynamic ' + inp.dynamicFactor : 'static ' + inp.staticFactor} × storage MBR)`,
        'Guards against over-bending in service. Bending below the MBR overstrains the sheaths and loads the interlocked profiles directly.'),
      chk('Outer sheath strain in service', sOp.strain * 100, (dyn ? inp.strainDynamic : inp.strainStatic), '%', `${STD_17J} §5.3.1 (allowable sheath strain)`,
        'Guards against cracking or crazing of the polymer sheath at the outside of the bend.'),
      chk('Outer sheath cyclic strain range', sOp.strainRange * 100, inp.strainRange, '%', 'Material qualification (cyclic strain); API 17J §5.3.1',
        'Guards against fatigue of the polymer sheath from repeated bending; exceedance leads to sheath cracking and annulus flooding.'),
    ];
    let lockUc = 0, lockMargin = 0;
    if (m.lockingRadius) {
      const ob = overbendCheck({ curvature: 1 / Math.min(inp.reelRadius, inp.opRadius), lockingRadius: m.lockingRadius });
      lockUc = ob.utilisation; lockMargin = ob.marginRadius;
      checks.push(chk(`Locking of ${m.lockingLayer} at the tightest bend`, m.lockingRadius, Math.min(inp.reelRadius, inp.opRadius), 'm', `${STD_17B} (over-bending: unlocking of interlocked layers)`,
        'Guards against the interlocked profile reaching its geometric lock. Beyond that point further bending is carried by the profile hooks, which can unlock or yield.'));
    }
    const gov = worst(checks);
    const loads = loadsFrom(inp);
    const r = run(() => axisymmetricResponse(S.stack, loads), 'Axisymmetric solution');
    const bend = bendingProperties(S.stack, { response: r, mu: inp.mu });
    const momentAt = (rad) => Math.min(bend.EIstick / rad, bend.frictionMoment + bend.EIslip / rad);
    const radii = linspace(0.6 * Math.min(m.storage, inp.reelRadius), Math.max(3 * m.dynamic, 1.2 * inp.opRadius), 40);
    const limit = (dyn ? inp.strainDynamic : inp.strainStatic);
    return {
      kpis: [
        kpi('mbrStorage', 'Storage MBR', m.storage, 'm', 'Smallest radius allowed on a reel or carousel: the larger of the factored locking radius and the sheath-strain radius.'),
        kpi('mbrStatic', 'Static service MBR', m.static, 'm', 'Smallest radius allowed for a pipe that does not flex in service.'),
        kpi('mbrDynamic', 'Dynamic service MBR', m.dynamic, 'm', 'Smallest radius allowed where the pipe flexes continuously (risers, jumpers).'),
        kpi('lockingRadius', 'Locking radius', m.lockingRadius ?? 0, 'm', `Radius at which ${m.lockingLayer ?? 'no layer'} locks (interlock gaps closed). Zero when no layer carries a locking strain.`),
        kpi('ucStorage', 'Storage utilisation', m.storage / inp.reelRadius, '-', 'Storage MBR over the reel radius. Above 1.0 the reel is too small for this pipe.', 1, 'max'),
        kpi('ucService', 'Service utilisation', mbrService / inp.opRadius, '-', 'Service MBR over the smallest radius reached in service. Above 1.0 the pipe is over-bent.', 1, 'max'),
        kpi('sheathStrain', 'Outer sheath strain in service', sOp.strain * 100, '%', 'Bending strain at the outside of the outer sheath at the smallest service radius.', limit, 'max'),
        kpi('sheathStrainReel', 'Outer sheath strain on the reel', sReel.strain * 100, '%', 'Bending strain of the outer sheath when wound on the reel.', inp.strainStatic, 'max'),
        kpi('lockUc', 'Locking utilisation at the tightest bend', lockUc, '-', 'Locking radius over the tightest applied radius; 1.0 means the interlocked layer is fully locked.', 1, 'max'),
        kpi('momentMbr', 'Bending moment at the service MBR', momentAt(mbrService) / 1e3, 'kN·m', 'Moment needed to bend the pressurised, tensioned pipe to its service MBR — the load a bend restrictor or stiffener must be designed for.'),
      ],
      checks,
      plots: [
        bar('Bend radii', ['Locking radius', 'Storage MBR', 'Static MBR', 'Dynamic MBR', 'Reel radius', 'Smallest service radius'], 'Radius', 'm',
          [{ name: 'Radius', y: [round(m.lockingRadius ?? 0), round(m.storage), round(m.static), round(m.dynamic), inp.reelRadius, inp.opRadius] }],
          'Capacity radii (first four bars) against the radii actually applied (last two). The reel radius must exceed the storage MBR and the service radius the static or dynamic MBR.'),
        plot('line', 'Outer sheath strain and bending moment vs bend radius', 'Bend radius', 'm', 'Sheath strain', '%', [
          { name: 'Outer sheath strain', x: radii.map((x) => round(x)), y: radii.map((x) => round(O.rOut / x * 100)) },
          { name: 'Allowable strain', x: [round(radii[0]), round(radii[radii.length - 1])], y: [limit, limit], dash: true },
          { name: 'Bending moment / 10 (kN·m)', x: radii.map((x) => round(x)), y: radii.map((x) => round(momentAt(x) / 1e4)) },
        ], 'Sheath strain rises as the radius tightens and crosses the allowable at the sheath-strain radius. The moment curve (scaled by 1/10) shows how hard the pipe resists being bent that far.'),
      ],
      tables: [tbl('Bend-radius summary', [col('item', 'Item'), col('value', 'Radius', 'm'), col('basis', 'Basis')], [
        { item: 'Locking radius', value: round(m.lockingRadius ?? 0), basis: m.lockingLayer ? `${m.lockingLayer}: outer radius / locking strain` : 'no interlocked layer with a locking strain' },
        { item: 'Sheath-strain radius, static', value: round(m.sheathStaticRadius ?? 0), basis: `outer radius / ${inp.strainStatic} %` },
        { item: 'Sheath-strain radius, dynamic', value: round(m.sheathDynamicRadius ?? 0), basis: `outer radius / ${inp.strainDynamic} %` },
        { item: 'Storage MBR', value: round(m.storage), basis: `max(${inp.storageFactor} × locking, sheath strain)` },
        { item: 'Static MBR', value: round(m.static), basis: `${inp.staticFactor} × storage` },
        { item: 'Dynamic MBR', value: round(m.dynamic), basis: `${inp.dynamicFactor} × storage` },
      ], 'How each minimum bend radius is built up. Use the storage value for reels and chutes, the static value for flowlines and the dynamic value for risers.')],
      notes: [`${gov.name} governs (${fmt(gov.uc)}). Margin to locking at the tightest bend: ${fmt(lockMargin)} m.`, `Service MBR (${inp.service}) ${fmt(mbrService)} m; storage MBR ${fmt(m.storage)} m.`, GENERIC],
      advice: [
        adv('engineering', sev(gov.uc), `Storage MBR ${fmt(m.storage)} m against a reel radius of ${fmt(inp.reelRadius)} m; ${inp.service} MBR ${fmt(mbrService)} m against ${fmt(inp.opRadius)} m in service.`,
          gov.uc > 1 ? `Provide a radius of at least ${fmt(Math.max(m.storage, mbrService))} m where the check fails: larger reel hub, bend restrictor locking at ≥ ${fmt(mbrService)} m, or a longer bend stiffener.` : `Margin of ${fmt((1 - gov.uc) * 100)} %. Specify the bend restrictor or stiffener to lock no tighter than ${fmt(mbrService)} m and confirm the locking strain from the profile drawings.`,
          `${STD_17J} §5.3.1; ${STD_17B}. ${CONFIRM}`),
        adv('economic', 'info', `The storage MBR sets the smallest reel: hub diameter ≥ ${fmt(2 * m.storage)} m for this pipe.`,
          inp.reelRadius > 1.3 * m.storage ? `The reel hub is ${fmt((inp.reelRadius / m.storage - 1) * 100)} % larger than needed; a smaller standard reel holds more pipe per trip and may save a transit.` : 'The reel is close to the limit: make sure the first wraps are laid without overbending — a damaged inner wrap means cutting out pipe and a re-termination.',
          'Installation logistics planning (reel capacity and number of transits) per DNV-ST-N001; AACE Class 3 schedule basis',
          `Each avoided round trip of a reel vessel (about 2 days) saves roughly ${fmt(2 * MGO_T_PER_DAY)} t of fuel ≈ ${fmt(2 * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening: ${MGO_T_PER_DAY} t/day, ${CO2_MGO} t CO₂ per t fuel).`),
      ],
      publish: { 'flexible.mbrStorage': m.storage, 'flexible.mbrDynamic': m.dynamic },
    };
  },
  verify: [{
    // Default stack: pressure armour outer radius = 76.2 + 7 + 8 + 8 = 99.2 mm, locking strain 7 % → 1.4171 m (governs over
    // the carcass, 83.2/0.09 = 0.924 m). Storage = max(1.1·1.4171, outer sheath radius 118.2 mm / 7.7 % = 1.5351) = 1.5589 m;
    // static = 1.0 × storage; dynamic = 1.5 × storage = 2.3383 m.
    name: 'Locking radius and API 17J bend-radius factors',
    source: 'API Spec 17J §5.3.1 (factors as recalled); gap-closure geometry, hand calculation',
    inputs: { storageFactor: 1.1, staticFactor: 1.0, dynamicFactor: 1.5, strainStatic: 7.7, strainDynamic: 7.7 },
    expect: [{ kpi: 'lockingRadius', value: 1.4171, tol: 0.001 }, { kpi: 'mbrStorage', value: 1.5589, tol: 0.001 }, { kpi: 'mbrDynamic', value: 2.3383, tol: 0.001 }],
  }],
};
