// Installation and deployment: reel storage, lay over a chute or wheel, catenary lay tension, tensioner grip and
// crushing, dynamic tension bounds from vessel heave, touchdown bend radius and campaign duration.
import { axisymmetricResponse, equivalentBeamProperties, minimumBendRadius } from '../../engine/xsection/flexible.js';
import { cableLay, layDynamicTension } from '../../engine/lines/pipelay.js';
import { stackInputs, stackFrom, massSplit, capacityOf, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, linspace, sev, worst, UNBONDED, GENERIC, STD_17J, STD_17B, CONFIRM, vesselCo2, CO2_STEEL } from './_common.js';

export default {
  id: 's14-installation',
  suite: 14,
  title: 'Installation and deployment',
  summary: 'Checks reel and chute radii, lay tension, tensioner grip and crushing, dynamic tension and touchdown bending for a flexible-lay operation, and estimates the lay duration.',
  covers: [14],
  fidelity: 'analytical',
  standards: ['API RP 17B (5th ed.) §10 (installation)', 'API Spec 17J (4th ed.) §5.3.1, Table 6', 'DNV-ST-N001 (marine operations)', 'DNV-RP-N103 (lifting through the wave zone, dynamic amplification)'],
  models: ['Nonlinear catenary lay (inextensible, flat seabed)', 'Capstan / sidewall line load T/R on the chute', 'Friction-grip tensioner statics', 'Quasi-dynamic tension bounds for harmonic top motion', 'Minimum bend-radius compliance on reel, chute and touchdown'],
  limits: 'Static catenary without current or bending stiffness; dynamic tension is given as a lower bound (quasi-static plus inertia) and an upper bound (drag-locked elastic stretch), not a time-domain prediction. Crush and squeeze limits are manufacturer data. Reel packing, pull-in and seabed friction are not modelled. ' + UNBONDED + ' ' + GENERIC,
  publishes: ['install.vesselDays', 'flexible.layTension'],
  searchTerms: 'flexible pipe installation reel lay chute tensioner crushing',
  inputs: [
    ...stackInputs(),
    { key: 'fillState', label: 'Pipe contents during lay', type: 'select', default: 'seawater', options: [{ value: 'seawater', label: 'Flooded with sea water' }, { value: 'empty', label: 'Empty (air filled)' }], group: 'Deployment', role: 'ic' },
    { key: 'wOverride', label: 'Submerged weight (0 = from the layer stack)', unit: 'N/m', type: 'number', default: 0, min: 0, max: 20000, group: 'Deployment', role: 'load', help: 'Data-sheet submerged weight for the lay condition (flooded or empty as selected).' },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 1000, min: 5, max: 3500, group: 'Deployment', role: 'environment', bind: 'site.depth' },
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Deployment', role: 'environment', bind: 'site.rhoW' },
    { key: 'chuteHeight', label: 'Chute departure height above water', unit: 'm', type: 'number', default: 15, min: 0, max: 60, group: 'Deployment', role: 'geometry' },
    { key: 'H', label: 'Bottom (horizontal) lay tension', unit: 'kN', type: 'number', default: 25, min: 0.5, max: 2000, group: 'Deployment', role: 'bc', design: true, help: 'Horizontal tension held at touchdown; it sets the lay-back and the touchdown bend radius H/w.' },
    { key: 'chuteRadius', label: 'Chute / lay-wheel radius', unit: 'm', type: 'number', default: 3.5, min: 0.5, max: 20, group: 'Handling equipment', role: 'bc', design: true },
    { key: 'reelRadius', label: 'Reel hub radius', unit: 'm', type: 'number', default: 2.2, min: 0.5, max: 20, group: 'Handling equipment', role: 'ic', help: 'Initial storage curvature of the pipe = 1 / hub radius.' },
    { key: 'nReel', label: 'Number of reeling operations', type: 'int', default: 2, min: 1, max: 10, group: 'Handling equipment', role: 'operation', help: 'Factory spooling, transpooling and re-reeling each bend and straighten the pipe once.' },
    { key: 'crushLimit', label: 'Allowable line load on the chute', unit: 'kN/m', type: 'number', default: 400, min: 5, max: 5000, group: 'Handling equipment', role: 'criteria', help: 'Manufacturer limit for radial load per metre of pipe over a curved support.' },
    { key: 'nTracks', label: 'Tensioner tracks', type: 'int', default: 4, min: 2, max: 4, group: 'Handling equipment', role: 'geometry' },
    { key: 'trackLength', label: 'Tensioner contact length', unit: 'm', type: 'number', default: 8, min: 0.5, max: 20, group: 'Handling equipment', role: 'geometry', design: true },
    { key: 'padMu', label: 'Pad / outer sheath friction', unit: '-', type: 'number', default: 0.12, min: 0.02, max: 0.5, group: 'Handling equipment', role: 'material', help: 'Wet polymer on polyurethane pads: about 0.07–0.15.' },
    { key: 'squeezeLimit', label: 'Allowable squeeze per track', unit: 'kN/m', type: 'number', default: 400, min: 5, max: 5000, group: 'Handling equipment', role: 'criteria', help: 'Manufacturer crushing limit per track and per metre for the number of tracks used.' },
    { key: 'gripSF', label: 'Safety factor on tensioner grip', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'fTa', label: 'Tensile armour utilisation factor for installation', unit: '-', type: 'number', default: 0.67, min: 0.3, max: 1, group: 'Acceptance criteria', role: 'criteria', help: 'API 17J Table 6: 0.67 for functional loads, 0.85 with environmental loads. ' + CONFIRM },
    { key: 'storageFactor', label: 'Storage MBR / locking radius', unit: '-', type: 'number', default: 1.1, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'heaveAmp', label: 'Chute vertical motion amplitude', unit: 'm', type: 'number', default: 1.0, min: 0, max: 8, group: 'Vessel motion', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Heave plus pitch-induced motion at the chute in the limiting sea state.' },
    { key: 'period', label: 'Motion period', unit: 's', type: 'number', default: 9, min: 3, max: 25, group: 'Vessel motion', role: 'bc' },
    { key: 'dynBound', label: 'Dynamic tension used for the checks', type: 'select', default: 'lower', options: [{ value: 'lower', label: 'Lower bound (quasi-static + inertia)' }, { value: 'upper', label: 'Upper bound (drag-locked, elastic)' }], group: 'Vessel motion', role: 'control',
      help: 'The true amplitude lies between the bounds; use the upper bound for deep water and short periods unless a dynamic analysis is available.' },
    { key: 'length', label: 'Length to lay', unit: 'km', type: 'number', default: 5, min: 0.05, max: 200, group: 'Campaign', role: 'operation', bind: 'pipeline.length' },
    { key: 'layRate', label: 'Lay rate', unit: 'm/h', type: 'number', default: 350, min: 10, max: 2000, group: 'Campaign', role: 'operation' },
    { key: 'downtime', label: 'Weather downtime', unit: '%', type: 'number', default: 15, min: 0, max: 80, group: 'Campaign', role: 'operation', uncertain: { dist: 'normal', cov: 0.3 } },
    { key: 'mobDays', label: 'Mobilisation, load-out and transit', unit: 'days', type: 'number', default: 8, min: 0, max: 90, group: 'Campaign', role: 'operation' },
  ],
  compute(inp, ctx) {
    const S = stackFrom(inp, ctx);
    const flooded = inp.fillState === 'seawater';
    const beam = run(() => equivalentBeamProperties(S.stack, { weight: { rhoContent: flooded ? inp.rhoW : 0, rhoSeawater: inp.rhoW } }), 'Section properties');
    const w = inp.wOverride > 0 ? inp.wOverride : beam.mass.submergedOperating;
    need(w > 0, `The ${flooded ? 'flooded' : 'empty'} pipe has a submerged weight of ${fmt(w)} N/m: it floats and cannot be laid as a catenary. Lay it flooded or add ballast.`);
    const h = inp.depth + inp.chuteHeight;
    const mbr = run(() => minimumBendRadius(S.stack, { storageFactor: inp.storageFactor }), 'Minimum bend radius');
    const cl = run(() => cableLay({ w, h, H: inp.H * 1e3, chuteRadius: inp.chuteRadius, mbr: mbr.storage, swpLimit: inp.crushLimit * 1e3 }), 'Catenary lay');
    const mLine = (flooded ? beam.mass.massSeawaterFilled : beam.mass.massEmpty) + inp.rhoW * Math.PI / 4 * S.od ** 2;
    const dyn = run(() => layDynamicTension({ staticTension: cl.topTension, w, EA: beam.EA, m: mLine, suspendedLength: cl.suspendedLength, topAngleDeg: cl.departureAngleDeg, heaveAmp: inp.heaveAmp, period: inp.period }), 'Dynamic tension');
    const dT = inp.dynBound === 'upper' ? dyn.upper : dyn.lower;
    const Tdes = cl.topTension + dT;
    // allowable tension: wire stress linear in tension, pipe at the surface with no pressure difference
    let allowT = Infinity;
    if (S.idx.ta.length) {
      const a = axisymmetricResponse(S.stack, { tension: 0 }), b = axisymmetricResponse(S.stack, { tension: 1e5 });
      for (const i of S.idx.ta) {
        const s1 = (b.layers[i].wireStress - a.layers[i].wireStress) / 1e5;
        if (s1 > 0) allowT = Math.min(allowT, (inp.fTa * capacityOf(S.stack.layers[i]) - a.layers[i].wireStress) / s1);
      }
    }
    need(Number.isFinite(allowT), 'The stack has no tensile armour with strength data: the allowable lay tension cannot be found.');
    const squeeze = Tdes * inp.gripSF / (inp.padMu * inp.nTracks * inp.trackLength); // N/m per track
    const swp = Tdes / inp.chuteRadius;
    const checks = [
      chk('Lay tension vs allowable tension', Tdes / 1e3, allowT / 1e3, 'kN', `${STD_17J} Table 6 (factor ${inp.fTa})`, 'Guards against overstressing the tensile armour at the top of the suspended span, where static weight and vessel motion add.'),
      chk('Reel: storage MBR vs hub radius', mbr.storage, inp.reelRadius, 'm', `${STD_17J} §5.3.1 (storage MBR)`, 'Guards against over-bending the first wraps on the reel hub, which would lock or unlock the interlocked layers.'),
      chk('Chute: storage MBR vs chute radius', mbr.storage, inp.chuteRadius, 'm', `${STD_17J} §5.3.1; ${STD_17B} §10`, 'Guards against over-bending under tension where the pipe passes over the chute or lay wheel.'),
      chk('Chute line load (crushing)', swp / 1e3, inp.crushLimit, 'kN/m', `${STD_17B} §10 (crushing over curved supports); manufacturer limit`, 'Guards against crushing the carcass and pressure armour by the radial reaction T/R on the chute.'),
      chk('Tensioner squeeze per track', squeeze / 1e3, inp.squeezeLimit, 'kN/m', `${STD_17B} §10 (tensioner crushing); friction-grip statics with factor ${inp.gripSF}`, 'The squeeze needed to hold the lay tension by friction must not exceed the crushing limit of the pipe; otherwise the pipe ovalises or slips.'),
      chk('Touchdown: storage MBR vs catenary radius', mbr.storage, cl.tdpBendRadius, 'm', `${STD_17J} §5.3.1`, 'Guards against over-bending at touchdown when too little bottom tension is held.'),
      chk('No slack at the chute (upper-bound dynamic amplitude)', dyn.upper / 1e3, cl.topTension / 1e3, 'kN', `${STD_17B} §8, §10 — tension must stay positive through the motion cycle; DNV-RP-N103 slack criterion in intent`, 'Guards against the suspended pipe going slack and snatching in each wave cycle, which risks compression, bird-caging and looping. Compression at touchdown itself needs a dynamic analysis (global-response module).'),
    ];
    const gov = worst(checks);
    const layDays = inp.length * 1e3 / inp.layRate / 24 / Math.max(1 - inp.downtime / 100, 0.05);
    const days = inp.mobDays + layDays;
    const a = cl.H / w;
    const xs = linspace(0, cl.layBack, 60);
    const Hs = linspace(Math.max(0.2 * inp.H, 0.5), 4 * inp.H, 30);
    const win = Hs.map((Hk) => cableLay({ w, h, H: Hk * 1e3 }));
    const Hmin = mbr.storage * w, Hmax = Math.max(0, (allowT - dT) - w * h);
    const ms = massSplit(S.stack);
    const cycles = 2 * inp.nReel + 2;
    return {
      kpis: [
        kpi('topTension', 'Static top tension', cl.topTension / 1e3, 'kN', 'Tension at the chute from the catenary: bottom tension plus submerged weight times the departure height.'),
        kpi('designTension', 'Design lay tension incl. vessel motion', Tdes / 1e3, 'kN', `Static top tension plus the ${inp.dynBound}-bound dynamic amplitude; the load the tensioner and armour must take.`, allowT / 1e3, 'max'),
        kpi('dafLower', 'Dynamic amplification, lower bound', dyn.dafLower, '-', 'Peak over static tension if the line follows the motion quasi-statically.'),
        kpi('dafUpper', 'Dynamic amplification, upper bound', dyn.dafUpper, '-', 'Peak over static tension if drag locks the catenary and the motion is taken by elastic stretch. A large gap between the bounds calls for a dynamic analysis.'),
        kpi('departureAngle', 'Departure angle from horizontal', cl.departureAngleDeg, 'deg', 'Angle of the pipe leaving the chute; close to 90° for near-vertical lay.'),
        kpi('layBack', 'Lay-back to touchdown', cl.layBack, 'm', 'Horizontal distance from the chute to the touchdown point.'),
        kpi('tdpRadius', 'Touchdown bend radius', cl.tdpBendRadius, 'm', 'Radius of curvature at touchdown, H/w. Must stay above the storage MBR.', mbr.storage, 'min'),
        kpi('swp', 'Line load on the chute', swp / 1e3, 'kN/m', 'Radial reaction per metre of pipe on the chute, T/R.', inp.crushLimit, 'max'),
        kpi('squeeze', 'Tensioner squeeze needed per track', squeeze / 1e3, 'kN/m', 'Normal load per metre and per track needed to hold the design tension by friction with the grip safety factor.', inp.squeezeLimit, 'max'),
        kpi('maxUc', 'Governing installation utilisation', gov.uc, '-', `Highest utilisation of the installation checks (${gov.name}).`, 1, 'max'),
        kpi('Hmin', 'Minimum bottom tension (touchdown MBR)', Hmin / 1e3, 'kN', 'Smallest bottom tension keeping the touchdown radius at the storage MBR.'),
        kpi('Hmax', 'Maximum bottom tension (armour)', Hmax / 1e3, 'kN', 'Largest bottom tension before the design top tension reaches the allowable. Hmin to Hmax is the lay window.'),
        kpi('bendCycles', 'Bending cycles before service', cycles, '-', 'Bend-and-straighten cycles from reeling operations and the pass over the chute; each takes the sheaths to the storage strain.'),
        kpi('layDays', 'Vessel days for the campaign', days, 'days', 'Mobilisation plus laying time including weather downtime.'),
      ],
      checks,
      plots: [
        plot('line', 'Suspended catenary during lay', 'Horizontal distance from touchdown', 'm', 'Height above seabed', 'm', [{ name: 'Pipe', x: xs.map((x) => round(x)), y: xs.map((x) => round(a * (Math.cosh(x / a) - 1))) }],
          'Shape of the suspended pipe from touchdown (left) up to the chute. The curvature is greatest at touchdown.'),
        plot('line', 'Lay window: top tension and touchdown radius vs bottom tension', 'Bottom tension', 'kN', 'Top tension (kN) / radius (m)', '-', [
          { name: 'Static top tension, kN', x: Hs.map((v) => round(v)), y: win.map((q) => round(q.topTension / 1e3)) },
          { name: 'Touchdown radius, m', x: Hs.map((v) => round(v)), y: win.map((q) => round(q.tdpBendRadius)) },
          { name: 'Storage MBR, m', x: [round(Hs[0]), round(Hs[Hs.length - 1])], y: [round(mbr.storage), round(mbr.storage)], dash: true },
          { name: 'This case', x: [inp.H], y: [round(cl.tdpBendRadius)] },
        ], 'More bottom tension opens the touchdown radius but raises the top tension. The touchdown-radius line must stay above the dashed MBR line; the operating point is marked.'),
      ],
      tables: [tbl('Installation summary', [col('item', 'Item'), col('value', 'Value'), col('unit', 'Unit')], [
        { item: 'Submerged weight during lay', value: round(w), unit: 'N/m' }, { item: 'Suspended length', value: round(cl.suspendedLength), unit: 'm' },
        { item: 'Dynamic tension amplitude, lower / upper bound', value: `${fmt(dyn.lower / 1e3)} / ${fmt(dyn.upper / 1e3)}`, unit: 'kN' }, { item: 'Allowable tension', value: round(allowT / 1e3), unit: 'kN' },
        { item: 'Storage MBR', value: round(mbr.storage), unit: 'm' }, { item: 'Sheath strain on the reel hub', value: round(S.stack.outerRadius / inp.reelRadius * 100), unit: '%' },
        { item: 'Laying time incl. weather', value: round(layDays), unit: 'days' }, { item: 'Pipe mass on the reel(s)', value: round(beam.mass.massEmpty * inp.length), unit: 't' },
      ], 'Key figures for the installation procedure and the vessel specification.')],
      notes: [`${gov.name} governs (${fmt(gov.uc)}).`, `Lay window on bottom tension: ${fmt(Hmin / 1e3)} to ${fmt(Hmax / 1e3)} kN; operating at ${fmt(inp.H)} kN.`, `Dynamic amplification between ${fmt(dyn.dafLower)} and ${fmt(dyn.dafUpper)}; the ${inp.dynBound} bound is used.`, GENERIC],
      advice: [
        adv('engineering', sev(gov.uc), `Governing installation utilisation ${fmt(gov.uc)} (${gov.name}); design lay tension ${fmt(Tdes / 1e3)} kN against ${fmt(allowT / 1e3)} kN allowable.`,
          gov.uc > 1 ? `Bring the failing item inside its limit: bottom tension between ${fmt(Hmin / 1e3)} and ${fmt(Hmax / 1e3)} kN, chute and hub radius ≥ ${fmt(mbr.storage)} m, tensioner contact length ≥ ${fmt(inp.trackLength * squeeze / (inp.squeezeLimit * 1e3))} m.` : `Write the lay window (${fmt(Hmin / 1e3)}–${fmt(Hmax / 1e3)} kN bottom tension) and the limiting chute motion into the installation procedure; confirm crush limits with the manufacturer's crushing analysis.`,
          `${STD_17B} §10; ${STD_17J} Table 6; DNV-ST-N001 (2021)`),
        adv('economic', 'info', `The campaign needs about ${fmt(days)} vessel days (${fmt(layDays)} laying at ${fmt(inp.layRate)} m/h with ${fmt(inp.downtime)} % weather downtime).`,
          `Each 5 points of weather downtime removed saves ${fmt(inp.length * 1e3 / inp.layRate / 24 * (1 / Math.max(1 - inp.downtime / 100, 0.05) - 1 / Math.max(1 - (inp.downtime - 5) / 100, 0.05)))} days; a wider lay window (higher limiting sea state) is worth more than a faster lay rate for short lines.`,
          'Marine-operation planning per DNV-ST-N001; schedule estimate at AACE Class 3–4',
          `${vesselCo2(days)} The pipe itself embodies about ${fmt(ms.steel * inp.length * CO2_STEEL)} t CO₂e in ${fmt(ms.steel * inp.length)} t of steel (screening factor ${CO2_STEEL} t CO₂e/t).`),
      ],
      publish: { 'install.vesselDays': days, 'flexible.layTension': Tdes },
      quantities: { pipe_km: inp.length, steel_t: ms.steel * inp.length, polymer_t: ms.polymer * inp.length, vessel_days: { flexlay: days } },
    };
  },
  verify: [{
    // Inextensible catenary: T_top = H + w·h, touchdown radius a = H/w. w = 500 N/m, h = 1000 + 20 m, H = 50 kN →
    // T_top = 50 + 0.5·1020 = 560 kN; a = 100 m; departure angle acos(H/T) = acos(50/560) = 84.88°.
    name: 'Catenary lay: top tension, touchdown radius and departure angle',
    source: 'Irvine, Cable Structures, §1.2 (catenary); hand calculation',
    inputs: { wOverride: 500, depth: 1000, chuteHeight: 20, H: 50 },
    expect: [{ kpi: 'topTension', value: 560, tol: 0.001 }, { kpi: 'tdpRadius', value: 100, tol: 0.001 }, { kpi: 'departureAngle', value: 84.88, tol: 0.001 }],
  }],
};
