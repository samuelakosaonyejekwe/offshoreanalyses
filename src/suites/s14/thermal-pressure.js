// Thermal and pressure-induced response: steady temperature through the layers, shutdown cooldown (transient
// radial conduction), thermal and pressure elongation and twist, and the end force when expansion is restrained.
import { axisymmetricResponse, thermalProfile, pressureThermalElongation } from '../../engine/xsection/flexible.js';
import { radialTransient, timeToTemperature } from '../../engine/flow/thermalpipe.js';
import { stackInputs, stackFrom, loadInputs, loadsFrom, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, thin, sev, worst, LIMIT_AXI, GENERIC, STD_17J, STD_17B, DEG, BAR } from './_common.js';

const K0 = 273.15;

export default {
  id: 's14-thermal-pressure',
  suite: 14,
  title: 'Thermal and pressure-induced response',
  summary: 'Gives the layer temperatures, cooldown after shutdown, and the elongation, twist and restrained end force caused by temperature and pressure.',
  covers: [12],
  fidelity: 'numerical',
  standards: ['API RP 17B (5th ed.) §6 (thermal analysis), §8', 'API Spec 17J (4th ed.) §5.3.1 (polymer temperature limits)', 'API TR 17TR2 (PA11 ageing)'],
  models: ['Steady radial conduction through a composite cylinder', 'Transient heat conduction (finite-volume, implicit)', 'Thermal expansion of helical and polymer layers', 'Pressure-induced deformation of the layer stack', 'Temperature-dependent polymer modulus (linear derating)'],
  limits: LIMIT_AXI + ' Conduction is axisymmetric with constant properties; helical layers use an effective radial conductivity. Polymer softening with temperature is a user-set linear derating of modulus, not a viscoelastic model; coupled thermo-mechanical finite-element analysis is outside the module. ' + GENERIC,
  publishes: ['flexible.Uvalue', 'flexible.thermalElongation'],
  searchTerms: 'flexible pipe thermal analysis U-value cooldown thermal expansion',
  inputs: [
    ...stackInputs(),
    ...loadInputs({ Te: 0 }).filter((i) => i.key !== 'Te'),
    { key: 'Tin', label: 'Bore fluid temperature', unit: '°C', type: 'number', default: 80, min: -40, max: 150, group: 'Thermal boundaries', role: 'bc', design: true, uncertain: { dist: 'normal', cov: 0.05 } },
    { key: 'Tsea', label: 'Sea-water temperature', unit: '°C', type: 'number', default: 4, min: -2, max: 35, group: 'Thermal boundaries', role: 'bc', bind: 'site.Tbed' },
    { key: 'hIn', label: 'Internal film coefficient', unit: 'W/m²K', type: 'number', default: 1000, min: 1, max: 1e9, group: 'Thermal boundaries', role: 'bc', help: 'Flowing oil 300–2000, gas 100–1000; lower when stagnant.' },
    { key: 'hOut', label: 'External film coefficient', unit: 'W/m²K', type: 'number', default: 500, min: 1, max: 1e9, group: 'Thermal boundaries', role: 'bc', help: 'Sea water with current 500–2000; still water about 100–200.' },
    { key: 'Tinstall', label: 'Installation (stress-free) temperature', unit: '°C', type: 'number', default: 10, min: -20, max: 50, group: 'Thermal boundaries', role: 'ic', help: 'Temperature of all layers when the pipe was terminated.' },
    { key: 'length', label: 'Pipe length', unit: 'm', type: 'number', default: 1500, min: 1, max: 50000, group: 'Pipe construction', role: 'geometry', bind: 'pipeline.length' },
    { key: 'ends', label: 'End restraint', type: 'select', default: 'free', options: [{ value: 'free', label: 'Free to expand (riser, flowline with expansion)' }, { value: 'fixed', label: 'Fully restrained (buried or anchored)' }], group: 'Thermal boundaries', role: 'bc' },
    { key: 'eDerate', label: 'Polymer modulus loss per °C above 20 °C', unit: '%/°C', type: 'number', default: 1, min: 0, max: 3, group: 'Thermal boundaries', role: 'material', help: 'Linear derating of the polymer modulus with sheath temperature (floor at 20 % of the cold value). From the material data sheet.' },
    { key: 'tMaxSheath', label: 'Maximum pressure-sheath temperature', unit: '°C', type: 'number', default: 90, min: 40, max: 150, group: 'Acceptance criteria', role: 'criteria', help: 'Service limit of the sheath polymer: about 60 °C for HDPE, 65–90 °C for PA11 depending on water cut and pH (API 17TR2), 130 °C for PVDF.' },
    { key: 'holdHours', label: 'Shutdown duration (no-touch time required)', unit: 'h', type: 'number', default: 1.5, min: 0.1, max: 240, group: 'Shutdown cooldown', role: 'operation' },
    { key: 'Tmin', label: 'Minimum allowable fluid temperature', unit: '°C', type: 'number', default: 20, min: -20, max: 80, group: 'Shutdown cooldown', role: 'criteria', help: 'Hydrate or wax appearance temperature at shut-in pressure.' },
    { key: 'rhoFluid', label: 'Shut-in fluid density', unit: 'kg/m³', type: 'number', default: 700, min: 1, max: 1500, group: 'Shutdown cooldown', role: 'ic' },
    { key: 'cpFluid', label: 'Shut-in fluid heat capacity', unit: 'J/kgK', type: 'number', default: 2200, min: 500, max: 5000, group: 'Shutdown cooldown', role: 'ic' },
    { key: 'nCells', label: 'Radial cells per layer', type: 'int', default: 4, min: 1, max: 32, group: 'Numerics', role: 'control' },
    { key: 'nSteps', label: 'Time steps', type: 'int', default: 200, min: 20, max: 4000, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nCells', label: 'Radial cells per layer', levels: [2, 4, 8], qoi: 'Tcool' },
  meshTime: { key: 'nSteps', label: 'Time steps', levels: [100, 200, 400], qoi: 'Tcool' },
  compute(inp, ctx) {
    const cold = stackFrom(inp, ctx);
    const th = run(() => thermalProfile(cold.stack, { TinC: inp.Tin, ToutC: inp.Tsea, hIn: inp.hIn, hOut: inp.hOut }), 'Steady conduction');
    const iS = cold.idx.sheath >= 0 ? cold.idx.sheath : cold.stack.sealInner;
    const Tsheath = th.interfaceTempsC[iS]; // inner face of the pressure sheath: its hottest point
    const eFac = Math.max(0.2, 1 - inp.eDerate / 100 * Math.max(0, Tsheath - 20));
    const S = stackFrom(inp, ctx, { polymerEFactor: eFac });
    const dT = th.layerMeanTempsC.map((T) => T - inp.Tinstall);
    const loads = loadsFrom({ ...inp, Te: 0 });
    const el = run(() => pressureThermalElongation(S.stack, { pInt: loads.pInt, pExt: loads.pExt, deltaT: dT, annulus: loads.annulus }), 'Pressure and thermal elongation');
    const fix = run(() => axisymmetricResponse(S.stack, { pInt: loads.pInt, pExt: loads.pExt, annulus: loads.annulus, pAnnulus: loads.pAnnulus, deltaT: dT, axialStrain: 0, twist: 0 }), 'Restrained solution');
    const fixed = inp.ends === 'fixed';
    // shutdown: contents and wall cool from the steady flowing profile
    const tEnd = inp.holdHours * 3600;
    const layers = S.stack.layers.map((L) => {
      const solid = L.kind === 'helical' ? Math.min(1, L.nA / Math.cos(L.alpha) / (Math.PI * (L.rOut ** 2 - L.rIn ** 2))) : 1;
      return { t: L.t, k: L.k, rho: L.rho * solid, cp: L.cp ?? (L.kind === 'helical' ? 500 : 2000) };
    });
    const tr = run(() => radialTransient({ Di: 2 * S.stack.innerRadius, contents: { rho: inp.rhoFluid, cp: inp.cpFluid }, layers, hi: inp.hIn, ho: inp.hOut, Tamb: inp.Tsea + K0, T0: inp.Tin + K0, mode: 'cooldown', initial: 'steady', tEnd, dt: tEnd / inp.nSteps, cellsPerLayer: inp.nCells, theta: 1, maxOutput: 300 }), 'Cooldown simulation');
    const Tcool = tr.Tcontents[tr.Tcontents.length - 1] - K0;
    const tHit = inp.Tmin < inp.Tin && inp.Tmin > inp.Tsea ? timeToTemperature({ t: tr.t, T: tr.Tcontents, threshold: inp.Tmin + K0 }) : (inp.Tmin >= inp.Tin ? 0 : Infinity);
    const cooldownH = Math.min(tHit, tEnd) / 3600;
    const reached = Number.isFinite(tHit) && tHit <= tEnd;
    const elong = el.strainTotal * inp.length;
    const checks = [
      chk('Pressure sheath temperature', Tsheath, inp.tMaxSheath, '°C', `${STD_17J} §5.3.1; API 17TR2 for PA11 ageing`,
        'Guards against thermal ageing and softening of the pressure sheath. Above the limit the polymer hydrolyses or creeps faster and the sheath life falls steeply.'),
      chk('Shutdown: fluid temperature margin at the end of the no-touch time', inp.Tmin - inp.Tsea, Math.max(Tcool - inp.Tsea, 1e-6), '°C above sea temperature', 'Flow-assurance no-touch time (project philosophy); transient conduction per Carslaw & Jaeger',
        'Guards against hydrate or wax formation during a shutdown: the fluid must still be above its minimum allowable temperature when the no-touch time ends, otherwise the line must be depressurised or inhibited first.'),
    ];
    const gov = worst(checks);
    const [tt, TT, To] = thin([tr.t, tr.Tcontents, tr.Touter], 200);
    const radii = [S.stack.innerRadius, ...S.stack.layers.map((L) => L.rOut)];
    return {
      kpis: [
        kpi('heatLoss', 'Steady heat loss', th.heatLoss, 'W/m', 'Heat flowing from the bore fluid to the sea per metre of pipe in steady operation.'),
        kpi('Uvalue', 'Overall U-value on the bore area', th.UvalueInner, 'W/m²K', 'Overall heat-transfer coefficient referred to the bore surface; the number flow assurance uses. Lower means better insulated.'),
        kpi('sheathTemp', 'Pressure sheath temperature', Tsheath, '°C', 'Temperature at the inner face of the pressure sheath in steady flow, to compare with the polymer service limit.', inp.tMaxSheath, 'max'),
        kpi('outerTemp', 'Outer surface temperature', th.interfaceTempsC[th.interfaceTempsC.length - 1], '°C', 'Temperature of the outer sheath surface in steady flow.'),
        kpi('strainThermal', 'Thermal axial strain (free pipe)', el.strainThermal * 100, '%', 'Axial strain of the free pipe from the temperature rise of each layer above installation temperature.'),
        kpi('strainPressure', 'Pressure axial strain (free pipe)', el.strainPressure * 100, '%', 'Axial strain of the free, capped pipe from internal and external pressure; positive is lengthening.'),
        kpi('elongation', 'Free elongation over the length', elong, 'm', 'Combined pressure and thermal length change of the unrestrained pipe. It must be absorbed by slack, a bend or the end connection.'),
        kpi('rotation', 'Free end rotation over the length', (el.twistPressure + el.twistThermal) * inp.length / DEG, 'deg', 'Rotation of one end relative to the other caused by pressure and temperature when the pipe is free to twist.'),
        kpi('restrainedForce', 'Effective axial force if fully restrained', fix.tensionEffective / 1e3, 'kN', 'Effective force that develops when expansion and twist are prevented; negative is compression, which can buckle a flowline or bird-cage the armour.'),
        kpi('restrainedTorque', 'Torque if fully restrained', fix.torque / 1e3, 'kN·m', 'Torque reacted by the end fittings when twist is prevented.'),
        kpi('eFactor', 'Polymer modulus factor at temperature', eFac, '-', 'Fraction of the cold polymer modulus used for the hot pipe (linear derating on sheath temperature).'),
        kpi('Tcool', 'Fluid temperature at the end of the shutdown', Tcool, '°C', 'Bore fluid temperature after the shutdown duration, from the transient conduction solution.', inp.Tmin, 'min'),
        kpi('cooldownTime', 'Cooldown time to the minimum temperature', cooldownH, 'h', 'Time for the stagnant fluid to cool to the minimum allowable temperature; equal to the shutdown duration when the limit is not reached within it.'),
      ],
      checks,
      plots: [
        plot('line', 'Steady temperature through the wall', 'Radius', 'mm', 'Temperature', '°C', [{ name: 'Temperature', x: radii.map((x) => round(x * 1e3)), y: th.interfaceTempsC.map((T) => round(T)) },
          { name: 'Sheath limit', x: [round(radii[0] * 1e3), round(radii[radii.length - 1] * 1e3)], y: [inp.tMaxSheath, inp.tMaxSheath], dash: true }],
        'Temperature at each layer boundary from the bore outwards. The steep drops are across the polymer layers, which provide almost all of the insulation.'),
        plot('line', 'Shutdown cooldown', 'Time', 'h', 'Temperature', '°C', [
          { name: 'Bore fluid', x: tt.map((t) => round(t / 3600)), y: TT.map((T) => round(T - K0)) },
          { name: 'Outer surface', x: tt.map((t) => round(t / 3600)), y: To.map((T) => round(T - K0)) },
          { name: 'Minimum allowable', x: [0, round(inp.holdHours)], y: [inp.Tmin, inp.Tmin], dash: true },
        ], 'Fluid temperature after flow stops. The crossing with the dashed line is the cooldown time available before hydrate or wax risk.'),
      ],
      tables: [tbl('Layer temperatures', [col('layer', 'Layer'), col('Tin', 'Inner face', '°C'), col('Tmean', 'Mean', '°C'), col('dT', 'Rise above installation', '°C')],
        S.stack.layers.map((L, i) => ({ layer: L.name, Tin: round(th.interfaceTempsC[i]), Tmean: round(th.layerMeanTempsC[i]), dT: round(dT[i]) })),
        'Steady temperature of every layer and its rise above the installation temperature, which drives thermal strain.')],
      notes: [`Sheath at ${fmt(Tsheath)} °C against a limit of ${inp.tMaxSheath} °C; U = ${fmt(th.UvalueInner)} W/m²K.`,
        `${fixed ? 'Fully restrained: effective force ' + fmt(fix.tensionEffective / 1e3) + ' kN' : 'Free ends: elongation ' + fmt(elong) + ' m over ' + fmt(inp.length) + ' m'} (thermal ${fmt(el.strainThermal * 100)} %, pressure ${fmt(el.strainPressure * 100)} %).`,
        reached ? `Fluid reaches ${inp.Tmin} °C after ${fmt(cooldownH)} h.` : `Fluid stays above ${inp.Tmin} °C for the ${fmt(inp.holdHours)} h shutdown.`, GENERIC],
      advice: [
        adv('engineering', sev(gov.uc), `Sheath temperature ${fmt(Tsheath)} °C (limit ${inp.tMaxSheath} °C); ${fixed ? 'restrained force ' + fmt(fix.tensionEffective / 1e3) + ' kN' : 'free elongation ' + fmt(elong) + ' m'}; fluid at ${fmt(Tcool)} °C after ${fmt(inp.holdHours)} h.`,
          Tsheath > inp.tMaxSheath ? `Select a sheath polymer rated above ${fmt(Tsheath)} °C (PVDF) or cool the fluid by ${fmt(Tsheath - inp.tMaxSheath)} °C.` : (fixed && fix.tensionEffective < 0 ? `The restrained pipe goes into ${fmt(-fix.tensionEffective / 1e3)} kN of compression: check radial instability and upheaval, or allow the ${fmt(elong)} m of expansion.` : `Provide for ${fmt(elong)} m of end movement in the tie-in design and confirm the polymer ageing life at ${fmt(Tsheath)} °C.`),
          `${STD_17B} §6; ${STD_17J} §5.3.1; API 17TR2`),
        adv('economic', reached ? 'warn' : 'info', reached ? `The no-touch time is ${fmt(cooldownH)} h; a longer shutdown needs depressurisation or inhibitor.` : `The pipe keeps the fluid above ${inp.Tmin} °C for the whole ${fmt(inp.holdHours)} h shutdown without intervention.`,
          reached ? `An insulation layer that halves the U-value (now ${fmt(th.UvalueInner)} W/m²K) would roughly double the no-touch time; weigh its cost against the methanol and deferred production of each restart.` : `Steady heat loss is ${fmt(th.heatLoss * inp.length / 1e3)} kW over ${fmt(inp.length)} m; no additional insulation is justified by cooldown alone.`,
          'Flow-assurance operating-cost assessment; whole-life costing per ISO 15686-5',
          `Heat lost to the sea is ${fmt(th.heatLoss * inp.length * 8760 / 1e6)} MWh per year over this length; where the fluid is heated, better insulation saves that energy directly, at the price of more polymer foam in the wall.`),
      ],
      publish: { 'flexible.Uvalue': th.UvalueInner, 'flexible.thermalElongation': elong },
    };
  },
  verify: [{
    // Single PA11 tube, bore 200 mm, 20 mm wall, k = 0.33 W/mK (library), films removed: q = 2πk·ΔT/ln(ro/ri) =
    // 2π·0.33·(60 − 4)/ln(0.12/0.10) = 116.11/0.18232 = 636.9 W/m. Linear profile in ln r: mean of the face
    // temperatures = 32 °C, so with α = 1.1e-4 /K and installation at 4 °C the free thermal strain is 1.1e-4·28 = 0.308 %.
    name: 'Conduction through a single cylinder and free thermal strain',
    source: 'Incropera et al., Fundamentals of Heat and Mass Transfer, §3.3 (radial conduction); ε = α·ΔT',
    inputs: { boreD: 200, matSheath: 'pa11', Tin: 60, Tsea: 4, hIn: 1e9, hOut: 1e9, Tinstall: 4, pInt: 0, depth: 0, pAnn: 0,
      layers: [{ name: 'Tube', type: 'pressureSheath', t: 20, alphaDeg: 0, lay: 0, fill: 1, wireW: 0, E: 0, yield: 0, uts: 0, lock: 0 }] },
    expect: [{ kpi: 'heatLoss', value: 636.9, tol: 0.003 }, { kpi: 'strainThermal', value: 0.308, tol: 0.003 }],
  }],
};
