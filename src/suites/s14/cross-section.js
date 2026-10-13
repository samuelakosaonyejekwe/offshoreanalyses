// Multilayer axisymmetric cross-section analysis of an unbonded flexible pipe: layer forces, wire and sheath
// stresses, interface contact pressures, API 17J Table 6 utilisations and the equivalent beam properties handed
// to global analysis.
import { axisymmetricResponse, layerUtilisations, equivalentBeamProperties } from '../../engine/xsection/flexible.js';
import { stackInputs, stackFrom, loadInputs, loadsFrom, factorInputs, factorsFrom, muInput, massSplit, kpi, chk, adv, bar, tbl, col, run, need, fmt, round, worst, sev,
  LIMIT_AXI, GENERIC, STD_17J, STD_17B, CO2_STEEL, BAR } from './_common.js';

export default {
  id: 's14-cross-section',
  suite: 14,
  title: 'Multilayer cross-section analysis',
  summary: 'Solves the layer-by-layer axisymmetric response of the pipe wall to tension, pressure, torque and temperature and checks every armour layer against its permissible utilisation.',
  covers: [2],
  fidelity: 'analytical',
  standards: ['API Spec 17J (4th ed.) §5.3, Table 6', 'API RP 17B (5th ed.) §8', 'ISO 13628-2'],
  models: ['Multilayer axisymmetric analysis with unilateral interlayer contact', 'Helical armour wire equilibrium (Féret & Bournazel 1987; Witz & Tan 1992)', 'Generalised Hooke\'s law for polymer layers', 'Effective tension and end-cap bookkeeping (Sparks)', 'Coupled axial–torsional stiffness'],
  limits: LIMIT_AXI + ' Permissible utilisation factors are recalled defaults to be confirmed against the governing edition of API 17J. ' + GENERIC,
  publishes: ['flexible.EA', 'flexible.EIslip', 'flexible.EIstick', 'flexible.GJ', 'flexible.submergedWeight', 'flexible.outerDiameter', 'flexible.massPerLength'],
  searchTerms: 'unbonded flexible pipe axisymmetric layer model tensile armour stress',
  inputs: [
    ...stackInputs(),
    ...loadInputs({ torque: true }),
    { key: 'twistFixed', label: 'Ends restrained against twist', type: 'bool', default: false, group: 'Loads and pressure boundaries', role: 'bc',
      help: 'Restrained: twist is zero and the reaction torque is reported (set applied torque to 0). Free: the applied torque is carried and the pipe twists.' },
    { key: 'deltaT', label: 'Uniform temperature rise from installation', unit: '°C', type: 'number', default: 0, min: -60, max: 150, group: 'Loads and pressure boundaries', role: 'ic',
      help: 'Temperature change of all layers relative to the stress-free (manufacturing) state.' },
    { key: 'rhoContent', label: 'Bore fluid density', unit: 'kg/m³', type: 'number', default: 700, min: 0, max: 2000, group: 'Internal fluid', role: 'load', help: 'Oil ≈ 700–850, gas ≈ 100–250, water 1025.' },
    muInput(),
    ...factorInputs(),
  ],
  compute(inp, ctx) {
    const S = stackFrom(inp, ctx);
    const loads = loadsFrom(inp);
    need(!(inp.twistFixed && Math.abs(loads.torque) > 0), 'Over-constrained torsion: the ends are restrained against twist and a torque is also applied. Set the torque to 0 (the reaction torque is reported) or free the ends.');
    const full = { ...loads, deltaT: inp.deltaT };
    if (inp.twistFixed) { delete full.torque; full.twist = 0; }
    const r = run(() => axisymmetricResponse(S.stack, full), 'Axisymmetric solution');
    const util = layerUtilisations(S.stack, r, { condition: inp.condition, factors: factorsFrom(inp) });
    const beam = run(() => equivalentBeamProperties(S.stack, { loads: full, mu: inp.mu, weight: { rhoContent: inp.rhoContent, rhoSeawater: inp.rhoW, annulusFlooded: inp.annulus === 'flooded' } }), 'Equivalent beam properties');
    const m = beam.mass, ms = massSplit(S.stack);

    const checks = util.layers.filter((u) => u.utilisation !== null).map((u) => chk(`${u.name}: wire stress (${inp.condition})`, u.stress / 1e6, u.factor * u.capacity / 1e6, 'MPa', `${STD_17J} (factor ${u.factor} × structural capacity)`,
      `Guards against rupture or yielding of the ${u.type === 'tensileArmour' ? 'tensile' : 'pressure'} armour wires: above the limit the layer no longer has the margin the specification requires against burst or tensile failure.`));
    S.stack.layers.forEach((L, i) => {
      if (L.kind === 'shell' && L.fy && (L.type === 'pressureSheath' || L.type === 'outerSheath')) {
        // deviatoric stress only: the confined sheath is close to hydrostatic compression, which does not yield it
        checks.push(chk(`${L.name}: von Mises stress`, r.layers[i].vonMises / 1e6, L.fy / 1e6, 'MPa', 'Polymer yield (material data sheet); API 17J §5.3.1 limits sheath strain',
          'Guards against yielding of the polymer sheath between its supports; a sheath that yields creeps into armour gaps and thins.'));
      }
    });
    need(checks.length > 0, 'No armour layer carries strength data (yield or tensile strength): enter them in the layer table.');
    const gov = worst(checks);
    const maxContact = Math.max(0, ...r.interfaces.map((f) => f.contactPressure));
    const open = r.interfaces.filter((f) => !f.closed).length;

    const names = S.stack.layers.map((L) => L.name);
    const rows = S.stack.layers.map((L, i) => {
      const a = r.layers[i], u = util.layers.find((x) => x.name === L.name);
      return { layer: L.name, type: L.type, R: round(L.R * 1e3), N: round(a.axialForce / 1e3), Mt: round(a.torque / 1e3), hoop: round(a.hoopForce / 1e3), u: round(a.radialDisp * 1e3),
        stress: round((a.wireStress ?? a.vonMises) / 1e6), sAx: a.sigmaAxial != null ? round(a.sigmaAxial / 1e6) : null, sHoop: a.sigmaHoop != null ? round(a.sigmaHoop / 1e6) : null,
        sRad: a.sigmaRadial != null ? round(a.sigmaRadial / 1e6) : null, shear: a.shear != null ? round(a.shear / 1e6) : null, cIn: round(a.contactIn / 1e6), cOut: round(a.contactOut / 1e6),
        uc: u && u.utilisation !== null ? round(u.utilisation, 3) : null };
    });
    const kpis = [
      kpi('maxUc', 'Governing layer utilisation', gov.uc, '-', `Highest stress utilisation of any layer (${gov.name}). Above 1.0 the layer exceeds its permissible stress for the selected load condition.`, 1, 'max'),
      kpi('wallForce', 'True wall force', r.wallForce / 1e3, 'kN', 'Sum of the axial forces in all layers: effective tension plus the pressure end-cap forces on the sealing layers. This is what the armour actually carries.'),
      kpi('axialStrain', 'Axial strain', r.axialStrain * 100, '%', 'Elongation per unit length under the applied loads. Flexible pipes typically stay below about 0.5 % in service.'),
      kpi('twist', 'Twist', r.twist / (Math.PI / 180), 'deg/m', 'Rotation per metre. A torque-balanced cross-wound armour pair keeps this close to zero under tension and pressure.'),
      kpi('reactionTorque', 'Torque carried', Math.abs(r.torque) < 1e-3 ? 0 : r.torque / 1e3, 'kN·m', 'Torque in the section: the applied value, or the reaction torque when the ends are restrained against twist. A torque below 1 mN·m is round-off of a torque-balanced section and is reported as 0.'),
      kpi('maxContact', 'Largest interface contact pressure', maxContact / 1e6, 'MPa', 'Highest pressure between adjacent layers. It sets interlayer friction, wear of the anti-wear tapes and creep of the sheath into armour gaps.'),
      kpi('openInterfaces', 'Open interfaces', open, '-', 'Number of layer interfaces that have lost contact. A gap next to a tensile armour layer means the wires are no longer radially supported.'),
      kpi('EA', 'Axial stiffness (torque-free ends)', beam.EA / 1e6, 'MN', 'Tangent axial stiffness about this load state; use it as EA in global riser analysis.'),
      kpi('GJ', 'Torsional stiffness (free to extend)', beam.GJ / 1e3, 'kN·m²', 'Tangent torsional stiffness about this load state.'),
      kpi('EIstick', 'Bending stiffness before slip', beam.EIstick / 1e3, 'kN·m²', 'Bending stiffness while friction holds the armour wires (very small curvatures only).'),
      kpi('EIslip', 'Bending stiffness after slip', beam.EIslip / 1e3, 'kN·m²', 'Bending stiffness once the wires slide — the value normally used for global analysis.'),
      kpi('wSub', 'Submerged weight, operating', m.submergedOperating, 'N/m', 'Weight in sea water with the bore full of the operating fluid; negative would mean the pipe floats.'),
      kpi('massAir', 'Mass in air, empty', m.massEmpty, 'kg/m', 'Dry mass per metre of the empty pipe, for reel and lifting weights.'),
      kpi('OD', 'Outer diameter', S.od * 1e3, 'mm', 'Outside diameter of the layer stack — the hydrodynamic and buoyancy diameter.'),
    ];
    const steelPerKm = ms.steel, polyPerKm = ms.polymer; // kg/m = t/km
    const armourRows = util.layers.filter((u) => u.utilisation !== null);
    const least = armourRows.reduce((a, c) => (c.utilisation < a.utilisation ? c : a), armourRows[0] || { utilisation: 0, name: '' });
    const advice = [
      adv('engineering', sev(gov.uc), `${gov.name} governs at a utilisation of ${fmt(gov.uc)} for the ${inp.condition} condition (${fmt(gov.demand)} against ${fmt(gov.capacity)} MPa).`,
        gov.uc > 1 ? `Raise the capacity of that layer by about ${fmt((gov.uc - 1) * 100)} % (thicker wire, higher fill or higher strength grade), or reduce the pressure/tension.` : `Reserve of ${fmt((1 - gov.uc) * 100)} % on the governing layer. Confirm the layer data against the manufacturer's design report and re-run for the abnormal and hydrotest conditions.`,
        `${STD_17J}; ${STD_17B}. Utilisation factors are recalled defaults — confirm against the governing edition.`),
      adv('economic', 'info', `The pipe wall contains about ${fmt(steelPerKm)} t of steel and ${fmt(polyPerKm)} t of polymer per kilometre.`,
        least.name && least.utilisation < 0.6 ? `${least.name} is only ${fmt(least.utilisation * 100)} % utilised: a ${fmt(Math.min(30, (1 - least.utilisation / 0.9) * 100))} % thinner wire would still meet the criterion — ask the manufacturer for an optimised structure before ordering.` : 'The armour is reasonably utilised; cost reduction is more likely from the design pressure or water-depth basis than from the layer sizes.',
        'Whole-life costing per ISO 15686-5; material quantity take-off at AACE Class 4 level',
        `${steelCo2PerKm(steelPerKm)}`),
      adv('engineering', open > 0 ? 'warn' : 'info', open > 0 ? `${open} interface(s) are open under this load case — layers are separating.` : 'All interfaces stay in contact under this load case.',
        open > 0 ? 'Check the radial-instability module: loss of contact next to the tensile armour is the precursor of bird-caging.' : 'Re-check with the lowest tension and bore pressure of the load-case matrix, where contact is most easily lost.', STD_17B),
    ];
    return {
      kpis, checks,
      plots: [
        bar('Stress by layer', names, 'Wire stress (metal) or von Mises (polymer)', 'MPa', [{ name: 'Stress', y: r.layers.map((a) => round((a.wireStress ?? a.vonMises) / 1e6)) }],
          'Stress in each layer from the bore outwards. The pressure armour carries hoop load, the tensile armour axial load; compare with the allowables in the check list.'),
        bar('Interface contact pressure', r.interfaces.map((_, j) => `${names[j]} / ${names[j + 1]}`), 'Contact pressure', 'MPa', [{ name: 'Contact', y: r.interfaces.map((f) => round(f.contactPressure / 1e6)) }],
          'Pressure between adjacent layers. Zero means the interface has opened; high values drive friction and wear.'),
        bar('Axial force share', names, 'Axial force', 'kN', [{ name: 'Axial force', y: r.layers.map((a) => round(a.axialForce / 1e3)) }],
          'How the true wall force is shared. Nearly all of it should sit in the tensile armour; a large share elsewhere indicates an unusual construction or input error.'),
      ],
      tables: [tbl('Layer results', [col('layer', 'Layer'), col('type', 'Type'), col('R', 'Mean radius', 'mm'), col('N', 'Axial force', 'kN'), col('Mt', 'Torque', 'kN·m'), col('hoop', 'Hoop force', 'kN/m'), col('u', 'Radial displacement', 'mm'),
        col('stress', 'Wire / von Mises stress', 'MPa'), col('sAx', 'Axial stress', 'MPa'), col('sHoop', 'Circumferential stress', 'MPa'), col('sRad', 'Radial stress', 'MPa'), col('shear', 'Shear stress', 'MPa'),
        col('cIn', 'Contact inner', 'MPa'), col('cOut', 'Contact outer', 'MPa'), col('uc', 'Utilisation', '-')], rows,
        'One row per layer, bore outwards. Metallic layers report the stress along the wire; polymer layers report axial, circumferential, radial and shear stress and the von Mises equivalent.')],
      notes: [`${gov.name} governs (${fmt(gov.uc)}).`, `End-cap force ${fmt(r.endCapForce / 1e3)} kN; external pressure ${fmt(loads.pExt / BAR)} bar; annulus ${inp.annulus}.`,
        `Contact iteration converged in ${r.iterations} pass(es); axial residual ${fmt(r.residuals.axial)} N.`, GENERIC],
      advice,
      publish: { 'flexible.EA': beam.EA, 'flexible.EIslip': beam.EIslip, 'flexible.EIstick': beam.EIstick, 'flexible.GJ': beam.GJ, 'flexible.submergedWeight': m.submergedOperating, 'flexible.outerDiameter': S.od, 'flexible.massPerLength': m.massEmpty },
    };
  },
  verify: [
    {
      // Sparks: wall force = Te + pi·Ai − pe·Ae on the sealing faces. Inner seal face r = 76.2 + 7 = 83.2 mm,
      // outer seal face r = 118.2 mm. pe = 1025·9.80665·1000 = 10.0518 MPa.
      // With the annulus at zero pressure: 500 000 + 20e6·π·0.0832² − 10.0518e6·π·0.1182² = 500 000 + 434 938 − 441 195
      // = 493.74 kN
      name: 'True wall force = effective tension + end-cap forces',
      source: 'Sparks, Fundamentals of Marine Riser Mechanics (2007) Eq. 2.5; hand calculation',
      inputs: { Te: 500, pInt: 200, depth: 1000, rhoW: 1025, pAnn: 0 },
      expect: [{ kpi: 'wallForce', value: 493.74, tol: 0.001 }],
    },
    {
      // Single polymer tube in pure tension is uniaxial: strain = T/(E·2πR·t), R = 105 mm, t = 10 mm, E = 1 GPa.
      // A = 2π·0.105·0.010 = 6.5973e-3 m²; strain = 100e3/(1e9·6.5973e-3) = 1.5158 %; EA = 6.5973 MN.
      name: 'Single polymer tube in tension: uniaxial Hooke\'s law',
      source: 'Hand calculation, ε = T/(EA)',
      inputs: { boreD: 200, Te: 100, pInt: 0, depth: 0, pAnn: 0, layers: [{ name: 'Tube', type: 'pressureSheath', t: 10, alphaDeg: 0, lay: 0, fill: 1, wireW: 0, E: 1, yield: 0, uts: 0, lock: 0 }] },
      expect: [{ kpi: 'axialStrain', value: 1.5158, tol: 0.002 }, { kpi: 'EA', value: 6.5973, tol: 0.002 }],
    },
  ],
};

function steelCo2PerKm(tPerKm) {
  return `Each kilometre embodies roughly ${fmt(tPerKm * CO2_STEEL)} t CO₂e in its steel layers alone (screening factor ${CO2_STEEL} t CO₂e per tonne — not a product declaration); every per cent of armour saved removes about ${fmt(tPerKm * CO2_STEEL / 100)} t CO₂e per kilometre.`;
}
