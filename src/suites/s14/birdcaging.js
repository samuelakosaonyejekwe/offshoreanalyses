// Radial instability and bird-caging of the tensile armour under axial compression, with the radial (elastic
// foundation) and lateral wire-buckling screens.
import { axisymmetricResponse } from '../../engine/xsection/flexible.js';
import { birdcagingCheck, wireFoundationBuckling, lateralBucklingScreen } from '../../engine/xsection/flexcollapse.js';
import { stackInputs, stackFrom, loadInputs, loadsFrom, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, linspace, sev, worst, LIMIT_AXI, GENERIC, STD_17B, CO2_MGO, MGO_T_PER_DAY, BAR, G } from './_common.js';

export default {
  id: 's14-birdcaging',
  suite: 14,
  title: 'Radial instability and bird-caging',
  summary: 'Screens the tensile armour for radial (bird-caging) and lateral wire buckling when the pipe goes into axial compression.',
  covers: [9],
  fidelity: 'screening',
  standards: ['API RP 17B (5th ed.) §8 (bird-caging, lateral buckling failure modes)', 'API Spec 17J (4th ed.) §5.3'],
  models: ['Radial equilibrium of compressed helical wires: outward pressure C·tan²α/(2πR²)', 'Hoop strength of high-strength tape and outer sheath', 'Wire on an elastic foundation, P = 2√(k·E·I) (Hetényi)', 'Euler screening of lateral wire buckling', 'Axisymmetric layer model with unilateral contact'],
  limits: LIMIT_AXI + ' Bird-caging is a strength check of the restraining tape and sheath, not a post-buckling analysis. Lateral wire buckling is an Euler screening with an assumed buckling length — not a published capacity model; cyclic curvature, friction decay and tape restraint are not represented. Use it to rank cases and flag those needing detailed analysis or test. ' + GENERIC,
  searchTerms: 'flexible pipe birdcaging lateral buckling tensile armour compression',
  inputs: [
    ...stackInputs(),
    ...loadInputs({ Te: -150, pInt: 1 }).map((i) => (i.key === 'Te' ? { ...i, label: 'Effective tension (negative = compression)', design: true, help: 'Lowest effective tension from global analysis, typically at touchdown in the near vessel position or during installation. Enter compression as a negative value.' } : i)),
    { key: 'tapeStrength', label: 'High-strength tape hoop strength', unit: 'MPa', type: 'number', default: 600, min: 0, max: 3000, group: 'Radial restraint', role: 'material', help: 'Tensile strength of the anti-bird-caging tape over the outer tensile armour (aramid or glass-reinforced). From supplier data.' },
    { key: 'tapeThk', label: 'High-strength tape total thickness', unit: 'mm', type: 'number', default: 1, min: 0, max: 10, group: 'Radial restraint', role: 'geometry', design: true, help: '0 = no high-strength tape; the outer sheath alone restrains the armour.' },
    { key: 'tapeE', label: 'High-strength tape modulus', unit: 'GPa', type: 'number', default: 20, min: 0.1, max: 150, group: 'Radial restraint', role: 'material' },
    { key: 'lengthFactor', label: 'Lateral buckling length / wire length per pitch', unit: '-', type: 'number', default: 1, min: 0.25, max: 4, group: 'Lateral buckling screen', role: 'control',
      help: 'Assumed buckling length of the wire in the cylinder surface as a fraction of one pitch of wire. This is an assumption of the screen (study 0.5–2); no published capacity model is implied.' },
    { key: 'sf', label: 'Safety factor on instability capacity', unit: '-', type: 'number', default: 1.5, min: 1, max: 4, group: 'Acceptance criteria', role: 'criteria', help: 'Project factor on the calculated critical compression; no normative value — agree it with the manufacturer.' },
  ],
  compute(inp, ctx) {
    const S = stackFrom(inp, ctx);
    need(S.idx.ta.length > 0, 'The stack has no tensile armour layer.');
    const loads = loadsFrom(inp);
    const r = run(() => axisymmetricResponse(S.stack, loads), 'Axisymmetric solution');
    const taL = S.idx.ta.map((i) => S.stack.layers[i]);
    const layerC = S.idx.ta.map((i) => Math.max(0, -r.layers[i].axialForce));
    const C = layerC.reduce((a, c) => a + c, 0);
    const outerTa = taL[taL.length - 1];
    const O = S.stack.layers[S.idx.outer];
    const bcArgs = (flooded, pExt) => ({ compression: C, layers: taL.map((L) => ({ R: L.R, alpha: L.alpha, weight: L.nA * L.E * Math.cos(L.alpha) ** 3 })),
      tape: inp.tapeThk > 0 && inp.tapeStrength > 0 ? { strength: inp.tapeStrength * 1e6, thickness: inp.tapeThk / 1e3, R: outerTa.rOut } : undefined,
      sheath: O.type === 'outerSheath' && O.fy ? { strength: O.fy, thickness: O.t, R: O.R } : undefined, pExt, pAnnulus: loads.pAnnulus, annulusFlooded: flooded });
    const flooded = inp.annulus === 'flooded';
    const bc = run(() => birdcagingCheck(bcArgs(flooded, loads.pExt)), 'Bird-caging check');
    const bcDry = birdcagingCheck(bcArgs(false, loads.pExt)), bcWet = birdcagingCheck(bcArgs(true, loads.pExt));
    need(bc.capacityPressure > 0, 'Nothing restrains the armour radially: add a high-strength tape or an outer sheath with a yield strength.');
    // elastic foundation: combined hoop stiffness of tape and sheath referred to the sheath
    const Et = (inp.tapeThk > 0 ? inp.tapeE * 1e9 * inp.tapeThk / 1e3 : 0) + O.E * O.t;
    const buck = taL.map((L, k) => {
      const n = L.n ?? Math.max(1, Math.round(L.fill * 2 * Math.PI * L.R * Math.cos(L.alpha) / (3 * L.t)));
      const b = L.wireWidth ?? 3 * L.t, t = L.wireThickness ?? L.t;
      const rad = wireFoundationBuckling({ E: L.E, In: b * t ** 3 / 12, restraintE: Et / O.t, restraintThickness: O.t, restraintR: O.R, R: L.R, alpha: L.alpha, n });
      const lat = lateralBucklingScreen({ E: L.E, It: t * b ** 3 / 12, R: L.R, alpha: L.alpha, n, lengthFactor: inp.lengthFactor, compression: layerC[k] });
      return { L, rad, lat, comp: layerC[k] };
    });
    const checks = [
      chk(`Bird-caging: outward armour pressure (${flooded ? 'flooded' : 'dry'} annulus)`, bc.demandPressure / 1e6, bc.capacityPressure / inp.sf / 1e6, 'MPa', `${STD_17B} (bird-caging); radial statics of compressed helices, safety factor ${inp.sf}`,
        'Guards against the compressed armour wires bursting outwards through the tape and sheath. A bird-cage is permanent: the armour loses its lay and the pipe must be replaced.'),
      ...buck.map((b) => chk(`${b.L.name}: radial wire buckling on elastic foundation`, b.comp / 1e3, b.rad.layerCompression / inp.sf / 1e3, 'kN', 'Hetényi, Beams on Elastic Foundation, Ch. VII (screening)',
        'Guards against individual wires buckling radially against the elastic restraint of tape and sheath, the local form of bird-caging.')),
    ];
    const gov = worst(checks);
    const latRatio = Math.max(...buck.map((b) => b.comp / b.lat.layerCompression));
    const open = r.interfaces.map((f, j) => ({ f, j })).filter(({ f, j }) => !f.closed && (S.idx.ta.includes(j) || S.idx.ta.includes(j + 1)));
    const gapMax = Math.max(0, ...r.interfaces.map((f) => f.gap));
    const Cs = linspace(0, Math.max(2 * C, 1.5 * bc.criticalCompression, 1e4), 21);
    const perUnit = bc.criticalCompression > 0 ? bc.capacityPressure / bc.criticalCompression : 0;
    const ds = linspace(0, Math.max(2 * inp.depth, 500), 21);
    const ucAt = (d, fl) => { const x = birdcagingCheck(bcArgs(fl, inp.rhoW * G * d)); return x.demandPressure / (x.capacityPressure / inp.sf); };
    const tapeReq = Math.max(0, (bc.demandPressure * inp.sf - bc.sheathCapacity - bc.pressureSupport) * outerTa.rOut / Math.max(inp.tapeStrength * 1e6, 1));
    return {
      kpis: [
        kpi('armourCompression', 'Compression carried by the tensile armour', C / 1e3, 'kN', 'True compressive force in the armour layers from the axisymmetric solution (effective compression plus the reverse end-cap effect of external pressure). Zero means the armour stays in tension.'),
        kpi('ucBirdcage', 'Bird-caging utilisation', bc.demandPressure / (bc.capacityPressure / inp.sf), '-', 'Outward pressure of the compressed wires over the factored hoop capacity of tape, sheath and (dry annulus only) external pressure. Above 1.0 bird-caging is predicted.', 1, 'max'),
        kpi('critCompression', 'Critical armour compression (unfactored)', bc.criticalCompression / 1e3, 'kN', 'Armour compression at which the restraining layers reach their hoop strength, before the safety factor.'),
        kpi('ucDry', 'Bird-caging utilisation, dry annulus', bcDry.demandPressure / (bcDry.capacityPressure / inp.sf), '-', 'Utilisation with an intact outer sheath, where external pressure clamps the armour.'),
        kpi('ucFlooded', 'Bird-caging utilisation, flooded annulus', bcWet.capacityPressure > 0 ? bcWet.demandPressure / (bcWet.capacityPressure / inp.sf) : 99, '-', 'Utilisation after the outer sheath is breached: external pressure no longer helps and only tape and sheath strength remain.'),
        kpi('ucRadial', 'Radial wire-buckling utilisation', Math.max(...buck.map((b) => b.comp / (b.rad.layerCompression / inp.sf))), '-', 'Layer compression over the factored elastic-foundation buckling load of its wires.', 1, 'max'),
        kpi('lateralRatio', 'Lateral buckling screening ratio', latRatio, '-', 'Layer compression over the friction-free Euler screening load for sideways wire buckling with the assumed buckling length. It is a ranking number, not a pass/fail criterion: values above 1 flag the case for a detailed lateral-buckling assessment or test, particularly with a flooded annulus and cyclic bending.'),
        kpi('wirePcr', 'Radial buckling load per wire, outer layer', buck[buck.length - 1].rad.Pcr / 1e3, 'kN', 'Critical compressive force of one wire of the outer armour on its elastic foundation, 2√(k·E·I).'),
        kpi('gapMax', 'Largest interface gap', gapMax * 1e3, 'mm', 'Largest radial gap that opens between layers under this load case; any gap next to the armour means lost radial support.'),
      ],
      checks,
      plots: [
        plot('line', 'Bird-caging utilisation vs armour compression', 'Armour compression', 'kN', 'Utilisation', '-', [
          { name: `${flooded ? 'Flooded' : 'Dry'} annulus`, x: Cs.map((c) => round(c / 1e3)), y: Cs.map((c) => round(c * perUnit / (bc.capacityPressure / inp.sf))) },
          { name: 'Limit', x: [0, round(Cs[Cs.length - 1] / 1e3)], y: [1, 1], dash: true },
          { name: 'This case', x: [round(C / 1e3)], y: [round(bc.demandPressure / (bc.capacityPressure / inp.sf))] },
        ], 'Utilisation grows in proportion to armour compression. The crossing with the limit line is the allowable compression for this annulus condition.'),
        plot('line', 'Bird-caging utilisation vs water depth', 'Water depth', 'm', 'Utilisation', '-', [
          { name: 'Dry annulus', x: ds.map((d) => round(d)), y: ds.map((d) => round(ucAt(d, false))) },
          { name: 'Flooded annulus', x: ds.map((d) => round(d)), y: ds.map((d) => round(ucAt(d, true))) },
          { name: 'Limit', x: [0, round(ds[ds.length - 1])], y: [1, 1], dash: true },
        ], 'For the same armour compression, depth helps a dry annulus (pressure clamps the armour) but not a flooded one. The gap between the lines is the value of an intact outer sheath.'),
      ],
      tables: [tbl('Wire buckling by layer', [col('layer', 'Layer'), col('comp', 'Layer compression', 'kN'), col('rad', 'Radial buckling capacity', 'kN'), col('half', 'Radial half-wavelength', 'mm'), col('lat', 'Lateral screen capacity', 'kN'), col('Lb', 'Assumed lateral buckling length', 'm')],
        buck.map((b) => ({ layer: b.L.name, comp: round(b.comp / 1e3), rad: round(b.rad.layerCompression / 1e3), half: round(b.rad.halfWavelength * 1e3), lat: round(b.lat.layerCompression / 1e3), Lb: round(b.lat.bucklingLength) })),
        'Unfactored layer capacities; the lateral screen is reported for ranking only and is not a pass/fail check. The radial value is the elastic-foundation bifurcation load; the lateral value is an Euler screen that depends directly on the assumed buckling length.')],
      notes: [C > 0 ? `The armour carries ${fmt(C / 1e3)} kN of compression; ${gov.name} governs at ${fmt(gov.uc)}.` : 'The armour stays in tension for this load case: no instability demand.',
        open.length ? `${open.length} interface(s) next to the tensile armour are open.` : 'All interfaces next to the tensile armour remain closed.',
        'The lateral buckling value is a screening number with an assumed buckling length; it is not a capacity prediction.', GENERIC],
      advice: [
        adv('engineering', sev(gov.uc), `Governing instability utilisation ${fmt(gov.uc)} (${gov.name}); armour compression ${fmt(C / 1e3)} kN against ${fmt(bc.criticalCompression / inp.sf / 1e3)} kN allowable for bird-caging.`,
          gov.uc > 1 ? (tapeReq > inp.tapeThk / 1e3 ? `Increase the high-strength tape to about ${fmt(tapeReq * 1e3)} mm, or remove the compression (more top tension, buoyancy near touchdown, slower lay).` : 'Remove the compression at source (more top tension, buoyancy near touchdown) and commission a detailed lateral-buckling assessment.') : 'Check the near-offset and installation cases for the lowest effective tension, and repeat with a flooded annulus — the critical combination for deep water.',
          `${STD_17B} (bird-caging and lateral buckling). The lateral screen is not a published capacity model.`),
        adv('engineering', latRatio > 1 ? 'warn' : 'info', `Lateral wire-buckling screening ratio ${fmt(latRatio)} (friction-free Euler column over ${fmt(buck[0].lat.bucklingLength)} m of wire).`,
          latRatio > 1 ? 'The screen is exceeded: this does not prove failure, because interlayer friction and tape restraint are ignored, but the case needs a dedicated lateral-buckling analysis or cyclic compression–bending test before acceptance — above all if the annulus can flood.' : 'The friction-free screen is not exceeded; no further lateral-buckling work is indicated for this case.',
          'API RP 17B §8 (lateral buckling of tensile armour); Østergaard, Lyckegaard & Andreasen (2012) Marine Structures 27 for the mechanisms not represented here'),
        adv('economic', 'info', `A flooded annulus changes the bird-caging utilisation from ${fmt(bcDry.demandPressure / (bcDry.capacityPressure / inp.sf))} to ${fmt(bcWet.capacityPressure > 0 ? bcWet.demandPressure / (bcWet.capacityPressure / inp.sf) : 99)}.`,
          'Avoiding compression operationally (weather limits on the near offset, minimum lay tension) costs waiting time; compare that downtime with the cost of extra tape in the pipe specification before ordering.',
          'Installation weather-window planning per DNV-ST-N001; value-of-information on annulus monitoring (DNV-RP-G101)',
          `A bird-caged section cannot be repaired in place: replacing it means manufacturing a new length and a campaign of roughly 10 vessel days ≈ ${fmt(10 * MGO_T_PER_DAY * CO2_MGO)} t CO₂ (screening: ${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂ per t).`),
      ],
    };
  },
  verify: [{
    // One armour layer, flooded annulus, no tape: capacity = sheath hoop strength σ_y·t/R, demand per unit compression =
    // tan²α/(2πR²). Stack: liner 5 mm on a 200 mm bore, armour 4 mm at 35° → R = 0.107 m; HDPE sheath 8 mm, R = 0.113 m,
    // σ_y = 22 MPa (library). Capacity = 22e6·0.008/0.113 = 1.5575 MPa; tan²35°/(2π·0.107²) = 0.49029/0.071936 = 6.8157 m⁻²
    // → C_cr = 1.5575e6/6.8157 = 228.5 kN.
    name: 'Critical armour compression for bird-caging, single layer restrained by the sheath',
    source: 'Radial statics of a compressed helix (Féret & Bournazel 1987); hand calculation',
    inputs: { boreD: 200, annulus: 'flooded', tapeThk: 0, matOuter: 'hdpe', layers: [
      { name: 'Liner', type: 'pressureSheath', t: 5, alphaDeg: 0, lay: 0, fill: 1, wireW: 0, E: 0, yield: 0, uts: 0, lock: 0 },
      { name: 'Armour', type: 'tensileArmour', t: 4, alphaDeg: 35, lay: 1, fill: 0.9, wireW: 12, E: 207, yield: 1100, uts: 1400, lock: 0 },
      { name: 'Outer sheath', type: 'outerSheath', t: 8, alphaDeg: 0, lay: 0, fill: 1, wireW: 0, E: 0, yield: 0, uts: 0, lock: 0 }] },
    expect: [{ kpi: 'critCompression', value: 228.5, tol: 0.003 }],
  }],
};
