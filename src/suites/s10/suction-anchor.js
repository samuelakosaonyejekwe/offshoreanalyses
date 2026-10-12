// Suction anchor (caisson) in clay: suction installation, horizontal and vertical holding capacity, inclined-load
// interaction at the padeye, optimal padeye depth and a lower-bound shell-buckling check on the suction.
import { suctionInstallClay, suctionPenetrationCurve, caissonVerticalCapacityClay, caissonHorizontalCapacityClay, inclinedInteraction, ringBucklingPressure, pumpFlow } from '../../engine/geotech/suction.js';
import { DEG, G, need, run, round, fmt, adv, check, worst, sev, linspace, steelCo2Text, vesselCo2Text } from './_common.js';
import { anchorLoadInputs, clayInputs, capacityFactorInputs, anchorLoad, LIMIT_ANCHOR } from './_anchor.js';

export default {
  id: 's10-suction-anchor',
  suite: 10,
  title: 'Suction anchor installation and capacity',
  summary: 'Checks that a suction caisson can be installed by under-pressure without plug failure or buckling and that it holds the inclined mooring load at the padeye.',
  covers: [10],
  fidelity: 'analytical',
  standards: ['DNV-RP-E303 (2017) Geotechnical design and installation of suction anchors in clay', 'API RP 2SK (3rd ed., 2005) §6.5, anchor safety factors', 'ISO 19901-4 (2016) §9', 'DNV-RP-C202 (2019) shell buckling (reference for detailed check)'],
  models: ['Suction installation in clay: wall friction, tip bearing, plug-heave limit (Houlsby & Byrne 2005)', 'Lateral bearing factor Np(z) (Murff & Hamilton 1993) with rigid-body rotation about the padeye', 'Reverse end bearing / plug pull-out vertical capacity (undrained shear strength method)', 'Inclined-load interaction (Supachawarote, Randolph & Gourvenec 2004)', 'Elastic ring buckling (Bresse–Bryan) as lower bound'],
  limits: `${LIMIT_ANCHOR} Clay only (no sand installation or drained tensile capacity here); the lateral solution ignores base shear and the buckling pressure is that of a long unstiffened ring — a conservative screen, not a DNV-RP-C202 code check. Padeye structural strength and torsion from line misalignment are not checked.`,
  publishes: ['mooring.anchorCapacity'],
  searchTerms: 'suction anchor caisson clay installation holding capacity padeye',
  inputs: [
    ...anchorLoadInputs({ H: 2750, V: 1300 }),
    { key: 'D', label: 'Caisson outer diameter', unit: 'm', type: 'number', default: 5, min: 1, max: 15, group: 'Caisson', role: 'geometry', design: true },
    { key: 'L', label: 'Skirt (embedded) length', unit: 'm', type: 'number', default: 15, min: 2, max: 60, group: 'Caisson', role: 'geometry', design: true, help: 'Aspect ratio L/D of 3–6 is usual in soft clay.' },
    { key: 'tw', label: 'Wall thickness', unit: 'mm', type: 'number', default: 35, min: 10, max: 120, group: 'Caisson', role: 'geometry', design: true, help: 'D/t of 120–200 is typical; thinner walls install more easily but buckle sooner.' },
    { key: 'steel', label: 'Steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Caisson', role: 'material' },
    { key: 'topMass', label: 'Top cap, padeye and stiffener mass', unit: 't', type: 'number', default: 25, min: 0, max: 500, group: 'Caisson', role: 'geometry' },
    { key: 'zPadeye', label: 'Padeye depth below mudline (0 = optimal)', unit: 'm', type: 'number', default: 0, min: 0, max: 60, group: 'Caisson', role: 'geometry',
      help: 'Attachment depth of the mooring line. At the optimal depth the caisson translates without rotating and the lateral capacity is greatest (about 2/3 of L in normally consolidated clay).' },
    { key: 'sealed', label: 'Top sealed in service', type: 'bool', default: true, group: 'Caisson', role: 'bc', help: 'Sealed: uplift mobilises reverse end bearing under the plug. Vented: wall friction inside and outside only.' },
    ...clayInputs(),
    { key: 'gammaSub', label: 'Submerged unit weight of soil', unit: 'kN/m³', type: 'number', default: 6, min: 2, max: 12, group: 'Soil', role: 'soil', bind: 'site.gammaSub' },
    { key: 'alphaInst', label: 'Adhesion factor during installation', unit: '-', type: 'number', default: 0.35, min: 0.05, max: 1, group: 'Soil', role: 'ic',
      help: 'Remoulded wall adhesion while penetrating ≈ 1 / sensitivity (0.25–0.5). Initial condition: the soil along the wall is disturbed by installation.' },
    { key: 'alphaSet', label: 'Adhesion factor after set-up', unit: '-', type: 'number', default: 0.65, min: 0.1, max: 1, group: 'Soil', role: 'soil', help: 'Long-term wall adhesion once the clay has reconsolidated (typically 0.6–0.7 after ~90 days).' },
    { key: 'NpConst', label: 'Lateral bearing factor Np (0 = Murff–Hamilton profile)', unit: '-', type: 'number', default: 0, min: 0, max: 14, group: 'Soil', role: 'soil', help: '0 uses Np rising from 2 at the mudline to 9 at depth. A constant 10–12 is a common quick estimate for a rough caisson.' },
    { key: 'penRate', label: 'Penetration rate', unit: 'mm/s', type: 'number', default: 1, min: 0.1, max: 10, group: 'Installation', role: 'operation' },
    ...capacityFactorInputs({ sf: 1.6 }),
    { key: 'sfLat', label: 'Required safety factor (0 = API table: lateral 1.6 / axial 2.0 intact permanent)', unit: '-', type: 'number', default: 0, min: 0, max: 5, group: 'Capacity factors', role: 'criteria',
      help: 'Overrides the factor applied to the inclined load. With 0 the lateral and axial API factors are interpolated on the load angle.' },
  ],
  compute(inp, ctx) {
    const ld = anchorLoad(inp);
    const mat = ctx.material(inp.steel);
    need(mat, `Steel '${inp.steel}' is not in the materials library.`);
    const D = inp.D, L = inp.L, t = inp.tw / 1e3, Di = D - 2 * t;
    need(t < D / 4, 'Wall thickness is not small against the diameter: check the units (mm).');
    const su0 = inp.su0 * 1e3, k = inp.suGrad * 1e3, gamma = inp.gammaSub * 1e3;
    need(su0 > 0 || k > 0, 'The clay has no strength: enter a mudline strength or a gradient.');
    const su = { su0, k };
    const steelMass = Math.PI * (D - t) * t * L * mat.rho + inp.topMass * 1e3;      // kg
    const Wsub = steelMass * G * (1 - 1025 / mat.rho);
    // --- installation
    const instPar = { Do: D, t, su, alphaO: inp.alphaInst, alphaI: inp.alphaInst, gamma, Vprime: Wsub };
    const curve = run(() => suctionPenetrationCurve({ soil: 'clay', L, n: 40, params: instPar }), 'Suction penetration curve');
    const full = run(() => suctionInstallClay({ ...instPar, h: L }), 'Suction installation');
    const pcr = ringBucklingPressure({ E: mat.E, nu: mat.nu ?? 0.3, t, D });
    const flow = pumpFlow({ Di, v: inp.penRate / 1e3 });
    // --- capacity
    const zP = inp.zPadeye > 0 ? Math.min(inp.zPadeye, L) : undefined;
    const hor = run(() => caissonHorizontalCapacityClay({ D, L, su0, k, NpConst: inp.NpConst > 0 ? inp.NpConst : undefined, zPadeye: zP }), 'Horizontal capacity');
    const ver = run(() => caissonVerticalCapacityClay({ Do: D, t, L, su, alphaO: inp.alphaSet, alphaI: inp.alphaSet, W: Wsub, gamma }), 'Vertical capacity');
    const Vmax = inp.sealed ? ver.sealed : ver.vented;
    const fac = inp.setUp * inp.cyclic;
    const Hd = hor.H * fac, Vd = Wsub + (Vmax - Wsub) * fac;                           // weight is not factored
    const inter = run(() => inclinedInteraction({ Hmax: Hd, Vmax: Vd, L, D, H: ld.H, V: ld.V, angle: ld.angle }), 'Inclined interaction');
    // API factors: 1.6 lateral, 2.0 axial (permanent intact); interpolate on the load angle unless overridden
    const tab = { permanent: { intact: [1.6, 2.0], damaged: [1.2, 1.5] }, mobile: { intact: [1.2, 1.5], damaged: [1.0, 1.2] } }[inp.service][inp.condition];
    const sf = inp.sfLat > 0 ? inp.sfLat : inp.sfRequired > 0 ? inp.sfRequired : tab[0] + (tab[1] - tab[0]) * ld.angle / (Math.PI / 2);
    const fs = inter.Tult / ld.T;
    const sDesign = full.s, sAllow = Math.min(full.sAllow, pcr);
    const checks = [
      check('Inclined holding capacity at the padeye', ld.T * sf / 1e3, inter.Tult / 1e3, 'kN', 'DNV-RP-E303 (2017) §4; API RP 2SK (3rd ed.) pile/suction anchor factors',
        'Guards against the caisson translating, rotating or pulling out under the inclined line load. Exceeding it means progressive displacement of the anchor and loss of line pretension.'),
      check('Suction against plug heave', sDesign / 1e3, full.sAllow / 1e3, 'kPa', 'Houlsby & Byrne (2005) plug-heave limit; DNV-RP-E303 (2017) §5',
        'Guards against the soil plug being sucked up inside the caisson before full penetration, which leaves the anchor short of its design depth.'),
      check('Suction against shell buckling (lower bound)', sDesign / 1e3, pcr / 1e3, 'kPa', 'Bresse–Bryan ring buckling; detailed check per DNV-RP-C202 (2019)',
        'Guards against the thin wall collapsing inwards under the installation under-pressure. The long-ring value ignores soil support and stiffeners, so a pass is safe and a fail calls for a code check.'),
    ];
    if (zP !== undefined) checks.push(check('Capacity lost to padeye position', hor.Hmax - hor.H, 0.1 * hor.Hmax, 'kN', 'DNV-RP-E303 (2017) §4.3 (optimal load attachment point)',
      'Compares the lateral capacity lost by attaching away from the optimal depth with a 10 % allowance. A misplaced padeye makes the caisson rotate and mobilise less soil.'));
    const g = worst(checks);
    const th = linspace(0, Math.PI / 2, 31);
    const env = th.map((a) => inclinedInteraction({ Hmax: Hd, Vmax: Vd, L, D, angle: a }).Tult);
    const steel_t = steelMass / 1e3;
    return {
      kpis: [
        { key: 'Hmax', label: 'Lateral capacity at optimal padeye', value: round(hor.Hmax / 1e3, 5), unit: 'kN', explain: 'Horizontal capacity when the line is attached at the centroid of soil resistance so the caisson translates without rotating.' },
        { key: 'Hpad', label: 'Lateral capacity at the chosen padeye', value: round(hor.H / 1e3, 5), unit: 'kN', explain: 'Horizontal capacity for the actual attachment depth; lower than the optimum because the caisson rotates.' },
        { key: 'zOpt', label: 'Optimal padeye depth', value: round(hor.zOptimal, 4), unit: 'm', explain: 'Depth of the centroid of lateral soil resistance; attach the line here (allowing for the embedded-line angle) for maximum capacity.' },
        { key: 'Vmax', label: 'Vertical pull-out capacity', value: round(Vmax / 1e3, 5), unit: 'kN', explain: 'Uplift capacity: with a sealed top the lesser of reverse end bearing and plug pull-out, vented the wall friction inside and outside, plus submerged weight.' },
        { key: 'Tult', label: 'Inclined capacity at the load angle', value: round(inter.Tult / 1e3, 5), unit: 'kN', explain: 'Ultimate line tension at the padeye for the actual inclination, from the H–V interaction envelope with set-up and cyclic factors.' },
        { key: 'fs', label: 'Factor of safety', value: round(fs, 4), unit: '-', limit: round(sf, 4), sense: 'min', explain: 'Inclined capacity over design load. API RP 2SK asks for 1.6 (lateral) to 2.0 (axial) for permanent intact moorings.' },
        { key: 'interaction', label: 'H–V interaction value', value: round(inter.f, 4), unit: '-', explain: 'Position of the unfactored load in the failure envelope: 1.0 is on the envelope, lower is inside.' },
        { key: 'selfWeight', label: 'Self-weight penetration', value: round(curve.selfWeight, 4), unit: 'm', explain: 'Depth the caisson sinks under its own submerged weight before pumping; it must be enough to seal the skirt tip (typically > 1 m).' },
        { key: 'suction', label: 'Suction required at full penetration', value: round(sDesign / 1e3, 5), unit: 'kPa', limit: round(sAllow / 1e3, 5), sense: 'max', explain: 'Under-pressure needed to push the skirt to its final depth against wall friction and tip resistance.' },
        { key: 'fsPlug', label: 'Margin against plug heave', value: round(Math.min(full.fsPlug, 99), 4), unit: '-', limit: 1.5, sense: 'min', explain: 'Allowable over required suction for plug stability at final depth; 1.5 is a customary minimum.' },
        { key: 'pBuckle', label: 'Ring buckling pressure (lower bound)', value: round(pcr / 1e3, 5), unit: 'kPa', explain: 'Elastic collapse pressure of a long unstiffened cylinder of this D/t. Ring stiffeners or soil support raise the real value.' },
        { key: 'steel', label: 'Caisson steel mass', value: round(steel_t, 4), unit: 't', explain: 'Skirt plus top cap and appurtenances; drives fabrication cost and the crane capacity of the installation vessel.' },
        { key: 'pump', label: 'Pump flow', value: round(flow.Q * 3600, 4), unit: 'm³/h', explain: 'Water to be evacuated to advance at the chosen penetration rate.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Suction vs penetration', x: { label: 'Required / allowable suction', unit: 'kPa' }, y: { label: 'Skirt tip depth', unit: 'm' },
          series: [{ name: 'Required suction', x: curve.s.map((v) => round(v / 1e3)), y: curve.h.map((v) => round(-v)) }, { name: 'Plug-heave limit', x: curve.sLimit.map((v) => round(v / 1e3)), y: curve.h.map((v) => round(-v)), dash: true },
            { name: 'Ring buckling (lower bound)', x: [pcr / 1e3, pcr / 1e3], y: [0, -L], dash: true }],
          explain: 'Under-pressure needed as the skirt goes down, against the two limits. Installation is feasible while the solid curve stays left of both dashed lines; where it starts at zero is the self-weight penetration.' },
        { type: 'line', title: 'H–V capacity envelope at the padeye', x: { label: 'Horizontal load', unit: 'kN' }, y: { label: 'Vertical load', unit: 'kN' },
          series: [{ name: 'Capacity envelope', x: th.map((a, i) => round(env[i] * Math.cos(a) / 1e3)), y: th.map((a, i) => round(env[i] * Math.sin(a) / 1e3)) },
            { name: 'Envelope / safety factor', x: th.map((a, i) => round(env[i] * Math.cos(a) / sf / 1e3)), y: th.map((a, i) => round(env[i] * Math.sin(a) / sf / 1e3)), dash: true }, { name: 'Design load', x: [round(ld.H / 1e3)], y: [round(ld.V / 1e3)] }],
          explain: 'Interaction diagram of horizontal and vertical capacity with the operating point. The design load must sit inside the dashed (factored) envelope; the nearly square shape shows that H and V interact weakly for long caissons.' },
      ],
      tables: [
        { title: 'Capacity components', columns: [{ key: 'q', label: 'Component' }, { key: 'v', label: 'Value', unit: 'kN' }],
          rows: [['Outside wall friction (set-up)', ver.frictionOut], ['Inside wall friction (set-up)', ver.frictionIn], ['Reverse end bearing', ver.endBearing], ['Sealed capacity: reverse bearing mechanism', ver.reverseBearing], ['Sealed capacity: plug pull-out mechanism', ver.plugPullout], ['Vented capacity', ver.vented], ['Submerged weight', Wsub], ['Lateral capacity, optimal padeye', hor.Hmax], ['Lateral capacity, chosen padeye', hor.H]].map(([q, v]) => ({ q, v: round(v / 1e3) })),
          explain: 'Break-down of the vertical and horizontal resistance before set-up and cyclic factors. The smaller of the two sealed mechanisms governs uplift.' },
        { title: 'Installation resistance at full penetration', columns: [{ key: 'q', label: 'Component' }, { key: 'v', label: 'Value', unit: 'kN' }],
          rows: [['Outside wall friction (remoulded)', full.outside], ['Inside wall friction (remoulded)', full.inside], ['Skirt tip bearing', full.tip], ['Total resistance', full.R], ['Submerged weight available', Wsub]].map(([q, v]) => ({ q, v: round(v / 1e3) })),
          explain: 'Soil resistance to penetration at final depth and the weight available to overcome it; the difference is supplied by suction over the caisson top area.' },
      ],
      notes: [`L/D = ${fmt(L / D)}, D/t = ${fmt(D / t)}; average lateral bearing factor Np = ${fmt(hor.NpAverage)}.`,
        `Load ${fmt(ld.T / 1e3)} kN at ${fmt(ld.angleDeg)}° at the padeye; required safety factor ${fmt(sf)} (${inp.sfLat > 0 || inp.sfRequired > 0 ? 'user value' : 'API lateral/axial factors interpolated on angle'}).`,
        inp.sealed ? `Sealed uplift governed by ${ver.reverseBearing <= ver.plugPullout ? 'reverse end bearing' : 'plug pull-out'}.` : 'Vented top: uplift relies on wall friction only.'],
      advice: [
        adv('engineering', sev(g.uc), `Factor of safety ${fmt(fs)} against ${fmt(sf)} required; installation needs ${fmt(sDesign / 1e3)} kPa of suction against ${fmt(full.sAllow / 1e3)} kPa (plug) and ${fmt(pcr / 1e3)} kPa (ring buckling). Governing: ${g.name} (${fmt(g.uc)}).`,
          g.uc > 1 ? (g.name.startsWith('Inclined') ? `Lengthen the skirt to about ${fmt(L * Math.sqrt(g.uc))} m (capacity grows roughly with L² in normally consolidated clay) or increase the diameter by ${fmt(100 * (g.uc - 1))} %.`
            : g.name.includes('buckling') ? `Increase the wall to about ${fmt(inp.tw * Math.cbrt(g.uc))} mm (ring buckling ∝ t³) or add ring stiffeners and verify to DNV-RP-C202.`
              : g.name.includes('plug') ? 'Reduce the required suction: add ballast for more self-weight penetration, or reduce the aspect ratio L/D.' : `Move the padeye to ${fmt(hor.zOptimal)} m below mudline.`)
            : `Place the padeye at ${fmt(hor.zOptimal)} m and confirm the installation adhesion factor ${inp.alphaInst} against the clay sensitivity; allow for a misalignment and tilt tolerance (typically ±7.5° and 5°).`,
          'DNV-RP-E303 (2017) §4–5; API RP 2SK (3rd ed., 2005) anchor safety factors (confirm the governing edition).'),
        adv('economic', 'info', `The caisson weighs ${fmt(steel_t)} t; capacity reserve is ${fmt(100 * (fs / sf - 1))} % over the required factor.`,
          fs > 1.15 * sf ? `A skirt about ${fmt(L * (1 - Math.sqrt(sf / fs)))} m shorter would still meet the factor, saving roughly ${fmt(steel_t * (1 - Math.sqrt(sf / fs)) * 0.8)} t per anchor and reducing lift weight; re-check installation suction after shortening.` : 'Capacity is matched to demand; savings would come from standardising one caisson size across all anchor locations to simplify fabrication and installation.',
          'Quantity-based Class 4 estimate per AACE RP 18R-97; installation by subsea crane vessel planned per DNV-ST-N001 (2021).',
          `${steelCo2Text(steel_t)} ${vesselCo2Text(1)} Suction anchors are removable by reversing the pumps, leaving no steel in the seabed at decommissioning.`),
      ],
      publish: { 'mooring.anchorCapacity': inter.Tult },
      quantities: { steel_t: round(steel_t), vessel_days: { suction_anchor_installation: 1 } },
    };
  },
  verify: [
    { // Uniform clay su = 20 kPa with a constant Np = 10: pure translation capacity H = Np·su·D·L = 10·20·5·15 = 15 000 kN,
      // and the centroid of a uniform reaction is at mid-depth, 7.5 m.
      name: 'Lateral capacity in uniform clay, H = Np·su·D·L',
      source: 'Randolph & Gourvenec (2011) Offshore Geotechnical Engineering, §9.3 (translational failure of a caisson)',
      inputs: { D: 5, L: 15, su0: 20, suGrad: 0, NpConst: 10, zPadeye: 0 },
      expect: [{ kpi: 'Hmax', value: 15000, tol: 0.005 }, { kpi: 'zOpt', value: 7.5, tol: 0.01 }] },
  ],
};
