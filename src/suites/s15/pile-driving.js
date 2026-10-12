// Pile driving: static resistance to driving, one-dimensional wave-equation drivability, driving stresses,
// refusal, self-weight penetration and pile-run screening, hammer selection and driving fatigue.
import { createProfile } from '../../engine/geotech/soilprofile.js';
import { GENERIC_HAMMERS, impactVelocity, srdProfile, drivabilityAnalysis, bearingGraph, selfWeightPenetration, pileSubmergedWeight } from '../../engine/geotech/driving.js';
import { snCurve, snConstants } from '../../engine/fatigue/sn.js';
import { G, need, run, round, fmt, kpi, check, adv, line, bar, tbl, worst, linspace, interp, steelCo2Text, vesselCo2Text, CONFIRM } from './_common.js';

const DEFAULT_LAYERS = [
  { type: 'sand', thickness: 8, gammaSub: 9.5, su: 0, suGrad: 0, phi: 32, Dr: 0.5, OCR: 1 },
  { type: 'clay', thickness: 22, gammaSub: 8, su: 40, suGrad: 2, phi: 0, Dr: 0, OCR: 2 },
  { type: 'sand', thickness: 60, gammaSub: 10, su: 0, suGrad: 0, phi: 37, Dr: 0.8, OCR: 1 },
];
const CASES = [{ value: 'coringLower', label: 'Coring, lower bound' }, { value: 'coringUpper', label: 'Coring, upper bound' }, { value: 'pluggedLower', label: 'Plugged, lower bound' }, { value: 'pluggedUpper', label: 'Plugged, upper bound' }];

export default {
  id: 's15-pile-driving', suite: 15, title: 'Pile driving and penetration',
  summary: 'Predicts the soil resistance to driving, blow counts, driving stresses and driving time of an open-ended pile under an impact hammer, flags refusal and pile run, and finds the smallest hammer that reaches the target penetration.',
  covers: [12], fidelity: 'numerical',
  standards: ['API RP 2A-WSD (2014) §6.10 and §12.5 (pile wall for driving, drivability)', 'API RP 2GEO (2011, add. 2014) §8 (axial pile capacity basis of the SRD)', 'DNV-ST-N001 (2021) Sec. 18 (offshore pile installation)', 'Smith (1960) Bearing capacity of piles by the wave equation, ASCE J. SMFD 86(SM4); Stevens, Wiltsie & Turton (1982) OTC 4205'],
  models: ['One-dimensional pile stress-wave propagation', 'Smith-type pile driving analysis (lumped masses, quake and damping)', 'Hammer–pile–soil interaction', 'Dynamic soil resistance (Smith damping)', 'Shaft friction and toe resistance (API unit values with Stevens SRD bounds)', 'Undrained soil strength profile', 'Self-weight penetration (vertical equilibrium)', 'Miner fatigue sum of driving stress cycles'],
  limits: 'Continuous driving without set-up, a rigid ram without cushion or anvil dynamics, no residual-stress analysis and Smith quake and damping at their classical values. The hammer list is generic (round-number energies, not manufacturer data) and the sand upper-bound SRD factors are the commonly quoted Stevens values; the clay upper-bound factor is 1.0 unless changed in the engine call. Vibratory driving, underwater noise, CPT-based SRD (Alm & Hamre) and pile tip buckling are outside the module.',
  publishes: ['install.pileBlows', 'install.pileDrivingHours', 'install.pileSRD'],
  searchTerms: 'pile drivability wave equation SRD blow count hammer selection refusal offshore',
  inputs: [
    { key: 'D', label: 'Pile outer diameter', unit: 'm', type: 'number', default: 2.134, min: 0.3, max: 12, group: 'Pile', role: 'geometry' },
    { key: 't', label: 'Wall thickness', unit: 'mm', type: 'number', default: 60, min: 8, max: 200, group: 'Pile', role: 'geometry', design: true },
    { key: 'Lpile', label: 'Pile length', unit: 'm', type: 'number', default: 75, min: 5, max: 200, group: 'Pile', role: 'geometry' },
    { key: 'pen', label: 'Target penetration', unit: 'm', type: 'number', default: 60, min: 2, max: 180, group: 'Pile', role: 'bc', design: true, help: 'Boundary condition of the installation: the penetration required by the in-place design.' },
    { key: 'mat', label: 'Pile steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Pile', role: 'material' },
    { key: 'layers', label: 'Soil layers from the mudline down', type: 'table', group: 'Soil', role: 'soil',
      columns: [{ key: 'type', label: 'Type (clay / sand)', type: 'text' }, { key: 'thickness', label: 'Thickness', unit: 'm', type: 'number' }, { key: 'gammaSub', label: 'γ′', unit: 'kN/m³', type: 'number' }, { key: 'su', label: 'su at top', unit: 'kPa', type: 'number' }, { key: 'suGrad', label: 'su gradient', unit: 'kPa/m', type: 'number' }, { key: 'phi', label: 'φ′', unit: 'deg', type: 'number' }, { key: 'Dr', label: 'Relative density', type: 'number' }, { key: 'OCR', label: 'OCR', type: 'number' }],
      default: DEFAULT_LAYERS, help: 'Clay rows use su, gradient and OCR; sand rows use φ′ and relative density (0–1). The default is a generic layered profile, not site data.' },
    { key: 'srdCase', label: 'Resistance-to-driving case', type: 'select', default: 'coringUpper', options: CASES, group: 'Soil', role: 'soil', help: 'Upper bounds govern refusal and hammer size; lower bounds govern pile run and minimum blow count.' },
    { key: 'hammer', label: 'Hammer', type: 'select', default: 'generic-1200', group: 'Hammer', role: 'bc', design: true, options: Object.entries(GENERIC_HAMMERS).map(([value, h]) => ({ value, label: h.label })), help: 'Generic hydraulic hammers for screening: replace with the data sheet of the hammer to be used.' },
    { key: 'eff', label: 'Hammer efficiency (energy at impact / rated)', unit: '-', type: 'number', default: 0.9, min: 0.3, max: 1, group: 'Hammer', role: 'bc', uncertain: { dist: 'normal', cov: 0.08 } },
    { key: 'mHammer', label: 'Hammer and sleeve mass resting on the pile', unit: 't', type: 'number', default: 250, min: 0, max: 2000, group: 'Hammer', role: 'load', help: 'For self-weight penetration and pile-run screening.' },
    { key: 'blowRate', label: 'Blow rate', unit: '1/min', type: 'number', default: 35, min: 5, max: 80, group: 'Hammer', role: 'operation' },
    { key: 'refusal', label: 'Refusal criterion', unit: 'blows/0.25 m', type: 'number', default: 250, min: 20, max: 2000, group: 'Criteria', role: 'criteria', help: `Contractual: commonly 250 blows per 0.25 m sustained (or 800 blows per 0.3 m once). ${CONFIRM}` },
    { key: 'sigFrac', label: 'Allowable driving stress / yield', unit: '-', type: 'number', default: 0.9, min: 0.5, max: 1, group: 'Criteria', role: 'criteria', help: `API RP 2A: dynamic stresses 80–90 % of yield. ${CONFIRM}` },
    { key: 'snName', label: 'S–N curve of the pile girth welds', type: 'select', default: 'D', group: 'Criteria', role: 'material', options: ['C1', 'D', 'E', 'F'].map((n) => ({ value: n, label: `DNV-RP-C203 curve ${n} (air)` })) },
    { key: 'nStation', label: 'Penetration stations analysed', unit: '-', type: 'int', default: 6, min: 3, max: 60, group: 'Numerics', role: 'control' },
    { key: 'nSeg', label: 'Pile segments in the wave-equation model', unit: '-', type: 'int', default: 24, min: 8, max: 200, group: 'Numerics', role: 'control', help: 'Blow counts change by about 5 % from 24 to 96 segments; use 48 or more for final work (each doubling takes about four times longer).' },
  ],
  mesh: { key: 'nSeg', label: 'Pile segments', levels: [12, 24, 48, 96], qoi: 'blowsEnd' },
  compute(inp, ctx) {
    // mm → m, kN/m³ → N/m³, kPa → Pa, t → kg.
    const steel = ctx.material(inp.mat); need(steel, `Steel '${inp.mat}' is not in the materials library.`);
    const D = inp.D, t = inp.t / 1e3;
    need(t < D / 2, 'Wall thickness must be smaller than the radius.');
    need(inp.pen <= inp.Lpile, 'The target penetration exceeds the pile length.');
    let top = 0;
    const layers = (inp.layers || []).map((l, i) => {
      const type = String(l.type || '').trim().toLowerCase(), th = Number(l.thickness);
      need(type === 'clay' || type === 'sand', `Soil layer ${i + 1}: type must be "clay" or "sand".`);
      need(th > 0, `Soil layer ${i + 1}: thickness must be positive.`);
      const base = { type, top, thickness: th, gammaSub: Number(l.gammaSub) * 1e3 };
      top += th;
      return type === 'clay' ? { ...base, suTop: Number(l.su) * 1e3, suGrad: Number(l.suGrad || 0) * 1e3, OCR: Math.max(Number(l.OCR) || 1, 1) } : { ...base, phiDeg: Number(l.phi), Dr: Math.min(Math.max(Number(l.Dr) || 0.5, 0.05), 1) };
    });
    need(layers.length > 0, 'Enter at least one soil layer.');
    need(top >= inp.pen, `The soil profile is ${fmt(top)} m deep but the target penetration is ${inp.pen} m: extend the last layer.`);
    const profile = run(() => createProfile(layers), 'Soil profile');
    const h0 = GENERIC_HAMMERS[inp.hammer], hammer = { ratedEnergy: h0.ratedEnergy, ramMass: h0.ramMass, efficiency: inp.eff };
    const pile = { length: inp.Lpile, D, t, E: steel.E || 2.1e11, rho: 7850, nSegments: inp.nSeg };
    const srd = run(() => srdProfile({ profile, D, t, Lmax: inp.pen, dz: 0.5 }), 'Resistance to driving');
    const depths = linspace(inp.pen / inp.nStation, inp.pen, inp.nStation);
    const dr = run(() => drivabilityAnalysis({ profile, hammer, pile, srdCase: inp.srdCase, depths, refusalBlowsPer025m: inp.refusal, dz: 0.5 }), 'Drivability analysis');
    const nD = dr.depth.length, capB = (v) => (Number.isFinite(v) ? Math.min(v, 9999) : 9999);
    const b025 = dr.blowsPer025m.map(capB), blowsEnd = b025[nD - 1], blowsMax = Math.max(...b025);
    const sigC = Math.max(...dr.maxCompressiveStress), sigT = Math.max(...dr.maxTensileStress), srdEnd = dr.srd[nD - 1];
    const As = Math.PI / 4 * (D * D - (D - 2 * t) ** 2), c = Math.sqrt(pile.E / pile.rho), Z = pile.rho * c * As, v0 = impactVelocity(hammer);
    const Wp = pileSubmergedWeight({ D, t, length: inp.Lpile }), swp = selfWeightPenetration({ z: Array.from(srd.z), srd: Array.from(srd.cases.coringLower.total), pileWeight: Wp, hammerWeight: inp.mHammer * 1e3 * G });
    const hours = dr.totalBlows / inp.blowRate / 60;
    // Driving fatigue: one stress cycle per blow at the range reported for each depth interval.
    const sn = snConstants(snCurve(inp.snName, 'air'));
    let dmg = 0; for (const b of dr.fatigueBlocks) { const S = b.stressRange; if (S > 0) dmg += b.blows / (S >= sn.sKnee || !sn.m2 ? sn.a1 * S ** -sn.m1 : sn.a2 * S ** -sn.m2); }
    // Hammer comparison at the final penetration (bearing graph at the final SRD with the shaft share of this profile).
    const shaftFrac = Math.min(Math.max(srd.cases[inp.srdCase].shaft[srd.z.length - 1] / srd.cases[inp.srdCase].total[srd.z.length - 1], 0), 1);
    const keys = Object.keys(GENERIC_HAMMERS), iSel = keys.indexOf(inp.hammer);
    const cmpKeys = keys.filter((k, i) => Math.abs(i - iSel) <= 1);
    const cmp = cmpKeys.map((k) => { const hh = GENERIC_HAMMERS[k]; const bg = run(() => bearingGraph({ hammer: { ratedEnergy: hh.ratedEnergy, ramMass: hh.ramMass, efficiency: inp.eff }, pile, srd: [srdEnd], shaftFraction: shaftFrac, penetration: inp.pen }), 'Bearing graph'); return { key: k, E: hh.ratedEnergy / 1e3, b: capB(bg.blowsPer025m[0]), sig: bg.maxCompressiveStress[0] }; });
    const okH = cmp.filter((q) => q.b <= inp.refusal && q.sig <= inp.sigFrac * steel.fy).sort((p, q) => p.E - q.E)[0];
    const Ereq = okH ? okH.E : cmp[cmp.length - 1].E;
    const checks = [
      check('Blow count against refusal', blowsMax, inp.refusal, 'blows/0.25 m', 'API RP 2A-WSD (2014) §12.5.6 (refusal definition); project specification', 'Above the refusal criterion driving must stop: the pile cannot be brought to its design penetration with this hammer without relief drilling or a larger hammer.'),
      check('Compressive driving stress', sigC / 1e6, inp.sigFrac * steel.fy / 1e6, 'MPa', 'API RP 2A-WSD (2014) §6.10 (dynamic stresses during driving 80–90 % of yield)', 'Peak stress wave in the pile wall. Overstress buckles or yields the pile head or tip and accumulates fatigue damage.'),
      check('Driving fatigue of girth welds', dmg, 0.3, 'Miner sum (budget 0.3)', 'DNV-RP-C203 (2021) §2.4; installation share of the fatigue budget by project practice', 'Every blow is a stress cycle. A hard drive can consume a significant share of the fatigue life of the pile welds before the structure is in service.'),
    ];
    if (swp.runRisk) checks.push(check('Pile run: resistance falls below pile + hammer weight', 1.2, 1, '- (run risk present)', 'Vertical equilibrium on the lower-bound SRD (pile-run screening)', 'A strong layer over a weaker one: when the tip breaks through, the pile and hammer can drop metres without a blow, damaging the hammer and its rigging.'));
    const gov = worst(checks), ok = gov.uc <= 1;
    const steel_t = As * 7850 * inp.Lpile / 1e3;
    const zS = Array.from(srd.z).map((z) => -round(z));
    return {
      kpis: [
        kpi('srdEnd', 'Resistance to driving at the target penetration', srdEnd / 1e6, 'MN', 'Static soil resistance the hammer must overcome at final penetration for the selected case.'),
        kpi('blowsEnd', 'Blow count at the target penetration', blowsEnd, 'blows/0.25 m', 'From the wave-equation analysis at final penetration. Compare with the refusal criterion; values capped at 9999.', inp.refusal, 'max'),
        kpi('blowsMax', 'Maximum blow count during the drive', blowsMax, 'blows/0.25 m', 'Largest blow count at any analysed station: a dense layer above the target can govern.', inp.refusal, 'max'),
        kpi('totalBlows', 'Total blows to the target (or to refusal)', dr.totalBlows, '-', 'Blow count integrated over penetration; with the blow rate it gives the net driving time and the fatigue cycles.'),
        kpi('hours', 'Net driving time', hours, 'h', 'Total blows over the blow rate, without stops for add-ons, hammer cooling or weather.'),
        kpi('sigC', 'Maximum compressive driving stress', sigC / 1e6, 'MPa', 'Peak compressive stress anywhere in the pile during the drive.', inp.sigFrac * steel.fy / 1e6, 'max'),
        kpi('sigT', 'Maximum tensile driving stress', sigT / 1e6, 'MPa', 'Reflected tension wave; highest in easy driving. Relevant for welds and connectors.'),
        kpi('enthru', 'Energy transferred to the pile at final penetration', dr.enthru[nD - 1] / 1e3, 'kJ', 'Maximum energy passing the pile head (ENTHRU). For a rigid ram without cushion on a long pile it is close to the ram kinetic energy.'),
        kpi('v0', 'Ram impact velocity', v0, 'm/s', '√(2·η·E/m) from the rated energy, efficiency and ram mass.'),
        kpi('sigImp', 'Analytic peak stress of a rigid ram on an elastic rod', pile.rho * c * v0 / 1e6, 'MPa', 'ρ·c·v0: the theoretical stress at the instant of impact; the lumped model overshoots it by up to about 20 % without a cushion.'),
        kpi('Z', 'Pile impedance ρ·c·A', Z / 1e6, 'MN·s/m', 'Matches hammer to pile: heavy rams suit high-impedance piles.'),
        kpi('swp', 'Self-weight penetration (pile + hammer)', swp.swpPileHammer ?? inp.pen, 'm', 'Depth at which the lower-bound resistance first equals the submerged pile weight plus the hammer; the full target if the weight is never resisted.'),
        kpi('Ereq', 'Smallest rated hammer energy that reaches the target', Ereq, 'kJ', 'Smallest hammer among those compared whose blow count and driving stress at final penetration are acceptable.'),
        kpi('damage', 'Driving fatigue damage', dmg, '-', 'Miner sum over all blows with the stress range at each depth interval (an upper estimate of the range).', 0.3, 'max'),
      ],
      checks,
      plots: [
        line('Resistance to driving against depth', 'SRD', 'MN', 'Penetration', 'm', CASES.map((q) => ({ name: q.label, x: Array.from(srd.cases[q.value].total).map((v) => round(v / 1e6)), y: zS, dash: q.value !== inp.srdCase })),
          'The four Stevens cases bound the resistance. The gap between coring and plugged cases shows how much the outcome depends on plug behaviour; the solid line is the case analysed.'),
        line('Blow count and driving stress against depth', 'Blow count (blows/0.25 m) · stress (MPa)', '', 'Penetration', 'm', [{ name: 'Blow count', x: b025.map((v) => round(v)), y: dr.depth.map((z) => -round(z)) }, { name: 'Compressive stress (MPa)', x: dr.maxCompressiveStress.map((v) => round(v / 1e6)), y: dr.depth.map((z) => -round(z)) }, { name: 'Refusal', x: [inp.refusal, inp.refusal], y: [0, -inp.pen], dash: true }],
          'Blow count rises with resistance; where it reaches the refusal line driving stops. Stress stays near the impact value and peaks in hard driving.'),
        bar('Hammer comparison at the target penetration', cmp.map((q) => `${q.E} kJ`), 'Blow count', 'blows/0.25 m', [{ name: 'Blow count', y: cmp.map((q) => round(q.b)) }, { name: 'Refusal', y: cmp.map(() => inp.refusal) }],
          'Blow count at final penetration for neighbouring hammer sizes. The smallest hammer below the refusal bar, with acceptable stress, is the economic choice.'),
      ],
      tables: [tbl('Drivability by penetration', [['z', 'Penetration', 'm'], ['srd', 'SRD', 'MN'], ['b', 'Blow count', 'blows/0.25 m'], ['sc', 'Compressive stress', 'MPa'], ['st', 'Tensile stress', 'MPa'], ['e', 'ENTHRU', 'kJ'], ['r', 'Refusal']],
        dr.depth.map((z, i) => ({ z: round(z), srd: round(dr.srd[i] / 1e6), b: round(b025[i]), sc: round(dr.maxCompressiveStress[i] / 1e6), st: round(dr.maxTensileStress[i] / 1e6), e: round(dr.enthru[i] / 1e3), r: dr.refusal[i] ? 'yes' : 'no' })),
        'One wave-equation blow analysed at each penetration station: the predicted driving record to compare with the pile-driving monitoring on site.')],
      notes: [`${h0.label} at ${inp.eff} efficiency: ${fmt(blowsEnd)} blows/0.25 m at ${inp.pen} m against SRD ${fmt(srdEnd / 1e6)} MN (${CASES.find((q) => q.value === inp.srdCase).label}).`, dr.refusalDepth != null ? `Refusal is predicted at ${fmt(dr.refusalDepth)} m.` : 'No refusal predicted.', swp.runRisk ? 'Pile-run risk: the lower-bound resistance drops below the pile + hammer weight below the self-weight penetration.' : `Self-weight penetration about ${fmt(swp.swpPileHammer ?? inp.pen)} m with the hammer on.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `The pile reaches ${inp.pen} m with ${fmt(blowsEnd)} blows/0.25 m, a peak driving stress of ${fmt(sigC / 1e6)} MPa and about ${fmt(dr.totalBlows)} blows (${fmt(hours)} h net).`,
          `Run the lower-bound case for pile run and the plugged upper bound for refusal, then confirm with the actual hammer data and ${Math.max(48, inp.nSeg)} or more segments. Monitor the drive (PDA) and compare with this table.`, 'API RP 2A-WSD (2014) §12.5; DNV-ST-N001 (2021) Sec. 18')
          : adv('engineering', 'fail', `"${gov.name}" is not met (${fmt(gov.uc)}): ${fmt(blowsMax)} blows/0.25 m maximum, ${fmt(sigC / 1e6)} MPa, fatigue damage ${fmt(dmg)}.`,
            gov.name.startsWith('Blow') ? `Use a hammer of at least ${fmt(Ereq)} kJ rated energy${okH ? '' : ' (none of the hammers compared reaches the target — consider relief drilling or a thicker pile tip)'}.` : gov.name.startsWith('Compressive') ? `Increase the wall to about ${fmt(Math.ceil(inp.t * gov.uc))} mm at the driving head, reduce the hammer energy setting, or add a cushion.` : gov.name.startsWith('Pile run') ? 'Hold the pile on the crane or gripper through the weak layer and lower the hammer energy before break-through.' : 'Reduce the stress range with a lower energy setting or improve the weld class (ground flush girth welds).',
            'API RP 2A-WSD (2014) §6.10, §12.5; DNV-ST-N001 (2021) Sec. 18'),
        adv('economic', 'info', `Net driving time is about ${fmt(hours)} h per pile with the ${h0.ratedEnergy / 1e3} kJ hammer; the smallest hammer compared that reaches the target is ${fmt(Ereq)} kJ.`,
          Ereq < h0.ratedEnergy / 1e3 ? 'A smaller hammer is sufficient: check its spread cost and handling weight against the longer driving time.' : 'Hammer size is set by the upper-bound resistance: a site-specific SRD study (CPT-based) often narrows the bounds enough to avoid mobilising the next hammer size or a drilling spread as contingency.',
          'Installation spread day rate × (handling + driving time) per pile; contingency for refusal priced as relief drilling (AACE RP 18R-97 estimate classes)',
          `Pile steel ≈ ${steelCo2Text(steel_t)} Driving time on the installation vessel: ${vesselCo2Text(hours / 24 + 0.3, 40)}`),
      ],
      publish: { 'install.pileBlows': blowsEnd, 'install.pileDrivingHours': hours, 'install.pileSRD': srdEnd },
      quantities: { steel_t: round(steel_t), vessel_days: { piling: round(hours / 24 + 0.3) } },
    };
  },
  verify: [
    { name: 'SRD in uniform clay, coring upper bound: Fp·∫α·su·dz·π·(Do + Di) + 9·su·A_annulus', source: 'Hand calculation with the API RP 2GEO alpha method (α = 0.5·ψ^-0.5 for ψ ≤ 1, 0.5·ψ^-0.25 for ψ > 1, α ≤ 1) and Fp = 0.5 (OCR 1, Stevens et al. 1982): su = 50 kPa, γ′ = 8 kN/m³, D = 2.134 m, t = 60 mm, 40 m → ∫f dz = 1603.3 kN/m, shaft = 0.5·1603.3·π·4.148 = 10 446 kN, toe = 9·50·0.3909 = 175.9 kN → 10.62 MN',
      inputs: { layers: [{ type: 'clay', thickness: 80, gammaSub: 8, su: 50, suGrad: 0, phi: 0, Dr: 0, OCR: 1 }], pen: 40, Lpile: 50, nStation: 3, nSeg: 12 },
      expect: [{ kpi: 'srdEnd', value: (0.5 * Math.PI * (2 * 2.134 - 0.12) * (0.5 * 50e3 ** 0.75 * 8000 ** 0.25 * 6.25 ** 1.25 / 1.25 + 0.5 * Math.sqrt(50e3 * 8000) * (2 / 3) * (25 ** 1.5 - 6.25 ** 1.5) + 50e3 * 15) + 9 * 50e3 * Math.PI / 4 * (2.134 ** 2 - 2.014 ** 2)) / 1e6, tol: 0.005 }] },
    { name: 'Ram impact velocity and rod impact stress: v0 = √(2·η·E/m), σ = ρ·c·v0', source: 'Energy balance and one-dimensional wave mechanics (Timoshenko & Goodier §168): 1200 kJ, η = 0.9, ram 60 t → v0 = 6.0 m/s; c = √(210e9/7850) = 5172.2 m/s → σ = 7850·5172.2·6.0 = 243.6 MPa',
      inputs: { nStation: 3, nSeg: 12 }, expect: [{ kpi: 'v0', value: 6.0, tol: 1e-9 }, { kpi: 'sigImp', value: 7850 * Math.sqrt(210e9 / 7850) * 6 / 1e6, tol: 1e-9 }] },
    { name: 'Energy transfer to a long pile by a rigid ram equals the ram kinetic energy', source: 'Wave mechanics: the ram decelerates on the pile impedance within M/Z = 3.8 ms, long before the toe reflection returns (2L/c = 29 ms), so ENTHRU → η·E = 1080 kJ less the lumped-head loss, which vanishes with refinement (within 10 % at 48 segments)',
      inputs: { nStation: 3, nSeg: 48 }, expect: [{ kpi: 'enthru', value: 1080, tol: 0.1 }] },
  ],
};
