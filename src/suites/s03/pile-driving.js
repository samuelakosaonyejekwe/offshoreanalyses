// Pile installation by impact driving: soil resistance to driving, one-dimensional wave-equation analysis of each
// blow, blow counts and driving stresses against penetration, refusal, self-weight penetration and pile-run
// screening, driving fatigue, and installation time per foundation.
import { GENERIC_HAMMERS, impactVelocity, srdProfile, drivabilityAnalysis, selfWeightPenetration, pileSubmergedWeight } from '../../engine/geotech/driving.js';
import { snCurve, snConstants } from '../../engine/fatigue/sn.js';
import { minerSum } from '../../engine/fatigue/damage.js';
import { G, need, run, fmt, round, adv, check, worst, sev, vesselCo2Text, tube, supportInputs, steelInputs, steelOf, soilInputs, soilProfile, soilLabel, LIMIT_HIFI, CONFIRM } from './_common.js';

const geo = supportInputs().filter((i) => ['depth', 'D', 't', 'Lp'].includes(i.key));
const HAMMERS = Object.entries(GENERIC_HAMMERS).map(([value, h]) => ({ value, label: h.label }));

export default {
  id: 's03-pile-driving', suite: 3, title: 'Pile driving and penetration',
  summary: 'Drivability of the monopile: soil resistance to driving, blow count, driving stresses and transferred energy against penetration, refusal and pile-run risk, fatigue consumed by driving and the installation time.',
  covers: [14], fidelity: 'numerical',
  standards: ['DNV-ST-0126 (2021) §7.6.6 and §8 (installation)', 'API RP 2A-WSD (21st ed.) §6.10 and §12.5 (drivability, driving stresses)', 'DNV-RP-C203 (2021) §2 (fatigue from pile driving)', 'DNV-ST-N001 (2021) marine operations', 'Smith (1960) ASCE SM4; Stevens, Wiltsie & Turton (1982) OTC 4205'],
  models: ['One-dimensional stress-wave propagation in the pile', 'Smith wave-equation pile-driving model (ram, cushion, helmet, lumped pile, elastic–plastic soil springs with damping)', 'Soil resistance to driving after Stevens et al. (coring and plugged, lower and upper bound)', 'Hammer–pile–soil interaction and transferred energy (ENTHRU)', 'Self-weight penetration and pile-run screening', 'Palmgren–Miner damage from driving stress cycles'],
  limits: `Continuous driving without set-up after delays, no residual-stress analysis and a rigid-ram hammer model with generic hammer data (not manufacturer data). The Alm & Hamre friction-fatigue SRD, vibratory driving and underwater noise are not included. Refusal and stress limits are contractual: confirm them for the project. ${LIMIT_HIFI}`,
  publishes: ['install.vesselDays'],
  searchTerms: 'monopile drivability wave equation SRD blow count driving fatigue installation',
  inputs: [
    ...geo,
    { key: 'Ltot', label: 'Total pile length', unit: 'm', type: 'number', default: 68, min: 15, max: 140, group: 'Support structure', role: 'geometry', help: 'Toe to pile top at the end of driving: embedded length + water depth + stick-up.' },
    ...steelInputs(),
    ...soilInputs(),
    { key: 'hammer', label: 'Hammer', type: 'select', default: 'generic-3000', options: HAMMERS, group: 'Hammer', role: 'bc', help: 'Generic hydraulic hammers defined by rated energy and ram mass (pile-head boundary condition). Replace with the contractor\'s hammer data.' },
    { key: 'eff', label: 'Hammer efficiency', unit: '-', type: 'number', default: 0.9, min: 0.3, max: 1, group: 'Hammer', role: 'bc', uncertain: { dist: 'normal', cov: 0.08 }, help: 'Kinetic energy at impact over rated energy; 0.85–0.95 for hydraulic hammers, lower under water or at reduced stroke.' },
    { key: 'energyFrac', label: 'Energy setting', unit: '%', type: 'number', default: 100, min: 10, max: 100, group: 'Hammer', role: 'operation', design: true, help: 'Reduced energy lowers driving stresses and noise in the easy part of the drive.' },
    { key: 'helmetT', label: 'Anvil and sleeve mass on the pile', unit: 't', type: 'number', default: 80, min: 0, max: 600, group: 'Hammer', role: 'bc' },
    { key: 'cushion', label: 'Cushion stiffness (0 = steel on steel)', unit: 'MN/m', type: 'number', default: 0, min: 0, max: 1e6, group: 'Hammer', role: 'bc' },
    { key: 'hammerT', label: 'Hammer and sleeve weight on the pile', unit: 't', type: 'number', default: 500, min: 0, max: 2500, group: 'Hammer', role: 'bc', help: 'Rests on the pile during self-weight penetration.' },
    { key: 'srdCase', label: 'Resistance case for the drive', type: 'select', default: 'coringUpper', group: 'Soil resistance to driving', role: 'soil',
      options: [{ value: 'coringLower', label: 'Coring, lower bound' }, { value: 'coringUpper', label: 'Coring, upper bound' }, { value: 'pluggedLower', label: 'Plugged, lower bound' }, { value: 'pluggedUpper', label: 'Plugged, upper bound' }],
      help: 'Large open piles core through the soil; the plugged cases bound the resistance of small piles in dense sand.' },
    { key: 'quake', label: 'Soil quake', unit: 'mm', type: 'number', default: 2.54, min: 0.5, max: 10, group: 'Soil resistance to driving', role: 'soil', help: 'Elastic displacement at which the static resistance is mobilised (Smith: 2.54 mm).' },
    { key: 'Js', label: 'Shaft damping', unit: 's/m', type: 'number', default: 0.16, min: 0, max: 1, group: 'Soil resistance to driving', role: 'soil', help: 'Smith: 0.16 s/m in sand, about 0.65 s/m in clay.' },
    { key: 'Jt', label: 'Toe damping', unit: 's/m', type: 'number', default: 0.49, min: 0, max: 1.5, group: 'Soil resistance to driving', role: 'soil' },
    { key: 'refusal', label: 'Refusal criterion', unit: 'blows / 0.25 m', type: 'number', default: 250, min: 20, max: 2000, group: 'Acceptance criteria', role: 'criteria', help: `Contractual limit; 250 blows per 0.25 m (≈ 300 per foot) is a common value. ${CONFIRM}` },
    { key: 'stressLim', label: 'Allowable driving stress / yield', unit: '-', type: 'number', default: 0.9, min: 0.5, max: 1, group: 'Acceptance criteria', role: 'criteria', help: 'API RP 2A limits dynamic driving stress to 80–90 % of yield.' },
    { key: 'scf', label: 'Stress concentration factor at girth welds', unit: '-', type: 'number', default: 1.1, min: 1, max: 3, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'dAllow', label: 'Fatigue damage budget for driving', unit: '-', type: 'number', default: 0.1, min: 0.01, max: 1, group: 'Acceptance criteria', role: 'criteria', help: 'Share of the Miner sum reserved for installation, before the design fatigue factor.' },
    { key: 'bpm', label: 'Blow rate', unit: 'blows/min', type: 'number', default: 35, min: 10, max: 60, group: 'Installation spread', role: 'operation' },
    { key: 'tFixed', label: 'Other vessel time per foundation', unit: 'h', type: 'number', default: 20, min: 0, max: 200, group: 'Installation spread', role: 'operation', help: 'Positioning, upending, gripper and noise-mitigation deployment, hammer handling, transition piece or flange work, transit.' },
    { key: 'wow', label: 'Weather downtime factor', unit: '-', type: 'number', default: 1.3, min: 1, max: 3, group: 'Installation spread', role: 'operation', help: 'Total time over working time for the season and limiting sea state of the operation.' },
    { key: 'nFnd', label: 'Number of foundations', type: 'int', default: 1, min: 1, max: 300, group: 'Installation spread', role: 'operation' },
    { key: 'nSeg', label: 'Pile segments in the wave equation', type: 'int', default: 50, min: 20, max: 400, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nSeg', label: 'Pile segments', levels: [25, 50, 100, 200], qoi: 'totalBlows' },
  compute(inp, ctx) {
    const steel = steelOf(inp, ctx);
    need(inp.Ltot >= inp.Lp, 'The total pile length cannot be shorter than the target penetration.');
    const sec = tube(inp.D, inp.t / 1e3);
    const profile = soilProfile(inp);
    need(profile.depth >= inp.Lp, 'The soil profile is shallower than the target penetration: extend the last layer.');
    const h0 = GENERIC_HAMMERS[inp.hammer];
    const hammer = { ramMass: h0.ramMass, ratedEnergy: h0.ratedEnergy * inp.energyFrac / 100, efficiency: inp.eff };
    const v0 = impactVelocity(hammer);
    const pile = { length: inp.Ltot, D: inp.D, t: inp.t / 1e3, E: steel.E, rho: steel.rho, nSegments: inp.nSeg };
    const srd = run(() => srdProfile({ profile, D: inp.D, t: inp.t / 1e3, Lmax: inp.Lp, dz: 0.5 }), 'Soil resistance to driving');
    const nD = Math.min(12, Math.max(4, Math.round(inp.Lp / 3)));
    const depths = Array.from({ length: nD }, (_, i) => round((inp.Lp * (i + 1)) / nD, 6));
    const dr = run(() => drivabilityAnalysis({ profile, hammer, pile, srdCase: inp.srdCase, depths, helmetMass: inp.helmetT * 1e3, cushionStiffness: inp.cushion > 0 ? inp.cushion * 1e6 : undefined,
      quake: { shaft: inp.quake / 1e3, toe: inp.quake / 1e3 }, damping: { shaftClay: inp.Js, shaftSand: inp.Js, toeClay: inp.Jt, toeSand: inp.Jt }, refusalBlowsPer025m: inp.refusal }), 'Drivability analysis');
    const last = depths.length - 1;
    const fin = (a) => Array.from(a, (v) => (Number.isFinite(v) ? v : inp.refusal * 4));
    const b025 = fin(dr.blowsPer025m);
    const sMax = Math.max(...dr.maxCompressiveStress), sTen = Math.max(...dr.maxTensileStress);
    const refused = dr.refusalDepth != null;
    // Self-weight penetration with the lower-bound resistance (deepest plausible), pile then pile + hammer.
    const Wp = pileSubmergedWeight({ D: inp.D, t: inp.t / 1e3, length: inp.Ltot, rho: steel.rho, rhoWater: 1025 });
    const sw = selfWeightPenetration({ z: Array.from(srd.z), srd: Array.from(srd.cases.coringLower.total), pileWeight: Wp, hammerWeight: inp.hammerT * 1e3 * G });
    // Driving fatigue at the girth welds (curve D in air, as-welded, thickness-corrected).
    const sn = snConstants(snCurve('D', 'air'), { thickness: sec.t, scf: inp.scf });
    const blocks = dr.fatigueBlocks || [];
    const dmg = blocks.length ? minerSum({ ranges: blocks.map((b) => b.stressRange), counts: blocks.map((b) => (Number.isFinite(b.blows) ? b.blows : 0)), sn }).damage : 0;
    const totalBlows = Number.isFinite(dr.totalBlows) ? dr.totalBlows : blocks.reduce((a, b) => a + (Number.isFinite(b.blows) ? b.blows : 0), 0);
    const tDrive = totalBlows / inp.bpm / 60; // h
    const perFnd = (tDrive + inp.tFixed) * inp.wow / 24, vdays = perFnd * inp.nFnd;
    const checks = [
      check('Blow count at final penetration against refusal', b025[last], inp.refusal, 'blows / 0.25 m', 'API RP 2A-WSD §12.5.6 (refusal definition); project installation specification', 'Guards against the pile stopping short of target penetration. At refusal the hammer cannot advance the pile and relief drilling or a larger hammer is needed.'),
      check('Maximum compressive driving stress', sMax / 1e6, inp.stressLim * steel.fy / 1e6, 'MPa', 'API RP 2A-WSD §6.10.5 (dynamic stresses ≤ 80–90 % of yield)', 'Guards the pile head and toe against yielding and local buckling during driving. Overstress causes mushrooming of the head or toe damage that can stop the drive.'),
      check('Fatigue damage consumed by driving', dmg, inp.dAllow, '-', 'DNV-RP-C203 (2021) §2.2; DNV-ST-0126 (2021) §4.11.2 (installation fatigue)', 'Every blow is one large stress cycle at the girth welds. Damage used during driving is no longer available for the service life.'),
    ];
    const gov = worst(checks), ok = gov.uc <= 1 && !refused;
    const cs = srd.cases;
    const zz = Array.from(srd.z);
    // Energy setting that brings the peak stress to the limit (stress ∝ impact velocity ∝ √energy).
    const eOpt = Math.min(100, inp.energyFrac * (inp.stressLim * steel.fy / sMax) ** 2);
    return {
      kpis: [
        { key: 'srdFinal', label: 'Soil resistance to driving at final penetration', unit: 'MN', value: dr.srd[last] / 1e6, explain: 'Static resistance the hammer must overcome at target depth for the chosen resistance case; compare with the lower and upper bounds in the chart.' },
        { key: 'blowsFinal', label: 'Blow count at final penetration', unit: 'blows / 0.25 m', value: b025[last], limit: inp.refusal, sense: 'max', explain: 'Hammer blows needed for the last quarter metre. Easy driving is below 50; above about 125 driving is hard and refusal is near.' },
        { key: 'totalBlows', label: 'Total number of blows', unit: '-', value: totalBlows, explain: 'Blows from the start of driving to target penetration (or to refusal). It sets driving time, fatigue damage and the cumulative noise dose.' },
        { key: 'stressMax', label: 'Maximum compressive driving stress', unit: 'MPa', value: sMax / 1e6, limit: inp.stressLim * steel.fy / 1e6, sense: 'max', explain: 'Largest compressive stress anywhere in the pile over the drive. A rigid ram without cushion overstates the first peak by up to about 20 % in this lumped model.' },
        { key: 'stressTen', label: 'Maximum tensile driving stress', unit: 'MPa', value: sTen / 1e6, explain: 'Largest tensile stress from the wave reflected at the toe in easy driving. It matters for bolted or grouted connections and for welds with defects.' },
        { key: 'stressFirst', label: 'Compressive stress in the first analysed stage', unit: 'MPa', value: dr.maxCompressiveStress[0] / 1e6, explain: 'Peak stress at shallow penetration where soil resistance is small: close to the theoretical impact stress ρ·c·v0 of a rigid ram on an elastic rod.' },
        { key: 'enthru', label: 'Energy transferred to the pile (final stage)', unit: 'kJ', value: dr.enthru[last] / 1e3, explain: 'Energy that actually passes the pile head per blow (ENTHRU). The ratio to the rated energy is the overall transfer efficiency to compare with driving records.' },
        { key: 'v0', label: 'Ram impact velocity', unit: 'm/s', value: v0, explain: 'Velocity of the ram at impact from the delivered kinetic energy, √(2ηE/m). It sets the peak driving stress.' },
        { key: 'swp', label: 'Self-weight penetration (pile + hammer)', unit: 'm', value: sw.swpPileHammer ?? inp.Lp, explain: 'Depth the pile sinks under its own weight and the hammer before any blow, with the lower-bound resistance. It must give enough lateral stability before the gripper releases.' },
        { key: 'refusalDepth', label: 'Refusal depth (0 = none)', unit: 'm', value: dr.refusalDepth ?? 0, explain: 'Penetration at which the blow count first exceeds the refusal criterion; zero means the pile reaches target depth.' },
        { key: 'fatDamage', label: 'Fatigue damage from driving', unit: '-', value: dmg, limit: inp.dAllow, sense: 'max', explain: 'Miner sum at the girth welds from all blows (curve D in air with thickness correction and the stated SCF), using the full compressive-plus-tensile stress range of each stage — an upper estimate.' },
        { key: 'driveHours', label: 'Net driving time', unit: 'h', value: tDrive, explain: 'Hammering time at the stated blow rate, without stops for noise mitigation checks or hammer cooling.' },
        { key: 'vesselDays', label: 'Installation vessel time', unit: 'days', value: vdays, explain: `Driving plus other time per foundation, times the weather factor and the number of foundations (${fmt(perFnd)} days each).` },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Blow count against penetration', x: { label: 'Blows per 0.25 m', unit: '-' }, y: { label: 'Penetration', unit: 'm' },
          series: [{ name: 'Predicted', x: b025, y: depths.map((v) => -v) }, { name: 'Refusal', x: [inp.refusal, inp.refusal], y: [0, -inp.Lp], dash: true }],
          explain: 'Driving record prediction. The curve must stay left of the refusal line to target depth; a sharp rise marks a hard layer where energy should be at maximum.' },
        { type: 'line', title: 'Soil resistance to driving', x: { label: 'Resistance', unit: 'MN' }, y: { label: 'Penetration', unit: 'm' },
          series: [['coringLower', 'Coring, lower'], ['coringUpper', 'Coring, upper'], ['pluggedLower', 'Plugged, lower'], ['pluggedUpper', 'Plugged, upper']].map(([k, n]) => ({ name: n, x: Array.from(cs[k].total, (v) => v / 1e6), y: zz.map((v) => -v), dash: k !== inp.srdCase })).concat([{ name: 'Pile + hammer weight', x: [(Wp + inp.hammerT * 1e3 * G) / 1e6, (Wp + inp.hammerT * 1e3 * G) / 1e6], y: [0, -inp.Lp], dash: true }]),
          explain: 'Static resistance during driving for the four bounding cases, with the weight of pile and hammer. Where the weight line exceeds the lower-bound curve the pile penetrates without blows; a fall of resistance below the weight at depth signals pile-run risk.' },
        { type: 'line', title: 'Driving stresses against penetration', x: { label: 'Stress', unit: 'MPa' }, y: { label: 'Penetration', unit: 'm' },
          series: [{ name: 'Maximum compression', x: Array.from(dr.maxCompressiveStress, (v) => v / 1e6), y: depths.map((v) => -v) }, { name: 'Maximum tension', x: Array.from(dr.maxTensileStress, (v) => v / 1e6), y: depths.map((v) => -v) }, { name: 'Allowable', x: [inp.stressLim * steel.fy / 1e6, inp.stressLim * steel.fy / 1e6], y: [0, -inp.Lp], dash: true }],
          explain: 'Peak stresses in the pile per blow. Compression grows with toe resistance; tension is largest in easy driving when the wave reflects from a free toe.' },
      ],
      tables: [
        { title: 'Drivability by penetration', columns: [{ key: 'z', label: 'Penetration', unit: 'm' }, { key: 'srd', label: 'SRD', unit: 'MN' }, { key: 'b', label: 'Blows / 0.25 m' }, { key: 'bm', label: 'Blows / m' }, { key: 'sc', label: 'Compression', unit: 'MPa' }, { key: 'st', label: 'Tension', unit: 'MPa' }, { key: 'e', label: 'ENTHRU', unit: 'kJ' }, { key: 'r', label: 'Refusal' }],
          rows: depths.map((z, i) => ({ z, srd: round(dr.srd[i] / 1e6), b: round(b025[i]), bm: round(4 * b025[i]), sc: round(dr.maxCompressiveStress[i] / 1e6), st: round(dr.maxTensileStress[i] / 1e6), e: round(dr.enthru[i] / 1e3), r: dr.refusal[i] ? 'yes' : 'no' })),
          explain: 'One wave-equation analysis per penetration stage: resistance, blow count, stresses and transferred energy, to compare with the pile-driving monitoring records.' },
      ],
      notes: [`Hammer: ${h0.label} at ${inp.energyFrac} % energy, efficiency ${inp.eff}; soil: ${soilLabel(inp)}; resistance case: ${inp.srdCase}.`,
        refused ? `Refusal is predicted at ${fmt(dr.refusalDepth)} m, short of the ${inp.Lp} m target.` : `The pile reaches ${inp.Lp} m with ${fmt(b025[last])} blows per 0.25 m at the end.`,
        sw.runRisk ? `Pile-run risk: resistance falls below the weight of pile and hammer between ${sw.runIntervals.map((r) => `${fmt(r.from)}–${fmt(r.to)} m`).join(', ')}.` : 'No pile-run interval is found with the lower-bound resistance.',
        `${gov.name} governs (utilisation ${fmt(gov.uc)}).`, 'Generic hammer data: repeat the analysis with the contractor\'s hammer model and site-calibrated resistance before the installation procedure is fixed.'],
      advice: [
        ok ? adv('engineering', sev(gov.uc), `The pile is drivable to ${inp.Lp} m: ${fmt(b025[last])} blows per 0.25 m at the end, ${fmt(totalBlows)} blows in total, peak stress ${fmt(sMax / 1e6)} MPa (${fmt(sMax / (inp.stressLim * steel.fy) * 100)} % of allowable), driving damage ${fmt(dmg)}.`,
          `Start at reduced energy and step up with resistance; ${sMax < inp.stressLim * steel.fy ? `the stress limit allows up to ${fmt(eOpt)} % energy` : 'keep energy as analysed'}. Carry the driving damage ${fmt(dmg)} into the fatigue budget and monitor blows with a pile-driving analyser on the first positions.`,
          'API RP 2A-WSD §6.10 and §12.5; DNV-ST-0126 (2021) §8; DNV-RP-C203 (2021) §2',
          `${fmt(totalBlows)} blows at ${fmt(hammer.ratedEnergy / 1e3)} kJ are the underwater-noise dose of this foundation: fewer, softer blows reduce both the noise-mitigation burden and the fatigue used. ${vesselCo2Text(perFnd)}`)
          : adv('engineering', 'fail', refused ? `Refusal at ${fmt(dr.refusalDepth)} m before the ${inp.Lp} m target with the ${h0.label}.` : `${gov.name}: ${fmt(gov.demand)} against ${fmt(gov.capacity)} ${gov.unit}.`,
            refused ? 'Select a larger hammer, allow a higher refusal criterion for the last metres, or plan relief drilling; check the upper-bound resistance against the hammer capacity.' : gov.name.includes('stress') ? `Reduce the energy setting to about ${fmt(eOpt)} % or add a cushion; thicken the top and toe cans if full energy is needed.` : `Reduce the number of high-stress blows (lower energy in easy driving), improve the weld class by grinding, or raise the driving share of the fatigue budget by thickening the can (now ${inp.t} mm).`,
            'API RP 2A-WSD §6.10.5 and §12.5.6; DNV-RP-C203 (2021) §2; DNV-ST-0126 (2021) §8',
            `A stalled pile costs vessel time: each extra day is ${vesselCo2Text(1)}`),
        adv('economic', 'info', `Installation takes about ${fmt(perFnd)} vessel days per foundation (${fmt(tDrive)} h of driving, ${inp.tFixed} h other work, weather factor ${inp.wow}), ${fmt(vdays)} days for ${inp.nFnd} foundation(s).`,
          'Driving is a small part of the cycle: the schedule is governed by weather windows and handling. Match the hammer to the upper-bound resistance so that no position needs a second mobilisation, and price refusal risk (drilling spread on standby) against one hammer size up.',
          'Marine-operation planning per DNV-ST-N001 (weather-restricted operations, α-factors); schedule-risk allowances per AACE RP 57R-09'),
      ],
      publish: { 'install.vesselDays': vdays },
      quantities: { steel_t: round(sec.A * steel.rho * inp.Ltot / 1e3 * inp.nFnd), vessel_days: { foundation_installation: round(vdays) } },
    };
  },
  verify: [
    { name: 'Ram impact velocity from the delivered energy: v0 = √(2ηE/m)', source: 'Energy balance (hand calculation): rated 3000 kJ, efficiency 0.9, ram 150 t → v0 = √(2 × 0.9 × 3.0e6 / 1.5e5) = 6.00 m/s', inputs: { hammer: 'generic-3000', eff: 0.9, energyFrac: 100 }, expect: [{ kpi: 'v0', value: 6.0, tol: 1e-6 }] },
    { name: 'Peak impact stress of a rigid ram on an elastic rod: σ = ρ·c·v0 = √(Eρ)·v0', source: 'One-dimensional wave theory (Timoshenko & Goodier §168): √(210e9 × 7850) × 6.0 m/s = 243.6 MPa. The lumped wave-equation model overshoots a sharp front by up to about 20 % without a cushion, so a tolerance of 20 % is applied',
      inputs: { hammer: 'generic-3000', eff: 0.9, energyFrac: 100, helmetT: 0, cushion: 0, nSeg: 100 }, expect: [{ kpi: 'stressFirst', value: 243.6, tol: 0.2 }] },
  ],
};
