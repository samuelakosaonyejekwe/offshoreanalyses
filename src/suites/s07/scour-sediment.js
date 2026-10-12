// Seabed mobility, sediment transport rates and local scour at bridge piers, towers and structure toes, with the
// scour-protection rock size and quantity.
import { criticalShearStress, settlingVelocitySoulsby, grainRoughness, bedShearCurrent, waveOrbital, bedShearWave, combinedShearSoulsby, mobilityCheck,
  meyerPeterMuller, engelundHansen, vanRijnBedLoad, vanRijnSuspendedLoad, vanRijnSimplified, soulsbyVanRijn, rouseProfile, isbashStoneSize } from '../../engine/coastal/sediment.js';
import { scourHec18, scourSheppardMelville, scourBreusers, pileScour, scourTimeScale, scourDevelopment, largeStructureScour, armourTakeoff } from '../../engine/geotech/scour.js';
import { G, DEG, CO2_ROCK, CO2_MGO, MGO_T_PER_DAY, need, run, finite, r, f, linspace, check, eng, eco, slrInputs, slrRows, SLR_BASIS } from './_common.js';

const safe = (fn, fallback = 0) => { try { const v = fn(); return Number.isFinite(v) ? v : fallback; } catch { return fallback; } };

export default {
  id: 's07-scour-sediment', suite: 7, title: 'Seabed scour and sediment transport',
  summary: 'Bed mobility, transport rates and the local scour depth at a pier, tower or wall toe, checked against foundation embedment, with the rock protection that stops it.',
  covers: [6], fidelity: 'analytical',
  standards: ['FHWA HEC-18 Evaluating Scour at Bridges (5th ed., 2012) ch. 7', 'FHWA HEC-23 (3rd ed., 2009) Design Guideline 11', 'Soulsby, Dynamics of Marine Sands (1997)', 'Sumer & Fredsøe, The Mechanics of Scour in the Marine Environment (2002)', 'DNV-ST-0126 (2021) §7.6', 'USACE CEM Part VI-5-6'],
  models: ['Shields sediment initiation', 'Soulsby wave–current bed shear stress', 'Meyer-Peter–Müller bed load', 'Engelund–Hansen total load', 'Van Rijn bed and suspended load', 'Soulsby–Van Rijn wave–current transport', 'Rouse suspension profile', 'HEC-18 (CSU) pier scour', 'Sheppard–Melville pier scour', 'Sumer–Fredsøe wave–current pile scour', 'Isbash rock sizing'],
  limits: 'Equilibrium scour in uniform non-cohesive sediment from empirical formulae for isolated piers; cohesive beds, layered or armoured beds, pile groups with complex caps, contraction scour and long-term bed-level change need site-specific assessment or physical modelling. Toe scour at walls is a rule-of-thumb coefficient, not a prediction.',
  publishes: ['marine.scourDepth', 'marine.scourRockD50'],
  searchTerms: 'bridge pier scour HEC-18 Sheppard Melville sediment transport Shields scour protection',
  inputs: [
    { key: 'a', label: 'Pier / tower width normal to the flow', unit: 'm', type: 'number', default: 3, min: 0.2, max: 40, group: 'Structure', role: 'geometry', design: true },
    { key: 'Lp', label: 'Pier length along the flow', unit: 'm', type: 'number', default: 9, min: 0.2, max: 120, group: 'Structure', role: 'geometry' },
    { key: 'nose', label: 'Nose shape', type: 'select', default: 'round', options: [{ value: 'round', label: 'Round nose' }, { value: 'circular', label: 'Circular cylinder' }, { value: 'square', label: 'Square nose' }, { value: 'sharp', label: 'Sharp nose' }, { value: 'group', label: 'Group of cylinders' }], group: 'Structure', role: 'geometry' },
    { key: 'skewDeg', label: 'Angle of attack of the flow', unit: 'deg', type: 'number', default: 5, min: 0, max: 45, group: 'Structure', role: 'bc', help: 'Skew multiplies scour of long piers quickly: 15° on a 3:1 pier adds about 50 %.' },
    { key: 'embed', label: 'Foundation embedment below the bed', unit: 'm', type: 'number', default: 20, min: 0.5, max: 120, group: 'Structure', role: 'ic', design: true, help: 'Existing or proposed pile / caisson embedment at the initial bed level.' },
    { key: 'allow', label: 'Allowable scour as a share of embedment', unit: '-', type: 'number', default: 0.3, min: 0.05, max: 0.8, group: 'Criteria', role: 'criteria', help: 'Project criterion from the foundation design: the embedment the capacity check can afford to lose.' },
    { key: 'h', label: 'Water depth', unit: 'm', type: 'number', default: 8, min: 0.5, max: 80, group: 'Flow and waves', role: 'bc', bind: 'site.depth' },
    { key: 'U', label: 'Depth-averaged current', unit: 'm/s', type: 'number', default: 1.4, min: 0.05, max: 6, group: 'Flow and waves', role: 'bc', bind: 'site.Uc100', uncertain: { dist: 'lognormal', cov: 0.15 } },
    { key: 'Hs', label: 'Significant wave height', unit: 'm', type: 'number', default: 2.5, min: 0, max: 15, group: 'Flow and waves', role: 'bc', bind: 'marine.HsToe', uncertain: { dist: 'lognormal', cov: 0.1 } },
    { key: 'Tp', label: 'Peak wave period', unit: 's', type: 'number', default: 9, min: 2, max: 22, group: 'Flow and waves', role: 'bc', bind: 'marine.TpToe' },
    { key: 'phiDeg', label: 'Angle between waves and current', unit: 'deg', type: 'number', default: 30, min: 0, max: 90, group: 'Flow and waves', role: 'bc' },
    ...slrInputs({ group: 'Flow and waves' }),
    { key: 'd50', label: 'Median grain size d50', unit: 'mm', type: 'number', default: 0.3, min: 0.063, max: 60, group: 'Seabed', role: 'soil', bind: 'site.d50', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'd90r', label: 'd90 / d50', unit: '-', type: 'number', default: 2, min: 1.1, max: 6, group: 'Seabed', role: 'soil' },
    { key: 'rhoS', label: 'Sediment grain density', unit: 'kg/m³', type: 'number', default: 2650, min: 1800, max: 3200, group: 'Seabed', role: 'soil' },
    { key: 'general', label: 'General / contraction scour and bed lowering', unit: 'm', type: 'number', default: 0.5, min: 0, max: 20, group: 'Seabed', role: 'ic', help: 'Long-term degradation plus contraction scour, from surveys or a river-engineering study.' },
    { key: 'K3', label: 'HEC-18 bed-condition factor K3', unit: '-', type: 'number', default: 1.1, min: 1.1, max: 1.3, group: 'Seabed', role: 'criteria' },
    { key: 'kToe', label: 'Wall-toe scour depth / wave height', unit: '-', type: 'number', default: 1.0, min: 0, max: 2, group: 'Criteria', role: 'criteria', help: 'Rule of thumb for reflective walls (USACE CEM VI-5-6): maximum toe scour about equal to the unbroken wave height at the wall. Calibrate against surveys.' },
    { key: 'rock', label: 'Scour-protection rock', type: 'select', source: 'materials:rock', default: 'rock-granite', group: 'Protection', role: 'material' },
    { key: 'Kv', label: 'Velocity amplification at the pier', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Protection', role: 'criteria', help: 'HEC-23: 1.5 round-nose, 1.7 square-nose.' },
    { key: 'extent', label: 'Protection radius / pier width', unit: '-', type: 'number', default: 2.5, min: 1, max: 6, group: 'Protection', role: 'geometry', design: true },
    { key: 'stormH', label: 'Duration of the design flow', unit: 'h', type: 'number', default: 12, min: 0.5, max: 720, group: 'Flow and waves', role: 'bc' },
  ],
  compute(inp, ctx) {
    const d = inp.d50 / 1000, rhoS = inp.rhoS, h = inp.h + inp.slr;
    const rk = ctx.material(inp.rock); need(rk && rk.rho > 1500, 'The selected rock has no density in the materials library.');
    const at = (hh) => {
      const crit = criticalShearStress({ d, rhoS });
      const z0 = grainRoughness(d).z0;
      const tc = bedShearCurrent({ U: inp.U, h: hh, z0 });
      const orb = inp.Hs > 0 ? waveOrbital({ H: inp.Hs / Math.SQRT2, T: inp.Tp, h: hh }) : { Uw: 0, A: 0 };
      const tw = orb.Uw > 1e-4 ? bedShearWave({ Uw: orb.Uw, T: inp.Tp, z0 }).tau : 0;
      const comb = tw > 0 ? combinedShearSoulsby({ tauC: tc.tau, tauW: tw, phiDeg: inp.phiDeg }) : { tauMean: tc.tau, tauMax: tc.tau };
      const mob = mobilityCheck({ tau: comb.tauMax, d, rhoS });
      // threshold velocity: Van Rijn's fit inside its 0.1–2 mm range, otherwise Shields with the log-law drag coefficient
      let sim; try { sim = vanRijnSimplified({ U: inp.U, h: hh, d50: d, d90: d * inp.d90r, rhoS }); } catch { sim = { Ucr: crit.uStarCr / Math.sqrt(tc.Cd) }; }
      const hec = scourHec18({ a: inp.a, y1: hh, V1: inp.U, nose: inp.nose, skew: inp.skewDeg * DEG, L: Math.min(inp.Lp, 12 * inp.a), K3: inp.K3 }).ys;
      // effective width for the other methods: projected width under skew
      const ae = inp.a * Math.cos(inp.skewDeg * DEG) + inp.Lp * Math.sin(inp.skewDeg * DEG);
      const sm = safe(() => scourSheppardMelville({ a: ae, y1: hh, V1: inp.U, Vc: sim.Ucr, D50: d }).ys, 0);
      const br = safe(() => scourBreusers({ D: ae, h: hh }), 0);
      const KC = orb.Uw * inp.Tp / ae, Ucw = inp.U / (inp.U + orb.Uw);
      const sf = pileScour({ KC: Math.max(KC, 1e-6), Ucw }).SoverD * ae;
      return { crit, z0, tc, orb, tw, comb, mob, sim, hec, sm, br, sf, KC, Ucw, ae };
    };
    const s = run(() => at(h), 'Scour and transport');
    const local = Math.max(s.hec, s.sm), total = local + inp.general, toe = inp.Hs > 0 && inp.kToe > 0 ? largeStructureScour({ Lref: inp.Hs, kS: inp.kToe }).S : 0;
    finite(local, 'Local scour depth');
    const mobile = s.mob.mobile;
    // transport rates (solid volume per metre width)
    const mpm = safe(() => meyerPeterMuller({ theta: s.tc.tau / ((rhoS - 1025) * G * d), d, rhoS }).qb);
    const eh = safe(() => engelundHansen({ U: inp.U, tau: s.tc.tau, d50: d, rhoS }).qt);
    const vrb = safe(() => vanRijnBedLoad({ U: inp.U, h, d50: d, d90: d * inp.d90r, rhoS }).qb);
    let vrsO = null; try { vrsO = vanRijnSuspendedLoad({ U: inp.U, h, d50: d, d90: d * inp.d90r, rhoS }); } catch { vrsO = null; }
    const vrs = vrsO && Number.isFinite(vrsO.qs) ? vrsO.qs : 0;
    const svr = safe(() => soulsbyVanRijn({ U: inp.U, Urms: s.orb.Uw, h, d50: d, d90: d * inp.d90r, rhoS }).qt);
    const ws = settlingVelocitySoulsby({ d, rhoS });
    const day = 86400;
    // development in time (current-dominated time scale; live bed only)
    const theta = s.mob.theta;
    const T = mobile ? safe(() => scourTimeScale({ mode: 'current', D: s.ae, theta: Math.max(s.tc.tau / ((rhoS - 1025) * G * d), 1e-3), d, s: rhoS / 1025, delta: h }).T, 0) : 0;
    const tEnd = inp.stormH * 3600, frac = T > 0 ? scourDevelopment(1, tEnd, T) : (mobile ? 1 : 0);
    // protection: Isbash with the amplified velocity (HEC-23 DG11), layer 3 D50 thick
    const delta = rk.rho / 1025 - 1;
    const Uamp = inp.Kv * (inp.U + 0.5 * s.orb.Uw);
    const D50r = isbashStoneSize({ U: Uamp, delta, C: 1.2 });
    const tRock = Math.max(3 * D50r, 0.3);
    const take = armourTakeoff({ D: inp.a, extentFactor: Math.max(inp.extent, 0.6), thickness: tRock, porosity: 0.4, rhoS: rk.rho });
    const rockT = take.mass / 1000, days = Math.max(1, rockT / 1500);
    const cap = inp.allow * inp.embed, ok = total <= cap;
    const scen = slrRows(inp).map((sc) => { const q = at(inp.h + sc.rise); return { name: sc.name, rise: sc.rise, Uw: r(q.orb.Uw, 3), theta: r(q.mob.ratio, 3), S: r(Math.max(q.hec, q.sm), 3) }; });
    const us = linspace(0.2 * inp.U, 2 * inp.U, 30);
    const zs = vrsO ? linspace(vrsO.a, h, 25) : [];

    const advice = [];
    advice.push(ok
      ? eng('ok', `Design scour is ${f(total)} m (${f(local)} m local + ${f(inp.general)} m general) against ${f(cap)} m allowed — ${f(total / inp.embed * 100, 2)} % of the ${f(inp.embed)} m embedment.`,
        'The foundation tolerates the predicted scour unprotected. Survey the bed after major floods or storms and set a trigger depth at the allowable value.', 'FHWA HEC-18 (2012) ch. 7 and ch. 10 (monitoring); DNV-ST-0126 (2021) §7.6',
        `Leaving the bed unprotected avoids about ${f(rockT, 3)} t of rock (≈ ${f(rockT * CO2_ROCK, 2)} t CO₂e, screening factor) and the disturbance of ${f(take.area, 3)} m² of seabed habitat per pier.`)
      : eng('fail', `Design scour is ${f(total)} m (${f(local)} m local + ${f(inp.general)} m general), ${f(total / cap, 2)} × the ${f(cap)} m the foundation can lose (${f(total / inp.embed * 100, 2)} % of embedment).`,
        `Either deepen the foundation by about ${f((total - cap) / inp.allow, 2)} m, or place scour protection: D50 ≥ ${f(D50r, 2)} m rock, ${f(tRock, 2)} m thick, out to ${f(take.outerRadius)} m from the pier axis on a filter layer.`,
        'FHWA HEC-18 (2012) ch. 7; HEC-23 (2009) DG11; Sumer & Fredsøe (2002)',
        `Protection needs about ${f(rockT, 3)} t of rock per pier (≈ ${f(rockT * CO2_ROCK, 2)} t CO₂e, screening factor) over ${f(take.area, 3)} m² of seabed; a streamlined or aligned pier reduces scour at source with no added material.`));
    advice.push(eng(mobile ? 'info' : 'ok', `The bed is ${mobile ? 'mobile' : 'stable'} in the design flow: Shields ratio θmax/θcr = ${f(s.mob.ratio)} (${mobile ? 'live-bed' : 'clear-water'} scour; critical velocity ${f(s.sim.Ucr, 2)} m/s). ${T > 0 ? `About ${f(frac * 100, 2)} % of equilibrium scour develops in the ${f(inp.stormH)} h event (time scale ${f(T / 3600, 2)} h).` : ''}`,
      mobile ? 'Treat equilibrium scour as reached in fine sand; in clear-water conditions (ratio < 1) the largest scour occurs near the threshold velocity, so check the lower flows as well.' : 'Check a range of flows up to the threshold: clear-water scour peaks at U ≈ Ucr.', 'Soulsby (1997) ch. 6; HEC-18 (2012) §7.2'));
    advice.push(eco('info', `Scour protection would take ${f(rockT, 3)} t of ${f(D50r * 1000, 2)} mm rock per pier and about ${f(days, 2)} vessel day(s) to place; the wall-toe rule of thumb gives ${f(toe, 2)} m of toe scour.`,
      'Compare three options on whole-life cost: deeper foundation (steel/concrete once), rock protection (placing plus periodic top-up), or monitoring with a remedial contract held ready. Risk-rank piers so that inspection goes where scour consumes most embedment.',
      'Whole-life costing per ISO 15686-5; risk-based inspection per HEC-18 ch. 10 / API RP 580 principles',
      `Placing burns roughly ${f(days * MGO_T_PER_DAY, 2)} t of marine gas oil (≈ ${f(days * MGO_T_PER_DAY * CO2_MGO, 2)} t CO₂) per pier; monitoring-led management defers that until the bed actually moves.`));
    const worst = scen[scen.length - 1];
    advice.push(eng('info', `Sea-level rise deepens the water: at +${f(worst.rise)} m the orbital velocity is ${f(worst.Uw, 2)} m/s and the local scour ${f(worst.S)} m.`, 'Scour at piers is not sensitive to the rise itself; watch instead for changed tidal currents, storminess and bed-level trends, and resurvey on the adaptation review cycle.', SLR_BASIS));

    return {
      kpis: [
        { key: 'scour', label: 'Design scour depth (local + general)', value: r(total), unit: 'm', limit: r(cap), sense: 'max', explain: 'Bed lowering at the foundation in the design event. It removes lateral and axial soil support, so it is compared with the embedment the foundation can afford to lose.' },
        { key: 'scourHec', label: 'Local scour — HEC-18 (CSU)', value: r(s.hec), unit: 'm', explain: 'Pier scour from the CSU equation: depends on pier width, depth, Froude number, nose shape and skew. Generally conservative for wide piers.' },
        { key: 'scourSM', label: 'Local scour — Sheppard–Melville', value: r(s.sm), unit: 'm', explain: 'Pier scour including the sediment size and the clear-water / live-bed transition; the HEC-18 alternative for wide piers. Zero when outside the method\'s range.' },
        { key: 'scourBr', label: 'Local scour — Breusers (current)', value: r(s.br), unit: 'm', explain: 'Classical 1.5 D tanh(h/D) steady-current estimate for a cylinder, for comparison.' },
        { key: 'scourWC', label: 'Local scour — waves plus current (Sumer–Fredsøe)', value: r(s.sf), unit: 'm', explain: 'Equilibrium scour at a slender pile under combined waves and current; governs for marine towers where waves dominate (KC number small).' },
        { key: 'toeScour', label: 'Wall-toe scour (rule of thumb)', value: r(toe), unit: 'm', explain: 'Indicative maximum scour in front of a reflective wall, about one wave height. Deduct it from the toe embedment in the stability check.' },
        { key: 'shields', label: 'Shields ratio θmax / θcr', value: r(s.mob.ratio), unit: '-', explain: 'Peak wave-plus-current bed shear stress over the threshold of motion. Above 1 the bed is mobile and scour is live-bed; below 1 the far-field bed is stable.' },
        { key: 'Ucr', label: 'Threshold current speed', value: r(s.sim.Ucr, 3), unit: 'm/s', explain: 'Depth-averaged velocity at which the bed sediment starts to move. Clear-water scour is largest when the flow is near this value.' },
        { key: 'qt', label: 'Sediment transport (wave + current, Soulsby–Van Rijn)', value: r(svr * day), unit: 'm³/m/day', explain: 'Total volume of sand carried past each metre of width per day in the design conditions. It sets how fast scour holes fill, trenches infill and the bed adjusts.' },
        { key: 'tScale', label: 'Scour time scale', value: r(T / 3600, 3), unit: 'h', explain: 'Time for scour to reach about 63 % of equilibrium. Short compared with a storm means equilibrium scour applies; zero when the bed is not mobile.' },
        { key: 'D50rock', label: 'Scour-protection rock D50', value: r(D50r, 3), unit: 'm', explain: 'Median stone size stable under the accelerated flow at the pier (Isbash with the amplification factor). Use the next standard grading above.' },
        { key: 'rockT', label: 'Protection rock per pier', value: r(rockT), unit: 't', explain: 'Mass of an armour layer three stones thick over the protected radius. Add a filter layer or geotextile beneath.' },
      ],
      checks: [
        check('Scour against allowable loss of embedment', total, cap, 'm', 'FHWA HEC-18 (2012) ch. 7; project foundation criterion', 'Guards against the foundation losing the soil it relies on for lateral and axial capacity. Scour is the leading cause of bridge failure and usually happens unseen during the flood.'),
        check('Wave-plus-current scour against allowable loss of embedment', s.sf + inp.general, cap, 'm', 'Sumer & Fredsøe (2002) §3; DNV-ST-0126 (2021) §7.6', 'The same embedment criterion with the marine (oscillatory flow) scour formula. For slender towers in waves this is the relevant mechanism; for wide river piers the HEC-18 check governs.'),
      ],
      plots: [
        { type: 'line', title: 'Local scour versus current speed', x: { label: 'Current', unit: 'm/s' }, y: { label: 'Scour depth', unit: 'm' },
          series: [{ name: 'HEC-18', x: us.map((u) => r(u, 3)), y: us.map((u) => r(scourHec18({ a: inp.a, y1: h, V1: u, nose: inp.nose, skew: inp.skewDeg * DEG, L: Math.min(inp.Lp, 12 * inp.a), K3: inp.K3 }).ys, 3)) }, { name: 'Sheppard–Melville', x: us.map((u) => r(u, 3)), y: us.map((u) => r(safe(() => scourSheppardMelville({ a: s.ae, y1: h, V1: u, Vc: s.sim.Ucr, D50: d }).ys), 3)) }, { name: 'Allowable (less general scour)', x: [us[0], us[29]], y: [r(Math.max(cap - inp.general, 0), 3), r(Math.max(cap - inp.general, 0), 3)], dash: true }],
          explain: 'How scour grows with the flood or tidal current. The Sheppard–Melville curve peaks near the threshold velocity (clear-water peak) and again in strong live-bed flow; the crossing with the dashed line is the current the foundation can stand.' },
        { type: 'bar', title: 'Transport rate by formula', categories: ['MPM bed load', 'Van Rijn bed load', 'Van Rijn suspended', 'Engelund–Hansen total', 'Soulsby–Van Rijn (waves + current)'], series: [{ name: 'm³/m/day', y: [mpm, vrb, vrs, eh, svr].map((q) => r(q * day, 3)) }],
          explain: 'The same flow through five transport formulae. A spread of a factor two to five is normal — it is the honest uncertainty of sediment transport prediction and the reason to calibrate against bathymetric change.' },
        { type: 'line', title: 'Scour development in time', x: { label: 'Time', unit: 'h' }, y: { label: 'Scour depth', unit: 'm' },
          series: [(() => { const ts = linspace(0, Math.max(inp.stormH, T > 0 ? 3 * T / 3600 : inp.stormH), 40); return { name: 'S(t)', x: ts.map((t) => r(t, 3)), y: ts.map((t) => r(local * (T > 0 ? scourDevelopment(1, t * 3600, T) : (mobile ? 1 : 0)), 3)) }; })()],
          explain: 'Exponential approach to equilibrium scour. If the curve is still rising at the end of the event the equilibrium depth is conservative for a single storm but not for a sequence of them.' },
        ...(vrsO ? [{ type: 'line', title: 'Suspended-sediment concentration profile', x: { label: 'Volume concentration', unit: '-' }, y: { label: 'Height above bed', unit: 'm' },
          series: [{ name: 'Rouse profile', x: zs.map((z) => r(safe(() => rouseProfile({ z, h, a: vrsO.a, ca: vrsO.ca, ws: vrsO.ws, uStar: vrsO.uStar }).c), 3)), y: zs.map((z) => r(z, 3)) }],
          explain: 'Suspended sand concentration from the reference level to the surface. A profile hugging the bed means bed-load dominated transport; a fuller profile means sand is carried far and intakes will ingest it.' }] : []),
      ],
      tables: [
        { title: 'Bed shear stress and mobility', columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['Current shear stress', s.tc.tau, 'Pa'], ['Wave shear stress (amplitude)', s.tw, 'Pa'], ['Combined maximum', s.comb.tauMax, 'Pa'], ['Critical shear stress', s.crit.tauCr, 'Pa'], ['Critical Shields number', s.crit.thetaCr, '-'], ['Dimensionless grain size D*', s.crit.Dstar, '-'], ['Near-bed orbital velocity (rms wave)', s.orb.Uw, 'm/s'], ['Settling velocity', ws, 'm/s'], ['Keulegan–Carpenter number', s.KC, '-'], ['Current ratio Uc/(Uc+Um)', s.Ucw, '-'], ['Effective pier width', s.ae, 'm']].map(([q, v, u]) => ({ q, v: r(v), u })),
          explain: 'Intermediate hydraulic quantities behind the mobility and scour results, for checking against field measurements of bed shear stress and sediment flux.' },
        { title: 'Sea-level-rise scenarios', columns: [{ key: 'name', label: 'Scenario' }, { key: 'rise', label: 'Rise', unit: 'm' }, { key: 'Uw', label: 'Orbital velocity', unit: 'm/s' }, { key: 'theta', label: 'Shields ratio', unit: '-' }, { key: 'S', label: 'Local scour', unit: 'm' }],
          rows: scen, explain: 'Bed mobility and local scour with the water deepened by each scenario. Deeper water weakens wave stirring at the bed and slightly changes the depth term of the scour formulae.' },
      ],
      notes: [`${s.hec >= s.sm ? 'HEC-18' : 'Sheppard–Melville'} governs the local scour; design value is the larger of the two plus general scour.`, `${mobile ? 'Live-bed' : 'Clear-water'} conditions; rock sized for ${f(Uamp, 2)} m/s at the pier.`],
      advice,
      publish: { 'marine.scourDepth': total, 'marine.scourRockD50': D50r },
      quantities: { rock_t: r(ok ? 0 : rockT), vessel_days: { 'scour protection': r(ok ? 0 : days, 2) } },
    };
  },
  verify: [
    { name: 'HEC-18 CSU equation ys = 2.0 K1 K2 K3 a (y1/a)^0.35 Fr^0.43', source: 'FHWA HEC-18 (2012) Eq. 7.1. Hand calculation: a = 3 m, y1 = 8 m, V1 = 1.2 m/s, round nose, no skew, K3 = 1.1: Fr = 1.2/√(9.80665 × 8) = 0.13548, (8/3)^0.35 = 1.4096, Fr^0.43 = 0.42335 → ys = 2 × 1.1 × 1.4096 × 0.42335 × 3 = 3.939 m',
      inputs: { a: 3, Lp: 3, nose: 'round', skewDeg: 0, h: 8, slr: 0, U: 1.2, K3: 1.1 }, expect: [{ kpi: 'scourHec', value: 3.939, tol: 0.002 }] },
    { name: 'Isbash rock size with HEC-23 amplification, D50 = (K V)²/(2 g Δ C²), C = 1.2', source: 'FHWA HEC-23 (2009) DG11 Eq. 11.1 (0.692 (K V)²/((Ss − 1) 2g)). Hand calculation, no waves: K V = 1.5 × 1.2 = 1.8 m/s, Δ = 2650/1025 − 1 = 1.5854 → D50 = 3.24/(2 × 9.80665 × 1.5854 × 1.44) = 0.0724 m',
      inputs: { U: 1.2, Hs: 0, Kv: 1.5, rock: 'rock-granite', slr: 0 }, expect: [{ kpi: 'D50rock', value: 0.0724, tol: 0.005 }] },
  ],
};
