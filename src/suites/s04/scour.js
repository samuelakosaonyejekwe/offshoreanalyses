// Seabed scour and erosion around a subsea structure: sediment mobility, scour depth, loss of support and protection.
import { criticalShields, bedShearStress, pileScour, scourTimeScale, scourDevelopment, largeStructureScour, rockSizeShields } from '../../engine/geotech/scour.js';
import { soulsbyVanRijn } from '../../engine/coastal/sediment.js';
import { nearBedWaveCurrent } from '../../engine/metocean/current.js';
import { G, DEG, NU_SEA, need, run, fmt, round, adv, check, vesselCo2Text, matInputs, waterInputs } from './_common.js';

export default {
  id: 's04-scour', suite: 4, title: 'Seabed scour and erosion',
  summary: 'Sediment mobility at the structure, equilibrium scour at mat edges and legs and its development in time, the resulting loss of bearing area, and the rock size and quantity for protection.',
  covers: [12], fidelity: 'screening',
  standards: ['DNV-RP-C212 (2019) §9 (scour)', 'DNV-ST-0126 (2021) §7.7', 'DNV-RP-F109 (2021) §3 (near-bed flow)', 'Sumer & Fredsøe (2002) The Mechanics of Scour in the Marine Environment', 'Soulsby (1997) Dynamics of Marine Sands', 'CIRIA C683 The Rock Manual §5'],
  models: ['Shields critical mobility (Soulsby)', 'Bed shear stress under waves and current (Soulsby)', 'Logarithmic near-bed velocity profile', 'Keulegan–Carpenter dependent pile scour (Sumer & Fredsøe)', 'Scour time scale and exponential development', 'Edge scour of large structures (user coefficient)', 'Soulsby–van Rijn total-load sediment transport', 'Critical-shear rock sizing'],
  limits: 'Non-cohesive seabed (sand, silt); erosion of clay is strength-controlled and is not predicted. Edge scour of the mat is a user-supplied non-dimensional depth because no general formula exists; leg scour assumes slender, isolated cylinders in live-bed conditions. Morphological modelling with an evolving bed and CFD are outside the module.',
  publishes: ['subsea.scourDepth', 'subsea.scourAreaLoss'],
  searchTerms: 'scour subsea structure mudmat edge scour protection rock size Shields',
  inputs: [
    ...waterInputs(),
    { key: 'd50', label: 'Median grain size d50', unit: 'mm', type: 'number', default: 0.2, min: 0.02, max: 100, group: 'Seabed', role: 'soil', bind: 'site.d50', uncertain: { dist: 'lognormal', cov: 0.3 } },
    { key: 'rhoS', label: 'Grain density', unit: 'kg/m³', type: 'number', default: 2650, min: 1500, max: 3500, group: 'Seabed', role: 'soil' },
    { key: 'repose', label: 'Angle of repose of the scour slope', unit: 'deg', type: 'number', default: 30, min: 15, max: 45, group: 'Seabed', role: 'soil' },
    { key: 'Uref', label: 'Near-bed current speed', unit: 'm/s', type: 'number', default: 0.7, min: 0, max: 3, group: 'Flow', role: 'bc', bind: 'site.UcBed100', uncertain: { dist: 'lognormal', cov: 0.2 } },
    { key: 'zRef', label: 'Reference height of the current', unit: 'm', type: 'number', default: 1, min: 0.1, max: 20, group: 'Flow', role: 'bc' },
    { key: 'Hs', label: 'Significant wave height', unit: 'm', type: 'number', default: 13, min: 0, max: 25, group: 'Flow', role: 'bc', bind: 'site.Hs100' },
    { key: 'Tp', label: 'Peak period', unit: 's', type: 'number', default: 15.5, min: 3, max: 25, group: 'Flow', role: 'bc', bind: 'site.Tp100' },
    { key: 'angDeg', label: 'Angle between waves and current', unit: 'deg', type: 'number', default: 0, min: 0, max: 90, group: 'Flow', role: 'bc' },
    { key: 'tStorm', label: 'Duration of the condition', unit: 'h', type: 'number', default: 12, min: 0.1, max: 1e5, group: 'Flow', role: 'bc', help: 'Storm duration, or a long time for persistent currents; scour develops exponentially towards equilibrium.' },
    { key: 'S0', label: 'Scour already present', unit: 'm', type: 'number', default: 0, min: 0, max: 20, group: 'Flow', role: 'ic', help: 'Initial bed lowering at the structure from surveys; scour develops from this state.' },
    ...matInputs(),
    { key: 'D', label: 'Skirt depth', unit: 'm', type: 'number', default: 0.5, min: 0, max: 10, group: 'Mudmat geometry', role: 'geometry', design: true },
    { key: 'Dleg', label: 'Leg / guide-post diameter at the bed', unit: 'm', type: 'number', default: 0.6, min: 0.05, max: 10, group: 'Mudmat geometry', role: 'geometry' },
    { key: 'kS', label: 'Edge scour depth / mat width', unit: '-', type: 'number', default: 0.04, min: 0, max: 0.5, group: 'Scour coefficients', role: 'criteria', uncertain: { dist: 'lognormal', cov: 0.5 },
      help: 'Non-dimensional equilibrium scour at the edge of a large low structure, from model tests or published guidance (Whitehouse 1998 reports values of a few per cent of the plan dimension for gravity bases). Indicative — confirm for the configuration.' },
    { key: 'ampl', label: 'Bed shear amplification at the structure', unit: '-', type: 'number', default: 2.5, min: 1, max: 10, group: 'Scour coefficients', role: 'criteria', help: 'About 4 at a circular pile in current, 2–4 at corners of large structures (Sumer & Fredsøe 2002; Whitehouse 1998).' },
    { key: 'nSigma', label: 'Standard deviations on leg scour', unit: '-', type: 'number', default: 0, min: 0, max: 3, group: 'Scour coefficients', role: 'criteria', help: 'S/D = 1.3 + nσ·0.7 in steady current (DNV-ST-0126 quotes 1.3 D as default).' },
    { key: 'lossAllow', label: 'Allowable loss of bearing area', unit: '%', type: 'number', default: 10, min: 0, max: 60, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'wProt', label: 'Width of the protection apron', unit: 'm', type: 'number', default: 3, min: 0.5, max: 30, group: 'Protection', role: 'geometry', design: true },
    { key: 'rhoRock', label: 'Rock density', unit: 'kg/m³', type: 'number', default: 2650, min: 2000, max: 3500, group: 'Protection', role: 'material' },
  ],
  compute(inp) {
    // SI: mm → m, h → s, deg → rad.
    const d = inp.d50 / 1e3, rho = inp.rhoW, h = inp.depth, s = inp.rhoS / rho;
    need(s > 1.05, 'Grain density must exceed the water density.');
    const z0 = d / 12; // grain roughness, Soulsby (1997)
    const H = inp.Hs, T = inp.Tp;
    const Uw = H > 0 ? run(() => nearBedWaveCurrent({ H, T, d: h, zBed: 0 }), 'Wave kinematics').Uw : 0;
    // Depth-averaged current consistent with the log law through the near-bed reference value.
    need(inp.zRef > z0, 'The reference height must be above the bed roughness.');
    const Ubar = inp.Uref * (Math.log(h / z0) - 1) / Math.log(inp.zRef / z0);
    const bs = run(() => bedShearStress({ d, h, U: Ubar, Uw, T, phi: inp.angDeg * DEG, rho }), 'Bed shear stress');
    const cr = criticalShields({ d, s, nu: NU_SEA });
    const denom = (inp.rhoS - rho) * G * d;
    const theta = bs.tauMax / denom, thetaAmp = inp.ampl * theta;
    const mobFar = theta / cr.thetaCr, mobNear = thetaAmp / cr.thetaCr;
    const liveBed = mobFar > 1, scours = mobNear > 1;
    // Leg scour (Sumer & Fredsøe): KC and current ratio at the bed.
    const KC = Uw > 0 ? Uw * T / inp.Dleg : Infinity;
    const Ucw = inp.Uref + Uw > 0 ? inp.Uref / (inp.Uref + Uw) : 0;
    const legEq = scours ? pileScour({ KC, Ucw, nSigma: inp.nSigma }).SoverD * inp.Dleg : 0;
    const edge = largeStructureScour({ Lref: inp.B, kS: scours ? inp.kS : 0, KC: Uw > 0 ? Uw * T / inp.B : undefined });
    const edgeEq = edge.S;
    // Time scale on the leg diameter (current mode when current dominates, wave mode otherwise).
    let Tsc = Infinity;
    if (scours && theta > 0) {
      try { Tsc = (Ucw >= 0.5 || !Number.isFinite(KC) ? scourTimeScale({ mode: 'current', D: inp.Dleg, theta, d, s, delta: Math.min(h, 10 * inp.Dleg) }) : scourTimeScale({ mode: 'waves', D: inp.Dleg, theta, d, s, KC })).T; } catch { Tsc = Infinity; }
    }
    const t = inp.tStorm * 3600;
    const dev = (Seq, tt) => (Seq <= inp.S0 ? Math.max(Seq, inp.S0) : inp.S0 + (Number.isFinite(Tsc) ? scourDevelopment(Seq - inp.S0, tt, Tsc) : 0));
    const legS = dev(legEq, t), edgeS = dev(edgeEq, t);
    const Sgov = Math.max(edgeS, 0);
    // Undermining: bed lowered below skirt tip removes support over a band set by the scour slope.
    const under = Math.max(Sgov - inp.D, 0) / Math.tan(inp.repose * DEG);
    const Be = Math.max(inp.B - 2 * under, 0), Le = Math.max(inp.L - 2 * under, 0);
    const loss = 1 - (Be * Le) / (inp.B * inp.L);
    const ucFactor = loss < 1 ? 1 / (1 - loss) : 99;
    // Transport rate and protection.
    let qt = 0;
    try { qt = soulsbyVanRijn({ U: Ubar, Urms: Uw / Math.SQRT2, h, d50: d, rhoS: inp.rhoS, rho, nu: NU_SEA }).qt; } catch { qt = 0; }
    let rock = { d50: 0 };
    let rockNote = null;
    if (Ubar > 0 || Uw > 0) {
      try { rock = rockSizeShields({ h, U: Ubar, Uw, T, phi: inp.angDeg * DEG, rho, amplification: inp.ampl, rhoS: inp.rhoRock, nu: NU_SEA }); }
      catch (e) { rock = { d50: 0 }; rockNote = `No statically stable rock size exists for this loading (${e.message}): use a dynamically stable design, mattresses or deeper skirts.`; }
    }
    const tLayer = Math.max(3 * rock.d50, 0.3);
    const apronArea = (inp.B + 2 * inp.wProt) * (inp.L + 2 * inp.wProt) - inp.B * inp.L;
    const rockVol = apronArea * tLayer, rockT = rockVol * (1 - 0.4) * inp.rhoRock / 1e3;
    const needProt = scours && (loss * 100 > inp.lossAllow || edgeEq > inp.D);
    const checks = [
      check('Loss of bearing area by undermining', loss * 100, inp.lossAllow || 1e-9, '%', 'DNV-RP-C212 (2019) §9; ISO 19901-4:2016 §7.7 (scour allowance in foundation design)', 'Guards against scour reaching below the skirt tips and removing soil support under the mat edge, which raises contact pressure and tilts the structure.'),
      check('Edge scour against skirt depth', edgeS, Math.max(inp.D, 1e-9), 'm', 'DNV-RP-C212 (2019) §9 (skirts to extend below predicted scour)', 'Skirts must reach below the scoured bed. If the bed drops below the tip, piping under the skirt starts and base suction and sliding resistance are lost.'),
    ];
    const ok = checks[0].uc <= 1 && checks[1].uc <= 1;
    const ts = [], se = [], sl = [];
    for (let i = 0; i <= 40; i++) { const tt = t * 3 * i / 40; ts.push(round(tt / 3600)); se.push(dev(edgeEq, tt)); sl.push(dev(legEq, tt)); }
    const us = [], mob = [];
    for (let i = 0; i <= 20; i++) { const u = 0.1 * i; const ub = u * (Math.log(h / z0) - 1) / Math.log(inp.zRef / z0); const b2 = bedShearStress({ d, h, U: ub, Uw, T, phi: inp.angDeg * DEG, rho }); us.push(round(u)); mob.push(inp.ampl * b2.tauMax / denom / cr.thetaCr); }
    const advice = [
      ok ? adv('engineering', scours ? 'info' : 'ok', scours ? `The bed is mobile at the structure (θα/θcr = ${fmt(mobNear)}); predicted edge scour ${fmt(edgeS)} m stays above the ${fmt(inp.D)} m skirt tip and bearing-area loss is ${fmt(loss * 100)} %. Leg scour ${fmt(legS)} m.`
        : `The bed is not mobile even with the flow amplification at the structure (θα/θcr = ${fmt(mobNear)}): no scour is predicted for this condition.`,
        scours ? `Include ${fmt(edgeS)} m of scour in the sliding and bearing checks, and plan a post-storm survey; the equilibrium value ${fmt(edgeEq)} m is approached after about ${Number.isFinite(Tsc) ? fmt(3 * Tsc / 3600) : '–'} h.` : 'Confirm the grain size and the long-term current statistics; re-run for the smaller grain fraction and for persistent tidal or loop currents.',
        'DNV-RP-C212 (2019) §9; DNV-ST-0126 (2021) §7.7; Sumer & Fredsøe (2002)',
        'No rock is required, avoiding quarrying and placement: ' + fmt(rockT) + ' t of rock and its vessel time would otherwise be added.')
        : adv('engineering', 'fail', `Edge scour ${fmt(edgeS)} m against ${fmt(inp.D)} m of skirt removes ${fmt(loss * 100)} % of the bearing area (limit ${inp.lossAllow} %), raising bearing utilisation by a factor of ${fmt(ucFactor)}.`,
          `Extend the skirts to at least ${fmt(edgeEq * 1.2)} m, or place a rock apron ${fmt(inp.wProt)} m wide of d50 ≥ ${fmt(rock.d50 * 1e3)} mm, ${fmt(tLayer)} m thick (${fmt(rockT)} t) on a filter layer.`,
          'DNV-RP-C212 (2019) §9; CIRIA C683 Rock Manual §5.2 and §5.4; Sumer & Fredsøe (2002)',
          `Deeper skirts add about ${fmt(2 * (inp.B + inp.L) * Math.max(edgeEq * 1.2 - inp.D, 0) * 0.02 * 7.85)} t of steel; the rock alternative is ${fmt(rockT)} t placed over ${fmt(apronArea)} m² of seabed habitat, plus ${vesselCo2Text(1.5)}`),
      adv('economic', needProt ? 'warn' : 'info', needProt ? `Protection would need ${fmt(rockT)} t of rock (${fmt(rockVol)} m³ bulk) and a rock-placement vessel for about 1–2 days.` : `No protection is needed now; a protection apron, if later required, would be ${fmt(rockT)} t of rock.`,
        needProt ? 'Compare skirt extension at fabrication (steel only, no offshore time) with rock placement after installation; skirts are usually cheaper when decided before fabrication.' : 'Adopt monitor-and-react: include bed-level measurement in routine ROV inspections and trigger protection only if scour exceeds half the skirt depth.',
        'Risk-based inspection planning per DNV-RP-G101; whole-life costing per ISO 15686-5'),
    ];
    return {
      kpis: [
        { key: 'edgeS', label: 'Scour depth at the mat edge', unit: 'm', value: edgeS, limit: inp.D, sense: 'max', explain: 'Predicted bed lowering along the mat edge after the stated duration. It should stay above the skirt tip.' },
        { key: 'edgeEq', label: 'Equilibrium edge scour', unit: 'm', value: edgeEq, explain: 'Edge scour approached if the condition persists: user coefficient × mat width when the bed is mobile.' },
        { key: 'legS', label: 'Scour depth at legs / posts', unit: 'm', value: legS, explain: 'Local scour around slender vertical members after the stated duration (Sumer & Fredsøe).' },
        { key: 'legEq', label: 'Equilibrium leg scour', unit: 'm', value: legEq, explain: 'Equilibrium local scour at a leg: 1.3 diameters in steady current, less in wave-dominated flow at low KC.' },
        { key: 'thetaCr', label: 'Critical Shields parameter', unit: '-', value: cr.thetaCr, explain: 'Threshold of grain motion for this sediment (Soulsby 1997). Fine sand is about 0.05; coarse material lower in terms of D*.' },
        { key: 'theta', label: 'Shields parameter (undisturbed)', unit: '-', value: theta, explain: 'Maximum wave-current bed shear stress made dimensionless with grain weight, away from the structure.' },
        { key: 'mobility', label: 'Mobility ratio at the structure', unit: '-', value: mobNear, explain: 'Amplified Shields parameter over critical. Above 1 scour develops (clear-water if the far-field ratio is below 1, live-bed otherwise); it is a trigger for the scour assessment, not a failure in itself.' },
        { key: 'areaLoss', label: 'Loss of bearing area', unit: '%', value: loss * 100, limit: inp.lossAllow, sense: 'max', explain: 'Share of the base left unsupported once the scour slope undercuts the skirt tips.' },
        { key: 'ucFactor', label: 'Increase factor on bearing utilisation', unit: '-', value: ucFactor, explain: 'Multiplier on contact pressure and bearing utilisation caused by the lost area. Apply it to the bearing-capacity result.' },
        { key: 'Tscour', label: 'Scour time scale', unit: 'h', value: Number.isFinite(Tsc) ? Tsc / 3600 : 0, explain: 'Time for scour to reach 63 % of equilibrium. Short compared with a storm means equilibrium is reached; zero when the bed is immobile.' },
        { key: 'rockD50', label: 'Stable rock size d50', unit: 'mm', value: rock.d50 * 1e3, explain: 'Smallest stone that stays put under the amplified bed shear stress (critical-shear criterion). Zero when no statically stable size exists for the loading.' },
        { key: 'rockT', label: 'Rock for a protection apron', unit: 't', value: rockT, explain: 'Mass of an apron of the given width and three stone diameters (at least 0.3 m) thick, 40 % porosity.' },
        { key: 'qt', label: 'Sediment transport rate', unit: 'm³/m/day', value: qt * 86400, explain: 'Total-load transport per metre width in the undisturbed flow (Soulsby–van Rijn). A high rate means scour holes form — and backfill — quickly.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Scour development in time', x: { label: 'Time', unit: 'h' }, y: { label: 'Scour depth', unit: 'm' },
          series: [{ name: 'Mat edge', x: ts, y: se }, { name: 'Leg / post', x: ts, y: sl }, { name: 'Skirt tip', x: [0, ts[ts.length - 1]], y: [inp.D, inp.D], dash: true }],
          explain: 'Growth of scour towards equilibrium from the present bed level. The edge curve must stay below the dashed skirt-tip line for the duration considered.' },
        { type: 'line', title: 'Mobility at the structure against current speed', x: { label: 'Near-bed current', unit: 'm/s' }, y: { label: 'θ·α / θcr', unit: '-' },
          series: [{ name: 'Mobility ratio', x: us, y: mob }, { name: 'Threshold', x: [0, 2], y: [1, 1], dash: true }, { name: 'Design condition', x: [inp.Uref], y: [mobNear] }],
          explain: 'How the mobility ratio grows with current at the present wave condition. The crossing with 1.0 is the current speed at which scour starts.' },
      ],
      tables: [{ title: 'Flow and sediment parameters', columns: [{ key: 'k', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
        rows: [['Depth-averaged current', Ubar, 'm/s'], ['Wave orbital velocity at the bed', Uw, 'm/s'], ['Current shear stress', bs.tauC, 'Pa'], ['Wave shear stress', bs.tauW, 'Pa'], ['Maximum combined shear stress', bs.tauMax, 'Pa'], ['Dimensionless grain size D*', cr.Dstar, '-'],
          ['Far-field mobility θ/θcr', mobFar, '-'], ['KC on the leg', Number.isFinite(KC) ? KC : 0, '-'], ['Current ratio Uc/(Uc+Uw)', Ucw, '-'], ['Edge scour regime', edge.regime, ''], ['Undermined band width', under, 'm'], ['Effective mat after scour', `${fmt(Be)} × ${fmt(Le)}`, 'm'],
          ['Protection layer thickness', tLayer, 'm'], ['Apron area', apronArea, 'm²']].map(([k, v, u]) => ({ k, v: typeof v === 'number' ? round(v) : v, u })),
        explain: 'Intermediate hydraulic and sediment quantities behind the scour prediction, for checking against a hand calculation or a survey.' }],
      notes: [scours ? (liveBed ? 'Live-bed conditions: the whole seabed is mobile.' : 'Clear-water scour: the bed is mobile only where the structure amplifies the flow.') : 'No sediment motion is predicted at the structure.', `Edge scour uses the coefficient kS = ${inp.kS} on the mat width (indicative).`, `Scour ${edgeS > inp.D ? 'reaches below' : 'stays above'} the skirt tip.`, ...(rockNote ? [rockNote] : [])],
      advice,
      publish: { 'subsea.scourDepth': edgeS, 'subsea.scourAreaLoss': loss },
      quantities: needProt ? { rock_t: round(rockT), vessel_days: { rock_placement: 1.5 } } : {},
    };
  },
  verify: [
    { name: 'Critical Shields parameter of 0.2 mm sand: θcr = 0.30/(1+1.2D*) + 0.055(1 − e^(−0.02D*))', source: 'Soulsby (1997) eq. 77 with D* = d[g(s−1)/ν²]^(1/3) = 4.121 (s = 2.65, ν = 1.36e-6): θcr = 0.05046 + 0.00435 = 0.0548',
      inputs: { d50: 0.2, rhoS: 2650, rhoW: 1000 }, expect: [{ kpi: 'thetaCr', value: 0.05481, tol: 0.002 }] },
    { name: 'Leg scour in steady current: S = 1.3·D', source: 'Sumer, Fredsøe & Christiansen (1992); DNV-ST-0126 §7.7: 1.3 × 0.5 m = 0.65 m',
      inputs: { Hs: 0, Uref: 1.0, Dleg: 0.5, nSigma: 0, d50: 0.2 }, expect: [{ kpi: 'legEq', value: 0.65, tol: 1e-6 }] },
  ],
};
