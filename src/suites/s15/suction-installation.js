// Suction foundation installation: self-weight penetration, required suction against depth, and the limits set by
// plug heave or piping, cavitation, pump capacity and shell buckling.
import { suctionPenetrationCurve, suctionInstallClay, suctionInstallSand, ringBucklingPressure } from '../../engine/geotech/suction.js';
import { G, DEG, need, run, round, fmt, kpi, check, adv, line, tbl, worst, steelCo2Text, vesselCo2Text, CONFIRM } from './_common.js';

export default {
  id: 's15-suction-installation', suite: 15, title: 'Suction foundation installation',
  summary: 'Predicts how far a suction bucket or caisson sinks under its own weight, the suction needed to reach full penetration, and whether plug heave or piping, cavitation, pump capacity or shell buckling stops the installation.',
  covers: [13], fidelity: 'analytical',
  standards: ['DNV-RP-E303 (2021) §5 (installation of suction anchors in clay)', 'DNV-ST-0126 (2021) §7.9 (suction bucket installation)', 'Houlsby & Byrne (2005) Design procedures for installation of suction caissons in clay / in sand, Proc. ICE Geotech. Eng. 158', 'DNV-RP-C202 (2021) (shell buckling — full check outside this module)'],
  models: ['Suction penetration resistance (wall friction + tip bearing)', 'Undrained soil strength profile with adhesion factor', 'Effective stress analysis with wall-friction stress enhancement in sand', 'Soil seepage: reduced tip resistance and critical (piping) suction in sand', 'Plug heave (reverse end bearing) in clay', 'Cavitation limit of the pump head', 'Elastic ring buckling under suction'],
  limits: 'Uniform or linearly increasing clay strength, or uniform sand, without layering, boulders, stiffeners or internal ring plates; installation coefficients (adhesion, K·tan δ, bearing factors, seepage factor) are literature defaults to be confirmed for the site. Buckling is the long-cylinder elastic ring formula, a lower-bound screen only; coupled seepage finite elements and large-deformation penetration analysis are outside the module.',
  publishes: ['install.suctionRequired', 'install.suctionSelfWeightPen'],
  searchTerms: 'suction caisson bucket installation required suction plug heave piping self-weight penetration',
  inputs: [
    { key: 'Do', label: 'Outer diameter', unit: 'm', type: 'number', default: 6, min: 1, max: 40, group: 'Foundation', role: 'geometry', design: true },
    { key: 'L', label: 'Skirt length (target penetration)', unit: 'm', type: 'number', default: 14, min: 0.5, max: 60, group: 'Foundation', role: 'bc', design: true, help: 'Boundary condition of the installation: the skirt must reach this depth for the in-place capacity.' },
    { key: 't', label: 'Skirt wall thickness', unit: 'mm', type: 'number', default: 50, min: 8, max: 150, group: 'Foundation', role: 'geometry', design: true },
    { key: 'Vp', label: 'Submerged weight incl. ballast and attached structure', unit: 'kN', type: 'number', default: 1500, min: 0, max: 2e5, group: 'Foundation', role: 'load', uncertain: { dist: 'normal', cov: 0.05 } },
    { key: 'mat', label: 'Skirt steel', type: 'select', source: 'materials:steel', default: 's355', group: 'Foundation', role: 'material' },
    { key: 'soil', label: 'Soil', type: 'select', default: 'clay', group: 'Soil', role: 'soil', options: [{ value: 'clay', label: 'Clay (undrained)' }, { value: 'sand', label: 'Sand (drained with seepage)' }] },
    { key: 'su0', label: 'Undrained shear strength at mudline', unit: 'kPa', type: 'number', default: 5, min: 0.2, max: 400, group: 'Soil', role: 'soil', bind: 'site.su0', uncertain: { dist: 'lognormal', cov: 0.25 }, showIf: { key: 'soil', in: ['clay'] } },
    { key: 'suGrad', label: 'Strength gradient', unit: 'kPa/m', type: 'number', default: 1.6, min: 0, max: 50, group: 'Soil', role: 'soil', bind: 'site.suGrad', uncertain: { dist: 'normal', cov: 0.2 }, showIf: { key: 'soil', in: ['clay'] } },
    { key: 'alpha', label: 'Wall adhesion factor during installation', unit: '-', type: 'number', default: 0.4, min: 0.1, max: 1, group: 'Soil', role: 'soil', uncertain: { dist: 'normal', cov: 0.2 }, showIf: { key: 'soil', in: ['clay'] }, help: `About 1/sensitivity: 0.25–0.5 for soft clays, lower with painted steel. ${CONFIRM}` },
    { key: 'Nc', label: 'Tip bearing factor Nc', unit: '-', type: 'number', default: 7.5, min: 5, max: 12, group: 'Soil', role: 'soil', showIf: { key: 'soil', in: ['clay'] }, help: `7.5 for a strip-like skirt tip (DNV-RP-E303), 9 for deep bearing. ${CONFIRM}` },
    { key: 'gammaSub', label: 'Submerged unit weight γ′', unit: 'kN/m³', type: 'number', default: 6, min: 2, max: 13, group: 'Soil', role: 'soil', bind: 'site.gammaSub', uncertain: { dist: 'normal', cov: 0.05 } },
    { key: 'phiDeg', label: 'Friction angle φ′', unit: 'deg', type: 'number', default: 36, min: 25, max: 45, group: 'Soil', role: 'soil', bind: 'site.phiDeg', uncertain: { dist: 'normal', cov: 0.05 }, showIf: { key: 'soil', in: ['sand'] } },
    { key: 'KtanD', label: 'K·tan δ on the skirt wall', unit: '-', type: 'number', default: 0.5, min: 0.1, max: 1.5, group: 'Soil', role: 'soil', showIf: { key: 'soil', in: ['sand'] }, help: `Lateral stress ratio times wall friction; about 0.4–0.8 for dense sand. ${CONFIRM}` },
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 300, min: 5, max: 3500, group: 'Site', role: 'environment', bind: 'site.depth' },
    { key: 'pumpDp', label: 'Pump maximum differential pressure', unit: 'kPa', type: 'number', default: 400, min: 10, max: 5000, group: 'Pump', role: 'bc', design: true },
    { key: 'pumpQ', label: 'Pump flow rate', unit: 'm³/h', type: 'number', default: 150, min: 1, max: 5000, group: 'Pump', role: 'bc' },
    { key: 'fsSoil', label: 'Safety factor on plug heave / piping', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Criteria', role: 'criteria', help: `Allowable suction = critical suction / this factor. ${CONFIRM}` },
    { key: 'fsBuck', label: 'Safety factor on shell buckling', unit: '-', type: 'number', default: 2.0, min: 1, max: 5, group: 'Criteria', role: 'criteria', help: 'Applied to the elastic ring buckling pressure, which ignores imperfections; confirm with a DNV-RP-C202 shell check.' },
    { key: 'swpMin', label: 'Minimum self-weight penetration for a seal', unit: 'm', type: 'number', default: 1.0, min: 0.1, max: 10, group: 'Criteria', role: 'criteria', help: 'The skirt tip must be embedded enough to hold suction without piping around the rim when pumping starts.' },
    { key: 'vLower', label: 'Lowering speed during self-weight penetration', unit: 'm/min', type: 'number', default: 0.5, min: 0.05, max: 10, group: 'Pump', role: 'ic', help: 'Initial condition: touchdown and self-weight penetration under crane control.' },
  ],
  compute(inp, ctx) {
    // mm → m, kPa → Pa, kN → N, kN/m³ → N/m³.
    const steel = ctx.material(inp.mat); need(steel, `Steel '${inp.mat}' is not in the materials library.`);
    const Do = inp.Do, t = inp.t / 1e3, L = inp.L, Vp = inp.Vp * 1e3, gam = inp.gammaSub * 1e3, clay = inp.soil === 'clay';
    need(t < Do / 4, 'The wall is too thick for a skirt of this diameter.');
    let params;
    if (clay) params = { Do, t, su: { su0: inp.su0 * 1e3, k: inp.suGrad * 1e3 }, alphaO: inp.alpha, alphaI: inp.alpha, Nc: inp.Nc, NcPlug: 9, gamma: gam, Vprime: Vp };
    else {
      const phi = inp.phiDeg * DEG, Nq = Math.exp(Math.PI * Math.tan(phi)) * Math.tan(Math.PI / 4 + phi / 2) ** 2, Ng = 1.5 * (Nq - 1) * Math.tan(phi);
      params = { Do, t, gamma: gam, KtanDeltaO: inp.KtanD, KtanDeltaI: inp.KtanD, Nq, Ngamma: Ng, Vprime: Vp };
    }
    const n = 56;
    const cv = run(() => suctionPenetrationCurve({ soil: inp.soil, L, n, params }), 'Suction penetration curve');
    const end = run(() => (clay ? suctionInstallClay({ ...params, h: L }) : suctionInstallSand({ ...params, h: L })), 'Suction at full penetration');
    const sReq = Math.max(end.s, 0), sCrit = clay ? end.sAllow : end.sCrit, R = clay ? end.R : end.R0;
    const swp = Math.min(cv.selfWeight, L);
    const sCav = 1025 * G * inp.depth + 101325 - 2300; // absolute pressure available before the water under the lid cavitates
    const pBuck = ringBucklingPressure({ E: steel.E || 2.1e11, t, D: Do - t });
    const sReqMax = Math.max(...cv.s, 0);
    // Deepest penetration reachable within every limit.
    const lim = (i) => Math.min(cv.sLimit[i] / inp.fsSoil, sCav, inp.pumpDp * 1e3, pBuck / inp.fsBuck);
    let hMax = swp; for (let i = 0; i < cv.h.length; i++) { if (cv.s[i] <= lim(i)) hMax = Math.max(hMax, cv.h[i]); else if (cv.h[i] > swp) break; }
    const Ain = Math.PI / 4 * (Do - 2 * t) ** 2;
    const hours = swp / inp.vLower / 60 + Math.max(L - swp, 0) * Ain / inp.pumpQ;
    const checks = [
      check(clay ? 'Required suction against plug heave' : 'Required suction against piping', sReq / 1e3, sCrit / inp.fsSoil / 1e3, 'kPa', clay ? 'DNV-RP-E303 (2021) §5 (allowable underpressure: reverse end bearing of the soil plug); Houlsby & Byrne (2005a)' : 'Houlsby & Byrne (2005b) critical suction γ′·h/(1 − a); DNV-ST-0126 §7.9', clay ? 'Too much suction lifts the soil plug inside the bucket instead of pulling the skirt down: penetration stops short and the lid lands on heaved soil.' : 'Too much suction liquefies the sand inside (piping): the plug loosens, heaves and the seal is lost.'),
      check('Required suction against cavitation', sReqMax / 1e3, sCav / 1e3, 'kPa', 'Hydrostatic plus atmospheric pressure less vapour pressure at the lid', 'The pump cannot lower the pressure under the lid below vapour pressure. In shallow water this caps the available suction.'),
      check('Required suction against pump capacity', sReqMax / 1e3, inp.pumpDp, 'kPa', 'Pump skid rating; DNV-ST-N001 (2021) Sec. 18 (suction installation equipment)', 'The pump must deliver the largest suction needed over the whole penetration.'),
      check('Suction against skirt buckling', sReqMax / 1e3, pBuck / inp.fsBuck / 1e3, 'kPa', 'Bresse–Bryan ring buckling 2E/(1 − ν²)·(t/D)³ as a screen; DNV-RP-C202 (2021) for the design check', 'External pressure on the unsupported skirt above the mudline can collapse the shell. Thin, large-diameter buckets are governed by this.'),
      check('Self-weight penetration for an initial seal', inp.swpMin, Math.max(swp, 1e-6), 'm (required / achieved)', 'DNV-RP-E303 (2021) §5; DNV-ST-0126 §7.9 (seal before pumping)', 'Without enough initial embedment water flows under the rim when the pump starts and no differential pressure can build.'),
    ];
    const gov = worst(checks), ok = gov.uc <= 1;
    const steel_t = (Math.PI * (Do - t) * t * L + Math.PI / 4 * Do * Do * 1.5 * t) * 7.85;
    const tBuck = t * (sReqMax * inp.fsBuck / pBuck) ** (1 / 3) * 1e3;
    return {
      kpis: [
        kpi('swp', 'Self-weight penetration', swp, 'm', 'Depth at which soil resistance equals the submerged weight: where pumping has to take over.', inp.swpMin, 'min'),
        kpi('R', 'Penetration resistance at full depth', R / 1e3, 'kN', clay ? 'Wall adhesion inside and outside plus tip bearing at the target penetration.' : 'Wall friction and tip resistance at the target penetration without suction (seepage reduces it during pumping).'),
        kpi('sReq', 'Required suction at full penetration', sReq / 1e3, 'kPa', 'Differential pressure across the lid needed to complete the penetration: (resistance − weight) over the lid area, with the seepage relief in sand.', sCrit / inp.fsSoil / 1e3, 'max'),
        kpi('sCrit', clay ? 'Critical suction for plug heave at full depth' : 'Critical suction for piping at full depth', sCrit / 1e3, 'kPa', 'Suction at which the soil plug fails. The allowable value is this divided by the safety factor.'),
        kpi('fsSoil', 'Safety factor against plug failure', sReq > 0 ? Math.min(sCrit / sReq, 99) : 99, '-', 'Critical over required suction at full penetration; 99 when no suction is needed.', inp.fsSoil, 'min'),
        kpi('sCav', 'Cavitation limit', sCav / 1e3, 'kPa', 'Largest differential pressure physically available at this water depth.'),
        kpi('pBuck', 'Elastic ring buckling pressure of the skirt', pBuck / 1e3, 'kPa', 'Lower-bound collapse pressure of the unstiffened shell, before the safety factor.'),
        kpi('hMax', 'Deepest penetration within all limits', hMax, 'm', 'Greatest depth at which the required suction stays below the soil, cavitation, pump and buckling limits.', L, 'min'),
        kpi('hours', 'Installation time (self-weight + pumping)', hours, 'h', 'Lowering through the self-weight phase plus the water volume to pump out at the pump flow rate.'),
        kpi('steel', 'Bucket steel (skirt and lid, approximate)', steel_t, 't', 'Skirt shell plus a lid of 1.5 × the skirt thickness; stiffeners excluded.'),
      ],
      checks,
      plots: [
        line('Suction against penetration', 'Suction', 'kPa', 'Penetration', 'm', [{ name: 'Required', x: cv.s.map((v) => round(Math.max(v, 0) / 1e3)), y: cv.h.map((v) => -round(v)) }, { name: `Allowable (${clay ? 'plug heave' : 'piping'} / ${inp.fsSoil})`, x: cv.sLimit.map((v) => round(v / inp.fsSoil / 1e3)), y: cv.h.map((v) => -round(v)), dash: true }, { name: 'Pump capacity', x: [inp.pumpDp, inp.pumpDp], y: [0, -L], dash: true }, { name: 'Buckling / safety factor', x: [pBuck / inp.fsBuck / 1e3, pBuck / inp.fsBuck / 1e3], y: [0, -L], dash: true }],
          'Required suction is zero down to the self-weight penetration and then grows. Installation is feasible while the required curve stays left of every limit; where it crosses one, penetration stops.'),
      ],
      tables: [tbl('Installation resistance and limits at the target penetration', [['k', 'Quantity'], ['v', 'Value'], ['u', 'Unit']], (clay ? [
        ['Outside wall adhesion', end.outside / 1e3, 'kN'], ['Inside wall adhesion', end.inside / 1e3, 'kN'], ['Tip bearing', end.tip / 1e3, 'kN'], ['Total resistance', R / 1e3, 'kN'], ['Submerged weight', Vp / 1e3, 'kN'], ['Required suction', sReq / 1e3, 'kPa'], ['Allowable suction', sCrit / inp.fsSoil / 1e3, 'kPa'],
      ] : [
        ['Outside wall friction', end.outside / 1e3, 'kN'], ['Inside wall friction', end.inside / 1e3, 'kN'], ['Tip resistance', end.tip / 1e3, 'kN'], ['Self-weight resistance', R / 1e3, 'kN'], ['Pore-pressure factor a', end.a, '-'], ['Required suction', sReq / 1e3, 'kPa'], ['Critical suction', sCrit / 1e3, 'kPa'],
      ]).map(([k, v, u]) => ({ k, v: round(v), u })), 'Breakdown of the soil resistance at full penetration and the suction balance, for comparison with the installation monitoring (pressure against depth).')],
      notes: [`Self-weight penetration ${fmt(swp)} m of ${L} m; required suction ${fmt(sReq / 1e3)} kPa at full depth against ${fmt(sCrit / inp.fsSoil / 1e3)} kPa allowable.`, `Cavitation limit ${fmt(sCav / 1e3)} kPa, pump ${inp.pumpDp} kPa, buckling ${fmt(pBuck / inp.fsBuck / 1e3)} kPa.`, `"${gov.name}" governs at ${fmt(gov.uc)}.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `The bucket reaches ${L} m: it sinks ${fmt(swp)} m under its own weight and needs ${fmt(sReq / 1e3)} kPa of suction, with a factor of ${fmt(sReq > 0 ? sCrit / sReq : 99)} against ${clay ? 'plug heave' : 'piping'}.`,
          `Repeat with upper- and lower-bound soil strength (the required suction scales almost linearly with ${clay ? 'the adhesion factor and su' : 'K·tan δ'}), and confirm the skirt with a DNV-RP-C202 buckling check at ${fmt(sReqMax / 1e3)} kPa including fabrication tolerances.`, clay ? 'DNV-RP-E303 (2021) §5' : 'DNV-ST-0126 (2021) §7.9; Houlsby & Byrne (2005b)')
          : adv('engineering', 'fail', `"${gov.name}" is not met (${fmt(gov.uc)}): required suction ${fmt(sReqMax / 1e3)} kPa; the bucket can be brought to ${fmt(hMax)} m of ${L} m.`,
            gov.name.includes('plug') || gov.name.includes('piping') ? `Reduce the aspect ratio (larger diameter, shorter skirt), add ballast to cut the suction demand, or accept ${fmt(hMax)} m and re-check the in-place capacity.` : gov.name.includes('buckling') ? `Increase the skirt wall to about ${fmt(Math.ceil(tBuck))} mm or add ring stiffeners.` : gov.name.includes('pump') ? `Use a pump rated for at least ${fmt(sReqMax / 1e3)} kPa.` : gov.name.includes('cavitation') ? 'Add ballast weight or reduce resistance (smaller wall, skirt-tip water injection): the water depth cannot supply more suction.' : 'Add weight or reduce the tip area so the skirt embeds enough to seal before pumping.',
            clay ? 'DNV-RP-E303 (2021) §5' : 'DNV-ST-0126 (2021) §7.9; Houlsby & Byrne (2005b)'),
        adv('economic', 'info', `Installation takes about ${fmt(hours)} h per bucket with a ${inp.pumpQ} m³/h pump and no hammer; the bucket contains about ${fmt(steel_t)} t of steel.`,
          'Suction installation avoids the piling spread and its noise mitigation: compare the heavier fabricated bucket against driven piles on vessel time per foundation, and note that the bucket can be removed by reversing the pump at end of life.',
          'Foundation concept comparison on installed cost per location (fabrication tonnes + vessel hours); AACE RP 18R-97 estimate classes',
          `Bucket steel ≈ ${steelCo2Text(steel_t)} Installation vessel time ≈ ${vesselCo2Text(hours / 24, 35)} No impact-piling noise is generated.`),
      ],
      publish: { 'install.suctionRequired': sReq, 'install.suctionSelfWeightPen': swp },
      quantities: { steel_t: round(steel_t), vessel_days: { suctionInstallation: round(hours / 24, 3) } },
    };
  },
  verify: [
    { name: 'Clay, uniform strength: R = α·su·π·(Do + Di)·h + (γ′·h + Nc·su)·π·D·t, s = (R − V′)/(π·Do²/4)', source: 'Houlsby & Byrne (2005a) hand calculation: Do = 6 m, t = 40 mm, h = 14 m, su = 20 kPa uniform, α = 0.4, Nc = 7.5, γ′ = 6 kN/m³, V′ = 1500 kN → adhesion 0.4·20·π·11.92·14 = 4194.1 kN, tip (84 + 150)·π·5.96·0.04 = 175.3 kN, R = 4369.3 kN, s = 2869.3/28.274 = 101.5 kPa',
      inputs: { su0: 20, suGrad: 0, t: 40 }, expect: [{ kpi: 'R', value: 0.4 * 20 * Math.PI * 11.92 * 14 + (84 + 150) * Math.PI * 5.96 * 0.04, tol: 0.002 }, { kpi: 'sReq', value: (0.4 * 20 * Math.PI * 11.92 * 14 + (84 + 150) * Math.PI * 5.96 * 0.04 - 1500) / (Math.PI * 9), tol: 0.01 }] },
    { name: 'Cavitation limit: ρ·g·h + p_atm − p_vapour', source: 'Hydrostatics: 1025·9.80665·300 + 101 325 − 2300 = 3114.6 kPa', inputs: {}, expect: [{ kpi: 'sCav', value: (1025 * 9.80665 * 300 + 101325 - 2300) / 1e3, tol: 1e-9 }] },
    { name: 'Ring buckling: p = 2E/(1 − ν²)·(t/D)³', source: 'Bresse–Bryan (Timoshenko & Gere, Theory of Elastic Stability §7.2): E = 210 GPa, t = 40 mm, mean D = 5.96 m → 2·210e9/0.91·(0.04/5.96)³ = 139.5 kPa', inputs: { t: 40 }, expect: [{ kpi: 'pBuck', value: 2 * 210e9 / 0.91 * (0.04 / 5.96) ** 3 / 1e3, tol: 1e-6 }] },
  ],
};
