// Suction-bucket foundation: self-weight penetration, required suction against penetration with its limits
// (plug heave in clay, piping in sand, cavitation, skirt buckling), pump flow and installation time, and the
// in-place tensile and moment capacity.
import { suctionInstallClay, suctionInstallSand, selfWeightPenetrationClay, selfWeightPenetrationSand, pumpFlow, ringBucklingPressure, caissonVerticalCapacityClay, tensileCapacitySand, bucketMomentCapacity } from '../../engine/geotech/suction.js';
import { bearingNq, bearingNgamma } from '../../engine/geotech/bearing.js';
import { c202ElasticBuckling } from '../../engine/structural/shells.js';
import { G, DEG, need, run, fmt, round, adv, check, worst, sev, steelCo2Text, steelInputs, steelOf, LIMIT_CODE, LIMIT_HIFI, CONFIRM } from './_common.js';

export default {
  id: 's03-suction-bucket-installation', suite: 3, title: 'Suction-bucket installation and capacity',
  summary: 'Penetration resistance and required suction of a suction bucket with depth, the limits set by plug heave or piping, cavitation and skirt buckling, pump flow and installation time, and the in-place tensile and moment capacity.',
  covers: [14, 8], fidelity: 'analytical',
  standards: ['DNV-ST-0126 (2021) §7.4 and §8 (suction bucket foundations, installation)', 'DNV-RP-E303 (2021) suction anchors in clay', 'Offshore Wind Accelerator (2019) Suction Installed Caisson Foundations — design guidelines', 'Houlsby & Byrne (2005a, b) Proc. ICE Geotechnical Engineering 158', 'DNV-RP-C202 (2019) §3.4 (external pressure)', 'Houlsby, Kelly & Byrne (2005) tensile capacity'],
  models: ['Penetration resistance from wall friction and tip bearing with effective-stress enhancement (Houlsby & Byrne)', 'Seepage-induced reduction of tip resistance and critical suction for piping in sand', 'Reverse end-bearing limit on suction (plug heave) in clay', 'Cavitation limit from water depth', 'Shell buckling of the skirt under external pressure', 'Darcy seepage and displaced-water pump flow', 'Uplift resistance: drained friction and rapid (suction-assisted) capacity; reverse end bearing in clay', 'Rigid-rotation moment capacity of the bucket'],
  limits: `Uniform sand or clay with linear strength increase; layered profiles, boulders, lid-stiffener resistance, skirt-tip buckling on stones and cyclic degradation of the in-place capacity are not modelled. The moment capacity is a screening value without base resistance. ${LIMIT_CODE} ${LIMIT_HIFI}`,
  searchTerms: 'suction bucket caisson installation required suction piping plug heave offshore wind',
  inputs: [
    { key: 'Do', label: 'Bucket outer diameter', unit: 'm', type: 'number', default: 10, min: 2, max: 40, group: 'Bucket geometry', role: 'geometry', design: true },
    { key: 'L', label: 'Skirt length', unit: 'm', type: 'number', default: 6, min: 1, max: 40, group: 'Bucket geometry', role: 'geometry', design: true, help: 'Target penetration of the skirt tip (installation end state).' },
    { key: 'tw', label: 'Skirt wall thickness', unit: 'mm', type: 'number', default: 30, min: 8, max: 120, group: 'Bucket geometry', role: 'geometry', design: true },
    ...steelInputs(),
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 35, min: 5, max: 300, group: 'Site', role: 'environment', bind: 'site.depth' },
    { key: 'soil', label: 'Soil', type: 'select', default: 'sand', group: 'Soil', role: 'soil', options: [{ value: 'sand', label: 'Sand (drained, seepage)' }, { value: 'clay', label: 'Clay (undrained)' }] },
    { key: 'phiDeg', label: 'Effective friction angle φ′', unit: 'deg', type: 'number', default: 36, min: 25, max: 45, group: 'Soil', role: 'soil', bind: 'site.phiDeg', showIf: { key: 'soil', in: ['sand'] }, uncertain: { dist: 'normal', cov: 0.06 } },
    { key: 'Ktan', label: 'Wall friction K·tanδ', unit: '-', type: 'number', default: 0.5, min: 0.1, max: 1.2, group: 'Soil', role: 'soil', showIf: { key: 'soil', in: ['sand'] }, uncertain: { dist: 'lognormal', cov: 0.25 }, help: 'Lateral stress coefficient times interface friction; back-calculated values are 0.3–0.9 (Houlsby & Byrne 2005b).' },
    { key: 'kPerm', label: 'Permeability of the sand', unit: 'm/s', type: 'number', default: 1e-4, min: 1e-7, max: 1e-2, group: 'Soil', role: 'soil', showIf: { key: 'soil', in: ['sand'] } },
    { key: 'kf', label: 'Permeability ratio inside / outside', unit: '-', type: 'number', default: 3, min: 1, max: 10, group: 'Soil', role: 'soil', showIf: { key: 'soil', in: ['sand'] }, help: 'Loosening of the plug under upward seepage raises the inside permeability and the piping limit (2–5 back-calculated).' },
    { key: 'su0', label: 'Undrained shear strength at mudline', unit: 'kPa', type: 'number', default: 15, min: 1, max: 400, group: 'Soil', role: 'soil', bind: 'site.su0', showIf: { key: 'soil', in: ['clay'] }, uncertain: { dist: 'lognormal', cov: 0.25 } },
    { key: 'suGrad', label: 'Strength gradient', unit: 'kPa/m', type: 'number', default: 2, min: 0, max: 30, group: 'Soil', role: 'soil', bind: 'site.suGrad', showIf: { key: 'soil', in: ['clay'] } },
    { key: 'alpha', label: 'Adhesion factor during installation', unit: '-', type: 'number', default: 0.5, min: 0.1, max: 1, group: 'Soil', role: 'soil', showIf: { key: 'soil', in: ['clay'] }, help: 'About the inverse of the clay sensitivity (remoulded strength on the wall).' },
    { key: 'gammaSub', label: 'Submerged unit weight γ′', unit: 'kN/m³', type: 'number', default: 10, min: 4, max: 13, group: 'Soil', role: 'soil', bind: 'site.gammaSub' },
    { key: 'Vp', label: 'Submerged weight on the bucket at touchdown', unit: 'MN', type: 'number', default: 4, min: 0.05, max: 100, group: 'Installation', role: 'bc', help: 'Share of the submerged structure weight (bucket, jacket leg, ballast) carried by this bucket when the crane slacks.' },
    { key: 'rate', label: 'Penetration rate', unit: 'mm/s', type: 'number', default: 0.5, min: 0.05, max: 5, group: 'Installation', role: 'operation', help: 'Controlled by pump capacity; 1–2 m/h is common.' },
    { key: 'fsSuction', label: 'Safety factor on the suction limits', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Acceptance criteria', role: 'criteria', help: `Applied to plug heave, piping, cavitation and buckling limits. ${CONFIRM}` },
    { key: 'tens', label: 'Design uplift load in place', unit: 'MN', type: 'number', default: 5, min: 0, max: 200, group: 'In-place loads', role: 'bc', bind: 'jacket.pileAxialMax', help: 'Factored tension on the bucket under a jacket leg in the design storm (zero for a mono-bucket in compression).' },
    { key: 'Hin', label: 'Design horizontal load in place', unit: 'MN', type: 'number', default: 1.5, min: 0, max: 100, group: 'In-place loads', role: 'bc', bind: 'jacket.pileShearMax' },
    { key: 'eIn', label: 'Height of that load above the mudline', unit: 'm', type: 'number', default: 2, min: 0, max: 150, group: 'In-place loads', role: 'bc', help: 'Moment-to-shear ratio at the bucket lid; large for a mono-bucket (tens of metres), small under a jacket leg.' },
    { key: 'gammaR', label: 'Resistance factor on in-place capacity', unit: '-', type: 'number', default: 1.25, min: 1, max: 2, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'tOther', label: 'Other vessel time per foundation', unit: 'h', type: 'number', default: 18, min: 0, max: 200, group: 'Installation', role: 'operation', help: 'Positioning, lowering, levelling, pump recovery and grouting under the lid.' },
    { key: 'nStep', label: 'Penetration steps', type: 'int', default: 60, min: 10, max: 800, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nStep', label: 'Penetration steps', levels: [30, 60, 120, 240], qoi: 'pumpVol' },
  compute(inp, ctx) {
    const steel = steelOf(inp, ctx);
    const Do = inp.Do, t = inp.tw / 1e3, Di = Do - 2 * t, L = inp.L, gam = inp.gammaSub * 1e3, Vp = inp.Vp * 1e6;
    need(t < Do / 20, 'The skirt wall is too thick for a thin-walled bucket (t must be below D/20).');
    const clay = inp.soil === 'clay';
    const Ai = (Math.PI * Di * Di) / 4;
    const su = { su0: inp.su0 * 1e3, k: inp.suGrad * 1e3 };
    const phi = inp.phiDeg * DEG;
    const sandP = { Do, t, gamma: gam, KtanDeltaO: inp.Ktan, KtanDeltaI: inp.Ktan, Nq: bearingNq(phi), Ngamma: bearingNgamma(phi), Vprime: Vp, kf: inp.kf };
    const clayP = { Do, t, su, alphaO: inp.alpha, alphaI: inp.alpha, gamma: gam, Vprime: Vp };
    const at = (h) => {
      if (clay) {
        const r = suctionInstallClay({ ...clayP, h });
        // Equilibrium V′ + s·πDi²/4 = R (Houlsby & Byrne 2005a Eq. 7, the suction acts on the inside lid area).
        return { R: r.R, s: Math.max((r.R - Vp) / Ai, 0), lim: r.sAllow, a: 0 };
      }
      const r = suctionInstallSand({ ...sandP, h });
      return { R: r.R0, s: Math.max(r.s, 0), lim: r.sCrit, a: r.a };
    };
    const swp = Math.min(run(() => (clay ? selfWeightPenetrationClay({ ...clayP, L }) : selfWeightPenetrationSand({ ...sandP, L })), 'Self-weight penetration'), L);
    // Other suction limits: cavitation at the lid and elastic buckling of the skirt under external pressure.
    const sCav = 1025 * G * inp.depth + 101325 - 2300;
    const free = Math.max(L - swp, 0.1 * L);
    const sBuck = Math.max(c202ElasticBuckling({ r: (Do - t) / 2, t, l: Math.max(free, 0.5), E: steel.E, nu: steel.nu }).fEhLateral * t / ((Do - t) / 2), ringBucklingPressure({ E: steel.E, nu: steel.nu, t, D: Do }));
    const n = inp.nStep, v = inp.rate / 1e3;
    const hs = [], ss = [], lims = [], Rs = [];
    let vol = 0, sMaxNeed = 0, ucSoil = 0, hCrit = 0, tSuction = 0;
    for (let i = 0; i <= n; i++) {
      const h = Math.max((L * i) / n, 1e-3);
      const q = at(h);
      hs.push(h); ss.push(q.s); lims.push(q.lim); Rs.push(q.R);
      if (q.s > sMaxNeed) sMaxNeed = q.s;
      if (q.s > 0 && q.lim > 0 && q.s / q.lim > ucSoil) { ucSoil = q.s / q.lim; hCrit = h; }
      if (i > 0 && h > swp) {
        const dh = Math.min(h - swp, L / n), dt = dh / v;
        const pf = pumpFlow({ Di, v, k: clay ? 0 : inp.kPerm, s: q.s, h, a: q.a });
        vol += pf.Q * dt; tSuction += dt;
      }
    }
    const fin = at(L);
    const pfEnd = pumpFlow({ Di, v, k: clay ? 0 : inp.kPerm, s: fin.s, h: L, a: fin.a });
    // In-place capacity.
    const Wb = Vp;
    let tensile = 0, tensileRapid = 0;
    if (clay) { const c = caissonVerticalCapacityClay({ Do, t, L, su, W: Wb, gamma: gam }); tensile = c.vented; tensileRapid = c.sealed; } else { const c = tensileCapacitySand({ Do, t, h: L, gamma: gam, KtanDeltaO: inp.Ktan, KtanDeltaI: inp.Ktan, W: Wb, waterDepth: inp.depth }); tensile = c.drained; tensileRapid = c.rapid ?? c.drained; }
    const mom = run(() => bucketMomentCapacity(clay ? { soil: 'clay', D: Do, L, e: inp.eIn, su0: su.su0, k: su.k } : { soil: 'sand', D: Do, L, e: inp.eIn, phi, gamma: gam }), 'Moment capacity');
    const fsS = inp.fsSuction;
    const checks = [
      check(clay ? 'Required suction against plug heave (reverse end bearing)' : 'Required suction against piping (critical suction)', fin.s / 1e3 * (ucSoil / Math.max(fin.s / Math.max(fin.lim, 1e-9), 1e-12)), fin.lim / fsS / 1e3, 'kPa', clay ? 'Houlsby & Byrne (2005a) §4; DNV-RP-E303 (2021) §5.3' : 'Houlsby & Byrne (2005b) §5; OWA suction caisson guidelines (2019) §6.4',
        clay ? 'Guards against the soil plug being sucked up inside the bucket instead of the skirt going down. Plug heave stops penetration short of target and leaves a void under the lid.' : 'Guards against the upward seepage loosening the sand inside until it boils. Piping channels break the seal, the pump draws water instead of advancing the bucket, and the soil inside is left loose.'),
      check('Required suction against cavitation at the lid', sMaxNeed / 1e3, sCav / fsS / 1e3, 'kPa', 'DNV-ST-0126 (2021) §8 (installation); absolute pressure at the pump inlet', 'The pump cannot lower the pressure under the lid below the vapour pressure. In shallow water this limits the suction available to finish the penetration.'),
      check('Required suction against skirt buckling', sMaxNeed / 1e3, sBuck / fsS / 1e3, 'kPa', 'DNV-RP-C202 (2019) §3.4 (lateral pressure); Bresse–Bryan ring buckling as lower bound', 'Guards the unsupported part of the thin skirt against collapsing inwards under the pressure difference. Buckling is most critical at the start of suction, when most of the skirt is still above the seabed.'),
      check('In-place uplift: design tension against drained tensile capacity', inp.tens, tensile / inp.gammaR / 1e6, 'MN', clay ? 'DNV-RP-E303 (2021) §4 (vented capacity as sustained-load bound)' : 'Houlsby, Kelly & Byrne (2005); OWA guidelines (2019) §5.5 (no sustained tension in sand)', 'Guards against the bucket being pulled out under sustained or repeated tension. Sand offers only wall friction once pore suction has drained; cyclic tension ratchets the bucket upwards.'),
      check('In-place lateral load against moment capacity', inp.Hin, mom.H / inp.gammaR / 1e6, 'MN', 'Rigid-rotation method (Broms 1964; Murff & Hamilton 1993); DNV-ST-0126 (2021) §7.4', 'Guards against the bucket rotating in the soil under the horizontal load applied at its lever arm. Screening value without base resistance.'),
    ];
    checks[0].uc = Math.min(ucSoil * fsS, 99); checks[0].demand = round(ucSoil * fin.lim / 1e3, 5);
    const gov = worst(checks), ok = gov.uc <= 1 && swp > 0.3;
    const tInstall = tSuction / 3600, vdays = (tInstall + inp.tOther) / 24;
    const mSteel = (Math.PI * (Do - t) * t * L + (Math.PI * Do * Do) / 4 * 0.04 * 1.5) * steel.rho / 1e3; // skirt + stiffened lid (40 mm plate × 1.5)
    // Skirt length that keeps the suction within the soil limit: scan shorter skirts.
    let Lok = L;
    if (ucSoil * fsS > 1) { for (let q = n; q >= 1; q--) { const a2 = at(hs[q]); if (a2.s <= a2.lim / fsS) { Lok = hs[q]; break; } Lok = hs[q]; } }
    return {
      kpis: [
        { key: 'swp', label: 'Self-weight penetration', unit: 'm', value: swp, limit: 0.3, sense: 'min', explain: 'Depth the skirt reaches under the submerged weight alone. It must be enough to seal the rim (typically 0.5–1 m) before the pumps can build suction.' },
        { key: 'sFinal', label: 'Required suction at full penetration', unit: 'kPa', value: fin.s / 1e3, explain: 'Pressure difference across the lid needed to bring the skirt tip to target depth. It is the pump duty point and the load case for the lid and skirt.' },
        { key: 'sLimit', label: clay ? 'Suction limit for plug heave at full penetration' : 'Critical suction for piping at full penetration', unit: 'kPa', value: fin.lim / 1e3, explain: clay ? 'Suction at which the clay plug would fail in reverse end bearing and heave into the bucket.' : 'Suction at which the upward hydraulic gradient inside the bucket reaches the critical value and the sand plug loosens.' },
        { key: 'ucSuction', label: 'Suction utilisation against the soil limit', unit: '-', value: ucSoil * fsS, limit: 1, sense: 'max', explain: `Largest ratio of required suction to allowable suction over the penetration (at ${fmt(hCrit)} m), including the safety factor.` },
        { key: 'sCav', label: 'Cavitation limit', unit: 'kPa', value: sCav / 1e3, explain: 'Largest pressure difference physically available in this water depth (hydrostatic plus atmospheric minus vapour pressure).' },
        { key: 'sBuckle', label: 'Skirt buckling pressure', unit: 'kPa', value: sBuck / 1e3, explain: 'Elastic buckling pressure of the unsupported skirt length after self-weight penetration (DNV-RP-C202 lateral pressure, not below the long-cylinder value).' },
        { key: 'Rfinal', label: 'Penetration resistance at full depth', unit: 'MN', value: fin.R / 1e6, explain: 'Wall friction plus tip resistance without suction. The suction replaces the weight that is missing to overcome it and, in sand, also reduces it through seepage.' },
        { key: 'pumpQ', label: 'Pump flow at full penetration', unit: 'm³/h', value: pfEnd.Q * 3600, explain: 'Water to remove from under the lid: displaced volume at the penetration rate plus seepage through the sand plug. Sizes the pump skid.' },
        { key: 'pumpVol', label: 'Total pumped volume', unit: 'm³', value: vol, explain: 'Water pumped from self-weight penetration to target depth; larger than the bucket volume in permeable sand because of seepage.' },
        { key: 'tInstall', label: 'Suction penetration time', unit: 'h', value: tInstall, explain: 'Pumping time at the stated penetration rate from self-weight penetration to target depth.' },
        { key: 'tensile', label: 'Drained tensile capacity in place', unit: 'MN', value: tensile / 1e6, explain: 'Sustained uplift capacity: wall friction plus weight (vented case in clay). Use it for mean and long-duration tension.' },
        { key: 'tensileRapid', label: 'Rapid-loading tensile capacity', unit: 'MN', value: tensileRapid / 1e6, explain: 'Capacity under short wave loading when passive suction develops under the lid (sealed case in clay, cavitation-limited in sand). Not to be relied on for repeated tension without a specific assessment.' },
        { key: 'Hult', label: 'Ultimate horizontal load at the stated lever arm', unit: 'MN', value: mom.H / 1e6, explain: 'Horizontal load at the given height above the mudline that rotates the bucket in the soil (rigid-rotation screening).' },
        { key: 'vesselDays', label: 'Installation vessel time per foundation', unit: 'days', value: vdays, explain: 'Suction penetration plus the other installation steps. No hammer is used, so there is no piling noise and the operation is reversible.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Required suction against penetration', x: { label: 'Suction', unit: 'kPa' }, y: { label: 'Tip penetration', unit: 'm' },
          series: [{ name: 'Required suction', x: ss.map((q) => q / 1e3), y: hs.map((h) => -h) }, { name: clay ? 'Plug-heave limit' : 'Piping limit', x: lims.map((q) => Math.min(q, 5 * Math.max(sMaxNeed, 1e3)) / 1e3), y: hs.map((h) => -h), dash: true }, { name: 'Cavitation limit', x: [Math.min(sCav, 5 * Math.max(sMaxNeed, 1e3)) / 1e3, Math.min(sCav, 5 * Math.max(sMaxNeed, 1e3)) / 1e3], y: [0, -L], dash: true }],
          explain: 'Suction needed at each depth beside the soil and cavitation limits. Suction starts where self-weight penetration ends; the required curve must stay left of both limits to target depth.' },
        { type: 'line', title: 'Penetration resistance against depth', x: { label: 'Resistance without suction', unit: 'MN' }, y: { label: 'Tip penetration', unit: 'm' },
          series: [{ name: 'Resistance', x: Rs.map((q) => q / 1e6), y: hs.map((h) => -h) }, { name: 'Submerged weight', x: [Vp / 1e6, Vp / 1e6], y: [0, -L], dash: true }],
          explain: 'Wall friction and tip resistance as the skirt goes down, with the available weight. Their crossing is the self-weight penetration; the gap below it must be closed by suction.' },
      ],
      tables: [
        { title: 'Installation by depth', columns: [{ key: 'h', label: 'Penetration', unit: 'm' }, { key: 'R', label: 'Resistance', unit: 'MN' }, { key: 's', label: 'Required suction', unit: 'kPa' }, { key: 'l', label: 'Soil limit', unit: 'kPa' }, { key: 'u', label: 's / limit' }],
          rows: hs.filter((_, i) => i % Math.max(1, Math.round(n / 15)) === 0 || i === n).map((h) => { const q = at(h); return { h: round(h), R: round(q.R / 1e6), s: round(q.s / 1e3), l: round(q.lim / 1e3), u: q.lim > 0 ? round(q.s / q.lim) : 0 }; }),
          explain: 'Predicted installation log: resistance, required suction and the soil limit at each depth, for comparison with the pump pressure record during installation.' },
      ],
      notes: [`${clay ? `Clay, su = ${inp.su0} + ${inp.suGrad}·z kPa, adhesion ${inp.alpha}` : `Sand, φ′ = ${inp.phiDeg}°, K·tanδ = ${inp.Ktan}, permeability ratio ${inp.kf}`}; L/D = ${fmt(L / Do)}.`,
        `Self-weight penetration ${fmt(swp)} m; suction of ${fmt(fin.s / 1e3)} kPa completes the penetration (${fmt(ucSoil * fsS * 100)} % of the allowable soil limit).`, `${gov.name} governs (utilisation ${fmt(gov.uc)}).`,
        clay ? 'In clay the required suction is taken on the inside lid area (πDi²/4).' : 'In sand seepage reduces the tip and inside resistance; without it the required suction would be far higher.'],
      advice: [
        ok ? adv('engineering', sev(gov.uc), `The bucket can be installed: self-weight penetration ${fmt(swp)} m, required suction ${fmt(fin.s / 1e3)} kPa against ${fmt(fin.lim / fsS / 1e3)} kPa allowable for the soil and ${fmt(Math.min(sCav, sBuck) / fsS / 1e3)} kPa for cavitation and buckling.`,
          `Specify pumps for ${fmt(pfEnd.Q * 3600 * 1.3)} m³/h at ${fmt(fin.s / 1e3 * 1.2)} kPa and monitor pressure against penetration; stop and cycle the suction if the pressure rises above the predicted curve by more than the safety margin.`,
          'DNV-ST-0126 (2021) §8; OWA suction caisson design guidelines (2019); Houlsby & Byrne (2005a, b)',
          `Suction installation needs no impact hammer: no piling noise for marine mammals and full removal at end of life by reversing the pumps. The bucket holds about ${fmt(mSteel)} t of steel — ${steelCo2Text(mSteel)}`)
          : adv('engineering', 'fail', swp <= 0.3 ? `Self-weight penetration is only ${fmt(swp)} m: the rim may not seal.` : `${gov.name}: ${fmt(gov.demand)} against ${fmt(gov.capacity)} ${gov.unit}.`,
            swp <= 0.3 ? 'Add ballast or reduce the wall and stiffener area in contact with the soil so that the skirt penetrates at least 0.5 m before pumping.' : gov.name.includes('piping') || gov.name.includes('plug') ? `Shorten the skirt to about ${fmt(Lok)} m or widen the bucket (lower L/D needs less suction), add ballast (now ${inp.Vp} MN), or use water injection at the tip to cut resistance.` : gov.name.includes('buckling') ? `Thicken the skirt (buckling pressure ∝ t^2.5) or add ring stiffeners; ${fmt(inp.tw * (gov.uc) ** 0.4)} mm would pass.` : gov.name.includes('cavitation') ? 'The water is too shallow for the suction needed: reduce resistance (larger diameter, shorter skirt, ballast).' : 'Increase diameter or skirt length for in-place capacity, or reduce the tension by widening the jacket footprint.',
            'OWA suction caisson design guidelines (2019) §6; Houlsby & Byrne (2005a, b); DNV-RP-C202 (2019) §3.4',
            `A wider, shorter bucket of the same capacity uses similar steel (now about ${fmt(mSteel)} t — ${steelCo2Text(mSteel)}) but installs with less suction; ballast is the lowest-carbon way to add driving force.`),
        adv('economic', 'info', `Installation takes about ${fmt(tInstall)} h of pumping and ${fmt(vdays)} vessel days per foundation; the bucket is about ${fmt(mSteel)} t of steel.`,
          'Compare with a driven pile on total installed cost: the bucket uses more fabricated steel per unit capacity but saves the hammer spread, noise mitigation and, at decommissioning, complete removal. Installation risk concentrates in soil variability, so budget CPTs at every bucket position.',
          'AACE RP 18R-97 Class 4 quantities; installation risk allowance per AACE RP 57R-09; decommissioning provisions in whole-life costing (ISO 15686-5)'),
      ],
      quantities: { steel_t: round(mSteel), vessel_days: { foundation_installation: round(vdays) } },
    };
  },
  verify: [
    { name: 'Required suction in uniform clay: s = (αo·su·πDo·h + αi·su·πDi·h + (γ′h + Nc·su)·πD·t − V′) / (πDi²/4)', source: 'Houlsby & Byrne (2005a) Eq. (7), hand calculation: Do = 10 m, t = 30 mm, h = 8 m, su = 20 kPa, α = 0.5, γ′ = 6 kN/m³, V′ = 1 MN → outside 2513.3 kN, inside 2498.2 kN, tip (48 + 180) kPa × π × 9.97 × 0.03 = 214.2 kN, R = 5225.7 kN; s = 4225.7 kN / 77.60 m² = 54.45 kPa',
      inputs: { soil: 'clay', Do: 10, L: 8, tw: 30, su0: 20, suGrad: 0, alpha: 0.5, gammaSub: 6, Vp: 1 }, expect: [{ kpi: 'sFinal', value: 54.45, tol: 0.002 }, { kpi: 'Rfinal', value: 5.2257, tol: 0.001 }] },
    { name: 'Cavitation limit: s = ρ·g·hw + p_atm − p_vap', source: 'Hydrostatics: 1025 × 9.80665 × 35 + 101325 − 2300 = 450.8 kPa', inputs: { depth: 35 }, expect: [{ kpi: 'sCav', value: 450.85, tol: 0.001 }] },
  ],
};
