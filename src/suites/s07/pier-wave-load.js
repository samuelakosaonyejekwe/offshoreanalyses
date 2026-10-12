// Wave and current loading on a bridge pier or offshore tower shaft: Morison force with fifth-order Stokes
// kinematics integrated to the instantaneous surface, diffraction check, breaking-wave slam and deck air gap.
import { stokes5Wave, stokes5Kinematics, cnoidalWave, cnoidalKinematics, wavelength, breakingLimits } from '../../engine/metocean/waves.js';
import { pileForceClosedForm, maccamyFuchs, memberForceHistory, slamForcePerLength, marineGrowth, keuleganCarpenter } from '../../engine/hydro/morison.js';
import { waveInDeckForce } from '../../engine/hydro/impact.js';
import { G, CO2_STEEL, CO2_CONC, need, run, finite, r, f, linspace, check, eng, eco, slrInputs, slrRows, SLR_BASIS } from './_common.js';

/**
 * Morison base shear and mudline moment history over one wave period (midpoint rule over the wetted height).
 * `kin(z, t)` returns { eta, u, a } — surface elevation, horizontal velocity including current and horizontal acceleration.
 */
function waveHistory({ T, d, kin, De, Cd, Cm, rho, nz, nt = 72 }) {
  const out = { t: [], F: [], M: [], Fmax: 0, Mmax: 0, tF: 0, prof: null };
  const A = Math.PI * De * De / 4;
  for (let it = 0; it < nt; it++) {
    const t = it * T / nt;
    const eta = kin(0, t).eta;
    const dz = (d + eta) / nz;
    let F = 0, M = 0; const prof = [];
    for (let i = 0; i < nz; i++) {
      const z = -d + (i + 0.5) * dz;
      const k = kin(Math.min(z, eta), t);
      const fz = 0.5 * rho * Cd * De * Math.abs(k.u) * k.u + rho * Cm * A * k.a;
      F += fz * dz; M += fz * dz * (z + d); prof.push([z, fz]);
    }
    out.t.push(t); out.F.push(F); out.M.push(M);
    if (Math.abs(F) > out.Fmax) { out.Fmax = Math.abs(F); out.tF = t; out.prof = prof; }
    out.Mmax = Math.max(out.Mmax, Math.abs(M));
  }
  return out;
}

export default {
  id: 's07-pier-wave-load', suite: 7, title: 'Wave and current load on a pier or tower',
  summary: 'Base shear and overturning moment on a bridge pier or tower shaft from the design wave and current, the breaking-wave slam and the air gap under the deck.',
  covers: [2], fidelity: 'reduced-order',
  standards: ['DNV-RP-C205 (2019) §3, §6, §7, §8', 'ISO 19901-1:2015 (metocean actions)', 'IEC 61400-3-1:2019 Annex D (breaking-wave slam, Wienke & Oumeraci)', 'API RP 2A-WSD (21st ed.) §C17.6 (wave-in-deck)', 'AASHTO Guide Specifications for Bridges Vulnerable to Coastal Storms (2008)'],
  models: ['Morison equation', 'Stokes fifth-order wave theory', 'Cnoidal wave theory', 'Airy wave theory with Wheeler stretching', 'MacCamy–Fuchs diffraction', 'Wave slamming (Wienke–Oumeraci form)', 'Wave-in-deck silhouette force', 'Wave breaking limits (Miche, McCowan)'],
  limits: 'A single vertical shaft of constant effective diameter in a long-crested regular design wave on a uniform current; pile groups, shielding, caps near the surface, 3-D diffraction around wide piers beyond MacCamy–Fuchs and vortex-induced vibration are not modelled. The slam load is an impulsive screening value: dynamic amplification must be taken from a structural dynamic analysis; CFD and model-test precision are not implied.',
  publishes: ['marine.pierShear', 'marine.pierMoment', 'marine.pierSlam', 'marine.crestElevation'],
  searchTerms: 'Morison wave force bridge pier breaking wave slam air gap Stokes fifth order',
  inputs: [
    { key: 'D', label: 'Shaft diameter / width', unit: 'm', type: 'number', default: 3, min: 0.2, max: 30, group: 'Structure', role: 'geometry', design: true },
    { key: 'mg', label: 'Marine growth thickness', unit: 'mm', type: 'number', default: 80, min: 0, max: 300, group: 'Structure', role: 'ic', help: 'Initial fouled condition assumed over the wetted height.' },
    { key: 'Cd', label: 'Drag coefficient', unit: '-', type: 'number', default: 1.05, min: 0, max: 2.5, group: 'Structure', role: 'criteria', uncertain: { dist: 'normal', cov: 0.15 }, help: '0.65 smooth, 1.05 rough circular cylinder (DNV-RP-C205 §6.7); about 2 for rectangular piers.' },
    { key: 'Cm', label: 'Inertia coefficient', unit: '-', type: 'number', default: 2.0, min: 1, max: 3, group: 'Structure', role: 'criteria', uncertain: { dist: 'normal', cov: 0.1 } },
    { key: 'deckZ', label: 'Deck soffit level above still water', unit: 'm', type: 'number', default: 10, min: 0.5, max: 80, group: 'Structure', role: 'geometry', design: true, help: 'Underside of the lowest deck beam or platform, above the design still-water level before sea-level rise.' },
    { key: 'deckW', label: 'Deck width exposed to the wave', unit: 'm', type: 'number', default: 12, min: 0.5, max: 200, group: 'Structure', role: 'geometry' },
    { key: 'deckH', label: 'Deck silhouette height', unit: 'm', type: 'number', default: 2.5, min: 0.2, max: 20, group: 'Structure', role: 'geometry' },
    { key: 'h', label: 'Water depth at design still-water level', unit: 'm', type: 'number', default: 18, min: 2, max: 150, group: 'Wave and current boundary', role: 'bc', bind: 'site.depth' },
    { key: 'H', label: 'Design (maximum) wave height', unit: 'm', type: 'number', default: 9, min: 0.1, max: 35, group: 'Wave and current boundary', role: 'bc', bind: 'site.Hmax100', uncertain: { dist: 'lognormal', cov: 0.1 } },
    { key: 'T', label: 'Associated wave period', unit: 's', type: 'number', default: 11, min: 3, max: 25, group: 'Wave and current boundary', role: 'bc', uncertain: { dist: 'normal', cov: 0.08 } },
    { key: 'U', label: 'In-line current', unit: 'm/s', type: 'number', default: 1.0, min: 0, max: 4, group: 'Wave and current boundary', role: 'bc', bind: 'site.Uc100' },
    ...slrInputs({ group: 'Wave and current boundary' }),
    { key: 'rhoW', label: 'Sea-water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Wave and current boundary', role: 'environment', bind: 'site.rhoW' },
    { key: 'slam', label: 'Breaking-wave slam', type: 'select', default: 'auto', options: [{ value: 'auto', label: 'Apply when the wave is near breaking' }, { value: 'yes', label: 'Always apply' }, { value: 'no', label: 'Do not apply' }], group: 'Impact', role: 'criteria' },
    { key: 'Cs', label: 'Slam coefficient', unit: '-', type: 'number', default: 6.283, min: 3, max: 6.3, group: 'Impact', role: 'criteria', help: '2π (Wienke & Oumeraci / IEC 61400-3-1) or π (von Kármán); 5.15 Campbell–Weynberg.' },
    { key: 'lambda', label: 'Curling factor λ', unit: '-', type: 'number', default: 0.4, min: 0.1, max: 0.6, group: 'Impact', role: 'criteria', uncertain: { dist: 'normal', cov: 0.2 }, help: 'Share of the crest height that hits as a vertical wall of water: 0.4–0.5 for plunging breakers.' },
    { key: 'CdDeck', label: 'Deck silhouette drag coefficient', unit: '-', type: 'number', default: 2.0, min: 1, max: 3, group: 'Impact', role: 'criteria' },
    { key: 'gapMin', label: 'Required air gap', unit: 'm', type: 'number', default: 1.5, min: 0, max: 5, group: 'Criteria', role: 'criteria', help: 'ISO 19902 / API practice asks for 1.5 m above the design crest.' },
    { key: 'Rshear', label: 'Design lateral resistance at the mudline', unit: 'kN', type: 'number', default: 4000, min: 10, max: 1e6, group: 'Criteria', role: 'criteria', help: 'From the foundation and section design; compared with the factored wave load.' },
    { key: 'gammaF', label: 'Load factor on environmental load', unit: '-', type: 'number', default: 1.35, min: 1, max: 2, group: 'Criteria', role: 'criteria' },
    { key: 'nz', label: 'Integration points over the wetted height', type: 'int', default: 40, min: 5, max: 640, group: 'Solver', role: 'control' },
  ],
  mesh: { key: 'nz', label: 'Integration points', levels: [10, 20, 40, 80], qoi: 'Fmax' },
  compute(inp) {
    const nz = Math.round(inp.nz);
    const De = marineGrowth({ D: inp.D, thickness: inp.mg / 1000 }).De;
    const state = (rise, nzs = nz, nt = 72) => {
      const d = inp.h + rise;
      const lim = breakingLimits({ T: inp.T, d });
      const Hb = Math.min(...Object.values(lim).filter((v) => Number.isFinite(v) && v > 0));
      const H = Math.min(inp.H, 0.98 * Hb); // a wave above the breaking limit has already broken: cap and treat as breaking
      let hist, theory = 'Stokes 5th order', crest, uCrest, c;
      let sol = null; try { sol = stokes5Wave({ H, d, T: inp.T, currentEulerian: inp.U }); } catch { sol = null; }
      if (sol && sol.valid) {
        hist = waveHistory({ T: sol.T, d, kin: (z, t) => { const k = stokes5Kinematics(sol, { x: 0, z, t, rho: inp.rhoW }); return { eta: k.eta, u: k.u, a: k.axTotal }; }, De, Cd: inp.Cd, Cm: inp.Cm, rho: inp.rhoW, nz: nzs, nt });
        crest = sol.crest; c = sol.c; uCrest = stokes5Kinematics(sol, { x: 0, z: sol.crest, t: 0 }).u;
      } else {
        // long waves in shallow water (Ursell number above the Stokes range): cnoidal theory; Airy with Wheeler stretching
        // only where no cnoidal solution exists for the period
        let cn = null; try { cn = cnoidalWave({ H, d, T: inp.T }); } catch { cn = null; }
        if (cn && Number.isFinite(cn.c) && Number.isFinite(cn.crest)) {
          theory = 'Cnoidal';
          hist = waveHistory({ T: cn.T, d, kin: (z, t) => { const k = cnoidalKinematics(cn, { x: 0, z, t, rho: inp.rhoW }); return { eta: k.eta, u: k.u + inp.U, a: k.ax }; }, De, Cd: inp.Cd, Cm: inp.Cm, rho: inp.rhoW, nz: nzs, nt });
          crest = cn.crest; c = cn.c; uCrest = cnoidalKinematics(cn, { x: 0, z: cn.crest, t: 0 }).u + inp.U;
        } else {
          theory = 'Airy with Wheeler stretching (outside Stokes and cnoidal ranges)';
          const m = memberForceHistory({ p1: [0, 0, -d], p2: [0, 0, 0.75 * H], D: De, Cd: inp.Cd, Cm: inp.Cm, rho: inp.rhoW, wave: { H, T: inp.T, depth: d }, current: inp.U > 0 ? { speed: inp.U } : undefined, stretching: 'wheeler', nSeg: Math.max(1, Math.round(nzs / 8)), about: [0, 0, -d] });
          hist = { t: Array.from(m.t), F: Array.from(m.Fh), M: Array.from(m.Mh), Fmax: m.FhMax, Mmax: m.MhMax, prof: null };
          const L = wavelength(inp.T, d); c = L / inp.T; crest = 0.5 * H * (1 + Math.PI * H / L); // second-order crest estimate
          uCrest = Math.PI * H / inp.T / Math.tanh(2 * Math.PI * d / L) + inp.U;
        }
      }
      return { d, H, Hb, hist, theory, crest, uCrest, c, L: wavelength(inp.T, d) };
    };
    const s = run(() => state(inp.slr), 'Wave load');
    const { hist } = s;
    finite(hist.Fmax, 'Wave force');
    const lin = pileForceClosedForm({ H: s.H, T: inp.T, depth: s.d, D: De, Cd: inp.Cd, Cm: inp.Cm, rho: inp.rhoW });
    const mf = maccamyFuchs({ H: s.H, T: inp.T, depth: s.d, D: De, rho: inp.rhoW });
    const DL = De / s.L, diffraction = DL > 0.2;
    const KC = keuleganCarpenter(Math.PI * s.H / inp.T / Math.tanh(2 * Math.PI * s.d / s.L), inp.T, De);
    const nearBreaking = s.H >= 0.9 * s.Hb;
    const slamOn = inp.slam === 'yes' || (inp.slam === 'auto' && nearBreaking);
    // slam: line force ½ ρ Cs D C² over the curling height λ ηb, acting near the crest (Wienke & Oumeraci 2005)
    const Fslam = slamOn ? slamForcePerLength({ D: De, v: s.c, Cs: inp.Cs, rho: inp.rhoW }) * inp.lambda * s.crest : 0;
    const zSlam = s.d + s.crest * (1 - inp.lambda / 2), tSlam = 13 * De / (32 * s.c);
    const gap = inp.deckZ - inp.slr - s.crest;
    const deck = gap < 0 ? waveInDeckForce({ crestElevation: s.crest, deckBottom: inp.deckZ - inp.slr, deckTop: inp.deckZ - inp.slr + inp.deckH, width: inp.deckW, u: s.uCrest, Cd: inp.CdDeck, rho: inp.rhoW }).force : 0;
    const Fdes = inp.gammaF * (hist.Fmax + deck), R = inp.Rshear * 1e3;
    const ok = Fdes <= R, gapOk = gap >= inp.gapMin;
    const scen = slrRows(inp).map((sc) => { const q = state(sc.rise, Math.max(10, Math.round(nz / 4)), 36); return { name: sc.name, rise: sc.rise, d: r(q.d, 3), F: r(q.hist.Fmax / 1e3), M: r(q.hist.Mmax / 1e3), crest: r(q.crest, 3), gap: r(inp.deckZ - sc.rise - q.crest, 3) }; });
    // required diameter change is not the lever here: report the load per metre of diameter for the designer
    const dFdD = hist.Fmax / De;
    const raise = Math.max(0, inp.gapMin - gap);

    const advice = [];
    advice.push(ok
      ? eng('ok', `Maximum wave-plus-current base shear is ${f(hist.Fmax / 1e3)} kN (overturning ${f(hist.Mmax / 1e6)} MN·m, ${s.theory}); factored ${f(Fdes / 1e3)} kN against ${f(inp.Rshear)} kN resistance.`,
        `Utilisation ${f(Fdes / R, 2)}. ${diffraction ? 'The shaft is in the diffraction regime — use the MacCamy–Fuchs value for the inertia part.' : `Drag-to-inertia ratio ${f(lin.FD / Math.max(lin.FI, 1), 2)}: ${lin.FD > lin.FI ? 'drag' : 'inertia'} governs`}; confirm Cd and Cm for the actual shape and roughness.`,
        'DNV-RP-C205 (2019) §6.2–6.7; ISO 19901-1:2015', `Marine growth of ${f(inp.mg)} mm is carrying ${f((1 - inp.D / De) * 100, 2)} % of the load width; a cleaning regime lets a slimmer shaft and foundation do the same job over the life.`)
      : eng('fail', `Factored base shear ${f(Fdes / 1e3)} kN exceeds the ${f(inp.Rshear)} kN lateral resistance (utilisation ${f(Fdes / R, 2)}); unfactored ${f(hist.Fmax / 1e3)} kN, overturning ${f(hist.Mmax / 1e6)} MN·m.`,
        `Increase the lateral resistance to at least ${f(Fdes / 1e3)} kN, or cut the load: a slimmer or streamlined shaft sheds about ${f(dFdD / 1e3)} kN per metre of diameter, and marine-growth control removes up to ${f((1 - inp.D / De) * 100, 2)} %.`,
        'DNV-RP-C205 (2019) §6.2–6.7; ISO 19901-1:2015', 'Reducing the hydrodynamic load at source saves foundation steel or concrete — roughly in proportion to the overturning moment — rather than adding material to resist it.'));
    if (slamOn) advice.push(eng('warn', `Breaking-wave slam: ${f(Fslam / 1e3)} kN acting ${f(zSlam)} m above the mudline for about ${f(tSlam * 1000, 2)} ms (H/Hb = ${f(s.H / s.Hb, 2)}), ${f(Fslam / Math.max(hist.Fmax, 1), 2)} × the Morison force.`,
      'Run a dynamic analysis with the slam time history; its duration is short against the first natural period, so the dynamic response, not the peak force, governs. Model tests are warranted for plunging breakers.', 'IEC 61400-3-1:2019 Annex D; Wienke & Oumeraci (2005); DNV-RP-C205 §8'));
    advice.push(eng(gapOk ? 'info' : 'fail', `Crest elevation ${f(s.crest)} m gives an air gap of ${f(gap, 3)} m under the deck after ${f(inp.slr)} m of sea-level rise (required ${f(inp.gapMin)} m)${deck > 0 ? `; the wave strikes the deck with ${f(deck / 1e3)} kN` : ''}.`,
      gapOk ? `The air gap is consumed at about ${f(inp.slr + gap - inp.gapMin, 2)} m of total sea-level rise — record it as the adaptation trigger.` : `Raise the deck soffit by ${f(raise, 2)} m, or design the deck and its bearings for wave-in-deck loads including uplift (not covered here).`,
      `ISO 19902:2020 §6.3 air gap; AASHTO coastal-storm guide (2008); ${SLR_BASIS}`, 'Setting the soffit for the end-of-life sea level now avoids jacking or replacing the superstructure later — the largest avoidable material and disruption item for a coastal bridge.'));
    advice.push(eco(ok ? 'info' : 'warn', `The overturning moment of ${f(hist.Mmax / 1e6)} MN·m sizes the foundation; it grows by ${f(((scen[scen.length - 1].M / Math.max(scen[0].M, 1e-9)) - 1) * 100, 2)} % between the first and last sea-level scenario.`,
      'Foundation quantity scales roughly with the moment: compare the cost of a slimmer shaft, growth control and a higher soffit against added foundation steel or concrete on a whole-life basis.', 'Whole-life costing per ISO 15686-5; AACE RP 18R-97 quantity-based estimating',
      `As a screening guide each 10 % of extra moment costs about 10 % more foundation material (${f(CO2_STEEL)} t CO₂e per t steel, ${f(CO2_CONC)} t per m³ concrete, screening factors).`));

    const prof = hist.prof || linspace(-s.d, 0, 20).map((z) => [z, 0]);
    return {
      kpis: [
        { key: 'Fmax', label: 'Maximum base shear (wave + current)', value: r(hist.Fmax / 1e3), unit: 'kN', explain: 'Largest in-line Morison force on the shaft during the design wave cycle, integrated from the bed to the instantaneous surface. Unfactored.' },
        { key: 'Mmax', label: 'Maximum mudline overturning moment', value: r(hist.Mmax / 1e3), unit: 'kN·m', explain: 'Largest moment of the distributed wave force about the bed. It sizes the foundation and the shaft section at the mudline.' },
        { key: 'Flin', label: 'Linear-theory base shear (closed form)', value: r(lin.Fmax / 1e3), unit: 'kN', explain: 'Morison force from Airy kinematics integrated to still water, no current. The gap to the nonlinear value shows the effect of crest kinematics and current.' },
        { key: 'FI', label: 'Inertia force amplitude (linear)', value: r(lin.FI / 1e3), unit: 'kN', explain: 'Inertia part of the linear Morison force. Dominant for large diameters and small waves; scales with D² and H.' },
        { key: 'Fmf', label: 'MacCamy–Fuchs diffraction force', value: r(mf.F / 1e3), unit: 'kN', explain: 'Linear diffraction force on the cylinder. Use instead of Morison inertia when D/L exceeds about 0.2, where the shaft scatters the wave.' },
        { key: 'DL', label: 'Diameter / wavelength', value: r(DL, 3), unit: '-', explain: 'Below about 0.2 the shaft is slender and Morison applies; above it diffraction reduces the effective inertia coefficient.' },
        { key: 'KC', label: 'Keulegan–Carpenter number', value: r(KC, 3), unit: '-', explain: 'Orbital excursion relative to diameter. Below about 5 inertia dominates; above about 25 drag dominates and Cd matters most.' },
        { key: 'crest', label: 'Design crest elevation above still water', value: r(s.crest), unit: 'm', explain: 'Height of the wave crest above the design still-water level from the nonlinear wave theory; higher than H/2 because crests are peaked.' },
        { key: 'gap', label: 'Air gap under the deck', value: r(gap, 3), unit: 'm', limit: inp.gapMin, sense: 'min', explain: 'Clearance between the design crest and the deck soffit after sea-level rise. Negative means the wave hits the deck.' },
        { key: 'Fslam', label: 'Breaking-wave slam force', value: r(Fslam / 1e3), unit: 'kN', explain: 'Impulsive force of a plunging breaker on the shaft near the crest, lasting tens of milliseconds. Zero when slam is not applied.' },
        { key: 'Fdeck', label: 'Wave-in-deck horizontal force', value: r(deck / 1e3), unit: 'kN', explain: 'Silhouette drag force on the deck when the crest exceeds the soffit. Zero with a positive air gap.' },
        { key: 'HHb', label: 'Wave height / breaking limit', value: r(s.H / s.Hb, 3), unit: '-', explain: 'Design wave relative to the largest wave the depth and period allow. Near 1 the wave breaks at the structure and slam must be considered.' },
      ],
      checks: [
        check('Factored base shear against lateral resistance', Fdes / 1e3, inp.Rshear, 'kN', 'ISO 19902:2020 / EN 1990 partial-factor format; DNV-RP-C205 §6', 'Guards against the shaft or its foundation yielding laterally under the extreme wave and current. Exceeded, the pier displaces or the foundation fails in lateral bearing.'),
        check('Air gap', inp.gapMin, Math.max(gap, 1e-3), 'm', 'ISO 19902:2020 §6.3; API RP 2A §C17.6', 'Guards against wave crests striking the deck, where horizontal and uplift loads rise by an order of magnitude and bridge decks are lifted off their bearings.'),
      ],
      plots: [
        { type: 'line', title: 'Base shear and overturning moment over one wave', x: { label: 'Time', unit: 's' }, y: { label: 'Shear (kN) / moment (10 kN·m)', unit: '' },
          series: [{ name: 'Base shear', x: hist.t.map((t) => r(t, 3)), y: hist.F.map((v) => r(v / 1e3, 3)) }, { name: 'Moment / 10', x: hist.t.map((t) => r(t, 3)), y: hist.M.map((v) => r(v / 1e4, 3)) }],
          explain: 'Force and moment through the wave cycle with the crest at t = 0. A peak at the crest indicates drag dominance; a peak a quarter-period earlier indicates inertia dominance.' },
        { type: 'line', title: 'Force per unit height at maximum base shear', x: { label: 'Force intensity', unit: 'kN/m' }, y: { label: 'Elevation to still water', unit: 'm' },
          series: [{ name: 'Morison line load', x: prof.map((p) => r(p[1] / 1e3, 3)), y: prof.map((p) => r(p[0], 3)) }],
          explain: 'Distribution of the hydrodynamic load over the shaft at the instant of maximum shear. The load concentrates near the surface, which is why the moment arm is long and the air gap matters.' },
        { type: 'line', title: 'Load versus sea-level rise', x: { label: 'Sea-level rise', unit: 'm' }, y: { label: 'Base shear (kN) / air gap (m × 100)', unit: '' },
          series: [{ name: 'Base shear', x: scen.map((q) => q.rise), y: scen.map((q) => q.F) }, { name: 'Air gap × 100', x: scen.map((q) => q.rise), y: scen.map((q) => r(q.gap * 100, 3)) }],
          explain: 'Wave load and remaining air gap as the sea rises. The load changes slowly; the air gap falls one-for-one and usually sets the adaptation date.' },
      ],
      tables: [
        { title: 'Sea-level-rise scenarios', columns: [{ key: 'name', label: 'Scenario' }, { key: 'rise', label: 'Rise', unit: 'm' }, { key: 'd', label: 'Depth', unit: 'm' }, { key: 'F', label: 'Base shear', unit: 'kN' }, { key: 'M', label: 'Moment', unit: 'kN·m' }, { key: 'crest', label: 'Crest', unit: 'm' }, { key: 'gap', label: 'Air gap', unit: 'm' }],
          rows: scen, explain: 'The design wave on the deeper water of each scenario. Use the row for the end of the design life for the foundation loads and read the trigger for deck raising from the air-gap column.' },
        { title: 'Load components', columns: [{ key: 'q', label: 'Component' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [['Drag amplitude (linear)', lin.FD / 1e3, 'kN'], ['Inertia amplitude (linear)', lin.FI / 1e3, 'kN'], ['MacCamy–Fuchs effective Cm', mf.CmEff, '-'], ['Effective diameter with growth', De, 'm'], ['Wavelength', s.L, 'm'], ['Wave celerity', s.c, 'm/s'], ['Crest particle velocity', s.uCrest, 'm/s'], ['Slam height λ ηb', slamOn ? inp.lambda * s.crest : 0, 'm'], ['Slam duration', slamOn ? tSlam * 1000 : 0, 'ms']].map(([q, v, u]) => ({ q, v: r(v), u })),
          explain: 'The pieces of the hydrodynamic load and the kinematic quantities behind them, for comparison with a hand check or a model-test report.' },
      ],
      notes: [`${s.theory}; design wave ${f(s.H)} m${s.H < inp.H ? ` (capped from ${f(inp.H)} m at the breaking limit)` : ''} in ${f(s.d)} m of water.`,
        `${lin.FD > lin.FI ? 'Drag' : 'Inertia'}-dominated loading (KC = ${f(KC)}); ${diffraction ? 'diffraction regime: prefer the MacCamy–Fuchs force.' : 'slender regime: Morison is applicable.'}`],
      advice,
      publish: { 'marine.pierShear': hist.Fmax, 'marine.pierMoment': hist.Mmax, 'marine.pierSlam': Fslam, 'marine.crestElevation': s.crest },
    };
  },
  verify: [
    { name: 'Small-amplitude inertia limit: F = ρ Cm (π D²/4) g (H/2) tanh(kd)', source: 'Morison et al. (1950); Dean & Dalrymple (1991) Eq. 8.33. Hand calculation: D = 2 m, H = 0.2 m, T = 8 s, d = 20 m, Cm = 2, Cd = 0: k = 0.07079 rad/m, tanh(kd) = 0.8886 → F = 1025 × 2 × π × 9.80665 × 0.1 × 0.8886 = 5.61 kN',
      inputs: { D: 2, mg: 0, Cd: 0, Cm: 2, H: 0.2, T: 8, h: 20, U: 0, slr: 0, slam: 'no', nz: 80 }, expect: [{ kpi: 'Fmax', value: 5.612, tol: 0.01 }, { kpi: 'FI', value: 5.612, tol: 0.003 }] },
  ],
};
