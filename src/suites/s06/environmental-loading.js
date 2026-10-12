// Wave, current, wind and vessel-induced loading on a dolphin: Morison loading on the piles integrated to the
// instantaneous surface, wind on the superstructure, and the share of moored-vessel wind / current force.
import { airyKinematics, stokes5Wave, stokes5Kinematics, selectWaveTheory } from '../../engine/metocean/waves.js';
import { currentPowerLaw } from '../../engine/metocean/current.js';
import { windForce } from '../../engine/metocean/wind.js';
import { pileForceClosedForm, keuleganCarpenter } from '../../engine/hydro/morison.js';
import { ocimfForce } from '../../engine/hydro/drift.js';
import { pileScour } from '../../engine/geotech/scour.js';
import { need, run, fmt, round, cap, adv, check, linspace, CONFIRM, NOT_HIFI, ECON_BASIS } from './_common.js';

/** Morison force and mudline moment on one vertical pile through a wave cycle. Returns arrays over phase. */
function morisonCycle({ H, T, d, De, Cd, Cm, Uc, rho, nSeg, theory, nPhase = 36 }) {
  const useS5 = theory === 'stokes5';
  const sol = useS5 ? run(() => stokes5Wave({ H, d, T }), 'Stokes fifth-order wave') : null;
  const out = { phase: [], F: [], M: [], eta: [], uMax: 0 };
  for (let ip = 0; ip < nPhase; ip++) {
    const t = T * ip / nPhase;
    const eta = useS5 ? stokes5Kinematics(sol, { x: 0, z: 0, t }).eta : airyKinematics({ H, T, d, x: 0, z: 0, t, clipAboveSurface: false }).eta;
    const top = eta, dz = (top + d) / nSeg;
    let F = 0, M = 0;
    for (let i = 0; i < nSeg; i++) {
      const z = -d + (i + 0.5) * dz;
      const k = useS5 ? stokes5Kinematics(sol, { x: 0, z: Math.min(z, eta - 1e-9), t }) : airyKinematics({ H, T, d, x: 0, z, t, stretching: 'wheeler', rho });
      // current: 1/7 power profile stretched with the surface (value at the still-water-equivalent elevation)
      const zc = Math.min(0, (z - eta) * d / (d + eta));
      const u = k.u + (Uc > 0 ? currentPowerLaw(zc, { v0: Uc, d }) : 0);
      const f = 0.5 * rho * Cd * De * u * Math.abs(u) + rho * Cm * Math.PI * De * De / 4 * k.ax;   // Morison, N/m
      F += f * dz; M += f * dz * (z + d);
      if (Math.abs(u) > out.uMax) out.uMax = Math.abs(u);
    }
    out.phase.push(360 * ip / nPhase); out.F.push(F); out.M.push(M); out.eta.push(eta);
  }
  return out;
}

export default {
  id: 's06-environmental-loading', suite: 6, title: 'Wave, current, wind and vessel-induced loading',
  summary: 'Extreme-storm wave, current and wind load on the dolphin and the operating load from a moored vessel, with base shear, overturning moment, wave clearance under the cap and storm stability.',
  covers: [8], fidelity: 'analytical',
  standards: ['DNV-RP-C205 (2019) §3, §4, §6', 'ISO 19901-1:2015 §8–9', 'API RP 2A-WSD 22nd ed. §5.3', 'OCIMF MEG4 (2018) Appendix A (vessel force definitions)', 'BS 6349-1-2:2016'],
  models: ['Morison equation', 'Airy wave theory with Wheeler stretching', 'Stokes fifth-order wave theory', 'Power-law current profile', 'Wind drag on the superstructure', 'OCIMF-form vessel wind and current forces', 'Pile scour (Sumer & Fredsøe)'],
  limits: 'Slender vertical piles without diffraction (diameter below 0.2 wavelength); piles load in phase with a group shielding factor; no wave slam on the cap (checked as clearance only) and no breaking-wave impact. Vessel coefficients are user data: the defaults are generic magnitudes, not OCIMF curves. Wave-drift and passing-ship forces enter as a specified additional load. ' + NOT_HIFI,
  publishes: ['dolphin.stormShear', 'dolphin.stormMoment', 'dolphin.vesselLoad', 'dolphin.scourDepth'],
  searchTerms: 'Morison wave load pile dolphin current wind vessel mooring force',
  inputs: [
    { key: 'nPiles', label: 'Number of piles', unit: '-', type: 'int', default: 4, min: 1, max: 60, group: 'Structure', role: 'geometry' },
    { key: 'D', label: 'Pile outer diameter', unit: 'm', type: 'number', default: 1.4, min: 0.2, max: 6, group: 'Structure', role: 'geometry', bind: 'dolphin.pileD', design: true },
    { key: 'growth', label: 'Marine growth thickness', unit: 'mm', type: 'number', default: 50, min: 0, max: 300, group: 'Structure', role: 'ic', help: 'Adds twice its thickness to the hydrodynamic diameter and roughens the surface.' },
    { key: 'shield', label: 'Group shielding / phase factor', unit: '-', type: 'number', default: 1.0, min: 0.4, max: 1.2, group: 'Structure', role: 'load', help: '1.0 sums all piles in phase (conservative for closely spaced groups).' },
    { key: 'capSoffit', label: 'Cap soffit above still-water level', unit: 'm', type: 'number', default: 4.5, min: 0, max: 30, group: 'Structure', role: 'geometry' },
    { key: 'Awind', label: 'Projected wind area of cap and equipment', unit: 'm²', type: 'number', default: 30, min: 0, max: 2000, group: 'Structure', role: 'geometry' },
    { key: 'zWind', label: 'Centre of wind area above still water', unit: 'm', type: 'number', default: 6.5, min: 0, max: 60, group: 'Structure', role: 'geometry' },
    { key: 'depth', label: 'Still-water depth', unit: 'm', type: 'number', default: 15, min: 2, max: 120, group: 'Site', role: 'environment', bind: 'site.depth' },
    { key: 'surge', label: 'Tide plus storm surge above datum', unit: 'm', type: 'number', default: 1.5, min: -5, max: 12, group: 'Site', role: 'bc', bind: 'site.surge100', help: 'Added to the depth for the storm case.' },
    { key: 'rhoW', label: 'Water density', unit: 'kg/m³', type: 'number', default: 1025, min: 990, max: 1050, group: 'Site', role: 'environment', bind: 'site.rhoW' },
    { key: 'Hmax', label: 'Extreme wave height', unit: 'm', type: 'number', default: 4.0, min: 0, max: 35, group: 'Extreme storm (no vessel)', role: 'bc', bind: 'site.Hmax100', uncertain: { dist: 'gumbel', cov: 0.12 } },
    { key: 'T', label: 'Associated wave period', unit: 's', type: 'number', default: 8, min: 2, max: 25, group: 'Extreme storm (no vessel)', role: 'bc', bind: 'site.Tp100' },
    { key: 'Uc', label: 'Surface current', unit: 'm/s', type: 'number', default: 1.2, min: 0, max: 5, group: 'Extreme storm (no vessel)', role: 'bc', bind: 'site.Uc100', uncertain: { dist: 'normal', cov: 0.15 } },
    { key: 'Uw', label: 'Wind speed at 10 m (3 s gust)', unit: 'm/s', type: 'number', default: 40, min: 0, max: 90, group: 'Extreme storm (no vessel)', role: 'bc', bind: 'site.Uw100' },
    { key: 'theory', label: 'Wave theory', type: 'select', default: 'auto', group: 'Hydrodynamics', role: 'control', options: [{ value: 'auto', label: 'Automatic (Airy or Stokes 5th by steepness and depth)' }, { value: 'airy', label: 'Airy with Wheeler stretching' }, { value: 'stokes5', label: 'Stokes fifth order' }] },
    { key: 'Cd', label: 'Drag coefficient', unit: '-', type: 'number', default: 1.05, min: 0.4, max: 2, group: 'Hydrodynamics', role: 'load', uncertain: { dist: 'normal', cov: 0.1 }, help: `0.65 smooth, 1.05 rough (marine growth) for circular members (API RP 2A / ISO 19902). ${CONFIRM}` },
    { key: 'Cm', label: 'Inertia coefficient', unit: '-', type: 'number', default: 1.2, min: 1, max: 2.2, group: 'Hydrodynamics', role: 'load', help: `1.6 smooth, 1.2 rough (API RP 2A / ISO 19902); 2.0 potential-flow value. ${CONFIRM}` },
    { key: 'CsWind', label: 'Wind shape coefficient', unit: '-', type: 'number', default: 1.5, min: 0.5, max: 2.5, group: 'Hydrodynamics', role: 'load', help: 'API RP 2A-WSD §5.3.2: 1.5 for beams and sides of buildings, 1.0 overall platform projected area, 0.5 cylinders.' },
    { key: 'nSeg', label: 'Integration segments over the water column', unit: '-', type: 'int', default: 40, min: 5, max: 640, group: 'Hydrodynamics', role: 'control' },
    { key: 'UwOp', label: 'Wind speed limit with a vessel alongside', unit: 'm/s', type: 'number', default: 25, min: 0, max: 60, group: 'Operating case (vessel moored)', role: 'bc', help: 'Departure criterion of the terminal; the vessel is assumed to leave before the extreme storm.' },
    { key: 'UcOp', label: 'Current with a vessel alongside', unit: 'm/s', type: 'number', default: 1.0, min: 0, max: 4, group: 'Operating case (vessel moored)', role: 'bc' },
    { key: 'Hop', label: 'Wave height with a vessel alongside', unit: 'm', type: 'number', default: 1.5, min: 0, max: 10, group: 'Operating case (vessel moored)', role: 'bc' },
    { key: 'LBP', label: 'Vessel length between perpendiculars', unit: 'm', type: 'number', default: 220, min: 20, max: 420, group: 'Operating case (vessel moored)', role: 'load' },
    { key: 'draft', label: 'Vessel draft', unit: 'm', type: 'number', default: 12.5, min: 1, max: 26, group: 'Operating case (vessel moored)', role: 'load' },
    { key: 'AL', label: 'Vessel broadside windage', unit: 'm²', type: 'number', default: 3200, min: 20, max: 20000, group: 'Operating case (vessel moored)', role: 'load' },
    { key: 'CyW', label: 'Vessel lateral wind coefficient (user data)', unit: '-', type: 'number', default: 0.9, min: 0, max: 2, group: 'Operating case (vessel moored)', role: 'load', help: 'Generic magnitude, not OCIMF data: enter the coefficient for the vessel, loading condition and wind angle.' },
    { key: 'CyC', label: 'Vessel lateral current coefficient (user data)', unit: '-', type: 'number', default: 1.2, min: 0, max: 6, group: 'Operating case (vessel moored)', role: 'load', help: 'Includes the shallow-water increase for the depth-to-draft ratio; from the coefficient source.' },
    { key: 'angC', label: 'Current angle off the bow', unit: 'deg', type: 'number', default: 10, min: 0, max: 90, group: 'Operating case (vessel moored)', role: 'bc', help: 'The lateral current force is taken as proportional to sin of this angle (user coefficient at beam-on).' },
    { key: 'Fextra', label: 'Additional vessel force (wave drift, passing ship)', unit: 'kN', type: 'number', default: 150, min: 0, max: 1e5, group: 'Operating case (vessel moored)', role: 'bc', bind: 'vessel.Fmean' },
    { key: 'share', label: 'Share of the vessel force on this dolphin', unit: '-', type: 'number', default: 0.5, min: 0.05, max: 1, group: 'Operating case (vessel moored)', role: 'bc', help: '0.5 for two breasting dolphins loaded evenly; more for yawed vessels.' },
    { key: 'zVessel', label: 'Vessel contact / line height above still water', unit: 'm', type: 'number', default: 3.5, min: -5, max: 30, group: 'Operating case (vessel moored)', role: 'bc' },
    { key: 'Hcap', label: 'Lateral capacity of the dolphin', unit: 'kN', type: 'number', default: 3500, min: 1, max: 1e6, group: 'Acceptance', role: 'criteria', bind: 'dolphin.lateralCapacity', help: 'Ultimate horizontal resistance from the pile-capacity or integrity module (sum over the piles).' },
    { key: 'fs', label: 'Required factor against lateral failure', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Acceptance', role: 'criteria', help: CONFIRM },
  ],
  mesh: { key: 'nSeg', label: 'Integration segments', levels: [10, 20, 40, 80], qoi: 'Fwave' },
  compute(inp) {
    const d = inp.depth + inp.surge, De = inp.D + 2 * inp.growth / 1e3, rho = inp.rhoW, n = inp.nPiles * inp.shield;
    need(d > 1, 'The water depth including surge must exceed 1 m.');
    const rec = inp.Hmax > 0 ? selectWaveTheory({ H: inp.Hmax, T: inp.T, d }) : { theory: 'airy', HoverHb: 0, ursell: 0, dOverL: 0, steepness: 0 };
    need(rec.theory !== 'breaking', `A ${fmt(inp.Hmax)} m, ${fmt(inp.T)} s wave breaks in ${fmt(d)} m of water (H/Hb = ${fmt(rec.HoverHb)}): use the depth-limited wave height (about ${fmt(inp.Hmax / rec.HoverHb)} m) or a longer period.`);
    let theory = inp.theory === 'auto' ? (rec.theory === 'airy' ? 'airy' : 'stokes5') : inp.theory;
    let s5note = '';
    if (theory === 'stokes5') { try { const s = stokes5Wave({ H: inp.Hmax, d, T: inp.T }); if (!s.valid) { if (inp.theory === 'auto') { theory = 'airy'; s5note = 'Stokes fifth-order theory is outside its validity range here (Ursell number above about 40): Airy kinematics with Wheeler stretching are used; a stream-function solution is recommended.'; } else s5note = 'Stokes fifth-order theory is outside its validity range here: treat the result with caution.'; } } catch { theory = 'airy'; s5note = 'The Stokes fifth-order solution failed; Airy kinematics with Wheeler stretching are used.'; } }
    const cyc = inp.Hmax > 0 ? morisonCycle({ H: inp.Hmax, T: inp.T, d, De, Cd: inp.Cd, Cm: inp.Cm, Uc: inp.Uc, rho, nSeg: inp.nSeg, theory })
      : morisonCycle({ H: 1e-6, T: inp.T, d, De, Cd: inp.Cd, Cm: inp.Cm, Uc: inp.Uc, rho, nSeg: inp.nSeg, theory: 'airy', nPhase: 2 });
    let iMax = 0; cyc.F.forEach((f, i) => { if (Math.abs(f) > Math.abs(cyc.F[iMax])) iMax = i; });
    const F1 = Math.abs(cyc.F[iMax]), M1 = Math.abs(cyc.M[iMax]);
    const crest = Math.max(...cyc.eta);
    const cf = inp.Hmax > 0 ? pileForceClosedForm({ H: inp.Hmax, T: inp.T, depth: d, D: De, Cd: inp.Cd, Cm: inp.Cm, rho }) : { Fmax: 0, Mmax: 0, FD: 0, FI: 0 };
    // current alone (for the split) by the same integration with a vanishing wave
    const curOnly = inp.Uc > 0 ? Math.abs(morisonCycle({ H: 1e-6, T: inp.T, d, De, Cd: inp.Cd, Cm: inp.Cm, Uc: inp.Uc, rho, nSeg: inp.nSeg, theory: 'airy', nPhase: 1 }).F[0]) : 0;
    // wind on the superstructure, speed referred to the centre of area with a 1/8 power law (API RP 2A gust profile, approximate)
    const Uz = inp.Uw * Math.max(inp.zWind / 10, 0.1) ** 0.125;
    const Fw = inp.Awind > 0 && inp.Uw > 0 ? windForce({ U: Uz, A: inp.Awind, Cs: inp.CsWind }).F : 0;
    const Fstorm = n * F1 + Fw, Mstorm = n * M1 + Fw * (d + inp.zWind);
    // operating case: vessel wind (beam-on) and current, share on this dolphin, plus operating wave on the piles
    const vw = inp.UwOp > 0 ? ocimfForce({ kind: 'wind', V: inp.UwOp, angleDeg: 90, Cx: 0, Cy: inp.CyW, AT: 1, AL: inp.AL, L: inp.LBP }).Fy : 0;
    const vc = inp.UcOp > 0 ? ocimfForce({ kind: 'current', V: inp.UcOp, angleDeg: inp.angC, Cx: 0, Cy: inp.CyC * Math.sin(inp.angC * Math.PI / 180), L: inp.LBP, T: inp.draft, rho }).Fy : 0;
    const Fves = inp.share * (Math.abs(vw) + Math.abs(vc) + inp.Fextra * 1e3);
    const dOp = inp.depth;
    let Fop1 = 0, Mop1 = 0;
    if (inp.Hop > 0 || inp.UcOp > 0) {
      const recOp = inp.Hop > 0 ? selectWaveTheory({ H: inp.Hop, T: Math.min(inp.T, 7), d: dOp }) : { theory: 'airy' };
      if (recOp.theory !== 'breaking') { const c2 = morisonCycle({ H: Math.max(inp.Hop, 1e-6), T: Math.min(inp.T, 7), d: dOp, De, Cd: inp.Cd, Cm: inp.Cm, Uc: inp.UcOp, rho, nSeg: Math.max(10, Math.round(inp.nSeg / 2)), theory: 'airy', nPhase: 18 }); let j = 0; c2.F.forEach((f, i) => { if (Math.abs(f) > Math.abs(c2.F[j])) j = i; }); Fop1 = Math.abs(c2.F[j]); Mop1 = Math.abs(c2.M[j]); }
    }
    const Fop = Fves + n * Fop1, Mop = Fves * (dOp + inp.zVessel) + n * Mop1;
    const Fgov = Math.max(Fstorm, Fop), govName = Fstorm >= Fop ? 'extreme storm without vessel' : 'operating limit with vessel moored';
    const KC = cyc.uMax > 0 && inp.Hmax > 0 ? keuleganCarpenter(Math.max(cyc.uMax - inp.Uc, 1e-6), inp.T, De) : 0;
    const Um = Math.max(cyc.uMax - inp.Uc, 0);
    const sc = pileScour({ KC: inp.Hmax > 0 ? Math.max(KC, 1e-6) : Infinity, Ucw: inp.Uc + Um > 0 ? inp.Uc / (inp.Uc + Um) : 0 });
    const scour = sc.SoverD * inp.D;
    const fsAch = inp.Hcap * 1e3 / Math.max(Fgov, 1);
    const clearance = inp.capSoffit - (crest + inp.surge);
    const checks = [
      { name: `Storm stability: lateral capacity over governing load (${govName})`, demand: round(inp.fs, 4), capacity: round(cap(fsAch), 4), unit: 'FS (required / achieved)', uc: Math.min(inp.fs / fsAch, 99), ref: 'ISO 19902:2020 §9 / API RP 2A-WSD §6.3.4 (overridable factor)', explain: 'Ultimate lateral resistance of the dolphin against the larger of the storm and moored-vessel loads. Below the required factor the piles hinge or the soil fails in the design event.' },
      check('Wave crest against the cap soffit', crest + inp.surge, inp.capSoffit, 'm', 'ISO 19902:2020 §A.9.4 (air gap); BS 6349-1-2:2016', 'Crest elevation including surge against the underside of the cap. If the crest reaches the cap, wave-in-deck slam and uplift occur — loads this module does not compute.'),
    ];
    const w = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = w.uc <= 1;
    // diameter sensitivity: drag ∝ D, inertia ∝ D²
    const dF = n * (cf.FD * 0.1 + cf.FI * 0.21);
    const zs = linspace(-d, 0, 21);
    const uCrest = zs.map((z) => { const k = airyKinematics({ H: Math.max(inp.Hmax, 1e-6), T: inp.T, d, z, t: 0, stretching: 'none' }); return round(k.u + (inp.Uc > 0 ? currentPowerLaw(z, { v0: inp.Uc, d }) : 0)); });
    const advice = [
      adv('engineering', ok ? 'ok' : 'fail', `Governing horizontal load ${fmt(Fgov / 1e3)} kN (${govName}); storm: ${fmt(n * F1 / 1e3)} kN wave and current on ${inp.nPiles} piles plus ${fmt(Fw / 1e3)} kN wind; vessel alongside: ${fmt(Fop / 1e3)} kN. Achieved factor ${fmt(cap(fsAch))} against ${fmt(inp.fs)}.`,
        ok ? `Carry both cases into the global strength check with their lever arms (${fmt(Mstorm / Math.max(Fstorm, 1))} m and ${fmt(Mop / Math.max(Fop, 1))} m above the mudline). Replace the generic vessel coefficients with project values before design.`
          : (w.name.startsWith('Wave crest') ? `Raise the cap soffit by at least ${fmt(-clearance + 0.5)} m or design the cap for wave-in-deck loading.` : `Provide a lateral capacity of at least ${fmt(inp.fs * Fgov / 1e3)} kN, or lower the wind limit for vessels alongside (vessel wind force scales with speed squared).`),
        'DNV-RP-C205 (2019) §6; ISO 19901-1:2015; API RP 2A-WSD §5.3; OCIMF MEG4 (2018)',
        `Wave load grows with pile diameter (drag ∝ D, inertia ∝ D²): 10 % more diameter adds about ${fmt(dF / 1e3)} kN of storm load. Cleaning marine growth of ${fmt(inp.growth)} mm restores ${fmt((1 - inp.D / De) * 100)} % of drag area without any new steel.`),
      adv('economic', 'info', `The departure wind limit of ${fmt(inp.UwOp)} m/s produces ${fmt(Math.abs(vw) * inp.share / 1e3)} kN on this dolphin; estimated scour ${fmt(scour)} m (${fmt(scour / inp.D)} D).`,
        'Trade the wind limit against berth availability: each 2 m/s on the departure criterion changes the vessel wind load by about 15–20 % and the hours of downtime per year from the wind climate. Budget scour protection or the extra free length in the pile design.',
        `Berth downtime assessment (PIANC WG153); ${ECON_BASIS}`),
    ];
    return {
      kpis: [
        { key: 'Fwave', label: 'Wave and current force on one pile', unit: 'kN', value: F1 / 1e3, explain: 'Peak in-line Morison force on a single pile through the wave cycle, integrated from the seabed to the instantaneous surface with current.' },
        { key: 'Mwave', label: 'Mudline moment on one pile', unit: 'kN·m', value: M1 / 1e3, explain: 'Overturning moment about the mudline at the instant of peak force on one pile.' },
        { key: 'FwaveLin', label: 'Closed-form linear wave force (no current)', unit: 'kN', value: cf.Fmax / 1e3, explain: 'Textbook Morison result for an Airy wave integrated to still-water level, without current. The difference from the value above is the effect of current, crest kinematics and wave nonlinearity.' },
        { key: 'Fstorm', label: 'Storm base shear on the dolphin', unit: 'kN', value: Fstorm / 1e3, explain: 'All piles (with the group factor) plus wind on the superstructure in the extreme storm, vessel absent.' },
        { key: 'Mstorm', label: 'Storm overturning moment at the mudline', unit: 'kN·m', value: Mstorm / 1e3, explain: 'Moment of the storm loads about the mudline; resisted by pile bending and push–pull.' },
        { key: 'Fwind', label: 'Wind force on the superstructure', unit: 'kN', value: Fw / 1e3, explain: 'Drag on cap, fender panel, hooks and walkway ends at the extreme gust speed.' },
        { key: 'Fvessel', label: 'Moored-vessel load on this dolphin', unit: 'kN', value: Fves / 1e3, explain: 'Share of the vessel wind, current and additional force carried by this dolphin at the operating limit.' },
        { key: 'Fop', label: 'Operating-case base shear', unit: 'kN', value: Fop / 1e3, explain: 'Vessel load plus the operating wave and current on the piles.' },
        { key: 'Fgov', label: 'Governing horizontal load', unit: 'kN', value: Fgov / 1e3, explain: 'Larger of the storm and operating cases: the environmental design load to compare with the berthing load.' },
        { key: 'fsLat', label: 'Factor against lateral failure', unit: '-', value: cap(fsAch), limit: inp.fs, sense: 'min', explain: 'Lateral capacity of the dolphin over the governing environmental load.' },
        { key: 'crest', label: 'Crest elevation above datum', unit: 'm', value: crest + inp.surge, limit: inp.capSoffit, sense: 'max', explain: 'Wave crest height including tide and surge. It must stay below the cap soffit to avoid wave-in-deck loading.' },
        { key: 'KC', label: 'Keulegan–Carpenter number', unit: '-', value: KC, explain: 'Orbital excursion over pile diameter. Below about 6 the load is inertia-dominated and there is little wave scour; above about 20 drag dominates.' },
        { key: 'scour', label: 'Estimated local scour depth', unit: 'm', value: scour, explain: 'Equilibrium scour at a slender pile under the storm waves and current (up to 1.3 D in steady current). It lengthens the free pile and should be carried into the pile analyses.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Wave force on one pile through the wave cycle', x: { label: 'Phase', unit: 'deg' }, y: { label: 'Force (kN) / surface (m × 10)', unit: '' },
          series: [{ name: 'In-line force (kN)', x: cyc.phase, y: cyc.F.map((f) => round(f / 1e3)) }, { name: 'Surface elevation × 10 (m)', x: cyc.phase, y: cyc.eta.map((e) => round(e * 10)) }],
          explain: 'Force history over one wave. A peak under the crest means drag (and current) dominates; a peak near the zero crossing means inertia dominates.' },
        { type: 'line', title: 'Horizontal velocity under the crest (wave plus current)', x: { label: 'Velocity', unit: 'm/s' }, y: { label: 'Elevation', unit: 'm' }, series: [{ name: 'u', x: uCrest, y: zs.map((z) => round(z)) }],
          explain: 'Particle velocity profile below still-water level at the crest phase (linear wave plus current). The upper few metres carry most of the drag.' },
        { type: 'bar', title: 'Horizontal load by source', categories: ['Storm: wave + current', 'Storm: current alone', 'Storm: wind on structure', 'Operating: vessel', 'Operating: wave + current'],
          series: [{ name: 'Force (kN)', y: [n * F1, n * curOnly, Fw, Fves, n * Fop1].map((v) => round(v / 1e3)) }], y: { label: 'Force', unit: 'kN' },
          explain: 'Contribution of each source to the two load cases. It shows whether the dolphin is governed by the storm on the bare structure or by the vessel at the operating limit.' },
      ],
      tables: [{ title: 'Load cases at the mudline', columns: [{ key: 'case', label: 'Case' }, { key: 'F', label: 'Base shear', unit: 'kN' }, { key: 'M', label: 'Overturning moment', unit: 'kN·m' }, { key: 'arm', label: 'Lever arm', unit: 'm' }],
        rows: [['Extreme storm, no vessel', Fstorm, Mstorm], ['Operating limit, vessel moored', Fop, Mop]].map(([name, F, M]) => ({ case: name, F: round(F / 1e3), M: round(M / 1e3), arm: round(M / Math.max(F, 1)) })),
        explain: 'Resultant horizontal force and moment about the mudline for each case, ready to apply in the structural and pile modules.' }],
      notes: [`Wave theory used: ${theory === 'stokes5' ? 'Stokes fifth order' : 'Airy with Wheeler stretching'} (chart recommendation: ${rec.theory}; H/Hb = ${fmt(rec.HoverHb)}).`, ...(s5note ? [s5note] : []), `${govName} governs.`, 'Vessel coefficients are user-supplied magnitudes (generic defaults).'],
      advice,
      publish: { 'dolphin.stormShear': Fstorm, 'dolphin.stormMoment': Mstorm, 'dolphin.vesselLoad': Fves, 'dolphin.scourDepth': scour },
      quantities: {},
    };
  },
  verify: [
    // Inertia-dominated pile in deep water: k = ω²/g = 0.2516 1/m, tanh(kd) = 1;
    // F_I = ρ·Cm·(πD²/4)·g·(H/2) = 1025·2·3.1416·9.80665·0.5 = 31.58 kN; F_D ≈ 1.26 kN < F_I/2 so Fmax = F_I;
    // lever arm d − 1/k = 50 − 3.975 = 46.03 m → M = 1453 kN·m.
    { name: 'Inertia force on a vertical pile, deep water: F = ρ·Cm·(πD²/4)·g·H/2', source: 'Morison et al. (1950) with Airy kinematics; Dean & Dalrymple (1991) §8 (hand calculation in the module source)',
      inputs: { nPiles: 1, D: 2, growth: 0, shield: 1, depth: 50, surge: 0, Hmax: 1, T: 4, Uc: 0, Uw: 0, Cd: 1, Cm: 2, theory: 'airy', nSeg: 160, capSoffit: 10 },
      expect: [{ kpi: 'Fwave', value: 31.58, tol: 0.01 }, { kpi: 'Mwave', value: 1453, tol: 0.015 }, { kpi: 'FwaveLin', value: 31.58, tol: 0.002 }] },
    // Steady current only on one pile: F = ½ρCdD∫u²dz with u = U(1 + z/d)^(1/7): ∫ = U²·d·7/9 → 0.5·1025·1·1·1²·20·7/9 = 7.972 kN.
    { name: 'Current drag with a 1/7 power profile: F = ½ρCdD·U²·d·7/9', source: 'Closed-form integral of the power-law profile (DNV-RP-C205 §4.1.4)',
      inputs: { nPiles: 1, D: 1, growth: 0, shield: 1, depth: 20, surge: 0, Hmax: 0, T: 8, Uc: 1, Uw: 0, Cd: 1, Cm: 2, nSeg: 160, capSoffit: 10 },
      expect: [{ kpi: 'Fwave', value: 7.972, tol: 0.01 }] },
  ],
};
