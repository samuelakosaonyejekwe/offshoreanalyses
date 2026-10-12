// Seismic response of a marine structure site: equivalent-linear ground response through the soil column,
// surface response spectrum, Mononobe–Okabe and Westergaard thrusts on a wall, Newmark sliding displacement
// and SPT-based liquefaction triggering.
import { equivalentLinearResponse, transferFunction, responseSpectrum, ariasIntensity } from '../../engine/geotech/siteresponse.js';
import { mononobeOkabe, westergaard, rankine } from '../../engine/geotech/earthpressure.js';
import { sptProfile, liquefactionPotentialIndex } from '../../engine/geotech/liquefaction.js';
import { newmarkDisplacement } from '../../engine/geotech/slope.js';
import { G, DEG, CO2_ROCK, need, run, finite, r, f, linspace, thin, check, sfCheck, eng, eco, slrInputs, SLR_BASIS } from './_common.js';

/** Deterministic synthetic accelerogram: twelve sines around the dominant frequency with fixed phases and a build-up / decay envelope. */
function synthetic(pga, fDom, dur, dt) {
  const n = Math.round(dur / dt), a = new Float64Array(n);
  let seed = 12345; const rnd = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
  const comps = Array.from({ length: 12 }, (_, i) => ({ f: fDom * (0.35 + 1.9 * i / 11), ph: 2 * Math.PI * rnd(), w: Math.exp(-(((i - 4) / 4) ** 2)) }));
  const t1 = 0.15 * dur, t2 = 0.5 * dur;
  let peak = 0;
  for (let i = 0; i < n; i++) {
    const t = i * dt, env = t < t1 ? (t / t1) ** 2 : t < t2 ? 1 : Math.exp(-3 * (t - t2) / (dur - t2));
    let s = 0; for (const c of comps) s += c.w * Math.sin(2 * Math.PI * c.f * t + c.ph);
    a[i] = env * s; peak = Math.max(peak, Math.abs(a[i]));
  }
  for (let i = 0; i < n; i++) a[i] *= pga / peak;
  return a;
}

export default {
  id: 's07-seismic-response', suite: 7, title: 'Seismic response of the site and wall',
  summary: 'Carries the bedrock motion to the surface, then checks a gravity wall for seismic earth and water thrust and sliding displacement, and the foundation soil for liquefaction.',
  covers: [8], fidelity: 'reduced-order',
  standards: ['EN 1998-5:2004 §7 and Annex E (retaining structures), Annex B (liquefaction)', 'PIANC WG 34 Seismic Design Guidelines for Port Structures (2001)', 'ISO 19901-2:2022 (seismic design procedures)', 'Youd et al. (2001) NCEER liquefaction procedure', 'Kramer, Geotechnical Earthquake Engineering (1996) ch. 7, 11'],
  models: ['Equivalent-linear 1-D site response (SHAKE)', 'Darendeli modulus reduction and damping', 'Elastic response spectrum', 'Mononobe–Okabe seismic earth pressure', 'Westergaard hydrodynamic pressure', 'Newmark sliding block', 'SPT liquefaction triggering (CSR / CRR)', 'Liquefaction potential index'],
  limits: 'Vertically propagating shear waves in horizontal layers (total stress, no pore-pressure build-up in the response), pseudo-static wall thrusts and a rigid sliding block; basin effects, wall flexibility, tilting, and post-liquefaction flow are not modelled. The built-in record is synthetic: site-specific spectrum-matched records and, for important structures, effective-stress dynamic finite-element analysis are needed and are outside the module.',
  publishes: ['marine.pgaSurface', 'marine.seismicKh', 'marine.seismicSlide'],
  searchTerms: 'seismic site response Mononobe-Okabe quay wall Newmark sliding liquefaction',
  inputs: [
    { key: 'record', label: 'Bedrock outcrop accelerogram (optional)', type: 'series', default: { x: [], y: [] }, xLabel: 'Time', xUnit: 's', yLabel: 'Acceleration', yUnit: 'm/s²', group: 'Ground motion at the foundation boundary', role: 'bc', help: 'Evenly sampled record. Leave empty to use the built-in synthetic record scaled to the peak acceleration below.' },
    { key: 'pga', label: 'Peak bedrock acceleration', unit: 'm/s²', type: 'number', default: 1.5, min: 0.1, max: 12, group: 'Ground motion at the foundation boundary', role: 'bc', bind: 'site.pga475', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'For the design return period (e.g. 475 years). An imported record is scaled to this value.' },
    { key: 'fDom', label: 'Dominant frequency of the synthetic record', unit: 'Hz', type: 'number', default: 1.5, min: 0.3, max: 10, group: 'Ground motion at the foundation boundary', role: 'bc' },
    { key: 'dur', label: 'Duration of the synthetic record', unit: 's', type: 'number', default: 20, min: 5, max: 80, group: 'Ground motion at the foundation boundary', role: 'bc' },
    { key: 'Mw', label: 'Moment magnitude', unit: '-', type: 'number', default: 7.0, min: 5, max: 9, group: 'Ground motion at the foundation boundary', role: 'bc' },
    { key: 'column', label: 'Soil column (top to bottom)', type: 'table', group: 'Ground', role: 'soil',
      columns: [{ key: 'name', label: 'Layer', type: 'text' }, { key: 'h', label: 'Thickness', unit: 'm', type: 'number' }, { key: 'rho', label: 'Density', unit: 'kg/m³', type: 'number' }, { key: 'Vs', label: 'Shear-wave velocity', unit: 'm/s', type: 'number' }, { key: 'PI', label: 'Plasticity index', unit: '%', type: 'number' }, { key: 'N160', label: '(N1)60', unit: '-', type: 'number' }, { key: 'FC', label: 'Fines', unit: '%', type: 'number' }],
      default: [{ name: 'Compacted fill', h: 5, rho: 1900, Vs: 180, PI: 0, N160: 24, FC: 8 }, { name: 'Dense sand', h: 8, rho: 1950, Vs: 260, PI: 0, N160: 28, FC: 5 }, { name: 'Stiff clay', h: 12, rho: 2000, Vs: 330, PI: 30, N160: 0, FC: 100 }],
      help: '(N1)60 = 0 marks a layer not susceptible to liquefaction (clay, rock).' },
    { key: 'rockVs', label: 'Bedrock shear-wave velocity', unit: 'm/s', type: 'number', default: 900, min: 300, max: 100000, group: 'Ground', role: 'bc' },
    { key: 'rockRho', label: 'Bedrock density', unit: 'kg/m³', type: 'number', default: 2400, min: 1800, max: 3000, group: 'Ground', role: 'bc' },
    { key: 'nonlinear', label: 'Strain-compatible (equivalent-linear) soil', type: 'bool', default: true, group: 'Ground', role: 'soil', help: 'Off: linear soil with 3 % damping.' },
    { key: 'Hw', label: 'Wall height', unit: 'm', type: 'number', default: 10, min: 1, max: 40, group: 'Wall', role: 'geometry' },
    { key: 'W', label: 'Effective (buoyant) weight of wall and fill on the base', unit: 'kN/m', type: 'number', default: 1500, min: 10, max: 50000, group: 'Wall', role: 'ic', design: true, help: 'Gravity state before the earthquake; from the static stability module.' },
    { key: 'mu', label: 'Base friction coefficient', unit: '-', type: 'number', default: 0.6, min: 0.2, max: 0.9, group: 'Wall', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 } },
    { key: 'phiDeg', label: 'Backfill friction angle', unit: 'deg', type: 'number', default: 36, min: 25, max: 45, group: 'Wall', role: 'soil', uncertain: { dist: 'normal', cov: 0.07 } },
    { key: 'gam', label: 'Backfill effective unit weight (average)', unit: 'kN/m³', type: 'number', default: 13, min: 6, max: 22, group: 'Wall', role: 'soil' },
    { key: 'hw', label: 'Water depth in front of the wall', unit: 'm', type: 'number', default: 7, min: 0, max: 40, group: 'Wall', role: 'bc', bind: 'marine.hToe' },
    ...slrInputs({ group: 'Wall' }).slice(0, 1),
    { key: 'khRatio', label: 'Seismic coefficient / (surface PGA/g)', unit: '-', type: 'number', default: 0.5, min: 0.2, max: 1, group: 'Criteria', role: 'criteria', help: 'EN 1998-5 Table 7.1: 1/r with r = 2 for gravity walls that can slide, 1 for walls that cannot move.' },
    { key: 'dAllowPct', label: 'Allowable sliding (% of wall height)', unit: '%', type: 'number', default: 1.5, min: 0.1, max: 10, group: 'Criteria', role: 'criteria', help: 'PIANC (2001) degree I serviceable: below 1.5 % of height for gravity quay walls.' },
    { key: 'Tstr', label: 'Natural period of the structure', unit: 's', type: 'number', default: 0.8, min: 0.05, max: 6, group: 'Structure', role: 'geometry' },
    { key: 'mass', label: 'Seismic mass of the structure', unit: 't', type: 'number', default: 5000, min: 1, max: 1e6, group: 'Structure', role: 'geometry' },
  ],
  compute(inp) {
    const col = (inp.column || []).map((l, i) => ({ name: String(l.name || `Layer ${i + 1}`), h: Number(l.h), rho: Number(l.rho), Vs: Number(l.Vs), PI: Number(l.PI) || 0, N160: Number(l.N160) || 0, FC: Number(l.FC) || 0 }));
    need(col.length >= 1, 'Give at least one soil layer above bedrock.');
    for (const l of col) need(l.h > 0 && l.rho > 1000 && l.Vs > 30, `Layer '${l.name}': thickness, density (> 1000 kg/m³) and shear-wave velocity (> 30 m/s) are required.`);
    // input motion
    let acc, dt;
    const rec = inp.record;
    if (rec && rec.x && rec.x.length >= 64) {
      dt = (rec.x[rec.x.length - 1] - rec.x[0]) / (rec.x.length - 1);
      need(dt > 0 && dt < 0.1, 'The accelerogram must be evenly sampled with a time step below 0.1 s.');
      const pk = Math.max(...rec.y.map((v) => Math.abs(v)));
      need(pk > 0, 'The accelerogram is all zeros.');
      acc = Float64Array.from(rec.y, (v) => v * inp.pga / pk);
      need(acc.length <= 16384, 'The accelerogram is longer than 16384 samples: decimate or trim it.');
    } else { dt = 0.01; acc = synthetic(inp.pga, inp.fDom, inp.dur, dt); }
    // layers: Darendeli curves at the mid-layer mean effective stress (K0 = 0.5, submerged)
    let sv = 0;
    const layers = col.map((l) => {
      const gEff = (l.rho - 1025) * G, sMid = sv + gEff * l.h / 2; sv += gEff * l.h;
      return { h: l.h, rho: l.rho, Vs: l.Vs, damping: 0.03, ...(inp.nonlinear ? { curve: { model: 'darendeli', sigmaM: Math.max(sMid * 2 / 3, 5000), PI: l.PI } } : {}) };
    });
    const bedrock = { rho: inp.rockRho, Vs: inp.rockVs, damping: 0.01 };
    const sr = run(() => equivalentLinearResponse({ layers, bedrock, acc, dt, input: 'outcrop' }), 'Site response');
    const pgaS = finite(sr.pgaSurface, 'Surface acceleration');
    // small-strain fundamental frequency from the transfer function relative to the layer base (rigid-base definition)
    const fs = linspace(0.05, 25, 1248), tf = transferFunction({ layers: layers.map((l) => ({ h: l.h, rho: l.rho, Vs: l.Vs, damping: 0.01 })), bedrock, freqs: fs, input: 'within' });
    let f0 = fs[0]; for (let i = 1; i < fs.length - 1; i++) { if (tf.amp[i] > tf.amp[i - 1] && tf.amp[i] >= tf.amp[i + 1]) { f0 = fs[i]; break; } }
    const periods = linspace(0.05, 4, 60);
    const spec = responseSpectrum({ acc: sr.surfaceAcc, dt, periods }), specIn = responseSpectrum({ acc, dt, periods });
    const SaT = responseSpectrum({ acc: sr.surfaceAcc, dt, periods: [inp.Tstr] }).PSA[0];
    const Vbase = SaT * inp.mass * 1000;
    const arias = ariasIntensity({ acc: sr.surfaceAcc, dt });

    // wall: pseudo-static thrusts with kh = ratio × surface PGA / g
    const kh = inp.khRatio * pgaS / G, phi = inp.phiDeg * DEG, gam = inp.gam * 1e3, H = inp.Hw, hw = Math.min(inp.hw + inp.slr, H), W = inp.W * 1e3;
    const KA = rankine(phi).Ka;
    const thrust = (k) => { const mo = mononobeOkabe({ phi, kh: k }); return { KAE: mo.KAE, PAE: 0.5 * gam * H * H * mo.KAE, Pwd: hw > 0 ? westergaard({ kh: k, h: hw, gammaW: 10055 }).force : 0 }; };
    const fsOf = (k) => { const t = thrust(k); return inp.mu * W / (t.PAE + k * W + t.Pwd); };
    const khMax = Math.tan(phi) * 0.98; // Mononobe–Okabe has no solution beyond tan φ for level backfill
    need(inp.mu * W > 0.5 * gam * H * H * KA, 'The wall slides under static earth pressure alone (μW is less than the active thrust): fix the static design first.');
    const khUse = Math.min(kh, khMax);
    const t = thrust(khUse), fsSeis = fsOf(khUse);
    let lo = 0, hi = khMax; if (fsOf(hi) >= 1) lo = hi; else for (let i = 0; i < 50; i++) { const mid = 0.5 * (lo + hi); if (fsOf(mid) >= 1) lo = mid; else hi = mid; }
    const ky = lo;
    // sliding block driven by the full surface motion against the yield acceleration
    const slide = ky * G < pgaS ? newmarkDisplacement({ acc: sr.surfaceAcc, dt, ky: Math.max(ky * G, 1e-3), bothPolarities: true }).displacement : 0;
    const dAllow = inp.dAllowPct / 100 * H;

    // liquefaction of the susceptible layers (submerged column), CSR from the surface PGA
    let z = 0; const pts = [];
    for (const l of col) { z += l.h; pts.push({ z: z - l.h / 2, N160: l.N160, FC: l.FC, gamma: l.rho * G, name: l.name, sus: l.N160 > 0 }); }
    const liq = run(() => sptProfile({ points: pts.map((p) => ({ z: p.z, N160: p.sus ? p.N160 : 60, FC: p.FC, gamma: p.gamma })), waterTableDepth: 0, amax: pgaS, Mw: inp.Mw }), 'Liquefaction');
    const liqRows = liq.map((q, i) => ({ name: pts[i].name, z: r(q.z, 3), CSR: r(q.CSR, 3), CRR: pts[i].sus ? r(q.CRR75, 3) : 0, FS: pts[i].sus ? r(Math.min(q.FS, 99), 3) : 99, sus: pts[i].sus ? 'yes' : 'no' }));
    const susIdx = pts.map((p, i) => (p.sus && p.z <= 20 ? i : -1)).filter((i) => i >= 0);
    const fsLiq = susIdx.length ? Math.min(...susIdx.map((i) => liq[i].FS)) : 99;
    const lpi = susIdx.length ? liquefactionPotentialIndex(susIdx.map((i) => liq[i].z), susIdx.map((i) => liq[i].FS)) : 0;
    const okSlide = slide <= dAllow, okLiq = fsLiq >= 1.25;
    // extra weight needed to bring the yield coefficient up to the design coefficient (per metre run)
    const dW = Math.max(0, (t.PAE + t.Pwd) / Math.max(inp.mu - khUse, 0.02) - W);

    const [pt, pa, ps] = thin([Array.from(sr.t), Array.from(acc), Array.from(sr.surfaceAcc)], 800);
    const advice = [
      okSlide
        ? eng(fsSeis >= 1 ? 'ok' : 'info', `Surface PGA is ${f(pgaS)} m/s² (× ${f(sr.amplification)} on bedrock). With kh = ${f(khUse, 3)} the pseudo-static sliding factor is ${f(fsSeis)}; yield coefficient ky = ${f(ky, 3)} and the Newmark sliding is ${f(slide * 1000)} mm against ${f(dAllow * 1000)} mm allowed.`,
          `The wall meets the displacement criterion${fsSeis < 1 ? ' although the pseudo-static factor is below 1 — the design relies on tolerable sliding' : ''}. Confirm with site-specific spectrum-matched records (at least three, or seven for mean response).`, 'EN 1998-5:2004 §7.3 and Annex E; PIANC WG 34 (2001) performance grades',
          'Accepting a small, repairable sliding displacement instead of a no-movement design avoids a much heavier wall and its concrete.')
        : eng('fail', `Newmark sliding of ${f(slide * 1000)} mm exceeds the ${f(dAllow * 1000)} mm allowed (${f(inp.dAllowPct)} % of height); surface PGA ${f(pgaS)} m/s², yield coefficient ky = ${f(ky, 3)}, pseudo-static factor ${f(fsSeis)} at kh = ${f(khUse, 3)}.`,
          `Raise the yield acceleration: about ${f(dW / 1e3)} kN/m more effective weight, a rougher base or shear key, lighter free-draining rock backfill (higher φ, lower thrust), or ground anchors.`, 'EN 1998-5:2004 §7.3 and Annex E; PIANC WG 34 (2001)',
          `Replacing the backfill wedge with rock fill reduces the seismic thrust without adding concrete; about ${f(0.5 * H * H * Math.tan(Math.PI / 4 - phi / 2) * 1.9, 3)} t of rock per metre run (≈ ${f(0.5 * H * H * Math.tan(Math.PI / 4 - phi / 2) * 1.9 * CO2_ROCK, 2)} t CO₂e, screening factor).`),
      eng(okLiq ? 'ok' : 'fail', okLiq ? `Liquefaction: minimum factor of safety ${f(Math.min(fsLiq, 99))} in the susceptible layers (LPI = ${f(lpi, 2)}).` : `Liquefaction is triggered: factor of safety ${f(fsLiq)} in '${pts[susIdx.find((i) => liq[i].FS === fsLiq)]?.name ?? 'the fill'}' (LPI = ${f(lpi, 2)}; above 15 is severe).`,
        okLiq ? 'No ground improvement is indicated for triggering; confirm with CPT-based assessment and check cyclic softening of the clay.' : 'Densify or drain the liquefiable layer (vibro-compaction, stone columns, deep soil mixing) behind and beneath the wall, or found below it; the wall checks above assume no liquefaction.', 'EN 1998-5:2004 §4.1.4 and Annex B; Youd et al. (2001); Iwasaki et al. (1982) LPI'),
      eng('info', `Site fundamental frequency ${f(f0)} Hz (period ${f(1 / f0)} s); structure period ${f(inp.Tstr)} s gives Sa = ${f(SaT)} m/s² and an elastic base shear of ${f(Vbase / 1e6)} MN. Sea-level rise of ${f(inp.slr)} m is included in the hydrodynamic thrust.`,
        Math.abs(inp.Tstr - 1 / f0) / (1 / f0) < 0.25 ? 'The structure period is close to the site period: expect resonance amplification and consider stiffening or base isolation.' : 'Use the surface spectrum for the structural check, with a behaviour factor appropriate to the structural system.', `ISO 19901-2:2022; EN 1998-1 §3.2; ${SLR_BASIS}`),
      eco(okSlide && okLiq ? 'info' : 'warn', `Design-level shaking leaves ${f(slide * 1000)} mm of permanent sliding; ${okLiq ? 'no liquefaction settlement is expected' : 'liquefaction would add settlement and lateral spreading not quantified here'}.`,
        'Choose the performance grade by consequence: compare the cost of ground improvement and extra wall weight now with the repair and downtime cost of re-aligning the cope, crane rails and services after the design earthquake.', 'PIANC WG 34 (2001) performance-based design; whole-life costing per ISO 15686-5',
        'Ground improvement by vibro-compaction uses little material compared with a heavier wall; avoided post-earthquake reconstruction is the dominant environmental benefit.'),
    ];

    return {
      kpis: [
        { key: 'pgaSurface', label: 'Peak ground acceleration at the surface', value: r(pgaS), unit: 'm/s²', explain: 'Peak of the motion at the top of the soil column after site response. It drives the wall inertia, earth thrust and liquefaction demand.' },
        { key: 'amp', label: 'Amplification bedrock → surface', value: r(sr.amplification), unit: '-', explain: 'Surface over bedrock peak acceleration. Soft sites amplify moderate shaking and can de-amplify strong shaking as the soil softens and damps.' },
        { key: 'f0', label: 'Fundamental frequency of the site', value: r(f0), unit: 'Hz', explain: 'First resonance of the soil column at small strain (about Vs/4H for a uniform layer). Structures with a similar frequency are shaken hardest.' },
        { key: 'SaT', label: 'Spectral acceleration at the structure period', value: r(SaT), unit: 'm/s²', explain: '5 %-damped pseudo-acceleration of the surface motion at the structure\'s period; times mass gives the elastic seismic base shear.' },
        { key: 'Vbase', label: 'Elastic seismic base shear', value: r(Vbase / 1e3), unit: 'kN', explain: 'Spectral acceleration times seismic mass, before any behaviour (ductility) factor. Applied through the foundation into the ground.' },
        { key: 'kh', label: 'Horizontal seismic coefficient', value: r(khUse, 3), unit: '-', explain: 'Pseudo-static coefficient used for the wall thrusts: the stated share of surface PGA/g.' },
        { key: 'KAE', label: 'Mononobe–Okabe active coefficient', value: r(t.KAE), unit: '-', explain: 'Seismic active earth-pressure coefficient. Compare with the static Rankine value to see how much the earthquake adds to the thrust.' },
        { key: 'PAE', label: 'Seismic active thrust', value: r(t.PAE / 1e3), unit: 'kN/m', explain: 'Total (static plus dynamic) earth thrust on the wall during shaking; the dynamic increment acts at about 0.6 H.' },
        { key: 'Pwd', label: 'Hydrodynamic thrust (Westergaard)', value: r(t.Pwd / 1e3), unit: 'kN/m', explain: 'Reduction of seaward water pressure as the wall moves, acting as an additional seaward force; grows with the square of water depth, so with sea-level rise.' },
        { key: 'fsSeis', label: 'Pseudo-static sliding factor', value: r(Math.min(fsSeis, 99)), unit: '-', explain: 'Friction over seismic thrusts and wall inertia at the design coefficient. Below 1 the wall slides during the strongest pulses; the displacement check then governs.' },
        { key: 'ky', label: 'Yield coefficient', value: r(ky, 3), unit: '-', explain: 'Seismic coefficient at which the wall just starts to slide. Shaking above ky·g moves the wall a little with each pulse.' },
        { key: 'slide', label: 'Permanent sliding displacement', value: r(slide * 1000), unit: 'mm', limit: r(dAllow * 1000), sense: 'max', explain: 'Newmark sliding-block displacement accumulated over the surface record. Compared with the serviceability limit for cope alignment and crane rails.' },
        { key: 'fsLiq', label: 'Minimum liquefaction factor of safety', value: r(Math.min(fsLiq, 99)), unit: '-', limit: 1.25, sense: 'min', explain: 'Cyclic resistance over cyclic stress in the weakest susceptible layer within 20 m. Below 1 liquefaction is triggered; EN 1998-5 asks for 1.25.' },
        { key: 'lpi', label: 'Liquefaction potential index', value: r(lpi, 3), unit: '-', explain: 'Depth-weighted severity of liquefaction over the top 20 m: below 5 minor, 5–15 moderate, above 15 severe surface effects.' },
      ],
      checks: [
        check('Sliding displacement', slide * 1000, dAllow * 1000, 'mm', 'PIANC WG 34 (2001) Table 4.1; EN 1998-5 §7.3.2.2', 'Guards against permanent seaward movement that misaligns the cope, crane rails and services. Beyond the limit the wall needs realignment or rebuilding.'),
        sfCheck('Liquefaction triggering', Math.min(fsLiq, 99), 1.25, 'EN 1998-5:2004 §4.1.4(11); Youd et al. (2001)', 'Guards against loss of strength of saturated loose sand. Liquefied backfill multiplies the thrust on the wall and removes foundation support.'),
      ],
      plots: [
        { type: 'line', title: 'Acceleration time histories', x: { label: 'Time', unit: 's' }, y: { label: 'Acceleration', unit: 'm/s²' },
          series: [{ name: 'Surface', x: pt.map((v) => r(v, 3)), y: ps.map((v) => r(v, 3)) }, { name: 'Bedrock outcrop', x: pt.map((v) => r(v, 3)), y: pa.map((v) => r(v, 3)), dash: true }, { name: 'Yield acceleration', x: [pt[0], pt[pt.length - 1]], y: [r(ky * G, 3), r(ky * G, 3)], dash: true }],
          explain: 'Input and surface motion with the wall\'s yield acceleration. Each excursion of the surface trace above the yield line makes the wall slide a little.' },
        { type: 'line', title: 'Response spectra (5 % damping)', x: { label: 'Period', unit: 's' }, y: { label: 'Pseudo-acceleration', unit: 'm/s²' },
          series: [{ name: 'Surface', x: periods.map((v) => r(v, 3)), y: Array.from(spec.PSA, (v) => r(v, 3)) }, { name: 'Bedrock', x: periods.map((v) => r(v, 3)), y: Array.from(specIn.PSA, (v) => r(v, 3)), dash: true }, { name: 'Structure period', x: [inp.Tstr, inp.Tstr], y: [0, r(SaT, 3)] }],
          explain: 'Spectral demand at bedrock and at the surface. The peak near the site period is where the soil column amplifies; a structure period on that peak is the worst case.' },
        { type: 'line', title: 'Site transfer function (small strain)', x: { label: 'Frequency', unit: 'Hz' }, y: { label: 'Surface / base amplitude', unit: '-' },
          series: [{ name: 'Amplification', x: fs.filter((_, i) => i % 4 === 0).map((v) => r(v, 3)), y: tf.amp.filter((_, i) => i % 4 === 0).map((v) => r(Math.min(v, 50), 3)) }],
          explain: 'Frequency-by-frequency amplification of the soil column. Peaks are the site resonances; the first one is the fundamental frequency.' },
        { type: 'line', title: 'Liquefaction factor of safety with depth', x: { label: 'Factor of safety', unit: '-' }, y: { label: 'Depth', unit: 'm' },
          series: [{ name: 'FS', x: liqRows.map((q) => Math.min(q.FS, 3)), y: liqRows.map((q) => -q.z) }, { name: 'FS = 1', x: [1, 1], y: [0, -z], dash: true }],
          explain: 'Factor of safety against liquefaction at the middle of each layer (capped at 3; non-susceptible layers shown at the cap). Points left of the dashed line liquefy.' },
      ],
      tables: [
        { title: 'Strain-compatible layer properties', columns: [{ key: 'name', label: 'Layer' }, { key: 'ratio', label: 'G / Gmax', unit: '-' }, { key: 'damp', label: 'Damping', unit: '%' }, { key: 'gam', label: 'Max shear strain', unit: '%' }],
          rows: sr.layers.map((l, i) => ({ name: col[i].name, ratio: r(l.ratio, 3), damp: r(l.damping * 100, 3), gam: r(l.gammaMax * 100, 3) })), explain: 'Converged modulus reduction, damping and peak strain of each layer. Strains above about 0.5–1 % are beyond the reliable range of equivalent-linear analysis.' },
        { title: 'Liquefaction assessment', columns: [{ key: 'name', label: 'Layer' }, { key: 'z', label: 'Depth', unit: 'm' }, { key: 'sus', label: 'Susceptible' }, { key: 'CSR', label: 'CSR', unit: '-' }, { key: 'CRR', label: 'CRR7.5', unit: '-' }, { key: 'FS', label: 'FS', unit: '-' }],
          rows: liqRows, explain: 'Cyclic stress ratio from the surface PGA against the cyclic resistance from (N1)60 for each layer. Use it to target ground improvement at the layers with FS below 1.25.' },
      ],
      notes: [`${rec && rec.x && rec.x.length >= 64 ? 'Imported record' : 'Built-in synthetic record'} scaled to ${f(inp.pga)} m/s²; site response converged in ${sr.iterations} iteration(s)${sr.converged ? '' : ' (not converged: treat with caution)'}; Arias intensity ${f(arias.arias)} m/s.`,
        `Static Rankine Ka = ${f(KA)}; seismic KAE = ${f(t.KAE)}${kh > khMax ? '; the coefficient was capped at tan φ, the limit of the Mononobe–Okabe solution' : ''}.`],
      advice,
      publish: { 'marine.pgaSurface': pgaS, 'marine.seismicKh': khUse, 'marine.seismicSlide': slide },
    };
  },
  verify: [
    { name: 'Mononobe–Okabe coefficient, φ = 30°, δ = β = θ = 0, kh = 0.2', source: 'Kramer (1996) Eq. 11.19, hand calculation: ψ = atan 0.2 = 11.31°, cos²(φ − ψ) = 0.8973, KAE = 0.8973/{cos²ψ [1 + √(sin 30° sin 18.69°/cos 11.31°)]²} = 0.8973/(0.9615 × 1.9719) = 0.4733. Linear soil over near-rigid rock with the wall coefficient set so that kh = 0.2',
      inputs: { phiDeg: 30, nonlinear: false, pga: 1.96133, khRatio: 1, column: [{ name: 'Rock-like', h: 1, rho: 2400, Vs: 3000, PI: 0, N160: 0, FC: 0 }], rockVs: 3000, rockRho: 2400, record: { x: Array.from({ length: 200 }, (_, i) => i * 0.01), y: Array.from({ length: 200 }, (_, i) => (i >= 50 && i < 100 ? 0.2 * 9.80665 * Math.sin(Math.PI * (i - 50) / 50) : 0)) }, W: 3000, hw: 0, slr: 0 },
      expect: [{ kpi: 'KAE', value: 0.4733, tol: 0.01 }] },
    { name: 'Uniform layer on stiff rock: fundamental frequency f0 = Vs/(4H)', source: 'Kramer (1996) Eq. 7.13. Vs = 200 m/s, H = 20 m → f0 = 2.5 Hz',
      inputs: { nonlinear: false, column: [{ name: 'Uniform', h: 20, rho: 1900, Vs: 200, PI: 0, N160: 0, FC: 0 }], rockVs: 50000 }, expect: [{ kpi: 'f0', value: 2.5, tol: 0.02 }] },
  ],
};
