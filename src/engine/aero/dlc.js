/**
 * Design-load-case driver on the reduced-order coupled model (coupled.js), with the naming of IEC 61400-3-1:2019
 * Table 2 (IEC 61400-1 Table 2) and a deliberately SIMPLIFIED SUBSET of the cases:
 *   1.2 power production, NTM, normal sea state            fatigue   (F)
 *   1.3 power production, ETM, normal sea state            ultimate  (N)
 *   1.4 power production, ECD                              ultimate  (N)
 *   1.6 power production, NTM, severe sea state            ultimate  (N)
 *   2.3 power production, EOG with loss of grid            ultimate  (A)
 *   5.1 emergency stop, NTM                                ultimate  (N)
 *   6.1 parked, 50-year EWM (turbulent) + extreme sea      ultimate  (N)
 *   6.2 parked, 50-year EWM, grid loss (large yaw errors)  ultimate  (A)
 *   7.1 parked with fault, 1-year EWM (optional)           ultimate  (A)
 * Simplifications (in addition to those of coupled.js): rotor-effective wind only, so direction changes enter as
 * the cosine of the misalignment; one wave direction aligned with the wind; no wind/wave misalignment, currents or
 * water-level variation unless given in the sea state; few seeds and short records unless the caller asks for more;
 * no statistical extrapolation (DLC 1.1) and no fault cases other than grid loss. Results are for concept screening
 * of support-structure loads, not for certification.
 */

import { simulateCoupled, rainflowCount, seriesStats } from './coupled.js';
import { operatingSchedule } from './control.js';
import { ntmSigma, etmSigma, ewm, eogGust, eogSpeed, ecd, turbulenceScale, rotorEffectiveTurbulence } from './turbulence.js';

function need(cond, name, msg) {
  if (!cond) throw new RangeError(`${name}: ${msg}`);
}

const D2R = Math.PI / 180;

/**
 * Partial safety factors for loads (IEC 61400-1 Table 3): normal (N) 1.35, abnormal (A) 1.1, transport and erection
 * (T) 1.5; fatigue (F) 1.0 (clause 7.6.3). Applied to characteristic extremes.
 */
export const LOAD_SAFETY_FACTORS = { N: 1.35, A: 1.1, T: 1.5, F: 1.0 };

/**
 * GENERIC wind-sea correlation for when no site scatter table is available: fetch-limited deep-water growth curves
 * of the Shore Protection Manual (US Army CERC 1984, Eqs. 3-33, 3-34),
 *   g Hs / U^2 = 1.6e-3 (g F / U^2)^(1/2),    g Tp / U = 0.2857 (g F / U^2)^(1/3),
 * capped by the fully-developed limits g Hs/U^2 = 0.243 and g Tp/U = 8.134. U is taken as the 10-m wind speed
 * derived from the hub-height speed by the power law (the manual's adjusted wind speed is not applied). Wind sea
 * only: no swell. Replace with the site's conditional Hs-Tp table for design.
 * @param {object} o
 * @param {number} o.Vhub hub-height wind speed [m/s] @param {number} o.hubHeight [m]
 * @param {number} [o.fetch=100000] fetch [m] @param {number} [o.shearExp=0.14] @param {number} [o.g=9.80665]
 * @returns {{Hs:number, Tp:number, U10:number, generic:boolean}} [m], [s], [m/s]
 */
export function genericSeaState(o) {
  need(o && o.Vhub >= 0 && o.hubHeight > 0, 'Vhub/hubHeight', 'must be valid');
  const g = o.g ?? 9.80665;
  const U = Math.max(0.1, o.Vhub * Math.pow(10 / o.hubHeight, o.shearExp ?? 0.14));
  const x = g * (o.fetch ?? 1e5) / (U * U);
  return { Hs: Math.min(0.243, 1.6e-3 * Math.sqrt(x)) * U * U / g, Tp: Math.min(8.134, 0.2857 * Math.cbrt(x)) * U / g, U10: U, generic: true };
}

function seaStateAt(site, U) {
  const ss = site.seaState;
  if (typeof ss === 'function') return ss(U);
  if (Array.isArray(ss) && ss.length) {
    const tab = ss.slice().sort((a, b) => a.U - b.U);
    if (U <= tab[0].U) return { Hs: tab[0].Hs, Tp: tab[0].Tp };
    if (U >= tab[tab.length - 1].U) return { Hs: tab[tab.length - 1].Hs, Tp: tab[tab.length - 1].Tp };
    let i = 0;
    while (tab[i + 1].U < U) i++;
    const w = (U - tab[i].U) / (tab[i + 1].U - tab[i].U);
    return { Hs: tab[i].Hs + w * (tab[i + 1].Hs - tab[i].Hs), Tp: tab[i].Tp + w * (tab[i + 1].Tp - tab[i].Tp) };
  }
  return genericSeaState({ Vhub: U, hubHeight: site.hubHeight, fetch: site.fetch, shearExp: site.shearExp });
}

/**
 * Wind-speed bins with Weibull probabilities and the associated sea state (joint wind-wave description for
 * fatigue): P(bin) = exp(-(lo/A)^k) - exp(-(hi/A)^k); Hs and Tp conditional on the bin-centre wind speed from the
 * user table `seaState` ([{ U, Hs, Tp }] or function) or else from genericSeaState (flagged generic).
 * @param {object} o
 * @param {number} o.weibullA scale [m/s] @param {number} o.weibullK shape [-]
 * @param {number} o.vIn @param {number} o.vOut cut-in and cut-out [m/s] @param {number} [o.binWidth=2] [m/s]
 * @param {number} o.hubHeight [m]
 * @param {Array<{U:number,Hs:number,Tp:number}>|function(number):{Hs:number,Tp:number}} [o.seaState]
 * @param {number} [o.fetch] @param {number} [o.shearExp]
 * @returns {Array<{U:number, lo:number, hi:number, prob:number, Hs:number, Tp:number, generic:boolean}>} bins; the
 *   probabilities sum to the fraction of time between cut-in and cut-out
 */
export function windWaveBins(o) {
  need(o && o.weibullA > 0 && o.weibullK > 0, 'weibullA/weibullK', 'must be > 0');
  need(o.vOut > o.vIn && o.vIn >= 0, 'vIn/vOut', 'vOut must exceed vIn >= 0');
  const bw = o.binWidth ?? 2;
  const cdf = (v) => 1 - Math.exp(-Math.pow(v / o.weibullA, o.weibullK));
  const bins = [];
  const nb = Math.max(1, Math.round((o.vOut - o.vIn) / bw));
  for (let i = 0; i < nb; i++) {
    const lo = o.vIn + (o.vOut - o.vIn) * i / nb; const hi = o.vIn + (o.vOut - o.vIn) * (i + 1) / nb;
    const U = 0.5 * (lo + hi);
    const ss = seaStateAt(o, U);
    bins.push({ U, lo, hi, prob: cdf(hi) - cdf(lo), Hs: ss.Hs, Tp: ss.Tp, generic: !!ss.generic });
  }
  return bins;
}

function steadyStart(model, U) {
  const cfg = model.controller;
  const N = cfg.gearboxRatio ?? 1;
  const s = operatingSchedule({ surface: model.surface, ratedPower: cfg.ratedPowerMech * (cfg.efficiency ?? 1), efficiency: cfg.efficiency ?? 1,
    ratedRotorSpeed: cfg.ratedGenSpeed / N, minRotorSpeed: cfg.minGenSpeed / N, gearboxRatio: N, rho: model.rho ?? 1.225,
    pitchFineDeg: cfg.pitchMin / D2R, speeds: [U], cutIn: 0, cutOut: 1e3 });
  return { omega: s.omega[0], pitchDeg: s.pitchDeg[0] };
}

function extremes(r, tT, dt) {
  const i0 = Math.round(tT / dt);
  const mm = seriesStats(r.mudlineMoment, i0); const tm = seriesStats(r.towerBaseMoment, i0);
  const ms = seriesStats(r.mudlineShear, i0); const dd = seriesStats(r.towerTopDisp, i0);
  return { mudlineMoment: mm.absMax, towerBaseMoment: tm.absMax, mudlineShear: ms.absMax, towerTopDisp: dd.absMax, maxRotorSpeed: seriesStats(r.rotorSpeed, i0).max };
}

function summarise(name, description, type, runs) {
  const psf = LOAD_SAFETY_FACTORS[type];
  const keys = ['mudlineMoment', 'towerBaseMoment', 'mudlineShear', 'towerTopDisp', 'maxRotorSpeed'];
  const char = {}; const governing = {};
  for (const k of keys) {
    char[k] = 0;
    for (const r of runs) if (r.extremes[k] > char[k]) { char[k] = r.extremes[k]; governing[k] = r.label; }
  }
  const design = {};
  for (const k of ['mudlineMoment', 'towerBaseMoment', 'mudlineShear']) design[k] = psf * char[k];
  return { dlc: name, description, type, psf, characteristic: char, design, governing, runs };
}

/**
 * Run a set of design load cases and tabulate ultimate and fatigue loads at the mudline and the tower base.
 * Power-production cases start from the steady operating point of the mean wind speed and discard `tTransient`.
 * @param {object} o
 * @param {object} o.model simulateCoupled arguments common to all cases: surface, turbine, controller, structure or
 *   support, rho, dt, actuator, ... (wind, waves, duration and events are set per case)
 * @param {object} o.site { hubHeight [m], rotorDiameter [m], waterDepth [m], Vref [m/s], Iref [-],
 *   Vave [m/s] (default 0.2 Vref), ratedWind, cutIn, cutOut [m/s], weibullA [m/s], weibullK [-],
 *   seaState (normal sea state: table [{ U, Hs, Tp }] or function; default generic), Hs50, Tp50, Hs1, Tp1
 *   (extreme sea states) [m, s], gamma (JONSWAP, default 3.3), fetch [m], shearExp }
 * @param {string[]} [o.cases=['1.2','1.3','1.4','1.6','2.3','5.1','6.1','6.2']] cases to run ('7.1' optional)
 * @param {number} [o.nSeeds=1] turbulence/wave seeds per condition (IEC asks for at least six 10-min seeds)
 * @param {number} [o.duration=600] simulated length of stochastic cases after the transient [s]
 * @param {number} [o.tTransient=30] discarded start [s]
 * @param {number} [o.seed=1] base seed
 * @param {number[]} [o.windSpeeds] mean wind speeds for the ultimate power-production cases [m/s];
 *   default rated-2, rated, rated+2, cutOut
 * @param {number} [o.binWidth=2] wind bin width for DLC 1.2 [m/s]
 * @param {number} [o.lifetimeYears=25] @param {number} [o.availability=1] fraction of time in power production
 * @param {number} [o.wohlerM=4] S-N slope for DELs @param {number} [o.nEq=1e7] equivalent cycle number
 * @param {function(object):number} [o.snDamage] hook ({ channel, ranges, counts }) -> Miner damage for lifetime-scaled
 *   cycle counts (ranges in N m of section moment: convert to stress inside the hook)
 * @param {function(Float64Array):object} [o.rainflow] cycle-counting hook; default coupled.js rainflowCount
 * @param {number|function(number):number} [o.parkedCt] thrust coefficient of the parked/feathered rotor referred to
 *   the swept area, as a number or function of yaw error [deg] (e.g. from bemSolve at omega = 0, pitch 90 deg);
 *   required for DLC 6.x/7.1
 * @param {number|function(number):number} [o.parkedFaultCt] as parkedCt for the fault condition of DLC 7.1
 * @param {number[]} [o.yawErrors61Deg=[-8, 0, 8]] @param {number[]} [o.yawErrors62Deg=[-180,-90,-30,30,90,180]]
 * @param {number[]} [o.gridLossTimes=[0, 2.45, 4, 5.25]] grid-loss instants after the start of the EOG [s]
 * @returns {{cases:object[], table:object[], fatigue:(object|null), psf:object, notes:string[]}} per-case summaries
 *   (characteristic and design extremes, governing run), a flat table (moments N m, shear N), DLC 1.2 fatigue
 *   results ({ bins, del: { mudlineMoment, towerBaseMoment } [N m], damage }) and remarks
 */
export function runDesignLoadCases(o) {
  need(o && o.model && o.model.surface && o.model.controller, 'model', 'surface and controller required');
  const site = o.site;
  need(site && site.hubHeight > 0 && site.rotorDiameter > 0, 'site', 'hubHeight and rotorDiameter required');
  need(site.Vref > 0 && site.Iref > 0, 'site.Vref/Iref', 'must be > 0');
  need(site.ratedWind > 0 && site.cutOut > site.cutIn, 'site.ratedWind/cutIn/cutOut', 'required');
  const model = o.model;
  const dt = model.dt;
  need(dt > 0, 'model.dt', 'must be > 0');
  const cases = o.cases ?? ['1.2', '1.3', '1.4', '1.6', '2.3', '5.1', '6.1', '6.2'];
  const nSeeds = o.nSeeds ?? 1;
  const dur = o.duration ?? 600;
  const tT = o.tTransient ?? 30;
  const seed0 = o.seed ?? 1;
  const D = site.rotorDiameter; const R = 0.5 * D;
  const lam1 = turbulenceScale(site.hubHeight);
  const Vave = site.Vave ?? 0.2 * site.Vref;
  const speeds = o.windSpeeds ?? [site.ratedWind - 2, site.ratedWind, site.ratedWind + 2, site.cutOut];
  const depth = site.waterDepth;
  const hasWaves = depth > 0 && (model.structure || model.support);
  const waveFor = (Hs, Tp, seed) => (hasWaves && Hs > 0 ? { Hs, Tp, gamma: site.gamma ?? 3.3, depth, seed, macCamyFuchs: !!site.macCamyFuchs } : null);
  const turb = (V, sigma, seed, length) => rotorEffectiveTurbulence({ Vhub: V, sigma1: sigma, hubHeight: site.hubHeight, R, duration: length, dt, seed }).u;
  const run = (args) => simulateCoupled({ ...model, ...args });
  const notes = [];
  const out = [];
  let fatigue = null;
  const ctOf = (c, yaw) => (typeof c === 'function' ? c(yaw) : c);

  for (const id of cases) {
    if (id === '1.2') {
      need(site.weibullA > 0 && site.weibullK > 0, 'site.weibullA/weibullK', 'required for DLC 1.2');
      const bins = windWaveBins({ weibullA: site.weibullA, weibullK: site.weibullK, vIn: site.cutIn, vOut: site.cutOut, binWidth: o.binWidth ?? 2,
        hubHeight: site.hubHeight, seaState: site.seaState, fetch: site.fetch, shearExp: site.shearExp });
      if (bins.some((b) => b.generic)) notes.push('DLC 1.2 uses the GENERIC fetch-limited wind-sea correlation: replace with the site scatter table.');
      const m = o.wohlerM ?? 4;
      const life = (o.lifetimeYears ?? 25) * 365.25 * 86400 * (o.availability ?? 1);
      const chans = ['mudlineMoment', 'towerBaseMoment'];
      const sum = { mudlineMoment: 0, towerBaseMoment: 0 };
      const all = { mudlineMoment: { ranges: [], counts: [] }, towerBaseMoment: { ranges: [], counts: [] } };
      const rf = o.rainflow || rainflowCount;
      const i0 = Math.round(tT / dt);
      bins.forEach((b, ib) => {
        b.del = { mudlineMoment: 0, towerBaseMoment: 0 };
        const acc = { mudlineMoment: 0, towerBaseMoment: 0 };
        for (let s = 0; s < nSeeds; s++) {
          const seed = seed0 + 1000 * ib + s;
          const r = run({ wind: turb(b.U, ntmSigma({ Vhub: b.U, Iref: site.Iref }), seed, dur + tT), waves: waveFor(b.Hs, b.Tp, seed + 500),
            duration: dur + tT, tTransient: tT, initial: steadyStart(model, b.U) });
          for (const ch of chans) {
            const c = rf(r[ch].subarray(i0));
            const scale = life * b.prob / (dur * nSeeds);
            for (let i = 0; i < c.ranges.length; i++) {
              acc[ch] += c.counts[i] * Math.pow(c.ranges[i], m) / (dur * nSeeds);
              if (o.snDamage) { all[ch].ranges.push(c.ranges[i]); all[ch].counts.push(c.counts[i] * scale); }
            }
          }
        }
        for (const ch of chans) { b.del[ch] = Math.pow(acc[ch] * 1, 1 / m); sum[ch] += acc[ch] * life * b.prob; }
      });
      const nEq = o.nEq ?? 1e7;
      fatigue = {
        bins, wohlerM: m, nEq, lifetimeSeconds: life,
        del: { mudlineMoment: Math.pow(sum.mudlineMoment / nEq, 1 / m), towerBaseMoment: Math.pow(sum.towerBaseMoment / nEq, 1 / m) },
        damage: o.snDamage ? { mudlineMoment: o.snDamage({ channel: 'mudlineMoment', ...all.mudlineMoment }), towerBaseMoment: o.snDamage({ channel: 'towerBaseMoment', ...all.towerBaseMoment }) } : null,
        psf: LOAD_SAFETY_FACTORS.F,
      };
      continue;
    }
    const runs = [];
    if (id === '1.3' || id === '1.6' || id === '5.1') {
      for (const V of speeds) {
        for (let s = 0; s < nSeeds; s++) {
          const seed = seed0 + Math.round(V * 10) + 7919 * s + (id === '1.3' ? 11 : id === '1.6' ? 22 : 33);
          const sigma = id === '1.3' ? etmSigma({ Vhub: V, Iref: site.Iref, Vave }) : ntmSigma({ Vhub: V, Iref: site.Iref });
          let ss = seaStateAt(site, V);
          if (id === '1.6') {
            ss = site.Hs50 > 0 ? { Hs: site.Hs50, Tp: site.Tp50 } : ss;
            if (!(site.Hs50 > 0)) notes.push('DLC 1.6: no severe sea state given (site.Hs50/Tp50); the normal sea state was used.');
          }
          const tStop = tT + 0.5 * dur;
          const r = run({ wind: turb(V, sigma, seed, dur + tT), waves: waveFor(ss.Hs, ss.Tp, seed + 500), duration: dur + tT, tTransient: tT,
            initial: steadyStart(model, V), events: id === '5.1' ? (t) => (t >= tStop ? { command: 'estop' } : null) : undefined });
          runs.push({ label: `V = ${V.toFixed(1)} m/s, seed ${s + 1}, Hs = ${ss.Hs.toFixed(2)} m`, V, extremes: extremes(r, tT, dt) });
        }
      }
      out.push(summarise(id, id === '1.3' ? 'Power production, ETM, NSS' : id === '1.6' ? 'Power production, NTM, SSS' : 'Emergency stop, NTM, NSS', 'N', runs));
    } else if (id === '1.4') {
      for (const V of [site.ratedWind - 2, site.ratedWind, site.ratedWind + 2]) {
        const g = ecd({ Vhub: V });
        const t0 = tT + 10;
        for (const sign of [1, -1]) {
          const ss = seaStateAt(site, V);
          const r = run({ wind: (t) => g.speed(t - t0) * Math.cos(sign * g.directionDeg(t - t0) * D2R), waves: waveFor(ss.Hs, ss.Tp, seed0 + 41),
            duration: t0 + 60, tTransient: tT, initial: steadyStart(model, V) });
          runs.push({ label: `V = ${V.toFixed(1)} m/s, direction change ${sign > 0 ? '+' : '-'}${g.thetaCgDeg.toFixed(0)} deg`, V, extremes: extremes(r, tT, dt) });
        }
      }
      out.push(summarise('1.4', 'Power production, ECD, NSS', 'N', runs));
    } else if (id === '2.3') {
      for (const V of [site.ratedWind - 2, site.ratedWind + 2, site.cutOut]) {
        const sigma = ntmSigma({ Vhub: V, Iref: site.Iref });
        const ve1 = ewm({ Vref: site.Vref, z: site.hubHeight, hubHeight: site.hubHeight }).Ve1;
        const vg = eogGust({ Vhub: V, sigma1: sigma, D, lambda1: lam1, Ve1: ve1 });
        const t0 = tT + 10;
        for (const tl of (o.gridLossTimes ?? [0, 2.45, 4, 5.25])) {
          const ss = seaStateAt(site, V);
          const r = run({ wind: (t) => eogSpeed(t - t0, { Vhub: V, Vgust: vg }), waves: waveFor(ss.Hs, ss.Tp, seed0 + 43),
            duration: t0 + 60, tTransient: tT, initial: steadyStart(model, V), events: (t) => (t >= t0 + tl ? { gridLoss: true } : null) });
          runs.push({ label: `V = ${V.toFixed(1)} m/s, Vgust = ${vg.toFixed(2)} m/s, grid loss at +${tl} s`, V, extremes: extremes(r, tT, dt) });
        }
      }
      out.push(summarise('2.3', 'Power production, EOG with loss of electrical network', 'A', runs));
    } else if (id === '6.1' || id === '6.2' || id === '7.1') {
      const ctSrc = id === '7.1' ? o.parkedFaultCt : o.parkedCt;
      need(ctSrc !== undefined, id === '7.1' ? 'parkedFaultCt' : 'parkedCt', `required for DLC ${id}`);
      const e = ewm({ Vref: site.Vref, z: site.hubHeight, hubHeight: site.hubHeight });
      const V = id === '7.1' ? e.V1 : e.V50;
      const Hs = id === '7.1' ? site.Hs1 : site.Hs50; const Tp = id === '7.1' ? site.Tp1 : site.Tp50;
      if (!(Hs > 0)) notes.push(`DLC ${id}: no extreme sea state given; run without waves.`);
      const yaws = id === '6.1' ? (o.yawErrors61Deg ?? [-8, 0, 8]) : id === '6.2' ? (o.yawErrors62Deg ?? [-180, -90, -30, 30, 90, 180]) : [0];
      for (const yaw of yaws) {
        for (let s = 0; s < nSeeds; s++) {
          const seed = seed0 + 6100 + 31 * s + Math.round(yaw);
          const r = run({ controller: null, wind: turb(V, 0.11 * V, seed, dur + tT), waves: Hs > 0 ? waveFor(Hs, Tp, seed + 500) : null,
            duration: dur + tT, tTransient: tT, initial: { omega: 0, pitchDeg: 90 }, aeroOverride: { ct: ctOf(ctSrc, yaw), cq: 0 }, towerWindDrag: true,
            shearExp: 0.11 });
          runs.push({ label: `V = ${V.toFixed(1)} m/s, yaw error ${yaw} deg, seed ${s + 1}`, V, extremes: extremes(r, tT, dt) });
        }
      }
      out.push(summarise(id, id === '6.1' ? 'Parked, EWM 50-year (turbulent), ESS' : id === '6.2' ? 'Parked, EWM 50-year, grid loss (yaw misalignment), ESS' : 'Parked with fault, EWM 1-year', id === '6.1' ? 'N' : 'A', runs));
    } else {
      throw new RangeError(`cases: unknown design load case "${id}"`);
    }
  }
  const table = out.map((c) => ({
    dlc: c.dlc, description: c.description, type: c.type, psf: c.psf,
    mudlineMomentChar: c.characteristic.mudlineMoment, mudlineMomentDesign: c.design.mudlineMoment,
    towerBaseMomentChar: c.characteristic.towerBaseMoment, towerBaseMomentDesign: c.design.towerBaseMoment,
    mudlineShearChar: c.characteristic.mudlineShear, mudlineShearDesign: c.design.mudlineShear,
    governing: c.governing.mudlineMoment,
  }));
  if (fatigue) {
    table.push({ dlc: '1.2', description: 'Power production, NTM, NSS (fatigue)', type: 'F', psf: 1.0,
      mudlineMomentDel: fatigue.del.mudlineMoment, towerBaseMomentDel: fatigue.del.towerBaseMoment, wohlerM: fatigue.wohlerM, nEq: fatigue.nEq });
  }
  if (nSeeds < 6) notes.push(`Only ${nSeeds} seed(s) per condition were run; IEC 61400-1 requires at least six 10-min realisations (or one continuous hour) for stochastic cases.`);
  return { cases: out, table, fatigue, psf: LOAD_SAFETY_FACTORS, notes: [...new Set(notes)] };
}
