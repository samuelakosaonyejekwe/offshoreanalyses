// Site characterisation: pull every live source for a coordinate and reduce it to the design values the suites
// bind to (docs/CONTRACT.md §3). Each derived number keeps its provenance and method so it can be audited.
import * as C from './connectors.js';
import { fitGumbel, returnValue, scatterDiagram } from '../engine/reliability/extreme.js';
import { getProject, publish, save, emit } from '../core/store.js';
import { SHELF_DEPTH, SEABED_CLASSES, classifySeabed, gradeSeabedSamples, reduceSeabedSamples, shomClass } from './seabed.js';
import { hydroSorting, combineSeabedEvidence, seabedConfidence } from './seabedevidence.js';
import { SCREENING_KEYS } from './soilscreening.js';

const clean = (a) => a.filter((v) => typeof v === 'number' && Number.isFinite(v));
const mean = (a) => a.reduce((s, v) => s + v, 0) / (a.length || 1);
// reduce-based extremes: spreading a 175 000-element hindcast into Math.max overflows the call stack
const maxOf = (a) => a.reduce((m, v) => (v > m ? v : m), -Infinity);
const minOf = (a) => a.reduce((m, v) => (v < m ? v : m), Infinity);
const quant = (a, p) => { const s = clean(a).sort((x, y) => x - y); return s.length ? s[Math.min(s.length - 1, Math.floor(p * s.length))] : NaN; };

/** Annual-maximum Gumbel extremes; falls back to the empirical maximum when the record is too short to fit. */
function extremes(times, values, label) {
  const pairs = times.map((t, i) => [t, values[i]]).filter(([, v]) => typeof v === 'number' && Number.isFinite(v));
  if (pairs.length < 24 * 300) return null;
  // Calendar year is read from the ISO string itself (no time-zone shift), and years with less than 80 % of
  // their hours are dropped: a stray first or last partial year would otherwise enter as a false low maximum.
  const byYear = new Map();
  for (const [t, v] of pairs) { const y = Number(String(t).slice(0, 4)); const e = byYear.get(y) || { n: 0, max: -Infinity }; e.n++; if (v > e.max) e.max = v; byYear.set(y, e); }
  const kept = [...byYear.entries()].filter(([, e]) => e.n >= 0.8 * 8760).sort((a, b) => a[0] - b[0]);
  const am = { years: kept.map(([y]) => y), maxima: kept.map(([, e]) => e.max) };
  const maxima = am.maxima;
  // fewer than five complete years cannot support a fit: report the observed maximum with a stated allowance
  if (maxima.length < 5) return { n: maxima.length, v1: quant(maxima, 0.5), v100: 1.25 * maxOf(maxima), method: `${label}: only ${maxima.length} complete years — 100-year value taken as 1.25 × the observed maximum, not an extreme-value fit` };
  // L-moments are the more robust estimator for the short samples that annual maxima give
  const fit = fitGumbel(maxima, { method: maxima.length >= 15 ? 'mle' : 'lmom' });
  const g = { dist: 'gumbel', loc: fit.loc, scale: fit.scale };
  // 1-year value: the Gumbel mode of the annual maxima (return period 1/(1 − e⁻¹) ≈ 1.58 years)
  return { n: maxima.length, v1: returnValue(g, 1.582), v10: returnValue(g, 10), v50: returnValue(g, 50), v100: returnValue(g, 100), fit: g, maxima, years: am.years,
    method: `${label}: Gumbel (${maxima.length >= 15 ? 'maximum likelihood' : 'L-moments'}) on ${maxima.length} annual maxima${maxima.length < 15 ? ' — short record, so the 100-year value carries wide uncertainty' : ''}` };
}

/** Maximum-likelihood Weibull fit for wind speed (Newton iteration on the shape equation). */
export function weibullFit(v) {
  const x = clean(v).filter((u) => u > 0);
  if (x.length < 100) return null;
  const lnx = x.map(Math.log), ml = mean(lnx);
  let k = 2;
  for (let it = 0; it < 60; it++) {
    let a = 0, b = 0, c = 0;
    for (let i = 0; i < x.length; i++) { const p = x[i] ** k; a += p * lnx[i]; b += p; c += p * lnx[i] * lnx[i]; }
    const f = a / b - 1 / k - ml, df = (c * b - a * a) / (b * b) + 1 / (k * k);
    const kn = k - f / df;
    if (!Number.isFinite(kn) || kn <= 0) break;
    if (Math.abs(kn - k) < 1e-7) { k = kn; break; }
    k = kn;
  }
  const A = (x.reduce((s, u) => s + u ** k, 0) / x.length) ** (1 / k);
  return { A, k, mean: mean(x) };
}

// Screening-level peak ground acceleration: each catalogued event is attenuated to the site with
// Campbell (1981) PGA[g] = 0.0159 e^(0.868 M) (R + 0.0606 e^(0.7 M))^(−1.09), R in km, and the annual exceedance
// rate of the resulting site PGA is extrapolated log-linearly to 1/475 per year.
function seismicScreen(fc, lat, lon) {
  const feats = fc?.features || [];
  if (!feats.length) return { n: 0, pga475: 0, note: 'No M ≥ 4.5 events within 300 km since 1970: low historical seismicity.' };
  const hav = (la1, lo1, la2, lo2) => { const r = Math.PI / 180, a = Math.sin((la2 - la1) * r / 2) ** 2 + Math.cos(la1 * r) * Math.cos(la2 * r) * Math.sin((lo2 - lo1) * r / 2) ** 2; return 12742 * Math.asin(Math.sqrt(a)); };
  const tMin = Math.min(...feats.map((f) => f.properties.time));
  const yrs = Math.max(1, (Date.now() - tMin) / (365.25 * 86400e3));
  const pga = feats.map((f) => {
    const [lo, la, dep] = f.geometry.coordinates;
    const R = Math.hypot(hav(lat, lon, la, lo), dep || 10), M = f.properties.mag;
    return 0.0159 * Math.exp(0.868 * M) * (R + 0.0606 * Math.exp(0.7 * M)) ** -1.09;
  }).sort((a, b) => b - a);
  const mags = feats.map((f) => f.properties.mag);
  // exceedance rate λ(a) from ranked site accelerations; fit ln λ = c0 + c1 ln a over the upper half
  const pts = pga.slice(0, Math.max(5, Math.ceil(pga.length / 2))).map((a, i) => [Math.log(a), Math.log((i + 1) / yrs)]);
  const mx = mean(pts.map((p) => p[0])), my = mean(pts.map((p) => p[1]));
  const c1 = pts.reduce((s, p) => s + (p[0] - mx) * (p[1] - my), 0) / (pts.reduce((s, p) => s + (p[0] - mx) ** 2, 0) || 1e-300);
  const c0 = my - c1 * mx;
  // a sparse catalogue cannot support an extrapolation: fall back to 1.5 × the largest historical site value,
  // and never let the fitted value run beyond 5 × that maximum
  let pga475 = feats.length >= 10 && c1 < 0 ? Math.exp((Math.log(1 / 475) - c0) / c1) : 1.5 * pga[0];
  pga475 = Math.min(pga475, 5 * pga[0], 1.5);
  return { n: feats.length, years: yrs, maxMag: maxOf(mags), maxSitePga: pga[0], pga475, note: 'Screening estimate from the historical catalogue with a single attenuation relation (with fewer than 10 events: 1.5 × the largest historical site acceleration) — not a probabilistic seismic hazard assessment. Use a site-specific PSHA or the national hazard map for design.' };
}

export { SHELF_DEPTH, SEABED_CLASSES, classifySeabed, classifyGrainSize, gradeSeabedSamples, reduceSeabedSamples } from './seabed.js';

/**
 * Where the open sources hold no survey map or sample: the class of the worldwide sediment map (read at the site,
 * or from the bundled grid), the sediment type of the deep-sea grid, and — on a shelf — the comparison with the
 * hydrodynamic setting that decides whether a grain size may be passed on. Pure, so it can be tested.
 * @param {{ depth: number, live?: { code: string, sheet: number|null }|null, grid?: { sediment: object|null, deepsea: object|null }|null, hydro?: object|null }} o
 * @returns {null|object} the `data.seabed` record without provenance; null when no source has a value
 */
export function worldSeabed({ depth, live = null, grid = null, hydro = null }) {
  const node = grid?.sediment || null, deep = depth >= 500 ? grid?.deepsea || null : null;
  const sc = shomClass(live?.code) || (node && shomClass(node.code));
  if (!sc && !deep) return null;
  const fromLive = !!shomClass(live?.code), shelf = depth < SHELF_DEPTH;
  // in deep water the composition (carbonate, opal) matters more to engineering than the grain-size name of a
  // fine sediment, so the deep-sea type names the class wherever the map shows mud, silt or clay
  const fine = !sc || ['mud', 'sandyMud', 'silt', 'clay'].includes(sc.cls);
  const cls = deep && fine && deep.cls !== 'lithogenous' ? deep.cls : sc ? sc.cls : deep.cls;
  const c = SEABED_CLASSES[cls];
  const sb = { cls, ...c, basis: sc ? (fromLive ? 'worldmap' : 'grid') : 'deepsea', described: sc?.name || null, definition: sc?.definition || null, sheet: fromLive ? live.sheet : null, composition: deep?.name || null, shelf, hydro,
    resolutionKm: !fromLive && node ? node.stepDeg * 111.2 : null, offsetKm: !fromLive && node?.shifted ? node.offsetKm : 0, d50Class: c.d50 };
  if (shelf) {
    // on a shelf sediment changes over a few kilometres: the mapped class needs the hydrodynamics to agree
    sb.evidence = combineSeabedEvidence({ cls, d50: c.d50, d50Range: c.d50Range, hydro });
    sb.d50 = sb.evidence.verdict === 'agree' ? sb.evidence.d50 : null;
    if (sb.evidence.range) sb.d50Range = sb.evidence.range;
  }
  return sb;
}

/**
 * Pull and reduce everything for the project site. Every step is independent: one source failing leaves the
 * rest intact and is listed in `issues`.
 */
export async function pullSite(lat, lon, { years = 20, onProgress, force = false } = {}) {
  const data = { pulledAt: new Date().toISOString(), issues: [], provenance: {} };
  const bus = { 'site.lat': lat, 'site.lon': lon, 'site.rhoW': 1025 };
  const step = async (label, frac, fn) => {
    onProgress?.(frac, label);
    try { await fn(); } catch (e) { data.issues.push(`${label}: ${e.message}`); }
  };

  await step('Water depth', 0.03, async () => {
    const r = await C.depthAt(lat, lon, { force });
    data.elevation = r.data.elevation;
    data.provenance.depth = { source: 'NOAA NCEI global DEM mosaic', at: r.at };
    if (r.data.elevation < 0) bus['site.depth'] = -r.data.elevation;
    else data.issues.push(`Location is on land (elevation ${r.data.elevation.toFixed(0)} m): marine design values are taken from the nearest sea cell.`);
  });
  await step('Place', 0.06, async () => {
    const r = await C.reversePlace(lat, lon);
    data.place = { name: r.data.display_name || r.data.name || '', country: r.data.address?.country || '', iso2: (r.data.address?.country_code || '').toUpperCase() };
  });
  await step('Sea state now and 8-day forecast', 0.1, async () => {
    const r = await C.marineNow(lat, lon, { force });
    data.now = { ...r.data.current, units: r.data.current_units };
    data.forecast = r.data.hourly;
    data.provenance.now = { source: 'Open-Meteo Marine', at: r.at };
    if (typeof r.data.current?.sea_surface_temperature === 'number') bus['site.Tsea'] = r.data.current.sea_surface_temperature;
  });
  await step('Wave hindcast', 0.15, async () => {
    const r = await C.marineHindcast(lat, lon, { years, force, onProgress: (f, m) => onProgress?.(0.15 + 0.3 * f, m) });
    const hs = r.data.wave_height, t = r.data.time;
    // the long hindcast often carries the mean period only: the peak period is then taken as 1.15 × mean
    // (between 1.11 for a JONSWAP sea with γ = 3.3 and 1.17 for a Pierson–Moskowitz sea, energy period basis)
    const tpRaw = r.data.wave_peak_period || [], tm = r.data.wave_period || [];
    const nPeak = tpRaw.filter((p) => typeof p === 'number' && p > 0).length, nMean = tm.filter((p) => typeof p === 'number' && p > 0).length;
    const nHs = hs.filter((v) => typeof v === 'number').length;   // the record may start years after the requested period
    const tpFrom = nPeak >= 0.5 * nHs ? 'peak' : nMean >= 0.5 * nHs ? 'mean' : 'none';
    const tp = tpFrom === 'peak' ? tpRaw : tpFrom === 'mean' ? tm.map((p) => (typeof p === 'number' && p > 0 ? 1.15 * p : null)) : hs.map(() => null);
    const ex = extremes(t, hs, 'Hs');
    if (!ex) throw new Error('No wave data at this location (too far inland?)');
    // associated peak period: Tp = a·√Hs fitted to the storms above the 99th percentile
    const thr = quant(hs, 0.99);
    const storm = hs.map((h, i) => [h, tp[i]]).filter(([h, p]) => h >= thr && typeof p === 'number' && p > 0);
    const a = storm.length ? mean(storm.map(([h, p]) => p / Math.sqrt(h))) : 4.4;
    const hsOk = clean(hs), tpOk = clean(tp);
    Object.assign(bus, { 'site.Hs1': ex.v1, 'site.Hs100': ex.v100, 'site.Tp1': a * Math.sqrt(ex.v1), 'site.Tp100': a * Math.sqrt(ex.v100), 'site.Hmax100': 1.86 * ex.v100 });
    const pairs = hs.map((h, i) => [h, tp[i]]).filter(([h, p]) => Number.isFinite(h) && Number.isFinite(p) && p > 0);
    let scatter = null;
    try { scatter = scatterDiagram(Float64Array.from(pairs.map((p) => p[0])), Float64Array.from(pairs.map((p) => p[1])), { hsBin: 0.5, tpBin: 1 }); } catch { /* optional */ }
    // monthly climate
    const monthly = Array.from({ length: 12 }, () => []);
    t.forEach((s, i) => { if (Number.isFinite(hs[i])) monthly[Number(s.slice(5, 7)) - 1].push(hs[i]); });
    data.waves = { years: r.years, n: hsOk.length, mean: mean(hsOk), p50: quant(hsOk, 0.5), p90: quant(hsOk, 0.9), p99: thr, max: maxOf(hsOk), tpMean: tpOk.length ? mean(tpOk) : null, tpFactor: a, tpFrom, extremes: ex, scatter,
      monthlyMean: monthly.map(mean), monthlyP90: monthly.map((m) => quant(m, 0.9)), exceed: [0.5, 1, 1.5, 2, 2.5, 3, 4, 5].map((h) => ({ hs: h, frac: hsOk.filter((v) => v > h).length / hsOk.length })) };
    data.provenance.waves = { source: `Open-Meteo Marine hindcast (${r.grid}), ${r.years} years hourly`, at: r.at, method: `${ex.method}; ${tpFrom === 'none' ? `Tp = ${a.toFixed(2)}·√Hs, a generic storm-sea relation, because the hindcast carries no wave period here — replace with measured periods, especially where swell dominates` : `Tp = ${a.toFixed(2)}·√Hs fitted to the storms above the 99th percentile${tpFrom === 'mean' ? ' (peak period taken as 1.15 × the hindcast mean period)' : ''}`}; Hmax = 1.86·Hs (Rayleigh, ~1000 waves). Global reanalysis wave models are known to under-predict storm peaks, typically by 10–20 %: confirm extremes against a site-specific hindcast or measurements` };
  });
  await step('Wind hindcast', 0.48, async () => {
    const r = await C.windHindcast(lat, lon, { years, force, onProgress: (f, m) => onProgress?.(0.48 + 0.3 * f, m) });
    const u10 = r.data.wind_speed_10m, u100 = r.data.wind_speed_100m, t = r.data.time;
    const ex = extremes(t, u10, 'U10');
    const wb = weibullFit(u100);
    if (ex) Object.assign(bus, { 'site.Uw1': ex.v1, 'site.Uw100': ex.v100 });
    if (wb) Object.assign(bus, { 'site.windWeibullA': wb.A, 'site.windWeibullK': wb.k, 'site.hubWindMean': wb.mean });
    const temp = clean(r.data.temperature_2m);
    if (temp.length) bus['site.Tair'] = mean(temp);
    const ghi = clean(r.data.shortwave_radiation);
    if (ghi.length) bus['site.ghiAnnual'] = mean(ghi) * 8.76; // mean W/m² × 8760 h / 1000
    // 12-sector wind rose at 100 m
    const rose = Array.from({ length: 12 }, () => ({ n: 0, sum: 0 }));
    const dir = r.data.wind_direction_100m;
    for (let i = 0; i < u100.length; i++) if (Number.isFinite(u100[i]) && Number.isFinite(dir[i])) { const s = Math.floor(((dir[i] + 15) % 360) / 30); rose[s].n++; rose[s].sum += u100[i]; }
    const total = rose.reduce((s, x) => s + x.n, 0) || 1;
    const monthly = Array.from({ length: 12 }, () => []);
    t.forEach((s, i) => { if (Number.isFinite(u100[i])) monthly[Number(s.slice(5, 7)) - 1].push(u100[i]); });
    data.wind = { years: r.years, mean10: mean(clean(u10)), mean100: wb?.mean, weibull: wb, extremes: ex, shear: Math.log((wb?.mean || 1) / (mean(clean(u10)) || 1)) / Math.log(10),
      rose: rose.map((s, i) => ({ dir: i * 30, freq: s.n / total, mean: s.n ? s.sum / s.n : 0 })), monthlyMean: monthly.map(mean), tAirMin: temp.length ? minOf(temp) : null, tAirMax: temp.length ? maxOf(temp) : null };
    data.provenance.wind = { source: `Open-Meteo Historical Weather (${r.grid}), ${r.years} years hourly`, at: r.at, method: `${ex?.method || 'n/a'} on hourly-mean wind at 10 m; Weibull by maximum likelihood at 100 m` };
  });
  await step('Currents, tide and sea temperature', 0.8, async () => {
    const r = await C.oceanRecord(lat, lon, { force });
    const hh = r.data.hourly || {};
    const uc = clean(hh.ocean_current_velocity || []).map((v) => v / 3.6); // source unit km/h
    const sl = clean(hh.sea_level_height_msl || []);
    const sst = clean(hh.sea_surface_temperature || []);
    if (uc.length) {
      const umax = maxOf(uc);
      // one year of record cannot support a 100-year fit: scale the observed maximum by a stated factor instead
      Object.assign(bus, { 'site.Uc1': umax, 'site.Uc100': 1.25 * umax, 'site.UcBed100': 1.25 * umax * 0.6 });
      data.current = { mean: mean(uc), p90: quant(uc, 0.9), max: umax };
    }
    if (sl.length) {
      const range = quant(sl, 0.995) - quant(sl, 0.005);
      bus['site.tideRange'] = range;
      bus['site.surge100'] = Math.max(0.3, 0.25 * range);
      data.tide = { range, hat: maxOf(sl), lat: minOf(sl) };
    }
    if (sst.length) { bus['site.Tsea'] = mean(sst); data.sst = { mean: mean(sst), min: minOf(sst), max: maxOf(sst) }; }
    data.provenance.ocean = { source: 'Open-Meteo Marine (Copernicus Marine models), last 12 months hourly', at: r.at, method: '1-yr current = record maximum; 100-yr = 1.25 × record maximum (short record — confirm with site measurements); near-bed = 0.6 × surface; tidal range = 99 % spread of total sea level; surge allowance = max(0.3 m, 25 % of range)' };
  });
  await step('Seismicity', 0.88, async () => {
    const r = await C.earthquakes(lat, lon, { force });
    const s = seismicScreen(r.data, lat, lon);
    data.seismic = s;
    bus['site.pga475'] = s.pga475 * 9.80665;
    data.provenance.seismic = { source: 'USGS earthquake catalogue, M ≥ 4.5 within 300 km since 1970', at: r.at, method: s.note };
  });
  if (bus['site.depth'] != null) await step('Seabed sediment', 0.9, async () => {
    // best first: mapped substrate (European seas), described samples (United States), measured grain size from a
    // national survey (Australia, Canada); where none of them covers the site, the worldwide sediment map
    const got = await Promise.allSettled([C.seabedSubstrate(lat, lon, { force }), C.seabedSamples(lat, lon, { force }), C.seabedRegional(lat, lon, { force }), C.seabedWorldMap(lat, lon, { force })]);
    const map = got[0].status === 'fulfilled' ? got[0].value.data : null;
    const red = got[1].status === 'fulfilled' ? reduceSeabedSamples(got[1].value.data, lat, lon) : null;
    const reg = got[2].status === 'fulfilled' && got[2].value ? got[2].value.data : null;
    const lab = reg ? reduceSeabedSamples(gradeSeabedSamples(reg.samples), lat, lon) : null;
    const mapCls = map && classifySeabed(map.original && classifySeabed(map.original) ? map.original : map.text);
    // which grains the waves and currents of the site move: shown with every source, decisive only for the world map.
    // Shelf depths only: the surface current says nothing about the bed of the deep sea
    const hydro = bus['site.depth'] >= SHELF_DEPTH ? null : hydroSorting({ depth: bus['site.depth'], waves: data.waves && { p90: data.waves.p90, tpFactor: data.waves.tpFactor, hs1: bus['site.Hs1'], tp1: bus['site.Tp1'] }, current: data.current });
    const SURFACE = 'Surface sediment only.', MID = 'median grain size is the measured typical value of the class (see the card for its range and source)';
    if (mapCls) {
      data.seabed = { cls: mapCls, ...SEABED_CLASSES[mapCls], basis: 'map', described: map.original || map.text, scale: map.scale, holder: [map.holder, map.country].filter(Boolean).join(', '), hydro };
      data.provenance.seabed = { source: `EMODnet Geology seabed substrate map, 1:${Number(map.scale).toLocaleString('en')}${map.holder ? ` (${map.holder})` : ''}`, at: got[0].value.at, method: `Folk class of the mapped polygon at the site; ${MID}. ${SURFACE}` };
    } else if (red) {
      data.seabed = { cls: red.cls, ...SEABED_CLASSES[red.cls], basis: 'samples', n: red.n, agree: red.agree, nearestKm: red.nearestKm, farthestKm: red.farthestKm, years: red.years, described: red.samples.slice(0, 3).map((x) => x.text).join(' · '), hydro };
      data.provenance.seabed = { source: `NOAA NOS seabed samples (${red.n} within ${red.farthestKm.toFixed(1)} km, surveyed ${red.years[0]}–${red.years[1]})`, at: got[1].value.at, method: `Commonest class among the nearest described bottom samples, weighted by distance; ${MID}. ${SURFACE}` };
    } else if (lab) {
      const carb = lab.samples.map((x) => x.carbonate).filter((v) => v != null), c = SEABED_CLASSES[lab.cls];
      data.seabed = { cls: lab.cls, ...c, basis: 'grainsize', survey: C.source(reg.provider).name, n: lab.n, agree: lab.agree, nearestKm: lab.nearestKm, farthestKm: lab.farthestKm, years: lab.years, carbonate: carb.length ? quant(carb, 0.5) : null, described: lab.samples.slice(0, 3).map((x) => x.text).join(' · '), hydro,
        d50: lab.d50 ?? c.d50, d50Range: lab.d50 ? lab.d50Range : c.d50Range, d50Measured: lab.d50 ? lab.nSized : 0 };
      data.provenance.seabed = { source: `${C.source(reg.provider).name} (${C.source(reg.provider).provider}; ${lab.n} analysed samples within ${lab.farthestKm.toFixed(1)} km)`, at: got[2].value.at,
        method: `Class from the measured gravel, sand and mud fractions (Folk scheme; gravel from 30 %, then the sand-to-mud ratio) of the nearest surface samples, commonest class weighted by distance; ${lab.d50 ? `median grain size is the median of the mean grain size measured on ${lab.nSized} of those samples` : MID}. Point samples, resolution set by their spacing. ${SURFACE}` };
    } else {
      const live = got[3].status === 'fulfilled' ? got[3].value.data : null;
      const grid = (await C.seabedGrid(lat, lon)).data;
      const failed = got.filter((g) => g.status === 'rejected');
      if (failed.length) data.issues.push(`Seabed sediment: ${failed.length} of the live sources did not answer (${failed[0].reason.message}); the result below does not include them.`);
      const sb = worldSeabed({ depth: bus['site.depth'], live, grid, hydro });
      if (!sb) { data.seabed = { cls: null, basis: 'none', hydro }; data.seabed.confidence = seabedConfidence(data.seabed).key; return; }
      data.seabed = sb;
      const SHOM = 'Shom, 2026. Worldwide seabed sediment map, 1:500 000, doi:10.17183/HOM_GEOL_SEDIM_MONDIALE (Open Licence 2.0)', DEEP = 'Diesing 2020, Deep-sea sediments of the global ocean, PANGAEA, doi:10.1594/PANGAEA.911692 (CC BY 4.0), bundled 0.1° grid';
      data.provenance.seabed = { source: [sb.described && `${SHOM}${sb.basis === 'grid' ? ', bundled 0.1° grid' : ''}`, sb.composition && DEEP].filter(Boolean).join('; '), at: sb.basis === 'worldmap' ? got[3].value.at : null,
        method: `${sb.described ? `Map class “${sb.described}” (${sb.definition}) ${sb.basis === 'worldmap' ? `of the polygon at the site${sb.sheet != null ? `, source map no. ${sb.sheet} of the compilation` : ''}` : `at the ${sb.offsetKm ? `closest grid node with data, ${sb.offsetKm.toFixed(0)} km away` : 'nearest grid node'} (the map service did not answer)`}. ` : ''}${sb.composition ? `Deep-sea sediment type “${sb.composition}” at the nearest node of a 10 km machine-learning prediction. ` : ''}The map is a compilation of published sediment maps generalised to 1:500 000 and of a 1972 world map in the deep sea: it shows the sediment province, not the ground at the site. ${sb.shelf
          ? `In less than ${SHELF_DEPTH} m of water a grain size is passed to the analyses only when the hydrodynamic check agrees with the map: ${sb.evidence.verdict === 'agree' ? 'it does' : sb.evidence.verdict === 'conflict' ? 'it contradicts the map, so none is passed' : 'it cannot decide here, so none is passed'}. ${sb.evidence.why} ${hydro ? hydro.text : ''}`
          : 'Median grain size is the measured typical value of the class.'} ${SURFACE}` };
    }
    data.seabed.confidence = seabedConfidence(data.seabed).key;
    if (data.seabed.d50) bus['site.d50'] = data.seabed.d50;
  });
  await step('Solar resource', 0.93, async () => {
    const r = await C.solarClimatology(lat, lon, { force });
    const p = r.data.properties?.parameter || {};
    if (p.ALLSKY_SFC_SW_DWN?.ANN > 0) { bus['site.ghiAnnual'] = p.ALLSKY_SFC_SW_DWN.ANN * 365; data.solar = { ghiDaily: p.ALLSKY_SFC_SW_DWN, t2m: p.T2M }; }
    data.provenance.solar = { source: 'NASA POWER climatology', at: r.at };
  });
  // seabed temperature: a documented rule of thumb, since no open global point service exposes it
  if (bus['site.depth'] != null && bus['site.Tsea'] != null) {
    const d = bus['site.depth'];
    bus['site.Tbed'] = d < 50 ? bus['site.Tsea'] - 1 : Math.max(2, 4 + (bus['site.Tsea'] - 4) * Math.exp(-(d - 50) / 400));
    data.provenance.tbed = { source: 'derived', method: 'Seabed temperature decays exponentially from the sea-surface mean towards 4 °C with a 400 m scale below the mixed layer (screening value — replace with measured CTD data).' };
  }
  onProgress?.(1, 'Done');

  const p = getProject();
  // a screening soil profile was chosen for one site and one kind of seabed: it does not survive a change of either
  const before = p.site.data?.seabed?.cls, moved = p.site.lat != null && (Math.abs(p.site.lat - lat) > 0.01 || Math.abs(p.site.lon - lon) > 0.01);
  if (moved || (before && before !== data.seabed?.cls)) for (const k of SCREENING_KEYS) if (p.bus[k]?.from === 'screening') { delete p.bus[k]; data.screeningDropped = true; }
  p.site = { ...p.site, lat, lon, name: p.site.name || data.place?.name || `${lat.toFixed(3)}, ${lon.toFixed(3)}`, pulled: data.pulledAt, data };
  save();
  // values the user entered by hand (measurements, site investigation) outrank any pulled estimate
  for (const k of Object.keys(bus)) if (p.bus[k]?.from === 'user') { data.kept = [...(data.kept || []), k]; delete bus[k]; }
  publish(bus, 'site');
  emit('site', p.site);
  return { data, bus };
}

/** Site values a user enters by hand go on the same bus, so measured data always overrides the pulled estimate. */
export function setSiteValues(values) { publish(values, 'user'); }

export const SITE_KEYS = [
  ['site.depth', 'Water depth', 'm'], ['site.Hs1', '1-year significant wave height', 'm'], ['site.Tp1', '1-year peak period', 's'],
  ['site.Hs100', '100-year significant wave height', 'm'], ['site.Tp100', '100-year peak period', 's'], ['site.Hmax100', '100-year maximum wave height', 'm'],
  ['site.Uw1', '1-year wind (1 h, 10 m)', 'm/s'], ['site.Uw100', '100-year wind (1 h, 10 m)', 'm/s'], ['site.Uc1', '1-year surface current', 'm/s'],
  ['site.Uc100', '100-year surface current', 'm/s'], ['site.UcBed100', '100-year near-bed current', 'm/s'], ['site.tideRange', 'Tidal range', 'm'],
  ['site.surge100', '100-year storm surge', 'm'], ['site.rhoW', 'Seawater density', 'kg/m³'], ['site.Tsea', 'Sea-surface temperature', '°C'],
  ['site.Tbed', 'Seabed temperature', '°C'], ['site.Tair', 'Mean air temperature', '°C'], ['site.pga475', 'Peak ground acceleration (475-yr)', 'm/s²'],
  ['site.windWeibullA', 'Weibull scale A (100 m)', 'm/s'], ['site.windWeibullK', 'Weibull shape k', '-'], ['site.hubWindMean', 'Mean wind at 100 m', 'm/s'],
  ['site.ghiAnnual', 'Global horizontal irradiation', 'kWh/m²/yr'], ['site.su0', 'Seabed undrained shear strength', 'Pa'], ['site.suGrad', 'Shear-strength gradient', 'Pa/m'],
  ['site.phiDeg', 'Sand friction angle', 'deg'], ['site.gammaSub', 'Submerged unit weight', 'N/m³'], ['site.d50', 'Median grain size', 'm'],
];
