// Site characterisation: pull every live source for a coordinate and reduce it to the design values the suites
// bind to (docs/CONTRACT.md §3). Each derived number keeps its provenance and method so it can be audited.
import * as C from './connectors.js';
import { fitGumbel, returnValue, scatterDiagram } from '../engine/reliability/extreme.js';
import { getProject, publish, save, emit } from '../core/store.js';

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

/** Water depth (m) taken as the edge of the continental shelf: shallower, the global lithology grid sets no grain size. */
export const SHELF_DEPTH = 200;
/**
 * Seabed-surface classes with a representative median grain size (Wentworth scale mid-points) and what each
 * usually means for foundations, scour and burial. Grain size is a property of the class; strength is not —
 * no strength is ever derived from a surface description.
 */
export const SEABED_CLASSES = {
  rock: { label: 'Rock or boulders', d50: null, means: 'Driven piles may refuse and trenching is slow: expect drilled or rock-socketed foundations and surface-laid, protected cables or pipelines.' },
  gravel: { label: 'Gravel or coarse-grained sediment', d50: 8e-3, means: 'Free-draining and resistant to scour; hard on ploughs and jetting tools, and suction installation may not seal.' },
  mixed: { label: 'Mixed sediment', d50: 1e-3, means: 'Variable over short distances: expect both granular and cohesive behaviour and plan the site investigation to map the boundaries.' },
  sand: { label: 'Sand', d50: 2.5e-4, means: 'Drained behaviour and good bearing, but mobile under waves and current: scour protection and burial depth need attention.' },
  muddySand: { label: 'Muddy sand', d50: 1.2e-4, means: 'Mostly granular with fines: partly drained response, moderate scour potential.' },
  sandyMud: { label: 'Sandy mud', d50: 4e-5, means: 'Mostly cohesive: undrained strength governs and consolidation settlement should be expected.' },
  mud: { label: 'Mud (silt and clay)', d50: 1.5e-5, means: 'Soft cohesive seabed: low bearing capacity near the surface, large embedment of pipelines and mudmats, favourable for suction and drag anchors.' },
  biogenic: { label: 'Shell, coral or other biogenic material', d50: 2e-3, means: 'Carbonate material crushes under load: pile shaft friction can be far below that of silica sand, and coral is usually an environmental constraint.' },
  // classes below occur only in the global lithology grid (continental margins and the deep sea)
  silt: { label: 'Silt', d50: 2e-5, means: 'Fine, low-plasticity sediment between drained and undrained behaviour: easily eroded and jetted, sensitive to cyclic loading, and slow to regain strength after disturbance around anchors and trenches.' },
  clay: { label: 'Clay (terrigenous or deep-sea pelagic clay)', d50: 2e-6, means: 'Very soft, highly plastic clay at the mudline with strength rising slowly with depth: suited to suction, drag-embedment and plate anchors; expect deep pipeline and mudmat embedment, long consolidation times and low axial pipe-soil resistance.' },
  calcOoze: { label: 'Calcareous ooze', d50: 2e-5, means: 'Soft carbonate sediment of hollow microfossil shells: the grains crush under load, so pile shaft friction and anchor capacity can be far below those of a silica soil of the same grading, and it is compressible and often lightly cemented in layers — design with carbonate-specific methods and laboratory tests.' },
  silOoze: { label: 'Siliceous ooze (diatoms, radiolarians or sponge spicules)', d50: 3e-5, means: 'Very porous, light sediment of fragile opal skeletons: unit weight is exceptionally low and water content high, the grains crush, and index-test correlations for ordinary clays do not apply — anchors and mudmats need site-specific strength testing.' },
  mixedOoze: { label: 'Mixed calcareous and siliceous ooze', d50: 2e-5, means: 'Soft biogenic sediment with crushable carbonate and opal grains: low unit weight, high compressibility and low pile or anchor friction; treat as a carbonate soil until tests show otherwise.' },
  volcanic: { label: 'Volcanic ash, sand or gravel', d50: 5e-4, means: 'Angular, vesicular grains that crush and may be altered to clay: friction angles look high in the laboratory but pile and anchor friction is unreliable, the ground is abrasive to trenching tools, and the area is likely to be seismically active.' },
  silMud: { label: 'Siliceous mud', d50: 1e-5, means: 'Soft cohesive mud rich in opal microfossils: lower unit weight and higher water content than ordinary clay of the same strength, large embedment of pipelines and mudmats, favourable for suction and drag anchors once the strength profile is measured.' },
  calcMud: { label: 'Fine-grained calcareous sediment', d50: 1e-5, means: 'Carbonate mud or marl: soft and compressible near the surface, with crushable grains and patchy cementation — pile friction is low and variable, and hard cemented layers can stop ploughs and suction-anchor skirts.' },
};
const SEABED_WORDS = [
  ['rock', /\b(rock|rk|r|rky|boulders?|blds?|bo|bedrock|ledge|lava|hard)\b/], ['biogenic', /\b(coral|co|crl|shells?|sh|oysters?|oys|maerl|sponge|weed|grass|kelp)\b/],
  ['gravel', /\b(gravel|gravelly|g|grv|pebbles?|p|cobbles?|cb|stones?|st|shingle|coarse-grained)\b/], ['mixed', /\bmixed\b/],
  ['sand', /\b(sand|sandy|sd|s)\b/], ['mud', /\b(mud|muddy|md|m|clay|cy|cl|silt|si|ooze|oz|marl|mrl)\b/],
];
/** Class of one seabed description: a Folk class name ("muddy Sand"), chart abbreviations ("SD,BRN,CRS") or plain words. */
export function classifySeabed(text) {
  const t = String(text || '').toLowerCase().replace(/[;/]/g, ',').trim();
  if (!t || /no data|unknown/.test(t)) return null;
  // Folk names: the noun (last word) is the main fraction, an adjective before it the secondary one
  if (/\bmud to muddy sand\b/.test(t)) return 'sandyMud';
  if (/\b(sandy|gravelly sandy) mud\b/.test(t)) return 'sandyMud';
  if (/\b(muddy|gravelly muddy) sand\b/.test(t)) return 'muddySand';
  if (/\b(sandy|muddy sandy|muddy) gravel\b/.test(t)) return 'gravel';
  if (/\bgravelly sand\b/.test(t)) return 'sand';
  if (/\bgravelly mud\b/.test(t)) return 'sandyMud';
  // otherwise the first recognised term of the description is the main constituent (chart convention)
  for (const part of t.split(',').map((x) => x.trim()).filter(Boolean)) {
    for (const [k, re] of SEABED_WORDS) if (re.test(part)) return k;
  }
  return null;
}
/** Class of a laboratory grain-size analysis: gravel, sand and mud fractions (any common scale), or failing those a Folk code such as "(g)mS". */
export function classifyGrainSize({ gravel, sand, mud, folk } = {}) {
  const g = Math.max(0, gravel || 0), s = Math.max(0, sand || 0), m = Math.max(0, mud || 0);
  if (g + s + m > 0) {
    if (g / (g + s + m) >= 0.3 || s + m === 0) return 'gravel';
    const r = s / (s + m);
    return r >= 0.9 ? 'sand' : r >= 0.5 ? 'muddySand' : r >= 0.1 ? 'sandyMud' : 'mud';
  }
  // Folk code: the capital is the main fraction, lower-case letters the modifiers, brackets "slightly"
  const code = /^((?:\(?[gms]\)?)*)([GSMZC])$/.exec(String(folk || '').trim());
  if (!code) return null;
  const mod = code[1].replace(/\([gms]\)/g, '');
  return code[2] === 'G' ? 'gravel' : code[2] === 'S' ? (mod.includes('m') ? 'muddySand' : 'sand') : /[sg]/.test(mod) ? 'sandyMud' : 'mud';
}
/** Grain-size analyses in the shape reduceSeabedSamples expects: a class and a one-line description each. */
export function gradeSeabedSamples(samples) {
  const pc = (v) => `${Math.round(v || 0)} %`;
  return samples.map((sm) => ({ ...sm, cls: classifyGrainSize(sm), text: sm.sand != null || sm.mud != null || sm.gravel != null ? `gravel ${pc(sm.gravel)}, sand ${pc(sm.sand)}, mud ${pc(sm.mud)}` : `Folk class ${sm.folk}` }));
}
/** Reduce nearby described samples to a class: the commonest class among the nearest few, with distances. */
export function reduceSeabedSamples(samples, lat, lon, { maxKm = 15, take = 7 } = {}) {
  const R = Math.PI / 180;
  const near = samples.map((sm) => ({ ...sm, cls: sm.cls || classifySeabed(sm.text), km: 6371 * Math.hypot((sm.lon - lon) * R * Math.cos(lat * R), (sm.lat - lat) * R) }))
    .filter((sm) => sm.cls && sm.km <= maxKm).sort((a, b) => a.km - b.km).slice(0, take);
  if (!near.length) return null;
  const votes = {};
  for (const sm of near) votes[sm.cls] = (votes[sm.cls] || 0) + 1 / (1 + sm.km);
  const cls = Object.entries(votes).sort((a, b) => b[1] - a[1])[0][0];
  return { cls, n: near.length, nearestKm: near[0].km, farthestKm: near[near.length - 1].km, agree: near.filter((sm) => sm.cls === cls).length, years: [Math.min(...near.map((sm) => sm.year || 9999)), Math.max(...near.map((sm) => sm.year || 0))], samples: near };
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
    const hs = r.data.wave_height, tp = r.data.wave_peak_period, t = r.data.time;
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
    data.waves = { years: r.years, n: hsOk.length, mean: mean(hsOk), p50: quant(hsOk, 0.5), p90: quant(hsOk, 0.9), p99: thr, max: maxOf(hsOk), tpMean: mean(tpOk), tpFactor: a, extremes: ex, scatter,
      monthlyMean: monthly.map(mean), monthlyP90: monthly.map((m) => quant(m, 0.9)), exceed: [0.5, 1, 1.5, 2, 2.5, 3, 4, 5].map((h) => ({ hs: h, frac: hsOk.filter((v) => v > h).length / hsOk.length })) };
    data.provenance.waves = { source: `Open-Meteo Marine hindcast (${r.grid}), ${r.years} years hourly`, at: r.at, method: `${ex.method}; Tp = ${a.toFixed(2)}·√Hs from storms above the 99th percentile; Hmax = 1.86·Hs (Rayleigh, ~1000 waves). Global reanalysis wave models are known to under-predict storm peaks, typically by 10–20 %: confirm extremes against a site-specific hindcast or measurements` };
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
    // national survey (Australia, Canada); where none of them covers the site, the bundled global lithology grid
    const got = await Promise.allSettled([C.seabedSubstrate(lat, lon, { force }), C.seabedSamples(lat, lon, { force }), C.seabedRegional(lat, lon, { force })]);
    const map = got[0].status === 'fulfilled' ? got[0].value.data : null;
    const red = got[1].status === 'fulfilled' ? reduceSeabedSamples(got[1].value.data, lat, lon) : null;
    const reg = got[2].status === 'fulfilled' && got[2].value ? got[2].value.data : null;
    const lab = reg ? reduceSeabedSamples(gradeSeabedSamples(reg.samples), lat, lon) : null;
    const mapCls = map && classifySeabed(map.original && classifySeabed(map.original) ? map.original : map.text);
    if (mapCls) {
      data.seabed = { cls: mapCls, ...SEABED_CLASSES[mapCls], basis: 'map', described: map.original || map.text, scale: map.scale, holder: [map.holder, map.country].filter(Boolean).join(', ') };
      data.provenance.seabed = { source: `EMODnet Geology seabed substrate map, 1:${Number(map.scale).toLocaleString('en')}${map.holder ? ` (${map.holder})` : ''}`, at: got[0].value.at, method: 'Folk class of the mapped polygon at the site; median grain size is the mid-point of the class. Surface sediment only.' };
    } else if (red) {
      data.seabed = { cls: red.cls, ...SEABED_CLASSES[red.cls], basis: 'samples', n: red.n, agree: red.agree, nearestKm: red.nearestKm, farthestKm: red.farthestKm, years: red.years, described: red.samples.slice(0, 3).map((x) => x.text).join(' · ') };
      data.provenance.seabed = { source: `NOAA NOS seabed samples (${red.n} within ${red.farthestKm.toFixed(1)} km, surveyed ${red.years[0]}–${red.years[1]})`, at: got[1].value.at, method: 'Commonest class among the nearest described bottom samples, weighted by distance; median grain size is the mid-point of the class. Surface sediment only.' };
    } else if (lab) {
      const carb = lab.samples.map((x) => x.carbonate).filter((v) => v != null);
      data.seabed = { cls: lab.cls, ...SEABED_CLASSES[lab.cls], basis: 'grainsize', survey: C.source(reg.provider).name, n: lab.n, agree: lab.agree, nearestKm: lab.nearestKm, farthestKm: lab.farthestKm, years: lab.years, carbonate: carb.length ? quant(carb, 0.5) : null, described: lab.samples.slice(0, 3).map((x) => x.text).join(' · ') };
      data.provenance.seabed = { source: `${C.source(reg.provider).name} (${C.source(reg.provider).provider}; ${lab.n} analysed samples within ${lab.farthestKm.toFixed(1)} km)`, at: got[2].value.at, method: 'Class from the measured gravel, sand and mud fractions (Folk scheme; gravel from 30 %, then the sand-to-mud ratio) of the nearest surface samples, commonest class weighted by distance; median grain size is the mid-point of the class. Point samples, resolution set by their spacing. Surface sediment only.' };
    } else {
      const grid = (await C.seabedGrid(lat, lon)).data;
      const failed = got.filter((g) => g.status === 'rejected');
      if (failed.length) data.issues.push(`Seabed sediment: ${failed.length} of the live sources did not answer (${failed[0].reason.message}); the result below does not include them.`);
      if (!grid) { data.seabed = { cls: null, basis: 'none' }; return; }
      const km = grid.stepDeg * 111.2, shelf = bus['site.depth'] < SHELF_DEPTH;
      data.seabed = { cls: grid.cls, ...SEABED_CLASSES[grid.cls], basis: 'grid', described: grid.name, resolutionDeg: grid.stepDeg, resolutionKm: km, offsetKm: grid.shifted ? grid.offsetKm : 0, shelf };
      data.provenance.seabed = { source: `Global seafloor lithology grid, ${grid.stepDeg}° (about ${km.toFixed(0)} km) — Dutkiewicz, Müller, O’Callaghan & Jónasson 2015, Geology 43(9), doi:10.1130/G36883.1; CC BY-SA 4.0, bundled with the app`,
        method: `Lithology “${grid.name}” at the ${grid.shifted ? `closest grid node with data (${grid.offsetKm.toFixed(0)} km from the site)` : 'nearest grid node'}. The grid is a machine-learning interpolation of about 14,500 widely spaced sample descriptions: it indicates the regional sediment province only and is least reliable on continental shelves, where sediment changes over short distances. Median grain size is a nominal value for the class${shelf ? ` and is not passed to the analyses in less than ${SHELF_DEPTH} m of water` : ''}. Surface sediment only.` };
      // on the shelf the interpolated class is too uncertain to drive scour and sediment-transport defaults
      if (shelf) return;
    }
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
