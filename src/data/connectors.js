// Live-data connectors. Every request goes straight from the user's browser to the public source
// (all listed services send CORS headers), so freshness never depends on any server of ours.
// Responses are cached in IndexedDB with a time-to-live and reused when the device is offline.
import { idb, emit } from '../core/store.js';

const HOUR = 3600e3, DAY = 24 * HOUR;

async function cached(key, ttl, loader, { force = false } = {}) {
  const hit = await idb.get('cache', key);
  const fresh = hit && Date.now() - hit.at < ttl;
  if (hit && fresh && !force) return { ...hit, stale: false, source: 'cache' };
  if (typeof navigator !== 'undefined' && navigator.onLine === false) {
    if (hit) return { ...hit, stale: true, source: 'offline-cache' };
    throw new Error('Offline and nothing cached for this request yet.');
  }
  try {
    const data = await loader();
    const rec = { at: Date.now(), data };
    await idb.set('cache', key, rec);
    emit('data', { key, at: rec.at });
    return { ...rec, stale: false, source: 'network' };
  } catch (err) {
    if (hit) return { ...hit, stale: true, source: 'stale-cache', error: err.message };
    throw err;
  }
}

async function getJson(url, { timeout = 30000 } = {}) {
  const ac = new AbortController();
  const t = setTimeout(() => ac.abort(), timeout);
  try {
    const r = await fetch(url, { signal: ac.signal, headers: { Accept: 'application/json' } });
    if (!r.ok) throw new Error(`${new URL(url).hostname} answered ${r.status}`);
    return await r.json();
  } catch (e) {
    throw new Error(e.name === 'AbortError' ? `${new URL(url).hostname} timed out` : e.message);
  } finally { clearTimeout(t); }
}
async function getText(url) { const r = await fetch(url); if (!r.ok) throw new Error(`${new URL(url).hostname} answered ${r.status}`); return r.text(); }
const q = (o) => Object.entries(o).filter(([, v]) => v != null).map(([k, v]) => `${k}=${encodeURIComponent(v)}`).join('&');
const r4 = (v) => Math.round(v * 1e4) / 1e4;
const iso = (d) => d.toISOString().slice(0, 10);

/** Catalogue shown on the Live data page; `suites` drives the per-suite wiring list. */
export const SOURCES = [
  { id: 'marine', name: 'Open-Meteo Marine', provider: 'Open-Meteo (ECMWF WAM, ERA5-Ocean, Copernicus Marine)', home: 'https://open-meteo.com/en/docs/marine-weather-api', gives: 'Wave height, period and direction (total, wind sea, swell), ocean currents, sea level incl. tide, sea-surface temperature — forecast and hindcast back to 1940', ttl: HOUR, licence: 'CC BY 4.0', suites: [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16] },
  { id: 'era5', name: 'Open-Meteo Historical Weather', provider: 'Open-Meteo (ECMWF ERA5 reanalysis)', home: 'https://open-meteo.com/en/docs/historical-weather-api', gives: 'Hourly wind at 10 m and 100 m, gusts, air temperature, pressure, solar radiation since 1940', ttl: 7 * DAY, licence: 'CC BY 4.0', suites: [1, 2, 3, 5, 6, 7, 11, 15, 16] },
  { id: 'forecast', name: 'Open-Meteo Forecast', provider: 'Open-Meteo (ECMWF, GFS, ICON …)', home: 'https://open-meteo.com/en/docs', gives: '16-day weather forecast for operations planning', ttl: HOUR, licence: 'CC BY 4.0', suites: [15, 16] },
  { id: 'climate', name: 'Open-Meteo Climate', provider: 'Open-Meteo (CMIP6 HighResMIP)', home: 'https://open-meteo.com/en/docs/climate-api', gives: 'Downscaled climate-model projections to 2050 for wind, temperature and radiation', ttl: 30 * DAY, licence: 'CC BY 4.0', suites: [7, 16, 17] },
  { id: 'bathy', name: 'NOAA NCEI global DEM mosaic', provider: 'NOAA National Centers for Environmental Information', home: 'https://www.ncei.noaa.gov/maps/bathymetry/', gives: 'Seabed elevation / water depth anywhere (best available of ETOPO, GEBCO-class and regional surveys), point and grid sampling', ttl: 90 * DAY, licence: 'Public domain', suites: [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16] },
  { id: 'emodnet', name: 'EMODnet Bathymetry', provider: 'European Marine Observation and Data Network', home: 'https://emodnet.ec.europa.eu/en/bathymetry', gives: 'High-resolution depth for European seas', ttl: 90 * DAY, licence: 'CC BY 4.0', suites: [3, 7, 8, 13, 16] },
  { id: 'elevation', name: 'Open-Meteo Elevation', provider: 'Open-Meteo (Copernicus GLO-90 DEM)', home: 'https://open-meteo.com/en/docs/elevation-api', gives: 'Land elevation for landfalls, quays and coastal works', ttl: 365 * DAY, licence: 'CC BY 4.0', suites: [5, 6, 7, 13] },
  { id: 'quakes', name: 'USGS Earthquake Catalog', provider: 'U.S. Geological Survey', home: 'https://earthquake.usgs.gov/fdsnws/event/1/', gives: 'Global earthquake catalogue (real time) for seismicity screening', ttl: 6 * HOUR, licence: 'Public domain', suites: [1, 2, 3, 5, 6, 7] },
  { id: 'power', name: 'NASA POWER', provider: 'NASA Langley Research Center', home: 'https://power.larc.nasa.gov/', gives: 'Solar irradiance and meteorological climatology for floating solar and hybrid systems', ttl: 30 * DAY, licence: 'Public domain', suites: [16] },
  { id: 'soilgrids', name: 'ISRIC SoilGrids', provider: 'ISRIC World Soil Information', home: 'https://soilgrids.org', gives: 'Onshore soil texture at landfalls and port sites (clay / sand / silt fractions, bulk density)', ttl: 365 * DAY, licence: 'CC BY 4.0', suites: [5, 6, 7, 13] },
  { id: 'place', name: 'OpenStreetMap Nominatim', provider: 'OpenStreetMap contributors', home: 'https://nominatim.org', gives: 'Place names, country and nearest coast context', ttl: 365 * DAY, licence: 'ODbL', suites: [17] },
  { id: 'geocode', name: 'Open-Meteo Geocoding', provider: 'Open-Meteo (GeoNames)', home: 'https://open-meteo.com/en/docs/geocoding-api', gives: 'Search any place name → coordinates', ttl: 365 * DAY, licence: 'CC BY 4.0', suites: [] },
  { id: 'fx', name: 'Frankfurter', provider: 'European Central Bank reference rates', home: 'https://frankfurter.dev', gives: 'Daily currency exchange rates', ttl: 12 * HOUR, licence: 'Open', suites: [17] },
  { id: 'worldbank', name: 'World Bank Open Data', provider: 'World Bank', home: 'https://data.worldbank.org', gives: 'Inflation, lending interest rate, GDP deflator and electricity indicators by country', ttl: 7 * DAY, licence: 'CC BY 4.0', suites: [17] },
  { id: 'commodities', name: 'IMF Primary Commodity Prices', provider: 'International Monetary Fund via DBnomics', home: 'https://db.nomics.world/IMF/PCPS', gives: 'Monthly Brent crude, natural gas, iron ore, copper, aluminium, coal prices', ttl: DAY, licence: 'IMF terms', suites: [8, 13, 17] },
  { id: 'treasury', name: 'U.S. Treasury Fiscal Data', provider: 'U.S. Department of the Treasury', home: 'https://fiscaldata.treasury.gov', gives: 'Average interest rates on marketable securities (risk-free rate proxy)', ttl: 7 * DAY, licence: 'Public domain', suites: [17] },
  { id: 'carbon', name: 'Carbon Intensity API', provider: 'National Grid ESO', home: 'https://carbonintensity.org.uk', gives: 'Live grid carbon intensity (Great Britain) for avoided-emission estimates', ttl: HOUR, licence: 'CC BY 4.0', suites: [16, 17] },
  { id: 'turbines', name: 'NREL turbine-models', provider: 'National Renewable Energy Laboratory (GitHub)', home: 'https://github.com/NREL/turbine-models', gives: 'Published power and thrust curves of reference offshore wind turbines', ttl: 30 * DAY, licence: 'BSD-3-Clause', suites: [3, 16] },
  { id: 'literature', name: 'OpenAlex', provider: 'OurResearch', home: 'https://openalex.org', gives: 'Newest peer-reviewed publications on any analysis topic, with open-access links', ttl: DAY, licence: 'CC0', suites: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17] },
  { id: 'crossref', name: 'Crossref', provider: 'Crossref', home: 'https://www.crossref.org', gives: 'DOI metadata for standards-related literature and reports', ttl: DAY, licence: 'Open', suites: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17] },
  { id: 'wiki', name: 'Wikipedia', provider: 'Wikimedia Foundation', home: 'https://en.wikipedia.org', gives: 'Plain-language definitions for the built-in glossary', ttl: 30 * DAY, licence: 'CC BY-SA 4.0', suites: [] },
];
export const source = (id) => SOURCES.find((s) => s.id === id);
export async function lastFetch(prefix) {
  const keys = (await idb.keys('cache')).filter((k) => String(k).startsWith(`${prefix}:`));
  let at = 0;
  for (const k of keys) { const r = await idb.get('cache', k); if (r?.at > at) at = r.at; }
  return at || null;
}

// ---------------- metocean ----------------
const MARINE = 'https://marine-api.open-meteo.com/v1/marine';
const MARINE_VARS = 'wave_height,wave_peak_period,wave_direction,wind_wave_height,swell_wave_height,swell_wave_peak_period,ocean_current_velocity,ocean_current_direction,sea_level_height_msl,sea_surface_temperature';

export const marineNow = (lat, lon, o) => cached(`marine:now:${r4(lat)},${r4(lon)}`, HOUR, () =>
  getJson(`${MARINE}?${q({ latitude: lat, longitude: lon, current: MARINE_VARS, hourly: MARINE_VARS, forecast_days: 8, timezone: 'GMT', cell_selection: 'sea' })}`), o);

/** Hourly wave hindcast, fetched in 5-year blocks so a slow link never loses the whole request. */
export async function marineHindcast(lat, lon, { years = 20, onProgress, force } = {}) {
  const endYear = new Date().getUTCFullYear() - 1;
  const out = { time: [], wave_height: [], wave_peak_period: [], wave_direction: [] };
  const blocks = [];
  for (let y = endYear - years + 1; y <= endYear; y += 5) blocks.push([y, Math.min(endYear, y + 4)]);
  let at = 0, stale = false, i = 0;
  for (const [a, b] of blocks) {
    const res = await cached(`marine:hind:${r4(lat)},${r4(lon)}:${a}-${b}`, 365 * DAY, () =>
      getJson(`${MARINE}?${q({ latitude: lat, longitude: lon, start_date: `${a}-01-01`, end_date: `${b}-12-31`, hourly: 'wave_height,wave_peak_period,wave_direction', timezone: 'GMT', cell_selection: 'sea' })}`, { timeout: 90000 }), { force });
    const hh = res.data.hourly || {};
    for (const k of Object.keys(out)) if (hh[k]) out[k].push(...hh[k]);
    at = Math.max(at, res.at); stale ||= res.stale;
    onProgress?.(++i / blocks.length, `Wave hindcast ${a}–${b}`);
  }
  return { at, stale, data: out, years, grid: 'ERA5-Ocean / ECMWF WAM' };
}

/** Recent ocean-current and sea-level record (the public archive holds the last few years). */
export const oceanRecord = (lat, lon, o) => {
  const end = new Date(Date.now() - 2 * DAY), start = new Date(end.getTime() - 365 * DAY);
  return cached(`marine:ocean:${r4(lat)},${r4(lon)}`, 30 * DAY, () =>
    getJson(`${MARINE}?${q({ latitude: lat, longitude: lon, start_date: iso(start), end_date: iso(end), hourly: 'ocean_current_velocity,ocean_current_direction,sea_level_height_msl,sea_surface_temperature', timezone: 'GMT', cell_selection: 'sea' })}`, { timeout: 90000 }), o);
};

export async function windHindcast(lat, lon, { years = 20, onProgress, force } = {}) {
  const endYear = new Date().getUTCFullYear() - 1;
  const out = { time: [], wind_speed_10m: [], wind_speed_100m: [], wind_direction_100m: [], wind_gusts_10m: [], temperature_2m: [], shortwave_radiation: [] };
  const blocks = [];
  for (let y = endYear - years + 1; y <= endYear; y += 5) blocks.push([y, Math.min(endYear, y + 4)]);
  let at = 0, stale = false, i = 0;
  for (const [a, b] of blocks) {
    const res = await cached(`era5:${r4(lat)},${r4(lon)}:${a}-${b}`, 365 * DAY, () =>
      getJson(`https://archive-api.open-meteo.com/v1/archive?${q({ latitude: lat, longitude: lon, start_date: `${a}-01-01`, end_date: `${b}-12-31`, hourly: Object.keys(out).filter((k) => k !== 'time').join(','), wind_speed_unit: 'ms', timezone: 'GMT' })}`, { timeout: 120000 }), { force });
    const hh = res.data.hourly || {};
    for (const k of Object.keys(out)) if (hh[k]) out[k].push(...hh[k]);
    at = Math.max(at, res.at); stale ||= res.stale;
    onProgress?.(++i / blocks.length, `Wind hindcast ${a}–${b}`);
  }
  return { at, stale, data: out, years, grid: 'ERA5' };
}

export const weatherForecast = (lat, lon, o) => cached(`forecast:${r4(lat)},${r4(lon)}`, HOUR, () =>
  getJson(`https://api.open-meteo.com/v1/forecast?${q({ latitude: lat, longitude: lon, hourly: 'wind_speed_10m,wind_gusts_10m,wind_direction_10m,temperature_2m,visibility,precipitation', wind_speed_unit: 'ms', forecast_days: 16, timezone: 'GMT' })}`), o);

export const climateProjection = (lat, lon, o) => cached(`climate:${r4(lat)},${r4(lon)}`, 30 * DAY, () =>
  getJson(`https://climate-api.open-meteo.com/v1/climate?${q({ latitude: lat, longitude: lon, start_date: '1990-01-01', end_date: '2050-12-31', models: 'MRI_AGCM3_2_S,EC_Earth3P_HR', daily: 'wind_speed_10m_mean,wind_speed_10m_max,temperature_2m_mean,shortwave_radiation_sum', wind_speed_unit: 'ms' })}`, { timeout: 120000 }), o);

// ---------------- seabed and land ----------------
const NCEI = 'https://gis.ngdc.noaa.gov/arcgis/rest/services/DEM_mosaics/DEM_global_mosaic/ImageServer';

export async function depthAt(lat, lon, o) {
  const res = await cached(`bathy:pt:${r4(lat)},${r4(lon)}`, 90 * DAY, async () => {
    const j = await getJson(`${NCEI}/identify?${q({ geometry: `${lon},${lat}`, geometryType: 'esriGeometryPoint', sr: 4326, returnGeometry: false, returnCatalogItems: false, f: 'json' })}`);
    const z = Number(j.value);
    if (!Number.isFinite(z)) throw new Error('No elevation value at this location');
    return { elevation: z };
  }, o);
  return res;
}

/** Bathymetry grid around a point: n × n samples over ±halfKm, returned as a CONTRACT §4 grid model (z up, metres). */
export async function bathyGrid(lat, lon, { halfKm = 10, n = 25, force } = {}) {
  return cached(`bathy:grid:${r4(lat)},${r4(lon)}:${halfKm}:${n}`, 90 * DAY, async () => {
    const dLat = halfKm / 111.32, dLon = halfKm / (111.32 * Math.cos(lat * Math.PI / 180));
    const z = new Array(n * n).fill(NaN);
    // the service caps the number of points per request, so sample one row at a time
    for (let j = 0; j < n; j++) {
      const pts = [];
      for (let i = 0; i < n; i++) pts.push([lon - dLon + 2 * dLon * i / (n - 1), lat - dLat + 2 * dLat * j / (n - 1)]);
      const geometry = JSON.stringify({ points: pts, spatialReference: { wkid: 4326 } });
      const res = await getJson(`${NCEI}/getSamples?${q({ geometry, geometryType: 'esriGeometryMultipoint', returnFirstValueOnly: true, f: 'json' })}`);
      for (const s of res.samples || []) z[j * n + s.locationId] = Number(s.value);
    }
    const dx = 2 * halfKm * 1000 / (n - 1);
    return { kind: 'grid', name: `Bathymetry ±${halfKm} km`, nx: n, ny: n, x0: -halfKm * 1000, y0: -halfKm * 1000, dx, dy: dx, z, crs: `local tangent plane at ${lat.toFixed(4)}, ${lon.toFixed(4)} (x east, y north)` };
  }, { force });
}

export const emodnetDepth = (lat, lon, o) => cached(`emodnet:${r4(lat)},${r4(lon)}`, 90 * DAY, () =>
  getJson(`https://rest.emodnet-bathymetry.eu/depth_sample?geom=${encodeURIComponent(`POINT(${lon} ${lat})`)}`), o);

export const landElevation = (lat, lon, o) => cached(`elevation:${r4(lat)},${r4(lon)}`, 365 * DAY, () =>
  getJson(`https://api.open-meteo.com/v1/elevation?${q({ latitude: lat, longitude: lon })}`), o);

export const soilGrids = (lat, lon, o) => cached(`soilgrids:${r4(lat)},${r4(lon)}`, 365 * DAY, () =>
  getJson(`https://rest.isric.org/soilgrids/v2.0/properties/query?lon=${lon}&lat=${lat}&property=clay&property=sand&property=silt&property=bdod&depth=0-5cm&depth=30-60cm&depth=100-200cm&value=mean`, { timeout: 60000 }), o);

export const earthquakes = (lat, lon, { radiusKm = 300, minMag = 4.5, since = '1970-01-01', force } = {}) => cached(`quakes:${r4(lat)},${r4(lon)}:${radiusKm}:${minMag}`, 6 * HOUR, () =>
  getJson(`https://earthquake.usgs.gov/fdsnws/event/1/query?${q({ format: 'geojson', latitude: lat, longitude: lon, maxradiuskm: radiusKm, minmagnitude: minMag, starttime: since, orderby: 'time', limit: 5000 })}`, { timeout: 60000 }), { force });

export const solarClimatology = (lat, lon, o) => cached(`power:${r4(lat)},${r4(lon)}`, 30 * DAY, () =>
  getJson(`https://power.larc.nasa.gov/api/temporal/climatology/point?${q({ parameters: 'ALLSKY_SFC_SW_DWN,T2M,WS10M,WS50M', community: 'RE', longitude: lon, latitude: lat, format: 'JSON' })}`, { timeout: 60000 }), o);

// ---------------- places ----------------
export const searchPlace = (name) => cached(`geocode:${name.toLowerCase()}`, 365 * DAY, () =>
  getJson(`https://geocoding-api.open-meteo.com/v1/search?${q({ name, count: 8, language: 'en' })}`));
export const reversePlace = (lat, lon) => cached(`place:${r4(lat)},${r4(lon)}`, 365 * DAY, () =>
  getJson(`https://nominatim.openstreetmap.org/reverse?${q({ lat, lon, format: 'json', zoom: 5 })}`));

// ---------------- economics ----------------
export const fxRates = (base = 'USD', o) => cached(`fx:${base}`, 12 * HOUR, () => getJson(`https://api.frankfurter.dev/v1/latest?base=${base}`), o);

export const WB_INDICATORS = { inflation: 'FP.CPI.TOTL.ZG', lendingRate: 'FR.INR.LEND', gdpDeflator: 'NY.GDP.DEFL.KD.ZG', realInterest: 'FR.INR.RINR' };
export const worldBank = (iso2, indicator, o) => cached(`worldbank:${iso2}:${indicator}`, 7 * DAY, async () => {
  const j = await getJson(`https://api.worldbank.org/v2/country/${iso2}/indicator/${indicator}?format=json&per_page=40`);
  return (j[1] || []).filter((r) => r.value != null).map((r) => ({ year: Number(r.date), value: r.value })).sort((a, b) => a.year - b.year);
}, o);

export const COMMODITIES = {
  brent: { code: 'POILBRE', label: 'Brent crude', unit: 'USD/bbl' },
  natgasEU: { code: 'PNGASEU', label: 'Natural gas, EU (TTF)', unit: 'USD/MMBtu' },
  natgasUS: { code: 'PNGASUS', label: 'Natural gas, US (Henry Hub)', unit: 'USD/MMBtu' },
  ironOre: { code: 'PIORECR', label: 'Iron ore (62 % Fe)', unit: 'USD/t' },
  copper: { code: 'PCOPP', label: 'Copper', unit: 'USD/t' },
  aluminium: { code: 'PALUM', label: 'Aluminium', unit: 'USD/t' },
  coal: { code: 'PCOALAU', label: 'Coal (Australia)', unit: 'USD/t' },
};
export const commodity = (key, o) => cached(`commodities:${key}`, DAY, async () => {
  const c = COMMODITIES[key];
  const j = await getJson(`https://api.db.nomics.world/v22/series/IMF/PCPS/M.W00.${c.code}.USD?observations=1`);
  const s = j.series?.docs?.[0];
  if (!s) throw new Error('Series not found');
  const rows = s.period.map((p, i) => ({ period: p, value: s.value[i] })).filter((r) => typeof r.value === 'number');
  return { ...c, rows: rows.slice(-240) };
}, o);

export const treasuryRates = (o) => cached('treasury:avg', 7 * DAY, () =>
  getJson('https://api.fiscaldata.treasury.gov/services/api/fiscal_service/v2/accounting/od/avg_interest_rates?sort=-record_date&page%5Bsize%5D=20'), o);
export const gridCarbon = (o) => cached('carbon:gb', HOUR, () => getJson('https://api.carbonintensity.org.uk/intensity'), o);

// ---------------- engineering references ----------------
export const turbineList = (o) => cached('turbines:list', 30 * DAY, async () => {
  const j = await getJson('https://api.github.com/repos/NREL/turbine-models/contents/turbine_models/data/Offshore');
  return j.filter((f) => f.name.endsWith('.csv')).map((f) => ({ name: f.name.replace('.csv', ''), url: f.download_url }));
}, o);
export const turbineCurve = (item) => cached(`turbines:${item.name}`, 30 * DAY, async () => {
  const txt = await getText(item.url);
  const lines = txt.trim().split(/\r?\n/);
  const head = lines[0].split(',').map((s) => s.trim());
  const col = (re) => head.findIndex((c) => re.test(c));
  const iu = col(/wind.*speed/i), ip = col(/^power/i), ict = col(/^ct|thrust.*coef/i), icp = col(/^cp/i);
  const rows = lines.slice(1).map((l) => l.split(',').map(Number)).filter((r) => Number.isFinite(r[iu]));
  return { name: item.name, head, U: rows.map((r) => r[iu]), P: ip >= 0 ? rows.map((r) => r[ip]) : null, Ct: ict >= 0 ? rows.map((r) => r[ict]) : null, Cp: icp >= 0 ? rows.map((r) => r[icp]) : null };
});

/** Newest literature on a topic (sorted by publication date), with an open-access link where one exists. */
export const literature = (query, { n = 8, force } = {}) => cached(`literature:${query.toLowerCase()}:${n}`, DAY, async () => {
  const j = await getJson(`https://api.openalex.org/works?${q({ search: query, 'per-page': n, sort: 'publication_date:desc', filter: 'type:article,has_abstract:true', select: 'title,publication_year,publication_date,doi,cited_by_count,primary_location,open_access,authorships' })}`);
  return (j.results || []).map((w) => ({ title: w.title, year: w.publication_year, date: w.publication_date, doi: w.doi, cited: w.cited_by_count, venue: w.primary_location?.source?.display_name || '', oa: w.open_access?.oa_url || null, authors: (w.authorships || []).slice(0, 3).map((a) => a.author?.display_name).filter(Boolean).join(', ') + ((w.authorships || []).length > 3 ? ' et al.' : '') }));
}, { force });

export const mostCited = (query, { n = 6 } = {}) => cached(`literature:top:${query.toLowerCase()}:${n}`, 7 * DAY, async () => {
  const j = await getJson(`https://api.openalex.org/works?${q({ search: query, 'per-page': n, sort: 'cited_by_count:desc', select: 'title,publication_year,doi,cited_by_count,primary_location,open_access' })}`);
  return (j.results || []).map((w) => ({ title: w.title, year: w.publication_year, doi: w.doi, cited: w.cited_by_count, venue: w.primary_location?.source?.display_name || '', oa: w.open_access?.oa_url || null }));
});

export const glossary = (term) => cached(`wiki:${term.toLowerCase()}`, 30 * DAY, () =>
  getJson(`https://en.wikipedia.org/api/rest_v1/page/summary/${encodeURIComponent(term.replace(/ /g, '_'))}`));

/** Probe every source once so the Live data page can show which are reachable right now. */
export async function health() {
  const probes = {
    marine: () => marineNow(57, 2, { force: true }), bathy: () => depthAt(57, 2, { force: true }), quakes: () => earthquakes(36, 140, { radiusKm: 100, minMag: 6, force: true }),
    fx: () => fxRates('USD', { force: true }), worldbank: () => worldBank('GB', WB_INDICATORS.inflation, { force: true }), commodities: () => commodity('brent', { force: true }),
    literature: () => literature('offshore structure', { n: 1, force: true }), forecast: () => weatherForecast(57, 2, { force: true }), carbon: () => gridCarbon({ force: true }),
    turbines: () => turbineList({ force: true }), treasury: () => treasuryRates({ force: true }), power: () => solarClimatology(4, 6, { force: true }),
  };
  const out = {};
  await Promise.all(Object.entries(probes).map(async ([id, fn]) => {
    const t0 = performance.now();
    try { const r = await fn(); out[id] = { ok: !r.stale, ms: Math.round(performance.now() - t0), at: r.at }; } catch (e) { out[id] = { ok: false, error: e.message }; }
  }));
  return out;
}
