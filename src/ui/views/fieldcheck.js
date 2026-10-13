// Field validation of the site-data layer: the modelled environment is compared with what instruments measured —
// sea level at the nearest tide gauge, and wave height, wind speed and sea temperature at the nearest buoy.
// These are independent physical measurements, so agreement here is validation in the proper sense, for the
// environmental inputs and for the place and period compared. It says nothing about any structural model.
import { h, clear, fmt, icon } from '../dom.js';
import { renderPlot } from '../plot.js';
import { getProject, save } from '../../core/store.js';
import { tideStations, tideGauge, iocStations, iocGauge, modelledSeaLevel, buoysNear, buoyRecord, modelledRecent } from '../../data/connectors.js';

const R = Math.PI / 180;
const km = (a, b) => { const x = (a.lon - b.lon) * R * Math.cos((a.lat + b.lat) * R / 2), y = (a.lat - b.lat) * R; return 6371 * Math.hypot(x, y); };
const mean = (a) => a.reduce((x, y) => x + y, 0) / a.length;

/** Agreement statistics for paired [measured, modelled] values. `removeMean` drops a datum offset first. */
export function compare(pairs, { removeMean = false } = {}) {
  const o = pairs.map((p) => p[0]), m = pairs.map((p) => p[1]);
  const mo = mean(o), mm = mean(m);
  const d = pairs.map((p) => (removeMean ? (p[1] - mm) - (p[0] - mo) : p[1] - p[0]));
  const so = Math.sqrt(mean(o.map((v) => (v - mo) ** 2))), sm = Math.sqrt(mean(m.map((v) => (v - mm) ** 2)));
  // the largest tenth of measured values: what matters for design is how the peaks are predicted
  const top = pairs.slice().sort((a, b) => b[0] - a[0]).slice(0, Math.max(3, Math.ceil(pairs.length / 10)));
  return { n: pairs.length, bias: removeMean ? 0 : mm - mo, offset: mm - mo, rmse: Math.sqrt(mean(d.map((x) => x * x))), corr: so && sm ? mean(pairs.map((p) => (p[0] - mo) * (p[1] - mm))) / (so * sm) : 0,
    meanObs: mo, meanModel: mm, maxObs: Math.max(...o), maxModel: Math.max(...m), rangeObs: Math.max(...o) - Math.min(...o), peakRatio: mean(top.map((p) => p[1])) / (mean(top.map((p) => p[0])) || 1e-12) };
}

function chart(title, unit, pairs, removeMean) {
  const st = compare(pairs, { removeMean });
  const step = Math.max(1, Math.floor(pairs.length / 400));
  const pick = pairs.filter((_, i) => i % step === 0);
  const xs = pick.map((_, i) => i * step);
  return renderPlot({ type: 'line', title, x: { label: 'Hours' }, y: { label: removeMean ? 'Level about the mean' : title.split(' — ')[0], unit }, series: [
    { name: 'Measured', x: xs, y: pick.map((p) => p[0] - (removeMean ? st.meanObs : 0)) }, { name: 'Modelled', x: xs, y: pick.map((p) => p[1] - (removeMean ? st.meanModel : 0)) }] });
}

async function seaLevel(site, out) {
  out.firstChild.textContent = 'Looking for the nearest tide gauge…';
  // two networks: NOAA (verified hourly data, United States) and IOC (real time, worldwide) — take whichever is nearer
  const lists = await Promise.allSettled([tideStations(), iocStations()]);
  const cands = [];
  if (lists[0].status === 'fulfilled') for (const a of lists[0].value.data) cands.push({ ...a, net: 'noaa', km: km(a, site) });
  if (lists[1].status === 'fulfilled') for (const a of lists[1].value.data) cands.push({ ...a, net: 'ioc', km: km(a, site) });
  cands.sort((a, b) => a.km - b.km);
  let lastErr = 'no gauge network answered';
  for (const g of cands.slice(0, 4)) {
    if (g.km > 400) break;
    try {
      out.firstChild.textContent = `Fetching the gauge at ${g.name} (${fmt(g.km, 3)} km away)…`;
      const [obs, mod] = await Promise.all([g.net === 'noaa' ? tideGauge(g.id, 30) : iocGauge(g.id, 20), modelledSeaLevel(g.lat, g.lon, 30)]);
      const m = new Map(mod.data.map((r) => [r.t.slice(0, 13), r.v]));
      const pairs = obs.data.map((r) => [r.v, m.get(r.t.slice(0, 13))]).filter((p) => Number.isFinite(p[0]) && typeof p[1] === 'number');
      if (pairs.length < 100) throw new Error('too few overlapping hours');
      const st = compare(pairs, { removeMean: true });
      if (st.rangeObs > 20) throw new Error('gauge record is not in metres or is corrupted');
      const good = st.corr > 0.9 && st.rmse < 0.15 * st.rangeObs;
      return { kind: 'Sea level', good, station: `${g.name}${g.country || g.state ? `, ${g.country || g.state}` : ''}`, network: g.net === 'noaa' ? 'NOAA CO-OPS' : 'IOC Sea Level Monitoring', km: g.km, ...st,
        text: `correlation ${fmt(st.corr, 3)}, root-mean-square difference ${fmt(st.rmse, 2)} m after removing the datum offset, measured range ${fmt(st.rangeObs, 3)} m`, fig: chart(`Sea level — ${g.name}`, 'm', pairs, true),
        caution: 'Treat the tidal range and surge values derived for this site with caution and prefer local measurements.' };
    } catch (e) { lastErr = `${g.name}: ${e.message}`; }
  }
  throw new Error(cands.length && cands[0].km <= 400 ? lastErr : `the nearest gauge is ${cands.length ? `${fmt(cands[0].km, 3)} km` : 'unknown distance'} away — too far to say anything about this site`);
}

async function buoy(site, out) {
  out.firstChild.textContent = 'Looking for a buoy reporting waves near the site…';
  // observing networks are dense off Europe and North America and sparse elsewhere: widen the search once
  let near = [];
  for (const halfDeg of [3, 9]) {
    try { near = (await buoysNear(site.lat, site.lon, { variable: 'wvht', halfDeg })).data.map((b) => ({ ...b, lat: b.latitude, lon: b.longitude })).map((b) => ({ ...b, km: km(b, site) })).sort((a, b) => a.km - b.km); } catch { near = []; }
    if (near.length) break;
  }
  if (!near.length) throw new Error('no buoy or ship reported wave height within about 1,000 km of this site in the last two days — observing networks are sparse in this region');
  const b = near[0];
  out.firstChild.textContent = `Fetching 30 days from platform ${b.platform_code} (${fmt(b.km, 3)} km away)…`;
  const [rec, mod] = await Promise.all([buoyRecord(b.platform_code), modelledRecent(b.lat, b.lon, 30)]);
  const m = new Map(mod.data.map((r) => [r.t, r]));
  const res = [];
  for (const [key, label, unit, note, tol] of [['wvht', 'Significant wave height', 'm', '', 0.25], ['windspd', 'Wind speed', 'm/s', ' (buoy anemometers sit a few metres above the sea; the model value is at 10 m, so a modelled value some 5–10 % higher is expected)', 0.25], ['sst', 'Sea-surface temperature', '°C', '', 0.1]]) {
    const pairs = rec.data.map((r) => [r[key], m.get(String(r.time).slice(0, 13))?.[key]]).filter((p) => Number.isFinite(p[0]) && typeof p[1] === 'number');
    if (pairs.length < 60) continue;
    const st = compare(pairs);
    const scale = key === 'sst' ? Math.max(1, st.rangeObs) : st.meanObs;
    const good = st.corr > (key === 'sst' ? 0.6 : 0.8) && st.rmse < tol * Math.max(scale, 1e-9) * (key === 'sst' ? 10 : 1.6) && Math.abs(st.peakRatio - 1) < 0.25;
    res.push({ kind: label, good, station: `platform ${b.platform_code} (${String(b.platform_type || '').toLowerCase()})`, network: 'NOAA OSMC', km: b.km, ...st,
      text: `correlation ${fmt(st.corr, 3)}, mean bias ${fmt(st.bias, 2)} ${unit}, root-mean-square difference ${fmt(st.rmse, 2)} ${unit}; on the largest tenth of measured values the model gives ${fmt(100 * st.peakRatio, 3)} % of what was measured${note}`,
      fig: chart(`${label} — platform ${b.platform_code}`, unit, pairs, false),
      caution: key === 'wvht' ? (st.peakRatio < 0.9 ? `The model under-predicts the highest waves here by about ${fmt(100 * (1 - st.peakRatio), 2)} %: raise the extreme wave heights derived for this site accordingly, or replace them with a site-specific hindcast.` : 'Differences are within what global wave models usually show; storms still deserve a site-specific hindcast.') : '' });
  }
  if (!res.length) throw new Error(`platform ${b.platform_code} has too few recent observations to compare`);
  return res;
}

/** Card for the site page. */
export function fieldChecks(site) {
  const out = h('div');
  const prev = getProject().site.data?.fieldChecks;
  const card = h('div.card', { style: { marginTop: '14px' } }, h('h3', 'Check the site data against field measurements'),
    h('p.muted', 'Compares the modelled environment with what instruments measured over the last weeks: sea level at the nearest tide gauge (about 1,300 gauges worldwide) and wave height, wind speed and sea temperature at the nearest reporting buoy. Real measurements, fetched live — the only kind of evidence that validates a model.'),
    h('div.row', h('button.btn.primary', { type: 'button', onclick: async (e) => {
      e.target.disabled = true; clear(out).append(h('small.muted', 'Starting…'));
      const results = [], problems = [];
      for (const fn of [seaLevel, buoy]) { try { const r = await fn(site, out); results.push(...[].concat(r)); } catch (err) { problems.push(err.message); } }
      clear(out);
      for (const r of results) out.append(h('div.banner', { class: r.good ? '' : 'warn' }, icon(r.good ? 'check' : 'warn'), h('div', h('b', `${r.kind}: ${r.good ? 'the model agrees with the measurements' : 'the model departs from the measurements'}. `),
        `${r.station}, ${fmt(r.km, 3)} km from the site (${r.network}); ${r.n} hours compared: ${r.text}. `, r.good ? '' : `${r.caution} `)));
      if (results.length) out.append(h('div.grid.c2', results.map((r) => r.fig)));
      for (const pr of problems) out.append(h('div.banner', icon('info'), `Not compared: ${pr}.`));
      out.append(h('p.why', 'What this establishes: how well the open ocean and weather models behind the site values reproduce measurements near this site in recent weeks. What it does not: the extremes used for design (weeks of data contain few storms), and anything about a structure. Buoys and gauges can themselves be faulty; a single poor comparison is a prompt to look, not a verdict.'));
      getProject().site.data.fieldChecks = { at: new Date().toISOString(), results: results.map(({ fig, ...r }) => r), problems }; save();
      e.target.disabled = false;
    } }, icon('check'), 'Compare with the nearest gauge and buoy')),
    prev ? h('p.muted', `Last comparison ${new Date(prev.at).toLocaleString()}: ${prev.results.map((r) => `${r.kind} ${r.good ? '✓' : '⚠'} (${r.station})`).join(' · ') || 'no station in range'}.`) : null, out);
  return card;
}
