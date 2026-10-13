// Seabed and soil on the site page: what the open sources say about the seabed surface and how sure that is,
// the user's own observation, published studies of the area, and the two honest routes to soil strength — an
// imported cone penetration test or borehole, or an explicit screening profile until site data arrive.
import { h, clear, fmt, icon, toast, table } from '../dom.js';
import { renderPlot } from '../plot.js';
import { getProject, save, publish, emit } from '../../core/store.js';
import { SEABED_CLASSES, classifyObservation } from '../../data/seabed.js';
import { seabedConfidence, CONFIDENCE } from '../../data/seabedevidence.js';
import { screeningProfile, investigationPlan, STRUCTURES, SCREENING_KEYS, SCREENING_FLAG } from '../../data/soilscreening.js';
import { seabedStudies, seaAreas } from '../../data/connectors.js';
import { parseCptTable, cptDesignValues, boreholeDesignValues, PRESSURE_UNITS, LENGTH_UNITS } from '../../engine/geotech/siteinvestigation.js';

const mmOf = (d) => `${fmt(d * 1000, 2)} mm`;
const range = (r) => (r ? `${mmOf(r[0])} to ${mmOf(r[1])}` : '');
const KEY_LABEL = { 'site.su0': ['Undrained shear strength at the mudline', 'kPa', 1e-3], 'site.suGrad': ['Strength gain with depth', 'kPa/m', 1e-3], 'site.phiDeg': ['Sand friction angle', '°', 1], 'site.gammaSub': ['Submerged unit weight', 'kN/m³', 1e-3] };
// checked October 2026 at 88 shelf sites whose seabed is known from the survey maps and samples wired into the app
const HOW_WELL = 'How well this works: at 88 shelf sites where survey maps or samples give the answer, the seabed was of the kind stated (sandy or muddy) at 33 of the 36 sites where the world map and the hydrodynamic check agreed, and at 33 of the 45 where the check could not decide and the map stood alone. Where the two contradicted each other (7 sites) the map was wrong at 4.';

function seabedText(sb) {
  const span = (x) => `${fmt(x.nearestKm, 2)}–${fmt(x.farthestKm, 2)} km away${x.years?.[0] <= x.years?.[1] ? `, sampled ${x.years[0]}–${x.years[1]}` : ''}`;
  const world = () => `${sb.described ? `The Shom worldwide seabed sediment map (1:500 000) shows “${sb.described}” (${sb.definition}) ${sb.basis === 'worldmap' ? 'at the site' : `at the ${sb.offsetKm ? `closest node of its bundled 0.1° copy with data, ${fmt(sb.offsetKm, 2)} km from the site` : 'nearest node of its bundled 0.1° copy (the map service did not answer)'}`}. ` : ''}${sb.composition ? `The deep-sea sediment grid of Diesing (2020) predicts “${sb.composition}” here. ` : ''}${sb.described && sb.composition ? 'Both are regional compilations: they show' : 'This is a regional compilation: it shows'} the sediment province, not the ground at the site, and no survey map or sample in the open sources covers this location. `;
  return ({
    map: () => `Mapped at the site as “${sb.described}” on the EMODnet Geology substrate map at 1:${Number(sb.scale).toLocaleString('en')}${sb.holder ? ` (${sb.holder})` : ''}. `,
    samples: () => `${sb.agree} of the ${sb.n} nearest described bottom samples (${fmt(sb.nearestKm, 2)}–${fmt(sb.farthestKm, 2)} km away, surveyed ${sb.years[0]}–${sb.years[1]}) are of this class; nearest descriptions: ${sb.described}. `,
    grainsize: () => `${sb.agree} of the ${sb.n} nearest laboratory grain-size analyses in the ${sb.survey} (${span(sb)}) fall in this class; nearest results: ${sb.described}. These are point samples, so the answer is as local as their spacing. ${sb.carbonate >= 50 ? `Carbonate content is about ${fmt(sb.carbonate, 2)} %: treat the ground as a carbonate soil, whose grains crush and give low pile friction. ` : ''}`,
    worldmap: world, grid: world, deepsea: world,
  }[sb.basis] || (() => ''))();
}

function grainText(sb) {
  if (sb.basis === 'grainsize' && sb.d50Measured) return `Median grain size ${mmOf(sb.d50)}: ${sb.d50Measured > 1 ? `the median of the mean grain size measured on ${sb.d50Measured} of those samples (${range(sb.d50Range)})` : 'the mean grain size measured on the one sample of this class that has it'}, used as the default for scour and sediment-transport analyses. `;
  if (sb.d50) return `Median grain size ${mmOf(sb.d50)} is passed to the scour and sediment-transport analyses: the typical measured value for this class, which ranges from ${range(sb.d50Range)} — it is not a measurement at the site, and the analyses use the single value, so try both ends of the range. `;
  if (sb.shelf && sb.evidence) return h('b', `No grain size is passed to the analyses: ${sb.evidence.verdict === 'conflict' ? 'the two lines of evidence contradict each other' : 'on a continental shelf a coarse map needs a second line of evidence, and the hydrodynamic check cannot decide here'}. ${sb.d50Range ? `For this class measured values range from ${range(sb.d50Range)}${sb.d50Class ? ` (typically ${mmOf(sb.d50Class)})` : ''}. ` : ''}Enter one from a survey, a chart or a sample. `);
  return `No grain size is set for this class${sb.d50Range ? ` (measured values range from ${range(sb.d50Range)})` : ''}. `;
}

/** The “Seabed sediment” card. */
function seabedCard(redraw) {
  const p = getProject(), sb = p.site.data?.seabed, obs = p.site.observation;
  const conf = obs ? { key: 'observed', ...CONFIDENCE.observed } : seabedConfidence(sb);
  const card = h('div.card', { style: { marginTop: '14px' } }, h('div.row.spread', h('h3', { style: { margin: 0 } }, 'Seabed sediment'), h('span', { class: `chip ${conf.chip}`, title: 'How sure the statement below is' }, `Confidence: ${conf.label}`)));
  if (obs) {
    const c = SEABED_CLASSES[obs.cls];
    card.append(h('p', h('b', `${c.label} — your own observation. `), `${obs.note ? `“${obs.note}”. ` : ''}Classed ${obs.how}${obs.d50 ? `; median grain size ${mmOf(obs.d50)}${obs.measured ? ' as you measured it' : ' (typical measured value of the class)'} is passed to the analyses as your own value` : ''}. `, c.means));
  }
  if (sb?.cls) {
    card.append(h('p', obs ? h('b', 'The open sources say: ') : null, h('b', `${sb.label}. `), seabedText(sb), grainText(sb), obs ? '' : sb.means));
    if (sb.hydro) card.append(h('p', h('b', 'What the waves and currents allow. '), sb.hydro.text, sb.evidence ? h('span', ' ', h('span', { class: `chip ${sb.evidence.verdict === 'agree' ? 'ok' : sb.evidence.verdict === 'conflict' ? 'fail' : 'warn'}` }, sb.evidence.verdict === 'agree' ? 'consistent with the map — still not a measurement' : sb.evidence.verdict === 'conflict' ? 'contradicts the map' : 'cannot decide'), ` ${sb.evidence.why}`) : null));
    if (sb.evidence?.check?.length) card.append(h('p', { style: { marginBottom: 0 } }, h('b', 'What could explain it, and what to check:')), h('ul.plain', sb.evidence.check.map((x) => h('li', x))));
    if (sb.evidence) card.append(h('p.why', HOW_WELL));
  } else if (!obs) {
    card.append(h('p', 'No open seabed-sediment map or sample covers this location, and the worldwide map holds no value here (land, or sea beyond 72° of latitude). Record your own observation below, or enter the soil values from a site investigation or a regional geological map in the table above.'));
  }
  card.append(h('p.why', h('b', 'What would raise the confidence: '), conf.raise));
  if (sb?.cls) card.append(h('details.fold', h('summary', 'Where the grain size comes from'), h('div', h('p', sb.d50Source), sb.hydro ? h('p.muted', `Hydrodynamic check: orbital velocity at the bed from linear wave theory, bed shear stress from Soulsby (1997) with waves and current in line, threshold of motion from Soulsby & Whitehouse (1997). Conditions exceeded 10 % of the time: Hs ${fmt(sb.hydro.frequent.Hs, 2)} m, Tp ${fmt(sb.hydro.frequent.Tp, 2)} s, depth-mean current ${fmt(sb.hydro.frequent.U, 2)} m/s, giving ${fmt(sb.hydro.frequent.tau, 2)} Pa on a fine-sand bed; 1-year storm: Hs ${fmt(sb.hydro.storm.Hs, 2)} m, Tp ${fmt(sb.hydro.storm.Tp, 2)} s, current ${fmt(sb.hydro.storm.U, 2)} m/s, giving ${fmt(sb.hydro.storm.tau, 2)} Pa.`) : null)));
  card.append(observationFold(redraw), studiesFold());
  card.append(h('p.why', 'This describes the top few centimetres of the seabed. It indicates what kind of ground to expect and sets a grain size; it does not give strength, layering or depth to rock, and no strength value is derived from it. Strength comes only from site data or from the screening profile you choose below.'));
  return card;
}

/** The user's own seabed observation: stored in the project, published as the user's value. */
function observationFold(redraw) {
  const p = getProject(), o = p.site.observation || {};
  const f = (label, key, attrs = {}) => { const el = h('input', { type: 'text', inputmode: 'decimal', value: o.raw?.[key] ?? '', 'aria-label': label, ...attrs }); return [h('label', label, el), el]; };
  const note = h('input', { type: 'text', value: o.note || '', placeholder: 'e.g. grab sample: grey fine sand with shell fragments', 'aria-label': 'Description of the seabed' });
  const [gL, g] = f('Gravel %', 'gravel'), [sL, s] = f('Sand %', 'sand'), [mL, m] = f('Mud %', 'mud'), [dL, d] = f('Measured median grain size d50, mm', 'd50mm');
  const msg = h('small.muted');
  const apply = () => {
    const raw = { gravel: g.value, sand: s.value, mud: m.value, d50mm: d.value };
    const r = classifyObservation({ text: note.value, ...Object.fromEntries(Object.entries(raw).map(([k, v]) => [k, Number(String(v).replace(',', '.'))])) });
    if (!r) { msg.textContent = 'Nothing recognisable yet: describe the sediment in plain words (sand, mud, gravel, rock, shells…), or give percentages or a grain size.'; return; }
    p.site.observation = { ...r, note: note.value.trim(), raw, at: new Date().toISOString() };
    save();
    if (r.d50) publish({ 'site.d50': r.d50 }, 'user');
    toast(`Recorded: ${SEABED_CLASSES[r.cls].label}${r.d50 ? ` — grain size ${mmOf(r.d50)} is now your own value` : ''}.`);
    emit('site', p.site); redraw();
  };
  const remove = () => { delete p.site.observation; save(); toast('Observation removed. A grain size you set stays until you change it or pull the site again.'); if (p.bus['site.d50']?.from === 'user') { delete p.bus['site.d50']; save(); emit('bus', p.bus); } redraw(); };
  return h('details.fold', { open: !!p.site.observation }, h('summary', 'Record your own observation of the seabed'),
    h('div', h('p.muted', 'A grab sample, a diver’s or camera description, a chart bottom note or a laboratory grading. Any one of the fields is enough; a measured grain size outranks percentages, and percentages outrank the description. Your observation outranks every pulled source.'),
      h('div.soilgrid', h('label', { style: { gridColumn: '1 / -1' } }, 'Description', note), gL, sL, mL, dL),
      h('div.row', h('button.btn.primary', { type: 'button', onclick: apply }, 'Record observation'), p.site.observation ? h('button.btn', { type: 'button', onclick: remove }, 'Remove') : null, msg)));
}

/** Published seabed and geotechnical studies for the area, fetched on demand. */
function studiesFold() {
  const s = getProject().site, country = s.data?.place?.country;
  const area = h('input', { type: 'text', placeholder: 'Sea area, basin or field name', 'aria-label': 'Sea area to search for', style: { minWidth: '220px' } });
  const out = h('div'), picks = h('span.row');
  const go = async () => {
    const a = area.value.trim();
    if (!a) { clear(out).append(h('small.muted', 'Type the name of the sea area, basin or field.')); return; }
    clear(out).append(h('small.muted', 'Searching…'));
    try {
      const { via, items } = (await seabedStudies(a)).data;
      clear(out);
      if (!items.length) { out.append(h('small.muted', `Nothing found for “${a}”. Try the name geologists use for the area — a sea, basin, delta or field.`)); return; }
      out.append(h('ul.plain', items.map((w) => h('li', w.doi ? h('a', { href: w.doi, target: '_blank', rel: 'noopener' }, w.title) : w.title, ` — ${w.authors || 'authors not listed'}${w.venue ? `, ${w.venue}` : ''}${w.year ? ` (${w.year})` : ''}; cited ${w.cited} times`, w.oa ? h('span', ' · ', h('a', { href: w.oa, target: '_blank', rel: 'noopener' }, 'free full text')) : null))),
        h('small.muted', `Source: ${via}. A search result is a lead to read, not evidence about your site until you have read it.`));
    } catch (e) { clear(out).append(h('small.muted', `Search not available: ${e.message}`)); }
  };
  area.addEventListener('keydown', (e) => { if (e.key === 'Enter') go(); });
  // the names of the sea areas at the site are looked up when the section is first opened
  let asked = false;
  const fold = h('details.fold', h('summary', 'Published seabed and geotechnical studies for this area'),
    h('div', h('p.muted', 'Surveys of many shelves are published even where no open map service exists. Search by the name of the sea area, basin or field; the names found for the site are offered below.'),
      h('div.row', area, h('button.btn', { type: 'button', onclick: go }, icon('search'), 'Find studies'), picks), out));
  fold.addEventListener('toggle', async () => {
    if (!fold.open || asked || s.lat == null) return;
    asked = true;
    let names = [];
    try { names = (await seaAreas(Number(s.lat), Number(s.lon))).data; } catch { /* the field stays free text */ }
    if (country) names.push(`${country} offshore`);
    if (!area.value && names[0]) area.value = names[0];
    clear(picks).append(names.slice(1).map((n) => h('button.btn.small', { type: 'button', onclick: () => { area.value = n; go(); } }, n)));
  });
  return fold;
}

const kv = (rows) => table([{ key: 'k', label: 'Quantity', left: true }, { key: 'v', label: 'Value', left: true }, { key: 'how', label: 'How it was obtained', left: true }], rows);

/** Import of a cone penetration test or a layered borehole description. */
function investigationFold(redraw) {
  const p = getProject(), saved = p.site.soilData || {};
  const text = h('textarea.soiltext', { 'aria-label': 'CPT table or borehole layers', placeholder: 'depth [m], qc [MPa], fs [kPa], u2 [kPa]\n0.5, 0.31, 6.2, 12\n1.0, 0.34, 6.9, 21\n…' }, saved.text || '');
  const kind = h('select', { 'aria-label': 'Kind of record' }, [['cpt', 'Cone penetration test (CPT / CPTU)'], ['borehole', 'Borehole layers']].map(([v, l]) => h('option', { value: v, selected: saved.kind === v }, l)));
  const sel = (label, opts, def) => { const el = h('select', { 'aria-label': label }, Object.keys(opts).map((u) => h('option', { value: u, selected: u === def }, u))); return [h('label', label, el), el]; };
  const u = saved.units || {};
  const [zL, zU] = sel('Depth unit', LENGTH_UNITS, u.depth || 'm'), [qL, qU] = sel('Cone resistance unit', PRESSURE_UNITS, u.qc || 'MPa'), [fL, fU] = sel('Sleeve friction unit', PRESSURE_UNITS, u.fs || 'kPa'), [uL, uU] = sel('Pore pressure unit', PRESSURE_UNITS, u.u2 || 'kPa');
  const nkt = h('input', { type: 'text', inputmode: 'decimal', value: saved.Nkt || 15, 'aria-label': 'Cone factor Nkt' }), area = h('input', { type: 'text', inputmode: 'decimal', value: saved.a || 0.8, 'aria-label': 'Cone net area ratio' });
  const out = h('div');
  const file = h('input', { type: 'file', accept: '.csv,.txt,.gef,.ags,.dat,.tsv', style: { display: 'none' }, onchange: async () => { if (file.files[0]) { text.value = await file.files[0].text(); run(); } } });
  const hint = h('p.muted');
  const setHint = () => { hint.textContent = kind.value === 'cpt' ? 'Paste or upload the table: depth below the mudline, cone resistance (qc, or qt if already corrected), sleeve friction fs and, if measured, pore pressure u2 referenced to the seabed. Comma, semicolon, tab or blank separated; GEF files and the SCPT group of AGS files are read directly. Units written in the header are used, otherwise those chosen here.' : 'One layer per line: top depth m, bottom depth m, soil (clay or sand), then for clay the undrained strength at top and bottom in kPa, for sand the friction angle in degrees, and last the submerged unit weight in kN/m³. Example: 0, 6, clay, 4, 14, 6.5  —  6, 15, sand, 34, 9.5'; };
  kind.addEventListener('change', setHint); setHint();

  function interpret() {
    if (kind.value === 'borehole') {
      const layers = text.value.split(/\r?\n/).map((l) => l.trim()).filter((l) => l && /^\d/.test(l)).map((l) => l.split(/[,;\t]+/).map((c) => c.trim())).map((c) => {
        const soil = /sand|gravel/i.test(c[2]) ? 'sand' : 'clay', n = c.slice(3).map(Number);
        return soil === 'clay' ? { top: Number(c[0]), bottom: Number(c[1]), soil, suTop: n[0] * 1e3, suBottom: (Number.isFinite(n[1]) && c.length > 5 ? n[1] : n[0]) * 1e3, gammaSub: (c.length > 5 ? n[2] : n[1]) * 1e3 } : { top: Number(c[0]), bottom: Number(c[1]), soil, phiDeg: n[0], gammaSub: n[1] * 1e3 };
      });
      const d = boreholeDesignValues(layers);
      return { d, layers, rows: [
        d.su0 != null && { k: KEY_LABEL['site.su0'][0], v: `${fmt(d.su0 / 1e3, 3)} kPa`, how: 'Straight line fitted by least squares through the strength profile of the clay layers.' },
        d.su0 != null && { k: KEY_LABEL['site.suGrad'][0], v: `${fmt(d.suGrad / 1e3, 3)} kPa/m`, how: 'Slope of the same line.' },
        d.phiDeg != null && { k: KEY_LABEL['site.phiDeg'][0], v: `${fmt(d.phiDeg, 3)}°`, how: 'The lowest friction angle among the sand layers (cautious).' },
        d.gammaSub != null && { k: KEY_LABEL['site.gammaSub'][0], v: `${fmt(d.gammaSub / 1e3, 3)} kN/m³`, how: `Thickness-weighted mean of the ${d.share.sand > d.share.clay ? 'sand' : 'clay'} layers, which make up most of the profile.` }].filter(Boolean),
        plots: d.suFit ? [renderPlot({ type: 'line', title: 'Undrained shear strength with depth', x: { label: 'su', unit: 'kPa' }, y: { label: 'Depth below mudline', unit: 'm', reverse: true },
          series: [...layers.filter((l) => l.soil === 'clay').map((l, i) => ({ name: `Clay layer ${i + 1}`, x: [l.suTop / 1e3, l.suBottom / 1e3], y: [l.top, l.bottom] })), { name: 'Design line', x: [d.su0 / 1e3, (d.su0 + d.suGrad * d.depth) / 1e3], y: [0, d.depth], dash: true }] })] : [], notes: [] };
    }
    const t = parseCptTable(text.value, { depth: zU.value, qc: qU.value, fs: fU.value, u2: uU.value });
    const N = Number(String(nkt.value).replace(',', '.')), a = Number(String(area.value).replace(',', '.'));
    if (!(N >= 5 && N <= 30)) throw new Error('The cone factor Nkt must lie between 5 and 30 (10 to 20 is usual).');
    if (!(a > 0.3 && a <= 1)) throw new Error('The cone net area ratio must lie between 0.3 and 1 (0.75 to 0.85 is usual).');
    const d = cptDesignValues({ ...t, Nkt: N, a });
    const P = d.points, pc = (v) => `${Math.round(v * 100)} %`;
    const rows = [
      d.su0 != null && { k: KEY_LABEL['site.su0'][0], v: `${fmt(d.su0 / 1e3, 3)} kPa (${fmt(d.suRange.su0[0] / 1e3, 2)}–${fmt(d.suRange.su0[1] / 1e3, 2)} for Nkt 20–10)`, how: `su = (qt − σv0)/Nkt with Nkt = ${N} at each reading that behaves as clay or silt (Lunne, Robertson & Powell 1997), then a straight line by least squares through ${d.suFit.n} readings (fit quality R² = ${fmt(d.suFit.r2, 2)})${d.suFit.throughOrigin ? '; the free fit gave a negative mudline strength, so the line is taken through zero' : ''}. Nkt must be calibrated against laboratory strengths: the bracket shows what 10 and 20 would give.` },
      d.su0 != null && { k: KEY_LABEL['site.suGrad'][0], v: `${fmt(d.suGrad / 1e3, 3)} kPa/m (${fmt(d.suRange.suGrad[0] / 1e3, 2)}–${fmt(d.suRange.suGrad[1] / 1e3, 2)})`, how: 'Slope of the same line.' },
      d.phiDeg != null && { k: KEY_LABEL['site.phiDeg'][0], v: `${fmt(d.phiDeg, 3)}° (readings from ${fmt(d.phiRange[0], 3)}° to ${fmt(d.phiRange[2], 3)}°, median ${fmt(d.phiRange[1], 3)}°)`, how: `φ′ = 17.6° + 11° log₁₀[(qt/pa)/√(σ′v0/pa)] at each reading that behaves as sand (Kulhawy & Mayne 1990); the value passed on is the lower quartile, so three quarters of the readings are stronger. Relative density about ${pc(d.Dr)} (Jamiolkowski et al. 2003).` },
      { k: KEY_LABEL['site.gammaSub'][0], v: `${fmt(d.gammaSub / 1e3, 3)} kN/m³`, how: `Total unit weight from cone resistance and friction ratio (Robertson & Cabal 2010, about ±1 kN/m³) less the unit weight of sea water, averaged over the ${d.share.sand > 0.5 ? 'sand' : 'clay and silt'} readings, which make up most of the sounding.` },
      { k: 'Soil behaviour over the sounding', v: `clay ${pc(d.share.clay)}, silt mixtures ${pc(d.share.silt)}, sand ${pc(d.share.sand)} of ${fmt(d.depth, 3)} m`, how: 'Soil behaviour type index Ic of Robertson (2009): above 2.95 clay, 2.6 to 2.95 silt mixtures, below 2.6 sand. It describes how the soil responds to the cone, not its grading.' },
      { k: 'Median grain size', v: 'not derived', how: 'A cone test does not measure grain size: take it from samples.' }].filter(Boolean);
    const zs = P.z, coh = zs.map((_, i) => (P.group[i] === 'sand' ? NaN : P.su[i] / 1e3));
    const plots = [renderPlot({ type: 'scatter', title: 'Cone resistance and soil behaviour', subtitle: `${t.format}; ${t.z.length} readings`, x: { label: 'qt', unit: 'MPa' }, y: { label: 'Depth below mudline', unit: 'm', reverse: true }, series: ['clay', 'silt', 'sand'].map((g) => ({ name: g === 'silt' ? 'silt mixtures' : g, x: zs.map((_, i) => (P.group[i] === g ? P.qt[i] / 1e6 : NaN)).filter(Number.isFinite), y: zs.filter((_, i) => P.group[i] === g) })).filter((s) => s.x.length) })];
    if (d.su0 != null) plots.push(renderPlot({ type: 'scatter', title: 'Undrained shear strength with depth', x: { label: 'su', unit: 'kPa' }, y: { label: 'Depth below mudline', unit: 'm', reverse: true }, series: [{ name: `From the cone, Nkt = ${N}`, x: coh.filter(Number.isFinite), y: zs.filter((_, i) => Number.isFinite(coh[i])) }, { name: 'Design line', x: [d.su0 / 1e3, (d.su0 + d.suGrad * d.depth) / 1e3], y: [0, d.depth], line: true }] }));
    return { d, rows, plots, notes: [...t.warnings, t.u2 && !t.qcIsQt ? `Cone resistance corrected for pore pressure: qt = qc + u2 (1 − a), a = ${a}.` : t.qcIsQt ? 'The table holds the corrected cone resistance qt, used as given.' : 'No pore pressure in the table: the measured cone resistance is used uncorrected, which under-reads soft clay slightly.'], units: t.units };
  }

  function run() {
    clear(out);
    if (!text.value.trim()) return;
    let r;
    try { r = interpret(); } catch (e) { out.append(h('div.banner.fail', icon('warn'), h('div', h('b', 'This record could not be read. '), e.message))); return; }
    p.site.soilData = { kind: kind.value, text: text.value, units: r.units || { depth: zU.value, qc: qU.value, fs: fU.value, u2: uU.value }, Nkt: Number(nkt.value) || 15, a: Number(area.value) || 0.8, at: new Date().toISOString() };
    save();
    const values = Object.fromEntries([['site.su0', r.d.su0], ['site.suGrad', r.d.suGrad], ['site.phiDeg', r.d.phiDeg], ['site.gammaSub', r.d.gammaSub]].filter(([, v]) => Number.isFinite(v)));
    out.append(kv(r.rows), r.notes.length ? h('ul.plain', r.notes.map((n) => h('li', n))) : null, h('div.grid.c2', { style: { marginTop: '10px' } }, r.plots),
      h('div.row', { style: { marginTop: '10px' } }, h('button.btn.primary', { type: 'button', onclick: () => {
        // measured data replace the screening profile wherever they overlap
        publish(values, 'user');
        toast('Soil values from your site data are now used in every suite.'); redraw();
      } }, icon('check'), 'Use these values in every suite'), h('small.muted', `Sets ${Object.keys(values).map((k) => KEY_LABEL[k][0].toLowerCase()).join(', ')} as your own values. One sounding describes one location: where several exist, enter the cautious envelope in the table above.`)));
  }
  [zU, qU, fU, uU, kind].forEach((el) => el.addEventListener('change', run));
  [nkt, area].forEach((el) => el.addEventListener('change', run));
  text.addEventListener('change', run);
  const fold = h('details.fold', { open: !!saved.text }, h('summary', 'From real data: import a cone penetration test or a borehole'),
    h('div', hint, h('div.soilgrid', h('label', 'Kind of record', kind), zL, qL, fL, uL, h('label', 'Cone factor Nkt (10–20 usual)', nkt), h('label', 'Cone net area ratio a', area)), text,
      h('div.row', { style: { marginTop: '8px' } }, file, h('button.btn', { type: 'button', onclick: () => file.click() }, icon('upload'), 'Upload a file'), h('button.btn.primary', { type: 'button', onclick: run }, 'Interpret'),
        saved.text ? h('button.btn', { type: 'button', onclick: () => { delete p.site.soilData; save(); redraw(); } }, 'Remove the stored record') : null), out));
  if (saved.text) run();
  return fold;
}

/** Published ranges and the opt-in screening profile. */
function screeningFold(redraw) {
  const p = getProject(), sb = p.site.data?.seabed, obs = p.site.observation;
  const cls = obs?.cls || sb?.cls || null, prof = screeningProfile(cls, { carbonate: sb?.carbonate });
  const active = SCREENING_KEYS.filter((k) => p.bus[k]?.from === 'screening');
  const body = h('div');
  if (!prof.offered) body.append(h('p', prof.why));
  else {
    body.append(h('p', h('b', `${prof.title}. `), `Published typical ranges for the kind of seabed identified (${SEABED_CLASSES[cls].label.toLowerCase()}):`),
      table([{ key: 'k', label: 'Quantity', left: true }, { key: 'v', label: 'Published range', left: true }], prof.ranges.map(([k, v]) => ({ k, v }))),
      h('p', h('b', 'The screening profile takes: '), prof.chosen), h('p', prof.caution));
    const taken = Object.keys(prof.values).filter((k) => p.bus[k]?.from === 'user');
    body.append(h('div.row', h('button.btn', { type: 'button', onclick: () => {
      for (const k of SCREENING_KEYS) if (p.bus[k]?.from === 'screening') delete p.bus[k];
      publish(Object.fromEntries(Object.entries(prof.values).filter(([k]) => p.bus[k]?.from !== 'user')), 'screening');
      toast('Screening soil profile in use: every input and result that depends on it is marked.'); redraw();
    } }, 'Use a screening soil profile until site data arrive'),
    active.length ? h('button.btn', { type: 'button', onclick: () => { for (const k of active) delete p.bus[k]; save(); emit('bus', p.bus); toast('Screening soil profile removed.'); redraw(); } }, 'Remove the screening profile') : null),
    taken.length ? h('small.muted', `Kept as you entered them, not replaced by the profile: ${taken.map((k) => KEY_LABEL[k][0].toLowerCase()).join(', ')}.`) : null);
  }
  if (active.length) body.prepend(h('div.banner.warn.screening', icon('warn'), h('div', h('b', 'A screening soil profile is in use. '), `${active.map((k) => `${KEY_LABEL[k][0]} ${fmt(p.bus[k].value * KEY_LABEL[k][2], 3)} ${KEY_LABEL[k][1]}`).join('; ')} — ${SCREENING_FLAG}. Inputs fed by it are marked in every analysis, and every result computed on it is labelled a screening result.`)));
  body.append(h('p.why', h('b', 'Basis: '), prof.refs.join(' ')));
  return h('details.fold', { open: active.length > 0 }, h('summary', 'Before any data exists: published ranges and a screening soil profile'), body);
}

/** What investigation is normally expected for the seabed and the structure in hand. */
function plannerFold() {
  const p = getProject(), sb = p.site.data?.seabed, cls = p.site.observation?.cls || sb?.cls || null;
  const pick = h('select', { 'aria-label': 'Type of structure' }, Object.entries(STRUCTURES).map(([k, s]) => h('option', { value: k, selected: p.site.structure === k }, s.label)));
  const out = h('div');
  const draw = () => {
    p.site.structure = pick.value; save();
    const pl = investigationPlan(cls, pick.value, { carbonate: sb?.carbonate });
    clear(out).append(h('ul.plain', h('li', pl.first), h('li', h('b', 'Depth of investigation: '), pl.depth, '.'), h('li', h('b', 'Where: '), pl.where, '.'), pl.testing.map((x) => h('li', x))),
      h('p.why', h('b', 'Planning guidance, not a specification. Basis: '), pl.basis));
  };
  pick.addEventListener('change', draw); draw();
  return h('details.fold', h('summary', 'Planning the site investigation'), h('div', h('div.row', h('label', 'Structure'), pick, h('small.muted', cls ? `Seabed expected: ${SEABED_CLASSES[cls].label.toLowerCase()}.` : 'Seabed type unknown: the general programme is shown.')), out));
}

/** Cards for the site page: seabed sediment, and soil strength for the analyses. */
export function soilCards(redraw) {
  const p = getProject();
  if (!p.site.data?.seabed && !p.site.observation && p.site.lat == null) return [];
  const origin = (k) => { const b = p.bus[k]; return !b ? 'not set' : b.from === 'user' ? 'your own value' : b.from === 'screening' ? SCREENING_FLAG : b.from; };
  const strength = h('div.card', { style: { marginTop: '14px' } }, h('h3', 'Soil strength for the analyses'),
    h('p', 'Foundation, anchor and pipe–soil analyses need strength, friction angle and unit weight. A description of the seabed surface never becomes a strength here: the values come from your site data, or — clearly marked — from a screening profile you choose.'),
    table([{ key: 'k', label: 'Quantity', left: true }, { key: 'v', label: 'Value now' }, { key: 'u', label: 'Unit', left: true }, { key: 'o', label: 'Origin', left: true }],
      SCREENING_KEYS.map((k) => ({ k: KEY_LABEL[k][0], v: p.bus[k] ? p.bus[k].value * KEY_LABEL[k][2] : '—', u: KEY_LABEL[k][1], o: origin(k) }))),
    investigationFold(redraw), screeningFold(redraw), plannerFold());
  return [seabedCard(redraw), strength];
}
