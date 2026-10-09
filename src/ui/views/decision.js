// Decision support: one place that reads every run of the project and turns it into ranked actions,
// a cost and carbon picture from the quantities the suites produced, and live market context.
import { h, fmt, icon, table, clear } from '../dom.js';
import { renderPlot } from '../plot.js';
import { getProject } from '../../core/store.js';
import { suiteIndex, suiteModules, SHORT } from '../../core/registry.js';
import { commodity, fxRates, COMMODITIES, worldBank, WB_INDICATORS, gridCarbon } from '../../data/connectors.js';
import { adviceList } from '../results.js';
import { pager } from '../shell.js';

// Embodied-carbon factors, kg CO2e per unit. Mid-range open-literature values (ICE database v3 order of
// magnitude); they are screening defaults — replace with product EPDs in the Economics suite for a real LCA.
const CARBON = { steel_t: [1900, 't structural steel'], concrete_m3: [300, 'm³ concrete'], rock_t: [10, 't rock'], cable_km: [60000, 'km cable'], pipe_km: [250000, 'km line pipe'], copper_t: [3800, 't copper'], aluminium_t: [8600, 't aluminium'], polymer_t: [2500, 't polymer'] };
const FUEL_T_PER_DAY = 20; // generic offshore construction vessel; 3.206 t CO2 per t MGO (IMO)

export default async function decision(root) {
  const p = getProject();
  const idx = await suiteIndex();
  const mods = new Map();
  await Promise.all(idx.map(async (s) => { for (const m of await suiteModules(s.n)) mods.set(m.id, m); }));
  const runs = Object.entries(p.runs).map(([id, r]) => ({ id, ...r, mod: mods.get(id) })).filter((r) => r.mod);

  root.append(h('h1', 'Decision support'), h('p.lede', 'Every analysis you have run, read together: what fails, what governs, what to change first, and what the design costs in money and carbon.'));
  if (!runs.length) {
    root.append(h('div.empty', h('p', 'Run some analyses first — this page fills itself from their results.'), h('a.btn.primary', { href: '#/suite/1' }, 'Open the first suite')), pager({ href: '#/suite/17', label: '17. Economics' }, { href: '#/verification', label: 'Verification & validation' }));
    return { title: 'Decision support', crumbs: [{ href: '#/decision', label: 'Decision support' }] };
  }

  const failing = runs.filter((r) => r.status === 'fail').sort((a, b) => (b.maxUC || 0) - (a.maxUC || 0));
  const near = runs.filter((r) => r.status !== 'fail' && r.maxUC > 0.9);
  const reserve = runs.filter((r) => r.status !== 'fail' && r.maxUC > 0 && r.maxUC < 0.5);
  root.append(h('div.kpis',
    h('div.kpi', h('span', 'Analyses run'), h('b', runs.length), h('small.muted', `${new Set(runs.map((r) => r.mod.suite)).size} of 17 suites touched`)),
    h('div.kpi', { class: failing.length ? 'fail' : 'ok' }, h('span', 'Criteria not met'), h('b', failing.length), h('small.muted', failing.length ? 'action needed' : 'none')),
    h('div.kpi', h('span', 'Within 10 % of a limit'), h('b', near.length), h('small.muted', 'little margin')),
    h('div.kpi', h('span', 'Below 50 % utilisation'), h('b', reserve.length), h('small.muted', 'candidates to optimise'))));

  // ranked actions
  const actions = [];
  for (const r of failing) actions.push({ severity: 'fail', text: `${r.mod.title} (suite ${r.mod.suite}) fails: utilisation ${fmt(r.maxUC, 3)}${r.governing ? ` on “${r.governing}”` : ''}.`, action: h('span', (r.advice || []).filter((a) => a.action).slice(0, 2).map((a) => a.action).join(' ') || 'Open the module and use its Advisor tab to find the smallest design change that satisfies every criterion.', ' ', h('a', { href: `#/m/${r.id}` }, 'Open module →')) });
  for (const r of near) actions.push({ severity: 'warn', text: `${r.mod.title}: utilisation ${fmt(r.maxUC, 3)} — less than 10 % margin${r.governing ? ` on “${r.governing}”` : ''}.`, action: h('span', 'Check the sensitivity and uncertainty of this result before freezing the design. ', h('a', { href: `#/m/${r.id}` }, 'Open module →')) });
  for (const r of reserve.slice(0, 5)) actions.push({ severity: 'info', text: `${r.mod.title}: utilisation only ${fmt(r.maxUC, 3)}.`, action: h('span', 'There is reserve here: the Advisor tab shows how far each design variable can be reduced, which saves material, cost and carbon. ', h('a', { href: `#/m/${r.id}` }, 'Open module →')) });
  const notRun = idx.filter((s) => !runs.some((r) => r.mod.suite === s.n));
  if (notRun.length && notRun.length < 17) actions.push({ severity: 'info', text: `Suites not yet exercised: ${notRun.map((s) => s.n).join(', ')}.`, action: 'Interfaces between suites (loads, motions, reactions, quantities) are only closed when both sides have been run.' });
  if (!p.site.pulled) actions.push({ severity: 'warn', text: 'No site has been pulled for this project.', action: h('span', 'Environmental inputs are still module defaults. ', h('a', { href: '#/site' }, 'Choose the site →')) });
  root.append(h('h2', { style: { marginTop: '22px' } }, 'What to do next, in order'), adviceList(actions));

  // utilisation overview
  const withUC = runs.filter((r) => r.maxUC > 0).sort((a, b) => b.maxUC - a.maxUC).slice(0, 24);
  if (withUC.length) root.append(h('div', { style: { marginTop: '14px' } }, renderPlot({ type: 'bar', title: 'Highest utilisation per analysis', categories: withUC.map((r) => r.mod.title), series: [{ name: 'Utilisation', y: withUC.map((r) => r.maxUC), colors: withUC.map((r) => (r.maxUC > 1 ? 'var(--critical)' : r.maxUC > 0.9 ? 'var(--warning)' : 'var(--s1)')) }], y: { label: 'Utilisation (limit = 1)' }, limit: 1 })));

  // quantities → cost & carbon
  const q = {}, vessel = {};
  for (const r of runs) for (const [k, v] of Object.entries(r.quantities || {})) {
    if (k === 'vessel_days' && v && typeof v === 'object') for (const [vk, vv] of Object.entries(v)) vessel[vk] = (vessel[vk] || 0) + (Number(vv) || 0);
    else if (typeof v === 'number' && Number.isFinite(v)) q[k] = (q[k] || 0) + v;
  }
  const vesselDays = Object.values(vessel).reduce((s, v) => s + v, 0);
  const carbonRows = Object.entries(q).filter(([k]) => CARBON[k]).map(([k, v]) => ({ item: CARBON[k][1], qty: v, factor: CARBON[k][0], t: v * CARBON[k][0] / 1000 }));
  if (vesselDays) carbonRows.push({ item: 'vessel days (marine operations)', qty: vesselDays, factor: FUEL_T_PER_DAY * 3206, t: vesselDays * FUEL_T_PER_DAY * 3.206 });
  const totalC = carbonRows.reduce((s, r) => s + r.t, 0);
  const sus = h('div.card', { style: { marginTop: '14px' } }, h('h2', icon('leaf'), ' Material, marine-operations and carbon footprint'));
  if (carbonRows.length) {
    sus.append(table([{ key: 'item', label: 'Quantity from the analyses', left: true }, { key: 'qty', label: 'Amount' }, { key: 'factor', label: 'Factor [kg CO₂e / unit]' }, { key: 't', label: 'Embodied [t CO₂e]' }], carbonRows),
      h('p', h('b', `About ${fmt(totalC, 3)} t CO₂e embodied`), ' in the quantities produced so far.', p.bus['energy.aep'] ? ` Against an annual energy yield of ${fmt(p.bus['energy.aep'].value / 1e9, 3)} GWh that is ${fmt(totalC * 1e6 / (p.bus['energy.aep'].value / 1e3 * 25), 3)} g CO₂e per kWh over 25 years.` : ''),
      h('p.muted', 'Factors are mid-range open-literature screening values, not product declarations. The Economics suite carries the editable inventory, transport, operation, end-of-life credits and carbon pricing.'));
    const big = carbonRows.slice().sort((a, b) => b.t - a.t)[0];
    sus.append(adviceList([{ severity: 'info', text: `${big.item} carries ${fmt(100 * big.t / totalC, 2)} % of the embodied carbon.`, action: 'Reducing it — through the reserve identified above, higher-strength material, recycled content or a shorter installation campaign — is the most effective sustainability lever in this project.' }]));
  } else sus.append(h('p.muted', 'No material or vessel quantities have been produced yet. Structural, pipeline, cable and installation modules report them when run.'));
  root.append(sus);

  // economics snapshot from the bus + live market context
  const eco = h('div.card', { style: { marginTop: '14px' } }, h('h2', icon('coin'), ' Economics'));
  const b = (k) => p.bus[k]?.value;
  if (b('econ.capex') != null || b('econ.npv') != null || b('econ.lcoe') != null) {
    eco.append(h('div.kpis', [['econ.capex', 'Capital expenditure', p.currency], ['econ.opexAnnual', 'Operating cost per year', p.currency], ['econ.npv', 'Net present value', p.currency], ['econ.lcoe', 'Levelised cost of energy', `${p.currency}/MWh`]].filter(([k]) => b(k) != null).map(([k, l, u]) => h('div.kpi', { class: k === 'econ.npv' ? (b(k) >= 0 ? 'ok' : 'fail') : '' }, h('span', l), h('b', fmt(b(k))), h('em', u)))));
    if (b('econ.npv') < 0) eco.append(adviceList([{ severity: 'fail', text: 'The project destroys value at the current assumptions (NPV below zero).', action: 'Use the sensitivity and break-even modules of the Economics suite to find the price, cost or yield at which it turns positive.' }]));
  } else eco.append(h('p.muted', 'Run the cash-flow or LCOE modules of the Economics suite to see value metrics here. They pick up quantities, energy yield and installation days from the other suites.'), h('a.btn', { href: '#/suite/17' }, 'Open Economics'));
  const market = h('div', { style: { marginTop: '12px' } }, h('small.muted', 'Loading live market context…'));
  eco.append(market);
  root.append(eco);
  (async () => {
    const out = h('div.grid.c2');
    const notes = [];
    for (const key of ['brent', 'ironOre', 'copper', 'natgasEU']) {
      try {
        const r = await commodity(key);
        const rows = r.data.rows.slice(-120);
        const last = rows[rows.length - 1], yearAgo = rows[Math.max(0, rows.length - 13)];
        const chg = (last.value - yearAgo.value) / yearAgo.value;
        out.append(renderPlot({ type: 'line', title: `${COMMODITIES[key].label} — ${fmt(last.value, 4)} ${COMMODITIES[key].unit} (${last.period})`, subtitle: `IMF Primary Commodity Prices · ${chg >= 0 ? '+' : ''}${fmt(100 * chg, 2)} % over 12 months${r.stale ? ' · cached' : ''}`, x: { label: 'Months (last 10 years)' }, y: { label: COMMODITIES[key].unit, zero: false }, series: [{ name: COMMODITIES[key].label, x: rows.map((_, i) => i - rows.length + 1), y: rows.map((x) => x.value) }] }));
        if (key === 'ironOre' && Math.abs(chg) > 0.15) notes.push({ severity: chg > 0 ? 'warn' : 'info', text: `Iron ore is ${chg > 0 ? 'up' : 'down'} ${fmt(100 * Math.abs(chg), 2)} % in a year.`, action: `Steel-heavy cost lines (${q.steel_t ? `${fmt(q.steel_t, 3)} t in this project` : 'jackets, piles, line pipe'}) should be re-priced before sanction; consider an escalation clause or earlier procurement.` });
        if (key === 'copper' && Math.abs(chg) > 0.15) notes.push({ severity: chg > 0 ? 'warn' : 'info', text: `Copper is ${chg > 0 ? 'up' : 'down'} ${fmt(100 * Math.abs(chg), 2)} % in a year.`, action: 'Re-check the copper-versus-aluminium conductor choice in the cable suite; aluminium conductors need a larger cross-section but decouple the cable price from copper.' });
        if (key === 'brent' && Math.abs(chg) > 0.15) notes.push({ severity: 'info', text: `Brent is ${chg > 0 ? 'up' : 'down'} ${fmt(100 * Math.abs(chg), 2)} % in a year.`, action: 'Vessel day rates and hydrocarbon revenue both move with it: rerun the price scenarios in the Economics suite.' });
      } catch { /* a missing series leaves the others intact */ }
    }
    const extra = [];
    try { const fx = await fxRates(p.currency === 'EUR' ? 'EUR' : 'USD'); extra.push(`Exchange rates (ECB, ${fx.data.date}): 1 ${fx.data.base} = ${['EUR', 'USD', 'GBP', 'NOK', 'NGN', 'CNY', 'BRL'].filter((c) => fx.data.rates[c]).map((c) => `${fmt(fx.data.rates[c], 4)} ${c}`).join(' · ')}.`); } catch { /* optional */ }
    const iso2 = p.site.data?.place?.iso2;
    if (iso2) { try { const inf = await worldBank(iso2, WB_INDICATORS.inflation); const l = inf.data[inf.data.length - 1]; if (l) extra.push(`Consumer-price inflation in ${p.site.data.place.country}: ${fmt(l.value, 3)} % (${l.year}, World Bank).`); } catch { /* optional */ } }
    try { const c = await gridCarbon(); const v = c.data.data?.[0]?.intensity; if (v) extra.push(`Great Britain grid carbon intensity right now: ${v.actual ?? v.forecast} g CO₂/kWh (${v.index}).`); } catch { /* optional */ }
    clear(market).append(h('h3', 'Live market context'), adviceList(notes), out, h('ul.prose', extra.map((e) => h('li', e))));
  })();

  root.append(h('h2', { style: { marginTop: '22px' } }, 'All runs'), table([{ key: 'suite', label: 'Suite', left: true }, { key: 'title', label: 'Analysis', left: true }, { key: 'status', label: 'Status', left: true }, { key: 'uc', label: 'Highest utilisation' }, { key: 'gov', label: 'Governing check', left: true }, { key: 'at', label: 'Run', left: true }],
    runs.sort((a, b) => a.mod.suite - b.mod.suite).map((r) => ({ suite: `${r.mod.suite}. ${SHORT[r.mod.suite]}`, title: h('a', { href: `#/m/${r.id}` }, r.mod.title), status: h('span', { class: `chip ${r.status === 'fail' ? 'fail' : 'ok'}` }, r.status === 'fail' ? '✗ not met' : '✓ met'), uc: r.maxUC || null, gov: r.governing || '', at: new Date(r.at).toLocaleString() }))),
  pager({ href: '#/suite/17', label: '17. Economics' }, { href: '#/verification', label: 'Verification & validation' }));
  return { title: 'Decision support', crumbs: [{ href: '#/decision', label: 'Decision support' }] };
}
