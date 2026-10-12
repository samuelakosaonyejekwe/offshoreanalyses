// Mooring line strength: component-by-component capacity along the line (chain, wire, fibre rope, connectors,
// fairlead / chain stopper) after corrosion and wear, against intact and damaged design tensions.
import { chainMBL, corrodedDiameter } from '../../engine/lines/lineprops.js';
import { checkApi2sk, checkDnvE301 } from '../../engine/lines/mooringsystem.js';
import { need, run, round, fmt, adv, check, worst, sev, linspace, componentProps, steelCo2Text, CONFIRM, LIMIT_PROPS } from './_common.js';

const COMPONENTS = [
  { name: 'Anchor shackle', type: 'connector', d: 0, grade: '', cap: 10500, rate: 0.2 },
  { name: 'Bottom chain', type: 'chain-studless', d: 95, grade: 'R4', cap: 0, rate: 0.3 },
  { name: 'Chain–wire socket', type: 'connector', d: 0, grade: '', cap: 10000, rate: 0.1 },
  { name: 'Spiral-strand wire', type: 'wire-spiral', d: 100, grade: '', cap: 0, rate: 0 },
  { name: 'Top chain', type: 'chain-studless', d: 95, grade: 'R4', cap: 0, rate: 0.4 },
  { name: 'Chain stopper / fairlead', type: 'connector', d: 0, grade: '', cap: 11000, rate: 0.1 },
];
const GRADE_C = { R3: 0.0223, R3S: 0.0249, R4: 0.0274, R4S: 0.0304, R5: 0.0320 };

export default {
  id: 's10-line-strength',
  suite: 10,
  title: 'Mooring line strength and failure',
  summary: 'Finds the weakest component of a mooring line at end of life and checks it against the intact and one-line-broken design tensions in both API and DNV formats.',
  covers: [7],
  fidelity: 'analytical',
  standards: ['DNV-OS-E301 (2018) Ch.2 Sec.2 [4] (ULS, ALS) and Sec.4 (corrosion allowance)', 'API RP 2SK (3rd ed., 2005) §6.3', 'IACS UR W22 / DNV-OS-E302 (chain break load)', 'ISO 19901-7 (2013) §10'],
  models: ['Catalogue break load of offshore chain BL = c·d²·(44 − 0.08d)', 'Wire and fibre rope break-load estimators', 'Corrosion / wear reduction of diameter over the service life', 'API safety-factor format', 'DNV partial-factor format with characteristic strength 0.95·MBS', 'Weakest-link (series system) line strength'],
  limits: `Static strength of a series system of components under axial tension. Connector and fairlead capacities are the certified values you enter; local bending at fairleads, out-of-plane bending, torsion and brittle fracture are not assessed. ${LIMIT_PROPS} Safety factors are overridable defaults.`,
  publishes: ['mooring.MBLendOfLife'],
  searchTerms: 'mooring line strength minimum breaking load corrosion allowance safety factor',
  inputs: [
    { key: 'components', label: 'Line components (anchor → fairlead)', type: 'table', group: 'Line components', role: 'material',
      columns: [{ key: 'name', label: 'Component', type: 'text' }, { key: 'type', label: 'Type', type: 'text' }, { key: 'd', label: 'Diameter', unit: 'mm', type: 'number' }, { key: 'grade', label: 'Chain grade', type: 'text' },
        { key: 'cap', label: 'Certified break load (0 = estimate)', unit: 'kN', type: 'number' }, { key: 'rate', label: 'Corrosion / wear on diameter', unit: 'mm/yr', type: 'number' }],
      default: COMPONENTS,
      help: 'Type: chain-studless, chain-studlink, wire-six, wire-spiral, polyester, hmpe, nylon, or connector (shackles, links, sockets, chain stoppers, fairleads — these need a certified break load). Corrosion rates per DNV-OS-E301: about 0.2–0.4 mm/yr at the splash zone and seabed, less mid-water; sheathed wire and fibre rope 0.' },
    { key: 'life', label: 'Design service life', unit: 'years', type: 'number', default: 25, min: 0, max: 60, group: 'Line components', role: 'criteria', bind: 'life.designYears' },
    { key: 'TmeanI', label: 'Mean tension, intact', unit: 'kN', type: 'number', default: 2220, min: 0, max: 60000, group: 'Design tensions', role: 'load', help: 'Characteristic mean tension of the governing line in the design storm (tension-and-offset module).' },
    { key: 'TmaxI', label: 'Maximum tension, intact', unit: 'kN', type: 'number', default: 3100, min: 1, max: 60000, group: 'Design tensions', role: 'load', bind: 'mooring.Tmax', uncertain: { dist: 'gumbel', cov: 0.12 } },
    { key: 'TmeanD', label: 'Mean tension, one line broken', unit: 'kN', type: 'number', default: 3250, min: 0, max: 60000, group: 'Design tensions', role: 'load', help: 'From the line-break module: mean tension of the most loaded remaining line.' },
    { key: 'TmaxD', label: 'Maximum tension, one line broken', unit: 'kN', type: 'number', default: 4400, min: 1, max: 60000, group: 'Design tensions', role: 'load', bind: 'mooring.TmaxDamaged' },
    { key: 'sfIntact', label: 'API safety factor, intact', unit: '-', type: 'number', default: 1.67, min: 1, max: 4, group: 'Acceptance criteria', role: 'criteria', help: `Dynamic analysis (2.0 quasi-static). ${CONFIRM}` },
    { key: 'sfDamaged', label: 'API safety factor, one line broken', unit: '-', type: 'number', default: 1.25, min: 1, max: 3, group: 'Acceptance criteria', role: 'criteria', help: `Dynamic analysis (1.43 quasi-static). ${CONFIRM}` },
    { key: 'gmU', label: 'DNV γmean, ULS', unit: '-', type: 'number', default: 1.1, min: 0.8, max: 3, group: 'Acceptance criteria', role: 'criteria', help: `Consequence class 1: 1.10; class 2: 1.40. ${CONFIRM}` },
    { key: 'gdU', label: 'DNV γdyn, ULS', unit: '-', type: 'number', default: 1.5, min: 0.8, max: 3.5, group: 'Acceptance criteria', role: 'criteria', help: `Consequence class 1: 1.50; class 2: 2.10. ${CONFIRM}` },
    { key: 'gmA', label: 'DNV γmean, ALS', unit: '-', type: 'number', default: 1.0, min: 0.8, max: 3, group: 'Acceptance criteria', role: 'criteria', help: `1.00 for both consequence classes. ${CONFIRM}` },
    { key: 'gdA', label: 'DNV γdyn, ALS', unit: '-', type: 'number', default: 1.1, min: 0.8, max: 3.5, group: 'Acceptance criteria', role: 'criteria', help: `Consequence class 1: 1.10; class 2: 1.25. ${CONFIRM}` },
    { key: 'strengthFactor', label: 'Characteristic strength / MBL', unit: '-', type: 'number', default: 0.95, min: 0.7, max: 1, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'startCondition', label: 'Present condition', type: 'select', default: 'new', options: [{ value: 'new', label: 'New line' }, { value: 'used', label: 'In service (years entered below)' }], group: 'Initial condition', role: 'ic',
      help: 'Initial condition of the line components: new, or already corroded by the years in service.' },
    { key: 'yearsUsed', label: 'Years already in service', unit: 'years', type: 'number', default: 0, min: 0, max: 60, group: 'Initial condition', role: 'ic', showIf: { key: 'startCondition', in: ['used'] } },
  ],
  compute(inp) {
    need(Array.isArray(inp.components) && inp.components.length > 0, 'The component table is empty: add at least one component.');
    need(inp.TmaxI >= inp.TmeanI && inp.TmaxD >= inp.TmeanD, 'A maximum tension is smaller than the corresponding mean tension: check the design tensions.');
    const used = inp.startCondition === 'used' ? inp.yearsUsed : 0;
    const mblAt = (r, years) => {
      const type = String(r.type || '').trim().toLowerCase();
      if (type === 'connector') {
        need(r.cap > 0, `Component '${r.name}': connectors, shackles and fairleads need a certified break load.`);
        // capacity of a forged connector reduced with the section lost by corrosion on an equivalent 100 mm section
        return r.cap * 1e3 * Math.max(0.05, 1 - 2 * (r.rate || 0) * years / 100);
      }
      need(r.d > 0, `Component '${r.name}': diameter must be positive.`);
      const d = r.d / 1e3;
      const dc = (r.rate || 0) > 0 && years > 0 ? run(() => corrodedDiameter({ d, rate: r.rate / 1e3, years }), `Corrosion of '${r.name}'`) : d;
      need(dc > 0.3 * d, `Component '${r.name}' loses more than 70 % of its diameter over ${years} years: reduce the corrosion rate or the life.`);
      if (type.startsWith('chain')) {
        const g = String(r.grade || 'R4').toUpperCase();
        need(GRADE_C[g], `Component '${r.name}': chain grade '${r.grade}' is not R3, R3S, R4, R4S or R5.`);
        const m = chainMBL({ d: dc, grade: g });
        return r.cap > 0 ? r.cap * 1e3 * m / chainMBL({ d, grade: g }) : m;
      }
      const p = componentProps({ type, d });
      return (r.cap > 0 ? r.cap * 1e3 : p.MBL) * (dc / d) ** 2;
    };
    const rows = inp.components.map((r) => {
      const now = mblAt(r, used), eol = mblAt(r, inp.life), fresh = mblAt(r, 0);
      const apiI = checkApi2sk({ tension: inp.TmaxI * 1e3, MBL: eol, condition: 'intact', requiredSF: inp.sfIntact });
      const apiD = checkApi2sk({ tension: inp.TmaxD * 1e3, MBL: eol, condition: 'damaged', requiredSF: inp.sfDamaged });
      const uls = checkDnvE301({ Tmean: inp.TmeanI * 1e3, Tdyn: (inp.TmaxI - inp.TmeanI) * 1e3, MBL: eol, limitState: 'ULS', strengthFactor: inp.strengthFactor, gammaMean: inp.gmU, gammaDyn: inp.gdU });
      const als = checkDnvE301({ Tmean: inp.TmeanD * 1e3, Tdyn: (inp.TmaxD - inp.TmeanD) * 1e3, MBL: eol, limitState: 'ALS', strengthFactor: inp.strengthFactor, gammaMean: inp.gmA, gammaDyn: inp.gdA });
      return { r, fresh, now, eol, apiI, apiD, uls, als, uc: Math.max(apiI.uc, apiD.uc, uls.uc, als.uc) };
    });
    const weak = rows.reduce((a, x) => (x.eol < a.eol ? x : a), rows[0]);
    const weakNew = rows.reduce((a, x) => (x.fresh < a.fresh ? x : a), rows[0]);
    const body = rows.filter((x) => String(x.r.type).toLowerCase() !== 'connector');
    const bodyMin = body.length ? Math.min(...body.map((x) => x.fresh)) : weakNew.fresh;
    const conn = rows.filter((x) => String(x.r.type).toLowerCase() === 'connector');
    const connMin = conn.length ? Math.min(...conn.map((x) => x.fresh)) : bodyMin;
    const checks = [
      check(`API intact: ${weak.r.name}`, inp.TmaxI * inp.sfIntact, weak.eol / 1e3, 'kN', 'API RP 2SK (3rd ed.) §6.3',
        'Guards against breaking the weakest component in the intact design storm at end of life: maximum tension times the safety factor against the corroded break load.'),
      check(`API one line broken: ${weak.r.name}`, inp.TmaxD * inp.sfDamaged, weak.eol / 1e3, 'kN', 'API RP 2SK (3rd ed.) §6.3 (damaged condition)',
        'Guards against a second line parting after one has failed; a second failure usually means loss of position.'),
      check(`DNV ULS: ${weak.r.name}`, weak.uls.Td / 1e3, weak.uls.Sc / 1e3, 'kN', 'DNV-OS-E301 Ch.2 Sec.2 [4.2]',
        'Ultimate limit state in partial-factor format: factored mean plus factored dynamic tension against the characteristic strength of the weakest component.'),
      check(`DNV ALS: ${weak.r.name}`, weak.als.Td / 1e3, weak.als.Sc / 1e3, 'kN', 'DNV-OS-E301 Ch.2 Sec.2 [4.3]',
        'Accidental limit state with one line missing: the remaining lines must survive the same environment with reduced factors.'),
      check('Connectors at least as strong as the line body', bodyMin / 1e3, connMin / 1e3, 'kN', 'DNV-OS-E301 Ch.2 Sec.4 (accessories: break load not less than the line); API RP 2SK §6.3',
        'Guards against a shackle, link or stopper being the hidden weak link: accessories should be certified to at least the break load of the chain or rope they join.'),
    ];
    const g = worst(checks);
    const yrs = linspace(0, Math.max(inp.life, 1) * 1.4, 15);
    const allow = Math.max(inp.TmaxI * inp.sfIntact, inp.TmaxD * inp.sfDamaged);
    // life at which the weakest component's capacity falls to the governing API demand
    let lifeToLimit = yrs[yrs.length - 1] * 3;
    for (let y = 0; y <= 150; y += 0.5) { let m; try { m = Math.min(...inp.components.map((r) => mblAt(r, y))); } catch (e) { lifeToLimit = y; break; } if (m / 1e3 < allow) { lifeToLimit = y; break; } }
    const chains = rows.filter((x) => String(x.r.type).toLowerCase().startsWith('chain'));
    return {
      kpis: [
        { key: 'mblNew', label: 'Line strength when new (weakest component)', value: round(weakNew.fresh / 1e3, 5), unit: 'kN', explain: 'Break load of the weakest component before any corrosion; the catalogue strength of the line as a series system.' },
        { key: 'mblEol', label: 'Line strength at end of life', value: round(weak.eol / 1e3, 5), unit: 'kN', explain: 'Break load of the weakest component after corrosion and wear over the design life; the strength the checks must use.' },
        { key: 'retention', label: 'Strength retained at end of life', value: round(100 * weak.eol / weakNew.fresh, 4), unit: '%', explain: 'End-of-life strength as a share of new strength. Chain typically retains 80–90 % after 20–25 years with 0.2–0.4 mm/yr corrosion.' },
        { key: 'sfIntact', label: 'Safety factor, intact', value: round(weak.eol / (inp.TmaxI * 1e3), 4), unit: '-', limit: inp.sfIntact, sense: 'min', explain: 'End-of-life break load over the intact maximum tension.' },
        { key: 'sfDamaged', label: 'Safety factor, one line broken', value: round(weak.eol / (inp.TmaxD * 1e3), 4), unit: '-', limit: inp.sfDamaged, sense: 'min', explain: 'End-of-life break load over the maximum tension with one line missing.' },
        { key: 'ucUls', label: 'DNV ULS utilisation', value: round(weak.uls.uc, 4), unit: '-', limit: 1, sense: 'max', explain: 'Factored intact tension over characteristic strength of the weakest component.' },
        { key: 'ucAls', label: 'DNV ALS utilisation', value: round(weak.als.uc, 4), unit: '-', limit: 1, sense: 'max', explain: 'Factored damaged-condition tension over characteristic strength of the weakest component.' },
        { key: 'utilMbl', label: 'Intact maximum tension / end-of-life MBL', value: round(100 * inp.TmaxI * 1e3 / weak.eol, 4), unit: '%', explain: 'Line utilisation against minimum breaking strength (60 % is the API intact limit for dynamic analysis).' },
        { key: 'lifeToLimit', label: 'Service life until strength governs', value: round(lifeToLimit, 4), unit: 'years', limit: inp.life, sense: 'min', explain: 'Years of corrosion after which the weakest component no longer meets the API safety factors; must exceed the design life.' },
      ],
      checks,
      plots: [
        { type: 'bar', title: 'Component break load: new, now, end of life', categories: rows.map((x) => x.r.name),
          series: [{ name: 'New (kN)', y: rows.map((x) => round(x.fresh / 1e3)) }, { name: 'Now (kN)', y: rows.map((x) => round(x.now / 1e3)) }, { name: 'End of life (kN)', y: rows.map((x) => round(x.eol / 1e3)) }, { name: 'Required (kN)', y: rows.map(() => round(allow)) }],
          explain: 'Strength of each component along the line at three ages against the larger of the factored intact and damaged tensions. The shortest end-of-life bar is the weak link.' },
        { type: 'line', title: 'Strength of the line over time', x: { label: 'Years in service', unit: 'years' }, y: { label: 'Break load', unit: 'kN' },
          series: [...rows.filter((x) => (x.r.rate || 0) > 0).map((x) => ({ name: x.r.name, x: yrs.map((y) => round(y)), y: yrs.map((y) => { try { return round(mblAt(x.r, y) / 1e3); } catch (e) { return 0; } }) })),
            { name: 'Required', x: [0, yrs[yrs.length - 1]], y: [allow, allow], dash: true }, { name: 'Design life', x: [inp.life, inp.life], y: [0, weakNew.fresh / 1e3], dash: true }],
          explain: 'Loss of break load of the corroding components with time. A curve crossing the required line before the design-life marker means that component must be upsized or replaced during the life.' },
        { type: 'bar', title: 'Utilisation by component and check', categories: rows.map((x) => x.r.name),
          series: [{ name: 'API intact', y: rows.map((x) => round(x.apiI.uc, 3)) }, { name: 'API damaged', y: rows.map((x) => round(x.apiD.uc, 3)) }, { name: 'DNV ULS', y: rows.map((x) => round(x.uls.uc, 3)) }, { name: 'DNV ALS', y: rows.map((x) => round(x.als.uc, 3)) }],
          explain: 'All four strength checks for every component at end of life. Anything above 1.0 fails; compare formats to see which code governs for this line.' },
      ],
      tables: [
        { title: 'Component strength and safety factors at end of life', columns: [{ key: 'name', label: 'Component' }, { key: 'type', label: 'Type' }, { key: 'fresh', label: 'MBL new', unit: 'kN' }, { key: 'eol', label: 'MBL end of life', unit: 'kN' }, { key: 'sfI', label: 'SF intact' }, { key: 'sfD', label: 'SF damaged' }, { key: 'uls', label: 'DNV ULS' }, { key: 'als', label: 'DNV ALS' }],
          rows: rows.map((x) => ({ name: x.r.name, type: x.r.type, fresh: round(x.fresh / 1e3), eol: round(x.eol / 1e3), sfI: round(x.apiI.safetyFactor, 3), sfD: round(x.apiD.safetyFactor, 3), uls: round(x.uls.uc, 3), als: round(x.als.uc, 3) })),
          explain: 'Per-component break loads and results of each check. The row with the lowest safety factors is the critical component to report and to target in inspections.' },
      ],
      notes: [`Weakest component at end of life: ${weak.r.name} (${fmt(weak.eol / 1e3)} kN, ${fmt(100 * weak.eol / weak.fresh)} % of new).`,
        `Governing check: ${g.name} with utilisation ${fmt(g.uc)}.`,
        'Chain break loads follow the IACS formula; rope values are indicative estimators unless a certified load is entered.'],
      advice: [
        adv('engineering', sev(g.uc), `Line strength falls from ${fmt(weakNew.fresh / 1e3)} kN new to ${fmt(weak.eol / 1e3)} kN after ${inp.life} years; safety factors ${fmt(weak.eol / (inp.TmaxI * 1e3))} intact and ${fmt(weak.eol / (inp.TmaxD * 1e3))} damaged. Governing: ${g.name} (${fmt(g.uc)}).`,
          g.uc > 1 ? (g.name.startsWith('Connectors') ? `Specify accessories certified to at least ${fmt(bodyMin / 1e3)} kN (the line-body break load); the weakest is now ${fmt(connMin / 1e3)} kN.` : `Raise the end-of-life strength of ${weak.r.name} by ${fmt(100 * (g.uc - 1))} %: ${String(weak.r.type).startsWith('chain') ? `about ${fmt(weak.r.d * Math.sqrt(g.uc))} mm diameter, or one grade up` : 'a larger size or higher-rated item'}, or add corrosion allowance only where the rate is high (splash zone, touchdown).`)
            : `Strength is adequate to ${fmt(lifeToLimit)} years of corrosion. Measure chain diameter at the splash zone and touchdown at each class renewal and compare with the assumed ${fmt(Math.max(...inp.components.map((r) => r.rate || 0)))} mm/yr.`,
          'DNV-OS-E301 (2018) Ch.2 Sec.2 [4] and Sec.4 [2] (corrosion allowance); API RP 2SK (3rd ed., 2005) §6.3; IACS UR W22. Confirm the factors against the governing edition.'),
        adv('economic', 'info', `Corrosion allowance costs steel: ${chains.length ? `each millimetre added to a ${fmt(chains[0].r.d)} mm chain adds ${fmt(100 * 2 / chains[0].r.d)} % mass` : 'rope and connector sizes step in catalogue increments'}; the governing utilisation is ${fmt(g.uc)}.`,
          g.uc < 0.8 ? `There is ${fmt(100 * (1 - g.uc))} % reserve at end of life: a lower grade or about ${fmt(100 * (1 - Math.sqrt(g.uc / 0.95)))} % smaller diameter would still pass strength (check fatigue before reducing).` : 'Replace-in-service of the top chain (accessible from the platform) is usually cheaper than oversizing the whole line for splash-zone corrosion; size the submerged segments for their own, lower, rate.',
          'Whole-life costing per ISO 15686-5; in-service inspection per API RP 2I (2008) and class renewal surveys.',
          `${steelCo2Text(chains.length ? 19.9e3 * (chains[0].r.d / 1e3) ** 2 * 100 / 1e3 : 10)} (per 100 m of the first chain segment). Segment-specific corrosion allowance avoids carrying unnecessary steel along the whole line.`),
      ],
      publish: { 'mooring.MBLendOfLife': weak.eol },
      quantities: chains.length ? { steel_t: round(chains.reduce((a, x) => a + 19.9e3 * (x.r.d / 1e3) ** 2 * 100, 0) / 1e3) } : {},
    };
  },
  verify: [
    { // IACS UR W22: BL = c·d²·(44 − 0.08·d) kN, d in mm. R4 (c = 0.0274), d = 120 mm: 0.0274·14400·34.4 = 13 573 kN.
      // After 25 years at 0.4 mm/yr: d = 110 mm → 0.0274·12100·35.2 = 11 670 kN (86.0 % retained).
      name: 'Chain break load new and corroded (R4, 120 mm)',
      source: 'IACS UR W22 Table of break loads; DNV-OS-E301 Ch.2 Sec.4 corrosion allowance (hand calculation)',
      inputs: { components: [{ name: 'chain', type: 'chain-studless', d: 120, grade: 'R4', cap: 0, rate: 0.4 }], life: 25, TmaxI: 5000, TmeanI: 3000, TmaxD: 6000, TmeanD: 4000 },
      expect: [{ kpi: 'mblNew', value: 13573, tol: 0.001 }, { kpi: 'mblEol', value: 11670, tol: 0.001 }, { kpi: 'sfIntact', value: 2.334, tol: 0.001 }] },
  ],
};
