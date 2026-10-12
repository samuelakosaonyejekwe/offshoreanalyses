// Crush, collapse and external pressure: hydrostatic collapse of depressurised tubes, external-pressure wall stress,
// and crushing of tubes by tensioner squeeze, chute line load and the radial pressure of the armour under tension.
import { tubeCollapsePressure, lameStresses } from '../../engine/xsection/tubes.js';
import { umbilicalCrushCheck, umbilicalLoadSharing } from '../../engine/xsection/umbilical.js';
import { sectionInputs, buildSection, adv, check, worst, sev, fmt, round, run, need, linspace, G, CO2_STEEL, LIMIT_RIGID, STD_UMB, CONFIRM, vesselCo2Text } from './_common.js';

export default {
  id: 's12-crush-collapse', suite: 12, title: 'Crush, collapse and external pressure',
  summary: 'Collapse utilisation of the tubes under external hydrostatic pressure and crush utilisation under tensioner, chute and armour contact loads.',
  covers: [11], fidelity: 'analytical',
  standards: ['DNV-ST-F101 (2021) §5.4.4 Eq. 5.10–5.13 (collapse with ovality)', `${STD_UMB} §7 (crush and collapse of tubes, installation loads)`, "Roark's Formulas for Stress and Strain, Table 9.2 (ring under opposed radial loads)"],
  models: ['External pressure buckling: elastic instability 2E(t/D)³/(1−ν²) combined with plastic collapse and ovality (Haagsma & Schaap cubic)', 'Lamé thick-wall stresses under external pressure', 'Crush resistance: ring bending under opposed line loads', 'Helical contact pressure from tension (layer interaction)'],
  limits: `${LIMIT_RIGID} Crush is an elastic first-yield ring criterion with a single load-spreading factor through fillers and sheath (conservative: ductile tubes tolerate local yielding); pad geometry, creep of polymer layers and nonlinear contact are not modelled. Qualification crush tests govern over this calculation.`,
  publishes: ['umbilical.allowSqueeze', 'umbilical.collapseUtilisation', 'umbilical.minChuteRadius'],
  searchTerms: 'umbilical steel tube collapse external pressure crush tensioner',
  inputs: [
    ...sectionInputs(),
    { key: 'depth', label: 'Water depth', unit: 'm', type: 'number', default: 350, min: 1, max: 4000, group: 'External pressure', role: 'bc', bind: 'site.depth', uncertain: { dist: 'normal', cov: 0.01 } },
    { key: 'pMin', label: 'Minimum internal pressure in the tubes', unit: 'bar', type: 'number', default: 0, min: 0, max: 2000, group: 'External pressure', role: 'bc', help: 'Absolute minimum during installation or after loss of pressure (0 = empty / vented tube: the collapse design case).' },
    { key: 'ovality', label: 'Tube ovality (Dmax − Dmin)/D', unit: '%', type: 'number', default: 0.5, min: 0.1, max: 5, group: 'External pressure', role: 'geometry', uncertain: { dist: 'lognormal', cov: 0.3 }, help: 'After manufacture and reeling; the code minimum of 0.5 % is enforced.' },
    { key: 'alphaFab', label: 'Fabrication factor', unit: '-', type: 'number', default: 0.93, min: 0.8, max: 1, group: 'External pressure', role: 'material', help: 'Seamless 1.0; seam-welded tube 0.93 (DNV-ST-F101 Table 5-5).' },
    { key: 'gammaM', label: 'Material resistance factor γm', unit: '-', type: 'number', default: 1.15, min: 1, max: 1.5, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'gammaSC', label: 'Safety-class resistance factor γSC', unit: '-', type: 'number', default: 1.14, min: 1, max: 1.5, group: 'Acceptance criteria', role: 'criteria', help: 'DNV-ST-F101 Table 5-2: low 1.04, medium 1.14, high 1.26.' },
    { key: 'T', label: 'Tension at the tensioner / chute', unit: 'kN', type: 'number', default: 120, min: 0, max: 5000, group: 'Crush loads', role: 'bc', bind: 'umbilical.layTension', uncertain: { dist: 'lognormal', cov: 0.1 } },
    { key: 'nTracks', label: 'Tensioner tracks', type: 'select', default: 4, options: [{ value: 2, label: '2-track' }, { value: 3, label: '3-track' }, { value: 4, label: '4-track' }], group: 'Crush loads', role: 'bc' },
    { key: 'squeeze', label: 'Squeeze load per track', unit: 'kN/m', type: 'number', default: 40, min: 0.5, max: 500, group: 'Crush loads', role: 'bc', design: true },
    { key: 'chuteR', label: 'Chute / sheave radius', unit: 'm', type: 'number', default: 9, min: 0.5, max: 30, group: 'Crush loads', role: 'bc', design: true },
    { key: 'spread', label: 'Elements sharing one pad load', unit: '-', type: 'number', default: 3, min: 1, max: 12, group: 'Crush loads', role: 'geometry', help: 'Number of elements of the layer that carry one pad or chute line load after spreading through sheath, armour and fillers (1 = no spreading, most conservative). Calibrate against the crush test.' },
    { key: 'padMu', label: 'Pad / sheath friction coefficient', unit: '-', type: 'number', default: 0.15, min: 0.03, max: 0.6, group: 'Crush loads', role: 'operation' },
    { key: 'gripSF', label: 'Grip safety factor', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'etaCrush', label: 'Allowable ring-bending stress / yield', unit: '-', type: 'number', default: 0.8, min: 0.3, max: 1, group: 'Acceptance criteria', role: 'criteria', help: `Installation (abnormal) utilisation applied to the first-yield crush criterion. ${CONFIRM}` },
    { key: 'lengthKm', label: 'Manufactured length', unit: 'km', type: 'number', default: 5, min: 0.05, max: 200, group: 'Quantities', role: 'operation', bind: 'cable.length' },
  ],
  compute(inp, ctx) {
    const b = buildSection(inp, ctx), s = b.sec;
    need(b.tubes.length > 0, 'No steel tube (hollow component with a yield stress) is defined: there is nothing to check for collapse or crushing.');
    const pe = inp.rhoW * G * inp.depth, pmin = inp.pMin * 1e5, gam = inp.gammaM * inp.gammaSC, f0 = inp.ovality / 100;
    const ls = umbilicalLoadSharing(s, { tension: inp.T * 1e3 });
    const rows = b.tubes.map((e) => {
      const c = e.component, D = c.od, t = (c.od - c.id) / 2, idx = s.elements.indexOf(e);
      const col = run(() => tubeCollapsePressure({ D, t, E: c.E, fy: c.yieldStress, ovality: f0, alphaFab: inp.alphaFab }), `Collapse of ${e.name}`);
      const vm = lameStresses({ pInt: pmin, pExt: pe, ri: c.id / 2, ro: c.od / 2 }).vonMises;
      const cr = run(() => umbilicalCrushCheck({ tubeOd: D, tubeWall: t, allowableStress: inp.etaCrush * c.yieldStress, padLoad: inp.squeeze * 1e3, nTracks: Number(inp.nTracks), elementsUnderPad: inp.spread, alpha: e.alpha,
        tension: inp.T * 1e3, chuteRadius: inp.chuteR, mu: inp.padMu, safetyFactor: inp.gripSF }), `Crush of ${e.name}`);
      const qc = ls.elements[idx].contactInner;       // N/m per unit element length from the helix tension of this and outer layers
      return { name: e.name, D, t, dt: D / t, col, vm, fy: c.yieldStress, cr, qc, ucCol: Math.max(pe - pmin, 0) / (col.pc / gam), ucVm: vm / c.yieldStress, ucContact: qc / cr.tubeLineLoadCapacity };
    });
    const checks = [];
    for (const r of rows) {
      checks.push(check(`${r.name}: hydrostatic collapse`, Math.max(pe - pmin, 0) / 1e5, r.col.pc / gam / 1e5, 'bar', 'DNV-ST-F101 (2021) §5.4.4 Eq. 5.10 (pe − pmin ≤ pc/(γm·γSC))', 'Net external pressure on the depressurised tube against its factored collapse pressure. Above 1.0 the tube flattens; a collapsed control line is lost and may propagate along the tube.'));
      checks.push(check(`${r.name}: tensioner crush`, r.cr.tubeLineLoadTensioner / 1e3, r.cr.tubeLineLoadCapacity / 1e3, 'kN/m', `${STD_UMB} §7 (crush load); Roark Table 9.2 ring, M = 0.318·P·R`, 'Line load squeezed into the tube by one tensioner track against the load at which ring bending reaches the allowable stress. Exceeding it ovalises the tube permanently.'));
      checks.push(check(`${r.name}: chute line load`, r.cr.utilisationChute * r.cr.tubeLineLoadCapacity / 1e3, r.cr.tubeLineLoadCapacity / 1e3, 'kN/m', `${STD_UMB} §7; sidewall load T/R`, 'Share of the chute sidewall load T/R reaching the tube against its ring-bending capacity. Governs the minimum chute radius at a given lay tension.'));
    }
    const gov = worst(checks);
    const wCol = rows.reduce((a, r) => (r.ucCol > a.ucCol ? r : a), rows[0]), wCr = rows.reduce((a, r) => (r.cr.utilisationTensioner > a.cr.utilisationTensioner ? r : a), rows[0]);
    const allowPad = Math.min(...rows.map((r) => r.cr.allowablePadLoad)), minChute = Math.max(...rows.map((r) => r.cr.minimumChuteRadius)), gripL = Math.max(...rows.map((r) => r.cr.requiredTensionerLength));
    const collapseDepth = wCol.col.pc / gam / (inp.rhoW * G);
    const depths = linspace(0, Math.max(inp.depth * 1.5, collapseDepth * 1.2), 30);
    const ovs = linspace(0.5, 5, 19);
    const steelTubesT = b.tubes.reduce((a, e) => a + e.count * e.component.mass / Math.cos(e.alpha), 0) * inp.lengthKm;
    // wall needed for the governing tube to reach collapse utilisation 0.9 (bisection on t at constant OD)
    let tReq = wCol.t;
    if (wCol.ucCol > 0.9) { let lo = wCol.t, hi = wCol.D / 2 * 0.95; for (let i = 0; i < 50; i++) { const mid = 0.5 * (lo + hi); const pc = tubeCollapsePressure({ D: wCol.D, t: mid, E: b.tubes[0].component.E, fy: wCol.fy, ovality: f0, alphaFab: inp.alphaFab }).pc; if (Math.max(pe - pmin, 0) / (pc / gam) > 0.9) lo = mid; else hi = mid; } tReq = hi; }
    return {
      kpis: [
        { key: 'pe', label: 'External hydrostatic pressure', value: round(pe / 1e5), unit: 'bar', explain: 'Sea-water pressure at the seabed acting on the tubes (the bundle is flooded between elements).' },
        { key: 'collapseUc', label: 'Collapse utilisation (governing tube)', value: round(wCol.ucCol), unit: '-', limit: 1, sense: 'max', explain: `Net external pressure over factored collapse pressure for "${wCol.name}" when depressurised. Above 1.0 the empty tube is predicted to collapse at this depth.` },
        { key: 'pc', label: 'Collapse pressure, governing tube', value: round(wCol.col.pc / 1e5), unit: 'bar', explain: 'Characteristic collapse pressure from the elastic–plastic interaction with ovality, before resistance factors.' },
        { key: 'pel', label: 'Elastic collapse pressure, governing tube', value: round(wCol.col.pel / 1e5), unit: 'bar', explain: 'Elastic instability pressure 2E(t/D)³/(1−ν²) of a perfectly round tube. For the thick small-bore tubes of umbilicals it is far above the plastic value, so yield governs.' },
        { key: 'pp', label: 'Plastic collapse pressure, governing tube', value: round(wCol.col.pp / 1e5), unit: 'bar', explain: 'Pressure at which the wall yields in hoop compression, fy·αfab·2t/D.' },
        { key: 'collapseDepth', label: 'Allowable depth, empty governing tube', value: round(collapseDepth), unit: 'm', explain: 'Water depth at which the depressurised governing tube reaches its factored collapse pressure — the depth rating of the tube when vented.' },
        { key: 'vmExt', label: 'Wall stress utilisation under external pressure', value: round(Math.max(...rows.map((r) => r.ucVm))), unit: '-', explain: 'Bore von Mises stress of the depressurised tube under hydrostatic pressure over yield (Lamé). A cross-check on collapse for thick walls.' },
        { key: 'crushUc', label: 'Tensioner crush utilisation', value: round(wCr.cr.utilisationTensioner), unit: '-', limit: 1, sense: 'max', explain: `Line load on "${wCr.name}" from one tensioner track over its ring-bending capacity at the allowable stress.` },
        { key: 'chuteUc', label: 'Chute crush utilisation', value: round(Math.max(...rows.map((r) => r.cr.utilisationChute))), unit: '-', limit: 1, sense: 'max', explain: 'Share of the chute sidewall load T/R on the worst tube over its ring-bending capacity.' },
        { key: 'padAllow', label: 'Allowable squeeze per track', value: round(allowPad / 1e3), unit: 'kN/m', explain: 'Largest tensioner squeeze per track that keeps every tube within the allowable ring-bending stress, with the stated load spreading. This is the crushing capacity to give the installation contractor.' },
        { key: 'minChute', label: 'Minimum chute radius for crush', value: round(minChute), unit: 'm', explain: 'Chute or sheave radius below which the sidewall load T/R crushes the worst tube at this tension (bending MBR is checked separately).' },
        { key: 'gripLength', label: 'Tensioner length needed for grip', value: round(gripL), unit: 'm', explain: 'Contact length needed to hold the tension by friction at the chosen squeeze and safety factor; shorter tensioners need more squeeze, which raises crush utilisation.' },
        { key: 'contactUc', label: 'Armour contact-load utilisation', value: round(Math.max(...rows.map((r) => r.ucContact))), unit: '-', explain: 'Radial line load squeezed onto the tubes by the tensioned helical layers over the tube ring-bending capacity; it shows how much crush capacity tension alone consumes.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Collapse utilisation vs water depth (tubes empty)', x: { label: 'Water depth', unit: 'm' }, y: { label: 'Utilisation', unit: '-' },
          series: [...rows.map((r) => ({ name: r.name, x: depths, y: depths.map((d) => Math.max(inp.rhoW * G * d - pmin, 0) / (r.col.pc / gam)) })), { name: 'Limit', x: [0, depths[depths.length - 1]], y: [1, 1], dash: true }, { name: 'Site depth', x: [inp.depth, inp.depth], y: [0, 1], dash: true }],
          explain: 'Each line rises linearly with depth; where it crosses 1.0 is the depth rating of that tube when depressurised. Tubes with the highest D/t cross first.' },
        { type: 'line', title: `Collapse pressure vs ovality — ${wCol.name}`, x: { label: 'Ovality', unit: '%' }, y: { label: 'Collapse pressure', unit: 'bar' },
          series: [{ name: 'pc', x: ovs, y: ovs.map((o) => tubeCollapsePressure({ D: wCol.D, t: wCol.t, E: b.tubes[0].component.E, fy: wCol.fy, ovality: o / 100, alphaFab: inp.alphaFab }).pc / 1e5) }, { name: 'Factored demand', x: [ovs[0], ovs[ovs.length - 1]], y: [(pe - pmin) * gam / 1e5, (pe - pmin) * gam / 1e5], dash: true }],
          explain: 'Collapse pressure falls as the tube becomes more oval — which is what tensioner squeeze and reeling do. The crossing with the dashed demand line is the ovality at which the vented tube would collapse at this depth.' },
        { type: 'bar', title: 'Crush utilisation by tube', categories: rows.map((r) => r.name), series: [{ name: 'Tensioner', y: rows.map((r) => round(r.cr.utilisationTensioner)) }, { name: 'Chute', y: rows.map((r) => round(r.cr.utilisationChute)) }, { name: 'Armour contact under tension', y: rows.map((r) => round(r.ucContact)) }],
          explain: 'Ring-bending utilisation of each tube under the three contact load sources. The tallest bar sets the allowable squeeze or the minimum chute radius.' },
      ],
      tables: [{ title: 'Tube collapse and crush', explain: 'One row per tube type: slenderness, elastic, plastic and combined collapse pressure, collapse utilisation at the site depth, line-load capacity of the tube ring and the crush utilisations.',
        columns: [{ key: 'name', label: 'Tube' }, { key: 'dt', label: 'D/t', unit: '-' }, { key: 'pel', label: 'p elastic', unit: 'bar' }, { key: 'pp', label: 'p plastic', unit: 'bar' }, { key: 'pc', label: 'p collapse', unit: 'bar' }, { key: 'ucCol', label: 'Collapse UC', unit: '-' }, { key: 'cap', label: 'Ring capacity', unit: 'kN/m' }, { key: 'ucT', label: 'Tensioner UC', unit: '-' }, { key: 'ucC', label: 'Chute UC', unit: '-' }],
        rows: rows.map((r) => ({ name: r.name, dt: round(r.dt, 3), pel: round(r.col.pel / 1e5), pp: round(r.col.pp / 1e5), pc: round(r.col.pc / 1e5), ucCol: round(r.ucCol), cap: round(r.cr.tubeLineLoadCapacity / 1e3), ucT: round(r.cr.utilisationTensioner), ucC: round(r.cr.utilisationChute) })) }],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        `External pressure ${fmt(pe / 1e5)} bar at ${inp.depth} m; the weakest tube in collapse (${wCol.name}, D/t = ${fmt(wCol.dt)}) is rated to ${fmt(collapseDepth)} m when empty.`,
        `Allowable squeeze ${fmt(allowPad / 1e3)} kN/m per track assumes each pad load is shared by ${inp.spread} elements — confirm with the manufacturer crush test.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}. Collapse utilisation ${fmt(wCol.ucCol)}, tensioner crush ${fmt(wCr.cr.utilisationTensioner)}, allowable squeeze ${fmt(allowPad / 1e3)} kN/m, minimum chute radius for crush ${fmt(minChute)} m.`,
          gov.uc > 1 ? (gov.name.includes('collapse') ? `Increase the wall of "${wCol.name}" to at least ${fmt(tReq * 1e3)} mm (collapse utilisation 0.9) or keep the tube pressurised / liquid-filled during installation and service.` : `Limit the squeeze to ${fmt(allowPad / 1e3)} kN/m per track and use a tensioner of at least ${fmt(inp.gripSF * inp.T * 1e3 / (inp.padMu * Number(inp.nTracks) * allowPad))} m contact length, or more tracks.`)
            : 'State the allowable squeeze, minimum chute radius and the requirement to keep tubes filled in the installation specification; confirm the load-spreading assumption with a full-section crush test.',
          `DNV-ST-F101 (2021) §5.4.4; ${STD_UMB} §7 and qualification crush test`),
        adv('economic', 'info', `Tubes are ${fmt(1 / Math.max(wCol.ucCol, 1e-6))} × stronger than needed against collapse at ${inp.depth} m; the allowable squeeze of ${fmt(allowPad / 1e3)} kN/m needs ${fmt(inp.gripSF * inp.T * 1e3 / (inp.padMu * Number(inp.nTracks) * allowPad))} m of ${inp.nTracks}-track tensioner at ${inp.T} kN.`,
          'A low allowable squeeze forces a longer or multi-track tensioner and can exclude smaller lay vessels; a modest increase in tube wall or filler stiffness that raises the crush capacity can widen the vessel market more cheaply than chartering a larger spread.',
          'Installation-vessel selection on capability and day-rate basis (DNV-ST-N001 planning; AACE RP 18R-97 Class 4)',
          `Steel tubes total ${fmt(steelTubesT)} t over ${inp.lengthKm} km ≈ ${fmt(steelTubesT * CO2_STEEL)} t CO₂e at the screening factor ${CO2_STEEL} t/t (super-duplex is higher). A crushed section found at the post-lay test means recovery and re-lay: ${vesselCo2Text(5)}`),
      ],
      publish: { 'umbilical.allowSqueeze': allowPad, 'umbilical.collapseUtilisation': wCol.ucCol, 'umbilical.minChuteRadius': minChute },
      quantities: { steel_t: steelTubesT },
    };
  },
  verify: [
    // Thin tube D = 20 mm, t = 0.5 mm (D/t = 40), E = 200 GPa, ν = 0.3, fy = 550 MPa, αfab = 1, f0 = 0.5 %:
    // pel = 2E(t/D)³/(1 − ν²) = 2 × 200e9 × 0.025³/0.91 = 6.868 MPa = 68.68 bar; pp = fy·2t/D = 27.5 MPa = 275 bar;
    // cubic (pc − pel)(pc² − pp²) = pc·pel·pp·f0·D/t solved independently by bisection: pc = 6.523 MPa = 65.23 bar.
    { name: 'Elastic and combined collapse pressure of a thin tube',
      source: 'Timoshenko & Gere, Theory of Elastic Stability §7.5 (pel); DNV-ST-F101 (2021) Eq. 5.10 cubic solved by hand',
      inputs: { sheathT: 0, layers: [], alphaFab: 1, ovality: 0.5, depth: 100, components: [{ name: 'Tube', kind: 'tube', layer: 0, count: 1, od: 20, wall: 0.5, E: 200, rho: 7800, fy: 550, p: 0, rhoF: 0 }] },
      expect: [{ kpi: 'pel', value: 68.68, tol: 0.002 }, { kpi: 'pp', value: 275, tol: 0.002 }, { kpi: 'pc', value: 65.23, tol: 0.003 }] },
    // Ring under two opposed line loads (Roark Table 9.2): M = P·R/π. First yield at σ = 6M/t² → P = σ·t²·π/(6R);
    // D = 20, t = 2 mm → R = 9 mm, σ = 0.8 × 550 = 440 MPa: P = 440e6 × 4e-6 × π/(6 × 0.009) = 102.4 kN/m.
    { name: 'Crush capacity of a tube ring under opposed line loads: P = π·σ·t²/(6·R)',
      source: "Roark's Formulas for Stress and Strain, Table 9.2 case 1; hand calculation",
      inputs: { sheathT: 0, layers: [], spread: 1, depth: 100, components: [{ name: 'Tube', kind: 'tube', layer: 0, count: 1, od: 20, wall: 2, E: 200, rho: 7800, fy: 550, p: 0, rhoF: 0 }] },
      expect: [{ kpi: 'padAllow', value: 102.4, tol: 0.003 }] },
  ],
};
