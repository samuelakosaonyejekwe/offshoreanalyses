// Wind loading of the topsides: Frøya profile and gust conversion, projected-area wind forces per deck level
// for eight directions, frame response with gravity, and the wind resultant at the jacket interface.
import { froyaWind, windForce, RHO_AIR, API_SHAPE_COEFFICIENTS } from '../../engine/metocean/wind.js';
import { deckInputs, solverInput, buildDeck, loadCase, solveDeck, memberForces, memberChecks, framePlot, interfaceResultant, equipmentOf, need, fmt, round, linspace, adv, check, steelCo2Text, CONFIRM, DECK_LIMITS, DEG } from './_common.js';

export default {
  id: 's02-wind-load', suite: 2, title: 'Wind load on the topsides',
  summary: 'Wind speed profile, wind forces on each deck level by direction, base shear and overturning moment at the jacket interface, and the frame utilisation and sway under wind with gravity.',
  covers: [11], fidelity: 'numerical',
  standards: ['ISO 19901-1:2015 §A.7.3 (Frøya wind profile and gust averaging)', 'API RP 2A-WSD (2014) §5.3.2 (wind force, shape coefficients)', 'DNV-RP-C205 (2019) §5.3 (wind loads)', 'ISO 19901-3:2014 §7.3 (environmental actions on topsides)', 'NORSOK N-003 (2017) §6.3'],
  models: ['Frøya (NPD) wind profile and gust factor', 'Projected-area wind pressure F = ½ρ Cs A U²', 'Static equilibrium', 'Matrix stiffness method (3-D direct stiffness)'],
  limits: `${DECK_LIMITS} Quasi-static wind on block projected areas per deck level with one shape coefficient: shielding between modules, lift on deck undersides, vortex shedding of slender items and dynamic gust response are not modelled. Wind-tunnel or CFD pressure distributions are outside this module — enter their resultants through the block areas and coefficients.`,
  publishes: ['topsides.windFx', 'topsides.windFy', 'topsides.windMy'],
  searchTerms: 'wind load offshore topsides Froya profile projected area shape coefficient base shear',
  inputs: [
    ...deckInputs(),
    { key: 'U0', label: '1-hour mean wind speed at 10 m', unit: 'm/s', type: 'number', default: 38, min: 5, max: 80, group: 'Wind', role: 'environment', bind: 'site.Uw100', uncertain: { dist: 'gumbel', cov: 0.1 }, help: 'Extreme (100-year) 1-hour mean at 10 m above sea level.' },
    { key: 'tAvg', label: 'Gust averaging period', type: 'select', default: '60', group: 'Wind', role: 'environment', options: [{ value: '3', label: '3 s gust (individual members, small items)' }, { value: '15', label: '15 s (modules smaller than 50 m)' }, { value: '60', label: '1 min (whole topsides, with extreme waves)' }, { value: '3600', label: '1 h mean' }],
      help: `API RP 2A / ISO 19901-1: shorter gusts for smaller structures. ${CONFIRM}` },
    { key: 'zInt', label: 'Jacket interface above mean sea level', unit: 'm', type: 'number', default: 18, min: 3, max: 80, group: 'Wind', role: 'bc', help: 'Elevation of the deck-leg interface; deck elevations are measured from it.' },
    { key: 'rhoAir', label: 'Air density', unit: 'kg/m³', type: 'number', default: RHO_AIR, min: 1, max: 1.5, group: 'Wind', role: 'environment' },
    { key: 'hBlock', label: 'Exposed height of each deck level', unit: 'm', type: 'number', default: 6, min: 0.5, max: 30, group: 'Exposed areas', role: 'geometry', help: 'Height of the block of equipment, modules and cladding standing on each deck that the wind sees.' },
    { key: 'solid', label: 'Solidity of the exposed block', unit: '-', type: 'number', default: 0.7, min: 0.05, max: 1, group: 'Exposed areas', role: 'geometry', uncertain: { dist: 'normal', cov: 0.15 }, help: 'Projected solid area over the enclosing rectangle: 1.0 for clad modules, 0.4–0.7 for open process decks.' },
    { key: 'Cs', label: 'Shape coefficient', unit: '-', type: 'number', default: API_SHAPE_COEFFICIENTS.overallPlatformProjectedArea, min: 0.3, max: 2.5, group: 'Exposed areas', role: 'load', help: `API RP 2A: 1.0 overall projected area of a platform, 1.5 beams and sides of buildings, 0.5 cylinders. ${CONFIRM}` },
    { key: 'psi', label: 'Fraction of live load acting with the wind', unit: '-', type: 'number', default: 0.75, min: 0, max: 1, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'gG', label: 'Partial factor on permanent and variable actions', unit: '-', type: 'number', default: 1.1, min: 0.8, max: 1.6, group: 'Acceptance criteria', role: 'criteria', help: `ISO 19902 extreme combination 1.1 G + 1.1 Q + 1.35 E. ${CONFIRM}` },
    { key: 'gE', label: 'Partial factor on the wind action', unit: '-', type: 'number', default: 1.35, min: 1, max: 1.8, group: 'Acceptance criteria', role: 'criteria', help: CONFIRM },
    { key: 'driftRatio', label: 'Sway limit (height / value)', unit: '-', type: 'number', default: 300, min: 100, max: 1000, group: 'Acceptance criteria', role: 'criteria' },
    solverInput(1),
  ],
  mesh: { key: 'nSub', label: 'Elements per member', levels: [1, 2, 4], qoi: 'sway' },
  compute(inp, ctx) {
    const deck = buildDeck(inp, ctx), eqs = equipmentOf(deck, inp), zs = deck.levelsZ, t = Number(inp.tAvg);
    // Wind force on each level block for wind along x (face Ly × hBlock) and along y (face Lx × hBlock), characteristic.
    const lev = zs.map((z, k) => {
      const zc = Math.max(inp.zInt + z + inp.hBlock / 2, 1), w = froyaWind({ U0: inp.U0, z: zc, t });
      const fx = windForce({ U: w.u, A: deck.Ly * inp.hBlock * inp.solid, Cs: inp.Cs, rho: inp.rhoAir }), fy = windForce({ U: w.u, A: deck.Lx * inp.hBlock * inp.solid, Cs: inp.Cs, rho: inp.rhoAir });
      return { k, z, zc, u: w.u, Uz: w.Uz, q: fx.q, Fx: fx.F, Fy: fy.F };
    });
    const Fx0 = lev.reduce((s, l) => s + l.Fx, 0), Fy0 = lev.reduce((s, l) => s + l.Fy, 0);
    const My0 = lev.reduce((s, l) => s + l.Fx * (l.z + inp.hBlock / 2), 0), Mx0 = lev.reduce((s, l) => s + l.Fy * (l.z + inp.hBlock / 2), 0);
    // Eight directions; each level force goes to its nodes with the couple of the block centroid above the deck.
    const dirs = linspace(0, 315, 8), cases = dirs.map((d, n) => {
      const c = Math.cos(d * DEG), s = Math.sin(d * DEG), extra = [];
      lev.forEach((l) => { const ns = deck.levelNodes[l.k], fx = inp.gE * l.Fx * c / ns.length, fy = inp.gE * l.Fy * s / ns.length; for (const node of ns) extra.push({ node, F: [fx, fy, 0, -fy * inp.hBlock / 2, fx * inp.hBlock / 2, 0] }); });
      return loadCase(deck, inp, { name: `d${n}`, gDead: inp.gG, gEquip: inp.gG, gLive: inp.gG * inp.psi, extra, eqs });
    });
    cases.push(loadCase(deck, inp, { name: 'grav', gDead: inp.gG, gEquip: inp.gG, gLive: inp.gG * inp.psi, eqs }));
    const { ref, sol } = solveDeck(deck, cases, inp.nSub), nN = ref.model.nodes.length, u0 = sol.results.grav.displacements;
    const res = dirs.map((d, n) => {
      const r = sol.results[`d${n}`], rows = memberChecks(deck, memberForces(deck, ref, r)), gov = rows.reduce((a, c) => (c.uc > a.uc ? c : a)), ir = interfaceResultant(deck, r);
      let sway = 0; for (let i = 0; i < nN; i++) sway = Math.max(sway, Math.hypot(r.displacements[6 * i] - u0[6 * i], r.displacements[6 * i + 1] - u0[6 * i + 1]));
      return { d, r, rows, gov, ir, sway, shear: Math.hypot(ir.R[0], ir.R[1]), otm: Math.hypot(ir.R[3], ir.R[4]), minLeg: Math.min(...ir.rows.map((x) => x.Rz)), maxLeg: Math.max(...ir.rows.map((x) => x.Rz)) };
    });
    const g = res.reduce((a, c) => (c.gov.uc > a.gov.uc ? c : a)), sw = res.reduce((a, c) => (c.sway > a.sway ? c : a));
    const ucGrav = Math.max(...memberChecks(deck, memberForces(deck, ref, sol.results.grav)).map((x) => x.uc));
    const H = zs[zs.length - 1], swayLim = H / inp.driftRatio, swayChar = sw.sway / inp.gE;
    const minLeg = Math.min(...res.map((x) => x.minLeg));
    need(Math.abs(res[0].ir.R[0] + inp.gE * Fx0) < 1e-3 * inp.gE * Fx0 + 1, 'Equilibrium check failed: the horizontal reactions do not balance the wind force.');
    const checks = [
      check(`Member strength with wind from ${g.d}°`, g.gov.uc, 1, 'interaction', g.gov.ref, `Highest utilisation of any member (${g.gov.group.toLowerCase()} ${g.gov.i + 1}) under factored gravity with the factored wind from the worst of eight directions. Above 1.0 the member yields or buckles in the design storm.`),
      check(`Sway under characteristic wind from ${sw.d}°`, swayChar * 1e3, swayLim * 1e3, 'mm', `Height / ${inp.driftRatio} (ISO 19901-3 §7.5 serviceability)`, 'Horizontal movement caused by the unfactored wind alone. Excess sway strains bridge landings, risers and piping crossing to adjacent structures.'),
    ];
    if (minLeg < 0) checks.push(check('Leg uplift under wind', -minLeg / 1e3, Math.max(...res.map((x) => x.maxLeg)) / 1e3, 'kN', 'ISO 19901-3 §9 (interface connections)', 'A deck leg goes into tension under overturning; its connection to the jacket must carry the uplift.'));
    const gov = checks.reduce((a, c) => (c.uc > a.uc ? c : a)), ok = gov.uc <= 1;
    const zp = linspace(Math.max(inp.zInt, 5), inp.zInt + H + inp.hBlock + 5, 25);
    const braceT = deck.orig.filter((o) => o.group === 'Brace').reduce((s, o) => s + o.sec.A * o.L * deck.steel.rho, 0) / 1e3;
    return {
      kpis: [
        { key: 'Fx', label: 'Wind force along x (characteristic)', unit: 'kN', value: Fx0 / 1e3, explain: 'Unfactored wind force on the topsides for wind blowing along x, summed over the deck levels.' },
        { key: 'Fy', label: 'Wind force along y (characteristic)', unit: 'kN', value: Fy0 / 1e3, explain: 'Unfactored wind force for wind along y (the broadside when the deck is longer in x).' },
        { key: 'My', label: 'Wind overturning moment at the interface (x wind)', unit: 'kN·m', value: My0 / 1e3, explain: 'Moment of the x-direction wind forces about the jacket interface: it changes the leg loads by ± M / leg spacing.' },
        { key: 'uTop', label: 'Design gust speed at the top level', unit: 'm/s', value: lev[lev.length - 1].u, explain: `Wind speed averaged over ${inp.tAvg} s at the centroid of the top block (${fmt(lev[lev.length - 1].zc)} m above sea level).` },
        { key: 'qTop', label: 'Wind pressure at the top level', unit: 'kPa', value: lev[lev.length - 1].q / 1e3, explain: 'Basic velocity pressure ½ρU² at the top block before the shape coefficient.' },
        { key: 'maxUc', label: 'Highest member utilisation with wind', unit: '-', value: round(g.gov.uc, 5), limit: 1, sense: 'max', explain: `Governing utilisation over all wind directions (gravity alone gives ${fmt(ucGrav)}): the difference is what the wind costs the frame.` },
        { key: 'sway', label: 'Sway under characteristic wind', unit: 'mm', value: swayChar * 1e3, limit: swayLim * 1e3, sense: 'max', explain: 'Largest horizontal displacement caused by the unfactored wind, over all directions.' },
        { key: 'baseShear', label: 'Factored base shear (worst direction)', unit: 'kN', value: Math.max(...res.map((x) => x.shear)) / 1e3, explain: 'Largest horizontal reaction at the interface under the factored wind.' },
        { key: 'minLeg', label: 'Smallest leg reaction', unit: 'kN', value: minLeg / 1e3, explain: 'Lowest vertical leg reaction over all directions; a negative value is uplift.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Wind speed profile', x: { label: 'Wind speed', unit: 'm/s' }, y: { label: 'Height above sea level', unit: 'm' }, series: [{ name: `${inp.tAvg} s gust`, x: zp.map((z) => round(froyaWind({ U0: inp.U0, z, t }).u)), y: zp.map((z) => round(z)) }, { name: '1-hour mean', x: zp.map((z) => round(froyaWind({ U0: inp.U0, z, t: 3600 }).Uz)), y: zp.map((z) => round(z)), dash: true }, { name: 'Deck block centroids', x: lev.map((l) => round(l.u)), y: lev.map((l) => round(l.zc)), mode: 'markers' }],
          explain: 'Frøya profile of the mean and the gust speed used for the load; markers show the speed taken for each deck block. The gust curve lies further above the mean near the sea surface where turbulence is highest.' },
        { type: 'polar', title: 'Response by wind direction', series: [{ name: 'Member utilisation', thetaDeg: [...dirs, 360], r: [...res.map((x) => round(x.gov.uc)), round(res[0].gov.uc)] }, { name: 'Base shear / largest', thetaDeg: [...dirs, 360], r: [...res.map((x) => round(x.shear / Math.max(...res.map((y) => y.shear)))), round(res[0].shear / Math.max(...res.map((y) => y.shear)))] }],
          explain: 'Governing member utilisation and normalised base shear for wind towards each direction (0° = +x). A lobe shows the weak direction of the bracing system.' },
        framePlot({ title: `Utilisation and deformed shape, wind from ${g.d}°`, ref, u: g.r.displacements, values: (() => { const v = new Array(nN).fill(0); g.rows.forEach((x) => { for (const n of ref.nodeMap[x.i]) v[n] = Math.max(v[n], x.uc); }); return v.map((x) => round(x)); })(), valueLabel: 'Utilisation', valueUnit: '-', explain: 'Frame coloured by member utilisation for the governing wind direction with gravity. Windward legs unload and leeward legs and braces pick up the overturning.' }),
        { type: 'bar', title: 'Wind force by deck level', categories: lev.map((l) => `Level ${l.k + 1}`), series: [{ name: 'Along x', y: lev.map((l) => round(l.Fx / 1e3)) }, { name: 'Along y', y: lev.map((l) => round(l.Fy / 1e3)) }], y: { label: 'Characteristic force', unit: 'kN' }, explain: 'Unfactored wind force on the exposed block of each deck. Higher decks see more wind per square metre and a longer lever arm to the interface.' },
      ],
      tables: [
        { title: 'Wind by deck level', columns: [{ key: 'level', label: 'Level' }, { key: 'zc', label: 'Centroid above sea', unit: 'm' }, { key: 'u', label: 'Gust speed', unit: 'm/s' }, { key: 'q', label: 'Pressure', unit: 'kPa' }, { key: 'Fx', label: 'Fx', unit: 'kN' }, { key: 'Fy', label: 'Fy', unit: 'kN' }],
          rows: lev.map((l) => ({ level: l.k + 1, zc: round(l.zc), u: round(l.u), q: round(l.q / 1e3), Fx: round(l.Fx / 1e3), Fy: round(l.Fy / 1e3) })), explain: 'Characteristic wind speed, velocity pressure and force for each deck block. These forces are the wind input for the in-place and interface analyses.' },
        { title: 'Results by direction', columns: [{ key: 'd', label: 'Direction', unit: 'deg' }, { key: 'shear', label: 'Base shear', unit: 'kN' }, { key: 'otm', label: 'Overturning moment', unit: 'kN·m' }, { key: 'minLeg', label: 'Min leg reaction', unit: 'kN' }, { key: 'maxLeg', label: 'Max leg reaction', unit: 'kN' }, { key: 'sway', label: 'Sway (factored)', unit: 'mm' }, { key: 'uc', label: 'Max UC' }],
          rows: res.map((x) => ({ d: x.d, shear: round(x.shear / 1e3), otm: round(x.otm / 1e3), minLeg: round(x.minLeg / 1e3), maxLeg: round(x.maxLeg / 1e3), sway: round(x.sway * 1e3), uc: round(x.gov.uc) })), explain: 'Factored interface resultants and frame response for each of the eight wind directions (overturning moment includes the eccentricity of gravity loads about the plan centre).' },
      ],
      notes: [`Wind from ${g.d}° governs strength at ${fmt(g.gov.uc)} (${g.gov.group.toLowerCase()} ${g.gov.i + 1}); gravity alone gives ${fmt(ucGrav)}.`, `Characteristic wind force ${fmt(Fx0 / 1e3)} kN along x and ${fmt(Fy0 / 1e3)} kN along y with the ${inp.tAvg} s gust.`],
      advice: [
        ok ? adv('engineering', gov.uc > 0.9 ? 'warn' : 'ok', `Wind adds ${fmt(Math.max(g.gov.uc - ucGrav, 0))} to the governing utilisation (${fmt(g.gov.uc)} with wind from ${g.d}°); sway ${fmt(swayChar * 1e3)} mm against ${fmt(swayLim * 1e3)} mm.`,
          'Pass the level forces to the in-place and jacket analyses, and check cladding, flare boom, crane boom rest and other slender items separately with the 3 s gust and their own shape coefficients.', 'ISO 19901-1:2015 §A.7.3; API RP 2A-WSD (2014) §5.3.2; ISO 19902:2020 §9.10',
          `Wind is carried mainly by the vertical bracing (${fmt(braceT)} t); lowering the exposed solidity with open grating and wind walls only where needed reduces the load at source without steel (${steelCo2Text(braceT)} is the bracing in place).`)
          : adv('engineering', 'fail', `${gov.name} is not met (utilisation ${fmt(gov.uc)}): base shear ${fmt(Math.max(...res.map((x) => x.shear)) / 1e3)} kN factored, sway ${fmt(swayChar * 1e3)} mm.`,
            gov.name.startsWith('Sway') ? 'Add or stiffen vertical bracing in the faces parallel to the governing wind direction; bracing is far more effective against sway than heavier legs.' : `Strengthen the ${g.gov.group.toLowerCase()} group by about ${fmt((g.gov.uc - 1) * 100, 2)} % in area or modulus, add a braced bay in the governing direction, or reduce the exposed area or solidity.`,
            'ISO 19902:2020 §9.10, §13; API RP 2A-WSD (2014) §5.3.2', `Strengthening by ${fmt((gov.uc - 1) * 100, 2)} % in the governing group costs in the order of ${steelCo2Text(Math.max(gov.uc - 1, 0.05) * Math.max(braceT, 0.1 * deck.frameSteelT))}`),
        adv('economic', 'info', `The wind resultant of ${fmt(Math.max(Fx0, Fy0) / 1e3)} kN with a ${fmt(My0 / Math.max(Fx0, 1))} m lever arm is passed to the jacket and its piles as well: topsides wind is typically 10–25 % of the global base shear of a fixed platform in a storm.`,
          'Fix the exposed-area estimate early (module envelopes, cladding extent) and carry a contingency on area as on weight; late growth of wind area is paid for in bracing, jacket and pile steel together.', 'ISO 19901-5:2021 (control of weight and area growth); AACE RP 18R-97'),
      ],
      publish: { 'topsides.windFx': Fx0, 'topsides.windFy': Fy0, 'topsides.windMy': My0 },
    };
  },
  verify: [{
    // One level at 10 m above sea level (interface 3 m, deck 4 m, block centroid 3 m above the deck), 1-hour mean
    // 40 m/s with no gust conversion: U = U0 at 10 m. Face 20 m × 6 m, solidity 1, Cs = 1:
    //   Fx = ½·1.225·40²·120 = 117.6 kN; along y the face is 30 m × 6 m → 176.4 kN; My = 117.6·(4 + 3) = 823.2 kN·m.
    name: 'Wind force on a block at the reference height: F = ½ρ Cs A U²', source: 'API RP 2A-WSD (2014) Eq. 5.6 (hand calculation)',
    inputs: { U0: 40, tAvg: '3600', zInt: 3, hBlock: 6, solid: 1, Cs: 1, rhoAir: 1.225, nx: 1, ny: 1, nLevels: 1, bracing: 'none', columns: 'corners', live: 0, equip: [{ name: 'Block', W: 500, x: 15, y: 10, level: 1, a: 0, b: 0, zc: 0 }] },
    expect: [{ kpi: 'Fx', value: 117.6, tol: 1e-3 }, { kpi: 'Fy', value: 176.4, tol: 1e-3 }, { kpi: 'My', value: 823.2, tol: 1e-3 }],
  }],
};
