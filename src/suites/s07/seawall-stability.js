// Gravity seawall / quay-type retaining wall: sliding, overturning and foundation bearing under active earth
// pressure, differential water pressure, uplift, surcharge and an optional line load, with sea-level rise.
import { coulombKa, rankine, pressureDiagram, gravityWallStability } from '../../engine/geotech/earthpressure.js';
import { bearingCapacity } from '../../engine/geotech/bearing.js';
import { G, DEG, CO2_CONC, need, run, finite, r, f, linspace, check, sfCheck, eng, eco, slrInputs, slrRows, SLR_BASIS } from './_common.js';

export default {
  id: 's07-seawall-stability', suite: 7, title: 'Seawall sliding, overturning and bearing',
  summary: 'External stability of a gravity seawall at low water: earth and water thrust, uplift, base pressure, bearing capacity and the base width the margins require.',
  covers: [5], fidelity: 'analytical',
  standards: ['EN 1997-1:2004 §6, §9 and Annex D', 'BS 6349-2:2019 (quay walls, maritime retaining structures)', 'USACE EM 1110-2-2502 Retaining and Flood Walls', 'CIRIA C683 Rock Manual (2007) §5.4'],
  models: ['Coulomb / Rankine active earth pressure', 'Hydrostatic pressure and uplift', 'Limit equilibrium: sliding resistance and overturning', 'Meyerhof effective width', 'Vesic bearing capacity with load inclination', 'Passive toe resistance'],
  limits: 'Rigid wall per metre run on a homogeneous drained foundation, lumped (global) safety factors, static loads only; the governing seaward case with low water in front and lagging groundwater behind is analysed. Deep-seated slip, settlement, seepage erosion and seismic thrust are covered by other modules; finite-element soil–structure analysis is outside the module.',
  publishes: ['marine.wallBaseWidth', 'marine.wallVertical', 'marine.wallHorizontal', 'marine.wallQmax'],
  searchTerms: 'gravity seawall sliding overturning bearing capacity retaining wall tidal lag',
  inputs: [
    { key: 'Hw', label: 'Wall height (base to cope)', unit: 'm', type: 'number', default: 8, min: 1, max: 30, group: 'Wall', role: 'geometry' },
    { key: 'B', label: 'Base width', unit: 'm', type: 'number', default: 6.5, min: 0.5, max: 40, group: 'Wall', role: 'geometry', design: true },
    { key: 'bt', label: 'Top width', unit: 'm', type: 'number', default: 2, min: 0.3, max: 40, group: 'Wall', role: 'geometry', design: true, help: 'Vertical seaward face; the back face batters from the top width to the base width.' },
    { key: 'conc', label: 'Wall concrete', type: 'select', source: 'materials:concrete', default: 'c35', group: 'Wall', role: 'material', help: 'Gives the mass-concrete density.' },
    { key: 'length', label: 'Length of wall', unit: 'm', type: 'number', default: 300, min: 1, max: 20000, group: 'Wall', role: 'geometry' },
    { key: 'phiDeg', label: 'Backfill friction angle φ′', unit: 'deg', type: 'number', default: 35, min: 20, max: 45, group: 'Backfill', role: 'soil', uncertain: { dist: 'normal', cov: 0.07 } },
    { key: 'gam', label: 'Backfill bulk unit weight', unit: 'kN/m³', type: 'number', default: 18, min: 12, max: 23, group: 'Backfill', role: 'soil' },
    { key: 'gamSat', label: 'Backfill saturated unit weight', unit: 'kN/m³', type: 'number', default: 20, min: 14, max: 24, group: 'Backfill', role: 'soil' },
    { key: 'deltaW', label: 'Wall friction angle δ', unit: 'deg', type: 'number', default: 0, min: 0, max: 30, group: 'Backfill', role: 'soil', help: '0 gives Rankine; up to ⅔ φ′ for rough concrete (Coulomb).' },
    { key: 'q', label: 'Surcharge behind the wall', unit: 'kPa', type: 'number', default: 15, min: 0, max: 200, group: 'Loads', role: 'bc', uncertain: { dist: 'lognormal', cov: 0.2 }, help: 'Traffic, storage or crane surcharge on the backfill.' },
    { key: 'Hx', label: 'Additional seaward line load', unit: 'kN/m', type: 'number', default: 20, min: 0, max: 2000, group: 'Loads', role: 'bc', help: 'Bollard pull, wave-trough suction beyond the stated front level, or debris load.' },
    { key: 'yx', label: 'Height of the line load above the base', unit: 'm', type: 'number', default: 8, min: 0, max: 40, group: 'Loads', role: 'bc' },
    { key: 'hs', label: 'Sea level above the base (low water / wave trough)', unit: 'm', type: 'number', default: 3, min: 0, max: 30, group: 'Water levels', role: 'bc', help: 'The low seaward level gives the governing seaward case.' },
    { key: 'lag', label: 'Groundwater lag above the sea level', unit: 'm', type: 'number', default: 1.5, min: 0, max: 10, group: 'Water levels', role: 'ic', help: 'Initial pore-water state of the backfill: tidal lag plus overtopping recharge. Good drainage keeps it below about ⅓ of the tidal range.' },
    ...slrInputs({ group: 'Water levels' }),
    { key: 'gamW', label: 'Unit weight of sea water', unit: 'kN/m³', type: 'number', default: 10.05, min: 9.7, max: 10.3, group: 'Water levels', role: 'environment' },
    { key: 'tanD', label: 'Base friction tan δb', unit: '-', type: 'number', default: 0.6, min: 0.2, max: 0.9, group: 'Foundation', role: 'bc', uncertain: { dist: 'normal', cov: 0.1 }, help: 'Cast-in-situ concrete on rubble or dense gravel ≈ 0.6; precast on levelled stone 0.5.' },
    { key: 'phiF', label: 'Foundation friction angle', unit: 'deg', type: 'number', default: 36, min: 0, max: 45, group: 'Foundation', role: 'soil', bind: 'site.phiDeg', uncertain: { dist: 'normal', cov: 0.07 } },
    { key: 'cF', label: 'Foundation cohesion', unit: 'kPa', type: 'number', default: 0, min: 0, max: 500, group: 'Foundation', role: 'soil' },
    { key: 'gamF', label: 'Foundation submerged unit weight', unit: 'kN/m³', type: 'number', default: 10, min: 4, max: 13, group: 'Foundation', role: 'soil' },
    { key: 'D', label: 'Embedment of the base at the toe', unit: 'm', type: 'number', default: 1.0, min: 0, max: 8, group: 'Foundation', role: 'geometry', design: true, help: 'Reduce by the predicted toe scour depth.' },
    { key: 'passive', label: 'Share of toe passive resistance relied on', unit: '-', type: 'number', default: 0, min: 0, max: 1, group: 'Foundation', role: 'criteria', help: '0 where scour or excavation can remove the toe soil.' },
    { key: 'sfS', label: 'Required sliding factor', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Criteria', role: 'criteria' },
    { key: 'sfO', label: 'Required overturning factor', unit: '-', type: 'number', default: 2.0, min: 1, max: 4, group: 'Criteria', role: 'criteria' },
    { key: 'sfB', label: 'Required bearing factor', unit: '-', type: 'number', default: 2.5, min: 1, max: 5, group: 'Criteria', role: 'criteria' },
  ],
  compute(inp, ctx) {
    need(inp.bt <= inp.B + 1e-9, 'The top width cannot exceed the base width.');
    need(inp.gamSat >= inp.gam, 'Saturated unit weight must be at least the bulk unit weight.');
    const mat = ctx.material(inp.conc);
    need(mat && mat.rho > 0, 'The selected concrete has no density in the materials library.');
    const gc = mat.rho * G, gw = inp.gamW * 1e3, g1 = inp.gam * 1e3, g2 = inp.gamSat * 1e3, phi = inp.phiDeg * DEG, dW = inp.deltaW * DEG;
    need(dW <= phi, 'Wall friction δ cannot exceed the backfill friction angle.');
    const Ka = dW > 0 ? coulombKa({ phi, delta: dW }) : rankine(phi).Ka;

    /** Stability for a base width B and a sea-level rise (both water levels move up together). */
    const state = (B, rise) => {
      const bt = Math.min(inp.bt, B), Hw = inp.Hw;
      const hs = Math.min(inp.hs + rise, Hw), hg = Math.min(hs + inp.lag, Hw);
      // wall: rectangle bt wide at the seaward face + triangle behind it; arms from the seaward toe
      const W1 = gc * bt * Hw, W2 = gc * 0.5 * (B - bt) * Hw;
      // backfill wedge resting on the battered back face (total unit weights; water enters through uplift)
      const k = (B - bt) / Hw;
      const aW = (y) => k * y * y / 2, mW = (y) => k * B * y * y / 2 - k * k * y * y * y / 6;
      const Vs = g2 * aW(hg) + g1 * (aW(Hw) - aW(hg)), Ms = g2 * mW(hg) + g1 * (mW(Hw) - mW(hg));
      const Vq = inp.q * 1e3 * (B - bt);
      const pd = pressureDiagram({ layers: [{ thickness: Hw, gamma: g1, gammaSat: g2, phi, c: 0, K: Ka * Math.cos(dW) }], state: 'active', waterTableDepth: Hw - hg, surcharge: inp.q * 1e3, gammaW: gw });
      // arm of the effective thrust alone: remove the water triangle from the total resultant
      const yEff = pd.Peff > 0 ? (pd.P * pd.armAboveBase - pd.Pwater * hg / 3) / pd.Peff : Hw / 3;
      const Pv = pd.Peff * Math.tan(dW);
      const Fsea = 0.5 * gw * hs * hs;
      const U = gw * 0.5 * (hs + hg) * B, xU = hs + hg > 0 ? B * (hs + 2 * hg) / (3 * (hs + hg)) : B / 2;
      const Kp = rankine(inp.phiF * DEG).Kp;
      const Pp = inp.passive * (0.5 * Kp * inp.gamF * 1e3 * inp.D * inp.D + 2 * inp.cF * 1e3 * Math.sqrt(Kp) * inp.D);
      const vertical = [{ V: W1, x: bt / 2 }, { V: W2, x: bt + (B - bt) / 3 }, { V: Vs, x: Vs > 0 ? Ms / Vs : B }, { V: Vq, x: (B + bt) / 2 }, { V: Pv, x: B }];
      const horizontal = [{ H: pd.Peff, y: yEff }, { H: pd.Pwater, y: hg / 3 }, { H: -Fsea, y: hs / 3 }, { H: inp.Hx * 1e3, y: inp.yx }];
      const st = gravityWallStability({ B, vertical, horizontal, U, xU, deltaBase: Math.atan(inp.tanD), cBase: 0, Pp });
      let qult = 0;
      const e = Math.abs(st.e);
      if (e < B / 2 && st.H < st.N * Math.tan(inp.phiF * DEG) + inp.cF * 1e3 * B + 1) {
        try { qult = bearingCapacity({ method: 'vesic', B, D: inp.D, c: inp.cF * 1e3, phi: inp.phiF * DEG, gamma: inp.gamF * 1e3, V: st.N, HB: Math.max(st.H, 0), eB: e }).qult; } catch { qult = 0; }
      }
      const qM = Number.isFinite(st.qMeyerhof) ? st.qMeyerhof : 1e12;
      return { st, pd, U, Fsea, Pv, W: W1 + W2, Vs, hs, hg, qult, fsB: qult > 0 ? qult / qM : 0, qM, yEff, area: bt * Hw + 0.5 * (B - bt) * Hw };
    };
    const s = run(() => state(inp.B, inp.slr), 'Wall stability (the net vertical load must stay positive — check water levels and wall weight)');
    const { st, pd } = s;
    finite(st.fsOverturning === Infinity ? 1 : st.fsOverturning, 'Overturning factor');
    const fsS = Math.min(st.fsSliding, 999), fsO = Math.min(st.fsOverturning, 999);
    const margin = (B) => { try { const q = state(B, inp.slr); return Math.min(q.st.fsSliding / inp.sfS, q.st.fsOverturning / inp.sfO, q.fsB / inp.sfB); } catch { return 0; } };
    let lo = Math.max(inp.bt, 0.5), hi = 60;
    if (margin(lo) >= 1) hi = lo; else for (let i = 0; i < 50; i++) { const mid = 0.5 * (lo + hi); if (margin(mid) < 1) lo = mid; else hi = mid; }
    const Breq = hi;
    const ratios = { sliding: fsS / inp.sfS, overturning: fsO / inp.sfO, bearing: s.fsB / inp.sfB };
    const gov = Object.entries(ratios).sort((a, b) => a[1] - b[1])[0][0];
    const ok = Math.min(...Object.values(ratios)) >= 1;
    const eccRatio = Math.abs(st.e) / inp.B;
    const scen = slrRows(inp).map((sc) => { try { const q = state(inp.B, sc.rise); return { name: sc.name, rise: sc.rise, U: r(q.U / 1e3), N: r(q.st.N / 1e3), fsS: r(Math.min(q.st.fsSliding, 999), 3), fsO: r(Math.min(q.st.fsOverturning, 999), 3), fsB: r(q.fsB, 3) }; } catch { return { name: sc.name, rise: sc.rise, U: 0, N: 0, fsS: 0, fsO: 0, fsB: 0 }; } });
    const bs = linspace(Math.max(inp.bt, 0.5 * Math.min(inp.B, Breq)), 1.6 * Math.max(inp.B, Breq), 30);
    const sweep = bs.map((b) => { try { const q = state(b, inp.slr); return [Math.min(q.st.fsSliding, 10), Math.min(q.st.fsOverturning, 10), Math.min(q.fsB, 10)]; } catch { return [0, 0, 0]; } });
    const conc = s.area * inp.length, dConc = (Breq - inp.B) * 0.5 * inp.Hw * inp.length;
    const qMax = Number.isFinite(st.qMax) ? st.qMax : 1e12;

    const advice = [];
    advice.push(ok
      ? eng('ok', `The ${f(inp.B)} m base meets all three margins: sliding ${f(fsS)}, overturning ${f(fsO)}, bearing ${f(s.fsB)} (${gov} is closest to its limit).`,
        `The minimum base is ${f(Breq)} m. Keep the reserve unless toe scour, groundwater lag (${f(inp.lag)} m assumed) and backfill strength are confirmed by investigation and monitoring.`,
        'EN 1997-1:2004 §9.7 with lumped factors per BS 6349-2:2019 / USACE EM 1110-2-2502 (confirm the factor format of the governing code)',
        `Trimming to the minimum base would save about ${f(Math.max(0, -dConc), 3)} m³ of concrete (≈ ${f(Math.max(0, -dConc) * CO2_CONC, 3)} t CO₂e, screening factor) over ${f(inp.length)} m.`)
      : eng('fail', `The ${f(inp.B)} m base fails in ${gov}: sliding ${f(fsS)} (need ${f(inp.sfS)}), overturning ${f(fsO)} (need ${f(inp.sfO)}), bearing ${f(s.fsB)} (need ${f(inp.sfB)}).`,
        `Widen the base to about ${f(Breq)} m, or cut the groundwater lag with weep holes and a free-draining backfill, add a shear key or deeper toe, or use lighter (rock) backfill behind the wall.`,
        'EN 1997-1:2004 §9.7; BS 6349-2:2019; USACE EM 1110-2-2502',
        `Widening adds about ${f(Math.max(0, dConc), 3)} m³ of concrete (≈ ${f(Math.max(0, dConc) * CO2_CONC, 3)} t CO₂e, screening factor); drainage that halves the lag often achieves the same margin with almost no material.`));
    if (eccRatio > 1 / 6) advice.push(eng('warn', `The base resultant lies outside the middle third (e/B = ${f(eccRatio, 2)}): the heel lifts off and the toe pressure rises to ${f(qMax / 1e3)} kPa.`, 'Widen the base or move weight towards the heel; keep the resultant inside the middle third for normal loading.', 'USACE EM 1110-2-2502 §4; BS 6349-2:2019'));
    const worst = scen[scen.length - 1];
    advice.push(eng(worst.fsS < inp.sfS || worst.fsB < inp.sfB ? 'warn' : 'info', `Sea-level rise raises sea and groundwater together: at +${f(worst.rise)} m (${worst.name}) uplift is ${f(worst.U)} kN/m, sliding factor ${f(worst.fsS)}, bearing factor ${f(worst.fsB)}.`,
      'Allow for the end-of-life groundwater level in the design, and detail drainage outlets (flap valves) so that they still discharge at the higher sea level.', SLR_BASIS));
    advice.push(eco('info', `The wall section is ${f(s.area, 3)} m² (${f(conc, 3)} m³ of concrete over ${f(inp.length)} m).`,
      `Each 0.5 m of base width adds ${f(0.25 * inp.Hw * inp.length, 3)} m³ of concrete; compare that with drainage measures, a shear key and geotechnical investigation that would justify a higher friction angle.`,
      'AACE RP 18R-97 quantity-based estimating; value-of-information for ground investigation',
      `Embodied carbon of the wall ≈ ${f(conc * CO2_CONC, 3)} t CO₂e (screening factor 0.3 t/m³); mass concrete tolerates high cement replacement, which cuts this substantially.`));

    return {
      kpis: [
        { key: 'fsSliding', label: 'Sliding safety factor', value: r(fsS), unit: '-', limit: inp.sfS, sense: 'min', explain: 'Base friction (plus any passive resistance relied on) over the net seaward horizontal load. Below 1 the wall translates seaward.' },
        { key: 'fsOverturning', label: 'Overturning safety factor', value: r(fsO), unit: '-', limit: inp.sfO, sense: 'min', explain: 'Restoring over overturning moments about the seaward toe, uplift counted as overturning. Below 1 the wall rotates seaward.' },
        { key: 'fsBearing', label: 'Bearing safety factor', value: r(s.fsB), unit: '-', limit: inp.sfB, sense: 'min', explain: 'Ultimate bearing pressure under the inclined, eccentric load over the Meyerhof contact pressure. Low values mean the toe punches into the foundation.' },
        { key: 'qMax', label: 'Maximum toe pressure', value: r(Math.min(qMax, 1e9) / 1e3), unit: 'kPa', explain: 'Peak effective contact pressure at the seaward toe from the linear (or triangular) base distribution; compare with what the bedding and foundation can take.' },
        { key: 'ecc', label: 'Eccentricity ratio e / B', value: r(eccRatio, 3), unit: '-', limit: 0.1667, sense: 'max', explain: 'Offset of the base resultant from the centre as a share of the width. Up to 1/6 the whole base is in compression; beyond it the heel lifts.' },
        { key: 'Pa', label: 'Effective active earth thrust', value: r(pd.Peff / 1e3), unit: 'kN/m', explain: 'Horizontal thrust of the backfill skeleton on the vertical plane through the heel, including surcharge.' },
        { key: 'Pw', label: 'Net water thrust (behind − in front)', value: r((pd.Pwater - s.Fsea) / 1e3), unit: 'kN/m', explain: 'Unbalanced hydrostatic force from groundwater standing above the sea level. It is created entirely by the lag and is the cheapest load to remove by drainage.' },
        { key: 'U', label: 'Uplift on the base', value: r(s.U / 1e3), unit: 'kN/m', explain: 'Water pressure under the base, varying from the sea level at the toe to the groundwater level at the heel. It reduces the weight available for friction.' },
        { key: 'N', label: 'Effective vertical load', value: r(st.N / 1e3), unit: 'kN/m', explain: 'Total weight of wall and fill on the base less uplift; the normal force that generates friction and bearing pressure.' },
        { key: 'H', label: 'Net horizontal load', value: r(st.H / 1e3), unit: 'kN/m', explain: 'Sum of earth thrust, unbalanced water thrust and the line load, acting seaward.' },
        { key: 'Breq', label: 'Base width for the required margins', value: r(Breq), unit: 'm', explain: 'Smallest base that meets the sliding, overturning and bearing factors in the design case; the strengthening target or the optimisation limit.' },
      ],
      checks: [
        sfCheck('Sliding on the base', st.fsSliding, inp.sfS, 'EN 1997-1 §6.5.3; USACE EM 1110-2-2502 §4', 'Guards against the wall being pushed seaward along its base. Movement opens the backfill, drops the cope and breaks services.'),
        sfCheck('Overturning about the toe', st.fsOverturning, inp.sfO, 'USACE EM 1110-2-2502 §4; BS 6349-2', 'Guards against rotation about the seaward toe; in soil foundations a bearing failure under the toe normally precedes it.'),
        sfCheck('Bearing capacity of the foundation', s.fsB > 0 ? s.fsB : 1e-6, inp.sfB, 'EN 1997-1 Annex D; Vesic (1975)', 'Guards against a shear failure of the foundation under the inclined, eccentric base load, which tilts the wall seaward.'),
        check('Resultant within the middle third', eccRatio, 1 / 6, '-', 'USACE EM 1110-2-2502; BS 6349-2', 'Keeps the whole base in compression. Outside the middle third the heel separates, water pressure enters the gap and toe pressure escalates.'),
      ],
      plots: [
        { type: 'line', title: 'Pressures behind the wall', x: { label: 'Pressure', unit: 'kPa' }, y: { label: 'Depth below cope', unit: 'm' },
          series: [{ name: 'Effective active pressure', x: pd.points.map((p) => r(p.sigmaH / 1e3, 3)), y: pd.points.map((p) => -r(p.z, 3)) }, { name: 'Pore-water pressure', x: pd.points.map((p) => r(p.u / 1e3, 3)), y: pd.points.map((p) => -r(p.z, 3)) }],
          explain: 'Active earth pressure and groundwater pressure on the vertical plane through the heel. The water triangle below the water table often rivals the earth pressure — the case for drainage.' },
        { type: 'line', title: 'Safety factors versus base width', x: { label: 'Base width', unit: 'm' }, y: { label: 'Achieved / required factor', unit: '-' },
          series: [{ name: 'Sliding', x: bs.map((b) => r(b, 3)), y: sweep.map((q) => r(q[0] / inp.sfS, 3)) }, { name: 'Overturning', x: bs.map((b) => r(b, 3)), y: sweep.map((q) => r(q[1] / inp.sfO, 3)) }, { name: 'Bearing', x: bs.map((b) => r(b, 3)), y: sweep.map((q) => r(q[2] / inp.sfB, 3)) }, { name: 'Limit', x: [bs[0], bs[29]], y: [1, 1], dash: true }],
          explain: 'Each margin normalised by its required factor. The mode whose curve crosses 1.0 last governs the width; bearing usually rises fastest because widening both spreads and centres the load.' },
        { type: 'line', title: 'Base contact pressure', x: { label: 'Distance from seaward toe', unit: 'm' }, y: { label: 'Effective pressure', unit: 'kPa' },
          series: [{ name: 'Contact pressure', x: eccRatio > 1 / 6 ? [0, r(3 * (inp.B / 2 - Math.abs(st.e)), 3), inp.B] : [0, inp.B], y: eccRatio > 1 / 6 ? [r(Math.min(qMax, 1e9) / 1e3, 3), 0, 0] : [r(st.qMax / 1e3, 3), r(st.qMin / 1e3, 3)] }],
          explain: 'Linear distribution of effective contact pressure under the base. A triangle that ends short of the heel shows loss of contact.' },
      ],
      tables: [
        { title: 'Load summary per metre run', columns: [{ key: 'l', label: 'Load' }, { key: 'v', label: 'Force', unit: 'kN/m' }, { key: 'd', label: 'Direction' }],
          rows: [['Wall self-weight', s.W, 'down'], ['Backfill on the back face', s.Vs, 'down'], ['Surcharge on the back face', inp.q * 1e3 * (inp.B - Math.min(inp.bt, inp.B)), 'down'], ['Wall-friction component', s.Pv, 'down'], ['Uplift', s.U, 'up'], ['Effective earth thrust', pd.Peff, 'seaward'], ['Groundwater thrust', pd.Pwater, 'seaward'], ['Sea-water thrust', s.Fsea, 'landward'], ['Line load', inp.Hx * 1e3, 'seaward']].map(([l, v, d]) => ({ l, v: r(v / 1e3), d })),
          explain: 'Every force in the equilibrium, for checking against an independent hand calculation and for seeing which load dominates.' },
        { title: 'Sea-level-rise scenarios', columns: [{ key: 'name', label: 'Scenario' }, { key: 'rise', label: 'Rise', unit: 'm' }, { key: 'U', label: 'Uplift', unit: 'kN/m' }, { key: 'N', label: 'Effective vertical', unit: 'kN/m' }, { key: 'fsS', label: 'Sliding', unit: '-' }, { key: 'fsO', label: 'Overturning', unit: '-' }, { key: 'fsB', label: 'Bearing', unit: '-' }],
          rows: scen, explain: 'The same wall with sea and groundwater raised by each scenario. Rising water adds uplift and removes friction; zero rows mean the wall is no longer held down.' },
      ],
      notes: [`${gov[0].toUpperCase()}${gov.slice(1)} governs. Ka = ${f(Ka)}; sea ${f(s.hs)} m and groundwater ${f(s.hg)} m above the base including ${f(inp.slr)} m of sea-level rise.`,
        'The wave-crest (landward) case is resisted by the backfill and is checked with the wave-pressure module.'],
      advice,
      publish: { 'marine.wallBaseWidth': inp.B, 'marine.wallVertical': st.N, 'marine.wallHorizontal': st.H, 'marine.wallQmax': Math.min(qMax, 1e9) },
      quantities: { concrete_m3: r(conc) },
    };
  },
  verify: [
    { name: 'Dry rectangular gravity wall with Rankine backfill, statics by hand',
      source: 'Hand calculation (Craig, Soil Mechanics, ch. 11 method): wall 3 × 6 m, ρc = 2400 → W = 423.65 kN/m; φ′ = 30°, γ = 18 kN/m³ → Ka = ⅓, Pa = ½ × ⅓ × 18 × 36 = 108 kN/m at 2 m; sliding = 0.6 × 423.65/108 = 2.354; overturning = 423.65 × 1.5/(108 × 2) = 2.942; e = 1.5 − 419.47/423.65 = 0.5099 m > B/6 → qmax = 2N/[3(B/2 − e)] = 285.2 kPa',
      inputs: { Hw: 6, B: 3, bt: 3, conc: 'c35', phiDeg: 30, gam: 18, gamSat: 20, deltaW: 0, q: 0, Hx: 0, hs: 0, lag: 0, slr: 0, tanD: 0.6, D: 0, passive: 0 },
      expect: [{ kpi: 'fsSliding', value: 2.354, tol: 0.002 }, { kpi: 'fsOverturning', value: 2.942, tol: 0.002 }, { kpi: 'qMax', value: 285.2, tol: 0.003 }, { kpi: 'Pa', value: 108, tol: 0.001 }] },
  ],
};
