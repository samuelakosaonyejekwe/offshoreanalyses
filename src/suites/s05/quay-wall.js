// Quay-wall earth pressures and retaining-structure stability: an anchored sheet-pile (or combi) wall by the
// free-earth support method, or a gravity (block or caisson) wall by statics about its toe, under backfill, cargo
// surcharge, tidal lag, bollard pull and pseudo-static earthquake loading, with the permissible cargo surcharge.
import { rankine, coulombKa, coulombKp, pressureDiagram, mononobeOkabe, gravityWallStability, bollardPull, sheetPileWall, deadmanAnchor } from '../../engine/geotech/earthpressure.js';
import { bearingCapacity } from '../../engine/geotech/bearing.js';
import { DEG, need, run, r, f, cap, check, sfCheck, eng, eco, linspace, maxWhere, steelCo2, concCo2, CONFIRM, ECON_BASIS } from './_common.js';

const GW = 10055; // unit weight of sea water, N/m³ (1025 kg/m³)

/** Earth-pressure coefficients (horizontal components) for the static and pseudo-static cases. */
function coefficients(inp, kh) {
  const phi = inp.phi * DEG, da = inp.deltaA * DEG, dp = inp.deltaP * DEG;
  if (kh > 0) {
    need(phi - Math.atan(kh) > 0.02, `A seismic coefficient of ${f(kh)} exceeds what a backfill of φ′ = ${f(inp.phi)}° can sustain (no active wedge solution): improve the fill or reduce the coefficient.`);
    const a = mononobeOkabe({ phi, delta: da, kh }), p = mononobeOkabe({ phi, delta: dp, kh });
    return { Ka: a.KAE * Math.cos(da), Kp: p.KPE * Math.cos(dp), KaTot: a.KAE };
  }
  const Ka = coulombKa({ phi, delta: da }), Kp = coulombKp({ phi, delta: dp });
  return { Ka: Ka * Math.cos(da), Kp: Kp * Math.cos(dp), KaTot: Ka };
}

export default {
  id: 's05-quay-wall', suite: 5, title: 'Quay-wall earth pressure and stability',
  summary: 'Earth and water pressures on a quay wall and its stability: embedment, anchor force and bending moment of an anchored sheet-pile wall, or sliding, overturning and bearing of a gravity wall, in service and in the design earthquake, with the permissible cargo surcharge behind the cope.',
  covers: [5, 12], fidelity: 'analytical',
  standards: ['EN 1997-1:2004 §9 (retaining structures), Annex C', 'BS 6349-2:2019 §6–7 (quay walls)', 'EAU 2012 Recommendations of the Committee for Waterfront Structures', 'EN 1998-5:2004 §7, Annex E (Mononobe–Okabe)', 'PIANC WG34 (2001) seismic design of port structures', 'EN 1993-5:2007 (steel piling)'],
  models: ['Rankine earth pressure theory', 'Coulomb earth pressure theory with wall friction', 'Mononobe–Okabe pseudo-static earth pressure', 'Westergaard hydrodynamic pressure', 'Free-earth support method for anchored walls', 'Static force and moment equilibrium of a gravity wall', 'Meyerhof bearing capacity with eccentric and inclined load', 'Deadman anchor resistance'],
  limits: 'Limit-equilibrium design in cohesionless backfill with a single set of soil parameters; layered or cohesive soils, wall flexibility (beyond a moment-reduction factor), staged construction, anchor prestress and deep-seated slip surfaces are not modelled. The sheet-pile solution takes equal water levels on both sides and adds the tidal lag as a separate check of the anchor; earthquake loading is pseudo-static. Soil–structure finite-element analysis (Mohr–Coulomb or Cam-Clay continuum) is outside the module.',
  publishes: ['wharf.quayAnchorForce', 'wharf.quaySurchargePermissible', 'wharf.quayThrust'],
  searchTerms: 'quay wall sheet pile anchored combi wall gravity wall earth pressure Mononobe-Okabe',
  inputs: [
    { key: 'wallType', label: 'Wall type', type: 'select', default: 'sheet', group: 'Wall', role: 'geometry', options: [{ value: 'sheet', label: 'Anchored sheet-pile or combi wall' }, { value: 'gravity', label: 'Gravity wall (blocks, caisson, mass concrete)' }] },
    { key: 'H', label: 'Retained height, cope to dredged level', unit: 'm', type: 'number', default: 14, min: 2, max: 35, group: 'Wall', role: 'geometry', design: true, help: 'Include the over-dredging tolerance and the scour allowance at the toe.' },
    { key: 'zw', label: 'Water level below the cope', unit: 'm', type: 'number', default: 3.5, min: 0, max: 35, group: 'Wall', role: 'bc', help: 'Low water in front of the wall. Set equal to the retained height for a wall in the dry.' },
    { key: 'lag', label: 'Tidal lag (ground water above the harbour level)', unit: 'm', type: 'number', default: 0.8, min: 0, max: 5, group: 'Wall', role: 'bc', help: 'Excess head behind the wall at low water; 0.5–1.5 m with weep holes and a free-draining fill, more without drainage.' },
    { key: 'gamma', label: 'Unit weight of the backfill above water', unit: 'kN/m³', type: 'number', default: 18, min: 12, max: 23, group: 'Backfill', role: 'soil' },
    { key: 'gammaSub', label: 'Submerged unit weight of the backfill', unit: 'kN/m³', type: 'number', default: 10, min: 5, max: 13, group: 'Backfill', role: 'soil', bind: 'site.gammaSub' },
    { key: 'phi', label: 'Friction angle of the backfill and founding soil', unit: 'deg', type: 'number', default: 35, min: 20, max: 45, group: 'Backfill', role: 'soil', bind: 'site.phiDeg', uncertain: { dist: 'normal', cov: 0.07 } },
    { key: 'deltaA', label: 'Wall friction on the active side', unit: 'deg', type: 'number', default: 20, min: 0, max: 35, group: 'Backfill', role: 'soil', help: 'Up to two thirds of φ′ for steel sheet piles and rough concrete; 0 gives the Rankine value.' },
    { key: 'deltaP', label: 'Wall friction on the passive side', unit: 'deg', type: 'number', default: 12, min: 0, max: 25, group: 'Backfill', role: 'soil', help: 'Kept to about one third of φ′: the plane-wedge Coulomb passive coefficient is unsafe at higher wall friction.' },
    { key: 'q', label: 'Cargo surcharge behind the cope', unit: 'kPa', type: 'number', default: 40, min: 0, max: 300, group: 'Actions', role: 'bc', design: true, uncertain: { dist: 'gumbel', cov: 0.2 }, help: 'Uniform storage and traffic load on the apron and yard: 20 kPa general cargo, 40–60 kPa containers, 100+ kPa bulk.' },
    { key: 'bollard', label: 'Bollard pull per metre of wall', unit: 'kN/m', type: 'number', default: 40, min: 0, max: 500, group: 'Actions', role: 'bc', help: 'Bollard load divided by the wall length that shares it (anchor wall panel or monolith length), acting at cope level.' },
    { key: 'kh', label: 'Horizontal seismic coefficient', unit: '-', type: 'number', default: 0.1, min: 0, max: 0.5, group: 'Actions', role: 'bc', help: 'Pseudo-static coefficient, typically half to two thirds of the peak ground acceleration in g for walls that can move a little (EN 1998-5 §7.3.2.2, r = 1.5–2).' },
    { key: 'anchorZ', label: 'Anchor level below the cope', unit: 'm', type: 'number', default: 2, min: 0.2, max: 15, group: 'Sheet-pile wall', role: 'geometry', showIf: { key: 'wallType', in: ['sheet'] } },
    { key: 'dProv', label: 'Embedment provided below dredged level', unit: 'm', type: 'number', default: 9, min: 0.5, max: 40, group: 'Sheet-pile wall', role: 'geometry', design: true, showIf: { key: 'wallType', in: ['sheet'] } },
    { key: 'Wel', label: 'Section modulus of the wall', unit: 'cm³/m', type: 'number', default: 9000, min: 300, max: 60000, group: 'Sheet-pile wall', role: 'geometry', design: true, showIf: { key: 'wallType', in: ['sheet'] }, help: 'About 1500–5000 cm³/m for Z-section sheet piles, 5000–20000 cm³/m for combined tube or H-pile walls. Deduct corrosion.' },
    { key: 'wallMass', label: 'Mass of the wall per square metre', unit: 'kg/m²', type: 'number', default: 280, min: 50, max: 800, group: 'Sheet-pile wall', role: 'geometry', showIf: { key: 'wallType', in: ['sheet'] }, help: 'From the section table, for the quantity estimate: about 100–190 kg/m² for Z-section sheet piles, 200–400 kg/m² for combined walls.' },
    { key: 'steel', label: 'Wall steel grade', type: 'select', source: 'materials:steel', default: 's355', group: 'Sheet-pile wall', role: 'material', showIf: { key: 'wallType', in: ['sheet'] } },
    { key: 'rowe', label: 'Moment reduction factor for wall flexibility', unit: '-', type: 'number', default: 1, min: 0.4, max: 1, group: 'Sheet-pile wall', role: 'criteria', showIf: { key: 'wallType', in: ['sheet'] }, help: "Rowe's chart value for flexible walls in dense sand; 1.0 for stiff combined walls and as a safe default." },
    { key: 'tieSpacing', label: 'Tie-rod spacing', unit: 'm', type: 'number', default: 2.5, min: 0.5, max: 8, group: 'Sheet-pile wall', role: 'geometry', showIf: { key: 'wallType', in: ['sheet'] } },
    { key: 'tieCap', label: 'Design resistance of one tie rod', unit: 'kN', type: 'number', default: 1800, min: 50, max: 20000, group: 'Sheet-pile wall', role: 'criteria', design: true, showIf: { key: 'wallType', in: ['sheet'] }, help: 'Thread area × design strength, with corrosion allowance; for example an M100 upset rod in S460 is about 2000 kN.' },
    { key: 'deadH', label: 'Depth to the base of the anchor wall', unit: 'm', type: 'number', default: 5.5, min: 1, max: 15, group: 'Sheet-pile wall', role: 'geometry', showIf: { key: 'wallType', in: ['sheet'] }, help: 'A continuous deadman (anchor wall) placed beyond the active wedge of the quay wall.' },
    { key: 'B', label: 'Base width of the gravity wall', unit: 'm', type: 'number', default: 13, min: 1, max: 40, group: 'Gravity wall', role: 'geometry', design: true, showIf: { key: 'wallType', in: ['gravity'] } },
    { key: 'gammaWall', label: 'Unit weight of the wall (average, with fill in caissons)', unit: 'kN/m³', type: 'number', default: 22, min: 15, max: 26, group: 'Gravity wall', role: 'material', showIf: { key: 'wallType', in: ['gravity'] } },
    { key: 'deltaB', label: 'Base friction angle', unit: 'deg', type: 'number', default: 30, min: 10, max: 45, group: 'Gravity wall', role: 'soil', showIf: { key: 'wallType', in: ['gravity'] }, help: 'Concrete on a rubble bed: about 30–35° (friction coefficient 0.6–0.7); precast blocks on blocks 26–30°.' },
    { key: 'fsS', label: 'Required factor of safety, sliding and passive resistance', unit: '-', type: 'number', default: 1.5, min: 1, max: 3, group: 'Acceptance', role: 'criteria', help: `Applied to sliding of a gravity wall and to the passive resistance of a sheet-pile wall. 1.5 static; the seismic value below applies in the earthquake case. ${CONFIRM}` },
    { key: 'fsO', label: 'Required factor of safety, overturning', unit: '-', type: 'number', default: 2.0, min: 1, max: 3, group: 'Acceptance', role: 'criteria', help: CONFIRM },
    { key: 'fsB', label: 'Required factor of safety, bearing', unit: '-', type: 'number', default: 2.5, min: 1, max: 4, group: 'Acceptance', role: 'criteria', help: CONFIRM },
    { key: 'fsE', label: 'Required factor of safety in the earthquake case', unit: '-', type: 'number', default: 1.1, min: 1, max: 2, group: 'Acceptance', role: 'criteria', help: `PIANC WG34 and EN 1998-5 accept values near 1.0–1.2 for pseudo-static checks. ${CONFIRM}` },
    { key: 'gM', label: 'γM0 on the wall steel', unit: '-', type: 'number', default: 1.0, min: 1, max: 1.3, group: 'Acceptance', role: 'criteria', help: `EN 1993-5 §5.1.1. ${CONFIRM}` },
    { key: 'gF', label: 'Load factor on bending moment and anchor force', unit: '-', type: 'number', default: 1.35, min: 1, max: 1.6, group: 'Acceptance', role: 'criteria', help: `Applied to the effects from the unfactored limit-equilibrium analysis (EN 1997-1 Design Approach 2*). ${CONFIRM}` },
    { key: 'nInt', label: 'Integration intervals along the wall', unit: '-', type: 'int', default: 2000, min: 100, max: 32000, group: 'Numerics', role: 'control' },
  ],
  mesh: { key: 'nInt', label: 'Integration intervals', levels: [250, 500, 1000, 2000], qoi: 'Mmax' },
  compute(inp, ctx) {
    const H = inp.H, zw = Math.min(inp.zw, H), g = inp.gamma * 1e3, gs = inp.gammaSub * 1e3, phi = inp.phi * DEG, sheet = inp.wallType === 'sheet';
    const cS = coefficients(inp, 0), cE = inp.kh > 0 ? coefficients(inp, inp.kh) : cS, rk = rankine(phi);
    // active thrust on the retained height (horizontal component), static and seismic, for reporting and for the gravity wall
    const thrust = (K, q) => run(() => pressureDiagram({ layers: [{ thickness: Math.max(zw, 1e-6), gamma: g, gammaSat: g, phi, K }, ...(H > zw ? [{ thickness: H - zw, gamma: gs, gammaSat: gs, phi, K }] : [])], state: 'active', surcharge: q, tensionCrack: false }), 'Earth pressure');
    const pS = thrust(cS.Ka, inp.q * 1e3), pE = thrust(cE.Ka, 0.5 * inp.q * 1e3);
    // net water pressure from the tidal lag: triangle over the lag, then constant to the dredged level (and to the toe for the sheet wall)
    const lag = Math.min(inp.lag, zw), hw = H - zw;
    const waterLag = (depthBelow) => ({ P: GW * lag * (lag / 2 + hw + depthBelow), y: (GW * lag * lag / 2 * (hw + depthBelow + lag / 3) + GW * lag * (hw + depthBelow) ** 2 / 2) / Math.max(GW * lag * (lag / 2 + hw + depthBelow), 1e-9) });
    let out;
    if (sheet) {
      need(inp.anchorZ < H, 'The anchor level must be above the dredged level.');
      const m = ctx?.material ? ctx.material(inp.steel) : null; need(m && m.fy > 0, `Steel grade '${inp.steel}' is not in the materials library.`);
      const solve = (c, q, fsP) => run(() => sheetPileWall({ H, gamma: g, gammaSub: gs, waterDepth: zw, Ka: c.Ka, Kp: c.Kp, fsPassive: fsP, surcharge: q, anchorDepth: inp.anchorZ, roweFactor: inp.rowe, n: Math.round(inp.nInt) }), 'Sheet-pile analysis (the wall may be too high for the passive resistance available)');
      const st = solve(cS, inp.q * 1e3, inp.fsS), eq = inp.kh > 0 ? solve(cE, 0.5 * inp.q * 1e3, inp.fsE) : st;
      // tidal lag and bollard pull go to the anchor: lag pressure as a beam between anchor and the centre of passive resistance, bollard directly
      const wl = waterLag(st.d0 / 2), span = H + st.d0 / 2 - inp.anchorZ, Tlag = wl.P * (wl.y) / Math.max(span, 0.1), Mlag = wl.P * span / 8;
      const T = st.anchorForce + Tlag + inp.bollard * 1e3, TE = eq.anchorForce + Tlag;
      const Mst = st.Mmax + Mlag, ME = eq.Mmax + Mlag;
      const MRd = inp.Wel * 1e-6 * m.fy / inp.gM, Td = Math.max(inp.gF * T, TE) * inp.tieSpacing;
      const dead = deadmanAnchor({ gamma: g, H: inp.deadH, Ka: rk.Ka, Kp: rk.Kp, safetyFactor: inp.fsS });
      const dReq = Math.max(st.d, eq.d);
      // permissible surcharge: embedment, wall bending and tie rod all within their limits
      const okQ = (qq) => { try { const s = solve(cS, qq, inp.fsS); return s.d <= inp.dProv && inp.gF * (s.Mmax + Mlag) <= MRd && inp.gF * (s.anchorForce + Tlag + inp.bollard * 1e3) * inp.tieSpacing <= inp.tieCap * 1e3 && inp.gF * (s.anchorForce + Tlag + inp.bollard * 1e3) <= dead.Tult / inp.fsS * inp.gF; } catch { return false; } };
      const qPerm = maxWhere(okQ, 0, 400e3, 24);
      const checks = [
        check('Embedment of the wall', dReq, inp.dProv, 'm', 'EN 1997-1:2004 §9.7.4; EAU 2012 R 55 (free-earth support)', 'Embedment needed for moment equilibrium about the anchor with the passive resistance divided by its factor of safety, against the embedment provided. Too little and the toe kicks out seaward.'),
        check('Bending of the wall', Math.max(inp.gF * Mst, ME) / 1e3, MRd / 1e3, 'kN·m/m', 'EN 1993-5:2007 §5.2.2', 'Largest bending moment between anchor and toe against the elastic resistance of the section. If exceeded the wall yields and bulges below the anchor.'),
        check('Tie rod', Td / 1e3, inp.tieCap, 'kN', 'EN 1993-5:2007 §7.2; EAU 2012 R 20', 'Anchor force per tie (earth, tidal lag and bollard pull, factored) against its design resistance. A broken tie unloads onto its neighbours and can unzip the anchorage.'),
        check('Anchor wall (deadman) resistance', T / 1e3, dead.Tallow / 1e3, 'kN/m', 'EN 1997-1:2004 §9.7; Terzaghi (1943) deadman theory', 'Anchor force against the net passive resistance of the anchor wall over its factor of safety. If exceeded the anchorage moves forward and the quay wall follows.'),
        ...(inp.kh > 0 ? [check('Embedment in the earthquake case', eq.d, inp.dProv, 'm', 'EN 1998-5:2004 §7.3.2, Annex E', 'Embedment needed with Mononobe–Okabe pressures: the active thrust rises and the passive resistance falls together.')] : []),
      ];
      const WelNeed = Math.max(inp.gF * Mst, ME) * inp.gM / m.fy * 1e6;
      out = { checks, st, eq, T, TE, Mst, ME, MRd, dReq, qPerm, dead, WelNeed, Mplot: Math.max(inp.gF * Mst, ME),
        kp: [
          { key: 'd', label: 'Required embedment', value: r(dReq), unit: 'm', limit: inp.dProv, sense: 'max', explain: 'Depth below dredged level needed for equilibrium with the factored passive resistance (larger of the static and earthquake cases), against the embedment provided.' },
          { key: 'T', label: 'Anchor force', value: r(T / 1e3), unit: 'kN/m', explain: `Unfactored force in the anchorage per metre of wall from earth pressure (${f(st.anchorForce / 1e3)}), tidal lag (${f(Tlag / 1e3)}) and bollard pull (${f(inp.bollard)}). It is what the tie rods and the anchor wall must hold.` },
          { key: 'Mmax', label: 'Largest bending moment in the wall', value: r(Mst / 1e3, 8), unit: 'kN·m/m', explain: `Unfactored span moment at about ${f(st.zMmax)} m below the cope, where the shear is zero. With the load factor it is compared with the section resistance.` },
          { key: 'WelNeed', label: 'Section modulus required', value: r(WelNeed), unit: 'cm³/m', limit: inp.Wel, sense: 'max', explain: 'Elastic section modulus that brings the wall bending check to unity, for comparing sheet-pile and combi-wall sections.' },
          { key: 'Tdead', label: 'Anchor wall resistance (allowable)', value: r(dead.Tallow / 1e3), unit: 'kN/m', explain: 'Net passive resistance of a continuous anchor wall of the depth entered, over the factor of safety.' },
          { key: 'TE', label: 'Anchor force in the earthquake case', value: r(TE / 1e3), unit: 'kN/m', explain: 'Anchor force with Mononobe–Okabe pressures and half the cargo surcharge.' },
        ] };
    } else {
      const B = inp.B, gwl = inp.gammaWall * 1e3;
      // wall weight with buoyancy below the harbour water level; vertical component of the active thrust at the heel
      const Wdry = gwl * B * zw, Wsub = (gwl - GW) * B * (H - zw), W = Wdry + Wsub;
      const tanD = Math.tan(inp.deltaA * DEG);
      const stab = (p, c, kh, qShare, fsReqB) => {
        const wl = waterLag(0), bp = bollardPull({ pull: inp.bollard * 1e3, spreadLength: 1, height: H });
        const horizontal = [{ H: p.Peff, y: p.armAboveBase }, { H: wl.P, y: wl.y }, ...(kh > 0 ? [] : [{ H: bp.H, y: H }])];
        if (kh > 0) { // wall inertia and hydrodynamic suction on the seaward face (Westergaard, 7/12·kh·γw·h² at 0.4 h)
          horizontal.push({ H: kh * gwl * B * H, y: H / 2 });
          if (hw > 0) horizontal.push({ H: (7 / 12) * kh * GW * hw * hw, y: 0.4 * hw });
        }
        const vertical = [{ V: W, x: B / 2 }, { V: p.Peff * tanD, x: B }];
        const Vtot = W + p.Peff * tanD, Htot = horizontal.reduce((s, q) => s + q.H, 0);
        const first = gravityWallStability({ B, vertical, horizontal, deltaBase: inp.deltaB * DEG });
        const Beff = Math.max(B - 2 * Math.max(first.e, 0), 0.05 * B);
        const bc = bearingCapacity({ method: 'meyerhof', B: Beff, L: Infinity, D: 0, c: 0, phi, gamma: gs, V: Math.max(Vtot, 1), HB: Math.min(Htot, 0.95 * Vtot * Math.tan(phi)) });
        return { ...gravityWallStability({ B, vertical, horizontal, deltaBase: inp.deltaB * DEG, qult: bc.qult }), qult: bc.qult, Htot, Vtot, qShare, fsReqB };
      };
      const sS = stab(pS, cS, 0, 1), sE = inp.kh > 0 ? stab(pE, cE, inp.kh, 0.5) : sS;
      const okQ = (qq) => { const s = stab(thrust(cS.Ka, qq), cS, 0, 1); return s.fsSliding >= inp.fsS && s.fsOverturning >= inp.fsO && s.fsBearing >= inp.fsB && s.e <= B / 6 + 1e-9; };
      const qPerm = maxWhere(okQ, 0, 400e3, 24);
      const checks = [
        sfCheck('Sliding on the base', sS.fsSliding, inp.fsS, 'EN 1997-1:2004 §9.7.3; BS 6349-2:2019 §6', 'Base friction over the horizontal thrust of earth, surcharge, tidal lag and bollard pull. Below the required factor the wall slides seaward on its bed.'),
        sfCheck('Overturning about the toe', sS.fsOverturning, inp.fsO, 'EN 1997-1:2004 §9.7.2; BS 6349-2:2019 §6', 'Restoring moment of the wall weight over the overturning moment of the thrusts. Low values mean the resultant is near the toe and the wall is close to tipping.'),
        sfCheck('Bearing under the toe', sS.fsBearing, inp.fsB, 'EN 1997-1:2004 §6.5.2, Annex D; Meyerhof (1953)', 'Bearing resistance of the bed on the effective width over the applied pressure, reduced for the inclined resultant. Failure shows as forward tilting with settlement of the toe.'),
        check('Eccentricity of the resultant (middle third)', Math.abs(sS.e), B / 6, 'm', 'EN 1997-1:2004 §6.5.4', 'Distance of the base resultant from the centre against one sixth of the width. Outside the middle third the heel lifts and the toe pressure climbs steeply.'),
        ...(inp.kh > 0 ? [sfCheck('Sliding in the earthquake case', sE.fsSliding, inp.fsE, 'EN 1998-5:2004 §7.3.2; PIANC WG34 (2001) §T7', 'Sliding with Mononobe–Okabe thrust, wall inertia and hydrodynamic suction. Below 1 the wall moves seaward during shaking; the displacement can be estimated with a sliding-block analysis.'), sfCheck('Overturning in the earthquake case', sE.fsOverturning, inp.fsE, 'EN 1998-5:2004 §7.3.2', 'Overturning about the toe with the seismic thrusts acting together.')] : []),
      ];
      // first estimate of the width that meets every factor: sliding scales with B, overturning and bearing roughly with B²
      const Bneed = B * Math.max(inp.fsS / sS.fsSliding, Math.sqrt(inp.fsO / sS.fsOverturning), Math.sqrt(inp.fsB / Math.max(sS.fsBearing, 0.01)), Math.sqrt(Math.abs(sS.e) / (B / 6)), inp.kh > 0 ? inp.fsE / sE.fsSliding : 0);
      out = { checks, sS, sE, qPerm, Bneed, W, Mplot: pS.Peff * pS.armAboveBase,
        kp: [
          { key: 'fsSlide', label: 'Factor of safety against sliding', value: r(cap(sS.fsSliding, 999)), unit: '-', limit: inp.fsS, sense: 'min', explain: 'Base friction over horizontal thrust in service. It is usually the governing check of a gravity quay wall.' },
          { key: 'fsOT', label: 'Factor of safety against overturning', value: r(cap(sS.fsOverturning, 999)), unit: '-', limit: inp.fsO, sense: 'min', explain: 'Restoring over overturning moment about the toe.' },
          { key: 'fsBear', label: 'Factor of safety in bearing', value: r(cap(sS.fsBearing, 999)), unit: '-', limit: inp.fsB, sense: 'min', explain: `Ultimate bearing pressure of the bed (${f(sS.qult / 1e3)} kPa on the effective width) over the applied pressure of ${f(sS.qMeyerhof / 1e3)} kPa.` },
          { key: 'qToe', label: 'Toe pressure', value: r(sS.qMax / 1e3), unit: 'kPa', explain: 'Largest base pressure under the seaward edge from the linear distribution. It sizes the rubble bed and its thickness.' },
          { key: 'ecc', label: 'Eccentricity of the base resultant', value: r(sS.e), unit: 'm', explain: `Offset of the resultant from the centre of the base (seaward positive). The middle third ends at ${f(B / 6)} m.` },
          { key: 'fsSlideE', label: 'Factor of safety against sliding, earthquake', value: r(cap(sE.fsSliding, 999)), unit: '-', limit: inp.fsE, sense: 'min', explain: 'Sliding with pseudo-static earth pressure, wall inertia and hydrodynamic suction.' },
          { key: 'Mmax', label: 'Overturning moment of the earth thrust', value: r(pS.Peff * pS.armAboveBase / 1e3, 6), unit: 'kN·m/m', explain: 'Moment of the horizontal earth and surcharge thrust about the base: the main overturning action on the wall.' },
        ] };
    }
    const checks = out.checks, w = checks.reduce((m, c) => (c.uc > m.uc ? c : m)), ok = w.uc <= 1;
    const zs = linspace(0, H, 29), pAct = (K, q, z) => K * (q + (z <= zw ? g * z : g * zw + gs * (z - zw)));
    const steelPerM = sheet ? inp.wallMass / 1e3 * (H + inp.dProv) : 0; // t per metre run of quay
    const concPerM = sheet ? 0 : inp.B * H;
    return {
      kpis: [
        { key: 'Ka', label: 'Active earth pressure coefficient (horizontal)', value: r(cS.Ka), unit: '-', explain: `Coulomb coefficient with wall friction, horizontal component; the Rankine value without wall friction is ${f(rk.Ka)}.` },
        { key: 'Kp', label: 'Passive earth pressure coefficient (horizontal)', value: r(cS.Kp), unit: '-', explain: 'Coulomb coefficient with the passive-side wall friction entered. It is divided by its factor of safety in the embedment calculation.' },
        { key: 'KAE', label: 'Seismic active coefficient (Mononobe–Okabe)', value: r(cE.KaTot), unit: '-', explain: 'Total active coefficient with the horizontal seismic coefficient; the rise over the static value is the dynamic earth pressure increment.' },
        { key: 'Pa', label: 'Active thrust on the retained height', value: r(pS.Peff / 1e3), unit: 'kN/m', explain: `Horizontal earth and surcharge thrust over the retained height, acting ${f(pS.armAboveBase)} m above dredged level.` },
        { key: 'PaE', label: 'Active thrust in the earthquake case', value: r(pE.Peff / 1e3), unit: 'kN/m', explain: 'Thrust with the Mononobe–Okabe coefficient and half the cargo surcharge.' },
        ...out.kp,
        { key: 'qPerm', label: 'Maximum permissible cargo surcharge', value: r(out.qPerm / 1e3), unit: 'kPa', explain: 'Largest uniform load behind the cope for which every static check stays satisfied (reported at 400 kPa when not reached). It is the apron and yard load rating set by the wall.' },
      ],
      checks,
      plots: [
        { type: 'line', title: 'Lateral pressure on the wall', x: { label: 'Pressure', unit: 'kPa' }, y: { label: 'Level below the cope', unit: 'm' },
          series: [{ name: 'Active earth pressure, static', x: zs.map((z) => r(pAct(cS.Ka, inp.q * 1e3, z) / 1e3)), y: zs.map((z) => r(-z)) }, ...(inp.kh > 0 ? [{ name: 'Active earth pressure, earthquake', x: zs.map((z) => r(pAct(cE.Ka, 0.5 * inp.q * 1e3, z) / 1e3)), y: zs.map((z) => r(-z)), dash: true }] : []), { name: 'Net water pressure (tidal lag)', x: zs.map((z) => r(GW * Math.min(Math.max(z - (zw - lag), 0), lag) / 1e3)), y: zs.map((z) => r(-z)) }],
          explain: 'Horizontal pressure against depth on the retained side. The kink is at the water table, where the soil becomes buoyant; the surcharge adds a constant strip and the tidal lag a net water pressure.' },
        { type: 'bar', title: 'Utilisation by check', categories: checks.map((c) => c.name), series: [{ name: 'Utilisation', y: checks.map((c) => r(c.uc, 3)) }], x: { label: 'Check', unit: '' }, y: { label: 'Utilisation', unit: '-' },
          explain: 'Demand over capacity (or required over achieved factor of safety) for each stability and strength check. The tallest bar sets the wall dimensions.' },
        { type: 'line', title: 'Governing utilisation against cargo surcharge', x: { label: 'Surcharge', unit: 'kPa' }, y: { label: 'Utilisation', unit: '-' },
          series: (() => { const qs = linspace(0, Math.max(2 * inp.q, 80), 9); const ucs = qs.map((qq) => { try { if (sheet) { const s = sheetPileWall({ H, gamma: g, gammaSub: gs, waterDepth: zw, Ka: cS.Ka, Kp: cS.Kp, fsPassive: inp.fsS, surcharge: qq * 1e3, anchorDepth: inp.anchorZ, roweFactor: inp.rowe, n: 400 }); return Math.max(s.d / inp.dProv, inp.gF * s.Mmax / out.MRd); } const p = thrust(cS.Ka, qq * 1e3); const fr = (out.W + p.Peff * Math.tan(inp.deltaA * DEG)) * Math.tan(inp.deltaB * DEG); return inp.fsS * (p.Peff + GW * lag * (lag / 2 + hw) + inp.bollard * 1e3) / fr; } catch { return 99; } }); return [{ name: sheet ? 'Embedment or bending' : 'Sliding', x: qs.map((v) => r(v)), y: ucs.map((v) => r(Math.min(v, 5), 3)) }, { name: 'Limit', x: [0, r(qs[qs.length - 1])], y: [1, 1], dash: true }]; })(),
          explain: 'How the governing check responds to the load placed behind the cope. Where the curve crosses 1 is the surcharge rating of the quay apron.' },
      ],
      tables: [{ title: 'Design cases', columns: [{ key: 'item', label: 'Quantity' }, { key: 'st', label: 'Service' }, { key: 'eq', label: 'Earthquake' }, { key: 'unit', label: 'Unit' }],
        rows: sheet ? [{ item: 'Active coefficient (horizontal)', st: r(cS.Ka), eq: r(cE.Ka), unit: '-' }, { item: 'Passive coefficient (horizontal)', st: r(cS.Kp), eq: r(cE.Kp), unit: '-' }, { item: 'Embedment required', st: r(out.st.d), eq: r(out.eq.d), unit: 'm' }, { item: 'Anchor force', st: r(out.T / 1e3), eq: r(out.TE / 1e3), unit: 'kN/m' }, { item: 'Largest moment', st: r(out.Mst / 1e3), eq: r(out.ME / 1e3), unit: 'kN·m/m' }]
          : [{ item: 'Active coefficient (horizontal)', st: r(cS.Ka), eq: r(cE.Ka), unit: '-' }, { item: 'Horizontal force', st: r(out.sS.Htot / 1e3), eq: r(out.sE.Htot / 1e3), unit: 'kN/m' }, { item: 'Vertical force', st: r(out.sS.Vtot / 1e3), eq: r(out.sE.Vtot / 1e3), unit: 'kN/m' }, { item: 'Sliding factor', st: r(cap(out.sS.fsSliding, 999)), eq: r(cap(out.sE.fsSliding, 999)), unit: '-' }, { item: 'Overturning factor', st: r(cap(out.sS.fsOverturning, 999)), eq: r(cap(out.sE.fsOverturning, 999)), unit: '-' }, { item: 'Toe pressure', st: r(out.sS.qMax / 1e3), eq: r(out.sE.qMax / 1e3), unit: 'kPa' }],
        explain: 'The service and earthquake cases side by side, with the coefficients used. The earthquake column uses half the cargo surcharge and the reduced factors of safety.' }],
      notes: [`${w.name} governs at ${f(w.uc)}. Water ${f(zw)} m below the cope with ${f(lag)} m of tidal lag; wall friction ${f(inp.deltaA)}° active and ${f(inp.deltaP)}° passive.`,
        sheet ? `Free-earth support: equilibrium embedment ${f(out.st.d0)} m, zero shear at ${f(out.st.zMmax)} m below the cope; moment resistance of the section ${f(out.MRd / 1e3)} kN·m/m.` : `Wall weight ${f(out.W / 1e3)} kN/m (buoyant below water); the width that would just meet the governing factors is about ${f(out.Bneed)} m.`],
      advice: [
        ok ? eng('ok', sheet ? `Embedment ${f(out.dReq)} m needed against ${f(inp.dProv)} m provided; anchor force ${f(out.T / 1e3)} kN/m; wall moment ${f(out.Mst / 1e3)} kN·m/m needs ${f(out.WelNeed)} cm³/m against ${f(inp.Wel)}.` : `Sliding ${f(cap(out.sS.fsSliding, 999))}, overturning ${f(cap(out.sS.fsOverturning, 999))}, bearing ${f(cap(out.sS.fsBearing, 999))}; toe pressure ${f(out.sS.qMax / 1e3)} kPa.`,
          `The apron can be rated at ${f(out.qPerm / 1e3)} kPa. Keep the drainage behind the wall working: each extra 0.5 m of tidal lag adds about ${f(GW * 0.5 * (hw + 0.25) / 1e3)} kN/m of thrust. Check overall stability on a deep slip surface separately.`, 'EN 1997-1:2004 §9; BS 6349-2:2019; EAU 2012',
          sheet ? `The wall uses about ${f(steelPerM)} t of steel per metre of quay (${steelCo2(steelPerM)}); a section of ${f(out.WelNeed)} cm³/m would suffice in bending.` : `The wall uses ${f(concPerM)} m³ of concrete per metre of quay (${concCo2(concPerM)}); caissons filled with dredged sand carry the same weight with a fraction of the cement.`)
          : eng('fail', `${w.name}: ${f(w.demand)} against ${f(w.capacity)} ${w.unit} (${f(w.uc)}).`,
            sheet ? (w.name.startsWith('Embedment') ? `Drive the wall to ${f(out.dReq)} m below dredged level, or limit the surcharge to ${f(out.qPerm / 1e3)} kPa.` : w.name.startsWith('Bending') ? `Use a section of at least ${f(out.WelNeed)} cm³/m, lower the anchor to reduce the span, or relieve the wall with a platform on piles.` : w.name.startsWith('Tie') ? `Close the tie spacing to about ${f(inp.tieSpacing / w.uc)} m or use larger rods.` : `Deepen the anchor wall (resistance grows with depth squared) or replace it with raking anchor piles.`)
              : `Widen the base to about ${f(out.Bneed)} m, improve the base friction with a rubble bed and shear key, or reduce the thrust with a lightweight or free-draining backfill; otherwise rate the apron at ${f(out.qPerm / 1e3)} kPa.`,
            'EN 1997-1:2004 §9.7; EN 1998-5:2004 §7; EAU 2012; PIANC WG34 (2001)',
            sheet ? `Each extra metre of wall length adds about ${f(steelPerM / (H + inp.dProv))} t of steel per metre of quay (${steelCo2(steelPerM / (H + inp.dProv))}); a relieving platform or better backfill reduces the thrust instead.` : `Widening to ${f(out.Bneed)} m adds ${f(Math.max(out.Bneed - inp.B, 0) * H)} m³ per metre of quay: ${concCo2(Math.max(out.Bneed - inp.B, 0) * H)}. Rock or slag backfill with a higher friction angle cuts the thrust without more concrete.`),
        eco(ok ? 'info' : 'warn', `The wall supports a cargo surcharge of ${f(out.qPerm / 1e3)} kPa against ${f(inp.q)} kPa applied; material is about ${sheet ? `${f(steelPerM)} t of steel` : `${f(concPerM)} m³ of concrete`} per metre of quay.`,
          'The surcharge rating is the commercial capacity of the apron. Where it is short, compare restricting the stacking zone to a setback behind the active wedge (no construction) with relieving the wall; dredging deeper for larger ships has to be paid for in embedment or base width first.', `${ECON_BASIS}; berth upgrade studies (PIANC WG 164)`),
      ],
      publish: { 'wharf.quayAnchorForce': sheet ? out.T : 0, 'wharf.quaySurchargePermissible': out.qPerm, 'wharf.quayThrust': pS.Peff },
      quantities: sheet ? { steel_t: r(steelPerM) } : { concrete_m3: r(concPerM) },
    };
  },
  verify: [
    // Rankine: dry sand φ′ = 30°, no wall friction, H = 10 m, γ = 18 kN/m³, no surcharge: Ka = 1/3, Pa = ½·Ka·γ·H² = 300 kN/m.
    // Gravity wall B = 6 m, γ = 24 kN/m³ in the dry: W = 1440 kN/m; base friction 30°: sliding FS = 1440·tan 30°/300 = 2.771;
    // overturning FS = (1440·3)/(300·10/3) = 4.32.
    { name: 'Rankine thrust and gravity-wall statics: Pa = ½KaγH², FS = W tan δ / Pa', source: 'Rankine (1857); Terzaghi, Peck & Mesri, Soil Mechanics in Engineering Practice §45–46 (hand calculation)',
      inputs: { wallType: 'gravity', H: 10, zw: 10, lag: 0, gamma: 18, phi: 30, deltaA: 0, deltaP: 0, q: 0, bollard: 0, kh: 0, B: 6, gammaWall: 24, deltaB: 30 },
      expect: [{ kpi: 'Ka', value: 0.3333, tol: 1e-3 }, { kpi: 'Pa', value: 300, tol: 1e-3 }, { kpi: 'fsSlide', value: 2.771, tol: 2e-3 }, { kpi: 'fsOT', value: 4.32, tol: 2e-3 }] },
    // Mononobe–Okabe, φ′ = 30°, δ = 0, kh = 0.1: ψ = atan 0.1 = 5.711°;
    // KAE = cos²(24.289°)/{cos²ψ·[1 + √(sin 30°·sin 24.289°/cos ψ)]²} = 0.8308/(0.9901·2.1160) = 0.3966.
    { name: 'Mononobe–Okabe coefficient for φ′ = 30°, kh = 0.1', source: 'Okabe (1926), Mononobe & Matsuo (1929); Kramer, Geotechnical Earthquake Engineering §11.6 (hand calculation)',
      inputs: { wallType: 'gravity', H: 10, zw: 10, lag: 0, phi: 30, deltaA: 0, deltaP: 0, q: 0, bollard: 0, kh: 0.1, B: 8 }, expect: [{ kpi: 'KAE', value: 0.3966, tol: 2e-3 }] },
  ],
};
