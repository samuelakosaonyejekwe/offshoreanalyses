// On-bottom stability: absolute lateral static stability and vertical stability in the DNV-RP-F109 format, with the
// submerged weight and concrete thickness required and the effect of penetration and trenching.
import { nearBedSpectralVelocity, designOscillation, currentOverPipe, initialPenetration, verticalStability, morisonCycle,
  absoluteLateralStability, requiredSubmergedWeight, requiredConcreteThickness } from '../../engine/pipeline/stability.js';
import { pipeInputs, pipeFrom, lengthInput, depthInput, num, sel, kpi, chk, adv, plot, tbl, col, run, need, fmt, round, linspace, sev, arr,
  ECON_BASIS, CO2_CONCRETE, concreteCo2, steelCo2 } from './_common.js';

const Z0 = { clay: 5e-6, 'fine-sand': 1e-5, 'medium-sand': 4e-5, 'coarse-sand': 1e-4, gravel: 3e-4, pebble: 2e-3, cobble: 1e-2 };
const G_ENV = 'Design storm';
const G_SOIL = 'Seabed';
const G_HYD = 'Hydrodynamic load model';

export default {
  id: 's08-on-bottom-stability',
  suite: 8,
  title: 'On-bottom stability',
  summary: 'Checks whether the pipe stays put on the seabed in the design storm and sizes the weight coating if it does not.',
  covers: [1],
  fidelity: 'analytical',
  standards: ['DNV-RP-F109 (2021) §3.4–3.6, Eq. 3.1, 3.38, 3.39', 'DNV RP E305 (1988) §5.3 (Morison route)'],
  models: ['JONSWAP spectrum transferred to the seabed (linear wave theory)', 'Morison equation with worst-phase search', 'Coulomb friction', 'Passive soil resistance (Verley & Sotberg / Verley & Lund)', 'Logarithmic current boundary layer'],
  limits: 'Absolute (zero-displacement) static criterion only. Peak loads come from a Morison worst-phase search with effective coefficients; the DNV-RP-F109 peak-load coefficient tables (C_Y*, C_Z*) are not built in — select the coefficient method and type them in to apply the RP tables. The safety factor γ_SC has no default in the engine and must be set from the F109 table for the region, soil and safety class. The generalised (allowable-displacement) method and 3-D CFD are outside this module.',
  publishes: ['pipeline.submergedWeight', 'quantities.concrete_m3'],
  searchTerms: 'subsea pipeline on-bottom stability DNV-RP-F109 concrete weight coating',
  inputs: [
    ...pipeInputs(),
    lengthInput(10, 'Pipe'),
    sel('condition', 'Pipe condition', 'operating', [['operating', 'Operation (product filled)'], ['empty', 'Installation (empty)'], ['flooded', 'Water filled']], 'Pipe', 'ic',
      'Contents state assumed during the storm. The empty installation case usually governs the weight coating.'),
    depthInput(80, G_ENV),
    num('Hs', 'Significant wave height', 'm', 8, 0, 20, G_ENV, 'bc', 'Design-storm sea state (100-year wave with 10-year current, or the reverse, in operation). Zero gives the current-only case.', { bind: 'site.Hs100', uncertain: { dist: 'gumbel', cov: 0.1 } }),
    num('Tp', 'Peak period', 's', 13, 3, 25, G_ENV, 'bc', 'Spectral peak period of the design sea state.', { bind: 'site.Tp100' }),
    num('RD', 'Directionality and spreading factor', '-', 1, 0.3, 1, G_ENV, 'bc', 'Reduction of the near-bed velocity normal to the pipe for wave heading and short-crestedness (F109 §3.4.4). 1.0 = long-crested waves normal to the pipe.'),
    num('duration', 'Storm duration', 'h', 3, 0.5, 24, G_ENV, 'bc', 'Duration over which the single design oscillation is extracted.'),
    num('Vref', 'Current speed at reference height', 'm/s', 0.5, 0, 3, G_ENV, 'bc', 'Steady current used with the design waves.', { bind: 'site.UcBed100', uncertain: { dist: 'gumbel', cov: 0.15 } }),
    sel('currentRef', 'Current definition', 'log', [['log', 'At a reference height (log boundary layer)'], ['pipe', 'Already the mean over the pipe height']], G_ENV, 'bc', 'How the current speed above is converted to the mean velocity over the pipe.'),
    num('zRef', 'Reference height above seabed', 'm', 1, 0.2, 20, G_ENV, 'bc', 'Height at which the current speed applies.', { showIf: { key: 'currentRef', in: ['log'] } }),
    sel('bed', 'Seabed roughness class', 'medium-sand', Object.keys(Z0).map((k) => [k, k.replace('-', ' ')]), G_ENV, 'soil', 'Sets the bottom roughness z0 of the current boundary layer (F109 Table 3-1).', { showIf: { key: 'currentRef', in: ['log'] } }),
    num('currentAngle', 'Angle between current and pipe axis', 'deg', 90, 0, 90, G_ENV, 'bc', '90° = current normal to the pipe.'),
    sel('soil', 'Soil model for passive resistance', 'sand', [['sand', 'Sand'], ['clay', 'Clay'], ['none', 'None (rock or friction only)']], G_SOIL, 'soil', 'Passive resistance from penetration per F109 §3.4.6.'),
    num('mu', 'Pipe–soil friction coefficient', '-', 0.6, 0.05, 1.5, G_SOIL, 'soil', 'F109: 0.6 on sand and rock, 0.2 on clay (concrete-coated pipe).', { uncertain: { dist: 'lognormal', cov: 0.15 } }),
    num('gammaSoil', 'Submerged unit weight of sand', 'N/m³', 9500, 4000, 13000, G_SOIL, 'soil', 'Used for sand passive resistance and initial penetration.', { bind: 'site.gammaSub', showIf: { key: 'soil', in: ['sand'] } }),
    num('su', 'Undrained shear strength of clay', 'kPa', 10, 0.5, 500, G_SOIL, 'soil', 'Strength at pipe level.', { bind: 'site.su0', showIf: { key: 'soil', in: ['clay'] } }),
    num('gammaDry', 'Dry unit weight of clay', 'N/m³', 18000, 10000, 22000, G_SOIL, 'soil', 'Enters the clay penetration parameter G_c.', { showIf: { key: 'soil', in: ['clay'] } }),
    sel('zpMode', 'Penetration', 'initial', [['initial', 'Initial penetration from the F109 expressions'], ['given', 'Specified'], ['zero', 'None']], G_SOIL, 'ic', 'Penetration reduces the hydrodynamic load and adds passive resistance.'),
    num('zpGiven', 'Specified penetration', 'mm', 30, 0, 800, G_SOIL, 'ic', 'Total penetration of the pipe invert below the seabed.', { showIf: { key: 'zpMode', in: ['given'] } }),
    num('trenchAngle', 'Trench side-slope angle (0 = no trench)', 'deg', 0, 0, 45, G_SOIL, 'geometry', 'Trench shielding per F109 Eq. 3.22–3.23 (recalled expressions — confirm). Applied for 5–45°.'),
    sel('method', 'Peak load method', 'morison', [['morison', 'Morison, worst phase (E305-type coefficients)'], ['coefficients', 'Peak coefficients C_Y*, C_Z* (from the F109 tables)']], G_HYD, 'control', 'The F109 peak-coefficient tables are not built in; read C_Y*, C_Z* for the K* and M* reported below and type them in.'),
    num('CY', 'Peak horizontal coefficient C_Y*', '-', 1.2, 0, 15, G_HYD, 'load', 'From DNV-RP-F109 Table 3-9 at K* and M*.', { showIf: { key: 'method', in: ['coefficients'] } }),
    num('CZ', 'Peak vertical coefficient C_Z*', '-', 1.0, 0, 15, G_HYD, 'load', 'From DNV-RP-F109 Table 3-10 at K* and M*.', { showIf: { key: 'method', in: ['coefficients'] } }),
    num('CD', 'Drag coefficient', '-', 0.7, 0.3, 2, G_HYD, 'load', 'Effective value for a pipe on the bed (E305: 0.7).', { showIf: { key: 'method', in: ['morison'] } }),
    num('CM', 'Inertia coefficient', '-', 3.29, 1, 5, G_HYD, 'load', 'Effective value for a pipe on the bed (E305: 3.29).', { showIf: { key: 'method', in: ['morison'] } }),
    num('CL', 'Lift coefficient', '-', 0.9, 0, 2, G_HYD, 'load', 'Effective value for a pipe on the bed (E305: 0.9).', { showIf: { key: 'method', in: ['morison'] } }),
    num('Fw', 'Calibration factor on Morison loads', '-', 1, 0.5, 2, G_HYD, 'load', 'E305 Fig. 5.12 factor F_w (function of K and M); 1.0 leaves the loads uncalibrated.', { showIf: { key: 'method', in: ['morison'] } }),
    num('gammaSC', 'Safety factor γ_SC', '-', 1.1, 1, 2, 'Acceptance criteria', 'criteria', 'No default exists: take it from DNV-RP-F109 Tables 3-5 to 3-8 for the region, soil type and safety class (roughly 1.0–1.5). The value shown is a placeholder.'),
    num('gammaW', 'Vertical stability factor γ_w', '-', 1.1, 1, 1.5, 'Acceptance criteria', 'criteria', 'F109 Eq. 3.1: specific gravity of the pipe must be at least γ_w.'),
  ],
  compute(inp, ctx) {
    const p = pipeFrom(inp, ctx);
    const rhoByCond = { operating: inp.rhoCont, empty: 0, flooded: p.rhoW };
    const ws = inp.condition === 'empty' ? p.wsEmpty : inp.condition === 'flooded' ? p.wsFlood : p.ws;
    need(inp.depth > 1, 'Water depth must be more than 1 m for the wave transfer to the seabed.');
    need(ws > 0, `The pipe is buoyant in the ${inp.condition} condition (submerged weight ${fmt(ws)} N/m): add weight coating or wall thickness before checking lateral stability.`);
    const D = p.Do;
    // design oscillation at the seabed
    let Ustar = 0, Tstar = inp.Tp, Us = 0, Tu = inp.Tp;
    if (inp.Hs > 0) {
      const nb = run(() => nearBedSpectralVelocity({ Hs: inp.Hs, Tp: inp.Tp, depth: inp.depth }), 'Near-bed wave velocity');
      const os = designOscillation({ Us: nb.Us, Tu: nb.Tu, durationS: inp.duration * 3600, RD: inp.RD, gamma: nb.gamma, Tn: nb.Tn });
      Ustar = os.Ustar; Tstar = os.Tstar; Us = os.UsReduced; Tu = nb.Tu;
    }
    const sinA = Math.sin(inp.currentAngle * Math.PI / 180);
    const Vstar = inp.currentRef === 'pipe' ? inp.Vref * sinA
      : (inp.Vref > 0 ? currentOverPipe({ Vref: inp.Vref, zRef: inp.zRef, z0: Z0[inp.bed], D, angleDeg: inp.currentAngle }) : 0);
    need(Ustar + Vstar > 0, 'Both the wave-induced velocity and the current are zero: there is no hydrodynamic load to check.');
    const soilArgs = inp.soil === 'sand' ? { soil: 'sand', gammaSoilSub: inp.gammaSoil } : inp.soil === 'clay' ? { soil: 'clay', su: inp.su * 1e3, gammaSoilDry: inp.gammaDry } : { soil: 'none' };
    let zp = 0;
    if (inp.zpMode === 'given') zp = inp.zpGiven / 1e3;
    else if (inp.zpMode === 'initial' && inp.soil !== 'none') zp = Math.min(run(() => initialPenetration({ ...soilArgs, FC: ws, D }), 'Initial penetration').zpi, 0.5 * D);
    const base = { D, mu: inp.mu, Ustar, Vstar, Tstar, gammaSC: inp.gammaSC, method: inp.method, CY: inp.CY, CZ: inp.CZ, Fw: inp.Fw, CD: inp.CD, CM: inp.CM, CL: inp.CL, rhoW: p.rhoW, zp, ...soilArgs };
    if (inp.trenchAngle >= 5) base.trenchAngleDeg = inp.trenchAngle;
    const res = run(() => absoluteLateralStability({ ...base, ws }), 'Lateral stability');
    const vert = verticalStability({ ws, b: p.buoy, gammaW: inp.gammaW });
    const req = run(() => requiredSubmergedWeight({ ...base }), 'Required submerged weight');
    // concrete thickness that satisfies both criteria for this condition (penetration held at the present value)
    const { D: _d, ...stab } = base;
    let conc = null;
    try {
      conc = requiredConcreteThickness({ pipe: { D: p.D, t: p.t, rhoSteel: p.rhoS, coatings: p.coatings.filter((c) => c.name !== 'concrete'), rhoContents: rhoByCond[inp.condition], rhoWater: p.rhoW },
        concreteDensity: p.concRho, concreteAbsorption: p.concAbs, stability: stab, gammaW: inp.gammaW, tMax: 0.25 });
    } catch { conc = null; }
    const cyc = morisonCycle({ D, Ustar, Vstar, Tstar, rhoW: p.rhoW, CD: inp.CD, CM: inp.CM, CL: inp.CL, nPhase: 72 });
    const r = res.reduction;
    const FYc = Array.from(cyc.FY, (v) => v * inp.Fw * r.rTotY), FZc = Array.from(cyc.FZ, (v) => v * inp.Fw * r.rTotZ);
    // sweep of the sea state
    const hs = linspace(0.25 * Math.max(inp.Hs, 2), 1.5 * Math.max(inp.Hs, 2), 11);
    const ucHs = hs.map((H) => {
      const nb = nearBedSpectralVelocity({ Hs: H, Tp: inp.Tp, depth: inp.depth, nOmega: 300 });
      const os = designOscillation({ Us: nb.Us, Tu: nb.Tu, durationS: inp.duration * 3600, RD: inp.RD, gamma: nb.gamma, Tn: nb.Tn });
      return absoluteLateralStability({ ...base, ws, Ustar: os.Ustar, Tstar: os.Tstar }).uc;
    });
    const alt = (extra) => { try { return absoluteLateralStability({ ...base, ws, ...extra }).uc; } catch { return NaN; } };
    const L = inp.Lkm * 1e3;
    const tConcReq = conc ? conc.thickness : NaN;
    const dConc = Number.isFinite(tConcReq) ? tConcReq - inp.tConc / 1e3 : NaN;
    const rIn = p.D / 2 + inp.tCoat / 1e3;
    const volAt = (tc) => Math.PI * ((rIn + tc) ** 2 - rIn ** 2) * L;
    const dVol = Number.isFinite(tConcReq) ? volAt(tConcReq) - volAt(inp.tConc / 1e3) : 0;
    const checks = [
      chk('Lateral stability (absolute static)', inp.gammaSC * (res.FY + inp.mu * res.FZ), inp.mu * ws + res.FR, 'N/m', 'DNV-RP-F109 Eq. 3.38',
        'Guards against the pipe sliding sideways in the design oscillation: factored drag, inertia and lift-reduced friction against Coulomb friction plus passive soil resistance. Exceeding it means the pipe moves in the design storm.'),
      chk('Lift-off (vertical load vs weight)', inp.gammaSC * res.ucLift * ws / inp.gammaSC, ws, 'N/m', 'DNV-RP-F109 Eq. 3.39',
        'Guards against the hydrodynamic lift exceeding the submerged weight; beyond it the pipe loses contact with the seabed and all friction.'),
      chk('Vertical stability (flotation)', inp.gammaW * p.buoy, ws + p.buoy, 'N/m', 'DNV-RP-F109 Eq. 3.1',
        'Guards against the pipe floating: the weight in the assessed condition must exceed the buoyancy by the factor γ_w.'),
    ];
    const uc = Math.max(res.uc, vert.uc);
    const table = [
      { option: 'As specified', detail: `${fmt(inp.tConc)} mm concrete, penetration ${fmt(zp * 1e3)} mm`, uc: round(res.uc, 3) },
      { option: 'Concrete thickness for UC = 1', detail: Number.isFinite(tConcReq) ? `${fmt(tConcReq * 1e3)} mm (specific gravity ${fmt(conc.sg)})` : 'not reachable within 250 mm', uc: Number.isFinite(tConcReq) ? 1 : 99 },
      { option: 'Trench with 30° side slopes', detail: 'Load reduction r_tr applied', uc: round(alt({ trenchAngleDeg: 30 }), 3) },
      { option: 'Penetration 0.25 D (self-burial or pre-lay)', detail: `${fmt(0.25 * D * 1e3)} mm`, uc: round(inp.soil === 'none' ? alt({ zp: 0.25 * D, soil: 'none' }) : alt({ zp: 0.25 * D }), 3) },
    ].filter((x) => Number.isFinite(x.uc));
    return {
      kpis: [
        kpi('ucLat', 'Lateral stability utilisation', res.ucLateral, '-', 'Factored driving load over resistance in the design oscillation (F109 Eq. 3.38). Above 1.0 the pipe is predicted to move; the value includes γ_SC.', 1, 'max'),
        kpi('sfLat', 'Lateral stability factor', inp.gammaSC / Math.max(res.ucLateral, 1e-9), '-', 'Unfactored resistance over driving load. It must be at least the safety factor γ_SC that the project takes from the F109 table.', inp.gammaSC, 'min'),
        kpi('ucLift', 'Lift utilisation', res.ucLift, '-', 'Factored peak lift over submerged weight (F109 Eq. 3.39). Above 1.0 the pipe lifts off the seabed.', 1, 'max'),
        kpi('sg', 'Specific gravity', vert.sg, '-', 'Weight in the assessed condition over buoyancy. Below γ_w the pipe can float or be lifted during installation.', inp.gammaW, 'min'),
        kpi('ws', 'Submerged weight (assessed condition)', ws, 'N/m', 'Weight in water of the pipe in the selected condition, the resisting quantity of the check.'),
        kpi('wsReq', 'Submerged weight required', req.ws, 'N/m', 'Submerged weight at which the lateral criterion is just met for this diameter, soil and storm.'),
        kpi('tConcReq', 'Concrete thickness required', Number.isFinite(tConcReq) ? tConcReq * 1e3 : 250, 'mm', 'Weight-coating thickness that satisfies both the lateral and the vertical criterion, with the larger diameter it creates taken into account. 250 means not reachable with concrete alone.'),
        kpi('Ustar', 'Design oscillation velocity U*', Ustar, 'm/s', 'Amplitude of the largest near-bed wave-induced velocity expected in the storm, normal to the pipe.'),
        kpi('Tstar', 'Design oscillation period T*', Tstar, 's', 'Period of the design oscillation; with U* it sets the inertia load and the Keulegan–Carpenter number.'),
        kpi('Vstar', 'Current over the pipe V*', Vstar, 'm/s', 'Mean steady current over the pipe height, normal to the pipe.'),
        kpi('Kstar', 'Keulegan–Carpenter number K*', res.K, '-', 'U*·T*/D. Use it with M* to read C_Y* and C_Z* from the F109 tables.'),
        kpi('Mstar', 'Current-to-wave ratio M*', Number.isFinite(res.M) ? res.M : 99, '-', 'V*/U*. 99 is shown for the current-only case.'),
        kpi('FYpeak', 'Peak horizontal load', res.FY, 'N/m', 'Horizontal hydrodynamic load at the governing phase, after penetration and trench reduction.'),
        kpi('FZpeak', 'Lift at the governing phase', res.FZ, 'N/m', 'Vertical hydrodynamic load at the phase that governs sliding; it removes friction.'),
        kpi('zp', 'Penetration used', zp * 1e3, 'mm', 'Pipe penetration into the seabed used for load reduction and passive resistance.'),
        kpi('FR', 'Passive soil resistance', res.FR, 'N/m', 'Resistance of the soil ridge ahead of the penetrated pipe, added to Coulomb friction.'),
      ],
      checks,
      plots: [
        plot('line', 'Hydrodynamic load through the design oscillation', 'Phase', 'deg', 'Load per unit length', 'N/m', [
          { name: 'Horizontal F_Y', x: arr(cyc.phaseDeg), y: arr(FYc) }, { name: 'Lift F_Z', x: arr(cyc.phaseDeg), y: arr(FZc) },
          { name: 'Sliding resistance μ(w_s − F_Z) + F_R', x: arr(cyc.phaseDeg), y: FZc.map((fz) => round(Math.max(0, inp.mu * (ws - fz)) + res.FR)) , dash: true },
        ], 'Morison loads over one design oscillation with the unfactored sliding resistance at each phase. The pipe is closest to sliding where the horizontal load approaches the dashed line.'),
        plot('line', 'Utilisation vs significant wave height', 'Hs', 'm', 'Utilisation', '-', [
          { name: 'Governing utilisation', x: arr(hs), y: arr(ucHs) }, { name: 'Limit', x: [round(hs[0]), round(hs[hs.length - 1])], y: [1, 1], dash: true },
        ], 'How the stability utilisation grows with the storm. The crossing of the limit line is the largest sea state the pipe survives without moving.'),
      ],
      tables: [tbl('Stabilisation options', [col('option', 'Option'), col('detail', 'Detail'), col('uc', 'Utilisation', '-')], table,
        'Each row is one way of achieving stability with its governing utilisation, so that weight coating, trenching and reliance on embedment can be compared.')],
      notes: [
        `Method: ${res.method}. Governing phase ${res.phaseDeg ?? 0}°; K* = ${fmt(res.K)}, M* = ${fmt(res.M)}.`,
        `Significant near-bed velocity ${fmt(Us)} m/s at ${fmt(Tu)} s; design oscillation ${fmt(Ustar)} m/s.`,
        'γ_SC is a project input from the DNV-RP-F109 tables; the peak-load coefficient tables of the RP are not built in.',
      ],
      advice: [
        adv('engineering', sev(uc), `Governing stability utilisation is ${fmt(uc)} in the ${inp.condition} condition (lateral ${fmt(res.ucLateral)}, lift ${fmt(res.ucLift)}, flotation ${fmt(vert.uc)}).`,
          uc > 1 ? (Number.isFinite(tConcReq) ? `Increase the concrete coating from ${fmt(inp.tConc)} to at least ${fmt(tConcReq * 1e3)} mm, or trench the line (see the options table).` : 'Concrete alone cannot stabilise this pipe: trench, rock-cover or anchor the section.')
            : `Stable with ${fmt((1 - uc) * 100)} % reserve. Confirm γ_SC from the F109 table for the region and repeat for the empty installation case and the 10-year/100-year combinations.`,
          'DNV-RP-F109 (2021) §3.6 absolute lateral static stability; γ_SC from Tables 3-5 to 3-8 (project input)'),
        adv('economic', 'info', Number.isFinite(dConc) ? `The required coating differs from the specified one by ${fmt(dConc * 1e3)} mm, i.e. ${fmt(dVol)} m³ of concrete over ${fmt(inp.Lkm)} km.` : 'No concrete thickness up to 250 mm stabilises the line; the cost driver becomes trenching or rock placement.',
          dConc < -0.005 ? `About ${fmt(-dVol)} m³ of weight coating could be removed, which also lowers lay tension and vessel requirements.` : dConc > 0 ? `Budget ${fmt(dVol)} m³ of additional coating, or compare with the vessel days of a trenching spread.` : 'The coating is matched to the storm; no saving is available without changing the route or burial.',
          ECON_BASIS,
          `${concreteCo2(dVol)} ${steelCo2(p.steelKgM * L / 1e3)}`),
      ],
      publish: { 'pipeline.submergedWeight': p.ws, 'quantities.concrete_m3': p.concM3M * L },
      quantities: { concrete_m3: round(p.concM3M * L), steel_t: round(p.steelKgM * L / 1e3), pipe_km: inp.Lkm },
    };
  },
  verify: [{
    // Current only (Hs = 0), bare 400 mm pipe, V = 1 m/s over the pipe, no penetration, friction only.
    // F_Y = ½·1025·0.4·0.7·1² = 143.5 N/m; F_Z = ½·1025·0.4·0.9·1² = 184.5 N/m.
    // UC = γ_SC (F_Y + μ F_Z)/(μ w_s) = 1.1·(143.5 + 0.6·184.5)/(0.6·1000) = 0.4660.
    name: 'Steady current on a bare pipe: Morison drag and lift with Coulomb friction',
    source: 'Hand calculation, Morison equation with C_D = 0.7, C_L = 0.9 (DNV RP E305 §5.3)',
    inputs: { D: 400, t: 20, tCoat: 0, tConc: 0, wsGiven: 1000, Hs: 0, Vref: 1, currentRef: 'pipe', soil: 'none', zpMode: 'zero', mu: 0.6, gammaSC: 1.1 },
    expect: [{ kpi: 'FYpeak', value: 143.5, tol: 0.005 }, { kpi: 'FZpeak', value: 184.5, tol: 0.005 }, { kpi: 'ucLat', value: 0.466, tol: 0.005 }],
  }],
};
