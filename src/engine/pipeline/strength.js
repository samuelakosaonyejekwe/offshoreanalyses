/**
 * Pipeline wall-thickness and strength checks: DNV-ST-F101 limit states (pressure containment, collapse,
 * propagating buckling, combined-loading local buckling), ASME B31.4 / B31.8 stress checks, thick-wall
 * stresses, bending strain, ovalisation and bend (elbow) factors.
 *
 * All quantities SI. Pressures are absolute or gauge consistently (only differences enter the checks).
 */

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function nonNeg(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}
function fin(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v))) throw new RangeError(`${name} must be a finite number`);
  return v;
}

const SC = { low: 0, medium: 1, high: 2 };
function scIndex(safetyClass) {
  const i = SC[String(safetyClass).toLowerCase()];
  if (i === undefined) throw new RangeError('safetyClass must be "low", "medium" or "high"');
  return i;
}

/**
 * Partial resistance factors of DNV-ST-F101.
 *
 * Reference: DNV-ST-F101 (2021) Table 5-1 (γm = 1.15 for SLS/ULS/ALS, 1.00 for FLS), Table 5-2
 * (γSC: pressure containment 1.046 / 1.138 / 1.308; other limit states 1.04 / 1.14 / 1.26), strain
 * resistance factor γε 2.0 / 2.5 / 3.3 (Table 5-8) and the test-pressure factors of Table 5-9
 * (α_mpt 1.000 / 1.088 / 1.251; α_spt 1.03 / 1.05 / 1.05) for safety class low / medium / high.
 *
 * @param {'low'|'medium'|'high'} safetyClass
 * @param {'SLS'|'ULS'|'ALS'|'FLS'} [limitState='ULS']
 * @returns {{gammaM:number, gammaScPc:number, gammaSc:number, gammaEps:number, alphaMpt:number, alphaSpt:number}} all [-]
 */
export function dnvResistanceFactors(safetyClass, limitState = 'ULS') {
  const i = scIndex(safetyClass);
  return {
    gammaM: String(limitState).toUpperCase() === 'FLS' ? 1.0 : 1.15,
    gammaScPc: [1.046, 1.138, 1.308][i],
    gammaSc: [1.04, 1.14, 1.26][i],
    gammaEps: [2.0, 2.5, 3.3][i],
    alphaMpt: [1.0, 1.088, 1.251][i],
    alphaSpt: [1.03, 1.05, 1.05][i],
  };
}

/**
 * Characteristic wall thicknesses t1 and t2.
 *
 * Reference: DNV-ST-F101 (2021) Table 5-6: t1 = t − t_fab (prior to operation), t − t_fab − t_corr (operation);
 * t2 = t (prior to operation), t − t_corr (operation). t1 is used where failure is governed by the weakest
 * section (pressure containment, collapse), t2 where it is governed by the average section.
 *
 * @param {object} o
 * @param {number} o.t nominal wall thickness [m]
 * @param {number} [o.tFab=0] fabrication (minus) tolerance as a thickness [m]
 * @param {number} [o.tFabFrac=0] fabrication tolerance as a fraction of t [-] (added to tFab)
 * @param {number} [o.tCorr=0] corrosion allowance [m]
 * @param {'operation'|'test'} [o.condition='operation'] 'test' = prior to operation (installation, system test)
 * @returns {{t1:number, t2:number, tFab:number}} [m]
 */
export function characteristicWallThickness({ t, tFab = 0, tFabFrac = 0, tCorr = 0, condition = 'operation' }) {
  pos('t', t); nonNeg('tFab', tFab); nonNeg('tFabFrac', tFabFrac); nonNeg('tCorr', tCorr);
  const fab = tFab + tFabFrac * t;
  const corr = condition === 'operation' ? tCorr : 0;
  const t1 = t - fab - corr;
  const t2 = t - corr;
  if (t1 <= 0) throw new RangeError('t is fully consumed by the fabrication tolerance and corrosion allowance');
  return { t1, t2, tFab: fab };
}

/**
 * Pressure containment (bursting) resistance p_b(t) = 2t/(D − t) · f_cb · 2/√3, f_cb = min(f_y, f_u/1.15).
 *
 * Reference: DNV-ST-F101 (2021) Eq. 5.8 and 5.9.
 * Validity: thin to moderately thick pipe (the mid-wall diameter form).
 *
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.t wall thickness to use (t1 or t2) [m]
 * @param {number} o.fy characteristic yield strength [Pa]
 * @param {number} o.fu characteristic tensile strength [Pa]
 * @returns {{pb:number, fcb:number}} [Pa]
 */
export function burstPressure({ D, t, fy, fu }) {
  pos('D', D); pos('t', t); pos('fy', fy); pos('fu', fu);
  if (t >= D / 2) throw new RangeError('t must be smaller than D/2');
  const fcb = Math.min(fy, fu / 1.15);
  return { pb: ((2 * t) / (D - t)) * fcb * (2 / Math.sqrt(3)), fcb };
}

/**
 * Local incidental pressure p_li = p_d·γ_inc + ρ_cont·g·(h_ref − h_l) and local test pressure likewise.
 *
 * Reference: DNV-ST-F101 (2021) Eq. 4.1 / 4.2. γ_inc is 1.10 for a typical pipeline control system and may
 * be reduced to 1.05 where the pressure safety system justifies it.
 *
 * @param {object} o
 * @param {number} o.pRef pressure at the reference elevation (design pressure, or system test pressure) [Pa]
 * @param {number} [o.gammaInc=1.1] incidental-to-design pressure ratio [-] (use 1 for a test pressure)
 * @param {number} [o.rhoCont=0] density of the contents (or test fluid) [kg/m³]
 * @param {number} [o.hRef=0] elevation of the reference point, positive upwards [m]
 * @param {number} [o.hLocal=0] elevation of the local point, positive upwards [m]
 * @param {number} [o.g=9.80665] [m/s²]
 * @returns {number} local incidental (or local test) pressure [Pa]
 */
export function localIncidentalPressure({ pRef, gammaInc = 1.1, rhoCont = 0, hRef = 0, hLocal = 0, g = 9.80665 }) {
  nonNeg('pRef', pRef); pos('gammaInc', gammaInc); nonNeg('rhoCont', rhoCont);
  return pRef * gammaInc + rhoCont * g * (hRef - hLocal);
}

/**
 * Mill pressure test p_mpt = k · 2 t_min/(D − t_min) · min(0.96 SMYS, 0.84 SMTS).
 *
 * Reference: DNV-ST-F101 (2021) Eq. 7.3. k = 1.0 for end-face sealing with no axial compression from the
 * test rig; lower values apply to other sealing arrangements (argument).
 *
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.tMin minimum wall thickness t − t_fab [m]
 * @param {number} o.smys [Pa]
 * @param {number} o.smts [Pa]
 * @param {number} [o.k=1] end-cap factor [-]
 * @returns {number} mill test pressure [Pa]
 */
export function millTestPressure({ D, tMin, smys, smts, k = 1 }) {
  pos('D', D); pos('tMin', tMin); pos('smys', smys); pos('smts', smts); pos('k', k);
  return k * ((2 * tMin) / (D - tMin)) * Math.min(0.96 * smys, 0.84 * smts);
}

/**
 * Pressure containment check in operation and in system test, with the mill-test compatibility checks.
 *
 *   operation:   p_li − p_e ≤ min( p_b(t1)/(γm γSC) ; p_lt/α_spt − p_e ; p_mpt α_U/α_mpt )
 *   system test: p_lt − p_e ≤ min( p_b(t1,test)/(γm γSC) ; p_mpt )
 *
 * Reference: DNV-ST-F101 (2021) Eq. 5.6 and 5.7 with Tables 5-1, 5-2, 5-9.
 * The second and third terms of the operation check are only evaluated when pLocalTest / pMill are given.
 *
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.t nominal wall thickness [m]
 * @param {number} [o.tFab=0] fabrication tolerance [m]
 * @param {number} [o.tFabFrac=0] fabrication tolerance as fraction of t [-]
 * @param {number} [o.tCorr=0] corrosion allowance [m]
 * @param {number} o.fy characteristic yield strength (de-rated, incl. α_U) [Pa]
 * @param {number} o.fu characteristic tensile strength (de-rated, incl. α_U) [Pa]
 * @param {number} o.pLocalInc local incidental pressure p_li [Pa]
 * @param {number} [o.pe=0] external pressure at the point [Pa]
 * @param {'low'|'medium'|'high'} [o.safetyClass='medium']
 * @param {number} [o.pLocalTest] local system test pressure p_lt [Pa]
 * @param {number} [o.pMill] mill test pressure p_mpt [Pa]
 * @param {number} [o.alphaU=0.96] material strength factor [-]
 * @param {number} [o.gammaM] override of γm [-]
 * @param {number} [o.gammaSc] override of γSC for pressure containment [-]
 * @param {number} [o.alphaSpt] override of the system test factor [-]
 * @param {number} [o.alphaMpt] override of the mill test factor [-]
 * @returns {{t1:number, t1Test:number, pb:number, pbDesign:number, demand:number, capacity:number, uc:number,
 *   governing:string, ucBurst:number, ucSystemTest:(number|null), ucMillTest:(number|null),
 *   testCheck:({demand:number, capacity:number, uc:number}|null)}}
 *   thicknesses [m]; pressures [Pa]; uc = demand/capacity [-]
 */
export function pressureContainment({ D, t, tFab = 0, tFabFrac = 0, tCorr = 0, fy, fu, pLocalInc, pe = 0,
  safetyClass = 'medium', pLocalTest, pMill, alphaU = 0.96, gammaM, gammaSc, alphaSpt, alphaMpt }) {
  const f = dnvResistanceFactors(safetyClass);
  const gm = gammaM === undefined ? f.gammaM : pos('gammaM', gammaM);
  const gsc = gammaSc === undefined ? f.gammaScPc : pos('gammaSc', gammaSc);
  const aSpt = alphaSpt === undefined ? f.alphaSpt : pos('alphaSpt', alphaSpt);
  const aMpt = alphaMpt === undefined ? f.alphaMpt : pos('alphaMpt', alphaMpt);
  fin('pLocalInc', pLocalInc); fin('pe', pe);
  const { t1 } = characteristicWallThickness({ t, tFab, tFabFrac, tCorr, condition: 'operation' });
  const t1Test = characteristicWallThickness({ t, tFab, tFabFrac, tCorr, condition: 'test' }).t1;
  const { pb } = burstPressure({ D, t: t1, fy, fu });
  const pbDesign = pb / (gm * gsc);
  const demand = pLocalInc - pe;
  const caps = [['burst', pbDesign]];
  if (pLocalTest !== undefined) caps.push(['system test', pos('pLocalTest', pLocalTest) / aSpt - pe]);
  if (pMill !== undefined) caps.push(['mill test', (pos('pMill', pMill) * alphaU) / aMpt]);
  let governing = caps[0];
  for (const c of caps) if (c[1] < governing[1]) governing = c;
  const ucOf = (name) => {
    const c = caps.find((x) => x[0] === name);
    return c ? demand / c[1] : null;
  };
  let testCheck = null;
  if (pLocalTest !== undefined) {
    const pbT = burstPressure({ D, t: t1Test, fy, fu }).pb / (gm * gsc);
    const capT = pMill !== undefined ? Math.min(pbT, pMill) : pbT;
    testCheck = { demand: pLocalTest - pe, capacity: capT, uc: (pLocalTest - pe) / capT };
  }
  return {
    t1, t1Test, pb, pbDesign, demand, capacity: governing[1], uc: demand / governing[1], governing: governing[0],
    ucBurst: demand / pbDesign, ucSystemTest: ucOf('system test'), ucMillTest: ucOf('mill test'), testCheck,
  };
}

/**
 * Required wall thickness for pressure containment (burst term only), inverted in closed form:
 * t1 = D·x/(1 + x), x = γm γSC (p_li − p_e)/(2 f_cb · 2/√3); nominal t = (t1 + t_corr + t_fab)/(1 − tFabFrac).
 *
 * Reference: DNV-ST-F101 (2021) Eq. 5.6, 5.8 rearranged.
 *
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.fy [Pa]
 * @param {number} o.fu [Pa]
 * @param {number} o.pLocalInc local incidental pressure [Pa]
 * @param {number} [o.pe=0] external pressure [Pa]
 * @param {'low'|'medium'|'high'} [o.safetyClass='medium']
 * @param {number} [o.tCorr=0] corrosion allowance [m]
 * @param {number} [o.tFab=0] fabrication tolerance [m]
 * @param {number} [o.tFabFrac=0] fabrication tolerance as fraction of the nominal thickness [-]
 * @param {number} [o.gammaM] override [-]
 * @param {number} [o.gammaSc] override [-]
 * @returns {{t1:number, tNominal:number}} required minimum t1 and nominal thickness [m]
 */
export function requiredWallThicknessBurst({ D, fy, fu, pLocalInc, pe = 0, safetyClass = 'medium', tCorr = 0,
  tFab = 0, tFabFrac = 0, gammaM, gammaSc }) {
  pos('D', D); pos('fy', fy); pos('fu', fu); nonNeg('tCorr', tCorr); nonNeg('tFab', tFab);
  if (!(tFabFrac >= 0 && tFabFrac < 1)) throw new RangeError('tFabFrac must be in [0, 1)');
  const f = dnvResistanceFactors(safetyClass);
  const gm = gammaM === undefined ? f.gammaM : gammaM;
  const gsc = gammaSc === undefined ? f.gammaScPc : gammaSc;
  const dp = pLocalInc - pe;
  if (!(dp > 0)) return { t1: 0, tNominal: (tCorr + tFab) / (1 - tFabFrac) };
  const fcb = Math.min(fy, fu / 1.15);
  const x = (gm * gsc * dp) / (2 * fcb * (2 / Math.sqrt(3)));
  const t1 = (D * x) / (1 + x);
  return { t1, tNominal: (t1 + tCorr + tFab) / (1 - tFabFrac) };
}

/**
 * Characteristic collapse pressure of a pipe under external pressure only:
 * (p_c − p_el)(p_c² − p_p²) = p_c p_el p_p f_0 D/t, with p_el = 2E(t/D)³/(1 − ν²), p_p = f_y α_fab 2t/D.
 * The physical root (0 < p_c < min(p_el, p_p)) is taken from the trigonometric solution of the cubic.
 *
 * Reference: DNV-ST-F101 (2021) Eq. 5.10–5.12; ovality f_0 = (D_max − D_min)/D is not to be taken below 0.005.
 * Validity: D/t up to about 45 (the calibration range of the code).
 *
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.t wall thickness (t1 for system collapse, t2 in combined-loading checks) [m]
 * @param {number} o.fy characteristic yield strength [Pa]
 * @param {number} [o.E=207e9] Young's modulus [Pa]
 * @param {number} [o.nu=0.3] Poisson's ratio [-]
 * @param {number} [o.alphaFab=1] fabrication factor [-]
 * @param {number} [o.f0=0.005] ovality [-]; values below 0.005 are raised to 0.005
 * @returns {{pc:number, pel:number, pp:number, f0:number, residual:number}}
 *   pressures [Pa]; residual of the cubic normalised by p_el p_p² [-]
 */
export function collapsePressure({ D, t, fy, E = 207e9, nu = 0.3, alphaFab = 1, f0 = 0.005 }) {
  pos('D', D); pos('t', t); pos('fy', fy); pos('E', E); pos('alphaFab', alphaFab); nonNeg('f0', f0);
  if (!(nu > 0 && nu < 0.5)) throw new RangeError('nu must be in (0, 0.5)');
  if (t >= D / 2) throw new RangeError('t must be smaller than D/2');
  const fo = Math.max(f0, 0.005);
  const pel = (2 * E * (t / D) ** 3) / (1 - nu * nu);
  const pp = fy * alphaFab * ((2 * t) / D);
  // p³ + b p² + c p + d = 0
  const b = -pel;
  const c = -(pp * pp + pp * pel * fo * (D / t));
  const d = pel * pp * pp;
  const u = (-b * b / 3 + c) / 3;
  const v = ((2 * b * b * b) / 27 - (b * c) / 3 + d) / 2;
  const arg = Math.max(-1, Math.min(1, -v / Math.sqrt(-u * u * u)));
  const phi = Math.acos(arg);
  const y = -2 * Math.sqrt(-u) * Math.cos(phi / 3 + Math.PI / 3);
  const pc = y - b / 3;
  const residual = ((pc - pel) * (pc * pc - pp * pp) - pc * pel * pp * fo * (D / t)) / (pel * pp * pp);
  return { pc, pel, pp, f0: fo, residual };
}

/**
 * System collapse check p_e − p_min ≤ p_c(t1)/(γm γSC).
 *
 * Reference: DNV-ST-F101 (2021) Eq. 5.10.
 *
 * @param {object} o
 * @param {number} o.D [m]
 * @param {number} o.t1 characteristic wall thickness t1 [m]
 * @param {number} o.fy [Pa]
 * @param {number} o.pe external pressure [Pa]
 * @param {number} [o.pMin=0] minimum internal pressure that can be sustained [Pa]
 * @param {'low'|'medium'|'high'} [o.safetyClass='medium']
 * @param {number} [o.E=207e9] [Pa]
 * @param {number} [o.nu=0.3] [-]
 * @param {number} [o.alphaFab=1] [-]
 * @param {number} [o.f0=0.005] [-]
 * @param {number} [o.gammaM] override [-]
 * @param {number} [o.gammaSc] override [-]
 * @returns {{pc:number, pel:number, pp:number, demand:number, capacity:number, uc:number}} pressures [Pa], uc [-]
 */
export function collapseCheck({ D, t1, fy, pe, pMin = 0, safetyClass = 'medium', E = 207e9, nu = 0.3, alphaFab = 1,
  f0 = 0.005, gammaM, gammaSc }) {
  const f = dnvResistanceFactors(safetyClass);
  const gm = gammaM === undefined ? f.gammaM : gammaM;
  const gsc = gammaSc === undefined ? f.gammaSc : gammaSc;
  const r = collapsePressure({ D, t: t1, fy, E, nu, alphaFab, f0 });
  const demand = fin('pe', pe) - pMin;
  const capacity = r.pc / (gm * gsc);
  return { pc: r.pc, pel: r.pel, pp: r.pp, demand, capacity, uc: demand / capacity };
}

/**
 * Required wall thickness for system collapse, by bisection on t (collapse pressure is monotonic in t).
 *
 * Reference: DNV-ST-F101 (2021) Eq. 5.10.
 *
 * @param {object} o same arguments as collapseCheck without t1
 * @returns {number} required t1 [m]
 */
export function requiredWallThicknessCollapse(o) {
  const demand = fin('pe', o.pe) - (o.pMin || 0);
  if (!(demand > 0)) return 0;
  let lo = o.D / 500;
  let hi = o.D * 0.49;
  if (collapseCheck({ ...o, t1: hi }).uc > 1) throw new RangeError('pe cannot be resisted by any wall thickness below D/2');
  for (let i = 0; i < 80; i++) {
    const mid = 0.5 * (lo + hi);
    if (collapseCheck({ ...o, t1: mid }).uc > 1) lo = mid; else hi = mid;
  }
  return hi;
}

/**
 * Propagating buckling: p_pr = 35 f_y α_fab (t2/D)^2.5 and the check p_e − p_min ≤ p_pr/(γm γSC).
 * When the check fails a buckle, once initiated, will run; buckle arrestors are then needed.
 *
 * Reference: DNV-ST-F101 (2021) Eq. 5.15, 5.16.
 * Validity: 15 < D/t2 < 45.
 *
 * @param {object} o
 * @param {number} o.D [m]
 * @param {number} o.t2 [m]
 * @param {number} o.fy [Pa]
 * @param {number} [o.alphaFab=1] [-]
 * @param {number} [o.pe=0] external pressure [Pa]
 * @param {number} [o.pMin=0] minimum internal pressure [Pa]
 * @param {'low'|'medium'|'high'} [o.safetyClass='medium']
 * @param {number} [o.gammaM] override [-]
 * @param {number} [o.gammaSc] override [-]
 * @returns {{ppr:number, demand:number, capacity:number, uc:number, arrestorsRequired:boolean, inValidityRange:boolean}}
 */
export function propagationPressure({ D, t2, fy, alphaFab = 1, pe = 0, pMin = 0, safetyClass = 'medium', gammaM, gammaSc }) {
  pos('D', D); pos('t2', t2); pos('fy', fy); pos('alphaFab', alphaFab);
  const f = dnvResistanceFactors(safetyClass);
  const gm = gammaM === undefined ? f.gammaM : gammaM;
  const gsc = gammaSc === undefined ? f.gammaSc : gammaSc;
  const ppr = 35 * fy * alphaFab * (t2 / D) ** 2.5;
  const capacity = ppr / (gm * gsc);
  const demand = pe - pMin;
  const r = D / t2;
  return { ppr, demand, capacity, uc: demand / capacity, arrestorsRequired: demand > capacity, inValidityRange: r > 15 && r < 45 };
}

/**
 * Integral buckle-arrestor crossover pressure
 * p_X = p_pr + (p_pr,BA − p_pr)·[1 − exp(−20 t2 L_BA / D²)] and the check p_e ≤ p_X/(1.1 γm γSC).
 *
 * Reference: DNV-ST-F101 (2021) Eq. 5.17, 5.18 (Torselletti et al.). The factor 1.1 on the resistance and
 * the coefficient 20 are as recalled from the standard — confirm against the edition in force.
 *
 * @param {object} o
 * @param {number} o.D pipe outer diameter [m]
 * @param {number} o.t2 pipe wall thickness [m]
 * @param {number} o.tBA arrestor wall thickness [m]
 * @param {number} o.LBA arrestor length [m]
 * @param {number} o.fy pipe yield strength [Pa]
 * @param {number} [o.fyBA=o.fy] arrestor yield strength [Pa]
 * @param {number} [o.alphaFab=1] pipe fabrication factor [-]
 * @param {number} [o.alphaFabBA=1] arrestor fabrication factor [-]
 * @param {number} [o.pe=0] external pressure [Pa]
 * @param {'low'|'medium'|'high'} [o.safetyClass='medium']
 * @param {number} [o.resistanceFactor=1.1] additional factor on the crossover resistance [-]
 * @returns {{ppr:number, pprBA:number, pX:number, capacity:number, uc:number}} pressures [Pa]
 */
export function buckleArrestorCapacity({ D, t2, tBA, LBA, fy, fyBA, alphaFab = 1, alphaFabBA = 1, pe = 0,
  safetyClass = 'medium', resistanceFactor = 1.1 }) {
  pos('tBA', tBA); pos('LBA', LBA);
  const f = dnvResistanceFactors(safetyClass);
  const ppr = 35 * pos('fy', fy) * alphaFab * (pos('t2', t2) / pos('D', D)) ** 2.5;
  const pprBA = 35 * (fyBA === undefined ? fy : fyBA) * alphaFabBA * (tBA / D) ** 2.5;
  const pX = ppr + (pprBA - ppr) * (1 - Math.exp((-20 * t2 * LBA) / (D * D)));
  const capacity = pX / (resistanceFactor * f.gammaM * f.gammaSc);
  return { ppr, pprBA, pX, capacity, uc: pe / capacity };
}

/**
 * Plastic capacities of a pipe: S_p = f_y π (D − t) t and M_p = f_y (D − t)² t.
 *
 * Reference: DNV-ST-F101 (2021) Eq. 5.20, 5.21.
 *
 * @param {object} o
 * @param {number} o.D [m]
 * @param {number} o.t [m]
 * @param {number} o.fy [Pa]
 * @returns {{Sp:number, Mp:number}} Sp [N], Mp [N·m]
 */
export function plasticCapacities({ D, t, fy }) {
  pos('D', D); pos('t', t); pos('fy', fy);
  return { Sp: fy * Math.PI * (D - t) * t, Mp: fy * (D - t) * (D - t) * t };
}

/**
 * Local buckling, load-controlled condition (combined bending, effective axial force and pressure).
 *
 * Internal overpressure (p_i > p_e):
 *   { γmγSC |M_Sd|/(α_c M_p) + [γmγSC S_Sd/(α_c S_p)]² }² + [α_p (p_i − p_e)/(α_c p_b(t2))]² ≤ 1
 * External overpressure (p_i ≤ p_e), with p_min the minimum internal pressure:
 *   { γmγSC |M_Sd|/(α_c M_p) + [γmγSC S_Sd/(α_c S_p)]² }² + [γmγSC (p_e − p_min)/p_c(t2)]² ≤ 1
 * α_c = (1 − β) + β f_u/f_y;  β = (60 − D/t2)/90;
 * α_p = 1 − β for (p_i − p_e)/p_b < 2/3, else 1 − 3β(1 − (p_i − p_e)/p_b).
 *
 * Reference: DNV-ST-F101 (2021) Eq. 5.19–5.24 and 5.28.
 * Validity: 15 ≤ D/t2 ≤ 45 and |S_Sd|/S_p < 0.4 (reported in `warnings`; the value is still returned).
 * M_Sd and S_Sd are design load effects (already including γF, γE, γc).
 *
 * @param {object} o
 * @param {number} o.D [m]
 * @param {number} o.t2 [m]
 * @param {number} o.fy [Pa]
 * @param {number} o.fu [Pa]
 * @param {number} o.MSd design bending moment [N·m]
 * @param {number} o.SSd design effective axial force, tension positive [N]
 * @param {number} [o.pi=0] internal pressure [Pa]
 * @param {number} [o.pe=0] external pressure [Pa]
 * @param {number} [o.pMin=o.pi] minimum internal pressure for the external-overpressure case [Pa]
 * @param {'low'|'medium'|'high'} [o.safetyClass='medium']
 * @param {number} [o.alphaFab=1] [-]
 * @param {number} [o.f0=0.005] ovality [-]
 * @param {number} [o.E=207e9] [Pa]
 * @param {number} [o.nu=0.3] [-]
 * @param {number} [o.gammaM] override [-]
 * @param {number} [o.gammaSc] override [-]
 * @returns {{mode:'internal'|'external', uc:number, momentTerm:number, axialTerm:number, pressureTerm:number,
 *   alphaC:number, alphaP:(number|null), beta:number, Sp:number, Mp:number, pb:number, pc:(number|null),
 *   allowableMoment:number, warnings:string[]}}
 *   uc is the left-hand side of the interaction equation (≤ 1 passes); allowableMoment [N·m] is the design
 *   moment that brings the interaction to 1 with the other load effects unchanged (0 if they alone exceed it)
 */
export function localBucklingLoadControlled({ D, t2, fy, fu, MSd, SSd, pi = 0, pe = 0, pMin, safetyClass = 'medium',
  alphaFab = 1, f0 = 0.005, E = 207e9, nu = 0.3, gammaM, gammaSc }) {
  fin('MSd', MSd); fin('SSd', SSd); fin('pi', pi); fin('pe', pe);
  const f = dnvResistanceFactors(safetyClass);
  const gm = gammaM === undefined ? f.gammaM : gammaM;
  const gsc = gammaSc === undefined ? f.gammaSc : gammaSc;
  const g = gm * gsc;
  const { Sp, Mp } = plasticCapacities({ D, t: t2, fy });
  const { pb } = burstPressure({ D, t: t2, fy, fu });
  const beta = (60 - D / t2) / 90;
  const alphaC = (1 - beta) + beta * (fu / fy);
  const warnings = [];
  const r = D / t2;
  if (r < 15 || r > 45) warnings.push('D/t2 is outside 15–45');
  if (Math.abs(SSd) / Sp >= 0.4) warnings.push('|SSd|/Sp is not below 0.4');
  const momentTerm = (g * Math.abs(MSd)) / (alphaC * Mp);
  const axialTerm = ((g * SSd) / (alphaC * Sp)) ** 2;
  let pressureTerm;
  let alphaP = null;
  let pc = null;
  let mode;
  if (pi > pe) {
    mode = 'internal';
    const q = (pi - pe) / pb;
    alphaP = q < 2 / 3 ? 1 - beta : 1 - 3 * beta * (1 - q);
    pressureTerm = (alphaP * q) / alphaC;
  } else {
    mode = 'external';
    const pmin = pMin === undefined ? pi : pMin;
    pc = collapsePressure({ D, t: t2, fy, E, nu, alphaFab, f0 }).pc;
    pressureTerm = (g * (pe - pmin)) / pc;
  }
  const uc = (momentTerm + axialTerm) ** 2 + pressureTerm ** 2;
  const rem = 1 - pressureTerm ** 2;
  const allowableMoment = rem > 0 ? Math.max(0, (Math.sqrt(rem) - axialTerm) * alphaC * Mp / g) : 0;
  return { mode, uc, momentTerm, axialTerm, pressureTerm, alphaC, alphaP, beta, Sp, Mp, pb, pc, allowableMoment, warnings };
}

/**
 * Local buckling, displacement-controlled condition (compressive strain from bending and axial effects).
 *
 * Internal overpressure:  ε_Sd ≤ ε_c(t2, p_min − p_e)/γ_ε,
 *   ε_c(t, Δp) = 0.78 (t/D − 0.01)(1 + 5.75 Δp/p_b(t)) α_h^−1.5 α_gw
 * External overpressure:  [ε_Sd/(ε_c(t2, 0)/γ_ε)]^0.8 + (p_e − p_min)/(p_c(t2)/(γm γSC)) ≤ 1
 *
 * Reference: DNV-ST-F101 (2021) Eq. 5.30–5.32; γ_ε from Table 5-8.
 * Validity: D/t2 ≤ 45, p_i ≥ p_e for the first form. α_h is the yield-to-tensile ratio (R_t0.5/R_m)max of
 * the pipe material (0.93 unless otherwise specified); α_gw is the girth-weld factor (1.0 for D/t2 ≤ 20,
 * decreasing for higher D/t2 — taken from the standard's figure, pass it explicitly).
 *
 * @param {object} o
 * @param {number} o.D [m]
 * @param {number} o.t2 [m]
 * @param {number} o.fy [Pa]
 * @param {number} o.fu [Pa]
 * @param {number} o.epsSd design compressive strain, positive [-]
 * @param {number} [o.pMin=0] minimum internal pressure sustained together with the strain [Pa]
 * @param {number} [o.pe=0] external pressure [Pa]
 * @param {number} [o.alphaH=0.93] [-]
 * @param {number} [o.alphaGw=1] [-]
 * @param {'low'|'medium'|'high'} [o.safetyClass='medium']
 * @param {number} [o.gammaEps] override of the strain resistance factor [-]
 * @param {number} [o.alphaFab=1] [-]
 * @param {number} [o.f0=0.005] [-]
 * @param {number} [o.E=207e9] [Pa]
 * @param {number} [o.nu=0.3] [-]
 * @returns {{mode:'internal'|'external', epsC:number, epsRd:number, uc:number, strainTerm:number,
 *   pressureTerm:number, pc:(number|null)}} strains [-]
 */
export function localBucklingDisplacementControlled({ D, t2, fy, fu, epsSd, pMin = 0, pe = 0, alphaH = 0.93, alphaGw = 1,
  safetyClass = 'medium', gammaEps, alphaFab = 1, f0 = 0.005, E = 207e9, nu = 0.3 }) {
  nonNeg('epsSd', epsSd); pos('alphaH', alphaH); pos('alphaGw', alphaGw);
  const f = dnvResistanceFactors(safetyClass);
  const ge = gammaEps === undefined ? f.gammaEps : pos('gammaEps', gammaEps);
  const { pb } = burstPressure({ D, t: t2, fy, fu });
  if (t2 / D <= 0.01) throw new RangeError('t2/D must exceed 0.01 for the strain capacity formula');
  const epsCof = (dp) => 0.78 * (t2 / D - 0.01) * (1 + (5.75 * dp) / pb) * alphaH ** -1.5 * alphaGw;
  if (pMin >= pe) {
    const epsC = epsCof(pMin - pe);
    const epsRd = epsC / ge;
    return { mode: 'internal', epsC, epsRd, uc: epsSd / epsRd, strainTerm: epsSd / epsRd, pressureTerm: 0, pc: null };
  }
  const epsC = epsCof(0);
  const epsRd = epsC / ge;
  const pc = collapsePressure({ D, t: t2, fy, E, nu, alphaFab, f0 }).pc;
  const strainTerm = (epsSd / epsRd) ** 0.8;
  const pressureTerm = (pe - pMin) / (pc / (f.gammaM * f.gammaSc));
  return { mode: 'external', epsC, epsRd, uc: strainTerm + pressureTerm, strainTerm, pressureTerm, pc };
}

/**
 * Thin-wall stresses and the equivalent-stress (allowable stress design) check
 * σ_e = √(σ_h² + σ_l² − σ_h σ_l + 3τ²) ≤ η f_y, evaluated at the tension and compression bending fibres.
 * σ_h = (p_i − p_e)(D − t)/(2t); σ_l = N/A_s ± M D/(2I), with the true wall force N.
 *
 * Reference: von Mises yield criterion; format of DNV-ST-F101 (2021) Sec. 5 (ASD alternative) and
 * ISO 13623 §6.4. The usage factor η is code- and location-dependent (0.72–1.0 in common pipeline codes)
 * and must be supplied.
 *
 * @param {object} o
 * @param {number} o.D [m]
 * @param {number} o.t [m]
 * @param {number} o.fy yield strength [Pa]
 * @param {number} o.eta usage factor [-]
 * @param {number} [o.pi=0] internal pressure [Pa]
 * @param {number} [o.pe=0] external pressure [Pa]
 * @param {number} [o.N=0] true axial wall force, tension positive [N]
 * @param {number} [o.M=0] bending moment [N·m]
 * @param {number} [o.T=0] torque [N·m]
 * @returns {{sigmaH:number, sigmaLTension:number, sigmaLCompression:number, tau:number, sigmaE:number, allowable:number, uc:number}} stresses [Pa]
 */
export function equivalentStressCheck({ D, t, fy, eta, pi = 0, pe = 0, N = 0, M = 0, T = 0 }) {
  pos('D', D); pos('t', t); pos('fy', fy); pos('eta', eta);
  const Di = D - 2 * t;
  const As = (Math.PI / 4) * (D * D - Di * Di);
  const I = (Math.PI / 64) * (D ** 4 - Di ** 4);
  const sigmaH = ((pi - pe) * (D - t)) / (2 * t);
  const sb = (Math.abs(M) * D) / (2 * I);
  const tau = (Math.abs(T) * D) / (4 * I);
  const vm = (sl) => Math.sqrt(sigmaH * sigmaH + sl * sl - sigmaH * sl + 3 * tau * tau);
  const a = N / As + sb;
  const b = N / As - sb;
  const sigmaE = Math.max(vm(a), vm(b));
  return { sigmaH, sigmaLTension: a, sigmaLCompression: b, tau, sigmaE, allowable: eta * fy, uc: sigmaE / (eta * fy) };
}

/**
 * ASME B31.8 basic design factor F by location class.
 *
 * Reference: ASME B31.8 Table 841.1.6-1: Class 1 Div. 1 0.80; Class 1 Div. 2 0.72; Class 2 0.60; Class 3 0.50;
 * Class 4 0.40. Offshore (Table A842.2.2-1): hoop 0.72 (pipeline) / 0.50 (riser and platform piping),
 * longitudinal 0.80, combined 0.90.
 *
 * @param {'1-1'|'1-2'|'2'|'3'|'4'|'offshore'|'offshore-riser'} locationClass
 * @returns {number} hoop design factor F [-]
 */
export function b318DesignFactor(locationClass) {
  const m = { '1-1': 0.8, '1-2': 0.72, 1: 0.72, 2: 0.6, 3: 0.5, 4: 0.4, offshore: 0.72, 'offshore-riser': 0.5 };
  const v = m[String(locationClass)];
  if (v === undefined) throw new RangeError('locationClass must be one of 1-1, 1-2, 2, 3, 4, offshore, offshore-riser');
  return v;
}

/**
 * ASME B31.8 temperature de-rating factor T for steel pipe.
 *
 * Reference: ASME B31.8 Table 841.1.8-1: 1.000 up to 121 °C (250 °F), 0.967 at 149 °C, 0.933 at 177 °C,
 * 0.900 at 204 °C, 0.867 at 232 °C; linear interpolation between.
 *
 * @param {number} tempC [°C], up to 232
 * @returns {number} T [-]
 */
export function b318TemperatureFactor(tempC) {
  fin('tempC', tempC);
  const xs = [250, 300, 350, 400, 450];
  const ys = [1.0, 0.967, 0.933, 0.9, 0.867];
  const tf = tempC * 1.8 + 32;
  if (tf <= xs[0]) return 1;
  if (tf > xs[4] + 1e-9) throw new RangeError('tempC exceeds 232 °C (450 °F)');
  for (let i = 1; i < 5; i++) if (tf <= xs[i] + 1e-9) return ys[i - 1] + ((ys[i] - ys[i - 1]) * (tf - xs[i - 1])) / (xs[i] - xs[i - 1]);
  return ys[4];
}

/**
 * ASME B31.4 / B31.8 hoop, longitudinal and combined stress checks.
 *
 *   hoop:          S_H = (P_i − P_e) D/(2t) ≤ F1 · S · E_j · T
 *   longitudinal:  restrained   S_L = ν S_H − E α ΔT + N_ext/A_s ± M/Z
 *                  unrestrained S_L = S_H/2 + N_ext/A_s ± M/Z                       |S_L| ≤ F2 · S · T
 *   combined:      Tresca 2·√(((S_L − S_H)/2)² + τ²)  or  von Mises √(S_H² − S_H S_L + S_L² + 3τ²) ≤ F3 · S · T
 *
 * Reference: ASME B31.8 §841.1.1, §833 and §A842.2 (offshore: F1 = 0.72, F2 = 0.80, F3 = 0.90);
 * ASME B31.4 §402, §403.2 and Table A402.3.5-1 (offshore: 0.72 / 0.80 / 0.90). Design factors are arguments
 * so that onshore location classes (b318DesignFactor) or riser values can be used.
 * Validity: thin-wall (D/t ≳ 20; the codes permit the mid-wall diameter for thicker pipe).
 *
 * @param {object} o
 * @param {number} o.D [m]
 * @param {number} o.t pressure-design wall thickness (nominal minus allowances) [m]
 * @param {number} o.smys [Pa]
 * @param {number} o.pi internal design pressure [Pa]
 * @param {number} [o.pe=0] external pressure [Pa]
 * @param {number} [o.F1=0.72] hoop design factor [-]
 * @param {number} [o.F2=0.80] longitudinal design factor [-]
 * @param {number} [o.F3=0.90] combined design factor [-]
 * @param {number} [o.jointFactor=1] longitudinal weld joint factor E_j [-]
 * @param {number} [o.tempFactor=1] temperature de-rating factor T [-]
 * @param {boolean} [o.restrained=true] fully axially restrained (buried / long seabed line) or free
 * @param {number} [o.deltaT=0] temperature rise above installation [K]
 * @param {number} [o.E=207e9] [Pa]
 * @param {number} [o.alpha=11.7e-6] [1/K]
 * @param {number} [o.nu=0.3] [-]
 * @param {number} [o.axialForce=0] additional external axial force, tension positive [N]
 * @param {number} [o.M=0] bending moment [N·m]
 * @param {number} [o.T=0] torque [N·m]
 * @param {'tresca'|'vonMises'} [o.criterion='tresca']
 * @returns {{SH:number, SLmax:number, SLmin:number, Scombined:number, allowHoop:number, allowLong:number,
 *   allowCombined:number, ucHoop:number, ucLong:number, ucCombined:number, tRequiredHoop:number}}
 *   stresses [Pa]; tRequiredHoop pressure-design thickness for the hoop criterion [m]
 */
export function asmeStressCheck({ D, t, smys, pi, pe = 0, F1 = 0.72, F2 = 0.8, F3 = 0.9, jointFactor = 1, tempFactor = 1,
  restrained = true, deltaT = 0, E = 207e9, alpha = 11.7e-6, nu = 0.3, axialForce = 0, M = 0, T = 0, criterion = 'tresca' }) {
  pos('D', D); pos('t', t); pos('smys', smys); fin('pi', pi); fin('pe', pe);
  pos('F1', F1); pos('F2', F2); pos('F3', F3); pos('jointFactor', jointFactor); pos('tempFactor', tempFactor);
  const Di = D - 2 * t;
  const As = (Math.PI / 4) * (D * D - Di * Di);
  const I = (Math.PI / 64) * (D ** 4 - Di ** 4);
  const SH = ((pi - pe) * D) / (2 * t);
  const sb = (Math.abs(M) * D) / (2 * I);
  const tau = (Math.abs(T) * D) / (4 * I);
  const base = (restrained ? nu * SH - E * alpha * deltaT : SH / 2) + axialForce / As;
  const SLmax = base + sb;
  const SLmin = base - sb;
  const comb = criterion === 'vonMises'
    ? (sl) => Math.sqrt(SH * SH - SH * sl + sl * sl + 3 * tau * tau)
    : (sl) => 2 * Math.sqrt(((sl - SH) / 2) ** 2 + tau * tau);
  if (criterion !== 'vonMises' && criterion !== 'tresca') throw new RangeError('criterion must be "tresca" or "vonMises"');
  const Scombined = Math.max(comb(SLmax), comb(SLmin));
  const allowHoop = F1 * smys * jointFactor * tempFactor;
  const allowLong = F2 * smys * tempFactor;
  const allowCombined = F3 * smys * tempFactor;
  return {
    SH, SLmax, SLmin, Scombined, allowHoop, allowLong, allowCombined,
    ucHoop: SH / allowHoop, ucLong: Math.max(Math.abs(SLmax), Math.abs(SLmin)) / allowLong, ucCombined: Scombined / allowCombined,
    tRequiredHoop: ((pi - pe) * D) / (2 * allowHoop),
  };
}

/**
 * Thick-wall (Lamé) stresses in a cylinder under internal and external pressure.
 * σ_h(r) = A + B/r², σ_r(r) = A − B/r², A = (p_i r_i² − p_e r_o²)/(r_o² − r_i²), B = (p_i − p_e) r_i² r_o²/(r_o² − r_i²);
 * capped-end axial stress σ_a = A.
 *
 * Reference: Lamé (1852); Timoshenko, Strength of Materials Part II.
 * Validity: linear elastic, plane cross-sections, any D/t.
 *
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.t wall thickness [m]
 * @param {number} [o.pi=0] internal pressure [Pa]
 * @param {number} [o.pe=0] external pressure [Pa]
 * @param {number} [o.r] radius at which to evaluate [m] (default: inner surface)
 * @returns {{sigmaH:number, sigmaR:number, sigmaAxialCapped:number, sigmaHInner:number, sigmaHOuter:number, vonMisesCapped:number}} [Pa]
 */
export function lameStresses({ D, t, pi = 0, pe = 0, r }) {
  pos('D', D); pos('t', t); fin('pi', pi); fin('pe', pe);
  const ro = D / 2;
  const ri = ro - t;
  if (!(ri > 0)) throw new RangeError('t must be smaller than D/2');
  const rr = r === undefined ? ri : r;
  if (!(rr >= ri - 1e-12 && rr <= ro + 1e-12)) throw new RangeError('r must lie within the wall');
  const den = ro * ro - ri * ri;
  const A = (pi * ri * ri - pe * ro * ro) / den;
  const B = ((pi - pe) * ri * ri * ro * ro) / den;
  const sh = A + B / (rr * rr);
  const sr = A - B / (rr * rr);
  const vm = Math.sqrt(0.5 * ((sh - sr) ** 2 + (sr - A) ** 2 + (A - sh) ** 2));
  return { sigmaH: sh, sigmaR: sr, sigmaAxialCapped: A, sigmaHInner: A + B / (ri * ri), sigmaHOuter: A + B / (ro * ro), vonMisesCapped: vm };
}

/**
 * Hydrostatic (strength) test pressure and the hoop stress it produces.
 *
 * Reference: ASME B31.4 §437.4.1 and ASME B31.8 §A847.2 (offshore): not less than 1.25 × the design pressure
 * / MAOP. Onshore B31.8 factors vary by location class (Table 841.3.2-1) and DNV-ST-F101 uses the system
 * test factor α_spt on the local incidental pressure — pass the applicable `factor`.
 *
 * @param {object} o
 * @param {number} o.pDesign design pressure / MAOP [Pa]
 * @param {number} [o.factor=1.25] test factor [-]
 * @param {number} [o.D] outer diameter [m] (optional, for the hoop stress)
 * @param {number} [o.t] wall thickness [m]
 * @param {number} [o.smys] [Pa]
 * @param {number} [o.pe=0] external pressure during test [Pa]
 * @returns {{pTest:number, hoopStress:(number|null), hoopFractionSmys:(number|null)}}
 */
export function hydrotestPressure({ pDesign, factor = 1.25, D, t, smys, pe = 0 }) {
  pos('pDesign', pDesign); pos('factor', factor);
  const pTest = factor * pDesign;
  let hoopStress = null;
  let frac = null;
  if (D !== undefined && t !== undefined) {
    hoopStress = ((pTest - pe) * pos('D', D)) / (2 * pos('t', t));
    if (smys !== undefined) frac = hoopStress / pos('smys', smys);
  }
  return { pTest, hoopStress, hoopFractionSmys: frac };
}

/**
 * Bending strain and bend radius of a pipe: ε = D/(2R); elastic limit radius R = E D/(2 σ_allow).
 *
 * Reference: Euler–Bernoulli kinematics.
 *
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} [o.R] bend radius [m] (gives the strain)
 * @param {number} [o.epsAllow] allowable bending strain [-] (gives the minimum radius)
 * @param {number} [o.sigmaAllow] allowable elastic bending stress [Pa] (gives the minimum elastic radius)
 * @param {number} [o.E=207e9] [Pa]
 * @returns {{strain:(number|null), minRadiusStrain:(number|null), minRadiusElastic:(number|null)}} radii [m]
 */
export function bendingStrainRadius({ D, R, epsAllow, sigmaAllow, E = 207e9 }) {
  pos('D', D);
  return {
    strain: R === undefined ? null : D / (2 * pos('R', R)),
    minRadiusStrain: epsAllow === undefined ? null : D / (2 * pos('epsAllow', epsAllow)),
    minRadiusElastic: sigmaAllow === undefined ? null : (E * D) / (2 * pos('sigmaAllow', sigmaAllow)),
  };
}

/**
 * Ovalisation of a thin elastic tube in pure bending (Brazier effect).
 * ζ = (1 − ν²) r⁴ κ²/t² is the relative flattening of the radius in the plane of bending;
 * ovality f = (D_max − D_min)/D = 2ζ; M = π E r³ t κ (1 − 3ζ/2); limit moment at ζ = 2/9:
 * M_Brazier = (2√2/9) π E r t²/√(1 − ν²).
 *
 * Reference: Brazier, L.G. (1927), Proc. R. Soc. A 116, 104–114.
 * Validity: thin wall, elastic, long tube; the first-order solution is meaningful for ζ ≤ 2/9.
 *
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.t [m]
 * @param {number} o.curvature κ [1/m]
 * @param {number} [o.E=207e9] [Pa]
 * @param {number} [o.nu=0.3] [-]
 * @returns {{zeta:number, ovality:number, moment:number, brazierMoment:number, criticalCurvature:number, beyondLimit:boolean}}
 *   moment [N·m]; criticalCurvature [1/m]
 */
export function brazierOvalisation({ D, t, curvature, E = 207e9, nu = 0.3 }) {
  pos('D', D); pos('t', t); nonNeg('curvature', curvature);
  const r = (D - t) / 2;
  const zeta = ((1 - nu * nu) * r ** 4 * curvature * curvature) / (t * t);
  return {
    zeta, ovality: 2 * zeta,
    moment: Math.PI * E * r ** 3 * t * curvature * (1 - 1.5 * zeta),
    brazierMoment: ((2 * Math.SQRT2) / 9) * Math.PI * E * r * t * t / Math.sqrt(1 - nu * nu),
    criticalCurvature: (Math.SQRT2 / 3) * t / (r * r * Math.sqrt(1 - nu * nu)),
    beyondLimit: zeta > 2 / 9,
  };
}

/**
 * Ovalisation (flattening) under bending strain combined with external pressure:
 * f_0' = [f_0 + c1 (1 + D/(c2 t)) (2 ε D/t)²]/(1 − p_e/p_c), to be kept below 3 % in general.
 *
 * Reference: DNV-ST-F101 (2021) Eq. 5.34 (form of Murphey & Langner). The coefficients c1 = 0.030 and
 * c2 = 120 are as recalled from the standard and are exposed so they can be set from the edition in force.
 *
 * @param {object} o
 * @param {number} o.D [m]
 * @param {number} o.t [m]
 * @param {number} o.strain bending strain [-]
 * @param {number} [o.f0=0.005] initial ovality [-]
 * @param {number} [o.pe=0] external overpressure [Pa]
 * @param {number} [o.pc] collapse pressure [Pa] (required when pe > 0)
 * @param {number} [o.c1=0.030] [-]
 * @param {number} [o.c2=120] [-]
 * @param {number} [o.limit=0.03] acceptance limit [-]
 * @returns {{ovality:number, limit:number, uc:number}}
 */
export function ovalisationBendingPressure({ D, t, strain, f0 = 0.005, pe = 0, pc, c1 = 0.03, c2 = 120, limit = 0.03 }) {
  pos('D', D); pos('t', t); nonNeg('strain', strain); nonNeg('f0', f0); nonNeg('pe', pe);
  let amp = 1;
  if (pe > 0) {
    pos('pc', pc);
    if (pe >= pc) throw new RangeError('pe must be below the collapse pressure pc');
    amp = 1 / (1 - pe / pc);
  }
  const ovality = (f0 + c1 * (1 + D / (c2 * t)) * ((2 * strain * D) / t) ** 2) * amp;
  return { ovality, limit, uc: ovality / limit };
}

/**
 * Pipe bend (elbow) flexibility characteristic, flexibility factor and stress intensification factor:
 * h = t R/r², k = 1.65/h, i = 0.9/h^(2/3) (not less than 1), r = mean radius of the pipe wall.
 *
 * Reference: ASME B31.4 Table 402.1-1 / ASME B31.8 Table E-1 (welding elbow or pipe bend).
 * Validity: no pressure stiffening correction, no flange correction.
 *
 * @param {object} o
 * @param {number} o.t nominal wall thickness [m]
 * @param {number} o.R bend radius to the pipe centreline [m]
 * @param {number} o.r mean radius of the pipe cross-section, (D − t)/2 [m]
 * @returns {{h:number, k:number, i:number}} all [-]
 */
export function elbowFactors({ t, R, r }) {
  pos('t', t); pos('R', R); pos('r', r);
  const h = (t * R) / (r * r);
  return { h, k: Math.max(1, 1.65 / h), i: Math.max(1, 0.9 / h ** (2 / 3)) };
}

/**
 * Wall thickness change and pressure hoop-stress factors of a pipe bend.
 *
 * Thinning estimate (kinematic, constant volume with the neutral axis on the centreline and no hoop strain):
 * t_extrados = t R/(R + r_o), t_intrados = t R/(R − r_o), r_o = D/2. Hoop stress intensification of a torus
 * under pressure (Lorenz factors): intrados (4R/D − 1)/(4R/D − 2), extrados (4R/D + 1)/(4R/D + 2).
 *
 * Reference: first-principles kinematics for the thinning; ASME B31.3 §304.2.1 for the Lorenz factors.
 * Validity: R > D/2; induction-bend thinning in practice depends on the process — treat as an estimate.
 *
 * @param {object} o
 * @param {number} o.D outer diameter [m]
 * @param {number} o.t wall thickness of the mother pipe [m]
 * @param {number} o.R bend centreline radius [m]
 * @returns {{tExtrados:number, tIntrados:number, thinningFraction:number, lorenzIntrados:number, lorenzExtrados:number}}
 */
export function bendWallThickness({ D, t, R }) {
  pos('D', D); pos('t', t); pos('R', R);
  if (R <= D / 2) throw new RangeError('R must exceed D/2');
  const ro = D / 2;
  const q = (4 * R) / D;
  return {
    tExtrados: (t * R) / (R + ro), tIntrados: (t * R) / (R - ro), thinningFraction: ro / (R + ro),
    lorenzIntrados: (q - 1) / (q - 2), lorenzExtrados: (q + 1) / (q + 2),
  };
}
