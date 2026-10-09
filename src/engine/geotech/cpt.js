// Cone penetration test (CPT/CPTU) interpretation and SPT corrections.
// All stresses in Pa. pa is the atmospheric reference pressure (default 100 kPa as used by the correlations).

const PA = 100e3;

function pos(name, v, allowZero = false) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < 0 || (!allowZero && v === 0)) {
    throw new RangeError(`${name}: must be a finite number ${allowZero ? '>= 0' : '> 0'} (got ${v})`);
  }
}

/**
 * Corrected cone resistance qt = qc + u2 (1 − a), with a the cone net area ratio.
 * Reference: Campanella, Gillespie & Robertson (1982); Lunne, Robertson & Powell (1997) CPT in Geotechnical
 * Practice. Validity: u2 measured behind the cone shoulder.
 * @param {{ qc: number, u2: number, a: number }} o qc [Pa] measured cone resistance; u2 [Pa] pore pressure at
 *   the shoulder (absolute above atmospheric, may be negative in dilative soil); a [-] net area ratio (0.5–0.9)
 * @returns {number} qt [Pa]
 */
export function correctedConeResistance({ qc, u2, a }) {
  pos('qc', qc, true);
  if (!Number.isFinite(u2)) throw new RangeError('u2: must be finite');
  if (!(a > 0 && a <= 1)) throw new RangeError('a: must be in (0, 1]');
  return qc + u2 * (1 - a);
}

/**
 * Basic CPTU ratios: friction ratio Rf = fs/qt · 100 %, pore pressure ratio Bq = (u2 − u0)/(qt − σv0),
 * normalised cone resistance Qt = (qt − σv0)/σ′v0 and normalised friction ratio Fr = fs/(qt − σv0) · 100 %.
 * Reference: Robertson (1990) Can. Geotech. J. 27(1) 151–158; Lunne et al. (1997). Validity: qt > σv0.
 * @param {{ qt: number, fs: number, sigmaV0: number, sigmaV0Eff: number, u2?: number, u0?: number }} o qt [Pa];
 *   fs [Pa] sleeve friction; sigmaV0 [Pa] total vertical stress; sigmaV0Eff [Pa] effective vertical stress;
 *   u2 [Pa]; u0 [Pa] in-situ pore pressure
 * @returns {{ Rf: number, Bq: number|null, Qt: number, Fr: number, qnet: number }} Rf [%]; Bq [-] (null without
 *   u2, u0); Qt [-]; Fr [%]; qnet [Pa]
 */
export function cptRatios({ qt, fs, sigmaV0, sigmaV0Eff, u2, u0 }) {
  pos('qt', qt); pos('fs', fs, true); pos('sigmaV0', sigmaV0, true); pos('sigmaV0Eff', sigmaV0Eff);
  const qnet = qt - sigmaV0;
  if (!(qnet > 0)) throw new RangeError('qt: must exceed the total vertical stress sigmaV0');
  return { Rf: (100 * fs) / qt, Bq: u2 !== undefined && u0 !== undefined ? (u2 - u0) / qnet : null, Qt: qnet / sigmaV0Eff, Fr: (100 * fs) / qnet, qnet };
}

/**
 * Soil behaviour type zone (Robertson normalised chart) from the index Ic:
 * Ic > 3.60 zone 2; 2.95–3.60 zone 3; 2.60–2.95 zone 4; 2.05–2.60 zone 5; 1.31–2.05 zone 6; < 1.31 zone 7.
 * Zones 1, 8 and 9 are not identified by Ic.
 * Reference: Robertson & Wride (1998) Can. Geotech. J. 35(3); Robertson (2009).
 * @param {number} Ic soil behaviour type index [-]
 * @returns {{ zone: number, description: string, group: string }} group 'clay' | 'silt' | 'sand' (coarse
 *   grouping used for layering: zones 2–3 clay, 4 silt, 5–7 sand)
 */
export function sbtZone(Ic) {
  if (!Number.isFinite(Ic)) throw new RangeError('Ic: must be finite');
  if (Ic > 3.6) return { zone: 2, description: 'Organic soils – clay', group: 'clay' };
  if (Ic > 2.95) return { zone: 3, description: 'Clays – silty clay to clay', group: 'clay' };
  if (Ic > 2.6) return { zone: 4, description: 'Silt mixtures – clayey silt to silty clay', group: 'silt' };
  if (Ic > 2.05) return { zone: 5, description: 'Sand mixtures – silty sand to sandy silt', group: 'sand' };
  if (Ic > 1.31) return { zone: 6, description: 'Sands – clean sand to silty sand', group: 'sand' };
  return { zone: 7, description: 'Gravelly sand to dense sand', group: 'sand' };
}

/**
 * Soil behaviour type index with the stress exponent iteration of Robertson (2009):
 *   Qtn = [(qt − σv0)/pa] (pa/σ′v0)^n, Fr = fs/(qt − σv0) · 100 %,
 *   Ic = sqrt[(3.47 − log10 Qtn)² + (log10 Fr + 1.22)²], n = 0.381 Ic + 0.05 (σ′v0/pa) − 0.15 <= 1.
 * Reference: Robertson (2009) Can. Geotech. J. 46(11) 1337–1355. Validity: uncemented, young soils.
 * @param {{ qt: number, fs: number, sigmaV0: number, sigmaV0Eff: number, pa?: number }} o qt, fs, sigmaV0,
 *   sigmaV0Eff [Pa]; pa [Pa] (default 100e3)
 * @returns {{ Ic: number, Qtn: number, Fr: number, n: number, iterations: number, zone: number,
 *   description: string, group: string }} Ic, Qtn, n [-]; Fr [%]
 */
export function soilBehaviourType({ qt, fs, sigmaV0, sigmaV0Eff, pa = PA }) {
  pos('qt', qt); pos('fs', fs); pos('sigmaV0', sigmaV0, true); pos('sigmaV0Eff', sigmaV0Eff); pos('pa', pa);
  const qnet = qt - sigmaV0;
  if (!(qnet > 0)) throw new RangeError('qt: must exceed the total vertical stress sigmaV0');
  const Fr = (100 * fs) / qnet;
  let n = 1, Ic = 0, Qtn = 0, it = 0;
  for (; it < 100; it++) {
    Qtn = (qnet / pa) * Math.pow(pa / sigmaV0Eff, n);
    Ic = Math.hypot(3.47 - Math.log10(Qtn), Math.log10(Fr) + 1.22);
    const nNew = Math.min(0.381 * Ic + 0.05 * (sigmaV0Eff / pa) - 0.15, 1);
    if (Math.abs(nNew - n) < 1e-6) { n = nNew; break; }
    n = nNew;
  }
  Qtn = (qnet / pa) * Math.pow(pa / sigmaV0Eff, n);
  Ic = Math.hypot(3.47 - Math.log10(Qtn), Math.log10(Fr) + 1.22);
  return { Ic, Qtn, Fr, n, iterations: it + 1, ...sbtZone(Ic) };
}

/**
 * Relative density of clean sand from cone resistance.
 * 'baldi': Dr = (1/C2) ln[qc / (C0 σ′^C1)] with qc and σ′ in kPa; C0 = 157, C1 = 0.55, C2 = 2.41 for normally
 *   consolidated Ticino sand with σ′ = σ′v0 (Baldi et al. 1986).
 * 'jamiolkowski': Dr = 0.268 ln[(qt/pa) / sqrt(σ′v0/pa)] − 0.675 (dry-sand calibration-chamber fit,
 *   Jamiolkowski et al. 2003).
 * The constants are arguments because they depend on sand compressibility and stress history.
 * Reference: Baldi et al. (1986) Proc. 4th Int. Geotech. Seminar, Singapore; Jamiolkowski, Lo Presti & Manassero
 * (2003) ASCE GSP 119. Validity: clean, uncemented, unaged silica sands; result not clipped to 0..1.
 * @param {{ qc: number, sigmaV0Eff: number, method?: string, C0?: number, C1?: number, C2?: number, a?: number,
 *   b?: number, pa?: number }} o qc [Pa] (qt for sand); sigmaV0Eff [Pa]; method 'baldi'|'jamiolkowski';
 *   C0, C1, C2 [-] Baldi constants; a, b [-] Jamiolkowski constants (0.268, 0.675); pa [Pa]
 * @returns {number} Dr [-] as a fraction
 */
export function relativeDensity({ qc, sigmaV0Eff, method = 'baldi', C0 = 157, C1 = 0.55, C2 = 2.41, a = 0.268, b = 0.675, pa = PA }) {
  pos('qc', qc); pos('sigmaV0Eff', sigmaV0Eff);
  if (method === 'baldi') return Math.log(qc / 1000 / (C0 * Math.pow(sigmaV0Eff / 1000, C1))) / C2;
  if (method === 'jamiolkowski') return a * Math.log(qc / pa / Math.sqrt(sigmaV0Eff / pa)) - b;
  throw new RangeError("method: must be 'baldi' or 'jamiolkowski'");
}

/**
 * Peak friction angle of sand from cone resistance.
 * 'robertson-campanella': tan φ′ = (1/2.68) [log10(qc/σ′v0) + 0.29] (Robertson & Campanella 1983);
 * 'kulhawy-mayne': φ′ = 17.6° + 11° log10[(qt/pa) / sqrt(σ′v0/pa)] (Kulhawy & Mayne 1990).
 * Reference: Robertson & Campanella (1983) Can. Geotech. J. 20(4); Kulhawy & Mayne (1990) EPRI EL-6800.
 * Validity: clean, uncemented, moderately compressible quartz sands.
 * @param {{ qc: number, sigmaV0Eff: number, method?: string, pa?: number }} o qc [Pa]; sigmaV0Eff [Pa];
 *   method 'kulhawy-mayne' (default) | 'robertson-campanella'; pa [Pa]
 * @returns {number} φ′ [deg]
 */
export function frictionAngleDeg({ qc, sigmaV0Eff, method = 'kulhawy-mayne', pa = PA }) {
  pos('qc', qc); pos('sigmaV0Eff', sigmaV0Eff);
  if (method === 'kulhawy-mayne') return 17.6 + 11 * Math.log10(qc / pa / Math.sqrt(sigmaV0Eff / pa));
  if (method === 'robertson-campanella') return (Math.atan((Math.log10(qc / sigmaV0Eff) + 0.29) / 2.68) * 180) / Math.PI;
  throw new RangeError("method: must be 'kulhawy-mayne' or 'robertson-campanella'");
}

/**
 * Undrained shear strength from net cone resistance: su = (qt − σv0) / Nkt.
 * Reference: Lunne, Robertson & Powell (1997). Validity: clays; Nkt is site specific, typically 10–18
 * (often 14 for a first estimate) and must be calibrated against laboratory strengths.
 * @param {{ qt: number, sigmaV0: number, Nkt: number }} o qt [Pa]; sigmaV0 [Pa]; Nkt [-] cone factor (required)
 * @returns {number} su [Pa]
 */
export function undrainedShearStrength({ qt, sigmaV0, Nkt }) {
  pos('qt', qt); pos('sigmaV0', sigmaV0, true); pos('Nkt', Nkt);
  return Math.max(qt - sigmaV0, 0) / Nkt;
}

/**
 * Overconsolidation ratio of clay: OCR = k (qt − σv0)/σ′v0.
 * Reference: Kulhawy & Mayne (1990) EPRI EL-6800 (k = 0.33; reported range 0.2–0.5). Validity: intact,
 * non-fissured clays.
 * @param {{ qt: number, sigmaV0: number, sigmaV0Eff: number, k?: number }} o stresses [Pa]; k [-] (default 0.33)
 * @returns {number} OCR [-] (not clipped; values below 1 indicate under-prediction or underconsolidation)
 */
export function overconsolidationRatio({ qt, sigmaV0, sigmaV0Eff, k = 0.33 }) {
  pos('qt', qt); pos('sigmaV0', sigmaV0, true); pos('sigmaV0Eff', sigmaV0Eff); pos('k', k);
  return (k * Math.max(qt - sigmaV0, 0)) / sigmaV0Eff;
}

/**
 * Small-strain shear modulus from cone resistance.
 * 'rix-stokoe' (sands): Gmax = A qc^0.25 σ′v0^0.375 with all stresses in kPa and A = 1634 (Rix & Stokoe 1991);
 * 'robertson' (all soils): Vs = [αvs (qt − σv0)/pa]^0.5 m/s with αvs = 10^(0.55 Ic + 1.68), Gmax = ρ Vs²
 *   (Robertson 2009).
 * Reference: Rix & Stokoe (1991) in Calibration Chamber Testing, Elsevier; Robertson (2009) Can. Geotech. J.
 * 46(11). Validity: uncemented Holocene–Pleistocene soils; A, the exponents and the αvs constants are arguments.
 * @param {{ method?: string, qc?: number, qt?: number, sigmaV0?: number, sigmaV0Eff?: number, Ic?: number,
 *   rho?: number, A?: number, nq?: number, ns?: number, c1?: number, c0?: number, pa?: number }} o method
 *   'rix-stokoe' (needs qc, sigmaV0Eff [Pa]) | 'robertson' (needs qt, sigmaV0 [Pa], Ic, rho [kg/m³] saturated
 *   density); A (1634), nq (0.25), ns (0.375); c1 (0.55), c0 (1.68); pa [Pa]
 * @returns {{ Gmax: number, Vs: number|null }} Gmax [Pa]; Vs [m/s] (Robertson method only)
 */
export function smallStrainShearModulus({ method = 'rix-stokoe', qc, qt, sigmaV0, sigmaV0Eff, Ic, rho, A = 1634, nq = 0.25, ns = 0.375, c1 = 0.55, c0 = 1.68, pa = PA }) {
  if (method === 'rix-stokoe') {
    pos('qc', qc); pos('sigmaV0Eff', sigmaV0Eff);
    return { Gmax: 1000 * A * Math.pow(qc / 1000, nq) * Math.pow(sigmaV0Eff / 1000, ns), Vs: null };
  }
  if (method === 'robertson') {
    pos('qt', qt); pos('sigmaV0', sigmaV0, true); pos('rho', rho);
    if (!Number.isFinite(Ic)) throw new RangeError('Ic: must be finite');
    const Vs = Math.sqrt((Math.pow(10, c1 * Ic + c0) * Math.max(qt - sigmaV0, 0)) / pa);
    return { Gmax: rho * Vs * Vs, Vs };
  }
  throw new RangeError("method: must be 'rix-stokoe' or 'robertson'");
}

/**
 * Total unit weight estimate: γ/γw = 0.27 log10(Rf) + 0.36 log10(qt/pa) + 1.236.
 * Reference: Robertson & Cabal (2010) Proc. CPT'10, Huntington Beach. Validity: Rf in %, approximate (about
 * ±1 kN/m³); not for cemented or very organic soils.
 * @param {{ qt: number, Rf: number, gammaW?: number, pa?: number }} o qt [Pa]; Rf [%] friction ratio (> 0);
 *   gammaW [N/m³] (default 9810, as in the source); pa [Pa]
 * @returns {number} total unit weight γ [N/m³]
 */
export function unitWeightEstimate({ qt, Rf, gammaW = 9810, pa = PA }) {
  pos('qt', qt); pos('Rf', Rf); pos('gammaW', gammaW);
  return gammaW * (0.27 * Math.log10(Rf) + 0.36 * Math.log10(qt / pa) + 1.236);
}

/**
 * Interpret a CPTU sounding into a layered soil profile for soilprofile.js.
 * At each reading: unit weight (Robertson & Cabal) → running σv0, σ′v0 (hydrostatic u0 from the mudline) →
 * Ic (Robertson 2009) → group clay / silt / sand by Ic thresholds. Consecutive readings of the same group
 * form a layer; layers thinner than minThickness are merged into the layer above. Layer properties are
 * averages of the point estimates: clay/silt su (with a linear least-squares gradient), OCR; sand φ′, Dr.
 * Reference: see the individual correlations. Validity: offshore sounding starting at the mudline, fully
 * saturated, hydrostatic; screening-level layering to be reviewed against borehole logs.
 * @param {{ z: number[], qc: number[], fs: number[], u2?: number[], a?: number, Nkt: number, gammaW?: number,
 *   minThickness?: number, sandThresholdIc?: number, clayThresholdIc?: number, drMethod?: string,
 *   phiMethod?: string, eps50?: number }} o z [m] depths below mudline, ascending, first > 0; qc, fs [Pa];
 *   u2 [Pa] pore pressure referenced to the water pressure at the mudline, i.e. cone zeroed at the seabed
 *   (optional — without it qt = qc); a [-] net area ratio (default 0.8); Nkt [-] (required);
 *   gammaW [N/m³] (default 10 052); minThickness [m] (default 1); sandThresholdIc (2.6), clayThresholdIc (2.95);
 *   eps50 [-] assigned to cohesive layers (default: not set)
 * @returns {{ layers: object[], points: object }} layers for createProfile (type, top, thickness, gammaSub,
 *   suTop, suGrad, OCR, phiDeg, Dr, Ic); points: per-reading arrays { z, qt, sigmaV0, sigmaV0Eff, Ic, Qtn, Fr,
 *   group, su, phiDeg, Dr, gamma }
 */
export function cptToLayers(o) {
  const { z, qc, fs, u2, a = 0.8, Nkt, gammaW = 10052, minThickness = 1, sandThresholdIc = 2.6, clayThresholdIc = 2.95, drMethod = 'baldi', phiMethod = 'kulhawy-mayne', eps50 } = o;
  if (!Array.isArray(z) || z.length < 2 || qc?.length !== z.length || fs?.length !== z.length) throw new RangeError('z: z, qc and fs must be arrays of equal length >= 2');
  pos('Nkt', Nkt); pos('minThickness', minThickness); pos('gammaW', gammaW);
  if (!(z[0] > 0)) throw new RangeError('z: first depth must be > 0');
  if (!(clayThresholdIc >= sandThresholdIc)) throw new RangeError('clayThresholdIc: must be >= sandThresholdIc');
  const n = z.length;
  const P = { z: z.slice(), qt: [], sigmaV0: [], sigmaV0Eff: [], Ic: [], Qtn: [], Fr: [], group: [], su: [], phiDeg: [], Dr: [], gamma: [] };
  let sv = 0, zPrev = 0;
  for (let i = 0; i < n; i++) {
    if (i > 0 && !(z[i] > z[i - 1])) throw new RangeError('z: must be strictly ascending');
    const u0 = gammaW * z[i];
    const qt = u2 ? correctedConeResistance({ qc: qc[i], u2: u2[i], a }) : qc[i];
    const fsi = Math.max(fs[i], 1e-3 * qt * 0.01); // guard against a zero sleeve reading (Rf floor 0.001 %)
    const gamma = Math.max(unitWeightEstimate({ qt, Rf: (100 * fsi) / qt, gammaW: 9810 }), gammaW * 1.1);
    sv += gamma * (z[i] - zPrev); zPrev = z[i];
    const se = sv - u0;
    const s = soilBehaviourType({ qt, fs: fsi, sigmaV0: sv, sigmaV0Eff: se });
    const group = s.Ic > clayThresholdIc ? 'clay' : s.Ic > sandThresholdIc ? 'silt' : 'sand';
    P.qt.push(qt); P.sigmaV0.push(sv); P.sigmaV0Eff.push(se); P.Ic.push(s.Ic); P.Qtn.push(s.Qtn); P.Fr.push(s.Fr); P.group.push(group); P.gamma.push(gamma);
    P.su.push(group === 'sand' ? NaN : undrainedShearStrength({ qt, sigmaV0: sv, Nkt }));
    P.phiDeg.push(group === 'sand' ? frictionAngleDeg({ qc: qt, sigmaV0Eff: se, method: phiMethod }) : NaN);
    P.Dr.push(group === 'sand' ? relativeDensity({ qc: qt, sigmaV0Eff: se, method: drMethod }) : NaN);
  }
  // runs of equal group; each reading represents the interval up to its depth
  let runs = [];
  for (let i = 0; i < n; i++) {
    const last = runs[runs.length - 1];
    if (last && last.group === P.group[i]) last.to = i; else runs.push({ group: P.group[i], from: i, to: i });
  }
  const topOf = (r) => (r.from === 0 ? 0 : z[r.from - 1]);
  const merged = [];
  for (const r of runs) {
    const prev = merged[merged.length - 1];
    if (prev && z[r.to] - topOf(r) < minThickness) prev.to = r.to; else if (prev && prev.group === r.group) prev.to = r.to; else merged.push({ ...r });
  }
  // a thin first run is absorbed by the run below it
  if (merged.length > 1 && z[merged[0].to] - 0 < minThickness) { merged[1].from = 0; merged.shift(); }
  runs = merged;
  const layers = runs.map((r) => {
    const top = topOf(r), bottom = z[r.to];
    const idx = [];
    for (let i = r.from; i <= r.to; i++) if (P.group[i] === r.group) idx.push(i);
    const mean = (arr) => idx.reduce((s, i) => s + arr[i], 0) / idx.length;
    const layer = { type: r.group, top, thickness: bottom - top, gammaSub: mean(P.gamma) - gammaW, Ic: mean(P.Ic) };
    if (r.group === 'sand') {
      layer.phiDeg = Math.min(Math.max(mean(P.phiDeg), 20), 45);
      layer.Dr = Math.min(Math.max(mean(P.Dr), 0), 1);
    } else {
      // least-squares line su = suTop + suGrad (z − top) through the readings of the layer
      const zm = mean(P.z), sm = mean(P.su);
      let sxx = 0, sxy = 0;
      for (const i of idx) { sxx += (z[i] - zm) ** 2; sxy += (z[i] - zm) * (P.su[i] - sm); }
      let grad = sxx > 0 ? sxy / sxx : 0, suTop = sm - grad * (zm - top);
      if (suTop < 0 || suTop + grad * (bottom - top) < 0) { grad = 0; suTop = sm; }
      layer.suTop = suTop; layer.suGrad = grad;
      layer.OCR = Math.max(1, idx.reduce((s, i) => s + overconsolidationRatio({ qt: P.qt[i], sigmaV0: P.sigmaV0[i], sigmaV0Eff: P.sigmaV0Eff[i] }), 0) / idx.length);
      if (eps50 !== undefined) layer.eps50 = eps50;
    }
    return layer;
  });
  return { layers, points: P };
}

/**
 * SPT energy and procedure correction: N60 = N (ER/60) CB CS CR, and overburden normalisation
 * (N1)60 = CN N60 with CN = (pa/σ′v0)^0.5 <= CNmax.
 * Reference: Skempton (1986) Géotechnique 36(3) for the N60 corrections; Liao & Whitman (1986) J. Geotech. Eng.
 * 112(3) for CN; CN <= 1.7 per Youd et al. (2001). Validity: sands; the factors CB (borehole diameter),
 * CS (sampler liner) and CR (rod length) are inputs from the cited tables.
 * @param {{ N: number, energyRatioPct?: number, CB?: number, CS?: number, CR?: number, sigmaV0Eff: number,
 *   CNmax?: number, pa?: number }} o N [blows/0.3 m] field blow count; energyRatioPct [%] hammer energy ratio
 *   (default 60); CB, CS, CR [-] (default 1); sigmaV0Eff [Pa]; CNmax [-] (default 1.7); pa [Pa]
 * @returns {{ N60: number, CN: number, N160: number }} corrected blow counts [blows/0.3 m]
 */
export function sptCorrection({ N, energyRatioPct = 60, CB = 1, CS = 1, CR = 1, sigmaV0Eff, CNmax = 1.7, pa = PA }) {
  pos('N', N, true); pos('energyRatioPct', energyRatioPct); pos('CB', CB); pos('CS', CS); pos('CR', CR); pos('sigmaV0Eff', sigmaV0Eff); pos('CNmax', CNmax);
  const N60 = N * (energyRatioPct / 60) * CB * CS * CR;
  const CN = Math.min(Math.sqrt(pa / sigmaV0Eff), CNmax);
  return { N60, CN, N160: CN * N60 };
}
