// Axial capacity of driven open- or closed-ended steel pipe piles (static, API main-text method).

import { propertiesAt, scourBaseDepth, verticalEffectiveStress } from './soilprofile.js';


function pos(name, v, allowZero = false) {
  if (typeof v !== 'number' || !Number.isFinite(v) || v < 0 || (!allowZero && v === 0)) {
    throw new RangeError(`${name}: must be a finite number ${allowZero ? '>= 0' : '> 0'} (got ${v})`);
  }
}

/**
 * Design parameters for cohesionless siliceous soil: shaft friction factor β [-], limiting unit skin friction
 * fLim [Pa], end-bearing factor Nq [-], limiting unit end bearing qLim [Pa].
 * β applies to open-ended piles driven unplugged; for full-displacement piles β may be taken 25 % higher.
 * Reference: API RP 2GEO (2011) Table 1 (API RP 2A-WSD 21st ed. Supplement 3, Table 6.4.3-1).
 * Not applicable (the standard gives no values): very loose and loose sand, loose sand-silt, medium-dense and
 * dense silt.
 */
export const API_SAND_TABLE = Object.freeze({
  'medium-dense sand-silt': Object.freeze({ beta: 0.29, fLim: 67e3, Nq: 12, qLim: 3e6 }),
  'medium-dense sand': Object.freeze({ beta: 0.37, fLim: 81e3, Nq: 20, qLim: 5e6 }),
  'dense sand-silt': Object.freeze({ beta: 0.37, fLim: 81e3, Nq: 20, qLim: 5e6 }),
  'dense sand': Object.freeze({ beta: 0.46, fLim: 96e3, Nq: 40, qLim: 10e6 }),
  'very dense sand-silt': Object.freeze({ beta: 0.46, fLim: 96e3, Nq: 40, qLim: 10e6 }),
  'very dense sand': Object.freeze({ beta: 0.56, fLim: 115e3, Nq: 50, qLim: 12e6 }),
});

// Densities for which the table gives no values, and the class substituted for them on request.
const LOOSE_CLASSES = ['very loose sand', 'loose sand', 'very loose sand-silt', 'loose sand-silt'];
const LOWEST_SAND_CLASS = 'medium-dense sand-silt';

/**
 * Density description from relative density: very loose 0–15 %, loose 15–35 %, medium dense 35–65 %,
 * dense 65–85 %, very dense 85–100 %.
 * Reference: API RP 2GEO (2011) Table 1 footnote (API RP 2A-WSD Table 6.4.3-1 notes).
 * @param {number} Dr relative density [-], 0..1
 * @returns {string} 'very loose' | 'loose' | 'medium-dense' | 'dense' | 'very dense'
 */
export function densityDescription(Dr) {
  if (!(Dr >= 0 && Dr <= 1)) throw new RangeError(`Dr: must be within 0..1 (got ${Dr})`);
  if (Dr < 0.15) return 'very loose';
  if (Dr < 0.35) return 'loose';
  if (Dr < 0.65) return 'medium-dense';
  if (Dr < 0.85) return 'dense';
  return 'very dense';
}

/**
 * API sand design parameters for a layer. Priority: explicit layer.beta/fLim/Nq/qLim, then layer.apiSand (a key
 * of API_SAND_TABLE), then the density description from layer.Dr with 'sand' or (for silt) 'sand-silt'.
 * Reference: API RP 2GEO (2011) Table 1. Validity: throws for densities the table does not cover (very loose and
 * loose sand, loose sand-silt), for which the standard asks for CPT-based methods or site-specific data.
 * With `allowLooseSand: true` such a layer is given the lowest tabulated class instead, 'medium-dense sand-silt'
 * (β = 0.29, fLim = 67 kPa, Nq = 12, qLim = 3 MPa), and the result carries `screening: true` with a `warning`:
 * the table does not support these values for loose material, the true resistance is probably lower, and the
 * figure serves only to let a screening calculation run to the end.
 * @param {object} layer layer record (soilprofile.js)
 * @param {{ allowLooseSand?: boolean }} [o] allowLooseSand: substitute the lowest tabulated class for a density
 *   the table does not cover (default false)
 * @returns {{ beta: number, fLim: number, Nq: number, qLim: number, category: string, screening?: boolean,
 *   warning?: string }} β [-], fLim [Pa], Nq [-], qLim [Pa], category label; screening and warning only when the
 *   lowest class was substituted
 */
export function apiSandParameters(layer, o = {}) {
  if (layer.beta !== undefined && layer.fLim !== undefined && layer.Nq !== undefined && layer.qLim !== undefined) {
    return { beta: layer.beta, fLim: layer.fLim, Nq: layer.Nq, qLim: layer.qLim, category: 'user' };
  }
  let key = layer.apiSand;
  if (key === undefined) {
    if (layer.Dr === undefined) throw new RangeError('layer: needs apiSand, Dr, or explicit beta/fLim/Nq/qLim for the API sand method');
    key = `${densityDescription(layer.Dr)} ${layer.type === 'silt' ? 'sand-silt' : 'sand'}`;
  }
  const row = API_SAND_TABLE[key];
  if (!row) {
    if (o.allowLooseSand && LOOSE_CLASSES.includes(key)) {
      return { ...API_SAND_TABLE[LOWEST_SAND_CLASS], category: `${LOWEST_SAND_CLASS} (substituted for ${key})`, screening: true,
        warning: `API RP 2GEO Table 1 has no design parameters for ${key}; the lowest tabulated class (${LOWEST_SAND_CLASS}) was used. `
          + 'Screening value only: the resistance of loose material is probably lower — use a CPT-based method or site-specific parameters for design.' };
    }
    throw new RangeError(`layer: API RP 2GEO Table 1 has no design parameters for '${key}'; supply beta, fLim, Nq and qLim`
      + (LOOSE_CLASSES.includes(key) ? ', or pass allowLooseSand: true for a screening value from the lowest tabulated class' : ''));
  }
  return { ...row, category: key };
}

/**
 * Adhesion factor for clay: α = 0.5 ψ^−0.5 for ψ <= 1, 0.5 ψ^−0.25 for ψ > 1, α <= 1, with ψ = su/σ′v.
 * Reference: API RP 2GEO (2011) §8.1 (API RP 2A-WSD §6.4.2). Validity: driven pipe piles in clay.
 * @param {number} psi strength ratio su/σ′v [-], > 0
 * @returns {number} α [-]
 */
export function clayAlpha(psi) {
  pos('psi', psi);
  return Math.min(1, psi <= 1 ? 0.5 / Math.sqrt(psi) : 0.5 * Math.pow(psi, -0.25));
}

/**
 * Unit skin friction f [Pa] and unit end bearing q [Pa] at a depth by the API main-text method:
 * clay f = α su, q = 9 su; sand f = min(β σ′v, fLim), q = min(Nq σ′v, qLim).
 * Silt is treated as clay when it has su, otherwise as sand-silt. Rock layers must supply fRock and qRock [Pa].
 * Above the scoured surface f = q = 0; σ′v includes the scour overburden reduction.
 * Reference: API RP 2GEO (2011) §8.1 (API RP 2A-WSD §6.4.2–6.4.3).
 * @param {object} profile soil profile
 * @param {number} z depth below original mudline [m]
 * @param {{ betaFactor?: number, allowLooseSand?: boolean }} [o] betaFactor [-] multiplier on β (1.25 for
 *   full-displacement piles; default 1); allowLooseSand: see apiSandParameters (default false)
 * @returns {{ f: number, q: number, type: string, sigmaV: number, alpha?: number, beta?: number,
 *   screening?: boolean, warning?: string }} f, q, sigmaV [Pa]; screening / warning when the lowest sand class was
 *   substituted for loose material
 */
export function unitResistanceAt(profile, z, o = {}) {
  const betaFactor = o.betaFactor ?? 1;
  pos('betaFactor', betaFactor);
  const pr = propertiesAt(profile, z);
  const sigmaV = verticalEffectiveStress(profile, z);
  if (z < scourBaseDepth(profile) - 1e-12) return { f: 0, q: 0, type: pr.type, sigmaV: 0 };
  if (pr.type === 'rock') {
    const l = pr.layer;
    if (l.fRock === undefined || l.qRock === undefined) throw new RangeError(`profile: rock layer ${pr.layerIndex} needs fRock and qRock`);
    return { f: l.fRock, q: l.qRock, type: 'rock', sigmaV };
  }
  if (pr.type === 'clay' || (pr.type === 'silt' && pr.su !== undefined)) {
    if (!(pr.su > 0)) return { f: 0, q: 0, type: 'clay', sigmaV, alpha: 1 };
    // at the soil surface σ′v = 0 → ψ = ∞ → α → 0 by the ψ > 1 branch
    const alpha = sigmaV > 0 ? clayAlpha(pr.su / sigmaV) : 0;
    return { f: alpha * pr.su, q: 9 * pr.su, type: 'clay', sigmaV, alpha };
  }
  const s = apiSandParameters(pr.layer, { allowLooseSand: o.allowLooseSand });
  const beta = s.beta * betaFactor;
  const out = { f: Math.min(beta * sigmaV, s.fLim), q: Math.min(s.Nq * sigmaV, s.qLim), type: 'sand', sigmaV, beta };
  if (s.screening) { out.screening = true; out.warning = s.warning; }
  return out;
}

const GX = [-Math.sqrt(0.6), 0, Math.sqrt(0.6)];
const GW = [5 / 9, 8 / 9, 5 / 9];

function slices(profile, L, dz) {
  const cuts = [0];
  const add = (z) => { if (z > cuts[cuts.length - 1] + 1e-12 && z < L - 1e-12) cuts.push(z); };
  const stops = profile.layers.map((l) => l.bottom);
  if (profile.scour) stops.push(profile.scour.general + profile.scour.local, profile.scour.general + profile.scour.reductionDepth);
  stops.sort((a, b) => a - b).forEach(add);
  cuts.push(L);
  const out = [];
  for (let i = 0; i < cuts.length - 1; i++) {
    const n = Math.max(1, Math.ceil((cuts[i + 1] - cuts[i]) / dz - 1e-9));
    const h = (cuts[i + 1] - cuts[i]) / n;
    for (let j = 0; j < n; j++) out.push([cuts[i] + j * h, cuts[i] + (j + 1) * h]);
  }
  return out;
}

/**
 * Ultimate axial capacity of a pipe pile at penetration L.
 *   Outer shaft friction Qfo = ∫ f π D dz; inner shaft friction Qfi = ∫ f π Di dz · internalFrictionFactor;
 *   plugged end bearing q π D²/4; unplugged (coring) end bearing q A_wall + Qfi; the smaller governs.
 *   Compression Qc = Qfo + Qp; tension Qt = tensionFactor · Qfo (pile and plug weight not included).
 * setupFactor multiplies shaft friction (ageing / set-up, or a reduction for driving disturbance); it is an
 * engineering input, 1 by default.
 * Reference: API RP 2GEO (2011) §8.1–8.2 (API RP 2A-WSD §6.4–6.5). Validity: static capacity of driven pipe
 * piles in siliceous soils; pile weight, cyclic degradation and group effects are not included. Loose and very
 * loose sand layers are outside the API table and throw; `allowLooseSand: true` lets the calculation run with the
 * lowest tabulated class for those layers and marks the result `screening: true` with `warnings` (not for design).
 * @param {{ profile: object, L: number, D: number, t: number, closedEnded?: boolean, plugMode?: string,
 *   internalFrictionFactor?: number, tensionFactor?: number, setupFactor?: number, betaFactor?: number,
 *   dz?: number }} o profile; L [m] penetration below original mudline; D [m] outer diameter; t [m] wall thickness;
 *   closedEnded [bool]; plugMode 'auto'|'plugged'|'unplugged' (default 'auto' = smaller end bearing);
 *   internalFrictionFactor [-] inner/outer unit friction ratio (default 1); tensionFactor [-] (default 1);
 *   setupFactor [-] (default 1); betaFactor [-] (default 1); dz [m] integration slice (default 0.25);
 *   allowLooseSand [bool] (default false)
 * @returns {object} { L [m], Qfo, Qfi, QpPlugged, QpUnplugged, Qp, Qc, Qt [N], plugged [bool], qTip [Pa],
 *   z [m] slice mid-depths, f [Pa] unit friction at mid-depths (incl. setupFactor), QfoCum [N] cumulative outer
 *   friction at slice bottoms, screening [bool] true when a loose layer was given the lowest tabulated class,
 *   warnings [string[]] }
 */
export function axialCapacity(o) {
  const { profile, L, D, t, closedEnded = false, plugMode = 'auto', internalFrictionFactor = 1, tensionFactor = 1, setupFactor = 1, betaFactor = 1, dz = 0.25,
    allowLooseSand = false } = o;
  const warnings = [];
  const at = (zz) => {
    const u = unitResistanceAt(profile, zz, { betaFactor, allowLooseSand });
    if (u.screening && !warnings.includes(u.warning)) warnings.push(u.warning);
    return u;
  };
  pos('L', L); pos('D', D); pos('t', t); pos('dz', dz);
  pos('internalFrictionFactor', internalFrictionFactor, true); pos('tensionFactor', tensionFactor, true); pos('setupFactor', setupFactor, true);
  if (t >= D / 2) throw new RangeError('t: must be less than D/2');
  if (L > profile.depth + 1e-9) throw new RangeError(`L: ${L} m exceeds the profile depth ${profile.depth} m`);
  if (!['auto', 'plugged', 'unplugged'].includes(plugMode)) throw new RangeError("plugMode: must be 'auto', 'plugged' or 'unplugged'");
  const Di = D - 2 * t;
  let integral = 0;
  const z = [], f = [], QfoCum = [];
  for (const [a, b] of slices(profile, L, dz)) {
    const h = b - a, mid = 0.5 * (a + b);
    let s = 0;
    for (let g = 0; g < 3; g++) s += GW[g] * at(mid + 0.5 * h * GX[g]).f;
    integral += 0.5 * h * s * setupFactor;
    z.push(mid); f.push(at(mid).f * setupFactor); QfoCum.push(integral * Math.PI * D);
  }
  const Qfo = integral * Math.PI * D;
  const Qfi = closedEnded ? 0 : integral * Math.PI * Di * internalFrictionFactor;
  const qTip = at(Math.min(L, profile.depth)).q;
  const QpPlugged = (qTip * Math.PI * D * D) / 4;
  const QpUnplugged = closedEnded ? QpPlugged : (qTip * Math.PI * (D * D - Di * Di)) / 4 + Qfi;
  let plugged;
  if (closedEnded || plugMode === 'plugged') plugged = true;
  else if (plugMode === 'unplugged') plugged = false;
  else plugged = QpPlugged <= QpUnplugged;
  const Qp = plugged ? QpPlugged : QpUnplugged;
  return { L, Qfo, Qfi, QpPlugged, QpUnplugged, Qp, Qc: Qfo + Qp, Qt: tensionFactor * Qfo, plugged, qTip, z, f, QfoCum,
    screening: warnings.length > 0, warnings };
}

/**
 * Capacity versus penetration curve.
 * Reference: see axialCapacity.
 * @param {object} o axialCapacity options without L, plus Lmin [m] (default D), Lmax [m] (default profile depth),
 *   step [m] (default 1)
 * @returns {{ L: number[], Qc: number[], Qt: number[], Qfo: number[], Qp: number[], plugged: boolean[],
 *   screening: boolean, warnings: string[] }} L [m]; capacities [N]; screening / warnings as axialCapacity (any
 *   penetration)
 */
export function capacityVsPenetration(o) {
  const Lmin = o.Lmin ?? o.D, Lmax = o.Lmax ?? o.profile.depth, step = o.step ?? 1;
  pos('step', step);
  if (!(Lmax > Lmin && Lmin > 0)) throw new RangeError('Lmax: must exceed Lmin > 0');
  const out = { L: [], Qc: [], Qt: [], Qfo: [], Qp: [], plugged: [], screening: false, warnings: [] };
  const n = Math.ceil((Lmax - Lmin) / step - 1e-9);
  for (let i = 0; i <= n; i++) {
    const L = Math.min(Lmin + i * step, Lmax);
    const c = axialCapacity({ ...o, L });
    out.L.push(L); out.Qc.push(c.Qc); out.Qt.push(c.Qt); out.Qfo.push(c.Qfo); out.Qp.push(c.Qp); out.plugged.push(c.plugged);
    if (c.screening) { out.screening = true; for (const w of c.warnings) if (!out.warnings.includes(w)) out.warnings.push(w); }
  }
  return out;
}

/**
 * Smallest penetration at which the design resistance reaches the load:
 *   design resistance = Q · resistanceFactor / safetyFactor, Q = Qc (compression) or Qt (tension).
 * Working-stress design uses safetyFactor (API RP 2A-WSD §6.3.4: 1.5 for extreme, 2.0 for operating conditions)
 * with resistanceFactor = 1; load-and-resistance-factor design uses resistanceFactor with factored loads and
 * safetyFactor = 1. Both are arguments; no default is implied by this function beyond 1.
 * The search steps down the profile and bisects in the first interval that satisfies the load, so a weaker
 * layer deeper down is not screened: review the full capacity curve.
 * Reference: API RP 2A-WSD §6.3.4; see axialCapacity. Validity: as axialCapacity.
 * @param {object} o axialCapacity options without L, plus load [N] (> 0), mode 'compression'|'tension',
 *   safetyFactor [-] (default 1), resistanceFactor [-] (default 1), Lmin [m], Lmax [m], step [m] (default 1),
 *   tol [m] (default 0.01)
 * @returns {{ L: number, capacity: number, designResistance: number, found: boolean, screening: boolean,
 *   warnings: string[] }} L [m] required penetration (Lmax if not found); capacity and designResistance [N] at L;
 *   screening / warnings as axialCapacity at that penetration
 */
export function requiredPenetration(o) {
  const { load, mode = 'compression', safetyFactor = 1, resistanceFactor = 1, step = 1, tol = 0.01 } = o;
  pos('load', load); pos('safetyFactor', safetyFactor); pos('resistanceFactor', resistanceFactor); pos('step', step); pos('tol', tol);
  if (mode !== 'compression' && mode !== 'tension') throw new RangeError("mode: must be 'compression' or 'tension'");
  const Lmin = o.Lmin ?? o.D, Lmax = o.Lmax ?? o.profile.depth;
  const R = (L) => { const c = axialCapacity({ ...o, L }); return ((mode === 'compression' ? c.Qc : c.Qt) * resistanceFactor) / safetyFactor; };
  const done = (L, found) => {
    const d = R(L); const c = axialCapacity({ ...o, L });
    return { L, capacity: (d * safetyFactor) / resistanceFactor, designResistance: d, found, screening: c.screening, warnings: c.warnings };
  };
  if (R(Lmin) >= load) return done(Lmin, true);
  let a = Lmin;
  while (a < Lmax - 1e-12) {
    let b = Math.min(a + step, Lmax);
    if (R(b) >= load) {
      while (b - a > tol) { const m = 0.5 * (a + b); if (R(m) >= load) b = m; else a = m; }
      return done(b, true);
    }
    a = b;
  }
  return done(Lmax, false);
}

/**
 * Unified CPT-based unit shaft friction for driven open-ended steel pipe piles in sand:
 *   f(z) = u · qc · (σ′v/pa)^a · Ar^b · [max((L − z)/D, v)]^(−c) · (tan δcv)^d · [min((L − z)/(D v), 1)]^e
 * with Ar = 1 − (Di/D)² and pa = 100 kPa. The coefficient set selects the method (see CPT_SHAFT_COEFFICIENTS).
 * Reference: API RP 2GEO (2011) Annex C (API RP 2A-WSD 21st ed. Supplement 3, commentary C6.4.3), unified form of the
 * simplified ICP-05, offshore UWA-05 and Fugro-05 methods. Validity: siliceous sands; qc averaged as the method
 * requires. End bearing by these methods is not provided here except cptEndBearingUWA.
 * @param {{ qc: number, sigmaV: number, z: number, L: number, D: number, t: number, deltaCvDeg: number,
 *   coeff: { a: number, b: number, c: number, d: number, e: number, u: number, v: number | string },
 *   pa?: number }} o qc [Pa] cone resistance at depth z; sigmaV [Pa] vertical effective stress; z [m] depth;
 *   L [m] pile penetration; D [m]; t [m]; deltaCvDeg [deg] constant-volume interface friction angle; coeff the
 *   coefficient set, v either a number or '4sqrtAr' / '2sqrtAr'; pa [Pa] (default 100e3)
 * @returns {number} unit shaft friction f [Pa]
 */
export function cptShaftFriction({ qc, sigmaV, z, L, D, t, deltaCvDeg, coeff, pa = 100e3 }) {
  pos('qc', qc, true); pos('sigmaV', sigmaV, true); pos('z', z, true); pos('L', L); pos('D', D); pos('t', t); pos('deltaCvDeg', deltaCvDeg); pos('pa', pa);
  if (z > L) throw new RangeError('z: must not exceed L');
  if (!coeff) throw new RangeError('coeff: coefficient set is required');
  const Di = D - 2 * t;
  const Ar = 1 - (Di / D) ** 2;
  let v = coeff.v;
  if (v === '4sqrtAr') v = 4 * Math.sqrt(Ar);
  else if (v === '2sqrtAr') v = 2 * Math.sqrt(Ar);
  pos('coeff.v', v);
  const h = (L - z) / D;
  return coeff.u * qc * Math.pow(sigmaV / pa, coeff.a) * Math.pow(Ar, coeff.b) * Math.pow(Math.max(h, v), -coeff.c)
    * Math.pow(Math.tan((deltaCvDeg * Math.PI) / 180), coeff.d) * Math.pow(Math.min(h / v, 1), coeff.e);
}

/**
 * Coefficient sets for cptShaftFriction (open-ended steel pipe piles): simplified ICP-05 and offshore UWA-05,
 * compression and tension.
 * Reference: API RP 2GEO (2011) Annex C, table of unit skin friction parameter values. The Fugro-05 and NGI-05
 * methods are not included (see docs/engine/geotech-piles.md).
 */
export const CPT_SHAFT_COEFFICIENTS = Object.freeze({
  'ICP-05 compression': Object.freeze({ a: 0.1, b: 0.2, c: 0.4, d: 1, e: 0, u: 0.023, v: '4sqrtAr' }),
  'ICP-05 tension': Object.freeze({ a: 0.1, b: 0.2, c: 0.4, d: 1, e: 0, u: 0.016, v: '4sqrtAr' }),
  'UWA-05 compression': Object.freeze({ a: 0, b: 0.3, c: 0.5, d: 1, e: 0, u: 0.030, v: 2 }),
  'UWA-05 tension': Object.freeze({ a: 0, b: 0.3, c: 0.5, d: 1, e: 0, u: 0.022, v: 2 }),
});

/**
 * Offshore UWA-05 unit end bearing on the gross pile area: q = qc,av (0.15 + 0.45 Ar), Ar = 1 − (Di/D)².
 * Reference: API RP 2GEO (2011) Annex C (offshore UWA-05); Lehane et al. (2005). Validity: open-ended driven
 * pipe piles in siliceous sand; qc,av is the average cone resistance over ±1.5 D about the tip.
 * @param {{ qcAvg: number, D: number, t: number }} o qcAvg [Pa]; D [m]; t [m]
 * @returns {{ q: number, Qp: number, Ar: number }} q [Pa] on gross area; Qp [N]; Ar [-]
 */
export function cptEndBearingUWA({ qcAvg, D, t }) {
  pos('qcAvg', qcAvg, true); pos('D', D); pos('t', t);
  const Ar = 1 - ((D - 2 * t) / D) ** 2;
  const q = qcAvg * (0.15 + 0.45 * Ar);
  return { q, Qp: (q * Math.PI * D * D) / 4, Ar };
}

/**
 * Block-failure capacity of a pile group in clay: shear on the perimeter of the enclosing block plus bearing
 * on its base, Qblock = 2 (B + W) L su,shaft + Nc su,base B W, compared with the sum of single-pile capacities.
 * Reference: Terzaghi & Peck (1967) block-failure mechanism; Nc = 9 for a deep square base (Skempton 1951).
 * Validity: closely spaced groups in clay; undrained; α = 1 on the block perimeter (soil-to-soil shear).
 * @param {{ B: number, W: number, L: number, suShaft: number, suBase: number, Nc?: number, nPiles: number,
 *   Qsingle: number }} o B, W [m] plan dimensions of the block (outside of outer piles); L [m] embedded length;
 *   suShaft [Pa] average su over L; suBase [Pa] su at the base; Nc [-] (default 9); nPiles [-]; Qsingle [N]
 * @returns {{ Qblock: number, Qsum: number, Qgroup: number, efficiency: number, governs: string }} capacities
 *   [N]; efficiency = Qgroup / Qsum [-]; governs 'block' | 'individual'
 */
export function groupBlockCapacityClay({ B, W, L, suShaft, suBase, Nc = 9, nPiles, Qsingle }) {
  pos('B', B); pos('W', W); pos('L', L); pos('suShaft', suShaft, true); pos('suBase', suBase, true); pos('Nc', Nc); pos('nPiles', nPiles); pos('Qsingle', Qsingle);
  const Qblock = 2 * (B + W) * L * suShaft + Nc * suBase * B * W;
  const Qsum = nPiles * Qsingle;
  const Qgroup = Math.min(Qblock, Qsum);
  return { Qblock, Qsum, Qgroup, efficiency: Qgroup / Qsum, governs: Qblock < Qsum ? 'block' : 'individual' };
}

