/**
 * Reinforced concrete sections and durability.
 *
 *  - Material laws of EN 1992-1-1:2004: strength classes (Table 3.1), parabola–rectangle concrete in compression
 *    (§3.1.7), bilinear reinforcing steel (§3.2.7).
 *  - Fibre (strip) integration of rectangular, T and circular sections in uniaxial bending with axial force:
 *    moment–curvature, ultimate moment, N–M interaction diagram (strain pivots of §6.1), cracking moment, elastic
 *    cracked-section properties.
 *  - Shear resistance without and with links (§6.2.2, §6.2.3), crack width (§7.3.4), deflection by interpolation
 *    between the uncracked and fully cracked states (§7.4.3).
 *  - Durability: chloride ingress by Fick's second law with an ageing diffusion coefficient, time to corrosion
 *    initiation, carbonation depth, corrosion propagation (Faraday's law), bar section loss, residual moment
 *    capacity with time; thickness loss of steel piles by exposure zone.
 *
 * Conventions: SI units (m, N, Pa, s). Section depth coordinate y is measured from the top fibre; strains and
 * stresses in the section analysis are compression-positive; axial force N is compression-positive; bending
 * moment M is positive when it compresses the top fibre and is taken about the centroid of the gross concrete
 * section. Concrete tensile strength is neglected in all strength calculations.
 */
import { brent } from '../math/roots.js';
import { erfc, erfinv } from '../math/special.js';

/** One Julian year in seconds, for converting durability times and rates. */
export const YEAR = 365.25 * 86400;

function pos(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v > 0)) throw new RangeError(`${name} must be a positive finite number`);
  return v;
}
function nonNeg(name, v) {
  if (!(typeof v === 'number' && Number.isFinite(v) && v >= 0)) throw new RangeError(`${name} must be a non-negative finite number`);
  return v;
}

/* ------------------------------------------------------------------ materials */

/**
 * Strength and deformation characteristics of normal-weight concrete, EN 1992-1-1 Table 3.1:
 * fcm = fck + 8 MPa; fctm = 0.30 fck^(2/3) (≤ C50/60) or 2.12 ln(1 + fcm/10) (> C50/60); fctk,0.05 = 0.7 fctm;
 * Ecm = 22 (fcm/10)^0.3 GPa; parabola–rectangle parameters εc2 = 2.0 ‰, εcu2 = 3.5 ‰, n = 2 (≤ C50/60) or
 * εc2 = 2.0 + 0.085 (fck − 50)^0.53 ‰, εcu2 = 2.6 + 35 [(90 − fck)/100]⁴ ‰, n = 1.4 + 23.4 [(90 − fck)/100]⁴.
 * @param {object} o
 * @param {number} o.fck characteristic cylinder strength [Pa], 12–90 MPa
 * @returns {{fck:number, fcm:number, fctm:number, fctk005:number, Ecm:number, epsC2:number, epsCu2:number, n:number}} strengths and modulus [Pa], strains [-]
 */
export function concreteProperties({ fck }) {
  pos('fck', fck);
  const f = fck / 1e6;
  if (f < 12 || f > 90) throw new RangeError('fck must lie between 12 and 90 MPa');
  const fcm = f + 8;
  const fctm = f <= 50 ? 0.3 * f ** (2 / 3) : 2.12 * Math.log(1 + fcm / 10);
  const hi = f > 50;
  const q = ((90 - f) / 100) ** 4;
  return {
    fck, fcm: fcm * 1e6, fctm: fctm * 1e6, fctk005: 0.7 * fctm * 1e6, Ecm: 22e9 * (fcm / 10) ** 0.3,
    epsC2: hi ? (2.0 + 0.085 * (f - 50) ** 0.53) * 1e-3 : 2.0e-3,
    epsCu2: hi ? (2.6 + 35 * q) * 1e-3 : 3.5e-3,
    n: hi ? 1.4 + 23.4 * q : 2,
  };
}

/**
 * Parabola–rectangle stress–strain relation for concrete in compression, EN 1992-1-1 Eq. (3.17), (3.18):
 * σc = fcd [1 − (1 − εc/εc2)^n] for 0 ≤ εc ≤ εc2, σc = fcd for εc2 ≤ εc ≤ εcu2; zero in tension.
 * @param {object} o
 * @param {number} o.strain compressive strain [-] (positive in compression)
 * @param {number} o.fcd design compressive strength [Pa]
 * @param {number} [o.epsC2=0.002] strain at reaching the maximum strength [-]
 * @param {number} [o.epsCu2=0.0035] ultimate strain [-] (strains beyond it are treated as the plateau; limit them in the analysis)
 * @param {number} [o.n=2] exponent [-]
 * @returns {number} compressive stress [Pa]
 */
export function concreteStress({ strain, fcd, epsC2 = 0.002, epsCu2 = 0.0035, n = 2 }) {
  pos('fcd', fcd); pos('epsC2', epsC2); pos('epsCu2', epsCu2); pos('n', n);
  if (!(strain > 0)) return 0;
  return strain >= epsC2 ? fcd : fcd * (1 - (1 - strain / epsC2) ** n);
}

/**
 * Design stress–strain relation of reinforcing steel, EN 1992-1-1 §3.2.7(2): elastic up to fyd, then either a
 * horizontal branch (k = 1, no strain limit needed) or an inclined branch reaching k·fyd at εuk and limited to εud.
 * Symmetric in tension and compression.
 * @param {object} o
 * @param {number} o.strain strain [-]
 * @param {number} o.fyd design yield strength [Pa]
 * @param {number} [o.Es=2.0e11] modulus of elasticity [Pa]
 * @param {number} [o.k=1] ratio (ft/fy)k [-]
 * @param {number} [o.epsUk=0.05] characteristic strain at maximum force [-] (used when k > 1)
 * @returns {number} stress [Pa], same sign as the strain
 */
export function steelStress({ strain, fyd, Es = 2.0e11, k = 1, epsUk = 0.05 }) {
  pos('fyd', fyd); pos('Es', Es);
  const e = Math.abs(strain); const ey = fyd / Es;
  let s;
  if (e <= ey) s = Es * e;
  else if (k > 1) s = fyd * (1 + (k - 1) * Math.min(1, (e - ey) / (epsUk - ey)));
  else s = fyd;
  return strain < 0 ? -s : s;
}

function resolveMaterial(m) {
  if (!m) throw new RangeError('material is required');
  const cp = concreteProperties({ fck: m.fck });
  const gammaC = m.gammaC === undefined ? 1.5 : pos('material.gammaC', m.gammaC);
  const gammaS = m.gammaS === undefined ? 1.15 : pos('material.gammaS', m.gammaS);
  const alphaCC = m.alphaCC === undefined ? 1 : pos('material.alphaCC', m.alphaCC);
  const fcd = m.fcd === undefined ? alphaCC * m.fck / gammaC : pos('material.fcd', m.fcd);
  const fyd = m.fyd === undefined ? pos('material.fyk', m.fyk) / gammaS : pos('material.fyd', m.fyd);
  const Es = m.Es === undefined ? 2.0e11 : pos('material.Es', m.Es);
  return {
    fcd, fyd, Es, k: m.k === undefined ? 1 : m.k, epsUk: m.epsUk === undefined ? 0.05 : m.epsUk,
    epsUd: m.epsUd === undefined ? Infinity : pos('material.epsUd', m.epsUd),
    epsC2: m.epsC2 === undefined ? cp.epsC2 : m.epsC2, epsCu2: m.epsCu2 === undefined ? cp.epsCu2 : m.epsCu2, n: m.n === undefined ? cp.n : m.n,
  };
}

/* ------------------------------------------------------------------ section geometry */

/**
 * Reinforced concrete cross-section for uniaxial bending, discretised into horizontal strips.
 *   shape 'rect': width b, depth h;  'T': flange bf × hf at the top over a web bw, total depth h;  'circle': diameter D.
 * Bars are lumped layers { depth, area, diameter? } with the depth measured from the top fibre; the bar diameter is only
 * needed for the corrosion functions. Concrete displaced by the bars is not deducted.
 * @param {object} o
 * @param {'rect'|'T'|'circle'} o.shape
 * @param {number} [o.b] width [m] (rect)
 * @param {number} [o.h] depth [m] (rect, T)
 * @param {number} [o.bf] flange width [m] (T)
 * @param {number} [o.hf] flange thickness [m] (T)
 * @param {number} [o.bw] web width [m] (T)
 * @param {number} [o.D] diameter [m] (circle)
 * @param {Array<{depth:number, area:number, diameter?:number}>} [o.bars=[]] reinforcement layers [m], [m²], [m]
 * @param {number} [o.nStrips=400] number of concrete strips
 * @returns {{shape:string, h:number, area:number, centroid:number, inertia:number, bars:object[], y:Float64Array, dA:Float64Array,
 *   dy:Float64Array, widthAt:function(number):number, dMax:number, As:number}} depth [m], gross area [m²], centroid depth from the top [m],
 *   gross second moment about the centroid [m⁴], strip mid-depths [m], areas [m²] and heights [m], width function [m], depth of the lowest
 *   bar layer [m], total steel area [m²]
 */
export function rcSection({ shape, b, h, bf, hf, bw, D, bars = [], nStrips = 400 }) {
  let depth; let widthAt;
  if (shape === 'rect') { pos('b', b); depth = pos('h', h); widthAt = () => b; }
  else if (shape === 'T') {
    pos('bf', bf); pos('hf', hf); pos('bw', bw); depth = pos('h', h);
    if (hf >= h) throw new RangeError('hf must be smaller than h');
    if (bw > bf) throw new RangeError('bw must not exceed bf');
    widthAt = (y) => (y <= hf ? bf : bw);
  } else if (shape === 'circle') {
    depth = pos('D', D);
    widthAt = (y) => 2 * Math.sqrt(Math.max(0, y * (D - y)));
  } else throw new RangeError('shape must be "rect", "T" or "circle"');
  if (!Number.isInteger(nStrips) || nStrips < 10) throw new RangeError('nStrips must be an integer >= 10');
  const ys = new Float64Array(nStrips); const dA = new Float64Array(nStrips); const dys = new Float64Array(nStrips);
  // a T-section is divided so that the flange/web junction falls on a strip boundary
  const bounds = shape === 'T' ? [0, hf, depth] : [0, depth];
  let k = 0;
  for (let s = 0; s + 1 < bounds.length; s++) {
    const span = bounds[s + 1] - bounds[s];
    const n = s + 2 === bounds.length ? nStrips - k : Math.max(1, Math.round(nStrips * span / depth));
    for (let i = 0; i < n; i++) {
      const y0 = bounds[s] + span * i / n; const y1 = bounds[s] + span * (i + 1) / n;
      ys[k] = 0.5 * (y0 + y1); dys[k] = y1 - y0;
      // circle: exact strip area from the segment formula so that the total area is exact
      if (shape === 'circle') {
        const R = D / 2; const seg = (y) => { const u = Math.min(1, Math.max(-1, (y - R) / R)); return R * R * (Math.asin(u) + u * Math.sqrt(1 - u * u)); };
        dA[k] = seg(y1) - seg(y0);
      } else dA[k] = widthAt(ys[k]) * (y1 - y0);
      k++;
    }
  }
  let area = 0; let s1 = 0;
  for (let i = 0; i < nStrips; i++) { area += dA[i]; s1 += dA[i] * ys[i]; }
  const centroid = s1 / area;
  let inertia = 0;
  for (let i = 0; i < nStrips; i++) inertia += dA[i] * ((ys[i] - centroid) ** 2 + dys[i] * dys[i] / 12);
  if (shape === 'circle') inertia = Math.PI * D ** 4 / 64;
  let dMax = 0; let As = 0;
  const bl = bars.map((br, i) => {
    if (!(br.depth >= 0 && br.depth <= depth)) throw new RangeError(`bars[${i}].depth must lie within the section`);
    nonNeg(`bars[${i}].area`, br.area);
    dMax = Math.max(dMax, br.depth); As += br.area;
    return { depth: br.depth, area: br.area, diameter: br.diameter };
  });
  return { shape, h: depth, area, centroid, inertia, bars: bl, y: ys, dA, dy: dys, widthAt, dMax: bl.length ? dMax : depth, As };
}

/**
 * Bar layers of a circular section with n equal bars evenly spaced on a circle (one bar at the bottom).
 * @param {object} o
 * @param {number} o.D section diameter [m]
 * @param {number} o.cover distance from the concrete surface to the bar centre [m]
 * @param {number} o.nBars number of bars (>= 4)
 * @param {number} o.barDiameter bar diameter [m]
 * @returns {Array<{depth:number, area:number, diameter:number}>} one entry per bar, depth from the top [m]
 */
export function circularBars({ D, cover, nBars, barDiameter }) {
  pos('D', D); pos('cover', cover); pos('barDiameter', barDiameter);
  if (!Number.isInteger(nBars) || nBars < 4) throw new RangeError('nBars must be an integer >= 4');
  const r = D / 2 - cover;
  if (!(r > 0)) throw new RangeError('cover must be smaller than D/2');
  const a = Math.PI * barDiameter * barDiameter / 4;
  const out = [];
  for (let i = 0; i < nBars; i++) out.push({ depth: D / 2 + r * Math.cos(2 * Math.PI * i / nBars), area: a, diameter: barDiameter });
  return out;
}

/* Section forces for the strain plane ε(y) = epsTop − kappa·y. Returns [N, M, extreme steel tensile strain (≥ 0)]. */
function integrateSection(sec, mat, epsTop, kappa, out) {
  let N = 0; let M = 0;
  const { y, dA } = sec; const yc = sec.centroid;
  const { fcd, epsC2, n } = mat;
  for (let i = 0; i < y.length; i++) {
    const e = epsTop - kappa * y[i];
    if (e <= 0) continue;
    const s = e >= epsC2 ? fcd : fcd * (1 - (1 - e / epsC2) ** n);
    const f = s * dA[i];
    N += f; M += f * (yc - y[i]);
  }
  let et = 0;
  for (const b of sec.bars) {
    const e = epsTop - kappa * b.depth;
    const f = steelStress({ strain: e, fyd: mat.fyd, Es: mat.Es, k: mat.k, epsUk: mat.epsUk }) * b.area;
    N += f; M += f * (yc - b.depth);
    if (-e > et) et = -e;
  }
  out[0] = N; out[1] = M; out[2] = et;
  return out;
}

/* ------------------------------------------------------------------ ultimate limit state: pivots */

/*
 * Ultimate strain distributions (EN 1992-1-1 Fig. 6.1), parametrised by u ∈ [0, 3] from pure compression to pure tension:
 *  u ∈ [0, 1]  rotation about point C at depth (1 − εc2/εcu2) h where the strain stays εc2: uniform εc2 → (εcu2 at the top, 0 at the bottom)
 *  u ∈ [1, 2]  top strain εcu2, strain at the lowest bar from its value at u = 1 down to −εT
 *  u ∈ [2, 3]  strain at the lowest bar −εT, top strain from εcu2 down to −εT (uniform tension)
 * εT is εud of the steel; with the horizontal top branch (no limit) a nominal 2 % is used, which only fixes the
 * parametrisation of the tension-controlled range (all tension steel has yielded there).
 */
function pivotState(sec, mat, u) {
  const h = sec.h; const d = sec.dMax;
  const eT = Number.isFinite(mat.epsUd) ? mat.epsUd : 0.02;
  if (u <= 1) {
    const yC = (1 - mat.epsC2 / mat.epsCu2) * h;
    const top = mat.epsC2 + (mat.epsCu2 - mat.epsC2) * u;
    return [top, (top - mat.epsC2) / yC];
  }
  if (u <= 2) {
    const ed0 = mat.epsCu2 * (1 - d / h);
    const ed = ed0 + (-eT - ed0) * (u - 1);
    return [mat.epsCu2, (mat.epsCu2 - ed) / d];
  }
  const top = mat.epsCu2 + (-eT - mat.epsCu2) * (u - 2);
  return [top, (top + eT) / d];
}

/**
 * N–M interaction diagram of a reinforced concrete section at the ultimate limit state: the strain distribution is
 * rotated through the pivots of EN 1992-1-1 §6.1 (Fig. 6.1) from uniform compression εc2 to uniform tension and the
 * stresses of the parabola–rectangle and bilinear steel laws are integrated over the section.
 * @param {object} o
 * @param {object} o.section result of rcSection
 * @param {object} o.material { fck, fyk [Pa], gammaC = 1.5, gammaS = 1.15, alphaCC = 1 (recommended value; national annexes may give 0.85),
 *   Es = 200 GPa, k = 1, epsUk = 0.05, epsUd = Infinity; fcd, fyd may be given directly }
 * @param {number} [o.nPoints=91] number of points along the diagram (multiple of 3 plus 1 gives equal sampling of the three ranges)
 * @returns {{N:Float64Array, M:Float64Array, epsTop:Float64Array, kappa:Float64Array, NRdMax:number, NRdMin:number, balanced:{N:number, M:number}}}
 *   axial resistance (compression +) [N], moment resistance for bending that compresses the top [N·m], strain state of each point
 *   [-], [1/m]; squash load and tension resistance [N]; point of maximum moment
 */
export function interactionDiagram({ section, material, nPoints = 91 }) {
  const mat = resolveMaterial(material);
  if (!Number.isInteger(nPoints) || nPoints < 4) throw new RangeError('nPoints must be an integer >= 4');
  const N = new Float64Array(nPoints); const M = new Float64Array(nPoints); const eT = new Float64Array(nPoints); const ka = new Float64Array(nPoints);
  const r = [0, 0, 0];
  let best = 0;
  for (let i = 0; i < nPoints; i++) {
    const [top, kappa] = pivotState(section, mat, 3 * i / (nPoints - 1));
    integrateSection(section, mat, top, kappa, r);
    N[i] = r[0]; M[i] = r[1]; eT[i] = top; ka[i] = kappa;
    if (M[i] > M[best]) best = i;
  }
  return { N, M, epsTop: eT, kappa: ka, NRdMax: N[0], NRdMin: N[nPoints - 1], balanced: { N: N[best], M: M[best] } };
}

/**
 * Ultimate moment of resistance of a reinforced concrete section for a given axial force, by solving the ultimate
 * strain distribution (pivots of EN 1992-1-1 §6.1) for axial equilibrium.
 * @param {object} o
 * @param {object} o.section result of rcSection
 * @param {object} o.material see interactionDiagram
 * @param {number} [o.N=0] design axial force [N], compression positive
 * @returns {{MRd:number, x:number, kappa:number, epsTop:number, epsSteel:number, steelYielded:boolean, xOverD:number}} moment resistance about the
 *   gross centroid [N·m], neutral axis depth [m] (Infinity when the whole section is in compression), curvature [1/m], top strain [-], largest
 *   steel tensile strain [-], whether the extreme tension steel has yielded, x/d [-]
 */
export function momentCapacity({ section, material, N = 0 }) {
  const mat = resolveMaterial(material);
  if (!Number.isFinite(N)) throw new RangeError('N must be finite');
  const r = [0, 0, 0];
  const g = (u) => { const [top, kappa] = pivotState(section, mat, u); return integrateSection(section, mat, top, kappa, r)[0] - N; };
  if (g(0) < 0) throw new RangeError('N exceeds the compressive resistance of the section');
  if (g(3) > 0) throw new RangeError('N exceeds the tensile resistance of the section');
  const u = brent(g, 0, 3, { tol: 1e-13 });
  const [top, kappa] = pivotState(section, mat, u);
  integrateSection(section, mat, top, kappa, r);
  const x = kappa > 0 ? top / kappa : Infinity;
  return { MRd: r[1], x: x > 0 ? x : 0, kappa, epsTop: top, epsSteel: r[2], steelYielded: r[2] >= mat.fyd / mat.Es, xOverD: x / section.dMax };
}

/**
 * Moment–curvature relation of a reinforced concrete section under constant axial force. For each curvature the
 * top strain is solved for axial equilibrium with the parabola–rectangle concrete law (no tension) and bilinear steel;
 * the curve ends when the top strain reaches εcu2 (or the steel reaches εud). Use gammaC = gammaS = 1 with mean strengths
 * for a best-estimate response, design values for the design curve.
 * @param {object} o
 * @param {object} o.section result of rcSection
 * @param {object} o.material see interactionDiagram
 * @param {number} [o.N=0] axial force [N], compression positive
 * @param {number} [o.nPoints=60] number of curvature steps up to the ultimate curvature
 * @returns {{kappa:Float64Array, M:Float64Array, epsTop:Float64Array, epsSteel:Float64Array, x:Float64Array, Mu:number, kappaU:number,
 *   My:number, kappaY:number, ductility:number}} curvature [1/m], moment [N·m], top strain [-], extreme steel tensile strain [-], neutral
 *   axis depth [m]; ultimate moment and curvature; moment and curvature at first yield of the extreme tension steel (NaN if it does
 *   not yield); curvature ductility κu/κy [-]
 */
export function momentCurvature({ section, material, N = 0, nPoints = 60 }) {
  const mat = resolveMaterial(material);
  if (!Number.isInteger(nPoints) || nPoints < 2) throw new RangeError('nPoints must be an integer >= 2');
  const ult = momentCapacity({ section, material, N });
  if (!(ult.kappa > 0)) throw new RangeError('momentCurvature: the section is fully in compression or tension at this axial force');
  const r = [0, 0, 0];
  const solve = (kappa) => {
    const f = (top) => integrateSection(section, mat, top, kappa, r)[0] - N;
    // bracket: fully tensile strain plane up to the crushing strain
    const lo = -0.1; const hi = mat.epsCu2 * 1.0000001;
    if (f(hi) < 0) return null;
    const top = brent(f, lo, hi, { tol: 1e-14 });
    integrateSection(section, mat, top, kappa, r);
    return { top, M: r[1], es: r[2] };
  };
  const ks = []; const Ms = []; const tops = []; const ess = []; const xs = [];
  const m0 = solve(0);
  ks.push(0); Ms.push(m0 ? m0.M : 0); tops.push(m0 ? m0.top : 0); ess.push(m0 ? m0.es : 0); xs.push(Infinity);
  const ey = mat.fyd / mat.Es;
  let kappaY = NaN; let My = NaN;
  for (let i = 1; i <= nPoints; i++) {
    const kappa = ult.kappa * i / nPoints;
    const s = i === nPoints ? { top: ult.epsTop, M: ult.MRd, es: ult.epsSteel } : solve(kappa);
    if (!s) break;
    if (Number.isNaN(kappaY) && s.es >= ey && section.bars.length) {
      // first yield between the previous point and this one
      const k0 = ks[ks.length - 1];
      kappaY = brent((k) => { const q = solve(k); return (q ? q.es : 1) - ey; }, k0 === 0 ? 1e-9 * kappa : k0, kappa, { tol: 1e-12 });
      My = solve(kappaY).M;
    }
    ks.push(kappa); Ms.push(s.M); tops.push(s.top); ess.push(s.es); xs.push(s.top > 0 ? s.top / kappa : 0);
  }
  return {
    kappa: Float64Array.from(ks), M: Float64Array.from(Ms), epsTop: Float64Array.from(tops), epsSteel: Float64Array.from(ess), x: Float64Array.from(xs),
    Mu: ult.MRd, kappaU: ult.kappa, My, kappaY, ductility: ult.kappa / kappaY,
  };
}

/* ------------------------------------------------------------------ serviceability: elastic sections */

/**
 * Cracking moment of the uncracked transformed section: Mcr = (fctm + N/Au) Iu/(h − yu), with the reinforcement
 * transformed by αe = Es/Ec (EN 1992-1-1 §7.1(2), §7.4.3).
 * @param {object} o
 * @param {object} o.section result of rcSection
 * @param {number} o.fctm mean tensile strength [Pa] (or fctm,fl, fct,eff as appropriate)
 * @param {number} o.Ec concrete modulus [Pa] (Ecm, or Ec,eff = Ecm/(1 + φ) for sustained load)
 * @param {number} [o.Es=2.0e11] steel modulus [Pa]
 * @param {number} [o.N=0] axial force [N], compression positive
 * @returns {{Mcr:number, area:number, centroid:number, inertia:number, alphaE:number}} cracking moment [N·m], transformed area [m²], centroid
 *   depth [m], second moment [m⁴] (concrete units), modular ratio [-]
 */
export function crackingMoment({ section, fctm, Ec, Es = 2.0e11, N = 0 }) {
  pos('fctm', fctm); pos('Ec', Ec); pos('Es', Es);
  const ae = Es / Ec;
  let A = section.area; let S = section.area * section.centroid;
  for (const b of section.bars) { A += (ae - 1) * b.area; S += (ae - 1) * b.area * b.depth; }
  const yu = S / A;
  let I = section.inertia + section.area * (section.centroid - yu) ** 2;
  for (const b of section.bars) I += (ae - 1) * b.area * (b.depth - yu) ** 2;
  return { Mcr: Math.max(0, (fctm + N / A) * I / (section.h - yu)), area: A, centroid: yu, inertia: I, alphaE: ae };
}

/**
 * Fully cracked elastic section in pure bending (linear concrete in compression, no tension, reinforcement transformed
 * by αe; bars inside the compression zone by αe − 1): neutral axis from the first-moment condition, second moment about
 * it, and the stresses for a given moment. For a singly reinforced rectangular section this is the classical
 * x = d [√((αρ)² + 2αρ) − αρ], Icr = b x³/3 + α As (d − x)².
 * @param {object} o
 * @param {object} o.section result of rcSection
 * @param {number} o.Ec concrete modulus [Pa]
 * @param {number} [o.Es=2.0e11] steel modulus [Pa]
 * @param {number} [o.M] bending moment [N·m] for the stresses
 * @returns {{x:number, inertia:number, alphaE:number, sigmaC:(number|undefined), sigmaS:(number|undefined)}} neutral axis depth [m], cracked
 *   second moment in concrete units [m⁴], modular ratio, extreme concrete compressive stress and tensile stress in the lowest bar layer [Pa]
 */
export function crackedSection({ section, Ec, Es = 2.0e11, M }) {
  pos('Ec', Ec); pos('Es', Es);
  if (!section.bars.length) throw new RangeError('section has no reinforcement');
  const ae = Es / Ec;
  const { y, dA } = section;
  const props = (x) => {
    let S = 0; let I = 0;
    for (let i = 0; i < y.length; i++) {
      // partial strip at the neutral axis
      const dy = section.dy[i];
      const top = y[i] - 0.5 * dy;
      if (top >= x) break;
      const frac = Math.min(1, (x - top) / dy);
      const ym = top + 0.5 * frac * dy;
      S += dA[i] * frac * (x - ym); I += dA[i] * frac * ((x - ym) ** 2 + (frac * dy) ** 2 / 12);
    }
    for (const b of section.bars) {
      const f = b.depth < x ? ae - 1 : ae;
      S += f * b.area * (x - b.depth); I += f * b.area * (x - b.depth) ** 2;
    }
    return [S, I];
  };
  const x = brent((v) => props(v)[0], 1e-9 * section.h, section.dMax, { tol: 1e-14 });
  const I = props(x)[1];
  const out = { x, inertia: I, alphaE: ae };
  if (M !== undefined) { out.sigmaC = M * x / I; out.sigmaS = ae * M * (section.dMax - x) / I; }
  return out;
}

/**
 * Modulus for sustained loading: Ec,eff = Ecm/(1 + φ(∞, t0)) (EN 1992-1-1 Eq. 7.20).
 * @param {object} o
 * @param {number} o.Ecm secant modulus [Pa]
 * @param {number} o.phi creep coefficient [-]
 * @returns {number} Ec,eff [Pa]
 */
export function effectiveModulus({ Ecm, phi }) {
  pos('Ecm', Ecm); nonNeg('phi', phi);
  return Ecm / (1 + phi);
}

/**
 * Deformation of a member that is expected to crack, by interpolation between the uncracked and fully cracked
 * states, EN 1992-1-1 Eq. (7.18), (7.19): α = ζ αII + (1 − ζ) αI, ζ = 1 − β (Mcr/M)² (ζ = 0 for M ≤ Mcr); β = 1.0 for a
 * single short-term loading, 0.5 for sustained or repeated loading. Works for any deformation parameter; given the two
 * flexural stiffnesses it returns the mean curvature, and with a span and deflection coefficient K the deflection
 * δ = K L² κ (K = 5/48 simply supported with uniform load, 1/12 central point load, 1/4 cantilever with uniform load,
 * 1/3 cantilever with end load — elastic beam theory with κ at the section of maximum moment).
 * @param {object} o
 * @param {number} o.M bending moment at the critical section [N·m]
 * @param {number} o.Mcr cracking moment [N·m]
 * @param {number} o.EIuncracked flexural stiffness of the uncracked section [N·m²]
 * @param {number} o.EIcracked flexural stiffness of the fully cracked section [N·m²]
 * @param {number} [o.beta=1] load duration factor [-]
 * @param {number} [o.L] span [m]
 * @param {number} [o.K=5/48] deflection coefficient [-]
 * @returns {{zeta:number, curvature:number, EIeffective:number, deflection:(number|undefined)}} distribution coefficient [-], mean curvature
 *   [1/m], equivalent stiffness M/κ [N·m²], deflection [m] when L is given
 */
export function crackedDeflection({ M, Mcr, EIuncracked, EIcracked, beta = 1, L, K = 5 / 48 }) {
  nonNeg('M', M); nonNeg('Mcr', Mcr); pos('EIuncracked', EIuncracked); pos('EIcracked', EIcracked);
  if (!(beta > 0 && beta <= 1)) throw new RangeError('beta must lie in (0, 1]');
  const zeta = M > Mcr ? 1 - beta * (Mcr / M) ** 2 : 0;
  const curvature = zeta * M / EIcracked + (1 - zeta) * M / EIuncracked;
  return {
    zeta, curvature, EIeffective: M > 0 ? M / curvature : EIuncracked,
    deflection: L === undefined ? undefined : K * pos('L', L) ** 2 * curvature,
  };
}

/* ------------------------------------------------------------------ shear, EN 1992-1-1 §6.2 */

/**
 * Design shear resistance of a member without shear reinforcement, EN 1992-1-1 Eq. (6.2.a), (6.2.b):
 * VRd,c = [CRd,c k (100 ρl fck)^(1/3) + k1 σcp] bw d ≥ (vmin + k1 σcp) bw d, with k = 1 + √(200/d) ≤ 2 (d in mm),
 * ρl = Asl/(bw d) ≤ 0.02, σcp = NEd/Ac < 0.2 fcd, CRd,c = 0.18/γc, k1 = 0.15, vmin = 0.035 k^(3/2) fck^(1/2) (MPa)
 * (recommended values of the standard).
 * @param {object} o
 * @param {number} o.bw smallest web width in the tensile area [m]
 * @param {number} o.d effective depth [m]
 * @param {number} o.fck [Pa]
 * @param {number} o.Asl area of anchored tensile reinforcement [m²]
 * @param {number} [o.NEd=0] axial force [N], compression positive
 * @param {number} [o.Ac] concrete area [m²] (needed when NEd ≠ 0)
 * @param {number} [o.gammaC=1.5]
 * @param {number} [o.alphaCC=1] for the limit σcp < 0.2 fcd
 * @returns {{VRdc:number, k:number, rhoL:number, sigmaCp:number, vMin:number, governedByMinimum:boolean}} resistance [N], size factor, ratio,
 *   mean compressive stress [Pa], vmin [Pa]
 */
export function shearResistanceNoLinks({ bw, d, fck, Asl, NEd = 0, Ac, gammaC = 1.5, alphaCC = 1 }) {
  pos('bw', bw); pos('d', d); pos('fck', fck); nonNeg('Asl', Asl); pos('gammaC', gammaC);
  const f = fck / 1e6;
  const k = Math.min(2, 1 + Math.sqrt(0.2 / d));
  const rhoL = Math.min(0.02, Asl / (bw * d));
  let sigmaCp = 0;
  if (NEd !== 0) sigmaCp = Math.min(NEd / pos('Ac', Ac), 0.2 * alphaCC * fck / gammaC);
  const v = (0.18 / gammaC) * k * (100 * rhoL * f) ** (1 / 3) * 1e6;
  const vMin = 0.035 * k ** 1.5 * Math.sqrt(f) * 1e6;
  const a = (v + 0.15 * sigmaCp) * bw * d; const b = (vMin + 0.15 * sigmaCp) * bw * d;
  return { VRdc: Math.max(a, b, 0), k, rhoL, sigmaCp, vMin, governedByMinimum: b > a };
}

/**
 * Shear resistance of a member with vertical shear reinforcement by the variable-angle truss model, EN 1992-1-1
 * Eq. (6.8), (6.9): VRd,s = (Asw/s) z fywd cot θ, VRd,max = αcw bw z ν1 fcd/(cot θ + tan θ), ν1 = 0.6 (1 − fck/250),
 * 1 ≤ cot θ ≤ 2.5, z = 0.9 d by default. When VEd is given the links required are evaluated with the flattest strut the
 * web crushing limit allows (θ from VEd = VRd,max when cot θ = 2.5 is not sufficient).
 * @param {object} o
 * @param {number} o.bw web width [m]
 * @param {number} o.d effective depth [m]
 * @param {number} [o.z=0.9 d] inner lever arm [m]
 * @param {number} o.fck [Pa]
 * @param {number} o.fywk characteristic yield strength of the links [Pa]
 * @param {number} [o.AswOverS] link area per unit length Asw/s [m²/m]
 * @param {number} [o.cotTheta=2.5] strut inclination, 1 ≤ cot θ ≤ 2.5
 * @param {number} [o.VEd] design shear force [N] for the required reinforcement
 * @param {number} [o.gammaC=1.5]
 * @param {number} [o.gammaS=1.15]
 * @param {number} [o.alphaCC=1]
 * @param {number} [o.alphaCw=1] coefficient for the state of stress in the compression chord [-] (1 for non-prestressed members)
 * @returns {{VRds:(number|undefined), VRdmax:number, VRd:(number|undefined), nu1:number, cotTheta:number, required:(object|undefined)}} resistances [N];
 *   `required` = { AswOverS [m²/m], cotTheta, adequate } for VEd (adequate = false when VEd exceeds VRd,max at cot θ = 1)
 */
export function shearResistanceLinks({ bw, d, z, fck, fywk, AswOverS, cotTheta = 2.5, VEd, gammaC = 1.5, gammaS = 1.15, alphaCC = 1, alphaCw = 1 }) {
  pos('bw', bw); pos('d', d); pos('fck', fck); pos('fywk', fywk);
  if (!(cotTheta >= 1 && cotTheta <= 2.5)) throw new RangeError('cotTheta must lie in [1, 2.5]');
  const zz = z === undefined ? 0.9 * d : pos('z', z);
  const fcd = alphaCC * fck / gammaC; const fywd = fywk / gammaS;
  const nu1 = 0.6 * (1 - fck / 250e6);
  const vmax = (c) => alphaCw * bw * zz * nu1 * fcd / (c + 1 / c);
  const out = { VRds: undefined, VRdmax: vmax(cotTheta), VRd: undefined, nu1, cotTheta };
  if (AswOverS !== undefined) {
    out.VRds = nonNeg('AswOverS', AswOverS) * zz * fywd * cotTheta;
    out.VRd = Math.min(out.VRds, out.VRdmax);
  }
  if (VEd !== undefined) {
    nonNeg('VEd', VEd);
    let c = 2.5; let adequate = true;
    if (VEd > vmax(2.5)) {
      if (VEd > vmax(1)) { c = 1; adequate = false; }
      // VEd = αcw bw z ν1 fcd sin θ cos θ = ½ (…) sin 2θ
      else { const th = 0.5 * Math.asin(2 * VEd / (alphaCw * bw * zz * nu1 * fcd)); c = 1 / Math.tan(th); }
    }
    out.required = { AswOverS: VEd / (zz * fywd * c), cotTheta: c, adequate };
  }
  return out;
}

/* ------------------------------------------------------------------ crack width, EN 1992-1-1 §7.3.4 */

/**
 * Characteristic crack width, EN 1992-1-1 Eq. (7.8)–(7.11): wk = sr,max (εsm − εcm),
 * εsm − εcm = [σs − kt fct,eff (1 + αe ρp,eff)/ρp,eff]/Es ≥ 0.6 σs/Es,
 * sr,max = k3 c + k1 k2 k4 φ/ρp,eff (bar spacing ≤ 5(c + φ/2)), otherwise sr,max = 1.3 (h − x);
 * ρp,eff = As/Ac,eff, Ac,eff = b · min[2.5 (h − d), (h − x)/3, h/2]. Recommended values k3 = 3.4, k4 = 0.425;
 * k1 = 0.8 high-bond bars (1.6 plain), k2 = 0.5 bending (1.0 pure tension), kt = 0.6 short-term, 0.4 long-term loading.
 * @param {object} o
 * @param {number} o.sigmaS stress in the tension reinforcement on the cracked section [Pa]
 * @param {number} o.As area of tension reinforcement within Ac,eff [m²]
 * @param {number} o.b width of the tension zone [m]
 * @param {number} o.h section depth [m]
 * @param {number} o.d effective depth [m]
 * @param {number} o.x neutral axis depth of the cracked section [m]
 * @param {number} o.cover cover to the longitudinal reinforcement c [m]
 * @param {number} o.barDiameter bar diameter φ [m] (equivalent diameter for mixed sizes)
 * @param {number} o.fctEff mean tensile strength effective at cracking [Pa] (fctm, or less for early-age cracking)
 * @param {number} o.Ecm concrete modulus [Pa]
 * @param {number} [o.Es=2.0e11] [Pa]
 * @param {number} [o.kt=0.4] load duration factor
 * @param {number} [o.k1=0.8] bond factor
 * @param {number} [o.k2=0.5] strain distribution factor
 * @param {number} [o.k3=3.4]
 * @param {number} [o.k4=0.425]
 * @param {number} [o.barSpacing] spacing of the bonded bars [m]; when it exceeds 5(c + φ/2), sr,max = 1.3(h − x)
 * @returns {{wk:number, srMax:number, strainDifference:number, rhoEff:number, hcEff:number, minimumStrainGoverns:boolean}} crack width [m],
 *   maximum crack spacing [m], εsm − εcm [-], ρp,eff [-], effective tension depth [m]
 */
export function crackWidth({ sigmaS, As, b, h, d, x, cover, barDiameter, fctEff, Ecm, Es = 2.0e11, kt = 0.4, k1 = 0.8, k2 = 0.5, k3 = 3.4, k4 = 0.425, barSpacing }) {
  nonNeg('sigmaS', sigmaS); pos('As', As); pos('b', b); pos('h', h); pos('d', d); nonNeg('x', x); nonNeg('cover', cover);
  pos('barDiameter', barDiameter); pos('fctEff', fctEff); pos('Ecm', Ecm); pos('Es', Es);
  if (!(d < h)) throw new RangeError('d must be smaller than h');
  const hcEff = Math.min(2.5 * (h - d), (h - x) / 3, h / 2);
  const rhoEff = As / (b * hcEff);
  const ae = Es / Ecm;
  const full = (sigmaS - kt * fctEff / rhoEff * (1 + ae * rhoEff)) / Es;
  const minimum = 0.6 * sigmaS / Es;
  const diff = Math.max(full, minimum);
  const wide = barSpacing !== undefined && barSpacing > 5 * (cover + barDiameter / 2);
  const srMax = wide ? 1.3 * (h - x) : k3 * cover + k1 * k2 * k4 * barDiameter / rhoEff;
  return { wk: srMax * diff, srMax, strainDifference: diff, rhoEff, hcEff, minimumStrainGoverns: minimum > full };
}

/* ------------------------------------------------------------------ durability */

/**
 * Chloride content at depth x and time t from the error-function solution of Fick's second law for a constant
 * surface concentration, with a time-dependent apparent diffusion coefficient
 *   C(x, t) = Ci + (Cs − Ci) erfc[(x − Δx)/(2 √(Da(t) t))],   Da(t) = D0 (t0/t)^m
 * (Collepardi et al. 1972; ageing form as in fib Bulletin 34 "Model Code for Service Life Design" (2006), chloride model,
 * where D0 = ke·DRCM,0·kt at the reference age t0 and m is the ageing exponent). m = 0 gives constant diffusivity.
 * Concentrations in any consistent unit (e.g. % by mass of binder).
 * @param {object} o
 * @param {number} o.x depth from the surface [m]
 * @param {number} o.t exposure time [s]
 * @param {number} o.Cs surface chloride content (at depth Δx)
 * @param {number} [o.Ci=0] initial chloride content
 * @param {number} o.D0 apparent diffusion coefficient at the reference time [m²/s]
 * @param {number} [o.t0=28 days] reference time [s]
 * @param {number} [o.m=0] ageing exponent [-], 0 ≤ m < 1
 * @param {number} [o.dx=0] depth of the convection zone Δx [m]
 * @returns {{C:number, Da:number}} chloride content at (x, t) and the apparent diffusion coefficient at t [m²/s]
 */
export function chlorideProfile({ x, t, Cs, Ci = 0, D0, t0 = 28 * 86400, m = 0, dx = 0 }) {
  nonNeg('x', x); pos('t', t); pos('D0', D0); pos('t0', t0); nonNeg('dx', dx);
  if (!(m >= 0 && m < 1)) throw new RangeError('m must lie in [0, 1)');
  if (!Number.isFinite(Cs) || !Number.isFinite(Ci)) throw new RangeError('Cs and Ci must be finite');
  const Da = D0 * (t0 / t) ** m;
  const arg = Math.max(0, x - dx) / (2 * Math.sqrt(Da * t));
  return { C: Ci + (Cs - Ci) * erfc(arg), Da };
}

/**
 * Time to corrosion initiation by chloride ingress: the time at which the chloride content at the reinforcement
 * reaches the critical content, from the model of chlorideProfile solved in closed form:
 * (c − Δx)/(2 √(D0 t0^m t^(1−m))) = erfc⁻¹[(Ccrit − Ci)/(Cs − Ci)].
 * @param {object} o
 * @param {number} o.cover concrete cover to the reinforcement [m]
 * @param {number} o.Cs surface chloride content
 * @param {number} o.Ccrit critical chloride content (same unit)
 * @param {number} [o.Ci=0] initial chloride content
 * @param {number} o.D0 apparent diffusion coefficient at t0 [m²/s]
 * @param {number} [o.t0=28 days] reference time [s]
 * @param {number} [o.m=0] ageing exponent [-]
 * @param {number} [o.dx=0] convection zone depth [m]
 * @returns {number} initiation time [s] (Infinity when the surface content does not exceed the critical content)
 */
export function chlorideInitiationTime({ cover, Cs, Ccrit, Ci = 0, D0, t0 = 28 * 86400, m = 0, dx = 0 }) {
  pos('cover', cover); pos('D0', D0); pos('t0', t0); nonNeg('dx', dx);
  if (!(m >= 0 && m < 1)) throw new RangeError('m must lie in [0, 1)');
  if (!(Ccrit > Ci)) throw new RangeError('Ccrit must exceed Ci');
  if (!(Cs > Ccrit)) return Infinity;
  if (cover <= dx) return 0;
  // erfc(z) = r  →  z = erfinv(1 − r)
  const z = erfinv(1 - (Ccrit - Ci) / (Cs - Ci));
  const q = ((cover - dx) / (2 * z)) ** 2;
  return (q / (D0 * t0 ** m)) ** (1 / (1 - m));
}

/**
 * Carbonation depth by the square-root-of-time law xc = K √t (diffusion-controlled front; the simplified form of
 * the fib Bulletin 34 carbonation model) and the time for the front to reach the reinforcement, t = (c/K)².
 * K depends on the concrete and exposure and is an input (1 mm/√year = 1.780e-7 m/√s).
 * @param {object} o
 * @param {number} o.K carbonation coefficient [m/√s]
 * @param {number} [o.t] exposure time [s]
 * @param {number} [o.cover] cover [m]
 * @returns {{depth:(number|undefined), initiationTime:(number|undefined)}} xc [m] at t; time to depassivation [s] for the cover
 */
export function carbonation({ K, t, cover }) {
  pos('K', K);
  return {
    depth: t === undefined ? undefined : K * Math.sqrt(nonNeg('t', t)),
    initiationTime: cover === undefined ? undefined : (pos('cover', cover) / K) ** 2,
  };
}

/**
 * Corrosion penetration rate from the corrosion current density by Faraday's law: dx/dt = icorr M/(z F ρ).
 * For iron (M = 55.845 g/mol, z = 2) 1 µA/cm² = 0.01 A/m² corresponds to 11.6 µm/year.
 * @param {object} o
 * @param {number} o.icorr corrosion current density [A/m²] (1 µA/cm² = 0.01 A/m²)
 * @param {number} [o.molarMass=0.055845] [kg/mol]
 * @param {number} [o.valence=2] electrons per atom [-]
 * @param {number} [o.rho=7850] density of the metal [kg/m³]
 * @returns {number} penetration rate [m/s]
 */
export function corrosionRate({ icorr, molarMass = 0.055845, valence = 2, rho = 7850 }) {
  nonNeg('icorr', icorr); pos('molarMass', molarMass); pos('valence', valence); pos('rho', rho);
  return icorr * molarMass / (valence * 96485.33212 * rho);
}

/**
 * Loss of bar section by uniform corrosion after initiation: penetration p = rate·(t − ti), remaining diameter
 * d = d0 − 2p, area ratio (d/d0)². The largest pit depth is reported as R·p for a pitting factor R supplied by the
 * caller (ratio of maximum to average penetration, a material/exposure parameter; R = 1 means no pitting estimate).
 * @param {object} o
 * @param {number} o.diameter original bar diameter d0 [m]
 * @param {number} o.rate penetration rate [m/s] (see corrosionRate)
 * @param {number} o.t age [s]
 * @param {number} [o.initiationTime=0] time of depassivation [s]
 * @param {number} [o.pittingFactor=1] R [-]
 * @returns {{penetration:number, diameter:number, areaRatio:number, lossFraction:number, pitDepth:number}} p [m], remaining diameter [m],
 *   remaining / original area [-], lost fraction [-], maximum pit depth [m] (capped at d0)
 */
export function barSectionLoss({ diameter, rate, t, initiationTime = 0, pittingFactor = 1 }) {
  pos('diameter', diameter); nonNeg('rate', rate); nonNeg('t', t); nonNeg('initiationTime', initiationTime); pos('pittingFactor', pittingFactor);
  const p = rate * Math.max(0, t - initiationTime);
  const d = Math.max(0, diameter - 2 * p);
  const ar = (d / diameter) ** 2;
  return { penetration: p, diameter: d, areaRatio: ar, lossFraction: 1 - ar, pitDepth: Math.min(diameter, pittingFactor * p) };
}

/**
 * Residual moment capacity of a corroding reinforced concrete section with time: bars that carry a `diameter` lose
 * section uniformly after the initiation time (barSectionLoss) and the ultimate moment is recomputed (momentCapacity).
 * Bond deterioration, cover spalling and loss of steel ductility are not modelled.
 * @param {object} o
 * @param {object} o.section result of rcSection (bars with `diameter` corrode)
 * @param {object} o.material see interactionDiagram
 * @param {number[]} o.times ages [s]
 * @param {number} o.initiationTime time to depassivation [s]
 * @param {number} o.rate penetration rate [m/s]
 * @param {number} [o.N=0] axial force [N], compression positive
 * @returns {Array<{t:number, MRd:number, areaRatio:number, capacityRatio:number}>} moment resistance [N·m], remaining bar area ratio [-] and
 *   capacity relative to the uncorroded section [-] at each age
 */
export function residualCapacityHistory({ section, material, times, initiationTime, rate, N = 0 }) {
  if (!Array.isArray(times) && !ArrayBuffer.isView(times)) throw new RangeError('times must be an array');
  const M0 = momentCapacity({ section, material, N }).MRd;
  const out = [];
  for (const t of times) {
    let a0 = 0; let a1 = 0;
    const bars = section.bars.map((b) => {
      if (!(b.diameter > 0)) return b;
      const r = barSectionLoss({ diameter: b.diameter, rate, t, initiationTime }).areaRatio;
      a0 += b.area; a1 += b.area * r;
      return { ...b, area: b.area * r };
    });
    const MRd = momentCapacity({ section: { ...section, bars }, material, N }).MRd;
    out.push({ t, MRd, areaRatio: a0 > 0 ? a1 / a0 : 1, capacityRatio: MRd / M0 });
  }
  return out;
}

/**
 * Wall thickness loss of a steel pile (or sheet pile) by exposure zone: loss = rate·max(0, t − protection life) on each
 * exposed face, remaining thickness and, for a tubular pile, the remaining area and section modulus (loss taken from
 * the outside). Corrosion rates per zone (atmospheric, splash, tidal, low water, immersion, buried) are inputs — take
 * them from the governing standard (e.g. EN 1993-5 Tables 4-1 and 4-2 give total losses for 5 to 100 years) or from
 * site measurements.
 * @param {object} o
 * @param {Array<{name:string, rate:number, protectionLife?:number, faces?:number}>} o.zones corrosion rate per exposed face [m/s], life of
 *   coating or cathodic protection before corrosion starts [s] (default 0), number of corroding faces (default 1)
 * @param {number} o.t design life [s]
 * @param {number} o.thickness original wall thickness [m]
 * @param {number} [o.D] outer diameter of a tubular pile [m]
 * @returns {Array<{name:string, loss:number, remaining:number, thicknessRatio:number, areaRatio:(number|undefined), modulusRatio:(number|undefined)}>}
 *   total thickness loss [m], remaining thickness [m] and ratios to the original [-] for each zone
 */
export function pileCorrosionLoss({ zones, t, thickness, D }) {
  if (!Array.isArray(zones) || zones.length === 0) throw new RangeError('zones must be a non-empty array');
  nonNeg('t', t); pos('thickness', thickness);
  if (D !== undefined) { pos('D', D); if (thickness > D / 2) throw new RangeError('thickness must not exceed D/2'); }
  const tube = (Do, tt) => { const di = Do - 2 * tt; return { A: Math.PI / 4 * (Do * Do - di * di), W: Math.PI / 32 * (Do ** 4 - di ** 4) / Do }; };
  const ref = D === undefined ? null : tube(D, thickness);
  return zones.map((z, i) => {
    nonNeg(`zones[${i}].rate`, z.rate);
    const faces = z.faces === undefined ? 1 : z.faces;
    const loss = Math.min(thickness, faces * z.rate * Math.max(0, t - (z.protectionLife || 0)));
    const remaining = thickness - loss;
    const out = { name: z.name, loss, remaining, thicknessRatio: remaining / thickness, areaRatio: undefined, modulusRatio: undefined };
    if (ref) {
      // one face: metal lost from the outside (outer diameter shrinks); two faces: shared between outside and inside
      const outer = faces >= 2 ? loss / 2 : loss;
      const now = remaining > 0 ? tube(D - 2 * outer, remaining) : { A: 0, W: 0 };
      out.areaRatio = now.A / ref.A; out.modulusRatio = now.W / ref.W;
    }
    return out;
  });
}
