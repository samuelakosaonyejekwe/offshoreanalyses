/**
 * Settlement of shallow foundations: elastic, strain-influence, consolidation (magnitude and rate).
 *
 * SI units: m, Pa, N, s, N/m^3; hydraulic conductivity m/s; cv, ch in m^2/s.
 *
 * References
 *   Boussinesq J. (1885); Newmark N.M. (1935) Simplified computation of vertical pressures in elastic
 *     foundations, Univ. Illinois Eng. Exp. Stn Circ. 24.
 *   Schleicher F. (1926); Timoshenko S. & Goodier J.N., Theory of Elasticity; Poulos H.G. & Davis E.H.
 *     (1974) Elastic Solutions for Soil and Rock Mechanics; Steinbrenner W. (1934).
 *   Schmertmann J.H., Hartman J.P. & Brown P.R. (1978) Improved strain influence factor diagrams, JGED 104.
 *   Janbu N. (1963) Soil compressibility as determined by oedometer and triaxial tests, ECSMFE Wiesbaden.
 *   Terzaghi K. (1943); Barron R.A. (1948) Trans. ASCE 113; Hansbo S. (1981); Carrillo N. (1942).
 *   Biot M.A. (1941) J. Appl. Phys. 12; Verruijt A. (2013) Theory and Problems of Poroelasticity.
 */
import { brent } from '../math/roots.js';
import { solveTridiagonal } from '../math/linalg.js';

const PI = Math.PI;

function need(name, v, min = -Infinity, max = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < min || v > max) {
    throw new RangeError(`${name} must be a number in [${min}, ${max}] (got ${v})`);
  }
  return v;
}

/* ------------------------------------------------------------------ stress increments */

/**
 * Vertical stress influence factor below the corner of a uniformly loaded rectangle on an elastic
 * half-space (Newmark 1935 integration of Boussinesq):
 *   I = (1/4π) [ 2mn√(m²+n²+1)/(m²+n²+m²n²+1) · (m²+n²+2)/(m²+n²+1) + atan( 2mn√(m²+n²+1)/(m²+n²+1−m²n²) ) ]
 * with m = B/z, n = L/z and the arctangent taken in (0, π).
 * @param {number} B side [m]
 * @param {number} L other side [m]
 * @param {number} z depth below the loaded surface [m], ≥ 0
 * @returns {number} Δσz / q [-] (0.25 at z = 0)
 */
export function boussinesqRectCorner(B, L, z) {
  need('B', B, 0);
  need('L', L, 0);
  need('z', z, 0);
  if (B === 0 || L === 0) return 0;
  if (z === 0) return 0.25;
  const m = B / z;
  const n = L / z;
  const s = m * m + n * n + 1;
  const r = 2 * m * n * Math.sqrt(s);
  return ((r / (s + m * m * n * n)) * ((s + 1) / s) + Math.atan2(r, s - m * m * n * n)) / (4 * PI);
}

/**
 * Vertical stress increment at any point below (or beside) a uniformly loaded rectangle, by superposition
 * of four corner solutions.
 * @param {object} o
 * @param {number} o.q applied pressure [Pa]
 * @param {number} o.B footing width along x [m]
 * @param {number} o.L footing length along y [m]
 * @param {number} o.z depth [m]
 * @param {number} [o.x=0] offset from the centre along x [m]
 * @param {number} [o.y=0] offset from the centre along y [m]
 * @returns {number} Δσz [Pa]
 */
export function stressBelowRect({ q, B, L, z, x = 0, y = 0 }) {
  need('B', B, 1e-12);
  need('L', L, 1e-12);
  need('z', z, 0);
  const part = (a, b) => Math.sign(a) * Math.sign(b) * boussinesqRectCorner(Math.abs(a), Math.abs(b), z);
  const x1 = B / 2 - x;
  const x2 = B / 2 + x;
  const y1 = L / 2 - y;
  const y2 = L / 2 + y;
  return q * (part(x1, y1) + part(x1, y2) + part(x2, y1) + part(x2, y2));
}

/**
 * Vertical stress increment on the axis of a uniformly loaded circle (Boussinesq):
 *   Δσz = q [1 − (1 + (R/z)²)^(−3/2)].
 * @param {number} q pressure [Pa]
 * @param {number} R radius [m]
 * @param {number} z depth [m]
 * @returns {number} Δσz [Pa]
 */
export function stressBelowCircleCentre(q, R, z) {
  need('R', R, 1e-12);
  need('z', z, 0);
  if (z === 0) return q;
  return q * (1 - (1 + (R / z) ** 2) ** -1.5);
}

/**
 * Boussinesq point load: Δσz = 3 Q z³ / (2π (r² + z²)^(5/2)).
 * @param {number} Q load [N]
 * @param {number} r horizontal offset [m]
 * @param {number} z depth [m], > 0
 * @returns {number} Δσz [Pa]
 */
export function stressPointLoad(Q, r, z) {
  need('r', r, 0);
  need('z', z, 1e-12);
  return (3 * Q * z ** 3) / (2 * PI * (r * r + z * z) ** 2.5);
}

/**
 * 2:1 (vertical : horizontal) load-spread approximation: Δσz = q B L / ((B + z)(L + z)).
 * @param {object} o q [Pa], B [m], L [m] (Infinity for a strip), z [m]
 * @returns {number} Δσz [Pa]
 */
export function stress2to1({ q, B, L = Infinity, z }) {
  need('B', B, 1e-12);
  need('z', z, 0);
  return Number.isFinite(L) ? (q * B * L) / ((B + z) * (L + z)) : (q * B) / (B + z);
}

/* ------------------------------------------------------------------ immediate settlement */

/**
 * Surface settlement influence factor at the corner of a flexible uniformly loaded rectangle on an elastic
 * half-space (Schleicher 1926; Timoshenko & Goodier):  s = q B (1 − ν²) I / E,
 *   I = (1/π) [ m ln((1 + √(1 + m²))/m) + ln(m + √(1 + m²)) ],  m = L/B.
 * (0.561 for a square.)
 * @param {number} m aspect ratio L/B [-], > 0
 * @returns {number} I [-]
 */
export function influenceRectCorner(m) {
  need('m', m, 1e-9);
  const r = Math.sqrt(1 + m * m);
  return (m * Math.log((1 + r) / m) + Math.log(m + r)) / PI;
}

/**
 * Steinbrenner (1934) factors for the corner of a flexible rectangle B × L on a layer of thickness H over a
 * rigid base:  s = q B (1 − ν²) [F1 + (1 − 2ν)/(1 − ν) F2] / E, with M = L/B, N = H/B:
 *   F1 = (1/π) { M ln[ (1 + √(M²+1)) √(M²+N²) / (M (1 + √(M²+N²+1))) ] + ln[ (M + √(M²+1)) √(1+N²) / (M + √(M²+N²+1)) ] }
 *   F2 = (N/2π) atan( M / (N √(M²+N²+1)) ).
 * @param {number} M L/B [-]
 * @param {number} N H/B [-]
 * @returns {{F1:number,F2:number}} [-]
 */
export function steinbrennerFactors(M, N) {
  need('M', M, 1e-9);
  need('N', N, 0);
  if (N === 0) return { F1: 0, F2: 0 };
  const a = Math.sqrt(M * M + 1);
  const b = Math.sqrt(M * M + N * N);
  const c = Math.sqrt(M * M + N * N + 1);
  const F1 = (M * Math.log(((1 + a) * b) / (M * (1 + c))) + Math.log(((M + a) * Math.sqrt(1 + N * N)) / (M + c))) / PI;
  const F2 = (N / (2 * PI)) * Math.atan(M / (N * c));
  return { F1, F2 };
}

/**
 * Immediate (elastic) settlement of a uniformly loaded footing on a homogeneous elastic half-space or
 * finite layer:  s = q B (1 − ν²) I / E.
 *   Rectangle, flexible: corner and centre from influenceRectCorner (centre = 2 × corner factor of the
 *     quarter rectangle); average by numerical integration of the corner solution over the footing
 *     (first principles; 0.95 for a square).
 *   Rectangle, rigid: 0.93 × flexible average (Schleicher 1926; Bowles 1996 §5.6) – an approximation.
 *   Circle (B = diameter): flexible centre I = 1, edge 2/π, average 8/(3π); rigid π/4 (exact, Boussinesq).
 *   Finite layer (H given, rectangle only): Steinbrenner factors replace the half-space factor.
 * @param {object} o
 * @param {number} o.q net applied pressure [Pa]
 * @param {number} o.B width or diameter [m]
 * @param {number} [o.L] length [m] (rectangle; ≥ B)
 * @param {number} o.E Young's modulus [Pa] (undrained Eu with ν = 0.5 for immediate settlement of clay)
 * @param {number} o.nu Poisson's ratio [-]
 * @param {'rect'|'circle'} [o.shape='rect']
 * @param {'centre'|'corner'|'edge'|'average'|'rigid'} [o.point='rigid'] ('edge' for circles only)
 * @param {number} [o.H=Infinity] compressible layer thickness below the base [m]
 * @returns {{s:number,I:number}} settlement [m]; influence factor [-]
 */
export function elasticSettlement({ q, B, L, E, nu, shape = 'rect', point = 'rigid', H = Infinity }) {
  need('B', B, 1e-12);
  need('E', E, 1e-12);
  need('nu', nu, 0, 0.5);
  let I;
  if (shape === 'circle') {
    I = { centre: 1, edge: 2 / PI, average: 8 / (3 * PI), rigid: PI / 4 }[point];
    if (I === undefined) throw new RangeError(`point '${point}' is not available for a circle`);
    if (Number.isFinite(H)) throw new RangeError('finite layer thickness H is implemented for rectangles only');
  } else if (shape === 'rect') {
    need('L', L, B);
    const cornerI = (b, l) => {
      // settlement at the corner of a b × l rectangle, expressed as a factor on the full width B
      if (b <= 0 || l <= 0) return 0;
      const short = Math.min(b, l);
      const long = Math.max(b, l);
      if (!Number.isFinite(H)) return (short / B) * influenceRectCorner(long / short);
      const f = steinbrennerFactors(long / short, H / short);
      return (short / B) * (f.F1 + ((1 - 2 * nu) / (1 - nu)) * f.F2);
    };
    const at = (x, y) => cornerI(B / 2 - x, L / 2 - y) + cornerI(B / 2 - x, L / 2 + y) + cornerI(B / 2 + x, L / 2 - y) + cornerI(B / 2 + x, L / 2 + y);
    if (point === 'corner') I = cornerI(B, L);
    else if (point === 'centre') I = at(0, 0);
    else if (point === 'average' || point === 'rigid') {
      // 16-point Gauss–Legendre on a quadrant (the integrand is smooth inside the footing)
      const gx = [0.0950125098376374, 0.2816035507792589, 0.4580167776572274, 0.6178762444026438, 0.755404408355003, 0.8656312023878318, 0.9445750230732326, 0.9894009349916499];
      const gw = [0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767, 0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541];
      const xs = [];
      const ws = [];
      for (let i = 0; i < 8; i++) { xs.push(0.5 - 0.5 * gx[i], 0.5 + 0.5 * gx[i]); ws.push(0.5 * gw[i], 0.5 * gw[i]); }
      let sum = 0;
      for (let i = 0; i < 16; i++) for (let j = 0; j < 16; j++) sum += ws[i] * ws[j] * at((xs[i] * B) / 2, (xs[j] * L) / 2);
      I = point === 'rigid' ? 0.93 * sum : sum;
    } else throw new RangeError(`point '${point}' is not available for a rectangle`);
  } else throw new RangeError(`shape must be 'rect' or 'circle' (got ${shape})`);
  return { s: (q * B * (1 - nu * nu) * I) / E, I };
}

/**
 * Rotation of a rigid footing under moment on an elastic half-space.
 *   Circle (exact, Borowicka 1943): θ = 3 (1 − ν) M / (8 G R³) = 6 (1 − ν²) M / (E D³).
 *   Rectangle: tan θ = M (1 − ν²) Iθ / (B² L E), with Iθ an argument; the default
 *   Iθ = 16 / [π (1 + 0.22 B/L)] is the expression quoted by Bowles (1996, eq. 5-17, after Taylor 1967) for
 *   a rigid base and should be treated as approximate.
 * @param {object} o
 * @param {number} o.M moment [N·m]
 * @param {number} o.B dimension in the plane of rotation, or diameter [m]
 * @param {number} [o.L] other dimension [m] (rectangle)
 * @param {number} o.E [Pa]
 * @param {number} o.nu [-]
 * @param {'rect'|'circle'} [o.shape='rect']
 * @param {number} [o.Itheta] rotation influence factor [-] (rectangle)
 * @returns {{theta:number,differential:number}} rotation [rad]; edge-to-edge differential settlement [m]
 */
export function footingRotation({ M, B, L, E, nu, shape = 'rect', Itheta }) {
  need('B', B, 1e-12);
  need('E', E, 1e-12);
  need('nu', nu, 0, 0.5);
  let theta;
  if (shape === 'circle') {
    const G = E / (2 * (1 + nu));
    theta = (3 * (1 - nu) * M) / (8 * G * (B / 2) ** 3);
  } else {
    need('L', L, 1e-12);
    const It = Itheta !== undefined ? need('Itheta', Itheta, 0) : 16 / (PI * (1 + (0.22 * B) / L));
    theta = Math.atan((M * (1 - nu * nu) * It) / (B * B * L * E));
  }
  return { theta, differential: Math.tan(theta) * B };
}

/**
 * Schmertmann, Hartman & Brown (1978) strain-influence settlement of footings on sand:
 *   s = C1 C2 Δp Σ (Iz / Es) Δz,
 *   C1 = 1 − 0.5 σ′v0/Δp ≥ 0.5,  C2 = 1 + 0.2 log10(t / 0.1 yr),
 *   Iz: axisymmetric (L/B = 1) 0.1 at z = 0, peak at B/2, zero at 2B; plane strain (L/B ≥ 10) 0.2 at z = 0,
 *   peak at B, zero at 4B; peak Izp = 0.5 + 0.1 √(Δp / σ′vp), σ′vp = initial effective stress at the depth
 *   of the peak. Intermediate L/B is interpolated linearly between the two cases (Schmertmann's advice).
 * Es is supplied per layer (Schmertmann: 2.5 qc axisymmetric, 3.5 qc plane strain).
 * @param {object} o
 * @param {number} o.B footing width [m]
 * @param {number} [o.L=o.B] length [m]
 * @param {number} o.dp net foundation pressure Δp [Pa]
 * @param {number} o.sigmaV0 initial effective overburden at foundation level [Pa]
 * @param {number} o.gamma effective unit weight below the base [N/m^3] (for σ′vp)
 * @param {{zTop:number,zBottom:number,Es:number}[]} o.layers depths below the base [m], Es [Pa]
 * @param {number} [o.timeYears=0.1] time for the creep factor C2 [yr] (≥ 0.1)
 * @param {number} [o.nSub=40] integration sub-intervals per layer
 * @returns {{s:number,C1:number,C2:number,Izp:number,sumIzEs:number}} s [m]; Σ(Iz/Es)Δz [m/Pa]
 */
export function schmertmannSettlement({ B, L, dp, sigmaV0, gamma, layers, timeYears = 0.1, nSub = 40 }) {
  need('B', B, 1e-12);
  const Lx = L === undefined ? B : need('L', L, B);
  need('dp', dp, 1e-12);
  need('sigmaV0', sigmaV0, 0);
  need('gamma', gamma, 0);
  need('timeYears', timeYears, 0.1);
  if (!Array.isArray(layers) || layers.length === 0) throw new RangeError('layers must be a non-empty array');
  const w = Math.min(Math.max((Lx / B - 1) / 9, 0), 1); // 0 axisymmetric … 1 plane strain
  const iz = (z) => {
    const one = (iz0, zp, zmax, izp) => (z <= zp ? iz0 + ((izp - iz0) * z) / zp : z >= zmax ? 0 : (izp * (zmax - z)) / (zmax - zp));
    const izpA = 0.5 + 0.1 * Math.sqrt(dp / (sigmaV0 + (gamma * B) / 2));
    const izpP = 0.5 + 0.1 * Math.sqrt(dp / (sigmaV0 + gamma * B));
    return (1 - w) * one(0.1, B / 2, 2 * B, izpA) + w * one(0.2, B, 4 * B, izpP);
  };
  let sum = 0;
  for (const ly of layers) {
    need('layer.Es', ly.Es, 1e-12);
    if (!(ly.zBottom > ly.zTop) || ly.zTop < 0) throw new RangeError('layers need 0 ≤ zTop < zBottom');
    const dz = (ly.zBottom - ly.zTop) / nSub;
    for (let i = 0; i < nSub; i++) sum += (iz(ly.zTop + (i + 0.5) * dz) / ly.Es) * dz;
  }
  const C1 = Math.max(1 - (0.5 * sigmaV0) / dp, 0.5);
  const C2 = 1 + 0.2 * Math.log10(timeYears / 0.1);
  const Izp = (1 - w) * (0.5 + 0.1 * Math.sqrt(dp / (sigmaV0 + (gamma * B) / 2))) + w * (0.5 + 0.1 * Math.sqrt(dp / (sigmaV0 + gamma * B)));
  return { s: C1 * C2 * dp * sum, C1, C2, Izp, sumIzEs: sum };
}

/**
 * Janbu (1963) tangent-modulus strain for a stress change σ′0 → σ′1 with M = m σr (σ′/σr)^(1−j):
 *   ε = [ (σ′1/σr)^j − (σ′0/σr)^j ] / (m j)   (j ≠ 0);   ε = ln(σ′1/σ′0) / m   (j = 0, NC clay).
 * Typical: j = 1 over-consolidated clay/rock (linear), 0.5 sand and silt, 0 normally consolidated clay.
 * @param {object} o
 * @param {number} o.sigma0 initial vertical effective stress [Pa]
 * @param {number} o.sigma1 final vertical effective stress [Pa]
 * @param {number} o.m modulus number [-]
 * @param {number} o.j stress exponent [-] (Janbu's a = 1 − j)
 * @param {number} [o.sigmaRef=100e3] reference stress [Pa]
 * @returns {number} vertical strain [-]
 */
export function janbuStrain({ sigma0, sigma1, m, j, sigmaRef = 100e3 }) {
  need('sigma0', sigma0, 1e-9);
  need('sigma1', sigma1, 1e-9);
  need('m', m, 1e-9);
  need('j', j, 0, 1);
  if (j === 0) return Math.log(sigma1 / sigma0) / m;
  return ((sigma1 / sigmaRef) ** j - (sigma0 / sigmaRef) ** j) / (m * j);
}

/* ------------------------------------------------------------------ consolidation magnitude */

/**
 * One-dimensional primary consolidation strain of a clay element (e–log σ′ model):
 *   σ′f ≤ σ′p:  ε = Cr/(1+e0) log10(σ′f/σ′0);
 *   σ′0 < σ′p < σ′f:  ε = Cr/(1+e0) log10(σ′p/σ′0) + Cc/(1+e0) log10(σ′f/σ′p);
 *   σ′0 ≥ σ′p (normally consolidated):  ε = Cc/(1+e0) log10(σ′f/σ′0).
 * @param {object} o sigma0, dSigma, sigmaP [Pa]; Cc, Cr, e0 [-]
 * @returns {number} vertical strain [-]
 */
export function oedometerStrain({ sigma0, dSigma, sigmaP, Cc, Cr, e0 }) {
  need('sigma0', sigma0, 1e-9);
  need('dSigma', dSigma, 0);
  need('Cc', Cc, 0);
  need('Cr', Cr, 0);
  need('e0', e0, 0);
  const sp = Math.max(sigmaP === undefined ? sigma0 : sigmaP, sigma0);
  const sf = sigma0 + dSigma;
  if (sf <= sp) return (Cr / (1 + e0)) * Math.log10(sf / sigma0);
  return (Cr / (1 + e0)) * Math.log10(sp / sigma0) + (Cc / (1 + e0)) * Math.log10(sf / sp);
}

/**
 * Primary consolidation settlement under a rectangular or circular footing by layer summation. Each layer
 * is split into sub-layers; the stress increment at mid-height comes from Boussinesq (rectangle centre by
 * corner integration, circle axis) or the 2:1 method.
 * @param {object} o
 * @param {number} o.q net foundation pressure [Pa]
 * @param {number} o.B width or diameter [m]
 * @param {number} [o.L=o.B] length [m]
 * @param {'rect'|'circle'} [o.shape='rect']
 * @param {'boussinesq'|'2:1'} [o.stressMethod='boussinesq']
 * @param {number} [o.sigmaV0Base=0] effective overburden at foundation level [Pa]
 * @param {{thickness:number,gamma:number,Cc:number,Cr:number,e0:number,ocr?:number,sigmaP?:number}[]} o.layers
 *   from the base downward; gamma effective unit weight [N/m^3]; either ocr [-] or sigmaP [Pa]
 * @param {number} [o.nSub=10] sub-layers per layer
 * @returns {{s:number,layers:{s:number,zMid:number}[],profile:{z:number,sigma0:number,dSigma:number,strain:number}[]}}
 *   s total settlement [m]
 */
export function consolidationSettlement({ q, B, L, shape = 'rect', stressMethod = 'boussinesq', sigmaV0Base = 0, layers, nSub = 10 }) {
  need('B', B, 1e-12);
  const Lx = L === undefined ? B : L;
  if (!Array.isArray(layers) || layers.length === 0) throw new RangeError('layers must be a non-empty array');
  const ds = (z) => {
    if (stressMethod === '2:1') return stress2to1({ q, B, L: Lx, z });
    if (shape === 'circle') return stressBelowCircleCentre(q, B / 2, z);
    return 4 * q * boussinesqRectCorner(B / 2, Lx / 2, z);
  };
  let zTop = 0;
  let sig = sigmaV0Base;
  let total = 0;
  const out = [];
  const profile = [];
  for (const ly of layers) {
    need('layer.thickness', ly.thickness, 1e-12);
    need('layer.gamma', ly.gamma, 0);
    const dz = ly.thickness / nSub;
    let sl = 0;
    for (let i = 0; i < nSub; i++) {
      const z = zTop + (i + 0.5) * dz;
      const s0 = sig + ly.gamma * (i + 0.5) * dz;
      const sp = ly.sigmaP !== undefined ? ly.sigmaP : (ly.ocr !== undefined ? ly.ocr : 1) * s0;
      const d = ds(z);
      const strain = oedometerStrain({ sigma0: s0, dSigma: d, sigmaP: sp, Cc: ly.Cc, Cr: ly.Cr, e0: ly.e0 });
      sl += strain * dz;
      profile.push({ z, sigma0: s0, dSigma: d, strain });
    }
    out.push({ s: sl, zMid: zTop + ly.thickness / 2 });
    total += sl;
    sig += ly.gamma * ly.thickness;
    zTop += ly.thickness;
  }
  return { s: total, layers: out, profile };
}

/**
 * Secondary compression: s = Cα/(1 + e0) · H · log10(t2/t1)  (Cα = Δe per log cycle of time).
 * @param {object} o Calpha [-], e0 [-], H layer thickness [m], t1 end of primary [s], t2 [s]
 * @returns {number} settlement [m]
 */
export function secondaryCompression({ Calpha, e0, H, t1, t2 }) {
  need('Calpha', Calpha, 0);
  need('e0', e0, 0);
  need('H', H, 0);
  need('t1', t1, 1e-12);
  need('t2', t2, 1e-12);
  return t2 <= t1 ? 0 : (Calpha / (1 + e0)) * H * Math.log10(t2 / t1);
}

/* ------------------------------------------------------------------ consolidation rate */

/**
 * Terzaghi average degree of consolidation for uniform initial excess pore pressure (exact series):
 *   U = 1 − Σ_{m=0}^{∞} (2/M²) exp(−M² Tv),  M = (2m + 1) π/2,  Tv = cv t / Hdr².
 * @param {number} Tv time factor [-], ≥ 0
 * @returns {number} U [-] in [0, 1]
 */
export function consolidationU(Tv) {
  need('Tv', Tv, 0);
  if (Tv === 0) return 0;
  if (Tv < 1e-4) return Math.sqrt((4 * Tv) / PI); // early-time asymptote, error < 1e-12 here
  let sum = 0;
  for (let m = 0; m < 200; m++) {
    const M = (2 * m + 1) * (PI / 2);
    const term = (2 / (M * M)) * Math.exp(-M * M * Tv);
    sum += term;
    if (term < 1e-17) break;
  }
  return Math.min(Math.max(1 - sum, 0), 1);
}

/**
 * Inverse of consolidationU: time factor for a given average degree of consolidation.
 * @param {number} U [-], 0 ≤ U < 1
 * @returns {number} Tv [-]
 */
export function consolidationTv(U) {
  need('U', U, 0, 1 - 1e-12);
  if (U === 0) return 0;
  if (U < 0.02) return (PI / 4) * U * U;
  return brent((t) => consolidationU(t) - U, 1e-5, 12, { tol: 1e-13 });
}

/**
 * Excess pore pressure ratio at normalised depth Z = z/Hdr (exact Terzaghi series):
 *   u/u0 = Σ (2/M) sin(M Z) exp(−M² Tv).
 * @param {number} Z z/Hdr [-] (0 at the drained face, 1 at the impervious face / mid-plane)
 * @param {number} Tv [-]
 * @returns {number} u/u0 [-]
 */
export function excessPorePressureRatio(Z, Tv) {
  need('Z', Z, 0, 2);
  need('Tv', Tv, 0);
  if (Tv === 0) return Z > 0 && Z < 2 ? 1 : 0;
  let sum = 0;
  for (let m = 0; m < 2000; m++) {
    const M = (2 * m + 1) * (PI / 2);
    const e = Math.exp(-M * M * Tv);
    sum += (2 / M) * Math.sin(M * Z) * e;
    if (e < 1e-17) break;
  }
  return sum;
}

/**
 * Radial consolidation to a vertical drain.
 *   Barron (1948) equal-strain, ideal drain:  Ur = 1 − exp(−8 Th / F(n)),  Th = ch t / de²,
 *     F(n) = n²/(n² − 1) ln n − (3n² − 1)/(4n²),  n = de/dw.
 *   Hansbo (1981) with smear and well resistance (used when s > 1 or qw is finite):
 *     F = ln(n/s) + (kh/ks) ln s − 0.75 + π z (2 l − z) kh / qw.
 * de = 1.05 × spacing (triangular) or 1.128 × spacing (square).
 * @param {object} o
 * @param {number} o.ch horizontal coefficient of consolidation [m^2/s]
 * @param {number} o.t time [s]
 * @param {number} o.de diameter of the drained cylinder [m]
 * @param {number} o.dw equivalent drain diameter [m]
 * @param {number} [o.s=1] smear ratio ds/dw [-]
 * @param {number} [o.khOverKs=1] permeability ratio undisturbed/smear [-]
 * @param {number} [o.qw=Infinity] drain discharge capacity [m^3/s]
 * @param {number} [o.kh] horizontal permeability [m/s] (needed with qw)
 * @param {number} [o.z=0] depth along the drain [m]
 * @param {number} [o.l=0] drain length (drainage path) [m]
 * @returns {{Ur:number,F:number,Th:number}} [-]
 */
export function radialConsolidation({ ch, t, de, dw, s = 1, khOverKs = 1, qw = Infinity, kh, z = 0, l = 0 }) {
  need('ch', ch, 0);
  need('t', t, 0);
  need('dw', dw, 1e-12);
  need('de', de, dw * (1 + 1e-9));
  need('s', s, 1);
  const n = de / dw;
  let F;
  if (s === 1 && !Number.isFinite(qw)) F = ((n * n) / (n * n - 1)) * Math.log(n) - (3 * n * n - 1) / (4 * n * n);
  else {
    F = Math.log(n / s) + khOverKs * Math.log(s) - 0.75;
    if (Number.isFinite(qw)) F += (PI * z * (2 * l - z) * need('kh', kh, 0)) / qw;
  }
  const Th = (ch * t) / (de * de);
  return { Ur: 1 - Math.exp((-8 * Th) / F), F, Th };
}

/**
 * Carrillo (1942) combination of vertical and radial drainage: U = 1 − (1 − Uv)(1 − Ur).
 * @param {number} Uv [-]
 * @param {number} Ur [-]
 * @returns {number} U [-]
 */
export function combinedConsolidation(Uv, Ur) {
  need('Uv', Uv, 0, 1);
  need('Ur', Ur, 0, 1);
  return 1 - (1 - Uv) * (1 - Ur);
}

/**
 * Time–settlement curve: s(t) = si + U(Tv) sc + secondary compression after tp.
 * @param {object} o
 * @param {number[]} o.times [s]
 * @param {number} [o.si=0] immediate settlement [m]
 * @param {number} o.sc final primary consolidation settlement [m]
 * @param {number} o.cv [m^2/s]
 * @param {number} o.Hdr drainage path length [m]
 * @param {number} [o.Calpha=0] [-]
 * @param {number} [o.e0=1] [-]
 * @param {number} [o.H=0] layer thickness for secondary compression [m]
 * @param {number} [o.tp] end-of-primary time [s] (default: time for U = 95 %)
 * @returns {{t:number[],s:number[],U:number[],tp:number}} settlements [m]
 */
export function timeSettlementCurve({ times, si = 0, sc, cv, Hdr, Calpha = 0, e0 = 1, H = 0, tp }) {
  need('cv', cv, 1e-18);
  need('Hdr', Hdr, 1e-12);
  const tP = tp !== undefined ? tp : (consolidationTv(0.95) * Hdr * Hdr) / cv;
  const s = [];
  const U = [];
  for (const t of times) {
    const u = consolidationU((cv * t) / (Hdr * Hdr));
    U.push(u);
    s.push(si + u * sc + (Calpha > 0 && t > tP ? secondaryCompression({ Calpha, e0, H, t1: tP, t2: t }) : 0));
  }
  return { t: Array.from(times), s, U, tp: tP };
}

/**
 * Finite-volume solution of one-dimensional consolidation in a layered column under time-dependent load:
 *   mv ∂u/∂t = ∂/∂z ( (k/γw) ∂u/∂z ) + mv dq/dt
 * Backward-Euler (θ = 1) or Crank–Nicolson (θ = 0.5) in time, harmonic-mean conductance at interfaces,
 * so layer contrasts in k and mv are handled conservatively. Second-order in space.
 * @param {object} o
 * @param {{thickness:number,k:number,mv:number}[]} o.layers top to bottom; k [m/s]; mv [1/Pa]
 * @param {(t:number)=>number} o.load applied total stress q(t) [Pa] (load(0) is applied instantaneously)
 * @param {number} o.tEnd [s]
 * @param {number} [o.nSteps=400]
 * @param {number} [o.nCells=40] cells per layer
 * @param {'drained'|'undrained'} [o.top='drained']
 * @param {'drained'|'undrained'} [o.bottom='undrained']
 * @param {number} [o.gammaW=9810] unit weight of water [N/m^3]
 * @param {number} [o.theta=1] time-integration parameter [-]
 * @returns {{t:Float64Array,settlement:Float64Array,U:Float64Array,z:Float64Array,u:Float64Array,sFinal:number}}
 *   settlement(t) [m]; U = settlement / final settlement under the current load [-]; z cell centres [m];
 *   u excess pore pressure at tEnd [Pa]; sFinal = Σ mv q(tEnd) h [m]
 */
export function consolidationFD({ layers, load, tEnd, nSteps = 400, nCells = 40, top = 'drained', bottom = 'undrained', gammaW = 9810, theta = 1 }) {
  if (!Array.isArray(layers) || layers.length === 0) throw new RangeError('layers must be a non-empty array');
  need('tEnd', tEnd, 1e-12);
  need('theta', theta, 0.5, 1);
  if (typeof load !== 'function') throw new RangeError('load must be a function of time');
  const N = layers.length * nCells;
  const dz = new Float64Array(N);
  const mv = new Float64Array(N);
  const kk = new Float64Array(N);
  const z = new Float64Array(N);
  let zc = 0;
  let idx = 0;
  for (const ly of layers) {
    need('layer.thickness', ly.thickness, 1e-12);
    need('layer.k', ly.k, 1e-30);
    need('layer.mv', ly.mv, 1e-30);
    const h = ly.thickness / nCells;
    for (let i = 0; i < nCells; i++) { dz[idx] = h; mv[idx] = ly.mv; kk[idx] = ly.k / gammaW; z[idx] = zc + (i + 0.5) * h; idx++; }
    zc += ly.thickness;
  }
  // conductances between cell centres (and to the boundaries)
  const c = new Float64Array(N + 1);
  for (let i = 1; i < N; i++) c[i] = 1 / (dz[i - 1] / (2 * kk[i - 1]) + dz[i] / (2 * kk[i]));
  c[0] = top === 'drained' ? (2 * kk[0]) / dz[0] : 0;
  c[N] = bottom === 'drained' ? (2 * kk[N - 1]) / dz[N - 1] : 0;
  const dt = tEnd / nSteps;
  const u = new Float64Array(N).fill(load(0));
  const lower = new Float64Array(N);
  const diag = new Float64Array(N);
  const upper = new Float64Array(N);
  const rhs = new Float64Array(N);
  for (let i = 0; i < N; i++) {
    const cap = (mv[i] * dz[i]) / dt;
    lower[i] = i > 0 ? -theta * c[i] : 0;
    upper[i] = i < N - 1 ? -theta * c[i + 1] : 0;
    diag[i] = cap + theta * (c[i] + c[i + 1]);
  }
  const t = new Float64Array(nSteps + 1);
  const settlement = new Float64Array(nSteps + 1);
  const U = new Float64Array(nSteps + 1);
  let sumMvH = 0;
  for (let i = 0; i < N; i++) sumMvH += mv[i] * dz[i];
  const settle = (q) => { let s = 0; for (let i = 0; i < N; i++) s += mv[i] * dz[i] * (q - u[i]); return s; };
  let qPrev = load(0);
  settlement[0] = settle(qPrev);
  for (let n = 1; n <= nSteps; n++) {
    const tn = n * dt;
    const qNow = load(tn);
    for (let i = 0; i < N; i++) {
      const cap = (mv[i] * dz[i]) / dt;
      const flux = (i > 0 ? c[i] * (u[i - 1] - u[i]) : -c[0] * u[0]) + (i < N - 1 ? c[i + 1] * (u[i + 1] - u[i]) : -c[N] * u[N - 1]);
      rhs[i] = cap * u[i] + (1 - theta) * flux + cap * (qNow - qPrev);
    }
    const sol = solveTridiagonal(lower, diag, upper, rhs);
    for (let i = 0; i < N; i++) u[i] = sol[i];
    qPrev = qNow;
    t[n] = tn;
    settlement[n] = settle(qNow);
    U[n] = qNow !== 0 ? settlement[n] / (sumMvH * qNow) : 0;
  }
  return { t, settlement, U, z, u, sFinal: sumMvH * qPrev };
}

/**
 * Biot poroelastic column under a suddenly applied surface load (Terzaghi's problem with compressible
 * constituents; Verruijt 2013 ch. 2). With confined drained modulus Kv = K + 4G/3, Biot coefficient α and
 * Biot modulus M (1/M = n/Kf + (α − n)/Ks):
 *   undrained response  p0 = q α M / (Kv + α² M);
 *   consolidation coefficient  c = (k/γw) M Kv / (Kv + α² M);
 *   p(z, t)/p0 = excessPorePressureRatio(z/h, c t / h²)  (drained top, impervious base at depth h);
 *   settlement  w(t) = q h / (Kv + α² M) + [q h/Kv − q h/(Kv + α² M)] U(Tv).
 * For incompressible grains and fluid (α = 1, M → ∞) this reduces to Terzaghi: p0 = q, c = cv = k Kv/γw.
 * @param {object} o
 * @param {number} o.q surface load [Pa]
 * @param {number} o.h column height (drainage path) [m]
 * @param {number} o.K drained bulk modulus [Pa]
 * @param {number} o.G shear modulus [Pa]
 * @param {number} o.k hydraulic conductivity [m/s]
 * @param {number} [o.alpha=1] Biot coefficient [-]
 * @param {number} [o.M=Infinity] Biot modulus [Pa]
 * @param {number} [o.gammaW=9810] [N/m^3]
 * @param {number} o.t time [s]
 * @param {number} [o.z] depth at which to return pore pressure [m] (default h, the base)
 * @returns {{p0:number,c:number,Tv:number,p:number,U:number,settlement:number,settlementImmediate:number,
 *   settlementFinal:number}} pressures [Pa]; c [m^2/s]; settlements [m]
 */
export function biotColumn({ q, h, K, G, k, alpha = 1, M = Infinity, gammaW = 9810, t, z }) {
  need('h', h, 1e-12);
  need('K', K, 1e-12);
  need('G', G, 0);
  need('k', k, 1e-30);
  need('alpha', alpha, 0, 1);
  need('t', t, 0);
  const Kv = K + (4 * G) / 3;
  const fin = Number.isFinite(M);
  const p0 = fin ? (q * alpha * M) / (Kv + alpha * alpha * M) : q / alpha;
  const c = fin ? ((k / gammaW) * M * Kv) / (Kv + alpha * alpha * M) : ((k / gammaW) * Kv) / (alpha * alpha);
  const Tv = (c * t) / (h * h);
  const zz = z === undefined ? h : need('z', z, 0, h);
  const s0 = fin ? (q * h) / (Kv + alpha * alpha * M) : 0;
  const sInf = (q * h) / Kv;
  const U = consolidationU(Tv);
  return { p0, c, Tv, p: p0 * excessPorePressureRatio(zz / h, Tv), U, settlement: s0 + (sInf - s0) * U, settlementImmediate: s0, settlementFinal: sInf };
}
