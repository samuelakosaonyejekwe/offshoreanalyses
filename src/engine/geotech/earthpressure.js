/**
 * Lateral earth pressure, retaining structures and sheet-pile walls.
 *
 * SI units: m, N (per metre run of wall: N/m), Pa, N/m^3, rad.
 * Angle conventions (all positive as drawn in the usual textbook figure): φ soil friction, δ wall friction,
 * β backfill slope above horizontal, θ wall back-face batter from vertical (positive when the face leans
 * towards the backfill, i.e. the retained wedge is larger).
 *
 * References
 *   Rankine W.J.M. (1857); Coulomb C.A. (1776); Jaky J. (1944); Mayne P.W. & Kulhawy F.H. (1982) K0–OCR.
 *   Okabe S. (1926); Mononobe N. & Matsuo H. (1929); Seed H.B. & Whitman R.V. (1970).
 *   Westergaard H.M. (1933) Water pressures on dams during earthquakes, Trans. ASCE 98.
 *   Terzaghi K. (1954) Anchored bulkheads, Trans. ASCE 119; Rowe P.W. (1952) Proc. ICE 1.
 *   Meyerhof G.G. (1953); BS 6349-2; EAU 2012; USS Steel Sheet Piling Design Manual (1984).
 */
import { brent } from '../math/roots.js';

const PI = Math.PI;

function need(name, v, min = -Infinity, max = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < min || v > max) {
    throw new RangeError(`${name} must be a number in [${min}, ${max}] (got ${v})`);
  }
  return v;
}

/**
 * Rankine coefficients for a vertical smooth wall and backfill sloping at β (pressure acts parallel to the
 * slope):  Ka,p = cos β [cos β ∓ √(cos²β − cos²φ)] / [cos β ± √(cos²β − cos²φ)].
 * For β = 0: Ka = tan²(π/4 − φ/2), Kp = tan²(π/4 + φ/2).
 * @param {number} phi friction angle [rad], 0 ≤ φ < 90°
 * @param {number} [beta=0] backfill slope [rad], 0 ≤ β ≤ φ
 * @returns {{Ka:number,Kp:number}} [-]
 */
export function rankine(phi, beta = 0) {
  need('phi', phi, 0, PI / 2 - 1e-9);
  need('beta', beta, 0, Math.max(phi, 0));
  const cb = Math.cos(beta);
  const r = Math.sqrt(Math.max(cb * cb - Math.cos(phi) ** 2, 0));
  return { Ka: (cb * (cb - r)) / (cb + r), Kp: (cb * (cb + r)) / (cb - r) };
}

/**
 * Coulomb active coefficient (resultant inclined at δ to the wall normal):
 *   Ka = cos²(φ − θ) / { cos²θ cos(δ + θ) [1 + √( sin(φ + δ) sin(φ − β) / (cos(δ + θ) cos(θ − β)) )]² }.
 * @param {object} o
 * @param {number} o.phi [rad]
 * @param {number} [o.delta=0] wall friction [rad], 0 ≤ δ ≤ φ
 * @param {number} [o.beta=0] backfill slope [rad], β ≤ φ
 * @param {number} [o.theta=0] wall batter from vertical [rad]
 * @returns {number} Ka [-]
 */
export function coulombKa({ phi, delta = 0, beta = 0, theta = 0 }) {
  need('phi', phi, 0, PI / 2 - 1e-9);
  need('delta', delta, 0, PI / 2 - 1e-9);
  need('beta', beta, -PI / 2, phi);
  need('theta', theta, -PI / 3, PI / 3);
  const root = Math.sqrt((Math.sin(phi + delta) * Math.sin(phi - beta)) / (Math.cos(delta + theta) * Math.cos(theta - beta)));
  return Math.cos(phi - theta) ** 2 / (Math.cos(theta) ** 2 * Math.cos(delta + theta) * (1 + root) ** 2);
}

/**
 * Coulomb passive coefficient (plane wedge; unconservative for δ > φ/3, where log-spiral values such as
 * Caquot & Kerisel should be used instead):
 *   Kp = cos²(φ + θ) / { cos²θ cos(δ − θ) [1 − √( sin(φ + δ) sin(φ + β) / (cos(δ − θ) cos(β − θ)) )]² }.
 * @param {object} o phi, delta, beta, theta [rad] as coulombKa
 * @returns {number} Kp [-]
 */
export function coulombKp({ phi, delta = 0, beta = 0, theta = 0 }) {
  need('phi', phi, 0, PI / 2 - 1e-9);
  need('delta', delta, 0, PI / 2 - 1e-9);
  need('theta', theta, -PI / 3, PI / 3);
  const root = Math.sqrt((Math.sin(phi + delta) * Math.sin(phi + beta)) / (Math.cos(delta - theta) * Math.cos(beta - theta)));
  if (!(root < 1)) throw new RangeError('Coulomb passive wedge is not defined for this combination of phi, delta and beta');
  return Math.cos(phi + theta) ** 2 / (Math.cos(theta) ** 2 * Math.cos(delta - theta) * (1 - root) ** 2);
}

/**
 * At-rest coefficient: Jaky (1944) K0 = 1 − sin φ′ for normally consolidated soil, extended for
 * over-consolidation by Mayne & Kulhawy (1982): K0 = (1 − sin φ′) OCR^(sin φ′).
 * @param {number} phi effective friction angle [rad]
 * @param {number} [ocr=1] over-consolidation ratio [-], ≥ 1
 * @returns {number} K0 [-]
 */
export function k0(phi, ocr = 1) {
  need('phi', phi, 0, PI / 2 - 1e-9);
  need('ocr', ocr, 1);
  const s = Math.sin(phi);
  return (1 - s) * ocr ** s;
}

/**
 * Depth of the tension crack behind a wall in c–φ soil (Rankine active): zc = (2c/√Ka − q) / γ ≥ 0.
 * @param {object} o c [Pa], gamma [N/m^3], Ka [-], q surcharge [Pa]
 * @returns {number} zc [m]
 */
export function tensionCrackDepth({ c, gamma, Ka, q = 0 }) {
  need('c', c, 0);
  need('gamma', gamma, 1e-9);
  need('Ka', Ka, 1e-9);
  return Math.max(((2 * c) / Math.sqrt(Ka) - q) / gamma, 0);
}

/**
 * Lateral pressure diagram on a vertical wall in layered soil with a water table and uniform surcharge:
 *   σ′h = K σ′v − 2c√K (active),  K σ′v + 2c√K (passive),  K0 σ′v (rest);  u = γw (z − zw).
 * K defaults to the Rankine value of each layer (K0 from Jaky for 'rest'); a layer may carry its own K
 * (e.g. Coulomb with wall friction, horizontal component). Tensile active stresses are discarded when
 * `tensionCrack` is true. The resultant and its line of action follow from exact integration of the
 * piecewise-linear diagram.
 * @param {object} o
 * @param {{thickness:number,gamma:number,gammaSat?:number,phi:number,c?:number,K?:number}[]} o.layers
 *   top to bottom; gamma bulk unit weight above the water table, gammaSat below [N/m^3]; φ [rad]; c [Pa]
 * @param {'active'|'passive'|'rest'} [o.state='active']
 * @param {number} [o.waterTableDepth=Infinity] depth of the water table below the top [m]
 * @param {number} [o.surcharge=0] uniform surface surcharge [Pa]
 * @param {number} [o.gammaW=9810] unit weight of water [N/m^3]
 * @param {boolean} [o.tensionCrack=true]
 * @returns {{points:{z:number,sigmaV:number,u:number,sigmaH:number}[],Peff:number,Pwater:number,P:number,
 *   zBar:number,armAboveBase:number,H:number}}
 *   points: effective vertical and horizontal stress and pore pressure [Pa] at each break (two entries at
 *   a layer interface); Peff, Pwater, P = total thrust [N/m]; zBar depth of the total resultant [m];
 *   armAboveBase [m]; H total height [m]
 */
export function pressureDiagram({ layers, state = 'active', waterTableDepth = Infinity, surcharge = 0, gammaW = 9810, tensionCrack = true }) {
  if (!Array.isArray(layers) || layers.length === 0) throw new RangeError('layers must be a non-empty array');
  need('surcharge', surcharge, 0);
  need('waterTableDepth', waterTableDepth, 0);
  const points = [];
  let z = 0;
  let sv = surcharge;
  let Peff = 0;
  let Meff = 0;
  const addSegment = (z1, z2, p1, p2) => {
    // integrate the (optionally clipped) linear segment p(z) between z1 and z2
    if (tensionCrack && state === 'active') {
      if (p1 <= 0 && p2 <= 0) return;
      if (p1 < 0) { z1 += ((z2 - z1) * -p1) / (p2 - p1); p1 = 0; } else if (p2 < 0) { z2 = z1 + ((z2 - z1) * p1) / (p1 - p2); p2 = 0; }
    }
    const h = z2 - z1;
    const F = 0.5 * (p1 + p2) * h;
    Peff += F;
    Meff += h * (p1 * (z1 + h / 3) + p2 * (z1 + (2 * h) / 3)) / 2;
  };
  for (const ly of layers) {
    need('layer.thickness', ly.thickness, 1e-12);
    need('layer.gamma', ly.gamma, 0);
    const c = ly.c || 0;
    let K = ly.K;
    if (K === undefined) K = state === 'rest' ? k0(ly.phi) : state === 'active' ? rankine(ly.phi).Ka : rankine(ly.phi).Kp;
    const sh = (s) => (state === 'active' ? K * s - 2 * c * Math.sqrt(K) : state === 'passive' ? K * s + 2 * c * Math.sqrt(K) : K * s);
    const zEnd = z + ly.thickness;
    const breaks = [z];
    if (waterTableDepth > z && waterTableDepth < zEnd) breaks.push(waterTableDepth);
    breaks.push(zEnd);
    for (let i = 0; i < breaks.length - 1; i++) {
      const za = breaks[i];
      const zb = breaks[i + 1];
      const wet = za >= waterTableDepth;
      const gEff = wet ? (ly.gammaSat !== undefined ? ly.gammaSat : ly.gamma) - gammaW : ly.gamma;
      const svB = sv + gEff * (zb - za);
      const pa = sh(sv);
      const pb = sh(svB);
      points.push({ z: za, sigmaV: sv, u: Math.max(za - waterTableDepth, 0) * gammaW, sigmaH: pa });
      points.push({ z: zb, sigmaV: svB, u: Math.max(zb - waterTableDepth, 0) * gammaW, sigmaH: pb });
      addSegment(za, zb, pa, pb);
      sv = svB;
    }
    z = zEnd;
  }
  const hw = Math.max(z - waterTableDepth, 0);
  const Pwater = 0.5 * gammaW * hw * hw;
  const Mwater = Pwater * (z - hw / 3);
  const P = Peff + Pwater;
  const zBar = P !== 0 ? (Meff + Mwater) / P : 0;
  return { points, Peff, Pwater, P, zBar, armAboveBase: z - zBar, H: z };
}

/**
 * Horizontal stress on a rigid wall from a strip surcharge of intensity q parallel to the wall (elastic
 * solution doubled for the rigid, non-yielding wall; Terzaghi 1954):
 *   σh = (2q/π) (β − sin β cos 2α),
 * β = angle subtended by the strip at the point, α = angle from the wall to the bisector.
 * @param {object} o
 * @param {number} o.q strip pressure [Pa]
 * @param {number} o.a distance from the wall to the near edge of the strip [m]
 * @param {number} o.b strip width [m]
 * @param {number} o.z depth below the surface [m]
 * @returns {number} σh [Pa]
 */
export function stripSurchargePressure({ q, a, b, z }) {
  need('a', a, 0);
  need('b', b, 1e-12);
  need('z', z, 0);
  if (z === 0) return 0;
  const t1 = Math.atan(a / z);
  const t2 = Math.atan((a + b) / z);
  const beta = t2 - t1;
  const alpha = 0.5 * (t1 + t2);
  return ((2 * q) / PI) * (beta - Math.sin(beta) * Math.cos(2 * alpha));
}

/**
 * Horizontal stress on a rigid wall of height H from a line load Q (N/m) at distance mH, at depth nH
 * (Boussinesq doubled, with Terzaghi's (1954) empirical modification for loads close to the wall):
 *   m > 0.4:  σh = (4Q/πH) m² n / (m² + n²)²;    m ≤ 0.4:  σh = (Q/H) 0.203 n / (0.16 + n²)².
 * @param {object} o Q [N/m], x distance from the wall [m], z depth [m], H wall height [m]
 * @returns {number} σh [Pa]
 */
export function lineSurchargePressure({ Q, x, z, H }) {
  need('x', x, 0);
  need('H', H, 1e-12);
  need('z', z, 0, H);
  const m = x / H;
  const n = z / H;
  if (m > 0.4) return ((4 * Q) / (PI * H)) * ((m * m * n) / (m * m + n * n) ** 2);
  return (Q / H) * ((0.203 * n) / (0.16 + n * n) ** 2);
}

/**
 * Mononobe–Okabe seismic earth-pressure coefficients (cohesionless, dry backfill):
 *   ψ = atan(kh / (1 − kv)),
 *   KAE = cos²(φ − θ − ψ) / { cos ψ cos²θ cos(δ + θ + ψ) [1 + √( sin(φ + δ) sin(φ − β − ψ) / (cos(δ + θ + ψ) cos(β − θ)) )]² }
 *   KPE = cos²(φ + θ − ψ) / { cos ψ cos²θ cos(δ − θ + ψ) [1 − √( sin(φ + δ) sin(φ + β − ψ) / (cos(δ − θ + ψ) cos(β − θ)) )]² }
 * Total thrust PAE = ½ γ H² (1 − kv) KAE. The active solution exists only for φ − β − ψ ≥ 0.
 * Reduces to the Coulomb coefficients at kh = kv = 0.
 * @param {object} o
 * @param {number} o.phi [rad]
 * @param {number} [o.delta=0] [rad]
 * @param {number} [o.beta=0] [rad]
 * @param {number} [o.theta=0] [rad]
 * @param {number} o.kh horizontal seismic coefficient [-]
 * @param {number} [o.kv=0] vertical seismic coefficient (positive downward inertia reduction) [-]
 * @returns {{KAE:number,KPE:number,psi:number}} [-,-,rad]
 */
export function mononobeOkabe({ phi, delta = 0, beta = 0, theta = 0, kh, kv = 0 }) {
  need('phi', phi, 0, PI / 2 - 1e-9);
  need('delta', delta, 0, PI / 2 - 1e-9);
  need('kh', kh, 0, 1);
  need('kv', kv, -0.9, 0.9);
  const psi = Math.atan(kh / (1 - kv));
  if (phi - beta - psi < 0) throw new RangeError('kh is too large: phi − beta − psi < 0, no Mononobe–Okabe active solution');
  const ra = Math.sqrt((Math.sin(phi + delta) * Math.sin(phi - beta - psi)) / (Math.cos(delta + theta + psi) * Math.cos(beta - theta)));
  const KAE = Math.cos(phi - theta - psi) ** 2 / (Math.cos(psi) * Math.cos(theta) ** 2 * Math.cos(delta + theta + psi) * (1 + ra) ** 2);
  const sp = Math.sin(phi + beta - psi);
  let KPE = NaN;
  if (sp >= 0) {
    const rp = Math.sqrt((Math.sin(phi + delta) * sp) / (Math.cos(delta - theta + psi) * Math.cos(beta - theta)));
    if (rp < 1) KPE = Math.cos(phi + theta - psi) ** 2 / (Math.cos(psi) * Math.cos(theta) ** 2 * Math.cos(delta - theta + psi) * (1 - rp) ** 2);
  }
  return { KAE, KPE, psi };
}

/**
 * Westergaard (1933) hydrodynamic water pressure on a rigid vertical face during an earthquake
 * (incompressible approximation):  p(z) = (7/8) kh γw √(h z),  total force (7/12) kh γw h² acting 0.4 h
 * above the base.
 * @param {object} o
 * @param {number} o.kh horizontal seismic coefficient [-]
 * @param {number} o.h water depth [m]
 * @param {number} [o.z] depth below the water surface for the local pressure [m] (default h)
 * @param {number} [o.gammaW=9810] [N/m^3]
 * @returns {{p:number,force:number,armAboveBase:number}} [Pa, N/m, m]
 */
export function westergaard({ kh, h, z, gammaW = 9810 }) {
  need('kh', kh, 0, 1);
  need('h', h, 0);
  const zz = z === undefined ? h : need('z', z, 0, h);
  return { p: (7 / 8) * kh * gammaW * Math.sqrt(h * zz), force: (7 / 12) * kh * gammaW * h * h, armAboveBase: 0.4 * h };
}

/**
 * Stability of a gravity quay wall / caisson per metre run, from statics about the toe.
 *   sliding      FS = [(ΣV − U) tan δb + cb B + Pp] / ΣH
 *   overturning  FS = ΣM_resisting / ΣM_overturning (uplift counted as overturning)
 *   bearing      e = B/2 − (ΣMr − ΣMo)/N, N = ΣV − U; Meyerhof q = N / (B − 2e); trapezoidal
 *                qmax,min = N/B (1 ± 6e/B) or triangular 2N / (3 (B/2 − e)) when e > B/6.
 * A global (deep-seated) factor of safety can be attached through `globalFs`, e.g. from slope.js.
 * @param {object} o
 * @param {number} o.B base width [m]
 * @param {{V:number,x:number}[]} o.vertical vertical loads [N/m] (down positive) with arms from the toe [m]
 * @param {{H:number,y:number}[]} o.horizontal horizontal loads [N/m] (towards the toe positive) with heights
 *   above the base [m]
 * @param {number} [o.U=0] uplift resultant on the base [N/m]
 * @param {number} [o.xU=o.B/2] arm of the uplift from the toe [m]
 * @param {number} o.deltaBase base friction angle [rad]
 * @param {number} [o.cBase=0] base adhesion [Pa]
 * @param {number} [o.Pp=0] passive resistance at the toe counted for sliding [N/m]
 * @param {number} [o.qult] ultimate bearing pressure on the effective width [Pa]
 * @param {number} [o.globalFs] externally computed global factor of safety [-]
 * @returns {{fsSliding:number,fsOverturning:number,e:number,N:number,H:number,qMeyerhof:number,qMax:number,
 *   qMin:number,fsBearing:number,inMiddleThird:boolean,globalFs:number|undefined}} forces [N/m], q [Pa]
 */
export function gravityWallStability({ B, vertical, horizontal, U = 0, xU, deltaBase, cBase = 0, Pp = 0, qult, globalFs }) {
  need('B', B, 1e-9);
  need('deltaBase', deltaBase, 0, PI / 2 - 1e-9);
  let V = 0;
  let Mr = 0;
  let Mo = 0;
  let H = 0;
  for (const v of vertical) { V += v.V; if (v.V >= 0) Mr += v.V * v.x; else Mo += -v.V * v.x; }
  for (const h of horizontal) { H += h.H; if (h.H >= 0) Mo += h.H * h.y; else Mr += -h.H * h.y; }
  Mo += U * (xU === undefined ? B / 2 : xU);
  const N = V - U;
  if (!(N > 0)) throw new RangeError('net vertical load on the base must be positive');
  const e = B / 2 - (Mr - Mo) / N;
  const Rs = N * Math.tan(deltaBase) + cBase * B + Pp;
  const eAbs = Math.abs(e);
  const within = eAbs <= B / 6;
  const qMax = eAbs >= B / 2 ? Infinity : within ? (N / B) * (1 + (6 * eAbs) / B) : (2 * N) / (3 * (B / 2 - eAbs));
  const qMin = within ? (N / B) * (1 - (6 * eAbs) / B) : 0;
  const qM = eAbs >= B / 2 ? Infinity : N / (B - 2 * eAbs);
  return {
    fsSliding: H > 0 ? Rs / H : Infinity, fsOverturning: Mo > 0 ? Mr / Mo : Infinity, e, N, H,
    qMeyerhof: qM, qMax, qMin, fsBearing: qult !== undefined ? qult / qM : NaN, inMiddleThird: within, globalFs,
  };
}

/**
 * Effect of a bollard (mooring line) pull on a quay wall, per metre run: the horizontal component is
 * distributed over an effective wall length and applied at cope level; the vertical component relieves the
 * wall weight.
 * @param {object} o
 * @param {number} o.pull line pull [N]
 * @param {number} [o.angle=0] line angle above horizontal [rad]
 * @param {number} o.spreadLength wall length sharing one bollard (monolith length or bollard spacing) [m]
 * @param {number} o.height height of the bollard above the wall base [m]
 * @returns {{H:number,Vup:number,M:number}} horizontal [N/m], uplift [N/m], overturning moment about the base [N·m/m]
 */
export function bollardPull({ pull, angle = 0, spreadLength, height }) {
  need('pull', pull, 0);
  need('spreadLength', spreadLength, 1e-9);
  need('height', height, 0);
  const H = (pull * Math.cos(angle)) / spreadLength;
  return { H, Vup: (pull * Math.sin(angle)) / spreadLength, M: H * height };
}

function wallPressures({ H, gamma, gammaSub, waterDepth, Ka, Kp, surcharge }) {
  // effective vertical stress at depth z below a surface located at depth z0 (z ≥ z0)
  const sv = (z, z0) => {
    if (z <= z0) return 0;
    const dry = Math.max(Math.min(z, waterDepth) - z0, 0);
    const wet = z - Math.max(z0, waterDepth);
    return gamma * dry + gammaSub * Math.max(wet, 0);
  };
  return {
    active: (z) => Ka * (surcharge + sv(z, 0)),
    passive: (z) => (z > H ? Kp * sv(z, H) : 0),
  };
}

/**
 * Sheet-pile wall in cohesionless soil by limit equilibrium (free-earth support), per metre run.
 *   Cantilever (no anchor): simplified method – embedment d0 from moment equilibrium of active and
 *   (factored) passive pressure about the toe, then increased by `embedmentIncrease` (20 % by convention)
 *   to develop the toe reaction.
 *   Anchored: free-earth support – embedment from moment equilibrium about the anchor level, anchor force
 *   from horizontal equilibrium, maximum bending moment at the point of zero shear. Rowe's (1952)
 *   moment-reduction factor for wall flexibility is applied when supplied as `roweFactor` (chart value).
 * Water is taken at the same level on both sides (no net water pressure). The passive coefficient is
 * divided by `fsPassive`.
 * @param {object} o
 * @param {number} o.H retained height (top to dredge level) [m]
 * @param {number} o.gamma bulk unit weight above water [N/m^3]
 * @param {number} [o.gammaSub=o.gamma] submerged unit weight below water [N/m^3]
 * @param {number} [o.waterDepth=Infinity] water level below the top of the wall, both sides [m]
 * @param {number} o.Ka active coefficient (horizontal component) [-]
 * @param {number} o.Kp passive coefficient (horizontal component) [-]
 * @param {number} [o.fsPassive=1] factor of safety on Kp [-]
 * @param {number} [o.surcharge=0] uniform surcharge on the retained side [Pa]
 * @param {number} [o.anchorDepth] depth of the anchor below the top [m]; omit for a cantilever wall
 * @param {number} [o.embedmentIncrease=1.2] multiplier on d0 for cantilever walls [-]
 * @param {number} [o.roweFactor=1] moment reduction factor [-], ≤ 1
 * @param {number} [o.n=4000] integration intervals
 * @returns {{d0:number,d:number,length:number,anchorForce:number,Mmax:number,zMmax:number,Pa:number,Pp:number}}
 *   d0 equilibrium embedment [m]; d design embedment [m]; total length [m]; anchor force [N/m]
 *   (0 for cantilever); Mmax maximum bending moment [N·m/m] at depth zMmax [m]; thrusts at d0 [N/m]
 */
export function sheetPileWall(o) {
  const { H, gamma, waterDepth = Infinity, Ka, Kp, fsPassive = 1, surcharge = 0, anchorDepth, embedmentIncrease = 1.2, roweFactor = 1, n = 4000 } = o;
  need('H', H, 1e-9);
  need('gamma', gamma, 1e-9);
  need('Ka', Ka, 1e-9);
  need('Kp', Kp, 1e-9);
  need('fsPassive', fsPassive, 1e-9);
  need('surcharge', surcharge, 0);
  need('roweFactor', roweFactor, 0.1, 1);
  const gammaSub = o.gammaSub === undefined ? gamma : need('gammaSub', o.gammaSub, 1e-9);
  const kp = Kp / fsPassive;
  if (kp <= Ka) throw new RangeError('Kp/fsPassive must exceed Ka for equilibrium to exist');
  const anchored = anchorDepth !== undefined;
  if (anchored) need('anchorDepth', anchorDepth, 0, H);
  const pr = wallPressures({ H, gamma, gammaSub, waterDepth, Ka, Kp: kp, surcharge });
  const net = (z) => pr.active(z) - pr.passive(z);
  // Simpson on each smooth range (breaks at the water table and the dredge level)
  const integ = (f, a, b, m = 200) => {
    if (b <= a) return 0;
    const h = (b - a) / m;
    let s = 0;
    for (let i = 0; i <= m; i++) s += (i === 0 || i === m ? 1 : i % 2 ? 4 : 2) * f(a + i * h);
    return (s * h) / 3;
  };
  const breaks = (Lw) => [0, ...[waterDepth, H].filter((x) => x > 0 && x < Lw).sort((a, b) => a - b), Lw];
  const momentExact = (d) => {
    const Lw = H + d;
    const pivot = anchored ? anchorDepth : Lw;
    const bs = breaks(Lw);
    let s = 0;
    for (let i = 0; i < bs.length - 1; i++) s += integ((z) => net(z) * (anchored ? z - pivot : pivot - z), bs[i], bs[i + 1]);
    return s;
  };
  let hi = H;
  let guard = 0;
  while (momentExact(hi) > 0 && guard++ < 60) hi *= 1.6;
  if (momentExact(hi) > 0) throw new Error('sheet-pile embedment did not converge: passive resistance insufficient');
  const d0 = brent(momentExact, 1e-9, hi, { tol: 1e-11 });
  const Lw0 = H + d0;
  const bs = breaks(Lw0);
  let Pa = 0;
  let Pp = 0;
  for (let i = 0; i < bs.length - 1; i++) { Pa += integ(pr.active, bs[i], bs[i + 1]); Pp += integ(pr.passive, bs[i], bs[i + 1]); }
  const T = anchored ? Pa - Pp : 0;
  // shear and moment by marching down the wall; the anchor contributes −T to the shear and −T (z − a) to
  // the moment below its level, added in closed form so that the jump is not smeared over a step
  const dz = Lw0 / n;
  let S = 0;
  let Ms = 0;
  let Mmax = 0;
  let zM = 0;
  let pPrev = net(0);
  for (let i = 1; i <= n; i++) {
    const z = i * dz;
    const p = net(z);
    const Sn = S + 0.5 * (pPrev + p) * dz;
    // Simpson-consistent moment increment for a shear varying quadratically within the step
    Ms += 0.5 * (S + Sn) * dz - ((p - pPrev) * dz * dz) / 12;
    S = Sn;
    pPrev = p;
    const M = Ms - (anchored && z > anchorDepth ? T * (z - anchorDepth) : 0);
    if (Math.abs(M) > Math.abs(Mmax)) { Mmax = M; zM = z; }
  }
  const d = anchored ? d0 : embedmentIncrease * d0;
  return { d0, d, length: H + d, anchorForce: T, Mmax: Math.abs(Mmax) * roweFactor, zMmax: zM, Pa, Pp };
}

/**
 * Ultimate resistance of a deadman anchor (anchor wall or block) in cohesionless soil near the surface,
 * per metre run for a continuous wall: Tult = ½ γ H² (Kp − Ka), valid when the wall height h is at least
 * about half the depth H to its base (Terzaghi 1943; Teng 1962). For a short block of width b the
 * three-dimensional end effect is applied as a multiplier `shapeFactor` on the net pressure (argument;
 * Ovesen's or Teng's charts), giving Tult = shapeFactor · ½ γ H² (Kp − Ka) · b.
 * The deadman must lie outside the active wedge of the wall for the full value to apply.
 * @param {object} o
 * @param {number} o.gamma effective unit weight [N/m^3]
 * @param {number} o.H depth to the base of the deadman [m]
 * @param {number} o.Ka [-]
 * @param {number} o.Kp [-]
 * @param {number} [o.b] block width [m]; omit for a continuous wall (result per metre)
 * @param {number} [o.shapeFactor=1] 3-D enhancement factor [-]
 * @param {number} [o.safetyFactor=1] [-]
 * @returns {{Tult:number,Tallow:number}} [N/m] for a continuous wall, [N] for a block
 */
export function deadmanAnchor({ gamma, H, Ka, Kp, b, shapeFactor = 1, safetyFactor = 1 }) {
  need('gamma', gamma, 0);
  need('H', H, 0);
  need('Ka', Ka, 0);
  need('Kp', Kp, Ka);
  need('safetyFactor', safetyFactor, 1e-9);
  const perM = 0.5 * gamma * H * H * (Kp - Ka);
  const Tult = b === undefined ? perM : shapeFactor * perM * need('b', b, 0);
  return { Tult, Tallow: Tult / safetyFactor };
}
