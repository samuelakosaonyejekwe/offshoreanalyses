/**
 * Subsea mudmat foundations: combined V–H–M–T capacity, skirt penetration, contact pressure, landing,
 * breakout and on-bottom stability.
 *
 * SI units: m, N, Pa, N/m^3, rad, kg, s.
 *
 * References
 *   API RP 2GEO (2011) 7 and Annex A; ISO 19901-4:2016 7; DNV-RP-C212 (2019) 5–6.
 *   Green A.P. (1954) The plastic yielding of metal junctions due to combined shear and pressure,
 *     J. Mech. Phys. Solids 2 (exact V–H interaction of a strip on Tresca material).
 *   Meyerhof G.G. (1953) The bearing capacity of foundations under eccentric and inclined loads.
 *   Bransby M.F. & Randolph M.F. (1998) Geotechnique 48(5); Gourvenec S. (2007) Geotechnique 57(9);
 *   Feng X., Randolph M.F., Gourvenec S. & Wallerand R. (2014) Geotechnique 64(1).
 *   DNV-RP-H103 / DNV-RP-N103 6 (landing on the seabed, retrieval).
 */
import { brent } from '../math/roots.js';
import { bearingCapacity, effectiveAreaRect } from './bearing.js';

const PI = Math.PI;

function need(name, v, min = -Infinity, max = Infinity) {
  if (typeof v !== 'number' || Number.isNaN(v) || v < min || v > max) {
    throw new RangeError(`${name} must be a number in [${min}, ${max}] (got ${v})`);
  }
  return v;
}

/**
 * Exact V–H interaction of a rough strip on uniform Tresca soil (Green 1954):
 *   V/(B su) = 1 + π/2 + acos(h) + √(1 − h²),  h = H/(B su),  for V/Vult ≥ 0.5;
 *   H = B su (pure sliding) for V/Vult < 0.5.
 * @param {number} h normalised horizontal load H/(A su) [-], 0 ≤ h ≤ 1
 * @returns {number} largest admissible V/(A su) [-] (2 + π at h = 0, 1 + π/2 at h = 1)
 */
export function greenVH(h) {
  need('h', h, 0, 1);
  return 1 + PI / 2 + Math.acos(h) + Math.sqrt(1 - h * h);
}

/**
 * Inverse of greenVH: limiting H/(A su) for a given vertical mobilisation v = V/Vult.
 * @param {number} v V/Vult [-], 0 ≤ v ≤ 1
 * @returns {number} h = H/(A su) [-]
 */
export function greenHofV(v) {
  need('v', v, 0, 1);
  if (v <= 0.5) return 1;
  if (v >= 1) return 0;
  const target = v * (2 + PI);
  let lo = 0;
  let hi = 1;
  for (let i = 0; i < 60; i++) {
    const mid = 0.5 * (lo + hi);
    if (greenVH(mid) > target) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}

/**
 * Plastic torsional capacity of a rectangular base in full contact on undrained soil:
 *   Tult = su ∫∫ r dA  (shear su mobilised tangentially everywhere), with the closed form
 *   ∫∫ r dA = (2/3) [ 2 a b d + a³ ln((b + d)/a) + b³ ln((a + d)/b) ],  a = B/2, b = L/2, d = √(a² + b²).
 * @param {object} o
 * @param {number} o.B width [m]
 * @param {number} o.L length [m]
 * @param {number} o.su interface undrained strength [Pa]
 * @returns {number} Tult [N·m]
 */
export function torsionalCapacityRect({ B, L, su }) {
  need('B', B, 1e-9);
  need('L', L, 1e-9);
  need('su', su, 0);
  const a = B / 2;
  const b = L / 2;
  const d = Math.hypot(a, b);
  const quarter = (2 * a * b * d + a ** 3 * Math.log((b + d) / a) + b ** 3 * Math.log((a + d) / b)) / 6;
  return 4 * su * quarter;
}

/**
 * Horizontal load equivalent to combined horizontal load and torsion on the effective area
 * (DNV-RP-C212 5.4 / DNV Classification Notes 30.4):  H′ = 2T/L′ + √(H² + (2T/L′)²).
 * @param {object} o
 * @param {number} o.H resultant horizontal load [N]
 * @param {number} o.T torsional moment about the vertical axis [N·m]
 * @param {number} o.Leff effective length L′ [m]
 * @returns {number} H′ [N]
 */
export function torsionEquivalentH({ H, T, Leff }) {
  need('H', H, 0);
  need('Leff', Leff, 1e-9);
  const a = (2 * Math.abs(T)) / Leff;
  return a + Math.hypot(H, a);
}

/**
 * Classical combined-loading check of a rectangular mudmat: Meyerhof effective area for the moments,
 * equivalent horizontal load for torsion, Brinch Hansen or Vesic factors for inclination and shape
 * (bearing.js), and a separate sliding check.
 * Undrained: uniform su (use the strength at about B′/4 below skirt-tip level for a moderately increasing
 * profile, or isoUndrainedCapacity in bearing.js with a Davis & Booker factor). Drained: c′–φ′.
 * @param {object} o
 * @param {number} o.B width [m]
 * @param {number} o.L length [m]
 * @param {number} o.V vertical load incl. submerged weight [N]
 * @param {number} [o.HB=0] horizontal load along B [N]
 * @param {number} [o.HL=0] horizontal load along L [N]
 * @param {number} [o.MB=0] moment producing eccentricity along B [N·m]
 * @param {number} [o.ML=0] moment producing eccentricity along L [N·m]
 * @param {number} [o.T=0] torsion [N·m]
 * @param {number} [o.su] undrained strength [Pa] (undrained case, phi = 0)
 * @param {number} [o.phi=0] friction angle [rad]
 * @param {number} [o.c=0] effective cohesion [Pa]
 * @param {number} o.gamma effective unit weight [N/m^3]
 * @param {number} [o.D=0] skirt depth (base taken at skirt tip) [m]
 * @param {number} [o.delta] interface friction angle for drained sliding [rad] (default φ)
 * @param {'hansen'|'vesic'} [o.method='hansen']
 * @param {number} [o.materialFactor=1] γm on strength-derived capacities [-]
 * @returns {{Vult:number,ucBearing:number,Hres:number,Hequiv:number,Hult:number,ucSliding:number,
 *   Beff:number,Leff:number,Aeff:number,eB:number,eL:number,qult:number}} forces [N], lengths [m], uc [-]
 */
export function mudmatClassicalCapacity(o) {
  const { B, L, V, HB = 0, HL = 0, MB = 0, ML = 0, T = 0, phi = 0, c = 0, gamma, D = 0, method = 'hansen', materialFactor = 1 } = o;
  need('B', B, 1e-9);
  need('L', L, 1e-9);
  need('V', V, 1e-9);
  need('materialFactor', materialFactor, 1e-9);
  const und = phi < 1e-9;
  const coh = und ? need('su', o.su, 1e-9) : c;
  const eB = Math.abs(MB) / V;
  const eL = Math.abs(ML) / V;
  if (eB >= B / 2 || eL >= L / 2) throw new RangeError('load eccentricity lies outside the mudmat (MB, ML too large for V)');
  const ea = effectiveAreaRect({ B, L, eB, eL });
  const Hres = Math.hypot(HB, HL);
  const Hequiv = T !== 0 ? torsionEquivalentH({ H: Hres, T, Leff: ea.Leff }) : Hres;
  const scale = Hres > 0 ? Hequiv / Hres : 0;
  // Torsion without a horizontal load is applied along the short effective side (most onerous for iγ, ic).
  const hb = Hres > 0 ? HB * scale : Hequiv;
  const hl = Hres > 0 ? HL * scale : 0;
  const Hult = (und ? coh * ea.Aeff : c * ea.Aeff + V * Math.tan(o.delta !== undefined ? o.delta : phi)) / materialFactor;
  const cap = bearingCapacity({ method, B, L, D, c: coh, phi, gamma, V, HB: hb, HL: hl, eB, eL, q: gamma * D });
  // Net capacity: overburden at skirt-tip level is balanced by the soil plug weight.
  const Vult = ((cap.qult - gamma * D) * cap.Aeff) / materialFactor + gamma * D * cap.Aeff;
  return {
    Vult, ucBearing: V / Vult, Hres, Hequiv, Hult, ucSliding: Hequiv / Hult,
    Beff: cap.Beff, Leff: cap.Leff, Aeff: cap.Aeff, eB, eL, qult: cap.qult,
  };
}

/**
 * Failure-envelope check for a rectangular mat on undrained clay.
 *
 * The polynomial envelopes of Gourvenec (2007) and Feng et al. (2014) are not reproduced. This function
 * implements a generic, documented interaction
 *     (h / h*)^a + (m / m*)^b + (t / t*)^c = 1
 * in which the uniaxial limits shrink with vertical mobilisation v = V/Vult:
 *     h* = greenHofV(v)                  – exact for a strip on uniform Tresca soil (Green 1954);
 *     m* = 4 v (1 − v)   (no tension)    – exact Meyerhof effective-width result, peak at v = 0.5;
 *     m* = 1 − v^p       (full tension, skirted; p is an argument);
 *     t* = h*                            – torsion is a pure base-shear mechanism.
 * h = H/Hult, m = M/Mult, t = T/Tult are normalised by the uniaxial capacities supplied by the caller
 * (from bearing.js, torsionalCapacityRect, or published factors: e.g. Mult ≈ 0.69 A B su for a rough strip
 * with full tension, Bransby & Randolph 1998; (2 + π) A B su / 8 without tension).
 * Exponents default to a = b = c = 2 – a neutral choice, not a published fit: calibrate against the
 * project's finite-element or published envelope before relying on it.
 * @param {object} o
 * @param {number} o.V [N]
 * @param {number} [o.H=0] resultant horizontal load [N]
 * @param {number} [o.M=0] resultant moment [N·m]
 * @param {number} [o.T=0] torsion [N·m]
 * @param {number} o.Vult uniaxial vertical capacity [N]
 * @param {number} o.Hult uniaxial horizontal capacity [N]
 * @param {number} o.Mult uniaxial (peak) moment capacity [N·m]
 * @param {number} [o.Tult=Infinity] uniaxial torsional capacity [N·m]
 * @param {boolean} [o.tension=false] true for skirted mats able to sustain transient suction
 * @param {number} [o.a=2]
 * @param {number} [o.b=2]
 * @param {number} [o.c=2]
 * @param {number} [o.p=2] exponent of the full-tension moment reduction [-]
 * @returns {{f:number,inside:boolean,v:number,hStar:number,mStar:number,loadFactor:number}}
 *   f envelope value (≤ 1 inside) [-]; loadFactor: factor on (H, M, T) at constant V that reaches the envelope
 */
export function mudmatEnvelope(o) {
  const { V, H = 0, M = 0, T = 0, Vult, Hult, Mult, Tult = Infinity, tension = false, a = 2, b = 2, c = 2, p = 2 } = o;
  need('V', V, 0);
  need('Vult', Vult, 1e-12);
  need('Hult', Hult, 1e-12);
  need('Mult', Mult, 1e-12);
  need('a', a, 0.5, 10);
  need('b', b, 0.5, 10);
  need('c', c, 0.5, 10);
  const v = V / Vult;
  if (v >= 1) return { f: Infinity, inside: false, v, hStar: 0, mStar: 0, loadFactor: 0 };
  const hStar = greenHofV(v);
  const mStar = tension ? 1 - v ** p : 4 * v * (1 - v);
  const h = Math.abs(H) / Hult;
  const m = Math.abs(M) / Mult;
  const t = Number.isFinite(Tult) ? Math.abs(T) / Tult : 0;
  const val = (s) => (h > 0 ? ((s * h) / hStar) ** a : 0) + (m > 0 ? (mStar > 0 ? ((s * m) / mStar) ** b : Infinity) : 0) + (t > 0 ? ((s * t) / hStar) ** c : 0);
  const f = val(1);
  let loadFactor = Infinity;
  if (f > 0 && Number.isFinite(f)) {
    let lo = 0;
    let hi = 1;
    while (val(hi) < 1 && hi < 1e9) hi *= 2;
    for (let i = 0; i < 80; i++) {
      const mid = 0.5 * (lo + hi);
      if (val(mid) < 1) lo = mid; else hi = mid;
    }
    loadFactor = 0.5 * (lo + hi);
  } else if (!Number.isFinite(f)) loadFactor = 0;
  return { f, inside: f <= 1, v, hStar, mStar, loadFactor };
}

/**
 * Penetration resistance of mudmat skirts in clay with su = su0 + k z:
 *   R(z) = α su,avg(z) · As(z) + (Nc su(z) + γ′ z) · Atip,
 *   As = z · wallLength · 2 (both faces),  Atip = wallLength · t.
 * Nc = 7.5 is the deep strip value commonly used for skirt tips (argument).
 * @param {object} o
 * @param {number} o.z tip penetration [m]
 * @param {number} o.wallLength total plan length of skirts incl. internal ribs [m]
 * @param {number} o.t skirt wall thickness [m]
 * @param {number} o.su0 mudline strength [Pa]
 * @param {number} [o.k=0] strength gradient [Pa/m]
 * @param {number} [o.alpha=1] adhesion factor during penetration (≈ 1/sensitivity) [-]
 * @param {number} [o.Nc=7.5] tip bearing factor [-]
 * @param {number} [o.gamma=0] effective unit weight [N/m^3]
 * @returns {{R:number,friction:number,tip:number}} [N]
 */
export function skirtPenetrationClay({ z, wallLength, t, su0, k = 0, alpha = 1, Nc = 7.5, gamma = 0 }) {
  need('z', z, 0);
  need('wallLength', wallLength, 0);
  need('t', t, 0);
  need('su0', su0, 0);
  need('alpha', alpha, 0, 1.5);
  const friction = alpha * (su0 + 0.5 * k * z) * z * wallLength * 2;
  const tip = (Nc * (su0 + k * z) + gamma * z) * wallLength * t;
  return { R: friction + tip, friction, tip };
}

/**
 * CPT-based skirt penetration resistance (DNV-RP-C212 6.2 / DNV CN 30.4 4.5):
 *   R(z) = kp(z) Ap qc(z) + As′ ∫0^z kf qc dz,
 * with Ap the tip area and As′ the wall area per unit depth (both faces). Empirical coefficients quoted by
 * DNV for North Sea conditions: clay kp = 0.4 / 0.6, kf = 0.03 / 0.05; sand kp = 0.3 / 0.6,
 * kf = 0.001 / 0.003 (most probable / highest expected). They are arguments here.
 * @param {object} o
 * @param {number} o.z tip penetration [m]
 * @param {(z:number)=>number} o.qc cone resistance profile [Pa] as a function of depth [m]
 * @param {number} o.wallLength total plan length of skirts [m]
 * @param {number} o.t wall thickness [m]
 * @param {number} o.kp tip coefficient [-]
 * @param {number} o.kf friction coefficient [-]
 * @param {number} [o.n=200] integration intervals
 * @returns {{R:number,friction:number,tip:number}} [N]
 */
export function skirtPenetrationCpt({ z, qc, wallLength, t, kp, kf, n = 200 }) {
  need('z', z, 0);
  need('wallLength', wallLength, 0);
  need('t', t, 0);
  need('kp', kp, 0, 2);
  need('kf', kf, 0, 1);
  if (typeof qc !== 'function') throw new RangeError('qc must be a function of depth');
  let integral = 0;
  const dz = z / n;
  for (let i = 0; i < n; i++) integral += 0.5 * (qc(i * dz) + qc((i + 1) * dz)) * dz;
  const friction = 2 * wallLength * kf * integral;
  const tip = kp * wallLength * t * qc(z);
  return { R: friction + tip, friction, tip };
}

/**
 * Self-weight penetration of skirts in clay: depth at which skirtPenetrationClay equals the available
 * submerged weight, limited to the skirt height.
 * @param {object} o as skirtPenetrationClay (without z) plus
 * @param {number} o.W submerged weight available for penetration [N]
 * @param {number} o.skirtHeight [m]
 * @returns {{z:number,full:boolean,Rfull:number,ballastRequired:number}} z [m]; resistance at full
 *   penetration [N]; extra weight needed to reach full penetration [N]
 */
export function selfWeightPenetration(o) {
  const { W, skirtHeight } = o;
  need('W', W, 0);
  need('skirtHeight', skirtHeight, 1e-9);
  const R = (z) => skirtPenetrationClay({ ...o, z }).R;
  const Rfull = R(skirtHeight);
  if (Rfull <= W) return { z: skirtHeight, full: true, Rfull, ballastRequired: 0 };
  if (R(0) >= W) return { z: 0, full: false, Rfull, ballastRequired: Rfull - W };
  const z = brent((x) => R(x) - W, 0, skirtHeight, { tol: 1e-10 });
  return { z, full: false, Rfull, ballastRequired: Rfull - W };
}

/**
 * Smallest mat (fixed aspect ratio L/B) for which the classical bearing and sliding utilisations are both
 * ≤ targetUc. Solved by bracketing and root-finding on B.
 * @param {object} o all arguments of mudmatClassicalCapacity except B and L, plus
 * @param {number} [o.aspect=1] L/B [-]
 * @param {number} [o.targetUc=1] [-]
 * @param {number} [o.Bmin=0.2] [m]
 * @param {number} [o.Bmax=100] [m]
 * @returns {{B:number,L:number,ucBearing:number,ucSliding:number}} [m, m, -, -]
 */
export function requiredMatSize(o) {
  const { aspect = 1, targetUc = 1, Bmin = 0.2, Bmax = 100 } = o;
  need('aspect', aspect, 0.05, 20);
  need('targetUc', targetUc, 1e-6);
  const uc = (B) => {
    try {
      const r = mudmatClassicalCapacity({ ...o, B, L: aspect * B });
      const u = Math.max(r.ucBearing, r.ucSliding);
      return Number.isFinite(u) && u > 0 ? u : 1e6;
    } catch (e) {
      if (e instanceof RangeError) return 1e6; // eccentricity outside the base or H beyond sliding: too small
      throw e;
    }
  };
  if (uc(Bmax) > targetUc) throw new RangeError('Bmax is too small: no mat within the search range satisfies the loads');
  let lo = Bmin;
  let hi = Bmax;
  if (uc(lo) <= targetUc) hi = lo;
  else {
    for (let i = 0; i < 70; i++) {
      const mid = 0.5 * (lo + hi);
      if (uc(mid) > targetUc) lo = mid; else hi = mid;
    }
  }
  const r = mudmatClassicalCapacity({ ...o, B: hi, L: aspect * hi });
  return { B: hi, L: aspect * hi, ucBearing: r.ucBearing, ucSliding: r.ucSliding };
}

/**
 * Contact-pressure distribution under a rigid rectangular mat on a tensionless Winkler bed, bi-axial
 * eccentricity, with lift-off. The pressure is the positive part of a plane p = p0 + gx·x + gy·y (origin
 * at the centre); the three plane parameters are found from vertical and moment equilibrium by Newton
 * iteration with numerical integration over the base (midpoint rule, n × n cells).
 * Full contact reduces to p = V/A (1 ± 6 ex/B ± 6 ey/L); uniaxial lift-off to pmax = 2V / (3 L (B/2 − ex)).
 * @param {object} o
 * @param {number} o.B width along x [m]
 * @param {number} o.L length along y [m]
 * @param {number} o.V vertical load [N], > 0
 * @param {number} [o.ex=0] eccentricity along x (= My/V) [m], |ex| < B/2
 * @param {number} [o.ey=0] eccentricity along y [m], |ey| < L/2
 * @param {number} [o.n=160] cells per side
 * @returns {{pMax:number,pMin:number,p0:number,gx:number,gy:number,contactFraction:number,
 *   corners:{pp:number,pm:number,mp:number,mm:number},fullContact:boolean}}
 *   pressures [Pa] (pMin is the plane value at the least-loaded corner, negative = lift-off there);
 *   gradients [Pa/m]; contact area ratio [-]; corners at (±B/2, ±L/2)
 */
export function contactPressure({ B, L, V, ex = 0, ey = 0, n = 160 }) {
  need('B', B, 1e-9);
  need('L', L, 1e-9);
  need('V', V, 1e-12);
  need('ex', Math.abs(ex), 0, 0.5 * B * (1 - 1e-6));
  need('ey', Math.abs(ey), 0, 0.5 * L * (1 - 1e-6));
  const A = B * L;
  let p0 = V / A;
  let gx = (12 * V * ex) / (B * B * B * L);
  let gy = (12 * V * ey) / (L * L * L * B);
  const cornersOf = () => ({
    pp: p0 + (gx * B) / 2 + (gy * L) / 2, pm: p0 + (gx * B) / 2 - (gy * L) / 2,
    mp: p0 - (gx * B) / 2 + (gy * L) / 2, mm: p0 - (gx * B) / 2 - (gy * L) / 2,
  });
  let cs = cornersOf();
  let frac = 1;
  if (Math.min(cs.pp, cs.pm, cs.mp, cs.mm) < 0) {
    const dx = B / n;
    const dy = L / n;
    const dA = dx * dy;
    // Work in normalised unknowns so the 3×3 system is well conditioned.
    for (let it = 0; it < 60; it++) {
      let f0 = -V;
      let f1 = -V * ex;
      let f2 = -V * ey;
      let j00 = 0, j01 = 0, j02 = 0, j11 = 0, j12 = 0, j22 = 0;
      for (let i = 0; i < n; i++) {
        const x = -B / 2 + (i + 0.5) * dx;
        for (let j = 0; j < n; j++) {
          const y = -L / 2 + (j + 0.5) * dy;
          const p = p0 + gx * x + gy * y;
          if (p > 0) {
            f0 += p * dA; f1 += p * x * dA; f2 += p * y * dA;
            j00 += dA; j01 += x * dA; j02 += y * dA; j11 += x * x * dA; j12 += x * y * dA; j22 += y * y * dA;
          }
        }
      }
      const det = j00 * (j11 * j22 - j12 * j12) - j01 * (j01 * j22 - j12 * j02) + j02 * (j01 * j12 - j11 * j02);
      if (!(Math.abs(det) > 0)) throw new Error('contact pressure iteration failed: no contact area');
      const d0 = (f0 * (j11 * j22 - j12 * j12) - j01 * (f1 * j22 - j12 * f2) + j02 * (f1 * j12 - j11 * f2)) / det;
      const d1 = (j00 * (f1 * j22 - f2 * j12) - f0 * (j01 * j22 - j12 * j02) + j02 * (j01 * f2 - f1 * j02)) / det;
      const d2 = (j00 * (j11 * f2 - j12 * f1) - j01 * (j01 * f2 - f1 * j02) + f0 * (j01 * j12 - j11 * j02)) / det;
      p0 -= d0; gx -= d1; gy -= d2;
      frac = j00 / A;
      if (Math.abs(f0) < 1e-10 * V && Math.abs(f1) < 1e-10 * V * B && Math.abs(f2) < 1e-10 * V * L) break;
    }
    cs = cornersOf();
  }
  const vals = [cs.pp, cs.pm, cs.mp, cs.mm];
  return { pMax: Math.max(...vals), pMin: Math.min(...vals), p0, gx, gy, contactFraction: frac, corners: cs, fullContact: Math.min(...vals) >= 0 };
}

/**
 * Landing of a structure on a soil bed idealised as a linear spring (energy balance, cf. DNV-RP-N103 6.2):
 *   ½ (m + ma) v² + W′ δ = ½ k δ²  →  δ = [W′ + √(W′² + k (m + ma) v²)] / k,  Fmax = k δ.
 * @param {object} o
 * @param {number} o.mass structural mass [kg]
 * @param {number} [o.addedMass=0] hydrodynamic added mass incl. entrapped water [kg]
 * @param {number} o.W submerged weight carried by the soil at landing (hook load released) [N]
 * @param {number} o.v landing velocity [m/s]
 * @param {number} o.k vertical soil stiffness [N/m]
 * @returns {{Fmax:number,delta:number,daf:number}} peak soil reaction [N]; penetration [m]; Fmax/W′ [-]
 */
export function landingImpact({ mass, addedMass = 0, W, v, k }) {
  need('mass', mass, 0);
  need('addedMass', addedMass, 0);
  need('W', W, 0);
  need('v', v, 0);
  need('k', k, 1e-12);
  const delta = (W + Math.sqrt(W * W + k * (mass + addedMass) * v * v)) / k;
  return { Fmax: k * delta, delta, daf: W > 0 ? (k * delta) / W : Infinity };
}

/**
 * Largest set-down velocity for which the landing reaction of landingImpact stays within the allowable
 * bearing load Qallow (ultimate capacity divided by the chosen factor):
 *   v = √[ (Q² − 2 Q W′) / (k (m + ma)) ],  zero when Q ≤ 2 W′ (static sudden release already governs).
 * @param {object} o
 * @param {number} o.mass [kg]
 * @param {number} [o.addedMass=0] [kg]
 * @param {number} o.W submerged weight at landing [N]
 * @param {number} o.k soil stiffness [N/m]
 * @param {number} o.Qallow allowable soil reaction [N]
 * @returns {number} limiting velocity [m/s]
 */
export function landingVelocityLimit({ mass, addedMass = 0, W, k, Qallow }) {
  need('mass', mass, 0);
  need('W', W, 0);
  need('k', k, 1e-12);
  need('Qallow', Qallow, 0);
  const s = Qallow * Qallow - 2 * Qallow * W;
  const M = mass + addedMass;
  if (s <= 0 || M <= 0) return 0;
  return Math.sqrt(s / (k * M));
}

/**
 * Retrieval (breakout) force of a mat from clay – simple documented model:
 *   F = W′ + Fskirt + β · Nc · su · A,
 * i.e. submerged weight, skirt wall friction α su,avg As, and a fraction β of the reverse end-bearing
 * resistance mobilised as base suction. β is an argument: 0 for a vented/slowly lifted mat, up to 1 for a
 * rapid lift of a sealed base (cf. DNV-RP-N103 6.3, which allows suction up to the bearing capacity).
 * @param {object} o
 * @param {number} o.W submerged weight [N]
 * @param {number} o.A base area [m^2]
 * @param {number} o.su undrained strength at base level [Pa]
 * @param {number} o.breakoutFactor β [-], 0–1
 * @param {number} [o.Nc=6.17] reverse bearing factor incl. shape [-]
 * @param {number} [o.skirtFriction=0] wall friction on skirts [N]
 * @returns {{F:number,suction:number,suctionPressure:number}} [N, N, Pa]
 */
export function breakoutForce({ W, A, su, breakoutFactor, Nc = 1.2 * (2 + PI), skirtFriction = 0 }) {
  need('W', W, 0);
  need('A', A, 0);
  need('su', su, 0);
  need('breakoutFactor', breakoutFactor, 0, 1);
  const suction = breakoutFactor * Nc * su * A;
  return { F: W + skirtFriction + suction, suction, suctionPressure: A > 0 ? suction / A : 0 };
}

/**
 * On-bottom stability of a structure under hydrodynamic drag, inertia and lift (rigid-body statics):
 *   sliding      FS = [μ (W′ − FL) + su A + Hpassive] / FH
 *   overturning  FS = (W′ − FL) (B/2) / (FH · arm + M0)   (about the downstream edge)
 *   uplift       FS = W′ / FL
 * @param {object} o
 * @param {number} o.W submerged weight [N]
 * @param {number} o.FH horizontal hydrodynamic load (drag + inertia) [N]
 * @param {number} [o.FL=0] hydrodynamic lift [N]
 * @param {number} [o.arm=0] height of FH above the base [m]
 * @param {number} [o.M0=0] additional overturning moment [N·m]
 * @param {number} o.B base width in the load direction [m]
 * @param {number} [o.mu=0] drained friction coefficient tan δ [-]
 * @param {number} [o.su=0] undrained interface strength [Pa]
 * @param {number} [o.A=0] base area [m^2]
 * @param {number} [o.Hpassive=0] skirt/embedment resistance [N]
 * @returns {{fsSliding:number,fsOverturning:number,fsUplift:number,slidingResistance:number,e:number}}
 *   [-,-,-,N,m]; e eccentricity of the base reaction [m]
 */
export function onBottomStability({ W, FH, FL = 0, arm = 0, M0 = 0, B, mu = 0, su = 0, A = 0, Hpassive = 0 }) {
  need('W', W, 0);
  need('FH', FH, 0);
  need('FL', FL, 0);
  need('B', B, 1e-9);
  need('mu', mu, 0);
  const N = Math.max(W - FL, 0);
  const R = mu * N + su * A + Hpassive;
  const Mo = FH * arm + M0;
  return {
    fsSliding: FH > 0 ? R / FH : Infinity,
    fsOverturning: Mo > 0 ? (N * B) / (2 * Mo) : Infinity,
    fsUplift: FL > 0 ? W / FL : Infinity,
    slidingResistance: R,
    e: N > 0 ? Mo / N : Infinity,
  };
}
